from fwat.measure.FwatPreOP import FwatPreOP
import numpy as np 
from mpi4py import MPI

class Tele_PreOP(FwatPreOP):
    def __init__(self, measure_type, iter, evtid, run_opt):
        super().__init__(measure_type, iter, evtid, run_opt)
    
        # read tref and t_inj for teleseismic events
        from .tele.tele import get_injection_time
        from .tele.tele import compute_ak135_time
        self.t_inj = get_injection_time(self.evtid)

        if "_" in self.evtid:
            phase = self.evtid.split('_')[0]
        else:
            phase = 'P'
        self.t_ref = compute_ak135_time(
            self.evla,self.evlo,self.evdp,
            self.stla,self.stlo,phase
        )
    
    def _sanity_check(self):
        super()._sanity_check()

        # make sure adjsrc_type in [2,'cross-conv']
        if self.adjsrc_type not in ['cross-conv','2']:
            if self.myrank == 0:
                print("only L2 norm (imeas = 2) and cross-conv adjoint source are supported!")
                print(f"adjsrc_type = {self.adjsrc_type}")
            exit(1)
        
        # make sure Z/R components are used
        if self.components != ['R','Z']:
            if self.myrank == 0:
                print("R/Z components must be used!")
                print(f"components = {self.components}")
            exit(1)

    def save_forward(self):
        import os 
        from obspy.io.sac import SACTrace
        from scipy.signal import convolve

        # get some vars
        evtid = self.evtid 
        ncomp = self.ncomp
        components = self.components
        dt_syn = self.dt_syn
        npt_syn = self.npt_syn
        t0_syn = self.t0_syn
        
        # init a sac header
        tr = SACTrace(
            evla=self.evla,evlo=self.evlo,
            evdp=self.evdp,stla=0.,
            stlo=0.,stel=0,lcalda=1,
            delta = dt_syn
        )

        outdir = f"{self.DATA_DIR}/{evtid}"
        if self.myrank == 0:
            print("Synthetic Observations ...")
            os.makedirs(outdir,exist_ok=True)
        MPI.COMM_WORLD.Barrier()
        
        # load stf for tele seismic events
        stf = np.zeros((ncomp,npt_syn))
        for ib in range(len(self.Tmax)):
            bandname = self._get_bandname(ib)
            for ic in range(self.ncomp):
                ch = self.components[ic]
                tr = SACTrace.read(f'src_rec/stf_{ch}.sac.{bandname}_{evtid}')
                stf[ic,:] += tr.data 
        
        # loop every station to save sac
        for ir in range(self.nsta_loc):
            i = self._istart + ir 
            for ic in range(ncomp):
                ch = components[ic]

                # load synthetics from npy
                code = self._get_station_code(i,ic)
                filename = f"{self.syndir}/OUTPUT_FILES/{code}.sem.npy"
                data =  np.load(filename)

                tr.data = convolve(data[:,1],stf[ic,:],'same') * dt_syn   
                tr.b =  self.t_inj - self.t_ref[i]

                # channel and others
                tr.kcmpnm = f"{self.chcode}{ch}"
                tr.knetwk = self.netwk[i]
                tr.kstnm = self.stnm[i]
                tr.stla = self.stla[i]
                tr.stlo = self.stlo[i]

                # save to sac
                filename = f"{outdir}/{code}.sac"
                tr.write(filename)

    def _process_all_seismograms(self,ib:int):
        """
        do preprocessing for all seismigrams at a given frequency band, save in big arrays

        Parameters
        ------------
        ib: int
            frequency band index
        """
        from obspy.io.sac import SACTrace
        from .utils import interpolate_syn,bandpass

        # frequency band
        freqmin = 1. / self.Tmax[ib]
        freqmax = 1. / self.Tmin[ib]

        # get frequency band
        out_dir = f"{self.syndir}/OUTPUT_FILES/"

        # get time window 
        dt_syn = self.dt_syn
        win_tb,win_te = self.pdict['TIME_WINDOW']
        npt2 = int((win_te + win_tb) / dt_syn)
        if npt2 // 2 * 2 != npt2:
            npt2 += 1
        
        # glob arrays
        ncomp = self.ncomp
        nsta = self.nsta 
        npt_syn = self.npt_syn
        glob_obs = np.zeros((nsta,ncomp,npt_syn))
        glob_syn = np.zeros((nsta,ncomp,npt_syn))

        # loop each station
        for ir in range(self.nsta_loc):
            i = ir + self._istart
            for ic in range(ncomp):
                name = self._get_station_code(i,ic)

                # read data
                syn_data = np.load(f"{out_dir}/{name}.sem.npy")[:,1]
                obs_tr = SACTrace.read(f"{self.DATA_DIR}/{self.evtid}/{name}.sac")
                t0_obs = obs_tr.b 
                dt_obs = obs_tr.delta 
                npt_obs = obs_tr.npts 
                obs_data = obs_tr.data

                # filter 
                obs_data = bandpass(obs_data,dt_obs,freqmin,freqmax)
                syn_data = bandpass(syn_data,dt_syn,freqmin,freqmax)

                # interpolate the obs/syn data to the same sampling of syn data
                t0_inp = self.t_ref[i] - win_tb
                u1 = interpolate_syn(obs_data,t0_obs + self.t_ref[i],dt_obs,npt_obs,t0_inp,dt_syn,npt2)
                w1 = interpolate_syn(syn_data,self.t_inj,dt_syn,npt_syn,t0_inp,dt_syn,npt2)
                u = interpolate_syn(u1,t0_inp,dt_syn,npt2,self.t_inj,dt_syn,npt_syn)
                w = interpolate_syn(w1,t0_inp,dt_syn,npt2,self.t_inj,dt_syn,npt_syn)     

                # save to  global array   
                glob_obs[i,ic,:] = u
                glob_syn[i,ic,:] = w

        # save everything in big arrays
        tmp = MPI.COMM_WORLD.allreduce(glob_obs); glob_obs = tmp * 1.
        tmp = MPI.COMM_WORLD.allreduce(glob_syn); glob_syn = tmp * 1.

        return glob_obs,glob_syn

    def _cal_adj_source_l2(self,ib:int):
        from obspy.io.sac import SACTrace
        from .utils import interpolate_syn
        from fwat.adjoint.l2_misfit import measure_adj_l2
        from scipy.signal import convolve,correlate
        from .tele.tele import compute_stf,get_average_amplitude
        import os 

        # get frequency band
        bandname = self._get_bandname(ib)
        freqmin = 1. / self.Tmax[ib]
        freqmax = 1. / self.Tmin[ib]
        out_dir = f"{self.syndir}/OUTPUT_FILES/"

        # get vars
        t0_syn = self.t0_syn
        dt_syn = self.dt_syn
        npt_syn = self.npt_syn
        t_inj = self.t_inj

        # get time window 
        win_tb,win_te = self.pdict['TIME_WINDOW']
        npt2 = int((win_te + win_tb) / dt_syn)
        if npt2 // 2 * 2 != npt2:
            npt2 += 1
        
        # allocate global arrays 
        nsta = self.nsta 
        nsta_loc = self.nsta_loc

        # misfits
        ncomp = self.ncomp
        tstart = np.zeros((nsta_loc))
        tend = np.zeros((nsta_loc))
        win_chi = np.zeros((nsta_loc,ncomp,20))
        tr_chi = np.zeros((nsta_loc,ncomp))
        am_chi = np.zeros((nsta_loc,ncomp))

        # glob arrays
        ncomp = self.ncomp
        glob_obs,glob_syn = self._process_all_seismograms(ib)

        # get source time function 
        stf_names = [f'{self.SRC_REC}/stf_{ch}.sac.' + f"{bandname}" + f"_{self.evtid}" for ch in self.components]
        has_stf_flag = True 
        for ic in range(ncomp):
            has_stf_flag = has_stf_flag and os.path.isfile(stf_names[ic])
        if has_stf_flag: # read stf from directory
            if self.myrank == 0:
                print(f"read stf from {self.SRC_REC} ...")
            stf = np.zeros((ncomp,npt_syn),'f4')
            for ic in range(ncomp):
                stf[ic,:] = SACTrace.read(stf_names[ic]).data
        else: # estimate stf
            if self.myrank == 0:
                print("estimation stf by PCA ...")
            stf = compute_stf(glob_syn,glob_obs,dt_syn,freqmin,freqmax,self.components)

            # save to SRC_REC
            sac_head = SACTrace(b=0,npts=npt_syn,delta=dt_syn)
            for ic in range(ncomp):
                sac_head.data = stf[ic,:]
                if self.myrank ==0 : sac_head.write(stf_names[ic])
        
        # get average amplitude
        ic = self.components.index('Z')
        avgamp = get_average_amplitude(glob_obs,ic)
        if self.myrank == 0: 
            print("average amplitude for data gather: %g\n" %(avgamp))
    
        # call measure_adj for misfits
        for ir in range(nsta_loc):
            i = ir + self._istart 
            for ic in range(self.ncomp):
                # convolve with stf
                tmp = dt_syn * convolve(glob_syn[i,ic,:],stf[ic,:],'same')

                # only keep data in the window
                t0_inp = self.t_ref[i] - win_tb
                tmp1 = interpolate_syn(tmp,t_inj,dt_syn,npt_syn,t0_inp,dt_syn,npt2)
                tmp = interpolate_syn(tmp1,t0_inp,dt_syn,npt2,t_inj,dt_syn,npt_syn)
                glob_syn[i,ic,:] = tmp

                # time window
                tstart[ir] = self.t_ref[i] - self.t_inj - win_tb
                tend[ir] = self.t_ref[i] - self.t_inj + win_te

                # measure

                tr_chi[ir,ic],am_chi[ir,ic],win_chi[ir,ic,:],adjsrc =  \
                    measure_adj_l2(
                        glob_obs[i,ic,:],glob_syn[i,ic,:],
                        t0_syn,dt_syn,npt_syn,
                        tstart[ir],tend[ir]
                    )

                # contributions on stf
                adjsrc /= avgamp 
                tmp = correlate(adjsrc,stf[ic,:],'same') * dt_syn 
                if np.max(np.abs(tmp)) > 0.:
                    adjsrc = tmp / np.max(np.abs(stf[ic,:]))

                # save adjoint source
                data = np.zeros((npt_syn,2))
                data[:,0] = t0_syn + np.arange(npt_syn) * dt_syn
                data[:,1] = adjsrc
                name = self._get_station_code(i,ic) + ".adj.sem.npy"
                np.save(f"{out_dir}/{bandname}/{name}",data)

                # save SAC obs and syn
                # init a sac header
                tr = SACTrace(
                    evla=self.evla,evlo=self.evlo,
                    evdp=self.evdp,stla=self.stla[i],
                    stlo=self.stlo[i],stel=0,lcalda=1,
                    delta = dt_syn,b=t0_syn,
                    t0 = self.t_ref[i] - t_inj,
                    knetwk=self.netwk[i],
                    kstnm = self.stnm[i],
                    kcmpnm = self.chcode + self.components[ic]
                )
                name = self._get_station_code(i,ic)
                tr.data = glob_obs[i,ic,:]
                tr.write(f"{out_dir}/{bandname}/{name}.sac.obs")
                tr.data = glob_syn[i,ic,:]
                tr.write(f"{out_dir}/{bandname}/{name}.sac.syn")
        
        # normalize misfit function 
        tr_chi *= 1. / avgamp**2 
        am_chi *= 1. / avgamp**2
        win_chi[:,:,14] *= 1. / avgamp**2

        self._print_measure_info(bandname,tstart,tend,tr_chi,am_chi,win_chi)

    def _cal_adj_source_user(self,ib:int):
        """
        compute adjoint source by using cross-convolution for a given frequency band

        Parameters
        ------------
        ib: int
            frequency band index
        """
        from obspy.io.sac import SACTrace
        from .tele.tele import get_average_amplitude
        from fwat.adjoint.cross_conv import measure_adj_cross_conv

        # get frequency band
        bandname = self._get_bandname(ib)
        out_dir = f"{self.syndir}/OUTPUT_FILES/"

        # get vars
        t0_syn = self.t0_syn
        dt_syn = self.dt_syn
        npt_syn = self.npt_syn
        t_inj = self.t_inj

        # get time window 
        win_tb,win_te = self.pdict['TIME_WINDOW']
        npt2 = int((win_te + win_tb) / dt_syn)
        if npt2 // 2 * 2 != npt2:
            npt2 += 1
        
        # allocate global arrays 
        nsta_loc = self.nsta_loc

        # misfits
        tstart = np.zeros((nsta_loc))
        tend = np.zeros((nsta_loc))
        win_chi = np.zeros((nsta_loc,1,20))
        tr_chi = np.zeros((nsta_loc,1))
        am_chi = np.zeros((nsta_loc,1))

        # glob arrays
        glob_obs,glob_syn = self._process_all_seismograms(ib)

        # get average amplitude of obs data
        ic = self.components.index('Z')
        avgamp = get_average_amplitude(glob_obs,ic)
        if self.myrank == 0: 
            print("average amplitude for data gather: %g\n" %(avgamp))

        # normalize 
        glob_obs /= avgamp
        glob_syn /= avgamp

        # measure
        for ir in range(nsta_loc):
            i = ir + self._istart

            # get  window info
            tstart[ir] = self.t_ref[i] - self.t_inj - win_tb
            tend[ir] = self.t_ref[i] - self.t_inj + win_te
    
            # call exponentiated phase function
            ic_z = self.components.index('Z')
            ic_r = self.components.index('R')
            tr_chi[ir, 0],am_chi[ir, 0],  \
            win_chi[ir, 0, :],adj_r,adj_z,cc1,cc2 = \
                measure_adj_cross_conv(
                    glob_obs[i, ic_z, :],
                    glob_syn[i, ic_z, :],
                    glob_obs[i, ic_r, :],
                    glob_syn[i, ic_r, :],
                    t0_syn,
                    dt_syn,
                    self.Tmin[ib],
                    self.Tmax[ib],
                    tstart[ir],
                    tend[ir],
                )
        
            # save adjoint source
            data = np.zeros((npt_syn,2))
            data[:,0] = t0_syn + np.arange(npt_syn) * dt_syn
            data[:,1] = adj_z
            name = self._get_station_code(i,ic_z) + ".adj.sem.npy"
            np.save(f"{out_dir}/{bandname}/{name}",data)
            data[:,1] = adj_r
            name = self._get_station_code(i,ic_r) + ".adj.sem.npy"
            np.save(f"{out_dir}/{bandname}/{name}",data)

            # save obs and synthetic data
            for ic in range(1): # only one component
                # init a sac header
                tr = SACTrace(
                    evla=self.evla,evlo=self.evlo,
                    evdp=self.evdp,stla=self.stla[i],
                    stlo=self.stlo[i],stel=0,lcalda=1,
                    delta = dt_syn,b=t0_syn,
                    t0 = self.t_ref[i] - t_inj,
                    knetwk=self.netwk[i],
                    kstnm = self.stnm[i],
                    kcmpnm = self.chcode + self.components[ic]
                )
                name = self._get_station_code(i,ic)
                tr.data = cc1
                tr.write(f"{out_dir}/{bandname}/{name}.sac.obs")
                tr.data = cc2
                tr.write(f"{out_dir}/{bandname}/{name}.sac.syn")

        self._print_measure_info(bandname,tstart,tend,tr_chi,am_chi,win_chi)
    
    def cal_adj_source(self,ib:int):
        import os 
        bandname = self._get_bandname(ib)
        if self.myrank == 0:
            print(f"preprocessing for band {bandname} ...")
            os.makedirs(f"{self.syndir}/OUTPUT_FILES/{bandname}",exist_ok=True)
        MPI.COMM_WORLD.Barrier()

        if self.adjsrc_type == '2':
            self._cal_adj_source_l2(ib)
        else:
            self._cal_adj_source_user(ib)
