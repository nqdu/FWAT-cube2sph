from .FwatPreOP import FwatPreOP
import numpy as np 
from mpi4py import MPI
from fwat.adjoint.MeasureStats import MeasureStats
    
def _splitting_intensity(Rsyn:np.ndarray,Tsyn:np.ndarray,dt_syn:float):
    """
    compute splitting intensity from R and T components
    
    Parameters
    ---------------
    Rsyn,Tsyn: np.ndarray
        R and T components of synthetics
    dt_syn: float
        time step of synthetics
    
    Returns
    ---------------
    si_syn: float
        splitting intensity of synthetics
    """

    from .utils import dif1
    dRsyn = dif1(Rsyn,dt_syn)
    norm_syn = np.sum(dRsyn**2)
    norm_syn = 1. / (norm_syn + 1.0e-30)
    si_syn = -2. * np.sum(dRsyn * Tsyn) * norm_syn 

    return float(si_syn)

def _splitting_intensity_adjsrc(Rsyn,Tsyn,dt_syn,si_obs = 0.,weight = 1.):
    """
    compute splitting intensity from R and T components and its adjoint source
    
    Parameters
    ---------------
    Rsyn,Tsyn: np.ndarray
        R and T components of synthetics
    dt_syn: float
        time step of synthetics
    si_obs: float
        splitting intensity of observations
    weight: float
        weight for this station
    
    Returns
    ---------------
    si_syn: float
        splitting intensity of synthetics
    adjsrc_R,adjsrc_T: np.ndarray
        adjoint source for R and T components"""
    from .utils import dif1
    dRsyn = dif1(Rsyn,dt_syn)
    norm_syn = np.sum(dRsyn**2)
    norm_syn = 1. / (norm_syn + 1.0e-30)
    si_syn = -2. * np.sum(dRsyn * Tsyn) * norm_syn 

    si_diff = (si_syn - si_obs) * weight
    dTcomp = dif1(Tsyn,dt_syn)
    ddRsyn = dif1(dRsyn,dt_syn)
    adjsrc_T = -2. * si_diff * dRsyn * norm_syn
    adjsrc_R = -2. * si_diff * (  
        2. * np.sum(dRsyn * Tsyn) * norm_syn**2 * ddRsyn - 
        dTcomp * norm_syn
    )
    return si_syn,adjsrc_R,adjsrc_T

class SKS_PreOP(FwatPreOP):
    def __init__(self, measure_type, iter, evtid, run_opt):
        super().__init__(measure_type, iter, evtid, run_opt)
    
        # read tref and t_inj for teleseismic events
        from .tele.tele import get_injection_time
        from .tele.tele import compute_ak135_time
        self.t_inj = get_injection_time(self.evtid)

        if "_" in self.evtid:
            phase = self.evtid.split('_')[0]
        else:
            phase = 'sks'
        self.t_ref = compute_ak135_time(
            self.evla,self.evlo,self.evdp,
            self.stla,self.stlo,phase
        )

    def _write_si_obs(self,si_syn:np.ndarray):
        # save si
        for irank in range(self.nprocs):
            if irank == self.myrank:
                outfile = f"{self.DATA_DIR}/{self.evtid}/si_obs.txt"

                # open file
                if irank == 0:
                    fio = open(outfile,"w")
                else:
                    fio = open(outfile,"a")
                
                for ir in range(self.nsta_loc):
                    i = ir + self._istart
                    name = f"{self.netwk[i]}.{self.stnm[i]}"
                    fio.write(f"{name} %f 1. \n" %(si_syn[ir]))
                
                # close
                fio.close()
                
            # sync
            MPI.COMM_WORLD.Barrier()

    def _sanity_check(self):
        super()._sanity_check()

        assert self.adjsrc_type in ["SI","cross-conv"], f"Invalid ADJSRC_TYPE: {self.adjsrc_type}"

        if self.components != ["R","T"]:
            print("For SKS FWI, both R and T components are required!")
            exit(1)

    def save_forward(self):
        import os 
        from obspy.io.sac import SACTrace

        # get some vars
        evtid = self.evtid 
        ncomp = self.ncomp
        components = self.components
        dt_syn = self.dt_syn
        npt_syn = self.npt_syn

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
        
        # loop every station to save sac
        for ir in range(self.nsta_loc):
            i = self._istart + ir 

            for ic in range(ncomp):
                ch = components[ic]

                # load synthetics from npy
                code = self._get_station_code(i,ic)
                filename = f"{self.syndir}/OUTPUT_FILES/{code}.sem.npy"
                data = self.seismogram[filename]
                #data =  np.load(filename)

                tr.b =  self.t_inj - self.t_ref[i]
                tr.data = data[:,1] * 1.

                # channel and others
                tr.kcmpnm = f"{self.chcode}{ch}"
                tr.knetwk = self.netwk[i]
                tr.kstnm = self.stnm[i]
                tr.stla = self.stla[i]
                tr.stlo = self.stlo[i]

                # save to sac
                filename = f"{outdir}/{code}.sac"
                tr.write(filename)

    def cal_adj_source_cc(self,ib:int):
        from obspy.io.sac import SACTrace
        from .utils import interpolate_syn
        from .utils import bandpass 
        from fwat.adjoint.cross_conv import measure_adj_cross_conv
        import os 

        bandname = self._get_bandname(ib)
        if self.myrank == 0:
            print(f"preprocessing for band {bandname} ...")
            #os.makedirs(f"{self.syndir}/OUTPUT_FILES/{bandname}",exist_ok=True)
        MPI.COMM_WORLD.Barrier()

        freqmin = 1. / self.Tmax[ib]
        freqmax = 1. / self.Tmin[ib]
        out_dir = f"{self.syndir}/OUTPUT_FILES"

        # get vars
        t0_syn = self.t0_syn
        dt_syn = self.dt_syn
        npt_syn = self.npt_syn
        t_inj = self.t_inj

        # get time window 
        win_tb,win_te = self.pdict['TIME_WINDOW']
        
        # allocate global arrays 
        nsta_loc = self.nsta_loc

        # misfits
        stats_list = []

        # save seismo_win headers
        self.seismo_win['dt'] = dt_syn
        self.seismo_win['t0'] = t0_syn

        # get normalize amplitude for syn/obs data
        amp_file = f"{self.SRC_REC}/sks.factor.{bandname}.{self.evtid}.txt"
        if os.path.exists(amp_file):
            if self.myrank == 0:
                print("read average amplitude for syn/obs data from file ...")
            
            with open(amp_file,'r') as f:
                avgamp,avgamp_syn = map(float,f.read().strip().split())
        else:
            # get average amplitude of obs data
            ic_r = self.components.index("R")

            maxamp_obs = np.zeros((self.nsta))
            maxamp_syn = np.zeros((self.nsta))
            for ir in range(nsta_loc):
                i = ir + self._istart
                # time window
                tstart = self.t_ref[i] - t_inj - win_tb
                tend = self.t_ref[i] - t_inj + win_te
                name = self._get_station_code(i,ic_r)

                # read syn data and preprocessing
                #sdata = np.load(f"{out_dir}/{name}.sem.npy")[:,1]
                sdata = self.seismogram[f"{out_dir}/{name}.sem.npy"][:,1] * 1.
                sdata = bandpass(sdata,dt_syn,freqmin,freqmax)

                # read obs data
                obs_tr = SACTrace.read(f"{self.DATA_DIR}/{self.evtid}/{name}.sac")
                # obs info
                t0_obs = obs_tr.b 
                dt_obs = obs_tr.delta 
                npt_obs = obs_tr.npts 

                # filter 
                obs_tr.data = bandpass(obs_tr.data,dt_obs,freqmin,freqmax)
                
                # interpoate obs data to same series of synthetics
                odata = interpolate_syn(obs_tr.data,t0_obs + self.t_ref[i],dt_obs,npt_obs,t_inj,dt_syn,npt_syn)

                # save to local array
                maxamp_obs[i] = np.max(np.abs(odata))
                maxamp_syn[i] = np.max(np.abs(sdata))
            
            # gather max amplitude
            maxamp_obs = MPI.COMM_WORLD.allreduce(maxamp_obs, op=MPI.SUM)
            maxamp_syn = MPI.COMM_WORLD.allreduce(maxamp_syn, op=MPI.SUM)

            # get average amplitude
            avgamp = np.median(maxamp_obs)
            avgamp_syn = np.median(maxamp_syn)
            min_factor = 0.2
            max_factor = 5.0
            valid_mask = (maxamp_obs > min_factor * avgamp) & \
                        (maxamp_obs < max_factor * avgamp)
            valid_mask1 = (maxamp_syn > min_factor * avgamp_syn) & \
                        (maxamp_syn < max_factor * avgamp_syn)
            if np.sum(valid_mask) > 0:
                avgamp = np.mean(maxamp_obs[valid_mask])
            if np.sum(valid_mask1) > 0:
                avgamp_syn = np.mean(maxamp_syn[valid_mask1])

            # save to file
            if self.myrank == 0:
                with open(amp_file,'w') as f:
                    f.write(f"{avgamp} {avgamp_syn}\n")
        MPI.COMM_WORLD.Barrier()

        # loop each station
        for ir in range(nsta_loc):
            i = ir + self._istart
            obs_data = np.zeros((2,npt_syn))
            syn_data = np.zeros((2,npt_syn))

            # time window
            tstart = self.t_ref[i] - t_inj - win_tb
            tend = self.t_ref[i] - t_inj + win_te

            for ic in range(self.ncomp):
                name = self._get_station_code(i,ic)

                # read syn data and preprocessing
                #sdata = np.load(f"{out_dir}/{name}.sem.npy")[:,1]
                sdata = self.seismogram[f"{out_dir}/{name}.sem.npy"][:,1] * 1.
                syn_data[ic,:] = bandpass(sdata,dt_syn,freqmin,freqmax)

                # read obs data
                obs_tr = SACTrace.read(f"{self.DATA_DIR}/{self.evtid}/{name}.sac")
                # obs info
                t0_obs = obs_tr.b 
                dt_obs = obs_tr.delta 
                npt_obs = obs_tr.npts 

                # filter 
                obs_tr.data = bandpass(obs_tr.data,dt_obs,freqmin,freqmax)
                
                # interpoate obs data to same series of synthetics
                odata = interpolate_syn(obs_tr.data,t0_obs + self.t_ref[i],dt_obs,npt_obs,t_inj,dt_syn,npt_syn)

                # save to global array   
                obs_data[ic,:] = odata
            
            # components
            ic_r = self.components.index("R")
            ic_t = self.components.index("T")

            # normalize syn and obs data by average amplitude
            obs_data *= 1. / avgamp
            syn_data *= 1. / avgamp_syn

            # compute misfit and adjoint source
            stats,adj_r,adj_t,cc1,cc2 = \
                measure_adj_cross_conv(
                    obs_data[ic_r,:],
                    syn_data[ic_r,:],
                    obs_data[ic_t,:],
                    syn_data[ic_t,:],
                    t0_syn,dt_syn,
                    tstart,tend + 1. / freqmax
                )
            stats.code = f"{self.netwk[i]}.{self.stnm[i]}.cross-conv"
            stats_list.append(stats)

            # filter adjoint source
            adj_r = bandpass(adj_r,dt_syn,freqmin,freqmax) / avgamp_syn
            adj_t = bandpass(adj_t,dt_syn,freqmin,freqmax) / avgamp_syn

            data = np.zeros((npt_syn,2))
            data[:,0] = t0_syn + np.arange(npt_syn) * dt_syn
            data[:,1] = adj_t
            outname = f"{out_dir}/{bandname}/{self.netwk[i]}.{self.stnm[i]}.{self.chcode}T.adj.sem.npy"
            #np.save(outname,data)
            self.seismogram_adj[outname] = data * 1.

            data[:,1] = adj_r
            outname = f"{out_dir}/{bandname}/{self.netwk[i]}.{self.stnm[i]}.{self.chcode}R.adj.sem.npy"
            #np.save(outname,data)
            self.seismogram_adj[outname] = data * 1.

            # save syn
            name = self._get_station_code(i,ic_r)
            self.seismo_win[f"{out_dir}/{bandname}/{name}.dat.syn"] = cc1.copy()
            self.seismo_win[f"{out_dir}/{bandname}/{name}.dat.obs"] = cc2.copy()
            if ir == 0:
                self.seismo_win['npts'] = len(cc1)

    
        self._print_measure_info(bandname,stats_list=stats_list)
    
    def cal_adj_source_si(self,ib:int):
        from obspy.io.sac import SACTrace
        from .utils import interpolate_syn
        from .utils import bandpass,taper_window
        import os 

        bandname = self._get_bandname(ib)
        if self.myrank == 0:
            print(f"preprocessing for band {bandname} ...")
            #os.makedirs(f"{self.syndir}/OUTPUT_FILES/{bandname}",exist_ok=True)
        MPI.COMM_WORLD.Barrier()

        freqmin = 1. / self.Tmax[ib]
        freqmax = 1. / self.Tmin[ib]
        out_dir = f"{self.syndir}/OUTPUT_FILES"

        # get vars
        t0_syn = self.t0_syn
        dt_syn = self.dt_syn
        npt_syn = self.npt_syn
        t_inj = self.t_inj

        # get time window 
        win_tb,win_te = self.pdict['TIME_WINDOW']
        
        # allocate global arrays 
        nsta_loc = self.nsta_loc

        # misfits
        stats_list = []

        # load si if it exists
        cal_obs_si = False
        filename = f"{self.DATA_DIR}/{self.evtid}/si_obs.txt"
        si_obs_data = np.zeros((nsta_loc))
        if os.path.exists(filename):
            tmp = np.loadtxt(filename,dtype=str,ndmin=2)
            si_obs_evt = {}
            for i in range(tmp.shape[0]):
                si_obs_evt[tmp[i,0]] = tmp[i,1:].astype(float)
        else:
            si_obs_evt = {}
            cal_obs_si = True

        # save seismo_win headers
        self.seismo_win['dt'] = dt_syn
        self.seismo_win['t0'] = t0_syn
        self.seismo_win['npts'] = npt_syn

        # loop each station
        for ir in range(nsta_loc):
            i = ir + self._istart
            obs_data = np.zeros((2,npt_syn))
            syn_data = np.zeros((2,npt_syn))

            # time window
            taper = np.zeros((npt_syn))
            tstart = self.t_ref[i] - t_inj - win_tb
            tend = self.t_ref[i] - t_inj + win_te
            lpt,rpt,taper0 = taper_window(0,dt_syn,npt_syn,tstart,tend)
            taper[lpt:rpt] = taper0 * 1.

            for ic in range(self.ncomp):
                name = self._get_station_code(i,ic)

                # read syn data and preprocessing
                #sdata = np.load(f"{out_dir}/{name}.sem.npy")[:,1]
                sdata = self.seismogram[f"{out_dir}/{name}.sem.npy"][:,1] * 1.
                sdata = bandpass(sdata,dt_syn,freqmin,freqmax)
                syn_data[ic,:] = sdata * taper 

                # obs data
                if cal_obs_si:
                    obs_tr = SACTrace.read(f"{self.DATA_DIR}/{self.evtid}/{name}.sac")
                    # obs info
                    t0_obs = obs_tr.b 
                    dt_obs = obs_tr.delta 
                    npt_obs = obs_tr.npts 

                    # filter 
                    obs_tr.data = bandpass(obs_tr.data,dt_obs,freqmin,freqmax)
                    
                    # interpoate obs data to same series of synthetics
                    odata = interpolate_syn(obs_tr.data,t0_obs + self.t_ref[i],dt_obs,npt_obs,t_inj,dt_syn,npt_syn)

                    # save to global array   
                    obs_data[ic,:] = odata * taper
                
                # save syn
                self.seismo_win[f"{out_dir}/{bandname}/{name}.dat.syn"] = syn_data[ic,:].copy()

                # save obs
                if cal_obs_si:
                    self.seismo_win[f"{out_dir}/{bandname}/{name}.dat.obs"] = obs_data[ic,:].copy()
                    #tr.write(f"{out_dir}/{bandname}/{name}.sac.obs")

            # components
            ic_r = self.components.index("R")
            ic_t = self.components.index("T")

            # compute si for obs if required
            weight = 1.
            if cal_obs_si:
                SI_obs = _splitting_intensity(
                    obs_data[ic_r,:],
                    obs_data[ic_t,:],
                    dt_syn
                )

                # save to it
                si_obs_data[ir] = SI_obs
            else:
                name = f"{self.netwk[i]}.{self.stnm[i]}"
                if name not in si_obs_evt:
                    print(f"cannot locate {name} in  {self.DATA_DIR}/{self.evtid}!")
                    exit(1)
                tmp = si_obs_evt[name]
                SI_obs = float(tmp[0])
                if len(tmp) > 1:
                    weight = 1. / float(tmp[1])

            # compute si for syn data and adjoint source
            SI_syn,adjsrc_R,adjsrc_T = _splitting_intensity_adjsrc(
                syn_data[ic_r,:],
                syn_data[ic_t,:],
                dt_syn,
                si_obs=SI_obs,
                weight=weight
            )

            # compute misfit function
            si_diff = (SI_syn - SI_obs) * weight
            L = 0.5 * si_diff**2
            stats = MeasureStats(
                adj_type = self.adjsrc_type,
                misfit = L,
                tstart = tstart,
                tend = tend,
                tr_chi = L,
                am_chi = L,
                tshift = 0.
            )
            stats.code = f"{self.netwk[i]}.{self.stnm[i]}"
            stats_list.append(stats)

            # filter adjoint source and taper it
            adjsrc_R = bandpass(adjsrc_R,dt_syn,freqmin,freqmax) * taper 
            adjsrc_T = bandpass(adjsrc_T,dt_syn,freqmin,freqmax) * taper 

            data = np.zeros((npt_syn,2))
            data[:,0] = t0_syn + np.arange(npt_syn) * dt_syn
            data[:,1] = adjsrc_T
            outname = f"{out_dir}/{bandname}/{self.netwk[i]}.{self.stnm[i]}.{self.chcode}T.adj.sem.npy"
            #np.save(outname,data)
            self.seismogram_adj[outname] = data * 1.

            data[:,1] = adjsrc_R
            outname = f"{out_dir}/{bandname}/{self.netwk[i]}.{self.stnm[i]}.{self.chcode}R.adj.sem.npy"
            #np.save(outname,data)
            self.seismogram_adj[outname] = data * 1.

        # save SI 
        if cal_obs_si: self._write_si_obs(si_obs_data)
    
        self._print_measure_info(bandname,stats_list=stats_list)

    
    def cal_adj_source(self, ib: int):
        """
        calculate adjoint source for SKS FWI
        Parameters
        ----------
        ib : int
            index of frequency band
        Returns
        -------
        None
        """
        if self.adjsrc_type == "SI":
            self.cal_adj_source_si(ib)
        elif self.adjsrc_type == "cross-conv":
            self.cal_adj_source_cc(ib)
        else:
            print(f"Invalid ADJSRC_TYPE: {self.adjsrc_type}")
            exit(1)
        