from fwat.measure.FwatPreOP import FwatPreOP
import numpy as np 
from mpi4py import MPI

class Noise_PreOP(FwatPreOP):
    def __init__(self, measure_type, iter, evtid, run_opt):
        super().__init__(measure_type, iter, evtid, run_opt)

    def save_forward(self):
        import os 
        from obspy.io.sac import SACTrace
        from .utils import cumtrapz1

        # get some vars
        evtid = self.evtid 
        ncomp = self.ncomp
        components = self.components
        dt_syn = self.dt_syn
        npt_syn = self.npt_syn
        t0_syn = self.t0_syn
        myrank = self.myrank
        
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
                data =  np.load(filename)

                if self.pdict['USE_EGF'] == False:
                    if ir == 0 and myrank == 0 and ic == 0:
                        print("EGF => CCF ...")
                    tr.data = -cumtrapz1(data[:,1],dt_syn)
                else:
                    tr.data = data[:,1] * 1.
                tr.b = t0_syn

                # channel and others
                tr.kcmpnm = f"{self.chcode}{ch}"
                tr.knetwk = self.netwk[i]
                tr.kstnm = self.stnm[i]
                tr.stla = self.stla[i]
                tr.stlo = self.stlo[i]

                # save to sac
                filename = f"{outdir}/{code}.sac"
                tr.write(filename)
    
    def cal_adj_source(self,ib:int):
        from obspy.io.sac import SACTrace
        from .utils import interpolate_syn,dif1
        from .utils import bandpass
        from .measure import measure_adj
        import os 
        bandname = self._get_bandname(ib)
        if self.myrank == 0:
            print(f"preprocessing for band {bandname} ...")
            os.makedirs(f"{self.syndir}/OUTPUT_FILES/{bandname}",exist_ok=True)
        MPI.COMM_WORLD.Barrier()

        freqmin = 1. / self.Tmax[ib]
        freqmax = 1. / self.Tmin[ib]
        out_dir = f"{self.syndir}/OUTPUT_FILES/"

        # new dt/nt for interpolate
        dt_inp = 0.01
        t0_inp = -10.

        # get vars
        t0_syn = self.t0_syn
        dt_syn = self.dt_syn
        npt_syn = self.npt_syn

        # group velocity window
        vmin_list = [x[0] for x in self.pdict['GROUPVEL_WIN']]
        vmax_list = [x[1] for x in self.pdict['GROUPVEL_WIN']]
        imeas = self.pdict['IMEAS']

        # allocate mpi jobs
        nsta = self.nsta 
        nsta_loc = self.nsta_loc

        # misfits
        ncomp = self.ncomp
        tstart = np.zeros((nsta_loc))
        tend = np.zeros((nsta_loc))
        win_chi = np.zeros((nsta_loc,ncomp,20))
        tr_chi = np.zeros((nsta_loc,ncomp))
        am_chi = np.zeros((nsta_loc,ncomp))

        # loop every station to do preprocessing
        for ir in range(nsta_loc):
            i = ir + self._istart
            for icomp in range(self.ncomp):
                name = self._get_station_code(i,icomp)
                obs_tr = SACTrace.read(f'{self.DATA_DIR}/{self.evtid}/{name}.sac')
                syn_tr = np.load(f"{out_dir}/{name}.sem.npy")[:,1]

                # obs info
                t0_obs = obs_tr.b 
                dt_obs = obs_tr.delta 
                npt_obs = obs_tr.npts 

                # new npts/t0 for interpolate
                t1_inp = t0_syn + (npt_syn - 1) * dt_syn
                npt1_inp = int((npt_obs - 1) * dt_obs / dt_inp)
                npt_cut = int((t1_inp - t0_inp) / dt_inp) + 1

                # interp obs/syn
                dat_inp1 = interpolate_syn(obs_tr.data,t0_obs,dt_obs,npt_obs,
                                        t0_obs + dt_inp,dt_inp,npt1_inp)
                syn_inp = interpolate_syn(syn_tr.data,t0_syn,dt_syn,npt_syn,
                                         t0_inp,dt_inp,npt_cut)

                # compute time derivative 
                if self.pdict['USE_EGF'] == False:
                    dat_inp1 = -dif1(dat_inp1,dt_inp)
                    if self.myrank == 0: print("CCFs => EGFs ...")
                
                
                # cut 
                dat_inp = interpolate_syn(dat_inp1,t0_obs + dt_inp,dt_inp,npt1_inp,
                                         t0_inp,dt_inp,npt_cut)

                # preprocess
                dat_inp = bandpass(dat_inp,dt_inp,freqmin,freqmax)
                syn_inp = bandpass(syn_inp,dt_inp,freqmin,freqmax)
        
                # normalize my amplitude of the window
                dist = obs_tr.dist
                win_b = np.floor((dist / vmax_list[ib] - self.Tmax[ib] / 2. + t0_inp) / dt_inp)
                win_e = np.floor((dist / vmin_list[ib] + self.Tmax[ib] / 2. - t0_inp) / dt_inp)
                win_b = int(max(win_b,0))
                win_e = int(min(win_e,len(dat_inp)-1))
                dat_inp *= np.max(np.abs(syn_inp[win_b:win_e])) / np.max(np.abs(dat_inp[win_b:win_e]))

                # compute time window
                tstart[ir] = dist / vmax_list[ib] - self.Tmax[ib] * 0.5 
                tend[ir] = dist / vmin_list[ib] + self.Tmax[ib] * 0.5 
                tstart[ir] = max(tstart[ir],t0_inp)
                tend[ir] = min(tend[ir],t0_obs+(npt_obs-1)*dt_obs)

                # compute misfits and adjoint source
                verbose = (self.myrank == 0) and (ir == 0) and (icomp == 0)
                tr_chi[ir,icomp],am_chi[ir,icomp],win_chi[ir,icomp,:],adjsrc =   \
                    measure_adj(t0_inp,dt_inp,npt_cut,t0_syn,dt_syn,npt_syn,
                                tstart[ir],tend[ir],imeas,self.Tmax[ib]*1.01,
                                self.Tmin[ib]*0.99,verbose,dat_inp,
                                syn_inp,)

                # save adjoint source
                data = np.zeros((npt_syn,2))
                data[:,0] = t0_syn + np.arange(npt_syn) * dt_syn
                data[:,1] = adjsrc
                name = self._get_station_code(i,icomp) + ".adj.sem.npy"
                np.save(f"{out_dir}/{bandname}/{name}",data)

                # save obs and syn as sac
                name = self._get_station_code(i,icomp) + ".sac"
                filename = f"{out_dir}/{bandname}/{name}.obs"
                obs_tr.delta = dt_inp
                obs_tr.b = t0_inp
                obs_tr.data = dat_inp
                obs_tr.write(filename)
                filename = f"{out_dir}/{bandname}/{name}.syn"
                obs_tr.data = syn_inp
                obs_tr.write(filename)
        
        # print info and save MEASUREMENTS file
        self._print_measure_info(bandname,tstart,tend,tr_chi,am_chi,win_chi)