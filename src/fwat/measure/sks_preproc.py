from .FwatPreOP import FwatPreOP
import numpy as np 
from mpi4py import MPI
    
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

    def save_forward(self):
        import os 
        from obspy.io.sac import SACTrace
        from .utils import taper_window,bandpass

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

        # init syn_data 
        syn_data = np.zeros((2,npt_syn))
        
        # loop every station to save sac
        for ir in range(self.nsta_loc):
            i = self._istart + ir 

            # time window
            win_tb,win_te = self.pdict['TIME_WINDOW']
            taper = np.zeros((npt_syn))
            tstart = self.t_ref[i] - self.t_inj - win_tb
            tend = self.t_ref[i] - self.t_inj + win_te
            lpt,rpt,taper0 = taper_window(0,dt_syn,npt_syn,tstart,tend)
            taper[lpt:rpt] = taper0 * 1.

            for ic in range(ncomp):
                ch = components[ic]

                # load synthetics from npy
                code = self._get_station_code(i,ic)
                filename = f"{self.syndir}/OUTPUT_FILES/{code}.sem.npy"
                data = self.seismogram[filename]
                #data =  np.load(filename)

                tr.b =  self.t_inj - self.t_ref[i]
                tr.data = data[:,1] * 1.

                # syn data
                syn_data[ic,:] = data[:,1] * taper

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
        from .utils import interpolate_syn
        from .utils import bandpass,taper_window
        import os 

        bandname = self._get_bandname(ib)
        if self.myrank == 0:
            print(f"preprocessing for band {bandname} ...")
            os.makedirs(f"{self.syndir}/OUTPUT_FILES/{bandname}",exist_ok=True)
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
        tstart = np.zeros((nsta_loc))
        tend = np.zeros((nsta_loc))
        win_chi = np.zeros((nsta_loc,1,20))
        tr_chi = np.zeros((nsta_loc,1))
        am_chi = np.zeros((nsta_loc,1))

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

        # loop each station
        for ir in range(nsta_loc):
            i = ir + self._istart
            obs_data = np.zeros((2,npt_syn))
            syn_data = np.zeros((2,npt_syn))

            # time window
            taper = np.zeros((npt_syn))
            tstart[ir] = self.t_ref[i] - t_inj - win_tb
            tend[ir] = self.t_ref[i] - t_inj + win_te
            lpt,rpt,taper0 = taper_window(0,dt_syn,npt_syn,tstart[ir],tend[ir])
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
                
                # save SAC
                tr = SACTrace(
                    evla=self.evla,evlo=self.evlo,
                    evdp=self.evdp,stla=self.stla[i],
                    stlo=self.stlo[i],stel=0,lcalda=1,
                    delta = dt_syn,b=t0_syn,
                    t0 = self.t_ref[i] - t_inj,
                    knetwk = self.netwk[i],
                    kstnm = self.stnm[i],
                    kcmpnm = self.chcode + self.components[ic]
                )
                
                # save syn
                tr.data = syn_data[ic,:]
                self.seismo_sac[f"{out_dir}/{bandname}/{name}.sac.syn"] = tr.copy()
                #tr.write(f"{out_dir}/{bandname}/{name}.sac.syn")

                # save obs
                if cal_obs_si:
                    tr.data = obs_data[ic,:]
                    self.seismo_sac[f"{out_dir}/{bandname}/{name}.sac.obs"] = tr.copy()
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
            tr_chi[ir,0] = L 
            am_chi[ir,0] = L
            win_chi[ir,0,6] = SI_syn
            win_chi[ir,0,7] = SI_obs

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
        
        # save measurement files
        self._print_measure_info(bandname,tstart,tend,tr_chi,am_chi,win_chi)