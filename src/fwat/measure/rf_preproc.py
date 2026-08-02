import numpy as np 
from mpi4py import MPI

from fwat.measure.FwatPreOP import FwatPreOP
from ..adjoint.MeasureStats import MeasureStats
from .tele.deconit import deconit,myconvolve,nextpow2
from .tele.deconit import gauss_filter,apply_gaussian
from .utils import taper_window

def _rf_adj_src_freq(rf_obs,rf_syn,synr,synz,
                dt,tshift,f0):
    
    # append to avoid aliasing
    nfft = nextpow2(len(rf_obs)*2)
    rf_o_flt = np.zeros((nfft))
    rf_o_flt[:len(rf_obs)] = rf_obs.copy()
    rf_s_flt = np.zeros((nfft))
    rf_s_flt[:len(rf_syn)] = rf_syn.copy()
    synr_flt = np.zeros((nfft))
    synr_flt[:len(synr)] = synr.copy()
    synz_flt = np.zeros((nfft))
    synz_flt[:len(synz)] = synz.copy()

    # gauss filter
    gauss = gauss_filter(nfft,dt,f0)

    # go to frequency domain
    RFo = np.fft.fft(rf_o_flt) * dt 
    RFs = np.fft.fft(rf_s_flt) * dt
    SYNr = np.fft.fft(synr_flt) * dt
    SYNz = np.fft.fft(synz_flt) * dt
    dR = RFs - RFo

    # shift factor 
    om = 2 * np.pi * np.fft.fftfreq(nfft,dt)
    shift_factor = np.exp(1j * om * tshift)

    # compute H componet
    deno = abs(SYNz * SYNz.conj())
    eps = np.max(np.abs(deno)) * 1.0e-4
    deno[deno < eps] = eps
    adj_h_fq = dR * SYNz / deno * shift_factor * gauss

    # Z component
    deno1 = np.abs(SYNz * SYNz.conj())**2 
    eps = np.max(np.abs(deno1)) * 1.0e-4
    deno1[deno1 < eps] = eps
    adj_z_fq = -dR * SYNr.conj() * SYNz**2 / deno1 * shift_factor * gauss

    # go back to time domain
    adj_r = np.fft.ifft(adj_h_fq).real[:len(rf_obs)] / dt
    adj_z = np.fft.ifft(adj_z_fq).real[:len(rf_obs)] / dt 

    return adj_r,adj_z

class RF_PreOP(FwatPreOP):
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
        
        # get rf parameters
        self._f0 = self.pdict['GAUSS_F0']
        self._maxit = self.pdict['MAXIT']
        self._minerr = self.pdict['MINERR']
        self._tshift = self.pdict['TSHIFT']

        # only one frequency band is required
        if len(self.Tmax) > 0:
            if self.myrank == 0:
                print("only One frequency band is permitted!")
                print(f"previous band Tmax = {self.Tmax}")
                print(f"previous band Tmin = {self.Tmin}")
            self.Tmax = [np.max(self.Tmax)]
            self.Tmin = [np.min(self.Tmin)]

            if self.myrank == 0:
                print(f"keep max Tmax = {self.Tmax}")
                print(f"keep min Tmin = {self.Tmin}")

    def _get_bandname(self,ib:int):
        return "F%2.1f" %(self._f0[ib])
    
    def _get_rf_code(self,i:int,if0:int):
        bandname = self._get_bandname(if0)
        code = f"{self.netwk[i]}.{self.stnm[i]}.{self.chcode}R.{bandname}"

        return code 
    
    def _sanity_check(self):
        super()._sanity_check()

        # make sure adjsrc_type = 2
        assert self.pdict['ADJSRC_TYPE'] == '2', f"For receiver function, the ADJSRC_TYPE should be 2"

    def save_forward(self):
        import os 
        from obspy.io.sac import SACTrace

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
            delta = dt_syn,
            b=-self._tshift,
            isynth = 'irldta'
        )

        # time window 
        win_tb,win_te = self.pdict['TIME_WINDOW']

        outdir = f"{self.DATA_DIR}/{evtid}"
        if self.myrank == 0:
            print("Synthetic Observations ...")
            os.makedirs(outdir,exist_ok=True)
        MPI.COMM_WORLD.Barrier()
        
        # loop every station to save sac
        for ir in range(self.nsta_loc):
            i = self._istart + ir 

            # add tags
            tr.knetwk = self.netwk[i]
            tr.kstnm = self.stnm[i]
            tr.stla = self.stla[i]
            tr.stlo = self.stlo[i]
            tr.kcmpnm = f"{self.chcode}R"

            # get syn_data
            syn_data = np.zeros((2,npt_syn))
            for ic in range(ncomp):
                # load synthetics from npy
                code = self._get_station_code(i,ic)
                filename = f"{self.syndir}/OUTPUT_FILES/{code}.sem.npy"
                data =  np.load(filename)
                
                # copy to syn_data
                syn_data[ic,:] = data[:,1] * 1.

            # compute rf
            for ib in range(len(self._f0)):
                f0 = self._f0[ib]
                name = self._get_rf_code(i,ib)

                # get window used for RF
                pad = max(5.0,2.0 / self._f0[ib])
                tb = self.t_ref[i] - self.t_inj - win_tb - pad
                te = self.t_ref[i] - self.t_inj + win_te + pad
                lpt,rpt,taper0 = taper_window(0,dt_syn,npt_syn,tb,te)
                taper_p = np.zeros((npt_syn))
                taper_p[lpt:rpt] = taper0

                # bandpass
                idx_r = self.components.index('R')
                idx_z = self.components.index('Z')
                R = syn_data[idx_r,:] * taper_p
                Z = syn_data[idx_z,:] * taper_p

                # compute rf
                rf = deconit(
                    R,Z,
                    dt_syn,
                    self._tshift,
                    f0,0,self._maxit
                )

                # save to sac
                tr.data = rf * 1.
                tr.user1 = f0
                name = self._get_rf_code(i,ib) + ".rf.sac"
                filename = f"{outdir}/{name}"
                tr.write(filename)

    def cal_adj_source(self,ib:int):
        """
        Calculate adjoint source for RF measurement.
        
        Parameters
        -----------
        ib: int
            index of frequency band
        """

        from obspy.io.sac import SACTrace
        from fwat.measure.utils import interpolate_syn,bandpass,taper_window
        from scipy.integrate import trapezoid
        import os 
        bandname = self._get_bandname(ib)
        if self.myrank == 0:
            print(f"preprocessing for band {bandname} ...")
            os.makedirs(f"{self.syndir}/OUTPUT_FILES/{bandname}",exist_ok=True)
        MPI.COMM_WORLD.Barrier()

        # get frequency band
        freqmin = 1. / self.Tmax[0]
        freqmax = 1. / self.Tmin[0]
        out_dir = f"{self.syndir}/OUTPUT_FILES"

        # get vars
        t0_syn = self.t0_syn
        dt_syn = self.dt_syn
        npt_syn = self.npt_syn

        # get time window 
        win_tb,win_te = self.pdict['TIME_WINDOW']
                
        # allocate global arrays 
        nsta_loc = self.nsta_loc

        # tstart/tend for each station
        stats_list = [MeasureStats(adj_type=self.adjsrc_type) for _ in range(nsta_loc)]
        tstart =  max(- win_tb,t0_syn)
        tend = min(win_te,t0_syn + dt_syn * npt_syn)

        # gauss filter
        gauss = gauss_filter(npt_syn,dt_syn,self._f0[ib])

        # save seismo_win headers
        self.seismo_win['dt'] = dt_syn
        self.seismo_win['t0'] = -self._tshift
        self.seismo_win['npts'] = npt_syn

        # get average amplitude of rf_obs
        amp = 0.
        for ir in range(nsta_loc):
            i = ir + self._istart
            rfname = self._get_rf_code(i,ib)
            obs_tr = SACTrace.read(f"{self.DATA_DIR}/{self.evtid}/{rfname}.rf.sac")
            amp += np.max(np.abs(obs_tr.data))

        # sync
        amp = MPI.COMM_WORLD.allreduce(amp,op=MPI.SUM)
        amp /= self.nsta
        
        # loop each station
        for ir in range(nsta_loc):
            i = ir + self._istart

            # get window used for RF
            pad = max(5.0,2.0 / self._f0[ib])
            tb = self.t_ref[i] - self.t_inj - win_tb - pad
            te = self.t_ref[i] - self.t_inj + win_te + pad
            lpt,rpt,taper0 = taper_window(0,dt_syn,npt_syn,tb,te)
            taper_p = np.zeros((npt_syn))
            taper_p[lpt:rpt] = taper0

            # read synthetic data
            syn_data = np.zeros((2,npt_syn))
            for ic in range(self.ncomp):
                name = self._get_station_code(i,ic)
                syn_data[ic,:] = self.seismogram[f"{out_dir}/{name}.sem.npy"][:,1]
                syn_data[ic,:] *= taper_p # taper synthetic data
            
            # compute rf 
            idx_r = self.components.index('R')
            idx_z = self.components.index('Z')
            rf_syn = deconit(
                syn_data[idx_r,:],
                syn_data[idx_z,:],
                dt_syn,
                self._tshift,
                self._f0[ib],0,
                self._maxit
            )

            # read obs data
            rfname = self._get_rf_code(i,ib)
            obs_tr = SACTrace.read(f"{self.DATA_DIR}/{self.evtid}/{rfname}.rf.sac")
            t0_obs = obs_tr.b 
            dt_obs = obs_tr.delta 
            npt_obs = obs_tr.npts 

            # interpolate rf_obs to time window
            rf_obs = interpolate_syn(
                obs_tr.data,
                t0_obs,dt_obs,npt_obs,
                -self._tshift,dt_syn,len(rf_syn)
            )

            # get time window 
            lpt,rpt,taper1 = taper_window(
                -self._tshift,dt_syn,npt_syn,
                -win_tb,win_te
            )
            taper = rf_syn * 0
            taper[lpt:rpt] = taper1 

            # taper rf 
            rf_obs *= taper
            rf_syn *= taper

            # compute adjoint source
            adj_r,adj_z  = \
            _rf_adj_src_freq(
                rf_obs*taper,rf_syn*taper,
                syn_data[idx_r,:],
                syn_data[idx_z,:],
                dt_syn,
                self._tshift,
                self._f0[ib]
            )

            # taper again and apply amp
            adj_r *= taper_p / amp**2 
            adj_z *= taper_p / amp**2

            # misfit 
            chi = 0.5 * trapezoid((rf_obs - rf_syn)**2 / amp**2, dx=dt_syn)
            stats = MeasureStats(
                adj_type=self.adjsrc_type,
                misfit=chi,
                tstart=tstart,
                tend=tend,
                code=f"{self.netwk[i]}.{self.stnm[i]}.BXR",
                tr_chi=chi,
                am_chi=chi
            )
            stats_list[ir] = stats

            # save adjoint source
            data = np.zeros((npt_syn,2))
            data[:,0] = t0_syn + np.arange(npt_syn) * dt_syn
            data[:,1] = adj_z
            outname = f"{out_dir}/{bandname}/{self.netwk[i]}.{self.stnm[i]}.{self.chcode}Z.adj.sem.npy"
            self.seismogram_adj[outname] = data.copy()
            #np.save(outname,data)
            data[:,1] = adj_r
            outname = f"{out_dir}/{bandname}/{self.netwk[i]}.{self.stnm[i]}.{self.chcode}R.adj.sem.npy"
            self.seismogram_adj[outname] = data.copy()
            #np.save(outname,data)

            # save obs and syn
            name = f"{self.netwk[i]}.{self.stnm[i]}.{self.chcode}R.rf.dat"
            self.seismo_win[f"{out_dir}/{bandname}/{name}.obs"] = rf_obs.copy()
            self.seismo_win[f"{out_dir}/{bandname}/{name}.syn"] = rf_syn.copy()
        
        # save measurement files
        self._print_measure_info(bandname,stats_list)