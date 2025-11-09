from fwat.measure.FwatPreOP import FwatPreOP
import numpy as np 
from mpi4py import MPI

def _rf_adj_src(rf_obs,rf_syn,synr,synz,
                dt,tshift,f0,maxit):
    from .tele.deconit import deconit,myconvolve,nextpow2
    from .tele.deconit import gauss_filter,apply_gaussian

    def shift_data(a,dt,t0):
        a1 = np.fft.fft(a)
        om = 2 * np.pi * np.fft.fftfreq(len(a),dt)
        a1 = a1 * np.exp(-1j * om * t0)

        return np.fft.ifft(a1).real
    
    # get a filte
    nft = nextpow2(len(rf_obs)* 2)
    synr_filt = np.pad(synr,(0,nft - len(synr)),'constant')
    synz_filt = np.pad(synz,(0,nft - len(synz)),'constant')
    rf_o_filt = np.pad(rf_obs,(0,nft - len(rf_obs)),'constant')
    rf_s_filt = np.pad(rf_syn,(0,nft - len(rf_syn)),'constant')

    # gaussian filter
    gauss = gauss_filter(nft,dt,f0)
    synr_filt = apply_gaussian(synr_filt,gauss,dt)
    synz_filt = apply_gaussian(synz_filt,gauss,dt)

    # normalize
    zrf = deconit(
        synz_filt,synz_filt,dt,0.,
        f0,0,maxit
    )
    amp = np.max(abs(zrf))
    dobs_norm = rf_o_filt / amp
    dsyn_norm = rf_s_filt / amp

    # get negative shifted data, then reverse
    r_rev = shift_data(synr_filt,dt,-tshift)[::-1]
    z_rev = shift_data(synz_filt,dt,-tshift)[::-1]
    
    # adjoint source
    dR = dsyn_norm - dobs_norm
    adj_r_filt = deconit(dR,z_rev,dt,0,f0,1,maxit)
    tmp = myconvolve(-dR,r_rev)
    tmp1 = myconvolve(z_rev,z_rev)
    adj_z_filt = deconit(tmp,tmp1,dt,0,f0,1,maxit)

    # filter
    adj_r_filt = apply_gaussian(adj_r_filt,gauss,dt)
    adj_z_filt = apply_gaussian(adj_z_filt,gauss,dt)

    # copy to original length
    nt = len(rf_obs)
    adj_r = adj_r_filt[:nt]
    adj_z = adj_z_filt[:nt]

    return adj_r,adj_z,amp

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
        from .tele.deconit import deconit
        from .utils import bandpass,taper_window

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
            b=-self._tshift
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
            freqmin = 1. / self.Tmax[0]
            freqmax = 1. / self.Tmin[0]
            for ib in range(len(self._f0)):
                f0 = self._f0[ib]
                name = self._get_rf_code(i,ib)

                # bandpass
                idx_r = self.components.index('R')
                idx_z = self.components.index('Z')
                R = bandpass (
                    syn_data[idx_r,:],dt_syn,
                    freqmin,freqmax
                )
                Z = bandpass (
                    syn_data[idx_z,:],dt_syn,
                    freqmin,freqmax
                )

                # time window 
                lpt,rpt,taper = taper_window(
                    0,dt_syn,npt_syn,
                    self.t_ref[i] - self.t_inj - win_tb,
                    self.t_ref[i] - self.t_inj + win_te
                )
                r_win = R[lpt:rpt] * taper
                z_win = Z[lpt:rpt] * taper

                # compute rf
                rf = deconit(
                    r_win,z_win,
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
        from obspy.io.sac import SACTrace
        from fwat.measure.utils import interpolate_syn,bandpass,taper_window
        from fwat.measure.tele.deconit import deconit
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
        out_dir = f"{self.syndir}/OUTPUT_FILES/"

        # get vars
        t0_syn = self.t0_syn
        dt_syn = self.dt_syn
        npt_syn = self.npt_syn

        # get time window 
        win_tb,win_te = self.pdict['TIME_WINDOW']
                
        # allocate global arrays 
        nsta_loc = self.nsta_loc

        # pre-allocate misfits
        tstart = np.zeros((nsta_loc))
        tend = np.zeros((nsta_loc))
        win_chi = np.zeros((nsta_loc,1,20))
        tr_chi = np.zeros((nsta_loc,1))
        am_chi = np.zeros((nsta_loc,1))
        tstart[:] = max(- win_tb,t0_syn)
        tend[:] = min(win_te,t0_syn + dt_syn * npt_syn)

        # loop each station
        for ir in range(nsta_loc):
            i = ir + self._istart
        
            # read synthetic data
            syn_data = np.zeros((2,npt_syn))
            for ic in range(self.ncomp):
                name = self._get_station_code(i,ic)
                syn_data[ic,:] = self.seismogram[f"{out_dir}/{name}.sem.npy"][:,1]
                #syn_data[ic,:] = np.load(f"{out_dir}/{name}.sem.npy")[:,1]

                # filter
                syn_data[ic,:] = bandpass(syn_data[ic,:],dt_syn,freqmin,freqmax)
            
            # get time window 
            taper = np.zeros((npt_syn))
            t1 = self.t_ref[i] - self.t_inj - win_tb 
            t2 = self.t_ref[i] - self.t_inj + win_te
            lpt,rpt,taper0 = taper_window(0,dt_syn,npt_syn,t1,t2)
            taper[lpt:rpt] = taper0 * 1.

            # fetch windowed data and compute receiver function
            idx_r = self.components.index('R')
            idx_z = self.components.index('Z')
            r_win = syn_data[idx_r,lpt:rpt] * taper[lpt:rpt]
            z_win = syn_data[idx_z,lpt:rpt] * taper[lpt:rpt]
            rf_syn = deconit(
                r_win,
                z_win,
                dt_syn,self._tshift,
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

            # compute adjoint source
            tp = self.t_ref[i] - self.t_inj
            adj_r_win,adj_z_win,amp  = _rf_adj_src(
                rf_obs,rf_syn,
                r_win,z_win,
                dt_syn,self._tshift,self._f0[ib],
                self._maxit
            )

            # taper  RF
            rf_obs *= taper / amp 
            rf_syn *= taper / amp 

            # misfit 
            chi = 0.5 * trapezoid((rf_obs - rf_syn)**2, dx=dt_syn)
            tr_chi[ir] = chi 
            am_chi[ir] = chi
            win_chi[ir,0,13-1] = 0.5 * sum( rf_obs**2 )
            win_chi[ir,0,14-1] = 0.5 * sum( rf_syn**2 )
            win_chi[ir,0,15-1] = chi
            win_chi[ir,0,20-1] = npt_syn * dt_syn

            # filter and taper adj source
            adj_z = np.zeros_like(npt_syn)
            adj_r = np.zeros_like(npt_syn)
            adj_z[lpt:rpt] = bandpass(adj_z_win,dt_syn,freqmin,freqmax) * taper0
            adj_r[lpt:rpt] = bandpass(adj_r_win,dt_syn,freqmin,freqmax) * taper0

            # save adjoint source
            data = np.zeros((npt_syn,2))
            data[:,0] = t0_syn + np.arange(npt_syn) * dt_syn
            data[:,1] = adj_r
            outname = f"{out_dir}/{bandname}/{self.netwk[i]}.{self.stnm[i]}.{self.chcode}T.adj.sem.npy"
            self.seismogram_adj[outname] = data.copy()
            #np.save(outname,data)
            data[:,1] = adj_r
            outname = f"{out_dir}/{bandname}/{self.netwk[i]}.{self.stnm[i]}.{self.chcode}R.adj.sem.npy"
            self.seismogram_adj[outname] = data.copy()
            #np.save(outname,data)

            # save obs and syn
            name = f"{self.netwk[i]}.{self.stnm[i]}.{self.chcode}R.rf.sac"
            tr = SACTrace(
                evla=self.evla,evlo=self.evlo,
                evdp=self.evdp,stla=self.stla[i],
                stlo=self.stlo[i],stel=0,lcalda=1,
                delta = dt_syn,b=-self._tshift,
                knetwk = self.netwk[i],
                kstnm = self.stnm[i],
                kcmpnm = self.chcode + 'R'
            )
            tr.data = rf_obs
            #tr.write()
            self.seismo_sac[f"{out_dir}/{bandname}/{name}.obs"] = tr.copy()
            tr.data = rf_syn
            #tr.write(f"{out_dir}/{bandname}/{name}.syn")
            self.seismo_sac[f"{out_dir}/{bandname}/{name}.syn"] = tr.copy()
        
        # save measurement files
        self._print_measure_info(bandname,tstart,tend,tr_chi,am_chi,win_chi)