import numpy as np 
from mpi4py import MPI

class FwatPreOP:
    def __init__(self,measure_type:str,iter:int,evtid:str,run_opt:int):
        # import packages
        from utils import read_params,alloc_mpi_jobs
        from utils import get_source_loc
        import h5py 
        comm = MPI.COMM_WORLD
        self.myrank = comm.Get_rank()
        self.nprocs = comm.Get_size()

        # private vars
        self.SRC_REC = "./src_rec/"
        self.SOLVER = "./solver/"
        self.DATA_DIR = "./fwat_data/"
        self.MISFIT = "./misfits/"
        
        # copy in input parameters
        assert(measure_type in ['sks','tele','noise'])
        self.meatype = measure_type
        self.iter = iter 
        self.evtid = evtid
        self.run_opt = run_opt

        # path to simulation directory
        mdir = "M%02d" %(iter)
        if self.run_opt == 2:
            mdir = mdir + ".ls"
        self.mod = mdir
        self.syndir:str = f"{self.SOLVER}/{mdir}/{evtid}/"

        # read FWAT params below
        # frequency band
        pdict = read_params(f"{self.syndir}/DATA/FWAT.PAR.yaml")['measure'][self.meatype]
        self.Tmin = [x[0] for x in pdict['FILTER_BANDS']]
        self.Tmax = [x[1] for x in pdict['FILTER_BANDS']]

        # components
        self.ncomp = len(pdict['COMPS'])
        self.components = pdict['COMPS']
        self.components = sorted(self.components)

        # channel code
        self.chcode = pdict['CH_CODE']

        # backup pdict for furthur usage
        self.pdict = pdict
        ### end ##########

        # read sourcer/receiver locations below
        # ---------------------
        # source loc
        sourcefile = f"{self.SRC_REC}/sources.dat.{measure_type}"
        self.evla,self.evlo,self.evdp = get_source_loc(evtid,sourcefile)

        # station info
        # read station coordinates
        stationfile = f'{self.SRC_REC}/STATIONS_{evtid}_globe'
        statxt = np.loadtxt(stationfile,dtype=str,ndmin=2)
        self.nsta = statxt.shape[0]
        self.netwk = statxt[:,1]
        self.stnm = statxt[:,0]
        self.stla = statxt[:,2].astype(float)
        self.stlo = statxt[:,3].astype(float)
        self.bazd = np.zeros((self.nsta))

        # read simulation info dt,t0,npts
        fio = h5py.File(f"{self.syndir}/OUTPUT_FILES/seismograms.h5","r")
        t = fio[list(fio.keys())[0]][:,0]
        fio.close()
        self.t0_syn = t[0]
        self.dt_syn = t[1] - t[0]
        self.npt_syn = len(t)

        # read tref and t_inj for teleseismic events
        if measure_type in ['tele','sks']:
            from tele.tele import get_injection_time
            from tele.tele import compute_ak135_time
            self.t_inj = get_injection_time(self.evtid)

            if "_" in self.evtid:
                phase = self.evtid.split('_')[0]
            else:
                if measure_type == "tele":
                    phase = 'P'
                elif measure_type == 'sks':
                    phase = 'sks'
            self.t_ref = compute_ak135_time(
                self.evla,self.evlo,self.evdp,
                statxt,phase
            )
        
        # allocate jobs for each proc 
        istart,iend = alloc_mpi_jobs(self.nsta,self.nprocs,self.myrank)
        self.nsta_loc = iend - istart + 1
        self._istart = istart 

        # sanity check 
        self._sanity_check()
    
    def _sanity_check(self):
        pass

    def _get_station_code(self,i:int,ic:int):
        """
        get station code network.staname.channel 

        Parameters
        ----------------
        i: int 
            id of station in self.netwk
        ic: int
            id of channel in self.components

        Returns
        --------------
        code: str
            station code network.staname.channel 
        """
        ch = self.components[ic]
        name = f"{self.netwk[i]}.{self.stnm[i]}.{self.chcode}{ch}"

        return name 

    def _rotate_XYZ_to_ZNE(self):
        from cube2sph_rotate import rotate_seismo_fwd

        # set parameters
        fn_matrix = f"{self.SRC_REC}/rot_{self.evtid}"
        from_dir = f"{self.syndir}/OUTPUT_FILES/"
        to_dir = f"{self.syndir}/OUTPUT_FILES/"
        from_template='${nt}.${sta}.BX${comp}.semd'
        to_template='${nt}.${sta}.BX${comp}.sem.npy'

        # rotate seismograms from XYZ to ZNE
        rotate_seismo_fwd(fn_matrix,from_dir,to_dir,from_template,to_template)
    
    def _rotate_ZNE_to_XYZ(self):
        from cube2sph_rotate import rotate_seismo_adj

        # set parameters
        fn_matrix = f"{self.SRC_REC}/rot_{self.evtid}"
        from_dir = f"{self.syndir}/SEM/"
        to_dir = f"{self.syndir}/SEM/"
        from_template='${nt}.${sta}.BX${comp}.adj.sem.npy'
        to_template='${nt}.${sta}.BX${comp}.adj'

        # rotate seismograms from XYZ to ZNE
        rotate_seismo_adj(fn_matrix,from_dir,to_dir,from_template,to_template)
        
    def save_forward(self):
        import os 
        from obspy.io.sac import SACTrace
        from utils import cumtrapz1

        # get some vars
        evtid = self.evtid 
        ncomp = self.ncomp
        components = self.components
        dt_syn = self.dt_syn
        npt_syn = self.npt_syn
        t0_syn = self.t0_syn
        myrank = self.myrank
        nprocs = self.nprocs

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
        stf = None
        if self.meatype == "tele":
            from scipy.signal import convolve
            stf = np.zeros((ncomp,npt_syn))

            # sum stf from all bands
            for ib in range(self.Tmax):
                bandname = 'T%03g_T%03g' %(self.Tmin[ib],self.Tmax[ib])
                for ic in range(self.ncomp):
                    ch = self.components[ic]
                    tr = SACTrace.read(f'src_rec/stf_{ch}.sac.{bandname}_{evtid}')
                    stf[ic,:] += tr.data 
        
        # loop every station to save sac
        for ir in range(self.nsta_loc):
            i = self._istart + ir 
            for ic in range(ncomp):
                ch = components[ic]

                # load synthetics from ascii
                code = self._get_station_code(i,ic)
                filename = f"{self.syndir}/OUTPUT_FILES/{code}.sem.npy"
                data =  np.load(filename)

                # special handling for each type
                if self.meatype == "tele":
                    # convolve
                    tr.data = convolve(data[:,1],stf[ic,:],'same') * dt_syn    
                    tr.b =  self.t_inj - self.t_ref[i]
                elif self.meatype == "sks":
                    tr.b =  self.t_inj - self.t_ref[i]
                    tr.data = data[:,1] * 1.
                elif self.meatype == "noise":
                    if self.pdict['USE_EGF'] == False:
                        if i == myrank and myrank == 0 and ic == 0:
                            print("EGF => CCF ...")
                        tr.data = -cumtrapz1(data[:,1],dt_syn)
                    tr.b = t0_syn
                else:
                    print(f"{self.meatype} is not implemeted!")
                    exit(1)

                # channel and others
                tr.kcmpnm = f"{self.chcode}{ch}"
                tr.knetwk = self.netwk[i]
                tr.kstnm = self.stnm[i]
                tr.stla = self.stla[i]
                tr.stlo = self.stlo[i]

                # save to sac
                filename = f"{outdir}/{code}.sac"
                tr.write(filename)

    def _get_obs_info(self):
        from obspy.io.sac import SACTrace
        name = self._get_station_code(0,0) + ".sac"
        hd = SACTrace.read(f"{self.DATA_DIR}/{self.evtid}/{name}",headonly=True)
        self.t0_obs = hd.b 
        self.dt_obs = hd.delta
        self.npt_obs = hd.npts 

    def _print_measure_info(self,bandname,tstart,tend,tr_chi,am_chi,window_chi):
        ncomp = tr_chi.shape[1]
        nsta_loc = self.nsta_loc
        imeas = self.pdict['IMEAS']

        # sync
        MPI.COMM_WORLD.Barrier()

        # loop each proc to print info on the screen
        # and save files
        for irank in range(self.nprocs):
            if irank == self.myrank:

                # open outfile
                outfile = f"{self.MISFIT}/{self.mod}_{self.evtid}_{bandname}_{self.meatype}_window_chi"
                if irank == 0:
                    fio = open(outfile,"w")
                else:
                    fio = open(outfile,"a")

                for ir in range(nsta_loc):
                    i = ir + self._istart 
                    for ic in range(ncomp):
                        name = self._get_station_code(i,ic)
                        print(f'{name}')
                        print("Measurement window No.  1")
                        print("start and end time of window: %f %f" %(tstart[ir],tend[ir]) )
                        print(f"adjoint source and chi value for imeas = {imeas}")
                        print("%e" %(window_chi[ir,ic,6]))
                        print("tr_chi = %e am_chi = %e" %(tr_chi[ir,ic],am_chi[ir,ic]))
                        print("")

                        ch = self.components[ic]
                        fio.write(f"{self.evtid} {self.stnm[i]} {self.netwk[i]} {self.chcode}{ch} 1 {imeas} ")
                        fio.write("%g %g " %(tstart[ir],tend[ir]))
                        for j in range(20):
                            fio.write("%g " %(window_chi[ir,ic,j]))
                        fio.write("%g %g 0. 0.\n" %(tr_chi[ir,ic],am_chi[ir,ic]))
                
                # close output file
                fio.close()

            # barrier
            MPI.COMM_WORLD.Barrier()


    def _cal_adjsrc_noise(self,ib:int,bandname:str):
        from obspy.io.sac import SACTrace
        from utils import interpolate_syn,dif1
        from utils import bandpass
        from measure import measure_adj

        freqmin = 1. / self.Tmax[ib]
        freqmax = 1. / self.Tmin[ib]
        out_dir = f"{self.syndir}/OUTPUT_FILES/"

        # new dt/nt for interpolate
        dt_inp = 0.01
        t0_inp = -10.

        # get vars
        t0_obs = self.t0_obs
        t0_syn = self.t0_syn
        dt_obs = self.dt_obs
        dt_syn = self.dt_syn
        npt_syn = self.npt_syn
        npt_obs = self.npt_obs

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
    
    def _cal_adjsrc_tele(self,ib:int,bandname:str):
        from obspy.io.sac import SACTrace
        from utils import interpolate_syn
        from utils import bandpass
        from measure import measure_adj
        from scipy.signal import convolve,correlate
        from tele.tele import compute_stf,get_average_amplitude
        import os 

        freqmin = 1. / self.Tmax[ib]
        freqmax = 1. / self.Tmin[ib]
        out_dir = f"{self.syndir}/OUTPUT_FILES/"

        # get vars
        t0_obs = self.t0_obs
        t0_syn = self.t0_syn
        dt_obs = self.dt_obs
        dt_syn = self.dt_syn
        npt_syn = self.npt_syn
        npt_obs = self.npt_obs
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
        glob_obs = np.zeros((nsta,ncomp,npt_syn))
        glob_syn = np.zeros((nsta,ncomp,npt_syn))

        # loop each station
        for ir in range(nsta_loc):
            i = ir + self._istart
            for ic in range(ncomp):
                name = self._get_station_code(i,ic)

                # read data
                syn_data = np.load(f"{out_dir}/{name}.sem.npy")[:,1]
                obs_tr = SACTrace.read(f"{self.DATA_DIR}/{self.evtid}/{name}.sac")
                obs_data = obs_tr.data

                # data processing
                obs_data = bandpass(obs_data,dt_obs,freqmin,freqmax)
                syn_data = bandpass(syn_data,dt_syn,freqmin,freqmax)

                # interpolate the obs/syn data to the same sampling of syn data
                t0_inp = self.t_ref[i] - win_tb
                u1 = interpolate_syn(obs_data,t0_obs + self.t_ref[i],dt_obs,npt_obs,t0_inp,dt_syn,npt2)
                w1 = interpolate_syn(syn_data,t_inj,dt_syn,npt_syn,t0_inp,dt_syn,npt2)
                u = interpolate_syn(u1,t0_inp,dt_syn,npt2,t_inj,dt_syn,npt_syn)
                w = interpolate_syn(w1,t0_inp,dt_syn,npt2,t_inj,dt_syn,npt_syn)     

                # save to  global array   
                glob_obs[i,ic,:] = u
                glob_syn[i,ic,:] = w 
        
        # all reduce
        comm = MPI.COMM_WORLD
        tmp = comm.allreduce(glob_obs); glob_obs = tmp * 1.
        tmp = comm.allreduce(glob_syn); glob_syn = tmp * 1.

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
                glob_syn[i,ic,:] = dt_syn * convolve(glob_syn[i,ic,:],stf[ic,:],'same')

                # time window
                tstart[ir] = self.t_ref[i] - self.t_inj - win_tb
                tend[ir] = self.t_ref[i] - self.t_inj + win_te

                # verboase
                verbose = (self.myrank == 0) and (ir == 0) and (ic == 0)
                
                tr_chi[ir,ic],am_chi[ir,ic],win_chi[ir,ic,:],adjsrc =  \
                    measure_adj(t0_syn,dt_syn,npt_syn,
                                t0_syn,dt_syn,npt_syn,
                                tstart[ir],tend[ir],2,
                                1.01*self.Tmax[ib],0.99*self.Tmin[ib],verbose,
                                glob_obs[i,ic,:],glob_syn[i,ic,:])

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
        tr_chi *= dt_syn / avgamp**2 
        am_chi *= dt_syn / avgamp**2
        win_chi[:,:,14] *= dt_syn / avgamp**2

        # print info and save MEASUREMENTS file
        self._print_measure_info(bandname,tstart,tend,tr_chi,am_chi,win_chi)
    
    def _cal_adjsrc_sks(self,ib:int,bandname:str):
        from obspy.io.sac import SACTrace
        from utils import interpolate_syn
        from utils import bandpass,dif1,taper_window

        freqmin = 1. / self.Tmax[ib]
        freqmax = 1. / self.Tmin[ib]
        out_dir = f"{self.syndir}/OUTPUT_FILES/"

        # get vars
        t0_obs = self.t0_obs
        t0_syn = self.t0_syn
        dt_obs = self.dt_obs
        dt_syn = self.dt_syn
        npt_syn = self.npt_syn
        npt_obs = self.npt_obs
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

        # loop each station
        for ir in range(nsta_loc):
            i = ir + self._istart
            obs_data = np.zeros((2,npt_syn))
            syn_data = np.zeros((2,npt_syn))

            # time window
            taper = np.zeros((npt_syn))
            tstart[ir] = self.t_ref[i] - t_inj - win_tb
            tend[ir] = self.t_ref[i] - t_inj + win_te
            lpt,rpt,taper0 = taper_window(0,dt_syn,tstart[ir],tend[ir])
            taper[lpt:rpt] = taper0 * 1.

            for ic in range(self.ncomp):
                name = self._get_station_code(i,ic)

                # read data
                sdata = np.load(f"{out_dir}/{name}.sem.npy")[:,1]
                obs_tr = SACTrace.read(f"{self.DATA_DIR}/{self.evtid}/{name}.sac")
                
                # interpoate obs data to same series of synthetics
                odata = interpolate_syn(obs_tr.data,t0_obs + self.t_ref[i],dt_obs,npt_obs,t_inj,dt_syn,npt_syn)

                # filter
                odata = bandpass(odata,dt_syn,freqmin,freqmax)
                sdata = bandpass(sdata,dt_syn,freqmin,freqmax)

                # save to global array   
                obs_data[ic,:] = odata * taper
                syn_data[ic,:] = sdata * taper 

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
                tr.data = obs_data[ic,:]
                tr.write(f"{out_dir}/{bandname}/{name}.sac.obs")
                tr.data = syn_data[ic,:]
                tr.write(f"{out_dir}/{bandname}/{name}.sac.syn")

            # compute normalization factor
            ic_r = self.components.index("R")
            ic_t = self.components.index("T")
            dRobs = dif1(obs_data[ic_r,:],dt_syn)
            dRsyn = dif1(syn_data[ic_r,:],dt_syn)
            norm_obs = np.sum(dRobs**2)
            norm_syn = np.sum(dRsyn**2)
            norm_obs = 1. / (norm_obs + 1.0e-30)
            norm_syn = 1. / (norm_syn + 1.0e-30)

            # compute si 
            Tobs = obs_data[ic_t,:]
            Tsyn = syn_data[ic_t,:]
            dTcomp = dif1(Tsyn,dt_syn)
            SI_obs = -2. * np.sum(dRobs * Tobs) * norm_obs 
            SI_syn = -2. * np.sum(dRsyn * Tsyn) * norm_syn 

            # compute misfit function
            L = 0.5 * (SI_syn - SI_obs)**2
            tr_chi[ir,0] = L 
            am_chi[ir,0] = L
            win_chi[ir,0,6] = SI_syn
            win_chi[ir,0,7] = SI_obs

            # adjoint source
            ddRsyn = dif1(dRsyn,dt_syn)
            adjsrc_T = -2. * (SI_syn - SI_obs) * dRsyn * norm_syn
            adjsrc_R = -2. * (SI_syn - SI_obs) * (  
                2. * np.sum(dRsyn * Tsyn) * norm_syn**2 * ddRsyn - 
                dTcomp * norm_syn
            )

            # filter adjoint source and taper it
            adjsrc_R = bandpass(adjsrc_R,dt_syn,freqmin,freqmax) * taper 
            adjsrc_T = bandpass(adjsrc_T,dt_syn,freqmin,freqmax) * taper 

            data = np.zeros((npt_syn,2))
            data[:,0] = t0_syn + np.arange(npt_syn) * dt_syn
            data[:,1] = adjsrc_T
            outname = f"{out_dir}/{bandname}/{self.netwk[i]}.{self.stnm[i]}.{self.chcode}T.adj.sem.npy"
            np.save(outname,data)
            data[:,1] = adjsrc_R
            outname = f"{out_dir}/{bandname}/{self.netwk[i]}.{self.stnm[i]}.{self.chcode}R.adj.sem.npy"
            np.save(outname,data)
        
        # save measurement files
        self._print_measure_info(bandname,tstart,tend,tr_chi,am_chi,win_chi)
        
    def cal_adj_source(self,ib:int):
        import os 
        bandname = 'T%03g_T%03g' %(self.Tmin[ib],self.Tmax[ib])
        if self.myrank == 0:
            print(f"preprocessing for band {bandname} ...")
            os.makedirs(f"{self.syndir}/OUTPUT_FILES/{bandname}",exist_ok=True)
        MPI.COMM_WORLD.Barrier()

        if self.meatype == "noise":
            self._cal_adjsrc_noise(ib,bandname)
        elif self.meatype == "tele":
            self._cal_adjsrc_tele(ib,bandname)
        else:
            self._cal_adjsrc_sks(ib,bandname)
    
    def sum_adj_source(self):
        import os 
        from utils import rotate_RT_to_EN
        nb = len(self.Tmax)
        if self.myrank == 0:
            print("\nsum adjoint source ...")
            os.makedirs(f"{self.syndir}/SEM",exist_ok=True)
        MPI.COMM_WORLD.Barrier()

        # loop each stations
        for ir in range(self.nsta_loc):
            i = ir + self._istart

            # init adjoint source
            adj = np.zeros((3,self.npt_syn))

            # check if rotation is requried
            need_rotate = False
            comps_read = ['E','N','Z']
            if 'R' in self.components or 'T' in self.components:
                need_rotate = True
                comps_read = ['R','T','Z']

            for ib in range(nb):
                bandname = 'T%03g_T%03g' %(self.Tmin[ib],self.Tmax[ib])
                input_dir = f"{self.syndir}/OUTPUT_FILES/{bandname}"

                for ic,ch in enumerate(comps_read):
                    name = f"{self.netwk[i]}.{self.stnm[i]}.{self.chcode}{ch}.adj.sem.npy"
                    if os.path.exists(f"{input_dir}/{name}"):
                        data = np.load(f"{input_dir}/{name}")
                        adj[ic,:] += data[:,1]
                
            # rotation if required
            if need_rotate:
                adj[0,:],adj[1,:] = rotate_RT_to_EN(adj[0,:],adj[1,:],self.bazd[i])

            # save to SEM/
            data = np.zeros((self.npt_syn,2))
            data[:,0] = self.t0_syn + self.dt_syn * np.arange(self.npt_syn)
            for ic,ch in enumerate(['E','N','Z']):
                name = f"{self.netwk[i]}.{self.stnm[i]}.{self.chcode}{ch}.adj.sem.npy"
                data[:,1] = adj[ic,:]
                np.save(f"{self.syndir}/SEM/{name}",data)
        
        # wait for jobs finish
        MPI.COMM_WORLD.Barrier()

        # rotate to XYZ
        self._rotate_ZNE_to_XYZ()

    def cleanup(self):
        from glob import glob 
        import os 
        import h5py 
        import shutil 
        from obspy.io.sac import SACTrace
        if self.myrank == 0:
            print("cleaning up ...")

        # first pack all obs/syn sacs to hdf5 
        for ib in range(len(self.Tmax)):
            bandname = 'T%03g_T%03g' %(self.Tmin[ib],self.Tmax[ib])
            for tag in ['obs','syn']:
                # all sac files
                sacfiles = glob(f"{self.syndir}/OUTPUT_FILES/{bandname}/*.sac.{tag}")
                fio = h5py.File(f"{self.syndir}/OUTPUT_FILES/seismogram.{tag}.{bandname}.h5","w")


                for i in range(len(sacfiles)):
                    tr = SACTrace.read(sacfiles[i])
                    if i == 0:
                        fio.attrs['dt'] = tr.delta  
                        fio.attrs['t0'] = tr.b 
                        fio.attrs['npts'] = tr.npts
                    
                    dsetname = tr.knetwk + "." + tr.kstnm + "." + tr.kcmpnm
                    fio.create_dataset(dsetname,shape=tr.data.shape,dtype='f4')
                    fio[dsetname][:] = tr.data 
        
            # clean bandname
            shutil.rmtree(f"{self.syndir}/OUTPUT_FILES/{bandname}",ignore_errors=True)

        # clean semd and sem.ascii
        filenames = glob(f"{self.syndir}/OUTPUT_FILES/*.semd")
        for f in filenames:
            os.remove(f)
        filenames = glob(f"{self.syndir}/OUTPUT_FILES/*.sem.npy")
        for f in filenames:
            os.remove(f)
        filenames = glob(f"{self.syndir}/SEM/*.sem.npy")
        for f in filenames:
            os.remove(f)

    def execute(self):
        from utils import rotate_EN_to_RT,cal_dist_baz

        # rotate seismograms from XYZ to ZNE
        self._rotate_XYZ_to_ZNE()

        # check if we need rotate to RT
        need_rt_rotate = False
        if 'R' in self.components or 'T' in self.components:
            need_rt_rotate = True 
            for i in range(self.nsta):
                _,self.bazd[i] = cal_dist_baz(self.evla,self.evlo,self.stla[i],self.stlo[i])

        # rotate to R/T if required
        if need_rt_rotate:
            for ir in range(self.nsta_loc):
                i = ir + self._istart
                temp_syn = np.zeros((2,self.npt_syn))
                data = np.zeros((self.npt_syn,2))
                for ic,ch in enumerate(['E','N']):
                    name = f"{self.netwk[i]}.{self.stnm[i]}.{self.chcode}{ch}.sem.npy"
                    filename = f"{self.syndir}/OUTPUT_FILES/{name}"
                    data[:,:] = np.load(filename)
                    temp_syn[ic,:] = data[:,1]

                # rotate to R/T
                ve = temp_syn[0,:]
                vn = temp_syn[1,:]
                vr,vt = rotate_EN_to_RT(ve,vn,self.bazd[i])
                temp_syn[0,:] = vr * 1.
                temp_syn[1,:] = vt * 1.

                # save to ascii
                for ic,ch in enumerate(['R','T']):
                    name = f"{self.netwk[i]}.{self.stnm[i]}.{self.chcode}{ch}.sem.npy"
                    filename = f"{self.syndir}/OUTPUT_FILES/{name}"
                    data[:,1] = temp_syn[ic,:] 
                    #data.tofile(filename)
                    np.save(filename,data)
            MPI.COMM_WORLD.Barrier()

        # save current synthetics as observation if required
        if self.run_opt == 1:
            self.save_forward()
            return 0
        
        # loop each frequency band to compute misfit/adjoint source 
        self._get_obs_info()
        for ib in range(len(self.Tmax)):
            self.cal_adj_source(ib)
            MPI.COMM_WORLD.Barrier()
        
        # sum adjoint source
        self.sum_adj_source()

        # pack SACs to h5
        if self.myrank == 0: self.cleanup()

        # sync
        MPI.COMM_WORLD.Barrier()

def main():
    import sys 
    if len(sys.argv) != 5:
        print("Usage: python run_preprocess.py measure_type iter evtid run_opt")
        exit(1)

    # get input parameter
    mtype = sys.argv[1]
    iter0 = int(sys.argv[2])
    evtid = sys.argv[3]
    run_opt = int(sys.argv[4])

    # init operator
    op = FwatPreOP(mtype,iter0,evtid,run_opt)

    # run
    op.execute()

    # finalize
    MPI.Finalize()

if __name__ == "__main__":
    main()