import numpy as np 
from mpi4py import MPI

class FwatPreOP:
    def __init__(self,measure_type:str,iter:int,evtid:str,run_opt:int):
        # import packages
        from utils import read_params
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
        self.stla = np.float64(statxt[:,2])
        self.stlo = np.float64(statxt[:,3])
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
        from cube2sph_rotate import rotate_seismo

        # set parameters
        fn_matrix = f"{self.SRC_REC}/rot_{self.evtid}"
        rotate = "XYZ->NEZ"
        from_dir = f"{self.syndir}/OUTPUT_FILES/"
        to_dir = f"{self.syndir}/OUTPUT_FILES/"
        from_template='${nt}.${sta}.BX${comp}.semd'
        to_template='${nt}.${sta}.BX${comp}.sem.ascii'

        # rotate seismograms from XYZ to ZNE
        rotate_seismo(fn_matrix,rotate,from_dir,
                      to_dir,from_template,to_template)
    
    def _rotate_ZNE_to_XYZ(self):
        from cube2sph_rotate import rotate_seismo

        # set parameters
        fn_matrix = f"{self.SRC_REC}/rot_{self.evtid}"
        rotate = "XYZ<-NEZ"
        from_dir = f"{self.syndir}/SEM/"
        to_dir = f"{self.syndir}/SEM/"
        from_template='${nt}.${sta}.BX${comp}.adj.sem.ascii'
        to_template='${nt}.${sta}.BX${comp}.adj'

        # rotate seismograms from XYZ to ZNE
        rotate_seismo(fn_matrix,rotate,from_dir,
                      to_dir,from_template,to_template,'ascii')
        
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
        data = np.zeros((npt_syn,2))
        data[:,0] = np.arange(npt_syn) * dt_syn + t0_syn 
        for ir in range(self.nsta_loc):
            i = ir + self._istart
            for ic in range(ncomp):
                ch = components[ic]

                # fetch synthetics 
                code = self._get_station_code(i,ic)
                iloc = self._syncmp.index(ch)
                data[:,1] = self.seismo[ir,iloc,:]

                # special handling for each type
                if self.meatype == "tele":
                    # convolve
                    tr.data = convolve(data[:,1],stf[ic,:],'same') * dt_syn    
                    tr.b =  self.t_inj - self.t_ref[i]
                elif self.meatype == "sks":
                    tr.b =  self.t_inj - self.t_ref[i]
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
                tr.knetwk = self.netwk[i],
                tr.kstnm = self.stnm[i],
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
        from utils import alloc_mpi_jobs

        ncomp = tr_chi.shape[1]
        istart,iend = alloc_mpi_jobs(self.nsta,self.nprocs,self.myrank)
        nsta_loc = iend - istart + 1
        imeas = self.pdict['IMEAS']

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
                    i = ir + istart 
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
        window = np.zeros((nsta_loc,ncomp,20))
        tr_chi = np.zeros((nsta_loc,ncomp))
        am_chi = np.zeros((nsta_loc,ncomp))

        # loop every station to do preprocessing
        for ir in range(nsta_loc):
            i = ir + self._istart
            for icomp in range(self.ncomp):
                name = self._get_station_code(i,icomp)
                obs_tr = SACTrace.read(f'{self.DATA_DIR}/{self.evtid}/{name}.sac')

                # synthetic data
                iloc = self._syncmp.index(self.components[icomp])
                syn_tr = self.seismo[ir,iloc,:] * 1.


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
                tr_chi[ir,icomp],am_chi[ir,icomp],window[ir,icomp,:],adjsrc =   \
                    measure_adj(t0_inp,dt_inp,npt_cut,t0_syn,dt_syn,npt_syn,
                                tstart[ir],tend[ir],imeas,self.Tmax[ib]*1.01,
                                self.Tmin[ib]*0.99,verbose,dat_inp,
                                syn_inp,)

                # save adjoint source
                data = np.zeros((npt_syn,2))
                data[:,0] = t0_syn + np.arange(npt_syn) * dt_syn
                data[:,1] = adjsrc
                name = self._get_station_code(i,icomp) + ".adj.sem.ascii"
                np.savetxt(f"{out_dir}/{bandname}/{name}",data,fmt='%g')

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
        self._print_measure_info(bandname,tstart,tend,tr_chi,am_chi,window)
    
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
        window = np.zeros((nsta_loc,ncomp,20))
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

                # syn data
                iloc = self._syncmp.index(self.components[ic])
                syn_data = self.seismo[ir,iloc,:] * 1. 

                # read obs data
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
        ic = 0
        for i in range(ncomp):
            if self.components[i] == 'Z':
                ic = i 
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
                
                tr_chi[ir,ic],am_chi[ir,ic],window[ir,ic,:],adjsrc =  \
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
                name = self._get_station_code(i,ic) + ".adj.sem.ascii"
                np.savetxt(f"{out_dir}/{bandname}/{name}",data,fmt='%g')

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

        # print info and save MEASUREMENTS file
        self._print_measure_info(bandname,tstart,tend,tr_chi,am_chi,window)
    
    def _cal_adjsrc_sks(self,ib:int,bandname:str):
        from obspy.io.sac import SACTrace
        from utils import interpolate_syn
        from utils import bandpass,cumtrapz1,dif1

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
        tstart = np.zeros((nsta_loc))
        tend = np.zeros((nsta_loc))
        window = np.zeros((nsta_loc,1,20))
        tr_chi = np.zeros((nsta_loc,1))
        am_chi = np.zeros((nsta_loc,1))

        # loop each station
        for ir in range(nsta_loc):
            i = ir + self._istart
            obs_data = np.zeros((2,npt_syn))
            syn_data = np.zeros((2,npt_syn))

            # time window
            tstart[ir] = self.t_ref[i] - t_inj - win_tb
            tend[ir] = self.t_ref[i] - t_inj + win_te

            for ic in range(self.ncomp):
                name = self._get_station_code(i,ic)

                # read data
                obs_tr = SACTrace.read(f"{self.DATA_DIR}/{self.evtid}/{name}.sac")
                odata = obs_tr.data 

                # synthetic
                iloc = self._syncmp.index(self.components[ic])
                sdata = self.seismo[ir,iloc,:] * 1.

                # data processing
                odata = bandpass(odata,dt_obs,freqmin,freqmax)
                sdata = bandpass(sdata,dt_syn,freqmin,freqmax)

                # interpolate the obs/syn data to the same sampling of syn data
                t0_inp = self.t_ref[i] - win_tb
                u1 = interpolate_syn(odata,t0_obs + self.t_ref[i],dt_obs,npt_obs,t0_inp,dt_syn,npt2)
                w1 = interpolate_syn(sdata,t_inj,dt_syn,npt_syn,t0_inp,dt_syn,npt2)
                u = interpolate_syn(u1,t0_inp,dt_syn,npt2,t_inj,dt_syn,npt_syn)
                w = interpolate_syn(w1,t0_inp,dt_syn,npt2,t_inj,dt_syn,npt_syn)     

                # save to global array   
                obs_data[ic,:] = u
                syn_data[ic,:] = w 

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
                tr.data = u
                tr.write(f"{out_dir}/{bandname}/{name}.sac.obs")
                tr.data = w
                tr.write(f"{out_dir}/{bandname}/{name}.sac.syn")

            # compute SI
            ic_r = self.components.index("R")
            ic_t = self.components.index("T")
            Robs_dot = dif1(obs_data[ic_r,:],dt_syn)
            norm = cumtrapz1(Robs_dot**2,dt_syn)[-1]
            if norm < 1.0e-15:
                norm = 0.
            else:
                norm = 1. / norm
            Tobs = obs_data[ic_t,:]
            Tsyn = syn_data[ic_t,:]
            SI_obs = -2. * cumtrapz1(Robs_dot*Tobs,dt_syn)[-1] * norm
            SI_syn = -2. * cumtrapz1(Robs_dot*Tsyn,dt_syn)[-1] * norm
            L = 0.5 * (SI_syn - SI_obs)**2
            tr_chi[ir,0] = L 
            am_chi[ir,0] = L
            window[ir,0,6] = SI_syn
            window[ir,0,7] = SI_obs

            # adjoint source
            adjsrc_T = -2. * norm * (SI_syn - SI_obs) * Robs_dot
            data = np.zeros((npt_syn,2))
            data[:,0] = t0_syn + np.arange(npt_syn) * dt_syn
            data[:,1] = adjsrc_T
            outname = f"{out_dir}/{bandname}/{self.netwk[i]}.{self.stnm[i]}.{self.chcode}T.adj.sem.ascii"
            np.savetxt(outname,data,fmt='%g')
        
        # save measurement files
        self._print_measure_info(bandname,tstart,tend,tr_chi,am_chi,window)
        
    def cal_adj_source(self,ib:int):
        import os 
        bandname = 'T%03g_T%03g' %(self.Tmin[ib],self.Tmax[ib])
        if self.myrank == 0:
            print(f"preprocessing for band {bandname} ...")
            os.makedirs(f"{self.syndir}/OUTPUT_FILES/{bandname}",exist_ok=True)

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
        for i in range(self.myrank,self.nsta,self.nprocs):
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
                    name = f"{self.netwk[i]}.{self.stnm[i]}.{self.chcode}{ch}.adj.sem.ascii"
                    if os.path.exists(f"{input_dir}/{name}"):
                        data = np.loadtxt(f"{input_dir}/{name}")
                        adj[ic,:] += data[:,1]
                
            # rotation if required
            if need_rotate:
                adj[0,:],adj[1,:] = rotate_RT_to_EN(adj[0,:],adj[1,:],self.bazd[i])

            # save to SEM/
            data = np.zeros((self.npt_syn,2))
            data[:,0] = self.t0_syn + self.dt_syn * np.arange(self.npt_syn)
            for ic,ch in enumerate(['E','N','Z']):
                name = f"{self.netwk[i]}.{self.stnm[i]}.{self.chcode}{ch}.adj.sem.ascii"
                data[:,1] = adj[ic,:]
                np.savetxt(f"{self.syndir}/SEM/{name}",data,fmt="%g")
        
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
                fio.attrs['dt'] = self.dt_syn 
                fio.attrs['t0'] = self.t0_syn
                fio.attrs['npts'] = self.npt_syn

                for i in range(len(sacfiles)):
                    tr = SACTrace.read(sacfiles[i])
                    dsetname = tr.knetwk + "." + tr.kstnm + "." + tr.kcmpnm
                    fio.create_dataset(dsetname,shape=tr.data.shape,dtype='f4')
                    fio[dsetname][:] = tr.data 
        
            # clean bandname
            shutil.rmtree(f"{self.syndir}/OUTPUT_FILES/{bandname}",ignore_errors=True)

        # clean semd and sem.ascii
        filenames = glob(f"{self.syndir}/OUTPUT_FILES/*.semd")
        for f in filenames:
            os.remove(f)
        filenames = glob(f"{self.syndir}/OUTPUT_FILES/*.sem.ascii")
        for f in filenames:
            os.remove(f)
        filenames = glob(f"{self.syndir}/SEM/*.sem.ascii")
        for f in filenames:
            os.remove(f)

    def execute(self):
        from utils import rotate_EN_to_RT,cal_dist_baz
        from utils import alloc_mpi_jobs

        # rotate seismograms from XYZ to ZNE
        self._rotate_XYZ_to_ZNE()

        # check if we need rotate to RT
        need_rt_rotate = False
        if 'R' in self.components or 'T' in self.components:
            need_rt_rotate = True 
            for i in range(self.nsta):
                _,self.bazd[i] = cal_dist_baz(self.evla,self.evlo,self.stla[i],self.stlo[i])
        
        # allocate syn arrays
        istart,iend = alloc_mpi_jobs(self.nsta,self.nprocs,self.myrank)
        self.nsta_loc = iend - istart + 1
        self._istart = istart 
        self.seismo = np.zeros((self.nsta_loc,3,self.npt_syn))
        self._syncmp = ['E','N','Z']

        # load synthetic data 
        for ir in range(self.nsta_loc):
            i = ir + self._istart
            for ic,ch in enumerate(self._syncmp):
                name = f"{self.netwk[i]}.{self.stnm[i]}.{self.chcode}{ch}.sem.ascii"
                filename = f"{self.syndir}/OUTPUT_FILES/{name}"
                data = np.loadtxt(filename)
                self.seismo[ir,ic,:] = data[:,1]

        # rotate to R/T if required
        if need_rt_rotate:
            for ir in range(self.nsta_loc):
                i = ir + self._istart

                # rotate to R/T
                ve = self.seismo[ir,0,:]
                vn = self.seismo[ir,1,:]
                vr,vt = rotate_EN_to_RT(ve,vn,self.bazd[i])
                self.seismo[ir,0,:] = vr * 1.
                self.seismo[ir,1,:] = vt * 1.
            self._syncmp = ['R','T','Z']

        # save current synthetics as observation if required
        if self.run_opt == 1:
            self.save_forward()
            return 0
        
        # loop each frequency band to compute misfit/adjoint source 
        self._get_obs_info()
        for ib in range(len(self.Tmax)):
            self.cal_adj_source(ib)
        
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

if __name__ == "__main__":
    main()