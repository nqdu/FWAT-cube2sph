import numpy as np 
from mpi4py import MPI 
from fwat.const import PARAM_FILE

def _get_rotate_matrix(azd,bazd,rotate_src=True,rotate_sta = True):
    from numpy import sin,cos
    R_s = np.eye(3,3)
    R_r = np.eye(3,3)

    if rotate_src:
        az = np.deg2rad(azd)
        R_s = np.array([
            [sin(az),cos(az),0],
            [cos(az),-sin(az),0],
            [0,0,1.]
        ])

    if rotate_sta:
        baz = np.deg2rad(bazd)
        R_r = np.array([
            [-sin(baz),-cos(baz),0],
            [-cos(baz),sin(baz),0],
            [0.,0.,1.]
        ])

    return R_s,R_r

class NoiseMC_PreOP():
    def __init__(self, measure_type:str, iter:int, evtid:str, run_opt:int):
        # import packages
        from .utils import read_params,cal_dist_az_baz
        from .utils import get_source_loc
        import h5py 
        comm = MPI.COMM_WORLD
        self.myrank = comm.Get_rank()
        self.nprocs = comm.Get_size()

        # private vars
        from fwat import const as fpath
        self.SRC_REC = fpath.SRC_REC
        self.SOLVER = fpath.SOLVER
        self.DATA_DIR = fpath.DATA_DIR
        self.MISFIT = fpath.MISFIT

        # backup
        self.meatype = measure_type
        self.iter = iter 
        self.evtid = evtid
        self.run_opt = run_opt

        # read FWAT params below
        # frequency band
        pdict = read_params(f"{PARAM_FILE}")['measure'][self.meatype]
        self.Tmin = [x[0] for x in pdict['FILTER_BANDS']]
        self.Tmax = [x[1] for x in pdict['FILTER_BANDS']]

        # channel code
        self.chcode = pdict['CH_CODE']

        # backup pdict for further usage
        self.pdict = pdict

        # adjoint source
        self.adjsrc_type = pdict['ADJSRC_TYPE']

        # path to simulation directory
        mdir = "M%02d" %(iter)
        if self.run_opt == 2:
            mdir = mdir + ".ls"
        self.mod = mdir

        # get all synthetic directories
        self._get_mc_channel()

        # adjoint source
        self.adjsrc_type = self.pdict['ADJSRC_TYPE']

        # read sourcer/receiver locations below
        # ---------------------
        # source loc
        sourcefile = f"{self.SRC_REC}/sources.dat.{measure_type}"
        self.evla,self.evlo,self.evdp = get_source_loc(f"{evtid}",sourcefile)

        # station info
        # read station coordinates
        self._get_station_info()

        # read simulation info dt,t0,npts
        fio = h5py.File(f"{self.syndirs[0]}/OUTPUT_FILES/seismograms.h5","r")
        t = fio[list(fio.keys())[0]][:,0]
        fio.close()
        self.t0_syn = t[0]
        self.dt_syn = t[1] - t[0]
        self.npt_syn = len(t)
    
    def _get_mc_channel(self):
        import os 

        # get source comps and receiver comps
        self.cc_comps:list[str] = sorted(self.pdict['CC_COMPS'])
        self.ncomp = len(self.cc_comps)

        # sanity check
        self._sanity_check()

        # init
        comps_temp = set()
        scomps = [cc[0] for cc in self.cc_comps]

        for ch in ['E','N','Z']:
            name = f"{self.SRC_REC}/STATIONS_{self.evtid}_{ch}"
            if ch in scomps and os.path.exists(name):
                comps_temp.add(ch)

        for ch in ['R','T']:
            name1 = f"{self.SRC_REC}/STATIONS_{self.evtid}_{ch}"
            if ch in scomps and os.path.exists(name1):
                comps_temp.add("E")
                comps_temp.add("N")

        # recheck all cc_comps
        cc_comps = []
        for ic in range(self.ncomp):
            flag = scomps[ic] in comps_temp
            flag1 = scomps[ic] in ['Z','E','N']
            flag2 = scomps[ic] in ['R','T']
            if flag1 and  flag:
                cc_comps.append(self.cc_comps[ic])

            if flag2 and 'E' in comps_temp and 'N' in comps_temp :
                cc_comps.append(self.cc_comps[ic])
        
        # reset cc_comps and ncomps
        self.cc_comps = sorted(cc_comps)
        self.ncomp = len(cc_comps)
        if self.myrank == 0: print(f"{self.evtid}: used cc_comps =  {cc_comps}\n")

        comps = sorted(list(comps_temp))
        self.scomp_syn = comps # only ZNE
        self.syndirs:list[str] = []
        for ic in range(len(comps)):
            syndir = f"{self.SOLVER}/{self.mod}/{self.evtid}_{comps[ic]}/"
            self.syndirs.append(syndir)

    def _sanity_check(self):
        # make sure [RT][NE] donnot co-exist
        scomps = [self.cc_comps[i][0] for i in range(self.ncomp)]
        rcomps = [self.cc_comps[i][1] for i in range(self.ncomp)]

        if ('N' in scomps or 'E' in scomps) and (('R' in scomps or 'T' in scomps)):
            print("[NE] and [RT] cannot both exist in source-side CC-components")
            exit(1)

        if ('N' in rcomps or 'E' in rcomps) and (('R' in rcomps or 'T' in rcomps)):
            print("[NE] and [RT] cannot both exist in receiver-side CC-components")
            exit(1)

        # check adjsrc_type
        if self.adjsrc_type not in [5,7,'exp_phase','cc_time']:
            if self.myrank == 0:
                print(f"only [5,7,'exp_phase','cc_time'] adjoint sources are supported!")
                print(f"adjsrc_type = {self.adjsrc_type}")
            exit(1)

    
    def _get_station_info(self):
        from .utils  import cal_dist_az_baz
        self.stainfo = {}
        for ic in range(self.ncomp):
            chs = self.cc_comps[ic][0]
            stationfile = f'{self.SRC_REC}/STATIONS_{self.evtid}_{chs}_globe'
            statxt = np.loadtxt(stationfile,dtype=str,ndmin=2)
            nsta = statxt.shape[0]
            netwk = statxt[:,1]
            stnm = statxt[:,0]
            stla = statxt[:,2].astype(float)
            stlo = statxt[:,3].astype(float)

            for ir in range(nsta):
                name = netwk[ir] + '.' + stnm[ir]
                if name not in self.stainfo:
                    _,az,baz = cal_dist_az_baz(self.evla,self.evlo,stla[ir],stlo[ir])
                    self.stainfo[name] = [stlo[ir],stla[ir],az,baz]
        
        # get all names
        self.sta_names: list[str] = sorted(self.stainfo.keys())

    def _rotate_XYZ_to_ZNE(self):
        from .cube2sph_rotate import rotate_seismo_fwd

        # loop each synthetic directory
        for ic in range(len(self.scomp_syn)):
            fn_matrix = f"{self.SRC_REC}/rot_{self.evtid}"
            from_dir = f"{self.syndirs[ic]}/OUTPUT_FILES/"
            to_dir = f"{self.syndirs[ic]}/OUTPUT_FILES/"
            from_template='${nt}.${sta}.BX${comp}.semd'
            to_template='${nt}.${sta}.BX${comp}.sem.npy'

            # rotate seismograms from XYZ to ZNE
            rotate_seismo_fwd(fn_matrix,from_dir,to_dir,from_template,to_template)

    def _rotate_ZNE_to_XYZ(self):
        from .cube2sph_rotate import rotate_seismo_adj

        # loop each synthetic directory
        for ic in range(len(self.scomp_syn)):
            fn_matrix = f"{self.SRC_REC}/rot_{self.evtid}"
            from_dir = f"{self.syndirs[ic]}/SEM/"
            to_dir = f"{self.syndirs[ic]}/SEM/"
            from_template='${nt}.${sta}.BX${comp}.adj.sem.npy'
            to_template='${nt}.${sta}.BX${comp}.adj'

            # rotate seismograms from XYZ to ZNE
            rotate_seismo_adj(fn_matrix,from_dir,to_dir,from_template,to_template)

    def _get_bandname(self,ib:int):
        return 'T%03g_T%03g' %(self.Tmin[ib],self.Tmax[ib])

    def _rotate_seismogram_to_RT(self,):
        from .utils import rotate_EN_to_RT
        from .utils import alloc_mpi_jobs

        # get receiver components
        rcomp = [self.cc_comps[i][1] for i in range(self.ncomp)]

        # allocate jobs
        keys = sorted(self.stainfo.keys())
        nsta = len(keys) 
        istart,iend = alloc_mpi_jobs(nsta,self.nprocs,self.myrank)

        # loop each synthetic directory
        for ie in range(len(self.scomp_syn)):
            need_rt_rotate = False
            if 'R' in rcomp or 'T' in rcomp:
                need_rt_rotate = True 

            # check if need rotate
            if not need_rt_rotate: continue

            # rotate
            for i in range(istart,iend+1):
                temp_syn = np.zeros((2,self.npt_syn))
                data = np.zeros((self.npt_syn,2))
                for ic,ch in enumerate(['E','N']):
                    name = f"{keys[i]}.{self.chcode}{ch}.sem.npy"
                    filename = f"{self.syndirs[ie]}/OUTPUT_FILES/{name}"
                    data[:,:] = np.load(filename)
                    temp_syn[ic,:] = data[:,1]

                # rotate to R/T
                ve = temp_syn[0,:]
                vn = temp_syn[1,:]
                bazd = self.stainfo[keys[i]][-1]
                vr,vt = rotate_EN_to_RT(ve,vn,bazd)
                temp_syn[0,:] = vr * 1.
                temp_syn[1,:] = vt * 1.

                # save to ascii
                for ic,ch in enumerate(['R','T']):
                    name = f"{keys[i]}.{self.chcode}{ch}.sem.npy"
                    filename = f"{self.syndirs[ie]}/OUTPUT_FILES/{name}"
                    data[:,1] = temp_syn[ic,:] 
                    #data.tofile(filename)
                    np.save(filename,data)

        # sync
        MPI.COMM_WORLD.Barrier()


    def save_forward(self):
        import os 
        from obspy.io.sac import SACTrace
        from fwat.measure.utils import cumtrapz1,alloc_mpi_jobs
        from fwat.measure.utils import rotate_EN_to_RT

        # get some vars
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

        if self.myrank == 0:
            print("Synthetic Observations ...")

        # loop every source component
        for ic in range(self.ncomp):
            chs = self.cc_comps[ic][0]
            chr = self.cc_comps[ic][1]
            
            # make dir
            outdir = f"{self.DATA_DIR}/{self.evtid}_{chs}"
            os.makedirs(outdir,exist_ok=True)

            # load stations
            statxt = np.loadtxt(f"{self.SRC_REC}/STATIONS_{self.evtid}_{chs}_globe",dtype=str,ndmin=2)
            sta_names = [statxt[i,1] + "." + statxt[i,0] for i in range(statxt.shape[0])]

            # allocate jobs
            istart,iend = alloc_mpi_jobs(len(sta_names),self.nprocs,self.myrank)

            # now loop to read stations
            for i in range(istart,iend+1):
                # get station code
                code = f"{sta_names[i]}.{self.chcode}{chr}"

                # load [NEZ]? components if required
                if chs in ['N','E','Z']:
                    syndir = f"{self.SOLVER}/{self.mod}/{self.evtid}_{chs}/"
                    filename = f"{syndir}/OUTPUT_FILES/{code}.sem.npy"
                    data = np.load(filename)
                    tr.data = data[:,1] * 1.

                else: # chs is in ['R','T']
                    
                    # load E/N data
                    data = np.zeros((2,npt_syn))
                    for ic0,ch0 in enumerate(['N','E']):
                        syndir = f"{self.SOLVER}/{self.mod}/{self.evtid}_{ch0}/"
                        filename = f"{syndir}/OUTPUT_FILES/{code}.sem.npy"
                        data[ic0,:] = np.load(filename)[:,1]
                    
                    # rotation
                    az = self.stainfo[sta_names[i]][2]
                    ve = data[0,:]
                    vn = data[1,:]
                    vr,vt = rotate_EN_to_RT(ve,vn,az)
                    vr = -vr 
                    vt = -vt 

                    # choose components
                    if chs == 'R':
                        tr.data = vr * 1. 
                    else:
                        tr.data = vt * 1.

                # save sac
                if self.pdict['USE_EGF'] == False:
                    if i == istart and myrank == 0:
                        print("EGF => CCF ...")
                    tr.data = -cumtrapz1(tr.data,dt_syn)
                tr.b = t0_syn

                # channel and others
                tr.kcmpnm = f"{self.chcode}{chr}"
                tr.knetwk = sta_names[i].split('.')[0]
                tr.kstnm = sta_names[i].split('.')[1]
                info = self.stainfo[sta_names[i]]
                tr.stla = info[1]
                tr.stlo = info[0]

                # save to sac
                filename = f"{outdir}/{code}.sac"
                tr.write(filename)
            
            # sync
            MPI.COMM_WORLD.Barrier()
    
    def cal_adj_source(self,ib:int):
        from obspy.io.sac import SACTrace
        from .utils import interpolate_syn,dif1
        from .utils import bandpass,alloc_mpi_jobs
        from .utils import rotate_EN_to_RT
        import os 

        # get vars
        npt_syn = self.npt_syn
        t0_syn = self.t0_syn
        dt_syn = self.dt_syn

        # new dt/nt for interpolate
        dt_inp = 0.01
        t0_inp = -10.
        t1_inp = t0_syn + (npt_syn - 1) * dt_syn
        npt_cut = int((t1_inp - t0_inp) / dt_inp) + 1

        # get snr
        snr_threshold:float = -1.
        if 'SNR_THRESHOLD' in self.pdict:
            snr_threshold = self.pdict['SNR_THRESHOLD'][ib]

        # band name
        bandname = self._get_bandname(ib)
        freqmin = 1. / self.Tmax[ib]
        freqmax = 1. / self.Tmin[ib]

        # group velocity window
        vmin_list = [x[0] for x in self.pdict['GROUPVEL_WIN']]
        vmax_list = [x[1] for x in self.pdict['GROUPVEL_WIN']]

        # log
        if self.myrank == 0:
            print(f"preprocessing for band {bandname} ...")
            for ic in range(len(self.scomp_syn)):
                outdir = f"{self.syndirs[ic]}/OUTPUT_FILES/{bandname}"
                os.makedirs(outdir,exist_ok=True)
        MPI.COMM_WORLD.Barrier()

        # allocate jobs
        nsta = len(self.sta_names)
        istart,iend = alloc_mpi_jobs(nsta,self.nprocs,self.myrank)
        nsta_loc = iend - istart + 1

        # misfits
        ncomp = self.ncomp
        tstart = np.zeros((nsta_loc))
        tend = np.zeros((nsta_loc))
        win_chi = np.zeros((nsta_loc,ncomp,20))
        tr_chi = np.zeros((nsta_loc,ncomp))
        am_chi = np.zeros((nsta_loc,ncomp))

        # components to read
        rcomp = [self.cc_comps[i][1] for i in range(self.ncomp)]
        if 'R' in rcomp or 'T' in rcomp:
            rcomp = ['R','T','Z']
        else:
            rcomp = ['E','N','Z']
        scomp = [self.cc_comps[i][0] for i in range(self.ncomp)]
        if 'R' in scomp or 'T' in scomp:
            scomp = ['R','T','Z']
        else:
            scomp = ['E','N','Z']

        # names for each cc comps
        if ib == 0:
            sta_names_comps:list[set] = []
            for ic in range(3):
                filename = f"{self.SRC_REC}/STATIONS_{self.evtid}_{scomp[ic]}_globe"
                if os.path.exists(filename):
                    # load stations
                    statxt = np.loadtxt(filename,dtype=str,ndmin=2)
                    sta_names = set([statxt[i,1] + "." + statxt[i,0] for i in range(statxt.shape[0])])
                else:
                    sta_names = set()
                sta_names_comps.append(sta_names)
            self.sta_names_comps = sta_names_comps
        
        # loop each station
        for ir in range(nsta_loc):
            i = ir + istart 

            # allocate space for adjoint source
            adj_src_all = np.zeros((3,3,npt_syn))

            # loop each CC-component
            for ic in range(self.ncomp):
                chs = self.cc_comps[ic][0]
                chr = self.cc_comps[ic][1]
                i_s = scomp.index(chs)
                i_r = rcomp.index(chr)

                # check if this station exists in this source component
                if self.sta_names[i] not in self.sta_names_comps[i_s]:
                    continue

                # load obs_data
                code = f"{self.sta_names[i]}.{self.chcode}{chr}"
                obs_tr = SACTrace.read(f"{self.DATA_DIR}/{self.evtid}_{chs}/{code}.sac")
                t0_obs = obs_tr.b 
                dt_obs = obs_tr.delta 
                npt_obs = obs_tr.npts
                npt1_inp = int((npt_obs - 1) * dt_obs / dt_inp)
                dist = obs_tr.dist 

                # load synthetic data
                if chs in ['N','E','Z']:
                    syndir = f"{self.SOLVER}/{self.mod}/{self.evtid}_{chs}/"
                    filename = f"{syndir}/OUTPUT_FILES/{code}.sem.npy"
                    syn_tr = np.load(filename)[:,1]
                else: # chs is in ['R','T']
                    # load E/N data
                    data = np.zeros((2,npt_syn))
                    for ic0,ch0 in enumerate(['N','E']):
                        syndir = f"{self.SOLVER}/{self.mod}/{self.evtid}_{ch0}/"
                        filename = f"{syndir}/OUTPUT_FILES/{code}.sem.npy"
                        data[ic0,:] = np.load(filename)[:,1]
                    
                    # rotation
                    az = self.stainfo[self.sta_names[i]][2]
                    ve = data[0,:]
                    vn = data[1,:]
                    vr,vt = rotate_EN_to_RT(ve,vn,az)
                    vr = -vr 
                    vt = -vt 

                    # choose components
                    if chs == 'R':
                        syn_tr = vr * 1. 
                    else:
                        syn_tr = vt * 1.

                # interp obs/syn
                dat_inp1 = interpolate_syn(obs_tr.data,t0_obs,dt_obs,npt_obs,
                                        t0_obs + dt_inp,dt_inp,npt1_inp)
                syn_inp = interpolate_syn(syn_tr,t0_syn,dt_syn,npt_syn,
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

                # find amplitude of the in the window to normalize
                dist = obs_tr.dist
                win_b = np.floor((dist / vmax_list[ib] - self.Tmax[ib] / 2. + t0_inp) / dt_inp)
                win_e = np.floor((dist / vmin_list[ib] + self.Tmax[ib] / 2. - t0_inp) / dt_inp)
                win_b = int(max(win_b,0))
                win_e = int(min(win_e,len(dat_inp)-1))
                amp = np.max(np.abs(dat_inp[win_b:win_e]))

                # check snr
                std_in_tail = np.std(dat_inp[win_e+1:])
                if std_in_tail < 1.0e-20:
                    snr = 1.0e20
                else:
                    snr = amp / std_in_tail
                if snr > snr_threshold:
                    dat_inp *= np.max(np.abs(syn_inp[win_b:win_e])) / amp
                else:
                    # set dat/syn to zero, it will contribute nothing
                    dat_inp *= 0. 
                    syn_inp *= 0.

                # compute time window
                tstart[ir] = dist / vmax_list[ib] - self.Tmax[ib] * 0.5 
                tend[ir] = dist / vmin_list[ib] + self.Tmax[ib] * 0.5 
                tstart[ir] = max(tstart[ir],t0_inp)
                tend[ir] = min(tend[ir],t0_obs+(npt_obs-1)*dt_obs)

                # compute misfits and adjoint source
                if self.adjsrc_type == 'exp_phase':
                    from fwat.adjoint.exp_phase_misfit import measure_adj_exphase
                    tr_chi[ir,ic],am_chi[ir,ic],win_chi[ir,ic,:],adjsrc =  \
                        measure_adj_exphase(dat_inp,syn_inp,
                                            t0_inp,dt_inp,npt_cut,
                                            self.Tmin[ib],self.Tmax[ib],
                                            tstart[ir],tend[ir]
                        )
                    # reinterpolate adjoint source
                    adjsrc = interpolate_syn(adjsrc,t0_inp,dt_inp,npt_cut,
                                             t0_syn,dt_syn,npt_syn)
                elif self.adjsrc_type == 'cc_time':
                    from fwat.adjoint.cc_misfit import measure_adj_cc
                    tr_chi[ir,ic],am_chi[ir,ic],win_chi[ir,ic,:],adjsrc =  \
                        measure_adj_cc(dat_inp,syn_inp,
                                       t0_inp,dt_inp,npt_cut,
                                       self.Tmin[ib],self.Tmax[ib],
                                       tstart[ir],tend[ir]
                        )
                    # reinterpolate adjoint source
                    adjsrc = interpolate_syn(adjsrc,t0_inp,dt_inp,npt_cut,
                                             t0_syn,dt_syn,npt_syn)
                else:
                    from .measure import measure_adj
                    verbose = (self.myrank == 0) and (ir == 0)
                    tr_chi[ir,ic],am_chi[ir,ic],win_chi[ir,ic,:],adjsrc =   \
                        measure_adj(t0_inp,dt_inp,npt_cut,t0_syn,dt_syn,npt_syn,
                                    tstart[ir],tend[ir],self.adjsrc_type,self.Tmax[ib]*1.01,
                                    self.Tmin[ib]*0.99,verbose,dat_inp,
                                    syn_inp)

                # save adjoint source to adj_src_all
                adj_src_all[i_s,i_r,:] = adjsrc.copy()

                # save obs and syn data as sac
                if chs in ['N','E','Z']:
                    syndir = f"{self.SOLVER}/{self.mod}/{self.evtid}_{chs}/"
                else:
                    if chs == 'R':
                        syndir = f"{self.SOLVER}/{self.mod}/{self.evtid}_E/"
                    else:
                        syndir = f"{self.SOLVER}/{self.mod}/{self.evtid}_N/"
                outdir = f"{syndir}/OUTPUT_FILES/{bandname}"
                obs_tr.delta = dt_inp 
                obs_tr.b = t0_inp 
                obs_tr.data = dat_inp 
                obs_tr.write(f"{outdir}/{code}.sac.obs")
                obs_tr.data = syn_inp
                obs_tr.write(f"{outdir}/{code}.sac.syn")
            # end for loop cc_comps

            # rotate adjoint source to ENZ-ENZ
            azd,bazd = self.stainfo[self.sta_names[i]][2:]
            rotate_src = (scomp == ['R','T','Z'])
            rotate_sta = (rcomp == ['R','T','Z'])
            R_s,R_r = _get_rotate_matrix(azd,bazd,rotate_src=rotate_src,rotate_sta=rotate_sta)
            adj_src_all = np.einsum("ip,jq,pqk-> ijk",R_s,R_r,adj_src_all)

            # save adjoint source
            rcomp = ['E','N','Z']
            data = np.zeros((npt_syn,2))
            data[:,0] = t0_syn + np.arange(npt_syn) * dt_syn
            for ic in range(len(self.scomp_syn)):
                i_s = rcomp.index(self.scomp_syn[ic])
                for i_r in range(3):
                    code = f"{self.sta_names[i]}.{self.chcode}{rcomp[i_r]}"
                    filename = f"{self.syndirs[ic]}/OUTPUT_FILES/{bandname}/{code}.adj.sem.npy"
                    data[:,1] = adj_src_all[i_s,i_r,:]
                    np.save(filename,data)

        # end for each station

        # print measure_adj information
        self._print_measure_info(bandname,tstart,tend,tr_chi,am_chi,win_chi)


    def _print_measure_info(self,bandname:str,tstart:np.ndarray,tend:np.ndarray,
                            tr_chi:np.ndarray,am_chi:np.ndarray,
                            window_chi:np.ndarray):
        """ 
        print measure_adj info and save to file 
        
        Parameters
        ----------
        bandname: str
            band name
        tstart: np.ndarray, shape(nsta_loc,)
            start time of measurement window
        tend: np.ndarray,shape(nsta_loc,)
            end time of measurement window
        tr_chi: np.ndarray,shape(nsta_loc,ncomp)
            travel-time chi
        am_chi: np.ndarray,shape(nsta_loc,ncomp)
            amplitude chi
        window_chi: np.ndarray,shape(nsta_loc,ncomp,20)
            window chi and adjoint source info  (20 values)
        """
        import os 
        from .utils import alloc_mpi_jobs 
        
        ncomp = tr_chi.shape[1]
        nsta_loc = len(tstart)

        # create directory
        if self.myrank == 0:
            os.makedirs(f"{self.MISFIT}/{self.mod}",exist_ok=True)

        # sync
        MPI.COMM_WORLD.Barrier()

        rcomp = [self.cc_comps[i][1] for i in range(self.ncomp)]
        if 'R' in rcomp or 'T' in rcomp:
            rcomp = ['R','T','Z']
        else:
            rcomp = ['E','N','Z']
        scomp = [self.cc_comps[i][0] for i in range(self.ncomp)]
        if 'R' in scomp or 'T' in scomp:
            scomp = ['R','T','Z']
        else:
            scomp = ['E','N','Z']

        # job id
        istart,_ = alloc_mpi_jobs(len(self.sta_names),self.nprocs,self.myrank)

        # loop each proc to print info on the screen
        # and save files
        for irank in range(self.nprocs):
            if irank == self.myrank:

                # open outfile
                outfile = f"{self.MISFIT}/{self.mod}/{self.evtid}_{bandname}_{self.meatype}_window_chi"
                if irank == 0:
                    fio = open(outfile,"w")
                else:
                    fio = open(outfile,"a")

                for ir in range(nsta_loc):
                    i = ir + istart
                    for ic in range(ncomp):
                        chs = self.cc_comps[ic][0]
                        chr = self.cc_comps[ic][1]
                        i_s = scomp.index(chs)
                        i_r = rcomp.index(chr)

                        # check if this station exists in this source component
                        if self.sta_names[i] not in self.sta_names_comps[i_s]:
                            continue

                        name = self.sta_names[i]
                        print(f'{name}.{self.chcode}{chr} CC comp = {self.cc_comps[ic]}')
                        print("Measurement window No.  1")
                        print("start and end time of window: %f %f" %(tstart[ir],tend[ir]) )
                        print(f"adjoint source and chi value for type = {self.adjsrc_type}")
                        print("%e" %(window_chi[ir,ic,6]))
                        print("tr_chi = %e am_chi = %e" %(tr_chi[ir,ic],am_chi[ir,ic]))
                        print("")

                        net,stnm = name.split('.')
                        fio.write(f"{self.evtid}_{chs} {net} {stnm} {self.chcode}{chr} 1 {self.adjsrc_type} ")
                        fio.write("%g %g " %(tstart[ir],tend[ir]))
                        for j in range(20):
                            fio.write("%g " %(window_chi[ir,ic,j]))
                        fio.write("%g %g 0. 0.\n" %(tr_chi[ir,ic],am_chi[ir,ic]))
                
                # close output file
                fio.close()

            # barrier
            MPI.COMM_WORLD.Barrier()

    def sum_adj_source(self):
        import os 
        import glob 

        nb = len(self.Tmax)
        if self.myrank == 0:
            print("\nsum adjoint source ...")
            for ic in range(len(self.scomp_syn)):
                os.makedirs(f"{self.syndirs[ic]}/SEM",exist_ok=True)
        MPI.COMM_WORLD.Barrier()

        # loop each syndir
        for ic in range(len(self.scomp_syn)):
            syndir = self.syndirs[ic]

            comps_read = ['E','N','Z']

            # find stations
            bandname = self._get_bandname(0)
            adj_files = glob.glob(f"{syndir}/OUTPUT_FILES/{bandname}/*.{self.chcode}Z.adj.sem.npy")
            nsta = len(adj_files)

            # loop each station
            for ir in range(nsta):
                name = adj_files[ir].split('/')[-1].split(f'.{self.chcode}Z.adj.sem.npy')[0]
                adj = np.zeros((3,self.npt_syn)) 

                for ib in range(nb):
                    bandname = self._get_bandname(ib) 
                    for ic in range(3):
                        data = np.load(f"{syndir}/OUTPUT_FILES/{bandname}/{name}.{self.chcode}{comps_read[ic]}.adj.sem.npy")
                        adj[ic,:] += data[:,1]
                
                # save to SEM/
                data = np.zeros((self.npt_syn,2)) 
                data[:,0] = self.t0_syn + self.dt_syn * np.arange(self.npt_syn)

                for ic,ch in enumerate(['E','N','Z']):
                    name1 = f"{name}.{self.chcode}{ch}.adj.sem.npy"
                    data[:,1] = adj[ic,:]
                    np.save(f"{syndir}/SEM/{name1}",data)
        
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

        # loop each syn dir
        for ic in range(len(self.scomp_syn)): 
            syndir = self.syndirs[ic]
        
            # first pack all obs/syn sacs to hdf5 
            for ib in range(len(self.Tmax)):
                bandname = self._get_bandname(ib)
                for tag in ['obs','syn']:

                    # all sac files
                    sacfiles = glob(f"{syndir}/OUTPUT_FILES/{bandname}/*.sac.{tag}")

                    # check if we need to save data
                    if len(sacfiles) == 0:
                        continue

                    # open h5file 
                    fio = h5py.File(f"{syndir}/OUTPUT_FILES/seismogram.{tag}.{bandname}.h5","w")
                    for i in range(len(sacfiles)):
                        tr = SACTrace.read(sacfiles[i])
                        if i == 0:
                            fio.attrs['dt'] = tr.delta  
                            fio.attrs['t0'] = tr.b 
                            fio.attrs['npts'] = tr.npts
                        
                        dsetname = tr.knetwk + "." + tr.kstnm + "." + tr.kcmpnm
                        fio.create_dataset(dsetname,shape=tr.data.shape,dtype='f4')
                        fio[dsetname][:] = tr.data 
                    
                    # close 
                    fio.close()
            
                # clean bandname
                shutil.rmtree(f"{syndir}/OUTPUT_FILES/{bandname}",ignore_errors=True)

            # clean semd and sem.ascii
            filenames = glob(f"{syndir}/OUTPUT_FILES/*.semd")
            for f in filenames:
                os.remove(f)
            filenames = glob(f"{syndir}/OUTPUT_FILES/*.sem.npy")
            for f in filenames:
                os.remove(f)
            filenames = glob(f"{syndir}/SEM/*.sem.npy")
            for f in filenames:
                os.remove(f)

    def execute(self):
        # rotate seismograms from XYZ to ZNE
        self._rotate_XYZ_to_ZNE()

        # rotate to RT
        self._rotate_seismogram_to_RT()

        # save current synthetics as observation if required
        if self.run_opt == 1:
            self.save_forward()
            return 0
    
        # loop each frequency band to compute misfit/adjoint source 
        for ib in range(len(self.Tmax)):
            self.cal_adj_source(ib)
            MPI.COMM_WORLD.Barrier()

        # sum adjoint source
        self.sum_adj_source()

        # pack SACs to h5
        if self.myrank ==0: self.cleanup()

        # sync
        MPI.COMM_WORLD.Barrier()
