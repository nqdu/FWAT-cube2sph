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

        # read FWAT params below
        # frequency band
        pdict = read_params(f"{PARAM_FILE}")['measure'][self.meatype]
        self.Tmin = [x[0] for x in pdict['FILTER_BANDS']]
        self.Tmax = [x[1] for x in pdict['FILTER_BANDS']]

        # 
        self.meatype = measure_type
        self.iter = iter 
        self.evtid = evtid
        self.run_opt = run_opt

        # get source comps and receiver comps
        self.cc_comps:list[str] = pdict['CC_COMPS']
        self.ncomp = len(self.cc_comps)

        # get all synthetic directories
        comps_temp = set()
        for ic in range(self.ncomp):
            if self.cc_comps[ic][0] == 'T' or self.cc_comps[ic][0] == 'R':
                comps_temp.add('N')
                comps_temp.add('E')
            else:
                comps_temp.add(self.cc_comps[ic][0])
        comps = sorted(list(comps_temp))
        self.scomp_syn = comps # only ZNE
        self.syndirs = []
        for ic in range(len(comps)):
            syndir = f"{self.SOLVER}/{mdir}/{self.evtid}_{comps[ic]}/"
            self.syndirs.append(syndir)

        # channel code
        self.chcode = pdict['CH_CODE']

        # backup pdict for furthur usage
        self.pdict = pdict

        # path to simulation directory
        mdir = "M%02d" %(iter)
        if self.run_opt == 2:
            mdir = mdir + ".ls"
        self.mod = mdir

        # read sourcer/receiver locations below
        # ---------------------
        # source loc
        sourcefile = f"{self.SRC_REC}/sources.dat.{measure_type}"
        self.evla,self.evlo,self.evdp = get_source_loc(f"{evtid}",sourcefile)

        # station info
        # read station coordinates
        self.stainfo = {}
        for ie in range(len(self.scomp)):
            chs = self.scomp[ie]
            stationfile = f'{self.SRC_REC}/STATIONS_{self.evtid}_{chs}_globe'
            statxt = np.loadtxt(stationfile,dtype=str,ndmin=2)
            nsta = statxt.shape[0]
            netwk = statxt[:,1]
            stnm = statxt[:,0]
            stla = statxt[:,2].astype(float)
            stlo = statxt[:,3].astype(float)

            for ir in range(nsta):
                name = netwk[ir] + '.' + stnm
                if name not in self.stainfo:
                    _,az,baz = cal_dist_az_baz(self.evla,self.evlo,stla[ir],stlo[ir])
                    self.stainfo[name] = [stlo[ir],stla[ir],az,baz]

        # read simulation info dt,t0,npts
        fio = h5py.File(f"{self.syndirs[0]}/OUTPUT_FILES/seismograms.h5","r")
        t = fio[list(fio.keys())[0]][:,0]
        fio.close()
        self.t0_syn = t[0]
        self.dt_syn = t[1] - t[0]
        self.npt_syn = len(t)

        # sanity check
        self._sanity_check()
    
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

    def _rotate_XYZ_to_ZNE(self):
        from .cube2sph_rotate import rotate_seismo_fwd

        # loop each synthetic directory
        for ie in range(len(self.scomp_syn)):
            ch = self.scomp_syn[ie]
            fn_matrix = f"{self.SRC_REC}/rot_{self.evtid}"
            from_dir = f"{self.syndirs[ie]}/OUTPUT_FILES/"
            to_dir = f"{self.syndirs[ie]}/OUTPUT_FILES/"
            from_template='${nt}.${sta}.BX${comp}.semd'
            to_template='${nt}.${sta}.BX${comp}.sem.npy'

            # rotate seismograms from XYZ to ZNE
            rotate_seismo_fwd(fn_matrix,from_dir,to_dir,from_template,to_template)

    def _rotate_ZNE_to_XYZ(self):
        from .cube2sph_rotate import rotate_seismo_adj

        # loop each synthetic directory
        for ie in range(len(self.scomp_syn)):
            ch = self.scomp_syn[ie]
            fn_matrix = f"{self.SRC_REC}/rot_{self.evtid}"
            from_dir = f"{self.syndirs[ie]}/SEM/"
            to_dir = f"{self.syndirs[ie]}/SEM/"
            from_template='${nt}.${sta}.BX${comp}.adj.sem.npy'
            to_template='${nt}.${sta}.BX${comp}.adj'

            # rotate seismograms from XYZ to ZNE
            rotate_seismo_adj(fn_matrix,from_dir,to_dir,from_template,to_template)


    def _print_measure_info(self,bandname,tstart,tend,tr_chi,am_chi,window_chi):
        import os 
        ncomp = tr_chi.shape[1]
        nsta_loc = self.nsta_loc
        imeas = self.pdict['IMEAS']

        # create directory
        if self.myrank == 0:
            os.makedirs(f"{self.MISFIT}/{self.mod}",exist_ok=True)

        # sync
        MPI.COMM_WORLD.Barrier()

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

    def _get_bandname(self,ib:int):
        return 'T%03g_T%03g' %(self.Tmin[ib],self.Tmax[ib])

    def sum_adj_source(self):
        import os 
        from .utils import rotate_RT_to_EN
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
                bandname = self._get_bandname(ib)
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
            bandname = self._get_bandname(ib)
            for tag in ['obs','syn']:

                # all sac files
                sacfiles = glob(f"{self.syndir}/OUTPUT_FILES/{bandname}/*.sac.{tag}")

                # check if we need to save data
                if len(sacfiles) == 0:
                    continue

                # open h5file 
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
                
                # close 
                fio.close()
        
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

    def _rotate_seismogram_to_RT(self,):
        from .utils import rotate_EN_to_RT
        from .utils import alloc_mpi_jobs

        # get receiver components
        rcomp = [self.cc_comps[i][1] for i in range(self.ncomp)]

        # loop each synthetic directory
        for ie in range(len(self.scomp_syn)):
            need_rt_rotate = False
            if 'R' in rcomp or 'T' in rcomp:
                need_rt_rotate = True 

            # check if need rotate
            if not need_rt_rotate: continue
            
            # load stations
            keys = sorted(self.stainfo.keys())
            nsta = len(keys)

            # allocate jobs 
            istart,iend = alloc_mpi_jobs(nsta,self.nprocs,self.myrank)

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
        from .utils import cumtrapz1,alloc_mpi_jobs
        from .utils import rotate_EN_to_RT

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
        for ie in range(self.ncomp):
            chs = self.cc_comps[ie][0]
            chr = self.cc_comps[ie][1]
            
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
                    filename = f"{self.syndirs[ie]}/OUTPUT_FILES/{code}.sem.npy"
                    data = np.load(filename)
                    tr.data = data[:,1] * 1.

                else: # chs is in ['R','T']
                    
                    # load E/N data
                    data = np.zeros((2,npt_syn))
                    for ic0,ch0 in enumerate(['N','E']):
                        code0 = f"{sta_names[i]}.{self.chcode}{ch0}"
                        syndir = f"{self.SOLVER}/{self.mod}/{self.evtid}_{ch0}/"
                        filename = f"{syndir}/OUTPUT_FILES/{code0}.sem.npy"
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
        from .measure import measure_adj
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

        # band name
        bandname = self._get_bandname(ib)
        freqmin = 1. / self.Tmax[ib]
        freqmax = 1. / self.Tmin[ib]

        # group velocity window
        vmin_list = [x[0] for x in self.pdict['GROUPVEL_WIN']]
        vmax_list = [x[1] for x in self.pdict['GROUPVEL_WIN']]
        imeas = self.pdict['IMEAS']

        # log
        if self.myrank == 0:
            print(f"preprocessing for band {bandname} ...")
            for ie in range(len(self.scomp_syn)):
                outdir = f"{self.syndirs[ie]}/OUTPUT_FILES/{bandname}"
                os.makedirs(outdir,exist_ok=True)
        MPI.COMM_WORLD.Barrier()

        # allocate jobs
        sta_names_all = sorted(self.stainfo.keys())
        nsta = len(sta_names_all)
        istart,iend = alloc_mpi_jobs(nsta,self.nprocs,self.myrank)
        nsta_loc = iend - istart 

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

        # load station names for each source comp
        sta_names_comps = []
        for ic in range(3):
            filename = f"{self.SRC_REC}/STATIONS_{self.evtid}_{scomp[ic]}_globe"
            if os.path.exists(filename):
                # load stations
                statxt = np.loadtxt(filename,dtype=str,ndmin=2)
                sta_names = set([statxt[i,1] + "." + statxt[i,0] for i in range(statxt.shape[0])])
            else:
                sta_names = set()
            sta_names_comps.append(sta_names)
        
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
                if sta_names_all[i] not in sta_names_comps[i_s]:
                    continue

                # load obs_data
                code = f"{sta_names_all[i]}.{self.chcode}{chr}"
                obs_tr = SACTrace.read(f"{self.DATA_DIR}/{self.evtid}_{chs}/{code}.sac")
                t0_obs = obs_tr.b 
                dt_obs = obs_tr.delta 
                npt_obs = obs_tr.npts
                npt1_inp = int((npt_obs - 1) * dt_obs / dt_inp)
                dist = obs_tr.dist 

                # load synthetic data
                if chs in ['N','E','Z']:
                    filename = f"{self.syndirs[ic]}/OUTPUT_FILES/{code}.sem.npy"
                    syn_tr = np.load(filename)[:,1]
                else: # chs is in ['R','T']
                    # load E/N data
                    data = np.zeros((2,npt_syn))
                    for ic0,ch0 in enumerate(['N','E']):
                        code0 = f"{sta_names[i]}.{self.chcode}{ch0}"
                        syndir = f"{self.SOLVER}/{self.mod}/{self.evtid}_{ch0}/"
                        filename = f"{syndir}/OUTPUT_FILES/{code0}.sem.npy"
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
                verbose = (self.myrank == 0) and (ir == 0)
                tr_chi[ir,ic],am_chi[ir,ic],win_chi[ir,ic,:],adjsrc =   \
                    measure_adj(t0_inp,dt_inp,npt_cut,t0_syn,dt_syn,npt_syn,
                                tstart[ir],tend[ir],imeas,self.Tmax[ib]*1.01,
                                self.Tmin[ib]*0.99,verbose,dat_inp,
                                syn_inp)

                # save adjoint source to adj_src_all
                adj_src_all[i_s,i_r,:] = adjsrc.copy()

                # save obs and syn data as sac
                outdir = f"{self.syndirs[i_s]}/OUTPUT_FILES/{bandname}"
                obs_tr.delta = dt_inp 
                obs_tr.b = t0_inp 
                obs_tr.data = dat_inp 
                obs_tr.write(f"{outdir}/{code}.sac.obs")
                obs_tr.data = syn_inp
                obs_tr.write(f"{outdir}/{code}.sac.syn")

        # rotate adjoint source to ENZ-ENZ
        azd,bazd = self.stainfo[sta_names_all[i]]
        rotate_src = (scomp == ['R','T','Z'])
        rotate_sta = (rcomp == ['R','T','Z'])
        R_s,R_r = _get_rotate_matrix(azd,bazd,rotate_src=rotate_src,rotate_sta=rotate_sta)
        adj_src_all = np.einsum("ip,jq,pqk-> ijk",R_s,R_r,adj_src_all)

        # save adjoint source
        data = np.zeros((npt_syn,2))
        data[:,0] = t0_syn + np.arange(npt_syn) * dt_syn
        for ic in range(len(self.scomp_syn)):
            for ic1 in range(len(self.scomp_syn)):
                code = f"{sta_names_all[i]}.{self.chcode}{self.scomp_syn[ic1]}"
                filename = f"{self.syndirs[ic]}/OUTPUT_FILES/{bandname}/{code}.adj.sem.npy"
                data[:,1] = adj_src_all[ic,ic1,:]
                np.save(filename,data)

        # print measure_adj information


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
