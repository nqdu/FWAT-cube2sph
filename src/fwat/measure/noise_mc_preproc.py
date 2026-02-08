import numpy as np 
from mpi4py import MPI 
import os 
from fwat.const import PARAM_FILE
from obspy.io.sac import SACTrace
from fwat.measure.utils import cumtrapz1,dif1,bandpass \
                                ,interpolate_syn,alloc_mpi_jobs
from fwat.adjoint.MeasureStats import MeasureStats

def _get_snr(data:np.ndarray,win_b:int,win_e:int):
    """
    Get signal-to-noise ratio (SNR) from data within a specified window.

    Parameters
    -----------
    data: np.ndarray
        Input data array.
    win_b,win_e: int
        Start and end indices of the signal window.

    Returns
    -----------
    snr: float
        Signal-to-noise ratio, if data is very clean, return 1.0e20.
    """
    # compute max amplitude in signal window
    amp = np.max(abs(data[win_b:win_e]))

    # get noise std
    n = 0 
    std_in_noise = 1.0e-20 
    if win_b > 0:
        std_in_noise += np.std(data[0:win_b])
        n += 1
    if win_e < len(data) - 1:
        std_in_noise += np.std(data[win_e+1:])
        n += 1 
    if n > 0:
        std_in_noise /= n
    
    if std_in_noise <= 1.0e-20:
        snr = 1.0e20
    else: 
        snr = amp / std_in_noise

    return snr

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

def _get_window_by_group_vel(dist:float,vmin:float,vmax:float,Tmax:float,
                            npt:int,t0_inp:float,dt_inp:float):

    # start/end index of the normalized window
    win_b = np.floor((dist / vmax - Tmax / 2. + t0_inp) / dt_inp)
    win_e = np.floor((dist / vmin + Tmax / 2. - t0_inp) / dt_inp)
    win_b = int(max(win_b,0))
    win_e = int(min(win_e,npt-1))


    # compute time window
    tstart = dist / vmax - Tmax * 0.5 
    tend = dist / vmin + Tmax * 0.5 
    tstart = max(tstart,t0_inp)
    tend = min(tend,t0_inp+(npt-1)*dt_inp)

    return win_b,win_e,tstart,tend

class NoiseMC_PreOP():
    def __init__(self, measure_type:str, iter:int, evtid:str, run_opt:int):
        # import packages
        from .utils import read_params
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
        self.adjsrc_type:str = str(self.pdict['ADJSRC_TYPE'])

        # path to simulation directory
        mdir = "M%02d" %(iter)
        if self.run_opt == 2:
            mdir = mdir + ".ls"
        self.mod = mdir

        # get all synthetic directories
        self._get_mc_channel()

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
        t = np.array(fio[list(fio.keys())[0]])[:,0]
        fio.close()
        self.t0_syn = t[0] * 1
        self.dt_syn = t[1] - t[0]
        self.npt_syn = len(t)

        # seismograms 
        self.seismogram: dict[str, np.ndarray] = {}
        self.seismogram_adj: dict[str, np.ndarray] = {}

        # save sac here 
        self.seismogram_sac = {}
    
    def _get_mc_channel(self):
        import os 

        # get source comps and receiver comps
        self.cc_comps:list[str] = sorted(self.pdict['CC_COMPS'])
        self.ncomp = len(self.cc_comps)

        # sanity check
        self._sanity_check()

        # init source comps 
        comps_temp = set()

        # check available source components
        for comp in self.cc_comps:
            filename = f"{self.SRC_REC}/STATIONS_{self.evtid}_{comp}_globe"
            exist_flag = os.path.exists(filename)
            if comp[0] == 'Z' and exist_flag:
                comps_temp.add('Z')
            elif comp[0] in ['R','T'] and exist_flag:
                comps_temp.add('E')
                comps_temp.add('N')
            else:
                print(f"Source component {comp} file {filename} does not exist!")
                exit(1)
        scomp_syn = sorted(list(comps_temp))
        
        # add syndirs
        self.syndirs:list[str] = []
        for ic in range(len(scomp_syn)):
            syndir = f"{self.SOLVER}/{self.mod}/{self.evtid}_{scomp_syn[ic]}/"
            self.syndirs.append(syndir)

        # add source component 
        self.scomp_syn = scomp_syn

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
        if self.adjsrc_type not in ['5','7','exp_phase','cc_time','cc_time_dd']:
            if self.myrank == 0:
                print(f"only [5,7,'exp_phase','cc_time','cc_time_dd'] adjoint sources are supported!")
                print(f"adjsrc_type = {self.adjsrc_type}")
            exit(1)

    
    def _get_station_info(self):
        from .utils  import cal_dist_az_baz
        self.stainfo = {}
        self.sta_names_comp:list[set[str]] = []
        for ic in range(self.ncomp):
            comp = self.cc_comps[ic]
            stationfile = f'{self.SRC_REC}/STATIONS_{self.evtid}_{comp}_globe'
            statxt = np.loadtxt(stationfile,dtype=str,ndmin=2)
            nsta = statxt.shape[0]
            netwk = statxt[:,1]
            stnm = statxt[:,0]
            stla = statxt[:,2].astype(float)
            stlo = statxt[:,3].astype(float)

            # save station name in each component
            names = set()
            for ir in range(nsta):
                name = netwk[ir] + '.' + stnm[ir]
                if name not in self.stainfo:
                    dist,az,baz = cal_dist_az_baz(self.evla,self.evlo,stla[ir],stlo[ir])
                    self.stainfo[name] = [stlo[ir],stla[ir],az,baz,dist/1000]
                
                # save names for this component
                names.add(name)
            self.sta_names_comp.append(names)
        
        # get all names used for all multi-channels 
        self.sta_names: list[str] = sorted(self.stainfo.keys())

    def _rotate_XYZ_to_ZNE(self):
        from .cube2sph_rotate import rotate_seismo_fwd

        # get rotation list
        istart,iend = alloc_mpi_jobs(len(self.sta_names),self.nprocs,self.myrank)
        rot_list = [self.sta_names[i] for i in range(istart,iend+1)]

        # loop each synthetic directory
        for ic in range(len(self.scomp_syn)):
            fn_matrix = f"{self.SRC_REC}/rot_{self.evtid}"
            from_dir = f"{self.syndirs[ic]}/OUTPUT_FILES/"
            to_dir = f"{self.syndirs[ic]}/OUTPUT_FILES/"
            from_template='${nt}.${sta}.BX${comp}.semd'
            to_template='${nt}.${sta}.BX${comp}.sem.npy'

            # rotate seismograms from XYZ to ZNE
            out = rotate_seismo_fwd(rot_list, fn_matrix, from_dir, to_dir, from_template, to_template)
            self.seismogram.update(out)

    def _rotate_ZNE_to_XYZ(self):
        from .cube2sph_rotate import rotate_seismo_adj

        # get rotation list
        istart,iend = alloc_mpi_jobs(len(self.sta_names),self.nprocs,self.myrank)
        rot_list = [self.sta_names[i] for i in range(istart,iend+1)]

        # loop each synthetic directory
        for ic in range(len(self.scomp_syn)):
            fn_matrix = f"{self.SRC_REC}/rot_{self.evtid}"
            from_dir = f"{self.syndirs[ic]}/SEM/"
            to_dir = f"{self.syndirs[ic]}/SEM/"
            from_template='${nt}.${sta}.BX${comp}.adj.sem.npy'
            to_template='${nt}.${sta}.BX${comp}.adj'

            # rotate seismograms from XYZ to ZNE
            rotate_seismo_adj(rot_list, fn_matrix, from_dir, to_dir, from_template, to_template)

    def _get_bandname(self,ib:int):
        return 'T%03g_T%03g' %(self.Tmin[ib],self.Tmax[ib])

    def _rotate_seismogram_to_RT(self,):
        from .utils import rotate_EN_to_RT

        # get receiver components
        rcomp = [self.cc_comps[i][1] for i in range(self.ncomp)]

        # allocate jobs
        keys = self.sta_names
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
                    data[:,:] = self.seismogram[filename]
                    temp_syn[ic,:] = data[:,1]

                # rotate to R/T
                ve = temp_syn[0,:]
                vn = temp_syn[1,:]
                bazd = self.stainfo[keys[i]][3]
                vr,vt = rotate_EN_to_RT(ve,vn,bazd)
                temp_syn[0,:] = vr * 1.
                temp_syn[1,:] = vt * 1.

                # save to ascii
                for ic,ch in enumerate(['R','T']):
                    name = f"{keys[i]}.{self.chcode}{ch}.sem.npy"
                    filename = f"{self.syndirs[ie]}/OUTPUT_FILES/{name}"
                    data[:,1] = temp_syn[ic,:] 
                    #data.tofile(filename)
                    self.seismogram[filename] = data.copy()
                    #np.save(filename,data)

        # sync
        MPI.COMM_WORLD.Barrier()

    def _get_RTZRTZ_seismogram(self,i:int):
        data = np.zeros((3,3,self.npt_syn))

        # loop each component to load in data 
        for i_s,chs in enumerate(['E','N','Z']):
            syndir = f"{self.SOLVER}/{self.mod}/{self.evtid}_{chs}/"
            for i_r,chr in enumerate(['E','N','Z']):
                filename = f"{syndir}/OUTPUT_FILES/{self.sta_names[i]}.{self.chcode}{chr}.sem.npy"
                if filename not in self.seismogram:
                    continue
                data[i_s,i_r,:] = self.seismogram[filename][:,1] * 1.
        
        # rotate to RTZ-RTZ if needed
        azd,bazd = self.stainfo[self.sta_names[i]][2:4]
        R_s,R_r = _get_rotate_matrix(azd,bazd)
        data = np.einsum("ip,jq,pqk-> ijk",R_s,R_r,data)

        return data
    
    def _rotate_mc_adj(self,i:int,adj_src_all:np.ndarray):
        # rotate adjoint source to ENZ-ENZ
        azd,bazd = self.stainfo[self.sta_names[i]][2:4]
        R_s,R_r = _get_rotate_matrix(azd,bazd)
        adj_src_all = np.einsum("ip,jq,pqk-> ijk",R_s.T,R_r.T,adj_src_all)

        return adj_src_all
    
    def _cal_adj_by_name(self,name,dat_inp,syn_inp,
                         t0_inp,dt_inp,npt_cut,
                         ib:int,tstart:float,tend:float):
        
        if self.adjsrc_type == 'exp_phase':
            from fwat.adjoint.exp_phase_misfit import measure_adj_exphase
            return measure_adj_exphase(
                        dat_inp,syn_inp,
                        t0_inp,dt_inp,npt_cut,
                        self.Tmin[ib],self.Tmax[ib],
                        tstart,tend)
        elif self.adjsrc_type == 'cc_time':
            from fwat.adjoint.cc_misfit import measure_adj_cc
            return measure_adj_cc(
                        dat_inp,syn_inp,
                        t0_inp,dt_inp,npt_cut,
                        self.Tmin[ib],self.Tmax[ib],
                        tstart,tend)
        else:
            from .measure import measure_adj
            imeas = int(self.adjsrc_type)
            verbose = (self.myrank == 0) and (name == self.sta_names[0])
            return measure_adj(
                        t0_inp,dt_inp,npt_cut,
                        t0_inp,dt_inp,npt_cut,
                        tstart,tend,imeas,
                        self.Tmax[ib]*1.01,
                        self.Tmin[ib]*0.99,
                        verbose,dat_inp,
                        syn_inp)
    
    def _cal_adj_by_name_dd(self,dat_inp_i,dat_inp_j,
                            syn_inp_i,syn_inp_j,
                            t0_inp,dt_inp,npt_cut,
                            ib:int,tstart_i:float,tend_i:float,
                            tstart_j:float,tend_j:float):
        from fwat.adjoint.cc_misfit import measure_adj_cc_dd
        return measure_adj_cc_dd(
                    dat_inp_i,syn_inp_i,
                    dat_inp_j,syn_inp_j,
                    t0_inp,dt_inp,npt_cut,
                    self.Tmin[ib],self.Tmax[ib],
                    tstart_i,tend_i,tstart_j,tend_j)
        
    def save_forward(self):
        # get some vars
        dt_syn = self.dt_syn
        npt_syn = self.npt_syn
        t0_syn = self.t0_syn
        myrank = self.myrank
        
        # init a sac header
        tr = SACTrace(
            evla=self.evla,evlo=self.evlo,
            evdp=self.evdp,stla=0.,
            stlo=0.,stel=0,lcalda=True,
            delta = dt_syn,
            b=t0_syn
        )

        if self.myrank == 0:
            print("Synthetic Observations ...")

        # make directories
        if myrank == 0:
            for ic in range(self.ncomp):
                comp = self.cc_comps[ic]
                outdir = f"{self.DATA_DIR}/{self.evtid}_{comp[0]}"
                os.makedirs(outdir,exist_ok=True)
        MPI.COMM_WORLD.Barrier()
        
        # note, all data in seismograms are in ENZ-ENZ coordinates
        # we should rotate them to RTZ-RTZ if needed

        # allocate jobs 
        nsta = len(self.sta_names)
        istart,iend = alloc_mpi_jobs(nsta,self.nprocs,self.myrank)
        nsta_loc = iend - istart + 1

        # loop each stations 
        print_EGF = True
        for ir in range(nsta_loc):
            i = ir + istart 

            # fetch seismogram in RTZ-RTZ
            data = self._get_RTZRTZ_seismogram(i)

            # now loop each component to save sac files
            for ic in range(self.ncomp):
                # check if this component exists
                if self.sta_names[i] not in self.sta_names_comp[ic]:
                    continue

                # out dir
                outdir = f"{self.DATA_DIR}/{self.evtid}_{self.cc_comps[ic][0]}"

                # get chs/chr
                chs = self.cc_comps[ic][0]
                chr = self.cc_comps[ic][1]
                i_r = ['R','T','Z'].index(chr)
                i_s = ['R','T','Z'].index(chs)
                tr.data = data[i_s,i_r,:] * 1.

                # convert to EGF if required
                if self.pdict['USE_EGF'] == False:
                    if print_EGF and myrank == 0:
                        print("EGF => CCF ...")
                        print_EGF = False
                    tr.data = -cumtrapz1(tr.data,dt_syn)

                # save sac 
                code = f"{self.sta_names[i]}.{self.chcode}{chr}"
                tr.kcmpnm = f"{self.chcode}{chr}"
                tr.knetwk = self.sta_names[i].split('.')[0]
                tr.kstnm = self.sta_names[i].split('.')[1]
                info = self.stainfo[self.sta_names[i]]
                tr.stla = info[1]
                tr.stlo = info[0]
                filename = f"{outdir}/{code}.sac"
                tr.write(filename)

        # sync
        MPI.COMM_WORLD.Barrier()



    def cal_adj_source(self,ib:int):
        """ 
        Calculate adjoint source for noise cross-correlation measurement. 

        Parameters
        -----------
        ib: int
            Index of frequency band.
        """

        # check if double difference is enabled
        if '_dd' in self.adjsrc_type:
            self.cal_adj_source_dd(ib)
            return 

        # get vars
        npt_syn = self.npt_syn
        t0_syn = self.t0_syn
        dt_syn = self.dt_syn

        # new dt/nt for interpolate
        dt_inp = 0.01
        t0_inp = -10.
        if dt_inp > dt_syn:
            dt_inp = dt_syn * 1.
            t0_inp = t0_syn * 1.
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
        stats_list = []

        # temp 
        RTZ = ['R','T','Z']

        # loop each station
        for ir in range(nsta_loc):
            i = ir + istart 

            # fetch RTZRTZ seismogram
            data = self._get_RTZRTZ_seismogram(i)

            # allocate space for adjoint source
            adj_src_all = np.zeros((3,3,npt_syn))

            # loop each CC-component
            for ic in range(self.ncomp):
                # check if this component exists
                if self.sta_names[i] not in self.sta_names_comp[ic]:
                    continue
                
                # get chs/chr
                chs = self.cc_comps[ic][0]
                chr = self.cc_comps[ic][1]
                i_s = RTZ.index(chs)
                i_r = RTZ.index(chr)

                # load obs_data
                code = f"{self.sta_names[i]}.{self.chcode}{chr}"
                obs_tr = SACTrace.read(f"{self.DATA_DIR}/{self.evtid}_{chs}/{code}.sac")
                t0_obs = obs_tr.b * 1.
                dt_obs = obs_tr.delta * 1.
                npt_obs = obs_tr.npts * 1
                npt1_inp = int((npt_obs - 1) * dt_obs / dt_inp)
                dist = obs_tr.dist 

                # filter obs/syn in the band
                obs_tr.data = bandpass(obs_tr.data,dt_obs,freqmin,freqmax)
                syn_inp = bandpass(data[i_s,i_r,:],dt_syn,freqmin,freqmax)

                # now we will resample data to (t0_inp,dt_inp,npt_cut)
                # interp obs/syn
                dat_inp1 = interpolate_syn(obs_tr.data,t0_obs,dt_obs,npt_obs,
                                        t0_obs + dt_inp,dt_inp,npt1_inp)
                syn_inp = interpolate_syn(syn_inp,t0_syn,dt_syn,npt_syn,
                                         t0_inp,dt_inp,npt_cut)
                
                # compute time derivative 
                if self.pdict['USE_EGF'] == False:
                    dat_inp1 = -dif1(dat_inp1,dt_inp)
                    if self.myrank == 0: print("CCFs => EGFs ...")
        
                # cut 
                dat_inp = interpolate_syn(dat_inp1,t0_obs + dt_inp,dt_inp,npt1_inp,
                                         t0_inp,dt_inp,npt_cut)

                # find amplitude of the in the window to normalize
                dist = obs_tr.dist
                win_b,win_e,tstart,tend = _get_window_by_group_vel(dist,vmin_list[ib],vmax_list[ib],self.Tmax[ib],npt_cut,t0_inp,dt_inp)
                amp = np.max(np.abs(dat_inp[win_b:win_e]))

                # check snr
                snr = _get_snr(dat_inp,win_b,win_e)

                # normalize
                if amp == 0.: amp = 1.
                dat_inp *= np.max(np.abs(syn_inp[win_b:win_e])) / amp

                # compute misfits and adjoint source
                stats,adjsrc =  \
                    self._cal_adj_by_name(
                        self.sta_names[i],dat_inp,syn_inp,
                        t0_inp,dt_inp,npt_cut,
                        ib,tstart,tend
                    )
                adjsrc = interpolate_syn(adjsrc,t0_inp,dt_inp,npt_cut,
                                            t0_syn,dt_syn,npt_syn)

                # make sure the snr > snr_threshold
                if snr < snr_threshold:
                    adjsrc[:] = 0.
                    stats.misfit = 0.
                    stats.tr_chi = 0.
                    stats.am_chi = 0.
                    stats.tshift = 0.
                
                # append stats
                stats.code = f"{self.evtid}_{chs} {self.sta_names[i]}.{self.chcode}{chr}"
                stats_list.append(stats)

                # save adjoint source to adj_src_all
                adj_src_all[i_s,i_r,:] = adjsrc.copy()

                # save obs and syn data as sac
                outdir = f"{self.SOLVER}/{self.mod}/{self.evtid}_{chs}/OUTPUT_FILES/{bandname}"
                obs_tr.delta = dt_inp 
                obs_tr.b = t0_inp 
                obs_tr.data = dat_inp * 1.
                self.seismogram_sac[f"{outdir}/{code}.sac.obs"] = obs_tr.copy()
                #obs_tr.write(f"{outdir}/{code}.sac.obs")
                obs_tr.data = syn_inp * 1.
                self.seismogram_sac[f"{outdir}/{code}.sac.syn"] = obs_tr.copy()
                #obs_tr.write(f"{outdir}/{code}.sac.syn")
            # end for loop cc_comps

            # rotate adjoint source to ENZ-ENZ
            adj_src_all = self._rotate_mc_adj(i,adj_src_all)

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
                    self.seismogram_adj[filename] = data.copy()
                    #np.save(filename,data)

        # end for each station

        # print measure_adj information
        self._print_measure_info(bandname,stats_list)

    def _alloc_shared_memory(self,npt_cut:int):
        # allocate shared memory for syn/obs data 
        nsta = len(self.sta_names)
        ncomp = self.ncomp

        # 1. Create a communicator for processes on the SAME physical node
            # This segregates the world into groups based on shared memory availability.
        comm = MPI.COMM_WORLD
        comm_node = comm.Split_type(MPI.COMM_TYPE_SHARED)
        node_rank = comm_node.Get_rank()

        # Create Inter-node communicator (Leader group)
        # Only Rank 0 of every node joins this communicator
        leader_color = 0 if node_rank == 0 else MPI.UNDEFINED
        lead_comm = comm.Split(leader_color, key=self.myrank)

        # 2. Calculate memory size (Only Node-Rank 0 allocates the actual bytes)
        itemsize = MPI.DOUBLE.Get_size()
        if node_rank == 0:
            nbytes = nsta * 9 * npt_cut * itemsize
        else:
            nbytes = 0
        
        # check if comm_node is intracomm 
        if not isinstance(comm_node, MPI.Intracomm):
            print("Error: comm_node is not an Intracomm. Shared memory may not be supported on this platform.")
            exit(1)

        # 3. Allocate the Shared Memory Window
        sh_syn_win = MPI.Win.Allocate_shared(nbytes, itemsize, comm=comm_node)
        sh_obs_win = MPI.Win.Allocate_shared(nbytes, itemsize, comm=comm_node)
        
        # Wait for initialization to finish
        comm_node.Barrier()

        return sh_syn_win, sh_obs_win,lead_comm, comm_node


    def cal_adj_source_dd(self,ib:int):
        """ 
        Calculate adjoint source for noise cross-correlation measurement, double-difference version.

        Parameters
        -----------
        ib: int
            Index of frequency band.
        """

        # get vars
        npt_syn = self.npt_syn
        t0_syn = self.t0_syn
        dt_syn = self.dt_syn

        # new dt/nt for interpolate
        dt_inp = 0.01
        t0_inp = -10.
        if dt_inp > dt_syn:
            dt_inp = dt_syn * 1.
            t0_inp = t0_syn * 1.
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

        # allocate shared memory for syn/obs data 
        sh_syn_win, sh_obs_win, lead_comm, comm_node = self._alloc_shared_memory(npt_cut)
        buf_syn,_ = sh_syn_win.Shared_query(0)
        buf_obs,_ = sh_obs_win.Shared_query(0)
        syn_data = np.ndarray((nsta,3,3,npt_cut), dtype=np.float64, buffer=buf_syn)
        obs_data = np.ndarray((nsta,3,3,npt_cut), dtype=np.float64, buffer=buf_obs)
        if comm_node.Get_rank() == 0:
            syn_data.fill(0.)
            obs_data.fill(0.)
        comm_node.Barrier()
        
        # misfits
        ncomp = self.ncomp
        stats_list = []

        # temp 
        RTZ = ['R','T','Z']

        # loop to fetch data
        for ir in range(nsta_loc):
            i = ir + istart 

            # fetch RTZRTZ seismogram
            data = self._get_RTZRTZ_seismogram(i)

            # loop each CC-component
            for ic in range(self.ncomp):
                # check if this component exists
                if self.sta_names[i] not in self.sta_names_comp[ic]:
                    continue
                
                # get chs/chr
                chs = self.cc_comps[ic][0]
                chr = self.cc_comps[ic][1]
                i_s = RTZ.index(chs)
                i_r = RTZ.index(chr)

                # load obs_data
                code = f"{self.sta_names[i]}.{self.chcode}{chr}"
                obs_tr = SACTrace.read(f"{self.DATA_DIR}/{self.evtid}_{chs}/{code}.sac")
                t0_obs = obs_tr.b * 1.
                dt_obs = obs_tr.delta * 1.
                npt_obs = obs_tr.npts * 1
                npt1_inp = int((npt_obs - 1) * dt_obs / dt_inp)

                # filter obs/syn in the band
                obs_tr.data = bandpass(obs_tr.data,dt_obs,freqmin,freqmax)
                syn_inp = bandpass(data[i_s,i_r,:],dt_syn,freqmin,freqmax)

                # now we will resample data to (t0_inp,dt_inp,npt_cut)
                # interp obs/syn
                dat_inp1 = interpolate_syn(obs_tr.data,t0_obs,dt_obs,npt_obs,
                                        t0_obs + dt_inp,dt_inp,npt1_inp)
                syn_inp = interpolate_syn(syn_inp,t0_syn,dt_syn,npt_syn,
                                         t0_inp,dt_inp,npt_cut)
                
                # compute time derivative 
                if self.pdict['USE_EGF'] == False:
                    dat_inp1 = -dif1(dat_inp1,dt_inp)
                    if self.myrank == 0: print("CCFs => EGFs ...")
        
                # cut 
                dat_inp = interpolate_syn(dat_inp1,t0_obs + dt_inp,dt_inp,npt1_inp,
                                         t0_inp,dt_inp,npt_cut)

                # save to shared memory
                syn_data[i,i_s,i_r,:] = syn_inp.copy()
                obs_data[i,i_s,i_r,:] = dat_inp.copy()

                # save obs and syn data as sac
                outdir = f"{self.SOLVER}/{self.mod}/{self.evtid}_{chs}/OUTPUT_FILES/{bandname}"
                obs_tr.delta = dt_inp 
                obs_tr.b = t0_inp 
                obs_tr.data = dat_inp * 1.
                self.seismogram_sac[f"{outdir}/{code}.sac.obs"] = obs_tr.copy()
                obs_tr.data = syn_inp * 1.
                self.seismogram_sac[f"{outdir}/{code}.sac.syn"] = obs_tr.copy()

        # sync to make sure all data are ready
        MPI.COMM_WORLD.Barrier()
        if comm_node.Get_rank() == 0:
            lead_comm.Allreduce(MPI.IN_PLACE, syn_data, op=MPI.SUM)
            lead_comm.Allreduce(MPI.IN_PLACE, obs_data, op=MPI.SUM)
        comm_node.Barrier()

        # adjoint sources 
        adj_src_all = np.zeros((nsta,3,3,npt_syn),dtype=float)

        # loop each (i,j) pair to compute adjoint source for double-difference measurement
        njobs = nsta * (nsta - 1) // 2
        kstart,kend = alloc_mpi_jobs(njobs,self.nprocs,self.myrank)
        for k in range(kstart,kend+1):
            # get (i,j) pair from k
            i = int(nsta - 2 - np.floor(np.sqrt(-8*k + 4*nsta*(nsta-1) - 7) / 2.0 - 0.5))
            row_start_k = i * nsta - (i * (i + 1)) // 2
            j = int(k - row_start_k + i + 1)

            if i >= nsta or j >= nsta:
                print(f"Error in calculating (i,j) from k: i={i}, j={j}, k={k}")
                exit(1)

            # allocate space for adjoint source
            adj_src_i = np.zeros((3,3,npt_syn))
            adj_src_j = np.zeros((3,3,npt_syn)) 

            # loop each CC-component
            for ic in range(self.ncomp):
                # check if this component exists
                if self.sta_names[i] not in self.sta_names_comp[ic] or \
                    self.sta_names[j] not in self.sta_names_comp[ic]:
                    continue

                # get chs/chr
                chs = self.cc_comps[ic][0]
                chr = self.cc_comps[ic][1]
                i_s = RTZ.index(chs)
                i_r = RTZ.index(chr)

                # get distances
                dist_i = self.stainfo[self.sta_names[i]][4]
                dist_j = self.stainfo[self.sta_names[j]][4]

                # get syn/obs data from shared memory
                syn_inp_i = syn_data[i,i_s,i_r,:].copy()
                dat_inp_i = obs_data[i,i_s,i_r,:].copy()
                syn_inp_j = syn_data[j,i_s,i_r,:].copy()
                dat_inp_j = obs_data[j,i_s,i_r,:].copy()

                # normalize by amplitude in the window
                # find amplitude of the in the window to normalize
                win_b_i,win_e_i,tstart_i,tend_i = \
                        _get_window_by_group_vel(
                            dist_i,vmin_list[ib],vmax_list[ib],self.Tmax[ib],
                            npt_cut,t0_inp,dt_inp)
                amp_i = np.max(np.abs(dat_inp_i[win_b_i:win_e_i]))
                if amp_i == 0.: amp_i = 1.
                dat_inp_i *= np.max(np.abs(syn_inp_i[win_b_i:win_e_i])) / amp_i

                win_b_j,win_e_j,tstart_j,tend_j = \
                        _get_window_by_group_vel(
                            dist_j,vmin_list[ib],vmax_list[ib],self.Tmax[ib],
                            npt_cut,t0_inp,dt_inp)
                amp_j = np.max(np.abs(dat_inp_j[win_b_j:win_e_j]))
                if amp_j == 0.: amp_j = 1.
                dat_inp_j *= np.max(np.abs(syn_inp_j[win_b_j:win_e_j])) / amp_j

                # compute adjoint source for i and j
                stats,adjsrc_i,adjsrc_j =  \
                    self._cal_adj_by_name_dd(
                        dat_inp_i,dat_inp_j,
                        syn_inp_i,syn_inp_j,
                        t0_inp,dt_inp,npt_cut,
                        ib,tstart_i,tend_i,
                        tstart_j,tend_j
                    )
                adjsrc_i = interpolate_syn(adjsrc_i,t0_inp,dt_inp,npt_cut,
                                        t0_syn,dt_syn,npt_syn)
                adjsrc_j = interpolate_syn(adjsrc_j,t0_inp,dt_inp,npt_cut,
                                        t0_syn,dt_syn,npt_syn)

                # check snr
                snr_i = _get_snr(dat_inp_i,win_b_i,win_e_i)
                snr_j = _get_snr(dat_inp_j,win_b_j,win_e_j)
                if snr_i < snr_threshold or snr_j < snr_threshold:
                    adjsrc_i[:] = 0.
                    adjsrc_j[:] = 0.
                    stats.misfit = 0.
                    stats.tr_chi = 0.
                    stats.am_chi = 0.
                    stats.tshift = 0.
                
                # append stats
                stats.code = f"{self.evtid}_{chs} {self.sta_names[i]}.{self.chcode}{chr}-{self.sta_names[j]}.{self.chcode}{chr}"
                stats_list.append(stats)

                # save adjoint source
                adj_src_i[i_s,i_r,:] = adjsrc_i.copy()
                adj_src_j[i_s,i_r,:] = adjsrc_j.copy()
            
            # end for loop cc_comps

            # rotate adjoint source to ENZ-ENZ
            adj_src_i = self._rotate_mc_adj(i,adj_src_i)
            adj_src_j = self._rotate_mc_adj(j,adj_src_j)

            # accumulate adjoint source
            adj_src_all[i,:,:,:] += adj_src_i.copy()
            adj_src_all[j,:,:,:] += adj_src_j.copy()

        # sync to make sure all adjoint sources are ready
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, adj_src_all, op=MPI.SUM)

        # free space for shared memory
        if lead_comm != MPI.COMM_NULL:
            lead_comm.Free()
        if comm_node != MPI.COMM_NULL:
            comm_node.Free()
        sh_syn_win.Free()
        sh_obs_win.Free()

        # now loop each station to save adjoint sources and obs/syn data
        for ir in range(nsta_loc):
            i = ir + istart

            # save adjoint source
            rcomp = ['E','N','Z']
            data = np.zeros((npt_syn,2))
            data[:,0] = t0_syn + np.arange(npt_syn) * dt_syn
            for ic in range(len(self.scomp_syn)):
                i_s = rcomp.index(self.scomp_syn[ic])
                for i_r in range(3):
                    code = f"{self.sta_names[i]}.{self.chcode}{rcomp[i_r]}"
                    filename = f"{self.syndirs[ic]}/OUTPUT_FILES/{bandname}/{code}.adj.sem.npy"
                    data[:,1] = adj_src_all[i,i_s,i_r,:]
                    self.seismogram_adj[filename] = data.copy()

        # print measure_adj information
        self._print_measure_info(bandname,stats_list)

    def _print_measure_info(self,bandname:str,stats_list:list[MeasureStats]):
        """ 
        print measure_adj info and save to file 
        
        Parameters
        ----------
        bandname: str
            band name
        stats_list: list[MeasureStats]
            list of MeasureStats for each station and component
        """
        import os 
        from .utils import alloc_mpi_jobs 

        # create directory
        if self.myrank == 0:
            os.makedirs(f"{self.MISFIT}/{self.mod}",exist_ok=True)

        # sync
        MPI.COMM_WORLD.Barrier()

        # temp comps
        RTZ = ['R','T','Z']

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
                
                # loop each measure and print info
                nmeasure = len(stats_list)
                for im in range(nmeasure):
                    stats = stats_list[im]

                    # print info
                    print(f"Measuring stats of {stats.code}")
                    print("start and end time of window: %f %f" %(stats.tstart,stats.tend) )
                    print(f"adjoint source and chi value for adjoint type = {self.adjsrc_type}")
                    print("tshift = %e" %(stats.tshift))
                    print("tr_chi = %e am_chi = %e" %(stats.tr_chi,stats.am_chi))
                    print("")

                    # write to file
                    # format: evtid code adj_type tstart tend tshift tr_chi am_chi misfit
                    code_split = stats.code.split()
                    fio.write(f"{code_split[0]} {code_split[1]} {stats.adj_type} ")
                    fio.write("%g %g " %(stats.tstart,stats.tend))
                    fio.write("%g %g %g %g\n" %(stats.tshift,stats.tr_chi,stats.am_chi,stats.misfit))
                
                # close output file
                fio.close()

            # barrier
            MPI.COMM_WORLD.Barrier()
    
    def sum_adj_source(self):
        import os 
        import glob 
        from fwat.measure.utils import alloc_mpi_jobs

        nb = len(self.Tmax)
        if self.myrank == 0:
            print("\nsum adjoint source ...")
            for ic in range(len(self.scomp_syn)):
                os.makedirs(f"{self.syndirs[ic]}/SEM",exist_ok=True)
        MPI.COMM_WORLD.Barrier()

        # all stations
        nsta = len(self.sta_names)
        istart,iend = alloc_mpi_jobs(nsta,self.nprocs,self.myrank)
        nsta_loc = iend - istart + 1

        # loop each syndir
        for ic in range(len(self.scomp_syn)):
            syndir = self.syndirs[ic]

            comps_read = ['E','N','Z']

            # loop each station
            for ir in range(nsta_loc):
                i = ir + istart
                name = self.sta_names[i]
                adj = np.zeros((3,self.npt_syn)) 

                for ib in range(nb):
                    bandname = self._get_bandname(ib) 
                    for ic in range(3):
                        #data = np.load(f"{syndir}/OUTPUT_FILES/{bandname}/{name}.{self.chcode}{comps_read[ic]}.adj.sem.npy")
                        key = f"{syndir}/OUTPUT_FILES/{bandname}/{name}.{self.chcode}{comps_read[ic]}.adj.sem.npy"
                        adj[ic,:] += self.seismogram_adj[key][:,1]
                
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
        import re 
        if self.myrank == 0:
            print("cleaning up ...")

        # save obs/syn
        for ic in range(self.ncomp):
            comp = self.cc_comps[ic]
            chs = comp[0]
            chs_enz = ''
            syndir = f"{self.SOLVER}/{self.mod}/{self.evtid}_{chs}"

            if chs == 'Z':
                chs_enz = 'Z'
            elif chs == 'R':
                chs_enz = 'E'
            else: # chs == 'T'
                chs_enz = 'N'   
            syndir1 = f"{self.SOLVER}/{self.mod}/{self.evtid}" + f"_{chs_enz}"

            # pack all obs/syn sacs to hdf5
            for ib in range(len(self.Tmax)):
                bandname = self._get_bandname(ib)
                for tag in ['obs','syn']:
                    # all sac file
                    pattern = re.compile(rf"^{syndir}/OUTPUT_FILES/{bandname}/.*\.sac\.{tag}")
                    sacfiles = [s for s in self.seismogram_sac.keys() if pattern.match(s)]

                    # check if we need to save data
                    nfiles = len(sacfiles)
                    nfiles_tot = MPI.COMM_WORLD.allreduce(nfiles,op=MPI.SUM)
                    if nfiles_tot == 0:
                        continue

                    # loop each proc
                    for irank in range(self.nprocs):
                        # open h5file 
                        if irank == self.myrank:
                            # open file 
                            if self.myrank == 0:
                                fio = h5py.File(f"{syndir1}/OUTPUT_FILES/seismogram_{chs}.{tag}.{bandname}.h5","w")
                            else:
                                fio = h5py.File(f"{syndir1}/OUTPUT_FILES/seismogram_{chs}.{tag}.{bandname}.h5","a")
                            
                            # loop each sac file
                            for i in range(len(sacfiles)):
                                tr = self.seismogram_sac[sacfiles[i]]
                                if i == 0 and irank == 0:
                                    fio.attrs['dt'] = tr.delta  
                                    fio.attrs['t0'] = tr.b 
                                    fio.attrs['npts'] = tr.npts
                                
                                dsetname = tr.knetwk + "." + tr.kstnm + "." + tr.kcmpnm
                                dst = fio.create_dataset(dsetname,shape=tr.data.shape,dtype='f4')
                                dst[:] = np.array(tr.data) 
                        
                            # close 
                            fio.close()
                        MPI.COMM_WORLD.Barrier()


        # loop each syn dir
        for ic in range(len(self.scomp_syn)): 
            syndir = self.syndirs[ic]

            if self.myrank == 0:
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
        self.cleanup()

        # sync
        MPI.COMM_WORLD.Barrier()