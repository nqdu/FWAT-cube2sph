from os import stat
import numpy as np 
from mpi4py import MPI
from fwat.adjoint.MeasureStats import MeasureStats
from fwat.const import PARAM_FILE

def create_shared_communicators(global_comm):
    """
    Splits the global communicator into:
      1. node_comm: Processes on the same physical node (shared memory).
      2. leader_comm: A communicator of only the 'Rank 0' from each node.
    """
    
    # Create intra-node communicator (all processes on the same machine)
    node_comm = global_comm.Split_type(MPI.COMM_TYPE_SHARED)
    node_rank = node_comm.Get_rank()
    
    # Identify if this process is the 'leader' of its node
    is_leader = (node_rank == 0)
    
    # Create inter-node communicator (only leaders talk to each other)
    # color=0 for leaders, MPI.UNDEFINED for workers (who get None)
    leader_comm = global_comm.Split(color=0 if is_leader else MPI.UNDEFINED, 
                                    key=global_comm.Get_rank())
    
    return node_comm, leader_comm

def create_shared_array(node_comm,shape, dtype = float):
    """
    Allocates a shared memory window and wraps it in a NumPy array.
    """
    
    # Calculate bytes required for the array
    # Only the leader actually asks for memory; others ask for 0.
    is_leader = (node_comm.Get_rank() == 0)
    itemsize = np.dtype(dtype).itemsize
    nbytes = (np.prod(shape) * itemsize) if is_leader else 0
    
    # Allocate the shared window
    win = MPI.Win.Allocate_shared(nbytes, itemsize, comm=node_comm)
    
    # Query the memory location. 
    # rank=0 in Shared_query always points to the leader's allocated block.
    buf, _ = win.Shared_query(0)
    
    # Create the NumPy wrapper
    array = np.ndarray(buffer=buf, dtype=dtype, shape=shape)
    
    # Initialize memory to zero (Leader only)
    if is_leader:
        array.fill(0)
    
    # Ensure initialization is complete before any process returns
    node_comm.Barrier()
    
    return win, array


class FwatPreOP:
    def __init__(self,measure_type:str,iter:int,evtid:str,run_opt:int):
        # import packages
        from .utils import read_params,alloc_mpi_jobs
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
        
        # copy in input parameters
        assert(measure_type in ['sks','tele','noise',"rf"])
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
        pdict = read_params(f"{PARAM_FILE}")['measure'][self.meatype]
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

        # adjoint source type 
        self.adjsrc_type:str = str(self.pdict['ADJSRC_TYPE'])

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
        t = np.asarray(fio[list(fio.keys())[0]][:,0])
        fio.close()
        self.t0_syn = t[0]
        self.dt_syn = t[1] - t[0]
        self.npt_syn = len(t)
        
        # allocate jobs for each proc 
        istart,iend = alloc_mpi_jobs(self.nsta,self.nprocs,self.myrank)
        self.nsta_loc = iend - istart + 1
        self._istart = istart 

        # seismograms here, only save seismograms for local stations
        self.seismogram : dict[str,np.ndarray] = {}
        self.seismogram_adj : dict[str,np.ndarray] = {}
        self.seismo_sac = {}

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
        from .cube2sph_rotate import rotate_seismo_fwd

        # get rotation list
        rot_list = [self.netwk[i] + "." + self.stnm[i] for i in range(self._istart,self._istart + self.nsta_loc)]

        # set parameters
        fn_matrix = f"{self.SRC_REC}/rot_{self.evtid}"
        from_dir = f"{self.syndir}/OUTPUT_FILES/"
        to_dir = f"{self.syndir}/OUTPUT_FILES/"
        from_template='${nt}.${sta}.BX${comp}.semd'
        to_template='${nt}.${sta}.BX${comp}.sem.npy'

        # rotate seismograms from XYZ to ZNE
        self.seismogram =  \
            rotate_seismo_fwd(rot_list,fn_matrix,from_dir,to_dir,from_template,to_template)
    
    def _rotate_ZNE_to_XYZ(self):
        from .cube2sph_rotate import rotate_seismo_adj

        # get rotation list
        rot_list = [self.netwk[i] + "." + self.stnm[i] for i in range(self._istart,self._istart + self.nsta_loc)]

        # set parameters
        fn_matrix = f"{self.SRC_REC}/rot_{self.evtid}"
        from_dir = f"{self.syndir}/SEM/"
        to_dir = f"{self.syndir}/SEM/"
        from_template='${nt}.${sta}.BX${comp}.adj.sem.npy'
        to_template='${nt}.${sta}.BX${comp}.adj'

        # rotate seismograms from XYZ to ZNE
        rotate_seismo_adj(rot_list,fn_matrix,from_dir,to_dir,from_template,to_template)

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
        import sys 

        # create directory
        if self.myrank == 0:
            os.makedirs(f"{self.MISFIT}/{self.mod}",exist_ok=True)

        # sync
        comm = MPI.COMM_WORLD
        comm.Barrier()

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
                    fio.write(f"{self.evtid} {stats.code} {stats.adj_type} ")
                    fio.write("%g %g " %(stats.tstart,stats.tend))
                    fio.write("%g %g %g %g\n" %(stats.tshift,stats.tr_chi,stats.am_chi,stats.misfit))
                
                # close output file
                fio.close()

            # barrier
            comm.Barrier()

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
                    if f"{input_dir}/{name}" in self.seismogram_adj:
                        adj[ic,:] += self.seismogram_adj[f"{input_dir}/{name}"][:,1]
                
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
        import re 
        import os 
        import h5py 
        if self.myrank == 0:
            print("cleaning up ...")

        # first pack all obs/syn sacs to hdf5 
        for ib in range(len(self.Tmax)):
            bandname = self._get_bandname(ib)
            for tag in ['obs','syn']:

                # all sac files
                pattern = re.compile(rf"^{self.syndir}/OUTPUT_FILES/{bandname}/.*\.sac\.{tag}")
                sacfiles = [s for s in self.seismo_sac.keys() if pattern.match(s)]

                # check if we need to save data
                nfiles = len(sacfiles)
                nfiles_tot = MPI.COMM_WORLD.allreduce(nfiles,op=MPI.SUM)
                if nfiles_tot == 0:
                    continue

                # loop each proc
                for irank in range(self.nprocs):
                    if irank == self.myrank:

                        # open h5file 
                        if irank == 0:
                            fio = h5py.File(f"{self.syndir}/OUTPUT_FILES/seismogram.{tag}.{bandname}.h5","w")
                        else:
                            fio = h5py.File(f"{self.syndir}/OUTPUT_FILES/seismogram.{tag}.{bandname}.h5","a")
                        for i in range(len(sacfiles)):
                            tr = self.seismo_sac[sacfiles[i]]
                            if i == 0 and irank == 0:
                                fio.attrs['dt'] = tr.delta  
                                fio.attrs['t0'] = tr.b 
                                fio.attrs['npts'] = tr.npts
                            
                            dsetname = tr.knetwk + "." + tr.kstnm + "." + tr.kcmpnm
                            fio.create_dataset(dsetname,shape=tr.data.shape,dtype='f4')
                            fio[dsetname][:] = tr.data 
                        
                        # close 
                        fio.close()

                    MPI.COMM_WORLD.Barrier()



        # clean semd and sem.ascii
        if self.myrank == 0:
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
        from .utils import cal_dist_baz,rotate_EN_to_RT

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
                    data[:,:] = self.seismogram[filename]
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
                    self.seismogram[filename] = data.copy()
        MPI.COMM_WORLD.Barrier()

    def save_forward(self):
        pass

    def cal_adj_source(self,ib:int):
        pass

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

        # cleanup temporary files 
        self.cleanup()

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