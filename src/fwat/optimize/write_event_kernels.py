import numpy as np  
import h5py 
from mpi4py import MPI
import glob

from fwat.optimize.model import FwatModel
from fwat.FortranIO import FortranIO

def run(argv):
    # mpi nprocs
    myrank = MPI.COMM_WORLD.Get_rank()
    nprocs = MPI.COMM_WORLD.Get_size()

    # get params
    MODEL_DIR = argv[0]
    KERNEL_DIR = argv[1]
    mdtype = argv[2]
    kltype = int(argv[3])

    # init FwatModel if required
    if mdtype == "iso":
        print("it's ISO model, vp,vs,rho are already existed")
        print("do nothing ...")
        exit(0)
    M = FwatModel(None,mdtype,kltype)

    # get name list
    mname_list = M.get_model_names()
    gname_list = M.get_grad_names()
    nmod = len(mname_list)

    # get how many files in the MODEL_DIR
    nfiles = len(glob.glob(f"{MODEL_DIR}/proc*_{mname_list[0]}.bin"))

    # print sth
    if myrank == 0:
        print("\nwrite kernels from cijkl ...")
        print(f"MODEL_DIR = {MODEL_DIR} KERNEL_DIR = {KERNEL_DIR}")
        print(f"nfiles = {nfiles}\n")

    # loop each proc 
    for irank in range(myrank,nfiles,nprocs):
        # get size of this model
        filename = MODEL_DIR + "/proc%06d_" %(irank) + mname_list[0] + ".bin"
        size = np.fromfile(filename,dtype='i4',count=1)[0] // 4

        # allocate space
        md = np.zeros((nmod,size))
        grad = np.zeros((nmod,size))

        # read base model
        for im in range(nmod):
            filename = MODEL_DIR + "/proc%06d_" %(irank) + mname_list[im] + ".bin"
            fio = FortranIO(filename,"r")
            md[im,:] = fio.read_record('f4')
            fio.close()

        # read base kernels
        for im in range(nmod):
            filename = f"{KERNEL_DIR}/{gname_list[im]}.h5"
            fh5 = h5py.File(filename,"r")
            grad[im,:] = fh5[str(irank)][:] * 1.
            fh5.close()
        
        # convert to plotting kernels
        _,direc = M.convert_kl(md,grad)
        direc = np.asarray(-direc,dtype='f4') # we need search direction

        # get search direction names
        dname_list = M.get_direc_names()
        nker = len(dname_list)

        # write user model
        for i in range(nker):
            filename = KERNEL_DIR + "/proc%06d_" %(irank) + dname_list[i] + ".bin"
            fio = FortranIO(filename,"w")
            fio.write_record(direc[i,...])
            fio.close()

def main():
    import sys
    if len(sys.argv) !=5 :
        print("need 5 parameters: MODEL_DIR KERNEL_DIR model_type kernel_type")
        print("example: python write_user_model.py M06 GRADIENT mdtype kltype")
        exit(1)
    
    run(sys.argv[1:])

if __name__ == "__main__":
    main()