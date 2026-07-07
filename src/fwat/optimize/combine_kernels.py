import numpy as np  
import h5py 
from mpi4py import MPI
import glob
import yaml

from fwat import FwatModel
from fwat.FortranIO import FortranIO
from fwat.const import PARAM_FILE


def gather_write_kls(irank:int,kl_opt:np.ndarray,gname_usr:list[str],KERNEL_DIR:str):
    # get search direction names
    nker = len(gname_usr)

    # write user model
    for i in range(nker):
        filename = f"{KERNEL_DIR}/{gname_usr[i]}.h5"
        if irank == 0:
            fh5 = h5py.File(filename,"w")
            dset = fh5.create_dataset(str(irank),data=kl_opt[i,:],dtype='f4')
            fh5.close()
            print("combine kernel %s ..."%(gname_usr[i]))
        else:
            fh5 = h5py.File(filename,"a")
            dset = fh5.create_dataset(str(irank),data=kl_opt[i,:],dtype='f4')
            fh5.close()

def run(argv):
    # mpi nprocs
    myrank = MPI.COMM_WORLD.Get_rank()
    nprocs = MPI.COMM_WORLD.Get_size()

    if len(argv) != 2:
        print("need 2 parameters: MODEL_DIR KERNEL_DIR")
        print("example: python combine_kernels.py M06 GRADIENT")
        exit(1)

    # get params
    MODEL_DIR = argv[0]
    KERNEL_DIR = argv[1]
    M = FwatModel()

    # get name list
    mname_list = M.model_names()
    gname_list = M.grad_names()
    nmod = len(mname_list)
    nkers = len(gname_list)

    # print sth
    if myrank == 0:
        print(f"combine kernels: MODEL_DIR = {MODEL_DIR} KERNEL_DIR = {KERNEL_DIR}")

    # loop each proc 
    # get size of this model
    filename = MODEL_DIR + "/proc%06d_" %(myrank) + mname_list[0] + ".bin"
    size = np.fromfile(filename,dtype='i4',count=1)[0] // 4

    # allocate space
    md = np.zeros((nmod,size))
    grad = np.zeros((nkers,size))

    # read base model
    for im in range(nmod):
        filename = MODEL_DIR + "/proc%06d_" %(myrank) + mname_list[im] + ".bin"
        fio = FortranIO(filename,"r")
        md[im,:] = fio.read_record('f4')
        fio.close()

    # read base kernels
    for im in range(nkers):
        filename = KERNEL_DIR + "/proc%06d_" %(myrank) + gname_list[im] + ".bin"
        fio = FortranIO(filename,"r")
        grad[im,:] = fio.read_record('f4')
        fio.close()

    # get search direction names
    gname_usr = M.grad_names(base=False)
    nkers_usr = len(gname_usr)
    kl_opt = np.zeros((nkers_usr + 1,size),dtype='f4') # another one for hess_kernel

    # read hess kernel
    filename = KERNEL_DIR + "/proc%06d_" %(myrank) + "hess_kernel.bin"
    fio = FortranIO(filename,"r")
    kl_opt[-1,:] = fio.read_record('f4')
    fio.close()
    
    # convert base kernels to optimization-space kernels
    _,kl_opt[:nkers_usr,:] = M.convert_kl(md,grad)

    # gather size of each proc
    size_list = MPI.COMM_WORLD.allgather(size)

    # output name list
    outnames = gname_usr + ["hess_kernel"]

    # write user model
    if myrank == 0:
        gather_write_kls(myrank, kl_opt, outnames, KERNEL_DIR)

        for irank in range(1,nprocs):
            size_rank = size_list[irank]
            kl_buf = np.zeros((nkers_usr + 1,size_rank),dtype='f4')
            MPI.COMM_WORLD.Recv(kl_buf,source=irank,tag=irank)
            gather_write_kls(irank, kl_buf, outnames, KERNEL_DIR)
    else:
        kl_buf = np.zeros((nkers_usr + 1,size),dtype='f4')
        kl_buf[:,:] = kl_opt[:,:]
        MPI.COMM_WORLD.Send(kl_buf, dest=0, tag=myrank)

    # sync
    MPI.COMM_WORLD.Barrier()

def main():
    import sys
    if len(sys.argv) !=5 :
        print("need 5 parameters: MODEL_DIR KERNEL_DIR model_type kernel_type")
        print("example: python write_user_model.py M06 GRADIENT mdtype kltype")
        exit(1)
    
    run(sys.argv[1:])

if __name__ == "__main__":
    main()