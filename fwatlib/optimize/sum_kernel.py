import sys 
import os 
from scipy.io import FortranFile 
import numpy as np 
from tools import *
from mpi4py import MPI
import h5py

def compute_zpred_hess(iter_cur):

    # mpi rank/size
    comm = MPI.COMM_WORLD
    myrank = comm.Get_rank()

    # read external_mesh.bin for zstore/ibool
    filename = './optimize/MODEL_M%02d'%(iter_cur) + '/proc%06d'%(myrank) + '_external_mesh.bin'
    f = FortranFile(filename)
    _ = f.read_ints('i4')[0] # nspec
    _ = f.read_ints('i4')[0] #gnlob
    _ = f.read_ints('i4')
    ibool = f.read_ints('i4') - 1
    xstore = f.read_reals('f4')
    ystore = f.read_reals('f4')
    zstore = f.read_reals('f4')
    f.close()

    # compute depth in spherical coordinates
    EARTH = 6371000
    rstore = np.sqrt((xstore/EARTH)**2 + (ystore/EARTH)**2 + (zstore/EARTH)**2) * EARTH
    zstore = rstore - EARTH

    # find max
    dmax_loc = np.max(abs(zstore))
    dmax = comm.allreduce(dmax_loc,MPI.MAX)

    # compute kernel
    zl = zstore[ibool]
    hess = zl * 0
    idx = zl >= 0.
    hess[idx] = 1.0e-8
    idx = np.logical_not(idx)
    hess[idx] = abs(zl[idx]) / dmax

    # get hess min/max
    maxh_loc = np.max(abs(hess))
    maxh = comm.allreduce(maxh_loc,MPI.MAX)
    if myrank == 0:
        print(f"hess maximum = {maxh}")

    return hess

def main():
    if len(sys.argv) !=4 :
        print("need 4 parameters: iter evtfile PRECOND_TYPE [zmin zmax]")
        print("example: python sum_kernels.py 0 src_rec/sources.dat [default,none,z_precond,z2_precond] ")
        exit(1)
    
    # get input 
    iter_cur = int(sys.argv[1])
    evtfile = sys.argv[2]
    PRECOND = sys.argv[3]
    assert(PRECOND in ['default','none','z_precond','z2_precond'])
    KERNEL_DIR = './optimize/SUM_KERNELS_M%02d'%(iter_cur)
    os.makedirs(KERNEL_DIR,exist_ok=True)
    
    # mpi rank/size
    comm = MPI.COMM_WORLD
    myrank = comm.Get_rank()

    # only do preconditioning for first iteration
    if iter_cur!= 0 and PRECOND == 'default':
        if myrank == 0:
            print(f"iteration = {iter_cur} > 0, skip preconditiononing")
        PRECOND = 'none' 

    # print info
    if myrank == 0:
        print(f"PRECOND_TYPE = {PRECOND} ")

    # read sourcenames
    srctxt = np.loadtxt(evtfile,dtype=str,ndmin=2)
    nevts = srctxt.shape[0]

    # open weights_file
    weight = np.ones((nevts))
    if os.path.exists("optimize/weight_kl.txt"):
        if myrank == 0:
            print("find weight_kl.txt file, sum weighted kernels ... ")
        weight = np.loadtxt("optimize/weight_kl.txt")

    # sum kernels
    grad_list = get_gradname_list()
    nker = len(grad_list)
    for i in range(nker):
        if myrank == 0: print(f'sum kernel {grad_list[i]}')
        kl = np.float32(0.)
        for ievt in range(nevts):
            # read kernel from h5
            setname = '/' + srctxt[ievt,0]
            filename = './solver/M%02d'%(iter_cur) + setname + '/GRADIENT/' + grad_list[i] + '.h5'
            fio = h5py.File(filename,"r")
            arr = fio[str(myrank)][:]
            fio.close()

            # sum kernel
            kl = kl + arr * weight[ievt]
        
        outname = KERNEL_DIR + "/proc%06d"%myrank + '_' + grad_list[i] + '.bin'
        f = FortranFile(outname,"w")
        kl = np.float32(kl)
        f.write_record(kl)
        f.close()

        if i == 0 and PRECOND == 'none':
            kl = np.float32(kl * 0 + 1.)
            outname = KERNEL_DIR + "/proc%06d"%myrank + "_hess_kernel" + '.bin'
            f = FortranFile(outname,"w")
            f.write_record(kl)
            f.close()
    
    if PRECOND == 'z_precond':
        outname = KERNEL_DIR + "/proc%06d"%myrank + "_hess_kernel" + '.bin'
        kl = compute_zpred_hess(iter_cur)
        f = FortranFile(outname,"w")
        kl = np.float32(kl)
        f.write_record(kl)
        f.close()

    if PRECOND == 'z2_precond':
        outname = KERNEL_DIR + "/proc%06d"%myrank + "_hess_kernel" + '.bin'
        kl = compute_zpred_hess(iter_cur)
        idx = np.logical_not(kl == 1.0e-8)
        kl[idx] = kl[idx]**2
        f = FortranFile(outname,"w")
        kl = np.float32(kl)
        f.write_record(kl)
        f.close()

    # sum hess if required
    if PRECOND == 'default':
        kl = np.float32(0.)
        for ievt in range(nevts):
            setname = '/' + srctxt[ievt,0]
            filename = './solver/M%02d'%(iter_cur) + setname + '/GRADIENT/' + 'hess_kernel.h5'

            fio = h5py.File(filename,'r')
            arr = fio[str(myrank)][:]
            fio.close()
                
            # get hessian norm
            s = np.sum(arr * arr)
            s_all = comm.allreduce(s,MPI.SUM)
            if myrank == 0:
                print(f'event {ievt + 1} of {nevts}, hessian norm = {np.sqrt(s_all)}')
            kl = kl + np.abs(arr)
        
        outname = KERNEL_DIR + "/proc%06d"%myrank + "_hess_kernel" + '.bin'
        
        # normalize kl 
        maxh_loc = np.max(np.abs(kl))
        maxh = comm.allreduce(maxh_loc,MPI.MAX)
        if myrank == 0:
            print(f'hessian maximum : {maxh}')
        if maxh < 1.0e-18:
            kl = kl * 0 + 1.
        else:
            kl = kl / maxh 
        
        # inverse 
        THRESHOLD_HESS = 1.0e-3
        idx = kl > THRESHOLD_HESS
        idx1 = np.logical_not(idx)
        kl[idx] = 1. / kl[idx]
        kl[idx1] = 1. / THRESHOLD_HESS

        f = FortranFile(outname,"w")
        kl = np.float32(kl)
        f.write_record(kl)
        f.close()

main()