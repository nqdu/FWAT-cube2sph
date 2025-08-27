import os 
from scipy.io import FortranFile 
import numpy as np 
from mpi4py import MPI
import h5py
import yaml
import glob 

from fwat.optimize.model import FwatModel
from fwat.const import OPT_DIR,SOLVER,PARAM_FILE

def compute_zpred_hess(iter_cur):

    # mpi rank/size
    comm = MPI.COMM_WORLD
    myrank = comm.Get_rank()

    # read external_mesh.bin for zstore/ibool
    filename = f'{OPT_DIR}/MODEL_M%02d'%(iter_cur) + '/proc%06d'%(myrank) + '_external_mesh.bin'
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

def run(argv):
    if len(argv) !=3 :
        print("Usage: fwat sum_kernel evtfile iter_cur MODEL ")
        print("example: fwat sum_kernel src_rec/sources.dat 0 M00.ls")
        exit(1)
    
    # get input 
    evtfile = argv[0]
    iter_cur = int(argv[1])
    MODEL = argv[2] 

    # get precond
    with open(f"{PARAM_FILE}","r") as fio:
        PRECOND = yaml.safe_load(fio)['optimize']['PRECOND_TYPE']

    assert(PRECOND in ['default','none','z_precond','z2_precond'])
    KERNEL_DIR = f'./{OPT_DIR}/SUM_KERNELS_{MODEL}'
    os.makedirs(KERNEL_DIR,exist_ok=True)
    
    # mpi rank/size
    comm = MPI.COMM_WORLD
    myrank = comm.Get_rank()

    # only do preconditioning for first iteration
    if iter_cur!= 0 and PRECOND == 'default' and 'ls' in MODEL:
        if myrank == 0:
            print(f"iteration = {iter_cur} > 0, skip preconditiononing")
        PRECOND = 'none' 

    # print info
    if myrank == 0:
        print(f"PRECOND_TYPE = {PRECOND} ")

    # init FwatModel
    M = FwatModel()

    # read sourcenames
    srctxt = np.loadtxt(evtfile,dtype=str,ndmin=2)
    nevts = srctxt.shape[0]

    # open weights_file
    weight = np.ones((nevts))
    if os.path.exists(f"{OPT_DIR}/weight_kl.txt"):
        if myrank == 0:
            print("find weight_kl.txt file, sum weighted kernels ... ")
        weight = np.loadtxt(f"{OPT_DIR}/weight_kl.txt")

    # sum kernels
    grad_list_base = M.get_grad_names(base=True)
    nmod = len(grad_list_base)
    for i in range(nmod):
        if myrank == 0: print(f'sum kernel {grad_list_base[i]}')
        kl = np.array([0])
        for ievt in range(nevts):
            # check if it's noise source
            filenames = glob.glob(f'./{SOLVER}/{MODEL}/{srctxt[ievt,0]}_[NEZRT]')
            filenames.append(f'./{SOLVER}/{MODEL}/{srctxt[ievt,0]}')

            # check if one of these files exists
            file_exists = False
            for j,f in enumerate(filenames):
                filename = f + '/GRADIENT/' + grad_list_base[i] + '.h5'
                if os.path.exists(filename):
                    file_exists = True 

            if not file_exists:
                if myrank == 0:
                    print(f"no kernel found for {srctxt[ievt,0]}, name = {grad_list_base[i]}")
                exit(1)
                
            for f in filenames:
                filename = f + '/GRADIENT/' + grad_list_base[i] + '.h5'
                if not os.path.exists(filename): continue

                # read 
                fio = h5py.File(filename,"r")
                arr = fio[str(myrank)][:]
                fio.close()

                # sum kernel
                kl = kl + arr * weight[ievt]
        
        # write out
        outname = KERNEL_DIR + "/proc%06d"%myrank + '_' + grad_list_base[i] + '.bin'
        f = FortranFile(outname,"w")
        kl = np.asarray(kl,dtype='f4')
        f.write_record(kl)
        f.close()

        if i == 0 and PRECOND == 'none':
            kl = kl * 0 + 1
            outname = KERNEL_DIR + "/proc%06d"%myrank + "_hess_kernel" + '.bin'
            kl = np.asarray(kl,dtype='f4')
            f = FortranFile(outname,"w")
            f.write_record(kl)
            f.close()
    
    if PRECOND == 'z_precond':
        outname = KERNEL_DIR + "/proc%06d"%myrank + "_hess_kernel" + '.bin'
        kl = compute_zpred_hess(iter_cur)
        f = FortranFile(outname,"w")
        kl = np.asarray(kl,dtype='f4')
        f.write_record(kl)
        f.close()

    if PRECOND == 'z2_precond':
        outname = KERNEL_DIR + "/proc%06d"%myrank + "_hess_kernel" + '.bin'
        kl = compute_zpred_hess(iter_cur)
        idx = np.logical_not(kl == 1.0e-8)
        kl[idx] = kl[idx]**2
        f = FortranFile(outname,"w")
        kl = np.asarray(kl,dtype='f4')
        f.write_record(kl)
        f.close()

    # sum hess if required
    if PRECOND == 'default':
        kl = np.array([0])
        for ievt in range(nevts):

            # check if it's noise source
            filenames = glob.glob(f'./{SOLVER}/{MODEL}/{srctxt[ievt,0]}_[NEZRT]')
            filenames.append(f'./{SOLVER}/{MODEL}/{srctxt[ievt,0]}')

            for f in filenames:
                filename = f + '/GRADIENT/hess_kernel.h5'
                if not os.path.exists(filename): continue

                # read 
                fio = h5py.File(filename,"r")
                arr = np.asarray(fio[str(myrank)][:],dtype='f4')
                fio.close()

                # sum kernel
                kl = kl + arr * weight[ievt]

                # filename = f'./{SOLVER}/{MODEL}' + setname + '/GRADIENT/' + 'hess_kernel.h5'

                # fio = h5py.File(filename,'r')
                # arr = fio[str(myrank)][:]
                # fio.close()
                
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
        kl = np.asarray(kl,dtype='f4')
        f.write_record(kl)
        f.close()


def main():
    import sys 
    if len(sys.argv) !=4 :
        print("need 3 parameters: evtfile iter_cur MODEL ")
        print("example: python sum_kernels.py src_rec/sources.dat 0 M00.ls")
        exit(1)

    # run
    run(sys.argv[1:])
    

if __name__ == "__main__":
    main()