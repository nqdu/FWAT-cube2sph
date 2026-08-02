import os
from sys import argv
import numpy as np 
from mpi4py import MPI
import h5py
import yaml
import glob 

from fwat import FwatModel
from fwat.const import OPT_DIR,SOLVER,PARAM_FILE,SRC_REC,NGLL
from fwat.FortranIO import FortranIO
from fwat.measure.cal_misfit import compute_misfit
from fwat.optimize.search_direction import compute_inner_dot
from fwat.optimize.libgll import get_gll_weights

def compute_zpred_hess(iter_cur):

    # mpi rank/size
    comm = MPI.COMM_WORLD
    myrank = comm.Get_rank()

    # read external_mesh.bin for zstore/ibool
    filename = f'{OPT_DIR}/MODEL_M%02d'%(iter_cur) + '/proc%06d'%(myrank) + '_external_mesh.bin'
    f = FortranIO(filename,"r")
    _ = f.read_record('i4')[0] # nspec
    _ = f.read_record('i4')[0] #gnlob
    _ = f.read_record('i4')
    ibool = f.read_record('i4') - 1
    xstore = f.read_record('f4')
    ystore = f.read_record('f4')
    zstore = f.read_record('f4')
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

def _get_database_size(MODEL:str) -> int:

    comm = MPI.COMM_WORLD
    myrank = comm.Get_rank()

    filename = f'{OPT_DIR}/MODEL_{MODEL}/proc%06d'%(myrank) + '_external_mesh.bin'
    f = FortranIO(filename,"r")
    nspec = f.read_record('i4')[0]
    f.close()

    return nspec * NGLL**3

def get_summed_kernel(MODEL:str,simu_type:str) -> np.ndarray:

    # mpi rank/size
    comm = MPI.COMM_WORLD
    myrank = comm.Get_rank()

    # init FwatModel
    M = FwatModel()
    grad_list_user = M.grad_names(base=False)

    # init grad_base
    nkers = len(grad_list_user)
    ksize = _get_database_size(MODEL)
    grad_user = np.zeros((nkers,ksize),dtype=float)

    # read source list
    srctxt = np.loadtxt(f'{SRC_REC}/sources.dat.{simu_type}',dtype=str,ndmin=2,usecols=[0])
    nevts = srctxt.shape[0]

    # loop all events and sum kernels
    for ievt in range(nevts):
        # check if it's noise source
        if simu_type == 'noise':
            filenames = glob.glob(f'./{SOLVER}/{MODEL}/{srctxt[ievt,0]}_[NEZRT]')
        else:
            filenames = []
        filenames.append(f'./{SOLVER}/{MODEL}/{srctxt[ievt,0]}')
        if myrank == 0:
            print(f"sum kernels for event {srctxt[ievt,0]} in simu type {simu_type}, filenames = {filenames}")
            print("nkernels = %d, ksize = %d" %(nkers,ksize))

        for i in range(nkers):                
            for f in filenames:
                filename = f + '/GRADIENT/' + grad_list_user[i] + '.h5'
                if not os.path.exists(filename): continue

                # read 
                fio = h5py.File(filename,"r")
                #arr = np.array(fio[str(myrank)])[:]
                arr = np.array(fio[str(myrank)][:])
                fio.close()

                # sum kernel
                grad_user[i,:] += arr

    return grad_user

def _get_grad_norm(MODEL:str,SIMU_TYPES:list[str],grad_user:np.ndarray) -> np.ndarray:
    """
    Compute the norm of the gradient for each simulation type.

    Parameters
    ----------
    MODEL : str
        current model name, e.g., M00.ls/M00
    SIMU_TYPES : list[str]
        List of simulation types, shape (nsims,).
    grad_user : np.ndarray
        User-defined gradient for this iteration, shape (nsims,nkers,ksize).

    Returns
    -------
    gnorm : np.ndarray
        Norm of the gradient for each simulation type, shape (nsims,)
    """
    # mpi rank
    comm = MPI.COMM_WORLD
    myrank = comm.Get_rank()

    # get 3-D GLL weights
    w = get_gll_weights()
    NGLL3 = NGLL**3
    wgll3d = np.zeros((NGLL3),dtype='f4')
    for k in range(NGLL):
        for j in range(NGLL):
            for i in range(NGLL):
                wgll3d[k*NGLL*NGLL+j*NGLL+i] = w[i] * w[j] * w[k]
    
    # read jacobian
    # read nspec/nglob/ibool
    iter_cur = int(MODEL.split('.')[0][1:])
    filename = f'{OPT_DIR}/MODEL_M%02d'%(iter_cur) + '/proc%06d'%(myrank) + '_external_mesh.bin'
    f = FortranIO(filename,"r")
    nspec = f.read_record('i4')[0]
    _ = f.read_record('i4')[0]
    f.read_record('i4')
    _= f.read_record('i4')
    for _ in range(3):
        f.read_record('f4')
    f.read_record('i4') # irregular_element_number
    f.read_record('f4')
    f.read_record('f4')
    for _ in range(9):
        f.read_record('f4')
    jaco = f.read_record('f4').reshape(nspec,NGLL3)
    f.close()

    # compute norm of the gradient for each simulation type, shape (nsims,)
    nsims = grad_user.shape[0]
    gnorm = np.zeros(nsims)
    for i in range(nsims):
        # compute |g|
        nkers = grad_user.shape[1]
        g0 = 0.
        for j in range(nkers):
            grad = grad_user[i,j,:].reshape(nspec,NGLL3)
            g0 += compute_inner_dot(grad,grad,jaco,wgll3d)
        g0_all = comm.allreduce(g0,MPI.SUM)
        gnorm[i] = np.sqrt(g0_all)

        # read events for this simu type
        if myrank == 0:
            s = SIMU_TYPES[i]
            print(f"MODEL {MODEL}: Original |g| for simu type {s} = {gnorm[i]}")
    
    return gnorm

def _get_normalized_weights(param:dict,MODEL:str, 
                        SIMU_TYPES:list[str], iter_cur:int,
                        gnorm:np.ndarray) -> tuple[np.ndarray,bool]:
    """
    compute adaptive weights for each simulation type based on the misfit of this iteration, only for multiple simu types, e.g., ["noise","tele"]
    chi = sum_i L_i(x) / L_i(x0) / |g_i(x0)| * w_i, and weights will be re-written as (w_i / L_i(x0) / |g_i(x0)|)

    Parameters
    ----------
    param : dict
        fwat parameters read from fwat.yaml
    MODEL : str
        current model name, e.g., M00.ls
    SIMU_TYPES : list[str]
        list of simulation types, e.g., ["noise","tele"]
    iter_cur : int
        current iteration id
    grad_user : np.ndarray
        user defined gradient for this iteration, shape (nsims,nkers,ksize)

    Returns
    -------
    weights_type : np.ndarray
        adaptive weights for each simulation type, shape (len(SIMU_TYPES),)
    success : bool
        whether the adaptive weights were successfully computed, if False, the caller should use the original weights from fwat.yaml
      
    """
    # mpi info 
    comm = MPI.COMM_WORLD
    myrank = comm.Get_rank()

    # sanity check
    nsim = len(SIMU_TYPES)
    iter_wts = param['simulation']['iter_wts']
    if iter_wts > iter_cur:
        if myrank == 0:
            print(f"ERROR! iter_wts = {iter_wts} > iter_cur = {iter_cur}, STOP PROGRAM!!!!")
        exit(1)
    if nsim == 1:
        weights_usr = np.ones((nsim))
        return weights_usr, False
    
    # read weights in parameter file
    weights_bak = np.array(param['simulation']['weights'])

    # find if we have weights.txt in OPT_DIR
    if os.path.exists(f'{OPT_DIR}/weights.txt'):
        weights_usr = np.loadtxt(f'{OPT_DIR}/weights.txt')
        # check if weights_usr has the same length as nsim
        if len(weights_usr) != nsim:
            if myrank == 0:
                print(f"ERROR! weights.txt found but length {len(weights_usr)} != nsim {nsim}, \
                     ignore weights.txt and use original weights from fwat.yaml")
    else:
        # not found
        weights_usr = weights_bak.copy()


    # return current weights_usr
    is_ls_model = 'ls' in MODEL
    if iter_wts < 0 or nsim == 1 or iter_wts < iter_cur or is_ls_model :
        return weights_usr, False

    # in this case iter_wts == iter_cur, we need to compute adaptive weights based on misfit of this iteration
    # reset weights_usr to that in parameter file, weights will be re-computed based on misfit of this iteration
    weights_usr = weights_bak.copy()

    # normalize type 
    norm_type = param['simulation']['norm_type']
    assert norm_type in ['gradient','misfit'], f"norm_type = {norm_type} not supported, should be gradient or misfit"
    if myrank == 0:
        print(f"\ncompute adaptive weights based on norm_type = {norm_type}")

    # init
    chi = np.ones(nsim)
    gnorminv = np.ones(nsim)

    if norm_type == 'misfit':
        # compute misfit for each simulation type
        for i,s in enumerate(SIMU_TYPES):
            chi[i],_ = compute_misfit(MODEL,simu_type=s)
    else:
        # normalize gnorm
        gnorm_max = np.max(gnorm)
        gnorminv = 1. / (gnorm / gnorm_max)
        gnorminv = gnorminv / np.sum(gnorminv)

    # update weights
    weights_type = weights_usr * gnorminv / chi 
    success = True

    return weights_type, success

def compute_hessian_kernel(iter_cur:int, SIMU_TYPES:list[str], weights_type:np.ndarray) -> np.ndarray:
    """
    compute hessian kernel by summing kernels from different simulation types, only for multiple simu types, e.g., ["noise","tele"]
    kl = sum_i w_i * |k_i|, where k_i is the kernel for simu type i, and w_i is the weight for simu type i

    Parameters
    ----------
    iter_cur : int
        current iteration id
    SIMU_TYPES : list[str]
        list of simulation types, e.g., ["noise","tele"]
    weights_type : np.ndarray
        adaptive weights for each simulation type, shape (len(SIMU_TYPES),)

    Returns
    -------
    kl : np.ndarray
        hessian kernel, shape (ksize,)
      
    """
    comm = MPI.COMM_WORLD
    myrank = comm.Get_rank()

    kl = np.array([0])
    for i in range(len(SIMU_TYPES)):
        srctxt = np.loadtxt(f'{SRC_REC}/sources.dat.{SIMU_TYPES[i]}',dtype=str,ndmin=2,usecols=[0])
        nevts = srctxt.shape[0]
        for ievt in range(nevts):
            # check if it's noise source
            if SIMU_TYPES[i] == 'noise':
                filenames = glob.glob(f'./{SOLVER}/M%02d'%(iter_cur) + f'/{srctxt[ievt,0]}_[NEZRT]')
            else:
                filenames = []
            filenames.append(f'./{SOLVER}/M%02d'%(iter_cur) + f'/{srctxt[ievt,0]}')

            for f in filenames:
                filename = f + '/GRADIENT/hess_kernel.h5'
                if not os.path.exists(filename): continue

                # read 
                fio = h5py.File(filename,"r")
                arr = np.array(fio[str(myrank)])[:].astype('f4')
                fio.close()
                
                # get hessian norm
                s = np.sum(arr * arr)
                s_all = comm.allreduce(s,MPI.SUM)
                if myrank == 0:
                    print(f'event {srctxt[ievt,0]}, hessian norm = {np.sqrt(s_all)}')
                kl = kl + np.abs(arr) * weights_type[i]

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

    return kl

def compute_preconditioner(iter_cur:int,PRECOND:str, SIMU_TYPES:list[str], weights_type:np.ndarray) -> np.ndarray:
    # this is a placeholder function for computing preconditioner, you can replace it with your own function, e.g., compute_zpred_hess or compute_hessian_kernel

    # base preconditioner on zpred_hess
    kl = compute_zpred_hess(iter_cur)
    if PRECOND == 'none':
        kl = kl * 0 + 1.
    elif PRECOND == 'z_precond':
        pass
    elif PRECOND == 'z2_precond':
        idx = np.logical_not(kl == 1.0e-8)
        kl[idx] = kl[idx]**2
    elif PRECOND == 'z_sqrt_precond':
        idx = np.logical_not(kl == 1.0e-8)
        kl[idx] = np.sqrt(kl[idx])
    elif PRECOND == 'default':
        kl = compute_hessian_kernel(iter_cur,SIMU_TYPES,weights_type)
    else:
        raise ValueError(f"PRECOND = {PRECOND} not supported")
    
    return kl

def run(argv):
    if len(argv) !=1 :
        print("Usage: fwat sum_kernel MODEL ")
        print("example: fwat sum_kernel M00.ls")
        exit(1)
    
    # get input 
    MODEL = argv[0]
    iter_cur = int(MODEL.split('.')[0][1:])

    # get fwat params
    with open(f"{PARAM_FILE}","r") as fio:
        param = yaml.safe_load(fio)
    PRECOND = param['optimize']['PRECOND_TYPE']
    SIMU_TYPES = param['simulation']['types']

    assert(PRECOND in ['default','none','z_precond','z2_precond','z_sqrt_precond']), \
        f"PRECOND_TYPE = {PRECOND} not supported "
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

    # get summed user defined kernels 
    grad_user_names = M.grad_names(base=False)
    ksize = _get_database_size(MODEL)
    nsims = len(SIMU_TYPES)
    nkers = len(grad_user_names)
    grad_user = np.zeros((nsims,nkers,ksize),dtype=float)
    for i,s in enumerate(SIMU_TYPES):
        grad_user[i,...] = get_summed_kernel(MODEL,s)

    # compute norm of the gradient for each simulation type
    gnorm = _get_grad_norm(MODEL,SIMU_TYPES,grad_user)

    # open weights_file
    weights_type, success = _get_normalized_weights(param,MODEL,SIMU_TYPES,iter_cur,gnorm)
    if success and myrank == 0:
        print(f"adaptive weights computed based on misfit of this iteration: {weights_type}")

        # save it to FWAT_OPT_DIR
        np.savetxt(f'{OPT_DIR}/weights.txt', weights_type)
    comm.Barrier()

    # print info about normalized grad
    if myrank == 0:
        for i,s in enumerate(SIMU_TYPES):
            print(f"MODEL {MODEL}: Weighted |g| for simu type {s} = {gnorm[i] * weights_type[i]}")
    comm.Barrier()

    # get weighted summed kernel
    grad_user_weighted = np.zeros((nkers,ksize),dtype=float)
    for i in range(nsims):
        grad_user_weighted[...] += grad_user[i,...] * weights_type[i]
    
    # write weighted summed kernel to file
    grad_user_names = M.grad_names(base=False)
    for i in range(nkers):
        outname = KERNEL_DIR + "/proc%06d"%myrank + '_' + grad_user_names[i] + '.bin'
        f = FortranIO(outname,"w")
        arr = np.asarray(grad_user_weighted[i,...],dtype='f4')
        f.write_record(arr)
        f.close()

    # compute preconditioner
    kl = compute_preconditioner(iter_cur,PRECOND,SIMU_TYPES,weights_type)
    outname = KERNEL_DIR + "/proc%06d"%myrank + "_hess_kernel" + '.bin'
    fio = FortranIO(outname,"w")
    kl = np.asarray(kl,dtype='f4')
    fio.write_record(kl)
    fio.close()


def main():
    import sys 
    if len(sys.argv) !=2 :
        print("need 1 parameter: MODEL ")
        print("example: python sum_kernels.py M00.ls")
        exit(1)

    # run
    run(sys.argv[1:])
    

if __name__ == "__main__":
    main()