import numpy as np 
import sys 
from scipy.io import FortranFile 
from mpi4py import MPI 
import os 
from tools import FwatModel
import h5py
import yaml

NGLL = 5; NGLL3 = NGLL**3

def get_model_grad(iter:int,nspec:int,M:FwatModel):
    comm = MPI.COMM_WORLD
    myrank = comm.Get_rank()

    # initialize
    grad_list = M.get_grad_names(base=True)
    mod_list = M.get_model_names()
    nmod = len(grad_list)
    mod_vec = np.zeros((nmod,nspec,NGLL3),'f4')
    ker_vec = np.zeros((nmod,nspec,NGLL3),'f4')

    # print log info
    KERNEL_DIR = './optimize/SUM_KERNELS_M%02d/'%(iter)
    MODEL_DIR = './optimize/MODEL_M%02d/'%(iter)
    if myrank == 0: print('reading ',KERNEL_DIR,MODEL_DIR)

    # reading model/kernels
    for i in range(nmod):
        # read kernel 
        filename = KERNEL_DIR + grad_list[i] + ".h5"
        fio = h5py.File(filename,"r")
        ker_vec[i,:,:] = fio[str(myrank)][:].reshape(nspec,NGLL3)
        fio.close()

        # read model
        filename = MODEL_DIR + "/proc%06d"%myrank + '_' + mod_list[i] + ".bin"
        f = FortranFile(filename)
        mod_vec[i,:,:] = f.read_reals('f4').reshape(nspec,NGLL3)
        f.close()

    # convert kernel and model to required type
    mod_vec1,ker_vec1 = M.convert_kl(mod_vec,ker_vec)
    ker_vec1 = np.float32(ker_vec1)

    # convert model to required optimzed type
    mod_vec1 = np.float32(M.get_used_model(mod_vec1))

    return mod_vec1,ker_vec1

def compute_inner_dot(a,b,weights,jaco):
    """
    compute inner product for mpi vector
    """
    # compute normalize vector 
    comm = MPI.COMM_WORLD
    amax = np.max(abs(a))
    bmax = np.max(abs(b))
    c1 = comm.allreduce(amax,MPI.MAX)
    c2 = comm.allreduce(bmax,MPI.MAX)

    # set threshold
    if c1 == 0.: c1 = 1.
    if c2 == 0.: c2 = 1.

    # compute inner product in this rank
    q = np.sum((a / c1) * (b / c2) * jaco * weights)
    q_sum = comm.allreduce(q,MPI.SUM)

    return q_sum * c1 * c2 

def get_lbfgs_direc(iter:int,paramfile:str,M:FwatModel):

    from libgll import get_gll_weights

    comm = MPI.COMM_WORLD
    myrank = comm.Get_rank()

    # get weights
    w = get_gll_weights()
    weights = np.zeros((NGLL3),dtype='f4')
    for k in range(NGLL):
        for j in range(NGLL):
            for i in range(NGLL):
                weights[k*NGLL*NGLL+j*NGLL+i] = w[i] * w[j] * w[k]
    
    # read jacobian
    # read nspec/nglob/ibool
    filename = './optimize/MODEL_M%02d'%(iter) + '/proc%06d'%(myrank) + '_external_mesh.bin'
    f = FortranFile(filename)
    nspec = f.read_ints('i4')[0]
    _ = f.read_ints('i4')[0]
    f.read_ints('i4')
    _= f.read_ints('i4')
    for _ in range(3):
        f.read_reals('f4')
    f.read_ints('i4') # irregular_element_number
    f.read_reals('f4')
    f.read_reals('f4')
    for _ in range(9):
        f.read_reals('f4')
    jaco = f.read_reals('f4').reshape(nspec,NGLL3)
    f.close()

    # alloc space for temporary vars 
    a = np.zeros((1000),dtype='f4'); p = a * np.float32(1.) 

    # read current gradient
    _,q_vec = get_model_grad(iter,nspec,M)
    grad_bak = q_vec.copy()

    # load lbfgs parameters
    with open(paramfile,"r") as f:
        pdict = yaml.safe_load(f)

    # current search direction is -grad
    iter_start:int = pdict['iter_start']

    # get istore
    MAX_HIS = pdict['MSTORE']
    iter_store = iter - MAX_HIS
    if iter_store <= iter_start: iter_store = iter_start

    # forward 
    if myrank == 0:
        print("**************************************")
        print("******* starting backward store *****")
        print("**************************************")
    
    for istore in range(iter-1,iter_store-1,-1):
        model1,grad1 = get_model_grad(istore+1,nspec,M)
        model0,grad0 = get_model_grad(istore,nspec,M)
        grad_diff = grad1 - grad0 
        model_diff = model1 - model0

        # compute p and a
        p_sum = compute_inner_dot(grad_diff,model_diff,weights,jaco)
        p[istore] = 1. / p_sum 
        a_sum = compute_inner_dot(model_diff,q_vec,weights,jaco)
        a[istore] = p[istore] * a_sum

        if myrank == 0: print('a,p = ',a[istore],p[istore])
        q_vec -= a[istore] * grad_diff
    
    istore = iter - 1
    model1,grad1 = get_model_grad(istore+1,nspec,M)
    model0,grad0 = get_model_grad(istore,nspec,M)
    grad_diff = grad1 - grad0 
    model_diff = model1 - model0

    # read hess from hdf5
    filename = './optimize/SUM_KERNELS_M%02d'%(iter)
    filename += "/hess_kernel.h5"
    if myrank == 0: print(f"reading hess from {filename}")
    fio = h5py.File(filename,"r")
    hess = fio[str(myrank)][:].reshape(nspec,NGLL3)
    fio.close()

    # get hessian max
    hess_max_loc = np.max(np.abs(hess))
    hess_max = comm.allreduce(hess_max_loc,MPI.MAX)
    if myrank == 0: print("hess maximum :  %g"% hess_max)
    
    #
    p_k_up_sum = compute_inner_dot(grad_diff,model_diff,weights,jaco)
    p_k_down_sum = compute_inner_dot(grad_diff,grad_diff,weights,jaco)
    p_k = p_k_up_sum / p_k_down_sum
    if myrank == 0: print("p_k = ",p_k) 
    r_vec = p_k * (hess * q_vec) 

    if myrank == 0:
        print("**************************************")
        print("******* starting forward store *****")
        print("**************************************")
    
    for istore in range(iter_store,iter):
        model1,grad1 = get_model_grad(istore+1,nspec,M)
        model0,grad0 = get_model_grad(istore,nspec,M)
        grad_diff = grad1 - grad0 
        model_diff = model1 - model0
        b_sum = compute_inner_dot(grad_diff,r_vec,weights,jaco)
        b = p[istore] * b_sum 

        if myrank == 0: print("a,b = ",a[istore],b)
        r_vec += model_diff * (a[istore] - b)
    
    # search direc
    direc = -r_vec

    # check the angle between search direction and negative grad
    grad_bak = - grad_bak
    grad_sum = compute_inner_dot(grad_bak,direc,weights,jaco)
    grad_norm = compute_inner_dot(grad_bak,grad_bak,weights,jaco)
    direc_norm = compute_inner_dot(direc,direc,weights,jaco)
    grad_norm = np.sqrt(grad_norm); direc_norm = np.sqrt(direc_norm)
    a = grad_sum / (direc_norm * grad_norm)
    theta = np.arccos(a) * 180 / np.pi 

    if theta <= 92.:
        if myrank == 0: print("The search direction is accepted!")
    else:
        if myrank == 0: 
            print("The search direction is not accepted!")
            print("clear previous information !")
            print("use negative grad as search direction")
            direc = -grad_bak

            # write new info 
            pdict['iter_start'] = iter
            with open(paramfile,"w") as f:
                yaml.safe_dump(pdict,f)

    comm.barrier()

    return direc

def get_sd_direction(iter:int,paramfile:str,M:FwatModel):
    comm = MPI.COMM_WORLD
    myrank = comm.Get_rank()

    # read nspec/nglob/ibool
    filename = './optimize/MODEL_M%02d'%(iter) + '/proc%06d'%(myrank) + '_external_mesh.bin'
    f = FortranFile(filename)
    nspec = f.read_ints('i4')[0]
    f.close()

    # get gradient
    _,q_vec = get_model_grad(iter,nspec,M)
    q_vec = np.float32(-q_vec)

    return q_vec


def get_search_direction(iter:int,paramfile:str):

    comm = MPI.COMM_WORLD
    myrank = comm.Get_rank()

    # load lbfgs parameters
    with open(paramfile,"r") as f:
        pdict = yaml.safe_load(f)

    # initialize FwatModel
    M = FwatModel()

    # current search direction is -grad
    iter_start:int = pdict['iter_start']
    if iter == iter_start:
        direc = get_sd_direction(iter,paramfile,M)
    else:
        direc = get_lbfgs_direc(iter,paramfile,M)

    # get min/max
    direc_list =  M.get_direc_names()
    if myrank == 0:
        print("")
    for i in range(len(direc_list)):
        min_vp = np.min(direc[i,...])
        max_vp = np.max(direc[i,...])
        min_all = comm.allreduce(min_vp,MPI.MIN)
        max_all = comm.allreduce(max_vp,MPI.MAX)
        if myrank == 0:
            print(f"search direction : {direc_list[i]} min/max = %g %g"%(min_all,max_all))
    
    # save search
    direc = np.float32(direc)
    ndirec = direc.shape[0]
    for i in range(ndirec):
        outname =  './optimize/SUM_KERNELS_M%02d'%(iter) + '/proc%06d'%(myrank) + '_' + direc_list[i] + ".bin"
        f = FortranFile(outname,'w')
        f.write_record(direc[i,:,:])
        f.close()
    
def main():
    if len(sys.argv) !=3 :
        print("need 2 parameters: iter paramfile")
        print("example: mpirun -np 160 python get_lbfgs_step_direc.py iter fwat_params/lbfgs.yaml")
        exit(1)

    # get mpi info
    comm = MPI.COMM_WORLD
    myrank = comm.Get_rank()
    
    # get input 
    iter = int(sys.argv[1])
    paramfile = sys.argv[2]
    KERNEL_DIR = './optimize/SUM_KERNELS_M%02d'%(iter)
    if myrank == 0:  
        os.makedirs(KERNEL_DIR,exist_ok=True)
        print(f"get L-BFGS search direction for iteration {iter}")
    
    # get search direction
    get_search_direction(iter,paramfile)

if __name__ == "__main__":
    main()
        
