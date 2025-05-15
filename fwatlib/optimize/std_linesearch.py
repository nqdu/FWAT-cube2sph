import numpy as np 
import sys 
from scipy.io import FortranFile 
from mpi4py import MPI 
import os 
from tools import *
from get_lbfgs_direc import compute_inner_dot
from libgll import get_gll_weights
import yaml
import h5py

NGLL = 5; NGLL3 = NGLL**3

def main():
    if len(sys.argv) !=5 :
        print("need 7 parameters: Model lbfgsfile fcost fcost1")
        print("example: python std_linesearch M06 lbfgsfile fcost fcost1")
        exit(1)

    # get mpi info
    comm = MPI.COMM_WORLD
    myrank = comm.Get_rank()

    # load parameters
    lbfgsfile = sys.argv[2]
    with open(lbfgsfile, 'r') as f:
        opt = yaml.safe_load(f)

    model_str = sys.argv[1]
    MODEL_DIR = "./optimize/MODEL_" +  model_str
    KERNEL_DIR = "./optimize/SUM_KERNELS_" + model_str
    fcost = float(sys.argv[3])
    fcost1 = float(sys.argv[4])
    if myrank == 0: print("\nChecking Wolfe Condtion ...")

    # get weights
    w = get_gll_weights()
    weights = np.zeros((NGLL3),dtype='f4')
    for k in range(NGLL):
        for j in range(NGLL):
            for i in range(NGLL):
                weights[k*NGLL*NGLL+j*NGLL+i] = w[i] * w[j] * w[k]

    # read jaco 
    # read nspec/nglob/ibool
    filename = MODEL_DIR + '/proc%06d'%(myrank) + '_external_mesh.bin'
    f = FortranFile(filename)
    nspec = f.read_ints('i4')[0]
    _ = f.read_ints('i4')[0]
    f.read_ints('i4')
    _ = f.read_ints('i4')
    for _ in range(3):
        f.read_reals('f4')
        f.read_ints('i4') # irregular_element_number
        f.read_reals('f4')
        f.read_reals('f4')
    for _ in range(9):
        f.read_reals('f4')
    jaco = f.read_reals('f4').reshape(nspec,NGLL3)
    f.close()

    # read search direction, grad, grad for linesearch
    grad_list = get_gradname_list()
    direc_list = get_direc_name_list()
    nker = len(grad_list)
    direc = np.zeros((nker,nspec,NGLL3),'f4')
    grad_sub = np.zeros((nker,nspec,NGLL3),'f4')
    grad_sub1 = np.zeros((nker,nspec,NGLL3),'f4')
    for i in range(nker):
        # read forward kernel 
        filename = f"{KERNEL_DIR}/{grad_list[i]}.h5"
        f = h5py.File(filename,"r")
        grad_sub[i,:] = f[str(myrank)][:].reshape(nspec,NGLL3)
        f.close()

        # kernel for line search
        filename = f"{KERNEL_DIR}.ls/{grad_list[i]}.h5"
        f = h5py.File(filename,"r")
        grad_sub1[i,:] = f[str(myrank)][:].reshape(nspec,NGLL3)
        f.close()

        # search direction
        filename = f"{KERNEL_DIR}/{direc_list[i]}.h5"
        f = h5py.File(filename,"r")
        direc[i,:] = f[str(myrank)][:].reshape(nspec,NGLL3)
        f.close()

    # initialize alpha_L/R if required
    if opt['iter_ls'] == 0:
        opt['alpha_L'] = 0. 
        opt['alpha_R'] = 0.

    # check if wolfe condition is satisfied
    # wolfe condition
    # f(x + dx) <= f(x) + m1 * alpha * p^T \nablaf
    # p^T grad(x + dx) >= m2 * p^T grad(x)
    m1 = opt['M1']
    m2 = opt['M2']
    q = 0. 
    q1 = 0.
    for i in range(nker):
        q += compute_inner_dot(grad_sub[i,:,:],direc[i,:,:],weights,jaco)
        q1 += compute_inner_dot(grad_sub1[i,:,:],direc[i,:,:],weights,jaco)

    # determine current step_fac
    with open("fwat_params/FWAT.PAR.yaml","r") as f:
        params = yaml.safe_load(f)['optimize']
    step_fac = opt['alpha']

    # check wolfe conditions
    cond1:bool = fcost1 <= (fcost + m1 * step_fac * q)
    cond2:bool = (q1 >= m2 * q)

    if myrank == 0:
        print(f"Wolfe condition : {cond1} {cond2}")

    if cond1 and cond2: # wolfe condition is satisfied 
        opt['flag'] = 'GRAD'
        opt['iter'] += 1
        opt['iter_ls'] = 0

    elif not cond1:
        opt['alpha_R'] = opt['alpha']
        new_alpha = (opt['alpha_L'] + opt['alpha_R'])/2.
        opt['alpha'] = new_alpha
        opt['iter_ls']  += 1

    elif cond1 and (not cond2):
        opt['alpha_L'] = opt['alpha']
        if opt['alpha_R'] != 0. :
            new_alpha = (opt['alpha_L']+opt['alpha_R'])/2.
        else:
            new_alpha = opt['FACTOR'] * opt['alpha']
            
        opt['alpha'] = new_alpha
        opt['iter_ls']  += 1
    
    # write over parameters
    if myrank == 0:
        with open(lbfgsfile, 'w') as f:
            yaml.safe_dump(opt, f)

    # wait until finish
    comm.barrier()


if __name__ == "__main__":
    main()
        