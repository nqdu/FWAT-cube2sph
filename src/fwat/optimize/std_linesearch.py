import numpy as np 
from mpi4py import MPI 
import yaml
import h5py

from typing import Final

from fwat.optimize.libgll import get_gll_weights
from fwat.optimize.model import FwatModel
from fwat.optimize.search_direction import compute_inner_dot
from fwat.const import PARAM_FILE,LBFGS_FILE,OPT_DIR,NGLL
from fwat.FortranIO import FortranIO

NGLL3 = NGLL**3

def run(argv):
    if len(argv) !=3 :
        print("Usage: fwat linesearch Model fcost fcost1")
        print("example: fwat linesearch M06 fcost fcost1")
        exit(1)

    # get mpi info
    comm = MPI.COMM_WORLD
    myrank = comm.Get_rank()

    # load parameters
    lbfgsfile:Final = LBFGS_FILE
    with open(lbfgsfile, 'r') as f:
        opt = yaml.safe_load(f)

    model_str = argv[0]
    MODEL_DIR = f"./{OPT_DIR}/MODEL_" +  model_str
    KERNEL_DIR = f"./{OPT_DIR}/SUM_KERNELS_" + model_str
    fcost = float(argv[1])
    fcost1 = float(argv[2])
    if myrank == 0: print("\nChecking Wolfe Condtion ...")

    # load FWAT Model
    M = FwatModel()

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
    f = FortranIO(filename,"r")
    nspec = f.read_record('i4')[0]
    _ = f.read_record('i4')[0]
    f.read_record('i4')
    _ = f.read_record('i4')
    for _ in range(3):
        f.read_record('f4')
        f.read_record('i4') # irregular_element_number
        f.read_record('f4')
        f.read_record('f4')
    for _ in range(9):
        f.read_record('f4')
    jaco = f.read_record('f4').reshape(nspec,NGLL3)
    f.close()

    # read search direction, grad, grad for linesearch
    grad_list_base = M.get_grad_names()
    direc_list = M.get_direc_names()
    mod_list = M.get_model_names()
    nker_base = len(grad_list_base)
    nker = len(direc_list)
    direc = np.zeros((nker,nspec,NGLL3),'f4')
    grad_sub_base = np.zeros((nker_base,nspec,NGLL3),'f4')
    grad_sub1_base = np.zeros((nker_base,nspec,NGLL3),'f4')
    mod_sub_base = np.zeros((nker_base,nspec,NGLL3),'f4')
    mod_sub1_base = np.zeros((nker_base,nspec,NGLL3),'f4')

    # read gradient and model
    for i in range(nker_base):
        # read forward kernel 
        filename = f"{KERNEL_DIR}/{grad_list_base[i]}.h5"
        f = h5py.File(filename,"r")
        grad_sub_base[i,:] = f[str(myrank)][:].reshape(nspec,NGLL3)
        f.close()

        # read model
        filename = f"{MODEL_DIR}/proc%06d_{mod_list[i]}.bin" %(myrank)
        fio = FortranIO(filename,"r")
        mod_sub_base[i,:] = fio.read_record('f4').reshape(nspec,NGLL3)
        fio.close()

        # kernel for line search
        filename = f"{KERNEL_DIR}.ls/{grad_list_base[i]}.h5"
        f = h5py.File(filename,"r")
        grad_sub1_base[i,:] = f[str(myrank)][:].reshape(nspec,NGLL3)
        f.close()

        # read model for line search
        filename = f"{MODEL_DIR}.ls/proc%06d_{mod_list[i]}.bin" %(myrank)
        fio = FortranIO(filename,"r")
        mod_sub1_base[i,:] = fio.read_record('f4').reshape(nspec,NGLL3)
        fio.close()
    
    # convert to user defined model
    _,grad_sub = M.convert_kl(mod_sub_base,grad_sub_base)
    _,grad_sub1 = M.convert_kl(mod_sub1_base,grad_sub1_base)

    # read search direction
    for i in range(nker):
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
    paramfile = f"{PARAM_FILE}"
    with open(f"{paramfile}","r") as f:
        params = yaml.safe_load(f)['optimize']
    step_fac_in_per = params['MAX_PER']
    step_fac = opt['alpha']
    if opt['iter'] == opt['iter_start'] and opt['iter_ls'] == 0:
        step_fac = -1
    d0 = np.max(abs(direc))
    direc_max = comm.allreduce(d0,MPI.MAX)

    if step_fac <= 0 or step_fac * direc_max > step_fac_in_per:
        step_fac = step_fac_in_per / direc_max
    opt['alpha'] = float(step_fac)

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

def main():
    import sys 
    if len(sys.argv) !=4 :
        print("need 3 parameters: Model fcost fcost1")
        print("example: python std_linesearch M06 fcost fcost1")
        exit(1)

    run(sys.argv[1:])


if __name__ == "__main__":
    main()
        