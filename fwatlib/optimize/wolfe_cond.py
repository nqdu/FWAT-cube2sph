import numpy as np 
import sys 
from scipy.io import FortranFile 
from tools import * 
from libgll import get_gll_weights

NGLL = 5
NGLL3 = NGLL**3

def main():
    if len(sys.argv) !=7 :
        print("need 7 parameters: MODEL simu_type chi chi1 step_fac nprocs")
        print("example: python get_lbfgs_step_fac.py M06 tele chi chi_next step_fac nprocs")
        exit(1)
    model_str = sys.argv[1]
    MODEL_DIR = "./optimize/MODEL_" +  model_str
    KERNEL_DIR = "./optimize/SUM_KERNELS_" + model_str
    simu_type = sys.argv[2]
    chi = float(sys.argv[3])
    chi1 = float(sys.argv[4])
    step_fac = float(sys.argv[5])
    nprocs = int(sys.argv[6])
    print("\nChecking Wolfe Condtion ...")

    # read event set
    #sta_list = np.loadtxt(f"src_rec/sources.dat.ls.{simu_type}")[:,0]

    # grad/direc list
    direc_list = get_direc_name_list()
    #grad_list = ['alpha_kernel_smooth','beta_kernel_smooth','rhop_kernel_smooth']
    grad_list = get_gradname_list()
    
    # get weights
    w = get_gll_weights()
    weights = np.zeros((NGLL3),dtype='f4')
    for k in range(NGLL):
        for j in range(NGLL):
            for i in range(NGLL):
                weights[k*NGLL*NGLL+j*NGLL+i] = w[i] * w[j] * w[k]

    # sum kernels for each param/rank
    g_list = np.zeros((nprocs))
    for irank in range(nprocs):
        # read nspec/nglob/ibool
        filename = MODEL_DIR + '/proc%06d'%(irank) + '_external_mesh.bin'
        f = FortranFile(filename)
        nspec = f.read_ints('i4')[0]
        nglob = f.read_ints('i4')[0]
        f.read_ints('i4')
        ibool = f.read_ints('i4') - 1
        for _ in range(3):
           f.read_reals('f4')
        f.read_ints('i4') # irregular_element_number
        f.read_reals('f4')
        f.read_reals('f4')
        for _ in range(9):
            f.read_reals('f4')
        jaco = f.read_reals('f4').reshape(nspec,NGLL3)
        f.close()

        # allocate space
        direc_vec = np.zeros((3,nspec,NGLL3),dtype='f4')
        grad_vec = np.zeros((3,nspec,NGLL3),dtype='f4')

        # read gradient/direction and compute inner product
        for iker in range(3):
            filename = KERNEL_DIR + "/proc%06d"%irank + '_' + grad_list[iker] + ".bin"
            f = FortranFile(filename)
            grad_vec[iker,:,:] = f.read_reals('f4').reshape(nspec,NGLL3)
            f.close()

            filename = KERNEL_DIR + "/proc%06d"%irank + '_' + direc_list[iker] + ".bin"
            f = FortranFile(filename)
            #print(grad.max())
            direc_vec[iker,:,:] = f.read_reals('f4').reshape(nspec,NGLL3)
            f.close()
        d1 = np.max(np.abs(direc_vec))
        g1 = np.max(np.abs(grad_vec))
        if d1 == 0. : d1 = 1.
        if g1 == 0. : g1 = 1.
        grad_vec /= g1; direc_vec /= d1 
        dd = np.sum(direc_vec * grad_vec * jaco * weights) * (g1 * d1)

        # compute grad.dot(direc)
        g_list[irank] = dd
    
    # check wolf condition
    g1 = np.sum(g_list)
    interp_sep = False
    backtrack = False
    Armijo_TOL = 0.01
    if chi1 <= chi + Armijo_TOL * step_fac * g1 :
        interp_sep = False 
    elif chi1 < chi:
        interp_sep = True
    else:
        backtrack = True 
    if backtrack:
        step_fac *= 0.5
        print("line search failed! please try step_fac = %g" %(step_fac))
        exit(1)
    
    # interp new step_fac if required
    if interp_sep:
        a = (chi1 - chi - g1 * step_fac) / step_fac**2
        b = g1; c = chi
        alpha = -b / (2. * a)
        chimin = a * alpha**2 + b * alpha + c 
        if alpha > step_fac or alpha < 0 or chimin < 0:
            alpha = step_fac 
            chimin = chi1 
        
        if chimin < 0.1 * chi1: chmin = chi1
    else:
        alpha = step_fac 
        chimin = chi1 

    # get optimal step_fac 
    step_fac = alpha 

    print('step_fac chi_min = %g %g'%(step_fac,chimin))


main()
