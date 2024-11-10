import numpy as np 
import sys 
from scipy.io import FortranFile 
from tools import *
MAX_PER = 0.02

def main():
    if len(sys.argv) !=5 :
        print("need 5 parameters: MODEL ,OUT_DIR,step_fac,nprocs")
        print("example: python get_lbfgs_step_fac.py M06 OUT_DIR step_fac nprocs")
        exit(1)
    model_str = sys.argv[1]
    MODEL_DIR = "./optimize/MODEL_" +  model_str
    KERNEL_DIR = "./optimize/SUM_KERNELS_" + model_str
    OUT_DIR = sys.argv[2]
    step_fac = float(sys.argv[3])
    nprocs = int(sys.argv[4])

    # print sth
    print("\n get L-BFGS step fac ...")
    print(f"MODEL_DIR = {MODEL_DIR}")
    print(f"KERNEL_DIR = {KERNEL_DIR}")
    print(f"OUTPUT_DIR = {OUT_DIR}")
    print(f"nprocs = {nprocs}")
    
    # get direc_max
    kername_list = get_direc_name_list()
    nker = len(kername_list)
    direc_max = np.float32(-1.0)
    array_size = np.zeros((nprocs),'i4')
    for myrank in range(nprocs):
        for i in range(nker):
            filename = KERNEL_DIR + '/proc%06d'%myrank + '_' + kername_list[i] + '.bin'
            f = FortranFile(filename)
            direc_vp = f.read_reals('f4')
            direc_max = max(direc_max,np.max(np.abs(direc_vp)))
            f.close()

            # get num_elments
            array_size[myrank] = direc_vp.size
    
    # correct step_fac
    step_fac_in_per = MAX_PER
    if step_fac <= 0 or step_fac * direc_max > step_fac_in_per:
        step_fac = step_fac_in_per / direc_max
    print("step_fac dmax relvar = %g %g %g"%(step_fac,direc_max,step_fac * direc_max))

    # write new model
    for myrank in range(nprocs):
        # get model list
        size = array_size[myrank]
        modname_list = get_model_name_list()
        vec = np.zeros((nker,size),'f4')

        # read model
        for i in range(nker):
            filename = MODEL_DIR + '/proc%06d'%myrank + '_' + modname_list[i] + '.bin'
            f = FortranFile(filename)
            vec[i,...] = f.read_reals('f4')
            f.close()

        # convert model if required
        vec = np.float32(convert_md(vec,False))

        # read search direction
        direc = np.zeros((nker,size),'f4')
        for i in range(nker):
            filename = KERNEL_DIR + '/proc%06d'%myrank + '_' + kername_list[i] + '.bin'
            f = FortranFile(filename)
            direc[i,...] = f.read_reals('f4')
            f.close()

        # new model
        vec = vec * np.exp(step_fac * direc)

        # convert back if required
        vec = np.float32(convert_md(vec,True))

        # write model
        for i in range(nker):
            filename = OUT_DIR + '/proc%06d'%myrank + '_' + modname_list[i] + '.bin'
            f = FortranFile(filename,"w")
            f.write_record(vec[i,...])
            f.close()

main()
