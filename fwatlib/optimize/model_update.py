import numpy as np 
import sys 
from scipy.io import FortranFile 
from tools import FwatModel
import yaml
import h5py

def main():
    if len(sys.argv) !=6 :
        print("need 5 parameters: MODEL OUT_DIR paramfile lbfgs nprocs")
        print("example: python model_update.py M06 OUT_DIR paramfile lbfgsnprocs")
        exit(1)
    model_str = sys.argv[1]
    MODEL_DIR = "./optimize/MODEL_" +  model_str + "/"
    KERNEL_DIR = "./optimize/SUM_KERNELS_" + model_str + "/"
    OUT_DIR = sys.argv[2]
    paramfile = sys.argv[3]
    lbfgsfile = sys.argv[4]
    nprocs = int(sys.argv[5])

    # read yaml 
    with open(paramfile,"r") as f:
        params = yaml.safe_load(f)['optimize']
    with open(lbfgsfile,"r") as f:
        opt = yaml.safe_load(f)
    M = FwatModel(paramfile)

    # print sth
    print("\n get L-BFGS step fac ...")
    print(f"MODEL_DIR = {MODEL_DIR}")
    print(f"KERNEL_DIR = {KERNEL_DIR}")
    print(f"OUTPUT_DIR = {OUT_DIR}")
    print(f"nprocs = {nprocs}")

    # get name list
    mname_list = M.get_model_names()
    dname_list = M.get_direc_names()
    nker = len(dname_list)
    nmod = len(mname_list)
    mname_user = dname_list.copy()
    for i in range(nker):
        mname_user[i] = mname_user[i][1:]
    
    # get direc_max
    direc_max = np.float32(-1.0)
    array_size = np.zeros((nprocs),'i4')
    for myrank in range(nprocs):
        for i in range(nker):
            filename = KERNEL_DIR + dname_list[i] + ".h5"
            f = h5py.File(filename,"r")
            direc_vp = f[str(myrank)][:]
            direc_max = max(direc_max,np.max(np.abs(direc_vp)))
            f.close()

            # get num_elments
            array_size[myrank] = direc_vp.size
    
    # correct step_fac
    step_fac = opt['alpha']
    step_fac_in_per = params['MAX_PER']

    # first lbfgs step_fac = step_fac_in_per
    if opt['iter'] == opt['iter_start']:
        step_fac = -1

    if step_fac <= 0 or step_fac * direc_max > step_fac_in_per:
        step_fac = step_fac_in_per / direc_max
    print("step_fac dmax relvar = %g %g %g"%(step_fac,direc_max,step_fac * direc_max))

    # write new model
    for myrank in range(nprocs):
        # get model list
        size = array_size[myrank]
        nmod = len(mname_list)
        vec0 = np.zeros((nmod,size),'f4')

        # read model
        for i in range(nmod):
            filename = MODEL_DIR + '/proc%06d'%myrank + '_' + mname_list[i] + '.bin'
            f = FortranFile(filename,"r")
            vec0[i,...] = f.read_reals('f4')
            f.close()

        # convert model if required
        vec = np.float32(M.convert_md(vec0,False))

        # read search direction
        direc = np.zeros((nker,size),'f4')
        for i in range(nker):
            filename = KERNEL_DIR + dname_list[i] + ".h5"
            f = h5py.File(filename,"r")
            direc[i,...] = f[str(myrank)][:]
            f.close()

        # new model
        vec = M.model_update(vec,step_fac * direc)

        # convert back if required
        vec0 = np.float32(M.convert_md(vec,True))

        # write base model
        for i in range(nmod):
            filename = OUT_DIR + '/proc%06d'%myrank + '_' + mname_list[i] + '.bin'
            f = FortranFile(filename,"w")
            f.write_record(vec0[i,...])
            f.close()
    
    # update opt and write over 
    if opt['flag'] != 'LS': 
        opt['flag'] = 'LS'
        opt['iter_ls'] = 0 
        with open(lbfgsfile,"w") as f:
            yaml.safe_dump(opt,f)

if __name__ == "__main__":
    main()
