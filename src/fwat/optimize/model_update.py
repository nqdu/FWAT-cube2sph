import numpy as np 
import yaml
import h5py
from mpi4py import MPI

from typing import Final

from fwat.optimize.model import FwatModel
from fwat.const import OPT_DIR,PARAM_FILE,LBFGS_FILE
from fwat.FortranIO import FortranIO


def run(argv):
    if len(argv) !=2 :
        print("Usage: fwat update MODEL OUT_DIR")
        print("example: fwat update M06 OUT_DIR")
        exit(1)
    model_str = argv[0]
    MODEL_DIR = f"{OPT_DIR}/MODEL_" +  model_str + "/"
    KERNEL_DIR = f"{OPT_DIR}/SUM_KERNELS_" + model_str + "/"
    OUT_DIR = argv[1]

    # mpi 
    comm = MPI.COMM_WORLD
    myrank = comm.Get_rank()

    # read yaml
    paramfile: Final = PARAM_FILE
    lbfgsfile: Final = LBFGS_FILE
    with open(paramfile,"r") as f:
        params = yaml.safe_load(f)['optimize']
    with open(lbfgsfile,"r") as f:
        opt = yaml.safe_load(f)
    M = FwatModel(paramfile)

    # print sth
    if myrank == 0:
        print("\n update model ...")
        print(f"MODEL_DIR = {MODEL_DIR}")
        print(f"KERNEL_DIR = {KERNEL_DIR}")
        print(f"OUTPUT_DIR = {OUT_DIR}")

    # get name list
    mname_list = M.get_model_names()
    dname_list = M.get_direc_names()
    nker = len(dname_list)
    nmod = len(mname_list)
    mname_user = dname_list.copy()
    for i in range(nker):
        mname_user[i] = mname_user[i][1:]
    
    # get direc_max
    direc_max0 = -1.
    array_size = 0
    for i in range(nker):
        filename = KERNEL_DIR + dname_list[i] + ".h5"

        # print 
        if myrank == 0:
            print(f"reading search direction from {filename}")
        
        # load
        f = h5py.File(filename,"r")
        direc_vp = f[str(myrank)][:]
        direc_max0 = max(direc_max0,np.max(np.abs(direc_vp)))
        f.close()

        # get num_elments
        array_size = direc_vp.size
    direc_max = comm.allreduce(direc_max0,MPI.MAX)
    
    # correct step_fac
    step_fac = opt['alpha']
    step_fac_in_per = params['MAX_PER']

    # first lbfgs step_fac = step_fac_in_per    
    if opt['iter'] == opt['iter_start'] + 1 and opt['iter_ls'] == 0:
        # init step_fac is 1
        step_fac = 1.

    # limit step_fac
    flag = (step_fac > 0) and (step_fac * direc_max > step_fac_in_per) and opt['iter_ls'] ==0 
    if step_fac <= 0 or flag:
        step_fac = step_fac_in_per / direc_max
    
    if myrank == 0:
        print("step_fac dmax relvar = %g %g %g"%(step_fac,direc_max,step_fac * direc_max))
    opt['alpha'] = float(step_fac) 

    # write new model
    # get model list
    size = array_size
    nmod = len(mname_list)
    vec0 = np.zeros((nmod,size),'f4')

    # read model
    for i in range(nmod):
        filename = MODEL_DIR + '/proc%06d'%myrank + '_' + mname_list[i] + '.bin'

        if myrank == 0:
            print(f'reading model from {filename}')

        f = FortranIO(filename,"r")
        vec0[i,...] = f.read_record('f4')
        f.close()

    # convert model if required
    vec = np.asarray(M.convert_model(vec0,False),dtype='f4')

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
    vec0 = np.asarray(M.convert_model(vec,True),dtype='f4')

    # write base model
    for i in range(nmod):
        filename = OUT_DIR + '/proc%06d'%myrank + '_' + mname_list[i] + '.bin'

        if myrank == 0:
            print(f'wrting new model to {filename}')
        
        f = FortranIO(filename,"w")
        f.write_record(vec0[i,...])
        f.close()

    # update opt and write over 
    if opt['flag'] != 'LS': 
        opt['flag'] = 'LS'
        opt['iter_ls'] = 0 

        if myrank == 0:
            with open(lbfgsfile,"w") as f:
                yaml.safe_dump(opt,f)
    
    # sync
    comm.Barrier()


def main():
    import sys 
    if len(sys.argv) !=3 :
        print("need 2 parameters: MODEL OUT_DIR")
        print("example: python model_update.py M06 OUT_DIR")
        exit(1)
    run(sys.argv[1:])

if __name__ == "__main__":
    main()
