import numpy as np 
import sys 
from scipy.io import FortranFile 
from tools import FwatModel

def main():
    if len(sys.argv) !=5 :
        print("need 4 parameters: MODEL_DIR paramfile nprocs")
        print("example: python write_user_model M06 mdtype kltype nprocs")
        exit(1)

    # get params
    MODEL_DIR = sys.argv[1]
    mdtype = sys.argv[2]
    kltype = int(sys.argv[3])
    nprocs = int(sys.argv[4])

    # init FwatModel if required
    if mdtype == "iso":
        print("it's ISO model, vp,vs,rho are already existed")
        print("do nothing ...")
        exit(0)
    M = FwatModel(None,mdtype,kltype)

    # print sth
    print("\n write user model from cijkl ...")
    print(f"MODEL_DIR = {MODEL_DIR}")
    print(f"nprocs = {nprocs}")

    # get name list
    mname_list = M.get_model_names()
    nmod = len(mname_list)

    # get size for each proc
    array_size = np.zeros((nprocs),'i4')
    for myrank in range(nprocs):
        filename = MODEL_DIR + "/proc%06d_" %(myrank) + mname_list[0] + ".bin"
        fio = FortranFile(filename,"r")
        x = fio.read_reals('f4') 
        array_size[myrank] = x.size
        fio.close()

    # loop each proc 
    for myrank in range(nprocs):
        size = array_size[myrank]
        md = np.zeros((nmod,size))

        # read base model
        for im in range(nmod):
            filename = MODEL_DIR + "/proc%06d_" %(myrank) + mname_list[im] + ".bin"
            fio = FortranFile(filename,"r")
            x = fio.read_reals('f4') 
            fio.close()
            md[im,:] = x * 1.
        
        # convert to user model
        md_new,all_names = M.convert_md_visual(md)
        md_new = np.float32(md_new)

        # write user model
        nker = len(all_names)
        for i in range(nker):
            filename = MODEL_DIR + "/proc%06d_" %(myrank) + all_names[i] + ".bin"
            fio = FortranFile(filename,"w")
            fio.write_record(md_new[i,...])
            fio.close()

if __name__ == "__main__":
    main()