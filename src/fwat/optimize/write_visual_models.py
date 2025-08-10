import numpy as np 
from scipy.io import FortranFile 
from fwat.optimize.model import FwatModel
import glob

def run(argv):
    if len(argv) !=3 :
        print("need 3 parameters: MODEL_DIR model_type kernel_type")
        print("model_type and kernel_type can be find in class FwatModel")
        print("example: fwat-model visual M06 dtti 1")
        exit(1)

    # get params
    MODEL_DIR = argv[0]
    mdtype = argv[1]
    kltype = int(argv[2])

    # init FwatModel if required
    if mdtype == "iso":
        print("it's ISO model, vp,vs,rho are already existed")
        print("do nothing ...")
        exit(0)
    M = FwatModel(None,mdtype,kltype)

    # get name list
    mname_list = M.get_model_names()
    nmod = len(mname_list)

    # get how many files in the MODEL_DIR
    nfiles = len(glob.glob(f"{MODEL_DIR}/proc*_{mname_list[0]}.bin"))

    # print sth
    print("write visual model from cijkl ...")
    print(f"MODEL_DIR = {MODEL_DIR}")
    print(f"nfiles = {nfiles}\n")

    # get size for each proc
    array_size = np.zeros((nfiles),'i4')
    for myrank in range(nfiles):
        filename = MODEL_DIR + "/proc%06d_" %(myrank) + mname_list[0] + ".bin"
        fio = FortranFile(filename,"r")
        x = fio.read_reals('f4') 
        array_size[myrank] = x.size
        fio.close()

    # loop each proc 
    for myrank in range(nfiles):
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
            print(f"working on {filename}")
            fio = FortranFile(filename,"w")
            fio.write_record(md_new[i,...])
            fio.close()


def main():
    import sys 
    if len(sys.argv) !=4 :
        print("need 3 parameters: MODEL_DIR model_type kernel_type")
        print("example: python write_user_model.py M06 mdtype kltype")
        exit(1)

    run(sys.argv[1:])
    
if __name__ == "__main__":
    main()