import sys 
import os 
from scipy.io import FortranFile 
import numpy as np 
from tools import *
import h5py

def main():
    if len(sys.argv) !=5 :
        print("need 4 parameters: WORK_DIR name nprocs combine=1,decompose=0")
        print("example: python bin2h5.py solver/M00/P21/GRADIENT beta_kernel 160 1")
        exit(1)

    WORK_DIR = sys.argv[1] + "/"
    name = sys.argv[2]
    nprocs = int(sys.argv[3])
    combordecomp = int(sys.argv[4])
    assert(combordecomp in [0,1])

    if combordecomp == 1:
        # create hdf5 file
        fio = h5py.File(f"{WORK_DIR}" + f"{name}.h5","w")
        
        for i in range(nprocs):
            procname = str(i)
            dset_name = "proc%06d"%(i) + "_" + name
            filename = WORK_DIR +  dset_name + ".bin"
            fio_f = FortranFile(filename,"r")
            a = fio_f.read_reals('f4')
            fio_f.close()

            # create a h5 database
            fio.create_dataset(procname,dtype='f4',shape = a.shape)
            fio[procname][:] = a[:]
        fio.close()

    else:
        # open hdf5 file
        fio = h5py.File(f"{WORK_DIR}" + f"{name}.h5","r")

        for i in range(nprocs):
            # get dataset
            procname = str(i)
            dset_name = "proc%06d"%(i) + "_" + name
            a = fio[procname][:]

            filename = WORK_DIR +  dset_name + ".bin"
            fio_f = FortranFile(filename,"w")
            fio_f.write_record(a)
            fio_f.close()

        fio.close()


    

if __name__ == "__main__":
    main()