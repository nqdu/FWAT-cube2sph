import numpy as np 
import h5py
from fwat.FortranIO import FortranIO

def run(argv):
    if len(argv) !=4 :
        print("Usage: fwat bin2h5 WORK_DIR name nprocs combine=1,decompose=0")
        print("example:  fwat bin2h5 solver/M00/P21/GRADIENT beta_kernel 160 1")
        exit(1)

    WORK_DIR = argv[0] + "/"
    name = argv[1]
    nprocs = int(argv[2])
    combordecomp = int(argv[3])
    assert(combordecomp in [0,1])

    if combordecomp == 1:
        # create hdf5 file
        fio = h5py.File(f"{WORK_DIR}" + f"{name}.h5","w")
        
        for i in range(nprocs):
            procname = str(i)
            dset_name = "proc%06d"%(i) + "_" + name
            filename = WORK_DIR +  dset_name + ".bin"
            fio_f = FortranIO(filename,"r")
            a = fio_f.read_record('f4')
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
            fio_f = FortranIO(filename,"w")
            fio_f.write_record(a)
            fio_f.close()

        fio.close()


def main():
    import sys
    if len(sys.argv) !=5 :
        print("need 4 parameters: WORK_DIR name nprocs combine=1,decompose=0")
        print("example: python bin2h5.py solver/M00/P21/GRADIENT beta_kernel 160 1")
        exit(1)

    run(sys.argv[1:])

    

if __name__ == "__main__":
    main()