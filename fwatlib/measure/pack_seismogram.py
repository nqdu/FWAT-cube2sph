import numpy as np 
import sys  
import h5py
# from tqdm import tqdm

def main():
    if len(sys.argv) <3:
        print("Usage: pack_seismo.py outname *list")
        exit(1)
    
    # load stations
    outname = sys.argv[1]
    seislist = sys.argv[2:]

    # open hdf5 file
    fio = h5py.File(outname,"w")

    for i in range(len(seislist)):
        filename = seislist[i]
        name = filename.split('/')[-1]
        data = np.loadtxt(filename)

        # create database
        fio.create_dataset(name,dtype='f4',shape=data.shape)
        fio[name][:,:] = data[:,:]
    
    fio.close()

if __name__ == "__main__":
    main()