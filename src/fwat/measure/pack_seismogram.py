import numpy as np 
import sys  
import h5py
import os 

def run(argv):

    if len(argv) <2:
        print("Usage: fwat pack outname(.h5) *list")
        exit(1)

    # load stations
    outname = argv[0]
    seislist = argv[1:]

    # open hdf5 file
    fout = h5py.File(outname,"w")

    for i in range(len(seislist)):
        filename = seislist[i]
        
        if 'all_seismograms.bin' in filename:
            # open dirname and find how many stations in it
            path = os.path.dirname(filename)
            ncomp = np.loadtxt(f"{path}/../DATA/STATIONS_FILTERED",ndmin=2,dtype=str).shape[0] * 3

            # get file size
            fsize = os.path.getsize(filename)
            nt = (fsize // ncomp - 512 - 8) // 16
            if ncomp * (nt * 2 * 8 + 512 + 8) != fsize:
                print("format error !\n")
                exit(1)
            # print(nt,ncomp)
            # print(fsize,ncomp * (nt * 2 * 4 + 512))
            
            with open(filename,"rb") as fio:
                for ir in range(ncomp):
                    # read name 
                    _ = np.fromfile(fio,'i4',count=1)[0]
                    name = np.fromfile(fio,dtype='S512',count=1)[0]
                    _ = np.fromfile(fio,'i4',count=1)[0]
                    name = str(name.decode('utf-8')).rstrip()

                    # read data
                    data = np.zeros((nt,2))
                    temp = np.fromfile(fio,'f4',count=nt*4).reshape((nt,4)) 
                    data[:,:] = temp[:,1:3] 
                    
                    fout.create_dataset(name,dtype='f4',shape=data.shape)
                        
                    fout[name][:,:] = data[:,:]

        elif 'all_seismograms.ascii' in filename:
            # get how many components
            ncomp = 0
            with open(filename,"r") as fio:
                for line in fio:
                    if '.semd' in line:
                        ncomp += 1
            
            # get no. of time step
            nt = 0
            with open(filename,"r") as fio:
                line = fio.readline()
                while line:
                    line = fio.readline()
                    if '.semd' in line:
                        break
                    nt = nt + 1
            
            # read seismograms
            with open(filename,"r") as fio:
                for ir in range(ncomp):
                    line = fio.readline()
                    name = line.split()[0]
                    data = np.zeros((nt,2),dtype='f4')
                    for it in range(nt):
                        line = fio.readline()
                        info = line.split()[:2]
                        data[it,:] = float(info[0]),float(info[1])
                    fout.create_dataset(name,dtype='f4',shape=data.shape)    
                    fout[name][:,:] = data[:,:]
        else:
            name = filename.split('/')[-1]
            data = np.loadtxt(filename)

            # create database
            fout.create_dataset(name,dtype='f4',shape=data.shape)
            fout[name][:,:] = data[:,:]
    
    fout.close()

def main():
    if len(sys.argv) <3:
        print("Usage: pack_seismo.py outname *list")
        exit(1)

    run(sys.argv[1:])    


if __name__ == "__main__":
    main()