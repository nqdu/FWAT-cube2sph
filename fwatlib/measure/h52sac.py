from obspy.io.sac import  SACTrace
import h5py 
import sys

def convert2sac(workdir,saclist):
    # open h5file 
    fio = h5py.File(workdir + "/seismograms.h5","w")

    # first get header 
    sachd = SACTrace.read(saclist[0],headonly=True)
    dt = sachd.delta 
    t0 = sachd.b 
    npts = sachd.npts 
    fio.attrs['dt'] = dt 
    fio.attrs['t0'] = t0 
    fio.attrs['npts'] = npts 

    for i in range(len(saclist)):
        filename = saclist[i]
        tr = SACTrace.read(filename)

        # get name
        dsetname = tr.knetwk + "." + tr.kstnm + "." + tr.kcmpnm

        # add data to fio 
        fio.create_dataset(dsetname,shape=sachd.data.shape,dtype='f4')
        fio[dsetname][:] = sachd.data 
    
    fio.close()

def main():
    if len(sys.argv) <= 3:
        print("Usage: python sac2h5.py 2sac=0/2h5=1 workdir *list")
        exit(1)

    # get input args
    flag = int(sys.argv[1])
    workdir = sys.argv[2]
    namelist = sys.argv[3:]