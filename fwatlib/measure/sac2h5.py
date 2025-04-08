from obspy.io.sac import  SACTrace
import h5py 
import sys

def convert2h5(outname,saclist):
    # open h5file 
    if '.h5' not in outname:
        outname = outname + '.h5'
        
    fio = h5py.File(outname,"w")
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
        fio.create_dataset(dsetname,shape=tr.data.shape,dtype='f4')
        fio[dsetname][:] = tr.data 
    
    fio.close()


def main():
    if len(sys.argv) < 3:
        print("Usage: python sac2h5.py outname *list")
        exit(1)

    # get input args
    outname = sys.argv[1]
    namelist = sys.argv[2:]

    # conversion ... 
    convert2h5(outname,namelist)

if __name__ == "__main__":
    main()