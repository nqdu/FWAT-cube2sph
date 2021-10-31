import numpy as np
import obspy 
from glob import glob 
import sys 

def main():
    argv =  sys.argv 
    if len(argv) !=3:
        print("please run this like:")
        print("./this datadir noiselevel")
        exit(1)
    
    # unpack parameters
    datadir = argv[1]
    noiselevel = float(argv[2])

    sourcedir = glob(f"{datadir}/*")
    for s in sourcedir:
        sacfiles = glob(f"{s}/*")
        for file in  sacfiles:
            st = obspy.read(f"{file}",format='sac')
            n = st[0].stats.npts 
            st[0].data = st[0].data * (1.0 + noiselevel * np.random.randn(n))
            st.write(f"{file}",format='sac')

            print(file)


main()