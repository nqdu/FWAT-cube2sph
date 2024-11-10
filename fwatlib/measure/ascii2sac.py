import numpy as np
import sys 
from obspy.io.sac import  SACTrace
from obspy import read_events

def main():
    if len(sys.argv) != 4:
        print("Usage: python ascii2sac.py stationfile sourcefile syndir")
        exit(1)
    
    # get params
    stationfile = sys.argv[1]
    sourcefile = sys.argv[2]
    syndir = sys.argv[3] + "/"

    # read station
    statxt = np.loadtxt(stationfile,dtype=str,ndmin=2)
    nsta = statxt.shape[0]

    # read event location
    f = open(sourcefile,"r")
    lines = f.readlines()
    f.close()
    source_dict = {}
    for l in lines:
        if ':' in l:
            info = l.split(':')
            if 'name' in info[0]:
                continue
            source_dict[info[0]] = float(info[1].split()[0])
    if 'latitude' in source_dict.keys():
        evla = source_dict['latitude']
        evlo = source_dict['longitude']
    else:
        evla = source_dict['latorUTM']
        evlo = source_dict['longorUTM']
    evdp = source_dict['depth']

    # loop around all stations
    for i in range(nsta):
        staname_base = statxt[i,1] + "." + statxt[i,0] + ".BX"
        stla = float(statxt[i,2])
        stlo = float(statxt[i,3])

        # read obs data, only header
        tr = SACTrace(evla=evla,evlo=evlo,evdp=evdp,stla=stla,
                    stlo=stlo,stel=0,lcalda=1,knetwk=statxt[i,1],kstnm=statxt[i,0])
        
        for c in ['N','E','Z']:
            staname = staname_base + c
            filename = f'{syndir}' + staname + ".sem.ascii"
            data = np.loadtxt(filename)
            dt = data[1,0] - data[0,0]
            t0 = data[0,0]

            # create sac trace 
            sac = tr.copy()
            if c == 'N':
                sac.cmpaz = 0
                sac.cmpinc = 90.
            elif c == 'E':
                sac.cmpaz = 90.
                sac.cmpinc = 90.
            else:
                sac.cmpaz = 0
                sac.cmpinc = 0.
            sac.kcmpnm = 'BH' + c
            sac.b = t0
            sac.data = data[:,1]
            sac.delta = dt

            # write otu
            filename = f'{syndir}' + staname + ".sac"
            sac.write(filename)

if __name__ == "__main__":
    main()