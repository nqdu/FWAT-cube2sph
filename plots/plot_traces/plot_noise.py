import numpy as np
import matplotlib.pyplot as plt 
import sys 
import os
from multiprocessing import Pool


import matplotlib as mpl
mpl.rcParams['lines.linewidth'] = 0.5
mpl.rcParams['font.size'] = 10
mpl.rcParams['xtick.labelsize']=10
mpl.rcParams['ytick.labelsize']=10
mpl.rcParams['axes.labelsize']=10
mpl.rcParams['legend.fontsize'] = 8
mpl.rcParams['legend.fontsize'] = 8
mpl.rcParams['savefig.bbox'] = 'tight'

def compute_epc_dist(evla,evlo,statxt):
    from obspy.geodetics import locations2degrees

    stla = np.float64(statxt[:,2])
    stlo = np.float64(statxt[:,3])

    dist = locations2degrees(evla,evlo,stla,stlo) * 6371 * np.pi / 180
    
    if isinstance(dist,np.float64):
        dist = np.array([dist])
    return dist


def plot_event(line:str,M1:str,solver:str,outdir:str,paramfile:str):
    import h5py 
    import yaml

    # get path
    info = line.split()
    evtid = info[0]
    path2last = f"{solver}/{M1}/{evtid}/OUTPUT_FILES/"

    # load paramfile
    with open(paramfile,"r") as f:
        pdict = yaml.safe_load(f)['measure']['noise']
    nbands = len(pdict['FILTER_BANDS'])

    # loop every bands
    for ib in range(nbands):
        Tmin,Tmax = pdict['FILTER_BANDS'][ib]
        vmin,vmax = pdict['GROUPVEL_WIN'][ib]
        band="T%03g_T%03g" %(Tmin,Tmax)

        # compute distance
        src_rec = f"{solver}/../src_rec/"
        statxt = np.loadtxt(f"{src_rec}/STATIONS_{evtid}_globe",dtype=str,ndmin=2)
        dist = compute_epc_dist(float(info[1]),float(info[2]),statxt)

        # sort distance
        idx = np.argsort(dist)
        dist = dist[idx]
        statxt[:,:] = statxt[idx,:]

        # open h5file
        fobs = h5py.File(f"{path2last}/seismogram.obs.{band}.h5","r")
        fsyn = h5py.File(f"{path2last}/seismogram.syn.{band}.h5","r")

        # make sure station with same number
        assert(len(fsyn.keys()) == len(fobs.keys()))

        # find all station names
        names = []
        for ista in range(statxt.shape[0]):
            names.append(statxt[ista,1] + '.' + statxt[ista,0] + '.BX')
        if ib == 0:
            print(f"Plotting {evtid}, no. of stations = {len(names)} ...")

        # create figures
        fig = plt.figure(1,figsize=(6,4))
        ax1=fig.add_subplot(1,1,1)

        # time vector
        t0 = fobs.attrs['t0']
        t = t0 + np.arange(fobs.attrs['npts']) * fobs.attrs['dt']

        # loop each station
        for i,myname in enumerate(names):
            name = myname  + 'Z'
            synz = fsyn[name][:]
            obsz = fobs[name][:]

            # scale them, as we only focus on phases
            synz /= np.max(abs(synz))
            obsz /= np.max(abs(obsz))

            # legend
            labelo = 'Observed'
            labels = 'Synthetic'
            if i != len(names) - 1:
                labelo = None
                labels = None
            
            #time window
            tstart = dist[i] / vmax - Tmax * 0.5 
            tend = dist[i] / vmin + Tmax * 0.5 
            tstart = max(tstart,t0)
            tend = min(tend,t[-1])

            d = i
            ax1.plot(t,obsz + d ,color='k',label=labelo)
            ax1.plot(t,synz + d,color='r',ls='--',label=labels)
            ax1.vlines(tstart,d-0.2,d+0.2,color='m')
            ax1.vlines(tend,d-0.2,d+0.2,color='m')
            tname = statxt[i,1] + '.' + statxt[i,0]
            ax1.text(t[-1],i,tname,size=5,color='k')
        
        #ax1.set_ylim(np.min(dist)*0.95,np.max(dist)*1.05)
        ax1.set_ylim(-1,len(names))
        ax1.yaxis.set_major_locator(mpl.ticker.MaxNLocator(integer=True))
        ax1.set_title(f"Model {M1} : {evtid} {band} ")
        ax1.set_xlabel("Time (s)")
        ax1.set_ylabel("Station No.")
        ax1.legend(loc='upper left')
        fig.savefig(f"{outdir}/{evtid}.{band}.jpg",dpi=300)
        fig.clear()

        # close h5file
        fsyn.close()
        fobs.close()



def main():
    if len(sys.argv) != 2:
        print("Usage: ./compare_tele.py model(M03)")
        exit(1)
    
    # set directory
    solver="../../solver/"
    seisdir="seismograms_noise/"
    sourcelist="../../src_rec/sources.dat.noise"
    paramfile = "../../fwat_params/FWAT.PAR.yaml"

    #### stop here

    # read model name
    M1 = sys.argv[1]

    # create directory
    seisdir = seisdir + "/" + M1
    os.makedirs(seisdir,exist_ok=True)

    # loop each line to plot 
    infile = open(sourcelist,"r")
    lines = infile.readlines()
    infile.close()

    # plot figures
    args = []
    for line in lines:
        args.append((line,M1,solver,seisdir,paramfile))
    pool = Pool(4)
    pool.starmap(plot_event,args)
    pool.close()
    pool.join()

if __name__ == "__main__":
    main()
