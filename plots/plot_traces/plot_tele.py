import numpy as np
import matplotlib.pyplot as plt 
import sys 
import os

import matplotlib as mpl
mpl.rcParams['lines.linewidth'] = 0.5
mpl.rcParams['font.size'] = 10
mpl.rcParams['xtick.labelsize']=10
mpl.rcParams['ytick.labelsize']=10
mpl.rcParams['axes.labelsize']=10
mpl.rcParams['legend.fontsize'] = 8
mpl.rcParams['legend.fontsize'] = 8
mpl.rcParams['savefig.bbox'] = 'tight'

def find_amp(names,comp,fio):
    amp = -1.
    for n in names:
        data = fio[f'{n}{comp}'][:]
        amp = max(np.max(abs(data)),amp)
    
    return amp

def compute_ak135_time(evla,evlo,evdp,statxt):
    from obspy.taup import TauPyModel

    nsta = statxt.shape[0]
    t_ref = np.zeros((nsta))

    # create taup model
    model = TauPyModel("ak135")
    
    for i in range(nsta):
        stla = float(statxt[i,2])
        stlo = float(statxt[i,3])
        t_ref[i] = model.get_travel_times_geo(evdp,evla,evlo,stla,stlo,['P'])[0].time
    
    return t_ref


def plot_event(line:str,M1:str,solver:str,outdir:str,band:str,
                twb:float,twe:float):
    import h5py 

    M0="M00"
    info = line.split()
    evtid = info[0]
    path2init = f"{solver}/{M0}/{evtid}/OUTPUT_FILES/"
    path2last = f"{solver}/{M1}/{evtid}/OUTPUT_FILES/"

    # compute ref travel time
    src_rec = f"{solver}/../src_rec/"
    t_inj_txt = np.loadtxt(f"{src_rec}/injection_time",dtype=str)
    statxt = np.loadtxt(f"{src_rec}/STATIONS_{evtid}_globe",dtype=str)
    t_inj = 0.
    for i in range(t_inj_txt.shape[0]):
        if t_inj_txt[i,0] == evtid:
            t_inj = float(t_inj_txt[i,1])
            break 
    tref = compute_ak135_time(float(info[1]),float(info[2]),
                            float(info[3]),statxt)
    tref = tref - t_inj 

    # open h5file
    f0syn = h5py.File(f"{path2init}/seismogram.syn.{band}.h5","r")
    fobs = h5py.File(f"{path2init}/seismogram.obs.{band}.h5","r")
    f1syn = h5py.File(f"{path2last}/seismogram.syn.{band}.h5","r")

    # make sure station with same number
    assert(len(f1syn.keys()) == len(fobs.keys())) 

    # find all station names
    names = []
    for ista in range(statxt.shape[0]):
        names.append(statxt[ista,1] + '.' + statxt[ista,0] + '.BX')
    print(f"Plotting {evtid}, no. of stations = {len(names)} ...")


    # create figures
    fig = plt.figure(1,figsize=(6.7,9.3))
    ax1=fig.add_subplot(2,2,1)
    ax2=fig.add_subplot(2,2,2)
    ax3=fig.add_subplot(2,2,3)
    ax4=fig.add_subplot(2,2,4)

    # get scaling factor
    sz = 1. / find_amp(names,'Z',fobs)
    sr = 1. / find_amp(names,'R',fobs)

    # time vector
    t = np.arange(fobs.attrs['npts']) * fobs.attrs['dt']

    # z component
    for i,myname in enumerate(names):
        name = myname  + 'Z'
        synz0 = f0syn[name][:]
        synz1 = f1syn[name][:]
        obsz = fobs[name][:]
        
        name = myname + 'R'
        synr0 = f0syn[name][:]
        synr1 = f1syn[name][:]
        obsr = fobs[name][:]
        
        # legend
        labelo = 'Observed'
        labels = 'Synthetic'
        if i != len(names) - 1:
            labelo = None
            labels = None

        ax1.plot(t,obsz * sz + i ,color='k',label=labelo)
        ax1.plot(t,synz0 * sz + i,color='r',label=labels)
        ax1.vlines(tref[i]-twb,i-0.2,i+0.2,color='m')
        ax1.vlines(tref[i]+twe,i-0.2,i+0.2,color='m')

        ax2.plot(t,obsz * sz + i ,color='k',label=labelo)
        ax2.plot(t,synz1 * sz + i,color='r',label=labels)
        ax2.vlines(tref[i]-twb,i-0.2,i+0.2,color='m')
        ax2.vlines(tref[i]+twe,i-0.2,i+0.2,color='m')

        ax3.plot(t,obsr * sr + i ,color='k',label=labelo)
        ax3.plot(t,synr0 * sr + i,color='r',label=labels)
        ax3.vlines(tref[i]-twb,i-0.2,i+0.2,color='m')
        ax3.vlines(tref[i]+twe,i-0.2,i+0.2,color='m')

        ax4.plot(t,obsr * sr + i ,color='k',label=labelo)
        ax4.plot(t,synr1 * sr + i,color='r',label=labels)
        ax4.vlines(tref[i]-twb,i-0.2,i+0.2,color='m')
        ax4.vlines(tref[i]+twe,i-0.2,i+0.2,color='m')
    
    ax1.set_ylim(-1,len(names))
    ax2.set_ylim(-1,len(names))
    ax3.set_ylim(-1,len(names))
    ax4.set_ylim(-1,len(names))
    ax1.set_title("BXZ before FWI")
    ax2.set_title("BXZ after FWI")
    ax3.set_title("BXR before FWI")
    ax4.set_title("BXR after FWI")
    ax3.set_xlabel("Time (s)")
    ax4.set_xlabel("Time (s)")
    ax1.legend(loc='upper right')
    fig.savefig(f"{outdir}/inv.{evtid}.jpg",dpi=300)
    fig.clear()

    # close h5file
    f0syn.close()
    f1syn.close()
    fobs.close()


def main():
    if len(sys.argv) != 2:
        print("Usage: ./compare_tele.py model(M03)")
        exit(1)
    
    # set directory
    path = "../"

    # ploting window
    twb=5
    twe=45

    # ploting band
    band="T005_T050"

    #### stop here 
    
    solver=f"{path}/solver/"
    seisdir="seismograms_comp/"
    sourcelist=f"{path}/src_rec/sources.dat.tele"

    # create directory
    os.makedirs(seisdir,exist_ok=True)

    # read model name
    M1 = sys.argv[1]


    # loop each line to plot 
    infile = open(sourcelist,"r")
    lines = infile.readlines()
    infile.close()
    for line in lines:
        # plot
        plot_event(line,M1,solver,seisdir,band,twb,twe)

if __name__ == "__main__":
    main()
