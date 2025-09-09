import numpy as np 
import matplotlib.pyplot as plt 
from fwat.measure.utils import bandpass,taper_window
from fwat.adjoint.cc_misfit import measure_adj_cc 
from fwat.measure.measure import measure_adj

def main():
    nt = 1024
    t = np.linspace(0,10,nt)
    dt = t[1]-t[0]
    syn = np.exp(-((t-4.5)/0.5)**2)
    obs = np.exp(-((t-5.)/0.5)**2) + 0.2 * np.exp(-((t-4)/0.5)**2) - 0.2 * np.exp(-((t-6)/0.5)**2)

    Tmin = 1.
    Tmax = 50.
    tstart = 3.
    tend = 8.

    tr,_,_,adj = measure_adj_cc(
        obs,syn,t[0],dt,len(t),Tmin,Tmax,
        tstart,tend)

    lpt,rpt,taper0 = taper_window(t[0],dt,len(t),tstart,tend,p=0.05)
    taper = adj * 0 
    taper[lpt:rpt] = taper0 

    # adjoint source by measure_adj
    tr_mt,_,_,adj_mt = measure_adj(0,dt,nt,0.,dt,nt,tstart,tend,5,
                         Tmax,Tmin,False,obs,syn)

    print("misfit from user cc: ",tr)
    print("misfit from measure adj: ",tr_mt)

    # plot 
    plt.figure(1)
    plt.plot(t,adj,label='adj')
    plt.plot(t,adj_mt ,label='adj-mt')
    plt.legend()
    plt.savefig("cc_phase.jpg",dpi=300)

if __name__ == "__main__":
    main()