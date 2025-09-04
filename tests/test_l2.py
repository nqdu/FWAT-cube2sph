import numpy as np 
import matplotlib.pyplot as plt 
from fwat.measure.utils import bandpass,taper_window
from fwat.adjoint.l2_misfit import measure_adj_l2
from fwat.measure.measure import measure_adj

def main():
    nt = 1024
    t = np.linspace(0,10,nt)
    dt = t[1]-t[0]
    syn = np.exp(-((t-5)/0.5)**2)
    obs = np.exp(-((t-6.5)/0.5)**2) + 0.2 * np.exp(-((t-4)/0.3)**2) - 0.2 * np.exp(-((t-6)/0.3)**2)

    Tmin = 1.
    Tmax = 50.
    tstart = 3.
    tend = 8.

    syn = bandpass(syn,dt,1./Tmax,1./Tmin)
    obs = bandpass(obs,dt,1./Tmax,1./Tmin)

    tr,_,_,adj = measure_adj_l2(
        obs,syn,t[0],dt,len(t),tstart,tend) 
    # adj = bandpass(adj,dt,1./Tmax,1./Tmin)

    lpt,rpt,taper0 = taper_window(t[0],dt,len(t),tstart,tend,p=0.1)
    taper = adj * 0 
    taper[lpt:rpt] = taper0 

    # adjoint source by measure_adj
    tr_mt,_,_,adj_mt = measure_adj(0,dt,nt,0.,dt,nt,tstart,tend,2,
                         Tmax,Tmin,False,obs,syn)

    print("misfit from user l2: ",tr)
    print("misfit from measure adj: ",tr_mt)

    # plot 
    plt.figure(1)
    plt.plot(t,adj,label='adj')
    plt.plot(t,adj_mt ,label='adj-mt')
    plt.legend()
    plt.savefig("l2.jpg",dpi=300)

if __name__ == "__main__":
    main()