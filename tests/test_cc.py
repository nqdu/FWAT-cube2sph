import numpy as np 
import matplotlib.pyplot as plt 
from fwat.measure.utils import bandpass,taper_window
from fwat.adjoint.cc_misfit import measure_adj_cc 

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

    _,_,_,adj = measure_adj_cc(
        obs,syn,t[0],dt,len(t),Tmin,Tmax,
        tstart,tend)

    lpt,rpt,taper0 = taper_window(t[0],dt,len(t),tstart,tend,p=0.1)
    taper = adj * 0 
    taper[lpt:rpt] = taper0 

    # adjoint source by fd
    syn1 = syn * 1.
    adj_fd = syn1 * 0 

    for it in range(lpt,rpt):
        syn0 = syn[it] * 1.
        dv = 1.0
        syn1[it] = syn0 + dv 
        tr1,_,_,_ = measure_adj_cc(
                    obs,syn1,t[0],dt,len(t),
                    Tmin,Tmax,
                    tstart,tend)

        syn1[it] = syn0 - dv 
        tr2,_,_,_ = measure_adj_cc(
                    obs,syn1,t[0],dt,len(t),
                    Tmin,Tmax,
                    tstart,tend)
        adj_fd[it] = (tr1 - tr2) / (2 * dv)
        syn1[it] = syn0 

    adj_fd = bandpass(adj_fd,dt,1./Tmax,1./Tmin) * taper
    print(adj_fd)

    # plot 
    plt.figure(1)
    plt.plot(t,adj,label='adj')
    plt.plot(t,adj_fd / dt ,label='adj-fd')
    # print(np.max(abs(adj)) / np.max(abs(adj_fd)))
    # plt.plot(t,adj/np.max(abs(adj)))
    # plt.plot(t,adj_fd/np.max(abs(adj_fd)))
    plt.legend()
    plt.savefig("cc_phase.jpg",dpi=300)

if __name__ == "__main__":
    main()