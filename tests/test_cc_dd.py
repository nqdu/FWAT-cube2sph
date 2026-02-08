import numpy as np 
import matplotlib.pyplot as plt 
from fwat.measure.measure import measure_adj
from fwat.measure.utils import bandpass,taper_window
from fwat.adjoint.cc_misfit import measure_adj_cc_dd

def main():
    nt = 1024
    t = np.linspace(0,10,nt)
    dt = t[1]-t[0]
    syn = np.exp(-((t-4.5)/0.5)**2)
    obs = np.exp(-((t-5.)/0.5)**2)
    syn2 = np.exp(-((t-4.7)/0.5)**2)
    obs2 = np.exp(-((t-5.6)/0.5)**2)

    Tmin = 1.   
    Tmax = 50.
    tstart = 3.
    tend = 8.

    stats,adj_i,adj_j = measure_adj_cc_dd(
        obs,syn,obs2,syn2,t[0],dt,len(t),Tmin,Tmax,
        tstart,tend,tstart,tend)

    lpt,rpt,taper0 = taper_window(t[0],dt,len(t),tstart,tend,p=0.05)
    taper_i = adj_i * 0 
    taper_i[lpt:rpt] = taper0 
    taper_j = adj_j * 0
    taper_j[lpt:rpt] = taper0

    print("misfit/deltat_dd from user cc: %g %g" %(stats.misfit,stats.tshift))

    # plot 
    plt.figure(1)
    plt.subplot(211)
    plt.plot(t,syn,label='syn')
    plt.plot(t,obs,label='obs') 
    plt.legend()
    plt.subplot(212)
    plt.plot(t,syn2,label='syn2')
    plt.plot(t,obs2,label='obs2')
    plt.legend() 

    plt.figure(2)
    plt.subplot(211)
    plt.plot(t,adj_i,label='adj_i')
    plt.legend()
    plt.subplot(212)
    plt.plot(t,adj_j,label='adj_j')
    plt.legend()
    plt.savefig("cc_phase.jpg",dpi=300)

    plt.show()

if __name__ == "__main__":
    main()