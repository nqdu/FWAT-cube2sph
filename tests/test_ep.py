import numpy as np 
import matplotlib.pyplot as plt 
from fwat.measure.utils import bandpass,taper_window
from fwat.measure.measure import measure_adj_exphase

def exphase_misfit(obs,syn,t0,dt,nt,
               min_period,max_period,
               tstart,tend,water=0.1,
               taper_ratio = 0.1):
    """
    Parameters
    ------------
    obs: np.ndarray
        observed data, shape(nt)
    syn: np.ndarray
        synthetic data, shape(nt)
    t0,dt,nt: float,float,int
        starttime/sampling/no.of points of adjoint source
    min/max_period: float
        minimum/maximum period used
    tstart,tend: float
        starttime/endtime of measurement window
    water: float
        waterlevel for synthetic envelope, default = 0.1
    taper_ratio: float
        taper of the window, default = 0.1

    Returns
    ----------------
    tr,am : float
        misfit function (am = tr)
    win: np.ndarray, shape(20)
        measure_adj window
    adj: np.ndarray
        adjoint source, shape(nt)
    """
    from scipy.signal import hilbert
    from fwat.measure.utils import bandpass,taper_window
    from fwat.measure.utils import interpolate_syn
    from scipy.integrate import trapezoid

    # make sure len(obs) == len(syn)
    assert len(obs) == len(syn), "Observed and synthetic data must have the same length"

    # get window info
    lpt, rpt, taper0 = taper_window(t0, dt, nt, tstart, tend, p=taper_ratio)

    # compute hilbert transform
    s = syn[lpt:rpt] * taper0 
    d = obs[lpt:rpt] * taper0 
    syn_a = hilbert(s)
    obs_a = hilbert(d)
    Hsyn = np.imag(syn_a)
    Hobs = np.imag(obs_a)
    Es = np.abs(syn_a)
    Ed = np.abs(obs_a)

    # determine waterlevel threshold
    w_s = water * np.max(Es)
    w_d = water * np.max(Ed)
    Es_wtr = Es + w_s 
    Ed_wtr = Ed + w_d 

    # get real/imag part
    dR = d / Ed_wtr - s / Es_wtr 
    dI = Hobs / Ed_wtr - Hsyn / Es_wtr

    # compute misfit
    misfit = trapezoid(dR**2,dx=dt) + trapezoid(dI**2,dx=dt)
    misfit = misfit * 0.5

    return misfit

def main():
    nt = 1024
    t = np.linspace(0,10,nt)
    dt = t[1]-t[0]
    syn = np.exp(-((t-5)/0.5)**2)
    obs = np.exp(-((t-5.5)/0.5)**2) + 0.2 * np.exp(-((t-4)/0.3)**2) - 0.2 * np.exp(-((t-6)/0.3)**2)

    Tmin = 1.
    Tmax = 50.
    tstart = 3.
    tend = 8.

    _,_,_,adj = measure_adj_exphase(
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
        dv = 1.0e-6
        syn1[it] = syn0 + dv 
        tr1 = exphase_misfit(
                    obs,syn1,t[0],dt,len(t),
                    Tmin,Tmax,
                    tstart,tend)

        syn1[it] = syn0 - dv 
        tr2 = exphase_misfit(
                    obs,syn1,t[0],dt,len(t),
                    Tmin,Tmax,
                    tstart,tend)
        adj_fd[it] = (tr1 - tr2) / (2 * dv)
        syn1[it] = syn0 

    adj_fd = bandpass(adj_fd,dt,1./Tmax,1./Tmin) * taper

    # plot 
    plt.figure(1)
    plt.plot(t,adj,label='adj')
    plt.plot(t,adj_fd / dt,label='adj-fd')
    # print(np.max(abs(adj)) / np.max(abs(adj_fd)))
    # plt.plot(t,adj/np.max(abs(adj)))
    # plt.plot(t,adj_fd/np.max(abs(adj_fd)))
    plt.legend()
    plt.savefig("exp_phase.jpg",dpi=300)

if __name__ == "__main__":
    main()