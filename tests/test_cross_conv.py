import numpy as np
import matplotlib.pyplot as plt 

from scipy.signal import convolve

from fwat.measure.measure import measure_adj_cross_conv
from fwat.measure.utils import taper_window,bandpass

def shift_data(a,dt,t0):
    a1 = np.fft.rfft(a)
    om = 2 * np.pi * np.fft.rfftfreq(len(a),dt)
    a1 = a1 * np.exp(-1j * om * t0)

    return np.fft.irfft(a1).real

def measure_adj_cross_conv_mis(
    obs_v,syn_v,obs_h,syn_h,
    t0,dt,min_period,max_period,
    tstart,tend,maxit:int):

    from fwat.measure.utils import bandpass,taper_window
    from fwat.measure.tele.deconit import time_decon
    from scipy.signal import convolve 
    from scipy.integrate import trapezoid

    # make sure len(obs) == len(syn)
    nt = len(obs_v)
    assert len(obs_v) == len(syn_v), "Observed and synthetic data must have the same length"
    assert len(obs_h) == len(obs_v), "Observed_h and observed_r must have the same length"

    # get window info
    lpt, rpt, taper0 = taper_window(t0, dt, nt, tstart, tend, p=0.1)

    # get windowed data
    vsyn = syn_v[lpt:rpt] * taper0 
    hsyn = syn_h[lpt:rpt] * taper0 
    vobs = obs_v[lpt:rpt] * taper0
    hobs = obs_h[lpt:rpt] * taper0 

    # compute convolution of h/v
    chi1 = convolve(vobs,hsyn,'same') * dt
    chi2 = convolve(vsyn,hobs,'same') * dt
    dchi = chi1 - chi2 

    # misfit function
    mis = 0.5 * trapezoid(dchi**2,dx=dt)

    return mis

def main():
    # source time function  
    nt = 1024
    t = np.linspace(0,10,nt)
    dt = t[1]-t[0]
    src = np.exp(-((t-5)/0.5)**2)
    src_t = np.exp(-((t-5)/0.5)**2) + 0.2 * np.exp(-((t-4)/0.3)**2) - 0.2 * np.exp(-((t-6)/0.3)**2)

    # h/v components
    t0 = 4.5
    h = (t > t0) * (t-t0) * np.exp(-(t -t0) / 0.2 )
    v = src * 0.5

    # observed data
    hobs = convolve(h,src_t,'same') * dt
    vobs = convolve(v,src_t,'same') * dt

    # syn data 
    hsyn = convolve(shift_data(h,dt,-0.5),src,'same') * dt
    vsyn = convolve(shift_data(v,dt,-0.5),src,'same') * dt

    # compute adjoint source 
    Tmin = 5
    Tmax = 50.
    tstart = 3.
    tend = 8.
    chi,_,_,adj_r,adj_z = measure_adj_cross_conv(vobs,vsyn,hobs,hsyn,0,dt,Tmin,Tmax,tstart,tend,120)
        # obs_v,syn_v,obs_h,syn_h,
        # t0,dt,min_period,max_period,
        # tstart,tend,maxit:int

    taper = np.zeros((nt))
    lpt,rpt,taper0 = taper_window(0,dt,len(t),tstart,tend,0.1)
    taper[lpt:rpt] = taper0 * 1.

    # fd approx
    adj_z_fd = adj_r * 0 
    adj_r_fd = adj_r * 0 
    zsyn1 = vsyn * 1.
    hsyn1 = hsyn * 1. 

    for it in range(nt):
        # if t[it] < 3  or t[it] > 8:
        #     continue 
        dv = 1.0e-8
        syn0 = vsyn[it] * 1. 
        zsyn1[it] = syn0 + dv
        chi1 = measure_adj_cross_conv_mis(vobs,zsyn1,hobs,hsyn1,0,dt,Tmin,Tmax,tstart,tend,120)

        zsyn1[it] = syn0 - dv
        chi2 = measure_adj_cross_conv_mis(vobs,zsyn1,hobs,hsyn1,0,dt,Tmin,Tmax,tstart,tend,120)

        adj_z_fd[it] = (chi1 - chi2) / (2. * dv)
        zsyn1[it] = syn0 * 1.

        syn0 = hsyn[it] * 1. 
        hsyn1[it] = syn0 + dv
        chi1 = measure_adj_cross_conv_mis(vobs,zsyn1,hobs,hsyn1,0,dt,Tmin,Tmax,tstart,tend,120)

        hsyn1[it] = syn0 - dv
        chi2 = measure_adj_cross_conv_mis(vobs,zsyn1,hobs,hsyn1,0,dt,Tmin,Tmax,tstart,tend,120)

        adj_r_fd[it] = (chi1 - chi2) / (2. * dv)
        hsyn1[it] = syn0 * 1.


    # taper 
    freqmin = 1. / Tmax 
    freqmax = 1. / Tmin
    adj_z_fd = bandpass(adj_z_fd,dt,freqmin,freqmax)
    adj_r_fd = bandpass(adj_r_fd,dt,freqmin,freqmax)
    adj_z_fd *= taper
    adj_r_fd *= taper

    plt.figure(1,figsize=(14,10))
    plt.subplot(211)
    plt.plot(t,adj_z)
    plt.plot(t,adj_z_fd / dt,label='fd')
    plt.legend()

    plt.subplot(212)
    plt.plot(t,adj_r)
    plt.plot(t,adj_r_fd / dt,label='fd')
    plt.legend()

    plt.savefig("cross-conv.jpg",dpi=300)

if __name__ == "__main__":
    main()