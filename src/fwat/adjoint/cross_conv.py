
import numpy as np 

from fwat.measure.utils import bandpass,taper_window
from fwat.measure.tele.deconit import time_decon
from scipy.signal import convolve 
from scipy.integrate import trapezoid

def measure_adj_cross_conv(
    obs_v,syn_v,obs_h,syn_h,
    t0,dt,min_period,max_period,
    tstart,tend,maxit:int=120,taper_ratio = 0.05):

    """
    Parameters
    ------------
    obs_v/h: np.ndarray
        observed data, R/Z component, shape(nt)
    syn_v/h: np.ndarray
        synthetic data, R/Z component, shape(nt)
    t0: float
        starttime of obs/syn
    dt: float
        time sampling interval
    min/max_period: float
        minimum/maximum period used
    tstart,tend: float
        starttime/endtime of measurement window
    maxit: int
        no. of iter-decon iterations 
    taper_ratio: float
        taper of the window, default = 0.1
    cc_vohs: np.ndarray
        cross-correlation between observed vertical and synthetic horizontal data
    cc_hovs: np.ndarray
        cross-correlation between synthetic vertical and observed horizontal data

    Returns
    ----------------
    tr,am : float
        misfit function (am = tr)
    win: np.ndarray, shape(20)
        measure_adj window
    adj_r,adj_z: np.ndarray
        adjoint source on R/Z component shape(nt)
    """

    # make sure len(obs) == len(syn)
    nt = len(obs_v)
    assert len(obs_v) == len(syn_v), "Observed and synthetic data must have the same length"
    assert len(obs_h) == len(obs_v), "Observed_h and observed_r must have the same length"

    # get window info
    lpt, rpt, taper0 = taper_window(t0, dt, nt, tstart, tend, p=taper_ratio)

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

    # reverse several arrays
    v_rev = vsyn[::-1]
    h_rev = hsyn[::-1]
    chi1_rev = chi1[::-1]
    chi2_rev = chi2[::-1]
    
    # adjoint source
    tmp_z = -convolve(dchi,chi2_rev,'same') * dt 
    tmp_r = convolve(dchi,chi1_rev,'same') * dt 
    adj_r_tp = time_decon(tmp_r,h_rev,dt)
    adj_z_tp = time_decon(tmp_z,v_rev,dt)

    # copy to global arrays
    adj_r = obs_h * 0 
    adj_z = obs_h * 0
    adj_r[lpt:rpt] = adj_r_tp
    adj_z[lpt:rpt] = adj_z_tp 

    # filter
    taper = adj_r * 0 
    taper[lpt:rpt] = taper0 
    adj_r = bandpass(adj_r,dt,1./max_period,1./min_period) * taper
    adj_z = bandpass(adj_z,dt,1./max_period,1./min_period) * taper

    # measure_adj arrays
    tr_chi = mis
    am_chi = mis 
    win_chi = np.zeros((20))
    win_chi[13-1] = 0.5 * np.sum( obs_v**2 + obs_h**2 )
    win_chi[14-1] = 0.5 * np.sum( syn_v**2 + syn_h**2 )
    win_chi[15-1] = tr_chi 
    win_chi[20-1] = nt * dt

    # also return the cross-conv for visualization
    cc1 = obs_v * 0.
    cc1[lpt:rpt] = chi1 * taper0 
    cc2 = obs_v * 0. 
    cc2[lpt:rpt] = chi2 * taper0

    return tr_chi,am_chi,win_chi,adj_r,adj_z,cc1,cc2
