
import numpy as np 

from fwat.measure.utils import bandpass,taper_window
from fwat.measure.tele.deconit import time_decon
from scipy.signal import convolve 
from scipy.integrate import trapezoid
from fwat.adjoint.MeasureStats import MeasureStats

def measure_adj_cross_conv(
    obs_v:np.ndarray,syn_v:np.ndarray,obs_h:np.ndarray,syn_h:np.ndarray,
    t0:float,dt:float,
    tstart:float,tend:float,
    taper_ratio:float = 0.05
):

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

    # filter
    vsyn = syn_v * 1.
    hsyn = syn_h * 1.
    vobs = obs_v * 1.
    hobs = obs_h * 1.

    # get window info
    lpt, rpt, taper0 = taper_window(t0, dt, nt*2-1, tstart, tend, p=taper_ratio)
    taper = np.zeros(nt*2-1)
    taper[lpt:rpt] = taper0 
    taper[:] = 1.

    # compute convolution of h/v
    chi1 = convolve(vobs,hsyn,'full') * dt
    chi2 = convolve(vsyn,hobs,'full') * dt
    dchi = (chi1 - chi2) * taper

    # misfit function
    mis = 0.5 * trapezoid(dchi**2,dx=dt)

    # reverse several arrays
    v_rev = vobs[::-1].copy()
    h_rev = hobs[::-1].copy()
    
    # adjoint source
    dchi_adj = dchi * taper 
    adj_z:np.ndarray = -convolve(dchi_adj,h_rev,'valid') * dt
    adj_r:np.ndarray = convolve(dchi_adj,v_rev,'valid') * dt

    # stats 
    stats = MeasureStats(
        adj_type="cross_conv",
        misfit=mis,
        tr_chi=mis,
        am_chi=mis,
        tstart = tstart,
        tend = tend
    )

    # also return the cross-conv for visualization
    cc1 = chi1 * taper 
    cc2 = chi2 * taper

    
    return stats,adj_z,adj_r,cc1,cc2