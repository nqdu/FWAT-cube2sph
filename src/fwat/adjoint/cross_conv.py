
import numpy as np 
from scipy.signal import convolve 
from scipy.integrate import trapezoid

from fwat.measure.utils import taper_window
from fwat.adjoint.MeasureStats import MeasureStats

def measure_adj_cross_conv(
    obs_v:np.ndarray,syn_v:np.ndarray,
    obs_h:np.ndarray,syn_h:np.ndarray,
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
    tstart,tend: float
        starttime/endtime of measurement window
    maxit: int
        no. of iter-decon iterations 
    taper_ratio: float
        taper of the window, default = 0.1

    Returns
    ----------------
    tr,am : float
        misfit function (am = tr)
    win: np.ndarray, shape(20)
        measure_adj window
    adj_r,adj_z: np.ndarray
        adjoint source on R/Z component shape(nt)
    cc_vohs: np.ndarray
        cross-correlation between observed vertical and synthetic horizontal data
    cc_hovs: np.ndarray
        cross-correlation between synthetic vertical and observed horizontal data
    """

    # make sure len(obs) == len(syn)
    nt = len(obs_v)
    assert len(obs_v) == len(syn_v), "Observed and synthetic data must have the same length"
    assert len(obs_h) == len(obs_v), "Observed_h and observed_r must have the same length"

    # taper 
    lpt,rpt,taper0 = taper_window(t0, dt, nt, tstart, tend, p=taper_ratio)

    # filter
    vsyn = syn_v[lpt:rpt] * taper0
    hsyn = syn_h[lpt:rpt] * taper0
    vobs = obs_v[lpt:rpt] * taper0
    hobs = obs_h[lpt:rpt] * taper0

    # compute convolution of h/v
    chi1 = convolve(vobs,hsyn,'full') * dt
    chi2 = convolve(vsyn,hobs,'full') * dt
    dchi = (chi1 - chi2)

    # misfit function
    mis = 0.5 * trapezoid(dchi**2,dx=dt)

    # reverse several arrays
    v_rev = vobs[::-1].copy()
    h_rev = hobs[::-1].copy()
    
    # adjoint source
    dchi_adj = dchi * 1.
    adj_z0:np.ndarray = -convolve(dchi_adj,h_rev,'valid') * dt
    adj_r0:np.ndarray = convolve(dchi_adj,v_rev,'valid') * dt

    adj_z = np.zeros_like(syn_v)
    adj_r = np.zeros_like(syn_v)
    adj_z[lpt:rpt] = adj_z0 * taper0
    adj_r[lpt:rpt] = adj_r0 * taper0

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
    cc1 = chi1 * 1.
    cc2 = chi2 * 1.

    return stats,adj_z,adj_r,cc1,cc2