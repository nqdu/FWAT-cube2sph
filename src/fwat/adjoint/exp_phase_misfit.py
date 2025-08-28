
import  numpy as np 
from scipy.signal import hilbert
from fwat.measure.utils import bandpass,taper_window
from scipy.integrate import trapezoid

def measure_adj_exphase(obs,syn,t0,dt,nt,
               min_period,max_period,
               tstart,tend,water=0.05,
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
    max_es = np.max(Es) 
    max_ed = np.max(Ed)
    w_s = water * max_es
    w_d = water * max_ed
    Es_wtr = Es + w_s
    Ed_wtr = Ed + w_d

    # get real/imag part
    dR = d / Ed_wtr - s / Es_wtr 
    dI = Hobs / Ed_wtr - Hsyn / Es_wtr

    # compute misfit
    misfit = trapezoid(dR**2,dx=dt) + trapezoid(dI**2,dx=dt)
    misfit = misfit * 0.5

    # adjoint source 
    Es_wtr_cubic = Es_wtr**3 
    adj_real = - (dR * Hsyn**2 / Es_wtr_cubic) - \
                np.imag(hilbert(dR * s * Hsyn / Es_wtr_cubic))
    
    adj_imag = dI * s * Hsyn / Es_wtr_cubic +  \
                np.imag(hilbert(dI * s**2 / Es_wtr_cubic))
    
    adjsrc = obs * 0 
    adjsrc[lpt:rpt] = (adj_imag + adj_real) * taper0 

    # filter
    taper = adjsrc * 0 
    taper[lpt:rpt] = taper0 
    adjsrc = bandpass(adjsrc,dt,1./max_period,1./min_period) * taper

    # measure_adj arrays
    tr_chi = misfit 
    am_chi = misfit 
    win_chi = np.zeros((20))
    win_chi[13-1] = 0.5 * np.sum( obs**2 )
    win_chi[14-1] = 0.5 * np.sum( syn**2 )
    win_chi[15-1] = tr_chi 
    win_chi[20-1] = nt * dt

    # reinterpolate adjoint source

    return tr_chi,am_chi,win_chi,adjsrc
