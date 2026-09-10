
import  numpy as np 
from scipy.signal import hilbert
from .MeasureStats import MeasureStats
from fwat.measure.utils import bandpass,taper_window
from scipy.integrate import trapezoid
from .cc_misfit import cc_measure

def measure_adj_exphase(obs,syn,t0,dt,nt,
               min_period,max_period,
               tstart,tend,water=0.05,
               taper_ratio = 0.1,
               tshift_min = -4.5,
               tshift_max = 4.5,
               dlna_min = -1.5,
               dlna_max = 1.5,
               cc_min = 0.8):
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
    tshift_min: float
        minimum cross-correlation time shift
    tshift_max: float
        maximum cross-correlation time shift
    dlna_min: float
        minimum log amplitude shift
    dlna_max: float
        maximum log amplitude shift
    cc_min: float
        minimum cross-correlation coefficient to do measurement

    Returns
    ----------------
    stats: MeasureStats
        misfit and other measurement info
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

    # data selection, same rules as cc_misfit
    tshift,dlna,_,_,cc_coef = cc_measure(d,s,dt,1.,0.5,False)
    misfit_flag = 1.
    if tshift < tshift_min or tshift > tshift_max or  \
        cc_coef < cc_min or dlna < dlna_min or dlna > dlna_max:
        misfit_flag = 0.

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
    misfit = misfit * 0.5 * misfit_flag

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
    adjsrc = bandpass(adjsrc,dt,1./max_period,1./min_period) * taper * misfit_flag

    # measure_adj arrays
    stats = MeasureStats(
        adj_type = "exp_phase",
        misfit=misfit,
        tstart=tstart,
        tend=tend,
        tr_chi=misfit,
        am_chi=misfit,
        tshift = tshift
    )

    return stats, adjsrc
