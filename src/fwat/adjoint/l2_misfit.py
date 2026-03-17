import numpy as np 
from .MeasureStats import MeasureStats
from fwat.measure.utils import taper_window
from scipy.integrate import trapezoid

def measure_adj_l2(
        obs,syn,t0,dt,nt,
        tstart,tend,
        taper_ratio = 0.05
):
    """
    Compute the L2 misfit between observed and synthetic waveforms.

    Parameters
    ----------
    obs : numpy.ndarray
        Observed waveform.
    syn : numpy.ndarray
        Synthetic waveform.
    t0 : float
        Start time of the window.
    dt : float
        Time step.
    nt : int
        Number of time steps.
    tstart : float
        Start time of the measurement window.
    tend : float
        End time of the measurement window.
    taper_ratio : float
        Taper ratio for the windowing function.

    Returns
    -------
    stats: MeasureStats
        stats class containing measurement info
    adj: np.ndarray
        adjoint source, shape(nt)
    """

    # get window info
    lpt, rpt, taper0 = taper_window(t0, dt, nt, tstart, tend, p=taper_ratio)
    taper = syn * 0
    taper[lpt:rpt] = taper0 

    # get windowed data
    s = syn[lpt:rpt] * taper0 
    d = obs[lpt:rpt] * taper0 

    # compute misfit 
    ydiff = s - d 
    misfit = trapezoid(ydiff**2, dx=dt) * 0.5 

    # adjoint source and filter 
    adjsrc = obs * 0. 
    adjsrc[lpt:rpt] = ydiff 

    stats = MeasureStats(
        adj_type='l2',
        misfit=misfit,
        tstart=tstart,
        tend=tend,
        code='',
        tr_chi=misfit,
        am_chi=misfit,
        tshift=0.
    )

    return stats, adjsrc