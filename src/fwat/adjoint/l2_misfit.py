import numpy as np 

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
    tr,am : float
        misfit function (am = tr if return_type = dt)
    win: np.ndarray, shape(20)
        measure_adj window
    adj: np.ndarray
        adjoint source, shape(nt)
    """
    from fwat.measure.utils import taper_window
    from scipy.integrate import trapezoid

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

    # measure_adj arrays
    tr_chi = misfit 
    am_chi = misfit 
    win_chi = np.zeros((20))
    win_chi[13-1] = 0.5 * np.sum( obs**2 )
    win_chi[14-1] = 0.5 * np.sum( syn**2 )
    win_chi[15-1] = tr_chi 
    win_chi[20-1] = nt * dt

    return tr_chi,am_chi,win_chi,adjsrc