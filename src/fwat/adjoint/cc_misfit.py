import numpy as np 
from numba import jit 

@jit(nopython=True)
def _cal_cc_correction(s,ishift,dlna):
    """
    Calculate the cross-correlation correction terms.

    Parameters
    ----------
    s : np.ndarray
        Input signal array.
    ishift : int
        Time shift for cross-correlation.
    dlna : float
        Logarithmic amplitude shift.

    Returns
    -------
    s_cc_dt : np.ndarray
        Cross-correlation correction for time.
    s_cc_dtdlna : np.ndarray
        Cross-correlation correction for amplitude.s
    """
    nlen = len(s)
    s_cc_dt = s * 0 
    s_cc_dtdlna = s * 0

    expdlna = np.exp(dlna)

    # correction
    for i in range(nlen):
        j = i - ishift

        if j >=0 and j < nlen:
            s_cc_dt[i] = s[j]
            s_cc_dtdlna[i] = expdlna * s[j]

    return s_cc_dt,s_cc_dtdlna


def _cal_cc_error(d,s,dt,ishift,dlna,dt_sigma_min,dlna_sigma_min):
    """
    Calculate the cross-correlation error terms.

    Parameters
    ----------
    d : np.ndarray
        Observed data.
    s : np.ndarray
        Synthetic data.
    dt : float
        Time sampling interval.
    ishift : int
        Time shift for cross-correlation.
    dlna : float
        Logarithmic amplitude shift.
    dt_sigma_min : float
        Minimum uncertainty for time shift.
    dlna_sigma_min : float
        Minimum uncertainty for amplitude shift.

    Returns
    -------
    sigma_dt : float
        Estimated uncertainty for time shift.
    sigma_dlna : float
        Estimated uncertainty for amplitude shift.
    """
    from fwat.measure.utils import dif1 

    # apply a scaling
    s_cc_dt,s_cc_dtdlna = _cal_cc_correction(s,ishift,dlna)

    # time derivative
    s_cc_vel = dif1(s_cc_dtdlna,dt)

    # The estimated error for dt and dlna with uncorrelation assumption
    sigma_dt_top = np.sum((d - s_cc_dtdlna)**2)
    sigma_dt_bot = np.sum(s_cc_vel**2)

    sigma_dlna_top = sigma_dt_top
    sigma_dlna_bot = np.sum(s_cc_dt**2)

    sigma_dt = np.sqrt(sigma_dt_top / sigma_dt_bot)
    sigma_dlna = np.sqrt(sigma_dlna_top / sigma_dlna_bot)
    
    # Check that errors do not go below the pre-defined threshold value
    if sigma_dt < dt_sigma_min or np.isnan(sigma_dt):
        sigma_dt = dt_sigma_min

    if sigma_dlna < dlna_sigma_min or np.isnan(sigma_dlna):
        sigma_dlna = dlna_sigma_min

    return sigma_dt, sigma_dlna

def _cc_shift(d,s,dt,dt_sigma_min,dlna_sigma_min,weight_by_uncertainty=True):
    """
    Calculate the cross-correlation shift.

    Parameters
    ----------
    d : np.ndarray
        Observed data.
    s : np.ndarray
        Synthetic data.
    dt : float
        Time sampling interval.
    dt_sigma_min : float
        Minimum uncertainty for time shift.
    dlna_sigma_min : float
        Minimum uncertainty for amplitude shift.
    weight_by_uncertainty : bool
        Whether to weight the shift by uncertainty.

    Returns
    -------
    tshift : float
        Estimated time shift.
    dlna : float
        Estimated logarithmic amplitude shift.
    sigma_dt : float
        Estimated uncertainty for time shift.
    sigma_dlna : float
        Estimated uncertainty for amplitude shift.
    """
    from scipy.signal import correlate
    
    # compute time shift between data and syn
    cc = correlate(d,s,'full')
    ishift = np.argmax(cc) - len(d) + 1
    tshift = ishift * dt 

    # compute dlna
    dlna = 0.5 * np.log(np.sum(d**2) / np.sum(s**2))

    # compute uncertainty
    if weight_by_uncertainty:
        sigma_dt,sigma_dlna = _cal_cc_error(d,s,dt,ishift,dlna,dt_sigma_min,dlna_sigma_min)
    else:
        sigma_dt = 1. 
        sigma_dlna = 1.

    return tshift,dlna,sigma_dt,sigma_dlna

def _cc_shift_dd(d1,s1,d2,s2,dt,dt_sigma_min,dlna_sigma_min):
    from scipy.signal import correlate

    # compute time shift between data and syn
    cc1 = correlate(d1,d2,'full')
    ishift_obs = np.argmax(cc1) - len(d1) + 1
    tshift_obs = ishift_obs * dt

    cc1 = correlate(s1,s2,'full')
    ishift_syn = np.argmax(cc1) - len(s1) + 1
    tshift_syn = ishift_syn * dt

    # overall shift
    ishift_dd = ishift_syn - ishift_obs
    tshift = ishift_dd * dt

    # compute dlna
    dlna_obs = 0.5 * np.log(np.sum(d1**2) / np.sum(d2**2))
    dlna_syn = 0.5 * np.log(np.sum(s2**2) / np.sum(s1**2))

    # uncertainties
    sigma_dt,sigma_dlna = _cal_cc_error(d1,d2,dt,ishift_obs,dlna_obs,dt_sigma_min,dlna_sigma_min)

    return tshift,tshift_obs,tshift_syn,dlna_obs,dlna_syn,sigma_dt,sigma_dlna

def measure_adj_cc(obs,syn,t0,dt,nt,
               min_period,max_period,
               tstart,tend,
               return_type = 'dt',
               taper_ratio = 0.05,
               tshift_min = -4.5,
               tshift_max = 4.5,
               weight_by_uncertainty = True,
               dt_sigma_min = 1.,
               dlna_sigma_min = 0.5):
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
    return_type: str
        type of return value, either 'dt' or 'am'
    taper_ratio: float
        taper of the window, default = 0.05
    tshift_min: float
        minimum time shift
    tshift_max: float
        maximum time shift
    weight_by_uncertainty: bool
        whether to weight the shift by uncertainty
    dt_sigma_min: float
        minimum uncertainty for time shift
    dlna_sigma_min: float
        minimum uncertainty for amplitude shift

    Returns
    ----------------
    tr,am : float
        misfit function (am = tr if return_type = dt)
    win: np.ndarray, shape(20)
        measure_adj window
    adj: np.ndarray
        adjoint source, shape(nt)
    """
    
    from fwat.measure.utils import bandpass,taper_window
    from fwat.measure.utils import dif1
    from scipy.integrate import trapezoid

    # make sure len(obs) == len(syn)
    assert len(obs) == len(syn), "Observed and synthetic data must have the same length"
    assert return_type in ['dt','am'], "return type should be one of dt,am"

    # get window info
    lpt, rpt, taper0 = taper_window(t0, dt, nt, tstart, tend, p=taper_ratio)
    taper = syn * 0
    taper[lpt:rpt] = taper0 

    # get windowed data
    s = syn[lpt:rpt] * taper0 
    d = obs[lpt:rpt] * taper0 

    # calculate time shift
    tshift,dlna,sigma_dt,sigma_dlna = _cc_shift(d,s,dt,dt_sigma_min,dlna_sigma_min,weight_by_uncertainty)

    # limit the tshift in [tshift_min,tshift_max]
    if tshift < tshift_min or tshift > tshift_max:
        tshift = 0.
        dlna = 0. 
        sigma_dt = 1.
        sigma_dlna = 1.

    # misfit
    misfit_p = 0.5 * (tshift / sigma_dt) **2 
    misfit_q = 0.5 * (dlna / sigma_dlna) ** 2  

    # calculate adjoint source for time
    dsdt = dif1(s,dt)
    nnorm = trapezoid(dsdt**2,dx=dt)
    fp = dsdt * tshift / nnorm / sigma_dt ** 2 

    # adjoint source for amplitude
    mnorm = trapezoid(s**2,dx=dt)
    fq = -1. * s * dlna / mnorm / sigma_dlna**2 

    # adjoint source
    adjsrc = obs * 0 
    tr_chi = misfit_p 
    am_chi = misfit_q 


    if return_type == 'dt':
        adjsrc[lpt:rpt] = fp * taper0
    else:
        adjsrc[lpt:rpt] = fq * taper0 
        tr_chi = am_chi 

    # filter adjoint source 
    adjsrc = bandpass(adjsrc,dt,1./max_period,1./min_period) * taper

    # measure_adj arrays
    win_chi = np.zeros((20))
    win_chi[6] = tshift
    win_chi[13-1] = 0.5 * np.sum( obs**2 )
    win_chi[14-1] = 0.5 * np.sum( syn**2 )
    win_chi[15-1] = tr_chi 
    win_chi[20-1] = nt * dt

    return tr_chi,am_chi,win_chi,adjsrc

def measure_adj_cc_dd(
        obs1,syn1,obs2,syn2,
        t0,dt,nt,
        min_period,max_period,
        tstart,tend,
        taper_ratio = 0.05,
        dt_sigma_min = 1.,
        dlna_sigma_min = 0.5):
    """
    Measure the double difference cross-correlation adjoint source

    Parameters
    ------------
    obs1: np.ndarray
        observed data for first pair, shape(nt)
    syn1: np.ndarray
        synthetic data for first pair, shape(nt)
    obs2: np.ndarray
        observed data for second pair, shape(nt)
    syn2: np.ndarray
        synthetic data for second pair, shape(nt)
    t0,dt,nt: float,float,int
        starttime/sampling/no.of points of adjoint source
    min/max_period: float
        minimum/maximum period used
    tstart,tend: float
        starttime/endtime of measurement window
    taper_ratio: float
        taper of the window, default = 0.05
    dt_sigma_min: float
        minimum uncertainty for time shift
    dlna_sigma_min: float
        minimum uncertainty for amplitude shift

    Returns
    ----------------
    tr,am : float
        misfit function
    win: np.ndarray, shape(20)
        measure_adj window
    adj: np.ndarray
        adjoint source, shape(nt)
    """

    from fwat.measure.utils import bandpass,taper_window
    from fwat.measure.utils import dif1
    from scipy.integrate import trapezoid

    # make sure len(obs) == len(syn)
    assert len(obs1) == len(syn1), "Observed and synthetic data must have the same length"

    # get window info
    lpt, rpt, taper0 = taper_window(t0, dt, nt, tstart, tend, p=taper_ratio)
    taper = syn1 * 0
    taper[lpt:rpt] = taper0 

    # get windowed data
    s1 = syn1[lpt:rpt] * taper0 
    d1 = obs1[lpt:rpt] * taper0 
    s2 = syn2[lpt:rpt] * taper0
    d2 = obs2[lpt:rpt] * taper0

    # calculate time shift
    tshift,_,_,dlna,_,sigma_dt,sigma_dlna = _cc_shift_dd(d1,s1,d2,s2,dt,dt_sigma_min,dlna_sigma_min)
    ishfit_dd = int(tshift / dt)

    # compute misfit 
    misfit_p = 0.5 * (tshift / sigma_dt) **2 
    misfit_q = 0.5 * (dlna / sigma_dlna) ** 2  

    # adjoint source pre-computing
    s1_cc_dt,_ = _cal_cc_correction(s1,-1. * ishfit_dd,0.)
    dsdt_cc1 = dif1(s1_cc_dt,dt)
    s2_cc_dt,_ = _cal_cc_correction(s2,ishfit_dd,0.)
    dsdt_cc2 = dif1(s2_cc_dt,dt)

    # norm
    ds1dt = dif1(s1,dt)
    nnorm = trapezoid(ds1dt*dsdt_cc2,dx=dt)
    fp_1 = -1 * dsdt_cc2 * tshift / nnorm / sigma_dt ** 2  # -1
    fp_2 = +1 * dsdt_cc1 * tshift / nnorm / sigma_dt ** 2  # +1

    # adjoint source
    adjsrc_1 = obs1 * 0 
    adjsrc_2 = obs2 * 0
    adjsrc_1[lpt:rpt] = fp_1 * taper0
    adjsrc_2[lpt:rpt] = fp_2 * taper0

    # filter
    adjsrc_1 = bandpass(adjsrc_1,dt,1./max_period,1./min_period) * taper
    adjsrc_2 = bandpass(adjsrc_2,dt,1./max_period,1./min_period) * taper

    #
    tr_chi = misfit_p
    am_chi = misfit_q
    win_chi = np.zeros((20))
    win_chi[6] = tshift
    win_chi[13-1] = 0.5 * np.sum( obs1**2 )
    win_chi[14-1] = 0.5 * np.sum( syn1**2 )
    win_chi[15-1] = tr_chi
    win_chi[20-1] = nt * dt

    return tr_chi,am_chi,win_chi,adjsrc_1,adjsrc_2
