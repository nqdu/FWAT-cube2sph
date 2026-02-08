import numpy as np 
from fwat.measure.utils import bandpass,taper_window,dif1
from scipy.integrate import trapezoid
from fwat.adjoint.cc_misfit import cc_measure
from numba import jit 

def _prepare_mtm_data(d:np.ndarray,obs:np.ndarray,
                      tshift:float,dlna:float,
                      dt:float,nt:int,t0:float,
                      tstart:float,tend:float):
    """
    Prepare the data for multitaper measurement by shifting and scaling the observed data.

    Parameters
    --------------
    d: np.ndarray
        windowed data to be shifted and scaled, shape(nwin)
    obs: np.ndarray
        original observed data, shape(nt)
    tshift: float
        time shift to apply to the observed data
    dlna: float
        log amplitude shift to apply to the observed data
    dt: float
        sampling interval
    nt: int
        number of points in the original observed data
    t0: float
        start time of the original observed data
    tstart: float
        start time of the measurement window
    tend: float
        end time of the measurement window

    Returns
    --------------
    d: np.ndarray
        shifted and scaled data, shape(nwin)
    is_mtm: bool
        whether the data is still suitable for multitaper measurement after shifting and scaling
    """

    lpt,rpt,taper0 = taper_window(t0, dt, nt, tstart, tend, p=0.05)
    ishift = int(tshift / dt)

    # left/right sample
    lpt1 = max(lpt + ishift, 0)
    rpt1 = min(rpt + ishift, nt)
    nlen_d = rpt1 - lpt1

    if nlen_d == len(d):
        d = obs[lpt1:rpt1] * taper0
        d *= np.exp(-dlna) 
        is_mtm = True
    else:
        is_mtm = False

    return d, is_mtm

def _search_freq_limit(fnum:int,i_ampmax:int,scale_wl:float,
                        s_spec:np.ndarray, find_max: bool = False):
    # search for frequency limit with amplitude above water level
    is_search = True

    if not find_max:
        nfreq_limit = 0
        for iw in range(i_ampmax-1,-1,-1):
            if iw < i_ampmax:
                if abs(s_spec[iw]) < scale_wl and is_search:
                    is_search = False
                    nfreq_limit = iw 
                
                if abs(s_spec[iw]) >= scale_wl * 10. and not is_search:
                    is_search = True
                    nfreq_limit = iw
    else:
        nfreq_limit = fnum - 1
        for iw in range(i_ampmax+1,fnum):
            if abs(s_spec[iw]) < scale_wl and is_search:
                is_search = False
                nfreq_limit = iw 
            
            if abs(s_spec[iw]) >= scale_wl * 10. and not is_search:
                is_search = True
                nfreq_limit = iw        

    return is_search, nfreq_limit

def _cal_freq_limits(syn:np.ndarray,
                     df:float,dt:float,nlen_f:int,min_period:float,
                     max_period:float,water_level: float,
                     ncyle_in_window: float,npi: float):
    # calclute frequency limits for multitaper measurement
    fnum = int(nlen_f / 2 + 1 )
    s_spec = np.fft.rfft(syn,nlen_f) * dt 

    # max amplitute 
    ampmax = np.max(np.abs(s_spec))
    i_ampmax = int(np.argmax(np.abs(s_spec)))

    # scale factor 
    scale_wl = water_level * ampmax 

    # frequency limits
    ifmin = int(1. / max_period / df)
    ifmax = int(1. / min_period / df)

    # get max frequency index with amplitude above water level
    is_search, nfreq_max = _search_freq_limit(fnum,i_ampmax,scale_wl,s_spec,find_max=True)
    nfreq_max = min(nfreq_max, ifmax,int(1./(2.*dt)/df)-1)

    # get min frequency index with amplitude above water level
    is_search, nfreq_min = _search_freq_limit(fnum,i_ampmax,scale_wl,s_spec)
    
    # limit nfreq_min to at least N cycles in the window
    nfreq_min = max(nfreq_min, ifmin, int(ncyle_in_window / (len(syn) * dt * df)) - 1)
    half_taper_width = npi / (4.0 * len(syn) * dt)
    chosen_bandwidth = (nfreq_max - nfreq_min) * df

    if chosen_bandwidth < half_taper_width:
        nfreq_min = None 
        nfreq_max = None 
        is_mtm = False
    else:
        is_mtm = True

    return is_mtm, nfreq_min, nfreq_max

def measure_adj_mt(obs,syn,t0,dt,nt,
               min_period,max_period,
               tstart,tend,
               return_type = 'dt',
               taper_ratio = 0.05,
               tshift_min = -4.5,
               tshift_max = 4.5,
               dlna_min = -1.5,
               dlna_max=1.5,
               cc_min = 0.8,
               weight_by_uncertainty = True,
               dt_sigma_min = 1.,
               dlna_sigma_min = 0.5,
            wtr: float = 0.02,npi: float = 2.5,
            dt_fac: float = 2.,err_fac: float = 2.5,
            dt_max_scale: float = 3.5,
            ncyle_in_window: float = 1.5,):
    """
    Measure the multitaper adjoint source    

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
    dlna_min: float
        minimum log amplitude shift
    dlna_max: float
        maximum log amplitude shift
    cc_min: float
        minimum cross-correlation coefficient to do measurement
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
    tshift,dlna,sigma_dt,sigma_dlna,cc_coef = cc_measure(d,s,dt,dt_sigma_min,dlna_sigma_min,weight_by_uncertainty)

    # iterate
    is_mtm = True 
    while is_mtm:
        if tshift < dt:
            is_mtm = False
        
        if not is_mtm:
            break

        # shift and scale observed data
        d, is_mtm = _prepare_mtm_data(d, obs, tshift, dlna, dt, nt, t0, tstart, tend)
        if not is_mtm:
            break

        # determine fft info 
        lnpt = 15 
        nlen_f = 2**lnpt
        freq = np.fft.rfftfreq(nlen_f, dt)
        df = freq[1] - freq[0]
        wvec = np.pi * 2 * freq

        # check for sufficient cycles in the window
        is_mtm, nfreq_min, nfreq_max = _cal_freq_limits(s, df, dt, nlen_f, min_period, max_period, wtr, ncyle_in_window, npi)
        if not is_mtm:
            break

        # determine taper bandwidth