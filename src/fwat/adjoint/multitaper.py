"""
Simplified line-by-line translation of Fortran mt_adj subroutine
Only supports imeas=5 (cross-correlation) and imeas=7 (multitaper)
Dynamic sizing: NPT = nextpow2(npts), NDIM = npts
Disabled features: USE_PHYSICAL_DISPERSION, DO_RAY_DENSITY_SOURCE, DISPLAY_DETAILS
"""

import numpy as np
from scipy.signal.windows import dpss
from scipy.signal import hilbert

from fwat.measure.utils import bandpass
from fwat.adjoint.MeasureStats import MeasureStats

# Constants (from ma_constants.f90)
PI = np.pi
TWOPI = 2.0 * np.pi
WTR = 0.02
NCHI = 20


def nextpow2(n):
    """Return the next power of 2 greater than or equal to n"""
    return int(2 ** np.ceil(np.log2(n)))


def measure_adj_mt(obs, syn, t0, dt, npts, 
                   tstart, tend,
                   min_period, max_period,
                    return_type = 'dt',
                   tshift_min=-4.5, tshift_max=4.5,
                   dlna_min=-1.5, dlna_max=1.5,
                   cc_min=0.8,
                   dt_sigma_min=1.0, dlna_sigma_min=0.5,
                   itaper=1, wtr=0.02, npi=2.5,
                   dt_fac=2.0, err_fac=2.5,
                   dt_max_scale=3.5,
                   ncycle_in_window=1.5,
                   error_type=1):
    """
    High-level interface for multitaper adjoint source computation
    
    Parameters
    ----------
    obs, syn : ndarray
        Observed and synthetic data (npts,)
    t0 : float
        Start time of the trace
    dt : float
        Sampling interval
    npts : int
        Number of points in trace
    tstart, tend : float
        Start and end time of measurement window
    min_period, max_period : float
        Period band for measurements (seconds)
    return_type : str
        Type of measurement to return ('dt' for time shift, 'am' for amplitude)
    tshift_min, tshift_max : float
        Min/max allowed time shift
    dlna_min, dlna_max : float
        Min/max allowed amplitude anomaly
    cc_min : float
        Minimum cross-correlation coefficient
    dt_sigma_min, dlna_sigma_min : float
        Minimum uncertainties
    itaper : int
        Taper type (1=MT, 2=cosine, 3=boxcar)
    wtr : float
        Water level threshold
    npi : float
        Time-bandwidth product for multitaper
    dt_fac, err_fac, dt_max_scale : float
        Factors for MT measurement validation
    ncycle_in_window : float
        Minimum cycles in window
    error_type : int
        Error type (0=none, 1=CC, 2=jackknife)

    Returns
    -------
    tr_chi : float
        Traveltime misfit
    am_chi : float
        Amplitude misfit
    window_chi : ndarray
        Measurement array (20,)
    adj_src : ndarray
        Adjoint source (npts,)
    """
    
    # Validate inputs
    assert len(obs) == npts and len(syn) == npts, "Data length mismatch"
    assert tstart >= t0 and tend <= t0 + (npts-1)*dt, "Window out of bounds"
    
    # Set imeas based on return_type if not specified
    if return_type == 'dt':
        imeas = 5  # Cross-correlation for traveltime
    else:
        imeas = 7  # Multitaper for amplitude
    
    # Set is_mtm based on imeas and itaper
    if imeas == 5:
        is_mtm = itaper if itaper in [2, 3] else 2  # Cosine or boxcar for CC
    elif imeas == 7:
        is_mtm = 1  # Multitaper required for MT
    else:
        is_mtm = 2  # Default to cosine taper for unsupported imeas
    
    # Window and taper data
    istart = int((tstart - t0) / dt)
    iend = int((tend - t0) / dt) + 1
    nlen = iend - istart
    
    # Apply cosine taper to window
    ipwr_t = 10
    taper = 1.0 - np.cos(PI * np.arange(nlen) / (nlen - 1))**ipwr_t
    
    dat_dtw = obs[istart:iend] * taper
    syn_dtw = syn[istart:iend] * taper
    
    # Compute cross-correlation measurements
    tshift, dlna, sigma_dt_cc, sigma_dlnA_cc, cc_max = compute_cc_measurement(
        syn_dtw, dat_dtw, dt, dt_sigma_min, dlna_sigma_min
    )
    
    # Check CC quality
    if cc_max < cc_min:
        # Return zero adjoint if CC is too low
        window_chi = np.zeros(NCHI)
        window_chi[10] = sigma_dt_cc
        window_chi[11] = sigma_dlnA_cc
        adj_src = np.zeros(npts)
        
        stats = MeasureStats(
            adj_type='CC',
            tshift=tshift
        )
        return stats, adj_src
    
    # Check if measurements are within bounds
    if not (tshift_min <= tshift <= tshift_max):
        print(f"Warning: tshift {tshift:.3f} outside [{tshift_min}, {tshift_max}]")
    if not (dlna_min <= dlna <= dlna_max):
        print(f"Warning: dlna {dlna:.3f} outside [{dlna_min}, {dlna_max}]")
    
    # Deconstruct data using CC measurements
    ishift = int(tshift / dt)
    dat_dtw_cc = np.zeros(nlen)
    if ishift >= 0:
        dat_dtw_cc[ishift:] = dat_dtw[:nlen-ishift] * np.exp(-dlna)
    else:
        dat_dtw_cc[:nlen+ishift] = dat_dtw[-ishift:] * np.exp(-dlna)
    
    # Setup frequency parameters
    NPT = nextpow2(npts)
    df = 1.0 / (NPT * dt)
    fnum = NPT // 2 + 1
    
    # Frequency limits
    i_left = max(1, int(ncycle_in_window / (nlen * dt * df)))
    i_right = min(fnum - 1, int(1.0 / (min_period * df)))
    
    # Initialize frequency-dependent measurements
    dtau_w = np.zeros(NPT)
    dlnA_w = np.zeros(NPT)
    err_dtau = np.ones(NPT) * sigma_dt_cc
    err_dlnA = np.ones(NPT) * sigma_dlnA_cc
    sigma_dt = sigma_dt_cc
    sigma_dlnA = sigma_dlnA_cc
    
    # For multitaper, compute frequency-dependent measurements
    if is_mtm == 1 and imeas == 7:
        trans_mtm, i_pmax_syn, i_right = compute_mt_transfer_function(
            syn_dtw, dat_dtw_cc, nlen, dt, npi, wtr, i_left, i_right
        )
        
        # Extract measurements from transfer function
        dtau_w, dlnA_w = extract_mt_measurements(
            trans_mtm, df, NPT, tshift, dlna, i_right
        )
        
        # Validate MT measurements
        use_mt = validate_mt_measurements(
            dtau_w, dlnA_w, i_left, i_right, dt, df, nlen, 
            dt_fac, dt_max_scale, tshift, ncycle_in_window, min_period
        )
        
        if not use_mt:
            # Fall back to CC
            #print("MT measurement failed validation, falling back to CC")
            imeas = 5
            is_mtm = 2
    
    # Initialize window_chi
    window_chi = np.zeros(NCHI)
    
    # Compute adjoint source if requested
    tr_adj_src, am_adj_src, tr_chi, am_chi, window_chi = mt_adj_simplified(
        istart, dat_dtw, syn_dtw, nlen, dt, tshift, dlna,
        sigma_dt_cc, sigma_dlnA_cc, dtau_w, dlnA_w, 
        err_dtau, err_dlnA, sigma_dt, sigma_dlnA,
        i_left, i_right, window_chi, npts,
        imeas=imeas, is_mtm=is_mtm, ERROR_TYPE=error_type, NPI=npi
    )
    
    # Extract adjoint source for appropriate measurement type
    misfit = 0.
    if return_type == 'dt':
        adj_src = tr_adj_src  # Traveltime adjoint
        misfit = tr_chi
    else:
        adj_src = am_adj_src  # Amplitude adjoint
        misfit = am_chi

    # bandpass adj_src to measurement band
    adj_src_bp = bandpass(adj_src, dt, 1./max_period, 1./min_period)
    taper_all = adj_src_bp * 0
    taper_all[istart:iend] = taper
    adj_src = taper_all * adj_src_bp

    # CREATE MeasureStats object for output
    adj_type = 'CC' if imeas == 5 else 'MT'
    stats = MeasureStats(
        adj_type=adj_type,
        misfit=misfit,
        tr_chi=tr_chi,
        am_chi=am_chi,
        tshift=tshift
    )


    return stats, adj_src


def compute_cc_measurement(syn, dat, dt, dt_sigma_min, dlna_sigma_min):
    """Compute cross-correlation time shift and amplitude measurements"""
    nlen = len(syn)
    
    # Cross-correlation
    cc = np.correlate(dat, syn, mode='full')
    lag = np.arange(-nlen+1, nlen)
    
    # Find maximum correlation
    imax = np.argmax(np.abs(cc))
    cc_max = cc[imax] / (np.sqrt(np.sum(dat**2) * np.sum(syn**2)) + 1e-10)
    
    # Time shift
    tshift = lag[imax] * dt
    
    # Amplitude measurement using Hilbert transform
    syn_env = np.abs(hilbert(syn))
    dat_env = np.abs(hilbert(dat))
    
    # Log amplitude ratio
    dlna = 0.5 * np.log(np.sum(dat_env**2) / (np.sum(syn_env**2) + 1e-10) + 1e-10)
    
    # Uncertainty estimates (simplified)
    Tmax = nlen * dt
    sigma_dt_cc = max(dt_sigma_min, np.sqrt(1.0 - cc_max**2) * Tmax / (PI * cc_max + 1e-10))
    sigma_dlnA_cc = max(dlna_sigma_min, np.sqrt(1.0 - cc_max**2) / (2.0 * cc_max + 1e-10))
    
    return tshift, dlna, sigma_dt_cc, sigma_dlnA_cc, cc_max


def compute_mt_transfer_function(syn, dat, nlen, dt, npi, wtr, i_left, i_right):
    """Compute multitaper transfer function"""
    NPT = nextpow2(len(syn) * 2)
    ntaper = int(npi * 2.0)
    
    # Get DPSS tapers
    tapers, eigens = dpss(nlen, npi, Kmax=ntaper, return_ratios=True)
    ntaper = tapers.shape[0]
    
    # Initialize transfer function terms
    fnum = NPT // 2 + 1
    top_mtm = np.zeros(fnum, dtype=complex)
    bot_mtm = np.zeros(fnum, dtype=complex)
    
    # Find peak frequency in synthetics
    syn_fft = np.fft.rfft(syn, n=NPT) * dt
    i_pmax_syn = np.argmax(np.abs(syn_fft))
    
    # Compute transfer function components
    for ictaper in range(ntaper):
        # Apply taper
        syn_h = np.zeros(NPT)
        dat_h = np.zeros(NPT)
        syn_h[:nlen] = syn * tapers[ictaper, :]
        dat_h[:nlen] = dat * tapers[ictaper, :]
        
        # FFT
        syn_ho = np.fft.rfft(syn_h, n=NPT) * dt
        dat_ho = np.fft.rfft(dat_h, n=NPT) * dt
        
        # Accumulate
        top_mtm += dat_ho * np.conj(syn_ho)
        bot_mtm += syn_ho * np.conj(syn_ho)
    
    # Water level
    ampmax = np.max(np.abs(bot_mtm))
    wtr_use = ampmax * wtr**2
    
    # Compute transfer function
    trans_mtm = np.zeros(fnum, dtype=complex)
    mask = np.abs(bot_mtm) > wtr_use
    trans_mtm[mask] = top_mtm[mask] / bot_mtm[mask]
    trans_mtm[~mask] = top_mtm[~mask] / (bot_mtm[~mask] + wtr_use)
    
    # Update i_right based on water level
    for i in range(i_pmax_syn, fnum):
        if np.abs(bot_mtm[i]) <= wtr_use:
            i_right = min(i, i_right)
            break
    
    return trans_mtm, i_pmax_syn, i_right


def extract_mt_measurements(trans_mtm, df, NPT, tshift_cc, dlna_cc, i_right):
    """Extract phase and amplitude from transfer function"""
    dtau_w = np.zeros(NPT)
    dlnA_w = np.zeros(NPT)
    
    # Phase unwrapping and conversion to time shift
    phi = np.angle(trans_mtm)
    
    # Unwrap phase
    phi_unwrap = np.unwrap(phi[:i_right+1])
    
    # Convert to time shift: dtau = -phi / (2*pi*f)
    for i in range(1, i_right+1):
        freq = i * df
        dtau_w[i] = -phi_unwrap[i] / (TWOPI * freq) + tshift_cc
    
    # Amplitude: dlnA = 0.5 * log(|trans|)
    dlnA_w[:i_right+1] = 0.5 * np.log(np.abs(trans_mtm[:i_right+1]) + 1e-10) + dlna_cc
    
    return dtau_w, dlnA_w


def validate_mt_measurements(dtau_w, dlnA_w, i_left, i_right, dt, df, nlen,
                             dt_fac, dt_max_scale, tshift, ncycle_in_window, min_period):
    """Validate multitaper measurements"""
    
    # Check if time shifts are reasonable
    for i in range(i_left, i_right+1):
        freq = i * df
        period = 1.0 / freq
        
        # Check 1: dt should not be larger than a fraction of the period
        if np.abs(dtau_w[i]) > period / dt_fac:
            return False
        
        # Check 2: dt should not be much larger than CC measurement
        if np.abs(dtau_w[i]) > dt_max_scale * abs(tshift):
            return False
    
    # Check 3: sufficient cycles in window
    if ncycle_in_window * min_period > nlen * dt:
        return False
    
    return True


def mt_adj_simplified(istart: int, 
                      dat_dtw: np.ndarray, 
                      syn_dtw: np.ndarray, 
                      nlen: int, 
                      dt: float, 
                      tshift: float, 
                      dlnA: float, 
                      sigma_dt_cc: float, 
                      sigma_dlnA_cc: float, 
                      dtau_w: np.ndarray, 
                      dlnA_w: np.ndarray, 
                      err_dtau: np.ndarray, 
                      err_dlnA: np.ndarray, 
                      sigma_dt: float, 
                      sigma_dlnA: float, 
                      i_left: int, 
                      i_right: int, 
                      window_chi: np.ndarray, 
                      npts_orig: int,
                      imeas: int = 7, 
                      is_mtm: int = 1, 
                      ERROR_TYPE: int = 1, 
                      NPI: float = 2.5):
    """
    Simplified line-by-line translation of Fortran mt_adj subroutine
    Only supports imeas=5 (CC) and imeas=7 (MT)
    
    Input:
        istart: starting index of the windowed portion in original trace
        dat_dtw: windowed data (nlen)
        syn_dtw: windowed synthetic (nlen)
        nlen: length of windowed data
        dt: sampling rate
        tshift: CC time shift measurement
        dlnA: CC amplitude measurement
        sigma_dt_cc: CC time shift uncertainty
        sigma_dlnA_cc: CC amplitude uncertainty
        dtau_w: frequency-dependent time shift (NPT)
        dlnA_w: frequency-dependent amplitude (NPT)
        err_dtau: error estimates for dtau (NPT)
        err_dlnA: error estimates for dlnA (NPT)
        sigma_dt: average time shift error
        sigma_dlnA: average amplitude error
        i_left: left frequency index
        i_right: right frequency index
        window_chi: measurement array (NCHI)
        npts_orig: original number of points in full trace (for output sizing)
        imeas: measurement type (5=CC, 7=MT)
        is_mtm: multitaper flag (1 for MT, 2 for cosine taper)
        ERROR_TYPE: error type (0=none, 1=CC, 2=jackknife)
        NPI: time-bandwidth product
    
    Output:
        tr_adj_src: time shift adjoint source (npts)
        am_adj_src: amplitude adjoint source (npts)
        tr_chi: time shift misfit
        am_chi: amplitude misfit
        window_chi: updated measurement array
    """
    
    # Dynamic sizing based on actual data length
    # Use original npts for output array sizing to match input obs/syn
    NPT = nextpow2(npts_orig)
    NDIM = npts_orig
    
    # Define post-processing time-domain taper (boxcar window)
    time_window = np.zeros(NPT)
    time_window[:nlen] = 1.0  # boxcar window
    
    # CROSS CORRELATION ADJOINT SOURCES
    # Initialize arrays
    ft_bar_t = np.zeros(NPT)
    fa_bar_t = np.zeros(NPT)
    syn_vtw = np.zeros(NPT)
    
    if imeas == 5:
        # Compute synthetic velocity (no physical dispersion)
        syn_vtw[1:nlen-1] = (syn_dtw[2:nlen] - syn_dtw[0:nlen-2]) / (2.0*dt)
        syn_vtw[0] = (syn_dtw[1] - syn_dtw[0]) / dt
        syn_vtw[nlen-1] = (syn_dtw[nlen-1] - syn_dtw[nlen-2]) / dt
        
        # Compute CC traveltime and amplitude banana-doughnut kernels
        Nnorm = dt * np.sum(syn_vtw[:nlen] * syn_vtw[:nlen])
        ft_bar_t[:nlen] = -syn_vtw[:nlen] / Nnorm
        
        Mnorm = dt * np.sum(syn_dtw[:nlen] * syn_dtw[:nlen])
        fa_bar_t[:nlen] = syn_dtw[:nlen] / Mnorm
    
    # FREQUENCY-DOMAIN TAPERS FOR MT ADJOINT SOURCES
    fp = np.zeros(NPT)
    fq = np.zeros(NPT)
    
    if is_mtm == 1 and imeas == 7:
        # Initialize water levels for err_dtau/dlnA division
        dtau_wtr = WTR * np.sum(np.abs(dtau_w[i_left:i_right+1])) / (i_right - i_left + 1)
        dlnA_wtr = WTR * np.sum(np.abs(dlnA_w[i_left:i_right+1])) / (i_right - i_left + 1)
        
        # Frequency-domain tapers (cosine taper)
        ipwr_w = 10
        w_taper = np.zeros(NPT)
        i_range = np.arange(i_left, i_right + 1)
        w_taper[i_left:i_right+1] = 1.0 - np.cos(PI * (i_range - i_left) / (i_right - i_left))**ipwr_w
        
        # Compute normalization factor for w_taper
        df = 1.0 / (NPT * dt)
        ffac = 2.0 * df * np.sum(w_taper[i_left:i_right+1])
        
        # wp_taper and wq_taper are modified frequency-domain tapers
        wp_taper = np.zeros(NPT)
        wq_taper = np.zeros(NPT)
        idx = np.arange(i_left, i_right + 1)
        
        if ERROR_TYPE == 0:
            # No error estimate, only adds normalization factor
            wp_taper[idx] = w_taper[idx] / ffac
            wq_taper[idx] = w_taper[idx] / ffac
            
        elif ERROR_TYPE == 1:
            # CC error as a constant for all freqs
            wp_taper[idx] = w_taper[idx] / ffac / (sigma_dt ** 2)
            wq_taper[idx] = w_taper[idx] / ffac / (sigma_dlnA ** 2)
            
        elif ERROR_TYPE == 2:
            # MT jack-knife error estimate
            err_t = err_dtau[idx].copy()
            mask_t = err_dtau[idx] < dtau_wtr
            err_t[mask_t] = err_t[mask_t] + dtau_wtr
            
            err_A = err_dlnA[idx].copy()
            mask_A = err_dlnA[idx] < dlnA_wtr
            err_A[mask_A] = err_A[mask_A] + dlnA_wtr
            
            wp_taper[idx] = w_taper[idx] / ffac / (err_t ** 2)
            wq_taper[idx] = w_taper[idx] / ffac / (err_A ** 2)
        
        # Allocate MT variables
        ntaper = int(NPI * 2.0)
        
        # Get the MT tapers
        tapers, eigens = dpss(nlen, NPI, Kmax=ntaper, return_ratios=True)
        ntaper = tapers.shape[0]
        # tapers has shape (ntaper, nlen), transpose for column-wise indexing
        tas = np.zeros((NPT, ntaper))
        for ictaper in range(ntaper):
            tas[:nlen, ictaper] = tapers[ictaper, :]
        
        # Initialize (use rfft output size)
        nfreq = NPT // 2 + 1
        d_bot_mtm = np.zeros(nfreq, dtype=complex)
        v_bot_mtm = np.zeros(nfreq, dtype=complex)
        
        # Compute the bot required to compute p_j's and q_j's
        syn_dtw_ho_all = np.zeros((nfreq, ntaper), dtype=complex)
        syn_vtw_ho_all = np.zeros((nfreq, ntaper), dtype=complex)
        
        for ictaper in range(ntaper):
            syn_dtw_h = np.zeros(NPT)
            syn_vtw_h = np.zeros(NPT)
            
            # Tapered synthetic displacement (no physical dispersion)
            syn_dtw_h[:nlen] = syn_dtw[:nlen] * tas[:nlen, ictaper]
            
            # Compute velocity of tapered syn
            syn_vtw_h[1:nlen-1] = (syn_dtw_h[2:nlen] - syn_dtw_h[0:nlen-2]) / (2.0*dt)
            syn_vtw_h[0] = (syn_dtw_h[1] - syn_dtw_h[0]) / dt
            syn_vtw_h[nlen-1] = (syn_dtw_h[nlen-1] - syn_dtw_h[nlen-2]) / dt
            
            # Apply RFFT to get complex spectra
            syn_dtw_ho_all[:, ictaper] = np.fft.rfft(syn_dtw_h, n=NPT) * dt
            syn_vtw_ho_all[:, ictaper] = np.fft.rfft(syn_vtw_h, n=NPT) * dt
            
            d_bot_mtm[:] = d_bot_mtm[:] + syn_dtw_ho_all[:, ictaper] * np.conj(syn_dtw_ho_all[:, ictaper])
            v_bot_mtm[:] = v_bot_mtm[:] + syn_vtw_ho_all[:, ictaper] * np.conj(syn_vtw_ho_all[:, ictaper])
        
        # Compute p_j, q_j, P_j, Q_j and adjoint source fp, fq
        for ictaper in range(ntaper):
            # Compute p_j(w) and q_j(w)
            pwc_adj = np.zeros(nfreq, dtype=complex)
            qwc_adj = np.zeros(nfreq, dtype=complex)
            
            idx = np.arange(1, min(i_right + 1, nfreq))  # Skip DC, stay within rfft bounds
            pwc_adj[idx] = syn_vtw_ho_all[idx, ictaper] / v_bot_mtm[idx]
            qwc_adj[idx] = -syn_dtw_ho_all[idx, ictaper] / d_bot_mtm[idx]
            
            # Compute P_j(w) and Q_j(w) - no ray density, add misfit measurement
            pwc_adj[:] = pwc_adj[:] * (dtau_w[:nfreq] + 0j) * (wp_taper[:nfreq] + 0j)
            qwc_adj[:] = qwc_adj[:] * (dlnA_w[:nfreq] + 0j) * (wq_taper[:nfreq] + 0j)
            
            # IRFFT into the time domain
            dtau_pj_t = np.fft.irfft(pwc_adj, n=NPT) / dt
            dlnA_qj_t = np.fft.irfft(qwc_adj, n=NPT) / dt
            
            # Create adjoint source: applies taper to time signal
            fp[:] = fp[:] + tas[:, ictaper] * dtau_pj_t[:]
            fq[:] = fq[:] + tas[:, ictaper] * dlnA_qj_t[:]
    
    # ASSEMBLE VARIOUS ADJOINT SOURCES
    tr_adj_src = np.zeros(NDIM)
    am_adj_src = np.zeros(NDIM)
    
    # Integrated waveform difference squared
    dat_win = dat_dtw[:nlen] * time_window[:nlen]
    syn_win = syn_dtw[:nlen] * time_window[:nlen]
    waveform_d2 = np.sum(dat_win**2)
    waveform_s2 = np.sum(syn_win**2)
    waveform_chi = np.sum((dat_win - syn_win)**2)
    
    # Compute traveltime and amplitude adjoint sources for imeas
    idx = np.arange(istart, istart + nlen)
    
    if imeas == 5:
        # Cross-correlation (no ray density)
        tr_adj_src[idx] = -(tshift / sigma_dt_cc**2) * ft_bar_t[:nlen] * time_window[:nlen]
        am_adj_src[idx] = -(dlnA / sigma_dlnA_cc**2) * fa_bar_t[:nlen] * time_window[:nlen]
    
    elif imeas == 7:
        # Multitaper
        tr_adj_src[idx] = fp[:nlen] * time_window[:nlen]
        am_adj_src[idx] = fq[:nlen] * time_window[:nlen]
    
    # COMPUTE MISFIT FUNCTION VALUE
    if is_mtm == 1 and imeas == 7:
        window_chi[0] = 0.5 * 2.0 * df * np.sum((dtau_w[:i_right+1])**2 * wp_taper[:i_right+1])
        window_chi[1] = 0.5 * 2.0 * df * np.sum((dlnA_w[:i_right+1])**2 * wq_taper[:i_right+1])
    
    window_chi[2] = 0.5 * (tshift / sigma_dt_cc)**2
    window_chi[3] = 0.5 * (dlnA / sigma_dlnA_cc)**2
    
    # CC/averaged mt measurement (no uncertainty estimates)
    if is_mtm == 1 and imeas == 7:
        window_chi[4] = np.sum(dtau_w[:i_right+1] * w_taper[:i_right+1]) / np.sum(w_taper[:i_right+1])
        window_chi[5] = np.sum(dlnA_w[:i_right+1] * w_taper[:i_right+1]) / np.sum(w_taper[:i_right+1])
    
    window_chi[6] = tshift
    window_chi[7] = dlnA
    
    # Estimated measurement uncertainties
    if is_mtm == 1 and imeas == 7:
        window_chi[8] = sigma_dt
        window_chi[9] = sigma_dlnA
    
    window_chi[10] = sigma_dt_cc
    window_chi[11] = sigma_dlnA_cc
    
    # For normalization, divide by duration of window
    window_chi[12] = 0.5 * waveform_d2
    window_chi[13] = 0.5 * waveform_s2
    window_chi[14] = 0.5 * waveform_chi
    window_chi[15] = nlen * dt
    
    # Compute tr_chi and am_chi based on imeas
    if imeas == 5:
        # Cross-correlation
        tr_chi = window_chi[2]
        am_chi = window_chi[3]
    elif imeas == 7:
        # Multitaper
        tr_chi = window_chi[0]
        am_chi = window_chi[1]
    else:
        tr_chi = 0.0
        am_chi = 0.0
    
    return tr_adj_src, am_adj_src, tr_chi, am_chi, window_chi
