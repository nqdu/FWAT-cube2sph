import numpy as np 
from typing import Final

def measure_adj(t0_inp,dt_inp,npt_inp,
                t0_syn,dt_syn,npt_syn,
                tstart,tend,imeas:int,
                tlong,tshort,verbose:bool,
                obs_data,syn_data,
                compute_adj_source = True,
                run_bandpass=False,display_details=False,
                output_measure_files = False,
                tshift_min = -4.5,tshift_max=4.5,
                dlna_min = -1.5,dlna_max=1.5,
                cc_min = 0.8,err_type = 1,
                dt_sigma_min=1.,dlna_sigma_min=0.5,
                itaper = 1,wtr=0.02,npi=2.5,
                dt_fac=2.,err_fac=2.5,
                dt_max_scale = 3.5,
                ncyle_in_window=1.5,
                use_physical_disp=False):
    """
    MEASURE ADJ wrapper
    """
    from .lib import libmeas
    
    if imeas == 5:
        itaper = 2

    tr,amp, window,adj =  \
        libmeas.measure(
            t0_inp,dt_inp,npt_inp,
            t0_syn,dt_syn,npt_syn,
            tstart,tend,imeas,
            tlong,tshort,verbose,
            obs_data,syn_data,run_bandpass,
            display_details,output_measure_files,
            compute_adj_source,tshift_min,
            tshift_max,dlna_min,dlna_max,
            cc_min,err_type,dt_sigma_min,
            dlna_sigma_min,itaper,wtr,npi,
            dt_fac,err_fac,dt_max_scale,
            ncyle_in_window,use_physical_disp
        )
    
    return tr,amp,window,adj

def measure_adj_exphase(obs,syn,t0,dt,nt,
               min_period,max_period,
               tstart,tend,water=0.1,
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
    from scipy.signal import hilbert
    from fwat.measure.utils import bandpass,taper_window
    from fwat.measure.utils import interpolate_syn
    from scipy.integrate import trapezoid

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
    w_s = water * np.max(Es)
    w_d = water * np.max(Ed)
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

def measure_adj_cross_conv(
    obs_v,syn_v,obs_h,syn_h,
    t0,dt,min_period,max_period,
    tstart,tend,maxit:int=120,taper_ratio = 0.1):

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

    Returns
    ----------------
    tr,am : float
        misfit function (am = tr)
    win: np.ndarray, shape(20)
        measure_adj window
    adj_r,adj_z: np.ndarray
        adjoint source on R/Z component shape(nt)
    """

    from fwat.measure.utils import bandpass,taper_window
    from fwat.measure.tele.deconit import time_decon
    from scipy.signal import convolve 
    from scipy.integrate import trapezoid

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

    return tr_chi,am_chi,win_chi,adj_r,adj_z

