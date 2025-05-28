import numpy as np 

def read_params(paramfile="fwat_params/FWAT.PAR.yaml"):
    import yaml
    with open(paramfile,"r") as f:
        pdict = yaml.safe_load(f)
    
    return pdict

def measure_adj_file(t0:float,dt:float,nt:int,Tmin:float,Tmax:float,imeas:int,ccode:str):
    itaper = 1
    if imeas == 5:
        itaper = 2
    template = """ \
  %f %f   %d  # tstart, DT, npts: time vector for simulations
                      %d  # imeas (1-8; see manual)
                     %s  # channel: BH or LH
      %f      %f  # TLONG and TSHORT: band-pass periods for records
                .false.  # RUN_BANDPASS: use band-pass on records
                .false.  # DISPLAY_DETAILS
                .false.  # OUTPUT_MEASUREMENT_FILES
                 .true.  # COMPUTE_ADJOINT_SOURCE
     -4.5000     4.5000  # TSHIFT_MIN; TSHIFT_MAX
     -1.5000     1.5000  # DLNA_MIN; DLNA_MAX
                  0.800  # CC_MIN
                      1  # ERROR_TYPE -- 0 none; 1 CC, MT-CC; 2 MT-jack-knife
                  1.000  # DT_SIGMA_MIN
                  0.500  # DLNA_SIGMA_MIN
                      %d  # ITAPER -- taper type: 1 multi-taper; 2 cosine; 3 boxcar
            0.020  2.50  # WTR, NPI (ntaper = 2*NPI)
                  2.000  # DT_FAC
                  2.500  # ERR_FAC
                  3.500  # DT_MAX_SCALE
                  1.500  # NCYCLE_IN_WINDOW
		 .false. # USE_PHYSICAL_DISPERSION """% (t0,dt,nt,imeas,ccode,Tmax,Tmin,itaper)
    
    return template

def sac_cos_taper(npts,p):
    if p == 0.0 or p == 1.0:
        frac = int(npts * p / 2.0)
    else:
        frac = int(npts * p / 2.0 + 0.5)

    idx1 = 0
    idx2 = frac - 1
    idx3 = npts - frac
    idx4 = npts - 1

    # Very small data lengths or small decimal taper percentages can result in
    # idx1 == idx2 and idx3 == idx4. This breaks the following calculations.
    if idx1 == idx2:
        idx2 += 1
    if idx3 == idx4:
        idx3 -= 1

    # the taper at idx1 and idx4 equals zero and
    # at idx2 and idx3 equals one
    cos_win = np.zeros(npts,'f4')
    cos_win[idx1:idx2 + 1] = np.cos(-(
        np.pi / 2.0 * (float(idx2) -
                        np.arange(idx1, idx2 + 1)) / (idx2 - idx1)))
    cos_win[idx2 + 1:idx3] = 1.0
    cos_win[idx3:idx4 + 1] = np.cos((
        np.pi / 2.0 * (float(idx3) -
                        np.arange(idx3, idx4 + 1)) / (idx4 - idx3)))

    # if indices are identical division by zero
    # causes NaN values in cos_win
    if idx1 == idx2:
        cos_win[idx1] = 0.0
    if idx3 == idx4:
        cos_win[idx3] = 0.0
    return cos_win

def interpolate_syn(data,t1,dt1,npt1,t2,dt2,npt2,taper=0.05):
    """
    interpolate data from (t1, dt1, npt1) to a new data (t2,dt2,npt2)

    data: np.ndarray
        input data, shape(npt1) with starttime t1 and interval dt1 
    taper: float
        taper input data if required, taper * len(data) * 2 is the window used 
    
    """
    # taper input data if required
    data1 = np.float32(data)
    if taper > 0.:
        cos_tp = sac_cos_taper(npt1,taper)
        data1 = data1 * cos_tp

    temp = np.zeros((npt2),'f4')
    time = t2 + np.arange(npt2) * dt2 
    idx = np.logical_and(time > t1,time < t1 + (npt1-1) * dt1)
    ii = np.int64((time[idx] - t1) / dt1)
    tt = time[idx] - (ii * dt1 + t1) 
    temp[idx] = (data1[ii+1] - data1[ii]) * tt / dt1 + data1[ii]
    
    return temp

def preprocess(u,dt,freqmin,freqmax):
    import obspy
    tr = obspy.Trace(data = u)
    tr.stats.delta = dt 

    tr.detrend("demean")
    tr.detrend("linear")
    tr.taper(0.05)
    tr.filter("bandpass",freqmin=freqmin,freqmax=freqmax,zerophase=True,corners=4)
    tr.detrend("demean")
    tr.detrend("linear")
    tr.taper(0.05)

    w = np.float32(tr.data)

    return w

# diff function,central difference 1-st order 
def dif1(data,dt):
    n = len(data)
    data1 = np.zeros((n))
    data1[1:n-1] = 1. / (2 * dt) * (data[2:n] - data[0:n-2])

    return data1

def cumtrapz1(data,dt):
    from scipy.integrate import cumtrapz
    n = len(data)
    data1 = np.zeros((n))
    data1[1:n] = cumtrapz(data) * dt 

    return data1

def get_sem_seismo_info(sacfile:str):
    from obspy.io.sac import SACTrace
    syn_z_hd = SACTrace.read(sacfile,headonly=True)
    npt = syn_z_hd.npts
    dt = syn_z_hd.delta
    t0 = syn_z_hd.b 

    return t0,dt,npt

def rotate_NE_to_RT(vn:np.ndarray,ve:np.ndarray,bazd):
    """
    rotate NE components to RT, by using bazd
    
    Parameters
    ---------------
    vn,ve: np.ndarrray
        N/E componenets
    bazd: float
        back azimuth, in deg

    Returns:
    vr,vt: np.ndarray
        R/T components
    """
    from numpy import sin,cos 
    baz = np.deg2rad(bazd)
    vr = -ve * sin(baz) - vn * cos(baz)
    vt = -ve * cos(baz) + vn * sin(baz)

    return vr,vt

def rotate_RT_to_NE(vr:np.ndarray,vt:np.ndarray,bazd):
    """
    rotate NE components to RT, by using bazd
    
    Parameters
    ---------------
    vr,vt: np.ndarrray
        R/T componenets
    bazd: float
        back azimuth, in deg

    Returns:
    vn,ve: np.ndarray
        N/E components
    """
    from numpy import sin,cos 
    baz = np.deg2rad(bazd)
    ve = -vr * sin(baz) - vt * cos(baz)
    vn = -vr * cos(baz) + vt * sin(baz)

    return vn,ve