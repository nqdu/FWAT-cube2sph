import numpy as np 
from fwat.const import PARAM_FILE

def alloc_mpi_jobs(ntasks,nprocs,myrank):
    sub_n = ntasks // nprocs
    num_larger_procs = ntasks - nprocs * sub_n
    startid = 0
    if myrank < num_larger_procs:
        sub_n = sub_n + 1
        startid = 0 + myrank * sub_n
    elif sub_n > 0 : 
        startid = 0 + num_larger_procs + myrank * sub_n
    else : #// this process has only zero elements
        startid = -1
        sub_n = 0
    
    endid = startid + sub_n - 1

    return startid,endid

def read_params(paramfile:str):
    import yaml
    with open(paramfile,"r") as f:
        pdict = yaml.safe_load(f)
    
    return pdict

def hann_taper(npts,p):
    from scipy.signal import windows
    max_len = int(npts * p)
    max_len = min(max_len, int(npts / 2))

    # get hann window
    wlen = max_len * 2 
    if wlen < npts:
        wlen += 1

    win = np.ones(npts)
    hanwindow = windows.hann(wlen)
    win[:max_len] = hanwindow[:max_len]
    win[npts-max_len:] = hanwindow[wlen-max_len:]

    return win

def sac_cos_taper(npts,p):
    frac = int(npts * p) 
    frac = min(frac,int(npts/2) - 1)

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

_TAPER_ENTRY_POINT = {
    "hann": hann_taper,
    "cos": sac_cos_taper
}

def interpolate_syn(data,t1,dt1,npt1,t2,dt2,npt2,max_percentage=0.05,type_='hann'):
    """
    interpolate data from (t1, dt1, npt1) to a new data (t2,dt2,npt2)

    data: np.ndarray
        input data, shape(npt1) with starttime t1 and interval dt1 
    max_percentage: float
        max_percentage of taper on one side, if required, taper * len(data) * 2 is the window used 
    
    """
    # taper input data if required
    data1 = data * 1.
    if max_percentage > 0.:
        func = _TAPER_ENTRY_POINT[type_]
        cos_tp = func(npt1,max_percentage)
        data1 = data1 * cos_tp

    temp = np.zeros((npt2))
    time = t2 + np.arange(npt2) * dt2 
    idx = np.logical_and(time > t1,time < t1 + (npt1-2) * dt1)
    ii = np.int64((time[idx] - t1) / dt1)
    tt = time[idx] - (ii * dt1 + t1) 
    temp[idx] = (data1[ii+1] - data1[ii]) * tt / dt1 + data1[ii]
    
    return temp


def bandpass(u,dt,freqmin,freqmax,max_percentage=0.05,type_='hann'):
    """
    Apply a SAC-like bandpass filter with pre-processing to a 1D signal.

    This function performs the following steps:
        1. Detrends and demeans the input signal to remove linear trends and DC offset.
        2. Applies a Hann taper to the signal to reduce edge effects.
        3. Applies a zero-phase Butterworth bandpass filter between `freqmin` and `freqmax`.
        4. taper again to reduce edge effects after filtering.

    Parameters
    ----------
    u : array_like
        Input 1D signal to be filtered.
    dt : float
        Sampling interval (in seconds).
    freqmin : float
        Lower frequency bound of the bandpass filter (in Hz).
    freqmax : float
        Upper frequency bound of the bandpass filter (in Hz).
    max_percentage : float
        Maximum percentage of taper to apply at both ends of the signal.
    type_ : str
        Type of tapering window to use ('hann' or 'cos').

    Returns
    -------
    u1 : ndarray
        The filtered signal after bandpass filtering and pre/post-processing.

    Notes
    -----
    - The function uses a 4th-order Butterworth filter implemented with second-order sections (SOS) for numerical stability.
    - The Hann taper is applied with a default fraction of 0.05 (5%) at both ends of the signal.
    - The function assumes the existence of a `hann_taper` function for tapering.
    """
    from scipy import signal 
    assert type_ in ['hann','cos'], "type_ should be one of hann/cos"

    # detrend/demean
    u1 = u - np.mean(u)
    u1 = signal.detrend(u1,type='linear')

    # taper 
    func = _TAPER_ENTRY_POINT[type_]
    win = func(len(u1),max_percentage)
    u1 = u1 * win

    # filter
    sos = signal.butter(4, [freqmin,freqmax], btype='bandpass', fs=1.0/dt, output='sos')
    u1 = signal.sosfiltfilt(sos, u1)

    # taper again
    win = func(len(u1),max_percentage * 0.5)
    u1 = u1 * win

    return u1

# diff function,central difference 1-st order 
def dif1(data,dt):
    n = len(data)
    data1 = np.zeros((n))
    data1[1:n-1] = 1. / (2 * dt) * (data[2:n] - data[0:n-2])

    return data1

def cumtrapz1(data,dt):
    from scipy.integrate import cumulative_trapezoid
    data1 = cumulative_trapezoid(data,dx=dt,initial=0.)

    return data1

def get_window_info(t0,dt,tstart,tend):
    """
    get window start/end sample index

    Parameters
    -----------
    t0,dt: float
        start time and sampling interval
    tstart,tend: float
        start/end time of this window

    Returns
    -----------
    left_pt,right_pt: int
        the window sample index [left_pt,right_pt)
    """
    assert tend >= tstart, "window is reversed!"

    nlen = int((tend - tstart) / dt) + 1
    left_pt = int(np.floor(tstart - t0) / dt)
    right_pt = left_pt + nlen

    return left_pt,right_pt 

def taper_window(t0,dt,nt,tstart,tend,p=0.05,type_='hann'):
    """
    get window start/end sample index

    Parameters
    -----------
    t0,dt: float
        start time and sampling interval
    nt: int
        no. of points in total
    tstart,tend: float
        start/end time of this window
    p: float
        ratio of data to be tapered on one side
    type_: str
        type of tapering window, either 'hann' or 'cos'

    Returns
    -----------
    left_pt,right_pt: int
        the window sample index [left_pt,right_pt)
    """
    assert tend >= tstart, "window is reversed!"
    assert type_ in ['hann','cos'], "type_ should be one of hann/cos"

    nlen = int((tend - tstart) / dt) + 1
    left_pt = int(np.floor(tstart - t0) / dt)
    if left_pt < 0:
        left_pt = 0
    right_pt = left_pt + nlen
    if right_pt > nt:
        right_pt = nt
    nlen = right_pt - left_pt

    func = _TAPER_ENTRY_POINT[type_]
    win = func(nlen,p)

    return left_pt,right_pt,win 

def get_source_loc(evtid:str,sourcefile:str):
    """
    get source location evla,evlo,evdp

    Parameters
    ----------------
    evtid: str
        source tag
    sourcefile: str
        source file, like sources.dat.tele

    Returns
    -------------
    evla,evlo,evdp: lat/lon/depth (in km) of this source
    """
    # load source file in txt 
    temp = np.loadtxt(sourcefile,dtype=str,ndmin=2)

    # loop every source to find travel time
    find_source = False
    evla,evlo,evdp = 0,0,0
    for i in range(temp.shape[0]):
        if temp[i,0] == evtid:
            evla,evlo,evdp = float(temp[i,1]),float(temp[i,2]),float(temp[i,3])
            find_source = True
    
    if not find_source:
        print(f"cannot find source, please check {evtid} and {sourcefile}")
        exit(1)
    
    return evla,evlo, evdp 

def get_simu_info(sacfile:str):
    from obspy.io.sac import SACTrace
    syn_z_hd = SACTrace.read(sacfile,headonly=True)
    npt = syn_z_hd.npts
    dt = syn_z_hd.delta
    t0 = syn_z_hd.b 

    return t0,dt,npt

def _geod2geoc(geographic_lat_deg, flattening=0.0033528106647474805):
    """
    Convert geographic (geodetic) latitude to geocentric latitude.
    
    Parameters:
        geographic_lat_deg : float
            Geographic latitude in degrees.
        flattening : float
            Flattening of the ellipsoid (default is WGS-84).
    
    Returns:
        geocentric_lat_deg : float
            Geocentric latitude in degrees.
    """
    phi = np.radians(geographic_lat_deg)
    geocentric_phi = np.arctan((1 - flattening) ** 2 * np.tan(phi))
    return np.degrees(geocentric_phi)

def cal_dist_az_baz(lat1:float,lon1:float,lat2:float,lon2:float):
    from obspy.geodetics import calc_vincenty_inverse

    dist_in_m,az,baz = calc_vincenty_inverse(_geod2geoc(lat1),lon1,_geod2geoc(lat2),lon2,6371000.,0.)
    
    return dist_in_m,az,baz

def cal_dist_baz(lat1:float,lon1:float,lat2:float,lon2:float):
    dist_in_m,_,baz = cal_dist_az_baz(lat1,lon1,lat2,lon2)

    return dist_in_m,baz

def rotate_EN_to_RT(ve:np.ndarray,vn:np.ndarray,bazd):
    """
    rotate NE components to RT, by using bazd
    
    Parameters
    ---------------
    ve,vn: np.ndarrray
        E/N componenets
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

def rotate_RT_to_EN(vr:np.ndarray,vt:np.ndarray,bazd):
    """
    rotate NE components to RT, by using bazd
    
    Parameters
    ---------------
    vr,vt: np.ndarrray
        R/T componenets
    bazd: float
        back azimuth, in deg

    Returns:
    ve,vn: np.ndarray
        E/N components
    """
    from numpy import sin,cos 
    baz = np.deg2rad(bazd)
    ve = -vr * sin(baz) - vt * cos(baz)
    vn = -vr * cos(baz) + vt * sin(baz)

    return ve,vn