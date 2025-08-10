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

def interpolate_syn(data,t1,dt1,npt1,t2,dt2,npt2,max_percentage=0.05):
    """
    interpolate data from (t1, dt1, npt1) to a new data (t2,dt2,npt2)

    data: np.ndarray
        input data, shape(npt1) with starttime t1 and interval dt1 
    max_percentage: float
        max_percentage of taper on one side, if required, taper * len(data) * 2 is the window used 
    
    """
    # taper input data if required
    data1 = np.float32(data)
    if max_percentage > 0.:
        cos_tp = sac_cos_taper(npt1,max_percentage*2)
        data1 = data1 * cos_tp

    temp = np.zeros((npt2))
    time = t2 + np.arange(npt2) * dt2 
    idx = np.logical_and(time > t1,time < t1 + (npt1-1) * dt1)
    ii = np.int64((time[idx] - t1) / dt1)
    tt = time[idx] - (ii * dt1 + t1) 
    temp[idx] = (data1[ii+1] - data1[ii]) * tt / dt1 + data1[ii]
    
    return temp


def bandpass(u,dt,freqmin,freqmax):
    import obspy
    tr = obspy.Trace(data = u * 1.)
    tr.stats.delta = dt 

    tr.detrend("demean")
    tr.detrend("linear")
    tr.taper(0.05)
    tr.filter("bandpass",freqmin=freqmin,freqmax=freqmax,zerophase=True,corners=4)
    tr.detrend("demean")
    tr.detrend("linear")
    tr.taper(0.05)

    return tr.data

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

def taper_window(t0,dt,tstart,tend,p=0.05):
    """
    get window start/end sample index

    Parameters
    -----------
    t0,dt: float
        start time and sampling interval
    tstart,tend: float
        start/end time of this window
    p: float
        ratio of data to be tapered on one side

    Returns
    -----------
    left_pt,right_pt: int
        the window sample index [left_pt,right_pt)
    """
    assert tend >= tstart, "window is reversed!"

    nlen = int((tend - tstart) / dt) + 1
    left_pt = int(np.floor(tstart - t0) / dt)
    right_pt = left_pt + nlen

    win = sac_cos_taper(nlen,p*2)

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