import numpy as np 

def read_fwat_params(paramfile="fwat_params/FWAT.PAR.yaml"):
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

def get_average_amplitude(glob_obs,icomp):
    """
    get the average amplitude for icomp-th component 

    Parameters:
        glob_obs: np.ndarray
            observation data in global, shape(nsta,ncomp,nt)
        icomp: int
            icomp-th componnet
    """
    nsta = glob_obs.shape[0]
    avgamp0 = 0.
    for i in range(nsta):
        avgamp0 += np.max(np.abs(glob_obs[i,0,:]))
    avgamp0 /= nsta 
    igood =0 
    avgamp = 0. 
    for i in range(nsta):
        if np.max(np.abs(glob_obs[i,0,:])) - avgamp0 < 0.2 * avgamp0:
            igood += 1 
            avgamp += np.max(np.abs(glob_obs[i,0,:])) 
    avgamp /= igood 

    return avgamp

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
