import numpy as np 

def nextpow2(n):
    i = 1
    while i < n:
        i = i * 2 
    
    return i

def gauss_filter(nt,dt,f0):
    freqs = np.fft.rfftfreq(nt,dt)
    out = np.exp(-0.25 * (2 * np.pi * freqs / f0) **2)

    return out 

def lowpass_filter(d,dt,freqmax):
    from scipy.signal import sosfiltfilt,butter 

    # get sos coefs
    sos = butter(4,freqmax * dt * 2,output='sos')
    out = sosfiltfilt(sos,d)

    return out

def apply_gaussian(data,gauss,dt):
    data_freq = np.fft.rfft(data) * dt 
    data_freq *= gauss 

    out = np.fft.irfft(data_freq) / dt 

    return out 

def mycorrelate(a,b):
    aft = np.fft.rfft(a)
    bft = np.fft.rfft(b)
    out = np.fft.irfft(aft * np.conjugate(bft))

    return out 

def myconvolve(a,b):
    aft = np.fft.rfft(a)
    bft = np.fft.rfft(b)
    out = np.fft.irfft(aft * bft)

    return out 

def deconit(u, w, dt, tshift, f0,ipart=1,maxiter=120,):
    """
    Time iterative deconvolution from cps330

    Parameters
    -----------
    u: np.ndarray
        numerator, shape(npts)
    w: np.ndarray
        denominator, shape(npts)
    tshift: float
        time shift of deconvoluted arrays
    f0: float
        gauss filter parameter
    ipart: int
        =1 full arrays
        =0 positive half
    maxiter: int
        max iteration
    
    Returns
    ----------
    rf: np.ndarray
        deconvoluted array, shape(npts)
    """
    nt = len(u)
    nft = nextpow2(nt * 2) # avoid aliasing

    # allocate space for new array 
    uflt = np.zeros((nft)); uflt[:nt] = u.copy()
    wflt = np.zeros((nft)); wflt[:nt] = w.copy()
    wcopy = wflt.copy()
    gauss = gauss_filter(nft,dt,f0)

    # filter 
    uflt = apply_gaussian(uflt,gauss,dt)
    wflt = apply_gaussian(wflt,gauss,dt)

    # init
    invpw = 1. / np.sum(wflt**2) / dt 
    invpu = 1. / np.sum(uflt**2) / dt
    P = np.zeros((nft))
    sumsq_i = 1.0
    minderr = 0.001
    d_error = 100 * invpw + minderr
    rflt = uflt.copy()

    # max lag
    max_lag = nft 
    if ipart == 0:
        max_lag = int(0.5 * nft)

    for i in range(maxiter):
        if abs(d_error) <= minderr : break; 

        cuw = mycorrelate(rflt,wflt) * dt 
        idx = np.argmax(abs(cuw[0:max_lag]))
        P[idx] += cuw[idx] * invpw / dt 
        temp = apply_gaussian(P,gauss,dt)
        rflt = uflt - myconvolve(temp,wcopy) * dt 

        # compute error 
        sumsq = np.sum(rflt**2) * dt * invpu
        d_error = 100 * (sumsq_i - sumsq)
        sumsq_i = sumsq
    
    # get rf and phase shift 
    P = apply_gaussian(P,gauss,dt)
    pft = np.fft.rfft(P) * dt 
    pft *= np.exp(-1j * np.fft.rfftfreq(nft,dt) * np.pi * 2 * tshift)
    rf = np.fft.irfft(pft)[:nt] / dt

    return rf 

def _time_decon_with_filt(u,w,dt,freqmax):
    from scipy.signal import fftconvolve
    nt = len(u)
    nt2 = nextpow2(nt)

    # allocate space for new array 
    uflt = np.zeros((nt2))
    wflt = np.zeros((nt2))
    uflt[:nt] = u.copy()
    wflt[:nt] = w.copy()
    wcopy = wflt.copy()

    # filter
    # wflt = lowpass_filter(wflt,dt,freqmax)
    # uflt = lowpass_filter(uflt,dt,freqmax)
    P = np.zeros((nt2))

    invpow = 1. / (dt * np.sum(wflt**2))
    invpowu = 1. / (dt * np.sum(uflt**2))
    sumsq_i = 1.
    sumsq = 50.
    minderr = 0.001 
    d_error = 100 * invpow + minderr

    rflt = uflt.copy()
    for i in range(120):
        if abs(d_error) < minderr or d_error < 0: break

        # find max
        cuw = fftconvolve(rflt,wflt[::-1],'full') * dt 
        idx = np.argmax(abs(cuw[nt2:]))
        amp = cuw[idx + nt2] * invpow  / dt
        P[idx] += amp
        temp = lowpass_filter(P,dt,freqmax)
        rflt = uflt - fftconvolve(temp,wcopy,'full')[:nt2] * dt 
        
        sumsq = np.sum(rflt**2) * dt * invpowu 
        d_error = 100 * (sumsq_i - sumsq)
        sumsq_i = sumsq

        print(i,d_error)

    rf = lowpass_filter(P,dt,freqmax)
    return rf[:nt] 


def time_decon(u,w,dt,freqmax=None):
    # load modules
    from scipy.signal import correlate,convolve
    nt = len(u)

    if freqmax != None:
        return _time_decon_with_filt(u,w,dt,freqmax)

    # allocate space for new array 
    uflt = u.copy()
    wflt = w.copy()
    rf = np.zeros((nt),'f4')
    invpow = 1. / (dt * np.sum(wflt**2))
    invpowu = 1. / (dt * np.sum(uflt**2))
    sumsq_i = 1.
    sumsq = 50.
    minderr = 0.001 
    d_error = 100 * invpow + minderr

    dobs = uflt * 1.
    src_one = np.zeros((nt))
    for i in range(200):
        if abs(d_error) < minderr: break

        #cuw = mycorrelate(dobs,wflt) * dt 
        cuw = correlate(dobs,wflt,'same') * dt 
        
        # fetch center
        n = len(cuw)
        n1 = int(10 / dt)
        idx = np.argmax(abs(cuw[:]))
        amp = cuw[idx] * invpow  / dt
        rf[idx] += amp
        src_one[idx] = amp 
        #dobs -= myconvolve(src_one,wflt) * dt 
        dobs -= convolve(src_one,wflt,'same') * dt 

        sumsq = np.sum(dobs**2) * dt * invpowu 
        d_error = 100 * (sumsq_i - sumsq)
        sumsq_i = sumsq

        src_one[idx] = 0.

    return rf 