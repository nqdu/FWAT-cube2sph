try:
    import cupy as cp
except:
    import numpy as cp

# import numpy as cp
import numpy as np 


def nextpow2(n):
    i = 1
    while i < n:
        i = i * 2 
    
    return i

def gauss_filter(nt,dt,f0,device='cpu'):
    if device == 'cpu':
        freqs = np.fft.rfftfreq(nt,dt)
        out = np.exp(-0.25 * (2 * np.pi * freqs / f0) **2)
    else:
        freqs = cp.fft.rfftfreq(nt,dt)
        out = cp.exp(-0.25 * (2 * cp.pi * freqs / f0) **2)

    return out 

def apply_gaussian(data,gauss,dt,device='cpu'):

    if device == 'cpu':
        data_freq = np.fft.rfft(data) * dt 
        data_freq *= gauss 

        out = np.fft.irfft(data_freq)/ dt 

    else:
        data_freq = cp.fft.rfft(data) * dt 
        data_freq *= gauss 
        out = cp.fft.irfft(data_freq) / dt
    return out 

def mycorrelate(a,b,device='cpu'):
    if device == 'cpu':
        aft = np.fft.rfft(a)
        bft = np.fft.rfft(b)
        out = np.fft.irfft(aft * np.conjugate(bft))
        return out
    else:
        aft = cp.fft.rfft(a)
        bft = cp.fft.rfft(b)
        out = cp.fft.irfft(aft * cp.conjugate(bft))

    return out 

def myconvolve(a,b,device='cpu'):
    if device == 'cpu':
        aft = np.fft.rfft(a)
        bft = np.fft.rfft(b)
        out = np.fft.irfft(aft * bft)
    else:
        aft = cp.fft.rfft(a)
        bft = cp.fft.rfft(b)
        out = cp.fft.irfft(aft * bft)
    return out 


def deconit(u, w, dt, tshift, f0,ipart=1,maxiter=2000):
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
    nft = nextpow2(nt*2) # avoid aliasing

    # device is gpu if cupy is used
    device = 'cpu'
    if cp is not np:
        device = 'gpu'

    # allocate space for new array 
    uflt = cp.zeros((nft)); uflt[:nt] = cp.asarray(u.copy())
    wflt = cp.zeros((nft)); wflt[:nt] = cp.asarray(w.copy())
    wcopy = cp.asarray(wflt.copy())
    gauss = gauss_filter(nft,dt,f0,device=device)

    # filter 
    uflt = apply_gaussian(uflt,gauss,dt,device=device)
    wflt = apply_gaussian(wflt,gauss,dt,device=device)
    # init
    invpw = 1. / cp.sum(wflt**2) / dt 
    invpu = 1. / cp.sum(uflt**2) / dt
    P = cp.zeros((nft))
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

        cuw = mycorrelate(rflt,wflt,device=device) * dt 
        idx = cp.argmax(abs(cuw[0:max_lag]))
        P[idx] += cuw[idx] * invpw / dt 
        temp = apply_gaussian(P,gauss,dt,device=device)
        rflt = uflt - myconvolve(temp,wcopy,device=device) * dt 

        # compute error 
        sumsq = cp.sum(rflt**2) * dt * invpu
        d_error = 100 * (sumsq_i - sumsq)
        sumsq_i = sumsq

    if d_error > minderr: print(d_error)
    
    # get rf and phase shift 
    P = apply_gaussian(P,gauss,dt,device=device)
    pft = cp.fft.rfft(P) * dt
    pft *= cp.exp(-1j * cp.fft.rfftfreq(nft,dt) * cp.pi * 2 * tshift)
    rf = cp.fft.irfft(pft)[:nt]/ dt

    # convert back to numpy array if using cupy
    if cp is not np:
        rf = cp.asnumpy(rf)

    return rf 


def time_decon(u,w,dt,freqmax=None):
    # load modules
    from scipy.signal import correlate,convolve
    nt = len(u)

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
    for i in range(1500):
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