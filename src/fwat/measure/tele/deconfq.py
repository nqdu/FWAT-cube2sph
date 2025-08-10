import numpy as np 

TSHIFT = 80. # tackle negative time

def compute_stf_freq(h,v,dt):
    n = len(h)
    n2 = n *2

    # compute fft
    h_freq = np.fft.rfft(h,n=n2)
    v_freq = np.fft.rfft(v,n=n2)
    w = 2. * np.pi * np.fft.rfftfreq(n2,dt)
    gauss = np.exp(-(w / (2 * 3.5))**2)

    # frequency domain deconvove
    vabs = np.abs(v_freq * v_freq.conj())
    max_v = np.max(vabs)
    fai = vabs * 1.
    water = 1.0e-3
    idx = vabs < water * max_v
    fai[idx] = water * max_v
    out = h_freq * v_freq.conj() / fai * gauss * np.exp(-1j * TSHIFT *w )

    return np.fft.irfft(out)[:n] / dt

def shiftconvolve(u,stf,dt):
    n = len(u)
    n2 = n * 2
    u_freq = np.fft.rfft(u,n=n2)
    stf_freq = np.fft.rfft(stf,n=n2)
    w = 2. * np.pi * np.fft.rfftfreq(n2,dt)

    # compute convolved arrays
    out = u_freq * stf_freq * np.exp(1j * TSHIFT * w) 

    return np.fft.irfft(out)[:n]  * dt