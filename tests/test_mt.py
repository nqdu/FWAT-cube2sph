import numpy as np 
import matplotlib.pyplot as plt 
from fwat.measure.utils import bandpass,taper_window
from fwat.adjoint.multitaper import measure_adj_mt
from fwat.measure.measure import measure_adj
import time 

def main():
    nt = 1024
    t = np.linspace(0,10,nt)
    dt = t[1]-t[0]
    Tmin = 1.
    Tmax = 50.
    tstart = 3.
    tend = 8.

    lpt,rpt,taper0 = taper_window(t[0],dt,len(t),tstart,tend,p=0.05)

    # Create dispersive signals using frequency-dependent group velocity
    freq = np.fft.rfftfreq(nt, dt)
    # Simulate dispersive wave with varying arrival times by frequency
    syn_fft = np.exp(-((freq - 0.2)/0.1)**2)
    obs_fft = np.exp(-((freq - 0.2)/0.1)**2)
    
    # Apply frequency-dependent phase shifts to create dispersion
    t_center_syn = 4.5
    t_center_obs = 5.0
    phase_syn = -2 * np.pi * freq * t_center_syn - 0.5 * (2 * np.pi * freq)**2
    phase_obs = -2 * np.pi * freq * t_center_obs - 0.3 * (2 * np.pi * freq)**2
    
    syn_fft = syn_fft * np.exp(1j * phase_syn)
    obs_fft = obs_fft * np.exp(1j * phase_obs)
    
    syn = np.fft.irfft(syn_fft, nt)
    obs = np.fft.irfft(obs_fft, nt)
    
    # Add noise to make signals more complex
    noise_level = 0.0
    syn += noise_level * np.random.randn(nt)
    obs += noise_level * np.random.randn(nt)
    taper = syn * 0
    taper[lpt:rpt] = taper0
    syn *= taper
    obs *= taper
    
    # Add multiple Gaussian peaks for complexity
    #syn += 0.5 * np.exp(-((t-4.5)/0.3)**2)
    #obs += 0.5 * np.exp(-((t-5.0)/0.3)**2)


    tic = time.time()
    stats,adj = measure_adj_mt(
        obs,syn,t[0],dt,len(t),
        tstart,tend,Tmin,Tmax)
    toc = time.time()

    print("elapsed time: ",toc-tic)

    # adjoint source by measure_adj
    tic = time.time()
    stats_mt,adj_mt = measure_adj(0,dt,nt,0.,dt,nt,tstart,tend,5,
                         Tmax,Tmin,False,obs,syn)
    toc = time.time()
    print("elapsed time: ",toc-tic)

    print("misfit from user cc: ",stats.misfit)
    print("misfit from measure adj: ",stats_mt.misfit)
    # plot 
    plt.figure(1)
    plt.subplot(211)
    plt.plot(t,adj,label='adj')
    plt.plot(t,adj_mt ,label='adj-mt',ls='--')
    plt.legend()

    # syn/obs
    plt.subplot(212)
    plt.plot(t,syn,label='syn')
    plt.plot(t,obs,label='obs')
    plt.legend()

    plt.savefig("cc_phase.jpg",dpi=300)

if __name__ == "__main__":
    main()