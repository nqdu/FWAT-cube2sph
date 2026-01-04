import numpy as np
import matplotlib.pyplot as plt

import yaml 
from fwat.adjoint.cross_conv import measure_adj_cross_conv
from fwat.measure.utils import bandpass,interpolate_syn
from pathlib import Path

def main():
    # load seismograms

    # get path 
    bin_path = Path(__file__).parent / "data"
    syn_r0 = np.load(bin_path / 'syn_r.npy')
    syn_z0 = np.load(bin_path / 'syn_z.npy')
    obs_r0 = np.load(bin_path / 'obs_r.npy')
    obs_z0 = np.load(bin_path / 'obs_z.npy')

    # get yaml info
    dt = 0.02
    win_tb = 5.
    win_te = 30.
    t0_fk = 35.815201
    t0 = -0.200000003
    freqmin = 0.02
    freqmax = 0.4
    tstart = -win_tb + t0_fk
    tend = win_te + t0_fk

    t = np.arange(len(syn_r0)) * dt + t0
    nt = len(t)
    dt1 = 0.1 # original dt for finite-difference seismograms
    nt1 = int((nt -1) * dt / dt1) + 1
    t1 = np.arange(nt1) * dt1 + t0

    # reset freqmax if beyond Nyquist
    nyq = 0.5 / dt1
    if freqmax > nyq:
        freqmax = nyq - 0.01
        print(f"Reset freqmax to {freqmax} Hz due to Nyquist limit.")

    # band pass 
    syn_z0 = interpolate_syn(syn_z0, t0, dt, nt, t0, dt1, nt1)
    syn_r0 = interpolate_syn(syn_r0, t0, dt, nt, t0, dt1, nt1)
    obs_z0 = interpolate_syn(obs_z0, t0, dt, nt, t0, dt1, nt1)
    obs_r0 = interpolate_syn(obs_r0, t0, dt, nt, t0, dt1, nt1)
    syn_z = bandpass(syn_z0, dt1, freqmin, freqmax)
    syn_r = bandpass(syn_r0, dt1, freqmin, freqmax)
    obs_z = bandpass(obs_z0, dt1, freqmin, freqmax)
    obs_r = bandpass(obs_r0, dt1, freqmin, freqmax)

    # add min_period to tend to avoid edge effect
    tend += 1. / freqmax

    # compute adjoint source
    chi,_,_,adj_z,adj_r,cc1,cc2 =  \
        measure_adj_cross_conv(
            obs_z,syn_z,
            obs_r,syn_r,
            t0,dt1,
            tstart,tend
        )
    adj_r = bandpass(adj_r,dt1,freqmin,freqmax)
    adj_z = bandpass(adj_z,dt1,freqmin,freqmax)

    # compute FD adjoint source for validation
    eps = 1e-6
    adj_r_fd = np.zeros_like(adj_r)
    adj_z_fd = np.zeros_like(adj_z)
    print("Computing finite-difference adjoint sources Z...")
    for i in range(len(syn_z)):
        syn_z_p = syn_z0.copy()
        syn_z_m = syn_z0.copy()
        syn_z_p[i] += eps
        syn_z_m[i] -= eps

        # bandpass amd interpolate
        syn_z_p = bandpass(syn_z_p,dt1, freqmin, freqmax)
        syn_z_m = bandpass(syn_z_m,dt1, freqmin, freqmax)
        chi_p,_,_,_,_,_,_ =  \
            measure_adj_cross_conv(
                obs_z,syn_z_p,
                obs_r,syn_r,
                t0,dt1,
                tstart,tend
            )
        chi_m,_,_,_,_,_,_ =  \
            measure_adj_cross_conv(
                obs_z,syn_z_m,
                obs_r,syn_r,
                t0,dt1,
                tstart,tend
            )
        adj_z_fd[i] = (chi_p - chi_m) / (2 * eps)

    print("Computing finite-difference adjoint sources R...")
    for i in range(len(syn_r)):
        syn_r_p = syn_r0.copy()
        syn_r_m = syn_r0.copy()
        syn_r_p[i] += eps
        syn_r_m[i] -= eps

        # bandpass
        syn_r_p = bandpass(syn_r_p,dt1, freqmin, freqmax)
        syn_r_m = bandpass(syn_r_m,dt1, freqmin, freqmax)
        chi_p,_,_,_,_,_,_ =  \
            measure_adj_cross_conv(
                obs_z,syn_z,
                obs_r,syn_r_p,
                t0,dt1,
                1./freqmax,1./freqmin,
            )
        chi_m,_,_,_,_,_,_ =  \
            measure_adj_cross_conv(
                obs_z,syn_z,
                obs_r,syn_r_m,
                t0,dt1,
                tstart,tend
            )
        adj_r_fd[i] = (chi_p - chi_m) / (2 * eps)  

    # bandpass adjoint sources for better visualization
    adj_z_fd = bandpass(adj_z_fd,dt1,freqmin,freqmax) / dt1
    adj_r_fd = bandpass(adj_r_fd,dt1,freqmin,freqmax) / dt1

    plt.figure(figsize=(12, 10))
    
    # Plot Z
    plt.subplot(211)
    plt.plot(t1, adj_z, 'k', linewidth=3, label='Analytical (Code)')
    plt.plot(t1, adj_z_fd, 'r--', linewidth=2, label='Finite Difference (Truth)')
    plt.title(f"Vertical Adjoint Source)")
    plt.ylabel("Amplitude")
    plt.legend()
    plt.grid(True, alpha=0.3)

    # Plot R
    plt.subplot(212)
    plt.plot(t1, adj_r, 'k', linewidth=3, label='Analytical (Code)')
    plt.plot(t1, adj_r_fd, 'r--', linewidth=2, label='Finite Difference (Truth)')
    plt.title(f"Radial Adjoint Source (Max Error: %)")
    plt.xlabel("Time (s)")
    plt.ylabel("Amplitude")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig("cc-adj.png", dpi=150)

    # plot vsyn vobs hsyn hobs and cc1 
    plt.figure(figsize=(12, 8))
    plt.subplot(311)
    plt.plot(t1, obs_z, 'b-', linewidth=2, label='Vobs')
    plt.plot(t1, syn_z, 'r--', linewidth=2, label='Vsyn')
    plt.axvline(x=tstart, color='k', linestyle='--', linewidth=1)
    plt.axvline(x=tend, color='k', linestyle='--', linewidth=1)
    plt.legend()

    plt.subplot(312)
    plt .plot(t1, syn_r, 'g--', linewidth=2, label='HSyn')
    plt.plot(t1, obs_r, 'g-', linewidth=2, label='Hobs')
    plt.legend()

    # cc1   and cc2 for visualization
    plt.subplot(313)
    t2 = t1[0] + np.arange(len(cc1)) * dt
    plt.plot(t2, cc1, 'b-', linewidth=2, label='Cross-Conv Vobs * Hsyn')
    plt.plot(t2, cc2, 'g--', linewidth=2, label='Cross-Conv Vsyn * Hobs')
    plt.title("Cross-Convolutions for Visualization")
    plt.xlabel("Time (s)")
    plt.ylabel("Amplitude")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig("cc-vs.png", dpi=150)

if __name__ == "__main__":
    main()