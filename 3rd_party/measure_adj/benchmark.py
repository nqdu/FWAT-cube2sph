import glob
import os
import sys
import time

import matplotlib.pyplot as plt
import numpy as np

from fwat.measure.measure import measure_adj


def _load_libmeas():
    candidates = glob.glob(os.path.join("lib/", "libmeas*.so"))
    if not candidates:
        raise RuntimeError("Could not find libmeas module in lib/")
    sys.path.insert(0, os.path.abspath("lib"))
    import libmeas  # type: ignore

    return libmeas


def make_gaussian_signals(npts: int, dt: float, shift: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    t0 = 0.0
    t = t0 + np.arange(npts) * dt
    t_center = 20.0
    sigma = 2.0
    obs = np.exp(-0.5 * ((t - t_center) / sigma) ** 2)
    syn = np.exp(-0.5 * ((t - (t_center + shift)) / sigma) ** 2)
    return t, obs, syn


def make_complex_signals(npts: int, dt: float, shift: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rng = np.random.default_rng(12345)
    t0 = 0.0
    t = t0 + np.arange(npts) * dt
    centers = np.array([12.0, 25.0, 42.0])
    sigmas = np.array([1.5, 3.0, 2.2])
    amps = np.array([1.0, 0.6, 0.8])

    obs = np.zeros_like(t)
    syn = np.zeros_like(t)
    for c, s, a in zip(centers, sigmas, amps):
        obs += a * np.exp(-0.5 * ((t - c) / s) ** 2)
        syn += (a * 0.95) * np.exp(-0.5 * ((t - (c + shift)) / (s * 1.05)) ** 2)

    f0 = 0.03
    f1 = 0.12
    k = (f1 - f0) / (t[-1] - t[0])
    chirp = np.sin(2.0 * np.pi * (f0 * t + 0.5 * k * t * t))
    obs += 0.2 * chirp * np.exp(-0.5 * ((t - 30.0) / 12.0) ** 2)
    syn += 0.18 * np.sin(2.0 * np.pi * (f0 * (t - shift) + 0.5 * k * (t - shift) ** 2)) * np.exp(
        -0.5 * ((t - 30.0) / 12.0) ** 2
    )

    obs += 0.01 * rng.standard_normal(size=t.shape)
    syn += 0.01 * rng.standard_normal(size=t.shape)

    obs *= np.hanning(len(obs))
    syn *= np.hanning(len(syn))
    return t, obs, syn


def run_measurements():
    libmeas = _load_libmeas()

    npts = 1027
    dt = 0.1
    tt = 0.0
    dtt = dt
    nn = npts
    t0 = 0.0
    tstart = 5.0
    tend = 60.0
    tlong = 60.0
    tshort = 5.0
    shift = 2
    t, obs, syn = make_complex_signals(npts, dt, shift)
    obs *= 0
    syn *= 0
    #t, obs, syn = make_gaussian_signals(npts, dt, shift)


    base_kwargs = dict(
        t0=t0,
        dt=dt,
        npts=npts,
        tt=tt,
        dtt=dtt,
        nn=nn,
        tstart=tstart,
        tend=tend,
        tlong=tlong,
        tshort=tshort,
        verbose=False,
        RUN_BANDPASS=False,
        DISPLAY_DETAILS=False,
        OUTPUT_MEASUREMENT_FILES=False,
        COMPUTE_ADJOINT_SOURCE=True,
        TSHIFT_MIN=-4.5,
        TSHIFT_MAX=4.5,
        DLNA_MIN=-1.5,
        DLNA_MAX=1.5,
        CC_MIN=0.1,
        ERROR_TYPE=1,
        DT_SIGMA_MIN=1,
        DLNA_SIGMA_MIN=0.5,
        ITAPER=1,
        WTR=0.02,
        NPI=2.5,
        DT_FAC=2,
        ERR_FAC=2.5,
        DT_MAX_SCALE=3.5,
        NCYCLE_IN_WINDOW=1.5,
        USE_PHYSICAL_DISPERSION=False,
    )

    out_dir = "plots"
    os.makedirs(out_dir, exist_ok=True)

    timing_rows = []
    for imeas in range(1,9):
        if imeas <=6 and imeas >=3:
            itaper = 2
        else:
            itaper = 1
        base_kwargs["ITAPER"] = itaper
        t_start = time.perf_counter()
        tr_chi_f, am_chi_f, window_chi_f, adj_f = libmeas.measure(
            t0,
            dt,
            npts,
            tt,
            dtt,
            nn,
            tstart,
            tend,
            imeas,
            tlong,
            tshort,
            False,
            obs,
            syn,
            False,
            False,
            False,
            True,
            base_kwargs["TSHIFT_MIN"],
            base_kwargs["TSHIFT_MAX"],
            base_kwargs["DLNA_MIN"],
            base_kwargs["DLNA_MAX"],
            base_kwargs["CC_MIN"],
            base_kwargs["ERROR_TYPE"],
            base_kwargs["DT_SIGMA_MIN"],
            base_kwargs["DLNA_SIGMA_MIN"],
            base_kwargs["ITAPER"],
            base_kwargs["WTR"],
            base_kwargs["NPI"],
            base_kwargs["DT_FAC"],
            base_kwargs["ERR_FAC"],
            base_kwargs["DT_MAX_SCALE"],
            base_kwargs["NCYCLE_IN_WINDOW"],
            False,
        )
        t_fortran = time.perf_counter() - t_start

        t_start = time.perf_counter()
        stats_p,adj_out = \
            measure_adj(  \
                t0,dt,npts,tt,dtt,nn,
                tstart,tend,imeas,tlong,tshort,
                verbose=False,obs_data=obs,syn_data=syn,
                itaper = itaper,tshift_min=base_kwargs["TSHIFT_MIN"],
                tshift_max=base_kwargs["TSHIFT_MAX"],
                dlna_min=base_kwargs["DLNA_MIN"],
                dlna_max=base_kwargs["DLNA_MAX"],
                cc_min=base_kwargs["CC_MIN"],
                err_type=int(base_kwargs["ERROR_TYPE"]),
                dt_sigma_min=base_kwargs["DT_SIGMA_MIN"],
                dlna_sigma_min=base_kwargs["DLNA_SIGMA_MIN"],
                wtr=base_kwargs["WTR"],
                npi=base_kwargs["NPI"],
                dt_fac=base_kwargs["DT_FAC"],
                err_fac=base_kwargs["ERR_FAC"],
                dt_max_scale=base_kwargs["DT_MAX_SCALE"],
                ncyle_in_window=base_kwargs["NCYCLE_IN_WINDOW"],
                use_physical_disp=bool(base_kwargs["USE_PHYSICAL_DISPERSION"]),
            )
        
        mtm_used = imeas in (7, 8)
        t_python = time.perf_counter() - t_start
    

        adj_f = np.array(adj_f)
        adj_out = np.array(adj_out)
        tr_chi_p = stats_p.tr_chi
        am_chi_p = stats_p.am_chi

        adj_diff = adj_out - adj_f
        adj_ref = np.max(np.abs(adj_f)) + 1e-12
        adj_rel = np.max(np.abs(adj_diff)) / adj_ref
        adj_rms = np.sqrt(np.mean(adj_diff**2)) / (np.sqrt(np.mean(adj_f**2)) + 1e-12)

        print(f"imeas={imeas}")
        print(f"  tr_chi fortran/python: {tr_chi_f} / {tr_chi_p}")
        print(f"  am_chi fortran/python: {am_chi_f} / {am_chi_p}")
        if imeas in (7, 8):
            print(f"  mtm used: {mtm_used}")
        print(f"  adj max abs diff: {np.max(np.abs(adj_diff))}")
        print(f"  fortran/python time (s): {t_fortran} / {t_python}")
        timing_rows.append((imeas, t_fortran, t_python))

        fig, axes = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
        axes[0].plot(t, adj_f, label="Fortran", linewidth=1.2)
        axes[0].plot(t, adj_out, label="Python", linewidth=1.0, linestyle="--")
        axes[0].set_title(f"Adjoint Source Comparison (imeas={imeas})")
        axes[0].legend()
        axes[0].grid(True, alpha=0.3)

        axes[1].plot(t, adj_diff, label="Python - Fortran", color="tab:red", linewidth=1.0)
        axes[1].set_xlabel("Time (s)")
        axes[1].set_ylabel("Difference")
        axes[1].grid(True, alpha=0.3)

        fig.tight_layout()
        fig.savefig(os.path.join(out_dir, f"adjoint_compare_imeas{imeas}.png"), dpi=150)
        plt.close(fig)

        if len(timing_rows) == 1:  # Only plot obs/syn for the first imeas to avoid too many plots
            fig = plt.figure(figsize=(10, 4))
            plt.plot(t,obs, label="Observed", color="black", linewidth=1.0)
            plt.plot(t,syn, label="Synthetic", color="tab:blue", linewidth=1.0, linestyle="--")
            plt.xlabel("Time (s)")
            plt.tight_layout()
            plt.legend()
            plt.savefig(os.path.join(out_dir, f"obs_syn_data.png"), dpi=150)

    print("\nTime Comparison (seconds)")
    print("imeas | fortran | python | speedup")
    print("----- | ------- | ------ | -------")
    for imeas, t_fortran, t_python in timing_rows:
        speedup = t_fortran / t_python if t_python > 0 else float("inf")
        print(f"{imeas:5d} | {t_fortran:7.6f} | {t_python:6.6f} | {speedup:7.2f}x")


if __name__ == "__main__":
    run_measurements()
