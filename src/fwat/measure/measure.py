from dataclasses import dataclass
from functools import lru_cache
from typing import Optional, Tuple, Final

import numpy as np
from scipy.signal.windows import dpss

from fwat.adjoint.MeasureStats import MeasureStats
from fwat.measure.utils import bandpass as bandpass_fwat
 

PI = np.pi
TWOPI = 2.0 * PI
CCI = 1j
LARGE_VAL: Final[float] = 1.0e8
PHASE_STEP: Final[float] = 1.5 * PI
TRBDNDW: Final[float] = 0.3
APARM: Final[float] = 30.0
IORD: Final[int] = 4
PASSES: Final[int] = 2
def next2pow2(n: int) -> int:
    if n <= 1:
        return 1
    return 1 << (n - 1).bit_length()


def select_nfft(n: int, lpt: int) -> int:
    return max(1 << lpt, next2pow2(n))


@dataclass
class MAParams:
    out_dir: str = "OUTPUT_FILES"
    tlong: float = 0.0
    tshort: float = 0.0
    wtr: float = 0.0
    npi: float = 0.0
    dt_fac: float = 0.0
    err_fac: float = 0.0
    dt_max_scale: float = 0.0
    ncycle_in_window: float = 0.0
    tshift_min: float = -LARGE_VAL
    tshift_max: float = LARGE_VAL
    dlna_min: float = -LARGE_VAL
    dlna_max: float = LARGE_VAL
    cc_min: float = 0.0
    dt_sigma_min: float = 0.0
    dlna_sigma_min: float = 0.0
    ntaper: int = 0
    ipwr_t: int = 10
    ipwr_w: int = 10
    error_type: int = 0
    imeas0: int = 1
    imeas: int = 1
    itaper: int = 1
    is_mtm0: int = 0
    is_mtm: int = 0
    display_details: bool = False
    output_measurement_files: bool = False
    run_bandpass: bool = False
    compute_adjoint_source: bool = True
    use_physical_dispersion: bool = False
    lpt: int = 15




def _rfft(x: np.ndarray, nfft: int, dt: float) -> np.ndarray:
    return np.fft.rfft(x, nfft) * dt


def _irfft(x: np.ndarray, nfft: int, dt: float) -> np.ndarray:
    return np.fft.irfft(x, nfft) / dt


def bandpass(x: np.ndarray, n: int, delta_t: float, f1: float, f2: float) -> None:
    if n <= 1:
        return
    # sos = butter(N=4, Wn=[f1, f2], btype='bandpass', fs=1.0/delta_t, output='sos')
    # x[:n] = sosfiltfilt(sos, x[:n])
    x[:n] = bandpass_fwat(x[:n], delta_t, f1, f2)


def xapiir_sub(
    data: np.ndarray,
    aproto: str,
    trbndw: float,
    a: float,
    iord: int,
    ftype: str,
    flo: float,
    fhi: float,
    ts: float,
    passes: int,
) -> np.ndarray:
    sn, sd, nsects = _design(iord, ftype, aproto, a, trbndw, flo, fhi, ts)
    zp = passes != 1
    return _apply_iir(data.astype(float, copy=True), zp, sn, sd, nsects)


def _apply_iir(data: np.ndarray, zp: bool, sn: np.ndarray, sd: np.ndarray, nsects: int) -> np.ndarray:
    if nsects <= 0:
        return data
    jptr = 0
    for _ in range(nsects):
        x1 = 0.0
        x2 = 0.0
        y1 = 0.0
        y2 = 0.0
        b0, b1, b2 = sn[jptr : jptr + 3]
        a1, a2 = sd[jptr + 1 : jptr + 3]
        for i in range(len(data)):
            out = b0 * data[i] + b1 * x1 + b2 * x2
            out = out - (a1 * y1 + a2 * y2)
            y2 = y1
            y1 = out
            x2 = x1
            x1 = data[i]
            data[i] = out
        jptr += 3

    if not zp:
        return data

    jptr = 0
    for _ in range(nsects):
        x1 = 0.0
        x2 = 0.0
        y1 = 0.0
        y2 = 0.0
        b0, b1, b2 = sn[jptr : jptr + 3]
        a1, a2 = sd[jptr + 1 : jptr + 3]
        for i in range(len(data) - 1, -1, -1):
            out = b0 * data[i] + b1 * x1 + b2 * x2
            out = out - (a1 * y1 + a2 * y2)
            y2 = y1
            y1 = out
            x2 = x1
            x1 = data[i]
            data[i] = out
        jptr += 3
    return data


def _design(
    iord: int,
    ftype: str,
    aproto: str,
    a: float,
    trbndw: float,
    fl: float,
    fh: float,
    ts: float,
) -> tuple[np.ndarray, np.ndarray, int]:
    rtype = [""] * 10
    p = np.zeros(10, dtype=np.complex128)
    z = np.zeros(10, dtype=np.complex128)
    nsects = 0
    if aproto == "BU":
        p, rtype, dcvalue, nsects = _buroots(iord)
    else:
        raise ValueError("Only BU prototype supported")

    sn = np.zeros(3 * 10)
    sd = np.zeros(3 * 10)

    if ftype == "BP":
        flw = _warp(fl * ts / 2.0, 2.0)
        fhw = _warp(fh * ts / 2.0, 2.0)
        sn, sd, nsects = _lptbp(p, z, rtype, dcvalue, nsects, flw, fhw)
    elif ftype == "LP":
        fhw = _warp(fh * ts / 2.0, 2.0)
        sn, sd = _lp(p, z, rtype, dcvalue, nsects)
        sn, sd = _cutoffs(sn, sd, nsects, fhw)
    elif ftype == "HP":
        flw = _warp(fl * ts / 2.0, 2.0)
        sn, sd = _lpthp(p, z, rtype, dcvalue, nsects)
        sn, sd = _cutoffs(sn, sd, nsects, flw)
    else:
        raise ValueError("Only BP/LP/HP supported")

    sn, sd = _bilin2(sn, sd, nsects)
    return sn, sd, nsects


def _buroots(iord: int) -> tuple[np.ndarray, list[str], float, int]:
    p = np.zeros(10, dtype=np.complex128)
    rtype = [""] * 10
    half = iord // 2
    nsects = 0
    if 2 * half < iord:
        p[0] = -1.0 + 0.0j
        rtype[0] = "SP"
        nsects = 1
    for k in range(1, half + 1):
        angle = PI * (0.5 + (2 * k - 1) / (2 * iord))
        nsects += 1
        p[nsects - 1] = complex(np.cos(angle), np.sin(angle))
        rtype[nsects - 1] = "CP"
    dcvalue = 1.0
    return p, rtype, dcvalue, nsects


def _warp(f: float, ts: float) -> float:
    twopi = 2.0 * PI
    angle = twopi * f * ts / 2.0
    return (2.0 * np.tan(angle) / ts) / twopi


def _lp(p: np.ndarray, z: np.ndarray, rtype: list[str], dcvalue: float, nsects: int) -> tuple[np.ndarray, np.ndarray]:
    sn = np.zeros(3 * 10)
    sd = np.zeros(3 * 10)
    iptr = 0
    for i in range(nsects):
        rt = rtype[i]
        if rt == "CPZ":
            scale = (p[i] * np.conj(p[i])).real / (z[i] * np.conj(z[i])).real
            sn[iptr] = (z[i] * np.conj(z[i])).real * scale
            sn[iptr + 1] = -2.0 * z[i].real * scale
            sn[iptr + 2] = 1.0 * scale
            sd[iptr] = (p[i] * np.conj(p[i])).real
            sd[iptr + 1] = -2.0 * p[i].real
            sd[iptr + 2] = 1.0
        elif rt == "CP":
            scale = (p[i] * np.conj(p[i])).real
            sn[iptr] = scale
            sn[iptr + 1] = 0.0
            sn[iptr + 2] = 0.0
            sd[iptr] = (p[i] * np.conj(p[i])).real
            sd[iptr + 1] = -2.0 * p[i].real
            sd[iptr + 2] = 1.0
        elif rt == "SP":
            scale = -p[i].real
            sn[iptr] = scale
            sn[iptr + 1] = 0.0
            sn[iptr + 2] = 0.0
            sd[iptr] = -p[i].real
            sd[iptr + 1] = 1.0
            sd[iptr + 2] = 0.0
        iptr += 3
    sn[0:3] = dcvalue * sn[0:3]
    return sn, sd


def _lptbp(
    p: np.ndarray,
    z: np.ndarray,
    rtype: list[str],
    dcvalue: float,
    nsects: int,
    fl: float,
    fh: float,
) -> tuple[np.ndarray, np.ndarray, int]:
    sn = np.zeros(3 * 10)
    sd = np.zeros(3 * 10)
    pi = PI
    twopi = 2.0 * pi
    aa = twopi * twopi * fl * fh
    bb = twopi * (fh - fl)

    n = nsects
    nsects = 0
    iptr = 0
    for i in range(n):
        rt = rtype[i]
        if rt == "CPZ":
            ctmp = (bb * z[i]) ** 2 - 4.0 * aa
            ctmp = complex(ctmp) ** 0.5
            z1 = 0.5 * (bb * z[i] + ctmp)
            z2 = 0.5 * (bb * z[i] - ctmp)

            ctmp = (bb * p[i]) ** 2 - 4.0 * aa
            ctmp = complex(ctmp) ** 0.5
            p1 = 0.5 * (bb * p[i] + ctmp)
            p2 = 0.5 * (bb * p[i] - ctmp)

            sn[iptr] = (z1 * np.conj(z1)).real
            sn[iptr + 1] = -2.0 * z1.real
            sn[iptr + 2] = 1.0
            sd[iptr] = (p1 * np.conj(p1)).real
            sd[iptr + 1] = -2.0 * p1.real
            sd[iptr + 2] = 1.0
            iptr += 3

            sn[iptr] = (z2 * np.conj(z2)).real
            sn[iptr + 1] = -2.0 * z2.real
            sn[iptr + 2] = 1.0
            sd[iptr] = (p2 * np.conj(p2)).real
            sd[iptr + 1] = -2.0 * p2.real
            sd[iptr + 2] = 1.0
            iptr += 3
            nsects += 2
        elif rt == "CP":
            ctmp = (bb * p[i]) ** 2 - 4.0 * aa
            ctmp = complex(ctmp) ** 0.5
            p1 = 0.5 * (bb * p[i] + ctmp)
            p2 = 0.5 * (bb * p[i] - ctmp)

            sn[iptr] = 0.0
            sn[iptr + 1] = bb
            sn[iptr + 2] = 0.0
            sd[iptr] = (p1 * np.conj(p1)).real
            sd[iptr + 1] = -2.0 * p1.real
            sd[iptr + 2] = 1.0
            iptr += 3

            sn[iptr] = 0.0
            sn[iptr + 1] = bb
            sn[iptr + 2] = 0.0
            sd[iptr] = (p2 * np.conj(p2)).real
            sd[iptr + 1] = -2.0 * p2.real
            sd[iptr + 2] = 1.0
            iptr += 3
            nsects += 2
        elif rt == "SP":
            sn[iptr] = 0.0
            sn[iptr + 1] = bb
            sn[iptr + 2] = 0.0
            sd[iptr] = aa
            sd[iptr + 1] = -bb * p[i].real
            sd[iptr + 2] = 1.0
            iptr += 3
            nsects += 1

    s = complex(0.0, np.sqrt(aa))
    h = complex(1.0, 0.0)
    iptr = 0
    for _ in range(nsects):
        num = (sn[iptr + 2] * s + sn[iptr + 1]) * s + sn[iptr]
        den = (sd[iptr + 2] * s + sd[iptr + 1]) * s + sd[iptr]
        h *= num / den
        iptr += 3
    scale = dcvalue / np.sqrt((h * np.conj(h)).real)
    sn[0:3] = sn[0:3] * scale
    return sn, sd, nsects


def _lpthp(
    p: np.ndarray,
    z: np.ndarray,
    rtype: list[str],
    dcvalue: float,
    nsects: int,
) -> tuple[np.ndarray, np.ndarray]:
    sn = np.zeros(3 * 10)
    sd = np.zeros(3 * 10)
    iptr = 0
    for i in range(nsects):
        rt = rtype[i]
        if rt == "CPZ":
            scale = (p[i] * np.conj(p[i])).real / (z[i] * np.conj(z[i])).real
            sn[iptr] = 1.0 * scale
            sn[iptr + 1] = -2.0 * z[i].real * scale
            sn[iptr + 2] = (z[i] * np.conj(z[i])).real * scale
            sd[iptr] = 1.0
            sd[iptr + 1] = -2.0 * p[i].real
            sd[iptr + 2] = (p[i] * np.conj(p[i])).real
        elif rt == "CP":
            scale = (p[i] * np.conj(p[i])).real
            sn[iptr] = 0.0
            sn[iptr + 1] = 0.0
            sn[iptr + 2] = scale
            sd[iptr] = 1.0
            sd[iptr + 1] = -2.0 * p[i].real
            sd[iptr + 2] = (p[i] * np.conj(p[i])).real
        elif rt == "SP":
            scale = -p[i].real
            sn[iptr] = 0.0
            sn[iptr + 1] = scale
            sn[iptr + 2] = 0.0
            sd[iptr] = 1.0
            sd[iptr + 1] = -p[i].real
            sd[iptr + 2] = 0.0
        iptr += 3
    sn[0:3] = sn[0:3] * dcvalue
    return sn, sd


def _cutoffs(sn: np.ndarray, sd: np.ndarray, nsects: int, f: float) -> tuple[np.ndarray, np.ndarray]:
    scale = 2.0 * PI * f
    iptr = 0
    for _ in range(nsects):
        sn[iptr + 1] = sn[iptr + 1] / scale
        sn[iptr + 2] = sn[iptr + 2] / (scale * scale)
        sd[iptr + 1] = sd[iptr + 1] / scale
        sd[iptr + 2] = sd[iptr + 2] / (scale * scale)
        iptr += 3
    return sn, sd


def _bilin2(sn: np.ndarray, sd: np.ndarray, nsects: int) -> tuple[np.ndarray, np.ndarray]:
    iptr = 0
    for _ in range(nsects):
        a0 = sd[iptr]
        a1 = sd[iptr + 1]
        a2 = sd[iptr + 2]
        scale = a2 + a1 + a0
        sd[iptr] = 1.0
        sd[iptr + 1] = 2.0 * (a0 - a2) / scale
        sd[iptr + 2] = (a2 - a1 + a0) / scale

        b0 = sn[iptr]
        b1 = sn[iptr + 1]
        b2 = sn[iptr + 2]
        sn[iptr] = (b2 + b1 + b0) / scale
        sn[iptr + 1] = 2.0 * (b0 - b2) / scale
        sn[iptr + 2] = (b2 - b1 + b0) / scale
        iptr += 3
    return sn, sd


def cc_measure_select(params: MAParams, tshift: float, dlnA: float, cc_max: float) -> Tuple[float, float, float]:
    if (
        (cc_max < params.cc_min)
        or (tshift < params.tshift_min)
        or (tshift > params.tshift_max)
        or (dlnA < params.dlna_min)
        or (dlnA > params.dlna_max)
    ):
        return 0.0, 0.0, cc_max
    return tshift, dlnA, cc_max


def interpolate_dat_and_syn(
    params: MAParams,
    data: np.ndarray,
    syn: np.ndarray,
    syn_phydisp: np.ndarray,
    tstart: float,
    tend: float,
    t0: float,
    dt: float,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, int, int]:
    nlen = int(np.floor((tend - tstart) / dt) + 1)
    if nlen <= 1:
        raise ValueError("Check the length of the data and syn arrays")
    nlen = max(1, nlen)
    istart = int(np.floor((tstart - t0) / dt))

    times = tstart + np.arange(nlen) * dt
    idx = np.floor((times - t0) / dt).astype(int)
    valid = (idx >= 0) & (idx < len(data) - 1)
    idx_safe = np.clip(idx, 0, len(data) - 2)
    t1 = (idx_safe * dt) + t0
    w = (times - t1) / dt

    dat_win = np.zeros(nlen)
    syn_win = np.zeros(nlen)
    dat_win[valid] = data[idx_safe[valid]] + (data[idx_safe[valid] + 1] - data[idx_safe[valid]]) * w[valid]
    syn_win[valid] = syn[idx_safe[valid]] + (syn[idx_safe[valid] + 1] - syn[idx_safe[valid]]) * w[valid]
    syn_win_phydisp = np.zeros(nlen)
    if params.use_physical_dispersion:
        syn_win_phydisp[valid] = syn_phydisp[idx_safe[valid]] + (
            syn_phydisp[idx_safe[valid] + 1] - syn_phydisp[idx_safe[valid]]
        ) * w[valid]

    return dat_win, syn_win, syn_win_phydisp, nlen, istart


def _compute_cc_core(
    syn: np.ndarray,
    data: np.ndarray,
    nlen: int,
    dt: float,
    tshift_min: float,
    tshift_max: float,
    lpt: int,
):
    ishift = 0
    cc_max = 0.0
    i_left = -int(nlen / 2.0)
    i_right = int(nlen / 2.0)

    syn0 = syn[:nlen]
    data0 = data[:nlen]
    norm_s = np.sqrt(np.sum(syn0 * syn0))
    if norm_s == 0.0:
        return ishift, ishift * dt, cc_max

    shifts = np.arange(i_left, i_right + 1, dtype=np.int64)

    nfft = select_nfft(2 * nlen - 1, lpt)
    data_fft = np.fft.rfft(data0, nfft)
    syn_rev_fft = np.fft.rfft(syn0[::-1], nfft)
    conv_full = np.fft.irfft(data_fft * syn_rev_fft, nfft)[: 2 * nlen - 1]
    raw_cc = conv_full[shifts + (nlen - 1)]

    data_sq = data0 * data0
    prefix = np.empty(nlen + 1, dtype=np.float64)
    prefix[0] = 0.0
    prefix[1:] = np.cumsum(data_sq)
    id_left = np.maximum(0, shifts)
    id_right = np.minimum(nlen, nlen + shifts)
    segment_power = prefix[id_right] - prefix[id_left]
    norms = norm_s * np.sqrt(segment_power)

    cc = np.full_like(raw_cc, -np.inf, dtype=np.float64)
    valid_norm = norms > 0.0
    cc[valid_norm] = raw_cc[valid_norm] / norms[valid_norm]

    shift_sec = shifts.astype(np.float64) * dt
    valid = (shift_sec >= tshift_min) & (shift_sec <= tshift_max) & (cc > 0.0)
    if np.any(valid):
        cc_masked = np.where(valid, cc, -np.inf)
        i_best = int(np.argmax(cc_masked))
        cc_max = float(cc_masked[i_best])
        ishift = int(shifts[i_best])

    tshift = ishift * dt
    return ishift, tshift, cc_max


def _compute_cc_with_params(params: MAParams, syn: np.ndarray, data: np.ndarray, nlen: int, dt: float) -> Tuple[int, float, float, float]:
    ishift, tshift, cc_max = _compute_cc_core(
        syn, data, nlen, dt, params.tshift_min, params.tshift_max, params.lpt
    )
    if np.sum(syn[:nlen] * syn[:nlen]) == 0.0:
        dlnA = 0.0
    else:
        dlnA = 0.5 * np.log(
            np.sum(data[:nlen] * data[:nlen]) / np.sum(syn[:nlen] * syn[:nlen])
        )
    return ishift, tshift, dlnA, cc_max


def compute_cc(*args) -> Tuple[int, float, float, float]:
    if len(args) == 5 and isinstance(args[0], MAParams):
        params, syn, data, nlen, dt = args
        return _compute_cc_with_params(params, syn, data, nlen, dt)
    if len(args) == 4:
        syn, data, nlen, dt = args
        params = MAParams()
        return _compute_cc_with_params(params, syn, data, nlen, dt)
    raise TypeError("compute_cc expects (params, syn, data, nlen, dt) or (syn, data, nlen, dt)")


def deconstruct_dat_cc(dat_dtw: np.ndarray, ishift: int, nlen: int, dlnA: float) -> np.ndarray:
    dat_dtw_cc = dat_dtw.copy()
    idx = np.arange(nlen)
    j = idx + ishift
    mask = (j >= 1) & (j <= nlen - 2)
    dat_dtw_cc[mask] = dat_dtw[j[mask]]
    if ishift < 0:
        fill_index = -ishift + 1
        if 0 <= fill_index < nlen:
            dat_dtw_cc[: -ishift + 1] = dat_dtw_cc[fill_index]
    if ishift > 0:
        fill_index = nlen - ishift - 2
        if 0 <= fill_index < nlen:
            dat_dtw_cc[nlen - ishift - 1 : nlen] = dat_dtw_cc[fill_index]
    dat_dtw_cc *= np.exp(-dlnA)
    return dat_dtw_cc


def reconstruct_syn_cc(syn_dtw: np.ndarray, ishift: int, nlen: int, dlnA: float) -> Tuple[np.ndarray, np.ndarray]:
    syn_dtw_cc_dt = syn_dtw.copy()
    idx = np.arange(nlen)
    j = idx - ishift
    mask = (j >= 1) & (j <= nlen - 2)
    syn_dtw_cc_dt[mask] = syn_dtw[j[mask]]
    if ishift > 0:
        if ishift + 2 < nlen:
            syn_dtw_cc_dt[: ishift + 1] = syn_dtw_cc_dt[ishift + 2]
    if ishift < 0:
        fill_index = nlen + ishift - 2
        if 0 <= fill_index < nlen:
            syn_dtw_cc_dt[nlen + ishift - 1 : nlen] = syn_dtw_cc_dt[fill_index]
    syn_dtw_cc = syn_dtw_cc_dt * np.exp(dlnA)
    return syn_dtw_cc, syn_dtw_cc_dt


def compute_average_error(
    params: MAParams,
    data_dtw: np.ndarray,
    syn_dtw_cc: np.ndarray,
    syn_dtw_cc_dt: np.ndarray,
    nlen: int,
    dt: float,
    sigma_dt: float,
    sigma_dlnA: float,
) -> Tuple[float, float]:
    syn_vtw_cc = np.zeros(nlen)
    syn_vtw_cc[1:-1] = (syn_dtw_cc[2:nlen] - syn_dtw_cc[: nlen - 2]) / (2.0 * dt)
    syn_vtw_cc[0] = (syn_dtw_cc[1] - syn_dtw_cc[0]) / dt
    syn_vtw_cc[-1] = (syn_dtw_cc[-1] - syn_dtw_cc[-2]) / dt

    diff = data_dtw[:nlen] - syn_dtw_cc[:nlen]
    sigma_dt_top = np.sum(diff * diff)
    sigma_dlnA_top = sigma_dt_top
    sigma_dt_bot = np.sum(syn_vtw_cc[:nlen] ** 2)
    sigma_dlnA_bot = np.sum(syn_dtw_cc_dt[:nlen] ** 2)
    sigma_dt = np.sqrt(sigma_dt_top / sigma_dt_bot) if sigma_dt_bot > 0 else 1.0
    sigma_dlnA = np.sqrt(sigma_dlnA_top / sigma_dlnA_bot) if sigma_dlnA_bot > 0 else 1.0

    if params.error_type == 0:
        sigma_dt = 1.0
        sigma_dlnA = 1.0
    else:
        if sigma_dt < params.dt_sigma_min:
            sigma_dt = params.dt_sigma_min
        if sigma_dlnA < params.dlna_sigma_min:
            sigma_dlnA = params.dlna_sigma_min

    return sigma_dt, sigma_dlnA


def write_average_meas(*_args, **_kwargs) -> None:
    return


def write_trans(
    trans: np.ndarray,
    wvec: np.ndarray,
    i_right: int,
    idf_new: int,
    df: float,
    tshift: float,
    dlnA: float,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, float, float]:
    fnum = len(trans)
    phi_wt = np.zeros(fnum)
    abs_wt = np.zeros(fnum)
    dtau_wt = np.zeros(fnum)
    dlnA_wt = np.zeros(fnum)

    if i_right > 0:
        phi = np.angle(trans[:i_right])
        phi = np.unwrap(phi, discont=PHASE_STEP)
        amp = np.abs(trans[:i_right])

        phi_wt[:i_right] = phi
        abs_wt[:i_right] = amp

        dtau_wt[0] = tshift
        if i_right > 1:
            dtau_wt[1:i_right] = (-1.0 / wvec[1:i_right]) * phi[1:i_right] + tshift

        dlnA_wt[:i_right] = dlnA
        nz = amp > 0.0
        dlnA_wt[:i_right][nz] = np.log(amp[nz]) + dlnA

    if i_right > 0:
        dtau_wa = float(np.sum(dtau_wt[:i_right]) / i_right)
        dlnA_wa = float(np.sum(dlnA_wt[:i_right]) / i_right)
    else:
        dtau_wa = float(tshift)
        dlnA_wa = float(dlnA)
    return phi_wt, abs_wt, dtau_wt, dlnA_wt, dtau_wa, dlnA_wa


def reconstruct_syn(
    syn_dtwo: np.ndarray,
    wvec: np.ndarray,
    dtau_wt: np.ndarray,
    dlnA_wt: np.ndarray,
    i_right: int,
    dt: float,
    nlen: int,
    nfft: int,
) -> Tuple[np.ndarray, np.ndarray]:
    fnum = len(syn_dtwo)
    wseis_recon = np.zeros(fnum, dtype=np.complex128)
    omegas = wvec[:i_right]
    wseis_recon[:i_right] = syn_dtwo[:i_right] * np.exp(dlnA_wt[:i_right]) * np.exp(-CCI * omegas * dtau_wt[:i_right])
    syn_dtw_mt = _irfft(wseis_recon, nfft, dt)

    wseis_recon = np.zeros(fnum, dtype=np.complex128)
    wseis_recon[:i_right] = syn_dtwo[:i_right] * np.exp(-CCI * omegas * dtau_wt[:i_right])
    syn_dtw_mt_dt = _irfft(wseis_recon, nfft, dt)

    return syn_dtw_mt[:nlen], syn_dtw_mt_dt[:nlen]


@lru_cache(maxsize=32)
def _costaper_cached(ipoint: int) -> tuple[float, ...]:
    idx = np.arange(1, ipoint + 1, dtype=np.float64)
    tas = 1.0 - np.cos(2.0 * PI * idx / ipoint)
    return tuple((tas / np.sqrt(1.5)).tolist())


def costaper(ipoint: int) -> np.ndarray:
    return np.asarray(_costaper_cached(ipoint), dtype=np.float64)


def boxcar(ipoint: int) -> np.ndarray:
    return np.ones(ipoint)


def _root(a: np.ndarray, bb: np.ndarray, w: np.ndarray, n: int, ik: int, u: float, el: float, elam: float) -> float:
    epsi = 1.0e-15
    epsi1 = 5.0e-15

    while True:
        elam = 0.5 * (u + el)
        if abs(u - el) <= 1.5 * epsi1:
            return elam

        an = a[1] - elam
        b = 0.0
        bn = -1.0 / an
        iag = 1 if an >= 0.0 else 0

        for i in range(2, n + 1):
            x = abs(bb[i - 1]) / epsi if an == 0.0 else w[i - 1] / an
            an = a[i] - elam - x
            if an == 0.0:
                an = epsi
            bm = b
            b = bn
            bn = ((a[i] - elam) * b - bm * x - 1.0) / an
            if an >= 0.0:
                iag += 1

        if iag == ik:
            el = elam
        else:
            u = elam

        dlam = 1.0 / bn
        if abs(dlam) <= epsi1:
            dlam = np.copysign(epsi1, dlam)
        elam = elam - dlam
        if elam < u and elam > el:
            if abs(u - el) <= 1.5 * epsi1:
                return elam
            continue


def tsturm(nt: int, n: int, a: np.ndarray, b: np.ndarray, w: np.ndarray, nev: int, ipar: int) -> Tuple[np.ndarray, np.ndarray]:
    epsi = 1.0e-15
    epsi1 = 5.0e-15
    ev = np.full(max(n, nev) + 1, -1.0)
    r = np.zeros((n + 1, 2 * nev + ipar + 2))

    if n <= 0 or nev <= 0:
        return r[1:], ev[1 : nev + 1]

    a1 = np.zeros(n + 1)
    b1 = np.zeros(n + 1)
    w1 = np.zeros(n + 1)
    a1[1 : n + 1] = a[:n]
    b1[1 : n + 1] = b[:n]
    w1[1 : n + 1] = w[:n]

    umeps = 1.0 - epsi
    u = 1.0
    for ik in range(nev):
        if ik > 0:
            u = ev[ik] * umeps
        el = min(ev[ik + 1], u)

        while True:
            elam = 0.5 * (u + el)
            if abs(u - el) <= epsi1:
                break

            iag = 0
            q = a1[1] - elam
            if q >= 0:
                iag += 1

            for i in range(2, n + 1):
                x = abs(b1[i - 1]) / epsi if q == 0.0 else w1[i - 1] / q
                q = a1[i] - elam - x
                if q >= 0.0:
                    iag += 1
                if iag > nev:
                    break

            if iag >= ik + 1:
                el = elam
            else:
                u = elam
                continue

            if iag == ik + 1:
                break
            upper = min(iag, nev)
            for i in range(ik + 2, upper + 1):
                ev[i] = elam

        el = elam
        elam = _root(a1, b1, w1, n, ik + 1, u, el, elam)
        ev[ik + 1] = elam

        jk = 2 * (ik + 1) + ipar - 1
        r[1, jk] = 1.0
        if n >= 2:
            r[2, jk] = -(a1[1] - ev[ik + 1]) / b1[1]
            ddot = 1.0 + r[2, jk] * r[2, jk]
            jm1 = 2
            for j in range(3, n + 1):
                r[j, jk] = -((a1[jm1] - ev[ik + 1]) * r[jm1, jk] + b1[j - 2] * r[j - 2, jk]) / b1[jm1]
                ddot += r[j, jk] * r[j, jk]
                jm1 = j
        else:
            ddot = 1.0

        rnorm = np.sqrt(nt / (2.0 * ddot))
        r[1 : n + 1, jk] = r[1 : n + 1, jk] * rnorm

    return r[1 : n + 1, :], ev[1 : nev + 1]


@lru_cache(maxsize=32)
def _staper_cached(nt: int, fw: float, nev: int) -> tuple[tuple[float, ...], ...]:
    if nt <= 0 or nev <= 0:
        return tuple(tuple() for _ in range(max(nt, 0)))
    tapers = dpss(nt, NW=fw, Kmax=nev, sym=False, norm=2)
    out = tapers.T
    return tuple(tuple(float(value) for value in row) for row in out)


def staper(nt: int, fw: float, nev: int) -> np.ndarray:
    return np.asarray(_staper_cached(nt, fw, nev), dtype=np.float64)


def mt_measure_select(
    params: MAParams,
    nlen: int,
    tshift: float,
    i_pmax_syn: int,
    dtau_w: np.ndarray,
    err_dt: np.ndarray,
    dt: float,
    i_left: int,
    i_right: int,
    fstart: float,
    fend: float,
) -> Tuple[int, int, float, float, bool]:
    use_trace = True
    nfft = select_nfft(nlen, params.lpt)
    df = 1.0 / (dt * nfft)
    f_pmax = df * i_pmax_syn
    t_pmax = 1.0 / f_pmax
    wlen = dt * nlen

    if params.ncycle_in_window * t_pmax > wlen:
        use_trace = False

    fstart = max(fstart, params.ncycle_in_window / wlen)
    fend = min(fend, 1.0 / (2.0 * dt))

    ntaper = int(params.npi * 2.0)
    if ntaper > 10:
        ntaper = 10
    if ntaper < 1:
        ntaper = 10
    if use_trace and fstart >= fend - ntaper * df:
        use_trace = False

    i_left_old = i_left
    i_right_old = i_right
    fvec = df * np.arange(nfft)
    fseg = fvec[i_left_old : i_right_old + 1]
    left_hits = np.nonzero(fseg > fstart)[0]
    if left_hits.size > 0:
        i_left = i_left_old + int(left_hits[0]) - 1

    right_hits = np.nonzero(fseg > fend)[0]
    if right_hits.size > 0:
        i_right = i_left_old + int(right_hits[0]) - 1

    fstart = (i_left) * df
    fend = (i_right) * df

    if abs(tshift) <= 1.01 * dt:
        dtau_w[:] = 0.0
        use_trace = False

    if use_trace and i_right >= i_left:
        fseg = fvec[i_left:i_right + 1]
        valid = fseg > 0.0
        if np.any(valid):
            dtau_seg = np.abs(dtau_w[i_left:i_right + 1])[valid]
            err_seg = err_dt[i_left:i_right + 1][valid]
            fseg = fseg[valid]
            if np.any(dtau_seg > 1.0 / (params.dt_fac * fseg)):
                use_trace = False
            elif np.any(err_seg > 1.0 / (params.err_fac * fseg)):
                use_trace = False
            elif np.any(dtau_seg > params.dt_max_scale * abs(tshift)):
                use_trace = False

    return i_left, i_right, fstart, fend, use_trace


def mt_measure(
    params: MAParams,
    dat_dt: np.ndarray,
    syn_dt: np.ndarray,
    syn_dt_phydisp: np.ndarray,
    t0: float,
    dt: float,
    npts: int,
    tstart: float,
    tend: float,
) -> Tuple:
    if tstart < t0 or tend > t0 + (npts - 1) * dt or tstart >= tend:
        raise ValueError("Check tstart and tend")

    dat_dtw, syn_dtw, syn_dtw_phydisp, nlen, istart = interpolate_dat_and_syn(
        params, dat_dt, syn_dt, syn_dt_phydisp, tstart, tend, t0, dt
    )

    sfac1 = (2.0 / float(nlen)) ** 2
    ipwr_t = 10
    idx = np.arange(nlen)
    fac = 1.0 - np.cos(PI * idx / (nlen - 1)) ** ipwr_t
    syn_dtw[:nlen] *= fac
    dat_dtw[:nlen] *= fac
    if params.use_physical_dispersion:
        syn_dtw_phydisp[:nlen] *= fac

    ishift, tshift, dlnA, cc_max = compute_cc(params, syn_dtw, dat_dtw, nlen, dt)
    dat_dtw_cc = deconstruct_dat_cc(dat_dtw, ishift, nlen, dlnA)
    syn_dtw_cc, syn_dtw_cc_dt = reconstruct_syn_cc(syn_dtw, ishift, nlen, dlnA)

    sigma_dt_cc = 1.0
    sigma_dlnA_cc = 1.0
    sigma_dt_cc, sigma_dlnA_cc = compute_average_error(
        params, dat_dtw, syn_dtw_cc, syn_dtw_cc_dt, nlen, dt, sigma_dt_cc, sigma_dlnA_cc
    )

    write_average_meas(None, 2, tshift, dlnA, sigma_dt_cc, sigma_dlnA_cc)

    if params.is_mtm == 0:
        return (
            istart,
            dat_dtw,
            syn_dtw,
            syn_dtw_phydisp,
            nlen,
            tshift,
            sigma_dt_cc,
            dlnA,
            sigma_dlnA_cc,
            cc_max,
            syn_dtw_cc,
            1,
            1,
            0,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
        )

    nfft = select_nfft(nlen, params.lpt)
    df = 1.0 / (nfft * dt)
    dw = TWOPI * df
    fnum = nfft // 2 + 1
    df_new = 1.0 / (tend - tstart)
    idf_new = int(df_new / df) if df > 0 else 1
    wvec = np.zeros(fnum)
    wvec = dw * np.arange(fnum)

    syn_dtwo = _rfft(syn_dtw, nfft, dt)
    dat_dtwo = _rfft(dat_dtw_cc, nfft, dt)

    i_pmax_dat = int(np.argmax(np.abs(dat_dtwo))) + 1
    i_pmax_syn = int(np.argmax(np.abs(syn_dtwo))) + 1
    ampmax_unw = float(np.max(np.abs(syn_dtwo)))

    wtr_use_unw = ampmax_unw * params.wtr
    i_right = fnum
    if fnum > i_pmax_syn:
        mags = np.abs(syn_dtwo)
        tail = mags[i_pmax_syn:fnum]
        if tail.size > 0:
            below = tail <= abs(wtr_use_unw)
            if np.any(below):
                i_right = i_pmax_syn + int(np.argmax(below)) + 1

    if params.is_mtm == 1:
        ntaper = int(params.npi * 2.0)
    else:
        ntaper = 1

    if params.is_mtm == 1:
        tas = staper(nlen, params.npi, ntaper)
    elif params.is_mtm == 2:
        tas = costaper(nlen)[:, None]
    else:
        tas = boxcar(nlen)[:, None]

    tas_view = tas[:nlen, :ntaper]
    syn_tapered = syn_dtw[:nlen, None] * tas_view
    dat_tapered = dat_dtw_cc[:nlen, None] * tas_view
    syn_dtw_ho_all = np.fft.rfft(syn_tapered, nfft, axis=0) * dt
    dat_dtw_ho_all = np.fft.rfft(dat_tapered, nfft, axis=0) * dt

    top_mtm = np.sum(dat_dtw_ho_all * np.conj(syn_dtw_ho_all), axis=1)
    bot_mtm = np.sum(syn_dtw_ho_all * np.conj(syn_dtw_ho_all), axis=1)
    trans_w = np.zeros(fnum, dtype=np.complex128)
    if params.is_mtm != 1:
        syn_dtw_ho = syn_dtw_ho_all[:, 0]
        dat_dtw_ho = dat_dtw_ho_all[:, 0]
        ampmax = float(np.max(np.abs(syn_dtw_ho)))
        wtr_use = ampmax * params.wtr
        denom = np.where(np.abs(syn_dtw_ho) > abs(wtr_use), syn_dtw_ho, syn_dtw_ho + wtr_use)
        trans_w = dat_dtw_ho / denom

    if params.is_mtm != 1:
        phi_wt, abs_wt, dtau_w, dlnA_w, dtau_wa, dlnA_wa = write_trans(
            trans_w, wvec, i_right, idf_new, df, tshift, dlnA
        )
        syn_dtw_mt, syn_dtw_mt_dt = reconstruct_syn(
            syn_dtwo, wvec, dtau_w, dlnA_w, i_right, dt, nlen, nfft
        )
        return (
            istart,
            dat_dtw,
            syn_dtw,
            syn_dtw_phydisp,
            nlen,
            tshift,
            sigma_dt_cc,
            dlnA,
            sigma_dlnA_cc,
            cc_max,
            syn_dtw_cc,
            i_pmax_dat,
            i_pmax_syn,
            i_right,
            trans_w,
            dtau_w,
            dlnA_w,
            sigma_dt_cc,
            sigma_dlnA_cc,
            syn_dtw_mt,
            None,
        )

    wtr_mtm = 1.0e-10
    ampmax = np.max(np.abs(bot_mtm))
    wtr_use = ampmax * wtr_mtm ** 2
    denom = np.where(np.abs(bot_mtm) > abs(wtr_use), bot_mtm, bot_mtm + wtr_use)
    trans_mtm = top_mtm / denom

    phi_wt, abs_wt, dtau_w, dlnA_w, dtau_wa, dlnA_wa = write_trans(
        trans_mtm, wvec, i_right, idf_new, df, tshift, dlnA
    )
    syn_dtw_mt, syn_dtw_mt_dt = reconstruct_syn(
        syn_dtwo, wvec, dtau_w, dlnA_w, i_right, dt, nlen, nfft
    )

    sigma_dt = sigma_dt_cc
    sigma_dlnA = sigma_dlnA_cc
    write_average_meas(None, 1, dtau_wa, dlnA_wa, sigma_dt, sigma_dlnA)

    err_dt = np.zeros(fnum)
    err_dlnA = np.zeros(fnum)
    if ntaper > 1:
        phi_mul = np.zeros((fnum, ntaper))
        abs_mul = np.zeros((fnum, ntaper))
        dtau_mul = np.zeros((fnum, ntaper))
        dlnA_mul = np.zeros((fnum, ntaper))
        top_leave = top_mtm[:, None] - dat_dtw_ho_all * np.conj(syn_dtw_ho_all)
        bot_leave = bot_mtm[:, None] - syn_dtw_ho_all * np.conj(syn_dtw_ho_all)
        for iom in range(ntaper):
            top_mtm_jk = top_leave[:, iom]
            bot_mtm_jk = bot_leave[:, iom]
            ampmax = np.max(np.abs(bot_mtm_jk))
            wtr_use = ampmax * wtr_mtm ** 2
            trans_mtm_jk = np.where(
                np.abs(bot_mtm_jk) > abs(wtr_use),
                top_mtm_jk / bot_mtm_jk,
                top_mtm_jk / (bot_mtm_jk + wtr_use),
            )
            phi_mul[:, iom], abs_mul[:, iom], dtau_mul[:, iom], dlnA_mul[:, iom], _, _ = write_trans(
                trans_mtm_jk, wvec, i_right, idf_new, df, tshift, dlnA
            )
        err_phi = np.zeros(i_right)
        err_abs = np.zeros(i_right)
        for i in range(i_right):
            eph_ave = 0.0
            edt_ave = 0.0
            eabs_ave = 0.0
            eabs2_ave = 0.0
            for iom in range(ntaper):
                eph_iom = ntaper * phi_wt[i] - (ntaper - 1) * phi_mul[i, iom]
                edt_iom = ntaper * dtau_w[i] - (ntaper - 1) * dtau_mul[i, iom]
                eabs_iom = ntaper * abs_wt[i] - (ntaper - 1) * abs_mul[i, iom]
                eabs2_iom = ntaper * dlnA_w[i] - (ntaper - 1) * dlnA_mul[i, iom]
                eph_ave += eph_iom
                edt_ave += edt_iom
                eabs_ave += eabs_iom
                eabs2_ave += eabs2_iom
            eph_ave /= ntaper
            edt_ave /= ntaper
            eabs_ave /= ntaper
            eabs2_ave /= ntaper
            for iom in range(ntaper):
                err_phi[i] += (phi_mul[i, iom] - eph_ave) ** 2
                err_dt[i] += (dtau_mul[i, iom] - edt_ave) ** 2
                err_abs[i] += (abs_mul[i, iom] - eabs_ave) ** 2
                err_dlnA[i] += (dlnA_mul[i, iom] - eabs2_ave) ** 2
            err_phi[i] = np.sqrt(err_phi[i] / (ntaper * (ntaper - 1)))
            err_dt[i] = np.sqrt(err_dt[i] / (ntaper * (ntaper - 1)))
            if i == 0:
                err_dt[i] = LARGE_VAL
            err_abs[i] = np.sqrt(err_abs[i] / (ntaper * (ntaper - 1)))
            err_dlnA[i] = np.sqrt(err_dlnA[i] / (ntaper * (ntaper - 1)))

    return (
        istart,
        dat_dtw,
        syn_dtw,
        syn_dtw_phydisp,
        nlen,
        tshift,
        sigma_dt_cc,
        dlnA,
        sigma_dlnA_cc,
        cc_max,
        syn_dtw_cc,
        i_pmax_dat,
        i_pmax_syn,
        i_right,
        trans_mtm,
        dtau_w,
        dlnA_w,
        sigma_dt,
        sigma_dlnA,
        syn_dtw_mt,
        err_dt,
    )


def mt_adj(
    params: MAParams,
    istart: int,
    dat_dtw: np.ndarray,
    syn_dtw: np.ndarray,
    syn_dtw_phydisp: np.ndarray,
    nlen: int,
    npts: int,
    dt: float,
    tshift: float,
    dlnA: float,
    sigma_dt_cc: float,
    sigma_dlnA_cc: float,
    dtau_w: np.ndarray,
    dlnA_w: np.ndarray,
    err_dtau: np.ndarray,
    err_dlnA: np.ndarray,
    sigma_dt: float,
    sigma_dlnA: float,
    i_left: int,
    i_right: int,
) -> Tuple[np.ndarray, float, np.ndarray, float, np.ndarray]:
    window_chi = np.zeros(3 * (5 - 1) + 8)
    tr_adj_src = np.zeros(npts)
    am_adj_src = np.zeros(npts)

    time_window = np.ones(nlen)

    if 3 <= params.imeas <= 6:
        if params.use_physical_dispersion:
            syn_vtw = np.zeros(nlen)
            syn_vtw[1:-1] = (syn_dtw_phydisp[2:nlen] - syn_dtw_phydisp[: nlen - 2]) / (2.0 * dt)
            syn_vtw[0] = (syn_dtw_phydisp[1] - syn_dtw_phydisp[0]) / dt
            syn_vtw[-1] = (syn_dtw_phydisp[-1] - syn_dtw_phydisp[-2]) / dt
            Nnorm = dt * np.sum(syn_vtw[:nlen] * syn_vtw[:nlen])
            ft_bar_t = -syn_vtw[:nlen] / Nnorm
            Mnorm = dt * np.sum(syn_dtw_phydisp[:nlen] * syn_dtw_phydisp[:nlen])
            fa_bar_t = syn_dtw_phydisp[:nlen] / Mnorm
        else:
            syn_vtw = np.zeros(nlen)
            syn_vtw[1:-1] = (syn_dtw[2:nlen] - syn_dtw[: nlen - 2]) / (2.0 * dt)
            syn_vtw[0] = (syn_dtw[1] - syn_dtw[0]) / dt
            syn_vtw[-1] = (syn_dtw[-1] - syn_dtw[-2]) / dt
            Nnorm = dt * np.sum(syn_vtw[:nlen] * syn_vtw[:nlen])
            ft_bar_t = -syn_vtw[:nlen] / Nnorm
            Mnorm = dt * np.sum(syn_dtw[:nlen] * syn_dtw[:nlen])
            fa_bar_t = syn_dtw[:nlen] / Mnorm
    else:
        ft_bar_t = np.zeros(nlen)
        fa_bar_t = np.zeros(nlen)

    fp = np.zeros(nlen)
    fq = np.zeros(nlen)
    if params.is_mtm == 1:
        dtau_wtr = params.wtr * np.sum(np.abs(dtau_w[i_left:i_right])) / (i_right - i_left)
        dlnA_wtr = params.wtr * np.sum(np.abs(dlnA_w[i_left:i_right])) / (i_right - i_left)
        w_taper = np.zeros(i_right)
        idx = np.arange(i_left, i_right)
        w_taper[i_left:i_right] = 1.0 - np.cos(PI * (idx - i_left) / (i_right - i_left)) ** params.ipwr_w
        nfft = select_nfft(nlen, params.lpt)
        df = 1.0 / (nfft * dt)
        ffac = 2.0 * df * np.sum(w_taper[i_left:i_right])
        wp_taper = np.zeros(i_right)
        wq_taper = np.zeros(i_right)
        if params.error_type == 0:
            wp_taper[i_left:i_right] = w_taper[i_left:i_right] / ffac
            wq_taper[i_left:i_right] = w_taper[i_left:i_right] / ffac
        elif params.error_type == 1:
            wp_taper[i_left:i_right] = w_taper[i_left:i_right] / ffac / (sigma_dt**2)
            wq_taper[i_left:i_right] = w_taper[i_left:i_right] / ffac / (sigma_dlnA**2)
        else:
            err_t = err_dtau[i_left:i_right]
            err_t = err_t + (err_t < dtau_wtr) * dtau_wtr
            err_A = err_dlnA[i_left:i_right]
            err_A = err_A + (err_A < dlnA_wtr) * dlnA_wtr
            wp_taper[i_left:i_right] = w_taper[i_left:i_right] / ffac / (err_t**2)
            wq_taper[i_left:i_right] = w_taper[i_left:i_right] / ffac / (err_A**2)

        ntaper = int(params.npi * 2.0)
        tas = staper(nlen, params.npi, ntaper)
        nfft = select_nfft(nlen, params.lpt)
        fnum = len(dtau_w)
        tas_view = tas[:nlen, :ntaper]
        syn_source = syn_dtw_phydisp[:nlen] if params.use_physical_dispersion else syn_dtw[:nlen]
        syn_dtw_h_all = syn_source[:, None] * tas_view
        syn_vtw_h_all = np.zeros((nlen, ntaper))
        syn_vtw_h_all[1:-1, :] = (syn_dtw_h_all[2:, :] - syn_dtw_h_all[:-2, :]) / (2.0 * dt)
        syn_vtw_h_all[0, :] = (syn_dtw_h_all[1, :] - syn_dtw_h_all[0, :]) / dt
        syn_vtw_h_all[-1, :] = (syn_dtw_h_all[-1, :] - syn_dtw_h_all[-2, :]) / dt

        syn_dtw_ho_all = np.fft.rfft(syn_dtw_h_all, nfft, axis=0) * dt
        syn_vtw_ho_all = np.fft.rfft(syn_vtw_h_all, nfft, axis=0) * dt
        d_bot_mtm = np.sum(syn_dtw_ho_all * np.conj(syn_dtw_ho_all), axis=1)
        v_bot_mtm = np.sum(syn_vtw_ho_all * np.conj(syn_vtw_ho_all), axis=1)

        pwc_adj = np.zeros((fnum, ntaper), dtype=np.complex128)
        qwc_adj = np.zeros((fnum, ntaper), dtype=np.complex128)
        pwc_adj[:i_right, :] = syn_vtw_ho_all[:i_right, :] / v_bot_mtm[:i_right, None]
        qwc_adj[:i_right, :] = -syn_dtw_ho_all[:i_right, :] / d_bot_mtm[:i_right, None]
        if params.error_type == 0:
            pwc_adj[:i_right, :] *= wp_taper[:i_right, None]
            qwc_adj[:i_right, :] *= wq_taper[:i_right, None]
        else:
            pwc_adj[:i_right, :] *= (dtau_w[:i_right] * wp_taper[:i_right])[:, None]
            qwc_adj[:i_right, :] *= (dlnA_w[:i_right] * wq_taper[:i_right])[:, None]

        dtau_pj_t = np.fft.irfft(pwc_adj, nfft, axis=0) / dt
        dlnA_qj_t = np.fft.irfft(qwc_adj, nfft, axis=0) / dt
        fp = np.sum(tas_view * dtau_pj_t[:nlen, :], axis=1)
        fq = np.sum(tas_view * dlnA_qj_t[:nlen, :], axis=1)

    waveform_d2 = np.sum((dat_dtw[:nlen] * time_window) ** 2)
    waveform_s2 = np.sum((syn_dtw[:nlen] * time_window) ** 2)
    waveform_chi = np.sum(((dat_dtw[:nlen] - syn_dtw[:nlen]) * time_window) ** 2)

    i1 = istart + np.arange(nlen)
    mask = (i1 >= 0) & (i1 < len(tr_adj_src))
    i1_valid = i1[mask]
    
    if params.imeas in (1, 2):
        tr_adj_src[i1_valid] = -dat_dtw[mask] / waveform_d2 * time_window[mask]
        am_adj_src[i1_valid] = (syn_dtw[mask] - dat_dtw[mask]) * time_window[mask]
    elif params.imeas in (3, 4):
        tr_adj_src[i1_valid] = ft_bar_t[mask] * time_window[mask]
        am_adj_src[i1_valid] = fa_bar_t[mask] * time_window[mask]
    elif params.imeas in (5, 6):
        tr_adj_src[i1_valid] = -(tshift / sigma_dt_cc**2) * ft_bar_t[mask] * time_window[mask]
        am_adj_src[i1_valid] = -(dlnA / sigma_dlnA_cc**2) * fa_bar_t[mask] * time_window[mask]
    elif params.imeas in (7, 8):
        tr_adj_src[i1_valid] = fp[mask] * time_window[mask]
        am_adj_src[i1_valid] = fq[mask] * time_window[mask]

    nfft = select_nfft(nlen, params.lpt)
    df = 1.0 / (nfft * dt)
    if params.is_mtm == 1:
        window_chi[0] = 0.5 * 2.0 * df * np.sum((dtau_w[:i_right]) ** 2 * wp_taper[:i_right])
        window_chi[1] = 0.5 * 2.0 * df * np.sum((dlnA_w[:i_right]) ** 2 * wq_taper[:i_right])
    window_chi[2] = 0.5 * (tshift / sigma_dt_cc) ** 2
    window_chi[3] = 0.5 * (dlnA / sigma_dlnA_cc) ** 2
    if params.is_mtm == 1:
        wsum = np.sum(w_taper[:i_right])
        if wsum > 0.0:
            window_chi[4] = np.sum(dtau_w[:i_right] * w_taper[:i_right]) / wsum
            window_chi[5] = np.sum(dlnA_w[:i_right] * w_taper[:i_right]) / wsum
    window_chi[6] = tshift
    window_chi[7] = dlnA
    if params.is_mtm == 1:
        window_chi[8] = sigma_dt
        window_chi[9] = sigma_dlnA
    window_chi[10] = sigma_dt_cc
    window_chi[11] = sigma_dlnA_cc
    window_chi[12] = 0.5 * waveform_d2
    window_chi[13] = 0.5 * waveform_s2
    window_chi[14] = 0.5 * waveform_chi
    window_chi[15] = nlen * dt

    if params.imeas <= 2:
        tr_chi = 0.5 * waveform_chi
        am_chi = 0.5 * waveform_chi
    elif 3 <= params.imeas <= 6:
        tr_chi = window_chi[2]
        am_chi = window_chi[3]
    else:
        tr_chi = window_chi[0]
        am_chi = window_chi[1]

    return tr_adj_src, tr_chi, am_adj_src, am_chi, window_chi


def interpolate_syn(syn: np.ndarray, t1: float, dt1: float, npt1: int, t2: float, dt2: float, npt2: int) -> None:
    if t1 == t2 and dt1 == dt2 and npt1 == npt2:
        return

    syn1 = np.zeros(npt2, dtype=syn.dtype)
    times = t2 + np.arange(npt2) * dt2
    mask = (times > t1) & (times < t1 + (npt1 - 1) * dt1)
    ii = np.floor((times[mask] - t1) / dt1).astype(np.int64)
    tt = times[mask] - (ii * dt1 + t1)
    syn1[mask] = (syn[ii + 1] - syn[ii]) * tt / dt1 + syn[ii]
    syn[:npt2] = syn1
    if npt1 > npt2:
        syn[npt2:npt1] = 0.0


def taper_start(syn: np.ndarray, npt: int, itmax: int) -> None:
    if 2 * itmax > npt:
        raise ValueError("Check taper_start of adjoint source")
    wt = TWOPI / (2.0 * (itmax - 1))
    for i in range(itmax):
        syn[i] = syn[i] * (0.5 * (1.0 - np.cos(wt * i)))


def measure_adj_impl(
    data_in: np.ndarray,
    syn_in: np.ndarray,
    npts: int,
    t0: float,
    dt: float,
    fstart0: float,
    fend0: float,
    tstart: float,
    tend: float,
    tt: float,
    dtt: float,
    nn: int,
    chan: str,
    params: MAParams,
) -> Tuple[np.ndarray, float, float, float, np.ndarray]:
    data = data_in.copy()
    syn = syn_in.copy()
    syn_phydisp = syn.copy()

    if params.run_bandpass:
        bandpass(data, npts, dt, fstart0, fend0)
        bandpass(syn, npts, dt, fstart0, fend0)
        if params.use_physical_dispersion:
            bandpass(syn_phydisp, npts, dt, fstart0, fend0)

    cmp = chan[2:3]
    chan_syn = f"{chan}{cmp}"
    _ = chan_syn

    nwin = 0
    all_chi = 0.0
    adj_syn_all = np.zeros_like(syn)
    recon_cc_all = np.zeros_like(syn)

    nwin += 1
    window_chi_out = np.zeros(20)
    adj_syn_out = np.zeros(npts)
    window_chi = np.zeros_like(window_chi_out)
    window_chi[16] = 0.5 * np.sum(data**2)
    window_chi[17] = 0.5 * np.sum(syn**2)
    window_chi[18] = 0.5 * np.sum((data - syn) ** 2)
    window_chi[19] = npts * dt

    (
        istart,
        data_dtw,
        syn_dtw,
        syn_dtw_phydisp,
        nlen,
        tshift,
        sigma_dt_cc,
        dlnA,
        sigma_dlnA_cc,
        cc_max,
        syn_dtw_cc,
        i_pmax_dat,
        i_pmax_syn,
        i_right,
        trans_mtm,
        dtau_w,
        dlnA_w,
        sigma_dt,
        sigma_dlnA,
        syn_dtw_mt,
        err_dt,
    ) = mt_measure(params, data, syn, syn_phydisp, t0, dt, npts, tstart, tend)

    i_left_adj = 0
    i_right_adj = min(i_right, nlen)

    if params.is_mtm == 1:
        fstart = fstart0
        fend = fend0
        i_left = 0
        i_right_idx = min(i_right, len(dtau_w))
        err_dt = err_dt if err_dt is not None else np.zeros_like(dtau_w)
        i_left, i_right_idx, fstart, fend, use_trace = mt_measure_select(
            params,
            nlen,
            tshift,
            i_pmax_syn,
            dtau_w,
            err_dt,
            dt,
            i_left,
            i_right_idx,
            fstart,
            fend,
        )
        i_left_adj = i_left
        i_right_adj = min(i_right_idx, nlen)
        if not use_trace:
            params.imeas = params.imeas0 - 2
            params.is_mtm = 3
            (
                istart,
                data_dtw,
                syn_dtw,
                syn_dtw_phydisp,
                nlen,
                tshift,
                sigma_dt_cc,
                dlnA,
                sigma_dlnA_cc,
                cc_max,
                syn_dtw_cc,
                i_pmax_dat,
                i_pmax_syn,
                i_right,
                trans_mtm,
                dtau_w,
                dlnA_w,
                sigma_dt,
                sigma_dlnA,
                syn_dtw_mt,
                err_dt,
            ) = mt_measure(params, data, syn, syn_phydisp, t0, dt, npts, tstart, tend)
            i_left_adj = 0
            i_right_adj = min(i_right, nlen)

    if params.imeas >= 5:
        tshift, dlnA, cc_max = cc_measure_select(params, tshift, dlnA, cc_max)

    tr_chi_out = 0.0
    am_chi_out = 0.0
    if params.compute_adjoint_source:
        tr_adj_src, tr_chi, am_adj_src, am_chi, window_chi_local = mt_adj(
            params,
            istart,
            data_dtw,
            syn_dtw,
            syn_dtw_phydisp,
            nlen,
            npts,
            dt,
            tshift,
            dlnA,
            sigma_dt_cc,
            sigma_dlnA_cc,
            dtau_w if dtau_w is not None else np.zeros(nlen),
            dlnA_w if dlnA_w is not None else np.zeros(nlen),
            err_dt if err_dt is not None else np.zeros(nlen),
            np.zeros(nlen),
            sigma_dt if sigma_dt is not None else sigma_dt_cc,
            sigma_dlnA if sigma_dlnA is not None else sigma_dlnA_cc,
            i_left_adj,
            i_right_adj,
        )
        if params.imeas % 2 == 1:
            adj_syn_all += tr_adj_src
            all_chi += tr_chi
        else:
            adj_syn_all += am_adj_src
            all_chi += am_chi
        if params.imeas >= 7 and syn_dtw_mt is not None:
            recon_cc_all[istart : istart + nlen] += syn_dtw_mt[:nlen]
        else:
            recon_cc_all[istart : istart + nlen] += syn_dtw_cc[:nlen]
        window_chi[: len(window_chi_local)] = window_chi_local
        window_chi[16] = 0.5 * np.sum(data**2)
        window_chi[17] = 0.5 * np.sum(syn**2)
        window_chi[18] = 0.5 * np.sum((data - syn) ** 2)
        window_chi[19] = npts * dt
        window_chi_out[:] = window_chi
        tr_chi_out = tr_chi
        am_chi_out = am_chi

    if params.compute_adjoint_source:
        bandpass(adj_syn_all, npts, dt, fstart0, fend0)
        interpolate_syn(adj_syn_all, t0, dt, npts, tt, dtt, nn)
        itmax = int(params.tshort / dtt)
        taper_start(adj_syn_all, nn, itmax)
        adj_syn_out[:] = adj_syn_all[:nn]

    return window_chi_out, tr_chi_out, am_chi_out,tshift, adj_syn_out



def measure_adj(t0_inp: float,dt_inp: float,npt_inp: int,
                t0_syn: float,dt_syn: float,npt_syn: int,
                tstart: float,tend: float,imeas:int,
                tlong: float,tshort: float,verbose:bool,
                obs_data: np.ndarray,syn_data: np.ndarray,
                compute_adj_source: bool = True,
                run_bandpass: bool = False,display_details: bool = False,
                output_measure_files: bool = False,
                tshift_min: float = -4.5,tshift_max: float = 4.5,
                dlna_min: float = -1.5,dlna_max: float = 1.5,
                cc_min: float = 0.8,err_type: int = 1,
                dt_sigma_min: float = 1.,dlna_sigma_min: float = 0.5,
                itaper: int = 1,wtr: float = 0.02,npi: float = 2.5,
                dt_fac: float = 2.,err_fac: float = 2.5,
                dt_max_scale: float = 3.5,
                ncyle_in_window: float = 1.5,
                lpt:int = 15,
                use_physical_disp: bool = False) -> tuple[MeasureStats,np.ndarray]:
    
    if imeas <=6 and imeas <=3:
        itaper = 2
    params = MAParams(
        tlong=tlong,
        tshort=tshort,
        wtr=wtr,
        npi=npi,
        dt_fac=dt_fac,
        err_fac=err_fac,
        dt_max_scale=dt_max_scale,
        ncycle_in_window=ncyle_in_window,
        tshift_min=tshift_min,
        tshift_max=tshift_max,
        dlna_min=dlna_min,
        dlna_max=dlna_max,
        cc_min=cc_min,
        dt_sigma_min=dt_sigma_min,
        dlna_sigma_min=dlna_sigma_min,
        ntaper=0,
        ipwr_t=10,
        ipwr_w=10,
        error_type=err_type,
        imeas0=imeas,
        imeas=imeas,
        itaper=itaper,
        is_mtm0=0,
        is_mtm=0,
        display_details=display_details,
        output_measurement_files=output_measure_files,
        run_bandpass=run_bandpass,
        compute_adjoint_source=compute_adj_source,
        use_physical_dispersion=use_physical_disp,
        lpt=lpt
    )

    if params.imeas in (1, 2):
        params.is_mtm0 = 0
    elif 3 <= params.imeas <= 6:
        if params.itaper == 1:
            raise ValueError("Change ITAPER to 2/3 for CC measurements")
        params.is_mtm0 = params.itaper
    elif params.imeas in (7, 8):
        params.is_mtm0 = 1
    else:
        raise ValueError("imeas must be 1-8")

    params.is_mtm = params.is_mtm0
    fstart0 = 1.0 / params.tlong
    fend0 = 1.0 / params.tshort

    # sanity check for input data
    empty_data = np.sum(obs_data**2) == 0.0
    if empty_data:
        # If the observed data is empty, we can skip the actual measurement and return a default misfit and zero adjoint source.
        obs_data = np.ones_like(syn_data)
        syn_data = np.ones_like(syn_data)
    
    _, tr_chi_out, am_chi_out, tshift, adj =  \
        measure_adj_impl(
        obs_data,
        syn_data,
        npt_inp,
        t0_inp,
        dt_inp,
        fstart0,
        fend0,
        tstart,
        tend,
        t0_syn,
        dt_syn,
        npt_syn,
        'BX',
        params,
    )

    misfit = 0.
    if imeas <= 2:
        misfit = tr_chi_out
    elif imeas % 2 == 1:
        misfit = tr_chi_out
    else:
        misfit = am_chi_out

    stats = MeasureStats(
        adj_type = str(params.imeas),
        misfit=misfit,
        tstart=tstart,
        tend=tend,
        tr_chi=tr_chi_out,
        am_chi=am_chi_out,
        tshift = tshift
    )

    return stats, adj
    
