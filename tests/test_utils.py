"""
pytest tests for all functions in fwat/measure/utils.py
"""
import os
import tempfile
import numpy as np
import pytest

from fwat.measure.utils import (
    alloc_mpi_jobs,
    read_params,
    hann_taper,
    sac_cos_taper,
    interpolate_syn,
    bandpass,
    dif1,
    cumtrapz1,
    get_window_info,
    taper_window,
    get_source_loc,
    _geod2geoc,
    cal_dist_az_baz,
    cal_dist_baz,
    rotate_EN_to_RT,
    rotate_RT_to_EN,
)


# ---------------------------------------------------------------------------
# alloc_mpi_jobs
# ---------------------------------------------------------------------------

class TestAllocMpiJobs:
    def test_even_distribution(self):
        """All ranks receive exactly the same number of tasks."""
        ntasks, nprocs = 12, 4
        ranges = [alloc_mpi_jobs(ntasks, nprocs, r) for r in range(nprocs)]
        # verify contiguous coverage [0, ntasks)
        covered = set()
        for start, end in ranges:
            covered.update(range(start, end + 1))
        assert covered == set(range(ntasks))

    def test_uneven_distribution(self):
        """ntasks not divisible by nprocs – all tasks still assigned once."""
        ntasks, nprocs = 10, 3
        ranges = [alloc_mpi_jobs(ntasks, nprocs, r) for r in range(nprocs)]
        covered = []
        for start, end in ranges:
            covered.extend(range(start, end + 1))
        assert sorted(covered) == list(range(ntasks))

    def test_more_procs_than_tasks(self):
        """Excess ranks should receive zero-length allocations (start=-1, end=-2).
        range(start, end+1) = range(-1,-1) is empty; nsta_loc = end-start+1 = 0."""
        ntasks, nprocs = 2, 5
        idle_count = 0
        for r in range(nprocs):
            start, end = alloc_mpi_jobs(ntasks, nprocs, r)
            if start == -1:
                idle_count += 1
                assert end == -2
        assert idle_count == nprocs - ntasks

    def test_single_proc(self):
        """Single process should get all tasks."""
        ntasks = 7
        start, end = alloc_mpi_jobs(ntasks, 1, 0)
        assert start == 0
        assert end == ntasks - 1

    def test_single_task(self):
        """A single task goes to rank 0 only."""
        start0, end0 = alloc_mpi_jobs(1, 4, 0)
        assert start0 == 0 and end0 == 0
        for r in range(1, 4):
            s, e = alloc_mpi_jobs(1, 4, r)
            assert s == -1 and e == -2


# ---------------------------------------------------------------------------
# read_params
# ---------------------------------------------------------------------------

class TestReadParams:
    def test_reads_yaml(self, tmp_path):
        yaml_content = "key1: value1\nkey2: 42\n"
        p = tmp_path / "params.yaml"
        p.write_text(yaml_content)
        result = read_params(str(p))
        assert result == {"key1": "value1", "key2": 42}

    def test_nested_yaml(self, tmp_path):
        yaml_content = "outer:\n  inner: 3.14\n"
        p = tmp_path / "nested.yaml"
        p.write_text(yaml_content)
        result = read_params(str(p))
        assert result["outer"]["inner"] == pytest.approx(3.14)


# ---------------------------------------------------------------------------
# hann_taper
# ---------------------------------------------------------------------------

class TestHannTaper:
    def test_shape(self):
        win = hann_taper(100, 0.05)
        assert win.shape == (100,)

    def test_interior_is_one(self):
        """Middle of the window (far from edges) should be 1."""
        win = hann_taper(200, 0.05)
        assert np.allclose(win[20:180], 1.0)

    def test_edges_taper(self):
        """Edges should be less than 1 when p > 0."""
        win = hann_taper(200, 0.1)
        assert win[0] < 1.0
        assert win[-1] < 1.0

    def test_large_taper(self):
        """p=0.5 tapers half on each side; window remains non-negative."""
        win = hann_taper(100, 0.5)
        assert np.all(win >= 0.0)

    def test_values_in_range(self):
        win = hann_taper(256, 0.1)
        assert np.all(win >= 0.0) and np.all(win <= 1.0)


# ---------------------------------------------------------------------------
# sac_cos_taper
# ---------------------------------------------------------------------------

class TestSacCosTaper:
    def test_shape(self):
        win = sac_cos_taper(128, 0.05)
        assert len(win) == 128

    def test_endpoints_zero(self):
        """First and last samples must be 0."""
        win = sac_cos_taper(200, 0.1)
        assert win[0] == pytest.approx(0.0, abs=1e-6)
        assert win[-1] == pytest.approx(0.0, abs=1e-6)

    def test_interior_is_one(self):
        """Interior samples (well away from edges) should equal 1."""
        win = sac_cos_taper(200, 0.05)
        assert np.allclose(win[15:185], 1.0, atol=1e-6)

    def test_values_in_range(self):
        win = sac_cos_taper(300, 0.1)
        assert np.all(win >= 0.0) and np.all(win <= 1.0)


# ---------------------------------------------------------------------------
# interpolate_syn
# ---------------------------------------------------------------------------

class TestInterpolateSyn:
    def _linear_signal(self, t0, dt, npt):
        t = t0 + np.arange(npt) * dt
        return t - t0  # simple ramp starting at 0

    def test_same_grid_returns_same_values(self):
        dt = 0.01
        npt = 200
        t0 = 0.0
        data = np.sin(2 * np.pi * 1.0 * (t0 + np.arange(npt) * dt))
        out = interpolate_syn(data, t0, dt, npt, t0, dt, npt, max_percentage=0.0)
        assert out.shape == (npt,)

    def test_upsampling_preserves_smooth_signal(self):
        """Upsample a low-frequency sinusoid and check interior fidelity."""
        dt1, npt1, t1 = 0.02, 500, 0.0
        dt2, npt2, t2 = 0.01, 1000, 0.0
        t_coarse = t1 + np.arange(npt1) * dt1
        data = np.sin(2 * np.pi * 0.5 * t_coarse)
        out = interpolate_syn(data, t1, dt1, npt1, t2, dt2, npt2, max_percentage=0.0)
        t_fine = t2 + np.arange(npt2) * dt2
        idx = np.logical_and(t_fine > t1 + 1.0,
                              t_fine < t1 + (npt1 - 2) * dt1 - 1.0)
        expected = np.sin(2 * np.pi * 0.5 * t_fine[idx])
        assert np.allclose(out[idx], expected, atol=1e-3)

    def test_outside_range_is_zero(self):
        """Samples outside [t1, t1+(npt1-2)*dt1] must remain zero."""
        dt1, npt1, t1 = 0.01, 100, 5.0
        dt2, npt2, t2 = 0.01, 500, 0.0
        data = np.ones(npt1)
        out = interpolate_syn(data, t1, dt1, npt1, t2, dt2, npt2, max_percentage=0.0)
        # samples before t1
        early_idx = t2 + np.arange(npt2) * dt2 < t1
        assert np.all(out[early_idx] == 0.0)


# ---------------------------------------------------------------------------
# bandpass
# ---------------------------------------------------------------------------

class TestBandpass:
    def _make_signal(self, dt=0.01, duration=30.0):
        t = np.arange(0, duration, dt)
        # low-freq (0.05 Hz) + target (1 Hz) + high-freq (20 Hz)
        sig = (np.sin(2 * np.pi * 0.05 * t)
               + np.sin(2 * np.pi * 1.0 * t)
               + np.sin(2 * np.pi * 20.0 * t))
        return t, sig

    def test_output_shape(self):
        t, sig = self._make_signal()
        out = bandpass(sig, 0.01, 0.5, 2.0)
        assert out.shape == sig.shape

    def test_pass_band_preserved(self):
        """A pure in-band sinusoid should survive with most of its energy."""
        dt = 0.01
        t = np.arange(0, 20.0, dt)
        sig = np.sin(2 * np.pi * 1.0 * t)  # 1 Hz, well inside [0.5, 2] Hz
        out = bandpass(sig, dt, 0.5, 2.0, max_percentage=0.05)
        # use interior to avoid taper effects
        interior = slice(200, len(t) - 200)
        energy_ratio = np.var(out[interior]) / np.var(sig[interior])
        assert energy_ratio > 0.5

    def test_stop_band_attenuated(self):
        """A pure out-of-band frequency should lose most of its energy."""
        dt = 0.005
        t = np.arange(0, 40.0, dt)
        nyq = 0.5 / dt
        high_f = nyq * 0.9  # well above the filter band
        sig = np.sin(2 * np.pi * high_f * t)
        out = bandpass(sig, dt, 0.2, 1.0, max_percentage=0.05)
        interior = slice(400, len(t) - 400)
        energy_ratio = np.var(out[interior]) / (np.var(sig[interior]) + 1e-30)
        assert energy_ratio < 0.01

    def test_both_taper_types(self):
        dt = 0.01
        t = np.arange(0, 20.0, dt)
        sig = np.sin(2 * np.pi * 1.0 * t)
        out_hann = bandpass(sig, dt, 0.5, 2.0, type_='hann')
        out_cos  = bandpass(sig, dt, 0.5, 2.0, type_='cos')
        assert out_hann.shape == sig.shape
        assert out_cos.shape  == sig.shape

    def test_invalid_taper_type_raises(self):
        sig = np.ones(100)
        with pytest.raises(AssertionError):
            bandpass(sig, 0.01, 0.5, 2.0, type_='invalid')


# ---------------------------------------------------------------------------
# dif1
# ---------------------------------------------------------------------------

class TestDif1:
    def test_derivative_of_linear(self):
        """d/dt (a*t) = a, so interior samples should equal a/dt * dt = a."""
        dt = 0.01
        t = np.arange(0, 1.0, dt)
        a = 3.0
        data = a * t
        deriv = dif1(data, dt)
        # interior points only (central difference is exact for linear)
        assert np.allclose(deriv[1:-1], a, rtol=1e-10)

    def test_derivative_of_constant_is_zero(self):
        dt = 0.1
        data = np.full(50, 5.0)
        deriv = dif1(data, dt)
        assert np.allclose(deriv[1:-1], 0.0, atol=1e-12)

    def test_output_shape(self):
        data = np.random.randn(100)
        deriv = dif1(data, 0.01)
        assert deriv.shape == data.shape

    def test_endpoints_are_zero(self):
        """Central difference scheme leaves endpoints at 0."""
        data = np.sin(np.linspace(0, np.pi, 50))
        deriv = dif1(data, 0.01)
        assert deriv[0] == 0.0 and deriv[-1] == 0.0


# ---------------------------------------------------------------------------
# cumtrapz1
# ---------------------------------------------------------------------------

class TestCumtrapz1:
    def test_integral_of_constant(self):
        """∫ c dt from 0 to T ≈ c * T."""
        dt = 0.1
        n = 100
        c = 2.5
        data = np.full(n, c)
        result = cumtrapz1(data, dt)
        assert result.shape == data.shape
        assert result[0] == pytest.approx(0.0)
        assert result[-1] == pytest.approx(c * (n - 1) * dt, rel=1e-6)

    def test_integral_of_sine_is_negative_cosine(self):
        """∫ sin(t) dt = 1 - cos(t), starting at 0."""
        dt = 0.001
        t = np.arange(0, np.pi, dt)
        data = np.sin(t)
        result = cumtrapz1(data, dt)
        expected = 1 - np.cos(t)
        assert np.allclose(result, expected, atol=1e-4)

    def test_output_shape_and_initial_zero(self):
        data = np.random.randn(200)
        result = cumtrapz1(data, 0.01)
        assert result.shape == data.shape
        assert result[0] == pytest.approx(0.0)


# ---------------------------------------------------------------------------
# get_window_info
# ---------------------------------------------------------------------------

class TestGetWindowInfo:
    def test_basic_window(self):
        t0, dt = 0.0, 0.01
        tstart, tend = 1.0, 2.0
        lpt, rpt = get_window_info(t0, dt, tstart, tend)
        assert lpt == pytest.approx(int((tstart - t0) / dt), abs=1)
        assert rpt > lpt

    def test_window_length(self):
        t0, dt = 0.0, 0.01
        tstart, tend = 2.0, 5.0
        lpt, rpt = get_window_info(t0, dt, tstart, tend)
        nlen_expected = int((tend - tstart) / dt) + 1
        assert rpt - lpt == nlen_expected

    def test_reversed_window_raises(self):
        with pytest.raises(AssertionError):
            get_window_info(0.0, 0.01, 5.0, 3.0)


# ---------------------------------------------------------------------------
# taper_window
# ---------------------------------------------------------------------------

class TestTaperWindow:
    def test_basic_output(self):
        t0, dt, nt = 0.0, 0.01, 1000
        tstart, tend = 2.0, 5.0
        lpt, rpt, win = taper_window(t0, dt, nt, tstart, tend, p=0.05)
        assert 0 <= lpt < rpt <= nt
        nlen = rpt - lpt
        assert win.shape == (nlen,)

    def test_taper_values_in_range(self):
        t0, dt, nt = 0.0, 0.01, 1000
        lpt, rpt, win = taper_window(t0, dt, nt, 1.0, 4.0, p=0.1)
        assert np.all(win >= 0.0) and np.all(win <= 1.0)

    def test_clamps_to_signal_bounds(self):
        """Windows that extend beyond the signal should be clamped."""
        t0, dt, nt = 0.0, 0.01, 100  # signal ends at t=0.99
        lpt, rpt, win = taper_window(t0, dt, nt, -5.0, 200.0)
        assert lpt == 0
        assert rpt == nt

    def test_hann_and_cos_types(self):
        t0, dt, nt = 0.0, 0.01, 500
        lpt1, rpt1, w1 = taper_window(t0, dt, nt, 1.0, 3.0, type_='hann')
        lpt2, rpt2, w2 = taper_window(t0, dt, nt, 1.0, 3.0, type_='cos')
        assert w1.shape == w2.shape

    def test_invalid_taper_type_raises(self):
        with pytest.raises(AssertionError):
            taper_window(0.0, 0.01, 100, 0.1, 0.5, type_='bogus')

    def test_reversed_window_raises(self):
        with pytest.raises(AssertionError):
            taper_window(0.0, 0.01, 100, 5.0, 1.0)


# ---------------------------------------------------------------------------
# get_source_loc
# ---------------------------------------------------------------------------

class TestGetSourceLoc:
    def _write_source_file(self, tmp_path, lines):
        p = tmp_path / "sources.dat"
        p.write_text("\n".join(lines) + "\n")
        return str(p)

    def test_finds_existing_source(self, tmp_path):
        lines = [
            "EVT001 10.5 20.3 15.0",
            "EVT002  5.0 30.0  8.0",
        ]
        sf = self._write_source_file(tmp_path, lines)
        lat, lon, dep = get_source_loc("EVT001", sf)
        assert lat == pytest.approx(10.5)
        assert lon == pytest.approx(20.3)
        assert dep == pytest.approx(15.0)

    def test_finds_second_source(self, tmp_path):
        lines = [
            "EVT001 10.5 20.3 15.0",
            "EVT002  5.0 30.0  8.0",
        ]
        sf = self._write_source_file(tmp_path, lines)
        lat, lon, dep = get_source_loc("EVT002", sf)
        assert lat == pytest.approx(5.0)
        assert lon == pytest.approx(30.0)
        assert dep == pytest.approx(8.0)

    def test_missing_source_exits(self, tmp_path, capsys):
        lines = ["EVT001 10.5 20.3 15.0"]
        sf = self._write_source_file(tmp_path, lines)
        with pytest.raises(SystemExit):
            get_source_loc("EVT999", sf)


# ---------------------------------------------------------------------------
# _geod2geoc
# ---------------------------------------------------------------------------

class TestGeod2Geoc:
    def test_equator(self):
        """Geographic equator maps to geocentric equator."""
        assert _geod2geoc(0.0) == pytest.approx(0.0, abs=1e-10)

    def test_pole(self):
        """Geographic pole (90°) maps to geocentric pole (90°)."""
        assert _geod2geoc(90.0) == pytest.approx(90.0, abs=1e-6)

    def test_south_pole(self):
        assert _geod2geoc(-90.0) == pytest.approx(-90.0, abs=1e-6)

    def test_geocentric_slightly_smaller_than_geographic(self):
        """For oblate spheroid, geocentric lat < geographic lat in mid-latitudes."""
        geo_lat = 45.0
        gc_lat = _geod2geoc(geo_lat)
        assert gc_lat < geo_lat

    def test_spherical_earth_no_change(self):
        """With flattening=0 (sphere), geocentric == geographic."""
        lat = 37.5
        assert _geod2geoc(lat, flattening=0.0) == pytest.approx(lat, rel=1e-6)


# ---------------------------------------------------------------------------
# cal_dist_az_baz
# ---------------------------------------------------------------------------

class TestCalDistAzBaz:
    def test_same_point_is_zero(self):
        dist, az, baz = cal_dist_az_baz(0.0, 0.0, 0.0, 0.0)
        assert dist == pytest.approx(0.0, abs=1.0)

    def test_equatorial_distance(self):
        """1 degree along equator ≈ 111,195 m."""
        dist, az, baz = cal_dist_az_baz(0.0, 0.0, 0.0, 1.0)
        assert dist == pytest.approx(111_195.0, rel=0.01)

    def test_azimuth_northward(self):
        """Point due north should have azimuth ≈ 0°."""
        dist, az, baz = cal_dist_az_baz(0.0, 0.0, 10.0, 0.0)
        assert az == pytest.approx(0.0, abs=1.0)

    def test_azimuth_southward(self):
        """Point due south should have azimuth ≈ 180°."""
        dist, az, baz = cal_dist_az_baz(10.0, 0.0, 0.0, 0.0)
        assert az == pytest.approx(180.0, abs=1.0)

    def test_baz_is_roughly_opposite_az(self):
        """Back azimuth should be roughly opposite to azimuth for short paths."""
        dist, az, baz = cal_dist_az_baz(0.0, 0.0, 0.0, 1.0)
        assert abs(baz - (az + 180) % 360) < 2.0


# ---------------------------------------------------------------------------
# cal_dist_baz
# ---------------------------------------------------------------------------

class TestCalDistBaz:
    def test_returns_dist_and_baz(self):
        dist, baz = cal_dist_baz(0.0, 0.0, 0.0, 1.0)
        dist_full, _, baz_full = cal_dist_az_baz(0.0, 0.0, 0.0, 1.0)
        assert dist == pytest.approx(dist_full)
        assert baz == pytest.approx(baz_full)


# ---------------------------------------------------------------------------
# rotate_EN_to_RT  /  rotate_RT_to_EN
# ---------------------------------------------------------------------------

class TestRotateComponents:
    def _random_signal(self, n=256):
        rng = np.random.default_rng(42)
        return rng.standard_normal(n), rng.standard_normal(n)

    def test_round_trip_EN_to_RT_and_back(self):
        """EN → RT → EN should recover the originals."""
        ve, vn = self._random_signal()
        baz = 123.4
        vr, vt = rotate_EN_to_RT(ve, vn, baz)
        ve2, vn2 = rotate_RT_to_EN(vr, vt, baz)
        assert np.allclose(ve, ve2, atol=1e-12)
        assert np.allclose(vn, vn2, atol=1e-12)

    def test_round_trip_RT_to_EN_and_back(self):
        """RT → EN → RT should recover the originals."""
        vr, vt = self._random_signal()
        baz = 55.0
        ve, vn = rotate_RT_to_EN(vr, vt, baz)
        vr2, vt2 = rotate_EN_to_RT(ve, vn, baz)
        assert np.allclose(vr, vr2, atol=1e-12)
        assert np.allclose(vt, vt2, atol=1e-12)

    def test_energy_preserved(self):
        """Rotation is orthogonal: total energy is conserved."""
        ve, vn = self._random_signal()
        baz = 210.0
        vr, vt = rotate_EN_to_RT(ve, vn, baz)
        assert np.dot(ve, ve) + np.dot(vn, vn) == pytest.approx(
            np.dot(vr, vr) + np.dot(vt, vt), rel=1e-10
        )

    def test_baz_zero(self):
        """At baz=0: R ← -N, T ← -E (from the formulae with sin0=0, cos0=1)."""
        ve = np.array([1.0, 0.0])
        vn = np.array([0.0, 1.0])
        vr, vt = rotate_EN_to_RT(ve, vn, 0.0)
        assert np.allclose(vr, -vn, atol=1e-12)  # -ve*0 - vn*1
        assert np.allclose(vt, -ve, atol=1e-12)  # -ve*1 + vn*0

    def test_output_shape_preserved(self):
        n = 512
        ve, vn = np.ones(n), np.zeros(n)
        vr, vt = rotate_EN_to_RT(ve, vn, 45.0)
        assert vr.shape == (n,) and vt.shape == (n,)
