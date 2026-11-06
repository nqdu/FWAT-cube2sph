import numpy as np
import pytest

from fwat.adjoint.caller import run_measurement
from fwat.adjoint.cc_misfit import measure_adj_cc, measure_adj_cc_dd
from fwat.adjoint.cross_conv import measure_adj_cross_conv
from fwat.adjoint.exp_phase_misfit import measure_adj_exphase
from fwat.adjoint.l2_misfit import measure_adj_l2
from fwat.adjoint.multitaper import measure_adj_mt
from fwat.measure.noise_mc_preproc import NoiseMC_PreOP


@pytest.fixture
def traces():
    dt = 0.1
    t = np.arange(256) * dt
    syn = np.exp(-((t - 9.0) / 1.0) ** 2)
    obs = np.exp(-((t - 9.4) / 1.0) ** 2)
    syn2 = np.exp(-((t - 8.8) / 1.2) ** 2)
    obs2 = np.exp(-((t - 9.3) / 1.2) ** 2)
    return obs, syn, obs2, syn2, dt


def params(dt, nt):
    return {
        "trace": {"t0": 0.0, "dt": dt, "nt": nt},
        "window": {"start": 6.0, "end": 13.0},
        "band": {"min_period": 1.0, "max_period": 20.0},
        "quality": {"check_excluded": False},
    }


@pytest.mark.parametrize("kind, function", [
    ("l2", measure_adj_l2),
    ("cc_time", measure_adj_cc),
    ("exp_phase", measure_adj_exphase),
    ("mt", measure_adj_mt),
])
def test_single_trace_dispatch(traces, kind, function):
    obs, syn, _, _, dt = traces
    inputs = params(dt, len(obs))
    if kind == "l2":
        del inputs["band"]
        del inputs["quality"]
    result = run_measurement(kind, {"obs": obs, "syn": syn}, inputs)
    expected_stats, expected_adj = function(obs, syn, inputs)
    assert result.stats.misfit == pytest.approx(expected_stats.misfit)
    np.testing.assert_allclose(result.adjoints["trace"], expected_adj)
    assert result.diagnostics == {}


def test_pair_dispatch(traces):
    obs, syn, obs2, syn2, dt = traces
    inputs = params(dt, len(obs))
    inputs["windows"] = {
        "first": {"start": 6.0, "end": 13.0},
        "second": {"start": 6.0, "end": 13.0},
    }
    del inputs["window"]
    inputs["quality"] = {}
    waveforms = {
        "obs_first": obs, "syn_first": syn,
        "obs_second": obs2, "syn_second": syn2,
    }
    result = run_measurement("cc_dd", waveforms, inputs)
    stats, adj_first, adj_second = measure_adj_cc_dd(
        obs, syn, obs2, syn2, inputs)
    assert result.stats.misfit == pytest.approx(stats.misfit)
    np.testing.assert_allclose(result.adjoints["first"], adj_first)
    np.testing.assert_allclose(result.adjoints["second"], adj_second)


def test_cross_conv_diagnostics(traces):
    obs, syn, obs2, syn2, dt = traces
    inputs = {
        "trace": {"t0": 0.0, "dt": dt},
        "window": {"start": 6.0, "end": 13.0},
    }
    waveforms = {
        "obs_first": obs, "syn_first": syn,
        "obs_second": obs2, "syn_second": syn2,
    }
    result = run_measurement("cross_conv", waveforms, inputs)
    stats, first, second, cc1, cc2 = measure_adj_cross_conv(
        obs, syn, obs2, syn2, inputs)
    assert result.stats.misfit == pytest.approx(stats.misfit)
    np.testing.assert_allclose(result.adjoints["first"], first)
    np.testing.assert_allclose(result.adjoints["second"], second)
    np.testing.assert_allclose(result.diagnostics["cross_conv_first"], cc1)
    np.testing.assert_allclose(result.diagnostics["cross_conv_second"], cc2)


def test_legacy_numeric_dispatch(traces):
    obs, syn, _, _, dt = traces
    inputs = params(dt, len(obs))
    inputs["options"] = {"verbose": False, "itaper": 2}
    result = run_measurement("5", {"obs": obs, "syn": syn}, inputs)
    assert result.adjoints["trace"].shape == obs.shape
    assert np.isfinite(result.stats.misfit)


def test_runner_does_not_handle_noise_snr(traces):
    obs, syn, _, _, dt = traces
    inputs = params(dt, len(obs))
    inputs["quality"].update({"snr": 2.0, "snr_min": 3.0})
    with pytest.raises(ValueError, match="unknown keys"):
        run_measurement("cc", {"obs": obs, "syn": syn}, inputs)


def test_noise_snr_exclusion_is_recorded_and_reused(tmp_path, traces):
    obs, syn, _, _, dt = traces
    op = NoiseMC_PreOP.__new__(NoiseMC_PreOP)
    op.adjsrc_type = "cc_time"
    op.Tmin = [1.0]
    op.Tmax = [20.0]
    op.pdict = {}
    op.MISFIT = str(tmp_path)
    op.mod = "M00"
    op.evtid = "EVENT"
    op.meatype = "noise"
    op.iter_start = 0
    op.check_excluded = True
    op._record_exclusions = {}
    op.myrank = 0
    op.nprocs = 1
    op.sta_names = ["AA.BB"]
    code = "EVENT_Z AA.BB.BXZ"
    bandname = op._get_bandname(0)
    excluded_codes = op.read_excluded_codes(bandname)

    stats, adj = op._cal_adj_by_name(
        "AA.BB", obs, syn, 0., dt, len(obs), 0, 6., 13.,
        code, excluded_codes, op._record_exclusions[bandname],
        snr=1.0, snr_threshold=3.0)
    stats.code = code
    assert stats.exclude
    assert stats.misfit == 0.0
    assert np.count_nonzero(adj) == 0

    op._print_measure_info(bandname, [stats])
    assert op.exclusion_file(bandname).read_text().splitlines() == [code]
    excluded_codes = op.read_excluded_codes(bandname)
    stats, adj = op._cal_adj_by_name(
        "AA.BB", obs, syn, 0., dt, len(obs), 0, 6., 13.,
        code, excluded_codes, False, snr=10.0, snr_threshold=3.0)
    assert stats.exclude
    assert np.count_nonzero(adj) == 0

    stats, adj = op._cal_adj_by_name(
        "CC.DD", obs, syn, 0., dt, len(obs), 0, 6., 13.,
        "EVENT_Z CC.DD.BXZ", excluded_codes, False,
        snr=1.0, snr_threshold=3.0)
    assert not stats.exclude
    assert stats.misfit > 0
    assert np.count_nonzero(adj) > 0


def test_noise_pair_snr_excludes_both_traces(traces):
    obs, syn, obs2, syn2, dt = traces
    op = NoiseMC_PreOP.__new__(NoiseMC_PreOP)
    op.Tmin = [1.0]
    op.Tmax = [20.0]
    code = "EVENT_Z AA.BB.BXZ-CC.DD.BXZ"
    stats, first, second = op._cal_adj_by_name_dd(
        obs, obs2, syn, syn2, 0., dt, len(obs), 0,
        6., 13., 6., 13., 10.0, 1.0, 3.0,
        code, {}, True)
    assert stats.exclude
    assert np.count_nonzero(first) == 0
    assert np.count_nonzero(second) == 0

    stats, first, second = op._cal_adj_by_name_dd(
        obs, obs2, syn, syn2, 0., dt, len(obs), 0,
        6., 13., 6., 13., 10.0, 1.0, 3.0,
        code, {}, False)
    assert not stats.exclude
    assert stats.misfit > 0
    assert np.count_nonzero(first) > 0
    assert np.count_nonzero(second) > 0

    stats, first, second = op._cal_adj_by_name_dd(
        obs, obs2, syn, syn2, 0., dt, len(obs), 0,
        6., 13., 6., 13., 10.0, 10.0, 3.0,
        code, {code: True}, False)
    assert stats.exclude
    assert np.count_nonzero(first) == 0
    assert np.count_nonzero(second) == 0
