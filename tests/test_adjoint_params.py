"""The grouped adjoint API must agree with existing positional calls."""

import numpy as np
import pytest

from fwat.adjoint.cc_misfit import (measure_adj_cc, measure_adj_cc_impl,
                                   measure_adj_cc_dd, measure_adj_cc_dd_impl)
from fwat.adjoint.cross_conv import (measure_adj_cross_conv,
                                     measure_adj_cross_conv_impl)
from fwat.adjoint.exp_phase_misfit import (measure_adj_exphase,
                                           measure_adj_exphase_impl)
from fwat.adjoint.l2_misfit import measure_adj_l2, measure_adj_l2_impl
from fwat.adjoint.multitaper import measure_adj_mt, measure_adj_mt_impl


@pytest.fixture
def signals():
    dt = 0.1
    t = np.arange(256) * dt
    syn = np.exp(-((t - 9.0) / 1.0) ** 2)
    obs = np.exp(-((t - 9.4) / 1.0) ** 2)
    syn2 = np.exp(-((t - 8.8) / 1.2) ** 2)
    obs2 = np.exp(-((t - 9.3) / 1.2) ** 2)
    return obs, syn, obs2, syn2, dt


def assert_same_result(grouped, legacy):
    assert grouped[0].misfit == pytest.approx(legacy[0].misfit)
    assert grouped[0].exclude == legacy[0].exclude
    for new, old in zip(grouped[1:], legacy[1:]):
        np.testing.assert_allclose(new, old)


def test_single_trace_methods(signals):
    obs, syn, _, _, dt = signals
    trace = {"t0": 0.0, "dt": dt, "nt": len(obs)}
    window = {"start": 6.0, "end": 13.0}
    band = {"min_period": 1.0, "max_period": 20.0}
    quality = {"cc_min": 1.1, "check_excluded": False}

    assert_same_result(
        measure_adj_l2(obs, syn, {"trace": trace, "window": window}),
        measure_adj_l2_impl(obs, syn, 0.0, dt, len(obs), 6.0, 13.0),
    )
    assert_same_result(
        measure_adj_l2(obs=obs, syn=syn,
                       params={"trace": trace, "window": window}),
        measure_adj_l2_impl(obs, syn, 0.0, dt, len(obs), 6.0, 13.0),
    )
    assert_same_result(
        measure_adj_cc(obs, syn, {"trace": trace, "window": window,
                                  "band": band, "quality": quality}),
        measure_adj_cc_impl(obs, syn, 0.0, dt, len(obs), 1.0, 20.0, 6.0, 13.0,
                       cc_min=1.1, check_excluded=False),
    )
    assert_same_result(
        measure_adj_exphase(obs, syn, params={"trace": trace, "window": window,
                                               "band": band, "quality": quality}),
        measure_adj_exphase_impl(obs, syn, 0.0, dt, len(obs), 1.0, 20.0, 6.0, 13.0,
                            cc_min=1.1, check_excluded=False),
    )
    assert_same_result(
        measure_adj_mt(obs, syn, {"trace": trace, "window": window,
                                  "band": band, "quality": quality,
                                  "options": {"measure_type": "cc"}}),
        measure_adj_mt_impl(obs, syn, 0.0, dt, len(obs), 6.0, 13.0, 1.0, 20.0,
                       measure_type="cc", cc_min=1.1, check_excluded=False),
    )
    assert_same_result(
        measure_adj_mt(obs, syn, {"trace": trace, "window": window,
                                  "band": band, "quality": quality,
                                  "options": {"measure_type": "mt"}}),
        measure_adj_mt_impl(obs, syn, 0.0, dt, len(obs), 6.0, 13.0, 1.0, 20.0,
                       measure_type="mt", cc_min=1.1, check_excluded=False),
    )


def test_two_trace_methods(signals):
    obs, syn, obs2, syn2, dt = signals
    trace = {"t0": 0.0, "dt": dt, "nt": len(obs)}
    windows = {
        "first": {"start": 6.0, "end": 13.0},
        "second": {"start": 6.0, "end": 13.0},
    }
    assert_same_result(
        measure_adj_cc_dd(
            obs, syn, obs2, syn2,
            {"trace": trace, "windows": windows,
             "band": {"min_period": 1.0, "max_period": 20.0}}),
        measure_adj_cc_dd_impl(obs, syn, obs2, syn2, 0.0, dt, len(obs),
                          1.0, 20.0, 6.0, 13.0, 6.0, 13.0),
    )
    assert_same_result(
        measure_adj_cross_conv(
            obs, syn, obs2, syn2,
            {"trace": {"t0": 0.0, "dt": dt},
             "window": {"start": 6.0, "end": 13.0}}),
        measure_adj_cross_conv_impl(obs, syn, obs2, syn2,
                               0.0, dt, 6.0, 13.0),
    )


def test_invalid_dictionary_keys(signals):
    obs, syn, _, _, dt = signals
    params = {
        "trace": {"t0": 0.0, "dt": dt, "nt": len(obs)},
        "window": {"start": 6.0, "end": 13.0},
        "quality": {"cc_mim": 0.8},
    }
    with pytest.raises(ValueError, match="unknown keys"):
        measure_adj_l2(obs, syn, params)
