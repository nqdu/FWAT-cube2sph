"""Select an adjoint measurement and present its outputs consistently."""

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

import numpy as np

from .MeasureStats import MeasureStats
from .cc_misfit import measure_adj_cc, measure_adj_cc_dd
from .cross_conv import measure_adj_cross_conv
from .exp_phase_misfit import measure_adj_exphase
from .l2_misfit import measure_adj_l2
from .multitaper import measure_adj_mt
from .params import expand_params


@dataclass
class MeasurementResult:
    stats: MeasureStats
    adjoints: dict[str, np.ndarray]
    diagnostics: dict[str, Any]


_SINGLE_TRACE = {
    "l2": measure_adj_l2,
    "cc": measure_adj_cc,
    "cc_time": measure_adj_cc,
    "exp_phase": measure_adj_exphase,
    "mt": measure_adj_mt,
}


def _legacy_measure(imeas: int, waveforms: Mapping, params: Mapping):
    """Adapt numeric measurement types handled by the older Fortran port."""
    from fwat.measure.measure import measure_adj

    values = expand_params(params, "legacy")
    t0 = values.pop("t0")
    dt = values.pop("dt")
    nt = values.pop("nt")
    tstart = values.pop("tstart")
    tend = values.pop("tend")
    min_period = values.pop("min_period")
    max_period = values.pop("max_period")
    verbose = values.pop("verbose", False)
    itaper = values.pop("itaper", 2 if imeas == 5 else 1)
    return measure_adj(
        t0, dt, nt, t0, dt, nt, tstart, tend, imeas,
        max_period * 1.01, min_period * 0.99, verbose,
        waveforms["obs"], waveforms["syn"], itaper=itaper, **values)


def run_measurement(kind: str, waveforms: Mapping,
                    params: Mapping) -> MeasurementResult:
    """Run one method using named waveforms and grouped parameters."""

    if kind in _SINGLE_TRACE:
        stats, adj = _SINGLE_TRACE[kind](
            waveforms["obs"], waveforms["syn"], params)
        result = MeasurementResult(stats, {"trace": adj}, {})
    elif kind in ("cc_dd", "cc_time_dd"):
        stats, adj_first, adj_second = measure_adj_cc_dd(
            waveforms["obs_first"], waveforms["syn_first"],
            waveforms["obs_second"], waveforms["syn_second"], params)
        result = MeasurementResult(
            stats, {"first": adj_first, "second": adj_second}, {})
    elif kind == "cross_conv":
        stats, adj_first, adj_second, cc1, cc2 = measure_adj_cross_conv(
            waveforms["obs_first"], waveforms["syn_first"],
            waveforms["obs_second"], waveforms["syn_second"], params)
        result = MeasurementResult(
            stats, {"first": adj_first, "second": adj_second},
            {"cross_conv_first": cc1, "cross_conv_second": cc2})
    elif str(kind).isdigit():
        stats, adj = _legacy_measure(int(kind), waveforms, params)
        result = MeasurementResult(stats, {"trace": adj}, {})
    else:
        raise ValueError(f"Unknown measurement kind: {kind}")

    return result
