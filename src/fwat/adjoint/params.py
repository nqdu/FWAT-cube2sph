"""Dictionary input for the public adjoint measurement functions.

Waveform arrays stay as explicit arguments. The dictionary groups trace,
window, period band, quality checks, and method-specific options. Omitted
optional keys retain the defaults of the underlying measurement function.
"""

from collections.abc import Mapping
from typing import Any, TypedDict, Union


class TraceParams(TypedDict):
    t0: float
    dt: float
    nt: int


class CrossTraceParams(TypedDict):
    t0: float
    dt: float


class WindowParams(TypedDict):
    start: float
    end: float


class BandParams(TypedDict):
    min_period: float
    max_period: float


class DoubleWindowParams(TypedDict):
    first: WindowParams
    second: WindowParams


class QualityParams(TypedDict, total=False):
    tshift_min: float
    tshift_max: float
    dlna_min: float
    dlna_max: float
    cc_min: float
    check_excluded: bool
    dt_sigma_min: float
    dlna_sigma_min: float


class AdjointParams(TypedDict, total=False):
    trace: Union[TraceParams, CrossTraceParams]
    window: WindowParams
    windows: DoubleWindowParams
    band: BandParams
    quality: QualityParams
    options: dict[str, Any]


_COMMON_QUALITY = {
    "tshift_min", "tshift_max", "dlna_min", "dlna_max", "cc_min",
    "check_excluded",
}
_SIGMA_QUALITY = {"dt_sigma_min", "dlna_sigma_min"}
_METHODS = {
    "cc": (True, False, _COMMON_QUALITY | _SIGMA_QUALITY,
           {"return_type", "taper_ratio", "weight_by_uncertainty"}),
    "cc_dd": (True, True, _SIGMA_QUALITY,
              {"tshift_obs_user", "save_measurement", "taper_ratio"}),
    "exphase": (True, False, _COMMON_QUALITY,
                {"water", "taper_ratio"}),
    "mt": (True, False, _COMMON_QUALITY | _SIGMA_QUALITY,
           {"measure_type", "return_type", "itaper", "wtr", "npi",
            "dt_fac", "err_fac", "dt_max_scale", "ncycle_in_window",
            "error_type"}),
    "l2": (False, False, set(), {"taper_ratio"}),
    "cross_conv": (False, False, set(), {"taper_ratio"}),
    "legacy": (True, False, _COMMON_QUALITY | _SIGMA_QUALITY,
               {"verbose", "itaper"}),
}


def _section(value: Any, name: str, required: set[str], allowed: set[str]) -> dict:
    if not isinstance(value, Mapping):
        raise TypeError(f"params['{name}'] must be a mapping")
    missing = required - value.keys()
    unknown = value.keys() - allowed
    if missing:
        raise ValueError(f"params['{name}'] missing keys: {sorted(missing)}")
    if unknown:
        raise ValueError(f"params['{name}'] unknown keys: {sorted(unknown)}")
    return dict(value)


def expand_params(params: AdjointParams, method: str) -> dict[str, Any]:
    """Validate grouped parameters and expand them to implementation keywords."""
    if not isinstance(params, Mapping):
        raise TypeError("params must be a mapping")
    needs_band, double_window, quality_keys, option_keys = _METHODS[method]
    required_sections = {"trace", "windows" if double_window else "window"}
    if needs_band:
        required_sections.add("band")
    allowed_sections = required_sections | {"quality", "options"}
    missing = required_sections - params.keys()
    unknown = params.keys() - allowed_sections
    if missing:
        raise ValueError(f"params missing sections: {sorted(missing)}")
    if unknown:
        raise ValueError(f"params unknown sections: {sorted(unknown)}")

    trace_keys = {"t0", "dt"} if method == "cross_conv" else {"t0", "dt", "nt"}
    trace = _section(params["trace"], "trace", trace_keys, trace_keys)
    if method == "mt":
        trace["npts"] = trace.pop("nt")
    result = trace

    if double_window:
        windows = _section(params["windows"], "windows", {"first", "second"},
                           {"first", "second"})
        for label, number in (("first", 1), ("second", 2)):
            window = _section(windows[label], f"windows.{label}",
                              {"start", "end"}, {"start", "end"})
            result[f"tstart{number}"] = window["start"]
            result[f"tend{number}"] = window["end"]
    else:
        window = _section(params["window"], "window", {"start", "end"},
                          {"start", "end"})
        result["tstart"] = window["start"]
        result["tend"] = window["end"]

    if needs_band:
        result.update(_section(params["band"], "band",
                               {"min_period", "max_period"},
                               {"min_period", "max_period"}))
    result.update(_section(params.get("quality", {}), "quality", set(), quality_keys))
    result.update(_section(params.get("options", {}), "options", set(), option_keys))
    return result
