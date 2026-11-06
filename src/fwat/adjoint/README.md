# Adjoint measurement parameters

The public functions in this directory accept waveform arrays followed by a
grouped `params` dictionary. The original positional equations are available
under names ending in `_impl` (for example, `measure_adj_cc_impl`). Optional
dictionary keys use the implementation's existing defaults.

```python
params = {
    "trace": {"t0": t0, "dt": dt, "nt": len(obs)},
    "window": {"start": tstart, "end": tend},
    "band": {"min_period": Tmin, "max_period": Tmax},
    "quality": {"cc_min": 0.8, "check_excluded": True},
    "options": {"return_type": "dt"},
}
stats, adj = measure_adj_cc(obs, syn, params)
```

`trace` and `window` are required for every method. `band` is required for
`measure_adj_cc`, `measure_adj_cc_dd`, `measure_adj_exphase`, and
`measure_adj_mt`. `measure_adj_cross_conv` uses only `t0` and `dt` in `trace`;
the other methods also require `nt`. The multitaper function maps `nt` to its
existing `npts` argument.

Double-difference calls use `windows` in place of `window`:

```python
params = {
    "trace": {"t0": t0, "dt": dt, "nt": len(obs_i)},
    "windows": {
        "first": {"start": tstart_i, "end": tend_i},
        "second": {"start": tstart_j, "end": tend_j},
    },
    "band": {"min_period": Tmin, "max_period": Tmax},
    "options": {"tshift_obs_user": dd_shift_ij},
}
stats, adj_i, adj_j = measure_adj_cc_dd(
    obs_i, syn_i, obs_j, syn_j, params)
```

`quality` contains the applicable shift, amplitude, correlation, and
uncertainty limits. `options` contains method-specific settings such as
`taper_ratio`, `measure_type`, or `npi`. Unknown sections and keys raise
`ValueError`, so misspelled settings do not silently use a default.

## Shared caller

Use `run_measurement` when the method is selected from configuration. It
returns named adjoints and optional diagnostics:

```python
from fwat.adjoint.caller import run_measurement

result = run_measurement(
    "cc_time", {"obs": obs, "syn": syn}, params)
stats = result.stats
adj = result.adjoints["trace"]
```

For `cc_dd` and `cross_conv`, provide `obs_first`, `syn_first`,
`obs_second`, and `syn_second`. Their adjoints are named `first` and `second`.
The cross-convolution arrays are in `result.diagnostics`. Numeric method names
such as `"5"` use the older `measure_adj` implementation.

The noise preprocessor checks SNR during its data quality control, before it
calls `run_measurement`. Low-SNR traces are marked excluded and recorded for
later iterations. The generic measurement caller does not inspect SNR.
