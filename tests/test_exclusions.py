import copy

import numpy as np
import pytest

from fwat.adjoint.MeasureStats import MeasureStats
from fwat.adjoint.cc_misfit import measure_adj_cc_impl
from fwat.adjoint.exp_phase_misfit import measure_adj_exphase_impl
from fwat.adjoint.multitaper import measure_adj_mt_impl
from fwat.measure.FwatPreOP import FwatPreOP
from fwat.measure.measure import measure_adj
from fwat.measure.noise_mc_preproc import NoiseMC_PreOP


@pytest.fixture
def traces():
    nt = 256
    dt = 0.1
    t = np.arange(nt) * dt
    syn = np.exp(-((t - 9.0) / 1.0) ** 2)
    obs = np.exp(-((t - 9.5) / 1.0) ** 2)
    return obs, syn, nt, dt


@pytest.mark.parametrize("measure", ["cc", "exp_phase", "mt", "legacy_mt"])
def test_quality_check_runs_only_when_requested(traces, measure):
    obs, syn, nt, dt = traces
    if measure == "cc":
        call = lambda check: measure_adj_cc_impl(
            obs, syn, 0., dt, nt, 1., 20., 6., 13.,
            cc_min=1.1, check_excluded=check)
    elif measure == "exp_phase":
        call = lambda check: measure_adj_exphase_impl(
            obs, syn, 0., dt, nt, 1., 20., 6., 13.,
            cc_min=1.1, check_excluded=check)
    elif measure == "mt":
        call = lambda check: measure_adj_mt_impl(
            obs, syn, 0., dt, nt, 6., 13., 1., 20.,
            cc_min=1.1, check_excluded=check)
    else:
        call = lambda check: measure_adj(
            0., dt, nt, 0., dt, nt, 6., 13., 7, 20., 1., False,
            obs, syn, cc_min=1.1, check_excluded=check)

    excluded, excluded_adj = call(True)
    assert excluded.exclude
    assert excluded.misfit == 0
    assert np.count_nonzero(excluded_adj) == 0

    included, included_adj = call(False)
    assert not included.exclude
    assert included.misfit > 0
    assert np.count_nonzero(included_adj) > 0


def test_later_numeric_cc_can_measure_shift_outside_initial_bounds():
    nt = 256
    dt = 0.1
    t = np.arange(nt) * dt
    syn = np.exp(-((t - 9.) / 1.) ** 2)
    obs = np.exp(-((t - 11.) / 1.) ** 2)
    stats, adj = measure_adj(
        0., dt, nt, 0., dt, nt, 6., 15., 5, 20., 1., False,
        obs, syn, tshift_min=-0.5, tshift_max=0.5,
        check_excluded=False, itaper=2)
    assert not stats.exclude
    assert stats.tshift == pytest.approx(2.)
    assert stats.misfit > 0
    assert np.count_nonzero(adj) > 0


def test_excluded_codes_are_saved_and_reused(tmp_path, traces, capsys):
    op = FwatPreOP.__new__(FwatPreOP)
    op.MISFIT = str(tmp_path)
    op.mod = "M00"
    op.evtid = "EVENT"
    op.meatype = "noise"
    op.iter_start = 0
    op.check_excluded = True
    op._record_exclusions = {}
    op.myrank = 0
    op.nprocs = 1
    op.adjsrc_type = "cc_time"

    bandname = "T001_T020"
    assert op.read_excluded_codes(bandname) == {}
    excluded_code = "EVENT_Z AA.BB.BXZ"
    accepted_code = "EVENT_Z CC.DD.BXZ"
    excluded = MeasureStats(adj_type="cc_time", code=excluded_code, exclude=True)
    assert copy.copy(excluded).exclude
    accepted = MeasureStats(adj_type="cc_time", code=accepted_code, misfit=1.)
    op._print_measure_info(bandname, [excluded, accepted])

    path = op.exclusion_file(bandname)
    assert path.read_text().splitlines() == [excluded_code]
    assert op.read_excluded_codes(bandname) == {excluded_code: True}
    assert not op._record_exclusions[bandname]
    output = (tmp_path / "M00" / "EVENT_T001_T020_noise_window_chi").read_text()
    assert accepted_code in output
    assert excluded_code not in output
    assert excluded_code not in capsys.readouterr().out

    noise = NoiseMC_PreOP.__new__(NoiseMC_PreOP)
    noise.adjsrc_type = "cc_time"
    obs, syn, nt, dt = traces
    stats, adj = noise._cal_adj_by_name(
        "AA.BB", obs, syn, 0., dt, nt, 0, 6., 13.,
        excluded_code, op.read_excluded_codes(bandname), False)
    assert stats.exclude and stats.misfit == 0
    assert np.count_nonzero(adj) == 0

    noise.MISFIT = str(tmp_path)
    noise.mod = "M00"
    noise.evtid = "EVENT"
    noise.meatype = "noise"
    noise.iter_start = 0
    noise.myrank = 0
    noise.nprocs = 1
    noise.sta_names = ["AA.BB"]
    noise._record_exclusions = {"T002_T020": True}
    noise._print_measure_info("T002_T020", [excluded, accepted])
    noise_path = noise.exclusion_file("T002_T020")
    assert noise_path.read_text().splitlines() == [excluded_code]
    noise_output = (tmp_path / "M00" / "EVENT_T002_T020_noise_window_chi").read_text()
    assert accepted_code in noise_output
    assert excluded_code not in noise_output
