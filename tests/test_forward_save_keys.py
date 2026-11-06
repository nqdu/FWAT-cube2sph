"""Forward SAC saves read rotated traces from their in-memory keys."""

import os
import numpy as np
import pytest
from obspy.io.sac import SACTrace

from fwat.measure.noise_mc_preproc import NoiseMC_PreOP
from fwat.measure.rf_preproc import RF_PreOP
from fwat.measure.sks_preproc import SKS_PreOP
from fwat.measure.tele_preproc import Tele_PreOP


def _forward_operator(cls, tmp_path):
    op = cls.__new__(cls)
    op.evtid = "EVENT"
    op.syndir = str(tmp_path / "synthetics")
    op.DATA_DIR = str(tmp_path / "observations")
    op.evla = 0.0
    op.evlo = 0.0
    op.evdp = 10.0
    op.netwk = ["AA"]
    op.stnm = ["BB"]
    op.stla = [1.0]
    op.stlo = [2.0]
    op.chcode = "BX"
    op.components = ["R", "Z"]
    op.ncomp = 2
    op.nsta_loc = 1
    op._istart = 0
    op.myrank = 0
    op.t0_syn = 0.0
    op.dt_syn = 0.1
    op.npt_syn = 128
    op.t_inj = 0.0
    op.t_ref = [0.0]
    op.seismogram = {}
    for component, amplitude in (("R", 1.0), ("Z", 2.0)):
        name = f"AA.BB.BX{component}.sem.npy"
        time = np.arange(op.npt_syn) * op.dt_syn
        op.seismogram[os.path.join(op.syndir, name)] = np.column_stack(
            (time, amplitude * np.ones(op.npt_syn)))
    return op


@pytest.mark.parametrize("cls", [Tele_PreOP, SKS_PreOP])
def test_forward_sac_uses_in_memory_trace(cls, tmp_path):
    op = _forward_operator(cls, tmp_path)
    if cls is Tele_PreOP:
        op.adjsrc_type = "cross-conv"

    op.save_forward()

    for component, amplitude in (("R", 1.0), ("Z", 2.0)):
        path = tmp_path / "observations" / "EVENT" / f"AA.BB.BX{component}.sac"
        trace = SACTrace.read(str(path))
        np.testing.assert_allclose(trace.data, amplitude)


def test_rf_forward_sac_uses_in_memory_trace(tmp_path, monkeypatch):
    op = _forward_operator(RF_PreOP, tmp_path)
    op._tshift = 1.0
    op._f0 = [1.0]
    op._maxit = 1
    op.pdict = {"TIME_WINDOW": [1.0, 1.0]}
    monkeypatch.setattr("fwat.measure.rf_preproc.deconit",
                        lambda radial, vertical, *args: radial.copy())

    op.save_forward()

    path = tmp_path / "observations" / "EVENT" / "AA.BB.BXR.F1.0.rf.sac"
    assert path.is_file()
    assert SACTrace.read(str(path)).data.shape == (op.npt_syn,)


def test_noise_forward_sac_uses_source_and_receiver_keys(tmp_path):
    op = NoiseMC_PreOP.__new__(NoiseMC_PreOP)
    op.evtid = "EVENT"
    op.mod = "M00"
    op.SOLVER = str(tmp_path / "solver")
    op.DATA_DIR = str(tmp_path / "observations")
    op.evla = 0.0
    op.evlo = 0.0
    op.evdp = 10.0
    op.sta_names = ["AA.BB"]
    op.sta_names_comp = [{"AA.BB"}]
    op.stainfo = {"AA.BB": [2.0, 1.0, 0.0, 0.0, 100.0]}
    op.cc_comps = ["ZZ"]
    op.ncomp = 1
    op.nprocs = 1
    op.myrank = 0
    op.chcode = "BX"
    op.t0_syn = 0.0
    op.dt_syn = 0.1
    op.npt_syn = 128
    op.pdict = {"USE_EGF": True}
    source_dir = f"{op.SOLVER}/{op.mod}/{op.evtid}_Z"
    name = "AA.BB.BXZ.sem.npy"
    op.seismogram = {
        os.path.join(source_dir, name): np.column_stack(
            (np.arange(op.npt_syn) * op.dt_syn, np.ones(op.npt_syn)))
    }

    op.save_forward()

    path = tmp_path / "observations" / "EVENT_Z" / "AA.BB.BXZ.sac"
    np.testing.assert_allclose(SACTrace.read(str(path)).data, 1.0)
