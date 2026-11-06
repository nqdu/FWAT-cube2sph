"""Summed adjoint traces reach rotation without temporary .npy files."""

from pathlib import Path
import os

import h5py
import numpy as np

from fwat.measure.FwatPreOP import FwatPreOP
from fwat.measure.noise_mc_preproc import NoiseMC_PreOP


def _rotation_file(tmp_path):
    (tmp_path / "rot_EVENT").write_text(
        "AA BB\n"
        "0 1 0\n"
        "-1 0 0\n"
        "0 0 1\n"
    )


def _trace(values):
    return np.column_stack((np.arange(len(values), dtype=float), values))


def _forward_input(syndir, scale):
    output_dir = Path(syndir) / "OUTPUT_FILES"
    output_dir.mkdir(parents=True)
    with h5py.File(output_dir / "seismograms.h5", "w") as stream:
        for component, values in {
            "X": [1., 2., 3.], "Y": [4., 5., 6.], "Z": [7., 8., 9.],
        }.items():
            stream.create_dataset(
                f"AA.BB.BX{component}.semd",
                data=_trace(scale * np.asarray(values)))


def _assert_forward_keys(seismogram, syndir, scale):
    expected = {
        "N": [4., 5., 6.], "E": [-1., -2., -3.], "Z": [7., 8., 9.],
    }
    for component, values in expected.items():
        key = os.path.join(syndir, f"AA.BB.BX{component}.sem.npy")
        np.testing.assert_allclose(seismogram[key][:, 1],
                                   scale * np.asarray(values))
        assert "OUTPUT_FILES" not in key


def _assert_rotated(syndir, east, north, vertical):
    sem = Path(syndir) / "SEM"
    for component, values in {
        "X": -np.asarray(east), "Y": north, "Z": vertical,
    }.items():
        trace = np.loadtxt(sem / f"AA.BB.BX{component}.adj")
        np.testing.assert_allclose(trace[:, 0], [0., 1., 2.])
        np.testing.assert_allclose(trace[:, 1], values)
    assert not list(sem.glob("*.npy"))


def test_base_forward_rotation_uses_logical_keys(tmp_path):
    _rotation_file(tmp_path)
    syndir = str(tmp_path / "base")
    _forward_input(syndir, 1.0)
    op = FwatPreOP.__new__(FwatPreOP)
    op.SRC_REC = str(tmp_path)
    op.evtid = "EVENT"
    op.syndir = syndir
    op.netwk = ["AA"]
    op.stnm = ["BB"]
    op._istart = 0
    op.nsta_loc = 1

    op._rotate_XYZ_to_ZNE()

    _assert_forward_keys(op.seismogram, syndir, 1.0)


def test_noise_forward_rotation_keeps_source_namespaces(tmp_path):
    _rotation_file(tmp_path)
    syndirs = [str(tmp_path / "noise_E"), str(tmp_path / "noise_N")]
    for syndir, scale in zip(syndirs, (1.0, 2.0)):
        _forward_input(syndir, scale)
    op = NoiseMC_PreOP.__new__(NoiseMC_PreOP)
    op.SRC_REC = str(tmp_path)
    op.evtid = "EVENT"
    op.syndirs = syndirs
    op.scomp_syn = ["E", "N"]
    op.sta_names = ["AA.BB"]
    op.nprocs = 1
    op.myrank = 0
    op.seismogram = {}

    op._rotate_XYZ_to_ZNE()

    for syndir, scale in zip(syndirs, (1.0, 2.0)):
        _assert_forward_keys(op.seismogram, syndir, scale)
    assert len(op.seismogram) == 6


def test_base_sum_rotates_from_memory(tmp_path):
    _rotation_file(tmp_path)
    op = FwatPreOP.__new__(FwatPreOP)
    op.myrank = 0
    op.nsta_loc = 1
    op._istart = 0
    op.netwk = ["AA"]
    op.stnm = ["BB"]
    op.chcode = "BX"
    op.components = ["E", "N", "Z"]
    op.bazd = np.zeros(1)
    op.npt_syn = 3
    op.t0_syn = 0.0
    op.dt_syn = 1.0
    op.Tmin = [1.0, 2.0]
    op.Tmax = [2.0, 3.0]
    op.SRC_REC = str(tmp_path)
    op.evtid = "EVENT"
    op.syndir = str(tmp_path / "base")
    op.seismogram_adj = {}

    for ib, scale in enumerate((1.0, 2.0)):
        band = op._get_bandname(ib)
        for component, values in {
            "E": [1., 2., 3.], "N": [4., 5., 6.], "Z": [7., 8., 9.],
        }.items():
            key = os.path.join(op.syndir,band,
                              f"AA.BB.BX{component}.adj.sem.npy")
            op.seismogram_adj[key] = _trace(scale * np.asarray(values))

    op.sum_adj_source()

    _assert_rotated(op.syndir, [3., 6., 9.], [12., 15., 18.],
                    [21., 24., 27.])
    assert len(op.seismogram_adj) == 9  # six band traces and three summed traces


def test_noise_sum_rotates_each_source_directory_from_memory(tmp_path):
    _rotation_file(tmp_path)
    op = NoiseMC_PreOP.__new__(NoiseMC_PreOP)
    op.myrank = 0
    op.nprocs = 1
    op.sta_names = ["AA.BB"]
    op.chcode = "BX"
    op.npt_syn = 3
    op.t0_syn = 0.0
    op.dt_syn = 1.0
    op.Tmin = [1.0]
    op.Tmax = [2.0]
    op.SRC_REC = str(tmp_path)
    op.evtid = "EVENT"
    op.scomp_syn = ["E", "N"]
    op.syndirs = [str(tmp_path / "noise_E"), str(tmp_path / "noise_N")]
    op.seismogram_adj = {}

    for syndir, scale in zip(op.syndirs, (1.0, 2.0)):
        band = op._get_bandname(0)
        for component, values in {
            "E": [1., 2., 3.], "N": [4., 5., 6.], "Z": [7., 8., 9.],
        }.items():
            key = os.path.join(syndir,band,
                              f"AA.BB.BX{component}.adj.sem.npy")
            op.seismogram_adj[key] = _trace(scale * np.asarray(values))

    op.sum_adj_source()

    _assert_rotated(op.syndirs[0], [1., 2., 3.], [4., 5., 6.],
                    [7., 8., 9.])
    _assert_rotated(op.syndirs[1], [2., 4., 6.], [8., 10., 12.],
                    [14., 16., 18.])
