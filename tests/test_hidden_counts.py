"""s1 build review (2026-11-01): no tool prints a test board's cell limits or a test file's cell count, and the
board's required keys are checked. A test board's limits can encode the held-out stage's cell count, so they are read
in code only (red-team review 2026-10-31)."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pytest

from vec.boards import Board
from vec.io import check_submission, make_submission, write_submission

GENES = [f"G{i}" for i in range(20)]
SECRET_MAX = 4321  # stands for a test board's max_cells: never printed


def _file(tmp_path, n, coords=True, name="pred.h5ad"):
    rng = np.random.default_rng(0)
    X = np.log1p(rng.poisson(2.0, (n, len(GENES))).astype(np.float32))
    a = make_submission(X, GENES, rng.normal(size=(n, 3)).astype(np.float32) if coords else None)
    path = tmp_path / name
    write_submission(a, path)
    return path


def _board(key="T2:heart:test_extrap", split="test_extrap", **kw):
    return Board(key=key, task="T2", genes=GENES, needs_coords=True, split=split, setting="heart", min_cells=10,
                 max_cells=SECRET_MAX, **kw)


def test_check_hides_the_limits_of_a_test_board(tmp_path):
    rep = check_submission(_file(tmp_path, 5), _board())
    assert not rep.ok and any("outside the board's limits" in e for e in rep.errors)
    assert not any(str(SECRET_MAX) in e or " 5 cells" in e for e in rep.errors)
    shown = check_submission(_file(tmp_path, 5), _board(key="T2:heart:val_extrap", split="val_extrap"))
    assert any(str(SECRET_MAX) in e for e in shown.errors)  # validation boards keep their numbers


def test_check_enforces_the_board_required_keys(tmp_path):
    path = _file(tmp_path, 50)
    ok = check_submission(path, _board(obsm_required=["spatial_3D"]))
    assert ok.ok, ok.errors
    rep = check_submission(path, _board(obs_required=["sample"], obsm_required=["spatial_3D", "X_other"]))
    assert "obs['sample'] is required by the board spec" in rep.errors
    assert "obsm['X_other'] is required by the board spec" in rep.errors


def test_vec_check_prints_max_cells_not_the_count(tmp_path, monkeypatch):
    from typer.testing import CliRunner

    from vec import cli

    path = _file(tmp_path, SECRET_MAX)
    monkeypatch.setattr(cli, "get_board", lambda key: _board())
    res = CliRunner().invoke(cli.app, ["check", str(path), "--board", "T2:heart:test_extrap"])
    assert res.exit_code == 0, res.output
    assert "max_cells cells" in res.output and str(SECRET_MAX) not in res.output


def test_vec_boards_hides_test_board_limits(tmp_path, monkeypatch):
    from typer.testing import CliRunner

    from vec import cli

    boards = {"T2:heart:test_extrap": _board(), "T2:heart:val_extrap": _board(key="T2:heart:val_extrap",
                                                                              split="val_extrap")}
    monkeypatch.setattr(cli, "load_boards", lambda: boards)
    res = CliRunner().invoke(cli.app, ["boards"])
    assert res.exit_code == 0, res.output
    assert res.output.count(str(SECRET_MAX)) == 1 and "hidden" in res.output


# ---------------------------------------------------------------------------------------------------------------------
# The T3 builder: seed and count (experiments/t3_decomp/build_submission.py)

T3_DIR = Path(__file__).resolve().parents[1] / "experiments" / "t3_decomp"


@pytest.fixture(scope="module")
def t3b():
    sys.path.insert(0, str(T3_DIR))
    import build_submission

    return build_submission


def test_t3_count_is_read_from_the_board(t3b):
    b = Board(key="T3:bcatenin", task="T3", genes=GENES, needs_coords=True, split="test", max_cells=SECRET_MAX)
    assert t3b.resolve_n_cells("max", b) == SECRET_MAX
    assert t3b.count_label(SECRET_MAX, b) == "max_cells" and t3b.count_label(7, b) == "not max_cells"
    with pytest.raises(t3b.Refused, match="no max_cells"):
        t3b.resolve_n_cells("max", Board(key="T3:bcatenin", task="T3", genes=GENES, needs_coords=True))
    with pytest.raises(t3b.Refused, match="'max' or an integer"):
        t3b.resolve_n_cells("lots", b)
    syn = Board(key="T3:synthetic", task="T3", genes=GENES, needs_coords=True, max_cells=3000)
    assert t3b.resolve_n_cells("250", syn) == 250 and t3b.count_label(250, syn) == 250


@pytest.mark.parametrize("seed", [0, 1, 999_999])
def test_t3_build_refuses_a_seed_below_the_offset(t3b, tmp_path, seed):
    b = Board(key="T3:synthetic", task="T3", genes=GENES, needs_coords=True, max_cells=3000)
    with pytest.raises(t3b.Refused, match="below 1000000"):
        t3b.build_submission(tmp_path / "d.json", "allowed", "s1", tmp_path / "wt.h5ad", b, 100,
                             tmp_path / "out.h5ad", seed=seed)


def test_t3_cli_defaults_are_the_offset_seed_and_max_cells(t3b):
    import inspect

    params = inspect.signature(t3b.build).parameters
    assert params["seed"].default.default == 1_000_001 and params["n_cells"].default.default == "max"
    assert inspect.signature(t3b.build_submission).parameters["seed"].default == 1_000_001


# ---------------------------------------------------------------------------------------------------------------------
# experiments/builds/build_test.py: dry runs are named, test builds are not, counts are hidden

BUILDS = Path(__file__).resolve().parents[1] / "experiments" / "builds"


@pytest.fixture(scope="module")
def bt():
    sys.path.insert(0, str(BUILDS))
    import build_test

    return build_test


@pytest.mark.parametrize("argv, msg", [
    (["t1", "--board", "T1:val", "--out", "runs/x/T1_s2.h5ad"], "must contain 'dry'"),
    (["t1", "--board", "T2:heart:val_extrap", "--out", "runs/x/T1_dry.h5ad"], "takes only"),
    (["t1", "--out", "runs/x/T1_dry.h5ad"], "must not contain 'dry'"),
    (["heart", "--entry", "s1", "--seed", "5", "--out", "runs/x/h.h5ad"], "below 1000000"),
    (["heart", "--entry", "s2_trim", "--seed", "5", "--out", "runs/x/h.h5ad"], "below 1000000"),
])
def test_build_test_refuses_misnamed_dry_runs_and_low_seeds(bt, argv, msg):
    with pytest.raises(SystemExit) as e:
        bt.main(argv)
    assert msg in str(e.value)


def test_build_record_hides_counts_on_a_test_board(bt):
    b = _board()
    rec = {"board": b.key, "n_cells": SECRET_MAX, "record": {"n": SECRET_MAX, "n_before": 2000, "n_after": 2321,
                                                              "n_jittered": 3, "t": 0.5}}
    out = bt._hide_counts(rec, b)
    assert out["n_cells"] == "max_cells" and out["record"] == {"t": 0.5}
    assert str(SECRET_MAX) not in json.dumps(out)
    with pytest.raises(SystemExit, match="does not hold the board's max_cells"):
        bt._hide_counts(dict(rec, n_cells=7), b)
    val = _board(key="T2:heart:val_extrap", split="val_extrap")
    assert bt._hide_counts(rec, val) is rec
