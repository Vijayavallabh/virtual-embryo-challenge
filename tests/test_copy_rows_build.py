"""Row E50 as specified: the builder's off-by-default --copy-rows (experiments/builds/build_test.py), on a fake board
and synthetic stages only. No board file, released stage or test count is read."""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import anndata
import numpy as np
import pytest
from test_t2_build import _labelled_pair

from vec import t2_build as tb
from vec.boards import Board
from vec.io import COORD_KEY

REPO = Path(__file__).resolve().parents[1]
OUT = "T2embryo_cr.h5ad"
BASE = ["embryo", "--method", "fgw_shape0", "--positions", "floor_exact"]


@pytest.fixture(scope="module")
def bt():
    sys.path.insert(0, str(REPO / "experiments" / "builds"))
    import build_test

    return build_test


@pytest.fixture(scope="module")
def world():
    """A released earlier bracket with exact copies inserted (every 7th cell repeated), the bracket after the
    duplicate filter (uns["qc"]["dropped_rows"]), a later bracket and the panel."""
    before, after, panel = _labelled_pair(np.random.default_rng(61))
    order, dropped = [], []
    for i in range(before.n_obs):
        order.append(i)
        if i % 7 == 3:
            dropped.append(len(order))
            order.append(i)
    released = before[order].copy()
    released.obs_names_make_unique()
    deduped = released[np.setdiff1d(np.arange(len(order)), dropped)].copy()
    deduped.uns["qc"] = {"drop_duplicates": True, "dropped_rows": np.asarray(dropped, np.int64)}
    return released, deduped, after, list(panel)


def _board(panel, **kw):
    base = dict(key="T2:embryo:test_interp", task="T2", genes=panel, needs_coords=True, setting="embryo",
                split="test", min_cells=10, max_cells=5000)
    return Board(**{**base, **kw})


def _raise(*a, **kw):
    raise AssertionError("a loader ran")


def _fake(bt, monkeypatch, world, board, dry=False):
    released, deduped, after, panel = world

    def load(path, genes=None, **kw):
        return deduped if Path(path).name == "E7.5.h5ad" else after

    def read(path, genes=None):
        assert Path(path).name == "E7.5.h5ad"
        return released if genes is None else released[:, list(genes)].copy()

    monkeypatch.setattr(bt, "_board", lambda which, override: (board, dry))
    monkeypatch.setattr(bt, "load_boards", lambda: {board.key: board})
    monkeypatch.setattr(bt, "load_stage", load)
    monkeypatch.setattr(bt, "read_stage", read)
    monkeypatch.setattr(bt, "_t2_count", lambda board, n: {"decision": "fake"})


def _floor(bt):
    return bt.floor_copy_coords("T2/E7.5.h5ad")


def _bits(C):
    return np.ascontiguousarray(np.asarray(C, dtype=np.float32)).view(np.uint32)


def _ints(obj):
    """Every integer anywhere in a record (bools excepted)."""
    if isinstance(obj, dict):
        return [x for v in obj.values() for x in _ints(v)]
    if isinstance(obj, list):
        return [x for v in obj for x in _ints(v)]
    return [obj] if isinstance(obj, int) and not isinstance(obj, bool) else []


# ------------------------------------------------------------------------------------------------ off by default
def test_off_by_default(bt):
    a = bt.build_parser().parse_args([*BASE, "--out", f"runs/x/{OUT}"])
    assert a.copy_rows is False and a.copy_rows_test_ok is False
    assert bt._copy_rows_refusals(a) is None
    board = _board(["g0"], max_cells=777)
    rec = {"n_cells": 777, "record": {"n": 777, "t": 0.5}, "entry": "fgw_shape0"}
    assert bt._hide_counts(rec, board) == {"n_cells": "max_cells", "record": {"t": 0.5}, "entry": "fgw_shape0"}
    val = Board(key="T2:embryo:val_interp", task="T2", genes=["g0"], needs_coords=True, split="val")
    assert bt._hide_counts(rec, val) == rec


@pytest.mark.parametrize("extra, match", [
    (["--row-order", "copy_parity"], "--copy-rows is refused"),
    (["--shape-interp", "overlay"], "--copy-rows is refused"),
    (["--herd"], "--copy-rows is refused"),
    (["--draw", "midpoint"], "--copy-rows is refused"),
    (["--move-coupling", "aligned"], "--copy-rows is refused"),
    (["--expression", "blend"], "--copy-rows is refused"),
])
def test_incompatible_flags_are_refused_before_anything_loads(bt, monkeypatch, extra, match):
    for name in ("_board", "load_stage", "read_stage", "load_boards"):
        monkeypatch.setattr(bt, name, _raise)
    a = bt.build_parser().parse_args([*BASE, "--copy-rows", "--copy-rows-test-ok", *extra, "--out", f"x/{OUT}"])
    with pytest.raises(SystemExit, match=match):
        bt.build_embryo(a)


def test_other_methods_and_the_ok_flag_alone_are_refused(bt, monkeypatch):
    for name in ("_board", "load_stage", "read_stage", "load_boards"):
        monkeypatch.setattr(bt, name, _raise)
    ap = bt.build_parser()
    with pytest.raises(SystemExit, match="--copy-rows-test-ok is for --copy-rows only"):
        bt.build_embryo(ap.parse_args([*BASE, "--copy-rows-test-ok", "--out", f"x/{OUT}"]))
    for argv in (["embryo", "--method", "fgw_shape0", "--copy-rows"],
                 ["embryo", "--method", "fgw", "--positions", "floor_exact", "--copy-rows"]):
        with pytest.raises(SystemExit, match="fgw_shape0 --positions floor_exact only"):
            bt.build_embryo(ap.parse_args([*argv, "--out", f"x/{OUT}"]))


def test_a_test_board_needs_the_ok_flag_before_anything_loads(bt, monkeypatch, world):
    board = _board(world[3])
    monkeypatch.setattr(bt, "_board", lambda which, override: (board, False))
    for name in ("load_stage", "read_stage", "_t2_count"):
        monkeypatch.setattr(bt, name, _raise)
    a = bt.build_parser().parse_args([*BASE, "--copy-rows", "--out", f"x/{OUT}"])
    with pytest.raises(SystemExit, match="--copy-rows-test-ok"):
        bt.build_embryo(a)
    # a dry-run board whose split says test is a test board too
    monkeypatch.setattr(bt, "_board", lambda which, override: (_board(world[3], key="T2:embryo:val_interp"), True))
    with pytest.raises(SystemExit, match="--copy-rows-test-ok"):
        bt.build_embryo(a)


# ------------------------------------------------------------------------------------------------ the build
@pytest.mark.parametrize("expression", ["mixture", "series"])
def test_the_test_build_is_the_floor_copy_and_its_record_holds_no_count(bt, monkeypatch, world, own_repo, tmp_path,
                                                                        capsys, expression):
    released, deduped, after, panel = world
    board = _board(panel)
    _fake(bt, monkeypatch, world, board)
    out = tmp_path / OUT
    bt.main([*BASE, "--expression", expression, "--copy-rows", "--copy-rows-test-ok", "--out", str(out)])
    a = anndata.read_h5ad(out)
    floor = _floor(bt)
    np.testing.assert_array_equal(_bits(a.obsm[COORD_KEY]), _bits(floor))
    rec = json.loads(out.with_suffix(".build.json").read_text())
    assert rec["n_cells"] == "floor_copy" and rec["copy_rows"] is True and rec["allow_below_truth"] is True
    assert rec["count_rule"] == bt.COPY_ROWS_COUNT_RULE and rec["size_rule"] == "not applied"
    assert rec["copy_rows_copy"] == "vec.evaluate.subsample(T2/E7.5.h5ad as released, 0.1, 0)"
    assert bt.COPY_ROWS_N not in rec and rec["record"]["expression"] == expression
    assert len(floor) not in _ints(rec)
    printed = json.loads(capsys.readouterr().out.strip().splitlines()[-1])
    assert printed["n_cells"] == "floor_copy" and len(floor) not in _ints(printed)
    if expression == "mixture":  # half the slots keep their own released row's expression
        Xr = np.asarray(released[:, panel].X.toarray(), np.float32)
        idx = np.sort(np.random.default_rng(0).choice(released.n_obs, len(floor), replace=False))
        Xf = a.X.toarray() if hasattr(a.X, "toarray") else np.asarray(a.X)
        own = sum(np.array_equal(Xf[i], Xr[j]) for i, j in enumerate(idx))
        assert own == int(round(0.5 * len(floor)))


def test_the_own_limits_check_and_the_wrapped_refusal_show_no_value(bt, monkeypatch, world, own_repo, tmp_path):
    panel = world[3]
    out = tmp_path / OUT
    _fake(bt, monkeypatch, world, _board(panel, min_cells=4000, max_cells=5000))  # limits that exclude the copy
    with pytest.raises(SystemExit) as e:
        bt.main([*BASE, "--copy-rows", "--copy-rows-test-ok", "--out", str(out)])
    assert "outside the board's limits" in str(e.value) and not re.search(r"\d", str(e.value))
    assert not out.exists()
    _fake(bt, monkeypatch, world, _board(panel))

    def refuse(*a, **kw):
        raise ValueError("T2:embryo:test_interp: 69 cells outside [10, 5000]")

    monkeypatch.setattr(bt.t2_build, "build_submission", refuse)
    with pytest.raises(SystemExit) as e:
        bt.main([*BASE, "--copy-rows", "--copy-rows-test-ok", "--out", str(out)])
    assert "build_submission refused the copy-rows file" in str(e.value) and not re.search(r"\d", str(e.value))


def test_the_bit_self_checks_fire_on_one_ulp(bt, monkeypatch, world, own_repo, tmp_path):
    panel = world[3]
    out = tmp_path / OUT
    _fake(bt, monkeypatch, world, _board(panel))
    orig = tb.copy_rows_cloud

    def nudged(*a, **kw):
        res = orig(*a, **kw)
        c = res.coords.astype(np.float32)
        c[3, 1] = np.nextafter(c[3, 1], np.float32(np.inf))
        res.coords = c.astype(np.float64)
        return res

    monkeypatch.setattr(bt.t2_build, "copy_rows_cloud", nudged)
    with pytest.raises(SystemExit, match="bit for bit"):
        bt.main([*BASE, "--copy-rows", "--copy-rows-test-ok", "--out", str(out)])
    assert not out.exists()
    monkeypatch.setattr(bt.t2_build, "copy_rows_cloud", orig)
    real = anndata.read_h5ad

    def read_back(path, *a, **kw):
        x = real(path, *a, **kw)
        c = np.asarray(x.obsm[COORD_KEY], np.float32).copy()
        c[0, 0] = np.nextafter(c[0, 0], np.float32(-np.inf))
        x.obsm[COORD_KEY] = c
        return x

    monkeypatch.setattr(anndata, "read_h5ad", read_back)
    with pytest.raises(SystemExit, match="written file"):
        bt.main([*BASE, "--copy-rows", "--copy-rows-test-ok", "--out", str(out)])
    assert not out.exists()


def test_the_launcher_list_names_this_file():
    sh = (REPO / "experiments" / "external" / "n54_step.sh").read_text()
    assert "tests/test_copy_rows_build.py" in sh
