"""Row E50 as specified: the builder's off-by-default --copy-rows (experiments/builds/build_test.py), on a fake board
and synthetic stages only. No board file, released stage or test count is read."""

from __future__ import annotations

import hashlib
import json
import re
import sys
from pathlib import Path
from types import SimpleNamespace

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


@pytest.mark.parametrize("expression", ["mixture", "series"])
def test_the_builder_path_calls_no_cloud_guard(bt, monkeypatch, world, own_repo, tmp_path, expression):
    """Row E50, the analyst's decision: the builder's copy-rows path has no cloud guard (the coordinates are checked
    bit for bit instead); a patched guard that raises is never reached and the file is written."""
    def never(*a, **kw):
        raise AssertionError("the cloud guard ran on the copy-rows path")

    _fake(bt, monkeypatch, world, _board(world[3]))
    monkeypatch.setattr(bt.t2_build, "cloud_guard", never)
    out = tmp_path / OUT
    bt.main([*BASE, "--expression", expression, "--copy-rows", "--copy-rows-test-ok", "--out", str(out)])
    assert out.exists()
    np.testing.assert_array_equal(_bits(anndata.read_h5ad(out).obsm[COORD_KEY]), _bits(_floor(bt)))


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


# ------------------------------------------------------------------------------------------------ row E52
# --copy-rows-size (off by default; the copy-rows file at the team's pre-registered size, the captain's order of
# 2026-11-21, decision log 2026-11-21 (v)). Synthetic stages and a synthetic size JSON only; the size file is not read.
SIZE_ARGS = ["--copy-rows", "--copy-rows-test-ok", "--copy-rows-size", "--copy-rows-size-test-ok"]
# The flag-off copy-rows file of the fixture world (mixture, the default seed), pinned on the builder before
# --copy-rows-size existed and unchanged after it: the flag-off path keeps its bytes and its record.
GOLDEN_CR_SHA256 = "3e126d6352eaaaf606fa6e75e66ada46e612381bcbdfcdaa9090fec18e6e048b"
GOLDEN_CR_KEYS = ["allow_below_truth", "board", "copy_rows", "copy_rows_copy", "count_decision", "count_rule",
                  "dry_run", "entry", "n_cells", "record", "seed", "size_rule", "target_rms"]


def test_the_flag_off_copy_rows_file_keeps_its_bytes_and_record(bt, monkeypatch, world, own_repo, tmp_path):
    _fake(bt, monkeypatch, world, _board(world[3]))
    out = tmp_path / OUT
    bt.main([*BASE, "--copy-rows", "--copy-rows-test-ok", "--out", str(out)])
    assert hashlib.sha256(out.read_bytes()).hexdigest() == GOLDEN_CR_SHA256
    rec = json.loads(out.with_suffix(".build.json").read_text())
    assert sorted(rec) == GOLDEN_CR_KEYS
    assert rec["size_rule"] == "not applied" and rec["target_rms"] is None and "copy_rows_size" not in rec
    assert rec["n_cells"] == "floor_copy" and rec["record"]["expression"] == "mixture"


# The exact-position file P's command (--method fgw_shape0 --positions floor_exact --expression mixture) on the fixture
# world with a synthetic size JSON (the bracket's own RMS) and a dry board, pinned with the builder of main before
# --copy-rows-size existed (b1430aa, run from a copy of that file) and unchanged after it. The cloud guard, a check that
# changes nothing it passes, is stubbed: this small synthetic world trips it.
GOLDEN_P_SHA256 = "ca775cc85b5f3b3efd9d5fe7376369b889d1819ea948b1f38ddf6783f36f98ff"


def test_the_exact_position_file_keeps_its_bytes(bt, monkeypatch, world, own_repo, tmp_path):
    board = _board(world[3], key="T2:embryo:val_interp", split="val")
    _fake(bt, monkeypatch, world, board, dry=True)
    monkeypatch.setattr(bt.t2_build, "cloud_guard", lambda *a, **kw: SimpleNamespace(ok=True, flags=[]))
    own = tb.rms_radius(np.asarray(world[1].obsm[COORD_KEY], dtype=np.float64)[:, :3])
    p = tmp_path / "size_embryo.json"
    p.write_text(json.dumps({"board": "embryo", "valid": True,
                             "preregistered": {"z_mode": "raw", "valid": True, "pred_rms_full": own}}))
    monkeypatch.setitem(bt.SIZES, "embryo", p)
    out = tmp_path / "T2embryo_s2x.h5ad"
    bt.main([*BASE, "--expression", "mixture", "--out", str(out)])
    assert hashlib.sha256(out.read_bytes()).hexdigest() == GOLDEN_P_SHA256


def test_the_size_flags_are_off_by_default(bt):
    a = bt.build_parser().parse_args([*BASE, "--out", f"runs/x/{OUT}"])
    assert a.copy_rows_size is False and a.copy_rows_size_test_ok is False
    a = bt.build_parser().parse_args([*BASE, "--copy-rows", "--out", f"runs/x/{OUT}"])
    assert a.copy_rows_size is False and a.copy_rows_size_test_ok is False
    assert bt._copy_rows_refusals(a) is None


def _no_loader(bt, monkeypatch):
    for name in ("_board", "load_stage", "read_stage", "load_boards", "_t2_count", "_copy_rows_size_target"):
        monkeypatch.setattr(bt, name, _raise)


@pytest.mark.parametrize("argv, match", [
    (["--copy-rows-size"], "--copy-rows-size is for --copy-rows only"),
    (["--copy-rows-size", "--copy-rows-size-test-ok"], "--copy-rows-size is for --copy-rows only"),
    (["--copy-rows", "--copy-rows-size-test-ok"], "--copy-rows-size-test-ok is for --copy-rows-size only"),
    (["--copy-rows-size-test-ok"], "--copy-rows-size-test-ok is for --copy-rows-size only"),
    (["--copy-rows", "--copy-rows-test-ok", "--copy-rows-size", "--copy-rows-size-test-ok", "--expression", "series"],
     "--copy-rows-size is for --expression mixture only"),
    (["--copy-rows", "--copy-rows-test-ok", "--copy-rows-size", "--copy-rows-size-test-ok", "--expression", "blend"],
     "--copy-rows-size is for --expression mixture only"),
])
def test_the_size_refusals_come_before_anything_loads_and_show_no_digit(bt, monkeypatch, argv, match):
    _no_loader(bt, monkeypatch)
    a = bt.build_parser().parse_args([*BASE, *argv, "--out", f"x/{OUT}"])
    with pytest.raises(SystemExit, match=re.escape(match)) as e:
        bt.build_embryo(a)
    assert not re.search(r"\d", str(e.value))


@pytest.mark.parametrize("drop, match", [("--copy-rows-size-test-ok", "--copy-rows-size-test-ok"),
                                         ("--copy-rows-test-ok", "--copy-rows-test-ok")])
@pytest.mark.parametrize("dry", [False, True])
def test_a_test_board_needs_both_yes_flags_before_anything_loads(bt, monkeypatch, world, drop, match, dry):
    board = _board(world[3]) if not dry else _board(world[3], key="T2:embryo:val_interp")  # split test either way
    _no_loader(bt, monkeypatch)
    monkeypatch.setattr(bt, "_board", lambda which, override: (board, dry))
    a = bt.build_parser().parse_args([*BASE, *[x for x in SIZE_ARGS if x != drop], "--out", f"x/{OUT}"])
    with pytest.raises(SystemExit) as e:
        bt.build_embryo(a)
    assert "needs the captain's yes" in str(e.value) and match in str(e.value)
    assert not re.search(r"\d", str(e.value))


def _size_json(tmp_path, target=3.25, valid=True, z_mode="raw"):
    p = tmp_path / "size_embryo_E7.75.json"
    p.write_text(json.dumps({"board": "embryo", "valid": valid,
                             "preregistered": {"z_mode": z_mode, "valid": valid, "pred_rms_full": target}}))
    return p


@pytest.mark.parametrize("kw", [{"valid": False}, {"z_mode": "xy_grow_z_hold"}, {"z_mode": "common_step"}])
def test_an_invalid_or_not_raw_size_json_is_refused_with_build_embryos_message(bt, monkeypatch, world, own_repo,
                                                                               tmp_path, kw):
    _fake(bt, monkeypatch, world, _board(world[3]))
    monkeypatch.setitem(bt.SIZES, "embryo", _size_json(tmp_path, **kw))
    out = tmp_path / OUT
    with pytest.raises(SystemExit, match="the embryo size JSON is not a valid raw-mode result"):
        bt.main([*BASE, *SIZE_ARGS, "--out", str(out)])
    assert not out.exists()


@pytest.mark.parametrize("target", [0.0, -1.0, float("inf"), float("nan")])
def test_a_target_that_is_not_finite_and_positive_is_refused_with_no_value(bt, monkeypatch, world, own_repo, tmp_path,
                                                                           target):
    _fake(bt, monkeypatch, world, _board(world[3]))
    monkeypatch.setitem(bt.SIZES, "embryo", _size_json(tmp_path, target=target))
    with pytest.raises(SystemExit) as e:
        bt.main([*BASE, *SIZE_ARGS, "--out", str(tmp_path / OUT)])
    assert "finite and positive" in str(e.value) and not re.search(r"\d", str(e.value))


@pytest.fixture
def sized_build(bt, monkeypatch, world, own_repo, tmp_path):
    """The flag-off copy-rows file and the sized one of the fixture world, at a synthetic pre-registered size."""
    _fake(bt, monkeypatch, world, _board(world[3]))
    size = _size_json(tmp_path, target=3.25)
    monkeypatch.setitem(bt.SIZES, "embryo", size)
    cr, cs = tmp_path / "cr" / OUT, tmp_path / "cs" / "T2embryo_crs.h5ad"
    bt.main([*BASE, "--copy-rows", "--copy-rows-test-ok", "--out", str(cr)])
    bt.main([*BASE, *SIZE_ARGS, "--out", str(cs)])
    return cr, cs, size


def test_the_sized_file_is_the_copy_rows_file_rescaled_and_its_record(bt, sized_build):
    cr, cs, size = sized_build
    assert hashlib.sha256(cr.read_bytes()).hexdigest() == GOLDEN_CR_SHA256  # the flag-off file, unchanged
    A, B = anndata.read_h5ad(cr), anndata.read_h5ad(cs)
    XA = A.X.toarray() if hasattr(A.X, "toarray") else np.asarray(A.X)
    XB = B.X.toarray() if hasattr(B.X, "toarray") else np.asarray(B.X)
    assert XA.dtype == XB.dtype and XA.tobytes() == XB.tobytes()  # C's expression, byte for byte
    assert list(A.obs_names) == list(B.obs_names) and list(A.var_names) == list(B.var_names)
    floor = _floor(bt)
    C64 = floor.astype(np.float64)
    mu = C64.mean(0)
    want = (mu + (C64 - mu) * (3.25 / tb.rms_radius(C64))).astype(np.float32)
    np.testing.assert_array_equal(_bits(B.obsm[COORD_KEY]), _bits(want))
    np.testing.assert_array_equal(_bits(B.obsm[COORD_KEY]), _bits(tb.copy_rows_sized(floor, 3.25)))
    assert abs(tb.rms_radius(np.asarray(B.obsm[COORD_KEY], np.float64)) / 3.25 - 1) <= 1e-6
    rec = json.loads(cs.with_suffix(".build.json").read_text())
    assert rec["copy_rows_size"] is True and rec["target_rms"] == 3.25
    assert rec["size_rule"] == ("applied: rescale_to_rms about the copy centroid to preregistered pred_rms_full "
                                "(raw mode)") == bt.COPY_ROWS_SIZE_RULE
    assert rec["size_sha256"] == hashlib.sha256(size.read_bytes()).hexdigest() and rec["size_file"].endswith(size.name)
    assert rec["n_cells"] == "floor_copy" and bt.COPY_ROWS_N not in rec and len(floor) not in _ints(rec)
    c_rec = json.loads(cr.with_suffix(".build.json").read_text())
    assert sorted(rec) == sorted([*c_rec, "copy_rows_size", "size_file", "size_sha256"])
    assert {k: v for k, v in rec.items() if k not in ("copy_rows_size", "size_rule", "target_rms", "size_file",
                                                     "size_sha256")} == {k: v for k, v in c_rec.items()
                                                                         if k not in ("size_rule", "target_rms")}


def test_the_parity_self_check_is_count_free_and_catches_a_mirror(bt, monkeypatch, world):
    _fake(bt, monkeypatch, world, _board(world[3]))
    floor = _floor(bt)
    S32 = tb.copy_rows_sized(floor, 7.5).astype(np.float32)
    assert bt._copy_rows_parity_ok(S32, floor)
    assert not bt._copy_rows_parity_ok(ev_mirror(S32), floor)
    assert not bt._copy_rows_parity_ok(S32[:-1], floor)
    assert bt._copy_rows_rms_ok(S32, 7.5) and not bt._copy_rows_rms_ok(S32, 7.5 * (1 + 3e-6))


def ev_mirror(C):
    M = np.array(C, copy=True)
    M[:, 2] = -M[:, 2]
    return M


@pytest.mark.parametrize("check", ["_copy_rows_rms_ok", "_copy_rows_parity_ok"])
def test_a_failing_self_check_refuses_before_writing_with_no_number(bt, monkeypatch, world, own_repo, tmp_path,
                                                                    check):
    _fake(bt, monkeypatch, world, _board(world[3]))
    monkeypatch.setitem(bt.SIZES, "embryo", _size_json(tmp_path))
    monkeypatch.setattr(bt, check, lambda *a, **kw: False)
    out = tmp_path / OUT
    with pytest.raises(SystemExit) as e:
        bt.main([*BASE, *SIZE_ARGS, "--out", str(out)])
    assert "--copy-rows-size" in str(e.value) and not re.search(r"\d", str(e.value))
    assert not out.exists()


def test_a_failing_read_back_removes_the_file_and_its_stamp(bt, monkeypatch, world, own_repo, tmp_path):
    from vec import provenance

    _fake(bt, monkeypatch, world, _board(world[3]))
    monkeypatch.setitem(bt.SIZES, "embryo", _size_json(tmp_path))
    real = anndata.read_h5ad

    def read_back(path, *a, **kw):
        x = real(path, *a, **kw)
        c = np.asarray(x.obsm[COORD_KEY], np.float32).copy()
        c[1, 2] = np.nextafter(c[1, 2], np.float32(np.inf))
        x.obsm[COORD_KEY] = c
        return x

    monkeypatch.setattr(anndata, "read_h5ad", read_back)
    out = tmp_path / OUT
    with pytest.raises(SystemExit) as e:
        bt.main([*BASE, *SIZE_ARGS, "--out", str(out)])
    assert "rescaled copy's" in str(e.value) and not re.search(r"\d", str(e.value))
    assert not out.exists() and not provenance.stamp_path(out).exists()


def test_the_sized_record_hides_counts_like_the_copy_rows_record(bt):
    board = _board(["g0"], max_cells=777)
    rec = {"n_cells": 69, bt.COPY_ROWS_N: 69, "copy_rows_size": True, "target_rms": 3.25,
           "record": {"n": 69, "n_before": 35, "n_after": 34, "n_jittered": 0, "t": 0.5}}
    out = bt._hide_counts(rec, board)
    assert out["n_cells"] == "floor_copy" and bt.COPY_ROWS_N not in out
    assert out["record"] == {"t": 0.5} and not set(bt.COUNT_FIELDS) & set(out["record"])
    with pytest.raises(SystemExit):
        bt._hide_counts({**rec, "n_cells": 70}, board)
