"""The occupancy frame check (experiments/coords/occupancy_frame_check.py) and the off-by-default row order
`--row-order copy_parity` (vec.t2_build.copy_parity_order, experiments/builds/build_test.py embryo), on synthetic
clouds, stages and a synthetic test board whose limits are never shown. Nothing is read from data/ and nothing is
scored: the scorer's functions are called on synthetic arrays only.

What is pinned: the check reports agreement for a cloud and its plain copy and disagreement for a mirrored cloud (both
scorer builds) and for a re-ordered cloud whose pre-fix parity differs (the pre-fix build; the main build cannot tell
them apart, by design); its overlap is the main build's occupancy_dice; the row order restores the copy's pre-fix parity
and keeps every cell's expression with its position, the pairwise distances and the scale; the builder's default path
writes the bytes it wrote before the option existed."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import re
import sys
from functools import partial
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse
from scipy.spatial.distance import pdist

from vec import evaluate as ev
from vec import t2_build as tb
from vec.boards import Board

REPO = Path(__file__).resolve().parents[1]
GENES = [f"g{i}" for i in range(12)]
SECRET_MAX = 400  # stands for the test board's max_cells: never printed, never in a record
SECRET_TRUTH = 120  # stands for its truth_cells
N_BEFORE, N_AFTER = 1500, 900
SEED = 1_000_001
# sha256 of the .h5ad `build_test.py embryo --method fgw_shape0 [--positions floor_exact]` wrote on `embryo_env` at the
# default seed, computed with the builder of main at d569eb3 (before --row-order was added)
PRE_OPTION_SHA256 = {"moved": "b05224fda4f1a72b40a3bb39c6de8de793b9e2b2f324c5aa27077500922273d2",
                     "floor_exact": "e01377a03d71dc462366bf939867d695ecc68a7e72139c29deacfcd6e8e4c5de"}


@pytest.fixture(scope="module")
def chk():
    spec = importlib.util.spec_from_file_location("occupancy_frame_check",
                                                  REPO / "experiments" / "coords" / "occupancy_frame_check.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def _cloud(n: int, seed: int) -> np.ndarray:
    """An anisotropic, rotated cloud with no mirror symmetry (a Gaussian body with two small lumps placed so that no
    reflection maps the cloud onto itself or onto a proper sign flip of itself), float32 as a file holds it."""
    rng = np.random.default_rng(seed)
    Z = rng.normal(size=(n, 3)) * [3.0, 1.6, 0.7]
    k1, k2 = n // 5, n // 8
    Z[:k1] = [3.5, 2.2, 1.0] + rng.normal(size=(k1, 3)) * 0.4
    Z[k1:k1 + k2] = [-2.5, 1.6, -0.9] + rng.normal(size=(k2, 3)) * 0.4
    Z = Z[rng.permutation(n)]
    th = 0.7
    R = np.array([[np.cos(th), -np.sin(th), 0.0], [np.sin(th), np.cos(th), 0.0], [0.0, 0.0, 1.0]])
    return (Z @ R.T + [10.0, -4.0, 2.0]).astype(np.float32)


def _parity_flip_order(C: np.ndarray, want_other_than: int, n_true: int, seed: int = 0) -> np.ndarray:
    """A row order of C whose pre-fix parity (n_true truth cells, metric seed `seed`) is not `want_other_than`."""
    rng = np.random.default_rng(99)
    for _ in range(200):
        o = rng.permutation(len(C))
        if ev.prefix_parity(C[o], n_true, seed) != want_other_than:
            return o
    raise AssertionError("no order flips the parity")


# ------------------------------------------------------------------------------------------------ the check
def test_a_cloud_and_its_plain_copy_agree(chk):
    stage = _cloud(4000, 1)
    copy = stage[np.sort(np.random.default_rng(5).choice(len(stage), 400, replace=False))]
    refs = [stage[np.sort(np.random.default_rng(50 + k).choice(len(stage), 400, replace=False))] for k in range(5)]
    f = chk.frame_facts(copy.copy(), copy, refs=refs)
    assert f["axes_agree"] and f["handedness_agree"] and f["coincide"]
    assert f["overlap"] == 1.0 and f["axis_cos"] == [1.0, 1.0, 1.0]
    # another draw of the same stage (a file on the stage's own positions) agrees too
    other = stage[np.sort(np.random.default_rng(6).choice(len(stage), 1200, replace=False))]
    g = chk.frame_facts(other, copy, refs=refs)
    assert g["axes_agree"] and g["handedness_agree"] and g["coincide"]
    assert set(g) == set(chk.MAIN_FIELDS)
    # translation and an isotropic rescale (the size rule) change nothing in the scorer's frame
    h = chk.frame_facts(((other - other.mean(0)) * 1.37 + 5.0).astype(np.float32), copy, refs=refs)
    assert h == g


def test_a_mirrored_cloud_disagrees_under_both_builds(chk):
    stage = _cloud(4000, 2)
    copy = stage[np.sort(np.random.default_rng(5).choice(len(stage), 400, replace=False))]
    refs = [stage[np.sort(np.random.default_rng(50 + k).choice(len(stage), 400, replace=False))] for k in range(5)]
    m = ev.mirror_z(copy)
    f = chk.frame_facts(m, copy, refs=refs)
    assert not f["handedness_agree"] and not f["coincide"] and f["overlap"] < f["overlap_mirror"]
    # under the pre-fix build a z-mirror has the opposite parity on the same rows, so it falls in the other mode
    p = chk.prefix_facts(m, copy, n_true=300)
    assert p["truth_rows_shared"] and not p["parity_agree"] and not p["same_mode_whatever_truth"]


def test_a_reordered_cloud_disagrees_only_where_the_convention_sees_order(chk):
    stage = _cloud(4000, 3)
    copy = stage[np.sort(np.random.default_rng(5).choice(len(stage), 400, replace=False))]
    refs = [stage[np.sort(np.random.default_rng(50 + k).choice(len(stage), 400, replace=False))] for k in range(5)]
    o = _parity_flip_order(copy, ev.prefix_parity(copy, 300), 300)
    re = copy[o]
    # main build: row order is absorbed (right-handed basis, proper flips), so the frames agree
    f = chk.frame_facts(re, copy, refs=refs)
    assert f["axes_agree"] and f["handedness_agree"] and f["coincide"] and f["overlap"] == 1.0
    # pre-fix build: the same rows in another order have the other parity, the other mode
    p = chk.prefix_facts(re, copy, n_true=300)
    assert p["truth_rows_shared"] and not p["parity_agree"]
    assert chk.prefix_facts(copy.copy(), copy, n_true=300)["parity_agree"]


def test_the_overlap_is_the_main_builds_occupancy_dice(chk):
    sm = ev._shape_metrics()
    stage = _cloud(3000, 4)
    copy = stage[np.sort(np.random.default_rng(1).choice(len(stage), 300, replace=False))]
    for P in (stage[:900], ev.mirror_z(stage[:900]), _cloud(700, 9)):
        with ev.scorer_build("main"):
            want = sm.occupancy_dice(P, copy, seed=0)[0]
            want_m = sm.occupancy_dice(ev.mirror_z(P), copy, seed=0)[0]
        f = chk.frame_facts(P, copy, refs=[])
        assert f["overlap"] == round(want, 4) and f["overlap_mirror"] == round(want_m, 4)
    # under the pre-fix build the scorer returns the proper or mirror value as the two parities say
    modes = set()
    for k in range(8):
        P = stage[:900][np.random.default_rng([3, k]).permutation(900)]
        with ev.scorer_build("prefix"):
            got = sm.occupancy_dice(P, copy, seed=0)[0]
        p = chk.prefix_facts(P, copy, n_true=10_000, overlap=True)
        assert p["overlap_prefix"] == round(got, 4)
        f = chk.frame_facts(P, copy, refs=[])
        modes.add(p["overlap_prefix"] == f["overlap"])
    assert modes == {True, False}  # both the proper and the mirror mode were met


def test_file_mode_reads_a_file_and_its_stage_and_prints_main_facts_only(chk, tmp_path, capsys):
    stage = _cloud(5000, 13)
    (tmp_path / "T2").mkdir()
    a = ad.AnnData(X=np.zeros((len(stage), 1), dtype=np.float32))
    a.obsm["spatial_3D"] = stage
    a.write_h5ad(tmp_path / "T2" / "E7.5.h5ad")
    f = ad.AnnData(X=np.zeros((1700, 1), dtype=np.float32))
    f.obsm["spatial_3D"] = ((stage[np.random.default_rng(2).permutation(5000)[:1700]] - 3.0) * 0.8).astype(np.float32)
    f.write_h5ad(tmp_path / "T2embryo_dry.h5ad")
    chk.main(["file", "--file", str(tmp_path / "T2embryo_dry.h5ad"), "--data-dir", str(tmp_path), "--n-ref", "8"])
    out = capsys.readouterr().out
    row = json.loads(out)
    assert row["coincide"] and row["prefix"] == chk.WITHHELD and row["build"] == "main"
    assert set(row) == {"mode", "file", "copy", "metric_seed", "build", "prefix", *chk.MAIN_FIELDS}
    assert not _shows(out, 1700, 5000, 500, 170)


def test_practice_mode_loops_tasks_and_seeds_with_both_row_orders(chk, tmp_path, capsys, monkeypatch):
    """The practice loop on stand-ins for the harness (its task list, board, inputs and candidate): one row per task,
    seed and row order, the bake copy of the reference at the seed, the pre-fix facts at the task's truth count."""
    from types import SimpleNamespace

    from vec import boards as vb
    from vec import harness, prep
    from vec import io as vio

    stage = _cloud(6000, 21)
    ref = ad.AnnData(X=sparse.csr_matrix(np.ones((len(stage), len(GENES)), dtype=np.float32)),
                     var=pd.DataFrame(index=GENES))
    ref.obsm["spatial_3D"] = stage
    t = SimpleNamespace(key="T2e/x", task="T2", setting="embryo", kind="interp", panel_board="P", reference="T2/r",
                        target="T2/t", inputs={"before": "T2/r"})
    board = Board(key="pseudo:x", task="T2", genes=list(GENES), needs_coords=True, split="pseudo", min_cells=50,
                  max_cells=2000, truth_cells=450)
    built = []

    def build(ctx):
        built.append(ctx.seed)
        rows = np.random.default_rng(ctx.seed).permutation(len(stage))[:2000]
        return ev.Stage(X=np.zeros((2000, 0), np.float32), coords=stage[rows])

    monkeypatch.setattr(harness, "load_tasks", lambda: {"T2:embryo:test_interp": [t]})
    monkeypatch.setattr(harness, "task_board", lambda t_, b, n_target, n_ref: board)
    monkeypatch.setattr(harness.TaskView, "of", classmethod(lambda cls, t_: "view"))
    monkeypatch.setitem(harness.CANDIDATES, "fgw_shape0x", SimpleNamespace(applies=lambda t_: True, build=build))
    monkeypatch.setattr(vb, "load_boards", lambda: {"P": board})
    monkeypatch.setattr(vio, "read_stage", lambda p, genes=None: ref)
    monkeypatch.setattr(prep, "load_stage", lambda p, genes=None, drop_duplicates=False: ref)
    monkeypatch.setattr(chk, "_cells_only", lambda p: 9000)
    chk.main(["practice", "--seeds", "0", "1", "--n-ref", "6", "--data-dir", str(tmp_path)])
    rows = [json.loads(x) for x in capsys.readouterr().out.splitlines()]
    assert built == [harness.BUILD_SEED_OFFSET, harness.BUILD_SEED_OFFSET + 1]
    assert [(r["seed"], r["row_order"]) for r in rows] == [(0, "as_built"), (0, "copy_parity"), (1, "as_built"),
                                                           (1, "copy_parity")]
    for r in rows:
        assert r["coincide"] and set(r["prefix"]) == {"parity_file", "parity_copy", "parity_agree",
                                                      "truth_rows_shared", "same_mode_whatever_truth",
                                                      "overlap_prefix"}
        if r["row_order"] == "copy_parity":
            assert r["prefix"]["same_mode_whatever_truth"]


def test_the_check_prints_no_count_size_or_score(chk):
    stage = _cloud(4000, 7)
    copy = stage[:400]
    f = chk.frame_facts(stage[1000:2200], copy, refs=[stage[2200:2600]])
    text = json.dumps(f)
    for v in ("1200", "400", "4000"):
        assert v not in text
    assert all(k in chk.MAIN_FIELDS for k in f)
    assert not any(w in k for k in f for w in ("score", "skill", "rms", "n_", "count", "size"))


# ------------------------------------------------------------------------------------------------ the row order
def test_truth_rows_shared_follows_the_counts_only():
    assert tb.truth_rows_shared(5000, 5000, 600)
    assert tb.truth_rows_shared(601, 9000, 600)  # the same route of numpy's choice
    assert not tb.truth_rows_shared(600, 5000, 600)  # a permutation route against a subsample route
    assert not tb.truth_rows_shared(5000, 300, 600)  # the copy below the truth: another common count


def test_copy_parity_order_matches_the_copy_and_keeps_every_cell(chk):
    stage = _cloud(4000, 11)
    file = stage[np.random.default_rng(3).permutation(len(stage))[:1500]]
    copy = stage[np.sort(np.random.default_rng(4).choice(len(stage), 400, replace=False))]
    X = np.random.default_rng(8).gamma(2.0, 1.0, size=(len(file), 5)).astype(np.float32)
    n_true = 300
    assert tb.truth_rows_shared(len(file), len(copy), n_true)
    want = ev.prefix_parity(copy, n_true)
    for s in (1, 2, 3, 4):
        bad = file[_parity_flip_order(file, want * -1, n_true)] if s % 2 else file  # either parity as built
        order, rec = tb.copy_parity_order(bad, copy, n_true, seed=s)
        assert ev.prefix_parity(bad[order], n_true) == want
        assert chk.prefix_facts(bad[order], copy, n_true)["same_mode_whatever_truth"]
        assert sorted(order.tolist()) == list(range(len(bad)))  # a permutation: no cell added, lost or repeated
        assert rec == {"row_order": "copy_parity", "row_order_stream": tb.COPY_PARITY_STREAM,
                       "row_order_metric_seed": 0, "row_order_matched": True}
    # the cells travel whole: expression row and position together; distances and scale unchanged
    order, _ = tb.copy_parity_order(file, copy, n_true, seed=7)
    pairs = {(tuple(c), tuple(x)) for c, x in zip(file.tolist(), X.tolist(), strict=True)}
    assert {(tuple(c), tuple(x)) for c, x in zip(file[order].tolist(), X[order].tolist(), strict=True)} == pairs
    np.testing.assert_array_equal(np.sort(pdist(file[order].astype(float))), np.sort(pdist(file.astype(float))))
    assert tb.rms_radius(file[order]) == pytest.approx(tb.rms_radius(file), rel=1e-12)
    # under the main build the order changes nothing the check sees beyond the draw
    f = chk.frame_facts(file[order], copy, refs=[])
    assert f["axes_agree"] and f["handedness_agree"]
    # deterministic in the seed; the record says nothing about the file as built
    np.testing.assert_array_equal(tb.copy_parity_order(file, copy, n_true, seed=7)[0], order)


def test_copy_parity_order_refuses_without_shared_truth_rows():
    stage = _cloud(3000, 12)
    with pytest.raises(ValueError, match="different truth rows"):
        tb.copy_parity_order(stage[:2000], stage[:200], 600, seed=1)  # the copy below the truth count
    with pytest.raises(ValueError, match="values not shown") as e:
        tb.copy_parity_order(stage[:2000], stage[:200], 600, seed=1)
    assert "600" not in str(e.value) and "200" not in str(e.value)


# ------------------------------------------------------------------------------------------------ the builder
def _stage(n: int, seed: int, shift: float) -> ad.AnnData:
    rng = np.random.default_rng(seed)
    types = rng.integers(0, 3, n)
    lam = np.full((n, len(GENES)), 0.5)
    for t in range(3):
        lam[types == t, 4 * t:4 * t + 4] += 3.0 + shift
    C = rng.poisson(lam).astype(np.int32)
    C[C.sum(1) == 0, 0] = 1
    X = np.log1p(1e4 * C / C.sum(1, keepdims=True)).astype(np.float32)
    a = ad.AnnData(X=sparse.csr_matrix(X), obs=pd.DataFrame({"celltype": [f"t{t}" for t in types]},
                                                             index=[f"c{seed}_{i}" for i in range(n)]),
                   var=pd.DataFrame(index=GENES))
    a.layers["counts"] = sparse.csr_matrix(C)
    centres = np.array([[-4.0, 0.0, 0.0], [0.0, 3.0, 0.0], [4.0, 0.0, 1.0]]) * (1.0 + shift)
    a.obsm["spatial_3D"] = (centres[types] + rng.normal(size=(n, 3)) * [1.5, 1.0, 0.5]).astype(np.float32)
    return a


def embryo_test_board() -> Board:
    return Board(key="T2:embryo:test_interp", task="T2", genes=list(GENES), needs_coords=True, split="test_interp",
                 setting="embryo", min_cells=50, max_cells=SECRET_MAX, truth_cells=SECRET_TRUTH)


def install(bt, setattr_, setitem_, tmp: Path, reads: list) -> None:
    """Point build_test's embryo builder at synthetic E7.5 and E8.0 stages, a synthetic test board and size JSON. The
    interpolation runs with 150 anchors (the builder's own default is 2000; the code path is the same)."""
    stages = {"E7.5.h5ad": _stage(N_BEFORE, 1, 0.0), "E8.0.h5ad": _stage(N_AFTER, 2, 0.3)}
    board = embryo_test_board()

    def load(path, genes=None, drop_duplicates=False):
        name = Path(path).name
        assert Path(path).parent.name == "T2" and genes == list(GENES) and drop_duplicates, path
        reads.append(name)
        return stages[name].copy()

    size = tmp / "size_embryo.json"
    size.write_text(json.dumps({"valid": True, "preregistered": {"z_mode": "raw", "pred_rms_full": 6.0}}))
    setattr_(bt, "load_stage", load)
    setattr_(bt, "load_boards", lambda: {board.key: board})
    setattr_(bt, "_t2_count", lambda b, n: {"decision": "synthetic"})
    setattr_(bt.t2_build, "interpolate_clouds", partial(tb.interpolate_clouds, n_anchor=150))
    setitem_(bt.SIZES, "embryo", size)
    if hasattr(bt, "read_stage"):  # the released stage the option's floor copy is drawn from
        setattr_(bt, "read_stage", lambda path, genes=None: (reads.append("released:" + Path(path).name),
                                                             stages[Path(path).name].copy())[1])


def _shows(text: str, *values: int) -> bool:
    """Whether any of the integers appears in the text as a number of its own (not as digits inside a float)."""
    return any(re.search(rf"(?<![\d.]){v}(?![\d.])", text) for v in values)


def sha(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


@pytest.fixture(scope="module")
def bt():
    sys.path.insert(0, str(REPO / "experiments" / "builds"))
    import build_test

    return build_test


@pytest.fixture
def embryo_env(bt, monkeypatch, own_repo, tmp_path):
    reads: list[str] = []
    install(bt, monkeypatch.setattr, monkeypatch.setitem, tmp_path, reads)
    return reads


ARMS = {"moved": ["--method", "fgw_shape0"], "floor_exact": ["--method", "fgw_shape0", "--positions", "floor_exact"]}


def _build(bt, tmp_path, name: str, *argv: str) -> tuple[Path, dict, str]:
    out = tmp_path / name / "T2embryo.h5ad"
    bt.main(["embryo", *argv, "--out", str(out)])
    text = out.with_suffix(".build.json").read_text()
    return out, json.loads(text), text


@pytest.mark.parametrize("arm", list(ARMS))
def test_the_default_path_writes_the_old_bytes(bt, embryo_env, tmp_path, capsys, arm):
    out, rec, text = _build(bt, tmp_path, "plain", *ARMS[arm])
    assert sha(out) == PRE_OPTION_SHA256[arm]
    assert embryo_env == ["E7.5.h5ad", "E8.0.h5ad"]
    assert not any(k.startswith("row_order") for k in rec["record"])
    assert not _shows(text + capsys.readouterr().out, SECRET_MAX)


def test_the_option_is_off_by_default_and_as_built_is_the_default(bt, embryo_env, tmp_path):
    ap = bt.build_parser()
    assert ap.parse_args(["embryo", *ARMS["floor_exact"], "--out", "runs/x/e.h5ad"]).row_order == "as_built"
    a, rec_a, text_a = _build(bt, tmp_path, "a", *ARMS["floor_exact"])
    b, rec_b, text_b = _build(bt, tmp_path, "b", *ARMS["floor_exact"], "--row-order", "as_built")
    assert a.read_bytes() == b.read_bytes() and text_a == text_b
    assert "released:E7.5.h5ad" not in embryo_env  # the default reads no released copy


def _file(path: Path) -> tuple[np.ndarray, np.ndarray]:
    a = ad.read_h5ad(path)
    return np.asarray(a.X.toarray() if sparse.issparse(a.X) else a.X), np.asarray(a.obsm["spatial_3D"])


@pytest.mark.parametrize("arm", list(ARMS))
def test_copy_parity_reorders_the_file_into_the_copys_parity_and_changes_nothing_else(bt, embryo_env, tmp_path,
                                                                                    capsys, arm, chk):
    base, _, _ = _build(bt, tmp_path, "base", *ARMS[arm])
    out, rec, text = _build(bt, tmp_path, "cp", *ARMS[arm], "--row-order", "copy_parity")
    assert "released:E7.5.h5ad" in embryo_env
    X0, C0 = _file(base)
    X1, C1 = _file(out)
    # the same cells: every expression row with its own position, nothing added, lost or moved
    def cells(X, C):
        return sorted(map(tuple, np.hstack([C.astype(np.float64), X.astype(np.float64)]).tolist()))

    assert cells(X1, C1) == cells(X0, C0) and not np.array_equal(C1, C0)
    np.testing.assert_array_equal(np.sort(pdist(C1.astype(float))), np.sort(pdist(C0.astype(float))))
    assert tb.rms_radius(C1) == pytest.approx(tb.rms_radius(C0), rel=1e-9)
    # the copy's pre-fix parity, against the same truth rows, at the board's truth count (read by code)
    released = _stage(N_BEFORE, 1, 0.0)
    copy = ev.subsample(ev.Stage(X=np.zeros((N_BEFORE, 0)), coords=np.asarray(released.obsm["spatial_3D"])), 0.1,
                        0).coords
    np.testing.assert_array_equal(bt.floor_copy_coords("T2/E7.5.h5ad"), copy)
    assert chk.prefix_facts(C1, copy, SECRET_TRUTH)["same_mode_whatever_truth"]
    # the record says what was done and holds no count, parity or number of tries
    r = rec["record"]
    assert r["row_order"] == "copy_parity" and r["row_order_matched"] is True
    assert r["row_order_copy"] == "vec.evaluate.subsample(T2/E7.5.h5ad as released, 0.1, 0)"
    assert not any("parity_" in k or "tries" in k for k in r)
    assert not _shows(text + capsys.readouterr().out, SECRET_MAX, SECRET_TRUTH)


def test_copy_parity_refuses_when_the_truth_rows_differ_and_shows_no_value(bt, embryo_env, tmp_path, capsys,
                                                                         monkeypatch):
    big = Board(**{**embryo_test_board().__dict__, "truth_cells": 390})  # the copy's 150 rows are below it
    monkeypatch.setattr(bt, "load_boards", lambda: {big.key: big})
    out = tmp_path / "r" / "T2embryo.h5ad"
    with pytest.raises(SystemExit, match="different truth rows") as e:
        bt.main(["embryo", *ARMS["floor_exact"], "--row-order", "copy_parity", "--out", str(out)])
    assert not out.exists()
    assert not _shows(str(e.value) + capsys.readouterr().out, 390, 150, SECRET_MAX)
