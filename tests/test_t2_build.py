"""vec.t2_build on synthetic data only (no competition data, no board anchors).

Size JSONs are produced by experiments/coords/size_procedure.py itself on synthetic stage geometries with random
growth, so the reader is tested against the procedure's real output schema; none of the numbers means anything.
"""

from __future__ import annotations

import hashlib
import importlib.util
import io
import sys
import tokenize
from pathlib import Path

import anndata as ad
import numpy as np
import pytest
from conftest import GENES, simulate_stage

from vec import t2_build as tb
from vec.boards import Board
from vec.evaluate import SCORER_BUILDS, Stage, score_board
from vec.io import COORD_KEY, check_submission, read_stage

REPO = Path(__file__).resolve().parents[1]


def t2_board(min_cells=10, max_cells=3000, truth_cells=None, key="T2:synthetic", task="T2", needs_coords=True):
    return Board(key=key, task=task, genes=GENES, needs_coords=needs_coords, min_cells=min_cells,
                 max_cells=max_cells, truth_cells=truth_cells)


def stack(n, seed=0, extent=(120.0, 80.0, 60.0), step=15.0):
    """A uniformly filled box cut into z sections, rows stored section by section in descending z (like the released
    MERFISH files)."""
    rng = np.random.default_rng(seed)
    C = rng.uniform(-1.0, 1.0, size=(n, 3)) * np.asarray(extent)
    C[:, 2] = np.round(C[:, 2] / step) * step
    return C[np.argsort(-C[:, 2], kind="stable")]


def random_rotation(seed):
    Q, R = np.linalg.qr(np.random.default_rng(seed).normal(size=(3, 3)))
    Q = Q @ np.diag(np.sign(np.diag(R)))
    if np.linalg.det(Q) < 0:
        Q[:, 0] *= -1
    return Q


@pytest.fixture(scope="module")
def sp():
    """experiments/coords/size_procedure.py, loaded by path."""
    spec = importlib.util.spec_from_file_location("coords_size_procedure",
                                                  REPO / "experiments" / "coords" / "size_procedure.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


# ---------------------------------------------------------------------------------------------------- cell count
def test_count_is_max_cells_with_a_record():
    n, rec = tb.n_cells(t2_board(min_cells=100, max_cells=3000, truth_cells=800), pool_size=5000)
    assert n == 3000
    assert rec["rule"] == "max_cells" and rec["floor"] == 800 and rec["floor_source"] == "truth_cells"
    assert rec["exceeds_pool"] is False and rec["repeated_rows"] == 0 and rec["fraction_of_max_cells"] == 1.0


def test_count_never_below_truth_cells_and_overrides_are_refused_not_clipped():
    board = t2_board(min_cells=100, max_cells=3000, truth_cells=800)
    assert tb.n_cells(board, 5000, override="0.5max")[0] == 1500
    assert tb.n_cells(board, 5000, override="1max")[0] == 3000
    assert tb.n_cells(board, 5000, override=800)[0] == 800  # exactly truth_cells is allowed
    for bad in (799, "0.25max", 3001):  # below truth_cells (750 cells for 0.25max) or above max_cells
        with pytest.raises(ValueError, match="refused"):
            tb.n_cells(board, 5000, override=bad)
    for bad in ("1.5max", "half", "0max", "2x"):
        with pytest.raises(ValueError):
            tb.n_cells(board, 5000, override=bad)
    for bad in (0.5, True):
        with pytest.raises(TypeError):
            tb.n_cells(board, 5000, override=bad)


def test_count_floor_is_min_cells_without_truth_cells_and_bad_boards_are_refused():
    n, rec = tb.n_cells(t2_board(min_cells=100, max_cells=3000), pool_size=5000, override="0.25max")
    assert n == 750 and rec["floor"] == 100 and rec["floor_source"] == "min_cells" and rec["truth_cells"] is None
    with pytest.raises(ValueError, match="no max_cells"):
        tb.n_cells(t2_board(max_cells=None), 5000)
    with pytest.raises(ValueError, match="inconsistent"):
        tb.n_cells(t2_board(min_cells=10, max_cells=500, truth_cells=800), 5000)
    with pytest.raises(ValueError):
        tb.n_cells(t2_board(), 0)


def test_count_above_the_pool_is_annotated():
    n, rec = tb.n_cells(t2_board(min_cells=100, max_cells=3000), pool_size=1200)
    assert n == 3000 and rec["exceeds_pool"] and rec["repeated_rows"] == 1800
    assert rec["x_pool"] == pytest.approx(2.5) and "2 or 3 times" in rec["note"]


# ---------------------------------------------------------------------------------------------------- sampling
def test_sampling_without_replacement_is_random_rows_in_random_order():
    pool = stack(2000)
    s = tb.sample_cells(pool, 500, seed=3)
    assert s.rows.shape == (500,) and np.unique(s.rows).size == 500 and s.rows.min() >= 0 and s.rows.max() < 2000
    assert not (np.all(np.diff(s.rows) > 0) or np.all(np.diff(s.rows) < 0))  # neither a block nor file order
    np.testing.assert_array_equal(s.coords, pool[s.rows])  # nothing to jitter
    assert s.record["with_replacement"] is False and s.record["n_jittered"] == 0 and s.record["max_shift"] == 0.0
    np.testing.assert_array_equal(tb.sample_cells(pool, 500, seed=3).rows, s.rows)  # deterministic per seed
    assert not np.array_equal(tb.sample_cells(pool, 500, seed=4).rows, s.rows)


def test_sampling_above_the_pool_repeats_rows_evenly_and_jitter_breaks_every_duplicate_within_the_bound():
    pool = stack(2000, seed=1)
    n = 5300
    s = tb.sample_cells(pool, n, seed=0)
    counts = np.bincount(s.rows, minlength=2000)
    assert set(np.unique(counts)) == {2, 3} and counts.sum() == n
    assert s.record["copies_per_row"] == [2, 3] and s.record["n_repeated_rows"] == n - 2000
    assert s.record["n_jittered"] == n - 2000 == s.record["dup_from_repeated_rows"]
    assert np.unique(s.coords, axis=0).shape[0] == n
    assert np.unique(s.coords.astype(np.float32), axis=0).shape[0] == n  # as the scorer reads them
    nn = tb.median_nn_distance(pool)
    assert s.record["median_nn"] == pytest.approx(nn)
    shift = np.linalg.norm(s.coords - pool[s.rows], axis=1)
    moved = shift > 0
    assert moved.sum() == n - 2000
    bound = tb.JITTER_NN_FRAC * nn
    assert shift.max() <= bound * (1 + 1e-12) and shift[moved].min() >= tb.JITTER_MIN_FRAC * bound * (1 - 1e-12)
    # the first occurrence of each row keeps its exact position
    _, first = np.unique(s.rows, return_index=True)
    np.testing.assert_array_equal(s.coords[first], pool[s.rows[first]])
    # isotropic: no axis is preferred
    d = (s.coords - pool[s.rows])[moved]
    assert np.all(np.abs(d.mean(0)) < 0.1 * bound)
    assert np.all(np.abs(np.abs(d).mean(0) / np.abs(d).mean() - 1) < 0.15)


def test_sampling_breaks_duplicate_positions_already_in_the_pool():
    pool = stack(1000, seed=2)
    pool[11] = pool[10]
    pool[500] = pool[10]
    pool[700] = pool[20] + np.array([0.0, 0.0, 1e-12])  # distinct in float64, one position as the scorer reads it
    assert np.unique(pool, axis=0).shape[0] == 998
    s = tb.sample_cells(pool, 1000, seed=0)
    assert s.record["dup_from_pool_positions"] == 3 and s.record["dup_from_repeated_rows"] == 0
    assert np.unique(s.coords.astype(np.float32), axis=0).shape[0] == 1000


def test_sampling_inputs_and_refusals():
    pool = stack(800, seed=3)
    a = ad.AnnData(X=np.zeros((800, len(GENES)), np.float32))
    a.obsm[COORD_KEY] = pool.astype(np.float32)
    st = Stage(X=np.zeros((800, len(GENES)), np.float32), coords=pool.astype(np.float32))
    rows = [tb.sample_cells(x, 300, seed=5).rows for x in (pool, a, st)]
    np.testing.assert_array_equal(rows[0], rows[1])
    np.testing.assert_array_equal(rows[0], rows[2])
    s = tb.sample_cells(pool, 2000, seed=1, jitter_frac=0.1)
    assert s.record["max_shift"] <= 0.1 * tb.median_nn_distance(pool) * (1 + 1e-12)
    with pytest.raises(ValueError, match="free budget"):
        tb.sample_cells(pool, 900, seed=0, jitter_frac=0.3)
    with pytest.raises(ValueError):
        tb.sample_cells(pool, 0, seed=0)
    with pytest.raises(ValueError):
        tb.sample_cells(pool, 10, seed=-1)
    with pytest.raises(ValueError, match="fewer than two distinct"):
        tb.sample_cells(np.zeros((1, 3)), 2, seed=0)
    with pytest.raises(KeyError):
        tb.sample_cells(ad.AnnData(X=np.zeros((5, 2), np.float32)), 2, seed=0)


# ---------------------------------------------------------------------------------------------------- coordinates
def test_centre_and_rms_radius():
    C = stack(500, seed=4) + np.array([1000.0, -50.0, 7.0])
    np.testing.assert_allclose(tb.centre(C).mean(0), 0.0, atol=1e-9)
    X = C - C.mean(0)
    assert tb.rms_radius(C) == pytest.approx(np.sqrt((X ** 2).sum(1).mean()))


def test_proper_rotation_rotates_about_the_centroid_and_refuses_improper_matrices():
    C = stack(400, seed=5) + 300.0
    R = random_rotation(0)
    D = tb.proper_rotation(C, R)
    np.testing.assert_allclose(D.mean(0), C.mean(0), atol=1e-9)
    np.testing.assert_allclose(np.linalg.norm(D[1:] - D[:-1], axis=1), np.linalg.norm(C[1:] - C[:-1], axis=1),
                               rtol=1e-12)
    assert tb.rms_radius(D) == pytest.approx(tb.rms_radius(C))
    np.testing.assert_allclose(tb.proper_rotation(C, R, about=np.zeros(3)), C @ R.T, atol=1e-9)
    mirror = np.diag([1.0, 1.0, -1.0])
    for improper in (mirror, R @ mirror, -np.eye(3)):
        with pytest.raises(ValueError, match="improper"):
            tb.proper_rotation(C, improper)
    with pytest.raises(ValueError, match="not orthonormal"):
        tb.proper_rotation(C, 2.0 * np.eye(3))
    with pytest.raises(ValueError):
        tb.check_rotation(np.eye(2))


def test_rescale_to_rms_is_proper_isotropic_and_matches_the_size_procedure(sp):
    C = stack(600, seed=6) + 40.0
    D = tb.rescale_to_rms(C, 321.0)
    assert tb.rms_radius(D) == pytest.approx(321.0)
    np.testing.assert_allclose(D.mean(0), C.mean(0), atol=1e-9)
    f = 321.0 / tb.rms_radius(C)
    np.testing.assert_allclose(D - D.mean(0), f * (C - C.mean(0)), atol=1e-9)
    np.testing.assert_allclose(D, sp.rescale_to_rms(C, 321.0), rtol=0, atol=1e-9)
    for bad in (0.0, -1.0, float("nan"), float("inf")):
        with pytest.raises(ValueError):
            tb.rescale_to_rms(C, bad)
    with pytest.raises(ValueError, match="no spread"):
        tb.rescale_to_rms(np.ones((5, 3)), 1.0)


# ---------------------------------------------------------------------------------------------------- size JSON
def _geom(sp, name, t, scale, sd, seed):
    """A synthetic stage: a core blob inside a wider shell, scaled by `scale`, z on 15-unit sections."""
    rng = np.random.default_rng(seed)
    core = rng.normal(size=(1500, 3)) * np.asarray(sd) * 0.6
    shell = rng.normal(size=(1500, 3)) * np.asarray(sd)
    C = np.vstack([core, shell]) * scale
    C[:, 2] = np.round(C[:, 2] / 15.0) * 15.0
    return sp.geometry_from_arrays(name, t, C, np.r_[np.ones(1500, bool), np.zeros(1500, bool)])


def _size_json(sp, board, geoms, target_t, path):
    """What size_procedure.run writes, built from synthetic geometries (no files read)."""
    rep = sp.heart_report(geoms, target_t) if board == "heart" else sp.embryo_report(geoms, target_t)
    res = {"board": board, "target_t": float(target_t), "timestamp": "2026-10-19T00:00:00Z", "valid": rep["valid"],
           "preregistered_settings": dict(sp.PREREGISTERED),
           "provenance": {"git_commit": "synthetic", "git_uncommitted": [], "size_procedure_sha256": "synthetic"},
           "inputs": [], **rep}
    path.write_text(sp.dumps(res))
    return res


def _heart(sp, tmp_path, sd, name="size_heart.json", seed=0):
    rng = np.random.default_rng(seed)
    scales = np.cumprod(1.0 + rng.uniform(0.05, 0.6, size=3))  # random growth
    geoms = [_geom(sp, f"s{i}", t, s, sd, seed + i) for i, (t, s) in enumerate(zip((8.75, 9.5, 10.5), scales,
                                                                                    strict=True))]
    path = tmp_path / name
    return _size_json(sp, "heart", geoms, 12.5, path), path


@pytest.mark.parametrize("hypothesis", ["H1", "H0", "H2"])
def test_size_raw_mode_is_an_isotropic_rescale_to_the_hypothesis_rms(sp, tmp_path, hypothesis):
    res, path = _heart(sp, tmp_path, sd=(30.0, 80.0, 100.0))  # core thin axis along x: auto resolves to raw
    assert res["preregistered"]["z_mode"] == "raw" and res["valid"] is True
    C = stack(900, seed=7)
    out, rec = tb.apply_size_procedure(C, path, hypothesis, expect_board="heart")
    target = res["preregistered"]["scope"][hypothesis]["pred_rms_full"]
    assert tb.rms_radius(out) == pytest.approx(target) and rec["target_rms_full"] == pytest.approx(target)
    np.testing.assert_allclose(out.mean(0), C.mean(0), atol=1e-9)
    np.testing.assert_allclose(out - out.mean(0), rec["scale"] * (C - C.mean(0)), atol=1e-9)
    assert rec["z_mode_recorded"] == rec["z_mode_used"] == "raw" and rec["z_held"] is False
    assert rec["size_json_sha256"] == hashlib.sha256(path.read_bytes()).hexdigest()
    assert rec["git_commit"] == "synthetic" and rec["warnings"] == []


def test_size_z_hold_mode_leaves_z_native(sp, tmp_path):
    res, path = _heart(sp, tmp_path, sd=(100.0, 80.0, 30.0))  # core thin axis along z: auto resolves to z hold
    pre = res["preregistered"]
    assert pre["z_mode"] == "xy_grow_z_hold"
    C = stack(900, seed=8)
    out, rec = tb.apply_size_procedure(C, path, "H2")
    np.testing.assert_array_equal(out[:, 2], C[:, 2])  # z untouched, bit for bit
    assert tb.rms_radius(out) == pytest.approx(pre["scope"]["H2"]["pred_rms_full"])
    assert rec["z_held"] is True and rec["target_rms_xy_full"] == pytest.approx(pre["scope"]["H2"]["pred_rms_xy_full"])
    out_xy, rec_xy = tb.apply_size_procedure(C, path, "H2", match="xy")
    X = out_xy - out_xy.mean(0)
    assert np.sqrt((X[:, :2] ** 2).sum(1).mean()) == pytest.approx(pre["scope"]["H2"]["pred_rms_xy_full"])
    np.testing.assert_array_equal(out_xy[:, 2], C[:, 2])
    # a logged alternative can be named for diagnosis; the record says which mode was used
    out_raw, rec_raw = tb.apply_size_procedure(C, path, "H1", z_mode="raw")
    assert rec_raw["z_mode_recorded"] == "xy_grow_z_hold" and rec_raw["z_mode_used"] == "raw"
    assert tb.rms_radius(out_raw) == pytest.approx(res["z_modes"]["raw"]["scope"]["H1"]["pred_rms_full"])


def test_size_embryo_has_no_hypothesis(sp, tmp_path):
    rng = np.random.default_rng(3)
    scales = np.cumprod(1.0 + rng.uniform(0.05, 0.6, size=4))
    geoms = [_geom(sp, f"e{i}", t, s, (60.0, 50.0, 90.0), 10 + i)
             for i, (t, s) in enumerate(zip((6.75, 7.25, 7.5, 8.0), scales, strict=True))]
    res = _size_json(sp, "embryo", geoms, 7.75, tmp_path / "size_embryo.json")
    C = stack(700, seed=9)
    out, rec = tb.apply_size_procedure(C, tmp_path / "size_embryo.json", expect_board="embryo")
    assert tb.rms_radius(out) == pytest.approx(res["preregistered"]["pred_rms_full"])
    assert rec["z_mode_used"] == "raw" and rec["reference"] == res["preregistered"]["bracket"]
    with pytest.raises(ValueError, match="no scope hypotheses"):
        tb.apply_size_procedure(C, tmp_path / "size_embryo.json", "H1")
    out_hold, _ = tb.apply_size_procedure(C, res, z_mode="xy_grow_z_hold")  # a dict works too
    np.testing.assert_array_equal(out_hold[:, 2], C[:, 2])
    assert tb.rms_radius(out_hold) == pytest.approx(res["z_modes"]["xy_grow_z_hold"]["pred_rms_full"])


def test_size_refusals(sp, tmp_path):
    res, path = _heart(sp, tmp_path, sd=(30.0, 80.0, 100.0))
    C = stack(300, seed=10)
    with pytest.raises(ValueError, match="scope hypothesis"):
        tb.apply_size_procedure(C, path)
    with pytest.raises(ValueError, match="scope hypothesis"):
        tb.apply_size_procedure(C, path, "H3")
    with pytest.raises(ValueError, match="diagnosis only"):
        tb.apply_size_procedure(C, path, "H1", z_mode="split")
    with pytest.raises(ValueError, match="cannot be applied"):
        tb.apply_size_procedure(C, path, "H1", z_mode="auto")
    with pytest.raises(ValueError, match="expected embryo"):
        tb.apply_size_procedure(C, path, "H1", expect_board="embryo")
    with pytest.raises(ValueError, match="only to the xy_grow_z_hold"):
        tb.apply_size_procedure(C, path, "H1", match="xy")
    with pytest.raises(ValueError, match="not a size procedure output"):
        tb.apply_size_procedure(C, res["preregistered"], "H1")
    # dry runs carry redacted forecasts: refused
    with pytest.raises(ValueError, match="redacted"):
        tb.apply_size_procedure(C, sp.redact_forecasts(res), "H1")
    # no growth rule applicable (no rate interval starts at or after E9.0): valid=false, refused
    early = [_geom(sp, "a", 8.25, 1.0, (30.0, 80.0, 100.0), 1), _geom(sp, "b", 8.75, 1.2, (30.0, 80.0, 100.0), 2)]
    invalid = _size_json(sp, "heart", early, 12.5, tmp_path / "invalid.json")
    assert invalid["valid"] is False
    with pytest.raises(ValueError, match="valid=false"):
        tb.apply_size_procedure(C, tmp_path / "invalid.json", "H1")
    # uncommitted procedure files are reported
    res["provenance"]["git_uncommitted"] = [" M experiments/coords/size_procedure.py"]
    assert "uncommitted" in tb.apply_size_procedure(C, res, "H1")[1]["warnings"][0]


# ---------------------------------------------------------------------------------------------------- guard
def test_guard_passes_random_subsets_and_flags_contiguous_blocks():
    pool = stack(20000, seed=11)
    for m in (1000, 4000):
        g = tb.cloud_guard(tb.sample_cells(pool, m, seed=1).coords, pool, seed=0)
        assert g.ok, g.flags
        assert 0.9 < g.cover_ratio < 1.1 and abs(np.log(g.rms_ratio)) < 0.03
        mid = (len(pool) - m) // 2
        for block in (pool[mid:mid + m], pool[:m]):  # mid-file slab and the first rows (what a[:n] submits)
            gb = tb.cloud_guard(block, pool, seed=0)
            assert not gb.ok and any(f.startswith("block") for f in gb.flags)
            assert gb.cover_ratio < 0.6
    assert tb.cloud_guard(pool[:1000], pool, seed=0).as_dict()["ok"] is False


def test_guard_handles_translation_rescaling_rotation_and_oversampling():
    pool = stack(20000, seed=12)
    sub = tb.sample_cells(pool, 3000, seed=2).coords
    assert tb.cloud_guard(tb.centre(sub) + 500.0, pool).ok  # translation drops out
    big = tb.rescale_to_rms(sub, 1.8 * tb.rms_radius(sub))
    assert any(f.startswith("rms") for f in tb.cloud_guard(big, pool).flags)
    g = tb.cloud_guard(big, pool, expected_rms_ratio=tb.rms_radius(big) / tb.rms_radius(pool))
    assert g.ok and g.expected_rms_source == "argument"
    rot = tb.proper_rotation(sub, random_rotation(1))
    assert not tb.cloud_guard(rot, pool).ok  # frame="source" assumes no rotation
    assert tb.cloud_guard(rot, pool, frame="own").ok
    assert not tb.cloud_guard(pool[:3000], pool, frame="own").ok
    over = tb.sample_cells(pool, 26000, seed=3).coords  # copy-like above the pool size
    g = tb.cloud_guard(over, pool)
    assert g.ok and g.n_pred == 26000
    with pytest.raises(ValueError, match="frame"):
        tb.cloud_guard(sub, pool, frame="pca")
    with pytest.raises(ValueError, match="at least"):
        tb.cloud_guard(sub[:3], pool)
    with pytest.raises(ValueError):
        tb.cloud_guard(sub, pool, expected_rms_ratio=0.0)


def test_guard_weighted_reference_for_composition_weighted_draws():
    left = stack(10000, seed=15, extent=(60.0, 80.0, 60.0)) + np.array([-70.0, 0.0, 0.0])
    right = stack(10000, seed=16, extent=(60.0, 80.0, 60.0)) + np.array([70.0, 0.0, 0.0])
    pool = np.vstack([left, right])  # two cell types in two regions
    w = np.r_[np.ones(10000), np.full(10000, 0.03)]  # the right-hand type strongly down-weighted
    pred = pool[np.random.default_rng(17).choice(20000, 2000, replace=False, p=w / w.sum())]
    assert not tb.cloud_guard(pred, pool).ok  # uniform draws are the wrong reference for a weighted draw
    g = tb.cloud_guard(pred, pool, source_weights=w)
    assert g.ok and g.weighted_reference, g.flags
    mid = 10000 - 1000
    assert not tb.cloud_guard(pool[mid:mid + 2000], pool, source_weights=np.ones(20000)).ok
    with pytest.raises(ValueError, match="source_weights"):
        tb.cloud_guard(pred, pool, source_weights=w[:-1])
    with pytest.raises(ValueError, match="source_weights"):
        tb.cloud_guard(pred, pool, source_weights=-w)


# ---------------------------------------------------------------------------------------------------- builder
def _cells(n, seed=0):
    rng = np.random.default_rng(seed)
    counts = rng.poisson(rng.gamma(2.0, 1.0, size=(1, len(GENES))), size=(n, len(GENES)))
    lib = counts.sum(1, keepdims=True).clip(min=1)
    return np.log1p(1e4 * counts / lib).astype(np.float32)


def _counts(n, n_genes, seed=0):
    """MERFISH-like raw counts: a few per gene, a median library of about 150, no empty cell."""
    rng = np.random.default_rng(seed)
    mu = rng.gamma(0.5, 300.0 / n_genes, size=(1, n_genes))
    counts = rng.poisson(mu * rng.lognormal(0.0, 0.4, size=(n, 1)), size=(n, n_genes)).astype(np.float64)
    counts[counts.sum(1) == 0, 0] = 1.0
    return counts


def _cp10k(counts):
    return 1e4 * counts / counts.sum(1, keepdims=True)


def test_builder_orders_genes_and_passes_the_format_check(tmp_path):
    board = t2_board(min_cells=10, max_cells=1000)
    X = _cells(200)
    order = np.random.default_rng(1).permutation(len(GENES))
    extra = np.ones((200, 2), np.float32)
    a = ad.AnnData(X=np.hstack([X[:, order], extra]))
    a.var_names = [GENES[i] for i in order] + ["Extra1", "Extra2"]
    coords = stack(200, seed=13)
    sub, rep = tb.build_submission(board, a, coords, out=tmp_path / "pred.h5ad")
    assert rep.ok and list(sub.var_names) == GENES
    np.testing.assert_array_equal(sub.X, X)
    assert sub.obsm[COORD_KEY].dtype == np.float32 and sub.obsm[COORD_KEY].shape == (200, 3)
    assert any("outside the board panel dropped" in w for w in rep.warnings)
    assert check_submission(tmp_path / "pred.h5ad", board).ok
    back = read_stage(tmp_path / "pred.h5ad", GENES)
    np.testing.assert_array_equal(np.asarray(back.X), X)
    # matrix with gene names, Stage with gene names, matrix already in panel order
    assert tb.build_submission(board, X[:, order], coords, genes=[GENES[i] for i in order])[1].ok
    assert tb.build_submission(board, Stage(X=X[:, order]), coords, genes=[GENES[i] for i in order])[1].ok
    _, rep = tb.build_submission(board, X, coords)
    assert rep.ok and any("panel order" in w for w in rep.warnings)


def test_builder_refusals(tmp_path):
    board = t2_board(min_cells=50, max_cells=300)
    X, coords = _cells(100), stack(100, seed=14)
    with pytest.raises(ValueError, match="no gene names"):
        tb.build_submission(board, X[:, :-1], coords)
    with pytest.raises(KeyError):
        tb.build_submission(board, X[:, :-1], coords, genes=GENES[:-1])
    with pytest.raises(ValueError, match="do not pass genes"):
        tb.build_submission(board, ad.AnnData(X=X), coords, genes=GENES)
    with pytest.raises(ValueError, match="negative"):
        tb.build_submission(board, X - 1.0, coords)
    bad = X.copy()
    bad[0, 0] = np.nan
    with pytest.raises(ValueError, match="NaN"):
        tb.build_submission(board, bad, coords)
    with pytest.raises(ValueError, match="coordinate rows"):
        tb.build_submission(board, X, coords[:99])
    with pytest.raises(ValueError, match="outside"):
        tb.build_submission(board, X[:40], coords[:40])
    with pytest.raises(ValueError, match="outside"):
        tb.build_submission(t2_board(min_cells=10, max_cells=99), X, coords)
    dup = coords.copy()
    dup[1] = dup[0]
    with pytest.raises(ValueError, match="repeat the exact"):
        tb.build_submission(board, X, dup)
    _, rep = tb.build_submission(board, X, dup, allow_duplicate_coords=True)
    assert rep.ok and any("exact coordinates" in w for w in rep.warnings)
    with pytest.raises(ValueError, match="takes no coordinates"):
        tb.build_submission(t2_board(task="T1", needs_coords=False), X, coords)
    counts = np.random.default_rng(0).poisson(3.0, size=(100, len(GENES))).astype(np.float32)
    _, rep = tb.build_submission(board, counts, coords)
    assert any("log1p CP10k" in w for w in rep.warnings)
    assert not list(tmp_path.iterdir())  # nothing written without `out`


# ---------------------------------------------------------------------------------------------------- end to end
def test_end_to_end_copy_like_submission_scores_under_both_builds(tmp_path, profiles):
    _, later = profiles
    centres = np.array([[0, 0, 0], [6, 0, 0], [0, 6, 0]], float)

    def draw(seed, n=4000):
        return simulate_stage(n, later, [0.3, 0.3, 0.4], np.random.default_rng(seed), coords_scale=3.0,
                              coords_shift=1.5 * centres)

    target, pool = draw(21), draw(22)
    base, _ = profiles
    ref = simulate_stage(4000, base, [0.5, 0.3, 0.2], np.random.default_rng(23), coords_scale=2.0,
                         coords_shift=centres)
    board = t2_board(min_cells=100, max_cells=1500, truth_cells=200)
    n, _ = tb.n_cells(board, pool.n)
    s = tb.sample_cells(pool, n, seed=0)
    assert tb.cloud_guard(s.coords, pool.coords, seed=0).ok
    tb.build_submission(board, pool.X[s.rows], s.coords, genes=GENES, out=tmp_path / "copy.h5ad")
    pred = Stage.from_anndata(read_stage(tmp_path / "copy.h5ad", GENES), GENES, need_coords=True)
    assert pred.n == 1500
    for build in SCORER_BUILDS:
        r = score_board(board, pred, target, ref, seed=0, build=build)
        assert r.extras["build"] == build and np.isfinite(r.score)
        assert r.score > 75.0, (build, r.group_skill)  # cells from the target's own distribution


# ---------------------------------------------------------------------------------------------------- no constants
# Every numeric literal allowed in vec/t2_build.py: indices and generic factors, the scorer's grid (16, 3.0), the
# jitter budget (0.25) and its lower end (0.5), guard thresholds (0.9, 0.03, 0.1), reference draws (5), the warning
# factor (2.0), the implied library band (0.1 and 2.0 x 1e4), the smallest cloud (4) and tolerances. A size, a count or
# a composition share of any stage fails.
ALLOWED_LITERALS = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 16.0, 0.25, 0.5, 0.9, 0.03, 0.1, 1e-6, 1e-9, 1e-12}


def test_module_holds_no_stage_constants():
    src = (REPO / "src" / "vec" / "t2_build.py").read_text()
    lits = {float(t.string) for t in tokenize.generate_tokens(io.StringIO(src).readline) if t.type == tokenize.NUMBER}
    assert lits <= ALLOWED_LITERALS, sorted(lits - ALLOWED_LITERALS)


def test_builder_enforces_the_truth_floor_unless_a_control_asks():
    board = t2_board(min_cells=10, max_cells=1000, truth_cells=300)
    X, coords = _cells(200), stack(200, seed=21)
    with pytest.raises(ValueError, match="below the board's truth_cells"):
        tb.build_submission(board, X, coords)
    assert tb.build_submission(board, X, coords, allow_below_truth=True)[1].ok


def test_no_raw_counts_warning_on_a_float32_single_gene_cell():
    board = t2_board(min_cells=10, max_cells=1000)
    X = _cells(50)
    X[0] = 0.0
    X[0, 0] = np.float32(np.log1p(1e4))  # a valid cell with all its counts in one gene, stored as float32
    _, rep = tb.build_submission(board, X, stack(50, seed=22))
    assert not any("raw counts" in w for w in rep.warnings)


def test_builder_warns_on_log1p_of_counts_never_scaled_to_1e4():
    board = t2_board(min_cells=10, max_cells=1000)
    counts, coords = _counts(300, len(GENES)), stack(300, seed=23)
    _, rep = tb.build_submission(board, np.log1p(counts), coords)  # fractional and below the ceiling: only size shows
    assert any("implied library" in w for w in rep.warnings), rep.warnings
    _, rep = tb.build_submission(board, np.log1p(_cp10k(counts)), coords)  # the same cells, scaled to 1e4 first
    assert not any("implied library" in w for w in rep.warnings), rep.warnings


def test_no_library_warning_on_released_like_cells_over_part_of_the_measured_panel():
    measured = [f"m{i}" for i in range(500)]
    board = Board(key="T2:synthetic", task="T2", genes=measured[:249] + measured[251:], needs_coords=True,
                  min_cells=10, max_cells=1000)  # 498 of the 500 measured genes
    X = np.log1p(_cp10k(_counts(300, 500, seed=1)))  # scaled over all 500 measured genes, as released
    _, rep = tb.build_submission(board, X, stack(300, seed=24), genes=measured)
    assert not any("implied library" in w or "raw counts" in w for w in rep.warnings), rep.warnings
    assert any("2 genes outside the board panel dropped" in w for w in rep.warnings)


def test_builder_warns_on_non_integer_values_above_the_log1p_ceiling():
    board = t2_board(min_cells=10, max_cells=1000)
    cp10k = _cp10k(_counts(300, len(GENES), seed=2))
    for X in (cp10k, np.log1p(100.0 * cp10k)):  # the log forgotten; log1p CPM
        nz = X[X > 0]
        assert not np.all(nz == np.rint(nz)) and X.max() > np.log1p(1e4)  # the max branch, not the whole-number one
        _, rep = tb.build_submission(board, X, stack(300, seed=25))
        raw = [w for w in rep.warnings if "raw counts" in w]
        assert raw and "whole numbers only: False" in raw[0], rep.warnings
