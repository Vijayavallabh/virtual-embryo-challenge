"""Row E57's step, RS (one run per side): vec.run_sort, its two harness arms (vec.harness.RS_CANDIDATES) and the
builder's off-by-default `--neighbourhood run_sorted` (experiments/builds/build_test.py embryo), on synthetic clouds and
synthetic stages only. Nothing is read from data/ or boards/, a synthetic TRIES.md stands in for the ledger, and nothing
is scored.

Pinned: the reassignment moves expression rows only, so every coordinate stays byte-identical in its order, every
row's expression is unchanged (the rows as a multiset) and every position keeps its pool; the midline is the plane that
mirrors the earlier bracket's rows onto their own labels; the swaps are linear_sum_assignment's on mirror-homologous
positions and move each bracket to its own side; the placebo matches RS's number of swaps per pool and side and leaves
the brackets' arrangement as it was; on O's command the builder writes O's file but for the expression rows' places;
every refusal comes before anything is read and shows no count; and every other embryo command gives what main's
builder before the option (BASE) gives, run in this same process from git's copy of build_test.py."""

from __future__ import annotations

import hashlib
import itertools
import json
import re
import sys
import types
from functools import partial
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import yaml
from conftest import git
from scipy import sparse
from scipy.optimize import linear_sum_assignment

from vec import harness
from vec import run_sort as rs
from vec import t2_build as tb
from vec.boards import Board
from vec.prep import LINEAGE_MAP

REPO = Path(__file__).resolve().parents[1]
BASE = "397b0e12996261446830d0507555e63e188302cf"  # main before --neighbourhood existed (the tip merged in, row R57)
PANEL = [f"g{i}" for i in range(12)]
EXTRA = ["x0", "x1", "x2"]  # genes outside the panel: part of each cell's library, never of a file
SECRET_MAX, SECRET_TRUTH = 400, 120  # stand for the test board's max_cells and truth_cells: never shown
# label: (lineage group, cells at E7.5, cells at E8.0, the blob's centre on the side x > 0)
WORLD = {"Epi": ("EpiblastStreak", 160, 140, (2.0, 0.5, 0.0)), "NeuralA": ("Neural", 140, 160, (1.5, 3.0, 0.5)),
         "Heart": ("CM", 120, 140, (2.5, -2.5, 0.0)), "ExM": ("ExtraEmbryonicMeso", 140, 120, (1.0, 0.0, 3.0)),
         "Gut": ("Endoderm", 120, 120, (0.5, -1.0, -2.0))}
DEPTH = {"E7.5": 60.0, "E8.0": 200.0}
FLOOR = ["--method", "fgw_shape0", "--positions", "floor_exact"]
OV = ["--shape-interp", "overlay", "--shape-interp-test-ok"]
O_ARGS = [*FLOOR, *OV]
RULING = "R60"
RS_ARGS = [*O_ARGS, "--neighbourhood", "run_sorted", "--captain-ruling", RULING]
W_ARGS = [*FLOOR, "--expression", "series_back", *OV, "--captain-ruling", "R48"]
RULING_ROW = (f"| {RULING} (the captain's ruling of 2026-11-25, ruling, no run: build T2embryo_s2x_ovrs.h5ad; the "
              f"embryo builder is to take `--neighbourhood run_sorted` on O's command only with `--captain-ruling "
              f"{RULING}`) | T2 embryo | synthetic |")
R48_ROW = ("| R48 (the captain's order of 2026-11-24, ruling, no run: \"Build T2embryo_s2x_ovsb.h5ad as named\"; the "
           "embryo builder is to take `--expression series_back` on O's command only with `--captain-ruling R48`; not "
           "a try) | T2 embryo | synthetic |")


# ------------------------------------------------------------------------------------------------ synthetic clouds
def rotation(seed: int) -> np.ndarray:
    Q, R = np.linalg.qr(np.random.default_rng(seed).normal(size=(3, 3)))
    Q = Q * np.sign(np.diag(R))
    return Q if np.linalg.det(Q) > 0 else -Q


def bilateral_cloud(seed: int, n: int = 1200, labels=("A", "B", "C", "D"), groups=("Neural", "CM", "Neural", "CM"),
                    turn: int | None = None, mix: float = 0.5):
    """Rows of two brackets on a cloud symmetric about the plane x = 0 (each label two mirror-image blobs), the brackets
    interleaved at random (as O's overlay interleaves two embryos), optionally turned by a random rotation. Returns
    coordinates, stage, labels, groups, expression rows (each row distinct) and the true normal."""
    rng = np.random.default_rng(seed)
    lab = rng.choice(len(labels), n)
    centre = np.array([[2.0, 0.5 * i, (-1) ** i * 1.0] for i in range(len(labels))])
    side = np.where(rng.random(n) < 0.5, -1.0, 1.0)
    C = centre[lab] * np.c_[side, np.ones(n), np.ones(n)] + rng.normal(scale=0.35, size=(n, 3))
    stage = (rng.random(n) < mix).astype(np.int64)
    X = rng.gamma(2.0, 1.0, size=(n, 6)).astype(np.float32)
    u = np.array([1.0, 0.0, 0.0])
    if turn is not None:
        Q = rotation(turn)
        C, u = C @ Q.T, Q @ u
    return C, stage, np.asarray(labels)[lab], np.asarray(groups)[lab], X, u


# ------------------------------------------------------------------------------------------------ the module
def test_the_pools_partition_the_lineage_maps_groups():
    groups = set(yaml.safe_load(Path(LINEAGE_MAP).read_text())["groups"])
    flat = [g for pool in rs.POOLS for g in pool]
    assert len(flat) == len(set(flat)) and set(flat) == groups
    assert rs.POOLS[0] == ("EpiblastStreak", "Neural", "Ectoderm", "NCC") and rs.POOLS[-1] == ("Unknown",)
    assert (rs.N_NORMALS, rs.MIN_LABEL_ROWS, rs.K) == (2000, 20, 15) == (2000, 20, harness.NBHD_K)
    with pytest.raises(rs.RunSortRefused, match="no pool"):
        rs.pool_index(["Neural", "Epiblast"])


def test_the_normals_are_a_fixed_fibonacci_hemisphere():
    U = rs.fibonacci_hemisphere()
    assert U.shape == (2000, 3) and np.allclose(np.linalg.norm(U, axis=1), 1.0) and (U[:, 2] > 0).all()
    np.testing.assert_array_equal(U, rs.fibonacci_hemisphere())
    # spread evenly: every direction of the hemisphere lies within 3 degrees of a normal, up to sign
    probe = rotation(1) @ np.random.default_rng(2).normal(size=(3, 500))
    probe /= np.linalg.norm(probe, axis=0)
    assert np.degrees(np.arccos(np.abs(U @ probe).max(0).min())) < 3.0


@pytest.mark.parametrize("turn", [None, 3, 4])
def test_the_midline_is_the_plane_that_mirrors_the_earlier_rows_onto_their_labels(turn):
    C, stage, labels, groups, _, u = bilateral_cloud(5, turn=turn)
    early = stage == 0
    found, contrast, scores = rs.midline(C[early], labels[early], C.mean(0))
    assert abs(found @ u) > 0.99 and contrast < 0.5
    assert scores.shape == (2000,) and float(scores.min()) / float(np.median(scores)) == contrast
    # a round cloud whose labels are mixed has no preferred plane: the best normal is about as good as the median one
    rng = np.random.default_rng(6)
    ball = rng.normal(size=(600, 3))
    _, flat, _ = rs.midline(ball, rng.choice(["a", "b"], 600), ball.mean(0))
    assert flat > 0.85


def test_the_score_by_hand_on_a_tiny_cloud():
    """Two labels of 20 rows: the score of a normal is the median over labels of the median distance from a mirrored row
    to its label's nearest row, the row itself included."""
    rng = np.random.default_rng(7)
    P = rng.normal(size=(40, 3))
    lab = np.array(["a"] * 20 + ["b"] * 20)
    c = P.mean(0)
    U = rs.fibonacci_hemisphere(50)
    got = rs.plane_scores(P, lab, c, U)
    for k in (0, 17, 49):
        per = []
        for name in ("a", "b"):
            Q = P[lab == name]
            M = rs.mirror(Q, c, U[k])
            per.append(np.median(np.sqrt(((M[:, None, :] - Q[None, :, :]) ** 2).sum(-1)).min(1)))
        assert got[k] == pytest.approx(np.median(per), rel=1e-12)
    with pytest.raises(rs.RunSortRefused, match="no label holds"):
        rs.plane_scores(P[:30], np.array(["a"] * 15 + ["b"] * 15), c, U)


def test_the_reassignment_keeps_every_coordinate_and_every_rows_expression():
    C, stage, labels, groups, X, _ = bilateral_cloud(11, turn=2)
    C32 = C.astype(np.float32)
    before = C32.tobytes()
    for mode in ("run_sorted", "run_sorted_same"):
        out = rs.sort_rows(C32, stage, labels, groups, mode)
        assert C32.tobytes() == before  # the step never touches the coordinates it reads
        assert sorted(out.src.tolist()) == list(range(len(C)))  # a permutation of the rows
        Y = X[out.src]
        # every row's expression is unchanged: the rows as a multiset, and each moved row is one whole row of X
        np.testing.assert_array_equal(np.sort(Y.view(np.dtype((np.void, Y.dtype.itemsize * Y.shape[1]))).ravel()),
                                      np.sort(X.view(np.dtype((np.void, X.dtype.itemsize * X.shape[1]))).ravel()))
        moved = np.flatnonzero(out.src != np.arange(len(C)))
        assert len(moved) > 0 and all((Y[i] == X[out.src[i]]).all() for i in moved)
        # each position keeps its pool, and labels and groups move with their rows
        pools = rs.pool_index(groups)
        np.testing.assert_array_equal(pools[out.src], pools)
        assert set(map(tuple, out.pairs.tolist())) == {(int(i), int(out.src[i])) for i in out.pairs[:, 0]}
    same = rs.sort_rows(C, stage, labels, groups, "as_built")
    np.testing.assert_array_equal(same.src, np.arange(len(C)))
    assert len(same.pairs) == 0


def test_rs_sends_each_bracket_to_its_own_side_and_raises_the_own_run_share():
    C, stage, labels, groups, _, _ = bilateral_cloud(12, turn=8)
    out = rs.sort_rows(C, stage, labels, groups, "run_sorted", later_contrast=True)
    s = (C - out.centre) @ out.normal
    after = out.stage_after
    np.testing.assert_array_equal(after, stage[out.src])
    for i, j in out.pairs:  # i on the later side held an earlier row, j on the earlier side a later row
        assert s[i] >= 0 > s[j] and stage[i] == 0 and stage[j] == 1 and after[i] == 1 and after[j] == 0
    pools = rs.pool_index(groups)
    for p, pool in enumerate(rs.POOLS):
        name, here = pool[0], pools == p
        wrong_before = min(((stage == 0) & (s >= 0) & here).sum(), ((stage == 1) & (s < 0) & here).sum())
        assert out.pools[name]["rs_swaps"] == out.pools[name]["swaps"] == wrong_before
        # afterwards one bracket has no row left on the wrong side in this pool (the leftovers are the other's)
        assert min(((after == 0) & (s >= 0) & here).sum(), ((after == 1) & (s < 0) & here).sum()) == 0
    o, r = rs.bracket_mix(C, stage), rs.bracket_mix(C, after)
    assert r["own_run_share"] > 0.85 > 0.6 > o["own_run_share"] and r["sd_f"] > 1.5 * o["sd_f"]
    assert out.record["contrast_earlier"] < 0.7 and out.record["contrast_later"] < 0.85
    assert out.record["median_homology_error"] < rs.reach(C)
    np.testing.assert_allclose(out.errors, np.linalg.norm(C[out.pairs[:, 1]] - rs.mirror(C[out.pairs[:, 0]],
                                                                                          out.centre, out.normal),
                                                           axis=1))
    assert rs.reach(C) == tb.neighbourhood_radius(C, 15)


def test_the_pairs_are_the_least_cost_assignment_on_the_mirror_cost():
    """Within a pool the swaps pair the wrong-side rows of the two brackets at the least total |x_j - mirror(x_i)|^2:
    checked against every assignment on a pool small enough to enumerate (4 wrong-side earlier rows, 5 later)."""
    rng = np.random.default_rng(13)
    C = rng.normal(size=(48, 3))
    C[:, 0] = np.abs(C[:, 0]) + 0.1
    C[24:] = C[:24] * [-1.0, 1.0, 1.0]  # rows 24-47 mirror rows 0-23: the centroid's x is 0
    stage = np.ones(48, int)
    stage[[0, 1, 2, 3]] = 0  # four earlier rows on x > 0 (the later side): wrong
    stage[24:43] = 0  # nineteen earlier rows on x < 0 (their side); rows 43-47 are later rows on x < 0: wrong
    labels = np.array(["a"] * 48)
    groups = np.array(["Blood"] * 48)
    plane = np.array([[1.0, 0.0, 0.0]])
    out = rs.sort_rows(C, stage, labels, groups, "run_sorted", normals=plane)
    A, B = np.array([0, 1, 2, 3]), np.arange(43, 48)
    assert len(out.pairs) == 4 and set(out.pairs[:, 0]) == set(A) and set(out.pairs[:, 1]) <= set(B)
    M = rs.mirror(C[A], out.centre, out.normal)
    cost = ((M[:, None, :] - C[B][None, :, :]) ** 2).sum(-1)
    best = min(itertools.permutations(range(5), 4), key=lambda p: sum(cost[i, p[i]] for i in range(4)))
    want = {(int(A[i]), int(B[best[i]])) for i in range(4)}
    assert set(map(tuple, out.pairs.tolist())) == want
    assert (out.errors ** 2).sum() == pytest.approx(sum(cost[i, best[i]] for i in range(4)), rel=1e-12)
    r, c = linear_sum_assignment(cost)
    assert {(int(A[i]), int(B[j])) for i, j in zip(r, c, strict=True)} == want


def test_the_placebo_matches_the_number_of_swaps_per_pool_and_side_and_leaves_the_arrangement():
    C, stage, labels, groups, _, _ = bilateral_cloud(14, turn=9, groups=("Neural", "CM", "Blood", "CM"))
    arm = rs.sort_rows(C, stage, labels, groups, "run_sorted")
    pl = rs.sort_rows(C, stage, labels, groups, "run_sorted_same")
    np.testing.assert_array_equal(pl.normal, arm.normal)  # the same plane
    s = (C - pl.centre) @ pl.normal
    for name, t in pl.pools.items():
        assert t["swaps"] == arm.pools[name]["swaps"] == t["rs_swaps"]
        assert t["placebo_earlier"] - t["placebo_later"] in (0, 1) and t["placebo_shortfall"] == 0
    pools = rs.pool_index(groups)
    for i, j in pl.pairs:  # two rows of one bracket and one pool, one on each side of the plane
        assert stage[i] == stage[j] and pools[i] == pools[j] and s[i] >= 0 > s[j]
    # per side: as many rows cross from each side as under RS
    assert (s[pl.pairs[:, 0]] >= 0).sum() == (s[arm.pairs[:, 0]] >= 0).sum() == len(arm.pairs)
    np.testing.assert_array_equal(pl.stage_after, stage)  # every position keeps its bracket: f is O's exactly
    assert rs.bracket_mix(C, pl.stage_after) == rs.bracket_mix(C, stage)
    assert pl.record["placebo_matches_swaps"] is True and pl.record["placebo_shortfall"] == 0
    assert pl.record["moved_share"] == arm.record["moved_share"]
    # each bracket's swaps are the cheapest pairs of its own assignment across the plane
    for b in (0, 1):
        plus = np.flatnonzero((stage == b) & (pools == 0) & (s >= 0))
        minus = np.flatnonzero((stage == b) & (pools == 0) & (s < 0))
        M = rs.mirror(C[plus], pl.centre, pl.normal)
        cost = ((M[:, None, :] - C[minus][None, :, :]) ** 2).sum(-1)
        r, c = linear_sum_assignment(cost)
        q = pl.pools["EpiblastStreak"]["placebo_earlier" if b == 0 else "placebo_later"]
        want = sorted(cost[r, c])[:q]
        got = sorted(float(((rs.mirror(C[[i]], pl.centre, pl.normal) - C[[j]]) ** 2).sum())
                     for i, j in pl.pairs if stage[i] == b and pools[i] == 0)
        np.testing.assert_allclose(got, want, rtol=1e-12)


def test_the_placebo_tops_up_from_the_other_bracket_and_records_a_shortfall():
    """A pool whose earlier bracket sits on one side only: its swaps all go to the later bracket. A pool where neither
    bracket has rows on both sides: no placebo swap is possible and the shortfall is recorded."""
    rng = np.random.default_rng(15)
    n = 200
    C = rng.normal(size=(n, 3))
    C[:, 0] = np.abs(C[:, 0]) + 0.2
    C[100:] = C[:100] * [-1.0, 1.0, 1.0]  # rows 0-99 at x > 0 (the later side), rows 100-199 their mirror images
    labels = np.array(["a", "b"] * (n // 2))
    groups = np.array(["Neural"] * n)
    stage = np.zeros(n, int)
    stage[100:150] = 1  # 50 later rows on x < 0: wrong; 90 earlier rows on x > 0: wrong
    stage[50:60] = 1  # 10 later rows on x > 0, so the later bracket has only 10 pairs across the plane
    plane = np.array([[1.0, 0.0, 0.0]])
    arm = rs.sort_rows(C, stage, labels, groups, "run_sorted", normals=plane)
    pl = rs.sort_rows(C, stage, labels, groups, "run_sorted_same", normals=plane)
    t = pl.pools["EpiblastStreak"]
    assert t["rs_swaps"] == arm.pools["EpiblastStreak"]["swaps"] == 50
    # ceil(50/2) = 25 each, but the later bracket has 10 pairs: the earlier one takes the other 40
    assert (t["placebo_earlier"], t["placebo_later"], t["placebo_shortfall"], t["swaps"]) == (40, 10, 0, 50)
    assert pl.record["placebo_matches_swaps"] is True
    # now no bracket has rows on both sides of the plane
    stage2 = np.r_[np.zeros(100, int), np.ones(100, int)]
    arm2 = rs.sort_rows(C, stage2, labels, groups, "run_sorted", normals=plane)
    pl2 = rs.sort_rows(C, stage2, labels, groups, "run_sorted_same", normals=plane)
    assert arm2.pools["EpiblastStreak"]["swaps"] == 100 and pl2.pools["EpiblastStreak"]["swaps"] == 0
    assert pl2.pools["EpiblastStreak"]["placebo_shortfall"] == 100
    assert pl2.record["placebo_matches_swaps"] is False and pl2.record["placebo_shortfall"] == 100


def test_the_step_is_deterministic_and_refuses_what_the_row_does_not_define():
    C, stage, labels, groups, _, _ = bilateral_cloud(16)
    a = rs.sort_rows(C, stage, labels, groups, "run_sorted")
    b = rs.sort_rows(C.copy(), stage.copy(), labels.copy(), groups.copy(), "run_sorted")
    np.testing.assert_array_equal(a.src, b.src)
    np.testing.assert_array_equal(a.normal, b.normal)
    with pytest.raises(rs.RunSortRefused, match="mode"):
        rs.sort_rows(C, stage, labels, groups, "shuffled")
    with pytest.raises(rs.RunSortRefused, match="no pool"):
        rs.sort_rows(C, stage, labels, np.where(groups == "CM", "Mesoderm", groups), "run_sorted")
    with pytest.raises(rs.RunSortRefused, match="stage"):
        rs.sort_rows(C, stage + 1, labels, groups, "run_sorted")
    with pytest.raises(rs.RunSortRefused, match="both brackets"):
        rs.sort_rows(C, np.zeros_like(stage), labels, groups, "run_sorted")
    bad = C.copy()
    bad[3, 1] = np.nan
    with pytest.raises(rs.RunSortRefused, match="finite"):
        rs.sort_rows(bad, stage, labels, groups, "run_sorted")
    with pytest.raises(rs.RunSortRefused, match="no label holds"):
        rs.sort_rows(C, stage, np.arange(len(C)).astype(str), groups, "run_sorted")  # every label one row


def test_the_mix_figures_are_n15s():
    C, stage, *_ = bilateral_cloud(17)
    from scipy.spatial import cKDTree

    nb = cKDTree(C).query(C, k=15)[1]
    f = stage[nb].mean(1)
    mix = rs.bracket_mix(C, stage)
    assert mix["sd_f"] == pytest.approx(float(f.std()), rel=1e-12)
    assert mix["own_run_share"] == pytest.approx(float(np.where(stage == 1, f, 1 - f).mean()), rel=1e-12)


# ------------------------------------------------------------------------------------------------ synthetic stages
def stage_file(name: str, seed: int, world=None) -> ad.AnnData:
    """A bracket as a released file holds it after the duplicate filter: Poisson counts over the panel and three genes
    outside it, X = log(1 + 1e4 c / N), obs celltype and group (lineage-map groups), each label two mirror-image blobs
    about x = 0."""
    rng = np.random.default_rng(seed)
    world = WORLD if world is None else world
    depth = DEPTH[name]
    C, labs, groups, P = [], [], [], []
    for k, (lab, (grp, n75, n80, centre)) in enumerate(world.items()):
        n = n80 if name == "E8.0" else n75
        prof = np.r_[np.ones(8), np.full(4, 0.3), np.ones(3)]
        prof[[k, (k + 5) % 12]] += 2.0
        if name == "E8.0":
            prof[:4] *= 2.0
        c = rng.poisson(prof / prof.sum() * depth, size=(n, len(prof)))
        c[c[:, :len(PANEL)].sum(1) == 0, 0] += 1
        C.append(c)
        labs += [lab] * n
        groups += [grp] * n
        side = np.where(rng.random(n) < 0.5, -1.0, 1.0)
        P.append(np.asarray(centre) * np.c_[side, np.ones(n), np.ones(n)] + rng.normal(size=(n, 3)) * 0.4)
    order = rng.permutation(sum(len(c) for c in C))
    C, P = np.vstack(C)[order], np.vstack(P)[order]
    labs, groups = np.asarray(labs)[order], np.asarray(groups)[order]
    N = C.sum(1)
    a = ad.AnnData(X=sparse.csr_matrix(np.log1p(1e4 * C / N[:, None]).astype(np.float32)),
                   obs=pd.DataFrame({"celltype": labs, "group": groups}, index=[f"{name}_{i}" for i in range(len(C))]),
                   var=pd.DataFrame(index=PANEL + EXTRA))
    a.layers["counts"] = sparse.csr_matrix(C.astype(np.int32))
    a.obsm["spatial_3D"] = np.asarray(P * (1.0 if name == "E7.5" else 1.3), dtype=np.float32)
    return a


@pytest.fixture(scope="module")
def stages():
    return {"E7.5.h5ad": stage_file("E7.5", 21), "E8.0.h5ad": stage_file("E8.0", 22)}


# ------------------------------------------------------------------------------------------------ the harness arms
VIEW = harness.TaskView(key="syn:E57", task="T2", kind="interp", setting="embryo", test_board="T2:embryo:test_interp",
                        roles=("before", "after"), role_days={"before": 7.5, "after": 8.0},
                        role_stems={"before": "E7.5", "after": "E8.0"}, reference_role="before", target_day=7.75)


@pytest.fixture
def ctx(stages, monkeypatch):
    monkeypatch.setattr(tb, "interpolate_clouds", partial(tb.interpolate_clouds, n_anchor=150))
    board = Board(key="pseudo:synthetic_embryo", task="T2", genes=list(PANEL), needs_coords=True, setting="embryo",
                  min_cells=10, max_cells=SECRET_MAX)
    inputs = {"before": stages["E7.5.h5ad"][:, PANEL].copy(), "after": stages["E8.0.h5ad"][:, PANEL].copy()}
    return harness.BuildContext(task=VIEW, board=board, inputs=inputs, seed=harness.BUILD_SEED_OFFSET + 2)


def _bytes(a) -> bytes:
    return np.ascontiguousarray(a).tobytes()


def _rows(X) -> np.ndarray:
    X = np.ascontiguousarray(X)
    return X[np.lexsort(X.T[::-1])]


def test_the_harness_arms_are_o_with_only_the_expression_rows_reassigned(ctx):
    o = harness.CANDIDATES["fgw_shape0x_ov"].build(ctx)
    again = harness._fgw_shape0x_ov_rs("as_built")(ctx)  # _overlay_mixture repeats O's path word for word
    assert _bytes(again.X) == _bytes(o.X) and _bytes(again.coords) == _bytes(o.coords)
    for name in harness.RS_CANDIDATES:
        arm = harness.CANDIDATES[name].build(ctx)
        assert arm.coords.dtype == o.coords.dtype and _bytes(arm.coords) == _bytes(o.coords)  # every coordinate
        assert arm.X.shape == o.X.shape and _bytes(_rows(arm.X)) == _bytes(_rows(o.X))  # every row's expression
        assert not np.array_equal(arm.X, o.X)  # rows moved
        c = harness.CANDIDATES[name]
        assert c.params["neighbourhood"] == harness.RS_CANDIDATES[name] and c.params["shape_form"] == "overlay"
        assert c.requires is harness.CANDIDATES["fgw_shape0x_ov"].requires and c.tasks == frozenset({"T2"})
    assert harness.RS_CANDIDATES == {"fgw_shape0x_ov_rs": "run_sorted", "fgw_shape0x_ov_rs_same": "run_sorted_same"}


def test_the_harness_step_reads_each_rows_own_label_and_group(ctx):
    before, after, res, X = harness._overlay_mixture(ctx)
    step = harness.rs_step(before, after, res, "run_sorted")
    lab = tb.drawn_labels(res, before.obs["celltype"], after.obs["celltype"])
    grp = tb.drawn_labels(res, before.obs["group"], after.obs["group"])
    want = rs.sort_rows(res.coords, res.stage, lab, grp, "run_sorted")
    np.testing.assert_array_equal(step.src, want.src)
    assert len(step.pairs) > 0 and step.record["contrast_earlier"] < 1.0


# ------------------------------------------------------------------------------------------------ the builder
def _board(key="T2:embryo:test_interp", split="test_interp") -> Board:
    return Board(key=key, task="T2", genes=list(PANEL), needs_coords=True, split=split, setting="embryo",
                 min_cells=50, max_cells=SECRET_MAX, truth_cells=SECRET_TRUTH)


def install(mod, setattr_, setitem_, shared: Path, reads: list, stages: dict, ledger_rows=(RULING_ROW, R48_ROW)):
    """Point a build_test module (this checkout's or BASE's) at the synthetic stages, test and validation boards, size
    JSON and ledger; the interpolation runs with 150 anchors (the code path is the builder's)."""
    boards = {b.key: b for b in (_board(), _board("T2:embryo:val_interp", "val_interp"))}

    def load(path, genes=None, drop_duplicates=False):
        name = Path(path).name
        assert Path(path).parent.name == "T2" and genes == PANEL and drop_duplicates, path
        reads.append(name)
        return stages[name].copy()

    def read(path, genes=None):
        reads.append("released:" + Path(path).name)
        a = stages[Path(path).name]
        return (a if genes is None else a[:, list(genes)]).copy()

    shared.mkdir(parents=True, exist_ok=True)
    size = shared / "size_embryo.json"
    size.write_text(json.dumps({"board": "embryo", "valid": True,
                                "preregistered": {"z_mode": "raw", "valid": True, "pred_rms_full": 6.0}}))
    ledger = shared / "TRIES.md"
    ledger.write_text("| ID | what |\n|---|---|\n" + "".join(r + "\n" for r in ledger_rows))
    setattr_(mod, "load_stage", load)
    setattr_(mod, "read_stage", read)
    setattr_(mod, "load_boards", lambda: boards)
    setattr_(mod, "_t2_count", lambda b, n: {"decision": "synthetic"})
    setattr_(mod.t2_build, "interpolate_clouds", partial(mod.t2_build.interpolate_clouds, n_anchor=150))
    setitem_(mod.SIZES, "embryo", size)
    setattr_(mod, "TRIES_FILE", ledger)


@pytest.fixture(scope="module")
def bt():
    sys.path.insert(0, str(REPO / "experiments" / "builds"))
    import build_test

    return build_test


@pytest.fixture
def env(bt, monkeypatch, own_repo, tmp_path, stages):
    reads: list[str] = []
    install(bt, monkeypatch.setattr, monkeypatch.setitem, tmp_path / "shared", reads, stages)
    return reads


def _shows(text: str, *values: int) -> bool:
    """Whether any of the integers appears in the text as a number of its own (not as digits inside a float)."""
    return any(re.search(rf"(?<![\d.]){v}(?![\d.])", text) for v in values)


def _build(mod, tmp_path: Path, name: str, *argv: str) -> tuple[Path, dict, str]:
    out = tmp_path / name / "T2embryo.h5ad"
    mod.main(["embryo", *argv, "--out", str(out)])
    text = out.with_suffix(".build.json").read_text()
    return out, json.loads(text), text


def _file(path: Path):
    a = ad.read_h5ad(path)
    X = np.asarray(a.X.toarray() if sparse.issparse(a.X) else a.X)
    return X, np.asarray(a.obsm["spatial_3D"]), list(a.obs_names), list(a.var_names)


def test_on_os_command_the_file_is_os_but_for_where_its_expression_rows_sit(bt, env, tmp_path, capsys):
    o, rec_o, _ = _build(bt, tmp_path, "o", *O_ARGS)
    w, rec_w, text = _build(bt, tmp_path, "rs", *RS_ARGS)
    Xo, Co, no, vo = _file(o)
    Xw, Cw, nw, vw = _file(w)
    assert Cw.dtype == Co.dtype and Cw.tobytes() == Co.tobytes()  # O's positions, bit for bit, in O's order
    assert nw == no and vw == vo and len(Xw) == len(Xo)  # O's names, genes and count
    assert _bytes(_rows(Xw)) == _bytes(_rows(Xo))  # every row's expression, unchanged
    assert (np.abs(Xw - Xo).max(1) > 0).mean() > 0.2  # rows moved
    assert rec_w["record"] == rec_o["record"] and rec_w["n_cells"] == rec_o["n_cells"] == "max_cells"
    assert rec_w["shape_interp"] == rec_o["shape_interp"]
    blk = rec_w["neighbourhood"]
    assert blk["captain_ruling"] == RULING and blk["file"] == "T2embryo_s2x_ovrs.h5ad"
    assert blk["command"] == bt.EMBRYO_RS_COMMAND.format(row=RULING) == "embryo " + " ".join(RS_ARGS)
    assert blk["o_command"] == "embryo " + " ".join(O_ARGS) and blk["rule"] == rs.RULE and "E57" in blk["read"]
    assert set(blk) == {"option", "file", "captain_ruling", "command", "o_command", "rule", "read", "normal_x",
                        "normal_y", "normal_z", "contrast_earlier"}
    assert all(not isinstance(v, (dict, list)) for v in blk.values())
    assert "neighbourhood" not in json.dumps(rec_o)
    assert not _shows(text + capsys.readouterr().out, SECRET_MAX, SECRET_TRUTH)
    w2, _, _ = _build(bt, tmp_path, "rs2", *RS_ARGS)  # built again: the same bytes
    assert w2.read_bytes() == w.read_bytes()


def test_the_option_runs_after_the_cloud_guard_and_moves_rows_only(bt, env, tmp_path, monkeypatch):
    calls = []
    inner_guard, inner_step = bt.t2_build.cloud_guard, bt.harness.rs_step

    def guard(*a, **kw):
        calls.append("guard")
        return inner_guard(*a, **kw)

    def step(before, after, res, mode):
        calls.append(("step", mode))
        coords = np.array(res.coords, copy=True)
        out = inner_step(before, after, res, mode)
        assert np.array_equal(coords, res.coords)
        return out

    monkeypatch.setattr(bt.t2_build, "cloud_guard", guard)
    monkeypatch.setattr(bt.harness, "rs_step", step)
    _build(bt, tmp_path, "rs", *RS_ARGS)
    assert calls == ["guard", ("step", "run_sorted")]
    calls.clear()
    _build(bt, tmp_path, "o", *O_ARGS)
    assert calls == ["guard"]


def _refused(bt, env, tmp_path, argv, match, name="T2embryo.h5ad"):
    out = tmp_path / "r" / name
    with pytest.raises(SystemExit, match=match) as e:
        bt.main(["embryo", *argv, "--out", str(out)])
    assert not out.exists() and env == []  # refused before any stage was read
    assert not _shows(str(e.value), SECRET_MAX, SECRET_TRUTH)
    return str(e.value)


@pytest.mark.parametrize("argv", [
    ["--method", "fgw", "--neighbourhood", "run_sorted", "--captain-ruling", RULING],
    ["--method", "fgw_shape0", *OV, "--neighbourhood", "run_sorted", "--captain-ruling", RULING],
    [*FLOOR, "--neighbourhood", "run_sorted", "--captain-ruling", RULING],  # no form
    [*FLOOR, "--shape-interp", "overlay", "--neighbourhood", "run_sorted", "--captain-ruling", RULING],  # no yes flag
    [*FLOOR, "--shape-interp", "stretch_ratio", "--shape-interp-test-ok", "--neighbourhood", "run_sorted",
     "--captain-ruling", RULING],
    [*RS_ARGS, "--expression", "series"],
    [*RS_ARGS, "--expression", "render_hi"],
    [*RS_ARGS, "--row-order", "copy_parity", "--row-order-test-ok"],
    [*RS_ARGS, "--herd"],
    [*RS_ARGS, "--draw", "midpoint"],
    [*RS_ARGS, "--move-coupling", "aligned"],
    [*RS_ARGS, "--copy-rows", "--copy-rows-test-ok"],
    [*RS_ARGS, "--seed", "1000002"],
])
def test_any_command_but_os_is_refused(bt, env, tmp_path, argv):
    _refused(bt, env, tmp_path, argv, "O's command")


def test_no_dry_run_and_never_with_series_back(bt, env, tmp_path):
    _refused(bt, env, tmp_path, [*RS_ARGS, "--board", "T2:embryo:val_interp"], "O's command", "T2embryo_dry.h5ad")
    _refused(bt, env, tmp_path, [*FLOOR, "--expression", "series_back", *OV, "--captain-ruling", "R48",
                                 "--neighbourhood", "run_sorted"], "O's command")


def test_the_ruling_is_needed_and_must_name_the_option(bt, env, tmp_path):
    _refused(bt, env, tmp_path, [*O_ARGS, "--neighbourhood", "run_sorted"], "needs --captain-ruling")
    for row in ("R48", "R61", "E57", "R60; true", "r60"):
        _refused(bt, env, tmp_path, [*O_ARGS, "--neighbourhood", "run_sorted", "--captain-ruling", row],
                 "no single row")


@pytest.mark.parametrize("ledger", [
    (R48_ROW,),  # no row R60
    (RULING_ROW, RULING_ROW),  # two rows R60
    (RULING_ROW.replace("`--neighbourhood run_sorted`", "--neighbourhood run_sorted"),),  # not in backticks
    (RULING_ROW.replace(f"`--captain-ruling {RULING}`", f"--captain-ruling {RULING}"),),
    (RULING_ROW.replace("ruling, no run", "ruling"),),
    (RULING_ROW.replace("(the captain's ruling", "(the AI assistant's ruling"),),
    (f"| {RULING} (the captain's ruling, ruling, no run) | `--neighbourhood run_sorted` "
     f"`--captain-ruling {RULING}` |",),
    (RULING_ROW.replace(f"| {RULING} ", f"| {RULING}a "),),
])
def test_the_ledger_must_hold_exactly_one_ruling_row_naming_the_option(bt, stages, monkeypatch, own_repo, tmp_path,
                                                                     ledger):
    reads: list[str] = []
    install(bt, monkeypatch.setattr, monkeypatch.setitem, tmp_path / "shared", reads, stages, ledger_rows=ledger)
    _refused(bt, reads, tmp_path, RS_ARGS, f"no single row {RULING}")


def test_the_ruling_check_is_the_builders_existing_one(bt, tmp_path):
    src = (REPO / "experiments" / "builds" / "build_test.py").read_text()
    assert "_captain_ruling_row(ruling, EMBRYO_RS_OPTION)" in src
    assert bt.EMBRYO_RS_OPTION == "`--neighbourhood run_sorted`"
    p = tmp_path / "TRIES.md"
    p.write_text("| ID | what |\n|---|---|\n" + RULING_ROW + "\n")
    assert bt._captain_ruling_row(RULING, bt.EMBRYO_RS_OPTION, p) is True
    assert bt._captain_ruling_row("R48", bt.EMBRYO_RS_OPTION, p) is False


def test_the_real_ledger_holds_ruling_row_r57_as_the_builder_needs_it(bt):
    """Row R57 on main, the captain's ruling that orders the file: the builder's own check passes on its first cell,
    read from the real experiments/harness/TRIES.md by code (a boolean; nothing of the row is printed), and the
    option's refusals take the file's command with --captain-ruling R57 before anything is read."""
    assert bt.TRIES_FILE == REPO / "experiments" / "harness" / "TRIES.md"
    assert bt._captain_ruling_row("R57", bt.EMBRYO_RS_OPTION) is True
    a = bt.build_parser().parse_args(["embryo", *O_ARGS, "--expression", "mixture", "--neighbourhood", "run_sorted",
                                      "--captain-ruling", "R57", "--out", "x.h5ad"])
    assert bt._run_sorted_refusals(a) == "R57"
    assert bt._captain_ruling_row("R48", bt.EMBRYO_RS_OPTION) is False  # R48 orders another file


def test_a_ruling_alone_is_refused_as_before(bt, env, tmp_path):
    for argv in ([*O_ARGS, "--captain-ruling", RULING], ["--method", "fgw", "--captain-ruling", RULING]):
        _refused(bt, env, tmp_path, argv, "--captain-ruling goes with --expression series_back only")


def test_the_option_is_off_by_default_named_in_the_parser_and_never_the_placebo(bt):
    a = bt.build_parser().parse_args(["embryo", *O_ARGS, "--out", "x.h5ad"])
    assert a.neighbourhood == "as_built" and a.captain_ruling is None
    assert bt.NEIGHBOURHOODS == ("as_built", "run_sorted") and "run_sorted_same" not in bt.NEIGHBOURHOODS
    with pytest.raises(SystemExit):
        bt.build_parser().parse_args(["embryo", *O_ARGS, "--neighbourhood", "run_sorted_same", "--out", "x.h5ad"])
    src = (REPO / "experiments" / "builds" / "build_test.py").read_text()
    assert "rs_ruling = _run_sorted_refusals(a)" in src and src.count('default="mixture"') == 1
    assert bt.EMBRYO_RS_FILE == "T2embryo_s2x_ovrs.h5ad"


# ------------------------------------------------------------------------------------------------ every other command
def _module(name: str, rel: str) -> types.ModuleType:
    """A module made from BASE's text of `rel`, registered under `name`."""
    mod = types.ModuleType(name)
    mod.__file__ = str(REPO / rel)
    sys.modules[name] = mod
    exec(compile(git(REPO, "show", f"{BASE}:{rel}") + "\n", str(REPO / rel), "exec"), mod.__dict__)
    return mod


@pytest.fixture(scope="module")
def base_bt(bt):
    """BASE's builder (main before the option); every module it imports is this checkout's, so a difference below is
    the option's."""
    mod = _module("_build_test_before_rs", "experiments/builds/build_test.py")
    assert not hasattr(mod, "EMBRYO_RS_FILE") and "neighbourhood" not in mod.__doc__
    return mod


DRY = "T2:embryo:val_interp"
# name: (argv, what main's builder does with it in this synthetic world: "built" or "refused")
COMMANDS = {
    "O": (O_ARGS, "built"),
    "W": (W_ARGS, "built"),
    "SO": ([*FLOOR, "--expression", "series", *OV], "built"),
    "O_hi": ([*FLOOR, "--expression", "render_hi", *OV], "built"),
    "P": (FLOOR, "built"),
    "s2": (["--method", "fgw_shape0"], "built"),
    "s1": (["--method", "fgw"], "built"),
    "stretch": ([*FLOOR, "--shape-interp", "stretch_ratio", "--shape-interp-test-ok"], "built"),
    "seed": ([*O_ARGS, "--seed", "1000003"], "built"),
    "copy_parity": ([*FLOOR, "--row-order", "copy_parity", "--row-order-test-ok"], "refused"),
    "ruling_alone": ([*O_ARGS, "--captain-ruling", "R48"], "refused"),
    "w_wrong_row": ([*FLOOR, "--expression", "series_back", *OV, "--captain-ruling", "R47"], "refused"),
    "form_needs_ok": ([*FLOOR, "--shape-interp", "overlay"], "refused"),
    "seed_low": ([*O_ARGS, "--seed", "5"], "refused"),
}


def _outcome(mod, tmp: Path, stages, which: str, argv: list[str], name: str):
    """What one builder does with one command: the h5ad's sha256 and its record's text, or ("refused", its message)."""
    with pytest.MonkeyPatch.context() as mp:
        install(mod, mp.setattr, mp.setitem, tmp / "shared", [], stages)
        out = tmp / which / name
        try:
            mod.main(["embryo", *argv, "--out", str(out)])
        except SystemExit as e:
            return "refused", str(e.code).replace(str(out), "<out>")
        return hashlib.sha256(out.read_bytes()).hexdigest(), out.with_suffix(".build.json").read_text()


@pytest.mark.parametrize("key", list(COMMANDS))
def test_every_other_command_does_what_main_did(bt, base_bt, own_repo, tmp_path, stages, key):
    argv, kind = COMMANDS[key]
    want = _outcome(base_bt, tmp_path, stages, "base", argv, "T2embryo.h5ad")
    got = _outcome(bt, tmp_path, stages, "new", argv, "T2embryo.h5ad")
    assert ("refused" if want[0] == "refused" else "built") == kind
    assert got[0] == want[0] and got[1] == want[1]  # the same bytes and record, or the same refusal


def test_a_dry_run_of_o_does_what_main_did(bt, base_bt, own_repo, tmp_path, stages):
    argv = [*FLOOR, "--shape-interp", "overlay", "--board", DRY]
    want = _outcome(base_bt, tmp_path, stages, "base", argv, "T2embryo_dry.h5ad")
    got = _outcome(bt, tmp_path, stages, "new", argv, "T2embryo_dry.h5ad")
    assert want[0] != "refused" and got == want
