"""Embryo arrangement batch (campaign/emb-arr): fgw_balanced, fgw_patchy and fgw_shuffle, fgw with one re-arrangement
of its drawn expression rows over its own positions (vec.t2_build.balanced_bracket_labels, patchy_bracket_labels,
relabel_rows, span_groups). Synthetic data only; nothing scored, no target read."""

from __future__ import annotations

import hashlib
from dataclasses import replace

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from conftest import GENES, simulate_stage

from vec import harness, t2_build
from vec.boards import Board
from vec.harness import CANDIDATES, NBHD_K, PseudoTask, load_tasks

ARMS = ("fgw_balanced", "fgw_patchy", "fgw_shuffle")
BRACKET_ARMS = ("fgw_balanced", "fgw_patchy")
BOTH = ("T2e/E7.25+E8.0->E7.5", "T2e/E6.75+E7.5->E7.25")


def _pool(n, rng, means, props, scale=1.0):
    centres = np.array([[0, 0, 0], [6, 0, 0], [0, 6, 0]], float) * scale
    st = simulate_stage(n, means, props, rng, coords_scale=1.5 * scale, coords_shift=centres)
    a = ad.AnnData(X=st.X, obs=pd.DataFrame({"celltype": st.celltype, "group": st.celltype},
                                            index=[f"c{i}" for i in range(n)]), var=pd.DataFrame(index=GENES))
    a.obsm["spatial_3D"] = st.coords.astype(np.float32)
    return a


def _task(**over) -> PseudoTask:
    base = PseudoTask(key="syn_arr", test_board="T2:embryo:test_interp", task="T2", setting="embryo", kind="interp",
                      target="T2/E7.5.h5ad", reference="T2/E7.25.h5ad",
                      inputs={"before": "T2/E7.25.h5ad", "after": "T2/E8.0.h5ad"}, twin_board="T2:embryo:val_interp")
    return replace(base, **over)


def _make_ctx(seed: int = 57, n_after: int = 900):
    rng = np.random.default_rng(seed)
    base = rng.gamma(2.0, 2.0, size=(3, len(GENES)))
    before = _pool(600, rng, base, [0.5, 0.3, 0.2])
    after = _pool(n_after, rng, base * 1.4, [0.3, 0.4, 0.3], scale=1.5)
    board = Board(key="pseudo:synthetic_embryo", task="T2", genes=list(GENES), needs_coords=True, setting="embryo",
                  min_cells=100, max_cells=500)
    return harness.BuildContext(task=harness.TaskView.of(_task()), board=board,
                                inputs={"before": before, "after": after}, seed=harness.BUILD_SEED_OFFSET + 4)


@pytest.fixture(scope="module")
def ctx():
    return _make_ctx()


def _sha(st) -> str:
    h = hashlib.sha256(np.ascontiguousarray(st.X, dtype=np.float32).tobytes())
    h.update(np.ascontiguousarray(st.coords, dtype=np.float32).tobytes())
    return h.hexdigest()


def _res(ctx, **kw):
    """fgw's own interpolation of `ctx` (the call harness._fgw makes), for each output cell's bracket."""
    before, after = ctx.inputs["before"], ctx.inputs["after"]
    g = ctx.board.genes
    Eb, Ea = (harness.dense(a[:, g].X).astype(np.float32) for a in (before, after))
    t = ctx.task
    frac = (t.target_day - t.role_days["before"]) / (t.role_days["after"] - t.role_days["before"])
    return harness.t2_build.interpolate_clouds(np.asarray(before.obsm["spatial_3D"], float)[:, :3],
                                               np.asarray(after.obsm["spatial_3D"], float)[:, :3], Eb, Ea, frac,
                                               harness.count_for(ctx.board, before.n_obs), ctx.seed,
                                               target_rms=harness.target_rms(ctx), **kw)


def _sorted_rows(X: np.ndarray) -> np.ndarray:
    return X[np.lexsort(X.T[::-1])]


def _share_spread(coords: np.ndarray, stage: np.ndarray, k: int = NBHD_K) -> float:
    """The SD over positions of the later bracket's share among each position's k nearest positions."""
    near = t2_build.cKDTree(coords).query(coords, k=k)[1]
    return float(stage[near].mean(1).std())


# ------------------------------------------------------------------------------------------- existing paths unchanged
# sha256 of X (float32 bytes) then coords (float32 bytes) of each candidate on `_make_ctx()`, taken on main (f63f964)
# before the arrangement step was added: the default path of `_fgw` is bit for bit what it was
PRE_EDIT_SHA256 = {
    "fgw": "37f9a5d4993a1ba9b22a5bcc522f42326f5054652dc102a2a9e706f39a0d5ac7",
    "fgw_aligned": "28d87f05ea8be86b498a9078635729b6a09ed6bfa1ba4375bb65344a28f2b8b5",
    "fgw_herd": "fabc738f35842cc8d511e44b6ff8d262eae425268576e8b0de2bb2647249fed5",
    "fgw_perm": "5779cf13d73024fa01c2aeb8799a4d20a70a4c3fb00c6fb61a2f1e096e4ac772",
    "fgw_shape0": "1eab0ac7232d4f747db90e9290daf57d9e898f40121549033395a8d88326d960",
    "fgw_late": "7c98eb721d88e76c15524ac0f4d9aacf2fd36762ee7ce276d738923c2af1af86",
}


@pytest.mark.parametrize("name", sorted(PRE_EDIT_SHA256))
def test_existing_fgw_candidates_are_bit_for_bit_as_before(ctx, name):
    assert _sha(CANDIDATES[name].build(ctx)) == PRE_EDIT_SHA256[name]


# ------------------------------------------------------------------------------------------------------- registration
def test_registration_matches_fgw_aligned():
    tasks = [t for ts in load_tasks().values() for t in ts]
    keys = {t.key for t in tasks}
    assert set(BOTH) <= keys
    for name in ARMS:
        c = CANDIDATES[name]
        assert c.name == name and c.tasks == CANDIDATES["fgw_aligned"].tasks == frozenset({"T2"})
        assert [c.applies(t) for t in tasks] == [CANDIDATES["fgw_aligned"].applies(t) for t in tasks]
        assert all(c.applies(t) for t in tasks if t.key in BOTH)
        assert not c.tuned_on and c.params["scope"] == "whole"
        assert c.params["stream"] == harness.ARRANGE_STREAMS[name.removeprefix("fgw_")]
    assert CANDIDATES["fgw_balanced"].params["k"] == CANDIDATES["fgw_shuffle"].params["k"] == NBHD_K == 15
    assert CANDIDATES["fgw_patchy"].params["kmeans_seed"] == t2_build.MIDPOINT_KMEANS_SEED
    streams = list(harness.ARRANGE_STREAMS.values())
    assert len(set(streams)) == 3 and harness.NBHD_STREAM not in streams and min(streams) > 9
    assert "fgw_potts" not in CANDIDATES  # blocked: no update rule on record


def test_refusals(ctx):
    with pytest.raises(ValueError, match="arrange must be one of"):
        harness._fgw(ctx, arrange="potts")
    with pytest.raises(ValueError, match="use one of them"):
        harness._fgw(ctx, arrange="shuffle", local_perm=True)
    with pytest.raises(ValueError, match="mode must be one of"):
        harness.bracket_arrangement(np.zeros((3, 3)), np.zeros(3, int), "potts", 0)
    C = np.random.default_rng(0).normal(size=(10, 3))
    s = np.r_[np.zeros(5, int), np.ones(5, int)]
    with pytest.raises(ValueError, match="number"):
        t2_build.relabel_rows(C, s, np.r_[np.zeros(4, int), np.ones(6, int)])
    with pytest.raises(ValueError, match="brackets 0"):
        t2_build.balanced_bracket_labels(C, s + 1, 15, np.random.default_rng(0))
    with pytest.raises(ValueError, match="one entry per position"):
        t2_build.patchy_bracket_labels(C, s[:-1], np.random.default_rng(0))


# ------------------------------------------------------------------------------- each arm differs only by its idea
@pytest.mark.parametrize("name", ARMS)
def test_arm_keeps_fgws_cells_positions_and_counts_and_moves_rows_only(ctx, name):
    base, arm = CANDIDATES["fgw"].build(ctx), CANDIDATES[name].build(ctx)
    assert arm.n == base.n == ctx.board.max_cells
    np.testing.assert_array_equal(arm.coords, base.coords)  # the same positions, in the same order
    np.testing.assert_array_equal(_sorted_rows(arm.X), _sorted_rows(base.X))  # the same drawn cells
    np.testing.assert_array_equal(np.sort(arm.X.sum(1)), np.sort(base.X.sum(1)))  # every cell's total kept
    assert not np.array_equal(arm.X, base.X)  # something moved
    res = _res(ctx)
    at = harness.bracket_arrangement(res.coords, res.stage, name.removeprefix("fgw_"), ctx.seed)
    assert np.array_equal(np.sort(at), np.arange(base.n))
    np.testing.assert_array_equal(arm.X, base.X[at])  # exactly the arrangement step, nothing else
    assert (res.stage[at] == 1).sum() == (res.stage == 1).sum()


def test_balanced_spreads_and_patchy_gathers_the_brackets(ctx):
    res = _res(ctx)
    C, s = res.coords, res.stage
    spread = _share_spread(C, s)
    for mode, cmp in (("balanced", np.less), ("patchy", np.greater)):
        at = harness.bracket_arrangement(C, s, mode, ctx.seed)
        assert cmp(_share_spread(C, s[at]), spread), mode
    # the shuffle ignores the brackets: every row ends within its own group of NBHD_K positions
    group = t2_build.span_groups(C, NBHD_K, np.random.default_rng([ctx.seed, harness.ARRANGE_STREAMS["shuffle"]]))
    at = harness.bracket_arrangement(C, s, "shuffle", ctx.seed)
    np.testing.assert_array_equal(group[at], group)
    sizes = np.bincount(group)
    assert sizes.max() == NBHD_K and (sizes == NBHD_K).mean() > 0.9


def test_bracket_arrangement_keeps_rows_of_unchanged_positions_and_swaps_the_rest():
    rng = np.random.default_rng(3)
    C = rng.normal(size=(400, 3))
    s = (rng.random(400) < 1 / 3).astype(np.int64)
    for mode in ("balanced", "patchy"):
        lab = (t2_build.balanced_bracket_labels(C, s, NBHD_K, np.random.default_rng([5, 21])) if mode == "balanced"
               else t2_build.patchy_bracket_labels(C, s, np.random.default_rng([5, 22])))
        assert (lab == 1).sum() == (s == 1).sum()
        at = t2_build.relabel_rows(C, s, lab)
        np.testing.assert_array_equal(s[at], lab)
        same = lab == s
        np.testing.assert_array_equal(at[same], np.flatnonzero(same))
        np.testing.assert_array_equal(at[at], np.arange(400))  # pairwise swaps
        np.testing.assert_array_equal(at, harness.bracket_arrangement(C, s, mode, 5))
    assert t2_build.patch_count(s) == 3 and t2_build.patch_count(np.r_[0, 1]) == 2
    lab = t2_build.patchy_bracket_labels(C, s, np.random.default_rng(1))
    from sklearn.cluster import KMeans

    patch = KMeans(n_clusters=3, init="k-means++", n_init=1, algorithm="lloyd",
                   random_state=t2_build.PATCH_KMEANS_SEED).fit(C).labels_
    mixed = [j for j in range(3) if 0 < lab[patch == j].mean() < 1]
    assert len(mixed) <= 1  # one bracket per patch, one patch split to hold the number exactly


# ------------------------------------------------------------------------- nothing to do: the arm is fgw bit for bit
def test_with_one_bracket_drawn_the_bracket_arms_are_fgw(ctx):
    base = harness._fgw(ctx, t_mix=1.0)
    for mode in ("balanced", "patchy"):
        assert _sha(harness._fgw(ctx, t_mix=1.0, arrange=mode)) == _sha(base)
    C = np.random.default_rng(0).normal(size=(50, 3))
    for s in (np.zeros(50, int), np.ones(50, int)):
        for mode in ("balanced", "patchy"):
            np.testing.assert_array_equal(harness.bracket_arrangement(C, s, mode, 1), np.arange(50))


def test_with_groups_of_one_the_shuffle_is_fgw(ctx, monkeypatch):
    C = np.random.default_rng(0).normal(size=(50, 3))
    np.testing.assert_array_equal(t2_build.span_shuffle(C, 1, np.random.default_rng(0)), np.arange(50))
    monkeypatch.setattr(harness, "NBHD_K", 1)
    assert _sha(CANDIDATES["fgw_shuffle"].build(ctx)) == PRE_EDIT_SHA256["fgw"]


def test_an_already_balanced_draw_is_left_as_it_is():
    """On a draw the greedy pass would itself produce, the balanced arm returns it unchanged."""
    C = np.random.default_rng(4).normal(size=(200, 3))
    s0 = (np.arange(200) % 3 == 0).astype(np.int64)
    lab = t2_build.balanced_bracket_labels(C, s0, NBHD_K, np.random.default_rng([9, 21]))
    again = t2_build.balanced_bracket_labels(C, lab, NBHD_K, np.random.default_rng([9, 21]))
    np.testing.assert_array_equal(again, lab)
    np.testing.assert_array_equal(harness.bracket_arrangement(C, lab, "balanced", 9), np.arange(200))
