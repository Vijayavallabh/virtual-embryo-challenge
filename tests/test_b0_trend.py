"""Heart candidate `b0_trend` (row H11 of experiments/harness/TRIES.md; both heart files by the captain's ruling of
2026-11-13, docs/PLAN.md decision log (w); ported from campaign/trend 8bb23bf): B0 with each lineage's per-gene velocity
the least-squares trend of log(mean CP10k + eps) over every input stage, read over B0's time budget; and the hook it
uses, `vec.t1_scope.momentum_shift(..., logf_by_type=...)`. Synthetic data only; nothing scored."""

from __future__ import annotations

import hashlib
import inspect
import sys
from dataclasses import replace
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from vec import core, harness, t1_scope
from vec.boards import Board
from vec.harness import B0_DAMPING, CANDIDATES, PseudoTask, load_tasks

# ---------------------------------------------------------------------------------------------------- momentum_shift


def _adata(C: np.ndarray, types: list[str], genes: list[str], prefix: str) -> ad.AnnData:
    a = ad.AnnData(X=np.log1p(1e4 * C / C.sum(1, keepdims=True).clip(1)).astype(np.float32),
                   obs=pd.DataFrame({"b0_type": types}, index=[f"{prefix}{i}" for i in range(len(C))]),
                   var=pd.DataFrame(index=genes))
    a.layers["counts"] = sparse.csr_matrix(C.astype(np.int32))
    return a


@pytest.fixture(scope="module")
def shift_inputs():
    """Types A and B at both stages; C only at the reference (skipped for too few start cells)."""
    rng = np.random.default_rng(21)
    genes = [f"g{i}" for i in range(50)]
    mu = rng.gamma(2.0, 2.0, 50)
    start = _adata(rng.poisson(mu, (300, 50)), ["A"] * 200 + ["B"] * 100, genes, "s")
    ref = _adata(rng.poisson(mu * np.linspace(0.5, 2.0, 50), (400, 50)), ["A"] * 250 + ["B"] * 120 + ["C"] * 30,
                 genes, "r")
    return start, ref, ref[np.arange(0, 400, 3)].copy(), genes


# sha256 of the output counts (int32 dense bytes) of `momentum_shift(start, ref, like, 0.5, key="b0_type",
# seed=[7, 3])` on `shift_inputs`, computed with the function as it was before `logf_by_type` was added
# (campaign/gen, 459916f; the same function as main's before this port)
PRE_CHANGE_SHA256 = "82b17858ce44b34cc83e4c92ae539f3cd34adee3f19c54f470a9a7565e2b6437"


def _sha(a: ad.AnnData) -> str:
    return hashlib.sha256(np.ascontiguousarray(a.layers["counts"].toarray(), dtype=np.int32).tobytes()).hexdigest()


def test_momentum_shift_without_the_mapping_is_unchanged(shift_inputs):
    """1. Default None is the function as before (pinned hash); an empty mapping, and a mapping that hands each type
    B0's own factor, give the same counts bit for bit (same stream, same factors)."""
    start, ref, like, _ = shift_inputs
    base, rep = t1_scope.momentum_shift(start, ref, like, 0.5, key="b0_type", seed=[7, 3])
    assert _sha(base) == PRE_CHANGE_SHA256
    empty, _ = t1_scope.momentum_shift(start, ref, like, 0.5, key="b0_type", seed=[7, 3], logf_by_type={})
    S, R = core.counts_of(start), core.counts_of(ref)
    own = {}
    for t in ("A", "B"):
        s, r = start.obs["b0_type"] == t, ref.obs["b0_type"] == t
        own[t] = 0.5 * (np.log(core.mean_cp10k(R[r.to_numpy()]) + 0.5) - np.log(core.mean_cp10k(S[s.to_numpy()]) + 0.5))
    same, rep2 = t1_scope.momentum_shift(start, ref, like, 0.5, key="b0_type", seed=[7, 3], logf_by_type=own)
    for out in (empty, same):
        np.testing.assert_array_equal(base.layers["counts"].toarray(), out.layers["counts"].toarray())
    assert rep.rows == rep2.rows


def test_the_mapping_replaces_only_named_types_and_never_revives_a_skipped_one(shift_inputs):
    """2. A's factor replaced by zero (its cells keep their counts), B keeps its own factor (its report row is B0's),
    and C, skipped for too few start cells, stays skipped although the mapping names it."""
    start, ref, like, genes = shift_inputs
    base, rep = t1_scope.momentum_shift(start, ref, like, 0.5, key="b0_type", seed=[7, 3])
    bump = np.linspace(-1.0, 1.0, len(genes))
    out, rep2 = t1_scope.momentum_shift(start, ref, like, 0.5, key="b0_type", seed=[7, 3],
                                        logf_by_type={"A": np.zeros(len(genes)), "C": bump})
    rows, rows2 = {r["type"]: r for r in rep.rows}, {r["type"]: r for r in rep2.rows}
    assert rows2["A"]["status"] == "shifted" and rows2["A"]["mean_abs_log_f"] == 0.0
    assert rows["A"]["mean_abs_log_f"] > 0
    assert {k: v for k, v in rows2["B"].items() if k != "median_depth_scale"} == {
        k: v for k, v in rows["B"].items() if k != "median_depth_scale"}
    assert rows2["C"]["status"].startswith("skipped") and rows2["C"] == rows["C"]
    t = like.obs["b0_type"].to_numpy()
    old, new = like.layers["counts"].toarray(), out.layers["counts"].toarray()
    np.testing.assert_array_equal(old[t == "A"], new[t == "A"])
    np.testing.assert_array_equal(old[t == "C"], new[t == "C"])
    assert not np.array_equal(old[t == "B"], new[t == "B"])
    with pytest.raises(ValueError, match="shape"):
        t1_scope.momentum_shift(start, ref, like, 0.5, key="b0_type", seed=[7, 3], logf_by_type={"A": np.zeros(3)})


# ---------------------------------------------------------------------------------------------------- heart fixtures
# (as campaign/trend's tests/test_h7_gmom.py, which main does not have)
GENES = [f"g{i}" for i in range(120)]
LINEAGES = np.array(["CM", "Endothelial", "Untracked", "Tiny"])


def _stage(n: int, seed: int, shift: float, tiny: int) -> ad.AnnData:
    """A T2 stage: CM and Endothelial lineages with a within-lineage state axis (a module of genes that a cell's state
    score moves), Untracked cells, and `tiny` cells of a lineage too small for the minimum at this stage. `shift`
    moves the CM states' module along the axis between stages, so the within-lineage change differs by state."""
    rng = np.random.default_rng(seed)
    types = rng.choice(3, size=n, p=[0.5, 0.35, 0.15])
    types[:tiny] = 3
    pi = np.random.default_rng(1).gamma(1.0, 1.0, len(GENES))
    pi /= pi.sum()
    a = np.zeros(len(GENES))
    a[:30] = np.random.default_rng(4).choice([-1.0, 1.0], 30) * 1.5
    s = rng.normal(size=n)
    lam = pi * np.exp(np.outer(s, a) + np.outer(shift * s * (types == 0), a) + 0.3 * (types == 1)[:, None])
    depth = np.maximum(np.round(np.exp(rng.normal(np.log(80.0), 0.4, n))), 5)[:, None]
    C = rng.poisson(lam / lam.sum(1, keepdims=True) * depth).astype(np.int32)
    C[C.sum(1) == 0, 0] = 1
    X = np.log1p(1e4 * C / C.sum(1, keepdims=True)).astype(np.float32)
    centres = np.array([[0.0, 0.0, 0.0], [5.0, 0.0, 0.0], [0.0, 5.0, 0.0], [5.0, 5.0, 0.0]])
    out = ad.AnnData(X=sparse.csr_matrix(X), obs=pd.DataFrame(
        {"celltype": [f"type{t}" for t in types], "merged": LINEAGES[types]},
        index=[f"c{seed}_{i}" for i in range(n)]), var=pd.DataFrame(index=GENES))
    out.layers["counts"] = sparse.csr_matrix(C)
    out.obsm["spatial_3D"] = (centres[types] + rng.normal(size=(n, 3)) * 1.5).astype(np.float32)
    return out


def _heart_task(**over) -> PseudoTask:
    kw = {"key": "syn_h7", "test_board": "T2:heart:test_extrap", "task": "T2", "kind": "extrap",
          "target": "T2/E10.5.h5ad", "reference": "T2/E9.5.h5ad", "setting": "heart",
          "inputs": {"ref": "T2/E9.5.h5ad", "start": "T2/E8.75.h5ad"}, "validation_board": "T2:heart:val_extrap"}
    return PseudoTask(**{**kw, **over})


BOARD = Board(key="pseudo:h7", task="T2", genes=list(GENES), needs_coords=True, setting="heart", min_cells=100,
              max_cells=700)

# ---------------------------------------------------------------------------------------------------- the trend
FLOW_INPUTS = {"ref": "T2/E9.5.h5ad", "start": "T2/E8.75.h5ad", "earlier": "T2/E8.5.h5ad"}
DAYS = {"earlier": 8.5, "start": 8.75, "ref": 9.5}


def _view(inputs=FLOW_INPUTS, **over):
    return harness.TaskView.of(_heart_task(inputs=inputs, **over))


def _ctx(inputs: dict, view=None, seed: int = harness.BUILD_SEED_OFFSET + 1) -> harness.BuildContext:
    return harness.BuildContext(task=view or _view(), board=BOARD, inputs=inputs, seed=seed)


def _trend(ctx: harness.BuildContext):
    """b0_trend's output and the record `_b0_t2` hands back (the build record's source)."""
    rec: dict = {}
    return harness._b0_t2(ctx, trend=True, trend_record=rec), rec


@pytest.fixture(scope="module")
def three():
    """Three heart-like stages: CM and Endothelial usable at all three, Tiny below min_cells at start (not shifted),
    Untracked cells."""
    return {"ref": _stage(1600, seed=2, shift=0.0, tiny=60), "start": _stage(1100, seed=3, shift=0.8, tiny=10),
            "earlier": _stage(900, seed=5, shift=1.6, tiny=10)}


def test_with_no_lineage_usable_before_start_b0_trend_is_b0_bit_for_bit(three):
    """3. The earlier stage holds no lineage at min_cells (its cells are relabelled to a lineage absent elsewhere and
    Untracked): every shifted lineage keeps B0's own factor, and the output is b0's bit for bit."""
    early = three["earlier"].copy()
    m = np.where(np.arange(early.n_obs) % 2 == 0, "Elsewhere", "Untracked")
    m[:5] = "CM"  # a few CM cells, below min_cells
    early.obs["merged"] = m
    inputs = {**three, "earlier": early}
    a, rec = _trend(_ctx(inputs))
    b = CANDIDATES["b0"].build(_ctx(inputs))
    np.testing.assert_array_equal(a.X, b.X)
    np.testing.assert_array_equal(a.coords, b.coords)
    np.testing.assert_array_equal(CANDIDATES["b0_trend"].build(_ctx(inputs)).X, b.X)
    assert rec["n_trend"] == 0
    assert {r["lineage"] for r in rec["lineages"]} == {"CM", "Endothelial"}
    assert all(r["own_factor"] and r["stages"] == ["start", "ref"] and r["weights"] == [1.0] for r in rec["lineages"])
    assert rec["median_cosine"] == pytest.approx(1.0) and rec["median_norm_ratio"] == pytest.approx(1.0)


@pytest.mark.parametrize("days", [[8.5, 8.75, 9.5], [8.3, 8.5, 8.75, 9.5]])
def test_a_linear_trend_gives_b0s_factor(days):
    """4. Pseudobulks with log(mean CP10k + eps) exactly linear in the day per gene: the trend factor is B0's (the last
    step of a linear trend is the trend). Exact linearity for every gene cannot be built from counts (CP10k sums to
    1e4 at every stage), so the stage pseudobulks are built directly."""
    rng = np.random.default_rng(3)
    a, b = rng.uniform(2.0, 4.0, 80), rng.uniform(-1.5, 1.5, 80)  # a + b (d - 9.5) > log(eps): means stay positive
    eps = 0.5
    means = [np.exp(a + b * (d - 9.5)) - eps for d in days]
    assert all((m > 0).all() for m in means)
    y = [np.log(m + eps) for m in means]
    budget = B0_DAMPING * (days[-1] - days[-2])
    trend = harness.trend_log_factor(y, days, budget)
    own = B0_DAMPING * (y[-1] - y[-2])
    np.testing.assert_allclose(trend, own, atol=1e-9, rtol=0)


def test_a_linear_trend_through_the_candidate_gives_b0s_factor(three, monkeypatch):
    """4 (wiring). The candidate's own path, with each stage's pseudobulk set to a linear-in-day design (vec.core.
    mean_cp10k replaced by the design, keyed on the stage's cell number): the mapping it hands momentum_shift equals
    B0_DAMPING x the last step to 1e-9, for each lineage with a trend."""
    rng = np.random.default_rng(4)
    G = len(BOARD.genes)
    coef = {lin: (rng.uniform(2.0, 4.0, G), rng.uniform(-1.5, 1.5, G)) for lin in ("CM", "Endothelial")}
    design_of = {}  # a stage's lineage cell number -> (lineage, day)
    for r, s in three.items():
        for lin in coef:
            n = int((np.asarray(s.obs["merged"]) == lin).sum())
            assert n not in design_of
            design_of[n] = (lin, DAYS[r])
    real = core.mean_cp10k

    def design(C):
        if C.shape[0] not in design_of:
            return real(C)
        lin, d = design_of[C.shape[0]]
        a, b = coef[lin]
        return np.exp(a + b * (d - 9.5)) - 0.5

    monkeypatch.setattr(core, "mean_cp10k", design)
    ctx = _ctx(three)
    out, rec = harness.b0_trend_factors(ctx, list(BOARD.genes), ["CM", "Endothelial", "Tiny", "Untracked@end"])
    assert set(out) == {"CM", "Endothelial"}
    for lin, logf in out.items():
        np.testing.assert_allclose(logf, B0_DAMPING * coef[lin][1] * (9.5 - 8.75), atol=1e-9, rtol=0)
        row = next(r for r in rec["lineages"] if r["lineage"] == lin)
        assert row["stages"] == ["earlier", "start", "ref"] and row["cosine"] == pytest.approx(1.0)


@pytest.mark.parametrize("days", [[8.5, 8.75, 9.5], [8.3, 8.5, 8.75, 9.5]])
def test_a_zig_zag_at_start_is_damped_by_the_trend(days):
    """5. A planted per-stage offset added to `start`'s log means (a middle stage: the second to last) on a linear
    trend: the trend factor is closer to the planted true slope (over B0's budget) than B0's factor is."""
    rng = np.random.default_rng(5)
    a, b = rng.normal(2.0, 1.0, 80), rng.normal(0.0, 1.5, 80)
    y = [a + b * d for d in days]
    y[-2] = y[-2] + rng.normal(0.0, 0.4, 80)
    budget = B0_DAMPING * (days[-1] - days[-2])
    truth = budget * b
    trend = harness.trend_log_factor(y, days, budget)
    own = B0_DAMPING * (y[-1] - y[-2])
    assert np.linalg.norm(trend - truth) < np.linalg.norm(own - truth)


def test_transition_weights():
    """6. The record's transition weights for gaps (0.2, 0.25, 0.75) and (0.2, 0.25), and for the test build's five
    stages (E8.25_late at 8.3, E8.5, E8.75, E9.5, E10.5); they sum to 1."""
    w3 = harness.trend_transition_weights(np.cumsum([0.0, 0.2, 0.25, 0.75]))
    w2 = harness.trend_transition_weights(np.cumsum([0.0, 0.2, 0.25]))
    w4 = harness.trend_transition_weights([8.3, 8.5, 8.75, 9.5, 10.5])
    np.testing.assert_allclose(np.round(w3, 4), [0.1119, 0.2192, 0.6689])
    np.testing.assert_allclose(np.round(w2, 4), [0.4262, 0.5738])
    np.testing.assert_allclose(np.round(w4, 4), [0.0500, 0.1095, 0.4118, 0.4287])
    assert w3.sum() == pytest.approx(1.0) and w2.sum() == pytest.approx(1.0) and w4.sum() == pytest.approx(1.0)
    assert harness.trend_transition_weights([8.75, 9.5]).tolist() == [1.0]


def test_trend_stages_reads_roles_and_days_only_and_refuses_bad_geometry():
    """The stages in day order from the view's roles and days; fewer than three, a reference that is not the latest,
    `start` not second to last, or equal days are refused."""
    v = _view(inputs={"ref": "T2/E9.5.h5ad", "start": "T2/E8.75.h5ad", "earlier": "T2/E8.5.h5ad",
                      "earliest": "T2/E8.25_late.h5ad"})
    assert harness.trend_stages(v) == [("earliest", 8.3), ("earlier", 8.5), ("start", 8.75), ("ref", 9.5)]
    assert harness.trend_stages(replace(v, target_day=float("nan"))) == harness.trend_stages(v)
    with pytest.raises(ValueError, match="three"):
        harness.trend_stages(_view(inputs={"ref": "T2/E9.5.h5ad", "start": "T2/E8.75.h5ad"}))
    with pytest.raises(ValueError, match="latest"):
        harness.trend_stages(replace(v, role_days={**v.role_days, "earliest": 9.6}))
    with pytest.raises(ValueError, match="latest"):
        harness.trend_stages(replace(v, role_days={**v.role_days, "earlier": 8.9}))
    with pytest.raises(ValueError, match="increase"):
        harness.trend_stages(replace(v, role_days={**v.role_days, "earlier": 8.75}))


def test_the_target_is_never_read_and_other_tasks_are_refused(three):
    """7. Built with the target's day NaN and no target among the inputs, under the read guard: nothing refused, and
    the target's day changes nothing. Two input stages, an interpolation task and T1 are refused, and so is the trend
    with depth, herd or cfm, or a record asked for without it."""
    blind = replace(_view(), target_day=float("nan"))
    assert set(three) == set(FLOW_INPUTS) and "T2/E10.5.h5ad" not in FLOW_INPUTS.values()
    with harness.guarded_reads() as g:
        a, rec = _trend(_ctx(three, view=blind))
    assert not g["refused"] and np.isfinite(a.X).all()
    b = CANDIDATES["b0_trend"].build(_ctx(three))
    np.testing.assert_array_equal(a.X, b.X)
    assert rec["stages"] == ["earlier", "start", "ref"] and rec["budget_days"] == pytest.approx(0.375)
    assert rec["n_trend"] == 2 and {r["lineage"] for r in rec["lineages"]} == {"CM", "Endothelial"}  # not Tiny
    two = {k: three[k] for k in ("ref", "start")}
    with pytest.raises(ValueError, match="three"):
        CANDIDATES["b0_trend"].build(_ctx(two, view=_view(inputs={"ref": "T2/E9.5.h5ad", "start": "T2/E8.75.h5ad"})))
    interp = replace(_view(), kind="interp")
    with pytest.raises(ValueError, match="heart extrapolation"):
        CANDIDATES["b0_trend"].build(_ctx(three, view=interp))
    with pytest.raises(ValueError, match="heart extrapolation"):
        CANDIDATES["b0_trend"].build(_ctx(three, view=replace(_view(), task="T1", setting=None)))
    with pytest.raises(ValueError, match="heart extrapolation"):
        CANDIDATES["b0_trend"].build(_ctx(three, view=replace(_view(), setting="embryo")))
    for kw in ({"depth": True}, {"herd": True}, {"cfm": True}):
        with pytest.raises(ValueError, match="depth, herd or cfm"):
            harness._b0_t2(_ctx(three), trend=True, **kw)
    with pytest.raises(ValueError, match="trend_record"):
        harness._b0_t2(_ctx(three), trend_record={})


def test_registration():
    """8. Applies to the two heart backtests only; b0 and b0_trend share a scope, so `decide` judges them alike."""
    tasks = {t.key: t for ts_ in load_tasks().values() for t in ts_}
    c = CANDIDATES["b0_trend"]
    assert {k for k, t in tasks.items() if c.applies(t)} == {"T2h/<=E9.5->E10.5", "T2h/<=E8.75->E10.5"}
    for k in ("T1/E9.5->E10.5", "T1/E9.5->E10.5:chamber", "T1/E8.5->E9.5", "T2e/E7.25+E8.0->E7.5",
              "T2e/E6.75+E7.5->E7.25"):
        assert not c.applies(tasks[k])
    assert c.tasks == frozenset({"T2"}) and c.derived == () and c.requires is harness._heart_extrap_with_start
    assert c.params["scope"] == CANDIDATES["b0"].params["scope"] == "chamber"
    assert c.params["damping"] == B0_DAMPING and c.params["rule"] == "vec.harness.B0_TREND_RULE"
    assert "b0_trend_sized" not in CANDIDATES


def test_cells_coordinates_and_rows_are_b0s(three, monkeypatch):
    """9. The same seed draws the same rows; coordinates and cell number equal b0's; Untracked cells are held flat in
    both, and a lineage with a trend moves differently."""
    seen = []
    real = harness.draw_copy

    def spy(ctx, herd=False):
        out = real(ctx, herd=herd)
        seen.append(out[1].copy())
        return out

    monkeypatch.setattr(harness, "draw_copy", spy)
    ctx = _ctx(three)
    a, b = CANDIDATES["b0_trend"].build(ctx), CANDIDATES["b0"].build(ctx)
    assert len(seen) == 2 and np.array_equal(seen[0], seen[1])
    assert a.n == b.n == BOARD.max_cells
    np.testing.assert_array_equal(a.coords, b.coords)
    merged = np.asarray(three["ref"].obs["merged"]).astype(str)[seen[0]]
    flat = (merged == "Untracked") | (merged == "Tiny")
    assert flat.any() and np.array_equal(a.X[flat], b.X[flat])
    assert not np.array_equal(a.X[merged == "CM"], b.X[merged == "CM"])


# sha256 of b0_trend's output X (float32 bytes) on `three` at BUILD_SEED_OFFSET + 1, as campaign/trend pins it
# (computed there at 28f25fd): the port builds the scored candidate's very output
B0_TREND_SHA256 = "8320700817f003ae9b5daff16ac5375fcb8442313f3454e69a048e7bfe87855c"


def test_b0_trend_output_is_the_campaign_branchs(three):
    st = CANDIDATES["b0_trend"].build(_ctx(three))
    assert hashlib.sha256(np.ascontiguousarray(st.X, dtype=np.float32).tobytes()).hexdigest() == B0_TREND_SHA256


# ---------------------------------------------------------------------------------------------------- review follow-up
# (the main session's review of 272b6d9: a lineage usable at only some earlier stages, the min_cells boundary, the
# mapping with the clip and held-flat genes, and the record's top-level weights)
FOUR_INPUTS = {"ref": "T2/E9.5.h5ad", "start": "T2/E8.75.h5ad", "earlier": "T2/E8.5.h5ad",
               "earliest": "T2/E8.25_late.h5ad"}
FOUR_DAYS = [8.3, 8.5, 8.75, 9.5]
TYPES = ["CM", "Endothelial", "Tiny", "Untracked@end"]


@pytest.fixture(scope="module")
def four(three):
    return {**three, "earliest": _stage(1000, seed=7, shift=2.2, tiny=10)}


def _with_lineage_cells(a: ad.AnnData, lineage: str, keep: int) -> ad.AnnData:
    """`a` with only the first `keep` cells of `lineage` left in it (the others relabelled to a lineage absent
    elsewhere)."""
    out = a.copy()
    m = np.asarray(out.obs["merged"]).astype(object)
    idx = np.flatnonzero(m == lineage)
    assert len(idx) > keep
    m[idx[keep:]] = "Elsewhere"
    out.obs["merged"] = m.astype(str)
    return out


def _log_means(a: ad.AnnData, lineage: str) -> np.ndarray:
    m = np.asarray(a.obs["merged"]).astype(str) == lineage
    return np.log(core.mean_cp10k(core.counts_of(a, list(BOARD.genes))[np.flatnonzero(m)]) + 0.5)


@pytest.fixture(scope="module")
def bt():
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "experiments" / "builds"))
    import build_test

    return build_test


def test_a_lineage_short_at_one_earlier_stage_takes_the_slope_over_its_usable_stages(four, bt):
    """Endothelial falls one cell below min_cells at `earlier` only: its factor is the budget x the least-squares slope
    over earliest, start and ref; its record lists those stages and their own transition weights; CM is unaffected.
    The build record's top-level weights come from all four stage days, each lineage's from its usable days."""
    n_min = harness._momentum_default("min_cells")
    view = _view(inputs=FOUR_INPUTS)
    partial = {**four, "earlier": _with_lineage_cells(four["earlier"], "Endothelial", n_min - 1)}
    out, rec = harness.b0_trend_factors(_ctx(partial, view=view), list(BOARD.genes), TYPES)
    full, rec_full = harness.b0_trend_factors(_ctx(four, view=view), list(BOARD.genes), TYPES)
    assert set(out) == set(full) == {"CM", "Endothelial"}
    budget = B0_DAMPING * (9.5 - 8.75)
    y = [_log_means(partial[r], "Endothelial") for r in ("earliest", "start", "ref")]
    want = harness.trend_log_factor(y, [8.3, 8.75, 9.5], budget)
    np.testing.assert_allclose(out["Endothelial"], want, atol=1e-12, rtol=0)
    own = B0_DAMPING * (y[2] - y[1])
    assert not np.allclose(out["Endothelial"], own) and not np.allclose(out["Endothelial"], full["Endothelial"])
    np.testing.assert_array_equal(out["CM"], full["CM"])
    rows = {r["lineage"]: r for r in rec["lineages"]}
    assert rows["Endothelial"]["stages"] == ["earliest", "start", "ref"]
    assert rows["Endothelial"]["days"] == [8.3, 8.75, 9.5] and not rows["Endothelial"]["own_factor"]
    w3 = harness.trend_transition_weights([8.3, 8.75, 9.5]).tolist()
    w4 = harness.trend_transition_weights(FOUR_DAYS).tolist()
    assert rows["Endothelial"]["weights"] == w3 and rows["CM"]["weights"] == w4
    assert rows["CM"] == next(r for r in rec_full["lineages"] if r["lineage"] == "CM")
    assert rec["stages"] == ["earliest", "earlier", "start", "ref"] and rec["days"] == FOUR_DAYS
    block = bt._trend_record(rec)  # what the build record shows
    assert block["weights"] == w4 and block["days"] == FOUR_DAYS
    shown = {r["lineage"]: r for r in block["lineages"]}
    assert shown["Endothelial"]["weights"] == w3 and shown["CM"]["weights"] == w4
    assert shown["Endothelial"]["stages"] == ["earliest", "start", "ref"]


@pytest.mark.parametrize("delta, usable", [(0, True), (-1, False)])
def test_the_min_cells_boundary_is_momentum_shifts_own(three, delta, usable):
    """An earlier stage with exactly min_cells cells of a lineage is usable, one cell fewer is not; min_cells is
    momentum_shift's own default (the value B0 runs with), and the record says so."""
    n_min = inspect.signature(t1_scope.momentum_shift).parameters["min_cells"].default
    assert harness._momentum_default("min_cells") == n_min
    inputs = {**three, "earlier": _with_lineage_cells(three["earlier"], "Endothelial", n_min + delta)}
    out, rec = harness.b0_trend_factors(_ctx(inputs), list(BOARD.genes), TYPES)
    assert rec["min_cells"] == n_min
    row = next(r for r in rec["lineages"] if r["lineage"] == "Endothelial")
    if usable:
        assert row["stages"] == ["earlier", "start", "ref"] and "Endothelial" in out
    else:
        assert row["stages"] == ["start", "ref"] and row["own_factor"] and "Endothelial" not in out
    assert next(r for r in rec["lineages"] if r["lineage"] == "CM")["stages"] == ["earlier", "start", "ref"]


def test_a_mapped_factor_is_clipped_and_keeps_held_flat_genes_at_one(shift_inputs):
    """With logf_by_type, the clip (max_abs_log_factor) and the held-flat zeros apply to the mapped factor as to the
    own one: a mapped factor of 3 clipped at log 2 is the mapped factor log 2, and a mapped factor that is nonzero on
    held-flat genes gives what it gives with those genes zeroed (counts and report)."""
    start, ref, like, genes = shift_inputs
    G, cap = len(genes), float(np.log(2.0))
    kw = {"key": "b0_type", "seed": [7, 3]}
    big, rep_big = t1_scope.momentum_shift(start, ref, like, 0.5, max_abs_log_factor=cap,
                                           logf_by_type={"A": np.full(G, 3.0)}, **kw)
    at_cap, rep_cap = t1_scope.momentum_shift(start, ref, like, 0.5, max_abs_log_factor=cap,
                                              logf_by_type={"A": np.full(G, cap)}, **kw)
    np.testing.assert_array_equal(big.layers["counts"].toarray(), at_cap.layers["counts"].toarray())
    assert rep_big.rows == rep_cap.rows
    row_a = next(r for r in rep_big.rows if r["type"] == "A")
    assert row_a["max_abs_log_f"] == pytest.approx(cap, abs=1e-15)
    held = ["g0", "g1", "g7"]
    bump = np.linspace(-1.0, 1.0, G) + 0.05
    zeroed = bump.copy()
    zeroed[[0, 1, 7]] = 0.0
    a, rep_a = t1_scope.momentum_shift(start, ref, like, 0.5, held, logf_by_type={"A": bump}, **kw)
    b, rep_b = t1_scope.momentum_shift(start, ref, like, 0.5, held, logf_by_type={"A": zeroed}, **kw)
    np.testing.assert_array_equal(a.layers["counts"].toarray(), b.layers["counts"].toarray())
    assert rep_a.rows == rep_b.rows
    row = next(r for r in rep_a.rows if r["type"] == "A")
    assert row["mean_abs_log_f"] == pytest.approx(float(np.abs(zeroed).mean()), abs=1e-15)


# ---------------------------------------------------------------------------------------------------- H13e, H13r, H13h
# Three variants of b0_trend (heart tries H13e, H13r, H13h; specification fixed by the main session on 2026-11-13): only
# the per-lineage, per-gene log factor differs. Synthetic data only; nothing scored.
ARMS = {"b0_trend_eb": "eb", "b0_trend_ts": "ts", "b0_trend_h": "horizon"}

# taken on this branch BEFORE the arms were added (main 7717366, the synthetic stages above, BUILD_SEED_OFFSET + 1):
# sha256 of the output X (float32) and coordinates (float32) of b0 and b0_trend, and of b0_trend's record from
# `_b0_t2(ctx, trend=True, trend_record=rec)` (json.dumps(rec, sort_keys=True))
PIN_B0_X = "7f9b4113111857bdc351c027f81b80bbde24153a02e07ff99c9cbd31f764e137"
PIN_COORDS = "1704b9c9251aad9a5e6df8d0fe2a5aab8a42c4647c58af60536c86915664f9ee"
PIN_TREND_X = {"three": B0_TREND_SHA256, "four": "2f176a19017393e88078138a572d5baf5bf7dad4589613aafd45fec6609916dd"}
PIN_TREND_REC = {"three": "019f1eb8c555321479c388b125eb8c30f885dbb3c85b7ab704e5c86bd6a36543",
                 "four": "0d010c3465b80d63dfa17ed284cfe03c6b135bd621e154a30a7e3459e81adb57"}


def _h(a) -> str:
    return hashlib.sha256(np.ascontiguousarray(a, dtype=np.float32).tobytes()).hexdigest()


@pytest.mark.parametrize("case", ["three", "four"])
def test_b0_and_b0_trend_are_unchanged_by_the_arms(three, four, case):
    """1. b0's and b0_trend's outputs (counts and coordinates) and b0_trend's record through `_b0_t2(ctx, trend=True,
    trend_record=rec)` (the test builder's call) are what they were before the arms were added."""
    import json

    inputs, view = (three, _view()) if case == "three" else (four, _view(inputs=FOUR_INPUTS))
    b = CANDIDATES["b0"].build(_ctx(inputs, view=view))
    a = CANDIDATES["b0_trend"].build(_ctx(inputs, view=view))
    assert b.n == a.n == 700
    assert _h(b.X) == PIN_B0_X and _h(b.coords) == PIN_COORDS
    assert _h(a.X) == PIN_TREND_X[case] and _h(a.coords) == PIN_COORDS
    rec: dict = {}
    c = harness._b0_t2(_ctx(inputs, view=view), trend=True, trend_record=rec)
    np.testing.assert_array_equal(c.X, a.X)
    np.testing.assert_array_equal(c.coords, a.coords)
    assert hashlib.sha256(json.dumps(rec, sort_keys=True).encode()).hexdigest() == PIN_TREND_REC[case]
    assert "form" not in rec and rec["rule"] == "vec.harness.B0_TREND_RULE"
    out, _ = harness.b0_trend_factors(_ctx(inputs, view=view), list(BOARD.genes), TYPES)
    out_ols, _ = harness.b0_trend_factors(_ctx(inputs, view=view), list(BOARD.genes), TYPES, form="ols")
    assert out.keys() == out_ols.keys() and all(np.array_equal(out[k], out_ols[k]) for k in out)
    with pytest.raises(ValueError, match="trend_form"):
        harness._b0_t2(_ctx(inputs, view=view), trend_form="eb")
    with pytest.raises(ValueError, match="unknown trend form"):
        harness.b0_trend_factors(_ctx(inputs, view=view), list(BOARD.genes), TYPES, form="huber")


# ---- eb
def _eb_by_hand(Y: np.ndarray, d: list[float]) -> tuple[list[float], list[list[int]]]:
    """The H13e weights computed with plain loops: the bins (gene indices) and their weights."""
    n, G = Y.shape
    dbar = sum(d) / n
    sxx = sum((x - dbar) ** 2 for x in d)
    b, rss, lev = [], [], []
    for g in range(G):
        ybar = sum(Y[k, g] for k in range(n)) / n
        bg = sum((d[k] - dbar) * Y[k, g] for k in range(n)) / sxx
        b.append(bg)
        rss.append(sum((Y[k, g] - ybar - bg * (d[k] - dbar)) ** 2 for k in range(n)))
        lev.append(float(np.mean(Y[:, g])))
    order = sorted(range(G), key=lambda g: (lev[g], g))
    sizes = [G // 10 + (1 if q < G % 10 else 0) for q in range(10)]
    bins, at = [], 0
    for s in sizes:
        bins.append(order[at:at + s])
        at += s
    weights = []
    for q in bins:
        v = sum(rss[g] for g in q) / (len(q) * (n - 2)) / sxx
        t2 = max(0.0, sum(b[g] ** 2 for g in q) / len(q) - v)
        weights.append(t2 / (t2 + v) if t2 + v > 0 else 0.0)
    return weights, bins


@pytest.mark.parametrize("days", [[8.5, 8.75, 9.5], [8.3, 8.5, 8.75, 9.5], [8.0, 8.3, 8.5, 8.75, 9.5]])
def test_eb_weights_follow_the_formula_bin_by_bin(days):
    """2. Log means with known noise (slopes grow with the level, the noise is constant): the weights equal the formula
    recomputed by hand bin by bin, the bins are of equal size up to one, and n = 3 (df = 1) works."""
    rng = np.random.default_rng(13)
    G = 95
    level = rng.uniform(0.0, 5.0, G)
    slope = rng.normal(0.0, 1.0, G) * level / 5.0
    Y = np.vstack([level + slope * (d - 9.5) + rng.normal(0.0, 0.3, G) for d in days])
    w_gene, w_bin = harness.trend_eb_weights(list(Y), days)
    want, bins = _eb_by_hand(Y, days)
    np.testing.assert_allclose(w_bin, want, rtol=1e-9, atol=1e-12)
    for q, idx in enumerate(bins):
        np.testing.assert_array_equal(w_gene[idx], w_bin[q])
    sizes = [len(q) for q in bins]
    assert max(sizes) - min(sizes) <= 1 and sum(sizes) == G
    assert 0.0 <= w_bin.min() and w_bin.max() <= 1.0 and len(set(np.round(w_bin, 12))) > 1


def test_eb_ties_are_broken_by_gene_order():
    """2. Every gene at the same level (exact arithmetic): the bins are consecutive blocks of the gene list, of sizes
    3, 3, 3, 2, ... for 23 genes, and each block's weight is the formula's on that block."""
    G = 23
    d = [0.0, 1.0, 2.0, 3.0]
    dc = np.array([-1.5, -0.5, 0.5, 1.5])
    e = np.array([1.0, -1.0, -1.0, 1.0])  # orthogonal to 1 and to dc: residual only
    s = np.arange(G) % 4 - 1.0
    r = (np.arange(G) * 7 % 5) / 4.0
    Y = 2.0 + np.outer(dc, s) + np.outer(e, r)
    assert (Y.mean(0) == 2.0).all()
    w_gene, w_bin = harness.trend_eb_weights(list(Y), d)
    want, bins = _eb_by_hand(Y, d)
    assert bins[0] == [0, 1, 2] and bins[3] == [9, 10] and bins[-1] == [21, 22]
    assert [len(q) for q in bins] == [3, 3, 3] + [2] * 7
    np.testing.assert_allclose(w_bin, want, rtol=1e-12, atol=1e-15)
    for q, idx in enumerate(bins):
        np.testing.assert_array_equal(w_gene[idx], w_bin[q])
    assert len(set(np.round(w_bin, 12))) > 1


@pytest.mark.parametrize("days", [[8.5, 8.75, 9.5], [8.3, 8.5, 8.75, 9.5]])
def test_eb_on_an_exact_line_keeps_the_ols_factor(days):
    """2. Stage means exactly on a line: every weight is 1 (to rounding) and the factor is b0_trend's."""
    rng = np.random.default_rng(3)
    a, b = rng.uniform(2.0, 4.0, 80), rng.uniform(-1.5, 1.5, 80)
    y = [a + b * (d - 9.5) for d in days]
    w_gene, w_bin = harness.trend_eb_weights(y, days)
    np.testing.assert_allclose(w_bin, 1.0, atol=1e-9)
    budget = B0_DAMPING * (days[-1] - days[-2])
    np.testing.assert_allclose(budget * w_gene * harness.trend_log_factor(y, days, 1.0),
                               harness.trend_log_factor(y, days, budget), atol=1e-9, rtol=0)


def test_eb_on_pure_noise_about_a_flat_line_shrinks_to_near_zero():
    """2. No trend, only noise: the slopes' spread is their sampling noise, so the weights are near 0."""
    rng = np.random.default_rng(17)
    G, days = 4000, [8.3, 8.5, 8.75, 9.5, 10.0]
    level = rng.uniform(0.0, 5.0, G)
    y = [level + rng.normal(0.0, 0.5, G) for _ in days]
    _, w_bin = harness.trend_eb_weights(y, days)
    assert w_bin.max() < 0.25 and np.median(w_bin) < 0.1


def test_eb_through_the_candidate_and_its_record(three):
    """2 (wiring). On the synthetic stages: each trend lineage's factor is budget x B_q(g) x b_g from its own log
    means; the record holds the ten weights, their median, and norm ratio and cosine against b0_trend's factor."""
    ctx = _ctx(three)
    out, rec = harness.b0_trend_factors(ctx, list(BOARD.genes), TYPES, form="eb")
    ols, _ = harness.b0_trend_factors(ctx, list(BOARD.genes), TYPES)
    assert set(out) == set(ols) == {"CM", "Endothelial"}
    budget = B0_DAMPING * (9.5 - 8.75)
    rows = {r["lineage"]: r for r in rec["lineages"]}
    for lin in out:
        y = [_log_means(three[r], lin) for r in ("earlier", "start", "ref")]
        w_gene, w_bin = harness.trend_eb_weights(y, [8.5, 8.75, 9.5])
        slope = harness.trend_log_factor(y, [8.5, 8.75, 9.5], 1.0)
        np.testing.assert_allclose(out[lin], budget * w_gene * slope, atol=1e-12, rtol=0)
        np.testing.assert_allclose(ols[lin], budget * slope, atol=1e-12, rtol=0)
        row = rows[lin]
        assert row["eb_weights"] == [float(w) for w in w_bin] and len(row["eb_weights"]) == 10
        assert row["eb_median_weight"] == pytest.approx(float(np.median(w_bin)))
        assert row["norm_ratio_ols"] == pytest.approx(np.linalg.norm(out[lin]) / np.linalg.norm(ols[lin]))
        assert row["cosine_ols"] == pytest.approx(harness._cosine(out[lin], ols[lin]))
    assert rec["rule"] == "vec.harness.B0_TREND_EB_RULE" and rec["form"] == "eb"
    assert rec["n_trend"] == rec["n_mapped"] == 2 and rec["budget_days"] == pytest.approx(budget)


def _only_start_and_ref(three) -> dict:
    early = three["earlier"].copy()
    early.obs["merged"] = np.where(np.arange(early.n_obs) % 2 == 0, "Elsewhere", "Untracked")
    return {**three, "earlier": early}


def test_eb_and_ts_keep_b0s_own_factor_where_only_start_and_ref_are_usable(three):
    """A lineage usable at only start and the reference is left out of the mapping, as in b0_trend."""
    inputs = _only_start_and_ref(three)
    for form in ("eb", "ts"):
        out, rec = harness.b0_trend_factors(_ctx(inputs), list(BOARD.genes), TYPES, form=form)
        assert out == {} and rec["n_trend"] == rec["n_mapped"] == 0
        assert all(r["own_factor"] and "cosine_ols" not in r for r in rec["lineages"])


# ---- ts
def test_ts_equals_ols_on_an_exact_line():
    """3. Stage means exactly on a line: the Theil-Sen slope is the least-squares one (and the line's)."""
    rng = np.random.default_rng(3)
    days = [8.3, 8.5, 8.75, 9.5]
    a, b = rng.uniform(2.0, 4.0, 80), rng.uniform(-1.5, 1.5, 80)
    y = [a + b * (d - 9.5) for d in days]
    np.testing.assert_allclose(harness.theil_sen_slope(y, days), harness.trend_log_factor(y, days, 1.0),
                               atol=1e-9, rtol=0)
    np.testing.assert_allclose(harness.theil_sen_slope(y, days), b, atol=1e-9, rtol=0)


def test_ts_ignores_one_outlying_stage_where_ols_does_not():
    """3. Five stages on a line, one (the second) shifted by a constant per-gene offset: the Theil-Sen factor stays at
    the clean slope (six of ten pair slopes are clean, so both middle values are) while the OLS factor moves by the
    offset's leverage, (d_1 - dbar) / Sxx x offset x budget."""
    rng = np.random.default_rng(23)
    days = [8.0, 8.3, 8.5, 8.75, 9.5]
    a, b = rng.uniform(2.0, 4.0, 80), rng.uniform(-1.5, 1.5, 80)
    y = [a + b * (d - 9.5) for d in days]
    off = rng.normal(0.0, 1.0, 80)
    y[1] = y[1] + off
    budget = B0_DAMPING * (9.5 - 8.75)
    ts = budget * harness.theil_sen_slope(y, days)
    ols = harness.trend_log_factor(y, days, budget)
    clean = budget * b
    assert np.abs(ts - clean).max() < 1e-9
    dc = np.asarray(days) - np.mean(days)
    np.testing.assert_allclose(ols - clean, budget * dc[1] / (dc ** 2).sum() * off, atol=1e-9, rtol=0)
    assert np.linalg.norm(ols - clean) > 0.1 * np.linalg.norm(clean)


def test_ts_with_three_stages_is_the_median_of_the_three_pair_slopes():
    rng = np.random.default_rng(29)
    days = [8.5, 8.75, 9.5]
    y = [rng.normal(2.0, 1.0, 50) for _ in days]
    pairs = np.vstack([(y[1] - y[0]) / 0.25, (y[2] - y[0]) / 1.0, (y[2] - y[1]) / 0.75])
    np.testing.assert_allclose(harness.theil_sen_slope(y, days), np.sort(pairs, axis=0)[1], atol=1e-12, rtol=0)


def test_ts_through_the_candidate_and_its_record(three):
    """3 (wiring). Each trend lineage's factor is the budget x its Theil-Sen slope; the record holds the cosine and
    norm ratio against b0_trend's factor."""
    ctx = _ctx(three)
    out, rec = harness.b0_trend_factors(ctx, list(BOARD.genes), TYPES, form="ts")
    ols, _ = harness.b0_trend_factors(ctx, list(BOARD.genes), TYPES)
    budget = B0_DAMPING * (9.5 - 8.75)
    rows = {r["lineage"]: r for r in rec["lineages"]}
    assert set(out) == {"CM", "Endothelial"}
    for lin in out:
        y = [_log_means(three[r], lin) for r in ("earlier", "start", "ref")]
        np.testing.assert_allclose(out[lin], budget * harness.theil_sen_slope(y, [8.5, 8.75, 9.5]), atol=1e-12,
                                   rtol=0)
        assert rows[lin]["cosine_ols"] == pytest.approx(harness._cosine(out[lin], ols[lin]))
        assert rows[lin]["norm_ratio_ols"] == pytest.approx(np.linalg.norm(out[lin]) / np.linalg.norm(ols[lin]))
    assert rec["rule"] == "vec.harness.B0_TREND_TS_RULE" and rec["form"] == "ts"


# ---- horizon
def test_horizon_scales_b0_trends_factor_by_the_horizon(three):
    """4. A trend lineage's factor is b0_trend's x (target_day - day_ref) / (day_ref - day_start); the record holds
    budget_h, b0_trend's budget and their ratio."""
    ctx = _ctx(three)  # target E10.5: target_day 10.5
    out, rec = harness.b0_trend_factors(ctx, list(BOARD.genes), TYPES, form="horizon")
    ols, _ = harness.b0_trend_factors(ctx, list(BOARD.genes), TYPES)
    scale = (10.5 - 9.5) / (9.5 - 8.75)
    assert set(out) == set(ols) == {"CM", "Endothelial"}
    for lin in out:
        np.testing.assert_allclose(out[lin], ols[lin] * scale, atol=1e-12, rtol=1e-12)
    assert rec["budget_h"] == pytest.approx(B0_DAMPING * 1.0)
    assert rec["budget_b0_trend"] == pytest.approx(B0_DAMPING * 0.75) == rec["budget_days"]
    assert rec["budget_ratio"] == pytest.approx(scale)
    assert rec["rule"] == "vec.harness.B0_TREND_H_RULE" and rec["n_trend"] == rec["n_mapped"] == 2


def test_horizon_maps_an_own_factor_lineage_scaled_the_same_way(three):
    """4. A lineage usable at only start and the reference IS in the mapping, with B0's own factor x the same ratio."""
    inputs = _only_start_and_ref(three)
    out, rec = harness.b0_trend_factors(_ctx(inputs), list(BOARD.genes), TYPES, form="horizon")
    scale = (10.5 - 9.5) / (9.5 - 8.75)
    assert set(out) == {"CM", "Endothelial"}
    for lin in out:
        own = B0_DAMPING * (_log_means(inputs["ref"], lin) - _log_means(inputs["start"], lin))
        np.testing.assert_allclose(out[lin], own * scale, atol=1e-12, rtol=0)
    assert rec["n_trend"] == 0 and rec["n_mapped"] == 2
    assert all(r["own_factor"] and r["norm_ratio"] == pytest.approx(scale) for r in rec["lineages"])


@pytest.mark.parametrize("day", [float("nan"), float("inf"), 9.5, 9.0])
def test_horizon_refuses_a_target_day_not_finite_or_not_after_the_reference(three, day):
    view = replace(_view(), target_day=day)
    with pytest.raises(ValueError, match="horizon"):
        harness.b0_trend_factors(_ctx(three, view=view), list(BOARD.genes), TYPES, form="horizon")
    with pytest.raises(ValueError, match="horizon"):
        CANDIDATES["b0_trend_h"].build(_ctx(three, view=view))


def test_eb_and_ts_never_need_the_target_day(three):
    """eb and ts read no target day: with it NaN, under the read guard, they build what they build with it."""
    blind = replace(_view(), target_day=float("nan"))
    for name in ("b0_trend_eb", "b0_trend_ts"):
        with harness.guarded_reads() as g:
            a = CANDIDATES[name].build(_ctx(three, view=blind))
        assert not g["refused"]
        np.testing.assert_array_equal(a.X, CANDIDATES[name].build(_ctx(three)).X)


# ---- each arm
@pytest.mark.parametrize("name", sorted(ARMS))
def test_each_arm_is_registered_like_b0_trend(name):
    """5. Registered with b0_trend's tasks and requirement; applies to exactly the tasks b0_trend applies to."""
    tasks = [t for ts_ in load_tasks().values() for t in ts_]
    c, base = CANDIDATES[name], CANDIDATES["b0_trend"]
    assert [t.key for t in tasks if c.applies(t)] == [t.key for t in tasks if base.applies(t)]
    assert c.tasks == base.tasks and c.derived == () and c.requires is harness._heart_extrap_with_start
    assert c.params["scope"] == "chamber" and c.params["damping"] == B0_DAMPING
    rule = {"eb": "B0_TREND_EB_RULE", "ts": "B0_TREND_TS_RULE", "horizon": "B0_TREND_H_RULE"}[ARMS[name]]
    assert c.params["rule"] == f"vec.harness.{rule}" and isinstance(getattr(harness, rule), str)


@pytest.mark.parametrize("name", sorted(ARMS))
def test_each_arm_refuses_what_b0_trend_refuses(three, name):
    """5. depth, herd and cfm raise with the arm's form; so do non-heart, interpolation and two-stage tasks."""
    form = ARMS[name]
    for kw in ({"depth": True}, {"herd": True}, {"cfm": True}):
        with pytest.raises(ValueError, match="depth, herd or cfm"):
            harness._b0_t2(_ctx(three), trend=True, trend_form=form, **kw)
    c = CANDIDATES[name]
    with pytest.raises(ValueError, match="heart extrapolation"):
        c.build(_ctx(three, view=replace(_view(), kind="interp")))
    with pytest.raises(ValueError, match="heart extrapolation"):
        c.build(_ctx(three, view=replace(_view(), setting="embryo")))
    two = {k: three[k] for k in ("ref", "start")}
    with pytest.raises(ValueError, match="three"):
        c.build(_ctx(two, view=_view(inputs={"ref": "T2/E9.5.h5ad", "start": "T2/E8.75.h5ad"})))


@pytest.mark.parametrize("name", sorted(ARMS))
def test_each_arm_leaves_its_record_and_keeps_b0s_rows_coordinates_and_cell_number(three, name, monkeypatch):
    """5. The arm's record (its rule, every b0_trend field) is left in TREND_ARM_RECORDS; the same seed draws the same
    rows as b0, with b0's coordinates and cell number; held-flat cells equal b0's and a trend lineage moves."""
    seen = []
    real = harness.draw_copy

    def spy(ctx, herd=False):
        out = real(ctx, herd=herd)
        seen.append(out[1].copy())
        return out

    monkeypatch.setattr(harness, "draw_copy", spy)
    harness.TREND_ARM_RECORDS.pop(name, None)
    ctx = _ctx(three)
    a, b = CANDIDATES[name].build(ctx), CANDIDATES["b0"].build(ctx)
    assert len(seen) == 2 and np.array_equal(seen[0], seen[1])
    assert a.n == b.n == BOARD.max_cells
    np.testing.assert_array_equal(a.coords, b.coords)
    merged = np.asarray(three["ref"].obs["merged"]).astype(str)[seen[0]]
    flat = (merged == "Untracked") | (merged == "Tiny")
    assert np.array_equal(a.X[flat], b.X[flat]) and not np.array_equal(a.X[merged == "CM"], b.X[merged == "CM"])
    rec = harness.TREND_ARM_RECORDS[name]
    _, base = harness.b0_trend_factors(ctx, list(BOARD.genes), TYPES)
    assert set(base) <= set(rec) and rec["form"] == ARMS[name] and rec["rule"] != base["rule"]
    assert all(set(r0) <= set(r1) for r0, r1 in zip(base["lineages"], rec["lineages"], strict=True))
    assert not np.array_equal(a.X, CANDIDATES["b0_trend"].build(ctx).X)
