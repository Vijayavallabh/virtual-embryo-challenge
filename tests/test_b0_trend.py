"""Heart candidate `b0_trend` (row H11 of experiments/harness/TRIES.md; both heart files by the captain's ruling of
2026-11-13, docs/PLAN.md decision log (w); ported from campaign/trend 8bb23bf): B0 with each lineage's per-gene velocity
the least-squares trend of log(mean CP10k + eps) over every input stage, read over B0's time budget; and the hook it
uses, `vec.t1_scope.momentum_shift(..., logf_by_type=...)`. Synthetic data only; nothing scored."""

from __future__ import annotations

import hashlib
from dataclasses import replace

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
