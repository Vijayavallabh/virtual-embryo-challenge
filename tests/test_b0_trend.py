"""Heart candidate `b0_trend` (campaign branch; specification of 2026-11-13): B0 with each lineage's per-gene velocity
the least-squares trend of log(mean CP10k + eps) over every input stage, read over B0's time budget; and the hook it
uses, `vec.t1_scope.momentum_shift(..., logf_by_type=...)`. Synthetic data only; no torch, nothing scored."""

from __future__ import annotations

import hashlib
from dataclasses import replace

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse
from test_h7_gmom import BOARD, _heart_task, _stage

from vec import core, harness, latent_extrap, t1_scope
from vec.harness import B0_DAMPING, CANDIDATES, load_tasks

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
# (campaign/gen, 459916f)
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


# ---------------------------------------------------------------------------------------------------- the trend
FLOW_INPUTS = {"ref": "T2/E9.5.h5ad", "start": "T2/E8.75.h5ad", "earlier": "T2/E8.5.h5ad"}
DAYS = {"earlier": 8.5, "start": 8.75, "ref": 9.5}


def _view(inputs=FLOW_INPUTS, **over):
    return harness.TaskView.of(_heart_task(inputs=inputs, **over))


def _ctx(inputs: dict, view=None, seed: int = harness.BUILD_SEED_OFFSET + 1) -> harness.BuildContext:
    return harness.BuildContext(task=view or _view(), board=BOARD, inputs=inputs, seed=seed)


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
    latent_extrap.RECORDS.clear()
    a = CANDIDATES["b0_trend"].build(_ctx(inputs))
    b = CANDIDATES["b0"].build(_ctx(inputs))
    np.testing.assert_array_equal(a.X, b.X)
    np.testing.assert_array_equal(a.coords, b.coords)
    rec = latent_extrap.RECORDS[-1]
    latent_extrap.RECORDS.clear()
    assert rec["candidate"] == "b0_trend" and rec["n_trend"] == 0
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
    """6. The record's transition weights for gaps (0.2, 0.25, 0.75) and (0.2, 0.25); they sum to 1."""
    w3 = harness.trend_transition_weights(np.cumsum([0.0, 0.2, 0.25, 0.75]))
    w2 = harness.trend_transition_weights(np.cumsum([0.0, 0.2, 0.25]))
    np.testing.assert_allclose(np.round(w3, 4), [0.1119, 0.2192, 0.6689])
    np.testing.assert_allclose(np.round(w2, 4), [0.4262, 0.5738])
    assert w3.sum() == pytest.approx(1.0) and w2.sum() == pytest.approx(1.0)
    assert harness.trend_transition_weights([8.75, 9.5]).tolist() == [1.0]


def test_the_target_is_never_read_and_other_tasks_are_refused(three):
    """7. Built with the target's day NaN and no target among the inputs, under the read guard: nothing refused, and
    the target's day changes nothing. Two input stages, an interpolation task and T1 are refused."""
    blind = replace(_view(), target_day=float("nan"))
    assert set(three) == set(FLOW_INPUTS) and "T2/E10.5.h5ad" not in FLOW_INPUTS.values()
    latent_extrap.RECORDS.clear()
    with harness.guarded_reads() as g:
        a = CANDIDATES["b0_trend"].build(_ctx(three, view=blind))
    assert not g["refused"] and np.isfinite(a.X).all()
    b = CANDIDATES["b0_trend"].build(_ctx(three))
    np.testing.assert_array_equal(a.X, b.X)
    rec = latent_extrap.RECORDS[0]
    latent_extrap.RECORDS.clear()
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
    latent_extrap.RECORDS.clear()
    assert len(seen) == 2 and np.array_equal(seen[0], seen[1])
    assert a.n == b.n == BOARD.max_cells
    np.testing.assert_array_equal(a.coords, b.coords)
    merged = np.asarray(three["ref"].obs["merged"]).astype(str)[seen[0]]
    flat = (merged == "Untracked") | (merged == "Tiny")
    assert flat.any() and np.array_equal(a.X[flat], b.X[flat])
    assert not np.array_equal(a.X[merged == "CM"], b.X[merged == "CM"])


# ---------------------------------------------------------------------------------------------------- b0_trend_sized
# sha256 of b0_trend's output X (float32 bytes) on `three` at BUILD_SEED_OFFSET + 1, computed at 28f25fd (before the
# sized arm was added)
B0_TREND_SHA256 = "8320700817f003ae9b5daff16ac5375fcb8442313f3454e69a048e7bfe87855c"


def test_b0_trend_output_is_unchanged_by_the_sized_arm(three):
    """(e) b0_trend's output on a fixed synthetic case is the one pinned before the sized arm."""
    st = CANDIDATES["b0_trend"].build(_ctx(three))
    latent_extrap.RECORDS.clear()
    assert hashlib.sha256(np.ascontiguousarray(st.X, dtype=np.float32).tobytes()).hexdigest() == B0_TREND_SHA256


def test_the_sized_factor_has_b0s_norm_and_the_trends_direction(three):
    """(a) For every lineage with a trend factor, the applied factor's L2 norm is B0's own and its cosine with the
    natural trend factor is 1 (both to 1e-12); the record keeps the natural direction's figures."""
    ctx = _ctx(three)
    genes = list(BOARD.genes)
    types = ["CM", "Endothelial", "Tiny", "Untracked@end"]
    nat, rec_n = harness.b0_trend_factors(ctx, genes, types)
    siz, rec_s = harness.b0_trend_factors(ctx, genes, types, sized=True)
    assert set(nat) == set(siz) == {"CM", "Endothelial"}
    for lin in nat:
        row_s = next(r for r in rec_s["lineages"] if r["lineage"] == lin)
        row_n = next(r for r in rec_n["lineages"] if r["lineage"] == lin)
        own_norm = np.linalg.norm(nat[lin]) / row_n["norm_ratio"]
        assert np.linalg.norm(siz[lin]) == pytest.approx(own_norm, abs=1e-12, rel=0)
        cos = siz[lin] @ nat[lin] / (np.linalg.norm(siz[lin]) * np.linalg.norm(nat[lin]))
        assert cos == pytest.approx(1.0, abs=1e-12)
        assert row_s["norm_ratio_applied"] == pytest.approx(1.0, abs=1e-12)
        assert row_n["norm_ratio_applied"] == row_n["norm_ratio"] != pytest.approx(1.0)
        assert row_s["cosine"] == row_n["cosine"] and row_s["norm_ratio"] == row_n["norm_ratio"]
    assert rec_s["sized"] is True and rec_n["sized"] is False
    assert rec_s["rule"] == "vec.harness.B0_TREND_SIZED_RULE" and rec_n["rule"] == "vec.harness.B0_TREND_RULE"
    assert rec_s["median_cosine"] == rec_n["median_cosine"]


def test_sized_with_no_lineage_usable_before_start_is_b0_bit_for_bit(three):
    """(b) As test 3, for the sized arm."""
    early = three["earlier"].copy()
    early.obs["merged"] = np.where(np.arange(early.n_obs) % 2 == 0, "Elsewhere", "Untracked")
    inputs = {**three, "earlier": early}
    a = CANDIDATES["b0_trend_sized"].build(_ctx(inputs))
    b = CANDIDATES["b0"].build(_ctx(inputs))
    latent_extrap.RECORDS.clear()
    np.testing.assert_array_equal(a.X, b.X)
    np.testing.assert_array_equal(a.coords, b.coords)


def test_sized_registration_and_refusals(three):
    """(c) As b0_trend: the two heart backtests only, b0's scope, the same refusals."""
    tasks = {t.key: t for ts_ in load_tasks().values() for t in ts_}
    c, n = CANDIDATES["b0_trend_sized"], CANDIDATES["b0_trend"]
    assert all(c.applies(t) == n.applies(t) for t in tasks.values())
    assert {k for k, t in tasks.items() if c.applies(t)} == {"T2h/<=E9.5->E10.5", "T2h/<=E8.75->E10.5"}
    assert c.tasks == n.tasks and c.derived == () and c.requires is n.requires
    assert c.params["scope"] == CANDIDATES["b0"].params["scope"] and c.params["damping"] == B0_DAMPING
    assert c.params["rule"] == "vec.harness.B0_TREND_SIZED_RULE"
    two = {k: three[k] for k in ("ref", "start")}
    with pytest.raises(ValueError, match="three"):
        c.build(_ctx(two, view=_view(inputs={"ref": "T2/E9.5.h5ad", "start": "T2/E8.75.h5ad"})))
    with pytest.raises(ValueError, match="heart extrapolation"):
        c.build(_ctx(three, view=replace(_view(), kind="interp")))
    with pytest.raises(ValueError, match="heart extrapolation"):
        c.build(_ctx(three, view=replace(_view(), task="T1", setting=None)))
    with pytest.raises(ValueError, match="trend must be"):
        harness._b0_t2(_ctx(three), trend="other")


def test_each_arm_leaves_its_own_record(three):
    """(d) Both arms on one task and seed: two records, each with its own candidate name; the sized arm's output
    differs from the natural one's and keeps its coordinates."""
    latent_extrap.RECORDS.clear()
    ctx = _ctx(three)
    with harness.guarded_reads() as g:
        a = CANDIDATES["b0_trend"].build(ctx)
        s = CANDIDATES["b0_trend_sized"].build(ctx)
    assert not g["refused"]
    recs = list(latent_extrap.RECORDS)
    latent_extrap.RECORDS.clear()
    assert [r["candidate"] for r in recs] == ["b0_trend", "b0_trend_sized"]
    assert [r["sized"] for r in recs] == [False, True] and {r["task"] for r in recs} == {ctx.task.key}
    assert not np.array_equal(a.X, s.X) and np.array_equal(a.coords, s.coords)


# ---------------------------------------------------------------------------------------------------- the split arms
# b0_trend_contrast (lineage-specific part) and b0_trend_shared (the unweighted mean over the shifted lineages). The
# synthetic `three` has exactly two shifted lineages (CM, Endothelial; Tiny is short at start), so the mean of two
# equal float64 arrays is that array and the mean of x and -x is 0, both exactly: cases 2 and 3 are exact comparisons.
# sha256 of the outputs X (float32 bytes) on `three` at BUILD_SEED_OFFSET + 1, computed at e9d96aa (before the split
# arms)
PRE_SPLIT_SHA256 = {"b0": "7f9b4113111857bdc351c027f81b80bbde24153a02e07ff99c9cbd31f764e137",
                    "b0_trend": "8320700817f003ae9b5daff16ac5375fcb8442313f3454e69a048e7bfe87855c",
                    "b0_trend_sized": "2eb66dab3ee33c6426a70ac41391c3cee0c44cf89cf53e58c5f8f0546d17c7fb"}
SPLIT_ARMS = ("b0_trend_contrast", "b0_trend_shared")
TYPES = ["CM", "Endothelial", "Tiny", "Untracked@end"]


def _patch_trend(monkeypatch, factors):
    """Replace the trend factor by the given arrays, in call order (one call per lineage with a trend, lineages in
    sorted order: CM then Endothelial), cycling over builds."""
    import itertools

    it = itertools.cycle(factors)
    monkeypatch.setattr(harness, "trend_log_factor", lambda log_means, days, budget: next(it).copy())


def test_outputs_before_the_split_arms_are_unchanged(three):
    """7. b0, b0_trend and b0_trend_sized give the outputs pinned before the split arms."""
    for name, sha in PRE_SPLIT_SHA256.items():
        st = CANDIDATES[name].build(_ctx(three))
        assert hashlib.sha256(np.ascontiguousarray(st.X, dtype=np.float32).tobytes()).hexdigest() == sha, name
    latent_extrap.RECORDS.clear()


def test_contrast_plus_shared_is_the_trend_and_contrasts_sum_to_zero(three):
    """1. contrast + shared = the trend factor per lineage, and the contrasts sum to zero per gene (1e-12); 4 (the
    general case): shared_share_of_total is in (0, 1) and equals n ||shared||^2 / sum_l ||trend_l||^2."""
    ctx, genes = _ctx(three), list(BOARD.genes)
    nat, _ = harness.b0_trend_factors(ctx, genes, TYPES)
    con, rec_c = harness.b0_trend_factors(ctx, genes, TYPES, split="contrast")
    sha, rec_s = harness.b0_trend_factors(ctx, genes, TYPES, split="shared")
    assert set(nat) == set(con) == set(sha) == {"CM", "Endothelial"}
    for lin in nat:
        np.testing.assert_allclose(con[lin] + sha[lin], nat[lin], atol=1e-12, rtol=0)
    np.testing.assert_allclose(sum(con.values()), 0.0, atol=1e-12)
    s = sha["CM"]
    share = 2 * float(s @ s) / sum(float(v @ v) for v in nat.values())
    assert 0 < rec_c["shared_share_of_total"] < 1
    assert rec_c["shared_share_of_total"] == pytest.approx(share, rel=1e-12)
    assert rec_c["n_shifted"] == 2 and rec_c["split"] == "contrast" and rec_s["split"] == "shared"
    assert rec_c["rule"] == "vec.harness.B0_TREND_CONTRAST_RULE" and rec_s["rule"] == "vec.harness.B0_TREND_SHARED_RULE"
    for row in rec_c["lineages"]:
        t, c = nat[row["lineage"]], con[row["lineage"]]
        assert row["norm_trend"] == pytest.approx(np.linalg.norm(t), rel=1e-12)
        assert row["norm_ratio_contrast"] == pytest.approx(np.linalg.norm(c) / np.linalg.norm(t), rel=1e-12)
        assert row["norm_ratio_shared"] == pytest.approx(np.linalg.norm(s) / np.linalg.norm(t), rel=1e-12)
        assert row["cos2_trend_shared"] == pytest.approx((t @ s) ** 2 / ((t @ t) * (s @ s)), rel=1e-12)
        assert "cosine" in row and "norm_ratio" in row  # the natural-trend fields are kept
    assert rec_c["median_norm_ratio_contrast"] == pytest.approx(
        np.median([r["norm_ratio_contrast"] for r in rec_c["lineages"]]))


def test_equal_factors_shared_is_b0_trend_and_contrast_is_the_copy(three, monkeypatch):
    """2 and 4. Every lineage given the same trend factor: b0_trend_shared is b0_trend bit for bit, b0_trend_contrast
    leaves the drawn cells' counts as they are (b0 at damping 0) bit for bit, and shared_share_of_total is 1."""
    x = np.random.default_rng(9).normal(0.0, 0.3, len(BOARD.genes))
    _patch_trend(monkeypatch, [x])
    ctx = _ctx(three)
    latent_extrap.RECORDS.clear()
    a, b = CANDIDATES["b0_trend_shared"].build(ctx), CANDIDATES["b0_trend"].build(ctx)
    np.testing.assert_array_equal(a.X, b.X)
    pool, rows, _, _, shifted = harness._b0_t2_shift(ctx, trend="contrast")
    np.testing.assert_array_equal(shifted.layers["counts"].toarray(), core.counts_of(pool).toarray()[rows])
    recs = {r["candidate"]: r for r in latent_extrap.RECORDS}
    latent_extrap.RECORDS.clear()
    assert recs["b0_trend_shared"]["shared_share_of_total"] == 1.0
    assert recs["b0_trend_contrast"]["shared_share_of_total"] == 1.0


def test_opposite_factors_contrast_is_b0_trend(three, monkeypatch):
    """3 and 4. Lineage factors that sum to zero per gene (x and -x): b0_trend_contrast is b0_trend bit for bit, and
    shared_share_of_total is 0."""
    x = np.random.default_rng(10).normal(0.0, 0.3, len(BOARD.genes))
    _patch_trend(monkeypatch, [x, -x])
    ctx = _ctx(three)
    latent_extrap.RECORDS.clear()
    a, b = CANDIDATES["b0_trend_contrast"].build(ctx), CANDIDATES["b0_trend"].build(ctx)
    np.testing.assert_array_equal(a.X, b.X)
    rec = latent_extrap.RECORDS[0]
    latent_extrap.RECORDS.clear()
    assert rec["candidate"] == "b0_trend_contrast" and rec["shared_share_of_total"] == 0.0


def test_one_shifted_lineage_raises_and_an_own_factor_lineage_is_in_the_split(three):
    """5. With Endothelial gone (relabelled Untracked at every stage) one lineage is shifted and both split arms raise
    (b0_trend still builds). With Endothelial short at the earlier stage it keeps B0's own factor, and it enters the
    mean and both mappings."""
    def relabel(a, lin, to):
        a = a.copy()
        m = np.asarray(a.obs["merged"]).astype(str)
        a.obs["merged"] = np.where(m == lin, to, m)
        return a

    one = {r: relabel(a, "Endothelial", "Untracked") for r, a in three.items()}
    for name in SPLIT_ARMS:
        with pytest.raises(ValueError, match="at least two shifted lineages"):
            CANDIDATES[name].build(_ctx(one))
    CANDIDATES["b0_trend"].build(_ctx(one))
    own_inputs = {**three, "earlier": relabel(three["earlier"], "Endothelial", "Elsewhere")}
    ctx, genes = _ctx(own_inputs), list(BOARD.genes)
    nat, rec = harness.b0_trend_factors(ctx, genes, TYPES)
    assert set(nat) == {"CM"}
    endo = next(r for r in rec["lineages"] if r["lineage"] == "Endothelial")
    assert endo["own_factor"]
    con, _ = harness.b0_trend_factors(ctx, genes, TYPES, split="contrast")
    sha, _ = harness.b0_trend_factors(ctx, genes, TYPES, split="shared")
    assert set(con) == set(sha) == {"CM", "Endothelial"}
    own = con["Endothelial"] + sha["Endothelial"]  # B0's own factor for Endothelial
    np.testing.assert_allclose(sha["CM"], (nat["CM"] + own) / 2, atol=1e-12, rtol=0)
    np.testing.assert_allclose(con["CM"], nat["CM"] - sha["CM"], atol=1e-12, rtol=0)
    latent_extrap.RECORDS.clear()


def test_split_registration_records_and_cells(three, monkeypatch):
    """6 and 8. Registered as b0_trend (tasks, requires, scope, damping), refused where it is; each arm leaves its own
    record; rows, coordinates and cell number are b0's for the same seed."""
    tasks = {t.key: t for ts_ in load_tasks().values() for t in ts_}
    n = CANDIDATES["b0_trend"]
    two = {k: three[k] for k in ("ref", "start")}
    for name in SPLIT_ARMS:
        c = CANDIDATES[name]
        assert all(c.applies(t) == n.applies(t) for t in tasks.values())
        assert c.tasks == n.tasks and c.requires is n.requires and c.derived == ()
        assert c.params["scope"] == n.params["scope"] and c.params["damping"] == B0_DAMPING
        assert c.params["rule"] == {"b0_trend_contrast": "vec.harness.B0_TREND_CONTRAST_RULE",
                                    "b0_trend_shared": "vec.harness.B0_TREND_SHARED_RULE"}[name]
        with pytest.raises(ValueError, match="three"):
            c.build(_ctx(two, view=_view(inputs={"ref": "T2/E9.5.h5ad", "start": "T2/E8.75.h5ad"})))
        with pytest.raises(ValueError, match="heart extrapolation"):
            c.build(_ctx(three, view=replace(_view(), kind="interp")))
        with pytest.raises(ValueError, match="heart extrapolation"):
            c.build(_ctx(three, view=replace(_view(), task="T1", setting=None)))
    seen = []
    real = harness.draw_copy

    def spy(ctx, herd=False):
        out = real(ctx, herd=herd)
        seen.append(out[1].copy())
        return out

    monkeypatch.setattr(harness, "draw_copy", spy)
    ctx = _ctx(three)
    latent_extrap.RECORDS.clear()
    with harness.guarded_reads() as g:
        built = {name: CANDIDATES[name].build(ctx) for name in ("b0", *SPLIT_ARMS)}
    assert not g["refused"]
    assert [r["candidate"] for r in latent_extrap.RECORDS] == list(SPLIT_ARMS)
    latent_extrap.RECORDS.clear()
    assert len(seen) == 3 and all(np.array_equal(s, seen[0]) for s in seen)
    for name in SPLIT_ARMS:
        assert built[name].n == built["b0"].n
        np.testing.assert_array_equal(built[name].coords, built["b0"].coords)
    assert not np.array_equal(built["b0_trend_contrast"].X, built["b0_trend_shared"].X)
