"""Heart arm `b0_trend_quad` (branch campaign/heart-e): b0_trend with each lineage's least-squares line replaced by the
least-squares quadratic over its usable stages, read as q(t_ref + b) - q(t_ref) over b0_trend's own budget b. Only the
per-lineage, per-gene log factor differs from b0_trend. Synthetic data only; nothing scored, no prediction compared
with a target."""

from __future__ import annotations

import hashlib
import json
from dataclasses import replace

import numpy as np
import pytest
from test_b0_trend import (
    BOARD,
    DAYS,
    FOUR_INPUTS,
    TYPES,
    _ctx,
    _log_means,
    _only_start_and_ref,
    _stage,
    _view,
    _with_lineage_cells,
)

from vec import core, harness
from vec.harness import B0_DAMPING, CANDIDATES, load_tasks

NAME = "b0_trend_quad"
BUDGET = B0_DAMPING * (9.5 - 8.75)


@pytest.fixture(scope="module")
def three():
    """test_b0_trend's `three` heart-like stages, built the same way."""
    return {"ref": _stage(1600, seed=2, shift=0.0, tiny=60), "start": _stage(1100, seed=3, shift=0.8, tiny=10),
            "earlier": _stage(900, seed=5, shift=1.6, tiny=10)}


@pytest.fixture(scope="module")
def four(three):
    """test_b0_trend's `four` stages, built the same way."""
    return {**three, "earliest": _stage(1000, seed=7, shift=2.2, tiny=10)}


def _h(a) -> str:
    return hashlib.sha256(np.ascontiguousarray(a, dtype=np.float32).tobytes()).hexdigest()


def _lagrange_change(y: list[np.ndarray], d: list[float], budget: float) -> np.ndarray:
    """The exact interpolating parabola through three points, by hand (Lagrange form, plain loops per gene):
    p(d[2] + budget) - p(d[2])."""
    G = len(y[0])
    out = np.empty(G)
    t1, t0 = d[2] + budget, d[2]
    for g in range(G):
        def p(t, g=g):
            s = 0.0
            for i in range(3):
                term = y[i][g]
                for j in range(3):
                    if j != i:
                        term *= (t - d[j]) / (d[i] - d[j])
                s += term
            return s
        out[g] = p(t1) - p(t0)
    return out


# ------------------------------------------------------------------------------------------------ the fitted change
def test_three_stages_give_the_exact_interpolating_parabolas_change():
    """2. Three stages: the factor is the exact parabola through them, read from the last day over the budget (by hand),
    and q(t_ref) is the reference's own value."""
    rng = np.random.default_rng(31)
    days = [8.5, 8.75, 9.5]
    y = [rng.normal(2.0, 1.0, 60) for _ in days]
    got = harness.trend_quad_log_factor(y, days, BUDGET, 9.5)
    np.testing.assert_allclose(got, _lagrange_change(y, days, BUDGET), rtol=1e-9, atol=1e-12)
    # a parabola with known coefficients: q(t) = a + b t + c t^2, so the change is b B + c ((t+B)^2 - t^2)
    a, b, c = rng.normal(size=60), rng.normal(size=60), rng.normal(size=60)
    yq = [a + b * t + c * t * t for t in days]
    want = b * BUDGET + c * ((9.5 + BUDGET) ** 2 - 9.5 ** 2)
    np.testing.assert_allclose(harness.trend_quad_log_factor(yq, days, BUDGET, 9.5), want, rtol=1e-8, atol=1e-9)


def test_four_stages_give_numpy_polyfits_change():
    """3. Four (and five) stages: the factor is numpy.polyfit's degree-2 fit's change over the budget."""
    rng = np.random.default_rng(37)
    for days in ([8.3, 8.5, 8.75, 9.5], [8.0, 8.3, 8.5, 8.75, 9.5]):
        Y = np.vstack([rng.normal(2.0, 1.0, 70) for _ in days])
        p = np.polyfit(days, Y, 2)
        want = np.polyval(p, 9.5 + BUDGET) - np.polyval(p, 9.5)
        got = harness.trend_quad_log_factor(list(Y), days, BUDGET, 9.5)
        np.testing.assert_allclose(got, want, rtol=1e-8, atol=1e-10)
        # with four or more stages q(t_ref) is a fit, not the reference's value: the arm reads the fit's change
        assert not np.allclose(np.polyval(p, 9.5), Y[-1])


def test_the_quadratic_needs_three_stages():
    with pytest.raises(ValueError, match="three"):
        harness.trend_quad_log_factor([np.ones(3), np.ones(3)], [8.75, 9.5], BUDGET, 9.5)


@pytest.mark.parametrize("days", [[8.5, 8.75, 9.5], [8.3, 8.5, 8.75, 9.5], [8.0, 8.3, 8.5, 8.75, 9.5]])
def test_collinear_means_give_b0_trends_factor_to_rounding(days):
    """4. Stage means exactly on a line in the day: the quadratic's change equals b0_trend's factor to 1e-9 relative
    (the curvature is zero only up to rounding; no exact equality is claimed)."""
    rng = np.random.default_rng(3)
    a, b = rng.uniform(2.0, 4.0, 80), rng.uniform(-1.5, 1.5, 80)
    y = [a + b * (d - 9.5) for d in days]
    quad = harness.trend_quad_log_factor(y, days, BUDGET, days[-1])
    ols = harness.trend_log_factor(y, days, BUDGET)
    np.testing.assert_allclose(quad, ols, rtol=1e-9, atol=0)


def test_collinear_means_through_the_candidate(three, monkeypatch):
    """4 (wiring and built counts). Each stage's pseudobulk set to a linear-in-day design (vec.core.mean_cp10k
    replaced, keyed on the stage's lineage cell number, as test_b0_trend does): the arm's mapping equals b0_trend's to
    1e-9 relative, and the built counts of the two candidates on this case are equal."""
    rng = np.random.default_rng(4)
    G = len(BOARD.genes)
    coef = {lin: (rng.uniform(2.0, 4.0, G), rng.uniform(-1.5, 1.5, G)) for lin in ("CM", "Endothelial")}
    design_of = {}
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
    quad, rec = harness.b0_trend_factors(_ctx(three), list(BOARD.genes), TYPES, form="quad")
    ols, _ = harness.b0_trend_factors(_ctx(three), list(BOARD.genes), TYPES)
    assert set(quad) == set(ols) == {"CM", "Endothelial"}
    for lin in quad:
        np.testing.assert_allclose(quad[lin], ols[lin], rtol=1e-9, atol=0)
    assert all(r["quad"] for r in rec["lineages"])
    a = CANDIDATES[NAME].build(_ctx(three))
    b = CANDIDATES["b0_trend"].build(_ctx(three))
    np.testing.assert_array_equal(a.X, b.X)
    np.testing.assert_array_equal(a.coords, b.coords)


def test_the_arm_reads_the_quadratic_on_the_synthetic_stages(three, four):
    """2, 3 (wiring). Through b0_trend_factors: on three stages each trend lineage's factor is the parabola's change by
    hand from its own log means; on four it is numpy.polyfit's; and it is not b0_trend's (the arm is not idle)."""
    out, rec = harness.b0_trend_factors(_ctx(three), list(BOARD.genes), TYPES, form="quad")
    ols, _ = harness.b0_trend_factors(_ctx(three), list(BOARD.genes), TYPES)
    assert set(out) == set(ols) == {"CM", "Endothelial"}
    for lin in out:
        y = [_log_means(three[r], lin) for r in ("earlier", "start", "ref")]
        np.testing.assert_allclose(out[lin], _lagrange_change(y, [8.5, 8.75, 9.5], BUDGET), rtol=1e-9, atol=1e-12)
        assert not np.allclose(out[lin], ols[lin])
    assert rec["rule"] == "vec.harness.B0_TREND_QUAD_RULE" and rec["n_trend"] == 2
    view = _view(inputs=FOUR_INPUTS)
    out4, _ = harness.b0_trend_factors(_ctx(four, view=view), list(BOARD.genes), TYPES, form="quad")
    for lin in out4:
        d = [8.3, 8.5, 8.75, 9.5]
        Y = np.vstack([_log_means(four[r], lin) for r in ("earliest", "earlier", "start", "ref")])
        p = np.polyfit(d, Y, 2)
        np.testing.assert_allclose(out4[lin], np.polyval(p, 9.5 + BUDGET) - np.polyval(p, 9.5), rtol=1e-8,
                                   atol=1e-10)


# ------------------------------------------------------------------------------------------------ the idle case
def test_a_lineage_with_two_usable_stages_keeps_b0_trends_factor(three):
    """5. Endothelial below min_cells at the earlier stage (usable at start and the reference only): it is left out of
    the mapping as in b0_trend (B0's own factor), with quad False; CM takes the quadratic."""
    inputs = {**three, "earlier": _with_lineage_cells(three["earlier"], "Endothelial", 5)}
    out, rec = harness.b0_trend_factors(_ctx(inputs), list(BOARD.genes), TYPES, form="quad")
    ols, rec_ols = harness.b0_trend_factors(_ctx(inputs), list(BOARD.genes), TYPES)
    assert set(out) == set(ols) == {"CM"}
    rows = {r["lineage"]: r for r in rec["lineages"]}
    assert rows["Endothelial"]["quad"] is False and rows["Endothelial"]["own_factor"]
    assert rows["CM"]["quad"] is True and not rows["CM"]["own_factor"]
    row_ols = next(r for r in rec_ols["lineages"] if r["lineage"] == "Endothelial")
    assert {k: v for k, v in rows["Endothelial"].items() if k != "quad"} == row_ols


def test_a_task_where_no_lineage_has_three_usable_stages_is_b0_trend_bit_for_bit(three):
    """5. No lineage usable at the earlier stage: the arm's output (counts and coordinates) is b0_trend's bit for bit,
    and so is b0's."""
    inputs = _only_start_and_ref(three)
    a = CANDIDATES[NAME].build(_ctx(inputs))
    b = CANDIDATES["b0_trend"].build(_ctx(inputs))
    c = CANDIDATES["b0"].build(_ctx(inputs))
    np.testing.assert_array_equal(a.X, b.X)
    np.testing.assert_array_equal(a.coords, b.coords)
    np.testing.assert_array_equal(a.X, c.X)
    rec = harness.TREND_ARM_RECORDS[NAME]
    assert rec["n_trend"] == 0 and all(r["quad"] is False and r["own_factor"] for r in rec["lineages"])


# ------------------------------------------------------------------------------------------------ shape and record
def test_cells_coordinates_and_number_are_b0_trends(three, monkeypatch):
    """6. The same seed draws the same rows as b0_trend, with its coordinates and cell number; held-flat cells equal
    b0_trend's and a quadratic lineage moves differently."""
    seen = []
    real = harness.draw_copy

    def spy(ctx, herd=False):
        out = real(ctx, herd=herd)
        seen.append(out[1].copy())
        return out

    monkeypatch.setattr(harness, "draw_copy", spy)
    a = CANDIDATES[NAME].build(_ctx(three))
    b = CANDIDATES["b0_trend"].build(_ctx(three))
    assert len(seen) == 2 and np.array_equal(seen[0], seen[1])
    assert a.n == b.n == BOARD.max_cells
    np.testing.assert_array_equal(a.coords, b.coords)
    merged = np.asarray(three["ref"].obs["merged"]).astype(str)[seen[0]]
    flat = (merged == "Untracked") | (merged == "Tiny")
    assert flat.any() and np.array_equal(a.X[flat], b.X[flat])
    assert not np.array_equal(a.X[merged == "CM"], b.X[merged == "CM"])


def test_same_seed_same_output(three, four):
    """7. Two builds at the same seed are equal bit for bit (counts, coordinates, record), on three and four stages."""
    for inputs, view in ((three, _view()), (four, _view(inputs=FOUR_INPUTS))):
        a = CANDIDATES[NAME].build(_ctx(inputs, view=view))
        ra = json.dumps(harness.TREND_ARM_RECORDS[NAME], sort_keys=True)
        b = CANDIDATES[NAME].build(_ctx(inputs, view=view))
        rb = json.dumps(harness.TREND_ARM_RECORDS[NAME], sort_keys=True)
        assert _h(a.X) == _h(b.X) and _h(a.coords) == _h(b.coords) and ra == rb


def test_the_record_gains_only_the_quad_flag(three, four):
    """The arm's record is b0_trend's with the per-lineage boolean "quad" added and the rule named; nothing else is
    added (top-level keys equal; each lineage row gains exactly "quad")."""
    for inputs, view in ((three, _view()), (four, _view(inputs=FOUR_INPUTS))):
        harness.TREND_ARM_RECORDS.pop(NAME, None)
        CANDIDATES[NAME].build(_ctx(inputs, view=view))
        rec = harness.TREND_ARM_RECORDS[NAME]
        base: dict = {}
        harness._b0_t2(_ctx(inputs, view=view), trend=True, trend_record=base)
        assert set(rec) == set(base) and rec["rule"] == "vec.harness.B0_TREND_QUAD_RULE"
        for r0, r1 in zip(base["lineages"], rec["lineages"], strict=True):
            assert set(r1) - set(r0) == {"quad"} and set(r0) <= set(r1)
            assert isinstance(r1["quad"], bool) and r1["quad"] == (not r1["own_factor"])
        for k in ("stages", "days", "budget_days", "eps", "min_cells", "n_trend"):
            assert rec[k] == base[k]


def test_registered_like_b0_trend():
    tasks = [t for ts_ in load_tasks().values() for t in ts_]
    c, base = CANDIDATES[NAME], CANDIDATES["b0_trend"]
    assert [t.key for t in tasks if c.applies(t)] == [t.key for t in tasks if base.applies(t)]
    assert c.tasks == base.tasks == frozenset({"T2"}) and c.derived == ()
    assert c.requires is base.requires is harness._heart_extrap_with_start
    assert c.params["scope"] == "chamber" and c.params["damping"] == B0_DAMPING
    assert c.params["rule"] == "vec.harness.B0_TREND_QUAD_RULE" and isinstance(harness.B0_TREND_QUAD_RULE, str)


def test_refuses_what_b0_trend_refuses_and_never_reads_the_target(three):
    for kw in ({"depth": True}, {"herd": True}, {"cfm": True}):
        with pytest.raises(ValueError, match="depth, herd or cfm"):
            harness._b0_t2(_ctx(three), trend=True, trend_form="quad", **kw)
    c = CANDIDATES[NAME]
    with pytest.raises(ValueError, match="heart extrapolation"):
        c.build(_ctx(three, view=replace(_view(), kind="interp")))
    with pytest.raises(ValueError, match="heart extrapolation"):
        c.build(_ctx(three, view=replace(_view(), setting="embryo")))
    two = {k: three[k] for k in ("ref", "start")}
    with pytest.raises(ValueError, match="three"):
        c.build(_ctx(two, view=_view(inputs={"ref": "T2/E9.5.h5ad", "start": "T2/E8.75.h5ad"})))
    blind = replace(_view(), target_day=float("nan"))
    with harness.guarded_reads() as g:
        a = c.build(_ctx(three, view=blind))
    assert not g["refused"]
    np.testing.assert_array_equal(a.X, c.build(_ctx(three)).X)
