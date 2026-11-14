"""The embryo batch E (campaign branch campaign/emb-e): a within-file mixture of the two read second-upload renders
(vec.t2_build.mix_rows, vec.harness.mix_half, the arm fgw_aligned_mix), the series render with a whole-cell multinomial
redraw (series_rendered_expression(..., kernel="multinomial_redraw"), the arm fgw_aligned_series_mn), a quadratic in
place of the trend's line (trend_factor / trend_rendered_expression(..., form="quadratic"), the arms fgw_quad_render
and fgw_aligned_quad), and the runner experiments/harness/n44_embryo_e_arms.py. Synthetic clouds and stages only; no
data file, committed result or scorer is read.

Tolerances stated here: the multinomial redraw's Monte Carlo mean against its expected counts, 5 standard errors of
that mean per entry plus the rounding of the total (|rint(S) - S| x w / S, which the draw's mean w rint(S) / S carries
by construction); against w rint(S) / S itself, 5 standard errors alone. The quadratic against numpy.polyfit and the
hand interpolation: 1e-9 relative. Collinear profiles: the quadratic's factor against the line's, 1e-9 relative (the
brief's)."""

from __future__ import annotations

import dataclasses
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy import sparse

from vec import harness
from vec import t2_build as tb

REPO = Path(__file__).resolve().parents[1]
HARNESS = REPO / "experiments" / "harness"
BASE = "f29a3eb"  # campaign/emb-d, the branch's base
ARMS = ("fgw_aligned_mix", "fgw_aligned_series_mn", "fgw_quad_render", "fgw_aligned_quad")
MN = tb.SERIES_KERNELS[2]


def _helpers(name: str):
    spec = importlib.util.spec_from_file_location(f"{name}_helpers_e", REPO / "tests" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def polar():
    """tests/test_emb_polar.py's base loader (vec.t2_build at a commit, read with git show)."""
    return _helpers("test_emb_polar")


@pytest.fixture(scope="module")
def emb_expr():
    """The synthetic embryo stages and contexts of tests/test_emb_expr.py (no data file)."""
    return _helpers("test_emb_expr")


@pytest.fixture
def base(tmp_path, polar):
    return polar._base_module(tmp_path, BASE)


def _spy_interp(monkeypatch):
    """Every interpolation the harness builds, captured as returned."""
    seen = []
    real = tb.interpolate_clouds

    def spy(*a, **kw):
        out = real(*a, **kw)
        seen.append(out)
        return out

    monkeypatch.setattr(harness.t2_build, "interpolate_clouds", spy)
    return seen


# ================================================================================ the mixture
# ------------------------------------------------------------------------------------------ M1. the same cells
def test_the_two_read_arms_hold_the_same_cells_positions_and_row_order(emb_expr, monkeypatch):
    c = emb_expr.ctx()
    seen = _spy_interp(monkeypatch)
    s, t = (harness.CANDIDATES[n].build(c) for n in ("fgw_aligned_series", "fgw_aligned_trend"))
    a, b = seen
    for f in ("stage", "rows", "coords", "pool_coords", "pool_weights"):
        assert np.array_equal(getattr(a, f), getattr(b, f)), f
    assert json.dumps(a.record, sort_keys=True, default=str) == json.dumps(b.record, sort_keys=True, default=str)
    assert np.array_equal(s.coords, t.coords) and s.n == t.n == emb_expr.BOARD.max_cells
    assert not np.array_equal(s.X, t.X)  # two different renders of the same cells


# ------------------------------------------------------------------------------------------ M2. the helper
def test_the_mixing_helper_is_s_with_an_all_true_mask_and_t_with_an_all_false_one():
    rng = np.random.default_rng(0)
    S = rng.gamma(1.0, 1.0, size=(9, 5)).astype(np.float32)
    T = rng.gamma(1.0, 1.0, size=(9, 5)).astype(np.float32)
    X, rec = tb.mix_rows(S, T, np.ones(9, dtype=bool))
    assert X.dtype == S.dtype and X.tobytes() == S.tobytes() and rec == {"mix_share_from_s": 1.0}
    X, rec = tb.mix_rows(S, T, np.zeros(9, dtype=bool))
    assert X.tobytes() == T.tobytes() and rec == {"mix_share_from_s": 0.0}
    m = np.array([1, 0, 0, 1, 1, 0, 1, 0, 0], dtype=bool)
    X, rec = tb.mix_rows(S, T, m)
    assert X[m].tobytes() == S[m].tobytes() and X[~m].tobytes() == T[~m].tobytes()
    assert rec == {"mix_share_from_s": 4 / 9}
    with pytest.raises(ValueError, match="boolean"):
        tb.mix_rows(S, T, m.astype(int))
    with pytest.raises(ValueError, match="boolean"):
        tb.mix_rows(S, T, m[:5])
    with pytest.raises(ValueError, match="one shape"):
        tb.mix_rows(S, T[:5], m)


def test_the_half_is_floor_n_over_2_rows_seeded_on_its_own_stream():
    for n in (0, 1, 7, 400, 401):
        m = harness.mix_half(n, 5)
        assert m.dtype == bool and m.shape == (n,) and int(m.sum()) == n // 2
        assert np.array_equal(m, harness.mix_half(n, 5))  # same seed, same half
    assert not np.array_equal(harness.mix_half(400, 5), harness.mix_half(400, 6))  # another seed, another half
    want = np.zeros(400, dtype=bool)
    want[np.random.default_rng([5, harness.SERIES_MIX_STREAM]).permutation(400)[:200]] = True
    assert np.array_equal(harness.mix_half(400, 5), want)
    # not stratified, and not the first rows: over many seeds every row is in the half about half the time
    share = np.mean([harness.mix_half(40, s) for s in range(2000)], axis=0)
    assert np.abs(share - 0.5).max() < 5 * np.sqrt(0.25 / 2000)


# ------------------------------------------------------------------------------------------ M3. the arm
def test_the_mixture_arm_is_series_on_its_half_and_trend_elsewhere(emb_expr, monkeypatch):
    c = emb_expr.ctx()
    seen = _spy_interp(monkeypatch)
    b = {n: harness.CANDIDATES[n].build(c) for n in ("fgw_aligned_mix", "fgw_aligned_series", "fgw_aligned_trend")}
    for r in seen[1:]:
        assert np.array_equal(r.rows, seen[0].rows) and np.array_equal(r.coords, seen[0].coords)
    m = harness.mix_half(b["fgw_aligned_mix"].n, c.seed)
    assert int(m.sum()) == b["fgw_aligned_mix"].n // 2
    X = b["fgw_aligned_mix"].X
    assert X[m].tobytes() == b["fgw_aligned_series"].X[m].tobytes()
    assert X[~m].tobytes() == b["fgw_aligned_trend"].X[~m].tobytes()
    assert b["fgw_aligned_mix"].coords.tobytes() == b["fgw_aligned_series"].coords.tobytes()
    # same seed, same output; another seed, another half
    again = harness.CANDIDATES["fgw_aligned_mix"].build(c)
    assert emb_expr.digest(again) == emb_expr.digest(b["fgw_aligned_mix"])
    assert not np.array_equal(harness.mix_half(400, c.seed), harness.mix_half(400, c.seed + 1))


def test_the_mixture_calls_fgw_with_the_aligned_moves(monkeypatch):
    seen = {}
    monkeypatch.setattr(harness, "_fgw", lambda ctx, **kw: seen.setdefault(ctx, kw))
    for name, kw in (("fgw_aligned_mix", {"move_coupling": "aligned", "expression": "series_trend_mix"}),
                     ("fgw_aligned_series_mn", {"move_coupling": "aligned", "expression": "series_mn"}),
                     ("fgw_quad_render", {"expression": "trend_quad"}),
                     ("fgw_aligned_quad", {"move_coupling": "aligned", "expression": "trend_quad"})):
        c = harness.CANDIDATES[name]
        assert c.tasks == frozenset({"T2"}) and c.name == name and not c.tuned_on and not c.externals
        assert not c.derived and c.build(name) == kw


# ------------------------------------------------------------------------------------------ M4. refusals
def test_the_mixture_applies_and_refuses_where_either_read_arm_does(emb_expr, monkeypatch):
    mix, ser, tr = (harness.CANDIDATES[n] for n in ("fgw_aligned_mix", "fgw_aligned_series", "fgw_aligned_trend"))
    two = emb_expr.task(roles=("before", "after"))
    assert ser.applies(two) and not tr.applies(two) and not mix.applies(two)
    three = emb_expr.task()
    assert ser.applies(three) and tr.applies(three) and mix.applies(three)
    heart = {t.key: t for ts in harness.load_tasks().values() for t in ts}["T2h/<=E9.5->E10.5"]
    assert not mix.applies(heart)
    # the trend refuses (the earlier sample holds no extraembryonic label), the series builds: the mixture refuses
    s = emb_expr.stages()
    s["earlier"].obs["celltype"] = s["earlier"].obs["celltype"].astype(str).str.replace("EXE-", "EMB-")
    c = emb_expr.ctx(inputs=s)
    ser.build(c)
    for cand in (tr, mix):
        with pytest.raises(ValueError, match="no extraembryonic label held by every sample"):
            cand.build(c)
    # both refuse (no label anywhere): the mixture refuses with the series' own refusal
    s = emb_expr.stages()
    for a in s.values():
        a.obs["celltype"] = a.obs["celltype"].astype(str).str.replace("EXE-", "EMB-")
    c = emb_expr.ctx(inputs=s)
    for cand in (ser, mix):
        with pytest.raises(ValueError, match="no matched extraembryonic label"):
            cand.build(c)


# ================================================================================ the multinomial redraw
def _render(emb_expr, monkeypatch, kernel=None, seed=(9, 7), redraw=(9, harness.SERIES_MN_STREAM), stages=None):
    """The series render on the synthetic before/after pair at the aligned positions, with the rendered cells' counts
    and library sizes captured as `_log_normalise` receives them and the expected counts as computed."""
    s = emb_expr.stages() if stages is None else stages
    c = emb_expr.ctx(inputs=s, roles=("before", "after"))
    res, _, _, frac = emb_expr._interp(c, "aligned")
    seen: dict = {}
    real_norm, real_exp, real_mn = tb._log_normalise, tb.series_expected_counts, tb._multinomial_series_redraw

    def norm(C, lib):
        seen["counts"], seen["lib"] = sparse.csr_matrix(C).toarray(), np.asarray(lib, dtype=np.float64)
        return real_norm(C, lib)

    def exp(*a):
        out = real_exp(*a)
        seen["expected"] = out
        return out

    def mn(expected, expected_other, rng):
        seen["rng_state"] = rng.bit_generator.state
        return real_mn(expected, expected_other, rng)

    monkeypatch.setattr(tb, "_log_normalise", norm)
    monkeypatch.setattr(tb, "series_expected_counts", exp)
    monkeypatch.setattr(tb, "_multinomial_series_redraw", mn)
    kw = {} if kernel is None else {"kernel": kernel, "repair_seed": list(redraw)}
    X, rec = tb.series_rendered_expression(s["before"], s["after"], res, list(emb_expr.PANEL), frac, list(seed), **kw)
    monkeypatch.undo()
    return X, rec, seen, res, s


# ------------------------------------------------------------------------------------------ N1. unset and exact_total
def test_kernel_unset_and_exact_total_are_the_bases_bit_for_bit(base, emb_expr):
    s = emb_expr.stages()
    c = emb_expr.ctx(inputs=s, roles=("before", "after"))
    res, _, _, frac = emb_expr._interp(c, "aligned")
    g = list(emb_expr.PANEL)
    X0, r0 = base.series_rendered_expression(s["before"], s["after"], res, g, frac, [9, 7])
    for kw in ({}, {"kernel": "independent"}):
        X1, r1 = tb.series_rendered_expression(s["before"], s["after"], res, g, frac, [9, 7], **kw)
        assert X0.tobytes() == X1.tobytes() and json.dumps(r0, sort_keys=True) == json.dumps(r1, sort_keys=True)
    xt = {"kernel": "exact_total", "repair_seed": [9, harness.SERIES_XT_STREAM]}
    X0, r0 = base.series_rendered_expression(s["before"], s["after"], res, g, frac, [9, 7], **xt)
    X1, r1 = tb.series_rendered_expression(s["before"], s["after"], res, g, frac, [9, 7], **xt)
    assert X0.tobytes() == X1.tobytes() and json.dumps(r0, sort_keys=True) == json.dumps(r1, sort_keys=True)
    # the built candidates: fgw_aligned_series is pinned by sha256 in tests/test_emb_expr.py; the new arm keeps its
    # coupling, moves, draws and cells
    built, xtb, mnb = (harness.CANDIDATES[n].build(c) for n in ("fgw_aligned_series", "fgw_aligned_series_xt",
                                                                 "fgw_aligned_series_mn"))
    assert emb_expr.digest(built) == emb_expr.PINS_AT_E709655[("before", "after")]["fgw_aligned_series"]
    assert built.coords.tobytes() == mnb.coords.tobytes() == xtb.coords.tobytes()
    assert not np.array_equal(built.X, mnb.X) and not np.array_equal(xtb.X, mnb.X)
    with pytest.raises(ValueError, match="repair_seed"):
        tb.series_rendered_expression(s["before"], s["after"], res, g, frac, [9, 7], kernel=MN)
    with pytest.raises(ValueError, match="kernel"):
        tb.series_rendered_expression(s["before"], s["after"], res, g, frac, [9, 7], kernel="multinomial")


# ------------------------------------------------------------------------------------------ N2. the totals
def test_every_rendered_cell_is_at_its_rounded_expected_totals(emb_expr, monkeypatch):
    X, rec, seen, res, _ = _render(emb_expr, monkeypatch, kernel=MN)
    expected, expected_other = seen["expected"]
    C, lib = seen["counts"], seen["lib"]
    assert len(C) == int((res.stage == 0).sum()) > 0  # the shallower bracket (before) is rendered
    np.testing.assert_array_equal(C.sum(1), np.rint(expected.sum(1)))
    np.testing.assert_array_equal(lib - C.sum(1), np.rint(expected_other))
    assert ((C == 0) | (expected > 0)).all()  # a gene whose expected count is 0 stays 0
    assert rec["mn_share_redrawn"] == pytest.approx(float((expected > 0).any(1).mean()))
    assert rec["mn_mean_rounded_expected_total"] == pytest.approx(float(np.rint(expected.sum(1)).mean()))
    _, rec0, _, _, _ = _render(emb_expr, monkeypatch)
    assert set(rec) - set(rec0) == {"mn_share_redrawn", "mn_mean_rounded_expected_total"}
    same = set(rec0) - {"rendered_median"}  # the rendered cells' median library follows the redraw
    assert {k: rec[k] for k in same} == {k: rec0[k] for k in same}
    assert rec["rendered_median"] == float(np.median(lib))


def test_the_redraw_by_hand():
    expected = np.array([[2.0, 1.0, 0.0, 1.2],  # total 4.2 -> 4; gene 2 cannot be drawn
                         [0.0, 0.0, 0.0, 0.0],  # nothing to draw: the cell draws nothing
                         [0.1, 0.2, 0.0, 0.1]])  # total 0.4 -> 0
    rng = np.random.default_rng(4)
    C, o, rec = tb._multinomial_series_redraw(expected, np.array([10.4, 3.0, 2.6]), rng)
    assert C.dtype == np.int64 and C.sum(1).tolist() == [4, 0, 0] and (C[:, 2] == 0).all()
    assert o.tolist() == [10.0, 3.0, 3.0]
    assert rec == {"mn_share_redrawn": pytest.approx(2 / 3), "mn_mean_rounded_expected_total": pytest.approx(4 / 3)}
    # one stream, cells in row order: the draws are those of the same generator called cell by cell
    g = np.random.default_rng(4)
    first = g.multinomial(4, expected[0, [0, 1, 3]] / 4.2)
    assert C[0, [0, 1, 3]].tolist() == first.tolist()


# ------------------------------------------------------------------------------------------ N3. Monte Carlo
def test_the_redraws_mean_over_many_seeds_matches_the_expected_counts(emb_expr, monkeypatch):
    _, _, seen, _, _ = _render(emb_expr, monkeypatch, kernel=MN)
    expected = seen["expected"][0][:6]
    reps = 4000
    draws = np.stack([tb._multinomial_series_redraw(expected, np.zeros(len(expected)),
                                                    np.random.default_rng([s, harness.SERIES_MN_STREAM]))[0]
                      for s in range(reps)]).astype(np.float64)
    m, se = draws.mean(0), draws.std(0) / np.sqrt(reps)
    S = expected.sum(1, keepdims=True)
    scaled = expected * np.rint(S) / S  # the multinomial's own mean
    assert (np.abs(m - scaled) <= 5.0 * se + 1e-12).all()
    assert (np.abs(m - expected) <= 5.0 * se + np.abs(np.rint(S) - S) * expected / S + 1e-12).all()


# ------------------------------------------------------------------------------------------ N4. the deeper bracket
def test_the_deeper_brackets_cells_are_bit_for_bit_as_before(emb_expr, monkeypatch):
    X0, _, _, res, _ = _render(emb_expr, monkeypatch)
    X1, _, _, _, _ = _render(emb_expr, monkeypatch, kernel=MN)
    hi = res.stage == 1
    assert hi.any() and X0[hi].tobytes() == X1[hi].tobytes()
    assert not np.array_equal(X0[~hi], X1[~hi])


# ------------------------------------------------------------------------------------------ N5. the stream
def test_the_redraw_stream_is_its_own(emb_expr, monkeypatch):
    from vec import t2_latent

    before = {*range(10), harness.NBHD_STREAM, harness.SERIES_XT_STREAM, t2_latent.EXPRESSION_STREAM,
              *tb.HERD_STREAMS, *tb.MIDPOINT_STREAMS}  # every stream a build used before this batch
    assert harness.SERIES_MIX_STREAM not in before
    used = before | {harness.SERIES_MIX_STREAM}
    assert harness.SERIES_MN_STREAM not in used
    seed = harness.BUILD_SEED_OFFSET + 4
    for new, others in ((harness.SERIES_MN_STREAM, used), (harness.SERIES_MIX_STREAM, before)):
        mine = np.random.default_rng([seed, new]).random(64)
        for k in others:
            assert not np.array_equal(mine, np.random.default_rng([seed, k]).random(64)), (new, k)
    # the render's own generator is [seed, 7]; the redraw receives a fresh [seed, SERIES_MN_STREAM]
    _, _, seen, _, _ = _render(emb_expr, monkeypatch, kernel=MN, seed=(seed, 7),
                               redraw=(seed, harness.SERIES_MN_STREAM))
    assert seen["rng_state"] == np.random.default_rng([seed, harness.SERIES_MN_STREAM]).bit_generator.state
    # the redraw's draws differ from the same redraw on the repair's stream
    a, _, _, _, _ = _render(emb_expr, monkeypatch, kernel=MN, redraw=(9, harness.SERIES_MN_STREAM))
    b, _, _, _, _ = _render(emb_expr, monkeypatch, kernel=MN, redraw=(9, harness.SERIES_XT_STREAM))
    assert not np.array_equal(a, b)
    # the harness passes the build seed and the redraw stream
    calls = []
    real = tb.series_rendered_expression

    def spy(*a, **kw):
        calls.append(kw)
        return real(*a, **kw)

    monkeypatch.setattr(harness.t2_build, "series_rendered_expression", spy)
    c = emb_expr.ctx(roles=("before", "after"))
    harness.CANDIDATES["fgw_aligned_series_mn"].build(c)
    assert calls == [{"kernel": MN, "repair_seed": [c.seed, harness.SERIES_MN_STREAM]}]


# ------------------------------------------------------------------------------------------ N6. reproducible
def test_the_redraw_is_reproducible(emb_expr, monkeypatch):
    a, ra, _, _, _ = _render(emb_expr, monkeypatch, kernel=MN)
    b, rb, _, _, _ = _render(emb_expr, monkeypatch, kernel=MN)
    assert a.tobytes() == b.tobytes() and ra == rb
    ctx = emb_expr.ctx(roles=("before", "after"))
    x, y = (harness.CANDIDATES["fgw_aligned_series_mn"].build(ctx) for _ in range(2))
    assert emb_expr.digest(x) == emb_expr.digest(y)


def test_the_mn_arm_is_registered_with_the_second_uploads_requirements():
    c, base_arm = harness.CANDIDATES["fgw_aligned_series_mn"], harness.CANDIDATES["fgw_aligned_series"]
    assert c.requires is base_arm.requires and c.params["kernel"] == MN and c.params["rule"] == tb.SERIES_RULE


# ================================================================================ the quadratic trend
def _prof(M):
    M = M.toarray().astype(float)
    return (M / np.maximum(M.sum(1, keepdims=True), 1)).mean(0)


def _samples(n, seed=3):
    rng = np.random.default_rng(seed)
    C = [sparse.csr_matrix(rng.poisson(rng.gamma(2.0, 2.5, 6), size=(200, 6)) + 1) for _ in range(n)]
    return C, [np.array(["EXE-A"] * 200)] * n


def _finish(fit, den):
    want = np.exp(fit) / (den + tb.SERIES_EPS)
    return np.clip(want / (den @ want / den.sum()), *tb.SERIES_CLIP)


# ------------------------------------------------------------------------------------------ Q1. unset is the line
def test_form_unset_is_the_line_bit_for_bit(base, emb_expr):
    C, names = _samples(3)
    days = [6.75, 7.25, 8.0]
    sam = list(zip(C, names, strict=True))
    old = base.trend_factor(C[1], names[1], sam, days, 7.5, ["EXE-A"])
    for kw in ({}, {"form": "line"}):
        assert tb.trend_factor(C[1], names[1], sam, days, 7.5, ["EXE-A"], **kw).tobytes() == old.tobytes()
    c = emb_expr.ctx()
    g = list(emb_expr.BOARD.genes)
    samples = [(a, c.task.role_days[r]) for r, a in c.inputs.items()]
    for mc in tb.MOVE_COUPLINGS:
        res = emb_expr._interp(c, mc)[0]
        X0, r0 = base.trend_rendered_expression(c.inputs["before"], c.inputs["after"], res, g, samples,
                                                c.task.target_day, [c.seed, 7])
        for kw in ({}, {"form": "line"}):
            X1, r1 = tb.trend_rendered_expression(c.inputs["before"], c.inputs["after"], res, g, samples,
                                                  c.task.target_day, [c.seed, 7], **kw)
            assert X0.tobytes() == X1.tobytes() and json.dumps(r0, sort_keys=True) == json.dumps(r1, sort_keys=True)
    # the built trend arms: fgw_trend_render pinned at e709655; fgw_aligned_trend is the base's render at the aligned
    # positions
    built = {n: harness.CANDIDATES[n].build(c) for n in ("fgw_trend_render", "fgw_aligned_trend", "fgw_aligned")}
    assert emb_expr.digest(built["fgw_trend_render"]) == emb_expr.PINS_AT_E709655[
        ("earlier", "before", "after")]["fgw_trend_render"]
    assert built["fgw_aligned_trend"].X.tobytes() == built["fgw_trend_render"].X.tobytes()
    assert built["fgw_aligned_trend"].coords.tobytes() == built["fgw_aligned"].coords.tobytes()
    with pytest.raises(ValueError, match="form"):
        tb.trend_factor(C[1], names[1], sam, days, 7.5, ["EXE-A"], form="cubic")


# ------------------------------------------------------------------------------------------ Q2. three samples
def test_three_samples_give_the_exact_interpolating_parabola():
    C, names = _samples(3)
    days = np.array([6.75, 7.25, 8.0])
    x = 7.5
    phi = tb.trend_factor(C[1], names[1], list(zip(C, names, strict=True)), days, x, ["EXE-A"], form="quadratic")
    Y = np.log(np.vstack([_prof(M) for M in C]) + tb.SERIES_EPS)
    lag = sum(Y[i] * np.prod([(x - days[j]) / (days[i] - days[j]) for j in range(3) if j != i]) for i in range(3))
    np.testing.assert_allclose(phi, _finish(lag, _prof(C[1])), rtol=1e-9)
    # the parabola passes through the samples: read at a sample's own day it gives that sample's profile
    for i in range(3):
        at = tb.trend_factor(C[i], names[i], list(zip(C, names, strict=True)), days, days[i], ["EXE-A"],
                             form="quadratic")
        np.testing.assert_allclose(at, np.clip(np.ones(6), *tb.SERIES_CLIP), rtol=1e-9)
    line = tb.trend_factor(C[1], names[1], list(zip(C, names, strict=True)), days, x, ["EXE-A"])
    assert np.abs(np.log(phi / line)).max() > 1e-3  # not the line on these samples


# ------------------------------------------------------------------------------------------ Q3. four samples
def test_four_samples_give_numpys_polyfit():
    C, names = _samples(4, seed=5)
    days = [6.75, 7.25, 7.5, 8.0]
    for x in (7.75, 8.5):
        phi = tb.trend_factor(C[2], names[2], list(zip(C, names, strict=True)), days, x, ["EXE-A"],
                              form="quadratic")
        Y = np.log(np.vstack([_prof(M) for M in C]) + tb.SERIES_EPS)
        fit = np.polyfit(days, Y, 2)
        np.testing.assert_allclose(phi, _finish(fit[0] * x * x + fit[1] * x + fit[2], _prof(C[2])), rtol=1e-9)


# ------------------------------------------------------------------------------------------ Q4. collinear profiles
def _collinear_profiles(monkeypatch, day_of_rows: dict[int, float]):
    """Every profile the trend reads replaced by exp(a + b day) - eps, its sample's day found by its number of rows:
    log(profile + eps) is exactly collinear in day up to the rounding of exp and log."""
    rng = np.random.default_rng(12)
    a, b = np.log(rng.gamma(2.0, 0.02, 24)), rng.normal(0.0, 0.4, 24)

    def prof(C, names, labels, w):
        return np.exp(a + b * day_of_rows[C.shape[0]]) - tb.SERIES_EPS

    monkeypatch.setattr(tb, "_label_profile", prof)


def test_collinear_profiles_give_the_lines_factor(emb_expr, monkeypatch):
    s = emb_expr.stages()
    c = emb_expr.ctx(inputs=s)
    rows = {a.n_obs: c.task.role_days[r] for r, a in c.inputs.items()}
    assert len(rows) == 3  # the three samples differ in size
    _collinear_profiles(monkeypatch, rows)
    names = [np.asarray(a.obs["celltype"]).astype(str) for a in c.inputs.values()]
    counts = [sparse.csr_matrix(a.layers["counts"]) for a in c.inputs.values()]
    days = list(rows.values())
    for i in range(3):
        sam = list(zip(counts, names, strict=True))
        line = tb.trend_factor(counts[i], names[i], sam, days, c.task.target_day, ["EXE-Endoderm"])
        quad = tb.trend_factor(counts[i], names[i], sam, days, c.task.target_day, ["EXE-Endoderm"], form="quadratic")
        np.testing.assert_allclose(quad, line, rtol=1e-9)
    # the built files on this case (reported, not required: no exact idle case is claimed)
    q, t = (harness.CANDIDATES[n].build(c) for n in ("fgw_quad_render", "fgw_trend_render"))
    assert q.coords.tobytes() == t.coords.tobytes()
    print(f"collinear case: built expression equal = {q.X.tobytes() == t.X.tobytes()}")


# ------------------------------------------------------------------------------------------ Q5. three days or refuse
def test_fewer_than_three_distinct_days_refuse(emb_expr, monkeypatch):
    C, names = _samples(3)
    with pytest.raises(ValueError, match="three distinct days"):
        tb.trend_factor(C[1], names[1], list(zip(C, names, strict=True)), [6.75, 7.25, 7.25], 7.5, ["EXE-A"],
                        form="quadratic")
    tb.trend_factor(C[1], names[1], list(zip(C, names, strict=True)), [6.75, 7.25, 7.25], 7.5, ["EXE-A"])  # the line
    # the builder's own refusal: the trend sample (a role the size rule does not read) given the before bracket's day
    c = emb_expr.trend_ctx()
    (extra,) = [r for r in c.inputs if r.startswith(harness.TREND_SAMPLE_PREFIX)]
    view = dataclasses.replace(c.task, role_days={**c.task.role_days, extra: c.task.role_days["before"]})
    c2 = dataclasses.replace(c, task=view)
    harness.CANDIDATES["fgw_trend_render"].build(c2)
    for name in ("fgw_quad_render", "fgw_aligned_quad"):
        with pytest.raises(ValueError, match="three distinct days"):
            harness.CANDIDATES[name].build(c2)
    # registration: the trend arms' requirements plus the three-day condition
    for name, tr in (("fgw_quad_render", "fgw_trend_render"), ("fgw_aligned_quad", "fgw_aligned_trend")):
        cand = harness.CANDIDATES[name]
        assert cand.applies(emb_expr.task()) and not cand.applies(emb_expr.task(roles=("before", "after")))
        assert cand.params["form"] == "quadratic" and cand.params["rule"] == tb.TREND_RENDER_RULE
        monkeypatch.setitem(harness.STAGE_DAYS, emb_expr.FILES["earlier"], harness.STAGE_DAYS[emb_expr.FILES["before"]])
        assert harness.CANDIDATES[tr].applies(emb_expr.task()) and not cand.applies(emb_expr.task())
        monkeypatch.undo()


# ------------------------------------------------------------------------------------------ Q6. the arms
@pytest.mark.parametrize("arm, base_arm", [("fgw_quad_render", "fgw_trend_render"),
                                           ("fgw_aligned_quad", "fgw_aligned_trend")])
def test_the_quadratic_arms_keep_their_bases_cells_and_are_reproducible(emb_expr, arm, base_arm):
    c = emb_expr.ctx()
    a, again, b = (harness.CANDIDATES[n].build(c) for n in (arm, arm, base_arm))
    assert emb_expr.digest(a) == emb_expr.digest(again)  # same seed, same output
    assert a.coords.tobytes() == b.coords.tobytes() and a.n == b.n
    assert not np.array_equal(a.X, b.X)  # the quadratic is not the line on these samples
    g = list(emb_expr.BOARD.genes)
    res = emb_expr._interp(c, "aligned" if "aligned" in arm else "fgw")[0]
    _, rec = tb.trend_rendered_expression(c.inputs["before"], c.inputs["after"], res, g,
                                          [(x, c.task.role_days[r]) for r, x in c.inputs.items()],
                                          c.task.target_day, [c.seed, 7], form="quadratic")
    assert rec["form"] == "quadratic"


def test_the_arms_apply_on_the_counted_backtests_and_the_nine_tasks(n44):
    tasks = {t.key: t for ts in harness.load_tasks().values() for t in ts}
    for key in ("T2e/E7.25+E8.0->E7.5", "T2e/E6.75+E7.5->E7.25"):
        assert all(harness.CANDIDATES[n].applies(tasks[key]) for n in ARMS)
    assert not any(harness.CANDIDATES[n].applies(tasks["T2h/<=E9.5->E10.5"]) for n in ARMS)
    nine = n44.all_tasks()
    assert set(nine) == {"half", "late", "new1", "new2", "e75", "e725", "x1", "x2", "x3"}
    for name in n44.SCORED:
        assert all(harness.CANDIDATES[name].applies(t) for t in nine.values()), name


# ================================================================================ the runner
@pytest.fixture(scope="module")
def n44():
    sys.path.insert(0, str(HARNESS))
    spec = importlib.util.spec_from_file_location("n44_embryo_e_arms", HARNESS / "n44_embryo_e_arms.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_the_runner_scores_the_floor_the_four_arms_and_the_four_rescored_arms(n44):
    assert n44.NEW == ARMS
    assert n44.RESCORED == ("fgw", "fgw_aligned_series", "fgw_trend_render", "fgw_aligned_trend")
    assert n44.ARMS == (harness.FLOOR, *ARMS, *n44.RESCORED)


def test_the_runners_task_file_loads_through_nine_tasks(tmp_path, monkeypatch, n44):
    calls = []

    def fake_run_task(t, cands, seeds, data_dir):  # synthetic scores: no data, no scorer
        calls.append((t, [c.name for c in cands], list(seeds)))
        rows = [{"seed": s, "scores": {c.name: 50.0 + 0.1 * s + i for i, c in enumerate(cands)}} for s in seeds]
        return {"task": {"key": t.key, "inputs": dict(t.inputs)}, "rows": rows,
                "provenance": {"commit": "c" * 40, "tree_clean": True}}

    monkeypatch.setattr(n44.harness, "run_task", fake_run_task)
    out = tmp_path / "runs"
    n44.main(["run", "--task", "x1", "--out", str(out), "--data", str(tmp_path / "no_data")])
    t = n44.all_tasks()["x1"]
    assert any(r.startswith(harness.TREND_SAMPLE_PREFIX) for r in t.inputs)  # the trend samples are handed over
    assert calls == [(t, list(n44.ARMS), list(n44.SEEDS))]  # one call with the floor and the eight arms
    (path,) = out.glob("T2e_*.json")
    data = json.loads(path.read_text())
    assert data["task_source"] == "experiments/harness/n44_embryo_e_arms.py" and data["not_buildable"] == {}
    assert data["sample_rule"]["samples"] == n44.n36.trend_sample_files(t.inputs["before"], t.inputs["after"],
                                                                         t.target)
    scores, why, meta = n44.nt._load([out])
    assert why == [] and meta == ["runs:ccccccc"]
    assert set(scores[t.key]) == set(n44.ARMS) and all(len(v) == 10 for v in scores[t.key].values())
    for name in n44.RESCORED:  # the rescored arms under their own names
        assert name in scores[t.key]
    with pytest.raises(SystemExit, match="exists"):
        n44.main(["run", "--task", "x1", "--out", str(out), "--data", str(tmp_path / "no_data")])
    assert len(calls) == 1
    with pytest.raises(SystemExit, match="unknown task"):
        n44.main(["run", "--task", "nope", "--out", str(out)])


def _fake_runner(score: dict[str, float]):
    """run_task's shape without data or scorer: each candidate is built (its build may raise) and given a constant."""
    calls = []

    def runner(t, cands, seeds, data_dir):
        calls.append([c.name for c in cands])
        rows = []
        for s in seeds:
            row = {"seed": s, "scores": {"floor": 50.0}, "skills": {"floor": {}}, "metrics": {"floor": {}},
                   "n_cells": {"floor": 1}, "groups": {"floor": {}}, "reads": {}}
            for c in cands:
                if c.name == "floor":
                    continue
                c.build(dataclasses.make_dataclass("Ctx", ["seed"])(seed=harness.BUILD_SEED_OFFSET + s))
                row["scores"][c.name] = score[c.name]
                for k in ("skills", "metrics", "groups"):
                    row[k][c.name] = {}
                row["n_cells"][c.name], row["reads"][c.name] = 1, []
            rows.append(row)
        return {"task": {"key": t.key, "inputs": dict(t.inputs)}, "rows": rows,
                "candidates": [c.name for c in cands], "params": {c.name: {} for c in cands},
                "provenance": {"commit": "a" * 40, "tree_clean": True}}

    return runner, calls


def _patch_builds(monkeypatch, n44, fail_on: set[str], key: str, seed_fails: int = 0):
    """Every scored arm's build replaced: the arms in fail_on raise the quadratic's refusal from seed_fails."""
    cands = harness.CANDIDATES

    def ok(ctx):
        return "built"

    def bad(ctx):
        if ctx.seed - harness.BUILD_SEED_OFFSET >= seed_fails:
            raise ValueError(f"{key}: a quadratic trend needs at least 3 distinct days")
        return "built"

    for x in n44.SCORED:
        monkeypatch.setitem(cands, x, dataclasses.replace(cands[x], build=bad if x in fail_on else ok))


def test_an_arm_the_builder_refuses_is_recorded_and_the_other_arms_are_still_scored(tmp_path, monkeypatch, n44):
    t = n44.all_tasks()["e725"]
    score = {n: 50.0 + i for i, n in enumerate(n44.SCORED)}
    _patch_builds(monkeypatch, n44, {"fgw_quad_render"}, t.key, seed_fails=2)
    runner, calls = _fake_runner(score)
    res = n44.score_task(t, tmp_path, runner=runner)
    assert calls == [["floor", *n44.SCORED], *[["floor", n] for n in n44.SCORED]]  # together first, then arm by arm
    assert set(res["not_buildable"]) == {"fgw_quad_render"}
    rec = res["not_buildable"]["fgw_quad_render"]
    assert rec["status"] == n44.n36.NOT_BUILDABLE and rec["seed"] == 2
    assert rec["message"] == "ValueError: <task>: a quadratic trend needs at least # distinct days"
    scored = {"floor", *n44.SCORED} - {"fgw_quad_render"}
    assert all(set(r["scores"]) == scored for r in res["rows"]) and set(res["candidates"]) == scored
    assert [r["seed"] for r in res["rows"]] == list(n44.SEEDS)
    monkeypatch.setattr(n44.harness, "run_task", runner)
    out = tmp_path / "runs"
    n44.main(["run", "--task", "e725", "--out", str(out), "--data", str(tmp_path / "no_data")])
    scores, why, _ = n44.nt._load([out])
    assert why == [] and set(scores[t.key]) == scored
    (path,) = out.glob("T2e_*.json")
    assert set(json.loads(path.read_text())["not_buildable"]) == {"fgw_quad_render"}
    # every arm failing still leaves the floor scored
    _patch_builds(monkeypatch, n44, set(n44.SCORED), t.key)
    runner, calls = _fake_runner(score)
    res = n44.score_task(t, tmp_path, runner=runner)
    assert set(res["not_buildable"]) == set(n44.SCORED) and calls[-1] == ["floor"]
    assert all(set(r["scores"]) == {"floor"} for r in res["rows"])
    # nothing fails: one run, nothing recorded
    _patch_builds(monkeypatch, n44, set(), t.key)
    runner, calls = _fake_runner(score)
    res = n44.score_task(t, tmp_path, runner=runner)
    assert calls == [["floor", *n44.SCORED]] and res["not_buildable"] == {}


def test_a_failure_outside_a_build_stops_the_task(tmp_path, n44):
    def runner(t, cands, seeds, data_dir):
        raise OSError("the data are not there")

    with pytest.raises(OSError, match="not there"):
        n44.score_task(n44.all_tasks()["x2"], tmp_path, runner=runner)
