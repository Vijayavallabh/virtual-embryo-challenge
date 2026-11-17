"""experiments/t1_m1/gmom_size_read.py (try T49 of experiments/harness/TRIES.md) on synthetic data: the fold's fit and
its three finishes, the arms of both forms and the reading. No released file, statistics file or score is read."""

from __future__ import annotations

import importlib.util
import json
import sys

import numpy as np
import pytest
from test_g_mom_verify import _small_inputs
from test_t1_state_transport import _gmom_ctx

from vec import harness, t1_scope
from vec import t1_state_transport as ts
from vec.boards import REPO_ROOT

SCRIPT = REPO_ROOT / "experiments" / "t1_m1" / "gmom_size_read.py"
GROUPS = {"V-CM": "CM", "Endocardium": "Endocardium", "Blood": "Blood", "Foregut": "other:Foregut"}


@pytest.fixture
def sm(monkeypatch, tmp_path):
    spec = importlib.util.spec_from_file_location("gmom_size_read_module", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    monkeypatch.setattr(mod, "RESULTS", tmp_path / "results")
    mod._FOLD.clear()
    harness._G_MOM_PREP.clear()
    return mod


def _world(seed=7):
    ref, start, _ = _small_inputs(seed)
    ctx = _gmom_ctx(ref, start)
    groups = {lab: GROUPS.get(lab, "x") for lab in ref.obs[t1_scope.LABEL_KEY].astype(str).unique()}
    rng = np.random.default_rng(5)
    dirs = {g: rng.normal(0, 0.4, ref.n_vars) for g in ("CM", "Endocardium")}
    return ref, start, ctx, groups, dirs, list(ctx.derived["t1_held_flat"])


def _classified(ctx):
    ref, start = harness._shallow(ctx.inputs["ref"]), harness._shallow(ctx.inputs["start"])
    src, start_src = ctx.task.role_stems["ref"], ctx.task.role_stems["start"]
    t1_scope.classify(ref, source=src)
    t1_scope.classify(start, source=start_src)
    return ref, src, start, start_src, harness.t1_chamber_shares(ref, src)


# ---------------------------------------------------------------------------------------------------------------------
# The fit and its finishes


def test_the_dampings_are_one_one_and_a_half_and_two_times_b0s(sm):
    assert sm.DAMPINGS == {"g_mom": 0.5, "g_mom_d075": 0.75, "g_mom_d1": 1.0, "m1_r_gmom": 0.5,
                           "m1_r_gmom_g075": 0.75, "m1_r_gmom_g1": 1.0}
    assert sm.LEVELS == (0.5, 0.75, 1.0) and sm.BASE_DAMPING == harness.B0_DAMPING == 0.5
    assert sm.REPRO_TOL == 1e-9


def test_the_finish_at_half_is_g_mom_factors_own_bit_for_bit_and_the_fit_is_made_once(sm, monkeypatch):
    fits = []
    real_fit = ts.fit_state_factors
    monkeypatch.setattr(ts, "fit_state_factors", lambda *a, **k: fits.append(1) or real_fit(*a, **k))
    _, _, ctx, groups, dirs, held = _world()
    fold = sm.fold_factors(ctx)
    assert len(fits) == 1
    assert fold["finished"][0.5] is harness._G_MOM_PREP["fold"]["factors"]
    # a fresh call of _g_mom_factors itself, its cache cleared, gives the same factors bit for bit
    harness._G_MOM_PREP.clear()
    direct = harness._g_mom_factors(ctx, *_classified(ctx))
    assert len(fits) == 2 and direct is not fold["finished"][0.5]
    assert sm.same_factors(direct, fold["finished"][0.5])
    harness._G_MOM_PREP.clear()
    sm._FOLD.clear()
    fits.clear()
    # one run: two seeds, every arm of both forms, one fit
    for seed in (1000000, 1000001):
        c = _gmom_ctx(ctx.inputs["ref"], ctx.inputs["start"], seed=seed)
        sm.form_a_arms(c)
        sm.form_b_arms(c, dirs, groups, 0.2, held, held)
    assert len(fits) == 1


def test_the_larger_finishes_are_finish_momentum_at_their_dampings_not_a_rescaling(sm):
    _, _, ctx, *_ = _world()
    fold = sm.fold_factors(ctx)
    raw = fold["raw"]
    for d in (0.75, 1.0):
        assert sm.same_factors(fold["finished"][d], ts.finish_momentum(raw, d))
        assert f"damped by {d:g}" in fold["finished"][d].record["finish"]
    revived = False
    for c, cf in fold["finished"][0.5].categories.items():
        big = fold["finished"][1.0].categories[c].logf
        revived |= bool(np.any((cf.logf == 0) & (big != 0)))
        doubled = ts._zero_and_clip(cf.logf * 2.0)
        assert not np.array_equal(doubled, big)  # not the 0.5 factors scaled up: entries zeroed at 0.5 come back
        assert np.abs(big).max() <= ts.CLIP_AT
    assert revived
    rec = sm.factor_records(fold["finished"], ctx.derived["t1_held_flat"])
    assert set(rec) == {"0.5", "0.75", "1"}
    for c in rec["0.5"]:
        shares = [rec[d][c]["share_not_zeroed"] for d in ("0.5", "0.75", "1")]
        assert shares == sorted(shares)
        assert rec["0.5"][c]["mean_abs_log_not_zeroed"] < rec["1"][c]["mean_abs_log_not_zeroed"]
        assert set(rec["1"][c]) == {"share_not_zeroed", "share_at_clip", "mean_abs_log_not_zeroed",
                                    "mean_abs_log_all_entries", "share_unzeroed_vs_0.5"}
        assert rec["0.5"][c]["share_unzeroed_vs_0.5"] == 0.0
        assert rec["0.5"][c]["mean_abs_log_all_entries"] < rec["1"][c]["mean_abs_log_all_entries"]


def test_a_build_that_reaches_another_fit_is_refused(sm):
    _, _, ctx, *_ = _world()
    fold = sm.fold_factors(ctx)
    with pytest.raises(RuntimeError, match="other than the fold's own"):
        sm._base_factors(fold, ts.finish_momentum(fold["raw"], 0.5))


def test_form_a_refuses_a_transform_handed_factors_other_than_the_folds_own(sm, monkeypatch):
    """A fresh finish at 0.5 is equal in value but another object: form A's transform must refuse it."""
    _, _, ctx, *_ = _world()
    fold = sm.fold_factors(ctx)
    fresh = ts.finish_momentum(fold["raw"], 0.5)
    assert sm.same_factors(fresh, fold["finished"][0.5]) and fresh is not fold["finished"][0.5]
    monkeypatch.setitem(harness._G_MOM_PREP, "fold", {**harness._G_MOM_PREP["fold"], "factors": fresh})
    with pytest.raises(RuntimeError, match="other than the fold's own"):
        sm.form_a_arms(ctx)


def _records_world():
    genes = ["a", "b", "h"]
    half = np.array([[ts.CLIP_AT, -ts.CLIP_AT, ts.CLIP_AT], [0.0, 0.06, 0.0]])
    one = np.array([[ts.CLIP_AT, -ts.CLIP_AT, ts.CLIP_AT], [0.08, 0.12, 0.2]])

    def sf(logf):
        cf = ts.CategoryFactors(category="cm", index={"c0": 0, "c1": 1}, sizes=np.array([3, 4]), logf=logf)
        return ts.StateFactors(genes=genes, metacells_on="later", categories={"cm": cf}, skipped={})

    return {0.5: sf(half), 1.0: sf(one)}


def test_the_factor_records_count_the_clip_exactly_and_skip_the_held_genes(sm):
    """Entries exactly at CLIP_AT count as at the clip; the held gene h (factor 1 when applied) is left out; the share
    newly not zeroed and the mean over all entries are reported."""
    rec = sm.factor_records(_records_world(), ["h"])
    half, one = rec["0.5"]["cm"], rec["1"]["cm"]
    assert half["share_at_clip"] == 0.5 and one["share_at_clip"] == 0.5
    assert half["share_not_zeroed"] == 0.75 and one["share_not_zeroed"] == 1.0
    assert half["share_unzeroed_vs_0.5"] == 0.0 and one["share_unzeroed_vs_0.5"] == 0.25
    assert half["mean_abs_log_all_entries"] == pytest.approx((2 * ts.CLIP_AT + 0.06) / 4)
    assert half["mean_abs_log_not_zeroed"] == pytest.approx((2 * ts.CLIP_AT + 0.06) / 3)
    assert one["mean_abs_log_all_entries"] == pytest.approx((2 * ts.CLIP_AT + 0.2) / 4)
    with_held = sm.factor_records(_records_world(), [])
    assert with_held["0.5"]["cm"]["share_at_clip"] == pytest.approx(3 / 6)
    assert with_held["1"]["cm"]["share_unzeroed_vs_0.5"] == pytest.approx(2 / 6)


def test_the_labels_say_no_keep_and_the_limits_are_the_docstrings(sm):
    """No text calls a pass a keep or form A the adoption rule's read; the limits are the docstring's, word for word,
    in the order written, and the reading writes them out."""
    assert "not a keep" in sm.PASSES and "passes the adoption rule" not in sm.PASSES
    assert "the adoption rule's read" not in sm.LABEL_T49
    assert "in both forms a labelled in-sample diagnostic: form A selects a constant on the cuts that kept g_mom" in (
        sm.LABEL_T49)
    doc = " ".join(sm.__doc__.split())
    assert len(sm.LIMITS) == 14 and all(x in doc for x in sm.LIMITS)
    assert sm.LIMITS[0].startswith("T49 chooses a value of one constant (B0's damping applied to g_mom)")
    assert sm.LIMITS[1].endswith("No outcome is read as 'g_mom was under-sized by a factor of about 3'.")
    assert sm.LIMITS[-1].endswith("and the change of sample kind in the truth.")
    assert "the adoption rule's read" not in doc and "So what g_mom delivers is about 0.3" not in doc


def test_the_reading_writes_the_limits(sm, tmp_path):
    out, t11, t42 = _write(sm, tmp_path, {k: _task(sm, {}) for k in sm.CUT_FILES})
    res = sm.read(out, t11, t42)
    assert res["limits"] == list(sm.LIMITS)
    text = (out / "t49.txt").read_text()
    assert all(f"- {x}" in text for x in sm.LIMITS) and "PASSES" not in text


def test_the_fold_is_refused_unless_g_mom_factors_finishes_once_at_half(sm, monkeypatch):
    _, _, ctx, *_ = _world()
    monkeypatch.setattr(harness, "B0_DAMPING", 0.6)
    with pytest.raises(RuntimeError, match="not once at 0.5"):
        sm.fold_factors(ctx)


# ---------------------------------------------------------------------------------------------------------------------
# The arms


def test_the_base_arms_are_the_registered_builds_bit_for_bit(sm):
    _, _, ctx, groups, dirs, held = _world()
    a = sm.form_a_arms(ctx)
    b = sm.form_b_arms(ctx, dirs, groups, 0.2, held, held)
    assert set(a) == {"g_mom", "g_mom_d075", "g_mom_d1"} and set(b) == {"m1_r_gmom", "m1_r_gmom_g075", "m1_r_gmom_g1"}
    assert harness._g_mom_factors is not None and harness._g_mom_factors.__name__ == "_g_mom_factors"  # restored
    harness._G_MOM_PREP.clear()
    gm = harness._g_mom(ctx)
    # the fit try T11's own build makes (inputs classified by b0_population) is the fold's, bit for bit
    assert sm.same_factors(harness._G_MOM_PREP["fold"]["factors"], sm._FOLD["fold"]["finished"][0.5])
    harness._G_MOM_PREP.clear()
    m1 = sm.gc.m1_r_arms(ctx, dirs, groups, 0.2, held, held)["m1_r_gmom"]
    assert a["g_mom"].X.dtype == gm.X.dtype and a["g_mom"].X.tobytes() == gm.X.tobytes()
    assert b["m1_r_gmom"].X.dtype == m1.X.dtype and b["m1_r_gmom"].X.tobytes() == m1.X.tobytes()
    for base, larger in (("g_mom", ("g_mom_d075", "g_mom_d1")), ("m1_r_gmom", ("m1_r_gmom_g075", "m1_r_gmom_g1"))):
        arms = {**a, **b}
        for x in larger:
            assert arms[x].X.shape == arms[base].X.shape and not np.array_equal(arms[x].X, arms[base].X)


def _spy_apply(monkeypatch):
    calls = []
    real = ts.apply_state_factors

    def spy(pop, sf, held_flat, seed, **kw):
        calls.append({"counts": pop.layers["counts"].toarray().tobytes() if hasattr(pop.layers["counts"], "toarray")
                      else np.asarray(pop.layers["counts"]).tobytes(),
                      "cells": list(pop.obs_names), "source": list(pop.obs["source_obs"].astype(str)),
                      "held": list(held_flat), "seed": list(seed), "rescale": kw.get("rescale"),
                      "no_edit": np.asarray(kw.get("no_edit")).tobytes(), "sf": sf})
        return real(pop, sf, held_flat, seed, **kw)

    monkeypatch.setattr(ts, "apply_state_factors", spy)
    return calls


@pytest.mark.parametrize("form", ["A", "B"])
def test_the_larger_arms_differ_from_their_base_in_the_finished_factors_only(sm, monkeypatch, form):
    _, _, ctx, groups, dirs, held = _world()
    fold = sm.fold_factors(ctx)
    calls = _spy_apply(monkeypatch)
    fits = []
    real_fit = ts.fit_state_factors
    monkeypatch.setattr(ts, "fit_state_factors", lambda *a, **k: fits.append(1) or real_fit(*a, **k))
    if form == "A":
        arms = sm.form_a_arms(ctx)
    else:
        arms = sm.form_b_arms(ctx, dirs, groups, 0.2, held, held)
    assert not fits  # the fold's fit is the one object every arm uses
    assert len(calls) == 3
    base = calls[0]
    for call, d in zip(calls, (0.5, 0.75, 1.0), strict=True):
        assert call["sf"] is fold["finished"][d]
        for key in ("counts", "cells", "source", "held", "seed", "rescale", "no_edit"):
            assert call[key] == base[key], f"a larger arm differs from its base in {key}"
    assert base["seed"] == [ctx.seed, 5] and base["rescale"] is False
    assert len({a.X.tobytes() for a in arms.values()}) == 3


# ---------------------------------------------------------------------------------------------------------------------
# The reading


def _arm(score, mmd_u=0.1):
    return {"score": score, "groups": {"distribution": score / 100},
            "metrics": {"mmd_u": mmd_u, "variogram": 0.01},
            "skills": {"de_score": score / 200, "de_direction": 0.4, "mmd_u": 0.3, "variogram": 0.2}}


B0, GM, M1 = 50.0, 52.0, 53.0
LARGER = {"g_mom_d075": "g_mom", "g_mom_d1": "g_mom", "m1_r_gmom_g075": "m1_r_gmom", "m1_r_gmom_g1": "m1_r_gmom"}


def _task(sm, gains, seeds=range(10), commit="c" * 40, clean=True, task=None, dampings=None, shift=None,
          identical=False):
    """A task file. gains: larger arm -> mean gain over its base, or a list of ten per-seed gains; shift: base arm ->
    an offset from its committed score. A larger arm's raw metric differs from its base's by 1e-9 unless identical."""
    shift = shift or {}
    rows = []
    for i, s in enumerate(seeds):
        arms = {"b0": _arm(B0 + shift.get("b0", 0.0)), "g_mom": _arm(GM + shift.get("g_mom", 0.0)),
                "m1_r_gmom": _arm(M1 + shift.get("m1_r_gmom", 0.0))}
        for a, base in LARGER.items():
            g = gains.get(a, 0.0)
            g = g[i] if isinstance(g, list) else g
            arms[a] = _arm(arms[base]["score"] + g, mmd_u=0.1 if identical else 0.1 + 1e-9)
        rows.append({"seed": s, "floor": 50.0, "arms": arms})
    return {"label": "synthetic", "task": task, "dampings": dampings or dict(sm.DAMPINGS), "seeds": list(seeds),
            "reproduction_tolerance": 1e-9, "factor_records": {"0.5": {"cm": {"share_not_zeroed": 0.2,
                                                                               "share_at_clip": 0.0,
                                                                               "mean_abs_log_not_zeroed": 0.07,
                                                                               "mean_abs_log_all_entries": 0.014,
                                                                               "share_unzeroed_vs_0.5": 0.0}}},
            "rows": rows, "stamp": {"commit": commit, "tree_clean": clean, "utc": "x", "inputs": {}}}


def _write(sm, tmp_path, cuts, drop=()):
    """cuts: cut key -> task record; T11's and T42's committed files, synthetic, beside them (drop: (kind, cut))."""
    out = tmp_path / "out"
    (out / "tasks").mkdir(parents=True)
    for k, rec in cuts.items():
        rec = {**rec, "task": k if rec.get("task") is None else rec["task"]}
        (out / "tasks" / sm.CUT_FILES[k]).write_text(json.dumps(rec))
    t11, t42 = tmp_path / "t11", tmp_path / "t42"
    t11.mkdir()
    t42.mkdir()
    for k, f in sm.CUT_FILES.items():
        if ("t11", k) not in drop:
            (t11 / f).write_text(json.dumps({"rows": [{"seed": s, "scores": {"floor": 50.0, "b0": B0, "g_mom": GM}}
                                                      for s in range(10)]}))
        if ("t42", k) not in drop:
            (t42 / f).write_text(json.dumps({"rows": [{"seed": s, "arms": {"b0": _arm(B0), "m1_r_gmom": _arm(M1)}}
                                                      for s in range(10)]}))
    return out, t11, t42


def _read(sm, tmp_path, gains=None, drop=(), **over):
    gains = gains or {}
    cuts = {k: _task(sm, gains) for k in sm.CUT_FILES}
    for k, kw in over.items():
        if kw is None:
            del cuts[getattr(sm.gc, k)]
        else:
            cuts[getattr(sm.gc, k)] = _task(sm, {**gains, **kw.pop("gains", {})}, **kw)
    out, t11, t42 = _write(sm, tmp_path, cuts, drop)
    return sm.read(out, t11, t42)["reading"]


@pytest.mark.parametrize("gains, kind, arm", [
    ({"g_mom_d075": 0.5, "g_mom_d1": 1.0}, "passes", "g_mom_d1"),
    ({"g_mom_d075": 1.0, "g_mom_d1": 1.2}, "passes", "g_mom_d075"),
    ({"g_mom_d075": 0.1, "g_mom_d1": 0.35}, "passes", "g_mom_d1"),
    ({"g_mom_d075": -0.1, "g_mom_d1": -0.2}, "loses", None),
    ({"g_mom_d075": 0.1, "g_mom_d1": -0.2}, "not separated", None),
    ({}, "not separated", None),
    ({"g_mom_d075": 0.5, "g_mom_d1": 1.0, "m1_r_gmom_g1": -0.31}, "vetoed", "g_mom_d1"),
    ({"g_mom_d075": 1.0, "g_mom_d1": 1.2, "m1_r_gmom_g1": -0.5}, "passes", "g_mom_d075"),
    ({"g_mom_d075": 1.0, "g_mom_d1": 1.2, "m1_r_gmom_g075": -0.4}, "vetoed", "g_mom_d075"),
    ({"g_mom_d075": 1.0, "g_mom_d1": 1.2, "m1_r_gmom_g075": -0.01}, "passes", "g_mom_d075"),
])
def test_each_statement_is_reached_by_a_table_built_for_it(sm, tmp_path, gains, kind, arm):
    r = _read(sm, tmp_path, gains)
    want_text = {"passes": f"{sm.PASSES}: {arm}", "vetoed": f"{sm.VETOED} ({arm})"}
    if kind in want_text:
        assert r["outcome"] == want_text[kind]
    else:
        assert r["outcome"].startswith({"loses": sm.LOSES, "not separated": sm.NOT_SEPARATED}[kind] + " (")
    assert r["valid"] and r["statement"]["arm"] == arm
    assert set(r["form_a"]) == {"A1", "A2", "A3"} and set(r["form_b_diagnostic"]) == {"B1", "B2", "B3"}
    a3 = r["form_a"]["A3"]
    assert a3["base"] == "g_mom_d075" and a3["change"] == "g_mom_d1"
    want = gains.get("g_mom_d1", 0.0) - gains.get("g_mom_d075", 0.0)
    cut = a3["cuts"][sm.gc.TASK_CHAMBER]
    assert cut["mean"] == pytest.approx(want) and cut["range"] == pytest.approx([want, want])
    assert cut["skills"]["de_score"] == pytest.approx(want / 200)
    assert cut["points"]["de_score"] == pytest.approx(25.0 * want / 200)
    assert r["mean_scores"][sm.gc.TASK_WHOLE]["g_mom_d1"] == pytest.approx(GM + gains.get("g_mom_d1", 0.0))
    assert r["factor_records"][sm.gc.TASK_WHOLE]["0.5"]["cm"]["share_not_zeroed"] == 0.2


def test_the_order_rule_between_the_two_larger_arms(sm, tmp_path):
    """Both pass against g_mom: g_mom_d1 only if it also passes against g_mom_d075; and g_mom_d1 alone, when
    g_mom_d075 fails against g_mom and g_mom_d1 fails against g_mom_d075, is still put forward on the last branch."""
    r = _read(sm, tmp_path, {"g_mom_d075": 0.5, "g_mom_d1": 1.0})
    assert r["statement"]["basis"] == "passes against g_mom and against g_mom_d075"
    r = _read(sm, tmp_path / "b", {"g_mom_d075": 1.0, "g_mom_d1": 1.2})
    assert r["form_a"]["A2"]["passes"] and not r["form_a"]["A3"]["passes"]
    assert r["statement"]["arm"] == "g_mom_d075" and r["statement"]["basis"] == "passes against g_mom"
    r = _read(sm, tmp_path / "c", {"g_mom_d075": 0.1, "g_mom_d1": 0.35})
    assert not r["form_a"]["A1"]["passes"] and r["form_a"]["A2"]["passes"] and not r["form_a"]["A3"]["passes"]
    assert r["statement"]["basis"] == "alone passes against g_mom"


def test_the_veto_reads_the_twin_on_both_cuts(sm, tmp_path):
    r = _read(sm, tmp_path, {"g_mom_d075": 0.5, "g_mom_d1": 1.0}, TASK_WHOLE={"gains": {"m1_r_gmom_g1": -0.31}})
    assert r["valid"] and r["outcome"].startswith(sm.VETOED) and r["statement"]["veto"]["loses"]
    r = _read(sm, tmp_path / "zero", {"g_mom_d075": 0.5, "g_mom_d1": 1.0, "m1_r_gmom_g1": 0.0})
    assert r["outcome"] == f"{sm.PASSES}: g_mom_d1" and not r["statement"]["veto"]["loses"]


@pytest.mark.parametrize("twin, vetoed", [(-0.01, False), (-0.05, False), (-0.3, False), (-0.31, True)])
def test_the_veto_needs_a_loss_beyond_the_adoption_rules_bound_not_a_noise_sign(sm, tmp_path, twin, vetoed):
    """A twin's mean of -0.01 is inside a ten-seed mean's noise (paired sd about 0.06 to 0.08); only a loss beyond
    vec.adoption.MAX_LOSS (0.3) on either cut vetoes."""
    assert sm.VETO_LOSS == 0.3
    r = _read(sm, tmp_path, {"g_mom_d075": 0.5, "g_mom_d1": 1.0}, TASK_CHAMBER={"gains": {"m1_r_gmom_g1": twin}})
    assert r["valid"] and r["statement"]["veto"]["loses"] is vetoed
    assert r["outcome"] == (f"{sm.VETOED} (g_mom_d1)" if vetoed else f"{sm.PASSES}: g_mom_d1")


def test_losing_needs_at_most_two_seeds_ahead_on_each_cut(sm, tmp_path):
    two_ahead = [0.05, 0.05] + [-0.3] * 8
    three_ahead = [0.06, 0.06, 0.06] + [-0.3] * 7
    r = _read(sm, tmp_path, {"g_mom_d075": two_ahead, "g_mom_d1": two_ahead})
    assert r["valid"] and r["outcome"] == f"{sm.LOSES} (A1 -0.23 / -0.23, A2 -0.23 / -0.23)"
    r = _read(sm, tmp_path / "more", {"g_mom_d075": two_ahead, "g_mom_d1": two_ahead},
              TASK_WHOLE={"gains": {"g_mom_d1": three_ahead}})
    assert r["valid"] and r["outcome"] == f"{sm.NOT_SEPARATED} (largest mean gain -0.19)"


def test_the_loses_and_not_separated_texts_carry_their_means(sm, tmp_path):
    """'loses' is a sign, not a size, and 'not separated' covers a steady gain below 0.3: both texts say the means."""
    r = _read(sm, tmp_path, {"g_mom_d075": -0.05, "g_mom_d1": -0.05}, TASK_WHOLE={"gains": {"g_mom_d1": -0.07}})
    assert r["outcome"] == f"{sm.LOSES} (A1 -0.05 / -0.05, A2 -0.05 / -0.07)"
    r = _read(sm, tmp_path / "steady", {"g_mom_d075": 0.25, "g_mom_d1": 0.25})
    assert r["form_a"]["A1"]["cuts"][sm.gc.TASK_WHOLE]["ahead"] == 10
    assert r["outcome"] == f"{sm.NOT_SEPARATED} (largest mean gain +0.25)"


def test_a_gain_of_exactly_point_three_does_not_pass(sm, tmp_path):
    r = _read(sm, tmp_path, {"g_mom_d075": 0.3, "g_mom_d1": 0.3})
    assert r["valid"] and not r["form_a"]["A1"]["passes"] and not r["form_a"]["A2"]["passes"]
    assert r["outcome"] == f"{sm.NOT_SEPARATED} (largest mean gain +0.30)"
    r = _read(sm, tmp_path / "more", {"g_mom_d1": 0.31})
    assert r["outcome"] == f"{sm.PASSES}: g_mom_d1"


def test_a_loss_beyond_point_three_on_one_cut_vetoes_the_rule(sm, tmp_path):
    r = _read(sm, tmp_path, {"g_mom_d1": 1.0}, TASK_WHOLE={"gains": {"g_mom_d1": -0.4}})
    assert r["valid"] and not r["form_a"]["A2"]["passes"] and r["outcome"].startswith(sm.NOT_SEPARATED)


@pytest.mark.parametrize("over, reason", [
    ({"TASK_WHOLE": None}, "missing cuts"),
    ({"TASK_CHAMBER": {"seeds": range(9)}}, "seeds"),
    ({"TASK_WHOLE": {"shift": {"b0": 1e-6}}}, "do not reproduce"),
    ({"TASK_WHOLE": {"shift": {"g_mom": 1e-6}}}, "do not reproduce"),
    ({"TASK_CHAMBER": {"shift": {"m1_r_gmom": 1e-6}}}, "do not reproduce"),
    ({"TASK_WHOLE": {"dampings": {"g_mom": 0.5, "g_mom_d075": 0.75, "g_mom_d1": 1.0, "m1_r_gmom": 0.5,
                                  "m1_r_gmom_g075": 0.75, "m1_r_gmom_g1": 0.5}}}, "wrong damping"),
    ({"TASK_CHAMBER": {"commit": "d" * 40}}, "2 commits"),
    ({"TASK_WHOLE": {"clean": False}}, "not clean"),
    ({"TASK_CHAMBER": {"task": "T1/E9.5->E10.5"}}, "the file holds task"),
    ({"TASK_CHAMBER": {"identical": True}, "TASK_WHOLE": {"identical": True}}, "exactly as their base"),
])
def test_the_void_cases(sm, tmp_path, over, reason):
    r = _read(sm, tmp_path, {"g_mom_d1": 1.0}, **over)
    assert not r["valid"] and r["outcome"] == "VOID"
    assert any(reason in x for x in r["validity"]["reasons"]), r["validity"]["reasons"]


@pytest.mark.parametrize("drop", [(("t11", "T1/E9.5->E10.5"),), (("t42", "T1/E9.5->E10.5:chamber"),)])
def test_a_missing_committed_file_voids(sm, tmp_path, drop):
    r = _read(sm, tmp_path, {"g_mom_d1": 1.0}, drop=drop)
    assert not r["valid"] and any("do not reproduce" in x for x in r["validity"]["reasons"])


def test_a_difference_within_the_tolerance_and_a_larger_arm_identical_on_one_cut_stay_valid(sm, tmp_path):
    r = _read(sm, tmp_path, {"g_mom_d1": 1.0}, TASK_WHOLE={"shift": {"g_mom": 1e-12, "m1_r_gmom": -1e-12}})
    assert r["valid"]
    assert r["validity"]["reproduction"][sm.gc.TASK_WHOLE]["largest_difference"]["g_mom"] <= 1e-9
    r = _read(sm, tmp_path / "one", {"g_mom_d1": 1.0}, TASK_WHOLE={"identical": True})
    assert r["valid"]


def test_the_reading_refuses_to_overwrite(sm, tmp_path):
    out, t11, t42 = _write(sm, tmp_path, {k: _task(sm, {}) for k in sm.CUT_FILES})
    sm.read(out, t11, t42)
    assert (out / "t49.json").exists() and (out / "t49.txt").exists()
    with pytest.raises(SystemExit, match="already holds"):
        sm.read(out, t11, t42)


def test_the_run_needs_its_row(sm, monkeypatch):
    monkeypatch.setattr(sm.n12, "tries_rows", lambda: {"T44", "T47"})
    monkeypatch.setattr(sm.gc, "n22_passed", lambda path: True)
    with pytest.raises(SystemExit, match="T49 has no row"):
        sm.main(["run", "--task", sm.gc.TASK_CHAMBER, "--out", "t49"])


def test_the_run_needs_n22(sm, monkeypatch):
    monkeypatch.setattr(sm.n12, "tries_rows", lambda: {"T49"})
    with pytest.raises(SystemExit, match="does not hold a passing N22"):
        sm.main(["run", "--task", sm.gc.TASK_CHAMBER, "--out", "t49", "--n22", "n22x"])


def test_the_run_refuses_an_existing_task_file(sm, monkeypatch):
    monkeypatch.setattr(sm.n12, "tries_rows", lambda: {"T49"})
    monkeypatch.setattr(sm.gc, "n22_passed", lambda path: True)
    tasks = sm.RESULTS / "t49" / "tasks"
    tasks.mkdir(parents=True)
    (tasks / sm.CUT_FILES[sm.gc.TASK_WHOLE]).write_text("{}")
    with pytest.raises(SystemExit, match="already holds"):
        sm.main(["run", "--task", sm.gc.TASK_WHOLE, "--out", "t49"])
