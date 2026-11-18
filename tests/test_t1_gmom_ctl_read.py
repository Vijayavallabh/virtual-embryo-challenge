"""experiments/t1_m1/gmom_ctl_read.py (try T59 of experiments/harness/TRIES.md) on synthetic data: the controls'
factors at each damping, the arms, the reading's statements, its void cases and the run's refusals. No released file,
statistics file or score is read."""

from __future__ import annotations

import importlib.util
import json
import sys

import numpy as np
import pytest
from test_g_mom_verify import _small_inputs
from test_t1_state_transport import _gmom_ctx

from vec import harness
from vec import t1_state_transport as ts
from vec.boards import REPO_ROOT
from vec.harness import CANDIDATES

SCRIPT = REPO_ROOT / "experiments" / "t1_m1" / "gmom_ctl_read.py"


@pytest.fixture
def cm(monkeypatch, tmp_path):
    spec = importlib.util.spec_from_file_location("gmom_ctl_read_module", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    monkeypatch.setattr(mod, "RESULTS", tmp_path / "results")
    mod._CTL.clear()
    mod.gs._FOLD.clear()
    harness._G_MOM_PREP.clear()
    harness._G_MOM_NULL_PREP.clear()
    return mod


def _ctx(seed=1000000):
    ref, start, _ = _small_inputs(7)
    return _gmom_ctx(ref, start, seed=seed)


# ---------------------------------------------------------------------------------------------------------------------
# The controls' factors and the arms


def test_the_arms_and_their_dampings(cm):
    assert cm.CUT_ARMS == {"g_mom": 0.5, "g_mom_d075": 0.75, "g_mom_d1": 1.0, "g_mom_null": 0.5,
                           "g_mom_null_d075": 0.75, "g_mom_null_d1": 1.0, "g_mom_perm": 0.5,
                           "g_mom_perm_d075": 0.75, "g_mom_perm_d1": 1.0}
    assert cm.WORLD_DAMPINGS == {"m1_r_gmom": 0.5, "m1_r_gmom_g075": 0.75, "m1_r_gmom_g1": 1.0}
    assert cm.WORLDS == ("lib3", "lib3_t67", "lib3_t33", "ref67", "ref50")
    assert cm.NULL_SEED == cm.PERM_SEED == 20261110 and cm.REPRO_TOL == 1e-9


def test_at_half_the_twin_and_the_permuted_arm_are_the_registered_candidates_builds_bit_for_bit(cm):
    ctx = _ctx()
    b0 = harness._b0_t1(harness.B0_DAMPING)
    arms = cm.cut_arms(ctx, b0)
    assert set(arms) == set(cm.CUT_ARM_NAMES)
    ctl = cm._CTL["fold"]
    assert ctl["g_mom_null"][0.5] is harness._G_MOM_NULL_PREP["fold"]["factors"]
    # the registered builds, every cache cleared (their own fits made afresh)
    harness._G_MOM_PREP.clear()
    harness._G_MOM_NULL_PREP.clear()
    reg = {"g_mom_null": CANDIDATES["g_mom_null"].build(ctx), "g_mom_perm": CANDIDATES["g_mom_perm"].build(ctx),
           "g_mom": harness._g_mom(ctx), "b0": b0(ctx)}
    for a, x in reg.items():
        assert arms[a].X.dtype == x.X.dtype and arms[a].X.tobytes() == x.X.tobytes(), a
    assert cm.gs.same_factors(harness._G_MOM_NULL_PREP["fold"]["factors"], ctl["g_mom_null"][0.5])


def _spy_apply(monkeypatch):
    calls = []
    real = ts.apply_state_factors

    def spy(pop, sf, held_flat, seed, **kw):
        c = pop.layers["counts"]
        calls.append({"counts": (c.toarray() if hasattr(c, "toarray") else np.asarray(c)).tobytes(),
                      "cells": list(pop.obs_names), "source": list(pop.obs["source_obs"].astype(str)),
                      "held": list(held_flat), "seed": list(seed), "rescale": kw.get("rescale"),
                      "no_edit": np.asarray(kw.get("no_edit")).tobytes(), "sf": sf})
        return real(pop, sf, held_flat, seed, **kw)

    monkeypatch.setattr(ts, "apply_state_factors", spy)
    return calls


def test_each_control_arm_is_the_earlier_trys_build_with_only_the_dampings_finish_changed(cm, monkeypatch):
    ctx = _ctx()
    fold = cm.gs.fold_factors(ctx)
    finishes = []
    real_fnd = ts.finish_no_development

    def spy_fnd(first, second, target, held_flat, reference_damping):
        finishes.append((first, second, target, list(held_flat), reference_damping))
        return real_fnd(first, second, target, held_flat, reference_damping)

    monkeypatch.setattr(ts, "finish_no_development", spy_fnd)
    ctl = cm.control_factors(ctx)
    monkeypatch.setattr(ts, "finish_no_development", real_fnd)
    # one registered call at 0.5 on g_mom's own factors (the spy of the script saw it), then one per larger damping,
    # on the same two raw directions, the target g_mom finished at that damping
    assert len(finishes) == 4
    first, second = finishes[0][0], finishes[0][1]
    assert finishes[0][2] is fold["finished"][0.5] and finishes[0][4] == 0.5
    for (f, s, target, held, ref_d), d in zip(finishes[2:], (0.75, 1.0), strict=True):
        assert f is first and s is second and held == finishes[0][3]
        assert target is fold["finished"][d] and ref_d == d
    free = np.array([g not in set(ctx.derived["t1_held_flat"]) for g in fold["finished"][0.5].genes])
    for d in (0.75, 1.0):
        assert cm.gs.same_factors(ctl["g_mom_null"][d], real_fnd(first, second, fold["finished"][d],
                                                                 finishes[0][3], d))
        assert cm.gs.same_factors(ctl["g_mom_perm"][d], ts.permute_factors(fold["finished"][d], cm.PERM_SEED))
        assert not cm.gs.same_factors(ctl["g_mom_null"][d], ctl["g_mom_null"][0.5])
        for c, cf in ctl["g_mom_null"][d].categories.items():  # the twin's size is g_mom's at that damping
            goal = fold["finished"][d].categories[c]
            assert ts._factor_rms(cf.logf, cf.sizes, free) == pytest.approx(
                ts._factor_rms(goal.logf, goal.sizes, free), rel=0.1)
    assert ctl["g_mom_perm"][0.75].record["permutation"] == ctl["g_mom_perm"][0.5].record["permutation"]
    # the arms: every control call shares g_mom's population, held genes, stream and no-edit mask; only sf differs
    calls = _spy_apply(monkeypatch)
    g = cm.gs.form_a_arms(ctx, {"g_mom": 1.0})
    arms = cm.control_arms(ctx)
    assert len(calls) == 7 and set(arms) == {f"{f}{s}" for f in ("g_mom_null", "g_mom_perm")
                                             for s in ("", "_d075", "_d1")}
    base = calls[0]
    want = [ctl[f][d] for f in ("g_mom_null", "g_mom_perm") for d in (0.5, 0.75, 1.0)]
    for call, sf in zip(calls[1:], want, strict=True):
        assert call["sf"] is sf
        for key in ("counts", "cells", "source", "held", "seed", "rescale", "no_edit"):
            assert call[key] == base[key], f"a control arm differs from g_mom in {key}"
    assert base["seed"] == [ctx.seed, 5] and base["rescale"] is False
    assert len({a.X.tobytes() for a in arms.values()} | {g["g_mom"].X.tobytes()}) == 7


def test_the_fits_are_made_once_per_fold(cm, monkeypatch):
    fits = []
    real_fit = ts.fit_state_factors
    monkeypatch.setattr(ts, "fit_state_factors", lambda *a, **k: fits.append(1) or real_fit(*a, **k))
    b0 = harness._b0_t1(harness.B0_DAMPING)
    ctx = _ctx()
    cm.cut_arms(ctx, b0)
    n = len(fits)
    assert n == 3  # g_mom's fit and the twin's two cross-fit directions
    cm.cut_arms(_gmom_ctx(ctx.inputs["ref"], ctx.inputs["start"], seed=1000001), b0)
    assert len(fits) == n


def test_the_run_stops_unless_a_fresh_finish_at_half_is_the_registered_twins(cm, monkeypatch):
    ctx = _ctx()
    cm.gs.fold_factors(ctx)  # the fold's own check passes first
    monkeypatch.setattr(cm.gs, "same_factors", lambda a, b: False)
    with pytest.raises(RuntimeError, match="not _g_mom_null_factors' own"):
        cm.control_factors(ctx)


def test_the_records_are_target_free_and_cover_the_three_families(cm):
    ctx = _ctx()
    ctl = cm.control_factors(ctx)
    rec = cm.control_records(ctl, ctx.derived["t1_held_flat"])
    assert set(rec["factors"]) == set(cm.FAMILIES) and set(rec["sizes"]) == set(cm.FAMILIES)
    for fam in cm.FAMILIES:
        assert set(rec["factors"][fam]) == {"0.5", "0.75", "1"}
    assert rec["permutation_same_at_every_damping"] is True
    assert rec["permutation"]["0.5"]["seed"] == cm.PERM_SEED
    assert set(rec["null_match"]["1"]) == set(ctl["g_mom_null"][1.0].categories)
    for d in ("0.5", "0.75", "1"):
        assert set(rec["null_match"][d]) == set(ctl["g_mom_null"][float(d)].categories)
        for m in rec["null_match"][d].values():  # the watched size match: g_mom's size inside the twin's range
            assert m["bracketed"] is True and m["damping"] is not None
            assert m["rms_at_smallest_damping"] <= m["rms_target"] <= m["rms_at_largest_damping"]
    assert cm.size_match_reasons("x", rec["null_match"]) == []
    for c, m in rec["null_match"]["1"].items():
        assert m["reference_damping"] == 1.0 and m["rms_target"] == pytest.approx(rec["sizes"]["g_mom"]["1"][c])


# ---------------------------------------------------------------------------------------------------------------------
# The reading on synthetic score tables


def _arm(score, mmd_u=0.1):
    return {"score": float(score), "groups": {"distribution": score / 100},
            "metrics": {"mmd_u": mmd_u, "variogram": 0.01},
            "skills": {"de_score": score / 200, "de_direction": 0.4, "mmd_u": 0.3, "variogram": 0.2}}


CUT_BASE = {"b0": 50.0, "g_mom": 52.0, "g_mom_null": 51.0, "g_mom_perm": 51.5}
WORLD_BASE = {"b0": 50.0, "m1_r_gmom": 53.0}


CH, WH = "T1/E9.5->E10.5:chamber", "T1/E9.5->E10.5"
CH_KEY = CH


def _gain(gains, arm, i):
    g = gains.get(arm, 0.0)
    return g[i] if isinstance(g, list) else g


def _match(**over):
    """A twin's size-match record that passes: per damping, one category matched inside its range."""
    out = {d: {"cm": {"damping": float(d) * 1.1, "relative_error": 1e-12, "rms_target": 0.1, "bracketed": True,
                      "rms_at_smallest_damping": 1e-4, "rms_at_largest_damping": 5.0}} for d in ("0.5", "0.75", "1")}
    for d, m in over.items():
        out[d]["cm"].update(m)
    return out


def _cut(cm, gains, seeds=range(10), commit="c" * 40, clean=True, task=None, shift=None, copy=None, dampings=None,
         **top):
    """A part-1 task file. gains: larger arm -> gain over its family's arm at 0.5 (a number or ten); shift: arm ->
    offset added after the committed files are made; copy: arm -> arm whose record it takes exactly."""
    shift, copy = shift or {}, copy or {}
    rows = []
    for i, s in enumerate(seeds):
        arms = {"b0": _arm(CUT_BASE["b0"])}
        for j, f in enumerate(cm.FAMILIES):
            arms[f] = _arm(CUT_BASE[f])
            for n, suffix in enumerate(("_d075", "_d1")):
                a = f + suffix
                arms[a] = _arm(CUT_BASE[f] + _gain(gains, a, i), mmd_u=0.1 + 1e-9 * (1 + n + 2 * j))
        for a, x in shift.items():
            arms[a] = {**arms[a], "score": arms[a]["score"] + x}
        for a, b in copy.items():
            arms[a] = dict(arms[b])
        rows.append({"seed": s, "floor": 50.0, "arms": arms})
    rec = {"label": "synthetic", "part": "cut", "task": task, "seeds": list(seeds),
           "dampings": dampings or dict(cm.CUT_ARMS), "null_seed": cm.NULL_SEED, "perm_seed": cm.PERM_SEED,
           "at_half": {f: {"used": "a" * 64, "registered_cache": "a" * 64} for f in cm.FAMILIES},
           "registered_rebuild": {"seed": 0, **{a: {"equal": True} for a in cm.T20_ARMS}},
           "records": {"permutation_same_at_every_damping": True, "permutation": {"0.5": {"seed": cm.PERM_SEED}},
                       "null_categories": {"0.5": ["cm"], "0.75": ["cm"], "1": ["cm"]}, "factors": {}, "sizes": {},
                       "null_match": _match()},
           "rows": rows, "stamp": {"commit": commit, "tree_clean": clean, "utc": "x", "inputs": {}}}
    rec.update(top)
    return rec


def _world(cm, k, gains, seeds=range(10), commit="c" * 40, clean=True, shift=None, thinning=None, sha_off=False,
           offset=None, copy=None, sha_same=None, **top):
    """A part-2 task file of world k. The thinned references' scores sit 0.5 above lib3's (their thinning is not
    idle) unless offset says otherwise."""
    shift = shift or {}
    off = (0.5 if k in cm.T54_WORLDS else 0.0) if offset is None else offset
    rows = []
    for i, s in enumerate(seeds):
        arms = {"b0": _arm(WORLD_BASE["b0"] + off), "m1_r_gmom": _arm(WORLD_BASE["m1_r_gmom"] + off)}
        for n, a in enumerate(("m1_r_gmom_g075", "m1_r_gmom_g1")):
            arms[a] = _arm(WORLD_BASE["m1_r_gmom"] + off + _gain(gains.get(k, {}), a, i), mmd_u=0.1 + 1e-9 * (n + 1))
        for a, x in shift.items():
            arms[a] = {**arms[a], "score": arms[a]["score"] + x}
        for a, b in (copy or {}).items():
            arms[a] = dict(arms[b])
        sha = {a: f"{a}-{s}" + ("x" if sha_off else "") for a in arms}
        if sha_same is not None and sha_same[1] == s:
            sha[sha_same[0]] = sha[cm.M1]
        rows.append({"seed": s, "floor": 50.0, "arms": arms, "sha256": sha})
    if k == cm.HEART:
        th = thinning
    else:
        spec = cm.THINNING_OF[k]
        th = {"p": spec["p"], "stream": list(spec["stream"]), "depth_ratio": spec["p"], "any_cell_emptied": False,
              "rebuild_largest_difference": 0.0}
        th.update(thinning or {})
    rec = {"label": "synthetic", "part": "world", "task": cm.TASK, "world": k,
           "truth": k if k in cm.T50_WORLDS else cm.HEART, "truth_rows": 100, "truth_contamination": 5,
           "fraction": 0.2, "seeds": list(seeds), "thinning": th, "dampings": dict(cm.WORLD_DAMPINGS), "rows": rows,
           "stamp": {"commit": commit, "tree_clean": clean, "utc": "x", "inputs": {}}}
    rec.update(top)
    return rec


def _committed_cut(cm, rec):
    rows = rec["rows"]
    t49 = {"rows": [{"seed": r["seed"], "arms": {a: {"score": r["arms"][a]["score"]} for a in cm.T49_ARMS}}
                    for r in rows]}
    summ = {a: [float(np.mean([r["arms"][a]["score"] for r in rows])),
                float(np.std([r["arms"][a]["score"] for r in rows], ddof=1))] for a in cm.T20_ARMS}
    dec = {a: {"mean_gain": float(np.mean([r["arms"][a]["score"] - r["arms"]["b0"]["score"] for r in rows])),
               "seeds_ahead": int(sum(r["arms"][a]["score"] > r["arms"]["b0"]["score"] for r in rows))}
           for a in cm.T20_ARMS}
    t20_rows = {"rows": [{"seed": r["seed"], "scores": {a: r["arms"][a]["score"] for a in ("b0", *cm.T20_ARMS)}}
                         for r in rows]}
    return t49, summ, dec, t20_rows


def _read(cm, tmp_path, gains=None, world_gains=None, cut_over=None, world_over=None, drop=(), t20_rows=False,
          t20_shift=None, t20_ahead_off=0):
    """Write both cuts and the five worlds (with overrides: key -> _cut / _world keyword arguments, or None to leave
    the file out) and the committed files made from the unshifted tables; read them."""
    gains, world_gains = gains or {}, world_gains or {}
    cut_over, world_over = cut_over or {}, world_over or {}
    out = tmp_path / "out"
    (out / "tasks").mkdir(parents=True)
    t49d, t50d, t54d, t20r = (tmp_path / x for x in ("t49", "t50", "t54", "t20rows"))
    for d in (t49d, t50d, t54d, t20r):
        d.mkdir()
    t20 = {"summaries": {}, "decisions": [{"base": "b0", "change": a, "tasks": {}} for a in cm.T20_ARMS]}
    for k, f in cm.CUT_FILES.items():
        kw = cut_over.get(k, {})
        clean = _cut(cm, {**gains, **(kw or {}).get("gains", {})}, task=k)
        t49, summ, dec, rows = _committed_cut(cm, clean)
        if ("t49", k) not in drop:
            (t49d / f).write_text(json.dumps(t49))
        if ("t20", k) not in drop:
            t20["summaries"][k] = summ
            for d in t20["decisions"]:
                d["tasks"][k] = dec[d["change"]]
        for r in rows["rows"]:
            for a, x in (t20_shift or {}).items():
                r["scores"][a] += x
        (t20r / f).write_text(json.dumps(rows))
        if kw is None:
            continue
        kw = {x: v for x, v in kw.items() if x != "gains"}
        rec = _cut(cm, {**gains, **cut_over.get(k, {}).get("gains", {})}, **{"task": k, **kw})
        (out / "tasks" / f).write_text(json.dumps(rec))
    t20["decisions"][0]["tasks"].get(CH_KEY, {"seeds_ahead": 0})["seeds_ahead"] -= t20_ahead_off
    t20_file = tmp_path / "decisions.json"
    t20_file.write_text(json.dumps(t20))
    for k in cm.WORLDS:
        kw = world_over.get(k, {})
        clean = _world(cm, k, world_gains)
        folder = t50d if k in cm.T50_WORLDS else t54d
        if ("world", k) not in drop:
            (folder / f"{k}.json").write_text(json.dumps({"rows": clean["rows"]}))
        if kw is None:
            continue
        (out / "tasks" / f"{k}.json").write_text(json.dumps(_world(cm, k, world_gains, **kw)))
    res = cm.read(out, t49d, t20_file, t50d, t54d, t20_tasks=t20r if t20_rows else None)
    return res, out


def _all(value):
    return {k: {"m1_r_gmom_g1": value, "m1_r_gmom_g075": value} for k in ("lib3", "lib3_t67", "lib3_t33", "ref67",
                                                                         "ref50")}


NEEDS_GAINS = {"g_mom_d1": 1.0, "g_mom_null_d1": 0.2, "g_mom_perm_d1": 0.1, "g_mom_d075": 0.6,
               "g_mom_null_d075": 0.1, "g_mom_perm_d075": 0.05}


def test_statement_one_needs_development_and_the_assignment_and_the_values(cm, tmp_path):
    res, out = _read(cm, tmp_path, NEEDS_GAINS, _all(0.5))
    r = res["reading"]
    assert r["valid"], r["validity"]["reasons"]
    assert r["statement_1"]["statement"] == cm.EXCEEDS and r["statement_1"]["holds"]
    c = r["reading_1"][cm.gc.TASK_CHAMBER]["1"]
    assert c["G"]["mean"] == pytest.approx(1.0) and c["Z"]["mean"] == pytest.approx(0.2)
    assert c["D"]["mean"] == pytest.approx(0.8) and c["Q"]["mean"] == pytest.approx(0.1)
    assert c["A"]["mean"] == pytest.approx(0.9) and c["D"]["label"] == c["A"]["label"] == "ahead"
    assert c["D"]["seeds_ahead"] == 10 and c["D"]["range"] == pytest.approx([0.8, 0.8])
    assert c["D"]["points"]["de"] == pytest.approx(c["G"]["points"]["de"] - c["Z"]["points"]["de"])
    assert c["G"]["points"]["de"] > 0 and c["G"]["points"]["mmd_u_variogram"] == pytest.approx(0.0)
    assert r["reading_1"][cm.gc.TASK_WHOLE]["0.75"]["D"]["mean"] == pytest.approx(0.5)
    assert r["statement_2"]["statement"] == cm.STAYS
    assert r["mean_scores"]["ref50"]["m1_r_gmom_g1"] == pytest.approx(54.0)
    assert r["mean_scores"][cm.gc.TASK_CHAMBER]["g_mom_null_d1"] == pytest.approx(51.2)
    text = (out / "t59.txt").read_text()
    assert cm.EXCEEDS in text and all(f"Limit: {x}" in text for x in cm.LIMITS)


@pytest.mark.parametrize("over, named", [
    ({"g_mom_null_d1": 0.9}, ["D on T1/E9.5->E10.5:chamber level", "D on T1/E9.5->E10.5 level"]),
    ({"g_mom_null_d1": 1.5}, ["D on T1/E9.5->E10.5:chamber behind", "D on T1/E9.5->E10.5 behind"]),
    ({"g_mom_perm_d1": 0.9}, ["A on T1/E9.5->E10.5:chamber level", "A on T1/E9.5->E10.5 level"]),
])
def test_statement_one_not_shown_names_each_reading_that_is_level_or_behind(cm, tmp_path, over, named):
    r = _read(cm, tmp_path, {**NEEDS_GAINS, **over}, _all(0.5))[0]["reading"]
    s = r["statement_1"]["statement"]
    assert s.startswith(cm.NOT_SHOWN + " (") and all(n in s for n in named)
    assert len(s[len(cm.NOT_SHOWN) + 2:-1].split("; ")) == len(named)


def test_statement_one_reads_one_cut_alone(cm, tmp_path):
    r = _read(cm, tmp_path, NEEDS_GAINS, _all(0.5),
              cut_over={cm.gc.TASK_WHOLE: {"gains": {"g_mom_perm_d1": 0.95}}})[0]["reading"]
    assert r["statement_1"]["statement"] == f"{cm.NOT_SHOWN} (A on T1/E9.5->E10.5 level)"


def test_statement_one_ignores_the_three_quarter_damping(cm, tmp_path):
    r = _read(cm, tmp_path, {**NEEDS_GAINS, "g_mom_null_d075": 2.0}, _all(0.5))[0]["reading"]
    assert r["reading_1"][cm.gc.TASK_CHAMBER]["0.75"]["D"]["label"] == "behind"
    assert r["statement_1"]["statement"] == cm.EXCEEDS


def test_statement_two_stays_positive_with_no_size_threshold(cm, tmp_path):
    r = _read(cm, tmp_path, NEEDS_GAINS, _all(0.05))[0]["reading"]
    assert r["statement_2"]["statement"] == cm.STAYS
    assert r["reading_2"]["lib3"]["1"]["label"] == "level"


def test_statement_two_turns_names_the_world_and_its_side(cm, tmp_path):
    wg = _all(0.5)
    wg["ref50"]["m1_r_gmom_g1"] = -0.5
    wg["lib3_t33"]["m1_r_gmom_g1"] = [0.1] * 7 + [-0.1] * 3
    r = _read(cm, tmp_path, NEEDS_GAINS, wg)[0]["reading"]
    assert r["statement_2"]["statement"] == (f"{cm.TURNS} (ref50, {cm.DEEPER}: mean -0.50, 0 of 10 seeds ahead, "
                                             "negative)")
    assert r["statement_2"]["neither"] == ["lib3_t33"]


def test_statement_two_not_separated_names_every_world_meeting_neither(cm, tmp_path):
    wg = _all(0.5)
    wg["lib3_t33"]["m1_r_gmom_g1"] = [0.1] * 7 + [-0.1] * 3  # mean 0.04, 7 ahead
    wg["ref67"]["m1_r_gmom_g1"] = -0.2  # behind in sign, level by the test
    wg["lib3"]["m1_r_gmom_g075"] = -1.0  # the 0.75 damping is not read by statement 2
    r = _read(cm, tmp_path, NEEDS_GAINS, wg)[0]["reading"]
    assert r["statement_2"]["statement"] == (
        f"{cm.NOT_SEPARATED} (lib3_t33, {cm.SHALLOWER}: mean +0.04, 7 of 10 seeds ahead, positive; "
        f"ref67, {cm.DEEPER}: mean -0.20, 0 of 10 seeds ahead, negative)")
    assert r["reading_2"]["lib3"]["0.75"]["label"] == "behind"


def test_seven_of_ten_ahead_does_not_stay_positive(cm, tmp_path):
    wg = _all(0.5)
    wg["lib3"]["m1_r_gmom_g1"] = [1.0] * 7 + [-0.1] * 3
    r = _read(cm, tmp_path, NEEDS_GAINS, wg)[0]["reading"]
    assert r["statement_2"]["statement"] == (f"{cm.NOT_SEPARATED} (lib3, {cm.SHALLOWER}: mean +0.67, 7 of 10 seeds "
                                             "ahead, positive)")


def test_statement_one_needs_g_itself_ahead(cm, tmp_path):
    # G level (+0.2) while D and A read ahead because the controls lose: 'not shown', naming G on both cuts
    gains = {**NEEDS_GAINS, "g_mom_d1": 0.2, "g_mom_null_d1": -0.6, "g_mom_perm_d1": -0.6}
    r = _read(cm, tmp_path, gains, _all(0.5))[0]["reading"]
    assert r["valid"], r["validity"]["reasons"]
    c = r["reading_1"][CH]["1"]
    assert c["D"]["label"] == c["A"]["label"] == "ahead" and c["G"]["label"] == "level"
    assert r["statement_1"]["statement"] == f"{cm.NOT_SHOWN} (G on {CH} level; G on {WH} level)"
    assert not r["statement_1"]["holds"]
    # no gain at all (G 0, the controls behind): not shown either
    gains = {**NEEDS_GAINS, "g_mom_d1": 0.0, "g_mom_null_d1": -0.6, "g_mom_perm_d1": -0.6}
    r = _read(cm, tmp_path / "zero", gains, _all(0.5))[0]["reading"]
    assert r["statement_1"]["statement"].startswith(cm.NOT_SHOWN)


def test_statement_one_says_exceeds_and_prints_the_controls_share_of_the_gain(cm, tmp_path):
    # the twin reproduces 70% of g_mom's gain: D (+0.42) is ahead, the words say 'exceeds', and Z/G 0.70 is printed
    gains = {**NEEDS_GAINS, "g_mom_d1": 1.4, "g_mom_null_d1": 0.98, "g_mom_perm_d1": 0.2}
    res, out = _read(cm, tmp_path, gains, _all(0.5))
    r = res["reading"]
    assert r["statement_1"]["statement"] == cm.EXCEEDS
    share = r["statement_1"]["control_share_of_gain"][CH]
    assert share["Z/G"] == pytest.approx(0.7) and share["Q/G"] == pytest.approx(0.2 / 1.4)
    assert r["size_match_largest_relative_error"][CH] == {"0.5": 1e-12, "0.75": 1e-12, "1": 1e-12}
    text = (out / "t59.txt").read_text()
    assert f"{CH} Z/G 0.70, Q/G 0.14" in text and "'exceeds' is not 'needs'" in text
    assert "largest relative error per damping" in text
    assert cm.statement_one({k: {"1": {q: {"label": "ahead", "mean": 0.0} for q in "GZDQA"}}
                             for k in cm.CUT_FILES})["control_share_of_gain"][CH]["Z/G"] != 0  # nan where G is 0


def test_statement_two_names_a_small_steady_reversal_with_its_sign(cm, tmp_path):
    wg = _all(0.5)
    wg["ref67"]["m1_r_gmom_g1"] = -0.25  # every seed behind, smaller than the 0.3 bar
    r = _read(cm, tmp_path, NEEDS_GAINS, wg)[0]["reading"]
    assert r["statement_2"]["statement"] == (f"{cm.NOT_SEPARATED} (ref67, {cm.DEEPER}: mean -0.25, 0 of 10 seeds "
                                             "ahead, negative)")
    assert r["statement_2"]["signs"]["ref67"] == "negative" and r["statement_2"]["signs"]["lib3"] == "positive"


@pytest.mark.parametrize("over, ok", [
    ({}, True),
    ({"1": {"damping": None}}, False),
    ({"0.75": {"bracketed": False}}, False),
    ({"0.5": {"bracketed": None}}, False),
    ({"1": {"damping": None, "rms_target": 0.0, "bracketed": False}}, True),  # g_mom has no size: both zero
])
def test_a_failed_size_match_voids(cm, tmp_path, over, ok):
    r = _read(cm, tmp_path, NEEDS_GAINS, _all(0.5),
              cut_over={WH: {"records": {"permutation_same_at_every_damping": True,
                                         "permutation": {"0.5": {"seed": cm.PERM_SEED}},
                                         "null_categories": {"0.5": ["cm"], "0.75": ["cm"], "1": ["cm"]},
                                         "null_match": _match(**over)}}})[0]["reading"]
    assert r["valid"] is ok, r["validity"]["reasons"]
    if not ok:
        assert any("size match failed" in x for x in r["validity"]["reasons"])


def test_an_unrecorded_size_match_voids(cm):
    m = _match()
    del m["0.75"]
    assert "not recorded at every damping" in cm.size_match_reasons("x", m)[0]
    assert "not recorded at every damping" in cm.size_match_reasons("x", None)[0]
    assert cm.size_match_reasons("x", _match()) == []


def test_the_size_match_watch_sees_a_match_at_an_end_of_the_range(cm):
    mr = np.array([[0.2, -0.2, 0.1], [-0.2, 0.2, -0.1]])
    sizes, free = np.array([1.0, 2.0]), np.array([True, True, False])
    calls = []
    with cm._capture_size_match(calls):
        _, inside = ts._finish_to_size(mr, sizes, free, 0.3)
        _, outside = ts._finish_to_size(mr, sizes, free, 100.0)  # above what the clip allows: stops at the top end
    assert ts._finish_to_size is not None and [c["bracketed"] for c in calls] == [True, False]
    # once every entry sits at the clip the size stops rising, so the bisection keeps the first damping that reached
    # the cap: the damping found is not at an end of the range, and only the bracket shows the failed match
    lo, hi = ts.SIZE_MATCH_DAMPING
    assert lo < outside["damping"] < hi and outside["relative_error"] > 0.5
    assert calls[1]["rms_at_largest_damping"] < 100.0 and inside["relative_error"] < 1e-6


def test_the_at_half_record_is_the_registered_builds_sha256(cm):
    ctx = _ctx()
    cm.cut_arms(ctx, harness._b0_t1(harness.B0_DAMPING))
    ctl = cm._CTL["fold"]
    rec = cm.at_half_record(ctl)
    assert rec["null_seed"] == rec["perm_seed"] == 20261110
    for f in cm.FAMILIES:
        assert rec["sha256"][f]["used"] == rec["sha256"][f]["registered_cache"]
        assert len(rec["sha256"][f]["used"]) == 64 and "same objects" in rec["sha256"][f]["note"]
    assert len({rec["sha256"][f]["used"] for f in cm.FAMILIES}) == 3
    assert cm.factors_sha256(ctl["g_mom"][1.0]) != rec["sha256"]["g_mom"]["used"]
    # a twin at 0.5 other than the registered one, or no registered build, is recorded as different
    swapped = {**ctl, "g_mom_null": {**ctl["g_mom_null"], 0.5: ctl["g_mom_null"][1.0]}}
    assert (cm.at_half_record(swapped)["sha256"]["g_mom_null"]["used"]
            != rec["sha256"]["g_mom_null"]["registered_cache"])
    harness._G_MOM_NULL_PREP.clear()
    assert cm.at_half_record(ctl)["sha256"]["g_mom_null"]["registered_cache"] is None


def test_in_a_thinned_reference_world_the_arms_are_built_on_the_thinned_inputs(cm, monkeypatch):
    """cmd_run_world with every loader stubbed: the inputs world_route returns, not the inputs as loaded, reach every
    seed's build (try T54's arms fit on the thinned reference)."""
    from types import SimpleNamespace as NS

    import anndata as ad

    def adata(tag):
        a = ad.AnnData(np.zeros((2, 1), dtype=np.float32))
        a.obs["celltype"] = ["x", "y"]
        a.uns["tag"] = tag
        return a

    released = {"ref": adata("released"), "start": adata("released")}
    thinned = {"ref": adata("thinned"), "start": adata("thinned")}
    t = NS(key=cm.TASK, panel_board="pb", reference="r_stage", target="t_stage",
           inputs={"ref": "ref.h5ad", "start": "start.h5ad"})
    routed = []

    def world_route(world, t_, data, genes, inputs, role, ref_stage, read_full):
        routed.append(world)
        assert inputs["ref"] is released["ref"]
        return thinned, NS(n=3), NS(n=4), 0, {"p": 0.5}

    seen = []

    def form_b_arms(ctx, directions, groups, fraction, held_m1, held_factors):
        seen.append(ctx.inputs)
        cm.gs._FOLD["fold"] = {"finished": {d: NS(record={"finish": f"per-category mean removed, damped by {d:g}, "
                                                                     "zero below log 1.05"}) for d in (0.5, 0.75, 1.0)}}
        return {a: f"pred-{a}" for a in cm.WORLD_ARMS}

    for mod, name, value in [
        (cm.n12, "tries_rows", lambda: {"T59"}), (cm.gc, "n22_passed", lambda p: True),
        (cm.gv, "tree_state", lambda: {}), (cm.gv, "stamp", lambda *a, **k: {"commit": "c", "tree_clean": True}),
        (harness, "load_tasks", lambda: {"T1:test": [t]}), (harness, "STAGE_DAYS", {"r_stage": 9.5, "t_stage": 10.5}),
        (cm, "load_boards", lambda: {"pb": NS(genes=["g0"])}), (cm, "read_stage", lambda *a, **k: None),
        (cm, "Stage", NS(from_anndata=lambda a, g: NS(n=3))),
        (cm, "load_stage", lambda path, lineage: released[path.stem]),
        (harness, "base_role", lambda t_: "ref"), (cm.pdr, "world_route", world_route),
        (harness, "task_board", lambda *a: "board"),
        (cm.ls, "all_anchors", lambda t_, b, tr, rs, seeds: ({s: NS(score=50.0, anchors=None) for s in seeds}, [])),
        (cm.feng_dirs, "feng_directions", lambda genes: ({}, {})), (cm.feng_dirs, "label_groups", lambda labels: {}),
        (cm.feng_dirs, "bracket_fraction", lambda a, b: 0.2), (cm.feng_dirs, "clean_stats_path", lambda: "x"),
        (cm.t1_scope, "held_flat_genes", lambda: []), (cm.hygiene, "load_union", lambda: {"file": "x"}),
        (cm.gc, "held_feng", lambda u, p: []), (harness.TaskView, "of", staticmethod(lambda t_: "view")),
        (harness, "_candidate_board", lambda b: "cboard"), (harness, "_b0_t1", lambda d: (lambda ctx: "b0")),
        (cm.gs, "form_b_arms", form_b_arms),
        (cm.er, "score_arms", lambda b, tr, rs, s, an, preds: {k: {"score": 1.0} for k in preds}),
        (cm.ob, "anchor_record", lambda a: {}), (cm.ls, "prediction_sha256", lambda p: "h"),
        (cm.ls, "sparsity", lambda p: 0.0), (cm.gs, "factor_records", lambda f, p: {}),
    ]:
        monkeypatch.setattr(mod, name, value)
    cm.main(["run", "--world", "ref67", "--out", "t59", "--seeds", "0,1"])
    assert routed == ["ref67"] and len(seen) == 2
    assert all(x is thinned for x in seen)
    rec = json.loads((cm.RESULTS / "t59" / "tasks" / "ref67.json").read_text())
    assert rec["thinning"] == {"p": 0.5} and rec["fraction"] == 0.2 and rec["truth"] == cm.HEART
    assert rec["dampings"] == cm.WORLD_DAMPINGS and rec["dampings_intended"] == cm.WORLD_DAMPINGS




@pytest.mark.parametrize("kw, reason", [
    ({"world_over": {"ref50": None}}, "missing cuts or worlds"),
    ({"cut_over": {WH: None}}, "missing cuts or worlds"),
    ({"cut_over": {CH: {"seeds": range(9)}}}, "seeds"),
    ({"world_over": {"lib3_t67": {"seeds": range(9)}}}, "seeds"),
    ({"cut_over": {CH: {"shift": {"b0": 1e-6}}}}, "do not reproduce"),
    ({"cut_over": {WH: {"shift": {"g_mom_d1": 1e-6}}}}, "do not reproduce"),
    ({"cut_over": {WH: {"shift": {"g_mom_null": 1e-6}}}}, "do not reproduce"),
    ({"cut_over": {CH: {"shift": {"g_mom_perm": -1e-6}}}}, "do not reproduce"),
    ({"world_over": {"ref67": {"shift": {"m1_r_gmom": 1e-6}}}}, "do not reproduce try T54's"),
    ({"world_over": {"lib3": {"shift": {"b0": 1e-6}}}}, "do not reproduce try T50's"),
    ({"drop": (("t49", CH),)}, "do not reproduce"),
    ({"drop": (("t20", WH),)}, "do not reproduce"),
    ({"drop": (("world", "ref50"),)}, "do not reproduce try T54's"),
    ({"cut_over": {CH: {"dampings": {**{a: 0.5 for a in ("g_mom", "g_mom_null", "g_mom_perm")},
                                     "g_mom_d075": 0.75, "g_mom_d1": 1.0, "g_mom_null_d075": 0.75,
                                     "g_mom_null_d1": 0.5, "g_mom_perm_d075": 0.75, "g_mom_perm_d1": 1.0}}}},
     "wrong damping"),
    ({"world_over": {"ref50": {"commit": "d" * 40}}}, "2 commits"),
    ({"cut_over": {WH: {"clean": False}}}, "not clean"),
    ({"cut_over": {CH: {"task": WH}}}, "the file holds task"),
    ({"cut_over": {CH: {"anchor_failures": [{"seed": 3, "reason": "x"}]}}}, "anchors"),
    ({"world_over": {"lib3_t33": {"anchor_failures": [{"seed": 3, "reason": "x"}]}}}, "anchors"),
    ({"world_over": {"lib3_t33": {"thinning": {"p": 0.5}}}}, "the file records p"),
    ({"world_over": {"ref50": {"thinning": {"depth_ratio": 0.52}}}}, "realised depth ratio"),
    ({"world_over": {"ref67": {"thinning": {"any_cell_emptied": True}}}}, "no count"),
    ({"world_over": {"lib3_t67": {"thinning": {"rebuild_largest_difference": 1e-3}}}}, "rebuild check"),
    ({"world_over": {"lib3": {"thinning": {"p": 1.0}}}}, "records a thinning"),
    ({"world_over": {"lib3_t67": {"sha_off": True}}}, "sha256"),
    ({"world_over": {"lib3_t33": {"truth_rows": 99}}}, "truth rows"),
    ({"world_over": {"ref50": {"truth": "ref50"}}}, "holds truth"),
    ({"world_over": {"ref67": {"world": "ref50"}}}, "holds world"),
    ({"cut_over": {CH: {"null_seed": 1}}}, "tries T23's and T21's"),
    ({"cut_over": {WH: {"at_half": {**{f: {"used": "a" * 64, "registered_cache": "a" * 64} for f in ("g_mom",
                                                                                                "g_mom_null")},
                                     "g_mom_perm": {"used": "a" * 64, "registered_cache": "b" * 64}}}}},
     "registered caches' sha256 on ['g_mom_perm']"),
    ({"cut_over": {CH: {"at_half": None}}}, "registered caches' sha256"),
    ({"cut_over": {CH: {"registered_rebuild": None}}}, "rebuilt with the caches cleared"),
    ({"cut_over": {WH: {"registered_rebuild": {"seed": 0, "g_mom_null": {"equal": True},
                                               "g_mom_perm": {"equal": False}}}}},
     "do not equal the arms bit for bit (or the rebuild is not recorded) on ['g_mom_perm']"),
    ({"world_over": {"lib3_t33": {"copy": {"m1_r_gmom_g1": "m1_r_gmom"}}}}, "the damping hook idle"),
    ({"world_over": {"ref50": {"sha_same": ("m1_r_gmom_g075", 4)}}}, "m1_r_gmom_g075 (sha256 on seed 4)"),
    ({"world_over": {"ref67": {"dampings": {"m1_r_gmom": 0.5, "m1_r_gmom_g075": 0.75, "m1_r_gmom_g1": None}}}},
     "wrong damping on ['m1_r_gmom_g1']"),
    ({"world_over": {"ref67": {"fraction": 0.4}}}, "fraction 0.4, not 0.2"),
    ({"world_over": {"lib3": {"fraction": None}}}, "fraction None"),
    ({"cut_over": {WH: {"copy": {"g_mom_null_d1": "g_mom_null"}}}}, "exactly as its twin"),
    ({"cut_over": {CH: {"copy": {"g_mom_perm_d075": "g_mom_d075"}}}}, "exactly as its twin"),
])
def test_the_void_cases(cm, tmp_path, kw, reason):
    r = _read(cm, tmp_path, NEEDS_GAINS, _all(0.5), **kw)[0]["reading"]
    assert not r["valid"] and r["outcome"] == "VOID"
    assert any(reason in x for x in r["validity"]["reasons"]), r["validity"]["reasons"]


def test_a_thinned_reference_scoring_as_lib3_on_every_seed_voids(cm, tmp_path):
    r = _read(cm, tmp_path, NEEDS_GAINS, _all(0.5))[0]["reading"]
    assert r["valid"]
    # the committed try T54 file is made from the same table, so only the idle check can catch it
    res, _ = _read(cm, tmp_path / "idle", NEEDS_GAINS, _all(0.5), world_over={"ref67": {"offset": 0.0}})
    reasons = res["reading"]["validity"]["reasons"]
    assert any("thinning idle" in x for x in reasons)


def test_a_record_off_voids(cm, tmp_path):
    rec = {"permutation_same_at_every_damping": False, "permutation": {"0.5": {"seed": 20261110}},
           "null_categories": {"0.5": ["cm"], "0.75": ["cm"], "1": ["cm"]}}
    r = _read(cm, tmp_path, NEEDS_GAINS, _all(0.5), cut_over={CH: {"records": rec}})[0]["reading"]
    assert any("same at every damping" in x for x in r["validity"]["reasons"])
    rec = {**rec, "permutation_same_at_every_damping": True, "null_categories": {"0.5": ["cm"], "1": ["cm", "en"]}}
    r = _read(cm, tmp_path / "b", NEEDS_GAINS, _all(0.5), cut_over={CH: {"records": rec}})[0]["reading"]
    assert any("twin's categories" in x for x in r["validity"]["reasons"])


def test_the_aggregate_reproduction_and_the_seed_for_seed_one(cm, tmp_path):
    r = _read(cm, tmp_path, NEEDS_GAINS, _all(0.5))[0]["reading"]
    rep = r["validity"]["reproduction"][CH]
    assert rep["holds"] and not rep["seed_for_seed_t21_t23"]
    assert set(rep["largest_difference"]) == {"b0", "g_mom", "g_mom_d075", "g_mom_d1", "g_mom_null", "g_mom_perm"}
    r = _read(cm, tmp_path / "rows", NEEDS_GAINS, _all(0.5), t20_rows=True)[0]["reading"]
    assert r["valid"] and r["validity"]["reproduction"][CH]["seed_for_seed_t21_t23"]
    r = _read(cm, tmp_path / "off", NEEDS_GAINS, _all(0.5), t20_rows=True,
              t20_shift={"g_mom_null": 1e-6})[0]["reading"]
    assert not r["valid"] and any("do not reproduce" in x for x in r["validity"]["reasons"])


def test_a_seeds_ahead_mismatch_against_the_committed_record_voids(cm, tmp_path):
    r = _read(cm, tmp_path, NEEDS_GAINS, _all(0.5))[0]["reading"]
    agg = r["validity"]["reproduction"][CH]["aggregates"]["g_mom_null"]
    assert agg["seeds_ahead"] == agg["seeds_ahead_committed"] == 10
    r = _read(cm, tmp_path / "off", NEEDS_GAINS, _all(0.5), t20_ahead_off=1)[0]["reading"]
    assert not r["valid"] and any("seeds ahead of b0 10, committed 9" in x
                                  for x in r["validity"]["reproduction"][CH]["absent"])


def test_the_reading_refuses_to_overwrite(cm, tmp_path):
    _, out = _read(cm, tmp_path, NEEDS_GAINS, _all(0.5))
    assert (out / "t59.json").exists() and (out / "t59.txt").exists()
    with pytest.raises(SystemExit, match="already holds"):
        cm.read(out)


def test_the_limits_are_the_docstrings_and_say_what_the_brief_says(cm):
    doc = " ".join(cm.__doc__.split())
    assert len(cm.LIMITS) == 17 and all(x in doc for x in cm.LIMITS)
    joined = " ".join(cm.LIMITS)
    for phrase in ("One E8.5 to E9.5 fit serves every seed and both cuts", "Twice the damping is not twice g_mom",
                   "A larger arm redraws more counts", "'Stays positive' has no size threshold",
                   "is not a deeper truth", "One released day, in-sample",
                   "Nothing here says how deep any hidden truth is",
                   "the damping is a size chosen on the cuts that kept g_mom",
                   "'Not shown' is not 'does not need development'", "means the twin or the permuted set gains more",
                   "on the cell-weighted root mean square only", "Statement 1 is not carried into the entry's form",
                   "Statement 2's bars are not symmetric", "Each thinned world is one draw",
                   "share their predictions bit for bit", "one draw of predictions",
                   "D compares unlike shapes at equal root mean square",
                   "The permuted set keeps g_mom's values exactly but not g_mom's size",
                   "grows linearly with the damping", "A partly reads that shift",
                   "is printed and bounds nothing", "the reader checks it before reading D",
                   "'not shown' is the expected outcome even when the effect is real and modest",
                   "printed whatever the label", "cannot fail by itself", "with the caches cleared"):
        assert phrase in joined
    assert "keep" not in cm.EXCEEDS and "adoption result" in cm.LABEL_T59
    assert "needs" not in cm.EXCEEDS and "exceeds" in cm.EXCEEDS and "is ahead" in cm.EXCEEDS
    assert cm.EXCEEDS == ("the larger damping's gain is ahead and exceeds by more than 0.3 on both cuts the gains of "
                          "the no-development twin (matched to g_mom's cell-weighted root mean square at each damping) "
                          "and of the permuted control (g_mom's factor values with the metacells reassigned)")
    assert "at g_mom's size" not in cm.EXCEEDS and f'"{cm.EXCEEDS}"' in doc


# ---------------------------------------------------------------------------------------------------------------------
# The run's refusals


@pytest.mark.parametrize("where", [["--cut", "chamber"], ["--world", "ref50"]])
def test_the_run_needs_its_row(cm, monkeypatch, where):
    monkeypatch.setattr(cm.n12, "tries_rows", lambda: {"T49", "T57"})
    monkeypatch.setattr(cm.gc, "n22_passed", lambda path: True)
    with pytest.raises(SystemExit, match="T59 has no row"):
        cm.main(["run", *where, "--out", "t59"])


@pytest.mark.parametrize("where", [["--cut", "whole"], ["--world", "lib3"]])
def test_the_run_needs_n22(cm, monkeypatch, where):
    monkeypatch.setattr(cm.n12, "tries_rows", lambda: {"T59"})
    with pytest.raises(SystemExit, match="does not hold a passing N22"):
        cm.main(["run", *where, "--out", "t59", "--n22", "n22x"])


@pytest.mark.parametrize("where, name", [(["--cut", "whole"], "T1_E9.5_E10.5.json"),
                                         (["--world", "lib3_t33"], "lib3_t33.json")])
def test_the_run_refuses_an_existing_task_file(cm, monkeypatch, where, name):
    monkeypatch.setattr(cm.n12, "tries_rows", lambda: {"T59"})
    monkeypatch.setattr(cm.gc, "n22_passed", lambda path: True)
    tasks = cm.RESULTS / "t59" / "tasks"
    tasks.mkdir(parents=True)
    (tasks / name).write_text("{}")
    with pytest.raises(SystemExit, match="already holds"):
        cm.main(["run", *where, "--out", "t59"])


@pytest.mark.parametrize("where", [[], ["--cut", "chamber", "--world", "lib3"]])
def test_the_run_takes_exactly_one_cut_or_world(cm, where):
    with pytest.raises(SystemExit, match="exactly one"):
        cm.main(["run", *where, "--out", "t59"])


def test_the_cut_and_world_names_are_the_briefs(cm):
    assert cm.CUTS == {"chamber": "T1/E9.5->E10.5:chamber", "whole": "T1/E9.5->E10.5"}
    with pytest.raises(SystemExit):
        cm.main(["run", "--world", "lib3_t33b", "--out", "t59"])


# ---------------------------------------------------------------------------------------------------------------------
# The second review's fixes


def test_the_permuted_sets_size_and_mean_shift_are_recorded_and_grow_with_the_damping(cm):
    ctx = _ctx()
    ctl = cm.control_factors(ctx)
    rec = cm.control_records(ctl, ctx.derived["t1_held_flat"])
    assert set(rec["mean_shift"]) == set(cm.FAMILIES)
    shift, size = rec["mean_shift"], rec["sizes"]
    for d in ("0.5", "0.75", "1"):
        for c in size["g_mom"][d]:
            # g_mom's mean is removed before the finish, the permuted set's is not: a larger shift and a larger size
            assert shift["g_mom_perm"][d][c] > 2 * shift["g_mom"][d][c]
            assert size["g_mom_perm"][d][c] > size["g_mom"][d][c]
    # on these stages, category cm: the shift is about 0.10 at 0.5 and doubles at 1.0 (linear in the damping)
    assert shift["g_mom_perm"]["0.5"]["cm"] == pytest.approx(0.104, abs=0.002)
    assert shift["g_mom_perm"]["1"]["cm"] == pytest.approx(2 * shift["g_mom_perm"]["0.5"]["cm"], rel=0.02)
    assert size["g_mom_perm"]["1"]["cm"] == pytest.approx(0.2994, abs=1e-3)
    assert size["g_mom"]["1"]["cm"] == pytest.approx(0.2640, abs=1e-3)
    # the record is the root mean square of the per-gene cell-weighted mean over the free genes, by hand
    cf = ctl["g_mom_perm"][1.0].categories["cm"]
    held = set(ctx.derived["t1_held_flat"])
    free = np.array([g not in held for g in ctl["g_mom_perm"][1.0].genes])
    mean = (cf.logf * (cf.sizes / cf.sizes.sum())[:, None]).sum(0)[free]
    assert shift["g_mom_perm"]["1"]["cm"] == pytest.approx(float(np.sqrt((mean ** 2).mean())), rel=1e-9)


def _records_with_shapes():
    def per(v):
        return {d: {"cm": v} for d in ("0.5", "0.75", "1")}
    return {"permutation_same_at_every_damping": True, "permutation": {"0.5": {"seed": 20261110}},
            "null_categories": {"0.5": ["cm"], "0.75": ["cm"], "1": ["cm"]}, "null_match": _match(),
            "factors": {"g_mom": per({"share_not_zeroed": 0.936, "share_at_clip": 0.0026,
                                      "mean_abs_log_not_zeroed": 0.1, "mean_abs_log_all_entries": 0.1}),
                        "g_mom_null": per({"share_not_zeroed": 0.59, "share_at_clip": 0.0256,
                                           "mean_abs_log_not_zeroed": 0.1, "mean_abs_log_all_entries": 0.1})},
            "sizes": {"g_mom": per(0.264), "g_mom_perm": per(0.2994), "g_mom_null": per(0.264)},
            "mean_shift": {"g_mom": per(0.0064), "g_mom_perm": per(0.208), "g_mom_null": per(0.03)}}


def test_beside_statement_one_the_shapes_the_shift_and_d_and_a_are_printed_whatever_the_label(cm, tmp_path):
    # D and A level (not shown): their means, the twin's shape and the permuted set's shift are printed all the same
    gains = {**NEEDS_GAINS, "g_mom_null_d1": 0.9, "g_mom_perm_d1": 0.9}
    rec = _records_with_shapes()
    res, out = _read(cm, tmp_path, gains, _all(0.5), cut_over={CH: {"records": rec}, WH: {"records": rec}})
    r = res["reading"]
    assert r["valid"], r["validity"]["reasons"]
    assert r["statement_1"]["statement"].startswith(cm.NOT_SHOWN)
    b = r["beside_statement_1"][CH]
    assert b["D_mean_at_1"] == pytest.approx(0.1) and b["A_mean_at_1"] == pytest.approx(0.1)
    assert b["twin_shape"]["1"]["cm"] == {"share_not_zeroed": 0.59, "share_at_clip": 0.0256,
                                          "damping_found": 1.1, "g_mom_share_not_zeroed": 0.936,
                                          "g_mom_share_at_clip": 0.0026}
    assert b["permuted_size"]["1"]["cm"] == {"rms": 0.2994, "g_mom_rms": 0.264, "mean_shift": 0.208,
                                             "g_mom_mean_shift": 0.0064}
    text = (out / "t59.txt").read_text()
    assert f"{CH} at damping 1, whatever the label: D mean +0.10, A mean +0.10" in text
    assert "cm not zeroed 0.590 (g_mom 0.936), at clip 0.0256 (g_mom 0.0026), damping found 1.1" in text
    assert "cm rms 0.2994 (g_mom 0.2640), mean shift 0.2080 (g_mom 0.0064)" in text
    assert "D compares unlike shapes" in text and "Q holds the shift besides the assignment" in text
    assert "bounds nothing; check it before reading D" in text


def test_the_text_prints_missing_shape_records_without_failing(cm, tmp_path):
    res, out = _read(cm, tmp_path, NEEDS_GAINS, _all(0.5))
    assert res["reading"]["valid"]
    assert res["reading"]["beside_statement_1"][CH]["twin_shape"]["1"]["cm"]["share_not_zeroed"] is None
    assert "cm not zeroed None (g_mom None)" in (out / "t59.txt").read_text()


# the 8-of-10 edges (a mutant '>' in place of '>=' must fail one of these)


def test_eight_of_ten_ahead_with_a_positive_mean_stays_positive(cm, tmp_path):
    wg = _all(0.5)
    wg["lib3"]["m1_r_gmom_g1"] = [0.5] * 8 + [-0.1] * 2  # mean +0.38, 8 of 10 ahead
    r = _read(cm, tmp_path, NEEDS_GAINS, wg)[0]["reading"]
    assert r["reading_2"]["lib3"]["1"]["seeds_ahead"] == 8
    assert r["statement_2"]["statement"] == cm.STAYS


def test_eight_of_ten_ahead_with_a_mean_at_or_below_zero_does_not_stay_positive(cm, tmp_path):
    wg = _all(0.5)
    wg["lib3"]["m1_r_gmom_g1"] = [0.05] * 8 + [-0.5] * 2  # mean -0.06, 8 of 10 ahead
    wg["ref67"]["m1_r_gmom_g1"] = [0.25] * 8 + [-1.0] * 2  # mean 0 exactly, 8 of 10 ahead
    r = _read(cm, tmp_path, NEEDS_GAINS, wg)[0]["reading"]
    assert r["statement_2"]["statement"].startswith(cm.NOT_SEPARATED)
    assert r["statement_2"]["neither"] == ["lib3", "ref67"]


def test_eight_of_ten_behind_with_a_mean_below_the_bar_turns_and_seven_does_not(cm, tmp_path):
    wg = _all(0.5)
    wg["ref50"]["m1_r_gmom_g1"] = [-0.5] * 8 + [0.1] * 2  # mean -0.38, 8 of 10 behind
    r = _read(cm, tmp_path, NEEDS_GAINS, wg)[0]["reading"]
    assert r["statement_2"]["statement"].startswith(cm.TURNS) and r["statement_2"]["behind"] == ["ref50"]
    wg["ref50"]["m1_r_gmom_g1"] = [-0.6] * 7 + [0.1] * 3  # mean -0.39, 7 of 10 behind
    r = _read(cm, tmp_path / "seven", NEEDS_GAINS, wg)[0]["reading"]
    assert r["statement_2"]["statement"].startswith(cm.NOT_SEPARATED) and r["statement_2"]["behind"] == []


@pytest.mark.parametrize("null_d1, label", [
    ([0.5] * 8 + [1.1] * 2, "ahead"),   # D: 0.5 on 8 seeds, -0.1 on 2: mean +0.38, 8 ahead
    ([0.5] * 7 + [1.1] * 3, "level"),   # 7 ahead, mean +0.32
    ([0.9] * 8 + [1.5] * 2, "level"),   # 8 ahead, mean -0.02
    ([1.5] * 8 + [0.9] * 2, "behind"),  # 8 behind, mean -0.38
    ([1.5] * 7 + [0.9] * 3, "level"),   # 7 behind, mean -0.32
])
def test_the_eight_of_ten_edge_of_statement_ones_d_label(cm, tmp_path, null_d1, label):
    r = _read(cm, tmp_path, {**NEEDS_GAINS, "g_mom_null_d1": null_d1}, _all(0.5))[0]["reading"]
    c = r["reading_1"][CH]["1"]["D"]
    assert c["label"] == label
    assert (r["statement_1"]["statement"] == cm.EXCEEDS) is (label == "ahead")


# the independent rebuild at 0.5, the built dampings and the world's idle hook


def test_the_rebuild_at_half_equals_the_arms_and_puts_the_caches_back(cm, monkeypatch):
    ctx = _ctx()
    preds = cm.cut_arms(ctx, harness._b0_t1(harness.B0_DAMPING))
    g_before = harness._G_MOM_PREP["fold"]["factors"]
    null_before = harness._G_MOM_NULL_PREP["fold"]["factors"]
    fits = []
    real_fit = ts.fit_state_factors
    monkeypatch.setattr(ts, "fit_state_factors", lambda *a, **k: fits.append(1) or real_fit(*a, **k))
    rec = cm.registered_rebuild(ctx, preds)
    assert len(fits) == 3  # the caches were cleared: g_mom's fit and the twin's two directions made again
    assert all(rec[a]["equal"] for a in cm.T20_ARMS)
    assert harness._G_MOM_PREP["fold"]["factors"] is g_before
    assert harness._G_MOM_NULL_PREP["fold"]["factors"] is null_before
    swapped = {**preds, "g_mom_perm": preds["g_mom_perm_d1"]}
    rec = cm.registered_rebuild(ctx, swapped)
    assert rec["g_mom_null"]["equal"] and not rec["g_mom_perm"]["equal"]


def test_the_dampings_are_read_from_the_built_factors(cm):
    ctx = _ctx()
    ctl = cm.control_factors(ctx)
    assert cm.built_cut_dampings(ctl) == cm.CUT_ARMS
    wrong = {**ctl, "g_mom_perm": {**ctl["g_mom_perm"], 1.0: ctl["g_mom_perm"][0.5]},
             "g_mom_null": {**ctl["g_mom_null"], 0.75: ctl["g_mom_null"][1.0]}}
    got = cm.built_cut_dampings(wrong)
    assert got["g_mom_perm_d1"] == 0.5 and got["g_mom_null_d075"] == 1.0
    assert cm._dampings_off(got, cm.CUT_ARMS) == ["g_mom_null_d075", "g_mom_perm_d1"]
    assert cm.built_world_dampings(cm.gs.fold_factors(ctx)["finished"]) == cm.WORLD_DAMPINGS
    assert cm.built_world_dampings({0.5: ctl["g_mom"][0.5]})["m1_r_gmom_g1"] is None
    empty = ts.StateFactors(genes=[], metacells_on="later", categories={}, skipped={})
    assert cm.finish_damping(empty) is None and cm.twin_damping(empty) is None


def test_the_cut_run_writes_the_built_dampings_and_the_rebuild(cm, monkeypatch):
    """cmd_run_cut with the loaders and the scorer stubbed, on the synthetic stages: the task file's dampings come from
    the built factors and the rebuild check is recorded once, on the first seed, and passes."""
    from types import SimpleNamespace as NS

    base = _ctx()
    ref, start = base.inputs["ref"], base.inputs["start"]
    t = NS(key=cm.gc.TASK_CHAMBER, panel_board="pb", reference="r_stage", target="t_stage", target_filter=None,
           inputs={"ref": "ref.h5ad", "start": "start.h5ad"})
    rebuilds = []
    real_rebuild = cm.registered_rebuild
    for mod, name, value in [
        (cm.n12, "tries_rows", lambda: {"T59"}), (cm.gc, "n22_passed", lambda p: True),
        (cm.gv, "tree_state", lambda: {}), (cm.gv, "stamp", lambda *a, **k: {"commit": "c", "tree_clean": True}),
        (harness, "load_tasks", lambda: {"T1:test": [t]}),
        (cm, "load_boards", lambda: {"pb": NS(genes=["g0"])}), (cm, "read_stage", lambda *a, **k: None),
        (cm, "Stage", NS(from_anndata=lambda a, g: NS(n=3))),
        (cm, "load_stage", lambda path, lineage: {"ref": ref, "start": start}[path.stem]),
        (harness, "task_board", lambda *a: "board"),
        (cm.ls, "all_anchors", lambda t_, b, tr, rs, seeds: ({s: NS(score=50.0, anchors=None) for s in seeds}, [])),
        (cm.t1_scope, "held_flat_genes", lambda: base.derived["t1_held_flat"]),
        (harness.TaskView, "of", staticmethod(lambda t_: base.task)),
        (harness, "_candidate_board", lambda b: base.board),
        (cm.er, "score_arms", lambda b, tr, rs, s, an, preds: {k: {"score": 1.0} for k in preds}),
        (cm, "registered_rebuild", lambda ctx, preds: rebuilds.append(ctx.seed) or real_rebuild(ctx, preds)),
    ]:
        monkeypatch.setattr(mod, name, value)
    cm.main(["run", "--cut", "chamber", "--out", "t59", "--seeds", "0,1"])
    rec = json.loads((cm.RESULTS / "t59" / "tasks" / cm.CUT_FILES[cm.gc.TASK_CHAMBER]).read_text())
    assert rebuilds == [harness.BUILD_SEED_OFFSET]
    assert rec["dampings"] == cm.CUT_ARMS and rec["dampings_intended"] == cm.CUT_ARMS
    assert rec["registered_rebuild"]["seed"] == 0
    assert all(rec["registered_rebuild"][a]["equal"] for a in cm.T20_ARMS)
    assert set(rec["records"]["mean_shift"]) == set(cm.FAMILIES)
    assert cm.record_reasons("x", rec) == []


def test_a_world_whose_larger_arm_is_m1_r_gmom_voids(cm):
    ok = _world(cm, "lib3", _all(0.5))
    assert cm.world_idle_reasons("lib3", ok) == []
    idle = _world(cm, "lib3", _all(0.5), copy={"m1_r_gmom_g075": "m1_r_gmom"})
    assert "m1_r_gmom_g075" in cm.world_idle_reasons("lib3", idle)[0]
    one = _world(cm, "lib3", _all(0.5), sha_same=("m1_r_gmom_g1", 0))
    assert "m1_r_gmom_g1 (sha256 on seed 0)" in cm.world_idle_reasons("lib3", one)[0]
