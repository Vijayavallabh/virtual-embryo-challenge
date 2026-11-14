"""experiments/t1_m1/eb_read.py (try T42 of experiments/harness/TRIES.md) on synthetic data: the shrunk directions, the
arms' twins and the reading. No statistics file, released file or score is read."""

from __future__ import annotations

import importlib.util
import json
import sys

import numpy as np
import pytest
from test_g_mom_verify import _small_inputs
from test_t1_state_transport import _gmom_ctx

from vec import harness, t1_scope
from vec.boards import REPO_ROOT

SCRIPT = REPO_ROOT / "experiments" / "t1_m1" / "eb_read.py"
GROUPS = {"V-CM": "CM", "Endocardium": "Endocardium", "Blood": "Blood", "Foregut": "other:Foregut"}


@pytest.fixture
def er(monkeypatch, tmp_path):
    spec = importlib.util.spec_from_file_location("eb_read_module", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    monkeypatch.setattr(mod, "RESULTS", tmp_path / "results")
    return mod


def _halves(d, noise, seed, agree=False):
    rng = np.random.default_rng(seed)
    if agree:
        return {"d_a": d.copy(), "d_b": d.copy(), "mean_cp": np.linspace(0.1, 40.0, d.size)}
    e = rng.normal(0, noise, d.size)
    return {"d_a": d + e, "d_b": d - e, "mean_cp": rng.uniform(0.1, 40.0, d.size)}


def _synthetic(n=300, agree=()):
    rng = np.random.default_rng(11)
    dirs = {g: rng.normal(0, 0.5, n) for g in ("CM", "Endocardium")}
    dirs["CM"][:7] = np.nan
    halves = {g: _halves(d, 0.25, i, agree=g in agree) for i, (g, d) in enumerate(dirs.items())}
    return dirs, halves


# ---------------------------------------------------------------------------------------------------------------------
# The shrunk directions


def test_the_decile_form_is_eb_shrinks_output_and_the_scalar_form_is_scale_times_d(er):
    import feng_dirs
    import n48_m1_direction_shrink as n48

    dirs, halves = _synthetic()
    eb, js, record = er.shrunk_directions(dirs, halves)
    assert set(eb) == set(js) == set(record) == set(dirs)
    for g, d in dirs.items():
        h = halves[g]
        want, factors = feng_dirs.eb_shrink(d, h["d_a"], h["d_b"], h["mean_cp"], harness.LP_BINS)
        scale = n48.js_scale(d, h["d_a"], h["d_b"])["scale"]
        assert np.array_equal(eb[g], want, equal_nan=True)
        assert np.array_equal(js[g], scale * d, equal_nan=True)
        assert record[g] == {"eb_factors": [float(f) for f in factors], "js_scale": float(scale)}
        assert 0 < scale < 1 and min(factors) < 1, "the synthetic halves disagree: something is shrunk"
        assert np.array_equal(np.isnan(eb[g]), np.isnan(d)) and np.array_equal(np.isnan(js[g]), np.isnan(d))


def test_a_group_whose_halves_agree_keeps_its_direction_bit_for_bit(er):
    dirs, halves = _synthetic(agree=("CM",))
    eb, js, record = er.shrunk_directions(dirs, halves)
    assert eb["CM"].tobytes() == dirs["CM"].tobytes() and js["CM"].tobytes() == dirs["CM"].tobytes()
    assert record["CM"] == {"eb_factors": [1.0] * harness.LP_BINS, "js_scale": 1.0}
    assert not np.array_equal(eb["Endocardium"], dirs["Endocardium"])


def test_a_direction_without_halves_is_refused(er):
    dirs, halves = _synthetic()
    del halves["Endocardium"]
    with pytest.raises(ValueError, match="cannot be shrunk"):
        er.shrunk_directions(dirs, halves)


def _arms_world():
    harness._G_MOM_PREP.clear()
    ref, start, _ = _small_inputs()
    ctx = _gmom_ctx(ref, start)
    groups = {lab: GROUPS.get(lab, "x") for lab in ref.obs[t1_scope.LABEL_KEY].astype(str).unique()}
    rng = np.random.default_rng(5)
    dirs = {g: rng.normal(0, 0.4, ref.n_vars) for g in ("CM", "Endocardium")}
    return ctx, groups, dirs, list(ctx.derived["t1_held_flat"])


def test_the_twins_differ_in_the_direction_only_and_held_genes_stay_held(er, monkeypatch):
    ctx, groups, dirs, held = _arms_world()
    halves = {g: _halves(d, 0.25, i) for i, (g, d) in enumerate(dirs.items())}
    eb, js, _ = er.shrunk_directions(dirs, halves)
    calls = []
    real = t1_scope.external_shift

    def spy(target_like, directions, fraction, group_of, held_flat=(), **kw):
        calls.append({"directions": directions, "fraction": fraction, "groups": dict(group_of),
                      "held": list(held_flat), "seed": kw.get("seed"), "cells": list(target_like.obs_names)})
        return real(target_like, directions, fraction, group_of, held_flat, **kw)

    monkeypatch.setattr(t1_scope, "external_shift", spy)
    arms = er.all_arms(ctx, dirs, eb, js, groups, 0.2, held, held)
    assert set(arms) == set(er.ARMS) - {"b0"}
    assert len(calls) == 6, "one direction step per M1 build"
    want = [dirs, dirs, eb, eb, js, js]
    for call, d in zip(calls, want, strict=True):
        assert set(call["directions"]) == set(d)
        assert all(np.array_equal(call["directions"][g], d[g], equal_nan=True) for g in d)
    for key in ("fraction", "groups", "held", "seed", "cells"):
        assert all(c[key] == calls[0][key] for c in calls), f"the twins differ in {key}"
    assert calls[0]["held"] and set(held) <= set(calls[0]["held"]), "held genes stay held"
    for base in ("m1_r", "m1_r_gmom"):
        assert arms[base].X.shape == arms[base + "_eb"].X.shape == arms[base + "_js"].X.shape
        assert not np.array_equal(arms[base].X, arms[base + "_eb"].X)
        assert not np.array_equal(arms[base].X, arms[base + "_js"].X)


def test_with_agreeing_halves_the_shrunk_arms_equal_their_twins_bit_for_bit(er):
    ctx, groups, dirs, held = _arms_world()
    halves = {g: _halves(d, 0.0, i, agree=True) for i, (g, d) in enumerate(dirs.items())}
    eb, js, _ = er.shrunk_directions(dirs, halves)
    arms = er.all_arms(ctx, dirs, eb, js, groups, 0.2, held, held)
    for base in ("m1_r", "m1_r_gmom"):
        for suffix in ("_eb", "_js"):
            assert arms[base + suffix].X.tobytes() == arms[base].X.tobytes()


# ---------------------------------------------------------------------------------------------------------------------
# The reading


N48_GROUPS = {g: {"eb_factors": [0.5 + 0.05 * i for i in range(10)], "js": {"scale": 0.9}}
              for g in ("CM", "Endocardium", "Endothelium", "Epicardium")}


def _arm(score):
    return {"score": score, "groups": {"distribution": score / 100},
            "metrics": {"mmd_u": 0.1, "variogram": 0.01},
            "skills": {"de_score": score / 200, "de_direction": 0.4, "mmd_u": 0.3, "variogram": 0.2}}


def _task(gains, b0=50.0, seeds=range(10), commit="c" * 40, clean=True, factor_shift=0.0):
    base = {"m1_r": 52.0, "m1_r_gmom": 53.0}
    rows = []
    for s in seeds:
        arms = {"b0": _arm(b0), **{a: _arm(v) for a, v in base.items()}}
        for a in ("m1_r_eb", "m1_r_gmom_eb", "m1_r_js", "m1_r_gmom_js"):
            arms[a] = _arm(base["m1_r_gmom" if "gmom" in a else "m1_r"] + gains.get(a, 0.0))
        rows.append({"seed": s, "floor": 50.0, "arms": arms})
    shrink = {g: {"eb_factors": [f + factor_shift for f in v["eb_factors"]], "js_scale": v["js"]["scale"]}
              for g, v in N48_GROUPS.items()}
    return {"label": "synthetic", "seeds": list(seeds), "rows": rows, "shrinkage": shrink,
            "stamp": {"commit": commit, "tree_clean": clean, "utc": "x", "inputs": {}}}


def _write(er, tmp_path, cuts):
    """cuts: cut key -> task record. T26's committed files and N48's file, synthetic, beside them."""
    out = tmp_path / "out"
    (out / "tasks").mkdir(parents=True)
    for k, rec in cuts.items():
        (out / "tasks" / er.CUT_FILES[k]).write_text(json.dumps(rec))
    t26 = tmp_path / "t26"
    t26.mkdir()
    for f in er.CUT_FILES.values():
        rows = [{"seed": s, "arms": {"b0": {"score": 50.0}, "m1_r": {"score": 52.0}, "m1_r_gmom": {"score": 53.0}}}
                for s in range(10)]
        (t26 / f).write_text(json.dumps({"rows": rows}))
    n48 = tmp_path / "n48.json"
    n48.write_text(json.dumps({"groups": N48_GROUPS}))
    return out, t26, n48


def _read(er, tmp_path, gains=None, **over):
    gains = gains or {}
    cuts = {k: _task(gains) for k in er.CUT_FILES}
    for k, kw in over.items():
        if kw is None:
            del cuts[getattr(er.gc, k)]
        else:
            cuts[getattr(er.gc, k)] = _task(gains, **kw)
    out, t26, n48 = _write(er, tmp_path, cuts)
    return er.read(out, t26, n48)["reading"]


@pytest.mark.parametrize("gains, verdict", [
    ({"m1_r_eb": 1.0, "m1_r_gmom_eb": 1.0}, "DECILE"),  # (a), (b), (e)
    ({"m1_r_eb": 1.0, "m1_r_gmom_eb": 1.0, "m1_r_js": 1.0, "m1_r_gmom_js": 1.0}, "SCALAR"),  # not (e), (c), (d)
    ({"m1_r_eb": 1.0, "m1_r_gmom_eb": 1.0, "m1_r_gmom_js": 0.8}, "DECILE"),  # (a), (b), not (e), not (c)
    ({"m1_r_js": 1.0, "m1_r_gmom_js": 1.0}, "SCALAR"),  # not (a), (c), (d)
    ({"m1_r_eb": 1.0, "m1_r_js": 1.0}, "NOTHING"),  # (a) and (c) only
    ({}, "NOTHING"),
])
def test_each_verdict_is_reached_by_a_table_built_for_it(er, tmp_path, gains, verdict):
    r = _read(er, tmp_path, gains)
    assert r["valid"] and r["outcome"] == getattr(er, verdict)
    comp = r["comparisons"]["b"]["cuts"][er.gc.TASK_CHAMBER]
    want = gains.get("m1_r_gmom_eb", 0.0)
    assert comp["mean"] == pytest.approx(want) and comp["ahead"] == (10 if want > 0 else 0)
    assert set(comp["skills"]) == {"de_score", "de_direction", "mmd_u", "variogram"}
    assert comp["skills"]["de_score"] == pytest.approx(want / 200)
    assert r["mean_scores"][er.gc.TASK_WHOLE]["m1_r_gmom_eb"] == pytest.approx(53.0 + want)


def test_a_gain_of_exactly_point_three_does_not_pass(er, tmp_path):
    r = _read(er, tmp_path, {"m1_r_eb": 0.3, "m1_r_gmom_eb": 0.3})
    assert r["valid"] and not r["comparisons"]["a"]["passes"] and not r["comparisons"]["b"]["passes"]
    assert r["outcome"] == er.NOTHING
    r = _read(er, tmp_path / "more", {"m1_r_eb": 0.31, "m1_r_gmom_eb": 0.31})
    assert r["comparisons"]["a"]["passes"] and r["outcome"] == er.DECILE


def test_a_loss_beyond_point_three_on_one_cut_vetoes(er, tmp_path):
    gains = {"m1_r_eb": 1.0, "m1_r_gmom_eb": 1.0}
    cuts = {er.gc.TASK_CHAMBER: _task(gains), er.gc.TASK_WHOLE: _task({"m1_r_eb": -0.4, "m1_r_gmom_eb": 1.0})}
    out, t26, n48 = _write(er, tmp_path, cuts)
    r = er.read(out, t26, n48)["reading"]
    assert not r["comparisons"]["a"]["passes"] and r["comparisons"]["b"]["passes"]
    assert r["outcome"] == er.NOTHING


@pytest.mark.parametrize("over, reason", [
    ({"TASK_WHOLE": None}, "missing cuts"),
    ({"TASK_CHAMBER": {"seeds": range(9)}}, "seeds"),
    ({"TASK_WHOLE": {"b0": 50.03}}, "b0 does not reproduce"),
    ({"TASK_CHAMBER": {"factor_shift": 1e-9}}, "factors differ"),
    ({"TASK_WHOLE": {"clean": False}}, "not clean"),
    ({"TASK_WHOLE": {"commit": "d" * 40}}, "2 commits"),
])
def test_the_void_cases(er, tmp_path, over, reason):
    r = _read(er, tmp_path, {"m1_r_eb": 1.0, "m1_r_gmom_eb": 1.0}, **over)
    assert not r["valid"] and r["outcome"] == "VOID"
    assert any(reason in x for x in r["validity"]["reasons"]), r["validity"]["reasons"]


def test_b0_within_the_tolerance_and_tiny_factor_noise_stay_valid(er, tmp_path):
    r = _read(er, tmp_path, {"m1_r_eb": 1.0, "m1_r_gmom_eb": 1.0},
              TASK_WHOLE={"b0": 50.01, "factor_shift": 1e-14})
    assert r["valid"] and r["outcome"] == er.DECILE
    assert r["validity"]["reproduction"][er.gc.TASK_WHOLE]["largest_difference"]["b0"] == pytest.approx(0.01)


def test_the_reading_refuses_to_overwrite_and_the_run_needs_n22(er, tmp_path):
    out, t26, n48 = _write(er, tmp_path, {k: _task({}) for k in er.CUT_FILES})
    er.read(out, t26, n48)
    assert (out / "t42.json").exists() and (out / "t42.txt").exists()
    with pytest.raises(SystemExit, match="already holds"):
        er.read(out, t26, n48)
    with pytest.raises(SystemExit, match="refused"):
        er.main(["run", "--task", er.gc.TASK_CHAMBER, "--out", "t42", "--n22", "n22x"])
