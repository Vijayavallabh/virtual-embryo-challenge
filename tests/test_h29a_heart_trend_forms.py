"""Try H29a (experiments/harness/h29a_heart_trend_forms.py) and its launcher (experiments/harness/run_h29a.sh):
synthetic arrays and hand-made tables only, and committed result files of earlier rows read for arithmetic. Nothing
here loads a competition file or scores a real cell."""

from __future__ import annotations

import ast
import importlib.util
import inspect
import json
import os
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from vec import harness
from vec.adoption import SEEDS

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "experiments" / "harness" / "h29a_heart_trend_forms.py"
LAUNCHER = REPO_ROOT / "experiments" / "harness" / "run_h29a.sh"
RESULTS = REPO_ROOT / "experiments" / "harness" / "results"


@pytest.fixture(scope="module")
def h29():
    spec = importlib.util.spec_from_file_location("h29a_heart_trend_forms_module", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def synth():
    spec = importlib.util.spec_from_file_location("b0_trend_tests_for_h29a", REPO_ROOT / "tests" / "test_b0_trend.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    inputs = {"ref": mod._stage(1600, seed=2, shift=0.0, tiny=60),
              "start": mod._stage(1100, seed=3, shift=0.8, tiny=10),
              "earlier": mod._stage(900, seed=5, shift=1.6, tiny=10)}
    for a in inputs.values():
        a.obs["core"] = np.asarray(a.obs["merged"]).astype(str) == "CM"
    return mod, inputs


def _pools(h29, inputs):
    whole = inputs["ref"]
    PC = np.asarray(whole.obsm["spatial_3D"], dtype=float)
    return h29.form_pools(whole, PC, np.asarray(whole.obs["core"], bool))


def _size(h29, inputs):
    sp = h29.h27.sp
    geoms = [sp.geometry_from_arrays(n, t, np.asarray(inputs[r].obsm["spatial_3D"], float),
                                     np.asarray(inputs[r].obs["core"], bool))
             for r, n, t in (("earlier", "E8.5", 8.5), ("start", "E8.75", 8.75), ("ref", "E9.5", 9.5))]
    return h29.h27.analogue_size_result(geoms, 10.5)


# ------------------------------------------------------------------------------------------------ row H28 unchanged
def test_row_h28s_committed_reading_is_reproduced_byte_for_byte(h29):
    """The H28 path after the parameterisation: its read on its committed parts and row H27's committed parts gives
    reading.json and reading.txt exactly as committed (arithmetic on committed files; no new reading)."""
    h28 = h29.h28
    out = RESULTS / "h28_heart_trend_strength"
    committed = json.loads((out / "reading.json").read_text())
    res = h28.reading(h28.load_parts(out / "parts"), h28.load_h27(), committed["commit"])
    assert json.dumps(res, indent=1, default=float) == (out / "reading.json").read_text()
    assert h28.text(res) == (out / "reading.txt").read_text()


def test_row_h28s_build_arms_defaults_are_its_own(h29, synth):
    """The new keywords default to row H28's form, strengths and trend_s0, and the default path builds row H28's arms
    and checks in row H28's order."""
    h28 = h29.h28
    sig = inspect.signature(h28.build_arms)
    assert sig.parameters["form"].default == "core_iso_half" == h28.FORM
    assert sig.parameters["strengths"].default == (1.0, 0.5, 0.25) == h28.STRENGTHS
    assert sig.parameters["zero"].default is True
    mod, inputs = synth
    pool = inputs["ref"]
    PC = np.asarray(pool.obsm["spatial_3D"], dtype=float)
    keep = h28.h26.pool_masks(PC, np.asarray(pool.obs["core"], bool))[h28.FORM]
    ins = {**inputs, "ref": pool[keep].copy()}
    out, extra = h28.build_arms(mod._ctx(ins, seed=harness.BUILD_SEED_OFFSET + 3), _size(h29, inputs),
                                whole_ref=inputs["ref"])
    assert list(out) == list(h28.BUILT) == ["control", "trend_s0", "trend_s1", "placebo_s1", "trend_s0.5",
                                            "placebo_s0.5", "trend_s0.25", "placebo_s0.25"]
    assert list(extra["checks"]) == ["s1_equals_harness_b0_trend", "s0_counts_equal_drawn",
                                     "s0_x_max_abs_diff_control", "same_draw", "pre_size_coords_equal",
                                     "post_size_coords_equal", "n_cells", "whole_pool"]
    assert all(extra["records"][k]["size"]["hypothesis"] == "core" for k in out)


# ------------------------------------------------------------------------------------------------ synthetic builds
def test_each_forms_arms_on_its_own_pool_and_size(h29, synth):
    mod, inputs = synth
    pools, info = _pools(h29, inputs)
    assert pools[h29.FULL] is inputs["ref"] and info[h29.FULL]["pool_cells"] == inputs["ref"].n_obs
    assert info[h29.HALF]["pool_cells"] == pools[h29.HALF].n_obs < inputs["ref"].n_obs
    size = _size(h29, inputs)
    for form, hyp in ((h29.FULL, "H1"), (h29.HALF, "H2")):
        assert h29.h27.HYP[form] == hyp
        ctx = mod._ctx({**inputs, "ref": pools[form]}, seed=harness.BUILD_SEED_OFFSET + 5)
        out, extra = h29.build_form(form, ctx, size, whole_ref=None if form == h29.FULL else inputs["ref"])
        assert list(out) == [h29.name(form, a) for a in h29.LOCAL]
        c = extra["checks"]
        for k in ("s1_equals_harness_b0_trend", "same_draw", "pre_size_coords_equal", "post_size_coords_equal"):
            assert c[k] is True, (form, k)
        assert "s0_counts_equal_drawn" not in c and ("whole_pool" in c) == (form == h29.HALF)
        assert len(set(c["n_cells"].values())) == 1 and c["n_cells"]["control"] == mod.BOARD.max_cells // 2
        for k, st in out.items():
            np.testing.assert_array_equal(st.coords, out[h29.name(form, "control")].coords)
            rec = extra["records"][k]
            assert rec["size"]["hypothesis"] == hyp and abs(rec["size"]["log_rms_after_vs_target"]) <= 1e-9
        tr = extra["records"][h29.name(form, "trend_s1")]["trend"]
        assert tr["strength"] == 1.0 and not tr["placebo"] and tr["ref_cells"] == pools[form].n_obs
        pl = extra["records"][h29.name(form, "placebo_s1")]["trend"]
        assert pl["placebo"] and not any(pl["fixed_points"].values()) and pl["shifted"] == tr["shifted"]
        assert not np.array_equal(out[h29.name(form, "trend_s1")].X, out[h29.name(form, "control")].X)
        assert not np.array_equal(out[h29.name(form, "trend_s1")].X, out[h29.name(form, "placebo_s1")].X)
        # the control is row H27's f arm on the same context: h27.build_forms's f sizing, bit for bit
        ref = h29.h27.build_forms(form, mod._ctx({**inputs, "ref": pools[form]}, seed=harness.BUILD_SEED_OFFSET + 5),
                                  {"even": 1.0}, size, {form: 1.0})["f"][0]
        np.testing.assert_array_equal(out[h29.name(form, "control")].X, ref.X)
        np.testing.assert_array_equal(out[h29.name(form, "control")].coords, ref.coords)


# ------------------------------------------------------------------------------------------------ hand-made tables
METRICS = ("de_score", "de_direction", "mmd_u", "variogram", "d2_shape", "occupancy_dice", "scale_log_ratio",
           "neighborhood_mmd")
POOL = {"full_iso_half": 1600, "half_iso_half": 900}


def _info(h29):
    t = {q: {"cells": 1000 + i, "core_share": 0.4, "xy_scope": 1.0} for i, q in enumerate(h29.TRUTHS)}
    return t, {f: {"pool_cells": POOL[f], "core_share": 0.5, "xy_scope": 0.2} for f in h29.FORMS}


def _rec(h29, key, nm):
    form, local = nm.split(":")
    if local == "control":
        return {"size": {}, "trend": {"strength": None}}
    return {"size": {}, "trend": {"strength": 1.0, "placebo": local.startswith("placebo"),
                                  "budget_days": h29.PRACTICE_BUDGET.get(key, 0.5), "ref_cells": POOL[form],
                                  "shifted": ["CM", "Endocardium"], "fixed_points": {"CM": 0, "Endocardium": 0}}}


def _part(h29, key, seeds, score=None, de=None, commit="c" * 40, clean=True, checks=None):
    """score: arm name -> per-seed array of board points added to 50 in every truth (or a dict world -> array);
    de: arm name -> skill added to de_score and de_direction (or a dict world -> value)."""
    score, de = score or {}, de or {}
    t_info, p_info = _info(h29)
    rows = []
    for s in seeds:
        sc, sk, gr = {}, {}, {}
        for q in h29.TRUTHS:
            w = "own" if q == h29.OWN else "CLEAN"
            sc[q], sk[q], gr[q] = {}, {}, {}
            for nm in h29.ALL:
                v = score.get(nm, 0.0)
                v = v.get(w, 0.0) if isinstance(v, dict) else v
                sc[q][nm] = 50.0 + float(np.broadcast_to(v, (10,))[s])
                d = de.get(nm, 0.0)
                d = d.get(w, 0.0) if isinstance(d, dict) else d
                sk[q][nm] = {m: 0.5 + (d if m in h29.DE_TERMS else 0.0) for m in METRICS}
                gr[q][nm] = {g: 0.5 for g in ("expression_change", "cell_state", "shape_growth",
                                               "local_organisation")}
        ck = {}
        for f in h29.FORMS:
            ck[f] = {"s1_equals_harness_b0_trend": True, "same_draw": True, "pre_size_coords_equal": True,
                     "post_size_coords_equal": True, "n_cells": {a: 350 for a in h29.LOCAL}}
            ck[f].update((checks or {}).get(f, {}))
        rows.append({"seed": s, "n_cells": {n: 350 for n in h29.BUILT}, "checks": ck, "scores": sc, "skills": sk,
                     "groups": gr, "records": {n: _rec(h29, key, n) for n in h29.BUILT}})
    return {"axis": h29.AXIS, "task": key, "rows": rows, "truth_info": t_info, "pred_info": p_info,
            "forecast_targets": {f: 100.0 + i for i, f in enumerate(("full_iso_half", "half_iso_half",
                                                                     "core_iso_half"))},
            "v7": {"ok": True, "log_diff": 0.0}, "provenance": {"commit": commit, "tree_clean": clean}}


def _parts(h29, **kw):
    return {k: [_part(h29, k, range(0, 5), **kw), _part(h29, k, range(5, 10), **kw)] for k in h29.TASK_KEYS}


def _h27(h29, off=None):
    """Row H27's committed shape: every form's f arm at 50 (+ off[form])."""
    off = off or {}
    t_info, p_info = _info(h29)
    rows = [{"seed": s, "scores": {q: {h29.h27_arm(f): 50.0 + off.get(f, 0.0) for f in h29.FORMS}
                                   for q in h29.TRUTHS}} for s in SEEDS]
    return {k: [{"rows": rows, "truth_info": t_info, "pred_info": p_info,
                 "forecast_targets": {f: 100.0 + i for i, f in enumerate(("full_iso_half", "half_iso_half",
                                                                          "core_iso_half"))}}]
            for k in h29.TASK_KEYS}


def _read(h29, parts=None, h27p=None, commit="c" * 40):
    return h29.reading(parts if parts is not None else _parts(h29), h27p if h27p is not None else _h27(h29), commit)


def _good(h29, form, gain=0.5, placebo_gain=0.0, de=0.01):
    return ({h29.name(form, "trend_s1"): np.full(10, gain), h29.name(form, "placebo_s1"): np.full(10, placebo_gain)},
            {h29.name(form, "trend_s1"): de})


def _both(h29, *specs):
    score, de = {}, {}
    for s, d in specs:
        score.update(s)
        de.update(d)
    return _parts(h29, score=score, de=de)


def test_a_clean_null_run_qualifies_neither(h29):
    res = _read(h29)
    assert not res["void"] and {c: v["outcome"] for c, v in res["candidates"].items()} == {
        "H-T1": "NOT QUALIFIED", "H-T1B": "NOT QUALIFIED"}
    assert res["guards"]["V4"][h29.FULL][h29.TASK_KEYS[0]] == {"mean": 0.0, "seed": 0.0}
    txt = h29.text(res)
    assert txt.splitlines()[0].startswith(f"Try H29a (commit {'c' * 40}). A diagnostic on one truth")
    assert h29.NOT_BUILT_LINE in txt and "OUTCOME H-T1: NOT QUALIFIED" in txt and "OUTCOME H-T1B: NOT QUALIFIED" in txt
    assert "row H27's committed full_iso_half:f" in txt and "row H27's committed half_iso_half:f" in txt
    assert all(ln.startswith(h29.TAG) for ln in txt.split("TABLES")[1].splitlines()[1:])
    assert len(txt.split("TABLES")[1].splitlines()[1:]) == 2 * 2 * 3 * 3  # forms x worlds x units x pairs


def test_h_t1_qualifies_on_s1_s2_s3_on_the_whole_pool(h29):
    res = _read(h29, _both(h29, _good(h29, h29.FULL)))
    assert res["candidates"]["H-T1"]["outcome"] == "QUALIFIED"
    assert res["candidates"]["H-T1B"]["outcome"] == "NOT QUALIFIED"  # its own S1 fails
    assert res["follows"]["H-T1"] == h29.FOLLOWS["H-T1", "QUALIFIED"]
    assert "says yes to it by name" in res["follows"]["H-T1"] and "never a point" in res["follows"]["H-T1"]
    assert res["statements"][h29.FULL] == {"S1": True, "S2": True, "S3": True, "qualified": True}


def test_h_t1b_needs_its_own_s1_to_s3_and_h_t1s_s2(h29):
    # its own three hold and H-T1's S2 holds (H-T1 itself failing S1): H-T1B qualifies, H-T1 does not
    full = ({h29.name(h29.FULL, "trend_s1"): np.full(10, 0.2), h29.name(h29.FULL, "placebo_s1"): np.full(10, -0.5)},
            {h29.name(h29.FULL, "trend_s1"): 0.01})
    res = _read(h29, _both(h29, full, _good(h29, h29.HALF)))
    assert res["statements"][h29.FULL] == {"S1": False, "S2": True, "S3": True, "qualified": False}
    assert res["candidates"]["H-T1"]["outcome"] == "NOT QUALIFIED"
    assert res["candidates"]["H-T1"]["why"] == "fails S1"
    assert res["candidates"]["H-T1B"]["outcome"] == "QUALIFIED"
    assert "K1's scope for this one file" in res["follows"]["H-T1B"]
    # its own three hold but H-T1's S2 fails on the whole pool (H28's confound): H-T1B does not qualify
    full = ({h29.name(h29.FULL, "trend_s1"): np.full(10, 0.5), h29.name(h29.FULL, "placebo_s1"): np.full(10, 0.4)},
            {h29.name(h29.FULL, "trend_s1"): 0.01})
    res = _read(h29, _both(h29, full, _good(h29, h29.HALF)))
    assert res["statements"][h29.HALF]["qualified"] and not res["statements"][h29.FULL]["S2"]
    assert res["candidates"]["H-T1B"] == {"outcome": "NOT QUALIFIED", "why": "fails H-T1's S2 on the whole pool"}
    assert res["follows"]["H-T1B"].startswith("H-T1B NOT QUALIFIED under row H29a (fails H-T1's S2 on the whole")
    # both qualify
    res = _read(h29, _both(h29, _good(h29, h29.FULL), _good(h29, h29.HALF)))
    assert [v["outcome"] for v in res["candidates"].values()] == ["QUALIFIED", "QUALIFIED"]


def test_the_bars_at_their_edges(h29):
    eight = np.array([1.0] * 8 + [-0.5] * 2)  # mean 0.7, 8 of 10
    seven = np.array([1.0] * 7 + [-0.1] * 3)  # mean 0.67, 7 of 10
    f = h29.name(h29.FULL, "trend_s1")
    for arm, ok in ((eight, True), (seven, False), (np.full(10, 0.3), False), (np.full(10, 0.3001), True)):
        res = _read(h29, _parts(h29, score={f: arm}, de={f: 0.01}))
        assert res["statements"][h29.FULL]["S1"] is ok and res["statements"][h29.FULL]["S2"] is ok
    # one world short: S1 is required in own and in CLEAN
    res = _read(h29, _parts(h29, score={f: {"own": np.full(10, 1.0), "CLEAN": np.full(10, 0.2)}}, de={f: 0.01}))
    assert not res["statements"][h29.FULL]["S1"]
    # S3: both DE terms above 0 in both worlds; zero, or one world below, fails
    for de in (0.0, {"own": 0.01, "CLEAN": -0.001}):
        res = _read(h29, _parts(h29, score={f: np.full(10, 1.0)}, de={f: de}))
        st = res["statements"][h29.FULL]
        assert st["S1"] and st["S2"] and not st["S3"] and not st["qualified"]
        assert res["candidates"]["H-T1"]["why"] == "fails S3"
    # S2 against the placebo: a placebo within the bar of the arm fails S2 only
    res = _read(h29, _both(h29, _good(h29, h29.FULL, placebo_gain=0.25)))
    assert res["statements"][h29.FULL] == {"S1": True, "S2": False, "S3": True, "qualified": False}


def test_every_combination_maps_to_one_outcome_per_candidate(h29):
    import itertools
    for bits in itertools.product([False, True], repeat=6):
        st = {h29.FULL: dict(zip(("S1", "S2", "S3"), bits[:3], strict=True)),
              h29.HALF: dict(zip(("S1", "S2", "S3"), bits[3:], strict=True))}
        for v in st.values():
            v["qualified"] = v["S1"] and v["S2"] and v["S3"]
        c = h29.candidates(st)
        assert (c["H-T1"]["outcome"] == "QUALIFIED") == all(bits[:3])
        assert (c["H-T1B"]["outcome"] == "QUALIFIED") == (all(bits[3:]) and bits[1])
        assert all(v["outcome"] in ("QUALIFIED", "NOT QUALIFIED") for v in c.values())
    assert {k for k in h29.FOLLOWS} == {(c, o) for c in ("H-T1", "H-T1B")
                                        for o in ("QUALIFIED", "NOT QUALIFIED", "VOID")}


def test_void_of_the_whole_read(h29):
    def why(p, **kw):
        r = _read(h29, p, **kw)
        assert r["void"] and [v["outcome"] for v in r["candidates"].values()] == ["VOID", "VOID"]
        assert h29.text(r).count("VOID (the whole read): ") == 1
        return r["why"]
    p = _parts(h29)
    del p[h29.TASK_KEYS[1]]
    assert any(w.startswith("V1: missing transition") for w in why(p))
    p = _parts(h29)
    p[h29.TASK_KEYS[0]].append(_part(h29, h29.TASK_KEYS[0], [3]))
    assert any("occur in more than one part" in w for w in why(p))
    p = _parts(h29)
    del p[h29.TASK_KEYS[0]][0]["rows"][0]["scores"]["own"][h29.name(h29.HALF, "placebo_s1")]
    assert any("missing from scores" in w for w in why(p))
    p = _parts(h29)
    del p[h29.TASK_KEYS[0]][0]["rows"][0]["checks"][h29.HALF]
    assert any("records or checks missing" in w for w in why(p))
    p = _parts(h29)
    p[h29.TASK_KEYS[0]][1]["provenance"]["commit"] = "d" * 40
    assert any("different commits" in w for w in why(p))
    assert any("dirty tree" in w for w in why(_parts(h29, clean=False)))
    assert any("row's commit" in w for w in why(_parts(h29), commit="c" * 39))  # a prefix is not the run commit
    assert any("row's commit" in w for w in why(_parts(h29), commit=""))
    p = _parts(h29)
    p[h29.NOT_BUILT] = [_part(h29, h29.NOT_BUILT, SEEDS)]
    assert any("other transitions" in w for w in why(p))
    p = _parts(h29)
    p[h29.TASK_KEYS[1]][0]["truth_info"]["thin_0"]["cells"] = 1
    assert any(w.startswith("V3:") and "truth thin_0" in w for w in why(p))
    p = _parts(h29)
    p[h29.TASK_KEYS[1]][0]["v7"]["ok"] = False
    assert any("reference cloud" in w for w in why(p))
    p = _parts(h29)
    del p[h29.TASK_KEYS[0]][0]["truth_info"]
    assert any("cannot be read" in w for w in why(p))


def test_a_form_void_voids_its_candidate_only(h29):
    # the half form's gate fails: H-T1 is still read and qualifies; H-T1B is void
    res = _read(h29, _both(h29, _good(h29, h29.FULL), _good(h29, h29.HALF)), h27p=_h27(h29, off={h29.HALF: 0.2}))
    assert not res["void"] and res["form_why"][h29.FULL] == []
    assert res["candidates"]["H-T1"]["outcome"] == "QUALIFIED" and res["candidates"]["H-T1B"]["outcome"] == "VOID"
    assert any(w.startswith("V4: half_iso_half ") for w in res["form_why"][h29.HALF])
    txt = h29.text(res)
    assert "VOID (half_iso_half, H-T1B): V4: half_iso_half" in txt and "GATE full_iso_half" in txt
    assert "S1 H-T1B" not in txt and h29.HALF not in res["tables"]
    # the full form's gate fails: H-T1 is void and H-T1B cannot qualify (H-T1's S2 cannot be read)
    res = _read(h29, _both(h29, _good(h29, h29.FULL), _good(h29, h29.HALF)), h27p=_h27(h29, off={h29.FULL: -0.2}))
    assert res["candidates"]["H-T1"]["outcome"] == "VOID"
    assert res["candidates"]["H-T1B"] == {"outcome": "NOT QUALIFIED",
                                          "why": "fails H-T1's S2 cannot be read (H-T1's form is void)"}


def test_the_gate_at_its_edges_per_form(h29):
    for form in h29.FORMS:
        assert not _read(h29, h27p=_h27(h29, off={form: 0.0999}))["form_why"][form]
        assert _read(h29, h27p=_h27(h29, off={form: 0.1001}))["form_why"][form]
        h = _h27(h29)
        h[h29.TASK_KEYS[0]][0]["rows"][0]["scores"]["own"][h29.h27_arm(form)] = 50.5001  # one seed beyond 0.5
        r = _read(h29, h27p=h)
        assert any("a seed is" in w for w in r["form_why"][form])
        assert not r["form_why"][[f for f in h29.FORMS if f != form][0]]


def test_each_form_void_on_its_identity_record_pool_and_forecast(h29):
    def form_why(p):
        r = _read(h29, p)
        assert not r["void"] and not r["form_why"][h29.FULL] and r["form_why"][h29.HALF]
        assert r["candidates"]["H-T1B"]["outcome"] == "VOID"
        return r["form_why"][h29.HALF]
    for k in ("s1_equals_harness_b0_trend", "same_draw", "pre_size_coords_equal", "post_size_coords_equal"):
        assert any(w.startswith("V5:") for w in form_why(_parts(h29, checks={h29.HALF: {k: False}}))), k
    assert any("cell numbers" in w for w in form_why(_parts(h29, checks={h29.HALF: {"n_cells": {
        "control": 350, "trend_s1": 349, "placebo_s1": 350}}})))
    p = _parts(h29)
    p[h29.TASK_KEYS[0]][0]["rows"][0]["n_cells"][h29.name(h29.HALF, "trend_s1")] = 349
    assert any("cell numbers" in w for w in form_why(p))
    p = _parts(h29)
    p[h29.TASK_KEYS[0]][0]["rows"][0]["skills"]["thin_0"][h29.name(h29.HALF, "trend_s1")]["d2_shape"] += 1e-8
    assert any("d2_shape is not the control's" in w for w in form_why(p))
    p = _parts(h29)
    p[h29.TASK_KEYS[1]][0]["rows"][0]["records"][h29.name(h29.HALF, "trend_s1")]["trend"]["budget_days"] = 0.375
    assert any("budget is" in w for w in form_why(p))
    p = _parts(h29)
    p[h29.TASK_KEYS[1]][0]["rows"][0]["records"][h29.name(h29.HALF, "placebo_s1")]["trend"]["ref_cells"] = 1600
    assert any("not the form's pool" in w for w in form_why(p))
    p = _parts(h29)
    p[h29.TASK_KEYS[1]][0]["rows"][0]["records"][h29.name(h29.HALF, "placebo_s1")]["trend"]["fixed_points"]["CM"] = 1
    assert any("own factor" in w for w in form_why(p))
    p = _parts(h29)
    p[h29.TASK_KEYS[1]][0]["rows"][0]["records"][h29.name(h29.HALF, "placebo_s1")]["trend"]["shifted"] = ["CM"]
    assert any("different lineages" in w for w in form_why(p))
    p = _parts(h29)
    p[h29.TASK_KEYS[1]][0]["rows"][0]["records"][h29.name(h29.HALF, "trend_s1")]["trend"]["strength"] = 0.5
    assert any("not strength 1" in w for w in form_why(p))
    p = _parts(h29)
    p[h29.TASK_KEYS[1]][0]["pred_info"][h29.HALF]["core_share"] = 0.51
    assert any("pool differs" in w for w in form_why(p))
    p = _parts(h29)
    p[h29.TASK_KEYS[1]][0]["forecast_targets"][h29.HALF] *= 1 + 1e-9
    assert any("forecast size differs" in w for w in form_why(p))


def test_folds_merge_the_two_backtests_seed_by_seed(h29):
    a = h29.name(h29.FULL, "trend_s1")
    p = {h29.TASK_KEYS[0]: [_part(h29, h29.TASK_KEYS[0], SEEDS, score={a: np.arange(10, dtype=float)})],
         h29.TASK_KEYS[1]: [_part(h29, h29.TASK_KEYS[1], SEEDS, score={a: -np.arange(10, dtype=float)})]}
    res = _read(h29, p)
    f = res["tables"][h29.FULL]["own"]["fold E10.5"][h29.pname(a, h29.name(h29.FULL, "control"))]["score"]
    assert f["mean"] == 0.0 and f["lo"] == 0.0 and f["hi"] == 0.0
    t = res["tables"][h29.FULL]["own"][h29.TASK_KEYS[0]][h29.pname(a, h29.name(h29.FULL, "control"))]["score"]
    assert t["mean"] == pytest.approx(4.5) and t["ahead"] == 9


# ------------------------------------------------------------------------------------------------ the script
def _prints(fn):
    tree = ast.parse(inspect.getsource(fn).lstrip())
    return [ast.unparse(n) for n in ast.walk(tree) if isinstance(n, ast.Call) and getattr(n.func, "id", "") == "print"]


def test_the_run_and_the_read_print_only_step_lines(h29):
    assert sorted(_prints(h29.run_transition)) == sorted(
        ['print(f"{key} seed {s} ({row[\'seconds\']} s) done", flush=True)', "print(f'wrote {path}', flush=True)"])
    for fn in (h29.cmd_run, h29.build_form, h29.form_pools, h29.reading, h29.text, h29.void_reasons):
        assert _prints(fn) == [], fn
    assert sorted(_prints(h29.cmd_read)) == sorted(["print(f\"OUTCOME {c}: {res['candidates'][c]['outcome']}\")",
                                                    "print(res['follows'][c])"])


def test_the_printed_what_follows_lines_hold_no_figure(h29):
    """The read prints the outcome and its fixed what-follows line: no digit but the statement names S1 to S3 and
    the names of the row, J6, K1, sha256 and the arms, in every outcome either candidate can reach."""
    for c in (_read(h29), _read(h29, _both(h29, _good(h29, h29.FULL), _good(h29, h29.HALF))),
              _read(h29, _parts(h29, clean=False)), _read(h29, h27p=_h27(h29, off={h29.FULL: -0.2})),
              _read(h29, h27p=_h27(h29, off={h29.HALF: 0.2}))):
        for line in c["follows"].values():
            names = r"H29a|H27|S[123]|J6|K1|sha256|E10\.5|H-T1B?|_iso_half:f|trend_s1"
            assert not re.search(r"\d", re.sub(names, "", line)), line


def test_the_script_names_no_test_stage_file_and_reads_no_test_score(h29):
    src = SCRIPT.read_text()
    for bad in ("results/p3", "/p3", "p3/", "size_heart_E12", "submissions", "leaderboard", "registry", "build_test",
                "boards/", "upload", "12.5", "manual_uploads", "index.json", "SIZES", "size_or_hold("):
        assert bad not in src, bad


def _refusing(h29, monkeypatch):
    def boom(*a, **k):
        raise AssertionError("something loaded")
    monkeypatch.setattr(h29.payoff, "heart_task", boom)
    monkeypatch.setattr(h29.h27.sp, "measure_stage", boom)
    monkeypatch.setattr(h29, "load_parts", boom)
    monkeypatch.setattr(h29, "load_h27", boom)


def test_each_run_refuses_before_anything_loads(h29, monkeypatch, tmp_path):
    _refusing(h29, monkeypatch)
    head = "c" * 40
    monkeypatch.setattr(h29, "tries_rows", lambda path=h29.TRIES: {"H28", "H29"})
    with pytest.raises(SystemExit, match="no row"):
        h29.run_transition(h29.TASK_KEYS[0], (0,), tmp_path, tmp_path, head)
    monkeypatch.setattr(h29, "tries_rows", lambda path=h29.TRIES: {"H28", "H29a"})
    monkeypatch.setattr(h29.provenance, "git_state", lambda repo, ignore=(): (head, ["src/x.py"]))
    with pytest.raises(SystemExit, match="not clean"):
        h29.run_transition(h29.TASK_KEYS[0], (0,), tmp_path, tmp_path, head)
    monkeypatch.setattr(h29.provenance, "git_state", lambda repo, ignore=(): ("unknown", None))
    with pytest.raises(SystemExit, match="not clean"):
        h29.run_transition(h29.TASK_KEYS[0], (0,), tmp_path, tmp_path, head)
    monkeypatch.setattr(h29.provenance, "git_state", lambda repo, ignore=(): (head, []))
    for other in ("d" * 40, head[:12], "", None):
        with pytest.raises(SystemExit, match="HEAD is not the run commit"):
            h29.run_transition(h29.TASK_KEYS[0], (0,), tmp_path, tmp_path, other)
    (tmp_path / h29.PART).write_text("{}")
    with pytest.raises(SystemExit, match="exists"):
        h29.run_transition(h29.TASK_KEYS[0], (0,), tmp_path, tmp_path, head)
    with pytest.raises(SystemExit, match="exists"):
        h29.main(["run", "--task", h29.TASK_KEYS[0], "--seeds", "0-1", "--out", str(tmp_path), "--commit", head])
    (tmp_path / h29.PART).unlink()
    with pytest.raises(SystemExit, match="not buildable"):
        h29.run_transition(h29.NOT_BUILT, (0,), tmp_path, tmp_path, head)
    with pytest.raises(AssertionError, match="something loaded"):  # past every refusal, the first load
        h29.run_transition(h29.TASK_KEYS[0], (0,), tmp_path, tmp_path, head)


def test_the_read_refuses_without_row_clean_tree_head_or_with_a_reading(h29, monkeypatch, tmp_path):
    _refusing(h29, monkeypatch)
    head = "c" * 40
    o = tmp_path / "o"
    args = ["read", "--runs", str(tmp_path), "--out", str(o), "--commit", head]
    monkeypatch.setattr(h29, "tries_rows", lambda path=h29.TRIES: {"H28"})
    with pytest.raises(SystemExit, match="no row"):
        h29.main(args)
    monkeypatch.setattr(h29, "tries_rows", lambda path=h29.TRIES: {"H28", "H29a"})
    monkeypatch.setattr(h29.provenance, "git_state", lambda repo, ignore=(): (head, ["x"]))
    with pytest.raises(SystemExit, match="not clean"):
        h29.main(args)
    monkeypatch.setattr(h29.provenance, "git_state", lambda repo, ignore=(): ("d" * 40, []))
    with pytest.raises(SystemExit, match="HEAD is not"):
        h29.main(args)
    monkeypatch.setattr(h29.provenance, "git_state", lambda repo, ignore=(): (head, []))
    o.mkdir()
    for f in ("reading.txt", "reading.json"):
        (o / f).write_text("x")
        with pytest.raises(SystemExit, match="already holds"):
            h29.main(args)
        (o / f).unlink()
    with pytest.raises(AssertionError, match="something loaded"):  # past every refusal, the first load
        h29.main(args)


def test_the_test_budget_and_j6_at_strength_one(h29):
    line = h29.h28.j6_line(h29.STRENGTH)
    assert "test budget 0.5 day" in line and "practice 0.375 day, ratio 1.33" in line
    assert "practice 0.125 day, ratio 4" in line and "first branch" in line
    assert h29.STRENGTHS == (1.0,) and h29.LOCAL == ("control", "trend_s1", "placebo_s1")
    assert h29.PERM_SEED == 20261120 and h29.BAR == 0.3 and h29.MIN_SEEDS == 8
    assert (h29.GUARD_MEAN, h29.GUARD_SEED) == (0.1, 0.5)


# ------------------------------------------------------------------------------------------------ the record
def test_the_gate_dry_runs_on_committed_results(h29):
    """Rule 9 J7: the gate's code and tolerances on committed reproduction data other than the pair it will read:
    row H27's rebuilt oracle twins of A and B against row H26's committed A and B, on the two E10.5 backtests and
    this row's four truths. It could have failed had the draw, the pools or the scorer moved between commits."""
    h27p = h29.load_h27()
    h26p = h29.h26.load_parts(RESULTS / "h26_heart_scope_forms" / "parts")
    for key in h29.TASK_KEYS:
        rows = h29.h26.merge_rows(h27p[key])[0]
        assert sorted(r["seed"] for r in rows) == list(SEEDS)
        for form in h29.FORMS:
            why, big = h29.h28.gate(key, rows, h29.h26.merge_rows(h26p[key])[0],
                                    mine=h29.h27.arm_name(form, "o"), theirs=form)
            assert why == [] and big == {"mean": 0.0, "seed": 0.0}, (key, form, why)
            assert all(h29.h27_arm(form) in r["scores"][q] for r in rows for q in h29.TRUTHS)
            assert h27p[key][0]["forecast_targets"][form] > 0 and form in h27p[key][0]["pred_info"]
    # a shifted copy of the same rows fails the gate (the guard can fail)
    key = h29.TASK_KEYS[0]
    rows = json.loads(json.dumps(h29.h26.merge_rows(h27p[key])[0]))
    for r in rows:
        r["scores"]["own"][h29.name(h29.HALF, "control")] = r["scores"]["own"][h29.h27_arm(h29.HALF)] + 0.2
        for q in h29.TRUTHS:
            r["scores"][q].setdefault(h29.name(h29.HALF, "control"), r["scores"][q][h29.h27_arm(h29.HALF)])
    why, _ = h29.gate(key, h29.HALF, rows, h29.h26.merge_rows(h27p[key])[0])
    assert why and all(w.startswith("V4: half_iso_half ") for w in why)


# ------------------------------------------------------------------------------------------------ the launcher
def test_the_launcher_runs_two_transitions_and_the_tests_first(h29):
    sh = LAUNCHER.read_text()
    assert ("tests/test_h29a_heart_trend_forms.py tests/test_h28_heart_trend_strength.py \\\n"
            "      tests/test_h27_heart_forecast_size.py") in sh
    assert all(f'"{k}:' in sh for k in h29.TASK_KEYS) and h29.NOT_BUILT not in sh
    assert "0-1 2-3 4-5 6-7 8-9" in sh
    for need in ("OMP_NUM_THREADS=8", "OPENBLAS_NUM_THREADS=8", "MKL_NUM_THREADS=8", "NUMPY_MADVISE_HUGEPAGE=0",
                 'TMPDIR="$S/tmp"', 'S="$HOME/v/scratch/h29a"', "status --porcelain", "python -I -B",
                 '--commit "$COMMIT"', "results/h29a_heart_trend_forms"):
        assert need in sh, need
    assert not re.search(r"(^|[\s;&|(])rm\s", sh) and "${" + "2:?" not in sh and "unlink" not in sh
    pat = dict(re.findall(r"^(DONE_LINE|WROTE_LINE)='([^']*)'$", sh, flags=re.M))
    assert set(pat) == {"DONE_LINE", "WROTE_LINE"} and "h29a" in pat["WROTE_LINE"]
    run = sh[sh.index("  run)\n    mkdir"):sh.index("  read)\n    if")]
    assert "grep -cvE" in run and not re.search(r"grep (?!-cvE)", run)
    read = sh[sh.index("  read)\n    if"):]
    assert '> "$S/read.out" 2>&1' in read


def _launch(tmp_path, *args, home=None):
    env = {k: v for k, v in os.environ.items() if not k.startswith("GIT_")}
    env.update(HOME=str(home or tmp_path / "home"), GIT_CEILING_DIRECTORIES=str(tmp_path))
    r = subprocess.run(["bash", str(LAUNCHER), *args], cwd=tmp_path / "repo", env=env, capture_output=True,
                       text=True, timeout=60)
    return r.returncode, r.stderr


@pytest.fixture()
def fake_repo(tmp_path):
    """A throwaway git repository with one commit (and a committed reading), under the test's own folder."""
    repo = tmp_path / "repo"
    out = repo / "experiments" / "harness" / "results" / "h29a_heart_trend_forms"
    out.mkdir(parents=True)
    (repo / "f.txt").write_text("x\n")
    env = {k: v for k, v in os.environ.items() if not k.startswith("GIT_")}
    env.update(HOME=str(tmp_path / "home"), GIT_CEILING_DIRECTORIES=str(tmp_path))
    (tmp_path / "home").mkdir()

    def git(*a):
        return subprocess.run(["git", "-c", "user.name=t", "-c", "user.email=t@t", "-c", "init.defaultBranch=x", *a],
                              cwd=repo, env=env, capture_output=True, text=True, check=True).stdout.strip()
    git("init", "-q")
    git("add", "f.txt")
    git("commit", "-q", "-m", "one")
    return repo, git("rev-parse", "HEAD"), git, out


def test_the_launcher_refusal_paths_under_a_temporary_home(tmp_path, fake_repo):
    repo, head, git, out = fake_repo
    S = tmp_path / "home" / "v" / "scratch" / "h29a"

    def exit_of(step):
        return (S / f"{step}.exit").read_text().strip()
    # no arguments, or one: usage, exit 2, usage.exit written
    assert _launch(tmp_path)[0] == 2 and exit_of("usage") == "2"
    assert _launch(tmp_path, "tests")[0] == 2
    # an unknown step: 97
    assert _launch(tmp_path, "build", head)[0] == 97 and exit_of("build") == "97"
    # run before tests, read before run: 96
    assert _launch(tmp_path, "run", head)[0] == 96 and exit_of("run") == "96"
    assert _launch(tmp_path, "read", head)[0] == 96 and exit_of("read") == "96"
    # tests exited 0 but at another commit: 96
    (S / "tests.exit").write_text("0\n")
    (S / "tests.head").write_text("d" * 40 + "\n")
    assert _launch(tmp_path, "run", head)[0] == 96
    # tests exited 1 at this commit: 96
    (S / "tests.exit").write_text("1\n")
    (S / "tests.head").write_text(head + "\n")
    assert _launch(tmp_path, "run", head)[0] == 96
    # an unclean tree: 98, before HEAD is recorded
    (repo / "g.txt").write_text("y\n")
    code, err = _launch(tmp_path, "tests", head)
    assert code == 98 and exit_of("tests") == "98" and "not clean" in err
    (repo / "g.txt").unlink()
    # another HEAD: 99
    code, err = _launch(tmp_path, "tests", "d" * 40)
    assert code == 99 and exit_of("tests") == "99" and "HEAD is not" in err
    assert (S / "tests.head").read_text().strip() == head  # the earlier head file stays as it was written
    # a reading already there: read refuses (95) and leaves it as it was
    (out / "reading.txt").write_text("the one reading\n")
    git("add", "-f", str(out / "reading.txt"))
    git("commit", "-q", "-m", "reading")
    head2 = git("rev-parse", "HEAD")
    (S / "run.exit").write_text("0\n")
    (S / "run.head").write_text(head2 + "\n")
    code, err = _launch(tmp_path, "read", head2)
    assert code == 95 and exit_of("read") == "95" and "already holds a reading" in err
    assert (out / "reading.txt").read_text() == "the one reading\n"
    assert (S / "read.head").read_text().strip() == head2
