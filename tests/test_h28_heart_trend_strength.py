"""Try H28 (experiments/harness/h28_heart_trend_strength.py): synthetic arrays and hand-made tables only, and committed
result files of earlier rows read for arithmetic. Nothing here loads a competition file or scores a real cell."""

from __future__ import annotations

import ast
import functools
import importlib.util
import inspect
import json
import re
import sys
from pathlib import Path

import numpy as np
import pytest

from vec import harness, t1_scope
from vec.adoption import SEEDS

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "experiments" / "harness" / "h28_heart_trend_strength.py"
RESULTS = REPO_ROOT / "experiments" / "harness" / "results"


@pytest.fixture(scope="module")
def h28():
    spec = importlib.util.spec_from_file_location("h28_heart_trend_strength_module", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


# ------------------------------------------------------------------------------------------------ synthetic stages
@pytest.fixture(scope="module")
def synth():
    spec = importlib.util.spec_from_file_location("b0_trend_tests_for_h28", REPO_ROOT / "tests" / "test_b0_trend.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    inputs = {"ref": mod._stage(1600, seed=2, shift=0.0, tiny=60),
              "start": mod._stage(1100, seed=3, shift=0.8, tiny=10),
              "earlier": mod._stage(900, seed=5, shift=1.6, tiny=10)}
    for a in inputs.values():
        a.obs["core"] = np.asarray(a.obs["merged"]).astype(str) == "CM"
    return mod, inputs


def _trimmed(h28, inputs):
    pool = inputs["ref"]
    PC = np.asarray(pool.obsm["spatial_3D"], dtype=float)
    keep = h28.h26.pool_masks(PC, np.asarray(pool.obs["core"], bool))[h28.FORM]
    return {**inputs, "ref": pool[keep].copy()}


def _size(h28, inputs):
    sp = h28.h27.sp
    geoms = [sp.geometry_from_arrays(n, t, np.asarray(inputs[r].obsm["spatial_3D"], float),
                                     np.asarray(inputs[r].obs["core"], bool))
             for r, n, t in (("earlier", "E8.5", 8.5), ("start", "E8.75", 8.75), ("ref", "E9.5", 9.5))]
    return h28.h27.analogue_size_result(geoms, 10.5)


def test_strength_one_is_the_harness_b0_trend_and_zero_keeps_every_count(h28, synth):
    mod, inputs = synth
    for ins in (inputs, _trimmed(h28, inputs)):
        ctx = mod._ctx(ins, seed=harness.BUILD_SEED_OFFSET + 4)
        st1, _, rows1, rec1 = h28.trend_expression(ctx, 1.0)
        ref = harness._b0_t2(mod._ctx(ins, seed=harness.BUILD_SEED_OFFSET + 4), trend=True)
        np.testing.assert_array_equal(st1.X, ref.X)
        np.testing.assert_array_equal(st1.coords, ref.coords)
        st0, c0, rows0, _ = h28.trend_expression(ctx, 0.0)
        pool, rows, X, coords = harness.draw_copy(ctx)
        np.testing.assert_array_equal(rows0, rows)
        assert (c0 != h28.counts_of(pool[rows], list(ctx.board.genes))).nnz == 0
        ctrl = harness._copy(ctx)
        np.testing.assert_array_equal(st0.coords, ctrl.coords)
        # the trend path re-renders X from the counts: equal to float32 rounding on these synthetic stages
        np.testing.assert_allclose(st0.X, ctrl.X, rtol=0, atol=1e-5)
        assert rec1["ref_cells"] == ins["ref"].n_obs and rec1["n_trend"] >= 1


def test_the_strength_scales_the_mapping_and_the_fallback_damping(h28, synth, monkeypatch):
    """Strength s multiplies b0_trend's mapping and momentum_shift's damping (B0's own factor for a lineage the
    mapping leaves out); the placebo maps every shifted lineage, fallbacks included."""
    mod, inputs = synth
    early = inputs["earlier"].copy()
    m = np.asarray(early.obs["merged"]).astype(str)
    early.obs["merged"] = np.where(m == "Endothelial", "Elsewhere", m)  # Endothelial usable at start and ref only
    ins = {**inputs, "earlier": early}
    seen = []
    real = t1_scope.momentum_shift

    @functools.wraps(real)  # b0_trend_factors reads momentum_shift's defaults through its signature
    def spy(*a, **k):
        seen.append((a[3], {key: np.array(v) for key, v in k["logf_by_type"].items()}))
        return real(*a, **k)

    monkeypatch.setattr(h28.t1_scope, "momentum_shift", spy)
    ctx = mod._ctx(ins, seed=harness.BUILD_SEED_OFFSET + 2)
    for s in (1.0, 0.5, 0.25, 0.0):
        h28.trend_expression(ctx, s)
    _, _, _, rec = h28.trend_expression(ctx, 0.5, placebo_arm=True)
    assert rec["fallback"] == ["Endothelial"]
    (d1, m1), (d5, m5), (d25, m25), (d0, m0), (dp, mp) = seen
    assert (d1, d5, d25, d0) == (0.5, 0.25, 0.125, 0.0)
    assert set(m1) == {"CM"}
    np.testing.assert_array_equal(m5["CM"], 0.5 * m1["CM"])
    np.testing.assert_array_equal(m25["CM"], 0.25 * m1["CM"])
    assert not np.any(m0["CM"])
    assert set(mp) == {"CM", "Endothelial"} and rec["placebo"] and rec["fixed_points"] == {"CM": 0, "Endothelial": 0}
    np.testing.assert_allclose(np.sort(mp["CM"]), np.sort(m5["CM"]), rtol=0, atol=0)
    assert not np.array_equal(mp["CM"], m5["CM"])
    # the fallback's placebo carries B0's own factor at strength 0.5, deranged
    own = harness.B0_DAMPING * 0.5 * (h28.lineage_log_means(ctx, "ref", list(ins["ref"].var_names), ["Endothelial"])
                                      ["Endothelial"]
                                      - h28.lineage_log_means(ctx, "start", list(ins["ref"].var_names),
                                                              ["Endothelial"])["Endothelial"])
    np.testing.assert_allclose(np.sort(mp["Endothelial"]), np.sort(own), rtol=1e-12, atol=1e-15)


def test_the_placebo_derangement_moves_every_gene_within_its_level_bin(h28):
    basis = np.random.default_rng(0).normal(size=500)
    p = h28.derangement(basis, "CM")
    assert sorted(p) == list(range(500)) and not np.any(p == np.arange(500))
    bins = np.array_split(np.argsort(basis, kind="stable"), h28.PERM_BINS)
    where = np.empty(500, int)
    for i, b in enumerate(bins):
        where[b] = i
    assert np.array_equal(where[p], where)
    np.testing.assert_array_equal(p, h28.derangement(basis, "CM"))
    assert not np.array_equal(p, h28.derangement(basis, "Endocardium"))


def test_every_arm_has_the_controls_cells_coordinates_and_half(h28, synth):
    mod, inputs = synth
    ins = _trimmed(h28, inputs)
    size = _size(h28, inputs)
    ctx = mod._ctx(ins, seed=harness.BUILD_SEED_OFFSET + 3)
    out, extra = h28.build_arms(ctx, size, whole_ref=inputs["ref"])
    c = extra["checks"]
    assert set(out) == set(h28.BUILT)
    for k in ("s1_equals_harness_b0_trend", "s0_counts_equal_drawn", "same_draw", "pre_size_coords_equal",
              "post_size_coords_equal"):
        assert c[k] is True, k
    assert len(set(c["n_cells"].values())) == 1 and c["n_cells"][h28.CONTROL] == mod.BOARD.max_cells // 2
    for name, st in out.items():
        np.testing.assert_array_equal(st.coords, out[h28.CONTROL].coords)
        assert abs(extra["records"][name]["size"]["log_rms_after_vs_target"]) <= 1e-9
    assert not np.array_equal(out[h28.arm(1.0)].X, out[h28.CONTROL].X)
    assert not np.array_equal(out[h28.arm(1.0)].X, out[h28.placebo(1.0)].X)
    assert "CM" in c["whole_pool"]


# ------------------------------------------------------------------------------------------------ hand-made tables
METRICS = ("de_score", "de_direction", "mmd_u", "variogram", "d2_shape", "occupancy_dice", "scale_log_ratio",
           "neighborhood_mmd")


def _info(h28):
    t = {q: {"cells": 1000 + i, "core_share": 0.4, "xy_scope": 1.0} for i, q in enumerate(h28.TRUTHS)}
    return t, {h28.FORM: {"pool_cells": 500, "core_share": 0.9, "xy_scope": 0.2}}


def _rec(h28, key, name):
    if name == h28.CONTROL:
        return {"size": {}, "trend": {"strength": None}}
    s = float(name.split("_s")[1])
    return {"size": {}, "trend": {"strength": s, "placebo": name.startswith("placebo"),
                                  "budget_days": h28.PRACTICE_BUDGET.get(key, 0.5), "ref_cells": 500,
                                  "shifted": ["CM", "Endocardium"], "fixed_points": {"CM": 0, "Endocardium": 0}}}


def _part(h28, key, seeds, score=None, de=None, commit="abc123", clean=True, checks=None):
    """score: name -> per-seed array of board points added to 50 in every truth (or a dict world -> array);
    de: name -> skill added to de_score and de_direction."""
    score, de = score or {}, de or {}
    t_info, p_info = _info(h28)
    rows = []
    for s in seeds:
        sc, sk, gr = {}, {}, {}
        for q in h28.TRUTHS:
            w = "own" if q == h28.OWN else "CLEAN"
            sc[q], sk[q], gr[q] = {}, {}, {}
            for name in h28.ALL:
                v = score.get(name, 0.0)
                v = v.get(w, 0.0) if isinstance(v, dict) else v
                sc[q][name] = 50.0 + float(np.broadcast_to(v, (10,))[s])
                sk[q][name] = {m: 0.5 + (de.get(name, 0.0) if m in ("de_score", "de_direction") else 0.0)
                               for m in METRICS}
                gr[q][name] = {g: 0.5 for g in ("expression_change", "cell_state", "shape_growth",
                                                 "local_organisation")}
        ck = {"s1_equals_harness_b0_trend": True, "s0_counts_equal_drawn": True, "same_draw": True,
              "pre_size_coords_equal": True, "post_size_coords_equal": True, "s0_x_max_abs_diff_control": 1e-7,
              "n_cells": {n: 350 for n in h28.BUILT}}
        ck.update(checks or {})
        rows.append({"seed": s, "n_cells": ck["n_cells"], "checks": ck, "scores": sc, "skills": sk, "groups": gr,
                     "records": {n: _rec(h28, key, n) for n in h28.BUILT}})
    return {"axis": h28.AXIS, "task": key, "rows": rows, "truth_info": t_info, "pred_info": p_info,
            "forecast_targets": {h28.FORM: 100.0}, "v7": {"ok": True, "log_diff": 0.0},
            "provenance": {"commit": commit, "tree_clean": clean}}


def _parts(h28, **kw):
    return {k: [_part(h28, k, range(0, 5), **kw), _part(h28, k, range(5, 10), **kw)] for k in h28.TASK_KEYS}


def _h27(h28, off=0.0):
    t_info, p_info = _info(h28)
    rows = [{"seed": s, "scores": {q: {h28.H27_ARM: 50.0 + off} for q in h28.TRUTHS}} for s in SEEDS]
    return {k: [{"rows": rows, "truth_info": t_info, "pred_info": p_info, "forecast_targets": {h28.FORM: 100.0}}]
            for k in h28.TASK_KEYS}


def _read(h28, parts=None, h27p=None, commit="abc"):
    return h28.reading(parts if parts is not None else _parts(h28), h27p if h27p is not None else _h27(h28), commit)


def test_a_clean_null_run_is_read_and_names_no_strength(h28):
    res = _read(h28)
    assert not res["void"] and res["outcome"] == "NONE" and res["output"] is None
    assert res["follows"] == h28.FOLLOWS["NONE"]
    assert res["guards"]["V4"][h28.TASK_KEYS[0]] == {"mean": 0.0, "seed": 0.0}
    txt = h28.text(res)
    assert txt.splitlines()[0].startswith("Try H28 (commit abc123). A diagnostic on one truth")
    assert h28.NOT_BUILT_LINE in txt and "OUTCOME: NONE" in txt
    assert all(ln.startswith(h28.TAG) for ln in txt.split("TABLES")[1].splitlines()[1:])


def _ahead(h28, s, gain=0.5, placebo_gain=0.0, de=0.01, over=None):
    score = {h28.arm(s): np.full(10, gain), h28.placebo(s): np.full(10, placebo_gain)}
    score.update(over or {})
    return _parts(h28, score=score, de={h28.arm(s): de})


def test_the_statements_and_their_eight_of_ten_edges(h28):
    # exactly 8 of 10 above zero with a mean above the bar passes; 7 of 10 does not; a mean of exactly 0.3 does not
    eight = np.array([1.0] * 8 + [-0.5] * 2)  # mean 0.7
    seven = np.array([1.0] * 7 + [-0.1] * 3)  # mean 0.67
    assert h28.passes(h28._stats(eight)) and not h28.passes(h28._stats(seven))
    assert not h28.passes(h28._stats(np.full(10, 0.3))) and h28.passes(h28._stats(np.full(10, 0.3001)))
    res = _read(h28, _parts(h28, score={h28.arm(1.0): eight, h28.placebo(1.0): np.zeros(10)},
                            de={h28.arm(1.0): 0.01}))
    assert res["outcome"] == "AHEAD" and res["output"] == 1.0
    st = res["statements"]["1"]
    assert st["ahead"] and st["over_placebo"] and st["on_de"] and st["qualified"]
    lines = res["lines"]
    assert lines[0].startswith("S1 strength 1 AHEAD of the copy on the full trim: TRUE (fold E10.5; own +0.70 "
                               "[-0.50, +1.00] 8/10")
    assert "bar: mean above +0.3 with at least 8 of 10 seeds above 0, in both worlds)" in lines[0]
    assert lines[1].startswith("S2 strength 1 AHEAD of its placebo: TRUE")
    assert lines[2].startswith("S3 strength 1 gains on the DE terms: TRUE")
    assert any(ln.startswith("S4 strength 0.5 AHEAD of strength 1: FALSE") for ln in lines)
    res7 = _read(h28, _parts(h28, score={h28.arm(1.0): seven}, de={h28.arm(1.0): 0.01}))
    assert res7["outcome"] == "NONE" and not res7["statements"]["1"]["ahead"]
    # one world short: ahead in own only is not ahead
    one = _read(h28, _parts(h28, score={h28.arm(1.0): {"own": np.full(10, 1.0), "CLEAN": np.full(10, 0.2)}},
                            de={h28.arm(1.0): 0.01}))
    assert one["outcome"] == "NONE"


def test_each_outcome_and_its_one_follows_line(h28):
    ahead = _read(h28, _ahead(h28, 1.0))
    assert ahead["outcome"] == "AHEAD" and ahead["follows"] == h28.FOLLOWS["AHEAD"].format(s=1.0, b=0.5)
    assert "ONE priced question" in ahead["follows"] and "amendment of 2026-11-18" in ahead["follows"]
    assert "range from zero" in ahead["follows"] and "never a point" in ahead["follows"]
    plc = _read(h28, _ahead(h28, 1.0, placebo_gain=0.4))
    assert plc["outcome"] == "PLACEBO" and plc["follows"] == h28.FOLLOWS["PLACEBO"]
    off = _read(h28, _ahead(h28, 1.0, de=-0.01))
    assert off["outcome"] == "OFF_DE" and off["follows"] == h28.FOLLOWS["OFF_DE"]
    assert _read(h28)["follows"] == h28.FOLLOWS["NONE"]
    void = _read(h28, h27p=_h27(h28, off=0.2))
    assert void["outcome"] == "VOID" and void["follows"].startswith("VOID: V4:")
    assert h28.text(void) == void["follows"] + "\n"
    assert set(h28.FOLLOWS) == {"VOID", "AHEAD", "PLACEBO", "OFF_DE", "NONE"}
    # every combination of the three statements over the three strengths maps to exactly one outcome
    import itertools
    for bits in itertools.product([False, True], repeat=9):
        st = {}
        for i, s in enumerate(h28.STRENGTHS):
            a, p, d = bits[3 * i:3 * i + 3]
            st[s] = {"ahead": a, "over_placebo": p, "on_de": d, "qualified": a and p and d,
                     "least_world_mean": 1.0 - i, "beats": {o: False for o in h28.STRENGTHS if o > s}}
        o = h28.outcome(st, False)
        assert o in h28.FOLLOWS and o != "VOID"
        assert (o == "AHEAD") == any(st[s]["qualified"] for s in h28.STRENGTHS)


def test_a_smaller_strength_is_named_only_if_it_beats_the_default(h28):
    # every strength qualified, 0.5 not beating 1 by the bar: the default, strength 1
    sc = {h28.arm(1.0): np.full(10, 0.6), h28.arm(0.5): np.full(10, 0.8), h28.arm(0.25): np.full(10, 0.7)}
    de = {h28.arm(s): 0.01 for s in h28.STRENGTHS}
    assert _read(h28, _parts(h28, score=sc, de=de))["output"] == 1.0
    # 0.5 beats 1 by more than the bar in both worlds: 0.5
    sc[h28.arm(0.5)] = np.full(10, 1.0)
    assert _read(h28, _parts(h28, score=sc, de=de))["output"] == 0.5
    # strength 1 not ahead: the next qualified strength is the default
    sc = {h28.arm(1.0): np.full(10, 0.1), h28.arm(0.5): np.full(10, 0.5)}
    r = _read(h28, _parts(h28, score=sc, de=de))
    assert r["output"] == 0.5 and r["follows"] == h28.FOLLOWS["AHEAD"].format(s=0.5, b=0.25)


def test_void_on_missing_repeated_mixed_or_dirty_parts(h28):
    p = _parts(h28)
    del p[h28.TASK_KEYS[1]]
    assert any(w.startswith("V1: missing transition") for w in _read(h28, p)["why"])
    p = _parts(h28)
    p[h28.TASK_KEYS[0]].append(_part(h28, h28.TASK_KEYS[0], [3]))
    assert any("occur in more than one part" in w for w in _read(h28, p)["why"])
    p = _parts(h28)
    del p[h28.TASK_KEYS[0]][0]["rows"][0]["scores"]["own"][h28.placebo(0.5)]
    assert any("missing from scores" in w for w in _read(h28, p)["why"])
    p = _parts(h28)
    p[h28.TASK_KEYS[0]][1]["provenance"]["commit"] = "zzz"
    assert any("different commits" in w for w in _read(h28, p)["why"])
    p = _parts(h28, clean=False)
    assert any("dirty tree" in w for w in _read(h28, p)["why"])
    assert any("row's commit" in w for w in _read(h28, commit="fff")["why"])
    assert any("row's commit" in w for w in _read(h28, commit="")["why"])
    p = _parts(h28)
    p[h28.NOT_BUILT] = [_part(h28, h28.NOT_BUILT, SEEDS)]
    assert any("other transitions" in w for w in _read(h28, p)["why"])


def test_void_when_a_truth_the_pool_or_the_forecast_differs_from_h27(h28):
    for edit in (lambda p: p["truth_info"]["thin_0"].update(cells=1),
                 lambda p: p["pred_info"][h28.FORM].update(core_share=0.91),
                 lambda p: p["forecast_targets"].update({h28.FORM: 100.0 * (1 + 1e-9)}),
                 lambda p: p["v7"].update(ok=False)):
        p = _parts(h28)
        edit(p[h28.TASK_KEYS[1]][0])
        r = _read(h28, p)
        assert r["void"] and any(w.startswith("V3:") for w in r["why"]), r["why"]


def test_the_gate_at_its_edges(h28):
    assert not _read(h28, h27p=_h27(h28, off=0.0999))["void"]
    assert _read(h28, h27p=_h27(h28, off=0.1001))["void"]
    h = _h27(h28)
    h[h28.TASK_KEYS[0]][0]["rows"][0]["scores"]["own"][h28.H27_ARM] = 50.5001  # one seed beyond 0.5, mean 0.05
    res = _read(h28, h27p=h)
    assert res["void"] and any("a seed is" in w for w in res["why"])


def test_void_on_every_identity_and_record_failure(h28):
    for k in ("s1_equals_harness_b0_trend", "s0_counts_equal_drawn", "same_draw", "pre_size_coords_equal",
              "post_size_coords_equal"):
        r = _read(h28, _parts(h28, checks={k: False}))
        assert r["void"] and any(w.startswith("V5:") for w in r["why"]), k
    p = _parts(h28)
    p[h28.TASK_KEYS[0]][0]["rows"][0]["n_cells"][h28.arm(0.5)] = 349
    assert any("cell numbers differ" in w for w in _read(h28, p)["why"])
    p = _parts(h28)
    p[h28.TASK_KEYS[0]][0]["rows"][0]["skills"]["thin_0"][h28.arm(0.25)]["d2_shape"] += 1e-8
    assert any("d2_shape is not the control's" in w for w in _read(h28, p)["why"])
    p = _parts(h28)
    p[h28.TASK_KEYS[1]][0]["rows"][0]["records"][h28.arm(1.0)]["trend"]["budget_days"] = 0.375
    assert any("budget is" in w for w in _read(h28, p)["why"])
    p = _parts(h28)
    p[h28.TASK_KEYS[1]][0]["rows"][0]["records"][h28.placebo(1.0)]["trend"]["ref_cells"] = 900
    assert any("not the trimmed pool" in w for w in _read(h28, p)["why"])
    p = _parts(h28)
    p[h28.TASK_KEYS[1]][0]["rows"][0]["records"][h28.placebo(1.0)]["trend"]["fixed_points"]["CM"] = 1
    assert any("own factor" in w for w in _read(h28, p)["why"])
    p = _parts(h28)
    p[h28.TASK_KEYS[1]][0]["rows"][0]["records"][h28.arm(0.25)]["trend"]["shifted"] = ["CM"]
    assert any("different lineages" in w for w in _read(h28, p)["why"])
    p = _parts(h28)
    del p[h28.TASK_KEYS[1]][0]["rows"][0]["checks"]
    assert any("records or checks missing" in w for w in _read(h28, p)["why"])


def test_a_malformed_part_voids_and_prints_only_the_void_line(h28):
    p = _parts(h28)
    del p[h28.TASK_KEYS[0]][0]["truth_info"]
    r = _read(h28, p)
    assert r["void"] and r["outcome"] == "VOID" and h28.text(r).startswith("VOID: ")
    assert "no rerun under it; the ruled file stands" in h28.text(r)


def test_folds_merge_the_two_backtests_seed_by_seed(h28):
    sc = {h28.arm(1.0): np.arange(10, dtype=float)}
    p = {h28.TASK_KEYS[0]: [_part(h28, h28.TASK_KEYS[0], SEEDS, score=sc)],
         h28.TASK_KEYS[1]: [_part(h28, h28.TASK_KEYS[1], SEEDS, score={h28.arm(1.0): -np.arange(10, dtype=float)})]}
    res = _read(h28, p)
    f = res["tables"]["own"]["fold E10.5"][h28.pname(h28.arm(1.0), h28.CONTROL)]["score"]
    assert f["mean"] == 0.0 and f["lo"] == 0.0 and f["hi"] == 0.0
    t = res["tables"]["own"][h28.TASK_KEYS[0]][h28.pname(h28.arm(1.0), h28.CONTROL)]["score"]
    assert t["mean"] == pytest.approx(4.5) and t["ahead"] == 9


# ------------------------------------------------------------------------------------------------ the script
def test_the_run_prints_no_score(h28):
    for fn in (h28.run_transition, h28.cmd_run, h28.build_arms, h28.trend_expression, h28.whole_pool_comparison):
        tree = ast.parse(inspect.getsource(fn).lstrip())
        prints = [ast.unparse(n) for n in ast.walk(tree) if isinstance(n, ast.Call) and getattr(n.func, "id", "") ==
                  "print"]
        for p in prints:
            assert p in ('print(f"{key} seed {s} ({row[\'seconds\']} s) done", flush=True)',
                         "print(f'wrote {path}', flush=True)"), p
        assert len(prints) == (2 if fn is h28.run_transition else 0), fn


def test_the_read_prints_only_the_outcome_and_what_follows(h28):
    tree = ast.parse(inspect.getsource(h28.cmd_read).lstrip())
    prints = [ast.unparse(n) for n in ast.walk(tree) if isinstance(n, ast.Call) and getattr(n.func, "id", "") ==
              "print"]
    assert len(prints) == 2 and prints[1] == "print(res['follows'])"
    assert "OUTCOME: " in prints[0] and "res['outcome']" in prints[0] and "res['output']" in prints[0]
    assert "tables" not in prints[0] and "lines" not in prints[0]


def test_the_script_names_no_test_stage_file_and_reads_no_test_score(h28):
    src = SCRIPT.read_text()
    for bad in ("results/p3", "/p3", "p3/", "size_heart_E12", "submissions", "leaderboard", "registry", "build_test",
                "boards/", "upload", "12.5", "46.4", "SIZES", "size_or_hold("):
        assert bad not in src, bad


def test_each_run_refuses_before_anything_loads(h28, monkeypatch, tmp_path):
    def boom(*a, **k):
        raise AssertionError("something loaded")

    monkeypatch.setattr(h28.payoff, "heart_task", boom)
    monkeypatch.setattr(h28.h27.sp, "measure_stage", boom)
    monkeypatch.setattr(h28, "tries_rows", lambda path=h28.TRIES: {"H27"})
    with pytest.raises(SystemExit, match="no row"):
        h28.run_transition(h28.TASK_KEYS[0], (0,), tmp_path, tmp_path)
    monkeypatch.setattr(h28, "tries_rows", lambda path=h28.TRIES: {"H27", "H28"})
    monkeypatch.setattr(h28.provenance, "git_state", lambda repo, ignore=(): ("c" * 40, ["src/x.py"]))
    with pytest.raises(SystemExit, match="not clean"):
        h28.run_transition(h28.TASK_KEYS[0], (0,), tmp_path, tmp_path)
    monkeypatch.setattr(h28.provenance, "git_state", lambda repo, ignore=(): ("unknown", None))
    with pytest.raises(SystemExit, match="not clean"):
        h28.run_transition(h28.TASK_KEYS[0], (0,), tmp_path, tmp_path)
    monkeypatch.setattr(h28.provenance, "git_state", lambda repo, ignore=(): ("c" * 40, []))
    (tmp_path / h28.PART).write_text("{}")
    with pytest.raises(SystemExit, match="exists"):
        h28.run_transition(h28.TASK_KEYS[0], (0,), tmp_path, tmp_path)
    (tmp_path / h28.PART).unlink()
    with pytest.raises(SystemExit, match="not buildable"):
        h28.run_transition(h28.NOT_BUILT, (0,), tmp_path, tmp_path)
    with pytest.raises(AssertionError, match="something loaded"):  # past every refusal, the first load
        h28.run_transition(h28.TASK_KEYS[0], (0,), tmp_path, tmp_path)


def test_the_read_refuses_without_the_row_or_with_a_reading(h28, monkeypatch, tmp_path):
    def boom(*a, **k):
        raise AssertionError("something loaded")

    monkeypatch.setattr(h28, "load_parts", boom)
    monkeypatch.setattr(h28, "load_h27", boom)
    monkeypatch.setattr(h28, "tries_rows", lambda path=h28.TRIES: {"H27"})
    with pytest.raises(SystemExit, match="no row"):
        h28.main(["read", "--runs", str(tmp_path), "--out", str(tmp_path / "o"), "--commit", "abc"])
    monkeypatch.setattr(h28, "tries_rows", lambda path=h28.TRIES: {"H27", "H28"})
    (tmp_path / "o").mkdir()
    (tmp_path / "o" / "reading.txt").write_text("x")
    with pytest.raises(SystemExit, match="already holds"):
        h28.main(["read", "--runs", str(tmp_path), "--out", str(tmp_path / "o"), "--commit", "abc"])


def test_the_test_and_practice_budgets_and_j6(h28):
    assert h28.TEST_BUDGET == 0.5
    line = h28.j6_line(0.25)
    assert "test budget 0.125 day" in line and "practice 0.09375 day, ratio 1.33" in line
    assert "practice 0.03125 day, ratio 4" in line and "first branch" in line


# ------------------------------------------------------------------------------------------------ the record
def test_the_gate_dry_runs_on_committed_results(h28):
    """Rule 9 J7: the gate's code and tolerances on committed reproduction data other than the pair it will read:
    row H27's rebuilt oracle twin of C against row H26's committed C, on the two E10.5 backtests and the truths of
    this row (distance 0.000; it could have failed had the draw, the pools or the scorer moved between commits)."""
    h27p = h28.load_h27()
    h26p = h28.h26.load_parts(RESULTS / "h26_heart_scope_forms" / "parts")
    for key in h28.TASK_KEYS:
        rows = h28.h26.merge_rows(h27p[key])[0]
        why, big = h28.gate(key, rows, h28.h26.merge_rows(h26p[key])[0], mine=h28.h27.arm_name(h28.FORM, "o"),
                            theirs=h28.FORM)
        assert why == [] and big == {"mean": 0.0, "seed": 0.0}
        assert sorted(r["seed"] for r in rows) == list(SEEDS)
        assert all(h28.H27_ARM in r["scores"][q] for r in rows for q in h28.TRUTHS)
        assert h27p[key][0]["forecast_targets"][h28.FORM] > 0
    # a shifted copy of the same rows fails the gate (the guard can fail)
    key = h28.TASK_KEYS[0]
    rows = json.loads(json.dumps(h28.h26.merge_rows(h27p[key])[0]))
    for r in rows:
        r["scores"]["own"][h28.CONTROL] = r["scores"]["own"][h28.H27_ARM] + 0.2
    why, _ = h28.gate(key, rows, h28.h26.merge_rows(h27p[key])[0], truths=(h28.OWN,))
    assert why


def test_the_noise_floor_of_rule_9_j5(h28):
    """J5 on one truth: the seed SE of the ruled file (row H27's C at f) minus the floor on fold E10.5, own and CLEAN,
    and of b0_trend minus the copy on row H25's committed backtests (whole pool)."""
    h27p = h28.load_h27()
    rows = {k: h28.h26.merge_rows(h27p[k])[0] for k in h28.TASK_KEYS}
    se = {}
    for w, ts in h28.WORLDS:
        d = h28.series(rows, h28.TASK_KEYS, ts, h28.H27_ARM) - h28.series(rows, h28.TASK_KEYS, ts, harness.FLOOR)
        se[w] = d.std(ddof=1) / np.sqrt(10)
    h25 = {}
    for p in sorted((RESULTS / "h25_heart_forms" / "tasks").glob("*.json")):
        r = json.loads(p.read_text())
        h25[r["task"]["key"]] = np.array([x["scores"]["b0_trend"] - x["scores"]["copy"]
                                          for x in sorted(r["rows"], key=lambda x: x["seed"])])
    fold = (h25["T2h/<=E9.5->E10.5"] + h25["T2h/<=E8.75->E10.5"]) / 2
    se25 = fold.std(ddof=1) / np.sqrt(10)
    print(f"\nH28 J5: seed SE own {se['own']:.3f}, CLEAN {se['CLEAN']:.3f}; b0_trend - copy (H25) fold mean "
          f"{fold.mean():.3f}, SE {se25:.3f}")
    assert round(se["own"], 3) == 0.111 and round(se["CLEAN"], 3) == 0.163
    assert round(fold.mean(), 3) == 1.105 and round(se25, 3) == 0.089
    assert 2 * max(se.values()) > h28.BAR  # the bar is under twice the floor: a diagnostic, as the row says


def _null_tables(h28, rng, gains, pgain=0.0, sd=0.5, rho=0.5):
    """Hand-made fold tables: per world, per seed, every arm and placebo's difference over the control as gain plus
    noise shared across arms (rho) and worlds (rho), with seed sd `sd` (SE about sd / sqrt(10))."""
    common = rng.normal(size=10)
    tab = {}
    for w, _ in h28.WORLDS:
        wc = np.sqrt(rho) * common + np.sqrt(1 - rho) * rng.normal(size=10)
        x = {}
        for s in h28.STRENGTHS:
            for name, g in ((h28.arm(s), gains[s]), (h28.placebo(s), pgain)):
                x[name] = g + sd * (np.sqrt(rho) * wc + np.sqrt(1 - rho) * rng.normal(size=10))
        fold = {}
        for a, b in h28.PAIRS:
            if b == h28.CONTROL and a in x:
                fold[h28.pname(a, b)] = {"score": h28._stats(x[a]), "de_score": 0.01, "de_direction": 0.01}
            elif a in x and b in x:
                fold[h28.pname(a, b)] = {"score": h28._stats(x[a] - x[b])}
        tab[w] = {"fold E10.5": fold}
    return tab


def test_the_output_rule_on_pure_noise(h28):
    """The challenger's point: the output rule on simulated arms with no effect, and with one equal effect at every
    strength (the winner's curse among nested arms). Seed sd 0.5 (SE about 0.16, the larger J5 figure), noise shared
    across arms and worlds at rho 0.5 and 0.8. Quoted in the row; the read does not use these figures."""
    rng = np.random.default_rng(20261120)
    n = 4000
    got = {}
    for rho in (0.5, 0.8):
        null = sum(h28.output(h28.statements(_null_tables(h28, rng, {s: 0.0 for s in h28.STRENGTHS}, rho=rho)))
                   is not None for _ in range(n)) / n
        eq = [h28.output(h28.statements(_null_tables(h28, rng, {s: 0.8 for s in h28.STRENGTHS}, rho=rho)))
              for _ in range(n)]
        got[rho] = (null, sum(o is not None and o < 1.0 for o in eq) / n, sum(o == 1.0 for o in eq) / n)
        print(f"\nH28 null (rho {rho}): a strength named in {got[rho][0]:.4f} of {n} pure-noise runs; with +0.8 at "
              f"every strength a smaller strength is named in {got[rho][1]:.4f}, strength 1 in {got[rho][2]:.4f}")
    assert all(v[0] < 0.01 for v in got.values())
    assert all(v[1] < 0.10 for v in got.values()) and got[0.8][1] < got[0.5][1]


def test_the_launcher_runs_two_transitions_and_the_tests_first(h28):
    sh = (REPO_ROOT / "experiments" / "harness" / "run_h28.sh").read_text()
    assert "tests/test_h28_heart_trend_strength.py tests/test_h27_heart_forecast_size.py" in sh
    assert all(f'"{k}:' in sh for k in h28.TASK_KEYS) and h28.NOT_BUILT not in sh
    assert "0-1 2-3 4-5 6-7 8-9" in sh
    for need in ("OMP_NUM_THREADS=8", "OPENBLAS_NUM_THREADS=8", "MKL_NUM_THREADS=8", "NUMPY_MADVISE_HUGEPAGE=0",
                 'TMPDIR="$S/tmp"', "S=~/v/scratch/h28", "status --porcelain", "python -I -B", "--commit"):
        assert need in sh, need
    pat = dict(re.findall(r"^(DONE_LINE|WROTE_LINE)='([^']*)'$", sh, flags=re.M))
    assert set(pat) == {"DONE_LINE", "WROTE_LINE"} and "h28" in pat["WROTE_LINE"]
    run = sh[sh.index("  run)"):sh.index("  read)")]
    assert "grep -cvE" in run and not re.search(r"grep (?!-cvE)", run)
    read = sh[sh.index("  read)"):]
    assert '> "$S/read.out" 2>&1' in read and "results/h28_heart_trend_strength" in read
