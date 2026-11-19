"""experiments/harness/h26_heart_scope_forms.py (try H26 of experiments/harness/TRIES.md) on synthetic arrays and
hand-made result tables. Nothing is read from data/ or runs/. Committed result files are read for arithmetic only:
the s2 hedge study's per-seed heart files (their truth and pool records, and the seed spreads of their paired
differences for the noise simulation) and row H25's and the earlier heart runs' committed rows (the guard's dry run
and the noise floor)."""

from __future__ import annotations

import importlib.util
import inspect
import json
import sys
from pathlib import Path

import numpy as np
import pytest

from vec import harness, t2_build
from vec.adoption import SEEDS

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "experiments" / "harness" / "h26_heart_scope_forms.py"
RESULTS = REPO_ROOT / "experiments" / "harness" / "results"


@pytest.fixture(scope="module")
def h26():
    spec = importlib.util.spec_from_file_location("h26_heart_scope_forms_module", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def _cloud_with_core(n=3000, seed=0):
    """A synthetic heart-like stack: a compact core near the origin and a wider shell (x widest, z thinnest)."""
    rng = np.random.default_rng(seed)
    core = rng.random(n) < 0.25
    sd = np.where(core[:, None], [[60.0, 50.0, 30.0]], [[300.0, 220.0, 60.0]])
    C = rng.normal(size=(n, 3)) * sd + np.array([500.0, -200.0, 40.0])
    return C, core


# ------------------------------------------------------------------------------------------------ truths and pools
def test_the_truths_are_the_studys_own_trims_at_its_seed(h26):
    TC, core = _cloud_with_core()
    m = h26.truth_masks(TC, core)
    assert tuple(m) == tuple(h26.payoff.HEART_TRUTHS) == h26.TRUTHS
    assert m["own"].all()
    for i, name in enumerate(h26.payoff.HEART_TRUTHS):
        if name == "own":
            continue
        mode, frac = name.split("_")
        want = h26.payoff.trim_mask(TC, core, float(frac), mode,
                                    np.random.default_rng([h26.payoff.TRUTH_SEED, i]))
        np.testing.assert_array_equal(m[name], want)
        assert m[name][core].all()  # core cells always kept: the five truths share one core
    assert "payoff.trim_mask" in inspect.getsource(h26.truth_masks)


def test_the_pools_use_the_builders_trim_which_equals_the_studys(h26):
    PC, core = _cloud_with_core(seed=3)
    m = h26.pool_masks(PC, core)
    assert m["full_iso_half"].all() and m["s1form_half"].all()
    for arm, frac in (("half_iso_half", 0.5), ("core_iso_half", 0.0)):
        np.testing.assert_array_equal(m[arm], t2_build.radial_trim_mask(PC, core, frac))
        np.testing.assert_array_equal(m[arm], h26.payoff.trim_mask(PC, core, frac, "radial",
                                                                  np.random.default_rng(0)))
        assert m[arm][core].all() and m[arm].sum() < len(PC)
    assert m["core_iso_half"].sum() <= m["half_iso_half"].sum()
    # B's pool is the builder's --entry s2_trim pool (build_test.py build_heart)
    assert "t2_build.radial_trim_mask(C, core, 0.5)" in (REPO_ROOT / "experiments" / "builds" /
                                                         "build_test.py").read_text()


def test_the_records_match_the_studys_layout(h26):
    C, core = _cloud_with_core(seed=4)
    keep = h26.pool_masks(C, core)["half_iso_half"]
    t, p = h26.truth_info(C, core, keep), h26.pool_info(C, core, keep)
    assert set(t) == {"cells", "core_share", "xy_scope"} and set(p) == {"pool_cells", "core_share", "xy_scope"}
    assert t["cells"] == p["pool_cells"] == int(keep.sum()) and t["core_share"] == float(core[keep].mean())


# ------------------------------------------------------------------------------------------------ the oracle size
def test_the_oracle_factors_come_from_the_cores_and_one_serves_every_truth(h26):
    TC, tcore = _cloud_with_core(seed=5)
    PC, pcore = _cloud_with_core(seed=6)
    TC = TC * 1.4
    f = h26.oracle_factors(TC, tcore, PC, pcore)
    assert f["label"] == "oracle size"
    assert f["even"] == pytest.approx(t2_build.rms_radius(TC[tcore]) / t2_build.rms_radius(PC[pcore]), rel=1e-14)
    assert f["in_plane"] == pytest.approx(h26.payoff.xy_rms(TC[tcore]) / h26.payoff.xy_rms(PC[pcore]), rel=1e-14)
    # the study's own line in payoff.run_heart
    assert "f_core = xy_rms(TC[tcore]) / xy_rms(PC[pcore])" in inspect.getsource(h26.payoff.run_heart)
    # every truth keeps the core, so a factor computed on any truth's cells is the same number
    for keep in h26.truth_masks(TC, tcore).values():
        g = h26.oracle_factors(TC[keep], tcore[keep], PC, pcore)
        assert g["even"] == f["even"] and g["in_plane"] == f["in_plane"]


def test_the_even_rescale_is_isotropic_about_the_centroid(h26):
    C, _ = _cloud_with_core(seed=7)
    out = h26.even_rescale(C, 1.37)
    mu = C.mean(0)
    np.testing.assert_allclose(out, mu + (C - mu) * 1.37, rtol=0, atol=1e-9)
    assert t2_build.rms_radius(out) == pytest.approx(1.37 * t2_build.rms_radius(C), rel=1e-12)
    assert "t2_build.rescale_to_rms" in inspect.getsource(h26.even_rescale)


def test_the_in_plane_rescale_holds_z_as_the_study_did(h26):
    C, _ = _cloud_with_core(seed=8)
    out = h26.in_plane_rescale(C, 1.37)
    np.testing.assert_array_equal(out[:, 2], C[:, 2])
    c = C[:, :2].mean(0)
    np.testing.assert_allclose(out[:, :2], c + (C[:, :2] - c) * 1.37, rtol=0, atol=1e-9)


def test_the_half_is_row_h25s(h26):
    X = np.arange(1001 * 3, dtype=float).reshape(1001, 3)
    C, _ = _cloud_with_core(n=1001, seed=9)
    f = {"even": 1.2, "in_plane": 1.3}
    for s in (0, 4, 9):
        rows = h26.h25.kept_rows(1001, s, 1)
        Xa, Ca = h26.form(X, C, "full_iso_half", f, s)
        np.testing.assert_array_equal(Xa, X[rows])
        np.testing.assert_allclose(Ca, h26.even_rescale(C, 1.2)[rows], atol=1e-12)
        Xs, Cs = h26.form(X, C, "s1form_half", f, s)
        np.testing.assert_array_equal(Xs, X[rows])
        np.testing.assert_array_equal(Cs[:, 2], C[rows, 2])
    assert h26.h25.HALF_SEED_BASE == 2 * harness.BUILD_SEED_OFFSET


# ------------------------------------------------------------------------------------------------ the transitions
def test_the_transitions_folds_and_expression_steps(h26):
    ts = {k: h26.payoff.heart_task(k) for k in h26.TASK_KEYS}
    assert h26.TASK_KEYS[0] == h26.payoff.HEART_TASK and h26.TASK_KEYS[1] in h26.payoff.HEART_EXTRA_TASKS
    assert {ts[k].target for k in h26.FOLDS[0]} == {"T2/E10.5.h5ad"}  # two transitions, one truth: one fold
    assert {ts[k].target for k in h26.FOLDS[1]} == {"T2/E8.5.h5ad"}
    assert sorted(k for f in h26.FOLDS for k in f) == sorted(h26.TASK_KEYS)
    harness_keys = {t.key for t in harness.load_tasks()["T2:heart:test_extrap"]}
    assert set(h26.HARNESS_KEYS) <= harness_keys and h26.TASK_KEYS[1] not in harness_keys
    for k in h26.HARNESS_KEYS:
        assert h26.s1_expression(k) == "b0_trend"
        harness._b0_trend_task_check(harness.TaskView.of(ts[k]))  # the trend applies there
        assert "ref" in ts[k].inputs and ts[k].inputs["ref"] == ts[k].reference
    assert h26.s1_expression(h26.TASK_KEYS[1]) == "copy" and set(ts[h26.TASK_KEYS[1]].inputs) == {"ref"}
    assert set(h26.STUDY_RUNS) == set(h26.TASK_KEYS)


# ------------------------------------------------------------------------------------------------ end to end, synthetic
@pytest.fixture(scope="module")
def synth():
    """The synthetic heart stages, board and context of tests/test_b0_trend.py (reused, not copied), with a core
    column (the CM lineage)."""
    spec = importlib.util.spec_from_file_location("b0_trend_tests_for_h26", REPO_ROOT / "tests" / "test_b0_trend.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    inputs = {"ref": mod._stage(1600, seed=2, shift=0.0, tiny=60),
              "start": mod._stage(1100, seed=3, shift=0.8, tiny=10),
              "earlier": mod._stage(900, seed=5, shift=1.6, tiny=10)}
    for a in inputs.values():
        a.obs["core"] = np.asarray(a.obs["merged"]).astype(str) == "CM"
    return mod, inputs


def test_full_iso_half_has_row_h25s_iso_half_cells_and_expression(h26, synth):
    """The identity behind the guard: on one context and seed, full_iso_half and row H25's iso_half are the same draw,
    the same half and the same expression; their positions differ by one even factor about the centroid."""
    mod, inputs = synth
    ref = np.asarray(inputs["ref"].obsm["spatial_3D"], dtype=float)
    size = h26.h25.oracle_size_result(1.6 * t2_build.rms_radius(ref), "E9.5", 9.5, t2_build.rms_radius(ref),
                                      float(ref[:, 2].std()), 10.5)
    f = {"even": 1.3, "in_plane": 1.4}
    for seed in (harness.BUILD_SEED_OFFSET + 1, harness.BUILD_SEED_OFFSET + 7):
        ctx = mod._ctx(inputs, seed=seed)
        theirs = h26.h25.make_arm("iso_half", size, {}).build(ctx)
        mine = h26.build_arm("full_iso_half", ctx, f, "b0_trend")
        np.testing.assert_array_equal(mine.X, theirs.X)
        a = mine.coords.astype(float) - mine.coords.astype(float).mean(0)
        b = theirs.coords.astype(float) - theirs.coords.astype(float).mean(0)
        k = np.sqrt((a ** 2).sum() / (b ** 2).sum())
        np.testing.assert_allclose(a, b * k, rtol=0, atol=1e-3 * np.abs(a).max())


def test_every_arm_builds_through_the_harness_and_the_checked_builder(h26, synth):
    mod, inputs = synth
    pool = inputs["ref"]
    PC = np.asarray(pool.obsm["spatial_3D"], dtype=float)
    pm = h26.pool_masks(PC, np.asarray(pool.obs["core"], bool))
    f = {"even": 1.3, "in_plane": 1.4}
    n = mod.BOARD.max_cells
    for seed in (harness.BUILD_SEED_OFFSET + 2,):
        st = {}
        for arm in h26.BUILT:
            p = pool if h26.TRIM[arm] is None else pool[pm[arm]].copy()
            ctx = mod._ctx({**inputs, "ref": p}, seed=seed)
            st[arm] = h26.build_arm(arm, ctx, f, "b0_trend")
        assert {a: s.n for a, s in st.items()} == {a: n // 2 for a in h26.BUILT}  # the harness's count, then the half
        s = seed - harness.BUILD_SEED_OFFSET
        rows = h26.h25.kept_rows(n, s, 1)
        trend = harness._b0_trend(mod._ctx(inputs, seed=seed))
        np.testing.assert_array_equal(st["s1form_half"].X, trend.X[rows])
        np.testing.assert_allclose(st["s1form_half"].coords[:, 2], trend.coords[rows, 2], atol=1e-3)
        copy = harness._copy(mod._ctx(inputs, seed=seed))
        np.testing.assert_array_equal(st["full_iso_half"].X, copy.X[rows])
        assert not np.array_equal(st["half_iso_half"].X, st["full_iso_half"].X)  # another pool, another draw
        # the even arms keep their pool's shape: centred positions are the copy's times one factor
        for arm in ("half_iso_half", "core_iso_half"):
            p = pool[pm[arm]].copy()
            c = harness._copy(mod._ctx({**inputs, "ref": p}, seed=seed)).coords.astype(float)[rows]
            a = st[arm].coords.astype(float)
            a, c = a - a.mean(0), c - c.mean(0)
            np.testing.assert_allclose(a, c * 1.3, rtol=0, atol=1e-2)


# ------------------------------------------------------------------------------------------------ hand-made tables
def _skills(**over):
    s = {m: 0.5 for m in ("de_score", "de_direction", "mmd_u", "variogram", "d2_shape", "occupancy_dice",
                          "scale_log_ratio", "neighborhood_mmd")}
    return {**s, **over}


def _groups(v=0.125):
    return {g: v for g in ("expression_change", "cell_state", "shape_growth", "local_organisation")}


def _part(h26, key, seeds, score=None, commit="abc123", clean=True, info=None):
    """A hand-made part. score: {truth: {arm: value or per-seed list (indexed by seed)}}; default 50."""
    score = score or {}
    rows = []
    for s in seeds:
        sc, sk, gr = {}, {}, {}
        for q in h26.TRUTHS:
            sc[q] = {}
            for a in h26.ALL:
                v = score.get(q, {}).get(a, 50.0)
                sc[q][a] = v[s] if isinstance(v, list) else v
            sk[q] = {a: _skills() for a in h26.ALL}
            gr[q] = {a: _groups() for a in h26.ALL}
        rows.append({"seed": s, "scores": sc, "skills": sk, "groups": gr})
    info = info or _info(h26)
    return {"axis": h26.AXIS, "task": key, "rows": rows, "provenance": {"commit": commit, "tree_clean": clean},
            **info}


def _info(h26):
    return {"truth_info": {q: {"cells": 100 + i, "core_share": 0.2 + 0.1 * i, "xy_scope": 0.1}
                           for i, q in enumerate(h26.TRUTHS)},
            "pred_info": {"full_iso_half": {"pool_cells": 900, "core_share": 0.23, "xy_scope": 0.2},
                          "half_iso_half": {"pool_cells": 600, "core_share": 0.3, "xy_scope": 0.1},
                          "core_iso_half": {"pool_cells": 300, "core_share": 1.0, "xy_scope": 0.0},
                          "s1form_half": {"pool_cells": 900, "core_share": 0.23, "xy_scope": 0.2}},
            "oracle_factors": {"label": "oracle size", "even": 1.1, "in_plane": 1.2}}


def _study(h26, score=None):
    """A hand-made study: the same truth and pool records (pred_info under the study's names) and rows."""
    i = _info(h26)
    pred = {"hold": i["pred_info"]["full_iso_half"], "half_radial": i["pred_info"]["half_iso_half"],
            "core_radial": i["pred_info"]["core_iso_half"]}
    out = {}
    for key in h26.TASK_KEYS:
        rows = []
        for s in SEEDS:
            sc = {q: {p: 50.0 for p in ("floor", "hold", "half_size", "half_radial", "half_thin", "core_radial")}
                  for q in h26.TRUTHS}
            for q, d in (score or {}).items():
                for p, v in d.items():
                    sc[q][p] = v[s] if isinstance(v, list) else v
            rows.append({"seed": s, "scores": sc})
        out[key] = {"rows": rows, "truth_info": i["truth_info"], "pred_info": pred,
                    "oracle_core_xy_factor": i["oracle_factors"]["in_plane"]}
    return out


def _h25(h26):
    return {k: {"rows": [{"seed": s, "skills": {"iso_half": _skills()}} for s in SEEDS]} for k in h26.HARNESS_KEYS}


def _parts(h26, score=None, **kw):
    """Five parts of two seeds per transition, as the run is launched."""
    return {k: [_part(h26, k, (2 * j, 2 * j + 1), score=(score or {}).get(k), **kw) for j in range(5)]
            for k in h26.TASK_KEYS}


def test_a_clean_run_is_read_and_a_tie_picks_the_file_in_force(h26):
    res = h26.reading(_parts(h26), _study(h26), _h25(h26), commit="abc")
    assert not res["void"], res["why"]
    assert res["pick"]["D3"] == 0.0 and res["pick"]["DJ4"] == 0.0 and res["pick"]["pick"].startswith("A")
    assert res["R1_carry_over"]["agreement"]["half_iso_half"] == {"agree": 15, "of": 15}
    assert "oracle size" in h26.text(res) and "PICK RULE" in h26.text(res)


def test_the_reading_is_void_on_missing_repeated_or_mixed_parts(h26):
    st, hh = _study(h26), _h25(h26)
    p = _parts(h26)
    del p[h26.TASK_KEYS[1]]
    assert h26.reading(p, st, hh)["void"]
    p = _parts(h26)
    p[h26.TASK_KEYS[0]] = p[h26.TASK_KEYS[0]][:4]  # seeds 8 and 9 missing
    assert h26.reading(p, st, hh)["void"]
    p = _parts(h26)
    p[h26.TASK_KEYS[0]].append(_part(h26, h26.TASK_KEYS[0], (9,)))  # seed 9 twice
    res = h26.reading(p, st, hh)
    assert res["void"] and any("more than one part" in w for w in res["why"])
    p = _parts(h26)
    del p[h26.TASK_KEYS[2]][1]["rows"][0]["scores"]["thin_0"]["core_iso_half"]
    assert h26.reading(p, st, hh)["void"]
    p = _parts(h26)
    del p[h26.TASK_KEYS[2]][1]["rows"][0]["groups"]["radial_0"]
    assert h26.reading(p, st, hh)["void"]
    p = _parts(h26)
    p[h26.TASK_KEYS[1]][3]["provenance"]["commit"] = "def456"
    assert h26.reading(p, st, hh)["void"]
    assert h26.reading(_parts(h26, clean=False), st, hh)["void"]
    res = h26.reading(_parts(h26), st, hh, commit="fff")
    assert res["void"] and any("row's commit" in w for w in res["why"])
    assert h26.reading(_parts(h26), st, hh, commit="")["void"]
    assert res["pick"].startswith("A") and "void" in res["follows"]


def test_the_reading_is_void_when_a_truth_a_pool_or_the_factor_differs_from_the_study(h26):
    hh = _h25(h26)
    st = _study(h26)
    st[h26.TASK_KEYS[0]]["truth_info"] = {**st[h26.TASK_KEYS[0]]["truth_info"],
                                          "radial_0": {"cells": 103, "core_share": 0.5, "xy_scope": 0.1}}
    res = h26.reading(_parts(h26), st, hh)
    assert res["void"] and any("truth radial_0 differs" in w for w in res["why"])
    st = _study(h26)
    st[h26.TASK_KEYS[1]]["pred_info"] = {**st[h26.TASK_KEYS[1]]["pred_info"],
                                         "half_radial": {"pool_cells": 601, "core_share": 0.3, "xy_scope": 0.1}}
    res = h26.reading(_parts(h26), st, hh)
    assert res["void"] and any("half_iso_half's pool differs" in w for w in res["why"])
    st = _study(h26)
    st[h26.TASK_KEYS[2]]["pred_info"] = {**st[h26.TASK_KEYS[2]]["pred_info"],
                                         "core_radial": {"pool_cells": 300, "core_share": 0.99, "xy_scope": 0.0}}
    assert h26.reading(_parts(h26), st, hh)["void"]
    st = _study(h26)
    st[h26.TASK_KEYS[2]]["oracle_core_xy_factor"] = 1.2 * (1 + 1e-9)
    res = h26.reading(_parts(h26), st, hh)
    assert res["void"] and any("in-plane oracle factor" in w for w in res["why"])
    st = _study(h26)
    st[h26.TASK_KEYS[2]]["oracle_core_xy_factor"] = 1.2 * (1 + 1e-13)  # float noise passes
    assert not h26.reading(_parts(h26), st, hh)["void"]


def test_the_reading_is_void_when_full_iso_half_is_not_row_h25s_iso_half(h26):
    st = _study(h26)
    hh = _h25(h26)
    hh[h26.HARNESS_KEYS[1]]["rows"][6]["skills"]["iso_half"] = _skills(mmd_u=0.5 + 2e-9)
    res = h26.reading(_parts(h26), st, hh)
    assert res["void"] and any("full_iso_half mmd_u differs from row H25's iso_half" in w for w in res["why"])
    hh = _h25(h26)
    hh[h26.HARNESS_KEYS[0]]["rows"][0]["skills"]["iso_half"] = _skills(de_score=0.5 + 5e-10)  # within 1e-9
    assert not h26.reading(_parts(h26), st, hh)["void"]
    hh = _h25(h26)
    hh[h26.HARNESS_KEYS[0]]["rows"][0]["skills"]["iso_half"] = _skills(d2_shape=0.9)  # size-dependent: not checked
    assert not h26.reading(_parts(h26), st, hh)["void"]
    hh = _h25(h26)
    del hh[h26.HARNESS_KEYS[0]]
    assert h26.reading(_parts(h26), st, hh)["void"]
    # the identity is checked on the truth as released only, and only on the two harness backtests
    assert h26.OWN == "own" and h26.TASK_KEYS[1] not in h26.HARNESS_KEYS


# ------------------------------------------------------------------------------------------------ the pick rule
def _const(h26, gain_by_key: dict[str, float], s1: float = 0.0) -> dict:
    """Parts where half_iso_half beats full_iso_half by a constant per transition on every truth, with s1form_half
    at 50 + s1 (below the arms when s1 < 0, so V is the arm's own score)."""
    return {k: {q: {"half_iso_half": 50.0 + g, "s1form_half": 50.0 + s1} for q in h26.TRUTHS}
            for k, g in gain_by_key.items()}


def test_the_pick_needs_both_d3_and_dj4_at_the_bar(h26):
    st, hh = _study(h26), _h25(h26)
    k0, k1, k2 = h26.TASK_KEYS  # E9.5->E10.5, E8.25_late->E8.5, E8.75->E10.5
    res = h26.reading(_parts(h26, score=_const(h26, {k0: 1.0, k1: 1.0, k2: 1.0}, s1=-10)), st, hh)
    assert res["pick"]["D3"] == pytest.approx(1.0) and res["pick"]["DJ4"] == pytest.approx(1.0)
    assert res["pick"]["pick"].startswith("B")
    res = h26.reading(_parts(h26, score=_const(h26, {k0: 0.99, k1: 0.99, k2: 0.99}, s1=-10)), st, hh)
    assert res["pick"]["pick"].startswith("A")
    # gains on the two E10.5 transitions only: D3 = 1.067 passes, DJ4 = 0.8 does not (one truth counts once)
    res = h26.reading(_parts(h26, score=_const(h26, {k0: 1.6, k1: 0.0, k2: 1.6}, s1=-10)), st, hh)
    assert res["pick"]["D3"] == pytest.approx(3.2 / 3) and res["pick"]["DJ4"] == pytest.approx(0.8)
    assert res["pick"]["pick"].startswith("A")
    # a gain on E8.5 alone: DJ4 = 1.0 passes, D3 = 0.667 does not
    res = h26.reading(_parts(h26, score=_const(h26, {k0: 0.0, k1: 2.0, k2: 0.0}, s1=-10)), st, hh)
    assert res["pick"]["DJ4"] == pytest.approx(1.0) and res["pick"]["D3"] == pytest.approx(2 / 3)
    assert res["pick"]["pick"].startswith("A")
    # the better of two: where s1form_half is above both arms, neither trim adds anything
    res = h26.reading(_parts(h26, score=_const(h26, {k0: 3.0, k1: 3.0, k2: 3.0}, s1=+5)), st, hh)
    assert res["pick"]["D3"] == 0.0 and res["pick"]["pick"].startswith("A")
    # a void run never picks B
    res = h26.reading(_parts(h26, score=_const(h26, {k0: 5.0, k1: 5.0, k2: 5.0}, s1=-10), clean=False), st, hh)
    assert res["void"] and res["pick"].startswith("A")


def test_c_is_reported_and_never_picked(h26):
    st, hh = _study(h26), _h25(h26)
    score = {k: {q: {"core_iso_half": 60.0, "s1form_half": 40.0} for q in h26.TRUTHS} for k in h26.TASK_KEYS}
    res = h26.reading(_parts(h26, score=score), st, hh)
    assert res["pick"]["pick"].startswith("A")
    assert res["R2_expected_V"][h26.TASK_KEYS[0]]["uniform"]["core_iso_half"] == pytest.approx(60.0)


def test_the_reported_statistics_decide_nothing(h26):
    """Plain-mean, per-prior and minimum-over-prior D3 and DJ4 are reported beside the pick and never move it."""
    st, hh = _study(h26), _h25(h26)
    # B ahead by 3 on every truth but the released one, where it is 2 behind; s1form_half far below
    score = {k: {q: {"half_iso_half": 50.0 + (-2.0 if q == "own" else 3.0), "s1form_half": 30.0}
                 for q in h26.TRUTHS} for k in h26.TASK_KEYS}
    res = h26.reading(_parts(h26, score=score), st, hh)
    p, rep = res["pick"], res["pick"]["reported_not_decided"]
    assert rep["V_per_prior"]["released"] == {"D3": pytest.approx(-2.0), "DJ4": pytest.approx(-2.0)}
    assert rep["V_per_prior"]["thin_family"]["D3"] == pytest.approx(3.0)
    assert rep["V_min_over_priors"]["D3"] == pytest.approx(-2.0)
    assert rep["plain_mean"]["D3"] == pytest.approx(p["D3"])  # s1 far below: V is the arm's own score
    assert p["D3"] >= 1.0 and p["pick"].startswith("B")  # the fixed rule picks B: the reported -2.0 does not decide
    assert "reported, not decided on" in h26.text(res) and "released D3 -2.000" in h26.text(res)


def test_the_carry_over_reference_figures_quoted_in_row_h26(h26):
    """Row H26 quotes, as the committed reference for its outcome, what the pick statistics would be if the study's
    own half_radial - hold carried over exactly (committed per-seed files, arithmetic only): on the plain mean and as
    V with s1 = hold, over the five priors and on the released target alone."""
    st = h26.load_study()
    pm, pv = {}, {}
    for key in h26.TASK_KEYS:
        rows = st[key]["rows"]
        for name, pr in h26.PRIORS.items():
            m = v = 0.0
            for q, w in pr.items():
                h = np.array([r["scores"][q]["hold"] for r in rows])
                b = np.array([r["scores"][q]["half_radial"] for r in rows])
                m += w * (b - h).mean()
                v += w * (np.maximum(h, b).mean() - h.mean())
            pm.setdefault(key, {})[name], pv.setdefault(key, {})[name] = m, v
    got = {"plain": h26.d3_dj4(pm)[:2], "V": h26.d3_dj4(pv)[:2],
           "plain released": h26.d3_dj4(pm, ("released",))[:2], "V released": h26.d3_dj4(pv, ("released",))[:2]}
    print("\nH26 carry-over reference (study half_radial - hold): "
          + "; ".join(f"{k} D3 {a:+.2f} DJ4 {b:+.2f}" for k, (a, b) in got.items()))
    want = {"plain": (3.08, 3.39), "V": (4.11, 4.19), "plain released": (-0.85, -0.20), "V released": (1.20, 1.41)}
    for k, (a, b) in want.items():
        assert got[k][0] == pytest.approx(a, abs=0.006) and got[k][1] == pytest.approx(b, abs=0.006), k


def test_the_priors_are_the_frozen_sections_five(h26):
    P = h26.PRIORS
    assert list(P) == ["uniform", "own_third", "released", "radial_family", "thin_family"]
    assert all(sum(w.values()) == pytest.approx(1.0) for w in P.values())
    assert P["own_third"]["own"] == pytest.approx(1 / 3) and P["own_third"]["thin_0"] == pytest.approx(1 / 6)
    assert set(P["radial_family"]) == {"own", "radial_0.5", "radial_0"} and set(P["thin_family"]) == {"thin_0.5",
                                                                                                        "thin_0"}
    assert set(P["released"]) == {"own"}


def test_the_expected_better_of_two_is_the_studys_hedge_arithmetic(h26, tmp_path):
    """E V under each prior equals payoff.hedge's (the study's own better-of-two) on the same hand-made table."""
    rng = np.random.default_rng(3)
    key = h26.TASK_KEYS[0]
    score = {q: {a: list(50 + rng.normal(0, 3, 10)) for a in h26.BUILT} for q in h26.TRUTHS}
    part = _part(h26, key, SEEDS, score=score)
    for j in range(5):  # payoff.hedge reads directories of per-seed parts
        d = tmp_path / f"{2 * j}-{2 * j + 1}"
        d.mkdir()
        sub = dict(part, rows=[r for r in part["rows"] if r["seed"] in (2 * j, 2 * j + 1)],
                   pred_values=list(h26.BUILT), truth_values=list(h26.TRUTHS))
        (d / "h.json").write_text(json.dumps(sub))
    extra = {k: list(v) for k, v in h26.PRIORS.items() if k not in ("uniform", "own_third")}
    hedge = h26.payoff.hedge(sorted(tmp_path.iterdir()), "s1form_half", None, extra)
    vals = {key: {q: {"V": {a: h26.better_of_two(part["rows"], q, a) for a in h26.BUILT}} for q in h26.TRUTHS}}
    ev = {name: {a: sum(w * vals[key][q]["V"][a] for q, w in pr.items()) for a in h26.BUILT}
          for name, pr in h26.PRIORS.items()}
    for name in h26.PRIORS:
        for a in h26.BUILT:
            assert ev[name][a] == pytest.approx(hedge["candidates"][a][name]["mean"], abs=1e-9)


def test_r1_agreement(h26):
    assert h26.agrees(2.0, 0.1) is True
    assert h26.agrees(-2.0, 0.5) is False  # signs differ and -2.0 is not near zero
    assert h26.agrees(-2.0, -0.5) is True
    assert h26.agrees(0.9, -0.9) is True  # both within 1.0 of zero
    assert h26.agrees(1.1, -0.2) is False


def test_r1_reads_the_study_against_hold(h26):
    st = _study(h26, score={"thin_0": {"half_radial": 55.0, "core_radial": 45.0}})
    score = {k: {"thin_0": {"half_iso_half": 53.0, "core_iso_half": 53.0}} for k in h26.TASK_KEYS}
    res = h26.reading(_parts(h26, score=score), st, _h25(h26))
    c = res["R1_carry_over"]["cells"][h26.TASK_KEYS[0]]["thin_0"]
    assert c["half_iso_half"]["new"]["mean"] == pytest.approx(3.0) and c["half_iso_half"]["study"]["mean"] == 5.0
    assert c["half_iso_half"]["agree"] and not c["core_iso_half"]["agree"]
    assert res["R1_carry_over"]["agreement"]["core_iso_half"] == {"agree": 12, "of": 15}


# ------------------------------------------------------------------------------------------------ noise, guards, J5
def _study_spreads(h26, how=max) -> dict[str, float]:
    """Per transition, the seed SD of the study's committed paired differences of a trim minus hold (the nearest
    committed analogue of a trimmed arm against the full one): `max` the largest over the five truths and both trims
    (half_radial, core_radial); `np.median` the median over the five truths of half_radial - hold (B's own analogue)."""
    st = h26.load_study()
    out = {}
    for key in h26.TASK_KEYS:
        rows = st[key]["rows"]
        trims = ("half_radial", "core_radial") if how is max else ("half_radial",)
        out[key] = float(how([float(np.std([r["scores"][q][p] - r["scores"][q]["hold"] for r in rows], ddof=1))
                              for q in h26.TRUTHS for p in trims]))
    return out


def _h25_stand_in(h26) -> tuple[dict[str, float], dict[str, float]]:
    """From row H25's committed per-seed scores: s1form_half minus iso_half (the stand-in against the file in force's
    form, the same pool, build seed and half rows), its mean (the level) and its seed SD, per E10.5 backtest. E8.25_late
    -> E8.5 has no H25 run: there the stand-in is the same copy as A with an in-plane instead of an even stretch, so it
    takes level 0 and the SMALLER of the two committed SDs (the closer copy, the larger premium)."""
    h25 = h26.load_h25()
    level, sd = {}, {}
    for key, r in h25.items():
        d = np.array([x["scores"]["s1form_half"] - x["scores"]["iso_half"] for x in r["rows"]])
        level[key], sd[key] = float(d.mean()), float(d.std(ddof=1))
    k1 = h26.TASK_KEYS[1]
    level[k1], sd[k1] = 0.0, min(sd.values())
    return level, sd


def _null_draws(h26, sd: dict[str, float], s1_sd: dict[str, float], s1_level: dict[str, float], draws: int,
                seed: int, share_truths: bool, share_e105: bool) -> dict[str, np.ndarray]:
    """Pure noise: half_iso_half (B) and full_iso_half (A) the same in distribution. Per transition, truth and seed,
    A and B are independent draws with SD sd / sqrt(2), so B - A has the committed paired spread sd; s1form_half is
    A's draw plus s1_level plus its own small term of SD s1_sd (it shares A's pool, build seed and half rows; H25's
    committed s1form_half and iso_half correlate 0.75 and 0.96 per seed). A component common to all three arms
    cancels in every difference, so it is left out. share_truths: one draw for all five truths of a transition (they
    are trims of one stage); share_e105: one draw for the two E10.5 transitions. Returns, per draw, the per-transition,
    per-prior differences B - A of E V and of the expected plain mean (arrays draws x 3 x 5)."""
    rng = np.random.default_rng(seed)
    W = np.array([[pr.get(q, 0.0) for q in h26.TRUTHS] for pr in h26.PRIORS.values()])  # priors x truths
    nt = 1 if share_truths else len(h26.TRUTHS)
    cache, dV, dM = {}, [], []
    for key in h26.TASK_KEYS:
        src = h26.TASK_KEYS[0] if (share_e105 and key == h26.TASK_KEYS[2]) else key
        if src not in cache:
            cache[src] = rng.normal(size=(3, draws, nt, len(SEEDS)))
        z = cache[src]
        A = sd[key] / np.sqrt(2) * z[0]
        B = sd[key] / np.sqrt(2) * z[1]
        S = A + s1_level[key] + s1_sd[key] * z[2]
        A, B, S = (np.broadcast_to(x, (draws, len(h26.TRUTHS), len(SEEDS))) for x in (A, B, S))
        dV.append((np.maximum(S, B).mean(-1) - np.maximum(S, A).mean(-1)) @ W.T)
        dM.append((B.mean(-1) - A.mean(-1)) @ W.T)
    return {"V": np.stack(dV, 1), "mean": np.stack(dM, 1)}


def _stats(per: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """D3 and DJ4 of per-draw arrays draws x transitions (TASK_KEYS order) x priors."""
    d3 = per.mean(axis=(1, 2))
    dj4 = 0.5 * (0.5 * (per[:, 0] + per[:, 2]).mean(1) + per[:, 1].mean(1))
    return d3, dj4


def test_the_vectorised_null_is_the_scripts_pick(h26):
    """The noise simulation's D3 and DJ4 equal h26.pick's on the same draws: its arithmetic is the rule's."""
    level, s1sd = _h25_stand_in(h26)
    sd = _study_spreads(h26)
    sim = _null_draws(h26, sd, s1sd, level, draws=5, seed=7, share_truths=False, share_e105=False)
    d3, dj4 = _stats(sim["V"])
    m3, mj4 = _stats(sim["mean"])
    for i in range(5):
        ev = {k: {n: {"full_iso_half": 0.0, "half_iso_half": float(sim["V"][i, j, p])}
                  for p, n in enumerate(h26.PRIORS)} for j, k in enumerate(h26.TASK_KEYS)}
        em = {k: {n: {"full_iso_half": 0.0, "half_iso_half": float(sim["mean"][i, j, p])}
                  for p, n in enumerate(h26.PRIORS)} for j, k in enumerate(h26.TASK_KEYS)}
        res = h26.pick(ev, True, em)
        assert res["D3"] == pytest.approx(d3[i], abs=1e-12) and res["DJ4"] == pytest.approx(dj4[i], abs=1e-12)
        rep = res["reported_not_decided"]
        assert rep["plain_mean"]["D3"] == pytest.approx(m3[i], abs=1e-12)
        assert rep["plain_mean"]["DJ4"] == pytest.approx(mj4[i], abs=1e-12)


def test_the_pick_rule_on_pure_noise(h26):
    """The fixed rule on pure noise (B and A the same in distribution), with s1form_half sharing A's draw as it does in
    the run (review of 2026-11-19, blocker 1). Spreads: the study's committed trim-minus-hold seed SDs, the largest
    over truths and trims and the median over truths of half_radial; the stand-in's level and SD from row H25's
    committed s1form_half - iso_half (and level, and far below, where the better of two is the arm's own score).
    The rates of three other rules the review proposed are printed for the main session; the pick uses none of them.
    Printed for row H26."""
    level_h25, s1sd = _h25_stand_in(h26)
    levels = {"H25's levels": level_h25, "level": {k: 0.0 for k in h26.TASK_KEYS},
              "far below": {k: -100.0 for k in h26.TASK_KEYS}}
    spreads = {"largest": _study_spreads(h26, max), "median": _study_spreads(h26, np.median)}
    rate = {}
    print("\nH26 null spreads (sd of paired differences, board points): " + "; ".join(
        f"{sn} " + ", ".join(f"{k} {v:.3f}" for k, v in sd.items()) for sn, sd in spreads.items()))
    print("H26 stand-in s1form_half - iso_half (row H25): " + ", ".join(
        f"{k} level {level_h25[k]:+.3f} sd {s1sd[k]:.3f}" for k in h26.TASK_KEYS))
    for sn, sd in spreads.items():
        for ln, lv in levels.items():
            for share, kw in (("independent truths", dict(share_truths=False, share_e105=False)),
                              ("shared truths", dict(share_truths=True, share_e105=False)),
                              ("shared truths and E10.5", dict(share_truths=True, share_e105=True))):
                sim = _null_draws(h26, sd, s1sd, lv, draws=4000, seed=20261119, **kw)
                d3, dj4 = _stats(sim["V"])
                m3, mj4 = _stats(sim["mean"])
                r3, rj4 = _stats(sim["V"][:, :, 2:3])  # the released target alone
                fixed = (d3 >= h26.PICK_BAR) & (dj4 >= h26.PICK_BAR)
                rate[(sn, ln, share)] = r = {
                    "null D3": float(d3.mean()), "fixed rule": float(fixed.mean()),
                    "+ plain mean >= 0": float((fixed & (m3 >= 0) & (mj4 >= 0)).mean()),
                    "+ plain mean >= 1": float((fixed & (m3 >= 1) & (mj4 >= 1)).mean()),
                    "+ released alone >= 1": float((fixed & (r3 >= 1) & (rj4 >= 1)).mean())}
                print(f"H26 null | spreads {sn} | s1 {ln} | {share}: " + ", ".join(
                    f"{k} {v:+.3f}" if k == "null D3" else f"{k} {v:.2%}" for k, v in r.items()))
    worst = max(r["fixed rule"] for r in rate.values())
    print(f"H26 the fixed rule's worst false pick of B on pure noise: {worst:.2%}")
    # with no premium (s1 far below), the rule holds its noise level
    assert all(r["fixed rule"] < 0.05 for (sn, ln, sh), r in rate.items() if ln == "far below")
    assert rate[("largest", "far below", "independent truths")]["fixed rule"] < 0.005
    # the premium: under pure noise E V favours B whenever s1form_half is within reach of the arms
    assert all(r["null D3"] > 0.2 for (sn, ln, sh), r in rate.items() if ln != "far below")
    assert all(abs(r["null D3"]) < 0.05 for (sn, ln, sh), r in rate.items() if ln == "far below")
    # and at the study's largest spreads it takes the fixed rule's false pick of B well above 5% (review, blocker 1)
    assert rate[("largest", "level", "shared truths and E10.5")]["fixed rule"] > 0.2
    # a plain-mean condition at the bar removes the premium: the rate falls to the no-premium level
    for (sn, _ln, sh), r in rate.items():
        assert r["+ plain mean >= 1"] <= rate[(sn, "far below", sh)]["fixed rule"] + 0.01


def test_the_study_parts_agree_on_every_record_the_guards_check(h26):
    """The guards' dry run (rule 9, J7) on committed reproduction data other than the pair the row reads: each
    transition's five study parts were five processes that rebuilt the truths, the pools and the oracle factor
    independently; every record the void conditions check agrees exactly across them."""
    st = h26.load_study()
    for key in h26.TASK_KEYS:
        parts = st[key]["parts"]
        assert len(parts) == 5 and st[key]["task"] == key
        assert sorted(r["seed"] for r in st[key]["rows"]) == list(SEEDS)
        for p in parts[1:]:
            assert p["truth_info"] == parts[0]["truth_info"]
            assert p["pred_info"] == parts[0]["pred_info"]
            assert p["oracle_core_xy_factor"] == parts[0]["oracle_core_xy_factor"]
            assert p["truth_seed"] == h26.payoff.TRUTH_SEED
        assert set(parts[0]["truth_info"]) == set(h26.TRUTHS)
        assert {"hold", "half_radial", "core_radial"} <= set(parts[0]["pred_info"])
        # the pools' records equal the truth-free part of the hold pool: the full pool is hold's
        assert parts[0]["pred_info"]["hold"]["pool_cells"] >= parts[0]["pred_info"]["half_radial"]["pool_cells"]


def test_the_expression_identity_on_the_committed_heart_runs(h26):
    """The guard's dry run (J7) for the H25 identity, on committed data other than the pair the row reads: (a) in row
    H25's own run, iso_half and copy_half (same draw, half and expression; another size) agree on the three
    expression-only skills to 1e-9 in every seed; (b) across every committed run of the two backtests, the copy's
    expression-only skills equal row H11's in every seed to 1e-9, so these skills reproduce across commits."""
    h25 = h26.load_h25()
    assert set(h25) == set(h26.HARNESS_KEYS)
    for r in h25.values():
        for row in r["rows"]:
            for m in h26.EXPRESSION_ONLY:
                assert abs(row["skills"]["iso_half"][m] - row["skills"]["copy_half"][m]) <= h26.TOL
    ref = h26.h25.load_dir(REPO_ROOT / h26.h25.GUARD_RUN)
    seen = 0
    for p in sorted(RESULTS.glob("*/tasks/T2h_*.json")):
        r = json.loads(p.read_text())
        key = r.get("task", {}).get("key")
        if key not in h26.HARNESS_KEYS or sorted(x["seed"] for x in r["rows"]) != list(SEEDS):
            continue
        if not all("copy" in x["skills"] for x in r["rows"]):
            continue
        mine = {x["seed"]: x["skills"]["copy"] for x in r["rows"]}
        for x in ref[key]["rows"]:
            for m in h26.EXPRESSION_ONLY:
                assert abs(mine[x["seed"]][m] - x["skills"]["copy"][m]) <= h26.TOL, (p, m)
        seen += 1
    print(f"\nH26 guard dry run: {seen} committed runs of the two backtests reproduce the copy's expression skills")
    assert seen >= 30


def test_the_noise_floor_of_rule_9_j5(h26):
    """J5's seed spread from committed results of the entry in force against the baseline on the counting tasks that
    have one: row H25's iso_half minus the copy on the two E10.5 backtests (one truth). Printed for row H26."""
    h25 = h26.load_h25()
    se = {}
    for key, r in h25.items():
        d = np.array([x["scores"]["iso_half"] - x["scores"]["copy"] for x in r["rows"]])
        se[key] = float(d.std(ddof=1) / np.sqrt(len(d)))
    sd = _study_spreads(h26)
    print("\nH26 J5 seed spread (SE of iso_half - copy): " + ", ".join(f"{k} {v:.3f}" for k, v in se.items())
          + "; the study's trimmed-minus-hold SE, largest over truths: "
          + ", ".join(f"{k} {v / np.sqrt(10):.3f}" for k, v in sd.items()))
    assert max(se.values()) < h26.PICK_BAR / 2
