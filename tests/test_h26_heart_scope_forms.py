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


def test_a_clean_run_is_read_decides_nothing_and_a_tie_leaves_the_file_in_force(h26):
    res = h26.reading(_parts(h26), _study(h26), _h25(h26), commit="abc")
    assert not res["void"], res["why"]
    assert "pick" not in res and not hasattr(h26, "pick") and not hasattr(h26, "PICK_BAR")
    assert res["decides"].startswith("nothing")
    assert res["conditional_outcome"]["value"] == "A"
    rep = res["reported_decides_nothing"]
    assert rep["label"] == "reported, decides nothing"
    assert rep["V"]["D3"] == 0.0 and rep["plain_mean"]["DJ4"] == 0.0
    assert res["R1_carry_over"]["agreement"]["half_iso_half"] == {"agree": 15, "of": 15}
    t = h26.text(res)
    assert t.splitlines()[0].startswith("This read decides nothing")
    assert "oracle size" in t and "pick" not in t.lower()
    assert ("conditional_outcome = A: used only under the organiser statement the row names; otherwise the file in "
            "force stays") in t


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
    assert res["conditional_outcome"]["value"] == "A" and "void" in res["follows"]
    assert h26.text(res).startswith("This read decides nothing")


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


# ---------------------------------------------------------------------------------------------- the conditional outcome
def _clean(h26, ba: dict[str, float], cb: dict[str, float] | None = None, s1: float = 0.0, own: float = 0.0) -> dict:
    """Parts where, on the CLEAN truths, B = A + ba[key] and C = B + cb[key] (A at 50), and on the truth as released B
    = A + own; s1form_half at 50 + s1 everywhere."""
    cb = cb or {}
    out = {}
    for k in h26.TASK_KEYS:
        out[k] = {}
        for q in h26.TRUTHS:
            g = ba.get(k, 0.0) if q in h26.CLEAN else (own if q == "own" else 0.0)
            h = cb.get(k, 0.0) if q in h26.CLEAN else 0.0
            out[k][q] = {"half_iso_half": 50.0 + g, "core_iso_half": 50.0 + g + h, "s1form_half": 50.0 + s1}
    return out


def _outcome(h26, **kw) -> str:
    res = h26.reading(_parts(h26, score=_clean(h26, **kw)), _study(h26), _h25(h26))
    return res["conditional_outcome"]["value"]


def test_clean_is_the_three_most_trimmed_truths(h26):
    assert h26.CLEAN == ("radial_0", "thin_0.5", "thin_0") and set(h26.CLEAN) < set(h26.TRUTHS)
    assert h26.WORLDS == {"CLEAN": h26.CLEAN, "own": ("own",)}
    st = h26.load_study()
    for key in h26.FOLDS[0]:  # E10.5 core shares 0.40, 0.85 and 1.00, against 0.23 as released
        cs = {q: st[key]["truth_info"][q]["core_share"] for q in h26.TRUTHS}
        assert [round(cs[q], 2) for q in h26.CLEAN] == [0.40, 0.85, 1.00] and round(cs["own"], 2) == 0.23


def test_the_conditional_outcome_rule(h26):
    k0, k1, k2 = h26.TASK_KEYS  # E9.5->E10.5, E8.25_late->E8.5, E8.75->E10.5
    every = lambda g: {k0: g, k1: g, k2: g}  # noqa: E731
    assert _outcome(h26, ba=every(1.0)) == "B"
    assert _outcome(h26, ba=every(0.99)) == "A"
    # the two E10.5 transitions are one fold: 2.0 and 0.0 there make 1.0, with 1.0 on E8.5 the mean is 1.0
    assert _outcome(h26, ba={k0: 2.0, k1: 1.0, k2: 0.0}) == "B"
    # a gain on the E10.5 fold alone: the mean of the folds is 1.0 but E8.5 is not above zero
    assert _outcome(h26, ba={k0: 2.0, k1: 0.0, k2: 2.0}) == "A"
    # a large gain on E8.5 with a loss on the E10.5 fold
    assert _outcome(h26, ba={k0: -0.5, k1: 4.0, k2: -0.5}) == "A"
    # only CLEAN counts: a gain on the truth as released alone names nothing
    assert _outcome(h26, ba=every(0.0), own=5.0) == "A"
    # C only on top of B's condition, and only with C-B at the bar on both folds
    assert _outcome(h26, ba=every(1.5), cb=every(1.0)) == "C (needs a size ruling)"
    assert _outcome(h26, ba=every(1.5), cb=every(0.99)) == "B"
    assert _outcome(h26, ba=every(1.5), cb={k0: 2.0, k1: 0.0, k2: 2.0}) == "B"
    assert _outcome(h26, ba=every(0.5), cb=every(3.0)) == "A"
    # a void run gives A
    res = h26.reading(_parts(h26, score=_clean(h26, ba=every(5.0)), clean=False), _study(h26), _h25(h26))
    assert res["void"] and res["conditional_outcome"]["value"] == "A"
    # the condition travels with the value
    res = h26.reading(_parts(h26, score=_clean(h26, ba=every(1.0))), _study(h26), _h25(h26))
    o = res["conditional_outcome"]
    assert o["condition"] == ("used only under the organiser statement the row names; otherwise the file in force "
                              "stays")
    assert o["inputs"]["B-A"] == pytest.approx(1.0) and o["inputs"]["B-A per fold"] == {"E10.5": pytest.approx(1.0),
                                                                                         "E8.5": pytest.approx(1.0)}


def test_the_stand_in_never_enters_the_primary_table_or_the_outcome(h26):
    """Plain means carry no draw premium: s1form_half's level moves V (secondary) and nothing else."""
    k0, k1, k2 = h26.TASK_KEYS
    ba = {k0: 1.2, k1: 1.1, k2: 1.0}
    lo = h26.reading(_parts(h26, score=_clean(h26, ba=ba, s1=-20.0)), _study(h26), _h25(h26))
    hi = h26.reading(_parts(h26, score=_clean(h26, ba=ba, s1=+0.5)), _study(h26), _h25(h26))
    assert lo["conditional_outcome"] == hi["conditional_outcome"] and lo["conditional_outcome"]["value"] == "B"
    for u in h26.UNITS:
        for q in h26.TRUTHS:
            assert lo["per_world"][u][q]["mean"] == hi["per_world"][u][q]["mean"]
            assert lo["per_world"][u][q]["diff"] == hi["per_world"][u][q]["diff"]
    assert lo["stated_worlds"]["CLEAN"]["fold-merged"]["diff"] == hi["stated_worlds"]["CLEAN"]["fold-merged"]["diff"]
    # V does move: with the stand-in at 50.5 it truncates A at 50.5 and B (51.0 to 51.2) is barely above it
    assert hi["per_world"][k0]["radial_0"]["V_secondary"]["diff"]["B-A"] == pytest.approx(0.7)
    assert lo["per_world"][k0]["radial_0"]["V_secondary"]["diff"]["B-A"] == pytest.approx(1.2)


def test_the_primary_table_per_world_and_fold(h26):
    """Per-seed arithmetic by hand: the E10.5 fold averages its two transitions seed by seed, the mean of the folds
    averages the two folds seed by seed, and min, max and seeds ahead are of those per-seed differences."""
    k0, k1, k2 = h26.TASK_KEYS
    s = np.arange(10, dtype=float)
    b = {k0: list(50 + s - 4.5), k1: list(50 + 2.0 + 0 * s), k2: list(50 - (s - 4.5))}  # k0 + k2 cancel by seed
    c = {k0: list(50 + 1.0 + 0 * s), k1: list(50 - 1.0 + 0 * s), k2: list(50 + 3.0 + 0 * s)}
    score = {k: {"thin_0": {"half_iso_half": b[k], "core_iso_half": c[k]}} for k in h26.TASK_KEYS}
    res = h26.reading(_parts(h26, score=score), _study(h26), _h25(h26))
    w = res["per_world"]
    assert w[k0]["thin_0"]["diff"]["B-A"] == {"mean": pytest.approx(0.0), "lo": -4.5, "hi": 4.5, "ahead": 5,
                                               "behind": 5}
    f = w["fold E10.5"]["thin_0"]
    assert f["diff"]["B-A"]["lo"] == f["diff"]["B-A"]["hi"] == pytest.approx(0.0)
    assert f["diff"]["C-A"]["mean"] == pytest.approx(2.0) and f["diff"]["C-B"]["mean"] == pytest.approx(2.0)
    m = w["mean of folds"]["thin_0"]
    assert m["diff"]["B-A"]["mean"] == pytest.approx(1.0) and m["diff"]["B-A"]["ahead"] == 10
    assert m["diff"]["C-B"]["mean"] == pytest.approx(0.5 * 2.0 + 0.5 * -3.0)
    assert m["mean"] == {"A": pytest.approx(50.0), "B": pytest.approx(51.0), "C": pytest.approx(50.5)}
    assert w["mean of folds"]["own"]["diff"]["C-A"]["mean"] == 0.0
    # CLEAN averages its three truths equally: one of three moved
    cl = res["stated_worlds"]["CLEAN"]
    assert cl["fold-merged"]["diff"]["B-A"]["mean"] == pytest.approx(1.0 / 3)
    assert cl["folds"]["E8.5"]["diff"]["B-A"]["mean"] == pytest.approx(2.0 / 3)
    assert res["stated_worlds"]["own"]["fold-merged"]["diff"]["B-A"]["mean"] == 0.0
    assert set(h26.UNITS) == {*h26.TASK_KEYS, "fold E10.5", "fold E8.5", "mean of folds"}


def test_the_reported_statistics_decide_nothing(h26):
    """D3 and DJ4 (on V and on the plain mean, per prior and at the minimum) are reported and never move the
    conditional outcome."""
    # B ahead by 3 on every truth but CLEAN, where it is level; the stand-in far below
    score = {k: {q: {"half_iso_half": 50.0 + (0.0 if q in h26.CLEAN else 3.0), "s1form_half": 30.0}
                 for q in h26.TRUTHS} for k in h26.TASK_KEYS}
    res = h26.reading(_parts(h26, score=score), _study(h26), _h25(h26))
    rep = res["reported_decides_nothing"]
    assert rep["V"]["per_prior"]["released"] == {"D3": pytest.approx(3.0), "DJ4": pytest.approx(3.0)}
    assert rep["V"]["per_prior"]["thin_family"]["D3"] == pytest.approx(0.0)
    assert rep["plain_mean"]["min_over_priors"]["D3"] == pytest.approx(0.0)
    assert rep["V"]["D3"] >= 1.0 and rep["plain_mean"]["DJ4"] >= 1.0
    assert res["conditional_outcome"]["value"] == "A"
    t = h26.text(res)
    assert "REPORTED, DECIDES NOTHING" in t and "released D3 +3.000" in t and "different draw" in t


def test_the_carry_over_reference_figures_quoted_in_row_h26(h26):
    """Row H26 quotes, as committed reference arithmetic (not a forecast of this run), what the statistics would be if
    the study's own half_radial - hold and core_radial - hold carried over exactly (committed per-seed files): the
    reported D3 and DJ4, and the stated worlds' fold-merged and per-fold plain-mean differences."""
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
    # the stated worlds under exact carry-over: A as hold, B as half_radial, C as core_radial
    twin = {"A": "hold", "B": "half_radial", "C": "core_radial"}

    def d(x, y, keys, truths):
        return float(np.mean([np.mean([r["scores"][q][twin[x]] - r["scores"][q][twin[y]] for r in st[k]["rows"]])
                              for k in keys for q in truths]))
    ref = {}
    for wn, ts in h26.WORLDS.items():
        for x, y in (("B", "A"), ("C", "A"), ("C", "B")):
            folds = [d(x, y, f, ts) for f in h26.FOLDS]
            ref[(wn, f"{x}-{y}")] = (float(np.mean(folds)), *folds)
    print("H26 stated worlds under exact carry-over (fold-merged; E10.5 fold, E8.5 fold): " + "; ".join(
        f"{w} {p} {m:+.2f} ({a:+.2f}, {b:+.2f})" for (w, p), (m, a, b) in ref.items()))
    ba, cb = ref[("CLEAN", "B-A")], ref[("CLEAN", "C-B")]
    print("H26 the conditional outcome's arithmetic on the carry-over reference: "
          + h26.outcome_from(ba[0], ba[1:], cb[0], cb[1:], True))


# ------------------------------------------------------------------------------------------------ noise
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


def _clean_null(h26, sd: dict[str, float], draws: int, seed: int, share_truths: bool, share_e105: bool) -> dict:
    """Pure noise for the conditional outcome: A, B and C equal in distribution. Per transition, CLEAN truth and seed,
    each arm is an independent draw of SD sd / sqrt(2), so every paired difference has the committed spread sd. A
    component common to the arms cancels and is left out. share_truths: one draw for the three CLEAN truths of a
    transition; share_e105: one draw for the two E10.5 transitions. The stand-in does not enter. Returns per draw the
    per-seed scores averaged over CLEAN, per transition (draws x 3 transitions x 3 arms x seeds)."""
    rng = np.random.default_rng(seed)
    nt = 1 if share_truths else len(h26.CLEAN)
    cache, out = {}, []
    for key in h26.TASK_KEYS:
        src = h26.TASK_KEYS[0] if (share_e105 and key == h26.TASK_KEYS[2]) else key
        if src not in cache:
            cache[src] = rng.normal(size=(draws, 3, nt, len(SEEDS)))
        out.append(sd[key] / np.sqrt(2) * cache[src].mean(2))  # CLEAN: the truths weighted equally
    return np.stack(out, 1)


def _clean_outcomes(h26, x: np.ndarray) -> dict:
    """The conditional outcome's arithmetic on draws x transitions x arms (A, B, C) x seeds, vectorised."""
    f = np.stack([0.5 * (x[:, 0] + x[:, 2]), x[:, 1]], 1)  # folds E10.5, E8.5: draws x 2 x arms x seeds
    ba_f = (f[:, :, 1] - f[:, :, 0]).mean(-1)
    cb_f = (f[:, :, 2] - f[:, :, 1]).mean(-1)
    ba, cb = ba_f.mean(1), cb_f.mean(1)
    b = (ba >= h26.OUTCOME_BAR) & (ba_f > 0).all(1)
    c = b & (cb >= h26.OUTCOME_BAR) & (cb_f > 0).all(1)
    return {"ba": ba, "ba_f": ba_f, "cb": cb, "cb_f": cb_f, "B": b & ~c, "C": c}


def test_the_vectorised_null_is_the_scripts_outcome(h26):
    """The noise simulation's arithmetic equals the script's reading on hand-made parts built from the same draws
    (each CLEAN truth given the draw, the other truths level), whatever the stand-in."""
    sd = {k: 6.0 for k in h26.TASK_KEYS}
    x = _clean_null(h26, sd, draws=12, seed=3, share_truths=True, share_e105=False)
    o = _clean_outcomes(h26, x)
    named = 0
    for i in range(len(x)):
        score = {k: {q: {"full_iso_half": list(50 + x[i, j, 0]), "half_iso_half": list(50 + x[i, j, 1]),
                         "core_iso_half": list(50 + x[i, j, 2]), "s1form_half": list(50 + x[i, j, 1] + 0.3)}
                     for q in h26.CLEAN} for j, k in enumerate(h26.TASK_KEYS)}
        res = h26.reading(_parts(h26, score=score), _study(h26), _h25(h26))
        inp = res["conditional_outcome"]["inputs"]
        assert inp["B-A"] == pytest.approx(o["ba"][i], abs=1e-9) and inp["C-B"] == pytest.approx(o["cb"][i], abs=1e-9)
        want = "C (needs a size ruling)" if o["C"][i] else ("B" if o["B"][i] else "A")
        assert res["conditional_outcome"]["value"] == want
        named += want != "A"
    assert 0 < named < len(x)  # both branches exercised at this spread


def test_the_conditional_outcome_on_pure_noise(h26):
    """A, B and C equal in distribution over CLEAN; the study's largest and median paired spreads; truths and the two
    E10.5 transitions independent or sharing draws. How often B or C is named. Printed for row H26."""
    spreads = {"largest": _study_spreads(h26, max), "median": _study_spreads(h26, np.median)}
    print("\nH26 null spreads (sd of paired differences, board points): " + "; ".join(
        f"{sn} " + ", ".join(f"{k} {v:.3f}" for k, v in sd.items()) for sn, sd in spreads.items()))
    rate = {}
    for sn, sd in spreads.items():
        for share, kw in (("independent truths", dict(share_truths=False, share_e105=False)),
                          ("shared truths", dict(share_truths=True, share_e105=False)),
                          ("shared truths and E10.5", dict(share_truths=True, share_e105=True))):
            o = _clean_outcomes(h26, _clean_null(h26, sd, draws=20000, seed=20261119, **kw))
            rate[(sn, share)] = r = {"mean B-A": float(o["ba"].mean()), "B named": float(o["B"].mean()),
                                     "C named": float(o["C"].mean()), "B or C": float((o["B"] | o["C"]).mean())}
            print(f"H26 conditional outcome on noise | spreads {sn} | {share}: mean B-A {r['mean B-A']:+.3f}, "
                  f"B named {r['B named']:.2%}, C named {r['C named']:.2%}, B or C {r['B or C']:.2%}")
    worst = max(r["B or C"] for r in rate.values())
    print(f"H26 the conditional outcome's worst false naming of B or C on pure noise: {worst:.2%}")
    # plain means carry no premium: the null difference centres on zero
    assert all(abs(r["mean B-A"]) < 0.03 for r in rate.values())
    assert all(r["B or C"] < 0.002 for (sn, _), r in rate.items() if sn == "median")
    assert worst < 0.10


def _v_null_draws(h26, sd: dict[str, float], s1_sd: dict[str, float], s1_level: dict[str, float], draws: int,
                  seed: int, share_truths: bool, share_e105: bool) -> dict[str, np.ndarray]:
    """The first fixer's null for V: B and A the same in distribution, s1form_half A's draw plus s1_level plus a small
    term of SD s1_sd (it shares A's pool, build seed and half rows). Per draw, the per-transition, per-prior B - A of
    E V and of the expected plain mean (draws x 3 x 5)."""
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


def _h25_stand_in(h26) -> tuple[dict[str, float], dict[str, float]]:
    """From row H25's committed per-seed scores: s1form_half minus iso_half, its mean and its seed SD, per E10.5
    backtest; E8.25_late -> E8.5 has no H25 run and takes level 0 and the smaller of the two SDs."""
    h25 = h26.load_h25()
    level, sd = {}, {}
    for key, r in h25.items():
        d = np.array([x["scores"]["s1form_half"] - x["scores"]["iso_half"] for x in r["rows"]])
        level[key], sd[key] = float(d.mean()), float(d.std(ddof=1))
    k1 = h26.TASK_KEYS[1]
    level[k1], sd[k1] = 0.0, min(sd.values())
    return level, sd


def _d3(per: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """D3 and DJ4 of per-draw arrays draws x transitions (TASK_KEYS order) x priors."""
    return per.mean(axis=(1, 2)), 0.5 * (0.5 * (per[:, 0] + per[:, 2]).mean(1) + per[:, 1].mean(1))


def test_the_v_null_arithmetic_is_the_scripts_reported_d3(h26):
    level, s1sd = _h25_stand_in(h26)
    sim = _v_null_draws(h26, _study_spreads(h26), s1sd, level, draws=5, seed=7, share_truths=False, share_e105=False)
    d3, dj4 = _d3(sim["V"])
    m3, mj4 = _d3(sim["mean"])
    for i in range(5):
        ev = {k: {n: {"full_iso_half": 0.0, "half_iso_half": float(sim["V"][i, j, p])}
                  for p, n in enumerate(h26.PRIORS)} for j, k in enumerate(h26.TASK_KEYS)}
        em = {k: {n: {"full_iso_half": 0.0, "half_iso_half": float(sim["mean"][i, j, p])}
                  for p, n in enumerate(h26.PRIORS)} for j, k in enumerate(h26.TASK_KEYS)}
        rep = h26.reported_d3_dj4(ev, em)
        assert rep["V"]["D3"] == pytest.approx(d3[i], abs=1e-12) and rep["V"]["DJ4"] == pytest.approx(dj4[i], abs=1e-12)
        assert rep["plain_mean"]["D3"] == pytest.approx(m3[i], abs=1e-12)
        assert rep["plain_mean"]["DJ4"] == pytest.approx(mj4[i], abs=1e-12)


def test_v_carries_a_draw_premium_on_pure_noise(h26):
    """Why V decides nothing: on pure noise (B and A the same in distribution), with s1form_half sharing A's draw as it
    does in the run, E V favours B with no scope effect, while the plain mean centres on zero. Printed for row H26."""
    level_h25, s1sd = _h25_stand_in(h26)
    levels = {"H25's levels": level_h25, "level": {k: 0.0 for k in h26.TASK_KEYS},
              "far below": {k: -100.0 for k in h26.TASK_KEYS}}
    spreads = {"largest": _study_spreads(h26, max), "median": _study_spreads(h26, np.median)}
    print("\nH26 stand-in s1form_half - iso_half (row H25): " + ", ".join(
        f"{k} level {level_h25[k]:+.3f} sd {s1sd[k]:.3f}" for k in h26.TASK_KEYS))
    null = {}
    for sn, sd in spreads.items():
        for ln, lv in levels.items():
            for share, kw in (("independent truths", dict(share_truths=False, share_e105=False)),
                              ("shared truths", dict(share_truths=True, share_e105=False)),
                              ("shared truths and E10.5", dict(share_truths=True, share_e105=True))):
                sim = _v_null_draws(h26, sd, s1sd, lv, draws=4000, seed=20261119, **kw)
                null[(sn, ln, share)] = (float(_d3(sim["V"])[0].mean()), float(_d3(sim["mean"])[0].mean()))
                print(f"H26 V premium | spreads {sn} | s1 {ln} | {share}: null D3 on V "
                      f"{null[(sn, ln, share)][0]:+.3f}, on the plain mean {null[(sn, ln, share)][1]:+.3f}")
    within = [v for (sn, ln, _), (v, _m) in null.items() if ln != "far below"]
    print(f"H26 null D3 on V with the stand-in within reach: median spreads "
          f"{min(v for (sn, ln, _), (v, _m) in null.items() if ln != 'far below' and sn == 'median'):+.2f} to "
          f"{max(v for (sn, ln, _), (v, _m) in null.items() if ln != 'far below' and sn == 'median'):+.2f}, largest "
          f"{min(v for (sn, ln, _), (v, _m) in null.items() if ln != 'far below' and sn == 'largest'):+.2f} to "
          f"{max(v for (sn, ln, _), (v, _m) in null.items() if ln != 'far below' and sn == 'largest'):+.2f}")
    assert all(v > 0.2 for v in within)
    assert all(abs(m) < 0.05 for (_v, m) in null.values())  # the plain mean carries no premium
    assert all(abs(v) < 0.05 for (sn, ln, _), (v, _m) in null.items() if ln == "far below")


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
    assert max(se.values()) < h26.OUTCOME_BAR / 2
