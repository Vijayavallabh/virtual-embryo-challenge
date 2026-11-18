"""experiments/harness/h25_heart_forms.py (try H25 of experiments/harness/TRIES.md) on synthetic arrays and
hand-made result tables, and the pinned shape metrics on a synthetic cloud. Nothing is read from data/ or runs/; the
committed copy scores of earlier heart runs (experiments/harness/results/*/tasks) are read for the guard's
arithmetic only, and no score of any other candidate is read."""

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
from vec.evaluate import Stage, panel_metrics

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "experiments" / "harness" / "h25_heart_forms.py"
RESULTS = REPO_ROOT / "experiments" / "harness" / "results"


@pytest.fixture(scope="module")
def h25():
    spec = importlib.util.spec_from_file_location("h25_heart_forms_module", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def _cloud(n=2000, seed=0, sd=(300.0, 200.0, 80.0)):
    """A synthetic anisotropic cloud (x widest, z thinnest, as a heart section stack)."""
    return np.random.default_rng(seed).normal(size=(n, 3)) * np.asarray(sd) + np.array([1000.0, -500.0, 40.0])


def _size(h25, target=500.0, sd_z=80.0):
    return h25.oracle_size_result(target, "E9.5", 9.5, 350.0, sd_z, 10.5)


# ------------------------------------------------------------------------------------------------ the size step
def test_the_procedure_has_no_figure_for_e10_5_on_either_backtest(h25):
    """At its pre-registered settings the size procedure applies no rule to either backtest's inputs (at most one
    stage at or after E9.0), so its result would be valid=False: the oracle is the only size."""
    for t in h25.tasks():
        days = [harness.STAGE_DAYS[f] for f in t.inputs.values()]
        assert h25.oracle_note(days) == []
    # the test layout has an interval starting at E9.5 and gets rules: the check is not empty by construction
    assert h25.oracle_note([8.3, 8.5, 8.75, 9.5, 10.5]) == ["latest_rate", "half_rate", "power_t"]


def test_the_raw_mode_is_an_isotropic_rescale_to_the_oracle_size(h25):
    C = _cloud()
    out, rec = h25.size_step(C, _size(h25), "raw")
    mu = C.mean(0)
    f = t2_build.rms_radius(out) / t2_build.rms_radius(C)
    assert t2_build.rms_radius(out) == pytest.approx(500.0, rel=1e-12)
    np.testing.assert_allclose(out, mu + (C - mu) * f, rtol=0, atol=1e-9)
    assert rec["label"] == "oracle size" and rec["z_mode_used"] == "raw" and rec["z_held"] is False


def test_the_recorded_mode_holds_z_and_reaches_the_same_size(h25):
    C = _cloud()
    out, rec = h25.size_step(C, _size(h25), "xy_grow_z_hold")
    assert rec["z_mode_used"] == "xy_grow_z_hold" and rec["z_held"] is True
    np.testing.assert_array_equal(out[:, 2], C[:, 2])
    assert t2_build.rms_radius(out) == pytest.approx(500.0, rel=1e-12)
    iso, _ = h25.size_step(C, _size(h25), "raw")
    # same size, different shape: the arms differ in form only
    assert t2_build.rms_radius(iso) == pytest.approx(t2_build.rms_radius(out), rel=1e-12)
    assert np.abs(iso - out).max() > 1.0


def test_the_size_step_is_the_builders_function(h25):
    """The size goes through t2_build.apply_size_procedure, which size_or_hold calls for a valid heart JSON."""
    src = inspect.getsource(h25.size_step)
    assert "t2_build.apply_size_procedure" in src
    assert "apply_size_procedure(coords, size_json, hypothesis" in inspect.getsource(t2_build.size_or_hold)


def test_the_oracle_refuses_a_size_that_cannot_hold_z(h25):
    with pytest.raises(ValueError):
        h25.oracle_size_result(70.0, "E9.5", 9.5, 350.0, 80.0, 10.5)
    with pytest.raises(ValueError):
        h25.oracle_size_result(float("nan"), "E9.5", 9.5, 350.0, 80.0, 10.5)


# ------------------------------------------------------------------------------------------------ the halving
def test_the_half_is_halve_pys_rows_at_the_documented_seed(h25):
    sys.path.insert(0, str(REPO_ROOT / "experiments" / "builds"))
    import halve

    n = 1237
    for s in (0, 1, 9):
        np.testing.assert_array_equal(h25.kept_rows(n, s, 1), halve.half_rows(n, h25.HALF_SEED_BASE + s))
        np.testing.assert_array_equal(h25.kept_rows(n, s, 0), np.arange(n))
    assert h25.halve.half_rows is halve.half_rows


def test_the_quarter_is_the_half_halved_and_nested_in_it(h25):
    n = 1237
    for s in SEEDS:
        half, quarter = h25.kept_rows(n, s, 1), h25.kept_rows(n, s, 2)
        assert len(half) == n // 2 and len(quarter) == (n // 2) // 2
        assert set(quarter) <= set(half) and np.all(np.diff(quarter) > 0)


def test_the_halving_seeds_avoid_every_build_and_scoring_stream(h25):
    seeds = {h25.HALF_SEED_BASE + s for s in SEEDS} | {h25.QUARTER_SEED_BASE + s for s in SEEDS}
    builds = {harness.BUILD_SEED_OFFSET + s for s in SEEDS}
    assert not seeds & builds and not seeds & set(SEEDS) and len(seeds) == 2 * len(SEEDS)
    assert min(seeds) >= harness.BUILD_SEED_OFFSET  # halve.py's floor
    sys.path.insert(0, str(REPO_ROOT / "experiments" / "builds"))
    import halve

    # the upload's seed is a build seed number; kept apart so that no seed number serves two purposes (for clarity:
    # sample_cells draws its rows from default_rng([seed, 0]), not halve.py's default_rng(seed))
    assert halve.DEFAULT_SEED in builds
    assert "default_rng([seed, 0])" in inspect.getsource(t2_build.sample_cells)


def test_the_halving_draw_is_uniform(h25):
    """Over the ten seeds every row is kept about half the time (a uniform half, not a block)."""
    n = 400
    kept = np.zeros(n)
    for s in range(200):
        kept[h25.kept_rows(n, s % 10 + 10 * (s // 10), 1)] += 1
    share = kept / 200
    assert abs(share.mean() - 0.5) < 1e-12 and share.min() > 0.3 and share.max() < 0.7


# ------------------------------------------------------------------------------------------------ the forms
def test_each_arm_can_differ_from_its_control(h25):
    """iso_half differs from copy_half (positions), copy_half from the copy (rows), copy_quarter from copy_half
    (rows), s1form from b0_trend's positions, s1form_half from s1form (rows), and iso_half from s1form_half's form."""
    rng = np.random.default_rng(1)
    C = _cloud(n=1001)
    X = rng.random((1001, 5))
    size = _size(h25)
    copy = h25.form(X, C, size, None, 0, 3)
    half = h25.form(X, C, size, None, 1, 3)
    quarter = h25.form(X, C, size, None, 2, 3)
    iso = h25.form(X, C, size, "raw", 1, 3)
    s1 = h25.form(X, C, size, "xy_grow_z_hold", 0, 3)
    s1h = h25.form(X, C, size, "xy_grow_z_hold", 1, 3)
    assert len(half[0]) == 500 and len(quarter[0]) == 250 and len(copy[0]) == 1001
    np.testing.assert_array_equal(iso[0], half[0])  # same cells and expression
    assert np.abs(iso[1] - half[1]).max() > 1.0  # other positions
    assert np.abs(s1[1] - C).max() > 1.0 and len(s1[0]) == 1001
    np.testing.assert_array_equal(s1h[1], s1[1][h25.kept_rows(1001, 3, 1)])
    assert np.abs(s1h[1] - iso[1]).max() > 1.0
    # iso_half is copy_half rescaled about copy's centroid
    mu = C.mean(0)
    f = iso[2]["rms_after"] / iso[2]["rms_before"]
    np.testing.assert_allclose(iso[1], mu + (half[1] - mu) * f, atol=1e-9)
    assert s1[2]["label"] == "oracle size"


def test_the_arms_apply_to_both_counted_backtests_and_use_the_builders_expression(h25):
    ts = h25.tasks()
    assert [t.key for t in ts] == list(h25.TASK_KEYS) and all(t.counted for t in ts)
    assert {t.target for t in ts} == {"T2/E10.5.h5ad"}  # one truth: one fold of evidence
    for name in h25.ARMS:
        arm = h25.make_arm(name, _size(h25), {})
        assert all(arm.applies(t) for t in ts)
    for name in h25.HARNESS_REFERENCES:
        assert all(harness.CANDIDATES[name].applies(t) for t in ts)
    assert all(h25.make_arm("copy_iso", _size(h25), {}).applies(t) for t in ts)
    assert h25.REFERENCES == ("copy", "b0_trend", "copy_iso") and "copy_iso" not in h25.ARMS  # not counted
    # the builder's --trend expression is the harness's b0_trend, the function the s1form arms build from
    build_src = REPO_ROOT / "experiments" / "builds" / "build_test.py"
    assert "harness._b0_t2(ctx, trend=True, trend_record=trend_rec)" in build_src.read_text()
    assert "_b0_t2(ctx, trend=True)" in inspect.getsource(harness._b0_trend)
    assert "harness._b0_trend" in inspect.getsource(h25.make_arm)


# ------------------------------------------------------------------------------------------------ the reading
def _skills(**over):
    s = {"de_score": 0.5, "de_direction": 0.5, "mmd_u": 0.5, "variogram": 0.5, "d2_shape": 0.5,
         "occupancy_dice": 0.5, "scale_log_ratio": 0.5, "neighborhood_mmd": 0.5}
    return {**s, **over}


def _groups():
    return {g: 0.125 for g in ("expression_change", "cell_state", "shape_growth", "local_organisation")}


def _result(h25, key, score=None, skills=None, commit="abc", clean=True):
    """A hand-made task result: score[arm] and skills[arm] are per-seed lists or constants."""
    score, skills = score or {}, skills or {}
    rows = []
    for k, s in enumerate(SEEDS):
        sc = {a: (score.get(a, 50.0)[k] if isinstance(score.get(a), list) else score.get(a, 50.0)) for a in h25.ALL}
        sk = {a: _skills(**{m: (v[k] if isinstance(v, list) else v) for m, v in skills.get(a, {}).items()})
              for a in h25.ALL}
        rows.append({"seed": s, "scores": sc, "skills": sk, "groups": {a: _groups() for a in h25.ALL}})
    return {"task": {"key": key}, "rows": rows, "provenance": {"commit": commit, "tree_clean": clean}}


def _both(h25, **kw):
    return {k: _result(h25, k, **kw) for k in h25.TASK_KEYS}


def _committed(h25, copy=50.0):
    return {k: _result(h25, k, score={"copy": copy}) for k in h25.TASK_KEYS}


def test_a_clean_run_is_read_and_a_tie_makes_only_statement_one_true(h25):
    res = h25.reading(_both(h25), _committed(h25))
    assert not res["void"]
    st = res["statements"]
    assert st["1"]["true"] and not st["2"]["true"] and not st["3"]["true"]
    assert st["2"]["follows"] == "nothing"


def test_statement_two_needs_the_bar_and_eight_seeds_on_both_backtests(h25):
    n = len(SEEDS)

    def half(v):  # copy_half and its rescaled twin iso_half move together, as the twin check requires
        return {"copy_half": {"neighborhood_mmd": v}, "iso_half": {"neighborhood_mmd": v}}

    lose = half([0.47] * n)  # -0.03 in every seed
    assert h25.reading(_both(h25, skills=lose), _committed(h25))["statements"]["2"]["true"]
    assert "priced question" in h25.reading(_both(h25, skills=lose), _committed(h25))["statements"]["2"]["follows"]
    small = half([0.485] * n)  # -0.015: under the bar
    assert not h25.reading(_both(h25, skills=small), _committed(h25))["statements"]["2"]["true"]
    seven = half([0.40] * 7 + [0.51] * 3)  # mean -0.067, 7 of 10 behind
    assert not h25.reading(_both(h25, skills=seven), _committed(h25))["statements"]["2"]["true"]
    one = _both(h25)
    one[h25.TASK_KEYS[0]] = _result(h25, h25.TASK_KEYS[0], skills=lose)  # one backtest only
    assert not h25.reading(one, _committed(h25))["statements"]["2"]["true"]


def test_statement_three_needs_the_bar_and_eight_seeds_on_both_backtests(h25):
    n = len(SEEDS)
    assert h25.reading(_both(h25, score={"s1form_half": [49.5] * n}), _committed(h25))["statements"]["3"]["true"]
    assert not h25.reading(_both(h25, score={"s1form_half": [49.8] * n}), _committed(h25))["statements"]["3"]["true"]
    st = h25.reading(_both(h25, score={"s1form_half": [49.5] * n}), _committed(h25))["statements"]["3"]
    assert "decides nothing" in st["follows"]


def test_statement_three_reads_against_the_sized_copy_not_the_copy(h25):
    """The oracle lifts every sized arm's size term (about +4.2 points over the unsized copy), so (3)'s control is
    copy_iso: an s1form_half 4 points above the copy but 0.5 below copy_iso makes (3) true, and one 0.5 below the copy
    but level with copy_iso makes it false."""
    n = len(SEEDS)
    above_copy = {"copy_iso": [54.2] * n, "s1form_half": [53.7] * n}
    res = h25.reading(_both(h25, score=above_copy), _committed(h25))
    assert res["statements"]["3"]["true"]
    assert res["minus_copy"][h25.TASK_KEYS[0]]["s1form_half"]["score"]["mean"] == pytest.approx(3.7)
    assert res["minus_copy_iso"][h25.TASK_KEYS[0]]["s1form_half"]["score"]["mean"] == pytest.approx(-0.5)
    level = {"copy_iso": [49.5] * n, "s1form_half": [49.5] * n}
    assert not h25.reading(_both(h25, score=level), _committed(h25))["statements"]["3"]["true"]
    assert set(res["minus_copy_iso"][h25.TASK_KEYS[0]]) == {"iso_half", "s1form", "s1form_half"}
    assert "not comparable" in h25.text(res)


def test_statement_one_fails_when_the_ruled_form_moves_a_shape_skill(h25):
    n = len(SEEDS)
    moved = {"iso_half": {"occupancy_dice": [0.45] * n}}  # -0.05 against copy_half
    assert not h25.reading(_both(h25, skills=moved), _committed(h25))["statements"]["1"]["true"]
    inside = {"iso_half": {"d2_shape": [0.515] * n}}  # +0.015
    assert h25.reading(_both(h25, skills=inside), _committed(h25))["statements"]["1"]["true"]


def _noise_rates(h25, sd_score, sd_skill, shared, draws, seed):
    """Share of pure-noise draws on which each statement comes out true. Every control (copy, copy_iso) is held level
    and every other arm gets per-seed noise, so a paired difference with a control has sd `sd_score` (board points) or
    `sd_skill`; shared=True draws one noise for both backtests (they share one truth, E10.5)."""
    rng = np.random.default_rng(seed)
    hits = {"1": 0, "2": 0, "3": 0}
    controls = ("copy", h25.SIZED_CONTROL)
    for _ in range(draws):
        results, noise = {}, None
        for key in h25.TASK_KEYS:
            if noise is None or not shared:
                noise = ({a: list(50.0 + rng.normal(0, sd_score, len(SEEDS))) for a in h25.ALL[1:]
                          if a not in controls},
                         {a: {m: list(0.5 + rng.normal(0, sd_skill, len(SEEDS)))
                              for m in ("d2_shape", "occupancy_dice", "neighborhood_mmd")}
                          for a in h25.ALL[1:] if a not in controls})
            results[key] = _result(h25, key, score=noise[0], skills=noise[1])
        st = h25.statements(results)
        for k in hits:
            hits[k] += st[k]["true"]
    return {k: v / draws for k, v in hits.items()}


def test_the_gates_on_pure_noise(h25):
    """Pure noise at the spread of an arm minus the copy on these backtests: board points sd 0.45 (row H23's committed
    arm-minus-copy spreads are 0.21 to 0.26 on E8.75->E10.5 and 0.39 to 0.46 on E9.5->E10.5); skills sd 0.01 (the
    same arms' neighborhood_mmd skill spreads are 0.001 to 0.005). Over 2000 draws, with the backtests independent and
    with one noise draw shared by both, statements (2) and (3) come out true in under 2.5%, and statement (1) true in
    most (it asks for agreement, so noise below its tolerance passes it)."""
    ind = _noise_rates(h25, 0.45, 0.01, shared=False, draws=2000, seed=20261118)
    cor = _noise_rates(h25, 0.45, 0.01, shared=True, draws=2000, seed=20261119)
    assert ind["2"] < 0.005 and ind["3"] < 0.005
    assert cor["2"] < 0.005 and cor["3"] < 0.025
    assert ind["1"] > 0.5 and cor["1"] > 0.5


def test_the_reading_is_void_when_the_copy_drifts_from_its_committed_scores(h25):
    n = len(SEEDS)
    assert h25.reading(_both(h25, score={"copy": [50.09] * n}), _committed(h25))["void"] is False
    assert h25.reading(_both(h25, score={"copy": [50.11] * n}), _committed(h25))["void"]  # mean 0.11
    spike = [50.0] * (n - 1) + [50.6]  # mean 0.06, one seed 0.6
    res = h25.reading(_both(h25, score={"copy": spike}), _committed(h25))
    assert res["void"] and any("a seed's copy score" in w for w in res["why"])
    assert res["statements"] is None
    inside = [50.0] * (n - 1) + [50.45]
    assert not h25.reading(_both(h25, score={"copy": inside}), _committed(h25))["void"]


def test_the_reading_is_void_on_missing_or_mixed_results(h25):
    r = _both(h25)
    del r[h25.TASK_KEYS[1]]
    assert h25.reading(r, _committed(h25))["void"]
    r = _both(h25)
    r[h25.TASK_KEYS[0]]["provenance"]["commit"] = "def"
    assert h25.reading(r, _committed(h25))["void"]
    r = _both(h25, clean=False)
    assert h25.reading(r, _committed(h25))["void"]
    r = _both(h25)
    r[h25.TASK_KEYS[0]]["rows"] = r[h25.TASK_KEYS[0]]["rows"][:9]
    assert h25.reading(r, _committed(h25))["void"]
    r = _both(h25)
    del r[h25.TASK_KEYS[0]]["rows"][2]["scores"]["iso_half"]
    assert h25.reading(r, _committed(h25))["void"]


def test_the_reading_is_void_when_twins_differ_on_expression(h25):
    n = len(SEEDS)
    bad = {"iso_half": {"mmd_u": [0.5] * (n - 1) + [0.51]}}
    res = h25.reading(_both(h25, skills=bad), _committed(h25))
    assert res["void"] and any("iso_half mmd_u differs from copy_half's" in w for w in res["why"])
    bad = {"s1form": {"de_score": [0.6] * n}}
    assert h25.reading(_both(h25, skills=bad), _committed(h25))["void"]
    bad = {"copy_iso": {"de_direction": [0.6] * n}}
    assert h25.reading(_both(h25, skills=bad), _committed(h25))["void"]


def test_the_reading_is_void_when_a_rescaled_twin_moves_its_neighbourhoods(h25):
    n = len(SEEDS)
    bad = {"copy_iso": {"neighborhood_mmd": [0.5] * (n - 1) + [0.49]}}
    res = h25.reading(_both(h25, skills=bad), _committed(h25))
    assert res["void"] and any("copy_iso neighborhood_mmd differs from copy's" in w for w in res["why"])
    bad = {"iso_half": {"neighborhood_mmd": [0.51] * n}}
    assert h25.reading(_both(h25, skills=bad), _committed(h25))["void"]
    rounding = {"iso_half": {"neighborhood_mmd": [0.5001] * n}}  # under NBHD_TWIN_TOL
    assert not h25.reading(_both(h25, skills=rounding), _committed(h25))["void"]
    moved = {"s1form": {"neighborhood_mmd": [0.45] * n}}  # an in-plane stretch may move it: not a twin check
    assert not h25.reading(_both(h25, skills=moved), _committed(h25))["void"]


def test_the_reading_is_void_unless_the_results_carry_the_rows_commit(h25):
    assert not h25.reading(_both(h25, commit="abc123def"), _committed(h25), commit="abc123")["void"]
    res = h25.reading(_both(h25, commit="abc123def"), _committed(h25), commit="fff000")
    assert res["void"] and any("row's commit" in w for w in res["why"])
    assert h25.reading(_both(h25), _committed(h25), commit="")["void"]


def test_the_guard_on_the_committed_heart_runs(h25):
    """The arithmetic quoted in row H25: the copy's scores in every committed run of the two heart backtests against
    row H11's run (the guard's reference). Earlier reruns reproduce it exactly, so the guard's limits are loose."""
    ref = h25.load_dir(REPO_ROOT / h25.GUARD_RUN)
    assert set(ref) == set(h25.TASK_KEYS)
    seen = 0
    for p in sorted(RESULTS.glob("*/tasks/T2h_*.json")):
        r = json.loads(p.read_text())
        key = r.get("task", {}).get("key")
        if key not in h25.TASK_KEYS or not all("copy" in row["scores"] for row in r["rows"]):
            continue
        if sorted(row["seed"] for row in r["rows"]) != list(SEEDS):
            continue
        assert h25.guard(r, ref[key]) == []
        d = h25.guard_distances(r, ref[key])
        assert d["mean"] == 0.0 and d["largest_seed"] == 0.0, p
        seen += 1
    assert seen >= 30  # 17 committed runs on each backtest at the time of writing, H11's own included


# ------------------------------------------------------------------------------------------------ the shape metrics
def _stage(C, rng):
    X = np.log1p(rng.poisson(2.0, size=(len(C), 12))).astype(np.float32)
    return Stage(X=X, coords=np.asarray(C, dtype=np.float32), celltype=np.array(["a"] * len(C)))


@pytest.fixture(scope="module")
def shape_values():
    """The pinned T2 panel (build main) as vec.evaluate.panel_metrics calls it, for one synthetic prediction under a
    uniform rescale and under an in-plane stretch with z held, against one synthetic truth."""
    rng = np.random.default_rng(5)
    truth_C = _cloud(n=600, seed=11)
    pred_C = _cloud(n=600, seed=12)
    mu = pred_C.mean(0)
    f = 1.7
    uniform = mu + (pred_C - mu) * f
    stretched = pred_C.copy()
    stretched[:, :2] = mu[:2] + (pred_C[:, :2] - mu[:2]) * f
    truth, ref = _stage(truth_C, rng), _stage(_cloud(n=600, seed=13), rng)
    X = _stage(pred_C, rng).X
    out = {}
    for name, C in (("base", pred_C), ("uniform", uniform), ("stretched", stretched)):
        pred = Stage(X=X, coords=np.asarray(C, dtype=np.float32), celltype=np.array(["a"] * len(C)))
        out[name] = panel_metrics("T2", pred, truth, ref, seed=0, build="main")
    return out, f


def test_d2_shape_and_occupancy_dice_ignore_a_uniform_rescale(shape_values):
    v, _ = shape_values
    for m in ("d2_shape", "occupancy_dice"):
        assert v["uniform"][m] == pytest.approx(v["base"][m], abs=1e-4), m


def test_neighborhood_mmd_ignores_a_uniform_rescale_within_the_twin_tolerance(h25, shape_values):
    """Why the reading may hold copy_iso to the copy, and iso_half to copy_half, on neighborhood_mmd."""
    v, _ = shape_values
    assert abs(v["uniform"]["neighborhood_mmd"] - v["base"]["neighborhood_mmd"]) <= h25.NBHD_TWIN_TOL / 10


def test_d2_shape_and_occupancy_dice_see_an_in_plane_stretch_with_z_held(shape_values):
    v, _ = shape_values
    for m in ("d2_shape", "occupancy_dice"):
        assert abs(v["stretched"][m] - v["base"][m]) > 1e-3, m


def test_scale_log_ratio_moves_by_the_log_of_the_factor(shape_values):
    v, f = shape_values
    assert v["uniform"]["scale_log_ratio"] - v["base"]["scale_log_ratio"] == pytest.approx(np.log(f), abs=1e-4)


# ------------------------------------------------------------------------------------------------ end to end, synthetic
@pytest.fixture(scope="module")
def synth():
    """The synthetic heart stages, board and context of tests/test_b0_trend.py (reused, not copied)."""
    spec = importlib.util.spec_from_file_location("b0_trend_tests_for_h25", REPO_ROOT / "tests" / "test_b0_trend.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    inputs = {"ref": mod._stage(1600, seed=2, shift=0.0, tiny=60),
              "start": mod._stage(1100, seed=3, shift=0.8, tiny=10),
              "earlier": mod._stage(900, seed=5, shift=1.6, tiny=10)}
    return mod, inputs


def test_every_arm_builds_through_the_harness_and_the_checked_builder(h25, synth):
    mod, inputs = synth
    ref = np.asarray(inputs["ref"].obsm["spatial_3D"], dtype=float)
    target = 1.6 * t2_build.rms_radius(ref)
    size = h25.oracle_size_result(target, "E9.5", 9.5, t2_build.rms_radius(ref), float(ref[:, 2].std()), 10.5)
    records: dict = {}
    n = mod.BOARD.max_cells
    for seed in (harness.BUILD_SEED_OFFSET + 1, harness.BUILD_SEED_OFFSET + 4):
        ctx = mod._ctx(inputs, seed=seed)
        st = {name: harness.CANDIDATES[name].build(ctx) for name in h25.HARNESS_REFERENCES}
        st |= {name: h25.make_arm(name, size, records).build(ctx) for name in ("copy_iso", *h25.ARMS)}
        assert {k: v.n for k, v in st.items()} == {"copy": n, "b0_trend": n, "copy_iso": n, "copy_half": n // 2,
                                                   "copy_quarter": (n // 2) // 2, "iso_half": n // 2, "s1form": n,
                                                   "s1form_half": n // 2}
        s = seed - harness.BUILD_SEED_OFFSET
        half, quarter = h25.kept_rows(n, s, 1), h25.kept_rows(n, s, 2)
        np.testing.assert_array_equal(st["copy_half"].X, st["copy"].X[half])
        np.testing.assert_array_equal(st["copy_half"].coords, st["copy"].coords[half])
        np.testing.assert_array_equal(st["copy_quarter"].X, st["copy"].X[quarter])
        np.testing.assert_array_equal(st["iso_half"].X, st["copy_half"].X)
        np.testing.assert_array_equal(st["copy_iso"].X, st["copy"].X)
        assert t2_build.rms_radius(st["copy_iso"].coords) == pytest.approx(target, rel=1e-5)
        np.testing.assert_allclose(st["iso_half"].coords, st["copy_iso"].coords[half], atol=1e-3)
        np.testing.assert_array_equal(st["s1form"].X, st["b0_trend"].X)
        np.testing.assert_array_equal(st["s1form"].coords[:, 2], st["b0_trend"].coords[:, 2])  # z held
        np.testing.assert_array_equal(st["s1form_half"].X, st["b0_trend"].X[half])
        np.testing.assert_array_equal(st["s1form_half"].coords, st["s1form"].coords[half])
        assert t2_build.rms_radius(st["s1form"].coords) == pytest.approx(target, rel=1e-5)
        assert t2_build.rms_radius(st["iso_half"].coords) == pytest.approx(target, rel=0.05)  # half of a sized copy
        assert not np.array_equal(st["iso_half"].coords, st["copy_half"].coords)
        assert not np.array_equal(st["s1form"].coords, st["copy"].coords)
    assert set(records) == {"copy_iso", "iso_half", "s1form", "s1form_half"}
    assert all(r["label"] == "oracle size" for arm in records.values() for r in arm.values())
    assert records["s1form"]["1"]["z_held"] and not records["iso_half"]["1"]["z_held"]
