"""experiments/t1_m1/gmom_combo.py (diagnostic N22 and try T26 of experiments/harness/TRIES.md) on synthetic data."""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pytest
from test_g_mom_verify import _small_inputs
from test_t1_state_transport import _gmom_ctx

from vec import harness, t1_scope
from vec import t1_state_transport as ts
from vec.boards import REPO_ROOT

SCRIPT = REPO_ROOT / "experiments" / "t1_m1" / "gmom_combo.py"
GROUPS = {"V-CM": "CM", "Endocardium": "Endocardium", "Blood": "Blood", "Foregut": "other:Foregut"}


@pytest.fixture
def gc(monkeypatch, tmp_path):
    spec = importlib.util.spec_from_file_location("gmom_combo_module", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    monkeypatch.setattr(mod, "RESULTS", tmp_path / "results")
    return mod


def _fit():
    """A small fold's cells, finished factors and R1 mask, as N22 uses them."""
    ref, start, _ = _small_inputs()
    factors = harness.g_mom_pair_factors(ref, "E9.5_RNA", start, "E8.5_RNA")[0]
    spec = t1_scope.load_spec()
    cls = harness._shallow(ref)
    t1_scope.classify(cls, source="E9.5_RNA", spec=spec)
    cats = cls.obs[t1_scope.CATEGORY_KEY].astype(str).to_numpy()
    pop = cls[np.isin(cats, list(factors.categories))].copy()
    pop.obs["source_obs"] = np.asarray(pop.obs_names).astype(str)
    if "counts" not in pop.layers:
        from vec.core import counts_of

        pop.layers["counts"] = counts_of(pop)
    labels = pop.obs[t1_scope.LABEL_KEY].astype(str).to_numpy()
    no_edit = np.isin(labels, sorted(ts.r1_unbracketed_labels(np.unique(labels))))
    return pop, no_edit, factors


def _direction(pop, size):
    rng = np.random.default_rng(3)
    return {g: rng.normal(0, size, pop.n_vars) for g in ("CM", "Endocardium")}


def test_n22_with_a_null_direction_reads_exactly_one(gc):
    pop, no_edit, factors = _fit()
    groups = {lab: GROUPS.get(lab, "x") for lab in pop.obs[t1_scope.LABEL_KEY].astype(str).unique()}
    moved = gc.moved_by_feng(pop, _direction(pop, 0.0), 0.2, groups, [], [1, 2])
    cat = "cm"
    r = gc.additivity_category(pop, moved, no_edit, factors, set(), cat)
    assert r["slope_without_m1"] == pytest.approx(r["slope_with_m1"])
    assert r["ratio"] == pytest.approx(1.0) and r["passes"]


def test_n22_a_real_direction_keeps_the_factors_slope_near_itself(gc):
    pop, no_edit, factors = _fit()
    groups = {lab: GROUPS.get(lab, "x") for lab in pop.obs[t1_scope.LABEL_KEY].astype(str).unique()}
    moved = gc.moved_by_feng(pop, _direction(pop, 0.3), 0.4, groups, [], [1, 2])
    for cat in ("cm", "endocardium"):
        r = gc.additivity_category(pop, moved, no_edit, factors, set(), cat)
        assert r["entries_in_the_regression"] > 0 and r["slope_without_m1"] > 0
        assert abs(r["ratio"] - 1) < 0.5


def test_moved_by_feng_leaves_the_uncovered_cells_and_the_row_order_alone(gc):
    pop, _, _ = _fit()
    labels = pop.obs[t1_scope.LABEL_KEY].astype(str).to_numpy()
    groups = {lab: GROUPS.get(lab, "x") for lab in np.unique(labels)}
    groups["Endocardium"] = "uncovered"
    moved = gc.moved_by_feng(pop, {"CM": _direction(pop, 0.5)["CM"]}, 1.0, groups, [], [1, 2])
    assert list(moved.obs_names) == list(pop.obs_names)
    other = labels != "V-CM"
    assert (moved.layers["counts"][other] != pop.layers["counts"][other]).nnz == 0
    assert (moved.layers["counts"][~other] != pop.layers["counts"][~other]).nnz > 0


def test_n22_reading_needs_every_rule_category_at_every_pair(gc):
    ok = {"passes": True}
    good = {p: {c: ok for c in gc.N22_CATEGORIES} for p in gc.N22_PAIRS}
    assert gc.n22_reading(good)["passes"]
    bad = {p: dict(v) for p, v in good.items()}
    bad["test"]["cm"] = {"passes": False}
    assert not gc.n22_reading(bad)["passes"]
    gone = {p: dict(v) for p, v in good.items()}
    del gone["backtest"]["endocardium"]
    r = gc.n22_reading(gone)
    assert not r["passes"] and r["missing"] == ["backtest:endocardium"]


def _arm(score, mmd, vario):
    return {"score": score, "groups": {"distribution": score / 10}, "metrics": {"mmd_u": mmd, "variogram": vario}}


def _rows(gc, key, share=1.0, drift=0.0, seeds=range(10)):
    t11 = gc.T11_RAW_GAIN[key]
    return [{"seed": s, "floor": 50.0,
             "arms": {"b0": _arm(50.0 + drift, 0.12, 0.008), "m1_feng": _arm(51.0, 0.12, 0.008),
                      "m1_r": _arm(51.0, 0.12, 0.008),
                      "m1_r_gmom": _arm(52.0, 0.12 - share * t11["mmd_u"], 0.008 - share * t11["variogram"])}}
            for s in seeds]


def _diag():
    return [{"seed": s, "scores": {"b0": 50.0, "m1_feng": 51.0}} for s in range(10)]


def _cuts(gc, share=1.0, drift=0.0):
    return {k: {"rows": _rows(gc, k, share, drift)} for k in gc.CUT_FILES}, {k: _diag() for k in gc.CUT_FILES}


def test_t26_reads_no_interference_at_t11s_gain_and_at_the_bar(gc):
    res, diag = _cuts(gc)
    assert gc.t26_reading(res, diag)["outcome"] == "NO INTERFERENCE"
    res, diag = _cuts(gc, share=0.51)
    assert gc.t26_reading(res, diag)["outcome"] == "NO INTERFERENCE"


def test_t26_below_half_of_t11s_gain_on_any_metric_or_cut_is_interference(gc):
    res, diag = _cuts(gc, share=0.4)
    assert gc.t26_reading(res, diag)["outcome"] == "INTERFERES"
    res, diag = _cuts(gc)
    for r in res[gc.TASK_WHOLE]["rows"]:
        r["arms"]["m1_r_gmom"]["metrics"]["variogram"] = r["arms"]["m1_r"]["metrics"]["variogram"]
    out = gc.t26_reading(res, diag)
    assert out["outcome"] == "INTERFERES" and out["cuts"][gc.TASK_WHOLE]["raw"]["variogram"]["below_bar"]
    assert not out["cuts"][gc.TASK_CHAMBER]["interferes"]


def test_t26_is_void_when_diag_v1_does_not_reproduce(gc):
    res, diag = _cuts(gc, drift=0.03)
    assert gc.t26_reading(res, diag)["outcome"] == "VOID"
    res, diag = _cuts(gc, drift=0.01)
    assert gc.t26_reading(res, diag)["outcome"] == "NO INTERFERENCE"


def test_t26_needs_both_cuts(gc):
    res, diag = _cuts(gc)
    del res[gc.TASK_WHOLE]
    with pytest.raises(ValueError, match="both cuts"):
        gc.t26_reading(res, diag)


def test_the_run_is_refused_without_a_passing_n22(gc, tmp_path):
    args = ["run", "--task", gc.TASK_CHAMBER, "--out", "t26", "--n22", "n22x"]
    with pytest.raises(SystemExit, match="refused"):
        gc.main(args)
    d = gc.RESULTS / "n22x"
    d.mkdir(parents=True)
    (d / "n22.json").write_text(json.dumps({"reading": {"passes": False}}))
    with pytest.raises(SystemExit, match="refused"):
        gc.main(args)
    (d / "n22.json").write_text(json.dumps({"reading": {"passes": True}}))
    assert gc.n22_passed(Path(d / "n22.json"))


def _dirs(ref):
    rng = np.random.default_rng(5)
    return {g: rng.normal(0, 0.4, ref.n_vars) for g in ("CM", "Endocardium")}


def test_the_factors_are_the_last_edit_on_m1_and_leave_r1s_cells_alone(gc):
    harness._G_MOM_PREP.clear()
    ref, start, _ = _small_inputs()
    ctx = _gmom_ctx(ref, start)
    groups = {lab: GROUPS.get(lab, "x") for lab in ref.obs[t1_scope.LABEL_KEY].astype(str).unique()}
    held = list(ctx.derived["t1_held_flat"])
    arms = gc.m1_r_arms(ctx, _dirs(ref), groups, 0.2, held, held)
    m1, on = arms["m1_r"], arms["m1_r_gmom"]
    assert m1.X.shape == on.X.shape and not np.array_equal(m1.X, on.X)
    twin = gc.m1_r_arms(ctx, _dirs(ref), groups, 0.2, held, held)
    assert np.array_equal(twin["m1_r_gmom"].X, on.X) and np.array_equal(twin["m1_r"].X, m1.X)
    same = (m1.X == on.X).all(axis=1)
    assert 0 < same.sum() < len(same)  # blood and non-chamber cells keep every count
