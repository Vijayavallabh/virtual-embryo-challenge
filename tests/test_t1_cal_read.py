"""experiments/t1_m1/cal_read.py (try T48 of experiments/harness/TRIES.md) on synthetic data: the calibration rule, the
arms and the reading. No released file, statistics file or score is read."""

from __future__ import annotations

import builtins
import importlib.util
import io
import json
import re
import sys

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import pytest
from scipy import sparse
from test_t1_step_read import _arms_world

from vec import t1_scope
from vec.boards import REPO_ROOT

SCRIPT = REPO_ROOT / "experiments" / "t1_m1" / "cal_read.py"
EPS, MIN_CELLS = 0.5, 20
GROUPS = {"V-CM": "CM", "A-CM": "CM", "Endocardium": "Endocardium", "Foregut": "other:Foregut"}
_SKIP_DIRS = {"boards", "data", "runs", "notes", "research", "reference", ".venv", ".venv-gpu", ".git", ".claude",
              "results"}


@pytest.fixture
def cr(monkeypatch, tmp_path):
    spec = importlib.util.spec_from_file_location("cal_read_module", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    monkeypatch.setattr(mod, "RESULTS", tmp_path / "results")
    return mod


# ---------------------------------------------------------------------------------------------------------------------
# A small synthetic pool: sparse counts, a known direction per group


def pool_world(seed=3, G=400, sizes=(("V-CM", 260), ("A-CM", 90), ("Endocardium", 120), ("Foregut", 80))):
    rng = np.random.default_rng(seed)
    genes = [f"g{i}" for i in range(G)]
    labels = [lab for lab, n in sizes for _ in range(n)]
    n = len(labels)
    lam = rng.gamma(0.4, 2.0, G)  # sparse: most genes below one count per cell
    C = rng.poisson(lam, size=(n, G)).astype(np.int64)
    C[:, 0] = 1  # a count of one in every cell: counts are recoverable from .X
    a = ad.AnnData(X=t1_scope.log1p_cp10k(sparse.csr_matrix(C)),
                   obs=pd.DataFrame({"celltype": labels}, index=[f"c{i}_1" for i in range(n)]),
                   var=pd.DataFrame(index=genes))
    dirs = {}
    for g in ("CM", "Endocardium"):
        d = rng.normal(0, 0.8, G)
        d[0] = 0.0
        d[rng.random(G) < 0.05] = np.nan
        dirs[g] = d
    held = ["g1", "g2", "g3"]
    return a, dirs, held


def _direct_share(cr, pool, dirs, held, f, seed=(48, 48)):
    """share_read's code by hand: external_shift on the covered cells, then per group share_stats over the moved genes
    with dense CP10k (share_read.calibration's form), combined by the groups' cells."""
    shr = cr.shr
    labels = pool.obs["celltype"].astype(str).to_numpy()
    cov = np.array([GROUPS.get(lab) in dirs for lab in labels])
    sub = pool[np.flatnonzero(cov)].copy()
    moved, _ = t1_scope.external_shift(sub, dirs, f, GROUPS, held, label_key="celltype", seed=list(seed))
    sl = labels[cov]
    not_held = ~np.isin(np.asarray(pool.var_names), held)
    shares, weights = {}, {}
    for g in dirs:
        m = np.array([GROUPS[lab] == g for lab in sl])
        weights[g] = float(m.sum())

        def cp(a, m=m):
            return {"cp": np.expm1(np.asarray(a.X[np.flatnonzero(m)].toarray(), np.float64)).mean(0)}

        st = shr.share_stats(shr.dl.lfc_ratio(cp(sub), cp(moved), EPS), dirs[g], np.isfinite(dirs[g]) & not_held,
                             None)
        shares[g] = st["share"]
    return shr.combine(shares, weights)[0]


def _cal(cr, pool, dirs, held, target, fraction=0.2, **kw):
    return cr.calibrated_fraction(pool, dirs, GROUPS, held, fraction, target, eps=EPS, min_cells=MIN_CELLS, **kw)


# ---------------------------------------------------------------------------------------------------------------------
# The rule


def test_delivered_is_share_reads_reading_of_external_shift_at_f(cr):
    pool, dirs, held = pool_world()
    for f in (0.2, 0.37, 1.0):
        got = cr.delivered(pool, dirs, GROUPS, held, f, eps=EPS, min_cells=MIN_CELLS)
        assert got == pytest.approx(_direct_share(cr, pool, dirs, held, f), rel=1e-10, abs=1e-12)
    # sparse counts and the pseudocount attenuate the step: it delivers less than its nominal size
    assert 0 < cr.delivered(pool, dirs, GROUPS, held, 0.2, eps=EPS, min_cells=MIN_CELLS) < 0.2


def test_delivered_rises_with_the_fraction_on_synthetic_data(cr):
    pool, dirs, held = pool_world()
    prep = cr.prepare(pool, dirs, GROUPS, held, EPS, MIN_CELLS)
    grid = np.linspace(0.2, 1.0, 17)
    d = [cr.delivered_on(prep, f)["share"] for f in grid]
    assert all(np.diff(d) > 0), d


def test_the_bisection_reaches_the_target_is_deterministic_and_counts_its_steps(cr):
    pool, dirs, held = pool_world()
    d02 = cr.delivered(pool, dirs, GROUPS, held, 0.2, eps=EPS, min_cells=MIN_CELLS)
    d10 = cr.delivered(pool, dirs, GROUPS, held, 1.0, eps=EPS, min_cells=MIN_CELLS)
    target = d02 + 0.4 * (d10 - d02)
    r = _cal(cr, pool, dirs, held, target)
    assert set(r) >= {"fraction", "delivered_at_base", "delivered", "steps"}
    assert 0.2 < r["fraction"] < 1.0 and r["steps"] >= 1
    assert abs(r["delivered"] - target) <= 1e-3 * abs(target)
    assert r["delivered_at_base"] == d02
    assert r["delivered"] == cr.delivered(pool, dirs, GROUPS, held, r["fraction"], eps=EPS, min_cells=MIN_CELLS)
    assert r["stream"] == [48, 48] and r["bracket"] == [0.2, 1.0]
    again = _cal(cr, pool, dirs, held, target)
    assert again == r


def test_a_target_equal_to_delivered_at_the_fraction_returns_the_fraction_itself(cr):
    pool, dirs, held = pool_world()
    d02 = cr.delivered(pool, dirs, GROUPS, held, 0.2, eps=EPS, min_cells=MIN_CELLS)
    r = _cal(cr, pool, dirs, held, d02)
    assert r["fraction"] == 0.2 and r["steps"] == 0 and r["delivered"] == d02


def test_the_rule_refuses_a_target_outside_its_bracket(cr):
    pool, dirs, held = pool_world()
    d02 = cr.delivered(pool, dirs, GROUPS, held, 0.2, eps=EPS, min_cells=MIN_CELLS)
    d10 = cr.delivered(pool, dirs, GROUPS, held, 1.0, eps=EPS, min_cells=MIN_CELLS)
    with pytest.raises(ValueError, match="outside the bracket"):
        _cal(cr, pool, dirs, held, d10 * 1.01)
    with pytest.raises(ValueError, match="outside the bracket"):
        _cal(cr, pool, dirs, held, d02 * 0.99)
    # a target at the top within the tolerance is inside: the answer is the upper end (which voids the read)
    r = _cal(cr, pool, dirs, held, d10)
    assert r["fraction"] == 1.0 and r["delivered"] == d10 and r["steps"] == 0


def test_the_rule_refuses_a_delivered_share_that_does_not_rise(cr):
    pool, dirs, held = pool_world()
    flat = {g: np.zeros_like(d) for g, d in dirs.items()}
    flat["CM"][5] = 1e-12
    with pytest.raises(ValueError, match="does not rise|not finite"):
        _cal(cr, pool, flat, held, 0.1)


def test_the_rule_refuses_when_a_jump_in_delivered_hides_the_target(cr, monkeypatch):
    pool, dirs, held = pool_world()
    monkeypatch.setattr(cr, "delivered_on", lambda prep, f, seed=None: {"share": 0.1 if f < 0.5 else 0.3})
    with pytest.raises(ValueError, match="did not reach"):
        _cal(cr, pool, dirs, held, 0.2)


def test_the_rule_reads_no_file_but_the_pool_it_is_given(cr, monkeypatch):
    pool, dirs, held = pool_world()
    d02 = cr.delivered(pool, dirs, GROUPS, held, 0.2, eps=EPS, min_cells=MIN_CELLS)
    opened = []

    def spy(*a, **k):
        opened.append(a[:1])
        raise AssertionError(f"the rule opened {a[:1]}")

    for mod, name in ((builtins, "open"), (io, "open"), (h5py, "File"), (ad, "read_h5ad"), (cr, "load_stage"),
                      (cr, "read_stage")):
        monkeypatch.setattr(mod, name, spy)
    monkeypatch.setattr(cr.Path, "read_text", spy)
    monkeypatch.setattr(cr.Path, "open", spy)
    r = _cal(cr, pool, dirs, held, d02 * 1.3)
    assert np.isfinite(r["fraction"]) and opened == []


def test_the_rule_uses_its_own_stream_which_no_build_step_uses(cr, monkeypatch):
    """Every stream index written as 'seed, <n>]' or as a '<NAME>_STREAM = <n>' constant in src/ and experiments/
    (*.py, cal_read.py left out): the calibration's index is none of them, and every external_shift call of the rule
    draws on [CAL_ROOT, CAL_STREAM]."""
    used, files = set(), 0
    for top in ("src", "experiments"):
        for path in (REPO_ROOT / top).rglob("*.py"):
            rel = path.relative_to(REPO_ROOT)
            if path.name == "cal_read.py" or _SKIP_DIRS & set(rel.parts[:-1]):
                continue
            files += 1
            text = path.read_text(errors="replace")
            used |= {int(n) for n in re.findall(r"seed, ?(\d+)\]", text)}
            used |= {int(n) for n in re.findall(r"_STREAM = (\d+)\b", text)}
    assert files > 10 and {1, 2, 4, 5, 15, 46} <= used, (files, sorted(used))
    assert cr.CAL_STREAM not in used, sorted(used)
    seeds = []
    real = t1_scope.external_shift

    def spy(*a, **k):
        seeds.append(k["seed"])
        return real(*a, **k)

    monkeypatch.setattr(cr.t1_scope, "external_shift", spy)
    pool, dirs, held = pool_world()
    d02 = cr.delivered(pool, dirs, GROUPS, held, 0.2, eps=EPS, min_cells=MIN_CELLS)
    _cal(cr, pool, dirs, held, d02 * 1.3)
    assert len(seeds) >= 3 and all(s == [cr.CAL_ROOT, cr.CAL_STREAM] for s in seeds)


def test_calibrate_records_both_targets(cr):
    pool, dirs, held = pool_world()
    L = 0.834
    rec = cr.calibrate(pool, dirs, GROUPS, held, 0.2, L, "f" * 64, eps=EPS, min_cells=MIN_CELLS)
    assert rec["targets"] == {"m1_r_gmom_calL": 0.2 * L, "m1_r_gmom_cal": 0.2}
    for a, t in rec["targets"].items():
        assert abs(rec["delivered"][a] - t) <= 1e-3 * t
    assert rec["fractions"]["m1_r_gmom_calL"] < rec["fractions"]["m1_r_gmom_cal"]
    assert rec["L"] == L and rec["t47_sha256"] == "f" * 64 and np.isfinite(rec["delivered_at_base"])


# ---------------------------------------------------------------------------------------------------------------------
# The arms


CAL_FR = {"m1_r_gmom_calL": 0.31, "m1_r_gmom_cal": 0.37}


def test_the_calibrated_arms_differ_from_the_base_build_in_the_fraction_only(cr, monkeypatch):
    calls = []

    def spy(ctx, directions, groups, fraction, held_m1, held_factors):
        calls.append({"ctx": ctx, "directions": directions, "groups": groups, "fraction": fraction,
                      "held_m1": held_m1, "held_factors": held_factors})
        return {"m1_r": ("m1_r", fraction), "m1_r_gmom": ("m1_r_gmom", fraction)}

    monkeypatch.setattr(cr.gc, "m1_r_arms", spy)
    ctx, dirs, groups, held_m1, held_f = object(), {"CM": np.zeros(3)}, {"V-CM": "CM"}, ["a", "b"], ["a"]
    arms = cr.cal_arms(ctx, dirs, groups, 0.2, held_m1, held_f, CAL_FR)
    assert set(arms) == {"m1_r_gmom", "m1_r_gmom_x15", "m1_r_gmom_x2", "m1_r_gmom_calL", "m1_r_gmom_cal"}
    assert [c["fraction"] for c in calls] == [0.2, 0.2 * 1.5, 0.2 * 2, 0.31, 0.37]
    base = calls[0]
    for c in calls:
        for key in ("ctx", "directions", "groups", "held_m1", "held_factors"):
            assert c[key] is base[key], f"an arm differs from the base build in {key}"
    assert arms["m1_r_gmom_calL"] == ("m1_r_gmom", 0.31) and arms["m1_r_gmom_cal"] == ("m1_r_gmom", 0.37)
    assert cr.arm_fractions(0.2, CAL_FR) == {"m1_r_gmom": 0.2, "m1_r_gmom_x15": 0.2 * 1.5, "m1_r_gmom_x2": 0.2 * 2,
                                             **CAL_FR}


def test_a_calibrated_arm_at_the_base_fraction_is_the_base_bit_for_bit(cr):
    """The idle case: a calibration that returns the fraction itself builds m1_r_gmom bit for bit."""
    ctx, groups, dirs, held = _arms_world()
    arms = cr.cal_arms(ctx, dirs, groups, 0.2, held, held, {"m1_r_gmom_calL": 0.2, "m1_r_gmom_cal": 0.33})
    assert arms["m1_r_gmom_calL"].X.tobytes() == arms["m1_r_gmom"].X.tobytes()
    assert not np.array_equal(arms["m1_r_gmom_cal"].X, arms["m1_r_gmom"].X)


# ---------------------------------------------------------------------------------------------------------------------
# The reading


L = 0.834


def _arm(score, mmd_u=0.1):
    return {"score": score, "groups": {"distribution": score / 100},
            "metrics": {"mmd_u": mmd_u, "variogram": 0.01},
            "skills": {"de_score": score / 200, "de_direction": 0.4, "mmd_u": 0.3, "variogram": 0.2}}


B0, BASE = 50.0, 53.0
X15, X2 = 0.9, 1.8


def _cal_record(fractions=None, delivered=None, L_=L):
    fr = fractions or {"m1_r_gmom_calL": 0.29, "m1_r_gmom_cal": 0.35}
    dv = delivered or {"m1_r_gmom_calL": 0.2 * L_, "m1_r_gmom_cal": 0.2}
    return {"L": L_, "t47": "x", "t47_sha256": "a" * 64, "base_fraction": 0.2,
            "targets": {"m1_r_gmom_calL": 0.2 * L_, "m1_r_gmom_cal": 0.2}, "fractions": fr,
            "delivered_at_base": 0.117, "delivered": dv, "steps": {"m1_r_gmom_calL": 9, "m1_r_gmom_cal": 10},
            "stream": [48, 48], "bracket": [0.2, 1.0], "rel_tol": 1e-3, "seconds": 1.0}


def _task(gains, seeds=range(10), commit="c" * 40, clean=True, task=None, fraction=0.2, fractions=None,
          base_shift=0.0, x_shift=0.0, identical=False, per_seed=None, cal=None):
    """A task file. gains: calibrated arm -> mean gain over m1_r_gmom; per_seed: arm -> ten per-seed gains. A
    calibrated arm's raw metric differs from the base's by 1e-9 unless `identical`."""
    cal = cal if cal is not None else _cal_record()
    rows = []
    for i, s in enumerate(seeds):
        arms = {"b0": _arm(B0), "m1_r_gmom": _arm(BASE + base_shift), "m1_r_gmom_x15": _arm(BASE + X15 + x_shift),
                "m1_r_gmom_x2": _arm(BASE + X2)}
        for a in ("m1_r_gmom_calL", "m1_r_gmom_cal"):
            g = per_seed[a][i] if per_seed and a in per_seed else gains.get(a, 0.0)
            arms[a] = _arm(BASE + base_shift + g, mmd_u=0.1 if identical else 0.1 + 1e-9)
        rows.append({"seed": s, "floor": 50.0, "arms": arms})
    fr = fractions or {"m1_r_gmom": fraction, "m1_r_gmom_x15": fraction * 1.5, "m1_r_gmom_x2": fraction * 2,
                       **(cal.get("fractions") or {})}
    return {"label": "synthetic", "task": task, "fraction": fraction, "fractions": fr, "calibration": cal,
            "seeds": list(seeds), "rows": rows, "stamp": {"commit": commit, "tree_clean": clean, "utc": "x",
                                                          "inputs": {}}}


def _write(cr, tmp_path, cuts):
    """cuts: cut key -> task record. T42's and T44's committed files, synthetic, beside them."""
    out = tmp_path / "out"
    (out / "tasks").mkdir(parents=True)
    for k, rec in cuts.items():
        rec = {**rec, "task": k if rec.get("task") is None else rec["task"]}
        (out / "tasks" / cr.CUT_FILES[k]).write_text(json.dumps(rec))
    t42, t44 = tmp_path / "t42", tmp_path / "t44"
    t42.mkdir()
    t44.mkdir()
    for f in cr.CUT_FILES.values():
        rows = [{"seed": s, "arms": {"b0": _arm(B0), "m1_r_gmom": _arm(BASE), "m1_r_gmom_js": _arm(BASE - 0.13)}}
                for s in range(10)]
        (t42 / f).write_text(json.dumps({"rows": rows}))
        rows = [{"seed": s, "arms": {"b0": _arm(B0), "m1_r_gmom": _arm(BASE), "m1_r_gmom_x15": _arm(BASE + X15),
                                     "m1_r_gmom_x2": _arm(BASE + X2)}} for s in range(10)]
        (t44 / f).write_text(json.dumps({"rows": rows}))
    return out, t42, t44


def _read(cr, tmp_path, gains=None, per_seed=None, **over):
    gains = gains or {}
    cuts = {k: _task(gains, per_seed=per_seed) for k in cr.CUT_FILES}
    for k, kw in over.items():
        if kw is None:
            del cuts[getattr(cr.gc, k)]
        else:
            cuts[getattr(cr.gc, k)] = _task(gains, **kw)
    out, t42, t44 = _write(cr, tmp_path, cuts)
    return cr.read(out, t42, t44)["reading"]


@pytest.mark.parametrize("gains, outcome, t48b", [
    ({"m1_r_gmom_calL": 1.0}, "the calibrated step passes on both cuts", False),
    ({"m1_r_gmom_calL": 1.0, "m1_r_gmom_cal": 1.2}, "the calibrated step passes on both cuts", True),
    ({"m1_r_gmom_calL": -0.1, "m1_r_gmom_cal": 1.2}, "the calibrated step loses", True),
    ({"m1_r_gmom_calL": 0.1}, "not separated", False),
    ({}, "not separated", False),
])
def test_each_statement_is_reached_by_a_table_built_for_it(cr, tmp_path, gains, outcome, t48b):
    r = _read(cr, tmp_path, gains)
    assert r["valid"] and r["outcome"] == outcome
    assert r["t48b_test"]["passes"] is t48b and r["t48b_test"]["decides"] is False
    for a in cr.CAL_ARMS:
        want = gains.get(a, 0.0)
        for other, off in (("m1_r_gmom", 0.0), ("m1_r_gmom_x15", X15), ("m1_r_gmom_x2", X2)):
            comp = r["comparisons"][a][other]
            assert comp["base"] == other and comp["change"] == a
            cut = comp["cuts"][cr.gc.TASK_CHAMBER]
            assert cut["mean"] == pytest.approx(want - off)
            assert cut["range"] == pytest.approx([want - off, want - off])
            assert set(cut["skills"]) == {"de_score", "de_direction", "mmd_u", "variogram"}
            assert cut["points_by_metric"]["de_score"] == pytest.approx(100 * 0.25 * (want - off) / 200)
            assert cut["points_by_metric"]["de_terms"] == pytest.approx(100 * 0.25 * (want - off) / 200)
            assert cut["points_by_metric"]["mmd_u"] == pytest.approx(0.0)
        assert r["mean_scores"][cr.gc.TASK_WHOLE][a] == pytest.approx(BASE + want)
    assert r["calibration_same_on_both_cuts"] is True


def test_losing_needs_at_most_two_seeds_ahead_on_each_cut(cr, tmp_path):
    two_ahead = [0.05, 0.05] + [-0.3] * 8
    three_ahead = [0.05, 0.05, 0.05] + [-0.3] * 7
    r = _read(cr, tmp_path, per_seed={"m1_r_gmom_calL": two_ahead})
    assert r["valid"] and r["outcome"] == cr.LOSES
    r = _read(cr, tmp_path / "more", {"m1_r_gmom_calL": -0.1},
              TASK_WHOLE={"per_seed": {"m1_r_gmom_calL": three_ahead}})
    assert r["valid"] and r["outcome"] == cr.NOT_SEPARATED


def test_a_gain_of_exactly_point_three_does_not_pass(cr, tmp_path):
    r = _read(cr, tmp_path, {"m1_r_gmom_calL": 0.3, "m1_r_gmom_cal": 0.3})
    assert r["valid"] and not r["comparisons"]["m1_r_gmom_calL"]["m1_r_gmom"]["passes"]
    assert r["outcome"] == cr.NOT_SEPARATED and r["t48b_test"]["passes"] is False
    r = _read(cr, tmp_path / "more", {"m1_r_gmom_calL": 0.31})
    assert r["outcome"] == cr.PASSES


def test_a_loss_beyond_point_three_on_one_cut_vetoes(cr, tmp_path):
    cuts = {cr.gc.TASK_CHAMBER: _task({"m1_r_gmom_calL": 1.0}), cr.gc.TASK_WHOLE: _task({"m1_r_gmom_calL": -0.4})}
    out, t42, t44 = _write(cr, tmp_path, cuts)
    r = cr.read(out, t42, t44)["reading"]
    assert r["valid"] and r["outcome"] == cr.NOT_SEPARATED


@pytest.mark.parametrize("over, reason", [
    ({"TASK_WHOLE": None}, "missing cuts"),
    ({"TASK_CHAMBER": {"seeds": range(9)}}, "seeds"),
    ({"TASK_WHOLE": {"base_shift": 1e-6}}, "do not reproduce T42"),
    ({"TASK_WHOLE": {"x_shift": 1e-6}}, "do not reproduce T44"),
    ({"TASK_WHOLE": {"fraction": 0.4}}, "fraction 0.4"),
    ({"TASK_WHOLE": {"fractions": {"m1_r_gmom": 0.2, "m1_r_gmom_x15": 0.3, "m1_r_gmom_x2": 0.4,
                                   "m1_r_gmom_calL": 0.29, "m1_r_gmom_cal": 0.36}}}, "wrong fraction"),
    ({"TASK_WHOLE": {"fractions": {"m1_r_gmom": 0.2, "m1_r_gmom_x15": 0.3, "m1_r_gmom_x2": 0.6,
                                   "m1_r_gmom_calL": 0.29, "m1_r_gmom_cal": 0.35}}}, "wrong fraction"),
    ({"TASK_CHAMBER": {"commit": "d" * 40}}, "2 commits"),
    ({"TASK_WHOLE": {"clean": False}}, "not clean"),
    ({"TASK_CHAMBER": {"task": "T1/E9.5->E10.5"}}, "the file holds task"),
    ({"TASK_WHOLE": {"cal": _cal_record(fractions={"m1_r_gmom_calL": 0.29, "m1_r_gmom_cal": 1.0})}}, "outside"),
    ({"TASK_WHOLE": {"cal": _cal_record(fractions={"m1_r_gmom_calL": 0.2, "m1_r_gmom_cal": 0.35})}}, "outside"),
    ({"TASK_WHOLE": {"cal": _cal_record(fractions={"m1_r_gmom_calL": 0.36, "m1_r_gmom_cal": 0.35})}},
     "not the smaller"),
    ({"TASK_WHOLE": {"cal": _cal_record(delivered={"m1_r_gmom_calL": 0.2 * L * 1.0011, "m1_r_gmom_cal": 0.2})}},
     "off its target"),
    ({"TASK_CHAMBER": {"cal": _cal_record(delivered={"m1_r_gmom_calL": 0.2 * L, "m1_r_gmom_cal": 0.1997})}},
     "off its target"),
    ({"TASK_CHAMBER": {"cal": {}}}, "lacks a fraction"),
    ({"TASK_CHAMBER": {"identical": True}, "TASK_WHOLE": {"identical": True}}, "exactly as m1_r_gmom"),
])
def test_the_void_cases(cr, tmp_path, over, reason):
    r = _read(cr, tmp_path, {"m1_r_gmom_calL": 1.0}, **over)
    assert not r["valid"] and r["outcome"] == "VOID"
    assert any(reason in x for x in r["validity"]["reasons"]), r["validity"]["reasons"]


def test_a_delivered_share_inside_the_tolerance_and_one_identical_cut_stay_valid(cr, tmp_path):
    cal = _cal_record(delivered={"m1_r_gmom_calL": 0.2 * L * 1.0009, "m1_r_gmom_cal": 0.2 * 0.9991})
    r = _read(cr, tmp_path, {"m1_r_gmom_calL": 1.0}, TASK_WHOLE={"identical": True, "cal": cal})
    assert r["valid"], r["validity"]["reasons"]
    assert r["calibration_same_on_both_cuts"] is False


def test_a_rescored_arm_off_by_one_e_minus_six_voids_b0_too(cr, tmp_path):
    cuts = {k: _task({"m1_r_gmom_calL": 1.0}) for k in cr.CUT_FILES}
    cuts[cr.gc.TASK_CHAMBER]["rows"][0]["arms"]["b0"]["score"] += 1e-6
    out, t42, t44 = _write(cr, tmp_path, cuts)
    r = cr.read(out, t42, t44)["reading"]
    worst = r["validity"]["reproduction"][cr.gc.TASK_CHAMBER]["t42"]["largest_difference"]
    assert not r["valid"] and worst["b0"] > 1e-9
    r = _read(cr, tmp_path / "ok", {"m1_r_gmom_calL": 1.0}, TASK_WHOLE={"base_shift": 1e-12, "x_shift": 1e-12})
    assert r["valid"]


def test_the_reading_refuses_to_overwrite(cr, tmp_path):
    out, t42, t44 = _write(cr, tmp_path, {k: _task({}) for k in cr.CUT_FILES})
    cr.read(out, t42, t44)
    assert (out / "t48.json").exists() and (out / "t48.txt").exists()
    with pytest.raises(SystemExit, match="already holds"):
        cr.read(out, t42, t44)


# ---------------------------------------------------------------------------------------------------------------------
# The run's refusals


def _pass_gates(cr, monkeypatch):
    monkeypatch.setattr(cr.n12, "tries_rows", lambda: {"T48"})
    monkeypatch.setattr(cr.gc, "n22_passed", lambda path: True)


def test_the_run_needs_its_row(cr, monkeypatch):
    monkeypatch.setattr(cr.n12, "tries_rows", lambda: {"T44", "T47"})
    monkeypatch.setattr(cr.gc, "n22_passed", lambda path: True)
    with pytest.raises(SystemExit, match="T48 has no row"):
        cr.main(["run", "--task", cr.gc.TASK_CHAMBER, "--out", "t48"])


def test_the_run_needs_n22(cr, monkeypatch):
    monkeypatch.setattr(cr.n12, "tries_rows", lambda: {"T48"})
    with pytest.raises(SystemExit, match="does not hold a passing N22"):
        cr.main(["run", "--task", cr.gc.TASK_CHAMBER, "--out", "t48", "--n22", "n22x"])


def test_the_run_refuses_an_existing_task_file(cr, monkeypatch):
    _pass_gates(cr, monkeypatch)
    tasks = cr.RESULTS / "t48" / "tasks"
    tasks.mkdir(parents=True)
    (tasks / cr.CUT_FILES[cr.gc.TASK_CHAMBER]).write_text("{}")
    with pytest.raises(SystemExit, match="already holds"):
        cr.main(["run", "--task", cr.gc.TASK_CHAMBER, "--out", "t48"])


def _t47(L_value=L, clean=True):
    return {"read": {"reading": {"L": {"value": L_value, "p5": 0.8, "p95": 0.9, "nonfinite_draws": 0}}},
            "stamp": {"commit": "c" * 40, "tree_clean": clean}}


@pytest.mark.parametrize("content, match", [
    (None, "is missing"),
    ("{not json", "is void"),
    (json.dumps({"read": {}}), "is void"),
    (json.dumps(_t47(clean=False)), "is void"),
    (json.dumps(_t47(float("nan"))), "not finite and in"),
    (json.dumps(_t47(1.5)), "not finite and in"),
    (json.dumps(_t47(0.5)), "not finite and in"),
    (json.dumps(_t47("0.8")), "not finite and in"),
])
def test_the_run_refuses_a_missing_or_void_t47(cr, monkeypatch, tmp_path, content, match):
    _pass_gates(cr, monkeypatch)
    path = tmp_path / "t47.json"
    if content is not None:
        path.write_text(content)
    monkeypatch.setattr(cr, "T47_JSON", path)
    loaded = []
    monkeypatch.setattr(cr.harness, "load_tasks", lambda *a, **k: loaded.append(1))
    with pytest.raises(SystemExit, match=match):
        cr.main(["run", "--task", cr.gc.TASK_CHAMBER, "--out", "t48"])
    assert loaded == []  # refused before any data is touched


def test_t47s_raw_L_and_sha256_are_read(cr, tmp_path):
    path = tmp_path / "t47.json"
    path.write_text(json.dumps(_t47(0.834)))
    L_, sha = cr.t47_scale(path)
    assert L_ == 0.834 and sha == cr.gv.sha256_file(path)
