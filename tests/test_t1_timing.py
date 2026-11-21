"""Try T67 (experiments/harness/TRIES.md): vec.t1_timing, experiments/t1_m1/timing_read.py, the late labelling
experiments/external/feng_relabel_late.py and the launcher experiments/t1_m1/timing.sh, on synthetic arrays, synthetic
statistics, synthetic stages and hand-made task tables only. No released file, statistics file, board file or score is
read; the launcher's refusal paths run in a throwaway git repository under a temporary home."""

from __future__ import annotations

import importlib.util
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
from test_t1_scope import HELD
from test_t1_state_transport import _gmom_ctx, classified_pair

from vec import t1_timing as tt
from vec.core import counts_of
from vec.t1_scope import load_spec

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "experiments" / "t1_m1" / "timing_read.py"
LATE = ROOT / "experiments" / "external" / "feng_relabel_late.py"
LAUNCHER = ROOT / "experiments" / "t1_m1" / "timing.sh"
GROUPS = ("CM", "Endocardium", "Endothelium", "Epicardium")
EPS = 0.5


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture
def tr(monkeypatch, tmp_path):
    mod = _load("timing_read_module", SCRIPT)
    monkeypatch.setattr(mod, "RESULTS", tmp_path / "results")
    return mod


@pytest.fixture
def frl():
    return _load("feng_relabel_late_module", LATE)


def _d(n=400, seed=1, nan_every=13):
    d = np.random.default_rng(seed).normal(0.0, 0.5, n)
    d[::nan_every] = np.nan
    return d


# ---------------------------------------------------------------------------------------------------------------------
# vec.t1_timing: the arithmetic


def test_the_hermite_coefficients_are_the_cubic_hermites_increment_over_the_test_span():
    def h(t, d, m0, m1):  # cubic Hermite on tau in [0, 1], slopes per day times the bracket's five days
        return ((3 * t**2 - 2 * t**3) * d + (t**3 - 2 * t**2 + t) * 5 * m0 + (t**3 - t**2) * 5 * m1)

    rng = np.random.default_rng(0)
    d, m1 = rng.normal(0, 1, 50), rng.normal(0, 0.2, 50)
    want = h(0.6, d, d / 5, m1) - h(0.2, d, d / 5, m1)
    assert np.allclose(tt.hermite_increment(d, m1), want, atol=1e-12)
    assert np.allclose(tt.hermite_increment(d, d / 5), 0.4 * d, atol=1e-12)  # a steady-rate gene: 0.4 D
    assert np.allclose(tt.hermite_increment(d, 0 * d), 0.512 * d)  # levelled off by E14.5
    assert np.allclose(tt.timing_shape(d, m1) * tt.TEST_FRACTION, tt.hermite_increment(d, m1), atol=1e-12)


def test_the_steady_rate_identity_is_bit_for_bit_with_nan_kept():
    for seed in range(5):
        d = _d(seed=seed)
        d[5] = -0.0
        out = tt.timing_direction(d, tt.steady_m1(d))
        assert out.tobytes() == d.tobytes() and tt.steady_identity(d)
    m1, fell = tt.timing_m1(_d(), np.full(400, 9.0), np.full(400, 9.0), np.full(400, 9.0), group_ok=False)
    assert tt.timing_direction(_d(), m1).tobytes() == _d().tobytes() and fell[np.isfinite(_d())].all()


def test_the_ols_slope():
    rng = np.random.default_rng(2)
    y = rng.normal(0, 1, (4, 30))
    want = np.polyfit(np.array(tt.SLOPE_DAYS), y, 1)[0]
    got = tt.ols_slope(y)
    assert np.allclose(got, want)
    y[2, 7] = np.nan
    assert np.isnan(tt.ols_slope(y)[7]) and np.isfinite(tt.ols_slope(y)[8])
    with pytest.raises(ValueError, match="one row per day"):
        tt.ols_slope(y[:3])


def test_the_limiter():
    d = np.array([1.0, 1.0, -1.0, -1.0, 1.0, 0.0, np.nan, 1.0, 2.0])
    m1 = np.array([-0.1, 0.9, 0.1, -0.9, 0.3, 0.2, 0.1, np.nan, 0.0])
    got = tt.limit(m1, d)
    assert got[:6].tolist() == [0.0, 0.6, 0.0, -0.6, 0.3, 0.0]
    assert np.isnan(got[6]) and np.isnan(got[7]) and got[8] == 0.0


def test_the_split_half_and_group_fallbacks():
    d = np.array([1.0, 1.0, 1.0, 1.0, np.nan])
    m1 = np.array([0.3, 0.3, 0.3, 0.9, 0.3])
    a = np.array([0.3, -0.1, 0.0, 0.9, 0.3])
    b = np.array([0.3, 0.2, 0.31, 0.9, 0.3])
    got, fell = tt.timing_m1(d, m1, a, b, group_ok=True)
    # gene 1: halves of different sign; gene 2: halves 0.31 apart, more than |m1| 0.3; gene 3: limited, not fallen
    assert fell.tolist() == [False, True, True, False, False]
    assert got[:4].tolist() == [0.3, 0.2, 0.2, 0.6] and np.isnan(got[4])
    nan_half, fell = tt.timing_m1(d, m1, np.array([np.nan, 0.3, 0.3, 0.3, 0.3]), b, group_ok=True)
    assert fell[0] and nan_half[0] == 0.2
    got, fell = tt.timing_m1(d, m1, a, b, group_ok=False)
    assert fell.tolist() == [True, True, True, True, False] and got[:4].tolist() == [0.2] * 4
    late = {"E15.5": 60, "E16.5": 49, "E17.5": 80}
    r = tt.group_timing("CM", d, m1, a, b, np.ones(5), late)
    assert r["group_ok"] is False and r["fallback_share"] == 1.0 and r["idle_ratio"] == 0.0
    assert r["direction"].tobytes() == d.tobytes() and r["perm_direction"].tobytes() == d.tobytes()
    assert tt.group_timing("CM", d, m1, a, b, np.ones(5), {**late, "E16.5": 50})["group_ok"] is True
    with pytest.raises(ValueError, match="no late cell count"):
        tt.group_timing("CM", d, m1, a, b, np.ones(5), {"E15.5": 60})


def test_the_norm_match_and_the_idle_ratio():
    d = _d()
    rng = np.random.default_rng(4)
    m1 = tt.limit(d / 5 + rng.normal(0, 0.05, d.size), d)
    out = tt.timing_direction(d, m1)
    ok = np.isfinite(d)
    assert np.array_equal(np.isfinite(out), ok)
    assert tt.l2(out) == pytest.approx(tt.l2(d), rel=1e-12)
    # the ratio is ||w_scaled - 0.4 D|| / ||0.4 D||, w = 0.512 D - 0.56 m1 rescaled to ||0.4 D||
    w = tt.hermite_increment(d, m1)
    w_scaled = w * tt.l2(0.4 * d) / tt.l2(w)
    want = tt.l2((w_scaled - 0.4 * d)[ok]) / tt.l2(0.4 * d)
    assert tt.idle_ratio(out, d) == pytest.approx(want, rel=1e-9)
    # the practice arm: the same shape at ||0.2 D|| is this direction at the practice fraction 0.2
    assert np.allclose((0.2 * out)[ok], (w * tt.l2(0.2 * d) / tt.l2(w))[ok])
    assert tt.is_idle({"CM": 0.01, "Endothelium": 0.049}) and not tt.is_idle({"CM": 0.01, "Endothelium": 0.05})
    assert not tt.is_idle({"CM": float("nan")})
    with pytest.raises(ValueError):
        tt.is_idle({})
    with pytest.raises(ValueError, match="no size"):
        tt.norm_match(np.zeros(3), np.ones(3))


def test_the_permuted_control_keeps_each_cells_deviations_and_breaks_the_gene_map():
    n = 2000
    rng = np.random.default_rng(7)
    d = rng.normal(0, 1, n)
    d[::50] = np.nan
    m1 = d / 5 + rng.normal(0, 0.1, n)
    level = rng.gamma(1.0, 5.0, n)
    perm = tt.permuted_m1(d, m1, level, "Endothelium")
    ok = np.flatnonzero(np.isfinite(d))
    cell = tt.equal_count_bins(np.abs(d[ok]), 20) * 5 + tt.equal_count_bins(level[ok], 5)
    assert np.unique(cell).size == 100
    raw = np.full(n, np.nan)  # the permutation before the limiter, from the same seed
    rng2 = np.random.default_rng([tt.PERM_SEED, tt.GROUP_ORDER.index("Endothelium")])
    delta = m1 - d / 5
    raw[ok] = delta[ok]
    for c in np.unique(cell):
        members = ok[cell == c]
        raw[members] = delta[members[rng2.permutation(members.size)]]
        assert np.allclose(np.sort(raw[members]), np.sort(delta[members]))  # each cell keeps its values
    assert np.array_equal(perm, tt.limit(d / 5 + raw, d), equal_nan=True)
    assert np.isnan(perm[::50]).all()
    assert abs(np.corrcoef(raw[ok], delta[ok])[0, 1]) < 0.2  # not the measured gene order
    assert np.array_equal(perm, tt.permuted_m1(d, m1, level, "Endothelium"), equal_nan=True)  # fixed seed
    assert not np.array_equal(perm, tt.permuted_m1(d, m1, level, "CM"), equal_nan=True)  # seeded per group
    with pytest.raises(ValueError, match="permutation seed"):
        tt.permuted_m1(d, m1, level, "Blood")


# ---------------------------------------------------------------------------------------------------------------------
# The timing kit on synthetic statistics


def _stats(stages, labels, genes, level, cells):
    """Statistics in the clean reduction's layout ('stage|zone|label|half'): for each stage, label and half, `cells`
    cells whose mean CP10k is exp(level(stage, label, half)) - EPS on every gene (so log(mean + EPS) is level)."""
    keys, n, s_cp = [], [], []
    for stage in stages:
        for lab in labels:
            for h in ("A", "B"):
                c = cells(stage, lab, h)
                if c == 0:
                    continue
                keys.append(f"{stage}|LV|{lab}|{h}")
                n.append(float(c))
                s_cp.append((np.exp(level(stage, lab, h)) - EPS) * c)
    s_cp = np.vstack(s_cp)
    return {"keys": np.array(keys, dtype=object), "genes": np.array(genes, dtype=object), "n": np.array(n),
            "S_cp": s_cp, "S_log": np.zeros_like(s_cp), "S_nz": np.zeros_like(s_cp), "depth_sum": np.ones(len(n))}


DAY = {"E9.5": 9.5, "E14.5": 14.5, "E15.5": 15.5, "E16.5": 16.5, "E17.5": 17.5}
LABELS = ("r_cm", "r_endocardial_ec", "r_vascular_ec", "r_epicardial")
CFG = {"metrics": {"pseudocount_cp10k": EPS},
       "external_labels": {"CM": [["r_cm"]], "Endocardium": [["r_endocardial_ec"]],
                           "Endothelium": [["r_endocardial_ec", "r_vascular_ec"]], "Epicardium": [["r_epicardial"]],
                           "Blood": [None]}}


def test_the_kit_reads_the_slope_on_the_groups_labels(tr):
    t1_genes = [f"G{i}" for i in range(40)]
    feng = [*reversed(t1_genes[:35]), "G35", "G35", "X1"]  # G35 twice (dropped), G36-G39 absent
    k = len(feng)
    slope = {lab: np.linspace(-0.2, 0.3, k) * (i + 1) for i, lab in enumerate(LABELS)}
    base = np.linspace(1.0, 3.0, k)

    def level(stage, lab, h):
        return base + slope[lab] * (DAY[stage] - 14.5) + (0.01 if h == "A" else 0.0)

    st_e = _stats(("E9.5", "E14.5"), LABELS, feng, level, lambda s, lab, h: 40)
    st_l = _stats(tt.LATE_STAGES, LABELS, feng, level,
                  lambda s, lab, h: 10 if lab == "r_epicardial" else 40)
    dirs = {g: _d(40, seed=i) for i, g in enumerate(GROUPS)}
    fb = {"Endocardium": {"takes": "Endothelium"}, "CM": {"takes": "CM"}}
    kit = tr.timing_kit(t1_genes, dirs, fb, st_e, st_l, CFG)
    pos = {g: i for i, g in enumerate(feng)}
    for g, lab in (("CM", "r_cm"), ("Epicardium", "r_epicardial")):
        got = kit[g]["m1_raw"]
        want = np.array([slope[lab][pos[x]] if x in pos and x != "G35" else np.nan for x in t1_genes])
        assert np.allclose(got, want, equal_nan=True)
    # Endocardium took Endothelium's direction, so its slope is read on Endothelium's labels (both labels pooled)
    assert kit["Endocardium"]["source_group"] == "Endothelium"
    assert kit["Endocardium"]["labels"] == ["r_endocardial_ec", "r_vascular_ec"]
    assert kit["Endocardium"]["cells"] == {"E14.5": 160, "E15.5": 160, "E16.5": 160, "E17.5": 160}
    assert kit["CM"]["group_ok"] and not kit["Epicardium"]["group_ok"]  # 20 Epicardium cells per late stage
    assert kit["Epicardium"]["direction"].tobytes() == dirs["Epicardium"].tobytes()
    dd = tr.arm_directions(dirs, kit)
    assert set(dd) == {"m1rc_gmom", "m1rc_timing_gmom", "m1rc_timing_perm_gmom", "m1rc_timing_steady"}
    rec = tr.kit_record(kit, dd)
    assert rec["steady_identity"] is True and set(rec["groups"]) == set(GROUPS)
    assert rec["direction_sha256"]["m1rc_gmom"] == rec["direction_sha256"]["m1rc_timing_steady"]
    with pytest.raises(ValueError, match="R2's guard"):
        tr.timing_kit(t1_genes, {"Blood": dirs["CM"]}, {}, st_e, st_l, CFG)


# ---------------------------------------------------------------------------------------------------------------------
# The arms on synthetic stages


def _world():
    start, ref = classified_pair(seed=5, n_scale=0.5)
    return ref, start


def _arm_dirs(tr, n):
    dirs = {g: _d(n, seed=i + 20, nan_every=17) for i, g in enumerate(GROUPS)}
    rng = np.random.default_rng(9)
    kit = {}
    for g, d in dirs.items():
        m1 = tt.limit(d / 5 + rng.normal(0, 0.3, n), d)
        kit[g] = {"direction": tt.timing_direction(d, m1),
                  "perm_direction": tt.timing_direction(d, tt.permuted_m1(d, m1, rng.gamma(1, 5, n), g))}
    return tr.arm_directions(dirs, kit)


def _held_m1(genes):
    return sorted({*HELD, genes[3], genes[7]} & set(genes))


def _dense(x):
    return np.asarray(x.todense()) if hasattr(x, "todense") else np.asarray(x)


def test_the_cut_arms(tr):
    ref, start = _world()
    ctx = _gmom_ctx(ref, start, seed=1000003)
    fits = tr.br.gmom_fits(ref, "E9.5_RNA", start, "E8.5_RNA")
    held_m1 = _held_m1(list(ref.var_names))
    dirs = _arm_dirs(tr, ref.n_vars)
    arms = tr.cut_arms(ctx, fits["file"], dirs, 0.2, held_m1, steady=True)
    assert tuple(arms) == (*tr.ARMS, tr.STEADY)
    X = {a: np.asarray(v.X) for a, v in arms.items()}
    t62 = tr.rc.cut_pair(ctx, fits["file"], dirs["m1rc_gmom"], 0.2, held_m1)
    assert X["m1rc_gmom"].tobytes() == t62["m1rc_gmom"].X.tobytes()  # the gate arm is T62's candidate
    assert X[tr.STEADY].tobytes() == X["m1rc_gmom"].tobytes()  # m1 := D / 5: bit for bit
    assert tr.prediction_sha(arms[tr.STEADY]) == tr.prediction_sha(arms["m1rc_gmom"])
    for a in ("m1rc_timing_gmom", "m1rc_timing_perm_gmom"):
        assert X[a].tobytes() != X["m1rc_gmom"].tobytes()
        alone = tr.rc.cut_pair(ctx, fits["file"], dirs[a], 0.2, held_m1)  # the candidate's build on its directions
        assert X[a].tobytes() == alone["m1rc_gmom"].X.tobytes()
    assert tuple(tr.cut_arms(ctx, fits["file"], dirs, 0.2, held_m1)) == tr.ARMS


def test_the_library_arms(tr):
    e95, e85 = _world()
    spec = load_spec()
    fb = {"E8.5": spec.fallbacks("E8.5_RNA"), "E9.5": spec.fallbacks("E9.5_RNA")}
    within = {"surrounding": ["Foregut"]}
    factors = tr.br.gmom_fits(e95, "E9.5_RNA", e85, "E8.5_RNA")["file"]
    shares = {"cm": 0.5, "endocardium": 0.2, "epicardium": 0.15, "blood": 0.05, "surrounding": 0.1}
    held_m1 = _held_m1(list(e95.var_names))
    dirs = _arm_dirs(tr, e95.n_vars)
    arms = tr.library_arms(e95, shares, factors, HELD, held_m1, 6, within, fb, dirs, 0.2, n_cells=300, steady=True)
    c = {a: _dense(counts_of(v)) for a, v in arms.items()}
    t62 = tr.rc.library_pair(e95, shares, factors, HELD, held_m1, 6, within, fb, dirs["m1rc_gmom"], 0.2, n_cells=300)
    assert np.array_equal(c["m1rc_gmom"], _dense(counts_of(t62["m1rc_gmom"])))
    assert np.array_equal(c[tr.STEADY], c["m1rc_gmom"])
    assert not np.array_equal(c["m1rc_timing_gmom"], c["m1rc_gmom"])
    assert not np.array_equal(c["m1rc_timing_gmom"], c["m1rc_timing_perm_gmom"])
    for a in ("m1rc_timing_gmom", "m1rc_timing_perm_gmom"):
        alone = tr.rc.library_pair(e95, shares, factors, HELD, held_m1, 6, within, fb, dirs[a], 0.2, n_cells=300)
        assert np.array_equal(c[a], _dense(counts_of(alone["m1rc_gmom"])))


# ---------------------------------------------------------------------------------------------------------------------
# The reading on hand-made tables

WEIGHTS = {"de_recovery": 0.25, "direction": 0.25, "distribution": 0.30, "covariance": 0.20}


def _arm(score, gain=0.0, share=0.8):
    """An arm of `score` whose gain over the base splits `share` into de_recovery plus direction (points)."""
    g = {"de_recovery": 0.5 + gain * share / 50, "direction": 0.5 + gain * share / 50,
         "distribution": 0.5 + gain * (1 - share) / 60, "covariance": 0.5 + gain * (1 - share) / 40}
    return {"score": float(score), "groups": g}


def _seq(v, s):
    return v[s] if isinstance(v, (list, tuple)) else v


def _late_fields(sha="a" * 64, **change):
    return {"stats_sha256": sha, "stats_sha256_hashed": sha, "yaml_sha256": TR_YAML, "authors_labels_read": False,
            "late_row": "N99", "usable": True, **change}


TR_YAML = "883e74c2c1ef5f01d5244555f8b3dcae24255568c731a6d53ee2788b7b677b49"
TR_STATS = "36d1f7edcbbf3d7a80c68b51d5624a94bb753668079fe3d69aed720c71baa7e9"


def _relabel(**change):
    return {"sha256": TR_STATS, "stats_sha256_hashed": TR_STATS, "yaml_sha256": TR_YAML, "route": "recomputed",
            "record_route": "recomputed", "authors_labels_read": False, **change}


def _task(lib=False, seeds=range(10), commit="c" * 40, clean=True, timing=0.5, perm=0.1, share=0.8, gate=0.0,
          relabel=None, late=None, steady=True, steady_equal=True, sha=True, inverted=False):
    rows = []
    for s in seeds:
        base = 50.0 + s / 100
        t, p = _seq(timing, s), _seq(perm, s)
        arms = {"m1rc_gmom": _arm(base + _seq(gate, s)), "m1rc_timing_gmom": _arm(base + t, t, share),
                "m1rc_timing_perm_gmom": _arm(base + p, p, share)}
        if lib:
            arms["copy_last"] = _arm(40.0)
        row = {"seed": s, "arms": arms, "sha256": {a: f"{a}-{s}" for a in arms} if sha else {}}
        if lib:
            row["inverted_anchors"] = ["x"] if inverted and s == 2 else []
        rows.append(row)
    if rows:
        rows[0]["steady_equal"] = steady_equal
    return {"rows": rows, "stamp": {"commit": commit, "tree_clean": clean}, "fraction": 0.2,
            "relabel": _relabel() if relabel is None else relabel, "late": _late_fields() if late is None else late,
            "timing": {"steady_identity": steady}}


def _loaded(chamber=None, whole=None, libs=None):
    out = {"chamber": _task(**(chamber or {})), "whole": _task(**(whole or {}))}
    for lib in ("1", "2", "3"):
        out[lib] = _task(lib=True, **((libs or {}).get(lib) or {}))
    return out


COMMITTED = {k: {s: 50.0 + s / 100 for s in range(10)} for k in ("chamber", "whole", "1", "2", "3")}
PREGATE = {"idle": False, "stamp": {"commit": "c" * 40, "tree_clean": True},
           "groups": {g: {"idle_ratio": 0.2, "perm_idle_ratio": 0.2, "fallback_share": 0.3, "group_ok": True}
                      for g in GROUPS}}
LATE_OK = {"stats_sha256": "a" * 64, "reasons": []}


def _read(tr, loaded, committed=COMMITTED, pregate=PREGATE, checks=None, late=LATE_OK, source=()):
    return tr.t67_reading(loaded, committed, pregate, checks or {"reasons": []}, late, source)


def test_the_constants(tr):
    assert tr.ROW == "T67" and tr.ARMS == ("m1rc_gmom", "m1rc_timing_gmom", "m1rc_timing_perm_gmom")
    assert tr.COUNTED == ("m1rc_timing_gmom", "m1rc_timing_perm_gmom") and tr.TRIES_COUNT == 10
    assert tr.LIB_ARMS[0] == "copy_last" and tr.GATE_TOL == 0.02 and tr.PRACTICE_FRACTION == 0.2
    assert (tr.CHAMBER_BAR, tr.WHOLE_BAR, tr.LIB_BAR, tr.LIB_FLOOR, tr.PERM_CHAMBER_BAR, tr.PERM_LIB_BAR,
            tr.MIN_SEEDS, tr.SHARE_BAR) == (0.3, 0.0, 0.3, -0.3, 0.2, 0.3, 8, 0.6)
    assert tr.rc.STATS_SHA256 == TR_STATS and tr.rc.YAML_SHA256 == TR_YAML == tr.frl.YAML_SHA256
    assert (tt.PERM_SEED, tt.PERM_D_BINS, tt.PERM_LEVEL_BINS, tt.IDLE_BAR, tt.MIN_LATE_CELLS) == (20261122, 20, 5,
                                                                                                0.05, 50)
    assert tt.SLOPE_STAGES == ("E14.5", "E15.5", "E16.5", "E17.5") and tr.frl.LATE_STAGES == tt.LATE_STAGES
    assert tr.source_guard() == [] and tr.rc.source_guard() == []


def test_the_default_table_is_ahead(tr):
    r = _read(tr, _loaded())
    assert r["outcome"] == "AHEAD" and r["ahead"] and r["i"]["holds"] and r["ii"]["holds"] and r["iii"]["holds"]
    assert r["vs_m1rc"]["chamber"]["mean"] == pytest.approx(0.5)
    assert r["vs_perm"]["chamber"]["mean"] == pytest.approx(0.4)
    assert r["iii"]["share"] == pytest.approx(0.8)
    assert tr.outcome_lines(r)[-1].startswith("What follows: one priced question")


def _p(mean, ahead=10, groups=None):
    return {"mean": mean, "ahead": ahead, "groups": groups or {g: 0.0 for g in WEIGHTS}}


def test_clause_i_boundaries(tr):
    ok = {"chamber": _p(0.3, 8), "whole": _p(0.01), "1": _p(0.3), "2": _p(0.3), "3": _p(0.3, 8)}
    assert tr.clause_i(ok)["holds"]
    for change, part in (({"chamber": _p(0.29, 10)}, "chamber"), ({"chamber": _p(0.5, 7)}, "chamber"),
                         ({"whole": _p(0.0)}, "whole"), ({"3": _p(0.3, 7)}, "library_3_seeds"),
                         ({"1": _p(0.29)}, "library_mean"), ({"1": _p(-0.31), "2": _p(1.0)}, "no_library_below_floor")):
        r = tr.clause_i({**ok, **change})
        assert not r["holds"] and not r["parts"][part], change
    assert tr.clause_i({**ok, "1": _p(-0.3), "2": _p(1.0)})["holds"]  # -0.3 is not below the floor


def test_clause_ii_boundaries(tr):
    ok = {"chamber": _p(0.2, 8), "whole": _p(-5.0), "1": _p(0.3), "2": _p(0.3), "3": _p(0.3, 0)}
    assert tr.clause_ii(ok)["holds"]  # the whole cut and library 3's seeds do not enter (ii)
    assert not tr.clause_ii({**ok, "chamber": _p(0.19, 10)})["holds"]
    assert not tr.clause_ii({**ok, "chamber": _p(0.5, 7)})["holds"]
    assert not tr.clause_ii({**ok, "2": _p(0.2)})["holds"]


def test_clause_iii(tr):
    def tasks(gain, de):
        return {k: _p(gain, groups={"de_recovery": de / 2, "direction": de / 2, "distribution": gain - de,
                                    "covariance": 0.0}) for k in ("chamber", "whole", "1", "2", "3")}

    assert tr.clause_iii(tasks(1.0, 0.6))["holds"] and not tr.clause_iii(tasks(1.0, 0.59))["holds"]
    assert not tr.clause_iii(tasks(-1.0, -0.9))["holds"]  # no positive gain: it cannot hold
    assert not tr.clause_iii(tasks(0.0, 0.0))["holds"]


@pytest.mark.parametrize("change", [
    {"chamber": {"timing": 0.25}},  # (i) chamber mean
    {"chamber": {"timing": [0.5] * 7 + [-0.1] * 3}},  # (i) chamber seeds
    {"whole": {"timing": -0.1, "perm": -0.6}},  # (i) whole cut
    {"libs": {"3": {"timing": [0.5] * 7 + [-0.1] * 3, "perm": -0.5}}},  # (i) library 3's seeds
    {"libs": {"1": {"timing": -0.4, "perm": -1.0}, "2": {"timing": 1.2, "perm": 0.0}}},  # (i) a library below -0.3
    {"chamber": {"perm": 0.35}},  # (ii) chamber against the control
    {"libs": {k: {"perm": 0.3} for k in "123"}},  # (ii) the libraries against the control
    {"chamber": {"share": 0.3}, "whole": {"share": 0.3}, "libs": {k: {"share": 0.3} for k in "123"}},  # (iii)
])
def test_every_failed_clause_reads_not_ahead(tr, change):
    r = _read(tr, _loaded(**change))
    assert r["outcome"] == tr.NOT_AHEAD_WORDS and not r["ahead"]
    assert "closes 'per-gene Hermite timing from Feng's post-window slope' on T1" in r["outcome"]


@pytest.mark.parametrize("case", [
    "missing task", "missing seed", "missing arm", "missing sha", "gate", "commit", "unclean", "late sha",
    "late yaml", "late authors", "authors", "stats", "yaml", "route", "steady", "steady build", "pregate missing",
    "pregate idle", "pregate commit", "checks", "late checks", "source", "fraction"])
def test_every_void_condition(tr, case):
    loaded, pregate, checks, late, source = _loaded(), PREGATE, {"reasons": []}, LATE_OK, ()
    if case == "missing task":
        loaded["2"] = None
    elif case == "missing seed":
        loaded["whole"] = _task(seeds=range(9))
    elif case == "missing arm":
        del loaded["3"]["rows"][4]["arms"]["copy_last"]
    elif case == "missing sha":
        loaded["chamber"] = _task(sha=False)
    elif case == "gate":
        loaded["1"] = _task(lib=True, gate=[0.0] * 9 + [0.021])
    elif case == "commit":
        loaded["whole"] = _task(commit="d" * 40)
    elif case == "unclean":
        loaded["chamber"] = _task(clean=False)
    elif case == "late sha":
        loaded["2"] = _task(lib=True, late=_late_fields(sha="b" * 64))
    elif case == "late yaml":
        loaded["2"] = _task(lib=True, late=_late_fields(yaml_sha256="0" * 64))
    elif case == "late authors":
        loaded["2"] = _task(lib=True, late=_late_fields(authors_labels_read=None))
    elif case == "authors":
        loaded["chamber"] = _task(relabel=_relabel(authors_labels_read=True))
    elif case == "stats":
        loaded["chamber"] = _task(relabel=_relabel(stats_sha256_hashed="0" * 64))
    elif case == "yaml":
        loaded["chamber"] = _task(relabel=_relabel(yaml_sha256="0" * 64))
    elif case == "route":
        loaded["chamber"] = _task(relabel=_relabel(route="none"))
    elif case == "steady":
        loaded["3"] = _task(lib=True, steady=False)
    elif case == "steady build":
        loaded["3"] = _task(lib=True, steady_equal=False)
    elif case == "pregate missing":
        pregate = None
    elif case == "pregate idle":
        pregate = {**PREGATE, "idle": True}
    elif case == "pregate commit":
        pregate = {**PREGATE, "stamp": {"commit": "e" * 40, "tree_clean": True}}
    elif case == "checks":
        checks = {"reasons": ["the settings file's sha256 is x"]}
    elif case == "late checks":
        late = {**LATE_OK, "reasons": ["row T67's line does not name the late statistics' sha256"]}
    elif case == "source":
        source = ["the source imports ['perm_read']"]
    elif case == "fraction":
        loaded["whole"] = {**_task(), "fraction": 0.4}
    r = _read(tr, loaded, pregate=pregate, checks=checks, late=late, source=source)
    assert r["outcome"] == "VOID" and r["reasons"], case
    assert tr.outcome_lines(r)[0].startswith("Outcome: VOID")


def test_a_gate_within_the_tolerance_passes(tr):
    loaded = _loaded()
    loaded["1"] = _task(lib=True, gate=[0.0] * 9 + [0.0199])
    assert _read(tr, loaded)["outcome"] == "AHEAD"
    assert _read(tr, _loaded(libs={"2": {"inverted": True}}))["outcome"] == "AHEAD"  # inverted anchors: reported only


def test_committed_t62_reads_the_gate_arm(tr, tmp_path):
    rows = [{"seed": s, "arms": {"m1rc_gmom": {"score": 51.0 + s}, "none_gmom": {"score": 0.0}}} for s in range(3)]
    (tmp_path / "lib1.json").write_text(json.dumps({"rows": rows}))
    assert tr.committed_t62("lib1.json", tmp_path) == {0: 51.0, 1: 52.0, 2: 53.0}


def test_the_read_prints_no_number_and_writes_its_files_once(tr, tmp_path, monkeypatch, capsys):
    tries = tmp_path / "TRIES.md"
    tries.write_text("| T67 (one row) | T1 | x |\n")
    monkeypatch.setattr(tr, "TRIES", tries)
    out = tmp_path / "results" / "t67"
    tasks = out / "tasks"
    tasks.mkdir(parents=True)
    for k, f in tr.task_files().items():
        (tasks / f).write_text(json.dumps(_task(lib=k in ("1", "2", "3"))))
    (tasks / "pregate.json").write_text(json.dumps(PREGATE))
    monkeypatch.setattr(tr, "committed_t62", lambda f, tasks=None: COMMITTED["chamber"])
    monkeypatch.setattr(tr.rc, "relabel_checks", lambda **k: {"reasons": []})
    monkeypatch.setattr(tr, "late_checks", lambda **k: LATE_OK)
    monkeypatch.setattr(tr.gv, "tree_state", lambda: {"commit": "c" * 40, "clean": True})
    res = tr.read(out)
    printed = capsys.readouterr().out
    assert res["reading"]["outcome"] == "AHEAD" and "Outcome: AHEAD" in printed
    assert not re.search(r"[+-]?\d+\.\d", printed)  # no score, no mean
    assert (out / "t67.json").exists() and "Tries counted: 10" in (out / "t67.txt").read_text()
    with pytest.raises(SystemExit, match="already holds"):
        tr.read(out)


# ---------------------------------------------------------------------------------------------------------------------
# Refusals, the late checks, the pre-gate


def test_the_runs_the_pregate_and_the_read_are_refused_without_the_row(tr, tmp_path, monkeypatch):
    tries = tmp_path / "TRIES.md"
    tries.write_text("| T63 (one row) | T1 | x |\n")
    monkeypatch.setattr(tr, "TRIES", tries)
    monkeypatch.setattr(tr.gv, "tree_state", lambda: pytest.fail("nothing before the row's refusal"))
    monkeypatch.setattr(tr.rc, "relabel_checks", lambda *a, **k: pytest.fail("nothing before the row's refusal"))
    for argv in (["pregate", "--out", "x"], ["run-cut", "--cut", "chamber", "--out", "x"],
                 ["run-lib", "--library", "3", "--out", "x"], ["read", "--out", "x"]):
        with pytest.raises(SystemExit, match="no row"):
            tr.main(argv)


def test_a_run_is_refused_in_order(tr, tmp_path, monkeypatch):
    tries = tmp_path / "TRIES.md"
    tries.write_text("| T67 (one row) | T1 | x |\n")
    monkeypatch.setattr(tr, "TRIES", tries)
    out = tmp_path / "out"
    out.mkdir()
    (out / "lib3.json").write_text("{}")
    with pytest.raises(SystemExit, match="already holds"):
        tr.refuse(out, ("lib3.json",))
    monkeypatch.setattr(tr.gv, "tree_state", lambda: {"commit": "c" * 40, "clean": False})
    with pytest.raises(SystemExit, match="not clean"):
        tr.refuse(out, ("lib1.json",))
    monkeypatch.setattr(tr.gv, "tree_state", lambda: {"commit": "c" * 40, "clean": True})
    monkeypatch.setattr(tr, "source_guard", lambda: ["the source names x"])
    with pytest.raises(SystemExit, match="the source names x"):
        tr.refuse(out, ("lib1.json",))
    monkeypatch.setattr(tr, "source_guard", lambda: [])
    monkeypatch.setattr(tr.rc, "relabel_checks", lambda *a, **k: {"passes": False, "reasons": ["the yaml"]})
    with pytest.raises(SystemExit, match="the yaml"):
        tr.refuse(out, ("lib1.json",))
    monkeypatch.setattr(tr.rc, "relabel_checks", lambda *a, **k: {"passes": True, "reasons": []})
    monkeypatch.setattr(tr, "late_checks", lambda: {"passes": False, "reasons": ["no late labelling manifest"]})
    with pytest.raises(SystemExit, match="no late labelling manifest"):
        tr.refuse(out, ("lib1.json",))
    monkeypatch.setattr(tr, "late_checks", lambda: {"passes": True, "reasons": []})
    monkeypatch.setattr(tr.frl, "x002_names_this_use", lambda: False)
    with pytest.raises(SystemExit, match="X002"):
        tr.refuse(out, ("lib1.json",))
    monkeypatch.setattr(tr.frl, "x002_names_this_use", lambda: True)
    with pytest.raises(SystemExit, match="pre-gate has not run"):
        tr.refuse(out, ("lib1.json",))
    assert tr.refuse(out, ("pregate.json",), need_pregate=False)[0]["clean"] is True
    pg = out / "pregate.json"
    for doc, why in (({"idle": True}, "idle"), ({"idle": False, "steady_identity": False}, "steady"),
                     ({"idle": False, "steady_identity": True, "stamp": {"commit": "d" * 40, "tree_clean": True}},
                      "this commit")):
        pg.write_text(json.dumps(doc))
        with pytest.raises(SystemExit, match=why):
            tr.refuse(out, ("lib1.json",))
    pg.write_text(json.dumps({"idle": False, "steady_identity": True,
                              "stamp": {"commit": "c" * 40, "tree_clean": True}}))
    assert tr.refuse(out, ("lib1.json",))[0]["commit"] == "c" * 40


def _late_files(tmp_path, frl_mod, **change):
    d = tmp_path / "late"
    d.mkdir(exist_ok=True)
    stats = d / frl_mod.STATS
    stats.write_bytes(b"late statistics")
    sha = frl_mod.fio.sha256_file(stats)
    man = {"source": "feng_relabel_late", "config": {"sha256": TR_YAML}, "authors_labels_read": False,
           "stages": ["E15.5", "E16.5", "E17.5"], "usable": True, "gate": {"row": "N99", "clean": True},
           "outputs": {frl_mod.STATS: {"sha256": sha}}, **change}
    path = d / "MANIFEST.json"
    path.write_text(json.dumps(man))
    return path, stats, sha


def test_the_late_checks(tr, tmp_path):
    path, stats, sha = _late_files(tmp_path, tr.frl)
    raw = path.read_bytes()

    def run(**kw):
        args = {"manifest": path, "stats_path": stats, "line": f"| T67 (x) | T1 | sha256 {sha} |",
                "record_at_head": raw, **kw}
        return tr.late_checks(**args)

    ok = run()
    assert ok["passes"] and ok["stats_sha256"] == sha == ok["stats_sha256_hashed"] and ok["late_row"] == "N99"
    assert not run(line="| T67 (x) | T1 | sha256 <stats sha> |")["passes"]  # the row must name the sha256
    assert not run(record_at_head=b"{}")["passes"]  # the RECORD at HEAD must be the manifest
    assert not run(stats_path=tmp_path / "absent.npz")["passes"]
    assert run(stats_path=tmp_path / "absent.npz", need_stats=False)["passes"]
    for k, change in enumerate(({"source": "feng_relabel"}, {"config": {"sha256": "0" * 64}},
                                {"authors_labels_read": True}, {"stages": ["E15.5", "E16.5"]}, {"usable": False},
                                {"stopped": True}, {"gate": {"row": "T67", "clean": True}},
                                {"gate": {"row": "N99", "clean": False}})):
        sub = tmp_path / f"c{k}"
        sub.mkdir()
        p2, s2, sha2 = _late_files(sub, tr.frl, **change)
        r = tr.late_checks(manifest=p2, stats_path=s2, line=f"| T67 | {sha2} |", record_at_head=p2.read_bytes())
        assert not r["passes"], change
    assert not tr.late_checks(manifest=tmp_path / "none.json", stats_path=stats, line=sha,
                              record_at_head=None)["passes"]


def test_the_pregate_prunes_an_idle_arm_and_prints_one_step_line(tr, tmp_path, monkeypatch, capsys):
    started = {"commit": "c" * 40, "clean": True}
    monkeypatch.setattr(tr, "refuse", lambda out, names, need_pregate=True: (started, {}, LATE_OK | {
        "stats_sha256_hashed": "a" * 64, "yaml_sha256": TR_YAML, "authors_labels_read": False, "late_row": "N99",
        "usable": True}))
    monkeypatch.setattr(tr, "get_board", lambda key: type("B", (), {"genes": ["G0", "G1"]}))
    monkeypatch.setattr(tr.rc, "_relabel_fields", lambda ext, checks: {})
    monkeypatch.setattr(tr.gv, "stamp", lambda s, inputs: {"commit": "c" * 40, "tree_clean": True, "utc": "now"})
    groups = {"CM": {"idle_ratio": 0.01, "perm_idle_ratio": 0.02, "fallback_share": 0.9, "group_ok": True,
                     "genes": 2, "cells": {"E14.5": 99}, "source_group": "CM"}}
    for idle, code in ((True, 4), (False, 0)):
        monkeypatch.setattr(tr, "load_kit", lambda genes, idle=idle: ({}, {}, {"groups": groups, "idle": idle,
                                                                             "steady_identity": True,
                                                                             "direction_sha256": {}}, {}))
        got = tr.main(["pregate", "--out", f"p{code}"])
        printed = capsys.readouterr().out
        assert got == code and printed.count("\n") == 1 and not re.search(r"\d", printed)
        saved = json.loads((tr.RESULTS / f"p{code}" / "tasks" / "pregate.json").read_text())
        assert saved["idle"] is idle and saved["bar"] == 0.05
    assert tr.pregate_refusal(tr.RESULTS / "p4" / "tasks" / "pregate.json", "c" * 40).startswith("the pre-gate found")
    assert tr.pregate_refusal(tr.RESULTS / "p0" / "tasks" / "pregate.json", "c" * 40) is None


def test_the_source_guard_refuses_other_routes(tr, tmp_path):
    for text in ("import perm_read as pr\n", 'x = fd.feng_directions(list(g))\n', 'g = "auth' + 'ors"\n',
                 "x = __imp" + "ort__('eb_read')\n"):
        bad = tmp_path / "bad.py"
        bad.write_text(text)
        assert tr.source_guard(bad), text


# ---------------------------------------------------------------------------------------------------------------------
# The late labelling (target-free; its refusals only)


def test_the_late_labelling_changes_only_the_stages(frl):
    cfg = frl.fr.load_config(frl.fr.CONFIG)
    late = frl.late_config(cfg)
    assert late["stages"] == ["E15.5", "E16.5", "E17.5"] and cfg["stages"] == ["E9.5", "E14.5"]
    assert {k: v for k, v in late.items() if k != "stages"} == {k: v for k, v in cfg.items() if k != "stages"}
    assert frl.fio.sha256_file(frl.fr.CONFIG) == frl.YAML_SHA256  # the frozen settings, unchanged
    text = LATE.read_text()
    assert "cell_type_spec" not in text.split('"""', 2)[2]  # the code never names the authors' column
    assert "fr.relabel(" in text and "feng2022_geo.group_stats(" in text


def test_the_x002_check_and_the_row_line(frl, tmp_path):
    reg = tmp_path / "EXTERNAL_DATA.md"
    reg.write_text("| X001 | feng_relabel_late |\n| X002 | Feng 2022 ... |\n")
    assert not frl.x002_names_this_use(reg)
    reg.write_text("| X002 | Feng 2022 ...; amended: data/derived/feng_relabel_late (E15.5-E17.5) |\n")
    assert frl.x002_names_this_use(reg)
    tries = tmp_path / "TRIES.md"
    tries.write_text("| N70 (target-free) | T1 | experiments/external/feng_relabel_late.py |\n| N7 | x |\n")
    assert "feng_relabel_late.py" in frl.row_line("N70", tries) and frl.row_line("N7", tries) == "| N7 | x |"
    assert frl.row_line("N71", tries) is None


def test_the_late_preflight_refuses_in_order(frl, tmp_path, monkeypatch):
    dst = tmp_path / "data" / "derived" / "feng_relabel_late"
    with pytest.raises(frl.Refused, match="not an N row"):
        frl.preflight("T67", dst)
    monkeypatch.setattr(frl, "YAML_SHA256", "0" * 64)
    with pytest.raises(frl.Refused, match="frozen"):
        frl.preflight("N70", dst)
    monkeypatch.setattr(frl, "YAML_SHA256", TR_YAML)
    head = {frl.fr.CONFIG: frl.fr.CONFIG.read_bytes(), frl.fr.RECORD: b"{}", frl.RECORD: None}
    monkeypatch.setattr(frl.fr, "committed_at_head", lambda p: head.get(p))
    with pytest.raises(frl.Refused, match="route recomputed"):
        frl.preflight("N70", dst)
    head[frl.fr.RECORD] = json.dumps({"config": {"sha256": TR_YAML}, "route": "recomputed"}).encode()
    monkeypatch.setattr(frl, "x002_names_this_use", lambda: False)
    with pytest.raises(frl.Refused, match="X002"):
        frl.preflight("N70", dst)
    monkeypatch.setattr(frl, "x002_names_this_use", lambda: True)
    monkeypatch.setattr(frl, "row_line", lambda row: "| N70 | T1 | another script |")
    with pytest.raises(frl.Refused, match="naming feng_relabel_late.py"):
        frl.preflight("N70", dst)
    monkeypatch.setattr(frl, "row_line", lambda row: "| N70 | T1 | feng_relabel_late.py |")
    from vec import provenance

    monkeypatch.setattr(provenance, "git_state", lambda repo: ("c" * 40, ["x.py"]))
    with pytest.raises(frl.Refused, match="not clean"):
        frl.preflight("N70", dst)
    monkeypatch.setattr(provenance, "git_state", lambda repo: ("c" * 40, []))
    monkeypatch.setattr(frl, "RECORD", tmp_path / "record" / "MANIFEST.json")
    head[frl.RECORD] = None
    dst.mkdir(parents=True)
    with pytest.raises(frl.Refused, match="already recorded or was started"):
        frl.preflight("N70", dst)
    gate = frl.preflight("N70", tmp_path / "fresh")
    assert gate == {"row": "N70", "commit": "c" * 40, "clean": True, "config_sha256": TR_YAML,
                    "first_record_route": "recomputed"}
    with pytest.raises(frl.Refused, match="no clean gate"):
        frl.run({}, tmp_path, tmp_path / "x", {})


# ---------------------------------------------------------------------------------------------------------------------
# The launcher


def test_the_launcher_deletes_nothing_and_opens_no_quote_in_a_usage_message():
    text = LAUNCHER.read_text()
    code = [line.split("#")[0] for line in text.splitlines()]
    assert not any(re.search(r"(^|[\s;&|(])(rm|rmdir|unlink|shred|truncate|mv)\s", line) for line in code)
    for msg in re.findall(r"\$\{\d+:\?([^}]*)\}", text):
        assert "'" not in msg and '"' not in msg, msg
    assert "OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 NUMPY_MADVISE_HUGEPAGE=0" in text
    assert "python -I -B" in text and "--basetemp=" in text and "^| T67 " in text
    listed = re.findall(r"tests/test_[a-z0-9_]+\.py", text)
    assert "tests/test_t1_timing.py" in listed and all((ROOT / f).exists() for f in listed)
    n54 = (ROOT / "experiments" / "external" / "n54_step.sh").read_text()
    assert "tests/test_t1_timing.py" in n54.split("  tests)")[1].split(";;")[0]


def _repo(tmp_path, row=True):
    """A throwaway git repository holding the launcher and a ledger (with or without row T67 and row N70)."""
    repo = tmp_path / "repo"
    (repo / "experiments" / "t1_m1").mkdir(parents=True)
    (repo / "experiments" / "harness").mkdir(parents=True)
    shutil.copy(LAUNCHER, repo / "experiments" / "t1_m1" / "timing.sh")
    (repo / "experiments" / "harness" / "TRIES.md").write_text(
        "| T62 | T1 | x |\n" + ("| T67 (one row) | T1 | x |\n| N70 (target-free) | T1 | x |\n" if row else ""))
    (repo / ".gitignore").write_text("data/\n")
    env = _env(tmp_path)
    for args in (["init", "-q"], ["add", "-A"], ["-c", "user.name=t", "-c", "user.email=t@t", "commit", "-qm", "x"]):
        subprocess.run(["git", *args], cwd=repo, env=env, check=True, capture_output=True)
    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=repo, env=env, capture_output=True, text=True)
    return repo, head.stdout.strip()


def _env(tmp_path):
    env = {k: v for k, v in os.environ.items() if k not in ("GIT_DIR", "GIT_WORK_TREE", "GIT_INDEX_FILE")}
    env["HOME"] = str(tmp_path / "home")
    env["GIT_CEILING_DIRECTORIES"] = str(tmp_path)
    return env


def _launch(tmp_path, repo, *argv):
    out = subprocess.run(["bash", "experiments/t1_m1/timing.sh", *argv], cwd=repo, env=_env(tmp_path),
                         capture_output=True, text=True)
    return out.returncode


def _exit(tmp_path, name):
    return (tmp_path / "home" / "v" / "scratch" / "t67" / f"{name}.exit").read_text().strip()


def _mark(tmp_path, name, code, head):
    s = tmp_path / "home" / "v" / "scratch" / "t67"
    s.mkdir(parents=True, exist_ok=True)
    (s / f"{name}.exit").write_text(f"{code}\n")
    (s / f"{name}.head").write_text(f"{head}\n")


def test_the_launchers_refusal_paths(tmp_path):
    repo, head = _repo(tmp_path)
    assert _launch(tmp_path, repo) == 2 and _exit(tmp_path, "usage") == "2"
    assert _launch(tmp_path, repo, "tests", head[:12]) == 2  # a prefix is never enough
    assert _launch(tmp_path, repo, "late", head, "T67") == 2  # the late chain takes an N row
    assert _launch(tmp_path, repo, "build", head) == 97 and _exit(tmp_path, "unknown") == "97"
    assert _launch(tmp_path, repo, "cut", "half", head) == 97
    assert _launch(tmp_path, repo, "pregate", "0" * 40) == 95 and _exit(tmp_path, "pregate") == "95"
    # the step before did not exit 0 at this HEAD
    assert _launch(tmp_path, repo, "pregate", head) == 92 and _exit(tmp_path, "pregate") == "92"
    assert _launch(tmp_path, repo, "late", head, "N70") == 92
    for step in (["cut", "chamber"], ["cut", "whole"], ["lib", "1"], ["lib", "2"], ["lib", "3"], ["read"]):
        assert _launch(tmp_path, repo, *step, head) == 92, step
    _mark(tmp_path, "tests", 0, "f" * 40)  # exit 0, but at another commit
    assert _launch(tmp_path, repo, "pregate", head) == 92
    _mark(tmp_path, "tests", 1, head)
    assert _launch(tmp_path, repo, "pregate", head) == 92
    _mark(tmp_path, "tests", 0, head)
    assert _launch(tmp_path, repo, "pregate", head) == 93  # the statistics are not in place
    _mark(tmp_path, "late_tests", 0, head)
    assert _launch(tmp_path, repo, "late", head, "N70") == 93
    data = repo / "data" / "derived"
    (data / "feng_relabel_c").mkdir(parents=True)
    (data / "feng_relabel_c" / "cd1_stats_relabel.npz").write_text("x")
    (data / "feng_relabel_late").mkdir()
    assert _launch(tmp_path, repo, "late", head, "N70") == 93  # its output folder exists
    (data / "feng_relabel_late" / "cd1_stats_relabel_late.npz").write_text("x")
    tasks = repo / "experiments" / "t1_m1" / "results" / "t67_timing" / "tasks"
    tasks.mkdir(parents=True)
    (tasks / "pregate.json").write_text("{}")  # untracked in the row's task folder: the tree stays clean
    assert _launch(tmp_path, repo, "pregate", head) == 91 and _exit(tmp_path, "pregate") == "91"
    _mark(tmp_path, "pregate", 4, head)  # idle: no run may start
    assert _launch(tmp_path, repo, "cut", "chamber", head) == 92
    _mark(tmp_path, "pregate", 0, head)
    (tasks / "lib2.json").write_text("{}")
    assert _launch(tmp_path, repo, "lib", "2", head) == 91 and _exit(tmp_path, "lib_2") == "91"
    for name in ("cut_chamber", "cut_whole", "lib_1", "lib_2", "lib_3"):
        _mark(tmp_path, name, 0, head)
    assert _launch(tmp_path, repo, "read", head) == 90  # task files missing
    for f in ("T1_E9.5_E10.5_chamber.json", "T1_E9.5_E10.5.json", "lib1.json", "lib3.json"):
        (tasks / f).write_text("{}")
    (tasks.parent / "t67.txt").write_text("x")
    assert _launch(tmp_path, repo, "read", head) == 91
    (repo / "stray.txt").write_text("x")  # an untracked file outside the task folder
    assert _launch(tmp_path, repo, "read", head) == 98 and _exit(tmp_path, "read") == "98"


def test_the_launcher_refuses_without_the_rows(tmp_path):
    repo, head = _repo(tmp_path, row=False)
    assert _launch(tmp_path, repo, "tests", head) == 96 and _exit(tmp_path, "tests") == "96"
    assert _launch(tmp_path, repo, "late-tests", head, "N70") == 96 and _exit(tmp_path, "late_tests") == "96"
