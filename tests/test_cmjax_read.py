"""Try T71 (experiments/harness/TRIES.md): experiments/t1_m1/cmjax_read.py and its launcher on synthetic stages,
hand-made task files and temporary git repositories only. No released file, reducer output, statistics file, board file
or score is read; the committed T62 and T63 task folders are checked for their file names only. It reuses row T64's and
T68's fixtures."""

from __future__ import annotations

import importlib.util
import inspect
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
from test_g_mom_verify import _small_inputs
from test_t1_scope import GENES, HELD
from test_t1_state_transport import _gmom_ctx

from vec import harness
from vec.core import counts_of
from vec.t1_scope import classify, load_spec, reweight

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "experiments" / "t1_m1" / "cmjax_read.py"
LAUNCHER = ROOT / "experiments" / "t1_m1" / "cmjax.sh"
SHARES = {"cm": 0.6, "endocardium": 0.25, "blood": 0.15}
DIRECTED = ("CM", "Endocardium", "Endothelium")
HELD_M1 = sorted({*HELD, GENES[3], GENES[7]})
WITHIN = {"surrounding": ["Foregut", "Epithelium"]}


@pytest.fixture
def cr(monkeypatch, tmp_path):
    spec = importlib.util.spec_from_file_location("cmjax_read_module", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    monkeypatch.setattr(mod, "RESULTS", tmp_path / "results")
    return mod


def _dense(x):
    return np.asarray(x.todense()) if hasattr(x, "todense") else np.asarray(x)


def _classified():
    ref, start, _ = _small_inputs()
    classify(ref, source="E9.5_RNA")
    classify(start, source="E8.5_RNA")
    return ref, start


def _directions(n_genes, groups=DIRECTED, size=0.4):
    rng = np.random.default_rng(11)
    return {g: rng.normal(0.0, size, n_genes) for g in groups}


def _fb():
    spec = load_spec()
    return {"E8.5": spec.fallbacks("E8.5_RNA"), "E9.5": spec.fallbacks("E9.5_RNA")}


def _cm_entries(cr, dirs, d_jax, genes):
    free = ~np.isin(np.asarray(genes, dtype=str), HELD_M1)
    level = np.random.default_rng(3).gamma(1.0, 50.0, len(genes))
    d_arm, _ = cr.blend(dirs["CM"], d_jax, free)
    plc0, _, _ = cr.placebo(dirs["CM"], d_jax, level, free)
    d_plc, _ = cr.blend(dirs["CM"], plc0, free)
    return d_arm, d_plc


def _lib(cr, e95, e85, d_jax=None, seed=4):
    factors = cr.br.gmom_fits(e95, "E9.5_RNA", e85, "E8.5_RNA")["file"]
    dirs = _directions(e95.n_vars)
    d_jax = dirs["CM"].copy() if d_jax is None else d_jax
    d_arm, d_plc = _cm_entries(cr, dirs, d_jax, list(e95.var_names))
    arms, rec = cr.library_arms(e95, SHARES, factors, HELD, HELD_M1, seed, WITHIN, _fb(), dirs, 0.2, d_arm, d_plc,
                                n_cells=300)
    base0, _ = reweight(e95, SHARES, 300, seed, within=WITHIN, fallback=_fb()["E9.5"])
    return arms, rec, base0, factors, dirs


def _cut(cr, d_jax=None, seed=1000002):
    ref, start, _ = _small_inputs()
    ctx = _gmom_ctx(ref, start, seed=seed)
    fits = cr.br.gmom_fits(ref, "E9.5_RNA", start, "E8.5_RNA")
    dirs = _directions(ref.n_vars)
    d_jax = dirs["CM"].copy() if d_jax is None else d_jax
    d_arm, d_plc = _cm_entries(cr, dirs, d_jax, list(ref.var_names))
    arms, rec = cr.cut_arms(ctx, fits["file"], dirs, 0.2, HELD_M1, d_arm, d_plc)
    return arms, rec, ctx, fits, dirs


def _cm_mask(cr, labels, dirs):
    groups, _, _ = cr.rc.direction_kit(np.unique(labels), dirs)
    return np.array([groups.get(lab) == "CM" for lab in labels], dtype=bool)


def _cut_draw_labels(ctx):
    ref = harness._shallow(ctx.inputs[ctx.task.reference_role])
    src = ctx.task.role_stems[ctx.task.reference_role]
    classify(ref, source=src)
    shares = harness.t1_chamber_shares(ref, src)
    pop, _ = reweight(ref, shares, harness.count_for(ctx.board, ref.n_obs), ctx.seed,
                      fallback=load_spec().fallbacks(src))
    return pop.obs["celltype"].astype(str).to_numpy()


# ---------------------------------------------------------------------------------------------------------------------
# The arms: equal to the base when D_jax is D_feng; only CM cells change; the bases are T62's and T63's


def test_with_d_jax_equal_to_d_feng_every_arm_is_its_base_on_both_routes(cr):
    e95, e85 = _classified()
    arms, rec, _, _, _ = _lib(cr, e95, e85)
    shas = {k: cr.prediction_sha(v.X) for k, v in arms.items()}
    assert shas[cr.ARM1] == shas[cr.PLC1] == shas[cr.BASE1] != shas[cr.BASE2] == shas[cr.ARM2] == shas[cr.PLC2]
    assert all(rec[n]["cm_cells_changed"] == 0 and rec[n]["cm_cells"] > 0 for n in cr.COUNTED)
    arms, _, _, _, _ = _cut(cr)
    shas = {k: cr.prediction_sha(v.X) for k, v in arms.items()}
    assert shas[cr.ARM1] == shas[cr.PLC1] == shas[cr.BASE1] and shas[cr.ARM2] == shas[cr.PLC2] == shas[cr.BASE2]


def test_the_blend_of_d_feng_with_itself_is_d_feng_bit_for_bit(cr):
    d = np.random.default_rng(1).normal(size=200)
    d[:5] = np.nan
    free = np.ones(200, dtype=bool)
    free[10:20] = False
    out, k = cr.blend(d, d.copy(), free)
    assert k == 1.0 and np.array_equal(out, d, equal_nan=True)
    plc, a, _ = cr.placebo(d, d.copy(), np.arange(200.0), free)
    assert a == pytest.approx(1.0) and np.allclose(plc, d, equal_nan=True, atol=1e-12)


def test_on_the_library_route_only_cm_cells_change(cr):
    e95, e85 = _classified()
    arms, rec, base0, _, dirs = _lib(cr, e95, e85, d_jax=np.random.default_rng(5).normal(0, 1.0, e95.n_vars))
    cm = _cm_mask(cr, base0.obs["celltype"].astype(str).to_numpy(), dirs)
    assert cm.any() and not cm.all()
    for base, names in ((cr.BASE1, (cr.ARM1, cr.PLC1)), (cr.BASE2, (cr.ARM2, cr.PLC2))):
        a = _dense(counts_of(arms[base]))
        for name in names:
            differs = (a != _dense(counts_of(arms[name]))).any(axis=1)
            assert differs.any() and not differs[~cm].any(), name
            assert rec[name]["non_cm_equal"] is True and rec[name]["cm_cells_changed"] == int(differs.sum())


def test_on_the_cut_route_only_cm_cells_change(cr):
    arms, rec, ctx, _, dirs = _cut(cr, d_jax=np.random.default_rng(6).normal(0, 1.0, len(GENES)))
    cm = _cm_mask(cr, _cut_draw_labels(ctx), dirs)
    assert cm.any()
    for base, names in ((cr.BASE1, (cr.ARM1, cr.PLC1)), (cr.BASE2, (cr.ARM2, cr.PLC2))):
        for name in names:
            differs = (_dense(arms[base].X) != _dense(arms[name].X)).any(axis=1)
            assert differs.any() and not differs[~cm].any(), name
            assert rec[name]["non_cm_equal"] is True


def test_a_cell_outside_the_cm_group_that_differs_stops_the_run(cr):
    labels = np.array(["V-CM", "Blood", "V-CM", "Endocardium"])
    groups = {"V-CM": "CM", "Endocardium": "Endocardium"}
    X = np.arange(12, dtype=np.float32).reshape(4, 3)
    Y = X.copy()
    Y[0, 1] += 1.0
    assert cr.cm_record(labels, groups, X, Y, cr.ARM1)["cm_cells_changed"] == 1
    Y[1, 0] += 1.0
    with pytest.raises(SystemExit, match="outside the CM group"):
        cr.cm_record(labels, groups, X, Y, cr.ARM1)


def test_the_bases_are_t62s_candidate_and_t63s_doubled_step(cr):
    arms, _, ctx, fits, dirs = _cut(cr, d_jax=np.random.default_rng(7).normal(size=len(GENES)), seed=1000003)
    x1 = cr.rc.cut_pair(ctx, fits["file"], dirs, 0.2, HELD_M1)["m1rc_gmom"]
    x2 = cr.rc.cut_pair(ctx, fits["file"], dirs, 0.4, HELD_M1)["m1rc_gmom"]  # T63's X2 is this call at 2/5
    assert arms[cr.BASE1].X.tobytes() == x1.X.tobytes() and arms[cr.BASE2].X.tobytes() == x2.X.tobytes()
    e95, e85 = _classified()
    arms, _, _, factors, dirs = _lib(cr, e95, e85, d_jax=np.random.default_rng(8).normal(size=e95.n_vars), seed=5)
    for base, f in ((cr.BASE1, 0.2), (cr.BASE2, 0.4)):
        t62 = cr.rc.library_pair(e95, SHARES, factors, HELD, HELD_M1, 5, WITHIN, _fb(), dirs, f, n_cells=300)
        assert _dense(arms[base].X).tobytes() == _dense(t62["m1rc_gmom"].X).tobytes()
    for fn in (cr.cut_arms, cr.library_arms):  # no truth reaches an arm
        assert not {"truth", "target", "e105"} & set(inspect.signature(fn).parameters)


def test_the_arm_dict_changes_the_cm_entry_alone(cr):
    dirs = _directions(30)
    d = np.ones(30)
    got = cr.arm_directions(dirs, d)
    assert set(got) == set(dirs) and np.array_equal(got["CM"], d)
    assert all(got[g] is dirs[g] for g in dirs if g != "CM")
    with pytest.raises(ValueError):
        cr.arm_directions({"Endocardium": d}, d)


# ---------------------------------------------------------------------------------------------------------------------
# The blend and the placebo


def test_the_blend_is_the_mean_rescaled_to_fengs_rms_on_the_free_finite_genes(cr):
    rng = np.random.default_rng(2)
    n = 300
    df, dj = rng.normal(size=n), rng.normal(0, 3.0, n)
    df[:4] = np.nan
    dj[4:9] = np.nan
    free = np.ones(n, dtype=bool)
    free[20:30] = False
    out, k = cr.blend(df, dj, free)
    both = np.isfinite(df) & np.isfinite(dj)
    raw = np.where(both, 0.5 * df + 0.5 * dj, df)
    assert np.allclose(out, raw * k, equal_nan=True) and np.array_equal(np.isnan(out), np.isnan(df))
    m = free & np.isfinite(out)
    assert np.sqrt(np.mean(out[m] ** 2)) == pytest.approx(np.sqrt(np.mean(df[m] ** 2)))
    assert np.allclose(out[4:9], df[4:9] * k)  # no JAX value: Feng's, then the common rescale


def test_the_placebo_keeps_the_shared_part_and_the_residual_values_within_bins(cr):
    rng = np.random.default_rng(4)
    n = 500
    df = rng.normal(size=n)
    dj = 0.7 * df + rng.normal(0, 0.5, n)
    dj[:6] = np.nan
    level = rng.gamma(1.0, 50.0, n)
    free = np.ones(n, dtype=bool)
    free[50:60] = False
    plc, a, perm = cr.placebo(df, dj, level, free)
    g = free & np.isfinite(dj)
    assert a == pytest.approx(float(np.sum(dj[g] * df[g]) / np.sum(df[g] ** 2)))
    r = dj - a * df
    assert np.array_equal(np.isnan(plc), np.isnan(dj)) and np.array_equal(plc[~g], dj[~g], equal_nan=True)
    assert np.allclose(plc[g], a * df[g] + perm[g])
    idx = np.flatnonzero(g)
    order = idx[np.argsort(level[idx], kind="stable")]
    for b in np.array_split(order, 10):
        assert np.array_equal(np.sort(perm[b]), np.sort(r[b]))
    assert not np.allclose(perm[g], r[g])
    again, _, _ = cr.placebo(df, dj, level, free)
    assert np.array_equal(again, plc, equal_nan=True) and cr.PLC_SEED == 20261123 and cr.PLC_BINS == 10
    other, _, _ = cr.placebo(df, dj, level, free, seed=1)
    assert not np.array_equal(other, plc, equal_nan=True)


def test_the_task_inputs_record_the_checks(cr):
    rng = np.random.default_rng(9)
    genes = [f"g{i}" for i in range(80)]
    arrays = {"symbol": np.array(genes[::-1] + ["dup", "dup"]), "direction": rng.normal(size=82),
              "m_early": rng.gamma(1.0, 10.0, 82)}
    dirs = {"CM": rng.normal(size=80), "Endocardium": rng.normal(size=80)}
    got = cr.jax_cm_inputs(dirs, arrays, genes, ["g1", "g2"])
    f = got["fields"]
    assert f["residual_same_values"] and f["plc_nan_where_jax_nan"] and f["g_genes"] == 78
    assert np.array_equal(got["d_jax"], arrays["direction"][:80][::-1])
    assert np.isnan(cr.on_t1(arrays, ["dup"])[0])  # a symbol the table holds twice gets no value


# ---------------------------------------------------------------------------------------------------------------------
# The pre-gate


def _manifest(cr, entering=50, assigned=60, high=3, **change):
    per = {}
    for e in cr.jc.ENDS:
        for n in cr.jc.ROW_EMBRYOS["practice"][e]:
            per[str(n)] = {"cm_assigned": assigned, "cm_entering": entering, "cm_skm_high": high,
                           "day": "E8.5" if e == "early" else "E16.0"}
    man = {"config": {"sha256": cr.jc.YAML_SHA256}, "gate": {"row": "T71", "clean": True, "mode": "practice"},
           "mode": "practice", "gene_universe": "gene_table", "columns_read": ["cell_id", "embryo_id", "day"],
           "authors_labels_read": False, "stats_npz_read": False,
           "embryos_read": sorted(n for e in cr.jc.ENDS for n in cr.jc.ROW_EMBRYOS["practice"][e]),
           "embryo_lists": {e: list(cr.jc.ROW_EMBRYOS["practice"][e]) for e in cr.jc.ENDS},
           "embryos_used": {e: list(cr.jc.ROW_EMBRYOS["practice"][e]) for e in cr.jc.ENDS}, "has_direction": True,
           "per_embryo": per, "delta_days": 7.5}
    man.update(change)
    return man


def _vectors(n=400, rho_noise=1.0):
    rng = np.random.default_rng(12)
    df = rng.normal(size=n)
    dj = df + rng.normal(0, rho_noise, n)
    cp = {"E9.5": np.full(n, 2.0), "E14.5": np.full(n, 0.5)}
    held = np.zeros(n, dtype=bool)
    held[:5] = True
    return dj, df, cp, held


def test_the_pregate_passes_and_reports(cr):
    dj, df, cp, held = _vectors()
    r = cr.pregate_reading(_manifest(cr), dj, df, cp, held)
    assert r["passes"] and r["failed"] == []
    assert r["ends"]["early"]["entering_total"] == 600 and r["ends"]["late"]["entering_total"] == 250
    assert r["late_skm_share"] == pytest.approx(0.05) and 0.2 <= r["spearman"] <= 0.95


@pytest.mark.parametrize("case,part", [("nuclei", "(a)"), ("purity", "(b)"), ("low", "(c)"), ("high", "(c)"),
                                       ("held", "(d)"), ("record", "record"), ("no_direction", "record"),
                                       ("expressed", "(c)")])
def test_the_pregate_fails(cr, case, part):
    dj, df, cp, held = _vectors()
    man = _manifest(cr)
    if case == "nuclei":
        man = _manifest(cr, entering=39)  # 5 x 39 = 195 at the late end
    elif case == "purity":
        man = _manifest(cr, high=7)  # 7 / 60 > 0.10
    elif case == "low":
        dj = np.random.default_rng(1).normal(size=len(df))
    elif case == "high":
        dj = df.copy()
    elif case == "held":
        dj = np.where(held, 100.0, dj)
    elif case == "record":
        man["gate"] = {"row": "T68", "clean": True, "mode": "practice"}
    elif case == "no_direction":
        man["has_direction"] = False
    elif case == "expressed":
        cp = {k: np.zeros(len(df)) for k in cp}  # no gene at the Feng level: no Spearman
    r = cr.pregate_reading(man, dj, df, cp, held)
    assert not r["passes"] and any(part in f for f in r["failed"]), (case, r["failed"])


def test_the_pregate_boundaries_that_pass(cr):
    dj, df, cp, held = _vectors()
    man = _manifest(cr, entering=40, assigned=60, high=6)  # 200 entering, share exactly 0.10
    r = cr.pregate_reading(man, dj, df, cp, held)
    assert r["passes"], r["failed"]
    used = {"early": list(range(1, 13)), "late": [59, 60, 61, 62]}  # only contributing embryos' nuclei count
    r = cr.pregate_reading(_manifest(cr, entering=50, embryos_used=used), dj, df, cp, held)
    assert r["ends"]["late"]["entering_total"] == 200 and r["passes"]


def test_the_pregate_prints_only_words_and_exits_3_on_a_failure(cr, tmp_path, monkeypatch, capsys):
    _ready(cr, tmp_path, monkeypatch)
    man = _manifest(cr, entering=10)
    genes = [f"g{i}" for i in range(40)]
    arrays = {"symbol": np.array(genes), "direction": np.linspace(-1, 1, 40), "m_early": np.ones(40),
              "m_late": np.ones(40)}
    monkeypatch.setattr(cr, "load_reducer", lambda *a: (man, arrays, "s" * 64))
    monkeypatch.setattr(cr, "get_board", lambda b: type("B", (), {"genes": genes})())
    monkeypatch.setattr(cr.rc, "recomputed_inputs", lambda g: ({"CM": np.linspace(-1, 1, len(g))}, {}))
    monkeypatch.setattr(cr.rc, "_held_m1", lambda p: (["g0"], {}))
    monkeypatch.setattr(cr, "feng_cm_means", lambda g: {"E9.5": np.ones(len(g)) * 2, "E14.5": np.ones(len(g))})
    monkeypatch.setattr(cr.gv, "stamp", lambda started, fs: {"commit": "c" * 40, "tree_clean": True, "utc": "x"})

    class A:
        out, commit = "pg", "c" * 40
    assert cr.cmd_pregate(A) == cr.PREGATE_FAIL == 3
    printed = capsys.readouterr().out
    assert printed.startswith("Pre-gate: FAILS") and "no try is counted" in printed
    assert not any(ch.isdigit() for ch in printed.replace("T71", ""))
    rec = json.loads((cr.RESULTS / "pg" / "tasks" / "pregate.json").read_text())
    assert rec["reading"]["passes"] is False and rec["manifest_sha256"] == "s" * 64
    with pytest.raises(SystemExit, match="already holds"):
        cr.cmd_pregate(A)


# ---------------------------------------------------------------------------------------------------------------------
# The reading on hand-made task files


def _arm(score):
    return {"score": score, "groups": {g: score / 100 for g in ("de_recovery", "direction", "distribution",
                                                               "covariance")}}


def _relabel(cr, **change):
    return {"sha256": cr.rc.STATS_SHA256, "stats_sha256_hashed": cr.rc.STATS_SHA256, "yaml_sha256": cr.rc.YAML_SHA256,
            "route": "recomputed", "record_route": "recomputed", "authors_labels_read": False, **change}


MSHA = "m" * 64


def _jfields(cr, **change):
    return {"manifest_sha256": MSHA, "yaml_sha256": cr.jc.YAML_SHA256, "mode": "practice",
            "embryos_read": sorted(n for e in cr.jc.ENDS for n in cr.jc.ROW_EMBRYOS["practice"][e]),
            "columns_read": ["cell_id", "embryo_id", "day"], "authors_labels_read": False, "stats_npz_read": False,
            "residual_same_values": True, "plc_nan_where_jax_nan": True, **change}


def _constants(cr, **change):
    return {"multipliers": {"x1": 1.0, "x2": 2.0}, "plc_seed": cr.PLC_SEED, "plc_bins": cr.PLC_BINS,
            "blend_weight": 0.5, **change}


def _task(cr, gain=(1.0, 1.0), plc=(0.0, 0.0), lib=False, commit="c" * 40, clean=True, seeds=range(10), **top):
    rows = []
    for s in seeds:
        b1, b2 = 50.0 + s / 100, 51.0 + s / 100
        arms = {cr.BASE1: _arm(b1), cr.ARM1: _arm(b1 + gain[0]), cr.PLC1: _arm(b1 + plc[0]),
                cr.BASE2: _arm(b2), cr.ARM2: _arm(b2 + gain[1]), cr.PLC2: _arm(b2 + plc[1])}
        sha = {a: f"{i}" * 64 for i, a in enumerate(arms)}
        rec = {n: {"non_cm_equal": True, "cm_cells": 30, "cm_cells_changed": 30} for n in cr.COUNTED}
        row = {"seed": s, "arms": arms, "sha256": sha, "records": rec}
        if lib:
            arms["copy_last"], sha["copy_last"] = _arm(40.0), "c" * 64
        rows.append(row)
    out = {"rows": rows, "stamp": {"commit": commit, "tree_clean": clean}, "feng_fraction": 0.2,
           "constants": _constants(cr), "relabel": _relabel(cr), "jax": _jfields(cr), "pregate_sha256": MSHA}
    out.update(top)
    return out


def _gates():
    return {"m1rc_gmom": {k: {s: 50.0 + s / 100 for s in range(10)} for k in ("chamber", "whole", "1", "2", "3")},
            "m1rc_gmom_x2": {k: {s: 51.0 + s / 100 for s in range(10)} for k in ("chamber", "whole", "1", "2", "3")}}


def _pregate(passes=True, commit="c" * 40, sha=MSHA, clean=True):
    return {"reading": {"passes": passes}, "stamp": {"commit": commit, "tree_clean": clean}, "manifest_sha256": sha}


def _loaded(cr, gain=(1.0, 1.0), plc=(0.0, 0.0)):
    out = {"chamber": _task(cr, gain, plc), "whole": _task(cr, gain, plc)}
    out.update({lib: _task(cr, gain, plc, lib=True) for lib in ("1", "2", "3")})
    return out


OK = {"reasons": []}


def _read(cr, loaded, gates=None, pregate=None, sha=MSHA, checks=OK, source=()):
    return cr.t71_reading(loaded, _gates() if gates is None else gates, checks,
                          _pregate() if pregate is None else pregate, sha, source)


def test_ahead_at_both_fractions(cr):
    r = _read(cr, _loaded(cr))
    assert r["outcome"] == "x1: AHEAD; x2: AHEAD" and r["passing"] == ["x1", "x2"]
    lines = cr.outcome_lines(r)
    assert lines[0] == "Outcome at x1: AHEAD" and all(b in lines for b in cr.BESIDE)
    assert "x1, x2" in lines[-1] and "default is no change" in lines[-1]


def test_ahead_at_one_fraction_asks_for_that_fraction_only(cr):
    r = _read(cr, _loaded(cr, gain=(0.1, 1.0)))
    assert r["passing"] == ["x2"] and "(x2)" in r["follows"]


def test_not_ahead_at_both_closes_the_line(cr):
    r = _read(cr, _loaded(cr, gain=(0.1, 0.1)))
    assert r["passing"] == [] and r["follows"] == cr.FOLLOWS_END
    assert "JAX nuclei CM co-direction" in r["follows"]


def _p(mean, ahead=10):
    return {"mean": mean, "ahead": ahead}


def _good(v=1.0):
    return {k: _p(v) for k in ("chamber", "whole", "1", "2", "3")}


@pytest.mark.parametrize("part,ab,ap", [
    ("(i) the chamber", {"chamber": _p(0.29)}, {}),
    ("(i) the chamber", {"chamber": _p(1.0, 7)}, {}),
    ("(i) the whole", {"whole": _p(0.0)}, {}),
    ("(i) the library mean", {"1": _p(0.2), "2": _p(0.2), "3": _p(0.4)}, {}),
    ("(i) library 3", {"3": _p(1.0, 7)}, {}),
    ("(i) libraries", {"2": _p(-0.31), "1": _p(1.0), "3": _p(1.0)}, {}),
    ("(ii) the chamber", {}, {"chamber": _p(0.19)}),
    ("(ii) the chamber", {}, {"chamber": _p(1.0, 7)}),
    ("(ii) the library", {}, {"1": _p(0.1), "2": _p(0.1), "3": _p(0.3)})])
def test_not_ahead_on_each_clause(cr, part, ab, ap):
    ahead, failed, _ = cr.fraction_statement({**_good(), **ab}, {**_good(), **ap})
    assert not ahead and any(f.startswith(part) for f in failed), failed


def test_the_boundaries_that_pass(cr):
    ab = {"chamber": _p(0.3, 8), "whole": _p(1e-9), "1": _p(0.6), "2": _p(-0.3), "3": _p(0.6, 8)}
    ap = {"chamber": _p(0.2, 8), "whole": _p(-5.0), "1": _p(0.2), "2": _p(0.2), "3": _p(0.2, 0)}
    ahead, failed, _ = cr.fraction_statement(ab, ap)
    assert ahead, failed


@pytest.mark.parametrize("case", ["missing_task", "missing_seed", "missing_arm", "missing_sha", "missing_copy_last",
                                  "two_commits", "unclean", "gate_x1", "gate_x2", "non_cm", "manifest", "yaml",
                                  "mode", "column", "authors", "embryo", "pregate_missing", "pregate_failed",
                                  "pregate_commit", "pregate_manifest", "pregate_unclean", "constants", "fraction",
                                  "stats", "route", "checks", "source", "residual", "read_commit"])
def test_every_void_reason_fires(cr, case):
    loaded, gates, pregate, sha, checks, source = _loaded(cr), _gates(), _pregate(), MSHA, OK, []
    read_commit = "c" * 40
    assert cr.t71_reading(loaded, gates, checks, pregate, sha, source, read_commit)["outcome"] != "VOID"
    if case == "missing_task":
        loaded["2"] = None
    elif case == "missing_seed":
        loaded["whole"] = _task(cr, seeds=range(9))
    elif case == "missing_arm":
        del loaded["chamber"]["rows"][2]["arms"][cr.PLC2]
    elif case == "missing_sha":
        del loaded["whole"]["rows"][5]["sha256"][cr.ARM1]
    elif case == "missing_copy_last":
        del loaded["3"]["rows"][0]["arms"]["copy_last"]
    elif case == "two_commits":
        loaded["1"] = _task(cr, lib=True, commit="d" * 40)
    elif case == "unclean":
        loaded["chamber"] = _task(cr, clean=False)
    elif case == "gate_x1":
        loaded["2"]["rows"][4]["arms"][cr.BASE1]["score"] += 0.0201
    elif case == "gate_x2":
        loaded["chamber"]["rows"][9]["arms"][cr.BASE2]["score"] -= 0.0201
    elif case == "non_cm":
        loaded["2"]["rows"][6]["records"][cr.PLC2]["non_cm_equal"] = False
    elif case == "manifest":
        loaded["whole"]["jax"]["manifest_sha256"] = "0" * 64
    elif case == "yaml":
        loaded["chamber"]["jax"]["yaml_sha256"] = "0" * 64
    elif case == "mode":
        loaded["chamber"]["jax"]["mode"] = "test"
    elif case == "column":
        loaded["3"]["jax"]["columns_read"] = ["cell_id", "embryo_id", "day", "celltype_update"]
    elif case == "authors":
        loaded["3"]["jax"]["authors_labels_read"] = True
    elif case == "embryo":
        loaded["1"]["jax"]["embryos_read"] = [*loaded["1"]["jax"]["embryos_read"], 22]
    elif case == "pregate_missing":
        pregate = None
    elif case == "pregate_failed":
        pregate = _pregate(passes=False)
    elif case == "pregate_commit":
        pregate = _pregate(commit="e" * 40)
    elif case == "pregate_manifest":
        pregate = _pregate(sha="0" * 64)
    elif case == "pregate_unclean":
        pregate = _pregate(clean=False)
    elif case == "constants":
        loaded["2"]["constants"] = _constants(cr, plc_seed=1)
    elif case == "fraction":
        loaded["chamber"]["feng_fraction"] = 0.4
    elif case == "stats":
        loaded["chamber"]["relabel"] = _relabel(cr, stats_sha256_hashed="0" * 64)
    elif case == "route":
        loaded["1"]["relabel"] = _relabel(cr, route="none")
    elif case == "checks":
        checks = {"reasons": ["the settings file's sha256 is x"]}
    elif case == "source":
        source = ["the source names a banned token"]
    elif case == "residual":
        loaded["2"]["jax"]["residual_same_values"] = False
    elif case == "read_commit":
        read_commit = "e" * 40
    r = cr.t71_reading(loaded, gates, checks, pregate, sha, source, read_commit)
    assert r["outcome"] == "VOID" and r["reasons"] and r["follows"] == cr.FOLLOWS_VOID, case
    lines = cr.outcome_lines(r)
    assert lines[0].startswith("Outcome: VOID") and not any(b in lines for b in cr.BESIDE)


def test_the_gate_passes_at_a_difference_of_0_02(cr):
    loaded = _loaded(cr)
    loaded["2"]["rows"][4]["arms"][cr.BASE1]["score"] += 0.02 - 1e-12
    loaded["3"]["rows"][1]["arms"][cr.BASE2]["score"] -= 0.02 - 1e-12
    assert _read(cr, loaded)["outcome"] != "VOID"
    assert cr.REPRO_TOL == 0.02 == cr.br.REPRO_TOL


def test_the_counts_and_constants(cr):
    assert cr.TRIES_COUNT == 20 and cr.COUNTED == ("cmjax_x1", "cmjax_plc_x1", "cmjax_x2", "cmjax_plc_x2")
    assert cr.GATED == ("m1rc_gmom", "m1rc_gmom_x2") and cr.LIB_ARMS[0] == "copy_last" and cr.ROW == "T71"
    assert (cr.MIN_END_NUCLEI, cr.MAX_SKM_SHARE, cr.SPEARMAN_LOW, cr.SPEARMAN_HIGH, cr.FENG_MIN_CP,
            cr.MAX_HELD_SHARE) == (200, 0.10, 0.2, 0.95, 1.0, 0.5)
    assert (cr.CHAMBER_BAR, cr.LIB_BAR, cr.LIB_FLOOR, cr.PLC_CHAMBER_BAR, cr.PLC_LIB_BAR, cr.MIN_AHEAD) == (
        0.3, 0.3, -0.3, 0.2, 0.2, 8)
    assert (cr.PRACTICE_FRACTION, cr.TEST_FRACTION, cr.BLEND_WEIGHT) == (0.2, 0.4, 0.5)
    assert cr.T62_TASKS == SCRIPT.parent / "results" / "t62_recomputed" / "tasks"
    assert cr.T63_TASKS == SCRIPT.parent / "results" / "t63_double_step" / "tasks"
    assert all((cr.T62_TASKS / f).exists() and (cr.T63_TASKS / f).exists() for f in cr.task_files().values())


def test_every_output_is_created_exclusively_and_the_source_guard(cr, tmp_path):
    p = tmp_path / "lib1.json"
    cr.write_new(p, "first")
    with pytest.raises(FileExistsError):
        cr.write_new(p, "second")
    text = SCRIPT.read_text()
    assert ".write_text(" not in text and text.count("write_new(") == 7  # the helper and six outputs
    assert cr.source_guard() == [] and cr.rc.source_guard(SCRIPT) == []
    assert "build_test.py" not in text.split('"""', 2)[2] and "boards/index" not in text
    for line in ("import double_step_read\n", "import perm_read\n", "import " + "n5" + "4_relabel_listing\n"):
        bad = tmp_path / "bad.py"
        bad.write_text(line)
        assert cr.source_guard(bad), line


# ---------------------------------------------------------------------------------------------------------------------
# Refusals and the printed lines


def _tries(tmp_path, row="T71"):
    t = tmp_path / "TRIES.md"
    t.write_text(f"| {row} (one row) | T1 | x |\n")
    return t


def _ready(cr, tmp_path, monkeypatch, commit="c" * 40, clean=True):
    monkeypatch.setattr(cr, "TRIES", _tries(tmp_path))
    monkeypatch.setattr(cr.gv, "tree_state", lambda: {"commit": commit, "clean": clean})
    monkeypatch.setattr(cr.rc, "relabel_checks", lambda *a, **k: {"passes": True, "reasons": []})
    man = tmp_path / "MANIFEST.json"
    man.write_text("{}")
    monkeypatch.setattr(cr, "reducer_paths", lambda: (man, tmp_path / "d.npz"))
    return cr.fio.sha256_file(man)


def test_the_runs_and_the_read_are_refused_without_the_row(cr, tmp_path, monkeypatch):
    monkeypatch.setattr(cr, "TRIES", _tries(tmp_path, "T68"))
    monkeypatch.setattr(cr.rc, "relabel_checks", lambda *a, **k: pytest.fail("nothing before the row's refusal"))
    monkeypatch.setattr(cr.gv, "tree_state", lambda: pytest.fail("nothing before the row's refusal"))
    with pytest.raises(SystemExit, match="no row"):
        cr.refuse(tmp_path / "out", "x.json", "c" * 40)
    with pytest.raises(SystemExit, match="no row"):
        cr.read(tmp_path / "out", "c" * 40)


def test_a_run_is_refused_on_another_head_an_unclean_tree_or_without_a_passed_pregate(cr, tmp_path, monkeypatch):
    sha = _ready(cr, tmp_path, monkeypatch)
    out = tmp_path / "out"
    with pytest.raises(SystemExit, match="pre-gate has not passed"):
        cr.refuse(out, "x.json", "c" * 40)
    out.mkdir()
    (out / "pregate.json").write_text(json.dumps(_pregate(commit="e" * 40, sha=sha)))
    with pytest.raises(SystemExit, match="another commit"):
        cr.refuse(out, "x.json", "c" * 40)
    (out / "pregate.json").write_text(json.dumps(_pregate(sha=sha)))
    assert cr.refuse(out, "x.json", "c" * 40)[0]["clean"] is True
    with pytest.raises(SystemExit, match="HEAD"):
        cr.refuse(out, "x.json", "d" * 40)
    (out / "x.json").write_text("{}")
    with pytest.raises(SystemExit, match="already holds"):
        cr.refuse(out, "x.json", "c" * 40)
    _ready(cr, tmp_path, monkeypatch, clean=False)
    with pytest.raises(SystemExit, match="not clean"):
        cr.refuse(out, "y.json", "c" * 40)


def test_the_read_prints_only_words_and_never_overwrites(cr, tmp_path, monkeypatch, capsys):
    sha = _ready(cr, tmp_path, monkeypatch)
    out = tmp_path / "results" / "t71"
    (out / "tasks").mkdir(parents=True)
    loaded = _loaded(cr, gain=(1.234567, 1.0))
    for k, f in cr.task_files().items():
        loaded[k]["jax"]["manifest_sha256"] = loaded[k]["pregate_sha256"] = sha
        (out / "tasks" / f).write_text(json.dumps(loaded[k]))
    (out / "tasks" / "pregate.json").write_text(json.dumps(_pregate(sha=sha)))
    gates = _gates()
    files = {v: k for k, v in cr.task_files().items()}
    monkeypatch.setattr(cr, "committed", lambda f, tasks, arm: gates[arm][files[f]])
    monkeypatch.setattr(cr.gv, "stamp", lambda started, fs: {"commit": "c" * 40, "tree_clean": True, "utc": "x"})
    res = cr.read(out, "c" * 40)
    assert res["reading"]["passing"] == ["x1", "x2"] and res["tries"] == 20
    printed = capsys.readouterr().out
    assert "Outcome at x1: AHEAD" in printed and "1.23" not in printed and "51." not in printed
    assert "+1.23" in (out / "t71.txt").read_text()
    with pytest.raises(SystemExit):
        cr.read(out, "c" * 40)


# ---------------------------------------------------------------------------------------------------------------------
# The launcher, traced under a temporary HOME in a temporary repository


def _git_env(home):
    env = {k: v for k, v in os.environ.items() if not k.startswith("GIT_")}
    env.update(HOME=str(home), GIT_CEILING_DIRECTORIES=str(home.parent), GIT_AUTHOR_NAME="t", GIT_AUTHOR_EMAIL="t@t",
               GIT_COMMITTER_NAME="t", GIT_COMMITTER_EMAIL="t@t")
    return env


def _repo(tmp_path, row="T71"):
    repo, home = tmp_path / "repo", tmp_path / "home"
    (repo / "experiments" / "t1_m1").mkdir(parents=True)
    (repo / "experiments" / "harness").mkdir(parents=True)
    home.mkdir()
    shutil.copy(LAUNCHER, repo / "experiments" / "t1_m1" / "cmjax.sh")
    (repo / "experiments" / "harness" / "TRIES.md").write_text(f"| {row} (one row) | T1 | x |\n")
    (repo / ".gitignore").write_text("data/\n")
    env = _git_env(home)
    for cmd in (["init", "-q"], ["add", "-A"], ["commit", "-q", "-m", "row"]):
        subprocess.run(["git", *cmd], cwd=repo, env=env, check=True)
    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=repo, env=env, check=True, capture_output=True,
                          text=True).stdout.strip()
    return repo, home, env, head


def _run(repo, env, *args, cwd=None):
    r = subprocess.run(["bash", str(repo / "experiments" / "t1_m1" / "cmjax.sh"), *args], cwd=cwd or repo, env=env,
                       capture_output=True, text=True)
    return r.returncode


def _s(home):
    return home / "v" / "scratch" / "t71"


def _refused(home, name):
    return (_s(home) / f"{name}.refused").read_text().strip()


def _no_record(home, name):
    return not (_s(home) / f"{name}.exit").exists() and not (_s(home) / f"{name}.head").exists()


def _passed(home, head, *names):
    for n in names:
        (_s(home) / f"{n}.exit").write_text("0\n")
        (_s(home) / f"{n}.head").write_text(head + "\n")


def test_the_launcher_refusals(tmp_path):
    repo, home, env, head = _repo(tmp_path)
    assert _run(repo, env) == 2 and _refused(home, "usage") == "2"
    assert _run(repo, env, "cut", head) == 2
    assert _run(repo, env, "bogus", head) == 97 and _refused(home, "bogus") == "97"
    assert _run(repo, env, "tests", "f" * 40) == 95 and _no_record(home, "tests")
    assert _run(repo, env, "tests", head[:12]) == 95 and _refused(home, "tests") == "95"  # a short hash refuses
    assert _run(repo, env, "tests", head.upper()) == 95
    assert _run(repo, env, "reduce", head) == 94 and _no_record(home, "reduce")
    _passed(home, head, "tests")
    assert _run(repo, env, "reduce", head) == 93
    (repo / "data" / "derived" / "feng_relabel_c").mkdir(parents=True)
    (repo / "data" / "derived" / "feng_relabel_c" / "cd1_stats_relabel.npz").write_text("x")
    assert _run(repo, env, "reduce", head) == 93
    (repo / "data" / "external" / "jax_qiu2024").mkdir(parents=True)
    (repo / "data" / "external" / "jax_qiu2024" / "MANIFEST.json").write_text("{}")
    (repo / "data" / "derived" / "jax_cm_t71_practice").mkdir()
    assert _run(repo, env, "reduce", head) == 91 and _refused(home, "reduce") == "91"
    assert _run(repo, env, "pregate", head) == 94
    _passed(home, head, "reduce")
    assert _run(repo, env, "pregate", head) == 93
    (repo / "data" / "derived" / "jax_cm_t71_practice" / "MANIFEST.json").write_text("{}")
    tasks = repo / "experiments" / "t1_m1" / "results" / "t71_cmjax" / "tasks"
    tasks.mkdir(parents=True)
    (tasks / "pregate.json").write_text("{}")
    assert _run(repo, env, "pregate", head) == 91
    assert _run(repo, env, "cut", "chamber", head) == 94
    (_s(home) / "pregate.exit").write_text("3\n")
    (_s(home) / "pregate.head").write_text(head + "\n")
    assert _run(repo, env, "lib", "2", head) == 94 and _refused(home, "lib_2") == "94"
    _passed(home, head, "pregate")
    (tasks / "T1_E9.5_E10.5_chamber.json").write_text("{}")
    assert _run(repo, env, "cut", "chamber", head) == 91 and _refused(home, "cut_chamber") == "91"
    assert _run(repo, env, "read", head) == 94
    _passed(home, head, "cut_chamber", "cut_whole", "lib_1", "lib_2", "lib_3")
    assert _run(repo, env, "cut", "chamber", head) == 91
    assert (_s(home) / "cut_chamber.exit").read_text().strip() == "0"  # a refused call never touches a finished exit
    (tasks.parent / "t71.json").write_text("{}")
    assert _run(repo, env, "read", head) == 91
    assert _run(repo, env, "read", head, cwd=repo / "experiments") == 91  # from a subfolder: the checkout's root
    (repo / "stray.txt").write_text("x")
    assert _run(repo, env, "read", head) == 98 and _refused(home, "read") == "98" and _no_record(home, "read")


def test_a_step_that_runs_marks_itself_running_then_writes_its_exit(tmp_path):
    repo, home, env, head = _repo(tmp_path)
    fake = tmp_path / "bin"
    fake.mkdir()
    seen = tmp_path / "seen.txt"
    (fake / "uv").write_text(f'#!/usr/bin/env bash\ncat "$HOME/v/scratch/t71/tests.exit" "$HOME/v/scratch/t71/'
                             f'tests.head" > "{seen}"\npwd >> "{seen}"\nexit 4\n')
    (fake / "uv").chmod(0o755)
    env = {**env, "PATH": f"{fake}:{env['PATH']}"}
    _s(home).mkdir(parents=True)
    (repo / "data" / "derived" / "feng_relabel_c").mkdir(parents=True)
    (repo / "data" / "derived" / "feng_relabel_c" / "cd1_stats_relabel.npz").write_text("x")
    assert _run(repo, env, "tests", head, cwd=repo / "experiments" / "t1_m1") == 4
    assert seen.read_text().split() == ["running", head, str(repo)]
    assert (_s(home) / "tests.exit").read_text().strip() == "4"


def test_the_launcher_refuses_without_the_row(tmp_path):
    repo, home, env, head = _repo(tmp_path, row="T68")
    assert _run(repo, env, "tests", head) == 96 and _refused(home, "tests") == "96"


def test_the_launcher_runs_the_tests_first_and_never_prints_a_log_or_removes_a_file():
    sh = LAUNCHER.read_text()
    tests = sh.split("  tests)\n")[1].split(";;")[0]
    assert "$PY -m pytest" in tests
    for f in ("test_cmjax_read.py", "test_jax_cm.py", "test_jax_blood.py", "test_jaxblood_read.py",
              "test_t62_recomputed_read.py", "test_t63_double_step_read.py", "test_core.py"):
        assert f"tests/{f}" in tests and (ROOT / "tests" / f).exists()
    assert "cmjax_read.py" in sh and "jax_cm.py --mode practice" in sh and "NUMPY_MADVISE_HUGEPAGE=0" in sh
    runs = [ln for ln in sh.splitlines() if "$PY" in ln and ".py" in ln] + [tests]
    assert all('> "$S/$NAME.out" 2> "$S/$NAME.err"' in r for r in runs)
    assert "rm " not in sh and "rmdir" not in sh and ":?" not in sh and "/v/vec-campaign" not in sh
    assert subprocess.run(["bash", "-n", str(LAUNCHER)]).returncode == 0
