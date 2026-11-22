"""Try T68 (experiments/harness/TRIES.md): experiments/t1_m1/jaxblood_read.py and its launcher on synthetic stages,
hand-made task files and temporary git repositories only. No released file, reducer output, statistics file, board file
or score is read; the committed T62 and T64 task folders are checked for their file names only. It reuses row T64's
fixtures."""

from __future__ import annotations

import importlib.util
import inspect
import json
import math
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

from vec import harness, t1_scope
from vec.core import counts_of
from vec.t1_scope import classify, load_spec, reweight

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "experiments" / "t1_m1" / "jaxblood_read.py"
LAUNCHER = ROOT / "experiments" / "t1_m1" / "jaxblood.sh"
SHARES = {"cm": 0.6, "endocardium": 0.25, "blood": 0.15}
DIRECTED = ("CM", "Endocardium", "Endothelium")
HELD_M1 = sorted({*HELD, GENES[3], GENES[7]})
WITHIN = {"surrounding": ["Foregut", "Epithelium"]}


@pytest.fixture
def jr(monkeypatch, tmp_path):
    spec = importlib.util.spec_from_file_location("jaxblood_read_module", SCRIPT)
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


def _jax(n_genes, size=3.0, seed=5):
    return np.random.default_rng(seed).normal(0.0, size, n_genes)


def _fb():
    spec = load_spec()
    return {"E8.5": spec.fallbacks("E8.5_RNA"), "E9.5": spec.fallbacks("E9.5_RNA")}


def _lib(jr, e95, e85, d, d_perm=None, seed=4):
    factors = jr.br.gmom_fits(e95, "E9.5_RNA", e85, "E8.5_RNA")["file"]
    dirs = _directions(e95.n_vars)
    arms, rec = jr.library_arms(e95, SHARES, factors, HELD, HELD_M1, seed, WITHIN, _fb(), dirs, 0.2, d,
                                d if d_perm is None else d_perm, n_cells=300)
    base0, _ = reweight(e95, SHARES, 300, seed, within=WITHIN, fallback=_fb()["E9.5"])
    return arms, rec, base0, factors, dirs


def _cut(jr, d, d_perm=None, seed=1000002):
    ref, start, _ = _small_inputs()
    ctx = _gmom_ctx(ref, start, seed=seed)
    fits = jr.br.gmom_fits(ref, "E9.5_RNA", start, "E8.5_RNA")
    dirs = _directions(ref.n_vars)
    arms, rec = jr.cut_arms(ctx, fits["file"], dirs, 0.2, HELD_M1, d, d if d_perm is None else d_perm)
    return arms, rec, ctx, fits, dirs


def _cut_draw_labels(ctx):
    ref = harness._shallow(ctx.inputs[ctx.task.reference_role])
    src = ctx.task.role_stems[ctx.task.reference_role]
    classify(ref, source=src)
    shares = harness.t1_chamber_shares(ref, src)
    pop, _ = reweight(ref, shares, harness.count_for(ctx.board, ref.n_obs), ctx.seed,
                      fallback=load_spec().fallbacks(src))
    return pop.obs["celltype"].astype(str).to_numpy()


# ---------------------------------------------------------------------------------------------------------------------
# (a) with f = 1 the arms are one file with the base; (b) only Blood cells change


def test_with_f_equal_to_one_the_arms_are_the_base_on_both_routes(jr):
    e95, e85 = _classified()
    arms, rec, _, _, _ = _lib(jr, e95, e85, np.zeros(e95.n_vars))
    shas = {k: jr.hmr.prediction_sha(v.X) for k, v in arms.items()}
    assert shas[jr.ARM] == shas[jr.BASE] == shas[jr.PERM]
    assert rec[jr.ARM]["blood_cells_changed"] == 0 and rec[jr.ARM]["blood_cells"] > 0
    arms, rec, _, _, _ = _cut(jr, np.full(len(GENES), np.nan))  # no JAX value anywhere: f = 1 too
    shas = {k: jr.hmr.prediction_sha(v.X) for k, v in arms.items()}
    assert shas[jr.ARM] == shas[jr.BASE] == shas[jr.PERM]


def test_on_the_library_route_only_blood_cells_change(jr):
    e95, e85 = _classified()
    d = _jax(e95.n_vars)
    arms, rec, base0, _, _ = _lib(jr, e95, e85, d, d_perm=d[::-1].copy())
    blood = base0.obs["celltype"].astype(str).to_numpy() == "Blood"
    a = _dense(counts_of(arms[jr.BASE]))
    for name in jr.COUNTED:
        differs = (a != _dense(counts_of(arms[name]))).any(axis=1)
        assert differs.any() and not differs[~blood].any(), name
        assert rec[name]["non_blood_equal"] is True and rec[name]["non_held_equal"] is True
        assert rec[name]["held_set"] == ["Blood"] and rec[name]["blood_cells_changed"] == int(differs.sum())


def test_on_the_cut_route_only_blood_cells_change(jr):
    arms, rec, ctx, _, _ = _cut(jr, _jax(len(GENES)))
    blood = _cut_draw_labels(ctx) == "Blood"
    for name in jr.COUNTED:
        differs = (_dense(arms[jr.BASE].X) != _dense(arms[name].X)).any(axis=1)
        assert differs.any() and not differs[~blood].any(), name
        assert rec[name]["non_blood_equal"] is True and rec[name]["blood_cells"] == int(blood.sum())
        assert rec[name]["shift"]["rows"][0]["type"] == "Blood"


def test_a_non_blood_cell_that_differs_stops_the_run(jr):
    e95, _ = _classified()
    pop, _ = reweight(e95, SHARES, 300, 4, within=WITHIN)
    X = _dense(pop.X).astype(np.float32)
    Y = X.copy()
    Y[np.flatnonzero(pop.obs["celltype"].astype(str).to_numpy() == "V-CM")[0], 0] += 1.0
    assert jr.blood_record(pop, X, X.copy(), jr.ARM)["non_blood_equal"] is True
    with pytest.raises(SystemExit, match="not labelled Blood"):
        jr.blood_record(pop, X, Y, jr.ARM)
    Z = X.copy()
    Z[np.flatnonzero(pop.obs["celltype"].astype(str).to_numpy() == "Blood")[0], 0] += 1.0
    assert jr.blood_record(pop, X, Z, jr.ARM)["blood_cells_changed"] == 1


# ---------------------------------------------------------------------------------------------------------------------
# (c) the move by hand: fraction, clip, held genes, keep_depth, stream


def test_the_blood_cells_are_external_shift_on_the_blood_rows(jr):
    e95, e85 = _classified()
    d = _jax(e95.n_vars, size=40.0)  # large: every free factor is clipped
    arms, rec, _, _, dirs = _lib(jr, e95, e85, d, seed=6)
    pop, groups, _, _ = jr.library_population(e95, SHARES, HELD_M1, 6, WITHIN, _fb(), dirs, 0.2, n_cells=300)
    blood = np.flatnonzero(pop.obs["celltype"].astype(str).to_numpy() == "Blood")
    want, rep = t1_scope.external_shift(pop[blood].copy(), {"Blood": d}, 1 / 7.5, {"Blood": "Blood"}, HELD_M1,
                                        max_abs_log_factor=math.log(3.0), keep_depth=True, seed=[6, 1])
    got = _dense(counts_of(arms[jr.ARM]))[blood]
    assert np.array_equal(got, _dense(counts_of(want)))
    assert rec[jr.ARM]["shift"]["rows"][0]["max_abs_log_f"] == pytest.approx(math.log(3.0))
    held = np.isin(np.asarray(e95.var_names), HELD_M1)
    assert np.array_equal(got[:, held], _dense(counts_of(pop))[blood][:, held])  # a held gene stays as drawn
    assert jr.FRACTION == 1 / 7.5 and jr.CLIP == math.log(3.0) and jr.JAX_STREAM == 1


def test_the_two_arms_share_one_stream_so_one_direction_gives_one_file(jr):
    e95, e85 = _classified()
    d = _jax(e95.n_vars)
    arms, _, _, _, _ = _lib(jr, e95, e85, d, d_perm=d.copy())
    assert jr.hmr.prediction_sha(arms[jr.ARM].X) == jr.hmr.prediction_sha(arms[jr.PERM].X)


def test_the_base_is_t62s_candidate_as_called_on_both_routes(jr):
    arms, _, ctx, fits, dirs = _cut(jr, _jax(len(GENES)), seed=1000003)
    t62 = jr.rc.cut_pair(ctx, fits["file"], dirs, 0.2, HELD_M1)
    assert arms[jr.BASE].X.tobytes() == t62["m1rc_gmom"].X.tobytes()
    e95, e85 = _classified()
    arms, _, _, factors, dirs = _lib(jr, e95, e85, _jax(e95.n_vars), seed=5)
    t62 = jr.rc.library_pair(e95, SHARES, factors, HELD, HELD_M1, 5, WITHIN, _fb(), dirs, 0.2, n_cells=300)
    assert _dense(arms[jr.BASE].X).tobytes() == _dense(t62["m1rc_gmom"].X).tobytes()
    for fn in (jr.cut_arms, jr.library_arms, jr.library_population):  # no truth reaches an arm
        assert not {"truth", "target", "e105"} & set(inspect.signature(fn).parameters)


# ---------------------------------------------------------------------------------------------------------------------
# (d) the permuted control and the gene map


def test_the_permutation_keeps_values_bins_held_and_unmapped_genes(jr):
    rng = np.random.default_rng(2)
    n = 400
    d, level = rng.normal(size=n), rng.gamma(1.0, 50.0, n)
    d[:7] = np.nan
    held = np.zeros(n, dtype=bool)
    held[10:20] = True
    p = jr.permuted(d, level, held)
    assert np.array_equal(np.isnan(p), np.isnan(d)) and np.array_equal(p[held], d[held])
    ok = np.isfinite(d) & ~held
    assert np.array_equal(np.sort(p[ok]), np.sort(d[ok])) and not np.array_equal(p[ok], d[ok])
    idx = np.flatnonzero(ok)
    order = idx[np.argsort(level[idx], kind="stable")]
    for b in np.array_split(order, 10):  # every value stays in its level bin
        assert np.array_equal(np.sort(p[b]), np.sort(d[b]))
    # one seed, one permutation (the array holds NaN where d has no value, so NaN compares equal here)
    assert np.array_equal(p, jr.permuted(d, level, held), equal_nan=True)
    assert jr.PERM_SEED == 20261122 and jr.PERM_BINS == 10
    assert not np.array_equal(p, jr.permuted(d, level, held, seed=1), equal_nan=True)


def test_the_gene_map_uses_unique_symbols_only(jr):
    arrays = {"symbol": np.array(["A", "B", "B", "C"]), "direction": np.array([1.0, 2.0, 3.0, 4.0])}
    got = jr.on_t1(arrays, ["C", "B", "A", "Z"])
    assert got[0] == 4.0 and np.isnan(got[1]) and got[2] == 1.0 and np.isnan(got[3])


# ---------------------------------------------------------------------------------------------------------------------
# (e) the pre-gate


def _counts(jr, early=40, late=60):
    """Assigned erythroid nuclei per embryo: 12 x 40 = 480 at the early end, 5 x 60 = 300 at the late end."""
    return {n: (early if e == "early" else late) for e in jr.jb.ENDS for n in jr.jb.ROW_EMBRYOS[e]}


def _manifest(jr, assigned=None, entering=None, **change):
    """The pre-gate decides on erythroid_assigned (the specification's count); erythroid_entering (nuclei at the UMI
    target) is reported only, set here below assigned by default."""
    assigned = assigned or _counts(jr)
    entering = entering if entering is not None else {n: max(k - 25, 0) for n, k in assigned.items()}
    per = {str(n): {"erythroid_assigned": k, "erythroid_entering": entering[n], "experimental_batch": ["run_4"],
                    "days": ["E8.5"]} for n, k in assigned.items()}
    man = {"config": {"sha256": jr.jb.YAML_SHA256}, "gate": {"row": "T68", "clean": True},
           "columns_read": list(jr.jb.CELL_WHITELIST), "authors_labels_read": False, "stats_npz_read": False,
           "embryos_read": sorted(assigned), "embryo_lists": {e: list(jr.jb.ROW_EMBRYOS[e]) for e in jr.jb.ENDS},
           "embryos_used": {e: list(jr.jb.ROW_EMBRYOS[e]) for e in jr.jb.ENDS}, "has_direction": True,
           "per_embryo": per}
    man.update(change)
    return man


def _held(n=50, k=5):
    h = np.zeros(n, dtype=bool)
    h[:k] = True
    return h


def test_the_pregate_passes_and_reports(jr):
    d = np.r_[np.full(5, 0.1), np.ones(45)]
    r = jr.pregate_reading(_manifest(jr), d, _held(), [f"g{i}" for i in range(50)])
    assert r["passes"] and r["failed"] == [] and r["held_share"] < 0.5
    assert r["ends"]["early"]["total"] == 480 and r["ends"]["late"]["total"] == 300
    assert r["ends"]["late"]["embryos_at_floor"] == 5 and r["ends"]["late"]["entering_total"] == 175
    assert r["logf_summary"]["moved_genes"] == 45 and r["ends"]["early"]["experimental_batch"] == ["run_4"]


def test_the_pregate_counts_assigned_nuclei_not_the_ones_at_the_umi_target(jr):
    """The specification's pre-gate counts nuclei assigned to the erythroid group; the nuclei that reach the UMI
    target and enter the direction are reported, never counted (blocker 3 of the review)."""
    d = np.r_[np.full(5, 0.1), np.ones(45)]
    few_entering = {n: 2 for e in jr.jb.ENDS for n in jr.jb.ROW_EMBRYOS[e]}
    r = jr.pregate_reading(_manifest(jr, entering=few_entering), d, _held(), [f"g{i}" for i in range(50)])
    assert r["passes"], r["failed"]
    assert r["ends"]["late"]["entering_total"] == 10 and r["ends"]["late"]["total"] == 300
    many_entering = {n: 500 for e in jr.jb.ENDS for n in jr.jb.ROW_EMBRYOS[e]}
    r = jr.pregate_reading(_manifest(jr, assigned=_counts(jr, late=59), entering=many_entering), d, _held(),
                           [f"g{i}" for i in range(50)])
    assert not r["passes"] and r["failed"] == ["fewer than 300 assigned erythroid nuclei at the late end"]


@pytest.mark.parametrize("case,part", [
    ("few_late", "fewer than 300 assigned erythroid nuclei at the late end"),
    ("few_embryos", "fewer than 3 embryos with at least 30 assigned erythroid nuclei at the late end"),
    ("held", "R2's held genes"),
    ("idle", "R2's held genes"),
    ("embryo", "outside the row's lists"),
    ("column", "outside the whitelist"),
    ("yaml", "another yaml"),
    ("no_direction", "wrote no direction"),
])
def test_the_pregate_fails(jr, case, part):
    assigned = _counts(jr)
    d, man = np.ones(50), None
    if case == "few_late":
        assigned.update({59: 60, 60: 60, 61: 60, 62: 60, 63: 59})
    elif case == "few_embryos":
        assigned.update({59: 29, 60: 29, 61: 29, 62: 200, 63: 200})
    elif case == "held":
        d = np.r_[np.full(5, 3.0), np.ones(45)]  # 45 / 90 on held genes: not below 0.5
    elif case == "idle":
        d = np.zeros(50)
    man = _manifest(jr, assigned)
    if case == "no_direction":
        man["has_direction"] = False
    if case == "embryo":
        man["embryos_read"] = [*man["embryos_read"], 22]
    elif case == "column":
        man["columns_read"] = [*man["columns_read"], "celltype_update"]
    elif case == "yaml":
        man["config"]["sha256"] = "0" * 64
    r = jr.pregate_reading(man, d, _held(), [f"g{i}" for i in range(50)])
    assert not r["passes"] and any(part in f for f in r["failed"]), r["failed"]
    assert len(r["failed"]) == 1, r["failed"]  # the fixture fails on the case's own reason only
    assert "FAILS" in jr.pregate_outcome(r)


def test_the_pregate_boundaries_that_pass(jr):
    assigned = {n: 23 for n in jr.jb.ROW_EMBRYOS["early"]}
    assigned.update({1: 30, 2: 30, 3: 30, 4: 26, 59: 30, 60: 30, 61: 30, 62: 210, 63: 0})
    d = np.r_[np.full(5, 2.999), np.ones(45)]
    r = jr.pregate_reading(_manifest(jr, assigned), d, _held(), [f"g{i}" for i in range(50)])
    assert r["passes"], r["failed"]  # 300 at each end, three embryos at 30, the share just below 0.5
    assert r["ends"]["early"]["total"] == 300 == r["ends"]["late"]["total"] and r["held_share"] < 0.5


# ---------------------------------------------------------------------------------------------------------------------
# The reading on hand-made task files


def _arm(score):
    return {"score": score, "groups": {g: score / 100 for g in ("de_recovery", "direction", "distribution",
                                                               "covariance")}}


def _relabel(jr, **change):
    return {"sha256": jr.rc.STATS_SHA256, "stats_sha256_hashed": jr.rc.STATS_SHA256, "yaml_sha256": jr.rc.YAML_SHA256,
            "route": "recomputed", "record_route": "recomputed", "authors_labels_read": False, **change}


MSHA = "m" * 64


def _jfields(jr, **change):
    return {"manifest_sha256": MSHA, "yaml_sha256": jr.jb.YAML_SHA256,
            "embryos_read": sorted(n for e in jr.jb.ENDS for n in jr.jb.ROW_EMBRYOS[e]),
            "columns_read": list(jr.jb.CELL_WHITELIST), "authors_labels_read": False, "stats_npz_read": False,
            "perm_same_values": True, **change}


def _constants(jr, **change):
    return {"fraction": jr.FRACTION, "clip": jr.CLIP, "jax_stream": jr.JAX_STREAM, "perm_seed": jr.PERM_SEED,
            "perm_bins": jr.PERM_BINS, **change}


def _task(jr, gain, perm=0.0, lib=False, per_seed=None, perm_seed=None, commit="c" * 40, clean=True,
          seeds=range(10), idle=False, inverted=False, **top):
    rows = []
    for s in seeds:
        g = per_seed[s] if per_seed is not None else gain
        q = perm_seed[s] if perm_seed is not None else perm
        base = 50.0 + s / 100
        arms = {jr.BASE: _arm(base), jr.ARM: _arm(base + g), jr.PERM: _arm(base + q)}
        sha = {jr.BASE: "a" * 64, jr.ARM: "a" * 64 if idle else "b" * 64, jr.PERM: "d" * 64}
        rec = {name: {"held_set": ["Blood"], "uncovered_set": ["Blood"], "non_held_equal": True,
                      "non_blood_equal": True, "blood_cells": 30, "blood_cells_changed": 0 if idle else 30}
               for name in jr.COUNTED}
        row = {"seed": s, "arms": arms, "sha256": sha, "records": rec}
        if lib:
            arms["copy_last"], sha["copy_last"] = _arm(40.0), "c" * 64
            row["inverted_anchors"] = ["x"] if inverted and s == 3 else []
        rows.append(row)
    out = {"rows": rows, "stamp": {"commit": commit, "tree_clean": clean}, "feng_fraction": 0.2,
           "constants": _constants(jr), "relabel": _relabel(jr), "jax": _jfields(jr), "pregate_sha256": MSHA}
    out.update(top)
    return out


def _t62():
    return {k: {s: 50.0 + s / 100 for s in range(10)} for k in ("chamber", "whole", "1", "2", "3")}


def _t64(gain=0.5):
    return {k: {s: 50.0 + s / 100 + gain for s in range(10)} for k in ("1", "2", "3")}


def _pregate(passes=True, commit="c" * 40, sha=MSHA, clean=True):
    return {"reading": {"passes": passes}, "stamp": {"commit": commit, "tree_clean": clean}, "manifest_sha256": sha}


def _loaded(jr, cuts=(0.0, 0.0), libs=(1.0, 1.0, 1.0), perm=0.2, **kw):
    out = {"chamber": _task(jr, cuts[0], perm), "whole": _task(jr, cuts[1], perm)}
    out.update({lib: _task(jr, g, perm, lib=True) for lib, g in zip(("1", "2", "3"), libs, strict=True)})
    out.update(kw)
    return out


OK = {"reasons": []}


def _read(jr, loaded, t64=None, pregate=None, sha=MSHA, checks=OK, source=()):
    return jr.t68_reading(loaded, _t62(), _t64() if t64 is None else t64, checks,
                          _pregate() if pregate is None else pregate, sha, source)


def test_ahead(jr):
    r = _read(jr, _loaded(jr))
    assert r["outcome"] == "AHEAD" and r["ahead"] and r["follows"] == jr.FOLLOWS_AHEAD and r["failed"] == []
    assert r["values"]["libraries_at_or_above_t64"] == ["1", "2", "3"]
    lines = jr.outcome_lines(r)
    assert lines[0] == "Outcome: AHEAD" and all(b in lines for b in jr.BESIDE) and "keep R1" in lines[-1]


def _p(mean, ahead=10):
    return {"mean": mean, "ahead": ahead}


GOOD = {"1": _p(1.0), "2": _p(1.0), "3": _p(1.0), "chamber": _p(0.0), "whole": _p(0.0)}
MEANS = {"1": 51.0, "2": 51.0, "3": 51.0}
T64M = {"1": 50.5, "2": 50.5, "3": 50.5}


@pytest.mark.parametrize("arm,ap,means,part", [
    ({"1": _p(0.29), "2": _p(0.3), "3": _p(0.3)}, {}, MEANS, "(i) the library-task mean"),
    ({"1": _p(1.0, 7)}, {}, MEANS, "(i) library 1 is ahead"),
    ({"3": _p(1.0, 7)}, {}, MEANS, "(i) library 3 is ahead"),
    ({}, {"1": _p(0.29), "2": _p(0.3), "3": _p(0.3)}, MEANS, "(ii) the library-task mean"),
    ({}, {"3": _p(1.0, 7)}, MEANS, "(ii) library 3 is ahead"),
    ({}, {}, {"1": 51.0, "2": 50.49, "3": 50.0}, "(iii)"),
    ({"chamber": _p(-0.0501)}, {}, MEANS, "(iv) the chamber cut"),
    ({"whole": _p(-0.06)}, {}, MEANS, "(iv) the whole cut"),
])
def test_not_ahead_on_each_clause(jr, arm, ap, means, part):
    ahead, failed, _ = jr.statement({**GOOD, **arm}, {**GOOD, **ap}, means, T64M)
    assert not ahead and len(failed) == 1 and failed[0].startswith(part), failed


def test_the_boundaries_that_pass(jr):
    # values whose library means are exactly the float 0.3 (0.9 / 3), so each clause sits on its nominal boundary
    arm = {"1": _p(0.5, 8), "2": _p(-0.5, 0), "3": _p(0.9, 8), "chamber": _p(-0.05), "whole": _p(-0.05)}
    ap = {"1": _p(0.0, 0), "2": _p(0.0, 0), "3": _p(0.9, 8)}
    ahead, failed, values = jr.statement(arm, ap, {"1": 50.5, "2": 50.5, "3": 0.0}, T64M)
    assert values["library_mean_d"] == 0.3 == values["library_mean_arm_minus_perm"]
    assert (ahead, failed) == (True, [])
    # library 2 behind on every seed and below T64 is allowed by the clauses as fixed
    below = {**arm, "3": _p(0.9 - 1e-12, 8)}  # just under the boundary fails (i): no tolerance is added
    assert jr.statement(below, ap, {"1": 50.5, "2": 50.5, "3": 0.0}, T64M)[1][0].startswith("(i)")


def test_not_ahead_end_to_end_closes_the_line(jr):
    r = _read(jr, _loaded(jr, perm=0.9))  # the arm is no better than its placebo: (ii) fails
    assert r["outcome"] == "NOT AHEAD" and r["follows"] == jr.FOLLOWS_END
    assert any(f.startswith("(ii)") for f in r["failed"])
    assert jr.outcome_lines(r)[-1] == jr.FOLLOWS_END and "closes" in jr.FOLLOWS_END


def test_clause_iii_reads_t64s_committed_means(jr):
    assert _read(jr, _loaded(jr), t64=_t64(gain=0.99))["outcome"] == "AHEAD"
    r = _read(jr, _loaded(jr), t64=_t64(gain=1.01))
    assert r["outcome"] == "NOT AHEAD" and r["failed"][0].startswith("(iii)")


def test_the_idle_case_is_named(jr):
    loaded = _loaded(jr)
    loaded["whole"] = _task(jr, 0.0, 0.2, idle=True)
    r = _read(jr, loaded)
    assert r["idle_tasks"] == ["whole"] and "the arm moved no cell on whole: nothing is priced there" in \
        jr.outcome_lines(r)


@pytest.mark.parametrize("case", ["missing_task", "missing_seed", "missing_arm", "missing_perm", "missing_sha",
                                  "missing_copy_last", "two_commits", "unclean", "gate", "non_blood", "non_held",
                                  "held_set", "manifest", "manifest_now", "yaml", "column", "authors", "embryo",
                                  "pregate_missing", "pregate_failed", "pregate_commit", "pregate_manifest",
                                  "constants", "fraction", "stats", "route", "checks", "source", "t64_seed",
                                  "perm_values", "held_not_blood", "pregate_unclean", "read_commit"])
def test_every_void_reason_fires(jr, case):
    loaded, t64, pregate, sha, checks, source = _loaded(jr), _t64(), _pregate(), MSHA, OK, []
    read_commit = "c" * 40
    assert jr.t68_reading(loaded, _t62(), t64, checks, pregate, sha, source, read_commit)["outcome"] != "VOID"
    if case == "held_not_blood":  # another R1-held label drawn: the arm would move Blood alone
        for name in jr.COUNTED:
            loaded["2"]["rows"][3]["records"][name].update(held_set=["Blood", "Macrophage"],
                                                           uncovered_set=["Blood", "Macrophage"])
    elif case == "pregate_unclean":
        pregate = _pregate(clean=False)
    elif case == "read_commit":
        read_commit = "e" * 40
    if case == "missing_task":
        loaded["2"] = None
    elif case == "missing_seed":
        loaded["whole"] = _task(jr, 0.0, seeds=range(9))
    elif case == "missing_arm":
        del loaded["chamber"]["rows"][2]["arms"][jr.ARM]
    elif case == "missing_perm":
        del loaded["1"]["rows"][2]["arms"][jr.PERM]
    elif case == "missing_sha":
        del loaded["whole"]["rows"][5]["sha256"][jr.PERM]
    elif case == "missing_copy_last":
        del loaded["3"]["rows"][0]["arms"]["copy_last"]
    elif case == "two_commits":
        loaded["1"] = _task(jr, 1.0, 0.2, lib=True, commit="d" * 40)
    elif case == "unclean":
        loaded["chamber"] = _task(jr, 0.0, 0.2, clean=False)
    elif case == "gate":
        loaded["2"]["rows"][4]["arms"][jr.BASE]["score"] += 0.0201
    elif case == "non_blood":
        loaded["2"]["rows"][6]["records"][jr.PERM]["non_blood_equal"] = False
    elif case == "non_held":
        loaded["3"]["rows"][1]["records"][jr.ARM]["non_held_equal"] = False
    elif case == "held_set":
        loaded["1"]["rows"][1]["records"][jr.ARM]["uncovered_set"] = ["Blood", "Foregut"]
    elif case == "manifest":
        loaded["whole"]["jax"]["manifest_sha256"] = "0" * 64
    elif case == "manifest_now":
        sha = "0" * 64
    elif case == "yaml":
        loaded["chamber"]["jax"]["yaml_sha256"] = "0" * 64
    elif case == "column":
        loaded["3"]["jax"]["columns_read"] = [*jr.jb.CELL_WHITELIST, "celltype_update"]
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
    elif case == "constants":
        loaded["2"]["constants"] = _constants(jr, perm_seed=1)
    elif case == "fraction":
        loaded["chamber"]["feng_fraction"] = 0.4
    elif case == "stats":
        loaded["chamber"]["relabel"] = _relabel(jr, stats_sha256_hashed="0" * 64)
    elif case == "route":
        loaded["1"]["relabel"] = _relabel(jr, route="none")
    elif case == "checks":
        checks = {"reasons": ["the settings file's sha256 is x"]}
    elif case == "source":
        source = ["the source names a banned token"]
    elif case == "t64_seed":
        del t64["2"][7]
    elif case == "perm_values":
        loaded["2"]["jax"]["perm_same_values"] = False
    r = jr.t68_reading(loaded, _t62(), t64, checks, pregate, sha, source, read_commit)
    assert r["outcome"] == "VOID" and r["reasons"] and r["follows"] == jr.FOLLOWS_VOID, case
    lines = jr.outcome_lines(r)
    assert lines[0].startswith("Outcome: VOID") and not any(b in lines for b in jr.BESIDE)


def test_the_gate_passes_at_a_difference_of_0_02(jr):
    loaded = _loaded(jr)
    loaded["2"]["rows"][4]["arms"][jr.BASE]["score"] += 0.02 - 1e-12
    assert _read(jr, loaded)["outcome"] != "VOID"
    assert jr.REPRO_TOL == 0.02 == jr.br.REPRO_TOL


def test_the_counts_and_constants(jr):
    assert jr.TRIES_COUNT == 10 and jr.COUNTED == ("m1rc_gmom_jaxblood", "m1rc_gmom_jaxperm")
    assert jr.ARMS == ("m1rc_gmom", "m1rc_gmom_jaxblood", "m1rc_gmom_jaxperm") and jr.LIB_ARMS[0] == "copy_last"
    assert jr.BASE == jr.rc.CANDIDATE and jr.T64_ARM == "m1rc_gmom_heldmom" and jr.ROW == "T68"
    assert (jr.MIN_END_NUCLEI, jr.MIN_EMBRYO_NUCLEI, jr.MIN_EMBRYOS, jr.MAX_HELD_SHARE) == (300, 30, 3, 0.5)
    assert (jr.LIB_BAR, jr.MIN_AHEAD, jr.MIN_T64_LIBS, jr.CUT_FLOOR) == (0.3, 8, 2, -0.05)
    assert jr.AHEAD_LIBS_ARM == ("1", "3") and jr.AHEAD_LIBS_PERM == ("3",)
    assert jr.T62_TASKS == SCRIPT.parent / "results" / "t62_recomputed" / "tasks"
    assert jr.T64_TASKS == SCRIPT.parent / "results" / "t64_heldmom" / "tasks"
    assert all((jr.T62_TASKS / f).exists() and (jr.T64_TASKS / f).exists() for f in jr.task_files().values())


def test_every_output_is_created_exclusively(jr, tmp_path):
    """Two copies of one step: the later finisher stops and the first's file stays (the review's lesser point)."""
    p = tmp_path / "lib1.json"
    jr.write_new(p, "first")
    with pytest.raises(FileExistsError):
        jr.write_new(p, "second")
    assert p.read_text() == "first"
    text = SCRIPT.read_text()
    assert ".write_text(" not in text and text.count("write_new(") == 7  # the helper and six outputs


def test_the_source_guard(jr, tmp_path):
    assert jr.source_guard() == [] and jr.rc.source_guard(SCRIPT) == []
    text = SCRIPT.read_text()
    assert "build_test.py" not in text.split('"""', 2)[2] and "boards/index" not in text  # it meets no test board
    bad = tmp_path / "bad.py"
    bad.write_text("import perm_read\n")
    assert jr.source_guard(bad)


# ---------------------------------------------------------------------------------------------------------------------
# Refusals and the printed lines


def _tries(tmp_path, row="T68"):
    t = tmp_path / "TRIES.md"
    t.write_text(f"| {row} (one row) | T1 | x |\n")
    return t


def test_the_runs_and_the_read_are_refused_without_the_row(jr, tmp_path, monkeypatch):
    monkeypatch.setattr(jr, "TRIES", _tries(tmp_path, "T64"))
    monkeypatch.setattr(jr.rc, "relabel_checks", lambda *a, **k: pytest.fail("nothing before the row's refusal"))
    monkeypatch.setattr(jr.gv, "tree_state", lambda: pytest.fail("nothing before the row's refusal"))
    with pytest.raises(SystemExit, match="no row"):
        jr.refuse(tmp_path / "out", "x.json", "c" * 40)
    with pytest.raises(SystemExit, match="no row"):
        jr.read(tmp_path / "out", "c" * 40)

    class A:
        out, commit = "o", "c" * 40
    with pytest.raises(SystemExit, match="no row"):
        jr.cmd_pregate(A)


def _ready(jr, tmp_path, monkeypatch, commit="c" * 40, clean=True):
    monkeypatch.setattr(jr, "TRIES", _tries(tmp_path))
    monkeypatch.setattr(jr.gv, "tree_state", lambda: {"commit": commit, "clean": clean})
    monkeypatch.setattr(jr.rc, "relabel_checks", lambda *a, **k: {"passes": True, "reasons": []})
    man = tmp_path / "MANIFEST.json"
    man.write_text("{}")
    monkeypatch.setattr(jr, "reducer_paths", lambda: (man, tmp_path / "d.npz"))
    return jr.fio.sha256_file(man)


def test_a_run_is_refused_on_another_head_an_unclean_tree_or_without_a_passed_pregate(jr, tmp_path, monkeypatch):
    sha = _ready(jr, tmp_path, monkeypatch)
    out = tmp_path / "out"
    with pytest.raises(SystemExit, match="pre-gate has not passed"):
        jr.refuse(out, "x.json", "c" * 40)
    out.mkdir()
    (out / "pregate.json").write_text(json.dumps(_pregate(passes=False, sha=sha)))
    with pytest.raises(SystemExit, match="pre-gate has not passed"):
        jr.refuse(out, "x.json", "c" * 40)
    (out / "pregate.json").write_text(json.dumps(_pregate(commit="e" * 40, sha=sha)))
    with pytest.raises(SystemExit, match="another commit"):
        jr.refuse(out, "x.json", "c" * 40)
    (out / "pregate.json").write_text(json.dumps(_pregate(sha="0" * 64)))
    with pytest.raises(SystemExit, match="another reducer manifest"):
        jr.refuse(out, "x.json", "c" * 40)
    (out / "pregate.json").write_text(json.dumps(_pregate(sha=sha)))
    assert jr.refuse(out, "x.json", "c" * 40)[0]["clean"] is True
    with pytest.raises(SystemExit, match="HEAD"):
        jr.refuse(out, "x.json", "d" * 40)
    (out / "x.json").write_text("{}")
    with pytest.raises(SystemExit, match="already holds"):
        jr.refuse(out, "x.json", "c" * 40)
    _ready(jr, tmp_path, monkeypatch, clean=False)
    with pytest.raises(SystemExit, match="not clean"):
        jr.refuse(out, "y.json", "c" * 40)


def test_the_read_prints_no_score_and_never_overwrites(jr, tmp_path, monkeypatch):
    sha = _ready(jr, tmp_path, monkeypatch)
    out = tmp_path / "results" / "t68"
    (out / "tasks").mkdir(parents=True)
    loaded = _loaded(jr, libs=(1.234567, 1.0, 1.0))
    for k, f in jr.task_files().items():
        loaded[k]["jax"]["manifest_sha256"] = loaded[k]["pregate_sha256"] = sha
        (out / "tasks" / f).write_text(json.dumps(loaded[k]))
    (out / "tasks" / "pregate.json").write_text(json.dumps(_pregate(sha=sha)))
    t62, t64 = _t62(), _t64()
    files = {v: k for k, v in jr.task_files().items()}
    monkeypatch.setattr(jr, "committed", lambda f, tasks, arm: t62[files[f]] if arm == jr.BASE else t64[files[f]])
    monkeypatch.setattr(jr.gv, "stamp", lambda started, fs: {"commit": "c" * 40, "tree_clean": True, "utc": "x"})
    res = jr.read(out, "c" * 40)
    assert res["reading"]["outcome"] == "AHEAD" and res["tries"] == 10
    text = (out / "t68.txt").read_text()
    assert "+1.23" in text and "Clause values" in text and "row T64's m1rc_gmom_heldmom" in text
    with pytest.raises(SystemExit):
        jr.read(out, "c" * 40)  # no second read into the same folder


def test_the_read_prints_only_words(jr, tmp_path, monkeypatch, capsys):
    sha = _ready(jr, tmp_path, monkeypatch)
    out = tmp_path / "results" / "t68w"
    (out / "tasks").mkdir(parents=True)
    loaded = _loaded(jr, libs=(1.234567, 1.0, 1.0))
    for k, f in jr.task_files().items():
        loaded[k]["jax"]["manifest_sha256"] = loaded[k]["pregate_sha256"] = sha
        (out / "tasks" / f).write_text(json.dumps(loaded[k]))
    (out / "tasks" / "pregate.json").write_text(json.dumps(_pregate(sha=sha)))
    t62, t64 = _t62(), _t64()
    files = {v: k for k, v in jr.task_files().items()}
    monkeypatch.setattr(jr, "committed", lambda f, tasks, arm: t62[files[f]] if arm == jr.BASE else t64[files[f]])
    monkeypatch.setattr(jr.gv, "stamp", lambda started, fs: {"commit": "c" * 40, "tree_clean": True, "utc": "x"})
    jr.read(out, "c" * 40)
    printed = capsys.readouterr().out
    assert "Outcome: AHEAD" in printed and "1.23" not in printed and "51." not in printed and "50." not in printed


def test_the_pregate_prints_only_words_and_exits_3_on_a_failure(jr, tmp_path, monkeypatch, capsys):
    _ready(jr, tmp_path, monkeypatch)
    man = _manifest(jr, {n: 10 for e in jr.jb.ENDS for n in jr.jb.ROW_EMBRYOS[e]})
    arrays = {"symbol": np.array(["Xist", "g1", "g2"]), "direction": np.array([0.1, 1.0, -1.0]),
              "m_early": np.array([1.0, 2.0, 3.0])}
    monkeypatch.setattr(jr, "load_reducer", lambda *a: (man, arrays, "s" * 64))
    monkeypatch.setattr(jr.gv, "stamp", lambda started, fs: {"commit": "c" * 40, "tree_clean": True, "utc": "x"})

    class A:
        out, commit = "pg", "c" * 40
    assert jr.cmd_pregate(A) == jr.PREGATE_FAIL == 3
    printed = capsys.readouterr().out
    assert printed.startswith("Pre-gate: FAILS") and "no try is counted" in printed
    assert not any(ch.isdigit() for ch in printed.replace("T68", ""))
    rec = json.loads((jr.RESULTS / "pg" / "tasks" / "pregate.json").read_text())
    assert rec["reading"]["passes"] is False and rec["manifest_sha256"] == "s" * 64
    with pytest.raises(SystemExit, match="already holds"):
        jr.cmd_pregate(A)


# ---------------------------------------------------------------------------------------------------------------------
# The launcher, traced under a temporary HOME in a temporary repository


def _git_env(home):
    env = {k: v for k, v in os.environ.items() if not k.startswith("GIT_")}
    env.update(HOME=str(home), GIT_CEILING_DIRECTORIES=str(home.parent), GIT_AUTHOR_NAME="t", GIT_AUTHOR_EMAIL="t@t",
               GIT_COMMITTER_NAME="t", GIT_COMMITTER_EMAIL="t@t")
    return env


def _repo(tmp_path, row="T68"):
    repo, home = tmp_path / "repo", tmp_path / "home"
    (repo / "experiments" / "t1_m1").mkdir(parents=True)
    (repo / "experiments" / "harness").mkdir(parents=True)
    home.mkdir()
    shutil.copy(LAUNCHER, repo / "experiments" / "t1_m1" / "jaxblood.sh")
    (repo / "experiments" / "harness" / "TRIES.md").write_text(f"| {row} (one row) | T1 | x |\n")
    (repo / ".gitignore").write_text("data/\n")
    env = _git_env(home)
    for cmd in (["init", "-q"], ["add", "-A"], ["commit", "-q", "-m", "row"]):
        subprocess.run(["git", *cmd], cwd=repo, env=env, check=True)
    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=repo, env=env, check=True, capture_output=True,
                          text=True).stdout.strip()
    return repo, home, env, head


def _run(repo, env, *args):
    r = subprocess.run(["bash", "experiments/t1_m1/jaxblood.sh", *args], cwd=repo, env=env, capture_output=True,
                       text=True)
    return r.returncode


def _exit(home, name):
    return (home / "v" / "scratch" / "t68" / f"{name}.exit").read_text().strip()


def _refused(home, name):
    return (home / "v" / "scratch" / "t68" / f"{name}.refused").read_text().strip()


def _no_record(home, name):
    s = home / "v" / "scratch" / "t68"
    return not (s / f"{name}.exit").exists() and not (s / f"{name}.head").exists()


def _passed(home, head, *names):
    s = home / "v" / "scratch" / "t68"
    for n in names:
        (s / f"{n}.exit").write_text("0\n")
        (s / f"{n}.head").write_text(head + "\n")


def test_the_launcher_refusals(tmp_path):
    repo, home, env, head = _repo(tmp_path)
    assert _run(repo, env) == 2 and _refused(home, "usage") == "2"
    assert _run(repo, env, "cut", head) == 2
    assert _run(repo, env, "bogus", head) == 97 and _refused(home, "bogus") == "97"
    assert _run(repo, env, "tests", "f" * 40) == 95 and _refused(home, "tests") == "95"
    assert _no_record(home, "tests")  # a refusal writes no exit or head record
    assert _run(repo, env, "reduce", head) == 94 and _refused(home, "reduce") == "94"  # no tests.exit
    assert _no_record(home, "reduce")
    _passed(home, head, "tests")
    assert _run(repo, env, "reduce", head) == 93  # the relabelling's output is not in place
    (repo / "data" / "derived" / "feng_relabel_c").mkdir(parents=True)
    (repo / "data" / "derived" / "feng_relabel_c" / "cd1_stats_relabel.npz").write_text("x")
    assert _run(repo, env, "reduce", head) == 93  # the JAX folder is not in place
    (repo / "data" / "external" / "jax_qiu2024").mkdir(parents=True)
    (repo / "data" / "external" / "jax_qiu2024" / "MANIFEST.json").write_text("{}")
    (repo / "data" / "derived" / "jax_blood_t68").mkdir()
    assert _run(repo, env, "reduce", head) == 91 and _refused(home, "reduce") == "91"  # never over a folder
    assert _run(repo, env, "pregate", head) == 94
    _passed(home, head, "reduce")
    assert _run(repo, env, "pregate", head) == 93  # no reducer manifest
    (repo / "data" / "derived" / "jax_blood_t68" / "MANIFEST.json").write_text("{}")
    tasks = repo / "experiments" / "t1_m1" / "results" / "t68_jaxblood" / "tasks"
    tasks.mkdir(parents=True)
    (tasks / "pregate.json").write_text("{}")  # untracked inside the row's folder: the tree stays clean
    assert _run(repo, env, "pregate", head) == 91
    assert _run(repo, env, "cut", "chamber", head) == 94  # the pre-gate has not exited 0
    (home / "v" / "scratch" / "t68" / "pregate.exit").write_text("3\n")
    (home / "v" / "scratch" / "t68" / "pregate.head").write_text(head + "\n")
    assert _run(repo, env, "lib", "2", head) == 94 and _refused(home, "lib_2") == "94"  # a failed pre-gate stops all
    _passed(home, head, "pregate")
    (tasks / "T1_E9.5_E10.5_chamber.json").write_text("{}")
    assert _run(repo, env, "cut", "chamber", head) == 91 and _refused(home, "cut_chamber") == "91"
    assert _run(repo, env, "read", head) == 94
    _passed(home, head, "cut_chamber", "cut_whole", "lib_1", "lib_2", "lib_3")
    # a second call of a step that succeeded is refused and leaves its record as it was (the review's lesser point)
    assert _run(repo, env, "cut", "chamber", head) == 91 and _exit(home, "cut_chamber") == "0"
    assert (home / "v" / "scratch" / "t68" / "cut_chamber.head").read_text().strip() == head
    (tasks.parent / "t68.json").write_text("{}")
    assert _run(repo, env, "read", head) == 91 and _refused(home, "read") == "91"  # past the 94 check: records hold
    (repo / "stray.txt").write_text("x")
    assert _run(repo, env, "read", head) == 98 and _refused(home, "read") == "98"
    assert _no_record(home, "read")


def test_a_step_that_runs_marks_itself_running_then_writes_its_exit(tmp_path):
    """A fake `uv` on PATH stands in for the step's Python: the launcher writes 'running' and HEAD before the step
    runs, so no later step can see an earlier exit 0 while it runs, then the step's own exit code."""
    repo, home, env, head = _repo(tmp_path)
    fake = tmp_path / "bin"
    fake.mkdir()
    seen = tmp_path / "seen.txt"
    (fake / "uv").write_text(f'#!/usr/bin/env bash\ncat "$HOME/v/scratch/t68/tests.exit" "$HOME/v/scratch/t68/'
                             f'tests.head" > "{seen}"\nexit 0\n')
    (fake / "uv").chmod(0o755)
    env = {**env, "PATH": f"{fake}:{env['PATH']}"}
    s = home / "v" / "scratch" / "t68"
    s.mkdir(parents=True)
    (repo / "data" / "derived" / "feng_relabel_c").mkdir(parents=True)
    (repo / "data" / "derived" / "feng_relabel_c" / "cd1_stats_relabel.npz").write_text("x")
    (s / "tests.exit").write_text("1\n")  # an earlier failed run at another commit
    (s / "tests.head").write_text("f" * 40 + "\n")
    assert _run(repo, env, "tests", head) == 0
    assert seen.read_text().split() == ["running", head]
    assert _exit(home, "tests") == "0" and (s / "tests.head").read_text().strip() == head


def test_the_launcher_refuses_without_the_row_and_on_another_heads_exit_file(tmp_path):
    repo, home, env, head = _repo(tmp_path, row="T64")
    assert _run(repo, env, "tests", head) == 96 and _refused(home, "tests") == "96"
    repo2, home2, env2, head2 = _repo(tmp_path / "b")
    s = home2 / "v" / "scratch" / "t68"
    s.mkdir(parents=True)
    (s / "tests.exit").write_text("0\n")
    (s / "tests.head").write_text("f" * 40 + "\n")  # exit 0, but at another commit
    assert _run(repo2, env2, "reduce", head2) == 94


def test_the_launcher_runs_the_tests_first_and_never_prints_a_log_or_removes_a_file():
    sh = LAUNCHER.read_text()
    tests = sh.split("  tests)\n")[1].split(";;")[0]  # the step body, not the line `  tests) PREV="" ;;`
    assert "$PY -m pytest" in tests
    for f in ("test_jaxblood_read.py", "test_jax_blood.py", "test_t64_heldmom_read.py", "test_t62_recomputed_read.py",
              "test_t1_bets.py", "test_t1_scope.py", "test_core.py", "test_g_mom_verify.py",
              "test_t1_state_transport.py"):
        assert f"tests/{f}" in tests and (ROOT / "tests" / f).exists()
    assert "jaxblood_read.py" in sh and "jax_blood.py" in sh and "NUMPY_MADVISE_HUGEPAGE=0" in sh
    runs = [ln for ln in sh.splitlines() if "$PY" in ln and ".py" in ln] + [tests]
    assert all('> "$S/$NAME.out" 2> "$S/$NAME.err"' in r for r in runs)  # no log is printed
    assert "rm " not in sh and "rmdir" not in sh and "${1:?" not in sh and "${2:?" not in sh
    assert subprocess.run(["bash", "-n", str(LAUNCHER)]).returncode == 0
