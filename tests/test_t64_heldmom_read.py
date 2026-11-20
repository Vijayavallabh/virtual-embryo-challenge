"""Try T64 (experiments/harness/TRIES.md): experiments/t1_m1/heldmom_read.py on synthetic stages and hand-made task
files only. No released file, statistics file, board file or score is read; the committed chamber config is read for its
labels and the committed T62 task folder for its file names. It reuses row T62's fixtures."""

from __future__ import annotations

import importlib.util
import inspect
import json
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy import sparse
from test_g_mom_verify import _small_inputs
from test_t1_scope import GENES, HELD
from test_t1_state_transport import _gmom_ctx

from vec import harness, t1_scope
from vec.boards import Board
from vec.core import counts_of
from vec.t1_scope import build_submission, classify, load_spec, reweight

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "experiments" / "t1_m1" / "heldmom_read.py"
LAUNCHER = ROOT / "experiments" / "t1_m1" / "t64_heldmom.sh"
BOARD = Board(key="pseudo:t64", task="T1", genes=list(GENES), needs_coords=False, min_cells=100, max_cells=600)
SHARES = {"cm": 0.6, "endocardium": 0.25, "blood": 0.15}
DIRECTED = ("CM", "Endocardium", "Endothelium")
HELD_M1 = sorted({*HELD, GENES[3], GENES[7]})
WITHIN = {"surrounding": ["Foregut", "Epithelium"]}


@pytest.fixture
def hm(monkeypatch, tmp_path):
    spec = importlib.util.spec_from_file_location("heldmom_read_module", SCRIPT)
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


def _without_blood(stage, keep=0):
    """The stage with all but `keep` of its Blood cells removed (the blood category then has fewer than 20 own
    cells)."""
    lab = stage.obs["celltype"].astype(str).to_numpy()
    blood = np.flatnonzero(lab == "Blood")
    rows = np.sort(np.concatenate([np.flatnonzero(lab != "Blood"), blood[:keep]]))
    return stage[rows].copy()


def _cut_draw_labels(ctx):
    """The labels of the draw build_submission makes on a harness cut, in row order (the same reweight call)."""
    ref = harness._shallow(ctx.inputs[ctx.task.reference_role])
    src = ctx.task.role_stems[ctx.task.reference_role]
    classify(ref, source=src)
    shares = harness.t1_chamber_shares(ref, src)
    pop, _ = reweight(ref, shares, harness.count_for(ctx.board, ref.n_obs), ctx.seed,
                      fallback=load_spec().fallbacks(src))
    return pop.obs["celltype"].astype(str).to_numpy()


def _lib(hm, e95, e85, seed=4, dirs=None, fb=None):
    factors = hm.br.gmom_fits(e95, "E9.5_RNA", e85, "E8.5_RNA")["file"]
    dirs = _directions(e95.n_vars) if dirs is None else dirs
    fb = _fb() if fb is None else fb
    arms, rec = hm.library_arms(e85, e95, SHARES, factors, HELD, HELD_M1, seed, WITHIN, fb, dirs, 0.2, n_cells=300)
    base0, _ = reweight(e95, SHARES, 300, seed, within=WITHIN, fallback=fb["E9.5"])
    return arms, rec, base0, factors, dirs


# ---------------------------------------------------------------------------------------------------------------------
# (a) only cells with held labels change, on both routes


def test_on_the_cut_route_only_held_cells_change(hm):
    ref, start, _ = _small_inputs()
    ctx = _gmom_ctx(ref, start, seed=1000002)
    fits = hm.br.gmom_fits(ref, "E9.5_RNA", start, "E8.5_RNA")
    arms, rec = hm.cut_arms(ctx, fits["file"], _directions(ref.n_vars), 0.2, HELD_M1)
    assert tuple(arms) == hm.ARMS
    labels = _cut_draw_labels(ctx)
    held = np.isin(labels, rec["held_set"])
    differs = (_dense(arms["m1rc_gmom"].X) != _dense(arms["m1rc_gmom_heldmom"].X)).any(axis=1)
    assert rec["held_set"] == ["Blood"] and rec["uncovered_set"] == ["Blood"] and rec["non_held_equal"] is True
    assert differs.any() and not differs[~held].any()
    assert int(differs.sum()) == rec["held_cells_changed"] and rec["held_labels_drawn"] == {
        "Blood": int((labels == "Blood").sum())}
    assert any(r["type"] == "blood" and r["status"] == "shifted" for r in rec["shift"])


def test_on_the_library_route_only_held_cells_change(hm):
    e95, e85 = _classified()
    arms, rec, base0, _, _ = _lib(hm, e95, e85)
    labels = base0.obs["celltype"].astype(str).to_numpy()
    held = labels == "Blood"
    a, b = (_dense(counts_of(arms[x])) for x in hm.ARMS)
    differs = (a != b).any(axis=1)
    assert differs.any() and not differs[~held].any() and rec["non_held_equal"] is True
    assert rec["held_set"] == rec["uncovered_set"] == ["Blood"]
    assert int(differs.sum()) == rec["held_cells_changed"]


def test_a_cell_outside_the_held_set_that_differs_stops_the_run(hm):
    e95, _ = _classified()
    pop, _ = reweight(e95, SHARES, 300, 4, within=WITHIN)
    X = _dense(pop.X).astype(np.float32)
    Y = X.copy()
    Y[np.flatnonzero(pop.obs["celltype"].astype(str).to_numpy() == "V-CM")[0], 0] += 1.0
    groups = hm.feng_dirs.label_groups(sorted(set(pop.obs["celltype"].astype(str))))
    dirs = _directions(e95.n_vars)
    assert hm.held_record(pop, groups, dirs, X, X.copy(), None)["non_held_equal"] is True
    with pytest.raises(SystemExit, match="outside the held set"):
        hm.held_record(pop, groups, dirs, X, Y, None)


# ---------------------------------------------------------------------------------------------------------------------
# (b) the held cells' values by hand


def _hand_momentum(C, cats, S, s_cats, R, r_cats, held_genes, seed):
    """B0's damped per-category momentum on counts C (cells x genes) by hand: per category, mean CP10k over its own
    cells at each stage, f = ((m_ref + 0.5) / (m_start + 0.5)) ** 0.5, f = 1 on held genes, the keep_depth constant,
    then from np.random.default_rng(seed) a Binomial over every entry with f < 1 in CSR order, then a Poisson top-up
    over every entry with f > 1; a cell left empty gets its counts back."""
    def mean_cp(M):
        M = np.asarray(M, dtype=np.float64)
        return (M / M.sum(axis=1, keepdims=True) * 1e4).mean(axis=0)

    F = {}
    for c in sorted(set(cats)):
        f = ((mean_cp(R[r_cats == c]) + 0.5) / (mean_cp(S[s_cats == c]) + 0.5)) ** 0.5
        f[held_genes] = 1.0
        F[c] = f
    M = sparse.csr_matrix(np.asarray(C, dtype=np.int64))
    M.sort_indices()
    row = np.repeat(np.arange(M.shape[0]), np.diff(M.indptr))
    c = M.data.astype(np.float64)
    f = np.array([F[cats[i]][j] for i, j in zip(row, M.indices, strict=True)])
    free = ~held_genes[M.indices]
    num = np.bincount(row, weights=c * free, minlength=M.shape[0])
    den = np.bincount(row, weights=c * f * free, minlength=M.shape[0])
    scale = np.divide(num, den, out=np.ones(M.shape[0]), where=den > 0)
    f = np.where(free, f * scale[row], 1.0)
    rng = np.random.default_rng(seed)
    new = M.data.copy()
    down, up = f < 1.0, f > 1.0
    new[down] = rng.binomial(M.data[down], f[down])
    new[up] = M.data[up] + rng.poisson((f[up] - 1.0) * M.data[up])
    dead = np.bincount(row, weights=new, minlength=M.shape[0]) == 0
    new[dead[row]] = M.data[dead[row]]
    return _dense(sparse.csr_matrix((new, M.indices, M.indptr), shape=M.shape))


def test_the_held_cells_values_are_the_hand_computation(hm):
    e95, e85 = _classified()
    arms, rec, base0, _, _ = _lib(hm, e95, e85, seed=4)
    labels = base0.obs["celltype"].astype(str).to_numpy()
    idx = np.flatnonzero(labels == "Blood")
    held_genes = np.isin(np.asarray(e95.var_names), HELD_M1)
    want = _hand_momentum(_dense(counts_of(base0))[idx], base0.obs["category"].astype(str).to_numpy()[idx],
                          _dense(counts_of(e85)), e85.obs["category"].astype(str).to_numpy(),
                          _dense(counts_of(e95)), e95.obs["category"].astype(str).to_numpy(), held_genes, [4, 1])
    got = _dense(counts_of(arms["m1rc_gmom_heldmom"]))[idx]
    assert np.array_equal(got, want)
    drawn = _dense(counts_of(base0))[idx]
    assert not np.array_equal(got, drawn)  # the momentum moved them
    assert np.array_equal(got[:, held_genes], drawn[:, held_genes])  # a held gene stays as drawn


def test_a_category_with_too_few_own_cells_keeps_its_counts(hm):
    e95, e85 = _classified()
    arms, rec, base0, _, _ = _lib(hm, e95, _without_blood(e85, keep=10))
    idx = base0.obs["celltype"].astype(str).to_numpy() == "Blood"
    assert np.array_equal(_dense(counts_of(arms["m1rc_gmom_heldmom"]))[idx], _dense(counts_of(base0))[idx])
    assert [r["status"] for r in rec["shift"] if r["type"] == "blood"] == ["skipped: fewer than 20 cells at a stage"]


def test_a_category_drawn_from_stand_ins_keeps_its_counts(hm):
    e95, e85 = _classified()
    fb = _fb()
    fb["E8.5"] = {**fb["E8.5"], "blood": ["Foregut"]}  # blood has no own E8.5 cells here: stand-ins stand for it
    arms, rec, base0, _, _ = _lib(hm, e95, _without_blood(e85), fb=fb)
    idx = base0.obs["celltype"].astype(str).to_numpy() == "Blood"
    assert np.array_equal(_dense(counts_of(arms["m1rc_gmom_heldmom"]))[idx], _dense(counts_of(base0))[idx])
    row = [r for r in rec["shift"] if r["type"] == "blood"]
    assert row and row[0]["status"].startswith("skipped: stand-in")


# ---------------------------------------------------------------------------------------------------------------------
# (c) the library route is build_submission's path; (d) the base is T62's candidate as called


def test_the_library_route_before_g_mom_is_build_submissions_path(hm):
    e95, e85 = _classified()
    fb, dirs = _fb(), _directions(e95.n_vars)
    pop, rep, groups, _, _ = hm.library_population(e85, e95, SHARES, HELD_M1, 4, WITHIN, fb, dirs, 0.2, n_cells=300)
    seen = []
    _, report = build_submission(e95, BOARD, SHARES, 300, 4, start=e85, damping=0.5, ref_source="E9.5_RNA",
                                 start_source="E8.5_RNA", within=WITHIN, fallback=fb["E9.5"], held_flat=HELD_M1,
                                 directions=dirs, direction_fraction=0.2, direction_groups=groups, min_cells=20,
                                 eps=0.5, post_cells=lambda p: seen.append(p) or p)
    assert list(pop.obs_names) == list(seen[0].obs_names)
    assert _dense(counts_of(pop)).tobytes() == _dense(counts_of(seen[0])).tobytes()
    assert _dense(pop.X).tobytes() == _dense(seen[0].X).tobytes()
    assert report["shift"]["rows"] == rep.as_dict()["rows"]


def test_the_base_is_t62s_candidate_as_called_on_both_routes(hm):
    ref, start, _ = _small_inputs()
    ctx = _gmom_ctx(ref, start, seed=1000003)
    fits = hm.br.gmom_fits(ref, "E9.5_RNA", start, "E8.5_RNA")
    dirs = _directions(ref.n_vars)
    arms, _ = hm.cut_arms(ctx, fits["file"], dirs, 0.2, HELD_M1)
    t62 = hm.rc.cut_pair(ctx, fits["file"], dirs, 0.2, HELD_M1)
    assert arms["m1rc_gmom"].X.tobytes() == t62["m1rc_gmom"].X.tobytes()
    e95, e85 = _classified()
    arms, _, _, factors, dirs = _lib(hm, e95, e85, seed=5)
    t62 = hm.rc.library_pair(e95, SHARES, factors, HELD, HELD_M1, 5, WITHIN, _fb(), dirs, 0.2, n_cells=300)
    assert _dense(arms["m1rc_gmom"].X).tobytes() == _dense(t62["m1rc_gmom"].X).tobytes()


def test_the_cut_arm_is_build_submission_with_start_and_damping_added(hm):
    ref, start, _ = _small_inputs()
    ctx = _gmom_ctx(ref, start, seed=1000004)
    fits = hm.br.gmom_fits(ref, "E9.5_RNA", start, "E8.5_RNA")
    dirs = _directions(ref.n_vars)
    arms, _ = hm.cut_arms(ctx, fits["file"], dirs, 0.2, HELD_M1)
    cref = harness._shallow(ref)
    classify(cref, source="E9.5_RNA")
    groups, covered, with_direction = hm.rc.direction_kit(cref.obs["celltype"].astype(str).unique(), dirs)
    sub, _ = build_submission(cref, ctx.board, harness.t1_chamber_shares(cref, "E9.5_RNA"),
                              harness.count_for(ctx.board, cref.n_obs), ctx.seed, ref_source="E9.5_RNA",
                              held_flat=HELD_M1, directions=dirs, direction_fraction=0.2, direction_groups=groups,
                              start=harness._shallow(start), start_source="E8.5_RNA", damping=0.5,
                              post_cells=hm.rc.rc_hook(fits["file"], ctx.seed, HELD, with_direction, covered))
    want = harness.finish(ctx, cref, hm.rc.dense(sub.X), None)
    assert arms["m1rc_gmom_heldmom"].X.tobytes() == want.X.tobytes()


# ---------------------------------------------------------------------------------------------------------------------
# (e) the held set


E95_CHAMBER = ("AVC-CM", "IFT-CM", "OFT/RV-CM", "SV-CM", "V-CM", "Endocardium", "BEC", "Proepicardium", "Blood")
NOT_E95 = ("LV-CM", "RV-CM", "Endothelium",  # E8.5 only
           "Epicardium", "Endocardium-1", "Endocardium-2", "EndoMT Endocardium", "Early VIC", "Macrophage", "OFT-CM")
RECOMPUTED = ("CM", "Endocardium", "Endothelium", "Epicardium")  # the four groups with a direction in N54's record


def test_the_held_set_is_the_uncovered_set_and_on_e95s_chamber_labels_it_is_blood(hm):
    import yaml

    cfg = yaml.safe_load((ROOT / "configs" / "t1_chamber_classes.yaml").read_text())
    chamber = {c for c, v in cfg["categories"].items() if v["scope"] == "chamber"}
    labels = {str(lab) for lab, cat in cfg["labels"].items() if cat in chamber}
    assert set(E95_CHAMBER) <= labels and labels == set(E95_CHAMBER) | set(NOT_E95)
    dirs = {g: np.zeros(3) for g in RECOMPUTED}
    for labs in (E95_CHAMBER, sorted(labels), sorted(cfg["labels"]), ("V-CM", "Blood", "Mystery cells")):
        groups = hm.feng_dirs.label_groups(sorted(labs))
        assert hm.held_set(labs, dirs) == hm.uncovered_set(labs, groups, dirs)
    assert hm.held_set(E95_CHAMBER, dirs) == {"Blood"}
    assert hm.held_set(sorted(labels), dirs) == {"Blood", "Macrophage", "Early VIC"}  # D's held cells at the test
    assert "Mystery cells" in hm.held_set(("V-CM", "Mystery cells"), dirs)  # a label of no group is held


# ---------------------------------------------------------------------------------------------------------------------
# (f) every held category skipped: one file


def test_with_every_held_category_skipped_the_two_arms_are_one_file(hm):
    e95, e85 = _classified()
    arms, rec, _, _, _ = _lib(hm, e95, _without_blood(e85, keep=5))
    assert hm.prediction_sha(arms["m1rc_gmom"].X) == hm.prediction_sha(arms["m1rc_gmom_heldmom"].X)
    assert rec["held_cells_changed"] == 0
    ref, start, _ = _small_inputs()
    start = _without_blood(start, keep=5)
    ctx = _gmom_ctx(ref, start, seed=1000005)
    fits = hm.br.gmom_fits(ref, "E9.5_RNA", start, "E8.5_RNA")
    arms, rec = hm.cut_arms(ctx, fits["file"], _directions(ref.n_vars), 0.2, HELD_M1)
    assert hm.prediction_sha(arms["m1rc_gmom"].X) == hm.prediction_sha(arms["m1rc_gmom_heldmom"].X)
    assert rec["held_cells_changed"] == 0 and rec["held_labels_drawn"]["Blood"] > 0


# ---------------------------------------------------------------------------------------------------------------------
# (g) the momentum's stream and fit inputs


def test_the_momentums_stream_and_fit_inputs(hm, monkeypatch):
    real, calls = t1_scope.momentum_shift, []

    def spy(start, ref, target_like, damping, held_flat=(), **kw):
        calls.append({"start": start, "ref": ref, "damping": damping, "seed": kw.get("seed"),
                      "held": list(held_flat)})
        return real(start, ref, target_like, damping, held_flat, **kw)

    monkeypatch.setattr(t1_scope, "momentum_shift", spy)
    e95, e85 = _classified()
    _lib(hm, e95, e85, seed=6)
    assert len(calls) == 1 and calls[0]["start"] is e85 and calls[0]["ref"] is e95
    assert calls[0]["seed"] == [6, 1] and calls[0]["damping"] == 0.5 and calls[0]["held"] == HELD_M1
    calls.clear()
    ref, start, target = _small_inputs()
    ctx = _gmom_ctx(ref, start, seed=1000006)
    fits = hm.br.gmom_fits(ref, "E9.5_RNA", start, "E8.5_RNA")
    hm.cut_arms(ctx, fits["file"], _directions(ref.n_vars), 0.2, HELD_M1)
    assert len(calls) == 1 and calls[0]["seed"] == [1000006, 1] and calls[0]["damping"] == 0.5
    assert list(calls[0]["start"].obs_names) == list(start.obs_names)
    assert list(calls[0]["ref"].obs_names) == list(ref.obs_names)
    assert list(calls[0]["ref"].obs_names) != list(target.obs_names)
    for fn in (hm.cut_arms, hm.library_arms, hm.library_population):  # no truth reaches an arm
        assert not {"truth", "target", "e105"} & set(inspect.signature(fn).parameters)
    assert hm.MOMENTUM_STREAM == 1 and hm.DAMPING == 0.5 == harness.B0_DAMPING


# ---------------------------------------------------------------------------------------------------------------------
# The reading on hand-made task files


def _arm(score):
    return {"score": score, "groups": {g: score / 100 for g in ("de_recovery", "direction", "distribution",
                                                               "covariance")}}


def _relabel(hm, **change):
    return {"sha256": hm.rc.STATS_SHA256, "stats_sha256_hashed": hm.rc.STATS_SHA256, "yaml_sha256": hm.rc.YAML_SHA256,
            "route": "recomputed", "record_route": "recomputed", "authors_labels_read": False, **change}


def _task(hm, gain, lib=False, per_seed=None, base=50.0, commit="c" * 40, clean=True, seeds=range(10), inverted=False,
          relabel=None, idle=False, **top):
    rows = []
    for s in seeds:
        g = per_seed[s] if per_seed is not None else gain
        arms = {"m1rc_gmom": _arm(base + s / 100), "m1rc_gmom_heldmom": _arm(base + s / 100 + g)}
        sha = {"m1rc_gmom": "a" * 64, "m1rc_gmom_heldmom": "a" * 64 if idle else "b" * 64}
        held = {"held_set": ["Blood"], "uncovered_set": ["Blood"], "held_labels_drawn": {"Blood": 30},
                "held_cells": 30, "held_cells_changed": 0 if idle else 30, "non_held_equal": True,
                "shift": [{"type": "blood", "status": "shifted", "n_start": 266, "n_ref": 213}]}
        row = {"seed": s, "arms": arms, "sha256": sha, "held": held}
        if lib:
            arms["copy_last"], sha["copy_last"] = _arm(40.0), "c" * 64
            row["inverted_anchors"] = ["x"] if inverted and s == 3 else []
        else:
            row["floor"] = 50.0
        rows.append(row)
    out = {"rows": rows, "stamp": {"commit": commit, "tree_clean": clean}, "fraction": 0.2, "damping": 0.5,
           "momentum_stream": 1, "relabel": relabel if relabel is not None else _relabel(hm),
           "shift_rows": rows[0]["held"]["shift"] if rows else None}
    out.update(top)
    return out


def _committed(base=50.0):
    return {k: {s: base + s / 100 for s in range(10)} for k in ("chamber", "whole", "1", "2", "3")}


def _loaded(hm, chamber=1.0, whole=0.5, libs=(0.0, 0.0, 0.0), **kw):
    out = {"chamber": _task(hm, chamber), "whole": _task(hm, whole)}
    out.update({lib: _task(hm, g, lib=True) for lib, g in zip(("1", "2", "3"), libs, strict=True)})
    out.update(kw)
    return out


OK = {"reasons": []}


def _p(mean, ahead=10):
    return {"mean": mean, "ahead": ahead}


# (h) the statement


def test_ahead(hm):
    r = hm.t64_reading(_loaded(hm), _committed(), OK)
    assert r["outcome"] == "D is AHEAD" and r["ahead"] and r["follows"] == hm.FOLLOWS_AHEAD
    assert r["j6_range"] == pytest.approx([0.0, 1.0]) and r["valid_libraries"] == ["1", "2", "3"]
    lines = hm.outcome_lines(r)
    assert lines[0] == "Outcome: D is AHEAD" and all(b in lines for b in hm.BESIDE)
    assert "QD" in lines[-1] and "keep R1" in lines[-1]


@pytest.mark.parametrize("tasks,valid,part", [
    ({"chamber": _p(0.3), "whole": _p(1.0)}, ["1", "2"], "chamber cut's mean"),
    ({"chamber": _p(1.0, 7), "whole": _p(1.0)}, ["1", "2"], "fewer than 8 of 10"),
    ({"chamber": _p(1.0), "whole": _p(0.0)}, ["1", "2"], "whole cut's mean"),
    ({"chamber": _p(1.0), "whole": _p(-0.01)}, ["1", "2"], "whole cut's mean"),
    ({"chamber": _p(1.0), "whole": _p(1.0), "1": _p(0.0), "2": _p(-0.31)}, ["1", "2"], "below -0.3 (library 2)"),
    ({"chamber": _p(1.0), "whole": _p(1.0), "3": _p(1.0)}, ["3"], "fewer than 2 library tasks"),
])
def test_not_ahead(hm, tasks, valid, part):
    tasks = {"1": _p(0.0), "2": _p(0.0), "3": _p(0.0), **tasks}
    ahead, failed = hm.statement(tasks, valid)
    assert not ahead and any(part in f for f in failed)


def test_the_boundaries_that_pass(hm):
    tasks = {"chamber": _p(0.3000001, 8), "whole": _p(1e-9), "1": _p(-0.3), "2": _p(0.0), "3": _p(-5.0)}
    assert hm.statement(tasks, ["1", "2"]) == (True, [])  # library 3 invalid: its -5 is ignored


def test_an_invalid_library_is_ignored_and_two_invalid_libraries_stop_it(hm):
    loaded = _loaded(hm, libs=(0.0, 0.0, -5.0))
    loaded["3"] = _task(hm, -5.0, lib=True, inverted=True)
    r = hm.t64_reading(loaded, _committed(), OK)
    assert r["outcome"] == "D is AHEAD" and r["valid_libraries"] == ["1", "2"]
    loaded["2"] = _task(hm, 0.0, lib=True, inverted=True)
    r = hm.t64_reading(loaded, _committed(), OK)
    assert r["outcome"] == "D is NOT AHEAD" and any("fewer than 2 library" in f for f in r["failed"])
    assert all(b in hm.outcome_lines(r) for b in hm.BESIDE) and hm.outcome_lines(r)[-1] == hm.FOLLOWS_END


def test_a_valid_library_below_the_floor_stops_it_end_to_end(hm):
    r = hm.t64_reading(_loaded(hm, libs=(0.0, -0.31, 0.0)), _committed(), OK)
    assert r["outcome"] == "D is NOT AHEAD" and r["follows"] == hm.FOLLOWS_END


def test_a_negative_chamber_mean_gives_the_j6_range_mean_to_zero(hm):
    r = hm.t64_reading(_loaded(hm, chamber=-0.5), _committed(), OK)
    assert r["outcome"] == "D is NOT AHEAD" and r["j6_range"] == pytest.approx([-0.5, 0.0])


def test_the_idle_case_is_named_and_read_as_it_is(hm):
    loaded = _loaded(hm, libs=(0.0, 0.0, 0.0))
    loaded["2"] = _task(hm, 0.0, lib=True, idle=True)
    r = hm.t64_reading(loaded, _committed(), OK)
    assert r["idle_tasks"] == ["2"] and r["outcome"] == "D is AHEAD"
    assert "the arm moved no cell on 2: nothing of D is priced there" in hm.outcome_lines(r)


def test_a_held_label_other_than_blood_is_named(hm):
    loaded = _loaded(hm)
    for x in loaded["3"]["rows"]:
        x["held"]["held_labels_drawn"]["Macrophage"] = 3
    r = hm.t64_reading(loaded, _committed(), OK)
    assert r["other_held_labels"] == {"3": ["Macrophage"]}
    assert "held labels other than Blood drawn on 3: Macrophage" in hm.outcome_lines(r)


def test_no_line_or_text_carries_advice(hm):
    for loaded in (_loaded(hm), _loaded(hm, chamber=0.0)):
        r = hm.t64_reading(loaded, _committed(), OK)
        res = {"label": hm.LABEL, "reading": r, "tries": 5, "stamp": {"commit": "e", "tree_clean": True, "utc": "x"}}
        assert "advice" not in hm.t64_text(res).lower()
    assert "advice" not in SCRIPT.read_text().lower()


# (i) every void reason


@pytest.mark.parametrize("case", ["missing_task", "missing_seed", "missing_arm", "missing_sha", "missing_copy_last",
                                  "two_commits", "unclean", "gate", "stats", "yaml", "route", "authors", "checks",
                                  "source", "fraction", "damping", "stream", "non_held", "held_set"])
def test_every_void_reason_fires(hm, case):
    loaded, committed, checks, source = _loaded(hm), _committed(), OK, []
    if case == "missing_task":
        loaded["1"] = None
    elif case == "missing_seed":
        loaded["whole"] = _task(hm, 0.5, seeds=range(9))
    elif case == "missing_arm":
        del loaded["chamber"]["rows"][2]["arms"]["m1rc_gmom_heldmom"]
    elif case == "missing_sha":
        del loaded["whole"]["rows"][5]["sha256"]["m1rc_gmom"]
    elif case == "missing_copy_last":
        del loaded["3"]["rows"][0]["arms"]["copy_last"]
    elif case == "two_commits":
        loaded["1"] = _task(hm, 0.0, lib=True, commit="d" * 40)
    elif case == "unclean":
        loaded["chamber"] = _task(hm, 1.0, clean=False)
    elif case == "gate":
        committed["2"][4] += 0.0201
    elif case == "stats":
        loaded["chamber"] = _task(hm, 1.0, relabel=_relabel(hm, stats_sha256_hashed="0" * 64))
    elif case == "yaml":
        loaded["whole"] = _task(hm, 0.5, relabel=_relabel(hm, yaml_sha256="0" * 64))
    elif case == "route":
        loaded["1"] = _task(hm, 0.0, lib=True, relabel=_relabel(hm, route="none"))
    elif case == "authors":
        loaded["2"] = _task(hm, 0.0, lib=True, relabel=_relabel(hm, authors_labels_read=True))
    elif case == "checks":
        checks = {"reasons": ["the settings file's sha256 is x"]}
    elif case == "source":
        source = ["a feng_directions call without the recomputed groups"]
    elif case == "fraction":
        loaded["chamber"] = _task(hm, 1.0, fraction=0.4)
    elif case == "damping":
        loaded["whole"] = _task(hm, 0.5, damping=1.0)
    elif case == "stream":
        loaded["3"] = _task(hm, 0.0, lib=True, momentum_stream=2)
    elif case == "non_held":
        loaded["2"]["rows"][6]["held"]["non_held_equal"] = False
    elif case == "held_set":
        loaded["1"]["rows"][1]["held"]["uncovered_set"] = ["Blood", "Foregut"]
    r = hm.t64_reading(loaded, committed, checks, source)
    assert r["outcome"] == "VOID" and r["reasons"] and r["follows"] == hm.FOLLOWS_END
    lines = hm.outcome_lines(r)
    assert lines[0].startswith("Outcome: VOID") and lines[-1] == hm.FOLLOWS_END
    assert not any(b in lines for b in hm.BESIDE)
    if case == "gate":
        assert "the gate failed" in lines[0] and "0.0201" in lines[0]


def test_the_gate_passes_at_a_difference_of_0_02(hm):
    loaded, committed = _loaded(hm), _committed()
    committed["2"][4] = 0.0
    loaded["2"]["rows"][4]["arms"]["m1rc_gmom"]["score"] = 0.02
    loaded["2"]["rows"][4]["arms"]["m1rc_gmom_heldmom"]["score"] = 0.02
    assert hm.t64_reading(loaded, committed, OK)["outcome"] != "VOID"
    loaded["2"]["rows"][4]["arms"]["m1rc_gmom"]["score"] = 0.0201
    assert hm.t64_reading(loaded, committed, OK)["outcome"] == "VOID"
    assert hm.REPRO_TOL == 0.02 == hm.br.REPRO_TOL


def test_committed_t62_reads_m1rc_gmom_per_seed(hm, tmp_path):
    (tmp_path / "lib1.json").write_text(json.dumps({"rows": [{"seed": s, "arms": {"m1rc_gmom": {"score": 60.0 + s}}}
                                                            for s in range(10)]}))
    assert hm.committed_t62("lib1.json", tmp_path) == {s: 60.0 + s for s in range(10)}
    assert hm.T62_TASKS == SCRIPT.parent / "results" / "t62_recomputed" / "tasks"
    assert all((hm.T62_TASKS / f).exists() for f in hm.task_files().values())  # the gate's committed files


# (j) refusals, counts, guards


def test_the_counts_and_constants(hm):
    assert hm.TRIES_COUNT == 5 and hm.TRY_RANGE == (413, 418) and hm.COUNTED == ("m1rc_gmom_heldmom",)
    assert hm.ARMS == ("m1rc_gmom", "m1rc_gmom_heldmom") and hm.LIB_ARMS[0] == "copy_last"
    assert hm.BASE == hm.rc.CANDIDATE and hm.GATED == ("m1rc_gmom",)
    assert hm.AHEAD_BAR == 0.3 and hm.MIN_AHEAD == 8 and hm.WHOLE_ABOVE == 0.0 and hm.LIB_FLOOR == -0.3
    assert hm.MIN_LIBRARIES == 2 and hm.PRACTICE_FRACTION == 0.2
    assert not hasattr(hm, "WHOLE_FLOOR") and "WHOLE_FLOOR" not in SCRIPT.read_text()


def test_the_run_and_the_read_are_refused_without_the_row(hm, tmp_path, monkeypatch):
    tries = tmp_path / "TRIES.md"
    tries.write_text("| T63 (one row) | T1 | x |\n")
    monkeypatch.setattr(hm, "TRIES", tries)
    monkeypatch.setattr(hm.rc, "relabel_checks", lambda *a, **k: pytest.fail("nothing before the row's refusal"))
    monkeypatch.setattr(hm.gv, "tree_state", lambda: pytest.fail("nothing before the row's refusal"))
    with pytest.raises(SystemExit, match="no row"):
        hm.refuse(tmp_path / "out", "x.json")
    with pytest.raises(SystemExit, match="no row"):
        hm.read(tmp_path / "out")


def test_the_run_is_refused_on_another_damping(hm, tmp_path, monkeypatch):
    tries = tmp_path / "TRIES.md"
    tries.write_text("| T64 (one row) | T1 | x |\n")
    monkeypatch.setattr(hm, "TRIES", tries)
    monkeypatch.setattr(hm.gv, "tree_state", lambda: {"commit": "e" * 40, "clean": True})
    monkeypatch.setattr(hm.rc, "relabel_checks", lambda *a, **k: {"passes": True, "reasons": []})
    assert hm.refuse(tmp_path / "out", "x.json")[0]["clean"] is True
    monkeypatch.setattr(hm, "DAMPING", 1.0)
    with pytest.raises(SystemExit, match="damping"):
        hm.refuse(tmp_path / "out", "x.json")
    monkeypatch.setattr(hm.gv, "tree_state", lambda: {"commit": "e" * 40, "clean": False})
    with pytest.raises(SystemExit, match="not clean"):
        hm.refuse(tmp_path / "out", "x.json")


def test_the_source_guard(hm, tmp_path):
    assert hm.source_guard() == [] and hm.rc.source_guard(SCRIPT) == []
    for text in ('x = fd.feng_directions(list(g))\n', 'x = fd.feng_directions(g, groups="none")\n',
                 "import perm_read\n", "from double_step_read import BANNED_IMPORTS\n",
                 "import n5" + "4_authors_agreement as a\n", "m = import" + "lib.import_module('x')\n"):
        bad = tmp_path / "bad.py"
        bad.write_text(text)
        assert hm.source_guard(bad), text


def test_the_read_prints_no_score(hm, tmp_path, monkeypatch, capsys):
    out = tmp_path / "results" / "t64"
    (out / "tasks").mkdir(parents=True)
    loaded = _loaded(hm, 1.234567, 0.5, (0.4, 0.4, 0.4))
    for k, f in hm.task_files().items():
        (out / "tasks" / f).write_text(json.dumps(loaded[k]))
    committed = _committed()
    tries = tmp_path / "TRIES.md"
    tries.write_text("| T64 (one row) | T1 | x |\n")
    monkeypatch.setattr(hm, "TRIES", tries)
    monkeypatch.setattr(hm, "committed_t62", lambda f: committed[{v: k for k, v in hm.task_files().items()}[f]])
    monkeypatch.setattr(hm.rc, "relabel_checks", lambda **k: {"reasons": [], "passes": True})
    monkeypatch.setattr(hm.gv, "stamp", lambda started, files: {"commit": "e" * 40, "tree_clean": True, "utc": "x"})
    monkeypatch.setattr(hm.gv, "tree_state", lambda: {"commit": "e" * 40, "clean": True})
    res = hm.read(out)
    printed = capsys.readouterr().out
    assert res["reading"]["outcome"] == "D is AHEAD" and res["tries"] == 5
    assert "1.23" not in printed and "51." not in printed and "50." not in printed and "Outcome" in printed
    assert all(b in printed for b in hm.BESIDE) and "advice" not in printed.lower()
    text = (out / "t64.txt").read_text()
    assert "+1.23" in text and "J6 range" in text and "QD, priced both ways" in text and hm.QD_FIXED[1] in text
    assert "momentum blood: shifted" in text and "T1 413 to 418" in text
    with pytest.raises(SystemExit):
        hm.read(out)  # no second read into the same folder


# ---------------------------------------------------------------------------------------------------------------------
# The launcher


def test_the_launcher_runs_the_tests_first_and_never_prints_a_log():
    sh = LAUNCHER.read_text()
    tests = sh.split("  tests)")[1].split(";;")[0]
    for f in ("test_t64_heldmom_read.py", "test_t62_recomputed_read.py", "test_t63_double_step_read.py",
              "test_t1_bets.py", "test_t1_scope.py", "test_core.py", "test_g_mom_verify.py",
              "test_t1_state_transport.py"):
        assert f"tests/{f}" in tests and (ROOT / "tests" / f).exists()
    for code in ("98", "96", "95", "93", "92", "91"):
        assert f"finish {code}" in sh
    assert "heldmom_read.py" in sh and "NUMPY_MADVISE_HUGEPAGE=0" in sh
    assert "OUT=t64_heldmom" in sh and '--out "$OUT"' in sh
    assert "feng_relabel_c/cd1_stats_relabel.npz" in sh and ".exit" in sh
    runs = [ln for ln in sh.splitlines() if "heldmom_read.py" in ln and "$PY" in ln] + [tests]
    assert len(runs) == 3 and all('> "$S/$NAME.out" 2> "$S/$NAME.err"' in r for r in runs)  # no log is printed
