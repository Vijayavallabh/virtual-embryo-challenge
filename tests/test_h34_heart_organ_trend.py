"""Row H34 (experiments/harness/h34_heart_organ_trend.py) and its launcher (experiments/harness/run_h34.sh): synthetic
arrays, hand-made tables and temporary files only. Nothing here loads a competition file, scores a real cell or reads
the repository's own ledger or a committed result; the launcher runs in throwaway git repositories under the test's
own folder, with a temporary home and a stand-in uv first on PATH."""

from __future__ import annotations

import ast
import copy
import hashlib
import importlib.util
import inspect
import json
import os
import re
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from vec import harness, prep, t1_scope, t2_build

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "experiments" / "harness" / "h34_heart_organ_trend.py"
LAUNCHER = REPO_ROOT / "experiments" / "harness" / "run_h34.sh"
COMMIT = "c" * 40
KEY1, KEY2 = "T2h/<=E9.5->E10.5", "T2h/<=E8.75->E10.5"
METRICS = ("de_score", "de_direction", "mmd_u", "variogram", "neighborhood_mmd", "d2_shape", "occupancy_dice",
           "scale_log_ratio")
GROUPS = ("expression_change", "cell_state", "shape_growth", "local_organisation")
LABELS = np.array(["V-CM", "BEC", "NCC", "Epi"])  # cm (the core), vascular endothelium, surrounding, epicardium


@pytest.fixture(scope="module")
def h34():
    spec = importlib.util.spec_from_file_location("h34_heart_organ_trend_module", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def _labelled(mod, a):
    merged = np.asarray(a.obs["merged"]).astype(str)
    a.obs["core"] = merged == "CM"
    kinds = {m: i for i, m in enumerate(mod.LINEAGES)}
    a.obs["celltype"] = LABELS[[kinds[m] for m in merged]]
    return a


@pytest.fixture(scope="module")
def synth():
    """test_b0_trend's synthetic T2 stages with T2 labels by lineage and the core (CM). The tiny lineage has enough
    cells at `start` and at the reference only: the trend mapping leaves it out and it takes B0's own factor, the
    damping's path through the strength."""
    spec = importlib.util.spec_from_file_location("b0_trend_tests_for_h34", REPO_ROOT / "tests" / "test_b0_trend.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    stages = {"T2/E9.5.h5ad": mod._stage(1600, seed=2, shift=0.0, tiny=160),
              "T2/E8.75.h5ad": mod._stage(1100, seed=3, shift=0.8, tiny=120),
              "T2/E8.5.h5ad": mod._stage(900, seed=5, shift=1.6, tiny=10),
              "T2/E8.25_late.h5ad": mod._stage(800, seed=9, shift=2.0, tiny=10),
              "T2/E10.5.h5ad": mod._stage(1400, seed=7, shift=-0.5, tiny=50)}
    for a in stages.values():
        _labelled(mod, a)
    return mod, stages


ROLES = {KEY1: {"ref": "T2/E9.5.h5ad", "start": "T2/E8.75.h5ad", "earlier": "T2/E8.5.h5ad"},
         KEY2: {"ref": "T2/E8.75.h5ad", "start": "T2/E8.5.h5ad", "earlier": "T2/E8.25_late.h5ad"}}
TARGET = "T2/E10.5.h5ad"


def _arrays(a):
    return (np.asarray(a.obsm["spatial_3D"], float), a.obs["celltype"].to_numpy(), np.asarray(a.obs["core"], bool))


def _geoms(h34, stages, files):
    sp = h34.h27.sp
    return [sp.geometry_from_arrays(harness.stem(f), day, _arrays(stages[f])[0], _arrays(stages[f])[2])
            for f, day in files]


@pytest.fixture(scope="module")
def dform(h34, synth):
    """D's form on the synthetic stages of the first transition: the pool, both ratios, D_f's size, a context maker."""
    mod, stages = synth
    inputs = {r: stages[f] for r, f in ROLES[KEY1].items()}
    whole = inputs["ref"]
    PC, labels, pcore = _arrays(whole)
    masks, _ = h34.h29b.organ_pools(PC, labels, pcore)
    ratios = {p: t2_build.organ_core_rms_ratio(PC, masks[p], pcore) for p in h34.h29b.POOLS}
    size = h34.h27.analogue_size_result(_geoms(h34, stages, [("T2/E8.5.h5ad", 8.5), ("T2/E8.75.h5ad", 8.75),
                                                              ("T2/E9.5.h5ad", 9.5)]), 10.5)
    ft = h34.h27.forecast_targets(size)
    ins = {**inputs, "ref": whole[masks["organ"]].copy()}

    def ctx(seed=4, view=None):
        return mod._ctx(ins, view=view, seed=harness.BUILD_SEED_OFFSET + seed)
    target = h34.h29b.d_targets(ft, {h34.h29b.CORE_FORM: float("nan")}, ratios)[h34.D_F]
    return SimpleNamespace(mod=mod, inputs=ins, whole=whole, masks=masks, ratios=ratios, ft=ft, ctx=ctx,
                           target=target)


# ------------------------------------------------------------------------------------------------ constants, imports
def test_the_fixed_constants_and_names(h34):
    assert h34.ROW == "H34" and h34.TASK_KEYS == (KEY1, KEY2) and tuple(h34.SEEDS) == tuple(range(10))
    assert h34.TRUTHS == ("own", "radial_0", "thin_0.5", "thin_0")
    assert h34.WORLDS == (("own", ("own",)), ("CLEAN", ("radial_0", "thin_0.5", "thin_0")))
    assert h34.STRENGTHS == (1.0, 1.5) and h34.CONTROL == "control" and h34.D_F == "organ:f"
    assert h34.COUNTED == ("trend_s1.0", "placebo_s1.0", "trend_s1.5", "placebo_s1.5")
    assert len(h34.COUNTED) * len(h34.TASK_KEYS) == 8  # eight tries; the control is a reference
    assert h34.BUILT == ("control", *h34.COUNTED) and h34.ALL == ("floor", *h34.BUILT)
    assert h34.PERM_SEED == 20261120 and h34.BAR == 0.3 and h34.MIN_SEEDS == 8
    assert h34.PRACTICE_BUDGET == {KEY1: 0.375, KEY2: 0.125} and h34.FILE == "OTS15" and h34.RULING_ROW == "R32"
    assert h34.GATE_EXIT == {"PASS": 0, "FAIL": 3} and h34.TOL == 1e-9 and h34.FACTOR_RTOL == 1e-12
    assert re.fullmatch(r"[0-9a-f]{64}", h34.H29B_READING_SHA256)
    assert h34.H29B_READING == REPO_ROOT / "experiments" / "harness" / "results" / "h29b_heart_organ" / "reading.json"
    assert h34.RULE_PAIRS == (("trend_s1.5", "trend_s1.0"), ("trend_s1.5", "control"))
    assert [w for w, _ in h34.G1_WORLDS] == ["own", "CLEAN", "truth own", "truth radial_0", "truth thin_0.5",
                                              "truth thin_0"]
    assert h34.G1_UNITS == (KEY1, KEY2, "fold E10.5") and set(h34.UNITS) == set(h34.G1_UNITS)
    assert all(u in h34.h29b.UNITS for u in h34.G1_UNITS)
    assert h34.RESULTS == REPO_ROOT / "experiments" / "harness" / "results" / "h34_heart_organ_trend"


def test_rows_h28s_and_h29bs_code_is_imported_not_copied(h34):
    h28, h29b = h34.h28, h34.h29b
    assert Path(h28.__file__).name == "h28_heart_trend_strength.py"
    assert Path(h29b.__file__).name == "h29b_heart_organ.py"
    assert h34.TASK_KEYS is h28.TASK_KEYS and h34.WORLDS is h28.WORLDS and h34.SEEDS is h28.SEEDS
    assert h34.UNITS is h28.UNITS and h34.BAR is h28.BAR and h34.pname is h28.pname
    assert h34.h27 is h29b.h27 and h34.h26 is h29b.h26 and h34.h25 is h29b.h25
    tree = ast.parse(SCRIPT.read_text())
    defined = {n.name for n in ast.walk(tree) if isinstance(n, ast.FunctionDef)}
    for theirs in ("trend_expression", "derangement", "lineage_log_means", "series", "passes", "size_to",
                   "organ_pools", "d_targets", "kept_rows", "truth_masks", "analogue_size_result", "_checked"):
        assert theirs not in defined, theirs
    src = SCRIPT.read_text()
    for used in ("h28.trend_expression(", "h28.series(", "h28.compare(", "h28.passes(", "h29b.organ_pools(",
                 "h29b.d_targets(", "h29b.size_to(", "h25.kept_rows(", "h29b._checked(", "h29b.series("):
        assert used in src, used


def test_vec_is_this_checkouts_own(h34):
    """The script puts its own checkout's src first, before any vec import, and refuses another checkout's vec."""
    body = ast.parse(SCRIPT.read_text()).body
    first_vec = next(i for i, n in enumerate(body) if isinstance(n, ast.ImportFrom) and (n.module or "")
                     .split(".")[0] == "vec")
    insert = next(i for i, n in enumerate(body) if isinstance(n, ast.Expr)
                  and ast.unparse(n).startswith("sys.path.insert(0, str(HERE / 'src'))"))
    assert insert < first_vec
    assert not [n for n in body[:insert] if isinstance(n, ast.Import | ast.ImportFrom)
                and any((getattr(n, "module", None) or a.name).split(".")[0] == "vec" for a in n.names)]
    assert h34.HERE == REPO_ROOT and h34.vec_home() == REPO_ROOT
    assert sys.path.index(str(REPO_ROOT / "src")) == 0 or Path(harness.__file__).is_relative_to(REPO_ROOT / "src")


def test_the_script_names_no_test_stage_and_reads_no_test_score():
    src = SCRIPT.read_text()
    for bad in ("results/p3", "/p3", "size_heart_E12", "submissions", "leaderboard", "registry", "build_test",
                "boards/", "12.5", "manual_uploads", "index.json", "test_extrap", "size_or_hold(", "headline",
                "halve.py"):
        assert bad not in src, bad
    # the one word of the test stage is in the row's fixed what-follows sentence: "NOT STOPPED uploads nothing"
    assert src.count("upload") == 1 and "NOT STOPPED uploads nothing" in src


# ------------------------------------------------------------------------------------------------ the arms (synthetic)
def test_the_control_is_row_h29bs_d_f_bit_for_bit(h34, dform):
    real = h34.h29b.d_targets(dform.ft, {h34.h29b.CORE_FORM: 3.0}, dform.ratios)
    assert real[h34.D_F] == dform.target  # D_f's size needs no target: the oracle control's field is not formed
    theirs, recs, draw = h34.h29b.build_pool_arms("organ", dform.ctx(), real)
    mine, rec, mdraw = h34.build_control(dform.ctx(), dform.target)
    np.testing.assert_array_equal(mine.X, theirs[h34.D_F].X)
    np.testing.assert_array_equal(mine.coords, theirs[h34.D_F].coords)
    assert rec == recs[h34.D_F] and mdraw == draw
    out, extra = h34.build_arms(dform.ctx(), dform.target)
    np.testing.assert_array_equal(out[h34.CONTROL].X, mine.X)
    np.testing.assert_array_equal(out[h34.CONTROL].coords, mine.coords)
    assert extra["records"][h34.CONTROL]["size"] == rec and extra["draw"] == draw
    assert extra["records"][h34.CONTROL]["trend"] == {"strength": None}


def test_every_arm_is_the_control_with_the_expression_step_before_the_size_step(h34, dform):
    out, extra = h34.build_arms(dform.ctx(), dform.target)
    assert list(out) == list(h34.BUILT) == list(extra["records"])
    c = extra["checks"]
    for k in ("s1_equals_harness_b0_trend", "same_draw", "pre_size_coords_equal", "post_size_coords_equal"):
        assert c[k] is True, k
    assert len(set(c["n_cells"].values())) == 1 and c["n_cells"]["control"] == dform.mod.BOARD.max_cells // 2
    assert set(c["changed_counts"]) == set(h34.COUNTED) and all(v > 0 for v in c["changed_counts"].values())
    for name, st in out.items():
        np.testing.assert_array_equal(st.coords, out["control"].coords)
        rec = extra["records"][name]["size"]
        assert abs(rec["log_rms_after_vs_target"]) <= 1e-9 and rec["target_rms_full"] == dform.target["target"]
    # each arm is row H28's step on D's context, then row H29b's size step and half: the step comes first
    shifted = set()
    for s in h34.STRENGTHS:
        for name, is_placebo in ((h34.arm(s), False), (h34.placebo(s), True)):
            st = h34.h28.trend_expression(dform.ctx(), s, placebo_arm=is_placebo)[0]
            want, _ = h34.size_and_half(dform.ctx(), st, dform.target, name)
            np.testing.assert_array_equal(out[name].X, want.X)
            tr = extra["records"][name]["trend"]
            assert tr["strength"] == s and tr["placebo"] is is_placebo
            assert tr["ref_cells"] == dform.inputs["ref"].n_obs < dform.whole.n_obs  # the reference is D's pool
            assert tr["budget_days"] == 0.375 and tr["effective_budget_days"] == s * 0.375
            assert not any(tr["fixed_points"].values()) and bool(tr["fixed_points"]) is is_placebo
            shifted.add(tuple(tr["shifted"]))
    assert len(shifted) == 1 and {"CM", "Endothelial", "Tiny"} == set(next(iter(shifted)))
    # strength 1 is the harness's own b0_trend on D's pool, sized and halved as the control
    b0, _ = h34.size_and_half(dform.ctx(), harness._b0_t2(dform.ctx(), trend=True), dform.target, "b0_trend")
    np.testing.assert_array_equal(out["trend_s1.0"].X, b0.X)
    for a, b in (("trend_s1.0", "control"), ("trend_s1.5", "trend_s1.0"), ("placebo_s1.5", "trend_s1.5")):
        assert not np.array_equal(out[a].X, out[b].X), (a, b)
    # a seed changes the draw; the same seed gives the same arms
    again, _ = h34.build_arms(dform.ctx(), dform.target)
    other, _ = h34.build_arms(dform.ctx(seed=5), dform.target)
    np.testing.assert_array_equal(again["trend_s1.5"].X, out["trend_s1.5"].X)
    assert not np.array_equal(other["control"].coords, out["control"].coords)


def test_with_the_trend_switched_off_the_arms_cells_and_counts_are_the_controls(h34, dform):
    out, extra = h34.build_arms(dform.ctx(), dform.target, strengths=(1.0, 0.0))
    off, ctrl = out[h34.arm(0.0)], out[h34.CONTROL]
    c = extra["checks"]
    # the counts: every drawn count as drawn, in the arm and in its placebo; strength 1 changes some
    assert c["changed_counts"][h34.arm(0.0)] == 0 == c["changed_counts"][h34.placebo(0.0)]
    assert c["changed_counts"][h34.arm(1.0)] > 0
    # the cells: the control's draw, positions (before and after the size step and half) and cell number
    assert c["same_draw"] is True and c["pre_size_coords_equal"] is True and c["post_size_coords_equal"] is True
    np.testing.assert_array_equal(off.coords, ctrl.coords)
    assert off.n == ctrl.n == c["n_cells"][h34.arm(0.0)]
    # the expression is the control's, re-rendered from the same counts
    np.testing.assert_allclose(off.X, ctrl.X, rtol=0, atol=1e-5)
    # and by the step itself: at strength 0 it hands back the drawn cells' counts
    _, counts, rows, _ = h34.h28.trend_expression(dform.ctx(), 0.0)
    pool, drawn_rows, _, _ = harness.draw_copy(dform.ctx())
    np.testing.assert_array_equal(rows, drawn_rows)
    assert (counts != h34.counts_of(pool[drawn_rows], list(dform.mod.BOARD.genes))).nnz == 0
    with pytest.raises(ValueError, match="strength 1 must be among"):
        h34.build_arms(dform.ctx(), dform.target, strengths=(1.5,))


# ------------------------------------------------------------------------------------------------ the gate's G2
def test_the_gate_reads_the_steps_arguments_and_builds_no_arm(h34, dform, monkeypatch):
    _, _, _, rec = h34.h28.trend_expression(dform.ctx(), 1.0)
    assert rec["fallback"] == ["Tiny"] and rec["n_trend"] == 2  # the mapping's two lineages; one through the damping
    own = t1_scope.momentum_shift
    pool, rows, _, _ = harness.draw_copy(dform.ctx())

    def boom(*a, **k):
        raise AssertionError("a count was touched or an arm was built")
    monkeypatch.setattr(t1_scope, "apply_factors", boom)  # the step is stopped before any count is touched
    monkeypatch.setattr(harness, "finish", boom)  # and it never reaches a stage: no arm is built
    lo, hi = h34.step_arguments(dform.ctx(), 1.0), h34.step_arguments(dform.ctx(), 1.5)
    assert t1_scope.momentum_shift is own
    (d_lo, m_lo), (d_hi, m_hi) = lo, hi
    assert (d_lo, d_hi) == (harness.B0_DAMPING, 1.5 * harness.B0_DAMPING) == (0.5, 0.75)
    assert sorted(m_lo) == sorted(m_hi) == ["CM", "Endothelial"]  # Tiny is left to the damping
    # the vectors are the harness's own b0_trend factors on D's context: the input stages and the draw only
    merged = np.asarray(pool.obs["merged"]).astype(str)[rows]
    types = sorted(set(np.where(merged == "Untracked", "Untracked@end", merged)))
    logf, _ = harness.b0_trend_factors(dform.ctx(), list(pool.var_names), types)
    assert sorted(logf) == sorted(m_lo)
    for k in m_lo:
        np.testing.assert_array_equal(m_lo[k], logf[k])
        np.testing.assert_array_equal(m_hi[k], 1.5 * logf[k])
        assert m_lo[k].shape == (len(dform.mod.GENES),) and np.abs(m_lo[k]).max() > 0.01
    assert h34.factor_identity(lo, hi) == {"ratio": 1.5, "lineages": ["CM", "Endothelial"], "damping": [0.5, 0.75],
                                           "ok": True, "why": None}
    assert h34.factor_identity(lo, lo)["ok"] is False  # strength 1 against itself is not the identity of 1.5


def test_reading_the_arguments_leaves_the_module_as_it_was_when_the_step_fails(h34, dform):
    own = t1_scope.momentum_shift
    two = dform.mod._view(inputs={"ref": "T2/E9.5.h5ad", "start": "T2/E8.75.h5ad"})  # b0_trend needs three stages
    with pytest.raises(ValueError, match="three input stages"):
        h34.step_arguments(dform.ctx(view=two), 1.0)
    assert t1_scope.momentum_shift is own


def test_the_factor_identity_by_hand_exactly_and_with_no_tolerance(h34):
    rng = np.random.default_rng(0)
    m = {"CM": rng.normal(0.0, 0.4, 30), "Endothelial": rng.normal(0.0, 0.4, 30)}

    def at(s, damping=0.5):
        return damping * s, {k: s * v for k, v in m.items()}
    ok = h34.factor_identity(at(1.0), at(1.5))
    assert ok == {"ratio": 1.5, "lineages": ["CM", "Endothelial"], "damping": [0.5, 0.75], "ok": True, "why": None}
    # one element one unit in the last place off: exactly means exactly
    d, hi = at(1.5)
    hi["CM"][7] = np.nextafter(hi["CM"][7], np.inf)
    bad = h34.factor_identity(at(1.0), (d, hi))
    assert bad["ok"] is False
    assert bad["why"] == "the log-factor vectors of ['CM'] are not 1.5 times the lower strength's"
    # a clip breaks the scale
    assert h34.factor_identity(at(1.0), (0.75, {k: np.clip(1.5 * v, -0.5, 0.5) for k, v in m.items()}))["ok"] is False
    # a damping that did not scale (the lineages the mapping leaves out would stay at strength 1)
    bad = h34.factor_identity(at(1.0), (0.5, at(1.5)[1]))
    assert bad["ok"] is False and bad["why"] == "the damping 0.5 is not 1.5 times the lower strength's 0.5"
    assert h34.factor_identity(at(1.0), (float(np.nextafter(0.75, 1.0)), at(1.5)[1]))["ok"] is False
    # a lineage in one mapping only; vectors of another length; a value that is not a number
    only = h34.factor_identity(at(1.0), (0.75, {"CM": 1.5 * m["CM"]}))
    assert only["ok"] is False and only["why"] == "the two strengths' mappings name different lineages"
    assert h34.factor_identity(at(1.0), (0.75, {k: 1.5 * v[:20] for k, v in m.items()}))["ok"] is False
    d, nan = at(1.5)
    nan["Endothelial"][0] = np.nan
    assert h34.factor_identity(at(1.0), (d, nan))["ok"] is False
    # no lineage in the mapping: the damping alone is held
    assert h34.factor_identity((0.5, {}), (0.75, {}))["ok"] is True
    assert h34.factor_identity((0.5, {}), (0.5, {}))["ok"] is False
    # another ratio
    assert h34.factor_identity(at(1.0), at(2.0), 2.0)["ok"] is True
    assert h34.factor_identity(at(1.0), at(2.0))["ok"] is False


# ------------------------------------------------------------------------------------------------ the loader
@pytest.fixture()
def fake_io(h34, synth, monkeypatch, tmp_path):
    """Every loader of rows H29b's and H34's scripts stood in by the synthetic stages, and the scorer by a recorder:
    nothing is read from disk and nothing is scored. Returns the recorder's log."""
    mod, stages = synth
    h29b, h27 = h34.h29b, h34.h27

    def stage_of(path, *a, **k):
        return stages["T2/" + Path(path).name].copy()

    def task(key):
        return mod._heart_task(key=key, inputs=ROLES[key], reference=ROLES[key]["ref"])
    log = []

    def score(board, pred, truth, ref, seed=0, anchors="computed", build=None, protocol=None):
        log.append({"truth_cells": board.truth_cells, "pred": pred, "truth": truth, "ref": ref, "seed": seed,
                    "anchors": anchors, "build": build, "protocol": protocol})
        return SimpleNamespace(score=50.0 + 1e-3 * len(log), skills={m: 0.5 for m in METRICS},
                               group_skill={g: 0.5 for g in GROUPS}, anchors={"of": [board.truth_cells, seed]})
    monkeypatch.setattr(h29b.payoff, "heart_task", task)
    monkeypatch.setattr(h27, "input_geometries", lambda t, data_dir: (_geoms(h34, stages, h27.input_files(t)), []))
    monkeypatch.setattr(prep, "add_lineage", lambda a, source=None: None)
    monkeypatch.setattr(harness, "task_board", lambda t, boards, n_target, n_ref: mod.BOARD)
    for m in (h34, h29b):
        monkeypatch.setattr(m, "load_boards", lambda: {"T2:heart:val_extrap": mod.BOARD})
        monkeypatch.setattr(m, "read_stage", stage_of)
        monkeypatch.setattr(m, "load_stage", stage_of)
        monkeypatch.setattr(m, "score_board", score)
    monkeypatch.setattr(h29b, "refuse_run", lambda *a, **k: None)
    monkeypatch.setattr(h34, "refuse_state", lambda commit, what="run": None)
    return log


def _call(log, score):
    """The scorer call that returned `score` (the recorder's scores name their call)."""
    return log[round((score - 50.0) / 1e-3) - 1]


def _same_call(a, b):
    for k in ("truth_cells", "seed", "anchors", "build", "protocol"):
        assert a[k] == b[k], k
    for k in ("pred", "truth", "ref"):
        np.testing.assert_array_equal(a[k].X, b[k].X)
        np.testing.assert_array_equal(a[k].coords, b[k].coords)


def test_the_gates_control_is_built_and_scored_exactly_as_row_h29b_builds_and_scores_d_f(h34, synth, fake_io,
                                                                                        monkeypatch, tmp_path):
    """Row H29b's own run step and this row's gate and run steps on the same stand-in stages: the control's cells,
    its size record, its draw, every scorer argument (board, truth, reference, seed, the floor and its anchors) and
    the target-free records are row H29b's for organ:f, bit for bit."""
    log, (_, stages) = fake_io, synth
    monkeypatch.setattr(h34, "refuse_run", lambda *a, **k: None)
    for key in h34.TASK_KEYS:
        theirs = h34.h29b.run_transition(key, (0, 1), tmp_path / "b" / h34.h26.safe(key), tmp_path, COMMIT,
                                         tmp_path / "g")
        n_theirs = len(log)
        tr = h34.load_transition(key, tmp_path)
        mine = [h34.gate_seed(tr, s, tmp_path) for s in (0, 1)]
        assert len(log) - n_theirs == 2 * len(h34.TRUTHS) * 2  # per seed and truth: the floor, the control
        ran = h34.run_transition(key, (0, 1), tmp_path / "r" / h34.h26.safe(key), tmp_path, COMMIT, tmp_path / "g")
        for i in (0, 1):
            their_row = theirs["rows"][i]
            for row, size_rec in ((mine[i], mine[i]["size_record"]),
                                  (ran["rows"][i], ran["rows"][i]["records"]["control"]["size"])):
                assert row["seed"] == their_row["seed"] == i and row["build_seed"] == their_row["build_seed"]
                assert size_rec == their_row["size_records"]["organ:f"]
                assert row["draw"] == their_row["draws"]["organ"]
                assert row["n_cells"]["control"] == their_row["n_cells"]["organ:f"]
                for q in h34.TRUTHS:
                    for mine_name, their_name in (("control", "organ:f"), ("floor", "floor")):
                        _same_call(_call(log, row["scores"][q][mine_name]),
                                   _call(log, their_row["scores"][q][their_name]))
                    # the control is scored against the floor's own anchors, at the truth's own cell number
                    arm_call = _call(log, row["scores"][q]["control"])
                    floor_call = _call(log, row["scores"][q]["floor"])
                    assert floor_call["anchors"] == "computed" and floor_call["truth_cells"] == tr.truths[q].n
                    assert arm_call["anchors"] == {"of": [tr.truths[q].n, i]}
                    assert (arm_call["build"], arm_call["protocol"]) == (harness.BUILD, harness.PROTOCOL)
            assert mine[i]["factors"]["ok"] is True and set(mine[i]["scores"]["own"]) == {"floor", "control"}
        f = tr.facts
        assert f["d_targets"]["organ:f"] == theirs["d_targets"]["organ:f"] == tr.target
        assert f["ratios"] == theirs["ratios"] and f["forecast_targets"] == theirs["forecast_targets"]
        assert f["pred_info"]["organ"] == theirs["pred_info"]["organ"] and f["v7"] == theirs["v7"]
        assert f["truth_info"] == {q: theirs["truth_info"][q] for q in h34.TRUTHS}
        assert f["organ_record"]["reference"] == theirs["organ_record"]["reference"]
        assert f["size_result"] == theirs["size_result"]
        for k, v in f.items():
            assert ran[k] == v, k  # the part carries the same records
        assert tr.inputs["ref"].n_obs == f["pred_info"]["organ"]["pool_cells"] < tr.ref_stage.n
        assert set(tr.inputs) == set(ROLES[key]) and tr.view.reference_role == "ref"
        for role in ("start", "earlier"):  # earlier stages whole, as the builder loads them
            assert tr.inputs[role].n_obs == stages[ROLES[key][role]].n_obs
        assert ran["arms"] == list(h34.BUILT) and ran["task"] == key and ran["axis"] == h34.AXIS
        assert [set(r["scores"][q]) for r in ran["rows"] for q in h34.TRUTHS] == [set(h34.ALL)] * 8
        assert ran["provenance"].keys() == {"commit", "tree_clean", "utc"}
        assert all(r["checks"]["s1_equals_harness_b0_trend"] for r in ran["rows"])


# ------------------------------------------------------------------------------------------------ hand-made tables
HI, LO = "trend_s1.5", "trend_s1.0"
POOL_CELLS, DRAWN = 900, 700
SIZE_REC = {"sizing": "f", "target_rms_full": 3.0, "base_rms": 2.0, "ratio": 1.5, "log_rms_after_vs_target": 0.0}
DRAW = {"pool_cells": POOL_CELLS, "drawn": DRAWN, "short": False, "repeated_rows_at_least": 0}


def _off(v, q, s) -> float:
    """An offset given as a number, ten per-seed numbers, or a dict world -> either."""
    if isinstance(v, dict):
        v = v.get("own" if q == "own" else "CLEAN", 0.0)
    return float(np.broadcast_to(np.asarray(v, dtype=float), (10,))[s])


def _control(h34, key, s, q) -> float:
    return 50.0 + 0.5 * h34.TASK_KEYS.index(key) + 0.01 * s + 0.001 * h34.TRUTHS.index(q)


def _trend(h34, key, name):
    if name == h34.CONTROL:
        return {"strength": None}
    kind, s = name.split("_s")
    s, b = float(s), h34.PRACTICE_BUDGET[key]
    return {"strength": s, "placebo": kind == "placebo", "budget_days": b, "effective_budget_days": s * b,
            "ref_cells": POOL_CELLS, "shifted": ["CM", "Endothelial"],
            "fixed_points": {"CM": 0, "Endothelial": 0} if kind == "placebo" else {}}


def _row(h34, key, s, score=None, de=None):
    """score: arm -> board points added to the control's score; de: arm -> skill added to de_score and de_direction."""
    score, de = score or {}, de or {}
    sc, sk, gr = {}, {}, {}
    for q in h34.TRUTHS:
        base = _control(h34, key, s, q)
        sc[q] = {"floor": 50.0, **{a: base + _off(score.get(a, 0.0), q, s) for a in h34.BUILT}}
        sk[q] = {a: {m: 0.5 + (_off(de.get(a, 0.0), q, s) if m in h34.DE_TERMS else 0.0) for m in METRICS}
                 for a in h34.ALL}
        gr[q] = {a: {g: 0.5 for g in GROUPS} for a in h34.ALL}
    checks = {"s1_equals_harness_b0_trend": True, "same_draw": True, "pre_size_coords_equal": True,
              "post_size_coords_equal": True, "changed_counts": {a: 9 for a in h34.COUNTED},
              "n_cells": {a: 350 for a in h34.BUILT}}
    return {"seed": s, "build_seed": harness.BUILD_SEED_OFFSET + s, "n_cells": {a: 350 for a in h34.BUILT},
            "checks": checks, "records": {a: {"size": dict(SIZE_REC), "trend": _trend(h34, key, a)}
                                          for a in h34.BUILT},
            "draw": dict(DRAW), "scores": sc, "skills": sk, "groups": gr}


def _facts(h34):
    organ = {"rule": "x", "label_config": "c.yaml", "label_config_sha256": "0" * 64,
             "kept_categories": ["cm", "endocardium"], "link_percentile": 99.0, "link_length": 2.5,
             "kept_share": 0.5, "chamber_share": 0.6, "core_share_of_kept": 0.7, "components_kept": 2,
             "components_total": 5}
    return {"truth_info": {q: {"cells": 1000 + i, "core_share": 0.4, "xy_scope": 1.0}
                           for i, q in enumerate(h34.TRUTHS)},
            "pred_info": {"organ": {"pool_cells": POOL_CELLS, "core_share": 0.7, "xy_scope": 0.1}},
            "organ_record": {"reference": organ}, "ratios": {"organ": 1.5, "label": 1.8},
            "d_targets": {"organ:f": {"target": 3.0, "base": 2.0, "ratio": 1.5, "source": "x"}},
            "forecast_targets": {"core_iso_half": 2.0}, "v7": {"ok": True, "log_diff": 0.0}}


def _part(h34, key, seeds, commit=COMMIT, clean=True, **kw):
    return {"axis": h34.AXIS, "task": key, **_facts(h34), "rows": [_row(h34, key, s, **kw) for s in seeds],
            "provenance": {"commit": commit, "tree_clean": clean, "utc": "x"}}


def _parts(h34, **kw):
    return {k: [_part(h34, k, range(0, 5), **kw), _part(h34, k, range(5, 10), **kw)] for k in h34.TASK_KEYS}


def _committed(h34, seeds=range(10)):
    """Row H29b's committed reading as the hand-made control makes it: D_f's plain means by row H29b's own
    arithmetic, and its target-free records."""
    rows = {k: [_row(h34, k, s) for s in seeds] for k in h34.TASK_KEYS}
    f = _facts(h34)
    tables = {w: {u: {"mean": {"organ:f": float(h34.h29b.series(rows, h34.h29b.UNITS[u], ts, "control").mean())}}
                  for u in h34.G1_UNITS} for w, ts in h34.G1_WORLDS}

    def per(v):
        return {k: copy.deepcopy(v) for k in h34.TASK_KEYS}
    short = {"organ": {"pool_cells": POOL_CELLS, "drawn": [DRAWN], "short_in_seeds": 0, "repeated_rows_at_least": 0}}
    return {"row": "H29b", "void": False, "outcome": "READ", "tables": tables, "truth_info": per(f["truth_info"]),
            "pool_info": per(f["pred_info"]), "organ_records": per(f["organ_record"]), "ratios": per(f["ratios"]),
            "d_targets": per(f["d_targets"]), "pools_short": per(short)}


def _gate(h34, verdict="PASS", commit=COMMIT, clean=True):
    return {"row": "H34", "axis": h34.GATE_AXIS, "verdict": verdict, "commit": commit,
            "why": [] if verdict == "PASS" else ["G1: x"], "distances": {},
            "provenance": {"commit": commit, "tree_clean": clean}}


def _read(h34, parts=None, committed="default", gate="default", commit=COMMIT, **kw):
    return h34.reading(_parts(h34, **kw) if parts is None else parts,
                       _committed(h34) if committed == "default" else committed,
                       _gate(h34) if gate == "default" else gate, commit)


NOT_STOPPED_LINES = ("Try H34: NOT STOPPED.",
                     "A screen, not evidence of development. What follows: step (3), the builder option under ruling "
                     "row R32, is coded and reviewed, then the identity builds (step 4), then the file (step 5), each "
                     "as the line's Build passage has it; NOT STOPPED uploads nothing, the file going in only under "
                     "the line's Rule and on the captain's word at the time.")
NOTHING_BUILT = "A screen, not evidence of development. What follows: nothing is built."
GOOD = {HI: 0.6, LO: 0.4}  # 1.5 above 1 and ahead of the copy by the bar, in every seed and both worlds


# ------------------------------------------------------------------------------------------------ the rule
def test_not_stopped_when_1_5_is_above_1_and_ahead_of_the_copy_in_both_worlds(h34):
    res = _read(h34, score=GOOD)
    assert res["why"] == [] and res["void"] is False and res["outcome"] == "NOT STOPPED" and res["clauses"] == []
    st = res["statements"]
    assert st["a_in_world"] == {"own": False, "CLEAN": False} and st["bar_met"] == {"own": True, "CLEAN": True}
    assert st["holds"] == {"a": False, "b": False} and st["stop"] is False
    assert (res["line"], res["follows"]) == NOT_STOPPED_LINES and res["notes"] == []
    for w in ("own", "CLEAN"):
        a = res["tables"][w]["fold E10.5"]["trend_s1.5 - trend_s1.0"]["score"]
        b = res["tables"][w]["fold E10.5"]["trend_s1.5 - control"]["score"]
        assert a["mean"] == pytest.approx(0.2) and a["ahead"] == 10 and a["behind"] == 0
        assert b["mean"] == pytest.approx(0.6) and b["ahead"] == 10 and b["behind"] == 0
    assert set(res["tables"]["own"]) == {KEY1, KEY2, "fold E10.5"}
    assert set(res["tables"]["own"][KEY1]) == {h34.pname(a, b) for a, b in h34.PAIRS} and len(h34.PAIRS) == 7


def _dips(n_below, rise=0.01):
    """trend_s1.5 minus trend_s1.0 per seed: 0.05 below zero in n_below seeds, `rise` above it in the others."""
    return np.where(np.arange(10) < n_below, -0.05, rise)


def test_clause_a_a_mean_below_zero_with_eight_seeds_below_in_either_world(h34):
    def run(own, clean):
        # trend_s1.5 stays 0.6 over the copy in both worlds, so clause (b) is silent; trend_s1.0 is it minus d
        return _read(h34, score={HI: 0.6, LO: {"own": 0.6 - np.asarray(own), "CLEAN": 0.6 - np.asarray(clean)}})
    up = np.full(10, 0.2)
    res = run(_dips(8), up)
    assert res["outcome"] == "STOPPED" and res["clauses"] == ["a"] and res["void"] is False
    assert res["statements"]["a_in_world"] == {"own": True, "CLEAN": False}
    assert res["statements"]["bar_met"] == {"own": True, "CLEAN": True}
    assert res["statements"]["holds"] == {"a": True, "b": False}
    assert res["line"] == "Try H34: STOPPED (a): OTS15 is not built." and res["follows"] == NOTHING_BUILT
    s = res["tables"]["own"]["fold E10.5"]["trend_s1.5 - trend_s1.0"]["score"]
    assert s["mean"] == pytest.approx(-0.038) and s["behind"] == 8 and s["ahead"] == 2
    assert run(up, _dips(8))["clauses"] == ["a"]  # in own or in CLEAN: CLEAN alone stops
    assert run(up, _dips(8))["statements"]["a_in_world"] == {"own": False, "CLEAN": True}
    assert run(_dips(10), _dips(10))["clauses"] == ["a"]
    assert run(_dips(9), up)["clauses"] == ["a"]
    assert run(_dips(7), up)["outcome"] == "NOT STOPPED"  # the mean is below 0, but 7 of 10 seeds are not 8
    assert run(_dips(8, rise=1.0), up)["outcome"] == "NOT STOPPED"  # 8 seeds below, but the mean is above 0
    # a difference of exactly 0 is neither above nor below 0
    tie = run(np.zeros(10), up)
    assert tie["outcome"] == "NOT STOPPED"
    s = tie["tables"]["own"]["fold E10.5"]["trend_s1.5 - trend_s1.0"]["score"]
    assert (s["mean"], s["ahead"], s["behind"]) == (0.0, 0, 0)
    assert run(up, up)["outcome"] == "NOT STOPPED"


def test_clause_b_holds_unless_the_s1_bar_is_met_in_own_and_also_in_clean(h34):
    def run(own, clean):
        # trend_s1.0 stays 0.1 under trend_s1.5 in both worlds, so clause (a) is silent
        own, clean = np.asarray(own, float), np.asarray(clean, float)
        return _read(h34, score={HI: {"own": own, "CLEAN": clean}, LO: {"own": own - 0.1, "CLEAN": clean - 0.1}})
    res = run(0.6, 0.29)
    assert res["outcome"] == "STOPPED" and res["clauses"] == ["b"]
    assert res["statements"]["bar_met"] == {"own": True, "CLEAN": False}
    assert res["statements"]["a_in_world"] == {"own": False, "CLEAN": False}
    assert res["statements"]["holds"] == {"a": False, "b": True}
    assert res["line"] == "Try H34: STOPPED (b): OTS15 is not built." and res["follows"] == NOTHING_BUILT
    assert run(0.29, 0.6)["clauses"] == ["b"]  # the bar failing in one world is enough for clause (b)
    assert run(0.29, 0.6)["statements"]["bar_met"] == {"own": False, "CLEAN": True}
    assert run(0.1, 0.1)["clauses"] == ["b"] and run(-0.5, -0.5)["clauses"] == ["b"]
    assert run(0.31, 0.31)["outcome"] == "NOT STOPPED"
    seven = np.where(np.arange(10) < 7, 1.0, -0.1)  # a mean of +0.67 with 7 of 10 seeds above 0
    eight = np.where(np.arange(10) < 8, 1.0, -0.1)
    assert run(seven, 0.6)["clauses"] == ["b"] and run(0.6, seven)["clauses"] == ["b"]
    assert run(eight, eight)["outcome"] == "NOT STOPPED"


def test_both_clauses_and_the_diagnostics_gate_nothing(h34):
    res = _read(h34, score={HI: 0.1, LO: 0.2})
    assert res["clauses"] == ["a", "b"] and res["line"] == "Try H34: STOPPED (a, b): OTS15 is not built."
    assert res["statements"]["holds"] == {"a": True, "b": True}
    # arm minus placebo and the DE terms: whatever they are, the outcome is as it was
    for plc, de in ((2.0, -0.3), (-2.0, 0.3), (0.6, 0.0)):
        extra = {"placebo_s1.5": plc, "placebo_s1.0": plc}
        ok = _read(h34, score={**GOOD, **extra}, de={HI: de, LO: de, "placebo_s1.5": -de})
        assert ok["outcome"] == "NOT STOPPED"
        d = ok["tables"]["own"]["fold E10.5"]
        assert d["trend_s1.5 - placebo_s1.5"]["score"]["mean"] == pytest.approx(0.6 - plc)
        assert d["trend_s1.5 - control"]["de_score"] == pytest.approx(de)
        assert d["trend_s1.5 - control"]["d2_shape"] == 0.0
        assert _read(h34, score={HI: 0.1, LO: 0.2, **extra}, de={HI: de, LO: de})["clauses"] == ["a", "b"]


def test_the_rules_two_bars_at_their_edges(h34):
    assert h34.below({"mean": -1e-9, "behind": 8}) and h34.below({"mean": -3.0, "behind": 10})
    assert not h34.below({"mean": 0.0, "behind": 10}) and not h34.below({"mean": -1.0, "behind": 7})
    assert not h34.below({"mean": 1e-9, "behind": 9})
    passes = h34.h28.passes  # row H28's S1 bar, in one world
    assert passes({"mean": 0.3000001, "ahead": 8}) and not passes({"mean": 0.3, "ahead": 10})
    assert not passes({"mean": 5.0, "ahead": 7})
    assert h34.outcome_line("VOID", ["v3: no gate record", "v4: missing transition x"]) == \
        "Try H34: VOID: v3: no gate record; v4: missing transition x; nothing is built."
    assert h34.outcome_line("STOPPED", clauses=["a"]) == "Try H34: STOPPED (a): OTS15 is not built."
    assert h34.outcome_line("NOT STOPPED") == "Try H34: NOT STOPPED."
    assert h34.follows_line("VOID") == h34.follows_line("STOPPED") == NOTHING_BUILT
    assert h34.follows_line("NOT STOPPED") == NOT_STOPPED_LINES[1]


def test_the_reading_and_its_text(h34):
    res = _read(h34, score=GOOD)
    json.dumps(res, default=float)
    assert res["row"] == "H34" and res["commit"] == COMMIT and res["gate"]["verdict"] == "PASS"
    assert res["prior"] == "No row on record had read the trend step on D's pool at any strength before this one."
    assert res["guards"]["G1"] == {"means": 18, "unequal": 0, "largest_abs_diff": 0.0}
    assert res["h29b_reading_sha256"] == h34.H29B_READING_SHA256
    assert res["budgets"] == {"1": {"test": 0.5, KEY1: 0.375, KEY2: 0.125},
                              "1.5": {"test": 0.75, KEY1: 0.5625, KEY2: 0.1875}}
    lines = h34.text(res).splitlines()
    assert tuple(lines[:2]) == NOT_STOPPED_LINES
    assert res["prior"] in lines and h34.RULE in lines
    txt = "\n".join(lines)
    for need in ("Row H34 (commit " + COMMIT, "no test score enters the row",
                 f"GATE (first, counting nothing): PASS at commit {COMMIT}.",
                 "Clause (a) in own: does not hold (trend_s1.5 minus trend_s1.0, fold E10.5: +0.20 [+0.20, +0.20] "
                 "10/10 above, 0/10 below", "Clause (a) in CLEAN: does not hold",
                 "The S1 bar in own: MET (trend_s1.5 minus the control D_f, fold E10.5: +0.60",
                 "The S1 bar in CLEAN: MET",
                 "Clause (a) holds in own or in CLEAN: FALSE. Clause (b) holds (the bar is not met in both worlds): "
                 "FALSE.",
                 "DIAGNOSTIC, gates nothing: trend_s1.0 minus placebo_s1.0, fold E10.5",
                 "DIAGNOSTIC, gates nothing: trend_s1.5 minus placebo_s1.5, fold E10.5",
                 "DIAGNOSTIC, gates nothing: the DE terms of trend_s1.0 minus control, fold E10.5",
                 "DIAGNOSTIC, gates nothing: the DE terms of trend_s1.5 minus control, fold E10.5",
                 "G1 repeated on the run's own control, against row H29b's committed D_f: 18 plain means compared, "
                 "0 unequal.", "Beside the read, gating nothing (the target-free records against row H29b's "
                 "committed ones): no difference.", "TABLES, practice diagnostics that gate nothing", "eight tries"):
        assert need in txt, need
    tagged = [ln for ln in lines if ln.startswith(h34.TAG)]
    assert len(tagged) == len(h34.WORLDS) * len(h34.UNITS) * len(h34.PAIRS) == 42
    stopped = h34.text(_read(h34, score={HI: 0.1, LO: 0.2}))
    assert stopped.startswith("Try H34: STOPPED (a, b): OTS15 is not built.\n" + NOTHING_BUILT + "\n")
    for need in ("Clause (a) in own: HOLDS", "Clause (a) in CLEAN: HOLDS", "The S1 bar in own: not met",
                 "Clause (a) holds in own or in CLEAN: TRUE. Clause (b) holds (the bar is not met in both worlds): "
                 "TRUE."):
        assert need in stopped, need
    # the rule's words are the eighth cell's
    for words in ("Clause (a) holds in a world if trend_s1.5 minus trend_s1.0 has a paired mean below 0 with at "
                  "least 8 of the 10 seeds below 0 in that world.",
                  "The S1 bar (row H28's) is met in a world if trend_s1.5 minus the control D_f has a paired mean "
                  "above +0.3 with at least 8 of the 10 seeds above 0 in that world.",
                  "Clause (b) holds unless the S1 bar is met in own and also in CLEAN; the bar failing in one world "
                  "is enough for clause (b).",
                  "STOPPED if clause (a) holds in own or in CLEAN, or clause (b) holds; NOT STOPPED otherwise.",
                  "VOID is decided first, and under VOID the rule is not read.",
                  "a difference of exactly 0 is neither above nor below 0"):
        assert words in h34.RULE, words


# ------------------------------------------------------------------------------------------------ the void
def test_every_void_condition(h34):
    def why(**kw):
        return _read(h34, **{"score": GOOD, **kw})["why"]

    def bad(edit, **kw):
        p = _parts(h34, score=GOOD)
        edit(p)
        return why(parts=p, **kw)

    def committed(edit):
        c = _committed(h34)
        edit(c)
        return why(committed=c)
    assert why() == []
    # v3: the gate fails, or its record is missing or of another commit than the run commit
    assert why(gate=None) == ["v3: no gate record"]
    assert why(gate={"verdict": "PASS", "commit": COMMIT}) == ["v3: no gate record"]  # another row's record
    assert why(gate=_gate(h34, "FAIL")) == ["v3: the gate said 'FAIL', not PASS"]
    assert why(gate=_gate(h34, commit="d" * 40)) == ["v3: the gate record is not of the run commit"]
    assert why(gate=_gate(h34, clean=False)) == ["v3: the gate record is not of the run commit"]
    # v4: a transition, seed, truth or arm is missing or a seed occurs twice
    assert bad(lambda p: p.pop(KEY2)) == [f"v4: missing transition {KEY2}"]
    assert why(parts={}) == [f"v4: missing transition {KEY1}", f"v4: missing transition {KEY2}"]
    assert any(w.startswith(f"v4: {KEY1}: seeds [0, 1, 2, 3, 4, 5, 6, 7, 8], not") for w in
               bad(lambda p: p[KEY1][1]["rows"].pop()))
    assert f"v4: {KEY1}: seeds [0] occur in more than one part" in bad(
        lambda p: p[KEY1][1]["rows"].append(p[KEY1][0]["rows"][0]))
    assert any("missing from scores: ['thin_0/control'" in w for w in
               bad(lambda p: p[KEY2][0]["rows"][0]["scores"].pop("thin_0")))
    assert f"v4: {KEY1} seed 3: missing from skills: ['own/placebo_s1.5']" in bad(
        lambda p: p[KEY1][0]["rows"][3]["skills"]["own"].pop("placebo_s1.5"))
    assert f"v4: {KEY1} seed 3: missing from groups: ['radial_0/floor']" in bad(
        lambda p: p[KEY1][0]["rows"][3]["groups"]["radial_0"].pop("floor"))
    assert f"v4: {KEY1} seed 3: missing from scores: ['own/control']" in bad(
        lambda p: p[KEY1][0]["rows"][3]["scores"]["own"].pop("control"))
    for gone in ("checks", "n_cells"):
        assert f"v4: {KEY1} seed 0: records, checks or cell numbers missing" in bad(
            lambda p, gone=gone: p[KEY1][0]["rows"][0].pop(gone)), gone
    assert f"v4: {KEY1} seed 0: records, checks or cell numbers missing" in bad(
        lambda p: p[KEY1][0]["rows"][0]["records"].pop("trend_s1.0"))
    assert "v4: parts of other transitions ['T2h/E8.25_late->E8.5']" in bad(
        lambda p: p.update({"T2h/E8.25_late->E8.5": p[KEY1]}))
    # v4: the parts come from more than one commit, from a commit other than the run commit or from an unclean tree
    assert "v4: the parts come from more than one commit" in bad(
        lambda p: p[KEY1][0]["provenance"].update(commit="d" * 40))
    assert why(parts=_parts(h34, clean=False, score=GOOD)) == ["v4: a part comes from an unclean tree"]
    assert bad(lambda p: p[KEY2][1]["provenance"].update(tree_clean=None)) == ["v4: a part comes from an unclean tree"]
    for other in ("d" * 40, COMMIT[:12], "", None):
        w = why(commit=other, gate=_gate(h34, commit=other))
        assert any(x.startswith("v4: a part comes from a commit other than the run commit") for x in w), other
    # v3 again: G1 on the run's own control, the control the rule's differences are paired with: bit for bit
    g1 = "v3: G1 on the run's own control: "
    w = bad(lambda p: p[KEY2][1]["rows"][2]["scores"]["thin_0"].update(
        control=p[KEY2][1]["rows"][2]["scores"]["thin_0"]["control"] * (1 + 1e-12)))
    assert len(w) == 4 and all(x.startswith(g1) and x.endswith(" from row H29b's committed D_f") for x in w)
    assert {x[len(g1):].split(": ")[0] for x in w} == {f"CLEAN {KEY2}", "CLEAN fold E10.5", f"truth thin_0 {KEY2}",
                                                       "truth thin_0 fold E10.5"}
    w = committed(lambda c: c["tables"]["own"]["fold E10.5"]["mean"].update(
        {"organ:f": c["tables"]["own"]["fold E10.5"]["mean"]["organ:f"] + 1e-12}))
    assert len(w) == 1 and w[0].startswith(g1 + "own fold E10.5: the control's mean is -")
    for void in (None, {"row": "H29b", "void": True}, {"row": "H27", "void": False, "outcome": "READ"},
                 {"row": "H29b", "void": False, "outcome": "VOID"}):
        assert why(committed=void) == ["v3: row H29b's committed reading is not at hand (missing, not the file its "
                                       "row records, or void)"]
    # v5: an arm is not the arm the third cell defines
    assert bad(lambda p: p[KEY1][0]["rows"][0]["checks"].update(s1_equals_harness_b0_trend=False)) == [
        f"v5: {KEY1} seed 0: trend_s1.0 before the size step is not the harness's b0_trend on the same context"]
    for k in ("same_draw", "pre_size_coords_equal", "post_size_coords_equal"):
        assert bad(lambda p, k=k: p[KEY2][1]["rows"][1]["checks"].pop(k)) == [
            f"v5: {KEY2} seed 6: {k.replace('_', ' ')} fails"]
    assert bad(lambda p: p[KEY1][0]["rows"][0]["n_cells"].update({"trend_s1.5": 351})) == [
        f"v5: {KEY1} seed 0: the arms' cell numbers are not the control's"]
    assert bad(lambda p: p[KEY1][0]["rows"][0]["n_cells"].pop("placebo_s1.0")) == [
        f"v5: {KEY1} seed 0: the arms' cell numbers are not the control's"]
    for term in ("d2_shape", "occupancy_dice", "scale_log_ratio"):  # more than 1e-9 from the control's
        assert bad(lambda p, term=term: p[KEY2][1]["rows"][4]["skills"]["thin_0.5"]["trend_s1.0"].update(
            {term: 0.5 + 2e-9})) == [f"v5: {KEY2} seed 9 thin_0.5: trend_s1.0's {term} is not the control's"]
        assert bad(lambda p, term=term: p[KEY2][1]["rows"][4]["skills"]["thin_0.5"]["trend_s1.0"].update(
            {term: 0.5 + 5e-10})) == []
    rec = "records"
    assert bad(lambda p: p[KEY1][0]["rows"][0][rec]["trend_s1.0"]["trend"].update(budget_days=0.5)) == [
        f"v5: {KEY1} seed 0: trend_s1.0's budget is 0.5, not 0.375"]
    assert bad(lambda p: p[KEY2][0]["rows"][0][rec]["placebo_s1.5"]["trend"].update(
        budget_days=0.125 * (1 + 1e-11))) != []  # relative 1e-12
    assert bad(lambda p: p[KEY2][0]["rows"][0][rec]["placebo_s1.5"]["trend"].update(
        budget_days=0.125 * (1 + 4e-13))) == []
    assert bad(lambda p: p[KEY1][0]["rows"][0][rec]["placebo_s1.5"]["trend"].update(ref_cells=POOL_CELLS + 5)) == [
        f"v5: {KEY1} seed 0: placebo_s1.5's reference is not D's pool"]
    assert bad(lambda p: p[KEY1][0]["rows"][0][rec]["trend_s1.5"]["trend"].update(shifted=["CM"])) == [
        f"v5: {KEY1} seed 0: the arms shift different lineages"]
    assert bad(lambda p: p[KEY1][0]["rows"][0][rec]["placebo_s1.5"]["trend"]["fixed_points"].update(CM=2)) == [
        f"v5: {KEY1} seed 0: placebo_s1.5 leaves a gene on its own factor"]
    # a malformed part or committed record is void, never a crash
    for edit in (lambda p: p[KEY1][0].pop("provenance"), lambda p: p[KEY1][0].pop("pred_info"),
                 lambda p: p[KEY2][0]["rows"][0]["records"]["trend_s1.5"].pop("trend"),
                 lambda p: p[KEY2][0].update(rows=[{"no": "seed"}])):
        w = bad(edit)
        assert len(w) == 1 and w[0].startswith("v4: a part or a committed record cannot be read ("), w
    assert why(committed={"row": "H29b", "void": False, "outcome": "READ"}) == [
        "v4: a part or a committed record cannot be read (KeyError: 'tables')"]


def test_what_the_row_does_not_fix_as_void_is_a_note_and_gates_nothing(h34):
    """The target-free records against row H29b's committed ones, the control's size record and an arm set that
    shifts no lineage: written beside the read, never a void condition and never part of the rule."""
    def read(edit=None, committed=None):
        p = _parts(h34, score=GOOD)
        if edit:
            edit(p)
        return _read(h34, parts=p, committed=committed or _committed(h34))

    def noted(edit, note):
        res = read(edit)
        assert res["why"] == [] and res["outcome"] == "NOT STOPPED" and res["notes"] == [note], res["notes"]
        assert note in h34.text(res)
    assert read()["notes"] == []
    noted(lambda p: p[KEY1][0]["truth_info"]["radial_0"].update(cells=999),
          f"note: {KEY1}: truth radial_0 differs from row H29b's (cells or core share)")
    noted(lambda p: p[KEY2][0]["truth_info"]["own"].update(core_share=0.41),
          f"note: {KEY2}: truth own differs from row H29b's (cells or core share)")
    noted(lambda p: p[KEY1][0]["pred_info"]["organ"].update(core_share=0.71),
          f"note: {KEY1}: D's pool differs from row H29b's (rows or core share)")
    for field, value in (("link_length", 2.5 + 1e-12), ("label_config_sha256", "1" * 64), ("components_kept", 3),
                         ("kept_categories", ["cm"]), ("kept_share", 0.51)):
        noted(lambda p, field=field, value=value: p[KEY1][0]["organ_record"]["reference"].update({field: value}),
              f"note: {KEY1}: D's record on the reference stage differs from row H29b's")
    for pool in ("organ", "label"):
        noted(lambda p, pool=pool: p[KEY2][0]["ratios"].update({pool: 1.5 + 1e-12}),
              f"note: {KEY2}: rms(pool) / rms(core) differs from row H29b's")
    for field in ("target", "base", "ratio"):
        noted(lambda p, field=field: p[KEY1][0]["d_targets"]["organ:f"].update({field: 3.0 + 1e-12}),
              f"note: {KEY1}: D_f's size differs from row H29b's")
    noted(lambda p: p[KEY1][0]["v7"].update(ok=False),
          f"note: {KEY1}: the reference cloud is not the cloud the size result measured")
    noted(lambda p: p[KEY2][0]["rows"][1]["draw"].update(drawn=DRAWN + 1),
          f"note: {KEY2}: the control's draw from D's pool differs from row H29b's")
    noted(lambda p: p[KEY1][0]["rows"][4]["records"]["control"]["size"].update(log_rms_after_vs_target=1e-6),
          f"note: {KEY1} seed 4: the control is not at D_f's size")
    noted(lambda p: p[KEY1][0]["rows"][4]["records"]["control"]["size"].update(target_rms_full=3.1),
          f"note: {KEY1} seed 4: the control's size is not D_f's field")

    def none_shifted(p):
        for a in h34.COUNTED:
            p[KEY1][0]["rows"][0]["records"][a]["trend"]["shifted"] = []
    noted(none_shifted, f"note: {KEY1} seed 0: the arms shift no lineage")
    # prose in row H29b's records is not compared: the rule's words and the config's path may be reworded
    assert read(lambda p: p[KEY1][0]["organ_record"]["reference"].update(rule="other", label_config="x"))["notes"] == []
    assert read(lambda p: p[KEY1][0]["d_targets"]["organ:f"].update(source="other words"))["notes"] == []
    # records that cannot be compared are a note too: the eighteen means are all G1 needs
    res = read(lambda p: p[KEY1][0].pop("v7"))
    assert res["why"] == [] and res["notes"] == ["note: the target-free records cannot be compared (KeyError: 'v7')"]
    bare = _committed(h34)
    for k in ("truth_info", "pool_info", "organ_records", "ratios", "d_targets", "pools_short"):
        bare.pop(k)
    res = read(committed=bare)
    assert res["why"] == [] and res["outcome"] == "NOT STOPPED" and len(res["notes"]) == 1
    assert res["notes"][0].startswith("note: the target-free records cannot be compared (KeyError")
    # the same notes at the gate
    rows, facts = _gate_rows(h34), {k: _facts(h34) for k in h34.TASK_KEYS}
    assert h34.notes(facts, rows, _committed(h34), lambda r: r["size_record"]) == []
    facts[KEY2]["ratios"]["label"] = 1.9
    rows[KEY1][2]["draw"]["short"] = True
    assert h34.notes(facts, rows, _committed(h34), lambda r: r["size_record"]) == [
        f"note: {KEY1}: the control's draw from D's pool differs from row H29b's",
        f"note: {KEY2}: rms(pool) / rms(core) differs from row H29b's"]
    assert h34.notes(facts, rows, None, lambda r: r["size_record"]) == [
        "note: row H29b's committed reading is not at hand; its target-free records are not compared"]


def test_void_wins_over_stop(h34):
    stop = {HI: 0.1, LO: 0.2}
    assert _read(h34, score=stop)["outcome"] == "STOPPED"
    dirty = _parts(h34, score=stop)
    dirty[KEY1][0]["provenance"]["tree_clean"] = False
    nudged = _parts(h34, score=stop)
    nudged[KEY1][0]["rows"][0]["scores"]["own"]["control"] += 1e-9
    for kw in ({"gate": _gate(h34, "FAIL")}, {"gate": None}, {"commit": "d" * 40}, {"committed": None},
               {"parts": dirty}, {"parts": nudged}):
        res = _read(h34, **{"score": stop, **kw})
        assert res["outcome"] == "VOID" and res["void"] is True and res["why"] and res["clauses"] == []
        assert "tables" not in res and "statements" not in res  # the rule is never applied to a void read
        assert res["line"] == "Try H34: VOID: " + "; ".join(res["why"]) + "; nothing is built."
        assert res["follows"] == NOTHING_BUILT
        txt = h34.text(res)
        assert txt.startswith(res["line"] + "\n" + NOTHING_BUILT + "\n")
        # under VOID the rule is not read: its words are there, no clause's finding, no table and no figure of an arm
        for absent in ("Try H34: STOPPED", "Clause (a) in ", "The S1 bar in ", "Clause (a) holds in own or in CLEAN: ",
                       "DIAGNOSTIC", h34.TAG, "TABLES"):
            assert absent not in txt, absent
        assert all(w in txt.splitlines() for w in res["why"]) and res["prior"] in txt and h34.RULE in txt
    res = _read(h34, score=stop, gate=_gate(h34, "FAIL"))
    assert res["line"] == "Try H34: VOID: v3: the gate said 'FAIL', not PASS; nothing is built."
    # every condition that holds is in the line
    res = _read(h34, parts=dirty, gate=None, score=stop)
    assert res["line"] == ("Try H34: VOID: v3: no gate record; v4: a part comes from an unclean tree; nothing is "
                           "built.")


# ------------------------------------------------------------------------------------------------ the gate's verdict
def _gate_rows(h34):
    rows = {}
    for k in h34.TASK_KEYS:
        rows[k] = []
        for s in range(10):
            sc = _row(h34, k, s)["scores"]
            rows[k].append({"seed": s, "scores": {q: {a: sc[q][a] for a in ("floor", "control")} for q in h34.TRUTHS},
                            "size_record": dict(SIZE_REC), "draw": dict(DRAW),
                            "factors": {"ratio": 1.5, "lineages": ["CM", "Endothelial"], "damping": [0.5, 0.75],
                                        "ok": True, "why": None}})
    return rows


def test_the_gates_two_checks(h34):
    why, dist = h34.gate_reasons(_gate_rows(h34), _committed(h34))
    assert why == [] and dist["G1"] == {"means": 18, "unequal": 0, "largest_abs_diff": 0.0}
    assert dist["G2"] == {"contexts": 20, "failed": 0, "lineages": ["CM", "Endothelial"]}
    # G1, exactly and with no tolerance: one seed's control score a part in 1e12 off
    rows = _gate_rows(h34)
    rows[KEY2][3]["scores"]["thin_0.5"]["control"] *= 1 + 1e-12
    why, dist = h34.gate_reasons(rows, _committed(h34))
    assert len(why) == 4 and all(w.startswith("G1: ") and "from row H29b's committed D_f" in w for w in why)
    assert dist["G1"]["unequal"] == 4 and 0 < dist["G1"]["largest_abs_diff"] < 1e-10
    assert dist["G2"]["failed"] == 0
    # G1: a committed mean a part in 1e12 off; each of the eighteen entries is held
    c = _committed(h34)
    c["tables"]["truth radial_0"][KEY1]["mean"]["organ:f"] *= 1 + 1e-12
    why, _ = h34.gate_reasons(_gate_rows(h34), c)
    assert len(why) == 1 and why[0].startswith(f"G1: truth radial_0 {KEY1}: the control's mean is ")
    seen = set()
    for w, _ in h34.G1_WORLDS:
        for u in h34.G1_UNITS:
            c = _committed(h34)
            c["tables"][w][u]["mean"]["organ:f"] = float(np.nextafter(c["tables"][w][u]["mean"]["organ:f"], 100.0))
            why, dist = h34.gate_reasons(_gate_rows(h34), c)
            assert len(why) == 1 and why[0].startswith(f"G1: {w} {u}: ") and dist["G1"]["unequal"] == 1
            seen.add((w, u))
    assert len(seen) == 18
    # G1 is the eighteen means and nothing else: the row's other figures of D_f's table do not enter
    c = _committed(h34)
    c["tables"]["own"][KEY1]["mean"]["core_iso_half:f"] = 1.0
    c["tables"]["truth organ"] = {"fold E10.5": {"mean": {"organ:f": 1.0}}}
    assert h34.gate_reasons(_gate_rows(h34), c)[0] == []
    # G2: one seed's context where the identity fails
    rows = _gate_rows(h34)
    rows[KEY1][7]["factors"] = {"ok": False, "lineages": ["CM"],
                                "why": "the damping 0.5 is not 1.5 times the lower strength's 0.5"}
    why, dist = h34.gate_reasons(rows, _committed(h34))
    assert why == [f"G2: {KEY1} seed 7: the damping 0.5 is not 1.5 times the lower strength's 0.5"]
    assert dist["G2"] == {"contexts": 20, "failed": 1, "lineages": ["CM", "Endothelial"]}
    assert dist["G1"]["unequal"] == 0
    rows = _gate_rows(h34)
    rows[KEY2][0]["factors"] = {}
    assert h34.gate_reasons(rows, _committed(h34))[0] == [f"G2: {KEY2} seed 0: None"]
    # both fail: both are named
    rows = _gate_rows(h34)
    rows[KEY1][0]["scores"]["own"]["control"] += 1e-9
    rows[KEY2][9]["factors"] = {"ok": False, "why": "the two strengths' mappings name different lineages"}
    why, _ = h34.gate_reasons(rows, _committed(h34))
    assert [w[:2] for w in why] == ["G1", "G1", "G1", "G1", "G2"]
    # both transitions and every seed, or the gate fails
    rows = _gate_rows(h34)
    rows.pop(KEY2)
    assert h34.gate_reasons(rows, _committed(h34))[0] == [f"G1: missing transition {KEY2}"]
    rows = _gate_rows(h34)
    rows[KEY1].pop()
    assert h34.gate_reasons(rows, _committed(h34))[0][0].startswith(f"G1: {KEY1}: seeds [0, 1, 2")
    for void in (None, {"row": "H29b", "void": True}):
        assert h34.gate_reasons(_gate_rows(h34), void)[0] == [
            "G1: row H29b's committed reading is not at hand (missing, not the file its row records, or void)"]


def test_row_h29bs_reading_is_the_file_its_row_records_and_no_other(h34, monkeypatch, tmp_path):
    """The reading is loaded only when its sha256 is the one named; any other file, a changed byte, a missing file
    or a file that is not JSON gives nothing, and G1 then cannot hold."""
    f = tmp_path / "reading.json"
    f.write_text(json.dumps(_committed(h34)))
    assert h34.load_h29b_reading(f) is None  # not the file the row records
    monkeypatch.setattr(h34, "H29B_READING_SHA256", hashlib.sha256(f.read_bytes()).hexdigest())
    assert h34.load_h29b_reading(f) == _committed(h34)
    monkeypatch.setattr(h34, "H29B_READING", f)
    assert h34.load_h29b_reading() == _committed(h34)  # the default path is the module's own
    f.write_text(json.dumps(_committed(h34)) + "\n")
    assert h34.load_h29b_reading(f) is None and h34.load_h29b_reading(tmp_path / "none.json") is None
    f.write_text("not json")
    monkeypatch.setattr(h34, "H29B_READING_SHA256", hashlib.sha256(f.read_bytes()).hexdigest())
    assert h34.load_h29b_reading(f) is None
    assert h34.gate_reasons(_gate_rows(h34), None)[0][0].startswith("G1: row H29b's committed reading is not at hand")


def test_the_gate_record_decides_whether_the_run_may_follow(h34, tmp_path):
    def passed(rec, commit=COMMIT):
        (tmp_path / "gate.json").write_text(rec if isinstance(rec, str) else json.dumps(rec))
        return h34.gate_passed(tmp_path, commit)
    assert passed(_gate(h34)) is True
    assert passed(_gate(h34, "FAIL")) is False and passed(_gate(h34, commit="d" * 40)) is False
    assert passed(_gate(h34, clean=False)) is False and passed(_gate(h34), commit=COMMIT[:12]) is False
    assert passed(_gate(h34), commit=None) is False and passed("not json") is False and passed([1]) is False
    assert passed({**_gate(h34), "axis": "h29b-heart-organ:gate"}) is False
    assert passed({**_gate(h34), "provenance": {"commit": "d" * 40, "tree_clean": True}}) is False
    assert h34.gate_passed(tmp_path / "none", COMMIT) is False
    dist = {"G1": {"means": 18, "unequal": 0, "largest_abs_diff": 0.0},
            "G2": {"contexts": 20, "failed": 0, "lineages": ["CM", "Endothelial"]}}
    txt = h34.gate_text({**_gate(h34), "distances": dist, "notes": []})
    for need in ("first, counting nothing", "No counted arm and no placebo is built or scored in it",
                 "the input stages and the draw only", "eighteen entries", "exactly, with no tolerance",
                 f"sha256 {h34.H29B_READING_SHA256}", "G1 found: 18 means compared, 0 unequal",
                 "G2 found: 20 contexts, 0 failed; the mapping's lineages: CM, Endothelial.",
                 "Beside the gate, gating nothing (the target-free records against row H29b's committed ones): no "
                 "difference."):
        assert need in txt, need
    assert txt.rstrip().endswith(h34.GATE_LINES["PASS"])
    noted = h34.gate_text({**_gate(h34), "distances": dist, "notes": ["note: x differs", "note: y differs"]})
    assert "committed ones): note: x differs; note: y differs" in noted and noted.rstrip().endswith(
        h34.GATE_LINES["PASS"])  # a note gates nothing
    failed = h34.gate_text({**_gate(h34, "FAIL"), "why": ["G1: a", "G2: b"]}).splitlines()
    assert failed[-4:] == [h34.GATE_LINES["FAIL"], "WHY: G1: a; G2: b",
                           "Try H34: VOID: G1: a; G2: b; nothing is built.", NOTHING_BUILT]
    assert "the run step refuses, no try is counted and the outcome is VOID" in h34.GATE_LINES["FAIL"]


# ------------------------------------------------------------------------------------------------ the steps, end to end
def _stub_scores(h34, offsets):
    def score_arms(tr, s, preds):
        sc, sk, gr = {}, {}, {}
        for q in h34.TRUTHS:
            base = _control(h34, tr.key, s, q)
            sc[q] = {"floor": 50.0, **{a: base + _off(offsets.get(a, 0.0), q, s) for a in preds}}
            sk[q] = {a: {m: 0.5 for m in METRICS} for a in ("floor", *preds)}
            gr[q] = {a: {g: 0.5 for g in GROUPS} for a in ("floor", *preds)}
        return sc, sk, gr
    return score_arms


def test_gate_run_and_read_end_to_end_on_stand_in_stages(h34, fake_io, monkeypatch, tmp_path, capsys):
    """The three steps through the command line on stand-in stages, two seeds in place of ten (the scores come from a
    stand-in scorer; the builds, the factors, the records and every check are the script's own)."""
    seeds = (0, 1)
    monkeypatch.setattr(h34, "SEEDS", seeds)
    monkeypatch.setattr(h34, "score_arms", _stub_scores(h34, GOOD))
    monkeypatch.setattr(h34, "git_state", lambda: {"commit": COMMIT, "tree_clean": True, "utc": "x"})
    # row H29b's committed reading as these stand-in stages make it: the stand-in control scores, the loader's records
    committed = _committed(h34, seeds)
    for key in h34.TASK_KEYS:
        tr = h34.load_transition(key, tmp_path)
        d = h34.build_control(h34.context(tr, 0), tr.target)[2]
        committed["truth_info"][key], committed["pool_info"][key] = tr.facts["truth_info"], tr.facts["pred_info"]
        committed["organ_records"][key], committed["ratios"][key] = tr.facts["organ_record"], tr.facts["ratios"]
        committed["d_targets"][key] = tr.facts["d_targets"]
        committed["pools_short"][key] = {"organ": {"pool_cells": d["pool_cells"], "drawn": [d["drawn"]],
                                                   "short_in_seeds": len(seeds) * d["short"],
                                                   "repeated_rows_at_least": d["repeated_rows_at_least"]}}
    reading = tmp_path / "h29b_reading.json"

    def commit_reading():
        """Write the stand-in reading and name it as the file row H29b's row records."""
        reading.write_text(json.dumps(committed, default=float))
        monkeypatch.setattr(h34, "H29B_READING_SHA256", hashlib.sha256(reading.read_bytes()).hexdigest())
    commit_reading()
    gate_dir, runs, out = tmp_path / "run" / "gate", tmp_path / "run", tmp_path / "out"
    common = ["--commit", COMMIT, "--data", str(tmp_path)]
    run = ["run", "--task", KEY1, "--seeds", "0-1", "--out", str(runs / "e95" / "0-1"), "--gate", str(gate_dir),
           *common]
    # the run refuses before any gate record exists
    with pytest.raises(SystemExit, match="has not said PASS"):
        h34.main(run)
    # the gate
    with pytest.raises(SystemExit) as e:
        h34.main(["gate", "--out", str(gate_dir), "--h29b", str(reading), *common])
    assert e.value.code == 0
    printed = capsys.readouterr().out.splitlines()
    assert printed[-1] == "GATE PASS" and len(printed) == 5
    assert all(re.fullmatch(r"T2h/[^ ]+ seed [01] \([0-9.]+ s\) done", ln) for ln in printed[:-1])
    g = json.loads((gate_dir / "gate.json").read_text())
    assert g["verdict"] == "PASS" and g["why"] == [] and g["commit"] == COMMIT and g["axis"] == h34.GATE_AXIS
    assert g["distances"]["G1"] == {"means": 18, "unequal": 0, "largest_abs_diff": 0.0}
    assert g["distances"]["G2"] == {"contexts": 4, "failed": 0, "lineages": ["CM", "Endothelial"]}
    assert g["notes"] == [] and g["h29b_reading_sha256"] == h34.H29B_READING_SHA256
    assert set(g["rows"]) == set(h34.TASK_KEYS) and [r["seed"] for r in g["rows"][KEY1]] == [0, 1]
    assert set(g["rows"][KEY1][0]["scores"]["own"]) == {"floor", "control"}  # only the control is built and scored
    assert g["rows"][KEY2][1]["factors"] == {"ratio": 1.5, "lineages": ["CM", "Endothelial"], "damping": [0.5, 0.75],
                                             "ok": True, "why": None}
    assert h34.GATE_LINES["PASS"] in (gate_dir / "gate.txt").read_text() and h34.gate_passed(gate_dir, COMMIT)
    with pytest.raises(SystemExit, match="already holds a gate record"):
        h34.main(["gate", "--out", str(gate_dir), "--h29b", str(reading), *common])
    # the run: two parts, one per transition
    h34.main(run)
    h34.main(["run", "--task", KEY2, "--seeds", "0-1", "--out", str(runs / "e875" / "0-1"), "--gate", str(gate_dir),
              *common])
    printed = capsys.readouterr().out.splitlines()
    assert [bool(re.fullmatch(r"T2h/[^ ]+ seed [01] \([0-9.]+ s\) done", ln)) for ln in printed] == \
        [True, True, False] * 2
    assert printed[2] == f"wrote {runs / 'e95' / '0-1' / 'h34.json'}"
    with pytest.raises(SystemExit, match="exists"):
        h34.main(run)
    parts = h34.load_parts(runs)
    assert set(parts) == set(h34.TASK_KEYS) and all(len(ps) == 1 for ps in parts.values())
    # the read
    read = ["read", "--runs", str(runs), "--gate", str(gate_dir), "--h29b", str(reading), "--out", str(out),
            "--commit", COMMIT]
    h34.main(read)
    res = json.loads((out / "reading.json").read_text())
    assert res["void"] is False and res["why"] == [] and res["notes"] == []
    assert res["guards"]["G1"] == {"means": 18, "unequal": 0, "largest_abs_diff": 0.0}
    # two seeds cannot reach the bar's eight: clause (b) holds here whatever the scores, and the read says so
    assert res["outcome"] == "STOPPED" and res["clauses"] == ["b"]
    assert capsys.readouterr().out == "Try H34: STOPPED (b): OTS15 is not built.\n" + NOTHING_BUILT + "\n"
    assert (out / "reading.txt").read_text().startswith("Try H34: STOPPED (b): OTS15 is not built.\n")
    assert res["trend_record"][KEY1]["ref_cells"] == committed["pool_info"][KEY1]["organ"]["pool_cells"]
    assert all(v > 0 for v in res["changed_counts"][KEY2].values())
    # a committed mean a part in 1e12 off: the gate fails, exits 3 and prints the void line in its fixed words
    committed["tables"]["CLEAN"]["fold E10.5"]["mean"]["organ:f"] *= 1 + 1e-12
    commit_reading()
    with pytest.raises(SystemExit) as e:
        h34.main(["gate", "--out", str(tmp_path / "gate2"), "--h29b", str(reading), *common])
    assert e.value.code == 3
    printed = capsys.readouterr().out.splitlines()
    assert printed[-3] == "GATE FAIL" and printed[-1] == NOTHING_BUILT
    assert printed[-2].startswith("Try H34: VOID: G1: CLEAN fold E10.5: the control's mean is ")
    assert printed[-2].endswith(" from row H29b's committed D_f; nothing is built.")
    g2 = json.loads((tmp_path / "gate2" / "gate.json").read_text())
    assert g2["verdict"] == "FAIL" and len(g2["why"]) == 1 and not h34.gate_passed(tmp_path / "gate2", COMMIT)
    with pytest.raises(SystemExit, match="has not said PASS"):
        h34.main([*run[:6], str(runs / "e95" / "again"), "--gate", str(tmp_path / "gate2"), *common])
    # and a read against that gate is void: the gate failed, and G1 fails on the run's own control too
    h34.main([*read[:4], str(tmp_path / "gate2"), *read[5:8], str(tmp_path / "out2"), *read[9:]])
    line = capsys.readouterr().out.splitlines()[0]
    assert line.startswith("Try H34: VOID: v3: the gate said 'FAIL', not PASS; v3: G1 on the run's own control: "
                           "CLEAN fold E10.5: the control's mean is ") and line.endswith("; nothing is built.")
    # a reading that is not the file row H29b's row records, or a missing one, fails the gate too
    for i, other in enumerate((reading.with_name("other.json"), tmp_path / "none.json")):
        if i == 0:
            other.write_text(reading.read_text() + "\n")
        with pytest.raises(SystemExit) as e:
            h34.main(["gate", "--out", str(tmp_path / f"gate{3 + i}"), "--h29b", str(other), *common])
        assert e.value.code == 3
        g3 = json.loads((tmp_path / f"gate{3 + i}" / "gate.json").read_text())
        assert g3["why"] == ["G1: row H29b's committed reading is not at hand (missing, not the file its row "
                             "records, or void)"]
        assert g3["notes"] == ["note: row H29b's committed reading is not at hand; its target-free records are not "
                               "compared"]
    capsys.readouterr()
    # tracked files or HEAD changing while the gate runs: no gate record at all
    monkeypatch.setattr(h34, "git_state", lambda: {"commit": COMMIT, "tree_clean": False, "utc": "x"})
    with pytest.raises(SystemExit, match="changed while the gate ran; no gate record is written"):
        h34.main(["gate", "--out", str(tmp_path / "gate5"), "--h29b", str(reading), *common])
    assert not (tmp_path / "gate5").exists()


def test_the_read_step_writes_one_reading_and_never_overwrites_it(h34, monkeypatch, tmp_path, capsys):
    monkeypatch.setattr(h34, "refuse_state", lambda commit, what="run": None)
    runs = tmp_path / "runs"
    for k, ps in _parts(h34, score=GOOD).items():
        for i, p in enumerate(ps):
            d = runs / h34.h26.safe(k) / str(i)
            d.mkdir(parents=True)
            (d / h34.PART).write_text(json.dumps(p))
    (runs / "other.json").write_text("{}")
    (runs / "gate").mkdir()
    (runs / "gate" / "gate.json").write_text(json.dumps(_gate(h34)))
    (tmp_path / "h29b.json").write_text(json.dumps(_committed(h34)))
    out = tmp_path / "out"
    args = ["read", "--runs", str(runs), "--gate", str(runs / "gate"), "--h29b", str(tmp_path / "h29b.json"),
            "--out", str(out), "--commit", COMMIT]
    # a reading that is not the file row H29b's row records: the read is void, and writes that
    h34.main([*args[:8], str(tmp_path / "void"), *args[9:]])
    assert capsys.readouterr().out.splitlines()[0] == (
        "Try H34: VOID: v3: row H29b's committed reading is not at hand (missing, not the file its row records, or "
        "void); nothing is built.")
    monkeypatch.setattr(h34, "H29B_READING_SHA256",
                        hashlib.sha256((tmp_path / "h29b.json").read_bytes()).hexdigest())
    h34.main(args)
    assert capsys.readouterr().out == "\n".join(NOT_STOPPED_LINES) + "\n"  # the two lines, and nothing else
    res = json.loads((out / "reading.json").read_text())
    assert res["outcome"] == "NOT STOPPED" and res["void"] is False
    assert (out / "reading.txt").read_text().splitlines()[:2] == list(NOT_STOPPED_LINES)
    before = {f: (out / f).read_bytes() for f in ("reading.json", "reading.txt")}
    with pytest.raises(SystemExit, match="already holds a reading: there is one read under row H34"):
        h34.main(args)
    for f in before:  # either file alone is a reading
        other = out / ("reading.txt" if f == "reading.json" else "reading.json")
        other.rename(tmp_path / "aside")
        with pytest.raises(SystemExit, match="already holds a reading"):
            h34.main(args)
        (tmp_path / "aside").rename(other)
    assert before == {f: (out / f).read_bytes() for f in before}
    with pytest.raises(FileExistsError):
        h34._write_new(out / "reading.txt", "another")
    assert before["reading.txt"] == (out / "reading.txt").read_bytes()
    assert 'd.add_argument("--out", default=str(RESULTS))' in SCRIPT.read_text()


# ------------------------------------------------------------------------------------------------ refusals and prints
def _ledger(tmp_path, rows):
    p = tmp_path / "TRIES.md"
    p.write_text("# Tries\n\n| ID | Board | What |\n|---|---|---|\n" + "".join(f"| {r} | heart | x |\n" for r in rows))
    return p


def _nothing_loads(h34, monkeypatch):
    def boom(*a, **k):
        raise AssertionError("something loaded")
    for name in ("load_transition", "load_parts", "load_h29b_reading"):
        monkeypatch.setattr(h34, name, boom)


def _calls(h34, tmp_path, commit=COMMIT):
    g = tmp_path / "gate"
    g.mkdir(exist_ok=True)
    (g / "gate.json").write_text(json.dumps(_gate(h34)))
    return {"gate": ["gate", "--out", str(tmp_path / "go"), "--commit", commit],
            "run": ["run", "--task", KEY1, "--seeds", "0-1", "--out", str(tmp_path / "o"), "--gate", str(g),
                    "--commit", commit],
            "read": ["read", "--runs", str(tmp_path / "runs"), "--gate", str(g), "--out", str(tmp_path / "r"),
                     "--commit", commit]}


def test_every_step_refuses_without_the_row_before_anything_loads(h34, monkeypatch, tmp_path):
    _nothing_loads(h34, monkeypatch)
    monkeypatch.setattr(h34, "git_state", lambda: {"commit": COMMIT, "tree_clean": True, "utc": "x"})
    calls = _calls(h34, tmp_path)
    # a ledger without the row: other rows, and rows that only mention it, are not row H34
    monkeypatch.setattr(h34, "TRIES", _ledger(tmp_path, ["H28 (one row)", "H29b (the form; see H34)", "H3", "H340",
                                                         "R32 (the captain's order; ruling, no run; H34)"]))
    assert h34.tries_rows() == {"H28", "H29b", "H3", "H340", "R32"}
    for step, args in calls.items():
        with pytest.raises(SystemExit, match="H34 is refused: it has no row in experiments/harness/TRIES.md") as e:
            h34.main(args)
        assert f"a row comes before its {step}" in str(e.value)
    assert not any((tmp_path / d).exists() for d in ("go", "o", "r"))
    # with the row, the same calls pass every refusal and reach their first load
    monkeypatch.setattr(h34, "TRIES", _ledger(tmp_path, ["H28 (one row)", "H34 (one row, one script; eight tries)"]))
    assert "H34" in h34.tries_rows()
    for args in calls.values():
        with pytest.raises(AssertionError, match="something loaded"):
            h34.main(args)


def test_the_other_refusals_of_every_step(h34, monkeypatch, tmp_path):
    _nothing_loads(h34, monkeypatch)
    monkeypatch.setattr(h34, "tries_rows", lambda path=None: {"H34"})
    state = {"commit": COMMIT, "tree_clean": True, "utc": "x"}
    monkeypatch.setattr(h34, "git_state", lambda: dict(state))
    calls = _calls(h34, tmp_path)
    # tracked files that differ from the commit
    for dirty in (False, None):
        state["tree_clean"] = dirty
        for args in calls.values():
            with pytest.raises(SystemExit, match="H34 is refused: tracked files differ from the commit"):
                h34.main(args)
    state["tree_clean"] = True
    # a HEAD other than the commit given
    state["commit"] = "d" * 40
    for args in calls.values():
        with pytest.raises(SystemExit, match="H34 is refused: HEAD is not the run commit"):
            h34.main(args)
    state["commit"] = COMMIT
    for other in (COMMIT[:12], ""):
        for args in _calls(h34, tmp_path, other).values():
            with pytest.raises(SystemExit, match="HEAD is not the run commit"):
                h34.main(args)
    # vec from another checkout
    monkeypatch.setattr(h34, "vec_home", lambda: tmp_path)
    for args in calls.values():
        with pytest.raises(SystemExit, match="H34 is refused: vec is imported from .*, not from this checkout"):
            h34.main(args)
    monkeypatch.setattr(h34, "vec_home", lambda: h34.HERE)
    # the run: the gate record of this commit must say PASS
    g = tmp_path / "gate" / "gate.json"
    for rec in (_gate(h34, "FAIL"), _gate(h34, commit="d" * 40), _gate(h34, clean=False),
                {"verdict": "PASS", "commit": COMMIT}):
        g.write_text(json.dumps(rec))
        with pytest.raises(SystemExit, match="H34 is refused: the gate of .* has not said PASS"):
            h34.main(calls["run"])
    g.unlink()
    with pytest.raises(SystemExit, match="has not said PASS"):
        h34.main(calls["run"])
    g.write_text(json.dumps(_gate(h34)))
    # the run: an existing output, another transition
    (tmp_path / "o").mkdir()
    (tmp_path / "o" / h34.PART).write_text("{}")
    with pytest.raises(SystemExit, match="H34 is refused: .*h34.json exists"):
        h34.main(calls["run"])
    for other in ("T2h/E8.25_late->E8.5", "T2h/other"):
        with pytest.raises(SystemExit, match="--task must be one of"):
            h34.main(["run", "--task", other, "--seeds", "0-1", "--out", str(tmp_path / "o2"), "--gate",
                      str(tmp_path / "gate"), "--commit", COMMIT])
    with pytest.raises(SystemExit, match="not a transition of row H34"):
        h34.run_transition("T2h/E8.25_late->E8.5", (0,), tmp_path / "o3", tmp_path, COMMIT, tmp_path / "gate")
    # the gate and the read: an existing record
    for f in ("gate.json", "gate.txt"):
        (tmp_path / "go").mkdir(exist_ok=True)
        (tmp_path / "go" / f).write_text("x")
        with pytest.raises(SystemExit, match="already holds a gate record"):
            h34.main(calls["gate"])
        (tmp_path / "go" / f).unlink()
    for f in ("reading.json", "reading.txt"):
        (tmp_path / "r").mkdir(exist_ok=True)
        (tmp_path / "r" / f).write_text("x")
        with pytest.raises(SystemExit, match="already holds a reading"):
            h34.main(calls["read"])
        (tmp_path / "r" / f).unlink()
    for args in (calls["gate"], calls["read"]):
        with pytest.raises(AssertionError, match="something loaded"):
            h34.main(args)


def test_the_state_of_the_checkout_is_read_from_the_scripts_own_location(h34, monkeypatch):
    seen = []

    def run(cmd, **kw):
        seen.append(cmd)
        return SimpleNamespace(stdout="" if "status" in cmd else COMMIT + "\n")
    monkeypatch.setattr(h34.h29b.subprocess, "run", run)
    state = h34.git_state()
    assert state["commit"] == COMMIT and state["tree_clean"] is True
    assert all(c[:3] == ["git", "-C", str(REPO_ROOT)] for c in seen)
    assert ["git", "-C", str(REPO_ROOT), "status", "--porcelain", "--untracked-files=no"] in seen


def _prints(fn):
    tree = ast.parse(inspect.getsource(fn).lstrip())
    return [ast.unparse(n) for n in ast.walk(tree) if isinstance(n, ast.Call) and getattr(n.func, "id", "") == "print"]


def test_each_step_prints_only_its_lines(h34):
    run = sorted(_prints(h34.run_transition))
    assert len(run) == 2 and "s) done" in run[0] and run[1] == "print(f'wrote {path}', flush=True)"
    gate = sorted(_prints(h34.run_gate))
    assert len(gate) == 4 and "s) done" in gate[1] and gate[0] == "print(f\"GATE {g['verdict']}\", flush=True)"
    assert gate[2:] == ["print(follows_line(VOID), flush=True)", "print(outcome_line(VOID, why), flush=True)"]
    assert sorted(_prints(h34.cmd_read)) == ["print(res['follows'])", "print(res['line'])"]
    own = [fn for name, fn in inspect.getmembers(h34, inspect.isfunction) if fn.__module__ == h34.__name__]
    assert len(own) > 30
    for fn in own:
        if fn not in (h34.run_transition, h34.run_gate, h34.cmd_read):
            assert _prints(fn) == [], fn.__name__
    # the fixed words hold names (rows, steps, the two strengths) and no figure
    names = r"H29b|H34|R32|G1|G2|step \(3\)|step 4|step 5|strength-1\.5|1\.5 times|strength-1 "
    for words in (*h34.FOLLOWS.values(), h34.SCREEN, h34.NO_PRIOR, *h34.GATE_LINES.values()):
        assert not re.search(r"\d", re.sub(names, "", words)), words


# ================================================================================================ the launcher
# A stand-in for uv: it records each call with the environment the launcher gave it and the run's exit files at that
# moment, and imitates each kind of call the launcher makes.
STUB_UV = r"""#!/usr/bin/env bash
S="$HOME/v/scratch/h34/$STUB_RUN"
threads="$OMP_NUM_THREADS|$OPENBLAS_NUM_THREADS|$MKL_NUM_THREADS"
echo "$UV_PROJECT_ENVIRONMENT|$threads|$NUMPY_MADVISE_HUGEPAGE|$TMPDIR|$*" >> "$HOME/uvcalls"
args=("$@")
[ "${args[0]}" = run ] || exit 90
args=("${args[@]:1}")
while [ "${args[0]}" = --frozen ] || [ "${args[0]}" = --no-sync ]; do args=("${args[@]:1}"); done
[ "${args[0]}" = python ] || exit 91
args=("${args[@]:1}")
while [ "${args[0]}" = -I ] || [ "${args[0]}" = -B ]; do args=("${args[@]:1}"); done
out=""; seeds=""; prev=""
for a in "${args[@]}"; do
  [ "$prev" = "--out" ] && out="$a"
  [ "$prev" = "--seeds" ] && seeds="$a"
  prev="$a"
done
state() { echo "$1 $(cat "$S/$1.exit" 2>/dev/null) $(cat "$S/done.exit" 2>/dev/null)" >> "$HOME/states"; }
case "${args[0]}" in
  -m) state tests; echo "some test output"; exit "${STUB_TESTS_EXIT:-0}" ;;
  experiments/harness/h34_heart_organ_trend.py)
    case "${args[1]}" in
      gate)
        state gate
        [ "${STUB_GATE_EXIT:-0}" = 1 ] && { echo "Traceback"; exit 1; }
        mkdir -p "$out"; echo '{"stub": "gate"}' > "$out/gate.json"
        echo "T2h/<=E9.5->E10.5 seed 0 (1.0 s) done"; echo "a library warning"
        if [ "${STUB_GATE_EXIT:-0}" = 3 ]; then
          echo "GATE FAIL"; echo "Try H34: VOID: G1: x; nothing is built."
          echo "A screen, not evidence of development. What follows: nothing is built."
        else
          echo "GATE PASS"
        fi
        exit "${STUB_GATE_EXIT:-0}" ;;
      run)
        state run
        if [ "${STUB_RUN_BAD:-}" = "$seeds" ]; then echo "boom"; exit 7; fi
        mkdir -p "$out"
        [ "${STUB_RUN_NOFILE:-}" = "$seeds" ] || echo '{"stub": "run"}' > "$out/h34.json"
        echo "T2h/x seed 0 (1.0 s) done"; echo "a library warning"; echo "wrote $out/h34.json"; exit 0 ;;
      read)
        state read
        [ "${STUB_READ_EXIT:-0}" = 0 ] || exit "$STUB_READ_EXIT"
        mkdir -p "$out"; echo "the reading" > "$out/reading.txt"; echo "{}" > "$out/reading.json"
        echo "Try H34: NOT STOPPED."; echo "A screen, not evidence of development. What follows: x"; exit 0 ;;
    esac ;;
esac
exit 92
"""
OUT_DIR = "experiments/harness/results/h34_heart_organ_trend"


def _env(tmp_path, extra=None):
    env = {k: v for k, v in os.environ.items() if not k.startswith(("GIT_", "H34_", "STUB_"))}
    env.update(HOME=str(tmp_path / "home"), GIT_CEILING_DIRECTORIES=str(tmp_path),
               PATH=f"{tmp_path / 'bin'}{os.pathsep}{os.environ['PATH']}")
    env.update(extra or {})
    return env


def _stand_in_env(path):
    """A folder that looks like a virtual environment to the launcher's check (an executable bin/python)."""
    (path / "bin").mkdir(parents=True)
    (path / "bin" / "python").write_text("#!/bin/sh\nexit 0\n")
    (path / "bin" / "python").chmod(0o755)
    return path


@pytest.fixture()
def repo(tmp_path):
    """A throwaway git repository with a ledger that holds a row H34, a stand-in for the script, for the virtual
    environment and for the data folder, a temporary home and the stand-in uv."""
    root = tmp_path / "repo"
    (root / "experiments" / "harness").mkdir(parents=True)
    (tmp_path / "home").mkdir()
    (tmp_path / "bin").mkdir()
    (tmp_path / "bin" / "uv").write_text(STUB_UV)
    (tmp_path / "bin" / "uv").chmod(0o755)
    (root / "experiments" / "harness" / "TRIES.md").write_text("| H29b (a row) | heart | x |\n"
                                                              "| H34 (the row) | heart | x |\n")
    (root / "experiments" / "harness" / "h34_heart_organ_trend.py").write_text("# a stand-in\n")
    (root / ".gitignore").write_text("runs/\n.venv/\ndata/\n")
    _stand_in_env(root / ".venv")
    (root / "data" / "raw" / "T2").mkdir(parents=True)

    def git(*a):
        return subprocess.run(["git", "-c", "user.name=t", "-c", "user.email=t@t", "-c", "init.defaultBranch=x", *a],
                              cwd=root, env=_env(tmp_path), capture_output=True, text=True, check=True).stdout.strip()
    git("init", "-q")
    git("add", ".")
    git("commit", "-q", "-m", "one")
    return SimpleNamespace(root=root, git=git, head=git("rev-parse", "HEAD"), tmp=tmp_path,
                           scratch=tmp_path / "home" / "v" / "scratch" / "h34")


def _launch(repo, *args, cwd=None, **extra):
    run = args[1] if len(args) > 1 else ""
    r = subprocess.run(["bash", str(LAUNCHER), *args], cwd=cwd or repo.root,
                       env=_env(repo.tmp, {"STUB_RUN": run, **extra}), capture_output=True, text=True, timeout=120)
    return r.returncode, r.stderr


def _uv_calls(repo):
    f = repo.tmp / "home" / "uvcalls"
    return f.read_text().splitlines() if f.exists() else []


def _exits(repo, run):
    return {p.name.removesuffix(".exit"): p.read_text().strip() for p in sorted((repo.scratch / run).glob("*.exit"))}


def test_the_launcher_text():
    sh = LAUNCHER.read_text()
    assert subprocess.run(["bash", "-n", str(LAUNCHER)], capture_output=True).returncode == 0
    assert os.access(LAUNCHER, os.X_OK)
    code = "\n".join(ln for ln in sh.splitlines() if not ln.lstrip().startswith("#"))
    assert not re.search(r"\brm\b|\brmdir\b|\bunlink\b|\bmv\b|-delete\b|\bshred\b", sh)  # it deletes nothing
    assert not re.search(r"\bscp\b|\brsync\b|\bcurl\b|\bwget\b|\bssh\b|vec submit|\bupload", sh)  # and sends nothing
    for m in re.finditer(r"\$\{[A-Za-z0-9_]+:\?([^}]*)\}", sh):
        assert "'" not in m.group(1) and '"' not in m.group(1), m.group(0)
    assert len(re.findall(r"\$\{[A-Za-z0-9_]+:\?", sh)) == 2
    assert 'PY="uv run --frozen --no-sync python -I -B"' in code
    assert "VENV=${H34_VENV:-$PWD/.venv}\n" in code and 'export UV_PROJECT_ENVIRONMENT="$VENV"\n' in code
    assert "DATA=${H34_DATA:-data/raw}\n" in code
    assert max(len(ln) for ln in sh.splitlines()) <= 120
    assert "export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 NUMPY_MADVISE_HUGEPAGE=0 " \
           'TMPDIR="$S/tmp"' in code
    assert 'if [ -n "${H34_DATA:-}" ]; then DATA_ARGS=(--data "$H34_DATA"); fi' in code
    order = [code.index(f"STEP={s}\n") for s in ("check", "tests", "gate", "run", "read")]
    assert order == sorted(order) and code.count("STEP=") == 5
    assert code.count("--untracked-files=no") == 1 and "git status --porcelain --untracked-files=no" in code
    # every output of a step goes to a file; nothing is printed but the launcher's own refusals
    assert code.count('> "$S/tests.log" 2>&1') == 1 and code.count('> "$S/gate.out" 2>&1') == 1
    assert code.count('> "$R/logs/$short-$seeds.log" 2>&1') == 1 and code.count('> "$S/read.out" 2>&1') == 1
    assert all(">&2" in ln for ln in code.splitlines() if re.search(r"\becho \"(refused|usage)", ln))
    for test in ("tests/test_h34_heart_organ_trend.py", "tests/test_h29b_heart_organ.py",
                 "tests/test_h27_heart_forecast_size.py"):
        assert test in code, test


def test_the_launcher_runs_the_four_steps_in_order(repo):
    code, err = _launch(repo, repo.head, "r1")
    assert (code, err) == (0, "")
    S, R = repo.scratch / "r1", f"runs/h34_{repo.head[:7]}_r1"
    assert _exits(repo, "r1") == {s: "0" for s in ("check", "tests", "gate", "run", "read", "done")}
    assert (S / "head").read_text().strip() == repo.head and (S / "output").read_text().strip() == R
    # while a step runs, its exit file and done.exit hold the word running
    states = (repo.tmp / "home" / "states").read_text().splitlines()
    assert states == ["tests running running", "gate running running", *["run running running"] * 10,
                      "read running running"]
    calls = _uv_calls(repo)
    env = f"{repo.root}/.venv|8|8|8|0|{S}/tmp|run --frozen --no-sync python -I -B "
    assert len(calls) == 13 and all(c.startswith(env) for c in calls)
    script = "experiments/harness/h34_heart_organ_trend.py"
    assert calls[0] == (env + f"-m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp={S}/tmp/pytest "
                        "tests/test_h34_heart_organ_trend.py tests/test_h29b_heart_organ.py "
                        "tests/test_h27_heart_forecast_size.py")
    assert calls[1] == env + f"{script} gate --out {R}/gate --commit {repo.head}"
    parts = ("0-1", "2-3", "4-5", "6-7", "8-9")
    for block, (key, short) in zip((calls[2:7], calls[7:12]), ((KEY1, "e95"), (KEY2, "e875")), strict=True):
        assert sorted(block) == sorted(env + f"{script} run --task {key} --seeds {p} --gate {R}/gate --commit "
                                       f"{repo.head} --out {R}/{short}/{p}" for p in parts)
    assert calls[12] == env + f"{script} read --runs {R} --gate {R}/gate --commit {repo.head} --out {OUT_DIR}"
    # the parts, their exit files and the counted lines (one line of each log is neither of the two expected kinds)
    assert len(sorted((repo.root / R).rglob("h34.json"))) == 10
    for short in ("e95", "e875"):
        for p in parts:
            assert (repo.root / R / "exit" / f"{short}-{p}").read_text().strip() == "0"
            assert (repo.root / R / "exit" / f"{short}-{p}.other_lines").read_text().strip() == "1"
    assert (repo.root / R / "exit" / "other_lines").read_text().strip() == "10"
    assert (S / "gate.other_lines").read_text().strip() == "1"
    assert (S / "gate.out").read_text().splitlines()[-1] == "GATE PASS"
    assert (S / "read.out").read_text().splitlines() == [
        "Try H34: NOT STOPPED.", "A screen, not evidence of development. What follows: x"]
    assert (S / "tests.log").read_text() == "some test output\n"
    assert (repo.root / OUT_DIR / "reading.txt").read_text() == "the reading\n"
    assert repo.git("status", "--porcelain", "--untracked-files=no") == ""
    # the run folder is never reused, and a second run at the same commit finds the reading: nothing is touched
    before = {p: p.read_bytes() for p in S.rglob("*") if p.is_file()}
    code, err = _launch(repo, repo.head, "r1")
    assert code == 91 and "the run folder" in err and "exists" in err
    assert before == {p: p.read_bytes() for p in S.rglob("*") if p.is_file()} and len(_uv_calls(repo)) == 13
    code, err = _launch(repo, repo.head, "r2")
    assert code == 95 and "already holds a reading" in err
    assert _exits(repo, "r2") == {"check": "95", "done": "95"} and len(_uv_calls(repo)) == 13
    assert (repo.root / OUT_DIR / "reading.txt").read_text() == "the reading\n"


def test_the_launchers_environment_overrides(repo):
    venv, data = repo.tmp / "shared env", repo.tmp / "data raw"
    # a named environment or data folder that is not there: refused before any step
    got, err = _launch(repo, repo.head, "v0", H34_VENV=str(venv))
    assert got == 99 and "no virtual environment at" in err and _exits(repo, "v0") == {"check": "99", "done": "99"}
    got, err = _launch(repo, repo.head, "d0", H34_DATA=str(data))
    assert got == 97 and "no T2 stages under" in err and _exits(repo, "d0") == {"check": "97", "done": "97"}
    assert _uv_calls(repo) == [] and not (repo.root / "runs").exists()
    _stand_in_env(venv)
    (data / "T2").mkdir(parents=True)
    code, err = _launch(repo, repo.head, "e1", H34_VENV=str(venv), H34_DATA=str(data))
    assert (code, err) == (0, "")
    calls = _uv_calls(repo)
    assert len(calls) == 13 and all(c.startswith(f"{venv}|8|8|8|0|") for c in calls)
    with_data = [c for c in calls if c.endswith(f" --data {data}")]
    assert len(with_data) == 11 and calls[0] not in with_data and calls[12] not in with_data  # the gate and the parts
    assert "--data" not in calls[0] and "--data" not in calls[12]


@pytest.mark.parametrize("extra, step, code, n_calls", [
    ({"STUB_TESTS_EXIT": "1"}, "tests", 1, 1),
    ({"STUB_TESTS_EXIT": "5"}, "tests", 5, 1),
    ({"STUB_GATE_EXIT": "3"}, "gate", 3, 2),
    ({"STUB_GATE_EXIT": "1"}, "gate", 1, 2),
    ({"STUB_RUN_BAD": "4-5"}, "run", 1, 7),
    ({"STUB_RUN_NOFILE": "8-9"}, "run", 1, 7),
    ({"STUB_READ_EXIT": "1"}, "read", 1, 13),
])
def test_the_launcher_stops_at_the_first_step_that_does_not_exit_0(repo, extra, step, code, n_calls):
    got, _ = _launch(repo, repo.head, "s1", **extra)
    assert got == code
    steps = ("check", "tests", "gate", "run", "read")
    want = {s: "0" for s in steps[:steps.index(step)]} | {step: str(code), "done": str(code)}
    assert _exits(repo, "s1") == want  # no later step wrote an exit file: none ran
    assert len(_uv_calls(repo)) == n_calls
    R = repo.root / f"runs/h34_{repo.head[:7]}_s1"
    if step == "gate" and code == 3:  # a failed gate: its fixed lines are kept, and counted as expected ones
        S = repo.scratch / "s1"
        assert (S / "gate.out").read_text().splitlines()[-3:] == [
            "GATE FAIL", "Try H34: VOID: G1: x; nothing is built.",
            "A screen, not evidence of development. What follows: nothing is built."]
        assert (S / "gate.other_lines").read_text().strip() == "1" and not (R / "logs").exists()
    if step == "run":  # the failed transition's five parts ran together; the second transition never started
        assert not (R / "e875").exists() and (R / "exit" / "other_lines").exists()
        assert not (repo.root / OUT_DIR).exists()


def test_the_launchers_own_refusals(repo):
    def refused(run, code, words, *args, **kw):
        got, err = _launch(repo, *(args or (repo.head, run)), **kw)
        assert got == code and words in err, (run, got, err)
        assert _exits(repo, run) == {"check": str(code), "done": str(code)}
    # usage, and a RUN that is not a plain name: nothing is written
    for args in ((), (repo.head,), (repo.head, ""), ("", "x")):
        got, err = _launch(repo, *args)
        assert got == 1 and "usage: run_h34.sh COMMIT RUN" in err
    for odd in ("a/b", "..", "a b", "r.1", "x;y"):
        got, err = _launch(repo, repo.head, odd)
        assert got == 90 and "letters, digits" in err
    assert not repo.scratch.exists() or not list(repo.scratch.iterdir())
    refused("h1", 96, "HEAD is not", "d" * 40, "h1")
    refused("h2", 96, "HEAD is not", repo.head[:12], "h2")
    (repo.root / ".gitignore").write_text("runs/\n# edited\n")
    refused("c1", 98, "tracked files differ from the commit")
    repo.git("checkout", "-q", ".gitignore")
    (repo.root / "untracked.txt").write_text("x\n")  # an untracked file is not a tracked change
    refused("n1", 94, "run from the root of the checkout", cwd=repo.root / "experiments")
    ledger = repo.root / "experiments" / "harness" / "TRIES.md"
    ledger.write_text("| H29b (a row) | heart | x |\n| H340 | heart | x |\n| R32 (see H34) | x |\n")
    repo.git("commit", "-q", "-am", "no row")
    head = repo.git("rev-parse", "HEAD")
    refused("t1", 93, "row H34 is not in experiments/harness/TRIES.md", head, "t1")
    ledger.write_text("| H34 (the row) | heart | x |\n")
    repo.git("commit", "-q", "-am", "the row")
    head = repo.git("rev-parse", "HEAD")
    (repo.root / "runs" / f"h34_{head[:7]}_o1").mkdir(parents=True)
    refused("o1", 92, "the output folder", head, "o1")
    assert _uv_calls(repo) == []  # no refusal reached a step
    for f in ("reading.json", "reading.txt"):
        (repo.root / OUT_DIR).mkdir(parents=True, exist_ok=True)
        (repo.root / OUT_DIR / f).write_text("the one reading\n")
        refused(f"k-{f[8:]}".replace(".", ""), 95, "already holds a reading", head, f"k-{f[8:]}".replace(".", ""))
        assert (repo.root / OUT_DIR / f).read_text() == "the one reading\n"
        (repo.root / OUT_DIR / f).rename(repo.tmp / f"aside-{f}")
    assert _uv_calls(repo) == []
    assert _launch(repo, head, "ok")[0] == 0  # and with nothing in the way it runs, the untracked file still there
    assert (repo.root / "untracked.txt").exists()
