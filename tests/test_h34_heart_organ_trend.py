"""Row H34 (experiments/harness/h34_heart_organ_trend.py) and its launcher (experiments/harness/run_h34.sh): synthetic
arrays, hand-made tables and temporary files only. Nothing here loads a competition file, scores a real cell or reads
the repository's own ledger or a committed result; the launcher runs in throwaway git repositories under the test's
own folder, with a temporary home and a stand-in uv first on PATH."""

from __future__ import annotations

import ast
import copy
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
    assert h34.GATE_EXIT == {"PASS": 0, "FAIL": 3} and h34.FACTOR_LOG_TOL == 1e-12
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
                "boards/", "upload", "12.5", "manual_uploads", "index.json", "test_extrap", "size_or_hold(",
                "headline", "halve.py"):
        assert bad not in src, bad


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
def test_the_gate_reads_the_steps_factors_and_builds_no_arm(h34, dform, monkeypatch):
    _, _, rows, rec = h34.h28.trend_expression(dform.ctx(), 1.0)
    fallback, shifted = rec["fallback"], rec["shifted"]
    assert fallback == ["Tiny"] and rec["n_trend"] == 2  # the mapping's two lineages and one through the damping
    own = t1_scope.apply_factors

    def boom(*a, **k):
        raise AssertionError("an arm was built")
    monkeypatch.setattr(harness, "finish", boom)  # the step never reaches the stage: no count changes
    f1, f15 = h34.step_factors(dform.ctx(), 1.0), h34.step_factors(dform.ctx(), 1.5)
    assert t1_scope.apply_factors is own
    pool = dform.inputs["ref"]
    merged = np.asarray(pool.obs["merged"]).astype(str)[rows]
    types = sorted(set(np.where(merged == "Untracked", "Untracked@end", merged)))
    assert f1.shape == f15.shape == (len(types), len(dform.mod.GENES))
    moved = [t for t, row in zip(types, f1, strict=True) if not (row == 1.0).all()]
    assert moved == shifted and "Untracked@end" not in moved
    got = h34.factor_identity(f1, f15)
    assert got["ok"] is True and got["why"] is None and got["max_abs_dev"] <= 1e-12 and got["ratio"] == 1.5
    assert got["moved_types"] == 3 and got["max_abs_log_factor"] > 0.01
    np.testing.assert_allclose(np.log(f15), 1.5 * np.log(f1), rtol=0, atol=1e-12)
    # the factors are B0's and the trend's own: strength 1 hands apply_factors the harness's b0_trend factors
    logf, _ = harness.b0_trend_factors(dform.ctx(), list(pool.var_names), types)
    for t in ("CM", "Endothelial"):
        np.testing.assert_allclose(np.log(f1[types.index(t)]), logf[t], rtol=0, atol=1e-12)
    # strength 1 against itself is not the identity of 1.5
    assert h34.factor_identity(f1, f1)["ok"] is False


def test_reading_the_factors_leaves_the_module_as_it_was_when_the_step_fails(h34, dform):
    own = t1_scope.apply_factors
    two = dform.mod._view(inputs={"ref": "T2/E9.5.h5ad", "start": "T2/E8.75.h5ad"})  # b0_trend needs three stages
    with pytest.raises(ValueError, match="three input stages"):
        h34.step_factors(dform.ctx(view=two), 1.0)
    assert t1_scope.apply_factors is own


def test_the_factor_identity_by_hand(h34):
    rng = np.random.default_rng(0)
    logf = rng.normal(0.0, 0.4, (4, 30))
    logf[3] = 0.0  # a type the step does not shift
    lo, hi = np.exp(logf), np.exp(1.5 * logf)
    ok = h34.factor_identity(lo, hi)
    assert ok["ok"] is True and ok["moved_types"] == 3 and ok["tol"] == 1e-12
    clipped = h34.factor_identity(lo, np.exp(np.clip(1.5 * logf, -0.5, 0.5)))  # a clip breaks the scale
    assert clipped["ok"] is False and "1.5 times" in clipped["why"]
    off = 1.5 * logf.copy()
    off[1, 7] += 1e-9
    assert h34.factor_identity(lo, np.exp(off))["ok"] is False
    off[1, 7] -= 1e-9 - 1e-14
    assert h34.factor_identity(lo, np.exp(off))["ok"] is True
    damping_only = 1.5 * logf.copy()
    damping_only[2] = logf[2]  # one lineage left at strength 1 (a damping that did not scale)
    assert h34.factor_identity(lo, np.exp(damping_only))["ok"] is False
    assert h34.factor_identity(lo, hi[:3])["ok"] is False and "shape" in h34.factor_identity(lo, hi[:3])["why"]
    assert h34.factor_identity(np.ones((0, 5)), np.ones((0, 5)))["ok"] is False
    still = h34.factor_identity(np.ones((3, 5)), np.ones((3, 5)))
    assert still["ok"] is False and "no factor moves" in still["why"]
    for bad in (np.nan, 0.0, -1.0, np.inf):
        x = hi.copy()
        x[0, 0] = bad
        assert h34.factor_identity(lo, x)["ok"] is False
    assert h34.factor_identity(lo, np.exp(2.0 * logf), 2.0)["ok"] is True


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
    monkeypatch.setattr(h34, "refuse_run", lambda *a, **k: None)
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
                                                                                        tmp_path):
    """Row H29b's own run step and this row's gate and run steps on the same stand-in stages: the control's cells,
    its size record, its draw, every scorer argument (board, truth, reference, seed, the floor and its anchors) and
    the target-free records are row H29b's for organ:f, bit for bit."""
    log, (_, stages) = fake_io, synth
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
                    arm_call, floor_call = _call(log, row["scores"][q]["control"]), _call(log, row["scores"][q]["floor"])
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
