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
                     "A screen, not evidence of development. What follows: the builder option under ruling row R32 "
                     "may be coded and reviewed.")
NOTHING_BUILT = "A screen, not evidence of development. What follows: nothing is built."
GOOD = {HI: 0.6, LO: 0.4}  # 1.5 above 1 and ahead of the copy by the bar, in every seed and both worlds


# ------------------------------------------------------------------------------------------------ the rule
def test_not_stopped_when_1_5_is_above_1_and_ahead_of_the_copy_in_both_worlds(h34):
    res = _read(h34, score=GOOD)
    assert res["why"] == [] and res["void"] is False and res["outcome"] == "NOT STOPPED" and res["clauses"] == []
    st = res["statements"]
    assert st["a_below"] == {"own": False, "CLEAN": False} and st["b_ahead"] == {"own": True, "CLEAN": True}
    assert st["fires"] == {"a": False, "b": False} and st["stop"] is False
    assert (res["line"], res["follows"]) == NOT_STOPPED_LINES
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
    assert res["statements"]["a_below"] == {"own": True, "CLEAN": False}
    assert res["statements"]["b_ahead"] == {"own": True, "CLEAN": True}
    assert res["line"] == "Try H34: STOPPED (a): OTS15 is not built." and res["follows"] == NOTHING_BUILT
    s = res["tables"]["own"]["fold E10.5"]["trend_s1.5 - trend_s1.0"]["score"]
    assert s["mean"] == pytest.approx(-0.038) and s["behind"] == 8 and s["ahead"] == 2
    assert run(up, _dips(8))["clauses"] == ["a"]  # EITHER world: CLEAN alone stops
    assert run(up, _dips(8))["statements"]["a_below"] == {"own": False, "CLEAN": True}
    assert run(_dips(10), _dips(10))["clauses"] == ["a"]
    assert run(_dips(9), up)["clauses"] == ["a"]
    assert run(_dips(7), up)["outcome"] == "NOT STOPPED"  # the mean is below 0, but 7 of 10 seeds are not 8
    assert run(_dips(8, rise=1.0), up)["outcome"] == "NOT STOPPED"  # 8 seeds below, but the mean is above 0
    assert run(np.zeros(10), up)["outcome"] == "NOT STOPPED"  # a tie is not below
    assert run(up, up)["outcome"] == "NOT STOPPED"


def test_clause_b_not_ahead_of_the_copy_by_the_s1_bar_in_both_worlds(h34):
    def run(own, clean):
        # trend_s1.0 stays 0.1 under trend_s1.5 in both worlds, so clause (a) is silent
        own, clean = np.asarray(own, float), np.asarray(clean, float)
        return _read(h34, score={HI: {"own": own, "CLEAN": clean}, LO: {"own": own - 0.1, "CLEAN": clean - 0.1}})
    res = run(0.6, 0.29)
    assert res["outcome"] == "STOPPED" and res["clauses"] == ["b"]
    assert res["statements"]["b_ahead"] == {"own": True, "CLEAN": False}
    assert res["statements"]["a_below"] == {"own": False, "CLEAN": False}
    assert res["line"] == "Try H34: STOPPED (b): OTS15 is not built." and res["follows"] == NOTHING_BUILT
    assert run(0.29, 0.6)["clauses"] == ["b"]  # the bar must hold in own and in CLEAN: one world failing stops
    assert run(0.29, 0.6)["statements"]["b_ahead"] == {"own": False, "CLEAN": True}
    assert run(0.1, 0.1)["clauses"] == ["b"] and run(-0.5, -0.5)["clauses"] == ["b"]
    assert run(0.31, 0.31)["outcome"] == "NOT STOPPED"
    seven = np.where(np.arange(10) < 7, 1.0, -0.1)  # a mean of +0.67 with 7 of 10 seeds above 0
    eight = np.where(np.arange(10) < 8, 1.0, -0.1)
    assert run(seven, 0.6)["clauses"] == ["b"] and run(0.6, seven)["clauses"] == ["b"]
    assert run(eight, eight)["outcome"] == "NOT STOPPED"


def test_both_clauses_and_the_diagnostics_gate_nothing(h34):
    res = _read(h34, score={HI: 0.1, LO: 0.2})
    assert res["clauses"] == ["a", "b"] and res["line"] == "Try H34: STOPPED (a, b): OTS15 is not built."
    assert res["statements"]["fires"] == {"a": True, "b": True}
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
    assert h34.outcome_line("VOID", ["G: no gate record", "V1: no part"]) == \
        "Try H34: VOID: G: no gate record; V1: no part; nothing is built."
    assert h34.outcome_line("STOPPED", clauses=["a"]) == "Try H34: STOPPED (a): OTS15 is not built."
    assert h34.outcome_line("NOT STOPPED") == "Try H34: NOT STOPPED."
    assert h34.follows_line("VOID") == h34.follows_line("STOPPED") == NOTHING_BUILT
    assert h34.follows_line("NOT STOPPED") == NOT_STOPPED_LINES[1]


def test_the_reading_and_its_text(h34):
    res = _read(h34, score=GOOD)
    json.dumps(res, default=float)
    assert res["row"] == "H34" and res["commit"] == COMMIT and res["gate"]["verdict"] == "PASS"
    assert res["prior"] == "No row on record had read the trend step on D's pool at any strength before this one."
    assert res["guards"]["V4"] == {"means": 18, "unequal": 0, "largest_abs_diff": 0.0}
    assert res["budgets"] == {"1": {"test": 0.5, KEY1: 0.375, KEY2: 0.125},
                              "1.5": {"test": 0.75, KEY1: 0.5625, KEY2: 0.1875}}
    lines = h34.text(res).splitlines()
    assert tuple(lines[:2]) == NOT_STOPPED_LINES
    assert res["prior"] in lines and h34.RULE in lines
    txt = "\n".join(lines)
    for need in ("Row H34 (commit " + COMMIT, "no test score enters the row", "GATE (first; not a try): PASS",
                 "(a) own: FALSE (trend_s1.5 minus trend_s1.0, fold E10.5: +0.20 [+0.20, +0.20] 10/10 above, "
                 "0/10 below", "(a) CLEAN: FALSE", "(a) fires: FALSE", "(b) own: the bar HOLDS (trend_s1.5 minus "
                 "the copy D_f, fold E10.5: +0.60", "(b) CLEAN: the bar HOLDS", "(b) fires: FALSE",
                 "DIAGNOSTIC, gates nothing: trend_s1.0 minus placebo_s1.0, fold E10.5",
                 "DIAGNOSTIC, gates nothing: trend_s1.5 minus placebo_s1.5, fold E10.5",
                 "DIAGNOSTIC, gates nothing: the DE terms of trend_s1.5 minus control, fold E10.5",
                 "DIAGNOSTIC, gates nothing: the DE terms of trend_s1.5 minus trend_s1.0, fold E10.5",
                 "18 plain means compared, 0 unequal", "eight tries"):
        assert need in txt, need
    tagged = [ln for ln in lines if ln.startswith(h34.TAG)]
    assert len(tagged) == len(h34.WORLDS) * len(h34.UNITS) * len(h34.PAIRS) == 42
    stopped = h34.text(_read(h34, score={HI: 0.1, LO: 0.2}))
    assert stopped.startswith("Try H34: STOPPED (a, b): OTS15 is not built.\n" + NOTHING_BUILT + "\n")
    for need in ("(a) own: TRUE", "(a) fires: TRUE", "(b) own: the bar FAILS", "(b) fires: TRUE"):
        assert need in stopped, need
    assert "in both worlds" in h34.RULE and "in EITHER world" in h34.RULE and "VOID wins over STOP" in h34.RULE
    assert "above +0.3 with at least 8 of 10 seeds above 0" in h34.RULE


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
    # G: the gate record
    assert why(gate=None) == ["G: no gate record"]
    assert why(gate={"verdict": "PASS", "commit": COMMIT}) == ["G: no gate record"]  # another row's record
    assert why(gate=_gate(h34, "FAIL")) == ["G: the gate said 'FAIL', not PASS"]
    assert why(gate=_gate(h34, commit="d" * 40)) == ["G: the gate record is not of the run commit"]
    assert why(gate=_gate(h34, clean=False)) == ["G: the gate record comes from a dirty tree"]
    # V1: a transition, seed, truth or arm missing; a seed twice; another transition's part
    assert bad(lambda p: p.pop(KEY2)) == [f"V1: missing transition {KEY2}"]
    assert why(parts={}) == [f"V1: missing transition {KEY1}", f"V1: missing transition {KEY2}"]
    assert any(w.startswith(f"V1: {KEY1}: seeds [0, 1, 2, 3, 4, 5, 6, 7, 8], not") for w in
               bad(lambda p: p[KEY1][1]["rows"].pop()))
    assert any("occur in more than one part" in w for w in
               bad(lambda p: p[KEY1][1]["rows"].append(p[KEY1][0]["rows"][0])))
    assert any("missing from scores: ['thin_0/control'" in w for w in
               bad(lambda p: p[KEY2][0]["rows"][0]["scores"].pop("thin_0")))
    assert any("missing from skills: ['own/placebo_s1.5']" in w for w in
               bad(lambda p: p[KEY1][0]["rows"][3]["skills"]["own"].pop("placebo_s1.5")))
    assert any("missing from groups: ['radial_0/floor']" in w for w in
               bad(lambda p: p[KEY1][0]["rows"][3]["groups"]["radial_0"].pop("floor")))
    for gone in ("draw", "checks", "n_cells"):
        assert any("records, checks, draw or cell numbers missing" in w for w in
                   bad(lambda p, gone=gone: p[KEY1][0]["rows"][0].pop(gone))), gone
    assert any("records, checks, draw or cell numbers missing" in w for w in
               bad(lambda p: p[KEY1][0]["rows"][0]["records"].pop("trend_s1.0")))
    assert any("parts of other transitions ['T2h/E8.25_late->E8.5']" in w for w in
               bad(lambda p: p.update({"T2h/E8.25_late->E8.5": p[KEY1]})))
    # V2: the commit and the tree
    assert "V2: the parts come from different commits" in bad(
        lambda p: p[KEY1][0]["provenance"].update(commit="d" * 40))
    assert why(parts=_parts(h34, clean=False, score=GOOD)) == ["V2: a part comes from a dirty tree"]
    assert "V2: a part comes from a dirty tree" in bad(lambda p: p[KEY2][1]["provenance"].update(tree_clean=None))
    for other in ("d" * 40, COMMIT[:12], "", None):
        w = why(commit=other, gate=_gate(h34, commit=other))
        assert any(x.startswith("V2: a part does not come from the row's run commit") for x in w), other
    # V3: the transition's target-free records against row H29b's committed reading
    for void in (None, {"row": "H29b", "void": True}, {"row": "H27", "void": False, "outcome": "READ"},
                 {"row": "H29b", "void": False, "outcome": "VOID"}):
        assert why(committed=void) == ["V3: row H29b's committed reading is void or unreadable"]
    assert f"V3: {KEY1}: truth radial_0 differs from row H29b's (cells or core share)" in bad(
        lambda p: p[KEY1][0]["truth_info"]["radial_0"].update(cells=999))
    assert f"V3: {KEY2}: truth own differs from row H29b's (cells or core share)" in bad(
        lambda p: p[KEY2][1]["truth_info"]["own"].update(core_share=0.41))
    assert f"V3: {KEY1}: D's pool differs from row H29b's (rows or core share)" in bad(
        lambda p: p[KEY1][0]["pred_info"]["organ"].update(pool_cells=POOL_CELLS + 1))
    for field, value in (("link_length", 2.5 + 1e-12), ("label_config_sha256", "1" * 64), ("components_kept", 3),
                         ("kept_categories", ["cm"]), ("kept_share", 0.51)):
        assert f"V3: {KEY1}: D's record on the reference stage differs from row H29b's" in bad(
            lambda p, field=field, value=value: p[KEY1][0]["organ_record"]["reference"].update({field: value})), field
    for pool in ("organ", "label"):
        assert f"V3: {KEY2}: rms(pool) / rms(core) differs from row H29b's" in bad(
            lambda p, pool=pool: p[KEY2][0]["ratios"].update({pool: 1.5 + 1e-12}))
    for field in ("target", "base", "ratio"):
        assert f"V3: {KEY1}: D_f's size differs from row H29b's" in bad(
            lambda p, field=field: p[KEY1][1]["d_targets"]["organ:f"].update({field: 3.0 + 1e-12})), field
    assert f"V3: {KEY1}: the reference cloud is not the cloud the size result measured" in bad(
        lambda p: p[KEY1][0]["v7"].update(ok=False))
    assert f"V3: {KEY2}: the control's draw from D's pool differs from row H29b's" in bad(
        lambda p: p[KEY2][0]["rows"][1]["draw"].update(drawn=DRAWN + 1))
    assert f"V3: {KEY1}: truth thin_0 differs from row H29b's (cells or core share)" in committed(
        lambda c: c["truth_info"][KEY1]["thin_0"].update(cells=1))
    # prose in row H29b's records is not compared: the rule's words and the config's path may be reworded
    assert bad(lambda p: p[KEY1][0]["organ_record"]["reference"].update(rule="other words", label_config="x")) == []
    assert bad(lambda p: p[KEY1][0]["d_targets"]["organ:f"].update(source="other words")) == []
    # V4: the run's own control against row H29b's committed D_f, bit for bit
    w = bad(lambda p: p[KEY2][1]["rows"][2]["scores"]["thin_0"].update(
        control=p[KEY2][1]["rows"][2]["scores"]["thin_0"]["control"] * (1 + 1e-12)))
    assert len(w) == 4 and all(x.startswith("V4: ") for x in w)
    assert {x.split(":")[1].strip() for x in w} == {f"CLEAN {KEY2}", "CLEAN fold E10.5", f"truth thin_0 {KEY2}",
                                                   "truth thin_0 fold E10.5"}
    w = committed(lambda c: c["tables"]["own"]["fold E10.5"]["mean"].update(
        {"organ:f": c["tables"]["own"]["fold E10.5"]["mean"]["organ:f"] + 1e-12}))
    assert len(w) == 1 and w[0].startswith("V4: own fold E10.5: the control's mean is -")
    # V5: the identities
    assert any("strength 1 is not the harness's b0_trend" in w for w in
               bad(lambda p: p[KEY1][0]["rows"][0]["checks"].update(s1_equals_harness_b0_trend=False)))
    for k in ("same_draw", "pre_size_coords_equal", "post_size_coords_equal"):
        assert f"V5: {KEY2} seed 6: {k.replace('_', ' ')} fails" in bad(
            lambda p, k=k: p[KEY2][1]["rows"][1]["checks"].pop(k))
    assert f"V5: {KEY1} seed 0: the arms' cell numbers differ" in bad(
        lambda p: p[KEY1][0]["rows"][0]["n_cells"].update({"trend_s1.5": 351}))
    assert f"V5: {KEY1} seed 0: the arms' cell numbers differ" in bad(
        lambda p: p[KEY1][0]["rows"][0]["n_cells"].pop("placebo_s1.0"))
    assert f"V5: {KEY1} seed 4: control is not at D_f's size" in bad(
        lambda p: p[KEY1][0]["rows"][4]["records"]["control"]["size"].update(log_rms_after_vs_target=1e-6))
    assert f"V5: {KEY1} seed 4: trend_s1.5's size is not D_f's field" in bad(
        lambda p: p[KEY1][0]["rows"][4]["records"]["trend_s1.5"]["size"].update(target_rms_full=3.1))
    assert f"V5: {KEY2} seed 9 thin_0.5: trend_s1.0's d2_shape is not the control's" in bad(
        lambda p: p[KEY2][1]["rows"][4]["skills"]["thin_0.5"]["trend_s1.0"].update(d2_shape=0.5 + 1e-6))
    # V6: the trend records
    rec = "records"
    assert any("trend_s1.5's record is not strength 1.5 trend" in w for w in
               bad(lambda p: p[KEY1][0]["rows"][0][rec]["trend_s1.5"]["trend"].update(strength=1.0)))
    assert any("placebo_s1.0's record is not strength 1 placebo" in w for w in
               bad(lambda p: p[KEY1][0]["rows"][0][rec]["placebo_s1.0"]["trend"].update(placebo=False)))
    assert any("trend_s1.0's budget is 0.5, not 0.375" in w for w in
               bad(lambda p: p[KEY1][0]["rows"][0][rec]["trend_s1.0"]["trend"].update(budget_days=0.5)))
    assert any("trend_s1.5's budget at its strength is not 1.5 times 0.125" in w for w in
               bad(lambda p: p[KEY2][0]["rows"][0][rec]["trend_s1.5"]["trend"].update(effective_budget_days=0.125)))
    assert any("placebo_s1.5's reference is not D's pool" in w for w in
               bad(lambda p: p[KEY1][0]["rows"][0][rec]["placebo_s1.5"]["trend"].update(ref_cells=POOL_CELLS + 5)))
    assert any("the arms shift different lineages, or none" in w for w in
               bad(lambda p: p[KEY1][0]["rows"][0][rec]["trend_s1.5"]["trend"].update(shifted=["CM"])))

    def none_shifted(p):
        for a in h34.COUNTED:
            p[KEY1][0]["rows"][0][rec][a]["trend"]["shifted"] = []
    assert f"V6: {KEY1} seed 0: the arms shift different lineages, or none" in bad(none_shifted)
    assert any("placebo_s1.5 leaves a gene on its own factor" in w for w in
               bad(lambda p: p[KEY1][0]["rows"][0][rec]["placebo_s1.5"]["trend"]["fixed_points"].update(CM=2)))
    # a malformed part is void, never a crash
    for edit in (lambda p: p[KEY1][0].pop("v7"), lambda p: p[KEY1][0].pop("provenance"),
                 lambda p: p[KEY2][0]["rows"][0]["records"]["trend_s1.5"].pop("trend"),
                 lambda p: p[KEY2][0].update(rows=[{"no": "seed"}])):
        w = bad(edit)
        assert len(w) == 1 and w[0].startswith("a part or a committed record cannot be read ("), w
    assert why(committed={"row": "H29b", "void": False, "outcome": "READ"})[0].startswith(
        "a part or a committed record cannot be read (KeyError")


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
        assert txt.startswith(res["line"] + "\n" + NOTHING_BUILT + "\n") and "STOPPED" not in txt
        assert all(w in txt.splitlines() for w in res["why"]) and res["prior"] in txt
    res = _read(h34, score=stop, gate=_gate(h34, "FAIL"))
    assert res["line"] == "Try H34: VOID: G: the gate said 'FAIL', not PASS; nothing is built."


# ------------------------------------------------------------------------------------------------ the gate's verdict
def _gate_rows(h34):
    rows = {}
    for k in h34.TASK_KEYS:
        rows[k] = []
        for s in range(10):
            sc = _row(h34, k, s)["scores"]
            rows[k].append({"seed": s, "scores": {q: {a: sc[q][a] for a in ("floor", "control")} for q in h34.TRUTHS},
                            "size_record": dict(SIZE_REC), "draw": dict(DRAW),
                            "factors": {"ok": True, "why": None, "max_abs_dev": 2e-16}})
    return rows


def test_the_gates_two_checks(h34):
    def facts():
        return {k: _facts(h34) for k in h34.TASK_KEYS}
    why, dist = h34.gate_reasons(_gate_rows(h34), facts(), _committed(h34))
    assert why == [] and dist["G1"] == {"means": 18, "unequal": 0, "largest_abs_diff": 0.0}
    assert dist["G2"] == {"draws": 20, "largest_abs_dev": 2e-16, "tol": 1e-12}
    # G1, bit for bit: one seed's control score a part in 1e12 off
    rows = _gate_rows(h34)
    rows[KEY2][3]["scores"]["thin_0.5"]["control"] *= 1 + 1e-12
    why, dist = h34.gate_reasons(rows, facts(), _committed(h34))
    assert len(why) == 4 and all(w.startswith("G1: ") and "from row H29b's committed D_f" in w for w in why)
    assert dist["G1"]["unequal"] == 4 and 0 < dist["G1"]["largest_abs_diff"] < 1e-10
    # G1: a committed mean a part in 1e12 off
    c = _committed(h34)
    c["tables"]["truth radial_0"][KEY1]["mean"]["organ:f"] *= 1 + 1e-12
    why, _ = h34.gate_reasons(_gate_rows(h34), facts(), c)
    assert len(why) == 1 and why[0].startswith(f"G1: truth radial_0 {KEY1}: the control's mean is ")
    # G1: the target-free records and the control's size and draw
    f = facts()
    f[KEY1]["ratios"]["organ"] = 1.5 + 1e-12
    assert h34.gate_reasons(_gate_rows(h34), f, _committed(h34))[0] == [
        f"G1: {KEY1}: rms(pool) / rms(core) differs from row H29b's"]
    f = facts()
    f[KEY2]["pred_info"]["organ"]["core_share"] = 0.71
    assert h34.gate_reasons(_gate_rows(h34), f, _committed(h34))[0] == [
        f"G1: {KEY2}: D's pool differs from row H29b's (rows or core share)"]
    rows = _gate_rows(h34)
    rows[KEY1][5]["size_record"]["log_rms_after_vs_target"] = 1e-6
    assert h34.gate_reasons(rows, facts(), _committed(h34))[0] == [
        f"G1: {KEY1} seed 5: the control is not at D_f's size"]
    rows = _gate_rows(h34)
    rows[KEY1][5]["draw"]["short"] = True
    assert h34.gate_reasons(rows, facts(), _committed(h34))[0] == [
        f"G1: {KEY1}: the control's draw from D's pool differs from row H29b's"]
    # G2: one draw whose factors are not 1.5 times strength 1's
    rows = _gate_rows(h34)
    rows[KEY1][7]["factors"] = {"ok": False, "why": "a log factor is 1.000e-03 from 1.5 times the lower strength's",
                                "max_abs_dev": 1e-3}
    why, dist = h34.gate_reasons(rows, facts(), _committed(h34))
    assert why == [f"G2: {KEY1} seed 7: the strength-1.5 log factors are not 1.5 times the strength-1 ones (a log "
                   "factor is 1.000e-03 from 1.5 times the lower strength's)"]
    assert dist["G2"]["largest_abs_dev"] == 1e-3 and dist["G1"]["unequal"] == 0
    rows = _gate_rows(h34)
    rows[KEY2][0]["factors"] = {}
    assert h34.gate_reasons(rows, facts(), _committed(h34))[0][0].startswith(f"G2: {KEY2} seed 0: ")
    # both transitions and every seed, or the gate fails
    rows = _gate_rows(h34)
    rows.pop(KEY2)
    assert h34.gate_reasons(rows, facts(), _committed(h34))[0] == [f"G1: missing transition {KEY2}"]
    rows = _gate_rows(h34)
    rows[KEY1].pop()
    assert h34.gate_reasons(rows, facts(), _committed(h34))[0][0].startswith(f"G1: {KEY1}: seeds [0, 1, 2")
    for void in (None, {"row": "H29b", "void": True}):
        assert h34.gate_reasons(_gate_rows(h34), facts(), void)[0][0] == \
            "G1: row H29b's committed reading is void or unreadable"


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
    txt = h34.gate_text({**_gate(h34), "distances": {"G1": {"means": 18, "unequal": 0, "largest_abs_diff": 0.0},
                                                     "G2": {"draws": 20, "largest_abs_dev": 3e-16, "tol": 1e-12}}})
    assert "not a try" in txt and "no arm is built" in txt and txt.rstrip().endswith(h34.GATE_LINES["PASS"])
    assert "18 means compared, 0 unequal" in txt and "over 20 draws: largest difference 3.000e-16" in txt
    failed = h34.gate_text({**_gate(h34, "FAIL"), "why": ["G1: a", "G2: b"]}).splitlines()
    assert failed[-4:] == [h34.GATE_LINES["FAIL"], "WHY: G1: a; G2: b",
                           "Try H34: VOID: G1: a; G2: b; nothing is built.", NOTHING_BUILT]


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
    reading.write_text(json.dumps(committed, default=float))
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
    assert g["distances"]["G2"]["draws"] == 4 and g["distances"]["G2"]["largest_abs_dev"] <= 1e-12
    assert set(g["rows"]) == set(h34.TASK_KEYS) and [r["seed"] for r in g["rows"][KEY1]] == [0, 1]
    assert set(g["rows"][KEY1][0]["scores"]["own"]) == {"floor", "control"}  # only the control is built and scored
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
    assert res["void"] is False and res["why"] == [] and res["guards"]["V4"]["unequal"] == 0
    # two seeds cannot reach the bar's eight: clause (b) fires here whatever the scores, and the read says so
    assert res["outcome"] == "STOPPED" and res["clauses"] == ["b"]
    assert capsys.readouterr().out == "Try H34: STOPPED (b): OTS15 is not built.\n" + NOTHING_BUILT + "\n"
    assert (out / "reading.txt").read_text().startswith("Try H34: STOPPED (b): OTS15 is not built.\n")
    assert res["trend_record"][KEY1]["ref_cells"] == committed["pool_info"][KEY1]["organ"]["pool_cells"]
    assert all(v > 0 for v in res["changed_counts"][KEY2].values())
    # a committed mean a part in 1e12 off: the gate fails, exits 3 and prints the void line in its fixed words
    committed["tables"]["CLEAN"]["fold E10.5"]["mean"]["organ:f"] *= 1 + 1e-12
    reading.write_text(json.dumps(committed, default=float))
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
    # and a read against that gate is void
    h34.main([*read[:4], str(tmp_path / "gate2"), *read[5:8], str(tmp_path / "out2"), *read[9:]])
    assert capsys.readouterr().out.startswith("Try H34: VOID: G: the gate said 'FAIL', not PASS; ")
    # a missing committed reading fails the gate too
    with pytest.raises(SystemExit) as e:
        h34.main(["gate", "--out", str(tmp_path / "gate3"), "--h29b", str(tmp_path / "none.json"), *common])
    assert e.value.code == 3
    assert json.loads((tmp_path / "gate3" / "gate.json").read_text())["why"] == [
        "G1: row H29b's committed reading is void or unreadable"]


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
    for words in (*h34.FOLLOWS.values(), h34.SCREEN, h34.NO_PRIOR, *h34.GATE_LINES.values()):
        assert not re.search(r"\d", re.sub(r"H29b|H34|R32|strength-1\.5|1\.5 times|strength-1 ", "", words)), words


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


@pytest.fixture()
def repo(tmp_path):
    """A throwaway git repository with a ledger that holds a row H34 and a stand-in for the script, a temporary home
    and the stand-in uv."""
    root = tmp_path / "repo"
    (root / "experiments" / "harness").mkdir(parents=True)
    (tmp_path / "home").mkdir()
    (tmp_path / "bin").mkdir()
    (tmp_path / "bin" / "uv").write_text(STUB_UV)
    (tmp_path / "bin" / "uv").chmod(0o755)
    (root / "experiments" / "harness" / "TRIES.md").write_text("| H29b (a row) | heart | x |\n"
                                                              "| H34 (the row) | heart | x |\n")
    (root / "experiments" / "harness" / "h34_heart_organ_trend.py").write_text("# a stand-in\n")
    (root / ".gitignore").write_text("runs/\n")

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
    assert 'export UV_PROJECT_ENVIRONMENT="${H34_VENV:-$PWD/.venv}"' in code
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
