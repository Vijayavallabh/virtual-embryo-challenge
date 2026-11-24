"""Row H34 (experiments/harness/h34_heart_organ_trend.py) and its launcher (experiments/harness/run_h34.sh): synthetic
arrays, hand-made tables and temporary files, and committed result files of rows H27 and H29b read for arithmetic only
(the gate's dry run). Nothing here loads a competition file, scores a real cell or reads the repository's own ledger;
the launcher runs in throwaway git repositories under the test's own folder, with a temporary home and a stand-in uv
first on PATH."""

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
    assert h34.FLOOR_PAIRS == tuple((a, "floor") for a in h34.BUILT) and len(h34.PAIRS) == 12
    assert [w for w, _ in h34.G1_WORLDS] == ["own", "CLEAN", "truth own", "truth radial_0", "truth thin_0.5",
                                              "truth thin_0"]
    assert h34.G1_UNITS == (KEY1, KEY2, "fold E10.5") and set(h34.UNITS) == set(h34.G1_UNITS)
    assert all(u in h34.h29b.UNITS for u in h34.G1_UNITS)
    assert h34.G1_ARMS == ("core_iso_half:f", "half_iso_half:f", "full_iso_half:f")  # C_f, B_f, A_f
    assert h34.STAT_FIELDS == ("mean", "se", "lo", "hi", "ahead") and h34.TARGET_FIELDS == ("target", "base", "ratio")
    assert len(h34.G1_ARMS) * len(h34.G1_WORLDS) * len(h34.G1_UNITS) * len(h34.STAT_FIELDS) == 270
    assert h34.RECORD_PROSE == ("rule", "label_config")
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
                 "h29b.d_targets(", "h29b.size_to(", "h25.kept_rows(", "h29b._checked(", "h29b.series(",
                 "h29b._stats(", "h29b.load_h27("):
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
        assert m_lo[k].dtype == m_hi[k].dtype == logf[k].dtype == np.float64
    assert h34.factor_identity(lo, hi) == {"ratio": 1.5, "lineages": ["CM", "Endothelial"], "damping": [0.5, 0.75],
                                           "ok": True, "why": None}
    assert h34.factor_identity(lo, lo)["ok"] is False  # strength 1 against itself is not the identity of 1.5
    # the arguments are taken at the call, as they are handed over: a copy of each vector in its own dtype
    handed = {"CM": np.float32([0.1, 0.2, 0.3])}

    def hands(ctx, s, placebo_arm=False):
        t1_scope.momentum_shift(None, None, None, 0.5 * s, (), key="b0_type", seed=[0, 3], logf_by_type=handed)
    monkeypatch.setattr(h34.h28, "trend_expression", hands)
    d, m = h34.step_arguments(dform.ctx(), 1.5)
    assert d == 0.75 and m["CM"].dtype == np.float32 and m["CM"] is not handed["CM"]
    np.testing.assert_array_equal(m["CM"], handed["CM"])
    monkeypatch.setattr(h34.h28, "trend_expression", lambda ctx, s, placebo_arm=False: None)
    with pytest.raises(RuntimeError, match="did not reach momentum_shift exactly once"):
        h34.step_arguments(dform.ctx(), 1.0)
    assert t1_scope.momentum_shift is own


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
    # every element finite: a NaN, and an infinity that 1.5 times itself would equal, both fail; so does a damping
    for bad_value in (np.nan, np.inf, -np.inf):
        (d_lo, m_lo), (d_hi, m_hi) = at(1.0), at(1.5)
        m_lo["Endothelial"][0], m_hi["Endothelial"][0] = bad_value, bad_value
        bad = h34.factor_identity((d_lo, m_lo), (d_hi, m_hi))
        assert bad["ok"] is False and bad["why"] == "a log factor or a damping is not finite", bad_value
        assert h34.factor_identity((bad_value, at(1.0)[1]), (bad_value, at(1.5)[1]))["ok"] is False
    # no lineage in the mapping: the damping alone is held
    assert h34.factor_identity((0.5, {}), (0.75, {}))["ok"] is True
    assert h34.factor_identity((0.5, {}), (0.5, {}))["ok"] is False
    # another ratio
    assert h34.factor_identity(at(1.0), at(2.0), 2.0)["ok"] is True
    assert h34.factor_identity(at(1.0), at(2.0))["ok"] is False
    # the identity is held in the vectors' own arithmetic: single-precision vectors scaled by the step's one product
    # pass, though that product is not the double-precision one
    m32 = {k: v.astype(np.float32) for k, v in m.items()}
    lo32, hi32 = (0.5, {k: 1.0 * v for k, v in m32.items()}), (0.75, {k: 1.5 * v for k, v in m32.items()})
    assert all(v.dtype == np.float32 for v in hi32[1].values()) and h34.factor_identity(lo32, hi32)["ok"] is True
    assert not np.array_equal(hi32[1]["CM"].astype(np.float64), 1.5 * m32["CM"].astype(np.float64))


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
G1_BASE = {"core_iso_half:f": 49.0, "half_iso_half:f": 49.4, "full_iso_half:f": 49.9}  # row H27's three controls
FACTORS = {"ratio": 1.5, "lineages": ["CM", "Endothelial"], "damping": [0.5, 0.75], "ok": True, "why": None}


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


def _gate_rows(h34, seeds=range(10)):
    """The gate's rows as hand-made: per transition and seed the control's and the floor's scores, its size record,
    its draw and G2's finding."""
    rows = {}
    for k in h34.TASK_KEYS:
        rows[k] = []
        for s in seeds:
            sc = _row(h34, k, s)["scores"]
            rows[k].append({"seed": s, "scores": {q: {a: sc[q][a] for a in ("floor", "control")} for q in h34.TRUTHS},
                            "size_record": dict(SIZE_REC), "draw": dict(DRAW), "factors": dict(FACTORS)})
    return rows


def _h27_rows(h34, seeds=range(10)):
    """Row H27's committed rows as hand-made: its three controls' scores per transition, seed and truth."""
    return {k: [{"seed": s, "scores": {q: {a: G1_BASE[a] + 0.3 * h34.TASK_KEYS.index(k) + 0.02 * ((7 * s) % 10)
                                           + 0.003 * h34.TRUTHS.index(q) for a in h34.G1_ARMS} for q in h34.TRUTHS}}
                for s in seeds] for k in h34.TASK_KEYS}


def _committed(h34, seeds=range(10)):
    """Row H29b's committed reading as the hand-made control and row H27's hand-made rows make it, by row H29b's own
    functions: D_f's plain means, its paired statistics against the three controls, and its target-free record."""
    h29b = h34.h29b
    rows, others = _gate_rows(h34, seeds), _h27_rows(h34, seeds)
    f = _facts(h34)
    tables = {}
    for w, ts in h34.G1_WORLDS:
        tables[w] = {}
        for u in h34.G1_UNITS:
            c = h29b.series(rows, h29b.UNITS[u], ts, "control")
            tables[w][u] = {"mean": {"organ:f": float(c.mean())},
                            "stated": {f"organ:f - {a}": h29b._stats(c - h29b.series(others, h29b.UNITS[u], ts, a))
                                       for a in h34.G1_ARMS}}

    def per(v):
        return {k: copy.deepcopy(v) for k in h34.TASK_KEYS}
    short = {"organ": {"pool_cells": POOL_CELLS, "drawn": [DRAWN], "short_in_seeds": 0, "repeated_rows_at_least": 0}}
    return {"row": "H29b", "void": False, "outcome": "READ", "tables": tables, "truth_info": per(f["truth_info"]),
            "pool_info": per(f["pred_info"]), "organ_records": per(f["organ_record"]), "ratios": per(f["ratios"]),
            "d_targets": per(f["d_targets"]), "pools_short": per(short)}


def _gate(h34, verdict="PASS", commit=COMMIT, clean=True, seeds=range(10)):
    return {"row": "H34", "axis": h34.GATE_AXIS, "verdict": verdict, "commit": commit,
            "why": [] if verdict == "PASS" else ["G1 (M): x"], "found": {}, "rows": _gate_rows(h34, seeds),
            "transitions": {k: _facts(h34) for k in h34.TASK_KEYS},
            "provenance": {"commit": commit, "tree_clean": clean}}


def _read(h34, parts=None, gate="default", commit=COMMIT, unreadable=(), **kw):
    return h34.reading(_parts(h34, **kw) if parts is None else parts, _gate(h34) if gate == "default" else gate,
                       commit, unreadable)


NOT_STOPPED_LINES = ("Try H34: NOT STOPPED.",
                     "A screen, not evidence of development. What follows: step (3), the builder option under ruling "
                     "row R32, is coded and reviewed, then the identity builds (step 4), then the file (step 5), each "
                     "as the line's Build passage has it; NOT STOPPED uploads nothing, the file going in only under "
                     "the line's Rule and on the captain's word at the time.")
NOTHING_BUILT = ("A screen, not evidence of development. What follows: nothing is built; the line ends there, with no "
                 "second read of OTS15's form under this or another row.")
GOOD = {HI: 0.6, LO: 0.4}  # 1.5 above 1 and ahead of the copy by the bar, in every seed and both worlds


# ------------------------------------------------------------------------------------------------ the rule
def test_not_stopped_when_1_5_is_above_1_and_ahead_of_the_copy_in_both_worlds(h34):
    res = _read(h34, score=GOOD)
    assert res["why"] == [] and res["void"] is False and res["outcome"] == "NOT STOPPED" and res["clauses"] == []
    st = res["statements"]
    assert st["a_in_world"] == {"own": False, "CLEAN": False} and st["bar_met"] == {"own": True, "CLEAN": True}
    assert st["holds"] == {"a": False, "b": False} and st["stop"] is False
    assert (res["line"], res["follows"]) == NOT_STOPPED_LINES
    for w in ("own", "CLEAN"):
        a = res["tables"][w]["fold E10.5"]["trend_s1.5 - trend_s1.0"]["score"]
        b = res["tables"][w]["fold E10.5"]["trend_s1.5 - control"]["score"]
        assert a["mean"] == pytest.approx(0.2) and a["ahead"] == 10 and a["behind"] == 0
        assert b["mean"] == pytest.approx(0.6) and b["ahead"] == 10 and b["behind"] == 0
    assert set(res["tables"]["own"]) == {KEY1, KEY2, "fold E10.5"}
    assert set(res["tables"]["own"][KEY1]) == {h34.pname(a, b) for a, b in h34.PAIRS}


def _dips(n_below, rise=0.01):
    """trend_s1.5 minus trend_s1.0 per seed: 0.05 below zero in n_below seeds, `rise` above it in the others."""
    return np.where(np.arange(10) < n_below, -0.05, rise)


def test_clause_a_a_mean_below_zero_with_eight_seeds_below_in_either_world(h34):
    def run(own, clean):
        # trend_s1.5 stays 0.6 over the copy in both worlds, so clause (b) is silent; trend_s1.0 is it minus d
        return _read(h34, score={HI: 0.6, LO: {"own": 0.6 - np.asarray(own), "CLEAN": 0.6 - np.asarray(clean)}})
    up = np.full(10, 0.2)
    res = run(_dips(8), up)
    assert res["outcome"] == "STOPPED" and res["clauses"] == ["a in own"] and res["void"] is False
    assert res["statements"]["a_in_world"] == {"own": True, "CLEAN": False}
    assert res["statements"]["bar_met"] == {"own": True, "CLEAN": True}
    assert res["statements"]["holds"] == {"a": True, "b": False}
    # the statement names every clause that holds and, for clause (a), its world or worlds
    assert res["line"] == "Try H34: STOPPED (a in own): OTS15 is not built." and res["follows"] == NOTHING_BUILT
    s = res["tables"]["own"]["fold E10.5"]["trend_s1.5 - trend_s1.0"]["score"]
    assert s["mean"] == pytest.approx(-0.038) and s["behind"] == 8 and s["ahead"] == 2
    res = run(up, _dips(8))  # in own or in CLEAN: CLEAN alone stops
    assert res["line"] == "Try H34: STOPPED (a in CLEAN): OTS15 is not built."
    assert res["statements"]["a_in_world"] == {"own": False, "CLEAN": True}
    assert run(_dips(10), _dips(10))["line"] == "Try H34: STOPPED (a in own and CLEAN): OTS15 is not built."
    assert run(_dips(9), up)["clauses"] == ["a in own"]
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
    assert res["clauses"] == ["a in own and CLEAN", "b"]
    assert res["line"] == "Try H34: STOPPED (a in own and CLEAN, b): OTS15 is not built."
    assert res["statements"]["holds"] == {"a": True, "b": True}
    res = _read(h34, score={HI: {"own": 0.1, "CLEAN": 0.6}, LO: {"own": 0.2, "CLEAN": 0.4}})
    assert res["line"] == "Try H34: STOPPED (a in own, b): OTS15 is not built."
    # arm minus placebo, the DE terms and the floor: whatever they are, the outcome is as it was
    for plc, de in ((2.0, -0.3), (-2.0, 0.3), (0.6, 0.0)):
        extra = {"placebo_s1.5": plc, "placebo_s1.0": plc}
        ok = _read(h34, score={**GOOD, **extra}, de={HI: de, LO: de, "placebo_s1.5": -de})
        assert ok["outcome"] == "NOT STOPPED"
        d = ok["tables"]["own"]["fold E10.5"]
        assert d["trend_s1.5 - placebo_s1.5"]["score"]["mean"] == pytest.approx(0.6 - plc)
        assert d["trend_s1.5 - control"]["de_score"] == pytest.approx(de)
        assert d["trend_s1.5 - control"]["d2_shape"] == 0.0
        assert _read(h34, score={HI: 0.1, LO: 0.2, **extra}, de={HI: de, LO: de})["clauses"] == [
            "a in own and CLEAN", "b"]
    # every arm below the floor, yet 1.5 above 1 and ahead of the copy: the floor gates nothing
    low = _read(h34, score=GOOD)
    d = low["tables"]["own"]["fold E10.5"]
    assert d["control - floor"]["score"]["mean"] == pytest.approx(0.25 + 0.045)
    assert d["trend_s1.5 - floor"]["score"]["mean"] == pytest.approx(0.25 + 0.045 + 0.6)


def test_the_rules_two_bars_at_their_edges(h34):
    assert h34.below({"mean": -1e-9, "behind": 8}) and h34.below({"mean": -3.0, "behind": 10})
    assert not h34.below({"mean": 0.0, "behind": 10}) and not h34.below({"mean": -1.0, "behind": 7})
    assert not h34.below({"mean": 1e-9, "behind": 9})
    passes = h34.h28.passes  # row H28's S1 bar, in one world; a paired mean of exactly +0.3 is not above +0.3
    assert passes({"mean": 0.3000001, "ahead": 8}) and not passes({"mean": 0.3, "ahead": 10})
    assert not passes({"mean": 5.0, "ahead": 7})
    assert h34.outcome_line("VOID", ["v3: no gate record", "v4: missing transition x"]) == \
        "Try H34: VOID: v3: no gate record; v4: missing transition x; nothing is built."
    assert h34.outcome_line("STOPPED", clauses=["a in own", "b"]) == \
        "Try H34: STOPPED (a in own, b): OTS15 is not built."
    assert h34.outcome_line("NOT STOPPED") == "Try H34: NOT STOPPED."
    assert h34.follows_line("VOID") == h34.follows_line("STOPPED") == NOTHING_BUILT
    assert h34.follows_line("NOT STOPPED") == NOT_STOPPED_LINES[1]
    assert h34.gate_void(["G1 (M): a", "G2: b"]) == ["v3: the gate fails (G1 (M): a; G2: b)"]


def test_the_control_in_the_rule_is_the_gates(h34):
    """The control's scores in the rule are the gate's, the ones G1 checked; everything else is the run step's."""
    parts, gate = _parts(h34, score=GOOD), _gate(h34)
    for r in gate["rows"][KEY1]:
        r["scores"]["own"]["control"] += 0.125  # (a read with such a gate is void; the function is shown alone)
    rows = h34.rule_rows(parts, gate)
    mine = {k: h34.h26.merge_rows(parts[k])[0] for k in h34.TASK_KEYS}
    for k in h34.TASK_KEYS:
        assert [r["seed"] for r in rows[k]] == list(range(10))
        for r, m, g in zip(rows[k], mine[k], gate["rows"][k], strict=True):
            for q in h34.TRUTHS:
                assert r["scores"][q]["control"] == g["scores"][q]["control"]
                assert {a: v for a, v in r["scores"][q].items() if a != "control"} == \
                    {a: v for a, v in m["scores"][q].items() if a != "control"}
            assert r["skills"] == m["skills"] and r["records"] == m["records"]
    assert rows[KEY1][0]["scores"]["own"]["control"] == mine[KEY1][0]["scores"]["own"]["control"] + 0.125
    assert _read(h34, parts=parts, gate=gate)["outcome"] == "VOID"
    assert _read(h34, score=GOOD)["control_scores"] == "the gate's, the ones G1 checked"


def test_the_reading_and_its_text(h34):
    res = _read(h34, score=GOOD)
    json.dumps(res, default=float)
    assert res["row"] == "H34" and res["commit"] == COMMIT and res["gate"]["verdict"] == "PASS"
    assert res["prior"] == "No row on record had read the trend step on D's pool at any strength before this one."
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
                 "TABLES, practice diagnostics that gate nothing", "eight tries"):
        assert need in txt, need
    for name in h34.BUILT:  # the control and each arm minus the floor, the published baseline
        assert f"DIAGNOSTIC, gates nothing: {name} minus the floor (the published baseline), fold E10.5: own " in txt
    tagged = [ln for ln in lines if ln.startswith(h34.TAG)]
    assert len(tagged) == len(h34.WORLDS) * len(h34.UNITS) * len(h34.PAIRS) == 72
    stopped = h34.text(_read(h34, score={HI: 0.1, LO: 0.2}))
    assert stopped.startswith("Try H34: STOPPED (a in own and CLEAN, b): OTS15 is not built.\n" + NOTHING_BUILT
                              + "\n")
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
                  "The control's scores in the rule are the gate's, the ones G1 checked.",
                  "VOID is decided first, and under VOID the rule is not read.",
                  "a difference of exactly 0 is neither above nor below 0, and a paired mean of exactly +0.3 is not "
                  "above +0.3"):
        assert words in h34.RULE, words


# ------------------------------------------------------------------------------------------------ the void
def test_every_void_condition(h34):
    def why(**kw):
        return _read(h34, **{"score": GOOD, **kw})["why"]

    def bad(edit, **kw):
        p = _parts(h34, score=GOOD)
        edit(p)
        return why(parts=p, **kw)

    def gated(edit):
        g = _gate(h34)
        edit(g)
        return why(gate=g)
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
    # v4: a score that is not a finite number, of an arm, of the control or of the floor
    for value in (float("nan"), float("inf"), float("-inf"), None, "50.0", True):
        for name in ("trend_s1.5", "control", "floor"):
            assert bad(lambda p, value=value, name=name: p[KEY1][0]["rows"][2]["scores"]["own"].update(
                {name: value})) == [f"v4: {KEY1} seed 2: a score is not a finite number: ['own/{name}']"], (value, name)
    assert bad(lambda p: p[KEY2][1]["rows"][0]["scores"]["thin_0"].update(floor=50)) == []  # a whole number is one
    # v4: a part that cannot be read or is of another transition
    assert why(unreadable=["/x/e95/0-1/h34.json"]) == ["v4: a part cannot be read: /x/e95/0-1/h34.json"]
    assert "v4: parts of other transitions ['T2h/E8.25_late->E8.5']" in bad(
        lambda p: p.update({"T2h/E8.25_late->E8.5": p[KEY1]}))
    for edit in (lambda p: p[KEY1][0].pop("provenance"), lambda p: p[KEY1][0].pop("pred_info"),
                 lambda p: p[KEY2][0]["rows"][0]["records"]["trend_s1.5"].pop("trend"),
                 lambda p: p[KEY2][0].update(rows=[{"no": "seed"}])):
        w = bad(edit)
        assert len(w) == 1 and w[0].startswith("v4: a part cannot be read ("), w
    # v4: the parts come from more than one commit, from a commit other than the run commit or from an unclean tree
    assert "v4: the parts come from more than one commit" in bad(
        lambda p: p[KEY1][0]["provenance"].update(commit="d" * 40))
    assert why(parts=_parts(h34, clean=False, score=GOOD)) == ["v4: a part comes from an unclean tree"]
    assert bad(lambda p: p[KEY2][1]["provenance"].update(tree_clean=None)) == ["v4: a part comes from an unclean tree"]
    for other in ("d" * 40, COMMIT[:12], "", None):
        w = why(commit=other, gate=_gate(h34, commit=other))
        assert any(x.startswith("v4: a part comes from a commit other than the run commit") for x in w), other
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
    # v5: an arm's record does not hold its own strength (1.0 or 1.5) and kind (trend or placebo)
    assert bad(lambda p: p[KEY1][0]["rows"][0][rec]["trend_s1.5"]["trend"].update(strength=1.0)) == [
        f"v5: {KEY1} seed 0: trend_s1.5's record does not hold strength 1.5 and kind trend"]
    assert bad(lambda p: p[KEY1][0]["rows"][0][rec]["placebo_s1.0"]["trend"].update(placebo=False)) == [
        f"v5: {KEY1} seed 0: placebo_s1.0's record does not hold strength 1.0 and kind placebo"]
    assert bad(lambda p: p[KEY1][0]["rows"][0][rec]["trend_s1.0"]["trend"].update(placebo=True)) == [
        f"v5: {KEY1} seed 0: trend_s1.0's record does not hold strength 1.0 and kind trend"]
    for gone in ("strength", "placebo"):
        assert bad(lambda p, gone=gone: p[KEY2][0]["rows"][1][rec]["placebo_s1.5"]["trend"].pop(gone)) == [
            f"v5: {KEY2} seed 1: placebo_s1.5's record does not hold strength 1.5 and kind placebo"], gone
    assert bad(lambda p: p[KEY1][0]["rows"][0][rec]["trend_s1.0"]["trend"].update(budget_days=0.5)) == [
        f"v5: {KEY1} seed 0: trend_s1.0's budget is 0.5, not 0.375"]
    assert bad(lambda p: p[KEY2][0]["rows"][0][rec]["placebo_s1.5"]["trend"].update(
        budget_days=0.125 * (1 + 1e-11))) != []  # relative 1e-12
    assert bad(lambda p: p[KEY2][0]["rows"][0][rec]["placebo_s1.5"]["trend"].update(
        budget_days=0.125 * (1 + 4e-13))) == []
    assert bad(lambda p: p[KEY1][0]["rows"][0][rec]["placebo_s1.5"]["trend"].update(ref_cells=POOL_CELLS + 5)) == [
        f"v5: {KEY1} seed 0: placebo_s1.5's reference is not D's pool"]
    assert bad(lambda p: p[KEY1][0]["rows"][0][rec]["trend_s1.5"]["trend"].update(shifted=["CM"])) == [
        f"v5: {KEY1} seed 0: the four arms do not shift the same lineages"]
    assert bad(lambda p: p[KEY1][0]["rows"][0][rec]["placebo_s1.5"]["trend"]["fixed_points"].update(CM=2)) == [
        f"v5: {KEY1} seed 0: placebo_s1.5 leaves a gene on its own factor"]
    # v5: the control or the floor as the run step scores them is not exactly the gate's
    assert bad(lambda p: p[KEY2][1]["rows"][2]["scores"]["thin_0"].update(
        control=p[KEY2][1]["rows"][2]["scores"]["thin_0"]["control"] * (1 + 1e-12))) == [
        f"v5: {KEY2} seed 7 thin_0: the run step's control score is not the gate's"]
    assert bad(lambda p: p[KEY1][0]["rows"][3]["scores"]["own"].update(floor=50.0 + 1e-12)) == [
        f"v5: {KEY1} seed 3 own: the run step's floor score is not the gate's"]
    assert gated(lambda g: g["rows"][KEY2][4]["scores"]["radial_0"].update(control=1.0)) == [
        f"v5: {KEY2} seed 4 radial_0: the run step's control score is not the gate's"]
    w = gated(lambda g: g["rows"][KEY1].pop())
    assert len(w) == 1 and w[0].startswith("v3: the gate record holds no score to hold the run step's to (KeyError")
    w = gated(lambda g: g.pop("rows"))
    assert len(w) == 1 and w[0].startswith("v3: the gate record holds no score to hold the run step's to (KeyError")


def test_void_wins_over_stop(h34):
    stop = {HI: 0.1, LO: 0.2}
    assert _read(h34, score=stop)["outcome"] == "STOPPED"
    dirty = _parts(h34, score=stop)
    dirty[KEY1][0]["provenance"]["tree_clean"] = False
    nudged = _parts(h34, score=stop)
    nudged[KEY1][0]["rows"][0]["scores"]["own"]["control"] += 1e-9
    for kw in ({"gate": _gate(h34, "FAIL")}, {"gate": None}, {"commit": "d" * 40}, {"unreadable": ["x/h34.json"]},
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
    # every void condition found is in the statement
    res = _read(h34, parts=dirty, gate=None, score=stop)
    assert res["line"] == ("Try H34: VOID: v3: no gate record; v4: a part comes from an unclean tree; nothing is "
                           "built.")


# ------------------------------------------------------------------------------------------------ the gate's verdict
def _gate_reasons(h34, rows=None, facts=None, committed="default", h27="default"):
    return h34.gate_reasons(_gate_rows(h34) if rows is None else rows,
                            {k: _facts(h34) for k in h34.TASK_KEYS} if facts is None else facts,
                            _committed(h34) if committed == "default" else committed,
                            _h27_rows(h34) if h27 == "default" else h27)


def test_the_gate_passes_only_if_g1_and_g2_both_hold(h34):
    why, found = _gate_reasons(h34)
    assert why == []
    assert found == {"R": {KEY1: [], KEY2: []}, "D": {KEY1: True, KEY2: True},
                     "M": {"means": 18, "unequal": [], "largest_abs_diff": 0.0},
                     "P": {"numbers": 270, "unequal": []},
                     "G2": {"contexts": 20, "failed": 0, "lineages": ["CM", "Endothelial"]}}
    # G2: one seed's context where the identity fails; G1 holds
    rows = _gate_rows(h34)
    rows[KEY1][7]["factors"] = {"ok": False, "lineages": ["CM"],
                                "why": "the damping 0.5 is not 1.5 times the lower strength's 0.5"}
    rows[KEY2][0]["factors"] = {}
    why, found = _gate_reasons(h34, rows=rows)
    assert why == [f"G2: 2 of 20 contexts fail (first: {KEY1} seed 7: the damping 0.5 is not 1.5 times the lower "
                   "strength's 0.5)"]
    assert found["G2"] == {"contexts": 20, "failed": 2, "lineages": ["CM", "Endothelial"]}
    assert found["M"]["unequal"] == [] and found["P"]["unequal"] == []
    # both fail: both are named
    rows[KEY1][0]["scores"]["own"]["control"] += 1e-9
    why, _ = _gate_reasons(h34, rows=rows)
    assert [w.split(":")[0] for w in why] == ["G1 (M)", "G1 (P)", "G2"]
    # the control on both transitions and every seed, the reading as the row names it, row H27's parts
    rows = _gate_rows(h34)
    rows.pop(KEY2)
    assert _gate_reasons(h34, rows=rows)[0] == [
        "G1: the control was not built and scored on both transitions and every seed"]
    rows = _gate_rows(h34)
    rows[KEY1].pop()
    assert _gate_reasons(h34, rows=rows)[0] == [
        "G1: the control was not built and scored on both transitions and every seed"]
    for void in (None, {"row": "H29b", "void": True}, {"row": "H27", "void": False, "outcome": "READ"}):
        why, found = _gate_reasons(h34, committed=void)
        assert why == ["G1: row H29b's reading.json is not at hand as the file its row records (missing, another "
                       "sha256, or void)"]
        assert found["M"] is None and found["G2"]["failed"] == 0
    why, found = _gate_reasons(h34, h27=None)
    assert why == ["G1 (P): row H27's committed parts are not at hand"] and found["P"] is None
    assert found["M"]["unequal"] == []


def test_g1_r_the_committed_record_of_d_f(h34):
    """(R): ratios.organ, the target, base and ratio of d_targets, pool_info.organ, and organ_records.reference but
    its two descriptions; each by exact equality."""
    def off(edit):
        facts = {k: _facts(h34) for k in h34.TASK_KEYS}
        edit(facts)
        why, found = _gate_reasons(h34, facts=facts)
        assert found["M"]["unequal"] == [] and found["D"] == {KEY1: True, KEY2: True}
        return why
    assert off(lambda f: f[KEY1]["ratios"].update(organ=1.5 + 1e-12)) == [
        f"G1 (R): {KEY1}: ratios.organ not row H29b's committed"]
    for field in ("target", "base", "ratio"):
        assert off(lambda f, field=field: f[KEY2]["d_targets"]["organ:f"].update(
            {field: f[KEY2]["d_targets"]["organ:f"][field] * (1 + 1e-12)})) == [
            f"G1 (R): {KEY2}: d_targets[organ:f].{field} not row H29b's committed"]
    for field, value in (("pool_cells", POOL_CELLS + 1), ("core_share", 0.7 + 1e-12), ("xy_scope", 0.1 + 1e-12)):
        assert off(lambda f, field=field, value=value: f[KEY1]["pred_info"]["organ"].update({field: value})) == [
            f"G1 (R): {KEY1}: pool_info.organ not row H29b's committed"], field
    for field, value in (("link_length", 2.5 + 1e-12), ("label_config_sha256", "1" * 64), ("components_kept", 3),
                         ("components_total", 6), ("kept_categories", ["cm"]), ("kept_share", 0.5 + 1e-12),
                         ("chamber_share", 0.61), ("core_share_of_kept", 0.71), ("link_percentile", 98.0),
                         ("a_new_field", 1)):
        assert off(lambda f, field=field, value=value: f[KEY2]["organ_record"]["reference"].update(
            {field: value})) == [f"G1 (R): {KEY2}: organ_records.reference not row H29b's committed"], field
    # its two descriptions are left out; so are the label pool's ratio, the truths and the size's wording
    assert off(lambda f: f[KEY1]["organ_record"]["reference"].update(rule="other words", label_config="x")) == []
    assert off(lambda f: f[KEY1]["ratios"].update(label=9.9)) == []
    assert off(lambda f: f[KEY1]["d_targets"]["organ:f"].update(source="other words")) == []
    assert off(lambda f: f[KEY1]["truth_info"]["own"].update(cells=1)) == []
    # several items on one transition are named together; a missing item cannot hold
    assert off(lambda f: (f[KEY1]["ratios"].update(organ=1.6), f[KEY1]["pred_info"]["organ"].update(
        pool_cells=1))) == [f"G1 (R): {KEY1}: ratios.organ, pool_info.organ not row H29b's committed"]
    committed = _committed(h34)
    committed["organ_records"][KEY1]["reference"].pop("link_length")
    assert _gate_reasons(h34, committed=committed)[0] == [
        f"G1 (R): {KEY1}: organ_records.reference not row H29b's committed"]


def test_g1_d_the_committed_draw_of_d_f(h34):
    """(D): pools_short.organ over the ten seeds, as row H29b's read forms it."""
    for field, value in (("pool_cells", POOL_CELLS + 1), ("drawn", DRAWN + 1), ("short", True),
                         ("repeated_rows_at_least", 3)):
        rows = _gate_rows(h34)
        rows[KEY1][0 if field == "pool_cells" else 6]["draw"][field] = value
        why, found = _gate_reasons(h34, rows=rows)
        assert why == [f"G1 (D): {KEY1}: the control's draw is not row H29b's committed pools_short.organ"], field
        assert found["D"] == {KEY1: False, KEY2: True}
    rows = _gate_rows(h34)
    assert h34.pool_draw(rows[KEY2]) == {"pool_cells": POOL_CELLS, "drawn": [DRAWN], "short_in_seeds": 0,
                                         "repeated_rows_at_least": 0}
    for s, r in enumerate(rows[KEY2]):  # as h29b.void_reasons forms it: the sorted counts, the seeds short, the most
        r["draw"] = {"pool_cells": 300, "drawn": 800 + s % 2, "short": s < 4, "repeated_rows_at_least": 500 + s}
    assert h34.pool_draw(rows[KEY2]) == {"pool_cells": 300, "drawn": [800, 801], "short_in_seeds": 4,
                                         "repeated_rows_at_least": 509}


def test_g1_m_and_p_the_committed_means_and_paired_statistics(h34):
    """(M) the eighteen means and (P) the 270 paired statistics, each by exact equality with no tolerance."""
    # one seed's control score a part in 1e12 off: the means and the paired statistics of its worlds and units
    rows = _gate_rows(h34)
    rows[KEY2][3]["scores"]["thin_0.5"]["control"] *= 1 + 1e-12
    why, found = _gate_reasons(h34, rows=rows)
    assert [(w, u) for w, u, _ in found["M"]["unequal"]] == [("CLEAN", KEY2), ("CLEAN", "fold E10.5"),
                                                             ("truth thin_0.5", KEY2),
                                                             ("truth thin_0.5", "fold E10.5")]
    assert 0 < found["M"]["largest_abs_diff"] < 1e-10
    assert why[0].startswith(f"G1 (M): 4 of 18 means are not row H29b's committed D_f means (first: CLEAN {KEY2}, ")
    assert len(why) == 2 and why[1].startswith(f"G1 (P): {len(found['P']['unequal'])} of 270 paired statistics are "
                                               f"not row H29b's committed (first: CLEAN {KEY2} organ:f - "
                                               "core_iso_half:f mean)")
    assert {(w, u) for w, u, _, _ in found["P"]["unequal"]} == {(w, u) for w, u, _ in found["M"]["unequal"]}
    assert {p for _, _, p, _ in found["P"]["unequal"]} == {f"organ:f - {a}" for a in h34.G1_ARMS}
    # each of the eighteen committed means, one unit in the last place off
    for w, _ in h34.G1_WORLDS:
        for u in h34.G1_UNITS:
            c = _committed(h34)
            c["tables"][w][u]["mean"]["organ:f"] = float(np.nextafter(c["tables"][w][u]["mean"]["organ:f"], 100.0))
            why, found = _gate_reasons(h34, committed=c)
            assert len(why) == 1 and why[0].startswith(f"G1 (M): 1 of 18 means are not row H29b's committed D_f "
                                                       f"means (first: {w} {u}, -")
            assert found["P"]["unequal"] == []
    # each of the 270 committed paired statistics, one unit in the last place (or one seed) off
    n = 0
    for w, _ in h34.G1_WORLDS:
        for u in h34.G1_UNITS:
            for a in h34.G1_ARMS:
                for k in h34.STAT_FIELDS:
                    c = _committed(h34)
                    s = c["tables"][w][u]["stated"][f"organ:f - {a}"]
                    s[k] = s[k] + 1 if k == "ahead" else float(np.nextafter(s[k], 100.0))
                    why, found = _gate_reasons(h34, committed=c)
                    assert found["P"]["unequal"] == [[w, u, f"organ:f - {a}", k]] and found["M"]["unequal"] == []
                    assert why == [f"G1 (P): 1 of 270 paired statistics are not row H29b's committed (first: {w} "
                                   f"{u} organ:f - {a} {k})"]
                    n += 1
    assert n == 270
    # two seeds exchanged on one transition: that transition's means are the same ten numbers; the paired
    # statistics, which pair by seed, are not
    rows = _gate_rows(h34)
    rows[KEY1][2]["scores"], rows[KEY1][5]["scores"] = rows[KEY1][5]["scores"], rows[KEY1][2]["scores"]
    why, found = _gate_reasons(h34, rows=rows)
    assert found["P"]["unequal"] and any(w.startswith("G1 (P): ") for w in why)
    # one of row H27's committed scores a part in 1e12 off: (P) alone
    others = _h27_rows(h34)
    others[KEY1][4]["scores"]["own"]["half_iso_half:f"] *= 1 + 1e-12
    why, found = _gate_reasons(h34, h27=others)
    assert found["M"]["unequal"] == [] and len(why) == 1 and why[0].startswith("G1 (P): ")
    assert {p for _, _, p, _ in found["P"]["unequal"]} == {"organ:f - half_iso_half:f"}
    # the figures of other arms and other worlds in row H29b's table do not enter
    c = _committed(h34)
    c["tables"]["own"][KEY1]["mean"]["core_iso_half:f"] = 1.0
    c["tables"]["own"][KEY1]["stated"]["label:f - organ:f"] = {"mean": 1.0}
    c["tables"]["truth organ"] = {"fold E10.5": {"mean": {"organ:f": 1.0}}}
    assert _gate_reasons(h34, committed=c)[0] == []


def test_g1_dry_runs_on_committed_results(h34):
    """Rule 9 J7, in the tests step, which every later step waits for: the script's own G1 arithmetic on committed
    results other than D_f's. Row H27's committed per-seed parts of A_f, B_f and C_f, through G1's means, give the
    means row H29b committed for those three controls, exactly (54 of 54 in the eighteen entries G1 reads); through
    G1's paired statistics they give the paired differences row H27 committed (216 of 216: mean, lo, hi, ahead). A
    copy with one seed's scores moved by one part in 1e12 fails both. Arithmetic on committed files only: nothing is
    built or scored."""
    committed = h34.load_h29b_reading()
    assert committed is not None and h34.committed_ok(committed)  # the file of the sha256 row H34's eighth cell names
    rows = h34.load_h27_rows()
    assert rows is not None and all([r["seed"] for r in rows[k]] == list(range(10)) for k in h34.TASK_KEYS)
    h27_reading = h34.h29b.load_h27_reading()
    letter = h34.h27.LETTER

    def m_off(parts):
        return [(a, w, u) for a in h34.G1_ARMS for w, _ in h34.G1_WORLDS for u in h34.G1_UNITS
                if h34.means(parts, a)[w][u] != committed["tables"][w][u]["mean"][a]]

    def p_off(parts):
        off = []
        for x, y in h34.h26.PAIRS:
            got = h34.paired(parts, f"{x}:f", parts, f"{y}:f")
            for w, _ in h34.G1_WORLDS:
                for u in h34.G1_UNITS:
                    want = h27_reading["tables"][w][u]["diff"]["f"][f"{letter[x]}-{letter[y]}"]
                    off += [(x, y, w, u, k) for k in ("mean", "lo", "hi", "ahead") if got[w][u][k] != want[k]]
        return off
    assert {f"{x}:f" for pair in h34.h26.PAIRS for x in pair} == set(h34.G1_ARMS)
    assert len(h34.G1_ARMS) * len(h34.G1_WORLDS) * len(h34.G1_UNITS) == 54
    assert len(h34.h26.PAIRS) * len(h34.G1_WORLDS) * len(h34.G1_UNITS) * 4 == 216
    assert m_off(rows) == [] and p_off(rows) == []
    # it can fail: one seed's scores a part in 1e12 off
    shifted = json.loads(json.dumps(rows))
    for q in h34.TRUTHS:
        for a in h34.G1_ARMS:
            shifted[KEY1][4]["scores"][q][a] *= 1 + 1e-12
    assert json.loads(json.dumps(rows)) == rows and m_off(shifted) and p_off(shifted)
    assert {u for _, _, u in m_off(shifted)} <= {KEY1, "fold E10.5"}  # the other transition's entries stand
    # and the paired statistics see what most means do not: two seeds exchanged on one transition
    swapped = json.loads(json.dumps(rows))
    swapped[KEY1][2]["scores"], swapped[KEY1][7]["scores"] = swapped[KEY1][7]["scores"], swapped[KEY1][2]["scores"]
    assert p_off(swapped) and len(m_off(swapped)) < len(m_off(shifted))


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
    assert h34.load_h27_rows(tmp_path / "none") is None  # row H27's parts that are not at hand give nothing either


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
    _, found = _gate_reasons(h34)
    txt = h34.gate_text({**_gate(h34), "found": found})
    for need in ("first, counting nothing", "no counted arm and no placebo is built or scored in it",
                 "by exact equality, with no tolerance", f"sha256 {h34.H29B_READING_SHA256}",
                 f"(R) found: {KEY1}: equal; {KEY2}: equal.", f"(D) found: {KEY1}: equal; {KEY2}: equal.",
                 "(M) found: 18 means compared, 0 unequal", "(P) found: 270 paired statistics compared, 0 unequal.",
                 "G2 found: 20 contexts, 0 failed; the mapping's lineages: CM, Endothelial."):
        assert need in txt, need
    assert txt.rstrip().endswith(h34.GATE_LINES["PASS"])
    failed = h34.gate_text({**_gate(h34, "FAIL"), "why": ["G1 (M): a", "G2: b"]}).splitlines()
    assert failed[-4:] == [h34.GATE_LINES["FAIL"], "WHY: G1 (M): a; G2: b",
                           "Try H34: VOID: v3: the gate fails (G1 (M): a; G2: b); nothing is built.", NOTHING_BUILT]
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


def _write_parts(root, parts):
    for k, ps in parts.items():
        for i, p in enumerate(ps):
            d = root / re.sub(r"[^A-Za-z0-9.]+", "_", k) / str(i)
            d.mkdir(parents=True)
            (d / "h34.json").write_text(json.dumps(p, default=float))


def test_gate_run_and_read_end_to_end_on_stand_in_stages(h34, fake_io, monkeypatch, tmp_path, capsys):
    """The three steps through the command line on stand-in stages, two seeds in place of ten (the scores come from a
    stand-in scorer; the builds, the factors, the records and every check are the script's own)."""
    seeds = (0, 1)
    monkeypatch.setattr(h34, "SEEDS", seeds)
    monkeypatch.setattr(h34, "score_arms", _stub_scores(h34, GOOD))
    monkeypatch.setattr(h34, "git_state", lambda: {"commit": COMMIT, "tree_clean": True, "utc": "x"})
    # row H27's parts and row H29b's reading as these stand-in stages make them: the stand-in control scores,
    # hand-made scores of the three controls, the loader's records
    h27 = tmp_path / "h27"
    for k, rows in _h27_rows(h34, seeds).items():
        (h27 / h34.h26.safe(k)).mkdir(parents=True)
        (h27 / h34.h26.safe(k) / "h27.json").write_text(json.dumps({"axis": h34.h27.AXIS, "task": k, "rows": rows}))
    committed = _committed(h34, seeds)
    for key in h34.TASK_KEYS:
        tr = h34.load_transition(key, tmp_path)
        d = h34.build_control(h34.context(tr, 0), tr.target)[2]
        committed["pool_info"][key], committed["organ_records"][key] = tr.facts["pred_info"], tr.facts["organ_record"]
        committed["ratios"][key], committed["d_targets"][key] = tr.facts["ratios"], tr.facts["d_targets"]
        committed["pools_short"][key] = {"organ": {"pool_cells": d["pool_cells"], "drawn": [d["drawn"]],
                                                   "short_in_seeds": len(seeds) * d["short"],
                                                   "repeated_rows_at_least": d["repeated_rows_at_least"]}}
    reading = tmp_path / "h29b_reading.json"

    def commit_reading():
        """Write the stand-in reading and name it as the file row H29b's row records."""
        reading.write_text(json.dumps(committed, default=float))
        monkeypatch.setattr(h34, "H29B_READING_SHA256", hashlib.sha256(reading.read_bytes()).hexdigest())
    commit_reading()
    out = tmp_path / "results"
    gate_dir, parts_dir = out / "gate", out / "parts"
    common = ["--commit", COMMIT, "--data", str(tmp_path)]

    def gate(folder, h29b=reading, parts=h27):
        with pytest.raises(SystemExit) as e:
            h34.main(["gate", "--out", str(folder), "--h29b", str(h29b), "--h27", str(parts), *common])
        return e.value.code, capsys.readouterr().out.splitlines()
    run = ["run", "--task", KEY1, "--seeds", "0-1", "--out", str(parts_dir / "e95" / "0-1"), "--gate", str(gate_dir),
           *common]
    # the run refuses before any gate record exists
    with pytest.raises(SystemExit, match="has not said PASS"):
        h34.main(run)
    # the gate
    code, printed = gate(gate_dir)
    assert code == 0 and printed[-1] == "GATE PASS" and len(printed) == 5
    assert all(re.fullmatch(r"T2h/[^ ]+ seed [01] \([0-9.]+ s\) done", ln) for ln in printed[:-1])
    g = json.loads((gate_dir / "gate.json").read_text())
    assert g["verdict"] == "PASS" and g["why"] == [] and g["commit"] == COMMIT and g["axis"] == h34.GATE_AXIS
    assert g["found"] == {"R": {KEY1: [], KEY2: []}, "D": {KEY1: True, KEY2: True},
                          "M": {"means": 18, "unequal": [], "largest_abs_diff": 0.0},
                          "P": {"numbers": 270, "unequal": []},
                          "G2": {"contexts": 4, "failed": 0, "lineages": ["CM", "Endothelial"]}}
    assert g["h29b_reading_sha256"] == h34.H29B_READING_SHA256 and g["provenance"]["commit"] == COMMIT
    assert set(g["rows"]) == set(h34.TASK_KEYS) and [r["seed"] for r in g["rows"][KEY1]] == [0, 1]
    assert set(g["rows"][KEY1][0]["scores"]["own"]) == {"floor", "control"}  # only the control is built and scored
    assert g["rows"][KEY2][1]["factors"] == FACTORS
    assert h34.GATE_LINES["PASS"] in (gate_dir / "gate.txt").read_text() and h34.gate_passed(gate_dir, COMMIT)
    with pytest.raises(SystemExit, match="already holds a gate record"):
        h34.main(["gate", "--out", str(gate_dir), "--h29b", str(reading), "--h27", str(h27), *common])
    # the run: two parts, one per transition; the four arms, the control once more, the floor's scores kept
    h34.main(run)
    h34.main(["run", "--task", KEY2, "--seeds", "0-1", "--out", str(parts_dir / "e875" / "0-1"), "--gate",
              str(gate_dir), *common])
    printed = capsys.readouterr().out.splitlines()
    assert [bool(re.fullmatch(r"T2h/[^ ]+ seed [01] \([0-9.]+ s\) done", ln)) for ln in printed] == \
        [True, True, False] * 2
    assert printed[2] == f"wrote {parts_dir / 'e95' / '0-1' / 'h34.json'}"
    with pytest.raises(SystemExit, match="exists"):
        h34.main(run)
    parts, unreadable = h34.load_parts(parts_dir)
    assert unreadable == [] and set(parts) == set(h34.TASK_KEYS) and all(len(ps) == 1 for ps in parts.values())
    assert set(parts[KEY1][0]["rows"][0]["scores"]["own"]) == set(h34.ALL)
    # the read
    read = ["read", "--runs", str(parts_dir), "--gate", str(gate_dir), "--out", str(out), "--commit", COMMIT]
    h34.main(read)
    res = json.loads((out / "reading.json").read_text())
    assert res["void"] is False and res["why"] == [] and res["control_scores"] == "the gate's, the ones G1 checked"
    # two seeds cannot reach the bar's eight: clause (b) holds here whatever the scores, and the read says so
    assert res["outcome"] == "STOPPED" and res["clauses"] == ["b"]
    assert capsys.readouterr().out == "Try H34: STOPPED (b): OTS15 is not built.\n" + NOTHING_BUILT + "\n"
    assert (out / "reading.txt").read_text().startswith("Try H34: STOPPED (b): OTS15 is not built.\n")
    assert res["trend_record"][KEY1]["ref_cells"] == committed["pool_info"][KEY1]["organ"]["pool_cells"]
    assert all(v > 0 for v in res["changed_counts"][KEY2].values())

    def read_into(name, runs=parts_dir, gate_folder=gate_dir):
        h34.main(["read", "--runs", str(runs), "--gate", str(gate_folder), "--out", str(tmp_path / name),
                  "--commit", COMMIT])
        return capsys.readouterr().out.splitlines()
    # v5: a run-part control score that is not the gate's; v4: a file named as a part that cannot be read
    moved = json.loads(json.dumps(parts))
    moved[KEY1][0]["rows"][1]["scores"]["radial_0"]["control"] += 1e-9
    _write_parts(tmp_path / "moved", moved)
    assert read_into("out_moved", runs=tmp_path / "moved")[0] == (
        f"Try H34: VOID: v5: {KEY1} seed 1 radial_0: the run step's control score is not the gate's; nothing is "
        "built.")
    _write_parts(tmp_path / "torn", parts)
    (tmp_path / "torn" / "x").mkdir()
    (tmp_path / "torn" / "x" / "h34.json").write_text("{not json")
    (tmp_path / "torn" / "y").mkdir()
    (tmp_path / "torn" / "y" / "h34.json").write_text(json.dumps({"axis": "h29b-heart-organ", "task": KEY1}))
    assert read_into("out_torn", runs=tmp_path / "torn")[0] == (
        f"Try H34: VOID: v4: a part cannot be read: {tmp_path / 'torn' / 'x' / 'h34.json'}; v4: a part cannot be "
        f"read: {tmp_path / 'torn' / 'y' / 'h34.json'}; nothing is built.")
    # a committed mean a part in 1e12 off: the gate fails, exits 3 and states the VOID outcome itself
    committed["tables"]["CLEAN"]["fold E10.5"]["mean"]["organ:f"] *= 1 + 1e-12
    commit_reading()
    code, printed = gate(tmp_path / "gate2")
    assert code == 3 and printed[-3] == "GATE FAIL" and printed[-1] == NOTHING_BUILT
    assert printed[-2].startswith("Try H34: VOID: v3: the gate fails (G1 (M): 1 of 18 means are not row H29b's "
                                  "committed D_f means (first: CLEAN fold E10.5, ")
    assert printed[-2].endswith(")); nothing is built.")
    g2 = json.loads((tmp_path / "gate2" / "gate.json").read_text())
    assert g2["verdict"] == "FAIL" and len(g2["why"]) == 1 and not h34.gate_passed(tmp_path / "gate2", COMMIT)
    with pytest.raises(SystemExit, match="has not said PASS"):
        h34.main([*run[:6], str(parts_dir / "e95" / "again"), "--gate", str(tmp_path / "gate2"), *common])
    assert read_into("out_gate2", gate_folder=tmp_path / "gate2")[0] == (
        "Try H34: VOID: v3: the gate said 'FAIL', not PASS; nothing is built.")
    # a reading that is not the file row H29b's row records, or a missing one; row H27's parts not at hand
    other = reading.with_name("other.json")
    other.write_text(reading.read_text() + "\n")
    for i, path in enumerate((other, tmp_path / "none.json")):
        code, printed = gate(tmp_path / f"gate{3 + i}", h29b=path)
        assert code == 3 and json.loads((tmp_path / f"gate{3 + i}" / "gate.json").read_text())["why"] == [
            "G1: row H29b's reading.json is not at hand as the file its row records (missing, another sha256, or "
            "void)"]
    committed["tables"]["CLEAN"]["fold E10.5"]["mean"]["organ:f"] /= 1 + 1e-12
    committed["tables"]["own"][KEY1]["stated"]["organ:f - full_iso_half:f"]["se"] *= 1 + 1e-12
    commit_reading()
    code, printed = gate(tmp_path / "gate5")
    assert code == 3 and printed[-2] == (
        "Try H34: VOID: v3: the gate fails (G1 (P): 1 of 270 paired statistics are not row H29b's committed (first: "
        f"own {KEY1} organ:f - full_iso_half:f se)); nothing is built.")
    code, printed = gate(tmp_path / "gate6", parts=tmp_path / "no_h27")
    assert code == 3 and "G1 (P): row H27's committed parts are not at hand" in printed[-2]
    # tracked files or HEAD changing while the gate runs: no gate record at all
    monkeypatch.setattr(h34, "git_state", lambda: {"commit": COMMIT, "tree_clean": False, "utc": "x"})
    with pytest.raises(SystemExit, match="changed while the gate ran; no gate record is written"):
        h34.main(["gate", "--out", str(tmp_path / "gate7"), "--h29b", str(reading), "--h27", str(h27), *common])
    assert not (tmp_path / "gate7").exists()


def test_the_read_step_writes_one_reading_and_never_overwrites_it(h34, monkeypatch, tmp_path, capsys):
    monkeypatch.setattr(h34, "refuse_state", lambda commit, what="run": None)
    runs = tmp_path / "results" / "parts"
    _write_parts(runs, _parts(h34, score=GOOD))
    (runs / "other.json").write_text("{}")
    gate_dir = tmp_path / "results" / "gate"
    gate_dir.mkdir()
    (gate_dir / "gate.json").write_text(json.dumps(_gate(h34)))
    out = tmp_path / "results"
    args = ["read", "--runs", str(runs), "--gate", str(gate_dir), "--out", str(out), "--commit", COMMIT]
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
    for name in ("load_transition", "load_parts", "load_h29b_reading", "load_h27_rows"):
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


def _asks(tmp_path, commit=COMMIT):
    """The launcher's question before it marks a step as started: the step's refusals and nothing else."""
    return {"gate": ["refusals", "--step", "gate", "--commit", commit, "--out", str(tmp_path / "go")],
            "run": ["refusals", "--step", "run", "--commit", commit, "--gate", str(tmp_path / "gate")],
            "read": ["refusals", "--step", "read", "--commit", commit, "--out", str(tmp_path / "r")]}


def test_every_step_refuses_without_the_row_before_anything_loads(h34, monkeypatch, tmp_path, capsys):
    _nothing_loads(h34, monkeypatch)
    monkeypatch.setattr(h34, "git_state", lambda: {"commit": COMMIT, "tree_clean": True, "utc": "x"})
    calls, asks = _calls(h34, tmp_path), _asks(tmp_path)
    # a ledger without the row: other rows, and rows that only mention it, are not row H34
    monkeypatch.setattr(h34, "TRIES", _ledger(tmp_path, ["H28 (one row)", "H29b (the form; see H34)", "H3", "H340",
                                                         "R32 (the captain's order; ruling, no run; H34)"]))
    assert h34.tries_rows() == {"H28", "H29b", "H3", "H340", "R32"}
    for step in calls:
        for args in (calls[step], asks[step]):
            with pytest.raises(SystemExit, match="H34 is refused: it has no row in experiments/harness/TRIES.md") as e:
                h34.main(args)
            assert f"a row comes before its {step}" in str(e.value)
    assert not any((tmp_path / d).exists() for d in ("go", "o", "r"))
    # with the row, the same calls pass every refusal and reach their first load; the question alone loads nothing,
    # prints nothing and writes nothing
    monkeypatch.setattr(h34, "TRIES", _ledger(tmp_path, ["H28 (one row)", "H34 (one row, one script; eight tries)"]))
    assert "H34" in h34.tries_rows()
    for step in calls:
        assert h34.main(asks[step]) is None
        with pytest.raises(AssertionError, match="something loaded"):
            h34.main(calls[step])
    assert capsys.readouterr().out == "" and not any((tmp_path / d).exists() for d in ("go", "o", "r"))


def test_the_other_refusals_of_every_step(h34, monkeypatch, tmp_path):
    _nothing_loads(h34, monkeypatch)
    monkeypatch.setattr(h34, "tries_rows", lambda path=None: {"H34"})
    state = {"commit": COMMIT, "tree_clean": True, "utc": "x"}
    monkeypatch.setattr(h34, "git_state", lambda: dict(state))
    calls, asks = _calls(h34, tmp_path), _asks(tmp_path)
    both = [*calls.values(), *asks.values()]
    # tracked files that differ from the commit
    for dirty in (False, None):
        state["tree_clean"] = dirty
        for args in both:
            with pytest.raises(SystemExit, match="H34 is refused: tracked files differ from the commit"):
                h34.main(args)
    state["tree_clean"] = True
    # a HEAD other than the commit given
    state["commit"] = "d" * 40
    for args in both:
        with pytest.raises(SystemExit, match="H34 is refused: HEAD is not the run commit"):
            h34.main(args)
    state["commit"] = COMMIT
    for other in (COMMIT[:12], ""):
        for args in (*_calls(h34, tmp_path, other).values(), *_asks(tmp_path, other).values()):
            with pytest.raises(SystemExit, match="HEAD is not the run commit"):
                h34.main(args)
    # vec from another checkout
    monkeypatch.setattr(h34, "vec_home", lambda: tmp_path)
    for args in both:
        with pytest.raises(SystemExit, match="H34 is refused: vec is imported from .*, not from this checkout"):
            h34.main(args)
    monkeypatch.setattr(h34, "vec_home", lambda: h34.HERE)
    # the run: the gate record of this commit must say PASS
    g = tmp_path / "gate" / "gate.json"
    for rec in (_gate(h34, "FAIL"), _gate(h34, commit="d" * 40), _gate(h34, clean=False),
                {"verdict": "PASS", "commit": COMMIT}):
        g.write_text(json.dumps(rec))
        for args in (calls["run"], asks["run"]):
            with pytest.raises(SystemExit, match="H34 is refused: the gate of .* has not said PASS"):
                h34.main(args)
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
        for args in (calls["gate"], asks["gate"]):
            with pytest.raises(SystemExit, match="already holds a gate record"):
                h34.main(args)
        (tmp_path / "go" / f).unlink()
    for f in ("reading.json", "reading.txt"):
        (tmp_path / "r").mkdir(exist_ok=True)
        (tmp_path / "r" / f).write_text("x")
        for args in (calls["read"], asks["read"]):
            with pytest.raises(SystemExit, match="already holds a reading"):
                h34.main(args)
        (tmp_path / "r" / f).unlink()
    for args in (calls["gate"], calls["read"]):
        with pytest.raises(AssertionError, match="something loaded"):
            h34.main(args)
    # the question names its step's folder
    with pytest.raises(SystemExit, match="refusals --step gate needs --out"):
        h34.main(["refusals", "--step", "gate", "--commit", COMMIT])
    with pytest.raises(SystemExit, match="refusals --step read needs --out"):
        h34.main(["refusals", "--step", "read", "--commit", COMMIT, "--gate", str(tmp_path / "gate")])
    with pytest.raises(SystemExit, match="refusals --step run needs --gate"):
        h34.main(["refusals", "--step", "run", "--commit", COMMIT, "--out", str(tmp_path / "o")])


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
    assert gate[2:] == ["print(follows_line(VOID), flush=True)",
                        "print(outcome_line(VOID, gate_void(why)), flush=True)"]
    assert sorted(_prints(h34.cmd_read)) == ["print(res['follows'])", "print(res['line'])"]
    own = [fn for name, fn in inspect.getmembers(h34, inspect.isfunction) if fn.__module__ == h34.__name__]
    assert len(own) > 30
    for fn in own:
        if fn not in (h34.run_transition, h34.run_gate, h34.cmd_read):
            assert _prints(fn) == [], fn.__name__
    # the fixed words hold names (rows, steps, the file, the two strengths) and no figure
    names = r"H29b|H34|R32|G1|G2|OTS15|step \(3\)|step 4|step 5|strength-1\.5|1\.5 times|strength-1 "
    for words in (*h34.FOLLOWS.values(), h34.SCREEN, h34.NO_PRIOR, *h34.GATE_LINES.values()):
        assert not re.search(r"\d", re.sub(names, "", words)), words


# ================================================================================================ the launcher
# A stand-in for uv: it records each call with the environment the launcher gave it and, at that moment, the call's
# exit file for the step, the call's done.exit and the step's once-only exit file; it imitates each kind of call the
# launcher makes.
STUB_UV = r"""#!/usr/bin/env bash
BASE="$HOME/v/scratch/h34"
S="$BASE/$STUB_RUN"
threads="$OMP_NUM_THREADS|$OPENBLAS_NUM_THREADS|$MKL_NUM_THREADS"
echo "$UV_PROJECT_ENVIRONMENT|$threads|$NUMPY_MADVISE_HUGEPAGE|$TMPDIR|$*" >> "$HOME/uvcalls"
args=("$@")
[ "${args[0]}" = run ] || exit 90
args=("${args[@]:1}")
while [ "${args[0]}" = --frozen ] || [ "${args[0]}" = --no-sync ]; do args=("${args[@]:1}"); done
[ "${args[0]}" = python ] || exit 91
args=("${args[@]:1}")
while [ "${args[0]}" = -I ] || [ "${args[0]}" = -B ]; do args=("${args[@]:1}"); done
out=""; seeds=""; step=""; prev=""
for a in "${args[@]}"; do
  [ "$prev" = "--out" ] && out="$a"
  [ "$prev" = "--seeds" ] && seeds="$a"
  [ "$prev" = "--step" ] && step="$a"
  prev="$a"
done
state() {
  once=$(cat "$BASE/$2.exit" 2>/dev/null || echo none)
  echo "$1 $2 $(cat "$S/$2.exit" 2>/dev/null) $(cat "$S/done.exit" 2>/dev/null) $once" >> "$HOME/states"
}
case "${args[0]}" in
  -m) state call tests; echo "some test output"; exit "${STUB_TESTS_EXIT:-0}" ;;
  experiments/harness/h34_heart_organ_trend.py)
    case "${args[1]}" in
      refusals)
        state ask "$step"
        if [ "${STUB_REFUSE:-}" = "$step" ]; then echo "H34 is refused: a reason"; exit 1; fi
        if [ "${STUB_RACE:-}" = "$step" ]; then echo running > "$BASE/$step.exit"; fi
        exit 0 ;;
      gate)
        state call gate
        if [ "${STUB_GATE_EXIT:-0}" = kill ]; then kill -9 "$PPID"; exit 137; fi
        if [ "${STUB_GATE_EXIT:-0}" = 1 ]; then echo "Traceback"; exit 1; fi
        mkdir -p "$out"; echo '{"stub": "gate"}' > "$out/gate.json"
        echo "T2h/<=E9.5->E10.5 seed 0 (1.0 s) done"; echo "a library warning"
        if [ "${STUB_GATE_EXIT:-0}" = 3 ]; then
          echo "GATE FAIL"; echo "Try H34: VOID: v3: the gate fails (G1 (M): x); nothing is built."
          echo "A screen, not evidence of development. What follows: nothing is built."
        else
          echo "GATE PASS"
        fi
        exit "${STUB_GATE_EXIT:-0}" ;;
      run)
        state call run
        if [ "${STUB_RUN_BAD:-}" = "$seeds" ]; then echo "boom"; exit 7; fi
        mkdir -p "$out"
        [ "${STUB_RUN_NOFILE:-}" = "$seeds" ] || echo '{"stub": "run"}' > "$out/h34.json"
        echo "T2h/x seed 0 (1.0 s) done"; echo "a library warning"; echo "wrote $out/h34.json"; exit 0 ;;
      read)
        state call read
        if [ "${STUB_READ_EXIT:-0}" = kill ]; then kill -9 "$PPID"; exit 137; fi
        [ "${STUB_READ_EXIT:-0}" = 0 ] || exit "$STUB_READ_EXIT"
        mkdir -p "$out"; echo "the reading" > "$out/reading.txt"; echo "{}" > "$out/reading.json"
        echo "Try H34: NOT STOPPED."; echo "A screen, not evidence of development. What follows: x"; exit 0 ;;
    esac ;;
esac
exit 92
"""
OUT_DIR = "experiments/harness/results/h34_heart_organ_trend"
STEPS = ("check", "tests", "gate", "run", "read")
ONE_CALL = ["tests", "ask gate", "gate", "ask run", *["run"] * 10, "ask read", "read"]  # the stand-in's calls


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


def _kinds(repo):
    """The stand-in's calls in order, by kind: tests, ask <step> (the script's refusals for a step), gate, run, read."""
    out = []
    for c in _uv_calls(repo):
        words = c.split("python -I -B ", 1)[1].split()
        out.append("tests" if words[0] == "-m" else f"ask {words[3]}" if words[1] == "refusals" else words[1])
    return out


def _exits(repo, run):
    """One call's exit files, in its own scratch folder."""
    return {p.name.removesuffix(".exit"): p.read_text().strip() for p in sorted((repo.scratch / run).glob("*.exit"))}


def _once(repo):
    """The once-only exit files: one for each gate, run or read ever started, whatever the call and the commit."""
    return {p.name.removesuffix(".exit"): p.read_text().strip() for p in sorted(repo.scratch.glob("*.exit"))}


def _next_commit(repo):
    (repo.root / "note.txt").write_text("x\n")
    repo.git("add", "note.txt")
    repo.git("commit", "-q", "-m", "two")
    return repo.git("rev-parse", "HEAD")


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
    at = {s: code.index(f"STEP={s}\n") for s in STEPS}
    assert list(at.values()) == sorted(at.values()) and code.count("STEP=") == 5
    assert code.count("--untracked-files=no") == 1 and "git status --porcelain --untracked-files=no" in code
    assert f"OUT={OUT_DIR}\n" in code  # the gate's record, the parts and the reading go to the results folder
    # every output of a step goes to a file; nothing is printed but the launcher's own refusals
    assert code.count('> "$S/tests.log" 2>&1') == 1 and code.count('> "$S/gate.out" 2>&1') == 1
    assert code.count('> "$S/logs/$short-$seeds.log" 2>&1') == 1 and code.count('> "$S/read.out" 2>&1') == 1
    assert code.count('> "$S/$STEP.refusals" 2>&1') == 1
    assert all(">&2" in ln for ln in code.splitlines() if re.search(r"\becho \"(refused|usage)", ln))
    for test in ("tests/test_h34_heart_organ_trend.py", "tests/test_h29b_heart_organ.py",
                 "tests/test_h27_heart_forecast_size.py"):
        assert test in code, test
    # the gate, the run and the read are each started through start_once, before the step's own call; the tests step,
    # which may be repeated, is not. The once-only file is made by a write that cannot replace a file
    assert re.findall(r"^start_once (--[a-z]+ \S+)$", code, flags=re.M) == [
        '--out "$OUT/gate"', '--gate "$OUT/gate"', '--out "$OUT"']
    assert code.count("start_once") == 4 and code.count("set -o noclobber; echo running > \"$BASE/$STEP.exit\"") == 1
    for step in ("gate", "run", "read"):
        assert at[step] < code.index("\nstart_once --", at[step]) < code.index(f'$PY "$SCRIPT" {step} ', at[step])
    assert code.index("STEP=tests\n") < code.index("-m pytest ") < code.index("\nstart_once --")
    # within start_once: the step's file (any content, any commit: nothing of it is read), this launcher's refusal,
    # the script's refusals, and only then the start
    body = code[code.index("start_once() {"):code.index("STEP=tests\n")]
    order = [body.index(s) for s in ('if [ -e "$BASE/$STEP.exit" ]; then', "finish 88", "finish 92",
                                     '"$SCRIPT" refusals --step ', "finish 89", "set -o noclobber", "OWNED=$STEP",
                                     '> "$BASE/$STEP.head"')]
    assert order == sorted(order) and "cat " not in body and "return" not in body


def test_the_launcher_runs_the_four_steps_in_order(repo):
    code, err = _launch(repo, repo.head, "r1")
    assert (code, err) == (0, "")
    S = repo.scratch / "r1"
    assert _exits(repo, "r1") == {s: "0" for s in (*STEPS, "done")}
    assert (S / "head").read_text().strip() == repo.head
    # the once-only files: one for each of the gate, the run and the read, none for the tests
    assert _once(repo) == {"gate": "0", "run": "0", "read": "0"}
    assert all((repo.scratch / f"{s}.head").read_text() == f"{repo.head} r1\n" for s in ("gate", "run", "read"))
    # while a step runs, its exit file and done.exit hold the word running; the script's refusals are asked before
    # the step's once-only file exists, and the step is called only once that file is there
    states = (repo.tmp / "home" / "states").read_text().splitlines()
    assert states == ["call tests running running none", "ask gate running running none",
                      "call gate running running running", "ask run running running none",
                      *["call run running running running"] * 10, "ask read running running none",
                      "call read running running running"]
    calls = _uv_calls(repo)
    env = f"{repo.root}/.venv|8|8|8|0|{S}/tmp|run --frozen --no-sync python -I -B "
    assert len(calls) == 16 and all(c.startswith(env) for c in calls) and _kinds(repo) == ONE_CALL
    script = "experiments/harness/h34_heart_organ_trend.py"
    gate = f"{OUT_DIR}/gate"
    assert calls[0] == (env + f"-m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp={S}/tmp/pytest "
                        "tests/test_h34_heart_organ_trend.py tests/test_h29b_heart_organ.py "
                        "tests/test_h27_heart_forecast_size.py")
    assert calls[1] == env + f"{script} refusals --step gate --commit {repo.head} --out {gate}"
    assert calls[2] == env + f"{script} gate --out {gate} --commit {repo.head}"
    assert calls[3] == env + f"{script} refusals --step run --commit {repo.head} --gate {gate}"
    parts = ("0-1", "2-3", "4-5", "6-7", "8-9")
    for block, (key, short) in zip((calls[4:9], calls[9:14]), ((KEY1, "e95"), (KEY2, "e875")), strict=True):
        assert sorted(block) == sorted(env + f"{script} run --task {key} --seeds {p} --gate {gate} --commit "
                                       f"{repo.head} --out {OUT_DIR}/parts/{short}/{p}" for p in parts)
    assert calls[14] == env + f"{script} refusals --step read --commit {repo.head} --out {OUT_DIR}"
    assert calls[15] == (env + f"{script} read --runs {OUT_DIR}/parts --gate {gate} --commit {repo.head} --out "
                         f"{OUT_DIR}")
    # the gate's record, the parts and the reading are in the results folder; the logs, the exit files and the
    # counted lines (one line of each log is neither of the expected kinds) in the call's scratch folder
    assert (repo.root / gate / "gate.json").exists()
    assert len(sorted((repo.root / OUT_DIR / "parts").rglob("h34.json"))) == 10
    for short in ("e95", "e875"):
        for p in parts:
            assert (repo.root / OUT_DIR / "parts" / short / p / "h34.json").exists()
            assert (S / "exit" / f"{short}-{p}").read_text().strip() == "0"
            assert (S / "exit" / f"{short}-{p}.other_lines").read_text().strip() == "1"
            assert len((S / "logs" / f"{short}-{p}.log").read_text().splitlines()) == 3
    assert (S / "exit" / "other_lines").read_text().strip() == "10"
    assert (S / "gate.other_lines").read_text().strip() == "1"
    assert (S / "gate.out").read_text().splitlines()[-1] == "GATE PASS"
    assert (S / "read.out").read_text().splitlines() == [
        "Try H34: NOT STOPPED.", "A screen, not evidence of development. What follows: x"]
    assert (S / "tests.log").read_text() == "some test output\n"
    assert all((S / f"{s}.refusals").read_text() == "" for s in ("gate", "run", "read"))
    assert (repo.root / OUT_DIR / "reading.txt").read_text() == "the reading\n"
    assert repo.git("status", "--porcelain", "--untracked-files=no") == ""
    # the run folder is never reused, and a later call at the same commit finds the reading: nothing is touched
    before = {p: p.read_bytes() for p in repo.scratch.rglob("*") if p.is_file()}
    code, err = _launch(repo, repo.head, "r1")
    assert code == 91 and "the run folder" in err and "exists" in err
    assert before == {p: p.read_bytes() for p in repo.scratch.rglob("*") if p.is_file()} and len(_uv_calls(repo)) == 16
    code, err = _launch(repo, repo.head, "r2")
    assert code == 95 and "already holds a reading" in err
    assert _exits(repo, "r2") == {"check": "95", "done": "95"} and len(_uv_calls(repo)) == 16
    assert (repo.root / OUT_DIR / "reading.txt").read_text() == "the reading\n"
    assert _once(repo) == {"gate": "0", "run": "0", "read": "0"}
    # with the reading out of the way a later call gets as far as the gate, which was started once already: at this
    # commit and at another one it is refused, and the script is neither asked nor called
    (repo.root / OUT_DIR / "reading.txt").rename(repo.tmp / "aside.txt")
    (repo.root / OUT_DIR / "reading.json").rename(repo.tmp / "aside.json")
    for i in range(2):
        code, err = _launch(repo, _next_commit(repo) if i else repo.head, f"r{3 + i}")
        assert code == 88 and "the gate step was started before" in err and "it is not started twice" in err
        assert _exits(repo, f"r{3 + i}") == {"check": "0", "tests": "0", "gate": "88", "done": "88"}
        assert _kinds(repo) == [*ONE_CALL, *["tests"] * (1 + i)]
    assert _once(repo) == {"gate": "0", "run": "0", "read": "0"} and not (repo.root / OUT_DIR / "reading.txt").exists()
    assert all((repo.scratch / f"{s}.head").read_text() == f"{repo.head} r1\n" for s in ("gate", "run", "read"))


def test_the_launchers_environment_overrides(repo):
    venv, data = repo.tmp / "shared env", repo.tmp / "data raw"
    # a named environment or data folder that is not there: refused before any step
    got, err = _launch(repo, repo.head, "v0", H34_VENV=str(venv))
    assert got == 99 and "no virtual environment at" in err and _exits(repo, "v0") == {"check": "99", "done": "99"}
    got, err = _launch(repo, repo.head, "d0", H34_DATA=str(data))
    assert got == 97 and "no T2 stages under" in err and _exits(repo, "d0") == {"check": "97", "done": "97"}
    assert _uv_calls(repo) == [] and _once(repo) == {} and not (repo.root / OUT_DIR).exists()
    _stand_in_env(venv)
    (data / "T2").mkdir(parents=True)
    code, err = _launch(repo, repo.head, "e1", H34_VENV=str(venv), H34_DATA=str(data))
    assert (code, err) == (0, "")
    calls = _uv_calls(repo)
    assert len(calls) == 16 and all(c.startswith(f"{venv}|8|8|8|0|") for c in calls)
    told = [k for k, c in zip(_kinds(repo), calls, strict=True) if c.endswith(f" --data {data}")]
    assert told == ["gate", *["run"] * 10]  # the gate and the parts load stages; nothing else is told the folder
    assert sum("--data" in c for c in calls) == 11


@pytest.mark.parametrize("extra, step, code, n_calls, once", [
    ({"STUB_TESTS_EXIT": "1"}, "tests", 1, 1, {}),
    ({"STUB_TESTS_EXIT": "5"}, "tests", 5, 1, {}),
    ({"STUB_REFUSE": "gate"}, "gate", 89, 2, {}),
    ({"STUB_GATE_EXIT": "3"}, "gate", 3, 3, {"gate": "3"}),
    ({"STUB_GATE_EXIT": "1"}, "gate", 1, 3, {"gate": "1"}),
    ({"STUB_REFUSE": "run"}, "run", 89, 4, {"gate": "0"}),
    ({"STUB_RUN_BAD": "4-5"}, "run", 1, 9, {"gate": "0", "run": "1"}),
    ({"STUB_RUN_NOFILE": "8-9"}, "run", 1, 9, {"gate": "0", "run": "1"}),
    ({"STUB_REFUSE": "read"}, "read", 89, 15, {"gate": "0", "run": "0"}),
    ({"STUB_READ_EXIT": "1"}, "read", 1, 16, {"gate": "0", "run": "0", "read": "1"}),
])
def test_the_launcher_stops_at_the_first_step_that_does_not_exit_0(repo, extra, step, code, n_calls, once):
    got, err = _launch(repo, repo.head, "s1", **extra)
    assert got == code
    want = {s: "0" for s in STEPS[:STEPS.index(step)]} | {step: str(code), "done": str(code)}
    assert _exits(repo, "s1") == want  # no later step wrote an exit file: none ran
    assert _kinds(repo) == ONE_CALL[:n_calls]
    # a step the script refuses is not started: it leaves no once-only file; a started step leaves its exit code
    assert _once(repo) == once
    S = repo.scratch / "s1"
    if code == 89:
        assert f"the script refuses the {step} step" in err and step not in _kinds(repo)
        assert (S / f"{step}.refusals").read_text() == "H34 is refused: a reason\n"
        assert not (repo.scratch / f"{step}.head").exists()
    if step == "gate" and code == 3:  # a failed gate: its fixed lines are kept, and counted as expected ones
        assert (S / "gate.out").read_text().splitlines()[-3:] == [
            "GATE FAIL", "Try H34: VOID: v3: the gate fails (G1 (M): x); nothing is built.",
            "A screen, not evidence of development. What follows: nothing is built."]
        assert (S / "gate.other_lines").read_text().strip() == "1" and not (S / "logs").exists()
    if step == "run" and code == 1:  # the failed transition's five parts ran together; the second never started
        assert not (repo.root / OUT_DIR / "parts" / "e875").exists() and (S / "exit" / "other_lines").exists()
    if step != "read" or code != 0:
        assert not (repo.root / OUT_DIR / "reading.txt").exists()


@pytest.mark.parametrize("extra, first, once", [
    ({"STUB_GATE_EXIT": "3"}, 3, {"gate": "3"}),  # the gate said FAIL
    ({"STUB_GATE_EXIT": "1"}, 1, {"gate": "1"}),  # the gate crashed: no verdict
    ({"STUB_GATE_EXIT": "kill"}, -9, {"gate": "running"}),  # the call was killed in the gate: no verdict
    ({"STUB_RUN_BAD": "4-5"}, 1, {"gate": "0", "run": "1"}),  # a part failed
    ({"STUB_RUN_NOFILE": "8-9"}, 1, {"gate": "0", "run": "1"}),  # a part wrote no file
    ({"STUB_READ_EXIT": "1"}, 1, {"gate": "0", "run": "0", "read": "1"}),  # the read crashed
    ({"STUB_READ_EXIT": "kill"}, -9, {"gate": "0", "run": "0", "read": "running"}),  # the call was killed in the read
])
def test_no_gate_run_or_read_is_started_a_second_time(repo, extra, first, once):
    """Row H34, no step twice: once a gate, run or read was started, however it ended (a FAIL, a crash, a killed
    call that left no verdict), no later call starts it again, at this commit or another, under any run name; the
    tests step, which may be repeated, runs again in every call."""
    assert _launch(repo, repo.head, "a1", **extra)[0] == first
    assert _once(repo) == once
    if first == -9:  # a killed call: nothing was written after the step's start
        assert _exits(repo, "a1")["done"] == "running" and list(once.values())[-1] == "running"
    heads = {s: (repo.scratch / f"{s}.head").read_text() for s in once}
    assert all(h == f"{repo.head} a1\n" for h in heads.values())
    n = len(_uv_calls(repo))
    for i in range(2):  # a later call at the same commit, then at another
        got, err = _launch(repo, _next_commit(repo) if i else repo.head, f"a{2 + i}")
        assert got == 88 and "the gate step was started before" in err and "it is not started twice" in err
        assert _exits(repo, f"a{2 + i}") == {"check": "0", "tests": "0", "gate": "88", "done": "88"}
        assert _kinds(repo)[n:] == ["tests"] * (1 + i)  # the script is neither asked about the step nor called
    # the first call's record stands
    assert _once(repo) == once and {s: (repo.scratch / f"{s}.head").read_text() for s in once} == heads


@pytest.mark.parametrize("step", ["gate", "run", "read"])
@pytest.mark.parametrize("holds", ["running", "0", "3", ""])
def test_each_of_the_three_steps_has_its_own_once_only_file(repo, step, holds):
    """A step whose exit file exists is refused whatever the file holds, exit 0 at this very commit included; the
    steps before it are started as usual, the ones after it never."""
    repo.scratch.mkdir(parents=True)
    (repo.scratch / f"{step}.exit").write_text(holds + "\n")
    (repo.scratch / f"{step}.head").write_text(f"{repo.head} elsewhere\n")
    got, err = _launch(repo, repo.head, "o1")
    assert got == 88 and f"the {step} step was started before" in err
    assert _exits(repo, "o1") == {s: "0" for s in STEPS[:STEPS.index(step)]} | {step: "88", "done": "88"}
    assert _kinds(repo) == ONE_CALL[:ONE_CALL.index(f"ask {step}")]  # the step is neither asked about nor called
    assert (repo.scratch / f"{step}.exit").read_text() == holds + "\n"
    assert (repo.scratch / f"{step}.head").read_text() == f"{repo.head} elsewhere\n"
    assert _once(repo) == {s: "0" for s in STEPS[2:STEPS.index(step)]} | {step: holds}


def test_a_refused_step_is_not_started_and_counts_nothing(repo):
    """A step is started once a call has passed the launcher's and the script's refusals: a refused call leaves no
    once-only file, so the step can still be started, once; the tests step may be repeated."""
    for run in ("t1", "t2"):  # the tests step may be repeated; a commit where it fails starts nothing
        assert _launch(repo, repo.head, run, STUB_TESTS_EXIT="1")[0] == 1
    assert _once(repo) == {} and _kinds(repo) == ["tests", "tests"]
    # the script refuses the gate: no gate was started
    got, err = _launch(repo, repo.head, "g1", STUB_REFUSE="gate")
    assert got == 89 and "the script refuses the gate step" in err
    assert _once(repo) == {} and _kinds(repo)[2:] == ["tests", "ask gate"]
    # so a later call starts it, once, and every step after it
    assert _launch(repo, repo.head, "g2") == (0, "")
    assert _kinds(repo)[4:] == ONE_CALL and _once(repo) == {"gate": "0", "run": "0", "read": "0"}
    assert all((repo.scratch / f"{s}.head").read_text() == f"{repo.head} g2\n" for s in ("gate", "run", "read"))
    kinds = _kinds(repo)
    assert (kinds.count("tests"), kinds.count("gate"), kinds.count("run"), kinds.count("read")) == (4, 1, 10, 1)


@pytest.mark.parametrize("how", ["parts", "script"])
def test_a_run_that_is_refused_after_the_gate_is_not_started(repo, how):
    """The run's refusals come before its start: parts this launcher did not make (the launcher's own refusal, before
    the script is asked) or the script's refusal leave no once-only file for the run, and no try is counted. The
    gate before it was started, so a later call stops there: stricter than the row, under which the refused run
    could still be started."""
    parts = repo.root / OUT_DIR / "parts"
    if how == "parts":
        parts.mkdir(parents=True)
        got, err = _launch(repo, repo.head, "p1")
        assert got == 92 and "made outside this launcher" in err
        assert _kinds(repo) == ["tests", "ask gate", "gate"]  # the script is not asked about the run
        parts.rmdir()
    else:
        got, err = _launch(repo, repo.head, "p1", STUB_REFUSE="run")
        assert got == 89 and "the script refuses the run step" in err
        assert _kinds(repo) == ["tests", "ask gate", "gate", "ask run"]
    assert _exits(repo, "p1") == {"check": "0", "tests": "0", "gate": "0", "run": str(got), "done": str(got)}
    assert _once(repo) == {"gate": "0"} and not (repo.scratch / "run.head").exists() and not parts.exists()
    n = len(_uv_calls(repo))
    got, err = _launch(repo, repo.head, "p2")
    assert got == 88 and "the gate step was started before" in err
    assert _once(repo) == {"gate": "0"} and _kinds(repo)[n:] == ["tests"]


def test_two_calls_cannot_both_start_a_step(repo):
    """The once-only file is made by a write that cannot replace a file: if another call makes it between this
    call's question to the script and its start, this call does not start the step and leaves that file alone."""
    got, err = _launch(repo, repo.head, "c1", STUB_RACE="gate")
    assert got == 88 and "the gate step was started before" in err
    assert _once(repo) == {"gate": "running"} and not (repo.scratch / "gate.head").exists()
    assert _exits(repo, "c1") == {"check": "0", "tests": "0", "gate": "88", "done": "88"}
    assert _kinds(repo) == ["tests", "ask gate"]  # the gate itself was never called


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
    assert _uv_calls(repo) == []  # no refusal reached a step
    for f in ("reading.json", "reading.txt"):
        (repo.root / OUT_DIR).mkdir(parents=True, exist_ok=True)
        (repo.root / OUT_DIR / f).write_text("the one reading\n")
        refused(f"k-{f[8:]}".replace(".", ""), 95, "already holds a reading", head, f"k-{f[8:]}".replace(".", ""))
        assert (repo.root / OUT_DIR / f).read_text() == "the one reading\n"
        (repo.root / OUT_DIR / f).rename(repo.tmp / f"aside-{f}")
    assert _uv_calls(repo) == [] and _once(repo) == {}  # none of these refusals started anything
    assert _launch(repo, head, "ok")[0] == 0  # and with nothing in the way it runs, the untracked file still there
    assert (repo.root / "untracked.txt").exists() and _once(repo) == {"gate": "0", "run": "0", "read": "0"}
