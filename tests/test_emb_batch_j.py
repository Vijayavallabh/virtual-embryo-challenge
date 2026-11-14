"""The embryo batch J (campaign branch campaign/emb-j): `fgw_damp_m_grp`, the movement-matched control fgw_damp_m with
the labels behind its two factors the brackets' obs["group"] in place of obs["celltype"] (so the partial coupling that
sets them is fgw_partial_grp's), and the runner experiments/harness/n50_embryo_j_arms.py. Synthetic stages only
(tests/test_emb_expr.py's); no data file, committed result or scorer is read.

Tolerances stated here: a factor against fgw_partial_grp's mean cell shift over fgw's, 1e-9 relative (the same
partial displacements; the after bracket's nearest-anchor indices are taken in the entry's frame against
fgw_partial_grp's own, the same up to the rounding of a rotation), as in tests/test_emb_batch_g.py; a cell-level ratio
against the factor it was scaled by, 1e-12 relative."""

from __future__ import annotations

import dataclasses
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pytest

from vec import harness
from vec import t2_build as tb

REPO = Path(__file__).resolve().parents[1]
HARNESS = REPO / "experiments" / "harness"
ARM = "fgw_damp_m_grp"
# sha256 (tests/test_emb_expr.digest: X float32, coords float32) of each existing arm on `_pair_stages()` with the
# default context of tests/test_emb_expr.ctx, computed before the merge on each arm's own branch: the tree of
# campaign/emb-g-fix at b33595e and the tree of campaign/emb-h at f8adf82, each extracted with `git archive` and
# imported first on sys.path with this environment (OPENBLAS/OMP/MKL threads 4). fgw and fgw_partial gave the same
# digest on both trees.
DIGESTS_BEFORE_MERGE = {
    "fgw": "4a53a67d0be4e5195b3ddf8e554b5ac2ac0b506bfb8808f747ab48cea1d8164f",  # b33595e and f8adf82
    "fgw_partial": "659fbae5513407f2f40ce6730e0418c3678304cb20eca8057c6f89da2fe6e711",  # b33595e and f8adf82
    "fgw_partial_grp": "3c6e99b0548240ff9e578ee1c6380417ebcbe0e0b01385616e1a952989afb77a",  # f8adf82
    "fgw_damp_m": "7a19fc7ec9a087d220d82e153435212ffaca2042e818e455d52061c87b8ac927",  # b33595e
    "fgw_damp": "e4c5d7ddfcfd6289d0aa5576efed3fc576de352a1e3714709e19f36cbe857be7",  # b33595e
    "fgw_part_rand": "e74d9f7ba4fca60e4d5d857d619b4ba7ff61fa773b91551412bc51cdab4322f7",  # b33595e
    "fgw_still": "9754173a91bc162a9d59e10eb87b282f75b8877731ca82c63185b139babdcac4",  # b33595e
}
COUPLING_KEYS = ("anchors_before", "anchors_after", "fgw", "moves", "rotation", "before_frame", "after_frame")


@pytest.fixture(scope="module")
def emb_expr():
    """The synthetic embryo stages and contexts of tests/test_emb_expr.py (no data file)."""
    spec = importlib.util.spec_from_file_location("test_emb_expr_helpers_j", REPO / "tests" / "test_emb_expr.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _rename_som(s: dict) -> dict:
    """The after bracket's SOM cells relabelled SOM2 (cell type only; their group stays ParaxialMesoderm)."""
    a = s["after"]
    a.obs["celltype"] = a.obs["celltype"].astype(str).replace({"SOM": "SOM2"}).to_numpy()
    return s


def _pair_stages(emb_expr) -> dict:
    """tests/test_emb_batch_h.py's pair: the cell types differ between the brackets (SOM / SOM2) and the groups agree
    in part (the after bracket's first 50 cells moved from group Endoderm to Ectoderm)."""
    s = _rename_som(emb_expr.stages())
    g = s["after"].obs["group"].astype(str).to_numpy().copy()
    g[:50] = "Ectoderm"
    s["after"].obs["group"] = g
    return s


def _ctx(emb_expr, s, seed=None):
    kw = {} if seed is None else {"seed": seed}
    return emb_expr.ctx(inputs=s, roles=("before", "after"), **kw)


def _calls(monkeypatch, coupling: bool = False) -> list:
    """Every call of interpolate_clouds as (keywords, result); the outer calls (no partial coupling) also return the
    coupling's internals (with_coupling, which changes no other output)."""
    seen = []
    real = tb.interpolate_clouds

    def spy(*a, **kw):
        if coupling and "coupling" not in kw:
            kw = {**kw, "with_coupling": True}
        res = real(*a, **kw)
        seen.append((kw, res))
        return res

    monkeypatch.setattr(tb, "interpolate_clouds", spy)
    return seen


# ------------------------------------------------------------------------------------------ 1. the existing arms
@pytest.mark.parametrize("name", sorted(DIGESTS_BEFORE_MERGE))
def test_the_seven_existing_arms_keep_their_files_from_their_own_branches(emb_expr, name):
    c = _ctx(emb_expr, _pair_stages(emb_expr))
    assert emb_expr.digest(harness.CANDIDATES[name].build(c)) == DIGESTS_BEFORE_MERGE[name]


# ------------------------------------------------------------------------------------------ 2. the factors
def test_the_twins_factors_are_fgw_partial_grps_cell_shift_over_the_entrys_and_differ_from_fgw_damp_ms(emb_expr):
    c = _ctx(emb_expr, _pair_stages(emb_expr))
    r_f, r_g, r_m, r_j = {}, {}, {}, {}
    harness._fgw(c, build_record=r_f)  # the entry, with its record
    harness._fgw_partial_grp(c, build_record=r_g)
    harness._fgw_damp_m(c, build_record=r_m)
    harness._fgw_damp_m_grp(c, build_record=r_j)
    fj, fm = r_j["anchor_move"]["factors"], r_m["anchor_move"]["factors"]
    for i in (0, 1):
        assert fj[i] == pytest.approx(r_g["shift"]["cell"][i] / r_f["shift"]["cell"][i], rel=1e-9, abs=0)
        assert abs(fj[i] - fm[i]) > 1e-6  # not fgw_damp_m's
    assert 0 < fj[0] < 1 and 0 < fj[1] < 1 and fj[0] != fj[1]
    # each bracket's cells move on average as far as under fgw_partial_grp
    np.testing.assert_allclose(r_j["shift"]["cell"], r_g["shift"]["cell"], rtol=1e-9, atol=0)
    # its share is the group overlap, fgw_partial_grp's; fgw_damp_m's is the celltype overlap
    assert r_j["anchor_move"]["s"] == r_g["partial"]["s"] != r_m["anchor_move"]["s"]
    assert r_j["anchor_move"]["anchor_moves"] == "damped_matched"


def test_the_twins_labels_are_the_groups_and_fgw_damp_m_still_gets_celltype(emb_expr, monkeypatch):
    s = _pair_stages(emb_expr)
    c = _ctx(emb_expr, s)
    seen = _calls(monkeypatch)
    harness.CANDIDATES["fgw_damp_m"].build(c)
    harness.CANDIDATES[ARM].build(c)
    outer = [kw for kw, _ in seen if "anchor_moves" in kw]
    inner = [kw for kw, _ in seen if "coupling" in kw]
    assert len(outer) == 2 and len(inner) == 2
    (km, kj), (im, ij) = outer, inner
    for side, i in (("before", 0), ("after", 1)):
        assert np.array_equal(np.asarray(km["anchor_labels"][i]), np.asarray(s[side].obs["celltype"]))
        assert np.array_equal(np.asarray(kj["anchor_labels"][i]), np.asarray(s[side].obs["group"]))
        # the partial coupling behind the factors: fgw_partial's labels for fgw_damp_m, fgw_partial_grp's for the twin
        assert np.array_equal(np.asarray(im["coupling_labels"][i]), np.asarray(s[side].obs["celltype"]))
        assert np.array_equal(np.asarray(ij["coupling_labels"][i]), np.asarray(s[side].obs["group"]))
    assert set(km) == set(kj) and "anchor_rng" not in kj and "coupling" not in kj
    for k in km:
        if k != "anchor_labels":
            assert json.dumps(km[k], default=str) == json.dumps(kj[k], default=str), k


# ------------------------------------------------------------------------------------------ 3. the displacements
def test_every_anchors_displacement_is_its_brackets_factor_times_the_entrys(emb_expr, monkeypatch):
    c = _ctx(emb_expr, _pair_stages(emb_expr))
    seen = []
    real = tb.controlled_displacements

    def spy(db, da, mode, s, rng, factors=None):
        out = real(db, da, mode, s, rng, factors=factors)
        seen.append((np.array(db), np.array(da), factors, out))
        return out

    calls = _calls(monkeypatch, coupling=True)
    monkeypatch.setattr(tb, "controlled_displacements", spy)
    rec: dict = {}
    harness._fgw_damp_m_grp(c, build_record=rec)
    (db0, da0, factors, (db, da)), = seen
    fb, fa = factors
    assert [fb, fa] == rec["anchor_move"]["factors"]
    assert np.array_equal(db, db0 * fb) and np.array_equal(da, da0 * fa)
    # db0, da0 are the entry's own displacements: its balanced barycentric partners in its own frame
    (_, res), = [(kw, r) for kw, r in calls if "anchor_moves" in kw]
    cp = res.coupling
    ib, ia, Tm, Zb, Za_b = cp["anchors_before"], cp["anchors_after"], cp["moves"], cp["before_frame"], cp["after_frame"]
    assert np.array_equal(db0, tb.barycentric(Tm, Za_b[ia]) - Zb[ib])
    assert np.array_equal(da0, tb.barycentric(Tm.T, Zb[ib]) - Za_b[ia])


def test_the_coupling_partners_and_rotation_are_the_entrys(emb_expr, monkeypatch):
    c = _ctx(emb_expr, _pair_stages(emb_expr))
    calls = _calls(monkeypatch, coupling=True)
    harness.CANDIDATES["fgw"].build(c)
    harness.CANDIDATES[ARM].build(c)
    outer = [r for kw, r in calls if "coupling" not in kw]
    assert len(outer) == 2
    entry, arm = outer
    for key in COUPLING_KEYS:
        assert np.array_equal(arm.coupling[key], entry.coupling[key]), key
    assert arm.record["rotation_deg"] == entry.record["rotation_deg"]
    fb, fa = arm.record["anchor_move"]["factors"]
    np.testing.assert_allclose(arm.coupling["shift_before"], fb * entry.coupling["shift_before"], rtol=1e-12, atol=0)
    np.testing.assert_allclose(arm.coupling["shift_after"], fa * entry.coupling["shift_after"], rtol=1e-12, atol=0)
    np.testing.assert_allclose(arm.movement["mean_cell_shift"], [fb * entry.movement["mean_cell_shift"][0],
                                                                 fa * entry.movement["mean_cell_shift"][1]],
                               rtol=1e-12, atol=0)


def test_the_drawn_cells_and_expression_rows_are_fgws(emb_expr, monkeypatch):
    c = _ctx(emb_expr, _pair_stages(emb_expr))
    fgw = harness.CANDIDATES["fgw"].build(c)
    calls = _calls(monkeypatch)
    harness.CANDIDATES["fgw"].build(c)
    arm = harness.CANDIDATES[ARM].build(c)
    entry_res, arm_res = [r for kw, r in calls if "coupling" not in kw]
    assert np.array_equal(arm_res.rows, entry_res.rows) and np.array_equal(arm_res.stage, entry_res.stage)
    assert np.array_equal(arm_res.pool_weights, entry_res.pool_weights)
    assert np.array_equal(arm.X, fgw.X) and arm.n == fgw.n
    assert not np.array_equal(arm.coords, fgw.coords)  # the idea is not idle here
    assert not np.array_equal(arm.coords, harness.CANDIDATES["fgw_damp_m"].build(c).coords)  # nor fgw_damp_m


# ------------------------------------------------------------------------------------------ 4. idle: overlap exactly 1
def test_with_groups_that_overlap_fully_the_twin_is_fgw_bit_for_bit_without_the_solver(emb_expr, monkeypatch):
    c = _ctx(emb_expr, _rename_som(emb_expr.stages()))  # cell types differ, groups in equal shares
    fgw = harness.CANDIDATES["fgw"].build(c)
    assert emb_expr.digest(harness.CANDIDATES["fgw_damp_m"].build(c)) != emb_expr.digest(fgw)  # not idle there

    def never(*a, **kw):
        raise AssertionError("the partial solver is not called at an overlap of 1")

    monkeypatch.setattr(tb, "partial_fgw_coupling", never)
    rec: dict = {}
    arm = harness._fgw_damp_m_grp(c, build_record=rec)
    assert emb_expr.digest(arm) == emb_expr.digest(fgw)
    assert np.array_equal(arm.X, fgw.X) and np.array_equal(arm.coords, fgw.coords)
    assert rec["anchor_move"]["factors"] == [1.0, 1.0] and rec["anchor_move"]["s"] == 1.0


# ------------------------------------------------------------------------------------------ 5. refusal
@pytest.mark.parametrize("side", ["before", "after"])
def test_a_bracket_without_group_labels_refuses(emb_expr, side):
    s = _pair_stages(emb_expr)
    del s[side].obs["group"]
    c = _ctx(emb_expr, s)
    with pytest.raises(ValueError, match="group"):
        harness.CANDIDATES[ARM].build(c)
    harness.CANDIDATES["fgw_damp_m"].build(c)  # fgw_damp_m reads celltype and still builds


# ------------------------------------------------------------------------------------------ 6. reproducible
def test_same_seed_same_output(emb_expr):
    c = _ctx(emb_expr, _pair_stages(emb_expr))
    a, b = harness._fgw_damp_m_grp(c), harness.CANDIDATES[ARM].build(c)
    assert emb_expr.digest(a) == emb_expr.digest(b)
    other = harness._fgw_damp_m_grp(dataclasses.replace(c, seed=c.seed + 1))
    assert emb_expr.digest(other) != emb_expr.digest(a)


# ------------------------------------------------------------------------------------------ 7. registration
def test_the_arm_is_registered_for_t2_with_fgw_damp_ms_requirements(monkeypatch):
    seen = {}
    monkeypatch.setattr(harness, "_fgw", lambda ctx, **kw: seen.setdefault(ctx, kw))
    c, m = harness.CANDIDATES[ARM], harness.CANDIDATES["fgw_damp_m"]
    assert c.name == ARM and c.tasks == m.tasks == frozenset({"T2"}) and c.requires is m.requires
    assert c.build("a") == {"anchor_moves": "damped_matched", "build_record": None, "anchor_label": "group"}
    assert m.build("b") == {"anchor_moves": "damped_matched", "build_record": None}
    assert c.params["label"] == "group" and m.params["label"] == "celltype"
    for k in ("anchor_moves", "alpha", "k", "num_iter_max", "scope"):
        assert c.params[k] == m.params[k], k


# ================================================================================ the runner
@pytest.fixture(scope="module")
def n50():
    sys.path.insert(0, str(HARNESS))
    spec = importlib.util.spec_from_file_location("n50_embryo_j_arms", HARNESS / "n50_embryo_j_arms.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_the_runner_scores_the_floor_fgw_the_arm_and_the_two_rescored_arms_on_the_nine_tasks(n50):
    assert n50.NEW == (ARM,) and n50.RESCORED == ("fgw_partial_grp", "fgw_damp_m")
    assert n50.ARMS == (harness.FLOOR, "fgw", "fgw_partial_grp", "fgw_damp_m", ARM)
    tasks = n50.nt.all_tasks()
    assert set(tasks) == {"half", "late", "new1", "new2", "e75", "e725", "x1", "x2", "x3"}
    for name in n50.ARMS[1:]:
        assert all(harness.CANDIDATES[name].applies(t) for t in tasks.values())


FAKE_SHIFT = {"fgw": {"anchor": [2.0, 4.0], "cell": [1.0, 0.0]}, "arm": {"anchor": [1.0, 1.0], "cell": [0.5, 0.0]}}
FAKE_RATIO = {"anchor": [0.5, 0.25], "cell": [0.5, None]}
SCORE = {"fgw": 51.0, "fgw_partial_grp": 52.0, "fgw_damp_m": 53.0, ARM: 54.0}
DAMPED = ("fgw_damp_m", ARM)


def _patch_builds(monkeypatch, n50, fail_on: set[str], key: str, seed_fails: int = 0):
    """Every guarded arm's build replaced (each writes a record and its shifts, the damped arms their factors): the
    arms in fail_on raise the builder's refusal from seed_fails."""
    cands = harness.CANDIDATES

    def make(name):
        def build(ctx, build_record=None):
            if name in fail_on and ctx.seed - harness.BUILD_SEED_OFFSET >= seed_fails:
                raise ValueError(f"{key}: anchor_moves='damped_matched' needs obs['group'] in both brackets")
            build_record["seen"] = 1
            if name in DAMPED:
                build_record["anchor_move"] = {"factors": [0.5, 0.25]}
            build_record["shift"] = FAKE_SHIFT["arm"]
            return "built"

        return build

    def fgw(ctx, build_record=None):
        build_record["shift"] = FAKE_SHIFT["fgw"]
        return "built"

    for x in n50.GUARDED:
        monkeypatch.setitem(cands, x, dataclasses.replace(cands[x], build=make(x)))
    monkeypatch.setitem(cands, "fgw", dataclasses.replace(cands["fgw"], build=fgw))


def _fake_runner(score: dict[str, float], ctx_of=None):
    """run_task's shape without data or scorer: each candidate is built (its build may raise) and given a constant.
    ctx_of(seed) gives the build's context (default: a context holding the seed only)."""
    calls = []
    if ctx_of is None:
        def ctx_of(seed):
            return dataclasses.make_dataclass("Ctx", ["seed"])(seed=seed)

    def runner(t, cands, seeds, data_dir):
        calls.append([c.name for c in cands])
        rows = []
        for s in seeds:
            row = {"seed": s, "scores": {"floor": 50.0}, "skills": {"floor": {}}, "metrics": {"floor": {}},
                   "n_cells": {"floor": 1}, "groups": {"floor": {}}, "reads": {}}
            for c in cands:
                if c.name == "floor":
                    continue
                c.build(ctx_of(harness.BUILD_SEED_OFFSET + s))
                row["scores"][c.name] = score[c.name]
                for k in ("skills", "metrics", "groups"):
                    row[k][c.name] = {}
                row["n_cells"][c.name], row["reads"][c.name] = 1, []
            rows.append(row)
        return {"task": {"key": t.key, "inputs": dict(t.inputs)}, "rows": rows,
                "candidates": [c.name for c in cands], "params": {c.name: {} for c in cands},
                "provenance": {"commit": "a" * 40, "tree_clean": True}}

    return runner, calls


def test_the_runners_task_file_loads_through_nine_tasks_with_the_rescored_arms_under_their_names(tmp_path,
                                                                                                 monkeypatch, n50):
    t = n50.nt.all_tasks()["e75"]
    _patch_builds(monkeypatch, n50, set(), t.key)
    runner, calls = _fake_runner(SCORE)
    monkeypatch.setattr(n50.harness, "run_task", runner)
    out = tmp_path / "runs"
    n50.main(["run", "--task", "e75", "--out", str(out), "--data", str(tmp_path / "no_data")])
    assert calls == [list(n50.ARMS)]
    (path,) = out.glob("T2e_*.json")
    data = json.loads(path.read_text())
    assert data["task_source"] == "experiments/harness/n50_embryo_j_arms.py" and data["not_buildable"] == {}
    # n46's shape: the arm's own fields, movement_ratio, and for the damped arms factors (fgw's and raw shifts not
    # written)
    want = {a: {str(s): {"seen": 1, "movement_ratio": FAKE_RATIO} for s in n50.SEEDS} for a in n50.GUARDED}
    for a in DAMPED:
        for s in n50.SEEDS:
            want[a][str(s)].update(anchor_move={"factors": [0.5, 0.25]}, factors=[0.5, 0.25])
    assert data["build_records"] == want
    scores, why, meta = n50.nt._load([out])
    assert why == [] and meta == ["runs:aaaaaaa"]
    assert set(scores[t.key]) == set(n50.ARMS) and {"fgw", "fgw_partial_grp", "fgw_damp_m"} <= set(scores[t.key])
    assert all(len(v) == 10 for v in scores[t.key].values())
    with pytest.raises(SystemExit, match="exists"):
        n50.main(["run", "--task", "e75", "--out", str(out), "--data", str(tmp_path / "no_data")])
    assert len(calls) == 1
    with pytest.raises(SystemExit, match="unknown task"):
        n50.main(["run", "--task", "nope", "--out", str(out)])


def test_an_arm_the_builder_refuses_is_recorded_and_the_others_are_still_scored(tmp_path, monkeypatch, n50):
    t = n50.nt.all_tasks()["e725"]
    _patch_builds(monkeypatch, n50, {ARM}, t.key, seed_fails=2)
    runner, calls = _fake_runner(SCORE)
    res = n50.score_task(t, tmp_path, runner=runner)
    assert calls == [list(n50.ARMS), ["floor", "fgw"], ["floor", "fgw_partial_grp"], ["floor", "fgw_damp_m"],
                     ["floor", ARM]]
    assert set(res["not_buildable"]) == {ARM}
    rec = res["not_buildable"][ARM]
    assert rec["status"] == n50.n36.NOT_BUILDABLE and rec["seed"] == 2
    assert rec["message"] == "ValueError: <task>: anchor_moves='damped_matched' needs obs['group'] in both brackets"
    scored = {"floor", "fgw", "fgw_partial_grp", "fgw_damp_m"}
    assert all(set(r["scores"]) == scored for r in res["rows"]) and set(res["candidates"]) == scored
    assert set(res["build_records"]) == {"fgw_partial_grp", "fgw_damp_m"}
    # a rescored arm refusing leaves the arm, the other rescored arm, fgw and the floor scored
    _patch_builds(monkeypatch, n50, {"fgw_partial_grp"}, t.key)
    runner, calls = _fake_runner(SCORE)
    res = n50.score_task(t, tmp_path, runner=runner)
    assert set(res["not_buildable"]) == {"fgw_partial_grp"}
    assert all(set(r["scores"]) == {"floor", "fgw", "fgw_damp_m", ARM} for r in res["rows"])
    assert set(res["build_records"]) == {"fgw_damp_m", ARM}
    # every guarded arm failing still leaves the floor and fgw scored
    _patch_builds(monkeypatch, n50, set(n50.GUARDED), t.key)
    runner, calls = _fake_runner(SCORE)
    res = n50.score_task(t, tmp_path, runner=runner)
    assert set(res["not_buildable"]) == set(n50.GUARDED)
    assert all(set(r["scores"]) == {"floor", "fgw"} for r in res["rows"])


def test_a_failure_outside_a_guarded_arms_build_stops_the_task(tmp_path, monkeypatch, n50):
    def runner(t, cands, seeds, data_dir):
        raise OSError("the data are not there")

    with pytest.raises(OSError, match="not there"):
        n50.score_task(n50.nt.all_tasks()["x2"], tmp_path, runner=runner)
    t = n50.nt.all_tasks()["x1"]
    _patch_builds(monkeypatch, n50, set(), t.key)

    def bad_fgw(ctx, build_record=None):
        raise ValueError("fgw fails")

    monkeypatch.setitem(harness.CANDIDATES, "fgw", dataclasses.replace(harness.CANDIDATES["fgw"], build=bad_fgw))
    runner, _ = _fake_runner(SCORE)
    with pytest.raises(ValueError, match="fgw fails"):
        n50.score_task(t, tmp_path, runner=runner)


def test_the_runner_writes_n46s_records_from_real_builds(tmp_path, emb_expr, n50):
    c = _ctx(emb_expr, _pair_stages(emb_expr))

    def ctx_of(seed):
        return dataclasses.replace(c, seed=seed)

    runner, _ = _fake_runner(SCORE, ctx_of=ctx_of)
    res = n50.score_task(n50.nt.all_tasks()["e75"], tmp_path, runner=runner, seeds=(0, 1))
    br = res["build_records"]
    assert set(br) == set(n50.GUARDED) and res["not_buildable"] == {}
    for arm in br:
        for seed in ("0", "1"):
            r = br[arm][seed]
            assert "shift" not in r and set(r["movement_ratio"]) == {"anchor", "cell"}
            assert all(len(v) == 2 and all(x is not None for x in v) for v in r["movement_ratio"].values())
            assert ("factors" in r) == (arm in DAMPED)
    for seed in ("0", "1"):
        j, g, m = br[ARM][seed], br["fgw_partial_grp"][seed], br["fgw_damp_m"][seed]
        assert set(g) == {"partial", "movement_ratio"}
        assert set(j) == set(m) == {"anchor_move", "movement_ratio", "factors"}
        np.testing.assert_allclose(j["movement_ratio"]["cell"], j["factors"], rtol=1e-12)
        np.testing.assert_allclose(j["movement_ratio"]["cell"], g["movement_ratio"]["cell"], rtol=1e-9)
        np.testing.assert_allclose(m["movement_ratio"]["cell"], m["factors"], rtol=1e-12)
        assert j["factors"] == j["anchor_move"]["factors"] and j["factors"] != m["factors"]
        assert j["anchor_move"]["s"] == g["partial"]["s"]
    json.dumps(br)  # floats and None only
