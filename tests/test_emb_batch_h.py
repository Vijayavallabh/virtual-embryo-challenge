"""The embryo batch H (campaign branch campaign/emb-h): `fgw_partial_grp`, fgw_partial with its transported share s
taken from the overlap of the brackets' anchors' obs["group"] distributions in place of obs["celltype"], and the runner
experiments/harness/n47_embryo_h_arms.py. Synthetic stages only (tests/test_emb_expr.py's); no data file, committed
result or scorer is read.

Tolerances stated here: the partial coupling's total mass against s and its marginals against the uniform masses,
1e-9 absolute, as in tests/test_emb_batch_f.py (POT's exact network simplex)."""

from __future__ import annotations

import dataclasses
import importlib.util
import json
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np
import pytest

from vec import harness
from vec import t2_build as tb

REPO = Path(__file__).resolve().parents[1]
HARNESS = REPO / "experiments" / "harness"
ARM = "fgw_partial_grp"
# sha256 of fgw_partial's (X float32, coords float32) and its build record on `_pair_stages()`, computed with the
# branch's base campaign/emb-f at fcc1f6b before this batch's edit (OPENBLAS/OMP/MKL threads 4)
FGW_PARTIAL_AT_BASE = "659fbae5513407f2f40ce6730e0418c3678304cb20eca8057c6f89da2fe6e711"
FGW_PARTIAL_RECORD_AT_BASE = {"partial": {"s": 0.6666666666666666, "transported_mass": 0.6666666666666663,
                                          "no_mass_share": [0.3145833333333333, 0.28939393939393937]}}


@pytest.fixture(scope="module")
def emb_expr():
    """The synthetic embryo stages and contexts of tests/test_emb_expr.py (no data file)."""
    spec = importlib.util.spec_from_file_location("test_emb_expr_helpers_h", REPO / "tests" / "test_emb_expr.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _rename_som(s: dict) -> dict:
    """The after bracket's SOM cells relabelled SOM2 (cell type only; their group stays ParaxialMesoderm): the cell
    types differ between the brackets while the groups agree."""
    a = s["after"]
    a.obs["celltype"] = a.obs["celltype"].astype(str).replace({"SOM": "SOM2"}).to_numpy()
    return s


def _pair_stages(emb_expr) -> dict:
    """The synthetic pair: cell types differ (SOM / SOM2) and the groups agree in part (the after bracket's first 50
    cells, all EXE-Endoderm, moved from group Endoderm to Ectoderm)."""
    s = _rename_som(emb_expr.stages())
    g = s["after"].obs["group"].astype(str).to_numpy().copy()
    g[:50] = "Ectoderm"
    s["after"].obs["group"] = g
    return s


def _ctx(emb_expr, s, seed=None):
    kw = {} if seed is None else {"seed": seed}
    return emb_expr.ctx(inputs=s, roles=("before", "after"), **kw)


class _Spy:
    """Records the partial solver's mass and coupling, the partners and the rotation's coupling of one build."""

    def __init__(self, monkeypatch):
        self.solver, self.partners, self.rot, self.interp = [], [], [], []
        real_p, real_pp = tb.partial_fgw_coupling, tb.partial_partners
        real_rot, real_i = tb.coupling_rotation, tb.interpolate_clouds

        def p(Fb, Fa, Zb, Za, alpha, m, it):
            T = real_p(Fb, Fa, Zb, Za, alpha, m, it)
            self.solver.append((m, T))
            return T

        def pp(T, own, other):
            out = real_pp(T, own, other)
            self.partners.append((np.array(T), np.array(own), np.array(other), out))
            return out

        def rot(a, b, T):
            self.rot.append(np.array(T))
            return real_rot(a, b, T)

        def interp(*a, **kw):
            res = real_i(*a, **kw)
            self.interp.append((kw, res))
            return res

        monkeypatch.setattr(tb, "partial_fgw_coupling", p)
        monkeypatch.setattr(tb, "partial_partners", pp)
        monkeypatch.setattr(tb, "coupling_rotation", rot)
        monkeypatch.setattr(tb, "interpolate_clouds", interp)


# ------------------------------------------------------------------------------------------ 1. fgw_partial unchanged
def test_fgw_partial_is_what_it_was_at_the_base_bit_for_bit(emb_expr):
    c = _ctx(emb_expr, _pair_stages(emb_expr))
    rec: dict = {}
    st = harness._fgw_partial(c, build_record=rec)
    assert emb_expr.digest(st) == FGW_PARTIAL_AT_BASE
    # campaign/emb-j: since the merge of campaign/emb-g every record also holds "shift" (fix brief M2); the rest is
    # the base's
    assert set(rec.pop("shift")) == {"anchor", "cell"}
    assert json.dumps(rec, sort_keys=True) == json.dumps(FGW_PARTIAL_RECORD_AT_BASE, sort_keys=True)
    assert emb_expr.digest(harness.CANDIDATES["fgw_partial"].build(c)) == FGW_PARTIAL_AT_BASE


def test_fgw_partial_still_passes_celltype(emb_expr, monkeypatch):
    seen = []
    real = tb.interpolate_clouds
    monkeypatch.setattr(tb, "interpolate_clouds", lambda *a, **kw: seen.append(kw) or real(*a, **kw))
    s = _pair_stages(emb_expr)
    c = _ctx(emb_expr, s)
    harness.CANDIDATES["fgw_partial"].build(c)
    harness.CANDIDATES[ARM].build(c)
    (kp, kg) = seen
    assert kp["coupling"] == kg["coupling"] == "partial"
    for side, i in (("before", 0), ("after", 1)):
        assert np.array_equal(np.asarray(kp["coupling_labels"][i]), np.asarray(s[side].obs["celltype"]))
        assert np.array_equal(np.asarray(kg["coupling_labels"][i]), np.asarray(s[side].obs["group"]))
    # every other argument of the call is fgw_partial's
    assert {k: v for k, v in kp.items() if k != "coupling_labels"}.keys() == \
        {k: v for k, v in kg.items() if k != "coupling_labels"}.keys()
    for k in kp:
        if k != "coupling_labels":
            assert json.dumps(kp[k], default=str) == json.dumps(kg[k], default=str), k


# ------------------------------------------------------------------------------------------ 2. s by hand
def test_the_group_overlap_by_hand_on_a_small_case():
    before = ["Endoderm", "Endoderm", "Mesoderm", "Ectoderm"]
    after = ["Endoderm", "Mesoderm", "Mesoderm", "Mesoderm", "Ectoderm"]
    # p_b: Endo 1/2, Meso 1/4, Ecto 1/4; p_a: Endo 1/5, Meso 3/5, Ecto 1/5; min: 1/5 + 1/4 + 1/5
    assert tb.label_overlap(before, after) == Fraction(13, 20)


def test_the_arms_s_is_the_overlap_of_the_anchors_groups_and_differs_from_fgw_partials(emb_expr, monkeypatch):
    s = _pair_stages(emb_expr)
    c = _ctx(emb_expr, s)
    spy = _Spy(monkeypatch)
    rec: dict = {}
    harness._fgw_partial_grp(c, build_record=rec)
    # by hand: every cell is an anchor here (480 and 660 cells, n_anchor 2000); before: three groups of 160 cells
    # (1/3 each); after: Endoderm 170, Ectoderm 50, ParaxialMesoderm 220, Mesoderm 220 of 660 cells;
    # min: 170/660 + 220/660 + 220/660 = 61/66
    assert s["before"].n_obs <= tb.INTERP_N_ANCHOR and s["after"].n_obs <= tb.INTERP_N_ANCHOR
    assert rec["partial"]["s"] == float(Fraction(61, 66))
    assert Fraction(61, 66) == tb.label_overlap(s["before"].obs["group"], s["after"].obs["group"])
    (m, _), = spy.solver
    assert m == float(Fraction(61, 66))
    # fgw_partial's s on the same pair: the celltype overlap, 2/3 by hand (SOM against SOM2)
    rec_p: dict = {}
    harness._fgw_partial(c, build_record=rec_p)
    assert rec_p["partial"]["s"] == float(Fraction(2, 3)) != rec["partial"]["s"]


def test_with_anchors_drawn_the_arms_s_is_the_overlap_of_the_anchors_groups(emb_expr, monkeypatch):
    """The share is read on the anchors, not the whole brackets: with fewer anchors than cells, s is label_overlap of
    the drawn anchors' groups (the anchors interpolate_clouds itself returns with with_coupling)."""
    s = _pair_stages(emb_expr)
    real = tb.interpolate_clouds
    seen = []

    def few(*a, **kw):
        kw = {**kw, "n_anchor": 100, "with_coupling": True}
        res = real(*a, **kw)
        seen.append(res)
        return res

    monkeypatch.setattr(tb, "interpolate_clouds", few)
    rec: dict = {}
    harness._fgw_partial_grp(_ctx(emb_expr, s), build_record=rec)
    (res,) = seen
    ib, ia = res.coupling["anchors_before"], res.coupling["anchors_after"]
    gb = np.asarray(s["before"].obs["group"]).astype(str)
    ga = np.asarray(s["after"].obs["group"]).astype(str)
    assert rec["partial"]["s"] == float(tb.label_overlap(gb[ib], ga[ia]))


# ------------------------------------------------------------------------------------------ 3. the coupling
def test_the_coupling_carries_mass_s_within_the_uniform_masses(emb_expr, monkeypatch):
    c = _ctx(emb_expr, _pair_stages(emb_expr))
    spy = _Spy(monkeypatch)
    rec: dict = {}
    harness._fgw_partial_grp(c, build_record=rec)
    (m, T), = spy.solver
    s = rec["partial"]["s"]
    n_b, n_a = T.shape
    assert m == s and 0.0 < s < 1.0
    assert abs(T.sum() - s) <= 1e-9 and rec["partial"]["transported_mass"] == float(T.sum())
    assert (T >= 0).all() and (T.sum(1) <= 1.0 / n_b + 1e-9).all() and (T.sum(0) <= 1.0 / n_a + 1e-9).all()
    assert rec["partial"]["no_mass_share"] == [float((T.sum(1) <= 0).mean()), float((T.sum(0) <= 0).mean())]
    (R_in,) = spy.rot  # the rotation is coupling_rotation on this coupling
    assert np.array_equal(R_in, T)


def test_an_anchor_without_mass_stays_and_a_matched_one_goes_to_its_row_barycentre(emb_expr, monkeypatch):
    c = _ctx(emb_expr, _pair_stages(emb_expr))
    spy = _Spy(monkeypatch)
    rec: dict = {}
    harness._fgw_partial_grp(c, build_record=rec)
    assert len(spy.partners) == 2 and max(rec["partial"]["no_mass_share"]) > 0.0
    (m, T), = spy.solver
    assert np.array_equal(spy.partners[0][0], T) and np.array_equal(spy.partners[1][0], T.T)
    for Tk, own, other, out in spy.partners:
        w = Tk.sum(1)
        none = w <= 0.0
        assert none.any()
        assert np.array_equal(out[none] - own[none], np.zeros_like(own[none]))  # displacement exactly 0
        np.testing.assert_allclose(out[~none], (Tk[~none] @ other) / w[~none, None], rtol=1e-12, atol=1e-12)


def test_the_drawn_cells_and_expression_rows_are_fgws(emb_expr, monkeypatch):
    c = _ctx(emb_expr, _pair_stages(emb_expr))
    fgw = harness.CANDIDATES["fgw"].build(c)
    res_f, _, _, _ = emb_expr._interp(c)
    spy = _Spy(monkeypatch)
    arm = harness._fgw_partial_grp(c)
    (_, res), = spy.interp
    assert np.array_equal(res.rows, res_f.rows) and np.array_equal(res.stage, res_f.stage)
    assert np.array_equal(arm.X, fgw.X) and arm.n == fgw.n
    assert not np.array_equal(arm.coords, fgw.coords)  # the idea is not idle here
    part = harness._fgw_partial(c)
    assert not np.array_equal(arm.coords, part.coords)  # and the arm is not fgw_partial


# ------------------------------------------------------------------------------------------ 4. idle: s exactly 1
def test_with_groups_that_overlap_fully_the_arm_is_fgw_bit_for_bit(emb_expr, monkeypatch):
    c = _ctx(emb_expr, _rename_som(emb_expr.stages()))  # cell types differ, groups in equal shares
    fgw = harness.CANDIDATES["fgw"].build(c)
    part = harness.CANDIDATES["fgw_partial"].build(c)
    assert emb_expr.digest(part) != emb_expr.digest(fgw)  # fgw_partial is not idle on this pair

    def never(*a, **kw):
        raise AssertionError("the partial solver is not called at s = 1")

    monkeypatch.setattr(tb, "partial_fgw_coupling", never)
    rec: dict = {}
    arm = harness._fgw_partial_grp(c, build_record=rec)
    assert rec["partial"] == {"s": 1.0, "transported_mass": None, "no_mass_share": [0.0, 0.0]}
    assert emb_expr.digest(arm) == emb_expr.digest(fgw)


# ------------------------------------------------------------------------------------------ 5. refusal
@pytest.mark.parametrize("side", ["before", "after"])
def test_a_bracket_without_group_labels_refuses(emb_expr, side):
    s = _pair_stages(emb_expr)
    del s[side].obs["group"]
    c = _ctx(emb_expr, s)
    with pytest.raises(ValueError, match="group"):
        harness.CANDIDATES[ARM].build(c)
    harness.CANDIDATES["fgw_partial"].build(c)  # fgw_partial reads celltype and still builds


# ------------------------------------------------------------------------------------------ 6. reproducible
def test_same_seed_same_output(emb_expr):
    c = _ctx(emb_expr, _pair_stages(emb_expr))
    a, b = harness._fgw_partial_grp(c), harness.CANDIDATES[ARM].build(c)
    assert emb_expr.digest(a) == emb_expr.digest(b)
    other = harness._fgw_partial_grp(dataclasses.replace(c, seed=c.seed + 1))
    assert emb_expr.digest(other) != emb_expr.digest(a)


# ------------------------------------------------------------------------------------------ 7. registration
def test_the_arm_is_registered_for_t2_with_fgw_partials_requirements(monkeypatch):
    seen = {}
    monkeypatch.setattr(harness, "_fgw", lambda ctx, **kw: seen.setdefault(ctx, kw))
    c = harness.CANDIDATES[ARM]
    assert c.name == ARM and c.tasks == harness.CANDIDATES["fgw_partial"].tasks == frozenset({"T2"})
    assert c.requires is harness.CANDIDATES["fgw_partial"].requires
    assert c.build("a") == {"coupling": "partial", "build_record": None, "coupling_label": "group"}
    assert harness.CANDIDATES["fgw_partial"].build("b") == {"coupling": "partial", "build_record": None}
    p, q = c.params, harness.CANDIDATES["fgw_partial"].params
    assert p["label"] == "group" and q["label"] == "celltype"
    for k in ("alpha", "k", "num_iter_max", "scope"):
        assert p[k] == q[k]


# ================================================================================ the runner
@pytest.fixture(scope="module")
def n47():
    sys.path.insert(0, str(HARNESS))
    spec = importlib.util.spec_from_file_location("n47_embryo_h_arms", HARNESS / "n47_embryo_h_arms.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_the_runner_scores_the_floor_fgw_fgw_partial_and_the_arm_on_the_nine_tasks(n47):
    assert n47.NEW == (ARM,) and n47.ARMS == (harness.FLOOR, "fgw", "fgw_partial", ARM)
    tasks = n47.nt.all_tasks()
    assert set(tasks) == {"half", "late", "new1", "new2", "e75", "e725", "x1", "x2", "x3"}
    for name in n47.ARMS[1:]:
        assert all(harness.CANDIDATES[name].applies(t) for t in tasks.values())


def _patch_builds(monkeypatch, n47, fail_on: set[str], key: str, seed_fails: int = 0):
    """fgw_partial's and the arm's builds replaced (each writes a record when given one): the arms in fail_on raise the
    builder's refusal from seed_fails."""
    cands = harness.CANDIDATES

    def ok(ctx, build_record=None):
        if build_record is not None:
            build_record["seen"] = 1
        return "built"

    def bad(ctx, build_record=None):
        if ctx.seed - harness.BUILD_SEED_OFFSET >= seed_fails:
            raise ValueError(f"{key}: the partial coupling needs obs['group'] in both brackets")
        if build_record is not None:
            build_record["seen"] = 1
        return "built"

    for x in n47.GUARDED:
        monkeypatch.setitem(cands, x, dataclasses.replace(cands[x], build=bad if x in fail_on else ok))
    monkeypatch.setitem(cands, "fgw", dataclasses.replace(cands["fgw"], build=lambda ctx: "built"))


def _fake_runner(score: dict[str, float]):
    """run_task's shape without data or scorer: each candidate is built (its build may raise) and given a constant."""
    calls = []

    def runner(t, cands, seeds, data_dir):
        calls.append([c.name for c in cands])
        rows = []
        for s in seeds:
            row = {"seed": s, "scores": {"floor": 50.0}, "skills": {"floor": {}}, "metrics": {"floor": {}},
                   "n_cells": {"floor": 1}, "groups": {"floor": {}}, "reads": {}}
            for c in cands:
                if c.name == "floor":
                    continue
                c.build(dataclasses.make_dataclass("Ctx", ["seed"])(seed=harness.BUILD_SEED_OFFSET + s))
                row["scores"][c.name] = score[c.name]
                for k in ("skills", "metrics", "groups"):
                    row[k][c.name] = {}
                row["n_cells"][c.name], row["reads"][c.name] = 1, []
            rows.append(row)
        return {"task": {"key": t.key, "inputs": dict(t.inputs)}, "rows": rows,
                "candidates": [c.name for c in cands], "params": {c.name: {} for c in cands},
                "provenance": {"commit": "a" * 40, "tree_clean": True}}

    return runner, calls


def test_the_runners_task_file_loads_through_nine_tasks_with_fgw_and_fgw_partial_under_their_names(
        tmp_path, monkeypatch, n47):
    score = {"fgw": 51.0, "fgw_partial": 52.0, ARM: 53.0}
    t = n47.nt.all_tasks()["e75"]
    _patch_builds(monkeypatch, n47, set(), t.key)
    runner, calls = _fake_runner(score)
    monkeypatch.setattr(n47.harness, "run_task", runner)
    out = tmp_path / "runs"
    n47.main(["run", "--task", "e75", "--out", str(out), "--data", str(tmp_path / "no_data")])
    assert calls == [list(n47.ARMS)]
    (path,) = out.glob("T2e_*.json")
    data = json.loads(path.read_text())
    assert data["task_source"] == "experiments/harness/n47_embryo_h_arms.py" and data["not_buildable"] == {}
    assert data["build_records"] == {ARM: {str(s): {"seen": 1} for s in n47.SEEDS}}  # the new arm's only
    scores, why, meta = n47.nt._load([out])
    assert why == [] and meta == ["runs:aaaaaaa"]
    assert set(scores[t.key]) == set(n47.ARMS) and {"fgw", "fgw_partial"} <= set(scores[t.key])
    assert all(len(v) == 10 for v in scores[t.key].values())
    with pytest.raises(SystemExit, match="exists"):
        n47.main(["run", "--task", "e75", "--out", str(out), "--data", str(tmp_path / "no_data")])
    assert len(calls) == 1
    with pytest.raises(SystemExit, match="unknown task"):
        n47.main(["run", "--task", "nope", "--out", str(out)])


def test_an_arm_the_builder_refuses_is_recorded_and_the_others_are_still_scored(tmp_path, monkeypatch, n47):
    t = n47.nt.all_tasks()["e725"]
    score = {"fgw": 51.0, "fgw_partial": 52.0, ARM: 0.0}
    _patch_builds(monkeypatch, n47, {ARM}, t.key, seed_fails=2)
    runner, calls = _fake_runner(score)
    res = n47.score_task(t, tmp_path, runner=runner)
    assert calls == [list(n47.ARMS), ["floor", "fgw"], ["floor", "fgw_partial"], ["floor", ARM]]
    assert set(res["not_buildable"]) == {ARM}
    rec = res["not_buildable"][ARM]
    assert rec["status"] == n47.n36.NOT_BUILDABLE and rec["seed"] == 2
    assert rec["message"] == "ValueError: <task>: the partial coupling needs obs['group'] in both brackets"
    scored = {"floor", "fgw", "fgw_partial"}
    assert all(set(r["scores"]) == scored for r in res["rows"]) and set(res["candidates"]) == scored
    assert res["build_records"] == {}
    # fgw_partial refusing leaves the arm, fgw and the floor scored, with the arm's records kept
    _patch_builds(monkeypatch, n47, {"fgw_partial"}, t.key)
    runner, calls = _fake_runner(score)
    res = n47.score_task(t, tmp_path, runner=runner)
    assert set(res["not_buildable"]) == {"fgw_partial"}
    assert all(set(r["scores"]) == {"floor", "fgw", ARM} for r in res["rows"])
    assert set(res["build_records"]) == {ARM}
    # both failing still leaves the floor and fgw scored
    _patch_builds(monkeypatch, n47, set(n47.GUARDED), t.key)
    runner, calls = _fake_runner(score)
    res = n47.score_task(t, tmp_path, runner=runner)
    assert set(res["not_buildable"]) == set(n47.GUARDED)
    assert all(set(r["scores"]) == {"floor", "fgw"} for r in res["rows"])


def test_a_failure_outside_a_guarded_arms_build_stops_the_task(tmp_path, monkeypatch, n47):
    def runner(t, cands, seeds, data_dir):
        raise OSError("the data are not there")

    with pytest.raises(OSError, match="not there"):
        n47.score_task(n47.nt.all_tasks()["x2"], tmp_path, runner=runner)
    t = n47.nt.all_tasks()["x1"]
    _patch_builds(monkeypatch, n47, set(), t.key)

    def bad_fgw(ctx):
        raise ValueError("fgw fails")

    monkeypatch.setitem(harness.CANDIDATES, "fgw", dataclasses.replace(harness.CANDIDATES["fgw"], build=bad_fgw))
    runner, _ = _fake_runner({"fgw": 1.0, "fgw_partial": 1.0, ARM: 1.0})
    with pytest.raises(ValueError, match="fgw fails"):
        n47.score_task(t, tmp_path, runner=runner)
