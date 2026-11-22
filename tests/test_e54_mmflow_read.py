"""Row E54 (the multi-stage flow's kill test, as specified): the read script experiments/coords/e54_mmflow_read.py and
the harness arm vec.harness.mmflow_candidate, on synthetic stages and hand-made per-seed records only. Nothing is read
from data/, no fit is made, no practice task is run and nothing is scored against a released stage.

Pinned: the folds (each task's exclusions by the rule equal the specification's table; seven task fits and the gate's
one more); the arm keeps C0's coordinates, count, bracket quotas and slots and changes only which real cells sit in
them (identity states give C0's file bit for bit); a refusal is recorded and C0's file stands in, any other error
propagates; the gate's rule and its counting; the pass line's clauses (a) to (e) with their boundaries; every void
condition; a reading is never overwritten."""

from __future__ import annotations

import copy
import sys
from dataclasses import replace
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from vec import harness
from vec import t2_mmflow as mf
from vec.boards import Board

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "experiments" / "coords"))
import e54_mmflow_read as e54  # noqa: E402

C = "c" * 40
TASKS = list(e54.nt.all_tasks())
SEEDS = list(range(10))
GROUPS = ("expression_change", "cell_state", "shape_growth", "local_organisation")
METRICS = ("de_score", "de_direction", "mmd_u", "variogram", "d2_shape", "occupancy_dice", "scale_log_ratio",
           "neighborhood_mmd")


# ------------------------------------------------------------------------------------------------ folds
def test_the_folds_are_the_specifications_seven_and_the_gates_one_more():
    assert e54.exclusion_problems() == []
    tf, gf = e54.task_folds(), e54.gate_folds()
    assert len(set(tf.values())) == 7 and len(e54.folds()) == 8
    assert set(gf.values()) - set(tf.values()) == {"x825_late"}
    assert tf["late"] == gf["T2/E8.0.h5ad"] and tf["e75"] == gf["T2/E7.5.h5ad"] and tf["x3"] == gf["T2/E8.5.h5ad"]
    for name, t in e54.nt.all_tasks().items():  # a fit never holds the target, and always holds both brackets
        left_out = e54.folds()[tf[name]]
        assert t.target in left_out and t.inputs["before"] not in left_out and t.inputs["after"] not in left_out
    assert e54.DAYS["T2/E8.25_late.h5ad"] == harness.STAGE_DAYS["T2/E8.25_late.h5ad"]


def test_the_gate_stages_and_rule():
    assert e54.GATE_STAGES == ("T2/E7.25.h5ad", "T2/E7.5.h5ad", "T2/E8.0.h5ad", "T2/E8.25_late.h5ad",
                               "T2/E8.5.h5ad")
    assert (e54.GATE_MUST, e54.GATE_NEED) == ("T2/E7.5.h5ad", 4)
    st = {f: {"beats": True, "flow_mean": 1.0, "mixture_mean": 2.0} for f in e54.GATE_STAGES}
    g = e54.gate_reading(st, C)
    assert g["pass"] and g["counted"] == ["E54a", "E54b", "E54c"] and "passes" in e54.gate_text(g)
    st["T2/E7.5.h5ad"]["beats"] = False
    g = e54.gate_reading(st, C)
    assert not g["pass"] and g["counted"] == ["E54a"] and "do not run" in g["outcome"]


def test_has_row_reads_the_ledger(tmp_path):
    p = tmp_path / "TRIES.md"
    p.write_text("| # | Board |\n|---|---|\n| E53 | embryo |\n")
    assert not e54.has_row(p)
    p.write_text("| # | Board |\n|---|---|\n| E54 | embryo |\n")
    assert e54.has_row(p)


def test_the_steps_refuse_until_the_register_names_torch_and_torchdiffeq(tmp_path, monkeypatch):
    """The review's lesser point 1: the X005 amendment is checked by code, in every step's refusal path."""
    p = tmp_path / "EXTERNAL_DATA.md"
    row = "| X005 | Open-source Python libraries: POT 0.9.7, numpy 2.5, scipy 1.18{} | MIT | code |\n"
    assert e54.register_problems(tmp_path / "missing.md") != []
    p.write_text(row.format(""))
    assert len(e54.register_problems(p)) == 2
    p.write_text(row.format("; torchdiffeq 0.2.5"))  # torchdiffeq alone does not name torch
    assert e54.register_problems(p) == ["the register's X005 row does not name torch: its amendment is committed "
                                        "before the run"]
    p.write_text(row.format("; pytorch-lightning, torch-geometric, torchdiffeq 0.2.5"))
    assert len(e54.register_problems(p)) == 1  # neither names torch itself
    p.write_text(row.format("; torch 2.14.1 and torchdiffeq 0.2.5 (row E54's flow; torchcfm installed, unused)"))
    assert e54.register_problems(p) == []
    p.write_text(row.format("; torch, torchdiffeq") + row.format(""))
    assert e54.register_problems(p) != []  # two X005 rows: refused
    monkeypatch.setattr(e54, "REGISTER", tmp_path / "missing.md")
    monkeypatch.setattr(e54, "_own_code", lambda: None)
    monkeypatch.setattr(e54, "has_row", lambda path=None: True)
    monkeypatch.setattr(e54, "repo_state", lambda repo: type("S", (), {"unclean": ()})())
    monkeypatch.setattr(e54, "commit_problems", lambda c: [])
    with pytest.raises(SystemExit, match="X005"):
        e54._refuse_unless_ready(C)
    p.write_text(row.format("; torch 2.14.1, torchdiffeq 0.2.5"))
    monkeypatch.setattr(e54, "REGISTER", p)
    e54._refuse_unless_ready(C)  # nothing else refuses here


def test_a_reading_is_never_overwritten(tmp_path):
    f = tmp_path / "r" / "reading.txt"
    e54._write_new(f, "one")
    with pytest.raises(FileExistsError):
        e54._write_new(f, "two")
    assert f.read_text() == "one"


def test_load_gate(tmp_path):
    assert e54.load_gate(tmp_path, C)[0] is None
    import json

    (tmp_path / "gate.json").write_text(json.dumps({"commit": C, "pass": True}))
    assert e54.load_gate(tmp_path, C)[1] == []
    assert e54.load_gate(tmp_path, "d" * 40)[1] == ["the gate reading is from another commit"]
    (tmp_path / "gate.json").write_text(json.dumps({"commit": C, "pass": False}))
    assert e54.load_gate(tmp_path, C)[1] == ["the gate did not pass: E54b and E54c are not run"]


# ------------------------------------------------------------------------------------------------ the harness arm
GENES = [f"g{i}" for i in range(15)]


def _pool(n, rng, shift, scale, prefix):
    types = rng.choice(2, n, p=[0.4, 0.6])
    means = np.array([np.r_[np.full(8, 3.0), np.full(7, 0.5)], np.r_[np.full(8, 0.5), np.full(7, 3.0)]]) * scale
    counts = rng.poisson(means[types]).astype(np.int32)
    counts[:, 0] += 1  # every cell has an expressed panel gene
    lib = counts.sum(1, keepdims=True) + rng.integers(20, 80, (n, 1))
    X = sparse.csr_matrix(np.log1p(1e4 * counts / lib).astype(np.float32))
    a = ad.AnnData(X=X, obs=pd.DataFrame({"celltype": np.where(types == 0, "EXE", "SOM"),
                                          "group": np.where(types == 0, "Endoderm", "Mesoderm")},
                                         index=[f"{prefix}{i}" for i in range(n)]),
                   var=pd.DataFrame(index=GENES))
    a.layers["counts"] = sparse.csr_matrix(counts)
    centres = np.array([[0, 0, 0], [5, 0, 0]], float)
    a.obsm["spatial_3D"] = (centres[types] + rng.normal(size=(n, 3)) * 1.5 + shift).astype(np.float32)
    return a


@pytest.fixture(scope="module")
def world():
    rng = np.random.default_rng(54)
    before, after = _pool(420, rng, 0.0, 1.0, "b"), _pool(520, rng, 1.0, 1.6, "a")
    board = Board(key="pseudo:synthetic_embryo", task="T2", genes=list(GENES), needs_coords=True, setting="embryo",
                  min_cells=60, max_cells=300)
    task = harness.PseudoTask(key="syn_e54", test_board="T2:embryo:test_interp", task="T2", setting="embryo",
                              kind="interp", target="T2/E7.5.h5ad", reference="T2/E7.25.h5ad",
                              inputs={"before": "T2/E7.25.h5ad", "after": "T2/E8.0.h5ad"},
                              twin_board="T2:embryo:val_interp")
    ctx = harness.BuildContext(task=harness.TaskView.of(task), board=board, inputs={"before": before, "after": after},
                               seed=harness.BUILD_SEED_OFFSET + 3)
    lb, la = rng.normal(size=(before.n_obs, 4)), rng.normal(size=(after.n_obs, 4))
    states = {"before_obs": np.asarray(before.obs_names, str), "after_obs": np.asarray(after.obs_names, str),
              "before_latent": lb, "after_latent": la, "before_snap": lb.copy(), "after_snap": la.copy(),
              "before_anti": lb.copy(), "after_anti": la.copy(), "days": np.array([7.25, 8.0, 7.5])}
    return task, ctx, states, rng


def _c0_res(ctx):
    b, a = ctx.inputs["before"], ctx.inputs["after"]
    g = ctx.board.genes
    Eb, Ea = (harness.dense(x[:, g].X).astype(np.float32) for x in (b, a))
    t = ctx.task
    frac = (t.target_day - t.role_days["before"]) / (t.role_days["after"] - t.role_days["before"])
    return harness.t2_build.interpolate_clouds(np.asarray(b.obsm["spatial_3D"], float), np.asarray(a.obsm["spatial_3D"],
                                               float), Eb, Ea, frac, harness.count_for(ctx.board, b.n_obs), ctx.seed,
                                               target_rms=harness.target_rms(ctx), t_pos=0.0,
                                               positions="floor_exact"), Eb, Ea


def test_identity_states_give_c0s_file_bit_for_bit(world):
    task, ctx, states, _ = world
    c0 = harness.CANDIDATES["fgw_shape0x"].build(ctx)
    facts = {}
    for name in e54.ARMS:
        arm = harness.mmflow_candidate(name, states, facts).build(ctx)
        np.testing.assert_array_equal(arm.X, c0.X)
        np.testing.assert_array_equal(arm.coords, c0.coords)
        rec = facts[ctx.seed]
        for br in ("before", "after"):
            assert rec[br]["library_ratio"] == 1.0 and rec[br]["kish_frac_arm"] == rec[br]["kish_frac_c0"]
            assert rec[br]["shares_arm"] == rec[br]["shares_c0"] and set(rec[br]["shares_c0"]) <= {"EXE", "SOM"}


def test_other_states_change_which_cells_sit_in_c0s_slots_and_nothing_else(world):
    task, ctx, states, rng = world
    pb, pa = rng.permutation(len(states["before_latent"])), rng.permutation(len(states["after_latent"]))
    moved = {**states, "before_snap": states["before_latent"][pb], "after_snap": states["after_latent"][pa],
             "before_anti": states["before_latent"][pb[pb]], "after_anti": states["after_latent"][pa[pa]]}
    c0 = harness.CANDIDATES["fgw_shape0x"].build(ctx)
    res, Eb, Ea = _c0_res(ctx)
    for name, (qb, qa) in (("fgw_mmflow_snap", (pb, pa)), ("fgw_mmflow_anti", (pb[pb], pa[pa]))):
        facts = {}
        arm = harness.mmflow_candidate(name, moved, facts).build(ctx)
        np.testing.assert_array_equal(arm.coords, c0.coords)  # the coordinates, jitter and order of C0's file
        assert arm.n == c0.n
        b = res.stage == 0
        want = np.empty_like(np.asarray(arm.X))
        want[b], want[~b] = Eb[qb[res.rows[b]]], Ea[qa[res.rows[~b]]]
        np.testing.assert_array_equal(arm.X, np.clip(want, 0.0, None))
        assert not np.array_equal(arm.X, c0.X)
        assert facts[ctx.seed]["before"]["library_ratio"] is not None


def test_the_arm_checks_its_states_belong_to_the_task(world):
    task, ctx, states, _ = world
    with pytest.raises(ValueError, match="lack"):
        harness.mmflow_candidate("fgw_mmflow_snap", {k: v for k, v in states.items() if k != "after_anti"})
    bad = {**states, "days": np.array([7.25, 8.0, 7.75])}
    with pytest.raises(ValueError, match="days"):
        harness.mmflow_candidate("fgw_mmflow_snap", bad).build(ctx)
    other = {**states, "before_obs": states["before_obs"][::-1]}
    with pytest.raises(ValueError, match="not the input's"):
        harness.mmflow_candidate("fgw_mmflow_snap", other).build(ctx)
    c = harness.mmflow_candidate("fgw_mmflow_anti", states)
    assert c.applies(task) and c.params["direction"] == "anti" and "fgw_mmflow_snap" not in harness.CANDIDATES
    heart = replace(task, key="syn_heart", kind="extrap", setting="heart", target="T2/E10.5.h5ad",
                    reference="T2/E9.5.h5ad", inputs={"ref": "T2/E9.5.h5ad"}, twin_board="T2:heart:val_extrap")
    assert not c.applies(heart)


def _seen(split="pseudo"):
    board = Board(key="pseudo:x", task="T2", genes=["g0"], needs_coords=True, split=split)
    return {"board": board, "stage": {(e54.C0, 2): "C0 file"}, "refused": {}}


def test_a_refusal_is_recorded_and_c0s_file_stands_in():
    def refuse(ctx):
        raise mf.MmflowRefused("one-to-one retrieval needs at least as many pool cells as predicted states")

    base = harness.Candidate("fgw_mmflow_snap", frozenset({"T2"}), (), refuse)
    seen = _seen()
    ctx = harness.BuildContext(task=None, board=None, inputs={}, seed=harness.BUILD_SEED_OFFSET + 2)
    assert e54.guarded(base, seen).build(ctx) == "C0 file"
    assert seen["refused"][("fgw_mmflow_snap", 2)].startswith("one-to-one")

    def other(ctx):
        raise ValueError("something else")

    with pytest.raises(ValueError, match="something else"):
        e54.guarded(replace(base, build=other), _seen()).build(ctx)
    with pytest.raises(RuntimeError, match="test board"):
        e54.guarded(base, _seen("test")).build(ctx)


# ------------------------------------------------------------------------------------------------ the reading
def _arm_facts(ratio=1.0, jit=False):
    br = {"library_ratio": ratio, "kish_frac_c0": 0.3, "kish_frac_arm": 0.3, "shares_c0": {"EXE": 0.5, "SOM": 0.5},
          "shares_arm": {"EXE": 0.4, "SOM": 0.6}}
    return {"n_jittered_any": jit, "before": dict(br), "after": dict(br)}


def _rec(task, seed, arm=52.0, anti=51.5, cs_arm=0.01, commit=C):
    names = (harness.FLOOR, e54.C0, e54.ARM, e54.ANTI)
    scores = {harness.FLOOR: 50.0, e54.C0: 51.0, e54.ARM: arm, e54.ANTI: anti}
    skills = {n: dict.fromkeys(METRICS, 0.4) for n in names}
    groups = {n: dict.fromkeys(GROUPS, 0.4) for n in names}
    groups[e54.ARM]["cell_state"] += cs_arm
    row = {"seed": seed, "scores": scores, "skills": skills, "groups": groups, "n_cells": dict.fromkeys(names, 200)}
    res = {"provenance": {"commit": commit, "tree_clean": True}, "build": "main", "protocol": "bake",
           "board": {"split": "pseudo"}, "rows": [row]}
    return {"task": task, "seed": seed, "run_commit": commit, "builds": {"main": res},
            "facts": {"board_split": "pseudo", "refused": dict.fromkeys(e54.ARMS),
                      "arm": {n: _arm_facts() for n in e54.ARMS}}}


def _recs(gain=None, anti=None, **kw):
    gain = gain or {}
    anti = anti or {}
    return {t: [_rec(t, s, arm=51.0 + gain.get(t, 1.0), anti=51.0 + anti.get(t, 0.0), **kw) for s in SEEDS]
            for t in TASKS}


GATE = {"pass": True, "commit": C, "beats": 5, "must": "T2/E7.5.h5ad", "stages": dict.fromkeys(e54.GATE_STAGES, {})}


def _read(recs, gate=GATE, why=()):
    return e54.reading(recs, gate, list(why), C)


def test_the_pass_line_holds_at_its_boundaries():
    res = _read(_recs())  # every task +1.0 over C0 exactly, the anti at C0
    assert not res["void"] and res["passes"], res["pass_line"]["clauses"]
    assert res["counted"] == ["E54a", "E54b", "E54c"] and "passes" in res["outcome"]
    p = res["pass_line"]
    assert p["merged_board"]["mean"] == 1.0 and p["test_like_board"]["mean"] == 1.0
    assert p["without_e725"]["mean"] == 1.0 and set(p["leave_one_truth_out"].values()) == {1.0}
    assert "holds" in e54.text(res) and res["e75_e725_apart"]["e725"]["arm"] == 1.0
    low = _read(_recs(gain=dict.fromkeys(TASKS, 0.999)))
    assert not low["passes"] and not low["pass_line"]["clauses"]["a"] and not low["pass_line"]["clauses"]["b"]


def test_a_fails_when_one_truth_carries_the_gain():
    gain = {t: 0.2 for t in TASKS}
    gain["e75"] = 10.0
    c = _read(_recs(gain=gain))["pass_line"]["clauses"]
    assert not c["a"] and c["e"]


def test_a_needs_four_truths_positive():
    gain = dict.fromkeys(TASKS, 3.0)
    gain.update(e75=-0.1, e725=-0.1)
    p = _read(_recs(gain=gain))["pass_line"]
    assert p["merged_board"]["mean"] > 1.0 and min(p["leave_one_truth_out"].values()) >= 0.5 and not p["clauses"]["a"]


def test_b_is_the_test_like_tasks_alone():
    gain = dict.fromkeys(TASKS, 3.0)
    gain.update(dict.fromkeys(e54.TEST_LIKE, 0.5))
    c = _read(_recs(gain=gain))["pass_line"]["clauses"]
    assert c["a"] and not c["b"]


def test_c_needs_the_margin_over_the_placebo_and_eight_seeds():
    res = _read(_recs(gain=dict.fromkeys(TASKS, 2.0), anti=dict.fromkeys(TASKS, 1.5)))
    assert not res["pass_line"]["clauses"]["c"] and res["pass_line"]["clauses"]["a"]
    recs = _recs(gain=dict.fromkeys(TASKS, 2.0), anti=dict.fromkeys(TASKS, 0.0))
    for t in TASKS:
        for s in (0, 1, 2):
            recs[t][s] = _rec(t, s, arm=53.0, anti=53.1)
    p = _read(recs)["pass_line"]
    assert p["arm_minus_anti"]["mean"] >= 0.7 and p["arm_minus_anti"]["ahead"] == 7 and not p["clauses"]["c"]


def test_d_is_the_cell_state_group():
    res = _read(_recs(cs_arm=-0.001))
    assert not res["pass_line"]["clauses"]["d"] and not res["passes"]
    assert _read(_recs(cs_arm=0.0))["pass_line"]["clauses"]["d"] is False  # above 0, not at it


def test_e_leaves_e725_out():
    gain = dict.fromkeys(TASKS, 0.6)
    gain["e725"] = 10.0
    c = _read(_recs(gain=gain))["pass_line"]["clauses"]
    assert not c["e"]


def _void(recs, **kw):
    res = _read(recs, **kw)
    assert res["void"] and "passes" not in res
    return " ".join(res["why"])


@pytest.mark.parametrize("arm", ["fgw_mmflow_snap", "fgw_mmflow_anti"])
@pytest.mark.parametrize("bracket", ["before", "after"])
def test_a_library_shift_voids(arm, bracket):
    recs = _recs()
    recs["x2"][4]["facts"]["arm"][arm][bracket]["library_ratio"] = 1.051
    assert "median library" in _void(recs)
    recs["x2"][4]["facts"]["arm"][arm][bracket]["library_ratio"] = 0.949
    assert "median library" in _void(recs)
    recs["x2"][4]["facts"]["arm"][arm][bracket]["library_ratio"] = 1.049
    assert not _read(recs)["void"]
    recs["x2"][4]["facts"]["arm"][arm][bracket]["library_ratio"] = None
    assert "unreadable" in _void(recs)


@pytest.mark.parametrize("term", ["d2_shape", "occupancy_dice", "scale_log_ratio"])
def test_a_shape_term_off_c0s_voids_unless_a_file_was_jittered(term):
    recs = _recs()
    recs["half"][0]["builds"]["main"]["rows"][0]["skills"][e54.ANTI][term] += 2e-9
    assert "shape terms" in _void(recs)
    recs["half"][0]["facts"]["arm"][e54.ANTI]["n_jittered_any"] = True
    assert not _read(recs)["void"]
    recs = _recs()
    recs["half"][0]["builds"]["main"]["rows"][0]["skills"][e54.ARM][term] += 5e-10
    assert not _read(recs)["void"]


def test_the_other_void_conditions():
    recs = _recs()
    recs["new1"][3]["builds"]["main"]["rows"][0]["n_cells"][e54.ARM] = 199
    assert "count differs" in _void(recs)
    recs = _recs()
    recs["e75"][9]["facts"]["refused"][e54.ARM] = "one-to-one retrieval needs ..."
    assert "missing arm" in _void(recs)
    recs = _recs()
    del recs["x3"][5]["builds"]["main"]["rows"][0]["scores"][e54.ANTI]
    assert "missing" in _void(recs)
    recs = _recs()
    recs["late"] = recs["late"][:9]
    assert "seeds are not" in _void(recs)
    recs = _recs()
    recs["late"][2]["run_commit"] = "d" * 40
    assert "another run commit" in _void(recs)
    recs = _recs()
    recs["late"][2]["builds"]["main"]["provenance"]["tree_clean"] = False
    assert "unclean" in _void(recs)
    recs = _recs()
    recs["x1"][0]["builds"]["main"]["board"]["split"] = "test"
    assert "test board" in _void(recs)
    recs = _recs()
    recs["x1"][0]["builds"]["main"]["protocol"] = "live"
    assert "bake" in _void(recs)
    assert "no passed gate" in _void(_recs(), gate=None)
    assert "no passed gate" in _void(_recs(), gate={**GATE, "pass": False})
    assert "a repro miss" in _void(_recs(), why=["a repro miss"])
    assert "VOID" in e54.text(_read(_recs(), gate=None))


def test_reported_tables_decide_nothing():
    recs = _recs()
    base = _read(recs)
    shifted = copy.deepcopy(recs)
    for t in TASKS:
        for r in shifted[t]:
            for n in e54.ARMS:
                for b in ("before", "after"):
                    r["facts"]["arm"][n][b].update(kish_frac_arm=0.01, shares_arm={"EXE": 1.0})
    other = _read(shifted)
    assert other["passes"] == base["passes"] and other["pass_line"] == base["pass_line"]
    rep = other["reported"]["half"][e54.ARM]["brackets"]["before"]
    assert rep["kish_frac_arm"] == pytest.approx(0.01) and rep["shares_arm"] == {"EXE": pytest.approx(1.0)}
