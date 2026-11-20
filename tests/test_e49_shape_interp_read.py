"""Row E49's run script (experiments/coords/e49_shape_interp_read.py) and its launcher, on synthetic arrays, hand-made
result records and the committed E47 score files only. Nothing is read from data/, nothing is built from a stage and
nothing is scored.

Pinned: the step-0 drop rule; the wrapped arms (a refusal is recorded and C0's file stands in its place, any other
error propagates); the pre-run gate passes E47's own scores and fails a seed moved by more than 0.01; the reading
applies the row's clauses P1 to P6 on the merged five-truth series, the void conditions, and the outcome-to-file rule
(at most two candidates, the ordering and ties, F1's block condition); the read refuses to overwrite and to run
without the row."""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from vec import evaluate as ev
from vec import harness
from vec import shape_interp as si
from vec.boards import Board

REPO = Path(__file__).resolve().parents[1]
SEEDS = list(range(10))
COMMIT = "a" * 40
ARMS = ("fgw_shape0x_sr", "fgw_shape0x_sg", "fgw_shape0x_ov")
CTLS = ("fgw_shape0x_sr_anti", "fgw_shape0x_sg_anti", "fgw_shape0x_sg_rot", "fgw_shape0x_ov_rot")
METRICS = ("de_score", "de_direction", "mmd_u", "variogram", "d2_shape", "occupancy_dice", "scale_log_ratio",
           "neighborhood_mmd")


@pytest.fixture(scope="module")
def e49():
    spec = importlib.util.spec_from_file_location("e49_shape_interp_read",
                                                  REPO / "experiments" / "coords" / "e49_shape_interp_read.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


# ------------------------------------------------------------------------------------------------ constants, row
def test_the_constants_name_the_row_its_arms_controls_and_gate(e49):
    assert e49.ARMS == {"fgw_shape0x_sr": ("fgw_shape0x_sr_anti",),
                        "fgw_shape0x_sg": ("fgw_shape0x_sg_anti", "fgw_shape0x_sg_rot"),
                        "fgw_shape0x_ov": ("fgw_shape0x_ov_rot",)}
    assert e49.TRY_IDS["fgw_shape0x_ov_rot"] == "E49g" and e49.C0 == "fgw_shape0x"
    assert e49.BUILDS == ("main", "prefix") and e49.TEST_LIKE == ("half", "late", "new1", "new2")
    assert (e49.GATE_TOL, e49.E47_RUN_COMMIT, e49.E47_FILES_COMMIT) == (0.01, "9967cf1", "d5ee9b5")
    assert (e49.BAR, e49.AHEAD, e49.EXPR_TOL, e49.SCALE_TOL) == (0.3, 8, 1e-9, 1e-6)
    assert e49.TIE_ORDER == ("fgw_shape0x_sg", "fgw_shape0x_sr", "fgw_shape0x_ov") and e49.MAX_CANDIDATES == 2
    assert e49.DESCRIBE_SEEDS == tuple(range(1_000_001, 1_000_011))
    assert all(e49.FORM_OF[a] in si.FORMS for a in e49.ARMS)


def test_a_run_is_refused_without_row_e49(e49, tmp_path, monkeypatch):
    tries = tmp_path / "TRIES.md"
    tries.write_text("| E48 (one try) | embryo |\n")
    monkeypatch.setattr(e49, "TRIES", tries)
    with pytest.raises(SystemExit, match="E49 has no row"):
        e49._refuse_unless_ready(COMMIT)
    for first in ("E49a to E49g (one row)", "E49 (one row)"):
        tries.write_text(f"| {first} | embryo |\n")
        assert e49.has_row(tries)


# ------------------------------------------------------------------------------------------------ step 0
def test_the_drop_rule(e49):
    tl = {"half": "1|2,3", "late": "1|2|3", "new1": "1|2,3", "new2": "1|2,3"}
    assert not e49.drop_decision("1|2,3", tl)["drop_f1"]
    d = e49.drop_decision("1,2,3", tl)
    assert d["drop_f1"] and d["f1_identity"]
    d = e49.drop_decision("1,2|3", tl)
    assert d["drop_f1"] and d["f1_unmatched"] and not d["f1_identity"]
    assert e49.dropped({"drop": d}) == ("fgw_shape0x_sr", "fgw_shape0x_sr_anti")
    assert e49.run_names(e49.dropped({"drop": d})) == ("fgw_shape0x_sg", "fgw_shape0x_sg_anti",
                                                         "fgw_shape0x_sg_rot", "fgw_shape0x_ov",
                                                         "fgw_shape0x_ov_rot")


def test_blocks_of_a_pair_ignore_place_and_size(e49):
    rng = np.random.default_rng(0)
    A = rng.normal(size=(8000, 3)) * [3.0, 1.0, 1.0]
    B = rng.normal(size=(8000, 3)) * [3.0, 2.0, 1.0]
    assert e49.blocks_of_pair(A, B) == "1|2,3" == e49.blocks_of_pair(A * 7 + 5, B)


# ------------------------------------------------------------------------------------------------ the wrapped arms
def _ctx(s=3):
    return SimpleNamespace(seed=harness.BUILD_SEED_OFFSET + s)


def _stage(seed):
    rng = np.random.default_rng(seed)
    return ev.Stage(X=rng.gamma(2.0, 1.0, size=(50, 4)).astype(np.float32), coords=rng.normal(size=(50, 3)))


def test_a_refused_form_is_recorded_and_c0s_file_stands_in(e49):
    board = Board(key="pseudo:T2e/x", task="T2", genes=["a"], needs_coords=True, split="pseudo", min_cells=1,
                  max_cells=50, truth_cells=10)
    c0 = harness.Candidate("fgw_shape0x", frozenset({"T2"}), (), lambda ctx: _stage(1), {})

    def refuse(ctx):
        raise si.ShapeRefused("stretch_registered: out of range")

    def guard(ctx):
        raise ValueError("pseudo: cloud guard flags ['coverage']")

    def other(ctx):
        raise ValueError("something else")

    with e49.captured() as seen:
        seen["board"] = board
        ctl = e49.observed(c0, seen).build(_ctx())
        for name, fn in (("fgw_shape0x_sg", refuse), ("fgw_shape0x_ov", guard)):
            out = e49.guarded(harness.Candidate(name, frozenset({"T2"}), (), fn, {}), seen).build(_ctx())
            np.testing.assert_array_equal(out.coords, ctl.coords)
        with pytest.raises(ValueError, match="something else"):
            e49.guarded(harness.Candidate("fgw_shape0x_sr", frozenset({"T2"}), (), other, {}), seen).build(_ctx())
        ok = e49.guarded(harness.Candidate("fgw_shape0x_sr", frozenset({"T2"}), (), lambda ctx: _stage(2), {}),
                         seen).build(_ctx())
    assert seen["refused"][("fgw_shape0x_sg", 3)].startswith("shape:")
    assert seen["refused"][("fgw_shape0x_ov", 3)] == "cloud guard"
    assert seen["refused"][("fgw_shape0x_sr", 3)] is None and ok.n == 50
    assert seen["x"][("fgw_shape0x_sg", 3)] == seen["x"][("fgw_shape0x", 3)]
    test = Board(**{**board.__dict__, "key": "T2:embryo:test_interp", "split": "test_interp"})
    with e49.captured() as seen:
        seen["board"] = test
        with pytest.raises(RuntimeError, match="test board"):
            e49.guarded(c0, seen).build(_ctx())


def test_the_spy_records_each_build_by_form_and_seed(e49, monkeypatch):
    calls = []

    def fake(*a, **kw):
        calls.append(kw.get("shape_form"))
        return SimpleNamespace(record={"n_jittered": 0, "floor_exact": {"fallback": False},
                                       "shape_interp": kw.get("shape_form"), "shape_interp_blocks": "1|2,3"})

    monkeypatch.setattr(e49.tb, "interpolate_clouds", fake)
    with e49.captured() as seen:
        e49.tb.interpolate_clouds(1, 2, 3, 4, 0.5, 100, 1_000_004, shape_form="stretch_ratio")
        e49.tb.interpolate_clouds(1, 2, 3, 4, 0.5, 100, 1_000_004)
    assert e49.tb.interpolate_clouds is fake  # the spy is removed afterwards
    assert seen["records"][("stretch_ratio", 1_000_004)]["shape_interp_blocks"] == "1|2,3"
    assert seen["records"][(None, 1_000_004)]["n_jittered"] == 0


def test_the_real_candidates_build_through_the_wrappers_on_a_synthetic_pair(e49):
    """The seven harness candidates on test_e47_floor_exact's synthetic embryo pair: each keeps C0's expression rows
    and count, moves positions, and leaves one interpolate_clouds record per form at the build seed."""
    from test_e47_floor_exact import _ctx as e47_ctx

    ctx = e47_ctx()
    s = ctx.seed - harness.BUILD_SEED_OFFSET
    with e49.captured() as seen:
        seen["board"] = ctx.board
        c0 = e49.observed(harness.CANDIDATES["fgw_shape0x"], seen).build(ctx)
        for n in e49.run_names():
            out = e49.guarded(harness.CANDIDATES[n], seen).build(ctx)
            assert seen["x"][(n, s)] == seen["x"][("fgw_shape0x", s)] and out.n == c0.n
            if seen["refused"][(n, s)] is None:
                assert not np.allclose(out.coords, c0.coords)
                assert seen["records"][(e49.FORM_OF[n], ctx.seed)]["shape_interp"] == e49.FORM_OF[n]
    assert seen["records"][(None, ctx.seed)]["fallback"] in (True, False)
    assert seen["records"][("stretch_ratio", ctx.seed)]["shape_interp_blocks"]


# ------------------------------------------------------------------------------------------------ the gate
def test_the_gate_passes_e47s_own_scores_and_fails_a_seed_moved_by_more_than_001(e49):
    key = e49.nt.all_tasks()["e75"].key
    old, why = e49.e47_scores(key)
    assert why == [] and set(e49.GATE_ARMS) <= set(old)
    g = e49.gate_check("e75", key, {c: old[c].copy() for c in e49.GATE_ARMS})
    assert g["pass"] and g["worst"] == {"mean": 0.0, "seed": 0.0}
    moved = {c: old[c].copy() for c in e49.GATE_ARMS}
    moved["fgw_shape0x"][4] += 0.011
    g = e49.gate_check("e75", key, moved)
    assert not g["pass"] and any("in a seed" in w for w in g["why"])
    assert not e49.gate_check("e75", "T2e/none", moved)["pass"]


def test_gate_files_are_required_at_the_run_commit(e49, tmp_path):
    assert len(e49.gate_problems(tmp_path, COMMIT)) == 9
    for name in e49.nt.all_tasks():
        p = e49._gate_path(tmp_path, name)
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(json.dumps({"pass": True, "commit": COMMIT}))
    assert e49.gate_problems(tmp_path, COMMIT) == []
    assert len(e49.gate_problems(tmp_path, "c" * 40)) == 9


# ------------------------------------------------------------------------------------------------ the reading
NAMES = ("floor", "fgw_shape0x", *ARMS, *CTLS)


def _row(seed: int, shift: dict, s_shift: dict, nb: dict) -> dict:
    """A seed's row: C0 at 47 board points; a file's board, gated shape group, shape terms and neighbourhood moved by
    the given shifts (S = the mean of d2_shape and occupancy_dice skills)."""
    scores, skills, groups = {}, {}, {}
    for n in NAMES:
        b = 50.0 if n == "floor" else 47.0 + 0.01 * seed + shift.get(n, 0.0)
        sk = {m: 0.5 for m in METRICS}
        sk["d2_shape"] = sk["occupancy_dice"] = 0.5 + s_shift.get(n, 0.0)
        sk["neighborhood_mmd"] = 0.4 + nb.get(n, 0.0)
        scores[n], skills[n] = b, sk
        groups[n] = {"shape_growth": 0.4 + s_shift.get(n, 0.0) * (2 / 3), "expression_change": 0.5,
                     "cell_state": 0.5, "local_organisation": sk["neighborhood_mmd"]}
    return {"seed": seed, "scores": scores, "skills": skills, "groups": groups}


def _facts(blocks="1|2,3", refused=None):
    names = ("fgw_shape0x", *ARMS, *CTLS)
    recs = {n: {"n_jittered": 0, "fallback": False} for n in names}
    recs["fgw_shape0x_sr"]["shape_interp_blocks"] = blocks
    return {"board_split": "pseudo", "x": dict.fromkeys(names, "x"), "n": dict.fromkeys(names, 100),
            "refused": {n: (refused or {}).get(n) for n in names[1:]}, "records": recs,
            "hands": {n: {"prefix_mode": "proper", "dice_proper": 0.8, "dice_mirror": 0.7, "dice_gap": 0.1}
                      for n in names},
            "files_identical_across_builds": True, "refusals_identical_across_builds": True}


def _recs(e49, shift=None, s_shift=None, nb=None, refused_on=None, seed_shift=None, blocks=None):
    out = {}
    for name in e49.nt.all_tasks():
        rs = []
        for s in SEEDS:
            sh = dict(shift or {})
            if seed_shift and s in seed_shift:
                sh.update(seed_shift[s])
            row = _row(s, sh, s_shift or {}, nb or {})
            builds = {b: {"provenance": {"commit": COMMIT, "tree_clean": True}, "build": b, "protocol": "bake",
                          "board": {"split": "pseudo"}, "rows": [json.loads(json.dumps(row))]} for b in e49.BUILDS}
            ref = (refused_on or {}).get((name, s))
            rs.append({"seed": s, "run_commit": COMMIT, "dropped": [], "builds": builds,
                       "facts": _facts((blocks or {}).get(name, "1|2,3"), ref)})
        out[name] = rs
    return out


DESCRIBE = {"drop": {"test_pair_blocks": "1|2,3", "test_like_blocks": {}, "f1_identity": False,
                     "f1_unmatched": False, "drop_f1": False}}
GOOD = {"fgw_shape0x_sg": 1.0, "fgw_shape0x_sr": 0.6, "fgw_shape0x_ov": 0.8}
GOOD_S = {"fgw_shape0x_sg": 0.03, "fgw_shape0x_sr": 0.02, "fgw_shape0x_ov": 0.02}


def test_three_passing_forms_give_two_candidates_by_merged_board(e49):
    res = e49.reading(_recs(e49, GOOD, GOOD_S), [], COMMIT, DESCRIBE)
    assert not res["void"]
    assert all(res["arms"][a]["passes"] for a in ARMS)
    cands = res["mapping"]["candidates"]
    assert [c["arm"] for c in cands] == ["fgw_shape0x_sg", "fgw_shape0x_ov"]
    assert "--shape-interp stretch_registered --shape-interp-test-ok" in cands[0]["command"]
    assert "--positions floor_exact" in cands[0]["command"]
    m = res["arms"]["fgw_shape0x_sg"]["merged"]["board"]
    assert m["mean"] == pytest.approx(1.0) and m["ahead"] == 10
    text = e49.text(res)
    assert "PASSES" in text and "VOID" not in text


def test_ties_go_to_the_test_like_board_then_the_fixed_order(e49):
    res = e49.reading(_recs(e49, dict.fromkeys(ARMS, 0.8), GOOD_S), [], COMMIT, DESCRIBE)
    assert [c["arm"] for c in res["mapping"]["candidates"]] == ["fgw_shape0x_sg", "fgw_shape0x_sr"]


@pytest.mark.parametrize("clause", ["P1_bar", "P1_seeds", "P2", "P4", "P5", "P6"])
def test_each_clause_can_stop_a_form(e49, clause):
    shift, s_shift, nb, refused, seed_shift = dict(GOOD), dict(GOOD_S), {}, None, None
    a = "fgw_shape0x_sg"
    if clause == "P1_bar":
        shift[a] = 0.25
    elif clause == "P1_seeds":
        seed_shift = {s: {a: -1.0} for s in (0, 1, 2)}  # mean still above 0.3, only 7 seeds ahead
    elif clause == "P2":
        s_shift[a] = -0.02
    elif clause == "P4":
        nb[a] = -0.01
    elif clause == "P5":
        s_shift.update({"fgw_shape0x_sg_rot": 0.05})
    else:
        refused = {("x2", 5): {"fgw_shape0x_sg_anti": "cloud guard"}}
    res = e49.reading(_recs(e49, shift, s_shift, nb, refused, seed_shift), [], COMMIT, DESCRIBE)
    assert not res["void"]
    r = res["arms"][a]
    assert not r["passes"]
    key = clause.split("_")[0]
    assert not r["clauses"][key]
    assert a not in [c["arm"] for c in res["mapping"]["candidates"]]


def test_nothing_passing_names_no_file(e49):
    res = e49.reading(_recs(e49), [], COMMIT, DESCRIBE)
    assert res["mapping"]["candidates"] == [] and "no candidate file" in res["mapping"]["notes"][0]


def test_f1_needs_its_test_pair_blocks_on_a_test_like_task(e49):
    shift = {"fgw_shape0x_sr": 1.0}
    s_shift = {"fgw_shape0x_sr": 0.02}
    blocks = dict.fromkeys(("half", "late", "new1", "new2"), "1|2|3")
    res = e49.reading(_recs(e49, shift, s_shift, blocks=blocks), [], COMMIT, DESCRIBE)
    assert res["arms"]["fgw_shape0x_sr"]["passes"]
    assert res["mapping"]["candidates"] == [] and "unpriced at the test" in res["mapping"]["notes"][0]
    res = e49.reading(_recs(e49, shift, s_shift), [], COMMIT, DESCRIBE)
    assert [c["arm"] for c in res["mapping"]["candidates"]] == ["fgw_shape0x_sr"]


def test_a_dropped_f1_is_neither_required_nor_read(e49):
    recs = _recs(e49, GOOD, GOOD_S)
    for rs in recs.values():
        for r in rs:
            r["dropped"] = ["fgw_shape0x_sr", "fgw_shape0x_sr_anti"]
            for b in r["builds"].values():
                for k in ("scores", "skills", "groups"):
                    for n in ("fgw_shape0x_sr", "fgw_shape0x_sr_anti"):
                        b["rows"][0][k].pop(n)
    d = {"drop": {**DESCRIBE["drop"], "drop_f1": True, "f1_identity": True}}
    res = e49.reading(recs, [], COMMIT, d)
    assert not res["void"] and res["arms"]["fgw_shape0x_sr"] == {"dropped": True}
    assert "dropped at step 0" in e49.text(res)


@pytest.mark.parametrize("breach", ["seed", "commit", "unclean", "expression", "count", "expr_skill", "scale",
                                    "builds", "refusal_builds", "gate", "test_board", "missing_arm", "describe"])
def test_the_void_conditions(e49, breach):
    recs = _recs(e49, GOOD, GOOD_S)
    why, describe = [], DESCRIBE
    r = recs["late"][3]
    f = r["facts"]
    if breach == "seed":
        recs["late"].pop(3)
    elif breach == "commit":
        r["builds"]["prefix"]["provenance"]["commit"] = "b" * 40
    elif breach == "unclean":
        r["builds"]["main"]["provenance"]["tree_clean"] = False
    elif breach == "expression":
        f["x"]["fgw_shape0x_ov"] = "y"
    elif breach == "count":
        f["n"]["fgw_shape0x_ov"] = 99
    elif breach == "expr_skill":
        r["builds"]["main"]["rows"][0]["skills"]["fgw_shape0x_sg"]["mmd_u"] += 1e-8
    elif breach == "scale":
        r["builds"]["main"]["rows"][0]["skills"]["fgw_shape0x_sg"]["scale_log_ratio"] += 1e-5
    elif breach == "builds":
        f["files_identical_across_builds"] = False
    elif breach == "refusal_builds":
        f["refusals_identical_across_builds"] = False
    elif breach == "gate":
        why = ["the gate did not pass on late at the run commit"]
    elif breach == "test_board":
        f["board_split"] = "test_interp"
    elif breach == "missing_arm":
        r["builds"]["main"]["rows"][0]["scores"].pop("fgw_shape0x_ov_rot")
    else:
        describe = None
    res = e49.reading(recs, why, COMMIT, describe)
    assert res["void"] and res["why"]
    assert e49.text(res).startswith(f"Row E49 (run commit {COMMIT}): VOID")


def test_a_jittered_seed_is_exempt_from_the_scale_check(e49):
    recs = _recs(e49, GOOD, GOOD_S)
    r = recs["late"][3]
    r["builds"]["main"]["rows"][0]["skills"]["fgw_shape0x_sg"]["scale_log_ratio"] += 1e-5
    r["facts"]["records"]["fgw_shape0x_sg"]["n_jittered"] = 2
    assert not e49.reading(recs, [], COMMIT, DESCRIBE)["void"]


def test_the_read_refuses_to_overwrite(e49, tmp_path):
    (tmp_path / "reading.txt").write_text("x")
    with pytest.raises(SystemExit, match="refusing to overwrite"):
        e49.main(["read", "--commit", COMMIT, "--runs", str(tmp_path), "--out", str(tmp_path)])


def test_the_merge_is_by_truth_then_seed(e49):
    """D_s is the mean over the five truths of each truth's task mean at seed s (E8.5 has three tasks)."""
    recs = _recs(e49)
    for r in recs["x2"]:
        for b in r["builds"].values():
            b["rows"][0]["scores"]["fgw_shape0x_ov"] += 3.0
    groups = e49.e47.truths(e49.nt.all_tasks())
    D = e49.merged_D(recs, groups, "fgw_shape0x_ov", "fgw_shape0x", "board")
    np.testing.assert_allclose(D, np.full(10, 3.0 / 3 / 5))


# ------------------------------------------------------------------------------------------------ the launcher
def test_the_launcher_caps_threads_runs_nine_processes_and_writes_exit_files():
    sh = (REPO / "experiments" / "coords" / "run_e49.sh").read_text()
    for need in ("OMP_NUM_THREADS=8", "OPENBLAS_NUM_THREADS=8", "MKL_NUM_THREADS=8", "NUMPY_MADVISE_HUGEPAGE=0",
                 "finish 98", "finish 99", '"$S/$STEP.exit"', "uv run --frozen --no-sync python -I -B",
                 "tests/test_e49_shape_interp_read.py", "tests/test_shape_interp.py", "worktree", "describe)"):
        assert need in sh
    tasks = sh.split("TASKS=(")[1].split(")")[0].split()
    assert sorted(tasks) == sorted(["half", "late", "new1", "new2", "e75", "e725", "x1", "x2", "x3"])
    assert len(tasks) <= 10 and "gate)" in sh and "run)" in sh and "read)" in sh
