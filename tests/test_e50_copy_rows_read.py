"""Row E50 as specified: the read script experiments/coords/e50_copy_rows_read.py (every branch of its fixed reading,
the gate, the void conditions) on synthetic per-seed records, and its launcher's refusal paths. Nothing here reads
data, runs a task or opens a result."""

from __future__ import annotations

import dataclasses
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from vec import evaluate as ev

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "experiments" / "coords"))
import e50_copy_rows_read as e50  # noqa: E402

C = "c" * 40
TASKS = list(e50.nt.all_tasks())
NAMES = (e50.harness.FLOOR, *e50.CONTROLS, e50.ARM)
POS = {"d2_shape": 0.1, "occupancy_dice": 0.6, "scale_log_ratio": 0.02, "sliced_wasserstein": 0.3}


def _skills(nbhd=0.8, shift=0.0):
    return {m: (nbhd if m == e50.NBHD else 0.4 + shift) for m in e50.METRICS}


def _row(seed, arm_board, arm_nbhd=0.8, unequal=None):
    skills = {n: _skills() for n in NAMES}
    skills[e50.ARM] = _skills(arm_nbhd, 0.05)
    groups = {n: ev.aggregate("T2", s)[1] for n, s in skills.items()}
    metrics = {n: {**POS, "mmd_u": 0.01 * i} for i, n in enumerate(NAMES)}
    if unequal:
        metrics[e50.ARM][unequal] = POS[unequal] + 1e-12
    scores = {e50.harness.FLOOR: 50.0, e50.PRIMARY: 51.0, e50.SECONDARY: 51.5, e50.ARM: arm_board}
    return {"seed": seed, "scores": scores, "skills": skills, "groups": groups, "metrics": metrics}


def _rec(task, seed, arm_board=52.0, refused=None, unequal=None, arm_nbhd=0.8, commit=C, coords=True):
    builds = {}
    for b in e50.BUILDS:
        term = unequal[1] if unequal and unequal[0] == b else None
        builds[b] = {"provenance": {"commit": commit, "tree_clean": True}, "build": b, "protocol": "bake",
                     "board": {"split": "pseudo"}, "rows": [_row(seed, arm_board, arm_nbhd, term)]}
    return {"task": task, "seed": seed, "run_commit": commit, "arm": e50.ARM, "builds": builds,
            "facts": {"board_split": "pseudo", "refused": refused, "coords_equal_floor": refused is None and coords,
                      "prefix_mode": {"arm": None if refused else "proper", "floor": "proper"},
                      "files_identical_across_builds": True, "refusals_identical_across_builds": True}}


def _recs(**over):
    return {t: [_rec(t, s, **over) for s in e50.SEEDS] for t in TASKS}


# ------------------------------------------------------------------------------------------------ reading (i)
def test_i_holds_and_the_file_is_built():
    res = e50.reading(_recs(), [], C)
    assert not res["void"] and res["i"]["holds"] and res["build_test_file"]
    assert res["i"]["built_task_seeds"] == len(TASKS) * len(e50.SEEDS) and res["outcome"].startswith("build")
    assert "HOLDS" in e50.text(res)


@pytest.mark.parametrize("build", ["main", "prefix"])
@pytest.mark.parametrize("term", ["d2_shape", "occupancy_dice", "scale_log_ratio"])
def test_i_fails_on_any_inequality_in_either_build(build, term):
    recs = _recs()
    recs["x2"][4] = _rec("x2", 4, unequal=(build, term))
    res = e50.reading(recs, [], C)
    assert not res["void"] and not res["i"]["holds"] and not res["build_test_file"]
    assert res["i"]["unequal"] == [{"task": "x2", "seed": 4, "build": build, "what": term}]
    assert res["outcome"].startswith("no test file") and "FAILS" in e50.text(res)


def test_i_fails_when_the_coordinates_are_not_the_floor_copy():
    recs = _recs()
    recs["half"][0] = _rec("half", 0, coords=False)
    assert not e50.reading(recs, [], C)["i"]["holds"]


def test_a_reported_term_and_a_mode_difference_decide_nothing():
    recs = _recs()
    recs["late"][2] = _rec("late", 2, unequal=("main", "sliced_wasserstein"))
    recs["late"][3]["facts"]["prefix_mode"]["arm"] = "mirror"
    res = e50.reading(recs, [], C)
    assert res["i"]["holds"] and len(res["i"]["reported_unequal"]) == 1 and len(res["i"]["prefix_mode_differs"]) == 1


# ------------------------------------------------------------------------------------------------ reading (ii)
def test_ii_lists_refusals_and_a_task_refused_everywhere_prices_nothing():
    recs = _recs()
    recs["e75"] = [_rec("e75", s, refused="count_outside_limits") for s in e50.SEEDS]
    recs["x1"][7] = _rec("x1", 7, refused="dropped_over_share", unequal=("main", "d2_shape"))  # not read: refused
    res = e50.reading(recs, [], C)
    assert res["i"]["holds"]  # a refused task-seed is left out of (i)
    assert res["ii"]["refused_at_every_seed"] == ["e75"]
    assert res["ii"]["refused"]["x1"] == {"7": "dropped_over_share"} and len(res["ii"]["refused"]["e75"]) == 10
    assert res["iv"]["e75"] is None
    board = res["iii"]["primary"]["main"]["board"]
    assert "e75" not in board["per_task"] and board["per_task"]["x1"]["seeds"] == 9
    assert "refused at every seed" in e50.text(res)


def test_refusal_categories():
    assert e50.refusal_category("copy_rows: more copy slots lack a deduplicated cell than the later bracket's share "
                                "(values not shown)") == "dropped_over_share"
    assert e50.refusal_category("pseudo: 9 cells outside [10, 20]") == "count_outside_limits"
    assert e50.refusal_category("syn: cloud guard flags ['rms']") is None
    assert e50.refusal_category("copy_rows: non-finite input (values not shown)") is None
    assert e50.refusal_category("copy_rows: fewer than three rows in the copy") == "copy_too_small"
    assert e50.refusal_category("copy_rows: the released rows and the duplicate filter's record do not line up "
                                "(values not shown)") is None
    assert e50.refusal_category("something else") is None
    # the analyst's decision: the row's three refusals and no other (an arm error outside them voids the run)
    assert {cat for _, cat in e50.REFUSALS} == {"dropped_over_share", "count_outside_limits", "copy_too_small"}


def test_a_refusal_outside_the_list_raises():
    """An arm error outside the listed refusals is raised (the run is void as a missing arm), never recorded."""
    seen = {"stage": {}, "file": {}, "refused": {}, "board": type("B", (), {"split": "val", "key": "T2:embryo:val"})()}

    def bad(ctx):
        raise ValueError("copy_rows: non-finite input (values not shown)")

    arm = e50.guarded(dataclasses.replace(e50.harness.CANDIDATES[e50.ARM], build=bad), seen)
    ctx = type("C", (), {"seed": e50.harness.BUILD_SEED_OFFSET + 3})()
    with pytest.raises(ValueError, match="non-finite"):
        arm.build(ctx)
    assert seen["refused"] == {}


# ------------------------------------------------------------------------------------------------ reading (iii), (iv)
def test_iii_by_hand():
    gain = {"half": 2.0, "late": 2.0, "new1": -1.0, "new2": 0.5, "e75": -3.0, "e725": 1.0, "x1": 0.0, "x2": 0.5,
            "x3": -0.5}
    recs = {t: [_rec(t, s, arm_board=51.0 + gain[t]) for s in e50.SEEDS] for t in TASKS}
    res = e50.reading(recs, [], C)
    for tier, base in (("primary", 51.0), ("secondary", 51.5)):
        b = res["iii"][tier]["main"]["board"]
        want = {t: g + 51.0 - base for t, g in gain.items()}
        assert {t: v["mean"] for t, v in b["per_task"].items()} == pytest.approx(want)
        merged = np.mean([np.mean([want[n] for n in names]) for names in res["truths"].values()])
        assert b["merged"] == pytest.approx(merged)
        assert b["D_s"] == pytest.approx([merged] * 10) and b["se"] == pytest.approx(0.0, abs=1e-12)
        pos, neg = [v for v in want.values() if v > 0], [v for v in want.values() if v < 0]
        assert b["gain_tasks"] == len(pos) and b["mean_gain"] == pytest.approx(np.mean(pos))
        assert b["loss_tasks"] == len(neg) and b["mean_loss"] == pytest.approx(np.mean(neg))
        assert b["worst_task"] == "e75" and b["best_task"] in ("half", "late")
        assert res["iii"][tier]["against"] == (e50.PRIMARY if tier == "primary" else e50.SECONDARY)
    p = res["iii"]["primary"]
    assert set(p["prefix"]) == set(p["main"]) == {"board", *(f"group:{g}" for g in e50.GROUP_NAMES),
                                                  *(f"skill:{m}" for m in e50.METRICS)}
    assert p["main"]["skill:mmd_u"]["merged"] == pytest.approx(0.05)
    assert p["main"]["board"]["seeds_ahead"] == (10 if p["main"]["board"]["merged"] > 0 else 0)


def test_iv_gate_bites_and_shape_gated_and_ungated():
    recs = _recs()
    for s in (0, 1, 2):
        recs["new2"][s] = _rec("new2", s, arm_nbhd=0.25)
    res = e50.reading(recs, [], C)
    iv = res["iv"]["new2"]["main"]
    assert iv["seeds_gate_bites"] == 3 and res["iv"]["half"]["main"]["seeds_gate_bites"] == 0
    ungated = ev.aggregate("T2", _skills(0.25, 0.05), gate="none")[1]["shape_growth"]
    gated = ev.aggregate("T2", _skills(0.25, 0.05))[1]["shape_growth"]
    assert gated < ungated
    assert iv["shape_ungated"] == pytest.approx(ungated)
    assert iv["shape_gated"] == pytest.approx((3 * gated + 7 * ungated) / 10)


# ------------------------------------------------------------------------------------------------ void
@pytest.mark.parametrize("breaker, needle", [
    (lambda r: r["x3"].pop(), "seeds are not"),
    (lambda r: r["x3"][0].update(run_commit="d" * 40), "another run commit"),
    (lambda r: r["x3"][0].update(arm="fgw_copyrows_other"), "another arm"),
    (lambda r: r["x3"][0]["builds"].pop("prefix"), "build prefix missing"),
    (lambda r: r["x3"][0]["builds"]["main"]["provenance"].update(tree_clean=False), "unclean"),
    (lambda r: r["x3"][0]["builds"]["main"].update(protocol="emulator"), "bake protocol"),
    (lambda r: r["x3"][0]["builds"]["main"]["board"].update(split="test"), "a test board"),
    (lambda r: r["x3"][0]["facts"].update(board_split="test_interp"), "a test board"),
    (lambda r: r["x3"][0]["builds"]["main"]["rows"][0]["scores"].pop(e50.ARM), "lacks an arm"),
    (lambda r: r["x3"][0]["facts"].update(files_identical_across_builds=False), "a file differs"),
    (lambda r: r["x3"][0]["facts"].update(refusals_identical_across_builds=False), "a refusal differs"),
])
def test_void_conditions(breaker, needle):
    recs = _recs()
    breaker(recs)
    res = e50.reading(recs, [], C)
    assert res["void"] and any(needle in w for w in res["why"]) and "VOID" in e50.text(res)


def test_a_gate_or_commit_problem_voids():
    res = e50.reading(_recs(), ["the gate did not pass on x1 at the run commit"], C)
    assert res["void"] and res["why"] == ["the gate did not pass on x1 at the run commit"]


# ------------------------------------------------------------------------------------------------ gate
def _gate(floor=50.0, d0=0.0, d1=0.0, dseed=0.0):
    old = {a: np.full(10, 55.0) for a in e50.GATE_ARMS}
    new = {e50.harness.FLOOR: np.full(10, floor), e50.PRIMARY: np.full(10, 55.0 + d0),
           e50.SECONDARY: np.full(10, 55.0 + d1)}
    new[e50.SECONDARY][3] += dseed
    return e50.gate_check("x1", "key", new, old, [])


def test_gate_branches():
    assert _gate()["pass"]
    assert _gate(d0=0.09, d1=-0.09)["pass"]
    assert not _gate(floor=50.1)["pass"]
    assert not _gate(d0=0.11)["pass"]
    assert not _gate(dseed=0.6)["pass"]  # within 0.1 in the mean, off by 0.6 in one seed
    assert _gate(dseed=0.45)["pass"]
    assert not e50.gate_check("x1", "key", {e50.harness.FLOOR: np.full(10, 50.0)}, {}, [])["pass"]
    assert not e50.gate_check("x1", "key", {}, None, ["no committed E47 file"])["pass"]


def test_has_row(tmp_path):
    p = tmp_path / "TRIES.md"
    p.write_text("| E49a to E49g (x) | embryo |\n")
    assert not e50.has_row(p)
    p.write_text("| E50 (one try) | embryo |\n")
    assert e50.has_row(p)


def test_same_bits():
    a = np.random.default_rng(0).normal(size=(5, 3)).astype(np.float32)
    b = a.copy()
    assert e50.same_bits(a, b) and e50.same_bits(a.astype(np.float64), b)
    b[2, 2] = np.nextafter(b[2, 2], np.float32(1e9))
    assert not e50.same_bits(a, b) and not e50.same_bits(a, a[::-1])


def test_only_the_arm_of_the_handed_over_form_runs():
    # the handed-over T2embryo_s2x.h5ad is fgw_shape0x built with --expression mixture (embryo2_2dccfee.txt)
    assert e50.HANDED_OVER_EXPRESSION == "mixture" and e50.ARM == "fgw_copyrows"
    assert e50.ARM == e50.ARMS_BY_EXPRESSION[e50.HANDED_OVER_EXPRESSION]
    assert e50.PRIMARY == "fgw_shape0x" and "fgw_copyrows_series" not in [
        c.name for c in e50.candidates({"stage": {}, "file": {}, "refused": {}})]
    assert set(e50.ARMS_BY_EXPRESSION.values()) == {"fgw_copyrows", "fgw_copyrows_series"}
    assert [c.name for c in e50.candidates({"stage": {}, "file": {}, "refused": {}})] == [
        "floor", "fgw_shape0x", "fgw_shape0x_series", e50.ARM]


# ------------------------------------------------------------------------------------------------ the launcher
SH = REPO / "experiments" / "coords" / "run_e50.sh"


def _launch(tmp_path, *args):
    env = {k: v for k, v in os.environ.items() if not k.startswith("GIT_")}
    env.update(HOME=str(tmp_path), GIT_CEILING_DIRECTORIES=str(tmp_path.parent))
    return subprocess.run(["bash", str(SH), *args], cwd=tmp_path, env=env, capture_output=True, text=True,
                          timeout=60)


def test_the_launcher_refuses_out_of_order_and_writes_an_exit_file(tmp_path):
    S = tmp_path / "v" / "scratch" / "e50"
    r = _launch(tmp_path)
    assert r.returncode == 2 and "usage" in r.stderr and (S / "usage.exit").read_text().strip() == "2"
    r = _launch(tmp_path, "gate", C)
    assert r.returncode == 96 and "step tests has not exited 0" in r.stderr
    assert (S / "gate.exit").read_text().strip() == "96"
    (S / "gate.exit").unlink()
    (S / "gate.exit").write_text("0\n")
    (S / "gate.head").write_text("e" * 40 + "\n")  # passed, but at another commit
    r = _launch(tmp_path, "run", C)
    assert r.returncode == 96 and (S / "run.exit").read_text().strip() == "96"
    r = _launch(tmp_path, "bogus", C)
    assert r.returncode == 97 and (S / "bogus.exit").read_text().strip() == "97"
    env = {k: v for k, v in os.environ.items() if not k.startswith("GIT_")}
    env.update(HOME=str(tmp_path), GIT_CEILING_DIRECTORIES=str(tmp_path.parent))
    trace = subprocess.run(["bash", "-x", str(SH), "read", C], cwd=tmp_path, env=env, capture_output=True, text=True,
                           timeout=60)  # the dry-run trace of a harmless path: it stops before git or Python
    assert trace.returncode == 96 and "finish 96" in trace.stderr and "uv run" not in trace.stderr


def test_the_launcher_has_no_unsafe_parameter_message_and_parses():
    text = SH.read_text()
    assert ":?" not in text
    assert subprocess.run(["bash", "-n", str(SH)], capture_output=True).returncode == 0
    for f in ("tests/test_copy_rows.py", "tests/test_copy_rows_build.py", "tests/test_e50_copy_rows_read.py",
              "tests/test_harness.py", "tests/test_hidden_counts.py"):  # the review's lesser point 5
        assert f in text and (REPO / f).exists()
