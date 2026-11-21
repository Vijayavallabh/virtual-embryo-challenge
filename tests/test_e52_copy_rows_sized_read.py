"""Row E52 as specified: the read script experiments/coords/e52_copy_rows_sized_read.py (its fixed reading, the gate,
the void conditions) on synthetic per-seed records and synthetic arrays, and its launcher's refusal paths. Nothing here
reads data, runs a task or opens a result."""

from __future__ import annotations

import copy
import dataclasses
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from vec import evaluate as ev

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "experiments" / "coords"))
import e52_copy_rows_sized_read as e52  # noqa: E402

C = "c" * 40
TASKS = list(e52.nt.all_tasks())
NAMES = (e52.harness.FLOOR, *e52.CONTROLS, e52.ARM)
T = 3.25


# ------------------------------------------------------------------------------------------------ synthetic arrays
@pytest.fixture(scope="module")
def arrays():
    rng = np.random.default_rng(11)
    floor = (rng.normal(size=(400, 3)) * np.array([4.0, 2.0, 0.8]) + 30.0).astype(np.float32)
    floor[7] = floor[3]  # an exact duplicate
    truth = (rng.normal(size=(500, 3)) * np.array([3.0, 2.5, 1.0])).astype(np.float32)
    X = rng.gamma(2.0, 1.0, size=(400, 6)).astype(np.float32)
    return floor, truth, X


def _metrics(arm_expr=None):
    m = {n: {"de_score": 0.3, "de_direction": 0.6, "mmd_u": 0.02, "variogram": 0.4, "sliced_wasserstein": 0.1,
             "scale_log_ratio": 0.05, "neighborhood_mmd": 0.03, "d2_shape": 0.01, "occupancy_dice": 0.8}
         for n in NAMES}
    if arm_expr:
        m[e52.ARM][arm_expr] += 1e-9
    return m


def _facts(arrays, build, arm_coords=None, arm_X=None, T_read=T, metrics=None):
    floor, truth, X = arrays
    arm = ev.Stage(X=X if arm_X is None else arm_X,
                   coords=e52.expected_coords(floor, T) if arm_coords is None else arm_coords)
    crow = ev.Stage(X=X, coords=floor)
    return e52.raw_facts(arm, crow, ev.Stage(X=X, coords=floor), ev.Stage(X=X[:1], coords=truth), T_read, 0, build,
                         _metrics() if metrics is None else metrics, with_reported=build == "main")


@pytest.fixture(scope="module")
def good_facts(arrays):
    return {b: _facts(arrays, b) for b in e52.BUILDS}


def _skills(nbhd=0.8, shift=0.0):
    return {m: (nbhd if m == e52.NBHD else 0.4 + shift) for m in e52.METRICS}


def _row(seed, arm_board):
    skills = {n: _skills() for n in NAMES}
    skills[e52.ARM] = _skills(0.8, 0.05)
    groups = {n: ev.aggregate("T2", s)[1] for n, s in skills.items()}
    scores = {e52.harness.FLOOR: 50.0, e52.P_ARM: 51.0, e52.C_ARM: 51.2, e52.ARM: arm_board}
    return {"seed": seed, "scores": scores, "skills": skills, "groups": groups, "metrics": _metrics()}


def _rec(task, seed, facts, arm_board=51.4, refused=None, refused_c=None, commit=C):
    builds = {b: {"provenance": {"commit": commit, "tree_clean": True}, "build": b, "protocol": "bake",
                  "board": {"split": "pseudo"}, "rows": [_row(seed, arm_board)]} for b in e52.BUILDS}
    return {"task": task, "seed": seed, "run_commit": commit, "arm": e52.ARM, "builds": builds,
            "facts": {"board_split": "pseudo", "refused": refused, "refused_c": refused_c,
                      "per_build": {} if refused else copy.deepcopy(facts),
                      "files_identical_across_builds": True, "refusals_identical_across_builds": True}}


def _recs(facts, **over):
    return {t: [_rec(t, s, facts, **over) for s in e52.SEEDS] for t in TASKS}


# ------------------------------------------------------------------------------------------------ reading (i)
def test_i_holds_on_records_built_to_the_arithmetic(good_facts):
    for b in e52.BUILDS:
        assert e52.i_failures(good_facts[b]) == [], b
    res = e52.reading(_recs(good_facts), [], C)
    assert not res["void"] and res["i"]["holds"] and res["build_test_file"] and res["outcome"].startswith("build")
    assert res["i"]["built_task_seeds"] == len(TASKS) * len(e52.SEEDS) and "HOLDS" in e52.text(res)
    assert good_facts["main"]["reported"]["parity_flips_all_counts"] == 0
    assert good_facts["main"]["reported"]["nn_sets_differ_vs_c"] <= 10  # ties only


def _fails_with(good_facts, bad: dict, build: str, what: str):
    facts = {**good_facts, build: bad}
    recs = _recs(good_facts)
    recs["x2"][4] = _rec("x2", 4, facts)
    res = e52.reading(recs, [], C)
    assert not res["void"] and not res["i"]["holds"] and not res["build_test_file"]
    assert {(f["task"], f["seed"], f["build"]) for f in res["i"]["failed"]} == {("x2", 4, build)}
    assert what in {f["what"] for f in res["i"]["failed"]}
    assert res["outcome"].startswith("no test file") and "FAILS" in e52.text(res)


@pytest.mark.parametrize("build", ["main", "prefix"])
def test_i_fails_on_one_changed_coordinate_bit(arrays, good_facts, build):
    C32 = e52.expected_coords(arrays[0], T)
    C32[5, 1] = np.nextafter(C32[5, 1], np.float32(np.inf))
    bad = _facts(arrays, build, arm_coords=C32)
    assert e52.i_failures(bad) == ["coordinates"]
    _fails_with(good_facts, bad, build, "coordinates")


@pytest.mark.parametrize("build", ["main", "prefix"])
def test_i_fails_on_a_parity_flip_at_the_scored_count(arrays, good_facts, build):
    M = e52.expected_coords(arrays[0], T)
    M[:, 2] = 2 * M[:, 2].mean() - M[:, 2]  # a mirror about the centroid: same size, the other parity
    bad = _facts(arrays, build, arm_coords=M)
    assert not bad["parity_equal_scored"] and "prefix_parity" in e52.i_failures(bad)
    flipped = {**good_facts[build], "parity_equal_scored": False}
    assert e52.i_failures(flipped) == ["prefix_parity"]
    _fails_with(good_facts, flipped, build, "prefix_parity")


@pytest.mark.parametrize("build", ["main", "prefix"])
def test_i_fails_on_an_x_byte_difference(arrays, good_facts, build):
    X = arrays[2].copy()
    X[10, 2] = np.nextafter(X[10, 2], np.float32(np.inf))
    bad = _facts(arrays, build, arm_X=X)
    assert e52.i_failures(bad) == ["x_bytes"]
    _fails_with(good_facts, bad, build, "x_bytes")


@pytest.mark.parametrize("build", ["main", "prefix"])
def test_i_fails_when_the_arm_target_is_not_the_task_role_target(arrays, good_facts, build):
    """An arm sized to another target (say the reference-augmented inputs' or the copy's own) than the T this read
    computes on the task's own role stages: the coordinates and the RMS check fail."""
    other = e52.expected_coords(arrays[0], T * (1 + 1e-4))
    bad = _facts(arrays, build, arm_coords=other)
    assert {"coordinates", "rms"} <= set(e52.i_failures(bad))
    _fails_with(good_facts, bad, build, "rms")


@pytest.mark.parametrize("term", list(e52.EXPRESSION_TERMS))
def test_i_fails_on_an_expression_term_unequal_to_fgw_copyrows(arrays, good_facts, term):
    bad = _facts(arrays, "main", metrics=_metrics(arm_expr=term))
    assert e52.i_failures(bad) == [term]
    _fails_with(good_facts, bad, "main", term)


@pytest.mark.parametrize("key, val, what", [("d2_diff", 2e-6, "d2_shape"), ("frame_dev", 2e-5, "frame"),
                                            ("rms_rel", 2e-6, "rms")])
def test_i_fails_past_each_tolerance(good_facts, key, val, what):
    assert e52.i_failures({**good_facts["main"], key: val}) == [what]


def test_i_occupancy_bound(good_facts):
    f = {**good_facts["main"], "dice_arm": 0.80, "dice_floor": 0.80, "v": 0, "occupied": 200}
    assert e52.i_failures(f) == []
    assert e52.i_failures({**f, "dice_arm": 0.80 + 0.9e-4}) == []
    assert e52.i_failures({**f, "dice_arm": 0.80 + 1.1e-4}) == ["occupancy_dice"]
    assert e52.i_failures({**f, "dice_arm": 0.80 + 0.02, "v": 1}) == []  # 4 / 200 + 1e-4 covers it
    assert e52.i_failures({**f, "dice_arm": 0.80 + 0.021, "v": 1}) == ["occupancy_dice"]


def test_read_target_uses_the_task_roles_only():
    calls = []

    class Fake:
        roles = ("before", "after")

    def target_rms(ctx):
        calls.append(set(ctx.inputs))
        return 2.5

    orig = e52.harness.target_rms
    e52.harness.target_rms = target_rms
    try:
        got = e52.read_target(Fake(), None, {"before": 1, "after": 2, e52.harness.REFERENCE_RELEASED: 3,
                                             e52.harness.REFERENCE_RELEASED_ROWS: 4})
    finally:
        e52.harness.target_rms = orig
    assert got == 2.5 and calls == [{"before", "after"}]


# ------------------------------------------------------------------------------------------------ (ii) and void
def test_ii_refusals_equal_to_c_and_left_out_of_i(good_facts):
    recs = _recs(good_facts)
    recs["e75"] = [_rec("e75", s, good_facts, refused="count_outside_limits", refused_c="count_outside_limits")
                   for s in e52.SEEDS]
    res = e52.reading(recs, [], C)
    assert not res["void"] and res["i"]["holds"] and res["ii"]["refused_at_every_seed"] == ["e75"]
    assert res["ii"]["equal_to_c"] and res["iv"]["e75"] is None


@pytest.mark.parametrize("arm, c", [("dropped_over_share", None), (None, "dropped_over_share"),
                                    ("copy_too_small", "dropped_over_share")])
def test_a_refusal_that_differs_from_cs_voids(good_facts, arm, c):
    recs = _recs(good_facts)
    recs["x1"][7] = _rec("x1", 7, good_facts, refused=arm, refused_c=c)
    res = e52.reading(recs, [], C)
    assert res["void"] and any("differs from fgw_copyrows" in w for w in res["why"])


def test_a_key_error_from_the_arm_is_raised_and_voids_the_run(good_facts):
    seen = {"stage": {}, "file": {}, "refused": {}, "board": type("B", (), {"split": "val", "key": "T2:embryo:val"})()}

    def bad(ctx):
        raise KeyError("reference_released")

    arm = e52.guarded(dataclasses.replace(e52.harness.CANDIDATES[e52.ARM], build=bad), seen)
    with pytest.raises(KeyError):
        arm.build(type("Ctx", (), {"seed": e52.harness.BUILD_SEED_OFFSET + 3})())
    assert seen["refused"] == {}
    recs = _recs(good_facts)
    recs["x3"] = [r for r in recs["x3"] if r["seed"] != 3]  # the seed whose arm raised is never written
    res = e52.reading(recs, [], C)
    assert res["void"] and any("seeds are not" in w for w in res["why"])


@pytest.mark.parametrize("breaker, needle", [
    (lambda r: r["x3"][0].update(run_commit="d" * 40), "another run commit"),
    (lambda r: r["x3"][0].update(arm="fgw_copyrows"), "another arm"),
    (lambda r: r["x3"][0]["builds"].pop("prefix"), "build prefix missing"),
    (lambda r: r["x3"][0]["builds"]["main"]["provenance"].update(tree_clean=False), "unclean"),
    (lambda r: r["x3"][0]["builds"]["main"].update(protocol="emulator"), "bake protocol"),
    (lambda r: r["x3"][0]["builds"]["main"]["board"].update(split="test"), "a test board"),
    (lambda r: r["x3"][0]["facts"].update(board_split="test_interp"), "a test board"),
    (lambda r: r["x3"][0]["builds"]["main"]["rows"][0]["scores"].pop(e52.C_ARM), "lacks an arm"),
    (lambda r: r["x3"][0]["facts"].update(files_identical_across_builds=False), "a file differs"),
    (lambda r: r["x3"][0]["facts"].update(refusals_identical_across_builds=False), "a refusal differs"),
    (lambda r: r["x3"][0]["facts"]["per_build"].pop("prefix"), "raw facts are missing"),
])
def test_void_conditions(good_facts, breaker, needle):
    recs = _recs(good_facts)
    breaker(recs)
    res = e52.reading(recs, [], C)
    assert res["void"] and any(needle in w for w in res["why"]) and "VOID" in e52.text(res)


def test_a_gate_or_commit_problem_voids(good_facts):
    res = e52.reading(_recs(good_facts), ["the gate did not pass on x1 at the run commit"], C)
    assert res["void"] and res["why"] == ["the gate did not pass on x1 at the run commit"]


# ------------------------------------------------------------------------------------------------ (iii), (iv)
def test_iii_against_p_and_c_over_five_truths_and_the_test_like_four(good_facts):
    gain = {"half": 2.0, "late": 1.0, "new1": -1.0, "new2": 0.5, "e75": -3.0, "e725": 1.0, "x1": 0.0, "x2": 0.5,
            "x3": -0.5}
    recs = {t: [_rec(t, s, good_facts, arm_board=51.0 + gain[t]) for s in e52.SEEDS] for t in TASKS}
    res = e52.reading(recs, [], C)
    for tag, base in (("vs_P", 51.0), ("vs_C", 51.2)):
        b = res["iii"][tag]["main"]["board"]
        want = {t: g + 51.0 - base for t, g in gain.items()}
        assert {t: v["mean"] for t, v in b["five_truths"]["per_task"].items()} == pytest.approx(want)
        merged = np.mean([np.mean([want[n] for n in names]) for names in res["truths"].values()])
        assert b["five_truths"]["merged"] == pytest.approx(merged)
        assert set(b["test_like"]["per_task"]) == set(e52.TEST_LIKE)
        tl = np.mean([np.mean([want[n] for n in names]) for names in res["test_like_truths"].values()])
        assert b["test_like"]["merged"] == pytest.approx(tl)
        assert b["five_truths"]["worst_task"] == "e75"
    assert res["iii"]["vs_P"]["against"] == e52.P_ARM and res["iii"]["vs_C"]["against"] == e52.C_ARM
    assert set(res["iii"]["vs_P"]["prefix"]) == {"board", *(f"group:{g}" for g in e52.GROUP_NAMES),
                                                 *(f"skill:{m}" for m in e52.METRICS)}
    assert res["iv"]["half"]["main"]["seeds_gate_bites"] == 0
    assert res["reported"]["half"]["main"][f"scale_log_ratio:{e52.ARM}"] == pytest.approx(0.05)


# ------------------------------------------------------------------------------------------------ gate
E50_GOOD = {b: {"per_task": {"x1": {"mean": 0.2}}, "D_s": [0.2] * 10} for b in ("main", "prefix")}


def _gate(e50_off=0.0, floor=50.0, p_off=0.0, build="main"):
    old = {e52.P_ARM: np.full(10, 55.0)}
    scores = {}
    for b in e52.BUILDS:
        p = np.full(10, 55.0 + (p_off if b == "main" else 0.0))
        c = p + 0.2 + (e50_off if b == build else 0.0)
        scores[b] = {e52.harness.FLOOR: np.full(10, floor), e52.P_ARM: p, e52.C_ARM: c}
    return e52.gate_check("x1", "key", scores, old, [], E50_GOOD, [])


def test_gate_branches():
    assert _gate()["pass"]
    assert _gate(e50_off=0.009)["pass"] and _gate(e50_off=-0.009, build="prefix")["pass"]
    for build in e52.BUILDS:
        g = _gate(e50_off=0.011, build=build)
        assert not g["pass"] and any("E50's per-task mean" in w for w in g["why"])
    assert not _gate(floor=50.1)["pass"]
    assert not _gate(p_off=0.11)["pass"] and _gate(p_off=0.09)["pass"]
    assert not e52.gate_check("x1", "key", {}, None, ["no committed E47 file"], E50_GOOD, [])["pass"]
    assert not e52.gate_check("x1", "key", {}, {}, [], None, ["E50's reading is missing"])["pass"]


def test_d_s_check():
    groups = {"A": ["a1", "a2"], "B": ["b1"]}
    gates = {n: {"diffs": {b: [0.2] * 10 for b in e52.BUILDS}} for n in ("a1", "a2", "b1")}
    assert e52.d_s_check(gates, E50_GOOD, groups) == []
    gates["b1"]["diffs"]["prefix"] = [0.2] * 9 + [0.32]  # D_s = (0.2 + 0.32) / 2: 0.06 off
    why = e52.d_s_check(gates, E50_GOOD, groups)
    assert len(why) == 1 and why[0].startswith("prefix")
    gates["b1"]["diffs"]["prefix"] = [0.2] * 9 + [0.29]  # 0.045 off
    assert e52.d_s_check(gates, E50_GOOD, groups) == []
    assert e52.d_s_check(gates, None, groups)


def test_e50_primary_reads_the_committed_reading_shape(tmp_path):
    import json

    p = tmp_path / "reading.json"
    good = {"void": False, "commit": e52.E50_RUN_COMMIT, "arm": "fgw_copyrows",
            "iii": {"primary": {"against": "fgw_shape0x", "main": {"board": {"per_task": {}, "D_s": []}},
                                "prefix": {"board": {"per_task": {}, "D_s": []}}}}}
    p.write_text(json.dumps(good))
    got, why = e52.e50_primary(p)
    assert why == [] and set(got) == {"main", "prefix"}
    p.write_text(json.dumps({**good, "commit": "0" * 40}))
    assert e52.e50_primary(p)[1]
    assert e52.e50_primary(tmp_path / "missing.json") == (None, ["E50's reading is missing"])


def test_has_row(tmp_path):
    p = tmp_path / "TRIES.md"
    p.write_text("| E50 (one try) | embryo |\n")
    assert not e52.has_row(p)
    p.write_text("| E52 (one try; the copy-rows file at the team's pre-registered size) | embryo |\n")
    assert e52.has_row(p)


def test_the_arms_and_controls():
    assert e52.ARM == "fgw_copyrows_sized" and e52.CONTROLS == ("fgw_shape0x", "fgw_copyrows")
    assert [c.name for c in e52.candidates({"stage": {}, "file": {}, "refused": {}})] == [
        "floor", "fgw_shape0x", "fgw_copyrows", "fgw_copyrows_sized"]
    assert {cat for _, cat in e52.REFUSALS} == {"dropped_over_share", "count_outside_limits", "copy_too_small"}


# ------------------------------------------------------------------------------------------------ the launcher
SH = REPO / "experiments" / "coords" / "run_e52.sh"


def _launch(tmp_path, *args, trace=False):
    env = {k: v for k, v in os.environ.items() if not k.startswith("GIT_")}
    env.update(HOME=str(tmp_path), GIT_CEILING_DIRECTORIES=str(tmp_path.parent))
    return subprocess.run(["bash", *(["-x"] if trace else []), str(SH), *args], cwd=tmp_path, env=env,
                          capture_output=True, text=True, timeout=60)


def test_the_launcher_refuses_out_of_order_and_writes_an_exit_file(tmp_path):
    S = tmp_path / "v" / "scratch" / "e52"
    r = _launch(tmp_path)
    assert r.returncode == 2 and "usage" in r.stderr and (S / "usage.exit").read_text().strip() == "2"
    r = _launch(tmp_path, "gate", C)
    assert r.returncode == 96 and "step tests has not exited 0" in r.stderr
    assert (S / "gate.exit").read_text().strip() == "96"
    (S / "gate.exit").write_text("0\n")
    (S / "gate.head").write_text("e" * 40 + "\n")  # passed, but at another commit
    r = _launch(tmp_path, "run", C)
    assert r.returncode == 96 and (S / "run.exit").read_text().strip() == "96"
    r = _launch(tmp_path, "bogus", C)
    assert r.returncode == 97 and (S / "bogus.exit").read_text().strip() == "97"
    trace = _launch(tmp_path, "read", C, trace=True)  # the dry-run trace of a harmless path: no git, no Python
    assert trace.returncode == 96 and "finish 96" in trace.stderr and "uv run" not in trace.stderr


def test_the_launcher_parses_deletes_nothing_and_lists_its_tests():
    text = SH.read_text()
    assert ":?" not in text
    assert subprocess.run(["bash", "-n", str(SH)], capture_output=True).returncode == 0
    code = [line.split("#")[0] for line in text.splitlines()]
    assert not any(__import__("re").search(r"(^|[\s;&|(])(rm|rmdir|unlink|shred|truncate)\s", line) for line in code)
    for f in ("tests/test_copy_rows_sized.py", "tests/test_e52_copy_rows_sized_read.py", "tests/test_copy_rows.py",
              "tests/test_copy_rows_build.py", "tests/test_harness.py", "tests/test_hidden_counts.py"):
        assert f in text and (REPO / f).exists()
    assert "SCRIPT=experiments/coords/e52_copy_rows_sized_read.py" in text and "R=runs/e52" in text


# ------------------------------------------------------------------------------------------------ fix round 2
def test_reading_i_compares_the_four_expression_terms_as_reported_and_the_x_bytes_raw(arrays, good_facts):
    """Row E52's wording of reading (i) for the four expression terms is the code's: the terms as the scorer reports
    them, the X bytes compared raw. Equal reported terms do not hide an X byte difference."""
    doc = " ".join(e52.__doc__.split())
    assert ("its de_score, de_direction, mmd_u and variogram, as the scorer reports them, equal fgw_copyrows's exactly "
            "(its X bytes compared raw)") in doc
    X = arrays[2].copy()
    X[0, 0] = np.nextafter(X[0, 0], np.float32(np.inf))
    bad = _facts(arrays, "main", arm_X=X)  # the reported terms equal C's, one X byte differs
    assert all(bad["expression_equal_c"].values()) and e52.i_failures(bad) == ["x_bytes"]


def test_raw_facts_computes_every_shape_quantity_inside_the_build(arrays, monkeypatch):
    """Lesser point 1: the scored parity and d2_shape, like the frame, the voxels and the Dice, are computed under the
    scorer build they are recorded for."""
    active, calls = [], []

    @e52.contextlib.contextmanager
    def scorer_build(build):
        active.append(build)
        try:
            yield
        finally:
            active.pop()

    orig_parity, orig_match = ev.svd_parity, ev.match_rows
    monkeypatch.setattr(ev, "scorer_build", scorer_build)
    monkeypatch.setattr(ev, "svd_parity", lambda Z: calls.append(("parity", tuple(active))) or orig_parity(Z))
    monkeypatch.setattr(ev, "match_rows", lambda *a: calls.append(("match", tuple(active))) or orig_match(*a))
    sm = ev._shape_metrics()
    orig_d2 = sm.d2_distance
    monkeypatch.setattr(sm, "d2_distance", lambda *a, **k: calls.append(("d2", tuple(active))) or orig_d2(*a, **k))
    _facts(arrays, "prefix")
    assert {k for k, _ in calls} == {"parity", "match", "d2"}
    assert all(where == ("prefix",) for _, where in calls), calls


def _gate_files(runs, commit=C, diff=0.2, passed=True):
    for name in TASKS:
        p = e52._gate_path(runs, name)
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(json.dumps({"pass": passed, "commit": commit, "diffs": {b: [diff] * 10 for b in e52.BUILDS}}))


@pytest.mark.parametrize("e50_d_s, passed, want", [(0.2, True, 0), (0.24, True, 0), (0.26, True, 3),
                                                     (0.2, False, 3)])
def test_gate_all_exits_3_on_a_d_s_miss_and_run_is_then_refused(tmp_path, monkeypatch, capsys, e50_d_s, passed,
                                                                 want):
    """Blocker 1: the nine-task gate, D_s included, is its own step (gate-all) that exits 3 on a miss; cmd_run keeps
    the same check as a second guard."""
    runs = tmp_path / "runs"
    _gate_files(runs, passed=passed)
    e50 = {b: {"per_task": {}, "D_s": [e50_d_s] * 10} for b in e52.BUILDS}
    monkeypatch.setattr(e52, "_refuse_unless_ready", lambda commit: None)
    monkeypatch.setattr(e52, "e50_primary", lambda path=None: (e50, []))
    with pytest.raises(SystemExit) as ex:
        e52.main(["gate-all", "--commit", C, "--runs", str(runs)])
    assert ex.value.code == want
    out = capsys.readouterr().out.strip().splitlines()
    assert out == [f"gate all {'pass' if want == 0 else 'fail'}"]
    if want:
        monkeypatch.setattr(e52, "one_seed", lambda *a: pytest.fail("cmd_run went past its gate check"))
        with pytest.raises(SystemExit) as ex:
            e52.main(["run", "--task", "x1", "--commit", C, "--runs", str(runs), "--data", str(tmp_path)])
        assert str(ex.value.code).startswith("refused (the pre-run gate)")


def test_gate_all_exits_1_when_it_may_not_run(tmp_path, monkeypatch):
    monkeypatch.setattr(e52, "_refuse_unless_ready", lambda commit: e52.sys.exit("refused: not here"))
    with pytest.raises(SystemExit) as ex:
        e52.main(["gate-all", "--commit", C, "--runs", str(tmp_path)])
    assert ex.value.code == "refused: not here"  # Python exits 1 on a message


FAKE_SCRIPT = """import os, sys
cmd = sys.argv[1]
if cmd == "gate":
    print("gate " + sys.argv[sys.argv.index("--task") + 1] + " pass")
    sys.exit(0)
if cmd == "gate-all":
    print("gate all " + ("pass" if os.environ["FAKE_GATE_ALL"] == "0" else "fail"))
    sys.exit(int(os.environ["FAKE_GATE_ALL"]))
sys.exit(0)
"""


def _fake_e52_checkout(tmp_path):
    """A temporary repository holding run_e52.sh and a stand-in read script (the per-task gates pass; gate-all exits
    FAKE_GATE_ALL), committed and clean, with the tests step recorded as passed at its HEAD."""
    from conftest import git

    repo = tmp_path / "checkout"
    (repo / "experiments" / "coords").mkdir(parents=True)
    (repo / "experiments" / "coords" / "run_e52.sh").write_text(SH.read_text())
    (repo / "experiments" / "coords" / "e52_copy_rows_sized_read.py").write_text(FAKE_SCRIPT)
    (repo / ".gitignore").write_text("runs/\n")
    git(repo, "init", "-q")
    git(repo, "add", "-A")
    git(repo, "commit", "-qm", "e52 stand-in")
    head = git(repo, "rev-parse", "HEAD")
    home = tmp_path / "home"
    S = home / "v" / "scratch" / "e52"
    S.mkdir(parents=True)
    (S / "tests.exit").write_text("0\n")
    (S / "tests.head").write_text(head + "\n")
    bin_ = tmp_path / "bin"
    bin_.mkdir()
    uv = bin_ / "uv"  # `uv run --frozen --no-sync python ARGS` -> this interpreter with ARGS
    uv.write_text(f'#!/bin/sh\nshift 4\nexec "{sys.executable}" "$@"\n')
    uv.chmod(0o755)
    return repo, head, home, S, bin_


@pytest.mark.parametrize("gate_all, want", [("3", 3), ("0", 0), ("1", 1), ("2", 1)])
def test_the_gate_step_ends_with_gate_all_and_a_d_s_miss_gives_3_and_refuses_run(tmp_path, gate_all, want):
    """Blocker 1 at the launcher: all nine task gates pass, then gate-all decides the gate step's exit; a D_s miss
    (gate-all 3) gives gate.exit 3 and the run step is refused (96)."""
    from conftest import git_env

    repo, head, home, S, bin_ = _fake_e52_checkout(tmp_path)
    env = {**git_env(), "HOME": str(home), "PATH": f"{bin_}:{os.environ['PATH']}", "FAKE_GATE_ALL": gate_all,
           "GIT_CEILING_DIRECTORIES": str(tmp_path)}

    def step(name):
        return subprocess.run(["bash", "experiments/coords/run_e52.sh", name, head], cwd=repo, env=env,
                              capture_output=True, text=True, timeout=120)

    r = step("gate")
    assert r.returncode == want, r.stderr
    assert (S / "gate.exit").read_text().strip() == str(want)
    assert (repo / "runs/e52/exit/gate-all").read_text().strip() == gate_all
    assert (repo / "runs/e52/exit/gate-all.other_lines").read_text().strip() == "0"
    assert all((repo / f"runs/e52/exit/gate-{t}").read_text().strip() == "0" for t in TASKS)
    r = step("run")
    if want == 0:
        assert r.returncode == 0 and (S / "run.exit").read_text().strip() == "0"
    else:
        assert r.returncode == 96 and "step gate has not exited 0" in r.stderr
        assert (S / "run.exit").read_text().strip() == "96"


def test_the_launcher_lists_the_e50_test_files_and_the_build_launcher_lists_the_read_tests():
    """Lesser point 3."""
    assert "tests/test_shape_interp_build.py" in SH.read_text()
    build = (REPO / "experiments" / "builds" / "embryo_copy_rows_sized.sh").read_text()
    assert "tests/test_e52_copy_rows_sized_read.py" in build
