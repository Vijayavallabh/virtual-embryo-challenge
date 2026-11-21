"""Row E53's kill-test script (experiments/coords/e53_geodesic_read.py) and its launcher (run_e53.sh), on synthetic
arrays, hand-made result records and committed practice files only (E47's score files, E48's gate files, E49's
reading). Nothing is read from data/, nothing is built from a stage and nothing is scored.

Pinned: the constants of the gate and of the pass line; the refusal without row E53; the target-free gate's clauses
g1 to g3 (the test pair decides, every pair is required at the run commit); J7's dry run; the reproduction gate (C0
against E47 by E49's own check, E49c against E49's committed per-task mean and SD, 0.01 each); the pass line (a) to (f)
each able to stop the arm, build main deciding; the void list; the read never overwrites; the launcher's order, exit
files and refusal paths under a temporary HOME, and that it removes nothing."""

from __future__ import annotations

import importlib.util
import json
import os
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
SEEDS = list(range(10))
COMMIT = "a" * 40
LAUNCHER = REPO / "experiments" / "coords" / "run_e53.sh"
C0, E49C, ARM, PLA = "fgw_shape0x", "fgw_shape0x_ov", "fgw_shape0x_geo", "fgw_shape0x_geo_rot"
METRICS = ("de_score", "de_direction", "mmd_u", "variogram", "d2_shape", "occupancy_dice", "scale_log_ratio",
           "neighborhood_mmd")


@pytest.fixture(scope="module")
def e53():
    spec = importlib.util.spec_from_file_location("e53_geodesic_read",
                                                  REPO / "experiments" / "coords" / "e53_geodesic_read.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


# ------------------------------------------------------------------------------------------------ constants, row
def test_the_constants_are_the_rows(e53):
    assert (e53.C0, e53.E49C, e53.ARM, e53.PLACEBO) == (C0, E49C, ARM, PLA)
    assert e53.TRY_IDS == {ARM: "E53a", PLA: "E53b"} and e53.ROW_IDS == ("E53a", "E53")
    assert e53.FORM_OF[ARM] == "geodesic_overlay" and e53.FORM_OF[PLA] == "geodesic_overlay_rot"
    assert e53.FORM_OF[E49C] == "overlay" and e53.BUILDS == ("main", "prefix")
    assert e53.TEST_LIKE == ("half", "late", "new1", "new2")
    assert (e53.A_MERGED, e53.A_TRUTHS, e53.A_LEAVE_ONE_OUT, e53.B_TEST_LIKE) == (1.0, 4, 0.5, 1.0)
    assert (e53.C_OVER_C0, e53.C_OVER_E49C, e53.D_LOCAL, e53.E_S, e53.E_AHEAD, e53.F_MAX_REFUSED) == (
        0.04, 0.02, 0.01, 0.03, 8, 5)
    assert (e53.G1_RATIO, e53.G1_DRAWS, e53.G2_JAC, e53.G3_DRAWS) == (0.5, 8, 0.05, 8)
    assert e53.DESCRIBE_SEEDS == tuple(range(1_000_001, 1_000_011)) and e53.TEST_T == 0.5
    assert e53.TEST_PAIR == {"before": "T2/E7.5.h5ad", "after": "T2/E8.0.h5ad"}
    assert e53.PAIRS == ("test", "half", "late", "new1", "new2", "e75", "e725", "x1", "x2", "x3")
    assert e53.REPRO_TOL == 0.01 and (e53.EXPR_TOL, e53.SCALE_TOL) == (1e-9, 1e-6)


def test_a_run_is_refused_without_row_e53(e53, tmp_path, monkeypatch):
    tries = tmp_path / "TRIES.md"
    tries.write_text("| E50 (one try) | embryo |\n| E51a and E51b | embryo |\n")
    monkeypatch.setattr(e53, "TRIES", tries)
    with pytest.raises(SystemExit, match="E53 has no row"):
        e53._refuse_unless_ready(COMMIT)
    for first in ("E53a and E53b (one row)", "E53 (one row)"):
        tries.write_text(f"| {first} | embryo |\n")
        assert e53.has_row(tries)


def test_the_task_fractions_follow_the_stage_days(e53):
    tasks = e53.nt.all_tasks()
    for name, t in tasks.items():
        f = e53.task_fraction(t)
        assert 0.0 < f < 1.0, name


# ------------------------------------------------------------------------------------------------ target-free gate
def _draw(ratio=0.3, jac=0.5, thin=0.3, lo=0.2, hi=0.4):
    return {"seed": 0, "data_ratio": ratio, "jac_min": jac, "thin_sd": thin, "thin_sd_before": lo,
            "thin_sd_after": hi, "thin_between": lo <= thin <= hi, "lbfgs_n_iter": 150, "p0_norm_mean": 0.1}


def test_the_gate_clauses(e53):
    good = [_draw() for _ in range(10)]
    assert e53.gate_clauses(good)["passes"]
    g1 = [_draw(ratio=0.6) for _ in range(3)] + [_draw() for _ in range(7)]
    r = e53.gate_clauses(g1)
    assert not r["clauses"]["g1"] and r["clauses"]["g2"] and r["clauses"]["g3"]
    assert e53.gate_clauses([_draw(ratio=0.6) for _ in range(2)] + [_draw() for _ in range(8)])["passes"]
    assert e53.gate_clauses([_draw(ratio=0.5)] * 10)["passes"]  # "at most 0.5"
    r = e53.gate_clauses([_draw(jac=0.05)] + [_draw() for _ in range(9)])  # "above 0.05", every draw
    assert not r["clauses"]["g2"] and r["clauses"]["g1"]
    r = e53.gate_clauses([_draw(thin=0.1) for _ in range(3)] + [_draw() for _ in range(7)])
    assert not r["clauses"]["g3"]
    assert e53.gate_clauses([_draw(thin=0.45, lo=0.5, hi=0.4)] * 10)["clauses"]["g3"]  # either order
    r = e53.gate_clauses([{"seed": 0, "refused": "x"}] + [_draw() for _ in range(9)])
    assert not r["clauses"]["g2"] and r["refused_draws"] == 1


def _describe(tmp_path, commit=COMMIT, test_draws=None, skip=()):
    for p in ("test", "half", "late", "new1", "new2", "e75", "e725", "x1", "x2", "x3"):
        if p in skip:
            continue
        draws = test_draws if (p == "test" and test_draws is not None) else [_draw() for _ in range(10)]
        draws = [{**d, "seed": 1_000_001 + i} for i, d in enumerate(draws)]
        f = tmp_path / "describe" / f"{p}.json"
        f.parent.mkdir(parents=True, exist_ok=True)
        f.write_text(json.dumps({"commit": commit, "t": 0.5, "draws": draws}))


def test_the_pregate_decides_on_the_test_pair_and_needs_every_pair(e53, tmp_path):
    _describe(tmp_path)
    assert e53.pregate(tmp_path, COMMIT)["pass"]
    assert not e53.pregate(tmp_path, "b" * 40)["pass"]
    bad = tmp_path / "bad"
    _describe(bad, test_draws=[_draw(jac=0.01)] + [_draw() for _ in range(9)])
    assert not e53.pregate(bad, COMMIT)["pass"]
    miss = tmp_path / "miss"
    _describe(miss, skip=("x3",))
    r = e53.pregate(miss, COMMIT)
    assert not r["pass"] and any("x3" in w for w in r["why"])
    task_bad = tmp_path / "taskbad"  # a task pair's values are reported, never deciding
    _describe(task_bad)
    p = task_bad / "describe" / "late.json"
    d = json.loads(p.read_text())
    d["draws"] = [{**x, "jac_min": 0.0} for x in d["draws"]]
    p.write_text(json.dumps(d))
    r = e53.pregate(task_bad, COMMIT)
    assert r["pass"] and not r["pairs"]["late"]["passes"]


def test_the_gate_values_come_from_the_forms_record(e53):
    rec = {"shape_interp_lddmm_data_ratio": 0.4, "shape_interp_jac_min": 0.2, "shape_interp_thin_sd": 0.3,
           "shape_interp_thin_sd_before": 0.2, "shape_interp_thin_sd_after": 0.35, "shape_interp_thin_between": True,
           "shape_interp_lddmm_lbfgs_n_iter": 150, "shape_interp_lddmm_p0_norm_mean": 0.01}
    v = e53.gate_values(rec)
    assert v["data_ratio"] == 0.4 and v["jac_min"] == 0.2 and v["thin_between"]


# ------------------------------------------------------------------------------------------------ J7, reproduction
def test_j7_reads_e48s_committed_rescore(e53, tmp_path):
    r = e53.j7()
    assert set(r) >= {"worst", "within", "tolerance", "why"} and r["tolerance"] == 0.01 and r["why"] == []
    for name in e53.nt.all_tasks():
        (tmp_path / f"{name}.json").write_text(json.dumps({"worst": {"mean": 0.0, "seed": 0.0}}))
    assert e53.j7(tmp_path)["within"]
    (tmp_path / "x3.json").write_text(json.dumps({"worst": {"mean": 0.002, "seed": 0.011}}))
    assert not e53.j7(tmp_path)["within"]
    (tmp_path / "x3.json").write_text(json.dumps({"worst": {"mean": None, "seed": None}}))
    assert not e53.j7(tmp_path)["within"] and e53.j7(tmp_path)["why"]


def test_e49c_is_checked_against_e49s_committed_mean_and_sd(e53, tmp_path):
    old, why = e53.e49c_committed("late")
    assert why == [] and set(old) >= {"mean", "sd"}
    rng = np.random.default_rng(0)
    d = rng.normal(size=10)
    ref = {"mean": float(d.mean()), "sd": float(d.std(ddof=1))}
    assert e53.e49c_check("late", d, ref, [])["pass"]
    assert not e53.e49c_check("late", d + 0.011, ref, [])["pass"]
    moved = d.copy()
    moved[0] += 0.2
    moved[1] -= 0.2  # same mean, another SD
    assert not e53.e49c_check("late", moved, ref, [])["pass"]
    assert not e53.e49c_check("late", d, None, ["no line"])["pass"]
    void = tmp_path / "r.json"
    void.write_text(json.dumps({"void": True}))
    assert e53.e49c_committed("late", void)[0] is None


def test_c0_is_checked_against_e47_by_e49s_own_gate(e53):
    key = e53.nt.all_tasks()["e75"].key
    old, why = e53.e49.e47_scores(key)
    assert why == []
    assert e53.e49.gate_check("e75", key, {c: old[c].copy() for c in (e53.harness.FLOOR, C0)})["pass"]


def test_gate_problems_need_the_pregate_j7_and_nine_gate_files(e53, tmp_path):
    assert len(e53.gate_problems(tmp_path, COMMIT)) == 11
    (tmp_path / "pregate.json").write_text(json.dumps({"commit": COMMIT, "pass": True}))
    (tmp_path / "j7.json").write_text(json.dumps({"commit": COMMIT, "within": False}))
    for name in e53.nt.all_tasks():
        p = e53._gate_path(tmp_path, name)
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(json.dumps({"pass": True, "commit": COMMIT}))
    assert e53.gate_problems(tmp_path, COMMIT) == []
    (tmp_path / "pregate.json").write_text(json.dumps({"commit": COMMIT, "pass": False}))
    assert e53.gate_problems(tmp_path, COMMIT) == ["the target-free gate failed"]


# ------------------------------------------------------------------------------------------------ the reading
NAMES = ("floor", C0, E49C, ARM, PLA)
GOOD = {"board": {ARM: 1.2}, "shape": {ARM: 0.05, E49C: 0.02}, "S": {ARM: 0.04}, "local": {}}


def _row(seed: int, sh: dict, task: str, tshift: dict, seed_s: dict) -> dict:
    scores, skills, groups = {}, {}, {}
    for n in NAMES:
        b = 50.0 if n == "floor" else 47.0 + 0.01 * seed + sh["board"].get(n, 0.0) + tshift.get(task, {}).get(n, 0.0)
        sk = dict.fromkeys(METRICS, 0.5)
        sk["d2_shape"] = sk["occupancy_dice"] = 0.5 + sh["S"].get(n, 0.0) + seed_s.get(seed, {}).get(n, 0.0)
        scores[n], skills[n] = b, sk
        groups[n] = {"shape_growth": 0.4 + sh["shape"].get(n, 0.0), "expression_change": 0.5, "cell_state": 0.5,
                     "local_organisation": 0.4 + sh["local"].get(n, 0.0)}
    return {"seed": seed, "scores": scores, "skills": skills, "groups": groups}


def _facts(refused=None):
    names = (C0, E49C, ARM, PLA)
    recs = {n: {"n_jittered": 0, "fallback": False} for n in names}
    recs[ARM].update({"shape_interp_lddmm_data_ratio": 0.3, "shape_interp_jac_min": 0.4, "shape_interp_thin_sd": 0.3})
    return {"board_split": "pseudo", "x": dict.fromkeys(names, "x"), "n": dict.fromkeys(names, 100),
            "refused": {n: (refused or {}).get(n) for n in names[1:]}, "records": recs,
            "hands": {n: {"prefix_mode": "proper", "dice_proper": 0.8, "dice_mirror": 0.7, "dice_gap": 0.1}
                      for n in names},
            "files_identical_across_builds": True, "refusals_identical_across_builds": True}


def _recs(e53, sh=None, tshift=None, seed_s=None, refused_on=None, prefix_board=None):
    sh = {k: dict(v) for k, v in (sh or GOOD).items()}
    out = {}
    for name in e53.nt.all_tasks():
        rs = []
        for s in SEEDS:
            row = _row(s, sh, name, tshift or {}, seed_s or {})
            builds = {b: {"provenance": {"commit": COMMIT, "tree_clean": True}, "build": b, "protocol": "bake",
                          "board": {"split": "pseudo"}, "rows": [json.loads(json.dumps(row))]} for b in e53.BUILDS}
            if prefix_board is not None:
                builds["prefix"]["rows"][0]["scores"][ARM] += prefix_board
            rs.append({"seed": s, "run_commit": COMMIT, "builds": builds,
                       "facts": _facts((refused_on or {}).get((name, s)))})
        out[name] = rs
    return out


def test_a_clean_pass_names_the_file_and_build_main_decides(e53):
    res = e53.reading(_recs(e53, prefix_board=-5.0), [], COMMIT)
    assert not res["void"] and res["passes"] and all(res["clauses"].values())
    assert res["D"]["board_e49c"]["mean"] == pytest.approx(1.2) and res["truths_positive"] == 5
    assert "priced question beside E49c" in res["outcome"]
    assert res["reported"]["prefix_board_vs_e49c"]["mean"] == pytest.approx(-3.8)
    t = e53.text(res)
    assert "(a) holds" in t and "(f) holds" in t and "VOID" not in t and "e75:" in t and "e725:" in t


def _fails(e53, key, **kw):
    res = e53.reading(_recs(e53, **kw), [], COMMIT)
    assert not res["void"]
    assert not res["clauses"][key] and not res["passes"] and res["outcome"].endswith("no file")
    return res


def test_a_fails_on_the_merged_bar(e53):
    sh = {**GOOD, "board": {ARM: 0.99}}
    _fails(e53, "a", sh=sh)


def test_a_fails_on_two_negative_truths(e53):
    sh = {**GOOD, "board": {}}
    tshift = {t: {ARM: 3.0} for t in ("half", "late", "new1", "x1", "new2", "x2", "x3")}
    tshift.update({"e75": {ARM: -1.0}, "e725": {ARM: -1.0}})
    res = _fails(e53, "a", sh=sh, tshift=tshift)
    assert res["D"]["board_e49c"]["mean"] >= 1.0 and min(res["leave_one_out_vs_e49c"].values()) >= 0.5
    assert res["truths_positive"] == 3


def test_a_fails_when_one_truth_carries_the_gain(e53):
    sh = {**GOOD, "board": {ARM: 0.1}}
    res = _fails(e53, "a", sh=sh, tshift={"e75": {ARM: 4.6}})
    assert res["D"]["board_e49c"]["mean"] >= 1.0 and res["truths_positive"] == 5


def test_b_fails_on_the_test_like_subset(e53):
    sh = {**GOOD, "board": {ARM: 0.5}}
    res = _fails(e53, "b", sh=sh, tshift={t: {ARM: 1.5} for t in ("x1", "x2", "x3", "e75", "e725")})
    assert res["clauses"]["a"]


def test_c_fails_against_c0_and_against_e49c(e53):
    _fails(e53, "c", sh={**GOOD, "shape": {ARM: 0.035}})
    _fails(e53, "c", sh={**GOOD, "shape": {ARM: 0.05, E49C: 0.04}})


def test_d_fails_a_local_loss_over_001(e53):
    _fails(e53, "d", sh={**GOOD, "local": {ARM: -0.011}})
    assert e53.reading(_recs(e53, sh={**GOOD, "local": {ARM: -0.009}}), [], COMMIT)["clauses"]["d"]


def test_e_fails_on_the_mean_and_on_the_seeds(e53):
    _fails(e53, "e", sh={**GOOD, "S": {ARM: 0.04, PLA: 0.015}})
    res = _fails(e53, "e", sh={**GOOD, "S": {ARM: 0.2}}, seed_s={s: {ARM: -0.21} for s in (0, 1, 2)})
    assert res["D"]["S_placebo"]["mean"] >= 0.03 and res["D"]["S_placebo"]["ahead"] == 7


def test_f_allows_five_refusals_of_ninety(e53):
    five = {(t, s): {ARM: "cloud guard"} for t, s in (("half", 0), ("late", 1), ("x1", 2), ("x2", 3), ("e75", 4))}
    res = e53.reading(_recs(e53, refused_on=five), [], COMMIT)
    assert res["clauses"]["f"] and res["arm_refused_task_seeds"] == 5
    six = {**five, ("e725", 5): {ARM: "shape: x"}}
    _fails(e53, "f", refused_on=six)
    placebo = {(t, 0): {PLA: "cloud guard"} for t in e53.nt.all_tasks()}  # the placebo's refusals are reported
    res = e53.reading(_recs(e53, refused_on=placebo), [], COMMIT)
    assert res["clauses"]["f"] and res["reported"]["refused_task_seeds"][PLA] == 9


@pytest.mark.parametrize("breach", ["seed", "commit", "unclean", "expression", "count", "expr_skill", "scale",
                                    "builds", "refusal_builds", "gate", "test_board", "missing_arm", "e49c_expr"])
def test_the_void_conditions(e53, breach):
    recs = _recs(e53)
    why = []
    r = recs["late"][3]
    f = r["facts"]
    if breach == "seed":
        recs["late"].pop(3)
    elif breach == "commit":
        r["builds"]["prefix"]["provenance"]["commit"] = "b" * 40
    elif breach == "unclean":
        r["builds"]["main"]["provenance"]["tree_clean"] = False
    elif breach == "expression":
        f["x"][ARM] = "y"
    elif breach == "count":
        f["n"][PLA] = 99
    elif breach == "expr_skill":
        r["builds"]["main"]["rows"][0]["skills"][ARM]["mmd_u"] += 1e-8
    elif breach == "scale":
        r["builds"]["main"]["rows"][0]["skills"][ARM]["scale_log_ratio"] += 1e-5
    elif breach == "builds":
        f["files_identical_across_builds"] = False
    elif breach == "refusal_builds":
        f["refusals_identical_across_builds"] = False
    elif breach == "gate":
        why = ["the gate did not pass on late at the run commit"]
    elif breach == "test_board":
        f["board_split"] = "test_interp"
    elif breach == "missing_arm":
        r["builds"]["main"]["rows"][0]["scores"].pop(PLA)
    else:
        f["x"][E49C] = "y"
    res = e53.reading(recs, why, COMMIT)
    assert res["void"] and res["why"]
    assert e53.text(res).startswith(f"Row E53 (run commit {COMMIT}): VOID")


def test_a_jittered_seed_is_exempt_from_the_scale_check(e53):
    recs = _recs(e53)
    r = recs["late"][3]
    r["builds"]["main"]["rows"][0]["skills"][ARM]["scale_log_ratio"] += 1e-5
    r["facts"]["records"][ARM]["n_jittered"] = 2
    assert not e53.reading(recs, [], COMMIT)["void"]


def test_the_read_refuses_to_overwrite(e53, tmp_path):
    (tmp_path / "reading.json").write_text("x")
    with pytest.raises(SystemExit, match="refusing to overwrite"):
        e53.main(["read", "--commit", COMMIT, "--runs", str(tmp_path), "--out", str(tmp_path)])


def test_the_wrappers_are_e49s_and_a_test_board_stops_them(e53):
    from types import SimpleNamespace

    from vec import evaluate as ev
    from vec.boards import Board

    board = Board(key="T2:embryo:test_interp", task="T2", genes=["a"], needs_coords=True, split="test_interp",
                  min_cells=1, max_cells=50, truth_cells=10)
    stage = ev.Stage(X=np.ones((5, 1), np.float32), coords=np.zeros((5, 3)))
    with e53.e49.captured() as seen:
        seen["board"] = board
        cands = e53.candidates(seen)
        assert [c.name for c in cands] == ["floor", C0, E49C, ARM, PLA]
        arm = cands[3]
        with pytest.raises(RuntimeError, match="test board"):
            arm.build(SimpleNamespace(seed=1_000_003))
    assert stage.n == 5


# ------------------------------------------------------------------------------------------------ the launcher
def test_the_launcher_orders_the_steps_caps_threads_and_lists_its_tests():
    sh = LAUNCHER.read_text()
    for need in ("OMP_NUM_THREADS=8", "OPENBLAS_NUM_THREADS=8", "MKL_NUM_THREADS=8", "NUMPY_MADVISE_HUGEPAGE=0",
                 "finish 98", "finish 99", "finish 96", "finish 95", "finish 97", '"$S/$STEP.exit"',
                 "uv run --frozen --no-sync python -I -B", "worktree", "single pregate", "single j7"):
        assert need in sh, need
    order = dict(re.findall(r"^  ([a-z0-9]+)\) PREV=\"([a-z0-9]*)\" ;;$", sh, flags=re.M))
    assert order == {"tests": "", "describe": "tests", "j7": "describe", "gate": "j7", "run": "gate", "read": "run"}
    tasks = sh.split("TASKS=(")[1].split(")")[0].split()
    assert sorted(tasks) == sorted(["half", "late", "new1", "new2", "e75", "e725", "x1", "x2", "x3"])
    assert 'PAIRS=(test "${TASKS[@]}")' in sh
    listed = re.findall(r"tests/test_[a-z0-9_]+\.py", sh)
    assert {"tests/test_lddmm.py", "tests/test_geodesic_overlay.py", "tests/test_e53_geodesic_read.py"} <= set(listed)
    assert all((REPO / f).exists() for f in listed)
    assert "parallel gate --task cpu" in sh and "parallel run --task gpu" in sh and "parallel describe --pair gpu" in sh


def test_the_launcher_removes_nothing_and_has_no_quote_in_a_usage_message():
    sh = LAUNCHER.read_text()
    code = [line.split("#")[0] for line in sh.splitlines()]
    assert not any(re.search(r"(^|[\s;&|(])(rm|rmdir|unlink|shred|truncate|mv)\s", line) for line in code)
    for msg in re.findall(r"\$\{\d+:\?([^}]*)\}", sh):
        assert "'" not in msg and '"' not in msg, msg


def _run(home: Path, *argv: str) -> subprocess.CompletedProcess:
    env = {k: v for k, v in os.environ.items() if k not in ("GIT_DIR", "GIT_WORK_TREE", "GIT_INDEX_FILE")}
    env["HOME"] = str(home)
    return subprocess.run(["bash", str(LAUNCHER), *argv], cwd=REPO, env=env, capture_output=True, text=True)


def test_the_launchers_refusal_paths_write_exit_files_and_start_nothing(tmp_path):
    scratch = tmp_path / "v" / "scratch" / "e53"
    out = _run(tmp_path)
    assert out.returncode == 2 and (scratch / "usage.exit").read_text().strip() == "2"
    out = _run(tmp_path, "bogus", "0" * 40)
    assert out.returncode == 97 and (scratch / "bogus.exit").read_text().strip() == "97"
    for step, prev in (("describe", "tests"), ("j7", "describe"), ("gate", "j7"), ("run", "gate"), ("read", "run")):
        out = _run(tmp_path, step, "0" * 40)
        assert out.returncode == 96 and f"step {prev} has not exited 0" in out.stderr, step
        assert (scratch / f"{step}.exit").read_text().strip() == "96"
    out = _run(tmp_path, "tests", "0" * 40)  # never the real HEAD: an unclean tree (98) or another HEAD (99)
    assert out.returncode in (98, 99) and (scratch / "tests.exit").read_text().strip() == str(out.returncode)
    assert not (scratch / "tests.log").exists()
    # a previous step that exited 0 at another commit does not open the next
    (scratch / "tests.exit").write_text("0\n")
    (scratch / "tests.head").write_text("1" * 40 + "\n")
    out = _run(tmp_path, "describe", "0" * 40)
    assert out.returncode == 96
    (scratch / "tests.head").write_text("0" * 40 + "\n")
    out = _run(tmp_path, "describe", "0" * 40)
    assert out.returncode in (98, 99)  # past the order check, stopped at the tree's own checks
    assert sorted(p.name for p in scratch.iterdir() if p.name.endswith(".log")) == []
