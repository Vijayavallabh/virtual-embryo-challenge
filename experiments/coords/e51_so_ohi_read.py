#!/usr/bin/env python
"""Row E51 of experiments/harness/TRIES.md (one row, one script; placed by fresh analysts before any run): the shared
practice read of two embryo candidates on the nine practice tasks, as the score-blind challenger's specification
(job_embryo_so_ohi) fixes it. No number here comes from a test score.

Arms (harness candidates, vec.harness.SO_OHI_CANDIDATES), both counted:
  SO    fgw_shape0x_series_ov  the overlay form's positions (E49c) with the series expression (E47b)
  O_hi  fgw_shape0x_ov_hi      the overlay form's positions with the depth top-up render_hi (E19c, E24b)
Controls, rescored and not counted: fgw_shape0x_series (E47b) and fgw_shape0x_ov (E49c). Reported only: fgw (for the
increment over the best of fgw, fgw_shape0x_series and fgw_shape0x_ov) and the floor. fgw_shape0_hi (E24b; committed
by E24 on half, late, new1, new2, e75, e725 and by E28 on x1, x2, x3) is run in the same call as O_hi's expression
reference, not counted, and decides only O_hi's void.

Path. harness.run_task on each task of nine_tasks.all_tasks(), seeds 0-9, bake protocol, build seed
BUILD_SEED_OFFSET + seed, under build "main" (primary) and "prefix" (reported only). No refusal is caught: an arm that
cannot build stops its task, and the missing seeds void the read.

Steps (each refuses unless row E51 is in TRIES.md, the tree is clean and HEAD is the run commit):
  gate --task NAME --commit SHA  the floor, fgw_shape0x, fgw_shape0x_series and fgw_shape0x_ov under main, ten seeds.
                                 Pass: fgw_shape0x_series within 0.01 board points of E47's committed scores per seed
                                 and per ten-seed mean; fgw_shape0x_ov's ten-seed mean minus fgw_shape0x's within 0.01
                                 of E49's committed per-task mean (the only committed E49 figure of that arm: E49's
                                 reading keeps no per-seed score). Exit 3 on a failure (stop; nothing is counted). A
                                 gate file of this commit is read again, not re-run.
  run --task NAME --commit SHA   the arms, controls and reported arms under both builds, seed by seed (resumable: a seed
                                 whose file exists is skipped); refused unless all nine gate files pass at this commit
  read --commit SHA              the reading into experiments/coords/results/e51_so_ohi/ (reading.json, reading.txt),
                                 never overwriting
The script prints only "<task> seed <s> done", "wrote <path>" and "gate <task> pass|fail". It reads no test board, test
file or test anchor, and prints no count.

Reading, fixed in row E51 before any output (applied by `reading`). J4 merge: five truths at equal weight; D_s(X) for
a difference X at seed s is the mean over the truths of each truth's task-mean difference at s; the merged mean is
mean_s D_s. The test-like series is the same over half, late, new1, new2 (three truths). Build main.
  SO passes only if all hold: (a) the test-like merged board against fgw_shape0x_series is above +0.3 with D_s above 0
     in at least 8 of 10 seeds; (b) the five-truth merged board against fgw_shape0x_series is above 0 with D_s above 0
     in at least 8 of 10 seeds; (c) no group is behind fgw_shape0x_series (E49's P4: merged mean below 0 with at least
     8 of 10 D_s below 0).
  O_hi passes only if all hold: (a) the five-truth merged board against fgw_shape0x_ov is above +0.3 with D_s above 0
     in at least 8 of 10 seeds; (b) the test-like merged board against fgw_shape0x_ov is above 0; (c) no group is
     behind fgw_shape0x_ov (E49's P4).
  A fail means no file. A pass only puts the question to the captain. If both pass, both go to him; the read does not
  pick between them.
Void, for both arms: a missing task or seed; a control or reported arm missing; another run commit or an unclean
tree; a gate file missing or failing; a test board. Void for SO only: SO missing; its four expression skills (de_score,
de_direction, mmd_u, variogram) off fgw_shape0x_series' by more than 1e-9; its d2_shape, occupancy_dice or
scale_log_ratio skill off fgw_shape0x_ov's by more than 1e-6. Void for O_hi only (SO's read stands): O_hi missing;
fgw_shape0_hi missing in any task, seed or build; O_hi's d2_shape, occupancy_dice or scale_log_ratio skill off
fgw_shape0x_ov's by more than 1e-6; O_hi's four expression skills (de_score, de_direction, mmd_u, variogram) off the
same run's fgw_shape0_hi's by more than 1e-9 in any task, seed or build (same drawn rows, same render_hi path, same
process and thread count); the same run's fgw_shape0_hi board score under main off its committed score (E24 for half,
late, new1, new2, e75, e725; E28 for x1, x2, x3) by more than 0.01 board points in any seed or in the ten-seed mean,
or that committed score missing (rule 9 J7, nine_tasks.TOL: the rendered arms move by up to about 1e-3 between runs
with the BLAS thread count, so the specification's cross-run 1e-9 on expression skills is replaced by these two
checks, a departure from the specification). fgw_shape0_hi is not in the pre-run gate.
Reported, deciding nothing: merged and per-truth means against both controls, neighborhood_mmd per task, the DEEP/OTHER
split, the per-truth increment over the best of fgw, fgw_shape0x_series and fgw_shape0x_ov, and, under main, the
largest difference between O_hi's four expression skills and the committed fgw_shape0_hi's (E24, E28), per task and
overall (with the same-run fgw_shape0_hi's worst-seed and mean distance from its committed board score); every group's
merged mean; the prefix build's board.
"""

from __future__ import annotations

import argparse
import contextlib
import dataclasses
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO / "experiments" / "t1_hygiene"))
sys.path.insert(0, str(REPO / "experiments" / "harness"))
sys.path.insert(0, str(REPO / "src"))  # this checkout's code, not the checkout a shared venv's vec.pth names
import e47_floor_exact as e47  # noqa: E402
import n12_hygiene_size as n12  # noqa: E402
import nine_tasks as nt  # noqa: E402

from vec import evaluate as ev  # noqa: E402
from vec import harness  # noqa: E402
from vec import shape_interp as si  # noqa: E402
from vec import t2_build as tb  # noqa: E402
from vec.adoption import SEEDS  # noqa: E402
from vec.boards import GROUPS  # noqa: E402
from vec.provenance import repo_state  # noqa: E402

ROW_IDS = ("E51", "E51a")  # the row's first word: "E51", or "E51a" if the analysts number its two tries
TRIES = REPO / "experiments" / "harness" / "TRIES.md"
SO, OHI = "fgw_shape0x_series_ov", "fgw_shape0x_ov_hi"
SERIES, OV, C0, FGW = "fgw_shape0x_series", "fgw_shape0x_ov", "fgw_shape0x", "fgw"
ARMS = (SO, OHI)
LABEL = {SO: "SO", OHI: "O_hi"}
CONTROL = {SO: SERIES, OHI: OV}  # the control each arm's pass rule reads against
CONTROLS = (SERIES, OV)
REPORTED = (FGW,)
HI = "fgw_shape0_hi"
HI_REF = (HI,)  # O_hi's expression reference: run, not counted, voids O_hi only
RUN_NAMES = (FGW, SERIES, OV, HI, SO, OHI)  # with the floor, in each run_task call
BEST_OF = (FGW, SERIES, OV)
GATE_RUN = (harness.FLOOR, C0, SERIES, OV)
BUILDS = ("main", "prefix")
TEST_LIKE = ("half", "late", "new1", "new2")
GROUP_NAMES = tuple(name for name, _, _ in GROUPS["T2"])
P4_GROUPS_E49 = ("expression_change", "cell_state", "local_organisation")  # reported beside, deciding nothing
EXPR_METRICS = e47.EXPR_METRICS  # de_score, de_direction, mmd_u, variogram
POS_METRICS = ("d2_shape", "occupancy_dice", "scale_log_ratio")
NBHD = "neighborhood_mmd"
EXPR_TOL, POS_TOL = 1e-9, 1e-6
BAR, AHEAD = 0.3, 8
GATE_TOL = 0.01  # board points (rule 9 J7; nine_tasks.TOL)
HI_TIE_TOL = GATE_TOL  # the same-run fgw_shape0_hi against its committed board score, per seed and ten-seed mean
FLOOR_SCORE = 50.0
E47_DIR = REPO / "experiments" / "harness" / "results" / "e47_floor_exact"
E47_RUN_COMMIT, E47_FILES_COMMIT = "9967cf1", "d5ee9b5"
E49_DIR = REPO / "experiments" / "coords" / "results" / "e49_shape_interp"
E49_RUN_COMMIT, E49_FILES_COMMIT = "d76ccd5", "53ba6cc"
HI_DIRS = ((REPO / "experiments" / "harness" / "results" / "e24_pricing" / "tasks", "802a0dd"),
           (REPO / "experiments" / "harness" / "results" / "e28_nine_tasks" / "tasks", "9702014"))
BASES = (E49_FILES_COMMIT,)
GUARD_FLAG = "--shape-interp-test-ok"
SCRIPT = "experiments/coords/e51_so_ohi_read.py"
RUNS = REPO / "runs" / "e51"
RESULTS = REPO / "experiments" / "coords" / "results" / "e51_so_ohi"


# ------------------------------------------------------------------------------------------------ checks
def _git(*args: str) -> subprocess.CompletedProcess:
    return subprocess.run(["git", "-C", str(REPO), *args], capture_output=True, text=True,
                          env={**os.environ, "GIT_OPTIONAL_LOCKS": "0"})


def has_row(path: Path | None = None) -> bool:
    rows = n12.tries_rows(TRIES if path is None else path)
    return any(r in rows for r in ROW_IDS)


def commit_problems(commit: str) -> list[str]:
    """Why this checkout may not run or read row E51 at `commit` (empty: it may)."""
    why = []
    head = _git("rev-parse", "HEAD").stdout.strip()
    if head != commit:
        why.append(f"HEAD is {head[:12]}, not the run commit {commit[:12]}")
    for base in BASES:
        if _git("merge-base", "--is-ancestor", base, "HEAD").returncode != 0:
            why.append(f"HEAD does not contain {base}")
    if GUARD_FLAG not in (REPO / "experiments" / "builds" / "build_test.py").read_text(encoding="utf-8"):
        why.append(f"the builder lacks the test-board guard {GUARD_FLAG}")
    for d, c in ((E47_DIR, E47_FILES_COMMIT), (E49_DIR, E49_FILES_COMMIT), *HI_DIRS):
        if _git("diff", "--quiet", c, "HEAD", "--", str(d.relative_to(REPO))).returncode != 0:
            why.append(f"{d.relative_to(REPO)} differs from {c}'s")
    return why


def _own_code() -> None:
    src = (REPO / "src").resolve()
    for mod in (ev, harness, tb, si):
        if src not in Path(mod.__file__).resolve().parents:
            sys.exit(f"{mod.__name__} is loaded from {mod.__file__}, not from this checkout's src/")


def _refuse_unless_ready(commit: str) -> None:
    _own_code()
    if not has_row():
        sys.exit("E51 has no row in experiments/harness/TRIES.md: a row is committed before its run")
    if repo_state(REPO).unclean != ():
        sys.exit("the tree is not clean: row E51 runs on a committed tree")
    why = commit_problems(commit)
    if why:
        sys.exit("refused: " + "; ".join(why))


def _task(name: str):
    tasks = nt.all_tasks()
    if name not in tasks:
        sys.exit(f"unknown task {name}; known: {sorted(tasks)}")
    t = tasks[name]
    if getattr(t, "target_filter", None) is not None:
        sys.exit(f"{name}: a task with a target filter; E51 does not run on one")
    return t


def _is_test_board(board) -> bool:
    return str(board.split).startswith("test") or "test" in str(board.key).split(":")[-1]


def _sha_file(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


# ------------------------------------------------------------------------------------------------ committed figures
def e47_scores(key: str) -> tuple[dict | None, list[str]]:
    """The committed E47 per-seed board scores of the task (task key), and why they cannot serve."""
    for p in sorted((E47_DIR / "tasks").glob("T2e_*.json")):
        r = json.loads(p.read_text())
        if r["task"]["key"] == key:
            why = []
            if not r["provenance"]["commit"].startswith(E47_RUN_COMMIT):
                why.append(f"{p.name}: E47's file is not from its run at {E47_RUN_COMMIT}")
            if not r["provenance"]["tree_clean"] or r["build"] != "main" or r["protocol"] != "bake":
                why.append(f"{p.name}: E47's file is not a clean main/bake run")
            return e47.per_seed(r)[0], why
    return None, [f"{key}: no committed E47 file"]


def e49_ov_mean(name: str, path: Path | None = None) -> tuple[float | None, list[str]]:
    """E49's committed per-task figure of fgw_shape0x_ov: its ten-seed mean board minus fgw_shape0x's, under main."""
    p = E49_DIR / "reading.json" if path is None else path
    if not p.exists():
        return None, ["no committed E49 reading"]
    r = json.loads(p.read_text())
    why = []
    if r.get("void") is not False or not str(r.get("commit", "")).startswith(E49_RUN_COMMIT):
        why.append(f"E49's reading is void or not from its run at {E49_RUN_COMMIT}")
    v = (r.get("tasks", {}).get(name, {}).get(OV) or {}).get("board", {}).get("mean")
    if v is None:
        why.append(f"{name}: no committed E49 figure for {OV}")
        return None, why
    return float(v), why


def committed_hi(dirs=None) -> tuple[dict[str, dict[str, np.ndarray]], list[str]]:
    """Task key -> metric -> the committed fgw_shape0_hi expression skills per seed (E24's six tasks, E28's three),
    and 'board' -> its committed board score per seed, in seed order."""
    out: dict[str, dict[str, np.ndarray]] = {}
    why: list[str] = []
    for d, _ in (HI_DIRS if dirs is None else dirs):
        for p in sorted(Path(d).glob("T2e_*.json")):
            r = json.loads(p.read_text())
            rows = sorted(r["rows"], key=lambda x: x["seed"])
            if HI not in rows[0].get("skills", {}):
                continue
            key = r["task"]["key"]
            if [x["seed"] for x in rows] != list(SEEDS) or r["build"] != "main" or r["protocol"] != "bake" \
                    or not r["provenance"]["tree_clean"]:
                why.append(f"{p.name}: the committed {HI} is not a clean main/bake run of seeds {list(SEEDS)}")
                continue
            if key in out:
                why.append(f"{key}: {HI} committed twice")
            out[key] = {m: np.array([x["skills"][HI][m] for x in rows], float) for m in EXPR_METRICS}
            out[key]["board"] = np.array([x["scores"][HI] for x in rows], float)
    return out, why


# ------------------------------------------------------------------------------------------------ gate
def _gate_path(runs: Path, name: str) -> Path:
    return Path(runs) / "gate" / f"{name}.json"


def gate_check(name: str, key: str, run_scores: dict, old: dict | None = None, why_old: list[str] | None = None,
               ov_mean: float | None = None, why_ov: list[str] | None = None) -> dict:
    """The pre-run gate on one task. Deciding: fgw_shape0x_series within GATE_TOL of E47's per seed and in the ten-seed
    mean; fgw_shape0x_ov's ten-seed mean minus fgw_shape0x's within GATE_TOL of E49's committed figure. Reported only:
    the floor's and fgw_shape0x's distance from 50 and from E47. `old`, `ov_mean` default to the committed figures."""
    if old is None and why_old is None:
        old, why_old = e47_scores(key)
    if ov_mean is None and why_ov is None:
        ov_mean, why_ov = e49_ov_mean(name)
    why = [*(why_old or []), *(why_ov or [])]
    worst: dict = {"series_mean": None, "series_seed": None, "ov_mean": None}
    n = {a: run_scores.get(a) for a in GATE_RUN}
    if old is not None:
        o = old.get(SERIES)
        if o is None or n[SERIES] is None or np.shape(o) != np.shape(n[SERIES]) or len(n[SERIES]) != len(SEEDS):
            why.append(f"{key}: {SERIES} missing or of another length")
        else:
            o, s = np.asarray(o, float), np.asarray(n[SERIES], float)
            worst["series_mean"], worst["series_seed"] = abs(float(s.mean() - o.mean())), float(np.abs(s - o).max())
            if not worst["series_mean"] <= GATE_TOL:
                why.append(f"{key}: {SERIES}'s ten-seed mean is {worst['series_mean']:.4f} off E47's")
            if not worst["series_seed"] <= GATE_TOL:
                why.append(f"{key}: {SERIES} is {worst['series_seed']:.4f} off E47's in a seed")
    if ov_mean is not None:
        if n[OV] is None or n[C0] is None or len(n[OV]) != len(SEEDS) or len(n[C0]) != len(SEEDS):
            why.append(f"{key}: {OV} or {C0} missing or of another length")
        else:
            d = float(np.mean(np.asarray(n[OV], float) - np.asarray(n[C0], float)))
            worst["ov_mean"] = abs(d - ov_mean)
            if not worst["ov_mean"] <= GATE_TOL:
                why.append(f"{key}: {OV} minus {C0} is {worst['ov_mean']:.4f} off E49's ten-seed mean")
    reported = {}
    if n[harness.FLOOR] is not None:
        reported["floor_off_50"] = float(np.abs(np.asarray(n[harness.FLOOR], float) - FLOOR_SCORE).max())
    if old is not None and old.get(C0) is not None and n[C0] is not None and np.shape(old[C0]) == np.shape(n[C0]):
        c = np.asarray(n[C0], float) - np.asarray(old[C0], float)
        reported["c0_vs_e47"] = {"mean": abs(float(c.mean())), "seed": float(np.abs(c).max())}
    return {"task": name, "key": key, "e47": {"run_commit": E47_RUN_COMMIT, "files_commit": E47_FILES_COMMIT},
            "e49": {"run_commit": E49_RUN_COMMIT, "files_commit": E49_FILES_COMMIT}, "tolerance": GATE_TOL,
            "worst": worst, "reported_deciding_nothing": reported, "why": why, "pass": not why}


def cmd_gate(a) -> None:
    _refuse_unless_ready(a.commit)
    t = _task(a.task)
    path = _gate_path(Path(a.runs), a.task)
    if path.exists():  # resumable: a gate file of this commit is read again, not re-run
        rec = json.loads(path.read_text())
        if rec.get("commit") != a.commit:
            sys.exit(f"{path} exists from another commit")
    else:
        res = harness.run_task(t, [harness.CANDIDATES[n] for n in GATE_RUN], seeds=SEEDS, data_dir=Path(a.data),
                               build="main", log=None)
        if str(res["board"].get("split", "")).startswith("test"):
            sys.exit("E51 met a test board")
        rec = {**gate_check(a.task, t.key, e47.per_seed(res)[0]), "commit": res["provenance"]["commit"],
               "tree_clean": res["provenance"]["tree_clean"]}
        rec["pass"] = bool(rec["pass"] and rec["commit"] == a.commit and rec["tree_clean"])
        path.parent.mkdir(parents=True, exist_ok=True)
        tmp = path.with_suffix(".tmp")
        tmp.write_text(json.dumps(rec, indent=1, default=float))
        os.replace(tmp, path)
    print(f"gate {a.task} {'pass' if rec['pass'] else 'fail'}", flush=True)
    sys.exit(0 if rec["pass"] else 3)


def gate_problems(runs: Path, commit: str) -> list[str]:
    why = []
    for name in nt.all_tasks():
        p = _gate_path(runs, name)
        if not p.exists():
            why.append(f"no gate file for {name}")
            continue
        g = json.loads(p.read_text())
        if not g.get("pass") or g.get("commit") != commit:
            why.append(f"the gate did not pass on {name} at the run commit")
    return why


# ------------------------------------------------------------------------------------------------ run
@contextlib.contextmanager
def captured():
    """While active, record the board harness.score_seeds is handed."""
    seen: dict = {}
    orig = harness.score_seeds

    def score_seeds(t, board, target, ref, inputs, candidates, *a, **kw):
        seen.update(board=board)
        return orig(t, board, target, ref, inputs, candidates, *a, **kw)

    harness.score_seeds = score_seeds
    try:
        yield seen
    finally:
        harness.score_seeds = orig


def guarded(base, seen):
    """A candidate that refuses to build on a test board."""
    def build(ctx):
        if "board" in seen and _is_test_board(seen["board"]):
            raise RuntimeError("E51 met a test board")
        return base.build(ctx)

    return dataclasses.replace(base, build=build)


def candidates(seen) -> list:
    return [harness.CANDIDATES[harness.FLOOR], *(guarded(harness.CANDIDATES[n], seen) for n in RUN_NAMES)]


def one_seed(t, s: int, data: Path) -> dict:
    out: dict = {"builds": {}, "facts": {}}
    for build in BUILDS:
        with captured() as seen:
            res = harness.run_task(t, candidates(seen), seeds=[s], data_dir=data, build=build, log=None)
        if _is_test_board(seen["board"]):
            raise RuntimeError("E51 met a test board")
        out["builds"][build] = res
        out["facts"][f"board_split_{build}"] = str(seen["board"].split)
    return out


def cmd_run(a) -> None:
    _refuse_unless_ready(a.commit)
    t = _task(a.task)
    runs = Path(a.runs)
    why = gate_problems(runs, a.commit)
    if why:
        sys.exit("refused (the pre-run gate): " + "; ".join(why))
    for c in RUN_NAMES:
        if not harness.CANDIDATES[c].applies(t):
            sys.exit(f"{a.task}: {c} does not apply")
    folder = runs / "run" / a.task
    folder.mkdir(parents=True, exist_ok=True)
    for s in SEEDS:
        path = folder / f"seed{s}.json"
        if path.exists():
            continue
        res = one_seed(t, s, Path(a.data))
        rec = {"task": a.task, "key": t.key, "seed": s, "run_commit": a.commit, "script": SCRIPT,
               "names": list(RUN_NAMES), **res}
        tmp = path.with_suffix(".tmp")
        tmp.write_text(json.dumps(rec, indent=1, default=float))
        os.replace(tmp, path)
        print(f"{a.task} seed {s} done", flush=True)
    print(f"wrote {folder}", flush=True)


# ------------------------------------------------------------------------------------------------ read
def _row(r: dict, build: str = "main") -> dict:
    return r["builds"][build]["rows"][0]


def _has(row: dict, n: str) -> bool:
    return all(n in row.get(k, {}) for k in ("scores", "skills", "groups"))


def _hi_boards(rs: list[dict], ref: dict | None) -> tuple[np.ndarray, np.ndarray] | None:
    """Under main, the same-run fgw_shape0_hi's board per seed and its committed one, or None where either is
    missing (a missing run arm is voided elsewhere)."""
    if ref is None or "board" not in ref:
        return None
    rows = [r["builds"].get("main", {}).get("rows", [{}])[0] for r in rs]
    if not all(HI in x.get("scores", {}) for x in rows):
        return None
    s, o = np.array([x["scores"][HI] for x in rows], float), np.asarray(ref["board"], float)
    return (s, o) if s.shape == o.shape else None


def hi_tie_problems(task: str, key: str, rs: list[dict], hi: dict) -> list[str]:
    """O_hi's link to the committed record: under main, the same-run fgw_shape0_hi's board within HI_TIE_TOL of its
    committed score (E24, E28) in every seed and in the ten-seed mean (rule 9 J7, nine_tasks.TOL). Voids O_hi only."""
    ref = hi.get(key)
    if ref is None or "board" not in ref:
        return [f"{task}: no committed {HI}"]
    if not all(_has(r["builds"].get("main", {}).get("rows", [{}])[0], HI) for r in rs):
        return []  # the missing arm already voids O_hi
    pair = _hi_boards(rs, ref)
    if pair is None:
        return [f"{task}: the committed {HI} is of another length"]
    s, o = pair
    why = []
    d_seed, d_mean = float(np.abs(s - o).max()), abs(float(s.mean() - o.mean()))
    if not d_seed <= HI_TIE_TOL:
        why.append(f"{task}: the same-run {HI} is {d_seed:.4f} off its committed score in a seed")
    if not d_mean <= HI_TIE_TOL:
        why.append(f"{task}: the same-run {HI} is {d_mean:.4f} off its committed score in the ten-seed mean")
    return why


def check_run(recs: dict[str, list[dict]], why_in: list[str], commit: str,
              hi: dict[str, dict[str, np.ndarray]], why_hi: list[str]) -> tuple[list[str], dict[str, list[str]]]:
    """The void conditions: (reasons that void both arms, reasons per arm)."""
    why = list(why_in)
    arm_why: dict[str, list[str]] = {a: [] for a in ARMS}
    arm_why[OHI].extend(why_hi)
    tasks = nt.all_tasks()
    for task in tasks:
        rs = recs.get(task, [])
        if [r["seed"] for r in rs] != list(SEEDS):
            why.append(f"{task}: seeds are not {list(SEEDS)}")
            continue
        key = tasks[task].key
        for r in rs:
            tag = f"{task} seed {r['seed']}"
            if r.get("run_commit") != commit:
                why.append(f"{tag}: another run commit")
            for b in BUILDS:
                res = r["builds"].get(b)
                if res is None:
                    why.append(f"{tag}: build {b} missing")
                    continue
                p = res["provenance"]
                if p["commit"] != commit or not p["tree_clean"]:
                    why.append(f"{tag}: build {b} at another commit or on an unclean tree")
                if res["build"] != b or res["protocol"] != "bake" or [x["seed"] for x in res["rows"]] != [r["seed"]]:
                    why.append(f"{tag}: build {b} is not the bake protocol at this seed and build")
                if str(res["board"].get("split", "")).startswith("test") \
                        or str(r["facts"].get(f"board_split_{b}", "")).startswith("test"):
                    why.append(f"{tag}: a test board")
                row = res["rows"][0]
                missing = [n for n in (harness.FLOOR, *CONTROLS, *REPORTED) if not _has(row, n)]
                if missing:
                    why.append(f"{tag}: build {b} lacks {', '.join(missing)}")
                    continue
                sk = row["skills"]
                if not _has(row, SO):
                    arm_why[SO].append(f"{tag}: build {b} lacks {SO}")
                else:
                    if not max(abs(sk[SO][m] - sk[SERIES][m]) for m in EXPR_METRICS) <= EXPR_TOL:
                        arm_why[SO].append(f"{tag}: {SO}'s expression skills differ from {SERIES}'s ({b})")
                    if not max(abs(sk[SO][m] - sk[OV][m]) for m in POS_METRICS) <= POS_TOL:
                        arm_why[SO].append(f"{tag}: {SO}'s position skills differ from {OV}'s ({b})")
                if not _has(row, HI):  # O_hi's expression reference: its absence voids O_hi only
                    arm_why[OHI].append(f"{tag}: build {b} lacks {HI}")
                if not _has(row, OHI):
                    arm_why[OHI].append(f"{tag}: build {b} lacks {OHI}")
                    continue
                if not max(abs(sk[OHI][m] - sk[OV][m]) for m in POS_METRICS) <= POS_TOL:
                    arm_why[OHI].append(f"{tag}: {OHI}'s position skills differ from {OV}'s ({b})")
                if _has(row, HI) and not max(abs(sk[OHI][m] - sk[HI][m]) for m in EXPR_METRICS) <= EXPR_TOL:
                    arm_why[OHI].append(f"{tag}: {OHI}'s expression skills differ from the same-run {HI}'s ({b})")
        arm_why[OHI].extend(hi_tie_problems(task, key, rs, hi))
    arm_why = {a: sorted(set(v), key=v.index) for a, v in arm_why.items()}
    return why, arm_why


def series(recs: list[dict], name: str, what: str, build: str = "main") -> np.ndarray:
    """Per seed, one quantity of one file: 'board', 'group:<g>' (gated, as the row holds it) or 'skill:<m>'."""
    rows = [_row(r, build) for r in recs]
    if what == "board":
        return np.array([x["scores"][name] for x in rows], float)
    if what.startswith("group:"):
        return np.array([x["groups"][name][what[6:]] for x in rows], float)
    if what.startswith("skill:"):
        return np.array([x["skills"][name][what[6:]] for x in rows], float)
    raise ValueError(what)


def merged_D(recs, groups: dict[str, list[str]], a: str, b: str, what: str, build: str = "main") -> np.ndarray:
    """D_s: the mean over truths of each truth's task-mean difference (a minus b) at each seed (J4)."""
    return e47.merged_seeds({n: series(recs[n], a, what, build) - series(recs[n], b, what, build)
                             for names in groups.values() for n in names}, groups)


def _stat(D: np.ndarray) -> dict:
    return {"mean": float(D.mean()), "se": float(D.std(ddof=1) / np.sqrt(len(D))), "ahead": int((D > 0).sum()),
            "behind": int((D < 0).sum())}


def behind(D: np.ndarray) -> bool:
    """E49's P4 definition: a merged mean below 0 with at least AHEAD of the ten D_s below 0."""
    return bool(D.mean() < 0 and (D < 0).sum() >= AHEAD)


def per_truth(recs, groups, name: str, build: str = "main") -> dict[str, float]:
    """Per truth, the mean board of one file (mean over the truth's tasks of the ten-seed mean)."""
    return {truth: float(np.mean([series(recs[n], name, "board", build).mean() for n in names]))
            for truth, names in groups.items()}


def increment_over_best(recs, groups, arm: str) -> dict:
    """Per truth, the arm's mean board minus the best of fgw, fgw_shape0x_series and fgw_shape0x_ov on that truth."""
    mine = per_truth(recs, groups, arm)
    others = {n: per_truth(recs, groups, n) for n in BEST_OF}
    out = {}
    for truth in groups:
        best = max(BEST_OF, key=lambda n: others[n][truth])
        out[truth] = {"increment": mine[truth] - others[best][truth], "best": best}
    return out


def expr_vs_committed_hi(recs, hi: dict) -> dict:
    """Reported, deciding nothing: under main, per task and overall, the largest |O_hi's expression skill minus the
    committed fgw_shape0_hi's| over seeds and the four expression metrics."""
    tasks = nt.all_tasks()
    per = {}
    for t, rs in recs.items():
        ref = hi.get(tasks[t].key)
        if ref is None:
            continue
        per[t] = float(max(np.abs(series(rs, OHI, f"skill:{m}") - np.asarray(ref[m], float)).max()
                           for m in EXPR_METRICS))
    return {"per_task": per, "overall": max(per.values()) if per else None}


def hi_board_vs_committed(recs, hi: dict) -> dict:
    """Per task, the same-run fgw_shape0_hi's worst-seed and ten-seed-mean distance from its committed board score
    under main (the deciding values of the tie, shown)."""
    tasks = nt.all_tasks()
    out = {}
    for t, rs in recs.items():
        pair = _hi_boards(rs, hi.get(tasks[t].key))
        if pair is not None:
            s, o = pair
            out[t] = {"seed": float(np.abs(s - o).max()), "mean": abs(float(s.mean() - o.mean()))}
    return out


def arm_reading(recs, groups, tl_groups, arm: str, hi: dict | None = None) -> dict:
    ctl = CONTROL[arm]
    full = merged_D(recs, groups, arm, ctl, "board")
    tl = merged_D(recs, tl_groups, arm, ctl, "board")
    grp = {g: merged_D(recs, groups, arm, ctl, f"group:{g}") for g in GROUP_NAMES}
    if arm == SO:
        clauses = {"a": bool(tl.mean() > BAR and (tl > 0).sum() >= AHEAD),
                   "b": bool(full.mean() > 0 and (full > 0).sum() >= AHEAD),
                   "c": not any(behind(D) for D in grp.values())}
    else:
        clauses = {"a": bool(full.mean() > BAR and (full > 0).sum() >= AHEAD),
                   "b": bool(tl.mean() > 0),
                   "c": not any(behind(D) for D in grp.values())}
    against = {}
    for c in CONTROLS:
        against[c] = {"merged_board": _stat(merged_D(recs, groups, arm, c, "board")),
                      "test_like_board": _stat(merged_D(recs, tl_groups, arm, c, "board")),
                      "prefix_board": _stat(merged_D(recs, groups, arm, c, "board", build="prefix")),
                      "per_truth": {t: v - per_truth(recs, groups, c)[t]
                                    for t, v in per_truth(recs, groups, arm).items()}}
    per_task = {}
    for t, rs in recs.items():
        per_task[t] = {"board_minus_control": float(np.mean(series(rs, arm, "board") - series(rs, ctl, "board"))),
                       NBHD: float(series(rs, arm, f"skill:{NBHD}").mean()),
                       **{f"{NBHD}_minus_{c}": float(np.mean(series(rs, arm, f"skill:{NBHD}")
                                                             - series(rs, c, f"skill:{NBHD}"))) for c in CONTROLS}}
    split = {k: float(np.mean([per_task[t]["board_minus_control"] for t in ks]))
             for k, ks in (("DEEP", nt.DEEP), ("OTHER", nt.OTHER))}
    return {"label": LABEL[arm], "control": ctl, "clauses": clauses, "passes": all(clauses.values()),
            "merged_board": _stat(full), "test_like_board": _stat(tl),
            "groups": {g: {**_stat(D), "is_behind": behind(D)} for g, D in grp.items()},
            "reported_deciding_nothing": {
                "behind_on_e49_p4_groups": [g for g in P4_GROUPS_E49 if behind(grp[g])],
                "against": against, "per_task": per_task, "deep_other": split,
                "increment_over_best": increment_over_best(recs, groups, arm),
                **({"expr_vs_committed_hi": expr_vs_committed_hi(recs, hi or {}),
                    "hi_board_vs_committed": hi_board_vs_committed(recs, hi or {})} if arm == OHI else {})}}


def outcome(arms: dict) -> str:
    words = []
    for a in ARMS:
        r = arms[a]
        if r.get("void"):
            words.append(f"{LABEL[a]}: void, read as nothing")
        elif r["passes"]:
            words.append(f"{LABEL[a]}: passes; the question goes to the captain")
        else:
            words.append(f"{LABEL[a]}: does not pass; no file")
    passing = [a for a in ARMS if not arms[a].get("void") and arms[a]["passes"]]
    if len(passing) == 2:
        words.append("both go to the captain; the read does not pick between them")
    return "; ".join(words)


def reading(recs: dict[str, list[dict]], why_in: list[str], commit: str,
            hi: dict[str, dict[str, np.ndarray]] | None = None, why_hi: list[str] | None = None) -> dict:
    if hi is None and why_hi is None:
        hi, why_hi = committed_hi()
    why, arm_why = check_run(recs, why_in, commit, hi or {}, list(why_hi or []))
    if why:
        return {"void": True, "why": why, "arm_why": arm_why, "commit": commit}
    tasks = nt.all_tasks()
    groups = e47.truths(tasks)
    tl_groups = e47.truths({n: tasks[n] for n in TEST_LIKE})
    arms = {a: ({"void": True, "why": arm_why[a], "label": LABEL[a]} if arm_why[a]
                else {"void": False, **arm_reading(recs, groups, tl_groups, a, hi)}) for a in ARMS}
    return {"void": False, "why": [], "commit": commit, "truths": groups, "test_like_truths": tl_groups,
            "arms": arms, "outcome": outcome(arms),
            "gate": {"e47_run_commit": E47_RUN_COMMIT, "e47_files_commit": E47_FILES_COMMIT,
                     "e49_run_commit": E49_RUN_COMMIT, "e49_files_commit": E49_FILES_COMMIT, "tolerance": GATE_TOL}}


def text(res: dict) -> str:
    if res["void"]:
        return f"Row E51 (run commit {res['commit']}): VOID: " + "; ".join(res["why"]) + "\n"
    L = [f"Row E51 (run commit {res['commit']}); the pre-run gate passed on all nine tasks.",
         "Readings fixed in row E51, build main (mean, SE, seeds ahead of 10):"]
    for a in ARMS:
        r = res["arms"][a]
        if r["void"]:
            L.append(f"  {LABEL[a]} {a}: VOID: " + "; ".join(r["why"]))
            continue
        L.append(f"  {r['label']} {a}: {'PASSES' if r['passes'] else 'does not pass'}; clauses "
                 + ", ".join(f"({k}) {'holds' if v else 'fails'}" for k, v in r["clauses"].items()))
        for k in ("merged_board", "test_like_board"):
            v = r[k]
            L.append(f"    {k} against {r['control']}: {v['mean']:+.3f} (SE {v['se']:.3f}, {v['ahead']} ahead)")
        L.append("    groups: " + "; ".join(f"{g} {v['mean']:+.4f} ({v['behind']} seeds behind"
                                            f"{'; BEHIND' if v['is_behind'] else ''})" for g, v in r["groups"].items()))
        rep = r["reported_deciding_nothing"]
        L.append("    reported, deciding nothing:")
        for c, v in rep["against"].items():
            L.append(f"      against {c}: merged {v['merged_board']['mean']:+.3f}, test-like "
                     f"{v['test_like_board']['mean']:+.3f}, prefix {v['prefix_board']['mean']:+.3f}; per truth "
                     + ", ".join(f"{t} {x:+.3f}" for t, x in v["per_truth"].items()))
        L.append("      increment over the best of fgw, fgw_shape0x_series, fgw_shape0x_ov per truth: "
                 + ", ".join(f"{t} {v['increment']:+.3f} ({v['best']})" for t, v in rep["increment_over_best"].items()))
        L.append(f"      DEEP/OTHER board against {r['control']}: {rep['deep_other']['DEEP']:+.3f}/"
                 f"{rep['deep_other']['OTHER']:+.3f}")
        L.append("      neighborhood_mmd skill per task: " + ", ".join(f"{t} {v[NBHD]:.3f}"
                                                                     for t, v in rep["per_task"].items()))
        L.append("      behind on E49's three P4 groups: " + (", ".join(rep["behind_on_e49_p4_groups"]) or "none"))
        if "expr_vs_committed_hi" in rep:
            e = rep["expr_vs_committed_hi"]
            L.append(f"      largest |expression skill - committed {HI}'s| (main): overall "
                     + ("none" if e["overall"] is None else f"{e['overall']:.2e}") + "; "
                     + ", ".join(f"{t} {v:.2e}" for t, v in e["per_task"].items()))
            L.append(f"      same-run {HI} board against its committed score (main), worst seed/mean: "
                     + ", ".join(f"{t} {v['seed']:.4f}/{v['mean']:.4f}"
                                 for t, v in rep["hi_board_vs_committed"].items()))
    L.append(f"Outcome: {res['outcome']}.")
    L.append("A fail means no file; a pass only puts the question to the captain. Nothing here shows E7.75's depth or "
             "profile, the overlay's geometry at the test pair, the live gate, scorer build or metric seed.")
    return "\n".join(L) + "\n"


def load_runs(runs: Path) -> dict[str, list[dict]]:
    out = {}
    for name in nt.all_tasks():
        folder = Path(runs) / "run" / name
        recs = [json.loads(p.read_text()) for p in sorted(folder.glob("seed*.json"))] if folder.exists() else []
        out[name] = sorted(recs, key=lambda r: r["seed"])
    return out


def provenance() -> dict:
    s = repo_state(REPO)
    return {"head": s.commit, "clean": s.unclean == (), "script_sha256": _sha_file(Path(__file__)),
            "harness_sha256": _sha_file(Path(harness.__file__)), "t2_build_sha256": _sha_file(Path(tb.__file__)),
            "shape_interp_sha256": _sha_file(Path(si.__file__)),
            "e49_reading_sha256": _sha_file(E49_DIR / "reading.json")}


def cmd_read(a) -> None:
    out = Path(a.out)
    files = [out / "reading.json", out / "reading.txt"]
    existing = [str(f) for f in files if f.exists()]
    if existing:
        sys.exit(f"refusing to overwrite: {', '.join(existing)}")
    _own_code()
    if not has_row():
        sys.exit("E51 has no row in experiments/harness/TRIES.md")
    runs = Path(a.runs)
    why = [*gate_problems(runs, a.commit), *commit_problems(a.commit)]
    if repo_state(REPO).unclean != ():
        why.append("the tree is not clean at the read")
    res = {**reading(load_runs(runs), why, a.commit), "reading_code": provenance()}
    out.mkdir(parents=True, exist_ok=True)
    files[0].write_text(json.dumps(res, indent=1, default=float))
    files[1].write_text(text(res))
    print(f"wrote {files[1]}", flush=True)


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("gate", "run"):
        p = sub.add_parser(name)
        p.add_argument("--task", required=True, help="a short name of nine_tasks.all_tasks()")
        p.add_argument("--commit", required=True, help="the run commit (full hash)")
        p.add_argument("--runs", default=str(RUNS))
        p.add_argument("--data", default=str(harness.DATA_DIR))
    r = sub.add_parser("read")
    r.add_argument("--commit", required=True)
    r.add_argument("--runs", default=str(RUNS))
    r.add_argument("--out", default=str(RESULTS))
    a = ap.parse_args(argv)
    {"gate": cmd_gate, "run": cmd_run, "read": cmd_read}[a.cmd](a)


if __name__ == "__main__":
    main()
