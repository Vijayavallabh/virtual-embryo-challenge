#!/usr/bin/env python
"""Row E58 of experiments/harness/TRIES.md: gate G1 of row R59's file T2embryo_s2x_ovt23.h5ad (Ot23, the captain's
named exception on the embryo's overlay form, bound to this gate), a counted practice read as the cold challenge of the
step-change round rewrites the design's G1 (EMBRYO_PORTFOLIO_challenge.md, s. 3 items 3-4 and s. 4). Its structure is
copied from experiments/coords/e55_overlay_stretch_read.py and changed where the challenge says. No number here comes
from a test score. It moves the overlay's form: the main session is recused on the embryo's shape forms and occupancy
term, and acts on the outcome line alone.

Arm, counted as one try when the run step starts: ARM = fgw_shape0x_ov_t23 = harness._fgw(ctx, t_pos=0.0,
positions='floor_exact', shape_form='overlay', t_mix=min(1, t + 1/6)), t the task's own bracket fraction: O's command
with its mix leaned a sixth toward the later bracket (at the test pair, t = 1/2, that is the builder's --mix-lean later,
t_mix = 2/3). Control, rescored and not counted: O = fgw_shape0x_ov (E49c). The arm is defined here, not in vec.harness.

Tasks (J4: tasks on one truth are merged into one first). Counted: E8.0 = half and late (one truth), E8.25_late = new1;
two truths at equal weight. Reported, deciding nothing: new2, x1, x2, x3. Not run at all: e75 and e725, barred for any
change to the mixing fraction (the ledger's prune list: "E7.25 and E7.5 cannot count for any change to the mixing
fraction"). Seeds 0-9, bake protocol, build seed BUILD_SEED_OFFSET + seed, scorer build main only.

The pass rule (PASS_RULE, fixed before any output): against O under build main on the two-truth merge, D_s the mean over
the two truths of each truth's task-mean difference at seed s; G1 passes only if P1 the merged board gain is at least
+0.3 with D_s above 0 in at least 8 of 10 seeds, P2 no group (expression_change, cell_state, local_organisation) is
behind on that merge by E49's definition (merged mean below 0 with at least 8 of 10 D_s below 0), and P3 the arm refused
on no task-seed. A merged gain below 0.3 in size is the idle case: STOP. shape_growth's merged change is reported beside
it for the captain (recused) and decides nothing. Outcomes, in e58_row's fixed words: PASS_TEXT, STOP_TEXT, IDLE_TEXT,
VOID_TEXT. A pass builds row R59's file (experiments/builds/embryo_ot23.sh, which reads row E58's outcome cell by code).

Steps (each refuses unless TRIES.md holds exactly one well-formed row E58 (e58_row), the tree is clean and HEAD is the
run commit):
  gate --task NAME --commit SHA   the floor and O under main, ten seeds (J7 guard): the floor against E47's committed
                                  scores within 0.01 board points per seed and per ten-seed mean, and O's ten-seed mean
                                  minus the floor's within 0.01 of E55's committed per-task figure (reading.json
                                  tasks[name].board_minus_floor.fgw_shape0x_ov, sha256 E55_READING_SHA256); exit 3 on a
                                  failure; writes <runs>/gate/<task>.json
  run --task NAME --commit SHA    the floor, O and the arm under main, seed by seed (resumable); refused unless all
                                  seven gate files pass at this commit
  read --commit SHA               the reading into experiments/coords/results/e58_mix_lean/ (reading.json, reading.txt),
                                  never overwriting
  j7 [--out FILE]                 the J7 dry run of the O guard on committed records only (not a try, no run): E55's
                                  committed O-minus-floor per task against E49's committed O-minus-C0 plus E47's
                                  committed C0-minus-floor, the same 0.01 tolerance; exit 0 pass, 3 fail
The script prints only "<task> seed <s> done", "wrote <path>", "gate <task> pass|fail" and "j7 pass|fail". It reads no
test board, test file or test anchor, and prints no count.
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
sys.path.insert(0, str(REPO / "experiments" / "harness"))
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(REPO / "src"))  # this checkout's code, not the checkout a shared venv's vec.pth names
import e47_floor_exact as e47  # noqa: E402
import e58_row as er  # noqa: E402
import nine_tasks as nt  # noqa: E402

from vec import harness  # noqa: E402
from vec import shape_interp as si  # noqa: E402
from vec import t2_build as tb  # noqa: E402
from vec.adoption import SEEDS  # noqa: E402
from vec.provenance import repo_state  # noqa: E402

TRIES = REPO / "experiments" / "harness" / "TRIES.md"
O = "fgw_shape0x_ov"  # noqa: E741  (the specification names the control O)
ARM = "fgw_shape0x_ov_t23"
TRY_ID = "E58"
MIX_LEAN_STEP = 1.0 / 6.0  # the lean: t_mix = min(1, t + 1/6); at t = 1/2 the same double as 2/3
COUNTED_TASKS = ("half", "late", "new1")
REPORTED_TASKS = ("new2", "x1", "x2", "x3")
BARRED_TASKS = ("e75", "e725")  # never run: barred for any change to the mixing fraction
TASKS = (*COUNTED_TASKS, *REPORTED_TASKS)
BUILD = "main"
OTHER_GROUPS = ("expression_change", "cell_state", "local_organisation")
GROUP_SHAPE = "shape_growth"
SHAPE_TERMS = ("d2_shape", "occupancy_dice")
BAR, AHEAD = 0.3, 8
T_MIX_TOL = 1e-12
GATE_TOL = 0.01  # board points: the floor against E47 per seed and per ten-seed mean (nine_tasks.TOL)
O_TOL = 0.01  # board points: O's ten-seed mean minus the floor's against E55's committed per-task figure
E47_DIR = REPO / "experiments" / "harness" / "results" / "e47_floor_exact"
E47_RUN_COMMIT = "9967cf1"
E47_SHA256 = {  # the committed E47 task files the gate reads (sha256 of each, at this script's writing)
    "T2e_E7.5_E8.5_E8.0.json": "e6ab844a3f8e8a2fc1d296a25ba75768f3c42c90a34ce9ac9d6c633442d36bc1",
    "T2e_E7.5_E8.25_late_E8.0.json": "59c345f2a64798c1362fe86b545efafed6ab6d8fa4d55a7254967fa7b3b0e041",
    "T2e_E7.5_E8.5_E8.25_late.json": "02c37ee2b5f987564429dde3ea98cde22a05b2ffbe240cbd15239399fac5013d",
    "T2e_E7.5_E8.75_E8.5.json": "878c5b5380d4e4046c94a711366bd258c0aea93b5eb06d868d3151eeb307a068",
    "T2e_E8.0_E8.75_E8.25_late.json": "03b0483396cd59cc6d4a013d9294a70ef06a9ec4ee27ba715024a5c8321aafb7",
    "T2e_E8.0_E8.75_E8.5.json": "c31044dcad6e19dd864c0032db415a7103dd27dd1f26901c2942ea92d9a154e5",
    "T2e_E8.25_late_E8.75_E8.5.json": "e40837de771777920abcb742ec22122f8345fd0e9de9d22c91b650544595d5d4"}
E55_READING = REPO / "experiments" / "coords" / "results" / "e55_overlay_stretch" / "reading.json"
E55_READING_SHA256 = "e500c9fa4161720191562598af76f8fb28bc12b0c3b0dfa9fdd3e38a38153ab5"  # as row E55's cell 11
E55_RUN_COMMIT = "4535f98c98ae1a2a8488185e8a99e4a2e8d92d86"
E49_READING = REPO / "experiments" / "coords" / "results" / "e49_shape_interp" / "reading.json"
E49_READING_SHA256 = "f4630e4583563d7c46dbb9ed39a97cc38a94d0f876dc31109cae8e7011e17115"
E49_RUN_COMMIT = "d76ccd5"
C0 = "fgw_shape0x"  # E47's and E49's reference arm, read from their committed records by the j7 dry run only
GATE_RUN = (harness.FLOOR, O)
SCRIPT = "experiments/coords/e58_mix_lean_read.py"
RUNS = REPO / "runs" / "e58"
RESULTS = REPO / "experiments" / "coords" / "results" / "e58_mix_lean"
REFUSAL = "cloud guard flags"
CANDIDATE_COMMAND = ("experiments/builds/build_test.py embryo --method fgw_shape0 --positions floor_exact "
                     "--expression mixture --shape-interp overlay --shape-interp-test-ok --mix-lean later "
                     "--captain-ruling R59 --out <run folder>/T2embryo_s2x_ovt23.h5ad "
                     "(experiments/builds/embryo_ot23.sh)")

CLAUSES = {
    "P1": "P1 the two-truth merged board gain against O is at least +0.3, with D_s above 0 in at least 8 of 10 seeds.",
    "P2": "P2 no group (expression_change, cell_state, local_organisation) is behind O on that merge by E49's "
          "definition (merged mean below 0 with at least 8 of 10 D_s below 0), with no tolerance.",
    "P3": "P3 the arm refused on no task-seed (a refusal makes it unreadable: STOP, not void).",
}
PASS_RULE = (
    "All under build main, against O (fgw_shape0x_ov), on the two-truth merge of the counted tasks (E8.0 = half and "
    "late, one truth; E8.25_late = new1; equal weight, J4), D_s the mean over the two truths of each truth's task-mean "
    "difference at seed s. " + " ".join(CLAUSES.values()) + " A merged gain below 0.3 in size is the idle case: STOP. "
    "shape_growth's merged change is reported beside it for the captain (recused) and decides nothing; new2, x1, x2 "
    "and x3 are reported and decide nothing; e75 and e725 are not run (barred for any change to the mixing fraction).")
VOID_CONDITIONS = (
    "a missing task, seed or arm (the floor, O and the arm in every row)",
    "a result from another run commit or an unclean tree, or not the bake protocol under build main at its seed",
    "a test board",
    "a gate file missing or failing at the run commit (the J7 guard)",
    "the arm's count off O's (the count rule is O's)",
    "the arm's recorded t_mix off min(1, t + 1/6), or O's off t (the file's own path)",
    "the arm's size target off O's (the size rule is O's)",
)
PRICE_LINES = (
    "J1: the run applies the file's own build path, O's command with the drawn cells' mix leaned a sixth toward the "
    "later bracket (t_mix = 2/3 at the test pair's t = 1/2), except the size step, which is the harness's task size "
    "rule (the file takes the pre-registered size; the shape terms ignore scale). The test's world (which run imaged "
    "E7.75) is not hosted by any counted task: a gain here prices only the world where the truth shares the later "
    "bracket's run, and nothing here prices the other worlds.",
    "J4: the counted merge weighs its two truths equally (half and late are one truth, E8.0); a world weight, where "
    "the captain's question uses one, takes the least favourable value of N23's range for the claim argued (0.4 for "
    "the later bracket's world).",
    "J5: diagnostic; the noise floor 6.54 exceeds twice any predicted gain.",
)


# ------------------------------------------------------------------------------------------------ the arm
def frac_of(t) -> float:
    """The task's bracket fraction t, as harness._fgw computes it from the task's view."""
    v = t if isinstance(t, harness.TaskView) else harness.TaskView.of(t)
    lo, hi = v.role_days["before"], v.role_days["after"]
    return (v.target_day - lo) / (hi - lo)


def arm_t_mix(frac: float) -> float:
    """The arm's later-bracket share: the task's own fraction leaned a sixth toward the later bracket, at most 1."""
    return min(1.0, float(frac) + MIX_LEAN_STEP)


def _arm_build(ctx: harness.BuildContext):
    return harness._fgw(ctx, t_pos=0.0, positions=tb.POSITIONS[1], shape_form="overlay",
                        t_mix=arm_t_mix(frac_of(ctx.task)))


ARM_CANDIDATE = harness.Candidate(
    ARM, frozenset({"T2"}), (), _arm_build,
    {"what": "fgw_shape0x_ov (the overlay form on the earlier bracket's own positions) with the drawn cells' mix "
             "leaned a sixth toward the later bracket, t_mix = min(1, t + 1/6) (row E58, gate G1 of row R59's file)",
     "positions": tb.POSITIONS[1], "stream": tb.FLOOR_EXACT_STREAM, "shape_form": "overlay",
     "t_mix": "min(1, t + 1/6)", "scope": "whole"}, requires=harness._interp)


# ------------------------------------------------------------------------------------------------ checks
def _git(*args: str) -> subprocess.CompletedProcess:
    return subprocess.run(["git", "-C", str(REPO), *args], capture_output=True, text=True,
                          env={**os.environ, "GIT_OPTIONAL_LOCKS": "0"})


def _sha_file(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def commit_problems(commit: str) -> list[str]:
    """Why this checkout may not run or read row E58 at `commit` (empty: it may)."""
    why = []
    head = _git("rev-parse", "HEAD").stdout.strip()
    if head != commit:
        why.append(f"HEAD is {head[:12]}, not the run commit {commit[:12]}")
    for name, sha in E47_SHA256.items():
        p = E47_DIR / "tasks" / name
        if not p.exists() or _sha_file(p) != sha:
            why.append(f"E47's committed {name} is missing or changed")
    if not E55_READING.exists() or _sha_file(E55_READING) != E55_READING_SHA256:
        why.append("E55's committed reading is missing or changed")
    return why


def _own_code() -> None:
    src = (REPO / "src").resolve()
    for mod in (harness, tb, si):
        if src not in Path(mod.__file__).resolve().parents:
            sys.exit(f"{mod.__name__} is loaded from {mod.__file__}, not from this checkout's src/")


def _refuse_unless_ready(commit: str) -> None:
    _own_code()
    why = er.row_problems(TRIES)
    if why:
        sys.exit("E58 has no well-formed row in experiments/harness/TRIES.md (a row is committed before its run): "
                 + "; ".join(why))
    if repo_state(REPO).unclean != ():
        sys.exit("the tree is not clean: row E58 runs on a committed tree")
    why = commit_problems(commit)
    if why:
        sys.exit("refused: " + "; ".join(why))


def _task(name: str):
    if name in BARRED_TASKS:
        sys.exit(f"{name} is barred for any change to the mixing fraction: row E58 does not run it")
    tasks = nt.all_tasks()
    if name not in TASKS or name not in tasks:
        sys.exit(f"unknown task {name}; row E58 runs {list(TASKS)}")
    t = tasks[name]
    if getattr(t, "target_filter", None) is not None:
        sys.exit(f"{name}: a task with a target filter; E58 does not run on one")
    return t


def _is_test_board(board) -> bool:
    return str(board.split).startswith("test") or "test" in str(board.key).split(":")[-1]


# ------------------------------------------------------------------------------------------------ the run's wrappers
@contextlib.contextmanager
def captured():
    """While active, record what harness.score_seeds is handed (the board) and each interpolate_clouds record by
    (who, build seed), who "arm" when the call passes a t_mix and "O" when it passes none: flat values only."""
    seen: dict = {"stage": {}, "refused": {}, "records": {}}
    orig_ss, orig_ic = harness.score_seeds, tb.interpolate_clouds

    def score_seeds(t, board, target, ref, inputs, candidates, *a, **kw):
        seen.update(board=board)
        return orig_ss(t, board, target, ref, inputs, candidates, *a, **kw)

    def interpolate_clouds(*a, **kw):
        res = orig_ic(*a, **kw)
        seed = int(kw["seed"] if "seed" in kw else a[6])
        who = O if kw.get("t_mix") is None else ARM
        r = res.record
        seen["records"][(who, seed)] = {"t_mix": float(r["t_mix"]), "t": float(r["t"]),
                                        "n_jittered": int(r.get("n_jittered", 0)),
                                        "target_rms": float(r["target_rms"])}
        return res

    harness.score_seeds, tb.interpolate_clouds = score_seeds, interpolate_clouds
    try:
        yield seen
    finally:
        harness.score_seeds, tb.interpolate_clouds = orig_ss, orig_ic


def observed(base, seen):
    """O unchanged; its file is kept as the stand-in for a refused arm."""
    def build(ctx):
        if _is_test_board(seen["board"]):
            raise RuntimeError("E58 met a test board")
        pred = base.build(ctx)
        seen["stage"][(base.name, ctx.seed - harness.BUILD_SEED_OFFSET)] = pred
        return pred

    return dataclasses.replace(base, build=build)


def guarded(base, seen):
    """The arm: its file, or on a refusal (ShapeRefused, the cloud guard) O's file of the seed."""
    def build(ctx):
        s = ctx.seed - harness.BUILD_SEED_OFFSET
        if _is_test_board(seen["board"]):
            raise RuntimeError("E58 met a test board")
        try:
            pred, why = base.build(ctx), None
        except si.ShapeRefused as e:
            pred, why = None, f"shape: {e}"
        except ValueError as e:
            if REFUSAL not in str(e):
                raise
            pred, why = None, "cloud guard"
        if pred is None:
            pred = seen["stage"][(O, s)]
        seen["refused"][(base.name, s)] = why
        seen["stage"][(base.name, s)] = pred
        return pred

    return dataclasses.replace(base, build=build)


def candidates(seen) -> list:
    return [harness.CANDIDATES[harness.FLOOR], observed(harness.CANDIDATES[O], seen), guarded(ARM_CANDIDATE, seen)]


def one_seed(t, s: int, data: Path) -> dict:
    """One task-seed under build main: run_task's result and E58's per-seed facts."""
    with captured() as seen:
        res = harness.run_task(t, candidates(seen), seeds=[s], data_dir=data, build=BUILD, log=None)
    board = seen["board"]
    if _is_test_board(board):
        raise RuntimeError("E58 met a test board")
    b = harness.BUILD_SEED_OFFSET + s
    facts = {"board_split": str(board.split), "frac": frac_of(t),
             "n": {n: int(len(seen["stage"][(n, s)].X)) for n in (O, ARM)},
             "records": {n: seen["records"].get((n, b), {}) for n in (O, ARM)},
             "refused": {ARM: seen["refused"].get((ARM, s))}}
    return {"builds": {BUILD: res}, "facts": facts}


# ------------------------------------------------------------------------------------------------ the J7 guard
def e55_o_per_task(path: Path | None = None, sha: str | None = E55_READING_SHA256) -> tuple[dict | None, list[str]]:
    """E55's committed per-task figure of O: its ten-seed mean board minus the floor's under main (reading.json
    tasks[name].board_minus_floor.fgw_shape0x_ov), from the reading of sha256 `sha` made at E55's run commit."""
    p = E55_READING if path is None else Path(path)
    if not p.exists():
        return None, ["no committed E55 reading"]
    if sha is not None and _sha_file(p) != sha:
        return None, [f"E55's reading is not the committed one of sha256 {sha[:12]}"]
    r = json.loads(p.read_text())
    why = []
    if r.get("void") is not False or r.get("commit") != E55_RUN_COMMIT:
        why.append(f"E55's reading is void or not from its run at {E55_RUN_COMMIT[:12]}")
    out = {}
    for name in TASKS:
        v = (((r.get("tasks") or {}).get(name) or {}).get("board_minus_floor") or {}).get(O)
        if v is None:
            why.append(f"{name}: no committed E55 figure for {O}")
        else:
            out[name] = float(v)
    return out, why


def e47_scores(key: str, tasks_dir: Path | None = None) -> tuple[dict | None, list[str]]:
    """The committed E47 per-seed board scores of the task (task key), and why they cannot serve."""
    for p in sorted(Path(tasks_dir or E47_DIR / "tasks").glob("T2e_*.json")):
        r = json.loads(p.read_text())
        if r["task"]["key"] == key:
            why = []
            if not r["provenance"]["commit"].startswith(E47_RUN_COMMIT):
                why.append(f"{p.name}: E47's file is not from its run at {E47_RUN_COMMIT}")
            if not r["provenance"]["tree_clean"] or r["build"] != "main" or r["protocol"] != "bake":
                why.append(f"{p.name}: E47's file is not a clean main/bake run")
            return e47.per_seed(r)[0], why
    return None, [f"{key}: no committed E47 file"]


def gate_check(name: str, key: str, run_scores: dict, o_ref: tuple[dict | None, list[str]] | None = None,
               floor_ref: tuple[dict | None, list[str]] | None = None) -> dict:
    """The J7 guard on one task: the floor against E47 within GATE_TOL per seed and per ten-seed mean, and O's ten-seed
    mean minus the floor's within O_TOL of E55's committed per-task figure. `o_ref` defaults to e55_o_per_task() and
    `floor_ref` to e47_scores(key)."""
    old, why = e47_scores(key) if floor_ref is None else floor_ref
    ref, why_ref = e55_o_per_task() if o_ref is None else o_ref
    why = [*why, *why_ref]
    worst = {"mean": None, "seed": None}
    fl = run_scores.get(harness.FLOOR)
    if old is not None:
        o = old.get(harness.FLOOR)
        if o is None or fl is None or np.shape(o) != np.shape(fl) or np.shape(fl) != (len(SEEDS),):
            why.append(f"{key}: the floor missing or of another length")
        else:
            o, n = np.asarray(o, float), np.asarray(fl, float)
            worst = {"mean": abs(float(n.mean() - o.mean())), "seed": float(np.abs(n - o).max())}
            if not worst["mean"] <= GATE_TOL:
                why.append(f"{key}: the floor's ten-seed mean is {worst['mean']:.4f} off E47's (tolerance {GATE_TOL})")
            if not worst["seed"] <= GATE_TOL:
                why.append(f"{key}: the floor is {worst['seed']:.4f} off E47's in a seed (tolerance {GATE_TOL})")
    o_minus_floor = None
    if run_scores.get(O) is not None and fl is not None \
            and np.shape(run_scores[O]) == np.shape(fl) == (len(SEEDS),):
        o_minus_floor = float(np.mean(np.asarray(run_scores[O], float) - np.asarray(fl, float)))
    else:
        why.append(f"{key}: {O} or the floor missing or of another length")
    o_off = None
    if o_minus_floor is not None:
        if ref is None or name not in ref:
            why.append(f"{name}: no committed E55 figure for {O}")
        else:
            o_off = abs(o_minus_floor - ref[name])
            if not o_off <= O_TOL:
                why.append(f"{name}: {O} minus the floor is {o_off:.4f} off E55's ten-seed mean (tolerance {O_TOL})")
    return {"task": name, "key": key, "e47_run_commit": E47_RUN_COMMIT, "e55_run_commit": E55_RUN_COMMIT,
            "e55_reading_sha256": E55_READING_SHA256, "tolerance": GATE_TOL, "o_tolerance": O_TOL, "worst": worst,
            "o_off_e55": o_off, "why": why, "pass": not why}


def _gate_path(runs: Path, name: str) -> Path:
    return Path(runs) / "gate" / f"{name}.json"


def cmd_gate(a) -> None:
    _refuse_unless_ready(a.commit)
    t = _task(a.task)
    path = _gate_path(Path(a.runs), a.task)
    if path.exists():
        sys.exit(f"{path} exists")
    res = harness.run_task(t, [harness.CANDIDATES[n] for n in GATE_RUN], seeds=SEEDS, data_dir=Path(a.data),
                           build=BUILD, log=None)
    if str(res["board"].get("split", "")).startswith("test"):
        sys.exit("E58 met a test board")
    rec = {**gate_check(a.task, t.key, e47.per_seed(res)[0]), "commit": res["provenance"]["commit"],
           "tree_clean": res["provenance"]["tree_clean"]}
    rec["pass"] = bool(rec["pass"] and rec["commit"] == a.commit and rec["tree_clean"])
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(rec, indent=1, default=float))
    print(f"gate {a.task} {'pass' if rec['pass'] else 'fail'}", flush=True)
    sys.exit(0 if rec["pass"] else 3)


def gate_problems(runs: Path, commit: str) -> list[str]:
    """Every gate file (the floor against E47, O against E55, per task) present and passing at the run commit."""
    why = []
    for name in TASKS:
        p = _gate_path(runs, name)
        if not p.exists():
            why.append(f"no gate file for {name}")
            continue
        g = json.loads(p.read_text())
        if not g.get("pass") or g.get("commit") != commit:
            why.append(f"the gate did not pass on {name} at the run commit")
    return why


def cmd_run(a) -> None:
    _refuse_unless_ready(a.commit)
    t = _task(a.task)
    runs = Path(a.runs)
    why = gate_problems(runs, a.commit)
    if why:
        sys.exit("refused (the pre-run gate): " + "; ".join(why))
    for c in (harness.CANDIDATES[O], ARM_CANDIDATE):
        if not c.applies(t):
            sys.exit(f"{a.task}: {c.name} does not apply")
    folder = runs / "run" / a.task
    folder.mkdir(parents=True, exist_ok=True)
    for s in SEEDS:
        path = folder / f"seed{s}.json"
        if path.exists():
            continue
        res = one_seed(t, s, Path(a.data))
        rec = {"task": a.task, "key": t.key, "seed": s, "run_commit": a.commit, "script": SCRIPT, **res}
        tmp = path.with_suffix(".tmp")
        tmp.write_text(json.dumps(rec, indent=1, default=float))
        os.replace(tmp, path)
        print(f"{a.task} seed {s} done", flush=True)
    print(f"wrote {folder}", flush=True)


# ------------------------------------------------------------------------------------------------ read
def series(recs: list[dict], name: str, what: str) -> np.ndarray:
    """Per seed, one quantity of one file under build main: 'board', 'group:<g>' (gated) or 'skill:<m>'."""
    rows = [r["builds"][BUILD]["rows"][0] for r in recs]
    if what == "board":
        return np.array([x["scores"][name] for x in rows], float)
    if what.startswith("group:"):
        return np.array([x["groups"][name][what[6:]] for x in rows], float)
    if what.startswith("skill:"):
        return np.array([x["skills"][name][what[6:]] for x in rows], float)
    raise ValueError(what)


def merged_D(recs: dict[str, list[dict]], groups: dict[str, list[str]], a: str, b: str, what: str) -> np.ndarray:
    """D_s: the mean over truths of each truth's task-mean difference (a minus b) at each seed."""
    return e47.merged_seeds({n: series(recs[n], a, what) - series(recs[n], b, what)
                             for names in groups.values() for n in names}, groups)


def _stat(D: np.ndarray) -> dict:
    D = np.asarray(D, float)
    return {"mean": float(D.mean()), "se": float(D.std(ddof=1) / np.sqrt(len(D))), "ahead": int((D > 0).sum()),
            "behind": int((D < 0).sum())}


def behind(D: np.ndarray) -> bool:
    """E49's P4 definition: a merged mean below 0 with at least AHEAD of the ten D_s below 0."""
    D = np.asarray(D, float)
    return bool(D.mean() < 0 and (D < 0).sum() >= AHEAD)


def decide(d: dict) -> dict:
    """THE PASS RULE (CLAUSES): `d` holds the per-seed D_s arrays board (the two-truth merged board against O) and
    groups ({group: D} for OTHER_GROUPS on the same merge), and arm_refused (the task-seeds the arm refused)."""
    D = np.asarray(d["board"], float)
    p = {"P1": bool(D.mean() >= BAR and (D > 0).sum() >= AHEAD),
         "P2": not any(behind(v) for v in d["groups"].values()) and set(d["groups"]) == set(OTHER_GROUPS),
         "P3": not d["arm_refused"]}
    idle = bool(abs(float(D.mean())) < BAR)
    passes = all(p.values()) and not idle
    return {"clauses": p, "idle": idle, "passes": passes,
            "outcome": "PASS" if passes else "STOP",
            "outcome_text": er.PASS_TEXT if passes else (er.IDLE_TEXT if idle else er.STOP_TEXT)}


def merges(tasks: dict | None = None) -> dict[str, dict[str, list[str]]]:
    """The merges fixed in row E58 (J4): counted (two truths), reported (the other four tasks on their truths) and all
    seven (three truths, new1 and x1 merged on E8.25_late), the last two deciding nothing."""
    tasks = nt.all_tasks() if tasks is None else tasks
    return {"counted": e47.truths({n: tasks[n] for n in COUNTED_TASKS}),
            "reported": e47.truths({n: tasks[n] for n in REPORTED_TASKS}),
            "all_seven": e47.truths({n: tasks[n] for n in TASKS})}


def check_run(recs: dict[str, list[dict]], why_in: list[str], commit: str) -> list[str]:
    """The void conditions on the per-seed records (task -> its records in seed order)."""
    why = list(why_in)
    tasks = nt.all_tasks()
    for task in TASKS:
        rs = recs.get(task, [])
        if [r.get("seed") for r in rs] != list(SEEDS):
            why.append(f"{task}: seeds are not {list(SEEDS)}")
            continue
        frac = frac_of(tasks[task])
        for r in rs:
            tag = f"{task} seed {r['seed']}"
            if r.get("run_commit") != commit:
                why.append(f"{tag}: another run commit")
            res = (r.get("builds") or {}).get(BUILD)
            if res is None:
                why.append(f"{tag}: build {BUILD} missing")
                continue
            p = res["provenance"]
            if p["commit"] != commit or not p["tree_clean"]:
                why.append(f"{tag}: another commit or an unclean tree")
            if res["build"] != BUILD or res["protocol"] != "bake" or [x["seed"] for x in res["rows"]] != [r["seed"]]:
                why.append(f"{tag}: not the bake protocol under build {BUILD} at this seed")
            if str(res["board"].get("split", "")).startswith("test"):
                why.append(f"{tag}: a test board")
            row = res["rows"][0]
            if any(n not in row["scores"] or n not in row["skills"] or n not in row.get("groups", {})
                   for n in (harness.FLOOR, O, ARM)):
                why.append(f"{tag}: an arm is missing")
            f = r["facts"]
            if str(f.get("board_split", "")).startswith("test"):
                why.append(f"{tag}: a test board")
            if f["n"].get(ARM) != f["n"].get(O):
                why.append(f"{tag}: the arm's count differs from O's")
            ro, ra = f["records"].get(O) or {}, f["records"].get(ARM) or {}
            if ro.get("t_mix") is None or abs(ro["t_mix"] - frac) > T_MIX_TOL:
                why.append(f"{tag}: O's recorded t_mix is not the task's t")
            if not f["refused"].get(ARM):  # a refused arm's file is O's; its own record is not read
                if ra.get("t_mix") is None or abs(ra["t_mix"] - arm_t_mix(frac)) > T_MIX_TOL:
                    why.append(f"{tag}: the arm's recorded t_mix is not min(1, t + 1/6)")
                if ro.get("target_rms") is None or ro.get("target_rms") != ra.get("target_rms"):
                    why.append(f"{tag}: the arm's size target differs from O's")
    return why


def _per_truth(recs, groups, a: str, b: str) -> dict:
    out = {truth: float(np.mean([np.mean(series(recs[n], a, "board") - series(recs[n], b, "board")) for n in names]))
           for truth, names in groups.items()}
    return {"per_truth": out, "gains": {k: v for k, v in out.items() if v > 0},
            "losses": {k: v for k, v in out.items() if v < 0}}


def reading(recs: dict[str, list[dict]], why_in: list[str], commit: str) -> dict:
    why = check_run(recs, why_in, commit)
    if why:
        return {"void": True, "why": why, "commit": commit, "outcome": "VOID", "outcome_text": er.VOID_TEXT,
                "passes": False}
    m = merges()
    refused = sum(bool(r["facts"]["refused"].get(ARM)) for t in TASKS for r in recs[t])
    d = {"board": merged_D(recs, m["counted"], ARM, O, "board"),
         "groups": {g: merged_D(recs, m["counted"], ARM, O, f"group:{g}") for g in OTHER_GROUPS},
         "arm_refused": refused}
    dec = decide(d)
    per_task, tasks = {}, nt.all_tasks()
    for t in TASKS:
        x = series(recs[t], ARM, "board") - series(recs[t], O, "board")
        frac = frac_of(tasks[t])
        per_task[t] = {"board": {"mean": float(x.mean()), "ahead": int((x > 0).sum()), "sd": float(x.std(ddof=1))},
                       "board_minus_floor": {n: float(np.mean(series(recs[t], n, "board")
                                                              - series(recs[t], harness.FLOOR, "board")))
                                             for n in (O, ARM)},
                       "t_mix": {"O": frac, "arm": arm_t_mix(frac)}, "counted": t in COUNTED_TASKS}
    sg = merged_D(recs, m["counted"], ARM, O, f"group:{GROUP_SHAPE}")
    return {"void": False, "why": [], "commit": commit, "row": TRY_ID, "arm": ARM, "control": O,
            "truths": m, "clauses": dec["clauses"], "clause_text": CLAUSES, "pass_rule": PASS_RULE,
            "idle": dec["idle"], "passes": dec["passes"], "outcome": dec["outcome"],
            "outcome_text": dec["outcome_text"], "arm_refused_task_seeds": refused,
            "merged": {"board": _stat(d["board"]),
                       "groups": {g: {**_stat(v), "is_behind": behind(v)} for g, v in d["groups"].items()}},
            "reported_deciding_nothing": {
                "shape_growth_counted_merge_recused": _stat(sg),
                "skills_counted_merge": {k: _stat(merged_D(recs, m["counted"], ARM, O, f"skill:{k}"))
                                         for k in SHAPE_TERMS},
                "reported_tasks_merge": _stat(merged_D(recs, m["reported"], ARM, O, "board")),
                "all_seven_merge": _stat(merged_D(recs, m["all_seven"], ARM, O, "board")),
                "against_floor_counted": {n: _stat(merged_D(recs, m["counted"], n, harness.FLOOR, "board"))
                                          for n in (O, ARM)}},
            "price": _per_truth(recs, m["counted"], ARM, O), "price_lines": PRICE_LINES, "tasks": per_task,
            "barred_not_run": list(BARRED_TASKS),
            "gate": {"e47_run_commit": E47_RUN_COMMIT, "tolerance": GATE_TOL, "e55_run_commit": E55_RUN_COMMIT,
                     "e55_reading_sha256": E55_READING_SHA256, "o_tolerance": O_TOL}}


def text(res: dict) -> str:
    if res["void"]:
        return (f"Row E58 (run commit {res['commit']}): VOID: " + "; ".join(res["why"]) + "\n" + er.VOID_TEXT + "\n")
    m, rep = res["merged"], res["reported_deciding_nothing"]
    L = [f"Row E58 (run commit {res['commit']}): gate G1 of row R59's file T2embryo_s2x_ovt23.h5ad; it moves the "
         "overlay's form; one counted try, the arm fgw_shape0x_ov_t23; diagnostic (J5, noise floor 6.54).",
         f"The pre-run gate passed on all seven tasks (the floor under main within {res['gate']['tolerance']} of E47, "
         f"run at {res['gate']['e47_run_commit']}; O minus the floor within {res['gate']['o_tolerance']} of E55's "
         "committed per-task ten-seed mean). e75 and e725 were not run (barred for any change to the mixing fraction).",
         "The pass rule, fixed in row E58 before any output, build main, against O on the two-truth merge (E8.0 = "
         "half, late; E8.25_late = new1); every clause must hold and the gain must not be idle:"]
    for k, v in res["clauses"].items():
        L.append(f"  {k} {'holds' if v else 'fails'}: {res['clause_text'][k]}")
    b = m["board"]
    L.append(f"  merged board {b['mean']:+.3f} (SE {b['se']:.3f}, {b['ahead']} of 10 seeds ahead)"
             + ("; idle (below 0.3 in size)" if res["idle"] else ""))
    L.append("  groups on that merge: " + "; ".join(f"{g} {v['mean']:+.4f} ({v['behind']} behind)"
                                                 for g, v in m["groups"].items()))
    if res["arm_refused_task_seeds"]:
        L.append(f"  the arm refused on {res['arm_refused_task_seeds']} task-seeds (O's file stood in): unreadable, "
                 "STOP, not void; the try stays counted")
    L.append("Reported, deciding nothing: shape_growth on the counted merge (the main session is recused) "
             f"{rep['shape_growth_counted_merge_recused']['mean']:+.4f}; per counted truth "
             + ", ".join(f"{k} {v:+.3f}" for k, v in res["price"]["per_truth"].items())
             + f"; the reported tasks' merge (new2, x1, x2, x3) {rep['reported_tasks_merge']['mean']:+.3f}; all seven "
             f"{rep['all_seven_merge']['mean']:+.3f}; against the floor (rule 5) O "
             f"{rep['against_floor_counted'][O]['mean']:+.3f}, the arm "
             f"{rep['against_floor_counted'][ARM]['mean']:+.3f}.")
    for t, v in res["tasks"].items():
        L.append(f"  {t}{' (counted)' if v['counted'] else ''}: the arm minus O {v['board']['mean']:+.3f} "
                 f"({v['board']['ahead']} ahead); t_mix O {v['t_mix']['O']:.4f}, the arm {v['t_mix']['arm']:.4f}")
    L.extend(res["price_lines"])
    if res["passes"]:
        L.append(f"  file: {CANDIDATE_COMMAND}")
    L.append(res["outcome_text"])
    L.append("Nothing here shows the test's world, E7.75's form, the union's form at the test pair, the live gate, "
             "build or metric seed.")
    return "\n".join(L) + "\n"


def load_runs(runs: Path) -> dict[str, list[dict]]:
    out = {}
    for name in TASKS:
        folder = Path(runs) / "run" / name
        recs = [json.loads(p.read_text()) for p in sorted(folder.glob("seed*.json"))] if folder.exists() else []
        out[name] = sorted(recs, key=lambda r: r["seed"])
    return out


def provenance() -> dict:
    s = repo_state(REPO)
    return {"head": s.commit, "clean": s.unclean == (), "script_sha256": _sha_file(Path(__file__)),
            "row_sha256": _sha_file(Path(er.__file__)), "harness_sha256": _sha_file(Path(harness.__file__)),
            "t2_build_sha256": _sha_file(Path(tb.__file__)), "shape_interp_sha256": _sha_file(Path(si.__file__))}


def cmd_read(a) -> None:
    out = Path(a.out)
    files = [out / "reading.json", out / "reading.txt"]
    existing = [str(f) for f in files if f.exists()]
    if existing:
        sys.exit(f"refusing to overwrite: {', '.join(existing)}")
    _own_code()
    why = er.row_problems(TRIES)
    if why:
        sys.exit("E58 has no well-formed row in experiments/harness/TRIES.md: " + "; ".join(why))
    runs = Path(a.runs)
    why = [*gate_problems(runs, a.commit), *commit_problems(a.commit)]
    if repo_state(REPO).unclean != ():
        why.append("the tree is not clean at the read")
    res = {**reading(load_runs(runs), why, a.commit), "reading_code": provenance()}
    out.mkdir(parents=True, exist_ok=True)
    files[0].write_text(json.dumps(res, indent=1, default=float))
    files[1].write_text(text(res))
    print(f"wrote {files[1]}", flush=True)


# ------------------------------------------------------------------------------------------------ J7 dry run
def j7_dry_run(e55_path: Path | None = None, e49_path: Path | None = None, e47_dir: Path | None = None,
               e55_sha: str | None = E55_READING_SHA256, e49_sha: str | None = E49_READING_SHA256) -> dict:
    """The O guard's tolerance tried on committed reproduction data other than the pair G1 reads (rule 9, J7), from
    committed records only (not a try): per task, E55's committed O-minus-floor ten-seed mean (its run at 4535f98)
    against E49's committed O-minus-C0 ten-seed mean (its run at d76ccd5) plus E47's committed C0-minus-floor ten-seed
    mean (its run at 9967cf1). The guard could fail there: three runs on three commits."""
    ref, why = e55_o_per_task(e55_path, e55_sha)
    p49 = E49_READING if e49_path is None else Path(e49_path)
    e49 = None
    if not p49.exists() or (e49_sha is not None and _sha_file(p49) != e49_sha):
        why = [*why, "E49's committed reading is missing or changed"]
    else:
        e49 = json.loads(p49.read_text())
        if e49.get("void") is not False or not str(e49.get("commit", "")).startswith(E49_RUN_COMMIT):
            why = [*why, f"E49's reading is void or not from its run at {E49_RUN_COMMIT}"]
    devs = {}
    tasks = nt.all_tasks()
    for name in TASKS:
        old, w = e47_scores(tasks[name].key, e47_dir)
        why = [*why, *w]
        o49 = (((e49 or {}).get("tasks") or {}).get(name) or {}).get(O, {}).get("board", {}).get("mean")
        if old is None or o49 is None or ref is None or name not in ref or C0 not in old \
                or harness.FLOOR not in old:
            why.append(f"{name}: a committed figure is missing")
            continue
        c0_minus_floor = float(np.mean(np.asarray(old[C0], float) - np.asarray(old[harness.FLOOR], float)))
        devs[name] = abs(ref[name] - (float(o49) + c0_minus_floor))
    worst = max(devs.values()) if devs else None
    ok = not why and len(devs) == len(TASKS) and worst is not None and worst <= O_TOL
    return {"pass": bool(ok), "tolerance": O_TOL, "worst": worst, "per_task": devs, "why": why,
            "records": {"e55": [E55_RUN_COMMIT, E55_READING_SHA256], "e49": [E49_RUN_COMMIT, E49_READING_SHA256],
                        "e47": E47_RUN_COMMIT},
            "what": "per task |E55's O minus the floor - (E49's O minus C0 + E47's C0 minus the floor)|, ten-seed "
                    "means under build main from committed records; not a try"}


def cmd_j7(a) -> None:
    res = j7_dry_run()
    if a.out:
        Path(a.out).parent.mkdir(parents=True, exist_ok=True)
        Path(a.out).write_text(json.dumps(res, indent=1, default=float))
    print(f"j7 {'pass' if res['pass'] else 'fail'}", flush=True)
    sys.exit(0 if res["pass"] else 3)


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("gate", "run"):
        p = sub.add_parser(name)
        p.add_argument("--task", required=True, help=f"one of {', '.join(TASKS)}")
        p.add_argument("--commit", required=True, help="the run commit (full hash)")
        p.add_argument("--runs", default=str(RUNS))
        p.add_argument("--data", default=str(harness.DATA_DIR))
    r = sub.add_parser("read")
    r.add_argument("--commit", required=True)
    r.add_argument("--runs", default=str(RUNS))
    r.add_argument("--out", default=str(RESULTS))
    j = sub.add_parser("j7")
    j.add_argument("--out", default=None, help="a JSON file for the dry run's per-task deviations")
    a = ap.parse_args(argv)
    {"gate": cmd_gate, "run": cmd_run, "read": cmd_read, "j7": cmd_j7}[a.cmd](a)


if __name__ == "__main__":
    main()
