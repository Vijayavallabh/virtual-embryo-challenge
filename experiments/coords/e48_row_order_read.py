#!/usr/bin/env python
"""Try E48 of experiments/harness/TRIES.md (committed in 1420c70, before this script): the pre-fix price of the floor
copy's parity on the embryo's floor-form file, on the nine practice tasks.

The arm. `fgw_shape0x_series` built per task and seed through the harness's own candidate at BUILD_SEED_OFFSET + seed,
then its rows re-ordered by vec.t2_build.copy_parity_order against the bake copy of the task's reference stage at that
seed (vec.evaluate.subsample(reference as released, 0.1, seed): the coordinates harness.score_seeds scores as floor),
with n_true the practice board's truth_cells and metric seed = the scoring seed. The same for `fgw_shape0x` (the row's
"fgw_shape0x in both orders"). Controls: both as built. A refusal of copy_parity_order (the copy and the file would get
different truth rows, an undefined parity, or no permutation found) leaves the file as built, as the builder would ship
it: the arm then equals its control on that task-seed.

Path. harness.run_task on each task of nine_tasks.all_tasks() (half, late, new1, new2, e75, e725, x1, x2, x3), seeds
0-9, bake protocol, once with build "main" and once with build "prefix" (the floor checked to score 50 under each).
Nothing of the task cut, the seeds or the read guard is re-implemented: the arms are harness Candidates whose build
wraps the base candidate's; a spy on harness.score_seeds and harness.subsample hands them the board, the reference
stage and the floor copy the run itself uses. A task with a target filter is refused.

Per task and seed it records, beside run_task's two results (full ranked panel and board score of every arm under
both builds): whether each file and the floor copy are handed the same truth rows (truth_rows_shared); whether the
option refused (no values); for every file and the floor copy, vec.evaluate.handedness against the fold truth half
(prefix mode, proper and mirror Dice; a counted read, under row E48).

Reading, fixed in row E48 before this script, applied by `read`:
  (i)   on every task-seed where the arm's file and the floor copy share their truth rows, the arm's prefix mode must
        equal the floor's; any difference voids the construction (not the run). Reported for both base candidates;
        the row's arm is fgw_shape0x_series_cp.
  (ii)  per task, the paired arm-minus-control occupancy skill under prefix and under main (mean, seed SD) and the
        paired board score difference, for both base candidates; every ranked metric's paired difference is reported.
        Primary figure: every seed, a refused task-seed counting as the control (difference 0), because the builder
        then ships the as-built file. Secondary: refused task-seeds left out, with their count. Only the primary is
        read under (ii).
  (iii) per task, the proper-minus-mirror Dice gap (the coin's stake), for the floor copy and every file.
Also reported: the number of shared task-seeds (if none, the read cannot price the option, and it says so); as built,
the share of shared seeds where the control's pre-fix parity agrees with the floor's (reading B.1, a fraction only);
the nine-task merged means (tasks sharing a truth merged, as E47) with their SE over seeds (reported, read by nothing).
Note (reviewer point 5): on the option's own rows, truth_rows_shared and the parity agreement re-test the option's gate;
only the prefix modes scored against the fold truth test the mechanism.

Void (the run is read as nothing; no task or seed is re-run selectively):
  - the pre-run gate did not pass on all nine tasks: the two controls under build main, all ten seeds, must reproduce
    the committed E47 scores (experiments/harness/results/e47_floor_exact, run at 9967cf1 on PrakashDGX_H2 in the
    campaign checkout ~/v/vec-campaign) within 0.1 board points in a task's ten-seed mean and 0.5 in any seed. This
    gate is the analyst's addition after the E47 guard (J7), not text of row E48; it runs before the arms start;
  - a result's commit differs from the run commit named at launch, or its tree is not clean; the run commit must
    contain 1420c70 (the row), 616984b and the test-board guard (--row-order-test-ok), and src/vec/t2_build.py,
    src/vec/evaluate.py, src/vec/harness.py and experiments/harness/nine_tasks.py must equal 616984b's apart from
    docstrings (compared as syntax trees);
  - a task, seed, build or arm is missing;
  - the copy handed to copy_parity_order is not bitwise the floor copy score_seeds scores at that seed, or the metric
    seed is not the scoring seed;
  - an arm changes anything but row order: its multiset of (expression row, position) pairs, and its count, differ
    from its control's;
  - the n_true handed to copy_parity_order is not the fold truth half's count;
  - a test board is met (the practice boards are "pseudo" or validation boards; --row-order copy_parity never runs on
    a test board here).
Nothing is selected or tuned on the output. The script prints no test count; it reads no test board, test file or test
anchor.

  gate  --task NAME --commit SHA   the controls under main against E47 (one task per process; start the nine together)
  run   --task NAME --commit SHA   the arms and controls, both builds, seed by seed (resumable: a seed whose file exists
                                   is skipped); refused unless all nine gate files pass at the same commit
  read  --commit SHA               the reading into experiments/coords/results/e48_row_order/ (reading.json, .txt)
Refused: a run while TRIES.md has no row E48, an unclean tree, another HEAD, or code not from this checkout.
"""

from __future__ import annotations

import argparse
import ast
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
from vec import t2_build as tb  # noqa: E402
from vec.adoption import SEEDS  # noqa: E402
from vec.provenance import repo_state  # noqa: E402

ROW = "E48"
TRIES = REPO / "experiments" / "harness" / "TRIES.md"
CONTROLS = ("fgw_shape0x_series", "fgw_shape0x")
ARMS = {"fgw_shape0x_series_cp": "fgw_shape0x_series", "fgw_shape0x_cp": "fgw_shape0x"}  # arm -> its control
ROW_ARM = "fgw_shape0x_series_cp"
FILES = (*CONTROLS, *ARMS)
BUILDS = ("main", "prefix")
OCC = "occupancy_dice"
BASES = ("1420c70", "616984b")
FROZEN_AT = "616984b"
FROZEN = ("src/vec/t2_build.py", "src/vec/evaluate.py", "src/vec/harness.py", "experiments/harness/nine_tasks.py")
GUARD_FLAG = "--row-order-test-ok"
E47_TASKS = REPO / "experiments" / "harness" / "results" / "e47_floor_exact" / "tasks"
E47_COMMIT = "9967cf1"
E47_WHERE = "PrakashDGX_H2, campaign checkout ~/v/vec-campaign"
GUARD_MEAN, GUARD_SEED = e47.GUARD_MEAN, e47.GUARD_SEED  # 0.1 and 0.5 board points
RUNS = REPO / "runs" / "e48"
RESULTS = REPO / "experiments" / "coords" / "results" / "e48_row_order"


# ------------------------------------------------------------------------------------------------ checks
def _git(*args: str) -> subprocess.CompletedProcess:
    return subprocess.run(["git", "-C", str(REPO), *args], capture_output=True, text=True,
                          env={**os.environ, "GIT_OPTIONAL_LOCKS": "0"})


def strip_docstrings(src: str) -> str:
    """The syntax tree of a module with every module, class and function docstring removed, dumped."""
    tree = ast.parse(src)
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)) and node.body:
            first = node.body[0]
            if isinstance(first, ast.Expr) and isinstance(first.value, ast.Constant) and isinstance(first.value.value,
                                                                                                    str):
                node.body = node.body[1:] or [ast.Pass()]
    return ast.dump(tree)


def frozen_changes(old_text=None, new_text=None) -> list[str]:
    """The frozen files whose code (docstrings aside) differs from FROZEN_AT's. old_text / new_text: path -> source,
    for tests; by default git show FROZEN_AT:path and the working file."""
    out = []
    for path in FROZEN:
        old = old_text(path) if old_text else _git("show", f"{FROZEN_AT}:{path}").stdout
        new = new_text(path) if new_text else (REPO / path).read_text(encoding="utf-8")
        try:
            if strip_docstrings(old) != strip_docstrings(new):
                out.append(path)
        except SyntaxError:
            out.append(path)
    return out


def commit_problems(commit: str) -> list[str]:
    """Why this checkout may not run or read try E48 at `commit` (empty: it may)."""
    why = []
    head = _git("rev-parse", "HEAD").stdout.strip()
    if head != commit:
        why.append(f"HEAD is {head[:12]}, not the run commit {commit[:12]}")
    for base in BASES:
        if _git("merge-base", "--is-ancestor", base, "HEAD").returncode != 0:
            why.append(f"HEAD does not contain {base}")
    changed = frozen_changes()
    if changed:
        why.append(f"code differs from {FROZEN_AT} (docstrings aside): {changed}")
    if GUARD_FLAG not in (REPO / "experiments" / "builds" / "build_test.py").read_text(encoding="utf-8"):
        why.append(f"the builder lacks the test-board guard {GUARD_FLAG}")
    return why


def _own_code() -> None:
    src = (REPO / "src").resolve()
    for mod in (ev, harness, tb):
        if src not in Path(mod.__file__).resolve().parents:
            sys.exit(f"{mod.__name__} is loaded from {mod.__file__}, not from this checkout's src/")


def _refuse_unless_ready(commit: str) -> None:
    _own_code()
    if ROW not in n12.tries_rows(TRIES):
        sys.exit(f"{ROW} has no row in experiments/harness/TRIES.md: a row is committed before its run")
    s = repo_state(REPO)
    if s.unclean != ():
        sys.exit("the tree is not clean: try E48 runs on a committed tree")
    why = commit_problems(commit)
    if why:
        sys.exit("refused: " + "; ".join(why))


def _task(name: str):
    tasks = nt.all_tasks()
    if name not in tasks:
        sys.exit(f"unknown task {name}; known: {sorted(tasks)}")
    t = tasks[name]
    if getattr(t, "target_filter", None) is not None:
        sys.exit(f"{name}: a task with a target filter; E48 does not run on one")
    return t


def _is_test_board(board) -> bool:
    return str(board.split).startswith("test") or "test" in str(board.key).split(":")[-1]


# ------------------------------------------------------------------------------------------------ the arms
def row_digests(stage) -> list[bytes]:
    """The multiset of a built file's (expression row, position) pairs, as sorted per-row digests (float32 bytes)."""
    X = np.ascontiguousarray(stage.X, dtype=np.float32)
    C = np.ascontiguousarray(stage.coords, dtype=np.float32)
    return sorted(hashlib.blake2b(X[i].tobytes() + C[i].tobytes(), digest_size=16).digest() for i in range(len(X)))


@contextlib.contextmanager
def captured():
    """While active, record what harness.score_seeds is handed (task, board, target, reference stage) and the floor
    copy it draws per seed (harness.subsample of that reference at SUBSAMPLE_FRAC), for the wrapped candidates."""
    seen: dict = {"floor": {}, "coords": {}, "digest": {}, "arm": {}}
    orig_ss, orig_sub = harness.score_seeds, harness.subsample

    def score_seeds(t, board, target, ref, inputs, candidates, *a, **kw):
        seen.update(task=t, board=board, target=target, ref=ref)
        return orig_ss(t, board, target, ref, inputs, candidates, *a, **kw)

    def subsample(stage, frac, seed):
        out = orig_sub(stage, frac, seed)
        if stage is seen.get("ref") and frac == ev.SUBSAMPLE_FRAC:
            seen["floor"][int(seed)] = np.array(out.coords, copy=True)
        return out

    harness.score_seeds, harness.subsample = score_seeds, subsample
    try:
        yield seen
    finally:
        harness.score_seeds, harness.subsample = orig_ss, orig_sub


def observed(base, seen):
    """The base candidate unchanged; its built file's coordinates and pair digests are kept for the checks."""
    def build(ctx):
        pred = base.build(ctx)
        s = ctx.seed - harness.BUILD_SEED_OFFSET
        seen["coords"][(base.name, s)] = np.array(pred.coords, copy=True)
        seen["digest"][(base.name, s)] = row_digests(pred)
        return pred

    return dataclasses.replace(base, build=build)


def reordered(name: str, base, seen):
    """The base candidate's file with its rows re-ordered by copy_parity_order against the floor copy of the seed."""
    def build(ctx):
        pred = base.build(ctx)
        s = ctx.seed - harness.BUILD_SEED_OFFSET
        board, ref = seen["board"], seen["ref"]
        if _is_test_board(board):
            raise RuntimeError("E48 met a test board")
        copy = ev.subsample(ref, ev.SUBSAMPLE_FRAC, s).coords
        floor = seen["floor"].get(s)
        rec = {"metric_seed": s, "scoring_seed_known": s in seen["floor"],
               "copy_is_floor": bool(floor is not None and copy.dtype == floor.dtype and copy.shape == floor.shape
                                     and copy.tobytes() == floor.tobytes()),
               "n_true_handed": None}
        try:
            order, _ = tb.copy_parity_order(pred.coords, copy, board.truth_cells, ctx.seed, metric_seed=s)
            out, rec["refused"] = pred.take(order), False
            rec["permutation"] = bool(np.array_equal(np.sort(order), np.arange(pred.n)))
        except ValueError:
            out, rec["refused"], rec["permutation"] = pred, True, True
        rec["n_true_handed"] = int(board.truth_cells)
        seen["arm"][(name, s)] = rec
        seen["coords"][(name, s)] = np.array(out.coords, copy=True)
        seen["digest"][(name, s)] = row_digests(out)
        return out

    params = {**dict(base.params), "row_order": tb.ROW_ORDERS[1], "row_order_stream": tb.COPY_PARITY_STREAM,
              "row_order_copy": "vec.evaluate.subsample(reference as released, 0.1, scoring seed)",
              "row_order_metric_seed": "the scoring seed", "try": ROW}
    return dataclasses.replace(base, name=name, build=build, params=params)


def candidates(seen) -> list:
    return [harness.CANDIDATES[harness.FLOOR], *(observed(harness.CANDIDATES[c], seen) for c in CONTROLS),
            *(reordered(a, harness.CANDIDATES[c], seen) for a, c in ARMS.items())]


def _hand(coords, truth, s: int) -> dict:
    h = ev.handedness(coords, truth.coords, seed=s, live_n_true=truth.n)
    return {k: h[k] for k in ("prefix_mode", "dice_proper", "dice_mirror", "dice_gap", "parity_pred")}


def one_seed(t, s: int, data: Path) -> dict:
    """Both builds of one task-seed: run_task's results and E48's per-seed facts."""
    out: dict = {"builds": {}, "facts": {}}
    digests = {}
    for build in BUILDS:
        with captured() as seen:
            res = harness.run_task(t, candidates(seen), seeds=[s], data_dir=data, build=build, log=None)
        board = seen["board"]
        if _is_test_board(board):
            raise RuntimeError("E48 met a test board")
        out["builds"][build] = res
        digests[build] = {n: seen["digest"][(n, s)] for n in FILES}
        if build != BUILDS[0]:
            continue
        truth = ev.working_copies(seen["target"], seen["ref"], s, ev.SUBSAMPLE_FRAC, harness.PROTOCOL)[0]
        floor = seen["floor"][s]
        hands = {harness.FLOOR: _hand(floor, truth, s), **{n: _hand(seen["coords"][(n, s)], truth, s) for n in FILES}}
        shared = {n: tb.truth_rows_shared(len(seen["coords"][(n, s)]), len(floor), truth.n, s) for n in FILES}
        arms = {}
        for a, c in ARMS.items():
            r = seen["arm"][(a, s)]
            arms[a] = {**{k: v for k, v in r.items() if k != "n_true_handed"},
                       "n_true_is_truth_half": r["n_true_handed"] == truth.n,
                       "pairs_equal_control": seen["digest"][(a, s)] == seen["digest"][(c, s)]}
        out["facts"] = {"board_split": str(board.split), "hands": hands, "truth_rows_shared": shared, "arms": arms}
    out["facts"]["files_identical_across_builds"] = all(digests[BUILDS[0]][n] == digests[b][n]
                                                        for b in BUILDS[1:] for n in FILES)
    return out


# ------------------------------------------------------------------------------------------------ steps
def _gate_path(runs: Path, name: str) -> Path:
    return runs / "gate" / f"{name}.json"


def e47_scores(key: str) -> tuple[dict | None, list[str]]:
    """The committed E47 per-seed board scores of the task (task key), and why they cannot serve."""
    for p in sorted(E47_TASKS.glob("T2e_*.json")):
        r = json.loads(p.read_text())
        if r["task"]["key"] == key:
            why = []
            if not r["provenance"]["commit"].startswith(E47_COMMIT):
                why.append(f"{p.name}: E47's file is not from {E47_COMMIT}")
            if not r["provenance"]["tree_clean"] or r["build"] != "main" or r["protocol"] != "bake":
                why.append(f"{p.name}: E47's file is not a clean main/bake run")
            return e47.per_seed(r)[0], why
    return None, [f"{key}: no committed E47 file"]


def gate_check(name: str, key: str, run_scores: dict) -> dict:
    old, why = e47_scores(key)
    worst = {"mean": None, "seed": None}
    if old is not None:
        worst, g = e47.guard({key: run_scores}, {key: old}, [key], arms=CONTROLS)
        why = [*why, *g]
    return {"task": name, "key": key, "e47": {"commit": E47_COMMIT, "where": E47_WHERE, "dir": str(
        E47_TASKS.relative_to(REPO))}, "tolerance": {"mean": GUARD_MEAN, "seed": GUARD_SEED}, "worst": worst,
            "why": why, "pass": not why}


def cmd_gate(a) -> None:
    _refuse_unless_ready(a.commit)
    t = _task(a.task)
    runs = Path(a.runs)
    path = _gate_path(runs, a.task)
    if path.exists():
        sys.exit(f"{path} exists")
    cands = [harness.CANDIDATES[harness.FLOOR], *(harness.CANDIDATES[c] for c in CONTROLS)]
    res = harness.run_task(t, cands, seeds=SEEDS, data_dir=Path(a.data), build="main", log=None)
    rec = {**gate_check(a.task, t.key, e47.per_seed(res)[0]), "commit": res["provenance"]["commit"],
           "tree_clean": res["provenance"]["tree_clean"]}
    rec["pass"] = bool(rec["pass"] and rec["commit"] == a.commit and rec["tree_clean"])
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(rec, indent=1, default=float))
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


def cmd_run(a) -> None:
    _refuse_unless_ready(a.commit)
    t = _task(a.task)
    runs = Path(a.runs)
    why = gate_problems(runs, a.commit)
    if why:
        sys.exit("refused (the pre-run gate): " + "; ".join(why))
    for c in (harness.FLOOR, *CONTROLS):
        if c != harness.FLOOR and not harness.CANDIDATES[c].applies(t):
            sys.exit(f"{a.task}: {c} does not apply")
    folder = runs / "run" / a.task
    folder.mkdir(parents=True, exist_ok=True)
    for s in SEEDS:
        path = folder / f"seed{s}.json"
        if path.exists():
            continue
        res = one_seed(t, s, Path(a.data))
        rec = {"task": a.task, "key": t.key, "seed": s, "run_commit": a.commit, "script": "experiments/coords/"
               "e48_row_order_read.py", **res}
        tmp = path.with_suffix(".tmp")
        tmp.write_text(json.dumps(rec, indent=1, default=float))
        os.replace(tmp, path)
        print(f"{a.task} seed {s} done", flush=True)
    print(f"wrote {folder}", flush=True)


# ------------------------------------------------------------------------------------------------ read
def _seed_values(recs: list[dict], build: str, name: str) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    rows = [r["builds"][build]["rows"][0] for r in recs]
    board = np.array([x["scores"][name] for x in rows], float)
    skills = {m: np.array([x["skills"][name][m] for x in rows], float) for m in rows[0]["skills"][name]}
    return board, skills


def _stat(d: np.ndarray) -> dict:
    return {"mean": float(d.mean()) if d.size else None, "sd": float(d.std(ddof=1)) if d.size > 1 else None,
            "n": int(d.size)}


def check_run(recs: dict[str, list[dict]], gates_why: list[str], commit: str) -> list[str]:
    """The void conditions on the per-seed records (task -> its records in seed order)."""
    why = list(gates_why)
    for name in nt.all_tasks():
        rs = recs.get(name, [])
        if [r["seed"] for r in rs] != list(SEEDS):
            why.append(f"{name}: seeds are not {list(SEEDS)}")
            continue
        for r in rs:
            tag = f"{name} seed {r['seed']}"
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
                row = res["rows"][0]
                if any(n not in row["scores"] or n not in row["skills"] for n in (harness.FLOOR, *FILES)):
                    why.append(f"{tag}: build {b} lacks an arm")
                if str(res["board"].get("split", "")).startswith("test"):
                    why.append(f"{tag}: a test board")
            f = r["facts"]
            if str(f.get("board_split", "")).startswith("test"):
                why.append(f"{tag}: a test board")
            for a, x in f["arms"].items():
                if not (x["copy_is_floor"] and x["scoring_seed_known"] and x["metric_seed"] == r["seed"]):
                    why.append(f"{tag}: {a}'s copy is not the scored floor copy, or its metric seed not the scoring "
                               "seed")
                if not (x["permutation"] and x["pairs_equal_control"]):
                    why.append(f"{tag}: {a} changes more than row order")
                if not x["n_true_is_truth_half"]:
                    why.append(f"{tag}: {a}'s n_true is not the truth half's count")
    return why


def reading(recs: dict[str, list[dict]], gates_why: list[str], commit: str) -> dict:
    why = check_run(recs, gates_why, commit)
    if why:
        return {"void": True, "why": why, "commit": commit}
    tasks = nt.all_tasks()
    groups = e47.truths(tasks)
    out_tasks, mode_fail, shared_n, refused_n = {}, {a: [] for a in ARMS}, {a: 0 for a in ARMS}, {a: 0 for a in ARMS}
    agree = {}
    for name in tasks:
        rs = recs[name]
        facts = [r["facts"] for r in rs]
        per: dict = {"paired": {}, "gap": {}, "shared_seeds": {}, "refused_seeds": {}}
        for a, c in ARMS.items():
            sh = [f["truth_rows_shared"][a] for f in facts]
            ref = [f["arms"][a]["refused"] for f in facts]
            per["shared_seeds"][a], per["refused_seeds"][a] = int(sum(sh)), int(sum(ref))
            shared_n[a] += int(sum(sh))
            refused_n[a] += int(sum(ref))
            for r, f, s_ in zip(rs, facts, sh, strict=True):
                if s_ and f["hands"][a]["prefix_mode"] != f["hands"][harness.FLOOR]["prefix_mode"]:
                    mode_fail[a].append(f"{name} seed {r['seed']}")
            keep = ~np.array(ref, bool)
            per["paired"][a] = {"control": c}
            for b in BUILDS:
                ba, ka = _seed_values(rs, b, a)
                bc, kc = _seed_values(rs, b, c)
                per["paired"][a][b] = {
                    "board": _stat(ba - bc), "occupancy_skill": _stat(ka[OCC] - kc[OCC]),
                    "skills": {m: _stat(ka[m] - kc[m]) for m in ka},
                    "secondary": {"board": _stat((ba - bc)[keep]), "occupancy_skill": _stat((ka[OCC] - kc[OCC])[keep]),
                                  "refused_left_out": int((~keep).sum())}}
        for n in (harness.FLOOR, *FILES):
            g = np.array([f["hands"][n]["dice_gap"] for f in facts], float)
            per["gap"][n] = {"dice_proper": float(np.mean([f["hands"][n]["dice_proper"] for f in facts])),
                             "dice_mirror": float(np.mean([f["hands"][n]["dice_mirror"] for f in facts])),
                             "gap": _stat(g)}
        for c in CONTROLS:
            sh = [f["truth_rows_shared"][c] for f in facts]
            ag = [f["hands"][c]["parity_pred"] == f["hands"][harness.FLOOR]["parity_pred"]
                  for f, s_ in zip(facts, sh, strict=True) if s_]
            agree.setdefault(c, {})[name] = None if not ag else round(float(np.mean(ag)), 3)
        out_tasks[name] = per
    merged = {}
    for a in ARMS:
        merged[a] = {}
        for b in BUILDS:
            d_board = {n: _seed_values(recs[n], b, a)[0] - _seed_values(recs[n], b, ARMS[a])[0] for n in tasks}
            d_occ = {n: _seed_values(recs[n], b, a)[1][OCC] - _seed_values(recs[n], b, ARMS[a])[1][OCC] for n in tasks}
            mb, mo = e47.merged_seeds(d_board, groups), e47.merged_seeds(d_occ, groups)
            merged[a][b] = {"board_mean": float(mb.mean()), "board_se": float(mb.std(ddof=1) / np.sqrt(len(SEEDS))),
                            "occupancy_mean": float(mo.mean()),
                            "occupancy_se": float(mo.std(ddof=1) / np.sqrt(len(SEEDS)))}
    rule_i = {a: {"shared_task_seeds": shared_n[a], "mode_differs_on": mode_fail[a],
                  "construction_void": bool(mode_fail[a]), "can_price": shared_n[a] > 0} for a in ARMS}
    return {"void": False, "why": [], "commit": commit, "row_arm": ROW_ARM, "truths": groups, "rule_i": rule_i,
            "refused_task_seeds": refused_n, "parity_agree_as_built": agree, "tasks": out_tasks, "merged": merged,
            "gate": {"e47_commit": E47_COMMIT, "where": E47_WHERE, "mean": GUARD_MEAN, "seed": GUARD_SEED}}


def text(res: dict) -> str:
    if res["void"]:
        return f"Try E48 (run commit {res['commit']}): VOID: " + "; ".join(res["why"]) + "\n"
    lines = [f"Try E48 (run commit {res['commit']}); the pre-run gate passed on all nine tasks (the controls under "
             f"main "
             f"within {res['gate']['mean']} / {res['gate']['seed']} of E47 at {res['gate']['e47_commit']}).",
             "Reading (i), fixed in row E48 (the row's arm is fgw_shape0x_series_cp):"]
    for a, r in res["rule_i"].items():
        if not r["can_price"]:
            lines.append(f"  {a}: no task-seed shares its truth rows with the floor copy; the read cannot price the "
                         "option")
        elif r["construction_void"]:
            lines.append(f"  {a}: the prefix mode differs from the floor's on {len(r['mode_differs_on'])} of "
                         f"{r['shared_task_seeds']} shared task-seeds ({', '.join(r['mode_differs_on'])}): the "
                         "construction is void")
        else:
            lines.append(f"  {a}: the prefix mode equals the floor's on all {r['shared_task_seeds']} shared task-seeds")
        lines.append(f"    refused task-seeds (the file as built ships): {res['refused_task_seeds'][a]} of "
                     f"{len(res['tasks']) * len(SEEDS)}")
    lines.append("Reading (ii), per task, arm minus control, paired by seed (primary: every seed, a refused seed "
                 "counts "
                 "as the control): board mean (sd) and occupancy skill mean (sd), main | prefix; secondary (refused "
                 "seeds left out) board mean, main | prefix:")
    for n, per in res["tasks"].items():
        for a, p in per["paired"].items():
            def f(x):
                return "  n/a  " if x["mean"] is None else f"{x['mean']:+.3f}" + (
                    f" ({x['sd']:.3f})" if x["sd"] is not None else "")
            lines.append(f"  {n:<5} {a:<22} board {f(p['main']['board'])} | {f(p['prefix']['board'])}   occupancy "
                         f"{f(p['main']['occupancy_skill'])} | {f(p['prefix']['occupancy_skill'])}   secondary "
                         f"{f(p['main']['secondary']['board'])} | {f(p['prefix']['secondary']['board'])} "
                         f"(left out {p['main']['secondary']['refused_left_out']})")
    lines.append("  every ranked metric's paired difference per task and build: reading.json, tasks.<task>.paired")
    lines.append("Reading (iii), per task, proper-minus-mirror Dice gap (mean over seeds):")
    for n, per in res["tasks"].items():
        lines.append(f"  {n:<5} " + ", ".join(f"{k} {v['gap']['mean']:+.4f}" for k, v in per["gap"].items()))
    lines.append("Reported, read by nothing: as built, the share of shared seeds where the control's pre-fix parity "
                 "agrees with the floor's (reading B.1):")
    for c, d in res["parity_agree_as_built"].items():
        lines.append(f"  {c}: " + ", ".join(f"{n} {'n/a' if v is None else v}" for n, v in d.items()))
    lines.append("Reported, read by nothing: nine-task merged paired means (tasks sharing a truth merged) with their "
                 "SE "
                 "over seeds:")
    for a, d in res["merged"].items():
        lines.append(f"  {a}: " + "; ".join(f"{b} board {v['board_mean']:+.3f} (SE {v['board_se']:.3f}), occupancy "
                                            f"{v['occupancy_mean']:+.4f} (SE {v['occupancy_se']:.4f})"
                                            for b, v in d.items()))
    lines.append("Note: on the option's own rows the shared-rows and parity checks re-test its gate; only the prefix "
                 "modes scored against the fold truth test the mechanism. Nothing here shows the live board's build, "
                 "metric seed, LAPACK, floor copy draw, row handling or whether the test counts share truth rows.")
    return "\n".join(lines) + "\n"


def load_runs(runs: Path) -> dict[str, list[dict]]:
    out = {}
    for name in nt.all_tasks():
        folder = runs / "run" / name
        recs = [json.loads(p.read_text()) for p in sorted(folder.glob("seed*.json"))] if folder.exists() else []
        out[name] = sorted(recs, key=lambda r: r["seed"])
    return out


def cmd_read(a) -> None:
    out = Path(a.out)
    files = [out / "reading.json", out / "reading.txt"]
    existing = [str(f) for f in files if f.exists()]
    if existing:
        sys.exit(f"refusing to overwrite: {', '.join(existing)}")
    _own_code()
    runs = Path(a.runs)
    why = [*gate_problems(runs, a.commit), *commit_problems(a.commit)]
    res = reading(load_runs(runs), why, a.commit)
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
