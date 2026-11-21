#!/usr/bin/env python
"""Row E50 of experiments/harness/TRIES.md (as specified; one try, placed by the main session before any run): the
copy-rows file, the investigator's Option A ("the second construction"), on the nine practice tasks. The captain's
order of 2026-11-21 ("reopen it [my ordeer]", docs/PLAN.md decision log 2026-11-21 (h)); specified by fresh analysts;
the main session is recused on this lever (decision log 2026-11-20 (i)). No number here comes from a test score.

Arm. The harness candidate ARM (fgw_copyrows_series, the series form: the exact-position file handed over is
fgw_shape0x_series; if the main session finds from that file's committed build record that it is the mixture form, it
sets HANDED_OVER_EXPRESSION to "mixture" in the row's commit, and the arm is fgw_copyrows). Only that arm counts. Per
task and scoring seed s: the rows, order, count and raw float32 coordinates of the bake floor copy
vec.evaluate.subsample(ref, 0.1, s), drawn from the very reference cells the harness's floor uses (vec.t2_build.
copy_rows_cloud, vec.harness._copy_rows). Controls: the floor, fgw_shape0x, fgw_shape0x_series.

Path. harness.run_task on each task of nine_tasks.all_tasks(), seeds 0-9, bake protocol, build seed BUILD_SEED_OFFSET
+ seed, under build "main" and under "prefix". A refusal of the arm (more dropped slots than the later bracket's
share, a count outside the task board's limits, the cloud guard) is caught and recorded by category; fgw_shape0x's
file of that seed stands in so the run goes on, and that task-seed is left out of every reading.

Steps (each refuses unless row E50 is in TRIES.md, the tree is clean and HEAD is the run commit):
  gate --task NAME --commit SHA  the floor scores 50, and fgw_shape0x and fgw_shape0x_series under main reproduce
                                 E47's committed scores (run at 9967cf1, committed in d5ee9b5) within 0.1 board points
                                 in the ten-seed mean and 0.5 in any seed; exit 3 on a failure (stop, report, count
                                 nothing). A gate file of this commit that exists is read again, not re-run.
  run --task NAME --commit SHA   the arm and controls under both builds, seed by seed (resumable: a seed whose file
                                 exists is skipped); refused unless all nine gate files pass at the same commit
  read --commit SHA              the reading into experiments/coords/results/e50_copy_rows/ (reading.json,
                                 reading.txt), never overwriting
The script prints only "<task> seed <s> done", "wrote <path>" and "gate <task> pass|fail". It reads no test board, test
file or test anchor, and prints no count.

Readings, fixed in row E50 before any output (applied by `reading`):
  (i)   on every task-seed where the arm builds, its d2_shape, occupancy_dice and scale_log_ratio equal the floor's
        exactly under main and under prefix; any inequality voids the construction and no test file is built;
  (ii)  refusals listed per task-seed; a task refused at every seed prices nothing;
  (iii) paired arm minus fgw_shape0x, per task and merged over the five truths, both builds: board score, each group,
        each metric, seeds ahead, SE, mean gain over gain tasks and mean loss over loss tasks, worst and best task; the
        same against fgw_shape0x_series as a secondary table;
  (iv)  per task, the arm's neighborhood_mmd skill and its shape group gated and ungated (where the assumed gate
        bites).
Outcome: if (i) holds, the test file is built under the captain's order with (iii) and (iv) beside it; (iii) and (iv)
decide nothing. Void: a missing task, seed or arm; mixed commits or an unclean tree; a gate miss; files differing
between the two builds' scoring; a refusal that differs between the builds; a test board.
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
from vec import t2_build as tb  # noqa: E402
from vec.adoption import SEEDS  # noqa: E402
from vec.boards import GROUPS, RANKED_METRICS  # noqa: E402
from vec.provenance import repo_state  # noqa: E402

ROW_IDS = ("E50",)
TRIES = REPO / "experiments" / "harness" / "TRIES.md"
HANDED_OVER_EXPRESSION = "series"  # the exact-position file handed over is fgw_shape0x_series (see the docstring)
ARMS_BY_EXPRESSION = {"series": "fgw_copyrows_series", "mixture": "fgw_copyrows"}
ARM = ARMS_BY_EXPRESSION[HANDED_OVER_EXPRESSION]
PRIMARY, SECONDARY = "fgw_shape0x", "fgw_shape0x_series"
CONTROLS = (PRIMARY, SECONDARY)
POSITION_TERMS = ("d2_shape", "occupancy_dice", "scale_log_ratio")
REPORTED_TERMS = ("sliced_wasserstein",)  # compared and reported where the panel holds it; deciding nothing
BUILDS = ("main", "prefix")
GROUP_SHAPE = "shape_growth"
GROUP_NAMES = tuple(name for name, _, _ in GROUPS["T2"])
METRICS = tuple(RANKED_METRICS["T2"])
NBHD = "neighborhood_mmd"
GATE_HALF = 0.5  # the assumed gate: shape * min(1, neighborhood skill / 0.5)
E47_DIR = REPO / "experiments" / "harness" / "results" / "e47_floor_exact"
E47_RUN_COMMIT, E47_FILES_COMMIT = "9967cf1", "d5ee9b5"
GATE_MEAN, GATE_SEED = 0.1, 0.5  # board points: a task's ten-seed mean, any one seed
GATE_ARMS = CONTROLS
FLOOR_SCORE = 50.0
BASES = ("35a21ac",)  # main when the copy-rows branch was cut (decision log 2026-11-21 (h))
GUARD_FLAG = "--copy-rows-test-ok"
SCRIPT = "experiments/coords/e50_copy_rows_read.py"
RUNS = REPO / "runs" / "e50"
RESULTS = REPO / "experiments" / "coords" / "results" / "e50_copy_rows"
REFUSALS = (("later bracket's share", "dropped_over_share"), ("cells outside [", "count_outside_limits"),
            ("cloud guard flags", "cloud_guard"), ("copy_rows:", "copy_rows_input"))


# ------------------------------------------------------------------------------------------------ checks
def _git(*args: str) -> subprocess.CompletedProcess:
    return subprocess.run(["git", "-C", str(REPO), *args], capture_output=True, text=True,
                          env={**os.environ, "GIT_OPTIONAL_LOCKS": "0"})


def has_row(path: Path | None = None) -> bool:
    rows = n12.tries_rows(TRIES if path is None else path)
    return any(r in rows for r in ROW_IDS)


def commit_problems(commit: str) -> list[str]:
    """Why this checkout may not run or read row E50 at `commit` (empty: it may)."""
    why = []
    head = _git("rev-parse", "HEAD").stdout.strip()
    if head != commit:
        why.append(f"HEAD is {head[:12]}, not the run commit {commit[:12]}")
    for base in BASES:
        if _git("merge-base", "--is-ancestor", base, "HEAD").returncode != 0:
            why.append(f"HEAD does not contain {base}")
    if GUARD_FLAG not in (REPO / "experiments" / "builds" / "build_test.py").read_text(encoding="utf-8"):
        why.append(f"the builder lacks the test-board guard {GUARD_FLAG}")
    if _git("diff", "--quiet", E47_FILES_COMMIT, "HEAD", "--", str(E47_DIR.relative_to(REPO))).returncode != 0:
        why.append(f"E47's committed results differ from {E47_FILES_COMMIT}'s")
    return why


def _own_code() -> None:
    src = (REPO / "src").resolve()
    for mod in (ev, harness, tb):
        if src not in Path(mod.__file__).resolve().parents:
            sys.exit(f"{mod.__name__} is loaded from {mod.__file__}, not from this checkout's src/")


def _refuse_unless_ready(commit: str) -> None:
    _own_code()
    if not has_row():
        sys.exit("E50 has no row in experiments/harness/TRIES.md: a row is committed before its run")
    if repo_state(REPO).unclean != ():
        sys.exit("the tree is not clean: row E50 runs on a committed tree")
    why = commit_problems(commit)
    if why:
        sys.exit("refused: " + "; ".join(why))


def _task(name: str):
    tasks = nt.all_tasks()
    if name not in tasks:
        sys.exit(f"unknown task {name}; known: {sorted(tasks)}")
    t = tasks[name]
    if getattr(t, "target_filter", None) is not None:
        sys.exit(f"{name}: a task with a target filter; E50 does not run on one")
    return t


def _is_test_board(board) -> bool:
    return str(board.split).startswith("test") or "test" in str(board.key).split(":")[-1]


def _sha_file(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def same_bits(C, R) -> bool:
    """Whether two coordinate arrays are the same as the scorer reads them: shape, row order and float32 bits."""
    A = np.ascontiguousarray(np.asarray(C)[:, :3], dtype=np.float32)
    B = np.ascontiguousarray(np.asarray(R)[:, :3], dtype=np.float32)
    return A.shape == B.shape and np.array_equal(A.view(np.uint32), B.view(np.uint32))


def refusal_category(message: str) -> str | None:
    """The category of an arm's refusal (no value kept), or None for an error that is not a refusal."""
    for needle, cat in REFUSALS:
        if needle in message:
            return cat
    return None


# ------------------------------------------------------------------------------------------------ the arm
def file_digest(stage) -> str:
    X = np.ascontiguousarray(stage.X, dtype=np.float32)
    C = np.ascontiguousarray(stage.coords, dtype=np.float32)
    return hashlib.sha256(np.asarray(X.shape, np.int64).tobytes() + X.tobytes() + C.tobytes()).hexdigest()


@contextlib.contextmanager
def captured():
    """While active, record what harness.score_seeds is handed (task, board, target, reference stage)."""
    seen: dict = {"stage": {}, "file": {}, "refused": {}}
    orig = harness.score_seeds

    def score_seeds(t, board, target, ref, inputs, candidates, *a, **kw):
        seen.update(task=t, board=board, target=target, ref=ref)
        return orig(t, board, target, ref, inputs, candidates, *a, **kw)

    harness.score_seeds = score_seeds
    try:
        yield seen
    finally:
        harness.score_seeds = orig


def _keep(seen, name: str, s: int, stage) -> None:
    seen["stage"][(name, s)] = stage
    seen["file"][(name, s)] = file_digest(stage)


def observed(base, seen):
    """A control unchanged; its file is kept (fgw_shape0x's also stands in for a refused arm)."""
    def build(ctx):
        pred = base.build(ctx)
        _keep(seen, base.name, ctx.seed - harness.BUILD_SEED_OFFSET, pred)
        return pred

    return dataclasses.replace(base, build=build)


def guarded(base, seen):
    """The arm: its file, or on a refusal fgw_shape0x's file of the seed, with the refusal's category recorded."""
    def build(ctx):
        s = ctx.seed - harness.BUILD_SEED_OFFSET
        if _is_test_board(seen["board"]):
            raise RuntimeError("E50 met a test board")
        try:
            pred, why = base.build(ctx), None
        except ValueError as e:
            why = refusal_category(str(e))
            if why is None:
                raise
            pred = seen["stage"][(PRIMARY, s)]
        seen["refused"][(base.name, s)] = why
        _keep(seen, base.name, s, pred)
        return pred

    return dataclasses.replace(base, build=build)


def candidates(seen) -> list:
    return [harness.CANDIDATES[harness.FLOOR], *(observed(harness.CANDIDATES[c], seen) for c in CONTROLS),
            guarded(harness.CANDIDATES[ARM], seen)]


def _mode(coords, truth, s: int) -> str:
    return ev.handedness(coords, truth.coords, seed=s, live_n_true=truth.n)["prefix_mode"]


def one_seed(t, s: int, data: Path) -> dict:
    """Both builds of one task-seed: run_task's results and E50's per-seed facts."""
    out: dict = {"builds": {}, "facts": {}}
    names = (*CONTROLS, ARM)
    files, refused = {}, {}
    for build in BUILDS:
        with captured() as seen:
            res = harness.run_task(t, candidates(seen), seeds=[s], data_dir=data, build=build, log=None)
        board = seen["board"]
        if _is_test_board(board):
            raise RuntimeError("E50 met a test board")
        out["builds"][build] = res
        files[build] = {n: seen["file"][(n, s)] for n in names}
        refused[build] = seen["refused"].get((ARM, s))
        if build != BUILDS[0]:
            continue
        floor = ev.subsample(seen["ref"], ev.SUBSAMPLE_FRAC, s)
        truth = ev.working_copies(seen["target"], seen["ref"], s, ev.SUBSAMPLE_FRAC, harness.PROTOCOL)[0]
        arm = seen["stage"][(ARM, s)]
        built = refused[build] is None
        out["facts"] = {"board_split": str(board.split), "refused": refused[build],
                        "coords_equal_floor": bool(built and same_bits(arm.coords, floor.coords)),
                        "prefix_mode": {"arm": _mode(arm.coords, truth, s) if built else None,
                                        "floor": _mode(floor.coords, truth, s)}}
    out["facts"]["files_identical_across_builds"] = all(files[BUILDS[0]][n] == files[b][n]
                                                        for b in BUILDS[1:] for n in names)
    out["facts"]["refusals_identical_across_builds"] = all(refused[BUILDS[0]] == refused[b] for b in BUILDS[1:])
    return out


# ------------------------------------------------------------------------------------------------ gate
def _gate_path(runs: Path, name: str) -> Path:
    return Path(runs) / "gate" / f"{name}.json"


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


def gate_check(name: str, key: str, run_scores: dict, old: dict | None = None,
               why_old: list[str] | None = None) -> dict:
    """The pre-run gate on one task: the floor at 50 in every seed, and each control within GATE_MEAN of E47's ten-seed
    mean and GATE_SEED of E47 in every seed. `old`/`why_old` default to the committed E47 scores of `key`."""
    if old is None and why_old is None:
        old, why_old = e47_scores(key)
    why = list(why_old or [])
    floor = run_scores.get(harness.FLOOR)
    if floor is None or len(floor) != len(SEEDS) or not np.all(np.abs(np.asarray(floor) - FLOOR_SCORE) <= 1e-9):
        why.append(f"{key}: the floor does not score {FLOOR_SCORE:g} in every seed")
    worst = {"mean": None, "seed": None}
    if old is not None:
        worst = {"mean": 0.0, "seed": 0.0}
        for arm in GATE_ARMS:
            o, n = old.get(arm), run_scores.get(arm)
            if o is None or n is None or np.shape(o) != np.shape(n):
                why.append(f"{key}: {arm} missing or of another length")
                continue
            o, n = np.asarray(o, float), np.asarray(n, float)
            dm, ds = abs(float(n.mean() - o.mean())), float(np.abs(n - o).max())
            worst = {"mean": max(worst["mean"], dm), "seed": max(worst["seed"], ds)}
            if not dm <= GATE_MEAN:
                why.append(f"{key}: {arm}'s ten-seed mean is {dm:.4f} off E47's (tolerance {GATE_MEAN})")
            if not ds <= GATE_SEED:
                why.append(f"{key}: {arm} is {ds:.4f} off E47's in a seed (tolerance {GATE_SEED})")
    return {"task": name, "key": key, "e47": {"run_commit": E47_RUN_COMMIT, "files_commit": E47_FILES_COMMIT,
                                               "dir": str(E47_DIR.relative_to(REPO))},
            "tolerance": {"mean": GATE_MEAN, "seed": GATE_SEED}, "worst": worst, "why": why, "pass": not why}


def cmd_gate(a) -> None:
    _refuse_unless_ready(a.commit)
    t = _task(a.task)
    path = _gate_path(Path(a.runs), a.task)
    if path.exists():  # resumable: a gate file of this commit is read again, not re-run
        rec = json.loads(path.read_text())
        if rec.get("commit") != a.commit:
            sys.exit(f"{path} exists from another commit")
    else:
        res = harness.run_task(t, [harness.CANDIDATES[harness.FLOOR], *(harness.CANDIDATES[c] for c in GATE_ARMS)],
                               seeds=SEEDS, data_dir=Path(a.data), build="main", log=None)
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


def cmd_run(a) -> None:
    _refuse_unless_ready(a.commit)
    t = _task(a.task)
    runs = Path(a.runs)
    why = gate_problems(runs, a.commit)
    if why:
        sys.exit("refused (the pre-run gate): " + "; ".join(why))
    for c in (*CONTROLS, ARM):
        if not harness.CANDIDATES[c].applies(t):
            sys.exit(f"{a.task}: {c} does not apply")
    folder = runs / "run" / a.task
    folder.mkdir(parents=True, exist_ok=True)
    for s in SEEDS:
        path = folder / f"seed{s}.json"
        if path.exists():
            continue
        res = one_seed(t, s, Path(a.data))
        rec = {"task": a.task, "key": t.key, "seed": s, "run_commit": a.commit, "script": SCRIPT, "arm": ARM, **res}
        tmp = path.with_suffix(".tmp")
        tmp.write_text(json.dumps(rec, indent=1, default=float))
        os.replace(tmp, path)
        print(f"{a.task} seed {s} done", flush=True)
    print(f"wrote {folder}", flush=True)


# ------------------------------------------------------------------------------------------------ read
def check_run(recs: dict[str, list[dict]], why_in: list[str], commit: str) -> list[str]:
    """The void conditions on the per-seed records (task -> its records in seed order)."""
    why = list(why_in)
    names = (harness.FLOOR, *CONTROLS, ARM)
    for task in nt.all_tasks():
        rs = recs.get(task, [])
        if [r["seed"] for r in rs] != list(SEEDS):
            why.append(f"{task}: seeds are not {list(SEEDS)}")
            continue
        for r in rs:
            tag = f"{task} seed {r['seed']}"
            if r.get("run_commit") != commit:
                why.append(f"{tag}: another run commit")
            if r.get("arm") != ARM:
                why.append(f"{tag}: run with another arm")
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
                if str(res["board"].get("split", "")).startswith("test"):
                    why.append(f"{tag}: a test board")
                row = res["rows"][0]
                if any(n not in row["scores"] or n not in row["skills"] or n not in row["metrics"] for n in names):
                    why.append(f"{tag}: build {b} lacks an arm")
            f = r["facts"]
            if str(f.get("board_split", "")).startswith("test"):
                why.append(f"{tag}: a test board")
            if not f.get("files_identical_across_builds"):
                why.append(f"{tag}: a file differs between the builds")
            if not f.get("refusals_identical_across_builds"):
                why.append(f"{tag}: a refusal differs between the builds")
    return why


def _row(r: dict, build: str) -> dict:
    return r["builds"][build]["rows"][0]


def built(r: dict) -> bool:
    return r["facts"].get("refused") is None


def reading_i(recs: dict[str, list[dict]]) -> dict:
    """Reading (i): on every task-seed where the arm builds, its position terms equal the floor's exactly under both
    builds (and its coordinates are the floor copy's, bit for bit). Reported beside it, deciding nothing: the terms
    in REPORTED_TERMS where the panel holds them, and the arm's and the floor's prefix modes."""
    unequal, reported, modes, n_built = [], [], [], 0
    for task, rs in recs.items():
        for r in rs:
            if not built(r):
                continue
            n_built += 1
            if not r["facts"].get("coords_equal_floor"):
                unequal.append({"task": task, "seed": r["seed"], "what": "coordinates"})
            pm = r["facts"].get("prefix_mode", {})
            if pm.get("arm") != pm.get("floor"):
                modes.append({"task": task, "seed": r["seed"]})
            for b in BUILDS:
                m = _row(r, b)["metrics"]
                for term in POSITION_TERMS:
                    if m[ARM].get(term) != m[harness.FLOOR].get(term):
                        unequal.append({"task": task, "seed": r["seed"], "build": b, "what": term})
                for term in REPORTED_TERMS:
                    if term in m[ARM] and term in m[harness.FLOOR] and m[ARM][term] != m[harness.FLOOR][term]:
                        reported.append({"task": task, "seed": r["seed"], "build": b, "what": term})
    return {"holds": bool(n_built > 0 and not unequal), "built_task_seeds": n_built, "unequal": unequal,
            "reported_unequal": reported, "prefix_mode_differs": modes}


def reading_ii(recs: dict[str, list[dict]]) -> dict:
    per = {task: {str(r["seed"]): r["facts"]["refused"] for r in rs if not built(r)} for task, rs in recs.items()}
    every = sorted(task for task, rs in recs.items() if rs and not any(built(r) for r in rs))
    return {"refused": {k: v for k, v in per.items() if v}, "refused_at_every_seed": every}


def value(row: dict, name: str, what: str) -> float:
    """One quantity of one file in one row: 'board', 'group:<g>' (gated, as the row holds it), 'ungated:<g>' or
    'skill:<m>'."""
    if what == "board":
        return float(row["scores"][name])
    if what.startswith("group:"):
        return float(row["groups"][name][what[6:]])
    if what.startswith("ungated:"):
        return float(ev.aggregate("T2", row["skills"][name], gate="none")[1][what[8:]])
    if what.startswith("skill:"):
        return float(row["skills"][name][what[6:]])
    raise ValueError(what)


def _diffs(rs: list[dict], a: str, b: str, what: str, build: str) -> dict[int, float]:
    """seed -> a minus b, on the seeds where the arm builds."""
    return {r["seed"]: value(_row(r, build), a, what) - value(_row(r, build), b, what) for r in rs if built(r)}


def paired(recs: dict[str, list[dict]], groups: dict[str, list[str]], a: str, b: str, what: str, build: str) -> dict:
    """Reading (iii) for one quantity: per task (mean, seeds ahead, seeds read), merged over the truths (a truth's
    value the mean over its tasks that have a built seed), D_s per seed (the same merge on that seed's built tasks),
    its SE, the mean gain over gain tasks, the mean loss over loss tasks, the worst and the best task."""
    per_task = {}
    for task, rs in recs.items():
        d = _diffs(rs, a, b, what, build)
        if d:
            v = np.array(list(d.values()), float)
            per_task[task] = {"mean": float(v.mean()), "ahead": int((v > 0).sum()), "behind": int((v < 0).sum()),
                              "seeds": len(v), "per_seed": {str(k): float(x) for k, x in sorted(d.items())}}
    truth_vals = [np.mean([per_task[n]["mean"] for n in names if n in per_task])
                  for names in groups.values() if any(n in per_task for n in names)]
    D = []
    for s in SEEDS:
        vals = [np.mean([per_task[n]["per_seed"][str(s)] for n in names
                         if n in per_task and str(s) in per_task[n]["per_seed"]])
                for names in groups.values()
                if any(n in per_task and str(s) in per_task[n]["per_seed"] for n in names)]
        if vals:
            D.append(float(np.mean(vals)))
    D = np.array(D, float)
    means = {t: v["mean"] for t, v in per_task.items()}
    gains = [v for v in means.values() if v > 0]
    losses = [v for v in means.values() if v < 0]
    return {"per_task": {t: {k: v[k] for k in ("mean", "ahead", "behind", "seeds")} for t, v in per_task.items()},
            "merged": float(np.mean(truth_vals)) if truth_vals else None,
            "D_s": D.tolist(), "se": float(D.std(ddof=1) / np.sqrt(len(D))) if len(D) > 1 else None,
            "seeds_ahead": int((D > 0).sum()), "seeds_behind": int((D < 0).sum()),
            "gain_tasks": len(gains), "mean_gain": float(np.mean(gains)) if gains else None,
            "loss_tasks": len(losses), "mean_loss": float(np.mean(losses)) if losses else None,
            "worst_task": min(means, key=means.get) if means else None,
            "best_task": max(means, key=means.get) if means else None}


def reading_iii(recs: dict[str, list[dict]], groups: dict[str, list[str]], against: str) -> dict:
    whats = ("board", *(f"group:{g}" for g in GROUP_NAMES), *(f"skill:{m}" for m in METRICS))
    return {b: {w: paired(recs, groups, ARM, against, w, b) for w in whats} for b in BUILDS}


def reading_iv(recs: dict[str, list[dict]]) -> dict:
    out = {}
    for task, rs in recs.items():
        rows = [r for r in rs if built(r)]
        if not rows:
            out[task] = None
            continue
        out[task] = {}
        for b in BUILDS:
            nb = [value(_row(r, b), ARM, f"skill:{NBHD}") for r in rows]
            out[task][b] = {"neighborhood_mmd_skill": float(np.mean(nb)),
                            "seeds_gate_bites": int(sum(x < GATE_HALF for x in nb)),
                            "shape_gated": float(np.mean([value(_row(r, b), ARM, f"group:{GROUP_SHAPE}")
                                                          for r in rows])),
                            "shape_ungated": float(np.mean([value(_row(r, b), ARM, f"ungated:{GROUP_SHAPE}")
                                                            for r in rows])),
                            "floor_shape": float(np.mean([value(_row(r, b), harness.FLOOR, f"group:{GROUP_SHAPE}")
                                                          for r in rows]))}
    return out


def reading(recs: dict[str, list[dict]], why_in: list[str], commit: str) -> dict:
    why = check_run(recs, why_in, commit)
    if why:
        return {"void": True, "why": why, "commit": commit, "arm": ARM}
    groups = e47.truths(nt.all_tasks())
    i = reading_i(recs)
    outcome = ("build: (i) holds; the test file is built under the captain's order, with (iii) and (iv) beside it "
               "(they decide nothing)" if i["holds"] else
               "no test file: (i) fails, which voids the construction")
    return {"void": False, "why": [], "commit": commit, "arm": ARM, "handed_over_expression": HANDED_OVER_EXPRESSION,
            "truths": groups, "i": i, "ii": reading_ii(recs),
            "iii": {"primary": {"against": PRIMARY, **reading_iii(recs, groups, PRIMARY)},
                    "secondary": {"against": SECONDARY, **reading_iii(recs, groups, SECONDARY)}},
            "iv": reading_iv(recs), "build_test_file": bool(i["holds"]), "outcome": outcome,
            "gate": {"e47_run_commit": E47_RUN_COMMIT, "e47_files_commit": E47_FILES_COMMIT,
                     "tolerance": {"mean": GATE_MEAN, "seed": GATE_SEED}}}


def _fmt(x, f="+.3f") -> str:
    return "n/a" if x is None else format(x, f)


def text(res: dict) -> str:
    if res["void"]:
        return f"Row E50 (run commit {res['commit']}): VOID: " + "; ".join(res["why"]) + "\n"
    i = res["i"]
    L = [f"Row E50 (run commit {res['commit']}), arm {res['arm']}; the pre-run gate passed on all nine tasks.",
         f"(i) position terms equal to the floor's on every built task-seed, both builds: "
         f"{'HOLDS' if i['holds'] else 'FAILS'}" + ("" if i["holds"] else f" ({len(i['unequal'])} inequalities)"),
         f"    reported: {len(i['reported_unequal'])} unequal reported terms; prefix mode differs from the floor's on "
         f"{len(i['prefix_mode_differs'])} task-seeds",
         "(ii) refusals per task-seed: " + (json.dumps(res["ii"]["refused"]) if res["ii"]["refused"] else "none")
         + "; refused at every seed: " + (", ".join(res["ii"]["refused_at_every_seed"]) or "none")]
    for tier in ("primary", "secondary"):
        t = res["iii"][tier]
        L.append(f"(iii) {tier}: arm minus {t['against']}, merged over five truths (mean, SE, seeds ahead):")
        for b in BUILDS:
            for w, v in t[b].items():
                L.append(f"    {b} {w}: {_fmt(v['merged'], '+.4f')} (SE {_fmt(v['se'], '.4f')}, {v['seeds_ahead']} "
                         f"ahead); gain tasks {v['gain_tasks']} mean {_fmt(v['mean_gain'], '+.4f')}; loss tasks "
                         f"{v['loss_tasks']} mean {_fmt(v['mean_loss'], '+.4f')}; worst {v['worst_task']}, best "
                         f"{v['best_task']}")
    L.append("(iv) per task (main): neighborhood_mmd skill, seeds where the assumed gate bites, shape gated / ungated:")
    for task, v in res["iv"].items():
        if v is None:
            L.append(f"    {task}: refused at every seed")
            continue
        m = v["main"]
        L.append(f"    {task}: {m['neighborhood_mmd_skill']:.3f}, {m['seeds_gate_bites']}, {m['shape_gated']:.4f} / "
                 f"{m['shape_ungated']:.4f}")
    L.append(f"Outcome: {res['outcome']}.")
    L.append("Nothing here shows that the live server draws its floor as subsample(E7.5 as released, 0.1, 0), scores "
             "rows as held, or baked its anchors with the scorer build it scores uploads with; nothing about E7.75.")
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
            "t2_build_sha256": _sha_file(Path(tb.__file__)), "harness_sha256": _sha_file(Path(harness.__file__))}


def cmd_read(a) -> None:
    out = Path(a.out)
    files = [out / "reading.json", out / "reading.txt"]
    existing = [str(f) for f in files if f.exists()]
    if existing:
        sys.exit(f"refusing to overwrite: {', '.join(existing)}")
    _own_code()
    if not has_row():
        sys.exit("E50 has no row in experiments/harness/TRIES.md")
    runs = Path(a.runs)
    why = [*gate_problems(runs, a.commit), *commit_problems(a.commit)]
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
