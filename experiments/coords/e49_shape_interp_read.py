#!/usr/bin/env python
"""Row E49 of experiments/harness/TRIES.md (as drafted; E49a to E49g, placed by the main session before any run): shape
interpolation between the released brackets on the embryo's exact-position file, on the nine practice tasks.

Arms (harness candidates, vec.harness.SHAPE_CANDIDATES; each is fgw_shape0x with one form of vec.shape_interp applied
to its exact-position cloud before the size rule):
  E49a fgw_shape0x_sr  F1 stretch_ratio        E49d fgw_shape0x_sr_anti  its control (the step mirrored)
  E49b fgw_shape0x_sg  F2 stretch_registered   E49e fgw_shape0x_sg_anti  its control (the covariance step mirrored)
                                               E49f fgw_shape0x_sg_rot   its control (same spectrum, random turn)
  E49c fgw_shape0x_ov  F3 overlay              E49g fgw_shape0x_ov_rot   its control (random turn for R)
C0 = fgw_shape0x as built and the floor are rescored and not counted.

Path. harness.run_task on each task of nine_tasks.all_tasks() (half, late, new1, new2, e75, e725, x1, x2, x3), seeds
0-9, bake protocol, build seed BUILD_SEED_OFFSET + seed, once with build "main" (primary) and once with "prefix"
(reported, deciding nothing). Nothing of the task cut, the seeds or the read guard is re-implemented: the arms and
controls are wrapped harness candidates; a refusal of a form (vec.shape_interp.ShapeRefused, or the harness's cloud
guard) is caught, recorded, and C0's file of that seed stands in its place, so the run goes on and that arm and its
controls are unreadable.

Steps (each refuses unless row E49 is in TRIES.md, the tree is clean and HEAD is the run commit):
  describe --commit SHA          step (0): target-free. E7.5 and E8.0 only (the test pair), with the validation
                                 board's panel and count rule, build seeds 1,000,001 to 1,000,010: F1's seed-free blocks
                                 and factors, F2's eigen-shares, guard and SDs, F3's thinnest SD; the seed-free blocks
                                 of every practice task's brackets (no target is read). Fixes whether F1 and its
                                 control are dropped: if F1 is the identity at the test pair, or its blocks there equal
                                 no test-like task's blocks. Prints no count.
  gate --task NAME --commit SHA  step (1): the floor and C0 under main, ten seeds, against E47's committed scores
                                 (experiments/harness/results/e47_floor_exact, run at 9967cf1, committed in d5ee9b5)
                                 within 0.01 board points per seed and per ten-seed mean; exit 3 on a failure
  run --task NAME --commit SHA   step (2): the arms and controls (less any dropped at step 0), both builds, seed by seed
                                 (resumable: a seed whose file exists is skipped); refused unless describe exists and
                                 all nine gate files pass at the same commit
  read --commit SHA              the reading into experiments/coords/results/e49_shape_interp/ (reading.json,
                                 reading.txt), never overwriting, with the reading code's provenance
Nothing is selected or tuned on the output. The script prints no count; it reads no test board, test file or test
anchor.

Reading, fixed in row E49 before any output (applied by `reading`). D_s(X) for a difference X (form minus reference)
at seed s is the mean over the five truths of each truth's task-mean difference at s (E8.0 = half, late; E8.25_late =
new1, x1; E8.5 = new2, x2, x3; E7.5 = e75; E7.25 = e725); the merged mean is mean_s D_s and its SE SD(D_s)/sqrt(10).
The test-like series is the same over half, late, new1, new2 (three truths). Under build main, form X passes only if:
  P1  merged board mean against C0 above 0.3, with D_s above 0 in at least 8 of 10 seeds;
  P2  merged gated shape_growth mean against C0 above 0;
  P3  test-like means against C0 above 0 on the board and on gated shape_growth;
  P4  behind on no other group (expression_change, cell_state, local_organisation): a group is behind if its merged
      mean is below 0 with at least 8 of 10 D_s below 0;
  P5  against each of its controls, the ungated shape mean S = (d2_shape skill + occupancy_dice skill) / 2 has a
      merged mean above 0 with at least 8 of 10 D_s above 0;
  P6  no refusal of X or of its controls on any task-seed (a refusal makes them unreadable).
A mixed outcome is a no-pass. Candidates: none if nothing passes; otherwise at most two, the passing forms with the
highest merged board mean (ties: the higher test-like board mean, then F2, F1, F3). F1 is a candidate only if its
blocks at the test pair (step 0) equal the blocks of at least one test-like task in this run.

Void (the run is read as nothing; no task or seed is re-run selectively): a missing task, seed, build or arm; a result
from another commit or an unclean tree; describe missing or from another commit; a gate file missing or failing; an
arm or control whose count or expression matrix (in row order) differs from C0's; an expression-only skill (de_score,
de_direction, mmd_u, variogram) off C0's by more than 1e-9; scale_log_ratio's skill off C0's by more than 1e-6 on a
task-seed where neither file was jittered (1e-6, not 1e-9: the file stores float32 coordinates, so the RMS the
scorer reads moves by rounding between two position sets); a file not identical across the two builds; a refusal that
differs between the builds; a test board.
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
from vec.provenance import repo_state  # noqa: E402

ROW_IDS = ("E49a", "E49")  # the row's first word: "E49a to E49g" as drafted, or "E49" as the main session numbers it
TRIES = REPO / "experiments" / "harness" / "TRIES.md"
C0 = "fgw_shape0x"
ARMS = {"fgw_shape0x_sr": ("fgw_shape0x_sr_anti",),
        "fgw_shape0x_sg": ("fgw_shape0x_sg_anti", "fgw_shape0x_sg_rot"),
        "fgw_shape0x_ov": ("fgw_shape0x_ov_rot",)}
TRY_IDS = {"fgw_shape0x_sr": "E49a", "fgw_shape0x_sg": "E49b", "fgw_shape0x_ov": "E49c",
           "fgw_shape0x_sr_anti": "E49d", "fgw_shape0x_sg_anti": "E49e", "fgw_shape0x_sg_rot": "E49f",
           "fgw_shape0x_ov_rot": "E49g"}
FORM_OF = dict(harness.SHAPE_CANDIDATES)
F1 = "fgw_shape0x_sr"
TIE_ORDER = ("fgw_shape0x_sg", "fgw_shape0x_sr", "fgw_shape0x_ov")
MAX_CANDIDATES = 2
BUILDS = ("main", "prefix")
TEST_LIKE = ("half", "late", "new1", "new2")
GROUP_SHAPE = "shape_growth"
OTHER_GROUPS = ("expression_change", "cell_state", "local_organisation")
SHAPE_TERMS = ("d2_shape", "occupancy_dice")
EXPR_METRICS = e47.EXPR_METRICS
EXPR_TOL = 1e-9
SCALE_TOL = 1e-6
BAR, AHEAD = 0.3, 8
E47_DIR = REPO / "experiments" / "harness" / "results" / "e47_floor_exact"
E47_RUN_COMMIT, E47_FILES_COMMIT = "9967cf1", "d5ee9b5"
GATE_TOL = 0.01  # board points, per seed and per ten-seed mean (nine_tasks.TOL)
GATE_ARMS = (harness.FLOOR, C0)
BASES = ("9cd4b6e",)  # the E48 branch this one merges (row order option, E48's script)
GUARD_FLAG = "--shape-interp-test-ok"
DESCRIBE_SEEDS = tuple(harness.BUILD_SEED_OFFSET + 1 + i for i in range(10))
DESCRIBE_BOARD = "T2:embryo:val_interp"
DESCRIBE_PAIR = {"before": "T2/E7.5.h5ad", "after": "T2/E8.0.h5ad"}
SCRIPT = "experiments/coords/e49_shape_interp_read.py"
RUNS = REPO / "runs" / "e49"
RESULTS = REPO / "experiments" / "coords" / "results" / "e49_shape_interp"
REFUSAL = "cloud guard flags"


# ------------------------------------------------------------------------------------------------ checks
def _git(*args: str) -> subprocess.CompletedProcess:
    return subprocess.run(["git", "-C", str(REPO), *args], capture_output=True, text=True,
                          env={**os.environ, "GIT_OPTIONAL_LOCKS": "0"})


def has_row(path: Path | None = None) -> bool:
    path = TRIES if path is None else path
    rows = n12.tries_rows(path)
    return any(r in rows for r in ROW_IDS)


def commit_problems(commit: str) -> list[str]:
    """Why this checkout may not run or read row E49 at `commit` (empty: it may)."""
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
    for mod in (ev, harness, tb, si):
        if src not in Path(mod.__file__).resolve().parents:
            sys.exit(f"{mod.__name__} is loaded from {mod.__file__}, not from this checkout's src/")


def _refuse_unless_ready(commit: str) -> None:
    _own_code()
    if not has_row():
        sys.exit("E49 has no row in experiments/harness/TRIES.md: a row is committed before its run")
    if repo_state(REPO).unclean != ():
        sys.exit("the tree is not clean: row E49 runs on a committed tree")
    why = commit_problems(commit)
    if why:
        sys.exit("refused: " + "; ".join(why))


def _task(name: str):
    tasks = nt.all_tasks()
    if name not in tasks:
        sys.exit(f"unknown task {name}; known: {sorted(tasks)}")
    t = tasks[name]
    if getattr(t, "target_filter", None) is not None:
        sys.exit(f"{name}: a task with a target filter; E49 does not run on one")
    return t


def _is_test_board(board) -> bool:
    return str(board.split).startswith("test") or "test" in str(board.key).split(":")[-1]


def _sha_file(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


# ------------------------------------------------------------------------------------------------ step 0
def _frame(C) -> np.ndarray:
    C = np.asarray(C, dtype=np.float64)[:, :3]
    return (C - C.mean(0)) / tb.rms_radius(C)


def blocks_of_pair(before, after) -> str:
    """F1's seed-free blocks for two bracket clouds (raw coordinates), as text ('1|2,3')."""
    return si.blocks_text(si.bracket_blocks(_frame(before), _frame(after)))


def drop_decision(test_blocks: str, test_like_blocks: dict[str, str]) -> dict:
    """Step 0's rule: F1 and its control are dropped if F1 is the identity at the test pair (one block of all three
    axes) or if its test-pair blocks equal no test-like task's blocks."""
    identity = test_blocks == "1,2,3"
    unmatched = test_blocks not in set(test_like_blocks.values())
    return {"test_pair_blocks": test_blocks, "test_like_blocks": dict(test_like_blocks), "f1_identity": identity,
            "f1_unmatched": unmatched, "drop_f1": bool(identity or unmatched)}


def dropped(describe: dict) -> tuple[str, ...]:
    return (F1, *ARMS[F1]) if describe["drop"]["drop_f1"] else ()


def cmd_describe(a) -> None:
    _refuse_unless_ready(a.commit)
    from vec.boards import load_boards
    from vec.io import COORD_KEY, dense
    from vec.prep import load_stage

    out = Path(a.runs) / "describe.json"
    if out.exists():
        sys.exit(f"{out} exists")
    board = load_boards()[DESCRIBE_BOARD]
    if _is_test_board(board):
        sys.exit("describe reads a validation board's panel and count rule only")
    genes, data = list(board.genes), Path(a.data)
    st = {r: load_stage(data / f, genes=genes, drop_duplicates=True) for r, f in DESCRIBE_PAIR.items()}
    Cb, Ca = (np.asarray(st[r].obsm[COORD_KEY], float)[:, :3] for r in ("before", "after"))
    Eb, Ea = (dense(st[r][:, genes].X).astype(np.float32) for r in ("before", "after"))
    n = harness.count_for(board, st["before"].n_obs)
    per_seed = {}
    for seed in DESCRIBE_SEEDS:
        per_seed[str(seed)] = {}
        for form in si.FORMS:
            try:
                res = tb.interpolate_clouds(Cb, Ca, Eb, Ea, 0.5, n, seed, target_rms=1.0, t_pos=0.0,
                                            positions=tb.POSITIONS[1], shape_form=form)
                per_seed[str(seed)][form] = {k: v for k, v in res.record.items() if k.startswith("shape_interp")}
            except si.ShapeRefused as e:
                per_seed[str(seed)][form] = {"refused": str(e)}
    test_blocks = blocks_of_pair(Cb, Ca)
    task_blocks = {}
    for name, t in nt.all_tasks().items():  # the brackets only: no target is read
        c = [np.asarray(load_stage(data / t.inputs[r], genes=genes, drop_duplicates=True).obsm[COORD_KEY], float)
             for r in ("before", "after")]
        task_blocks[name] = blocks_of_pair(*c)
    rec = {"commit": a.commit, "script": SCRIPT, "board_panel_and_count_rule": DESCRIBE_BOARD,
           "pair": DESCRIBE_PAIR, "t": 0.5, "seeds": list(DESCRIBE_SEEDS), "per_seed": per_seed,
           "task_blocks": task_blocks,
           "drop": drop_decision(test_blocks, {k: task_blocks[k] for k in TEST_LIKE})}
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(rec, indent=1, default=float))
    print(f"wrote {out}", flush=True)


def load_describe(runs: Path, commit: str) -> tuple[dict | None, list[str]]:
    p = Path(runs) / "describe.json"
    if not p.exists():
        return None, ["no describe.json (step 0)"]
    d = json.loads(p.read_text())
    return d, ([] if d.get("commit") == commit else ["describe.json is from another commit"])


# ------------------------------------------------------------------------------------------------ the arms
def x_digest(stage) -> str:
    """The expression matrix in row order (float32 bytes) and its count, as one digest."""
    X = np.ascontiguousarray(stage.X, dtype=np.float32)
    return hashlib.sha256(np.asarray(X.shape, np.int64).tobytes() + X.tobytes()).hexdigest()


def file_digest(stage) -> str:
    X = np.ascontiguousarray(stage.X, dtype=np.float32)
    C = np.ascontiguousarray(stage.coords, dtype=np.float32)
    return hashlib.sha256(X.tobytes() + C.tobytes()).hexdigest()


@contextlib.contextmanager
def captured():
    """While active, record what harness.score_seeds is handed (task, board, target, reference stage) and every
    interpolate_clouds record by (shape form, build seed): flat values only, no count."""
    seen: dict = {"stage": {}, "x": {}, "file": {}, "refused": {}, "records": {}}
    orig_ss, orig_ic = harness.score_seeds, tb.interpolate_clouds

    def score_seeds(t, board, target, ref, inputs, candidates, *a, **kw):
        seen.update(task=t, board=board, target=target, ref=ref)
        return orig_ss(t, board, target, ref, inputs, candidates, *a, **kw)

    def interpolate_clouds(*a, **kw):
        res = orig_ic(*a, **kw)
        seed = int(kw["seed"] if "seed" in kw else a[6])
        r = res.record
        seen["records"][(kw.get("shape_form"), seed)] = {
            "n_jittered": int(r.get("n_jittered", 0)),
            "fallback": bool((r.get("floor_exact") or {}).get("fallback", False)),
            **{k: v for k, v in r.items() if k.startswith("shape_interp") and not isinstance(v, (dict, list))}}
        return res

    harness.score_seeds, tb.interpolate_clouds = score_seeds, interpolate_clouds
    try:
        yield seen
    finally:
        harness.score_seeds, tb.interpolate_clouds = orig_ss, orig_ic


def _keep(seen, name: str, s: int, stage) -> None:
    seen["stage"][(name, s)] = stage
    seen["x"][(name, s)] = x_digest(stage)
    seen["file"][(name, s)] = file_digest(stage)


def observed(base, seen):
    """C0 unchanged; its file is kept for the checks and as the stand-in for a refused form."""
    def build(ctx):
        pred = base.build(ctx)
        _keep(seen, base.name, ctx.seed - harness.BUILD_SEED_OFFSET, pred)
        return pred

    return dataclasses.replace(base, build=build)


def guarded(base, seen):
    """An arm or control: its file, or on a refusal (ShapeRefused, the cloud guard) C0's file of the seed, recorded."""
    def build(ctx):
        s = ctx.seed - harness.BUILD_SEED_OFFSET
        if _is_test_board(seen["board"]):
            raise RuntimeError("E49 met a test board")
        try:
            pred, why = base.build(ctx), None
        except si.ShapeRefused as e:
            pred, why = None, f"shape: {e}"
        except ValueError as e:
            if REFUSAL not in str(e):
                raise
            pred, why = None, "cloud guard"
        if pred is None:
            pred = seen["stage"][(C0, s)]
        seen["refused"][(base.name, s)] = why
        _keep(seen, base.name, s, pred)
        return pred

    return dataclasses.replace(base, build=build)


def run_names(drop: tuple[str, ...] = ()) -> tuple[str, ...]:
    return tuple(n for a, cs in ARMS.items() for n in (a, *cs) if n not in drop)


def candidates(seen, drop: tuple[str, ...] = ()) -> list:
    return [harness.CANDIDATES[harness.FLOOR], observed(harness.CANDIDATES[C0], seen),
            *(guarded(harness.CANDIDATES[n], seen) for n in run_names(drop))]


def _hand(coords, truth, s: int) -> dict:
    h = ev.handedness(coords, truth.coords, seed=s, live_n_true=truth.n)
    return {k: h[k] for k in ("prefix_mode", "dice_proper", "dice_mirror", "dice_gap")}


def one_seed(t, s: int, data: Path, drop: tuple[str, ...] = ()) -> dict:
    """Both builds of one task-seed: run_task's results and E49's per-seed facts."""
    out: dict = {"builds": {}, "facts": {}}
    names = (C0, *run_names(drop))
    files, refused = {}, {}
    for build in BUILDS:
        with captured() as seen:
            res = harness.run_task(t, candidates(seen, drop), seeds=[s], data_dir=data, build=build, log=None)
        board = seen["board"]
        if _is_test_board(board):
            raise RuntimeError("E49 met a test board")
        out["builds"][build] = res
        files[build] = {n: seen["file"][(n, s)] for n in names}
        refused[build] = {n: seen["refused"].get((n, s)) for n in names if n != C0}
        if build != BUILDS[0]:
            continue
        b = harness.BUILD_SEED_OFFSET + s
        recs = {n: seen["records"].get((FORM_OF.get(n), b), {}) for n in names}
        truth = ev.working_copies(seen["target"], seen["ref"], s, ev.SUBSAMPLE_FRAC, harness.PROTOCOL)[0]
        out["facts"] = {"board_split": str(board.split),
                        "x": {n: seen["x"][(n, s)] for n in names},
                        "n": {n: int(len(seen["stage"][(n, s)].X)) for n in names},
                        "refused": refused[build], "records": recs,
                        "hands": {n: _hand(seen["stage"][(n, s)].coords, truth, s) for n in names}}
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


def gate_check(name: str, key: str, run_scores: dict) -> dict:
    old, why = e47_scores(key)
    worst = {"mean": None, "seed": None}
    if old is not None:
        worst = {"mean": 0.0, "seed": 0.0}
        for arm in GATE_ARMS:
            o, n = old.get(arm), run_scores.get(arm)
            if o is None or n is None or o.shape != n.shape:
                why.append(f"{key}: {arm} missing or of another length")
                continue
            dm, ds = abs(float(n.mean() - o.mean())), float(np.abs(n - o).max())
            worst = {"mean": max(worst["mean"], dm), "seed": max(worst["seed"], ds)}
            if not dm <= GATE_TOL:
                why.append(f"{key}: {arm}'s ten-seed mean is {dm:.4f} off E47's (tolerance {GATE_TOL})")
            if not ds <= GATE_TOL:
                why.append(f"{key}: {arm} is {ds:.4f} off E47's in a seed (tolerance {GATE_TOL})")
    return {"task": name, "key": key, "e47": {"run_commit": E47_RUN_COMMIT, "files_commit": E47_FILES_COMMIT,
                                               "dir": str(E47_DIR.relative_to(REPO))},
            "tolerance": GATE_TOL, "worst": worst, "why": why, "pass": not why}


def cmd_gate(a) -> None:
    _refuse_unless_ready(a.commit)
    t = _task(a.task)
    path = _gate_path(Path(a.runs), a.task)
    if path.exists():
        sys.exit(f"{path} exists")
    res = harness.run_task(t, [harness.CANDIDATES[harness.FLOOR], harness.CANDIDATES[C0]], seeds=SEEDS,
                           data_dir=Path(a.data), build="main", log=None)
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
    describe, why = load_describe(runs, a.commit)
    why = [*why, *gate_problems(runs, a.commit)]
    if why:
        sys.exit("refused (step 0 and the pre-run gate): " + "; ".join(why))
    drop = dropped(describe)
    for c in (C0, *run_names(drop)):
        if not harness.CANDIDATES[c].applies(t):
            sys.exit(f"{a.task}: {c} does not apply")
    folder = runs / "run" / a.task
    folder.mkdir(parents=True, exist_ok=True)
    for s in SEEDS:
        path = folder / f"seed{s}.json"
        if path.exists():
            continue
        res = one_seed(t, s, Path(a.data), drop)
        rec = {"task": a.task, "key": t.key, "seed": s, "run_commit": a.commit, "script": SCRIPT, "dropped": list(drop),
               **res}
        tmp = path.with_suffix(".tmp")
        tmp.write_text(json.dumps(rec, indent=1, default=float))
        os.replace(tmp, path)
        print(f"{a.task} seed {s} done", flush=True)
    print(f"wrote {folder}", flush=True)


# ------------------------------------------------------------------------------------------------ read
def _rows(recs: list[dict], build: str = "main") -> list[dict]:
    return [r["builds"][build]["rows"][0] for r in recs]


def series(recs: list[dict], name: str, what: str, build: str = "main") -> np.ndarray:
    """Per seed, one quantity of one file: 'board', 'group:<g>', 'skill:<m>', 'S' (the ungated shape mean) or
    'mult' (the board score under the multiplicative gate)."""
    rows = _rows(recs, build)
    if what == "board":
        return np.array([x["scores"][name] for x in rows], float)
    if what.startswith("group:"):
        return np.array([x["groups"][name][what[6:]] for x in rows], float)
    if what.startswith("skill:"):
        return np.array([x["skills"][name][what[6:]] for x in rows], float)
    if what == "S":
        return np.array([np.mean([x["skills"][name][m] for m in SHAPE_TERMS]) for x in rows], float)
    if what == "mult":
        return np.array([ev.aggregate("T2", x["skills"][name], gate="multiplicative")[0] for x in rows], float)
    raise ValueError(what)


def merged_D(recs: dict[str, list[dict]], groups: dict[str, list[str]], a: str, b: str, what: str,
             build: str = "main") -> np.ndarray:
    """D_s: the mean over truths of each truth's task-mean difference (a minus b) at each seed."""
    return e47.merged_seeds({n: series(recs[n], a, what, build) - series(recs[n], b, what, build)
                             for names in groups.values() for n in names}, groups)


def _stat(D: np.ndarray) -> dict:
    return {"mean": float(D.mean()), "se": float(D.std(ddof=1) / np.sqrt(len(D))), "ahead": int((D > 0).sum()),
            "behind": int((D < 0).sum())}


def check_run(recs: dict[str, list[dict]], why_in: list[str], commit: str, drop: tuple[str, ...]) -> list[str]:
    """The void conditions on the per-seed records (task -> its records in seed order)."""
    why = list(why_in)
    names = (C0, *run_names(drop))
    for task in nt.all_tasks():
        rs = recs.get(task, [])
        if [r["seed"] for r in rs] != list(SEEDS):
            why.append(f"{task}: seeds are not {list(SEEDS)}")
            continue
        for r in rs:
            tag = f"{task} seed {r['seed']}"
            if r.get("run_commit") != commit:
                why.append(f"{tag}: another run commit")
            if sorted(r.get("dropped", [])) != sorted(drop):
                why.append(f"{tag}: run with another set of dropped arms")
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
                if any(n not in row["scores"] or n not in row["skills"] for n in (harness.FLOOR, *names)):
                    why.append(f"{tag}: build {b} lacks an arm")
                    continue
                for n in names[1:]:
                    d = max(abs(row["skills"][n][m] - row["skills"][C0][m]) for m in EXPR_METRICS)
                    if not d <= EXPR_TOL:
                        why.append(f"{tag}: {n}'s expression-only skills differ from C0's ({b})")
            f = r["facts"]
            if str(f.get("board_split", "")).startswith("test"):
                why.append(f"{tag}: a test board")
            if not f.get("files_identical_across_builds") or not f.get("refusals_identical_across_builds"):
                why.append(f"{tag}: a file or a refusal differs between the builds")
            jit = {n: f["records"].get(n, {}).get("n_jittered", 0) > 0 for n in names}
            for n in names[1:]:
                if f["n"][n] != f["n"][C0] or f["x"][n] != f["x"][C0]:
                    why.append(f"{tag}: {n}'s count or expression matrix differs from C0's")
                if not (jit[n] or jit[C0]) and "main" in r["builds"]:
                    row = r["builds"]["main"]["rows"][0]
                    if n in row["skills"] and C0 in row["skills"]:
                        d = abs(row["skills"][n]["scale_log_ratio"] - row["skills"][C0]["scale_log_ratio"])
                        if not d <= SCALE_TOL:
                            why.append(f"{tag}: {n}'s scale_log_ratio skill is off C0's by {d:.2e}")
    return why


def _per_truth(recs, groups, a: str, b: str) -> dict:
    out = {}
    for truth, names in groups.items():
        out[truth] = float(np.mean([np.mean(series(recs[n], a, "board") - series(recs[n], b, "board"))
                                    for n in names]))
    gains = [v for v in out.values() if v > 0]
    losses = [v for v in out.values() if v < 0]
    return {"per_truth": out, "gain_mean": float(np.mean(gains)) if gains else 0.0, "n_gain": len(gains),
            "loss_mean": float(np.mean(losses)) if losses else 0.0, "n_loss": len(losses)}


def arm_reading(recs, groups, tl_groups, arm: str, refused: dict[str, int]) -> dict:
    ctl = ARMS[arm]
    d = {"board": merged_D(recs, groups, arm, C0, "board"),
         "shape_growth": merged_D(recs, groups, arm, C0, f"group:{GROUP_SHAPE}"),
         "S": merged_D(recs, groups, arm, C0, "S"),
         **{t: merged_D(recs, groups, arm, C0, f"skill:{t}") for t in (*SHAPE_TERMS, "scale_log_ratio",
                                                                         "neighborhood_mmd")},
         **{g: merged_D(recs, groups, arm, C0, f"group:{g}") for g in OTHER_GROUPS}}
    tl = {"board": merged_D(recs, tl_groups, arm, C0, "board"),
          "shape_growth": merged_D(recs, tl_groups, arm, C0, f"group:{GROUP_SHAPE}")}
    vs = {c: merged_D(recs, groups, arm, c, "S") for c in ctl}
    unreadable = [n for n in (arm, *ctl) if refused.get(n, 0)]
    p = {"P1": bool(d["board"].mean() > BAR and (d["board"] > 0).sum() >= AHEAD),
         "P2": bool(d["shape_growth"].mean() > 0),
         "P3": bool(tl["board"].mean() > 0 and tl["shape_growth"].mean() > 0),
         "P4": bool(not any(d[g].mean() < 0 and (d[g] < 0).sum() >= AHEAD for g in OTHER_GROUPS)),
         "P5": bool(all(v.mean() > 0 and (v > 0).sum() >= AHEAD for v in vs.values())),
         "P6": not unreadable}
    return {"try": TRY_IDS[arm], "form": FORM_OF[arm], "controls": list(ctl), "unreadable_because_refused": unreadable,
            "clauses": p, "passes": all(p.values()),
            "merged": {k: _stat(v) for k, v in d.items()}, "test_like": {k: _stat(v) for k, v in tl.items()},
            "against_controls_S": {c: _stat(v) for c, v in vs.items()},
            "prefix_board": _stat(merged_D(recs, groups, arm, C0, "board", build="prefix")),
            "multiplicative_board": _stat(merged_D(recs, groups, arm, C0, "mult")),
            "price": _per_truth(recs, groups, arm, C0)}


def mapping(arms: dict[str, dict], describe: dict, run_blocks: dict[str, set]) -> dict:
    """The outcome-to-file rule, fixed in the row."""
    passing = [a for a, r in arms.items() if r.get("passes")]
    notes = []
    if F1 in passing:
        tl = {b for n in TEST_LIKE for b in run_blocks.get(n, set())}
        tp = describe["drop"]["test_pair_blocks"]
        if tp not in tl or any(len(run_blocks.get(n, set())) != 1 for n in TEST_LIKE):
            passing.remove(F1)
            notes.append("F1 passed but its test-pair blocks equal no test-like task's blocks in this run: unpriced at "
                         "the test, not built")
    ranked = sorted(passing, key=lambda a: (-arms[a]["merged"]["board"]["mean"], -arms[a]["test_like"]["board"]["mean"],
                                            TIE_ORDER.index(a)))[:MAX_CANDIDATES]
    cands = [{"arm": a, "form": FORM_OF[a],
              "command": ("experiments/builds/build_test.py embryo --method fgw_shape0 --positions floor_exact "
                          f"--expression mixture --shape-interp {FORM_OF[a]} --shape-interp-test-ok "
                          f"--out runs/builds/T2embryo_{FORM_OF[a]}.h5ad"),
              "price": arms[a]["price"], "test_like_board": arms[a]["test_like"]["board"]["mean"],
              "merged_board": arms[a]["merged"]["board"]["mean"]} for a in ranked]
    if not cands:
        notes.append("no candidate file: the exact-position file and the series twin stand as handed over")
    return {"candidates": cands, "notes": notes,
            "row_order": "copy_parity is offered beside a shape-only candidate only by E48's own reading and the "
                         "captain's ruling of 2026-11-20 (m), as a separate file labelled unpriced under prefix (J1)"}


def reading(recs: dict[str, list[dict]], why_in: list[str], commit: str, describe: dict | None) -> dict:
    if describe is None:
        return {"void": True, "why": [*why_in, "no describe.json (step 0)"], "commit": commit}
    drop = dropped(describe)
    why = check_run(recs, why_in, commit, drop)
    if why:
        return {"void": True, "why": why, "commit": commit, "dropped": list(drop)}
    tasks = nt.all_tasks()
    groups = e47.truths(tasks)
    tl_groups = e47.truths({n: tasks[n] for n in TEST_LIKE})
    names = run_names(drop)
    refused = {n: sum(bool(r["facts"]["refused"].get(n)) for rs in recs.values() for r in rs) for n in names}
    run_blocks = {t: {r["facts"]["records"].get(F1, {}).get("shape_interp_blocks") for r in rs} - {None}
                  for t, rs in recs.items()}
    arms = {a: arm_reading(recs, groups, tl_groups, a, refused) if a not in drop else {"dropped": True}
            for a in ARMS}
    per_task = {}
    for t, rs in recs.items():
        per_task[t] = {}
        for n in names:
            d = series(rs, n, "board") - series(rs, C0, "board")
            per_task[t][n] = {"board": {"mean": float(d.mean()), "ahead": int((d > 0).sum()),
                                        "sd": float(d.std(ddof=1))},
                              **{m: float(np.mean(series(rs, n, f"skill:{m}") - series(rs, C0, f"skill:{m}")))
                                 for m in (*SHAPE_TERMS, "neighborhood_mmd")},
                              "dice_proper": float(np.mean([r["facts"]["hands"][n]["dice_proper"] for r in rs])),
                              "dice_mirror": float(np.mean([r["facts"]["hands"][n]["dice_mirror"] for r in rs])),
                              "thin_sd": [r["facts"]["records"].get(n, {}).get("shape_interp_thin_sd") for r in rs]}
    split = {a: {k: float(np.mean([per_task[t][a]["board"]["mean"] for t in ks]))
                 for k, ks in (("DEEP", nt.DEEP), ("OTHER", nt.OTHER))} for a in ARMS if a not in drop}
    return {"void": False, "why": [], "commit": commit, "dropped": list(drop), "describe_drop": describe["drop"],
            "truths": groups, "test_like_truths": tl_groups, "refused_task_seeds": refused,
            "f1_blocks_per_task": {t: sorted(v) for t, v in run_blocks.items()}, "arms": arms, "tasks": per_task,
            "deep_other": split, "mapping": mapping(arms, describe, run_blocks),
            "gate": {"e47_run_commit": E47_RUN_COMMIT, "e47_files_commit": E47_FILES_COMMIT, "tolerance": GATE_TOL}}


def text(res: dict) -> str:
    if res["void"]:
        return f"Row E49 (run commit {res['commit']}): VOID: " + "; ".join(res["why"]) + "\n"
    L = [f"Row E49 (run commit {res['commit']}); the pre-run gate passed on all nine tasks (the floor and C0 under "
         f"main within {res['gate']['tolerance']} of E47, run at {res['gate']['e47_run_commit']}, committed in "
         f"{res['gate']['e47_files_commit']}).",
         f"Step 0 (target-free, test pair): F1 blocks {res['describe_drop']['test_pair_blocks']}; dropped before the "
         f"read: {', '.join(res['dropped']) or 'none'}.",
         "Readings fixed in row E49, build main, merged over five truths (mean, SE, seeds ahead of 10):"]
    for a, r in res["arms"].items():
        if r.get("dropped"):
            L.append(f"  {a}: dropped at step 0, not run, not counted")
            continue
        m = r["merged"]
        L.append(f"  {r['try']} {a} ({r['form']}): {'PASSES' if r['passes'] else 'does not pass'}; clauses "
                 + ", ".join(f"{k} {'holds' if v else 'fails'}" for k, v in r["clauses"].items()))
        L.append(f"    board {m['board']['mean']:+.3f} (SE {m['board']['se']:.3f}, {m['board']['ahead']} ahead); gated "
                 f"shape_growth {m['shape_growth']['mean']:+.4f}; ungated S {m['S']['mean']:+.4f}; d2_shape "
                 f"{m['d2_shape']['mean']:+.4f}; occupancy_dice {m['occupancy_dice']['mean']:+.4f}; "
                 f"neighborhood_mmd {m['neighborhood_mmd']['mean']:+.4f}")
        L.append(f"    test-like: board {r['test_like']['board']['mean']:+.3f}, gated shape_growth "
                 f"{r['test_like']['shape_growth']['mean']:+.4f}")
        for c, v in r["against_controls_S"].items():
            L.append(f"    against {TRY_IDS[c]} {c}: ungated S {v['mean']:+.4f} (SE {v['se']:.4f}, {v['ahead']} ahead)")
        if r["unreadable_because_refused"]:
            L.append(f"    unreadable (refused on some task-seed): {', '.join(r['unreadable_because_refused'])}")
        L.append(f"    reported, deciding nothing: prefix board {r['prefix_board']['mean']:+.3f}; multiplicative-gate "
                 f"board {r['multiplicative_board']['mean']:+.3f}; per truth "
                 + ", ".join(f"{k} {v:+.3f}" for k, v in r["price"]["per_truth"].items()))
    L.append("Candidate files (the outcome-to-file rule fixed in the row):")
    for c in res["mapping"]["candidates"]:
        L.append(f"  {c['arm']}: {c['command']}")
    for n in res["mapping"]["notes"]:
        L.append(f"  {n}")
    L.append(f"  {res['mapping']['row_order']}")
    L.append("Reported, deciding nothing: refused task-seeds " + ", ".join(f"{k} {v}" for k, v in
                                                                         res["refused_task_seeds"].items())
             + "; F1 blocks per task " + ", ".join(f"{k} {v}" for k, v in res["f1_blocks_per_task"].items())
             + "; DEEP/OTHER board split " + "; ".join(f"{a} {v['DEEP']:+.3f}/{v['OTHER']:+.3f}"
                                                       for a, v in res["deep_other"].items()))
    L.append("Each candidate goes to the captain as a priced question (all tries are diagnostics, J5). Nothing here "
             "shows E7.75's shape, the live gate, build, metric seed or LAPACK signs.")
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
            "shape_interp_sha256": _sha_file(Path(si.__file__))}


def cmd_read(a) -> None:
    out = Path(a.out)
    files = [out / "reading.json", out / "reading.txt"]
    existing = [str(f) for f in files if f.exists()]
    if existing:
        sys.exit(f"refusing to overwrite: {', '.join(existing)}")
    _own_code()
    if not has_row():
        sys.exit("E49 has no row in experiments/harness/TRIES.md")
    runs = Path(a.runs)
    describe, why = load_describe(runs, a.commit)
    why = [*why, *gate_problems(runs, a.commit), *commit_problems(a.commit)]
    res = {**reading(load_runs(runs), why, a.commit, describe), "reading_code": provenance()}
    out.mkdir(parents=True, exist_ok=True)
    files[0].write_text(json.dumps(res, indent=1, default=float))
    files[1].write_text(text(res))
    print(f"wrote {files[1]}", flush=True)


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    d = sub.add_parser("describe")
    d.add_argument("--commit", required=True)
    d.add_argument("--runs", default=str(RUNS))
    d.add_argument("--data", default=str(harness.DATA_DIR))
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
    {"describe": cmd_describe, "gate": cmd_gate, "run": cmd_run, "read": cmd_read}[a.cmd](a)


if __name__ == "__main__":
    main()
