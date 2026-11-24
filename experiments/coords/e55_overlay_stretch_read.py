#!/usr/bin/env python
"""Row E55 of experiments/harness/TRIES.md (one row, placed by the main session before any run): the practice read of
OSR, the overlay's union stretched along its own principal axes to the brackets' interpolated profile, on the nine
practice tasks, as the cold challenger's specification (round_1124_t2, embryo_osr_attack.json) fixes it. Cloned from
experiments/coords/e49_shape_interp_read.py and changed only where that specification says. No number here comes from
a test score. Diagnostic (J5); the composition was chosen after reading E49 on the same nine tasks (selected on).

Arm and placebos (harness candidates, vec.harness.OSR_CANDIDATES; each is fgw_shape0x_ov with one composed form of
vec.shape_interp applied before the size rule), all three counted when the run step starts:
  E55a fgw_shape0x_ov_sr       overlay_stretch       the arm
  E55b fgw_shape0x_ov_sr_anti  overlay_stretch_anti  its placebo (the target from log_interp_sds(..., anti=True))
  E55c fgw_shape0x_ov_sr_axes  overlay_stretch_axes  its placebo (the arm's own f along a random frame, AXES_STREAM)
Control, rescored and not counted: O = fgw_shape0x_ov (E49c). The floor is run beside it, deciding nothing in the run.

Path. harness.run_task on each task of nine_tasks.all_tasks() (half, late, new1, new2, e75, e725, x1, x2, x3), seeds
0-9, bake protocol, build seed BUILD_SEED_OFFSET + seed, once with build "main" (primary) and once with "prefix"
(reported only). A refusal of a form (vec.shape_interp.ShapeRefused, or the harness's cloud guard) is caught and
recorded, and O's file of that seed stands in its place.

Step 0, the gate record (the specification author's ruling 4). The target-free gate already ran at commit d145d7882e41
(branch campaign/emb-osr-gate, experiments/coords/osr_gate.py, merged here); its folder (the ten pair files and
gate.json) is the read's --gate-record (default experiments/coords/results/osr_gate_g1), and the sha256 of its
gate.json, quoted in the row, its --gate-sha256 (default GATE_SHA256). The launcher reruns osr_gate.py pair for the ten
pairs at the run commit into <runs>/gate0 and osr_gate.py read there. The run and the read refuse or are void unless:
the record's gate.json has that sha256 and says PASS, its pair files are from d145d78; the recompute's pair files are
from the run commit and its gate.json says PASS (no stop clause S1 to S5); every pair's ten per-seed g equal the
record's within 1e-9. P7's task set is fixed from the record's median g.

Steps (each refuses unless TRIES.md holds a row starting "| E55 ", the tree is clean and HEAD is the run commit):
  gate --task NAME --commit SHA  the floor, fgw_shape0x (C0) and O under main, ten seeds: the floor and C0 against E47's
                                 committed scores within 0.01 board points per seed and per ten-seed mean (E49's gate);
                                 O's ten-seed mean minus C0's within 0.01 of E49's committed per-task figure (ruling 2:
                                 E49 the only reference, no per-seed check of O); exit 3 on a failure
  run --task NAME --commit SHA   O, the arm and its placebos, both builds, seed by seed (resumable); refused unless the
  --gate-record DIR              gate record and its recompute hold (step 0) and all nine gate files pass at this commit
  --gate-sha256 HEX
  read --commit SHA              the reading into experiments/coords/results/e55_overlay_stretch/ (reading.json,
  --gate-record DIR              reading.txt), never overwriting
  --gate-sha256 HEX
The script prints only "<task> seed <s> done", "wrote <path>" and "gate <task> pass|fail". It reads no test board,
test file or test anchor, and prints no count.

Refusals (ruling 6): a refusal of the arm or of either placebo on any task-seed makes the arm unreadable; the outcome
is STOP (no file), not VOID; O's file stands in; the three tries stay counted. P7 (ruling 5): the J4 merge over the
truths that keep a four-truth task of recorded |median g| >= 0.05, at equal weight; above 0, no seeds-ahead count.
Second set of rulings: the test's range is bounded by the J4 merge of the small-g test-like tasks (half, new1, new2)
against O, 'no gain at the test' when at or below 0 (small_g_range); VOID wins over a STOP (the void conditions are
checked before any clause); a step-0 failure makes the run step refuse to start (nothing counted), a void found after
it starts leaves the three tries counted; a refusal of O, C0 or the floor is not caught (the task stops, the missing
seeds void the read; nothing stands in for a control); the 1e-9 tolerance of the gate recompute is fixed (G_TOL, no
option changes it); after a VOID nothing follows under this row and a rerun needs a new row (VOID_TEXT).
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
import nine_tasks as nt  # noqa: E402
import osr_gate as og  # noqa: E402  (the gate's own computation, merged from campaign/emb-osr-gate)

from vec import evaluate as ev  # noqa: E402
from vec import harness  # noqa: E402
from vec import shape_interp as si  # noqa: E402
from vec import t2_build as tb  # noqa: E402
from vec.adoption import SEEDS  # noqa: E402
from vec.provenance import repo_state  # noqa: E402

ROW_PREFIX = "| E55 "  # the row's line in TRIES.md starts with this
TRIES = REPO / "experiments" / "harness" / "TRIES.md"
ROW_FIRST_LINE = ("diagnostic (J5, noise floor 6.54); composition chosen after reading E49 on the same nine tasks "
                  "(selected on); no J8 use")
C0 = "fgw_shape0x"
O = "fgw_shape0x_ov"  # noqa: E741  (the specification names the control O)
ARM = "fgw_shape0x_ov_sr"
PLACEBOS = ("fgw_shape0x_ov_sr_anti", "fgw_shape0x_ov_sr_axes")
COUNTED = (ARM, *PLACEBOS)
TRY_IDS = {ARM: "E55a", PLACEBOS[0]: "E55b", PLACEBOS[1]: "E55c"}
FORM_OF = {**harness.OSR_CANDIDATES, O: "overlay"}
RUN_NAMES = (O, *COUNTED)  # with the floor, in each run_task call
BUILDS = ("main", "prefix")
TEST_LIKE = ("half", "late", "new1", "new2")
E725_TASK = "e725"  # the four-truth merge leaves out its truth, E7.25 (rule 2's known feature)
GROUP_SHAPE = "shape_growth"
OTHER_GROUPS = ("expression_change", "cell_state", "local_organisation")
SHAPE_TERMS = ("d2_shape", "occupancy_dice")
EXPR_METRICS = e47.EXPR_METRICS
EXPR_TOL = 1e-9
SCALE_TOL = 1e-6
BAR, AHEAD = 0.3, 8
G_MIN = 0.05  # P7: tasks whose gate |median g| is at least this
G_TOL = 1e-9  # ruling 4: the recompute's per-seed g against the record's
GATE_COMMIT = "d145d7882e4113b450c643a5c23ca225f7c57ff0"  # the target-free gate's run (branch campaign/emb-osr-gate)
GATE_RECORD = REPO / "experiments" / "coords" / "results" / "osr_gate_g1"  # its ten pair files and gate.json
GATE_SHA256 = "7c0fc6b6398092d54f6e4318bd8c79c1e092c4ec7e1edbbb6b86fcbede7b4d67"  # of that gate.json, as the row
E47_DIR = REPO / "experiments" / "harness" / "results" / "e47_floor_exact"
E47_RUN_COMMIT, E47_FILES_COMMIT = "9967cf1", "d5ee9b5"
E49_DIR = REPO / "experiments" / "coords" / "results" / "e49_shape_interp"
E49_RUN_COMMIT, E49_FILES_COMMIT = "d76ccd5", "53ba6cc"
GATE_TOL = 0.01  # board points, per seed and per ten-seed mean (nine_tasks.TOL)
O_TOL = 0.01  # board points, the O guard per task: the ten-seed mean of O minus C0 against E49's (ruling 2)
GATE_RUN = (harness.FLOOR, C0, O)
GATE_ARMS = (harness.FLOOR, C0)
BASES = (E47_FILES_COMMIT, E49_FILES_COMMIT, GATE_COMMIT)
GUARD_FLAG = "--shape-interp-test-ok"
SCRIPT = "experiments/coords/e55_overlay_stretch_read.py"
RUNS = REPO / "runs" / "e55"
RESULTS = REPO / "experiments" / "coords" / "results" / "e55_overlay_stretch"
REFUSAL = "cloud guard flags"
LABEL = "diagnostic, selected on the nine tasks"
CANDIDATE_COMMAND = ("experiments/builds/build_test.py embryo --method fgw_shape0 --positions floor_exact "
                     "--expression mixture --shape-interp overlay_stretch --shape-interp-test-ok "
                     "--out runs/builds/T2embryo_overlay_stretch.h5ad")

# THE PASS RULE, word for word from the specification's practice_read.pass_rule
PASS_RULE = (
    "All under build main, against O. P1 four-truth merged board gain above +0.3, D_s above 0 in at least 8 of 10. "
    "P1b five-truth merged board gain above +0.3, D_s above 0 in at least 8 of 10. P2 four-truth merged gated "
    "shape_growth above 0. P3 test-like board and gated shape_growth above 0. P4 behind on no other group "
    "(expression_change, cell_state, local_organisation) on the five-truth or four-truth merge by E49's definition "
    "(merged mean below 0 with at least 8 of 10 D_s below 0), with no tolerance. P5 ungated S = (d2_shape skill + "
    "occupancy_dice skill)/2 of the arm above each placebo's on the four-truth merge, at least 8 of 10 seeds ahead. P6 "
    "no refusal of the arm. P7 the four-truth merged board gain is above 0 when restricted to tasks whose gate "
    "|median g| is at least 0.05. Void as E49: count or expression matrix off O's; expression-only skills off O's by "
    "more than 1e-9; scale skill off O's by more than 1e-6 where neither file jittered; mixed commits; an unclean "
    "tree; a test board; the gate record missing or from another commit. Pass = one priced naming question to the "
    "captain labelled 'diagnostic, selected on the nine tasks', with gains and losses per truth and the four-truth and "
    "test-like means. Anything else = no file.")
CLAUSES = {
    "P1": "P1 four-truth merged board gain above +0.3, D_s above 0 in at least 8 of 10.",
    "P1b": "P1b five-truth merged board gain above +0.3, D_s above 0 in at least 8 of 10.",
    "P2": "P2 four-truth merged gated shape_growth above 0.",
    "P3": "P3 test-like board and gated shape_growth above 0.",
    "P4": "P4 behind on no other group (expression_change, cell_state, local_organisation) on the five-truth or "
          "four-truth merge by E49's definition (merged mean below 0 with at least 8 of 10 D_s below 0), with no "
          "tolerance.",
    "P5": "P5 ungated S = (d2_shape skill + occupancy_dice skill)/2 of the arm above each placebo's on the four-truth "
          "merge, at least 8 of 10 seeds ahead.",
    "P6": "P6 no refusal of the arm.",
    "P7": "P7 the four-truth merged board gain is above 0 when restricted to tasks whose gate |median g| is at least "
          "0.05.",
}
VOID_CONDITIONS = (
    "count or expression matrix off O's",
    "expression-only skills off O's by more than 1e-9",
    "scale skill off O's by more than 1e-6 where neither file jittered",
    "mixed commits",
    "an unclean tree",
    "a test board",
    # 'the gate record missing or from another commit', replaced by the author's ruling 4:
    "the gate record missing, any stop clause S1-S5 firing there, or any pair's ten per-seed g values differing from "
    "the d145d78 record by more than 1e-9",
)
VOID_AS_E49 = ("a missing task, seed, build or arm", "a gate file missing or failing (the O guard included)",
               "a file not identical across the two builds", "a refusal that differs between the builds")
VOID_TEXT = ("VOID = read as nothing and no file. No task or seed is re-run selectively, and nothing further follows "
             "under this row. A rerun needs a new row placed before it, and its tries count again. Tries counted at "
             "the run step's start stay counted (rule 4).")
RULINGS = {  # the specification author's rulings on its open points (embryo_osr_open_points.json), as coded here
    "O guard": "per task, the same-run ten-seed mean of O minus C0 under main within 0.01 of E49's committed "
               "tasks[name]['fgw_shape0x_ov'].board.mean; no per-seed check of O; C0 and the floor against E47 per "
               "seed and per ten-seed mean within 0.01. E49 is the only reference; E51's pass is the J7 dry run.",
    "gate record": "accepted by its commit d145d7882e41 and the sha256 of its gate.json; recomputed at the run commit "
                   "with osr_gate.py's own computation and count; void if missing, if S1-S5 fire there, or if any g "
                   "differs by more than 1e-9.",
    "P7": "the J4 merge over the truths keeping a four-truth task of recorded gate |median g| >= 0.05, at equal "
          "weight; mean_s D_s > 0, no seeds-ahead count.",
    "refusal": "a refusal of the arm or of either placebo on any task-seed makes the arm unreadable: STOP (no "
               "file), not VOID; O's file stands in; the three tries stay counted.",
    "count": "the gate's count reading on the twin-board tasks changes nothing in the read; P7's g on x1 to x3 was "
             "measured at a count other than the one the read scores.",
}


# ------------------------------------------------------------------------------------------------ checks
def _git(*args: str) -> subprocess.CompletedProcess:
    return subprocess.run(["git", "-C", str(REPO), *args], capture_output=True, text=True,
                          env={**os.environ, "GIT_OPTIONAL_LOCKS": "0"})


def has_row(path: Path | None = None) -> bool:
    """Whether TRIES.md holds a row whose line starts with "| E55 "."""
    path = TRIES if path is None else path
    if not Path(path).exists():
        return False
    return any(ln.startswith(ROW_PREFIX) for ln in Path(path).read_text(encoding="utf-8").splitlines())


def commit_problems(commit: str) -> list[str]:
    """Why this checkout may not run or read row E55 at `commit` (empty: it may)."""
    why = []
    head = _git("rev-parse", "HEAD").stdout.strip()
    if head != commit:
        why.append(f"HEAD is {head[:12]}, not the run commit {commit[:12]}")
    for base in BASES:
        if _git("merge-base", "--is-ancestor", base, "HEAD").returncode != 0:
            why.append(f"HEAD does not contain {base}")
    if GUARD_FLAG not in (REPO / "experiments" / "builds" / "build_test.py").read_text(encoding="utf-8"):
        why.append(f"the builder lacks the test-board guard {GUARD_FLAG}")
    for d, c in ((E47_DIR, E47_FILES_COMMIT), (E49_DIR, E49_FILES_COMMIT)):
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
        sys.exit("E55 has no row in experiments/harness/TRIES.md: a row is committed before its run")
    if repo_state(REPO).unclean != ():
        sys.exit("the tree is not clean: row E55 runs on a committed tree")
    why = commit_problems(commit)
    if why:
        sys.exit("refused: " + "; ".join(why))


def _task(name: str):
    tasks = nt.all_tasks()
    if name not in tasks:
        sys.exit(f"unknown task {name}; known: {sorted(tasks)}")
    t = tasks[name]
    if getattr(t, "target_filter", None) is not None:
        sys.exit(f"{name}: a task with a target filter; E55 does not run on one")
    return t


def _is_test_board(board) -> bool:
    return str(board.split).startswith("test") or "test" in str(board.key).split(":")[-1]


def _sha_file(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


# ------------------------------------------------------------------------------------------------ the gate record
def _gate_ok_sha(path: Path, sha: str) -> bool:
    return Path(path).exists() and _sha_file(Path(path)) == sha


def gate0_problems(record: Path | None, sha: str | None, recompute: Path, commit: str) -> tuple[list[str], dict | None]:
    """Ruling 4: the earlier gate record (folder `record`, run at GATE_COMMIT, its gate.json of sha256 `sha`) and its
    recompute at the run commit (folder `recompute`, written by osr_gate.py pair for the ten pairs at the run commit).
    Void reasons: the record missing, its gate.json not of that sha256, a pair file of the record from another commit
    than GATE_COMMIT; the recompute missing or from another commit than the run commit; a stop clause S1 to S5 firing
    in the record or in the recompute; any pair's ten per-seed g differing from the record by more than G_TOL.
    Returns (reasons, the record's gate.json) where the record's gate.json is read (P7's task set is fixed from it)."""
    why: list[str] = []
    if record is None or sha is None or not (Path(record) / "gate.json").exists():
        return ["the gate record (gate.json of the d145d78 run) is missing"], None
    record = Path(record)
    if not _gate_ok_sha(record / "gate.json", sha):
        return [f"the gate record's gate.json is not the one of sha256 {sha[:12]}"], None
    rec = json.loads((record / "gate.json").read_text())
    if rec.get("verdict") != "PASS" or rec.get("stopped_by"):
        why.append("a stop clause fires in the gate record: " + ", ".join(rec.get("stopped_by") or ["verdict"]))
    if rec.get("commits") != [GATE_COMMIT]:
        why.append(f"the gate record is not from {GATE_COMMIT[:12]} alone")
    for name in og.NAMES:
        old_p, new_p = record / f"{name}.json", Path(recompute) / f"{name}.json"
        if not old_p.exists():
            why.append(f"the gate record lacks {name}")
            continue
        if not new_p.exists():
            why.append(f"the gate recompute at the run commit lacks {name}")
            continue
        old, new = json.loads(old_p.read_text()), json.loads(new_p.read_text())
        if old.get("commit") != GATE_COMMIT:
            why.append(f"the gate record's {name} is from another commit")
        if new.get("commit") != commit:
            why.append(f"the gate recompute's {name} is from another commit than the run commit")
        po, pn = old.get("per_seed") or {}, new.get("per_seed") or {}
        if sorted(po) != sorted(pn) or sorted(po) != sorted(str(x) for x in og.SEEDS):
            why.append(f"{name}: the recompute's seeds are not the record's")
            continue
        for sd in po:
            go, gn = po[sd].get("g"), pn[sd].get("g")
            if go is None or gn is None or not abs(float(go) - float(gn)) <= G_TOL:
                why.append(f"{name} seed {sd}: g differs from the record by more than {G_TOL}")
    tp = record / f"{og.TEST}.json"
    rec["test_pair_g"] = ([v.get("g") for v in (json.loads(tp.read_text()).get("per_seed") or {}).values()
                           if v.get("g") is not None] if tp.exists() else [])
    if (Path(recompute) / "gate.json").exists():
        rc = json.loads((Path(recompute) / "gate.json").read_text())
        if rc.get("verdict") != "PASS" or rc.get("stopped_by"):
            why.append("a stop clause fires in the gate recompute: " + ", ".join(rc.get("stopped_by") or ["verdict"]))
    else:
        why.append("the gate recompute's gate.json is missing")
    return why, rec


def p7_tasks(median_g: dict[str, float], groups4: dict[str, list[str]]) -> tuple[list[str], list[str]]:
    """Ruling 3: the four-truth tasks whose recorded gate |median g| is at least G_MIN, and those below."""
    four = [n for names in groups4.values() for n in names]
    at = sorted(n for n in four if median_g.get(n) is not None and abs(float(median_g[n])) >= G_MIN)
    return at, sorted(n for n in four if n not in at)


def tp_statement(median_g: dict[str, float], per_seed_g: list[float]) -> dict:
    """The recorded gate's statement about the test pair (anything_else_the_row_must_say), deciding nothing."""
    tp = median_g.get(og.TEST)
    if tp is None:
        return {"text": "the gate record holds no test-pair median g"}
    near = sorted(TEST_LIKE, key=lambda n: abs(abs(float(median_g[n])) - abs(float(tp)))
                  if median_g.get(n) is not None else float("inf"))[:2]
    below = sum(abs(float(g)) < G_MIN for g in per_seed_g)
    small = [n for n in TEST_LIKE if median_g.get(n) is not None and abs(float(median_g[n])) < G_MIN]
    text = (f"The test pair's median g is {float(tp):.4f}, and {below} of ten seeds are below {G_MIN}"
            + (f": the test pair falls in the class of tasks that P7 sets aside as below {G_MIN}"
               if abs(float(tp)) < G_MIN else "")
            + f". The test's step is nearest in size to {near[0]} ({float(median_g[near[0]]):.3f}) and {near[1]} "
            f"({float(median_g[near[1]]):.3f}); a practice gain, if any, applies at the test only as a range from zero "
            "up to the gain on the small-g test-like tasks (" + (", ".join(small) or "none") + ").")
    return {"test_pair_median_g": float(tp), "seeds_below": below, "nearest": near, "small_g_test_like": small,
            "below_p7": bool(abs(float(tp)) < G_MIN), "text": text}


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
    """O unchanged; its file is kept for the checks and as the stand-in for a refused form."""
    def build(ctx):
        if _is_test_board(seen["board"]):
            raise RuntimeError("E55 met a test board")
        pred = base.build(ctx)
        _keep(seen, base.name, ctx.seed - harness.BUILD_SEED_OFFSET, pred)
        return pred

    return dataclasses.replace(base, build=build)


def guarded(base, seen):
    """The arm or a placebo: its file, or on a refusal (ShapeRefused, the cloud guard) O's file of the seed."""
    def build(ctx):
        s = ctx.seed - harness.BUILD_SEED_OFFSET
        if _is_test_board(seen["board"]):
            raise RuntimeError("E55 met a test board")
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
        _keep(seen, base.name, s, pred)
        return pred

    return dataclasses.replace(base, build=build)


def candidates(seen) -> list:
    return [harness.CANDIDATES[harness.FLOOR], observed(harness.CANDIDATES[O], seen),
            *(guarded(harness.CANDIDATES[n], seen) for n in COUNTED)]


def _hand(coords, truth, s: int) -> dict:
    h = ev.handedness(coords, truth.coords, seed=s, live_n_true=truth.n)
    return {k: h[k] for k in ("prefix_mode", "dice_proper", "dice_mirror", "dice_gap")}


def one_seed(t, s: int, data: Path) -> dict:
    """Both builds of one task-seed: run_task's results and E55's per-seed facts."""
    out: dict = {"builds": {}, "facts": {}}
    files, refused = {}, {}
    for build in BUILDS:
        with captured() as seen:
            res = harness.run_task(t, candidates(seen), seeds=[s], data_dir=data, build=build, log=None)
        board = seen["board"]
        if _is_test_board(board):
            raise RuntimeError("E55 met a test board")
        out["builds"][build] = res
        files[build] = {n: seen["file"][(n, s)] for n in RUN_NAMES}
        refused[build] = {n: seen["refused"].get((n, s)) for n in COUNTED}
        if build != BUILDS[0]:
            continue
        b = harness.BUILD_SEED_OFFSET + s
        recs = {n: seen["records"].get((FORM_OF[n], b), {}) for n in RUN_NAMES}
        truth = ev.working_copies(seen["target"], seen["ref"], s, ev.SUBSAMPLE_FRAC, harness.PROTOCOL)[0]
        out["facts"] = {"board_split": str(board.split),
                        "x": {n: seen["x"][(n, s)] for n in RUN_NAMES},
                        "n": {n: int(len(seen["stage"][(n, s)].X)) for n in RUN_NAMES},
                        "refused": refused[build], "records": recs,
                        "hands": {n: _hand(seen["stage"][(n, s)].coords, truth, s) for n in RUN_NAMES}}
    out["facts"]["files_identical_across_builds"] = all(files[BUILDS[0]][n] == files[b][n]
                                                        for b in BUILDS[1:] for n in RUN_NAMES)
    out["facts"]["refusals_identical_across_builds"] = all(refused[BUILDS[0]] == refused[b] for b in BUILDS[1:])
    return out


# ------------------------------------------------------------------------------------------------ the O guard
def e49_o_per_task(path: Path | None = None) -> tuple[dict | None, list[str]]:
    """Ruling 2: E49's committed per-task figure of O, its ten-seed mean board minus C0's under main
    (reading.json tasks[name]['fgw_shape0x_ov'].board.mean). E49 is the only reference."""
    p = E49_DIR / "reading.json" if path is None else Path(path)
    if not p.exists():
        return None, ["no committed E49 reading"]
    r = json.loads(p.read_text())
    why = []
    if r.get("void") is not False or not str(r.get("commit", "")).startswith(E49_RUN_COMMIT):
        why.append(f"E49's reading is void or not from its run at {E49_RUN_COMMIT}")
    out = {}
    for name in nt.all_tasks():
        v = (((r.get("tasks") or {}).get(name) or {}).get(O) or {}).get("board", {}).get("mean")
        if v is None:
            why.append(f"{name}: no committed E49 figure for {O}")
        else:
            out[name] = float(v)
    return out, why


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


def gate_check(name: str, key: str, run_scores: dict, o_ref: tuple[dict | None, list[str]] | None = None) -> dict:
    """E49's gate on one task (the floor and C0 against E47 within GATE_TOL per seed and ten-seed mean), and the O guard
    (ruling 2): O's ten-seed mean minus C0's within O_TOL of E49's committed per-task figure. `o_ref` defaults to
    e49_o_per_task()."""
    old, why = e47_scores(key)
    ref, why_ref = e49_o_per_task() if o_ref is None else o_ref
    why = [*why, *why_ref]
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
            if not dm <= GATE_TOL:
                why.append(f"{key}: {arm}'s ten-seed mean is {dm:.4f} off E47's (tolerance {GATE_TOL})")
            if not ds <= GATE_TOL:
                why.append(f"{key}: {arm} is {ds:.4f} off E47's in a seed (tolerance {GATE_TOL})")
    o_minus_c0 = None
    if run_scores.get(O) is not None and run_scores.get(C0) is not None \
            and np.shape(run_scores[O]) == np.shape(run_scores[C0]) == (len(SEEDS),):
        o_minus_c0 = float(np.mean(np.asarray(run_scores[O], float) - np.asarray(run_scores[C0], float)))
    else:
        why.append(f"{key}: {O} or {C0} missing or of another length")
    o_off = None
    if o_minus_c0 is not None:
        if ref is None or name not in ref:
            why.append(f"{name}: no committed E49 figure for {O}")
        else:
            o_off = abs(o_minus_c0 - ref[name])
            if not o_off <= O_TOL:
                why.append(f"{name}: {O} minus {C0} is {o_off:.4f} off E49's ten-seed mean (tolerance {O_TOL})")
    return {"task": name, "key": key, "e47": {"run_commit": E47_RUN_COMMIT, "files_commit": E47_FILES_COMMIT,
                                               "dir": str(E47_DIR.relative_to(REPO))},
            "e49": {"run_commit": E49_RUN_COMMIT, "files_commit": E49_FILES_COMMIT}, "tolerance": GATE_TOL,
            "o_tolerance": O_TOL, "worst": worst, "o_off_e49": o_off, "why": why, "pass": not why}


def cmd_gate(a) -> None:
    _refuse_unless_ready(a.commit)
    t = _task(a.task)
    path = _gate_path(Path(a.runs), a.task)
    if path.exists():
        sys.exit(f"{path} exists")
    res = harness.run_task(t, [harness.CANDIDATES[n] for n in GATE_RUN], seeds=SEEDS, data_dir=Path(a.data),
                           build="main", log=None)
    if str(res["board"].get("split", "")).startswith("test"):
        sys.exit("E55 met a test board")
    rec = {**gate_check(a.task, t.key, e47.per_seed(res)[0]), "commit": res["provenance"]["commit"],
           "tree_clean": res["provenance"]["tree_clean"]}
    rec["pass"] = bool(rec["pass"] and rec["commit"] == a.commit and rec["tree_clean"])
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(rec, indent=1, default=float))
    print(f"gate {a.task} {'pass' if rec['pass'] else 'fail'}", flush=True)
    sys.exit(0 if rec["pass"] else 3)


def gate_problems(runs: Path, commit: str) -> list[str]:
    """Every gate file (the floor and C0 against E47, O against E49 per task) present and passing at the run commit."""
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


def _recompute_dir(runs: Path) -> Path:
    return Path(runs) / "gate0"


def cmd_run(a) -> None:
    _refuse_unless_ready(a.commit)
    t = _task(a.task)
    runs = Path(a.runs)
    why, _ = gate0_problems(a.gate_record, a.gate_sha256, _recompute_dir(runs), a.commit)
    why = [*why, *gate_problems(runs, a.commit)]
    if why:
        sys.exit("refused (the gate record and the pre-run gate): " + "; ".join(why))
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
        rec = {"task": a.task, "key": t.key, "seed": s, "run_commit": a.commit, "script": SCRIPT, **res}
        tmp = path.with_suffix(".tmp")
        tmp.write_text(json.dumps(rec, indent=1, default=float))
        os.replace(tmp, path)
        print(f"{a.task} seed {s} done", flush=True)
    print(f"wrote {folder}", flush=True)


# ------------------------------------------------------------------------------------------------ read
def _rows(recs: list[dict], build: str = "main") -> list[dict]:
    return [r["builds"][build]["rows"][0] for r in recs]


def series(recs: list[dict], name: str, what: str, build: str = "main") -> np.ndarray:
    """Per seed, one quantity of one file: 'board', 'group:<g>' (gated), 'skill:<m>' or 'S' (the ungated shape mean)."""
    rows = _rows(recs, build)
    if what == "board":
        return np.array([x["scores"][name] for x in rows], float)
    if what.startswith("group:"):
        return np.array([x["groups"][name][what[6:]] for x in rows], float)
    if what.startswith("skill:"):
        return np.array([x["skills"][name][what[6:]] for x in rows], float)
    if what == "S":
        return np.array([np.mean([x["skills"][name][m] for m in SHAPE_TERMS]) for x in rows], float)
    raise ValueError(what)


def merged_D(recs: dict[str, list[dict]], groups: dict[str, list[str]], a: str, b: str, what: str,
             build: str = "main") -> np.ndarray:
    """D_s: the mean over truths of each truth's task-mean difference (a minus b) at each seed."""
    return e47.merged_seeds({n: series(recs[n], a, what, build) - series(recs[n], b, what, build)
                             for names in groups.values() for n in names}, groups)


def _stat(D: np.ndarray) -> dict:
    return {"mean": float(D.mean()), "se": float(D.std(ddof=1) / np.sqrt(len(D))), "ahead": int((D > 0).sum()),
            "behind": int((D < 0).sum())}


def behind(D: np.ndarray) -> bool:
    """E49's P4 definition: a merged mean below 0 with at least AHEAD of the ten D_s below 0."""
    D = np.asarray(D, float)
    return bool(D.mean() < 0 and (D < 0).sum() >= AHEAD)


def _gain(D) -> bool:
    """Above +BAR with D_s above 0 in at least AHEAD of 10."""
    D = np.asarray(D, float)
    return bool(D.mean() > BAR and (D > 0).sum() >= AHEAD)


def decide(d: dict) -> dict:
    """THE PASS RULE (CLAUSES), all under build main, against O. `d` holds per-seed D_s arrays: board4, board5, sg4
    (four-truth gated shape_growth), board_tl and sg_tl (test-like), groups5 and groups4 ({group: D}), S_vs ({placebo:
    the arm's S minus the placebo's on the four-truth merge}), board4_g (the four-truth board restricted to tasks whose
    gate |median g| >= G_MIN, or None if no task qualifies); and counts: arm_refused (task-seeds the arm refused) and
    placebo_refused ({placebo: task-seeds})."""
    p = {"P1": _gain(d["board4"]),
         "P1b": _gain(d["board5"]),
         "P2": bool(np.mean(d["sg4"]) > 0),
         "P3": bool(np.mean(d["board_tl"]) > 0 and np.mean(d["sg_tl"]) > 0),
         "P4": not any(behind(D) for grp in (d["groups5"], d["groups4"]) for D in grp.values()),
         "P5": bool(all(np.mean(D) > 0 and (np.asarray(D) > 0).sum() >= AHEAD and not d["placebo_refused"].get(c)
                        for c, D in d["S_vs"].items()) and set(d["S_vs"]) == set(PLACEBOS)),
         "P6": not d["arm_refused"],
         "P7": bool(d["board4_g"] is not None and np.mean(d["board4_g"]) > 0)}
    return {"clauses": p, "passes": all(p.values())}


def restricted_groups(groups: dict[str, list[str]], keep) -> dict[str, list[str]]:
    """The merge's truths restricted to the tasks in `keep` (a truth left without a task drops out)."""
    out = {t: [n for n in names if n in keep] for t, names in groups.items()}
    return {t: v for t, v in out.items() if v}


def check_run(recs: dict[str, list[dict]], why_in: list[str], commit: str) -> list[str]:
    """The void conditions on the per-seed records (task -> its records in seed order)."""
    why = list(why_in)
    for task in nt.all_tasks():
        rs = recs.get(task, [])
        if [r["seed"] for r in rs] != list(SEEDS):
            why.append(f"{task}: seeds are not {list(SEEDS)}")
            continue
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
                if str(res["board"].get("split", "")).startswith("test"):
                    why.append(f"{tag}: a test board")
                row = res["rows"][0]
                if any(n not in row["scores"] or n not in row["skills"] or n not in row.get("groups", {})
                       for n in (harness.FLOOR, *RUN_NAMES)):
                    why.append(f"{tag}: build {b} lacks an arm")
                    continue
                for n in COUNTED:
                    d = max(abs(row["skills"][n][m] - row["skills"][O][m]) for m in EXPR_METRICS)
                    if not d <= EXPR_TOL:
                        why.append(f"{tag}: {n}'s expression-only skills differ from O's ({b})")
            f = r["facts"]
            if str(f.get("board_split", "")).startswith("test"):
                why.append(f"{tag}: a test board")
            if not f.get("files_identical_across_builds") or not f.get("refusals_identical_across_builds"):
                why.append(f"{tag}: a file or a refusal differs between the builds")
            jit = {n: f["records"].get(n, {}).get("n_jittered", 0) > 0 for n in RUN_NAMES}
            for n in COUNTED:
                if f["n"][n] != f["n"][O] or f["x"][n] != f["x"][O]:
                    why.append(f"{tag}: {n}'s count or expression matrix differs from O's")
                if not (jit[n] or jit[O]) and "main" in r["builds"]:
                    row = r["builds"]["main"]["rows"][0]
                    if n in row["skills"] and O in row["skills"]:
                        d = abs(row["skills"][n]["scale_log_ratio"] - row["skills"][O]["scale_log_ratio"])
                        if not d <= SCALE_TOL:
                            why.append(f"{tag}: {n}'s scale_log_ratio skill is off O's by {d:.2e}")
    return why


def _per_truth(recs, groups, a: str, b: str) -> dict:
    out = {}
    for truth, names in groups.items():
        out[truth] = float(np.mean([np.mean(series(recs[n], a, "board") - series(recs[n], b, "board"))
                                    for n in names]))
    gains = {k: v for k, v in out.items() if v > 0}
    losses = {k: v for k, v in out.items() if v < 0}
    return {"per_truth": out, "gains": gains, "losses": losses,
            "gain_mean": float(np.mean(list(gains.values()))) if gains else 0.0,
            "loss_mean": float(np.mean(list(losses.values()))) if losses else 0.0}


def merges(tasks: dict | None = None) -> dict[str, dict[str, list[str]]]:
    """The merges fixed in the specification: five-truth, four-truth (without E7.25) and test-like."""
    tasks = nt.all_tasks() if tasks is None else tasks
    five = e47.truths(tasks)
    e725 = tasks[E725_TASK].target
    return {"five": five, "four": {k: v for k, v in five.items() if k != e725},
            "test_like": e47.truths({n: tasks[n] for n in TEST_LIKE})}


def reading(recs: dict[str, list[dict]], why_in: list[str], commit: str, gate: dict | None) -> dict:
    """`gate`: the record's gate.json (P7's task set and the test-pair statement are read from it), or None."""
    if gate is None:
        return {"void": True, "why": [*why_in, "the gate record is missing"], "commit": commit}
    why = check_run(recs, why_in, commit)
    if why:
        return {"void": True, "why": why, "commit": commit}
    m = merges()
    refused = {n: sum(bool(r["facts"]["refused"].get(n)) for rs in recs.values() for r in rs) for n in COUNTED}
    median_g = {k: v for k, v in gate["median_g"].items() if k != og.TEST}
    g_tasks, below = p7_tasks(median_g, m["four"])
    g_groups = restricted_groups(m["four"], g_tasks)
    d = {"board4": merged_D(recs, m["four"], ARM, O, "board"),
         "board5": merged_D(recs, m["five"], ARM, O, "board"),
         "sg4": merged_D(recs, m["four"], ARM, O, f"group:{GROUP_SHAPE}"),
         "board_tl": merged_D(recs, m["test_like"], ARM, O, "board"),
         "sg_tl": merged_D(recs, m["test_like"], ARM, O, f"group:{GROUP_SHAPE}"),
         "groups5": {g: merged_D(recs, m["five"], ARM, O, f"group:{g}") for g in OTHER_GROUPS},
         "groups4": {g: merged_D(recs, m["four"], ARM, O, f"group:{g}") for g in OTHER_GROUPS},
         "S_vs": {c: merged_D(recs, m["four"], ARM, c, "S") for c in PLACEBOS},
         "board4_g": merged_D(recs, g_groups, ARM, O, "board") if g_groups else None,
         "arm_refused": refused[ARM], "placebo_refused": {c: refused[c] for c in PLACEBOS}}
    dec = decide(d)
    below_groups = restricted_groups(m["four"], below)
    per_task = {}
    for t, rs in recs.items():
        per_task[t] = {}
        for n in COUNTED:
            x = series(rs, n, "board") - series(rs, O, "board")
            per_task[t][n] = {"board": {"mean": float(x.mean()), "ahead": int((x > 0).sum()),
                                        "sd": float(x.std(ddof=1))},
                              **{k: float(np.mean(series(rs, n, f"skill:{k}") - series(rs, O, f"skill:{k}")))
                                 for k in (*SHAPE_TERMS, "neighborhood_mmd")},
                              "dice_proper": float(np.mean([r["facts"]["hands"][n]["dice_proper"] for r in rs])),
                              "dice_mirror": float(np.mean([r["facts"]["hands"][n]["dice_mirror"] for r in rs]))}
        per_task[t]["run_median_g"] = _median([r["facts"]["records"].get(ARM, {}).get("shape_interp_g") for r in rs])
        per_task[t]["gate_median_g"] = median_g.get(t)
    price = _per_truth(recs, m["five"], ARM, O)
    small = [n for n in TEST_LIKE if median_g.get(n) is not None and abs(float(median_g[n])) < G_MIN]
    statement = tp_statement(gate["median_g"], gate.get("test_pair_g", []))
    small_groups = e47.truths({n: nt.all_tasks()[n] for n in small}) if small else {}
    statement.update(small_g_range(recs, small_groups))
    unreadable = [n for n in COUNTED if refused[n]]
    res = {"void": False, "why": [], "commit": commit, "truths": m, "refused_task_seeds": refused,
           "unreadable_because_refused": unreadable, "outcome": ("PASS" if dec["passes"] else "STOP"),
           "test_pair": statement, "rulings": RULINGS,
           "clauses": dec["clauses"], "passes": dec["passes"], "clause_text": CLAUSES,
           "merged": {"board4": _stat(d["board4"]), "board5": _stat(d["board5"]), "sg4": _stat(d["sg4"]),
                      "board_tl": _stat(d["board_tl"]), "sg_tl": _stat(d["sg_tl"]),
                      "groups5": {g: {**_stat(v), "is_behind": behind(v)} for g, v in d["groups5"].items()},
                      "groups4": {g: {**_stat(v), "is_behind": behind(v)} for g, v in d["groups4"].items()},
                      "S_vs": {c: _stat(v) for c, v in d["S_vs"].items()},
                      "board4_g": None if d["board4_g"] is None else _stat(d["board4_g"])},
           "p7": {"tasks_at_or_above": g_tasks, "tasks_below": below,
                  "board4_below": (_stat(merged_D(recs, below_groups, ARM, O, "board")) if below_groups else None)},
           "reported_deciding_nothing": {
               "prefix_board5": _stat(merged_D(recs, m["five"], ARM, O, "board", build="prefix")),
               "prefix_board4": _stat(merged_D(recs, m["four"], ARM, O, "board", build="prefix")),
               "placebos_vs_O_board4": {c: _stat(merged_D(recs, m["four"], c, O, "board")) for c in PLACEBOS},
               "test_pair_median_g": gate["median_g"].get(og.TEST),
               "deep_other": {k: float(np.mean([per_task[t][ARM]["board"]["mean"] for t in ks]))
                              for k, ks in (("DEEP", nt.DEEP), ("OTHER", nt.OTHER))}},
           "tasks": per_task, "price": price,
           "gate": {"e47_run_commit": E47_RUN_COMMIT, "e47_files_commit": E47_FILES_COMMIT, "tolerance": GATE_TOL,
                    "e49_run_commit": E49_RUN_COMMIT, "e49_files_commit": E49_FILES_COMMIT, "o_tolerance": O_TOL}}
    if dec["passes"]:
        res["question"] = {"label": LABEL, "command": CANDIDATE_COMMAND, "per_truth": price["per_truth"],
                           "gains": price["gains"], "losses": price["losses"],
                           "four_truth_mean": float(d["board4"].mean()), "test_like_mean": float(d["board_tl"].mean())}
    else:
        res["question"] = None
    return res


def small_g_range(recs, small_groups: dict[str, list[str]]) -> dict:
    """Ruling 1 (second set): the upper end of the test's range, the J4 merge of the small-g test-like tasks over their
    truths at equal weight, the arm against O under main (mean over seeds of D_s); at or below 0 the statement is 'no
    gain at the test'. The per-task figures are reported beside it and decide nothing."""
    if not small_groups:
        return {"gain_on_small_g_test_like": None, "no_gain_at_the_test": True, "small_g_per_task": {}}
    D = merged_D(recs, small_groups, ARM, O, "board")
    per = {n: float(np.mean(series(recs[n], ARM, "board") - series(recs[n], O, "board")))
           for names in small_groups.values() for n in names}
    return {"gain_on_small_g_test_like": _stat(D), "no_gain_at_the_test": bool(D.mean() <= 0),
            "small_g_per_task": per}


def _median(vals) -> float | None:
    v = [float(x) for x in vals if x is not None]
    return float(np.median(v)) if v else None


def text(res: dict) -> str:
    if res["void"]:
        return (f"Row E55 (run commit {res['commit']}): VOID: " + "; ".join(res["why"]) + "\n" + VOID_TEXT + "\n")
    m = res["merged"]
    L = [f"Row E55 (run commit {res['commit']}); {ROW_FIRST_LINE}.",
         f"The pre-run gate passed on all nine tasks (the floor and C0 under main within {res['gate']['tolerance']} of "
         f"E47, run at {res['gate']['e47_run_commit']}); O minus C0 reproduced E49's per-task ten-seed mean within "
         f"{res['gate']['o_tolerance']} on every task; the gate record of {GATE_COMMIT[:12]} was reproduced at the run "
         f"commit within {G_TOL} in every g.",
         "The pass rule, fixed in row E55 before any output, build main, against O (mean, SE, seeds ahead of 10); all "
         "eight clauses must hold:"]
    for k, v in res["clauses"].items():
        L.append(f"  {k} {'holds' if v else 'fails'}: {res['clause_text'][k]}")
    L.append(f"  four-truth board {m['board4']['mean']:+.3f} (SE {m['board4']['se']:.3f}, {m['board4']['ahead']} "
             f"ahead); five-truth board {m['board5']['mean']:+.3f} (SE {m['board5']['se']:.3f}, "
             f"{m['board5']['ahead']} ahead); four-truth gated shape_growth {m['sg4']['mean']:+.4f}")
    L.append(f"  test-like board {m['board_tl']['mean']:+.3f}, gated shape_growth {m['sg_tl']['mean']:+.4f}")
    for c, v in m["S_vs"].items():
        L.append(f"  against {TRY_IDS[c]} {c}: ungated S {v['mean']:+.4f} (SE {v['se']:.4f}, {v['ahead']} ahead)")
    L.append("  groups (five-truth / four-truth): " + "; ".join(
        f"{g} {m['groups5'][g]['mean']:+.4f} / {m['groups4'][g]['mean']:+.4f}" for g in OTHER_GROUPS))
    p7 = res["p7"]
    L.append(f"  P7 tasks with gate |median g| >= {G_MIN}: {', '.join(p7['tasks_at_or_above']) or 'none'}; board "
             + ("none" if m["board4_g"] is None else f"{m['board4_g']['mean']:+.3f}")
             + f"; reported beside, the tasks below ({', '.join(p7['tasks_below']) or 'none'}): "
             + ("none" if p7["board4_below"] is None else f"{p7['board4_below']['mean']:+.3f}"))
    if res["unreadable_because_refused"]:
        L.append("  the arm is unreadable (a refusal of " + ", ".join(res["unreadable_because_refused"])
                 + " on some task-seed; O's file stood in): STOP, not void; the three tries stay counted")
    tp = res["test_pair"]
    L.append("The recorded gate at the test pair, deciding nothing: " + tp["text"]
             + (" No gain at the test." if tp.get("no_gain_at_the_test", True) else
                f" The upper end of that range, the J4 merge of those tasks against O: "
                f"{tp['gain_on_small_g_test_like']['mean']:+.3f}.")
             + ("" if not tp.get("small_g_per_task") else " Per task, deciding nothing: "
                + ", ".join(f"{k} {v:+.3f}" for k, v in tp["small_g_per_task"].items()) + "."))
    rep = res["reported_deciding_nothing"]
    L.append(f"Reported, deciding nothing: prefix board five-truth {rep['prefix_board5']['mean']:+.3f}, four-truth "
             f"{rep['prefix_board4']['mean']:+.3f}; per truth "
             + ", ".join(f"{k} {v:+.3f}" for k, v in res["price"]["per_truth"].items())
             + "; refused task-seeds " + ", ".join(f"{k} {v}" for k, v in res["refused_task_seeds"].items()))
    if res["passes"]:
        q = res["question"]
        L.append(f"PASS: one priced naming question to the captain, labelled '{q['label']}' (selected on the nine "
                 "tasks: the composition was chosen after reading E49 on them; a diagnostic, never 'kept').")
        L.append("  gains per truth: " + (", ".join(f"{k} {v:+.3f}" for k, v in q["gains"].items()) or "none")
                 + "; losses per truth: " + (", ".join(f"{k} {v:+.3f}" for k, v in q["losses"].items()) or "none")
                 + f"; four-truth mean {q['four_truth_mean']:+.3f}; test-like mean {q['test_like_mean']:+.3f}")
        L.append(f"  file, only on the captain's naming: {q['command']}")
        L.append("  at the test only as a range from zero (the test pair's median g is below P7's "
                 f"{G_MIN}).")
    else:
        L.append("STOP: no pass, no file.")
    L.append("Nothing here shows E7.75's shape, the live gate, build, metric seed or LAPACK signs.")
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
            "shape_interp_sha256": _sha_file(Path(si.__file__))}


def cmd_read(a) -> None:
    out = Path(a.out)
    files = [out / "reading.json", out / "reading.txt"]
    existing = [str(f) for f in files if f.exists()]
    if existing:
        sys.exit(f"refusing to overwrite: {', '.join(existing)}")
    _own_code()
    if not has_row():
        sys.exit("E55 has no row in experiments/harness/TRIES.md")
    runs = Path(a.runs)
    why, gate = gate0_problems(a.gate_record, a.gate_sha256, _recompute_dir(runs), a.commit)
    why = [*why, *gate_problems(runs, a.commit), *commit_problems(a.commit)]
    if repo_state(REPO).unclean != ():
        why.append("the tree is not clean at the read")
    res = {**reading(load_runs(runs), why, a.commit, gate), "reading_code": provenance(),
           "gate_record": {"commit": GATE_COMMIT, "gate_json_sha256": a.gate_sha256, "folder": str(a.gate_record)}}
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
        if name == "run":
            _gate_args(p)
    r = sub.add_parser("read")
    r.add_argument("--commit", required=True)
    r.add_argument("--runs", default=str(RUNS))
    r.add_argument("--out", default=str(RESULTS))
    _gate_args(r)
    a = ap.parse_args(argv)
    {"gate": cmd_gate, "run": cmd_run, "read": cmd_read}[a.cmd](a)


def _gate_args(p) -> None:
    p.add_argument("--gate-record", default=str(GATE_RECORD),
                   help=f"the target-free gate's run folder at {GATE_COMMIT[:12]} (pair files, gate.json)")
    p.add_argument("--gate-sha256", default=GATE_SHA256, help="the sha256 of that run's gate.json (the row's)")


if __name__ == "__main__":
    main()
