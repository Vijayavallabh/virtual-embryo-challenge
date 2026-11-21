#!/usr/bin/env python
"""Row E52 of experiments/harness/TRIES.md (as specified; one try, placed by the main session before any run): the
copy-rows file at the team's pre-registered size, on the nine practice tasks. The captain's order of 2026-11-21
(docs/PLAN.md decision log 2026-11-21 (v)); specified by a fresh analyst blind to test scores and amended by a cold
challenger blind to test scores; the main session is recused on this lever. Cut from e50_copy_rows_read.py and fixed
before any run. No number here comes from a test score.

Arm. The harness candidate ARM, fgw_copyrows_sized: per task and scoring seed s, fgw_copyrows exactly (rows, order,
count, mixture expression, streams [build seed, 1], [build seed, 4], [build seed, 51]), its float32 coordinates then
float32(rescale_to_rms(float64 of the bake floor copy subsample(ref, 0.1, s), T)), scaled about the copy's own
centroid, T = harness.target_rms on the task's own role stages only (fgw_shape0x's target). Controls, rescored and not
counted: the floor, fgw_copyrows (C, E50's arm) and fgw_shape0x (P).

Path. harness.run_task on each task of nine_tasks.all_tasks(), seeds 0-9, bake protocol, build seed BUILD_SEED_OFFSET
+ seed, under build "main" and under "prefix". A refusal of the arm or of fgw_copyrows in E50's three categories is
caught and recorded; fgw_shape0x's file of that seed stands in so the run goes on, and that task-seed is left out of
every reading. Any other error (a KeyError included) is raised: the seed is not written and the run is void.

Steps (each refuses unless row E52 is in TRIES.md, the tree is clean and HEAD is the run commit):
  gate --task NAME --commit SHA  the floor, fgw_shape0x and fgw_copyrows under both builds, ten seeds: the floor scores
                                 50 in every seed; fgw_shape0x under main is within 0.1 board points of E47's
                                 committed ten-seed mean and 0.5 in any seed; fgw_copyrows minus fgw_shape0x has each
                                 build's ten-seed mean within 0.01 of E50's committed reading (iii) primary board for
                                 the task. Exit 3 on a failure (stop, report, count nothing). A gate file of this commit
                                 that exists is read again, not re-run. The run step also needs every D_s (the nine
                                 tasks' paired difference merged over the five truths, per seed and build) within 0.05
                                 of E50's.
  run --task NAME --commit SHA   the arm and controls under both builds, seed by seed (resumable: a seed whose file
                                 exists is skipped); refused unless the gate holds on all nine tasks at the commit
  read --commit SHA              the reading into experiments/coords/results/e52_copy_rows_sized/ (reading.json,
                                 reading.txt), never overwriting
The script prints only "<task> seed <s> done", "wrote <path>" and "gate <task> pass|fail". It reads no test board, test
file, test anchor or size file, and prints no count.

Readings, fixed in row E52 before any output (applied by `reading`); raw values from the in-memory builds:
  (i)   on every task-seed where the arm builds, under main and under prefix: its X bytes equal fgw_copyrows's; its
        coordinates equal float32(rescale_to_rms(float64 floor copy, T)) bit for bit, with T computed by this read on
        the task's own role stages, and their RMS is within 1e-6 relative of T; its prefix parity on the matched rows
        at the scored count and metric seed equals the floor's; its d2_shape is within 1e-6 of the floor's; its
        canonical frame is within 1e-5 of the floor's up to a proper flip; its occupancy_dice equals the floor's except
        by at most 4v over (the arm's plus the truth's occupied voxels), plus 1e-4, v the matched rows whose voxel
        changes, maximised over the four proper flips; its de_score, de_direction, mmd_u and variogram equal
        fgw_copyrows's exactly; any failure voids the construction and no test file is built;
  (ii)  refusals in E50's three categories, equal to fgw_copyrows's at every task-seed;
  (iii) paired arm minus fgw_shape0x and arm minus fgw_copyrows, per task, merged over the five truths and over the
        test-like four, both builds: board, each group, each metric, seeds ahead, SE, mean gain and mean loss, worst
        and best task;
  (iv)  per task, the arm's neighborhood_mmd skill and its shape group gated and ungated.
Reported, deciding nothing: scale_log_ratio against fgw_shape0x; neighborhood_mmd against fgw_copyrows with the number
of cells whose 15-neighbour set differs; prefix parity for every count from 4 to n; sliced_wasserstein against the
floor.
Outcome: if (i) holds, the test file is built under the captain's order; (iii) and (iv) decide nothing. Void: a missing
task, seed or arm; an arm error outside the listed refusals (a KeyError included); a refusal that differs from
fgw_copyrows's or between the builds; mixed commits or an unclean tree; a gate miss; files differing between the two
builds' scoring; a test board.
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

ROW_IDS = ("E52",)
TRIES = REPO / "experiments" / "harness" / "TRIES.md"
ARM = "fgw_copyrows_sized"
C_ARM, P_ARM = "fgw_copyrows", "fgw_shape0x"  # the copy-rows file C and the exact-position file P
CONTROLS = (P_ARM, C_ARM)  # P first: its file stands in for a refused seed
AGAINST = (P_ARM, C_ARM)  # reading (iii): arm minus P, arm minus C
BUILDS = ("main", "prefix")
TEST_LIKE = ("half", "late", "new1", "new2")
GROUP_SHAPE = "shape_growth"
GROUP_NAMES = tuple(name for name, _, _ in GROUPS["T2"])
METRICS = tuple(RANKED_METRICS["T2"])
NBHD = "neighborhood_mmd"
GATE_HALF = 0.5  # the assumed gate: shape * min(1, neighborhood skill / 0.5)
EXPRESSION_TERMS = ("de_score", "de_direction", "mmd_u", "variogram")  # equal to fgw_copyrows's exactly
RMS_TOL, D2_TOL, FRAME_TOL, DICE_SLACK = 1e-6, 1e-6, 1e-5, 1e-4
DICE_V_FACTOR = 4
PARITY_FROM = 4  # three rows have rank two: their parity is arbitrary
NN_K = 15  # reported: cells whose 15-neighbour set differs from fgw_copyrows's
N_GRID, EXTENT = 16, 3.0  # the scorer's occupancy_dice grid (common.shape_metrics defaults)
E47_DIR = REPO / "experiments" / "harness" / "results" / "e47_floor_exact"
E47_RUN_COMMIT, E47_FILES_COMMIT = "9967cf1", "d5ee9b5"
GATE_MEAN, GATE_SEED = 0.1, 0.5  # board points against E47: a task's ten-seed mean, any one seed
E50_READING = REPO / "experiments" / "coords" / "results" / "e50_copy_rows" / "reading.json"
E50_RUN_COMMIT = "5fcf6b9ad187e2c5377900426f6dd749c4a1f9d9"  # row E50's run commit (its reading's commit field)
E50_FILES_COMMIT = "59ad7fa"  # the commit that holds E50's reading
E50_MEAN_TOL, E50_DS_TOL = 0.01, 0.05  # fgw_copyrows minus fgw_shape0x against E50's (iii) primary board
FLOOR_SCORE = 50.0
BASES = ("b1430aa",)  # main when the branch campaign/embryo-copy-rows-size was cut
GUARD_FLAG = "--copy-rows-size-test-ok"
SCRIPT = "experiments/coords/e52_copy_rows_sized_read.py"
RUNS = REPO / "runs" / "e52"
RESULTS = REPO / "experiments" / "coords" / "results" / "e52_copy_rows_sized"
REFUSALS = (("later bracket's share", "dropped_over_share"), ("cells outside [", "count_outside_limits"),
            ("fewer than three rows in the copy", "copy_too_small"))  # E50's three categories


# ------------------------------------------------------------------------------------------------ checks
def _git(*args: str) -> subprocess.CompletedProcess:
    return subprocess.run(["git", "-C", str(REPO), *args], capture_output=True, text=True,
                          env={**os.environ, "GIT_OPTIONAL_LOCKS": "0"})


def has_row(path: Path | None = None) -> bool:
    rows = n12.tries_rows(TRIES if path is None else path)
    return any(r in rows for r in ROW_IDS)


def commit_problems(commit: str) -> list[str]:
    """Why this checkout may not run or read row E52 at `commit` (empty: it may)."""
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
    if _git("diff", "--quiet", E50_FILES_COMMIT, "HEAD", "--", str(E50_READING.relative_to(REPO))).returncode != 0:
        why.append(f"E50's committed reading differs from {E50_FILES_COMMIT}'s")
    return why


def _own_code() -> None:
    src = (REPO / "src").resolve()
    for mod in (ev, harness, tb):
        if src not in Path(mod.__file__).resolve().parents:
            sys.exit(f"{mod.__name__} is loaded from {mod.__file__}, not from this checkout's src/")


def _refuse_unless_ready(commit: str) -> None:
    _own_code()
    if not has_row():
        sys.exit("E52 has no row in experiments/harness/TRIES.md: a row is committed before its run")
    if repo_state(REPO).unclean != ():
        sys.exit("the tree is not clean: row E52 runs on a committed tree")
    why = commit_problems(commit)
    if why:
        sys.exit("refused: " + "; ".join(why))


def _task(name: str):
    tasks = nt.all_tasks()
    if name not in tasks:
        sys.exit(f"unknown task {name}; known: {sorted(tasks)}")
    t = tasks[name]
    if getattr(t, "target_filter", None) is not None:
        sys.exit(f"{name}: a task with a target filter; E52 does not run on one")
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


def same_x(A, B) -> bool:
    """Whether two expression matrices are the same bytes as float32 (shape and every byte)."""
    a = np.ascontiguousarray(np.asarray(A), dtype=np.float32)
    b = np.ascontiguousarray(np.asarray(B), dtype=np.float32)
    return a.shape == b.shape and a.tobytes() == b.tobytes()


def refusal_category(message: str) -> str | None:
    """The category of a refusal (no value kept), or None for an error that is not a refusal."""
    for needle, cat in REFUSALS:
        if needle in message:
            return cat
    return None


# ------------------------------------------------------------------------------------------------ raw quantities
def read_target(view, board, inputs) -> float:
    """T as this read computes it, apart from the arm's code path: harness.target_rms on a context holding only the
    task's own role stages (fgw_shape0x's target; never the reference cells run_task adds for the copy-rows arms)."""
    ctx = harness.BuildContext(task=view, board=board, inputs={r: inputs[r] for r in view.roles}, seed=0)
    return float(harness.target_rms(ctx))


def expected_coords(floor_coords, T: float) -> np.ndarray:
    """float32(rescale_to_rms(float64 of the floor copy, T)), the arithmetic written out (row E52, construction 6)."""
    C64 = np.asarray(floor_coords, dtype=np.float32).astype(np.float64)
    mu = C64.mean(0)
    return (mu + (C64 - mu) * (T / tb.rms_radius(C64))).astype(np.float32)


def _occ(Z: np.ndarray) -> np.ndarray:
    """The scorer's voxel index of each row (common.shape_metrics.occupancy_dice's occ, kept per row)."""
    return np.floor((np.clip(Z, -EXTENT, EXTENT - 1e-9) + EXTENT) / (2 * EXTENT) * N_GRID).astype(int)


def frame_and_voxels(arm_coords, floor_coords, truth_coords, seed: int) -> dict:
    """Under the scorer build in force: the canonical frames of the arm's and the floor's matched rows (the same rows:
    _match_n depends on the counts only), the proper flip g that aligns them and the largest frame deviation after it;
    v, the most matched rows whose voxel differs between Z_arm * f and Z_floor * g * f over the four proper flips f;
    the smallest (arm's plus truth's occupied voxels) over the flips; and the raw occupancy_dice of each file."""
    sm = ev._shape_metrics()
    A, B = ev.match_rows(arm_coords, truth_coords, seed)
    F, _ = ev.match_rows(floor_coords, truth_coords, seed)
    Za, Zf, Zb = sm._canonicalise(A)[0], sm._canonicalise(F)[0], sm._canonicalise(B)[0]
    flips = [np.asarray(f, float) for f in sm._PROPER_FLIPS]
    devs = [float(np.abs(Za - Zf * g).max()) for g in flips]
    g = flips[int(np.argmin(devs))]
    ob = {tuple(x) for x in _occ(Zb)}
    v, den = 0, None
    for f in flips:
        va, vf = _occ(Za * f), _occ(Zf * g * f)
        v = max(v, int((va != vf).any(1).sum()))
        d = len({tuple(x) for x in va}) + len(ob)
        den = d if den is None else min(den, d)
    return {"frame_dev": min(devs), "v": v, "occupied": den,
            "dice_arm": float(sm.occupancy_dice(arm_coords, truth_coords, seed=seed)[0]),
            "dice_floor": float(sm.occupancy_dice(floor_coords, truth_coords, seed=seed)[0])}


def dice_bound(v: int, occupied: int) -> float:
    return DICE_V_FACTOR * v / max(occupied, 1) + DICE_SLACK


def parity_flips(arm_coords, floor_coords, seed: int, n_max: int | None = None) -> list[int]:
    """The counts n from PARITY_FROM to the arm's (or n_max) at which the SVD parity of the rows _match_n takes at
    metric seed `seed` differs between the arm and the floor copy (count-free: the truth's rows play no part)."""
    n_max = len(arm_coords) if n_max is None else n_max
    out = []
    for n in range(PARITY_FROM, n_max + 1):
        z = np.zeros((n, 3))
        if ev.svd_parity(ev.match_rows(arm_coords, z, seed)[0]) != ev.svd_parity(ev.match_rows(floor_coords, z,
                                                                                               seed)[0]):
            out.append(n)
    return out


def nn_sets_differ(A, B, k: int = NN_K) -> int:
    """Reported: the cells whose k-neighbour set (by index, self left out) differs between two clouds of one row
    order."""
    from scipy.spatial import cKDTree

    def sets(C):
        C = np.asarray(C, dtype=np.float64)
        _, idx = cKDTree(C).query(C, k=min(k + 1, len(C)))
        return [frozenset([j for j in row if j != i][:k]) for i, row in enumerate(idx)]

    return int(sum(a != b for a, b in zip(sets(A), sets(B), strict=True)))


def raw_facts(arm, crow, floor, truth, T: float, seed: int, build: str, metrics: dict, with_reported: bool) -> dict:
    """Reading (i)'s raw quantities for one task-seed under one build, from the built arrays (never panel values,
    except the four expression terms, compared as the scorer reports them)."""
    rms = tb.rms_radius(np.asarray(arm.coords, dtype=np.float32).astype(np.float64))
    with ev.scorer_build(build):
        fv = frame_and_voxels(arm.coords, floor.coords, truth.coords, seed)
    sm = ev._shape_metrics()
    A, _ = ev.match_rows(arm.coords, truth.coords, seed)
    F, _ = ev.match_rows(floor.coords, truth.coords, seed)
    out = {"x_equal_c": same_x(arm.X, crow.X),
           "coords_equal_expected": same_bits(arm.coords, expected_coords(floor.coords, T)),
           "rms_rel": abs(rms / T - 1.0), "target_read": T,
           "parity_equal_scored": ev.svd_parity(A) == ev.svd_parity(F),
           "d2_diff": abs(float(sm.d2_distance(arm.coords, truth.coords, seed=seed))
                          - float(sm.d2_distance(floor.coords, truth.coords, seed=seed))),
           **fv, "expression_equal_c": {m: metrics[ARM].get(m) == metrics[C_ARM].get(m) for m in EXPRESSION_TERMS}}
    if with_reported:
        out["reported"] = {"parity_flips_all_counts": len(parity_flips(arm.coords, floor.coords, seed)),
                           "nn_sets_differ_vs_c": nn_sets_differ(arm.coords, crow.coords)}
    return out


def i_failures(f: dict) -> list[str]:
    """Which of reading (i)'s checks one build's raw facts fail (empty: all hold)."""
    bad = []
    if not f.get("x_equal_c"):
        bad.append("x_bytes")
    if not f.get("coords_equal_expected"):
        bad.append("coordinates")
    if not (f.get("rms_rel") is not None and f["rms_rel"] <= RMS_TOL):
        bad.append("rms")
    if not f.get("parity_equal_scored"):
        bad.append("prefix_parity")
    if not (f.get("d2_diff") is not None and f["d2_diff"] <= D2_TOL):
        bad.append("d2_shape")
    if not (f.get("frame_dev") is not None and f["frame_dev"] <= FRAME_TOL):
        bad.append("frame")
    if not (f.get("dice_arm") is not None and f.get("dice_floor") is not None
            and abs(f["dice_arm"] - f["dice_floor"]) <= dice_bound(f.get("v", 0), f.get("occupied", 0))):
        bad.append("occupancy_dice")
    expr = f.get("expression_equal_c") or {}
    bad += [m for m in EXPRESSION_TERMS if expr.get(m) is not True]
    return bad


# ------------------------------------------------------------------------------------------------ the arms
def file_digest(stage) -> str:
    X = np.ascontiguousarray(stage.X, dtype=np.float32)
    C = np.ascontiguousarray(stage.coords, dtype=np.float32)
    return hashlib.sha256(np.asarray(X.shape, np.int64).tobytes() + X.tobytes() + C.tobytes()).hexdigest()


@contextlib.contextmanager
def captured():
    """While active, record what harness.score_seeds is handed (task, board, target, reference stage, role inputs)."""
    seen: dict = {"stage": {}, "file": {}, "refused": {}}
    orig = harness.score_seeds

    def score_seeds(t, board, target, ref, inputs, candidates, *a, **kw):
        seen.update(task=t, board=board, target=target, ref=ref, inputs=inputs)
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
    """A control unchanged; its file is kept (fgw_shape0x's also stands in for a refused seed)."""
    def build(ctx):
        pred = base.build(ctx)
        _keep(seen, base.name, ctx.seed - harness.BUILD_SEED_OFFSET, pred)
        return pred

    return dataclasses.replace(base, build=build)


def guarded(base, seen):
    """The arm or fgw_copyrows: its file, or on a refusal in E50's categories fgw_shape0x's file of the seed, with the
    category recorded. Any other error (a KeyError included) is raised and voids the run."""
    def build(ctx):
        s = ctx.seed - harness.BUILD_SEED_OFFSET
        if _is_test_board(seen["board"]):
            raise RuntimeError("E52 met a test board")
        try:
            pred, why = base.build(ctx), None
        except ValueError as e:
            why = refusal_category(str(e))
            if why is None:
                raise
            pred = seen["stage"][(P_ARM, s)]
        seen["refused"][(base.name, s)] = why
        _keep(seen, base.name, s, pred)
        return pred

    return dataclasses.replace(base, build=build)


def candidates(seen) -> list:
    return [harness.CANDIDATES[harness.FLOOR], observed(harness.CANDIDATES[P_ARM], seen),
            guarded(harness.CANDIDATES[C_ARM], seen), guarded(harness.CANDIDATES[ARM], seen)]


def one_seed(t, s: int, data: Path) -> dict:
    """Both builds of one task-seed: run_task's results and E52's raw per-seed facts."""
    out: dict = {"builds": {}, "facts": {"per_build": {}}}
    names = (*CONTROLS, ARM)
    files, refused = {}, {}
    view = harness.TaskView.of(t)
    for k, build in enumerate(BUILDS):
        with captured() as seen:
            res = harness.run_task(t, candidates(seen), seeds=[s], data_dir=data, build=build, log=None)
        board = seen["board"]
        if _is_test_board(board):
            raise RuntimeError("E52 met a test board")
        out["builds"][build] = res
        files[build] = {n: seen["file"][(n, s)] for n in names}
        refused[build] = {n: seen["refused"].get((n, s)) for n in (C_ARM, ARM)}
        if k == 0:
            out["facts"]["board_split"] = str(board.split)
        if refused[build][ARM] is not None:
            continue
        floor = ev.subsample(seen["ref"], ev.SUBSAMPLE_FRAC, s)
        truth = ev.working_copies(seen["target"], seen["ref"], s, ev.SUBSAMPLE_FRAC, harness.PROTOCOL)[0]
        T = read_target(view, harness._candidate_board(board), seen["inputs"])
        out["facts"]["per_build"][build] = raw_facts(seen["stage"][(ARM, s)], seen["stage"][(C_ARM, s)], floor, truth,
                                                     T, s, build, res["rows"][0]["metrics"], with_reported=k == 0)
    out["facts"]["refused"] = refused[BUILDS[0]][ARM]
    out["facts"]["refused_c"] = refused[BUILDS[0]][C_ARM]
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


def e50_primary(path: Path | None = None) -> tuple[dict | None, list[str]]:
    """E50's committed reading (iii) primary (fgw_copyrows minus fgw_shape0x) per build, and why it cannot serve."""
    p = E50_READING if path is None else path
    if not Path(p).exists():
        return None, ["E50's reading is missing"]
    r = json.loads(Path(p).read_text())
    why = []
    if r.get("void") is not False or r.get("commit") != E50_RUN_COMMIT or r.get("arm") != C_ARM:
        why.append("E50's reading is void, from another commit or of another arm")
    prim = r.get("iii", {}).get("primary", {})
    if prim.get("against") != P_ARM or any(b not in prim for b in BUILDS):
        why.append("E50's reading holds no primary table against fgw_shape0x under both builds")
        return None, why
    return {b: prim[b]["board"] for b in BUILDS}, why


def gate_check(name: str, key: str, scores: dict, old: dict | None = None, why_old: list[str] | None = None,
               e50: dict | None = None, why_e50: list[str] | None = None) -> dict:
    """The pre-run gate on one task. scores: build -> name -> per-seed board scores (ten seeds). The floor at 50 in
    every seed under both builds; fgw_shape0x under main within GATE_MEAN of E47's ten-seed mean and GATE_SEED in every
    seed; fgw_copyrows minus fgw_shape0x, per build, within E50_MEAN_TOL of E50's (iii) primary per-task mean. The
    per-seed differences are kept for the nine-task D_s check (gate_problems)."""
    if old is None and why_old is None:
        old, why_old = e47_scores(key)
    if e50 is None and why_e50 is None:
        e50, why_e50 = e50_primary()
    why = [*(why_old or []), *(why_e50 or [])]
    for b in BUILDS:
        floor = scores.get(b, {}).get(harness.FLOOR)
        if floor is None or len(floor) != len(SEEDS) or not np.all(np.abs(np.asarray(floor) - FLOOR_SCORE) <= 1e-9):
            why.append(f"{key}: the floor does not score {FLOOR_SCORE:g} in every seed under {b}")
    worst = {"e47_mean": None, "e47_seed": None, "e50_mean": {}}
    if old is not None:
        o, n = old.get(P_ARM), scores.get("main", {}).get(P_ARM)
        if o is None or n is None or np.shape(o) != np.shape(n):
            why.append(f"{key}: {P_ARM} missing or of another length")
        else:
            o, n = np.asarray(o, float), np.asarray(n, float)
            dm, ds = abs(float(n.mean() - o.mean())), float(np.abs(n - o).max())
            worst.update(e47_mean=dm, e47_seed=ds)
            if not dm <= GATE_MEAN:
                why.append(f"{key}: {P_ARM}'s ten-seed mean is {dm:.4f} off E47's (tolerance {GATE_MEAN})")
            if not ds <= GATE_SEED:
                why.append(f"{key}: {P_ARM} is {ds:.4f} off E47's in a seed (tolerance {GATE_SEED})")
    diffs = {}
    for b in BUILDS:
        c, p = scores.get(b, {}).get(C_ARM), scores.get(b, {}).get(P_ARM)
        if c is None or p is None or len(c) != len(SEEDS) or len(p) != len(SEEDS):
            why.append(f"{key}: {C_ARM} or {P_ARM} missing under {b}")
            continue
        d = np.asarray(c, float) - np.asarray(p, float)
        diffs[b] = d.tolist()
        want = None if e50 is None else e50.get(b, {}).get("per_task", {}).get(name, {}).get("mean")
        if want is None:
            why.append(f"{key}: E50 holds no {b} per-task mean for {name}")
            continue
        dm = abs(float(d.mean()) - float(want))
        worst["e50_mean"][b] = dm
        if not dm <= E50_MEAN_TOL:
            why.append(f"{key}: {C_ARM} minus {P_ARM} under {b} is {dm:.4f} off E50's per-task mean "
                       f"(tolerance {E50_MEAN_TOL})")
    return {"task": name, "key": key, "e47": {"run_commit": E47_RUN_COMMIT, "files_commit": E47_FILES_COMMIT},
            "e50": {"run_commit": E50_RUN_COMMIT, "files_commit": E50_FILES_COMMIT},
            "tolerance": {"e47_mean": GATE_MEAN, "e47_seed": GATE_SEED, "e50_mean": E50_MEAN_TOL,
                          "e50_d_s": E50_DS_TOL},
            "worst": worst, "diffs": diffs, "why": why, "pass": not why}


def d_s_check(gates: dict[str, dict], e50: dict | None, groups: dict[str, list[str]]) -> list[str]:
    """fgw_copyrows minus fgw_shape0x merged over the five truths per seed (D_s), per build, against E50's D_s."""
    if e50 is None:
        return ["E50's reading cannot serve for the D_s check"]
    why = []
    for b in BUILDS:
        if any(b not in g.get("diffs", {}) for g in gates.values()):
            why.append(f"a gate file lacks the {b} differences")
            continue
        D = e47.merged_seeds({n: np.asarray(g["diffs"][b], float) for n, g in gates.items()}, groups)
        want = e50.get(b, {}).get("D_s")
        if want is None or len(want) != len(D):
            why.append(f"E50 holds no {b} D_s of ten seeds")
            continue
        off = float(np.abs(D - np.asarray(want, float)).max())
        if not off <= E50_DS_TOL:
            why.append(f"{b}: a D_s is {off:.4f} off E50's (tolerance {E50_DS_TOL})")
    return why


def cmd_gate(a) -> None:
    _refuse_unless_ready(a.commit)
    t = _task(a.task)
    path = _gate_path(Path(a.runs), a.task)
    if path.exists():  # resumable: a gate file of this commit is read again, not re-run
        rec = json.loads(path.read_text())
        if rec.get("commit") != a.commit:
            sys.exit(f"{path} exists from another commit")
    else:
        cands = [harness.CANDIDATES[n] for n in (harness.FLOOR, P_ARM, C_ARM)]
        res = {b: harness.run_task(t, cands, seeds=SEEDS, data_dir=Path(a.data), build=b, log=None) for b in BUILDS}
        if any(str(r["board"].get("split", "")).startswith("test") for r in res.values()):
            sys.exit("E52 met a test board")
        rec = {**gate_check(a.task, t.key, {b: e47.per_seed(r)[0] for b, r in res.items()}),
               "commit": res["main"]["provenance"]["commit"],
               "tree_clean": all(r["provenance"]["tree_clean"] and r["provenance"]["commit"] == a.commit
                                 for r in res.values())}
        rec["pass"] = bool(rec["pass"] and rec["commit"] == a.commit and rec["tree_clean"])
        path.parent.mkdir(parents=True, exist_ok=True)
        tmp = path.with_suffix(".tmp")
        tmp.write_text(json.dumps(rec, indent=1, default=float))
        os.replace(tmp, path)
    print(f"gate {a.task} {'pass' if rec['pass'] else 'fail'}", flush=True)
    sys.exit(0 if rec["pass"] else 3)


def gate_problems(runs: Path, commit: str) -> list[str]:
    why, gates = [], {}
    for name in nt.all_tasks():
        p = _gate_path(runs, name)
        if not p.exists():
            why.append(f"no gate file for {name}")
            continue
        g = json.loads(p.read_text())
        if not g.get("pass") or g.get("commit") != commit:
            why.append(f"the gate did not pass on {name} at the run commit")
        gates[name] = g
    if not why:
        e50, why_e50 = e50_primary()
        why += why_e50 + d_s_check(gates, e50, e47.truths(nt.all_tasks()))
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
            if f.get("refused") != f.get("refused_c"):
                why.append(f"{tag}: the arm's refusal differs from {C_ARM}'s")
            if f.get("refused") is None and any(b not in f.get("per_build", {}) for b in BUILDS):
                why.append(f"{tag}: the arm built but its raw facts are missing")
    return why


def _row(r: dict, build: str) -> dict:
    return r["builds"][build]["rows"][0]


def built(r: dict) -> bool:
    return r["facts"].get("refused") is None


def reading_i(recs: dict[str, list[dict]]) -> dict:
    """Reading (i) on every built task-seed under both builds, from the raw facts. Reported beside it, deciding
    nothing: counts with a prefix parity flip (4 to n), cells whose 15-neighbour set differs from fgw_copyrows's,
    sliced_wasserstein unequal to the floor's (panel values)."""
    failed, n_built = [], 0
    flips, nn, sw = [], [], []
    for task, rs in recs.items():
        for r in rs:
            if not built(r):
                continue
            n_built += 1
            for b in BUILDS:
                f = r["facts"]["per_build"][b]
                for what in i_failures(f):
                    failed.append({"task": task, "seed": r["seed"], "build": b, "what": what})
                rep = f.get("reported")
                if rep is not None:
                    if rep.get("parity_flips_all_counts"):
                        flips.append({"task": task, "seed": r["seed"], "counts": rep["parity_flips_all_counts"]})
                    nn.append({"task": task, "seed": r["seed"], "cells": rep.get("nn_sets_differ_vs_c")})
                m = _row(r, b)["metrics"]
                if "sliced_wasserstein" in m[ARM] and m[ARM]["sliced_wasserstein"] != m[harness.FLOOR].get(
                        "sliced_wasserstein"):
                    sw.append({"task": task, "seed": r["seed"], "build": b})
    return {"holds": bool(n_built > 0 and not failed), "built_task_seeds": n_built, "failed": failed,
            "reported": {"parity_flip_counts": flips, "nn_sets_differ_vs_c": nn, "sliced_wasserstein_unequal": sw}}


def reading_ii(recs: dict[str, list[dict]]) -> dict:
    per = {task: {str(r["seed"]): r["facts"]["refused"] for r in rs if not built(r)} for task, rs in recs.items()}
    every = sorted(task for task, rs in recs.items() if rs and not any(built(r) for r in rs))
    return {"refused": {k: v for k, v in per.items() if v}, "refused_at_every_seed": every,
            "equal_to_c": all(r["facts"].get("refused") == r["facts"].get("refused_c")
                              for rs in recs.values() for r in rs)}


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


def paired(recs: dict[str, list[dict]], groups: dict[str, list[str]], a: str, b: str, what: str, build: str) -> dict:
    """Reading (iii) for one quantity over the tasks in `groups`: per task (mean, seeds ahead, seeds read), merged over
    the truths (a truth's value the mean over its tasks that have a built seed), D_s per seed, its SE, the mean gain
    over gain tasks, the mean loss over loss tasks, the worst and the best task."""
    tasks = {n for names in groups.values() for n in names}
    per_task = {}
    for task, rs in recs.items():
        if task not in tasks:
            continue
        d = {r["seed"]: value(_row(r, build), a, what) - value(_row(r, build), b, what) for r in rs if built(r)}
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


def reading_iii(recs: dict[str, list[dict]], groups: dict, test_like: dict, against: str) -> dict:
    whats = ("board", *(f"group:{g}" for g in GROUP_NAMES), *(f"skill:{m}" for m in METRICS))
    return {b: {w: {"five_truths": paired(recs, groups, ARM, against, w, b),
                    "test_like": paired(recs, test_like, ARM, against, w, b)} for w in whats} for b in BUILDS}


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
                                                            for r in rows]))}
    return out


def reported(recs: dict[str, list[dict]]) -> dict:
    """Deciding nothing, per task and build over the built seeds: the scale_log_ratio metric of the arm and of
    fgw_shape0x, and the neighborhood_mmd metric of the arm and of fgw_copyrows (panel values)."""
    out = {}
    for task, rs in recs.items():
        rows = [r for r in rs if built(r)]
        if not rows:
            continue
        out[task] = {b: {f"{m}:{n}": float(np.mean([_row(r, b)["metrics"][n][m] for r in rows]))
                         for m, ns in (("scale_log_ratio", (ARM, P_ARM)), (NBHD, (ARM, C_ARM))) for n in ns}
                     for b in BUILDS}
    return out


def reading(recs: dict[str, list[dict]], why_in: list[str], commit: str) -> dict:
    why = check_run(recs, why_in, commit)
    if why:
        return {"void": True, "why": why, "commit": commit, "arm": ARM}
    tasks = nt.all_tasks()
    groups = e47.truths(tasks)
    test_like = e47.truths({n: tasks[n] for n in TEST_LIKE})
    i = reading_i(recs)
    outcome = ("build: (i) holds; the test file is built under the captain's order ((iii) and (iv) decide nothing)"
               if i["holds"] else "no test file: (i) fails, which voids the construction")
    return {"void": False, "why": [], "commit": commit, "arm": ARM, "truths": groups, "test_like_truths": test_like,
            "i": i, "ii": reading_ii(recs),
            "iii": {tag: {"against": against, **reading_iii(recs, groups, test_like, against)}
                    for tag, against in (("vs_P", P_ARM), ("vs_C", C_ARM))},
            "iv": reading_iv(recs), "reported": reported(recs), "build_test_file": bool(i["holds"]),
            "outcome": outcome,
            "gate": {"e47_run_commit": E47_RUN_COMMIT, "e47_files_commit": E47_FILES_COMMIT,
                     "e50_run_commit": E50_RUN_COMMIT, "e50_files_commit": E50_FILES_COMMIT,
                     "tolerance": {"e47_mean": GATE_MEAN, "e47_seed": GATE_SEED, "e50_mean": E50_MEAN_TOL,
                                   "e50_d_s": E50_DS_TOL}}}


def _fmt(x, f="+.3f") -> str:
    return "n/a" if x is None else format(x, f)


def text(res: dict) -> str:
    if res["void"]:
        return f"Row E52 (run commit {res['commit']}): VOID: " + "; ".join(res["why"]) + "\n"
    i = res["i"]
    L = [f"Row E52 (run commit {res['commit']}), arm {res['arm']}; the pre-run gate passed on all nine tasks.",
         f"(i) on every built task-seed, both builds: {'HOLDS' if i['holds'] else 'FAILS'}"
         + ("" if i["holds"] else f" ({len(i['failed'])} failed checks: "
            + ", ".join(sorted({f['what'] for f in i['failed']})) + ")"),
         "    reported: task-seeds with a prefix parity flip at some count: "
         f"{len(i['reported']['parity_flip_counts'])}; sliced_wasserstein unequal to the floor's: "
         f"{len(i['reported']['sliced_wasserstein_unequal'])}",
         "(ii) refusals per task-seed: " + (json.dumps(res["ii"]["refused"]) if res["ii"]["refused"] else "none")
         + "; refused at every seed: " + (", ".join(res["ii"]["refused_at_every_seed"]) or "none")]
    for tag in ("vs_P", "vs_C"):
        t = res["iii"][tag]
        L.append(f"(iii) arm minus {t['against']} (mean, SE, seeds ahead), five truths / test-like four:")
        for b in BUILDS:
            for w, v in t[b].items():
                a, tl = v["five_truths"], v["test_like"]
                L.append(f"    {b} {w}: {_fmt(a['merged'], '+.4f')} (SE {_fmt(a['se'], '.4f')}, {a['seeds_ahead']} "
                         f"ahead; gain {a['gain_tasks']} {_fmt(a['mean_gain'], '+.4f')}, loss {a['loss_tasks']} "
                         f"{_fmt(a['mean_loss'], '+.4f')}; worst {a['worst_task']}, best {a['best_task']}) / "
                         f"{_fmt(tl['merged'], '+.4f')} (SE {_fmt(tl['se'], '.4f')}, {tl['seeds_ahead']} ahead)")
    L.append("(iv) per task (main): neighborhood_mmd skill, seeds where the assumed gate bites, shape gated / ungated:")
    for task, v in res["iv"].items():
        if v is None:
            L.append(f"    {task}: refused at every seed")
            continue
        m = v["main"]
        L.append(f"    {task}: {m['neighborhood_mmd_skill']:.3f}, {m['seeds_gate_bites']}, {m['shape_gated']:.4f} / "
                 f"{m['shape_ungated']:.4f}")
    L.append(f"Outcome: {res['outcome']}.")
    L.append("Nothing here shows that the live floor is subsample(E7.5 as released, 0.1, 0), that rows are scored as "
             "held, or that anchors and uploads share a scorer build and numeric library; nothing about E7.75; the "
             "practice size is the rule read on each task's own stages, not the test's size file.")
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
        sys.exit("E52 has no row in experiments/harness/TRIES.md")
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
