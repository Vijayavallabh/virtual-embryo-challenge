#!/usr/bin/env python
"""Row E59 of experiments/harness/TRIES.md: the gate of row R71's files T2embryo_s2x_ovbf.h5ad (BF_O) and
T2embryo_s2xs_ovbf.h5ad (BF_SO), the embryo's consensus-form files BF on O and on SO only (the captain's order of
2026-11-25). BF is EMBRYO_MOON2_design.md's P1 with every fix of its cold challenge (EMBRYO_MOON2_challenge.md): no step
depends on the cell count (fix 3), rule 7's with-and-without controls for the five extra specimens and for the flip
step (fix 4, clause (e)), the E8.0 truth's clause (fix 5), G0's idle case and reports (fix 6), the pool moved by the
same field (fix 7), one ruling flag (fix 8), O, SO and the two controls only (fix 9). Its structure follows
experiments/coords/e55_overlay_stretch_read.py and e58_mix_lean_read.py. No number here comes from a test score. It
moves the overlay's form, the occupancy term and the z extent: the main session is recused on the embryo's shape forms
and occupancy term and reads the outcome line alone (`python3 -I -B experiments/coords/bf_row.py outcome ...`).

Arms (vec.shape_interp, row R71's forms; the prior is loaded here, outside the harness's read guard, and handed to the
build; on every task, as at the test pair, it is the stages of BAYES_PRIOR_STAGES that do not lie strictly between the
task's brackets, so never the truth, J2), four counted tries when the run step starts:
  BF_O   fgw_shape0x_ov_bf         O's build with shape_form overlay_bayes (the file T2embryo_s2x_ovbf.h5ad's path)
  BF_SO  fgw_shape0x_series_ov_bf  SO's build (expression series) with shape_form overlay_bayes (T2embryo_s2xs_ovbf)
  BF2_O  fgw_shape0x_ov_bf2        overlay_bayes_brackets: the prior is the task's two brackets only (rule 7: without
                                   the five extra specimens)
  BFns_O fgw_shape0x_ov_bfns       overlay_bayes_noflip: no flip copies, each specimen turned onto the union's frame by
                                   the best of the four proper flips (rule 7: without the flip step)
Rescored, not counted: O fgw_shape0x_ov, SO fgw_shape0x_series_ov, and s2xs fgw_shape0x_series (main only, for its
guard against E47 and SO's construction checks). Every file is built once per task-seed under build main; O, SO,
BF_O and BF_SO are scored again under build prefix from the same in-memory files (clause (d)).

Steps (each refuses unless TRIES.md holds exactly one well-formed row E59 (bf_row), the tree is clean, HEAD is the run
commit and the committed reference records are the ones named below):
  g0 --seed S --commit SHA --out DIR   G0 at one build seed (target-free, not a try): the test pair's inputs (E7.5,
                                       E8.0, t = 1/2) with the validation board's panel and count rule, the test's
                                       prior
  g0read --commit SHA --out DIR        G0's verdict over the ten seeds into DIR/g0.json; exit 0 PASS, 4 STOP, 1 an error
  run --task NAME --commit SHA         the nine tasks, one process each, seed by seed (resumable); refused unless G0
                                       passed at this commit
  read --commit SHA                    the reading into experiments/coords/results/bf_gate/ (reading.json, reading.txt),
                                       never overwriting; the outcome line is the text's last line
  j7 [--out FILE]                      the J7 dry run of the O guard on committed records only (not a try, no run):
                                       E58's committed O-minus-floor per task against E55's; exit 0 pass, 3 fail
The script prints only "<task> seed <s> done", "wrote <path>", "g0 <seed> done", "g0 pass|stop" and "j7 pass|fail". It
reads no test board, test file or test anchor, and prints no count.
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
import traceback
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO / "experiments" / "harness"))
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(REPO / "src"))  # this checkout's code, not the checkout a shared venv's vec.pth names
import bf_row as br  # noqa: E402
import e47_floor_exact as e47  # noqa: E402
import nine_tasks as nt  # noqa: E402

from vec import harness  # noqa: E402
from vec import shape_interp as si  # noqa: E402
from vec import t2_build as tb  # noqa: E402
from vec.adoption import SEEDS  # noqa: E402
from vec.provenance import repo_state  # noqa: E402

TRIES = REPO / "experiments" / "harness" / "TRIES.md"
TRY_ID = br.ROW_ID
O = "fgw_shape0x_ov"  # noqa: E741  (the specification names the base O)
SO = "fgw_shape0x_series_ov"
S2XS = "fgw_shape0x_series"
BF_O, BF_SO = "fgw_shape0x_ov_bf", "fgw_shape0x_series_ov_bf"
BF2_O, BFNS_O = "fgw_shape0x_ov_bf2", "fgw_shape0x_ov_bfns"
ARMS = (BF_O, BF_SO, BF2_O, BFNS_O)  # the four counted tries
FILE_ARMS = (BF_O, BF_SO)  # the arms whose files row R71 orders
CONTROLS_E = (BF2_O, BFNS_O)  # clause (e)'s controls
BASE_OF = {BF_O: O, BF_SO: SO, BF2_O: O, BFNS_O: O}
FILE_OF = {BF_O: br.FILE_O, BF_SO: br.FILE_SO}
LABEL = {BF_O: "BF_O", BF_SO: "BF_SO", BF2_O: "BF2_O", BFNS_O: "BFns_O", O: "O", SO: "SO", S2XS: "s2xs"}
FORM_OF = {BF_O: "overlay_bayes", BF_SO: "overlay_bayes", BF2_O: "overlay_bayes_brackets",
           BFNS_O: "overlay_bayes_noflip"}
EXPRESSION_OF = {BF_O: "mixture", BF_SO: "series", BF2_O: "mixture", BFNS_O: "mixture"}
MAIN_NAMES = (O, SO, S2XS, *ARMS)  # with the floor, in the main pass
PREFIX_NAMES = (O, SO, BF_O, BF_SO)  # with the floor, scored again under prefix from the same files
BUILDS = ("main", "prefix")
TASKS = ("half", "late", "new1", "new2", "e75", "e725", "x1", "x2", "x3")
TEST_LIKE = ("half", "late", "new1", "new2")
E80_TASKS = ("half", "late")  # the truth collected like E7.75 (clause (f)), chosen by its protocol, by no score
E80_TRUTH = "T2/E8.0.h5ad"
EXPR_METRICS = ("de_score", "de_direction", "mmd_u", "variogram")
POS_METRICS = ("d2_shape", "occupancy_dice", "scale_log_ratio")
SHAPE_TERMS = ("d2_shape", "occupancy_dice")
NBHD = "neighborhood_mmd"
OTHER_GROUPS = ("expression_change", "cell_state", "local_organisation")
BAR_A, AHEAD, BAR_B, TRUTHS_B = 0.5, 8, 0.06, 4
NB_MERGED, NB_TRUTH, PREFIX_TOL, IDLE = -0.01, -0.02, 0.3, 0.3
EXPR_TOL, POS_TOL, RMS_TOL, GUARD_TOL = 1e-9, 1e-6, 1e-6, 0.01
G0_SEEDS = tuple(harness.BUILD_SEED_OFFSET + 1 + i for i in range(10))
G0_SHARE, G0_IDLE_DICE = 0.7, 0.9
G0_BOARD = "T2:embryo:val_interp"
G0_PAIR = {"before": "T2/E7.5.h5ad", "after": "T2/E8.0.h5ad"}
G0_T = 0.5
E47_DIR = REPO / "experiments" / "harness" / "results" / "e47_floor_exact" / "tasks"
E47_RUN_COMMIT = "9967cf1"
E47_SHA256 = {  # the committed E47 task files (sha256 of each at this script's writing)
    "T2e_E6.75_E7.5_E7.25.json": "e20f487c36c2fc706cdab027ab8a6ce72f082d6cf454fbeb63d932c9232a1010",
    "T2e_E7.25_E8.0_E7.5.json": "d926df49113c0d070fbf29cf09a488a908386bf319960bc270e1b47d4885168d",
    "T2e_E7.5_E8.25_late_E8.0.json": "59c345f2a64798c1362fe86b545efafed6ab6d8fa4d55a7254967fa7b3b0e041",
    "T2e_E7.5_E8.5_E8.0.json": "e6ab844a3f8e8a2fc1d296a25ba75768f3c42c90a34ce9ac9d6c633442d36bc1",
    "T2e_E7.5_E8.5_E8.25_late.json": "02c37ee2b5f987564429dde3ea98cde22a05b2ffbe240cbd15239399fac5013d",
    "T2e_E7.5_E8.75_E8.5.json": "878c5b5380d4e4046c94a711366bd258c0aea93b5eb06d868d3151eeb307a068",
    "T2e_E8.0_E8.75_E8.25_late.json": "03b0483396cd59cc6d4a013d9294a70ef06a9ec4ee27ba715024a5c8321aafb7",
    "T2e_E8.0_E8.75_E8.5.json": "c31044dcad6e19dd864c0032db415a7103dd27dd1f26901c2942ea92d9a154e5",
    "T2e_E8.25_late_E8.75_E8.5.json": "e40837de771777920abcb742ec22122f8345fd0e9de9d22c91b650544595d5d4"}
COMMITTED = {  # name -> (path, sha256, run commit), as rows E55, E51 and E58 record them
    "e55": (REPO / "experiments" / "coords" / "results" / "e55_overlay_stretch" / "reading.json",
            "e500c9fa4161720191562598af76f8fb28bc12b0c3b0dfa9fdd3e38a38153ab5",
            "4535f98c98ae1a2a8488185e8a99e4a2e8d92d86"),
    "e51": (REPO / "experiments" / "coords" / "results" / "e51_so_ohi" / "reading.json",
            "4c093dbb5f75dba05ba314eda5a160c292db249d3550e3edd3c9f1064ea5b894",
            "d0f0d6e54c5af45b3cd4038fc41354a73358e46c"),
    "e58": (REPO / "experiments" / "coords" / "results" / "e58_mix_lean" / "reading.json",
            "bd5da2ac3ab48e5a76c24b8494f29c975d4afe196400645b4c420e583ba322f7",
            "0a6d4ce8afd43bef2af5e43820ccd22e93cf8596"),
}
SCRIPT = "experiments/coords/bf_read.py"
RUNS = REPO / "runs" / "bf_gate"
RESULTS = REPO / "experiments" / "coords" / "results" / "bf_gate"
REFUSAL = "cloud guard flags"

CLAUSES = {
    "a": "(a) the five-truth merged gated board gain against its base is at least +0.5, with D_s above 0 in at least 8 "
         "of 10 seeds.",
    "b": "(b) the ungated d2_shape plus occupancy_dice skill (their sum) gains at least +0.06 merged, and above 0 on "
         "at least 4 of the 5 truths (the genericity claim).",
    "c": "(c) neighborhood_mmd's skill changes by at least -0.01 merged, and by no less than -0.02 on any truth.",
    "d": "(d) under build prefix the five-truth merged board gain is within 0.3 of (a)'s (the coin does not carry it).",
    "e": "(e) BF_O's ungated d2_shape plus occupancy_dice is above each control's (BF2_O, BFns_O) on the five-truth "
         "merge, with at least 8 of 10 D_s ahead (rule 7: the five extra specimens and the flip step each earn their "
         "place); a failed (e) is STOP for every arm.",
    "f": "(f) on the E8.0 truth (half and late merged; chosen by its collection protocol, by no score) the ungated "
         "d2_shape plus occupancy_dice and the gated board are both above 0.",
}
PASS_RULE = (
    "Per file arm against its own base (BF_O against O, BF_SO against SO), under build main unless named, on the J4 "
    "merge of the nine tasks into five truths at equal weight (E8.0 = half, late; E8.25_late = new1, x1; E8.5 = new2, "
    "x2, x3; E7.5 = e75; E7.25 = e725), D_s the mean over truths of each truth's task-mean difference at seed s. An "
    "arm passes only if all of " + " ".join(CLAUSES.values()) + " hold and its merged board gain is not idle (below "
    "0.3 in size: STOP). Its file is built only if it passes; neither file if (e) fails.")
VOID_CONDITIONS = (
    "a missing task, seed or file (the floor and every name of the main pass, the floor and O, SO, BF_O, BF_SO under "
    "prefix)",
    "a result from another run commit or an unclean tree, not the bake protocol, or a test board",
    "G0's record missing or not PASS at the run commit",
    "any refusal (an arm, a control or a base)",
    "an arm's count or expression matrix off its base's, or an expression term (de_score, de_direction, mmd_u, "
    "variogram) off its base's by more than 1e-9 in skill",
    "an arm's size target off its base's, or its RMS radius off that target by more than 1e-6 relative where it did "
    "not jitter",
    "an arm's recorded form or prior off its specification (the prior holds every stage of BAYES_PRIOR_STAGES that "
    "does not lie strictly between the task's brackets, and so never the truth)",
    "SO's guard, row E51's own void checks: SO's de_score, de_direction, mmd_u or variogram skill off s2xs's by more "
    "than 1e-9, or its d2_shape, occupancy_dice or scale_log_ratio skill off O's by more than 1e-6, in any task, seed "
    "or build",
    "the guard (J7): O's ten-seed mean minus the floor's off E55's committed per-task figure by more than 0.01; s2xs "
    "off E47's committed per-seed scores by more than 0.01 (per seed or ten-seed mean). The floor is not checked: the "
    "harness holds it at 50 on every seed or raises",
)
PRICE_LINES = (
    "J1: the run applies each file's own build path end to end except the size step (the harness's task size rule; the "
    "file takes the pre-registered size, which the form never moves); the prior holds, on every task as at the test "
    "pair, the released stages that do not lie strictly between the brackets (so never the truth, J2).",
    "J4: five truths at equal weight; clause (f) reads the one truth collected like E7.75.",
    "J5: a diagnostic; the noise floor 6.54 exceeds twice any predicted gain. J6: no time quantity sets a step "
    "(the union's mix follows each task's own t).",
)


# ------------------------------------------------------------------------------------------------ checks
def _git(*args: str) -> subprocess.CompletedProcess:
    return subprocess.run(["git", "-C", str(REPO), *args], capture_output=True, text=True,
                          env={**os.environ, "GIT_OPTIONAL_LOCKS": "0"})


def _sha_file(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def reference_problems() -> list[str]:
    """The committed reference records the guard reads: present and unchanged."""
    why = []
    for name, sha in E47_SHA256.items():
        p = E47_DIR / name
        if not p.exists() or _sha_file(p) != sha:
            why.append(f"E47's committed {name} is missing or changed")
    for name, (path, sha, _) in COMMITTED.items():
        if not path.exists() or _sha_file(path) != sha:
            why.append(f"{name.upper()}'s committed reading is missing or changed")
    return why


def commit_problems(commit: str) -> list[str]:
    """Why this checkout may not run or read row E59 at `commit` (empty: it may)."""
    why = []
    head = _git("rev-parse", "HEAD").stdout.strip()
    if head != commit:
        why.append(f"HEAD is {head[:12]}, not the run commit {commit[:12]}")
    return why + reference_problems()


def _own_code() -> None:
    src = (REPO / "src").resolve()
    for mod in (harness, tb, si):
        if src not in Path(mod.__file__).resolve().parents:
            sys.exit(f"{mod.__name__} is loaded from {mod.__file__}, not from this checkout's src/")


def _refuse_unless_ready(commit: str) -> None:
    _own_code()
    why = br.row_problems(TRIES)
    if why:
        sys.exit("E59 has no well-formed row in experiments/harness/TRIES.md (a row is committed before its run): "
                 + "; ".join(why))
    if repo_state(REPO).unclean != ():
        sys.exit("the tree is not clean: row E59 runs on a committed tree")
    why = commit_problems(commit)
    if why:
        sys.exit("refused: " + "; ".join(why))


def _task(name: str):
    tasks = nt.all_tasks()
    if name not in TASKS or name not in tasks:
        sys.exit(f"unknown task {name}; row E59 runs {list(TASKS)}")
    t = tasks[name]
    if getattr(t, "target_filter", None) is not None:
        sys.exit(f"{name}: a task with a target filter; E59 does not run on one")
    return t


def _is_test_board(board) -> bool:
    return str(board.split).startswith("test") or "test" in str(board.key).split(":")[-1]


# ------------------------------------------------------------------------------------------------ the prior
def stage_day(stage: str) -> float:
    return float(harness.STAGE_DAYS[f"T2/{stage}.h5ad"])


def prior_stages(lo: float, hi: float, truth: str | None = None) -> tuple[str, ...]:
    """The stages of BAYES_PRIOR_STAGES that do not lie strictly between the brackets' days lo and hi (the review's fix
    B1): at the test pair (E7.5, E8.0) all seven, since no released stage lies between them; on a practice task the
    truth, which lies between its brackets, leaves with every other stage there, so no practice prior holds a stage the
    test's cannot. `truth` (a target file such as T2/E8.0.h5ad) is checked to be out (J2)."""
    lo, hi = float(lo), float(hi)
    out = tuple(s for s in si.BAYES_PRIOR_STAGES if not lo < stage_day(s) < hi)
    if truth is not None and harness.stem(truth) in out:
        raise harness.LeakageError(f"{harness.stem(truth)}: the truth entered its prior")
    return out


def task_prior_stages(t) -> tuple[str, ...]:
    """A practice task's prior stages: prior_stages at its brackets' days, its truth checked out."""
    return prior_stages(harness.STAGE_DAYS[t.inputs["before"]], harness.STAGE_DAYS[t.inputs["after"]], t.target)


def pair_prior_stages() -> tuple[str, ...]:
    """The test pair's prior stages (G0 and the builder): all seven."""
    return prior_stages(harness.STAGE_DAYS[G0_PAIR["before"]], harness.STAGE_DAYS[G0_PAIR["after"]])


def load_prior(stages, data: Path) -> si.BayesPrior:
    """The prior's released specimens, each loaded as a prediction input is (synthetic records dropped), coordinates
    only."""
    from vec.io import COORD_KEY
    from vec.prep import load_stage

    coords = {}
    for s in stages:
        a = load_stage(Path(data) / "T2" / f"{s}.h5ad", lineage=False, drop_duplicates=True)
        coords[s] = np.asarray(a.obsm[COORD_KEY], dtype=float)[:, :3]
    return si.bayes_prior(coords)


# ------------------------------------------------------------------------------------------------ G0 (target-free)
def g0_inputs(data: Path) -> dict:
    """The test pair's brackets (coordinates, panel expression), t and count, with the validation board's panel and
    count rule (as experiments/coords/osr_gate.py's test pair). Never a target, a test board or a hidden count."""
    from dataclasses import replace

    from vec.boards import load_boards
    from vec.io import COORD_KEY, dense
    from vec.prep import load_stage

    board = replace(load_boards()[G0_BOARD], anchors={})
    if _is_test_board(board) or board.counts_hidden:
        raise PermissionError("G0 reads a validation board's panel and count rule only")
    genes = list(board.genes)
    st = {r: load_stage(Path(data) / f, genes=genes, drop_duplicates=True) for r, f in G0_PAIR.items()}
    Cb, Ca = (np.asarray(st[r].obsm[COORD_KEY], float)[:, :3] for r in ("before", "after"))
    Eb, Ea = (dense(st[r][:, genes].X).astype(np.float32) for r in ("before", "after"))
    return {"Cb": Cb, "Ca": Ca, "Eb": Eb, "Ea": Ea, "t": G0_T, "n": harness.count_for(board, st["before"].n_obs)}


def idle_dice(U, q) -> float:
    """The challenge's idle case: the 16-cube Dice between U's occupied voxels (U canonicalised) and the same number of
    q's highest-mass voxels (q summed from the 48-cube to the 16-cube; ties by voxel id)."""
    occ = si.occupancy(si.canonical_frame(U)[0])
    f = si.BAYES_GRID_N // si.OCC_GRID_N
    G = si.OCC_GRID_N
    q16 = np.asarray(q, float).reshape(G, f, G, f, G, f).sum(axis=(1, 3, 5)).ravel()
    top = np.sort(np.argsort(-q16, kind="stable")[:len(occ)])
    return si.dice(occ, top)


def z_extent(C) -> float:
    z = np.asarray(C, float)[:, 2]
    return float(np.percentile(z, 99) - np.percentile(z, 1))


def g0_seed(p: dict, seed: int, prior: si.BayesPrior) -> dict:
    """One build seed of G0: O's union and BF's file at the test pair (size rule at RMS 1, the shape terms ignore
    scale); the neighbour share, BF's cloud guard, the idle Dice; reported only: BF's z extent over O's (1st to 99th
    percentile of z) and the angle between the union's thinnest axis and z."""
    kw = {"target_rms": 1.0, "t_pos": 0.0, "positions": tb.POSITIONS[1]}
    ro = tb.interpolate_clouds(p["Cb"], p["Ca"], p["Eb"], p["Ea"], p["t"], p["n"], seed, shape_form="overlay", **kw)
    rb = tb.interpolate_clouds(p["Cb"], p["Ca"], p["Eb"], p["Ea"], p["t"], p["n"], seed, shape_form="overlay_bayes",
                               shape_prior=prior, **kw)
    guard = tb.cloud_guard(rb.coords, rb.pool_coords, seed=seed, source_weights=rb.pool_weights)
    q, _ = si.consensus_density(prior)
    Vt = si.canonical_frame(ro.coords)[2]
    return {"share": si.neighbour_share(ro.coords, rb.coords), "guard_ok": bool(guard.ok),
            "guard_flags": list(guard.flags), "idle_dice": idle_dice(ro.coords, q),
            "z_extent_ratio": z_extent(rb.coords) / z_extent(ro.coords),
            "angle_thin_axis_z_deg": float(np.degrees(np.arccos(min(1.0, abs(float(Vt[2][2])))))),
            "prior": "|".join(prior.names)}


def _g0_path(out: Path, seed: int) -> Path:
    return Path(out) / f"g0_{seed}.json"


def cmd_g0(a) -> int:
    _refuse_unless_ready(a.commit)
    seed = int(a.seed)
    if seed not in G0_SEEDS:
        sys.exit(f"G0 runs the build seeds {G0_SEEDS[0]} to {G0_SEEDS[-1]} only")
    path = _g0_path(Path(a.out), seed)
    if path.exists():
        sys.exit(f"{path} exists")
    rec = {"seed": seed, "commit": a.commit, "script": SCRIPT}
    try:
        p = g0_inputs(Path(a.data))
        rec.update(g0_seed(p, seed, load_prior(pair_prior_stages(), Path(a.data))))
    except Exception as e:  # recorded: a refusal or an error at G0 is a STOP
        rec["error"] = f"{type(e).__name__}: {e}"
        rec["traceback"] = traceback.format_exc(limit=5)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(rec, indent=1, default=float))
    print(f"g0 {seed} done", flush=True)
    return 0


def g0_verdict(recs: dict[int, dict], commit: str) -> dict:
    """STOP if a seed is missing, from another commit or in error, if any seed's median neighbour share is under 0.7,
    if BF's cloud guard refuses on any seed, or if the median over seeds of the idle Dice is above 0.9 (the idle case).
    The z extent ratio and the angle are reported and decide nothing."""
    why = []
    for s in G0_SEEDS:
        r = recs.get(s)
        if r is None:
            why.append(f"seed {s}: no record")
        elif r.get("commit") != commit:
            why.append(f"seed {s}: another commit")
        elif "error" in r:
            why.append(f"seed {s}: {r['error']}")
    if why:
        return {"verdict": "STOP", "why": why, "commit": commit}
    shares = [recs[s]["share"] for s in G0_SEEDS]
    dices = [recs[s]["idle_dice"] for s in G0_SEEDS]
    stops = []
    if min(shares) < G0_SHARE:
        stops.append(f"the median neighbour share is under {G0_SHARE} on {sum(x < G0_SHARE for x in shares)} seeds")
    if not all(recs[s]["guard_ok"] for s in G0_SEEDS):
        stops.append("BF's cloud guard refuses on a seed")
    idle = bool(np.median(dices) > G0_IDLE_DICE)
    if idle:
        stops.append(f"the idle case: the union's median Dice against q's top voxels is above {G0_IDLE_DICE}")
    return {"verdict": "STOP" if stops else "PASS", "why": stops, "commit": commit, "idle": idle,
            "per_seed": {str(s): {k: recs[s][k] for k in ("share", "guard_ok", "guard_flags", "idle_dice",
                                                         "z_extent_ratio", "angle_thin_axis_z_deg")}
                         for s in G0_SEEDS},
            "median_share": float(np.median(shares)), "median_idle_dice": float(np.median(dices)),
            "reported_deciding_nothing": {
                "median_z_extent_ratio": float(np.median([recs[s]["z_extent_ratio"] for s in G0_SEEDS])),
                "median_angle_thin_axis_z_deg": float(np.median([recs[s]["angle_thin_axis_z_deg"] for s in G0_SEEDS]))},
            "prior": recs[G0_SEEDS[0]]["prior"]}


def load_g0(out: Path) -> dict[int, dict]:
    recs = {}
    for s in G0_SEEDS:
        p = _g0_path(out, s)
        if p.exists():
            recs[s] = json.loads(p.read_text())
    return recs


def cmd_g0read(a) -> int:
    _refuse_unless_ready(a.commit)
    out = Path(a.out)
    path = out / "g0.json"
    if path.exists():
        sys.exit(f"{path} exists")
    res = g0_verdict(load_g0(out), a.commit)
    path.write_text(json.dumps(res, indent=1, default=float))
    print(f"g0 {'pass' if res['verdict'] == 'PASS' else 'stop'}", flush=True)
    return 0 if res["verdict"] == "PASS" else 4


def g0_problems(g0_dir: Path, commit: str) -> list[str]:
    p = Path(g0_dir) / "g0.json"
    if not p.exists():
        return ["no G0 record"]
    g = json.loads(p.read_text())
    if g.get("verdict") != "PASS" or g.get("commit") != commit:
        return ["G0 did not pass at the run commit"]
    return []


# ------------------------------------------------------------------------------------------------ the arms
def arm_candidate(name: str, holder: dict) -> harness.Candidate:
    """One counted arm: O's (or SO's) build with the arm's form; `holder["prior"]` is the task's prior (not used by
    BF2_O, whose prior is the task's two brackets)."""
    form, expression = FORM_OF[name], EXPRESSION_OF[name]

    def build(ctx):
        if form == "overlay_bayes_brackets":
            return harness._fgw(ctx, t_pos=0.0, expression=expression, positions=tb.POSITIONS[1], shape_form=form)
        inner = tb.interpolate_clouds  # harness._fgw stays as it is: the prior goes in at its one interpolation

        def with_prior(*a, **kw):
            return inner(*a, **kw, shape_prior=holder["prior"])

        tb.interpolate_clouds = with_prior
        try:
            return harness._fgw(ctx, t_pos=0.0, expression=expression, positions=tb.POSITIONS[1], shape_form=form)
        finally:
            tb.interpolate_clouds = inner

    return harness.Candidate(
        name, frozenset({"T2"}), (), build,
        {"what": f"{LABEL[BASE_OF[name]]} with the shape form {form} (row E59, the gate of row R71's files)",
         "positions": tb.POSITIONS[1], "stream": tb.FLOOR_EXACT_STREAM, "shape_form": form, "expression": expression,
         "scope": "whole"}, requires=harness._interp)


@contextlib.contextmanager
def captured():
    """While active, record what harness.score_seeds is handed (the board) and each interpolate_clouds record under the
    file being built (flat values only, no count)."""
    seen: dict = {"stage": {}, "refused": {}, "records": {}, "current": None}
    orig_ss, orig_ic = harness.score_seeds, tb.interpolate_clouds

    def score_seeds(t, board, target, ref, inputs, candidates, *a, **kw):
        seen.update(board=board)
        return orig_ss(t, board, target, ref, inputs, candidates, *a, **kw)

    def interpolate_clouds(*a, **kw):
        res = orig_ic(*a, **kw)
        r = res.record
        seen["records"][seen["current"]] = {
            "target_rms": float(r["target_rms"]), "n_jittered": int(r.get("n_jittered", 0)), "t": float(r["t"]),
            "t_mix": float(r["t_mix"]),
            **{k: v for k, v in r.items() if k.startswith("shape_interp") and not isinstance(v, (dict, list))}}
        return res

    harness.score_seeds, tb.interpolate_clouds = score_seeds, interpolate_clouds
    try:
        yield seen
    finally:
        harness.score_seeds, tb.interpolate_clouds = orig_ss, orig_ic


def observed(base, seen):
    """A base unchanged; its file is kept for the checks, the prefix pass and as a refused arm's stand-in."""
    def build(ctx):
        if _is_test_board(seen["board"]):
            raise RuntimeError("E59 met a test board")
        seen["current"] = base.name
        pred = base.build(ctx)
        seen["stage"][(base.name, ctx.seed - harness.BUILD_SEED_OFFSET)] = pred
        return pred

    return dataclasses.replace(base, build=build)


def guarded(base, seen):
    """An arm: its file, or on a refusal (ShapeRefused, the cloud guard) its base's file of the seed (a refusal makes
    the read VOID; the run goes on)."""
    def build(ctx):
        s = ctx.seed - harness.BUILD_SEED_OFFSET
        if _is_test_board(seen["board"]):
            raise RuntimeError("E59 met a test board")
        seen["current"] = base.name
        try:
            pred, why = base.build(ctx), None
        except si.ShapeRefused as e:
            pred, why = None, f"shape: {e}"
        except ValueError as e:
            if REFUSAL not in str(e):
                raise
            pred, why = None, "cloud guard"
        if pred is None:
            pred = seen["stage"][(BASE_OF[base.name], s)]
        seen["refused"][(base.name, s)] = why
        seen["stage"][(base.name, s)] = pred
        return pred

    return dataclasses.replace(base, build=build)


def replayed(base, stages: dict):
    """The prefix pass: the main pass's very file of the seed, scored again (no rebuild)."""
    def build(ctx):
        return stages[(base.name, ctx.seed - harness.BUILD_SEED_OFFSET)]

    return dataclasses.replace(base, build=build)


def main_candidates(seen, holder: dict) -> list:
    bases = [observed(harness.CANDIDATES[n], seen) for n in (O, SO, S2XS)]
    return [harness.CANDIDATES[harness.FLOOR], *bases, *(guarded(arm_candidate(n, holder), seen) for n in ARMS)]


def prefix_candidates(stages: dict, holder: dict) -> list:
    by_name = {**{n: harness.CANDIDATES[n] for n in (O, SO)}, **{n: arm_candidate(n, holder) for n in FILE_ARMS}}
    return [harness.CANDIDATES[harness.FLOOR], *(replayed(by_name[n], stages) for n in PREFIX_NAMES)]


def x_digest(stage) -> str:
    X = np.ascontiguousarray(stage.X, dtype=np.float32)
    return hashlib.sha256(np.asarray(X.shape, np.int64).tobytes() + X.tobytes()).hexdigest()


def rms_off(stage, target: float) -> float:
    C = np.asarray(stage.coords, dtype=np.float64)
    return abs(float(np.sqrt(((C - C.mean(0)) ** 2).sum(1).mean())) / float(target) - 1.0)


def one_seed(t, s: int, data: Path, prior: si.BayesPrior) -> dict:
    """One task-seed: every file built and scored under main, O, SO, BF_O and BF_SO scored again under prefix."""
    holder = {"prior": prior}
    out: dict = {"builds": {}, "facts": {}}
    with captured() as seen:
        out["builds"]["main"] = harness.run_task(t, main_candidates(seen, holder), seeds=[s], data_dir=data,
                                                 build="main", log=None)
    if _is_test_board(seen["board"]):
        raise RuntimeError("E59 met a test board")
    stages = dict(seen["stage"])
    with captured() as seen2:
        out["builds"]["prefix"] = harness.run_task(t, prefix_candidates(stages, holder), seeds=[s], data_dir=data,
                                                   build="prefix", log=None)
    if _is_test_board(seen2["board"]):
        raise RuntimeError("E59 met a test board")
    recs = seen["records"]
    out["facts"] = {
        "board_split": str(seen["board"].split), "prior": "|".join(prior.names),
        "n_equal": {n: int(len(stages[(n, s)].X)) == int(len(stages[(BASE_OF[n], s)].X)) for n in ARMS},
        "x": {n: x_digest(stages[(n, s)]) for n in MAIN_NAMES},
        "refused": {n: seen["refused"].get((n, s)) for n in ARMS},
        "records": {n: recs.get(n, {}) for n in MAIN_NAMES},
        "rms_off": {n: (rms_off(stages[(n, s)], recs[n]["target_rms"])
                        if n in recs and not recs[n]["n_jittered"] and not seen["refused"].get((n, s)) else None)
                    for n in ARMS}}
    return out


def task_prior(t, data: Path) -> si.BayesPrior:
    return load_prior(task_prior_stages(t), data)


def cmd_run(a) -> None:
    _refuse_unless_ready(a.commit)
    t = _task(a.task)
    runs = Path(a.runs)
    why = g0_problems(runs / "g0", a.commit)
    if why:
        sys.exit("refused (G0): " + "; ".join(why))
    holder = {"prior": None}
    for c in (*(harness.CANDIDATES[n] for n in (O, SO, S2XS)), *(arm_candidate(n, holder) for n in ARMS)):
        if not c.applies(t):
            sys.exit(f"{a.task}: {c.name} does not apply")
    prior = task_prior(t, Path(a.data))
    folder = runs / "run" / a.task
    folder.mkdir(parents=True, exist_ok=True)
    for s in SEEDS:
        path = folder / f"seed{s}.json"
        if path.exists():
            continue
        res = one_seed(t, s, Path(a.data), prior)
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
    """Per seed, one quantity of one file: 'board', 'group:<g>' (gated), 'skill:<m>' or 'shape_sum' (the ungated
    d2_shape plus occupancy_dice skill)."""
    rows = _rows(recs, build)
    if what == "board":
        return np.array([x["scores"][name] for x in rows], float)
    if what.startswith("group:"):
        return np.array([x["groups"][name][what[6:]] for x in rows], float)
    if what.startswith("skill:"):
        return np.array([x["skills"][name][what[6:]] for x in rows], float)
    if what == "shape_sum":
        return np.array([sum(x["skills"][name][m] for m in SHAPE_TERMS) for x in rows], float)
    raise ValueError(what)


def merged_D(recs, groups: dict[str, list[str]], a: str, b: str, what: str, build: str = "main") -> np.ndarray:
    """D_s: the mean over truths of each truth's task-mean difference (a minus b) at each seed."""
    return e47.merged_seeds({n: series(recs[n], a, what, build) - series(recs[n], b, what, build)
                             for names in groups.values() for n in names}, groups)


def per_truth(recs, groups: dict[str, list[str]], a: str, b: str, what: str, build: str = "main") -> dict:
    """Per truth, the mean over its tasks of the ten-seed mean difference (a minus b)."""
    return {truth: float(np.mean([np.mean(series(recs[n], a, what, build) - series(recs[n], b, what, build))
                                  for n in names])) for truth, names in groups.items()}


def _stat(D) -> dict:
    D = np.asarray(D, float)
    return {"mean": float(D.mean()), "se": float(D.std(ddof=1) / np.sqrt(len(D))), "ahead": int((D > 0).sum()),
            "behind": int((D < 0).sum())}


def merges(tasks: dict | None = None) -> dict[str, dict[str, list[str]]]:
    """The merges fixed in row E59 (J4): five truths (all nine tasks), the E8.0 truth (clause (f)); reported only: the
    test-like tasks on their three truths and the four truths without E7.25."""
    tasks = nt.all_tasks() if tasks is None else tasks
    return {"five": e47.truths({n: tasks[n] for n in TASKS}),
            "e80": e47.truths({n: tasks[n] for n in E80_TASKS}),
            "test_like": e47.truths({n: tasks[n] for n in TEST_LIKE}),
            "four": e47.truths({n: tasks[n] for n in TASKS if n != "e725"})}


def arm_figures(recs, m: dict, arm: str) -> dict:
    """The figures the pass rule reads for one file arm against its base."""
    base = BASE_OF[arm]
    five = m["five"]
    return {"board": merged_D(recs, five, arm, base, "board"),
            "board_prefix": merged_D(recs, five, arm, base, "board", build="prefix"),
            "shape_sum": merged_D(recs, five, arm, base, "shape_sum"),
            "shape_sum_truths": per_truth(recs, five, arm, base, "shape_sum"),
            "nbhd": merged_D(recs, five, arm, base, f"skill:{NBHD}"),
            "nbhd_truths": per_truth(recs, five, arm, base, f"skill:{NBHD}"),
            "e80_board": per_truth(recs, m["e80"], arm, base, "board")[E80_TRUTH],
            "e80_shape_sum": per_truth(recs, m["e80"], arm, base, "shape_sum")[E80_TRUTH]}


def decide_arm(f: dict, e_holds: bool) -> dict:
    """THE PASS RULE for one file arm (CLAUSES) on its figures `f` (arm_figures) and clause (e)'s verdict."""
    D, Dp, S, NB = (np.asarray(f[k], float) for k in ("board", "board_prefix", "shape_sum", "nbhd"))
    c = {"a": bool(D.mean() >= BAR_A and (D > 0).sum() >= AHEAD),
         "b": bool(S.mean() >= BAR_B and sum(v > 0 for v in f["shape_sum_truths"].values()) >= TRUTHS_B),
         "c": bool(NB.mean() >= NB_MERGED and min(f["nbhd_truths"].values()) >= NB_TRUTH),
         "d": bool(abs(float(Dp.mean()) - float(D.mean())) <= PREFIX_TOL),
         "e": bool(e_holds),
         "f": bool(f["e80_board"] > 0 and f["e80_shape_sum"] > 0)}
    idle = bool(abs(float(D.mean())) < IDLE)
    return {"clauses": c, "idle": idle, "passes": all(c.values()) and not idle}


def clause_e(recs, m: dict) -> dict:
    """(e): BF_O's ungated d2_shape plus occupancy_dice against each control on the five-truth merge."""
    out = {}
    for ctl in CONTROLS_E:
        D = merged_D(recs, m["five"], BF_O, ctl, "shape_sum")
        out[ctl] = {**_stat(D), "holds": bool(D.mean() > 0 and (D > 0).sum() >= AHEAD)}
    return {"controls": out, "holds": all(v["holds"] for v in out.values())}


def outcome_of(passes: dict[str, bool], idles: dict[str, bool]) -> str:
    if passes[BF_O] and passes[BF_SO]:
        return "PASS_BOTH"
    if passes[BF_O]:
        return "PASS_O"
    if passes[BF_SO]:
        return "PASS_SO"
    return "IDLE" if all(idles.values()) else "STOP"


def _e47_by_key(e47_dir: Path | None = None) -> dict[str, dict]:
    out = {}
    for p in sorted(Path(e47_dir or E47_DIR).glob("T2e_*.json")):
        r = json.loads(p.read_text())
        out[r["task"]["key"]] = r
    return out


def committed_figures() -> tuple[dict, list[str]]:
    """The guard's committed references: E47's per-seed scores by task key, E55's O-minus-floor and E51's SO-minus-s2xs
    per task (the committed readings of the sha256 named in COMMITTED)."""
    why = reference_problems()
    e47_runs = _e47_by_key()
    ref = {"e47": {}, "o": {}, "so": {}}
    for key, r in e47_runs.items():
        if not r["provenance"]["commit"].startswith(E47_RUN_COMMIT) or not r["provenance"]["tree_clean"] \
                or r["build"] != "main" or r["protocol"] != "bake":
            why.append(f"{key}: E47's file is not its clean main/bake run at {E47_RUN_COMMIT}")
        ref["e47"][key] = e47.per_seed(r)[0]
    r55 = json.loads(COMMITTED["e55"][0].read_text()) if COMMITTED["e55"][0].exists() else {}
    r51 = json.loads(COMMITTED["e51"][0].read_text()) if COMMITTED["e51"][0].exists() else {}
    if r55.get("void") is not False or r55.get("commit") != COMMITTED["e55"][2]:
        why.append("E55's reading is void or not from its run")
    if r51.get("void") is not False or r51.get("commit") != COMMITTED["e51"][2]:
        why.append("E51's reading is void or not from its run")
    for name in TASKS:
        v = (((r55.get("tasks") or {}).get(name) or {}).get("board_minus_floor") or {}).get(O)
        w = ((((((r51.get("arms") or {}).get(SO) or {}).get("reported_deciding_nothing") or {}).get("per_task") or {})
              .get(name) or {}).get("board_minus_control"))
        if v is None or w is None:
            why.append(f"{name}: a committed guard figure is missing")
        else:
            ref["o"][name], ref["so"][name] = float(v), float(w)
    return ref, why


def guard_problems(recs, ref: dict, tasks: dict | None = None) -> tuple[list[str], dict, dict]:
    """The guard (J7) on the run's main pass, per task, each within GUARD_TOL board points: s2xs against E47's committed
    per-seed scores (per seed and ten-seed mean) and O's ten-seed mean minus the floor's against E55's committed
    per-task figure. The floor is not checked: harness.score_seeds holds it at 50 on every seed or raises, so a check
    could not fail. SO's guard is row E51's construction checks (check_run); SO's ten-seed mean minus s2xs's against
    E51's committed per-task figure is returned third, reported and deciding nothing (the review's closure B)."""
    tasks = nt.all_tasks() if tasks is None else tasks
    why, worst, so_report = [], {}, {}
    for name in TASKS:
        key = tasks[name].key
        old = ref["e47"].get(key)
        fl, o, so, s2 = (series(recs[name], n, "board") for n in (harness.FLOOR, O, SO, S2XS))
        devs = {}
        if old is None or S2XS not in old:
            why.append(f"{name}: no committed E47 scores")
        else:
            ref_n = np.asarray(old[S2XS], float)
            if ref_n.shape != s2.shape:
                why.append(f"{name}: E47's {S2XS} is of another length")
            else:
                devs[f"{S2XS}_seed"] = float(np.abs(s2 - ref_n).max())
                devs[f"{S2XS}_mean"] = abs(float(s2.mean() - ref_n.mean()))
        if name in ref["o"]:
            devs["o_minus_floor"] = abs(float(np.mean(o - fl)) - ref["o"][name])
        else:
            why.append(f"{name}: no committed E55 figure for O")
        for k, v in devs.items():
            if not v <= GUARD_TOL:
                why.append(f"{name}: the guard's {k} is {v:.4f} off its committed figure (tolerance {GUARD_TOL})")
        worst[name] = devs
        if name in ref.get("so", {}):
            so_report[name] = abs(float(np.mean(so - s2)) - ref["so"][name])
    return why, worst, so_report


def check_run(recs: dict[str, list[dict]], why_in: list[str], commit: str) -> list[str]:
    """The void conditions on the per-seed records (task -> its records in seed order), the guard aside."""
    why = list(why_in)
    tasks = nt.all_tasks()
    for task in TASKS:
        rs = recs.get(task, [])
        if [r.get("seed") for r in rs] != list(SEEDS):
            why.append(f"{task}: seeds are not {list(SEEDS)}")
            continue
        want_prior = "|".join(task_prior_stages(tasks[task]))
        for r in rs:
            tag = f"{task} seed {r['seed']}"
            if r.get("run_commit") != commit:
                why.append(f"{tag}: another run commit")
            for b, names in (("main", MAIN_NAMES), ("prefix", PREFIX_NAMES)):
                res = (r.get("builds") or {}).get(b)
                if res is None:
                    why.append(f"{tag}: build {b} missing")
                    continue
                p = res["provenance"]
                if p["commit"] != commit or not p["tree_clean"]:
                    why.append(f"{tag}: another commit or an unclean tree")
                if res["build"] != b or res["protocol"] != "bake" or [x["seed"] for x in res["rows"]] != [r["seed"]]:
                    why.append(f"{tag}: not the bake protocol under build {b} at this seed")
                if str(res["board"].get("split", "")).startswith("test"):
                    why.append(f"{tag}: a test board")
                row = res["rows"][0]
                if any(n not in row["scores"] or n not in row["skills"] or n not in row.get("groups", {})
                       for n in (harness.FLOOR, *names)):
                    why.append(f"{tag}: a file is missing under {b}")
            f = r.get("facts") or {}
            if str(f.get("board_split", "")).startswith("test"):
                why.append(f"{tag}: a test board")
            if f.get("prior") != want_prior:
                why.append(f"{tag}: the prior is not the task's ({want_prior}): the truth must stay out")
            for n in ARMS:
                base = BASE_OF[n]
                if (f.get("refused") or {}).get(n):
                    why.append(f"{tag}: {LABEL[n]} refused ({f['refused'][n]})")
                    continue
                if not (f.get("n_equal") or {}).get(n):
                    why.append(f"{tag}: {LABEL[n]}'s count differs from its base's")
                if (f.get("x") or {}).get(n) != (f.get("x") or {}).get(base):
                    why.append(f"{tag}: {LABEL[n]}'s expression matrix differs from its base's")
                rec, rec_b = (f.get("records") or {}).get(n) or {}, (f.get("records") or {}).get(base) or {}
                if rec.get("shape_interp") != FORM_OF[n] or rec_b.get("shape_interp") != "overlay":
                    why.append(f"{tag}: {LABEL[n]}'s or its base's recorded form is not its specification's")
                want = "before|after" if n == BF2_O else want_prior
                if rec.get("shape_interp_prior") != want:
                    why.append(f"{tag}: {LABEL[n]}'s recorded prior is not {want}")
                if rec.get("target_rms") is None or rec.get("target_rms") != rec_b.get("target_rms"):
                    why.append(f"{tag}: {LABEL[n]}'s size target differs from its base's")
                off = (f.get("rms_off") or {}).get(n)
                if off is not None and not off <= RMS_TOL:
                    why.append(f"{tag}: {LABEL[n]}'s RMS radius is {off:.2e} off its target")
                if "main" in (r.get("builds") or {}):
                    row = r["builds"]["main"]["rows"][0]
                    for m in EXPR_METRICS:
                        x, y = row["skills"].get(n, {}).get(m), row["skills"].get(base, {}).get(m)
                        if x is None or y is None or not abs(float(x) - float(y)) <= EXPR_TOL:
                            why.append(f"{tag}: {LABEL[n]}'s {m} is off its base's")
            why += so_construction_problems(r, tag)
    return why


def so_construction_problems(r: dict, tag: str) -> list[str]:
    """SO's guard, row E51's own void checks (which E51's committed run at d0f0d6e applied on all nine tasks and passed,
    not void): SO's de_score, de_direction, mmd_u and variogram skills within 1e-9 of s2xs's under build main, and its
    d2_shape, occupancy_dice and scale_log_ratio skills within 1e-6 of O's under each build."""
    why = []
    for b, pairs in (("main", ((EXPR_METRICS, S2XS, EXPR_TOL), (POS_METRICS, O, POS_TOL))),
                     ("prefix", ((POS_METRICS, O, POS_TOL),))):
        res = (r.get("builds") or {}).get(b)
        if res is None:
            continue
        skills = res["rows"][0]["skills"]
        for metrics, other, tol in pairs:
            for m in metrics:
                x, y = skills.get(SO, {}).get(m), skills.get(other, {}).get(m)
                if x is None or y is None or not abs(float(x) - float(y)) <= tol:
                    why.append(f"{tag}: SO's {m} is off {LABEL[other]}'s under {b} (row E51's construction check)")
    return why


def reading(recs: dict[str, list[dict]], why_in: list[str], commit: str, g0: dict | None,
            ref: tuple[dict, list[str]] | None = None, run_started: bool = False) -> dict:
    """The reading: G0's STOP, VOID, or the per-arm verdicts and the fixed outcome. The four tries count from the run
    step's start (rule 4): `run_started` (the launcher's run step began) or any seed record."""
    if g0 is not None and g0.get("verdict") == "STOP" and g0.get("commit") == commit:
        return {"void": False, "why": list(g0.get("why", [])), "commit": commit, "row": TRY_ID, "g0": g0,
                "outcome": "G0_STOP", "outcome_text": br.G0_STOP_TEXT, "passes": False, "files_passing": [],
                "tries_counted": 0}
    why = check_run(recs, why_in, commit)
    if g0 is None or g0.get("verdict") != "PASS" or g0.get("commit") != commit:
        why.append("G0's record is missing or not PASS at the run commit")
    ref_figs, ref_why = committed_figures() if ref is None else ref
    worst, so_report = {}, {}
    if not why:
        gw, worst, so_report = guard_problems(recs, ref_figs)
        why += [*ref_why, *gw]
    started = bool(run_started) or any(recs.get(t) for t in TASKS)
    if why:
        return {"void": True, "why": why, "commit": commit, "row": TRY_ID, "outcome": "VOID",
                "outcome_text": br.VOID_TEXT, "passes": False, "files_passing": [],
                "tries_counted": len(ARMS) if started else 0, "guard": worst}
    m = merges()
    e = clause_e(recs, m)
    arms, figs = {}, {}
    for arm in FILE_ARMS:
        figs[arm] = arm_figures(recs, m, arm)
        arms[arm] = decide_arm(figs[arm], e["holds"])
    key = outcome_of({a: arms[a]["passes"] for a in FILE_ARMS}, {a: arms[a]["idle"] for a in FILE_ARMS})
    files = list(br.PASS_FILES.get(key, ()))
    per_task = {}
    for t in TASKS:
        per_task[t] = {LABEL[a]: {"board": float(np.mean(series(recs[t], a, "board")
                                                         - series(recs[t], BASE_OF[a], "board"))),
                                  "shape_sum": float(np.mean(series(recs[t], a, "shape_sum")
                                                             - series(recs[t], BASE_OF[a], "shape_sum")))}
                       for a in ARMS}
    report = {
        LABEL[a]: {"merged_board": _stat(figs[a]["board"]), "prefix_board": _stat(figs[a]["board_prefix"]),
                   "shape_sum": _stat(figs[a]["shape_sum"]), "shape_sum_truths": figs[a]["shape_sum_truths"],
                   "nbhd": _stat(figs[a]["nbhd"]), "nbhd_truths": figs[a]["nbhd_truths"],
                   "e80_board": figs[a]["e80_board"], "e80_shape_sum": figs[a]["e80_shape_sum"],
                   "board_truths": per_truth(recs, m["five"], a, BASE_OF[a], "board"),
                   "groups": {g: _stat(merged_D(recs, m["five"], a, BASE_OF[a], f"group:{g}"))
                              for g in (*OTHER_GROUPS, "shape_growth")},
                   "test_like_board": _stat(merged_D(recs, m["test_like"], a, BASE_OF[a], "board")),
                   "four_truth_board": _stat(merged_D(recs, m["four"], a, BASE_OF[a], "board")),
                   "against_floor": _stat(merged_D(recs, m["five"], a, harness.FLOOR, "board"))}
        for a in FILE_ARMS}
    report["controls"] = {LABEL[c]: {"board_vs_O": _stat(merged_D(recs, m["five"], c, O, "board")),
                                     "shape_sum_vs_O": _stat(merged_D(recs, m["five"], c, O, "shape_sum"))}
                          for c in CONTROLS_E}
    report["against_floor"] = {LABEL[n]: _stat(merged_D(recs, m["five"], n, harness.FLOOR, "board"))
                               for n in (O, SO, *ARMS)}
    report["so_minus_s2xs_off_e51"] = so_report  # per task, deciding nothing (SO's guard is E51's construction checks)
    return {"void": False, "why": [], "commit": commit, "row": TRY_ID, "truths": m,
            "clause_text": CLAUSES, "pass_rule": PASS_RULE, "void_conditions": list(VOID_CONDITIONS),
            "arms": {LABEL[a]: arms[a] for a in FILE_ARMS}, "clause_e": {LABEL[c]: v for c, v in e["controls"].items()},
            "outcome": key, "outcome_text": br.OUTCOME_TEXTS[key], "passes": bool(files), "files_passing": files,
            "tries_counted": len(ARMS) if started else 0, "reported_deciding_nothing": report, "per_task": per_task,
            "price_lines": list(PRICE_LINES), "g0": g0, "guard": worst}


def text(res: dict) -> str:
    """The reading as text; its last line is the outcome line, in the row's fixed words."""
    L = [f"Row {TRY_ID} (run commit {res['commit']}): the gate of row {br.RULING_ROW}'s files {br.FILE_O} (BF_O) and "
         f"{br.FILE_SO} (BF_SO); it moves the overlay's form, the occupancy term and the z extent; "
         f"{res.get('tries_counted', 0)} counted tries; a diagnostic (J5, noise floor 6.54)."]
    if res["outcome"] == "G0_STOP":
        L.append("G0 (target-free, at the test pair's inputs, before any try) said STOP: " + "; ".join(res["why"]))
    elif res["void"]:
        L.append("VOID: " + "; ".join(res["why"]))
    else:
        L.append("The pass rule, fixed in row E59 before any output: " + res["pass_rule"])
        for a, v in res["arms"].items():
            L.append(f"  {a}: " + ", ".join(f"({k}) {'holds' if h else 'fails'}" for k, h in v["clauses"].items())
                     + ("; idle (below 0.3 in size)" if v["idle"] else "") + f"; {'PASSES' if v['passes'] else 'STOP'}")
        rep = res["reported_deciding_nothing"]
        for a in ("BF_O", "BF_SO"):
            r = rep[a]
            L.append(f"  {a} against its base: merged board {r['merged_board']['mean']:+.3f} (SE "
                     f"{r['merged_board']['se']:.3f}, {r['merged_board']['ahead']} of 10 ahead), prefix "
                     f"{r['prefix_board']['mean']:+.3f}; d2+Dice {r['shape_sum']['mean']:+.4f} ("
                     + ", ".join(f"{k.split('/')[-1]} {v:+.4f}" for k, v in r["shape_sum_truths"].items())
                     + f"); neighbourhood {r['nbhd']['mean']:+.4f}; E8.0 board {r['e80_board']:+.3f}, d2+Dice "
                     f"{r['e80_shape_sum']:+.4f}; test-like board {r['test_like_board']['mean']:+.3f}; four-truth "
                     f"board {r['four_truth_board']['mean']:+.3f}; against the floor {r['against_floor']['mean']:+.3f}")
        for c, v in res["clause_e"].items():
            L.append(f"  (e) BF_O's d2+Dice minus {c}'s {v['mean']:+.4f} ({v['ahead']} of 10 ahead): "
                     f"{'holds' if v['holds'] else 'fails'}")
        L.extend(res["price_lines"])
    L.append("Nothing here shows E7.75's form, the union's form at the test count, the live gate, scorer build or "
             "metric seed.")
    L.append(res["outcome_text"])
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
            "row_sha256": _sha_file(Path(br.__file__)), "harness_sha256": _sha_file(Path(harness.__file__)),
            "t2_build_sha256": _sha_file(Path(tb.__file__)), "shape_interp_sha256": _sha_file(Path(si.__file__))}


def cmd_read(a) -> None:
    out = Path(a.out)
    files = [out / "reading.json", out / "reading.txt"]
    existing = [str(f) for f in files if f.exists()]
    if existing:
        sys.exit(f"refusing to overwrite: {', '.join(existing)}")
    _own_code()
    why = br.row_problems(TRIES)
    if why:
        sys.exit("E59 has no well-formed row in experiments/harness/TRIES.md: " + "; ".join(why))
    runs = Path(a.runs)
    g0_path = runs / "g0" / "g0.json"
    g0 = json.loads(g0_path.read_text()) if g0_path.exists() else None
    why = commit_problems(a.commit)
    if repo_state(REPO).unclean != ():
        why.append("the tree is not clean at the read")
    if why:  # a read on another commit or tree reads nothing
        g0 = None if g0 is None or g0.get("verdict") != "STOP" else {**g0, "verdict": "UNREAD"}
    res = {**reading(load_runs(runs), why, a.commit, g0, run_started=a.run_started), "reading_code": provenance()}
    out.mkdir(parents=True, exist_ok=True)
    files[0].write_text(json.dumps(res, indent=1, default=float))
    files[1].write_text(text(res))
    print(f"wrote {files[1]}", flush=True)


# ------------------------------------------------------------------------------------------------ J7 dry run
def j7_dry_run(e55: Path | None = None, e58: Path | None = None, check_sha: bool = True) -> dict:
    """The O guard's tolerance tried on committed reproduction data other than the pair the gate reads (rule 9, J7),
    from committed records only (not a try): per task E58 ran, E58's committed O-minus-floor ten-seed mean (its run at
    0a6d4ce) against E55's (its run at 4535f98). The guard could fail there: two runs on two commits."""
    why = []
    paths = {"e55": Path(e55) if e55 else COMMITTED["e55"][0], "e58": Path(e58) if e58 else COMMITTED["e58"][0]}
    recs = {}
    for k, p in paths.items():
        if not p.exists() or (check_sha and _sha_file(p) != COMMITTED[k][1]):
            why.append(f"{k.upper()}'s committed reading is missing or changed")
            continue
        recs[k] = json.loads(p.read_text())
        if recs[k].get("void") is not False or recs[k].get("commit") != COMMITTED[k][2]:
            why.append(f"{k.upper()}'s reading is void or not from its run")
    devs = {}
    for name in sorted((recs.get("e58") or {}).get("tasks") or {}):
        a = (((recs.get("e58") or {}).get("tasks") or {}).get(name) or {}).get("board_minus_floor", {}).get(O)
        b = (((recs.get("e55") or {}).get("tasks") or {}).get(name) or {}).get("board_minus_floor", {}).get(O)
        if a is None or b is None:
            why.append(f"{name}: a committed figure is missing")
            continue
        devs[name] = abs(float(a) - float(b))
    worst = max(devs.values()) if devs else None
    ok = not why and len(devs) >= 7 and worst is not None and worst <= GUARD_TOL
    return {"pass": bool(ok), "tolerance": GUARD_TOL, "worst": worst, "per_task": devs, "why": why,
            "what": "per task |E58's O minus the floor - E55's O minus the floor|, ten-seed means under build main, "
                    "from committed records; not a try"}


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
    g = sub.add_parser("g0")
    g.add_argument("--seed", required=True, type=int)
    g.add_argument("--commit", required=True)
    g.add_argument("--out", default=str(RUNS / "g0"))
    g.add_argument("--data", default=str(harness.DATA_DIR))
    gr = sub.add_parser("g0read")
    gr.add_argument("--commit", required=True)
    gr.add_argument("--out", default=str(RUNS / "g0"))
    p = sub.add_parser("run")
    p.add_argument("--task", required=True, help=f"one of {', '.join(TASKS)}")
    p.add_argument("--commit", required=True, help="the run commit (full hash)")
    p.add_argument("--runs", default=str(RUNS))
    p.add_argument("--data", default=str(harness.DATA_DIR))
    r = sub.add_parser("read")
    r.add_argument("--commit", required=True)
    r.add_argument("--runs", default=str(RUNS))
    r.add_argument("--out", default=str(RESULTS))
    r.add_argument("--run-started", action="store_true",
                   help="the launcher's run step began (the four tries count from its start, rule 4)")
    j = sub.add_parser("j7")
    j.add_argument("--out", default=None, help="a JSON file for the dry run's per-task deviations")
    a = ap.parse_args(argv)
    code = {"g0": cmd_g0, "g0read": cmd_g0read, "run": cmd_run, "read": cmd_read, "j7": cmd_j7}[a.cmd](a)
    sys.exit(code or 0)


if __name__ == "__main__":
    main()
