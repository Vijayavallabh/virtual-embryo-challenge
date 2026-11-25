#!/usr/bin/env python
"""Row E57 of experiments/harness/TRIES.md (one row, one script): "RS: one run per side", the score-blind ERROR
designer's proposal P1 of the embryo design round of 2026-11-24/25 (round_1125_gap), as its cold challenger fixed it
(EMBRYO_ERROR_challenge.md, sections 1 to 7; the ledger row text of its section 7 is this script's specification). No
number here comes from a test score. All three tries are DIAGNOSTICS under rule 9 J5.

  E57a  step 0, a counted reading by the strict reading of rows N59 to N61 as row E56 applied it (one embryo try,
        counted when step 0 starts). The test pair decides; the nine practice pairs are reported.
  E57b  fgw_shape0x_ov_rs       RS, the arm (vec.harness.RS_CANDIDATES; vec.run_sort mode run_sorted)
  E57c  fgw_shape0x_ov_rs_same  its placebo (mode run_sorted_same)
E57b and E57c are counted when the run step starts, and never run unless step 0 passed. Control, rescored and not
counted: O = fgw_shape0x_ov (E49c). The step itself (the arm and the placebo) is vec.run_sort's, whose docstring gives
it word for word with every reading the row left open.

Step 0 (E57a), build seeds 1,000,000-1,000,009, on O's cloud (vec.t2_build.interpolate_clouds with t_pos 0,
floor_exact positions and the overlay form, at size 1: every figure below is scale-free), each figure per build seed and
then the median over the ten seeds (the own-run share: the mean over the seeds). The test pair is E7.5 and E8.0 at
t = 1/2 with the panel and count rule of T2:embryo:val_interp (never a test board); a practice pair is a task of
experiments/harness/nine_tasks.py with its own t, panel and count rule, the target term of a twin board's count left out
as experiments/coords/osr_gate.py leaves it (no target is opened, loaded or stat'ed). In the row's words:
  (i)   contrast (best score over median) at most 0.70 for the earlier bracket's rows and 0.85 for the later's at that
        plane;
  (ii)  median SD(f) of E57b at least 1.5 x O's (N15's f) and its mean own-run share at least 0.85;
  (iii) E57b's median homology error at most 1 x the median 15-nearest-cell reach, and E57c's at most E57b's.
  Any fail: STOP; E57b and E57c never run. A test pair that cannot be read in full also reads STOP.

The run: the nine tasks of nine_tasks.all_tasks(), seeds 0-9, bake protocol, build seed 1,000,000 + seed, scorer build
main (prefix and the multiplicative gate reported, deciding nothing). Merge (J4): five truths at equal weight (E8.0 =
half, late; E8.25_late = new1, x1; E8.5 = new2, x2, x3; E7.5 = e75; E7.25 = e725); D_s = seed s's merged difference.
  Gate: O within 0.01 board points of its committed per-task means (E55) before the arms: per task, the same-run
        ten-seed mean of O minus the floor under main against E55's committed tasks[name].board_minus_floor
        ['fgw_shape0x_ov'] (E55's reading.json, pinned by its sha256 E55_READING_SHA256 and its run commit).
  Void: a missing task, seed or arm; mixed commits or an unclean tree; an arm whose coordinates, count or sorted
        expression rows differ from O's; raw d2_shape, occupancy_dice or scale_log_ratio more than 1e-9 off O's; raw
        de_score or de_direction more than 1e-4 off O's. In the row's words, the void clause is fixed at 1e-4 for
        these two: "the absolute difference in de_score and in de_direction against O's in the same task and seed,
        which covers float32 summation-order noise" (the scorer's pseudobulk is np.mean(X, 0) of float32 rows, whose
        rounding follows the row order; on synthetic arrays a re-assignment of rows moved de_direction by up to 1.7e-5,
        tests/test_e57_rs_read.py). DE_TOL holds it; the row's first cell names it as DE_TOL_TEXT.
  E57b passes only if all hold:
    (a) merged board minus O at least +0.30, with D_s > 0 in at least 8 of 10 seeds;
    (b) local_organisation, plus the gated part of shape_growth where O's neighbourhood skill is below 0.5, carries at
        least 80% of (a), and cell_state's merged change is within 3 seed SEs of 0;
    (c) against E57c, local_organisation merged at least +0.01 with D_s > 0 in at least 8 of 10 seeds;
    (d) new2 and e75 merged board at least -0.30;
    (e) local_organisation merged above 0 in each bracket-run class (half, late, new1; x1, x2, x3).
  |E57b - O| under 0.30 merged reads STOP. A pass goes to the captain as a priced question with both readings (per
  truth, per class, half and x3 apart); no file without his ruling.

Readings where the row left a choice (declared, fixed before any data were read):
  * (b) in board points: 100 x 0.25 x (the arm's local_organisation group skill minus O's), plus, on each task-seed
    where O's neighborhood_mmd skill is below 0.5, 100 x 0.25 x (the arm's gated shape_growth group skill minus O's);
    merged as (a); the clause holds when that merged mean is at least 0.8 x (a)'s merged mean. cell_state's change is
    its group skill's; the SE is the SD of the ten D_s (ddof 1) over sqrt(10).
  * (c) and (e) read the local_organisation group skill (neighborhood_mmd's skill, weight 1); a class is merged by
    truth as the five-truth merge (class 1: E8.0 = half, late; E8.25_late = new1; class 2: E8.25_late = x1; E8.5 = x2,
    x3); (d) merges new2 (E8.5) and e75 (E7.5) at equal weight.
  * "at least" is >=, "above" and "D_s > 0" are strict; the idle case is |merged mean| < 0.30.
  * Sorted expression rows: the rows of X (float32) in lexicographic order, compared byte for byte; coordinates
    byte for byte in file order (float32); both from the files the scorer reads, under both builds; the 1e-9 and
    1e-4 checks under both builds, each skill against O's of the same task, seed and build.
  * The multiplicative gate: the board recomputed from each file's skills by vec.evaluate.aggregate(..., gate=
    "multiplicative") (shape times neighbourhood skill), reported only.
  * Per class (reported with a pass, and always): bracket-run class 1 (half, late, new1), class 2 (x1, x2, x3), the
    own-profile pair (new2, e75) and e725, each merged by truth; half and x3 (the full dose) apart.

Steps (each refuses unless TRIES.md holds exactly one row whose first cell starts "E57 (" and names E57a, this script,
its launcher and the de tolerance as DE_TOL_TEXT, the tracked tree is clean, vec loads from this checkout's src/ and
HEAD is the run commit):
  step0 --name NAME --commit SHA       one pair (`test` or a task), ten seeds, into <runs>/step0/NAME.json
  step0-read --commit SHA --out DIR    the ten pair files into DIR/reading.json and reading.txt (never overwritten);
                                       exit 0 on PASS, 3 on STOP
  gate --task NAME --commit SHA        the floor and O under main, ten seeds, against E55; exit 3 on a failure
  run --task NAME --commit SHA         the floor, O, the arm and the placebo, both builds, seed by seed (resumable);
                                       refused unless step 0's reading says PASS at this commit and all nine gate
                                       files pass at this commit
  read --commit SHA --out DIR          the reading (never overwritten), VOID when step 0, the gate or the run is not
                                       complete and clean at the run commit
The script prints only "<name> seed <s> done", "wrote <path>", "gate <task> pass|fail" and "step0 PASS|STOP". It reads
no test board, test file or test anchor, and prints and writes no cell count.
"""

from __future__ import annotations

import argparse
import contextlib
import dataclasses
import hashlib
import json
import math
import os
import subprocess
import sys
import traceback
from dataclasses import replace
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO / "experiments" / "harness"))
sys.path.insert(0, str(REPO / "src"))  # this checkout's code, not the checkout a shared venv's vec.pth names
import e47_floor_exact as e47  # noqa: E402
import nine_tasks as nt  # noqa: E402

from vec import evaluate as ev  # noqa: E402
from vec import harness  # noqa: E402
from vec import run_sort as rs  # noqa: E402
from vec import t2_build as tb  # noqa: E402
from vec.adoption import SEEDS, keep_change  # noqa: E402
from vec.boards import GROUPS  # noqa: E402
from vec.provenance import repo_state  # noqa: E402

ROW = "E57"
SCRIPT = "experiments/coords/e57_rs_read.py"
LAUNCHER = "experiments/coords/run_e57.sh"
TRIES = REPO / "experiments" / "harness" / "TRIES.md"
ROW_LINE = f"| {ROW} "  # the row's line in TRIES.md starts with this; its first cell starts with "E57 ("
DE_TOL_TEXT = "de_score and de_direction more than 1e-4 off O's"  # the void clause's tolerance, as the row names it
# the first cell must name step 0's ID, the script, its launcher and the de tolerance
ROW_NAMES = ("E57a", SCRIPT, LAUNCHER, DE_TOL_TEXT)
O = "fgw_shape0x_ov"  # noqa: E741  (the row names the control O)
ARM = "fgw_shape0x_ov_rs"
PLACEBO = "fgw_shape0x_ov_rs_same"
COUNTED = (ARM, PLACEBO)
TRY_IDS = {"step0": "E57a", ARM: "E57b", PLACEBO: "E57c"}
MODE_OF = {O: "as_built", ARM: "run_sorted", PLACEBO: "run_sorted_same"}
RUN_NAMES = (O, ARM, PLACEBO)  # with the floor, in each run_task call
BUILDS = ("main", "prefix")
FILE = "T2embryo_s2x_ovrs.h5ad"
CANDIDATE_COMMAND = ("experiments/builds/build_test.py embryo --method fgw_shape0 --positions floor_exact "
                     "--expression mixture --shape-interp overlay --shape-interp-test-ok --neighbourhood run_sorted "
                     f"--captain-ruling <ROW> --out runs/builds/{FILE}")
LABEL = "diagnostic (J5)"
# step 0
TEST = "test"
TASK_NAMES = ("half", "late", "new1", "new2", "e75", "e725", "x1", "x2", "x3")
NAMES = (TEST, *TASK_NAMES)
STEP0_SEEDS = tuple(harness.BUILD_SEED_OFFSET + s for s in SEEDS)  # 1,000,000 to 1,000,009
TEST_BOARD = "T2:embryo:val_interp"  # the test pair's panel and count rule: a validation board, never a test board
TEST_PAIR = {"before": "T2/E7.5.h5ad", "after": "T2/E8.0.h5ad"}
TEST_T = 0.5
STEP0_RMS = 1.0  # the size rule's scale: every step-0 figure is scale-free
CONTRAST_EARLIER, CONTRAST_LATER = 0.70, 0.85  # (i)
SD_F_RATIO, OWN_RUN = 1.5, 0.85  # (ii)
HOMOLOGY_REACH = 1.0  # (iii)
OUTCOME_PASS = "E57a step 0: PASS"
# the gate
E55_DIR = REPO / "experiments" / "coords" / "results" / "e55_overlay_stretch"
E55_READING_SHA256 = "e500c9fa4161720191562598af76f8fb28bc12b0c3b0dfa9fdd3e38a38153ab5"  # its reading.json
E55_RUN_COMMIT = "4535f98c98ae1a2a8488185e8a99e4a2e8d92d86"
O_TOL = 0.01  # board points: O's ten-seed mean minus the floor's against E55's, per task
# the read
SHAPE_METRICS = ("d2_shape", "occupancy_dice", "scale_log_ratio")
DE_METRICS = ("de_score", "de_direction")
SHAPE_TOL = 1e-9  # raw d2_shape, occupancy_dice and scale_log_ratio against O's: the coordinates are O's bytes
DE_TOL = 1e-4  # raw de_score and de_direction against O's, same task and seed: float32 summation-order noise
BAR, AHEAD = 0.30, 8  # (a) and the idle case
CARRY = 0.8  # (b)
CELL_SE = 3.0  # (b)
PLACEBO_BAR = 0.01  # (c)
OWN_PROFILE_FLOOR = -0.30  # (d)
NB_GATE = 0.5  # (b): where O's neighbourhood skill is below this, the gated part of shape_growth counts
G = "local_organisation"
WEIGHT = {name: w for name, w, _ in GROUPS["T2"]}
CLASSES = {"class 1 (half, late, new1)": ("half", "late", "new1"), "class 2 (x1, x2, x3)": ("x1", "x2", "x3")}
OTHER_CLASSES = {"own profile (new2, e75)": ("new2", "e75"), "e725": ("e725",)}
SCRIPT_RUNS = REPO / "runs" / "e57"
RESULTS = REPO / "experiments" / "coords" / "results" / "e57_rs"
STEP0_RESULTS = REPO / "experiments" / "coords" / "results" / "e57a_step0"
CLAUSES = {
    "a": "(a) merged board minus O at least +0.30, with D_s > 0 in at least 8 of 10 seeds",
    "b": "(b) local_organisation, plus the gated part of shape_growth where O's neighbourhood skill is below 0.5, "
         "carries at least 80% of (a), and cell_state's merged change is within 3 seed SEs of 0",
    "c": "(c) against E57c, local_organisation merged at least +0.01 with D_s > 0 in at least 8 of 10 seeds",
    "d": "(d) new2 and e75 merged board at least -0.30",
    "e": "(e) local_organisation merged above 0 in each bracket-run class (half, late, new1; x1, x2, x3)",
}
STEP0_CLAUSES = {
    "i": "(i) contrast (best score over median) at most 0.70 for the earlier bracket's rows and 0.85 for the later's "
         "at that plane",
    "ii": "(ii) median SD(f) of E57b at least 1.5 x O's (N15's f) and its mean own-run share at least 0.85",
    "iii": "(iii) E57b's median homology error at most 1 x the median 15-nearest-cell reach, and E57c's at most "
           "E57b's",
}
VOID_TEXT = ("VOID = read as nothing and no file. No task or seed is re-run selectively, and nothing further follows "
             "under this row. A rerun needs a new row placed before it, and its tries count again. Tries counted at "
             "the start of step 0 and of the run step stay counted (rule 4).")
STOP_TEXT = "STOP: no file; the line closes and the handed-over embryo files stand."


# ------------------------------------------------------------------------------------------------ hooks (tests swap)
def _load(path: Path, genes: list[str]):
    from vec.prep import load_stage

    return load_stage(path, genes=genes, drop_duplicates=True)


def _rows_as_released(path: Path, genes: list[str]) -> int:
    from vec.io import read_stage

    return int(read_stage(path, genes).n_obs)


def _boards():
    from vec.boards import load_boards

    return load_boards()


def _tasks() -> dict:
    return nt.all_tasks()


# ------------------------------------------------------------------------------------------------ checks
def _git(*args: str) -> subprocess.CompletedProcess:
    return subprocess.run(["git", "-C", str(REPO), *args], capture_output=True, text=True,
                          env={**os.environ, "GIT_OPTIONAL_LOCKS": "0"})


def row_ok(path: Path | None = None) -> bool:
    """Whether TRIES.md holds exactly one row E57 whose first cell starts "E57 (" and names ROW_NAMES."""
    p = TRIES if path is None else Path(path)
    if not p.exists():
        return False
    lines = [ln for ln in p.read_text(encoding="utf-8").splitlines() if ln.startswith(ROW_LINE)]
    if len(lines) != 1:
        return False
    first = lines[0][2:].split("|", 1)[0]
    return first.startswith(f"{ROW} (") and all(name in first for name in ROW_NAMES)


def commit_problems(commit: str) -> list[str]:
    head = _git("rev-parse", "HEAD").stdout.strip()
    return [] if head == commit else [f"HEAD is {head[:12] or 'unknown'}, not the run commit {commit[:12]}"]


def own_code_problems() -> list[str]:
    src = (REPO / "src").resolve()
    return [f"{m.__name__} is loaded from {m.__file__}, not this checkout's src/" for m in (ev, harness, rs, tb)
            if src not in Path(m.__file__).resolve().parents]


def ready_problems(commit: str) -> list[str]:
    why = own_code_problems()
    if not row_ok():
        why.append(f"experiments/harness/TRIES.md holds no single row {ROW} whose first cell starts \"{ROW} (\" and "
                   f"names {', '.join(ROW_NAMES)}: the row is committed before its run")
    if repo_state(REPO).unclean != ():
        why.append("the tracked tree is not clean")
    return [*why, *commit_problems(commit)]


def _refuse_unless_ready(commit: str) -> None:
    why = ready_problems(commit)
    if why:
        sys.exit("refused: " + "; ".join(why))


def _is_test_board(board) -> bool:
    return str(board.split).startswith("test") or "test" in str(board.key).split(":")[-1]


def _sha_file(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _median(vals) -> float | None:
    v = [float(x) for x in vals if x is not None and np.isfinite(x)]
    return float(np.median(v)) if v else None


# ------------------------------------------------------------------------------------------------ step 0 (E57a)
def pair_inputs(name: str, data: Path) -> dict:
    """The pair's two brackets as the build loads them (never a target), its t and its count."""
    from vec.io import COORD_KEY, dense

    boards = _boards()
    if name == TEST:
        board, files, t = boards[TEST_BOARD], dict(TEST_PAIR), TEST_T
        panel_board = TEST_BOARD
        count_rule = f"harness.count_for({TEST_BOARD}, before's rows after drop_duplicates)"
    else:
        task = _tasks()[name]
        panel_board = task.panel_board
        files = {r: task.inputs[r] for r in ("before", "after")}
        view = harness.TaskView.of(task)
        lo, hi = view.role_days["before"], view.role_days["after"]
        t = (view.target_day - lo) / (hi - lo)
        if task.validation_board:
            board = replace(boards[task.validation_board], anchors={})
            count_rule = f"harness.count_for({task.validation_board} as published, before's rows)"
        else:
            twin = boards[task.twin_board]
            n_ref = _rows_as_released(data / task.reference, list(twin.genes))
            pc = harness.pseudo_counts(0, n_ref)  # the target term left out, as osr_gate.py
            board = replace(twin, key=f"pseudo:{task.key}", split="pseudo", anchors={}, max_cells=pc["max_cells"],
                            min_cells=1, truth_cells=None, ref_cells=pc["ref_cells"])
            count_rule = (f"harness.count_for(pseudo board on {task.twin_board}, before's rows), max_cells = "
                          "max(5000, 3 ref_cells): the target term left out")
    if _is_test_board(board) or board.counts_hidden:
        raise PermissionError("step 0 reads a validation or pseudo board's panel and count rule only")
    genes = list(board.genes)
    st = {r: _load(data / f, genes) for r, f in files.items()}
    for r, a in st.items():
        missing = [g for g in genes if g not in a.var_names]
        if missing:
            raise ValueError(f"{name}: input {r} lacks {len(missing)} panel genes")
    return {"before": st["before"], "after": st["after"],
            "Cb": np.asarray(st["before"].obsm[COORD_KEY], float)[:, :3],
            "Ca": np.asarray(st["after"].obsm[COORD_KEY], float)[:, :3],
            "Eb": dense(st["before"][:, genes].X).astype(np.float32),
            "Ea": dense(st["after"][:, genes].X).astype(np.float32),
            "t": float(t), "n": harness.count_for(board, st["before"].n_obs), "files": files,
            "panel_board": panel_board, "count_rule": count_rule}


def seed_figures(p: dict, seed: int) -> dict:
    """One pair-seed: O's cloud, the arm's and the placebo's step on it, and step 0's figures (no count)."""
    res = tb.interpolate_clouds(p["Cb"], p["Ca"], p["Eb"], p["Ea"], p["t"], p["n"], seed, target_rms=STEP0_RMS,
                                t_pos=0.0, positions=tb.POSITIONS[1], shape_form="overlay")
    labels = tb.drawn_labels(res, p["before"].obs["celltype"], p["after"].obs["celltype"])
    groups = tb.drawn_labels(res, p["before"].obs["group"], p["after"].obs["group"])
    arm = rs.sort_rows(res.coords, res.stage, labels, groups, MODE_OF[ARM], later_contrast=True)
    pl = rs.sort_rows(res.coords, res.stage, labels, groups, MODE_OF[PLACEBO])
    mix = {k: rs.bracket_mix(res.coords, s) for k, s in (("O", res.stage), ("RS", arm.stage_after),
                                                         ("placebo", pl.stage_after))}
    return {"contrast_earlier": arm.record["contrast_earlier"], "contrast_later": arm.record["contrast_later"],
            "normal": [float(x) for x in arm.normal],
            "sd_f": {k: v["sd_f"] for k, v in mix.items()},
            "own_run_share": {k: v["own_run_share"] for k, v in mix.items()},
            "homology_error": {"RS": arm.record["median_homology_error"],
                               "placebo": pl.record["median_homology_error"]},
            "reach": rs.reach(res.coords),
            "placebo_matches_swaps": bool(pl.record["placebo_matches_swaps"]),
            "placebo_short": bool(pl.record["placebo_shortfall"])}


def run_pair(name: str, data: Path, out: Path, commit: str) -> int:
    rec: dict = {"name": name, "commit": commit, "script": SCRIPT, "seeds": list(STEP0_SEEDS), "row": ROW,
                 "try": TRY_IDS["step0"]}
    code = 0
    try:
        p = pair_inputs(name, Path(data))
        rec.update(t=p["t"], files=p["files"], panel_board=p["panel_board"], count_rule=p["count_rule"])
        per_seed = {}
        for seed in STEP0_SEEDS:
            try:
                per_seed[str(seed)] = seed_figures(p, seed)
            except Exception as e:  # recorded for the seed; on the test pair it reads STOP
                per_seed[str(seed)] = {"error": f"{type(e).__name__}: {e}"}
                code = 3
        rec["per_seed"] = per_seed
    except Exception as e:
        rec["error"] = f"{type(e).__name__}: {e}"
        rec["traceback"] = traceback.format_exc(limit=5)
        code = 1
    (out / f"{name}.json").write_text(json.dumps(rec, indent=1))
    return code


def pair_summary(r: dict) -> dict:
    """A pair file's medians over its seeds (the own-run share: the mean), and whether every seed was read."""
    per = r.get("per_seed") or {}
    ok = [v for v in per.values() if "error" not in v]
    complete = (sorted(per) == sorted(str(s) for s in STEP0_SEEDS) and len(ok) == len(STEP0_SEEDS)
                and "error" not in r)

    def med(get):
        return _median([get(v) for v in ok])

    def mean(vals):  # correctly rounded, so ten equal shares average to that share
        return math.fsum(vals) / len(vals) if vals else None

    return {"complete": complete,
            "contrast_earlier": med(lambda v: v["contrast_earlier"]),
            "contrast_later": med(lambda v: v["contrast_later"]),
            "sd_f": {k: med(lambda v, k=k: v["sd_f"][k]) for k in ("O", "RS", "placebo")},
            "own_run_share": {k: mean([v["own_run_share"][k] for v in ok]) for k in ("O", "RS", "placebo")},
            "mean_own_run_share_rs": mean([v["own_run_share"]["RS"] for v in ok]),
            "homology_error": {k: med(lambda v, k=k: v["homology_error"][k]) for k in ("RS", "placebo")},
            "reach": med(lambda v: v["reach"]),
            "placebo_matches_swaps": bool(ok) and all(v["placebo_matches_swaps"] for v in ok),
            "errors": sorted({v["error"] for v in per.values() if "error" in v} | ({r["error"]} if "error" in r
                                                                                    else set()))}


def step0_decide(s: dict) -> dict:
    """The row's step-0 clauses on the test pair's summary: which hold."""
    def le(x, bar):
        return x is not None and bar is not None and x <= bar

    i = le(s["contrast_earlier"], CONTRAST_EARLIER) and le(s["contrast_later"], CONTRAST_LATER)
    sd_o, sd_rs = s["sd_f"]["O"], s["sd_f"]["RS"]
    ii = (sd_o is not None and sd_rs is not None and sd_rs >= SD_F_RATIO * sd_o
          and s["mean_own_run_share_rs"] is not None and s["mean_own_run_share_rs"] >= OWN_RUN)
    h_rs, h_pl, reach = s["homology_error"]["RS"], s["homology_error"]["placebo"], s["reach"]
    iii = (h_rs is not None and reach is not None and h_rs <= HOMOLOGY_REACH * reach and le(h_pl, h_rs))
    return {"i": bool(i), "ii": bool(ii), "iii": bool(iii)}


def step0_reading(out: Path) -> dict:
    """The ten pair files: the test pair decides (i) to (iii); the nine practice pairs are reported."""
    pairs, commits, missing = {}, set(), []
    for name in NAMES:
        p = Path(out) / f"{name}.json"
        try:
            r = json.loads(p.read_text())
        except Exception as e:
            missing.append(f"{name}: no readable pair file ({type(e).__name__})")
            continue
        commits.add(r.get("commit"))
        pairs[name] = pair_summary(r)
    fails, notes = [], []
    clauses = {"i": False, "ii": False, "iii": False}
    tp = pairs.get(TEST)
    if tp is None or not tp["complete"]:
        fails.append("incomplete")
        notes.append("the test pair was not read in full: " + "; ".join(
            (tp or {}).get("errors") or [m for m in missing if m.startswith(TEST)] or ["missing seeds"]))
    else:
        clauses = step0_decide(tp)
        fails += [k for k in ("i", "ii", "iii") if not clauses[k]]
    outcome = OUTCOME_PASS if not fails else f"E57a step 0: STOP ({', '.join(fails)})"
    return {"outcome": outcome, "pass": not fails, "clauses": clauses, "clause_text": STEP0_CLAUSES,
            "test_pair": tp, "practice_pairs_reported": {n: pairs.get(n) for n in TASK_NAMES},
            "unreadable": missing, "notes": notes, "pair_commits": sorted(str(c) for c in commits),
            "seeds": list(STEP0_SEEDS), "test_pair_board": TEST_BOARD, "script": SCRIPT}


def _fmt(x, f=".3f") -> str:
    return "none" if x is None else format(x, f)


def step0_text(res: dict) -> str:
    L = [res["outcome"],
         f"Row {ROW}, step 0 ({TRY_IDS['step0']}; build seeds {STEP0_SEEDS[0]} to {STEP0_SEEDS[-1]}; the test pair "
         f"E7.5 + E8.0 at t = 1/2 on {TEST_BOARD}'s panel and count rule decides; medians over the seeds):"]
    for k, text in res["clause_text"].items():
        L.append(f"  {k} {'holds' if res['clauses'][k] else 'fails'}: {text}")
    for n in res["notes"]:
        L.append("  " + n)
    rows = [(TEST, res["test_pair"])] + list(res["practice_pairs_reported"].items())
    L.append("Per pair (the test pair decides; the nine practice pairs are reported, deciding nothing): contrast "
             "earlier/later; SD(f) O/RS/placebo; own-run share O/RS/placebo; homology error RS/placebo; reach")
    for name, s in rows:
        if s is None:
            L.append(f"  {name:<5} missing")
            continue
        L.append(f"  {name:<5} {_fmt(s['contrast_earlier'])}/{_fmt(s['contrast_later'])}; "
                 + "/".join(_fmt(s["sd_f"][k]) for k in ("O", "RS", "placebo")) + "; "
                 + "/".join(_fmt(s["own_run_share"][k]) for k in ("O", "RS", "placebo")) + "; "
                 + "/".join(_fmt(s["homology_error"][k], ".4f") for k in ("RS", "placebo"))
                 + f"; {_fmt(s['reach'], '.4f')}"
                 + ("" if s["complete"] else " (incomplete)")
                 + ("" if s["placebo_matches_swaps"] else " (placebo swap count short)"))
    L.append("Any fail: STOP; E57b and E57c never run." if not res["pass"]
             else "PASS: the run step may start (E57b and E57c are counted when it starts).")
    p = res.get("provenance", {})
    L.append(f"Code: commit {p.get('commit')}, clean {p.get('clean')}, script sha256 {p.get('script_sha256')}, "
             f"run_sort sha256 {p.get('run_sort_sha256')}")
    return "\n".join(L) + "\n"


def provenance() -> dict:
    s = repo_state(REPO)
    return {"commit": s.commit, "clean": s.unclean == (), "script_sha256": _sha_file(Path(__file__)),
            "run_sort_sha256": _sha_file(Path(rs.__file__)), "harness_sha256": _sha_file(Path(harness.__file__)),
            "t2_build_sha256": _sha_file(Path(tb.__file__))}


def step0_passed(folder: Path, commit: str) -> list[str]:
    """Why step 0's reading in `folder` does not let the run start at `commit` (empty: it does)."""
    p = Path(folder) / "reading.json"
    if not p.exists() or not (Path(folder) / "reading.txt").exists():
        return [f"no step-0 reading in {folder}"]
    r = json.loads(p.read_text())
    first = (Path(folder) / "reading.txt").read_text().split("\n", 1)[0]
    pr = r.get("provenance") or {}
    why = []
    if r.get("outcome") != OUTCOME_PASS or first != OUTCOME_PASS:
        why.append("step 0 did not pass")
    if pr.get("commit") != commit or pr.get("clean") is not True:
        why.append("step 0's reading is not from the run commit on a clean tree")
    if r.get("pair_commits") != [commit]:
        why.append("step 0's pair files are not all from the run commit")
    return why


def cmd_step0(a) -> int:
    out = Path(a.runs) / "step0"
    out.mkdir(parents=True, exist_ok=True)
    if (out / f"{a.name}.json").exists():
        print(f"refused: {out / (a.name + '.json')} exists", file=sys.stderr)
        return 2
    ex = out / f"{a.name}.exit"
    ex.write_text("running\n")
    why = ready_problems(a.commit)
    if why:
        print("refused: " + "; ".join(why), file=sys.stderr)
        code = 2
    else:
        code = run_pair(a.name, Path(a.data), out, a.commit)
    ex.write_text(f"{code}\n")
    return code


def cmd_step0_read(a) -> int:
    out = Path(a.out)
    files = [out / "reading.json", out / "reading.txt"]
    existing = [str(f) for f in files if f.exists()]
    if existing:
        print(f"refusing to overwrite: {', '.join(existing)}", file=sys.stderr)
        return 95
    _refuse_unless_ready(a.commit)
    res = {**step0_reading(Path(a.runs) / "step0"), "provenance": provenance(), "row": ROW}
    out.mkdir(parents=True, exist_ok=True)
    files[0].write_text(json.dumps(res, indent=1, default=float))
    files[1].write_text(step0_text(res))
    print(f"wrote {files[1]}", flush=True)
    print("step0 " + ("PASS" if res["pass"] else "STOP"), flush=True)
    return 0 if res["pass"] else 3


# ------------------------------------------------------------------------------------------------ the gate
def e55_o_per_task(path: Path | None = None) -> tuple[dict | None, list[str]]:
    """E55's committed per-task figure of O: the ten-seed mean of O minus the floor under main."""
    p = E55_DIR / "reading.json" if path is None else Path(path)
    if not p.exists():
        return None, ["no committed E55 reading"]
    if _sha_file(p) != E55_READING_SHA256:
        return None, [f"E55's reading.json is not the one of sha256 {E55_READING_SHA256[:12]}"]
    r = json.loads(p.read_text())
    why = []
    if r.get("void") is not False or r.get("commit") != E55_RUN_COMMIT:
        why.append(f"E55's reading is void or not from its run at {E55_RUN_COMMIT[:12]}")
    out = {}
    for name in TASK_NAMES:
        v = (((r.get("tasks") or {}).get(name) or {}).get("board_minus_floor") or {}).get(O)
        if v is None:
            why.append(f"{name}: no committed E55 figure for {O}")
        else:
            out[name] = float(v)
    return out, why


def gate_check(name: str, run_scores: dict, ref: tuple[dict | None, list[str]] | None = None) -> dict:
    """O's ten-seed mean minus the floor's (this run, main) within O_TOL of E55's committed per-task figure."""
    ref_map, why = e55_o_per_task() if ref is None else ref
    why = list(why)
    o, f = run_scores.get(O), run_scores.get(harness.FLOOR)
    off = None
    if o is None or f is None or np.shape(o) != np.shape(f) or np.shape(o) != (len(SEEDS),):
        why.append(f"{name}: {O} or the floor missing or of another length")
    elif ref_map is None or name not in ref_map:
        why.append(f"{name}: no committed E55 figure for {O}")
    else:
        mean = float(np.mean(np.asarray(o, float) - np.asarray(f, float)))
        off = abs(mean - ref_map[name])
        if not off <= O_TOL:
            why.append(f"{name}: {O} minus the floor is {off:.4f} off E55's ten-seed mean (tolerance {O_TOL})")
    return {"task": name, "tolerance": O_TOL, "o_off_e55": off, "e55_run_commit": E55_RUN_COMMIT,
            "e55_reading_sha256": E55_READING_SHA256, "why": why, "pass": not why}


def _task(name: str):
    tasks = _tasks()
    if name not in tasks:
        sys.exit(f"unknown task {name}; known: {sorted(tasks)}")
    t = tasks[name]
    if getattr(t, "target_filter", None) is not None:
        sys.exit(f"{name}: a task with a target filter; E57 does not run on one")
    return t


def _gate_path(runs: Path, name: str) -> Path:
    return Path(runs) / "gate" / f"{name}.json"


def cmd_gate(a) -> int:
    _refuse_unless_ready(a.commit)
    t = _task(a.task)
    path = _gate_path(Path(a.runs), a.task)
    if path.exists():
        sys.exit(f"{path} exists")
    res = harness.run_task(t, [harness.CANDIDATES[harness.FLOOR], harness.CANDIDATES[O]], seeds=SEEDS,
                           data_dir=Path(a.data), build="main", log=None)
    if str(res["board"].get("split", "")).startswith("test"):
        sys.exit("E57 met a test board")
    rec = {**gate_check(a.task, e47.per_seed(res)[0]), "commit": res["provenance"]["commit"],
           "tree_clean": res["provenance"]["tree_clean"]}
    rec["pass"] = bool(rec["pass"] and rec["commit"] == a.commit and rec["tree_clean"])
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(rec, indent=1, default=float))
    print(f"gate {a.task} {'pass' if rec['pass'] else 'fail'}", flush=True)
    return 0 if rec["pass"] else 3


def gate_problems(runs: Path, commit: str) -> list[str]:
    why = []
    for name in TASK_NAMES:
        p = _gate_path(runs, name)
        if not p.exists():
            why.append(f"no gate file for {name}")
            continue
        g = json.loads(p.read_text())
        if not g.get("pass") or g.get("commit") != commit:
            why.append(f"the gate did not pass on {name} at the run commit")
    return why


# ------------------------------------------------------------------------------------------------ the run
def coords_digest(stage) -> str:
    C = np.ascontiguousarray(stage.coords, dtype=np.float32)
    return hashlib.sha256(np.asarray(C.shape, np.int64).tobytes() + C.tobytes()).hexdigest()


def sorted_rows_digest(stage) -> str:
    """The expression rows as a multiset: X (float32) with its rows in lexicographic order, as one digest."""
    X = np.ascontiguousarray(stage.X, dtype=np.float32)
    order = np.lexsort(X.T[::-1]) if X.size else np.arange(len(X))
    S = np.ascontiguousarray(X[order])
    return hashlib.sha256(np.asarray(S.shape, np.int64).tobytes() + S.tobytes()).hexdigest()


def file_digest(stage) -> str:
    X = np.ascontiguousarray(stage.X, dtype=np.float32)
    C = np.ascontiguousarray(stage.coords, dtype=np.float32)
    return hashlib.sha256(X.tobytes() + C.tobytes()).hexdigest()


@contextlib.contextmanager
def captured():
    """While active, record what harness.score_seeds is handed (task, board) and each vec.run_sort.sort_rows call's
    record by mode, with the bracket mix before and after it (shares only, no count)."""
    seen: dict = {"steps": {}}
    orig_ss, orig_sr = harness.score_seeds, rs.sort_rows

    def score_seeds(t, board, target, ref, inputs, candidates, *a, **kw):
        seen.update(task=t, board=board)
        return orig_ss(t, board, target, ref, inputs, candidates, *a, **kw)

    def sort_rows(coords, stage, labels, groups, mode=rs.MODES[1], **kw):
        out = orig_sr(coords, stage, labels, groups, mode=mode, **kw)
        before, after = rs.bracket_mix(coords, stage), rs.bracket_mix(coords, out.stage_after)
        seen["steps"][mode] = {**{k: v for k, v in out.record.items() if k != "rule"},
                               "sd_f_o": before["sd_f"], "sd_f": after["sd_f"],
                               "own_run_share_o": before["own_run_share"], "own_run_share": after["own_run_share"]}
        return out

    harness.score_seeds, rs.sort_rows = score_seeds, sort_rows
    try:
        yield seen
    finally:
        harness.score_seeds, rs.sort_rows = orig_ss, orig_sr


def _keeping(base, seen, s: int):
    """The candidate as harness.run_task builds it, its finished file's digests kept."""
    def build(ctx):
        if _is_test_board(seen["board"]):
            raise RuntimeError("E57 met a test board")
        pred = base.build(ctx)
        seen.setdefault("digests", {})[base.name] = {"coords": coords_digest(pred), "rows": sorted_rows_digest(pred),
                                                     "file": file_digest(pred), "n": int(len(pred.X))}
        return pred

    return dataclasses.replace(base, build=build)


def one_seed(t, s: int, data: Path) -> dict:
    """Both builds of one task-seed: run_task's results and E57's per-seed facts (no count but equalities)."""
    out: dict = {"builds": {}, "facts": {}}
    digests, steps = {}, {}
    for build in BUILDS:
        with captured() as seen:
            cands = [harness.CANDIDATES[harness.FLOOR],
                     *(_keeping(harness.CANDIDATES[n], seen, s) for n in RUN_NAMES)]
            res = harness.run_task(t, cands, seeds=[s], data_dir=data, build=build, log=None)
        if _is_test_board(seen["board"]):
            raise RuntimeError("E57 met a test board")
        out["builds"][build] = res
        digests[build] = seen.get("digests", {})
        steps[build] = seen["steps"]
    d = digests[BUILDS[0]]
    out["facts"] = {
        "board_split": str(out["builds"][BUILDS[0]]["board"].get("split")),
        "same_coords": {n: all(digests[b][n]["coords"] == digests[b][O]["coords"] for b in BUILDS) for n in COUNTED},
        "same_count": {n: all(digests[b][n]["n"] == digests[b][O]["n"] for b in BUILDS) for n in COUNTED},
        "same_sorted_rows": {n: all(digests[b][n]["rows"] == digests[b][O]["rows"] for b in BUILDS) for n in COUNTED},
        "rows_moved": {n: d[n]["file"] != d[O]["file"] for n in COUNTED},
        "files_identical_across_builds": all(digests[BUILDS[0]][n]["file"] == digests[b][n]["file"]
                                             for b in BUILDS[1:] for n in RUN_NAMES),
        "steps": steps[BUILDS[0]]}
    return out


def cmd_run(a) -> int:
    _refuse_unless_ready(a.commit)
    t = _task(a.task)
    runs = Path(a.runs)
    why = [*step0_passed(Path(a.step0), a.commit), *gate_problems(runs, a.commit)]
    if why:
        sys.exit("refused (step 0 and the pre-run gate): " + "; ".join(why))
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
    return 0


# ------------------------------------------------------------------------------------------------ the read
def _row(r: dict, build: str = "main") -> dict:
    return r["builds"][build]["rows"][0]


def per_seed_value(r: dict, name: str, what: str, build: str = "main") -> float:
    """One file's quantity at one task-seed: 'board', 'mult' (the multiplicative gate's board), 'group:<g>' (gated)
    or 'skill:<m>'."""
    x = _row(r, build)
    if what == "board":
        return float(x["scores"][name])
    if what == "mult":
        return float(ev.aggregate("T2", x["skills"][name], gate="multiplicative")[0])
    if what.startswith("group:"):
        return float(x["groups"][name][what[6:]])
    if what.startswith("skill:"):
        return float(x["skills"][name][what[6:]])
    raise ValueError(what)


def series(recs: list[dict], name: str, what: str, build: str = "main") -> np.ndarray:
    return np.array([per_seed_value(r, name, what, build) for r in recs], float)


def carried(recs: list[dict], a: str = ARM, b: str = O) -> np.ndarray:
    """(b)'s carried part per seed, in board points: the local_organisation group's change, plus the gated
    shape_growth group's change where b's neighborhood_mmd skill is below NB_GATE."""
    lo = series(recs, a, f"group:{G}") - series(recs, b, f"group:{G}")
    sg = series(recs, a, "group:shape_growth") - series(recs, b, "group:shape_growth")
    gated = series(recs, b, "skill:neighborhood_mmd") < NB_GATE
    return 100.0 * WEIGHT[G] * lo + np.where(gated, 100.0 * WEIGHT["shape_growth"] * sg, 0.0)


def merged_D(recs: dict[str, list[dict]], groups: dict[str, list[str]], a: str, b: str, what: str,
             build: str = "main") -> np.ndarray:
    """D_s: the mean over truths of each truth's task-mean difference (a minus b) at each seed."""
    if what == "carried":
        per = {n: carried(recs[n], a, b) for names in groups.values() for n in names}
    else:
        per = {n: series(recs[n], a, what, build) - series(recs[n], b, what, build)
               for names in groups.values() for n in names}
    return e47.merged_seeds(per, groups)


def _stat(D) -> dict:
    D = np.asarray(D, float)
    return {"mean": float(D.mean()), "se": float(D.std(ddof=1) / np.sqrt(len(D))), "ahead": int((D > 0).sum()),
            "behind": int((D < 0).sum())}


def merges(tasks: dict | None = None) -> dict[str, dict[str, list[str]]]:
    """The merges the row fixes (five truths), its classes and the parts reported apart, each merged by truth."""
    tasks = _tasks() if tasks is None else tasks

    def sub(names):
        return e47.truths({n: tasks[n] for n in names})

    return {"five": e47.truths(tasks), "own_profile_d": sub(("new2", "e75")),
            **{k: sub(v) for k, v in {**CLASSES, **OTHER_CLASSES}.items()},
            "half": sub(("half",)), "x3": sub(("x3",))}


def decide(d: dict) -> dict:
    """The row's clauses on per-seed D_s arrays: board (five-truth board, arm minus O), carried ((b)'s part), cell
    (cell_state group), lo_pl (local_organisation, arm minus the placebo), own ((d)'s merge), classes ({class: the
    local_organisation D_s, arm minus O})."""
    board = np.asarray(d["board"], float)
    cell = np.asarray(d["cell"], float)
    se = float(cell.std(ddof=1) / np.sqrt(len(cell)))
    lo_pl = np.asarray(d["lo_pl"], float)
    c = {"a": bool(board.mean() >= BAR and (board > 0).sum() >= AHEAD),
         "b": bool(np.mean(d["carried"]) >= CARRY * board.mean() and abs(cell.mean()) <= CELL_SE * se),
         "c": bool(lo_pl.mean() >= PLACEBO_BAR and (lo_pl > 0).sum() >= AHEAD),
         "d": bool(np.mean(d["own"]) >= OWN_PROFILE_FLOOR),
         "e": bool(all(np.mean(v) > 0 for v in d["classes"].values()) and len(d["classes"]) == len(CLASSES))}
    idle = bool(abs(board.mean()) < BAR)
    return {"clauses": c, "idle": idle, "passes": bool(all(c.values()) and not idle)}


def check_run(recs: dict[str, list[dict]], why_in: list[str], commit: str) -> list[str]:
    """The row's void conditions on the per-seed records (task -> its records in seed order)."""
    why = list(why_in)
    for task in TASK_NAMES:
        rs_ = recs.get(task, [])
        if [r.get("seed") for r in rs_] != list(SEEDS):
            why.append(f"{task}: seeds are not {list(SEEDS)} (a missing task or seed)")
            continue
        for r in rs_:
            tag = f"{task} seed {r['seed']}"
            if r.get("run_commit") != commit:
                why.append(f"{tag}: another run commit")
            f = r.get("facts") or {}
            for b in BUILDS:
                res = (r.get("builds") or {}).get(b)
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
                    for ms, tol, said in ((SHAPE_METRICS, SHAPE_TOL, "1e-9"), (DE_METRICS, DE_TOL, "1e-4")):
                        dev = max(abs(row["skills"][n][m] - row["skills"][O][m]) for m in ms)
                        if not dev <= tol:
                            why.append(f"{tag}: {n}'s raw {', '.join(ms)} more than {said} off O's ({b})")
            for n in COUNTED:
                if not (f.get("same_coords") or {}).get(n):
                    why.append(f"{tag}: {n}'s coordinates differ from O's")
                if not (f.get("same_count") or {}).get(n):
                    why.append(f"{tag}: {n}'s count differs from O's")
                if not (f.get("same_sorted_rows") or {}).get(n):
                    why.append(f"{tag}: {n}'s sorted expression rows differ from O's")
    return why


def _per_part(recs, groups, a: str, b: str, what: str = "board") -> float:
    return float(np.mean(merged_D(recs, groups, a, b, what)))


def reading(recs: dict[str, list[dict]], why_in: list[str], commit: str) -> dict:
    why = check_run(recs, why_in, commit)
    if why:
        return {"void": True, "why": why, "commit": commit, "outcome": f"{ROW}: VOID"}
    m = merges()
    d = {"board": merged_D(recs, m["five"], ARM, O, "board"),
         "carried": merged_D(recs, m["five"], ARM, O, "carried"),
         "cell": merged_D(recs, m["five"], ARM, O, "group:cell_state"),
         "lo_pl": merged_D(recs, m["five"], ARM, PLACEBO, f"group:{G}"),
         "own": merged_D(recs, m["own_profile_d"], ARM, O, "board"),
         "classes": {k: merged_D(recs, m[k], ARM, O, f"group:{G}") for k in CLASSES}}
    dec = decide(d)
    per_truth = {truth: float(np.mean([np.mean(series(recs[n], ARM, "board") - series(recs[n], O, "board"))
                                       for n in names])) for truth, names in m["five"].items()}
    per_class = {k: {"board": _per_part(recs, m[k], ARM, O), G: _per_part(recs, m[k], ARM, O, f"group:{G}")}
                 for k in (*CLASSES, *OTHER_CLASSES)}
    apart = {k: {"board": _per_part(recs, m[k], ARM, O), G: _per_part(recs, m[k], ARM, O, f"group:{G}")}
             for k in ("half", "x3")}
    per_task = {}
    for t, rs_ in recs.items():
        x = series(rs_, ARM, "board") - series(rs_, O, "board")
        steps = [r["facts"].get("steps") or {} for r in rs_]
        per_task[t] = {
            "board": {"mean": float(x.mean()), "ahead": int((x > 0).sum()), "sd": float(x.std(ddof=1))},
            G: float(np.mean(series(rs_, ARM, f"group:{G}") - series(rs_, O, f"group:{G}"))),
            "placebo_board": float(np.mean(series(rs_, PLACEBO, "board") - series(rs_, O, "board"))),
            "o_neighbourhood_skill": float(np.mean(series(rs_, O, "skill:neighborhood_mmd"))),
            "board_minus_floor": {n: float(np.mean(series(rs_, n, "board") - series(rs_, harness.FLOOR, "board")))
                                  for n in RUN_NAMES},
            "step": {mode: {k: _median([s.get(mode, {}).get(k) for s in steps])
                            for k in ("contrast_earlier", "moved_share", "median_homology_error", "sd_f_o", "sd_f",
                                      "own_run_share_o", "own_run_share")}
                     for mode in ("run_sorted", "run_sorted_same")},
            "placebo_matches_swaps": all(bool((s.get("run_sorted_same") or {}).get("placebo_matches_swaps"))
                                         for s in steps),
            "files_identical_across_builds": all(r["facts"].get("files_identical_across_builds") for r in rs_)}
    adoption = keep_change({n: {s: (float(series(recs[n], O, "board")[i]), float(series(recs[n], ARM, "board")[i]))
                                for i, s in enumerate(SEEDS)} for n in ("e75", "e725")})
    res = {"void": False, "why": [], "commit": commit, "outcome": f"{ROW}: {'PASS' if dec['passes'] else 'STOP'}",
           "label": LABEL, "clauses": dec["clauses"], "idle": dec["idle"], "passes": dec["passes"],
           "clause_text": CLAUSES, "truths": m["five"],
           "merged": {"board": _stat(d["board"]), "carried": _stat(d["carried"]), "cell_state": _stat(d["cell"]),
                      "local_organisation_vs_placebo": _stat(d["lo_pl"]), "own_profile_board": _stat(d["own"]),
                      "classes_local_organisation": {k: _stat(v) for k, v in d["classes"].items()}},
           "per_truth": per_truth, "per_class": per_class, "half_and_x3": apart, "tasks": per_task,
           "reported_deciding_nothing": {
               "prefix_board": _stat(merged_D(recs, m["five"], ARM, O, "board", build="prefix")),
               "multiplicative_gate_board": _stat(merged_D(recs, m["five"], ARM, O, "mult")),
               "placebo_vs_O_board": _stat(merged_D(recs, m["five"], PLACEBO, O, "board")),
               "against_floor": {n: _stat(merged_D(recs, m["five"], n, harness.FLOOR, "board")) for n in RUN_NAMES},
               "adoption_rule_e75_e725": {"keep": adoption["keep"], "reason": adoption["reason"],
                                          "tasks": adoption["tasks"]}}}
    if dec["passes"]:
        gains = {k: v for k, v in per_truth.items() if v > 0}
        losses = {k: v for k, v in per_truth.items() if v < 0}
        res["question"] = {"label": LABEL, "command": CANDIDATE_COMMAND, "file": FILE, "per_truth": per_truth,
                           "gains": gains, "losses": losses, "per_class": per_class, "half_and_x3": apart,
                           "merged_board": float(d["board"].mean())}
    else:
        res["question"] = None
    return res


def text(res: dict) -> str:
    if res["void"]:
        return (f"{res['outcome']}\nRow {ROW} (run commit {res['commit']}): " + "; ".join(res["why"]) + "\n"
                + VOID_TEXT + "\n")
    m = res["merged"]
    L = [res["outcome"],
         f"Row {ROW} (run commit {res['commit']}); {LABEL}. E57b {ARM} against O {O}, build main, five truths merged "
         "at equal weight (mean, SE, seeds ahead of 10). E57b passes only if all hold:"]
    for k, v in res["clauses"].items():
        L.append(f"  {k} {'holds' if v else 'fails'}: {res['clause_text'][k]}")
    L.append(f"  board {m['board']['mean']:+.3f} (SE {m['board']['se']:.3f}, {m['board']['ahead']} ahead); carried "
             f"{m['carried']['mean']:+.3f}; cell_state {m['cell_state']['mean']:+.5f} (SE {m['cell_state']['se']:.5f})")
    L.append(f"  against E57c {PLACEBO}: local_organisation {m['local_organisation_vs_placebo']['mean']:+.4f} "
             f"(SE {m['local_organisation_vs_placebo']['se']:.4f}, {m['local_organisation_vs_placebo']['ahead']} "
             f"ahead); new2 and e75 board {m['own_profile_board']['mean']:+.3f}")
    L.append("  local_organisation per bracket-run class: " + "; ".join(
        f"{k} {v['mean']:+.4f}" for k, v in m["classes_local_organisation"].items()))
    if res["idle"]:
        L.append("  |E57b - O| is under 0.30 merged: the idle case reads STOP")
    L.append("Per truth (board): " + ", ".join(f"{k} {v:+.3f}" for k, v in res["per_truth"].items()))
    L.append("Per class (board / local_organisation): " + "; ".join(
        f"{k} {v['board']:+.3f} / {v[G]:+.4f}" for k, v in res["per_class"].items()))
    L.append("Half and x3 apart (the full dose; board / local_organisation): " + "; ".join(
        f"{k} {v['board']:+.3f} / {v[G]:+.4f}" for k, v in res["half_and_x3"].items()))
    rep = res["reported_deciding_nothing"]
    L.append(f"Reported, deciding nothing: prefix board {rep['prefix_board']['mean']:+.3f}; multiplicative gate board "
             f"{rep['multiplicative_gate_board']['mean']:+.3f}; placebo minus O board "
             f"{rep['placebo_vs_O_board']['mean']:+.3f}; minus the floor "
             + ", ".join(f"{n} {v['mean']:+.3f}" for n, v in rep["against_floor"].items())
             + f"; the adoption rule on e75 and e725: {rep['adoption_rule_e75_e725']['reason']}")
    for t, v in res["tasks"].items():
        st = v["step"]["run_sorted"]
        L.append(f"  {t}: board {v['board']['mean']:+.3f} ({v['board']['ahead']} ahead), {G} {v[G]:+.4f}, O's "
                 f"neighbourhood skill {v['o_neighbourhood_skill']:.3f}; step: contrast "
                 f"{_fmt(st['contrast_earlier'])}, moved {_fmt(st['moved_share'])}, own-run share "
                 f"{_fmt(st['own_run_share_o'])} -> "
                 f"{_fmt(st['own_run_share'])}; placebo swap count matched {v['placebo_matches_swaps']}")
    if res["passes"]:
        q = res["question"]
        L.append(f"PASS: a priced question to the captain with both readings, labelled '{q['label']}'; no file without "
                 "his ruling.")
        L.append("  gains per truth: " + (", ".join(f"{k} {v:+.3f}" for k, v in q["gains"].items()) or "none")
                 + "; losses per truth: " + (", ".join(f"{k} {v:+.3f}" for k, v in q["losses"].items()) or "none"))
        L.append(f"  file, only on his ruling: {q['command']}")
    else:
        L.append(STOP_TEXT)
    L.append("Nothing here shows which run E7.75 shares, the full dose beyond half and x3, R and the midline at the "
             "test pair, the live gate or E7.75's density.")
    return "\n".join(L) + "\n"


def load_runs(runs: Path) -> dict[str, list[dict]]:
    out = {}
    for name in TASK_NAMES:
        folder = Path(runs) / "run" / name
        recs = [json.loads(p.read_text()) for p in sorted(folder.glob("seed*.json"))] if folder.exists() else []
        out[name] = sorted(recs, key=lambda r: r["seed"])
    return out


def cmd_read(a) -> int:
    out = Path(a.out)
    files = [out / "reading.json", out / "reading.txt"]
    existing = [str(f) for f in files if f.exists()]
    if existing:
        print(f"refusing to overwrite: {', '.join(existing)}", file=sys.stderr)
        return 95
    why = own_code_problems()
    if why:
        sys.exit("refused: " + "; ".join(why))
    if not row_ok():
        sys.exit(f"refused: no single row {ROW} as the read needs it")
    runs = Path(a.runs)
    why = [*step0_passed(Path(a.step0), a.commit), *gate_problems(runs, a.commit), *commit_problems(a.commit)]
    if repo_state(REPO).unclean != ():
        why.append("the tree is not clean at the read")
    res = {**reading(load_runs(runs), why, a.commit), "reading_code": provenance(), "row": ROW}
    out.mkdir(parents=True, exist_ok=True)
    files[0].write_text(json.dumps(res, indent=1, default=float))
    files[1].write_text(text(res))
    print(f"wrote {files[1]}", flush=True)
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    p0 = sub.add_parser("step0")
    p0.add_argument("--name", choices=NAMES, required=True)
    p0.add_argument("--commit", required=True)
    p0.add_argument("--runs", default=str(SCRIPT_RUNS))
    p0.add_argument("--data", default=str(harness.DATA_DIR))
    r0 = sub.add_parser("step0-read")
    r0.add_argument("--commit", required=True)
    r0.add_argument("--runs", default=str(SCRIPT_RUNS))
    r0.add_argument("--out", default=str(STEP0_RESULTS))
    for name in ("gate", "run"):
        p = sub.add_parser(name)
        p.add_argument("--task", required=True, choices=TASK_NAMES)
        p.add_argument("--commit", required=True)
        p.add_argument("--runs", default=str(SCRIPT_RUNS))
        p.add_argument("--data", default=str(harness.DATA_DIR))
        if name == "run":
            p.add_argument("--step0", default=str(STEP0_RESULTS))
    r = sub.add_parser("read")
    r.add_argument("--commit", required=True)
    r.add_argument("--runs", default=str(SCRIPT_RUNS))
    r.add_argument("--out", default=str(RESULTS))
    r.add_argument("--step0", default=str(STEP0_RESULTS))
    a = ap.parse_args(argv)
    return {"step0": cmd_step0, "step0-read": cmd_step0_read, "gate": cmd_gate, "run": cmd_run,
            "read": cmd_read}[a.cmd](a)


if __name__ == "__main__":
    sys.exit(main())
