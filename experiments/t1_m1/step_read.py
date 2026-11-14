#!/usr/bin/env python
"""Try T44 (experiments/harness/TRIES.md): is M1's step too small or too large on the E10.5 cuts? A labelled read on
try T42's route (experiments/t1_m1/eb_read.py), with M1's Feng step made larger.

M1 applies Feng's E9.5 -> E14.5 change in proportion to days: a fifth of the bracket on a one-day practice cut, two
fifths at the test. Try T42 measured the smaller side on these cuts (its scalar control, a step shrunk to 0.86 to 0.99
of its size, lost about 0.13 on every seed); the larger side has never been read. This read puts two larger steps to
the same cuts on the same route.

  run    One E10.5 cut, 10 seeds, on the harness's draws and scoring exactly as T42's run; refused unless the committed
         N22 passed, and unless experiments/harness/TRIES.md holds row T44 (a row is committed before its run). Arms per
         seed, all from one cell draw: b0 (the platform gate) and m1_r_gmom (M1 with R1, R2 and T11's factors: the
         entry's form) as T42 builds them; m1_r_gmom_x15 and m1_r_gmom_x2, the same build (gmom_combo.m1_r_arms) with
         the direction's fraction multiplied by 1.5 and by 2 (0.3 and 0.4 of Feng's bracket on these cuts) and nothing
         else changed: the same held genes, groups, composition, count, streams and g_mom factors. The two multipliers
         are the only numbers: 2 is the test's own fraction, 1.5 the point between. The task file holds every arm's
         fraction.
  read   T44's reading of both cuts' files (rule below, fixed in the row before the run).

What the arms share and do not: the scaled arms reach external_shift with the same directions, groups, held genes,
cells and stream ([seed, 2]) as the base; only `fraction` differs. Downstream of it, as in T42's twins, each covered
cell's keep-depth constant follows the step, and with it which counts are drawn up (Poisson) or down (Binomial) and so
every redraw of the covered cells on that stream; g_mom's factors (stream [seed, 5]) then act on different counts. A
per-seed difference therefore includes resampling noise the arms do not share.

NOT an adoption result: Feng's E9.5 cells lie inside these pseudo-tasks' external windows, the cuts are ones earlier T1
choices were made on, and the truth of a one-day cut differs from its reference in the kind of sample as well as in
age. The result can only go to the team captain.

The reading. Void unless both cuts hold the ten seeds at one commit on a clean tree; b0 and m1_r_gmom equal, seed for
seed, their committed scores in experiments/t1_m1/results/t42_m1_eb/tasks/ to 1e-9; each file's task is its cut, its
base fraction 0.2 and each arm's fraction its multiplier times that; and neither scaled arm scores exactly as
m1_r_gmom on every seed of both cuts. The test is the adoption rule's (vec.adoption.keep_change, paired by seed, both
cuts) for each scaled arm against m1_r_gmom. Verdict: "a larger step scores better on these cuts" if either arm passes
(naming which); "a larger step loses" if both arms' mean differences are below zero on both cuts with at most 2 of 10
seeds ahead on each; otherwise "the step's size is not separated on these cuts". The read puts nothing forward by
itself. Reported: per cut every arm's mean board score, each arm's mean paired difference with seeds ahead and per-seed
range, every metric's mean skill difference, and, read from T42's committed files, the scalar control's difference
(m1_r_gmom_js against m1_r_gmom), the response on the other side of the step.

Usage: on the DGX campaign checkout, repo root, both cuts at one commit, committing nothing until both task files
exist, and on the machine T42 ran on, so that b0 and m1_r_gmom reproduce. The environment variables are
docs/COMPUTE.md's for repository Python on the DGX.
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \
      uv run --frozen --no-sync python -I -B experiments/t1_m1/step_read.py run --out t44_m1_step \
      --task "T1/E9.5->E10.5:chamber"
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \
      uv run --frozen --no-sync python -I -B experiments/t1_m1/step_read.py run --out t44_m1_step \
      --task "T1/E9.5->E10.5"
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \
      uv run --frozen --no-sync python -I -B experiments/t1_m1/step_read.py read --out t44_m1_step
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import time
from collections.abc import Mapping, Sequence
from pathlib import Path

import numpy as np

from vec import harness, t1_scope
from vec import t1_state_transport as ts
from vec.adoption import SEEDS
from vec.boards import REPO_ROOT, load_boards
from vec.evaluate import Stage
from vec.io import read_stage
from vec.prep import load_stage

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
T42_TASKS = RESULTS / "t42_m1_eb" / "tasks"
sys.path.insert(0, str(HERE))
import eb_read as er  # noqa: E402
import feng_dirs  # noqa: E402
import gmom_combo as gc  # noqa: E402

gv, hygiene, n12 = gc.gv, gc.hygiene, er.n12
CUT_FILES = gc.CUT_FILES
BASE = "m1_r_gmom"
MULTIPLIERS = {"m1_r_gmom_x15": 1.5, "m1_r_gmom_x2": 2.0}  # the only numbers: the point between, the test's fraction
ARMS = ("b0", BASE, *MULTIPLIERS)
ROW = "T44"
FRACTION = er.FRACTION  # 0.2, Feng's step on an E9.5 -> E10.5 cut
FRACTION_TOL = 1e-12
REPRO_TOL = 1e-9  # b0 and m1_r_gmom against T42's committed scores: same machine and code path
REPRO_ARMS = ("b0", BASE)
SCALAR_CONTROL = "m1_r_gmom_js"  # T42's uniform smaller step on the entry's form, reported from its committed files
LOSES_MAX_AHEAD = 2
BETTER = "a larger step scores better on these cuts"
LOSES = "a larger step loses"
NOT_SEPARATED = "the step's size is not separated on these cuts"
LABEL_T44 = ("try T44, a labelled in-sample diagnostic, never an adoption result: M1 (with R1, R2 and g_mom's factors, "
             "the entry's form) with its Feng step multiplied by 1.5 and by 2, on T42's route; Feng's E9.5 cells lie "
             "inside these pseudo-tasks' external windows, the cuts are ones earlier T1 choices were made on, and the "
             "truth of a one-day cut differs from its reference in the kind of sample as well as in age; the result "
             "can only go to the team captain")


# ---------------------------------------------------------------------------------------------------------------------
# The run


def arm_fractions(fraction: float, multipliers: Mapping[str, float] = MULTIPLIERS) -> dict[str, float]:
    return {BASE: float(fraction), **{a: float(fraction) * m for a, m in multipliers.items()}}


def step_arms(ctx: harness.BuildContext, directions: Mapping, groups: Mapping, fraction: float,
              held_m1: Sequence[str], held_factors: Sequence[str],
              multipliers: Mapping[str, float] = MULTIPLIERS) -> dict[str, Stage]:
    """m1_r_gmom as T42 builds it (gmom_combo.m1_r_arms at `fraction`) and, per multiplier, the same call with the
    fraction multiplied and every other argument the same object."""
    fr = arm_fractions(fraction, multipliers)
    return {a: gc.m1_r_arms(ctx, directions, groups, f, held_m1, held_factors)[BASE] for a, f in fr.items()}


def cmd_run(a) -> None:
    if ROW not in n12.tries_rows():
        sys.exit(f"T44 is refused: {ROW} has no row in experiments/harness/TRIES.md (a row comes before its run)")
    n22 = RESULTS / a.n22 / "n22.json"
    if not gc.n22_passed(n22):
        sys.exit(f"T44 is refused: {gv._rel(n22)} does not hold a passing N22 (T26's route needs N22 passed)")
    out = RESULTS / a.out
    name = re.sub(r"[^A-Za-z0-9.]+", "_", a.task).strip("_") + ".json"
    gv.refuse_existing(out / "tasks", (name,))
    started = gv.tree_state()
    seeds = gv.parse_seeds(a.seeds)
    data = Path(a.data)
    t = {x.key: x for x in harness.load_tasks()["T1:test"]}[a.task]
    boards = load_boards()
    genes = boards[t.panel_board].genes
    target, _ = (harness.filter_target(t, data, genes) if t.target_filter is not None else
                 (Stage.from_anndata(read_stage(data / t.target, genes), genes), None))
    ref_stage = Stage.from_anndata(read_stage(data / t.reference, genes), genes)
    board = harness.task_board(t, boards, target.n, ref_stage.n)
    inputs = {r: load_stage(data / f, lineage=False) for r, f in t.inputs.items()}
    pool = inputs[harness.base_role(t)]
    directions, ext = feng_dirs.feng_directions(list(pool.var_names))
    groups = feng_dirs.label_groups(pool.obs["celltype"].astype(str).unique())
    fraction = feng_dirs.bracket_fraction(harness.STAGE_DAYS[t.reference], harness.STAGE_DAYS[t.target])
    preregistered = t1_scope.held_flat_genes()
    union = hygiene.load_union()
    held_m1 = gc.held_feng(union, preregistered)
    view, cboard = harness.TaskView.of(t), harness._candidate_board(board)
    b0 = harness._b0_t1(harness.B0_DAMPING)
    rows = []
    for s in seeds:
        t0 = time.time()
        f = gv.floor_anchors(t, board, target, ref_stage, s)
        ctx = harness.BuildContext(task=view, board=cboard, inputs=inputs, seed=harness.BUILD_SEED_OFFSET + s,
                                   derived={"t1_held_flat": preregistered})
        preds = {"b0": b0(ctx), **step_arms(ctx, directions, groups, fraction, held_m1, preregistered)}
        arms = er.score_arms(board, target, ref_stage, s, f.anchors, preds)
        rows.append({"seed": s, "floor": float(f.score), "arms": arms, "seconds": round(time.time() - t0, 1)})
        print(f"{t.key} seed {s}: " + ", ".join(f"{k} {v['score']:.2f}" for k, v in arms.items()), flush=True)
    res = {"label": LABEL_T44, "task": t.key, "fraction": fraction, "fractions": arm_fractions(fraction),
           "multipliers": dict(MULTIPLIERS), "external": ext, "seeds": list(seeds), "rows": rows,
           "stamp": gv.stamp(started, [data / t.target, data / t.reference, *[data / f for f in t.inputs.values()],
                                       t1_scope.HELD_FLAT_GENES, ts.R1_CONFIG, REPO_ROOT / union["file"],
                                       feng_dirs.clean_stats_path(), n22])}
    (out / "tasks").mkdir(parents=True, exist_ok=True)
    (out / "tasks" / name).write_text(json.dumps(res, indent=1, default=float) + "\n")


# ---------------------------------------------------------------------------------------------------------------------
# The reading


def reproduction(rows: Sequence[Mapping], t42_rows: Sequence[Mapping]) -> dict:
    """Per seed, b0 and m1_r_gmom against T42's committed scores; holds when every seed is there and each is within
    REPRO_TOL."""
    by_seed = {r["seed"]: r for r in t42_rows}
    missing = [r["seed"] for r in rows if r["seed"] not in by_seed]
    worst = {a: 0.0 for a in REPRO_ARMS}
    absent = []
    for r in rows:
        if r["seed"] not in by_seed:
            continue
        for arm in REPRO_ARMS:
            if arm not in r["arms"] or arm not in by_seed[r["seed"]]["arms"]:
                absent.append(f"{r['seed']}:{arm}")
                continue
            worst[arm] = max(worst[arm], abs(r["arms"][arm]["score"] - by_seed[r["seed"]]["arms"][arm]["score"]))
    holds = bool(t42_rows) and not missing and not absent and all(w <= REPRO_TOL for w in worst.values())
    return {"largest_difference": worst, "missing_seeds": missing, "absent": absent, "holds": holds}


def validity(results: Mapping[str, Mapping], t42: Mapping[str, Mapping]) -> dict:
    """The void conditions. results and t42: cut -> task file (a missing cut is absent)."""
    reasons = []
    missing = [k for k in CUT_FILES if k not in results]
    if missing:
        reasons.append(f"missing cuts {missing}")
    want_fr = arm_fractions(FRACTION)
    for k, r in results.items():
        seeds = [row["seed"] for row in r["rows"]]
        if sorted(seeds) != list(SEEDS) or list(r.get("seeds", seeds)) != list(SEEDS):
            reasons.append(f"{k}: seeds {seeds}, not {list(SEEDS)}")
        absent = sorted({a for row in r["rows"] for a in ARMS if a not in row["arms"]})
        if absent:
            reasons.append(f"{k}: arms {absent} missing")
        if r.get("task") != k:
            reasons.append(f"{k}: the file holds task {r.get('task')!r}")
        if not isinstance(r.get("fraction"), int | float) or abs(r["fraction"] - FRACTION) > FRACTION_TOL:
            reasons.append(f"{k}: fraction {r.get('fraction')!r}, not {FRACTION:g}")
        got_fr = r.get("fractions") or {}
        bad = [a for a, f in want_fr.items()
               if not isinstance(got_fr.get(a), int | float) or abs(got_fr[a] - f) > FRACTION_TOL]
        if bad:
            reasons.append(f"{k}: arm fractions {got_fr}, not {want_fr} (wrong fraction on {bad})")
    commits = {r["stamp"]["commit"] for r in results.values()}
    unclean = sorted(k for k, r in results.items() if r["stamp"]["tree_clean"] is not True)
    if unclean:
        reasons.append(f"tree not clean on {unclean}")
    if len(commits) > 1:
        reasons.append(f"the cuts ran at {len(commits)} commits")
    repro = {}
    for k, r in results.items():
        repro[k] = reproduction(r["rows"], t42.get(k, {}).get("rows", []))
        if not repro[k]["holds"]:
            reasons.append(f"{k}: b0 and m1_r_gmom do not reproduce T42's committed scores within {REPRO_TOL:g} per "
                           "seed")
    complete = all(a in row["arms"] for r in results.values() for row in r["rows"] for a in ARMS)
    if results and complete:
        same = [a for a in MULTIPLIERS
                if all(row["arms"][a] == row["arms"][BASE] for r in results.values() for row in r["rows"])]
        if same:
            reasons.append(f"scaled arms score exactly as {BASE} on every seed (built at the base fraction?): {same}")
    return {"valid": not reasons, "reasons": reasons, "reproduction": repro,
            "commits": sorted(c for c in commits if c)}


def compare(results: Mapping[str, Mapping], base: str, change: str) -> dict:
    """eb_read.compare (the adoption rule, paired by seed, both cuts; the mean paired differences) with each cut's
    per-seed range of the difference added."""
    out = er.compare(results, base, change)
    for k, r in results.items():
        d = [row["arms"][change]["score"] - row["arms"][base]["score"] for row in r["rows"]]
        out["cuts"][k]["range"] = [float(min(d)), float(max(d))]
    return out


def verdict(comps: Mapping[str, Mapping]) -> str:
    passing = [a for a in MULTIPLIERS if comps[a]["passes"]]
    if passing:
        return f"{BETTER}: {', '.join(passing)}"
    if all(c["mean"] < 0 and c["ahead"] <= LOSES_MAX_AHEAD for a in MULTIPLIERS for c in comps[a]["cuts"].values()):
        return LOSES
    return NOT_SEPARATED


def t44_reading(results: Mapping[str, Mapping], t42: Mapping[str, Mapping]) -> dict:
    v = validity(results, t42)
    out = {"valid": v["valid"], "validity": v}
    if not v["valid"]:
        out["outcome"] = "VOID"
        return out
    comps = {a: compare(results, BASE, a) for a in MULTIPLIERS}
    out["comparisons"] = comps
    out["outcome"] = verdict(comps)
    out["mean_scores"] = {k: {a: float(np.mean([row["arms"][a]["score"] for row in r["rows"]])) for a in ARMS}
                          for k, r in results.items()}
    if all(k in t42 for k in CUT_FILES) and all(SCALAR_CONTROL in row["arms"] and BASE in row["arms"]
                                                 for r in t42.values() for row in r["rows"]):
        out["scalar_control_t42"] = compare({k: t42[k] for k in CUT_FILES}, BASE, SCALAR_CONTROL)
    return out


def t44_text(res: Mapping) -> str:
    r = res["reading"]
    lines = [res["label"], f"Verdict: {r['outcome']}"]
    v = r["validity"]
    if not r["valid"]:
        lines.append("Void: " + "; ".join(v["reasons"]))
    for k, g in v["reproduction"].items():
        w = g["largest_difference"]
        lines.append(f"{k}: largest per-seed difference from T42's committed scores: b0 {w['b0']:.3g}, {BASE} "
                     f"{w[BASE]:.3g} (gate {REPRO_TOL:g})")
    for k, ms in r.get("mean_scores", {}).items():
        lines.append(f"mean scores {k}: " + ", ".join(f"{a} {x:.2f}" for a, x in ms.items()))

    def block(name: str, comp: Mapping) -> None:
        lines.append(f"{name}: {comp['change']} against {comp['base']}: "
                     f"{'PASSES' if comp['passes'] else 'does not pass'} ({comp['rule']['reason']})")
        for k, cut in comp["cuts"].items():
            sk = ", ".join(f"{m} {x:+.4f}" for m, x in cut["skills"].items())
            lines.append(f"    {k}: {cut['mean']:+.2f} points, ahead {cut['ahead']} of {cut['of']}, per seed "
                         f"{cut['range'][0]:+.2f} to {cut['range'][1]:+.2f}; skills: {sk}")

    for a, comp in r.get("comparisons", {}).items():
        block(f"x{MULTIPLIERS[a]:g}", comp)
    if "scalar_control_t42" in r:
        block("T42's scalar control, the smaller side (committed files, reported)", r["scalar_control_t42"])
    lines.append("Limit: in-sample labelled diagnostic; never adoption evidence; puts nothing forward by itself; it "
                 "can only go to the team captain")
    lines.append(f"Commit {res['stamp']['commit']}, tree clean {res['stamp']['tree_clean']}, UTC {res['stamp']['utc']}")
    return "\n".join(lines) + "\n"


def read(out: Path, t42_tasks: Path = T42_TASKS) -> dict:
    """Read the cuts present in out/tasks against T42's committed task files; write t44.{json,txt}."""
    gv.refuse_existing(out, ("t44.json", "t44.txt"))
    started = gv.tree_state()
    present = {k: out / "tasks" / f for k, f in CUT_FILES.items() if (out / "tasks" / f).exists()}
    results = {k: json.loads(p.read_text()) for k, p in present.items()}
    t42_paths = {k: t42_tasks / f for k, f in CUT_FILES.items() if (t42_tasks / f).exists()}
    t42 = {k: json.loads(p.read_text()) for k, p in t42_paths.items()}
    res = {"label": LABEL_T44, "reading": t44_reading(results, t42),
           "stamps": {k: v["stamp"] for k, v in results.items()}}
    res["stamp"] = gv.stamp(started, [*present.values(), *t42_paths.values()])
    gv.write_outputs(out, "t44", res, t44_text(res))
    return res


def cmd_read(a) -> None:
    res = read(RESULTS / a.out)
    if res["reading"]["outcome"] == "VOID":
        sys.exit("T44's reading is void and goes to the captain")


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("run", help="one cut of try T44")
    p.add_argument("--task", required=True, choices=list(CUT_FILES))
    p.add_argument("--out", required=True)
    p.add_argument("--n22", default="n22_additivity")
    p.add_argument("--seeds", default=",".join(map(str, SEEDS)))
    p.add_argument("--data", default=str(harness.DATA_DIR))
    p = sub.add_parser("read", help="try T44's reading")
    p.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    {"run": cmd_run, "read": cmd_read}[a.cmd](a)


if __name__ == "__main__":
    main()
