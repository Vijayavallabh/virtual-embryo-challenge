#!/usr/bin/env python
"""Try T42 (experiments/harness/TRIES.md): M1 with its Feng direction shrunk by reliability, a labelled read on T26's
route (experiments/t1_m1/gmom_combo.py), with the uniform James-Stein scale as its control.

N48 (target-free, committed) opened the decile form by its stop (a factor of deciles 4 to 10 below 0.9 in each of CM and
Endocardium) and closed the scalar James-Stein form (CM's scale 0.989 is not below 0.9). Its decile factors, deciles 4
to 10: CM 0.000, 0.203, 0.744, 0.932, 0.966, 0.982, 0.994; Endocardium 1, 1, 0.104, 0.422, 0.665, 0.793, 0.930;
Endothelium 1, 1, 0.111, 0.524, 0.697, 0.825, 0.932; Epicardium 1, 1, 0.333, 0.566, 0.667, 0.791, 0.925. A factor of 0
(max(msq - noise, 0)) sets the step of the first expressed decile to 0, CM's included. Those genes are not held:
like every free gene they are multiplied by the cell's depth constant and redrawn. CM holds most cells and is one of
the two groups the backtests read, so a pass may come from deleting the step of CM's weakly expressed genes (factor
0 in decile 4, 0.2 in decile 5), not from shrinking the small groups' noisy steps; this read cannot tell these
apart. The scalar
scale enters only as the uniform control of the decile form, never as an outcome (N48 closed it). The shrinkage uses
Feng's two halves only (no E10.5 information), so the paired difference between M1 with and without it is a fair
question to put to the two E10.5 cuts.

  run    One E10.5 cut, 10 seeds, on the harness's draws and scoring exactly as T26's run; refused unless the committed
         N22 passed. Arms per seed, all from the same draw and streams: b0 (the gate), m1_r and m1_r_gmom as T26 builds
         them; m1_r_eb and m1_r_gmom_eb, the same two builds with every group's direction replaced by
         feng_dirs.eb_shrink's output (try T18's rule and constants); m1_r_js and m1_r_gmom_js, every group's direction
         times that group's scalar n48_m1_direction_shrink.js_scale. Nothing else differs from the unshrunk twin.
         Refused unless experiments/harness/TRIES.md holds row T42 (a row is committed before its run).
  read   T42's reading of both cuts' files (rule below, fixed in the row before the run).

NOT an adoption result: Feng's E9.5 cells lie inside these pseudo-tasks' external windows and the cuts are ones earlier
T1 choices were made on. The result can only go to the team captain as a request.

The reading. Void unless both cuts hold the ten seeds at one commit on a clean tree, b0 is within 0.02 points per seed
of T26's committed scores on both cuts, the factors and scalar scales the run used equal N48's committed ones to 1e-12,
the statistics file's sha256 equals N48's, each file's task is its cut and its fraction 0.2, and no shrunk arm scores
exactly as its unshrunk twin (the same record) on every seed of both cuts (a build from unshrunk directions would pass
the factor gate, since the factors are recorded apart from the build). The test, the adoption rule's
(vec.adoption.keep_change, paired by seed, both cuts), on (a) m1_r_eb vs m1_r, (b) m1_r_gmom_eb vs m1_r_gmom, (c)
m1_r_js vs m1_r, (d) m1_r_gmom_js vs m1_r_gmom, the lead (e) m1_r_gmom_eb vs m1_r_gmom_js, and its reverse (f)
m1_r_gmom_js vs m1_r_gmom_eb, a veto on the decile. Verdict in order: (a), (b) and (e) pass and (f) does not: the
decile shrinkage is put forward. (a) and (b) pass otherwise: the decile shrinkage is not put forward; (f) passing: its
uniform control beats it; else (c) and (d) passing: its uniform control does as well; else: it is not shown to beat a
uniform step of the same projection. Not both (a) and (b): nothing is put forward; when (c) and (d) pass and none of
(a), (b), (e) does, only the closed scalar control reaches the rule. The scalar scale is never put forward (N48
closed it).

Usage: on the DGX campaign checkout, repo root, both cuts at one commit, committing nothing until both task files
exist (a commit between the cuts voids the read; the untracked task files leave the tree clean), and on the machine
T26 ran on, so that b0 reproduces. The environment variables are docs/COMPUTE.md's for repository Python on the DGX.
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \
      uv run --frozen --no-sync python -I -B experiments/t1_m1/eb_read.py run --out t42_m1_eb \
      --task "T1/E9.5->E10.5:chamber"
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \
      uv run --frozen --no-sync python -I -B experiments/t1_m1/eb_read.py run --out t42_m1_eb --task "T1/E9.5->E10.5"
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \
      uv run --frozen --no-sync python -I -B experiments/t1_m1/eb_read.py read --out t42_m1_eb
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
from vec.adoption import SEEDS, keep_change
from vec.boards import REPO_ROOT, load_boards
from vec.evaluate import Stage, score_board
from vec.io import read_stage
from vec.prep import load_stage

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
T26_TASKS = RESULTS / "t26_m1_gmom" / "tasks"
N48_FILE = RESULTS / "n48_m1_direction_shrink" / "n48.json"
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(REPO_ROOT / "experiments" / "t1_direction"))
sys.path.insert(0, str(REPO_ROOT / "experiments" / "external"))
import dirlib as dl  # noqa: E402
import feng_dirs  # noqa: E402
import filterio as fio  # noqa: E402
import gmom_combo as gc  # noqa: E402
import n12_hygiene_size as n12  # noqa: E402  (on the path through n48's import)
import n48_m1_direction_shrink as n48  # noqa: E402

gv, hygiene = gc.gv, gc.hygiene
CUT_FILES = gc.CUT_FILES
BASES = ("m1_r", "m1_r_gmom")
ARMS = ("b0", "m1_r", "m1_r_gmom", "m1_r_eb", "m1_r_gmom_eb", "m1_r_js", "m1_r_gmom_js")
COMPARISONS = {"a": ("m1_r", "m1_r_eb"), "b": ("m1_r_gmom", "m1_r_gmom_eb"), "c": ("m1_r", "m1_r_js"),
               "d": ("m1_r_gmom", "m1_r_gmom_js"), "e": ("m1_r_gmom_js", "m1_r_gmom_eb"),
               "f": ("m1_r_gmom_eb", "m1_r_gmom_js")}  # (base, change); (f) is (e) reversed, a veto on the decile
TWINS = (("m1_r", "m1_r_eb"), ("m1_r", "m1_r_js"), ("m1_r_gmom", "m1_r_gmom_eb"), ("m1_r_gmom", "m1_r_gmom_js"))
ROW = "T42"
FRACTION = 0.2  # Feng's step on an E9.5 -> E10.5 cut, feng_dirs.bracket_fraction(9.5, 10.5)
FACTOR_TOL = 1e-12  # the run's factors against N48's committed ones
REPRO_TOL = gc.REPRO_TOL  # 0.02 points per seed, b0 against T26's committed scores
DECILE = "the decile shrinkage is put forward"
DECILE_BEATEN = "the decile shrinkage is not put forward: its uniform control (closed by N48) beats it by the rule"
DECILE_AS_WELL = "the decile shrinkage is not put forward: its uniform control (closed by N48) does as well"
DECILE_NOT_SHOWN = ("the decile shrinkage is not put forward: it is not shown to beat a uniform step of the same "
                    "projection")
SCALAR_ONLY = "nothing is put forward: only the closed scalar control reaches the rule"
NOTHING = "nothing is put forward: closed by this reading"
LABEL_T42 = ("try T42, a labelled in-sample diagnostic, not an adoption result: M1 (with R1, R2, and with g_mom's "
             "factors) with its Feng direction shrunk by reliability (decile and scalar forms), on T26's route; Feng's "
             "E9.5 cells lie inside these pseudo-tasks' external windows and the cuts are ones earlier T1 choices were "
             "made on; the result can only go to the captain as a request")


# ---------------------------------------------------------------------------------------------------------------------
# The shrunk directions


def shrunk_directions(directions: Mapping[str, np.ndarray], halves: Mapping[str, Mapping[str, np.ndarray]]
                      ) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray], dict[str, dict]]:
    """Per group of `directions`: the decile form, feng_dirs.eb_shrink(d, d_a, d_b, mean_cp, LP_BINS)[0]; the scalar
    form, d times n48.js_scale(d, d_a, d_b)["scale"]; and the record of the factors used. NaN stays NaN. A group
    without halves, or whose scalar scale cannot be computed, is refused: no direction may go unshrunk silently."""
    eb, js, record = {}, {}, {}
    for g, d in directions.items():
        if g not in halves:
            raise ValueError(f"{g}: Feng's halves give no split for this direction: it cannot be shrunk")
        h = halves[g]
        d = np.asarray(d, dtype=np.float64)
        out, factors = feng_dirs.eb_shrink(d, h["d_a"], h["d_b"], h["mean_cp"], harness.LP_BINS)
        scale = n48.js_scale(d, h["d_a"], h["d_b"])["scale"]
        if scale is None:
            raise ValueError(f"{g}: no gene has a direction and both halves: no scalar scale")
        eb[g], js[g] = out, d * scale
        record[g] = {"eb_factors": [float(f) for f in factors], "js_scale": float(scale)}
    return eb, js, record


def feng_halves(genes: Sequence[str], groups) -> tuple[dict, str]:
    """feng_dirs.direction_halves on the clean statistics M1's direction is read from (manifest-checked), on `genes`,
    for `groups`, at M1's bracket; and the statistics file's sha256."""
    path = feng_dirs.clean_stats_path()
    entry = fio.verify_manifest_output(path)
    st = dl.load_stats(path, fields=dl.STAT_FIELDS)
    return feng_dirs.direction_halves(st, list(genes), feng_dirs.load_config(), feng_dirs.BRACKET,
                                      groups=set(groups)), entry["sha256"]


def all_arms(ctx: harness.BuildContext, directions: Mapping, eb: Mapping, js: Mapping, groups: Mapping,
             fraction: float, held_m1: Sequence[str], held_factors: Sequence[str]) -> dict[str, Stage]:
    """The six M1 arms of one seed: gmom_combo.m1_r_arms with M1's direction, with its decile-shrunk version and with
    its scalar-scaled version. The same context, groups, fraction, held genes and streams for each."""
    out = dict(gc.m1_r_arms(ctx, directions, groups, fraction, held_m1, held_factors))
    for suffix, dirs in (("_eb", eb), ("_js", js)):
        out.update({k + suffix: v for k, v in gc.m1_r_arms(ctx, dirs, groups, fraction, held_m1, held_factors).items()})
    return out


def score_arms(board, truth: Stage, ref_stage: Stage, seed: int, anchors, arms: Mapping[str, Stage]) -> dict:
    """g_mom_verify.score_arms's record (score, group skills, the metric panel) with each metric's skill added, from one
    score_board call per arm (T26's files hold no per-metric skill, which this reading reports)."""
    out = {}
    for name, pred in arms.items():
        r = score_board(board, pred, truth, ref_stage, seed=seed, anchors=anchors, build=harness.BUILD,
                        protocol=harness.PROTOCOL)
        out[name] = {"score": float(r.score), "groups": {g: float(v) for g, v in r.group_skill.items()},
                     "metrics": {m: float(r.metrics[m]) for m in gv.PANEL_METRICS if m in r.metrics},
                     "skills": {m: float(v) for m, v in r.skills.items()}}
    return out


def cmd_run(a) -> None:
    if ROW not in n12.tries_rows():
        sys.exit(f"T42 is refused: {ROW} has no row in experiments/harness/TRIES.md (a row comes before its run)")
    n22 = RESULTS / a.n22 / "n22.json"
    if not gc.n22_passed(n22):
        sys.exit(f"T42 is refused: {gv._rel(n22)} does not hold a passing N22 (T26's route needs N22 passed)")
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
    halves, stats_sha = feng_halves(list(pool.var_names), directions)
    if stats_sha != ext["sha256"]:
        sys.exit("the statistics the halves came from are not the ones M1's direction came from")
    eb, js, shrinkage = shrunk_directions(directions, halves)
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
        preds = {"b0": b0(ctx), **all_arms(ctx, directions, eb, js, groups, fraction, held_m1, preregistered)}
        arms = score_arms(board, target, ref_stage, s, f.anchors, preds)
        rows.append({"seed": s, "floor": float(f.score), "arms": arms, "seconds": round(time.time() - t0, 1)})
        print(f"{t.key} seed {s}: " + ", ".join(f"{k} {v['score']:.2f}" for k, v in arms.items()), flush=True)
    res = {"label": LABEL_T42, "task": t.key, "fraction": fraction, "external": ext, "seeds": list(seeds),
           "shrinkage": shrinkage, "feng_stats_sha256": stats_sha, "rows": rows,
           "stamp": gv.stamp(started, [data / t.target, data / t.reference, *[data / f for f in t.inputs.values()],
                                       t1_scope.HELD_FLAT_GENES, ts.R1_CONFIG, REPO_ROOT / union["file"],
                                       feng_dirs.clean_stats_path(), n22])}
    (out / "tasks").mkdir(parents=True, exist_ok=True)
    (out / "tasks" / name).write_text(json.dumps(res, indent=1, default=float) + "\n")


# ---------------------------------------------------------------------------------------------------------------------
# The reading


def validity(results: Mapping[str, Mapping], t26: Mapping[str, Mapping], n48_res: Mapping) -> dict:
    """The void conditions. results and t26: cut -> task file (a missing cut is absent); n48_res: N48's n48.json."""
    n48_groups, n48_sha = n48_res["groups"], n48_res["inputs"]["feng_stats_sha256"]
    reasons = []
    missing = [k for k in CUT_FILES if k not in results]
    if missing:
        reasons.append(f"missing cuts {missing}")
    for k, r in results.items():
        seeds = [row["seed"] for row in r["rows"]]
        if sorted(seeds) != list(SEEDS) or list(r.get("seeds", seeds)) != list(SEEDS):
            reasons.append(f"{k}: seeds {seeds}, not {list(SEEDS)}")
        absent = sorted({a for row in r["rows"] for a in ARMS if a not in row["arms"]})
        if absent:
            reasons.append(f"{k}: arms {absent} missing")
        if r.get("task") != k:
            reasons.append(f"{k}: the file holds task {r.get('task')!r}")
        if not isinstance(r.get("fraction"), int | float) or abs(r["fraction"] - FRACTION) > FACTOR_TOL:
            reasons.append(f"{k}: fraction {r.get('fraction')!r}, not {FRACTION:g}")
        if r.get("feng_stats_sha256") != n48_sha:
            reasons.append(f"{k}: the statistics' sha256 is not N48's")
    commits = {r["stamp"]["commit"] for r in results.values()}
    unclean = sorted(k for k, r in results.items() if r["stamp"]["tree_clean"] is not True)
    if unclean:
        reasons.append(f"tree not clean on {unclean}")
    if len(commits) > 1:
        reasons.append(f"the cuts ran at {len(commits)} commits")
    repro, worst_factor = {}, 0.0
    for k, r in results.items():
        repro[k] = reproduction(r["rows"], t26.get(k, {}).get("rows", []))
        if not repro[k]["b0_holds"]:
            reasons.append(f"{k}: b0 does not reproduce T26's committed scores within {REPRO_TOL:g} per seed")
        got = r.get("shrinkage", {})
        if set(got) != set(n48_groups):
            reasons.append(f"{k}: shrunk groups {sorted(got)}, N48's {sorted(n48_groups)}")
        for g in set(got) & set(n48_groups):
            want_eb, want_js = n48_groups[g]["eb_factors"], n48_groups[g]["js"]["scale"]
            if len(got[g]["eb_factors"]) != len(want_eb) or want_js is None:
                worst_factor = float("inf")
                continue
            diffs = [abs(x - y) for x, y in zip(got[g]["eb_factors"], want_eb, strict=True)]
            worst_factor = max(worst_factor, *diffs, abs(got[g]["js_scale"] - want_js))
    if results and worst_factor > FACTOR_TOL:
        reasons.append(f"the factors differ from N48's by {worst_factor:.3g} (> {FACTOR_TOL:g})")
    complete = all(a in row["arms"] for r in results.values() for row in r["rows"] for a in ARMS)
    if results and complete:
        same = [f"{c} = {b}" for b, c in TWINS
                if all(row["arms"][c] == row["arms"][b] for r in results.values() for row in r["rows"])]
        if same:
            reasons.append(f"shrunk arms score exactly as their twins on every seed (built unshrunk?): {same}")
    return {"valid": not reasons, "reasons": reasons, "reproduction": repro, "largest_factor_difference": worst_factor,
            "commits": sorted(c for c in commits if c)}


def reproduction(rows: Sequence[Mapping], t26_rows: Sequence[Mapping]) -> dict:
    """Per seed against T26's committed scores: b0 (the gate) and, reported only, m1_r and m1_r_gmom."""
    by_seed = {r["seed"]: r for r in t26_rows}
    missing = [r["seed"] for r in rows if r["seed"] not in by_seed]
    worst = {a: 0.0 for a in ("b0", *BASES)}
    for r in rows:
        if r["seed"] in by_seed:
            for arm in worst:
                if arm in r["arms"] and arm in by_seed[r["seed"]]["arms"]:
                    old = by_seed[r["seed"]]["arms"][arm]["score"]
                    worst[arm] = max(worst[arm], abs(r["arms"][arm]["score"] - old))
    return {"largest_difference": worst, "missing_seeds": missing,
            "b0_holds": bool(t26_rows) and not missing and worst["b0"] <= REPRO_TOL}


def _mean_diff(rows: Sequence[Mapping], base: str, change: str, field: str) -> dict:
    """The mean paired difference, change minus base, of every entry of `field` that every row holds for both arms."""
    keys = sorted(set.intersection(*[set(row["arms"][a].get(field, {})) for row in rows for a in (base, change)]))
    return {m: float(np.mean([row["arms"][change][field][m] - row["arms"][base][field][m] for row in rows]))
            for m in keys}


def compare(results: Mapping[str, Mapping], base: str, change: str) -> dict:
    """The adoption rule on change against base, paired by seed on both cuts, and the reported paired differences."""
    paired = {k: {row["seed"]: (row["arms"][base]["score"], row["arms"][change]["score"]) for row in r["rows"]}
              for k, r in results.items()}
    rule = keep_change(paired)
    cuts = {}
    for k, r in results.items():
        rows = r["rows"]
        d = np.array([row["arms"][change]["score"] - row["arms"][base]["score"] for row in rows])

        cuts[k] = {"mean": float(d.mean()), "ahead": int((d > 0).sum()), "of": len(d),
                   **{f: _mean_diff(rows, base, change, f) for f in ("skills", "groups", "metrics")}}
    return {"base": base, "change": change, "passes": bool(rule["keep"]), "rule": rule, "cuts": cuts}


def verdict(p: Mapping[str, bool]) -> str:
    """The decile goes forward only when (a), (b) and (e) pass and (f) does not; the scalar is never put forward, and
    'only the scalar control reaches the rule' is said only when (c) and (d) pass and none of (a), (b), (e) does."""
    if p["a"] and p["b"]:
        if p["e"] and not p["f"]:
            return DECILE
        if p["f"]:
            return DECILE_BEATEN
        return DECILE_AS_WELL if p["c"] and p["d"] else DECILE_NOT_SHOWN
    return SCALAR_ONLY if p["c"] and p["d"] and not (p["a"] or p["b"] or p["e"]) else NOTHING


def t42_reading(results: Mapping[str, Mapping], t26: Mapping[str, Mapping], n48_res: Mapping) -> dict:
    v = validity(results, t26, n48_res)
    out = {"valid": v["valid"], "validity": v}
    if not v["valid"]:
        out["outcome"] = "VOID"
        return out
    comps = {c: compare(results, *pair) for c, pair in COMPARISONS.items()}
    out["comparisons"] = comps
    out["outcome"] = verdict({c: comps[c]["passes"] for c in comps})
    out["mean_scores"] = {k: {a: float(np.mean([row["arms"][a]["score"] for row in r["rows"]])) for a in ARMS}
                          for k, r in results.items()}
    return out


def t42_text(res: Mapping) -> str:
    r = res["reading"]
    lines = [res["label"], f"Verdict: {r['outcome']}"]
    v = r["validity"]
    if not r["valid"]:
        lines.append("Void: " + "; ".join(v["reasons"]))
    for k, g in v["reproduction"].items():
        w = g["largest_difference"]
        lines.append(f"{k}: largest per-seed difference from T26's committed scores: b0 {w['b0']:.4f} (gate), m1_r "
                     f"{w['m1_r']:.4f}, m1_r_gmom {w['m1_r_gmom']:.4f} (reported)")
    lines.append(f"Largest difference of the factors from N48's: {v['largest_factor_difference']:.3g}")
    for c, comp in r.get("comparisons", {}).items():
        lines.append(f"({c}) {comp['change']} against {comp['base']}: {'PASSES' if comp['passes'] else 'does not pass'}"
                     f" ({comp['rule']['reason']})")
        for k, cut in comp["cuts"].items():
            sk = ", ".join(f"{m} {x:+.4f}" for m, x in cut["skills"].items())
            lines.append(f"    {k}: {cut['mean']:+.2f} points, ahead {cut['ahead']} of {cut['of']}; skills: {sk}")
    for k, ms in r.get("mean_scores", {}).items():
        lines.append(f"mean scores {k}: " + ", ".join(f"{a} {x:.2f}" for a, x in ms.items()))
    lines.append("Limit: in-sample labelled diagnostic; never adoption evidence; a request to the captain at most")
    lines.append(f"Commit {res['stamp']['commit']}, tree clean {res['stamp']['tree_clean']}, UTC {res['stamp']['utc']}")
    return "\n".join(lines) + "\n"


def read(out: Path, t26_tasks: Path = T26_TASKS, n48_file: Path = N48_FILE) -> dict:
    """Read the cuts present in out/tasks against T26's committed task files and N48's factors; write t42.{json,txt}."""
    gv.refuse_existing(out, ("t42.json", "t42.txt"))
    started = gv.tree_state()
    present = {k: out / "tasks" / f for k, f in CUT_FILES.items() if (out / "tasks" / f).exists()}
    results = {k: json.loads(p.read_text()) for k, p in present.items()}
    t26 = {k: json.loads((t26_tasks / f).read_text()) for k, f in CUT_FILES.items() if (t26_tasks / f).exists()}
    n48_res = json.loads(Path(n48_file).read_text())
    res = {"label": LABEL_T42, "reading": t42_reading(results, t26, n48_res),
           "stamps": {k: v["stamp"] for k, v in results.items()}}
    res["stamp"] = gv.stamp(started, [*present.values(), *[t26_tasks / f for f in CUT_FILES.values()
                                                          if (t26_tasks / f).exists()], Path(n48_file)])
    gv.write_outputs(out, "t42", res, t42_text(res))
    return res


def cmd_read(a) -> None:
    res = read(RESULTS / a.out)
    if res["reading"]["outcome"] == "VOID":
        sys.exit("T42's reading is void and goes to the captain")


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("run", help="one cut of try T42")
    p.add_argument("--task", required=True, choices=list(CUT_FILES))
    p.add_argument("--out", required=True)
    p.add_argument("--n22", default="n22_additivity")
    p.add_argument("--seeds", default=",".join(map(str, SEEDS)))
    p.add_argument("--data", default=str(harness.DATA_DIR))
    p = sub.add_parser("read", help="try T42's reading")
    p.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    {"run": cmd_run, "read": cmd_read}[a.cmd](a)


if __name__ == "__main__":
    main()
