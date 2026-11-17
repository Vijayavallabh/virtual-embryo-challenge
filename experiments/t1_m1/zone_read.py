#!/usr/bin/env python
"""Try T53 (experiments/harness/TRIES.md): M1's Feng direction with the heart zones held fixed between the bracket's two
ends. A labelled in-sample read on the routes of try T44 (experiments/t1_m1/step_read.py, the whole cut) and try T46
(experiments/t1_m1/lib_step_read.py, the chamber cut pooled and library by library).

M1 moves each covered cell along d, per cell group the log ratio of mean CP10k between Feng's E9.5 and E14.5 hearts
(feng_dirs.directions_from_stats). dirlib.feng_sel pools a stage's cells of the group's Feng labels over every heart
zone, and Feng's capture is not the same at the two ends (the ledger: mostly LV-zone cells at E9.5, zone-balanced at
E14.5), so a pooled direction holds, beside the change over five days, a change in which zones the cells come from.
This read puts two zone-fixed forms of the direction to the same cuts.

Directions (per group on the entry's gene map, the config's pseudocount and min_cells_ext; the statistics' keys are
stage|zone|cell type|half and both halves are pooled, as feng_sel pools them; a key without exactly four parts is
refused):
  within-zone  the log ratio (dirlib.lfc_ratio) of the group's mean CP10k over its Feng labels' cells of zone z at
               E14.5 and at E9.5; defined where the group has at least min_cells_ext cells in z at both ends.
  zlv          the cardiomyocyte group (CM) takes its within-zone direction in LV; every other group keeps the entry's
               pooled direction. Void (the run refuses) if CM has no LV direction.
  zbal         every group takes the unweighted mean, over the zones where it has a within-zone direction, of those
               directions (per gene, over the zones where the gene's value is finite; NaN where none is); a group with
               no such zone keeps the entry's pooled direction, and the task file says so.
The task file records, per form and group: the zones used; the slope through the origin and the Pearson correlation of
the form's direction on the entry's over the moved genes (finite in both, not in M1's held list); the form's squared
size over the entry's on those genes; and the sha256 of each direction. Ratios and names only: no cell count is
written or printed. No Feng stage other than the bracket's two ends is read (share_read.feng_inputs,
feng_dirs.bracket_stats).

  run    One truth, 10 seeds. pooled, lib1, lib2, lib3: the chamber cut T1/E9.5->E10.5:chamber exactly as T46's run
         makes them (lib_step_read.load_truth); whole: the cut T1/E9.5->E10.5 with its pooled target exactly as T44's
         run loads it. Arms per seed, all from one cell draw: b0, m1_r_gmom (fraction 0.2), m1_r_gmom_x15 and
         m1_r_gmom_x2 (0.3 and 0.4) as T44 and T46 build them, and m1_r_gmom_zlv, m1_r_gmom_zlv_x2, m1_r_gmom_zbal and
         m1_r_gmom_zbal_x2: the builds of m1_r_gmom and m1_r_gmom_x2 (gmom_combo.m1_r_arms) with the form's directions
         in place of the entry's and nothing else changed (the same groups, held genes, R1, R2, composition, count,
         streams and g_mom factors). Refused on an existing task file, unless the committed N22 passed, and unless
         experiments/harness/TRIES.md holds row T53.
  read   T53's reading of the five truths' files (below, fixed in the row before the run). Writes t53.json and t53.txt
         and refuses to overwrite them.

The reading. Void unless the five truths hold the ten seeds at one commit on a clean tree (tree_clean exactly true); no
file records an anchor failure; every arm is on every row with its fraction; on the four chamber truths b0, m1_r_gmom,
m1_r_gmom_x15 and m1_r_gmom_x2 equal, seed for seed, their committed scores in
experiments/t1_m1/results/t46_m1_lib_step/tasks/, and on whole theirs in experiments/t1_m1/results/t44_m1_step/tasks/
(the whole cut's file), each to 1e-9; no new arm scores exactly as the entry's arm of its fraction on every seed of a
truth; each new arm's recorded directions differ from the entry's in at least one group (zlv: in CM and in no other).
Tests. On the two cuts (pooled and whole): the adoption rule's test (vec.adoption.keep_change, both cuts) of each new
arm against the pooled arm of its fraction: "passes"; with the sign reversed "loses"; else "neither". On each library
truth: T46's single-task test (ahead, behind, level). The statement, per form: "zone-fixed scores better" if the form
passes at both fractions (_zlv against m1_r_gmom and _zlv_x2 against m1_r_gmom_x2), "zone-fixed scores worse" if it
loses at both, otherwise "no difference at the bar", with the four comparisons' numbers and each library truth's labels.
Reported: per truth every arm's mean board score; each new arm's paired difference from the pooled arm of its fraction
with seeds ahead, per-seed range and sd; every metric's mean skill difference and the parts on the two de terms and on
mmd_u and the variogram; beside each new arm, per group (the route gives no covered-cell weights), the recorded slope of
its direction on the entry's and, from the uniform arms' means on the same truth, the gain a uniform multiple of the
entry's step equal to that slope (times 2 for an _x2 arm) would have by linear interpolation between the multiples 1,
1.5 and 2, "what size alone would give", read by nothing and not computed outside 1 to 2. The read puts nothing forward
by itself.

NOT an adoption result: a labelled in-sample diagnostic, as T44 to T46 are; Feng's E9.5 cells lie inside these
pseudo-tasks' external windows and the cuts are ones earlier T1 choices were made on. The result can only go to the team
captain.

Limits, repeated in the text:
- One released day, in-sample, on cuts earlier choices were made on.
- A zone-fixed direction differs from the pooled one in size as well as in content, and larger steps score better on
  these cuts (tries T44 to T46), so a gain or loss can be size: the slope and the interpolated figure are printed for
  that and decide nothing.
- Within-zone directions rest on fewer cells than the pooled one and are noisier.
- Which zone mix the hidden test's cells have is not known from released data, and neither form claims to match it.

Usage: on the DGX campaign checkout, repo root, all five truths at one commit, committing nothing until the five task
files exist, and on the machine T44 and T46 ran on, so that the rescored arms reproduce. The environment variables are
docs/COMPUTE.md's for repository Python on the DGX.
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \
      uv run --frozen --no-sync python -I -B experiments/t1_m1/zone_read.py run --truth pooled --out t53_m1_zone
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \
      uv run --frozen --no-sync python -I -B experiments/t1_m1/zone_read.py run --truth lib1 --out t53_m1_zone
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \
      uv run --frozen --no-sync python -I -B experiments/t1_m1/zone_read.py run --truth lib2 --out t53_m1_zone
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \
      uv run --frozen --no-sync python -I -B experiments/t1_m1/zone_read.py run --truth lib3 --out t53_m1_zone
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \
      uv run --frozen --no-sync python -I -B experiments/t1_m1/zone_read.py run --truth whole --out t53_m1_zone
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \
      uv run --frozen --no-sync python -I -B experiments/t1_m1/zone_read.py read --out t53_m1_zone
"""

from __future__ import annotations

import argparse
import hashlib
import json
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
sys.path.insert(0, str(HERE))
import eb_read as er  # noqa: E402
import feng_dirs  # noqa: E402
import gmom_combo as gc  # noqa: E402
import lib_step_read as lsr  # noqa: E402
import share_read  # noqa: E402
import step_read as sr  # noqa: E402

dl, gv, hygiene, n12 = feng_dirs.dl, gc.gv, gc.hygiene, er.n12
ROW = "T53"
STEM = "t53"
BRACKET = feng_dirs.BRACKET  # ("E9.5", "E14.5")
CM, LV = "CM", "LV"  # the cardiomyocyte group and the zone zlv fixes
FORMS = ("zlv", "zbal")
TASK_CHAMBER, TASK_WHOLE = gc.TASK_CHAMBER, gc.TASK_WHOLE
WHOLE = "whole"
CHAMBER_TRUTHS = lsr.TRUTHS  # pooled, lib1, lib2, lib3
LIBS = tuple(lsr.LIBRARIES)
TRUTHS = (*CHAMBER_TRUTHS, WHOLE)
CUTS = (lsr.POOLED, WHOLE)  # the adoption rule's two cuts
TASK_OF = {**dict.fromkeys(CHAMBER_TRUTHS, TASK_CHAMBER), WHOLE: TASK_WHOLE}
BASE, X15, X2 = sr.BASE, "m1_r_gmom_x15", "m1_r_gmom_x2"
UNIFORM_MULTIPLIERS = {X15: 1.5, X2: 2.0}
UNIFORM = {BASE: 1.0, **UNIFORM_MULTIPLIERS}  # the uniform arms: multiple of the entry's step
NEW_ARMS = {f"{BASE}_{f}{s}": (f, p) for f in FORMS for s, p in (("", BASE), ("_x2", X2))}  # arm -> (form, pooled arm)
ARMS = ("b0", BASE, X15, X2, *NEW_ARMS)
REPRO_ARMS = ("b0", BASE, X15, X2)
FRACTION = sr.FRACTION  # 0.2
FRACTION_TOL = sr.FRACTION_TOL
REPRO_TOL = sr.REPRO_TOL  # 1e-9
COMMITTED = {**{k: RESULTS / "t46_m1_lib_step" / "tasks" / lsr.TASK_FILE.format(truth=k) for k in CHAMBER_TRUTHS},
             WHOLE: RESULTS / "t44_m1_step" / "tasks" / gc.CUT_FILES[TASK_WHOLE]}
TASK_FILE = "{truth}.json"
PASSES, LOSES, NEITHER = "passes", "loses", "neither"
AHEAD, BEHIND, LEVEL = "ahead", "behind", "level"
BETTER = "zone-fixed scores better"
WORSE = "zone-fixed scores worse"
NO_DIFFERENCE = "no difference at the bar"
LABEL_T53 = ("try T53, a labelled in-sample diagnostic, never an adoption result: M1 (with R1, R2 and g_mom's factors, "
             "the entry's form) with its Feng direction zone-fixed (zlv: the cardiomyocyte group's LV-zone direction; "
             "zbal: every group's unweighted mean of its within-zone directions) at the step and twice the step, on "
             "T44's and T46's routes (both E10.5 cuts, the chamber cut also library by library); Feng's E9.5 cells lie "
             "inside these pseudo-tasks' external windows and the cuts are ones earlier T1 choices were made on; the "
             "result can only go to the team captain")
LIMITS = (
    "One released day, in-sample, on cuts earlier choices were made on.",
    "A zone-fixed direction differs from the pooled one in size as well as in content, and larger steps score better "
    "on these cuts (tries T44 to T46), so a gain or loss can be size: the slope and the interpolated figure are "
    "printed for that and decide nothing.",
    "Within-zone directions rest on fewer cells than the pooled one and are noisier.",
    "Which zone mix the hidden test's cells have is not known from released data, and neither form claims to match it.",
)


# ---------------------------------------------------------------------------------------------------------------------
# The directions


def key_parts(st: Mapping) -> list[list[str]]:
    """Each statistics key split on '|'; a key that is not stage|zone|cell type|half (four parts) is refused."""
    parts = [str(k).split("|") for k in st["keys"]]
    bad = [str(k) for k, p in zip(st["keys"], parts, strict=True) if len(p) != 4]
    if bad:
        raise ValueError(f"statistics keys without the four parts stage|zone|cell type|half: {bad[:3]}")
    return parts


def zones_of(st: Mapping) -> list[str]:
    """The zones the statistics' keys name, sorted."""
    return sorted({p[1] for p in key_parts(st)})


def zone_sel(st: Mapping, stage: str, zone: str, labels, halves: Sequence[str] = ("A", "B")) -> dict | None:
    """dirlib.feng_sel restricted to one zone: the stage's cells of `labels` in `zone`, both halves pooled."""
    ls = set(labels)
    return dl.agg(st, lambda p: p[0] == stage and p[1] == zone and p[2] in ls and p[3] in halves)


def within_zone_directions(st: Mapping, genes, cfg: Mapping, groups,
                           bracket: tuple[str, str] = BRACKET) -> dict[str, dict[str, np.ndarray]]:
    """group -> zone -> the log ratio of the group's mean CP10k over its Feng labels' cells of the zone, bracket[1]
    against bracket[0], on `genes` (NaN where unmapped); a zone where the group has under min_cells_ext cells at either
    end gives none. `groups`: the entry's direction groups."""
    zones = zones_of(st)
    idx = dl.sym_index(st["genes"], list(genes))
    eps, mc = cfg["metrics"]["pseudocount_cp10k"], cfg["metrics"]["min_cells_ext"]
    out: dict[str, dict[str, np.ndarray]] = {}
    for g in sorted(groups):
        fl = cfg["external_labels"][g][0]
        per = {}
        for z in zones:
            a, b = zone_sel(st, bracket[0], z, fl), zone_sel(st, bracket[1], z, fl)
            if a is None or b is None or min(a["n"], b["n"]) < mc:
                continue
            per[z] = dl.on_t1(dl.lfc_ratio(a, b, eps), idx)
        out[g] = per
    return out


def zlv_directions(entry: Mapping[str, np.ndarray], within: Mapping) -> tuple[dict, dict]:
    """CM's LV direction, every other group the entry's own array; and the zones used per group ([] = pooled)."""
    if LV not in within.get(CM, {}):
        raise ValueError(f"T53 is void: the cardiomyocyte group {CM} has no {LV} direction")
    out = dict(entry)
    out[CM] = within[CM][LV]
    return out, {g: ([LV] if g == CM else []) for g in entry}


def zone_mean(vectors: Sequence[np.ndarray]) -> np.ndarray:
    """Per gene, the unweighted mean over the vectors whose value is finite there; NaN where none is."""
    V = np.vstack([np.asarray(v, dtype=np.float64) for v in vectors])
    fin = np.isfinite(V)
    cnt = fin.sum(axis=0)
    total = np.where(fin, V, 0.0).sum(axis=0)
    return np.where(cnt > 0, total / np.maximum(cnt, 1), np.nan)


def zbal_directions(entry: Mapping[str, np.ndarray], within: Mapping) -> tuple[dict, dict]:
    """Every group the unweighted mean of its within-zone directions (zone_mean); a group with none keeps the entry's
    own array (zones used [])."""
    out, used = {}, {}
    for g, v in entry.items():
        zs = sorted(within.get(g, {}))
        out[g] = zone_mean([within[g][z] for z in zs]) if zs else v
        used[g] = zs
    return out, used


def form_directions(entry: Mapping, within: Mapping) -> dict[str, tuple[dict, dict]]:
    return {"zlv": zlv_directions(entry, within), "zbal": zbal_directions(entry, within)}


def direction_sha256(v) -> str:
    x = np.ascontiguousarray(np.asarray(v, dtype=np.float64))
    h = hashlib.sha256()
    h.update(f"{x.dtype.str}{x.shape}".encode())
    h.update(x.tobytes())
    return h.hexdigest()


def _num(x) -> float | None:
    return float(x) if x is not None and np.isfinite(x) else None


def compare_directions(form: np.ndarray, entry: np.ndarray, genes: Sequence[str], held: Sequence[str]) -> dict:
    """Over the moved genes (finite in both, not held): the slope through the origin and the Pearson correlation of the
    form's direction on the entry's (share_read.origin_slope), and the form's squared size over the entry's."""
    held_set = set(map(str, held))
    f, e = np.asarray(form, dtype=np.float64), np.asarray(entry, dtype=np.float64)
    moved = np.isfinite(f) & np.isfinite(e) & np.array([str(g) not in held_set for g in genes], dtype=bool)
    slope, r = share_read.origin_slope(f[moved], e[moved])
    ee = float(e[moved] @ e[moved])
    return {"slope": _num(slope), "pearson": _num(r),
            "size_ratio": _num(float(f[moved] @ f[moved]) / ee) if ee > 0 else None}


def direction_record(entry: Mapping, forms: Mapping[str, tuple[dict, dict]], genes, held, zones: Sequence[str]) -> dict:
    """The task file's direction record: names, ratios and hashes only."""
    rec = {"zones": list(zones), "entry": {g: {"sha256": direction_sha256(entry[g])} for g in sorted(entry)}}
    for form, (dirs, used) in forms.items():
        rec[form] = {g: {"zones": list(used[g]), "pooled": not used[g], **compare_directions(dirs[g], entry[g], genes,
                                                                                            held),
                         "sha256": direction_sha256(dirs[g])} for g in sorted(dirs)}
    rec["zbal_fallback"] = sorted(g for g, v in rec["zbal"].items() if v["pooled"])
    return rec


# ---------------------------------------------------------------------------------------------------------------------
# The run


def arm_fractions(fraction: float) -> dict[str, float]:
    fr = sr.arm_fractions(fraction, UNIFORM_MULTIPLIERS)
    return {**fr, **{a: fr[p] for a, (_, p) in NEW_ARMS.items()}}


def load_truth(t: harness.PseudoTask, data: Path, genes, truth: str) -> Stage:
    """whole: the cut's target exactly as T44's run loads it; a chamber truth: as T46's run makes it."""
    if truth == WHOLE:
        return (harness.filter_target(t, data, genes)[0] if t.target_filter is not None else
                Stage.from_anndata(read_stage(data / t.target, genes), genes))
    return lsr.load_truth(t, data, genes, truth)


def build_seed(ctx: harness.BuildContext, b0, entry: Mapping, forms: Mapping[str, Mapping], groups: Mapping,
               fraction: float, held_m1: Sequence[str], held_factors: Sequence[str]) -> dict[str, Stage]:
    """One seed's arms: b0, m1_r_gmom and the uniform arms (step_read.step_arms), then each new arm, the build of the
    pooled arm of its fraction with the form's directions in place of the entry's and every other argument the same."""
    preds = {"b0": b0(ctx), **sr.step_arms(ctx, entry, groups, fraction, held_m1, held_factors,
                                           multipliers=UNIFORM_MULTIPLIERS)}
    fr = arm_fractions(fraction)
    for arm, (form, _) in NEW_ARMS.items():
        preds[arm] = gc.m1_r_arms(ctx, forms[form], groups, fr[arm], held_m1, held_factors)[BASE]
    return preds


def cmd_run(a) -> None:
    if ROW not in n12.tries_rows():
        sys.exit(f"T53 is refused: {ROW} has no row in experiments/harness/TRIES.md (a row comes before its run)")
    n22 = RESULTS / a.n22 / "n22.json"
    if not gc.n22_passed(n22):
        sys.exit(f"T53 is refused: {gv._rel(n22)} does not hold a passing N22 (T26's route needs N22 passed)")
    out = RESULTS / a.out
    name = TASK_FILE.format(truth=a.truth)
    gv.refuse_existing(out / "tasks", (name,))
    started = gv.tree_state()
    seeds = gv.parse_seeds(a.seeds)
    data = Path(a.data)
    t = {x.key: x for x in harness.load_tasks()["T1:test"]}[TASK_OF[a.truth]]
    boards = load_boards()
    genes = boards[t.panel_board].genes
    truth = load_truth(t, data, genes, a.truth)
    ref_stage = Stage.from_anndata(read_stage(data / t.reference, genes), genes)
    board = harness.task_board(t, boards, truth.n, ref_stage.n)
    stamp_inputs = [data / t.target, data / t.reference, *[data / f for f in t.inputs.values()]]
    res = {"label": LABEL_T53, "task": t.key, "truth": a.truth, "seeds": list(seeds)}
    floors, failures = lsr.all_anchors(t, board, truth, ref_stage, seeds)
    if failures:
        res.update({"anchor_failures": failures, "rows": [], "stamp": gv.stamp(started, stamp_inputs)})
        (out / "tasks").mkdir(parents=True, exist_ok=True)
        (out / "tasks" / name).write_text(json.dumps(res, indent=1, default=float) + "\n")
        sys.exit(f"T53 is refused on {a.truth}: the floor's anchors fail on seeds {[f['seed'] for f in failures]} "
                 "(recorded in the task file; T53 is void)")
    inputs = {r: load_stage(data / f, lineage=False) for r, f in t.inputs.items()}
    pool = inputs[harness.base_role(t)]
    pool_genes = list(pool.var_names)
    cfg = feng_dirs.load_config()
    st, src = share_read.feng_inputs()
    entry, _, _ = feng_dirs.directions_from_stats(st, pool_genes, cfg, BRACKET)
    try:
        forms = form_directions(entry, within_zone_directions(st, pool_genes, cfg, entry, BRACKET))
    except ValueError as e:
        sys.exit(f"T53 is refused: {e}")
    groups = feng_dirs.label_groups(pool.obs["celltype"].astype(str).unique())
    fraction = feng_dirs.bracket_fraction(harness.STAGE_DAYS[t.reference], harness.STAGE_DAYS[t.target])
    preregistered = t1_scope.held_flat_genes()
    union = hygiene.load_union()
    held_m1 = gc.held_feng(union, preregistered)
    directions = direction_record(entry, forms, pool_genes, held_m1, zones_of(st))
    form_dirs = {f: d for f, (d, _) in forms.items()}
    view, cboard = harness.TaskView.of(t), harness._candidate_board(board)
    b0 = harness._b0_t1(harness.B0_DAMPING)
    rows = []
    for s in seeds:
        t0 = time.time()
        f = floors[s]
        ctx = harness.BuildContext(task=view, board=cboard, inputs=inputs, seed=harness.BUILD_SEED_OFFSET + s,
                                   derived={"t1_held_flat": preregistered})
        preds = build_seed(ctx, b0, entry, form_dirs, groups, fraction, held_m1, preregistered)
        arms = er.score_arms(board, truth, ref_stage, s, f.anchors, preds)
        rows.append({"seed": s, "floor": float(f.score), "arms": arms,
                     "sha256": {k: lsr.prediction_sha256(p) for k, p in preds.items()},
                     "seconds": round(time.time() - t0, 1)})
        print(f"{a.truth} seed {s}: " + ", ".join(f"{k} {v['score']:.2f}" for k, v in arms.items()), flush=True)
    external = {**src, "bracket": list(BRACKET), "pseudocount_cp10k": cfg["metrics"]["pseudocount_cp10k"],
                "min_cells_ext": cfg["metrics"]["min_cells_ext"], "groups": sorted(entry)}
    res.update({"fraction": fraction, "fractions": arm_fractions(fraction), "new_arms": {a_: list(v) for a_, v in
                                                                                       NEW_ARMS.items()},
                "external": external, "directions": directions, "rows": rows,
                "stamp": gv.stamp(started, [*stamp_inputs, t1_scope.HELD_FLAT_GENES, ts.R1_CONFIG,
                                            REPO_ROOT / union["file"], feng_dirs.clean_stats_path(), n22])})
    (out / "tasks").mkdir(parents=True, exist_ok=True)
    (out / "tasks" / name).write_text(json.dumps(res, indent=1, default=float) + "\n")


# ---------------------------------------------------------------------------------------------------------------------
# The reading


def reproduction(rows: Sequence[Mapping], committed: Sequence[Mapping] | None) -> dict:
    """Per seed, each of REPRO_ARMS against its committed score; holds when every seed and arm is there and each is
    within REPRO_TOL."""
    by_seed = {r["seed"]: r for r in committed or []}
    worst = dict.fromkeys(REPRO_ARMS, 0.0)
    absent = []
    for r in rows:
        old = by_seed.get(r["seed"])
        for arm in REPRO_ARMS:
            if old is None or arm not in old["arms"] or arm not in r["arms"]:
                absent.append(f"{r['seed']}:{arm}")
                continue
            worst[arm] = max(worst[arm], abs(r["arms"][arm]["score"] - old["arms"][arm]["score"]))
    holds = bool(rows) and bool(committed) and not absent and all(w <= REPRO_TOL for w in worst.values())
    return {"largest_difference": worst, "absent": absent, "holds": holds}


def differing_groups(directions: Mapping, form: str) -> list[str] | None:
    """The groups whose recorded direction for `form` differs from the entry's (by sha256); None without a record."""
    try:
        entry, rec = directions["entry"], directions[form]
        return sorted(g for g in set(entry) | set(rec)
                      if (entry.get(g) or {}).get("sha256") != (rec.get(g) or {}).get("sha256"))
    except (KeyError, TypeError):
        return None


def validity(results: Mapping[str, Mapping], committed: Mapping[str, Sequence[Mapping]]) -> dict:
    """The void conditions. results: truth -> task file; committed: truth -> the committed file's rows (absent where the
    file is missing)."""
    reasons = []
    missing = [k for k in TRUTHS if k not in results]
    if missing:
        reasons.append(f"missing truths {missing}")
    want_fr = arm_fractions(FRACTION)
    for k, r in results.items():
        if r.get("anchor_failures"):
            reasons.append(f"{k}: the floor's anchors invert or go missing on seeds "
                           f"{[f['seed'] for f in r['anchor_failures']]}")
        rows = r.get("rows") or []
        seeds = [row["seed"] for row in rows]
        if sorted(seeds) != list(SEEDS) or list(r.get("seeds", seeds)) != list(SEEDS):
            reasons.append(f"{k}: seeds {seeds}, not {list(SEEDS)}")
        absent = sorted({a for row in rows for a in ARMS if a not in row["arms"]})
        if absent:
            reasons.append(f"{k}: arms {absent} missing")
        if r.get("task") != TASK_OF[k]:
            reasons.append(f"{k}: the file holds task {r.get('task')!r}")
        if r.get("truth") != k:
            reasons.append(f"{k}: the file holds truth {r.get('truth')!r}")
        if not isinstance(r.get("fraction"), int | float) or abs(r["fraction"] - FRACTION) > FRACTION_TOL:
            reasons.append(f"{k}: fraction {r.get('fraction')!r}, not {FRACTION:g}")
        got_fr = r.get("fractions") or {}
        bad = [a for a, f in want_fr.items()
               if not isinstance(got_fr.get(a), int | float) or abs(got_fr[a] - f) > FRACTION_TOL]
        if bad:
            reasons.append(f"{k}: wrong fraction on {bad}")
        rep = reproduction(rows, committed.get(k))
        if not rep["holds"]:
            reasons.append(f"{k}: {', '.join(REPRO_ARMS)} do not reproduce their committed scores within {REPRO_TOL:g} "
                           "per seed")
        if rows and not absent:
            same = [a for a, (_, p) in NEW_ARMS.items()
                    if all(row["arms"][a]["score"] == row["arms"][p]["score"] for row in rows)]
            if same:
                reasons.append(f"{k}: new arms score exactly as the pooled arm of their fraction on every seed: {same}")
        for form in FORMS:
            diff = differing_groups(r.get("directions") or {}, form)
            if diff is None:
                reasons.append(f"{k}: no direction record for {form}")
            elif form == "zlv" and diff != [CM]:
                reasons.append(f"{k}: zlv's directions differ from the entry's in {diff}, not in {CM} alone")
            elif not diff:
                reasons.append(f"{k}: {form}'s directions are the entry's in every group")
    commits = {(r.get("stamp") or {}).get("commit") for r in results.values()}
    unclean = sorted(k for k, r in results.items() if (r.get("stamp") or {}).get("tree_clean") is not True)
    if unclean:
        reasons.append(f"tree not clean on {unclean}")
    if len(commits) > 1:
        reasons.append(f"the truths ran at {len(commits)} commits")
    return {"valid": not reasons, "reasons": reasons, "commits": sorted(c for c in commits if c)}


def cut_test(results: Mapping[str, Mapping], arm: str, pooled: str) -> dict:
    """The adoption rule on both cuts (step_read.compare): 'passes'; with the sign reversed 'loses'; else 'neither'."""
    cuts = {k: results[k] for k in CUTS}
    fwd, rev = sr.compare(cuts, pooled, arm), sr.compare(cuts, arm, pooled)
    test = PASSES if fwd["passes"] else LOSES if rev["passes"] else NEITHER
    return {"base": pooled, "change": arm, "test": test, "reason": fwd["rule"]["reason"],
            "reversed_reason": rev["rule"]["reason"],
            "cuts": {k: {"mean": c["mean"], "ahead": c["ahead"], "of": c["of"], "range": c["range"]}
                     for k, c in fwd["cuts"].items()}}


def library_label(c: Mapping) -> str:
    return AHEAD if c["ahead"] else BEHIND if c["behind"] else LEVEL


def form_statement(tests: Sequence[str]) -> str:
    if all(t == PASSES for t in tests):
        return BETTER
    if all(t == LOSES for t in tests):
        return WORSE
    return NO_DIFFERENCE


def size_alone(means: Mapping[str, float], multiple: float | None, pooled: str) -> float | None:
    """The gain over the pooled arm of the new arm's fraction that a uniform multiple of the entry's step would have,
    by linear interpolation of the uniform arms' mean scores between the multiples 1, 1.5 and 2; None outside."""
    if multiple is None or not np.isfinite(multiple) or not 1.0 <= multiple <= 2.0:
        return None
    xs = list(UNIFORM.values())
    ys = [means[a] for a in UNIFORM]
    return float(np.interp(multiple, xs, ys)) - float(means[pooled])


def t53_reading(results: Mapping[str, Mapping], committed: Mapping[str, Sequence[Mapping]]) -> dict:
    v = validity(results, committed)
    out = {"valid": v["valid"], "validity": v}
    if not v["valid"]:
        out["outcome"] = "VOID"
        return out
    means = {k: {a: float(np.mean([row["arms"][a]["score"] for row in results[k]["rows"]])) for a in ARMS}
             for k in TRUTHS}
    cut_tests = {a: cut_test(results, a, p) for a, (_, p) in NEW_ARMS.items()}
    comps = {k: {a: lsr.truth_comparison(results[k], k, a, base=p) for a, (_, p) in NEW_ARMS.items()} for k in TRUTHS}
    labels = {k: {a: library_label(comps[k][a]) for a in NEW_ARMS} for k in LIBS}
    statements = {}
    for form in FORMS:
        arms = [a for a, (f, _) in NEW_ARMS.items() if f == form]
        statements[form] = {"statement": form_statement([cut_tests[a]["test"] for a in arms]),
                            "comparisons": {a: cut_tests[a] for a in arms},
                            "library_labels": {k: {a: labels[k][a] for a in arms} for k in LIBS}}
    size = {}
    for k in TRUTHS:
        rec = results[k]["directions"]
        size[k] = {a: {g: {"slope": rec[f][g]["slope"],
                           "what_size_alone_would_give": size_alone(
                               means[k], None if rec[f][g]["slope"] is None else rec[f][g]["slope"] * UNIFORM[p], p)}
                       for g in sorted(rec[f])} for a, (f, p) in NEW_ARMS.items()}
    out.update({"outcome": {f: s["statement"] for f, s in statements.items()}, "statements": statements,
                "mean_scores": means, "comparisons": comps, "size_alone": size,
                "directions": {f: results[lsr.POOLED]["directions"][f] for f in FORMS},
                "zbal_fallback": results[lsr.POOLED]["directions"].get("zbal_fallback", [])})
    return out


def _f(x, fmt="+.2f") -> str:
    return "n/a" if x is None else format(x, fmt)


def t53_text(res: Mapping) -> str:
    """The reading as text; no cell count is printed."""
    r = res["reading"]
    lines = [res["label"]]
    if not r["valid"]:
        lines.append("Outcome: VOID")
        lines.append("Void: " + "; ".join(r["validity"]["reasons"]))
    else:
        for form, s in r["statements"].items():
            lines.append(f"{form}: {s['statement']}")
            for a, c in s["comparisons"].items():
                lines.append(f"    {a} against {c['base']} on both cuts: {c['test']} ({c['reason']}); "
                             + "; ".join(f"{k} {x['mean']:+.2f} points, ahead {x['ahead']} of {x['of']}, per seed "
                                         f"{x['range'][0]:+.2f} to {x['range'][1]:+.2f}" for k, x in c["cuts"].items()))
            lines.append("    library truths (single-task test against the pooled arm of the same fraction): "
                         + "; ".join(f"{k}: " + ", ".join(f"{a} {lab}" for a, lab in ls_.items())
                                     for k, ls_ in s["library_labels"].items()))
        if r["zbal_fallback"]:
            lines.append(f"zbal keeps the entry's pooled direction for {r['zbal_fallback']} (no within-zone direction)")
        for form, rec in r["directions"].items():
            lines.append(f"{form}'s directions on the entry's (moved genes): " + "; ".join(
                f"{g}: zones {v['zones'] or ['pooled']}, slope {_f(v['slope'], '.3f')}, pearson "
                f"{_f(v['pearson'], '.3f')}, squared size ratio {_f(v['size_ratio'], '.3f')}" for g, v in rec.items()))
        for k in TRUTHS:
            lines.append(f"truth {k}: mean scores " + ", ".join(f"{a} {x:.2f}" for a, x in r["mean_scores"][k].items()))
            for a, c in r["comparisons"][k].items():
                lines.append("    " + lsr._comparison_line(c))
                lines.append(f"    {a}: what size alone would give (reported, read by nothing): " + ", ".join(
                    f"{g} slope {_f(v['slope'], '.3f')} -> {_f(v['what_size_alone_would_give'])}"
                    for g, v in r["size_alone"][k][a].items()))
    lines.extend(f"Limit: {x}" for x in LIMITS)
    lines.append("Limit: in-sample labelled diagnostic; never adoption evidence; puts nothing forward by itself; it "
                 "can only go to the team captain")
    lines.append(f"Commit {res['stamp']['commit']}, tree clean {res['stamp']['tree_clean']}, UTC {res['stamp']['utc']}")
    return "\n".join(lines) + "\n"


def read(out: Path, committed_files: Mapping[str, Path] | None = None) -> dict:
    """Read the truths present in out/tasks against the committed T46 and T44 files; write t53.{json,txt}."""
    gv.refuse_existing(out, (f"{STEM}.json", f"{STEM}.txt"))
    started = gv.tree_state()
    files = COMMITTED if committed_files is None else committed_files
    present = {k: out / "tasks" / TASK_FILE.format(truth=k) for k in TRUTHS
               if (out / "tasks" / TASK_FILE.format(truth=k)).exists()}
    results = {k: json.loads(p.read_text()) for k, p in present.items()}
    found = {k: Path(p) for k, p in files.items() if Path(p).exists()}
    committed = {k: json.loads(p.read_text())["rows"] for k, p in found.items()}
    res = {"label": LABEL_T53, "reading": t53_reading(results, committed),
           "stamps": {k: v.get("stamp") for k, v in results.items()}}
    res["stamp"] = gv.stamp(started, [*present.values(), *found.values()])
    gv.write_outputs(out, STEM, res, t53_text(res))
    return res


def cmd_read(a) -> None:
    res = read(RESULTS / a.out)
    if res["reading"]["outcome"] == "VOID":
        sys.exit("T53's reading is void and goes to the captain")


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("run", help="one truth of try T53")
    p.add_argument("--truth", required=True, choices=list(TRUTHS))
    p.add_argument("--out", required=True)
    p.add_argument("--n22", default="n22_additivity")
    p.add_argument("--seeds", default=",".join(map(str, SEEDS)))
    p.add_argument("--data", default=str(harness.DATA_DIR))
    p = sub.add_parser("read", help="try T53's reading")
    p.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    {"run": cmd_run, "read": cmd_read}[a.cmd](a)


if __name__ == "__main__":
    main()
