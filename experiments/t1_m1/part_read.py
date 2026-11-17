#!/usr/bin/env python
"""Try T55 (experiments/harness/TRIES.md): M1's Feng step with the part all cell groups share and the group-specific
parts sized apart. A labelled in-sample read on the routes of try T44 (experiments/t1_m1/step_read.py, the whole cut)
and try T46 (experiments/t1_m1/lib_step_read.py, the chamber cut pooled and library by library). Built to brief AI.

M1 moves each covered cell along d_g, per cell group the log ratio of mean CP10k between Feng's E9.5 and E14.5 hearts
(feng_dirs.directions_from_stats). Try T51 (split_read.py) split each d_g into c, the unweighted mean of the groups'
directions on M, and s_g = d_g - c. This read asks whether the content of the step matters at a fixed amount of
movement: whether the group-specific parts carry any of the score, and whether a step along c alone scores better or
worse than the entry's at equal movement.

Directions (per group, from the entry's directions alone):
  moved    a group's moved genes: those where d_g is finite and that are not in M1's held list (zone_read.moved_mask).
  M        split_read.common_mask of the groups' moved genes; c and s_g by split_read.split_vectors (c the unweighted
           mean over every direction group on M). The run refuses unless every group's moved genes are M (then it
           stops before any build and writes no task file; try T51's record says the real statistics meet this).
  form     (a, b): per group the direction a_g x c + b x s_g on M and the entry's own values off M (there d_g is NaN
           or the gene is held, so M1 does not move it).
  a_n,g    sqrt(sum(d_g^2) / sum(c^2)) over M: the common part alone at the entry's squared size.
  a_p,g    sum(d_g^2) / sum(c d_g) over M: the common part alone at the entry's projection on d_g. The run refuses
           where sum(c d_g) is not positive in a group (then sum(c^2) is positive wherever it runs).
  The task file records per group a_n,g, a_p,g, the energy parts (split_read.energy_parts: sum(c^2), sum(s_g^2) and
  2 sum(c s_g), each over sum(d_g^2)), the largest deviation of the form (1, 1) from d_g on M, and for every form its
  slope through the origin on the entry's direction over the moved genes, its squared size over the entry's and its
  sha256 (zone_read.compare_directions, zone_read.direction_sha256). Ratios, names and hashes only.

  run    One truth, 10 seeds. pooled, lib1, lib2, lib3: the chamber cut exactly as T46's run makes them; whole: the
         other cut exactly as T44's run loads it (zone_read.load_truth). Arms per seed, all from one cell draw: b0,
         m1_r_gmom (fraction 0.2), m1_r_gmom_x15, m1_r_gmom_x2 as T44 and T46 build them (step_read.step_arms); and
         seven new arms, each the build of m1_r_gmom (0.2) or m1_r_gmom_x2 (0.4) (gmom_combo.m1_r_arms) with a form's
         directions in place of the entry's and nothing else changed: m1_r_gmom_c1 (1, 0), m1_r_gmom_s05 (1, 0.5),
         m1_r_gmom_s15 (1, 1.5), m1_r_gmom_cn and m1_r_gmom_cn_x2 (a_n,g, 0), m1_r_gmom_cp and m1_r_gmom_cp_x2 (a_p,g,
         0). Refused on an existing task file, unless the committed N22 passed, unless experiments/harness/TRIES.md
         holds row T55, and where the directions refuse (above).
  read   T55's reading of the five truths' files (below, fixed in the row before the run). Writes t55.json and t55.txt
         and refuses to overwrite them.

The reading. Void unless: the five truths hold the ten seeds with every arm and its fraction at one commit on a clean
tree; no file records an anchor failure; on the four chamber truths b0, m1_r_gmom, m1_r_gmom_x15 and m1_r_gmom_x2 equal,
seed for seed, their committed scores in experiments/t1_m1/results/t46_m1_lib_step/tasks/, and on whole theirs in
experiments/t1_m1/results/t44_m1_step/tasks/, each to 1e-9; every recorded direction's sha256 and every a_n,g and a_p,g
are the same in the five files; the form (1, 1) rebuilt from c and s_g equals the entry's direction to 1e-12 in every
group (recorded by the run); no new arm scores exactly as the pooled arm of its fraction on every seed of a truth.
Tests. On the two cuts (pooled and whole) the adoption rule's (vec.adoption.keep_change: passes; with the sign
reversed, loses; else neither); on a library truth T46's single-task test (ahead, behind, level).
Statement 1, the group-specific parts at the entry's size: m1_r_gmom_c1 against m1_r_gmom: "dropping them scores
better" if it passes, "they carry part of the score" if it loses, otherwise "no difference at the bar"; _s05 and _s15
against m1_r_gmom with the same labels.
Statement 2, content at equal movement: m1_r_gmom_cn against m1_r_gmom and m1_r_gmom_cn_x2 against m1_r_gmom_x2: "the
common part alone scores better at equal movement" if both pass, "the common part alone scores worse at equal
movement" if both lose, "no difference at the bar" if both are neither, otherwise "mixed" with both labels.
Statement 3, content at equal projection, the same form and read by nothing: m1_r_gmom_cp and m1_r_gmom_cp_x2, and
beside each, per group, what the uniform arms give at the multiple with the same squared size (sqrt of the form's
squared size ratio, twice that for _x2; zone_read.size_alone, linear between the multiples 1, 1.5 and 2; "outside the
uniform arms" where it is not between them).
Reported: per truth every arm's mean score; each comparison's mean paired difference, seeds ahead, per-seed range and
library labels; the parts on the two de terms and on mmd_u and the variogram. The read puts nothing forward by itself.

NOT an adoption result: a labelled in-sample diagnostic, as T44 to T54 are; the result can only go to the team
captain. The limits are LIMITS below, repeated in the text.

Usage: on the DGX campaign checkout, repo root, all five truths at one commit, committing nothing until the five task
files exist, and on the machine T44 and T46 ran on, so that the rescored arms reproduce.
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \
      uv run --frozen --no-sync python -I -B experiments/t1_m1/part_read.py run --truth pooled --out t55_m1_part
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \
      uv run --frozen --no-sync python -I -B experiments/t1_m1/part_read.py run --truth lib1 --out t55_m1_part
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \
      uv run --frozen --no-sync python -I -B experiments/t1_m1/part_read.py run --truth lib2 --out t55_m1_part
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \
      uv run --frozen --no-sync python -I -B experiments/t1_m1/part_read.py run --truth lib3 --out t55_m1_part
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \
      uv run --frozen --no-sync python -I -B experiments/t1_m1/part_read.py run --truth whole --out t55_m1_part
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \
      uv run --frozen --no-sync python -I -B experiments/t1_m1/part_read.py read --out t55_m1_part
"""

from __future__ import annotations

import argparse
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
import split_read  # noqa: E402
import zone_read as zr  # noqa: E402

er, feng_dirs, gc, lsr, share_read, sr = zr.er, zr.feng_dirs, zr.gc, zr.lsr, zr.share_read, zr.sr
gv, hygiene, n12 = zr.gv, zr.hygiene, zr.n12
ROW = "T55"
STEM = "t55"
BRACKET = feng_dirs.BRACKET  # ("E9.5", "E14.5")
N, P = "a_n", "a_p"  # the per-group multiples
# form -> (a, b): a a number or the per-group multiple's name
FORMS = {"c1": (1.0, 0.0), "s05": (1.0, 0.5), "s15": (1.0, 1.5), "cn": (N, 0.0), "cp": (P, 0.0)}
REBUILD_TOL = 1e-12  # the form (1, 1) must equal the entry's direction on M to this
TASK_CHAMBER, TASK_WHOLE = zr.TASK_CHAMBER, zr.TASK_WHOLE
WHOLE = zr.WHOLE
CHAMBER_TRUTHS, LIBS, TRUTHS, CUTS, TASK_OF = zr.CHAMBER_TRUTHS, zr.LIBS, zr.TRUTHS, zr.CUTS, zr.TASK_OF
BASE, X15, X2 = zr.BASE, zr.X15, zr.X2
UNIFORM_MULTIPLIERS, UNIFORM = zr.UNIFORM_MULTIPLIERS, zr.UNIFORM
# new arm -> (its form, the pooled arm of its fraction)
NEW_ARMS = {f"{BASE}_c1": ("c1", BASE), f"{BASE}_s05": ("s05", BASE), f"{BASE}_s15": ("s15", BASE),
            f"{BASE}_cn": ("cn", BASE), f"{BASE}_cn_x2": ("cn", X2), f"{BASE}_cp": ("cp", BASE),
            f"{BASE}_cp_x2": ("cp", X2)}
ARMS = ("b0", BASE, X15, X2, *NEW_ARMS)
REPRO_ARMS = zr.REPRO_ARMS  # b0, m1_r_gmom, m1_r_gmom_x15, m1_r_gmom_x2
FRACTION, FRACTION_TOL, REPRO_TOL = zr.FRACTION, zr.FRACTION_TOL, zr.REPRO_TOL
COMMITTED = zr.COMMITTED
TASK_FILE = zr.TASK_FILE
PASSES, LOSES, NEITHER = zr.PASSES, zr.LOSES, zr.NEITHER
NO_DIFFERENCE = zr.NO_DIFFERENCE  # "no difference at the bar"
S1_ARMS = (f"{BASE}_c1", f"{BASE}_s05", f"{BASE}_s15")
S1_LABELS = {PASSES: "dropping them scores better", LOSES: "they carry part of the score", NEITHER: NO_DIFFERENCE}
S2_ARMS = (f"{BASE}_cn", f"{BASE}_cn_x2")
S3_ARMS = (f"{BASE}_cp", f"{BASE}_cp_x2")
BETTER_MOVE = "the common part alone scores better at equal movement"
WORSE_MOVE = "the common part alone scores worse at equal movement"
BETTER_PROJ = "the common part alone scores better at equal projection"
WORSE_PROJ = "the common part alone scores worse at equal projection"
MIXED = "mixed"
LABEL_T55 = ("try T55, a labelled in-sample diagnostic, never an adoption result: M1 (with R1, R2 and g_mom's factors, "
             "the entry's form) with its Feng direction split into the part all cell groups share (c) and the "
             "group-specific parts (s_g) as try T51 makes them, the parts sized apart at the step and twice the step, "
             "on T44's and T46's routes (both E10.5 cuts, the chamber cut also library by library); the cuts are ones "
             "earlier T1 choices were made on; the result can only go to the team captain")
CHOICES = (
    "Off M a form keeps the entry's own values (NaN or a held gene, which M1 does not move), so that the forms differ "
    "from the entry on M only.",
    "A group's moved genes are zone_read.moved_mask's (d_g finite and the gene not in M1's held list), the same "
    "definition as try T51's per-group mask; M is split_read.common_mask of them and c, s_g split_read.split_vectors' "
    "(c the unweighted mean over every direction group).",
    "Where a group's moved genes are not M, or sum(c d_g) is not positive in a group, the run stops before any build "
    "and writes no task file; the read is then void with the truth missing.",
    "Statement 2's 'worse' reads 'the common part alone scores worse at equal movement'; 'mixed' is printed with the "
    "two comparisons' labels. Statement 3 uses the same wording with 'equal projection' and is read by nothing.",
    "Statement 3's same-size multiple is per group: sqrt of the form's recorded squared size ratio for m1_r_gmom_cp "
    "and twice it for m1_r_gmom_cp_x2, the uniform arms' gain over the pooled arm of the arm's fraction by "
    "zone_read.size_alone (linear between 1, 1.5 and 2), 'outside the uniform arms' elsewhere.",
    "a_n,g, a_p,g and the (1, 1) deviation are compared between the five files exactly (one computation on the same "
    "inputs); the (1, 1) check is the largest absolute difference from d_g on M per group.",
    "Labels on each library truth: T46's single-task test of each arm against its comparison arm (ahead, behind, "
    "level).",
)
LIMITS = (
    "One released day, in-sample, on cuts earlier choices were made on.",
    "'Equal movement' is equal squared size of the direction over the moved genes, not equal counts moved or equal "
    "change in expression.",
    "c is a mean of four directions in which the endocardial cells count twice, and it holds a quarter of each group's "
    "own direction.",
    "The released day's low share along the specific parts (try T51) may be a slower pace of what separates cell types "
    "on this one day and says nothing about the two days of the test.",
    "Larger steps score better on these cuts and how far depends on the truth's depth (tries T50 and T54), so an arm "
    "that moves more than its comparison arm can gain by size alone, and no outcome sets a size at the test.",
    "The two fractions share the seeds, the draw and the directions, so agreement at both is one result read twice.",
)


# ---------------------------------------------------------------------------------------------------------------------
# The directions


def moved_masks(entry: Mapping[str, np.ndarray], genes: Sequence[str], held: Sequence[str]) -> dict[str, np.ndarray]:
    """Per group, the moved genes (zone_read.moved_mask of the entry's direction on itself)."""
    return {g: zr.moved_mask(np.asarray(entry[g], np.float64), np.asarray(entry[g], np.float64), genes, held)
            for g in sorted(entry)}


def split(entry: Mapping[str, np.ndarray], genes: Sequence[str], held: Sequence[str]) -> dict:
    """M, c and s_g as try T51 makes them; refused (ValueError) unless every group's moved genes are M."""
    masks = moved_masks(entry, genes, held)
    M = split_read.common_mask(masks)
    off = [g for g, m in masks.items() if not np.array_equal(m, M)]
    if off:
        raise ValueError(f"T55 is refused: the moved genes of {off} are not M (the genes moved in every group)")
    sv = split_read.split_vectors({g: np.asarray(entry[g], np.float64) for g in entry}, M)
    return {"M": M, "c": sv["c"], "s": sv["s"]}


def multiples(entry: Mapping[str, np.ndarray], sp: Mapping) -> dict[str, dict[str, float]]:
    """Per group a_n,g = sqrt(sum(d_g^2) / sum(c^2)) and a_p,g = sum(d_g^2) / sum(c d_g) over M; refused (ValueError)
    where sum(c d_g) is not positive."""
    M, c = sp["M"], sp["c"][sp["M"]]
    out, bad = {}, []
    for g in sorted(entry):
        d = np.asarray(entry[g], np.float64)[M]
        dd, cd, cc = float(d @ d), float(c @ d), float(c @ c)
        if not (np.isfinite(cd) and cd > 0):
            bad.append(g)
            continue
        out[g] = {N: float(np.sqrt(dd / cc)), P: dd / cd}
    if bad:
        raise ValueError(f"T55 is refused: sum(c d_g) is not positive in {bad}")
    return out


def form_direction(d: np.ndarray, c: np.ndarray, s: np.ndarray, M: np.ndarray, a: float, b: float) -> np.ndarray:
    """a c + b s on M, d's own values off M."""
    out = np.array(d, dtype=np.float64, copy=True)
    out[M] = a * c[M] + b * s[M]
    return out


def all_directions(entry: Mapping[str, np.ndarray], genes: Sequence[str], held: Sequence[str]) -> tuple[dict, dict]:
    """Every form's directions (form -> group -> array) and the task file's direction record. Raises ValueError where
    the run refuses."""
    sp = split(entry, genes, held)
    mult = multiples(entry, sp)
    M, c = sp["M"], sp["c"]
    forms: dict[str, dict[str, np.ndarray]] = {}
    for f, (a, b) in FORMS.items():
        forms[f] = {g: form_direction(entry[g], c, sp["s"][g], M, mult[g][a] if isinstance(a, str) else a, b)
                    for g in sorted(entry)}
    rec: dict = {"entry": {}, "forms": {f: list(ab) for f, ab in FORMS.items()}}
    for g in sorted(entry):
        d = np.asarray(entry[g], np.float64)
        one = form_direction(d, c, sp["s"][g], M, 1.0, 1.0)
        dev = float(np.max(np.abs(one[M] - d[M]))) if M.any() else 0.0
        rec["entry"][g] = {"sha256": zr.direction_sha256(d), N: mult[g][N], P: mult[g][P],
                           **split_read.energy_parts(d, c, sp["s"][g], M), "rebuild_11_max_abs": dev}
    for f in FORMS:
        rec[f] = {g: {**zr.compare_directions(forms[f][g], entry[g], genes, held),
                      "sha256": zr.direction_sha256(forms[f][g])} for g in sorted(entry)}
    return forms, rec


# ---------------------------------------------------------------------------------------------------------------------
# The run


def arm_fractions(fraction: float) -> dict[str, float]:
    fr = sr.arm_fractions(fraction, UNIFORM_MULTIPLIERS)
    return {**fr, **{a: fr[p] for a, (_, p) in NEW_ARMS.items()}}


def build_seed(ctx: harness.BuildContext, b0, entry: Mapping, forms: Mapping[str, Mapping], groups: Mapping,
               fraction: float, held_m1: Sequence[str], held_factors: Sequence[str]) -> dict[str, Stage]:
    """One seed's arms: b0, m1_r_gmom and the uniform arms (step_read.step_arms), then each new arm, the build of the
    pooled arm of its fraction with its form's directions in place of the entry's and every other argument the same."""
    preds = {"b0": b0(ctx), **sr.step_arms(ctx, entry, groups, fraction, held_m1, held_factors,
                                           multipliers=UNIFORM_MULTIPLIERS)}
    fr = arm_fractions(fraction)
    for arm, (key, _) in NEW_ARMS.items():
        preds[arm] = gc.m1_r_arms(ctx, forms[key], groups, fr[arm], held_m1, held_factors)[BASE]
    return preds


def cmd_run(a) -> None:
    if ROW not in n12.tries_rows():
        sys.exit(f"T55 is refused: {ROW} has no row in experiments/harness/TRIES.md (a row comes before its run)")
    n22 = RESULTS / a.n22 / "n22.json"
    if not gc.n22_passed(n22):
        sys.exit(f"T55 is refused: {gv._rel(n22)} does not hold a passing N22 (T26's route needs N22 passed)")
    out = RESULTS / a.out
    name = TASK_FILE.format(truth=a.truth)
    gv.refuse_existing(out / "tasks", (name,))
    started = gv.tree_state()
    seeds = gv.parse_seeds(a.seeds)
    data = Path(a.data)
    t = {x.key: x for x in harness.load_tasks()["T1:test"]}[TASK_OF[a.truth]]
    boards = load_boards()
    genes = boards[t.panel_board].genes
    truth = zr.load_truth(t, data, genes, a.truth)
    ref_stage = Stage.from_anndata(read_stage(data / t.reference, genes), genes)
    board = harness.task_board(t, boards, truth.n, ref_stage.n)
    stamp_inputs = [data / t.target, data / t.reference, *[data / f for f in t.inputs.values()]]
    res = {"label": LABEL_T55, "task": t.key, "truth": a.truth, "seeds": list(seeds)}
    floors, failures = lsr.all_anchors(t, board, truth, ref_stage, seeds)
    if failures:
        res.update({"anchor_failures": failures, "rows": [], "stamp": gv.stamp(started, stamp_inputs)})
        zr.write_task_file(out, name, res)
        sys.exit(f"T55 is refused on {a.truth}: the floor's anchors fail on seeds {[f['seed'] for f in failures]} "
                 "(recorded in the task file; T55 is void)")
    inputs = {r: load_stage(data / f, lineage=False) for r, f in t.inputs.items()}
    pool = inputs[harness.base_role(t)]
    pool_genes = list(pool.var_names)
    cfg = feng_dirs.load_config()
    st, src = share_read.feng_inputs()
    entry, _, _ = feng_dirs.directions_from_stats(st, pool_genes, cfg, BRACKET)
    groups = feng_dirs.label_groups(pool.obs["celltype"].astype(str).unique())
    fraction = feng_dirs.bracket_fraction(harness.STAGE_DAYS[t.reference], harness.STAGE_DAYS[t.target])
    preregistered = t1_scope.held_flat_genes()
    union = hygiene.load_union()
    held_m1 = gc.held_feng(union, preregistered)
    try:
        form_dirs, directions = all_directions(entry, pool_genes, held_m1)
    except ValueError as e:
        sys.exit(str(e))
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
    res.update({"fraction": fraction, "fractions": arm_fractions(fraction),
                "new_arms": {a_: list(v) for a_, v in NEW_ARMS.items()}, "external": external,
                "directions": directions, "rows": rows,
                "stamp": gv.stamp(started, [*stamp_inputs, t1_scope.HELD_FLAT_GENES, ts.R1_CONFIG,
                                            REPO_ROOT / union["file"], feng_dirs.clean_stats_path(), n22])})
    zr.write_task_file(out, name, res)


# ---------------------------------------------------------------------------------------------------------------------
# The reading


def direction_values(directions: Mapping) -> dict | None:
    """Every recorded sha256, a_n,g and a_p,g ((what, group) -> value); None without a full record."""
    try:
        out = {}
        for g, v in directions["entry"].items():
            out[("entry", g)] = v["sha256"]
            out[(N, g)] = v[N]
            out[(P, g)] = v[P]
        for f in FORMS:
            out.update({(f, g): v["sha256"] for g, v in zr._groups(directions[f]).items()})
        return out
    except (KeyError, TypeError, AttributeError):
        return None


def rebuild_off(directions: Mapping) -> list[str] | None:
    """The groups whose recorded (1, 1) deviation is missing or above REBUILD_TOL; None without a record."""
    try:
        ent = directions["entry"]
        if not ent:
            return None
        return sorted(g for g, v in ent.items()
                      if not isinstance(v.get("rebuild_11_max_abs"), int | float)
                      or not v["rebuild_11_max_abs"] <= REBUILD_TOL)
    except (KeyError, TypeError, AttributeError):
        return None


def validity(results: Mapping[str, Mapping], committed: Mapping[str, Sequence[Mapping]]) -> dict:
    """The void conditions. results: truth -> task file; committed: truth -> the committed file's rows."""
    reasons = []
    missing = [k for k in TRUTHS if k not in results]
    if missing:
        reasons.append(f"missing truths {missing}")
    want_fr = arm_fractions(FRACTION)
    values = {}
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
        rep = zr.reproduction(rows, committed.get(k))
        if not rep["holds"]:
            reasons.append(f"{k}: {', '.join(REPRO_ARMS)} do not reproduce their committed scores within {REPRO_TOL:g} "
                           "per seed")
        if rows and not absent:
            same = [a for a, (_, p) in NEW_ARMS.items()
                    if all(row["arms"][a]["score"] == row["arms"][p]["score"] for row in rows)]
            if same:
                reasons.append(f"{k}: new arms score exactly as the pooled arm of their fraction on every seed: {same}")
        directions = r.get("directions") or {}
        values[k] = direction_values(directions)
        if values[k] is None:
            reasons.append(f"{k}: no full direction record")
        off = rebuild_off(directions)
        if off is None:
            reasons.append(f"{k}: no record of the form (1, 1)")
        elif off:
            reasons.append(f"{k}: the form (1, 1) is not the entry's direction within {REBUILD_TOL:g} in {off}")
    known = {k: v for k, v in values.items() if v is not None}
    if len({json.dumps(sorted((f"{a}|{b}", x) for (a, b), x in v.items())) for v in known.values()}) > 1:
        reasons.append(f"the recorded directions' sha256 or multiples differ between the truths' files {sorted(known)}")
    commits = {(r.get("stamp") or {}).get("commit") for r in results.values()}
    unclean = sorted(k for k, r in results.items() if (r.get("stamp") or {}).get("tree_clean") is not True)
    if unclean:
        reasons.append(f"tree not clean on {unclean}")
    if len(commits) > 1:
        reasons.append(f"the truths ran at {len(commits)} commits")
    return {"valid": not reasons, "reasons": reasons, "commits": sorted(c for c in commits if c)}


def pair_statement(tests: Sequence[str], arms: Sequence[str], better: str, worse: str) -> str:
    """Both pass: better; both lose: worse; both neither: no difference at the bar; otherwise mixed with both labels."""
    if all(t == PASSES for t in tests):
        return better
    if all(t == LOSES for t in tests):
        return worse
    if all(t == NEITHER for t in tests):
        return NO_DIFFERENCE
    return f"{MIXED} (" + "; ".join(f"{a}: {t}" for a, t in zip(arms, tests, strict=True)) + ")"


def same_size(means: Mapping[str, float], size_ratio: float | None, pooled: str) -> dict:
    """The uniform arms' figure at the multiple with the same squared size as the arm: sqrt(size_ratio) times the
    pooled arm's multiple (1 or 2), linear between 1, 1.5 and 2 (zone_read.size_alone), as a gain over the pooled arm;
    None outside, with the reason."""
    m = (UNIFORM[pooled] * float(np.sqrt(size_ratio))
         if isinstance(size_ratio, int | float) and np.isfinite(size_ratio) and size_ratio >= 0 else None)
    return {"multiple": m, "what_size_alone_would_give": zr.size_alone(means, m, pooled),
            "why_not": zr.multiple_reason(m, "the same-size multiple")}


def comparison_pairs() -> dict[str, tuple[str, str]]:
    """Every comparison the read reports: name -> (arm, base arm)."""
    return {a: (a, p) for a, (_, p) in NEW_ARMS.items()}


def t55_reading(results: Mapping[str, Mapping], committed: Mapping[str, Sequence[Mapping]]) -> dict:
    v = validity(results, committed)
    out = {"valid": v["valid"], "validity": v}
    if not v["valid"]:
        out["outcome"] = "VOID"
        return out
    means = {k: {a: float(np.mean([row["arms"][a]["score"] for row in results[k]["rows"]])) for a in ARMS}
             for k in TRUTHS}
    pairs = comparison_pairs()
    cut_tests = {n: zr.cut_test(results, a, b) for n, (a, b) in pairs.items()}
    comps = {k: {n: lsr.truth_comparison(results[k], k, a, base=b) for n, (a, b) in pairs.items()} for k in TRUTHS}
    labels = {k: {n: zr.library_label(comps[k][n]) for n in pairs} for k in LIBS}
    rec = results[lsr.POOLED]["directions"]
    s1 = {a: {"statement": S1_LABELS[cut_tests[a]["test"]], "comparison": cut_tests[a],
              "library_labels": {k: labels[k][a] for k in LIBS}} for a in S1_ARMS}
    s2 = {"statement": pair_statement([cut_tests[a]["test"] for a in S2_ARMS], S2_ARMS, BETTER_MOVE, WORSE_MOVE),
          "comparisons": {a: cut_tests[a] for a in S2_ARMS},
          "library_labels": {k: {a: labels[k][a] for a in S2_ARMS} for k in LIBS}}
    s3 = {"statement": pair_statement([cut_tests[a]["test"] for a in S3_ARMS], S3_ARMS, BETTER_PROJ, WORSE_PROJ),
          "read_by_nothing": True, "comparisons": {a: cut_tests[a] for a in S3_ARMS},
          "library_labels": {k: {a: labels[k][a] for a in S3_ARMS} for k in LIBS},
          "same_size": {a: {k: {g: same_size(means[k], x.get("size_ratio"), NEW_ARMS[a][1])
                                for g, x in sorted(zr._groups(rec[NEW_ARMS[a][0]]).items())}
                            for k in TRUTHS} for a in S3_ARMS}}
    out.update({"outcome": {"statement_1": {a: s["statement"] for a, s in s1.items()},
                            "statement_2": s2["statement"], "statement_3": s3["statement"]},
                "statement_1": s1, "statement_2": s2, "statement_3": s3, "mean_scores": means,
                "comparisons": comps, "cut_tests": cut_tests, "directions": rec})
    return out


def t55_text(res: Mapping) -> str:
    """The reading as text; no cell count is printed."""
    r = res["reading"]
    lines = [res["label"], "Choices the brief left open:"]
    lines.extend(f"  - {x}" for x in CHOICES)
    if not r["valid"]:
        lines.append("Outcome: VOID")
        lines.append("Void: " + "; ".join(r["validity"]["reasons"]))
    else:
        lines.append("Statement 1, the group-specific parts at the entry's size:")
        for a, s in r["statement_1"].items():
            lines.append(f"  {a}: {s['statement']}")
            lines.append("    " + zr._cut_line(s["comparison"]) + "; library truths "
                         + ", ".join(f"{k} {lab}" for k, lab in s["library_labels"].items()))
        for name, key in (("Statement 2, content at equal movement", "statement_2"),
                          ("Statement 3, content at equal projection (read by nothing)", "statement_3")):
            s = r[key]
            lines.append(f"{name}: {s['statement']}")
            for a, c in s["comparisons"].items():
                lines.append("    " + zr._cut_line(c))
                if key == "statement_3":
                    lines.append("      what the uniform arms give at the same squared size (read by nothing): "
                                 + "; ".join(f"{k} " + ", ".join(
                                     f"{g} x{zr._f(v['multiple'], '.3f')} -> "
                                     + (zr._f(v["what_size_alone_would_give"]) if v["why_not"] is None
                                        else v["why_not"])
                                     for g, v in per.items()) for k, per in s["same_size"][a].items()))
            lines.append("    library truths (single-task test against the pooled arm of the same fraction): "
                         + "; ".join(f"{k}: " + ", ".join(f"{a} {lab}" for a, lab in ls_.items())
                                     for k, ls_ in s["library_labels"].items()))
        d = r["directions"]
        lines.append("Per group (over M; ratios only): " + "; ".join(
            f"{g}: a_n {v[N]:.4f}, a_p {v[P]:.4f}, sum(c^2)/sum(d^2) {v['e_c']:.4f}, sum(s^2)/sum(d^2) {v['e_s']:.4f}, "
            f"2 sum(c s)/sum(d^2) {v['cross']:+.4f}, (1, 1) off by {v['rebuild_11_max_abs']:.2e}"
            for g, v in d["entry"].items()))
        for f in FORMS:
            lines.append(f"form {f} {tuple(FORMS[f])} on the entry's direction: " + "; ".join(
                f"{g} slope {zr._f(v['slope'], '.4f')}, squared size ratio {zr._f(v['size_ratio'], '.4f')}"
                for g, v in zr._groups(d[f]).items()))
        for k in TRUTHS:
            lines.append(f"truth {k}: mean scores " + ", ".join(f"{a} {x:.2f}" for a, x in r["mean_scores"][k].items()))
            for c in r["comparisons"][k].values():
                lines.append("    " + lsr._comparison_line(c))
    lines.extend(f"Limit: {x}" for x in LIMITS)
    lines.append("Limit: in-sample labelled diagnostic; never adoption evidence; puts nothing forward by itself; it "
                 "can only go to the team captain")
    lines.append(f"Commit {res['stamp']['commit']}, tree clean {res['stamp']['tree_clean']}, UTC {res['stamp']['utc']}")
    return "\n".join(lines) + "\n"


def read(out: Path, committed_files: Mapping[str, Path] | None = None) -> dict:
    """Read the truths present in out/tasks against the committed T46 and T44 files; write t55.{json,txt}."""
    gv.refuse_existing(out, (f"{STEM}.json", f"{STEM}.txt"))
    started = gv.tree_state()
    files = COMMITTED if committed_files is None else committed_files
    present = {k: out / "tasks" / TASK_FILE.format(truth=k) for k in TRUTHS
               if (out / "tasks" / TASK_FILE.format(truth=k)).exists()}
    results = {k: json.loads(p.read_text()) for k, p in present.items()}
    found = {k: Path(p) for k, p in files.items() if Path(p).exists()}
    committed = {k: json.loads(p.read_text())["rows"] for k, p in found.items()}
    res = {"label": LABEL_T55, "choices": list(CHOICES), "limits": list(LIMITS),
           "reading": t55_reading(results, committed), "stamps": {k: v.get("stamp") for k, v in results.items()}}
    res["stamp"] = gv.stamp(started, [*present.values(), *found.values()])
    gv.write_outputs(out, STEM, res, t55_text(res))
    return res


def cmd_read(a) -> None:
    res = read(RESULTS / a.out)
    if res["reading"]["outcome"] == "VOID":
        sys.exit("T55's reading is void and goes to the captain")


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("run", help="one truth of try T55")
    p.add_argument("--truth", required=True, choices=list(TRUTHS))
    p.add_argument("--out", required=True)
    p.add_argument("--n22", default="n22_additivity")
    p.add_argument("--seeds", default=",".join(map(str, SEEDS)))
    p.add_argument("--data", default=str(harness.DATA_DIR))
    p = sub.add_parser("read", help="try T55's reading")
    p.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    {"run": cmd_run, "read": cmd_read}[a.cmd](a)


if __name__ == "__main__":
    main()
