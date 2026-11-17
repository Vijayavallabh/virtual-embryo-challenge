#!/usr/bin/env python
"""Try T57 (experiments/harness/TRIES.md): the common-part form (cn) and the form that keeps the group-specific parts
(dc), against the entry's direction at 1, 2, 3 and 4 times the practice step, in worlds made from libraries 3 and 2. A
labelled in-sample read on the routes of tries T50 (experiments/t1_m1/over_read.py), T54
(experiments/t1_m1/deep_read.py), T55 (experiments/t1_m1/part_read.py) and T56 (experiments/t1_m1/part_depth_read.py).
Built to brief AL.

Try T55 found that M1's step along the common part c of Feng's direction alone, at the entry's squared size (cn: per
group a_n x c), scores above the entry's step on the released day, and try T56 that it stays ahead where the truth is
shallower or deeper than the reference, at one and two times the step. A cold check found that cn raises the
coefficient on c while deleting the part of each group's direction orthogonal to c, which the released day follows.
This read adds the arm that keeps that orthogonal part and only raises the coefficient on c (dc), worlds built from
library 2, and every form at 1, 2, 3 and 4 times the step.

Forms, per group, on the common mask M (try T55's: part_read.split and part_read.multiples), from the directions alone:
  cn   a_n x c on M, the entry's own values off M (try T55's form, part_read.all_directions, unchanged).
  dc   a_n x c + o_g on M, where p_g = sum(c d_g) / sum(c c) is the entry's own coefficient on c and o_g = d_g - p_g x c
       its part orthogonal to c; written d_g + (a_n - p_g) x c, so that a_n = p_g gives d_g bit for bit; the entry's
       own values off M. The run refuses if a_n is below p_g in a group, and records per group p_g, a_n, the form's
       squared size over the entry's (zone_read.compare_directions) and its sha256.

  run    One world of the chamber cut T1/E9.5->E10.5:chamber, 10 seeds. Worlds (six):
           lib3        the truth library 3 as try T46 makes it (over_read.truth_and_contamination), the inputs as
                       they are;
           lib2        the truth library 2 as try T46 makes it, the inputs as they are;
           ref67       try T54's worlds of those names (deep_read.world_reference: the reference thinned by try T54's
           ref50       own code, p and stream, before anything else is made), the truth library 3;
           ref67_lib2  the same thinned references (deep_read.world_reference called under the name ref67 or ref50, so
           ref50_lib2  the same p and stream), the truth library 2.
         Arms per seed, all from one cell draw: b0; the entry's arms m1_r_gmom, m1_r_gmom_x2, _x3, _x4 and, as the
         plain-step control of statement 3 (added after review), _x5, _x6, _x8 (try T50's names and builds:
         step_read.step_arms on feng_dirs.feng_directions with over_read.MULTIPLIERS' values); cn and dc at the first
         four fractions (m1_r_gmom_cn, _cn_x2, _cn_x3, _cn_x4; m1_r_gmom_dc, _dc_x2, _dc_x3, _dc_x4), each the build of
         the entry's arm of its fraction (gmom_combo.m1_r_arms) with the form's directions in place of the entry's and
         every other argument the same. Refused on an existing task file, unless the committed N22 passed, unless
         experiments/harness/TRIES.md holds row T57, and where the forms refuse.
  read   T57's reading of the six worlds' files (below, fixed here and repeated in the ledger row before the run).
         Writes t57.json and t57.txt and refuses to overwrite them.

The reading (brief AL). Void unless: the six worlds hold the ten seeds with every arm and its fraction at one commit
on a clean tree; no file records an anchor failure; where a committed score exists it is reproduced, seed for seed and
to 1e-9: b0 and the entry's arms in lib3 and lib2 (x1 to x8: try T50's lib3.json and lib2.json) and in ref67 and ref50
(x1 to x6: try T54's files); m1_r_gmom_cn and _cn_x2 in lib3, ref67 and ref50 (try T56's files) and in lib2 (try
T55's lib2.json); the thinned worlds pass try T54's checks (their p and stream, depth ratio within 0.01 of p, no cell
emptied, the rebuild check); in ref67_lib2 and ref50_lib2 the thinning record's p and stream equal ref67's and
ref50's; per group the entry's directions as tries T50 and T54 load them and as try T55 loads them have the same
sha256; the forms' recorded directions and multiples are the same in the six files, and try T55's record (cn's among
it) equals try T55's lib3 file; in every file dc's a_n per group equals the entry's recorded a_n, p_g is at most a_n
and dc's squared size over the entry's is at least 1 (added after review); every file records the truth name the route
passed to try T46's truth maker and the truth matrix's sha256, and the worlds of one library have the same truth
sha256 (added after review); no form's arm scores exactly as the entry's arm of its fraction on every seed of a world.
A cell is a world and a multiple (1, 2, 3, 4): 24 cells. Per form and cell, the paired difference of the form's arm from
the entry's arm of the same fraction, with its mean, seeds ahead, per-seed range and its label by try T46's single-task
test ("ahead": a mean paired gain above 0.3 with at least 8 of 10 seeds ahead; "behind" the same reversed; else
"level").
Statement 1, per form: cn "ahead of the entry's direction in every cell read", dc "the entry plus a further step along c
is ahead of the entry in every cell read" (changed after review: dc minus the entry is exactly a further step along c of
m x fraction x (a_n - p_g) per group, so its sign reads depth along c, not the form), if in all 24 cells the difference
has a positive mean with at least 8 of 10 seeds ahead; "behind somewhere" if in some cell it is behind by try T46's
test, naming the cells; otherwise "not separated", naming every cell that meets neither condition. Each statement 1
carries its size (added after review): the number of the 24 cells ahead by try T46's test and the smallest and largest
cell means.
Statement 2, per cell: dc minus cn with its label by try T46's test; the summary counts the cells where dc is ahead,
behind and level, and names them.
Statement 3, beyond a larger step, per world and form: the form's highest mean over its four multiples against the
entry's highest mean over 1, 2, 3, 4, 5, 6 and 8 times (changed after review from the entry's four: the forms move along
c as far as the entry at m x a_n/p_g, up to about 5.9 at x4), paired by seed at those two multiples: mean, seeds ahead,
and "above the best plain step read" where the mean is positive with at least 8 of 10 seeds ahead (otherwise "not
above"); "not readable" where the entry's best is at the edge of its range (x8) and the multiple at which the entry
would move along c as far as the form at its best (the form's multiple x the largest a_n/p_g over the groups) lies
beyond 8.
Reported and read by nothing: every arm's mean score per world (x5 to x8 included); each difference's parts on the
two de terms and on mmd_u and the variogram, and its raw-metric differences; the forms' per-group records.

NOT an adoption result: a labelled in-sample diagnostic, as tries T44 to T56 are; the result can only go to the team
captain, and no outcome is evidence on Task 1's read depth, which stays closed. The limits are LIMITS below, repeated in
the text.

Usage: on the DGX campaign checkout, repo root, all six worlds at one commit, committing nothing until the six task
files exist, and on the machine tries T50, T54, T55 and T56 ran on, so that the rescored arms reproduce. The environment
variables are docs/COMPUTE.md's for repository Python on the DGX.
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \\
    uv run --frozen --no-sync python -I -B experiments/t1_m1/part_range_read.py run --world lib3 --out t57_range
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \\
    uv run --frozen --no-sync python -I -B experiments/t1_m1/part_range_read.py run --world lib2 --out t57_range
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \\
    uv run --frozen --no-sync python -I -B experiments/t1_m1/part_range_read.py run --world ref67 --out t57_range
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \\
    uv run --frozen --no-sync python -I -B experiments/t1_m1/part_range_read.py run --world ref50 --out t57_range
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \\
    uv run --frozen --no-sync python -I -B experiments/t1_m1/part_range_read.py run --world ref67_lib2 --out t57_range
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \\
    uv run --frozen --no-sync python -I -B experiments/t1_m1/part_range_read.py run --world ref50_lib2 --out t57_range
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \\
    uv run --frozen --no-sync python -I -B experiments/t1_m1/part_range_read.py read --out t57_range
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from collections.abc import Callable, Mapping, Sequence
from pathlib import Path

import anndata as ad
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
import deep_read as dr  # noqa: E402
import over_read as ob  # noqa: E402
import part_depth_read as pdr  # noqa: E402
import part_read as pr  # noqa: E402

ls, sr, gc, er, zr = ob.ls, ob.sr, ob.gc, ob.er, pr.zr
gv, hygiene, n12, feng_dirs = ob.gv, ob.hygiene, ob.n12, ob.feng_dirs
TASK = ob.TASK  # T1/E9.5->E10.5:chamber
LIB3, LIB2 = ob.HEART, "lib2"
ROW = "T57"
STEM = "t57"
WORLDS = (LIB3, LIB2, "ref67", "ref50", "ref67_lib2", "ref50_lib2")
REFERENCE_OF = {"ref67": "ref67", "ref50": "ref50", "ref67_lib2": "ref67", "ref50_lib2": "ref50"}  # try T54's world
TWIN = {"ref67_lib2": "ref67", "ref50_lib2": "ref50"}  # a library-2 world and its library-3 twin
TRUTH_OF = {LIB3: LIB3, LIB2: LIB2, "ref67": LIB3, "ref50": LIB3, "ref67_lib2": LIB2, "ref50_lib2": LIB2}
THINNING_OF = {k: dr.THINNED[v] for k, v in REFERENCE_OF.items()}
BASE, X2, X3, X4 = sr.BASE, ob.X2, ob.X3, ob.X4  # m1_r_gmom (0.2), _x2 (0.4), _x3 (0.6), _x4 (0.8)
X5, X6, X8 = ob.X5, ob.X6, ob.X8  # try T50's arms at 1.0, 1.2 and 1.6: statement 3's plain-step control only
MULTIPLES = (1, 2, 3, 4)  # the cells' multiples and the forms' arms
ENTRY_AT = {1: BASE, 2: X2, 3: X3, 4: X4}  # multiple -> the entry's arm of a cell
ENTRY_RANGE = {**ENTRY_AT, 5: X5, 6: X6, 8: X8}  # multiple -> the entry's arm, statement 3's plain-step range
ENTRY_EDGE = max(ENTRY_RANGE)
ENTRY_MULTIPLIERS = {a: float(ob.MULTIPLIERS[a]) for a in (X2, X3, X4, X5, X6, X8)}  # try T50's values
CN, DC = "cn", "dc"
FORMS = (CN, DC)


def form_arm(form: str, m: int) -> str:
    """The form's arm at the multiple m: m1_r_gmom_cn, m1_r_gmom_cn_x2, ..."""
    return f"{BASE}_{form}" + ("" if m == 1 else f"_x{m}")


FORM_ARMS = {form_arm(f, m): (f, ENTRY_AT[m]) for f in FORMS for m in MULTIPLES}  # arm -> (form, entry arm)
ENTRY_ARMS = tuple(ENTRY_RANGE.values())
ARMS = ("b0", *ENTRY_ARMS, *FORM_ARMS)
REPRO_ENTRY = {"t50": ("b0", *ENTRY_ARMS), "t54": ("b0", BASE, X2, X3, X4, X5, X6)}  # try T54 has no x8
REPRO_CN = (form_arm(CN, 1), form_arm(CN, 2))
CELLS = tuple((k, m) for k in WORLDS for m in MULTIPLES)
FRACTION, FRACTION_TOL, REPRO_TOL = sr.FRACTION, sr.FRACTION_TOL, sr.REPRO_TOL  # 0.2, its tolerance, 1e-9
DEPTH_TOL, REBUILD_TOL = ob.DEPTH_TOL, ob.REBUILD_TOL  # 0.01, 1e-5 (try T54's)
HOLD_AHEAD = 8  # seeds ahead of 10
T50_TASKS = RESULTS / "t50_m1_over_step" / "tasks"
T54_TASKS = RESULTS / "t54_m1_deep_step" / "tasks"
T55_TASKS = RESULTS / "t55_m1_part" / "tasks"
T56_TASKS = RESULTS / "t56_cn_depth" / "tasks"
ENTRY_SOURCE = {LIB3: "t50", LIB2: "t50", "ref67": "t54", "ref50": "t54"}  # the entry's arms' committed scores
CN_SOURCE = {LIB3: "t56", "ref67": "t56", "ref50": "t56", LIB2: "t55"}  # cn's and cn_x2's committed scores
TASK_FILE = "{world}.json"
AHEAD_ALL = "ahead of the entry's direction in every cell read"
AHEAD_ALL_DC = "the entry plus a further step along c is ahead of the entry in every cell read"
AHEAD_ALL_OF = {"cn": AHEAD_ALL, "dc": AHEAD_ALL_DC}
QUANTITY = {"cn": "cn minus the entry's arm of its fraction: the deleted orthogonal parts and the raised coefficient "
                  "on c together",
            "dc": "dc minus the entry's arm of its fraction: exactly the entry plus a further step along c of m x "
                  "fraction x (a_n - p_g) per group, a reading of depth along c, not of the form"}
BEHIND_SOMEWHERE = "behind somewhere"
NOT_SEPARATED = "not separated"
ABOVE = "above the best plain step read"
NOT_ABOVE = "not above"
NOT_READABLE = "not readable"
LABEL_T57 = ("try T57, a labelled in-sample diagnostic, never an adoption result: M1 (with R1, R2 and g_mom's factors, "
             "the entry's form) with its Feng step along the common part of Feng's direction alone at the entry's "
             "squared size (cn) and along the entry's direction with its coefficient on the common part raised to "
             "that size's (dc), at 1, 2, 3 and 4 times the step, against the entry's step of the same fraction, with "
             "library 3's and library 2's truths and the reference as released and thinned to two thirds and a half "
             "(try T54's draws) on the chamber cut, with the entry's step also at 5, 6 and 8 times as statement 3's "
             "plain-step control; the result can only go to the team captain, and no outcome is "
             "evidence on Task 1's read depth, which stays closed")
LIMITS = (
    "One released day, in-sample, the chamber cut, two libraries.",
    "Each difference is taken inside one world, in that world's own points, and is not compared across worlds.",
    "A thinned reference is one draw for all seeds and also changes the floor, the de terms' reference, B0's shift and "
    "g_mom's fit.",
    "'In every cell read' covers 1 to 4 times the practice step and the six worlds and nothing else.",
    "Both forms move further along c than the entry does (at the multiple m, as far as the entry at m x a_n/p_g, up to "
    "about 1.5 m), so a difference in a cell where larger plain steps still gain can be that larger movement and not "
    "the form. Statement 3 compares a form's best over 1 to 4 times with the entry's best over 1, 2, 3, 4, 5, 6 and 8 "
    "times, a range that holds the multiple at which the entry moves along c as far as the form; it is ruled not "
    "readable where that multiple lies beyond 8 while the entry's best is at 8. Within that range a larger plain step "
    "cannot give its 'above'; beyond 8 times and between the multiples read, nothing here says.",
    "cn deletes parts of the groups' directions that the released day follows, and what those parts carry over the two "
    "days of the test no released day shows. A cn lead needs a truth that follows those parts less than the released "
    "day does; no world here varies how much of them the truth carries (all six share one released reference and one "
    "direction set), so cn's statement 1 cannot speak to a test day that follows them more or less.",
    "Nothing here says how deep any hidden truth is, and no statement may be used with a figure for it.",
    "On M the directions obey two identities: dc - cn = o_g exactly (the same coefficient a_n on c), so statement 2 is "
    "the one clean reading of the deleted orthogonal parts at a fixed coefficient on c; and "
    "dc - entry = (a_n - p_g) x c exactly, so dc against the entry at the multiple m is the entry plus a further "
    "step along c of m x fraction x (a_n - p_g): a reading of depth along c, not of the form. dc's squared size over "
    "the entry's is 2 - (p_g/a_n)^2 (recorded per group).",
    "Statement 1's bar (a positive mean with at least 8 of 10 seeds ahead) has no size threshold: the arms share one "
    "cell draw, so a consistent tiny shift passes it in every cell while every cell is 'level' by try T46's test; the "
    "number of cells ahead by that test and the smallest and largest cell means are printed with it and go with it.",
    "A world split (a form ahead where plain steps still gain and behind where they have turned) is evidence of step "
    "size, the story of tries T44 to T50, not of the form.",
    "The library-2 thinned worlds have no draw control: try T54's second draw (ref67b) is not repeated with library "
    "2's truth, and the two thinned draws are shared between the libraries.",
    "The 24 cells share the seeds, the draw and the directions within a world, and the thinned worlds share their two "
    "draws between the libraries: agreement across cells is not 24 confirmations.",
    "Statement 3 picks each side's best multiple on the same ten seeds it then compares: both maxima are chosen in "
    "sample, and only the multiples read are searched (the forms' four, the entry's seven).",
    "6 of the 24 cells (lib3, ref67 and ref50 at 1 and 2 times the step) repeat tries T50, T54 and T56 by rule and add "
    "no evidence for cn there; the new evidence is dc, the multiples 3 and 4, and the library-2 worlds.",
)
CHOICES = (
    "dc is computed as d_g + (a_n - p_g) x c on M, algebraically a_n x c + o_g, so that where a_n equals p_g the form "
    "is the entry's direction bit for bit; off M it keeps the entry's own values (NaN or a held gene), as try T55's "
    "forms do; a_n is try T55's (part_read.multiples) and M, c try T55's (part_read.split).",
    "dc's record per group: p_g, a_n, sha256, and zone_read.compare_directions on the entry (slope, pearson, the "
    "squared size over the entry's over the moved genes, which are M).",
    "The directions are loaded by composing the calls of try T56 (feng_dirs.feng_directions for the entry's arms; "
    "share_read.feng_inputs and feng_dirs.directions_from_stats for the forms, as try T55 loads them), with try T56's "
    "record of whether the two routes agree by sha256 per group; the read voids where they do not (try T56's "
    "amendment, kept so that 'differs in the directions only' holds).",
    "The reproduction of the entry's arms includes b0 (as tries T56 and T54 reproduce it): b0 and the four entry arms "
    "against try T50's lib3 and lib2 files and try T54's ref67 and ref50 files.",
    "cn's committed scores in lib3 are taken from try T56's lib3 file (as the brief says), not try T55's (which T56's "
    "own read already matched).",
    "ref67_lib2 and ref50_lib2 call deep_read.world_reference under the name ref67 or ref50, so their thinning record "
    "names that world; the read checks p and stream against both try T54's table and the twin file's record, and "
    "applies try T54's depth, emptied-cell and rebuild checks to them.",
    "No check that lib2-truth worlds have the same predictions as their lib3 twins: the board's cell count follows "
    "the truth's size (harness.task_board), so predictions can differ between truths; nothing is added to the void "
    "list beyond the brief's but b0 and the agreement record.",
    "Statement 1 counts 'ahead' as a mean above zero and seeds strictly above zero (at least 8), as tries T46 to T56 "
    "count them; 'behind somewhere' wins over 'not separated' when both apply, and the cells meeting neither are still "
    "recorded and printed.",
    "Statement 3's best multiple is the one with the highest mean score (the smallest on a tie), among the four for "
    "the form and among the seven for the entry, inside each world; the comparison is try T50's (over_read.comparison) "
    "of those two arms paired by seed.",
    "A cell is named '<world> x<multiple>'.",
    "dc's refusal where a_n is below p_g cannot fire: a_n = sqrt(sum d_g^2 / sum c^2) >= |sum c d_g| / sum c^2 = |p_g| "
    "by Cauchy-Schwarz, with equality only where d_g is parallel to c on M; it is kept only as the brief's guard.",
    "Statement 3's c-matched multiple is the form's best multiple times the largest a_n/p_g over the groups (the "
    "farthest a group moves), from the dc record of the world's file (infinite where p_g is not positive); the entry's "
    "range adds try T50's x5, x6 and x8 (x8 is not in try T54's files, so it is not reproduced in ref67 and ref50).",
    "The truth's sha256 is that of its matrix as the scorer reads it (lib_step_read.prediction_sha256 on the truth "
    "Stage); the worlds of one library must agree on it, the thinned worlds thinning only the reference.",
    "The run's --out folder for the command lines is t57_range.",
)


# ---------------------------------------------------------------------------------------------------------------------
# The forms


def dc_direction(d: np.ndarray, c: np.ndarray, M: np.ndarray, a_n: float) -> tuple[np.ndarray, float]:
    """The dc form of one group's direction and p_g: on M, d + (a_n - p_g) c = a_n c + o_g with p_g = sum(c d)/sum(c c)
    and o_g = d - p_g c; off M, d's own values."""
    d = np.asarray(d, dtype=np.float64)
    cM, dM = np.asarray(c, dtype=np.float64)[M], d[M]
    p = float(cM @ dM) / float(cM @ cM)
    out = np.array(d, dtype=np.float64, copy=True)
    out[M] = dM + (a_n - p) * cM
    return out, p


def dc_forms(entry: Mapping[str, np.ndarray], genes: Sequence[str], held: Sequence[str]) -> tuple[dict, dict]:
    """dc per group (group -> array) and its record; M, c and a_n by try T55's functions. Raises ValueError where try
    T55's split or multiples refuse, and where a_n is below p_g in a group."""
    sp = pr.split(entry, genes, held)
    mult = pr.multiples(entry, sp)
    forms, rec, low = {}, {}, []
    for g in sorted(entry):
        d = np.asarray(entry[g], np.float64)
        a_n = mult[g][pr.N]
        f, p = dc_direction(d, sp["c"], sp["M"], a_n)
        if a_n < p:
            low.append(g)
        forms[g] = f
        rec[g] = {"p_g": p, "a_n": a_n, "sha256": zr.direction_sha256(f), **zr.compare_directions(f, d, genes, held)}
    if low:
        raise ValueError(f"T57 is refused: a_n is below p_g in {low} (dc would lower the coefficient on c there)")
    return forms, rec


def world_directions(pool_genes: Sequence[str], held_m1: Sequence[str]) -> tuple[dict, dict, dict, dict, dict]:
    """The entry's directions as tries T50 and T54 load them; the forms (cn by try T55's part_read.all_directions, dc
    by dc_forms) on the entry's directions as try T55 loads them; try T55's record, dc's record, and the external record
    (try T50's, try T55's source, and per group whether the two routes' entry directions agree by sha256). Depends on
    the pool's gene names and the held list only. Raises ValueError where a form refuses."""
    pool_genes = list(pool_genes)
    directions, ext = feng_dirs.feng_directions(pool_genes)
    cfg = feng_dirs.load_config()
    st, src = pr.share_read.feng_inputs()
    entry, _, _ = feng_dirs.directions_from_stats(st, pool_genes, cfg, pr.BRACKET)
    t55_forms, rec = pr.all_directions(entry, pool_genes, held_m1)
    dc, dc_rec = dc_forms(entry, pool_genes, held_m1)
    agree = {g: g in directions and g in entry
             and zr.direction_sha256(directions[g]) == zr.direction_sha256(entry[g])
             for g in sorted({*directions, *entry})}
    return directions, {CN: t55_forms[CN], DC: dc}, rec, dc_rec, {**ext, "t55_source": src,
                                                                  "entry_agrees_with_t55s": agree}


# ---------------------------------------------------------------------------------------------------------------------
# The worlds and the run


def world_route(world: str, t, data: Path, genes: Sequence[str], inputs: Mapping[str, ad.AnnData], role: str,
                ref_stage: Stage, read_full: Callable[[], ad.AnnData]
                ) -> tuple[Mapping[str, ad.AnnData], Stage, Stage, int, dict | None, str]:
    """The world's inputs, scorer's reference stage, truth, contamination count, thinning record and the truth name
    passed to try T46's truth maker. A thinned world: the reference thinned first by deep_read.world_reference under
    try T54's world name (ref67 or ref50); then the truth (library 3 or 2) by over_read.truth_and_contamination, which
    is try T46's for these names."""
    if world not in WORLDS:
        raise ValueError(f"unknown world {world!r}; one of {list(WORLDS)}")
    thinning = None
    if world in REFERENCE_OF:
        inputs, ref_stage, thinning = dr.world_reference(REFERENCE_OF[world], inputs, role, ref_stage, genes,
                                                         read_full)
    truth_name = TRUTH_OF[world]
    truth, contamination, _ = ob.truth_and_contamination(t, data, genes, truth_name)
    return inputs, ref_stage, truth, contamination, thinning, truth_name


def arm_fractions(fraction: float) -> dict[str, float]:
    fr = sr.arm_fractions(fraction, ENTRY_MULTIPLIERS)
    return {**fr, **{a: fr[p] for a, (_, p) in FORM_ARMS.items()}}


def build_seed(ctx: harness.BuildContext, b0, directions: Mapping, forms: Mapping[str, Mapping], groups: Mapping,
               fraction: float, held_m1: Sequence[str], held_factors: Sequence[str]) -> dict[str, Stage]:
    """One seed's arms: b0 and the entry's seven arms (x1 to x6 and x8) as tries T50 and T54 build them
    (step_read.step_arms with try T50's multipliers), then each
    form's arm, the build of the entry's arm of its fraction (gmom_combo.m1_r_arms) with the form's directions in place
    of the entry's and every other argument the same."""
    preds = {"b0": b0(ctx), **sr.step_arms(ctx, directions, groups, fraction, held_m1, held_factors,
                                           multipliers=ENTRY_MULTIPLIERS)}
    fr = arm_fractions(fraction)
    for arm, (form, _) in FORM_ARMS.items():
        preds[arm] = gc.m1_r_arms(ctx, forms[form], groups, fr[arm], held_m1, held_factors)[BASE]
    return preds


def seed_rows(view, cboard, board, truth: Stage, ref_stage: Stage, inputs: Mapping[str, ad.AnnData], floors: Mapping,
              seeds: Sequence[int], directions: Mapping, forms: Mapping, groups: Mapping, fraction: float,
              held_m1: Sequence[str], preregistered: Sequence[str], label: str = "") -> list[dict]:
    """Per seed: the arms built from the world's inputs and scored against the world's truth with the world's reference
    stage and that seed's anchors, recorded as tries T50, T54 and T56 record them."""
    b0 = harness._b0_t1(harness.B0_DAMPING)
    rows = []
    for s in seeds:
        t0 = time.time()
        f = floors[s]
        ctx = harness.BuildContext(task=view, board=cboard, inputs=inputs, seed=harness.BUILD_SEED_OFFSET + s,
                                   derived={"t1_held_flat": preregistered})
        preds = build_seed(ctx, b0, directions, forms, groups, fraction, held_m1, preregistered)
        arms = er.score_arms(board, truth, ref_stage, s, f.anchors, preds)
        rows.append({"seed": s, "floor": float(f.score), "anchors": ob.anchor_record(f.anchors), "arms": arms,
                     "sha256": {k: ls.prediction_sha256(p) for k, p in preds.items()},
                     "sparsity": {k: ls.sparsity(p) for k, p in preds.items()}, "seconds": round(time.time() - t0, 1)})
        print(f"{label} seed {s}: " + ", ".join(f"{k} {v['score']:.2f}" for k, v in arms.items()), flush=True)
    return rows


def cmd_run(a) -> None:
    if ROW not in n12.tries_rows():
        sys.exit(f"T57 is refused: {ROW} has no row in experiments/harness/TRIES.md (a row comes before its run)")
    n22 = RESULTS / a.n22 / "n22.json"
    if not gc.n22_passed(n22):
        sys.exit(f"T57 is refused: {gv._rel(n22)} does not hold a passing N22 (T26's route needs N22 passed)")
    out = RESULTS / a.out
    name = TASK_FILE.format(world=a.world)
    gv.refuse_existing(out / "tasks", (name,))
    started = gv.tree_state()
    seeds = gv.parse_seeds(a.seeds)
    data = Path(a.data)
    t = {x.key: x for x in harness.load_tasks()["T1:test"]}[TASK]
    boards = load_boards()
    genes = boards[t.panel_board].genes
    ref_stage = Stage.from_anndata(read_stage(data / t.reference, genes), genes)
    inputs = {r: load_stage(data / f, lineage=False) for r, f in t.inputs.items()}
    role = harness.base_role(t)
    try:  # a thinned reference comes before anything else is made (try T54's order)
        inputs, ref_stage, truth, contamination, thinning, truth_name = world_route(
            a.world, t, data, genes, inputs, role, ref_stage, lambda: read_stage(data / t.reference, None))
    except ValueError as e:
        sys.exit(f"T57 is refused on {a.world}: {e}")
    board = harness.task_board(t, boards, truth.n, ref_stage.n)
    stamp_inputs = [data / t.target, data / t.reference, *[data / f for f in t.inputs.values()]]
    res = {"label": LABEL_T57, "task": t.key, "world": a.world, "truth": truth_name,
           "truth_sha256": ls.prediction_sha256(truth),
           "reference_of": REFERENCE_OF.get(a.world), "truth_rows": int(truth.n),
           "truth_contamination": contamination, "seeds": list(seeds), "thinning": thinning}
    floors, failures = ls.all_anchors(t, board, truth, ref_stage, seeds)
    if failures:
        res.update({"anchor_failures": failures, "rows": [], "stamp": gv.stamp(started, stamp_inputs)})
        zr.write_task_file(out, name, res)
        sys.exit(f"T57 is refused on {a.world}: the floor's anchors fail on seeds {[f['seed'] for f in failures]} "
                 "(recorded in the task file; T57 is void)")
    pool = inputs[role]
    preregistered = t1_scope.held_flat_genes()
    union = hygiene.load_union()
    held_m1 = gc.held_feng(union, preregistered)
    try:
        directions, forms, rec, dc_rec, ext = world_directions(list(pool.var_names), held_m1)
    except ValueError as e:
        sys.exit(str(e))
    groups = feng_dirs.label_groups(pool.obs["celltype"].astype(str).unique())
    fraction = feng_dirs.bracket_fraction(harness.STAGE_DAYS[t.reference], harness.STAGE_DAYS[t.target])
    view, cboard = harness.TaskView.of(t), harness._candidate_board(board)
    rows = seed_rows(view, cboard, board, truth, ref_stage, inputs, floors, seeds, directions, forms, groups, fraction,
                     held_m1, preregistered, label=a.world)
    res.update({"fraction": fraction, "fractions": arm_fractions(fraction), "multipliers": dict(ENTRY_MULTIPLIERS),
                "form_arms": {k: list(v) for k, v in FORM_ARMS.items()}, "external": ext, "directions": rec,
                "dc_directions": dc_rec, "rows": rows,
                "stamp": gv.stamp(started, [*stamp_inputs, t1_scope.HELD_FLAT_GENES, ts.R1_CONFIG,
                                            REPO_ROOT / union["file"], feng_dirs.clean_stats_path(), n22])})
    zr.write_task_file(out, name, res)


# ---------------------------------------------------------------------------------------------------------------------
# The reading


def _rows(r: Mapping) -> list:
    return ob._rows(r)


def committed_files(t50_tasks: Path = T50_TASKS, t54_tasks: Path = T54_TASKS, t55_tasks: Path = T55_TASKS,
                    t56_tasks: Path = T56_TASKS) -> dict[str, dict[str, Path]]:
    """The committed task files the read reproduces: 'entry' (b0 and the entry's arms: x1 to x8 in try T50's lib3 and
    lib2, x1 to x6 in try T54's ref67 and ref50), 'cn' (cn and cn_x2: try T56's lib3, ref67 and ref50, try T55's
    lib2) and 't55_lib3' (try T55's direction record)."""
    folder = {"t50": Path(t50_tasks), "t54": Path(t54_tasks), "t55": Path(t55_tasks), "t56": Path(t56_tasks)}
    return {"entry": {k: folder[s] / TASK_FILE.format(world=k) for k, s in ENTRY_SOURCE.items()},
            "cn": {k: folder[s] / TASK_FILE.format(world=k) for k, s in CN_SOURCE.items()},
            "t55_lib3": {LIB3: folder["t55"] / TASK_FILE.format(world=LIB3)}}


def dc_values(rec) -> dict | None:
    """dc's recorded sha256, p_g and a_n per group ((what, group) -> value); None without a full record."""
    try:
        if not rec:
            return None
        return {(x, g): v[x] for g, v in zr._groups(rec).items() for x in ("sha256", "p_g", "a_n")}
    except (KeyError, TypeError, AttributeError):
        return None


def dc_record_off(dc_rec, directions) -> list[str]:
    """dc's record against try T55's in the same file: per group dc's a_n equals the entry's recorded a_n, p_g is at
    most a_n (Cauchy-Schwarz), and dc's squared size over the entry's is at least 1 (it is 2 - (p_g/a_n)^2)."""
    try:
        entry = directions["entry"]
        out = []
        for g, v in zr._groups(dc_rec).items():
            a_n, p_g, size = v["a_n"], v["p_g"], v.get("size_ratio")
            if g not in entry or entry[g][pr.N] != a_n:
                out.append(f"dc's a_n in {g} ({a_n!r}) is not the entry's recorded a_n "
                           f"({(entry.get(g) or {}).get(pr.N)!r})")
            if not p_g <= a_n:
                out.append(f"dc's p_g in {g} ({p_g!r}) is above its a_n ({a_n!r})")
            if not isinstance(size, int | float) or not size >= 1.0 - 1e-12:
                out.append(f"dc's squared size over the entry's in {g} ({size!r}) is not at least 1")
        return out
    except (KeyError, TypeError, AttributeError) as e:
        return [f"dc's record cannot be checked against the entry's ({e!r})"]


def _key(v: Mapping) -> str:
    return json.dumps(sorted((f"{a}|{b}", x) for (a, b), x in v.items()))


def _thinning_reasons(k: str, th, results: Mapping[str, Mapping]) -> list[str]:
    """Try T54's checks on a thinned world's record, and for a library-2 world its twin's p and stream."""
    reasons = []
    spec, th = THINNING_OF[k], th or {}
    if th.get("p") != spec["p"] or list(th.get("stream") or []) != list(spec["stream"]):
        reasons.append(f"{k}: the file records p {th.get('p')!r} and stream {th.get('stream')!r}, not {spec['p']!r} "
                       f"and {spec['stream']}")
    ratio = th.get("depth_ratio")
    if not isinstance(ratio, int | float) or not abs(ratio - spec["p"]) <= DEPTH_TOL:
        reasons.append(f"{k}: realised depth ratio {ratio!r} is not within {DEPTH_TOL:g} of p {spec['p']:.4f}")
    if th.get("any_cell_emptied") is not False:
        reasons.append(f"{k}: the thinning left a cell with no count (or did not record whether it did): "
                       f"any_cell_emptied {th.get('any_cell_emptied')!r}")
    rebuild = th.get("rebuild_largest_difference")
    if not isinstance(rebuild, int | float) or not rebuild <= REBUILD_TOL:
        reasons.append(f"{k}: the rebuild check's largest difference {rebuild!r} is not at most {REBUILD_TOL:g}")
    if k in TWIN:
        twin = results.get(TWIN[k])
        tt = (twin or {}).get("thinning") or {}
        if twin is None or th.get("p") != tt.get("p") or list(th.get("stream") or []) != list(tt.get("stream") or []):
            reasons.append(f"{k}: the thinning record's p and stream ({th.get('p')!r}, {th.get('stream')!r}) are not "
                           f"{TWIN[k]}'s ({tt.get('p')!r}, {tt.get('stream')!r})")
    return reasons


def validity(results: Mapping[str, Mapping], committed: Mapping[str, Mapping[str, Mapping]]) -> dict:
    """The void conditions. results: world -> task file; committed: 'entry' and 'cn' (world -> committed file) and
    't55_lib3' (lib3 -> try T55's lib3 file), a missing file absent."""
    reasons = []
    missing = [k for k in WORLDS if k not in results]
    if missing:
        reasons.append(f"missing worlds {missing}")
    want_fr = arm_fractions(FRACTION)
    repro, values, dcv = {}, {}, {}
    for k, r in results.items():
        rows = _rows(r)
        if r.get("anchor_failures"):
            reasons.append(f"{k}: the floor's anchors invert or go missing on seeds "
                           f"{[f['seed'] for f in r['anchor_failures']]}")
        seeds = [row["seed"] for row in rows]
        if sorted(seeds) != list(SEEDS) or list(r.get("seeds", seeds)) != list(SEEDS):
            reasons.append(f"{k}: seeds {seeds}, not {list(SEEDS)}")
        absent = sorted({a for row in rows for a in ARMS if a not in row["arms"]})
        if absent:
            reasons.append(f"{k}: arms {absent} missing")
        if r.get("task") != TASK:
            reasons.append(f"{k}: the file holds task {r.get('task')!r}")
        if r.get("world") != k:
            reasons.append(f"{k}: the file holds world {r.get('world')!r}")
        if k in TRUTH_OF and r.get("truth") != TRUTH_OF[k]:
            reasons.append(f"{k}: the file holds truth {r.get('truth')!r}, not {TRUTH_OF[k]}")
        if not isinstance(r.get("fraction"), int | float) or abs(r["fraction"] - FRACTION) > FRACTION_TOL:
            reasons.append(f"{k}: fraction {r.get('fraction')!r}, not {FRACTION:g}")
        got_fr = r.get("fractions") or {}
        bad = [a for a, f in want_fr.items()
               if not isinstance(got_fr.get(a), int | float) or abs(got_fr[a] - f) > FRACTION_TOL]
        if bad:
            reasons.append(f"{k}: arm fractions {got_fr}, not {want_fr} (wrong fraction on {bad})")
        if any(not isinstance((row.get("sha256") or {}).get(a), str) for row in rows for a in ARMS):
            reasons.append(f"{k}: sha256 not recorded for every arm")
        # the committed scores
        repro[k] = {}
        if k in ENTRY_SOURCE:
            want = REPRO_ENTRY[ENTRY_SOURCE[k]]
            rep = pdr.reproduction(rows, committed.get("entry", {}).get(k), want)
            repro[k]["entry_arms"] = rep
            if not rep["holds"]:
                reasons.append(f"{k}: {', '.join(want)} do not reproduce try T{ENTRY_SOURCE[k][1:]}'s "
                               f"committed scores within {REPRO_TOL:g} per seed")
        if k in CN_SOURCE:
            rep = pdr.reproduction(rows, committed.get("cn", {}).get(k), REPRO_CN)
            repro[k]["cn_arms"] = rep
            if not rep["holds"]:
                reasons.append(f"{k}: {', '.join(REPRO_CN)} do not reproduce try T{CN_SOURCE[k][1:]}'s "
                               f"committed scores within {REPRO_TOL:g} per seed")
        # the thinning
        th = r.get("thinning")
        if k in REFERENCE_OF:
            reasons.extend(_thinning_reasons(k, th, results))
        elif th is not None:
            reasons.append(f"{k}: the file records a thinning")
        # the entry's directions as tries T50 and T54 load them and as try T55 loads them
        agree = (r.get("external") or {}).get("entry_agrees_with_t55s")
        if not isinstance(agree, Mapping) or not agree or not all(x is True for x in agree.values()):
            off = sorted(g for g, x in agree.items() if x is not True) if isinstance(agree, Mapping) else agree
            reasons.append(f"{k}: the entry's directions as tries T50 and T54 load them and as try T55 loads them do "
                           f"not have the same sha256 in every group (or it is not recorded): {off!r}")
        # a form's arm built at the entry's arm (idle form, wrong directions or wrong fraction)
        if rows and not absent:
            same = [a for a, (_, p) in FORM_ARMS.items()
                    if all(row["arms"][a]["score"] == row["arms"][p]["score"] for row in rows)]
            if same:
                reasons.append(f"{k}: form arms score exactly as the entry's arm of their fraction on every seed: "
                               f"{same}")
        values[k] = pr.direction_values(r.get("directions") or {})
        if values[k] is None:
            reasons.append(f"{k}: no full direction record")
        dcv[k] = dc_values(r.get("dc_directions"))
        if dcv[k] is None:
            reasons.append(f"{k}: no full record of the form dc")
        else:
            reasons.extend(f"{k}: {x}" for x in dc_record_off(r.get("dc_directions"), r.get("directions")))
        if not isinstance(r.get("truth_sha256"), str):
            reasons.append(f"{k}: the truth matrix's sha256 is not recorded")
    known = {k: _key(v) for k, v in values.items() if v is not None}
    if len(set(known.values())) > 1:
        reasons.append(f"the recorded directions' sha256 or multiples differ between the worlds' files {sorted(known)}")
    known_dc = {k: _key(v) for k, v in dcv.items() if v is not None}
    if len(set(known_dc.values())) > 1:
        reasons.append(f"dc's recorded sha256, p_g or a_n differ between the worlds' files {sorted(known_dc)}")
    t55 = pr.direction_values(((committed.get("t55_lib3") or {}).get(LIB3) or {}).get("directions") or {})
    if t55 is None:
        reasons.append("try T55's lib3 file holds no full direction record")
    else:
        off = sorted(k for k, v in known.items() if v != _key(t55))
        if off:
            reasons.append(f"the recorded directions' sha256 or multiples differ from try T55's lib3 file in {off}")
    for lib in (LIB3, LIB2):
        hashes = {k: r.get("truth_sha256") for k, r in results.items()
                  if TRUTH_OF.get(k) == lib and isinstance(r.get("truth_sha256"), str)}
        if len(set(hashes.values())) > 1:
            reasons.append(f"the truth matrix's sha256 differs between the worlds of {lib} {sorted(hashes)}")
    commits = {(r.get("stamp") or {}).get("commit") for r in results.values()}
    unclean = sorted(k for k, r in results.items() if (r.get("stamp") or {}).get("tree_clean") is not True)
    if unclean:
        reasons.append(f"tree not clean on {unclean}")
    if len(commits) > 1:
        reasons.append(f"the worlds ran at {len(commits)} commits")
    return {"valid": not reasons, "reasons": reasons, "reproduction": repro,
            "commits": sorted(c for c in commits if c)}


def cell_name(world: str, m: int) -> str:
    return f"{world} x{m}"


def comparison(r: Mapping, world: str, change: str, base: str) -> dict:
    """Try T50's comparison (try T46's single-task test on one world, the reversed test, the points by part and the raw
    metric differences) with the label ahead / behind / level."""
    return ob.comparison(r, world, change, base)


def holds(c: Mapping) -> bool:
    """A positive mean with at least HOLD_AHEAD of 10 seeds ahead."""
    return c["mean"] > 0 and c["seeds_ahead"] >= HOLD_AHEAD


def statement_size(cells: Mapping[str, Mapping]) -> str:
    """Statement 1's size: the cells ahead by try T46's test and the smallest and largest cell means."""
    ahead = sum(c["label"] == "ahead" for c in cells.values())
    means = [c["mean"] for c in cells.values()]
    return (f"ahead by try T46's test in {ahead} of {len(cells)} cells, cell means {min(means):+.2f} to "
            f"{max(means):+.2f}")


def statement_1(cells: Mapping[str, Mapping], form: str = CN) -> dict:
    """cells: cell name -> the form's difference from the entry. AHEAD_ALL_OF[form] when every cell holds;
    BEHIND_SOMEWHERE naming the cells behind by try T46's test; otherwise NOT_SEPARATED naming every cell that meets
    neither. The statement carries its size (statement_size) so that it is never quoted without it."""
    behind = [n for n, c in cells.items() if c["label"] == "behind"]
    neither = [n for n, c in cells.items() if not holds(c) and c["label"] != "behind"]
    if all(holds(c) for c in cells.values()):
        said = AHEAD_ALL_OF[form]
    elif behind:
        said = f"{BEHIND_SOMEWHERE} ({', '.join(behind)})"
    else:
        said = f"{NOT_SEPARATED} ({', '.join(neither)})"
    size = statement_size(cells)
    return {"statement": f"{said}; {size}", "said": said, "size": size,
            "ahead_by_test": [n for n, c in cells.items() if c["label"] == "ahead"], "behind": behind,
            "neither": neither, "quantity": QUANTITY[form]}


def statement_2(cells: Mapping[str, Mapping]) -> dict:
    """cells: cell name -> dc minus cn. The cells where dc is ahead, behind and level by try T46's test, counted and
    named."""
    by = {lab: [n for n, c in cells.items() if c["label"] == lab] for lab in ("ahead", "behind", "level")}
    summary = "; ".join(f"{lab} in {len(ns)} of {len(cells)} cells" + (f" ({', '.join(ns)})" if ns else "")
                        for lab, ns in by.items())
    return {"summary": summary, "counts": {lab: len(ns) for lab, ns in by.items()}, "cells": by}


def best_multiple(mean_scores: Mapping[str, float], arm_at: Callable[[int], str],
                  multiples: Sequence[int] = MULTIPLES) -> int:
    """The multiple with the highest mean score among `multiples` (the smallest on a tie)."""
    return max(multiples, key=lambda m: (mean_scores[arm_at(m)], -m))


def c_ratio(dc_rec: Mapping) -> float:
    """The largest a_n/p_g over the groups: how much further along c than the entry a form moves at one multiple
    (both forms put a_n on c, the entry p_g); infinite where a p_g is not positive."""
    return max((v["a_n"] / v["p_g"] if v["p_g"] > 0 else float("inf")) for v in zr._groups(dc_rec).values())


def statement_3(r: Mapping, world: str, mean_scores: Mapping[str, float]) -> dict:
    """Per form, the form's best multiple (1 to 4) against the entry's best over ENTRY_RANGE (1 to 6 and 8), paired by
    seed. NOT_READABLE where the entry's best is at the edge of its range and the c-matched multiple (the form's
    multiple x c_ratio) lies beyond it: there a larger plain step than any read could still give the form's lead."""
    e = best_multiple(mean_scores, ENTRY_RANGE.__getitem__, tuple(ENTRY_RANGE))
    ratio = c_ratio(r["dc_directions"])
    out = {}
    for f in FORMS:
        m = best_multiple(mean_scores, lambda x, f=f: form_arm(f, x))
        matched = m * ratio
        c = comparison(r, world, form_arm(f, m), ENTRY_RANGE[e])
        if e == ENTRY_EDGE and matched > ENTRY_EDGE:
            said = (f"{NOT_READABLE}: the entry's best is at the edge of the range read (x{ENTRY_EDGE}) and the form "
                    f"at x{m} moves along c as far as the entry at x{matched:.2f}")
        else:
            said = ABOVE if holds(c) else NOT_ABOVE
        out[f] = {"form_multiple": m, "entry_multiple": e, "c_matched_multiple": matched, "comparison": c,
                  "statement": said}
    return out


def t57_reading(results: Mapping[str, Mapping], committed: Mapping[str, Mapping[str, Mapping]]) -> dict:
    v = validity(results, committed)
    out = {"valid": v["valid"], "validity": v}
    if not v["valid"]:
        out["outcome"] = "VOID"
        return out
    diffs = {f: {cell_name(k, m): comparison(results[k], k, form_arm(f, m), ENTRY_AT[m]) for k, m in CELLS}
             for f in FORMS}
    s1 = {f: statement_1(diffs[f], f) for f in FORMS}
    dc_cn = {cell_name(k, m): comparison(results[k], k, form_arm(DC, m), form_arm(CN, m)) for k, m in CELLS}
    s2 = statement_2(dc_cn)
    means = {k: {a: float(np.mean([row["arms"][a]["score"] for row in _rows(results[k])])) for a in ARMS}
             for k in WORLDS}
    s3 = {k: statement_3(results[k], k, means[k]) for k in WORLDS}
    out.update({"outcome": {"statement_1": {f: s["statement"] for f, s in s1.items()}, "statement_2": s2["summary"],
                            "statement_3": {k: {f: x["statement"] for f, x in s3[k].items()} for k in WORLDS}},
                "statement_1": s1, "statement_2": s2, "statement_3": s3, "differences": diffs,
                "dc_minus_cn": dc_cn, "mean_scores": means,
                "thinning": {k: {x: (results[k].get("thinning") or {}).get(x) for x in
                                 ("p", "stream", "depth_ratio", "rebuild_largest_difference", "any_cell_emptied")}
                             for k in REFERENCE_OF},
                "directions": results[LIB3]["directions"], "dc_directions": results[LIB3]["dc_directions"]})
    return out


def _raw_line(c: Mapping) -> str:
    return ", ".join(f"{m} {c['metrics'].get(m, float('nan')):+.4g}" for m in ls.WEIGHTS)


def _difference_line(c: Mapping) -> str:
    return (f"{c['change']} minus {c['base']}: {c['mean']:+.2f} points, ahead {c['seeds_ahead']} of {c['of']}, per "
            f"seed {c['range'][0]:+.2f} to {c['range'][1]:+.2f}, test: {c['label']}; de terms "
            f"{c['points']['de']:+.2f}, mmd_u and variogram {c['points']['mmd_u_variogram']:+.2f} (read by nothing); "
            f"raw metrics {_raw_line(c)} (read by nothing)")


def t57_text(res: Mapping) -> str:
    """The reading as text; no cell count is printed."""
    r = res["reading"]
    lines = [res["label"], f"Outcome: {'VOID' if not r['valid'] else 'see the three statements'}"]
    v = r["validity"]
    if not r["valid"]:
        lines.append("Void: " + "; ".join(v["reasons"]))
    for k, rep in v["reproduction"].items():
        for what, x in rep.items():
            lines.append(f"{k}: {what.replace('_', ' ')}, largest per-seed difference from the committed scores: "
                         + ", ".join(f"{a} {d:.3g}" for a, d in x["largest_difference"].items())
                         + f" (gate {REPRO_TOL:g})")
    if r["valid"]:
        for f in FORMS:
            s = r["statement_1"][f]
            lines.append(f"Statement 1, {f} ({s['quantity']}): {s['statement']}")
            if s["behind"] and s["neither"]:
                lines.append("    also meeting neither condition: " + ", ".join(s["neither"]))
        lines.append(f"Statement 2, dc minus cn per cell: {r['statement_2']['summary']}")
        lines.append("Statement 3, beyond a larger step (the form's best of x1 to x4 against the entry's best of x1 to "
                     "x6 and x8, paired by seed; the entry's multiple that moves along c as far as the form):")
        for k in WORLDS:
            for f, x in r["statement_3"][k].items():
                c = x["comparison"]
                lines.append(f"    {k}, {f} at x{x['form_multiple']} against the entry at x{x['entry_multiple']}: "
                             f"{c['mean']:+.2f} points, ahead {c['seeds_ahead']} of {c['of']}, c-matched multiple "
                             f"x{x['c_matched_multiple']:.2f}: {x['statement']}")
        lines.append("Per form and cell, the form minus the entry's arm of the same fraction (try T46's single-task "
                     "test inside the world):")
        for f in FORMS:
            for n, c in r["differences"][f].items():
                lines.append(f"    {f}, {n}: " + _difference_line(c))
        lines.append("Per cell, dc minus cn:")
        for n, c in r["dc_minus_cn"].items():
            lines.append(f"    {n}: " + _difference_line(c))
        lines.append("Reported and read by nothing:")
        for k in WORLDS:
            lines.append(f"    world {k}: mean scores "
                         + ", ".join(f"{a} {x:.2f}" for a, x in r["mean_scores"][k].items()))
        for k, th in r["thinning"].items():
            lines.append(f"    world {k}: the reference of try T54's {REFERENCE_OF[k]}, realised depth ratio "
                         f"{th['depth_ratio']:.4f} (p {THINNING_OF[k]['p']:.4f}, stream {th['stream']}), rebuild "
                         f"check's largest difference {th['rebuild_largest_difference']:.3g}, any cell left with no "
                         f"count: {th['any_cell_emptied']}")
        d, dc = r["directions"], r["dc_directions"]
        lines.append("    the forms per group (the same in the six files; try T55's equal to its lib3 file): "
                     + "; ".join(
            f"{g}: a_n {e[pr.N]:.4f}, p_g {dc[g]['p_g']:.4f}, dc's squared size over the entry's "
            f"{zr._f(dc[g].get('size_ratio'), '.4f')}" for g, e in d["entry"].items() if g in dc))
    lines.append("Choices the brief left open:")
    lines.extend(f"  - {x}" for x in CHOICES)
    lines.extend(f"Limit: {x}" for x in LIMITS)
    lines.append("Limit: in-sample labelled diagnostic; never adoption evidence; puts nothing forward by itself; it "
                 "can only go to the team captain; no outcome is evidence on Task 1's read depth, which stays closed")
    lines.append(f"Commit {res['stamp']['commit']}, tree clean {res['stamp']['tree_clean']}, UTC {res['stamp']['utc']}")
    return "\n".join(lines) + "\n"


def read(out: Path, t50_tasks: Path = T50_TASKS, t54_tasks: Path = T54_TASKS, t55_tasks: Path = T55_TASKS,
         t56_tasks: Path = T56_TASKS) -> dict:
    """Read the worlds present in out/tasks against the committed files of tries T50, T54, T55 and T56; write
    t57.{json,txt}."""
    gv.refuse_existing(out, (f"{STEM}.json", f"{STEM}.txt"))
    started = gv.tree_state()
    present = {k: out / "tasks" / TASK_FILE.format(world=k) for k in WORLDS
               if (out / "tasks" / TASK_FILE.format(world=k)).exists()}
    results = {k: json.loads(p.read_text()) for k, p in present.items()}
    files = committed_files(t50_tasks, t54_tasks, t55_tasks, t56_tasks)
    found = {what: {k: p for k, p in by.items() if p.exists()} for what, by in files.items()}
    committed = {what: {k: json.loads(p.read_text()) for k, p in by.items()} for what, by in found.items()}
    res = {"label": LABEL_T57, "choices": list(CHOICES), "limits": list(LIMITS),
           "reading": t57_reading(results, committed), "stamps": {k: v.get("stamp") for k, v in results.items()}}
    res["stamp"] = gv.stamp(started, [*present.values(), *(p for by in found.values() for p in by.values())])
    gv.write_outputs(out, STEM, res, t57_text(res))
    return res


def cmd_read(a) -> None:
    res = read(RESULTS / a.out)
    if res["reading"]["outcome"] == "VOID":
        sys.exit("T57's reading is void and goes to the captain")


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("run", help="one world of try T57")
    p.add_argument("--world", required=True, choices=list(WORLDS))
    p.add_argument("--out", required=True)
    p.add_argument("--n22", default="n22_additivity")
    p.add_argument("--seeds", default=",".join(map(str, SEEDS)))
    p.add_argument("--data", default=str(harness.DATA_DIR))
    p = sub.add_parser("read", help="try T57's reading")
    p.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    {"run": cmd_run, "read": cmd_read}[a.cmd](a)


if __name__ == "__main__":
    main()
