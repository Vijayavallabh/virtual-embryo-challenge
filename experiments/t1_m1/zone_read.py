#!/usr/bin/env python
"""Try T53 (experiments/harness/TRIES.md): M1's Feng direction with the heart zones held fixed between the bracket's two
ends. A labelled in-sample read on the routes of try T44 (experiments/t1_m1/step_read.py, the whole cut) and try T46
(experiments/t1_m1/lib_step_read.py, the chamber cut pooled and library by library). Built to brief AG as amended after
its first review (2026-11-17): zbal dropped, zw added, slope-matched forms and their permuted controls added.

M1 moves each covered cell along d, per cell group the log ratio of mean CP10k between Feng's E9.5 and E14.5 hearts
(feng_dirs.directions_from_stats). dirlib.feng_sel pools a stage's cells of the group's Feng labels over every heart
zone, and Feng's capture is not the same at the two ends (the ledger: mostly LV-zone cells at E9.5, zone-balanced at
E14.5), so a pooled direction holds, beside the change over five days, a change in which zones the cells come from.

Directions (per group on the entry's gene map, the config's pseudocount and min_cells_ext; the statistics' keys are
stage|zone|cell type|half and both halves are pooled, as feng_sel pools them; a key without exactly four parts is
refused):
  within-zone  the log ratio (dirlib.lfc_ratio) of the group's mean CP10k over its Feng labels' cells of zone z at
               E14.5 and at E9.5; defined where the group has at least min_cells_ext cells in z at both ends.
  zlv          the cardiomyocyte group (CM) takes its within-zone direction in LV; every other group keeps the entry's
               pooled direction. Void (the run refuses) if CM has no LV direction.
  zw           the form that differs from the entry by the zone mix only. Per group: w_z, the share of the group's E9.5
               cells (its Feng labels, both halves) in zone z; the E14.5 mean CP10k is sum_z w_z x the group's mean
               CP10k in zone z at E14.5, over the zones where the group has at least min_cells_ext cells at E14.5, the
               weights renormalised over those zones; the E9.5 mean is the entry's pooled one; the direction is
               dirlib.lfc_ratio of the two. The zones dropped and the share they held are recorded; a group that loses
               more than half (strictly) of its E9.5 weight this way keeps the entry's direction, and the file says so.
  F_sm         slope-matched: for a form F and each group it changes (direction sha256 not the entry's), F_g / s_g, s_g
               the slope through the origin of F_g on the entry's e_g over the moved genes. The run refuses if a changed
               group's s_g is not above 0.2. Unchanged groups keep the entry's own array.
  F_sm_perm    the control: per changed group, e_g + P(F_g^sm - e_g) on the moved genes (the entry's e_g elsewhere), P
               one permutation of the residual's values within T46's two strata (the board's scored genes, the other
               moved genes), groups in sorted order, the scored stratum first; a fresh generator per control and build
               seed on the stream [build seed, 53], which no other step uses.
The moved genes are those finite in both directions and not in M1's held list. The task file records, per form and
group: the zones used; the slope through the origin and the Pearson correlation of the form's direction on the entry's
over the moved genes; the form's squared size over the entry's on those genes; the sha256 of each direction; for zw the
zones dropped and their share of the group's E9.5 cells; for a slope-matched form the divisor s_g. Ratios, shares and
names only: no cell count is written or printed. No Feng stage other than the bracket's two ends is read
(share_read.feng_inputs, feng_dirs.bracket_stats).

  run    One truth, 10 seeds. pooled, lib1, lib2, lib3: the chamber cut T1/E9.5->E10.5:chamber exactly as T46's run
         makes them (lib_step_read.load_truth); whole: the cut T1/E9.5->E10.5 with its pooled target exactly as T44's
         run loads it. Arms per seed, all from one cell draw: b0, m1_r_gmom (fraction 0.2), m1_r_gmom_x15 and
         m1_r_gmom_x2 (0.3 and 0.4) as T44 and T46 build them; and ten new arms, each the build of m1_r_gmom or
         m1_r_gmom_x2 (gmom_combo.m1_r_arms) with other directions in place of the entry's and nothing else changed:
         m1_r_gmom_zlv, _zlv_x2, _zw, _zw_x2 (raw forms, 0.2 and 0.4); _zlv_sm, _zlv_sm_x2, _zw_sm, _zw_sm_x2
         (slope-matched, 0.2 and 0.4); _zlv_sm_perm, _zw_sm_perm (the controls, 0.2). Refused on an existing task file,
         unless the committed N22 passed, unless experiments/harness/TRIES.md holds row T53, and when a changed group's
         slope is not above 0.2.
  read   T53's reading of the five truths' files (below, fixed in the row before the run). Writes t53.json and t53.txt
         and refuses to overwrite them.

The reading. Void unless the five truths hold the ten seeds at one commit on a clean tree (tree_clean exactly true); no
file records an anchor failure; every arm is on every row with its fraction; on the four chamber truths b0, m1_r_gmom,
m1_r_gmom_x15 and m1_r_gmom_x2 equal, seed for seed, their committed scores in
experiments/t1_m1/results/t46_m1_lib_step/tasks/, and on whole theirs in experiments/t1_m1/results/t44_m1_step/tasks/
(the whole cut's file), each to 1e-9; no new arm scores exactly as the entry's arm of its fraction on every seed of a
truth; zlv and zlv_sm differ from the entry in CM and in no other group, zw and zw_sm in at least one group; every
recorded direction's sha256 (entry, forms, slope-matched forms; per group) is the same in the five files; each
slope-matched direction's recorded slope on the entry's is within 1e-9 of 1 in every changed group; no control's
prediction sha256 equals its real arm's on every seed of a truth.
Tests. On the two cuts (pooled and whole): the adoption rule's test (vec.adoption.keep_change, both cuts) of an arm
against its comparison arm: "passes"; with the sign reversed "loses"; else "neither". On each library truth: T46's
single-task test (ahead, behind, level).
The statement, per form (zlv, zw), on the slope-matched arms: "the zone-fixed content scores better at equal size" if
_F_sm passes against m1_r_gmom and _F_sm_x2 passes against m1_r_gmom_x2; "the zone-fixed content scores worse at equal
size" if both lose; otherwise "no difference at the bar"; with the four comparisons' numbers and each library truth's
labels. Where a form reads better, its gene statement: _F_sm against _F_sm_perm on the two cuts: "the gain needs the
form's genes" if it passes, otherwise "a scrambled departure gains as well or is not separated"; where it does not, the
comparison is reported and carries no statement. The raw arms are reported against the pooled arm of their fraction
with the recorded slope and the size-alone figure (the gain a uniform multiple of the entry's step equal to the slope,
twice the slope for an _x2 arm, would have by linear interpolation between the uniform arms' multiples 1, 1.5 and 2;
not computed outside 1 to 2, the text says why); they carry no statement. Reported: per truth every arm's mean board
score; each comparison's paired difference with seeds ahead, per-seed range and sd, every metric's mean skill difference
and the parts on the two de terms and on mmd_u and the variogram. The read puts nothing forward by itself.

NOT an adoption result: a labelled in-sample diagnostic, as T44 to T46 are; Feng's E9.5 cells lie inside these
pseudo-tasks' external windows and the cuts are ones earlier T1 choices were made on. The result can only go to the team
captain. The limits are LIMITS below, repeated in the text.

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
FORMS = ("zlv", "zw")  # the raw forms
SM = {f: f"{f}_sm" for f in FORMS}  # their slope-matched versions
PERM = {f: f"{f}_sm_perm" for f in FORMS}  # the controls of the slope-matched versions
RECORDED = (*FORMS, *SM.values())  # the directions the task file records (the controls depend on the seed)
ONLY_CM = ("zlv", SM["zlv"])  # these must differ from the entry in CM and in no other group
MIN_SLOPE = 0.2  # the run refuses a changed group whose slope on the entry's is not above this
SLOPE_TOL = 1e-9  # a slope-matched direction's recorded slope must be this close to 1
MAX_DROPPED = 0.5  # zw: a group that loses more than this share of its E9.5 weight keeps the entry's direction
CONTROL_STREAM = 53  # [build seed, 53]: the controls' permutations; a stream no other step uses
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
# new arm -> (its directions' key, the pooled arm of its fraction)
NEW_ARMS = {
    **{f"{BASE}_{f}{s}": (f, p) for f in FORMS for s, p in (("", BASE), ("_x2", X2))},
    **{f"{BASE}_{SM[f]}{s}": (SM[f], p) for f in FORMS for s, p in (("", BASE), ("_x2", X2))},
    **{f"{BASE}_{PERM[f]}": (PERM[f], BASE) for f in FORMS},
}
RAW_ARMS = tuple(a for a, (d, _) in NEW_ARMS.items() if d in FORMS)
SM_ARM = {f: (f"{BASE}_{SM[f]}", f"{BASE}_{SM[f]}_x2") for f in FORMS}
PERM_ARM = {f: f"{BASE}_{PERM[f]}" for f in FORMS}
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
BETTER = "the zone-fixed content scores better at equal size"
WORSE = "the zone-fixed content scores worse at equal size"
NO_DIFFERENCE = "no difference at the bar"
NEEDS_GENES = "the gain needs the form's genes"
NOT_SEPARATED = "a scrambled departure gains as well or is not separated"
LABEL_T53 = ("try T53, a labelled in-sample diagnostic, never an adoption result: M1 (with R1, R2 and g_mom's factors, "
             "the entry's form) with its Feng direction zone-fixed (zlv: the cardiomyocyte group's LV-zone direction; "
             "zw: every group's E14.5 mean taken at its E9.5 zone shares), raw and slope-matched to the entry's, at "
             "the step and twice the step, with permuted controls of the slope-matched forms, on T44's and T46's "
             "routes (both E10.5 cuts, the chamber cut also library by library); Feng's E9.5 cells lie inside these "
             "pseudo-tasks' external windows and the cuts are ones earlier T1 choices were made on; the result can "
             "only go to the team captain")
CHOICES = (
    "The size-alone figure's multiple is the recorded slope for an arm at the entry's fraction and twice the slope for "
    "an _x2 arm; it is not computed outside the uniform arms' multiples 1 to 2.",
    "Slopes and the size-alone figure are per group: the route gives no covered-cell weights.",
    "Labels on each library truth: T46's single-task test of each arm against its comparison arm (ahead, behind, "
    "level).",
    "A statistics key without exactly the four parts stage|zone|cell type|half is refused (the run stops).",
    "A form changes a group when its direction's sha256 differs from the entry's; slope-matching and the controls "
    "touch changed groups only, and the slope-matched direction is the whole vector divided by s_g.",
    "zw weighs the zones holding any of the group's E9.5 cells; a group keeps the entry's direction when the dropped "
    "share is strictly above 0.5.",
    "A control equals the entry's direction outside the moved genes; its permutation runs on the reference role's gene "
    "order with the board's panel as the scored stratum (as T46's placebo), a fresh generator per control on "
    "[build seed, 53].",
    "The gene statement's test is the adoption rule on both cuts with the control as the base arm.",
)
LIMITS = (
    "One released day, in-sample, on cuts earlier choices were made on.",
    "A zone-fixed direction differs from the pooled one in size as well as in content, and larger steps score better "
    "on these cuts (tries T44 to T46), so a raw arm's gain or loss can be size: the slope and the interpolated figure "
    "are printed for that and decide nothing.",
    "Within-zone directions rest on fewer cells than the pooled one and are noisier.",
    "Which zone mix the hidden test's cells have is not known from released data, and neither form claims to match it.",
    "The two fractions share the seeds, the draw and the direction, so agreement at both is one result read twice.",
    "'No difference at the bar' does not show that the zone mix is irrelevant, and the bar's power at the noise of "
    "within-zone means is not stated.",
    "zlv is the ventricular cardiomyocyte direction of one zone and drops atrial cells as well as fixing the mix.",
    "zw holds each group's zone mix at its E9.5 shares and rests at E14.5 mostly on that group's LV cells.",
    "A slope-matched direction is not the form as estimated but that form rescaled, and a request could only name what "
    "was read.",
    "The cuts' own truth may hold a real change of atrial share, so neither outcome says which direction suits the "
    "hidden test.",
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


def zw_directions(st: Mapping, entry: Mapping[str, np.ndarray], genes, cfg: Mapping,
                  bracket: tuple[str, str] = BRACKET) -> tuple[dict, dict, dict]:
    """Per entry group: the E14.5 mean CP10k at the group's E9.5 zone shares (over the zones with at least
    min_cells_ext of its cells at E14.5, the shares renormalised over them) against the entry's pooled E9.5 mean, as
    dirlib.lfc_ratio on `genes`. A group whose dropped zones held more than MAX_DROPPED of its E9.5 cells keeps the
    entry's own array. Returns the directions, the zones used ([] = the entry's) and per group the zones dropped, their
    share and whether the group fell back (shares only, no count)."""
    zones = zones_of(st)
    idx = dl.sym_index(st["genes"], list(genes))
    eps, mc = cfg["metrics"]["pseudocount_cp10k"], cfg["metrics"]["min_cells_ext"]
    out, used, info = {}, {}, {}
    for g in sorted(entry):
        fl = cfg["external_labels"][g][0]
        a = dl.feng_sel(st, bracket[0], fl)  # the entry's pooled E9.5 mean
        kept, dropped = [], []
        for z in zones:
            az = zone_sel(st, bracket[0], z, fl)
            if az is None:
                continue
            w = az["n"] / a["n"]
            bz = zone_sel(st, bracket[1], z, fl)
            (kept if bz is not None and bz["n"] >= mc else dropped).append((z, w, bz))
        lost = float(sum(w for _, w, _ in dropped))
        fallback = lost > MAX_DROPPED or not kept
        if fallback:
            out[g], used[g] = entry[g], []
        else:
            total = sum(w for _, w, _ in kept)
            cp = sum(w * bz["cp"] for _, w, bz in kept) / total
            out[g] = dl.on_t1(dl.lfc_ratio(a, {"cp": cp}, eps), idx)
            used[g] = [z for z, _, _ in kept]
        info[g] = {"dropped": [z for z, _, _ in dropped], "dropped_share": lost, "fallback": fallback}
    return out, used, info


def direction_sha256(v) -> str:
    x = np.ascontiguousarray(np.asarray(v, dtype=np.float64))
    h = hashlib.sha256()
    h.update(f"{x.dtype.str}{x.shape}".encode())
    h.update(x.tobytes())
    return h.hexdigest()


def changed_groups(entry: Mapping, dirs: Mapping) -> list[str]:
    """The groups whose direction's sha256 is not the entry's, sorted."""
    return sorted(g for g in dirs if direction_sha256(dirs[g]) != direction_sha256(entry[g]))


def _num(x) -> float | None:
    return float(x) if x is not None and np.isfinite(x) else None


def moved_mask(form: np.ndarray, entry: np.ndarray, genes: Sequence[str], held: Sequence[str]) -> np.ndarray:
    """The moved genes: finite in both directions and not in M1's held list."""
    held_set = set(map(str, held))
    return (np.isfinite(form) & np.isfinite(entry) & np.array([str(g) not in held_set for g in genes], dtype=bool))


def compare_directions(form: np.ndarray, entry: np.ndarray, genes: Sequence[str], held: Sequence[str]) -> dict:
    """Over the moved genes: the slope through the origin and the Pearson correlation of the form's direction on the
    entry's (share_read.origin_slope), and the form's squared size over the entry's."""
    f, e = np.asarray(form, dtype=np.float64), np.asarray(entry, dtype=np.float64)
    moved = moved_mask(f, e, genes, held)
    slope, r = share_read.origin_slope(f[moved], e[moved])
    ee = float(e[moved] @ e[moved])
    return {"slope": _num(slope), "pearson": _num(r),
            "size_ratio": _num(float(f[moved] @ f[moved]) / ee) if ee > 0 else None}


def slope_matched(entry: Mapping[str, np.ndarray], dirs: Mapping[str, np.ndarray], genes: Sequence[str],
                  held: Sequence[str]) -> tuple[dict, dict]:
    """Each changed group's direction divided by its slope through the origin on the entry's over the moved genes; the
    unchanged groups keep the entry's own array. Refused (ValueError) where a changed group's slope is not above
    MIN_SLOPE. Returns the directions and the divisors."""
    out, divisors = dict(entry), {}
    for g in changed_groups(entry, dirs):
        s = compare_directions(dirs[g], entry[g], genes, held)["slope"]
        if s is None or not s > MIN_SLOPE:
            raise ValueError(f"group {g}'s slope on the entry's direction is {s}, not above {MIN_SLOPE:g}")
        out[g] = np.asarray(dirs[g], dtype=np.float64) / s
        divisors[g] = s
    return out, divisors


def control_directions(entry: Mapping[str, np.ndarray], sm: Mapping[str, np.ndarray], genes: Sequence[str],
                       held: Sequence[str], scored: Sequence[str], seed) -> dict[str, np.ndarray]:
    """Per changed group: the entry's direction plus the residual F_sm - e permuted within the moved genes in two strata
    (the board's scored genes `scored`, then the other moved genes); the entry's direction elsewhere. Groups in sorted
    order, one generator default_rng([*seed, CONTROL_STREAM]); unchanged groups keep the entry's own array."""
    scored_set = set(map(str, scored))
    on_board = np.array([str(g) in scored_set for g in genes], dtype=bool)
    rng = np.random.default_rng([*np.atleast_1d(seed).astype(int).tolist(), CONTROL_STREAM])
    out = dict(entry)
    for g in changed_groups(entry, sm):
        e = np.asarray(entry[g], dtype=np.float64)
        f = np.asarray(sm[g], dtype=np.float64)
        moved = moved_mask(f, e, genes, held)
        r = np.where(moved, f - e, 0.0)
        v = e.copy()
        for stratum in (moved & on_board, moved & ~on_board):
            ix = np.flatnonzero(stratum)
            v[ix] = e[ix] + r[ix][rng.permutation(len(ix))]
        out[g] = v
    return out


def seed_controls(ctx: harness.BuildContext, entry: Mapping, forms: Mapping, held_m1: Sequence[str]) -> dict:
    """The controls of one build seed, on the reference role's gene order, M1's held genes and the board's scored
    genes (the order and strata T46's placebo uses)."""
    genes = list(ctx.inputs[ctx.task.reference_role].var_names)
    return {PERM[f]: control_directions(entry, forms[SM[f]], genes, held_m1, list(ctx.board.genes), ctx.seed)
            for f in FORMS}


def direction_record(entry: Mapping, forms: Mapping[str, tuple[dict, dict]], genes, held, zones: Sequence[str],
                     zw_info: Mapping | None = None, divisors: Mapping | None = None) -> dict:
    """The task file's direction record: names, ratios, shares and hashes only. forms: key -> (directions, zones
    used)."""
    rec = {"zones": list(zones), "entry": {g: {"sha256": direction_sha256(entry[g])} for g in sorted(entry)}}
    for form, (dirs, used) in forms.items():
        rec[form] = {g: {"zones": list(used[g]), "pooled": not used[g],
                         **compare_directions(dirs[g], entry[g], genes, held), "sha256": direction_sha256(dirs[g])}
                     for g in sorted(dirs)}
        rec.setdefault("changed", {})[form] = changed_groups(entry, dirs)
    for g, v in (zw_info or {}).items():
        rec["zw"][g].update(v)
    for form, d in (divisors or {}).items():
        for g, s in d.items():
            rec[form][g]["divisor"] = float(s)
    if "zw" in rec:
        rec["zw_fallback"] = sorted(g for g, v in rec["zw"].items() if v.get("fallback"))
    return rec


def all_directions(st: Mapping, entry: Mapping, genes, cfg: Mapping, held: Sequence[str],
                   bracket: tuple[str, str] = BRACKET) -> tuple[dict, dict]:
    """The raw and slope-matched forms (key -> directions) and their record. Raises ValueError where the run refuses."""
    within = within_zone_directions(st, genes, cfg, entry, bracket)
    zlv, zlv_used = zlv_directions(entry, within)
    zw, zw_used, zw_info = zw_directions(st, entry, genes, cfg, bracket)
    raw = {"zlv": (zlv, zlv_used), "zw": (zw, zw_used)}
    sm, divisors = {}, {}
    for f, (dirs, used) in raw.items():
        d, s = slope_matched(entry, dirs, genes, held)
        sm[SM[f]] = (d, used)
        divisors[SM[f]] = s
    forms = {**raw, **sm}
    rec = direction_record(entry, forms, genes, held, zones_of(st), zw_info, divisors)
    return {k: d for k, (d, _) in forms.items()}, rec


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
    pooled arm of its fraction with its directions (a form, a slope-matched form or this seed's control) in place of
    the entry's and every other argument the same."""
    preds = {"b0": b0(ctx), **sr.step_arms(ctx, entry, groups, fraction, held_m1, held_factors,
                                           multipliers=UNIFORM_MULTIPLIERS)}
    dirs = {**forms, **seed_controls(ctx, entry, forms, held_m1)}
    fr = arm_fractions(fraction)
    for arm, (key, _) in NEW_ARMS.items():
        preds[arm] = gc.m1_r_arms(ctx, dirs[key], groups, fr[arm], held_m1, held_factors)[BASE]
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
    groups = feng_dirs.label_groups(pool.obs["celltype"].astype(str).unique())
    fraction = feng_dirs.bracket_fraction(harness.STAGE_DAYS[t.reference], harness.STAGE_DAYS[t.target])
    preregistered = t1_scope.held_flat_genes()
    union = hygiene.load_union()
    held_m1 = gc.held_feng(union, preregistered)
    try:
        form_dirs, directions = all_directions(st, entry, pool_genes, cfg, held_m1, BRACKET)
    except ValueError as e:
        sys.exit(f"T53 is refused: {e}")
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
                "external": external, "directions": directions, "control_stream": CONTROL_STREAM, "rows": rows,
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


def _groups(rec: Mapping) -> dict:
    """A form's per-group records (a Mapping check, so that a malformed record reads as absent)."""
    if not isinstance(rec, Mapping):
        raise TypeError("not a direction record")
    return dict(rec)


def differing_groups(directions: Mapping, form: str) -> list[str] | None:
    """The groups whose recorded direction for `form` differs from the entry's (by sha256); None without a record."""
    try:
        entry, rec = directions["entry"], _groups(directions[form])
        return sorted(g for g in set(entry) | set(rec)
                      if (entry.get(g) or {}).get("sha256") != (rec.get(g) or {}).get("sha256"))
    except (KeyError, TypeError, AttributeError):
        return None


def direction_hashes(directions: Mapping) -> dict | None:
    """(direction, group) -> sha256 over the entry and every recorded form; None without a full record."""
    try:
        out = {("entry", g): v["sha256"] for g, v in directions["entry"].items()}
        for form in RECORDED:
            out.update({(form, g): v["sha256"] for g, v in _groups(directions[form]).items()})
        return out
    except (KeyError, TypeError, AttributeError):
        return None


def validity(results: Mapping[str, Mapping], committed: Mapping[str, Sequence[Mapping]]) -> dict:
    """The void conditions. results: truth -> task file; committed: truth -> the committed file's rows (absent where the
    file is missing)."""
    reasons = []
    missing = [k for k in TRUTHS if k not in results]
    if missing:
        reasons.append(f"missing truths {missing}")
    want_fr = arm_fractions(FRACTION)
    hashes = {}
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
            same_pred = [PERM_ARM[f] for f in FORMS
                         if all((row.get("sha256") or {}).get(PERM_ARM[f]) is not None
                                and row["sha256"][PERM_ARM[f]] == row["sha256"].get(SM_ARM[f][0]) for row in rows)]
            if same_pred:
                reasons.append(f"{k}: controls with their real arm's prediction sha256 on every seed: {same_pred}")
        directions = r.get("directions") or {}
        for form in RECORDED:
            diff = differing_groups(directions, form)
            if diff is None:
                reasons.append(f"{k}: no direction record for {form}")
            elif form in ONLY_CM and diff != [CM]:
                reasons.append(f"{k}: {form}'s directions differ from the entry's in {diff}, not in {CM} alone")
            elif not diff:
                reasons.append(f"{k}: {form}'s directions are the entry's in every group")
            elif form in SM.values():
                off = [g for g in diff if not isinstance(_groups(directions[form])[g].get("slope"), int | float)
                       or abs(_groups(directions[form])[g]["slope"] - 1.0) > SLOPE_TOL]
                if off:
                    reasons.append(f"{k}: {form}'s recorded slope on the entry's is not within {SLOPE_TOL:g} of 1 "
                                   f"in {off}")
        hashes[k] = direction_hashes(directions)
    known = {k: h for k, h in hashes.items() if h is not None}
    if len({json.dumps(sorted((f"{a}|{b}", s) for (a, b), s in h.items())) for h in known.values()}) > 1:
        reasons.append(f"the recorded directions' sha256 differ between the truths' files {sorted(known)}")
    commits = {(r.get("stamp") or {}).get("commit") for r in results.values()}
    unclean = sorted(k for k, r in results.items() if (r.get("stamp") or {}).get("tree_clean") is not True)
    if unclean:
        reasons.append(f"tree not clean on {unclean}")
    if len(commits) > 1:
        reasons.append(f"the truths ran at {len(commits)} commits")
    return {"valid": not reasons, "reasons": reasons, "commits": sorted(c for c in commits if c)}


def cut_test(results: Mapping[str, Mapping], arm: str, base: str) -> dict:
    """The adoption rule on both cuts (step_read.compare): 'passes'; with the sign reversed 'loses'; else 'neither'."""
    cuts = {k: results[k] for k in CUTS}
    fwd, rev = sr.compare(cuts, base, arm), sr.compare(cuts, arm, base)
    test = PASSES if fwd["passes"] else LOSES if rev["passes"] else NEITHER
    return {"base": base, "change": arm, "test": test, "reason": fwd["rule"]["reason"],
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


def size_alone_reason(slope: float | None, pooled: str) -> str | None:
    """Why the size-alone figure is not defined (None where it is)."""
    what = "slope" if pooled == BASE else "twice the slope"
    if slope is None:
        return "slope not defined"
    m = slope * UNIFORM[pooled]
    if m < 1.0:
        return f"{what} below 1: outside the uniform arms 1 to 2"
    if m > 2.0:
        return f"{what} above 2: outside the uniform arms 1 to 2"
    return None


def comparison_pairs() -> dict[str, tuple[str, str]]:
    """Every comparison the read reports: name -> (arm, base arm)."""
    pairs = {a: (a, p) for a, (_, p) in NEW_ARMS.items() if a not in PERM_ARM.values()}
    pairs.update({f"{SM_ARM[f][0]}_vs_perm": (SM_ARM[f][0], PERM_ARM[f]) for f in FORMS})
    pairs.update({f"{PERM_ARM[f]}_vs_pooled": (PERM_ARM[f], BASE) for f in FORMS})
    return pairs


def t53_reading(results: Mapping[str, Mapping], committed: Mapping[str, Sequence[Mapping]]) -> dict:
    v = validity(results, committed)
    out = {"valid": v["valid"], "validity": v}
    if not v["valid"]:
        out["outcome"] = "VOID"
        return out
    means = {k: {a: float(np.mean([row["arms"][a]["score"] for row in results[k]["rows"]])) for a in ARMS}
             for k in TRUTHS}
    pairs = comparison_pairs()
    cut_tests = {n: cut_test(results, a, b) for n, (a, b) in pairs.items()}
    comps = {k: {n: lsr.truth_comparison(results[k], k, a, base=b) for n, (a, b) in pairs.items()} for k in TRUTHS}
    labels = {k: {n: library_label(comps[k][n]) for n in pairs} for k in LIBS}
    statements = {}
    for form in FORMS:
        arms = SM_ARM[form]
        st = form_statement([cut_tests[a]["test"] for a in arms])
        gene_name = f"{arms[0]}_vs_perm"
        gene = cut_tests[gene_name]
        statements[form] = {
            "statement": st, "comparisons": {a: cut_tests[a] for a in arms},
            "library_labels": {k: {a: labels[k][a] for a in arms} for k in LIBS},
            "gene_statement": (NEEDS_GENES if gene["test"] == PASSES else NOT_SEPARATED) if st == BETTER else None,
            "gene_comparison": gene, "gene_library_labels": {k: labels[k][gene_name] for k in LIBS}}
    rec = results[lsr.POOLED]["directions"]
    raw = {}
    for a in RAW_ARMS:
        form, p = NEW_ARMS[a]
        groups = sorted(_groups(rec[form]))
        raw[a] = {"comparison": cut_tests[a], "library_labels": {k: labels[k][a] for k in LIBS},
                  "size_alone": {k: {g: {"slope": rec[form][g]["slope"],
                                         "what_size_alone_would_give": size_alone(
                                             means[k], None if rec[form][g]["slope"] is None
                                             else rec[form][g]["slope"] * UNIFORM[p], p),
                                         "why_not": size_alone_reason(rec[form][g]["slope"], p)}
                                     for g in groups} for k in TRUTHS}}
    out.update({"outcome": {f: s["statement"] for f, s in statements.items()}, "statements": statements,
                "raw_arms": raw, "mean_scores": means, "comparisons": comps, "cut_tests": cut_tests,
                "directions": {f: rec[f] for f in RECORDED}, "zw_fallback": rec.get("zw_fallback", [])})
    return out


def _f(x, fmt="+.2f") -> str:
    return "n/a" if x is None else format(x, fmt)


def _cut_line(c: Mapping) -> str:
    return (f"{c['change']} against {c['base']} on both cuts: {c['test']} ({c['reason']}); "
            + "; ".join(f"{k} {x['mean']:+.2f} points, ahead {x['ahead']} of {x['of']}, per seed "
                        f"{x['range'][0]:+.2f} to {x['range'][1]:+.2f}" for k, x in c["cuts"].items()))


def t53_text(res: Mapping) -> str:
    """The reading as text; no cell count is printed."""
    r = res["reading"]
    lines = [res["label"], "Choices the brief left open:"]
    lines.extend(f"  - {x}" for x in CHOICES)
    if not r["valid"]:
        lines.append("Outcome: VOID")
        lines.append("Void: " + "; ".join(r["validity"]["reasons"]))
    else:
        for form, s in r["statements"].items():
            lines.append(f"{form} (slope-matched): {s['statement']}")
            for c in s["comparisons"].values():
                lines.append("    " + _cut_line(c))
            lines.append("    library truths (single-task test against the pooled arm of the same fraction): "
                         + "; ".join(f"{k}: " + ", ".join(f"{a} {lab}" for a, lab in ls_.items())
                                     for k, ls_ in s["library_labels"].items()))
            gs = s["gene_statement"] or "no statement (the form does not read better)"
            lines.append(f"    {form}'s genes: {gs}"
                         f"; {_cut_line(s['gene_comparison'])}; library truths "
                         + ", ".join(f"{k} {lab}" for k, lab in s["gene_library_labels"].items()))
        lines.append("Raw arms (no statement; the slope and what size alone would give are reported, read by nothing):")
        for x in r["raw_arms"].values():
            lines.append("    " + _cut_line(x["comparison"]) + "; library truths "
                         + ", ".join(f"{k} {lab}" for k, lab in x["library_labels"].items()))
            for k, per in x["size_alone"].items():
                lines.append(f"      {k}: " + "; ".join(
                    f"{g} slope {_f(v['slope'], '.3f')} -> "
                    + (_f(v["what_size_alone_would_give"]) if v["why_not"] is None else v["why_not"])
                    for g, v in per.items()))
        if r["zw_fallback"]:
            lines.append(f"zw keeps the entry's pooled direction for {r['zw_fallback']} (more than half the E9.5 "
                         "weight in zones under the minimum at E14.5)")
        for form, rec in r["directions"].items():
            lines.append(f"{form}'s directions on the entry's (moved genes): " + "; ".join(
                f"{g}: zones {v['zones'] or ['pooled']}, slope {_f(v['slope'], '.6f')}, pearson "
                f"{_f(v['pearson'], '.3f')}, squared size ratio {_f(v['size_ratio'], '.3f')}"
                + (f", dropped {v['dropped']} holding {v['dropped_share']:.3f} of E9.5" if v.get("dropped") else "")
                + (f", divided by {v['divisor']:.4f}" if "divisor" in v else "")
                for g, v in _groups(rec).items()))
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
    """Read the truths present in out/tasks against the committed T46 and T44 files; write t53.{json,txt}."""
    gv.refuse_existing(out, (f"{STEM}.json", f"{STEM}.txt"))
    started = gv.tree_state()
    files = COMMITTED if committed_files is None else committed_files
    present = {k: out / "tasks" / TASK_FILE.format(truth=k) for k in TRUTHS
               if (out / "tasks" / TASK_FILE.format(truth=k)).exists()}
    results = {k: json.loads(p.read_text()) for k, p in present.items()}
    found = {k: Path(p) for k, p in files.items() if Path(p).exists()}
    committed = {k: json.loads(p.read_text())["rows"] for k, p in found.items()}
    res = {"label": LABEL_T53, "choices": list(CHOICES), "reading": t53_reading(results, committed),
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
