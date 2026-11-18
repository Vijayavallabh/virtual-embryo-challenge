#!/usr/bin/env python
"""Try T60 (experiments/harness/TRIES.md): the reachability read. How much of a perfect per-type expression shift can
the entry's external direction (Feng's E9.5 -> E14.5 change, M1's) reach on a one-day task at the best of several
sizes, and how much lies outside it?

Readout O1 showed that at an oracle mix a perfect per-type shift adds several board points on a one-day task
(shift_1 over oracle). Tries T47 and T51 measured how far along Feng's direction a released day moves; they did not
measure what lies outside the direction (decision log 2026-11-18 (t)). This read compares the oracle shift with the part
of it that lies along the direction, each at the best of the same five multiples of its own size.

Geometry. Source (and floor) E9.5, as rows T28 and T46. Each E10.5 library L is split once into halves A and B by T27's
rule (`split_library`: one permutation of the library's rows at SPLIT_SEED, the first half is A). The truth is half A
cut to the chamber categories, every cell as released (contamination 0, as row T28). Every arm starts from one draw per
seed of E9.5's real cells at half A's own chamber mix (`vec.t1_scope.reweight`, the oracle mix of rows O1 and T28).
The oracle shift is fitted on half B (no cell of the truth), E9.5 -> half B, per chamber category (O1's shift_1 route:
`vec.t1_scope.momentum_shift` with the pre-registered held list, stream [seed, 1]). Library 3 (the dissected heart the
uploaded file is drawn from) is the reading; libraries 1 and 2 (embryo portions, the kind E9.5 is) are its check, since
on library 3 the oracle shift also carries the change of sample kind. 10 seeds, scored with the board emulator as
`experiments/t1_scope/pseudo_task.py` scores (T1:val's genes and limits, 10% working copies, computed anchors, build
main), every arm against one seed's anchors.

Arms, all on the same drawn cells (k runs over MULTIPLES, 0.5, 1, 2, 3 and 4; the arm at k = 1 has the bare name):
  copy_last         a uniform draw of E9.5 (reference, not counted)
  oracle_mix        the base at the truth's mix, no expression change (reference, not counted)
  oracle_shift_x<k> oracle_mix plus the per-type shift fitted on half B at damping k (counted)
  proj_shift_x<k>   oracle_mix plus k times, per type, the least-squares projection of the oracle shift's per-gene log
                    factor (damping 1) onto Feng's direction for that type, over the genes on which M1 moves that type
                    (finite direction, not in R2's held set); every other gene unmoved. Applied as oracle_shift is
                    applied (the same call, damping 1, pseudocount, depth rule and stream, its factor replaced through
                    momentum_shift's logf_by_type), with one difference: the call's held list is M1's (R2's), so that
                    R2's genes are left exactly as they are, outside the depth constant too, as external_shift leaves
                    them; oracle_shift keeps pseudo_task's list. A type whose cells have no Feng group (R1's held types)
                    gets no change. Its size is chosen with the truth (counted)
  proj_perm         the control: proj_shift's size with each type's direction permuted across the moved genes within
                    deciles of the type's mean CP10k at E9.5 (`vec.harness.t3_decile_permutation`, one permutation per
                    type at [PERM_SEED, the type's index in the config's category order]): the same values on the same
                    genes at proj_shift's coefficient, so the same applied size (row T46's placebo) (counted)
  proj_perm_fit     the same permuted direction at its own least-squares coefficient (reported only, counted)
  entry_step        oracle_mix plus M1 at the entry's practice step (0.2 of Feng's bracket) with R1, R2 and g_mom's
                    factors (E8.5 -> E9.5, `vec.harness.g_mom_pair_factors`): the arm m1_r_gmom's form on this base
                    (`gmom_combo.moved_by_feng`, stream [seed, 2]; g_mom stream [seed, 5]; reported only, counted)

The guard, before any arm of the split is built: on row T28's own task (the whole library L cut to the chamber
categories, the same E9.5 source, E8.5 as momentum start) this script rebuilds T28's copy_last, oracle and momentum_0.5
(`guard_arms`, unit-tested equal to pseudo_task's own candidates) and scores them with its own scoring path. Each must
match T28's committed scores (experiments/t1_scope/results/t28_E9.5_E10.5_lib<L>/scores.csv): ten-seed mean within
GUARD_MEAN_TOL, every seed within GUARD_SEED_TOL. A library whose guard fails writes its task file with the guard and no
split arm, and the run exits non-zero. The run is refused before anything is loaded unless the tree is clean by
vec.provenance's rule.

  run    One library: the guard, then the split arms; writes results/<out>/tasks/lib<L>.json. Refused unless
         experiments/harness/TRIES.md holds row T60, on an unclean tree, and on an existing task file. Prints no score:
         per seed only its time, and the guard's verdict with its largest difference from T28's committed scores.
  read   T60's reading of the three task files (rule below, fixed in the row before the run); writes t60.{json,txt}.

The reading, paired by seed. On each library the best oracle multiple and the best projection multiple are the arms of
each family with the highest ten-seed mean (the smaller multiple on a tie). A = best oracle - best projection ("what a
new direction could buy at most on this task"). A library's test of statement (1) "A new direction has room" passes
when A's mean is above A_BAR with at least MIN_AHEAD of 10 seeds above 0. Library 3's state: "true" if the test passes
and the best projection multiple is below the largest; "not settled" if it passes with the best projection at the
largest multiple (still rising along the direction), or fails with the best oracle at the largest multiple (the oracle
still rising); "false" otherwise. Statement (1) is true only if library 3's state is true and libraries 1 and 2 are
valid and true as well; if library 3 is true and a portion library is not, it reads "library 3 only" (on library 3 the
room includes the change of sample kind, which the test's heart-to-heart pair does not have). C = proj_shift - proj_perm
("what the direction's identity reaches beyond a permutation of it at the same applied size"); statement (2) "The
direction's identity matters" is true if C's mean is above C_BAR with at least MIN_AHEAD of 10 seeds above 0, on
library 3. Reported: A at the least-squares size (oracle_shift - proj_shift), C against proj_perm_fit, B = proj_shift -
entry_step and every arm against oracle_mix. Void unless library 3's file holds the ten seeds with every arm, its guard
holds, its tree was clean at one commit, every row's anchors are uninverted, oracle_shift differs from oracle_mix and
proj_perm from proj_shift on at least one seed, and every moved type has a finite coefficient. A library 1 or 2 that is
void or ran at another commit than library 3 is named and reported, and the reading stands.

NOT an adoption result: the oracle shift is fitted on the truth's own library (in-sample by construction), the task is
one day and the test two, the sizes are chosen with the truth, Feng's E9.5 cells lie inside this task's external window,
and one truth is one fold of evidence. The result can only go to the team captain.

Usage: on the DGX campaign checkout, repo root, the three libraries at one commit, committing nothing until the three
task files exist (docs/COMPUTE.md):
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 \\
      uv run --frozen --no-sync python -I -B experiments/t1_m1/reach_read.py run --out t60_reach --library 3
  (the same with --library 1 and --library 2), then, at the commit that holds the three task files,
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 \\
      uv run --frozen --no-sync python -I -B experiments/t1_m1/reach_read.py read --out t60_reach
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from collections.abc import Mapping, Sequence
from pathlib import Path

import numpy as np

from vec import harness, t1_scope
from vec import t1_state_transport as ts
from vec.adoption import SEEDS
from vec.boards import GROUPS, RANKED_METRICS, get_board
from vec.core import counts_of, mean_cp10k
from vec.evaluate import DEFAULT_BUILD, SUBSAMPLE_FRAC, score_board
from vec.prep import load_stage

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]  # this checkout (vec.boards.REPO_ROOT is the installed package's)
RESULTS = HERE / "results"
T28_RESULTS = ROOT / "experiments" / "t1_scope" / "results"
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "experiments" / "t1_scope"))
import eb_read as er  # noqa: E402
import feng_dirs  # noqa: E402
import gmom_combo as gc  # noqa: E402
import pseudo_task as pt  # noqa: E402

gv, hygiene, n12 = gc.gv, gc.hygiene, er.n12

ROW = "T60"
TASK = "T1"
BOARD = "T1:val"  # rows O1 and T28: the gene order and cell limits
FILES = {"E8.5": "T1/E8.5_RNA.h5ad", "E9.5": "T1/E9.5_RNA.h5ad", "E10.5": "T1/E10.5_RNA.h5ad"}
LIBRARIES = ("1", "2", "3")
READ_LIBRARY = "3"  # the dissected heart; libraries 1 and 2 are reported beside it
SPLIT_SEED = 20261112  # row T27's split of library 3, applied to each library
PERM_SEED = 20261119  # proj_perm: one permutation per type, [PERM_SEED, type index in the config's category order]
PERM_BINS = harness.LP_BINS  # deciles, K2b's bins
N_CELLS = 5000  # rows O1 and T28
DAMPING = 1.0  # the oracle shift at full size (O1's shift_1)
EPS = 0.5  # momentum_shift's pseudocount, as pseudo_task
MIN_CELLS = 20
CONTAMINATION_LABELS = ("Foregut", "Epithelium")  # pseudo_task's default; unused at contamination 0 but passed alike
FRACTION = feng_dirs.bracket_fraction(9.5, 10.5)  # 0.2: M1's step on a one-day task
MOMENTUM_DAMPING = 0.5  # row T28's momentum_0.5, a guard arm only
GUARD_ARMS = ("copy_last", "oracle", "momentum_0.5")
GUARD_MEAN_TOL = 0.1  # board points, ten-seed mean
GUARD_SEED_TOL = 0.5  # board points, every seed
MULTIPLES = (0.5, 1.0, 2.0, 3.0, 4.0)  # each family's sizes, multiples of its own (damping 1; the fitted coefficient)


def scaled(name: str, k: float) -> str:
    """The arm of family `name` at multiple k: the bare name at 1, name_x<k> otherwise."""
    return name if k == 1.0 else f"{name}_x{k:g}"


ORACLE_ARMS = tuple(scaled("oracle_shift", k) for k in MULTIPLES)
PROJ_ARMS = tuple(scaled("proj_shift", k) for k in MULTIPLES)
REFERENCE_ARMS = ("copy_last", "oracle_mix")
COUNTED_ARMS = (*ORACLE_ARMS, *PROJ_ARMS, "proj_perm", "proj_perm_fit", "entry_step")
ARMS = (*REFERENCE_ARMS, *COUNTED_ARMS)
A_BAR = 1.5
C_BAR = 0.3
MIN_AHEAD = 8
STATEMENT_1 = "A new direction has room"
STATEMENT_2 = "The direction's identity matters"
PORTION_LIBRARIES = tuple(x for x in LIBRARIES if x != READ_LIBRARY)
# A itself is the best oracle multiple against the best projection multiple (best_arm); these are fixed pairs
READINGS = {"C": ("proj_shift", "proj_perm"), "A_least_squares": ("oracle_shift", "proj_shift"),
            "C_fit": ("proj_shift", "proj_perm_fit"), "B": ("proj_shift", "entry_step")}
LABEL = ("try T60, a labelled in-sample diagnostic, never an adoption result: how much of a perfect per-type shift "
         "Feng's direction can reach on a one-day task at the best of five sizes (E9.5 -> half of an E10.5 library, "
         "the oracle shift fitted on the library's other half); the sizes are chosen with the truth, the task is one "
         "day and the test two, and one truth is one fold of evidence; the result can only go to the team captain")
LIMITS = (
    "The oracle shift is fitted on the truth's own library (its other half): in-sample by construction. It carries "
    "that library's own signature (dissection, capture, ambient RNA), which no developmental direction could predict, "
    "so the room A measures can be non-developmental.",
    "Both families' sizes are chosen with the truth (the best of five multiples of each): the arms are optimistic "
    "bounds, never an estimate of what a step rule would deliver. The best of five sizes is not the best of every "
    "size: when the best projection is the largest multiple a passing test reads 'not settled', and when the best "
    "oracle is the largest a failing test reads 'not settled'.",
    "The projection's coefficient is the least-squares fit to the oracle's log factor over every moved gene, most of "
    "which the board barely scores; the multiples, not that fit, set the size A uses. The coefficient is in the units "
    "of M1's fraction (0.2 at the practice step) and is reported beside it; no size read here is a step setting.",
    "The task is one day (E9.5 -> E10.5) and the test two (E10.5 -> E12.5): the horizon ratio is 2 on every counting "
    "task, so nothing here is a test figure (rule 9, J6). The room outside the direction may be larger at two days.",
    "One truth is one fold of evidence: library 3 is one released sample, halved (rule 9, J5).",
    "On library 3 the oracle shift runs from E9.5's embryo portions to a dissected heart, so A there includes the "
    "change of sample kind (row T27: +8.24 points between the kinds), which the test's heart-to-heart pair does not "
    "have; statement (1) therefore needs the portion libraries 1 and 2 as well.",
    "The oracle is fitted on half a library: noise in its small types (epicardium, blood) lowers its score while the "
    "projection's one coefficient per type averages that noise away, so A can be understated; the multiple 0.5 covers "
    "part of this.",
    "A counts everything outside the direction's reach: genes M1 holds (R2), genes with no Feng value, types with no "
    "Feng group (R1, blood among them) and the part of each moved type's shift at right angles to the direction. Any "
    "source applied as M1 is would hold R2's genes and blood too. A does not say that any source could reach the rest.",
    "C compares the real direction with one permutation of it at the same applied size: one draw of the control.",
    "B mixes step size, g_mom's within-type factors and the route of application (per label through external_shift "
    "against per category through momentum_shift); it is reported and decides nothing.",
)


# ---------------------------------------------------------------------------------------------------------------------
# The geometry


def split_library(obs_names, library: str, seed: int = SPLIT_SEED) -> tuple[np.ndarray, np.ndarray]:
    """Row T27's split, for any library: the library's rows (in file order), permuted once at `seed`; the first half
    (rounded down) is A, the rest B; each returned sorted."""
    rows = np.flatnonzero(pt.library_of(obs_names) == library)
    perm = np.random.default_rng(seed).permutation(rows)
    return np.sort(perm[: len(perm) // 2]), np.sort(perm[len(perm) // 2:])


def chamber_rows(stage, rows: np.ndarray, chamber: Sequence[str]) -> np.ndarray:
    cats = stage.obs[t1_scope.CATEGORY_KEY].astype(str).to_numpy()
    return rows[np.isin(cats[rows], list(chamber))]


# ---------------------------------------------------------------------------------------------------------------------
# The projection


def type_groups(target_like, directions: Mapping, groups: Mapping[str, str],
                key: str = t1_scope.CATEGORY_KEY) -> dict[str, str | None]:
    """Each type (obs[key]) of `target_like` -> the Feng group whose direction its cells take under M1, or None when
    none of its cells has one (R1's held types). A type whose cells map to more than one group with a direction, or
    only some of them to one, is refused: momentum_shift moves a type by one factor."""
    types = target_like.obs[key].astype(str).to_numpy()
    labels = target_like.obs[t1_scope.LABEL_KEY].astype(str).to_numpy()
    out = {}
    for t in sorted(set(types)):
        gs = {groups.get(lab) if groups.get(lab) in directions else None for lab in set(labels[types == t])}
        if len(gs) != 1:
            raise ValueError(f"type {t!r}: its cells take {sorted(map(str, gs))} under M1, not one direction")
        out[t] = gs.pop()
    return out


def moved_genes(direction: np.ndarray, held: np.ndarray) -> np.ndarray:
    """The genes on which M1 moves a type: a finite direction value and not held (external_shift's zeros)."""
    return np.isfinite(np.asarray(direction, dtype=np.float64)) & ~np.asarray(held, dtype=bool)


def project(logf: np.ndarray, direction: np.ndarray, moved: np.ndarray) -> tuple[np.ndarray, dict]:
    """The least-squares multiple of `direction` that fits `logf` over the `moved` genes, as a full-length log factor
    (that multiple times the direction on the moved genes, 0 elsewhere), and a record: the coefficient, the share of
    logf's squared size on the moved genes that the fit explains, and logf's squared size on and off the moved genes."""
    logf = np.asarray(logf, dtype=np.float64)
    d = np.where(moved, np.nan_to_num(np.asarray(direction, dtype=np.float64), nan=0.0), 0.0)
    dd = float(d[moved] @ d[moved])
    coef = float(logf[moved] @ d[moved] / dd) if dd > 0 else 0.0
    out = coef * d
    on = float(logf[moved] @ logf[moved])
    dc, lc = d[moved] - d[moved].mean(), logf[moved] - logf[moved].mean()  # with a free constant (reported only)
    rec = {"coefficient": coef, "explained_share": float(coef * coef * dd / on) if on > 0 else float("nan"),
           "coefficient_centred": float(lc @ dc / (dc @ dc)) if float(dc @ dc) > 0 else 0.0,
           "size2_on_moved": on, "size2_off_moved": float(logf[~moved] @ logf[~moved]),
           "n_moved_genes": int(moved.sum())}
    return out, rec


def type_basis(start, t: str, key: str = t1_scope.CATEGORY_KEY) -> np.ndarray:
    """The type's mean CP10k over its own cells at the source stage (the permuted control's deciles)."""
    rows = np.flatnonzero(start.obs[key].astype(str).to_numpy() == t)
    if not len(rows):
        raise ValueError(f"type {t!r} has no own cells at the source stage")
    return np.asarray(mean_cp10k(counts_of(start)[rows]), dtype=np.float64).ravel()


def permuted_direction(direction: np.ndarray, basis: np.ndarray, moved: np.ndarray, seed) -> np.ndarray:
    """The direction's values on the moved genes permuted within deciles of `basis` over those genes (K2b's
    procedure); every other gene keeps its value. The same multiset of values on the same genes."""
    d = np.asarray(direction, dtype=np.float64).copy()
    idx = np.flatnonzero(moved)
    d[idx] = harness.t3_decile_permutation(d[idx], np.asarray(basis, dtype=np.float64)[idx], seed, PERM_BINS)
    return d


def type_seed(t: str, order: Sequence[str]) -> list[int]:
    return [PERM_SEED, list(order).index(t)]


# ---------------------------------------------------------------------------------------------------------------------
# The arms


def capture_shift(start, ref, target_like, held: Sequence[str], seed, **kw) -> tuple[object, object, dict]:
    """momentum_shift(start -> ref) on target_like at DAMPING, exactly as pseudo_task's shift_1 (stream `seed`), with an
    identity filter that records each shifted type's log factor as it is applied (after the held zeros). Returns the
    population, the report and type -> log factor."""
    seen: list[np.ndarray] = []

    def keep(logf, counts, held_mask):
        seen.append(np.array(logf, dtype=np.float64))
        return logf, None

    pop, rep = t1_scope.momentum_shift(start, ref, target_like, DAMPING, held, seed=seed, logf_filter=keep, **kw)
    shifted = [r["type"] for r in rep.rows if r["status"] == "shifted"]
    if len(shifted) != len(seen):
        raise RuntimeError("the shift's filter saw a different number of types than it shifted")
    return pop, rep, dict(zip(shifted, seen, strict=True))


def projections(logf_by_type: Mapping[str, np.ndarray], tgroups: Mapping[str, str | None], directions: Mapping,
                held_m1: np.ndarray, start, order: Sequence[str]) -> tuple[dict, dict, dict, dict]:
    """Per shifted type: the log factors of proj_shift, proj_perm (the permuted direction at proj_shift's coefficient:
    the same applied size) and proj_perm_fit (the permuted direction at its own least-squares coefficient), and the
    records. A type with no Feng group gets zeros (no change) in all three."""
    proj, perm, perm_fit, rec = {}, {}, {}, {}
    for t, logf in logf_by_type.items():
        g = tgroups.get(t)
        if g is None:
            proj[t], perm[t], perm_fit[t] = np.zeros_like(logf), np.zeros_like(logf), np.zeros_like(logf)
            rec[t] = {"group": None, "coefficient": None, "coefficient_perm": None}
            continue
        d = np.asarray(directions[g], dtype=np.float64)
        moved = moved_genes(d, held_m1)
        proj[t], r = project(logf, d, moved)
        dp = permuted_direction(d, type_basis(start, t), moved, type_seed(t, order))
        perm_fit[t], rp = project(logf, dp, moved)
        perm[t] = r["coefficient"] * np.where(moved, np.nan_to_num(dp, nan=0.0), 0.0)
        rec[t] = {"group": g, **r, "coefficient_perm": rp["coefficient"],
                  "explained_share_perm": rp["explained_share"]}
    return proj, perm, perm_fit, rec


def entry_arm(base, directions: Mapping, groups: Mapping[str, str], fraction: float, held_m1: Sequence[str],
              held_factors: Sequence[str], factors, seed: int):
    """M1 with R1 and R2 at `fraction` on `base` (`gmom_combo.moved_by_feng`, stream [seed, 2]), then g_mom's factors as
    the last edit (stream [seed, 5], no rescale, R1's held labels unedited): `gmom_combo.m1_r_arms`'s m1_r_gmom on this
    base. Refused if g_mom would edit a label Feng's direction does not cover (gmom_hook's check)."""
    moved = gc.moved_by_feng(base, directions, fraction, groups, held_m1, [seed, 2])
    labels = moved.obs[t1_scope.LABEL_KEY].astype(str).to_numpy()
    no_edit = np.isin(labels, sorted(ts.r1_unbracketed_labels(np.unique(labels))))
    cats = moved.obs[t1_scope.CATEGORY_KEY].astype(str).to_numpy()
    covered = {lab for lab, g in groups.items() if g in directions}
    extra = set(labels[~no_edit & np.isin(cats, list(factors.categories))]) - covered
    if extra:
        raise ValueError(f"g_mom would edit labels Feng's direction does not cover: {sorted(extra)}")
    return ts.apply_state_factors(moved, factors, list(held_factors), [seed, 5], rescale=False, no_edit=no_edit)[0]


def split_arms(base, start, sib, directions: Mapping, groups: Mapping[str, str], held_pre: Sequence[str],
               held_m1: Sequence[str], factors, seed: int, order: Sequence[str],
               shift_kw: Mapping) -> tuple[dict, dict]:
    """Every split arm of one seed on `base` (oracle_mix), and the projection records. The oracle family is
    momentum_shift at damping k and holds `held_pre` (pseudo_task's list); the projection family applies k times the
    projection at damping 1 and holds `held_m1` (M1's, R2's), as external_shift does. Every shifted arm uses the same
    stream [seed, 1]."""
    genes = list(base.var_names)
    held_m1 = sorted(set(held_m1))
    held_m1_mask = np.isin(genes, held_m1)
    shifted, _, logf = capture_shift(start, sib, base, held_pre, [seed, 1], **shift_kw)
    tgroups = type_groups(base, directions, groups)
    proj, perm, perm_fit, rec = projections(logf, tgroups, directions, held_m1_mask, start, order)

    def along(factors_by_type):
        return t1_scope.momentum_shift(start, sib, base, DAMPING, held_m1, seed=[seed, 1],
                                       logf_by_type=factors_by_type, **shift_kw)[0]

    arms = {"oracle_mix": base}
    for k in MULTIPLES:
        arms[scaled("oracle_shift", k)] = shifted if k == DAMPING else t1_scope.momentum_shift(
            start, sib, base, k, held_pre, seed=[seed, 1], **shift_kw)[0]
    for k in MULTIPLES:
        arms[scaled("proj_shift", k)] = along({t: k * v for t, v in proj.items()})
    arms["proj_perm"] = along(perm)
    arms["proj_perm_fit"] = along(perm_fit)
    arms["entry_step"] = entry_arm(base, directions, groups, FRACTION, held_m1, held_pre, factors, seed)
    return arms, rec


def guard_arms(mstart, ref, whole_shares: Mapping[str, float], within: Mapping, fb: Mapping, held_pre: Sequence[str],
               seed: int, n_cells: int = N_CELLS) -> dict:
    """Row T28's three arms of one seed on the whole library's shares, as pseudo_task's `candidates` builds them:
    copy_last (stream [seed, 99]), oracle (reweight at `seed`) and momentum_0.5 (mstart -> ref on the oracle's cells,
    stream [seed, 2]). fb: stage name -> fallbacks for "E8.5" and "E9.5"."""
    pick = np.sort(np.random.default_rng([seed, 99]).choice(ref.n_obs, n_cells, replace=False))
    oracle, _ = t1_scope.reweight(ref, whole_shares, n_cells, seed, within=within, fallback=fb["E9.5"])
    mom, _ = t1_scope.momentum_shift(mstart, ref, oracle, MOMENTUM_DAMPING, held_pre, within=within,
                                     start_fallback=fb["E8.5"], ref_fallback=fb["E9.5"], seed=[seed, 2],
                                     min_cells=MIN_CELLS, eps=EPS)
    return {"copy_last": pt.to_stage(ref[pick]), "oracle": pt.to_stage(oracle), "momentum_0.5": pt.to_stage(mom)}


# ---------------------------------------------------------------------------------------------------------------------
# Scoring and the guard


def score_arms(board, arms: Mapping, truth, ref_stage, seed: int) -> dict:
    """Every arm against one seed's anchors (computed once, on the first arm's call, as pseudo_task computes them: the
    reference working copy as floor, the truth's other half as ceiling): score, group skills, skills, metrics."""
    out, anchors, inverted = {}, "computed", None
    for name, pred in arms.items():
        r = score_board(board, pred, truth, ref_stage, seed=seed, anchors=anchors, frac=SUBSAMPLE_FRAC,
                        build=DEFAULT_BUILD)
        if anchors == "computed":
            anchors, inverted = r.anchors, sorted(r.extras.get("inverted_anchors") or [])
        out[name] = {"score": float(r.score), "groups": {g: float(v) for g, v in r.group_skill.items()},
                     "skills": {m: float(v) for m, v in r.skills.items()},
                     "metrics": {m: float(r.metrics[m]) for m in RANKED_METRICS[TASK] if m in r.metrics}}
    return {"arms": out, "inverted_anchors": inverted}


def committed_scores(path: Path, arms: Sequence[str] = GUARD_ARMS) -> dict[str, dict[int, float]]:
    """arm -> seed -> score from a pseudo_task scores.csv."""
    out: dict[str, dict[int, float]] = {a: {} for a in arms}
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            if row["candidate"] in out:
                out[row["candidate"]][int(row["seed"])] = float(row["score"])
    return out


def guard(fresh: Mapping[str, Mapping[int, float]], committed: Mapping[str, Mapping[int, float]],
          seeds: Sequence[int] = SEEDS, mean_tol: float = GUARD_MEAN_TOL, seed_tol: float = GUARD_SEED_TOL) -> dict:
    """Per arm, the ten-seed mean difference and the largest per-seed difference of `fresh` from `committed`; holds
    when every arm has every seed on both sides, every |mean difference| <= mean_tol and every |seed difference| <=
    seed_tol."""
    arms, missing = {}, []
    for a in committed:
        have = [s for s in seeds if s in fresh.get(a, {}) and s in committed[a]]
        missing += [f"{a}:{s}" for s in seeds if s not in have]
        d = np.array([fresh[a][s] - committed[a][s] for s in have], dtype=np.float64)
        arms[a] = {"mean_difference": float(d.mean()) if d.size else float("nan"),
                   "largest_seed_difference": float(np.abs(d).max()) if d.size else float("nan")}
    holds = bool(committed) and not missing and all(
        abs(v["mean_difference"]) <= mean_tol and v["largest_seed_difference"] <= seed_tol for v in arms.values())
    return {"arms": arms, "missing": missing, "holds": holds, "mean_tol": mean_tol, "seed_tol": seed_tol}


def t28_path(library: str) -> Path:
    return T28_RESULTS / f"t28_E9.5_E10.5_lib{library}" / "scores.csv"


# ---------------------------------------------------------------------------------------------------------------------
# The run


def cmd_run(a) -> None:
    if ROW not in n12.tries_rows():
        sys.exit(f"T60 is refused: {ROW} has no row in experiments/harness/TRIES.md (a row comes before its run)")
    out = RESULTS / a.out / "tasks"
    name = f"lib{a.library}.json"
    gv.refuse_existing(out, (name,))
    started = gv.tree_state()
    if started["clean"] is not True:  # an unclean stamp would void the library after its tries: refuse before them
        sys.exit("T60 is refused: the tree is not clean by vec.provenance's rule (git status --untracked-files=all, "
                 "outside results/, runs/ and experiments/**/results/)")
    seeds = gv.parse_seeds(a.seeds)
    data = Path(a.data)
    board = get_board(BOARD)
    if not board.min_cells <= N_CELLS <= (board.max_cells or N_CELLS):
        sys.exit(f"{N_CELLS} cells lie outside {BOARD}'s limits")
    spec = t1_scope.load_spec()
    chamber = spec.in_scope("chamber")
    stage = {k: load_stage(data / f, lineage=False) for k, f in FILES.items()}
    for k, s in stage.items():
        if list(s.var_names) != board.genes:
            sys.exit(f"{k}: genes differ from {BOARD}'s gene order")
        t1_scope.classify(s, source=f"{k}_RNA", spec=spec)
    e85, e95, e105 = stage["E8.5"], stage["E9.5"], stage["E10.5"]
    within = {"surrounding": list(CONTAMINATION_LABELS)}
    fb = {k: spec.fallbacks(f"{k}_RNA") for k in FILES}
    held_pre = t1_scope.held_flat_genes()
    ref_stage = pt.to_stage(e95)
    t0 = time.time()

    # The guard: row T28's task and arms, rebuilt and scored here
    whole, _ = pt.make_truth(e105, a.library, chamber, 0.0, list(CONTAMINATION_LABELS), 0)
    whole_shares = t1_scope.load_shares(t1_scope.composition(whole), spec)
    whole_stage = pt.to_stage(whole)
    fresh: dict[str, dict[int, float]] = {g: {} for g in GUARD_ARMS}
    for s in seeds:
        arms = guard_arms(e85, e95, whole_shares, within, fb, held_pre, s)
        res = score_arms(board, arms, whole_stage, ref_stage, s)
        for g in GUARD_ARMS:
            fresh[g][s] = res["arms"][g]["score"]
    committed = committed_scores(t28_path(a.library))
    g = guard(fresh, committed, seeds)
    worst = max(v["largest_seed_difference"] for v in g["arms"].values())
    print(f"library {a.library}: guard {'holds' if g['holds'] else 'FAILS'} (largest per-seed difference from T28's "
          f"committed scores {worst:.3g}; {time.time() - t0:.0f} s)", flush=True)
    record = {"label": LABEL, "row": ROW, "library": a.library, "seeds": list(seeds), "guard": g,
              "constants": {"split_seed": SPLIT_SEED, "perm_seed": PERM_SEED, "perm_bins": PERM_BINS,
                            "n_cells": N_CELLS, "damping": DAMPING, "eps": EPS, "min_cells": MIN_CELLS,
                            "fraction": FRACTION, "multiples": list(MULTIPLES), "board": BOARD,
                            "frac": SUBSAMPLE_FRAC, "build": DEFAULT_BUILD}}
    inputs = [data / f for f in FILES.values()] + [t1_scope.HELD_FLAT_GENES, t28_path(a.library)]
    if not g["holds"]:
        record.update(rows=None, stamp=gv.stamp(started, inputs))
        out.mkdir(parents=True, exist_ok=True)
        (out / name).write_text(json.dumps(record, indent=1, default=float) + "\n")
        sys.exit(f"library {a.library}: the guard fails; no split arm is built and the library is void")

    # The split
    half_a, half_b = split_library(e105.obs_names, a.library)
    truth = e105[chamber_rows(e105, half_a, chamber)].copy()
    sib = e105[half_b].copy()
    shares = t1_scope.load_shares(t1_scope.composition(truth), spec)
    truth_stage = pt.to_stage(truth)
    directions, ext = feng_dirs.feng_directions(list(e95.var_names))
    groups = feng_dirs.label_groups(e95.obs[t1_scope.LABEL_KEY].astype(str).unique())
    union = hygiene.load_union()
    held_m1 = gc.held_feng(union, held_pre)
    factors = harness.g_mom_pair_factors(e95, "E9.5_RNA", e85, "E8.5_RNA")[0]
    order = list(spec.categories)
    shift_kw = {"within": within, "start_fallback": fb["E9.5"], "ref_fallback": fb["E10.5"], "min_cells": MIN_CELLS,
                "eps": EPS}
    rows, coef = [], None
    for s in seeds:
        t1 = time.time()
        base, _ = t1_scope.reweight(e95, shares, N_CELLS, s, within=within, fallback=fb["E9.5"])
        pick = np.sort(np.random.default_rng([s, 99]).choice(e95.n_obs, N_CELLS, replace=False))
        arms, rec = split_arms(base, e95, sib, directions, groups, held_pre, held_m1, factors, s, order, shift_kw)
        if coef is None:
            coef = rec
        elif json.dumps(rec, sort_keys=True, default=float) != json.dumps(coef, sort_keys=True, default=float):
            sys.exit(f"library {a.library} seed {s}: the projection records differ from the first seed's")
        stages = {"copy_last": pt.to_stage(e95[pick]), **{k: pt.to_stage(v) for k, v in arms.items()}}
        res = score_arms(board, stages, truth_stage, ref_stage, s)
        rows.append({"seed": s, **res, "seconds": round(time.time() - t1, 1)})
        print(f"library {a.library} seed {s}: {time.time() - t1:.0f} s", flush=True)
    record.update(rows=rows, projection=coef, external=ext, arms={"reference": list(REFERENCE_ARMS),
                                                                  "counted": list(COUNTED_ARMS)},
                  truth={"n_rows": int(truth.n_obs), "n_half_b": int(sib.n_obs)},
                  stamp=gv.stamp(started, [*inputs, feng_dirs.clean_stats_path(), ROOT / union["file"]]))
    out.mkdir(parents=True, exist_ok=True)
    (out / name).write_text(json.dumps(record, indent=1, default=float) + "\n")


# ---------------------------------------------------------------------------------------------------------------------
# The reading

GROUP_WEIGHTS = {g: w for g, w, _ in GROUPS[TASK]}
DE_GROUPS, STATE_GROUPS = ("de_recovery", "direction"), ("distribution", "covariance")


def paired(rows: Sequence[Mapping], change: str, base: str) -> dict:
    """change - base per seed: mean, per-seed sd (ddof 1), seeds above 0, range, and the board points of each group
    (100 x weight x the mean group-skill difference) with the two de groups' and the two state groups' sums."""
    d = np.array([r["arms"][change]["score"] - r["arms"][base]["score"] for r in rows], dtype=np.float64)
    parts = {g: float(100.0 * w * np.mean([r["arms"][change]["groups"][g] - r["arms"][base]["groups"][g]
                                           for r in rows])) for g, w in GROUP_WEIGHTS.items()}
    return {"change": change, "base": base, "mean": float(d.mean()), "sd": float(d.std(ddof=1)) if d.size > 1 else 0.0,
            "ahead": int((d > 0).sum()), "of": int(d.size), "range": [float(d.min()), float(d.max())],
            "per_seed": d.tolist(), "groups": parts, "de_part": sum(parts[g] for g in DE_GROUPS),
            "state_part": sum(parts[g] for g in STATE_GROUPS)}


def statement(comp: Mapping, bar: float, min_ahead: int = MIN_AHEAD) -> bool:
    return bool(comp["mean"] > bar and comp["ahead"] >= min_ahead)


def best_arm(rows: Sequence[Mapping], family: Sequence[str]) -> tuple[str, dict[str, float]]:
    """The arm of `family` (in multiple order) with the highest ten-seed mean score, the smaller multiple on a tie,
    and every arm's mean."""
    means = {a: float(np.mean([r["arms"][a]["score"] for r in rows])) for a in family}
    top = max(means.values())
    return next(a for a in family if means[a] == top), means


def room_state(test: bool, best_oracle: str, best_proj: str) -> str:
    """One library's state of statement (1) from its test (A above the bar with enough seeds ahead) and where each
    family's best size lies: a passing test with the projection still rising at its largest multiple, or a failing
    test with the oracle still rising at its largest, is not settled."""
    if test:
        return "not settled" if best_proj == PROJ_ARMS[-1] else "true"
    return "not settled" if best_oracle == ORACLE_ARMS[-1] else "false"


def library_validity(r: Mapping | None, seeds: Sequence[int] = SEEDS) -> list[str]:
    """The void conditions of one library's file."""
    if r is None:
        return ["no task file"]
    reasons = []
    if not r.get("guard", {}).get("holds"):
        reasons.append("its guard against T28's committed scores does not hold")
    rows = r.get("rows")
    if not rows:
        return reasons + ["no split rows"]
    got = [row["seed"] for row in rows]
    if sorted(got) != list(seeds):
        reasons.append(f"seeds {got}, not {list(seeds)}")
    absent = sorted({arm for row in rows for arm in ARMS if arm not in row["arms"]})
    if absent:
        reasons.append(f"arms {absent} missing")
        return reasons
    if r.get("stamp", {}).get("tree_clean") is not True:
        reasons.append("its tree was not clean at one commit")
    inverted = sorted({m for row in rows for m in (row.get("inverted_anchors") or [])})
    if inverted:
        reasons.append(f"inverted anchors on {inverted}")
    for change, base in (("oracle_shift", "oracle_mix"), ("proj_perm", "proj_shift")):
        if all(row["arms"][change]["score"] == row["arms"][base]["score"] for row in rows):
            reasons.append(f"{change} scores exactly as {base} on every seed")
    proj = r.get("projection") or {}
    bad = [t for t, v in proj.items() if v.get("group") is not None
           and not (isinstance(v.get("coefficient"), int | float) and np.isfinite(v["coefficient"]))]
    if bad:
        reasons.append(f"no finite projection coefficient for {bad}")
    return reasons


def library_reading(r: Mapping) -> dict:
    """One valid library: the best of each family, A (best against best), the fixed pairs, every arm against the mix,
    the test and state of statement (1) and statement (2)'s test."""
    rows = r["rows"]
    best_o, means_o = best_arm(rows, ORACLE_ARMS)
    best_p, means_p = best_arm(rows, PROJ_ARMS)
    comps = {"A": paired(rows, best_o, best_p), **{k: paired(rows, *pair) for k, pair in READINGS.items()}}
    for arm in COUNTED_ARMS:
        comps[f"{arm}_over_mix"] = paired(rows, arm, "oracle_mix")
    test = statement(comps["A"], A_BAR)
    oracle_gain = comps[f"{best_o}_over_mix"]["mean"]
    return {"best_oracle": best_o, "best_projection": best_p, "family_means": {**means_o, **means_p},
            "comparisons": comps, "test_1": test, "state_1": room_state(test, best_o, best_p),
            "test_2": statement(comps["C"], C_BAR),
            "mean_scores": {a: float(np.mean([x["arms"][a]["score"] for x in rows])) for a in ARMS},
            "reached_share": comps[f"{best_p}_over_mix"]["mean"] / oracle_gain if oracle_gain > 0 else float("nan"),
            "projection": r["projection"]}


def t60_reading(results: Mapping[str, Mapping | None]) -> dict:
    """results: library -> task file (or None). Void unless library 3 is valid. A library 1 or 2 that is void or ran at
    another commit than library 3 is named and reported; statement (1) cannot then be true."""
    libs = {}
    for lib in LIBRARIES:
        reasons = library_validity(results.get(lib))
        libs[lib] = {"valid": not reasons, "reasons": reasons}
    if not libs[READ_LIBRARY]["valid"]:
        return {"valid": False, "reasons": [f"library {READ_LIBRARY}: " + "; ".join(libs[READ_LIBRARY]["reasons"])],
                "libraries": libs, "outcome": "VOID"}
    commit = results[READ_LIBRARY]["stamp"]["commit"]
    for lib in PORTION_LIBRARIES:
        if libs[lib]["valid"] and results[lib]["stamp"]["commit"] != commit:
            libs[lib] = {"valid": False, "reasons": [f"ran at commit {results[lib]['stamp']['commit']}, not library "
                                                     f"{READ_LIBRARY}'s {commit}"]}
    for lib in LIBRARIES:
        if libs[lib]["valid"]:
            libs[lib].update(library_reading(results[lib]))
    state3 = libs[READ_LIBRARY]["state_1"]
    if state3 != "true":
        s1 = state3
    elif all(libs[lib]["valid"] and libs[lib]["state_1"] == "true" for lib in PORTION_LIBRARIES):
        s1 = "true"
    else:
        s1 = f"library {READ_LIBRARY} only"
    s2 = "true" if libs[READ_LIBRARY]["test_2"] else "false"
    statements = {STATEMENT_1: s1, STATEMENT_2: s2}
    return {"valid": True, "reasons": [], "libraries": libs, "statements": statements,
            "outcome": "; ".join(f"{k}: {v}" for k, v in statements.items())}


def _cmp_line(name: str, c: Mapping) -> str:
    g = ", ".join(f"{k} {v:+.2f}" for k, v in c["groups"].items())
    return (f"    {name} ({c['change']} - {c['base']}): {c['mean']:+.2f} points, sd {c['sd']:.2f}, ahead {c['ahead']} "
            f"of {c['of']}, per seed {c['range'][0]:+.2f} to {c['range'][1]:+.2f}; de terms {c['de_part']:+.2f}, "
            f"mmd_u and variogram {c['state_part']:+.2f} ({g})")


def t60_text(res: Mapping) -> str:
    r = res["reading"]
    lines = [res["label"], f"Outcome: {r['outcome']}"]
    if not r["valid"]:
        lines.append("Void: " + "; ".join(r["reasons"]))
    lines.append(f"Statements (fixed in row T60): ({STATEMENT_1}) A, the best of the oracle multiples against the best "
                 f"of the projection multiples, above {A_BAR} on its ten-seed mean with at least {MIN_AHEAD} of 10 "
                 f"seeds ahead, on library {READ_LIBRARY} and on libraries {' and '.join(PORTION_LIBRARIES)}, the best "
                 "projection below the largest multiple; 'not settled' when a family's best is its largest multiple "
                 f"on the side that could change the answer; 'library {READ_LIBRARY} only' when a portion library does "
                 f"not show it. ({STATEMENT_2}) C above {C_BAR} with at least {MIN_AHEAD} of 10 ahead, on library "
                 f"{READ_LIBRARY}")
    for lib, v in r["libraries"].items():
        tag = "the reading" if lib == READ_LIBRARY else "the check on the change of kind"
        if not v["valid"]:
            lines.append(f"library {lib} ({tag}): void: " + "; ".join(v["reasons"]))
            continue
        lines.append(f"library {lib} ({tag}): statement (1)'s test {'passes' if v['test_1'] else 'fails'}, state "
                     f"{v['state_1']}; statement (2)'s test {'passes' if v['test_2'] else 'fails'}; best oracle "
                     f"{v['best_oracle']}, best projection {v['best_projection']}")
        lines.append("    mean scores: " + ", ".join(f"{a} {x:.2f}" for a, x in v["mean_scores"].items()))
        for k in ("A", "A_least_squares", "C", "C_fit", "B"):
            lines.append(_cmp_line(k, v["comparisons"][k]))
        for arm in COUNTED_ARMS:
            lines.append(_cmp_line(f"{arm} over the mix", v["comparisons"][f"{arm}_over_mix"]))
        lines.append("    share of the best oracle's gain over the mix that the best projection reaches: "
                     f"{v['reached_share']:.3f}")
        for t, p in v["projection"].items():
            if p.get("group") is None:
                lines.append(f"    type {t}: no Feng group (R1): unmoved in the projection arms")
                continue
            lines.append(f"    type {t} ({p['group']}): coefficient {p['coefficient']:.3f} (with a free constant "
                         f"{p['coefficient_centred']:.3f}; permuted {p['coefficient_perm']:.3f}; M1's step "
                         f"{FRACTION:.1f}); share of the oracle log shift's squared size on the moved genes along the "
                         f"direction {p['explained_share']:.3f} (permuted {p['explained_share_perm']:.3f}); squared "
                         f"size on the moved genes {p['size2_on_moved']:.2f}, off them {p['size2_off_moved']:.2f}")
    lines.append("What the read cannot show (written before the run):")
    lines += [f"  - {x}" for x in LIMITS]
    lines.append(f"Commit {res['stamp']['commit']}, tree clean {res['stamp']['tree_clean']}, UTC {res['stamp']['utc']}")
    return "\n".join(lines) + "\n"


def read(out: Path) -> dict:
    gv.refuse_existing(out, ("t60.json", "t60.txt"))
    started = gv.tree_state()
    paths = {lib: out / "tasks" / f"lib{lib}.json" for lib in LIBRARIES}
    results = {lib: json.loads(p.read_text()) if p.exists() else None for lib, p in paths.items()}
    res = {"label": LABEL, "reading": t60_reading(results),
           "stamps": {lib: r["stamp"] for lib, r in results.items() if r is not None}}
    res["stamp"] = gv.stamp(started, [p for p in paths.values() if p.exists()])
    gv.write_outputs(out, "t60", res, t60_text(res))
    return res


def cmd_read(a) -> None:
    res = read(RESULTS / a.out)
    if res["reading"]["outcome"] == "VOID":
        sys.exit("T60's reading is void and goes to the captain")


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("run", help="one library of try T60")
    p.add_argument("--library", required=True, choices=list(LIBRARIES))
    p.add_argument("--out", required=True)
    p.add_argument("--seeds", default=",".join(map(str, SEEDS)))
    p.add_argument("--data", default=str(harness.DATA_DIR))
    p = sub.add_parser("read", help="try T60's reading")
    p.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    {"run": cmd_run, "read": cmd_read}[a.cmd](a)


if __name__ == "__main__":
    main()
