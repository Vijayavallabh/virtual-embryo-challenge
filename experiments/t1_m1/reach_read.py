#!/usr/bin/env python
"""Try T60 (experiments/harness/TRIES.md): the reachability read. How much of a perfect per-type expression shift can
ANY size of the entry's external direction (Feng's E9.5 -> E14.5 change, M1's) reach on a one-day task, and how much
lies outside it?

Readout O1 showed that at an oracle mix a perfect per-type shift adds several board points on a one-day task
(shift_1 over oracle). Tries T47 and T51 measured how far along Feng's direction a released day moves; they did not
measure what lies outside the direction (decision log 2026-11-18 (t)). This read puts a bound on that: per type, the
part of the oracle shift that lies along the direction, at the size that fits it best, scored beside the whole oracle
shift.

Geometry. Source (and floor) E9.5, as rows T28 and T46. Each E10.5 library L is split once into halves A and B by T27's
rule (`split_library`: one permutation of the library's rows at SPLIT_SEED, the first half is A). The truth is half A
cut to the chamber categories, every cell as released (contamination 0, as row T28). Every arm starts from one draw per
seed of E9.5's real cells at half A's own chamber mix (`vec.t1_scope.reweight`, the oracle mix of rows O1 and T28).
The oracle shift is fitted on half B (no cell of the truth), E9.5 -> half B, per chamber category, at full size
(damping 1, O1's shift_1 route: `vec.t1_scope.momentum_shift` with the pre-registered held list, stream [seed, 1]).
Library 3 (the dissected heart the uploaded file is drawn from) is the reading; libraries 1 and 2 are reported beside
it. 10 seeds, scored with the board emulator as `experiments/t1_scope/pseudo_task.py` scores (T1:val's genes and
limits, 10% working copies, computed anchors, build main), every arm against one seed's anchors.

Arms, all on the same drawn cells:
  copy_last     a uniform draw of E9.5 (reference, not counted)
  oracle_mix    the base at the truth's mix, no expression change (reference, not counted)
  oracle_shift  oracle_mix plus the per-type shift fitted on half B at full size (counted)
  proj_shift    oracle_mix plus, per type, the least-squares projection of the oracle shift's per-gene log factor onto
                Feng's direction for that type, over the genes on which M1 moves that type (finite direction, not in
                R2's held set); every other gene unmoved. Applied as oracle_shift is applied (the same call, damping,
                pseudocount, depth rule and stream, its factor replaced through momentum_shift's logf_by_type), with
                one difference: the call's held list is M1's (R2's), so that R2's genes are left exactly as they are,
                outside the depth constant too, as external_shift leaves them; oracle_shift keeps pseudo_task's list. A
                type whose cells have no Feng group (R1's held types) gets no change. The most any size of the entry's
                direction can reach for that type on this task: an upper bound that uses the truth to choose the size
                (counted)
  proj_perm     the control: proj_shift with each type's direction permuted across the moved genes within deciles of
                the type's mean CP10k at E9.5 (`vec.harness.t3_decile_permutation`, one permutation per type at
                [PERM_SEED, the type's index in the config's category order]), projected and applied the same way: the
                same values and gene set, so the same size profile; its own fitted size (counted)
  entry_step    oracle_mix plus M1 at the entry's practice step (0.2 of Feng's bracket) with R1, R2 and g_mom's
                factors (E8.5 -> E9.5, `vec.harness.g_mom_pair_factors`): the arm m1_r_gmom's form on this base
                (`gmom_combo.moved_by_feng`, stream [seed, 2]; g_mom stream [seed, 5]; reported only, counted)

The guard, before any arm of the split is built: on row T28's own task (the whole library L cut to the chamber
categories, the same E9.5 source, E8.5 as momentum start) this script rebuilds T28's copy_last, oracle and momentum_0.5
with the same functions and streams as T28's script and scores them with its own scoring path. Each must match T28's
committed scores (experiments/t1_scope/results/t28_E9.5_E10.5_lib<L>/scores.csv): ten-seed mean within GUARD_MEAN_TOL,
every seed within GUARD_SEED_TOL. A library whose guard fails writes its task file with the guard and no split arm, and
the run exits non-zero.

  run    One library: the guard, then the split arms; writes results/<out>/tasks/lib<L>.json. Refused unless
         experiments/harness/TRIES.md holds row T60, and on an existing task file. Prints no score: per seed only its
         time, and the guard's verdict with its largest difference from T28's committed scores.
  read   T60's reading of the three task files (rule below, fixed in the row before the run); writes t60.{json,txt}.

The reading, on library 3, paired by seed. A = oracle_shift - proj_shift ("what a new direction could buy at most on
this task"); C = proj_shift - proj_perm ("what the direction's identity reaches beyond a random direction of the same
size profile"); B = proj_shift - entry_step, reported only (a matter of step size, a closed lever). Statement (1) "A new
direction has room": true if A's mean is above A_BAR with at least MIN_AHEAD of 10 seeds above 0. Statement (2) "The
direction's identity matters": true if C's mean is above C_BAR with at least MIN_AHEAD of 10 seeds above 0. Void unless
library 3's file holds the ten seeds with every arm, its guard holds, its tree was clean at one commit, every row's
anchors are uninverted, the coefficients are the same on every seed, oracle_shift differs from oracle_mix and proj_perm
from proj_shift on at least one seed. Libraries 1 and 2 are read the same way and reported; a void one is named.

NOT an adoption result: the oracle shift is fitted on the truth's own library (in-sample by construction), the task is
one day and the test two, the projection uses the truth to choose the size, Feng's E9.5 cells lie inside this task's
external window, and one truth is one fold of evidence. The result can only go to the team captain.

Usage: on the DGX campaign checkout, repo root, the three libraries at one commit, committing nothing until the three
task files exist (docs/COMPUTE.md):
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 \
      uv run --frozen --no-sync python -I -B experiments/t1_m1/reach_read.py run --out t60_reach --library 3
  (the same with --library 1 and --library 2), then, at the commit that holds the three task files,
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 \
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
REFERENCE_ARMS = ("copy_last", "oracle_mix")
COUNTED_ARMS = ("oracle_shift", "proj_shift", "proj_perm", "entry_step")
ARMS = (*REFERENCE_ARMS, *COUNTED_ARMS)
A_BAR = 1.5
C_BAR = 0.3
MIN_AHEAD = 8
STATEMENT_1 = "A new direction has room"
STATEMENT_2 = "The direction's identity matters"
READINGS = {"A": ("oracle_shift", "proj_shift"), "C": ("proj_shift", "proj_perm"), "B": ("proj_shift", "entry_step")}
LABEL = ("try T60, a labelled in-sample diagnostic, never an adoption result: how much of a perfect per-type shift any "
         "size of Feng's direction can reach on a one-day task (E9.5 -> half of an E10.5 library, the oracle shift "
         "fitted on the library's other half); the projection uses the truth to choose the size, the task is one day "
         "and the test two, and one truth is one fold of evidence; the result can only go to the team captain")
LIMITS = (
    "The oracle shift is fitted on the truth's own library (its other half): in-sample by construction.",
    "The projection's size is chosen with the truth: proj_shift is an upper bound for any size of the direction, never "
    "an estimate of what a step rule would deliver.",
    "The size is the least-squares fit to the oracle's log factor, not the size that scores best: rows T46 and T50 "
    "found the score along the direction still rising past the share the released day was measured to hold, so a "
    "better-scoring size along the direction can exist, and A can overstate the room outside it by that much. The "
    "coefficient is in the units of M1's fraction (0.2 at the practice step) and is reported beside it.",
    "The task is one day (E9.5 -> E10.5) and the test two (E10.5 -> E12.5): the horizon ratio is 2 on every counting "
    "task, so nothing here is a test figure (rule 9, J6).",
    "One truth is one fold of evidence: library 3 is one released sample, halved; libraries 1 and 2 are other samples "
    "of another kind and are reported only (rule 9, J5).",
    "A counts everything outside the direction's reach: genes M1 holds (R2), genes with no Feng value, types with no "
    "Feng group (R1, blood among them) and the part of each moved type's shift at right angles to the direction. It "
    "does not say that any source could reach it.",
    "C compares the best fit of the real direction with the best fit of one permutation of it: one draw of the "
    "control, the same values on the same genes, not the same applied size.",
    "B is a matter of step size, a closed lever; it is reported and decides nothing.",
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
    rec = {"coefficient": coef, "explained_share": float(coef * coef * dd / on) if on > 0 else float("nan"),
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
                held_m1: np.ndarray, start, order: Sequence[str]) -> tuple[dict, dict, dict]:
    """Per shifted type: proj_shift's and proj_perm's log factors and their records. A type with no Feng group gets
    zeros (no change) in both."""
    proj, perm, rec = {}, {}, {}
    for t, logf in logf_by_type.items():
        g = tgroups.get(t)
        if g is None:
            proj[t] = perm[t] = np.zeros_like(logf)
            rec[t] = {"group": None, "coefficient": None, "coefficient_perm": None}
            continue
        d = np.asarray(directions[g], dtype=np.float64)
        moved = moved_genes(d, held_m1)
        proj[t], r = project(logf, d, moved)
        dp = permuted_direction(d, type_basis(start, t), moved, type_seed(t, order))
        perm[t], rp = project(logf, dp, moved)
        rec[t] = {"group": g, **r, "coefficient_perm": rp["coefficient"],
                  "explained_share_perm": rp["explained_share"]}
    return proj, perm, rec


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
    """The five arms of one seed on `base` (oracle_mix), and the projection records. oracle_shift holds `held_pre`
    (pseudo_task's list); the projection arms hold `held_m1` (M1's, R2's), as external_shift does."""
    genes = list(base.var_names)
    held_m1 = sorted(set(held_m1))
    held_m1_mask = np.isin(genes, held_m1)
    shifted, _, logf = capture_shift(start, sib, base, held_pre, [seed, 1], **shift_kw)
    tgroups = type_groups(base, directions, groups)
    proj, perm, rec = projections(logf, tgroups, directions, held_m1_mask, start, order)
    arms = {"oracle_mix": base, "oracle_shift": shifted}
    arms["proj_shift"] = t1_scope.momentum_shift(start, sib, base, DAMPING, held_m1, seed=[seed, 1], logf_by_type=proj,
                                                 **shift_kw)[0]
    arms["proj_perm"] = t1_scope.momentum_shift(start, sib, base, DAMPING, held_m1, seed=[seed, 1], logf_by_type=perm,
                                                **shift_kw)[0]
    arms["entry_step"] = entry_arm(base, directions, groups, FRACTION, held_m1, held_pre, factors, seed)
    return arms, rec


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
        pick = np.sort(np.random.default_rng([s, 99]).choice(e95.n_obs, N_CELLS, replace=False))
        oracle, _ = t1_scope.reweight(e95, whole_shares, N_CELLS, s, within=within, fallback=fb["E9.5"])
        mom, _ = t1_scope.momentum_shift(e85, e95, oracle, MOMENTUM_DAMPING, held_pre, within=within,
                                         start_fallback=fb["E8.5"], ref_fallback=fb["E9.5"], seed=[s, 2],
                                         min_cells=MIN_CELLS, eps=EPS)
        arms = {"copy_last": pt.to_stage(e95[pick]), "oracle": pt.to_stage(oracle), "momentum_0.5": pt.to_stage(mom)}
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
                            "fraction": FRACTION, "board": BOARD, "frac": SUBSAMPLE_FRAC, "build": DEFAULT_BUILD}}
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


def t60_reading(results: Mapping[str, Mapping | None]) -> dict:
    """results: library -> task file (or None). Void unless library 3 is valid and the valid files share one commit."""
    libs = {}
    for lib in LIBRARIES:
        reasons = library_validity(results.get(lib))
        libs[lib] = {"valid": not reasons, "reasons": reasons}
    commits = {results[lib]["stamp"]["commit"] for lib in LIBRARIES if libs[lib]["valid"]}
    reasons = [f"library {READ_LIBRARY}: " + "; ".join(libs[READ_LIBRARY]["reasons"])] if not libs[READ_LIBRARY][
        "valid"] else []
    if len(commits) > 1:
        reasons.append(f"the valid libraries ran at {len(commits)} commits")
    out: dict = {"valid": not reasons, "reasons": reasons, "libraries": libs}
    if reasons:
        out["outcome"] = "VOID"
        return out
    for lib in LIBRARIES:
        if not libs[lib]["valid"]:
            continue
        rows = results[lib]["rows"]
        comps = {k: paired(rows, *pair) for k, pair in READINGS.items()}
        comps["oracle_shift_over_mix"] = paired(rows, "oracle_shift", "oracle_mix")
        comps["proj_shift_over_mix"] = paired(rows, "proj_shift", "oracle_mix")
        comps["proj_perm_over_mix"] = paired(rows, "proj_perm", "oracle_mix")
        comps["entry_step_over_mix"] = paired(rows, "entry_step", "oracle_mix")
        libs[lib].update(comparisons=comps, statements={STATEMENT_1: statement(comps["A"], A_BAR),
                                                        STATEMENT_2: statement(comps["C"], C_BAR)},
                         mean_scores={a: float(np.mean([r["arms"][a]["score"] for r in rows])) for a in ARMS},
                         projection=results[lib]["projection"])
        oracle_gain = comps["oracle_shift_over_mix"]["mean"]
        libs[lib]["reached_share"] = (comps["proj_shift_over_mix"]["mean"] / oracle_gain if oracle_gain > 0
                                      else float("nan"))
    out["statements"] = libs[READ_LIBRARY]["statements"]
    out["outcome"] = "; ".join(f"{k}: {'true' if v else 'false'}" for k, v in out["statements"].items())
    return out


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
    lines.append(f"Statements (library {READ_LIBRARY}, fixed in row T60): ({STATEMENT_1}) A above {A_BAR} with at "
                 f"least {MIN_AHEAD} of 10 seeds ahead; ({STATEMENT_2}) C above {C_BAR} with at least {MIN_AHEAD} "
                 "of 10")
    for lib, v in r["libraries"].items():
        tag = "the reading" if lib == READ_LIBRARY else "reported"
        if not v["valid"]:
            lines.append(f"library {lib} ({tag}): void: " + "; ".join(v["reasons"]))
            continue
        st = ", ".join(f"{k}: {'true' if b else 'false'}" for k, b in v["statements"].items())
        lines.append(f"library {lib} ({tag}): {st}")
        lines.append("    mean scores: " + ", ".join(f"{a} {x:.2f}" for a, x in v["mean_scores"].items()))
        for k in ("A", "C", "B", "oracle_shift_over_mix", "proj_shift_over_mix", "proj_perm_over_mix",
                  "entry_step_over_mix"):
            lines.append(_cmp_line(k, v["comparisons"][k]))
        lines.append("    share of the oracle shift's gain over the mix that proj_shift reaches: "
                     f"{v['reached_share']:.3f}")
        for t, p in v["projection"].items():
            if p.get("group") is None:
                lines.append(f"    type {t}: no Feng group (R1): unmoved in proj_shift and proj_perm")
                continue
            lines.append(f"    type {t} ({p['group']}): coefficient {p['coefficient']:.3f} (permuted "
                         f"{p['coefficient_perm']:.3f}); share of the oracle log shift's squared size on the moved "
                         f"genes along the direction {p['explained_share']:.3f} (permuted "
                         f"{p['explained_share_perm']:.3f}); squared size on the moved genes "
                         f"{p['size2_on_moved']:.2f}, off them {p['size2_off_moved']:.2f}")
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
