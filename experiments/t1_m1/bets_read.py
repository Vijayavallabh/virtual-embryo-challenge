#!/usr/bin/env python
"""Try T61 (experiments/harness/TRIES.md): released-data bets on T1's last file, read on the practice tasks.

T1's last file as built (row R8; docs/PLAN.md decision log 2026-11-19 (h), (i)) is the first upload's design with the
outside atlas's step taken out: library 3's cells at their own chamber mix, every cell at its E10.5 expression, then
g_mom's within-type factors on cardiomyocytes and endocardium. No practice run has scored that path. This script scores
its practice analogue, `none_gmom`, beside released-data variants of it, with no outside file in any arm (no Feng file,
no list fixed by counting Feng cells), on the harness's draws and scoring.

Tasks. The harness's two E9.5 -> E10.5 cuts of configs/pseudo_tasks.yaml (the chamber cut is primary; the whole cut is
the same truth in a wider scope), with the start stage E8.5, as row T26 ran them (bake protocol, build main, build seed
1,000,000 + seed, seeds 0-9). Beside them, row T60's three library tasks as that row ran them (source and floor E9.5;
each E10.5 library split once at 20261112; the truth is half A cut to the chamber categories; T1:val's genes and
limits; 5,000 cells; `pseudo_task.py`'s scoring). Every task scores one truth, E10.5 (rule 9: one truth).

The g_mom step of every arm with `_gmom` (and of none_gmom, gmom_s1, share_trend_* and advance_*) is the file's own:
the factors of `vec.harness.g_mom_pair_factors` fitted on the task's start and reference (E8.5 -> E9.5; at test E9.5 ->
E10.5), applied last through `vec.t1_state_transport.apply_state_factors` (stream [seed, 5], no rescale, the
pre-registered held genes) on the cells of `build_test.G_MOM_EDIT_CATEGORIES` (cm, endocardium; row R3b) whose label the
readout config gives a group (`build_test._labels_without_group`, the mask of row R8), as `build_test._t1_g_mom`'s hook
does under --feng-groups none (tests pin the two together).

Arms on each E10.5 cut (b0 is the reproduction gate, not counted; the other ten are counted):
  b0               B0: the chamber copy plus the damped (0.5) E8.5 -> E9.5 momentum on every cell (`harness._b0_t1`)
  copy             the chamber copy, no change (`harness._b0_t1(0.0)`'s cells)
  none_gmom        copy + g_mom: the practice analogue of the built last file
  b0_gmom          b0 + g_mom
  mom_d0.25_gmom   the momentum at damping 0.25 + g_mom
  mom_d1_gmom      the momentum at damping 1.0 + g_mom
  gmom_s1          copy + g_mom's factors finished at damping 1.0 (`t1_state_transport.finish_momentum`)
  share_trend_l<l> g_mom on the chamber copy drawn at own x (E9.5 / E8.5) ** l, renormalised (`vec.t1_bets`), l 0.5, 1
  advance_t<t>     g_mom on the chamber copy drawn inside each category by z-quintile quotas in proportion to exp(t z)
                   (`vec.t1_bets.advance_reweight`), t 0.5, 1
On each library task: the references copy_last and oracle_mix (row T60's, the gate, not counted); on T60's oracle base
none_gmom, b0, b0_gmom, mom_d0.25_gmom, mom_d1_gmom, gmom_s1, advance_t0.5 and advance_t1; on E9.5's own chamber mix
own_copy, own_none_gmom, share_trend_l0.5 and share_trend_l1 (a trend applied to the truth's own mix would be no rule);
twelve counted per library.

  run-cut  One E10.5 cut: the gate's arm and the ten arms; writes results/<out>/tasks/<cut>.json. Refused unless row
           T61 is in experiments/harness/TRIES.md, on an unclean tree, and on an existing task file. Prints no score.
  run-lib  One library task: the references and the twelve arms; writes results/<out>/tasks/lib<L>.json. Same
           refusals. Prints no score.
  read     T61's statements, applied to the task files before anything is printed (the row's rule, below); writes
           t61.json and t61.txt and prints the outcome lines only.

The reading (fixed in row T61 before the run). The gate: on each cut b0 within REPRO_TOL points of row T11's committed
b0 on every seed (experiments/harness/results/t11_g_mom/tasks/); on each library copy_last and oracle_mix within
REPRO_TOL of row T60's committed scores on every seed. A cut that fails its gate is void; a void chamber cut voids the
read. "Ahead of none_gmom": on the chamber cut the arm's ten-seed paired mean over none_gmom is above AHEAD_BAR with at
least MIN_AHEAD seeds above 0, and on the whole cut its paired mean is above 0. The output arm: the ahead arm with the
highest chamber-cut mean (ties: ARM_ORDER). Reading S (T1's first upload stands): the output arm. Reading D (it is
discounted): the output arm if on every valid library task (at least MIN_LIBRARIES of them) its paired mean over its
base's none_gmom is at least -AHEAD_BAR; otherwise none. No arm ahead: none in both readings. Named beside, deciding
nothing: whether the output is one of DEPTH_LINKED (then the best ahead arm outside them too), and each composition
arm's drawn median library size over its base's (docs/DATA_POLICY.md 2026-11-04 (f)(3): beyond 5% goes to the
captain).

NOT an adoption result: every task scores one truth (E10.5), the tasks span one day and the test two, the advance and
share_trend arms change composition (the post-exposure block: only the captain's ruling puts one in a file), and the
library tasks' expression arms sit on an oracle mix (a bound).

Usage: on a DGX campaign checkout, repo root, at the row's commit, committing nothing until the five task files exist:
  export NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
  uv run --frozen --no-sync python -I -B experiments/t1_m1/bets_read.py run-cut --out t61_bets --cut chamber
  uv run --frozen --no-sync python -I -B experiments/t1_m1/bets_read.py run-cut --out t61_bets --cut whole
  uv run --frozen --no-sync python -I -B experiments/t1_m1/bets_read.py run-lib --out t61_bets --library 3
  (and --library 1, --library 2), then, at the commit that holds the five task files,
  uv run --frozen --no-sync python -I -B experiments/t1_m1/bets_read.py read --out t61_bets
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

from vec import harness, t1_bets, t1_scope
from vec import t1_state_transport as ts
from vec.adoption import SEEDS
from vec.boards import GROUPS, get_board, load_boards
from vec.evaluate import Stage
from vec.io import dense, read_stage
from vec.prep import load_stage

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]  # this checkout (vec.boards.REPO_ROOT is the installed package's)
RESULTS = HERE / "results"
TRIES = ROOT / "experiments" / "harness" / "TRIES.md"
T11_TASKS = ROOT / "experiments" / "harness" / "results" / "t11_g_mom" / "tasks"
T60_TASKS = HERE / "results" / "t60_reach" / "tasks"
for _p in (HERE, ROOT / "experiments" / "t1_scope", ROOT / "experiments" / "builds",
           ROOT / "experiments" / "t1_hygiene", ROOT / "experiments" / "t1_strong_null"):
    sys.path.insert(0, str(_p))
import build_test as bt  # noqa: E402  (the file's g_mom mask and categories)
import g_mom_verify as gv  # noqa: E402
import pseudo_task as pt  # noqa: E402
import reach_read as rr  # noqa: E402  (row T60's geometry and scoring; nothing of it that reads Feng is called)

ROW = "T61"
ROW_ID = re.compile(r"[A-Z][0-9]+[a-z]?")
CUTS = {"chamber": "T1/E9.5->E10.5:chamber", "whole": "T1/E9.5->E10.5"}
CUT_FILES = {"chamber": "T1_E9.5_E10.5_chamber.json", "whole": "T1_E9.5_E10.5.json"}  # also T11's committed names
LIBRARIES = rr.LIBRARIES
REPRO_TOL = 0.02  # board points per seed (row T26's tolerance)
AHEAD_BAR = 0.3  # board points, the adoption rule's figure
MIN_AHEAD = 8
MIN_LIBRARIES = 2
B0_DAMPING = harness.B0_DAMPING  # 0.5
MOM_DAMPINGS = {"mom_d0.25_gmom": 0.25, "b0_gmom": B0_DAMPING, "mom_d1_gmom": 1.0}
GMOM_STRONG = 1.0  # gmom_s1: g_mom's factors finished at this damping in place of B0_DAMPING
TREND_ARMS = {f"share_trend_l{x:g}": x for x in t1_bets.LAMBDAS}
ADVANCE_ARMS = {f"advance_t{x:g}": x for x in t1_bets.TAUS}
COMPOSITION_ARMS = (*TREND_ARMS, *ADVANCE_ARMS)  # the post-exposure block: composition (J1a)
DEPTH_LINKED = ("gmom_s1", "mom_d1_gmom")  # larger sizes whose practice gain tracks the truth's depth (T46, T59)
LIBRARY_SHIFT_TOL = 0.05  # docs/DATA_POLICY.md 2026-11-04 (f)(3): a drawn median library shift beyond 5% -> captain
GMOM_STREAM = 5
# The candidate arms in the order that breaks a tie (the fewest changes to the built file first)
ARM_ORDER = ("copy", "gmom_s1", "mom_d0.25_gmom", "b0_gmom", "mom_d1_gmom", "b0", *TREND_ARMS, *ADVANCE_ARMS)
CUT_GATE = "b0"
CUT_ARMS = ("b0", "copy", "none_gmom", "b0_gmom", "mom_d0.25_gmom", "mom_d1_gmom", "gmom_s1", *TREND_ARMS,
            *ADVANCE_ARMS)
CUT_COUNTED = tuple(a for a in CUT_ARMS if a != CUT_GATE)
LIB_REFERENCES = ("copy_last", "oracle_mix")
LIB_ORACLE_ARMS = ("none_gmom", "b0", "b0_gmom", "mom_d0.25_gmom", "mom_d1_gmom", "gmom_s1", *ADVANCE_ARMS)
LIB_OWN_ARMS = ("own_copy", "own_none_gmom", *TREND_ARMS)
LIB_ARMS = (*LIB_REFERENCES, *LIB_ORACLE_ARMS, *LIB_OWN_ARMS)
LIB_COUNTED = (*LIB_ORACLE_ARMS, *LIB_OWN_ARMS)
# On a library task each candidate arm is read against the none_gmom of its own base; the copy is the oracle mix there
LIB_PAIR = {**{a: (a, "none_gmom") for a in ("b0", "b0_gmom", "mom_d0.25_gmom", "mom_d1_gmom", "gmom_s1",
                                              *ADVANCE_ARMS)},
            "copy": ("oracle_mix", "none_gmom"), **{a: (a, "own_none_gmom") for a in TREND_ARMS}}
TRIES_COUNT = len(CUT_COUNTED) * len(CUTS) + len(LIB_COUNTED) * len(LIBRARIES)  # 56
LABEL = ("try T61, a labelled in-sample diagnostic, never an adoption result: released-data variants of T1's last file "
         "(the first upload's design without the outside atlas's step) against its practice analogue none_gmom, on "
         "the E9.5 -> E10.5 cuts and row T60's library tasks; every task scores one truth (E10.5) over one day, the "
         "test spans two; the composition arms touch the post-exposure block; the result can only go to the captain")


# ---------------------------------------------------------------------------------------------------------------------
# The file's g_mom step


def file_hook(factors: ts.StateFactors, seed: int, held_flat: Sequence[str]):
    """The g_mom step of T1's last file (`build_test._t1_g_mom`'s hook with feng_direction False), on given factors:
    the cells of G_MOM_EDIT_CATEGORIES whose label the readout config gives a group are edited, every other cell keeps
    its counts; stream [seed, 5], no rescale, `held_flat` held. A category of the edit with no factors stops the run."""
    missing = sorted(set(bt.G_MOM_EDIT_CATEGORIES) - set(factors.categories))
    if missing:
        raise ValueError(f"g_mom has no factors for a category R3b edits: {missing}")

    def hook(pop):
        labels = pop.obs[t1_scope.LABEL_KEY].astype(str).to_numpy()
        cats = pop.obs[t1_scope.CATEGORY_KEY].astype(str).to_numpy()
        held_labels = bt._labels_without_group(np.unique(labels))
        no_edit = np.isin(labels, sorted(held_labels)) | ~np.isin(cats, bt.G_MOM_EDIT_CATEGORIES)
        return ts.apply_state_factors(pop, factors, list(held_flat), [seed, GMOM_STREAM], rescale=False,
                                      no_edit=no_edit)[0]

    return hook


def gmom_fits(ref, src: str, start, start_src: str) -> dict[str, ts.StateFactors]:
    """g_mom's factors on one pair of stages, fitted once: finished at B0's damping (the file's) and at GMOM_STRONG."""
    finished, raw = harness.g_mom_pair_factors(ref, src, start, start_src)
    return {"file": finished, "strong": ts.finish_momentum(raw, GMOM_STRONG)}


# ---------------------------------------------------------------------------------------------------------------------
# The arms on the harness's E10.5 cuts


def cut_arms(ctx: harness.BuildContext, fits: Mapping[str, ts.StateFactors]) -> tuple[dict[str, Stage], dict]:
    """Every arm of one seed on a harness cut (CUT_ARMS), each finished as the harness finishes a candidate, and the
    drawn populations' median library sizes over the copy's (docs/DATA_POLICY.md 2026-11-04 (f)(3); reported only)."""
    role = ctx.task.reference_role
    ref, src = harness._shallow(ctx.inputs[role]), ctx.task.role_stems[role]
    start, start_src = harness._shallow(ctx.inputs["start"]), ctx.task.role_stems["start"]
    spec = t1_scope.load_spec()
    t1_scope.classify(ref, source=src, spec=spec)
    t1_scope.classify(start, source=start_src, spec=spec)
    shares = harness.t1_chamber_shares(ref, src)
    n = harness.count_for(ctx.board, ref.n_obs)
    held = list(ctx.derived["t1_held_flat"])
    hook, strong = file_hook(fits["file"], ctx.seed, held), file_hook(fits["strong"], ctx.seed, held)
    momentum = {"start": start, "start_source": start_src}

    def build(sh, **kw) -> Stage:
        sub, _ = t1_scope.build_submission(ref, ctx.board, sh, n, ctx.seed, held_flat=held, ref_source=src, **kw)
        return harness.finish(ctx, ref, dense(sub.X), None)

    arms = {"b0": harness._b0_t1(B0_DAMPING)(ctx), "copy": build(shares), "none_gmom": build(shares, post_cells=hook)}
    for name in ("b0_gmom", "mom_d0.25_gmom", "mom_d1_gmom"):
        arms[name] = build(shares, damping=MOM_DAMPINGS[name], post_cells=hook, **momentum)
    arms["gmom_s1"] = build(shares, post_cells=strong)

    def drawn(sh):  # the draw build_submission makes at these shares (the same call), for its median library
        return t1_scope.reweight(ref, sh, n, ctx.seed, fallback=spec.fallbacks(src), spec=spec)[0]

    depth = {"copy": t1_bets.median_library(drawn(shares))}
    later, earlier = t1_bets.chamber_counts(ref, src, spec), t1_bets.chamber_counts(start, start_src, spec)
    for name, lam in TREND_ARMS.items():
        trended = t1_bets.trend_shares(shares, later, earlier, lam)[0]
        arms[name] = build(trended, post_cells=hook)
        depth[name] = t1_bets.median_library(drawn(trended))
    cats = [c for c, s in shares.items() if s > 0]
    axes, _ = t1_bets.momentum_axes(start, ref, cats, held, start_fallback=spec.fallbacks(start_src),
                                    later_fallback=spec.fallbacks(src))
    z = t1_bets.advance_scores(ref, axes)
    for name, tau in ADVANCE_ARMS.items():
        pop, _ = t1_bets.advance_reweight(ref, shares, n, ctx.seed, z, tau, fallback=spec.fallbacks(src), spec=spec)
        arms[name] = harness.finish(ctx, ref, dense(t1_bets.assemble(pop, ctx.board, hook).X), None)
        depth[name] = t1_bets.median_library(pop)
    return {a: arms[a] for a in CUT_ARMS}, {a: v / depth["copy"] for a, v in depth.items() if a != "copy"}


# ---------------------------------------------------------------------------------------------------------------------
# The arms on row T60's library tasks


def library_arms(e85, e95, half_shares: Mapping[str, float], fits: Mapping[str, ts.StateFactors], held: Sequence[str],
                 seed: int, within: Mapping, fb: Mapping[str, Mapping], z: np.ndarray,
                 n_cells: int = rr.N_CELLS) -> dict:
    """Every arm of one seed on a library task (LIB_ARMS but copy_last, which the caller draws), as AnnData, and each
    composition arm's drawn median library size over its base's (reported only). e85, e95: classified; half_shares:
    the truth's own chamber mix (row T60's oracle base); fb: stage name -> fallbacks."""
    hook, strong = file_hook(fits["file"], seed, held), file_hook(fits["strong"], seed, held)
    shift_kw = {"within": within, "start_fallback": fb["E8.5"], "ref_fallback": fb["E9.5"], "min_cells": rr.MIN_CELLS,
                "eps": rr.EPS}
    base, _ = t1_scope.reweight(e95, half_shares, n_cells, seed, within=within, fallback=fb["E9.5"])
    arms = {"oracle_mix": base, "none_gmom": hook(base)}
    for name, d in MOM_DAMPINGS.items():  # B0's damping is b0_gmom's: b0 is the same moved cells without g_mom
        moved, _ = t1_scope.momentum_shift(e85, e95, base, d, list(held), seed=[seed, 1], **shift_kw)
        arms[name] = hook(moved)
        if d == B0_DAMPING:
            arms["b0"] = moved
    arms["gmom_s1"] = strong(base)
    med = t1_bets.median_library
    depth = {}
    for name, tau in ADVANCE_ARMS.items():
        pop, _ = t1_bets.advance_reweight(e95, half_shares, n_cells, seed, z, tau, within=within, fallback=fb["E9.5"])
        arms[name] = hook(pop)
        depth[name] = med(pop) / med(base)
    own = harness.t1_chamber_shares(e95, "E9.5_RNA")
    own_copy, _ = t1_scope.reweight(e95, own, n_cells, seed, within=within, fallback=fb["E9.5"])
    arms["own_copy"], arms["own_none_gmom"] = own_copy, hook(own_copy)
    later, earlier = t1_bets.chamber_counts(e95, "E9.5_RNA"), t1_bets.chamber_counts(e85, "E8.5_RNA")
    for name, lam in TREND_ARMS.items():
        drawn, _ = t1_scope.reweight(e95, t1_bets.trend_shares(own, later, earlier, lam)[0], n_cells, seed,
                                     within=within, fallback=fb["E9.5"])
        arms[name] = hook(drawn)
        depth[name] = med(drawn) / med(own_copy)
    return {a: arms[a] for a in LIB_ARMS if a != "copy_last"}, depth


# ---------------------------------------------------------------------------------------------------------------------
# Runs


def tries_rows(path: Path = TRIES) -> set[str]:
    """The row IDs of experiments/harness/TRIES.md's tables (each table row's first word, where it is a row ID)."""
    words = (ln[2:].split("|", 1)[0].split() for ln in path.read_text(encoding="utf-8").splitlines() if ln[:2] == "| ")
    return {w[0] for w in words if w and ROW_ID.fullmatch(w[0])}


def refuse(out: Path, name: str) -> dict:
    """The refusals before anything is loaded: the row, an existing task file, an unclean tree. Returns the start
    state."""
    if ROW not in tries_rows():
        sys.exit(f"{ROW} is refused: it has no row in experiments/harness/TRIES.md (a row comes before its run)")
    gv.refuse_existing(out, (name,))
    started = gv.tree_state()
    if started["clean"] is not True:
        sys.exit(f"{ROW} is refused: the tree is not clean by vec.provenance's rule")
    return started


def cmd_run_cut(a) -> None:
    out = RESULTS / a.out / "tasks"
    started = refuse(out, CUT_FILES[a.cut])
    seeds = gv.parse_seeds(a.seeds)
    data = Path(a.data)
    t = {x.key: x for x in harness.load_tasks()["T1:test"]}[CUTS[a.cut]]
    boards = load_boards()
    genes = boards[t.panel_board].genes
    target, _ = (harness.filter_target(t, data, genes) if t.target_filter is not None else
                 (Stage.from_anndata(read_stage(data / t.target, genes), genes), None))
    ref_stage = Stage.from_anndata(read_stage(data / t.reference, genes), genes)
    board = harness.task_board(t, boards, target.n, ref_stage.n)
    inputs = {r: load_stage(data / f, lineage=False) for r, f in t.inputs.items()}
    preregistered = t1_scope.held_flat_genes()
    view, cboard = harness.TaskView.of(t), harness._candidate_board(board)
    ref_role = harness.base_role(t)
    fits = gmom_fits(inputs[ref_role], view.role_stems[ref_role], inputs["start"], view.role_stems["start"])
    rows = []
    for s in seeds:
        t0 = time.time()
        f = gv.floor_anchors(t, board, target, ref_stage, s)
        ctx = harness.BuildContext(task=view, board=cboard, inputs=inputs, seed=harness.BUILD_SEED_OFFSET + s,
                                   derived={"t1_held_flat": preregistered})
        preds, depth = cut_arms(ctx, fits)
        arms = gv.score_arms(board, target, ref_stage, s, f.anchors, preds)
        rows.append({"seed": s, "floor": float(f.score), "arms": arms, "median_library_over_base": depth,
                     "seconds": round(time.time() - t0, 1)})
        print(f"{t.key} seed {s}: {time.time() - t0:.0f} s", flush=True)  # no score
    res = {"label": LABEL, "row": ROW, "task": t.key, "seeds": list(seeds), "rows": rows,
           "constants": {"gmom_strong": GMOM_STRONG, "dampings": MOM_DAMPINGS, "lambdas": list(t1_bets.LAMBDAS),
                         "taus": list(t1_bets.TAUS), "edit_categories": list(bt.G_MOM_EDIT_CATEGORIES)},
           "stamp": gv.stamp(started, [data / t.target, data / t.reference, *[data / f for f in t.inputs.values()],
                                       t1_scope.HELD_FLAT_GENES, t1_scope.CHAMBER_CLASSES, ts.R1_CONFIG])}
    out.mkdir(parents=True, exist_ok=True)
    (out / CUT_FILES[a.cut]).write_text(json.dumps(res, indent=1, default=float) + "\n")


def cmd_run_lib(a) -> None:
    out = RESULTS / a.out / "tasks"
    name = f"lib{a.library}.json"
    started = refuse(out, name)
    seeds = gv.parse_seeds(a.seeds)
    data = Path(a.data)
    board = get_board(rr.BOARD)
    spec = t1_scope.load_spec()
    chamber = spec.in_scope("chamber")
    stage = {k: load_stage(data / f, lineage=False) for k, f in rr.FILES.items()}
    for k, s in stage.items():
        if list(s.var_names) != board.genes:
            sys.exit(f"{k}: genes differ from {rr.BOARD}'s gene order")
        t1_scope.classify(s, source=f"{k}_RNA", spec=spec)
    e85, e95, e105 = stage["E8.5"], stage["E9.5"], stage["E10.5"]
    within = {"surrounding": list(rr.CONTAMINATION_LABELS)}
    fb = {k: spec.fallbacks(f"{k}_RNA") for k in rr.FILES}
    held = t1_scope.held_flat_genes()
    ref_stage = pt.to_stage(e95)
    half_a, _ = rr.split_library(e105.obs_names, a.library)
    truth = e105[rr.chamber_rows(e105, half_a, chamber)].copy()
    shares = t1_scope.load_shares(t1_scope.composition(truth), spec)
    truth_stage = pt.to_stage(truth)
    fits = gmom_fits(e95, "E9.5_RNA", e85, "E8.5_RNA")
    cats = [c for c, s in harness.t1_chamber_shares(e95, "E9.5_RNA").items() if s > 0]
    axes, _ = t1_bets.momentum_axes(e85, e95, cats, held, within=within, start_fallback=fb["E8.5"],
                                    later_fallback=fb["E9.5"])
    z = t1_bets.advance_scores(e95, axes)
    rows = []
    for s in seeds:
        t0 = time.time()
        pick = np.sort(np.random.default_rng([s, 99]).choice(e95.n_obs, rr.N_CELLS, replace=False))
        arms, depth = library_arms(e85, e95, shares, fits, held, s, within, fb, z)
        stages = {"copy_last": pt.to_stage(e95[pick]), **{k: pt.to_stage(v) for k, v in arms.items()}}
        rows.append({"seed": s, **rr.score_arms(board, stages, truth_stage, ref_stage, s),
                     "median_library_over_base": depth, "seconds": round(time.time() - t0, 1)})
        print(f"library {a.library} seed {s}: {time.time() - t0:.0f} s", flush=True)  # no score
    res = {"label": LABEL, "row": ROW, "library": a.library, "seeds": list(seeds), "rows": rows,
           "constants": {"split_seed": rr.SPLIT_SEED, "n_cells": rr.N_CELLS, "board": rr.BOARD,
                         "gmom_strong": GMOM_STRONG, "dampings": MOM_DAMPINGS},
           "stamp": gv.stamp(started, [*(data / f for f in rr.FILES.values()), t1_scope.HELD_FLAT_GENES,
                                       t1_scope.CHAMBER_CLASSES, ts.R1_CONFIG])}
    out.mkdir(parents=True, exist_ok=True)
    (out / name).write_text(json.dumps(res, indent=1, default=float) + "\n")


# ---------------------------------------------------------------------------------------------------------------------
# The reading


GROUP_WEIGHTS = {g: w for g, w, _ in GROUPS["T1"]}


def score_of(row: Mapping, arm: str) -> float:
    return float(row["arms"][arm]["score"])


def paired(rows: Sequence[Mapping], change: str, base: str) -> dict:
    """change - base per seed, board points: mean, per-seed sd, seeds above 0, range, and the points of each scored
    group (100 x weight x the mean group-skill difference)."""
    d = np.array([score_of(r, change) - score_of(r, base) for r in rows], dtype=np.float64)
    parts = {g: float(100.0 * w * np.mean([r["arms"][change]["groups"][g] - r["arms"][base]["groups"][g]
                                           for r in rows])) for g, w in GROUP_WEIGHTS.items()}
    return {"change": change, "base": base, "mean": float(d.mean()),
            "sd": float(d.std(ddof=1)) if d.size > 1 else 0.0, "ahead": int((d > 0).sum()), "of": int(d.size),
            "range": [float(d.min()), float(d.max())], "groups": parts}


def gate(rows: Sequence[Mapping], committed: Mapping[int, Mapping[str, float]], arms: Sequence[str],
         seeds: Sequence[int] = SEEDS, tol: float = REPRO_TOL) -> dict:
    """Per arm, the largest per-seed |fresh - committed|; holds when every seed is on both sides and every difference
    is at most tol."""
    fresh = {int(r["seed"]): r for r in rows}
    missing = [s for s in seeds if s not in fresh or s not in committed]
    worst = {a: max((abs(score_of(fresh[s], a) - committed[s][a]) for s in seeds if s not in missing),
                    default=float("nan")) for a in arms}
    return {"largest_difference": worst, "missing_seeds": missing, "tolerance": tol,
            "holds": not missing and all(np.isfinite(v) and v <= tol for v in worst.values())}


def committed_t11(cut: str) -> dict[int, dict[str, float]]:
    rows = json.loads((T11_TASKS / CUT_FILES[cut]).read_text())["rows"]
    return {int(r["seed"]): {CUT_GATE: float(r["scores"][CUT_GATE])} for r in rows}


def committed_t60(library: str) -> dict[int, dict[str, float]]:
    rows = json.loads((T60_TASKS / f"lib{library}.json").read_text())["rows"]
    return {int(r["seed"]): {a: float(r["arms"][a]["score"]) for a in LIB_REFERENCES} for r in rows}


def validity(r: Mapping | None, arms: Sequence[str], seeds: Sequence[int] = SEEDS) -> list[str]:
    if r is None:
        return ["no task file"]
    rows = r.get("rows") or []
    reasons = []
    if sorted(int(x["seed"]) for x in rows) != list(seeds):
        reasons.append(f"seeds {[x['seed'] for x in rows]}, not {list(seeds)}")
    absent = sorted({a for x in rows for a in arms if a not in x["arms"]})
    if absent:
        return reasons + [f"arms {absent} missing"]
    if r.get("stamp", {}).get("tree_clean") is not True:
        reasons.append("its tree was not clean at one commit")
    if rows and all(x["arms"]["none_gmom"]["score"] == x["arms"]["copy" if "copy" in x["arms"] else "oracle_mix"]
                    ["score"] for x in rows):
        reasons.append("none_gmom scores exactly as its copy on every seed: g_mom's step did nothing")
    if rows and "inverted_anchors" in rows[0] and any(x.get("inverted_anchors") for x in rows):
        reasons.append("inverted anchors")
    return reasons


def cut_reading(r: Mapping, cut: str, committed: Mapping) -> dict:
    reasons = validity(r, CUT_ARMS)
    g = gate(r["rows"], committed, (CUT_GATE,)) if r is not None and not reasons else None
    if g is not None and not g["holds"]:
        reasons.append("its gate against row T11's committed b0 does not hold")
    if reasons:
        return {"valid": False, "reasons": reasons, "gate": g}
    rows = r["rows"]
    comps = {a: paired(rows, a, "none_gmom") for a in ARM_ORDER}
    reported = {"gmom_on_copy": paired(rows, "none_gmom", "copy"), "gmom_on_b0": paired(rows, "b0_gmom", "b0"),
                "none_gmom_over_floor": {"mean": float(np.mean([score_of(x, "none_gmom") - x["floor"] for x in rows]))}}
    return {"valid": True, "reasons": [], "gate": g, "commit": r["stamp"]["commit"], "against_none_gmom": comps,
            "reported": reported, "mean_scores": {a: float(np.mean([score_of(x, a) for x in rows])) for a in CUT_ARMS},
            "median_library_over_base": mean_depth(rows)}


def mean_depth(rows: Sequence[Mapping]) -> dict[str, float]:
    """Each composition arm's drawn median library size over its base's, mean over the seeds that record it."""
    got: dict[str, list[float]] = {}
    for x in rows:
        for a, v in (x.get("median_library_over_base") or {}).items():
            got.setdefault(a, []).append(float(v))
    return {a: float(np.mean(v)) for a, v in got.items()}


def library_reading(r: Mapping | None, library: str, committed: Mapping) -> dict:
    reasons = validity(r, LIB_ARMS)
    g = gate(r["rows"], committed, LIB_REFERENCES) if r is not None and not reasons else None
    if g is not None and not g["holds"]:
        reasons.append("its gate against row T60's committed copy_last and oracle_mix does not hold")
    if reasons:
        return {"valid": False, "reasons": reasons, "gate": g}
    rows = r["rows"]
    comps = {a: paired(rows, *LIB_PAIR[a]) for a in ARM_ORDER}
    return {"valid": True, "reasons": [], "gate": g, "commit": r["stamp"]["commit"], "against_none_gmom": comps,
            "mean_scores": {a: float(np.mean([score_of(x, a) for x in rows])) for a in LIB_ARMS},
            "median_library_over_base": mean_depth(rows)}


def ahead(chamber: Mapping, whole: Mapping, arm: str) -> bool:
    c, w = chamber["against_none_gmom"][arm], whole["against_none_gmom"][arm]
    return bool(c["mean"] > AHEAD_BAR and c["ahead"] >= MIN_AHEAD and w["mean"] > 0)


def t61_reading(cuts: Mapping[str, Mapping | None], libs: Mapping[str, Mapping | None],
                committed_cuts: Mapping[str, Mapping], committed_libs: Mapping[str, Mapping]) -> dict:
    """The statements of row T61, applied to the task files: the gates, the arms ahead of none_gmom, the output arm,
    and readings S and D."""
    cr = {c: cut_reading(cuts.get(c), c, committed_cuts[c]) for c in CUTS}
    lr = {lib: library_reading(libs.get(lib), lib, committed_libs[lib]) for lib in LIBRARIES}
    if not cr["chamber"]["valid"]:
        return {"outcome": "VOID", "cuts": cr, "libraries": lr, "reasons": ["the chamber cut is void"]}
    commit = cr["chamber"]["commit"]
    for v in (*cr.values(), *lr.values()):
        if v["valid"] and v["commit"] != commit:
            v.update(valid=False, reasons=[f"ran at commit {v['commit']}, not the chamber cut's {commit}"])
    ahead_arms = ([a for a in ARM_ORDER if ahead(cr["chamber"], cr["whole"], a)] if cr["whole"]["valid"] else [])
    best = (max(ahead_arms, key=lambda a: (cr["chamber"]["against_none_gmom"][a]["mean"], -ARM_ORDER.index(a)))
            if ahead_arms else None)
    valid_libs = [lib for lib in LIBRARIES if lr[lib]["valid"]]
    d_check = None
    if best is not None:
        lows = {lib: lr[lib]["against_none_gmom"][best]["mean"] for lib in valid_libs}
        d_check = {"libraries": lows, "enough": len(valid_libs) >= MIN_LIBRARIES,
                   "passes": len(valid_libs) >= MIN_LIBRARIES and all(v >= -AHEAD_BAR for v in lows.values())}
    out_s = best or "none"
    out_d = best if d_check and d_check["passes"] else "none"
    # The best ahead arm outside the two sizes whose practice gain is known to track the truth's depth (rows T46, T59;
    # decision log 2026-11-18 (f)): given beside the output, so the captain has a depth-free option to rule on
    free = [a for a in ahead_arms if a not in DEPTH_LINKED]
    depth_free = (max(free, key=lambda a: (cr["chamber"]["against_none_gmom"][a]["mean"], -ARM_ORDER.index(a)))
                  if free else "none")
    shift = cr["chamber"]["median_library_over_base"].get(best) if best else None
    ranges = {}
    for a in ARM_ORDER:
        vals = [cr[c]["against_none_gmom"][a]["mean"] for c in CUTS if cr[c]["valid"]]
        vals += [lr[lib]["against_none_gmom"][a]["mean"] for lib in valid_libs]
        ranges[a] = {"analogue_range": [min(vals), max(vals)],
                     "j6_range": [0.0, max(0.0, cr["chamber"]["against_none_gmom"][a]["mean"])]}
    return {"outcome": f"S: {out_s}; D: {out_d}", "reading_S": out_s, "reading_D": out_d, "ahead": ahead_arms,
            "depth_linked_output": best in DEPTH_LINKED if best else False, "depth_free_ahead": depth_free,
            "median_library_shift": None if shift is None else {"ratio": shift,
                                                                "beyond": bool(abs(shift - 1) > LIBRARY_SHIFT_TOL)},
            "d_check": d_check, "whole_cut_valid": cr["whole"]["valid"], "valid_libraries": valid_libs,
            "composition_arm": best in COMPOSITION_ARMS if best else False, "ranges": ranges, "cuts": cr,
            "libraries": lr}


def _line(name: str, c: Mapping) -> str:
    g = ", ".join(f"{k} {v:+.2f}" for k, v in c["groups"].items())
    return (f"    {name} ({c['change']} - {c['base']}): {c['mean']:+.2f} points, sd {c['sd']:.2f}, ahead "
            f"{c['ahead']} of {c['of']}, per seed {c['range'][0]:+.2f} to {c['range'][1]:+.2f} ({g})")


def outcome_lines(r: Mapping) -> list[str]:
    """The statements' results, the only lines printed to the terminal."""
    if r["outcome"] == "VOID":
        return ["Outcome: VOID (" + "; ".join(r["reasons"]) + ")"]
    lines = [f"Outcome: reading S (T1's first upload stands): {r['reading_S']}; reading D (it is discounted): "
             f"{r['reading_D']}", f"Arms ahead of none_gmom: {', '.join(r['ahead']) or 'none'}"]
    if not r["whole_cut_valid"]:
        lines.append("The whole cut is void: no arm can be ahead (its sign check cannot be made)")
    if r["reading_S"] != "none" and r["composition_arm"]:
        lines.append("The output arm changes composition: the post-exposure block; only the captain's ruling puts it "
                     "in a file (J1a: a range, never a point)")
    if r["depth_linked_output"]:
        lines.append("The output arm is a larger size whose practice gain tracks the truth's depth (T46, T59; decision "
                     "log 2026-11-18 (f)): released data alone or decline. The best ahead arm outside such sizes: "
                     f"{r['depth_free_ahead']}")
    if r["median_library_shift"] and r["median_library_shift"]["beyond"]:
        lines.append("The output arm's drawn cells shift the median library size beyond 5% of its base's "
                     "(docs/DATA_POLICY.md 2026-11-04 (f)(3)): the rule goes to the captain")
    return lines


def t61_text(res: Mapping) -> str:
    r = res["reading"]
    lines = [res["label"], *outcome_lines(r)]
    for c, v in r["cuts"].items():
        lines.append(f"cut {c}: " + ("valid" if v["valid"] else "void: " + "; ".join(v["reasons"])))
        if v.get("gate"):
            lines.append(f"    gate: largest per-seed difference from T11's b0 {v['gate']['largest_difference']}")
        if v["valid"]:
            lines.append("    mean scores: " + ", ".join(f"{a} {x:.2f}" for a, x in v["mean_scores"].items()))
            lines += [_line(a, cmp) for a, cmp in v["against_none_gmom"].items()]
            lines += [_line(k, cmp) for k, cmp in v["reported"].items() if "change" in cmp]
    for lib, v in r["libraries"].items():
        lines.append(f"library {lib}: " + ("valid" if v["valid"] else "void: " + "; ".join(v["reasons"])))
        if v["valid"]:
            lines.append("    mean scores: " + ", ".join(f"{a} {x:.2f}" for a, x in v["mean_scores"].items()))
            lines += [_line(a, cmp) for a, cmp in v["against_none_gmom"].items()]
    if r["outcome"] != "VOID":
        lines.append("Ranges (J1a: the analogue readings' paired means, cuts and valid libraries; J6: 0 to the chamber "
                     "cut's gain):")
        lines += [f"    {a}: analogue {v['analogue_range'][0]:+.2f} to {v['analogue_range'][1]:+.2f}; J6 "
                  f"{v['j6_range'][0]:+.2f} to {v['j6_range'][1]:+.2f}" for a, v in r["ranges"].items()]
    lines.append(f"Commit {res['stamp']['commit']}, tree clean {res['stamp']['tree_clean']}, UTC {res['stamp']['utc']}")
    return "\n".join(lines) + "\n"


def read(out: Path) -> dict:
    gv.refuse_existing(out, ("t61.json", "t61.txt"))
    started = gv.tree_state()
    tasks = out / "tasks"
    paths = {**{c: tasks / f for c, f in CUT_FILES.items()}, **{lib: tasks / f"lib{lib}.json" for lib in LIBRARIES}}
    loaded = {k: json.loads(p.read_text()) if p.exists() else None for k, p in paths.items()}
    reading = t61_reading({c: loaded[c] for c in CUTS}, {lib: loaded[lib] for lib in LIBRARIES},
                          {c: committed_t11(c) for c in CUTS}, {lib: committed_t60(lib) for lib in LIBRARIES})
    res = {"label": LABEL, "row": ROW, "reading": reading, "tries": TRIES_COUNT}
    res["stamp"] = gv.stamp(started, [p for p in paths.values() if p.exists()])
    out.mkdir(parents=True, exist_ok=True)
    (out / "t61.json").write_text(json.dumps(res, indent=1, default=float) + "\n")
    (out / "t61.txt").write_text(t61_text(res))
    print("\n".join(outcome_lines(reading)) + f"\n(numbers in {out / 't61.txt'})")
    return res


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("run-cut", help="one E10.5 cut of try T61")
    p.add_argument("--cut", required=True, choices=list(CUTS))
    p.add_argument("--out", required=True)
    p.add_argument("--seeds", default=",".join(map(str, SEEDS)))
    p.add_argument("--data", default=str(harness.DATA_DIR))
    p = sub.add_parser("run-lib", help="one library task of try T61")
    p.add_argument("--library", required=True, choices=list(LIBRARIES))
    p.add_argument("--out", required=True)
    p.add_argument("--seeds", default=",".join(map(str, SEEDS)))
    p.add_argument("--data", default=str(harness.DATA_DIR))
    p = sub.add_parser("read", help="try T61's reading")
    p.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    if a.cmd == "read":
        read(RESULTS / a.out)
    else:
        {"run-cut": cmd_run_cut, "run-lib": cmd_run_lib}[a.cmd](a)


if __name__ == "__main__":
    main()
