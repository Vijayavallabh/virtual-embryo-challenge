#!/usr/bin/env python
"""Try T62 (experiments/harness/TRIES.md): the recomputed-groups file for T1's last upload, read on the practice tasks.

The candidate for T1's last upload is `build_test.py t1 --g-mom --base lib3 --feng-groups recomputed`: the first
upload's design with the outside atlas's direction formed on the cell groups row N54 recomputed from the atlas's
permitted cells (`feng_dirs.feng_directions(genes, groups="recomputed")`). No practice run has scored it. This script
scores its practice analogue against the ruled file's (row R8's) practice analogue, on try T61's five tasks, draws,
seeds, protocol and scoring. It imports bets_read.py (T61) and gmom_combo.py (T26) and copies neither's arithmetic.

Tasks: T61's five. The two E9.5 -> E10.5 cuts of configs/pseudo_tasks.yaml (chamber, primary; whole), start E8.5, as
T61 and T26 ran them; row T60's three library tasks as T61 ran them. Every task scores one truth, E10.5.

Arms per seed, from one cell draw:
  none_gmom   the control, not counted: built exactly as T61 builds its none_gmom (`bets_read.file_hook` on
              `bets_read.gmom_fits`'s "file" factors; on the libraries on T60's oracle base)
  m1rc_gmom   the candidate, counted: try T26's m1_r_gmom construction with the recomputed direction. The same draw;
              the direction of feng_directions(genes, groups="recomputed") at the practice fraction
              (`feng_dirs.bracket_fraction` for the task's one-day step, 1/5 of the bracket) on the cells whose label's
              group (the readout config's regexes, `feng_dirs.label_groups`) has a recomputed direction after the
              yaml's fallback, holding R2's genes (the pre-registered list joined with R2's union, `gmom_combo.
              held_feng`, as the builder's `_t1_held`); on the cuts through `t1_scope.build_submission` (stream
              [seed, 2]) as T26 did, on the libraries through `gmom_combo.moved_by_feng` (the same stream) as row T60's
              entry_step did. Then g_mom as none_gmom's (the same factors, stream [seed, 5], no rescale, the
              pre-registered held genes, cardiomyocytes and endocardium only), with the builder's mask under
              recomputed (`build_test._labels_without_direction` on the groups with a direction, the with_direction
              argument of `_t1_g_mom`).
  m1rc_gmom_epi_held
              the second counted arm: m1rc_gmom with the Epicardium group's step taken out. Its cells are treated as
              the builder treats a group without a direction under recomputed: not moved by the direction (they keep
              their counts as drawn, R1's hold) and outside g_mom's covered labels and mask (with_direction without
              Epicardium). Every other cell is m1rc_gmom's, byte for byte: the same draw, the same direction step on
              stream [seed, 2] and the same g_mom draw on stream [seed, 5] (`hold_group`, tested). The direction is
              not removed before the step itself, because `vec.core.apply_factors` draws one stream over every moved
              cell: removing a group there would also redraw the other groups' cells (tested), and the paired
              difference would then carry sampling noise from cells the question is not about.
  A group without a recomputed direction, and a label with no group: its cells are not moved by the direction and
  get no momentum (R1: they keep their counts as drawn, the E9.5 cells' own counts on these tasks), and g_mom's mask
  holds them back too, even inside cardiomyocytes and endocardium. On the configs as committed, with CM, Endocardium
  and Endothelium directed, that mask holds back the same labels as none_gmom's (Early VIC alone; tested).
On each library task the reference copy_last (row T60's, not counted) is scored first, as T61 scored it.

  reliability  Target-free: no T1 cell is loaded or scored, no truth is read, no authors' group is used. From the
           recomputed statistics file alone (through feng_dirs' checks) and the T1 gene names (T1:val's gene order),
           for each recomputed group with a direction: the cells in each half at each end of the bracket; the Spearman
           correlation of the two halves' log ratios (`feng_dirs.direction_halves`, the readout's settings) over all
           mapped genes and within each expression decile (`feng_dirs.eb_shrink`'s basis and binning); eb_shrink's
           ten factors and their mean (the share of the direction's variance the halves agree on). Writes
           tasks/reliability.json and tasks/reliability.txt; prints only their path. It decides nothing.
  run-cut  One E10.5 cut. Refused before anything loads unless row T62 is in TRIES.md, the tree is clean, the task file
           does not exist, and the relabelling's checks pass (`relabel_checks`). Prints no score.
  run-lib  One library task; the same refusals. Prints no score.
  read     T62's statement, applied before anything is printed; writes t62.json and t62.txt; prints the outcome only.

The statement (fixed in row T62 before the run, T61's bar). Void (nothing follows, no rerun under the row): a task, seed
or arm missing; the five task files not at one commit or a tree not clean; on any task and seed none_gmom off T61's
committed none_gmom by more than 1e-9; the recomputed statistics file's sha256 not STATS_SHA256, or the committed record
not saying route recomputed and authors_labels_read false; the settings file's sha256 not YAML_SHA256; this script's
source calling feng_directions with any groups other than "recomputed" (`source_guard`).
"The recomputed step is AHEAD": on the chamber cut the ten-seed paired mean of m1rc_gmom minus none_gmom above +0.3 with
at least 8 of 10 seeds above 0, and on the whole cut the paired mean at least -0.3. Reading S (the first upload
stands): AHEAD or not. Reading D (it is discounted): AHEAD, at least two library tasks valid, and no valid library task
with a paired mean below -0.3. A library task is valid unless a seed has inverted anchors (T61's library rule).
"The epicardium step PAYS": on library 3's task the ten-seed paired mean of m1rc_gmom minus m1rc_gmom_epi_held above
+0.3 with at least 8 of 10 seeds above 0; "COSTS": below -0.3 with at least 8 of 10 seeds below 0; otherwise "within
noise". The same paired difference is reported on every task.
Reported for every task: the paired means, their spread, seeds ahead, the per-group skill differences and points, and
the arms' mean scores; J6's range from zero to the chamber cut's gain.
What follows: nothing is built or changed. The advice fixed before the run: the recomputed file if AHEAD under both
readings, the ruled file otherwise; if the recomputed file is advised and the epicardium step COSTS, the question also
names a third file (the recomputed groups with the epicardium group held), which needs an off-by-default builder
option, a check build and the captain's ruling.

NOT an adoption result: in-sample labelled diagnostic on one truth; the atlas's E9.5 cells lie at the tasks' own
reference stage; the labeller's threshold was moved on the released E9.5 (row N55); the practice fraction is half the
test's.

Usage: on the DGX campaign checkout ~/v/vec-campaign-d (it alone holds data/derived/feng_relabel_c/), repo root, at the
row's commit, committing nothing until the five task files exist:
  export NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
  uv run --frozen --no-sync python -I -B experiments/t1_m1/recomputed_read.py reliability --out t62_recomputed
  uv run --frozen --no-sync python -I -B experiments/t1_m1/recomputed_read.py run-cut --out t62_recomputed --cut chamber
  (and --cut whole, run-lib --library 1, 2, 3), then, at the commit that holds the five task files,
  uv run --frozen --no-sync python -I -B experiments/t1_m1/recomputed_read.py read --out t62_recomputed
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
import time
import warnings
from collections.abc import Mapping, Sequence
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

from vec import harness, t1_scope
from vec import t1_state_transport as ts
from vec.adoption import SEEDS
from vec.boards import GROUPS, get_board, load_boards
from vec.evaluate import Stage
from vec.io import dense, read_stage
from vec.prep import load_stage

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RESULTS = HERE / "results"
TRIES = ROOT / "experiments" / "harness" / "TRIES.md"
sys.path.insert(0, str(HERE))
import bets_read as br  # noqa: E402  (T61: its draws, hook, fits, scoring helpers; also puts the builder on the path)
import gmom_combo as gc  # noqa: E402  (T26: moved_by_feng, held_feng)

bt, gv, pt, rr = br.bt, br.gv, br.pt, br.rr
feng_dirs = bt.feng_dirs
hygiene = gc.hygiene

ROW = "T62"
CUTS, CUT_FILES, LIBRARIES = br.CUTS, br.CUT_FILES, br.LIBRARIES
T61_TASKS = RESULTS / "t61_bets" / "tasks"
CONTROL, CANDIDATE, EPI_HELD = "none_gmom", "m1rc_gmom", "m1rc_gmom_epi_held"
ARMS = (CONTROL, CANDIDATE, EPI_HELD)
LIB_REFERENCE = "copy_last"  # scored first, as T61 scored it; not counted
LIB_ARMS = (LIB_REFERENCE, *ARMS)
COUNTED = (CANDIDATE, EPI_HELD)
TRIES_COUNT = len(COUNTED) * (len(CUTS) + len(LIBRARIES))  # 10: T1 358 to 368
EPI_GROUP = "Epicardium"  # the group whose step the second counted arm takes out
EPI_TASK = "3"  # library 3, where the epicardium group has its test-like weight
EPI_BAR = 0.3
MIN_EPI = 8
DECILES = 10  # eb_shrink's bins (harness.LP_BINS)
RELIABILITY_FILES = ("reliability.json", "reliability.txt")
GROUPS_OPTION = "recomputed"
REPRO_TOL = 1e-9  # none_gmom against T61's committed none_gmom, per task and seed
AHEAD_BAR = 0.3
MIN_AHEAD = 8
WHOLE_FLOOR = -0.3  # the whole cut's paired mean must be at least this
LIB_FLOOR = -0.3  # under D no valid library task's paired mean below this
MIN_LIBRARIES = 2
PRACTICE_FRACTION = 0.2  # one day of Feng's E9.5 -> E14.5 bracket (the test's two days: 0.4)
STATS_SHA256 = "36d1f7edcbbf3d7a80c68b51d5624a94bb753668079fe3d69aed720c71baa7e9"
YAML_SHA256 = "883e74c2c1ef5f01d5244555f8b3dcae24255568c731a6d53ee2788b7b677b49"
STATS_NAME = "cd1_stats_relabel.npz"
RECORD = ROOT / "experiments" / "external" / "results" / "feng_relabel" / "MANIFEST.json"
GROUP_WEIGHTS = {g: w for g, w, _ in GROUPS["T1"]}
LABEL = ("try T62, a labelled in-sample diagnostic, never an adoption result: T1's recomputed-groups file "
         "(`--feng-groups recomputed`) against the ruled file (row R8) by their practice analogues m1rc_gmom and "
         "none_gmom, and the epicardium group's step by m1rc_gmom against m1rc_gmom_epi_held, on T61's five tasks; "
         "every task scores one truth (E10.5) over one day, the test spans two")
RELIABILITY_LABEL = ("try T62's reliability listing: target-free, from the recomputed statistics file's own two halves "
                     "of each group's cells; no T1 cell, no truth, no authors' group; it decides nothing")


# ---------------------------------------------------------------------------------------------------------------------
# Guards on the relabelling and on this source


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def relabel_checks(record: Path = RECORD, yaml_path: Path | None = None, stats_path: Path | None = None,
                   need_stats: bool = True) -> dict:
    """The relabelling's guards, read from committed files (the record, the frozen yaml) and, when need_stats, the
    recomputed statistics file's sha256 (its bytes are hashed, not loaded). Returns the values and the failed checks."""
    yaml_path = feng_dirs.RELABEL_CONFIG if yaml_path is None else Path(yaml_path)
    reasons = []
    rec = json.loads(Path(record).read_text()) if Path(record).exists() else {}
    if not rec:
        reasons.append(f"no committed relabelling record at {record}")
    if rec.get("route") != GROUPS_OPTION:
        reasons.append(f"the committed record's route is {rec.get('route')!r}, not 'recomputed'")
    if rec.get("authors_labels_read") is not False:
        reasons.append(f"the committed record's authors_labels_read is {rec.get('authors_labels_read')!r}, not false")
    rec_stats = (rec.get("outputs") or {}).get(STATS_NAME, {}).get("sha256")
    if rec_stats != STATS_SHA256:
        reasons.append(f"the record's statistics sha256 is {rec_stats}, not the row's {STATS_SHA256}")
    if (rec.get("config") or {}).get("sha256") != YAML_SHA256:
        reasons.append("the record was written from another settings file")
    yaml_sha = sha256_file(yaml_path) if yaml_path.exists() else None
    if yaml_sha != YAML_SHA256:
        reasons.append(f"the settings file's sha256 is {yaml_sha}, not the row's {YAML_SHA256}")
    stats_sha = None
    if need_stats:
        path = feng_dirs.recomputed_stats_path() if stats_path is None else Path(stats_path)
        stats_sha = sha256_file(path) if path.exists() else None
        if stats_sha != STATS_SHA256:
            reasons.append(f"the recomputed statistics file's sha256 is {stats_sha}, not the row's {STATS_SHA256}")
    return {"route": rec.get("route"), "authors_labels_read": rec.get("authors_labels_read"),
            "record_stats_sha256": rec_stats, "yaml_sha256": yaml_sha, "stats_sha256": stats_sha,
            "passes": not reasons, "reasons": reasons}


def _calls(text: str, name: str = "feng_directions") -> list[str]:
    """Every `name(...)` in `text` with its balanced arguments (a call cut off by the end of the text is returned as
    far as it goes, so it fails the check)."""
    out = []
    for m in re.finditer(re.escape(name) + r"\(", text):
        depth, i = 1, m.end()
        while i < len(text) and depth:
            depth += {"(": 1, ")": -1}.get(text[i], 0)
            i += 1
        out.append(text[m.start():i])
    return out


def source_guard(path: Path | None = None) -> list[str]:
    """This script's source never asks for any groups but the recomputed ones: every call of feng_directions names
    groups="recomputed", and nothing reaches the authors' groups, their statistics or row N54's listing. Returns the
    failures (the tokens are spelled in pieces here so that this function does not find itself)."""
    text = Path(__file__ if path is None else path).read_text(encoding="utf-8")
    bad = [f"a feng_directions call without groups=\"recomputed\": {c}" for c in _calls(text)
           if 'groups="recomputed"' not in c]
    q = '"'
    for token in (q + "auth" + "ors" + q, "'" + "auth" + "ors" + "'", "FENG_GROUPS" + "[0]", "clean_" + "stats_path",
                  "directions_" + "from_stats", "n54" + "_relabel_listing", "n54" + ".json", "n54" + "_authors",
                  "R1_FENG" + "_GROUPS", "r1_unbracketed" + "_labels"):
        if token in text:
            bad.append(f"the source names {token}")
    return bad


def recomputed_inputs(genes: Sequence[str]) -> tuple[dict, dict]:
    """The recomputed directions on `genes` and their record, through feng_dirs' own checks; refused unless the route
    is recomputed and Blood has no direction (the builder's guards)."""
    directions, ext = feng_dirs.feng_directions(list(genes), groups="recomputed")
    if ext.get("route") != GROUPS_OPTION:
        raise ValueError(f"the relabelling routes to {ext.get('route')!r}, not 'recomputed'")
    if "Blood" in directions:
        raise ValueError("Blood has a recomputed direction (R2's guard)")
    return directions, ext


# ---------------------------------------------------------------------------------------------------------------------
# The candidate's g_mom step and arms


def rc_hook(factors: ts.StateFactors, seed: int, held_flat: Sequence[str], with_direction: Sequence[str],
            covered: set[str]):
    """g_mom as none_gmom's (`bets_read.file_hook`: the same factors, stream [seed, 5], no rescale, `held_flat`,
    G_MOM_EDIT_CATEGORIES only) with the builder's mask under --feng-groups recomputed: the labels whose group has no
    recomputed direction (`build_test._labels_without_direction`) keep their counts; the builder's covered check."""
    missing = sorted(set(bt.G_MOM_EDIT_CATEGORIES) - set(factors.categories))
    if missing:
        raise ValueError(f"g_mom has no factors for a category R3b edits: {missing}")

    def hook(pop):
        labels = pop.obs[t1_scope.LABEL_KEY].astype(str).to_numpy()
        cats = pop.obs[t1_scope.CATEGORY_KEY].astype(str).to_numpy()
        held_labels = bt._labels_without_direction(np.unique(labels), with_direction)
        no_edit = np.isin(labels, sorted(held_labels)) | ~np.isin(cats, bt.G_MOM_EDIT_CATEGORIES)
        edited = set(labels[~no_edit & np.isin(cats, list(factors.categories))])
        if not edited <= set(covered):
            raise ValueError(f"g_mom would edit labels the recomputed direction does not cover: "
                             f"{sorted(edited - set(covered))}")
        return ts.apply_state_factors(pop, factors, list(held_flat), [seed, br.GMOM_STREAM], rescale=False,
                                      no_edit=no_edit)[0]

    return hook


def direction_kit(labels, directions: Mapping) -> tuple[dict, set[str], list[str]]:
    """The label -> group map (the readout config's regexes), the labels covered by a recomputed direction (the
    builder's `_g_mom_covered` under recomputed) and the groups with a direction."""
    groups = feng_dirs.label_groups(sorted({str(x) for x in labels}))
    return groups, bt._g_mom_covered(groups, directions, GROUPS_OPTION), sorted(directions)


def without_epicardium(directions: Mapping) -> dict:
    """The recomputed directions with the Epicardium group's taken out (the second counted arm's groups)."""
    return {g: v for g, v in directions.items() if g != EPI_GROUP}


def hold_group(moved, drawn, groups: Mapping[str, str], group: str = EPI_GROUP):
    """`moved` (after the direction step) with the cells whose label maps to `group` given back their counts as drawn
    (`drawn`, the same draw before the step: R1's hold, as the builder leaves a group without a direction), through
    the builder's own row stacking (`t1_scope._stack_rows`). Every other cell is `moved`'s, byte for byte."""
    if list(moved.obs_names) != list(drawn.obs_names):
        raise ValueError("the moved and drawn populations are not the same draw")
    labels = moved.obs[t1_scope.LABEL_KEY].astype(str).to_numpy()
    held = np.array([groups.get(lab) == group for lab in labels], dtype=bool)
    if not held.any():
        return moved
    if held.all():
        return drawn.copy()
    keep, back = np.flatnonzero(~held), np.flatnonzero(held)
    return t1_scope._stack_rows([(keep, moved[keep].copy()), (back, drawn[back].copy())], moved.n_obs)


def cut_pair(ctx: harness.BuildContext, factors: ts.StateFactors, directions: Mapping, fraction: float,
             held_m1: Sequence[str]) -> dict[str, Stage]:
    """none_gmom (T61's), m1rc_gmom (T26's m1_r_gmom with the recomputed direction, the builder's mask) and
    m1rc_gmom_epi_held (m1rc_gmom with the Epicardium group's cells held as drawn) from one draw on a harness cut,
    each finished as the harness finishes a candidate."""
    role = ctx.task.reference_role
    ref, src = harness._shallow(ctx.inputs[role]), ctx.task.role_stems[role]
    spec = t1_scope.load_spec()
    t1_scope.classify(ref, source=src, spec=spec)
    shares = harness.t1_chamber_shares(ref, src)
    n = harness.count_for(ctx.board, ref.n_obs)
    held = list(ctx.derived["t1_held_flat"])
    labels = ref.obs[t1_scope.LABEL_KEY].astype(str).unique()
    groups, covered, with_direction = direction_kit(labels, directions)
    _, covered_h, with_direction_h = direction_kit(labels, without_epicardium(directions))
    control_hook, drawn = br.file_hook(factors, ctx.seed, held), []

    def control(pop):
        drawn.append(pop)  # the draw itself: no direction and no momentum come before post_cells on this arm
        return control_hook(pop)

    def build(**kw) -> Stage:
        sub, _ = t1_scope.build_submission(ref, ctx.board, shares, n, ctx.seed, ref_source=src, **kw)
        return harness.finish(ctx, ref, dense(sub.X), None)

    moved_kw = {"held_flat": list(held_m1), "directions": directions, "direction_fraction": fraction,
                "direction_groups": groups}
    epi_hook = rc_hook(factors, ctx.seed, held, with_direction_h, covered_h)
    arms = {CONTROL: build(held_flat=held, post_cells=control)}  # first: it hands the draw to the third arm
    arms[CANDIDATE] = build(**moved_kw, post_cells=rc_hook(factors, ctx.seed, held, with_direction, covered))
    arms[EPI_HELD] = build(**moved_kw, post_cells=lambda moved: epi_hook(hold_group(moved, drawn[0], groups)))
    return arms


def library_pair(e95, half_shares: Mapping[str, float], factors: ts.StateFactors, held: Sequence[str],
                 held_m1: Sequence[str], seed: int, within: Mapping, fb: Mapping[str, Mapping], directions: Mapping,
                 fraction: float, n_cells: int = rr.N_CELLS) -> dict:
    """none_gmom (T61's, on T60's oracle base), m1rc_gmom on the same draw (row T60's entry_step route:
    `gmom_combo.moved_by_feng`, stream [seed, 2], then the candidate's g_mom) and m1rc_gmom_epi_held (the same, with
    the Epicardium group's cells held as drawn), as AnnData."""
    base, _ = t1_scope.reweight(e95, half_shares, n_cells, seed, within=within, fallback=fb["E9.5"])
    labels = e95.obs[t1_scope.LABEL_KEY].astype(str).unique()
    groups, covered, with_direction = direction_kit(labels, directions)
    _, covered_h, with_direction_h = direction_kit(labels, without_epicardium(directions))
    moved = gc.moved_by_feng(base, directions, fraction, groups, list(held_m1), [seed, 2])
    return {CONTROL: br.file_hook(factors, seed, held)(base),
            CANDIDATE: rc_hook(factors, seed, held, with_direction, covered)(moved),
            EPI_HELD: rc_hook(factors, seed, held, with_direction_h, covered_h)(hold_group(moved, base, groups))}


# ---------------------------------------------------------------------------------------------------------------------
# The reliability listing (target-free)


def _spearman(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) < 3:
        return float("nan")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # a constant half gives NaN, reported as such
        return float(spearmanr(x, y).statistic)


def reliability_listing(st: Mapping, genes: Sequence[str], cfg: Mapping, directions: Mapping,
                        fallbacks: Mapping | None = None, bracket: tuple[str, str] = feng_dirs.BRACKET,
                        bins: int = DECILES) -> dict:
    """Per group with a recomputed direction, from the statistics' own two halves of the group's cells (`cfg`: the
    readout's settings with the recomputed group map, `feng_dirs.groups_config`): the cells in each half at each end of
    the bracket; the Spearman correlation of the halves' log ratios over every mapped gene and within each expression
    decile (eb_shrink's basis and binning); eb_shrink's factors and their mean. A group whose direction came from the
    yaml's fallback has no halves of its own and is listed as such. Reads no T1 cell and no truth."""
    halves = feng_dirs.direction_halves(st, list(genes), cfg, bracket, groups=set(directions))
    out = {}
    for g in sorted(directions):
        if g not in halves:
            out[g] = {"own_direction": False, "fallback": (fallbacks or {}).get(g)}
            continue
        labels = cfg["external_labels"][g][0]
        cells = {stage: {h: int(feng_dirs.dl.feng_sel(st, stage, labels, halves=(h,))["n"]) for h in ("A", "B")}
                 for stage in bracket}
        h = halves[g]
        d, d_a, d_b, basis = (np.asarray(x, dtype=np.float64)
                              for x in (directions[g], h["d_a"], h["d_b"], h["mean_cp"]))
        idx = np.flatnonzero(np.isfinite(d) & np.isfinite(d_a) & np.isfinite(d_b) & np.isfinite(basis))
        runs = np.array_split(idx[np.argsort(basis[idx], kind="stable")], bins)
        _, factors = feng_dirs.eb_shrink(d, d_a, d_b, basis, bins)
        deciles = [{"decile": i + 1, "genes": int(r.size),
                    "mean_cp": [float(basis[r].min()), float(basis[r].max())] if r.size else None,
                    "spearman": _spearman(d_a[r], d_b[r]), "eb_factor": float(f)}
                   for i, (r, f) in enumerate(zip(runs, factors, strict=True))]
        out[g] = {"own_direction": True, "labels": list(labels), "cells": cells, "mapped_genes": int(idx.size),
                  "spearman_all": _spearman(d_a[idx], d_b[idx]), "deciles": deciles,
                  "eb_factors": [float(f) for f in factors], "eb_factor_mean": float(np.mean(factors))}
    return out


def reliability_text(res: Mapping) -> str:
    lines = [res["label"], f"Bracket {' -> '.join(res['bracket'])}; genes: {res['genes_from']}'s order "
             f"({res['of_genes']}); groups with a direction: {', '.join(res['with_direction'])}"]
    for g, r in res["groups"].items():
        if not r["own_direction"]:
            lines.append(f"{g}: no direction of its own (the yaml's fallback: {r['fallback']}); no halves to read")
            continue
        ends = "; ".join(f"{stage}: half A {v['A']}, half B {v['B']}" for stage, v in r["cells"].items())
        lines.append(f"{g}: cells {ends}; mapped genes {r['mapped_genes']}; halves' Spearman over all mapped genes "
                     f"{r['spearman_all']:.3f}; eb_shrink factor mean {r['eb_factor_mean']:.3f}")
        for x in r["deciles"]:
            lines.append(f"  decile {x['decile']:2d} ({x['genes']} genes): Spearman {x['spearman']:.3f}, eb factor "
                         f"{x['eb_factor']:.3f}")
    lines.append(f"Commit {res['stamp']['commit']}, tree clean {res['stamp']['tree_clean']}, UTC {res['stamp']['utc']}")
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------------------------------------------------
# Runs


def refuse(out: Path, name: str | Sequence[str]) -> tuple[dict, dict]:
    """The refusals before anything loads: the row, an existing output file, an unclean tree, the source guard and the
    relabelling's guards. Returns the start state and the relabelling's checks."""
    if ROW not in br.tries_rows(TRIES):
        sys.exit(f"{ROW} is refused: it has no row in experiments/harness/TRIES.md (a row comes before its run)")
    gv.refuse_existing(out, (name,) if isinstance(name, str) else tuple(name))
    started = gv.tree_state()
    if started["clean"] is not True:
        sys.exit(f"{ROW} is refused: the tree is not clean by vec.provenance's rule")
    bad = source_guard()
    if bad:
        sys.exit(f"{ROW} is refused: {bad}")
    checks = relabel_checks()
    if not checks["passes"]:
        sys.exit(f"{ROW} is refused: {checks['reasons']}")
    return started, checks


def _relabel_fields(ext: Mapping, checks: Mapping) -> dict:
    return {"file": ext["file"], "sha256": ext["sha256"], "stats_sha256_hashed": checks["stats_sha256"],
            "yaml_sha256": checks["yaml_sha256"], "route": ext["route"], "record_route": checks["route"],
            "authors_labels_read": checks["authors_labels_read"], "with_direction": ext["with_direction"],
            "fallbacks": ext["fallbacks"], "relabel_gate": ext.get("relabel_gate"),
            "epicardium_directed": EPI_GROUP in ext["with_direction"],
            "with_direction_epi_held": [g for g in ext["with_direction"] if g != EPI_GROUP]}


def cmd_reliability(a) -> None:
    """Row T62's reliability listing: the recomputed statistics file and the T1 gene names only; no T1 cell, no
    truth."""
    out = RESULTS / a.out / "tasks"
    started, checks = refuse(out, RELIABILITY_FILES)
    genes = list(get_board(rr.BOARD).genes)  # the gene names alone (the libraries' gene order)
    directions, ext = recomputed_inputs(genes)
    path = feng_dirs.recomputed_stats_path()
    feng_dirs.check_recomputed(path)
    st = feng_dirs.dl.load_stats(path)
    feng_dirs.check_recomputed_keys(st)
    cfg = feng_dirs.groups_config(feng_dirs.load_config(), feng_dirs.load_relabel_config())
    res = {"label": RELIABILITY_LABEL, "row": ROW, "bracket": list(feng_dirs.BRACKET), "genes_from": rr.BOARD,
           "of_genes": len(genes), "with_direction": sorted(directions), "relabel": _relabel_fields(ext, checks),
           "settings": {"pseudocount_cp10k": ext["pseudocount_cp10k"], "min_cells_ext": ext["min_cells_ext"],
                        "deciles": DECILES},
           "groups": reliability_listing(st, genes, cfg, directions, ext["fallbacks"]),
           "stamp": gv.stamp(started, [path, feng_dirs.RELABEL_CONFIG, feng_dirs.CONFIG, RECORD])}
    out.mkdir(parents=True, exist_ok=True)
    (out / RELIABILITY_FILES[0]).write_text(json.dumps(res, indent=1, default=float) + "\n")
    (out / RELIABILITY_FILES[1]).write_text(reliability_text(res))
    print(f"reliability listing written: {out / RELIABILITY_FILES[1]}")  # the path only


def _held_m1(preregistered: Sequence[str]) -> tuple[list[str], dict]:
    union = hygiene.load_union()
    return gc.held_feng(union, preregistered), union


def cmd_run_cut(a) -> None:
    out = RESULTS / a.out / "tasks"
    started, checks = refuse(out, CUT_FILES[a.cut])
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
    directions, ext = recomputed_inputs(list(inputs[ref_role].var_names))
    fraction = feng_dirs.bracket_fraction(harness.STAGE_DAYS[t.reference], harness.STAGE_DAYS[t.target])
    if abs(fraction - PRACTICE_FRACTION) > 1e-12:
        sys.exit(f"{t.key}: the practice fraction is {fraction}, not {PRACTICE_FRACTION}")
    held_m1, union = _held_m1(preregistered)
    factors = br.gmom_fits(inputs[ref_role], view.role_stems[ref_role], inputs["start"],
                           view.role_stems["start"])["file"]
    rows = []
    for s in seeds:
        t0 = time.time()
        f = gv.floor_anchors(t, board, target, ref_stage, s)
        ctx = harness.BuildContext(task=view, board=cboard, inputs=inputs, seed=harness.BUILD_SEED_OFFSET + s,
                                   derived={"t1_held_flat": preregistered})
        arms = gv.score_arms(board, target, ref_stage, s, f.anchors, cut_pair(ctx, factors, directions, fraction,
                                                                              held_m1))
        rows.append({"seed": s, "floor": float(f.score), "arms": arms, "seconds": round(time.time() - t0, 1)})
        print(f"{t.key} seed {s}: {time.time() - t0:.0f} s", flush=True)  # no score
    res = {"label": LABEL, "row": ROW, "task": t.key, "seeds": list(seeds), "rows": rows, "fraction": fraction,
           "relabel": _relabel_fields(ext, checks), "constants": {"edit_categories": list(bt.G_MOM_EDIT_CATEGORIES),
                                                                 "held_m1_genes": len(held_m1)},
           "stamp": gv.stamp(started, [data / t.target, data / t.reference, *[data / f for f in t.inputs.values()],
                                       t1_scope.HELD_FLAT_GENES, t1_scope.CHAMBER_CLASSES, ROOT / ext["file"],
                                       feng_dirs.RELABEL_CONFIG, RECORD, ROOT / union["file"]])}
    out.mkdir(parents=True, exist_ok=True)
    (out / CUT_FILES[a.cut]).write_text(json.dumps(res, indent=1, default=float) + "\n")


def cmd_run_lib(a) -> None:
    out = RESULTS / a.out / "tasks"
    name = f"lib{a.library}.json"
    started, checks = refuse(out, name)
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
    factors = br.gmom_fits(e95, "E9.5_RNA", e85, "E8.5_RNA")["file"]
    directions, ext = recomputed_inputs(list(e95.var_names))
    fraction = feng_dirs.bracket_fraction(9.5, 10.5)
    if abs(fraction - PRACTICE_FRACTION) > 1e-12:
        sys.exit(f"library {a.library}: the practice fraction is {fraction}, not {PRACTICE_FRACTION}")
    held_m1, union = _held_m1(held)
    rows = []
    for s in seeds:
        t0 = time.time()
        pick = np.sort(np.random.default_rng([s, 99]).choice(e95.n_obs, rr.N_CELLS, replace=False))
        arms = library_pair(e95, shares, factors, held, held_m1, s, within, fb, directions, fraction)
        stages = {LIB_REFERENCE: pt.to_stage(e95[pick]), **{k: pt.to_stage(v) for k, v in arms.items()}}
        rows.append({"seed": s, **rr.score_arms(board, stages, truth_stage, ref_stage, s),
                     "seconds": round(time.time() - t0, 1)})
        print(f"library {a.library} seed {s}: {time.time() - t0:.0f} s", flush=True)  # no score
    res = {"label": LABEL, "row": ROW, "library": a.library, "seeds": list(seeds), "rows": rows, "fraction": fraction,
           "relabel": _relabel_fields(ext, checks),
           "constants": {"split_seed": rr.SPLIT_SEED, "n_cells": rr.N_CELLS, "board": rr.BOARD,
                         "edit_categories": list(bt.G_MOM_EDIT_CATEGORIES), "held_m1_genes": len(held_m1)},
           "stamp": gv.stamp(started, [*(data / f for f in rr.FILES.values()), t1_scope.HELD_FLAT_GENES,
                                       t1_scope.CHAMBER_CLASSES, ROOT / ext["file"], feng_dirs.RELABEL_CONFIG, RECORD,
                                       ROOT / union["file"]])}
    out.mkdir(parents=True, exist_ok=True)
    (out / name).write_text(json.dumps(res, indent=1, default=float) + "\n")


# ---------------------------------------------------------------------------------------------------------------------
# The reading


def paired(rows: Sequence[Mapping]) -> dict:
    """m1rc_gmom - none_gmom per seed (`bets_read.paired`: mean, sd, seeds above 0, range, points per group) with the
    raw per-group skill differences and the two arms' mean scores."""
    p = br.paired(rows, CANDIDATE, CONTROL)
    p["skill_differences"] = {g: float(np.mean([r["arms"][CANDIDATE]["groups"][g] - r["arms"][CONTROL]["groups"][g]
                                                for r in rows])) for g in GROUP_WEIGHTS}
    p["mean_scores"] = {a: float(np.mean([br.score_of(r, a) for r in rows])) for a in ARMS}
    return p


def epi_paired(rows: Sequence[Mapping]) -> dict:
    """m1rc_gmom - m1rc_gmom_epi_held per seed (`bets_read.paired`) with the seeds below 0."""
    p = br.paired(rows, CANDIDATE, EPI_HELD)
    p["behind"] = int(sum(br.score_of(r, CANDIDATE) - br.score_of(r, EPI_HELD) < 0 for r in rows))
    return p


def epi_statement(p: Mapping) -> str:
    """Row T62's epicardium statement on library 3's paired difference."""
    if p["mean"] > EPI_BAR and p["ahead"] >= MIN_EPI:
        return "PAYS"
    if p["mean"] < -EPI_BAR and p["behind"] >= MIN_EPI:
        return "COSTS"
    return "within noise"


def committed_t61(name: str, tasks: Path | None = None) -> dict[int, float]:
    rows = json.loads(((T61_TASKS if tasks is None else tasks) / name).read_text())["rows"]
    return {int(r["seed"]): float(r["arms"][CONTROL]["score"]) for r in rows}


def task_files() -> dict[str, str]:
    return {**CUT_FILES, **{lib: f"lib{lib}.json" for lib in LIBRARIES}}


def guards(loaded: Mapping[str, Mapping | None], committed: Mapping[str, Mapping[int, float]], checks: Mapping,
           source_bad: Sequence[str], seeds: Sequence[int] = SEEDS) -> list[str]:
    """Every void condition of row T62; an empty list means none fired."""
    reasons = []
    commits = set()
    for k, r in loaded.items():
        if r is None:
            reasons.append(f"{k}: no task file")
            continue
        rows = r.get("rows") or []
        if sorted(int(x["seed"]) for x in rows) != list(seeds):
            reasons.append(f"{k}: seeds {[x['seed'] for x in rows]}, not {list(seeds)}")
        want = LIB_ARMS if k in LIBRARIES else ARMS
        absent = sorted({a for x in rows for a in want if a not in x["arms"]})
        if absent:
            reasons.append(f"{k}: arms {absent} missing")
            continue
        stamp = r.get("stamp") or {}
        commits.add(stamp.get("commit"))
        if stamp.get("tree_clean") is not True:
            reasons.append(f"{k}: its tree was not clean at one commit")
        ref = committed.get(k, {})
        worst = max((abs(br.score_of(x, CONTROL) - ref[int(x["seed"])]) if int(x["seed"]) in ref else float("inf")
                     for x in rows), default=float("inf"))
        if not worst <= REPRO_TOL:
            reasons.append(f"{k}: none_gmom differs from T61's committed none_gmom by {worst:.3g} (tolerance "
                           f"{REPRO_TOL:g})")
        rel = r.get("relabel") or {}
        if rel.get("sha256") != STATS_SHA256 or rel.get("stats_sha256_hashed") != STATS_SHA256:
            reasons.append(f"{k}: ran on another recomputed statistics file")
        if rel.get("yaml_sha256") != YAML_SHA256:
            reasons.append(f"{k}: ran under another settings file")
        if rel.get("route") != GROUPS_OPTION or rel.get("record_route") != GROUPS_OPTION:
            reasons.append(f"{k}: the route was not recomputed")
        if rel.get("authors_labels_read") is not False:
            reasons.append(f"{k}: the record did not say authors_labels_read false")
    if len(commits) > 1:
        reasons.append(f"the task files ran at more than one commit: {sorted(map(str, commits))}")
    reasons += [f"at the read: {x}" for x in checks.get("reasons", [])]
    reasons += [f"source: {x}" for x in source_bad]
    return reasons


def library_valid(r: Mapping) -> bool:
    return not any(x.get("inverted_anchors") for x in r["rows"])


def t62_reading(loaded: Mapping[str, Mapping | None], committed: Mapping[str, Mapping[int, float]], checks: Mapping,
                source_bad: Sequence[str] = ()) -> dict:
    """Row T62's statement, applied to the task files."""
    reasons = guards(loaded, committed, checks, source_bad)
    if reasons:
        return {"outcome": "VOID", "reasons": reasons}
    tasks = {k: paired(r["rows"]) for k, r in loaded.items()}
    c, w = tasks["chamber"], tasks["whole"]
    ahead = bool(c["mean"] > AHEAD_BAR and c["ahead"] >= MIN_AHEAD and w["mean"] >= WHOLE_FLOOR)
    valid = [lib for lib in LIBRARIES if library_valid(loaded[lib])]
    low = [lib for lib in valid if tasks[lib]["mean"] < LIB_FLOOR]
    d = bool(ahead and len(valid) >= MIN_LIBRARIES and not low)
    word = {True: "AHEAD", False: "NOT AHEAD"}
    means = [t["mean"] for t in tasks.values()]
    epi = {k: epi_paired(r["rows"]) for k, r in loaded.items()}
    epi_word = epi_statement(epi[EPI_TASK])
    recomputed = bool(ahead and d)
    return {"outcome": f"S: {word[ahead]}; D: {word[d]}; epicardium step: {epi_word}", "reading_S": word[ahead],
            "reading_D": word[d], "ahead": ahead, "valid_libraries": valid, "libraries_below_floor": low,
            "tasks": tasks, "j6_range": sorted([0.0, c["mean"]]), "task_range_reported": [min(means), max(means)],
            "epicardium_step": epi_word, "epicardium_task_valid": EPI_TASK in valid, "epicardium": epi,
            "epicardium_not_directed": [k for k, r in loaded.items()
                                        if (r.get("relabel") or {}).get("epicardium_directed") is not True],
            "advice_fixed_before_the_run": "the recomputed file" if recomputed else "the ruled file (row R8)",
            "third_file_named": bool(recomputed and epi_word == "COSTS")}


def outcome_lines(r: Mapping) -> list[str]:
    if r["outcome"] == "VOID":
        return ["Outcome: VOID (" + "; ".join(r["reasons"]) + ")"]
    lines = [f"Outcome: reading S (T1's first upload stands): the recomputed step is {r['reading_S']}; reading D (it "
             f"is discounted): {r['reading_D']}",
             f"The epicardium step (m1rc_gmom - m1rc_gmom_epi_held on library {EPI_TASK}, a valid library: "
             f"{'yes' if r['epicardium_task_valid'] else 'no'}): {r['epicardium_step']}",
             f"The advice fixed before the run: {r['advice_fixed_before_the_run']}; the question names a third file "
             f"(the recomputed groups with the epicardium group held): {'yes' if r['third_file_named'] else 'no'}"]
    if r["epicardium_not_directed"]:
        lines.append(f"The Epicardium group had no recomputed direction on: {', '.join(r['epicardium_not_directed'])} "
                     f"(the two candidate arms are the same file there)")
    return lines


def t62_text(res: Mapping) -> str:
    r = res["reading"]
    lines = [res["label"], *outcome_lines(r)]
    if r["outcome"] != "VOID":
        for k, p in r["tasks"].items():
            g = ", ".join(f"{x} {v:+.4f} skill ({p['groups'][x]:+.2f} points)"
                          for x, v in p["skill_differences"].items())
            lines.append(f"{k}: m1rc_gmom - none_gmom {p['mean']:+.2f} points, sd {p['sd']:.2f}, ahead {p['ahead']} of "
                         f"{p['of']}, per seed {p['range'][0]:+.2f} to {p['range'][1]:+.2f}; mean scores none_gmom "
                         f"{p['mean_scores'][CONTROL]:.2f}, m1rc_gmom {p['mean_scores'][CANDIDATE]:.2f}, "
                         f"m1rc_gmom_epi_held {p['mean_scores'][EPI_HELD]:.2f}; by group: {g}")
        for k, e in r["epicardium"].items():
            g = ", ".join(f"{x} {v:+.2f} points" for x, v in e["groups"].items())
            lines.append(f"{k}: m1rc_gmom - m1rc_gmom_epi_held {e['mean']:+.2f} points, sd {e['sd']:.2f}, ahead "
                         f"{e['ahead']} and behind {e['behind']} of {e['of']}, per seed {e['range'][0]:+.2f} to "
                         f"{e['range'][1]:+.2f}; by group: {g}")
        lines.append(f"Valid libraries: {', '.join(r['valid_libraries']) or 'none'}; below {LIB_FLOOR}: "
                     f"{', '.join(r['libraries_below_floor']) or 'none'}")
        lines.append(f"J6 range (0 to the chamber cut's gain; a negative chamber gain is reported as [gain, 0] and "
                     f"cannot be AHEAD): {r['j6_range'][0]:+.2f} to {r['j6_range'][1]:+.2f}; the five tasks' paired "
                     f"means, reported only: {r['task_range_reported'][0]:+.2f} to {r['task_range_reported'][1]:+.2f}")
    lines.append(f"Commit {res['stamp']['commit']}, tree clean {res['stamp']['tree_clean']}, UTC {res['stamp']['utc']}")
    return "\n".join(lines) + "\n"


def read(out: Path) -> dict:
    gv.refuse_existing(out, ("t62.json", "t62.txt"))
    started = gv.tree_state()
    paths = {k: out / "tasks" / f for k, f in task_files().items()}
    loaded = {k: json.loads(p.read_text()) if p.exists() else None for k, p in paths.items()}
    committed = {k: committed_t61(f) for k, f in task_files().items()}
    checks = relabel_checks(need_stats=False)  # the stats file's hash is in each task file, checked by `guards`
    reading = t62_reading(loaded, committed, checks, source_guard())
    res = {"label": LABEL, "row": ROW, "reading": reading, "tries": TRIES_COUNT}
    res["stamp"] = gv.stamp(started, [p for p in paths.values() if p.exists()])
    out.mkdir(parents=True, exist_ok=True)
    (out / "t62.json").write_text(json.dumps(res, indent=1, default=float) + "\n")
    (out / "t62.txt").write_text(t62_text(res))
    print("\n".join(outcome_lines(reading)) + f"\n(numbers in {out / 't62.txt'})")
    return res


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("run-cut", help="one E10.5 cut of try T62")
    p.add_argument("--cut", required=True, choices=list(CUTS))
    p.add_argument("--out", required=True)
    p.add_argument("--seeds", default=",".join(map(str, SEEDS)))
    p.add_argument("--data", default=str(harness.DATA_DIR))
    p = sub.add_parser("run-lib", help="one library task of try T62")
    p.add_argument("--library", required=True, choices=list(LIBRARIES))
    p.add_argument("--out", required=True)
    p.add_argument("--seeds", default=",".join(map(str, SEEDS)))
    p.add_argument("--data", default=str(harness.DATA_DIR))
    p = sub.add_parser("reliability", help="try T62's target-free reliability listing of the recomputed groups")
    p.add_argument("--out", required=True)
    p = sub.add_parser("read", help="try T62's reading")
    p.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    if a.cmd == "read":
        read(RESULTS / a.out)
    else:
        {"run-cut": cmd_run_cut, "run-lib": cmd_run_lib, "reliability": cmd_reliability}[a.cmd](a)


if __name__ == "__main__":
    main()
