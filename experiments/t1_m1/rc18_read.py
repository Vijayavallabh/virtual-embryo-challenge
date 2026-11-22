#!/usr/bin/env python
"""Try T69 (experiments/harness/TRIES.md): T1 candidate rc18, the outside heart atlas's direction over E9.5 to E18 (the
chord), read on row T62's five practice tasks before any score of it exists.

R13's file moves each recomputed Feng group along F, the E9.5 -> E14.5 log ratio of mean CP10k, at 2/5 of that
bracket. rc18 replaces F by D18, the E9.5 -> E18 log ratio, at that bracket's calendar share, 2/8.5 at the test and
1/8.5 on these one-day tasks. In log ratios D18 = F + G exactly, G the E14.5 -> E18 late-maturation term. The E18 cells
are row N58's labelling (experiments/external/feng_relabel_e18.py: N54's control sets, N54's E9.5 and E14.5 rows copied
unchanged). The directions come from `feng_dirs.e18_late_term`, built from the same pieces as
`feng_dirs.recomputed_e18_directions`; no builder option exists before an AHEAD reading and the captain's yes (the T1
order rule, decision log 2026-11-21 (y)). This script imports recomputed_read.py (T62) and, through it,
bets_read.py (T61) and gmom_combo.py (T26), and copies none of their arithmetic; it imports double_step_read.py (T63)
for its source guard only.

Inputs, all target-free: the recomputed statistics of row N54's run (sha256 rc.STATS_SHA256), row N58's E18 statistics
(data/derived/feng_relabel_e18/, its manifest committed at HEAD as its RECORD, its sha256 named in row T69), the frozen
settings file (sha256 rc.YAML_SHA256) and the T1 gene names.

Tasks, draws and seeds: row T62's five (the two E9.5 -> E10.5 cuts, chamber primary, start E8.5, build seed 1,000,000 +
seed; T60's three library tasks, split 20261112, copy_last scored first), seeds 0 to 9, bake protocol, build main. One
truth, E10.5: a labelled in-sample diagnostic.

Arms per seed, from the same build seed and streams (every arm is T62's candidate build on its own directions and
fraction; only the directions and the fraction differ):
  m1rc_gmom   gate, not counted: T62's candidate (F at 1/5), rebuilt; within GATE_TOL of T62's committed value on every
              task and seed
  rc18        counted: D18 at 1/8.5
  rc18_sm     counted, the size control: per group F scaled to rc18's applied L2 log-factor size (the norm of D18 / 8.5
              over the genes the step moves: finite, not held), applied at 1/8.5
  rc18_perm   counted, the permuted control: F + G', G' each group's G permuted across genes within 10 equal-count bins
              of the group's mean CP10k averaged over E14.5 and E18 (row T15's procedure: vec.harness.
              t3_decile_permutation, seed [20261110, the group's index in vec.harness.LP_GROUPS]), at 1/8.5

  run-cut  One E10.5 cut. Refused before anything loads unless row T69 is in TRIES.md, the output file does not exist,
           the tree is clean, the source guards pass, row N54's checks and row N58's checks pass (its E9.5 and E14.5
           rows identical to N54's byte for byte) and X002 names the E18 use. Prints no score.
  run-lib  One library task; the same refusals. Prints no score.
  read     T69's statement, applied before anything is printed; writes t69.json and t69.txt; prints the outcome words.

The statement, fixed in row T69 before the run (the challenger's text). VOID: a gate fails (m1rc_gmom off T62's
committed m1rc_gmom by more than 0.02 on a task and seed; the E9.5 and E14.5 statistics not N54's byte for byte), or a
condition of `guards`. STOP (the line ends): library 3's paired mean of rc18 minus m1rc_gmom below -0.3. AHEAD only if
rc18 minus m1rc_gmom has a chamber mean above +0.3 with at least 8 of 10 seeds ahead, a whole-cut mean above 0, at
least two valid libraries and no valid library mean below -0.3, and rc18 beats both rc18_sm and rc18_perm on the
chamber cut in at least 8 of 10 seeds (the same seed beating both). Otherwise NOT AHEAD. J6: horizon ratio 2, the gain
shown as 0 to the practice value.

Usage: on a DGX checkout holding data/derived/feng_relabel_c/ and data/derived/feng_relabel_e18/, through
experiments/t1_m1/rc18.sh (one step per call).
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
from vec.adoption import SEEDS
from vec.boards import get_board, load_boards
from vec.evaluate import Stage
from vec.io import dense, read_stage
from vec.prep import load_stage

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RESULTS = HERE / "results"
TRIES = ROOT / "experiments" / "harness" / "TRIES.md"
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "experiments" / "external"))
import double_step_read as ds  # noqa: E402  (T63: its source guard only)
import feng_relabel_e18 as fre  # noqa: E402  (row N58's names and record)
import recomputed_read as rc  # noqa: E402  (T62: its arms' pieces, guards and checks)

br, gc = rc.br, rc.gc
bt, gv, pt, rr = rc.bt, rc.gv, rc.pt, rc.rr
feng_dirs = rc.feng_dirs

ROW = "T69"
CUTS, CUT_FILES, LIBRARIES = rc.CUTS, rc.CUT_FILES, rc.LIBRARIES
T62_TASKS = RESULTS / "t62_recomputed" / "tasks"
CANDIDATE = rc.CANDIDATE  # m1rc_gmom: the gate
RC18, SM, PERM = "rc18", "rc18_sm", "rc18_perm"
ARMS = (CANDIDATE, RC18, SM, PERM)
COUNTED = (RC18, SM, PERM)
LIB_REFERENCE = rc.LIB_REFERENCE  # copy_last, scored first
LIB_ARMS = (LIB_REFERENCE, *ARMS)
TRIES_COUNT = len(COUNTED) * (len(CUTS) + len(LIBRARIES))  # 15
GATE_TOL = 0.02
PRACTICE_FRACTION = rc.PRACTICE_FRACTION  # 0.2: F's one-day share of E9.5 -> E14.5
RC18_FRACTION = feng_dirs.bracket_fraction(9.5, 10.5, feng_dirs.E18_BRACKET_DAYS)  # 1/8.5
FRACTION_TOL = 1e-12
FRACTIONS = {CANDIDATE: PRACTICE_FRACTION, RC18: RC18_FRACTION, SM: RC18_FRACTION, PERM: RC18_FRACTION}
OUT_FILES = ("t69.json", "t69.txt")
CHAMBER_BAR = 0.3  # above
MIN_SEEDS = 8
WHOLE_BAR = 0.0  # above
LIB_FLOOR = -0.3  # no valid library below
MIN_LIBRARIES = 2
STOP_TASK = "3"
STOP_FLOOR = -0.3  # library 3's mean below this ends the line
PERM_SEED, PERM_BINS, PERM_GROUPS = harness.LP_PERM_SEED, harness.LP_BINS, harness.LP_GROUPS
E18_PATH = feng_dirs.e18_stats_path()
E18_MANIFEST = E18_PATH.parent / "MANIFEST.json"
E18_RECORD = fre.RECORD

LABEL = ("try T69, a labelled in-sample diagnostic, never an adoption result: T1 candidate rc18, the outside atlas's "
         "direction over E9.5 to E18 at 1/8.5 (rc18), against m1rc_gmom (T62's candidate, R13's practice analogue), "
         "against its size control (rc18_sm) and its permuted control (rc18_perm), on row T62's five tasks; every task "
         "scores one truth (E10.5) over one day; the test spans two")
AHEAD_WORDS = "AHEAD"
NOT_AHEAD_WORDS = "NOT AHEAD: this ends rc18 (the atlas direction over E9.5 to E18) on T1"
STOP_WORDS = ("STOP: library 3's rc18 minus m1rc_gmom is below -0.3; this ends rc18 (the atlas direction over E9.5 to "
              "E18) on T1")


# ---------------------------------------------------------------------------------------------------------------------
# Guards


def source_guard(path: Path | None = None) -> list[str]:
    """T63's guard on this file (recomputed_read's guard, no import of a script forming the atlas direction on another
    route, no dynamic import)."""
    return ds.source_guard(Path(__file__ if path is None else path))


def row_line(row: str = ROW, path: Path | None = None) -> str | None:
    return fre.row_line(row, TRIES if path is None else path)


def e18_checks(manifest: Path | None = None, stats_path: Path | None = None, line: str | None = None,
               need_stats: bool = True, record_at_head: bytes | None = None) -> dict:
    """Row N58's guards: its manifest (written by feng_relabel_e18 from the frozen settings, E18 alone with N54's rows
    copied, those rows identical, authors_labels_read false, usable, not stopped, a clean gate of an N row), its copy
    committed at HEAD as the RECORD byte for byte, the statistics file's sha256 (hashed when need_stats) equal to the
    manifest's, and that sha256 named in row T69's line. Returns the values and the failed checks."""
    manifest = E18_MANIFEST if manifest is None else Path(manifest)
    reasons = []
    raw = manifest.read_bytes() if manifest.exists() else None
    man = json.loads(raw) if raw else {}
    if not man:
        reasons.append(f"no E18 labelling manifest at {manifest}")
    head = fre.fr.committed_at_head(E18_RECORD) if record_at_head is None else record_at_head
    if raw is None or head != raw:
        reasons.append("the E18 labelling's manifest is not committed at HEAD as its RECORD, byte for byte")
    if man.get("source") != fre.SOURCE:
        reasons.append(f"the E18 manifest's source is {man.get('source')!r}, not {fre.SOURCE!r}")
    if (man.get("config") or {}).get("sha256") != rc.YAML_SHA256:
        reasons.append("the E18 labelling ran under another settings file")
    if man.get("authors_labels_read") is not False:
        reasons.append(f"the E18 manifest's authors_labels_read is {man.get('authors_labels_read')!r}, not false")
    if man.get("stages") != list(fre.STAGES) or man.get("stages_copied") != list(fre.COPIED_STAGES):
        reasons.append(f"the E18 manifest's stages are {man.get('stages')} (copied {man.get('stages_copied')})")
    if (man.get("n54_stats") or {}).get("sha256") != rc.STATS_SHA256:
        reasons.append("the E18 file's copied rows are not from row N54's statistics")
    if (man.get("copied_rows") or {}).get("identical") is not True:
        reasons.append("the E18 file's E9.5 and E14.5 rows were not N54's byte for byte")
    if man.get("stopped") or man.get("usable") is not True:
        reasons.append("the E18 labelling stopped or failed the frozen yaml's unassigned cap (usable not true)")
    gate = man.get("gate") or {}
    if gate.get("clean") is not True or not fre.N_ROW.fullmatch(str(gate.get("row"))):
        reasons.append("the E18 labelling has no clean gate of an N row")
    sha = ((man.get("outputs") or {}).get(fre.STATS) or {}).get("sha256")
    hashed = None
    if need_stats:
        path = E18_PATH if stats_path is None else Path(stats_path)
        hashed = rc.sha256_file(path) if path.exists() else None
        if sha is None or hashed != sha:
            reasons.append(f"the E18 statistics file's sha256 is {hashed}, not the manifest's {sha}")
    line = row_line() if line is None else line
    if not sha or line is None or sha not in line:
        reasons.append(f"row {ROW}'s line does not name the E18 statistics' sha256 {sha}")
    return {"stats_sha256": sha, "stats_sha256_hashed": hashed, "yaml_sha256": (man.get("config") or {}).get("sha256"),
            "authors_labels_read": man.get("authors_labels_read"), "e18_row": gate.get("row"),
            "usable": man.get("usable"), "passes": not reasons, "reasons": reasons}


def refuse(out: Path, names: Sequence[str]) -> tuple[dict, dict, dict]:
    """The refusals before anything loads: the row, an existing output file, an unclean tree, the source guards, row
    N54's checks, row N58's checks and X002's entry for the E18 use."""
    if ROW not in br.tries_rows(TRIES):
        sys.exit(f"{ROW} is refused: it has no row in experiments/harness/TRIES.md (a row comes before its run)")
    gv.refuse_existing(out, tuple(names))
    started = gv.tree_state()
    if started["clean"] is not True:
        sys.exit(f"{ROW} is refused: the tree is not clean by vec.provenance's rule")
    bad = source_guard() + [f"recomputed_read.py: {x}" for x in rc.source_guard()]
    if bad:
        sys.exit(f"{ROW} is refused: {bad}")
    checks = rc.relabel_checks()
    if not checks["passes"]:
        sys.exit(f"{ROW} is refused: {checks['reasons']}")
    e18 = e18_checks()
    if not e18["passes"]:
        sys.exit(f"{ROW} is refused: {e18['reasons']}")
    if not fre.x002_names_this_use():
        sys.exit(f"{ROW} is refused: EXTERNAL_DATA.md's X002 row does not name {fre.X002_TOKEN!r}")
    return started, checks, e18


# ---------------------------------------------------------------------------------------------------------------------
# The directions (target-free)


def moved_mask(genes: Sequence[str], held: Sequence[str], *dirs: np.ndarray) -> np.ndarray:
    """The genes a direction step moves: finite in every one of `dirs` and not held flat (t1_scope.external_shift sets
    a NaN gene's and a held gene's log factor to 0)."""
    held_set = set(map(str, held))
    ok = np.array([g not in held_set for g in genes], dtype=bool)
    for d in dirs:
        ok &= np.isfinite(np.asarray(d, dtype=np.float64))
    return ok


def l2(v, mask: np.ndarray) -> float:
    v = np.asarray(v, dtype=np.float64)[mask]
    return float(np.sqrt(np.sum(v * v)))


def size_control(f, d18, mask: np.ndarray) -> np.ndarray:
    """rc18_sm's direction for one group: F times ||D18|| / ||F|| over `mask` (the moved genes), so that at rc18's
    fraction its applied log-factor size is rc18's. NaN stays NaN. A group whose F has no size is refused."""
    f = np.asarray(f, dtype=np.float64)
    size = l2(f, mask)
    if not size > 0:
        raise ValueError("F has no size on the moved genes: it cannot be scaled to rc18's")
    return f * (l2(d18, mask) / size)


def permuted_late(g, basis, group: str) -> np.ndarray:
    """rc18_perm's G' for one group: G's values permuted across the genes with finite G and basis, within PERM_BINS
    equal-count bins of the basis (vec.harness.t3_decile_permutation, row T15's procedure), seed [PERM_SEED, the
    group's index in PERM_GROUPS]. Every other gene keeps its own value (NaN stays NaN)."""
    if group not in PERM_GROUPS:
        raise ValueError(f"{group}: not one of {PERM_GROUPS}: it has no permutation seed (R2's guard)")
    g, basis = np.asarray(g, dtype=np.float64), np.asarray(basis, dtype=np.float64)
    out = g.copy()
    ok = np.isfinite(g) & np.isfinite(basis)
    out[ok] = harness.t3_decile_permutation(g[ok], basis[ok], [PERM_SEED, PERM_GROUPS.index(group)], PERM_BINS)
    return out


def arm_directions(late: Mapping, genes: Sequence[str], held: Sequence[str]) -> tuple[dict, dict]:
    """arm -> group -> the direction the build applies (at FRACTIONS[arm]), from `feng_dirs.e18_late_term`'s F, D18, G
    and basis; and per group the sizes (applied L2 over the moved genes) of every arm."""
    F, D18, G, basis = late["F"], late["D18"], late["G"], late["basis"]
    if set(F) != set(D18) or set(F) != set(G):
        raise ValueError(f"the groups of F {sorted(F)}, D18 {sorted(D18)} and G {sorted(G)} differ")
    dirs = {CANDIDATE: dict(F), RC18: dict(D18), SM: {}, PERM: {}}
    sizes = {}
    for grp in sorted(F):
        mask = moved_mask(genes, held, F[grp], D18[grp])
        dirs[SM][grp] = size_control(F[grp], D18[grp], mask)
        dirs[PERM][grp] = F[grp] + permuted_late(G[grp], basis[grp], grp)
        sizes[grp] = {a: FRACTIONS[a] * l2(dirs[a][grp], mask) for a in ARMS}
        sizes[grp]["moved_genes"] = int(mask.sum())
    return dirs, sizes


def direction_sha(dirs: Mapping[str, np.ndarray]) -> dict[str, str]:
    return {g: hashlib.sha256(np.ascontiguousarray(d, dtype=np.float64).tobytes()).hexdigest()
            for g, d in sorted(dirs.items())}


def load_kit(genes: Sequence[str], held: Sequence[str]) -> tuple[dict, dict, dict]:
    """From the files, through every check: the arms' directions on `genes`, their record (sizes, sha256 per arm and
    group, the E18 cells, the identity of the copied rows) and row R13's recomputed record. F from
    `rc.recomputed_inputs` must equal `feng_dirs.e18_late_term`'s F bit for bit."""
    F, ext = rc.recomputed_inputs(list(genes))
    late = feng_dirs.e18_late_term(list(genes))
    if direction_sha(F) != direction_sha(late["F"]):
        raise ValueError("F from row T62's route and from the E18 bundle differ")
    dirs, sizes = arm_directions(late, genes, held)
    record = {"sizes": sizes, "direction_sha256": {a: direction_sha(d) for a, d in dirs.items()},
              "e18_sha256": late["e18_sha256"], "e18_gate": late["e18_gate"], "cells": late["cells"],
              "copied_rows_identical": late["identical"].get("identical") is True,
              "fractions": dict(FRACTIONS)}
    return dirs, record, ext


# ---------------------------------------------------------------------------------------------------------------------
# The arms


def prediction_sha(stage) -> str:
    return ds.prediction_sha(stage)


def cut_arms(ctx: harness.BuildContext, factors, dirs: Mapping[str, Mapping], held_m1: Sequence[str],
             fractions: Mapping[str, float] = FRACTIONS) -> dict[str, Stage]:
    """Every arm of one seed on a harness cut, each T62's candidate build (`t1_scope.build_submission` with the
    candidate's g_mom hook and mask) on its own directions at its own fraction, finished as the harness finishes a
    candidate."""
    role = ctx.task.reference_role
    ref, src = harness._shallow(ctx.inputs[role]), ctx.task.role_stems[role]
    spec = t1_scope.load_spec()
    t1_scope.classify(ref, source=src, spec=spec)
    shares = harness.t1_chamber_shares(ref, src)
    n = harness.count_for(ctx.board, ref.n_obs)
    held = list(ctx.derived["t1_held_flat"])
    labels = ref.obs[t1_scope.LABEL_KEY].astype(str).unique()
    groups, covered, with_direction = rc.direction_kit(labels, dirs[CANDIDATE])
    hook = rc.rc_hook(factors, ctx.seed, held, with_direction, covered)

    def build(d, f) -> Stage:
        sub, _ = t1_scope.build_submission(ref, ctx.board, shares, n, ctx.seed, ref_source=src, post_cells=hook,
                                           held_flat=list(held_m1), directions=d, direction_fraction=f,
                                           direction_groups=groups)
        return harness.finish(ctx, ref, dense(sub.X), None)

    return {a: build(dirs[a], fractions[a]) for a in ARMS}


def library_arms(e95, half_shares: Mapping[str, float], factors, held: Sequence[str], held_m1: Sequence[str],
                 seed: int, within: Mapping, fb: Mapping[str, Mapping], dirs: Mapping[str, Mapping],
                 fractions: Mapping[str, float] = FRACTIONS, n_cells: int = rr.N_CELLS) -> dict:
    """Every arm of one seed on a library task, from T60's oracle draw (recomputed_read.library_pair's route:
    `gmom_combo.moved_by_feng` on stream [seed, 2], then the candidate's g_mom), as AnnData."""
    base, _ = t1_scope.reweight(e95, half_shares, n_cells, seed, within=within, fallback=fb["E9.5"])
    labels = e95.obs[t1_scope.LABEL_KEY].astype(str).unique()
    groups, covered, with_direction = rc.direction_kit(labels, dirs[CANDIDATE])
    hook = rc.rc_hook(factors, seed, held, with_direction, covered)
    return {a: hook(gc.moved_by_feng(base, dirs[a], fractions[a], groups, list(held_m1), [seed, 2])) for a in ARMS}


# ---------------------------------------------------------------------------------------------------------------------
# Runs


def _check_fraction(fraction: float, want: float, where: str) -> None:
    if abs(fraction - want) > FRACTION_TOL:
        sys.exit(f"{where}: the practice fraction is {fraction}, not {want}")


def _e18_fields(e18: Mapping) -> dict:
    return {k: e18[k] for k in ("stats_sha256", "stats_sha256_hashed", "yaml_sha256", "authors_labels_read",
                                "e18_row", "usable")}


def cmd_run_cut(a) -> None:
    out = RESULTS / a.out / "tasks"
    started, checks, e18 = refuse(out, (CUT_FILES[a.cut],))
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
    days = (harness.STAGE_DAYS[t.reference], harness.STAGE_DAYS[t.target])
    fractions = {CANDIDATE: feng_dirs.bracket_fraction(*days),
                 **{x: feng_dirs.bracket_fraction(*days, feng_dirs.E18_BRACKET_DAYS) for x in COUNTED}}
    for arm, f in fractions.items():
        _check_fraction(f, FRACTIONS[arm], f"{t.key} {arm}")
    held_m1, union = rc._held_m1(preregistered)
    dirs, record, ext = load_kit(list(inputs[ref_role].var_names), held_m1)
    factors = br.gmom_fits(inputs[ref_role], view.role_stems[ref_role], inputs["start"],
                           view.role_stems["start"])["file"]
    rows = []
    for s in seeds:
        t0 = time.time()
        f = gv.floor_anchors(t, board, target, ref_stage, s)
        ctx = harness.BuildContext(task=view, board=cboard, inputs=inputs, seed=harness.BUILD_SEED_OFFSET + s,
                                   derived={"t1_held_flat": preregistered})
        arms = cut_arms(ctx, factors, dirs, held_m1, fractions)
        rows.append({"seed": s, "floor": float(f.score),
                     "arms": gv.score_arms(board, target, ref_stage, s, f.anchors, arms),
                     "sha256": {k: prediction_sha(v) for k, v in arms.items()}, "seconds": round(time.time() - t0, 1)})
        del arms
        print(f"{t.key} seed {s}: {time.time() - t0:.0f} s", flush=True)  # no score
    res = {"label": LABEL, "row": ROW, "task": t.key, "seeds": list(seeds), "rows": rows, "fractions": fractions,
           "relabel": rc._relabel_fields(ext, checks), "e18": _e18_fields(e18), "kit": record,
           "constants": {"edit_categories": list(bt.G_MOM_EDIT_CATEGORIES), "held_m1_genes": len(held_m1)},
           "stamp": gv.stamp(started, [data / t.target, data / t.reference, *[data / f for f in t.inputs.values()],
                                       t1_scope.HELD_FLAT_GENES, t1_scope.CHAMBER_CLASSES, ROOT / ext["file"],
                                       E18_PATH, E18_MANIFEST, feng_dirs.RELABEL_CONFIG, rc.RECORD, E18_RECORD,
                                       ROOT / union["file"]])}
    out.mkdir(parents=True, exist_ok=True)
    (out / CUT_FILES[a.cut]).write_text(json.dumps(res, indent=1, default=float) + "\n")


def cmd_run_lib(a) -> None:
    out = RESULTS / a.out / "tasks"
    name = f"lib{a.library}.json"
    started, checks, e18 = refuse(out, (name,))
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
    fractions = {CANDIDATE: feng_dirs.bracket_fraction(9.5, 10.5),
                 **{x: feng_dirs.bracket_fraction(9.5, 10.5, feng_dirs.E18_BRACKET_DAYS) for x in COUNTED}}
    for arm, f in fractions.items():
        _check_fraction(f, FRACTIONS[arm], f"library {a.library} {arm}")
    held_m1, union = rc._held_m1(held)
    dirs, record, ext = load_kit(list(e95.var_names), held_m1)
    rows = []
    for s in seeds:
        t0 = time.time()
        pick = np.sort(np.random.default_rng([s, 99]).choice(e95.n_obs, rr.N_CELLS, replace=False))
        arms = library_arms(e95, shares, factors, held, held_m1, s, within, fb, dirs, fractions)
        stages = {LIB_REFERENCE: pt.to_stage(e95[pick]), **{k: pt.to_stage(v) for k, v in arms.items()}}
        rows.append({"seed": s, **rr.score_arms(board, stages, truth_stage, ref_stage, s),
                     "sha256": {k: prediction_sha(v) for k, v in stages.items()},
                     "seconds": round(time.time() - t0, 1)})
        del arms, stages
        print(f"library {a.library} seed {s}: {time.time() - t0:.0f} s", flush=True)  # no score
    res = {"label": LABEL, "row": ROW, "library": a.library, "seeds": list(seeds), "rows": rows,
           "fractions": fractions, "relabel": rc._relabel_fields(ext, checks), "e18": _e18_fields(e18), "kit": record,
           "constants": {"split_seed": rr.SPLIT_SEED, "n_cells": rr.N_CELLS, "board": rr.BOARD,
                         "edit_categories": list(bt.G_MOM_EDIT_CATEGORIES), "held_m1_genes": len(held_m1)},
           "stamp": gv.stamp(started, [*(data / f for f in rr.FILES.values()), t1_scope.HELD_FLAT_GENES,
                                       t1_scope.CHAMBER_CLASSES, ROOT / ext["file"], E18_PATH, E18_MANIFEST,
                                       feng_dirs.RELABEL_CONFIG, rc.RECORD, E18_RECORD, ROOT / union["file"]])}
    out.mkdir(parents=True, exist_ok=True)
    (out / name).write_text(json.dumps(res, indent=1, default=float) + "\n")


# ---------------------------------------------------------------------------------------------------------------------
# The reading


def task_files() -> dict[str, str]:
    return rc.task_files()


def committed_t62(name: str, tasks: Path | None = None) -> dict[int, float]:
    """seed -> m1rc_gmom's score in row T62's committed task file."""
    rows = json.loads(((T62_TASKS if tasks is None else tasks) / name).read_text())["rows"]
    return {int(r["seed"]): float(r["arms"][CANDIDATE]["score"]) for r in rows}


def guards(loaded: Mapping[str, Mapping | None], committed: Mapping[str, Mapping[int, float]], checks: Mapping,
           e18: Mapping, source_bad: Sequence[str], seeds: Sequence[int] = SEEDS) -> list[str]:
    """Every void condition of row T69 (its two gates among them); an empty list means none fired."""
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
        absent = sorted({a for x in rows for a in want if a not in x["arms"] or a not in (x.get("sha256") or {})})
        if absent or not rows:
            reasons.append(f"{k}: arms {absent} missing")
            continue
        stamp = r.get("stamp") or {}
        commits.add(stamp.get("commit"))
        if stamp.get("tree_clean") is not True:
            reasons.append(f"{k}: its tree was not clean at one commit")
        ref = committed.get(k, {})
        worst = max((abs(br.score_of(x, CANDIDATE) - ref[int(x["seed"])]) if int(x["seed"]) in ref else float("inf")
                     for x in rows), default=float("inf"))
        if not worst <= GATE_TOL:
            reasons.append(f"{k}: the gate fails: m1rc_gmom differs from row T62's committed m1rc_gmom by "
                           f"{worst:.3g} (tolerance {GATE_TOL:g})")
        fr_ = r.get("fractions") or {}
        if any(not abs(float(fr_.get(a, float("nan"))) - FRACTIONS[a]) <= FRACTION_TOL for a in ARMS):
            reasons.append(f"{k}: the practice fractions are {fr_}, not {FRACTIONS}")
        kit = r.get("kit") or {}
        if kit.get("copied_rows_identical") is not True:
            reasons.append(f"{k}: the gate fails: the E9.5 and E14.5 statistics were not N54's byte for byte")
        rel = r.get("relabel") or {}
        if rel.get("sha256") != rc.STATS_SHA256 or rel.get("stats_sha256_hashed") != rc.STATS_SHA256:
            reasons.append(f"{k}: ran on another recomputed statistics file")
        if rel.get("yaml_sha256") != rc.YAML_SHA256:
            reasons.append(f"{k}: ran under another settings file")
        if rel.get("route") != rc.GROUPS_OPTION or rel.get("record_route") != rc.GROUPS_OPTION:
            reasons.append(f"{k}: the route was not recomputed")
        if rel.get("authors_labels_read") is not False:
            reasons.append(f"{k}: the relabelling record did not say authors_labels_read false")
        e = r.get("e18") or {}
        want_sha = e18.get("stats_sha256")
        if not want_sha or e.get("stats_sha256") != want_sha or e.get("stats_sha256_hashed") != want_sha:
            reasons.append(f"{k}: ran on E18 statistics other than the row's")
        if kit.get("e18_sha256") != want_sha:
            reasons.append(f"{k}: its directions came from E18 statistics other than the row's")
        if e.get("yaml_sha256") != rc.YAML_SHA256:
            reasons.append(f"{k}: the E18 labelling ran under another settings file")
        if e.get("authors_labels_read") is not False:
            reasons.append(f"{k}: the E18 record did not say authors_labels_read false")
    if len(commits) > 1:
        reasons.append(f"the task files ran at more than one commit: {sorted(map(str, commits))}")
    reasons += [f"at the read: {x}" for x in checks.get("reasons", [])]
    reasons += [f"at the read: {x}" for x in e18.get("reasons", [])]
    reasons += [f"source: {x}" for x in source_bad]
    return reasons


def beats_both(rows: Sequence[Mapping]) -> dict:
    """On one task, the seeds where rc18 scores above rc18_sm and above rc18_perm (both on the same seed), and each
    count alone (reported)."""
    s = {a: np.array([br.score_of(r, a) for r in rows]) for a in (RC18, SM, PERM)}
    over_sm, over_perm = s[RC18] > s[SM], s[RC18] > s[PERM]
    return {"both": int((over_sm & over_perm).sum()), "over_sm": int(over_sm.sum()),
            "over_perm": int(over_perm.sum()), "of": len(rows)}


def t69_reading(loaded: Mapping[str, Mapping | None], committed: Mapping, checks: Mapping, e18: Mapping,
                source_bad: Sequence[str] = ()) -> dict:
    """Row T69's statement, applied to the task files."""
    reasons = guards(loaded, committed, checks, e18, source_bad)
    if reasons:
        return {"outcome": "VOID", "reasons": reasons}
    rows = {k: r["rows"] for k, r in loaded.items()}
    vs = {k: br.paired(v, RC18, CANDIDATE) for k, v in rows.items()}
    valid = [lib for lib in LIBRARIES if rc.library_valid(loaded[lib])]
    low = [lib for lib in valid if vs[lib]["mean"] < LIB_FLOOR]
    both = beats_both(rows["chamber"])
    parts = {"chamber": bool(vs["chamber"]["mean"] > CHAMBER_BAR and vs["chamber"]["ahead"] >= MIN_SEEDS),
             "whole": bool(vs["whole"]["mean"] > WHOLE_BAR),
             "two_valid_libraries": len(valid) >= MIN_LIBRARIES,
             "no_valid_library_below_floor": not low,
             "beats_both_controls": both["both"] >= MIN_SEEDS}
    stop = bool(vs[STOP_TASK]["mean"] < STOP_FLOOR)
    ahead = bool(all(parts.values()) and not stop)
    outcome = STOP_WORDS if stop else AHEAD_WORDS if ahead else NOT_AHEAD_WORDS
    reported = {"sm_minus_m1rc": {k: br.paired(v, SM, CANDIDATE) for k, v in rows.items()},
                "perm_minus_m1rc": {k: br.paired(v, PERM, CANDIDATE) for k, v in rows.items()},
                "rc18_minus_sm": {k: br.paired(v, RC18, SM) for k, v in rows.items()},
                "rc18_minus_perm": {k: br.paired(v, RC18, PERM) for k, v in rows.items()},
                "beats_both_by_task": {k: beats_both(v) for k, v in rows.items()},
                "mean_scores": {k: {a: float(np.mean([br.score_of(x, a) for x in v]))
                                    for a in (LIB_ARMS if k in LIBRARIES else ARMS)} for k, v in rows.items()},
                "inverted_anchor_libraries": [lib for lib in LIBRARIES if lib not in valid],
                "sizes": (loaded["chamber"].get("kit") or {}).get("sizes")}
    return {"outcome": outcome, "ahead": ahead, "stop": stop, "parts": parts, "valid_libraries": valid,
            "libraries_below_floor": low, "chamber_beats_both": both, "vs_m1rc": vs, "reported": reported,
            "j6_range": sorted([0.0, vs["chamber"]["mean"]])}


def outcome_lines(r: Mapping) -> list[str]:
    if r["outcome"] == "VOID":
        return ["Outcome: VOID (" + "; ".join(r["reasons"]) + "); no rerun under row T69 and its tries count"]
    lines = [f"Outcome: {r['outcome']}", "Parts: " + "; ".join(f"{k} {'holds' if v else 'fails'}"
                                                             for k, v in r["parts"].items())]
    lines.append("What follows: one priced question to the captain (rc18 beside R13's file)" if r["ahead"] else
                 "What follows: nothing is built; no question, no rerun under another row")
    return lines


def _line(name: str, p: Mapping) -> str:
    g = ", ".join(f"{x} {v:+.2f}" for x, v in p["groups"].items())
    return (f"  {name}: {p['mean']:+.2f} points, sd {p['sd']:.2f}, ahead {p['ahead']} of {p['of']}, per seed "
            f"{p['range'][0]:+.2f} to {p['range'][1]:+.2f}; by group: {g}")


def t69_text(res: Mapping) -> str:
    r = res["reading"]
    lines = [res["label"], *outcome_lines(r)]
    if r["outcome"] != "VOID":
        rep = r["reported"]
        for title, per in (("rc18 - m1rc_gmom", r["vs_m1rc"]), ("reported: rc18 - rc18_sm", rep["rc18_minus_sm"]),
                           ("reported: rc18 - rc18_perm", rep["rc18_minus_perm"]),
                           ("reported: rc18_sm - m1rc_gmom", rep["sm_minus_m1rc"]),
                           ("reported: rc18_perm - m1rc_gmom", rep["perm_minus_m1rc"])):
            lines.append(title)
            lines += [_line(k, p) for k, p in per.items()]
        b = r["chamber_beats_both"]
        lines.append(f"Chamber cut: rc18 above both controls on {b['both']} of {b['of']} seeds (above rc18_sm "
                     f"{b['over_sm']}, above rc18_perm {b['over_perm']})")
        lines.append(f"Valid libraries: {', '.join(r['valid_libraries']) or 'none'}; below {LIB_FLOOR}: "
                     f"{', '.join(r['libraries_below_floor']) or 'none'}; library {STOP_TASK} stop: "
                     f"{'fired' if r['stop'] else 'not fired'}")
        lines.append("Mean scores: " + "; ".join(f"{k}: " + ", ".join(f"{a} {v:.2f}" for a, v in m.items())
                                                 for k, m in rep["mean_scores"].items()))
        lines.append(f"Applied step sizes per group (target-free): {rep['sizes']}")
        lines.append(f"J6 range (horizon ratio 2; 0 to the chamber cut's gain): {r['j6_range'][0]:+.2f} to "
                     f"{r['j6_range'][1]:+.2f}")
    lines.append(f"Tries counted: {res['tries']}")
    lines.append(f"Commit {res['stamp']['commit']}, tree clean {res['stamp']['tree_clean']}, UTC {res['stamp']['utc']}")
    return "\n".join(lines) + "\n"


def read(out: Path) -> dict:
    if ROW not in br.tries_rows(TRIES):
        sys.exit(f"{ROW} is refused: it has no row in experiments/harness/TRIES.md")
    gv.refuse_existing(out, OUT_FILES)
    started = gv.tree_state()
    paths = {k: out / "tasks" / f for k, f in task_files().items()}
    loaded = {k: json.loads(p.read_text()) if p.exists() else None for k, p in paths.items()}
    committed = {k: committed_t62(f) for k, f in task_files().items()}
    checks = rc.relabel_checks(need_stats=False)  # the stats files' hashes are in each task file
    e18 = e18_checks(need_stats=False)
    bad = source_guard() + [f"recomputed_read.py: {x}" for x in rc.source_guard()]
    reading = t69_reading(loaded, committed, checks, e18, bad)
    res = {"label": LABEL, "row": ROW, "reading": reading, "tries": TRIES_COUNT}
    res["stamp"] = gv.stamp(started, [p for p in paths.values() if p.exists()])
    out.mkdir(parents=True, exist_ok=True)
    (out / OUT_FILES[0]).write_text(json.dumps(res, indent=1, default=float) + "\n")
    (out / OUT_FILES[1]).write_text(t69_text(res))
    print("\n".join(outcome_lines(reading)) + f"\n(numbers in {out / OUT_FILES[1]})")
    return res


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("run-cut", help="one E10.5 cut of try T69")
    p.add_argument("--cut", required=True, choices=list(CUTS))
    p.add_argument("--out", required=True)
    p.add_argument("--seeds", default=",".join(map(str, SEEDS)))
    p.add_argument("--data", default=str(harness.DATA_DIR))
    p = sub.add_parser("run-lib", help="one library task of try T69")
    p.add_argument("--library", required=True, choices=list(LIBRARIES))
    p.add_argument("--out", required=True)
    p.add_argument("--seeds", default=",".join(map(str, SEEDS)))
    p.add_argument("--data", default=str(harness.DATA_DIR))
    p = sub.add_parser("read", help="try T69's reading")
    p.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    if a.cmd == "read":
        read(RESULTS / a.out)
        return 0
    {"run-cut": cmd_run_cut, "run-lib": cmd_run_lib}[a.cmd](a)
    return 0


if __name__ == "__main__":
    sys.exit(main())
