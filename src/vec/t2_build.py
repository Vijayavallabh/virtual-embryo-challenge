"""T2 submission utilities: the submitted cell count, random sampling, coordinate helpers, a cloud guard, a builder.

What each piece enforces (docs/PLAN.md sections 3.6, 4.2 and 4.3; experiments/density/REPORT.md and
experiments/coords/REPORT.md, both measured on released data and independently reviewed):

- `n_cells`: submit max_cells, never below truth_cells when the board publishes it, always within [min_cells,
  max_cells] (the limits are inclusive). Below 1x truth_cells every seed loses (occupancy_dice thins the truth to the
  submitted count); above it imperfect predictions keep gaining. Keying the count to max_cells means truth_cells only
  bounds it, never chooses it. The release-day density rerun may replace max_cells by a lower "<f>max" level under the
  rule pre-registered in PLAN 3.6; that level is passed as `override`. `count_override` applies that rule (approved
  answer 7 of docs/RELEASE_RULES_SUMMARY.md, 2026-10-21) to one board's rerun scores and returns the override;
  experiments/density/summarise.py --decide writes it per test board to count_decision.json, and
  `n_cells_from_decision` reads that file for a test board (`read_count_decision`) and passes its override to n_cells,
  or, when that level falls below the board's lower limit, the next passing count up: the passing level with the
  fewest cells that the board allows, else max_cells (COUNT_FLOOR_RULE, approved by the user the same day, who left the
  reading of "next" to the AI assistant).
- `sample_cells`: random rows only. Released files store rows section by section, so a contiguous block of rows is a
  slab of sections (16-31 points lost at 1x on a near-perfect prediction, 13-18 at 8x). Rows repeat only when the
  count exceeds the pool; exact duplicate coordinates are then broken by isotropic jitter of at most 0.25 x the median
  nearest-neighbour distance, the budget measured to cost nothing. The rows are drawn by the shared core
  (vec.core.draw_balanced, docs/PLAN.md section 5; importable from here under its old name `_balanced_rows`); the
  count rule, the jitter, the coordinate helpers and the guard are this module's, the T2 head.
- Coordinates stay native: `centre`, `proper_rotation` (det +1 only: a mirrored copy of a correct shape loses up to
  0.12 Dice under the pinned scorer, so improper matrices are refused, never repaired), `rescale_to_rms` (an isotropic
  scaling about the centroid) and `apply_size_procedure`, which reads the pre-registered size procedure's output JSON
  (experiments/coords/size_procedure.py) and rescales in the z mode recorded there. z is never rescaled across stages
  except as that recorded mode says. The heart pipeline calls `size_or_hold` (approved answer 3): heart cells must come
  from E10.5, or the size rule goes back to the captain, and the caller's reference cloud must be the cloud the size
  JSON measured as E10.5, valid or not; then apply_size_procedure for a valid JSON, and for a heart JSON with
  valid=false the fallback (the cells keep E10.5's size).
- `cloud_guard`: RMS radius and occupied-voxel coverage on the scorer's 16^3 canonical grid, against same-count random
  subsets of the source. A random 1x subset covers only 81-83% of its source's occupied voxels and a contiguous block
  21-25%, so the guard compares with the random reference, never with 100%.
- `build_submission`: X in the board's gene order and obsm["spatial_3D"], written and passed through
  vec.io.check_submission, with a warning when X does not look like log1p CP10k (values above log1p(1e4), whole
  numbers only, or a median implied library size outside [0.1, 2] x 1e4).

No size, composition or class list of any stage is written here: every number about a target arrives as an argument
(the board spec, the size procedure's JSON, the reference cloud of the size fallback); the approved rules' own
thresholds are not properties of a stage. tests/test_t2_build.py checks the module's numeric literals against an
allow-list.

    pool = load_stage(path, genes=board.genes, drop_duplicates=True)     # prediction input: synthetic records dropped
    n, count_record = n_cells_from_decision(board, pool.n_obs, "count_decision.json")  # rule 7, or COUNT_FLOOR_RULE
    s = sample_cells(pool, n, seed=0)
    coords, size_record = size_or_hold(s.coords, "size_heart_E12.5.json", "H1", source_stage="E10.5",
                                       reference_cloud=pool, expect_board="heart")
    guard = cloud_guard(s.coords, pool.obsm["spatial_3D"], seed=0)        # same frame and scale as the pool
    adata, report = build_submission(board, pool.X[s.rows], coords, genes=list(pool.var_names), out="s1.h5ad")
"""

from __future__ import annotations

import hashlib
import json
import re
import tempfile
from collections.abc import Mapping
from dataclasses import asdict, dataclass, field
from pathlib import Path

import anndata as ad
import numpy as np
from scipy import sparse
from scipy.spatial import cKDTree

from vec.boards import Board
from vec.core import barycentric, coupling_rotation, draw_balanced, fgw_coupling, herd_rows, nearest_anchors
from vec.io import (
    COORD_KEY,
    TARGET_LIBRARY_SIZE,
    CheckReport,
    check_submission,
    dense,
    make_submission,
    write_submission,
)
from vec.provenance import stamp_path

JITTER_NN_FRAC = 0.25  # free jitter budget, as a fraction of the median nearest-neighbour distance (density)
JITTER_MIN_FRAC = 0.5  # a broken duplicate moves 0.5-1 x the budget: clearly apart, never beyond the budget
GRID_N = 16  # veckit occupancy_dice: 16^3 voxels spanning [-GRID_EXTENT, GRID_EXTENT] RMS radii
GRID_EXTENT = 3.0
VOXEL_EDGE_EPS = 1e-9  # veckit clips canonical coordinates just below the upper grid edge
PROPER_FLIPS = ((1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1))  # axis sign flips with det +1, as in veckit
GUARD_COVER_RATIO = 0.9  # flag below 0.9 x the random reference's coverage (contiguous blocks: about 0.3 at 1x)
GUARD_RMS_LOG_TOL = 0.03  # flag |log(RMS ratio / expected ratio)| above this (head blocks: 0.64-0.94 of the RMS)
GUARD_OUTSIDE_TOL = 0.1  # flag when the share of occupied voxels outside the source's exceeds the reference's by this
GUARD_REF_DRAWS = 5
GUARD_FRAMES = ("source", "own")
MIN_CLOUD = 4  # fewest points for a PCA frame
ROTATION_TOL = 1e-6
DEGENERATE_RMS = 1e-12
NATIVE_Z_WARN_FACTOR = 2.0  # z-hold mode: warn when the cloud's z SD and the reference stage's differ by more
LIBRARY_BAND = (0.1, 2.0)  # builder: warn when the median implied library size leaves this band, in units of 1e4
SIZE_BOARDS = ("heart", "embryo")
SIZE_HYPOTHESES = ("H1", "H0", "H2")
APPLICABLE_Z_MODES = ("raw", "common_step", "xy_grow_z_hold")  # "split" is logged by the procedure for diagnosis only
SIZE_MATCH = ("full", "xy")
REDACTED_PREFIX = "[redacted"
_MAX_LABEL = re.compile(r"^(?P<f>\d+(?:\.\d+)?|\.\d+)max$")  # density rerun level names: "0.75max", "1max"

# Release-day rules approved by the user on 2026-10-21 (docs/RELEASE_RULES_SUMMARY.md, approved answers 3 and 7).
RULES_APPROVED = "2026-10-21"
# Answer 3 (`size_or_hold`): heart cells must come from this stage, and keep its size when the size JSON is invalid.
SIZE_FALLBACK_STAGE = "E10.5"
SIZE_FALLBACK_RULE = ("approved rule 3: if the size procedure returns valid=false while E10.5 is usable, keep E10.5's "
                      "size, so the scale term sits at the floor instead of the build failing; s1 and s2 still differ "
                      "by scope. If E10.5 is barred (A20 fails), the heart plan is rebuilt from E9.5 and the size rule "
                      "comes back to the captain")
SIZE_CROSSCHECK_LOG_TOL = 1e-6  # float noise: the cloud the procedure measured gives its RMS radius; others refused
REASON_KEYS = ("reason", "reasons", "invalid_reason", "note", "notes")  # reason fields a size JSON may carry
# Answer 7, the count rule (`count_override`), on the density rerun (experiments/density/run_density.py).
COUNT_RULE_CONDITION = "shuffle50"  # the intermediate-quality condition: half the expression rows shuffled
COUNT_RULE_BASE = "1max"  # max_cells
COUNT_RULE_LEVELS = ("0.25max", "0.5max", "0.75max")  # the lower levels the rerun tries (its MAX_FRACS below 1)
COUNT_RULE_SEEDS = tuple(range(10))  # the 10 release seeds, 0-9 (docs/P3_RUNBOOK.md step 4.5)
COUNT_RULE_MIN_GAIN = 0.3  # board points: a level's paired mean gain over max_cells must exceed this
COUNT_RULE_MIN_AHEAD = 8  # of the 10 seeds, how many the level must beat max_cells in (difference > 0)
# each test board's kind -> the validation target (experiments/density/run_density.py TARGETS) whose rerun decides it
COUNT_RULE_TARGETS = {"heart": "heart_E10.5", "embryo": "embryo_E7.5"}
COUNT_RULE = ("approved rule 7: a lower level replaces max_cells on a board only if, on that board's validation "
              "target, its paired mean gain over max_cells exceeds 0.3 points and it beats max_cells in at least 8 of "
              "the 10 seeds. If several pass, the largest mean gain wins (ties go to the higher count). Heart and "
              "embryo are decided separately, each for its own test board. If the rerun cannot run, keep max_cells")
# The case answer 7 leaves open (`n_cells_from_decision`). The user approved the recommendation in the quoted words on
# 2026-10-21 and left the reading of "next" to the AI assistant ("anything that would improve my chances of winning");
# a decision panel chose the next count up (experiments/density/floor_reading/README.md).
COUNT_FLOOR_RULE = ("the case rule 7 leaves open, approved by the user on 2026-10-21 as 'the next passing count the "
                    "board allows, or max_cells if none does', read as the next count up (the reading the user left "
                    "to the AI assistant): if the winning level gives a test board fewer cells than its lower limit "
                    "(truth_cells where published and at least min_cells, else min_cells; a bound only), use the level "
                    "that passed rule 7 with the fewest cells at or above that limit, whatever its gain; if none, keep "
                    "max_cells")


# ---------------------------------------------------------------------------------------------------------------------
# Inputs
def _as_coords(C, name: str = "coords") -> np.ndarray:
    C = np.asarray(C, dtype=np.float64)
    if C.ndim != 2 or C.shape[1] != 3 or C.shape[0] == 0:
        raise ValueError(f"{name} must be a non-empty (n, 3) array, got shape {C.shape}")
    if not np.isfinite(C).all():
        raise ValueError(f"{name} contains NaN or infinite values")
    return C


def _coords_of(stage) -> np.ndarray:
    """(n, 3) float64 coordinates of an AnnData (obsm["spatial_3D"], first three columns as the scorer reads them), a
    vec.evaluate.Stage, or a plain array."""
    if isinstance(stage, ad.AnnData):
        if COORD_KEY not in stage.obsm:
            raise KeyError(f"obsm[{COORD_KEY!r}] is required")
        return _as_coords(np.asarray(stage.obsm[COORD_KEY], dtype=np.float64)[:, :3])
    if not isinstance(stage, np.ndarray) and hasattr(stage, "coords"):
        if stage.coords is None:
            raise ValueError("the stage has no coordinates")
        return _as_coords(np.asarray(stage.coords, dtype=np.float64)[:, :3])
    return _as_coords(stage)


def _whole(v, name: str, minimum: int) -> int:
    if isinstance(v, bool) or not isinstance(v, int | np.integer) or v < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}, got {v!r}")
    return int(v)


def rms_radius(coords) -> float:
    """sqrt(mean squared distance to the centroid) over all three coordinates: the scorer's RMS radius
    (veckit scale_log_ratio and the size procedure use the same definition)."""
    C = _as_coords(coords)
    X = C - C.mean(0)
    return float(np.sqrt((X ** 2).sum(1).mean()))


def median_nn_distance(coords) -> float:
    """Median distance from each distinct position to its nearest other distinct position (exact repeats count once,
    so copies do not pull it to zero); nan for fewer than two distinct positions."""
    U = np.unique(_as_coords(coords), axis=0)
    if U.shape[0] < 2:
        return float("nan")
    d, _ = cKDTree(U).query(U, k=2)
    return float(np.median(d[:, 1]))


# ---------------------------------------------------------------------------------------------------------------------
# Cell count
def _override_count(override, max_cells: int) -> int:
    if isinstance(override, str):
        m = _MAX_LABEL.match(override)
        if not m:
            raise ValueError(f"override {override!r}: use an integer count or a '<f>max' label such as '0.75max'")
        f = float(m["f"])
        if not 0 < f <= 1:
            raise ValueError(f"override {override!r}: the fraction of max_cells must be in (0, 1]")
        return int(round(f * max_cells))
    if isinstance(override, bool) or not isinstance(override, int | np.integer):
        raise TypeError(f"override must be an integer count or a '<f>max' label, got {type(override).__name__}")
    return int(override)


def _limits(board: Board) -> tuple[int, int, int | None, int, str]:
    """A board's (max_cells, min_cells, truth_cells, floor, floor source): the floor is truth_cells when the board
    publishes it (at least min_cells), else min_cells. Refused: a board without max_cells, a floor above max_cells."""
    if board.max_cells is None:
        raise ValueError(f"{board.key}: the board publishes no max_cells; the count rule needs it")
    max_cells = _whole(board.max_cells, "max_cells", 1)
    min_cells = _whole(board.min_cells, "min_cells", 1)
    truth = None if board.truth_cells is None else _whole(board.truth_cells, "truth_cells", 1)
    if truth is not None and truth >= min_cells:
        floor, floor_source = truth, "truth_cells"
    else:
        floor, floor_source = min_cells, "min_cells"
    if floor > max_cells:
        raise ValueError(f"{board.key}: inconsistent limits, {floor_source} {floor} > max_cells {max_cells}")
    return max_cells, min_cells, truth, floor, floor_source


def n_cells(board: Board, pool_size: int, override: int | str | None = None) -> tuple[int, dict]:
    """The number of cells to submit on a board, and a record of why.

    Rule (docs/PLAN.md 3.6): max_cells, never below truth_cells when the board publishes it, always within
    [min_cells, max_cells]. `override` replaces max_cells only as the release-day density rule allows: an integer
    count, or a "<f>max" label (round(f x max_cells), the density rerun's level names, which carry over from the
    validation boards to the test boards). An override outside [floor, max_cells] is refused, never clipped.
    `pool_size` (distinct cells available) only annotates the record: a count above it is drawn by `sample_cells`
    with rows repeated as evenly as possible, a regime the density experiment did not test (it drew at most 44% of a
    pool).
    """
    max_cells, min_cells, truth, floor, floor_source = _limits(board)
    pool_size = _whole(pool_size, "pool_size", 1)
    if override is None:
        n, rule = max_cells, "max_cells"
    else:
        n, rule = _override_count(override, max_cells), f"override {override!r}"
        if not floor <= n <= max_cells:
            raise ValueError(f"{board.key}: override {override!r} gives {n} cells, outside [{floor} ({floor_source}), "
                             f"{max_cells} (max_cells)]; refused, not clipped")
    record = {
        "board": board.key, "n": n, "rule": rule, "max_cells": max_cells, "min_cells": min_cells,
        "truth_cells": truth, "floor": floor, "floor_source": floor_source,
        "fraction_of_max_cells": n / max_cells, "pool_size": pool_size, "exceeds_pool": n > pool_size,
        "repeated_rows": max(0, n - pool_size), "x_pool": n / pool_size,
    }
    if n > pool_size:
        record["note"] = (f"count above the pool: {n - pool_size} rows repeat (each pool row {n // pool_size} or "
                          f"{-(-n // pool_size)} times) and their copies are jittered apart; untested by the density "
                          "experiment")
    return n, record


def _seed(s) -> int:
    if isinstance(s, bool) or not isinstance(s, int | np.integer):
        raise TypeError(f"seeds must be integers, got {s!r}")
    return int(s)


def _score_or_none(v, level: str, seed: int) -> float | None:
    """A finite score, or None for a missing one (None, NaN or infinite)."""
    if v is None:
        return None
    if isinstance(v, bool) or not isinstance(v, int | float | np.integer | np.floating):
        raise TypeError(f"{level} seed {seed}: a score must be a number (None or NaN when missing), got {v!r}")
    v = float(v)
    return v if np.isfinite(v) else None


def count_override(scores: Mapping) -> tuple[str | None, dict]:
    """Approved rule 7: the count level for one test board, from the density rerun on that board's validation target.

    scores: {level: {seed: score}}, the emulator board scores of the rerun's intermediate-quality condition (shuffle50
    in experiments/density/run_density.py: random real cells, half their expression rows shuffled) at the max_cells
    levels, keyed by the rerun's level labels ("1max" is max_cells; "0.25max", "0.5max" and "0.75max" the lower levels)
    and by seed. Decide each board separately (COUNT_RULE_TARGETS): heart_E10.5 for the heart test board, embryo_E7.5
    for the embryo one.

    The rule (docs/RELEASE_RULES_SUMMARY.md, approved answer 7, 2026-10-21): a lower level replaces max_cells only if
    its paired mean gain over max_cells (the mean over seeds of its score minus the 1max score of the same seed)
    exceeds 0.3 points and it beats max_cells (a difference above 0) in at least 8 of the 10 seeds. If several pass, the
    largest mean gain wins, an exact tie going to the higher count. If the rerun cannot run, keep max_cells: the rule
    needs the 10 release seeds (0-9) paired between 1max and every lower level, so a missing level, seed or 1max row
    returns None with the gap recorded. A None, NaN or infinite score is a missing seed. Other seeds and levels are
    ignored and listed; a label that is not "<f>max", a seed that is not an integer or a score that is not a number is
    refused.

    Returns (override, record): override is None (keep max_cells) or the winning level's label, which
    `n_cells_from_decision` passes to `n_cells` on the matching test board; a level that gives that board fewer cells
    than its lower limit gives way to the passing level with the fewest cells that the board allows (the next count
    up), else max_cells (COUNT_FLOOR_RULE), and n_cells refuses any count outside a board's limits, never clips it.
    The record holds each level's paired seeds, differences, mean gain, seeds ahead and verdict, the passing levels and
    the decision.
    """
    if not isinstance(scores, Mapping):
        raise TypeError("scores must map each count level to {seed: score}")
    for level, per_seed in scores.items():
        if not isinstance(level, str) or not _MAX_LABEL.match(level):
            raise ValueError(f"count level {level!r}: use the density rerun's '<f>max' labels, such as '0.75max'")
        if not isinstance(per_seed, Mapping):
            raise TypeError(f"{level}: scores must map each seed to a score")
    release = set(COUNT_RULE_SEEDS)
    finite: dict[str, dict[int, float]] = {}
    ignored_seeds: set[int] = set()
    for level in (COUNT_RULE_BASE, *COUNT_RULE_LEVELS):
        finite[level] = {}
        for seed, score in scores.get(level, {}).items():
            seed = _seed(seed)
            if seed not in release:
                ignored_seeds.add(seed)
            elif (v := _score_or_none(score, level, seed)) is not None:
                finite[level][seed] = v
    base = finite[COUNT_RULE_BASE]
    gaps = []
    for level in (COUNT_RULE_BASE, *COUNT_RULE_LEVELS):
        if not finite[level]:
            gaps.append(f"no {level} scores")
        elif missing := sorted(release - set(finite[level])):
            gaps.append(f"{level} lacks seeds {missing}")
    complete = not gaps
    per_level = {}
    for level in COUNT_RULE_LEVELS:
        paired = sorted(set(finite[level]) & set(base))
        diffs = [finite[level][s] - base[s] for s in paired]
        gain = float(np.mean(diffs)) if diffs else None
        ahead = sum(d > 0 for d in diffs)
        passes = bool(gain > COUNT_RULE_MIN_GAIN and ahead >= COUNT_RULE_MIN_AHEAD) if complete else None
        per_level[level] = {"paired_seeds": paired, "differences": diffs, "mean_gain": gain, "seeds_ahead": ahead,
                            "passes": passes}
    passing = [level for level in COUNT_RULE_LEVELS if per_level[level]["passes"]]
    override = (max(passing, key=lambda lv: (per_level[lv]["mean_gain"], float(_MAX_LABEL.match(lv)["f"])))
                if passing else None)
    n_seeds = len(COUNT_RULE_SEEDS)
    if not complete:
        decision = (f"keep max_cells: the rerun is incomplete ({'; '.join(gaps)}), and approved rule 7 keeps max_cells "
                    "when the rerun cannot run")
    elif override is None:
        decision = (f"keep max_cells: no lower level has a paired mean gain above {COUNT_RULE_MIN_GAIN:g} points with "
                    f"at least {COUNT_RULE_MIN_AHEAD} of the {n_seeds} seeds ahead")
    else:
        win = per_level[override]
        decision = (f"use {override}: paired mean gain {win['mean_gain']:+.3f} points over max_cells, "
                    f"{win['seeds_ahead']} of the {n_seeds} seeds ahead")
        if len(passing) > 1:
            decision += f"; the largest mean gain of the passing levels {passing} (ties go to the higher count)"
    record = {
        "rule": COUNT_RULE, "approved": RULES_APPROVED, "condition": COUNT_RULE_CONDITION, "base": COUNT_RULE_BASE,
        "levels": list(COUNT_RULE_LEVELS), "seeds": list(COUNT_RULE_SEEDS), "min_mean_gain": COUNT_RULE_MIN_GAIN,
        "min_seeds_ahead": COUNT_RULE_MIN_AHEAD, "complete": complete, "incomplete": gaps, "base_seeds": sorted(base),
        "per_level": per_level, "passing": passing, "override": override, "decision": decision,
        "ignored_levels": sorted(set(scores) - {COUNT_RULE_BASE, *COUNT_RULE_LEVELS}),
        "ignored_seeds": sorted(ignored_seeds),
    }
    return override, record


def _t2_kind(board: Board, kind: str | None) -> str:
    """A T2 test board's kind ("heart" or "embryo"): the one its key names; `kind`, when given, must agree, and is
    needed only for a key that names neither."""
    if board.task != "T2":
        raise ValueError(f"{board.key}: approved rule 7 sets the count of T2 boards only")
    if kind is not None and kind not in COUNT_RULE_TARGETS:
        raise ValueError(f"kind must be one of {sorted(COUNT_RULE_TARGETS)}, got {kind!r}")
    named = {p for p in board.key.split(":") if p in COUNT_RULE_TARGETS}
    if len(named) > 1 or (kind is not None and named and named != {kind}):
        raise ValueError(f"{board.key}: the key names {sorted(named)}, which does not agree with kind {kind!r}")
    if not named and kind is None:
        raise ValueError(f"{board.key}: the key names neither heart nor embryo; pass kind='heart' or kind='embryo'")
    return kind if kind is not None else named.pop()


def read_count_decision(decision_json, board: Board, kind: str | None = None) -> tuple[str | None, dict]:
    """The count level that approved rule 7's decision file gives a T2 test board, and a record of where it came from.

    decision_json: the count_decision.json written by experiments/density/summarise.py --decide, or its dict. The
    board's kind ("heart" or "embryo") is the one its key names (`kind`, when given, must agree; it is needed only for
    a key that names neither), so a heart board never takes the embryo decision. Refused: a file that is not a decision
    of this rule (another rule text or approval date), a kind the file does not hold, a decision whose validation target
    is not the kind's (COUNT_RULE_TARGETS), an override that is neither None (keep max_cells) nor one of the rule's
    lower levels, and for an override, a record of passing levels it does not head (`_passing_ranking`). Returns
    (override, record): the record holds the file's path and sha256 (None for a dict), the kind, the validation target,
    the override, the passing levels in rule order (`ranking`, empty for None; COUNT_FLOOR_RULE chooses among them by
    count), the decision text and the git state the decision recorded.
    """
    kind = _t2_kind(board, kind)
    if isinstance(decision_json, Mapping):
        d, path, sha = decision_json, None, None
    else:
        path = Path(decision_json)
        raw = path.read_bytes()
        d, sha = json.loads(raw), hashlib.sha256(raw).hexdigest()
    if not isinstance(d, Mapping) or not isinstance(d.get("boards"), Mapping):
        raise ValueError("not a count decision: pass the count_decision.json written by "
                         "experiments/density/summarise.py --decide")
    if d.get("rule") != COUNT_RULE or d.get("approved") != RULES_APPROVED:
        raise ValueError(f"the decision was made under another rule text (approved {d.get('approved')!r}), not "
                         f"approved rule 7 of {RULES_APPROVED} as this code applies it")
    b = d["boards"].get(kind)
    if not isinstance(b, Mapping) or "override" not in b:
        raise ValueError(f"the decision file holds no decision for the {kind} board")
    if b.get("validation_target") != COUNT_RULE_TARGETS[kind]:
        raise ValueError(f"the {kind} decision comes from {b.get('validation_target')!r}, not "
                         f"{COUNT_RULE_TARGETS[kind]}: each test board is decided by its own validation target")
    override = b["override"]
    if override is not None and override not in COUNT_RULE_LEVELS:
        raise ValueError(f"the {kind} decision's override {override!r} is neither None (keep max_cells) nor one of the "
                         f"rule's lower levels {COUNT_RULE_LEVELS}")
    record = {"decision_file": None if path is None else str(path), "decision_sha256": sha, "kind": kind,
              "validation_target": b["validation_target"], "override": override,
              "ranking": [] if override is None else _passing_ranking(b, override, kind), "decision": b.get("decision"),
              "decision_git_commit": d.get("git_commit"),
              "decision_allowed_uncommitted": d.get("allow_uncommitted")}
    return override, record


def _finite_number(v) -> float | None:
    """v as a float when it is a finite real number (Python or numpy, not a bool), else None."""
    if isinstance(v, bool | np.bool_) or not isinstance(v, int | float | np.integer | np.floating):
        return None
    try:
        v = float(v)
    except OverflowError:  # an integer too large for a float
        return None
    return v if np.isfinite(v) else None


def _passing_ranking(b: Mapping, override: str, kind: str) -> list[str]:
    """One board's passing levels in approved rule 7's order (largest mean gain first, an exact tie going to the higher
    count), read from the count_override record the decision holds; the decided override must come first. Every lower
    level's verdict is checked against its own recorded numbers, both ways: it passes exactly when its mean gain exceeds
    COUNT_RULE_MIN_GAIN with at least COUNT_RULE_MIN_AHEAD of the release seeds ahead. Refused: no such record or an
    incomplete one, a level without a finite mean gain, a whole seed count and a verdict, a verdict its numbers
    contradict, a passing list that is not the passing levels in the rule's level order, and an override that is not
    the first."""
    rec = b.get("record")
    per_level = rec.get("per_level") if isinstance(rec, Mapping) else None
    if not isinstance(per_level, Mapping) or not isinstance(rec.get("passing"), list):
        raise ValueError(f"the {kind} decision holds no record of its passing levels (decide mode writes one)")
    complete = rec.get("complete")
    if not (isinstance(complete, bool | np.bool_) and complete):
        raise ValueError(f"the {kind} decision names a level although its record says the rerun is incomplete")
    n_seeds, gains = len(COUNT_RULE_SEEDS), {}
    for lv in COUNT_RULE_LEVELS:
        entry = per_level[lv] if isinstance(per_level.get(lv), Mapping) else {}
        g, ahead, verdict = _finite_number(entry.get("mean_gain")), entry.get("seeds_ahead"), entry.get("passes")
        whole = not isinstance(ahead, bool | np.bool_) and isinstance(ahead, int | np.integer) and 0 <= ahead <= n_seeds
        if g is None or not whole or not isinstance(verdict, bool | np.bool_):
            raise ValueError(f"the {kind} decision's record has no mean gain, seed count (0-{n_seeds}) and verdict for "
                             f"{lv}: {entry.get('mean_gain')!r}, {ahead!r}, {verdict!r}")
        if bool(verdict) != (g > COUNT_RULE_MIN_GAIN and ahead >= COUNT_RULE_MIN_AHEAD):
            raise ValueError(f"the {kind} decision's verdict for {lv} (passes {verdict!r}) disagrees with its numbers "
                             f"(mean gain {g:+.3f}, {ahead} of {n_seeds} seeds ahead): rule 7 passes a mean gain above "
                             f"{COUNT_RULE_MIN_GAIN:g} with at least {COUNT_RULE_MIN_AHEAD} seeds ahead")
        if verdict:
            gains[lv] = g
    passing = [lv for lv in COUNT_RULE_LEVELS if lv in gains]
    if rec["passing"] != passing:
        raise ValueError(f"the {kind} decision lists {rec['passing']!r} as passing, but its per-level verdicts pass "
                         f"{passing}")
    ranking = sorted(passing, key=lambda lv: (gains[lv], float(_MAX_LABEL.match(lv)["f"])), reverse=True)
    if not ranking or ranking[0] != override:
        raise ValueError(f"the {kind} decision's override {override!r} is not the first of its passing levels "
                         f"{ranking} in approved rule 7's order (largest mean gain, ties to the higher count)")
    return ranking


def n_cells_from_decision(board: Board, pool_size: int, decision_json, kind: str | None = None) -> tuple[int, dict]:
    """`n_cells` with the level that approved rule 7's decision file gives this T2 test board (`read_count_decision`),
    so the count is applied by code, never copied by hand. The count record gains `count_decision` (the decision file's
    path and sha256, its validation target, override, passing levels in rule order and decision).

    A decided level that gives the test board fewer cells than its lower limit (truth_cells where published and at
    least min_cells, else min_cells) gives way under COUNT_FLOOR_RULE to the next passing count up: the passing level
    with the fewest cells that the board allows, whatever its gain, else max_cells. The user approved "the next passing
    count the board allows, or max_cells if none does" on 2026-10-21, beside rule 7, and left the reading of "next" to
    the AI assistant ("anything that would improve my chances of winning"); a decision panel found the next count up
    the better bet (experiments/density/floor_reading/README.md). The record then gains `count_floor` (the passing
    levels below the limit with their counts, the level used and why). Under the documented limits only 0.25max can
    fall below. No level is clipped to the limit; other refusals of n_cells (no max_cells, inconsistent limits) pass
    unchanged."""
    override, decision = read_count_decision(decision_json, board, kind)
    level, below = override, []
    if override is not None:
        max_cells, _, _, floor, floor_source = _limits(board)
        counts = {lv: _override_count(lv, max_cells) for lv in decision["ranking"]}
        if counts[override] < floor:  # the next count up: the allowed passing level with the fewest cells
            allowed = sorted((lv for lv in counts if counts[lv] >= floor),
                             key=lambda lv: (counts[lv], float(_MAX_LABEL.match(lv)["f"])))
            level = allowed[0] if allowed else None
            below = [{"level": lv, "n": counts[lv]} for lv in decision["ranking"] if counts[lv] < floor]
    n, record = n_cells(board, pool_size, override=level)
    record["count_decision"] = decision
    if below:
        used = (f"use {level}, the passing level with the fewest cells the board allows" if level is not None
                else "no passing level is allowed, so keep max_cells")
        record["count_floor"] = {
            "rule": COUNT_FLOOR_RULE, "floor": floor, "floor_source": floor_source, "below_floor": below,
            "level": level, "decision": f"{override} gives {below[0]['n']} cells, below {floor_source} {floor}: {used}"}
    return n, record


# ---------------------------------------------------------------------------------------------------------------------
# Sampling
_balanced_rows = draw_balanced  # its name here before it moved to vec.core (feat-020)


def _ball_shifts(m: int, r_max: float, rng: np.random.Generator) -> np.ndarray:
    """m isotropic displacements: uniform directions, lengths uniform in [JITTER_MIN_FRAC, 1] x r_max."""
    d = rng.normal(size=(m, 3))
    norm = np.linalg.norm(d, axis=1)
    while (bad := norm < DEGENERATE_RMS).any():  # a zero-length draw has probability zero; redraw it anyway
        d[bad] = rng.normal(size=(int(bad.sum()), 3))
        norm = np.linalg.norm(d, axis=1)
    return d / norm[:, None] * (r_max * rng.uniform(JITTER_MIN_FRAC, 1.0, size=m))[:, None]


@dataclass
class CellSample:
    """Cells drawn from a pool: `rows[i]` is the pool row of submitted cell i (in submission order) and `coords` its
    position after duplicate breaking. Take expression as `pool.X[sample.rows]`; `record` says what was done."""

    rows: np.ndarray
    coords: np.ndarray
    record: dict


def sample_cells(stage, n: int, seed: int, jitter_frac: float = JITTER_NN_FRAC, herd_expr=None) -> CellSample:
    """Draw n cells from a pool at random rows, then break exact duplicate coordinates.

    stage: an AnnData with obsm["spatial_3D"] (e.g. `vec.prep.load_stage(..., drop_duplicates=True)`), a
    vec.evaluate.Stage, or an (N, 3) array. Rows are random, never contiguous, and the output order is random too, so
    any later prefix of the file is still a random subset. Rows repeat only when n > N, each pool row then appearing
    n // N or n // N + 1 times. Every cell whose coordinates repeat an earlier cell's exactly (a repeated row, or a
    duplicate position in the pool) moves by an isotropic displacement of length 0.5-1 x jitter_frac x the median
    nearest-neighbour distance among the sample's distinct positions; the first occurrence stays exact. jitter_frac
    may not exceed 0.25, the budget measured to be free (experiments/density). The jitter is 3-D, so a moved cell can
    sit off its section's z level by at most that length.

    Draws use the streams [seed, 0] (rows) and [seed, 1] (jitter), so a jitter change never moves the rows.

    herd_expr (cells x genes of the pool, e.g. its log1p CP10k panel; T2 try H4 of experiments/harness/TRIES.md): the
    random rows are then herded towards the pool's mean expression (`vec.core.herd_rows`, stream [seed, 2]) before
    positions are read, so each cell keeps its own position; without it the draw is unchanged.
    """
    C = _coords_of(stage)
    N = C.shape[0]
    n = _whole(n, "n", 1)
    seed = _whole(seed, "seed", 0)
    if not 0 < jitter_frac <= JITTER_NN_FRAC:
        raise ValueError(f"jitter_frac must be in (0, {JITTER_NN_FRAC}] (the free budget), got {jitter_frac}")
    rows = draw_balanced(N, n, np.random.default_rng([seed, 0]))
    herd_rec = None
    if herd_expr is not None:
        rows, herd_rec = herd_rows(herd_expr, np.arange(N), rows, np.random.default_rng([seed, 2]))
    out = C[rows].copy()
    # first occurrence of each distinct position as the scorer reads it (float32: float-noise copies count as repeats)
    _, first = np.unique(out.astype(np.float32), axis=0, return_index=True)
    dup = np.ones(n, dtype=bool)
    dup[first] = False
    repeat_row = np.ones(n, dtype=bool)
    repeat_row[np.unique(rows, return_index=True)[1]] = False
    nn = median_nn_distance(out)  # over distinct positions
    r_max = jitter_frac * nn if np.isfinite(nn) else float("nan")
    shift = np.zeros(n)
    if dup.any():
        if not r_max > 0:
            raise ValueError("cannot break duplicate coordinates: fewer than two distinct positions in the sample")
        d = _ball_shifts(int(dup.sum()), r_max, np.random.default_rng([seed, 1]))
        out[dup] += d
        shift[dup] = np.linalg.norm(d, axis=1)
    n_left = n - np.unique(out.astype(np.float32), axis=0).shape[0]
    if n_left:
        raise RuntimeError(f"{n_left} duplicate coordinates remain after jitter (float32); the cloud is degenerate")
    k = n // N
    record = {
        "n": n, "pool_size": N, "seed": seed, "order": "random", "with_replacement": n > N,
        "copies_per_row": [k, k + 1] if n % N else [k, k],
        "n_distinct_rows": int(np.unique(rows).size), "n_repeated_rows": int(repeat_row.sum()),
        "n_duplicate_coords": int(dup.sum()), "dup_from_repeated_rows": int((dup & repeat_row).sum()),
        "dup_from_pool_positions": int((dup & ~repeat_row).sum()), "n_jittered": int(dup.sum()),
        "median_nn": nn, "jitter_frac": float(jitter_frac), "jitter_radius_max": r_max,
        "jitter_radius_min": JITTER_MIN_FRAC * r_max,
        "max_shift": float(shift.max()), "mean_shift_jittered": float(shift[dup].mean()) if dup.any() else 0.0,
    }
    if herd_rec is not None:
        record["herd"] = herd_rec
    return CellSample(rows=rows, coords=out, record=record)


# ---------------------------------------------------------------------------------------------------------------------
# Coordinates
def centre(coords) -> np.ndarray:
    """The cloud minus its centroid (a float64 copy)."""
    C = _as_coords(coords)
    return C - C.mean(0)


def check_rotation(R, tol: float = ROTATION_TOL) -> np.ndarray:
    """R as a float64 3x3 proper rotation, or ValueError. R must be orthonormal (R^T R = I within tol) with det +1.
    An improper matrix (det -1, a reflection possibly combined with a rotation) would mirror the cloud and is refused,
    never repaired: the released frames share one handedness and a mirrored copy of a correct shape loses up to 0.12
    Dice under the pinned scorer (experiments/coords/REPORT.md)."""
    R = np.asarray(R, dtype=np.float64)
    if R.shape != (3, 3) or not np.isfinite(R).all():
        raise ValueError(f"a rotation must be a finite 3x3 matrix, got shape {R.shape}")
    if not np.allclose(R.T @ R, np.eye(3), atol=tol):
        raise ValueError("matrix is not orthonormal: a scaling or shear is not a rotation (use rescale_to_rms for "
                         "size)")
    det = float(np.linalg.det(R))
    if det < 0:
        raise ValueError(f"improper rotation (det {det:+.6f}) would mirror the cloud: refused (proper rotations only)")
    return R


def proper_rotation(coords, R, about=None) -> np.ndarray:
    """Rotate a cloud by the proper rotation R (column-vector convention, x -> R x) about `about` (default: the
    centroid). Improper or non-orthonormal matrices are refused (`check_rotation`)."""
    C = _as_coords(coords)
    R = check_rotation(R)
    p = C.mean(0) if about is None else np.asarray(about, dtype=np.float64).reshape(3)
    return p + (C - p) @ R.T


def rescale_to_rms(coords, rms: float) -> np.ndarray:
    """Scale a cloud about its centroid so that its 3-D RMS radius equals `rms`: isotropic and proper (no rotation, no
    mirror). The same operation as experiments/coords/size_procedure.rescale_to_rms."""
    C = _as_coords(coords)
    target = float(rms)
    if not np.isfinite(target) or target <= 0:
        raise ValueError(f"target RMS must be finite and positive, got {rms!r}")
    current = rms_radius(C)
    if current < DEGENERATE_RMS:
        raise ValueError("the cloud has no spread: it cannot be rescaled")
    mu = C.mean(0)
    return mu + (C - mu) * (target / current)


def xy_scope(coords, core) -> float:
    """The size procedure's in-plane scope: log(in-plane RMS of the whole cloud / in-plane RMS of its core cells)."""
    C = _as_coords(coords)[:, :2]
    core = np.asarray(core, dtype=bool)

    def xy_rms(X: np.ndarray) -> float:
        return float(np.sqrt(((X - X.mean(0)) ** 2).sum(1).mean()))

    return float(np.log(xy_rms(C) / xy_rms(C[core])))


def radial_trim_mask(coords, core, frac: float) -> np.ndarray:
    """The heart s2's scope trim (docs/S2_RULES.md): cells kept when surrounding (non-core) cells are removed farthest
    in plane from the core centroid first, until the in-plane scope is at most `frac` of the cloud's own (the fewest
    removed that reach it, by bisection on the removal order). Core cells are always kept. The same rule as the s2
    hedge payoffs' radial trim (experiments/s2_hedge/payoff.py trim_mask, mode "radial"); tests pin the two together."""
    C = _as_coords(coords)
    core = np.asarray(core, dtype=bool)
    if core.shape != (len(C),) or not core.any():
        raise ValueError("core must be a boolean mask over the cells with at least one core cell")
    if not 0.0 <= float(frac) <= 1.0:
        raise ValueError(f"frac must lie in [0, 1], got {frac!r}")
    target = float(frac) * xy_scope(C, core)
    shell = np.flatnonzero(~core)
    d = np.linalg.norm(C[shell, :2] - C[core, :2].mean(0), axis=1)
    order = shell[np.argsort(-d, kind="stable")]

    def keep_after(k: int) -> np.ndarray:
        keep = np.ones(len(C), dtype=bool)
        keep[order[:k]] = False
        return keep

    lo, hi = 0, len(order)
    while lo < hi:
        mid = (lo + hi) // 2
        m = keep_after(mid)
        if xy_scope(C[m], core[m]) <= target:
            hi = mid
        else:
            lo = mid + 1
    return keep_after(lo)


def _scale_xy_hold_z(C: np.ndarray, target_full: float | None, target_xy: float | None) -> tuple[np.ndarray, dict]:
    """Scale x and y about the centroid, z untouched, so that the 3-D RMS radius equals target_full, or (when
    target_full is None) the in-plane RMS radius equals target_xy."""
    mu = C.mean(0)
    X = C - mu
    rms_xy = float(np.sqrt((X[:, :2] ** 2).sum(1).mean()))
    sd_z = float(np.sqrt((X[:, 2] ** 2).mean()))
    if rms_xy < DEGENERATE_RMS:
        raise ValueError("the cloud has no in-plane spread: it cannot be rescaled in x and y")
    if target_full is not None:
        if target_full <= sd_z:
            raise ValueError(f"target RMS {target_full:.6g} is not above the cloud's held z SD {sd_z:.6g}: z cannot be "
                             "held")
        target_xy = float(np.sqrt(target_full ** 2 - sd_z ** 2))
    s = target_xy / rms_xy
    out = C.copy()
    out[:, :2] = mu[:2] + X[:, :2] * s
    return out, {"scale_xy": s, "rms_xy_before": rms_xy, "rms_xy_after": rms_xy * s, "sd_z_cloud": sd_z}


def _forecast(v, name: str) -> float:
    """A positive finite number from the size JSON, or ValueError (dry-run JSONs carry a redaction marker instead)."""
    if isinstance(v, str) and v.startswith(REDACTED_PREFIX):
        raise ValueError(f"{name} is redacted: a dry-run JSON (a code check, not a prediction) cannot size a "
                         "submission")
    if isinstance(v, bool) or not isinstance(v, int | float) or not np.isfinite(v) or v <= 0:
        raise ValueError(f"{name} must be a positive number, got {v!r}")
    return float(v)


def load_size_result(size_json) -> tuple[dict, str | None]:
    """The size procedure's output (a path to its JSON, or the dict) and the sha256 of the file (None for a dict)."""
    if isinstance(size_json, dict):
        return size_json, None
    raw = Path(size_json).read_bytes()
    return json.loads(raw), hashlib.sha256(raw).hexdigest()


def apply_size_procedure(coords, size_json, hypothesis: str | None = None, z_mode: str | None = None,
                         match: str = "full", expect_board: str | None = None) -> tuple[np.ndarray, dict]:
    """Rescale a cloud to the size the pre-registered size procedure predicts, in the z mode it recorded.

    size_json: the JSON written by experiments/coords/size_procedure.py (`run`, CLI `--out`), or its dict. Refused:
    a result with valid=false (for the heart, size_or_hold applies approved rule 3's fallback, PLAN 4.3 item 5),
    redacted dry-run forecasts, unknown boards, modes or hypotheses. hypothesis: the heart scope hypothesis ("H1"
    trend for s1, "H2" half for s2, "H0" hold); required for the heart, refused for the embryo (which has none).
    z_mode: None applies the pre-registered mode recorded in the JSON; a logged alternative may be named for
    diagnosis.

    - "raw" / "common_step": isotropic `rescale_to_rms` about the centroid to the predicted full-cloud RMS (in
      common_step the target is in the reference stage's z step, so the cloud must be in that convention).
    - "xy_grow_z_hold": z stays native and untouched; x and y scale about the centroid so that the 3-D RMS equals the
      predicted full-cloud RMS (match="full", what scale_log_ratio scores), or so that the in-plane RMS equals the
      procedure's in-plane forecast (match="xy"). The two agree when the cloud's z SD equals the z SD the procedure
      held (the reference stage's for the heart, the lower bracketing stage's for the embryo); the record logs both.
    - "split" is refused: the procedure logs it for diagnosis and does not record the full cloud's xy/z split.

    Returns the new coordinates and a record (sizes, factors, recorded and used z mode, the JSON's provenance and
    sha256, warnings). The size numbers come only from the JSON.
    """
    size, sha = load_size_result(size_json)
    missing = [k for k in ("board", "valid", "preregistered") if k not in size]
    if missing:
        raise ValueError(f"not a size procedure output (missing {missing}): pass the JSON written by "
                         "experiments/coords/size_procedure.py")
    board = size["board"]
    if board not in SIZE_BOARDS:
        raise ValueError(f"unknown size board {board!r}")
    if expect_board is not None and board != expect_board:
        raise ValueError(f"size JSON is for the {board} board, expected {expect_board}")
    if size["valid"] is not True:
        raise ValueError("the size result has valid=false: no applicable growth rule, so its numbers must not be used")
    if match not in SIZE_MATCH:
        raise ValueError(f"match must be one of {SIZE_MATCH}, got {match!r}")
    pre = size["preregistered"]
    recorded = pre.get("z_mode")
    mode = recorded if z_mode is None else z_mode
    if mode == "split":
        raise ValueError("z mode 'split' is logged for diagnosis only; the JSON does not record the full cloud's "
                         "in-plane/z split")
    if mode not in APPLICABLE_Z_MODES:
        raise ValueError(f"z mode {mode!r} cannot be applied; use one of {APPLICABLE_Z_MODES}")
    res = pre if mode == recorded else size.get("z_modes", {}).get(mode)
    if res is None:
        raise ValueError(f"the size JSON holds no {mode!r} result")
    if res.get("valid") is not True:
        raise ValueError(f"the {mode!r} result has valid=false")
    if match == "xy" and mode != "xy_grow_z_hold":
        raise ValueError("match='xy' applies only to the xy_grow_z_hold mode")

    stages = {g.get("name"): g for g in res.get("stages", [])}
    if board == "heart":
        if hypothesis not in SIZE_HYPOTHESES:
            raise ValueError(f"heart: pass the scope hypothesis, one of {SIZE_HYPOTHESES} (s1 = H1, s2 = H2), got "
                             f"{hypothesis!r}")
        scope = res.get("scope", {}).get(hypothesis)
        if scope is None:
            raise ValueError(f"the size JSON holds no scope hypothesis {hypothesis!r}")
        target_full = _forecast(scope.get("pred_rms_full"), f"scope.{hypothesis}.pred_rms_full")
        target_xy = (_forecast(scope.get("pred_rms_xy_full"), f"scope.{hypothesis}.pred_rms_xy_full")
                     if mode == "xy_grow_z_hold" else None)
        reference = res.get("reference")
        z_reference = reference
    else:
        if hypothesis is not None:
            raise ValueError("embryo: the embryo procedure has no scope hypotheses; pass hypothesis=None")
        target_full = _forecast(res.get("pred_rms_full"), "pred_rms_full")
        target_xy = None
        bracket = res.get("bracket") or [None]
        reference, z_reference = bracket, bracket[0]

    C = _as_coords(coords)
    before = rms_radius(C)
    rec: dict = {
        "board": board, "target_t": size.get("target_t"), "size_timestamp": size.get("timestamp"),
        "size_json_sha256": sha, "z_mode_recorded": recorded, "z_mode_used": mode, "hypothesis": hypothesis,
        "match": match, "target_rms_full": target_full, "reference": reference, "rms_before": before,
        "warnings": [],
    }
    prov = size.get("provenance")
    if prov is None:
        rec["warnings"].append("the size JSON carries no provenance (git commit, input checksums)")
    else:
        rec["git_commit"] = prov.get("git_commit")
        rec["size_procedure_sha256"] = prov.get("size_procedure_sha256")
        if prov.get("git_uncommitted"):
            rec["warnings"].append(f"the size JSON was produced with uncommitted procedure files: "
                                   f"{prov['git_uncommitted']}")
    if board == "heart":
        rec["scale_vs_reference_full"] = scope.get("scale_vs_reference_full")
    ref_geom = stages.get(reference) if isinstance(reference, str) else None
    if ref_geom is not None and isinstance(ref_geom.get("rms_full"), int | float):
        rec["log_rms_before_vs_reference"] = float(np.log(before / ref_geom["rms_full"]))

    if mode in ("raw", "common_step"):
        out = rescale_to_rms(C, target_full)
        rec.update(z_held=False, scale=target_full / before)
        if mode == "common_step":
            rec["warnings"].append("common_step: the target is in the last stage's z step (one-pitch hypothesis); the "
                                   "cloud must be in that convention")
    else:
        zg = stages.get(z_reference)
        if zg is None or not isinstance(zg.get("sd_z_full"), int | float):
            raise ValueError(f"the size JSON holds no z SD for the stage {z_reference!r} whose z the procedure held")
        sd_z_ref = float(zg["sd_z_full"])
        if target_xy is None:  # embryo: in-plane forecast implied by the held z SD
            if target_full <= sd_z_ref:
                raise ValueError("inconsistent size JSON: the full-cloud forecast is not above the held z SD")
            target_xy = float(np.sqrt(target_full ** 2 - sd_z_ref ** 2))
        out, info = _scale_xy_hold_z(C, target_full if match == "full" else None, target_xy)
        rec.update(z_held=True, target_rms_xy_full=target_xy, sd_z_held_by_procedure=sd_z_ref, **info,
                   log_rms_xy_after_vs_forecast=float(np.log(info["rms_xy_after"] / target_xy)),
                   log_sd_z_cloud_vs_held=float(np.log(info["sd_z_cloud"] / sd_z_ref)))
        if abs(rec["log_sd_z_cloud_vs_held"]) > np.log(NATIVE_Z_WARN_FACTOR):
            rec["warnings"].append(f"the cloud's z SD ({info['sd_z_cloud']:.4g}) differs from the held z SD "
                                   f"({sd_z_ref:.4g}) by more than {NATIVE_Z_WARN_FACTOR:g}x: is z native to "
                                   f"{z_reference}?")
    rec["rms_after"] = rms_radius(out)
    rec["log_rms_after_vs_target"] = float(np.log(rec["rms_after"] / target_full))
    return out, rec


def _positive(v) -> bool:
    """A finite positive number (not a bool)."""
    return not isinstance(v, bool) and isinstance(v, int | float) and bool(np.isfinite(v)) and v > 0


def _redacted(obj) -> bool:
    """Whether a size JSON carries a dry-run redaction marker anywhere."""
    if isinstance(obj, str):
        return obj.startswith(REDACTED_PREFIX)
    if isinstance(obj, dict):
        return any(_redacted(v) for v in obj.values())
    if isinstance(obj, list | tuple):
        return any(_redacted(v) for v in obj)
    return False


def _header_error(size: dict, expect_board: str | None) -> str | None:
    """The first refusal apply_size_procedure makes of a size JSON's header, in its words, or None."""
    missing = [k for k in ("board", "valid", "preregistered") if k not in size]
    if missing:
        return (f"not a size procedure output (missing {missing}): pass the JSON written by "
                "experiments/coords/size_procedure.py")
    if size["board"] not in SIZE_BOARDS:
        return f"unknown size board {size['board']!r}"
    if expect_board is not None and size["board"] != expect_board:
        return f"size JSON is for the {size['board']} board, expected {expect_board}"
    return None


def _stage_times(size: dict) -> dict:
    """The stages of a size JSON's pre-registered result, as {name: nominal E-day}."""
    pre = size.get("preregistered") if isinstance(size.get("preregistered"), dict) else {}
    return {str(g.get("name")): g.get("t") for g in pre.get("stages") or [] if isinstance(g, dict)}


def _invalid_reason(size: dict) -> dict:
    """Why a size JSON says valid=false, in its own fields: any reason or note fields at its top level or in its
    pre-registered result, the growth series it marks invalid with their applicable rules (none: no rate interval
    starts at or after min_rate_start_t), and the stages it had (name: nominal E-day). No forecast is copied."""
    pre = size.get("preregistered") if isinstance(size.get("preregistered"), dict) else {}
    why = {f"{where}.{k}": src[k] for where, src in (("top", size), ("preregistered", pre)) for k in REASON_KEYS
           if k in src}
    growth = pre.get("growth") if isinstance(pre.get("growth"), dict) else {}
    why["invalid_growth"] = {tag: {"applicable": g.get("applicable")} for tag, g in growth.items()
                             if isinstance(g, dict) and g.get("valid") is not True}
    why["min_rate_start_t"] = pre.get("min_start_t")
    why["stages"] = _stage_times(size)
    return why


def _source_measurement(size: dict, ref: np.ndarray) -> tuple[dict, float]:
    """The heart size JSON's own measurement of E10.5 (a stage of its pre-registered result at t = 10.5) that the
    caller's reference cloud matches, and the reference cloud's RMS radius. The match is the RMS radius within float
    noise (SIZE_CROSSCHECK_LOG_TOL): the cloud the procedure measured gives exactly its RMS radius.

    Raised, with no size applied: the procedure did not measure E10.5, or the reference cloud is not the cloud it
    measured. For a valid JSON that leaves the heart cells' stage unshown (approved rule 3 sizes heart cells from E10.5
    only), and for an invalid one the reference's size is not E10.5's."""
    valid = size.get("valid") is True
    pre = size["preregistered"] if isinstance(size["preregistered"], dict) else {}
    t_source = float(SIZE_FALLBACK_STAGE.removeprefix("E"))
    measured = [g for g in pre.get("stages") or [] if isinstance(g, dict) and g.get("t") == t_source
                and _positive(g.get("rms_full"))]
    if not measured and valid:
        raise ValueError(f"the heart size result is valid but the procedure did not measure {SIZE_FALLBACK_STAGE} "
                         f"(its stages: {_stage_times(size)}): approved rule 3 ({RULES_APPROVED}) sizes heart cells "
                         f"from {SIZE_FALLBACK_STAGE} only, and under the pre-registered settings every valid "
                         f"release-day heart JSON measured {SIZE_FALLBACK_STAGE} (its last two stages are E9.5 and "
                         f"{SIZE_FALLBACK_STAGE}), so this is not the release-day JSON. If {SIZE_FALLBACK_STAGE} is "
                         "barred (assumption A20), the size rule comes back to the captain; no size is applied")
    if not measured:
        raise ValueError(f"the heart size result has valid=false and the procedure did not measure "
                         f"{SIZE_FALLBACK_STAGE} (its stages: {_stage_times(size)}): approved rule 3 "
                         f"keeps {SIZE_FALLBACK_STAGE}'s size only while {SIZE_FALLBACK_STAGE} is usable, and the "
                         "release-day run (experiments/coords/run_p3.sh) measures it then. A JSON without it has the "
                         "layout of E10.5 barred (assumption A20), whose size rule comes back to the captain, or is "
                         "not the release-day JSON; no size is applied")
    target = rms_radius(ref)
    same = [g for g in measured if target > 0 and abs(np.log(target / g["rms_full"])) <= SIZE_CROSSCHECK_LOG_TOL]
    if not same:
        seen = "; ".join(f"{g.get('name')}: {g.get('n_cells')} cells, RMS radius {g['rms_full']:.6g}"
                         for g in measured)
        full = (f"the full {SIZE_FALLBACK_STAGE} cloud (every cell, duplicated records dropped, no scope filter: "
                "vec.prep.load_stage(..., drop_duplicates=True))")
        if valid:
            why = (f"so the heart cells are not shown to come from {SIZE_FALLBACK_STAGE}: approved rule 3 "
                   f"({RULES_APPROVED}) sizes heart cells from {SIZE_FALLBACK_STAGE} only, and for cells from another "
                   "stage (E10.5 barred, assumption A20: the heart plan is rebuilt from E9.5) the size rule comes back "
                   f"to the captain. If the cells are {SIZE_FALLBACK_STAGE}'s, pass {full}, as the fallback for an "
                   "invalid JSON needs")
        else:
            why = f"so its size is not {SIZE_FALLBACK_STAGE}'s: pass {full}"
        raise ValueError(f"the reference cloud ({ref.shape[0]} cells, RMS radius {target:.6g}) is not the cloud the "
                         f"size procedure measured as {SIZE_FALLBACK_STAGE} ({seen}), {why}; no size is applied")
    return same[0], target


def size_or_hold(coords, size_json, hypothesis: str | None = None, *, source_stage: str | None = None,
                 reference_cloud=None, z_mode: str | None = None, match: str = "full",
                 expect_board: str | None = None) -> tuple[np.ndarray, dict]:
    """The T2 heart pipeline's size step under approved rule 3 (docs/RELEASE_RULES_SUMMARY.md, approved answer 3,
    2026-10-21): `apply_size_procedure`, or the rule's fallback for an invalid heart size JSON.

    Heart cells must come from E10.5. For cells from any other stage no size is applied, whatever the JSON says: if
    E10.5 is barred (assumption A20 fails), the heart plan is rebuilt from E9.5 and the size rule comes back to the
    captain. (Under the pre-registered settings every valid heart JSON measured E10.5: a valid result needs its last
    two stages at or after E9.0, which among the released heart stages are E9.5 and E10.5.) The caller states the
    stage, and for every heart JSON, valid or not, its reference cloud must be the cloud the JSON measured as E10.5
    (its stage at t = 10.5, RMS radius within float noise): a stale source_stage="E10.5" passed with another stage's
    cloud (the heart plan rebuilt from E9.5) gets no size.

    For E10.5 cells, a heart JSON with valid=true returns exactly what apply_size_procedure returns for the same
    arguments once the reference cloud matches (checked after apply_size_procedure's own refusals), and so does an
    embryo JSON whose valid flag is not false. A heart JSON with valid=false keeps E10.5's size: the cloud is rescaled
    isotropically about its centroid (`rescale_to_rms`) to the RMS radius of `reference_cloud`, the full E10.5 cloud,
    so the scale term sits at the floor instead of the build failing. The scope hypothesis does not change that size:
    s1 (H1) and s2 (H2) both keep E10.5's size and differ by scope, i.e. by which cells the pipeline includes; the
    hypothesis is recorded. The fallback needs the procedure's own measurement of E10.5 (the release-day run,
    experiments/coords/run_p3.sh, measures E10.5 whenever it is usable): a JSON without one has the layout of E10.5
    barred, or is not the release-day JSON, and is refused; so is a reference cloud whose RMS radius differs from that
    measurement beyond float noise, which is not E10.5's size. No forecast of the invalid JSON is used or copied: the
    size is the reference cloud's, checked against the JSON's measurement of E10.5; the record copies that measurement
    and the JSON's descriptive fields (target time, timestamp, recorded z mode, why it is invalid, its stage names and
    times).

    source_stage: the stage the prediction's cells were drawn from, as the caller states it (never inferred): "E10.5".
    reference_cloud: that stage's full cloud as the size procedure measured it (every cell, duplicated records dropped,
    no scope filter: the pool from `vec.prep.load_stage(..., drop_duplicates=True)`), as an AnnData, a
    vec.evaluate.Stage or an (n, 3) array. A heart JSON needs both, valid or not, and the reference is checked against
    the JSON's measurement of E10.5 either way. The rehearsal JSON (run_p3.sh --rehearsal) has redacted forecasts and
    is refused as redacted before that check, so a rehearsal catches a missing argument and a wrong stage name,
    hypothesis or board, but not a wrong reference cloud: a reference that is scope-filtered, rescaled or another
    stage's first meets the check on release day, where it is refused, not applied (the tests run the check on
    unredacted JSONs). The embryo ignores both.

    Raised, with no size applied: heart cells from a stage other than E10.5 (approved rule 3 returns the size rule to
    the captain); a heart JSON, valid or not, that holds no measurement of E10.5, or whose measurement the reference
    cloud does not match; an embryo JSON with valid=false (no fallback was approved for the embryo); and the refusals
    of apply_size_procedure (it decides every JSON whose valid flag is not false, before the reference check: not the
    procedure's output, an unknown board, mode or hypothesis, a board other than expect_board, redacted dry-run
    forecasts). For valid=false the same checks run here, and a redaction marker anywhere refuses the JSON. The
    fallback rescales isotropically, so it takes no z_mode and only match="full".

    Returns the new coordinates and a record. The fallback's record: `fallback` True, the rule and its approval date,
    the JSON's own account of why it is invalid (`fallback_reason`), the hypothesis, the sizes (target_rms_full, the
    reference cloud's RMS radius; rms_before; scale; rms_after), the JSON's measurement of E10.5 that the reference
    matched, the JSON's sha256 and provenance, and warnings.
    """
    size, sha = load_size_result(size_json)
    header = _header_error(size, expect_board) if isinstance(size, dict) else "not a size procedure output"
    ref = None
    if header is None and size["board"] == "heart":
        if not isinstance(source_stage, str) or not source_stage:
            raise ValueError(f"heart: state the stage the cells were drawn from (source_stage, e.g. "
                             f"{SIZE_FALLBACK_STAGE!r}); approved rule 3 needs it if the size JSON is invalid")
        if reference_cloud is None:
            raise ValueError("heart: pass the full cloud of the source stage (reference_cloud); approved rule 3 keeps "
                             "its size if the size JSON is invalid")
        ref = _coords_of(reference_cloud)
        if source_stage != SIZE_FALLBACK_STAGE:
            hint = (f" If these cells are {SIZE_FALLBACK_STAGE}'s, pass source_stage={SIZE_FALLBACK_STAGE!r} (the "
                    "stage, not a file name)." if SIZE_FALLBACK_STAGE in source_stage else "")
            raise ValueError(f"the heart cells come from {source_stage!r}, not {SIZE_FALLBACK_STAGE}: approved rule 3 "
                             f"({RULES_APPROVED}) sizes heart cells from {SIZE_FALLBACK_STAGE} only. For cells from "
                             "another stage (E10.5 barred, assumption A20: the heart plan is rebuilt from E9.5) the "
                             "size rule comes back to the captain, valid size JSON or not (a valid one measured "
                             f"E10.5). No size is applied.{hint}")
    if not (isinstance(size, dict) and size.get("valid") is False):
        out, rec = apply_size_procedure(coords, size_json, hypothesis, z_mode=z_mode, match=match,
                                        expect_board=expect_board)
        if ref is not None:  # a valid heart JSON, past apply_size_procedure's own refusals
            _source_measurement(size, ref)  # the cells' stage: the reference must be the cloud measured as E10.5
        return out, rec

    if header is not None:
        raise ValueError(header)
    board = size["board"]
    if board == "embryo":
        raise ValueError("the embryo size result has valid=false: no fallback was approved for the embryo (approved "
                         "rule 3 covers the heart only), so its numbers must not be used")
    if _redacted(size):
        raise ValueError("the size JSON is redacted: a dry-run JSON (a code check, not a prediction) cannot size a "
                         "submission, not even through the fallback")
    if hypothesis not in SIZE_HYPOTHESES:
        raise ValueError(f"heart: pass the scope hypothesis, one of {SIZE_HYPOTHESES} (s1 = H1, s2 = H2), got "
                         f"{hypothesis!r}")
    if z_mode is not None or match != "full":
        raise ValueError(f"the fallback keeps {SIZE_FALLBACK_STAGE}'s size by an isotropic rescale: it takes no z_mode "
                         f"and only match='full' (got z_mode={z_mode!r}, match={match!r})")
    pre = size["preregistered"] if isinstance(size["preregistered"], dict) else {}
    source, target = _source_measurement(size, ref)

    C = _as_coords(coords)
    before = rms_radius(C)
    out = rescale_to_rms(C, target)
    rec: dict = {
        "board": board, "target_t": size.get("target_t"), "size_timestamp": size.get("timestamp"),
        "size_json_sha256": sha, "fallback": True, "fallback_rule": SIZE_FALLBACK_RULE,
        "fallback_approved": RULES_APPROVED, "fallback_reason": _invalid_reason(size), "source_stage": source_stage,
        "hypothesis": hypothesis, "hypothesis_changes_size": False, "z_mode_recorded": pre.get("z_mode"),
        "z_mode_used": "raw", "z_held": False, "match": match, "target_rms_full": target,
        "target_source": f"RMS radius of the full {SIZE_FALLBACK_STAGE} cloud (reference_cloud, {ref.shape[0]} cells)",
        "size_json_stage_of_source": source.get("name"), "size_json_rms_full_of_source": float(source["rms_full"]),
        "rms_before": before, "scale": target / before, "warnings": [],
    }
    prov = size.get("provenance")
    if prov is None:
        rec["warnings"].append("the size JSON carries no provenance (git commit, input checksums)")
    else:
        rec["git_commit"] = prov.get("git_commit")
        rec["size_procedure_sha256"] = prov.get("size_procedure_sha256")
        if prov.get("git_uncommitted"):
            rec["warnings"].append(f"the size JSON was produced with uncommitted procedure files: "
                                   f"{prov['git_uncommitted']}")
    rec["rms_after"] = rms_radius(out)
    rec["log_rms_after_vs_target"] = float(np.log(rec["rms_after"] / target))
    return out, rec


# ---------------------------------------------------------------------------------------------------------------------
# Guard
def _right_handed_basis(C: np.ndarray) -> np.ndarray:
    """PCA axes of a cloud (rows, by decreasing variance) made right-handed (det +1), as the pinned scorer does."""
    X = C - C.mean(0)
    _, _, Vt = np.linalg.svd(X, full_matrices=False)
    if np.linalg.det(Vt) < 0:
        Vt[-1] *= -1
    return Vt


def _canonical(C: np.ndarray, basis: np.ndarray, scale: float | None = None) -> np.ndarray:
    """Centred on the cloud's own centroid, rotated by `basis`, divided by `scale` (default: the cloud's own RMS)."""
    X = C - C.mean(0)
    return X @ basis.T / (float(np.sqrt((X ** 2).sum(1).mean())) if scale is None else scale)


def _voxel_ids(Z: np.ndarray) -> np.ndarray:
    """Occupied voxels of canonical coordinates on veckit's occupancy_dice grid, as sorted unique integer ids."""
    v = np.floor((np.clip(Z, -GRID_EXTENT, GRID_EXTENT - VOXEL_EDGE_EPS) + GRID_EXTENT) / (2 * GRID_EXTENT)
                 * GRID_N).astype(np.int64)
    return np.unique((v[:, 0] * GRID_N + v[:, 1]) * GRID_N + v[:, 2])


@dataclass
class GuardReport:
    """`cloud_guard` result. Coverage = share of the source's occupied canonical voxels that the cloud occupies;
    outside = share of the cloud's occupied voxels that the source leaves empty. `ref_*` are means over same-count
    random draws from the source. `flags` is empty when the cloud passes."""

    n_pred: int
    n_source: int
    frame: str
    seed: int
    n_ref: int
    source_voxels: int
    rms_ratio: float
    expected_rms_ratio: float
    expected_rms_source: str
    ref_rms_ratio: float
    voxel_cover: float
    ref_voxel_cover: float
    ref_voxel_cover_sd: float
    cover_ratio: float
    voxel_outside: float
    ref_voxel_outside: float
    weighted_reference: bool = False
    flags: list[str] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return not self.flags

    def as_dict(self) -> dict:
        return asdict(self) | {"ok": self.ok}


def cloud_guard(pred_coords, source_coords, seed: int = 0, expected_rms_ratio: float | None = None,
                frame: str = "source", n_ref: int = GUARD_REF_DRAWS, source_weights=None) -> GuardReport:
    """Check a submitted cloud against the cloud it was drawn from, on the scorer's 16^3 canonical grid.

    The source is canonicalised as veckit does (centre, right-handed PCA basis, divide by its RMS radius) and its
    occupied voxels are the target set. frame="source" (default) centres the cloud and each reference draw on its own
    centroid, rotates them by the source's basis and divides the draws by the source's RMS radius and the cloud by
    that times expected / reference RMS ratio, so a translated copy or an intended isotropic rescale (passed as
    `expected_rms_ratio`) compares like a reference draw, while a slab of sections keeps its true footprint; the basis
    is the source's, so nearly degenerate principal axes cannot swap.
    It assumes the cloud was not rotated relative to the source. frame="own" gives the cloud its own basis and RMS
    and takes the best of the four proper sign flips, as occupancy_dice does (for rotated or reshaped clouds; scale
    then drops out). The same measures are taken on n_ref random draws of the same count from the source (rows
    repeated evenly when the cloud is larger than the source, as `sample_cells` does). Pass as source the pool the
    cells were drawn from (after any scope filter) in the submission's frame: the pool itself before a size step, or
    the pool put through the same size step. For a composition-weighted draw pass the per-row sampling weights as
    `source_weights`: the reference draws then use them (without replacement while the weighted rows suffice).

    Flags: "block" when coverage < 0.9 x the random reference's (a contiguous block of rows, i.e. a slab of sections,
    covers about 0.3 x at 1x and 0.55 x at 8x on released stages; under frame="source" a rotated or reshaped cloud also
    lands here); "rms" when |log(RMS ratio / expected)| > 0.03, the expected ratio being `expected_rms_ratio` or else
    the random reference's (about 1); "outside" when the share of occupied voxels outside the source's exceeds the
    reference's by more than 0.1.
    """
    P = _as_coords(pred_coords, "pred_coords")
    S = _as_coords(source_coords, "source_coords")
    if frame not in GUARD_FRAMES:
        raise ValueError(f"frame must be one of {GUARD_FRAMES}, got {frame!r}")
    if min(P.shape[0], S.shape[0]) < MIN_CLOUD:
        raise ValueError(f"the guard needs at least {MIN_CLOUD} points in each cloud")
    seed = _whole(seed, "seed", 0)
    n_ref = _whole(n_ref, "n_ref", 1)
    if expected_rms_ratio is not None and not (np.isfinite(expected_rms_ratio) and expected_rms_ratio > 0):
        raise ValueError(f"expected_rms_ratio must be finite and positive, got {expected_rms_ratio!r}")
    p = None
    if source_weights is not None:
        w = np.asarray(source_weights, dtype=np.float64)
        if w.shape != (S.shape[0],) or not np.isfinite(w).all() or (w < 0).any() or w.sum() <= 0:
            raise ValueError("source_weights must be one finite, non-negative weight per source row, not all zero")
        p = w / w.sum()
    rms_s, rms_p = rms_radius(S), rms_radius(P)
    if rms_p < DEGENERATE_RMS or rms_s < DEGENERATE_RMS:
        raise ValueError("a cloud without spread cannot be canonicalised")
    basis = _right_handed_basis(S)
    vox_s = _voxel_ids(_canonical(S, basis))

    def measure(Q: np.ndarray, ratio: float) -> tuple[float, float]:
        if frame == "source":
            Zs = [_canonical(Q, basis, rms_s * ratio)]
        else:
            Z = _canonical(Q, _right_handed_basis(Q))
            Zs = [Z * np.asarray(f, dtype=np.float64) for f in PROPER_FLIPS]
        best = (-1.0, 1.0)
        for Z in Zs:
            v = _voxel_ids(Z)
            inter = np.intersect1d(v, vox_s, assume_unique=True).size
            if inter / vox_s.size > best[0]:
                best = (inter / vox_s.size, 1.0 - inter / v.size)
        return best

    rng = np.random.default_rng([seed, 2])
    refs = []
    for _ in range(n_ref):
        if p is None:
            Q = S[draw_balanced(S.shape[0], P.shape[0], rng)]
        else:
            Q = S[rng.choice(S.shape[0], P.shape[0], replace=P.shape[0] > int((p > 0).sum()), p=p)]
        refs.append((*measure(Q, 1.0), rms_radius(Q) / rms_s))
    ref = np.asarray(refs)
    ref_cover, ref_outside, ref_rms = (float(x) for x in ref.mean(0))
    if expected_rms_ratio is None:
        expected, expected_source = ref_rms, "random reference"
    else:
        expected, expected_source = float(expected_rms_ratio), "argument"
    cover, outside = measure(P, expected / ref_rms)  # back to the scale of the reference draws
    rms_ratio = rms_p / rms_s
    cover_ratio = cover / ref_cover if ref_cover > 0 else float("nan")
    flags = []
    if not cover_ratio >= GUARD_COVER_RATIO:
        flags.append(f"block: covers {cover:.1%} of the source's occupied canonical voxels against {ref_cover:.1%} for "
                     f"random subsets of the same count (ratio {cover_ratio:.2f} < {GUARD_COVER_RATIO}): a contiguous "
                     "block of rows (a slab of sections), or a rotated or reshaped cloud, looks like this")
    if abs(np.log(rms_ratio / expected)) > GUARD_RMS_LOG_TOL:
        flags.append(f"rms: RMS radius is {rms_ratio:.3f} x the source's, expected {expected:.3f} ({expected_source}; "
                     f"tolerance {GUARD_RMS_LOG_TOL} in log)")
    if outside - ref_outside > GUARD_OUTSIDE_TOL:
        flags.append(f"outside: {outside:.1%} of occupied voxels lie outside the source's occupied region against "
                     f"{ref_outside:.1%} for random subsets")
    return GuardReport(
        n_pred=int(P.shape[0]), n_source=int(S.shape[0]), frame=frame, seed=seed, n_ref=n_ref,
        source_voxels=int(vox_s.size), rms_ratio=float(rms_ratio), expected_rms_ratio=float(expected),
        expected_rms_source=expected_source, ref_rms_ratio=ref_rms, voxel_cover=float(cover),
        ref_voxel_cover=ref_cover, ref_voxel_cover_sd=float(ref[:, 0].std(ddof=1)) if n_ref > 1 else 0.0,
        cover_ratio=float(cover_ratio), voxel_outside=float(outside), ref_voxel_outside=ref_outside,
        weighted_reference=p is not None, flags=flags)


# ---------------------------------------------------------------------------------------------------------------------
# Builder
def _panel_matrix(board: Board, cells, genes: list[str] | None) -> tuple[np.ndarray, list[str]]:
    """cells x board genes, float32, from an AnnData (columns by var_names) or a matrix / vec.evaluate.Stage whose
    columns are named by `genes` (or, with genes=None, already in the board's panel order)."""
    notes = []
    if isinstance(cells, ad.AnnData):
        if genes is not None:
            raise ValueError("an AnnData carries its gene names in var_names; do not pass genes")
        names, M = [str(g) for g in cells.var_names], cells.X
    else:
        names = None if genes is None else [str(g) for g in genes]
        M = cells.X if not (sparse.issparse(cells) or isinstance(cells, np.ndarray)) and hasattr(cells, "X") else cells
    if not sparse.issparse(M):
        M = np.asarray(M)
    if M.ndim != 2:
        raise ValueError(f"expression must be a cells x genes matrix, got shape {M.shape}")
    if names is None:
        if M.shape[1] != board.n_genes:
            raise ValueError(f"{M.shape[1]} columns and no gene names: pass genes= (the board panel has "
                             f"{board.n_genes} genes)")
        notes.append("columns taken as already in the board's panel order (no gene names given)")
        X = dense(M)
    else:
        if len(names) != M.shape[1]:
            raise ValueError(f"{len(names)} gene names for {M.shape[1]} columns")
        if len(set(names)) != len(names):
            raise ValueError("gene names are not unique")
        pos = {g: i for i, g in enumerate(names)}
        missing = [g for g in board.genes if g not in pos]
        if missing:
            raise KeyError(f"{len(missing)} board genes missing from the cells, e.g. {missing[:5]}")
        if len(names) > board.n_genes:
            notes.append(f"{len(names) - board.n_genes} genes outside the board panel dropped")
        X = dense(M[:, [pos[g] for g in board.genes]])
    return np.asarray(X, dtype=np.float32), notes


def build_submission(board: Board, cells, coords, genes: list[str] | None = None, out: str | Path | None = None,
                     allow_duplicate_coords: bool = False,
                     allow_below_truth: bool = False) -> tuple[ad.AnnData, CheckReport]:
    """A submission AnnData for a coordinate board: X (log1p CP10k, never renormalised) in the board's gene order,
    obsm["spatial_3D"] float32, obs names cell_0..cell_{n-1}, nothing else.

    cells: an AnnData (columns found by var_names), or a matrix or vec.evaluate.Stage with column names in `genes`.
    The file is written to `out` (or a temporary file) and must pass vec.io.check_submission; a rejected file is
    deleted and the call raises. Refused before writing: a board without coordinates, a count outside [min_cells,
    max_cells] or below truth_cells (unless allow_below_truth, for a control), negative or non-finite expression (never
    clipped here), coordinates that are not (n, 3) and finite, and exact duplicate positions (as the scorer reads them,
    float32) unless allow_duplicate_coords. Returns the AnnData and the check report, with warnings added when X does
    not look like log1p CP10k: values above log1p(1e4) or whole numbers only (raw counts), or a median implied library
    size (expm1 of X summed per cell) outside [0.1, 2] x 1e4 (log1p of counts never scaled to 1e4 per cell; the low
    end leaves room for a board panel that is part of the measured one).
    """
    if not board.needs_coords:
        raise ValueError(f"{board.key}: this board takes no coordinates; build its file with vec.io.make_submission")
    X, notes = _panel_matrix(board, cells, genes)
    C = _as_coords(coords)
    n = X.shape[0]
    if C.shape[0] != n:
        raise ValueError(f"{C.shape[0]} coordinate rows for {n} cells")
    if n < board.min_cells or (board.max_cells is not None and n > board.max_cells):
        raise ValueError(f"{board.key}: {n} cells outside [{board.min_cells}, {board.max_cells}]")
    if board.truth_cells and n < board.truth_cells and not allow_below_truth:
        raise ValueError(f"{board.key}: {n} cells is below the board's truth_cells; the count rule never goes below it "
                         "(use n_cells), or pass allow_below_truth=True for a control")
    if not np.isfinite(X).all():
        raise ValueError("expression contains NaN or infinite values")
    if (X < 0).any():
        raise ValueError("expression contains negative values (centred upstream?)")
    C32 = C.astype(np.float32)
    if not np.isfinite(C32).all():
        raise ValueError("coordinates overflow float32")
    n_dup = n - np.unique(C32, axis=0).shape[0]
    if n_dup and not allow_duplicate_coords:
        raise ValueError(f"{n_dup} cells repeat the exact (float32) coordinates of another cell: break them with "
                         "sample_cells' jitter, or pass allow_duplicate_coords=True")
    a = make_submission(X, board.genes, C32)
    if out is None:
        with tempfile.TemporaryDirectory() as tmp:
            rep = check_submission(write_submission(a, Path(tmp) / "submission.h5ad", stamp=False), board)
    else:
        path = write_submission(a, out)
        rep = check_submission(path, board)
        if not rep.ok:
            path.unlink()
            stamp_path(path).unlink(missing_ok=True)
    if not rep.ok:
        raise ValueError(f"{board.key}: vec.io.check_submission rejected the file: {rep.errors}")
    # float32 storage can put a valid single-gene cell a hair above the float64 log1p(1e4); allow that much
    vmax, ceiling = (float(X.max()) if X.size else 0.0), float(np.log1p(TARGET_LIBRARY_SIZE)) * (1 + 1e-6)
    nz = X[X > 0]
    whole = bool(nz.size) and bool(np.all(nz == np.rint(nz)))
    if vmax > ceiling or whole:
        rep.warnings.append(f"X does not look like log1p CP10k (max {vmax:.3g} against the largest possible "
                            f"{ceiling:.3g}; whole numbers only: {whole}): raw counts?")
    # log1p of counts never scaled to 1e4 per cell is fractional and stays below the ceiling: only its size shows
    with np.errstate(over="ignore"):
        lib = float(np.median(np.expm1(X).sum(1, dtype=np.float64)))
    lo, hi = (f * TARGET_LIBRARY_SIZE for f in LIBRARY_BAND)
    if not lo <= lib <= hi:
        rep.warnings.append(f"median implied library size {lib:.3g} (expm1 of X summed per cell) is outside "
                            f"[{lo:.3g}, {hi:.3g}]: log1p of counts never scaled to 1e4 per cell?")
    rep.warnings.extend(notes)
    return a, rep


# ---------------------------------------------------------------------------------------------------------------------
# Interpolation between two staged clouds (feat-010, the T2 embryo pipeline; docs/PLAN.md 4.2)
# Fixed on 2026-11-01, before any interpolation result: the plan's setting where it names one, else sizes chosen for
# run time, not for score.
INTERP_ALPHA = 0.5  # FGW structure weight: the middle of the plan's 0.3-0.6
INTERP_N_ANCHOR = 2000  # cells per stage in the coupling (POT's exact solver)
INTERP_N_PC = 30  # expression features: principal components of the log1p panel expression of both stages' anchors
INTERP_K = 8  # a cell moves by the mean displacement of its 8 nearest anchors of its own stage
INTERP_MAX_ITER = 200  # conditional-gradient iterations of the coupling
INTERP_EMD_ITER = 10_000_000  # network-simplex iterations of each step's exact transport (POT's default stops short)
PARTNER_CHUNK = 4096  # cells per block when drawing coupled partners (memory only; the draw does not depend on it)
HERD_STREAMS = (6, 7)  # random streams of the herded draws (before, after), apart from the interpolation's 0-5


def pchip_rms(days, rms, target_day: float) -> float:
    """The pre-registered embryo size rule in raw mode (experiments/coords/size_procedure.py, `_embryo_interp`): PCHIP
    of log full-cloud RMS radius against the stage day through the given stages, read at target_day, which must lie
    strictly inside them."""
    from scipy.interpolate import PchipInterpolator

    d = np.asarray(days, dtype=np.float64)
    r = np.asarray(rms, dtype=np.float64)
    order = np.argsort(d)
    d, r = d[order], r[order]
    if len(np.unique(d)) != len(d):
        raise ValueError("two stages share a day")
    if not d[0] < target_day < d[-1]:
        raise ValueError(f"target day {target_day} is not strictly inside the stages {d.tolist()}")
    if not (np.isfinite(r).all() and (r > 0).all()):
        raise ValueError("RMS radii must be finite and positive")
    return float(np.exp(PchipInterpolator(d, np.log(r), extrapolate=False)(target_day)))


@dataclass
class Interpolation:
    """`interpolate_clouds` result: output cell i is row `rows[i]` of stage `stage[i]` (0 before, 1 after), at
    `coords[i]` in the before stage's frame. `pool_coords` and `pool_weights`: where every cell of both stages lands
    (before rows, then after rows) and its chance of being drawn, the source for the cloud guard. `partners` (only
    with `with_partners=True`): for each output cell, a row of the other stage coupled to it (`blend_expression`)."""

    stage: np.ndarray
    rows: np.ndarray
    coords: np.ndarray
    pool_coords: np.ndarray
    pool_weights: np.ndarray
    record: dict
    partners: np.ndarray | None = None


def _joint_pcs(Eb: np.ndarray, Ea: np.ndarray, n_pc: int) -> tuple[np.ndarray, np.ndarray]:
    X = np.vstack([Eb, Ea]).astype(np.float64)
    X -= X.mean(0)
    _, _, Vt = np.linalg.svd(X, full_matrices=False)
    P = X @ Vt[: min(n_pc, Vt.shape[0])].T
    return P[: len(Eb)], P[len(Eb):]


def interpolate_clouds(before_coords, after_coords, before_expr, after_expr, t: float, n: int, seed: int, *,
                       target_rms: float, alpha: float = INTERP_ALPHA, n_anchor: int = INTERP_N_ANCHOR,
                       n_pc: int = INTERP_N_PC, k: int = INTERP_K, max_iter: int = INTERP_MAX_ITER,
                       emd_max_iter: int = INTERP_EMD_ITER, t_pos: float | None = None,
                       with_partners: bool = False, herd: bool = False) -> Interpolation:
    """Displacement interpolation between two stages at time fraction t (0 = before, 1 = after).

    1. Each cloud is centred and divided by its RMS radius (no PCA axes, so no sign choice).
    2. A fused Gromov-Wasserstein coupling (vec.core) of n_anchor random cells per stage: features are the first n_pc
       principal components of both stages' anchors' log1p panel expression, structure the within-cloud distances.
       The proper rotation that maps before onto after under the coupling brings the after cloud into the before
       stage's frame (a rotation, never a mirror image).
    3. Each anchor's displacement to its coupled partners' mean; every cell of either stage moves by the mean
       displacement of its k nearest anchors of its own stage: a before cell t of the way towards after, an after cell
       1 - t of the way back (McCann interpolation, carried to every cell).
    4. n cells drawn: round((1 - t) n) before cells and the rest after cells, random rows (vec.core.draw_balanced);
       each keeps its real expression, so the population is a mixture of real cells of the two stages.
    5. The drawn cloud is scaled about its centroid to RMS radius target_rms (the size rule) and placed at the before
       stage's centroid; exact repeats (a row drawn twice) move by the free jitter, as in `sample_cells`.

    t_pos (default: t) moves the cells by a different fraction than the mixture uses. 0 gives the before stage's
    geometry, every after cell carried onto its partners in the before cloud; 1 gives the after stage's geometry.
    These are the embryo s2 shape variants (`docs/S2_RULES.md`). Left unset, the result is bit-identical to t_pos = t.

    with_partners: also pick, for every output cell, a coupled cell of the other stage (one of its k nearest anchors at
    random, then a partner anchor from that anchor's coupling row or column), from its own random stream, so the rest
    of the result is unchanged.

    herd (embryo try E9 of experiments/harness/TRIES.md): each stage's drawn rows are herded towards that stage's mean
    expression (`vec.core.herd_rows` on before_expr / after_expr, streams HERD_STREAMS) before positions are
    read; the coupling and the moves are unchanged.
    """
    Cb, Ca = _as_coords(before_coords, "before_coords"), _as_coords(after_coords, "after_coords")
    Eb, Ea = np.asarray(before_expr, dtype=np.float64), np.asarray(after_expr, dtype=np.float64)
    if Eb.shape[0] != Cb.shape[0] or Ea.shape[0] != Ca.shape[0] or Eb.shape[1] != Ea.shape[1]:
        raise ValueError("expression must have one row per cell and the same genes in both stages")
    if not 0.0 <= t <= 1.0:
        raise ValueError(f"t must lie in [0, 1], got {t}")
    tp = t if t_pos is None else float(t_pos)
    if not 0.0 <= tp <= 1.0:
        raise ValueError(f"t_pos must lie in [0, 1], got {t_pos}")
    if not (np.isfinite(target_rms) and target_rms > 0):
        raise ValueError(f"target_rms must be finite and positive, got {target_rms}")
    n, seed = _whole(n, "n", 1), _whole(seed, "seed", 0)
    cb, sb, sa = Cb.mean(0), rms_radius(Cb), rms_radius(Ca)
    Zb, Za = (Cb - cb) / sb, (Ca - Ca.mean(0)) / sa
    rng = np.random.default_rng([seed, 4])
    ib = np.sort(rng.choice(len(Zb), min(n_anchor, len(Zb)), replace=False))
    ia = np.sort(rng.choice(len(Za), min(n_anchor, len(Za)), replace=False))
    Fb, Fa = _joint_pcs(Eb[ib], Ea[ia], n_pc)
    T = fgw_coupling(Fb, Fa, Zb[ib], Za[ia], alpha, max_iter, emd_max_iter)
    R = coupling_rotation(Zb[ib], Za[ia], T)
    Za_b = Za @ R
    db = barycentric(T, Za_b[ia]) - Zb[ib]
    da = barycentric(T.T, Zb[ib]) - Za_b[ia]
    near_b, near_a = nearest_anchors(Zb, Zb[ib], k)[1], nearest_anchors(Za_b, Za_b[ia], k)[1]
    pos_b = Zb + tp * db[near_b].mean(1)
    pos_a = Za_b + (1.0 - tp) * da[near_a].mean(1)
    n_b = int(round((1.0 - t) * n))
    n_a = n - n_b
    empty = np.array([], dtype=np.int64)
    rows_b = draw_balanced(len(Zb), n_b, np.random.default_rng([seed, 0])) if n_b else empty
    rows_a = draw_balanced(len(Za), n_a, np.random.default_rng([seed, 1])) if n_a else empty
    herd_rec = None
    if herd:
        herd_rec = {}
        if n_b:
            rows_b, herd_rec["before"] = herd_rows(Eb, np.arange(len(Zb)), rows_b,
                                                   np.random.default_rng([seed, HERD_STREAMS[0]]))
        if n_a:
            rows_a, herd_rec["after"] = herd_rows(Ea, np.arange(len(Za)), rows_a,
                                                  np.random.default_rng([seed, HERD_STREAMS[1]]))
    stage = np.r_[np.zeros(n_b, dtype=np.int64), np.ones(n_a, dtype=np.int64)]
    rows = np.r_[rows_b, rows_a].astype(np.int64)
    perm = np.random.default_rng([seed, 2]).permutation(n)
    stage, rows = stage[perm], rows[perm]
    Z = np.empty((n, 3))
    Z[stage == 0], Z[stage == 1] = pos_b[rows[stage == 0]], pos_a[rows[stage == 1]]
    zc = Z.mean(0)
    scale = target_rms / rms_radius(Z)
    coords = cb + (Z - zc) * scale
    pool = cb + (np.vstack([pos_b, pos_a]) - zc) * scale
    weights = np.r_[np.full(len(Zb), (1.0 - t) / len(Zb)), np.full(len(Za), t / len(Za))]
    _, first = np.unique(coords.astype(np.float32), axis=0, return_index=True)
    dup = np.ones(n, dtype=bool)
    dup[first] = False
    if dup.any():
        coords[dup] += _ball_shifts(int(dup.sum()), JITTER_NN_FRAC * median_nn_distance(coords),
                                    np.random.default_rng([seed, 3]))
    angle = float(np.degrees(np.arccos(np.clip((np.trace(R) - 1.0) / 2.0, -1.0, 1.0))))
    record = {"t": float(t), "t_pos": float(tp), "n": n, "n_before": n_b, "n_after": n_a, "alpha": float(alpha),
              "n_anchor": [len(ib), len(ia)], "n_pc": int(Fb.shape[1]), "k": int(k), "max_iter": int(max_iter),
              "rotation_deg": angle,
              "coupling_nonzero": int((T > 0).sum()), "target_rms": float(target_rms),
              "rms_before": sb, "rms_after": sa, "scale_from_canonical": float(scale),
              "mean_anchor_shift": [float(np.linalg.norm(db, axis=1).mean()), float(np.linalg.norm(da, axis=1).mean())],
              "n_jittered": int(dup.sum())}
    if herd_rec is not None:
        record["herd"] = herd_rec
    partners = None
    if with_partners:
        prng = np.random.default_rng([seed, 5])
        partners = np.empty(n, dtype=np.int64)
        for s, near, other, P in ((0, near_b, ia, T), (1, near_a, ib, T.T)):
            m = stage == s
            if m.any():
                anchor = near[rows[m], prng.integers(0, near.shape[1], int(m.sum()))]
                partners[m] = other[_sample_from_rows(P, anchor, prng)]
        record["partners"] = "one coupled cell of the other stage per output cell (stream [seed, 5])"
    return Interpolation(stage=stage, rows=rows, coords=coords, pool_coords=pool, pool_weights=weights, record=record,
                         partners=partners)


def _sample_from_rows(P: np.ndarray, which: np.ndarray, rng: np.random.Generator,
                      chunk: int = PARTNER_CHUNK) -> np.ndarray:
    """For each entry r of `which`, a column of P drawn with probability P[r, j] / P[r].sum()."""
    cum = np.cumsum(np.asarray(P, dtype=np.float64), axis=1)
    if (cum[:, -1] <= 0).any():
        raise ValueError("a row of the coupling carries no mass")
    cum /= cum[:, -1:]
    out = np.empty(len(which), dtype=np.int64)
    u = rng.random(len(which))
    for lo in range(0, len(which), chunk):
        c = cum[which[lo:lo + chunk]]
        out[lo:lo + chunk] = np.minimum((c < u[lo:lo + chunk, None]).sum(1), P.shape[1] - 1)
    return out


# ---------------------------------------------------------------------------------------------------------------------
# Detection depth of an interpolated stage (2026-11-03; embryo tries E7 and E8, experiments/harness/TRIES.md, each
# registered before its run). The released MERFISH stages differ several-fold in counts per cell (median totals:
# E6.75 30, E7.25 35, E7.5 62, E8.0 212), and the scorer reads log(1 + 1e4 c / N), whose zeros and levels follow the
# depth. A mixture of whole cells of two brackets carries both brackets' depths. The interpolated stage's depth is
# predicted log-linearly between the brackets' median depths at the bracket fraction t; a deeper bracket's cells are
# binomially thinned to it (a shallower bracket cannot be deepened and keeps its cells as they are).
DEPTH_RULE = ("depth: the interpolated stage's median counts per cell = exp((1 - t) log d_before + t log d_after), "
              "d = a bracket's median counts per cell over all its cells; cells of a bracket deeper than that are "
              "binomially thinned by target / d (panel genes and the rest of the library alike), a shallower "
              "bracket's cells are left as they are")
DEPTH_ZERO_RETRIES = 20  # a thinned cell left with no panel counts is thinned again (a cell without counts is not kept)


def library_sizes(X, counts) -> np.ndarray:
    """Each cell's library size N in `.X = log(1 + 1e4 c / N)`. N counts every gene of the cell's file, also genes
    outside the panel, so it is read from the cell's expressed panel genes: N = 1e4 c / expm1(x) (the mean over them,
    rounded; they agree to the rounding of X)."""
    Xc, C = sparse.csr_matrix(X, dtype=np.float64), sparse.csr_matrix(counts, dtype=np.float64)
    for M in (Xc, C):
        M.eliminate_zeros()
        M.sort_indices()
    if Xc.shape != C.shape or not (np.array_equal(Xc.indptr, C.indptr) and np.array_equal(Xc.indices, C.indices)):
        raise ValueError("expression and counts must have the same non-zero entries")
    per_row = np.diff(C.indptr)
    if (per_row == 0).any():
        raise ValueError(f"{int((per_row == 0).sum())} cells have no expressed panel gene")
    n = TARGET_LIBRARY_SIZE * C.data / np.expm1(Xc.data)
    lib = np.rint(np.add.reduceat(n, C.indptr[:-1]) / per_row)
    if (lib < np.asarray(C.sum(1)).ravel() - 0.5).any():
        raise ValueError("a recovered library size is below the cell's panel counts")
    return lib


def depth_target(d_before: float, d_after: float, t: float) -> float:
    """The interpolated stage's median counts per cell: log-linear between the brackets' medians at fraction t."""
    if not (d_before > 0 and d_after > 0 and 0.0 <= t <= 1.0):
        raise ValueError("depths must be positive and t in [0, 1]")
    return float(np.exp((1.0 - t) * np.log(d_before) + t * np.log(d_after)))


def _binomial_rows(C: sparse.csr_matrix, other: np.ndarray, w: np.ndarray,
                   rng: np.random.Generator) -> tuple[sparse.csr_matrix, np.ndarray]:
    """Each row's panel counts and its counts outside the panel, binomially thinned by the row's weight in [0, 1]."""
    data = rng.binomial(C.data.astype(np.int64), np.repeat(w, np.diff(C.indptr)))
    out = sparse.csr_matrix((data, C.indices.copy(), C.indptr.copy()), shape=C.shape)
    return out, rng.binomial(np.asarray(other, dtype=np.int64), w)


def _log_normalise(C: sparse.csr_matrix, lib: np.ndarray) -> np.ndarray:
    """log(1 + 1e4 c / N) as the released files store it, dense float32."""
    X = sparse.diags(TARGET_LIBRARY_SIZE / np.asarray(lib, dtype=np.float64)) @ sparse.csr_matrix(C, dtype=np.float64)
    X.data = np.log1p(X.data)
    return np.asarray(X.toarray(), dtype=np.float32)


def _stage_counts(a: ad.AnnData, genes) -> tuple[np.ndarray, sparse.csr_matrix, np.ndarray, np.ndarray]:
    """A bracket's panel expression (dense float32), panel counts, library sizes and counts outside the panel."""
    sub = a[:, list(genes)]
    C = sparse.csr_matrix(sub.layers["counts"], dtype=np.int64)
    C.eliminate_zeros()
    C.sort_indices()
    lib = library_sizes(sub.X, C)
    return dense(sub.X).astype(np.float32), C, lib, lib - np.asarray(C.sum(1)).ravel()


def _keep_some_counts(make, n_rows: int, rng: np.random.Generator) -> tuple[sparse.csr_matrix, np.ndarray, int]:
    """make(rows, rng) -> (panel counts, library sizes) for those rows; rows left without panel counts are made again,
    up to DEPTH_ZERO_RETRIES times. Returns the counts, the library sizes and how many rows needed another try."""
    C, lib = make(np.arange(n_rows), rng)
    C, lib = sparse.csr_matrix(C), np.array(lib, dtype=np.float64)
    retried = 0
    for _ in range(DEPTH_ZERO_RETRIES):
        empty = np.flatnonzero(np.asarray(C.sum(1)).ravel() == 0)
        if not len(empty):
            break
        retried += len(empty)
        Cn, ln = make(empty, rng)
        keep = np.ones(n_rows, dtype=bool)
        keep[empty] = False
        order = np.r_[np.flatnonzero(keep), empty]  # the original row of each stacked row
        C = sparse.vstack([C[keep], sparse.csr_matrix(Cn)]).tocsr()[np.argsort(order)]
        lib[empty] = ln
    if (np.asarray(C.sum(1)).ravel() == 0).any():
        raise ValueError("thinning left a cell without panel counts after every retry")
    return C, lib, retried


def depth_matched_expression(before: ad.AnnData, after: ad.AnnData, interp: Interpolation, genes, t: float,
                             seed) -> tuple[np.ndarray, dict]:
    """The interpolation's cells with each bracket's cells thinned to the predicted depth (DEPTH_RULE): expression
    (cells x genes, float32, log(1 + 1e4 c / N)) and a record. A bracket at or below the predicted depth keeps its
    rows exactly as released."""
    rng = np.random.default_rng(seed)
    parts = [_stage_counts(a, genes) for a in (before, after)]
    d = [float(np.median(p[2])) for p in parts]
    target = depth_target(d[0], d[1], t)
    out = np.empty((len(interp.rows), len(genes)), dtype=np.float32)
    rec = {"rule": DEPTH_RULE, "t": float(t), "median_depth": {"before": d[0], "after": d[1]}, "target": target,
           "factor": {}, "thinned": {}, "retried": {}}
    for s, name in ((0, "before"), (1, "after")):
        m = interp.stage == s
        E, C, lib, other = parts[s]
        f = min(1.0, target / d[s])
        rec["factor"][name] = f
        rows = interp.rows[m]
        if f >= 1.0 or not m.any():
            out[m] = E[rows]
            rec["thinned"][name], rec["retried"][name] = 0, 0
            continue
        w = np.full(len(rows), f)

        def make(idx, g, C=C, other=other, rows=rows, w=w):
            Cn, on = _binomial_rows(C[rows[idx]], other[rows[idx]], w[idx], g)
            return Cn, np.asarray(Cn.sum(1)).ravel() + on

        Cn, ln, retried = _keep_some_counts(make, len(rows), rng)
        out[m] = _log_normalise(Cn, ln)
        rec["thinned"][name], rec["retried"][name] = int(m.sum()), retried
    return out, rec


# Detection depth of an extrapolated stage (2026-11-03; heart try H3, experiments/harness/TRIES.md, registered before
# its run). Counts per cell fall steadily through heart development in the released series (median counts per cell:
# E8.25_late 233, E8.5 166, E8.75 130, E9.5 75, E10.5 56).
DEPTH_EXTRAP_RULE = ("depth (extrapolation): the target's median counts per cell = d_ref * exp(damping * r * "
                     "(target_day - ref_day)), r = (log d_ref - log d_start) / (ref_day - start_day), the latest log "
                     "rate per day, at B0's damping; the reference's cells are binomially thinned by target / d_ref "
                     "when the target lies below d_ref (counts are never added)")


def median_depth(a: ad.AnnData, genes) -> float:
    """A stage's median counts per cell (each cell's library size, `library_sizes`), over all its cells."""
    sub = a[:, list(genes)]
    return float(np.median(library_sizes(sub.X, sub.layers["counts"])))


def extrapolated_depth(d_start: float, d_ref: float, start_day: float, ref_day: float, target_day: float,
                       damping: float) -> float:
    """DEPTH_EXTRAP_RULE's target depth."""
    if not (d_start > 0 and d_ref > 0 and ref_day > start_day and target_day >= ref_day and np.isfinite(damping)):
        raise ValueError("depths must be positive, the start before the reference and the target not before it")
    rate = (np.log(d_ref) - np.log(d_start)) / (ref_day - start_day)
    return float(d_ref * np.exp(damping * rate * (target_day - ref_day)))


def thin_to_depth(a: ad.AnnData, genes, d_ref: float, target: float, seed) -> tuple[np.ndarray, dict]:
    """a's cells binomially thinned by min(1, target / d_ref), panel genes and the rest of each library alike, cells
    left without panel counts thinned again (as in `depth_matched_expression`): expression (cells x genes, float32,
    log(1 + 1e4 c / N)) and a record. At a factor of 1 the cells come back exactly as given."""
    E, C, lib, other = _stage_counts(a, genes)
    f = min(1.0, target / d_ref) if d_ref > 0 else 1.0
    rec = {"rule": DEPTH_EXTRAP_RULE, "reference_depth": float(d_ref), "target": float(target), "factor": f}
    if f >= 1.0:
        return E, {**rec, "thinned": 0, "retried": 0}
    rng = np.random.default_rng(seed)
    w = np.full(C.shape[0], f)

    def make(idx, g):
        Cn, on = _binomial_rows(C[idx], other[idx], w[idx], g)
        return Cn, np.asarray(Cn.sum(1)).ravel() + on

    Cn, ln, retried = _keep_some_counts(make, C.shape[0], rng)
    return _log_normalise(Cn, ln), {**rec, "thinned": int(C.shape[0]), "retried": retried}


def blend_expression(before: ad.AnnData, after: ad.AnnData, interp: Interpolation, genes, t: float,
                     seed) -> tuple[np.ndarray, dict]:
    """Each output cell's counts become Binomial(c_b, w_b) + Binomial(c_a, w_a) of itself and its coupled partner
    (`interp.partners`, from `interpolate_clouds(..., with_partners=True)`), with w_b = min(1, (1 - t) D / d_b) and
    w_a = min(1, t D / d_a): at the predicted depth D (DEPTH_RULE) the expected counts split (1 - t) : t between the
    earlier and the later cell. Then log(1 + 1e4 c / N). Returns the expression and a record."""
    if interp.partners is None:
        raise ValueError("blend_expression needs interpolate_clouds(..., with_partners=True)")
    rng = np.random.default_rng(seed)
    (Eb, Cb, lb, ob), (Ea, Ca, la, oa) = (_stage_counts(a, genes) for a in (before, after))
    db, da = float(np.median(lb)), float(np.median(la))
    target = depth_target(db, da, t)
    wb, wa = min(1.0, (1.0 - t) * target / db), min(1.0, t * target / da)
    is_b = interp.stage == 0
    rb = np.where(is_b, interp.rows, interp.partners)
    ra = np.where(is_b, interp.partners, interp.rows)

    def make(idx, g):
        B, obn = _binomial_rows(Cb[rb[idx]], ob[rb[idx]], np.full(len(idx), wb), g)
        A, oan = _binomial_rows(Ca[ra[idx]], oa[ra[idx]], np.full(len(idx), wa), g)
        S = (B + A).tocsr()
        S.eliminate_zeros()
        return S, np.asarray(S.sum(1)).ravel() + obn + oan

    C, lib, retried = _keep_some_counts(make, len(interp.rows), rng)
    rec = {"rule": "blend: " + DEPTH_RULE, "t": float(t), "median_depth": {"before": db, "after": da},
           "target": target, "weights": {"before": wb, "after": wa}, "retried": retried,
           "expected_share_before": float(wb * db / (wb * db + wa * da))}
    return _log_normalise(C, lib), rec

