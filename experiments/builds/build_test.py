#!/usr/bin/env python
"""Test-board candidate files, built with the harness's own builders in the test geometry (P3_RUNBOOK section 5).

Each entry is exactly the candidate the harness judged, run on the test's stages, then the test-only steps the
pseudo-tasks cannot run (the pre-registered sizes). Nothing is uploaded: the file is written with its build stamp and
checked locally (`vec check`); a server format check is a separate `vec submit ... --from-dgx` call from the local
machine, and an official upload needs the captain's go-ahead.

  t1      M1 (docs/PLAN.md 4.1; the captain's rulings of 2026-11-01): E10.5 cells re-weighted to E10.5's own chamber
          mix (or --cm-share, the cardiomyocyte share of heart cells), Feng's E9.5 -> E14.5 direction on the cells of
          Feng's groups at 2/5 of the bracket (E10.5 -> E12.5); the rest keep their E10.5 expression (R1 of
          experiments/harness/TRIES.md, the captain's ruling of 2026-11-05; --no-hold-unbracketed gives them B0's
          damped E9.5 -> E10.5 momentum, M1 as built before R1). R2 of TRIES.md also holds R2's gene modules flat in
          Feng's direction (the union N12 wrote, pinned in configs/t1_hygiene.yaml): on by default with R1 since the
          build commit after N12 (docs/PLAN.md decision log 2026-11-06 (k)); --no-feng-hygiene rebuilds the s1 of
          88b057f (M1 with R1, before R2). T12a and T12b (TRIES.md; the captain's rulings of 2026-11-10 (a) and (d)):
          --strong-null adds the strong-null down-tail prior after M1's steps (vec.t1_scope.apply_down_prior, genes by
          the reference's level alone), off by default and refused until configs/t1_strong_null.yaml says approved
          (a later build commit of the main session, after the stop rule and the captain's ruling). --g-mom (default
          off; only by the captain's ruling on try T11 of TRIES.md; refused together with --strong-null, whose order
          against it no row fixes) adds g_mom's within-type factors: fitted on E9.5 -> E10.5 and applied last, after
          Feng's direction, on the cells g_mom edits (metacell factors through apply_factors with keep_depth, R1's held
          labels left alone), so no share, cell count or scope changes. Row R3b (the captain's ruling of 2026-11-10
          (p)): the factors edit cardiomyocytes and endocardium only; epicardium and vascular endothelium keep M1's
          counts. --feng-groups none (default authors, every file as built; only the captain's ruling puts it in an
          upload) builds one file: the first upload's command (t1 --g-mom --base lib3) with Feng's direction taken
          out and nothing put in its place. Every cell keeps its E10.5 expression before g_mom's edit (R1 on every
          label), and no Feng file, R2 list or other outside data is read. g_mom's mask is the readout config's own
          (experiments/t1_direction/config.yaml): the factors edit the cm and endocardium cells whose label its
          regexes give a group, which on the labels configs/t1_chamber_classes.yaml places in those categories is
          every such cell but Early VIC, the same cells as in the entry. vec.t1_state_transport.R1_FENG_GROUPS, a
          list fixed in N12 by counting Feng cells, is not read under none
          (tested), and the hook's covered check is a self-check. Refused, before anything is read, without
          --g-mom and with --base pooled or mid, --cm-share, --blood-share, --lib3-pair, --herd,
          --no-hold-unbracketed, --feng-hygiene or --strong-null. The build record then carries a feng_groups field
          saying the direction step was removed, and its fraction is null. --released-bet share_trend | advance with
          --bet-strength (try T61 of experiments/harness/TRIES.md; off by default, --feng-groups none only): that file
          with its cells drawn at the chamber mix moved along the E9.5 -> E10.5-portions trend, or drawn inside each
          category towards the released momentum axis (vec.t1_bets). Both change composition, a lever the
          post-exposure block closes: refused unless --captain-ruling names the captain's ruling row for the arm.
          After a permitted-only relabelling whose record routes to recomputed, none is also refused unless
          --feng-none-ruling names a later decision-log entry of the captain (naming N54 and containing "sets the
          N54 route aside", not the run's own ruling, made after the record's commit) that sets the route aside;
          before any relabelling, none builds as before.
          --feng-groups recomputed (off unless the captain rules): the same one command, t1 --g-mom --base lib3, with
          Feng's direction formed on the groups experiments/external/feng_relabel.py recomputes from the permitted
          cells' raw counts (register rows X002 and X010), R1 and R2 as in the first upload. Refused, before
          anything is read, with any other base or option of that list (and with --no-feng-hygiene), and unless the
          labeller's record is committed at HEAD with the route recomputed; then refused unless the recomputed file
          passes feng_dirs.check_recomputed. R1's hold and g_mom's mask follow the groups with a recomputed direction
          after the fallback: a group without one moves no cell and g_mom leaves it alone. The record's feng_groups
          field names the file, its sha256, the yaml's sha256, the route, the fallbacks and the gate.
          The candidates fixed in docs/PLAN.md decision log 2026-11-20 (i) (off by default; no ruling; only the
          captain's ruling puts one in an upload), each a single change to that one file and refused on any other
          command or seed and on a dry run: --g-mom-damping 1.0 (candidate A: g_mom's factors finished at damping 1.0
          in place of B0_DAMPING 0.5, as practice arm gmom_s1 of experiments/t1_m1/bets_read.py finishes them) and
          --feng-fraction-mult 1.5 (candidate B: the direction step at 1.5 times the calendar fraction 2/5, as
          experiments/t1_m1/double_step_read.py multiplies it); both together are candidate C. Their defaults (0.5,
          1.0) build today's bytes; given at all, the record's t1_candidate field names both values.
          --feng-fraction-mult 2.0 --captain-ruling R20 (row R20 of experiments/harness/TRIES.md, the captain's
          order of 2026-11-21, docs/PLAN.md decision log 2026-11-21 (x)): the same file with the direction step at 2.0
          times the calendar fraction, nothing else changed; refused without that flag, with any other row, unless
          row R20 is the captain's order naming the option, with --g-mom-damping other than 0.5, with --released-bet,
          on any other command or seed and on a dry run. Its record's t1_candidate field names row R20.
          --feng-fraction-mult 2.5 or 3.0 --captain-ruling R22 (row R22 of experiments/harness/TRIES.md, the captain's
          order of 2026-11-22, docs/PLAN.md decision log 2026-11-22 (j): the first search line): the same file with the
          direction step at 2.5 or 3.0 times the calendar fraction, nothing else changed; refused as the 2.0 value is,
          with R22 in place of R20 (and 2.0 refused with R22). Its record's t1_candidate field names row R22.
          --g-mom-damping 1.5 or 2.0 --captain-ruling R22 (row R22, the second T1 search line of docs/PLAN.md decision
          log 2026-11-22): the same file with g_mom's within-type damping at 1.5 or 2.0, nothing else changed;
          refused without that flag, with any other row, unless row R22 is the captain's order naming the option,
          with --feng-fraction-mult other than 1.0, with --released-bet, on any other command or seed and on a dry
          run. Its record's t1_candidate field names row R22; --g-mom-damping 1.0 (candidate A) needs no ruling, as
          before
  heart   B0 (copy of E10.5 plus damped E9.5 -> E10.5 momentum in count space), sized by the pre-registered heart
          procedure: --entry s1 = scope hypothesis H1; s2_trim = the s2 of docs/S2_RULES.md (H2 as a trim to half
          E10.5's scope); s2 = H2 as first built (size only, superseded) (approved rule 3, `size_or_hold`). --trend
          (both heart files since the captain's ruling of 2026-11-13, docs/PLAN.md decision log (w); row H11 of
          experiments/harness/TRIES.md): b0_trend, each lineage's momentum the least-squares trend over E8.25_late to
          E10.5 in place of the E9.5 -> E10.5 step alone. --expression copy (off by default; --entry s1 only, refused
          with --trend, --depth, --herd and --cfm): the copy baseline's cells and expression (the draw B0 starts
          from, `vec.harness._copy`), no momentum. --z-mode raw (off by default): the size step's isotropic mode,
          every axis scaled by one factor to the hypothesis's target RMS, in place of the mode the size JSON
          recorded. The captain's direction for the heart's second upload (2026-11-18, docs/PLAN.md decision log
          2026-11-18 (k)): heart --entry s1 --expression copy --z-mode raw. The scope bet (no ruling; off unless all
          three flags are given; scope is a lever of the post-exposure block, so only the captain's ruling puts it in
          an upload): heart --entry s2_trim --expression copy --z-mode raw, the same two fixes on the half trim (H2,
          docs/S2_RULES.md), sized to the raw mode's own logged H2 figure. The count is s1's: vec.t2_build.n_cells
          does not read the pool's size (it only annotates the record), so the trimmed pool gets the count the server
          refused on 2026-11-18, and the file is halved with halve.py as the ruled file is. The full trim (no ruling;
          off unless all three flags are given; refused with --trend, --depth, --herd and --cfm): heart --entry
          core_trim --expression copy --z-mode raw, the same draw from E10.5 under radial_trim_mask at fraction 0,
          sized to the size JSON's raw-mode core forecast (z_modes.raw.pred_rms_core, the cloud at the core's own
          scope); its record and the half trim's carry scope_option (scope rule, kept-share field, size field,
          halving)
  embryo  --method pseudobulk_shift | fgw | fgw_shape0 | fgw_shape1: E7.5 cells moved by the E7.25 -> E7.5 change, or
          the FGW interpolation between E7.5 and E8.0 at t = 1/2 (fgw_shape0, the s2 of docs/S2_RULES.md, keeps the
          mixture in E7.5's geometry); sized by the pre-registered embryo procedure (size_embryo_E7.75.json).
          --positions floor_exact (off by default; with --method fgw_shape0 only; tries E47a and E47b of
          experiments/harness/TRIES.md): every cell on a real E7.5 position, the harness's fgw_shape0x, or with
          --expression series its fgw_shape0x_series
          --row-order copy_parity (off by default; no try, no ruling; branch campaign/embryo-occupancy): the same
          cells, expression and positions, the rows re-ordered (vec.t2_build.copy_parity_order) so that the file's
          pre-fix SVD parity on the rows the scorer's _match_n takes is the floor copy's (vec.evaluate.subsample of
          E7.5 as released, 0.1, seed 0, the bake protocol); refused unless the copy and the file get the same truth
          rows (a gate on the test board's truth_cells, read by code; nothing about it is printed). On a test board it
          is refused unless --row-order-test-ok (the captain's yes) is also given: a refusal there, or a success,
          shows one bit of the test count. It assumes the server scores the prediction's rows as the file holds them
          (file order, full count, no subsample or sort); score_h5ad's `_prediction` does this in the pinned copy,
          but the live server's handling is not shown. The as-built uploaded files' pre-fix parities are never
          computed, by any script, and the returned occupancy of a copy_parity upload is never compared with an
          earlier upload's to infer a mode (limit 2)
          --shape-interp stretch_ratio | stretch_registered | overlay (off by default: none; no try read, no ruling;
          row E49 as drafted, branch campaign/embryo-shape-interp): with --method fgw_shape0 --positions floor_exact
          only, one form of vec.shape_interp applied to the exact-position cloud between E7.5 and E8.0 at t = 1/2,
          before the size rule (so the pre-registered size and the count rule are untouched); the record gains a
          shape_interp block. On a test board it is refused unless --shape-interp-test-ok (the captain's yes) is also
          given. With --row-order copy_parity the order is fixed: form, size, jitter, cloud guard, then the row order
          on the final coordinates
          --copy-rows (off by default; the investigator's Option A, the captain's order of 2026-11-21, docs/PLAN.md
          decision log 2026-11-21 (h); row E50 as specified): with --method fgw_shape0 --positions floor_exact and
          --expression mixture or series only, the file takes the floor copy's own rows, order, count and float32
          coordinates (floor_copy_coords) and the expression of vec.t2_build.copy_rows_cloud; no size rule, centring,
          jitter or row-order step. It departs from the count rule (the record says "floor_copy" for n_cells on a test
          board, and its count_rule says departed) and may hold fewer cells than the truth (allow_below_truth). The
          limits are checked here, and every refusal on this path shows no value; the coordinates are checked bit for
          bit before writing and on the file read back. On a test board it is refused unless --copy-rows-test-ok (the
          captain's yes) is also given, before anything loads
          --copy-rows-size (off by default; the copy-rows file at the team's pre-registered size, the captain's order
          of 2026-11-21, docs/PLAN.md decision log 2026-11-21 (v); row E52 as specified): with --copy-rows and
          --expression mixture only. After the floor copy's bit check, the coordinates are scaled in float64 about the
          copy's own centroid to the embryo size JSON's preregistered.pred_rms_full (raw mode only, read as for the
          exact-position file; vec.t2_build.copy_rows_sized) and stored as float32; rows, order, count and expression
          stay the copy-rows file's, and nothing follows but writing (no jitter, centring, cloud guard or row-order
          step). Self-checks, each refusing with no value: the RMS within 1e-6 relative of the target, the prefix
          parity of the matched rows at metric seed 0 equal to the floor copy's for every count from 4 to the copy's,
          and the file read back equal to the float32 coordinates bit for bit (else the file and its stamp are
          removed). On a test board it needs --copy-rows-size-test-ok as well as --copy-rows-test-ok, before anything
          loads

The board is the test board from boards/ (after the intake, `vec boards update --from-file`); before that, --board
may name a validation board of the same task for a dry run, and the script refuses a dry run whose output file name
does not contain "dry" (such a file must never be uploaded). The build seed is never 0 (runbook section 5): default
vec.harness.BUILD_SEED_OFFSET + 1. T2 counts go through approved rule 7's committed decision
(`t2_build.n_cells_from_decision`, runbook section 5), checked against the count the harness's builder draws. On a
test board no cell count is printed or stored: the record says "max_cells" when the file holds the board's max_cells
(the count rule) and refuses otherwise (red-team review 2026-10-31; s1 build review 2026-11-01).

Usage (DGX, repo root, BLAS threads capped):
  uv run python experiments/builds/build_test.py t1 --out runs/builds/T1_s1.h5ad
  uv run python experiments/builds/build_test.py heart --entry s1 --trend --out runs/builds/T2heart_s1.h5ad
  uv run python experiments/builds/build_test.py embryo --method fgw --out runs/builds/T2embryo_s1.h5ad
  uv run python experiments/builds/build_test.py t1 --cm-share 0.60 --out runs/builds/T1_s2.h5ad
The commands of the entries as they stand are in experiments/builds/build_all.sh and docs/S2_RULES.md (since
2026-11-13 T1's two files are `t1 --g-mom --base lib3` and `t1 --g-mom --base mid --cm-share 0.60`, and both heart
files carry --trend: the captain's ruling of 2026-11-13, docs/PLAN.md decision log (w); row H11).
  uv run python experiments/builds/build_test.py heart --entry s2_trim --trend --out runs/builds/T2heart_s2.h5ad
  uv run python experiments/builds/build_test.py embryo --method fgw_shape0 --out runs/builds/T2embryo_s2.h5ad
Without --trend the heart builder gives the heart file as it was before that ruling (B0's E9.5 -> E10.5 step alone),
kept for reproducing the earlier reference set (builds_38616fe).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

from vec import harness, provenance, shape_interp, t1_bets, t1_scope, t2_build
from vec.boards import REPO_ROOT, load_boards
from vec.evaluate import SUBSAMPLE_FRAC, Stage, subsample
from vec.io import COORD_KEY, dense, read_stage
from vec.prep import load_stage

sys.path.insert(0, str(REPO_ROOT / "experiments" / "t1_m1"))
sys.path.insert(0, str(REPO_ROOT / "experiments" / "t1_hygiene"))
sys.path.insert(0, str(REPO_ROOT / "experiments" / "t1_strong_null"))
import feng_dirs  # noqa: E402
import hygiene  # noqa: E402
import strong_null  # noqa: E402

TEST_BOARDS = {"t1": "T1:test", "heart": "T2:heart:test_extrap", "embryo": "T2:embryo:test_interp"}
DRY_BOARDS = {"t1": {"T1:val"}, "heart": {"T2:heart:val_extrap", "T2:heart:val_interp"},
              "embryo": {"T2:embryo:val_interp"}}
COUNT_DECISION = REPO_ROOT / "experiments/density/results/release/count_decision.json"
COUNT_FIELDS = ("n", "n_before", "n_after", "n_jittered")  # the interpolation record's counts
SIZES = {"heart": REPO_ROOT / "experiments/coords/results/p3/size_heart_E12.5.json",
         "embryo": REPO_ROOT / "experiments/coords/results/p3/size_embryo_E7.75.json"}
DATA = harness.DATA_DIR
# fgw and its s2 shape variants (docs/S2_RULES.md): the position fraction passed to interpolate_clouds
EMBRYO_T_POS = {"fgw": None, "fgw_shape0": 0.0, "fgw_shape1": 1.0}
# R3b of experiments/harness/TRIES.md (the captain's ruling of 2026-11-10 (p)): the only categories g_mom's factors edit
# at test. These are the two any backtest read (T25, T26); epicardium and vascular endothelium have factors but no read.
G_MOM_EDIT_CATEGORIES = ("cm", "endocardium")
# --feng-groups (off by default): the cell groups Feng's direction is formed on. "authors" is Feng's own labels, every
# file as built. "none" is no Feng direction on any group: the entry with that one step taken out and nothing put in
# its place: no Feng file is read, and neither is any list that was fixed by counting Feng's cells (docs/PLAN.md
# decision log 2026-11-19 (a) to (e): the authors' labels were made with every stage pooled, the in-window ones
# included, which an organiser's condition of 2026-11-16 excludes). Only the captain's ruling puts it in an upload.
# "recomputed" is the first upload's design on the groups experiments/external/feng_relabel.py recomputes from the
# permitted cells' raw counts alone (register rows X002 and X010; the captain's "utilize the atlases", decision log
# 2026-11-19 (k)): it builds one file, the first upload's command, and only when the labeller's committed record routes
# to it. Only the captain's ruling puts it in an upload.
FENG_GROUPS = feng_dirs.FENG_GROUPS
FENG_GROUPS_DEFAULT = FENG_GROUPS[0]
# R6 of experiments/harness/TRIES.md (2026-11-13: the AI assistant's decision under the captain's delegation,
# docs/PLAN.md decision log 2026-11-13 (d)): the pool T1's cells are drawn from. "pooled": every E10.5 cell (the entry
# until then, and still T1's second upload). "lib3": E10.5's dissected-heart library alone, whose cell names end in
# this suffix, at that library's own chamber mix (tries T27 to T31). Feng's direction, R1, R2 and the fit of g_mom
# are untouched.
# R7 (the captain's ruling of 2026-11-13, docs/PLAN.md decision log 2026-11-13 (k); try T33): "mid", the base half-way
# between the two, for T1's second upload: every library-3 cell and, per chamber category, a seeded random fraction
# s / (1 + s) of the other libraries' cells (s: library 3's share of that category's cells in the file), so that
# library 3 supplies (s + 1) / 2 of each category, at that pool's own chamber mix.
T1_BASES = ("pooled", "lib3", "mid")
T1_DISSECTED_LIBRARY = "3"
T1_MID_STREAM = 21  # the mid base's draw of other-library cells: default_rng([seed, 21]), apart from every build stream
T1_MID_CM_SHARE = 0.60  # the one cardiomyocyte share the mid base is ruled for (docs/S2_RULES.md: T1's s2)
# b0_trend's earlier input stages for the heart (--trend; row H11, the captain's ruling of 2026-11-13, docs/PLAN.md
# decision log (w)): role -> (file, day as vec.harness.STAGE_DAYS gives it); with E9.5 ("start") and E10.5 ("ref")
HEART_TREND_STAGES = {role: (f, harness.STAGE_DAYS[f]) for role, f in
                      (("e8_75", "T2/E8.75.h5ad"), ("e8_5", "T2/E8.5.h5ad"), ("e8_25_late", "T2/E8.25_late.h5ad"))}
# The heart's scope bet (campaign/heart-iso-trim; no ruling; only the captain's ruling puts it in an upload): the ruled
# second file's two fixes (the copy's expression, the size step's raw mode) on the pre-registered trimmed scope, H2 as
# the half trim of docs/S2_RULES.md, sized by the size file's own logged raw-mode figure for H2. Scope is a lever of the
# post-exposure block (docs/DATA_POLICY.md 2026-11-04): the option is off unless all three flags are given.
HEART_ISO_TRIM = "heart --entry s2_trim --expression copy --z-mode raw"
# The full trim (campaign/heart-scope-options; no ruling; off unless all three flags are given, only the captain's
# ruling puts it in an upload): the same two fixes on E10.5 trimmed by vec.t2_build.radial_trim_mask at fraction 0
# (surrounding cells removed farthest in plane first until the in-plane scope is at most zero: the core and any
# surrounding cells inside its in-plane spread), sized to the raw mode's own core forecast
# (z_modes.raw.pred_rms_core of the size JSON: the cloud at the core's own scope, full-over-core ratio one).
HEART_CORE_TRIM = "heart --entry core_trim --expression copy --z-mode raw"
# Per scope option: its trim fraction, the size step's hypothesis, the size JSON field its size comes from
HEART_SCOPE_OPTIONS = {
    "s2_trim": {"frac": 0.5, "hypothesis": "H2", "size_source": "z_modes.raw.scope.H2.pred_rms_full",
                "rule": "surrounding cells removed farthest in plane from the core first, to half the in-plane scope"},
    "core_trim": {"frac": 0.0, "hypothesis": t2_build.SIZE_CORE_SCOPE, "size_source": "z_modes.raw.pred_rms_core",
                  "rule": "surrounding cells removed farthest in plane from the core first, to an in-plane scope of "
                          "at most zero (the core's own)"},
}
# --shape-interp (row E49 as drafted; no ruling): off by default; the fixed order of the steps when it is given
SHAPE_INTERP_OFF = "none"
SHAPE_INTERP_ORDER = ("floor_exact positions, the shape form, the size rule's RMS rescale, the duplicate jitter, the "
                      "cloud guard, then --row-order copy_parity if given (on the final coordinates)")
# --copy-rows (off by default; the investigator's Option A, the captain's order of 2026-11-21, docs/PLAN.md decision log
# 2026-11-21 (h); row E50 as specified): the floor copy's own rows, order, count and float32 coordinates with the
# exact-position file's expression rule (vec.t2_build.copy_rows_cloud)
COPY_ROWS_EXPRESSIONS = ("mixture", "series")
COPY_ROWS_N = "_copy_rows_n"  # the floor copy's count, handed to _hide_counts and never written
COPY_ROWS_COUNT_RULE = "departed: the floor copy's count (the captain's order of 2026-11-21)"
# --copy-rows-size (off by default; row E52 as specified, the captain's order of 2026-11-21, docs/PLAN.md decision log
# 2026-11-21 (v)): the copy-rows file (mixture only) rescaled about the copy's centroid to the pre-registered size
COPY_ROWS_SIZE_EXPRESSION = "mixture"
COPY_ROWS_SIZE_RULE = "applied: rescale_to_rms about the copy centroid to preregistered pred_rms_full (raw mode)"
COPY_ROWS_SIZE_RMS_TOL = 1e-6  # the RMS self-check, relative
COPY_ROWS_PARITY_FROM = 4  # the prefix parity self-check's smallest count (three rows have rank two: arbitrary)
COPY_ROWS_PARITY_SEED = 0  # its metric seed
HEART_HALVING = ("experiments/builds/halve.py <file> <half>: a uniform random half of the rows (halve.DEFAULT_SEED, "
                 "order kept), as for the ruled heart2 file; the halving is a separate step, not this build")


def _view(key: str, task: str, kind: str, setting: str | None, test_board: str, roles: dict[str, tuple[str, float]],
          reference_role: str, target_day: float) -> harness.TaskView:
    """A TaskView for the test geometry (the harness's view of a pseudo-task, with the test's stages)."""
    return harness.TaskView(key=key, task=task, kind=kind, setting=setting, test_board=test_board,
                            roles=tuple(roles), role_days={r: d for r, (_, d) in roles.items()},
                            role_stems={r: harness.stem(f) for r, (f, _) in roles.items()},
                            reference_role=reference_role, target_day=target_day)


def _board(which: str, override: str | None):
    boards = load_boards()
    key = override or TEST_BOARDS[which]
    if key not in boards:
        sys.exit(f"{key} is not in boards/: run the board intake first, or pass --board <validation board> for a dry "
                 "run")
    dry = key != TEST_BOARDS[which]
    return boards[key], dry


def _t2_count(board, pool_n: int) -> dict:
    """Approved rule 7 through its committed decision file (runbook section 5), checked against the count the harness's
    builder draws (harness.count_for on the candidate board); returns the decision's provenance, no counts."""
    n_rule = harness.count_for(harness._candidate_board(board), pool_n)
    n_dec, rec = t2_build.n_cells_from_decision(board, pool_n, COUNT_DECISION)
    if n_dec != n_rule:
        sys.exit("the committed count decision gives another count than the builder's count rule (values not shown); "
                 "build through n_cells_from_decision before any upload")
    d = rec["count_decision"]
    return {k: d.get(k) for k in ("decision_file", "decision_sha256", "kind", "validation_target", "override",
                                  "decision")}


def _hide_counts(rec: dict, board) -> dict:
    """On a test board, the record without cell counts: n_cells becomes "max_cells" (the count rule), and anything
    else refuses. T1's realised composition is left out (below)."""
    copy_n = rec.get(COPY_ROWS_N)  # --copy-rows only: the floor copy's own count, checked and then dropped
    if copy_n is not None:
        rec = {k: v for k, v in rec.items() if k != COPY_ROWS_N}
    if not board.counts_hidden:
        return rec
    if copy_n is not None:  # off by default: the copy-rows file holds the floor copy's count, not max_cells
        if rec["n_cells"] != copy_n:
            sys.exit("the built file does not hold the floor copy's count (values not shown)")
        out = dict(rec, n_cells="floor_copy")
        if isinstance(out.get("record"), dict):
            out["record"] = {k: v for k, v in out["record"].items() if k not in COUNT_FIELDS}
        return out
    if rec["n_cells"] != board.max_cells:
        sys.exit("the built file does not hold the board's max_cells, the approved count (values not shown)")
    out = dict(rec, n_cells="max_cells")
    # Stripped, not rounded (docs/PLAN.md decision log 2026-11-06 (g)): T1's realised composition is quota / n_cells per
    # category (t1_scope.composition of the drawn cells), so at full precision its denominators give the closed test
    # count. Rounding would not make it safe: the quotas are vec.core.quotas of the target shares, which the record
    # keeps at full precision ("shares"), so a rounded realised share still moves with n_cells wherever a target sits
    # within 1/n_cells of a rounding boundary. And it adds nothing: the quotas meet the targets to within 1/n_cells by
    # construction. The h5ad is untouched: this is the record written beside it.
    out.pop("realised", None)
    if isinstance(out.get("record"), dict):
        out["record"] = {k: v for k, v in out["record"].items() if k not in COUNT_FIELDS}
    return out


def _check_seed(seed: int) -> int:
    if seed < harness.BUILD_SEED_OFFSET:
        sys.exit(f"build seed {seed} is below {harness.BUILD_SEED_OFFSET} (runbook section 5: never seed 0)")
    return seed


def _t1_momentum(hold: bool, load=load_stage) -> dict:
    """build_submission's momentum arguments for T1. R1 (experiments/harness/TRIES.md; the captain's ruling of
    2026-11-05): with hold, cells without a Feng group keep their E10.5 expression, so no momentum is built and E9.5
    is not loaded; without it they get B0's damped E9.5 -> E10.5 momentum (M1 as built before R1)."""
    if hold:
        return {}
    return {"start": load(DATA / "T1/E9.5_RNA.h5ad", lineage=False), "damping": harness.B0_DAMPING,
            "start_source": "E9.5_RNA"}


def _feng_hygiene_on(flag: bool | None, hold_unbracketed: bool) -> bool:
    """R2's build flag (experiments/harness/TRIES.md R2). Since the build commit after N12 (docs/PLAN.md decision log
    2026-11-06 (k)), --feng-hygiene not given (None) follows --hold-unbracketed: on with R1, the default, and off with
    --no-hold-unbracketed, which keeps rebuilding the archived s1 of ae6d6b3 byte-identically. --no-feng-hygiene
    rebuilds the s1 of 88b057f (M1 with R1, before R2). An explicit --feng-hygiene with --no-hold-unbracketed is
    refused (R2 departs from M1 with R1)."""
    if flag and not hold_unbracketed:
        sys.exit("--feng-hygiene is refused with --no-hold-unbracketed: R2 departs from M1 with R1 "
                 "(experiments/harness/TRIES.md R2)")
    return bool(hold_unbracketed) if flag is None else bool(flag)


def _no_blood(direction_groups) -> None:
    """R2's guard, on group names only: --feng-hygiene is refused while Blood is among Feng's direction groups."""
    if "Blood" in set(direction_groups):
        sys.exit("--feng-hygiene is refused: Blood is among Feng's direction groups (experiments/harness/TRIES.md R2)")


def _feng_directions(genes, feng_groups: str = FENG_GROUPS_DEFAULT):
    """feng_dirs.feng_directions on the groups --feng-groups names. The default makes exactly the call the entries were
    built with (no keyword), so the default build is unchanged; "none" asks for no direction and reads no Feng file;
    "recomputed" asks for the permitted-only groups of experiments/external/feng_relabel.py (refused there unless that
    file is written through the labeller's gate, hashed in its manifest, usable, made from the yaml as it stands and
    its manifest committed as the labeller's RECORD), and is refused here unless the yaml's terminal rule routes that
    relabelling to recomputed."""
    if feng_groups == FENG_GROUPS_DEFAULT:
        return feng_dirs.feng_directions(genes)
    if feng_groups != "recomputed":
        return feng_dirs.feng_directions(genes, groups=feng_groups)
    import feng_relabel

    try:
        directions, ext = feng_dirs.feng_directions(genes, groups=feng_groups)
    except (ValueError, FileNotFoundError, StopIteration, feng_relabel.Refused) as e:
        sys.exit(f"--feng-groups recomputed is refused: {e}")
    if ext.get("route") != "recomputed":
        sys.exit(f"--feng-groups recomputed is refused: feng_relabel.yaml's terminal rule routes this relabelling to "
                 f"'{ext.get('route')}' (fallbacks {ext.get('fallbacks')}); the rule's file is --feng-groups "
                 f"{ext.get('route')}")
    return directions, ext


def _g_mom_covered(groups: dict, directions, feng_groups: str = FENG_GROUPS_DEFAULT) -> set[str]:
    """The labels g_mom's hook checks its edits against. Default: the labels the build moved along Feng's direction
    (as before). Under --feng-groups none: the labels the readout config's regexes give a group (`_labels_without_group`
    names the others), the same rule as the hook's own mask under none, so the check is a self-check there.
    vec.t1_state_transport.R1_FENG_GROUPS, a list fixed in N12 by counting Feng cells, is not read. Under recomputed:
    the default rule, the labels whose group has a recomputed direction after the fallback (`directions` holds only
    those groups), which is also the hook's own mask there (`_labels_without_direction`)."""
    if feng_groups == "none":
        return {lab for lab, group in groups.items() if not group.startswith(NO_GROUP)}
    return {lab for lab, group in groups.items() if group in directions}


def _labels_without_direction(labels, with_direction) -> set[str]:
    """Under --feng-groups recomputed, the labels g_mom holds back inside the categories it edits: those whose group
    (the readout config's regexes, `feng_dirs.label_groups`) has no recomputed direction after the fallback, labels
    with no group included. It takes the place of R1's mask (vec.t1_state_transport.r1_unbracketed_labels), whose list
    R1_FENG_GROUPS was fixed in N12 by counting cells under the authors' labels, so R1's mask in the build (cells of a
    group without a direction keep their E10.5 counts) and g_mom's mask follow the same groups."""
    groups = feng_dirs.label_groups(sorted({str(x) for x in labels}))
    return {lab for lab, group in groups.items() if group not in set(with_direction)}


NO_GROUP = "other:"  # what the readout config's harmonisation gives a label that none of its regexes matches


def _labels_without_group(labels) -> set[str]:
    """The labels the readout config's regexes map to no harmonised group (`feng_dirs.label_groups`: "other:..."). Under
    --feng-groups none this is the set g_mom holds back inside the categories it edits: the readout config's own
    rule, in place of R1's mask (vec.t1_state_transport.r1_unbracketed_labels), which reads R1_FENG_GROUPS, a list
    fixed in N12 by counting Feng cells. On every label configs/t1_chamber_classes.yaml can place in cm or
    endocardium the two masks agree: both hold back Early VIC and nothing else (tests/test_t1_feng_groups_none.py
    pins this on the committed configs)."""
    groups = feng_dirs.label_groups(sorted({str(x) for x in labels}))
    return {lab for lab, group in groups.items() if group.startswith(NO_GROUP)}


FENG_NONE_COMMAND = "t1 --g-mom --base lib3"  # the first upload's command, the one design --feng-groups none builds


def _feng_none_refusals(a) -> None:
    """--feng-groups none takes Feng's direction out and puts nothing in its place, and builds one file: the first
    upload's command (FENG_NONE_COMMAND). Refused before anything is read: without --g-mom, and with --base pooled or
    mid, --cm-share, --blood-share, --lib3-pair or --herd, which no row prices with it; then the options that act on or
    beside that direction: --no-hold-unbracketed (B0's momentum on every cell, a file no row prices), an explicit
    --feng-hygiene (R2 holds genes flat in a direction that is not there) and --strong-null (its prior acts on the
    cells Feng's direction moves)."""
    if (not a.g_mom or a.base != "lib3" or a.cm_share is not None or a.blood_share is not None or a.lib3_pair
            or a.herd):
        sys.exit(f"--feng-groups none builds one file: the first upload's command ({FENG_NONE_COMMAND}) with Feng's "
                 "direction taken out; it is refused with --base pooled or mid, --cm-share, --blood-share, "
                 "--lib3-pair, --herd or without --g-mom, which no row prices with it")
    if not a.hold_unbracketed:
        sys.exit("--feng-groups none is refused with --no-hold-unbracketed: with no Feng group every cell would get "
                 "B0's momentum, a file no row of experiments/harness/TRIES.md prices")
    if a.feng_hygiene:
        sys.exit("--feng-groups none is refused with --feng-hygiene: R2 holds genes flat in Feng's direction, and "
                 "this file has none")
    if a.strong_null:
        sys.exit("--feng-groups none is refused with --strong-null: its prior acts on the cells Feng's direction moves")
    _none_after_relabelling(getattr(a, "feng_none_ruling", None))


def _relabel_record() -> dict | None:
    """The permitted-only relabelling's record (experiments/external/feng_relabel.py RECORD), as in the tree or else
    as committed at HEAD; None while no relabelling exists. Its route is the yaml's terminal rule's outcome."""
    import feng_relabel

    raw = feng_relabel.RECORD.read_bytes() if feng_relabel.RECORD.exists() else feng_relabel.committed_at_head(
        feng_relabel.RECORD)
    return None if raw is None else json.loads(raw)


def _none_after_relabelling(ruling: str | None) -> None:
    """The final review's lesser point 1, for --feng-groups none. Before any relabelling exists, none builds as on main
    (and --feng-none-ruling is refused: it has nothing to override). After a relabelling whose recorded route is
    recomputed, none is refused unless --feng-none-ruling names a later docs/PLAN.md decision-log entry of the captain,
    naming N54 and containing the fixed phrase "sets the N54 route aside" (feng_relabel.NONE_RULING_PHRASE): not the
    ruling that authorised the run (the record's gate), dated on or after the record's commit, and not yet in
    docs/PLAN.md at that commit, which must be readable there (feng_relabel.later_ruling; the review of the port, lesser
    point 2; the final check of the fix round, items E and F). A switch to none made after N54's listing (which compares
    with the authors' groups) would otherwise be a choice by similarity to the authors' labels."""
    import feng_relabel

    rec = _relabel_record()
    if rec is None or rec.get("route") != "recomputed":
        if ruling is not None:
            sys.exit("--feng-none-ruling is refused: no relabelling routes to recomputed, so --feng-groups none needs "
                     "no ruling")
        return
    why = ("no --feng-none-ruling" if ruling is None else
           feng_relabel.later_ruling(ruling, feng_relabel.RECORD, feng_relabel.ROW,
                                     excluded=[(rec.get("gate") or {}).get("ruling")],
                                     phrase=feng_relabel.NONE_RULING_PHRASE))
    if why is not None:
        sys.exit("--feng-groups none is refused: the permitted-only relabelling's record routes T1's file to "
                 "recomputed (feng_relabel.yaml's terminal rule); only a later ruling of the captain, a docs/PLAN.md "
                 "decision-log entry naming N54 and containing \"sets the N54 route aside\", given as "
                 "--feng-none-ruling, other than the run's own ruling, dated on or after the record's commit and not "
                 f"yet in docs/PLAN.md then, sets that aside ({why})")


def _feng_recomputed_refusals(a) -> None:
    """--feng-groups recomputed builds one file: the first upload's command (FENG_NONE_COMMAND, t1 --g-mom --base lib3)
    on the recomputed groups, with R1 and R2 as in that upload. Refused before anything is read: without --g-mom, with
    --base pooled or mid, --cm-share, --blood-share, --lib3-pair, --herd, --no-hold-unbracketed, --no-feng-hygiene,
    --strong-null or --feng-none-ruling; then unless the labeller's record is committed at HEAD with the route
    recomputed (the yaml's terminal rule)."""
    import feng_relabel

    if (not a.g_mom or a.base != "lib3" or a.cm_share is not None or a.blood_share is not None or a.lib3_pair
            or a.herd or not a.hold_unbracketed or a.feng_hygiene is False or a.strong_null
            or getattr(a, "feng_none_ruling", None) is not None):
        sys.exit(f"--feng-groups recomputed builds one file: the first upload's command ({FENG_NONE_COMMAND}) on the "
                 "recomputed groups; it is refused with --base pooled or mid, --cm-share, --blood-share, --lib3-pair, "
                 "--herd, --no-hold-unbracketed, --no-feng-hygiene, --strong-null, --feng-none-ruling or without "
                 "--g-mom")
    head = feng_relabel.committed_at_head(feng_relabel.RECORD)
    if head is None:
        sys.exit("--feng-groups recomputed is refused: no permitted-only relabelling is committed "
                 f"({feng_relabel.RECORD.relative_to(REPO_ROOT)} at HEAD)")
    route = json.loads(head).get("route")
    if route != "recomputed":
        sys.exit(f"--feng-groups recomputed is refused: the committed relabelling's route is '{route}' "
                 "(feng_relabel.yaml's terminal rule); the rule's file is --feng-groups none")


# --released-bet (try T61 of experiments/harness/TRIES.md; off by default): a variant of T1's last file (the command
# FENG_NONE_COMMAND with --feng-groups none, row R8) that changes which cells are drawn. Both bets touch the
# post-exposure block (docs/DATA_POLICY.md 2026-11-04 (d): T1's composition, except by its legal route): each is
# refused unless --captain-ruling names a row of experiments/harness/TRIES.md that is the captain's ruling and names
# the arm. Read like for like: the trend and the axes run from E9.5 (embryo portions) to E10.5's portion libraries
# (every library but the dissected heart, T1_DISSECTED_LIBRARY), never to library 3, so no step is fitted on library
# 3's cells alone (row R6).
#   share_trend  library 3's own chamber mix x (E10.5 portions' share / E9.5's share) ** strength, renormalised
#                (vec.t1_bets.trend_shares; strength in t1_bets.LAMBDAS)
#   advance      library 3's cells drawn inside each category by quotas over five z-quantile bins, in proportion to
#                each bin's sum of exp(strength z), z the cell's projection on its category's E9.5 -> E10.5-portions
#                momentum axis (vec.t1_bets.advance_reweight; quotas, never importance weights: docs/DATA_POLICY.md
#                2026-11-04 (f)(2); strength in t1_bets.TAUS)
RELEASED_BETS = {"share_trend": ("l", t1_bets.LAMBDAS), "advance": ("t", t1_bets.TAUS)}
TRIES_FILE = Path(__file__).resolve().parents[2] / "experiments" / "harness" / "TRIES.md"


def _bet_arm(bet: str, strength: float) -> str:
    """The arm's name in row T61 (e.g. share_trend_l0.5, advance_t1)."""
    return f"{bet}_{RELEASED_BETS[bet][0]}{strength:g}"


def _captain_ruling_names(row: str, arm: str, path: Path | None = None) -> bool:
    """Whether experiments/harness/TRIES.md (TRIES_FILE) holds a ruling row `row` (its ID starts with R) whose first
    cell holds "(the captain's ruling" and whose text puts the arm in with the dedicated phrase "rules `<arm>` in". A
    row that merely mentions the arm, refuses it or awaits a ruling does not pass."""
    if not row.startswith("R"):
        return False
    for ln in (path or TRIES_FILE).read_text(encoding="utf-8").splitlines():
        if not ln.startswith(f"| {row} "):
            continue
        first = ln[2:].split("|", 1)[0]
        if "(the captain's ruling" in first and f"rules `{arm}` in" in ln:
            return True
    return False


def _bet_refusals(a) -> None:
    """--released-bet's refusals, before anything is read: only with --feng-groups none (whose own refusals then hold
    the rest of the command to the last file's); the strength must be one row T61 reads; and --captain-ruling must name
    the captain's ruling row for that arm. --captain-ruling or --bet-strength without --released-bet is refused too."""
    bet = getattr(a, "released_bet", None)
    if bet is None:
        ruling = getattr(a, "captain_ruling", None)
        if ruling is not None and (_ruled_fraction_ordered(a) or _ruled_damping_ordered(a)):  # rows R20, R22:
            return  # _t1_candidate refused every other use
        if ruling is not None or getattr(a, "bet_strength", None) is not None:
            sys.exit("--captain-ruling and --bet-strength go with --released-bet only")
        return
    if getattr(a, "feng_groups", FENG_GROUPS_DEFAULT) != "none":
        sys.exit("--released-bet is a variant of T1's last file: it needs --feng-groups none "
                 f"({FENG_NONE_COMMAND} --feng-groups none)")
    allowed = RELEASED_BETS[bet][1]
    if a.bet_strength is None or a.bet_strength not in allowed:
        sys.exit(f"--released-bet {bet} needs --bet-strength, one of {list(allowed)} (row T61's arms)")
    arm = _bet_arm(bet, a.bet_strength)
    if not a.captain_ruling or not _captain_ruling_names(a.captain_ruling, arm):
        sys.exit(f"--released-bet {bet} is refused: it changes T1's composition, a lever closed by the post-exposure "
                 "block, and only the captain's own ruling opens it; --captain-ruling must name a row R<n> of "
                 "experiments/harness/TRIES.md whose first cell holds \"(the captain's ruling\" and whose text says "
                 f"\"rules `{arm}` in\"")


# The T1 candidates of docs/PLAN.md decision log 2026-11-20 (i), each a single change to the recomputed last file
# (FENG_NONE_COMMAND with --feng-groups recomputed). Off by default (None: the option not given, today's path and record
# exactly). The first value of each tuple is the default and today's behaviour; the second is the candidate's.
G_MOM_DAMPINGS = (harness.B0_DAMPING, 1.0)  # A: g_mom's factors finished at 1.0 (practice arm gmom_s1, rows T49, T61)
FENG_FRACTION_MULTS = (1.0, 1.5)  # B: the direction at 1.5 x the calendar fraction (rows T50, T58, T63)
T1_CANDIDATE_COMMAND = f"{FENG_NONE_COMMAND} --feng-groups recomputed"
# Row R20 of experiments/harness/TRIES.md (the captain's order of 2026-11-21, docs/PLAN.md decision log 2026-11-21 (x),
# the amendment of docs/S2_RULES.md that names R20): one more file, T1_CANDIDATE_COMMAND with the direction at 2.0 x the
# calendar fraction and nothing else changed. Kept apart from FENG_FRACTION_MULTS, whose two values need no ruling: this
# one is taken only with --captain-ruling R20, and only while that row is the captain's order naming the option.
FENG_FRACTION_MULT_R20 = 2.0
T1_R20_ROW = "R20"
T1_R20_OPTION = "`--feng-fraction-mult 2.0`"
# Row R22 of experiments/harness/TRIES.md (the captain's order of 2026-11-22, docs/PLAN.md decision log 2026-11-22 (j),
# the amendment of docs/S2_RULES.md "a score-guided search of round steps may run"): the first search line, the next
# half steps of the same scale, 2.5 and 3.0 x the calendar fraction on T1_CANDIDATE_COMMAND with nothing else changed.
# Each is taken only with --captain-ruling R22, and only while that row is the captain's order naming the option.
FENG_FRACTION_MULTS_R22 = (2.5, 3.0)
T1_R22_ROW = "R22"
T1_R22_OPTIONS = {2.5: "`--feng-fraction-mult 2.5`", 3.0: "`--feng-fraction-mult 3.0`"}
# Every value taken only under a ruling row: value -> (the row, the option as the row's first cell names it, the
# record's order text). R20's entry is the 2.0 file exactly as before.
FENG_FRACTION_MULTS_RULED = {
    FENG_FRACTION_MULT_R20: (T1_R20_ROW, T1_R20_OPTION,
                             f"experiments/harness/TRIES.md row {T1_R20_ROW}, the captain's order of 2026-11-21 "
                             "(docs/PLAN.md decision log 2026-11-21 (x)): the direction at 2.0 x the calendar "
                             "fraction, nothing else changed"),
    **{v: (T1_R22_ROW, T1_R22_OPTIONS[v],
           f"experiments/harness/TRIES.md row {T1_R22_ROW}, the captain's order of 2026-11-22 (docs/PLAN.md decision "
           f"log 2026-11-22 (j)), the first search line: the direction at {v:.1f} x the calendar fraction, nothing "
           "else changed") for v in FENG_FRACTION_MULTS_R22},
}
# Row R22's second T1 search line (the lines added in docs/PLAN.md decision log 2026-11-22): g_mom's within-type
# damping at 1.5 and at 2.0 on T1_CANDIDATE_COMMAND, nothing else changed (fraction multiplier 1.0). Kept apart from
# G_MOM_DAMPINGS, whose two values need no ruling (1.0 is candidate A, as before): each is taken only with
# --captain-ruling R22, and only while that row is the captain's order naming the option.
G_MOM_DAMPINGS_R22 = (1.5, 2.0)
T1_R22_DAMPING_OPTIONS = {1.5: "`--g-mom-damping 1.5`", 2.0: "`--g-mom-damping 2.0`"}
G_MOM_DAMPING_ORDER_R22 = ("experiments/harness/TRIES.md row {row}, the captain's order of 2026-11-22 (docs/PLAN.md "
                           "decision log 2026-11-22), the second T1 search line: g_mom's within-type damping at "
                           "{v:.1f}, nothing else changed")


def _captain_order_names(row: str, option: str, path: Path | None = None) -> bool:
    """Whether experiments/harness/TRIES.md (TRIES_FILE) holds exactly one row `row`, a ruling row: its first cell
    starts with "<row> (the captain's order" and says "ruling, no run", and that first cell names `option` (in
    backticks) and "`--captain-ruling <row>`". A row that only mentions the option elsewhere does not pass."""
    lines = [ln for ln in (path or TRIES_FILE).read_text(encoding="utf-8").splitlines() if ln.startswith(f"| {row} ")]
    if len(lines) != 1:
        return False
    first = lines[0][2:].split("|", 1)[0]
    return (first.startswith(f"{row} (the captain's order") and "ruling, no run" in first and option in first
            and f"`--captain-ruling {row}`" in first)


def _ruled_fraction_ordered(a) -> bool:
    """The commands the ruling rows order, as _t1_candidate holds them: --feng-fraction-mult 2.0 with --captain-ruling
    R20, or 2.5 or 3.0 with --captain-ruling R22, no --released-bet or --bet-strength, and the value's row the
    captain's order naming its option."""
    ruled = FENG_FRACTION_MULTS_RULED.get(getattr(a, "feng_fraction_mult", None))
    return (ruled is not None and getattr(a, "captain_ruling", None) == ruled[0]
            and getattr(a, "released_bet", None) is None and getattr(a, "bet_strength", None) is None
            and _captain_order_names(ruled[0], ruled[1]))


def _ruled_damping_ordered(a) -> bool:
    """The commands row R22's second line orders, as _t1_candidate holds them: --g-mom-damping 1.5 or 2.0 with
    --captain-ruling R22, --feng-fraction-mult absent or 1.0, no --released-bet or --bet-strength, and row R22 the
    captain's order naming the option."""
    damping = getattr(a, "g_mom_damping", None)
    return (damping in T1_R22_DAMPING_OPTIONS and getattr(a, "captain_ruling", None) == T1_R22_ROW
            and getattr(a, "feng_fraction_mult", None) in (None, FENG_FRACTION_MULTS[0])
            and getattr(a, "released_bet", None) is None and getattr(a, "bet_strength", None) is None
            and _captain_order_names(T1_R22_ROW, T1_R22_DAMPING_OPTIONS[damping]))


def _t1_ruled_damping(a, damping: float) -> dict:
    """--g-mom-damping 1.5 or 2.0 (row R22's second T1 search line): refused before anything is read without
    --captain-ruling R22, with any other row, unless row R22 is the captain's order naming the option, with
    --released-bet or --bet-strength, with --feng-fraction-mult other than 1.0, and on any command but
    T1_CANDIDATE_COMMAND at the default seed on the test board. Returns the record's candidate values, naming the row."""
    row, option = T1_R22_ROW, T1_R22_DAMPING_OPTIONS[damping]
    mult, ruling = getattr(a, "feng_fraction_mult", None), getattr(a, "captain_ruling", None)
    if ruling is None:
        sys.exit(f"--g-mom-damping takes {G_MOM_DAMPINGS[0]:g} (the default) or {G_MOM_DAMPINGS[1]:g} only; "
                 f"{damping:g} is the captain's order of row {row} alone and needs --captain-ruling {row}")
    if ruling != row:
        sys.exit(f"--g-mom-damping {damping:g} is refused with --captain-ruling {ruling}: only row {row} orders it")
    if not _captain_order_names(row, option):
        sys.exit(f"--g-mom-damping {damping:g} is refused: experiments/harness/TRIES.md has no single row {row} whose "
                 f"first cell is the captain's order (\"ruling, no run\") naming {option} and "
                 f"`--captain-ruling {row}`")
    if getattr(a, "released_bet", None) is not None or getattr(a, "bet_strength", None) is not None:
        sys.exit(f"--g-mom-damping {damping:g} is refused with --released-bet and --bet-strength: row {row}'s line "
                 "is one file with one change")
    if mult is not None and mult != FENG_FRACTION_MULTS[0]:
        sys.exit(f"--g-mom-damping {damping:g} is refused with --feng-fraction-mult other than "
                 f"{FENG_FRACTION_MULTS[0]:g}: row {row}'s line is one file with one change")
    if (getattr(a, "feng_groups", FENG_GROUPS_DEFAULT) != "recomputed" or a.seed != harness.BUILD_SEED_OFFSET + 1
            or a.board is not None):
        sys.exit(f"--g-mom-damping {damping:g} changes T1's recomputed last file only: it needs the command "
                 f"{T1_CANDIDATE_COMMAND} at the default seed on the test board (row {row})")
    return {"g_mom_damping": float(damping), "feng_fraction_mult": FENG_FRACTION_MULTS[0], "captain_ruling": row,
            "order": G_MOM_DAMPING_ORDER_R22.format(row=row, v=float(damping))}


def _t1_ruled(a, mult: float) -> dict:
    """--feng-fraction-mult 2.0 (row R20's file) or 2.5 or 3.0 (row R22's first search line): refused before anything
    is read without --captain-ruling naming the value's row, with any other row, unless that row is the captain's order
    naming the option, with --released-bet or --bet-strength, with --g-mom-damping other than the default, and on any
    command but T1_CANDIDATE_COMMAND at the default seed on the test board. Returns the record's candidate values,
    naming the row. At 2.0 the messages and the record are row R20's as before."""
    row, option, order = FENG_FRACTION_MULTS_RULED[mult]
    damping, ruling = getattr(a, "g_mom_damping", None), getattr(a, "captain_ruling", None)
    if ruling is None:
        sys.exit(f"--feng-fraction-mult takes {FENG_FRACTION_MULTS[0]:g} (the default) or {FENG_FRACTION_MULTS[1]:g} "
                 f"only; {mult:g} is the captain's order of row {row} alone and needs --captain-ruling {row}")
    if ruling != row:
        sys.exit(f"--feng-fraction-mult {mult:g} is refused with --captain-ruling {ruling}: only row {row} orders it")
    if not _captain_order_names(row, option):
        sys.exit(f"--feng-fraction-mult {mult:g} is refused: experiments/harness/TRIES.md has no single row {row} "
                 f"whose first cell is the captain's order (\"ruling, no run\") naming {option} and "
                 f"`--captain-ruling {row}`")
    if getattr(a, "released_bet", None) is not None or getattr(a, "bet_strength", None) is not None:
        sys.exit(f"--feng-fraction-mult {mult:g} is refused with --released-bet and --bet-strength: row {row} is one "
                 "file with one change")
    if damping is not None and damping != G_MOM_DAMPINGS[0]:
        sys.exit(f"--feng-fraction-mult {mult:g} is refused with --g-mom-damping other than {G_MOM_DAMPINGS[0]:g}: "
                 f"row {row} is one file with one change")
    if (getattr(a, "feng_groups", FENG_GROUPS_DEFAULT) != "recomputed" or a.seed != harness.BUILD_SEED_OFFSET + 1
            or a.board is not None):
        sys.exit(f"--feng-fraction-mult {mult:g} changes T1's recomputed last file only: it needs the command "
                 f"{T1_CANDIDATE_COMMAND} at the default seed on the test board (row {row})")
    return {"g_mom_damping": G_MOM_DAMPINGS[0], "feng_fraction_mult": float(mult), "captain_ruling": row,
            "order": order}


def _t1_candidate(a) -> dict | None:
    """--g-mom-damping and --feng-fraction-mult: refused before anything is read with a value outside their two, on any
    command but T1_CANDIDATE_COMMAND (--feng-groups recomputed, whose own refusals then hold every other option to the
    first upload's), with a seed other than the default and on a dry run. Returns None when neither is given (today's
    file and record), else the two values as the record names them. --feng-fraction-mult 2.0 is row R20's file alone,
    2.5 and 3.0 row R22's line alone (_t1_ruled); --g-mom-damping 1.5 and 2.0 row R22's second T1 line alone
    (_t1_ruled_damping); with the flag --captain-ruling, any other value is refused here or by _bet_refusals."""
    damping, mult = getattr(a, "g_mom_damping", None), getattr(a, "feng_fraction_mult", None)
    if damping is None and mult is None:
        return None
    if mult in FENG_FRACTION_MULTS_RULED:
        return _t1_ruled(a, mult)
    if damping in T1_R22_DAMPING_OPTIONS:
        return _t1_ruled_damping(a, damping)
    if damping is not None and damping not in G_MOM_DAMPINGS:
        sys.exit(f"--g-mom-damping takes {G_MOM_DAMPINGS[0]:g} (the default) or {G_MOM_DAMPINGS[1]:g} only "
                 f"({G_MOM_DAMPINGS_R22[0]:g} and {G_MOM_DAMPINGS_R22[1]:g} only with --captain-ruling {T1_R22_ROW})")
    if mult is not None and mult not in FENG_FRACTION_MULTS:
        sys.exit(f"--feng-fraction-mult takes {FENG_FRACTION_MULTS[0]:g} (the default) or {FENG_FRACTION_MULTS[1]:g} "
                 f"only ({FENG_FRACTION_MULT_R20:g} only with --captain-ruling {T1_R20_ROW}; "
                 f"{FENG_FRACTION_MULTS_R22[0]:g} and {FENG_FRACTION_MULTS_R22[1]:g} only with --captain-ruling "
                 f"{T1_R22_ROW})")
    if (getattr(a, "feng_groups", FENG_GROUPS_DEFAULT) != "recomputed" or a.seed != harness.BUILD_SEED_OFFSET + 1
            or a.board is not None):
        sys.exit("--g-mom-damping and --feng-fraction-mult change T1's recomputed last file only: they need the "
                 f"command {T1_CANDIDATE_COMMAND} at the default seed on the test board (docs/PLAN.md decision log "
                 "2026-11-20 (i))")
    return {"g_mom_damping": G_MOM_DAMPINGS[0] if damping is None else float(damping),
            "feng_fraction_mult": FENG_FRACTION_MULTS[0] if mult is None else float(mult),
            "order": "docs/PLAN.md decision log 2026-11-20 (i): A = damping 1.0, B = 1.5 x the calendar fraction, C = "
                     "both; no ruling"}


def _feng_group_names(feng_groups: str = FENG_GROUPS_DEFAULT) -> list[str]:
    """The names of Feng's direction groups from Feng's statistics alone, before any T1 file is read. The groups
    feng_dirs.feng_directions returns do not depend on the gene list it is given, so it is called on one placeholder
    name (an empty list fails in dirlib.sym_index, whose empty index is a float array) and only the names are kept.
    The default makes exactly the call as before; recomputed asks for the recomputed groups' names."""
    if feng_groups == FENG_GROUPS_DEFAULT:
        return sorted(feng_dirs.feng_directions([""])[0])
    return sorted(_feng_directions([""], feng_groups)[0])


def _t1_union(hygiene_on: bool, feng_groups: str = FENG_GROUPS_DEFAULT) -> dict | None:
    """R2's union for --feng-hygiene, after R2's two refusals that need no T1 file, both made before anything else is
    read (R2: no E10.5 file is read before the build): hygiene.load_union refuses while configs/t1_hygiene.yaml pins no
    union_sha256 (until the build commit after N12) or pins another file's; then the Blood guard on the names of Feng's
    direction groups (on the groups --feng-groups names). Off: None, and nothing is read."""
    if not hygiene_on:
        return None
    try:
        union = hygiene.load_union()
    except hygiene.Refused as e:
        sys.exit(f"--feng-hygiene is refused: {e}")
    _no_blood(_feng_group_names() if feng_groups == FENG_GROUPS_DEFAULT else _feng_group_names(feng_groups))
    return union


def _t1_down_prior(on: bool, hold_unbracketed: bool, hygiene_on: bool) -> t1_scope.DownPrior | None:
    """The option of T12a and T12b (experiments/harness/TRIES.md; the captain's rulings of 2026-11-10 (a) and (d)): the
    strong-null down-tail prior after M1's steps. Off (the default): None, and nothing is read. On: refused unless the
    build is M1 with R1 and R2 (T12 departs from it) and unless configs/t1_strong_null.yaml says approved
    (`strong_null.check_approved`), both before anything else is read. The genes are ordered by the reference's level
    alone (`DownPriorTable.level_only`, the selection the backtests measured); the prior's cells are the ones M1 moves
    along Feng's direction (`covered_labels` None); `t1_scope.apply_down_prior` holds R2's and the pre-registered genes
    exactly as M1's own step does."""
    if not on:
        return None
    if not (hold_unbracketed and hygiene_on):
        sys.exit("--strong-null is refused without R1 and R2 (--no-hold-unbracketed, --no-feng-hygiene): T12 departs "
                 "from M1 with R1 and R2 (experiments/harness/TRIES.md T12a, T12b)")
    try:
        strong_null.check_approved()
    except strong_null.Refused as e:
        sys.exit(f"--strong-null is refused: {e}")
    return t1_scope.DownPrior(table=t1_scope.DownPriorTable.level_only())


def _t1_held(union: dict | None, direction_groups, build_genes) -> tuple[list[str], dict | None]:
    """The held-flat genes build_t1 passes to build_submission, and the build record's feng_hygiene field. Off (no
    union): the pre-registered list exactly as before R2, so the file is byte-identical, and None. On (R2): that list
    joined with the union (_t1_union), after the Blood guard again on the groups the build uses; the field records the
    union's file, sha256 and gene count, the genes per module, the commit and the sources' sha256 (docs/PLAN.md 4.1,
    item 4), and, by name, the union genes absent from the build's genes (N12 mapped the union on E9.5's names; a name
    E10.5 lacks cannot be held and changes nothing)."""
    held = t1_scope.held_flat_genes()
    if union is None:
        return held, None
    _no_blood(direction_groups)
    present = {str(g) for g in build_genes}
    absent = [g for g in union["genes"] if g not in present]
    commit, _ = provenance.git_state(provenance.REPO_ROOT)
    return sorted(set(held) | set(union["genes"])), {
        "file": union["file"], "sha256": union["sha256"], "commit": commit, "genes": len(union["genes"]),
        "union_genes_absent_from_build": {"count": len(absent), "genes": absent}, "modules": union["modules"],
        "sources": union["sources"]}


def _t1_g_mom(ref, categories: list[str], seed: int, covered_labels: set[str], feng_direction: bool = True,
              with_direction=None, damping: float | None = None):
    """The test build's g_mom step (--g-mom, default off; T1 try T11 of experiments/harness/TRIES.md; only the captain's
    ruling puts it in an upload): E9.5 is loaded here and nowhere else, and the factors are fitted exactly as the
    candidate's (`vec.harness.g_mom_pair_factors`: start E9.5, reference E10.5, the chamber categories of E10.5's own
    mix, metacells of the E10.5 cells, damped by B0_DAMPING, zeroed and clipped). Returns the `post_cells` hook of
    `t1_scope.build_submission` and a record (categories only, no cell number) that the hook completes when it runs.
    The hook is the build's last count edit, after R1's hold and after Feng's direction (the order of the backtests:
    g_mom was judged on top of B0's per-type shift, and at test Feng's direction takes that place). It edits the cells
    of `vec.t1_state_transport.apply_state_factors`: a cell's source cell must sit in a fitted metacell; cells of R1's
    held labels (no Feng group: `r1_unbracketed_labels`, the candidate's own mask) keep the counts they have, so R1
    stands; R3b: the fit is the candidate's, over all `categories`, and only cells of G_MOM_EDIT_CATEGORIES are edited,
    the mask added to R1's (a category the fit skipped stops the build); the pre-registered held-flat genes are the
    candidate's own list (R2's union is held in Feng's direction only, by its ruling; if the captain wants R2 held
    in this step too, that is a new variant with its own row, since
    T11 was judged without it); each edited cell keeps its expected depth (`apply_factors`, keep_depth), the T1 depth
    lever stays closed; the cells, their number, shares and scope are untouched. The stream is [seed, 5], the
    candidate's. covered_labels: the labels the build moved along Feng's direction (`label_groups` mapped to a group
    with a direction); the hook checks that every label it edits is one of them, so R1's mask and the build's own mask
    of Feng-covered labels cannot disagree silently. feng_direction False (--feng-groups none): the same fit, stream
    and edit; the held labels are `_labels_without_group`'s, the readout config's own rule, so R1_FENG_GROUPS is not
    read (the same cells on every label configs/t1_chamber_classes.yaml places in those categories); the record
    says so (applied last, with no Feng direction in the file), and covered_labels is `_g_mom_covered`'s, so the
    check is a self-check. with_direction (--feng-groups recomputed: the groups with a recomputed direction after
    feng_relabel.yaml's fallback): the same fit, stream and edit; the held labels are `_labels_without_direction`'s,
    the labels whose group has no recomputed direction, so g_mom's mask and R1's hold in the build follow the same
    groups and R1_FENG_GROUPS (fixed under the authors' labels) is not read; the record says so. damping (candidate A,
    --g-mom-damping; None: today's factors and record): the same raw fit finished by
    `vec.t1_state_transport.finish_momentum` at that damping, as practice arm gmom_s1 (experiments/t1_m1/bets_read.py
    gmom_fits) finishes it; at B0_DAMPING that is the call g_mom_pair_factors makes itself."""
    from vec import t1_state_transport as ts

    start = load_stage(DATA / "T1/E9.5_RNA.h5ad", lineage=False)
    factors, raw = harness.g_mom_pair_factors(ref, "E10.5_RNA", start, "E9.5_RNA", categories=categories)
    if damping is not None and damping != harness.B0_DAMPING:  # candidate A: the same fit, finished at `damping`
        factors = ts.finish_momentum(raw, damping)
    missing = sorted(set(G_MOM_EDIT_CATEGORIES) - set(factors.categories))
    if missing:
        raise ValueError(f"--g-mom has no factors for a category R3b edits: {missing}")
    held_flat = t1_scope.held_flat_genes()
    record = {"rule": "try T11's within-type factors at test: OT within each chamber category from E9.5 to E10.5, "
                      "metacells of the E10.5 cells, the per-category mean removed, damped by B0_DAMPING, zeroed and "
                      "clipped; applied last (after Feng's direction); R1's held labels unedited; keep_depth; R3b: "
                      "only the categories in edit_categories are edited; the factors are not scaled to the two-day "
                      "step",
              "edit_categories": list(G_MOM_EDIT_CATEGORIES),
              "categories_with_factors": sorted(factors.categories),
              "categories_without_factors": sorted(factors.skipped), "stream": [int(seed), 5],
              "held_flat": "the candidate's pre-registered list (R2's union is held in Feng's direction only)"}
    if not feng_direction:  # --feng-groups none: written only then, so every other file's record is as before
        record["rule"] = record["rule"].replace("applied last (after Feng's direction)",
                                                "applied last (no Feng direction in this file)")
        record["rule"] = record["rule"].replace("R1's held labels unedited",
                                                "labels the readout config's regexes give no group unedited")
        record["held_flat"] = "the candidate's pre-registered list"
        record["mask"] = ("the cells of edit_categories whose label the readout config's regexes give a group (the "
                          "readout config's own rule; the list of groups with a Feng direction is not read)")
    if with_direction is not None:  # --feng-groups recomputed: written only then
        record["rule"] = record["rule"].replace("applied last (after Feng's direction)",
                                                "applied last (after Feng's direction on the recomputed groups)")
        record["rule"] = record["rule"].replace("R1's held labels unedited",
                                                "labels of groups without a recomputed direction unedited")
        record["mask"] = ("the cells of edit_categories whose label's group has a recomputed direction after the "
                          f"fallback ({sorted(with_direction)}); R1_FENG_GROUPS is not read")
    if damping is not None:  # written only with --g-mom-damping, so every other file's record is as before
        record["damping"] = float(damping)
        if damping in T1_R22_DAMPING_OPTIONS:  # row R22's second T1 search line; candidate A's text below as before
            record["rule"] = record["rule"].replace("damped by B0_DAMPING",
                                                    f"damped by {float(damping):g} in place of B0_DAMPING (row "
                                                    f"{T1_R22_ROW}'s second T1 search line, docs/PLAN.md decision log "
                                                    "2026-11-22)")
        elif damping != harness.B0_DAMPING:
            record["rule"] = record["rule"].replace("damped by B0_DAMPING",
                                                    f"damped by {float(damping):g} in place of B0_DAMPING (candidate "
                                                    "A, docs/PLAN.md decision log 2026-11-20 (i))")

    def hook(pop):
        labels = pop.obs[t1_scope.LABEL_KEY].astype(str).to_numpy()
        cats = pop.obs[t1_scope.CATEGORY_KEY].astype(str).to_numpy()
        if with_direction is not None:
            held_labels = _labels_without_direction(np.unique(labels), with_direction)
        else:
            held_labels = (ts.r1_unbracketed_labels(np.unique(labels)) if feng_direction
                           else _labels_without_group(np.unique(labels)))
        no_edit = np.isin(labels, sorted(held_labels)) | ~np.isin(cats, G_MOM_EDIT_CATEGORIES)
        edited = set(labels[~no_edit & np.isin(cats, list(factors.categories))])
        if not edited <= set(covered_labels):
            extra = sorted(edited - set(covered_labels))
            what = "Feng's direction does not cover" if feng_direction else "the readout config gives no group"
            raise ValueError(f"--g-mom would edit labels {what}: {extra}")
        moved, rec = ts.apply_state_factors(pop, factors, held_flat, [seed, 5], rescale=False, no_edit=no_edit)
        record["median_depth_scale"] = rec["median_depth_scale"]  # no cell number
        return moved

    return hook, record


def _t1_mid_rows(ref, in3: np.ndarray, seed: int) -> tuple[np.ndarray, dict]:
    """The mid base's pool (row R7; try T33's arm as a builder gets it): every library-3 cell of the reference, and in
    each chamber category a random fraction s / (1 + s) of the other libraries' cells, s being library 3's share of that
    category's cells in the file. Categories in sorted order, one generator, default_rng([seed, T1_MID_STREAM]). The
    other libraries' cells outside the chamber categories are left out (no draw takes them). Returns the kept rows in
    file order and a record of shares only (no cell number)."""
    spec = t1_scope.load_spec()
    b = harness._shallow(ref)
    t1_scope.classify(b, source="E10.5_RNA", spec=spec)
    cat = np.asarray(b.obs[t1_scope.CATEGORY_KEY]).astype(str)
    rng = np.random.default_rng([seed, T1_MID_STREAM])
    keep = [np.flatnonzero(in3)]
    record = {}
    for c in sorted(spec.in_scope("chamber")):
        rows3, rows_o = np.flatnonzero(in3 & (cat == c)), np.flatnonzero(~in3 & (cat == c))
        if len(rows3) + len(rows_o) == 0:
            continue
        if len(rows3) == 0:
            sys.exit(f"--base mid: library 3 holds no cell of the chamber category {c}, so the mid base is undefined")
        s = len(rows3) / (len(rows3) + len(rows_o))
        k = int(round(len(rows_o) * s / (1.0 + s)))
        kept = np.sort(rng.choice(rows_o, size=k, replace=False)) if k else np.array([], dtype=np.int64)
        keep.append(kept)
        record[c] = {"library3_share_in_file": s, "other_fraction_kept": k / len(rows_o) if len(rows_o) else 0.0,
                     "library3_share_in_pool": len(rows3) / (len(rows3) + k)}
    return np.sort(np.concatenate(keep)), record


def _t1_trend(ref, portions: np.ndarray, shares: dict, lam: float) -> tuple[dict, dict]:
    """--released-bet share_trend: the pool's own chamber shares moved along the like-for-like trend from E9.5 to
    E10.5's portion libraries (`portions`, a mask of the reference's rows), `t1_bets.trend_shares` at `lam`. E9.5 is
    loaded here for its labels' counts. Returns the shares and a record (shares and ratios, no count)."""
    start = load_stage(DATA / "T1/E9.5_RNA.h5ad", lineage=False)
    later = t1_bets.chamber_counts(ref[np.flatnonzero(portions)], "E10.5_RNA")
    new, rec = t1_bets.trend_shares(shares, later, t1_bets.chamber_counts(start, "E9.5_RNA"), lam)
    return t1_scope.load_shares(new), {"rule": "library 3's own chamber mix x (E10.5 portion libraries' share / E9.5's "
                                               "share) ** lam, renormalised; ratio 1 below the trend's min_cells",
                                       **rec}


def _t1_advance(ref, pool, portions: np.ndarray, shares: dict, n: int, seed: int, held: list[str], board, hook,
                tau: float) -> tuple:
    """--released-bet advance: the pool's cells drawn at `shares`, each category's quota split over z-quantile bins in
    proportion to exp(tau z) (`t1_bets.advance_reweight`), z on each category's momentum axis from E9.5 to E10.5's
    portion libraries (`portions`), the held genes at 0; then g_mom's hook and the board's genes
    (`t1_bets.assemble`), as `t1_scope.build_submission` does on the no-shift path (at tau 0 the same file, tested).
    Returns the submission, a report with the realised composition (stripped from the record on a test board), a
    record and the drawn population (before the hook), for its median library size."""
    if n < board.min_cells or (board.max_cells is not None and n > board.max_cells):
        raise ValueError("n_cells is outside the board's limits")
    spec = t1_scope.load_spec()
    start = harness._shallow(load_stage(DATA / "T1/E9.5_RNA.h5ad", lineage=False))
    later = ref[np.flatnonzero(portions)].copy()
    t1_scope.classify(start, source="E9.5_RNA", spec=spec)
    t1_scope.classify(later, source="E10.5_RNA", spec=spec)
    t1_scope.classify(pool, source="E10.5_RNA", spec=spec)  # in place, as build_submission classifies its ref
    axes, axes_rec = t1_bets.momentum_axes(start, later, [c for c, s in shares.items() if s > 0], held,
                                           start_fallback=spec.fallbacks("E9.5_RNA"),
                                           later_fallback=spec.fallbacks("E10.5_RNA"))
    z = t1_bets.advance_scores(pool, axes)
    pop, _ = t1_bets.advance_reweight(pool, shares, n, seed, z, tau, fallback=spec.fallbacks("E10.5_RNA"), spec=spec)
    sub = t1_bets.assemble(pop, board, hook)
    return sub, {"realised": t1_scope.composition(pop)}, {
        "rule": "library 3's cells drawn inside each chamber category by quotas over five z-quantile bins in "
                "proportion to each bin's sum of exp(tau z), z the cell's projection on its category's E9.5 -> "
                "E10.5-portion-libraries momentum axis, standardised within the category; category quotas at library "
                "3's own chamber mix; then g_mom", "tau": float(tau), "axes": axes_rec}, pop


def _bet_median_library(pool, base_shares: dict, bet_pop, n: int, seed: int) -> dict:
    """docs/DATA_POLICY.md 2026-11-04 (f)(3) for --released-bet: the bet's drawn cells' median library size over the
    none file's draw at the same seed (`reweight` on the same pool at the pool's own chamber mix, as build_submission
    draws it on that file's path). A ratio of two medians of released cells, no count."""
    spec = t1_scope.load_spec()
    t1_scope.classify(pool, source="E10.5_RNA", spec=spec)  # in place and idempotent, as build_submission does
    base, _ = t1_scope.reweight(pool, base_shares, n, seed, fallback=spec.fallbacks("E10.5_RNA"), spec=spec)
    return {"median_library_over_none_draw": t1_bets.median_library(bet_pop) / t1_bets.median_library(base),
            "against": "the none file's draw at the same seed (library 3's cells at their own chamber mix)"}


def build_t1(a) -> dict:
    feng_groups = getattr(a, "feng_groups", FENG_GROUPS_DEFAULT)
    candidate = _t1_candidate(a)  # None unless --g-mom-damping or --feng-fraction-mult is given; refusals read nothing
    _bet_refusals(a)  # nothing at all unless --released-bet (or one of its two companions) is given
    if feng_groups == "none":
        _feng_none_refusals(a)
    elif feng_groups == "recomputed":
        _feng_recomputed_refusals(a)
    elif getattr(a, "feng_none_ruling", None) is not None:
        sys.exit("--feng-none-ruling is refused without --feng-groups none")
    if a.g_mom and a.strong_null:
        sys.exit("--g-mom is refused with --strong-null: the two count edits have no ruled order and were never tested "
                 "together; a new row in experiments/harness/TRIES.md and the captain's ruling fix it first")
    if a.g_mom and not a.hold_unbracketed:
        sys.exit("--g-mom is refused with --no-hold-unbracketed: g_mom's factors at test add to M1 with R1 and R2 as "
                 "they are (experiments/harness/TRIES.md, try T11)")
    if a.lib3_pair and (a.base != "lib3" or a.cm_share is None or a.blood_share is not None):
        sys.exit("--lib3-pair names the second file of the pair try T32 priced (library 3's cells with the "
                 "cardiomyocyte share moved): it needs --base lib3 with --cm-share, and no --blood-share")
    if a.base == "lib3" and (a.cm_share is not None or a.blood_share is not None) and not a.lib3_pair:
        sys.exit("--base lib3 is refused with --cm-share or --blood-share: the pair priced in try T30 keeps the "
                 "second upload on the pooled base (experiments/harness/TRIES.md, row R6). --lib3-pair builds try "
                 "T32's second file instead; it is no entry, and it is uploaded only by the captain's own ruling")
    if a.base == "lib3" and (not a.hold_unbracketed or a.strong_null or a.herd):
        sys.exit("--base lib3 is refused with --no-hold-unbracketed, --strong-null or --herd: each would fit or pull a "
                 "step on library 3's cells alone, which no try has read (experiments/harness/TRIES.md, row R6)")
    if a.base == "mid" and (a.cm_share != T1_MID_CM_SHARE or a.blood_share is not None or a.lib3_pair
                            or not a.hold_unbracketed or a.strong_null or a.herd):
        sys.exit("--base mid builds one file, T1's second upload by the captain's ruling of 2026-11-13 (row R7 of "
                 f"experiments/harness/TRIES.md; try T33): it needs --cm-share {T1_MID_CM_SHARE:.2f}, and it is "
                 "refused with another share, --blood-share, --lib3-pair, --no-hold-unbracketed, --strong-null or "
                 "--herd, which no try has read on that base")
    # R2's refusals come before anything is read. Under none R2 is off: it holds genes flat in Feng's direction, and
    # there is none (an explicit --feng-hygiene, the one case _feng_hygiene_on could refuse, was refused above)
    hygiene_on = False if feng_groups == "none" else _feng_hygiene_on(a.feng_hygiene, a.hold_unbracketed)
    prior = _t1_down_prior(a.strong_null, a.hold_unbracketed, hygiene_on)  # T12's, likewise
    union = _t1_union(hygiene_on) if feng_groups == FENG_GROUPS_DEFAULT else _t1_union(hygiene_on, feng_groups)
    board, dry = _board("t1", a.board)
    ref = load_stage(DATA / "T1/E10.5_RNA.h5ad", lineage=False)
    momentum = _t1_momentum(a.hold_unbracketed)
    shares = harness.t1_chamber_shares(harness._shallow(ref), "E10.5_RNA")
    fit_categories = [c for c, s in shares.items() if s > 0]  # g_mom's fit, as the candidate's: E10.5's own mix
    pool = ref
    mid_record = None
    if a.base in ("lib3", "mid"):
        names = [str(x) for x in ref.obs_names]
        library = np.array([x.rsplit("_", 1)[-1] if "_" in x else "" for x in names])
        kinds = set(library.tolist())
        if not all(k.isdigit() for k in kinds) or T1_DISSECTED_LIBRARY not in kinds or len(kinds) < 2:
            sys.exit(f"--base {a.base}: the reference's cell names do not separate library 3 from the others (every "
                     "name must end in _<library number>, and library 3 must be one of at least two)")
        in3 = library == T1_DISSECTED_LIBRARY
        if a.base == "lib3":
            pool = ref[np.flatnonzero(in3)].copy()
        else:
            rows, mid_record = _t1_mid_rows(ref, in3, a.seed)
            pool = ref[rows].copy()
        shares = harness.t1_chamber_shares(harness._shallow(pool), "E10.5_RNA")  # the pool's own chamber mix
    if a.cm_share is not None:
        shares = t1_scope.with_share(shares, "cm", a.cm_share)
    if a.blood_share is not None:  # blood's share of all chamber cells, the other categories in proportion
        shares = t1_scope.with_share(shares, "blood", a.blood_share, within=list(shares))
    directions, ext = _feng_directions(list(ref.var_names), feng_groups)
    held, feng_hygiene = _t1_held(union, directions, ref.var_names)
    groups = feng_dirs.label_groups(ref.obs["celltype"].astype(str).unique())
    calendar = feng_dirs.bracket_fraction(10.5, 12.5)
    # candidate B (--feng-fraction-mult): the multiplier times the calendar fraction, as double_step_read.py forms its
    # doubled arm (MULTIPLIER x fraction); without the option the calendar fraction itself, today's value
    fraction = calendar if candidate is None else candidate["feng_fraction_mult"] * calendar
    n = harness.count_for(board, ref.n_obs)
    extra = {} if prior is None else {"down_prior": prior}
    within_type, g_mom_record = ({}, None)
    if a.g_mom:
        covered = _g_mom_covered(groups, directions, feng_groups)  # the build's Feng-covered labels (none: self-check)
        if feng_groups == "recomputed":  # g_mom's mask follows the groups with a recomputed direction
            hook, g_mom_record = _t1_g_mom(ref, fit_categories, a.seed, covered, with_direction=sorted(directions),
                                           damping=None if candidate is None else candidate["g_mom_damping"])
        else:
            hook, g_mom_record = _t1_g_mom(ref, fit_categories, a.seed, covered,
                                           feng_direction=feng_groups != "none")
        within_type = {"post_cells": hook}
    bet = getattr(a, "released_bet", None)
    bet_record = bet_pop = None
    own_shares = shares
    if bet == "share_trend":  # refused above unless --feng-groups none, so the pool is library 3's (row R8's command)
        shares, bet_record = _t1_trend(ref, library != T1_DISSECTED_LIBRARY, shares, a.bet_strength)
    if bet == "advance":
        sub, rep, bet_record, bet_pop = _t1_advance(ref, pool, library != T1_DISSECTED_LIBRARY, shares, n, a.seed,
                                                    held, board, within_type["post_cells"], a.bet_strength)
    else:
        sub, rep = t1_scope.build_submission(pool, board, shares, n, a.seed, held_flat=held,
                                             ref_source="E10.5_RNA", directions=directions,
                                             direction_fraction=fraction, direction_groups=groups, herd=a.herd,
                                             **momentum, **extra, **within_type)
    if bet == "share_trend":  # the draw build_submission made at the trended shares (the same call), for (f)(3)
        spec = t1_scope.load_spec()
        bet_pop = t1_scope.reweight(pool, shares, n, a.seed, fallback=spec.fallbacks("E10.5_RNA"), spec=spec)[0]
    if bet_record is not None:
        bet_record = {**bet_record, **_bet_median_library(pool, own_shares, bet_pop, n, a.seed)}
    strong = None
    if prior is not None:  # T12's part (d): the file is written only if the target-free checks pass
        strong = {"selection": "the reference's level alone (vec.t1_scope.DownPriorTable.level_only), no table",
                  "constants": {"budget": prior.budget, "band": prior.band, "f_min": prior.f_min, "tol": prior.tol,
                                "mmd_bound": prior.mmd_bound},
                  **rep["down_prior"]}
        if not strong["checks"]["passed"]:
            sys.exit(f"the strong-null prior's target-free checks failed: {strong['checks']}")
    check = t1_scope.write_checked(sub, a.out, board)
    rec = {"entry": "M1", "board": board.key, "dry_run": dry, "n_cells": n, "seed": a.seed, "fraction": fraction,
           "herd": bool(a.herd), "hold_unbracketed": bool(a.hold_unbracketed), "feng_hygiene": feng_hygiene,
           "base": a.base,
           "cm_share_of_heart": t1_scope.share_of(shares, "cm", t1_scope.load_spec().heart),
           "blood_share_of_chamber": float(shares.get("blood", 0.0)),
           "shares": shares, "external": ext, "check_ok": check.ok, "realised": rep["realised"]}
    if a.lib3_pair:  # written only for that file, so every other file's record is as before
        rec["lib3_pair"] = True
    if feng_groups == "none":  # likewise: no direction step ran, so it has no fraction; no Feng file was read
        rec["fraction"] = None
        rec["feng_groups"] = {"groups": feng_groups, "file": ext["file"], "sha256": ext["sha256"],
                              "direction_step": "removed: no cell is moved along any direction; \"entry\": \"M1\" "
                                                "names the design this file comes from"}
        if getattr(a, "feng_none_ruling", None) is not None:  # only after a relabelling routed to recomputed
            rec["feng_groups"]["ruling_over_the_relabelling_route"] = a.feng_none_ruling
    if feng_groups == "recomputed":  # likewise: written only then; "external" names the same file
        rec["feng_groups"] = {"groups": feng_groups, "file": ext["file"], "sha256": ext["sha256"],
                              "relabel_config_sha256": ext["relabel_config_sha256"], "route": ext["route"],
                              "with_direction": ext["with_direction"], "fallbacks": ext["fallbacks"],
                              "relabel_gate": ext["relabel_gate"],
                              "groups_without_direction": "no cell of theirs moves: R1 holds them at their E10.5 "
                                                          "counts and g_mom's mask holds them back"}
    if candidate is not None:  # written only with --g-mom-damping or --feng-fraction-mult ("fraction" is the applied)
        rec["t1_candidate"] = {**candidate, "calendar_fraction": calendar, "applied_fraction": fraction,
                               "command": T1_CANDIDATE_COMMAND}
    if mid_record is not None:  # likewise: shares only, no cell number
        rec["mid_base"] = {"stream": [int(a.seed), T1_MID_STREAM], "categories": mid_record}
    if a.g_mom:
        rec["g_mom"] = g_mom_record
    if bet_record is not None:  # written only with --released-bet, so every other file's record is as before
        rec["released_bet"] = {"arm": _bet_arm(bet, a.bet_strength), "captain_ruling": a.captain_ruling,
                               **bet_record}
    return rec if strong is None else {**rec, "strong_null": strong}


def build_heart(a) -> dict:
    board, dry = _board("heart", a.board)
    genes = board.genes
    if a.trend and (a.entry == "s2" or a.depth or a.herd or a.cfm):
        sys.exit("--trend is refused with --entry s2, --depth, --herd and --cfm (no try read those combinations)")
    copy = a.expression == "copy"
    # The scope bet (HEART_ISO_TRIM; off unless asked for): --entry s2_trim with --expression copy --z-mode raw, the one
    # trimmed copy file. Every other use of --expression copy outside s1 is refused, as before.
    iso_trim = copy and a.entry == "s2_trim" and a.z_mode == "raw"
    if copy and ((a.entry not in ("s1", "core_trim") and not iso_trim) or a.trend or a.depth or a.herd or a.cfm):
        sys.exit("--expression copy is refused with --entry s2, --trend, --depth, --herd and --cfm, and with --entry "
                 "s2_trim unless --z-mode raw is given (it is the copy baseline's draw from --entry s1's pool, the "
                 f"full E10.5, with no step; or {HEART_ISO_TRIM}: the same draw from the half trim, sized in the raw "
                 "mode for H2)")
    # The full trim (HEART_CORE_TRIM; off unless asked for): --entry core_trim builds that one combination only
    core_trim = a.entry == "core_trim"
    if core_trim and (not (copy and a.z_mode == "raw") or a.trend or a.depth or a.herd or a.cfm):
        sys.exit(f"--entry core_trim is refused except as {HEART_CORE_TRIM} (the copy draw from the full trim, sized "
                 "in the raw mode to the core's own forecast), and with --trend, --depth, --herd and --cfm")
    scope_option = HEART_SCOPE_OPTIONS[a.entry] if iso_trim or core_trim else None
    ref = load_stage(DATA / "T2/E10.5.h5ad", genes=genes, drop_duplicates=True)
    start = load_stage(DATA / "T2/E9.5.h5ad", genes=genes, drop_duplicates=True)
    pool, trim = ref, None
    if a.entry == "s2_trim":  # docs/S2_RULES.md (amended 2026-11-02): H2 built as a trim to half E10.5's scope
        core = np.asarray(ref.obs["core"], dtype=bool)
        C = np.asarray(ref.obsm[COORD_KEY], dtype=float)[:, :3]
        keep = t2_build.radial_trim_mask(C, core, 0.5)
        pool = ref[keep].copy()
        trim = {"rule": "surrounding cells removed farthest in plane from the core first, to half the in-plane scope",
                "kept_share": float(keep.mean()), "xy_scope_before": t2_build.xy_scope(C, core),
                "xy_scope_after": t2_build.xy_scope(C[keep], core[keep])}
    elif core_trim:  # the full trim (no ruling): radial_trim_mask at fraction 0
        core = np.asarray(ref.obs["core"], dtype=bool)
        C = np.asarray(ref.obsm[COORD_KEY], dtype=float)[:, :3]
        keep = t2_build.radial_trim_mask(C, core, HEART_SCOPE_OPTIONS["core_trim"]["frac"])
        pool = ref[keep].copy()
        trim = {"rule": HEART_SCOPE_OPTIONS["core_trim"]["rule"], "frac": HEART_SCOPE_OPTIONS["core_trim"]["frac"],
                "kept_share": float(keep.mean()), "xy_scope_before": t2_build.xy_scope(C, core),
                "xy_scope_after": t2_build.xy_scope(C[keep], core[keep])}
    roles = {"ref": ("T2/E10.5.h5ad", 10.5), "start": ("T2/E9.5.h5ad", 9.5)}
    inputs = {"ref": pool, "start": start}
    if a.trend:  # b0_trend (row H11): the earlier stages on the board's panel, loaded as the two above, never trimmed
        for role, (f, day) in HEART_TREND_STAGES.items():
            inputs[role] = load_stage(DATA / f, genes=genes, drop_duplicates=True)
            roles[role] = (f, day)
    view = _view("test:T2h", "T2", "extrap", "heart", board.key, roles, "ref", 12.5)
    count = _t2_count(board, pool.n_obs)
    ctx = harness.BuildContext(task=view, board=harness._candidate_board(board), inputs=inputs, seed=a.seed)
    trend_rec: dict = {}
    if copy:  # the copy baseline: B0's own draw (draw_copy, same context, seed and count) through finish, no step
        stage = harness._copy(ctx)
    elif a.trend:
        stage = harness._b0_t2(ctx, trend=True, trend_record=trend_rec)  # b0_trend, the harness's B0 otherwise
    else:
        stage = harness._b0_t2(ctx, depth=a.depth, herd=a.herd, cfm=a.cfm)  # the harness's B0 with its guard and checks
    hyp = {"s1": "H1", "s2": "H2", "s2_trim": "H2", "core_trim": t2_build.SIZE_CORE_SCOPE}[a.entry]
    z_mode = None if a.z_mode == "recorded" else a.z_mode  # None: the mode the size JSON recorded
    # the size step checks its reference against the full E10.5 cloud it measured; it rescales the (trimmed) cells
    if z_mode is None:
        coords, size = t2_build.size_or_hold(stage.coords, SIZES["heart"], hyp, source_stage="E10.5",
                                             reference_cloud=ref, expect_board="heart")
    else:
        coords, size = t2_build.size_or_hold(stage.coords, SIZES["heart"], hyp, source_stage="E10.5",
                                             reference_cloud=ref, z_mode=z_mode, expect_board="heart")
    _, report = t2_build.build_submission(board, stage.X, coords, genes=list(genes), out=a.out)
    if not report.ok:
        raise SystemExit(f"the built file fails its check: {report.errors}")
    entry = f"copy {a.entry} ({hyp})" if copy else f"B0 {a.entry} ({hyp})"
    return {"entry": entry, "board": board.key, "dry_run": dry, "n_cells": int(stage.n),
            "seed": a.seed, "count_decision": count, "trim": trim, "depth_step": bool(a.depth),
            **({"cfm_step": t2_build.CFM_RULE} if a.cfm else {}),
            **({"trend": _trend_record(trend_rec)} if a.trend else {}),
            **({"expression": "copy"} if copy else {}),
            **({"z_mode": z_mode} if z_mode is not None else {}),
            **({"scope_option": _scope_option_record(a.entry, scope_option)} if scope_option is not None else {}),
            "size": {k: size[k] for k in size if not isinstance(size[k], (dict, list))}}


def _scope_option_record(entry: str, opt: dict) -> dict:
    """What a scope option's build record names (no ruling; no cell number): its scope rule and trim fraction, the
    record field that holds its kept share (a share of the released E10.5, never a count), the size JSON field its size
    comes from, and the halving that follows as a separate step."""
    return {"command": HEART_ISO_TRIM if entry == "s2_trim" else HEART_CORE_TRIM, "ruling": None,
            "scope_rule": opt["rule"], "trim_fraction": opt["frac"], "kept_share_field": "trim.kept_share",
            "size_hypothesis": opt["hypothesis"], "size_source_field": opt["size_source"], "halving": HEART_HALVING}


TREND_LINEAGE_FIELDS = ("lineage", "stages", "days", "own_factor", "weights", "cosine", "norm_ratio",
                        "mean_abs_logf_trend", "mean_abs_logf_b0")


def _trend_record(rec: dict) -> dict:
    """b0_trend's record (vec.harness.b0_trend_factors) for the build record: the rule, the stage roles and days, the
    transition weights, and per lineage the usable stages, cosine and size ratio (L2 norm, trend / B0) against B0's own
    factor and the two mean absolute log factors, plus the medians. Built from named fields only, so no cell number of
    any stage or pool enters it (the harness's per-stage cell numbers and min_cells are left out), on a dry run as on
    the test board."""
    return {"rule": "vec.harness.B0_TREND_RULE", "stages": list(rec["stages"]), "days": list(rec["days"]),
            "weights": [float(w) for w in harness.trend_transition_weights(rec["days"])],
            "budget_days": rec["budget_days"],
            "lineages": [{k: row[k] for k in TREND_LINEAGE_FIELDS} for row in rec["lineages"]],
            "median_cosine": rec["median_cosine"], "median_norm_ratio": rec["median_norm_ratio"]}


def floor_copy_coords(stage_file: str) -> np.ndarray:
    """The floor copy's coordinates as the emulator's bake protocol draws it on a board whose floor is a copy of
    `stage_file` (a T2 path under data/raw): vec.evaluate.subsample of the stage as released (no duplicate drop), at
    SUBSAMPLE_FRAC and seed t2_build.COPY_PARITY_COPY_SEED, rows in released order, float32 as a file holds them."""
    a = read_stage(DATA / stage_file)
    C = np.asarray(a.obsm[COORD_KEY], dtype=np.float32)[:, :3]
    return subsample(Stage(X=np.zeros((len(C), 0), dtype=np.float32), coords=C), SUBSAMPLE_FRAC,
                     t2_build.COPY_PARITY_COPY_SEED).coords


def _copy_rows_refusals(a) -> None:
    """--copy-rows' refusals, before any board or stage is read (off by default: nothing here acts without it)."""
    copy_rows, ok = getattr(a, "copy_rows", False), getattr(a, "copy_rows_test_ok", False)
    sized, sized_ok = getattr(a, "copy_rows_size", False), getattr(a, "copy_rows_size_test_ok", False)
    if sized_ok and not sized:  # row E52's option (off by default)
        raise SystemExit("--copy-rows-size-test-ok is for --copy-rows-size only")
    if sized and not copy_rows:
        raise SystemExit("--copy-rows-size is for --copy-rows only")
    if sized and a.expression != COPY_ROWS_SIZE_EXPRESSION:
        raise SystemExit("--copy-rows-size is for --expression mixture only")
    if ok and not copy_rows:
        raise SystemExit("--copy-rows-test-ok is for --copy-rows only")
    if not copy_rows:
        return
    if a.method != "fgw_shape0" or a.positions != t2_build.POSITIONS[1]:
        raise SystemExit("--copy-rows is for --method fgw_shape0 --positions floor_exact only")
    if (a.row_order != t2_build.ROW_ORDERS[0] or a.shape_interp != SHAPE_INTERP_OFF or a.herd
            or a.draw != t2_build.DRAWS[0] or a.move_coupling != t2_build.MOVE_COUPLINGS[0]
            or a.expression not in COPY_ROWS_EXPRESSIONS):
        raise SystemExit("--copy-rows is refused with --row-order copy_parity, --shape-interp, --herd, --draw "
                         "midpoint, --move-coupling aligned or an --expression other than mixture or series")


def _same_bits(C, ref32: np.ndarray) -> bool:
    """Whether coordinates are ref32 exactly: shape, row order and float32 bits."""
    C32 = np.ascontiguousarray(np.asarray(C)[:, :3], dtype=np.float32)
    R = np.ascontiguousarray(ref32, dtype=np.float32)
    return C32.shape == R.shape and np.array_equal(C32.view(np.uint32), R.view(np.uint32))


def _copy_rows_size_target() -> tuple[float, str, str]:
    """--copy-rows-size: the embryo size JSON's pre-registered full-cloud RMS, read by code with build_embryo's checks
    and message (valid, raw mode), and the file's path and sha256. Refusals show no value."""
    path = SIZES["embryo"]
    raw = Path(path).read_bytes()
    size = json.loads(raw)
    if not size.get("valid") or size["preregistered"]["z_mode"] != "raw":
        raise SystemExit("the embryo size JSON is not a valid raw-mode result")
    target = float(size["preregistered"]["pred_rms_full"])  # exactly as build_embryo reads it for the positions file
    if not np.isfinite(target) or target <= 0:
        raise SystemExit("the embryo size JSON's target is not finite and positive (values not shown)")
    try:
        shown = str(Path(path).resolve().relative_to(REPO_ROOT.resolve()))
    except ValueError:
        shown = str(path)
    return target, shown, hashlib.sha256(raw).hexdigest()


def _copy_rows_rms_ok(S32: np.ndarray, target: float) -> bool:
    """The RMS self-check: the stored float32 coordinates' RMS radius within COPY_ROWS_SIZE_RMS_TOL of the target."""
    return bool(abs(t2_build.rms_radius(np.asarray(S32, dtype=np.float64)) / target - 1.0) <= COPY_ROWS_SIZE_RMS_TOL)


def _copy_rows_parity_ok(S32: np.ndarray, floor32: np.ndarray) -> bool:
    """The prefix parity self-check, count-free: for every count n from COPY_ROWS_PARITY_FROM to the copy's own, the
    SVD parity of the rows the scorer's _match_n takes at metric seed COPY_ROWS_PARITY_SEED from the stored
    coordinates equals the floor copy's (the rows depend on the counts only)."""
    from vec.evaluate import match_rows, svd_parity

    S, F = np.asarray(S32, dtype=np.float32), np.asarray(floor32, dtype=np.float32)
    if S.shape != F.shape:
        return False
    for n in range(COPY_ROWS_PARITY_FROM, len(S) + 1):
        zeros = np.zeros((n, 3))
        if (svd_parity(match_rows(S, zeros, COPY_ROWS_PARITY_SEED)[0])
                != svd_parity(match_rows(F, zeros, COPY_ROWS_PARITY_SEED)[0])):
            return False
    return True


def _build_copy_rows(a, board, dry: bool, genes, files: dict, inputs: dict, count: dict) -> dict:
    """--copy-rows (row E50 as specified; the captain's order of 2026-11-21): the floor copy's rows, order, count and
    float32 coordinates (floor_copy_coords of the earlier bracket as released), the expression of the handed-over
    form's rule (mixture or series). No size rule, centring, jitter or row-order step. The limits are checked here and
    every refusal on this path shows no value; the coordinates are checked bit for bit before writing and on the file
    read back."""
    import anndata as ad

    stage_file = files["before"][0]
    before, after = inputs["before"], inputs["after"]
    released = read_stage(DATA / stage_file, list(genes))
    floor = floor_copy_coords(stage_file)
    Eb, Ea = dense(before[:, genes].X).astype(np.float32), dense(after[:, genes].X).astype(np.float32)
    try:
        dropped = before.uns.get("qc", {}).get("dropped_rows", [])
        res = t2_build.copy_rows_cloud(np.asarray(released.obsm[COORD_KEY])[:, :3], dropped, Eb, Ea,
                                       dense(released.X).astype(np.float32), 0.5, t2_build.COPY_PARITY_COPY_SEED,
                                       a.seed)
    except ValueError as e:
        raise SystemExit(f"--copy-rows refused: {e}") from None
    if not _same_bits(res.coords, floor):
        raise SystemExit("--copy-rows: the coordinates are not the floor copy's, bit for bit (values not shown)")
    sized = getattr(a, "copy_rows_size", False)  # row E52 (off by default): the rescale, then only writing
    coords, written = res.coords, floor
    if sized:
        target, size_file, size_sha = _copy_rows_size_target()
        try:
            coords = t2_build.copy_rows_sized(res.coords, target)
        except ValueError:
            raise SystemExit("--copy-rows-size refused: the copy cannot be rescaled (values not shown)") from None
        written = coords.astype(np.float32)  # what the file stores: round to nearest
        if not _copy_rows_rms_ok(written, target):
            raise SystemExit("--copy-rows-size: the stored coordinates' RMS is not the target's (values not shown)")
        if not _copy_rows_parity_ok(written, floor):
            raise SystemExit("--copy-rows-size: a prefix parity differs from the floor copy's (values not shown)")
    expr_rec = None
    if a.expression == "mixture":
        X = np.empty((len(res.rows), len(genes)), dtype=np.float32)
        X[res.stage == 0], X[res.stage == 1] = Eb[res.rows[res.stage == 0]], Ea[res.rows[res.stage == 1]]
    else:  # series, the handed-over exact-position file's own rule
        X, expr_rec = t2_build.series_rendered_expression(before, after, res, genes, 0.5, [a.seed, 7])
    if len(X) < board.min_cells or (board.max_cells is not None and len(X) > board.max_cells):
        raise SystemExit("the copy-rows count is outside the board's limits (values not shown)")
    try:
        _, report = t2_build.build_submission(board, X, coords, genes=list(genes), out=a.out,
                                              allow_duplicate_coords=True, allow_below_truth=True)
    except ValueError:
        raise SystemExit("build_submission refused the copy-rows file (values not shown)") from None
    if not report.ok:
        raise SystemExit("the built copy-rows file fails its check (values not shown)")
    back = ad.read_h5ad(a.out)
    if not _same_bits(back.obsm[COORD_KEY], written):
        Path(a.out).unlink(missing_ok=True)
        provenance.stamp_path(Path(a.out)).unlink(missing_ok=True)
        if sized:
            raise SystemExit("--copy-rows-size: the written file's coordinates are not the rescaled copy's, bit for "
                             "bit; the file was removed (values not shown)")
        raise SystemExit("--copy-rows: the written file's coordinates are not the floor copy's, bit for bit; the file "
                         "was removed (values not shown)")
    rec = dict(res.record)
    if expr_rec is not None:
        rec.update({f"series_{k}": v for k, v in expr_rec.items() if not isinstance(v, (dict, list))})
    rec["expression"] = a.expression
    out = {"entry": a.method, "board": board.key, "dry_run": dry, "n_cells": int(len(X)), "seed": a.seed,
           "count_decision": count, "target_rms": None,
           "record": {k: v for k, v in rec.items() if not isinstance(v, (dict, list))},
           "copy_rows": True,
           "copy_rows_copy": f"vec.evaluate.subsample({stage_file} as released, {SUBSAMPLE_FRAC}, "
                             f"{t2_build.COPY_PARITY_COPY_SEED})",
           "count_rule": COPY_ROWS_COUNT_RULE, "size_rule": "not applied", "allow_below_truth": True,
           COPY_ROWS_N: int(len(floor))}
    if sized:  # row E52 only: the record says the size step and where its target came from
        out.update({"copy_rows_size": True, "size_rule": COPY_ROWS_SIZE_RULE, "target_rms": target,
                    "size_file": size_file, "size_sha256": size_sha})
    return out


def build_embryo(a) -> dict:
    _copy_rows_refusals(a)
    if a.positions != t2_build.POSITIONS[0] and (a.method != "fgw_shape0" or a.herd
                                                 or a.draw != t2_build.DRAWS[0]):
        raise SystemExit(f"--positions {a.positions} is for --method fgw_shape0, without --herd and --draw midpoint")
    if a.row_order_test_ok and a.row_order != t2_build.ROW_ORDERS[1]:
        raise SystemExit("--row-order-test-ok is for --row-order copy_parity only")
    shape = None if a.shape_interp == SHAPE_INTERP_OFF else a.shape_interp
    if a.shape_interp_test_ok and shape is None:
        raise SystemExit("--shape-interp-test-ok is for --shape-interp with a form only")
    if shape is not None and (a.method != "fgw_shape0" or a.positions != t2_build.POSITIONS[1]):
        raise SystemExit("--shape-interp is for --method fgw_shape0 --positions floor_exact only")
    board, dry = _board("embryo", a.board)
    if (a.row_order != t2_build.ROW_ORDERS[0] and (not dry or str(board.split).startswith("test"))
            and not a.row_order_test_ok):  # before anything loads: nothing about the test count is computed
        raise SystemExit("--row-order copy_parity on a test board needs the captain's yes (--row-order-test-ok); a "
                         "refusal there shows one bit of the test count")
    if (shape is not None and (not dry or str(board.split).startswith("test"))
            and not a.shape_interp_test_ok):  # before anything loads
        raise SystemExit("--shape-interp on a test board needs the captain's yes (--shape-interp-test-ok)")
    if (getattr(a, "copy_rows", False) and (not dry or str(board.split).startswith("test"))
            and not a.copy_rows_test_ok):  # before anything loads: a refusal or a pass there shows one bit
        raise SystemExit("--copy-rows on a test board needs the captain's yes (--copy-rows-test-ok); a refusal or a "
                         "pass there shows one bit about the board's limits")
    if (getattr(a, "copy_rows_size", False) and (not dry or str(board.split).startswith("test"))
            and not a.copy_rows_size_test_ok):  # row E52's option: before anything loads
        raise SystemExit("--copy-rows-size on a test board needs the captain's yes (--copy-rows-size-test-ok)")
    genes = board.genes
    files = {"before": ("T2/E7.5.h5ad", 7.5), "after": ("T2/E8.0.h5ad", 8.0)}
    if a.method == "pseudobulk_shift":  # only it uses the stage before the bracket (E7.25 lacks two panel genes)
        files["earlier"] = ("T2/E7.25.h5ad", 7.25)
    inputs = {r: load_stage(DATA / f, genes=genes, drop_duplicates=True) for r, (f, _) in files.items()}
    count = _t2_count(board, inputs["before"].n_obs)  # provenance; under --copy-rows it does not set the count
    if getattr(a, "copy_rows", False):  # off by default
        return _build_copy_rows(a, board, dry, genes, files, inputs, count)
    view = _view("test:T2e", "T2", "interp", "embryo", board.key, files, "before", 7.75)
    ctx = harness.BuildContext(task=view, board=harness._candidate_board(board), inputs=inputs, seed=a.seed)
    size = json.loads(SIZES["embryo"].read_text())
    if not size.get("valid") or size["preregistered"]["z_mode"] != "raw":
        raise SystemExit("the embryo size JSON is not a valid raw-mode result")
    target = float(size["preregistered"]["pred_rms_full"])  # what apply_size_procedure applies in raw mode
    if a.method == "pseudobulk_shift":
        stage = harness._pseudobulk_shift(ctx)
        coords, rec = t2_build.apply_size_procedure(stage.coords, SIZES["embryo"], expect_board="embryo")
        X = stage.X
    else:
        before, after = inputs["before"], inputs["after"]
        Eb, Ea = dense(before[:, genes].X).astype(np.float32), dense(after[:, genes].X).astype(np.float32)
        res = t2_build.interpolate_clouds(np.asarray(before.obsm[COORD_KEY], float)[:, :3],
                                          np.asarray(after.obsm[COORD_KEY], float)[:, :3], Eb, Ea, 0.5,
                                          harness.count_for(board, before.n_obs), a.seed, target_rms=target,
                                          t_pos=EMBRYO_T_POS[a.method],
                                          with_partners=a.expression in ("blend", "latent"), herd=a.herd,
                                          move_coupling=a.move_coupling, draw=a.draw,
                                          log_depth=(None if a.draw == t2_build.DRAWS[0]
                                                     else (t2_build.log_panel_counts(before, genes),
                                                           t2_build.log_panel_counts(after, genes))),
                                          with_coupling=a.expression == "latent",
                                          **({} if shape is None else {"shape_form": shape}), positions=a.positions)
        expr_rec = None
        if a.expression == "mixture":
            X = np.empty((len(res.rows), len(genes)), dtype=np.float32)
            X[res.stage == 0], X[res.stage == 1] = Eb[res.rows[res.stage == 0]], Ea[res.rows[res.stage == 1]]
        elif a.expression == "depth":  # the harness's fgw_depth (try E7, experiments/harness/TRIES.md)
            X, expr_rec = t2_build.depth_matched_expression(before, after, res, genes, 0.5, [a.seed, 7])
        elif a.expression == "latent":  # the harness's fgw_latent (try E14): torch and scvi-tools load only here
            from vec import t2_latent

            X, expr_rec = t2_latent.latent_midpoint_counts(t2_latent.task_model(before, after, genes), before, after,
                                                           res, genes, 0.5, [a.seed, t2_latent.EXPRESSION_STREAM])
        elif a.expression == "series":  # the harness's series rendering (tries E24a-E24c, E25a-E25b)
            X, expr_rec = t2_build.series_rendered_expression(before, after, res, genes, 0.5, [a.seed, 7])
        elif a.expression == "render_hi":  # the harness's depth top-up to the deeper bracket (tries E19c and E24b)
            X, expr_rec = t2_build.depth_rendered_expression(before, after, res, genes, 0.5, [a.seed, 7], mode="hi")
        else:  # the harness's fgw_blend (try E8)
            X, expr_rec = t2_build.blend_expression(before, after, res, genes, 0.5, [a.seed, 7])
        guard = t2_build.cloud_guard(res.coords, res.pool_coords, seed=a.seed, source_weights=res.pool_weights)
        if not guard.ok:
            raise SystemExit(f"cloud guard flags {guard.flags}")
        coords, rec = res.coords, res.record
        if expr_rec is not None:
            prefix = a.expression if a.expression in ("latent", "series", "render_hi") else "depth"
            rec = {**rec, "expression": a.expression,
                   **{f"{prefix}_{k}": v for k, v in expr_rec.items() if not isinstance(v, (dict, list))
                      and (prefix, k) != ("latent", "retried")}}  # a test build's record holds no cell number
    if a.row_order != t2_build.ROW_ORDERS[0]:  # off by default: the rows re-ordered, nothing else changed
        if not board.truth_cells:
            raise SystemExit("--row-order copy_parity needs the board's truth_cells (read by code)")
        try:
            order, order_rec = t2_build.copy_parity_order(coords, floor_copy_coords(files["before"][0]),
                                                          board.truth_cells, a.seed)
        except ValueError as e:
            raise SystemExit(f"--row-order copy_parity refused: {e}") from None
        X, coords = np.asarray(X)[order], np.asarray(coords)[order]
        rec = {**rec, **order_rec, "row_order_copy": f"vec.evaluate.subsample({files['before'][0]} as released, "
               f"{SUBSAMPLE_FRAC}, {t2_build.COPY_PARITY_COPY_SEED})"}
    _, report = t2_build.build_submission(board, X, coords, genes=list(genes), out=a.out)
    if not report.ok:
        raise SystemExit(f"the built file fails its check: {report.errors}")
    out = {"entry": a.method, "board": board.key, "dry_run": dry, "n_cells": int(len(X)), "seed": a.seed,
           "count_decision": count, "target_rms": target,
           "record": {k: v for k, v in rec.items() if not isinstance(v, (dict, list))}}
    if shape is not None:  # written only when the option is given: flat values, no cell number
        out["shape_interp"] = {"form": shape, "module": "vec.shape_interp", "ruling": None, "row": "E49 (as drafted)",
                               "order": SHAPE_INTERP_ORDER,
                               **{k.removeprefix("shape_interp_"): v for k, v in rec.items()
                                  if k.startswith("shape_interp_") and not isinstance(v, (dict, list))}}
    return out


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="which", required=True)
    for name in ("t1", "heart", "embryo"):
        p = sub.add_parser(name)
        p.add_argument("--out", required=True)
        p.add_argument("--board", help="a validation board for a dry run (default: the test board)")
        p.add_argument("--seed", type=int, default=harness.BUILD_SEED_OFFSET + 1)
        if name == "t1":
            p.add_argument("--cm-share", type=float, help="cardiomyocyte share of heart cells (default: E10.5's own)")
            p.add_argument("--blood-share", type=float, help="blood share of chamber cells (default: E10.5's own)")
            p.add_argument("--herd", action="store_true",
                           help="herd each category's draw towards its pool's mean expression (T1 try T5 of "
                                "experiments/harness/TRIES.md; off: the entry as adopted so far)")
            p.add_argument("--hold-unbracketed", action=argparse.BooleanOptionalAction, default=True,
                           help="R1 (experiments/harness/TRIES.md; the captain's ruling of 2026-11-05): cells without "
                                "a Feng group keep their E10.5 expression (default); --no-hold-unbracketed gives "
                                "them B0's damped E9.5 -> E10.5 momentum, M1 as built before R1")
            p.add_argument("--feng-hygiene", action=argparse.BooleanOptionalAction, default=None,
                           help="R2 (experiments/harness/TRIES.md): also hold R2's gene modules flat in Feng's "
                                "direction (the union N12 wrote, pinned in configs/t1_hygiene.yaml; refused with "
                                "--no-hold-unbracketed or while Blood has a Feng direction). Not given: on with "
                                "--hold-unbracketed (the default), off with --no-hold-unbracketed")
            p.add_argument("--strong-null", action=argparse.BooleanOptionalAction, default=False,
                           help="T12a, T12b (experiments/harness/TRIES.md; the captain's rulings of 2026-11-10 (a) and "
                                "(d)): add the strong-null down-tail prior after M1's steps (off by default; refused "
                                "until configs/t1_strong_null.yaml says approved, and without R1 and R2)")
            p.add_argument("--base", choices=T1_BASES, default=T1_BASES[0],
                           help="the pool the cells are drawn from: pooled = every E10.5 cell (the default); lib3 = "
                                "E10.5's dissected-heart library alone at its own chamber mix (row R6 of "
                                "experiments/harness/TRIES.md; T1's first upload); mid = every library-3 cell and, "
                                "per chamber category, a seeded fraction of the other libraries' cells that puts "
                                "library 3's share half-way to 1 (row R7; T1's second upload, with --cm-share). The "
                                "count, Feng's direction, R1, R2 and g_mom's fit are the same in all three")
            p.add_argument("--lib3-pair", action="store_true",
                           help="with --base lib3 and --cm-share: build the second file of the pair try T32 priced "
                                "(experiments/harness/TRIES.md), library 3's cells with the cardiomyocyte share "
                                "moved. Off by default and no entry: T1's second upload stays on the pooled base "
                                "unless the captain rules")
            p.add_argument("--g-mom", action="store_true",
                           help="add g_mom's within-type factors at test (T1 try T11 of experiments/harness/TRIES.md; "
                                "fitted on E9.5 -> E10.5, applied last, after Feng's direction, on cardiomyocytes and "
                                "endocardium only (R3b), R1's held labels left alone, each cell's expected depth "
                                "kept, no share, cell count or scope changed; refused with --no-hold-unbracketed; "
                                "off: the entry as built, "
                                "and only the captain's ruling turns it on for an upload)")
            p.add_argument("--feng-groups", choices=FENG_GROUPS, default=FENG_GROUPS_DEFAULT,
                           help="the cell groups Feng's direction is formed on: authors = Feng's own labels (the "
                                "default, every file as built); none = one file, the first upload's command "
                                f"({FENG_NONE_COMMAND}) with Feng's direction taken out and nothing put in its place: "
                                "every cell keeps its E10.5 expression before g_mom's edit and no Feng file, R2 list "
                                "or other outside data is read; g_mom's factors edit the cm and endocardium cells "
                                "whose label the readout config's regexes give a group (all but Early VIC, as in "
                                "the entry; the list of groups with a Feng direction is not read). Refused "
                                "without --g-mom and with --base pooled or mid, --cm-share, --blood-share, "
                                "--lib3-pair, --herd, --no-hold-unbracketed, --feng-hygiene or --strong-null. The "
                                "build record then says the direction step was removed. recomputed = the same one "
                                "command on the groups experiments/external/feng_relabel.py recomputes from the "
                                "permitted cells alone, refused unless its committed record routes to recomputed. "
                                "Only the captain's ruling puts either in an upload")
            p.add_argument("--feng-none-ruling", default=None,
                           help="with --feng-groups none only, and only after a permitted-only relabelling whose "
                                "record routes to recomputed: a later docs/PLAN.md decision-log entry of the captain "
                                "(naming N54 and containing \"sets the N54 route aside\"; not the run's own "
                                "ruling; added after the record's commit) that sets the route aside")
            p.add_argument("--g-mom-damping", type=float, default=None,
                           help="candidate A of docs/PLAN.md decision log 2026-11-20 (i) (off by default, no ruling): "
                                f"g_mom's factors finished at 1.0 in place of {harness.B0_DAMPING:g} (practice arm "
                                "gmom_s1, rows T49 and T61). Takes 0.5 (the default, today's file) or 1.0; only with "
                                f"{T1_CANDIDATE_COMMAND} at the default seed on the test board. "
                                f"{G_MOM_DAMPINGS_R22[0]:g} and {G_MOM_DAMPINGS_R22[1]:g} are row {T1_R22_ROW}'s "
                                "second T1 search line (the captain's order of 2026-11-22) and need --captain-ruling "
                                f"{T1_R22_ROW}, with no other change")
            p.add_argument("--feng-fraction-mult", type=float, default=None,
                           help="candidate B of docs/PLAN.md decision log 2026-11-20 (i) (off by default, no ruling): "
                                "the outside atlas's direction at this multiple of the calendar fraction (2/5). Takes "
                                "1.0 (the default, today's file) or 1.5; only with the command "
                                f"{T1_CANDIDATE_COMMAND} at the default seed on the test board. "
                                f"{FENG_FRACTION_MULT_R20:g} is row {T1_R20_ROW}'s file (the captain's order of "
                                f"2026-11-21) and needs --captain-ruling {T1_R20_ROW}, with no other change; "
                                f"{FENG_FRACTION_MULTS_R22[0]:g} and {FENG_FRACTION_MULTS_R22[1]:g} are row "
                                f"{T1_R22_ROW}'s first search line (the captain's order of 2026-11-22) and need "
                                f"--captain-ruling {T1_R22_ROW}, with no other change")
            p.add_argument("--released-bet", choices=tuple(RELEASED_BETS), default=None,
                           help="try T61 of experiments/harness/TRIES.md (off by default): a variant of T1's last file "
                                "(--feng-groups none only) that changes which cells are drawn: share_trend (the "
                                "chamber mix moved along the E9.5 -> E10.5-portions trend) or advance (cells drawn "
                                "inside each category towards the released momentum axis). Composition, a lever "
                                "closed by the post-exposure block: refused unless --captain-ruling names the "
                                "captain's ruling row for the arm")
            p.add_argument("--bet-strength", type=float, default=None,
                           help="with --released-bet: lambda (share_trend) or tau (advance), one of row T61's values")
            p.add_argument("--captain-ruling", default=None,
                           help="with --released-bet: the ID of the row of experiments/harness/TRIES.md that is the "
                                "captain's ruling for the arm; with --feng-fraction-mult "
                                f"{FENG_FRACTION_MULT_R20:g}: {T1_R20_ROW}, the captain's order for that file; with "
                                f"{FENG_FRACTION_MULTS_R22[0]:g} or {FENG_FRACTION_MULTS_R22[1]:g}: {T1_R22_ROW}, the "
                                "captain's order for the first search line; with --g-mom-damping "
                                f"{G_MOM_DAMPINGS_R22[0]:g} or {G_MOM_DAMPINGS_R22[1]:g}: {T1_R22_ROW}, the captain's "
                                "order for the second T1 search line")
        if name == "heart":
            p.add_argument("--entry", choices=("s1", "s2", "s2_trim", "core_trim"), required=True,
                           help="s1 = H1; s2 = H2 as first built (size only); s2_trim = H2 as a trim to half the "
                                "scope, the s2 of docs/S2_RULES.md; core_trim = the full trim (no ruling; only as "
                                f"{HEART_CORE_TRIM}, sized to the raw mode's own core forecast)")
            p.add_argument("--herd", action="store_true",
                           help="herd B0's draw towards the pool's mean expression (try H4; off: the entry so far)")
            p.add_argument("--depth", action="store_true",
                           help="thin B0's cells to the extrapolated E12.5 depth (heart try H3 of "
                                "experiments/harness/TRIES.md; off: the entry as adopted so far)")
            p.add_argument("--cfm", action="store_true",
                           help="filter each lineage's momentum by its own within-stage covariance (heart try H6 of "
                                "experiments/harness/TRIES.md; off: the entry as adopted so far)")
            p.add_argument("--trend", action="store_true",
                           help="b0_trend: each lineage's momentum the least-squares trend over E8.25_late, E8.5, "
                                "E8.75, E9.5 and E10.5, over B0's time budget (the captain's ruling of 2026-11-13, "
                                "docs/PLAN.md decision log (w); ledger row H11 of experiments/harness/TRIES.md); "
                                "with --entry s1 or s2_trim, refused with --entry s2, --depth, --herd and --cfm")
            p.add_argument("--expression", choices=("b0", "copy"), default="b0",
                           help="b0 = B0's momentum step (the default, the entry as built); copy = the copy baseline's "
                                "cells and expression, B0's own draw with no step (the captain's direction for the "
                                "heart's second upload, docs/PLAN.md decision log 2026-11-18 (k)); --entry s1, or "
                                "--entry s2_trim with --z-mode raw (the scope bet, the same draw from the half trim "
                                "sized for H2; no ruling, scope is a closed lever, only the captain's ruling puts it "
                                "in an upload), or --entry core_trim with --z-mode raw (the full trim, likewise no "
                                "ruling); refused with --entry s2, --trend, --depth, --herd and --cfm")
            p.add_argument("--z-mode", choices=("recorded", "raw"), default="recorded",
                           help="the size step's z mode: recorded = the mode the size JSON pre-registered (the "
                                "default); raw = the isotropic mode the JSON also logs, one factor on all three axes "
                                "to that raw mode's own target RMS for the hypothesis (z_modes.raw, not the recorded "
                                "mode's pre-registered target); refused when the size JSON is not valid (its "
                                "fallback is already isotropic: build with recorded then)")
        if name == "embryo":
            p.add_argument("--method", choices=("pseudobulk_shift", *EMBRYO_T_POS), required=True,
                           help="fgw = s1; fgw_shape0 = the s2 of docs/S2_RULES.md (E7.5's geometry)")
            p.add_argument("--herd", action="store_true",
                           help="herd each bracket's draw towards its mean expression (try E9; off: the entry so far)")
            p.add_argument("--expression", choices=("mixture", "depth", "blend", "latent", "series", "render_hi"),
                           default="mixture",
                           help="fgw methods: the cells as released (the entry), or tries E7, E8, E14, E24 and E19c "
                                "of experiments/harness/TRIES.md: thinned to the predicted depth; blends; a count "
                                "model's latent interpolation (needs the gpu dependency group); the shallower "
                                "bracket's cells in the deeper bracket's run profile; every cell topped up to the "
                                "deeper bracket's depth. Neither of the last two is an entry until the captain rules")
            p.add_argument("--move-coupling", choices=t2_build.MOVE_COUPLINGS, default=t2_build.MOVE_COUPLINGS[0],
                           help="fgw methods: which coupling sets the moves: fgw (the entry) or aligned (try E10 of "
                                "experiments/harness/TRIES.md: an exact transport in the rotated common frame)")
            p.add_argument("--draw", choices=t2_build.DRAWS, default=t2_build.DRAWS[0],
                           help="fgw methods: which real cells of each bracket are drawn: uniform (the entry) or "
                                "midpoint (try E13 of experiments/harness/TRIES.md: tilted towards the local "
                                "displacement midpoint, every cell unedited at its moved position)")
            p.add_argument("--positions", choices=t2_build.POSITIONS, default=t2_build.POSITIONS[0],
                           help="--method fgw_shape0 only: where the cells sit: moved (the default, the entry) or "
                                "floor_exact (tries E47a and E47b of experiments/harness/TRIES.md: every cell on a "
                                "real E7.5 position, the E8.0 cells on unused ones by an exact expression assignment; "
                                "with --expression series the harness's fgw_shape0x_series). No entry until the "
                                "captain rules; refused with --herd and --draw midpoint")
            p.add_argument("--row-order", choices=t2_build.ROW_ORDERS, default=t2_build.ROW_ORDERS[0],
                           help="as_built (the default, every file as built) or copy_parity: the rows re-ordered so "
                                "the file's pre-fix SVD parity on the scorer's matched rows is the floor copy's "
                                "(vec.t2_build.copy_parity_order; branch campaign/embryo-occupancy). Same cells, "
                                "expression and positions. No try and no ruling: only the captain's ruling puts it "
                                "in an upload. It assumes the server scores the prediction's rows as the file holds "
                                "them (file order, full count, no subsample or sort); score_h5ad's _prediction does "
                                "this in the pinned copy, but the live server's handling is not shown")
            p.add_argument("--shape-interp", choices=(SHAPE_INTERP_OFF, *shape_interp.FORMS),
                           default=SHAPE_INTERP_OFF,
                           help="--method fgw_shape0 --positions floor_exact only: none (the default, every file as "
                                "built) or one shape form of vec.shape_interp between E7.5 and E8.0 at t = 1/2 on the "
                                "exact-position cloud, before the size rule (row E49 as drafted). No try read and no "
                                "ruling: only the captain's ruling puts it in an upload")
            p.add_argument("--shape-interp-test-ok", action="store_true",
                           help="--shape-interp on the test board: the captain's yes (off by default; without it the "
                                "option refuses on a test board before anything loads)")
            p.add_argument("--copy-rows", action="store_true",
                           help="--method fgw_shape0 --positions floor_exact only, with --expression mixture or "
                                "series (off by default; the captain's order of 2026-11-21, docs/PLAN.md decision log "
                                "2026-11-21 (h); row E50 as specified): the file takes the floor copy's own rows, "
                                "order, count and float32 coordinates (vec.evaluate.subsample of E7.5 as released, "
                                "0.1, seed 0), half its slots keep their own cell's expression and the rest take E8.0 "
                                "cells by an exact assignment (vec.t2_build.copy_rows_cloud). It departs from the "
                                "count rule and the size rule. Refused with --row-order copy_parity, --shape-interp, "
                                "--herd, --draw midpoint and --move-coupling aligned")
            p.add_argument("--copy-rows-test-ok", action="store_true",
                           help="--copy-rows on the test board: the captain's yes (off by default; without it the "
                                "option refuses on a test board before anything loads, because a refusal or a pass "
                                "there shows one bit about the board's limits)")
            p.add_argument("--copy-rows-size", action="store_true",
                           help="with --copy-rows and --expression mixture only (off by default; the captain's order "
                                "of 2026-11-21, docs/PLAN.md decision log 2026-11-21 (v); row E52 as specified): the "
                                "copy-rows file with its coordinates scaled about the copy's own centroid to the "
                                "embryo size JSON's pre-registered full-cloud RMS (raw mode), the last step before "
                                "writing (vec.t2_build.copy_rows_sized)")
            p.add_argument("--copy-rows-size-test-ok", action="store_true",
                           help="--copy-rows-size on the test board: the captain's yes (off by default; a test build "
                                "needs it and --copy-rows-test-ok, before anything loads)")
            p.add_argument("--row-order-test-ok", action="store_true",
                           help="--row-order copy_parity on the test board: the captain's yes (off by default; "
                                "without it the option refuses on a test board before anything loads, because a "
                                "refusal or a success there shows one bit of the test count)")
    return ap


def main(argv=None) -> None:
    a = build_parser().parse_args(argv)
    _check_seed(a.seed)
    if a.board not in (None, TEST_BOARDS[a.which]):  # a dry run
        if a.board not in DRY_BOARDS[a.which]:
            sys.exit(f"--board {a.board}: a {a.which} dry run takes only {sorted(DRY_BOARDS[a.which])}")
        if "dry" not in Path(a.out).name.lower():
            sys.exit(f"{a.out}: a dry run's file name must contain 'dry' (the file must never be uploaded)")
    elif "dry" in Path(a.out).name.lower():
        sys.exit(f"{a.out}: a test build's file name must not contain 'dry'")
    if getattr(a, "strong_null", False) and Path(a.out).exists():
        sys.exit(f"{a.out} exists: a build with --strong-null never overwrites a file (T12's row)")
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    rec = {"t1": build_t1, "heart": build_heart, "embryo": build_embryo}[a.which](a)
    rec = _hide_counts(rec, load_boards()[rec["board"]])
    Path(a.out).with_suffix(".build.json").write_text(json.dumps(rec, indent=1, default=float))
    print(json.dumps({k: v for k, v in rec.items() if k not in ("shares", "external", "realised", "strong_null")},
                     default=float))


if __name__ == "__main__":
    main()
