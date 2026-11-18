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
          counts. --g-mom-fit lib3 (with --g-mom; default whole; no entry, only a captain's ruling would use it): the
          same fit on E10.5's dissected-heart library alone
  heart   B0 (copy of E10.5 plus damped E9.5 -> E10.5 momentum in count space), sized by the pre-registered heart
          procedure: --entry s1 = scope hypothesis H1; s2_trim = the s2 of docs/S2_RULES.md (H2 as a trim to half
          E10.5's scope); s2 = H2 as first built (size only, superseded) (approved rule 3, `size_or_hold`). --trend
          (both heart files since the captain's ruling of 2026-11-13, docs/PLAN.md decision log (w); row H11 of
          experiments/harness/TRIES.md): b0_trend, each lineage's momentum the least-squares trend over E8.25_late to
          E10.5 in place of the E9.5 -> E10.5 step alone
  embryo  --method pseudobulk_shift | fgw | fgw_shape0 | fgw_shape1: E7.5 cells moved by the E7.25 -> E7.5 change, or
          the FGW interpolation between E7.5 and E8.0 at t = 1/2 (fgw_shape0, the s2 of docs/S2_RULES.md, keeps the
          mixture in E7.5's geometry); sized by the pre-registered embryo procedure (size_embryo_E7.75.json)

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
import json
import sys
from pathlib import Path

import numpy as np

from vec import harness, provenance, t1_scope, t2_build
from vec.boards import REPO_ROOT, load_boards
from vec.io import COORD_KEY, dense
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
    if not board.counts_hidden:
        return rec
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
    if isinstance(out.get("g_mom"), dict) and "edited_share" in out["g_mom"]:  # the drawn cells' share (--g-mom-fit):
        # its denominator follows from the count, as the realised composition's does; the pool's share stays
        out["g_mom"] = {k: v for k, v in out["g_mom"].items() if k != "edited_share"}
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


def _feng_group_names() -> list[str]:
    """The names of Feng's direction groups from Feng's statistics alone, before any T1 file is read. The groups
    feng_dirs.feng_directions returns do not depend on the gene list it is given, so it is called on one placeholder
    name (an empty list fails in dirlib.sym_index, whose empty index is a float array) and only the names are kept."""
    return sorted(feng_dirs.feng_directions([""])[0])


def _t1_union(hygiene_on: bool) -> dict | None:
    """R2's union for --feng-hygiene, after R2's two refusals that need no T1 file, both made before anything else is
    read (R2: no E10.5 file is read before the build): hygiene.load_union refuses while configs/t1_hygiene.yaml pins no
    union_sha256 (until the build commit after N12) or pins another file's; then the Blood guard on the names of Feng's
    direction groups. Off: None, and nothing is read."""
    if not hygiene_on:
        return None
    try:
        union = hygiene.load_union()
    except hygiene.Refused as e:
        sys.exit(f"--feng-hygiene is refused: {e}")
    _no_blood(_feng_group_names())
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


def _t1_library3(ref, what: str, need_others: bool) -> np.ndarray:
    """Which reference cells are E10.5's dissected-heart library (names ending in _T1_DISSECTED_LIBRARY), the selection
    of the lib3 base. Every name must end in _<library number> and library 3 must be among them; need_others (the
    bases): library 3 must be one of at least two. Otherwise the build stops."""
    names = [str(x) for x in ref.obs_names]
    library = np.array([x.rsplit("_", 1)[-1] if "_" in x else "" for x in names])
    kinds = set(library.tolist())
    if not all(k.isdigit() for k in kinds) or T1_DISSECTED_LIBRARY not in kinds or (need_others and len(kinds) < 2):
        sys.exit(f"{what}: the reference's cell names do not separate library 3 from the others (every name must end "
                 "in _<library number>, and library 3 must be " + ("one of at least two)" if need_others
                                                                   else "among them)"))
    return library == T1_DISSECTED_LIBRARY


def _edited_shares(cats: np.ndarray, n_edited: dict[str, int], weights: dict[str, float] | None = None) -> dict:
    """Per category of G_MOM_EDIT_CATEGORIES and over both ("all"): the share of the cells in that category that g_mom's
    hook edits (None where the category holds no cell). Shares only, no cell number. Without `weights` (the drawn
    cells), "all" is the share over both categories' cells as they stand. With `weights` (the build's target shares;
    the pool), "all" is sum_c w_c x share_c / sum_c w_c over the categories with a share and w_c > 0: the expected share
    of the drawn file's R3b cells that the hook edits, if each category's cells are drawn uniformly from the pool and
    no stand-ins are used; the pool's own count weights would not be that expectation once the build's shares differ
    from the pool's mix (--cm-share)."""
    out = {}
    for c in G_MOM_EDIT_CATEGORIES:
        total = int((cats == c).sum())
        out[c] = n_edited.get(c, 0) / total if total else None
    if weights is None:
        total = int(np.isin(cats, list(G_MOM_EDIT_CATEGORIES)).sum())
        out["all"] = sum(n_edited.get(c, 0) for c in G_MOM_EDIT_CATEGORIES) / total if total else None
    else:
        w = {c: float(weights.get(c, 0.0)) for c in G_MOM_EDIT_CATEGORIES if out[c] is not None}
        w_sum = sum(v for v in w.values() if v > 0)
        out["all"] = sum(v * out[c] for c, v in w.items() if v > 0) / w_sum if w_sum > 0 else None
    return {"all": out.pop("all"), **out}


G_MOM_FITS = ("whole", "lib3")  # --g-mom-fit: the reference g_mom's factors are fitted on
_G_MOM_DOSE_READING = (
    "Of the lib3 files only one on the lib3 base keeps the step's full dose, as the whole fit gives it. On a base "
    "other than lib3 the drawn cells whose source cell sits in no fitted metacell keep their counts, so the step's "
    "dose on the file falls to edited_share_pool (with lib3: on the mid base at most about (1 + s) / 2 of a category's"
    " cells, on the pooled base at most s, s being library 3's share of that category in the file, R1's held labels "
    "apart): such a file changes both what the factors carry and how many cells get any factor, and is not the same "
    "step without the collection part; a comparison of files names the base of each. ")
_G_MOM_NO_EVIDENCE = (
    "No fold, backtest or listing reads this option: it is a rehearsable code path with no evidence about the score, "
    "whose change it makes is unknown in sign; every figure this record gives (the option, the metacell counts, the "
    "shares) is a property of the construction, not a test")
G_MOM_LIB3_READING = (
    "the lib3 fit differs from the whole fit in its reference's size as well as its library make-up (fewer metacells, "
    "its own PCA, HVGs, OT marginals and per-category mean). A size-matched control (the same fit on random subsets of"
    " E10.5 with library 3's cell number per category, several draws with their spread, in a target-free listing; not "
    "built here) would remove size only: what it leaves is the part that goes with library make-up (states the "
    "libraries hold in different amounts, and any offset), net of size, which is not the same as the collection part. "
    "A library-share listing on this fit is zero by construction and is not a test that the collection part is gone. A"
    " uniform offset common to library 3 is about 1 to 3% of the factors' squared size (cardiomyocytes 1.1 to 2.8%, "
    "endocardium 1.2 to 1.5%; docs/PLAN.md decision log 2026-11-18 (i), on main at 3182da3). The rest of the "
    "library-linked part is a difference between states the libraries hold in different amounts; released stages "
    "cannot split it into development, identity or collection. The lib3 fit's leading axis has |cos| 0.81 "
    "(cardiomyocytes) and 0.75 (endocardium) with the fit as built's library slope (a libraries-1-and-2 fit gives 0.81"
    " and 0.76, the fit as built's own leading axis 0.98 and 0.81; no null for that cosine was computed, so it is an "
    "axis agreement only; experiments/t1_m1/results/gmom_test_fit/single_kind.txt, on main at 55d862c). "
    + _G_MOM_DOSE_READING + _G_MOM_NO_EVIDENCE)
EDITED_SHARE_POOL_ASSUMES = (
    "each category's drawn cells are drawn uniformly from the pool's cells of that category, and no stand-in is drawn "
    "(a stand-in's category differs from its source cell's, so it gets no factor yet counts in the pool share)")


def _t1_g_mom(ref, categories: list[str], seed: int, covered_labels: set[str], fit: str = G_MOM_FITS[0],
              pool=None, shares: dict[str, float] | None = None):
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
    of Feng-covered labels cannot disagree silently.

    fit (--g-mom-fit): the option the caller chose; `ref` is the reference the fit runs on (the caller passes the whole
    E10.5 object for "whole" and its library-3 rows for "lib3"), with the same start, function and `categories`. With
    "whole" the record and the hook are exactly as before the option (no field added), so every output of the default
    build, its .build.json included, is as on main. With "lib3" the record adds the option, the number of metacells
    per fitted category, a fixed reading (G_MOM_LIB3_READING: what the option's figures can and cannot show), and the
    share of the cells in R3b's edit categories that the hook edits: "edited_share_pool" over `pool` (the build's base,
    when given: its cells whose name sits in a fitted metacell and whose label R1 does not hold; its "all" weights the
    categories by `shares`, the build's target shares, see `_edited_shares`), and "edited_share" over the drawn cells,
    which `_hide_counts` strips on a test board (its denominator follows from the board's count). A drawn cell whose
    source cell is outside the fit keeps its counts and is not edited. The fit's other inputs follow from the smaller
    reference (fewer metacells, its own PCA, OT and per-category mean), so the lib3 fit differs from the whole one in
    more than library make-up. "edited_share_pool" is the expected share only under EDITED_SHARE_POOL_ASSUMES
    (uniform draws within a category, no stand-in), which the record states beside it ("edited_share_pool_assumes").
    No size control is built here: a lib3-against-size-matched comparison leaves what goes with library make-up
    (states the libraries hold in different amounts, and any offset), net of size, which is not the same as the
    collection part, and needs several draws read in a target-free listing (G_MOM_LIB3_READING)."""
    from vec import t1_state_transport as ts

    start = load_stage(DATA / "T1/E9.5_RNA.h5ad", lineage=False)
    factors, _ = harness.g_mom_pair_factors(ref, "E10.5_RNA", start, "E9.5_RNA", categories=categories)
    missing = sorted(set(G_MOM_EDIT_CATEGORIES) - set(factors.categories))
    if missing:
        raise ValueError(f"--g-mom has no factors for a category R3b edits: {missing}")
    held_flat = t1_scope.held_flat_genes()
    whole = fit == G_MOM_FITS[0]
    fit_record = {} if whole else {
        "fit": fit,
        "fit_reference": (f"E10.5's library-{T1_DISSECTED_LIBRARY} cells only (names ending in "
                          f"_{T1_DISSECTED_LIBRARY}, the lib3 base's selection); the fit's other inputs follow from "
                          "the smaller reference (fewer metacells, its own PCA, OT and per-category mean)"),
        "metacells": {c: int(factors.categories[c].n_metacells) for c in sorted(factors.categories)},
        "reading": G_MOM_LIB3_READING}
    if pool is not None and not whole:
        b = harness._shallow(pool)
        t1_scope.classify(b, source="E10.5_RNA", spec=t1_scope.load_spec())
        p_cats = np.asarray(b.obs[t1_scope.CATEGORY_KEY]).astype(str)
        p_labels = np.asarray(b.obs[t1_scope.LABEL_KEY]).astype(str)
        p_held = np.isin(p_labels, sorted(ts.r1_unbracketed_labels(np.unique(p_labels))))
        p_names = np.asarray(b.obs_names).astype(str)
        fitted = {c: (p_cats == c) & ~p_held & np.isin(p_names, list(factors.categories[c].index))
                  for c in G_MOM_EDIT_CATEGORIES}
        fit_record["edited_share_pool"] = _edited_shares(p_cats, {c: int(m.sum()) for c, m in fitted.items()},
                                                         weights=shares or {})
        fit_record["edited_share_pool_assumes"] = EDITED_SHARE_POOL_ASSUMES
    record = {"rule": "try T11's within-type factors at test: OT within each chamber category from E9.5 to E10.5, "
                      "metacells of the E10.5 cells, the per-category mean removed, damped by B0_DAMPING, zeroed and "
                      "clipped; applied last (after Feng's direction); R1's held labels unedited; keep_depth; R3b: "
                      "only the categories in edit_categories are edited; the factors are not scaled to the two-day "
                      "step",
              "edit_categories": list(G_MOM_EDIT_CATEGORIES),
              "categories_with_factors": sorted(factors.categories),
              "categories_without_factors": sorted(factors.skipped), "stream": [int(seed), 5],
              "held_flat": "the candidate's pre-registered list (R2's union is held in Feng's direction only)",
              **fit_record}

    def hook(pop):
        labels = pop.obs[t1_scope.LABEL_KEY].astype(str).to_numpy()
        cats = pop.obs[t1_scope.CATEGORY_KEY].astype(str).to_numpy()
        no_edit = (np.isin(labels, sorted(ts.r1_unbracketed_labels(np.unique(labels))))
                   | ~np.isin(cats, G_MOM_EDIT_CATEGORIES))
        edited = set(labels[~no_edit & np.isin(cats, list(factors.categories))])
        if not edited <= set(covered_labels):
            extra = sorted(edited - set(covered_labels))
            raise ValueError(f"--g-mom would edit labels Feng's direction does not cover: {extra}")
        moved, rec = ts.apply_state_factors(pop, factors, held_flat, [seed, 5], rescale=False, no_edit=no_edit)
        record["median_depth_scale"] = rec["median_depth_scale"]  # no cell number
        if not whole:  # the default's record stays as before the option
            record["edited_share"] = _edited_shares(cats, {c: v["n_edited"] for c, v in rec["by_category"].items()})
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


def build_t1(a) -> dict:
    if a.g_mom_fit != G_MOM_FITS[0] and not a.g_mom:
        sys.exit(f"--g-mom-fit {a.g_mom_fit} is refused without --g-mom: it chooses the reference g_mom's factors are "
                 "fitted on")
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
    hygiene_on = _feng_hygiene_on(a.feng_hygiene, a.hold_unbracketed)  # R2's refusals come before anything is read
    prior = _t1_down_prior(a.strong_null, a.hold_unbracketed, hygiene_on)  # T12's, likewise
    union = _t1_union(hygiene_on)
    board, dry = _board("t1", a.board)
    ref = load_stage(DATA / "T1/E10.5_RNA.h5ad", lineage=False)
    momentum = _t1_momentum(a.hold_unbracketed)
    shares = harness.t1_chamber_shares(harness._shallow(ref), "E10.5_RNA")
    fit_categories = [c for c, s in shares.items() if s > 0]  # g_mom's fit, as the candidate's: E10.5's own mix
    pool = ref
    mid_record = None
    if a.base in ("lib3", "mid"):
        in3 = _t1_library3(ref, f"--base {a.base}", need_others=True)
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
    directions, ext = feng_dirs.feng_directions(list(ref.var_names))
    held, feng_hygiene = _t1_held(union, directions, ref.var_names)
    groups = feng_dirs.label_groups(ref.obs["celltype"].astype(str).unique())
    fraction = feng_dirs.bracket_fraction(10.5, 12.5)
    n = harness.count_for(board, ref.n_obs)
    extra = {} if prior is None else {"down_prior": prior}
    within_type, g_mom_record = ({}, None)
    if a.g_mom:
        covered = {lab for lab, group in groups.items() if group in directions}  # the build's Feng-covered labels
        fit_ref = ref  # --g-mom-fit whole: the fit's reference is the whole E10.5 object, as before the option
        if a.g_mom_fit == "lib3":  # the same fit on E10.5's library-3 cells alone (the lib3 base's selection)
            fit_ref = ref[np.flatnonzero(_t1_library3(ref, "--g-mom-fit lib3", need_others=False))].copy()
        hook, g_mom_record = _t1_g_mom(fit_ref, fit_categories, a.seed, covered, fit=a.g_mom_fit, pool=pool,
                                       shares=shares)
        within_type = {"post_cells": hook}
    sub, rep = t1_scope.build_submission(pool, board, shares, n, a.seed, held_flat=held,
                                         ref_source="E10.5_RNA", directions=directions, direction_fraction=fraction,
                                         direction_groups=groups, herd=a.herd, **momentum, **extra, **within_type)
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
    if mid_record is not None:  # likewise: shares only, no cell number
        rec["mid_base"] = {"stream": [int(a.seed), T1_MID_STREAM], "categories": mid_record}
    if a.g_mom:
        rec["g_mom"] = g_mom_record
    return rec if strong is None else {**rec, "strong_null": strong}


def build_heart(a) -> dict:
    board, dry = _board("heart", a.board)
    genes = board.genes
    if a.trend and (a.entry == "s2" or a.depth or a.herd or a.cfm):
        sys.exit("--trend is refused with --entry s2, --depth, --herd and --cfm (no try read those combinations)")
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
    if a.trend:
        stage = harness._b0_t2(ctx, trend=True, trend_record=trend_rec)  # b0_trend, the harness's B0 otherwise
    else:
        stage = harness._b0_t2(ctx, depth=a.depth, herd=a.herd, cfm=a.cfm)  # the harness's B0 with its guard and checks
    hyp = {"s1": "H1", "s2": "H2", "s2_trim": "H2"}[a.entry]
    # the size step checks its reference against the full E10.5 cloud it measured; it rescales the (trimmed) cells
    coords, size = t2_build.size_or_hold(stage.coords, SIZES["heart"], hyp, source_stage="E10.5", reference_cloud=ref,
                                         expect_board="heart")
    _, report = t2_build.build_submission(board, stage.X, coords, genes=list(genes), out=a.out)
    if not report.ok:
        raise SystemExit(f"the built file fails its check: {report.errors}")
    return {"entry": f"B0 {a.entry} ({hyp})", "board": board.key, "dry_run": dry, "n_cells": int(stage.n),
            "seed": a.seed, "count_decision": count, "trim": trim, "depth_step": bool(a.depth),
            **({"cfm_step": t2_build.CFM_RULE} if a.cfm else {}),
            **({"trend": _trend_record(trend_rec)} if a.trend else {}),
            "size": {k: size[k] for k in size if not isinstance(size[k], (dict, list))}}


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


def build_embryo(a) -> dict:
    board, dry = _board("embryo", a.board)
    genes = board.genes
    files = {"before": ("T2/E7.5.h5ad", 7.5), "after": ("T2/E8.0.h5ad", 8.0)}
    if a.method == "pseudobulk_shift":  # only it uses the stage before the bracket (E7.25 lacks two panel genes)
        files["earlier"] = ("T2/E7.25.h5ad", 7.25)
    inputs = {r: load_stage(DATA / f, genes=genes, drop_duplicates=True) for r, (f, _) in files.items()}
    count = _t2_count(board, inputs["before"].n_obs)
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
                                          with_coupling=a.expression == "latent")
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
    _, report = t2_build.build_submission(board, X, coords, genes=list(genes), out=a.out)
    if not report.ok:
        raise SystemExit(f"the built file fails its check: {report.errors}")
    return {"entry": a.method, "board": board.key, "dry_run": dry, "n_cells": int(len(X)), "seed": a.seed,
            "count_decision": count, "target_rms": target,
            "record": {k: v for k, v in rec.items() if not isinstance(v, (dict, list))}}


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
            p.add_argument("--g-mom-fit", choices=G_MOM_FITS, default=G_MOM_FITS[0],
                           help="with --g-mom: the reference g_mom's factors are fitted on (start E9.5, the same "
                                "categories): whole = every E10.5 cell (the default, the entry as built); lib3 = "
                                "E10.5's dissected-heart library alone (names ending in _3, the lib3 base's "
                                "selection), so no contrast between metacells of different library make-up enters "
                                "the factors by construction. A uniform offset common to library 3 is about 1 to 3%% "
                                "of the factors' squared size (cardiomyocytes 1.1 to 2.8%%, endocardium 1.2 to 1.5%%; "
                                "decision log 2026-11-18 (i), on main at 3182da3); the rest of the library-linked "
                                "part is a difference between states the libraries hold in different amounts, which "
                                "released stages cannot split into development, identity or collection. The lib3 "
                                "fit's leading axis has |cos| 0.81 (cardiomyocytes) and 0.75 (endocardium) with the "
                                "fit as built's library slope (a libraries-1-and-2 fit gives 0.81 and 0.76; no null "
                                "for that cosine was computed, so it is an axis agreement only; single_kind.txt, on "
                                "main at 55d862c). The fit's other inputs follow from the smaller "
                                "reference (fewer metacells, its own PCA, OT and per-category mean), so a lib3 file "
                                "differs from the whole one in more than the library part; a size-matched control "
                                "(not built here; several draws, in a target-free listing) would remove size only, "
                                "and what it leaves is what goes with library make-up (states the libraries hold in "
                                "different amounts, and any offset), net of size, which is not the same as the "
                                "collection part. "
                                "Of the "
                                "lib3 files only one on the lib3 base keeps the step's full dose, as the whole fit "
                                "gives it. With lib3 and a base other "
                                "than lib3 (pooled, mid), the drawn "
                                "cells whose source cell is of another library sit in no fitted metacell and keep "
                                "their counts, a cut in the step's dose (on mid at most about (1 + s) / 2 of a "
                                "category's cells edited, on pooled at most s, s its library-3 share in the file, "
                                "R1's held labels apart): the record's "
                                "g_mom field gives the share of R3b's cells the step edits (edited_share_pool, its "
                                "'all' weighted by the build's shares, under the assumption it states; edited_share "
                                "on a dry run). No fold, backtest or listing reads the option: it is a rehearsable "
                                "code path with no evidence about the score, its "
                                "effect on the score is unknown in sign, and a library-share listing on a lib3 fit is "
                                "zero by construction and not a test. With whole the record is as before the option. "
                                "Refused without --g-mom; no entry uses it, and only the captain's ruling puts it in "
                                "an upload")
        if name == "heart":
            p.add_argument("--entry", choices=("s1", "s2", "s2_trim"), required=True,
                           help="s1 = H1; s2 = H2 as first built (size only); s2_trim = H2 as a trim to half the "
                                "scope, the s2 of docs/S2_RULES.md")
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
