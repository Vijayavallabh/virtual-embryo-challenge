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
          counts
  heart   B0 (copy of E10.5 plus damped E9.5 -> E10.5 momentum in count space), sized by the pre-registered heart
          procedure: --entry s1 = scope hypothesis H1; s2_trim = the s2 of docs/S2_RULES.md (H2 as a trim to half
          E10.5's scope); s2 = H2 as first built (size only, superseded) (approved rule 3, `size_or_hold`)
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
  uv run python experiments/builds/build_test.py heart --entry s1 --out runs/builds/T2heart_s1.h5ad
  uv run python experiments/builds/build_test.py embryo --method fgw --out runs/builds/T2embryo_s1.h5ad
  uv run python experiments/builds/build_test.py t1 --cm-share 0.60 --out runs/builds/T1_s2.h5ad
  uv run python experiments/builds/build_test.py heart --entry s2_trim --out runs/builds/T2heart_s2.h5ad
  uv run python experiments/builds/build_test.py embryo --method fgw_shape0 --out runs/builds/T2embryo_s2.h5ad
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


def _t1_g_mom(ref, categories: list[str], seed: int, covered_labels: set[str]):
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
    of Feng-covered labels cannot disagree silently."""
    from vec import t1_state_transport as ts

    start = load_stage(DATA / "T1/E9.5_RNA.h5ad", lineage=False)
    factors, _ = harness.g_mom_pair_factors(ref, "E10.5_RNA", start, "E9.5_RNA", categories=categories)
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
        return moved

    return hook, record


def build_t1(a) -> dict:
    if a.g_mom and a.strong_null:
        sys.exit("--g-mom is refused with --strong-null: the two count edits have no ruled order and were never tested "
                 "together; a new row in experiments/harness/TRIES.md and the captain's ruling fix it first")
    if a.g_mom and not a.hold_unbracketed:
        sys.exit("--g-mom is refused with --no-hold-unbracketed: g_mom's factors at test add to M1 with R1 and R2 as "
                 "they are (experiments/harness/TRIES.md, try T11)")
    hygiene_on = _feng_hygiene_on(a.feng_hygiene, a.hold_unbracketed)  # R2's refusals come before anything is read
    prior = _t1_down_prior(a.strong_null, a.hold_unbracketed, hygiene_on)  # T12's, likewise
    union = _t1_union(hygiene_on)
    board, dry = _board("t1", a.board)
    ref = load_stage(DATA / "T1/E10.5_RNA.h5ad", lineage=False)
    momentum = _t1_momentum(a.hold_unbracketed)
    shares = harness.t1_chamber_shares(harness._shallow(ref), "E10.5_RNA")
    fit_categories = [c for c, s in shares.items() if s > 0]  # g_mom's fit, as the candidate's: E10.5's own mix
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
        hook, g_mom_record = _t1_g_mom(ref, fit_categories, a.seed, covered)
        within_type = {"post_cells": hook}
    sub, rep = t1_scope.build_submission(ref, board, shares, n, a.seed, held_flat=held,
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
           "cm_share_of_heart": t1_scope.share_of(shares, "cm", t1_scope.load_spec().heart),
           "blood_share_of_chamber": float(shares.get("blood", 0.0)),
           "shares": shares, "external": ext, "check_ok": check.ok, "realised": rep["realised"]}
    if a.g_mom:
        rec["g_mom"] = g_mom_record
    return rec if strong is None else {**rec, "strong_null": strong}


def build_heart(a) -> dict:
    board, dry = _board("heart", a.board)
    genes = board.genes
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
    view = _view("test:T2h", "T2", "extrap", "heart", board.key,
                 {"ref": ("T2/E10.5.h5ad", 10.5), "start": ("T2/E9.5.h5ad", 9.5)}, "ref", 12.5)
    count = _t2_count(board, pool.n_obs)
    ctx = harness.BuildContext(task=view, board=harness._candidate_board(board), inputs={"ref": pool, "start": start},
                               seed=a.seed)
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
            "size": {k: size[k] for k in size if not isinstance(size[k], (dict, list))}}


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
            p.add_argument("--g-mom", action="store_true",
                           help="add g_mom's within-type factors at test (T1 try T11 of experiments/harness/TRIES.md; "
                                "fitted on E9.5 -> E10.5, applied last, after Feng's direction, on cardiomyocytes and "
                                "endocardium only (R3b), R1's held labels left alone, each cell's expected depth "
                                "kept, no share, cell count or scope changed; refused with --no-hold-unbracketed; "
                                "off: the entry as built, "
                                "and only the captain's ruling turns it on for an upload)")
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
