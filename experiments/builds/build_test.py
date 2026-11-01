#!/usr/bin/env python
"""Test-board candidate files, built with the harness's own builders in the test geometry (P3_RUNBOOK section 5).

Each entry is exactly the candidate the harness judged, run on the test's stages, then the test-only steps the
pseudo-tasks cannot run (the pre-registered sizes). Nothing is uploaded: the file is written with its build stamp and
checked locally (`vec check`); a server format check is a separate `vec submit ... --from-dgx` call from the local
machine, and an official upload needs the captain's go-ahead.

  t1      M1 (docs/PLAN.md 4.1; the captain's rulings of 2026-11-01): E10.5 cells re-weighted to E10.5's own chamber
          mix (or --cm-share, the cardiomyocyte share of heart cells), Feng's E9.5 -> E14.5 direction on the cells of
          Feng's groups at 2/5 of the bracket (E10.5 -> E12.5), B0's damped E9.5 -> E10.5 momentum on the rest
  heart   B0 (copy of E10.5 plus damped E9.5 -> E10.5 momentum in count space), sized by the pre-registered heart
          procedure: --entry s1 = scope hypothesis H1, s2 = H2 (approved rule 3, `size_or_hold`)
  embryo  --method pseudobulk_shift | fgw: E7.5 cells moved by the E7.25 -> E7.5 change, or the FGW interpolation
          between E7.5 and E8.0 at t = 1/2; sized by the pre-registered embryo procedure (size_embryo_E7.75.json)

The board is the test board from boards/ (after the intake, `vec boards update --from-file`); before that, --board
may name a validation board for a dry run (the output name says so, and the file must never be uploaded). The build
seed is never 0 (runbook section 5): default vec.harness.BUILD_SEED_OFFSET + 1.

Usage (DGX, repo root, BLAS threads capped):
  uv run python experiments/builds/build_test.py t1 --out runs/builds/T1_s1.h5ad
  uv run python experiments/builds/build_test.py heart --entry s1 --out runs/builds/T2heart_s1.h5ad
  uv run python experiments/builds/build_test.py embryo --method fgw --out runs/builds/T2embryo_s1.h5ad
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

from vec import harness, t1_scope, t2_build
from vec.boards import REPO_ROOT, load_boards
from vec.io import COORD_KEY, dense
from vec.prep import load_stage

sys.path.insert(0, str(REPO_ROOT / "experiments" / "t1_m1"))
import feng_dirs  # noqa: E402

TEST_BOARDS = {"t1": "T1:test", "heart": "T2:heart:test_extrap", "embryo": "T2:embryo:test_interp"}
SIZES = {"heart": REPO_ROOT / "experiments/coords/results/p3/size_heart_E12.5.json",
         "embryo": REPO_ROOT / "experiments/coords/results/p3/size_embryo_E7.75.json"}
DATA = harness.DATA_DIR


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


def _check_seed(seed: int) -> int:
    if seed < harness.BUILD_SEED_OFFSET:
        sys.exit(f"build seed {seed} is below {harness.BUILD_SEED_OFFSET} (runbook section 5: never seed 0)")
    return seed


def build_t1(a) -> dict:
    board, dry = _board("t1", a.board)
    ref = load_stage(DATA / "T1/E10.5_RNA.h5ad", lineage=False)
    start = load_stage(DATA / "T1/E9.5_RNA.h5ad", lineage=False)
    shares = harness.t1_chamber_shares(harness._shallow(ref), "E10.5_RNA")
    if a.cm_share is not None:
        shares = t1_scope.with_share(shares, "cm", a.cm_share)
    directions, ext = feng_dirs.feng_directions(list(ref.var_names))
    groups = feng_dirs.label_groups(ref.obs["celltype"].astype(str).unique())
    fraction = feng_dirs.bracket_fraction(10.5, 12.5)
    n = harness.count_for(board, ref.n_obs)
    sub, rep = t1_scope.build_submission(ref, board, shares, n, a.seed, start=start, damping=harness.B0_DAMPING,
                                         held_flat=t1_scope.held_flat_genes(), ref_source="E10.5_RNA",
                                         start_source="E9.5_RNA", directions=directions, direction_fraction=fraction,
                                         direction_groups=groups)
    check = t1_scope.write_checked(sub, a.out, board)
    return {"entry": "M1", "board": board.key, "dry_run": dry, "n_cells": n, "seed": a.seed, "fraction": fraction,
            "cm_share_of_heart": t1_scope.share_of(shares, "cm", t1_scope.load_spec().heart),
            "shares": shares, "external": ext, "check_ok": check.ok, "realised": rep["realised"]}


def build_heart(a) -> dict:
    board, dry = _board("heart", a.board)
    genes = board.genes
    ref = load_stage(DATA / "T2/E10.5.h5ad", genes=genes, drop_duplicates=True)
    start = load_stage(DATA / "T2/E9.5.h5ad", genes=genes, drop_duplicates=True)
    view = _view("test:T2h", "T2", "extrap", "heart", board.key,
                 {"ref": ("T2/E10.5.h5ad", 10.5), "start": ("T2/E9.5.h5ad", 9.5)}, "ref", 12.5)
    ctx = harness.BuildContext(task=view, board=harness._candidate_board(board), inputs={"ref": ref, "start": start},
                               seed=a.seed)
    stage = harness._b0_t2(ctx)  # the harness's B0, guard and checks included, in E10.5's frame and size
    hyp = {"s1": "H1", "s2": "H2"}[a.entry]
    coords, size = t2_build.size_or_hold(stage.coords, SIZES["heart"], hyp, source_stage="E10.5", reference_cloud=ref,
                                         expect_board="heart")
    _, report = t2_build.build_submission(board, stage.X, coords, genes=list(genes), out=a.out)
    if not report.ok:
        raise SystemExit(f"the built file fails its check: {report.errors}")
    return {"entry": f"B0 {a.entry} ({hyp})", "board": board.key, "dry_run": dry, "n_cells": int(stage.n),
            "seed": a.seed, "size": {k: size[k] for k in size if not isinstance(size[k], (dict, list))}}


def build_embryo(a) -> dict:
    board, dry = _board("embryo", a.board)
    genes = board.genes
    files = {"before": ("T2/E7.5.h5ad", 7.5), "after": ("T2/E8.0.h5ad", 8.0), "earlier": ("T2/E7.25.h5ad", 7.25)}
    inputs = {r: load_stage(DATA / f, genes=genes, drop_duplicates=True) for r, (f, _) in files.items()}
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
                                          harness.count_for(board, before.n_obs), a.seed, target_rms=target)
        X = np.empty((len(res.rows), len(genes)), dtype=np.float32)
        X[res.stage == 0], X[res.stage == 1] = Eb[res.rows[res.stage == 0]], Ea[res.rows[res.stage == 1]]
        guard = t2_build.cloud_guard(res.coords, res.pool_coords, seed=a.seed, source_weights=res.pool_weights)
        if not guard.ok:
            raise SystemExit(f"cloud guard flags {guard.flags}")
        coords, rec = res.coords, res.record
    _, report = t2_build.build_submission(board, X, coords, genes=list(genes), out=a.out)
    if not report.ok:
        raise SystemExit(f"the built file fails its check: {report.errors}")
    return {"entry": a.method, "board": board.key, "dry_run": dry, "n_cells": int(len(X)), "seed": a.seed,
            "target_rms": target, "record": {k: v for k, v in rec.items() if not isinstance(v, (dict, list))}}


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="which", required=True)
    for name in ("t1", "heart", "embryo"):
        p = sub.add_parser(name)
        p.add_argument("--out", required=True)
        p.add_argument("--board", help="a validation board for a dry run (default: the test board)")
        p.add_argument("--seed", type=int, default=harness.BUILD_SEED_OFFSET + 1)
        if name == "t1":
            p.add_argument("--cm-share", type=float, help="cardiomyocyte share of heart cells (default: E10.5's own)")
        if name == "heart":
            p.add_argument("--entry", choices=("s1", "s2"), required=True)
        if name == "embryo":
            p.add_argument("--method", choices=("pseudobulk_shift", "fgw"), required=True)
    a = ap.parse_args(argv)
    _check_seed(a.seed)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    rec = {"t1": build_t1, "heart": build_heart, "embryo": build_embryo}[a.which](a)
    Path(a.out).with_suffix(".build.json").write_text(json.dumps(rec, indent=1, default=float))
    print(json.dumps({k: v for k, v in rec.items() if k not in ("shares", "external", "realised")}, default=float))


if __name__ == "__main__":
    main()
