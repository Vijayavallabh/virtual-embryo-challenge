#!/usr/bin/env python
"""Target-free occupancy frame check (branch campaign/embryo-occupancy): does a built-form embryo file sit in the
scorer's occupancy frame where a plain copy of the stage whose positions it carries sits?

The scorer's occupancy pipeline, as the emulator runs it (vec.evaluate; the pinned veckit build "main" and the pre-fix
build "prefix" it documents): `_match_n` draws both clouds to a common count at the metric seed (the prediction
first); each cloud is centred, rotated onto its own PCA axes (decreasing variance, signs as LAPACK returns them; "main"
then makes the basis right-handed, "prefix" does not) and divided by its RMS radius; each is binned on a 16^3 grid of
[-3, 3) RMS radii; the prediction's four proper sign flips are tried and the best Dice kept. Translation and size drop
out; axis order, handedness and (under "prefix") the row order do not.

Two modes. Neither reads a target's coordinates or expression, a test anchor or a score, and nothing printed is a cell
count, a size or a score.

  file      a built file given by path (a dry run or a test build) against the plain copy of the released stage whose
            positions it carries (default T2/E7.5.h5ad): vec.evaluate.subsample(stage as released, 0.1, --copy-seed),
            the bake protocol's floor copy (seed 0 on the live board). Prints the MAIN-build frame facts only. The
            pre-fix parity of a test file is never computed here: printed beside that file's returned score it is the
            route to the hidden truth half's parity that docs/PLAN.md decision log 2026-11-18 (k) bars.
  practice  the builder's practice path on released stages: for each embryo pseudo-task of configs/pseudo_tasks.yaml
            and seed, the harness's fgw_shape0x (the second upload's position form, vec.harness) at its build seed,
            against the bake copy of the task's reference stage at that seed. Prints the main-build facts and the
            pre-fix facts (each cloud's parity on the rows `_match_n` would take against the task's truth half, read
            from the target file's cell count only; whether the two are handed the same truth rows), for the file as
            built and with `--row-order copy_parity` applied (vec.t2_build.copy_parity_order).

Main-build facts (file vs copy, on the rows `_match_n` hands over at --metric-seed, the file as the prediction):
  axes_agree        each of the file's three right-handed PCA axes is closest to the copy's axis of the same rank and
                    every |cos| >= AXIS_COS_MIN (a swap of near-equal axes, mechanism 3 of the brief, shows here)
  axis_cos          |cos| between the axes of the same rank
  axis_signs        sign of each cos: the proper flip the scorer's search must use (product +1 under main)
  overlap           best proper-flip Dice of the two occupancy sets: the main build's occupancy_dice(file, copy)
  overlap_mirror    the same for the file's z-mirror: what a mirror-mode comparison would give (the coin's stake here)
  handedness_agree  overlap >= overlap_mirror
  ref_overlap_lo/mean/sd  overlap of REF_DRAWS other bake-style copies (seeds REF_SEED_BASE + k) against the copy: the
                    spread of a plain copy's own draw-to-draw variance (the challenger's sampling mechanism)
  coincide          axes_agree and handedness_agree and overlap >= ref_overlap_mean - REF_SD_MULT x ref_overlap_sd
Pre-fix facts (practice only): parity_file, parity_copy, parity_agree, truth_rows_shared, same_mode_whatever_truth
(shared and agree: the file then falls in the copy's mode whichever parity the truth half has), overlap_prefix (the
pre-fix occupancy_dice(file, copy): overlap or overlap_mirror by the two clouds' parities on their matched rows).

Usage (DGX campaign checkout, repo root; BLAS threads capped; heavy only through the practice builds):
  uv run --frozen --no-sync python -I -B experiments/coords/occupancy_frame_check.py file --file <built.h5ad>
  uv run --frozen --no-sync python -I -B experiments/coords/occupancy_frame_check.py practice --seeds 0 1 2
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from vec import evaluate as ev  # noqa: E402
from vec import t2_build as tb  # noqa: E402

AXIS_COS_MIN = 0.95  # fixed before any output: an axis turned by more than about 18 degrees does not agree
REF_DRAWS = 10
REF_SD_MULT = 3.0  # fixed before any output: a file below the copies' mean by more than 3 of their sd does not coincide
REF_SEED_BASE = 4900  # seeds of the reference copies: apart from the copy's own seed and the practice seeds 0-9
COPY_STAGE = "T2/E7.5.h5ad"
MAIN_FIELDS = ("axes_agree", "axis_cos", "axis_signs", "overlap", "overlap_mirror", "handedness_agree",
               "ref_overlap_lo", "ref_overlap_mean", "ref_overlap_sd", "coincide")
WITHHELD = ("withheld for a file given by path: its pre-fix parity beside its returned score would solve the hidden "
            "truth half's parity (docs/PLAN.md decision log 2026-11-18 (k))")


# ------------------------------------------------------------------------------------------------ the scorer's frame
def basis(C: np.ndarray, build: str = "main") -> np.ndarray:
    """The PCA basis the scorer builds for these rows (rows of the result, by decreasing variance), with the array
    operations of vec.evaluate._canonicalise_prefix; "main" adds the right-handed flip of the third axis."""
    X = np.asarray(C, dtype=np.float64)
    X = X - X.mean(0)
    _, _, Vt = np.linalg.svd(X - X.mean(0), full_matrices=False)
    if build == "main" and np.linalg.det(Vt) < 0:
        Vt[-1] *= -1
    return Vt


def _voxels(C: np.ndarray, V: np.ndarray) -> np.ndarray:
    return tb._voxel_ids(tb._canonical(np.asarray(C, dtype=np.float64), V))


def best_dice(P: np.ndarray, Vp: np.ndarray, Q: np.ndarray, Vq: np.ndarray) -> float:
    """occupancy_dice's comparison of matched rows P (prediction) and Q (truth), each in its own frame: the best Dice
    over the prediction's proper sign flips, rounded as the panel rounds it."""
    vq = _voxels(Q, Vq)
    best = 0.0
    for f in tb.PROPER_FLIPS:
        vp = _voxels(P, Vp * np.asarray(f, dtype=np.float64)[:, None])
        best = max(best, 2.0 * np.intersect1d(vp, vq, assume_unique=True).size / (vp.size + vq.size))
    return round(float(best), 4)


def _read(C) -> np.ndarray:
    return tb._as_read(C)


def frame_facts(file_coords, copy_coords, refs=(), metric_seed: int = 0) -> dict:
    """The main-build facts (module docstring) of a file against a plain copy; refs: other copies of the same stage."""
    F, Q = _read(file_coords), _read(copy_coords)
    A, B = ev.match_rows(F, Q, metric_seed)
    Va, Vb = basis(A), basis(B)
    M = Va @ Vb.T
    absM = np.abs(M)
    order_ok = bool((absM.argmax(1) == np.arange(3)).all() and (absM.argmax(0) == np.arange(3)).all())
    cos = [round(float(c), 4) for c in np.abs(np.diag(M))]
    overlap = best_dice(A, Va, B, Vb)
    Am = ev.mirror_z(A)
    overlap_m = best_dice(Am, basis(Am), B, Vb)
    ref = []
    for R in refs:
        Ra, Rb = ev.match_rows(_read(R), Q, metric_seed)
        ref.append(best_dice(Ra, basis(Ra), Rb, basis(Rb)))
    mean = float(np.mean(ref)) if ref else None
    sd = float(np.std(ref, ddof=1)) if len(ref) > 1 else 0.0
    axes_agree = order_ok and min(cos) >= AXIS_COS_MIN
    hand = overlap >= overlap_m
    return {"axes_agree": axes_agree, "axis_cos": cos, "axis_signs": [int(s) for s in np.sign(np.diag(M))],
            "overlap": overlap, "overlap_mirror": overlap_m, "handedness_agree": bool(hand),
            "ref_overlap_lo": min(ref) if ref else None, "ref_overlap_mean": None if mean is None else round(mean, 4),
            "ref_overlap_sd": round(sd, 4) if ref else None,
            "coincide": bool(axes_agree and hand and (mean is None or overlap >= mean - REF_SD_MULT * sd))}


def prefix_facts(file_coords, copy_coords, n_true: int, metric_seed: int = 0, overlap: bool = False) -> dict:
    """The pre-fix facts (module docstring) for released stages only. n_true: the truth half's cell count, read by
    code; it is not printed."""
    F, Q = _read(file_coords), _read(copy_coords)
    pf, pc = ev.prefix_parity(F, n_true, metric_seed), ev.prefix_parity(Q, n_true, metric_seed)
    shared = tb.truth_rows_shared(len(F), len(Q), n_true, metric_seed)
    out = {"parity_file": pf, "parity_copy": pc, "parity_agree": pf == pc, "truth_rows_shared": shared,
           "same_mode_whatever_truth": bool(shared and pf == pc)}
    if overlap:
        A, B = ev.match_rows(F, Q, metric_seed)
        Va, Vb = basis(A, "prefix"), basis(B, "prefix")
        out["overlap_prefix"] = best_dice(A, Va, B, Vb)
    return out


# ------------------------------------------------------------------------------------------------ inputs
def _coords_only(C: np.ndarray) -> ev.Stage:
    return ev.Stage(X=np.zeros((len(C), 0), dtype=np.float32), coords=np.asarray(C, dtype=np.float32)[:, :3])


def released_copies(stage_path: Path, copy_seed: int, n_ref: int) -> tuple[np.ndarray, list[np.ndarray]]:
    """The bake copy of a released stage (no duplicate drop) and n_ref other copies drawn the same way."""
    import anndata as ad

    from vec.io import COORD_KEY

    a = ad.read_h5ad(stage_path)
    st = _coords_only(np.asarray(a.obsm[COORD_KEY], dtype=np.float32))
    copy = ev.subsample(st, ev.SUBSAMPLE_FRAC, copy_seed).coords
    return copy, [ev.subsample(st, ev.SUBSAMPLE_FRAC, REF_SEED_BASE + k).coords for k in range(n_ref)]


def _cells_only(path: Path) -> int:
    """A file's cell count from its obs index alone (no coordinates or expression read)."""
    import h5py

    with h5py.File(path, "r") as f:
        obs = f["obs"]
        key = obs.attrs.get("_index", "_index")
        return int(obs[key].shape[0])


def run_file(a) -> None:
    import anndata as ad

    from vec.io import COORD_KEY

    data = Path(a.data_dir)
    F = np.asarray(ad.read_h5ad(a.file).obsm[COORD_KEY], dtype=np.float32)[:, :3]
    copy, refs = released_copies(data / a.stage, a.copy_seed, a.n_ref)
    out = {"mode": "file", "file": Path(a.file).name, "copy": f"subsample({a.stage} as released, 0.1, {a.copy_seed})",
           "metric_seed": a.metric_seed, "build": "main", **frame_facts(F, copy, refs, a.metric_seed),
           "prefix": WITHHELD}
    print(json.dumps(out))


def run_practice(a) -> None:
    from vec import harness
    from vec.boards import load_boards
    from vec.io import read_stage
    from vec.prep import load_stage

    data = Path(a.data_dir)
    boards = load_boards()
    tasks = {t.key: t for ts in harness.load_tasks().values() for t in ts}
    keys = a.task or [k for k, t in tasks.items() if t.task == "T2" and t.setting == "embryo" and t.kind == "interp"]
    cand = harness.CANDIDATES["fgw_shape0x"]
    for key in keys:
        t = tasks[key]
        if not cand.applies(t):
            raise SystemExit(f"{key}: fgw_shape0x does not apply")
        genes = boards[t.panel_board].genes
        ref = ev.Stage.from_anndata(read_stage(data / t.reference, genes), genes, need_coords=True)
        board = harness.task_board(t, boards, _cells_only(data / t.target), ref.n)
        inputs = {r: load_stage(data / f, genes=genes, drop_duplicates=True) for r, f in t.inputs.items()}
        view, cboard = harness.TaskView.of(t), harness._candidate_board(board)
        refs = [ev.subsample(ref, ev.SUBSAMPLE_FRAC, REF_SEED_BASE + k).coords for k in range(a.n_ref)]
        for s in a.seeds:
            ctx = harness.BuildContext(task=view, board=cboard, inputs=inputs, seed=harness.BUILD_SEED_OFFSET + s)
            pred = cand.build(ctx)
            copy = ev.subsample(ref, ev.SUBSAMPLE_FRAC, s).coords  # the bake protocol's floor copy at seed s
            arms = {"as_built": pred.coords}
            try:
                order, _ = tb.copy_parity_order(pred.coords, copy, board.truth_cells, ctx.seed, metric_seed=s)
                arms["copy_parity"] = pred.coords[order]
            except ValueError as e:
                arms["copy_parity"] = None
                refusal = str(e)
            for arm, C in arms.items():
                row = {"mode": "practice", "task": key, "seed": s, "candidate": "fgw_shape0x", "row_order": arm,
                       "metric_seed": s}
                if C is None:
                    row["refused"] = refusal
                else:
                    row.update(frame_facts(C, copy, refs, s))
                    row["prefix"] = prefix_facts(C, copy, board.truth_cells, s, overlap=True)
                print(json.dumps(row))


def build_parser() -> argparse.ArgumentParser:
    from vec.harness import DATA_DIR

    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="mode", required=True)
    f = sub.add_parser("file", help="a built file against the plain copy of its stage (main-build facts only)")
    f.add_argument("--file", required=True)
    f.add_argument("--stage", default=COPY_STAGE, help="the released stage whose positions the file carries")
    f.add_argument("--copy-seed", type=int, default=0)
    p = sub.add_parser("practice", help="the builder's practice path on released stages")
    p.add_argument("--task", action="append", help="a pseudo-task key (default: every embryo interpolation task)")
    p.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2])
    for q in (f, p):
        q.add_argument("--data-dir", default=str(DATA_DIR))
        q.add_argument("--n-ref", type=int, default=REF_DRAWS)
    f.add_argument("--metric-seed", type=int, default=0)
    return ap


def main(argv=None) -> None:
    a = build_parser().parse_args(argv)
    {"file": run_file, "practice": run_practice}[a.mode](a)


if __name__ == "__main__":
    main()
