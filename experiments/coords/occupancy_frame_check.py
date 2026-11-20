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
            (or, with --nine, of the nine practice tasks of experiments/harness/nine_tasks.py) and seed, the
            harness's fgw_shape0x (the second upload's position form; its series twin has the same rows and
            positions) at its build seed, against the bake copy of the task's reference stage at that seed. Prints
            the main-build facts and the pre-fix facts (each cloud's parity on the rows `_match_n` would take
            against the task's truth half, read from the target file's cell count only; whether the two are handed
            the same truth rows), for the file as built and with `--row-order copy_parity` applied
            (vec.t2_build.copy_parity_order).

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

A third mode, added after the `file` output of 2026-11-20 fell under reading A.2 (below), before any output of its own:

  file_axes the same file and the same copy (metric seed 0, main build): whether plain copies of the stage turn their
            second and third axes as far as the file does. Released stage and the built file only; no score, no
            target, no count printed; no pre-fix parity of the file; whether the floor_exact fallback fires is not
            computed. Fields:
              ref_cos              per reference copy R (AXES_REF_DRAWS copies at seeds REF_SEED_BASE + k), the three
                                   diagonal |cos| of basis(match_rows(R, copy)) against the copy's, as frame_facts does
              ref_axes_agree_n     how many of them pass frame_facts' axes_agree test (axis order and AXIS_COS_MIN)
              file_axis_cos, file_min23   the file's |cos| and min(cos2, cos3)
              file_rank23          how many copies have min(cos2, cos3) <= file_min23 (0 to AXES_REF_DRAWS)
              ref_min23_3rd_lowest the third-lowest copy value of min(cos2, cos3)
              sv_ratio_copy, sv_ratio_stage   [s2/s1, s3/s2] of the copy's matched rows and of the whole stage
              file_min23_by_metric_seed   file_min23 at metric seeds AXES_METRIC_SEEDS
              dedup_draw_min23     min(cos2, cos3) of DEDUP_DRAWS uniform draws without replacement from the stage's
                                   de-duplicated positions (load_stage(..., drop_duplicates=True)) at the file's own
                                   count, or all of them if the file holds more (used in code, never printed; seeds
                                   DEDUP_SEED_BASE + k), each through match_rows against the copy; and its median
              turn_cost, ref_turn_cost   best proper-flip Dice of the matched rows canonicalised in the copy's basis,
                                   minus the same in their own basis (the scorer's value): the grid's price of the turn
                                   (positive: the turn costs Dice); the same for each reference copy
              overlap, ref_overlap_p2   the file's main-build overlap and the copies' 2nd percentile of theirs
            Readings, fixed before any output (the analysts' hand-back of 2026-11-20, as the challenger changed it):
              S1  file_rank23 >= S1_MIN_RANK: stage noise; every copy carries a turn like it; A.2 fired on an
                  uncalibrated threshold; Option A is not written.
              S2  file_rank23 < S1_MIN_RANK and the median of dedup_draw_min23 <= ref_min23_3rd_lowest: the turn comes
                  from de-duplication or count; a property of the file's position source.
              S3  file_rank23 < S1_MIN_RANK and that median above ref_min23_3rd_lowest: the turn comes from the file's
                  own rows; file_min23_by_metric_seed is read by the fixed text (the draw if it reaches the copies'
                  range at most seeds).
              option_a_eligible: S2 or S3, and a cost on the grid: turn_cost above the copies' TURN_PCT[1]th
                  percentile or overlap below their 2nd percentile, at metric seed 0. (The challenger wrote "turn_cost
                  below the 2nd percentile"; with turn_cost defined as above, a cost is the high side, so the code
                  reads the cost side; both turn_cost booleans are printed.) Without that cost: the turn is the
                  file's but costs nothing on the grid against the copy; Option A is not written.
            No reading says anything about the hidden truth's frame, and none is set beside a returned score.

Readings fixed before any output of `file` and `practice` (the analysts' hand-back of 2026-11-20; B.2 as the reviewer
reworded it, A.7 and A.8 added by the reviewer):
  A.1 coincide true: under main the file sits in the copy's frame; the one construction is --row-order copy_parity,
      a priced question after the practice read.
  A.2 axes_agree false, handedness_agree true: the file's draw orders or turns its principal axes differently from
      the copy's: mechanism 3 is active. The scorer searches signs only, so this costs occupancy under both builds.
      What follows: report it. The row-order option does not fix it. The only fix the code allows is giving the file
      the copy's own rows (the investigator's Option A), which breaks the count rule and needs a ruling. Nothing is
      built until he rules.
  A.3 handedness_agree false: a builder defect; stop, find it on released data, upload nothing on this term.
  A.4 axes and handedness agree, overlap below ref_overlap_mean - 3 sd: duplicate drop, jitter or the fallback; check
      by code only whether the fallback fires; report the shortfall; no construction on this branch.
  A.5 overlap - overlap_mirror within one ref_overlap_sd: the coin cannot be large; the option is not worth a ruling.
  A.6 first-upload file: coincide false is expected; a contrast only.
  A.7 the dry build coincides but the DGX test build does not, or the reverse: the difference is the count or the
      board, not the method; report it, and attribute nothing until the same code path is run at both counts on
      released data.
  A.8 the local second-upload file and the DGX test build of the same command give different main facts: the uploaded
      bytes are not this build; stop and reconcile before any reading.
  B.1 as_built rows: parity_agree true on about half the seeds where truth_rows_shared; far from half: report it.
  B.2 copy_parity rows: same_mode_whatever_truth must be true wherever truth_rows_shared; this recomputes the option's
      own gate and checks only the float32 path. Whether the mode really matches the floor's under the pre-fix scorer
      is shown only by E48's scoring against the fold truth.
  B.3 refused rows or truth_rows_shared false: those tasks cannot price the option; if every task refuses, reported.
  B.4 main facts should have the same spread in both row orders; beyond the copies' sd: report it.
  B.5 overlap versus overlap_mirror per task: the coin's stake on released geometries.
None of these says which mode or frame the hidden truth half has; none may be set beside a returned test score.

Usage (DGX campaign checkout, repo root; BLAS threads capped; heavy only through the practice builds):
  uv run --frozen --no-sync python -I -B experiments/coords/occupancy_frame_check.py file --file <built.h5ad>
  uv run --frozen --no-sync python -I -B experiments/coords/occupancy_frame_check.py file_axes --file <built.h5ad>
  uv run --frozen --no-sync python -I -B experiments/coords/occupancy_frame_check.py practice --nine --seeds 0 1 2
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
MIN_REF = 2  # coincide needs a spread (reviewer point 4)
# file_axes, fixed before its output (the analysts' hand-back of 2026-11-20)
AXES_REF_DRAWS = 50  # reference copies at seeds REF_SEED_BASE + k, k < 50 (4900-4949)
DEDUP_DRAWS = 20
DEDUP_SEED_BASE = 4950  # 4950-4969
AXES_METRIC_SEEDS = tuple(range(10))
S1_MIN_RANK = 3  # at least 3 of the 50 copies (6%) turned as far or further: stage noise
TURN_PCT = (2.0, 98.0)
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


def _axes(Va: np.ndarray, Vb: np.ndarray) -> tuple[bool, list[float]]:
    """Whether each axis of basis Va is nearest Vb's axis of the same rank (by row and by column), and the three
    diagonal |cos|, rounded to 4 decimals."""
    absM = np.abs(Va @ Vb.T)
    order_ok = bool((absM.argmax(1) == np.arange(3)).all() and (absM.argmax(0) == np.arange(3)).all())
    return order_ok, [round(float(c), 4) for c in np.diag(absM)]


def frame_facts(file_coords, copy_coords, refs=(), metric_seed: int = 0) -> dict:
    """The main-build facts (module docstring) of a file against a plain copy; refs: other copies of the same stage."""
    F, Q = _read(file_coords), _read(copy_coords)
    A, B = ev.match_rows(F, Q, metric_seed)
    Va, Vb = basis(A), basis(B)
    M = Va @ Vb.T
    order_ok, cos = _axes(Va, Vb)
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


def _min23(cos) -> float:
    return float(min(cos[1], cos[2]))


def turn_cost(A: np.ndarray, B: np.ndarray) -> float:
    """The grid's price of A's turn against B (matched rows): best proper-flip Dice with A canonicalised in B's basis
    minus the same in A's own basis (what the scorer compares). Positive: the turn costs Dice."""
    Va, Vb = basis(A), basis(B)
    return round(best_dice(A, Vb, B, Vb) - best_dice(A, Va, B, Vb), 4)


def _sv_ratio(C: np.ndarray) -> list[float]:
    X = np.asarray(C, dtype=np.float64)
    X = X - X.mean(0)
    s = np.linalg.svd(X, compute_uv=False)
    return [round(float(s[1] / s[0]), 4), round(float(s[2] / s[1]), 4)]


def file_axes_facts(file_coords, copy_coords, stage_coords, refs, dedup_draws, metric_seed: int = 0,
                    metric_seeds=AXES_METRIC_SEEDS) -> dict:
    """The `file_axes` facts and readings (module docstring). refs: the reference copies (AXES_REF_DRAWS of them for
    the fixed readings); dedup_draws: the de-duplicated stage draws. Nothing returned is a count, a size or a score."""
    if len(refs) < S1_MIN_RANK:
        raise ValueError(f"file_axes needs at least {S1_MIN_RANK} reference copies")
    F, Q = _read(file_coords), _read(copy_coords)
    A, B = ev.match_rows(F, Q, metric_seed)
    order_ok, cos = _axes(basis(A), basis(B))
    overlap = best_dice(A, basis(A), B, basis(B))
    tc = turn_cost(A, B)
    ref_cos, ref_ok, ref_overlap, ref_tc = [], 0, [], []
    for R in refs:
        Ra, Rb = ev.match_rows(_read(R), Q, metric_seed)
        ok, c = _axes(basis(Ra), basis(Rb))
        ref_cos.append(c)
        ref_ok += int(ok and min(c) >= AXIS_COS_MIN)
        ref_overlap.append(best_dice(Ra, basis(Ra), Rb, basis(Rb)))
        ref_tc.append(turn_cost(Ra, Rb))
    f23 = _min23(cos)
    r23 = sorted(_min23(c) for c in ref_cos)
    rank = int(sum(v <= f23 for v in r23))
    third = r23[S1_MIN_RANK - 1]
    by_seed = []
    for ms in metric_seeds:
        As, Bs = ev.match_rows(F, Q, ms)
        by_seed.append(_min23(_axes(basis(As), basis(Bs))[1]))
    dd = []
    for D in dedup_draws:
        Da, Db = ev.match_rows(_read(D), Q, metric_seed)
        dd.append(_min23(_axes(basis(Da), basis(Db))[1]))
    med = float(np.median(dd)) if dd else None
    if rank >= S1_MIN_RANK:
        reading = "S1"
    elif med is None:
        reading = "S2 or S3 (no de-duplicated draws)"
    else:
        reading = "S2" if med <= third else "S3"
    tc_lo, tc_hi = (float(v) for v in np.percentile(ref_tc, TURN_PCT))
    ov_p2 = float(np.percentile(ref_overlap, TURN_PCT[0]))
    cost = bool(tc > tc_hi or overlap < ov_p2)
    return {"file_axis_cos": cos, "file_axes_order_ok": order_ok, "file_min23": f23, "file_rank23": rank,
            "ref_cos": ref_cos, "ref_axes_agree_n": ref_ok, "ref_min23_3rd_lowest": third,
            "sv_ratio_copy": _sv_ratio(B), "sv_ratio_stage": _sv_ratio(_read(stage_coords)),
            "file_min23_by_metric_seed": by_seed,
            "file_min23_seeds_above_ref_3rd_lowest": int(sum(v > third for v in by_seed)),
            "dedup_draw_min23": dd, "dedup_draw_min23_median": None if med is None else round(med, 4),
            "turn_cost": tc, "ref_turn_cost": ref_tc, "ref_turn_cost_p2": round(tc_lo, 4),
            "ref_turn_cost_p98": round(tc_hi, 4), "turn_cost_within_2_98": bool(tc_lo <= tc <= tc_hi),
            "turn_cost_below_ref_p2": bool(tc < tc_lo), "turn_cost_above_ref_p98": bool(tc > tc_hi),
            "overlap": overlap, "ref_overlap_p2": round(ov_p2, 4), "overlap_below_ref_p2": bool(overlap < ov_p2),
            "reading": reading, "grid_cost": cost, "option_a_eligible": bool(reading in ("S2", "S3") and cost)}


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


def released_axes_inputs(stage_path: Path, copy_seed: int, n_file: int):
    """For file_axes: the bake copy of the released stage, the stage's own coordinates, AXES_REF_DRAWS reference
    copies (seeds REF_SEED_BASE + k) and DEDUP_DRAWS uniform draws without replacement (seeds DEDUP_SEED_BASE + k)
    from the de-duplicated stage at n_file rows, or all of them if it holds fewer (n_file is used here only)."""
    import anndata as ad

    from vec.io import COORD_KEY
    from vec.prep import load_stage

    a = ad.read_h5ad(stage_path)
    st = _coords_only(np.asarray(a.obsm[COORD_KEY], dtype=np.float32))
    copy = ev.subsample(st, ev.SUBSAMPLE_FRAC, copy_seed).coords
    refs = [ev.subsample(st, ev.SUBSAMPLE_FRAC, REF_SEED_BASE + k).coords for k in range(AXES_REF_DRAWS)]
    dedup = np.asarray(load_stage(stage_path, drop_duplicates=True).obsm[COORD_KEY], dtype=np.float32)[:, :3]
    n = min(int(n_file), len(dedup))
    draws = [dedup[np.random.default_rng(DEDUP_SEED_BASE + k).choice(len(dedup), n, replace=False)]
             for k in range(DEDUP_DRAWS)]
    return copy, st.coords, refs, draws


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


def run_file_axes(a) -> None:
    import anndata as ad

    from vec.io import COORD_KEY

    data = Path(a.data_dir)
    F = np.asarray(ad.read_h5ad(a.file).obsm[COORD_KEY], dtype=np.float32)[:, :3]
    copy, stage, refs, draws = released_axes_inputs(data / a.stage, a.copy_seed, len(F))
    out = {"mode": "file_axes", "file": Path(a.file).name,
           "copy": f"subsample({a.stage} as released, 0.1, {a.copy_seed})", "metric_seed": 0, "build": "main",
           "ref_seeds": f"{REF_SEED_BASE}-{REF_SEED_BASE + AXES_REF_DRAWS - 1}",
           "dedup_seeds": f"{DEDUP_SEED_BASE}-{DEDUP_SEED_BASE + DEDUP_DRAWS - 1}",
           **file_axes_facts(F, copy, stage, refs, draws, 0), "prefix": WITHHELD,
           "fallback": "not computed (a one-bit function of a test count)"}
    print(json.dumps(out))


def run_practice(a) -> None:
    from vec import harness
    from vec.boards import load_boards
    from vec.io import read_stage
    from vec.prep import load_stage

    data = Path(a.data_dir)
    boards = load_boards()
    if a.nine:  # the nine embryo practice tasks of tries E28 and E47 (experiments/harness/nine_tasks.py), by short name
        sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "harness"))
        import nine_tasks

        tasks = dict(nine_tasks.all_tasks())
        keys = a.task or list(tasks)
    else:
        tasks = {t.key: t for ts in harness.load_tasks().values() for t in ts}
        keys = a.task or [k for k, t in tasks.items()
                          if t.task == "T2" and t.setting == "embryo" and t.kind == "interp"]
    cand = harness.CANDIDATES["fgw_shape0x"]
    for key in keys:
        t = tasks[key]
        if getattr(t, "target_filter", None) is not None:
            raise SystemExit(f"{key}: a within-type task; the check does not cut stages as harness.run_task does")
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
    p.add_argument("--task", action="append",
                   help="a pseudo-task key, or with --nine a short name (default: every task of the set)")
    p.add_argument("--nine", action="store_true",
                   help="the nine embryo practice tasks of experiments/harness/nine_tasks.py (half, late, new1, new2, "
                        "e75, e725, x1, x2, x3) in place of configs/pseudo_tasks.yaml's embryo tasks")
    p.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2])
    for q in (f, p):
        q.add_argument("--data-dir", default=str(DATA_DIR))
        q.add_argument("--n-ref", type=int, default=REF_DRAWS)
    f.add_argument("--metric-seed", type=int, default=0)
    x = sub.add_parser("file_axes", help="whether plain copies turn their axes as far as the file (main build, metric "
                                         "seed 0, fixed seeds and readings S1-S3)")
    x.add_argument("--file", required=True)
    x.add_argument("--stage", default=COPY_STAGE, help="the released stage whose positions the file carries")
    x.add_argument("--copy-seed", type=int, default=0)
    x.add_argument("--data-dir", default=str(DATA_DIR))
    return ap


def main(argv=None) -> None:
    a = build_parser().parse_args(argv)
    if getattr(a, "n_ref", MIN_REF) < MIN_REF:
        raise SystemExit(f"--n-ref must be at least {MIN_REF} (coincide needs a spread)")
    {"file": run_file, "file_axes": run_file_axes, "practice": run_practice}[a.mode](a)


if __name__ == "__main__":
    main()
