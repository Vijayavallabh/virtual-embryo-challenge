"""Density experiment: how submitted cell count, sampling pattern and coordinate jitter move the ranked T2 metrics.

A released stage plays the held-out target. score_board builds its usual working copy (10% of the target, seed s),
splits it into truth and ceiling halves, and uses the preceding released stage as the DE reference and copy_last
floor. Predictions are drawn from the target cells that are NOT in that working copy, so a prediction differs from
the truth only in how many cells it has, how they were picked and whether their coordinates or expression rows were
perturbed.

Count levels, in two bases (n is capped at the pool size):
  <r>x     r x truth_cells, r in RATIOS. truth_cells = round(0.1 N) // 2 is computed from the target, never read from
           boards/index.json.
  <f>max   f x max_cells, f in MAX_FRACS. For a release target (the board's own held-out stage) max_cells is the
           board's published cap; for a pseudo-target it is the documented cap rule
           max(CAP_FLOOR, CAP_MULT x max(truth_cells, ref_cells)) from docs/CHALLENGE.md. --max-cells overrides both.

Conditions per (target, seed):
  random    random subsets of the pool at every level of both bases (nested: one permutation, prefixes)
  shuffle50 the same random subsets with 50% of their expression rows permuted among themselves: coordinates, the
            cell multiset and the pseudobulk are kept, local organisation is damaged. An intermediate-quality
            prediction, from the independent review (review/review_run.py; same RNG recipe on the <r>x levels, so
            those rows reproduce the review's)
  block     a contiguous run of pool rows (rows are stored in z-section order), one centre per count level, drawn
            uniformly over the positions where that block fits (BLOCK_DRAW; the first version nested every block
            around one centre chosen so the 8x block fit, which pinned the reference blocks to one slab)
  head      the first n rows of the pool in file order (what a row cap such as a[:n] submits)
  jitter    the random subsets at 1x and 4x with isotropic Gaussian jitter of RMS length f x median NN distance,
            f in {0.25, 0.5, 1, 2}; unit "own" = the submitted cloud's own median NN distance, "truth" = the truth
            half's (a fixed physical length, so 4x/"truth" and 1x/"own" move cells by about the same distance)
  refstage  random subsets of the REFERENCE stage (cells outside the reference working copy, i.e. outside the
            floor) at every level of both bases, coordinates untouched: the count sweep for a floor-quality
            prediction
  refblock / refhead  contiguous row blocks (per-level uniform centre / first rows) of that reference pool: the
            contiguous-rows penalty for a floor-quality prediction

Usage (DGX, repo root):
  uv run python experiments/density/run_density.py --target heart_E8.75 --seeds 0 1 2 3 4 --out runs/density/raw_v2
  uv run python experiments/density/run_density.py --target heart_E9.5 --seeds 0 --patterns block refblock
Rows are appended to <out>/<target>_s<seed>.jsonl as each call finishes; finished conditions are skipped on rerun.
Release targets whose files are not there yet are skipped with a message (exit code 0).
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
from sklearn.neighbors import NearestNeighbors

from vec.boards import RANKED_METRICS, get_board
from vec.evaluate import SUBSAMPLE_FRAC, Stage, score_board, split_half, subsample
from vec.prep import load_stage

TARGETS = {
    # Pseudo-targets: released training stages standing in for a held-out stage. The board supplies only its gene
    # panel; its published max_cells and truth_cells belong to its own held-out stage and are not used here.
    "heart_E9.5": dict(board="T2:heart:val_extrap", target="data/raw/T2/E9.5.h5ad", ref="data/raw/T2/E8.75.h5ad"),
    "heart_E8.75": dict(board="T2:heart:val_extrap", target="data/raw/T2/E8.75.h5ad",
                        ref="data/raw/T2/E8.25_late.h5ad"),
    "embryo_E8.0": dict(board="T2:embryo:val_interp", target="data/raw/T2/E8.0.h5ad", ref="data/raw/T2/E7.25.h5ad"),
    # Release targets: the validation answers (released 2026-10-20) ARE these boards' held-out stages, so the board's
    # own max_cells applies. Skipped with a message until the files exist.
    "heart_E10.5": dict(board="T2:heart:val_extrap", target="data/raw/T2/E10.5.h5ad", ref="data/raw/T2/E9.5.h5ad",
                        release=True),
    "embryo_E7.5": dict(board="T2:embryo:val_interp", target="data/raw/T2/E7.5.h5ad", ref="data/raw/T2/E7.25.h5ad",
                        release=True),
}
RATIOS = [0.25, 0.5, 1, 2, 4, 8]
MAX_FRACS = [0.25, 0.5, 0.75, 1.0]
CAP_MULT, CAP_FLOOR = 3, 5000
BLOCK_DRAW = "per_level_uniform"
SHUFFLE_FRAC = 0.5
# (ratio, jitter fraction, unit of the median NN distance)
JITTERS = [(r, f, unit) for r, unit in ((1, "own"), (4, "own"), (4, "truth")) for f in (0.25, 0.5, 1.0, 2.0)]
EXTRA_METRICS = ["energy_distance", "count_log_ratio", "morans_I_agreement", "pseudobulk_pearson", "pb_rel_err",
                 "variance_ratio", "sliced_wasserstein", "_dice_voxel_over_nn", "_n_up", "_n_dn", "_de_raw",
                 "_de_chance"]


def emulator_indices(n_target: int, n_ref: int, seed: int, frac: float = SUBSAMPLE_FRAC):
    """Row indices of score_board's target working copy, its truth/ceiling halves and the reference working copy,
    obtained by running vec.evaluate.subsample and split_half (same seeds as score_board) on columns of row
    indices."""
    rows = Stage(X=np.arange(n_target, dtype=np.int64)[:, None])
    wc = subsample(rows, frac, seed)
    truth, other = split_half(wc, seed=0)
    ref_wc = subsample(Stage(X=np.arange(n_ref, dtype=np.int64)[:, None]), frac, seed + 1)
    return wc.X[:, 0], truth.X[:, 0], other.X[:, 0], ref_wc.X[:, 0]


def rule_max_cells(truth_cells: int, ref_cells: int) -> int:
    """The documented cap rule (docs/CHALLENGE.md): 3 x max(truth half, reference working copy), at least 5,000."""
    return max(CAP_FLOOR, CAP_MULT * max(truth_cells, ref_cells))


def median_nn(coords: np.ndarray) -> float:
    d, _ = NearestNeighbors(n_neighbors=2).fit(coords).kneighbors(coords)
    return float(np.median(d[:, 1]))


def count_levels(truth_cells: int, max_cells: int, n_cap: int) -> list[tuple[str, float, str, int, float]]:
    """(basis, level, label suffix, n, nominal n / truth_cells) for every count level, n capped at n_cap."""
    out = [("truth", r, f"{r}x", min(int(round(r * truth_cells)), n_cap), float(r)) for r in RATIOS]
    out += [("max", f, f"{f:g}max", min(int(round(f * max_cells)), n_cap), f * max_cells / truth_cells)
            for f in MAX_FRACS]
    return out


def block_start(seed: int, stream: int, k: int, n_pool: int, n: int) -> int:
    """Start of an n-row block drawn uniformly over the n_pool - n + 1 positions where it fits."""
    return int(np.random.default_rng([seed, stream, k]).integers(0, n_pool - n + 1))


def conditions(pool: np.ndarray, ref_pool: np.ndarray, truth_cells: int, max_cells: int, seed: int):
    """Yield (label, meta, source stage, row indices, perturbation or None). Indices are rows of the full source
    stage ("target" or "ref"); the perturbation is ("jitter", f, unit, k) or ("shuffle", frac, rng index)."""
    n_pool, n_ref = len(pool), len(ref_pool)
    tgt = count_levels(truth_cells, max_cells, n_pool)
    ref = count_levels(truth_cells, max_cells, n_ref)

    def meta(pattern, basis, level, ratio, **kw):
        return dict(pattern=pattern, basis=basis, level=level, ratio=ratio, **kw)

    perm = np.random.default_rng([seed, 1]).permutation(n_pool)
    rand = {sfx: pool[np.sort(perm[:n])] for _, _, sfx, n, _ in tgt}
    for basis, lv, sfx, _, ratio in tgt:
        yield f"random_{sfx}", meta("random", basis, lv, ratio), "target", rand[sfx], None
    for i, (basis, lv, sfx, _, ratio) in enumerate(tgt):
        # i = index in RATIOS for the <r>x levels, as in the review; the <f>max levels use separate RNG streams.
        ri = i if basis == "truth" else 100 + i
        yield (f"shuffle50_{sfx}", meta("shuffle50", basis, lv, ratio, shuffle_frac=SHUFFLE_FRAC), "target",
               rand[sfx], ("shuffle", SHUFFLE_FRAC, ri))
    truth_levels = [x for x in tgt if x[0] == "truth"]
    for k, (basis, lv, sfx, n, ratio) in enumerate(truth_levels):
        start = block_start(seed, 2, k, n_pool, n)
        yield (f"block_{sfx}", meta("block", basis, lv, ratio, block_start=start, block_draw=BLOCK_DRAW), "target",
               pool[start:start + n], None)
    for basis, lv, sfx, n, ratio in truth_levels:
        yield f"head_{sfx}", meta("head", basis, lv, ratio), "target", pool[:n], None
    for k, (r, f, unit) in enumerate(JITTERS):
        yield (f"jitter_{r}x_{f}{unit}", meta(f"jitter_{unit}", "truth", r, float(r), jitter_f=f, jitter_unit=unit),
               "target", rand[f"{r}x"], ("jitter", f, unit, k))

    ref_perm = np.random.default_rng([seed, 4]).permutation(n_ref)
    for basis, lv, sfx, n, ratio in ref:
        yield f"refstage_{sfx}", meta("refstage", basis, lv, ratio), "ref", ref_pool[np.sort(ref_perm[:n])], None
    ref_truth_levels = [x for x in ref if x[0] == "truth"]
    for k, (basis, lv, sfx, n, ratio) in enumerate(ref_truth_levels):
        start = block_start(seed, 5, k, n_ref, n)
        yield (f"refblock_{sfx}", meta("refblock", basis, lv, ratio, block_start=start, block_draw=BLOCK_DRAW),
               "ref", ref_pool[start:start + n], None)
    for basis, lv, sfx, n, ratio in ref_truth_levels:
        yield f"refhead_{sfx}", meta("refhead", basis, lv, ratio), "ref", ref_pool[:n], None


def perturb(pred: Stage, spec, seed: int, nn_own: float, nn_truth: float) -> tuple[Stage, float]:
    """Apply a jitter or expression-shuffle perturbation; returns the new stage and the jitter RMS (0 if none)."""
    if spec is None:
        return pred, 0.0
    if spec[0] == "jitter":
        _, f, unit, k = spec
        rms = f * (nn_own if unit == "own" else nn_truth)
        noise = np.random.default_rng([seed, 3, k]).normal(size=pred.coords.shape) * (rms / np.sqrt(3))
        return Stage(X=pred.X, coords=(pred.coords + noise).astype(np.float32), celltype=pred.celltype), rms
    _, frac, ri = spec
    rng = np.random.default_rng([seed, 7, ri])
    sel = rng.choice(pred.n, int(round(frac * pred.n)), replace=False)
    X = pred.X.copy()
    X[sel] = X[sel[rng.permutation(len(sel))]]
    return Stage(X=X, coords=pred.coords, celltype=None), 0.0


def resolve(path: str) -> tuple[Path | None, str]:
    """The file itself, or the single <stem>*.h5ad next to where it should be (release file names are not known
    in advance). Returns (None, reason) when there is no unique match."""
    p = Path(path)
    if p.exists():
        return p, ""
    matches = sorted(p.parent.glob(f"{p.stem}*.h5ad")) if p.parent.exists() else []
    if len(matches) == 1:
        return matches[0], f"{p} not found; using {matches[0]}"
    if matches:
        return None, (f"{p} not found and several candidates exist ({', '.join(map(str, matches))}); "
                      "pass --target-file/--ref-file")
    return None, f"{p} not found"


def finite(v):
    try:
        v = float(v)
    except (TypeError, ValueError):
        return None
    return v if np.isfinite(v) else None


def done_conditions(path: Path) -> set[str]:
    """Conditions already in a JSONL file. Refuses files holding block rows from the first, nested-centre draw, so
    the two block designs never mix under one label."""
    if not path.exists():
        return set()
    rows = [json.loads(line) for line in path.read_text().splitlines() if line.strip()]
    stale = [r["condition"] for r in rows if r.get("pattern") in ("block", "refblock")
             and r.get("block_draw") != BLOCK_DRAW]
    if stale:
        raise SystemExit(f"{path} holds {len(stale)} block rows from the old nested-centre draw (e.g. {stale[0]}); "
                         "write to a fresh --out directory such as runs/density/raw_v2")
    return {r["condition"] for r in rows}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", required=True, choices=sorted(TARGETS))
    ap.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2, 3, 4])
    ap.add_argument("--out", default="runs/density/raw_v2")
    ap.add_argument("--only", nargs="*", default=None, help="run only these condition labels (timing/debug)")
    ap.add_argument("--patterns", nargs="*", default=None,
                    help="run only these patterns (random shuffle50 block head jitter refstage refblock refhead)")
    ap.add_argument("--max-cells", type=int, default=None, help="override max_cells for the <f>max levels")
    ap.add_argument("--target-file", default=None, help="override the target .h5ad path")
    ap.add_argument("--ref-file", default=None, help="override the reference .h5ad path")
    args = ap.parse_args()

    spec = TARGETS[args.target]
    tpath, tnote = resolve(args.target_file or spec["target"])
    rpath, rnote = resolve(args.ref_file or spec["ref"])
    if tpath is None or rpath is None:
        when = " (validation answers are released on 2026-10-20)" if spec.get("release") else ""
        print(f"[{args.target}] skipped: {tnote if tpath is None else rnote}{when}", flush=True)
        return
    for note in (tnote, rnote):
        if note:
            print(f"[{args.target}] {note}", flush=True)

    board = get_board(spec["board"])
    t0 = time.time()
    target = Stage.from_anndata(load_stage(str(tpath), genes=board.genes, lineage=False), board.genes,
                                need_coords=True)
    ref = Stage.from_anndata(load_stage(str(rpath), genes=board.genes, lineage=False), board.genes,
                             need_coords=True)
    print(f"[{args.target}] board={board.key} genes={board.n_genes} target={target.n} ({tpath}) ref={ref.n} "
          f"({rpath}) loaded in {time.time() - t0:.1f}s", flush=True)
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    for seed in args.seeds:
        wc_idx, truth_idx, _, ref_wc_idx = emulator_indices(target.n, ref.n, seed)
        # The bookkeeping must reproduce score_board's own working copies and truth half exactly.
        wc = subsample(target, SUBSAMPLE_FRAC, seed)
        assert np.array_equal(target.take(wc_idx).X, wc.X) and np.array_equal(target.take(wc_idx).coords, wc.coords)
        truth = split_half(wc, seed=0)[0]
        assert np.array_equal(target.take(truth_idx).coords, truth.coords)
        assert np.array_equal(ref.take(ref_wc_idx).X, subsample(ref, SUBSAMPLE_FRAC, seed + 1).X)
        pool = np.setdiff1d(np.arange(target.n), wc_idx)
        ref_pool = np.setdiff1d(np.arange(ref.n), ref_wc_idx)
        truth_cells, ref_cells = len(truth_idx), len(ref_wc_idx)
        rule = rule_max_cells(truth_cells, ref_cells)
        if args.max_cells:
            max_cells, max_source = args.max_cells, "cli"
        elif spec.get("release") and board.max_cells:
            max_cells, max_source = int(board.max_cells), "board"
        else:
            max_cells, max_source = rule, "rule"
        print(f"[{args.target} s{seed}] truth_cells={truth_cells} ref_cells={ref_cells} max_cells={max_cells} "
              f"({max_source}; rule gives {rule}; {max_cells / truth_cells:.2f}x truth) pool={len(pool)} "
              f"ref_pool={len(ref_pool)}", flush=True)
        if spec.get("release"):
            # Release targets only: the emulated 10% working copy should give the board's published truth half.
            print(f"[{args.target} s{seed}] board spec truth_cells={board.truth_cells} max_cells={board.max_cells}",
                  flush=True)
        nn_truth = median_nn(truth.coords)

        path = out_dir / f"{args.target}_s{seed}.jsonl"
        done = done_conditions(path)
        for label, meta, source, idx, pert in conditions(pool, ref_pool, truth_cells, max_cells, seed):
            if label in done or (args.only and label not in args.only):
                continue
            if args.patterns and not {meta["pattern"], meta["pattern"].split("_")[0]} & set(args.patterns):
                continue
            assert len(np.intersect1d(idx, wc_idx if source == "target" else ref_wc_idx)) == 0
            pred = (target if source == "target" else ref).take(idx)
            nn_own = median_nn(pred.coords)
            pred, jitter_rms = perturb(pred, pert, seed, nn_own, nn_truth)
            t1 = time.time()
            res = score_board(board, pred, target, ref, seed=seed, anchors="computed")
            row = dict(
                target=args.target, board=board.key, seed=seed, condition=label, source=source, **meta,
                n_pred=pred.n, truth_cells=truth_cells, ref_cells=res.extras["ref_cells"],
                max_cells=max_cells, max_cells_source=max_source,
                pool_cells=len(pool) if source == "target" else len(ref_pool),
                n_ratio=pred.n / truth_cells, nn_own=nn_own, nn_truth=nn_truth, jitter_rms=jitter_rms,
                score=res.score, group_skill=res.group_skill,
                skills={m: finite(res.skills.get(m)) for m in RANKED_METRICS["T2"]},
                metrics={m: finite(res.metrics.get(m)) for m in RANKED_METRICS["T2"] + EXTRA_METRICS},
                floor={m: finite(res.anchors.get(m, {}).get("floor")) for m in RANKED_METRICS["T2"]},
                ceiling={m: finite(res.anchors.get(m, {}).get("ceiling")) for m in RANKED_METRICS["T2"]},
                ceiling_extra={m: finite(res.extras["ceiling_metrics"].get(m)) for m in EXTRA_METRICS},
                inverted_anchors=res.extras["inverted_anchors"], seconds=time.time() - t1,
            )
            with path.open("a") as fh:
                fh.write(json.dumps(row) + "\n")
            nb, dice = (row["skills"][m] for m in ("neighborhood_mmd", "occupancy_dice"))
            print(f"[{args.target} s{seed}] {label:24s} n={pred.n:6d} score={res.score:6.2f} nb={nb} dice={dice} "
                  f"({row['seconds']:.1f}s)", flush=True)
    print(f"[{args.target}] done in {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
