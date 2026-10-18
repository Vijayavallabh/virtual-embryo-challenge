# Density experiment: submitted cell count, sampling pattern and jitter vs T2 metrics

2026-10-18. Local emulator only: no submissions, no external data, no anchors from `boards/index.json`. Code in
`experiments/density/`, raw outputs on the DGX in `runs/density/`, small results in `experiments/density/results/`.


_Saved by the lead on 2026-10-19 from the experiment agent's returned text (the subagent harness blocks agents from writing report files). The experiment text is unedited except where marked; the independent review follows it, and the lead's decisions are at the end._

## Answers

1. **Optimal count (assumption A13).** For a prediction drawn from the target's own distribution, the board score climbs
   steeply up to 1x truth_cells and is flat from 1x to 8x. Going below 1x is expensive: 0.5x costs 1.8 to 4.7 points and
   0.25x costs 5.1 to 9.5 points, and all 15 seeds lose at both counts. From 1x to 8x every mean paired change is within
   ±0.8 points, inside seed noise. A floor-quality prediction (cells of the preceding stage) shows the same cliff below
   1x (−1.0 to −5.5 points) and gains a little from more cells (+0.4 to +1.6 points at 8x versus 1x). Raw
   `neighborhood_mmd` is U-shaped, with its minimum at 2x on all three targets. **A13 passes, but with a modified
   window: the near-optimal range is 1x to 4x, not 0.5x to 2x; about 2x is the best single choice, and 0.5x is clearly
   suboptimal.**
   Confidence: High that counts below 1x lose points. High that 1x to 4x is a plateau. Medium that 2x is the best
   single count. Medium-low for the 8x end on real boards, where the sign depends on prediction quality.
2. **Contiguous-rows penalty.** Rows in every released T2 file are sorted by z section (Spearman(row, z) = −0.999 to
   −1.000). At 1x, a contiguous block costs a near-perfect prediction 16 to 28 points when it comes from mid-file and
   21 to 48 points when it is the first rows of the file. Even at 8x, where the block holds 44% of the rows, it still
   costs 12 to 21 points. For a floor-quality prediction the loss at 1x is 3 to 7 points (mid-file) and 3.5 to 18
   points (first rows). Most of the loss is in the shape group. Confidence: High.
3. **Jitter.** Isotropic Gaussian jitter with RMS equal to f times the median nearest-neighbour (NN) distance,
   submitted at 1x density:

   | f | mean change (points) | seeds that lose |
   |---|---|---|
   | 0.25 | +0.06 (nothing) | 6/15 |
   | 0.5 | −0.14 | 12/15 |
   | 1 | −0.29 | 13/15 |
   | 2 | −1.46 | 15/15 |

   One NN step is about 20 to 22 coordinate units, or 6 to 9% of the RMS radius. A denser submission tolerates the same
   physical displacement better. At 4x, jitter of up to 1 truth-half NN distance costs nothing (−0.01 to +0.16
   points), and 2 NN costs −0.96. Most of the loss comes through `occupancy_dice`, then `neighborhood_mmd`.
   Confidence: High for the shape of the curve. Medium for the exact sizes below 1 NN, which are at the level of
   paired noise.

## Setup

| Target (stands in for the held-out stage) | DE reference and copy_last floor | Board panel | truth_cells | Pool (target cells outside the working copy) |
|---|---|---|---|---|
| heart E9.5 (53,742 cells) | E8.75 | `T2:heart:val_extrap` (500 genes) | 2,687 | 48,368 |
| heart E8.75 (24,826) | E8.25_late | `T2:heart:val_extrap` (500) | 1,241 | 22,343 |
| embryo E8.0 (31,671; subset to the 498-gene panel) | E7.25 | `T2:embryo:val_interp` (498) | 1,583 | 28,504 |

- **Scoring.** Each call is `vec.evaluate.score_board(board, pred, target, ref, seed=s, anchors="computed")` with the
  default settings: 10% working copies, truth and ceiling halves from `split_half(seed=0)`, floor-normalised structure
  gate, no probe. Seeds are s = 0 to 4. Stages are loaded with `vec.prep.load_stage(path, genes=board.genes)` and
  converted with `Stage.from_anndata`.
- **Bookkeeping.** `run_density.py` runs `vec.evaluate.subsample` and `split_half` on a column of row indices with
  score_board's seeds, which recovers the rows of the target working copy (seed s), the truth half, and the reference
  working copy (seed s+1). Assertions check that these rows reproduce score_board's own working copies exactly. All
  predictions are drawn from rows outside the working copy, so they never share a cell with the truth or ceiling
  halves (also asserted per call). truth_cells is computed from the target (round(0.1 × N) // 2), not taken from
  `index.json`.
- **Conditions per (target, seed)**, with ratio r ∈ {0.25, 0.5, 1, 2, 4, 8} × truth_cells. No cap was needed: the
  largest draw, 8x, is 44% of the pool.
  - (a) `random`: random subsets of the pool. The subsets are nested: one permutation per seed, with prefixes of it
    as the counts.
  - (b) `block`: contiguous runs of pool rows in file order, nested around one random centre per seed (chosen so the
    8x block fits). `head`: the first n pool rows, which is what `a[:n]` submits.
  - (c) `jitter`: the same random subsets at 1x and 4x, with isotropic Gaussian noise of RMS length f × d,
    f ∈ {0.25, 0.5, 1, 2}. `own` means d is the submitted cloud's own median NN distance. `truth` means d is the truth
    half's, a fixed physical length, so 4x/`truth` moves cells as far as 1x/`own`. The task asked for 0.25 and 0.5; 1
    and 2 were added to find where jitter starts to hurt.
  - Supplementary: `refstage`, random subsets of the reference stage's pool (cells outside the floor working copy)
    at the same counts, with coordinates untouched. This is the count sweep for a floor-quality prediction. `refblock`
    and `refhead` are contiguous blocks of that pool. For embryo E8.0 the reference pool (11,965 cells) is smaller than
    8x, so the three reference conditions at 8x are the same set of cells.
- **Size of the run.** 720 score_board calls (3 targets × 5 seeds × 48 conditions), each 4 to 8 s for three metric
  panels. The run used 15 CPU processes with 2 BLAS threads each, about 10 minutes of wall time over two passes. No GPU
  was used.

## Commands

```bash
# local -> DGX code sync (before each DGX step)
cd /mnt/md0/IITM/BackUp/Home/vijayavallabh && rsync -az --exclude .venv --exclude __pycache__ --exclude data \
  --exclude runs --exclude .ruff_cache --exclude .pytest_cache -e ssh virtual-embryo-challenge/ \
  PrakashDGX_H2:v/virtual-embryo-challenge/

# DGX, in ~/v/virtual-embryo-challenge, with export PATH=$HOME/.local/bin:$PATH
# smoke tests (timing and one condition of each pattern; outputs deleted afterwards)
OMP_NUM_THREADS=2 uv run python experiments/density/run_density.py --target heart_E8.75 --seeds 0 \
  --out runs/density/timing --only random_1x random_8x
OMP_NUM_THREADS=2 uv run python experiments/density/run_density.py --target embryo_E8.0 --seeds 1 \
  --out runs/density/timing --only random_0.25x block_2x head_1x jitter_4x_0.5truth jitter_1x_2.0own refstage_8x
# full run: one process per (target, seed); resumable (finished conditions are skipped). It ran twice: the first pass
# without refblock/refhead, the second adding them. Per-call results are appended to runs/density/raw/<target>_s<seed>.jsonl
tmux new -d -s vec-density 'bash experiments/density/run_all.sh'   # logs: runs/logs/density_<target>_s<seed>.log
# structure-only diagnostics (row order, block coverage, NN spacing)
OMP_NUM_THREADS=2 uv run python experiments/density/diag_order.py --out runs/density/diag_order.csv
# summary CSVs, printed tables, plots
OMP_NUM_THREADS=2 uv run python experiments/density/summarise.py --raw runs/density/raw --out runs/density/summary \
  > runs/density/summary_tables.txt

# local: bring back the small results, then build the markdown tables used below
scp 'PrakashDGX_H2:v/virtual-embryo-challenge/runs/density/summary/*' experiments/density/results/
scp PrakashDGX_H2:v/virtual-embryo-challenge/runs/density/{diag_order.csv,summary_tables.txt} experiments/density/results/
python3 -I experiments/density/report_tables.py experiments/density/results
```

## Results

All scores are on the emulator's 0 to 100 board scale: mean (sd) over 5 seeds. "Δ" is a paired difference within a
seed, with the same working copy and truth half. `[k/5]` counts the seeds where Δ < 0. Skills are clipped at 1, so a
prediction at or above the split-half ceiling scores 1 on that metric.

### 0. Row order and coverage (`results/diag_order.csv`, structure only)

| Stage | Cells | z levels | Spearman(row, z) | RMS radius | Median NN, full stage | Median NN, truth half |
|---|---|---|---|---|---|---|
| E9.5 | 53,742 | 98 | −0.999 | 335 | 9.0 | 21.5 |
| E8.75 | 24,826 | 35 | −1.000 | 217 | 8.7 | 20.1 |
| E8.25_late | 58,716 | 38 | −1.000 | 354 | 8.7 | - |
| E8.0 | 31,671 | 30 | −0.999 | 327 | 8.4 | 19.9 |
| E7.25 | 13,295 | 40 | −1.000 | 147 | 7.2 | - |

- Rows are written section by section, in descending z.
- At 1x, a mid-file block covers 6-10% of the z levels and 21-25% of the target's occupied voxels (16³ canonical
  grid, as used by `occupancy_dice`). A random subset of the same size covers 94-100% of z levels and 81-83% of voxels.
- At 8x, a block covers 43-48% of z levels and 50-56% of voxels.
- A head block's RMS radius is only 0.64-0.94 of the full cloud's at 1x.
- The median NN distance of random subsets falls from 21.4 to 11.0 between 1x and 8x on E9.5. That factor of 1.95 is
  close to 8^(1/3) = 2, so for NN spacing the clouds behave as 3-D point sets.

### 1. Count sweep, near-perfect predictions (random subsets of the target)

| target | cells/truth | n | score | Δ vs 1x [seeds<0] | de | dir | mmd_u | vario | d2 | dice | scale | nbhd |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| heart_E9.5 | 0.25x | 672 | 94.32 (0.96) | -5.10 (1.11) [5/5] | 0.951 (0.027) | 0.973 (0.006) | 1.000 (0.000) | 0.996 (0.004) | 0.952 (0.038) | 0.662 (0.074) | 0.992 (0.008) | 0.944 (0.009) |
| heart_E9.5 | 0.5x | 1344 | 97.65 (1.14) | -1.77 (1.15) [5/5] | 0.973 (0.013) | 0.993 (0.003) | 0.999 (0.000) | 0.999 (0.002) | 0.982 (0.028) | 0.849 (0.096) | 0.995 (0.007) | 0.982 (0.001) |
| heart_E9.5 | 1x | 2687 | 99.42 (0.31) | +0.00 (0.00) [0/5] | 0.985 (0.033) | 1.000 (0.000) | 0.999 (0.002) | 1.000 (0.000) | 0.981 (0.021) | 0.980 (0.030) | 0.995 (0.008) | 1.000 (0.001) |
| heart_E9.5 | 2x | 5374 | 99.41 (0.17) | -0.01 (0.38) [3/5] | 0.997 (0.007) | 1.000 (0.000) | 0.999 (0.001) | 0.999 (0.003) | 0.979 (0.025) | 0.962 (0.028) | 0.996 (0.008) | 1.000 (0.000) |
| heart_E9.5 | 4x | 10748 | 99.44 (0.44) | +0.02 (0.43) [1/5] | 0.997 (0.007) | 1.000 (0.000) | 0.999 (0.002) | 0.999 (0.001) | 0.971 (0.035) | 0.973 (0.041) | 0.999 (0.002) | 0.999 (0.002) |
| heart_E9.5 | 8x | 21496 | 99.19 (0.64) | -0.23 (0.76) [3/5] | 0.997 (0.007) | 1.000 (0.000) | 0.998 (0.002) | 0.999 (0.002) | 0.979 (0.028) | 0.968 (0.044) | 1.000 (0.000) | 0.988 (0.002) |
| heart_E8.75 | 0.25x | 310 | 89.18 (0.79) | -9.52 (1.20) [5/5] | 0.923 (0.020) | 0.958 (0.016) | 0.996 (0.006) | 0.973 (0.005) | 0.981 (0.023) | 0.377 (0.030) | 0.980 (0.022) | 0.861 (0.016) |
| heart_E8.75 | 0.5x | 620 | 94.05 (0.98) | -4.66 (0.73) [5/5] | 0.962 (0.029) | 0.981 (0.011) | 0.997 (0.005) | 0.989 (0.004) | 0.983 (0.022) | 0.572 (0.073) | 0.994 (0.009) | 0.947 (0.009) |
| heart_E8.75 | 1x | 1241 | 98.71 (1.31) | +0.00 (0.00) [0/5] | 0.978 (0.020) | 0.998 (0.004) | 0.997 (0.002) | 0.998 (0.003) | 0.996 (0.008) | 0.919 (0.148) | 0.994 (0.005) | 0.993 (0.006) |
| heart_E8.75 | 2x | 2482 | 99.51 (0.33) | +0.80 (1.14) [1/5] | 0.993 (0.010) | 1.000 (0.000) | 0.996 (0.002) | 1.000 (0.001) | 0.985 (0.026) | 0.976 (0.030) | 0.997 (0.003) | 1.000 (0.000) |
| heart_E8.75 | 4x | 4964 | 98.69 (0.95) | -0.02 (0.97) [1/5] | 0.993 (0.010) | 1.000 (0.000) | 0.998 (0.002) | 0.998 (0.002) | 0.984 (0.023) | 0.884 (0.114) | 0.998 (0.002) | 0.997 (0.005) |
| heart_E8.75 | 8x | 9928 | 98.67 (1.47) | -0.04 (0.29) [2/5] | 0.993 (0.010) | 1.000 (0.000) | 0.996 (0.005) | 0.999 (0.001) | 0.991 (0.021) | 0.911 (0.153) | 0.999 (0.002) | 0.986 (0.010) |
| embryo_E8.0 | 0.25x | 396 | 93.77 (0.54) | -5.42 (0.45) [5/5] | 0.969 (0.021) | 0.981 (0.005) | 0.999 (0.001) | 0.994 (0.002) | 0.972 (0.053) | 0.460 (0.096) | 0.994 (0.010) | 0.970 (0.008) |
| embryo_E8.0 | 0.5x | 792 | 96.99 (0.54) | -2.20 (0.73) [5/5] | 0.983 (0.016) | 0.994 (0.004) | 0.999 (0.001) | 0.998 (0.001) | 0.986 (0.020) | 0.721 (0.094) | 0.997 (0.005) | 0.991 (0.005) |
| embryo_E8.0 | 1x | 1583 | 99.18 (0.54) | +0.00 (0.00) [0/5] | 0.993 (0.015) | 1.000 (0.000) | 1.000 (0.000) | 1.000 (0.000) | 0.987 (0.024) | 0.929 (0.072) | 0.997 (0.007) | 1.000 (0.000) |
| embryo_E8.0 | 2x | 3166 | 99.05 (1.27) | -0.13 (1.14) [2/5] | 0.998 (0.006) | 1.000 (0.000) | 1.000 (0.000) | 1.000 (0.000) | 0.998 (0.005) | 0.895 (0.145) | 0.998 (0.005) | 1.000 (0.000) |
| embryo_E8.0 | 4x | 6332 | 99.08 (0.97) | -0.10 (0.79) [2/5] | 1.000 (0.000) | 1.000 (0.000) | 1.000 (0.000) | 1.000 (0.000) | 0.982 (0.029) | 0.911 (0.118) | 0.998 (0.005) | 1.000 (0.000) |
| embryo_E8.0 | 8x | 12664 | 99.21 (0.70) | +0.03 (0.72) [3/5] | 1.000 (0.000) | 1.000 (0.000) | 1.000 (0.000) | 1.000 (0.000) | 0.983 (0.017) | 0.937 (0.079) | 0.999 (0.003) | 0.996 (0.001) |

Pooled over the three targets, the change against 1x is −6.68 at 0.25x, −2.87 at 0.5x, +0.22 at 2x, −0.03 at 4x and
−0.08 at 8x. Below 1x, most of the loss comes from `occupancy_dice`. veckit subsamples both clouds to min(n_pred,
n_true), so a short submission also thins the truth. At 1x the voxel edge is only 4 to 6 truth NN spacings, below the
code's own threshold of 10 for trusting Dice, so Dice behaves as a sampling-density statistic. Over-smoothed
neighbourhoods (`neighborhood_mmd`) and pseudobulk noise (`de_score`, `de_direction`) add to the loss.

Raw `neighborhood_mmd` divided by the split-half ceiling value (mean of 5 seeds):

| target | nbhd floor/ceiling | 0.25x | 0.5x | 1x | 2x | 4x | 8x |
|---|---|---|---|---|---|---|---|
| heart_E9.5 | 73 | 5.46 | 2.38 | 0.94 | 0.66 | 0.94 | 1.94 |
| heart_E8.75 | 31 | 6.09 | 2.76 | 1.21 | 0.66 | 0.90 | 1.47 |
| embryo_E8.0 | 139 | 5.64 | 2.45 | 0.86 | 0.57 | 0.92 | 1.63 |

The U-shape and its minimum at 2x repeat on all three targets. On these pseudo-boards the floor is 31 to 139 times the
ceiling, so the rise at 8x is only about 1% of the floor-to-ceiling gap (skill 0.986 to 0.996). See caveat 1 for
boards with a narrower gap.

A possible mechanism, not tested directly (Medium-low confidence). Nearby cells share most of their 15 neighbours, so
their neighbourhood profiles are near-duplicates, and those pairs inflate the within-sample kernel term of the
unbiased MMD. veckit draws at most 2,000 profiles, so the share of such pairs among predicted profiles falls roughly as
1/n_pred: a denser submission carries less of this positive bias. By 8x, neighbourhoods have become physically smaller
than the truth's (a smoothing mismatch), and that effect wins. The same bias argument would explain why floor-quality
predictions keep improving up to 8x (section 2).

### 2. Count sweep, floor-quality predictions (random subsets of the reference stage)

| target | cells/truth | score | Δ vs refstage 1x | nbhd skill | dice skill | de skill |
|---|---|---|---|---|---|---|
| heart_E9.5 | 0.25x | 47.94 (1.23) | -2.56 | 0.469 (0.010) | 0.448 (0.037) | 0.460 (0.012) |
| heart_E9.5 | 0.5x | 49.42 (1.09) | -1.08 | 0.494 (0.009) | 0.493 (0.012) | 0.457 (0.021) |
| heart_E9.5 | 1x | 50.50 (1.15) | +0.00 | 0.513 (0.009) | 0.495 (0.026) | 0.470 (0.022) |
| heart_E9.5 | 2x | 50.89 (1.17) | +0.39 | 0.525 (0.010) | 0.507 (0.024) | 0.468 (0.022) |
| heart_E9.5 | 4x | 50.54 (1.35) | +0.04 | 0.528 (0.011) | 0.506 (0.038) | 0.467 (0.032) |
| heart_E9.5 | 8x | 50.89 (1.00) | +0.39 | 0.534 (0.007) | 0.513 (0.019) | 0.467 (0.034) |
| heart_E8.75 | 0.25x | 42.82 (1.56) | -5.47 | 0.400 (0.014) | 0.304 (0.031) | 0.435 (0.014) |
| heart_E8.75 | 0.5x | 46.16 (0.99) | -2.14 | 0.442 (0.007) | 0.409 (0.019) | 0.451 (0.014) |
| heart_E8.75 | 1x | 48.30 (1.00) | +0.00 | 0.471 (0.008) | 0.498 (0.043) | 0.449 (0.014) |
| heart_E8.75 | 2x | 48.39 (0.65) | +0.10 | 0.487 (0.002) | 0.478 (0.036) | 0.436 (0.012) |
| heart_E8.75 | 4x | 49.20 (0.19) | +0.90 | 0.498 (0.007) | 0.491 (0.074) | 0.440 (0.015) |
| heart_E8.75 | 8x | 49.94 (0.50) | +1.64 | 0.508 (0.004) | 0.499 (0.045) | 0.449 (0.020) |
| embryo_E8.0 | 0.25x | 47.31 (0.37) | -2.16 | 0.485 (0.004) | 0.383 (0.032) | 0.405 (0.020) |
| embryo_E8.0 | 0.5x | 48.47 (1.02) | -1.00 | 0.493 (0.005) | 0.471 (0.026) | 0.408 (0.025) |
| embryo_E8.0 | 1x | 49.47 (0.79) | +0.00 | 0.500 (0.006) | 0.536 (0.014) | 0.407 (0.027) |
| embryo_E8.0 | 2x | 49.78 (0.74) | +0.31 | 0.507 (0.005) | 0.524 (0.027) | 0.416 (0.025) |
| embryo_E8.0 | 4x | 50.14 (0.92) | +0.67 | 0.512 (0.006) | 0.536 (0.015) | 0.418 (0.029) |
| embryo_E8.0 | 8x | 50.20 (0.76) | +0.73 | 0.516 (0.004) | 0.541 (0.023) | 0.418 (0.025) |

All 15 seeds lose at 0.25x and at 0.5x. Above 1x the changes are positive on average, and fewer seeds lose: 5/15 at
2x, 5/15 at 4x and 3/15 at 8x. The gain runs mostly through `neighborhood_mmd` skill, which rises by 0.016 to 0.037
from 1x to 8x. Unlike the near-perfect case, nothing turns down by 8x.

Side observation. A "no change" prediction built from fresh reference cells scores below the floor on `de_score`
(skill 0.41 to 0.47, against exactly 0.5 for the floor). The floor avoids this only because it is the server's own
reference working copy, which gives dp = 0. Any other cells bring that copy's sampling noise into dp, and the noise
then gets ranked.

### 3. Contiguous blocks

Board score per condition, with the paired Δ against random sampling of the same source at the same count in brackets
(mean ± sd). The block and head columns compare with random; the refblock and refhead columns compare with refstage.

| target | cells/truth | random | block (Δ same n) | head (Δ same n) | refstage | refblock (Δ) | refhead (Δ) |
|---|---|---|---|---|---|---|---|
| heart_E9.5 | 0.25x | 94.3 | 67.0 (-27.3 ± 5.9) | 45.4 (-49.0 ± 2.0) | 47.9 | 39.5 (-8.4 ± 1.4) | 31.8 (-16.1 ± 0.5) |
| heart_E9.5 | 0.5x | 97.7 | 72.3 (-25.4 ± 3.6) | 50.4 (-47.3 ± 1.7) | 49.4 | 42.7 (-6.7 ± 2.2) | 31.6 (-17.8 ± 0.7) |
| heart_E9.5 | 1x | 99.4 | 73.9 (-25.5 ± 2.3) | 53.0 (-46.5 ± 0.8) | 50.5 | 44.6 (-5.9 ± 1.2) | 32.4 (-18.1 ± 0.8) |
| heart_E9.5 | 2x | 99.4 | 75.8 (-23.6 ± 2.6) | 55.8 (-43.6 ± 0.8) | 50.9 | 45.4 (-5.5 ± 0.7) | 36.1 (-14.8 ± 0.8) |
| heart_E9.5 | 4x | 99.4 | 76.9 (-22.5 ± 2.3) | 62.4 (-37.0 ± 0.7) | 50.5 | 52.3 (+1.8 ± 1.2) | 45.0 (-5.5 ± 0.6) |
| heart_E9.5 | 8x | 99.2 | 81.3 (-17.9 ± 2.3) | 78.9 (-20.3 ± 1.0) | 50.9 | 50.4 (-0.5 ± 1.2) | 49.2 (-1.7 ± 0.4) |
| heart_E8.75 | 0.25x | 89.2 | 61.2 (-28.0 ± 5.6) | 43.7 (-45.5 ± 1.1) | 42.8 | 33.9 (-8.9 ± 5.3) | 33.6 (-9.2 ± 1.7) |
| heart_E8.75 | 0.5x | 94.0 | 67.6 (-26.5 ± 8.2) | 50.6 (-43.4 ± 1.0) | 46.2 | 38.4 (-7.7 ± 9.7) | 41.4 (-4.8 ± 1.4) |
| heart_E8.75 | 1x | 98.7 | 70.9 (-27.8 ± 4.0) | 50.6 (-48.1 ± 1.4) | 48.3 | 41.0 (-7.3 ± 5.2) | 40.2 (-8.1 ± 1.1) |
| heart_E8.75 | 2x | 99.5 | 73.3 (-26.2 ± 0.8) | 53.9 (-45.6 ± 0.9) | 48.4 | 43.6 (-4.8 ± 5.4) | 40.4 (-7.9 ± 1.0) |
| heart_E8.75 | 4x | 98.7 | 74.8 (-23.9 ± 1.6) | 60.8 (-37.9 ± 1.3) | 49.2 | 44.5 (-4.7 ± 5.0) | 39.1 (-10.0 ± 0.6) |
| heart_E8.75 | 8x | 98.7 | 80.5 (-18.2 ± 2.0) | 77.9 (-20.7 ± 1.1) | 49.9 | 44.6 (-5.3 ± 6.0) | 39.6 (-10.4 ± 0.5) |
| embryo_E8.0 | 0.25x | 93.8 | 75.9 (-17.9 ± 4.4) | 68.6 (-25.2 ± 1.0) | 47.3 | 42.6 (-4.7 ± 0.6) | 43.9 (-3.4 ± 0.1) |
| embryo_E8.0 | 0.5x | 97.0 | 82.2 (-14.8 ± 2.4) | 74.9 (-22.1 ± 1.0) | 48.5 | 44.3 (-4.1 ± 0.9) | 45.3 (-3.2 ± 0.8) |
| embryo_E8.0 | 1x | 99.2 | 83.4 (-15.8 ± 1.4) | 77.9 (-21.3 ± 0.8) | 49.5 | 46.2 (-3.3 ± 0.6) | 45.9 (-3.5 ± 0.5) |
| embryo_E8.0 | 2x | 99.1 | 84.5 (-14.5 ± 1.5) | 78.6 (-20.4 ± 1.9) | 49.8 | 48.3 (-1.4 ± 0.5) | 46.5 (-3.2 ± 0.3) |
| embryo_E8.0 | 4x | 99.1 | 85.5 (-13.6 ± 1.1) | 81.3 (-17.8 ± 1.5) | 50.1 | 49.8 (-0.4 ± 0.5) | 51.4 (+1.3 ± 0.4) |
| embryo_E8.0 | 8x | 99.2 | 86.8 (-12.4 ± 1.1) | 85.0 (-14.2 ± 0.7) | 50.2 | 50.2 (+0.0 ± 0.0) | 50.2 (+0.0 ± 0.0) |

Pooled at 1x: block −23.1, head −38.6, refblock −5.5, refhead −9.9. Pooled at 8x: block −16.2, head −18.4.

Where the loss comes from (mean skills):

| target | condition | score | de | dir | mmd_u | vario | d2 | dice | scale | nbhd | shape group |
|---|---|---|---|---|---|---|---|---|---|---|---|
| heart_E9.5 | random_1x | 99.4 | 0.99 | 1.00 | 1.00 | 1.00 | 0.98 | 0.98 | 1.00 | 1.00 | 0.99 |
| heart_E9.5 | block_1x | 73.9 | 0.85 | 0.92 | 0.82 | 0.88 | 0.26 | 0.27 | 0.78 | 0.79 | 0.44 |
| heart_E9.5 | head_1x | 53.0 | 0.64 | 0.68 | 0.54 | 0.64 | 0.32 | 0.28 | 0.49 | 0.51 | 0.36 |
| heart_E9.5 | block_8x | 81.3 | 0.89 | 0.95 | 0.90 | 0.94 | 0.38 | 0.42 | 0.84 | 0.87 | 0.55 |
| heart_E9.5 | head_8x | 78.9 | 0.81 | 0.89 | 0.88 | 0.92 | 0.53 | 0.44 | 0.71 | 0.85 | 0.56 |
| heart_E8.75 | random_1x | 98.7 | 0.98 | 1.00 | 1.00 | 1.00 | 1.00 | 0.92 | 0.99 | 0.99 | 0.97 |
| heart_E8.75 | block_1x | 70.9 | 0.82 | 0.92 | 0.82 | 0.87 | 0.24 | 0.23 | 0.69 | 0.75 | 0.38 |
| heart_E8.75 | head_1x | 50.6 | 0.68 | 0.70 | 0.50 | 0.61 | 0.25 | 0.23 | 0.67 | 0.45 | 0.34 |
| heart_E8.75 | block_8x | 80.5 | 0.88 | 0.96 | 0.91 | 0.94 | 0.40 | 0.47 | 0.75 | 0.84 | 0.54 |
| heart_E8.75 | head_8x | 77.9 | 0.84 | 0.88 | 0.83 | 0.89 | 0.63 | 0.45 | 0.71 | 0.81 | 0.59 |
| embryo_E8.0 | random_1x | 99.2 | 0.99 | 1.00 | 1.00 | 1.00 | 0.99 | 0.93 | 1.00 | 1.00 | 0.97 |
| embryo_E8.0 | block_1x | 83.4 | 0.89 | 0.97 | 0.98 | 0.98 | 0.24 | 0.26 | 0.90 | 0.96 | 0.46 |
| embryo_E8.0 | head_1x | 77.9 | 0.76 | 0.91 | 0.93 | 0.95 | 0.18 | 0.25 | 0.94 | 0.89 | 0.46 |
| embryo_E8.0 | block_8x | 86.8 | 0.96 | 0.99 | 0.99 | 0.99 | 0.35 | 0.32 | 0.93 | 0.97 | 0.53 |
| embryo_E8.0 | head_8x | 85.0 | 0.90 | 0.97 | 0.98 | 0.99 | 0.26 | 0.28 | 0.98 | 0.97 | 0.51 |

- **Shape.** A slab of sections has the wrong form. `d2_shape` and `occupancy_dice` skills fall to 0.18-0.32, and
  the shape group drops to 0.34-0.46 from 0.97-0.99, which alone is about 13 to 16 board points.
- **Composition and pseudobulk.** On the heart targets the slab's mix of cell types is regional, which costs
  `mmd_u`, `variogram` and the DE terms.
- **Structure gate.** It bites only where nbhd skill drops below 0.5 (for example heart E8.75 head 1x at 0.45, and
  heart E9.5 head 0.25x and 0.5x).
- **Lucky crops.** Occasionally a slab beats random sampling because its shape or size happens to match the target
  better: E9.5 refblock at 4x (+1.8, d2 skill 0.71) and E8.0 refhead at 4x (+1.3). On the live board this kind of
  coincidence could make a block look like an improvement.
- **Consistency.** The size of the penalty fits the community report of 48.94 → 38.6 for a contiguous versus random
  E9.5 sample on the live heart-extrap board (`docs/research/community_digest.md`).

### 4. Jitter

Paired Δ against the very same cells unjittered: mean (sd) [seeds < 0]. The column "jitter RMS at f=1" gives the NN
unit in coordinate units.

| target | density / unit | jitter RMS at f=1 | f=0.25 | f=0.5 | f=1.0 | f=2.0 |
|---|---|---|---|---|---|---|
| heart_E9.5 | 1x / own | 21.7 | -0.02 (0.06) [3/5] | -0.04 (0.27) [4/5] | -0.21 (0.45) [4/5] | -1.46 (0.57) [5/5] |
| heart_E9.5 | 4x / own | 14.1 | +0.02 (0.12) [2/5] | +0.12 (0.17) [1/5] | -0.36 (0.51) [4/5] | -0.32 (0.57) [4/5] |
| heart_E9.5 | 4x / truth | 21.6 | -0.09 (0.16) [3/5] | -0.00 (0.44) [3/5] | -0.05 (0.31) [2/5] | -1.11 (0.63) [5/5] |
| heart_E8.75 | 1x / own | 20.1 | +0.13 (0.24) [2/5] | -0.28 (0.30) [4/5] | -0.18 (0.10) [5/5] | -1.75 (0.57) [5/5] |
| heart_E8.75 | 4x / own | 14.0 | -0.17 (0.40) [2/5] | +0.25 (0.28) [0/5] | +0.02 (0.30) [3/5] | -0.41 (0.26) [5/5] |
| heart_E8.75 | 4x / truth | 20.1 | +0.09 (0.15) [1/5] | +0.40 (0.71) [1/5] | +0.55 (0.69) [1/5] | -0.42 (0.47) [4/5] |
| embryo_E8.0 | 1x / own | 20.4 | +0.06 (0.07) [1/5] | -0.11 (0.49) [4/5] | -0.48 (0.41) [4/5] | -1.17 (0.87) [5/5] |
| embryo_E8.0 | 4x / own | 13.8 | -0.12 (0.13) [4/5] | -0.07 (0.38) [3/5] | +0.03 (0.23) [2/5] | -0.42 (0.62) [3/5] |
| embryo_E8.0 | 4x / truth | 20.5 | -0.03 (0.18) [2/5] | +0.09 (0.32) [2/5] | -0.38 (0.53) [3/5] | -1.36 (1.06) [5/5] |

Pooled changes and losing seeds:

| Density / unit | f=0.25 | f=0.5 | f=1 | f=2 |
|---|---|---|---|---|
| 1x / own | +0.06 (6/15) | −0.14 (12/15) | −0.29 (13/15) | −1.46 (15/15) |
| 4x / own | −0.09 (8/15) | +0.10 (4/15) | −0.10 (9/15) | −0.38 (12/15) |
| 4x / truth | −0.01 (6/15) | +0.16 (6/15) | +0.04 (6/15) | −0.96 (14/15) |

| target | condition | Δ score | Δ dice skill | Δ nbhd skill | Δ d2 skill | Δ scale skill | nbhd raw/ceiling |
|---|---|---|---|---|---|---|---|
| heart_E9.5 | jitter_1x_1.0own | -0.21 | -0.024 | -0.002 | +0.004 | +0.000 | 1.14 |
| heart_E9.5 | jitter_1x_2.0own | -1.46 | -0.133 | -0.013 | +0.000 | -0.006 | 1.95 |
| heart_E9.5 | jitter_4x_2.0truth | -1.11 | -0.134 | +0.001 | +0.001 | -0.003 | 0.68 |
| heart_E8.75 | jitter_1x_1.0own | -0.18 | +0.009 | -0.009 | -0.003 | -0.001 | 1.47 |
| heart_E8.75 | jitter_1x_2.0own | -1.75 | -0.062 | -0.038 | -0.022 | -0.013 | 2.45 |
| heart_E8.75 | jitter_4x_2.0truth | -0.42 | -0.007 | -0.002 | -0.014 | -0.024 | 1.03 |
| embryo_E8.0 | jitter_1x_1.0own | -0.48 | -0.054 | -0.001 | -0.003 | +0.000 | 1.02 |
| embryo_E8.0 | jitter_1x_2.0own | -1.17 | -0.100 | -0.004 | -0.027 | -0.000 | 1.69 |
| embryo_E8.0 | jitter_4x_2.0truth | -1.36 | -0.125 | +0.000 | -0.039 | +0.001 | 0.64 |

- At 1x, jitter of 2 NN roughly doubles raw `neighborhood_mmd` (1.7 to 2.4 times the ceiling, up from 0.9 to 1.2).
- At 4x, jitter lowers raw `neighborhood_mmd` instead (0.94 → 0.50 times the ceiling at 2 own-NN on E9.5). Blurring
  positions smooths the dense cloud's too-local neighbourhoods toward the truth's 5%-sample scale, but skill
  clipping hides the gain.
- `occupancy_dice` takes the steady loss in every case.

### 5. Computed anchors (mean over seeds; for translating these effects to other boards)

| target | ref (floor) | truth_cells | de floor / ceiling | dir floor / ceiling | mmd_u floor / ceiling | vario floor / ceiling | d2 floor / ceiling | dice floor / ceiling | scale floor / ceiling | nbhd floor / ceiling |
|---|---|---|---|---|---|---|---|---|---|---|
| heart_E9.5 | E8.75 | 2687 | 0 / 0.9296 | 0 / 0.9695 | 0.05872 / 0.000122 | 0.05831 / 0.000784 | 0.05053 / 0.005132 | 0.7806 / 0.9297 | -0.4306 / 0.0061 | 0.1379 / 0.001882 |
| heart_E8.75 | E8.25_late | 1241 | 0 / 0.9219 | 0 / 0.965 | 0.05503 / 0.000558 | 0.0558 / 0.0008368 | 0.04262 / 0.005172 | 0.7647 / 0.8821 | 0.4906 / 0.00428 | 0.1123 / 0.003616 |
| embryo_E8.0 | E7.25 | 1583 | 0 / 0.9214 | 0 / 0.9848 | 0.2444 / 0.000576 | 0.1942 / 0.000921 | 0.03233 / 0.00554 | 0.7579 / 0.9032 | -0.798 / -0.00488 | 0.4615 / 0.003322 |

None of the 720 calls had inverted anchors. On these pseudo-boards the emulator score varies by 0.3 to 1.3 points
across seeds for a near-perfect prediction at 1x, and Dice is the noisiest term.

Complete per-condition mean (sd) tables for the score and all eight skills, including every jitter and reference
condition, are in `results/summary_tables.txt`. Machine-readable versions are in `results/summary_by_condition.csv`
(also raw metric values, gap fractions and paired deltas) and `results/results_long.csv` (one row per call).

## Conclusions

| # | Conclusion | Confidence |
|---|---|---|
| C1 | Never submit fewer cells than truth_cells on a T2 board. 0.5x costs 1.0-4.7 points and 0.25x costs 2.2-9.5 points, for near-perfect and floor-quality predictions alike (all 15 seeds lose at each count, for both prediction types). Most of the loss is `occupancy_dice` thinning the truth to the submitted count, plus over-smoothed neighbourhoods. | High |
| C2 | Between 1x and 4x the count is a plateau for near-perfect predictions (all paired changes within ±0.8 points, within noise), and floor-quality predictions improve slowly with count. About 2x truth_cells is the best single choice: never worse than 1x, the raw `neighborhood_mmd` minimum on 3/3 targets, and a margin against misjudging truth_cells. | Medium |
| C3 | At 8x, near-perfect predictions lose about 1% of the `neighborhood_mmd` gap (raw value 1.5-1.9x the ceiling), while floor-quality ones gain (+0.4 to +1.6). The net effect on a real board depends on prediction quality and on how close that board's floor sits to its ceiling. | Medium-low |
| C4 | A13 passes with a modification: matching the truth half is near-optimal, but the window is 1x to 4x (best about 2x), not 0.5x to 2x. | High (the cliff below 1x), Medium (the window's upper end) |
| C5 | Rows are stored in z-section order. A contiguous block costs a near-perfect prediction 16-28 points (mid-file) to 21-48 (first rows) at 1x, and still 12-21 at 8x. A floor-quality prediction loses 3-7 (mid-file) to 3.5-18 (first rows) at 1x. Random sampling over all rows is mandatory, and a block that happens to score higher (lucky crops at 4x) is a coincidence, not an improvement. | High |
| C6 | Jitter up to 0.25 NN is free. 0.5-1 NN costs about 0.1-0.5 points at 1x density, 2 NN costs 1.2-1.8. At 4x density the same physical displacement costs less (nothing up to 1 truth-NN, −1.0 at 2). One truth-half NN is about 20-22 units, or 6-9% of the RMS radius at these stages. | High (shape), Medium (sizes below 1 NN) |

## Caveats

1. **Pseudo-boards, not the real boards.** The targets are released training stages; the validation answers and test
   inputs arrive at P3. On these boards `neighborhood_mmd`'s floor is 31-139 times its ceiling, which dilutes count
   effects in skill units. On a board whose floor sits closer to its ceiling, the same raw changes cost more points.
   As illustrative arithmetic only: at 8x the raw value sits 0.5 to 0.9 ceiling units above the ceiling, which would
   cost about 1.2 to 2.4 points if floor/ceiling were 10.
2. **Only the two extremes of prediction quality were tested:** the same stage (perfect distribution) and the preceding
   stage (floor quality). Real predictions lie in between. The recommendation of at least 1x and about 2 to 4x holds at
   both extremes; the sign at 8x does not.
3. **Skills are clipped at 1**, so differences near the ceiling are invisible in scores. Raw values and gap fractions
   are in the CSVs. The structure gate is the emulated floor-normalised form; it never bit in the random or jitter
   conditions.
4. **Emulator seeds and working copies are not the server's.** Relative effects should transfer; absolute scores
   should not be read as leaderboard predictions.
5. **The jitter is isotropic in 3-D**, but z is a section index on a per-stage unknown scale. Real position errors are
   anisotropic, and real cells never sit between sections.
6. **Truth halves here have 1,241 to 2,687 cells.** For a much smaller truth half (the embryo validation board had
   [redacted 2026-10-19: A]), `mmd_u` (a 2,000-cell draw) and `variogram` (1,500 cells) add estimator noise at small absolute counts. That
   should strengthen the case for at least 2x, but it was not tested.
7. **Not tested:** blue-noise or even spacing (the PLAN says it helps Dice), duplicated cells, deformations, and the T1
   and T3 boards, which have no spatial terms.
8. **The first pass's per-process logs were overwritten** when the run resumed to add refblock and refhead.
   `runs/logs/density_*.log` therefore shows only the second pass; every call's result and timing is in
   `runs/density/raw/*.jsonl`.
9. **The U-shape mechanism** (overlapping-neighbourhood bias in the unbiased MMD) is an interpretation and was not
   tested directly.

## Recommendations for docs/PLAN.md (and ASSUMPTIONS.md)

1. **§3, principle 6 (spatial hygiene).** Replace "Submitted density near the truth half's density (neighborhood_mmd
   uses 15 nearest submitted cells)" with:
   > Submit about 2x truth_cells. 1x to 4x is a plateau. Never go below 1x: on emulator pseudo-boards 0.5x cost 1-5
   > points and 0.25x cost 2-10, because occupancy_dice thins the truth to the submitted count and neighbourhoods
   > over-smooth.

   Keep "never contiguous rows" and give its price:
   > Rows are stored in z-section order; a contiguous block costs 16-28 points (mid-file) or 21-48 (first rows) at 1x
   > on a near-perfect prediction, and still 12-21 at 8x.

   Add:
   > Jitter of at most 0.25 x the median NN distance (to break duplicate coordinates) is free. Keep per-cell position
   > noise below about 1 NN spacing (6-9% of the RMS radius at 5% sampling); 2 NN costs 1.2-1.8 points at 1x.
2. **§4.2, step 6.** Change to "Cell count about 2x the test truth half (never below 1x; within [min_cells,
   max_cells])", keeping the optional blue-noise clause.
3. **§4.3, pipeline step 1.** Change "sampled at random down to about the test truth-half size" to "sampled at random
   to about 2x the test truth-half size (never below 1x; within max_cells)".
4. **§7 timeline.** Mark the density experiment done (2026-10-18, `experiments/density/REPORT.md`) and add a P3
   follow-up. When the validation answers land, add E10.5 (reference E9.5) and E7.5 (reference E7.25) to `TARGETS` in
   `run_density.py`, one line each, and rerun with `run_all.sh`; it takes about 10 minutes on CPU. Those boards'
   floor/ceiling gaps will be narrower than here, so they give the realistic size of the 4x to 8x effects before the
   test count is fixed.
5. **`ASSUMPTIONS.md`, A13.** Set the status to "passed with modification: near-optimal window 1-4x truth_cells, best
   about 2x, never below 1x; recheck on P3 validation targets", with evidence `experiments/density/REPORT.md`. Edit
   the test text from 0.25x-4x to 0.25x-8x and from "optimum within 0.5-2x" to "optimum within 1-4x".
6. **§6 infrastructure (pseudo-task harness).**
   - Compare candidates on paired seeds, at least 5. The emulator's seed sd is 0.3 to 1.3 points even for a perfect
     prediction, so differences under about 1 point are noise.
   - Near the ceiling, compare raw metric values or gap fractions as well as the clipped skills (`summarise.py` shows
     one way to do it).
   - Add a guard to the submission builder: the submitted cloud's RMS radius and occupied-voxel fraction (16³
     canonical grid) must match the intended source cloud within a few percent. At 1x, a contiguous block shows an
     RMS ratio of 0.64 to 0.96 and covers 15 to 25% of the occupied voxels; a random subset sits at an RMS ratio of
     about 1.0 and covers 81 to 83%.

## Files

Local, in `experiments/density/`:
- `run_density.py`: the experiment (bookkeeping, conditions, score_board calls, JSONL output, resumable).
- `run_all.sh`: the DGX launcher, one process per (target, seed).
- `diag_order.py`: structure-only diagnostics (row order, block coverage, NN spacing).
- `summarise.py`: long CSV, per-condition mean/sd, paired deltas, pooled table and plots.
- `report_tables.py`: the markdown tables in this report, built from the summary CSV (stdlib only).
- `results/`:
  - `results_long.csv` (720 rows), `summary_by_condition.csv`, `summary_pooled.csv`, `summary_tables.txt`,
    `diag_order.csv`
  - `score_vs_count.png`, `gap_fraction_vs_count_random.png`, `jitter.png`

DGX, in `~/v/virtual-embryo-challenge/`: `runs/density/raw/*.jsonl` (15 files, 48 calls each),
`runs/density/summary/`, `runs/density/diag_order.csv`, `runs/density/summary_tables.txt` and
`runs/logs/density_*.log`. No competition data left the DGX, and nothing was committed to git.

## Independent review (2026-10-19)

An independent agent re-ran the analysis on the DGX with its own code. Its verdict, unedited:

The density experiment's numbers are right, and most of its conclusions hold. Its main recommendation, 'submit about 2x truth_cells, window 1x-4x', does not survive more seeds and more realistic prediction quality.

**What I checked**
- **Reproduction:** with independently written bookkeeping, 315 stored (target, seed, condition) rows reproduce bit-exactly.
- **Re-aggregation:** recomputing all 720 raw rows matches every table I checked.
- **Leakage:** none. Prediction cells never overlap the truth or ceiling halves, and no pool cell shares coordinates with a working-copy cell.
- **Setup:** panels, the floor, computed anchors and seeds are all wired correctly.

**What holds up**
- **Below 1x:** never submit fewer cells than truth_cells. All 20 seeds lose at 0.5x on every target, including a 665-cell truth target (-7.1 points) and an ~8k truth half (-0.9 points).
- **Contiguous rows:** the penalty is large (16-28 points mid-file, 21-48 for the first rows at 1x; still 12-21 at 8x).
- **Jitter:** about 0.25 NN is free; 2 NN costs about 1.2-1.7 points at 1x.
- **Denser jitter tolerance:** a denser submission tolerates the same physical jitter better on the heart targets, but not significantly on the embryo targets.

**What does not hold up**
- **The 2x recommendation.** I ran 20 seeds, added the small-truth target E7.25 (reference E6.75) and a prediction of intermediate quality (50% of expression rows permuted). Near-perfect predictions are flat from 2x to 8x, and 1x is only about 0.3 points worse. Imperfect predictions improve steadily with count:
  - floor-quality: 2x is worse than 4x by 0.35 ±0.17 points and worse than 8x by 0.56 ±0.23;
  - intermediate quality: 1x->8x gains +1.63 points pooled, and only 2 of 60 seeds lose;
  - with an ~8k truth half (score_board frac=0.3, [redacted 2026-10-19: A]), intermediate quality still gains +0.96 at 3x, which is that board's max_cells.
- **The raw neighborhood_mmd minimum at 2x.** It holds only for the truth sizes tested (1.2-2.7k cells). It moves to 4x for a 665-cell truth and to 1x for an ~8k truth, and it never shows in the score because skills are clipped at 1.
- **A13.** 'Matching truth-half density is near-optimal' fails for realistic, imperfect predictions.

**Smaller issues**
- The refblock centre is effectively fixed on heart_E9.5 and identical on embryo_E8.0, so the 'lucky crops' are single slabs, not seed luck.
- E6.75 was never checked for row order, and E9.5 and E6.75 are not strictly section-ordered.
- The fresh-reference de_score finding reverses on E7.25.
- C2's 'all paired changes within ±0.8' and 'never worse than 1x' are worded too strongly.
- The voxel guard in recommendation 6 is mis-specified.

**Data rules:** no violations. The lead should confirm whether keying the submission count to index.json truth_cells is acceptable; it is a board-spec field, not an anchor.

**Recommended changes to docs/PLAN.md**
- **§3.6, §4.2 step 6, §4.3 step 1:** 'never below 1x truth_cells; prefer the largest count allowed up to about 4-8x truth_cells, capped at max_cells (the heart extrapolation board caps at 3x)', rather than 'about 2x'.
- **A13:** mark as 'fails as stated: count below 1x is heavily penalised and above 1x is a lever for imperfect predictions'. Recheck on the P3 E10.5 and E7.5 targets.
- **Contiguous rows and jitter:** keep the report's rules and their costs.

**Files**
- Scripts, raw results and logs on the DGX: ~/v/virtual-embryo-challenge/runs/density_review/ (raw/ holds 1,460 calls, raw_big/ 200 calls, plus analysis.txt and analysis_big.txt), with logs at runs/logs/density_review_*.log. They were run with tmux sessions vec-density-review and vec-density-review-big via review_all.sh and review_big_all.sh, then review_analyse.py, review_reagg.py and review_structure.py.
- Local copies of the scripts and analysis files: /tmp/claude-1001/-mnt-md0-IITM-BackUp-Home-vijayavallabh/07d440b9-3ae8-4421-8a08-a046b42f307d/scratchpad/density_review/

### Checks

| # | Item | Verdict | Notes |
|---|---|---|---|
| 1 | Stored scores can be reproduced with independently written bookkeeping | confirmed | My own index derivation (direct numpy RNG, checked for equality against vec.evaluate.subsample/split_half on the real stages) re-scored 315 stored (target, seed, condition) rows covering 21 condition types, from both the first and second pass (random, refstage, block, head, refblock, refhead, jitter). Max /score diff/ = 0, max /skill diff/ = 0, and n_pred is identical. The DGX copies of run_density.py, summarise.py, diag_order.py, report_tables.py and run_all.sh have the same md5 as the local ones. |
| 2 | Report tables match the raw JSONL (720 rows) | confirmed | I re-aggregated runs/density/raw/*.jsonl without summarise.py. I found 720 unique rows, no duplicates and no inverted anchors. Every number I checked matches: the count sweep (0.25x -5.10/-9.52/-5.42, 0.5x -1.77/-4.66/-2.20), refstage deltas, losing-seed counts (5/15 at 2x, 5/15 at 4x, 3/15 at 8x), block/head/refblock/refhead deltas (pooled at 1x: -23.1/-38.6/-5.5/-9.9), all jitter cells and pooled jitter rows, nbhd raw/ceiling (minimum at 2x in all 15 seeds), the anchors table, and dice voxel/NN at 1x of 5.8/4.0/6.0. |
| 3 | No leakage between truth/ceiling cells and prediction cells | confirmed | The pool excludes the whole 10% working copy, so prediction cells never overlap the truth or ceiling halves, and this is asserted on every call. I added a content-level check: 0 pool cells have coordinates identical to any working-copy cell, for targets and references alike. No released WT T2 file contains duplicate coordinates. |
| 4 | Panels, floor, anchors and seeds wired correctly | confirmed | Heart targets use the 500-gene T2:heart:val_extrap panel. embryo_E8.0 uses the 498-gene embryo panel, with E8.0 subset to it. anchors='computed'; the floor is score_board's own reference working copy (seed s+1); split_half uses seed 0 as in score_board; one process per (target, seed), seeds 0-4. The board is used only for its gene list. truth_cells is computed as round(0.1N)//2. |
| 5 | Cliff below 1x (all seeds lose at 0.25x and 0.5x, for near-perfect and floor-quality predictions) | confirmed | With 20 seeds, near-perfect 0.5x vs 1x is -1.70/-4.26/-2.69 (20/20 seeds lose on each target), and floor-quality 0.5x is -1.09/-1.97/-0.79 (58/60 lose). On a small-truth target (E7.25 from E6.75, truth 665): -7.13. With an ~8k truth half (E9.5, score_board frac=0.3): -0.87 (10/10). The cliff is robust and shrinks as the truth half grows. |
| 6 | 1x-4x plateau for near-perfect predictions | confirmed | With 20 seeds, pooled change vs 1x is 2x +0.27 ±0.23, 4x +0.27 ±0.23, 8x +0.14 ±0.29; 2x minus 4x = -0.00 ±0.21. Wording problems: C2's 'all paired changes within ±0.8 points' holds only for target means, since per-seed paired changes run from -1.69 to +2.74. 'Never worse than 1x' is false in the report's own data: 6/15 seeds lose at 2x, and target means are -0.01 and -0.13. |
| 7 | About 2x truth_cells is the best single count; A13 window 1x-4x | refuted | This only holds for near-perfect predictions, where 2x, 4x and 8x tie. For imperfect predictions the score keeps rising with count. Floor-quality (refstage, 20 seeds): 2x is worse than 4x by -0.35 ±0.17 (48/60 seeds) and worse than 8x by -0.56 ±0.23. A new intermediate-quality condition (random subset with 50% of expression rows permuted, nbhd skill 0.59-0.86): 1x->2x +0.69, ->4x +1.17, ->8x +1.63 pooled (2/60 seeds lose at 8x). Small truth (E7.25): even near-perfect predictions gain at 4x (+1.17 ±0.40, 0/20). At ~8k truth: intermediate 3x vs 1x +0.96 ±0.29 (0/10). The evidence favours the largest count allowed up to about 8x, not 2x. |
| 8 | Raw neighborhood_mmd is U-shaped with its minimum at 2x on all targets | partially-confirmed | True for the three tested truth sizes (1,241-2,687 cells; every seed). The minimum moves to 4x for a 665-cell truth (E7.25) and to 1x for an ~8k truth (E9.5 frac 0.3: 1.29x ceiling at 1x, 1.58 at 2x, 3.04 at 3x). The 2x minimum is tied to the truth sizes tested (near the 2,000-profile cap). It is also invisible in score because skills are clipped at 1. |
| 9 | Floor-quality predictions gain with count (+0.4 to +1.6 at 8x) | confirmed | Stronger with 20 seeds: 8x vs 1x is +0.86/+1.94/+0.53, and +3.01 on E7.25. |
| 10 | Rows in every released T2 file are sorted by z section | partially-confirmed | Spearman(row, z) is -0.9994 to -0.9996 in all six files, including E6.75, which diag_order.py never checked. Only E7.25, E8.0, E8.25_late and E8.75 are strictly section by section (z-runs = z-levels). E9.5 has 3,484 runs over 98 levels and E6.75 has 988 over 65. The T3 KO file was not checked. |
| 11 | Contiguous-block penalty (1x: 16-28 points mid-file, 21-48 points first rows; 8x: 12-21 points) | confirmed | Reproduced exactly; all 15 seeds lose at every count. Caveat: because the centre must fit the 8x block, the refblock centre is nearly fixed on heart_E9.5 (block_start 9,425-10,105 of 22,343) and identical on embryo_E8.0 (5,191 in every seed). The 'lucky crops' (E9.5 refblock 4x +1.8 with 0/5 seeds losing; E8.0 refhead 4x +1.3) are therefore one fixed slab, not seed luck. |
| 12 | Jitter curve (free up to 0.25 NN; 2 NN costs 1.2-1.8 points at 1x) | confirmed | With 20 seeds at 2 NN, 1x/own costs -1.35/-1.68/-1.19 (59/60 seeds lose), and -2.93 on the small-truth E7.25. |
| 13 | A denser submission tolerates the same physical jitter better | partially-confirmed | With 20 seeds, loss(4x/truth, f=2) minus loss(1x/own, f=2) is +0.46 ±0.28 (E9.5) and +0.82 ±0.45 (E8.75), but not significant on embryo E8.0 (+0.28 ±0.38) or E7.25 (+0.55 ±0.72). The 5-seed E8.0 reversal (-1.36 vs -1.17) was noise. |
| 14 | Side finding: fresh reference cells score below the floor on de_score | partially-confirmed | de skill is 0.41-0.45 on the three targets, but 0.50-0.53 on E7.25 (reference E6.75). de_direction skill is above 0.5 everywhere (0.51-0.63). The effect depends on the stage pair. |
| 15 | Mechanism of the loss below 1x (Dice count-matching; voxel edge 4-6 truth NN) | confirmed | veckit git main (cc3859a, 0.1.1) _match_n subsamples both clouds to min(n_pred, n_true). At 0.25x the shape group carries most of the loss (-2.9/-4.8/-4.1 points vs nbhd -1.4/-3.3/-0.8 and expression -0.8/-1.2/-0.5). voxel/NN is 5.8/4.0/6.0 at 1x and 8.1 at an ~8k truth. |
| 16 | Emulator seed noise 0.3-1.3 points for a perfect 1x prediction | partially-confirmed | With 20 seeds the random_1x sd is 0.40/1.13/1.47, so the range reaches about 1.5. |
| 17 | Recommendations transfer to the real heart boards | unverified | The experiment's truth halves are 1.2-2.7k cells. The heart val_extrap spec has [redacted 2026-10-19: A] (=3x truth), so the 1x-4x window partly exceeds max_cells there. My frac=0.3 emulation ([redacted 2026-10-19: A]) suggests: near-perfect predictions are flat over 1x-3x, imperfect ones still gain up to 3x. |
| 18 | Compliance with data rules and run hygiene | confirmed | The experiment used only released T2 files and the pinned git-main veckit through vec.evaluate. No index.json anchors are used. Local results are metric summaries only (no h5ad or matrices). experiments/density is untracked in git. DGX writes stayed in runs/density and runs/logs; the tmux session was vec-density; CPU only. My own review wrote only to ~/v/virtual-embryo-challenge/runs/density_review and runs/logs/density_review_*.log; tmux sessions vec-density-review and vec-density-review-big have exited. |

### Code issues found by the review

1. run_density.py conditions(): the block centre is drawn from [n_max//2, n_pool-(n_max-n_max//2)] so that the 8x block fits. For the reference pools this range collapses: heart_E9.5 refblock centre lies in [10748, 11595] (block_start 9,425-10,105 across seeds), and embryo_E8.0 refblock is fixed at 5,982 (block_start 5,191 in every seed). refblock is therefore effectively deterministic, its seed sd understates position variability, and the 'lucky crop' consistency (0/5 seeds losing) reflects one slab. Draw a centre per ratio, or cap the nest at the largest block that leaves room.
2. The 'random' centre for target blocks is also restricted to the middle of the file (for example E9.5 block_start 17,128-35,726 of 48,368). That is fine for 'mid-file', but it is not a uniformly random position.
3. diag_order.py STAGES omits E6.75 (and the T3 KO), yet the report says 'every released T2 file'. E9.5 and E6.75 are not strictly section-ordered (3,484 and 988 z-runs). Block coverage numbers come from a single seed-0 mid-pool block.
4. Report wording, C2: 'all paired changes within ±0.8 points' holds only for means (per-seed range -1.69 to +2.74), and 'never worse than 1x' is contradicted (6/15 seeds lose at 2x; target means -0.01 and -0.13).
5. Report wording: '44% of the rows' at 8x is 44% of the pool (40% of the file).
6. Recommendation 6 (guard) says the occupied-voxel fraction must match the source cloud 'within a few percent', but a random 1x subset itself covers only 81-83% of the source's voxels. The guard needs a same-count random reference, and RMS alone misses mid-file blocks (ratio 0.93-0.96).
7. Recommendation 4 asserts that P3 boards will have narrower floor/ceiling gaps, with no supporting data. Larger truth halves shrink the split-half ceiling: refstage nbhd/ceiling is about 264 at an 8k truth vs 73-79 at 2.7k.
8. summarise.py: seeds_dscore_neg is computed against random_1x for reference-drawn conditions too, which is meaningless there. The refstage seed counts in the text (5/15, 5/15, 3/15) were computed elsewhere but are correct.
9. None of these affect the stored numbers: every stored row I re-scored reproduces bit-exactly.

### Compliance notes from the review

1. No rule violations found. The experiment uses only released files and anchors='computed', computes truth_cells from the stage, never opened veckit 0.1.2 or PR #6 (the scorer is git main cc3859a), copied no competition data locally and committed nothing.
2. Note: caveat 6 quotes the embryo validation board's truth_cells ([redacted 2026-10-19: A]) from boards/index.json. truth_cells is not an anchor and is not used to set any number, but it is an index.json field that encodes a held-out stage property.
3. Note for the lead: the recommendations set the test submission count as a multiple of the board's published truth_cells from index.json. truth_cells is a board-spec field, not an anchor, and PLAN §4.2/§4.3 already rely on it, but it encodes the held-out stage's cell count. Confirm this is acceptable under the team rule, or key the count to max_cells instead (max_cells >= 3x truth by the board rule).
4. Note: recommendation 4's claim about narrower floor/ceiling gaps on P3 boards has no stated source. Make sure it is not based on index.json anchors.
5. Disclosed: the first-pass per-process logs were overwritten. Per-call results survive in runs/density/raw.

## Lead's decisions (2026-10-19)

- **Cell count: max_cells, never below truth_cells.** The review's 20-seed runs, with a small-truth target, an 8k-cell
  truth half and an intermediate-quality prediction, show imperfect predictions gaining with count, so the experiment's
  "about 2x" recommendation is withdrawn. `ASSUMPTIONS.md` A13 is revised.
- **Random rows only; jitter at most 0.25 NN.**
- **Compliance:** keying the count to max_cells means truth_cells does not choose it. The earlier organiser question on
  truth_cells stands.
- **Applied by the follow-up run** (see "Revisions" below once added): the corrected block-centre draw, the
  intermediate-quality condition, count levels as fractions of max_cells, release-day targets E10.5 and E7.5.
