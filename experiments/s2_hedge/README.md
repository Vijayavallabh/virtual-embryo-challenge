# s2 hedge payoffs (sensitivity analysis for the s2 pre-registration)

**Not an adoption result and not a choice.** Nothing here picks an entry. The s2 rules are pre-registered on 17-23 Nov
with the captain's sign-off (`docs/PLAN.md` sections 7-8). These runs measure what a second upload gains where the
test's tissue mix is not the one s1 assumes, and what it costs where it is.

Each board allows 2 scored uploads, and the captain marks the better-scoring one as final, so s2 is worth most when it
differs from s1 along the axis of greatest uncertainty. The method: re-sample a released pseudo-task's truth to a grid
of tissue mixes, build the prediction at a grid of mixes, and score every prediction against every truth with the
harness's scoring (bake protocol, scorer build main, 10% working copies, each truth's floor and ceiling computed per
seed and shared by the predictions, build seed 1,000,000 + seed). Settings were fixed before each run
(`payoff.py`); truth re-sampling uses one fixed seed per grid value.

## T1 (`T1/E9.5->E10.5:chamber`; 10 seeds; commit cf476aa)

Both one-axis runs build the prediction from E9.5's own chamber mix (blood 2.6% of chamber cells, cardiomyocytes 66%
of heart cells) and move one share. B0 at its own mix against the unchanged truth equals the harness's B0 cell for
cell (d2_v1). Mean board scores (floor 50):

Blood's share of chamber cells (`results/t1_blood/`):

| prediction \ truth | own (22.9%) | 2.8% | 7.6% | 14.3% | 33.3% |
|---|---|---|---|---|---|
| own (2.6%) | 49.68 | 54.15 | 52.98 | 51.53 | 48.47 |
| 3% | 50.03 | 54.38 | 53.22 | 51.84 | 48.80 |
| 8% | 53.52 | 55.89 | 55.47 | 54.91 | 52.72 |
| 15% | 58.35 | 56.49 | 57.14 | 58.24 | 58.40 |
| 35% | 62.27 | 50.34 | 52.70 | 57.51 | 67.20 |

Cardiomyocyte share of heart cells (`results/t1_cm/`; the column is the truth's share of all cells):

| prediction \ truth | own (0.191) | 0.108 | 0.288 | 0.418 | 0.541 |
|---|---|---|---|---|---|
| own (66%) | 49.68 | 48.41 | 52.57 | 55.83 | 59.19 |
| 15% | 55.14 | 58.39 | 51.39 | 46.78 | 44.32 |
| 40% | 55.94 | 56.34 | 56.12 | 54.35 | 53.40 |
| 58% | 52.35 | 51.41 | 54.63 | 56.42 | 58.43 |
| 75% | 46.73 | 45.32 | 49.77 | 54.12 | 58.86 |

Reading: composition moves T1 by up to 12-19 points on this task, far more than any expression step (M1's direction
added 1.29). On the cardiomyocyte axis the matched prediction wins every column, and a mismatch of 0.6 costs about 14
points. **Both axes are confounded:** every blood row also carries E9.5's cardiomyocyte excess (adding blood dilutes
it, which is why 15% blood wins even against a 3%-blood truth), and every cardiomyocyte row carries E9.5's blood
shortfall. At test, s1 and s2 share E10.5's own mix and differ only on the hedged shares, so these tables are context,
not the decision basis.

`t1-grid` removes the confound (10 seeds; commit 90b25d3; `results/t1_grid/`). Every prediction starts from the
truth's own chamber mix (an oracle base) and sets two shares:
- blood: own (22.9%), 5% or 40% of chamber cells;
- cardiomyocytes: 15%, own (26.5%), 45%, 60% or 75% of heart cells.

The truths use the same 15-point grid; labels read `b<blood>_c<cardiomyocytes>`. The full 15 x 15 table is in
`results/t1_grid/payoff.txt`. The blood-own slice (mean board score; s1's analogue is `bown_cown`):

| prediction \ truth | c0.15 | cown | c0.45 | c0.6 | c0.75 |
|---|---|---|---|---|---|
| bown_c0.15 | 71.03 | 68.06 | 62.66 | 57.13 | 53.33 |
| bown_cown (s1) | 71.15 | 69.92 | 66.72 | 61.66 | 57.77 |
| bown_c0.45 | 66.90 | 68.34 | 70.62 | 68.44 | 65.64 |
| bown_c0.6 | 61.56 | 64.15 | 69.43 | 71.22 | 71.05 |
| bown_c0.75 | 55.81 | 58.52 | 64.81 | 69.38 | 72.68 |

The blood slice at the own cardiomyocyte share costs s1 much less. Against 5% blood it scores 64.04, where the matched
prediction scores 66.07; against 40% it scores 70.13, against 74.65.

The table below is the hedge (`results/t1_grid/hedge_bown_cown.*`): the expected better of s1 and each s2, paired by
seed, with each candidate's gain over s1 alone. Two priors come from `hedge`: uniform over the 15 truths, and one
third on own. Three more were computed from the same runs:
- blood own, cardiomyocytes in {own, 0.45, 0.6, 0.75};
- cardiomyocytes at 0.45 or above;
- cardiomyocytes at own, with blood varying.

| s2 | uniform | 1/3 on own | blood own, cm own..0.75 | cm >= 0.45 | cm own, blood varies |
|---|---|---|---|---|---|
| bown_c0.6 | 68.28 (+4.74) | 68.75 (+3.39) | +6.40 | +7.90 | +0.00 |
| bown_c0.75 | 67.89 (+4.35) | 68.47 (+3.11) | +5.68 | +7.25 | +0.00 |
| bown_c0.45 | 66.91 (+3.37) | 67.77 (+2.41) | +4.64 | +5.62 | +0.03 |
| b0.4_c0.6 | 66.50 (+2.96) | 67.48 (+2.12) | +3.16 | +4.87 | +0.24 |
| b0.4_cown | 65.03 (+1.49) | 66.43 (+1.07) | +0.03 | +1.59 | +1.51 |
| b0.05_cown | 63.99 (+0.45) | 65.68 (+0.32) | +0.00 | +0.21 | +0.69 |
| none (s1 alone) | 63.54 | 65.36 | 64.02 | 60.14 | 68.03 |

Reading:
- **The cardiomyocyte share is the lever.** A mismatch there costs up to 12-19 points, while blood costs s1 2-5.
- **bown_c0.6 is the best s2** (blood own, cardiomyocytes at 60% of heart cells) under every prior that lets the
  cardiomyocyte share move.
- **A blood hedge wins only if the cardiomyocyte share is known to stay at E10.5's,** and then by 1.5 points at most.
  This replaces the expectation in `docs/PLAN.md` 4.1 that the hedge belongs on blood.
- **Where the s2 lands:** the legal E14.5 datasets put the cardiomyocyte share at 40-77% of heart cells (a different
  assay, so the levels do not transfer), against E10.5's 26.5%. An s2 at 60% covers that upper range, and s1 keeps the
  same-assay estimate.

Caveats:
- An oracle base: every row starts from the truth's own chamber mix.
- One transition, E9.5 -> E10.5.
- The priors are stated assumptions.

## T2 heart (`T2h/<=E9.5->E10.5`; `heart-scope`; 10 seeds; commit 90b25d3; `results/heart_scope/`)

Prompted by a side note (2026-11-01): the planned heart s2 (scope hypothesis H2) uses the same cells as s1 and is only
9% smaller, while the E10.5 heart sample is mostly surrounding tissue (cardiomyocytes 15% of cells; the core, CM plus
endocardium, 23%). If E12.5 is a cleaner heart, both uploads miss the same way. Truths: E10.5 as released, and E10.5
with surrounding cells removed until its in-plane scope (log of whole-cloud over core RMS, the size procedure's scope)
is half or none of its own, farthest from the core first ("radial") or at random ("thin"). Predictions: B0 from E9.5
held ("hold", the s1 analogue), shrunk in plane by half its scope ("half_size", today's size-only H2), or trimmed to
half its scope or to none before the draw ("half_radial", "half_thin", "core_radial"); every prediction is scaled in
plane to the truth's core size (an oracle core size), so rows differ in scope only.

Mean board scores (floor 50). Brackets give the gain over hold and the number of seeds, out of 10, where the row
beats hold. The column headings give each truth's core share:

| prediction \ truth | own (23%) | radial_0.5 (29%) | radial_0 (40%) | thin_0.5 (85%) | thin_0 (100%) |
|---|---|---|---|---|---|
| hold (s1) | 50.22 | 52.09 | 47.46 | 52.05 | 43.63 |
| half_size (H2 today) | 49.93 [-0.30, 0] | 50.87 [-1.22, 0] | 50.20 [+2.74, 10] | 53.59 [+1.55, 10] | 44.73 [+1.10, 10] |
| half_radial | 44.34 [-5.88, 0] | 46.60 [-5.49, 0] | 47.77 [+0.31, 6] | 55.99 [+3.95, 10] | 49.40 [+5.77, 10] |
| half_thin | 27.23 [-23.00, 0] | 28.34 [-23.75, 0] | 30.33 [-17.13, 0] | 55.07 [+3.02, 10] | 51.20 [+7.57, 10] |
| core_radial | 38.93 [-11.29, 0] | 39.96 [-12.12, 0] | 42.75 [-4.71, 0] | 54.03 [+1.98, 10] | 54.47 [+10.84, 10] |

The hedge (`results/heart_scope/hedge_hold.*`) takes the better of s1 and the s2 for each seed and truth, then
averages over the truths. The two family columns average only over the truths named; s1 alone scores 49.92 on the
radial family and 47.84 on the thin family.

| s2 | uniform over 5 truths | 1/3 on own | radial family (own, radial_*) | thin family (thin_*) | worst truth |
|---|---|---|---|---|---|
| core_radial | 51.65 ± 0.32 | 51.42 ± 0.33 | 49.92 | 54.25 | radial_0 (47.46) |
| half_thin | 51.21 ± 0.42 | 51.04 ± 0.38 | 49.92 | 53.13 | radial_0 (47.46) |
| half_radial | 51.12 ± 0.29 | 50.97 ± 0.29 | 50.07 | 52.70 | radial_0 (47.89) |
| half_size (H2 today) | 50.17 ± 0.43 | 50.18 ± 0.37 | 50.84 | 49.16 | thin_0 (44.73) |
| none (s1 alone) | 49.09 ± 0.47 | 49.28 ± 0.40 | 49.92 | 47.84 | thin_0 (43.63) |

Reading:
- **H2 today is a weak hedge.** As s2 it adds about 1 point (+1.08 uniform, +0.90 with a third on own). It stays close
  to s1, and its worst case is still a clean heart (thin_0, 44.73).
- **A core-trimmed s2 is the strongest hedge here,** adding 2.1 to 2.6 points (core_radial). It covers the truths
  where s1 is weakest: against a clean heart it gains 10.8 points on expression and shape alike (mean skill: de_score
  0.59 against 0.52, mmd_u 0.53 against 0.44, d2_shape 0.61 against 0.40).
- **It would be a bad s1.** Against E10.5 as released it loses 11.3 points, mostly on the expression-spread metrics
  (mmd_u 0.30 against 0.47, neighborhood_mmd 0.36 against 0.52).
- **The choice depends on the prior.**
  - On the radial family only H2 helps (+0.92); core_radial never beats hold there.
  - On the thin family core_radial adds +6.41 against H2's +1.32.
  - half_radial sits between them (+0.15 and +4.86).
  - The released stages lean radial: the core share fell from 37% at E9.5 to 23% at E10.5.
- **Random thinning (half_thin)** collapses on the wide truths (down 17 to 24 points) and is dominated.

Caveats:
- The truths are synthetic trims of E10.5, and the priors are stated assumptions, not estimates.
- Every row uses an oracle core size; at test the size comes from the pre-registered procedure.
- E9.5 -> E10.5 is one transition.

A heart-focused s2 would be a shape operation, which `docs/PLAN.md` 4.3 step 5 admits only if it gains on both
released transitions; on E9.5 -> E10.5 the core share fell (37% to 23%), so it would need the captain's ruling. These
numbers inform that ruling; they change nothing in the plan.

## Runs

The first launch of `t1-grid` and `heart-scope` (2026-11-01 13:32 UTC, 46089fe) stalled in kernel memory compaction
(`docs/COMPUTE.md`) and was stopped after 6.9 hours. By then 2 heart seeds (0 and 5) had finished; their log lines are
kept in `runs/s2_hedge/logs/stuck_2026-11-01/` on the DGX as a check that the rerun reproduces them. The relaunch adds
`NUMPY_MADVISE_HUGEPAGE=0`. The heart scope runs as 5 processes of 2 seeds each. Each grid process now scores each
prediction straight after building it, instead of holding all 15 (about 20 GB per process). The scores are the same,
because every builder and scorer call is seeded on its own. That claim comes from reading the code (no global RNG in
`vec` or the scorer modules), not from a rerun, since the stalled grid finished no seed.
