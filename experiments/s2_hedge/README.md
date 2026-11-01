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

`t1-grid` (running) removes the confound: every prediction starts from the truth's own chamber mix (an oracle base) and
sets blood (own, 5%, 40% of chamber cells) and cardiomyocytes (15%, own, 45%, 60%, 75% of heart cells); truths use the
same 15-point grid.

## T2 heart (`T2h/<=E9.5->E10.5`; `heart-scope`, running)

Prompted by a side note (2026-11-01): the planned heart s2 (scope hypothesis H2) uses the same cells as s1 and is only
9% smaller, while the E10.5 heart sample is mostly surrounding tissue (cardiomyocytes 15% of cells; the core, CM plus
endocardium, 23%). If E12.5 is a cleaner heart, both uploads miss the same way. Truths: E10.5 as released, and E10.5
with surrounding cells removed until its in-plane scope (log of whole-cloud over core RMS, the size procedure's scope)
is half or none of its own, farthest from the core first ("radial") or at random ("thin"). Predictions: B0 from E9.5
held ("hold", the s1 analogue), shrunk in plane by half its scope ("half_size", today's size-only H2), or trimmed to
half its scope or to none before the draw ("half_radial", "half_thin", "core_radial"); every prediction is scaled in
plane to the truth's core size (an oracle core size), so rows differ in scope only.

A heart-focused s2 would be a shape operation, which `docs/PLAN.md` 4.3 step 5 admits only if it gains on both
released transitions; on E9.5 -> E10.5 the core share fell (37% to 23%), so it would need the captain's ruling.
