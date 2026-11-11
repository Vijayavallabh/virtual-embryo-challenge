# s2 hedge payoffs (sensitivity analysis for the s2 pre-registration)

**Not an adoption result and not a choice.** Nothing here picks an entry. The s2 rules are in `docs/S2_RULES.md`,
written by the AI assistant under the captain's delegation of 2026-11-02 and frozen on 2026-11-03. These runs measure what a second upload gains where the
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

Blood's share of chamber cells (`results/t1_blood/`; the column is the truth's share of all cells, which include 5%
surrounding cells):

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
added 1.29). On the cardiomyocyte axis the matched prediction wins three of the four columns that have one (in the 75% column the
own 66% row edges it, 59.19 against 58.86), and a mismatch of 0.6 costs about 14 points. **Both axes are confounded:** every blood row also carries E9.5's cardiomyocyte excess (adding blood dilutes
it, which is why 15% blood wins even against a 3%-blood truth), and every cardiomyocyte row carries E9.5's blood
shortfall. At test, s1 and s2 share E10.5's own mix and differ only on the hedged shares, so these tables are context,
not the decision basis.

`t1-grid` removes the confound (10 seeds; commit 90b25d3; `results/t1_grid/`). Every prediction starts from the
truth's own chamber mix (an oracle base) and sets two shares:
- blood: own (24.1% of chamber cells; 22.9% of all truth cells), 5% or 40% of chamber cells;
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
- **The cardiomyocyte share is the lever.** With blood at its own share, a mismatch there costs s1 up to 15 points
  (3.9, 9.6 and 14.9 against truths at 45%, 60% and 75% of heart cells), while blood costs s1 2-5.
- **bown_c0.6 is the best s2** (blood own, cardiomyocytes at 60% of heart cells) under every prior that lets the
  cardiomyocyte share move.
- **A blood hedge wins only if the cardiomyocyte share is known to stay at E10.5's,** and then by 1.5 points (at most
  1.6 under any prior here). This reverses the direction the captain agreed on 2026-11-01 (a hedge on the blood share;
  `docs/PLAN.md` 4.1, decision log 2026-11-01 (d)). It was decided by the AI assistant under his delegation
  (`docs/S2_RULES.md`).
- **Where the s2 lands:** Feng's two strains at E14.5 put the cardiomyocyte share at 42-77% of heart cells (a different
  assay, so the levels do not transfer), against E10.5's 26.5%. An s2 at 60% covers that upper range, and s1 keeps the
  same-assay estimate. (Restated 2026-11-05: a whole-heart figure was withdrawn and its value redacted, since its
  sample is coded Theiler stage 22, which the staging rule keeps out; the s2 stays inside the remaining range, so
  nothing changes; `docs/DATA_POLICY.md` 2026-11-05 (c) and (e).)

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
  - On the radial family H2 helps most (+0.92; half_radial +0.15); core_radial never beats hold there.
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
numbers inform that ruling; they change nothing in the plan. (2026-11-02: under the captain's delegation they were used
to set the heart s2 to the half trim, treated as a scope hedge; `docs/S2_RULES.md`.)

### The other heart transitions (`heart-scope --task`, commit dcf1779)

Two more transitions, each with the same five truths and predictions:
- **T2h/E8.25_late->E8.5** is step 5's second released transition, defined in `payoff.py` as a diagnostic. Its
  predictions are copies, because B0's momentum needs an earlier stage.
- **T2h/<=E8.75->E10.5** is the harness's 1.75-day pseudo-task, built with B0.

Step 5 (`results/heart_scope*/step5.txt`):

| transition | core_radial d2 | core_radial Dice | half_radial d2 | half_radial Dice | step 5 |
|---|---|---|---|---|---|
| E9.5 -> E10.5 | +0.416 | +0.041 | +0.153 | +0.034 | passes |
| E8.25_late -> E8.5 | -0.021 | +0.112 | -0.036 | +0.029 | fails |
| E8.75 -> E10.5 (not a step 5 transition) | +0.141 | -0.013 | +0.043 | -0.019 | fails |

The hedge on board score (`results/heart_scope*/hedge_hold.txt`). Each cell is the expected better of hold and the s2,
uniform over the five truths:

| s2 | E9.5 -> E10.5 | E8.25_late -> E8.5 | E8.75 -> E10.5 |
|---|---|---|---|
| core_radial | 51.65 | 63.59 | 60.26 |
| half_radial | 51.12 | 56.84 | 60.02 |
| half_size (H2 as built) | 50.17 | 53.45 | 54.50 |
| none | 49.09 | 51.61 | 52.73 |

On the released target alone, the half trim beats hold in two of the three transitions: +1.75 on E8.25_late -> E8.5
(paired mean of s2 minus s1, ahead in 9 of 10 seeds; +2.04 as the better of the two) and +1.57 on E8.75 -> E10.5,
against -5.88 on E9.5 -> E10.5. Step 5's shape bar fails. As an s2, the half trim beats H2 as built on every transition
under the uniform, one-third-on-own and thin-family priors. On E9.5 -> E10.5 it ties on the released target (50.22)
and trails under the radial family (50.07 against 50.84). `docs/S2_RULES.md` records the rule, its outcome and the
amendment that follows from this.

## Runs

The first launch of `t1-grid` and `heart-scope` (2026-11-01 13:32 UTC, 46089fe) stalled in kernel memory compaction
(`docs/COMPUTE.md`) and was stopped after 6.9 hours. By then 2 heart seeds (0 and 5) had finished; their log lines are
kept in `runs/s2_hedge/logs/stuck_2026-11-01/` on the DGX as a check that the rerun reproduces them. The relaunch adds
`NUMPY_MADVISE_HUGEPAGE=0`. The heart scope runs as 5 processes of 2 seeds each. Each grid process now scores each
prediction straight after building it, instead of holding all 15 (about 20 GB per process). The scores are the same,
because every builder and scorer call is seeded on its own. That claim comes from reading the code (no global RNG in
`vec` or the scorer modules), not from a rerun, since the stalled grid finished no seed.

## T2 embryo (`experiments/harness/results/emb_shape_v1/`; 10 seeds; commit 90ddaaa)

No candidate built before 2026-11-02 beats fgw in any seed of the embryo pseudo-tasks it runs on (`emb_fgw_v1`;
pseudobulk_shift and psb_size run on E7.5 only), so none of them adds anything as an s2. The shape variants keep fgw's coupling and mixture of real cells and take positions from
one bracket's geometry, at the size rule (`interpolate_clouds(t_pos=0 | 1)`; the run's fgw equals emb_fgw_v1's):

| | fgw | fgw_shape0 (earlier bracket's geometry) | fgw_shape1 (later bracket's geometry) |
|---|---|---|---|
| E7.5 from (E7.25, E8.0), t = 1/3 | 61.12 | 62.95, ahead in 10 of 10 | 59.47, ahead in 0 of 10 |
| E7.25 from (E6.75, E7.5), t = 2/3 | 58.78 | 53.67, ahead in 0 | 57.49, ahead in 0 |
| better of fgw and it, mean of both | 59.95 | 60.87 | 59.95 |

The E7.5 gain is shape alone (d2 +0.067, Dice +0.109); expression is unchanged. Under the rule fixed before the run,
s2 = fgw_shape0 (E7.5's geometry at test).

## Reproducing

The per-seed run files behind the T1 grid and heart tables are committed in `results/per_seed/`. They were copied
from the compute machine's `runs/s2_hedge/` with their provenance: T1 grid at 90b25d3, heart transitions at 90b25d3
and dcf1779, clean trees. The one-axis T1 tables keep only their cf476aa summaries, and the embryo rows are in
`experiments/harness/results/`.
The extra priors quoted here and in `docs/S2_RULES.md` are uniform over the truths each one lists:

```bash
P="uv run python experiments/s2_hedge/payoff.py"
$P hedge --runs experiments/s2_hedge/results/per_seed/t1-grid/* --s1 bown_cown --out experiments/s2_hedge/results/t1_grid \
  --prior blood_own_cm_own_up=bown_cown,bown_c0.45,bown_c0.6,bown_c0.75 \
  --prior cm_0.45_up=bown_c0.45,bown_c0.6,bown_c0.75,b0.05_c0.45,b0.05_c0.6,b0.05_c0.75,b0.4_c0.45,b0.4_c0.6,b0.4_c0.75 \
  --prior cm_own_blood_varies=bown_cown,b0.05_cown,b0.4_cown
$P hedge --runs experiments/s2_hedge/results/per_seed/heart-scope/* --s1 hold --out experiments/s2_hedge/results/heart_scope \
  --prior released=own --prior radial_family=own,radial_0.5,radial_0 --prior thin_family=thin_0.5,thin_0
```

Run these from the repo root. The same heart command reproduces `heart_scope_e825` and `heart_scope_e875`, with
`experiments/s2_hedge/results/per_seed/heart-scope-e825/*` and `.../heart-scope-e875/*` as the runs and the matching `--out`.
`summarise` and `step5` take the same run folders.

## The Task 1 statement rule (`statement_rule.py`, `results/t1_grid/statement_rule.*`; 2026-11-11)

Prices the rule fixed before the screened read of the organisers' Task 1 page statement (`docs/PLAN.md` decision log
2026-11-11 (c)), on the committed per-seed T1 grid above: if the statement gave a mix's shares, the pair the rule leaves
against today's pair (s1 `bown_cown`, s2 `bown_c0.6`), the better of two paired by seed. Gains, uniform over the 15
mixes / a third on own / blood at its own share:

| Statement gives | Rule | Gain |
|---|---|---|
| cardiomyocyte share (exact conversion) | the nearer upload moves to it | +0.51 / +0.36 / +0.53 |
| cardiomyocyte share of all cells, blood unstated (snapped to the grid) | as above, converted with E10.5's own blood | +0.38 / +0.27 / +0.53 |
| blood share | both uploads move to it | +2.30 / +1.64 / 0.00 |
| both | both take the blood share, the nearer one the cardiomyocyte share | +2.90 / +2.08 / +0.53 |
| not taken: cardiomyocyte share, s2 blood at 40% | both move to it, s2 with blood 40% | +1.95 / +1.39 / +0.53 |

Blood 0 and one-sided bounds lie off the grid and are unpriced. Outcome: the read found no such statement, so nothing
changed (decision log 2026-11-11 (d)). Reproduce from the repo root: `python -I -B experiments/s2_hedge/statement_rule.py`.

## Per-set marks (`bundle_sim.py`, `results/bundle_sim/`)

What the s2 set gains if one final mark covers a whole set of boards instead of each board. These are exact
expectations over the committed per-seed payoffs, with the boards independent, the T2 boards averaged and T3 left
out. Each cell gives the gain over s1 alone, per-board marks / per-set marks:

| Weight on each board's released-like truth | s2s of `docs/S2_RULES.md` | milder set (heart H2 as built, embryo redrawn) | T1's s2 alone |
|---|---|---|---|
| even | +7.63 / +5.54 | +5.52 / +5.16 | +4.74 / +4.74 |
| 1/3 | +5.97 / +3.81 | +4.07 / +3.64 | +3.39 / +3.39 |
| 0.6 | +4.00 / +2.09 | +2.51 / +2.14 | +2.03 / +2.03 |
| 1 | +1.06 / +0.00 | +0.17 / +0.00 | +0.00 / +0.00 |

With per-set marks the bold s2s keep a third to a half less value. They are never more than 0.05 behind the milder
set, and with per-board marks they lead by 1.5-2.1, so one set is kept.
