# Tries ledger: every change tried for an entry, per test board

Started 2026-11-03 by the AI assistant under the captain's delegation of 2026-11-02, after the captain asked for
continuous iteration with no fixed timeline and a side agent warned that endless tweaking against the same few
backtests lets lucky fits through even when the adoption rule passes. The rule's 10 seeds rule out draw noise, not a
change that only suits the particular released stages a backtest uses.

From 2026-11-04 the ledger also carries the captain's guiding principle (quoted in `AGENTS.md`): rules 5-8 below and
the headroom tables of where the entries stand against the published baselines, written by the AI assistant under the
same delegation. Rows written before that date stay as they were.

## Rules

1. Every change tried for an entry gets a row here **before its run**: the biological or physical reason, the metrics
   where its gain should show, and what would make a gain a fit to one backtest rather than a real effect. The row is
   committed before the run starts; the row names that commit.
2. A change is kept only if all three hold:
   - the adoption rule passes (`vec.adoption`, `docs/PLAN.md` section 3, item 8);
   - the gain sits on the metrics its reason predicts;
   - no known feature of one backtest explains it. Example: the E7.25 embryo target sits much closer to E6.75 than
     its day says (2026-11-03 diagnostic below), so a change that moves the mixture toward the earlier bracket gains
     there and says nothing about E7.75, where the bracket fraction is 1/2.
3. Changes without a reason that holds at test are pruned before they run, with the reason written down.
4. Every try counts, diagnostics included; each value of a swept parameter is one try. The counts go into every
   adoption report.
5. **The baseline and the question.** Every row, written before its run, names its board's published baseline
   (copy_last on T1, the heart and the embryo; wt_identity on T3) and the entry or candidate it changes. It answers the
   captain's question: what meaningful change beyond the baseline does it predict, and what held-out evidence would
   show it? The change names its facet: which states emerge, their proportions, their expression or, on T2, their
   spatial organisation. The evidence is named in advance: the backtests or folds that count for it under the
   adoption rule, none it was tuned or selected on (on T3 the cross-knockout folds; the Gata4 replicate fold never
   counts). The row also states the gain it expects there; a change expected to gain less than the adoption rule's
   0.3 points on its counting backtests is pruned before it runs (rule 3), since it could pass only by luck. Its
   outcome reports the result against the floor as well as against the entry. No test-board score is evidence.
6. **Weakest open group first.** Every row names the scored group it targets by its `vec.boards.GROUPS` id, its weight
   and its line in the headroom tables ("Where the entries stand against the published baselines"). A board's next
   try targets its weakest group (the most points left on the counting backtests; `docs/PLAN.md` section 3, item 9)
   whose lever is open; a try on another group says why the weaker ones are closed or have no reason that holds at
   test.
   - **Ties.** When a board's counting backtests name different weakest groups, those groups count as tied; across
     boards, gaps within the spread of a board's stand-ins count as tied. The row says which it targets first and why.
   - **Closed levers.** A lever is closed by the post-exposure block (`docs/DATA_POLICY.md`, 2026-11-04, (b), (c):
     cell numbers, sizes and scope on every board; T1's composition apart from births by the legal route; T3's
     composition), by a ruling of the captain (the heart's composition trend, dropped on 2026-11-01, `docs/PLAN.md`
     4.3, pipeline step 2), by a standing working answer (no new cross-modal input enters an upload until the
     organisers answer, A15) or by this ledger after its tries (the depth, blend, herding and coupling lines; the
     captain can reopen one). A lever fixed by a frozen or pre-registered rule moves only as that rule allows or by the
     captain's own ruling: on T3 a factor change that passes the adoption rule on the cross-knockout folds goes to him
     (`docs/PLAN.md` decision log 2026-11-01 (d)), and `experiments/t3_decomp/rules.yaml` is not edited. A closed
     lever stays closed whatever its headroom.
7. **Developmental reasons.** A developmental reason makes a change testable: it names the metrics the change should
   move and on which backtests. It never keeps a step the adoption rule rejects, and a gain off the predicted metrics
   is not kept (rule 2). An extra part of a change (a parameter, a step, a data source) needs its own held-out gain,
   measured with and without it, each run counted (rule 4); without that gain it is dropped.
8. **The tables follow the entries.** After every kept change or new entry, the script is rerun from committed
   results only (`python3 -I experiments/harness/headroom.py`), its output in `results/headroom/` is committed, and
   the headroom tables are copied from `headroom.txt`, naming the commit they were built at. An entry the script does
   not cover yet is added to its board list first. No official test score enters the tables.

## Counts per board

Before this ledger (from `results/registry.csv` and the decision files; floor and copy are reference rows, not tries):

| Board | Tries | Which |
|---|---|---|
| T1:test | 4 | b0, b0_scope, pseudobulk_shift (harness); m1_feng (labelled diagnostic, `experiments/t1_m1/results/diag_v1/`) |
| T2:embryo:test_interp | 6 | b0, pseudobulk_shift, psb_size, fgw (harness); fgw_shape0, fgw_shape1 (the s2 rule's run, `emb_shape_v1`) |
| T2:heart:test_extrap | 2 | b0, pseudobulk_shift (harness). The s2 hedge payoffs (heart scope) are a sensitivity analysis, not tries |
| T3:bcatenin | 1 | shift_transfer (the organisers' baseline); the T3 entries come from the pre-registered rules in `experiments/t3_decomp/` |

Since this ledger: see "Tries since 2026-11-03" and "Tries since 2026-11-04 (rules 5-8)" below. Running totals: T1 6, embryo 9, heart 3, T3 1 (H4 and E9 were pruned before their run; E10 was pruned by its gate and not counted; E11 counts).

## Where the entries stand against the published baselines

Built on 2026-11-04 by `experiments/harness/headroom.py` from the committed results at fcd8964 (stdlib only:
`python3 -I experiments/harness/headroom.py`; output in `results/headroom/headroom.txt` and `headroom.json`, with the
sha256 of its 19 inputs; the script and its output are added by the commit that adds this section). It reads
committed practice-task results only, never a target, a test-board file or a cell count, and it runs no candidate, so
it is not a try. Its self-checks: every stored board score recomputed from its skills to 1.4e-14 points, and every
mean used within 0.0005 of `results/registry.csv`.

Each cell gives a scored group's points against the published baseline, 100 x weight x (skill - 0.5), seed-paired,
then its points left to the split-half ceiling, 100 x weight x (1 - skill), on one backtest the adoption rule counts
for the entry (10 seeds). The baseline scores 50 on every seed of every backtest. * marks a stand-in, used where no
backtest can run the entry. "copy" is a fresh draw of the reference stage's real cells at the count rule, the
baseline as we would submit it: it is not the scorer's own working copy, so it can score below 50. The columns are
the `vec.boards.GROUPS` ids, which `headroom.txt` prints: T1 de_recovery, direction, distribution, covariance; T2
expression_change, cell_state, shape_growth, local_organisation (neighbourhoods); T3 response_identity,
response_direction, response_magnitude, mutant_state.

**T1:test** (counts 1 in Overall; baseline copy_last)

| s1 entry; backtest | Score (copy) | DE recovery (0.25) | Direction (0.25) | Distribution (0.30) | Covariance (0.20) |
|---|---|---|---|---|---|
| M1*: `m1_feng`, B0 plus Feng's direction at 1/5, a labelled diagnostic in-sample for the readout that chose Feng (`experiments/t1_m1/results/diag_v1/`, 401e3eb); E9.5 -> E10.5, chamber cut | 50.97 (48.06) | -0.50 / 13.00 | +3.01 / 9.49 | -0.08 / 15.08 | -1.45 / 11.45 |

- Weakest group: distribution (mmd_u), 15.08 points left, at the baseline. All of M1's gain over the floor is on
  direction; B0 alone (M1 without Feng's direction) scores 49.68.
- Open composition lever: births by the legal route (`docs/DATA_POLICY.md`, 2026-11-04 (c)); the rest of T1's
  composition is closed. Direction stays open under the block, but no backtest can judge a change to Feng's
  direction: Feng's E9.5 lies inside every T1 window.
- Shown, not counted: the whole-E10.5 truth, 45.21 (copy 47.98). Not scored by the harness: the s2 (M1 at 60%
  cardiomyocytes), whose evidence is the oracle-base grid in `experiments/s2_hedge/`.

**T2 boards** (each counts 0.5 in Overall; baseline copy_last). Groups weigh 0.25 each; the structure gate was 0.99
to 1.00 in every row, so the shape column is ungated in effect.

| Board; entry; backtest | Score (copy) | Expression change | Cell state | Shape and growth | Neighbourhoods |
|---|---|---|---|---|---|
| Heart s1, H1*: `b0`, H1 without the pre-registered size (`results/d2_v1/`, eae34eb); E10.5 from E9.5 | 50.43 (49.96) | +0.38 / 12.12 | -0.32 / 12.82 | -0.07 / 12.57 | +0.44 / 12.06 |
| same; E10.5 from <= E8.75 | 51.45 (50.55) | -0.46 / 12.96 | +0.44 / 12.06 | +0.19 / 12.31 | +1.28 / 11.22 |
| Embryo s1, fgw (`results/emb_fgw_v1/`, bd2e433); E7.5, t = 1/3 | 61.12 (50.58) | +3.54 / 8.96 | +4.21 / 8.29 | -0.13 / 12.63 | +3.49 / 9.01 |
| same; E7.25, t = 2/3 | 58.78 (53.27) | +4.81 / 7.69 | -2.22 / 14.72 | +3.64 / 8.86 | +2.55 / 9.95 |
| Embryo s2, fgw_shape0 (`results/emb_shape_v1/`, 90ddaaa); E7.5 | 62.95 (50.58) | +3.54 / 8.96 | +4.21 / 8.29 | +1.34 / 11.16 | +3.86 / 8.64 |
| same; E7.25 | 53.67 (53.27) | +4.81 / 7.69 | -2.22 / 14.72 | -0.49 / 12.99 | +1.56 / 10.94 |

- Weakest group, mean of the two backtests: heart, expression change (12.54 left, -0.04 against the baseline; every
  heart metric is within 1.3 points of it), a near-tie: cell state and shape leave 12.44 each, and the two backtests
  name different weakest groups (cell state on E10.5 from E9.5, 12.82; expression change from <= E8.75, 12.96), so
  they count as tied (rule 6); embryo s1, cell state (11.50 left, +1.00; the gap is E7.25's, the known feature in rule
  2), tied with shape (10.74), which E7.5 names (12.63); embryo s2, shape and growth (12.08 left), tied with cell state
  (11.50) in the same way.
- No open lever with a reason that holds at test on either board (`docs/PLAN.md` decision log 2026-11-04). Sizes,
  scope and cell numbers are blocked; the heart's composition trend was dropped by the captain (`docs/PLAN.md` 4.3,
  pipeline step 2); the depth (E7, H3), blend (E8), herding (H4, E9) and coupling (E11) lines are closed; the clock
  warp was pruned. An external heart direction or heart births would need cross-modal input, which no upload takes
  until the organisers answer (A15).
- The heart stand-in keeps the reference's coordinates, so its shape and scale skills are a copy's: they say nothing
  about the pre-registered size, which no backtest can score (A8). Not scored by the harness: the heart s2 (the half
  trim), whose evidence in `experiments/s2_hedge/` used an oracle core size.

**T3:bcatenin** (counts 1 in Overall; baseline wt_identity). No backtest runs the T3 entries: every fold is one they
were selected or fitted on. The stand-ins are rows of the pre-registered decomposition on the two Gata4 replicates,
scored against WT E8.75 by the T3 script's own scorer (`experiments/t3_decomp/results/gata4/gata4_records.csv`, 10
seeds); the cells are means over the two replicates, and the WT column is WT cells resampled (fac:none).

| Entry; stand-in | Score rep1 / rep2 (WT cells) | Response identity (0.30) | Direction (0.25) | Magnitude (0.25) | Mutant state (0.20) |
|---|---|---|---|---|---|
| s1*, B: `transfer:mab21l2_all`, the Mab21l2 factor on WT E8.75 cells, the one reading of a factor carried across knockouts (D11b, which failed its bar of 60 on rep1) | 59.40 / 62.15 (52.39 / 50.61) | -0.18 / 15.18 | +1.21 / 11.29 | +9.82 / 2.68 | -0.08 / 10.08 |
| s1*, A: `xrep:F+depth`, the Gata4 factor and depth borrowed from the other replicate (one knockout, so a fold the adoption rule never counts; no composition term) | 74.50 / 72.15 | +3.11 / 11.89 | +4.45 / 8.05 | +11.47 / 1.03 | +4.30 / 5.70 |
| s2*: `blind:G2_bcat_central_d50`, WT cells with s2's five states depleted (to 0.5, against s2's 0.61), no factor | 53.62 / 48.13 | -1.67 / 16.67 | +0.51 / 11.99 | +1.58 / 10.92 | +0.45 / 9.55 |

- Weakest group: response identity (de_score) on every stand-in, 15.18 points left across knockouts (11.89 on stand-in
  A).
- Open lever, through the captain only (T3's s1 is his ruling): a per-gene factor change that passes the adoption rule
  on the two cross-knockout folds. It would depart from the frozen T3 rules' outcome, so only his own ruling makes it,
  and `experiments/t3_decomp/rules.yaml` is not edited. Without his ruling both T3 uploads stay as the frozen rules
  build them; the composition is closed.
- Baseline rows on those folds (`results/baselines_v1/`; Gata4 -> Mab21l2, Mab21l2 -> Gata4 rep2): wt_identity 50 by
  construction; our copy of WT cells 50.74 and 48.59; the organisers' shift_transfer 47.05 and 42.18, against 69.75
  on the replicate fold.

**Across boards**, in Overall units (s1 or its stand-ins; `headroom.txt`, "Points left per scored group in Overall
units"): the largest gaps are T3 response identity 15.18 (11.89 with stand-in A), T1 distribution 15.08, T1 DE
recovery 13.00, T1 covariance 11.45, T3 direction 11.29 (8.05), T3 mutant state 10.08 (5.70) and T1 direction 9.49;
then the heart's four groups (expression change 6.27, cell state 6.22, shape 6.22, neighbourhoods 5.82) and embryo
cell state 5.75. The first two are tied (rule 6): T3's lies within its stand-ins' spread (11.89 to 15.18), and both
are a 0.30-weight metric at the floor. Points left in all ("Points left in Overall units per board"): T1 49.03, T3
39.23 (26.67 with stand-in A), heart 24.53, embryo 20.02.

Read with care:
- Points left are measured to the split-half ceiling, which no prediction reaches, and the skill map is 1/(1 + x), so
  equal points left are not equal reach. Near the floor (the M1 and H1 stand-ins, T3's s2) the ranking mostly
  restates the weights.
- These are backtests, not the test: T1 on E9.5 -> E10.5 (test E10.5 -> E12.5), whose stand-in draws E9.5's chamber
  mix against E10.5's chamber-cut truth; the heart 1 and 1.75 days ahead (test 2); the embryo at t = 1/3 and 2/3
  (test 1/2). Whether a backtest's weakest group is the test's is assumption A22.
- Earlier tries already aimed at weak groups (a check, not a change to their rows): T5 at B0's DE recovery on T1
  (de_score skill 0.450 on the chamber cut), H3 at the heart's cell state (-0.32 on the 1-day backtest), E7 and E8 at
  fgw's cell state (weak on E7.25 only), E11 at fgw's shape on E7.5 (its one group below the baseline there).

### What it means for the order of the next tries

1. **T1 births by the legal route.** The one open composition lever on T1's weakest groups (direction stays open
   under the block, but no backtest can judge a change to Feng's direction: item 5), and the one line that would
   predict states no released stage has (on the embryo, E8.0's real cells already bring in the states E7.5 lacks).
   Its design is a separate workflow (wf_77e506c8-014). Under rules 5-8 and `docs/DATA_POLICY.md` 2026-11-04 (c), its
   row names which states emerge (present in legal post-E14.5 chamber data, absent at E10.5) and their shares, both
   taken only from E10.5 and those data by a rule fixed before the run (the cardiomyocyte and blood shares stay as the
   s1 and s2 rules set them), the groups it should move (distribution and covariance first) and its held-out check,
   with the E9.5 chamber-mix confound of the chamber-cut backtest named in advance (as in the pruned exact top-up
   below). If that confound would decide the verdict, as it would have for the exact top-up, the row says so before
   the run; how the try is then judged is for the births design to settle, not this file.
2. **T3 response identity.** The largest gap in Overall units on the cross-knockout stand-in (15.18); on stand-in A,
   s1's own factor, it is 11.89, behind T1 distribution (15.08) and T1 DE recovery (13.00). So it ties with T1
   distribution (rule 6), and T1 goes first: its lever has a design under way and a date (the births decision by
   2026-11-10), and no T3 candidate exists yet. A factor change must pass the adoption rule on the two cross-knockout
   folds and then goes to the captain; the composition is closed. One reading the principle asks for, which changes
   no upload, is s1's own lead component (the Gata4 factor, without composition) on the Gata4 -> Mab21l2 fold: the
   first reading of it on a knockout it was not fitted on (A10). It needs a T3 harness candidate and counts as a try.
   Its row says in advance that this fold's WT reference is E9.5, not the E8.75 the test shares with Gata4, and what
   each outcome means: a pass is the first cross-knockout support for s1's factor; a loss under the adoption rule's
   test counts against its transfer (A10) and goes to the captain with T3's s1, which changes only by his ruling.
3. **Embryo: stop polishing.** It has the smallest gap (20.02 Overall points left) and the most tries (9), no group
   has an open lever with a reason that holds at test, and its one measured test risk, the flattened s1 geometry, is
   hedged by s2.
4. **Heart: no try.** Its entry is copy_last at the count rule plus a small count-space momentum, within 1.3 points of
   the baseline on every metric, and the pre-registered size, which no backtest can score. No lever on it has a
   reason that holds at test (above).
5. **T1 direction: no other try.** Feng's E9.5 lies inside every T1 window, and both T1 backtests that count for a
   same-scope change (the whole and the chamber-cut E10.5 truths) carry E9.5's chamber mix, so no check exists for a
   change to Feng's direction.

Once the heart and embryo s1 are uploaded, the frozen s2 rules build their s2 as s1, at s1's commit
(`docs/S2_RULES.md`), so from then on no try reaches those two boards unless the captain changes a rule.

## Tries since 2026-11-03

Order of decisions, fixed before the run: E7 against fgw (the s1 entry); then E8 against whichever of fgw and E7 is
the entry after that decision. Both on the two embryo backtests, 10 seeds, by `experiments/harness/run.py decide`.

| # | Board | Candidate | Reason | Gain should show on | A fit, not an effect, if | Registered at | Outcome |
|---|---|---|---|---|---|---|---|
| E7 | embryo | `fgw_depth` | Detection depth differs several-fold between the released embryo stages (median counts per cell E6.75 30, E7.25 35, E7.5 62, E8.0 212; `results/diag_emb_mix/`), and the scorer reads log(1 + 1e4 c / N), whose zeros and levels follow depth. fgw's mixture carries both brackets' depths (at test half its cells at 62 and half at 212). The target's depth is predicted log-linearly between the brackets' medians at the bracket fraction (`vec.t2_build.DEPTH_RULE`; at test about 115), and cells of a deeper bracket are binomially thinned to it, the counts a shallower acquisition of the same cell would give. Positions and draws are fgw's | variogram and mmd_u above all, then neighborhood_mmd; de_score and de_direction may move either way (depth moves pseudobulk); the shape terms identical to fgw's | it gains on one task only; or the predicted depth happens to sit near one target's depth while the rule misses the other (log-linear gives 64 against E7.5's 62 and 49 against E7.25's 35), so the gain should hold on both | f68836d | **Not kept** (`results/emb_depth_v1/`, 10 seeds): E7.5 -4.18, ahead in 0 of 10 (veto: mmd_u 0.671 -> 0.598, variogram 0.664 -> 0.602, neighborhood_mmd 0.640 -> 0.557); E7.25 +0.74, ahead in 10 of 10. Why: a mixture drawn (1 - t) : t already has the target's average log depth (35^(2/3) x 212^(1/3) = 64 against E7.5's 62), so thinning the deeper bracket pulled the population's average below the target; on E7.25 the mixture's average (49) sat above the truth (35), and thinning moved it toward it. What the scorer rewards is the population's average depth, which the mixture already gets right at the nominal fraction |
| E8 | embryo | `fgw_blend` | Development is continuous: an E7.75 cell is an intermediate state between its E7.5 ancestor's and its E8.0 descendant's, which a mixture of whole cells of the two brackets does not contain. Each drawn cell's counts become a binomial blend of itself and a partner from the other bracket chosen through fgw's coupling, splitting the expected counts (1 - t) : t at E7's predicted depth (a shallow side that cannot supply its share is capped at all its counts) | mmd_u, variogram and neighborhood_mmd beyond E7's; de_score and de_direction about unchanged; the shape terms identical to fgw's | it gains on one task only; or the gain comes from the blend's narrower within-state spread suiting one target (a blend averages two cells, so it carries about half a real cell's cell-to-cell variation) | f68836d | **Not kept**: E7.5 -2.23 and E7.25 -3.19, ahead in 0 of 10 on each (mmd_u 0.671 -> 0.638 and 0.497 -> 0.429; neighborhood_mmd falls on both). Blended cells are too narrow a population: the targets' cells are discrete real states, which a mixture of real cells matches better than averages of two |

| H3 | heart | `b0_depth` | (Restated after E7 and E8 and before this try's run: what the scorer rewards is the population's average log depth; E7 failed because it moved that average away from the target, and gained where it moved it toward the truth.) Counts per cell fall steadily through heart development in the released series (median counts per cell: E8.25_late 233, E8.5 166, E8.75 130, E9.5 75, E10.5 56; `results/diag_depth/`), and the scorer reads log(1 + 1e4 c / N), whose zeros and levels follow depth. B0 keeps E10.5's depth. The target's depth is extrapolated at B0's own damping: half the latest log rate per day, carried over the horizon (`vec.t2_build.DEPTH_EXTRAP_RULE`; at test 56 x (56/75)^(0.5 x 2) = 42); B0's cells are binomially thinned to it, never deepened. Leave-one-out on the depth series alone, each stage from the two before it: E8.75 134 (observed 130), E9.5 90 (75), E10.5 52 (56) | variogram, mmd_u, neighborhood_mmd; de_score and de_direction may move either way; the shape terms identical to B0's | it gains on one task only. Caveat written before the run: the survey that found the trend also showed both backtests' target depth (E10.5's 56), so the form of the rule (per-day rate, B0's damping 0.5, the size procedure's convention) was written knowing it, and E7's outcome was seen before this row was restated; the damping is B0's existing value, not fitted, and the leave-one-out above is the check that does not use a scored target | 71b18f9 | **Not kept** (`results/heart_depth_v1/`): the 1-day backtest -3.92, ahead in 0 of 10 (veto: mmd_u 0.466 -> 0.401, neighborhood_mmd 0.517 -> 0.458), although the extrapolated depth (52) was close to the truth's (56); the 1.75-day backtest +0.45, ahead in 9 of 10 (de_score and variogram up, mmd_u down). Across E7 and H3, thinning lowered mmd_u and neighborhood_mmd in three of four cases, whether or not it matched the target's depth: binomially thinned cells are not like a shallower acquisition's cells to these metrics. The depth line is closed (3 tries, none kept) |
| T5 | T1 | `b0_herd` | Sampling noise in a prediction's pseudobulk costs de_score: the metric ranks each gene's change of the mean log expression, the floor's change is exactly zero, and a fresh random draw of the reference's own cells scores -0.52 raw (skill 0.38) on E9.5 -> E10.5 (`results/t1_de_v1/`). Herding is a variance reduction with no biological content: within each category's quota, a drawn cell is swapped for an undrawn cell of the same category when that brings the sample's per-gene mean (log1p CP10k, all genes) closer to the category pool's mean (`vec.core.herd_rows`, random proposals, 20 per drawn cell, its own random stream so the draws it starts from are B0's); the same holds at test, since the scorer is the same | de_score above all, then de_direction; mmd_u and variogram about unchanged or slightly better (a more representative sample) | it gains on one task only; or it lowers mmd_u or variogram (herding toward the mean could favour typical cells and thin the tails) | 2bf32c4 | **Not kept** (`results/t1_herd_v1/`): -0.03 on E9.5 -> E10.5 and -0.01 on the chamber cut, ahead in 4 and 5 of 10 seeds. Why: most of the noise in a prediction's change comes from the scorer's own reference working copy (10% of the reference stage: 1,706 cells of E9.5 against B0's 5,118), which herding our draw cannot remove |
| T6 | T1 | `copy_herd` (diagnostic) | The same herding on a plain copy of the reference: measures how much of the copy's de_score handicap is sampling noise. Reported, not a candidate for an entry | de_score | - | 2bf32c4 | Small: +0.06, +0.21 and +0.17 on the three T1 backtests. The copy's handicap is the reference's noise, not ours |
| H4 | heart | `b0_herd` on the heart board | T5's reason on the T2 heart board, whose de_score also starts below the floor for a plain copy (skill 0.45-0.46 on both backtests): B0's draw of real cells is herded towards the pool's mean panel expression before B0's shift (`vec.t2_build.sample_cells(..., herd_expr=...)`, its own random stream); positions are the herded rows' own, jitter only for repeats, as before | de_score and de_direction; mmd_u, variogram and neighborhood_mmd about unchanged; the shape terms about unchanged (a different random-like subset of the same cloud) | it gains on one task only; or the cloud guard or a shape term moves (herding should not depend on position) | 47cd272 | **Pruned before its run** after T5: the reason does not hold, since the scorer's reference working copy carries most of the noise on T2 too (10% of the reference stage), and de metrics weigh half as much on T2 as on T1 |
| E9 | embryo | `fgw_herd` | The same on the embryo entry: each bracket's draw in fgw is herded towards that bracket's mean panel expression before positions are assigned (`interpolate_clouds(..., herd=True)`, its own random stream) | de_score and de_direction; the rest about unchanged | as H4 | 47cd272 | **Pruned before its run**, as H4 |
| E10 | embryo | `fgw_aligned` | fgw's coupling is an exact one-to-one pairing of its 2,000 anchors per stage: with Euclidean within-cloud distances the square-loss Gromov-Wasserstein term is concave on the couplings, so POT's conditional gradient ends on a vertex. Neither its features (panel expression at 30-212 counts per cell) nor the within-cloud distances can tell an embryo's left side from its right where expression is left-right symmetric (no lateralised marker before about E8.0), so partners can sit at random or mirrored across that direction, and McCann moves then pull both brackets towards its mid-plane: the moved cloud flattens there (wholly at t = 1/2 with mirrored partners), and so does the s2's mapped E8.0 half. The proper rotation is mostly unaffected (the other two axes and det +1 fix it). Change: keep FGW's rotation; take the moves (still the mean displacement of the 8 nearest anchors) from an exact transport of the same anchors on (1 - a) M / max M + a \|x - y R\|^2 / max, in the before stage's frame, a = INTERP_ALPHA = 0.5 reused (`vec.core.aligned_coupling`; `interpolate_clouds(..., move_coupling="aligned")`); in a common frame position tells the side. Synthetic left-right symmetric pair (tests/test_t2_build.py, AI-assistant scratch): FGW's partners have whitened cross-correlation about 0 along the symmetric axis, the aligned ones 0.4-0.7 (300-500 anchors), and the moved cloud's thinnest SD returns to the brackets'. **Gate, target-free, fixed here before it runs** (decided by the AI assistant under the captain's delegation of 2026-11-02, as the ideas workflow proposed): `experiments/harness/diag_emb_coupling.py` on the three bracket pairs at build seeds 1,000,000-1,000,009. A draw has the defect if fgw's smallest whitened anchor-partner cross-correlation is below 0.3 or its moved cloud's thinnest SD over RMS radius is below 0.9 x the thinner bracket's; the fix works if the aligned partners' cross-correlation is at least 0.6 and their merged-label agreement at most 0.05 below fgw's. The try runs only if draws with both reach 6 of 10 on each backtest pair and 3 of 10 (or the s1 seed's draw) on the test pair. A defect at test that the backtests cannot show: pruned, reported to the captain as a test-only risk. The diagnostic reads no target and no score, so it is not a try | E7.5: d2_shape up in at least 8 of 10 seeds, occupancy_dice flat or up; E7.25: d2 and Dice up in the seeds whose draws show the defect, about unchanged in the rest; de_score, de_direction, mmd_u and variogram identical to fgw's (the same cells); scale_log_ratio identical up to the scorer's rounding; neighborhood_mmd within +-0.01 skill on both | the gain sits in Dice alone, or only in E7.25's seeds 3, 4, 5, 8 and 9 (that block's anisotropy: E7.25's truth is thinner along PC3 than either bracket); or mean_anchor_shift falls well below fgw's while d2 moves toward fgw_shape0's (shorter moves across the sectioning-plane change, not the mechanism); or the gain is confined to draws where the rotation is a half-turn off (rotation and label agreement reported per seed). Caveats written before the run: the idea came after E7.25's seed split was seen; a = 0.5 is reused in a new role and no other value is tried; the backtests sit at t = 1/3 and 2/3, so a pass says nothing about t = 1/2 beyond the diagnostic | 5baba46 | **Pruned by its gate, not run, not counted** (`results/diag_emb_coupling/`, at 5baba46, clean). The defect showed in 10 of 10 draws on the test pair and on E7.5, and in 5 of 10 on E7.25 (seeds 3, 4, 5, 8 and 9, the seeds of fgw's earlier bimodal split there: mirrored partners), short of the 6 the gate needs on each backtest pair; the fix worked in every draw of all three pairs. Re-registered as E11 (below) |
| E11 | embryo | `fgw_aligned` (as E10) | E10's change, re-registered after its gate measured the defect on the real pairs (target-free, `results/diag_emb_coupling/`). **Test pair, every draw (10 of 10, the s1 seed included):** FGW's partners are mirrored across the symmetric direction (smallest whitened cross-correlation -0.26 to -0.50) and the moved cloud's thinnest axis is 0.21-0.25 of its RMS radius, against 0.36 for E8.0 and 0.53 for E7.5 (0.58-0.70 of the thinner bracket's): a cloud flatter than either bracket, which no sectioning of E7.75 gives; only 44-47% of each sampled cell's 15 nearest neighbours stay among them. **E7.5 backtest, every draw:** partners at random along that direction (cross-correlation about 0), thinnest axis 0.70-0.76 of the thinner bracket's. **E7.25, seeds 3, 4, 5, 8, 9:** mirrored (-0.55 to -0.66; 0.77-0.80). The aligned moves fix every draw: cross-correlation 0.71-0.87, thinnest axis 0.99-1.29 of the thinner bracket's (at test 0.43, between the brackets, as displacement interpolation should be), 65-75% of neighbourhoods kept; merged-label agreement 0.01-0.03 lower. Spending a try after the gate said no: decided by the AI assistant under the captain's delegation of 2026-11-02, because the gate was a rule for saving tries the backtests cannot judge, while E7.5 shows the defect in every draw and the test pair carries it in every draw; the adoption rule and rule 2 are unchanged. The captain can overrule | E7.5: d2_shape up in at least 8 of 10 seeds, occupancy_dice flat or up, neighborhood_mmd up or within +-0.01 skill; E7.25: d2 up and Dice down in seeds 3, 4, 5, 8 and 9 (E7.25's truth is thinner along its third axis than either bracket, so the flattened mode wins Dice there), small changes in the rest; de_score, de_direction, mmd_u and variogram identical to fgw's; scale_log_ratio identical up to the scorer's rounding | as E10. Written before the run: if E7.5 passes and E7.25 does not (without a veto), the change is not kept and the test-only risk goes to the captain with both results; a pass carried by Dice, or by E7.25's mirrored seeds alone, is not kept | b705279 | **Not kept** (`results/emb_aligned_v1/`, 10 seeds, b705279, clean): E7.5 +2.13, ahead in 10 of 10, all on shape as predicted (d2 skill +0.129, Dice +0.157; neighborhood_mmd -0.010, a little past the predicted band in 5 seeds; expression metrics identical); E7.25 -1.59, ahead in 1 of 10 (veto). Why E7.25 loses: the position term shortens the moves (mean anchor shift 0.48 against fgw's 0.62 where FGW was right), so d2 falls in the five unmirrored seeds (-0.13 to -0.19), and in the five mirrored ones the restored thickness costs Dice (-0.07 to -0.32) because E7.25's truth is itself thinner along its third axis than either bracket; neighborhood_mmd falls in every seed (-0.007 to -0.035). Correction to this row: "flatter than either bracket, which no sectioning of E7.75 gives" was too strong; E7.25 is such a stage. Pair check, the s2 rule's better-of-two criterion with s2 = fgw_shape0: s1 = fgw 60.87, s1 = fgw_aligned 60.24, so the s2 already covers fgw's flattened geometry on E7.5-like targets. The s1 stays fgw; at test its moved cloud is flattened in every draw (a known risk, hedged by the s2, which takes E7.5's geometry) |

## Tries since 2026-11-04 (rules 5-8)

Each row is committed before its run, as rule 1 says, and fills every column (rules 5-7). None yet.

| # | Board | Candidate | Departs from (baseline; entry) | Targets (GROUPS id, weight; table line) | Reason (the hypothesis) | Predicted change beyond the baseline (facet; metrics; expected gain) | Held-out evidence (counting backtests or folds, none tuned on) | A fit, not an effect, if | Registered at | Outcome (against the floor and the entry) |
|---|---|---|---|---|---|---|---|---|---|---|

## Pruned before running

- **Partial or unbalanced FGW coupling** (the plan's "mild unbalancing"). Its reason would be tissue that one bracket
  sampled and the other did not. But the embryo shape loss on the E7.5 backtest comes from the sectioning plane: z is
  the long axis at E6.75 and E7.25 and the thin axis at E8.0 (`experiments/coords/REPORT.md`), and the stack length is
  set by acquisition. Leaving mass unmatched does not address that; the s2 shape variant already hedges it.
- **A developmental-clock warp of the bracket fraction t.** On the backtests the targets sit off their nominal
  fractions in opposite directions (E7.25 near E6.75; E7.5 about midway between E7.25 and E8.0 by mmd_u), and at test
  the nominal fraction is 1/2. Nothing measured on released stages says which way E7.75 would sit.
- **Sweeps of FGW alpha, k or the anchor count.** No physical reason for another value than the plan's.

From the ideas workflow of 2026-11-03 (seven read-only proposers, three critics per idea, a ranking; results not
committed, their synthetic toys in the AI assistant's scratch). E10 above was its one survivor worth a try.

- **fgw_tissue (embryo): couple the brackets on 15-NN tissue profiles.** Tissue profiles make partners consistent but
  cannot set their left-right sign: expression is symmetric and FGW is mirror-blind. In a synthetic test it mirrored in
  about half the draws and lost on average with a round cross-section.
- **fgw_pairs (embryo): one coupled partner per cell.** The FGW coupling is already one-to-one, so this is exactly
  k = 1, a value of the pruned k sweep; in a synthetic test it raised raw neighborhood_mmd by 18-29%.
- **b0_smooth (heart): B0's top-up from kNN-smoothed expected counts.** B0's gain over the copy on the 1.75-day
  backtest lies on the same variogram-vs-variance line as pure thinning, that is, it comes from the own-count
  top-up's detection loss, which this change removes; heart depth falls at every released stage, so the test target
  is likely shallower still. About 3-5% to pass, the largest build of the lot.
- **The same exact top-up on T1.** The own-count top-up realises about half the intended rise of a sparse gene, but
  the expected test gain is 0.1-0.4 points, and both counted T1 backtests would judge it backwards through their
  composition mismatch (B0 draws E9.5's chamber mix, cardiomyocytes about 64% of cells, against about 19% in the
  chamber-cut truth), so either verdict would be that confound.
- **T3: a background- or spillover-aware factor.** It needs a new T3 harness candidate, a threshold for "background
  level" and an amendment of rules note 5 (hash-pinned); its sign cannot be checked on released data without tuning
  on a counted fold. About +0.3 expected, with an SD near 1.5.
- **A mirror-only correction of FGW's partners** (keep FGW's moves; put each partner on its anchor's side of the
  symmetric direction), considered after E11. It would avoid E11's shortened moves, but its effect on E7.25 would sit
  in the same five mirrored draws, where restored thickness already cost Dice, so it faces E11's veto for the same
  reason, and on E7.5 its gain would be E11's. Not tried; the coupling line is closed for the s1 unless the captain
  rules otherwise.
- **T3: an absolute knockout-depth rule; the factor measured against the target's WT.** Both would gain on the two
  cross-knockout folds and change nothing at test, where the knockouts share the WT reference: lucky-fit traps.

## Diagnostics behind these choices

2026-11-03, the embryo coupling (`experiments/harness/results/diag_emb_coupling/`, at 5baba46, clean; reads no target
and no score; every pair on T2:embryo:val_interp's panel, since the test board is not in boards/ yet). Per pair, over
build seeds 1,000,000-1,000,009: the smallest whitened anchor-partner cross-correlation of the entry's moves (fgw) and
the aligned ones, and the moved cloud's thinnest principal SD over its RMS radius at the pair's fraction (brackets in
brackets):

| Pair (t) | fgw cross-correlation | aligned | fgw thinnest axis | aligned | brackets' thinnest | 15-NN kept, fgw / aligned |
|---|---|---|---|---|---|---|
| E6.75, E7.5 -> E7.25 (2/3) | 0.59-0.66 in seeds 0, 1, 2, 6, 7; -0.55 to -0.66 in 3, 4, 5, 8, 9 | 0.84-0.87 | 0.45-0.46; 0.29-0.30 | 0.48-0.49 | 0.38, 0.53 | 0.65-0.67; 0.57-0.58 / 0.74-0.75 |
| E7.25, E8.0 -> E7.5 (1/3) | -0.06 to 0.09 | 0.71-0.73 | 0.25-0.27 | 0.35-0.36 | 0.38, 0.36 | 0.59-0.61 / 0.70-0.72 |
| E7.5, E8.0 -> E7.75 (1/2, test) | -0.26 to -0.50 | 0.73-0.75 | 0.21-0.25 | 0.43-0.44 | 0.53, 0.36 | 0.44-0.47 / 0.65-0.66 |

The s1 seed's test draw (1,000,001): cross-correlation -0.31, thinnest axis 0.67 of the thinner bracket's (aligned:
0.73 and 1.20). At t_pos = 0 (the s2's geometry) the test pool's thinnest axis is 0.47-0.51 with fgw's moves and
0.51-0.53 with the aligned ones (E7.5: 0.53); the mapped E8.0 half alone 0.33-0.43 against 0.48-0.51.

2026-11-03, the embryo stages (`experiments/harness/results/diag_emb_mix/`; DGX, released stages only, the bake
protocol's truth half). mmd_u and variogram of plain copies and of mixtures against each embryo backtest target:

| Target | Earlier bracket copied | Later bracket copied | Mixture at the nominal t | Mixture at 1/2 |
|---|---|---|---|---|
| E7.25 (E6.75, E7.5; t = 2/3) | 0.0206, 0.0088 | 0.0541, 0.0547 | 0.0264, 0.0271 | 0.0152, 0.0139 |
| E7.5 (E7.25, E8.0; t = 1/3) | 0.0858, 0.0529 | 0.0829, 0.1020 | 0.0433, 0.0262 | 0.0380, 0.0304 |

Lower is better (each cell: mmd_u, variogram). E7.25 is far closer to E6.75 than its day says, which is why fgw's
variogram falls on that backtest (it draws 2/3 of its cells from E7.5); E7.5 sits about midway by mmd_u. Mixtures
beat either bracket alone on both targets.

Detection depth, median counts per cell (duplicates dropped; N = the library each released value was normalised by):
E6.75 30, E7.25 35, E7.5 62, E8.0 212; E10.5 (heart) 56. Part of E7.25's closeness to E6.75 is depth (35 against 30
and 62), and the test's brackets are the furthest apart in depth of any (62 and 212). The table above was first run
from a scratch copy of this script on 2026-11-03 at 1093aae; `results/diag_emb_mix/` holds the rerun at the commit
that registered E7 and E8.
