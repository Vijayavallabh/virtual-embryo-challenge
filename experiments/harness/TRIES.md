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
     there and says nothing about E7.75, where the bracket fraction is 1/2. A second (added 2026-11-04 from the
     rebuild design, wf_237f8d8f-e7c): by mmd_u both embryo backtests prefer equal mixing to their nominal
     fraction (0.0152 against 0.0264 on E7.25, 0.0380 against 0.0433 on E7.5; the 2026-11-03 table below), so a
     change that pulls either backtest's mixture toward 1/2 gains on both and says nothing at test, where the
     nominal fraction is already 1/2.
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
   - **Closed levers.** A lever is closed by the post-exposure block (`docs/DATA_POLICY.md`, 2026-11-04, (b) to (d):
     cell numbers, sizes and scope on every board; T3's composition; T1's composition except by the legal route of
     (d)), by a ruling of the captain (the heart's composition trend, dropped on 2026-11-01, `docs/PLAN.md`
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

Since this ledger: see "Tries since 2026-11-03" and "Tries since 2026-11-04 (rules 5-8)" below. Running totals: T1 7, embryo 9, heart 3, T3 3 (H4 and E9 were pruned before their run; E10 was pruned by its gate and not counted; E11 counts; K2a and K2b, registered at 945de7d, ran at 0409d18 and count; T1's seventh is a campaign designer's unregistered comparison against E10.5 of 2026-11-05, counted under rule 4: "Pruned before running", the campaign design of 2026-11-05). N1 and N2 read no target and no score, so they are not tries, and R1 is a ruling request with no run.

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
- Open composition lever: T1's shares within the chamber scope, births included, by the legal route
  (`docs/DATA_POLICY.md`, 2026-11-04 (c), (d)); the births design was pruned on 2026-11-04 ("Pruned before running").
  Direction stays open under the block, but no backtest can judge a change to Feng's direction: Feng's E9.5 lies
  inside every T1 window.
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

1. **T1 composition by the legal route.** Births by the legal route were designed and pruned on 2026-11-04
   (wf_77e506c8-014; "Pruned before running"). `docs/DATA_POLICY.md` 2026-11-04 (d) then reopened T1's shares within
   the chamber scope by the legal route, and (f) says how a share rule is built (quotas at the frozen count, never
   importance weights; T1 depth closed). The rebuild design of 2026-11-04 (wf_237f8d8f-e7c) registered no
   composition try: no legal share rule has a counted check, since both counted T1 tasks draw E9.5's chamber mix
   (cardiomyocytes about 64% of cells, against about 19% in the chamber-cut truth: the pruned exact top-up below),
   and the committed oracle-base grid (`experiments/s2_hedge/`) puts the current s1/s2 pair within about 1.5 points
   of the best pair under every prior it states. A share change can enter only by the captain's ruling after he
   writes a prior (item 6). What remains on T1 is R1, a ruling request: cells without a Feng group keep their E10.5
   expression (row below). A later composition row still follows rules 5-8 and names that confound in advance.
2. **T3 response identity.** The largest gap in Overall units on the cross-knockout stand-in (15.18); on stand-in A,
   s1's own factor, it is 11.89, behind T1 distribution (15.08) and T1 DE recovery (13.00). So it ties with T1
   distribution (rule 6), and T1 went first; the rebuild design found no T1 try (item 1) and no T3 candidate with a
   reason that holds at test ("Pruned before running"). A factor change must pass the adoption rule on the two
   cross-knockout folds and then goes to the captain; the composition is closed. The one reading the principle asks
   for, s1's own lead component (the Gata4 factor, without composition) on the Gata4 -> Mab21l2 fold, the first
   reading of it on a knockout it was not fitted on (A10), is registered as K2a with a permuted-factor control K2b
   (rows below; decided by the AI assistant under the delegation of 2026-11-02, the captain can overrule). It changes
   no upload: T3's s1 changes only by the captain's ruling, whatever the outcome.
3. **Embryo: stop polishing.** It has the smallest gap (20.02 Overall points left) and the most tries (9), no group
   has an open lever with a reason that holds at test, and its one measured test risk, the flattened s1 geometry, is
   hedged by s2.
4. **Heart: no try.** Its entry is copy_last at the count rule plus a small count-space momentum, within 1.3 points of
   the baseline on every metric, and the pre-registered size, which no backtest can score. No lever on it has a
   reason that holds at test (above).
5. **T1 direction: no other try.** Feng's E9.5 lies inside every T1 window, and both T1 backtests that count for a
   same-scope change (the whole and the chamber-cut E10.5 truths) carry E9.5's chamber mix, so no check exists for a
   change to Feng's direction.
6. **After the first wave.** No further try is registered without a new reason that holds at test, written in its
   row before the run. The rebuild design's second-wave items run only at the captain's request: H5, a
   placebo-controlled reading of B0's momentum on both heart backtests (information only); D3, pricing T1's
   cardiomyocyte tail, which changes a frozen s2 rule and runs only after the captain writes a prior (under the
   stated priors it loses, 68.28 against 67.89 with a uniform prior; break-even needs about 25% of the prior on a
   cardiomyocyte share of at least 0.8); T3's s2 moved to the barred branch's s1 (a frozen rule; it would override
   D11b's failed bar after the fact; not recommended); a zone-fixed Feng cardiomyocyte direction (about 0 +- 0.3;
   not recommended). Scenarios kept as risks and not built (closed levers): T1 depth and the oesophagus-junction
   cells; the heart's core-focused s2 (scope; a frozen rule); stack-limited embryo geometry (z extent; a frozen
   rule); T3's composition magnitude (a frozen rule).

The captain's upload gate (2026-11-05; `docs/PLAN.md` section 8): no first upload on any board until our estimated
leaderboard total is above 209 (about 184 today) or his deadline exception opens (the s1 methods freeze at 2026-11-27
00:00 UTC; from then the first round may go in, each upload only on his go-ahead at the time), so every entry stays
open to tries until then. The gate is about uploads only: it is never evidence for a try, and a try's row still names
its own expected gain on its counting backtests.

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

Each row is committed before its run, as rule 1 says, and fills every column (rules 5-7). The first wave comes
from the rebuild design of 2026-11-04 (wf_237f8d8f-e7c: its synthesis's rows, with the compliance audit's fix to
N2 and K2b's permutation seed fixed by the AI assistant). N1 and N2 are diagnostics that read no target and no
score, so they are not tries (as E10's gate); R1 is a ruling request with no run; K2a and K2b are counted tries.
A row's registration commit is named with its outcome.

| # | Board | Candidate | Departs from (baseline; entry) | Targets (GROUPS id, weight; table line) | Reason (the hypothesis) | Predicted change beyond the baseline (facet; metrics; expected gain) | Held-out evidence (counting backtests or folds, none tuned on) | A fit, not an effect, if | Registered at | Outcome (against the floor and the entry) |
|---|---|---|---|---|---|---|---|---|---|---|
| N1 (diagnostic, not a try) | embryo | `diag_emb_coupling.py` rerun on the test pair with the test panel (497 genes) | copy_last; s1 fgw and s2 fgw_shape0, unchanged | shape_growth (0.25): s1 on E7.5 -0.13 / 12.63; s2 12.08 left; a check only | Every embryo backtest and diagnostic so far used the 498-gene validation panel; the test panel is new (D1, ae6d6b3). The recorded test-only risk the s2 hedges (FGW partners mirrored in every test draw, the moved cloud's thinnest axis 0.21-0.25 of its RMS radius; the s2 pool 0.47-0.51) was measured on the old panel | None to any entry. Expected: the recorded test-pair values hold within 0.03. Expected gain 0 | Reads only E7.5 (unlocked, A20) and E8.0, the test's inputs, at build seeds 1,000,000-1,000,009 and the s1 seed 1,000,001; prints no cell count; no target and no score, so not a try (as E10's gate) | Not applicable. Fixed now: if the s1 seed's thinnest axis or the s2 pool's leaves its recorded range by more than 0.05, the captain gets a test-only risk note; nothing is rebuilt on this result | 945de7d | **Hold** (`results/diag_emb_coupling_testpanel/`, at 0409d18, clean; target-free, not a try). On the test panel (497 genes) every value recorded on the validation panel holds at two decimals, departure 0.00 on all five: fgw's smallest cross-correlation -0.28 to -0.50 (recorded -0.26 to -0.50), the moved cloud's thinnest axis 0.21-0.25 of its RMS radius (0.21-0.25), the s1 seed's draw -0.31 and 0.67 of the thinner bracket's (-0.31 and 0.67), the s2 pool 0.47-0.51 (0.47-0.51; the s1 seed's draw alone 0.48). Reported, not judged: partners mirrored in 10 of 10 draws; 15-NN kept 0.44-0.47 (0.44-0.47); thinnest axis 0.58-0.69 of the thinner bracket's (E11: 0.58-0.70). No risk note: the test-only risk the s2 hedges is as recorded, and nothing is rebuilt |
| N2 (diagnostic, not a try) | T1 | size of M1's momentum on cells without a Feng group | copy_last; s1 M1 (s2 as s1) | distribution (mmd_u, 0.30): -0.08 / 15.08; de_recovery (0.25): -0.50 / 13.00; a check only | Gives ruling R1 the size of what it removes: M1 moves Blood, Macrophage and Early VIC (about a quarter of s1's cells at E10.5's chamber mix) by half the E9.5 -> E10.5 per-category change | None to any entry. Reports per category the genes whose factor moves by more than 10% and 25%, the mean absolute log factor, and the share of s1's and s2's cells touched, computed from E10.5's chamber mix and the s2 rule (not from built files) and printed to two decimals; it prints no cell count. Expected gain 0 | Reads released E9.5 and the unlocked E10.5 as the test build's inputs, as build_t1 does; no target, no score; not a try | Not applicable; its numbers go to the captain with R1 and decide nothing else | 945de7d | **Reported** (`experiments/t1_m1/results/n2_momentum_size/`, at 0409d18, clean; target-free, not a try). The momentum moves blood (Blood and Macrophage, the whole category) and endocardium's Early VIC cells (0.20 of that category at E10.5; EndoMT Endocardium, Endocardium-1 and Endocardium-2 follow Feng's Endocardium direction), by log f = 0.5 x [log(m_E10.5 + 0.5) - log(m_E9.5 + 0.5)] of the category's mean CP10k, no clip, the 23 held-flat genes at f = 1. Blood: of the 21,020 genes expressed in the moved E10.5 cells, the factor moves 2,375 by more than 10% (11.3%; 796 up, 1,579 down) and 218 by more than 25% (1.0%); mean absolute log factor 0.038; those genes hold 0.56 and 0.23 of the cells' CP10k. Early VIC: 1,370 of 20,622 (6.6%; 888 up, 482 down) and 186 (0.9%; 181 up); 0.029; 0.26 and 0.07. With apply_factors' keep-depth constant s (median 0.91 in blood, 0.96 in Early VIC) the applied change f s moves 7,443 blood genes by more than 10% (7,265 of them down) and 494 by more than 25%; 1,785 and 144 in Early VIC. Share of cells the momentum moves: s1 0.31 (blood 0.24, Early VIC 0.07), s2 0.28 (0.24 and 0.04), against the row's "about a quarter" of s1. Decides nothing; the numbers go to the captain with R1 |
| R1 (ruling, no run) | T1 | `m1_hold_unbracketed`: cells without a Feng group keep their E10.5 counts | copy_last; s1 M1; s2 follows as s1, and a changed s1 sends its s2 back to the captain (DATA_POLICY 2026-11-04 (d)) | distribution (mmd_u, 0.30): -0.08 / 15.08, the weakest open group; then de_recovery (0.25): -0.50 / 13.00 | M1 extrapolates half the E9.5 -> E10.5 per-category change for Blood, Macrophage and Early VIC. That step spans T1's sampling change (an embryo portion at E9.5, a dissected heart at E10.5). On the one transition measured, the same momentum had no sign skill on Blood: sign agreement 0.50 full and 0.55 cross-fitted, Spearman -0.06 and -0.05, de_direction_full -0.17 and -0.02 (pre-registered readout, experiments/t1_direction/results/E10.5/per_type.csv). Early VIC's category change contains the appearance of EndoMT and Early VIC states, a composition change carried as expression. A move with no measured sign skill adds distance in expectation | Expression only: these cells are predicted at their E10.5 expression, as the baseline predicts them; no state, share, cell number, size or scope changes. Metrics: mmd_u (the blood cluster), then de_score and de_direction (fewer blood genes among the top-ranked changes); the variogram may move either way. Expected at test +0.1 to +0.3 board points (range -0.3 to +0.5), in s1 and s2 | None that counts: both counted T1 tasks draw E9.5's chamber mix (blood 2.6% of chamber cells), so the expected change there is below 0.1 and rule 5 prunes it as a try. It enters only by the captain's ruling, as Feng's direction did (decision log 2026-11-01 (f)). Its evidence is the readout rows above (in-sample for that readout) and N2's sizes | Not applicable (no run). Fixed now: no test score confirms or reverses it; if the captain declines, M1 stays as built | this row's commit (named with its outcome) | |
| K2a | T3 | `t3_factor_xko`: s1's lead component (Gata4's pooled non-Mesp1 factor, exponent 1) on deduplicated WT E9.5 cells at WT composition, no depth, no composition | wt_identity; s1 (allowed branch: factor, Gata4 per-label depth, central composition at 0.38x), unchanged whatever the outcome | response_identity (de_score, 0.30): 15.18 left across knockouts, 11.89 on stand-in A; then response_direction (0.25): 11.29 (8.05) | A10: s1's factor, measured in compartments a Mesp1-Cre deletion cannot touch, is processing shared across knockouts, so it carries to a knockout it was not fitted on; no unseen knockout has tested it. Built with the T3 script's own primitives (gene_factor, eps 1e-5, clip 0.25-4; apply_gene_factor with the kNN-smoothed top-up profile, rules note 5). Depth is off because the knockouts' depth ratios flip between folds (0.52-0.63x against 1.27x), the pruned lucky-fit trap | Expression: the mutant's per-gene processing response predicted from another knockout's factor. Metrics: de_score and de_direction above K2b (gene identity); severity may sit at the r2 gate in both arms. Expected: both arms near or below 50 (the reverse fold's analogue, Mab21l2's factor against its own WT, scored 52.04 and 46.24). The reading is K2a minus K2b. Expected Overall gain 0 (no upload changes) | The cross-knockout fold T3/Gata4->Mab21l2 (configs/pseudo_tasks.yaml: Mab21l2 KO against WT E9.5), 10 seeds, never used by the T3 rules; reported against wt_identity (50), our WT copy (50.74) and shift_transfer (47.05). Written now: this fold's WT reference is E9.5, not the E8.75 the test shares with Gata4, so the factor crosses a WT run whose large per-gene file difference the test does not have. The fold is harsher than the test, and a fail is weak evidence against A10. Outcomes: if K2a beats K2b by the adoption rule's test (paired mean above 0.3, ahead in at least 8 of 10 seeds), that is the first cross-knockout support for s1's factor; otherwise there is no support, and it goes to the captain with T3's s1, which changes only by his ruling. No T3 file is rebuilt either way | The margin sits on severity or mutant state alone; or more than half the de_direction margin comes from 10 genes; or the r2 gate differs between the arms in more than 2 seeds, so the cliff decides rather than gene identity | 945de7d | **Not kept: no cross-knockout support** (`results/t3_k2_v1/`, 10 seeds, at 0409d18, clean; run by the AI assistant in the first-wave workflow wf_b4b013b3-a46). K2a 50.08 +- 0.47 against K2b 49.90 +- 1.80: paired +0.18, ahead in 9 of 10 seeds, short of the adoption rule's mean above 0.3, so this fold gives s1's factor no support (A10; the fold is harsher than the test, so this is weak evidence against it). Against the baselines: wt_identity (50) +0.08, ahead in 5 of 10; our WT copy (50.74) -0.65, ahead in 4 of 10 (a veto); shift_transfer (47.05) +3.04, ahead in 10 of 10 (copy and shift_transfer reproduce `baselines_v1` on every seed). Against K2b by group: response identity +0.27 and direction +0.38 (de_score 0.064 against 0.032, de_direction 0.007 against -0.054, the predicted side), magnitude -0.53, mutant state +0.06; against wt_identity: identity +0.54, direction +0.04, mutant state -0.50 (mmd_u and variogram skills 0.49 and 0.45). The fit checks: the margin sits on severity or mutant state alone, not met (severity works against K2a); more than half the de_direction margin from 10 genes, met (their share is 1.59 of a +0.061 margin; the other genes net negative); the r2 gate differing in more than 2 seeds, not met (1 seed: seed 3, where K2b's severity is scored; severity sits at the r2 gate in both arms on the other 9, as predicted, and the parts reproduce the panel on every seed). It goes to the captain with T3's s1, which changes only by his ruling; no T3 file is rebuilt and no upload changes. Added after the workflow's review, reported and not judged: (1) K2a's factor equals s1's factor file exactly on real data (max abs diff 0.0 over its 500 genes; the file's sha256 is the one `gata4_decisions.json` records; `results/t3_k2_v1/factor_check.txt`). (2) The miss rests on seed 3, where K2b's severity is scored and K2a's is r2-gated (K2a minus K2b -4.36). On the other nine seeds K2a leads K2b by +0.69, ahead in 9 of 9, so the factor's gene identity carries some signal and the verdict turns on condition 3's cliff in one seed. Against our WT copy K2a stays behind on those nine seeds too (-0.11, ahead in 4 of 9), so on this fold the factor does not beat a WT copy whichever way seed 3 is read. (3) Condition 2 is a weak test at margins this small: the per-gene de_direction parts are +0.515 and -0.454 around the +0.061 net, so almost any small margin meets it (a design note for later rows). No seed is excluded and the verdict stands |
| K2b | T3 | `t3_factor_xko_perm` (control): K2a's factor permuted across genes within deciles of WT E8.75 mean expression, one permutation, seed 20261104, fixed in this row | as K2a | as K2a | Separates the factor's gene identity from its magnitude and offset, which a comparison with wt_identity alone cannot | None by design. Expected: about K2a's level on severity, and below it on de_score and de_direction if A10's gene-identity part holds | as K2a; counts as a try (rule 4), so T3 goes from 1 to 3 | as K2a | 945de7d | Run with K2a and counted (`results/t3_k2_v1/`, at 0409d18, clean; T3 now 3): 49.90 +- 1.80; against wt_identity -0.10, ahead in 1 of 10; our WT copy -0.84, ahead in 1 of 10 (a veto); shift_transfer +2.85, ahead in 10 of 10. Below K2a on de_score and de_direction, as predicted if A10's gene-identity part holds; above it on severity (skill 0.521 against 0.500), all from seed 3, where its severity is scored rather than gated and it scores 54.95 (48.87-49.79 on the other nine). The K2a-K2b reading is in K2a's row |

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
- **T1 births by the legal route (`b0_births`; never registered).** Designed 2026-11-04 by the AI assistant's design
  workflow (a data scout, a code scout, a designer and three critics); design only, no T1 score of any kind. Pruned by
  the AI assistant under the captain's delegation of 2026-11-02; the captain can overrule. No developmental knowledge
  after E9.5 up to and including E13.5 entered any class, share or state: all came from Feng 2022's E14.5 cells and the
  released E8.5 and E9.5 stages. The strongest rule turned some of the latest stage's epicardium-category cells into
  fibroblast-like cells, shifted by Feng's within-E14.5 fibroblast-like-over-epicardial ratio, at t x (Feng's E14.5
  share - the latest stage's share) of heart cells, t = (target day - latest day) / (14.5 - latest day); immune was never
  built. No result of its run could put births into T1. (1) No rule written before any run makes a birth on the counted
  tasks (`docs/DATA_POLICY.md`, 2026-11-04 (c)). The scout's rule, fixed in its RULES block before its run (argmax and
  marker score both under 0.2 x the type's share of all Feng E14.5 cells; uncommitted, script saved 08:13 UTC, outputs
  08:14), finds only immune absent at E9.5, and immune cannot be built (4), so `b0_births` would equal `b0` on both.
  Fibroblast-like comes in only through a confident-call test (p >= 0.9) and a non-blood denominator chosen after that
  run with its E9.5 calls in view: a defensible change (the argmax calls were forced, mean probability 0.18 on SV-CM)
  but a later one. The scout's evidence also used Feng's own E9.5 hearts, which lie inside both counted windows,
  (8.5, 11.5]. Reading "before any run" as "before any scored run" would be a reading of `docs/DATA_POLICY.md`, which
  the delegation leaves to the captain. (2) Under the later rule a keep is still out of reach. As scorer sensitivities
  only, the digest's gap-fill 1 (an E8.5 -> E9.5 pseudo-task on a board near 75) gives +0.36 points per percentage
  point of emergent share with perfect states and 33-41% of that with nearest-progenitor proxies; the skill map
  (skill = d_floor / (d_floor + d), `docs/CHALLENGE.md`) makes the same change worth skill^2 as much, 0.35-0.44 of it
  at b0's 44.1 and 49.7. The rule's 1.9% of heart cells at t = 0.2 then gains about 0.25-0.5 on the chamber cut with
  perfect states and 0.07-0.22 with the states it builds, and less on the whole; a keep needs more than +0.3 and 8 of
  10 seeds ahead on both. The two counted tasks are nearly one result: the same predictions (inputs E9.5 and E8.5,
  build seeds 1,000,000 + s, 5,118 cells) scored against two cuts of one E10.5 library. A pass would exceed what the
  born cells can carry, and the born direction has no in-scope landing (target-free, E9.5 and Feng E14.5 only): it
  correlates best with NCC-derived (0.32), then Endocardium (0.26), NCC and BEC (0.25), and moves the cells 15-16% of
  the way to NCC-derived, 10-12% to Endocardium, BEC and ST, 2% to V-CM. No score pattern could tell such a gain from
  the scope and composition gaps the counted tasks already carry (the `configs/pseudo_tasks.yaml` notes; the T1 top-up
  entry above), so rule 2's third condition fails before the run. (3) The born state is a weak proxy: in the scorer's
  log1p pseudobulk the shifted cells reach 0.30-0.44 of Feng's within-E14.5 change; count-space edits keep zeros, so
  the 30 genes it raises most stay detected in 30% of cells against 82% in Feng's fibroblast-like; inside the chamber
  scope the nearest released cells (epicardium) correlate 0.20 with it (0.30 with chamber z). (4) Immune cannot be
  built: no released chamber cell detects more than 9 of its 30 markers (its own cells' 5th percentile is about 20),
  its nearest Feng type is blood (r = 0.07, 7 Feng E14.5 cells), and the s1 and s2 rules fix the blood share. (5) At
  test it is worth about zero or a little less (about -0.2, range -2 to +0.6). At t = 0.5 it adds 3.8-4.7% of heart
  cells: +0.3 to +0.6 if E12.5 holds that share, -0.5 to -2 if half or less (gap-fill 1: births at double the true
  share cost 4.4, against a reference it does not state). The confident-call test caught only 28-63% of E9.5's present
  non-cardiomyocyte nuclei (Proepicardium as epicardial: 244 of 883), so at E10.5 it would call fibroblast-like absent
  even at 3-7% of heart cells and add births on top. The share has no legal anchor: at E9.5, where both datasets are
  legal, the nuclei shares of the same types run 0.26x to 6.8x Feng's, and the three share bases move k by up to 3x.
  Births would go into both s1 and s2, so the s2 hedge does not cover them. Births are not added at test either; the
  births line is closed for T1 unless the captain rules otherwise or a new reason comes with its rule written before
  any run (`docs/DATA_POLICY.md` 2026-11-04 (d), which came after this entry), and T1's try count stays 6. Evidence
  on the DGX, not committed: `~/v/scratch/births/scout/` (`scout.py` sha256 387fda22, `scout_b2.py` 461f0959, post hoc; the per-cell
  calls), `critic_luckyfit/align2.py` 56e81d3c, `critic_value/reach.py` 0f7e21b3.

From the rebuild design of 2026-11-04 (wf_237f8d8f-e7c: five board designers, an adversarial critic each, a synthesis
and a compliance audit; design only, no scored run; the design texts are not committed). Pruned by the AI assistant
under the captain's delegation of 2026-11-02; the captain can overrule. Several diagnostics below were pruned because
they would read a counted task's truth and spend tries without making anything adoptable.

- **T1: a within-type maturity tilt** (importance re-weighting of real cells along a stage axis). It is composition
  under (d), not an expression change, and it fails on its own terms: its stage classifier separates the stages
  almost perfectly (AUC 0.991-1.000), the weights' effective sample is 0.12-0.28 of the pool, only 32-42% of the tilt
  is realised when the backtests draw 5,118 of 8,185 cells, and its Feng variant ranks cardiomyocytes by depth
  (rho -0.53).
- **T1: bracket transport of chamber sub-states, and Feng-internal held-out (FIB) checks.** Rule 2 fails before the
  run: Feng's capture changes from 69% LV-zone cells at E9.5 to zone-balanced at E14.5, so Feng-based growth and FIB
  brackets read capture, not development; item 5 applies as well.
- **T1: pricing the cardiomyocyte tail as a try.** A captain question only (item 6).
- **T1: a noise-aware, subtype-resolved Feng direction.** A split-half check of Feng's direction finds nothing to
  shrink, E9.5 has only 11 atrial cardiomyocytes, and its top-up is a variant of the pruned exact top-up.
- **T1: whole-ladder flow with growth.** 10-14 days of work, no check that reaches the test, and its critic's
  estimate at or below -0.5.
- **T1: a ladder-consistency gate and progress calibration of Feng's direction.** Its gate and clock read T1's
  sampling change (an embryo portion at E9.5, a dissected heart at E10.5).
- **T1: oracle ladders and decompositions on the counted tasks (eight diagnostics).** In-sample on a transition
  already read five times; they count as tries and can make nothing adoptable. Two had decision thresholds pinned to
  the field's validation totals, a use of leaderboard figures beyond calibration (`docs/DATA_POLICY.md` 2026-11-04
  (g)).
- **T1: PLAN item 3's decrease gate.** The form a backtest could test is not the form used at test.
- **T1: a composition-matched truth filter.** It judges only momentum variants that R1 makes moot.
- **T1: Feng's Endothelium group as a mislabelled identity.** Checked, no change: the group is 70 endocardial cells
  at E9.5 and 107 + 51 at E14.5, but on the pre-registered readout Feng's Endothelium direction agrees with the BEC
  change (sign agreement 0.813, Spearman 0.51, de_direction_full 0.26 against a 0.55 ceiling, where momentum scored
  -0.40).
- **Heart: Feng's direction, and the scRNA ladder's direction, on the heart panel.** Closed by A15 (no new
  cross-modal input enters an upload until the organisers answer); the target-free concordance check is null
  (cardiomyocyte Spearman 0.05-0.10 between Feng and MERFISH), and per-type momentum carries over only for
  cardiomyocytes.
- **Heart: removing the reference file's own transient per-gene deviation.** Leverage at the reference file is
  0.80-0.91, so the estimator is degenerate.
- **Heart: within-core state shares moved along Feng's change.** Closed by the captain's composition ruling of
  2026-11-01, and its sign was wrong.
- **Heart: two gap-anatomy diagnostics.** They would select on the counting backtests' truth.
- **Embryo: a depth-shaped draw of real bracket cells.** The depth line is closed; a simulation gives about +0.05 to
  +0.2 at test.
- **Embryo: exponential-growth state shares.** Share total variation against the current interpolation is 0.041,
  0.038 and 0.014 (about +0.05). As proposed, its check used an external atlas spanning E8.5, a held-out heart
  validation stage strictly inside a heart interpolation bracket, for a step in the shared core, before any window
  check: a window breach on paper (nothing was fetched or run; `docs/DATA_POLICY.md` 2026-11-04 (g)).
- **Embryo: s1 = fgw_aligned re-read on a test-like task.** It would re-decide E11, which its veto settled.
- **Embryo: displacement interpolation in a count-model latent space.** 7-10 days of work for about +0.15.
- **Embryo: a closability oracle, an attribution split, a third pseudo-task (heart E8.5) and a per-metric
  leaderboard read.** The first two read the counted targets and spend tries; the third would loosen the adoption
  threshold with a pair unlike the test; the fourth lies outside the screened route (totals and ranks only).
- **T3: a Mesp1-lineage timing shift on a WT-ladder basis.** The estimator fails a known-lag check, and the
  committed timing oracle is +0.17, +1.42 and -0.76.
- **T3: a same-reference transplant of Gata4's non-Mesp1 cells.** Within-label coverage is closed with T3's
  composition, and the build would apply the factor twice.
- **T3: a two-scenario factor exponent.** Its rho_F cannot be identified on released data.
- **T3: a per-gene Mesp1-lineage contrast shared by the released knockouts.** Nested in the timing shift above.
- **T3: three diagnostics on the cross-knockout folds** (s1's recipe with depth and gene identity separated, the
  knockout-slide variance net of the WT-run difference, an oracle bound of the Mesp1 response). They would spend
  both folds; K2a and K2b keep the one reading the principle asks for, on one fold.

From the DGX campaign design of 2026-11-05 (wf_7d1d7fe1-e6e: four program designers, an adversarial critic each, a
compliance audit and a synthesis; design only, no scored run; pruned by the AI assistant under the captain's delegation
of 2026-11-02, and he can overrule). Everything it proposed adds about +0.2 Overall in expectation after the critics,
and nothing passes the adoption rule on counting backtests it was not tuned on: the limit is legal evidence, not
compute (`docs/PLAN.md` decision log 2026-11-05 (b)).

- **T1: V3, Feng's E9.5 -> E14.5 epicardial-to-endothelial composition rate at f = 2/5, checked on a late-fetal family
  (Feng CD1 with Feulner 2026).** It re-proposes the pruned Feng-based growth without answering why it was pruned:
  Feng's capture moves from mostly LV-zone cells at E9.5 to zone-balanced at E14.5, and a late-fetal, whole-cell family
  cannot test that capture change or the transfer from whole cells to nuclei. The crosswalk of Feng's fibroblast-like
  group alone moves the test shift 2.7-fold, and the proposed power gate would preview the family's outcomes. Revival
  only at the captain's request, with a registered target-free capture check, every constant, the family under its own
  config key, a Theiler-to-E-day table and an X002 amendment committed before any gate. Before any row existed its
  designer compared the rule, one day forward from E9.5, with E10.5's committed non-cardiomyocyte split: an
  unregistered in-sample check that counts as T1's seventh try (rule 4). Its figures are carried nowhere, and E10.5 is
  never cited for V3 (`docs/DATA_POLICY.md` 2026-11-05 (b)).
- **T1: a hierarchical composition meta-model with scANVI.** Nothing identifies the (9.5, 13.5] gap, so it reduces to
  V3 with more parameters, and training scANVI on E10.5's labels reads a counted target.
- **T1: late-fetal pseudo-tasks (XP) as counting evidence; per-cell latent transport (JL-OT a and b); decoded
  sampling.** Feulner's TS22 sample is excluded by our staging rule (`docs/ORGANIZER_QUESTIONS.md` Q7's working
  answer; the X003 precedent for JAX embryos 50-53): a filter drops it at ingest, and it could be kept only with a
  documented per-sample age in (13.5, 14.5] and an organiser question citing Q7 sent first, since the captain cannot
  waive an organiser window. Every other legal late anchor is at or after birth, and XP0's gene-coverage check read the
  targets. JL-OT's counted expectation is 0.16-0.2, under rule 5's 0.3: an oracle per-type shift moved mmd_u only
  +0.12 and +0.30 points (`experiments/t1_scope/results/pseudo_E8.5_E9.5_lib1` and `lib2`); balanced transport would
  impose the anchor's capture-driven atrial:ventricular mix, a composition change with no quota rule
  (`docs/DATA_POLICY.md` 2026-11-04 (f)); and GPU non-determinism would break the freeze-morning hash check. Until a
  revival meets the audit's conditions (per-stage files with a loader test, a register row for Feulner with filter
  code, an X002 amendment per new Feng use confirmed by the captain, an independent leakage review of the harness
  extension, one late-fetal family approved as an adoption-rule reading), nothing from Feng enters a model beyond the
  ruled per-type direction. The two T1 designers saw per-stage cell counts and batch ids of Feng E17.5 and P0; a
  revival's rows disclose it.
- **T3: QC20, the knockout pipeline's 20-count and 20-gene cell filter in s1.** It claims more than the committed
  within-label ceiling (`experiments/t3_decomp/results/gata4/gata4_summary.csv`: at most about 1.4 and 0.8 T3
  points); its synthetic evidence was not built the way s1 is; its pseudobulk shift tracks pb (Spearman +0.97), so it
  is a depth increase; and its within-label redraws change T3's composition: two closed levers. Its check would read
  the Gata4 replicate, a fold that never counts. T3 stays at 3 tries, and nobody inspects E10.5's cell filter.
- **T3: additions to existing entries, not new prunes.** A Mesp1-lineage delay belongs to "a Mesp1-lineage timing
  shift" above (the WT axes are run differences). Generic cycle and stress programmes: no hand-typed gene list or
  pathway set enters a model without an `EXTERNAL_DATA.md` row, `docs/PLAN.md` 4.4 provenance excluding Ctnnb1
  perturbations, and A15. Shrinkage belongs to "a two-scenario factor exponent"; ensembles are covered by D3's
  factorial, D11b and K2a; a mean of runs is D10's x0 = 0.
- **Embryo.** E7.25 and E7.5 cannot count for any change to the mixing fraction: a designer derived each target's
  best fraction from the committed mixtures in `diag_emb_mix.json` (not quoted anywhere). The heart-E8.5
  "gate-estimate" diagnostic is dropped: it is the pruned third pseudo-task, and it would make the gate the reason for
  a try. The embryo s2 changed to fgw with aligned moves is a ruling request (R2, `docs/PLAN.md` decision log
  2026-11-05 (b)), not a row; its diagnostic N3 is registered only if the captain asks to consider it.
- **T1 and heart: population or share choices from qualitative E9.5-E13.5 background.** Closed by
  `docs/DATA_POLICY.md` 2026-11-04 (d) and decision log 2026-11-04 (e); struck.

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
