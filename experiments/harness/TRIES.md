# Tries ledger: every change tried for an entry, per test board

Started 2026-11-03 by the AI assistant under the captain's delegation of 2026-11-02, after the captain asked for
continuous iteration with no fixed timeline and a side agent warned that endless tweaking against the same few
backtests lets lucky fits through even when the adoption rule passes. The rule's 10 seeds rule out draw noise, not a
change that only suits the particular released stages a backtest uses.

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

## Counts per board

Before this ledger (from `results/registry.csv` and the decision files; floor and copy are reference rows, not tries):

| Board | Tries | Which |
|---|---|---|
| T1:test | 4 | b0, b0_scope, pseudobulk_shift (harness); m1_feng (labelled diagnostic, `experiments/t1_m1/results/diag_v1/`) |
| T2:embryo:test_interp | 6 | b0, pseudobulk_shift, psb_size, fgw (harness); fgw_shape0, fgw_shape1 (the s2 rule's run, `emb_shape_v1`) |
| T2:heart:test_extrap | 2 | b0, pseudobulk_shift (harness). The s2 hedge payoffs (heart scope) are a sensitivity analysis, not tries |
| T3:bcatenin | 1 | shift_transfer (the organisers' baseline); the T3 entries come from the pre-registered rules in `experiments/t3_decomp/` |

Since this ledger: see the table below. Running totals: T1 6, embryo 8, heart 3, T3 1 (H4 and E9 were pruned before their run; E10 is registered and counts only if its gate lets it run).

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
| E10 | embryo | `fgw_aligned` | fgw's coupling is an exact one-to-one pairing of its 2,000 anchors per stage: with Euclidean within-cloud distances the square-loss Gromov-Wasserstein term is concave on the couplings, so POT's conditional gradient ends on a vertex. Neither its features (panel expression at 30-212 counts per cell) nor the within-cloud distances can tell an embryo's left side from its right where expression is left-right symmetric (no lateralised marker before about E8.0), so partners can sit at random or mirrored across that direction, and McCann moves then pull both brackets towards its mid-plane: the moved cloud flattens there (wholly at t = 1/2 with mirrored partners), and so does the s2's mapped E8.0 half. The proper rotation is mostly unaffected (the other two axes and det +1 fix it). Change: keep FGW's rotation; take the moves (still the mean displacement of the 8 nearest anchors) from an exact transport of the same anchors on (1 - a) M / max M + a \|x - y R\|^2 / max, in the before stage's frame, a = INTERP_ALPHA = 0.5 reused (`vec.core.aligned_coupling`; `interpolate_clouds(..., move_coupling="aligned")`); in a common frame position tells the side. Synthetic left-right symmetric pair (tests/test_t2_build.py, AI-assistant scratch): FGW's partners have whitened cross-correlation about 0 along the symmetric axis, the aligned ones 0.4-0.7 (300-500 anchors), and the moved cloud's thinnest SD returns to the brackets'. **Gate, target-free, fixed here before it runs** (decided by the AI assistant under the captain's delegation of 2026-11-02, as the ideas workflow proposed): `experiments/harness/diag_emb_coupling.py` on the three bracket pairs at build seeds 1,000,000-1,000,009. A draw has the defect if fgw's smallest whitened anchor-partner cross-correlation is below 0.3 or its moved cloud's thinnest SD over RMS radius is below 0.9 x the thinner bracket's; the fix works if the aligned partners' cross-correlation is at least 0.6 and their merged-label agreement at most 0.05 below fgw's. The try runs only if draws with both reach 6 of 10 on each backtest pair and 3 of 10 (or the s1 seed's draw) on the test pair. A defect at test that the backtests cannot show: pruned, reported to the captain as a test-only risk. The diagnostic reads no target and no score, so it is not a try | E7.5: d2_shape up in at least 8 of 10 seeds, occupancy_dice flat or up; E7.25: d2 and Dice up in the seeds whose draws show the defect, about unchanged in the rest; de_score, de_direction, mmd_u and variogram identical to fgw's (the same cells); scale_log_ratio identical up to the scorer's rounding; neighborhood_mmd within +-0.01 skill on both | the gain sits in Dice alone, or only in E7.25's seeds 3, 4, 5, 8 and 9 (that block's anisotropy: E7.25's truth is thinner along PC3 than either bracket); or mean_anchor_shift falls well below fgw's while d2 moves toward fgw_shape0's (shorter moves across the sectioning-plane change, not the mechanism); or the gain is confined to draws where the rotation is a half-turn off (rotation and label agreement reported per seed). Caveats written before the run: the idea came after E7.25's seed split was seen; a = 0.5 is reused in a new role and no other value is tried; the backtests sit at t = 1/3 and 2/3, so a pass says nothing about t = 1/2 beyond the diagnostic | (this row's commit) | Pending: runs only if its gate passes |

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
- **T3: an absolute knockout-depth rule; the factor measured against the target's WT.** Both would gain on the two
  cross-knockout folds and change nothing at test, where the knockouts share the WT reference: lucky-fit traps.

## Diagnostics behind these choices

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
