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

Since this ledger: see the table below. Running totals: T1 4, embryo 8, heart 3, T3 1.

## Tries since 2026-11-03

Order of decisions, fixed before the run: E7 against fgw (the s1 entry); then E8 against whichever of fgw and E7 is
the entry after that decision. Both on the two embryo backtests, 10 seeds, by `experiments/harness/run.py decide`.

| # | Board | Candidate | Reason | Gain should show on | A fit, not an effect, if | Registered at | Outcome |
|---|---|---|---|---|---|---|---|
| E7 | embryo | `fgw_depth` | Detection depth differs several-fold between the released embryo stages (median counts per cell E6.75 30, E7.25 35, E7.5 62, E8.0 212; `results/diag_emb_mix/`), and the scorer reads log(1 + 1e4 c / N), whose zeros and levels follow depth. fgw's mixture carries both brackets' depths (at test half its cells at 62 and half at 212). The target's depth is predicted log-linearly between the brackets' medians at the bracket fraction (`vec.t2_build.DEPTH_RULE`; at test about 115), and cells of a deeper bracket are binomially thinned to it, the counts a shallower acquisition of the same cell would give. Positions and draws are fgw's | variogram and mmd_u above all, then neighborhood_mmd; de_score and de_direction may move either way (depth moves pseudobulk); the shape terms identical to fgw's | it gains on one task only; or the predicted depth happens to sit near one target's depth while the rule misses the other (log-linear gives 64 against E7.5's 62 and 49 against E7.25's 35), so the gain should hold on both | f68836d | **Not kept** (`results/emb_depth_v1/`, 10 seeds): E7.5 -4.18, ahead in 0 of 10 (veto: mmd_u 0.671 -> 0.598, variogram 0.664 -> 0.602, neighborhood_mmd 0.640 -> 0.557); E7.25 +0.74, ahead in 10 of 10. Why: a mixture drawn (1 - t) : t already has the target's average log depth (35^(2/3) x 212^(1/3) = 64 against E7.5's 62), so thinning the deeper bracket pulled the population's average below the target; on E7.25 the mixture's average (49) sat above the truth (35), and thinning moved it toward it. What the scorer rewards is the population's average depth, which the mixture already gets right at the nominal fraction |
| E8 | embryo | `fgw_blend` | Development is continuous: an E7.75 cell is an intermediate state between its E7.5 ancestor's and its E8.0 descendant's, which a mixture of whole cells of the two brackets does not contain. Each drawn cell's counts become a binomial blend of itself and a partner from the other bracket chosen through fgw's coupling, splitting the expected counts (1 - t) : t at E7's predicted depth (a shallow side that cannot supply its share is capped at all its counts) | mmd_u, variogram and neighborhood_mmd beyond E7's; de_score and de_direction about unchanged; the shape terms identical to fgw's | it gains on one task only; or the gain comes from the blend's narrower within-state spread suiting one target (a blend averages two cells, so it carries about half a real cell's cell-to-cell variation) | f68836d | **Not kept**: E7.5 -2.23 and E7.25 -3.19, ahead in 0 of 10 on each (mmd_u 0.671 -> 0.638 and 0.497 -> 0.429; neighborhood_mmd falls on both). Blended cells are too narrow a population: the targets' cells are discrete real states, which a mixture of real cells matches better than averages of two |

| H3 | heart | `b0_depth` | (Restated after E7 and E8 and before this try's run: what the scorer rewards is the population's average log depth; E7 failed because it moved that average away from the target, and gained where it moved it toward the truth.) Counts per cell fall steadily through heart development in the released series (median counts per cell: E8.25_late 233, E8.5 166, E8.75 130, E9.5 75, E10.5 56; `results/diag_depth/`), and the scorer reads log(1 + 1e4 c / N), whose zeros and levels follow depth. B0 keeps E10.5's depth. The target's depth is extrapolated at B0's own damping: half the latest log rate per day, carried over the horizon (`vec.t2_build.DEPTH_EXTRAP_RULE`; at test 56 x (56/75)^(0.5 x 2) = 42); B0's cells are binomially thinned to it, never deepened. Leave-one-out on the depth series alone, each stage from the two before it: E8.75 134 (observed 130), E9.5 90 (75), E10.5 52 (56) | variogram, mmd_u, neighborhood_mmd; de_score and de_direction may move either way; the shape terms identical to B0's | it gains on one task only. Caveat written before the run: the survey that found the trend also showed both backtests' target depth (E10.5's 56), so the form of the rule (per-day rate, B0's damping 0.5, the size procedure's convention) was written knowing it, and E7's outcome was seen before this row was restated; the damping is B0's existing value, not fitted, and the leave-one-out above is the check that does not use a scored target | this commit | pending |

## Pruned before running

- **Partial or unbalanced FGW coupling** (the plan's "mild unbalancing"). Its reason would be tissue that one bracket
  sampled and the other did not. But the embryo shape loss on the E7.5 backtest comes from the sectioning plane: z is
  the long axis at E6.75 and E7.25 and the thin axis at E8.0 (`experiments/coords/REPORT.md`), and the stack length is
  set by acquisition. Leaving mass unmatched does not address that; the s2 shape variant already hedges it.
- **A developmental-clock warp of the bracket fraction t.** On the backtests the targets sit off their nominal
  fractions in opposite directions (E7.25 near E6.75; E7.5 about midway between E7.25 and E8.0 by mmd_u), and at test
  the nominal fraction is 1/2. Nothing measured on released stages says which way E7.75 would sit.
- **Sweeps of FGW alpha, k or the anchor count.** No physical reason for another value than the plan's.

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
