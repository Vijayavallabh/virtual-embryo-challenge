# Non-mimetic levers for a top-3 place: screened cross-disciplinary research

> Date 2026-11-05 | AS_OF 2026-11-05 | Standard mode.
> Sources: 131 registry entries, 128 distinct works, 66 hosts, 91% academic or official; arxiv.org hosts make up 26%
> once merged, just over the 25% gate.
> Eight specialists, each notes file screened and audited. Three target-free diagnostics (N3-N5), run at 1edec56 on a
> clean DGX tree.
> Counter-reviewed by four adversaries (50 issues; section 8), and revised.

Written by the AI assistant for the captain's request of 2026-11-05: research relevant topics across disciplines from
first principles, then develop non-mimetic solutions that maximise the chance of a Human-track top-3 place.

No agent saw a test target or a test-board limit. The frozen count rules read the limits in code and print none. No try
was run. The research agents' compliance logs record four slips, listed in `docs/DATA_POLICY.md` 2026-11-05 (k):
- a memory file opened;
- a scorer docstring read that may hold a hidden-target count;
- one source with in-window results, used for its method only;
- toy runs on the loaded local machine.

The audited notes are in `notes/` (task-a to task-h) and the citation registry in `notes/registry.md` ([n] below).
The diagnostics are in `diagnostics.py`, with stamped outputs in `results/` (commit, clean tree and input sha256s).

## Executive summary

Eight specialists worked the problem from eight fields:
- statistical decision theory;
- statistical physics and systems biology;
- the provenance of pathway resources;
- computational statistics;
- contest economics;
- dynamical systems;
- tissue physics;
- assay measurement science.

**The bottom line holds after the counter-review: nothing found moves our total into its own top-3 region.** Reaching
our 90th-99th percentiles of reading A (about 204-208) takes 12-17 points over its mean of 191.5. The open levers come
to roughly +1 to +3 Overall in expectation, and several compete for the same upload slot. No leaderboard figure is used
as a threshold here; the screened calibration stays in `docs/PLAN.md`'s decision log as calibration only.

**Two readings of T3 give two pictures of the top-3 chance.**
- *The gate estimate* (`experiments/gate/estimate.py`, reading A) gives T3 a mean of 64.6. Half its weight is on an
  in-sample stand-in at 73.3, the Gata4 replicate fold, which never counts. Under it, T3 holds about 84% of our total's
  variance, and the top-3 chance is a few percent.
- *Task h's measurement-based payoffs* put s1's best case near 66 and the pair's mean near 58.7. Under them the chance
  of reaching our 95th percentile is under 1% (counter-review).
- Either way, the chance rests on s1's T3 and T1 bets landing near their best case.

**The sharpest finding is a quantification of a bet already on record:** s1's T3 severity rests on depth.
- The facts:
  - the scorer's severity gate rises with any common, one-signed part of a prediction;
  - Gata4's shallow depth is a property of its imaging run (it shows in cells the knockout cannot touch);
  - s1's depth step carries 43% of its predicted change's energy (N4).
- If the beta-catenin file is not as shallow as Gata4's, s1's severity, a quarter of T3, falls to its floor.
- Even if it is as shallow, s1 clears the gate only if its factor also agrees gene by gene. That has been seen only
  in-sample.
- A Gata4-free backup with a factor and no depth would cover part of the losing side. It is worth about +0.3 T3 under
  the default reading of a floored severity (skill 0.5). It is worth several points if a floored severity scores 0,
  which our working answer gives a 0.35 chance.
- So it is not settled. Section 3 lists what would decide it, all of it target-free, before T3's first upload on
  2026-11-27.

**The other open moves are captain's rulings and small:**
- an embryo backup leaning toward E7.5, which also hedges the known risk that E7.75 was imaged at E7.5-like depth
  (about +0.6 Overall from a toy; it amends a frozen s2 rule after the fact, as R2 would have);
- a skewed T1 bet on a "strong-null" test regime (+1.5 to +2 T1 if right, -0.2 to -0.45 if not);
- smaller T1 realisation and hygiene changes, two of which need gene sets that the open A15 question may close.

**Closed or pruned, by measurement or by argument:**
- the T1 tail decoder: N3 finds 5% of s1's bottom 400 genes infeasible;
- decoders for T3 and the heart (N5);
- knockout prediction from wild-type covariation (expected +0.02 to +0.05);
- clean TCF/LEF motif priors (about 0);
- time warps and step-size tuning on T1;
- within-upload mixtures;
- latent Gaussian-OT decoding, which is a composition change in disguise;
- sorting-physics extrapolation on T2.

## Contents
1. Where the top-3 chance lives
2. What the scorer rewards, from first principles
3. T3: the depth bet inside s1, and the backup question
4. T3: knockout responses from wild-type variation and clean pathway priors
5. T1: realisation, tails and dynamics
6. T2: the embryo backup, the heart's momentum, and the tails
7. Ruling sheets for the captain, by deadline
8. Counter-review
9. Key findings
10. Limitations and future directions
11. References

---

## 1. Where the top-3 chance lives

**A trailing team in a small-winner contest should add variance and cut correlation with the field.** This holds when
the winners are fewer than half the contestants [66]. In winner-take-all forecasting contests, honest forecasting is
not an equilibrium: contestants move away from the common prior [67]. Rank-based selection rewards exaggerated reports
[3][4]. In simulated tournaments the usual winner is a good forecaster whose errors fell the right way [71]. Field data
show trailing players take more risk, with mixed evidence that it pays [76][77]; [76] and [77] were read only as search
summaries.

**Our numbers, under two readings of T3.**

| | Reading A (gate estimate) | Task h's T3 payoffs |
|---|---|---|
| E[T3] for today's pair | 64.6 (stand-ins 73.3 / 60.8 / 50.1 at weights 1/2, 1/4, 1/4) | about 58.7 |
| T3 SD; share of Var(Overall) | 9.5; about 84% | about 7.5; about 77% |
| T3 outcome that matters for our 95th percentile | about 73-81 | s1's best case (about 66) lies below it |

- *Board SDs in reading A:* T1 3.2, heart 4.5 and embryo 2.6. T1's share is a lower bound, because reading A leaves
  out several unmeasured T1 risks.
- *Value of variance for a laggard:* certainty equivalent = mean + kappa x variance, with kappa = (c - mean)/(2
  sigma^2). That is 0.020, 0.038, 0.056 and 0.073 Overall points per unit of variance at c = 196, 200, 204 and 208,
  all within our own range. A zero-mean bet with SD 4 is then worth +0.3 to +1.1. It is worth about +1 only near our
  95th-99th percentile, and the arrays' thin right tail makes the normal formula an overstatement.
- *Task a's z-rule* took its threshold from a leaderboard figure, so it is struck as a decision input
  (`docs/DATA_POLICY.md`, Scores).

**Rank, not level.** The top-3 line is set by the field's third total, which moves with the hidden targets as ours
does. Task e's mimicry discount lambda is the share of the field's common shock that our s1 shares. Near lambda = 1,
only the second upload's gain over the first creates rank variance, so bold, uncorrelated backups gain value. A
synthetic bold T3 backup lifted the chance from 3.4% to 8.2% at lambda = 1, against 7.6% to 8.6% at lambda = 0.

Our T3 s1 transfers Gata4's delta, which is what informed teams would also do (task e, table F). If the beta-catenin
file is not Gata4-shallow, that hits much of the field at once, which raises lambda for T3. lambda cannot be read
without competitor material, so it is the captain's judgement. The rank-aware pricing tool (rule K3, task e) is
target-free to compute, but acting on it needs a ruling.

**Confidence:** Medium.

**Basis:** classical contest results [66][67][71] and committed arrays. Top-end probabilities are fragile under both
readings.

**Opposing interpretation:** validation-to-test shifts can preserve rankings [73]. Then lambda is near 0 and only mean
gains matter. Against that, the overfitting meta-analysis found large gaps exactly under non-i.i.d. shifts [72], which
is our case.

---

## 2. What the scorer rewards, from first principles

Task a derived the expected-skill-optimal submission for each pinned metric (veckit cc3859a; formulas only).

**mmd_u, variogram and neighbourhood MMD are kernel-score divergences.** The expected squared MMD is minimised by the
posterior predictive mixture [1][16][102][103]. The scorer bends this in two ways: it refits the basis and bandwidth on
each truth, and its skill map's floor distance shrinks when the truth sits near the baseline. That map acts as a
relative error, whose Bayes act is shaded toward the reference [2].
- *Within one upload,* mixing over continuous shift sizes buys nothing in a toy.
- *Across two uploads,* a two-point pair gains +0.08 to +0.11 skill. T1's frozen pair already has this form. The
  embryo pair does not: it shares one mixture (section 6).

**de_score** fills the predicted top n_up and bottom n_dn by dp, against a null ranked by reference level. The
probability ranking principle ranks by P(gene in the truth's list) [6][8][9], shrinking noisy genes [5]. A truth
down-gene needs pb(reference) >= 0.25, so a predicted decrease of a lower gene is a guaranteed miss if it reaches the
bottom tail. N3 and N5 measure this on every upload file (sections 5 and 6).

**de_direction** is a rank partial correlation, invariant to magnitude. Predicted and true deltas share the scorer's
10% reference-copy noise, so the optimal scale is the posterior mean. Inflation buys nothing; strong shading costs
0.02-0.15 in r.

**severity_slope (T3)** has an exact gate:

  r^2 > 0.01 if and only if cos^2(dp, dt) > 0.01 + 0.99 cos^2(dp, 1) over the truth's DE genes.

A predicted common shift raises the bar. A toy with 60 DE genes found:
- a predictor carrying the full true shift with small gene-specific agreement passed the gate in only 0-5% of draws;
- scaling the shift by that agreement passed in 70-100%.
Whether a floored severity scores skill 0.5 (the live board) or 0 (the Evaluation page's "undefined") is open. Our
working answer gives 0.5 a 0.65 chance (`docs/ORGANIZER_QUESTIONS.md`).

**Confidence:** High for the derivations and the gate identity; Medium for toy sizes.

**Opposing interpretation:** the refitted kernel and truth-dependent floors break exactness, so sizes must come from
diagnostics or registered tries.

---

## 3. T3: the depth bet inside s1, and the backup question

**What was known.** s1's knockout-depth term was already logged as "a bet, not a measurement" (`docs/PLAN.md` decision
log 2026-11-04 (c)). The ledger keeps depth out of cross-knockout candidates because "the knockouts' depth ratios flip
between folds" (`experiments/harness/TRIES.md`, K2a). What is new is the size.

**What the assays say (task h).**
- MERFISH detection efficiency depends on the imaging run and on local spot density [117].
- Panel-library normalisation can flip fold-change signs [118].
- Per-cell depth tracks tissue structure as well as technique [119].
- Segmentation choices alone can nearly double the cells called [122] and confound differential expression [121].

**What the released files say** (committed figures):
- Gata4's per-label depth ratios against WT E8.75 are 0.40-0.73 in non-Mesp1 labels as well as Mesp1 labels, so the
  deficit is a property of the run.
- The Gata4 stacks hold 1.30-1.48 times WT's objects per level at about 0.8 times the counts. That is either finer
  segmentation or more tissue captured less densely; task h's structure audit D2 would decide which.
- Every released knockout file is shallower in absolute terms than WT E8.75 (Gata4 82 and 68, Mab21l2 95, against 130
  median counts). Against WT E8.75 every one would sit at 0.52-0.73x: Gata4-shallow.

**What the scorer does with it** (N4, `results/t3-gate.json`; our own builds against WT E8.75 as released):
- The rebuild of s1 at the upload seed matches the archived file (pseudobulk correlation 1.0).
- The depth step carries 43% of s1's dp energy (an energy ratio), with rank correlation -0.945 with WT expression level.
- The gate bar depends on which genes the truth calls DE. It is 0.22 over all genes, 0.31 over down-feasible genes and
  0.66 over s1's own most-changed genes. On a real truth's DE genes, the committed offset share of the s1-like
  `xrep:F+depth` on Gata4 is 0.48-0.52.
- Without the depth step, the bar is 0.03.

**So s1's severity needs two things at test:**
1. a matching depth deficit (U2);
2. gene-specific agreement of its factor.

Both are needed: on Gata4, the target's own depth without a factor (`fac:depth`) floored in 20 of 20 seeds, and only
in-sample factors passed. Cross-knockout factors have shown no such agreement (K2a; A10 killed).

**The unknowns.**
- U1: was the beta-catenin file imaged apart from the WT reference? Task h's prior is about 0.85.
- U2: does it share Gata4's shallow depth? A defensible range is 0.4-0.9, not a coin flip. The absolute depths above
  point to a shallow file; the runs' different mechanisms allow a fresh draw.
- U3: does the shared stage offset carry? The barred factor's carry to Gata4 is mostly the WT E9.5-against-E8.75 axis:
  without it (`vs95`) it scored 52.04 and 46.24.
- The floor reading: does a floored severity score 0.5 or 0?

**The backup question (lever 1).** Make T3's s2 the barred branch's s1: the Mab21l2 factor against WT E8.75, no depth,
the same central composition. It keeps the Gata4-free insurance.

| Pricing | Expected T3 gain of the swap | Source |
|---|---|---|
| Task h's payoffs, before U3 | +1.5 to +2.1 over U2 = 0.55-0.65 (58.65 to 60.43 at 0.6) | task h; the arithmetic review corrected the figures |
| With U3 priced and the barred file valued unconditionally (about 49 + 9 x P(its gate passes), so 53-55) | about +0.3 (-0.2 to +0.9) | T3 adjudicator |
| The review of 2026-11-05, coupled model | -0.22 to +0.11 | `experiments/gate/results/t3_swap_review_2026-11-05.md` |
| If a floored severity scores 0 (chance 0.35) | several points positive (the critic's +2.6 to +6, before the adjudicator's lower price for the barred file) | completeness critic |

Today's s2 bets on A11's composition, which is killed. That composition alone raises its own gate bar to 0.35, which
is why s2's stand-in floored severity in 14 of 20 seeds on the released Gata4 replicates. Under a floor = 0 reading,
each floored upload loses another 12.5 T3 points, so today's s2 would average roughly 42 instead of 51
(illustrative: 51 - 0.7 x 12.5). That is what makes the swap robust there.

A side agent raised the same point to the captain on 2026-11-05. Today's s2 is T3's first upload, so this reading
weighs on the decision due before 2026-11-27.

**What still binds** (the review of 2026-11-05 and this counter-review):
- The swap would override D11b's failed pre-registered bar (59.40 < 60) after the fact; that bar, through rules.yaml
  note 6, is what made the refused branch the s2.
- It has no counting fold. Its anchors include the replicate fold, which never counts.
- U3 is unmeasured.
- The ruling switches must be rewritten, though they can be.
- The never-scored file matters little: the composition's own dp is about an eighth of the factor's.

**What would decide it, all target-free, before T3's first upload (T3's s2 goes first, from 2026-11-27):**
1. Task h's D5 and D6. D5 asks whether the component shared by s1's factor and D11b's lies along the WT developmental
   axis (staging) or across it (a technical anchor). D6 places WT E8.75 on the WT ladder.
2. D1 and D2, the run fingerprint and areal capture audits, which bear on U2.
3. The rank-aware pricing (rule K3) under both floor readings.

Decision rules are committed before any of them runs. K3a, a counted cross-knockout try registered in TRIES first, is
the only counted reading of the mechanism. Recommended default until then: keep today's s2.

**The checks, run 2026-11-05 at the captain's request** (TRIES N6-N10, registered with their rules before any run;
`results/d1.json`, `d2.json`, `d5d6.json`, `pricing.json`):
- **D1, run fingerprint.** The Gata4 replicates are on slide M164 (batches 1 and 2) and Mab21l2 on slide M125 (batch 0).
  They are two independent knockout runs.
- **D2, areal capture.** Object density is about the same in every file. Gata4's shortfall against WT E8.75 is lower
  capture per area (counts per area 0.56-0.69x), not finer segmentation. Every released knockout file is shallower
  than WT E8.75 (0.52-0.73x). Under the rule, p_U2 falls in the middle band, 0.7.
- **D5, the shared component.** 63% of the two factors' shared part lies across the WT developmental axis (80% on the
  secondary axis), which gives the high band, p_anchor 0.65. But the shared part is weak in size: the two factors
  correlate only 0.15, and the barred factor correlates 0.42 with the developmental axis itself.
- **D6, staging.** WT E8.75 sits at its nominal place on the ladder (0.400 against 0.404). Half of the known-lag check
  failed, so by its rule D6 shifts nothing.
- **The pricing at (0.7, 0.65).** Swapping gains +0.79 T3 under the half reading, +2.47 under the zero reading and
  +1.38 under the mixture. It raises the top-3 chance at our 90th and 95th percentiles by +1.9 to +2.5 points if the
  field shares our s1's fortunes (lambda 1), and by about 0 if it does not (lambda 0).
- **The rule's verdict: recommend the swap.** It holds at p_anchor 0.5 (+0.94) and fails at 0.35 (+0.48).
- **Caveats.**
  - D5's rule reads the split of the shared part, not its size. And D6's raw placement (W875 on its nominal stage)
    suggests Gata4's "older" look was its own litter's, which argues for a p_anchor nearer 0.5 than 0.65; the swap
    still clears the bar there.
  - Most of the gain is the zero-credit reading, which the organisers have now been asked about (fifth email,
    2026-11-05). Under "half" alone, the gain is +0.4 to +0.8 at p_anchor 0.5-0.65.

**Confidence:**
- High that s1's severity hinges on depth and on its factor's agreement.
- Low on every payoff figure.
- Medium on the recommendation.

**Opposing interpretation:** beta-catenin and Gata4 share lab, stage, driver and WT reference, and every released
knockout file is shallow. U2 may then be likely enough (break-even about 0.77-0.9) that today's pair is right as it
stands.

---

## 4. T3: knockout responses from wild-type variation and clean pathway priors

**Fluctuation-response (task b).** A wild-type covariance predicts the subspace and the gene-wise susceptibility of a
response, and its direction once one scale per perturbation is fitted [17][19][20][22][27]. Zero-shot sign and size
need exogenous variation of the regulator, which a snapshot lacks. Further limits:
- at sparse detection, slopes are attenuated or, after correction, confounded by cell state [23];
- network inference given interventions did not beat observational methods at scale [25];
- common metrics reward systematic variation [29];
- MERFISH admixture creates spurious co-expression [28].

On the released knockouts:
- the within-type residual's projection shares are -0.15 to +0.10 of dt;
- the oracle ceiling is +3.8 to +4.5 T3 points;
- the shared factor plus depth carry 0.78-0.80 of the decomposed Gata4 delta.

Task b puts the realistic gain at about +0.2 if it works. With a chance of passing near 0.1, the expected value is
+0.02 to +0.05 Overall. Not worth a try now.

**Clean pathway priors (task c).**
- *Excluded by provenance* [32][33][37][38]:
  - MSigDB's Wnt/beta-catenin hallmark;
  - perturbation footprints;
  - curated regulons of the held-out regulators.
- *Clean* [39][41][42]: sequence-derived matrices (JASPAR CORE, HOCOMOCO v12).
- *A judgement call* [44][45]: binding atlases.
- *Size:* a motif-only prior is worth -0.4 to +0.2 T3 points. A binding-informed one is worth +0.2 to +0.8 only under
  optimistic accuracy, and it loses on an entry that already carries signal.
- *Main value, defensive:* a provenance gate that keeps the common tools' default resources, which are the excluded
  ones, out of any build.

**Confidence:** Medium-High.

**Opposing interpretation:** CIPHER's strong fits [17] suggest more is possible. But its magnitudes are fitted and its
signal sits in about three global modes, which here are the factor and depth s1 already carries.

---

## 5. T1: realisation, tails and dynamics

**The tail decoder, pruned by measurement** (N3, `results/t1-tail.json`):
- None of s1's bottom 100 genes by dp is infeasible. 5 of the bottom 200 are, 20 of the bottom 400 and 87 of the
  bottom 800; all but two come from the chamber re-weighting.
- T1's s2 has fewer: 5 of its bottom 400 (N5).
- The released E9.5 -> E10.5 transitions had 382-426 truth down-genes, an in-sample property of a counted truth.
- At that depth, replacing the wasted slots is worth about +0.05 to +0.25 T1 points at a hit rate of 0.1-0.4, at most
  about +0.6. All of it would need a top-up of scope-driven genes, with a variogram cost.
- Pruned. The pruning holds only if the test's n_dn stays near 400.

**The count kernel for increases (task d).**
- Binomial thinning realises decreases exactly [51][52], and UMI counts fit a Poisson measurement model [53].
- `vec.core.apply_factors` adds increases as c + Poisson((f - 1)c). That never raises detection, and for sparse genes
  it realises about half an intended rise.
- An empirical-Bayes posterior-predictive top-up is exact in law; T3's frozen kNN kernel already is near-exact.
- But the ledger already pruned "the same exact top-up on T1", because the counted backtests judge it through the
  cardiomyocyte-share confound. Its expected value is +0.1 (range -0.3 to +0.5).
- It is a revival of a pruned line: a captain's blind bet only. It would also reach the heart's B0 through the shared
  core, unless confined to T1 by a flag.

**Dynamics (task f).** With a sampling break at E10.5 (the dissected heart), today's M1 after R1 is the state-space
optimum:
- level = the E10.5 cells;
- slope = Feng's within-dataset direction;
- horizon = the bracket fraction.

Damping is a trend-persistence probability [83], and the Theta benchmark is a half-damped drift [84][85]. Closed:
- moderately curved developmental clocks move the bracket fraction by about +-0.04 from 2/5; strongly curved ones reach
  about 0.32, still inside the range where T1 is flat (0.28-0.5 of the bracket);
- step size (flat to +-0.2 points);
- cell-level Wasserstein momentum (paths cross at s = 2);
- splines, which overshoot [90][91].

**Open T1 bets, all captain's rulings:**
- *L3, a strong-null prior in the down tail (task a).* The E9.5 -> E10.5 transition was "strong-null": mostly highly
  expressed genes going down. If the test regime matches, a rule fixed before any run gives +1.5 to +2 T1. If not, it
  gives -0.2 to -0.45 plus a distributional cost. Expected value is positive if P(strong null) is above about
  0.2-0.3. No counted check exists. It is the one skewed T1 bet, the kind a laggard should prefer.
- *A maturation quota on cycling cells (task f).* Expected +0.2 to +0.5; -0.3 to -1.5 if the sign is wrong.
  Conditions:
  - its rule is committed before any read of E10.5's within-category shares (DATA_POLICY 2026-11-04 (f)(1));
  - the drawn median library size is reported, and a shift above 5% goes to the captain ((f)(3)), as task f expects;
  - the cell-cycle gene set needs an EXTERNAL_DATA row and PLAN 4.4 provenance;
  - the open question whether A15's gene-set line covers T1 sets needs the captain's reading.
  Signature-based growth rates depend on the database [89]; growth and differentiation cannot be separated from
  snapshots [87]. [87] carries in-window results; only its method point is used here.
- *Feng direction hygiene (task h):* hold flat dissociation-stress, ambient and nuclear-detectability modules
  [120][125][126][127][128][129]. Expected +0.1 to +0.5. It needs gene-set rows and the same A15 reading.
- *Shading Feng's step toward the reference (task a L5, +0.3 to +0.5)* conflicts with any bolder-step lever. The bold
  T1 backup of task e (60% cardiomyocytes plus an exaggerated step) is worth only about +-0.2, because T1 is flat in
  step size, and it repeats D3, which lowered our chance at 8 of 9 levels. It is not recommended.

**Confidence:** High for the tail result and the closed dynamics levers; Medium for the bets.

---

## 6. T2: the embryo backup, the heart's momentum, and the tails

**Arrangement levers are spent (task g).**
- Permuting rows over positions moves only the neighbourhood term, which in effect is a normalised pseudobulk-shift
  metric.
- Content dominates arrangement (`experiments/density/REPORT.md`).
- Balanced, coherent or patchy assignment changed nothing or lost points in a toy.
- Sorting kinetics are slow and model-dependent [101][106][107]. Gromov-Wasserstein is blind to mirror images [114].

**An embryo backup leaning toward E7.5 (lever 2).**
- The skill denominator shrinks when the truth sits near copy_last, so a prediction at the nominal fraction 1/2 loses
  most when E7.75 sits close to E7.5.
- In closed form, skill(t; t') = t'^2/(t'^2 + (t - t')^2).
- A backup at t = 1/4 gains +0.10 skill per mixture-sensitive group (+2.4 points per 25-point group) on average over a
  uniform t' in [0.2, 0.8], positive only for t' below 0.375. The toy's central figure is +1.2 embryo points (+0.6
  Overall).
- Today's s2 shares s1's mixture and differs only in geometry, so the "near E7.5" case is unhedged on half the board.
- *A second reason, not priced in the toy:* a t = 1/4 mixture is three-quarters E7.5 cells (depth 62 against E8.0's
  212). It therefore partly hedges the logged, unhedged risk that E7.75 was imaged at E7.5-like depth (`docs/PLAN.md`
  decision log 2026-11-04 (c)).
- *Against it:*
  - it amends a frozen s2 rule after the fact, the ground on which the captain declined R2;
  - E7.25's closeness to E6.75 is barred as evidence for a mixing-fraction change (TRIES rule 2; DATA_POLICY 2026-11-05
    (b)(v)), so the scenario prior must come from the captain without the backtest positions.
- A target-free table, s1 and the declared backup scored against pseudo-truths built from the test pair at
  t' = 0.1-0.9, would show the asymmetry before he rules.

**The heart's momentum crosses the sampling break (task f).**
- At test, B0 moves each lineage by half of E9.5 -> E10.5, an increment across the E10.5 dissection change.
- The counted heart backtests use pre-break increments, so they never scored the estimator used at test.
- A hold, or cardiomyocyte-only momentum, removes that untested estimator: about +-0.25 Overall, expected about +0.05.

**Coherence (task e, lever 3).** For a laggard, a positive correlation between board outcomes raises the top-3 chance
at no cost to the mean. T1's s2 (cardiomyocyte-rich) and the heart's s2 (heart-focused trim) already bet on the same
"more chamber-dominated E12.5" world. Keep any future change on that side. This needs no ruling.

**Tails on T2 (N5).**
- *Heart:* 1-3 of its bottom 40 genes are infeasible, and 13-17 of its bottom 80.
- *Embryo:* 1 of the bottom 40, but 31 of the bottom 80. An embryo tail decoder would matter only if the truth's
  n_dn is 80 or more.
- Either decoder would need count edits to real cells on a half-weight board, so neither goes on the menu now.

**Confidence:** Medium, from a toy and a closed form.

**Opposing interpretation:** if E7.75 sits at or past its nominal fraction, as E7.5 did, the backup adds nothing and
dilutes today's geometry edge.

---

## 7. Ruling sheets for the captain, by deadline

Each sheet gives what changes, the recommendation, the deciding prior and break-even, the loss if wrong, the rule
overridden, the deadline and the default if he says nothing.

**A. T3 s2: today's refused s1, or the barred s1 (lever 1).**
- *What changes:* the T3 s2 file (`t3_build/barred_s1.h5ad`, already built and format-checked).
- *Recommendation (after the checks, by the rule fixed before they ran):* swap: T3's s2 becomes the barred branch's s1
  (`t3_build/barred_s1.h5ad`). It is Gata4-free, already built and format-checked. If the organisers answer "half"
  before 2026-11-27, the case weakens to +0.4 to +0.8 and becomes a close call. If they answer "zero", it strengthens
  to about +2.5.
- *Deciding priors:* P(U2), P(the stage anchor carries) and P(floor = 0). Under the default floor reading the swap
  breaks even near P(U2) of 0.77-0.9, or near P(anchor carries) of 0.35. It is robustly positive under the 0 reading:
  today's s2 floored severity in 14 of 20 practice seeds and would average roughly 42 rather than 51. Our working
  answer gives the 0 reading a 0.35 chance. A one-line question to the organisers on this reading is worth sending
  (draft on request).
- *Loss if wrong:* about 1-3 T3 points.
- *Overrides:* D11b's failed bar, through rules.yaml note 6, and the T3 ruling switches.
- *Deadline:* before T3's first upload, from 2026-11-27 00:00 UTC.
- *Default:* today's s2.

**B. Embryo s2 at mixture fraction 1/4 (lever 2).**
- *What changes:* the embryo s2 file. Geometry, size rule and count are unchanged.
- *Recommendation:* build the target-free pseudo-truth table, then decide. The AI assistant leans yes, because of the
  depth hedge.
- *Deciding prior:* P(E7.75 near E7.5 or E7.5-like in depth). Break-even is about 0.15-0.2.
- *Loss if wrong:* about 1.5 embryo points (0.75 Overall).
- *Overrides:* the frozen embryo s2 rule, after the fact (R2's precedent).
- *Deadline:* the embryo's second upload, by 2026-12-01 06:00 UTC.
- *Default:* today's s2.

**C. T1 strong-null down-tail prior (task a L3).**
- *What changes:* T1 s1's expression, so its s2 returns to the captain.
- *Recommendation:* the captain's call, as a skewed bet.
- *Deciding prior:* P(the test regime is strong-null like E9.5 -> E10.5). Break-even is about 0.2-0.3.
- *Loss if wrong:* about -0.2 to -0.45 T1, plus some distributional cost.
- *Overrides:* M1, which is his ruling.
- *Deadline:* the s1 freeze, 2026-11-27 00:00 UTC.
- *Default:* no change.

**D. T1 maturation quota, and Feng direction hygiene.**
- *What changes:* T1 s1's within-category shares, or Feng's direction.
- *Recommendation:* only after the captain rules on whether A15's gene-set line covers T1 gene sets. If it does not,
  register the rules first; if it does, both are closed.
- *Loss if wrong:* up to about -1.5 T1 for the quota.
- *Deadline:* the s1 freeze.
- *Default:* no change.

**E. The exact up-kernel (revives a pruned line), and the heart's momentum hold.**
- *Expected value:* +0.1 and +0.05.
- *Recommendation:* decline both, as below the line.
- *Default:* no change.

**Decision tools, target-free, to run next:** D5/D6 and D1/D2 for sheet A, the rank-aware pricing for A and B, and the
embryo pseudo-truth table for B. About 2-3 days in all, well before 2026-11-27.

---

## 8. Counter-review

Four independent adversaries reviewed the first draft (wf_3eef9333-d9d): arithmetic and traceability, a T3
adjudicator, compliance and governance, and completeness. They raised 50 issues: 1 critical, 13 high, 20 medium and 16
low. Their probabilities that adopting lever 1 as first framed would be a mistake were 0.6, 0.55, 0.55 and 0.4. All
four held that the bottom line stands and that the numbers did not.

**Key controversies:**
- **Lever 1's value.** The draft took task h's +1.8 "before the U3 risk" and dropped the qualifier.
  - Pricing U3 and valuing the barred file unconditionally gives about +0.3 (adjudicator).
  - The floor = 0 reading, which the draft never considered, makes the swap robustly positive (completeness critic).
  - Both are now in section 3, and the recommendation is to measure before ruling.
- **Two T3 models.** Section 1 used reading A's T3, which rests half on an in-sample stand-in, while section 3 used
  task h's payoffs. Both readings are now shown side by side. The gate estimate's T3 weights are flagged for the
  captain: they make T3's upper tail optimistic.
- **A leaderboard threshold.** Task a's z-rule and the "15-25 points" headline used a leaderboard figure as a
  threshold. Both are struck, and the distance is now stated against our own percentiles.

**Corrections made:**
- The depth bet is called a quantification of a known bet.
- The gate bar is quoted as 0.22-0.66, with the committed 0.48-0.52 on a real truth.
- The claim that s1 passes easily with matching depth is withdrawn: it also needs factor agreement.
- U2 is stated as 0.4-0.9.
- The exact up-kernel is relabelled a revival of a pruned line.
- The gene-set levers are gated on A15.
- E7.25 is removed as evidence.
- The embryo's idealised gain is corrected to +0.10 skill.
- The tail decoder is restated as +0.05 to +0.25.
- Citation fixes:
  - [73] for rank preservation, with [72] as the opposing case;
  - [25] for the intervention claim;
  - [78] dropped from the field-data sentence.
- The process slips are recorded in DATA_POLICY 2026-11-05 (k).
- The diagnostics are registered (TRIES N3-N5) and rerun stamped at 1edec56.
- The tails are measured on every upload file (N5).

**Not changed:** the eight notes themselves, which stay as audited. Where the report and a note differ, the report
says why.

## 9. Key findings

1. **Nothing found moves our total into its top-3 region.** The open levers add roughly +1 to +3 Overall in
   expectation, against the 12-17 points between reading A's mean and its 90th-99th percentiles.
2. **s1's T3 severity rests on two bets: a matching depth deficit and gene-specific agreement of its factor.** The
   depth step is 43% of s1's change; without matching depth, the gate bar is 0.22-0.66. The Gata4-free backup question
   turns on U2, U3 and the floor reading, all measurable or arguable before T3's first upload (section 3).
3. **The top-3 chance is a few percent at most.** Under the measurement-based T3 payoffs it is under 1% at our 95th
   percentile. It rests on s1's bets.
4. **The skill map's floor normalisation punishes boldness when the truth sits near the baseline.** That favours a
   near-baseline embryo backup, which also hedges a logged depth risk (section 6).
5. **Measurement closed three lines:** the T1 tail decoder (N3), the T3 and heart decoders (N5), and the T3 covariance
   programme (task b).

## 10. Limitations and future directions

**Limitations:**
- Payoffs come from toys, closed forms or judgement, except the three diagnostics.
- Several sources were read only at abstract or search-summary level: [4][66][69][70][75][76][77][81][82]; [30] was not
  read.
- The research could not use competitor material, the leaderboard, in-window literature or gene-level knockout
  literature, by design.
- The toy scripts live in the session scratchpad. Their sha256s are listed in `notes/toys.sha256`; they were not
  screened, so they are not committed.

**Future directions:** the decision tools of section 7, in deadline order:
1. D5/D6 and D1/D2, then the rank-aware pricing (sheet A, before 2026-11-27).
2. The embryo pseudo-truth table (sheet B).
3. The A15 reading (sheet D).

## 11. References

The full registry, with source type, accessibility, date, authority and task, is `notes/registry.md`: 131 entries and
128 distinct works ([3]=[68], [4]=[69] and [28]=[121] are the same works). The numbers cited above are its global
numbers. Main sources:
- [1] Steinwart & Ziegel, strictly proper kernel scores.
- [2] Gneiting, Making and evaluating point forecasts.
- [3] Witkowski et al., Incentive-compatible forecasting competitions.
- [4] Lichtendahl & Winkler 2007.
- [5] Stephens, ash.
- [6][8][9] the probability ranking principle.
- [16] Masnadi-Shirazi, kernel scoring rules.
- [17] CIPHER.
- [19][20] Kaneko; Sato et al.
- [22] Furusawa & Kaneko.
- [23] Gupta et al.
- [25] CausalBench.
- [27] Lehner.
- [28] Mitchel et al.
- [29] Systema.
- [32][33] MSigDB.
- [37][38] CollecTRI and DoRothEA.
- [39][41][42] JASPAR and HOCOMOCO.
- [44][45] ChIP-Atlas.
- [51][52] Gerard; Neufeld et al.
- [53] Sarkar & Stephens.
- [66] Gaba, Tsetlin & Winkler.
- [67] Ottaviani & Sorensen.
- [71] Aldous.
- [72] Roelofs et al.
- [73] Recht et al.
- [76][77] risk-taking in contests.
- [83][84][85] damped trend and Theta.
- [87] TIMO (method point only; in-window results not used).
- [89] TIGON.
- [90][91] Wasserstein splines.
- [101][106][107] cell-sorting physics.
- [102][103] kernel scores and the MMD test.
- [114] Memoli.
- [117] MERFISH with expansion.
- [118] Atta et al.
- [119] Bhuva et al.
- [120] Janssen et al.
- [121] Mitchel et al.
- [122] Baysor.
- [125][126] SoupX and DecontX.
- [127][128][129] dissociation and nuclei effects.
