---
task_id: a
role: Statistical decision theorist (proper scoring rules, kernel methods, ranking under uncertainty)
status: complete
sources_found: 16
---
## Sources
[1] Steinwart & Ziegel, Strictly proper kernel scores and characteristic kernels on compact spaces (Appl. Comput. Harmon. Anal.) | https://arxiv.org/pdf/1712.05279 | Source-Type: academic | Accessibility: public | As Of: 2021-03 | Authority: 9/10
[2] Gneiting, Making and Evaluating Point Forecasts (JASA) | https://ar5iv.arxiv.org/html/0912.0902 | Source-Type: academic | Accessibility: public | As Of: 2011-06 | Authority: 10/10
[3] Witkowski, Freeman, Vaughan, Pennock & Krause, Incentive-Compatible Forecasting Competitions (AAAI 2018) | https://object.cloud.sdsc.edu/v1/AUTH_da4962d3368042ac8337e2dfdd3e7bf3/ml-papers-txt/AAAI/2018/incentivecompatible_forecasting_competitions__5b0f617b.txt | Source-Type: academic | Accessibility: public | As Of: 2018-02 | Authority: 8/10
[4] Lichtendahl & Winkler, Probability Elicitation, Scoring Rules, and Competition Among Forecasters (Management Science 53(11)), abstract record | https://ideas.repec.org/a/inm/ormnsc/v53y2007i11p1745-1755.html | Source-Type: academic | Accessibility: semi-public | As Of: 2007-11 | Authority: 8/10
[5] Stephens, False discovery rates: a new deal (Biostatistics 18(2)) | https://pmc.ncbi.nlm.nih.gov/articles/PMC5379932 | Source-Type: academic | Accessibility: public | As Of: 2017-04 | Authority: 9/10
[6] Fairness in Ranking under Disparate Uncertainty (states the PRP top-k justification) | https://arxiv.org/pdf/2309.01610 | Source-Type: academic | Accessibility: public | As Of: 2023-09 | Authority: 6/10
[7] Wang & Zhu, Portfolio Theory of Information Retrieval (SIGIR 2009) | https://lintool.github.io/robust04-analysis-papers/p115-wang.pdf | Source-Type: academic | Accessibility: public | As Of: 2009-07 | Authority: 7/10
[8] York University EECS 6339 lecture notes, probabilistic retrieval (PRP minimises Bayes risk under 0/1 loss) | https://wiki.eecs.yorku.ca/course_archive/2013-14/W/6339/_media/11prob.pdf | Source-Type: other | Accessibility: public | As Of: 2014-01 | Authority: 5/10
[9] Gordon & Lenk, When is the probability ranking principle suboptimal? (catalogue record with abstract) | https://vufind-test.katalog.k.utb.cz/SummonRecord/FETCH-LOGICAL-c5111-8fab5ca1244d0cbc116178c13d083383eb3ee4cfb7a49012d31fcf3aa0702d683 | Source-Type: academic | Accessibility: semi-public | As Of: 1992-01 | Authority: 6/10
[10] scoringrules documentation, variogram score (vs_ensemble) | https://scoringrules.readthedocs.io/en/latest/generated/scoringrules.vs_ensemble.html | Source-Type: other | Accessibility: public | As Of: 2026-11 | Authority: 5/10
[11] Pic, Dombry, Naveau & Taillardat, Proper Scoring Rules for Multivariate Probabilistic Forecasts based on Aggregation and Transformation | https://arxiv.org/pdf/2407.00650 | Source-Type: academic | Accessibility: public | As Of: 2024-07 | Authority: 8/10
[12] Nowozin, Optimal Decisions from Probabilistic Models: the Intersection-over-Union Case (CVPR 2014) | https://www.microsoft.com/en-us/research/?p=166454 | Source-Type: academic | Accessibility: public | As Of: 2014-06 | Authority: 8/10
[13] Lichtendahl, Grushka-Cockayne & Winkler, Is It Better to Average Probabilities or Quantiles? (Management Science 59(7)) | https://pubsonline.informs.org/doi/fpi/10.1287/mnsc.1120.1667 | Source-Type: academic | Accessibility: semi-public | As Of: 2013-07 | Authority: 8/10
[14] Regression dilution (encyclopedia entry) | https://Www.Wikipedia.org/wiki/Disattenuation | Source-Type: other | Accessibility: public | As Of: 2026-11 | Authority: 4/10
[15] Ryvkin et al., Gambling in Contests (contest risk-taking literature, read at search-summary level) | https://drupal.uts.edu.au/sites/default/files/2022-06/Ryvkin%20et%20al%20Gambling_in_Contests.pdf | Source-Type: academic | Accessibility: public | As Of: 2022-06 | Authority: 5/10
[16] Masnadi-Shirazi, Strictly Proper Kernel Scoring Rules and Divergences with an Application to Kernel Two-Sample Hypothesis Testing | https://arxiv.org/pdf/1704.02578 | Source-Type: academic | Accessibility: public | As Of: 2017-04 | Authority: 6/10

## Findings
1. A kernel score's divergence is half the squared distance of the mean embeddings, ||Phi(P) - Phi(Q)||^2 = 2[ integral S_k(Q,x) dP - integral S_k(P,x) dP ], so expected squared MMD is minimised by the forecaster's own predictive distribution, strictly so exactly when the kernel is characteristic [1][16].
2. Adding a bounded kernel to a characteristic kernel keeps it characteristic (Lemma 3.3), so an average of RBF bandwidths such as veckit's is a strictly proper kernel score on its embedding, yet kernel metrics can still rate as close two distributions that are far apart in total variation (Theorem 3.7) [1].
3. The optimal point forecast is the Bayes rule of the scoring function actually used, and relative-error scores have Bayes rules shifted toward small values (beta-medians; the absolute percentage error "tends to support severe underforecasts"), so a forecast built for the wrong functional can be heavily penalised [2].
4. Ranking by calibrated marginal probability maximises the expected number of relevant items in every top-k prefix, and this optimality requires calibrated, certain and independently judged probabilities [6][8][9].
5. A mean-variance (portfolio) ranking that adds a tunable risk term to the probability of relevance beats the probability ranking principle most clearly on early-precision measures [7].
6. With heterogeneous standard errors, adaptive shrinkage re-orders effects relative to their raw estimates or z-scores, giving poorly measured observations a higher local false sign rate at the same p-value [5].
7. Winner-take-all selection by summed proper score is not incentive compatible: a forecaster raises her chance of ranking first by reporting more extreme probabilities, at a cost in expected score [3][4].
8. The variogram score is proper but not strictly proper, cannot see a bias shared by all components, and at order p = 0.5 is most sensitive to small deviations [10][11].
9. Decisions that maximise expected intersection-over-union must use the joint predictive model rather than per-element marginals, and the gain depends on calibrated probabilities [12].
10. Averaging quantiles gives the same mean as averaging probabilities with lower variance, a sharper combined forecast [13].

## First-principles analysis
Synthetic checks (numpy only; no data, no repo code, no score): /tmp/claude-1001/-mnt-md0-IITM-BackUp-Home-vijayavallabh/07d440b9-3ae8-4421-8a08-a046b42f307d/scratchpad/nonmimetic_research_2026-11-05/sims/sim1_mixture_vs_point.py, sim2b_de_score_prp.py, sim3_dedir_sev.py and sim4_offset_top3.py in the same folder. Every number marked "toy" comes from them and illustrates a mechanism, not a size at test.

### 0. Notation and the shared-noise structure
- R is the reference stage (T1 and heart E10.5, embryo E7.5, T3 WT E8.75) and T the target. Under the bake protocol `vec.evaluate` reproduces, the server draws R_w (10% of R), which serves as both the DE reference and the floor. It splits a 10% copy of T into the truth half T_h and the ceiling half T_c. Our file P (count rule) is scored whole, and each metric subsamples it itself.
- pb is the per-gene mean of log1p CP10k. dt = pb(T_h) - pb(R_w) = Delta + eta + eps, and dp = pb(P) - pb(R_w) = delta + xi + eps. Here eps = pb(R) - pb(R_w) is the scorer's reference-copy noise, which dp and dt share; eta is the truth-half noise, xi our own build noise and delta our intended change.
- Skill uses x = d / d_f, where d_f = |floor - ceiling| and d is our signed distance from the ceiling, and s = 1/(1 + x). The map is convex and decreasing for x > -1 and clipped at 1. By Jensen, for a fixed E[x], spread in x raises E[s]. The de_score and de_direction floors are exactly 0, because copy_last's dp is 0.
- A1 (assumption): the live server follows the bake protocol (all 66 validation anchors are reproduced exactly, METHOD_SUMMARY.md).

### 1. mmd_u
**Pinned formula.**
- PCA to 30 components is fitted on T_h, and up to 2000 cells are drawn from each side.
- gamma0 = 1 / median of the non-zero squared distances among the truth's drawn cells.
- k = mean over s in {1/4, 1/2, 1, 2, 4} of exp(-gamma0 s ||a - b||^2).
- The estimator is the unbiased U-statistic (diagonals dropped). Lower is better.

**(a) The Bayes act for E[MMD^2].** For a fixed kernel,
E_Pi ||mu_P - mu_Q||^2 = ||mu_P - E_Pi mu_Q||^2 + E_Pi ||mu_Q - E_Pi mu_Q||^2.
The minimiser is therefore the mean embedding of the posterior predictive mixture P_bar = E_Pi[Q]. This is kernel-score propriety [1][16].

Three qualifications:
- The basis and bandwidth are fitted on T_h, so the result is exact only conditional on the kernel (A2: the truth-fitted kernel varies little across plausible targets).
- The U-statistic scores the population our cells sample. Duplicated rows add about dup_share (1 - E k)/N, and a herded, low-discrepancy set subtracts about (1 - E k)/(N - 1). Both are negligible at our counts, and the count is closed anyway.
- MMD can be blind to large changes in total variation (Theorem 3.7 [1]).

**(b) What the skill map changes.**
- *Discrete scenarios* (e.g. two compositions, prior 1/2), with r = MMD^2(Q1, Q2)/d_f:
  E[s](w) = 0.5/(1 + (1 - w)^2 r) + 0.5/(1 + w^2 r).
  - The equal mixture is a local maximum only when r < 4/3 (second derivative proportional to 3r/4 - 1), and it beats a pure point bet only when r < 2.
  - At r = 3 the best act is w = 0.07, scoring 0.632 against 0.571 for the mixture (toy, sim1 a).
- *Continuous uncertainty about how far cells moved.* The toy uses truth N(theta v, I) in 30 PCs with theta ~ N(m, tau^2), the floor at theta = 0, and veckit's bandwidth mix.
  - The exact posterior predictive (variance inflated along v) lowers E[MMD^2] by only 2-5% and never raises E[s].
  - A point at 0.8 of the expected move beats a point at the expected move: 0.800 against 0.784 at tau/m = 0.5, and 0.627 against 0.620 at tau/m = 1 (toy, sim1 b).
  - The reason: d_f(theta) is small exactly when the true move is small, so the map acts as a relative error weighted by 1/d_f(theta). Its Bayes act is shaded toward the reference, the same mechanism by which relative-error scores favour low forecasts [2].
- *Better of two uploads.* Two points at m +/- (0.5 to 0.8) tau raise E[max s] by 0.08 to 0.11 over one point (0.787 to 0.875 at tau/m = 0.5; 0.629 to 0.742 at tau/m = 1). A pair of mixtures gains nothing (toy, sim1 c).
  - For equal-mean uploads with correlation r_12, E[max] - mu = sd sqrt((1 - r_12)/pi).

**(c) Our setting.**
- On T1 the mixture over composition scenarios is just the posterior-mean share vector. With two uploads the right object is a two-point quantisation, which the frozen s1/s2 pair already is (the committed oracle-base grid puts it within about 1.5 points of the best pair). T3's composition is closed.
- Hedging gain, summarised:
  - within one upload: about 0 for continuous shift uncertainty, and +0.04 to +0.05 skill for close discrete scenarios (r <= 1; about +1.2 to +1.5 T1 points at mmd_u's weight 0.30), turning negative when r > 2;
  - across two uploads: +0.08 to +0.11 skill (about +2.5 to +3.3 T1 points), but that pair is fixed by docs/S2_RULES.md.
- Not worth building: within-upload hedges over shift sizes (for example drawing Feng's fraction from a distribution).
- What survives is mild shading of uncertain moves toward the reference (lever L5), limited by de_direction (section 3).

### 2. de_score
**Pinned formula.**
- Truth sets: U = {dt >= 0.25} and D = {dt <= -0.25}; the code's own comment says the BH step is inert. n_up = |U| and n_dn = |D|.
- The prediction contributes its top n_up and bottom n_dn genes by dp, and raw = (|top ∩ U| + |bottom ∩ D|)/(n_up + n_dn).
- chance = max of raw for the rankings by +pb(R_w) and by -pb(R_w), so the orientation is chosen in hindsight.
- score = (raw - chance)/(1 - chance), and 0 if std(dp) < 0.01 std(dt).
- Skill = c/(2c - v), where c is the ceiling half's score and v ours; the slope at v = 0 is 1/(4c).

**(a) What our act moves.** chance and c depend only on the truth, so our act moves raw alone. The factor 1/(1 - chance) gives each hit more weight when the null is strong (1.6-1.7x on E9.5 -> E10.5).

**(b) The PRP ranking.** Given (n_up, n_dn), linearity gives
E[raw] = (sum over top of p^U_g + sum over bottom of p^D_g)/(n_up + n_dn),
with p^U = P(dt >= 0.25) and p^D = P(dt <= -0.25). One ranking by q = p^U - p^D serves every cutoff, so unknown n_up and n_dn do not matter. This is the probability ranking principle, optimal for every prefix when the probabilities are calibrated and certain [6][8][9].

**(c) Down-feasibility.** p^D_g = 0 whenever pb(R_w)_g < 0.25, because pb(T_h) >= 0.
- This rules out most of T1's 32,285 genes, and about 310-320 of a 500-gene panel (the digest counts 178-188 feasible).
- A decrease of such a gene gives a small negative dp. It fills a bottom slot whenever n_dn exceeds the number of feasible genes with larger predicted decreases. Every such slot is a guaranteed miss.
- At T1 test, scope removal makes this likely. Markers of the minor non-chamber tissues have pb(E10.5) below 0.25 (a tissue at a few % of cells), fall toward 0 in s1, and so sit in s1's bottom tail. They fall in the truth too, but by less than 0.25, so they can never be hits.

**(d) Heteroscedastic uncertainty.** With Delta_g ~ N(m_g, s_g^2), p^U = Phi((m - 0.25)/s). Ordering by m puts noisy genes with large |m| into the tails. Reliability re-ordering, as in ash, removes them [5].

**(e) Released evidence** (committed `experiments/harness/results/t1_de_v1/de_audit.txt`, released stages only).
- On E9.5 -> E10.5 chance is 0.39-0.41 (up/down 69/382 whole, 132/426 chamber): the changed genes are mostly highly expressed genes going down.
- On E8.5 -> E9.5 chance is 0.02 (31/30).
- A copy names 0.11 and B0 0.20-0.28, so both score below 0 on E9.5 -> E10.5.

**(f) Toy sizes** (sim2b: 3,000 genes; de_score units):

| Variant | Weak-null world | Strong-null world |
|---|---|---|
| Feasibility lift alone | +0.02 | +0.13 |
| Heteroscedastic PRP | +0.21 | +0.31 |

- An expression prior P(down | pb_ref, z-bin), calibrated on an independent transition of the *same* regime, adds +0.19 to +0.22 beyond PRP in the strong-null world.
- Calibrated on the *other* regime, it costs 0.02-0.05.
- The payoff is asymmetric: a large gain if the regime transfers, a small loss if it does not.

**(g) Realisation in count space without hurting the other metrics.** dp is our cells' mean log1p minus pb(R_w), so an order can only be imposed by per-gene count edits. Use shrink-only edits:
1. Compress the predicted decreases of infeasible genes, keeping their order, into a band just below 0 and above the feasible decreasers' tail.
2. Shrink unreliable genes' log-factors by their reliability.

Both move genes back toward the reference.
- On T1 these are low-expressed genes with negligible loadings on the truth's top 30 PCs, so mmd_u should not move.
- Their variogram entries return to the reference's sparsity.
- de_direction keeps their sign and, if anything, ranks them more plausibly, since their true change is bounded by pb_ref < 0.25.

Two cases are not free:
- Scope-driven drops come from absent cells, not from factors. Lifting them needs a low-level top-up of those genes in the retained cells, and the variogram is sensitive near zero (sqrt is steep at 0). That cost must be measured, or those genes left alone.
- Promoting genes into a tail (the expression prior) needs amplified decreases of highly expressed genes, which load on the PCs. That costs mmd_u and the variogram if wrong.

**(h) Notes on our entries.**
- CP10k closure already gives Feng-group cells a broad small decrease of their other highly expressed genes (log of the per-cell constant; N2 found 0.91-0.96 for momentum), a mild version of the strong-null pattern.
- With the shared eps, keep delta at the posterior-mean scale (section 3).

### 3. de_direction
**Pinned formula.** rp, rt and rr are the ranks of dp, dt and pb(R_w) over all genes, and
r = (r_pt - r_pr r_tr)/sqrt((1 - r_pr^2)(1 - r_tr^2)).
It is 0 if std(dp) < 0.01 std(dt) or |r_pr| >= 1 - 1e-6, and clipped to [-1, 1].

**What it measures.** Agreement in the order of change that the reference level cannot explain. It is invariant to every strictly increasing map of dp, so magnitude-free, but not to delta's scale relative to eps, which the scorer adds after our choice.

**The optimal act.** Without ties, rank vectors have fixed norms, so E[<rp_c, rt_c>] is linear in rp. By the rearrangement inequality, the best order is by E[rank(dt)_g]. In the partial version it is the expected residual rank given rr; any rr component in our ranks is partialled out, neither rewarded nor penalised beyond the collinearity guard.

**Shared noise.** For Pearson correlation, corr(lambda m + eps + xi, Delta + eta + eps) is maximised at
lambda* = 1 + V_xi/V_eps
(A3: m is the posterior mean and the errors are independent).
- The rank-partial version follows it with a flat optimum: the toy's best lambda is 1.0, 1.0, 1.5 and 2.0 for V_xi/V_eps = 0, 0, 0.33 and 1 (sim3 part 1).
- Shrinking to lambda = 0.5 costs 0.02-0.04 in r, and to 0.25 costs 0.12-0.15.
- Our builds draw the count rule from the whole pool with balanced repeats, so V_xi << V_eps (the scorer's 10% copy) and lambda* is about 1.
- The rule: realise the posterior-mean scale. Inflation buys nothing, and strong shading costs. This also explains why herding (T5, T6) found nothing.
- A copy's +0.10 to +0.22 on the released T1 backtests comes from this eps channel and is never evidence of direction.

### 4. severity_slope (T3)
**Pinned formula.** Over the truth's DE genes:
beta = dt.dp/||dt||^2, and r2 = 1 - ||dp - beta dt||^2/||dp - mean(dp)||^2.
- The value is log beta, scored as |log beta| against the floor |log 1e-3| = 6.91 (wt_identity) and the ceiling's |log beta_c|.
- The worst value applies if beta <= 0, r2 <= 0.01, or std(dp) < 0.01 std(dt).
- When it passes, skill is about 6.91/(6.91 + |log beta|).

**(a) Exact gate geometry.** With uncentred cosines over the DE genes,
r2 > 0.01 <=> cos^2(dp, dt) > 0.01 + 0.99 cos^2(dp, 1).
(The reason: the fit has no intercept, but r2's denominator is centred.)
Write dp = a 1 + u and dt = b 1 + v, with u, v orthogonal to 1.
- The offset that maximises the margin is a* = b (u.v)/||v||^2: the truth's common shift, scaled by the centred regression coefficient.
- At a*, the gate is cos^2(u, v) > 0.01, i.e. a centred |cos| above 0.1.
- At a = 0, it needs cos^2(u, v) > 0.01 (1 + n b^2/||v||^2).
- With u.v = 0 it fails, whatever the offsets.

So a predicted common shift (a depth ratio, a one-signed factor) without proportional gene-specific agreement sends severity to the floor.

**(b) Magnitude.** With dp = kappa phi and rho = cos(phi, Delta), beta = kappa rho ||phi||/||Delta||.
- At the posterior-mean scale (OLS of Delta on phi), beta = rho^2. The scorer regresses predicted on observed, the reverse of a calibration regression, so shrinkage is punished twice [14].
- At the correct norm, beta = rho; beta = 1 needs a further factor of 1/rho.
- For rho = 0.3, skill is 0.74 at the posterior mean, 0.85 at the correct norm and 1.0 at beta = 1.

**(c) Under uncertainty about rho.** s(l) = 6.91/(6.91 + |l|) is nearly linear in |l|, so E[s] peaks near the posterior median of log beta = 0 (the absolute-error Bayes rule [2]). The gate and the mutant_state group (0.20, which prefers kappa small) bound it.

**(d) Shared eps.** A near-WT prediction built from the whole WT passes the gate in 51-97% of toy draws (skill 0.59-0.77) as |eps|/|Delta| on DE genes goes from 0.15 to 0.4 (sim3 part 2). An uninformative added response dilutes that channel, and a very large kappa removes it.

**(e) Toy with a common shift** (sim4 part A: 60 DE genes, common-shift share 0.3-0.7 of the response energy).

| How the predictor handles the common shift | Small shared noise, centred agreement rho <= 0.2 | Large shared noise (0.4) |
|---|---|---|
| Carries the full true shift | skill 0.50, gate passed in 0-5% of draws | 0.81-0.95 at the posterior-mean scale; 0.50-0.60 when inflated x3 |
| Shift scaled by rho | 0.63-0.88, passes 70-100% | — |
| No shift | 0.60-0.75 | — |

- Across all 12 settings, the rho-scaled shift was best in 7 and never scored below 0.67.
- The full shift sat at or below 0.53 in 11 of its 24 cells.

**(f) Rule for global offsets.** Shrink any common component of the predicted response over plausible DE genes to rho_hat times its size, where rho_hat is the expected centred agreement.
- A5 (assumption, to check target-free): s1's Gata4 per-type depth is mostly such a component. In MERFISH most panel genes are sparse, and thinning to depth ratio p lowers their mean log1p roughly in proportion to pb, one-signed.
- Its de_score and de_direction signal is largely an expression-level pattern, which the null and the partialling neutralise, so its value sits in severity and mutant_state.

### 5. variogram
**Pinned formula.** Up to 1500 cells per side; 20,000 random gene pairs i != j, drawn after the cell subsamples from the same generator. gamma(i,j) = mean over cells of |x_i - x_j|^0.5, and VS = mean over pairs of (gamma_P - gamma_T)^2.

**What it measures.** Per pair, a 0.5-power dispersion of the two genes' joint log1p distribution.
- For pairs with a silent gene, it reduces to E[x^0.5] of the other gene. That is a marginal dominated by detection, since sqrt is steep at 0.
- It is proper but not strictly proper and blind to biases shared by all components [10][11]. A constant added to every gene of a cell cancels.

**The Bayes act.** Each gamma is linear in the cell distribution, so the posterior mean of gamma_T is the gamma of the posterior predictive mixture. The squared error adds our estimator's variance, which the scorer's 1500-cell subsample fixes, so there is no lever at a closed count.

**Edits.**
- Preserve: resampled real cells; count-space thinning or own-count top-up (zeros stay zero); per-cell uniform shifts.
- Damage:
  - per-gene additive log shifts on every cell (they fill zeros);
  - smoothing, averaging or blending cells (shrinks |x_i - x_j|);
  - depth changes that move detection away from the target's;
  - own-count top-up, which under-realises a rising gene's detection.

The skill-map effects of section 1 apply unchanged.

### 6. T2 under fixed count, size and scope
Board = 25[0.5 s(de_score) + 0.5 s(de_direction)] + 25[0.6 s(mmd_u) + 0.4 s(variogram)] + 25 mean(s(d2), s(dice), s(|scale|)) g + 25 s(nbhd), with g = min(1, s(nbhd)/0.5).
- **Expression change.** As on T1, on 500 genes with about 180-190 down-feasible; the smaller reference copy makes the eps channel stronger.
  - The embryo mixture realises dp = 0.5 (pb(E8.0) - pb(E7.5)) + eps exactly, with no factor to re-order without editing cells.
  - The heart's count-space B0 can take the tail decoder (L1).
- **Cell state.** As in section 1. Both embryo uploads share one expression mixture, so better-of-two hedges only the geometry (a frozen S2 rule).
- **Shape.**
  - d2 is the W1 distance of the RMS-normalised pair-distance distributions, i.e. the integrated absolute difference of their quantile functions. Its Bayes act is the pointwise median of the target's quantile function (Vincentisation logic [13]).
  - Dice: expected-Dice decisions include voxels in order of occupancy probability and need calibrated probabilities [12].
  - Scale is closed.
- **Local organisation.** The Bayes act is a posterior-predictive mixture of 15-NN profile distributions: real tissue with its own neighbours. The gate bites only below the floor.
- Nothing new survives the closures here except the heart tail decoder (in L1).

### 7. A top-3 objective
**The calculus.** P(top-3) is about Phi((mu - t)/sigma) (A4: normal, roughly independent board errors). With mu about 184, t about 200-209 and sigma about 10-14, z = (t - mu)/sigma is about 1.3-2.1 and P is about 2-9% (sim4 part B).
- An independent lever with mean g and sd s is worth g + z s^2/(2 sigma) mean points.
- With z/sigma about 0.1-0.2, s = 4 adds 0.9-1.7 and s = 6 adds 1.9-3.8.
- A zero-mean bet with sd 4 moves P about as much as a sure +1.

**When to trade mean for variance.** All four must hold:
1. The variance comes from genuine uncertainty about the target, not from added noise. Noise raises E[d] and lowers E[s].
2. The payoff is skewed: large in a plausible regime, small otherwise.
3. The bet is independent of our other bets or positively correlated with them.
4. It touches no frozen rule without the captain.

Competitive-forecasting theory agrees: winner-take-all rewards exaggerated reports [3][4]. Contests show trailing players take more risk, with mixed evidence that it pays [15].

**Within a board.** Better-of-two already favours negatively correlated pairs: E[max] - mu = sd sqrt((1 - r_12)/pi). A hedge pair doubles the chance of clearing a high bar compared with a duplicate.

**Limits.** The adoption rule is a mean rule and cannot value variance, so variance-seeking bets are the captain's blind bets. Uncertainty about t and about our own mu widens the effective sigma and lowers z, so do not over-weight variance.

## Candidate levers
Ranked by expected value to a top-3 objective (the mean gain plus the variance value of section 7).

**L1. T1 tail decoder: down-feasibility and reliability ordering of dp's tails (a realisation step)**
- **Board and group:** T1:test, de_recovery (de_score, 0.25), with an effect on direction (0.25). The same step serves heart expression_change (de_score 0.125) and, by the captain's ruling only, T3 response_identity.
- **Mechanism (section 2 c, g):** compress the predicted decreases of genes with pb(E10.5) < 0.25 + margin into a band just below 0, keeping their order. The margin is 2 sd of a 10% copy's pb, computed from E10.5's own cells. Optionally shrink log-factors whose sampling SE is large.
  - Factor-driven drops are compressed for free.
  - Scope-driven drops come from absent cells. They either get a measured low-level top-up in retained cells or are left alone.
  - Composition, depth, counts and Feng's direction are untouched.
- **Expected effect:** the toy gives +0.02 (weak null) to +0.13 (strong null) de_score units from feasibility alone. At ceiling c of about 0.7, that is +0.2 to +1.2 T1 points.
  - If the diagnostic finds 30-50% of the bottom tail infeasible, a hit rate of 0.2-0.4 on the replacement genes gives up to about +2.
  - de_direction moves between -0.2 and +0.2 points, and mmd_u should not move.
  - The gain is zero if no infeasible gene sits in the tail.
- **Legality:**
  - Changes expression realisation only. No composition, depth, count, size or scope change, no cross-modal input, no external data and no gene list; it reads only the build's own inputs.
  - Each scored run is a registered try.
  - M1 is the captain's ruling, so the step reaches s1, and s2 "as s1", only by his ruling or as a core builder step. Under S2_RULES a builder change for every entry carries over.
- **Validation route:**
  1. A target-free sizing diagnostic, not a try. For the built s1 against E10.5, and for B0's backtest builds against E9.5, count the genes with pb(reference) < 0.25 among the bottom K of pb(build) - pb(reference), K = 50, 100, 200, 400, 800. Split them by cause, scope against Feng or momentum, using the scope-only build. Prune if the count is near zero.
  2. A registered try on the counting backtests: B0 against B0 + decoder on T1/E9.5->E10.5 and T1/E9.5->E10.5:chamber (both counted, 10 seeds). Predicted gain on de_score, with de_direction within +/-0.2 and the distribution metrics flat. For the heart, the same comparison on both heart backtests.
- **Effort:** 1-2 days.
- **Main risk:**
  - The two counted T1 tasks score the same predictions against two cuts of one library, so they are nearly one result.
  - A pass on B0 may not carry to M1, whose tails come from scope removal and Feng.
  - Top-up for scope-driven genes may cost variogram.

**L2. T3 severity: offset discipline (remove or rho-scale the predicted common shift; s1's depth step first)**
- **Board and group:** T3:bcatenin, response_magnitude (severity_slope, 0.25), with side effects on mutant_state (0.20).
- **Mechanism (section 4 a, e, f):** a common shift over the DE genes passes the gate only in proportion to the centred agreement. s1 applies Gata4's per-type depth, which is mostly a one-signed shift in log1p (A5).
  - Knockout depth ratios flip between released folds.
  - K2a, without depth, already sat at the r2 gate in 9 of 10 seeds.
  - So under a sceptical prior the posterior-mean offset is near 0.
- **Expected effect:** conditional and high-stakes.
  - When the factor's centred agreement at test is about 0.1-0.4 and depth does not transfer, severity skill rises from about 0.50 (gated) to 0.60-0.75: +2.5 to +6 T3 points.
  - When bcatenin carries a depth shift matching Gata4's, removal loses 1-2.5 points, plus any mmd_u gain depth brings.
  - When the factor's agreement is about 0, both variants are gated unless the shared-noise channel carries the gate.
  - In expectation, roughly +0.5 to +2 T3 points if the captain's P(Gata4's depth transfers) is 0.5 or below, and negative above about 0.6. The sd is about 2-3 points, which adds roughly +0.2 to +0.8 of variance value.
- **Legality:**
  - Changes the output of frozen T3 rules, so it needs the captain's ruling; rules.yaml is not edited.
  - No composition, count, size or scope change, and no gene list, knockout literature, external data or cross-modal input.
- **Validation route:**
  1. A target-free diagnostic, not a try. For each built T3 file, decompose dp_proxy = pb(file) - pb(WT E8.75) into common and centred parts over plausible DE-candidate sets: feasible genes, and genes above the 80th percentile of |dp|. Report cos^2(dp, 1) and the required cos^2(dp, dt) = 0.01 + 0.99 cos^2(dp, 1), with and without the depth step.
  2. No counting fold can judge it: the test shares WT E8.75 with Gata4, and fold-spending diagnostics were pruned. Captain's blind bet only.
- **Effort:** 0.5 day for the diagnostic, 1 day to rebuild if ruled.
- **Main risk:** depth may be the dominant, correctly signed and gene-patterned part of bcatenin's response (same lab, stage, driver and WT reference). Then removing it costs severity and mmd_u.

**L3. T1 expression-level prior in the down tail (a blind bet)**
- **Board and group:** T1:test, de_recovery (0.25).
- **Mechanism (section 2 f):** PRP with P(down | pb_ref, predicted change), calibrated on both released T1 transitions by a rule committed before any run. E8.5->E9.5 is the weak-null regime and E9.5->E10.5 the strong one. It is realised by small decreases of highly expressed genes the model calls unchanged.
- **Expected effect:** if the test regime matches E9.5->E10.5, the toy gives +0.19 to +0.22 de_score units, about +1.5 to +2 T1 points. If not, -0.02 to -0.05, about -0.2 to -0.45, plus an unknown mmd_u and variogram cost from thinning highly expressed genes. The EV is positive if P(strong null at test) is above roughly 0.2-0.3, with an sd of about 1-1.5 points.
- **Legality:**
  - Uses released and unlocked data only, writes no gene list, and its rule is committed before any run.
  - Changes M1's expression, so it needs the captain.
- **Validation route:** none that counts. On the counting backtests the only legal calibration source is E8.5->E9.5, the other regime, so a backtest says nothing about its test value. Captain's blind bet only.
- **Effort:** 2 days.
- **Main risk:** E10.5 -> E12.5, heart to heart, may resemble E8.5 -> E9.5. The amplified decreases then cost the distributional terms.

**L4. Lever selection for a top-3 objective (a meta lever)**
- **Board and group:** all boards; no build change.
- **Mechanism:** value each blind bet at g + z s^2/(2 sigma), with z about 1.3-2.1 and sigma about 10-14. Prefer skewed bets that hinge on a regime (L2, L3) over symmetric noise, and keep s1/s2 negatively correlated.
- **Expected effect:** changes which blind bets the captain takes; worth about +1 to +2 mean-equivalent points on a top-3 objective.
- **Legality:** arithmetic on our own committed estimates (docs/PLAN.md section 8); no score inference.
- **Validation route:** target-free.
- **Effort:** 0.5 day.
- **Main risk:** mu, sigma and the threshold are themselves uncertain, so variance could be overvalued if we are closer than we think.

**L5. Shade uncertain moves toward the reference (the skill map's relative error)**
- **Board and group:** T1 distribution and covariance (Feng's fraction); heart cell_state (B0's damping).
- **Mechanism (section 1 b):** d_f(theta) shrinks with the true move, so the E[skill]-optimal point is shaded by roughly 1/(1 + c (tau/m)^2).
- **Expected effect:** the toy gives +0.016 skill at tau/m = 0.5, about +0.3 to +0.5 T1 points and under 0.3 heart board points, partly offset by de_direction below lambda = 1.
- **Legality:** T1's direction is a captain's ruling, and no T1 backtest can judge Feng's fraction (the window). Heart damping is a parameter of the shared B0.
- **Validation route:** the heart counting backtests could judge it, but the expected gain is under 0.3, so rule 5 prunes it. On T1: captain's blind bet only.
- **Effort:** 0.5 day.
- **Main risk:** under-moving when the true move is large.

**L6. T3 response scale (factor exponent) chosen by expected skill (closed as a try)**
- **Board and group:** T3, all four groups.
- **Mechanism:**
  - Severity wants beta near 1: kappa about rho^-2 times the posterior-mean scale when there is no common shift.
  - mutant_state wants kappa small.
  - The de metrics are scale-free except through the eps channel.
  - The gate wants offset discipline (L2).
- **Status:** the exponent line was pruned ("its rho_F cannot be identified on released data"), so this gives the captain the structure only. Closed; captain's blind bet only.
- **Effort:** none now.
- **Main risk:** no estimate of rho exists before test.

**Closed, or not worth building (for the record)**
- Within-upload mixtures over shift sizes or fractions: they cut E[MMD^2] by 2-5% and never raise E[skill] in the toy.
- Herded or low-discrepancy cell sets: the U-statistic moves by about 1/N, and the count is closed.
- The variogram's estimator variance: set by the scorer's subsample, and the count is closed.
- A second embryo expression hedge in s2: the S2 rule is frozen, and the embryo backtests cannot count for a mixing-fraction change.
- T1 depth, T3 composition, and cell numbers, sizes and scope everywhere: closed.

## Deep Read Notes
- [1] Steinwart & Ziegel, read via the ar5iv HTML rendering of arXiv 1712.05279, first 100k of 174k characters. Definition 1.1: S_k(P,x) = -int k(w,x) dP(w) + 1/2 int int k dP dP, negatively oriented. Theorem 1.1: ||Phi(P) - Phi(Q)||^2 = 2(int S_k(Q,x) dP - int S_k(P,x) dP); proper on M_1^k, strictly proper iff the embedding is injective (characteristic, for bounded kernels). Lemma 3.3: adding a kernel keeps injectivity. Theorem 3.7 and Corollary 3.8: distributions at TV distance 2 can be within epsilon in RKHS norm. For us: the posterior-mixture Bayes act, and MMD's blind spots.
- [2] Gneiting, ar5iv 0912.0902, first 100k of 133k characters. The Bayes rule x_hat = argmin_x E_F S(x, Y); consistency and elicitability; Bregman functions for the mean, GPL scores for quantiles (Theorem 3.3); S_beta(x, y) = |1 - (y/x)^beta| has the beta-median (density proportional to y^beta f(y)) as Bayes rule, and APE favours underforecasts; Osband's revelation principle (Theorem 2.6). For us: a normaliser that depends on the truth makes the board map a relative error, whose Bayes act shades toward small moves, and |log beta| points to the median of log beta.
- [3] Witkowski et al., full text. Choosing the top summed quadratic score is not truthful. Example 1: all believe 0.9, and a 0.95 report wins whenever the event occurs. The ELF event-lottery mechanism is strictly truthful for independent events, and correlated beliefs break it (Example 4). For us: under a top-k objective exaggeration is rational, and our four boards share latent scenarios.
- [5] Stephens, full body text (100k of 126k characters). The lfsr is the probability of a sign error, the minimum of the two one-sided posterior tails, and is more robust than lfdr. With heterogeneous SEs, ash re-orders significance relative to p or z, and poorly measured effects get a higher lfsr at equal p. A unimodal prior shrinks toward the mode, and lfsr is conservative in simulations. For us: reliability-ordered tails and shrink-only realisation.

## Gaps
- **Searched for, not found:**
  1. A treatment of skill maps of the form d_f/(d_f + d) with a d_f that depends on the truth. The relative-error and beta-median reading is my inference from [2].
  2. Bayes acts for rank-partial (Spearman) targets; I derived them from the rearrangement inequality.
  3. The expected-overlap-optimal ranking when the truth sets the cutoffs and the score is measured against a null chosen in hindsight.
  4. Calibration literature on predicted-on-observed ("reverse") regression slopes, beyond generic regression dilution [14].
  5. An expected-Dice decision rule, as distinct from IoU [12].
  6. The energy-score mixture result stated as a theorem; it follows from propriety, but no source states it directly.
- **Not read in full:**
  - Lichtendahl & Winkler 2007 (abstract only).
  - The contest literature (search summary only).
  - The last 40% of [1], which holds the family-specific kernel results.
- **All sizes are toys.** The real values must come from the target-free diagnostics and registered tries named in L1 and L2.
- **Counter-claims against L1 (the top lever):**
  - At test, n_dn may be small enough that the bottom tail never reaches infeasible genes, making the gain zero.
  - Most infeasible tail genes may be scope-driven. Their lift needs a top-up whose variogram cost (sqrt is steep at 0, over many sparse genes) can cancel the de_score gain.
  - The counted T1 pair is nearly one result, and B0's tails differ from M1's.
  - The compressed decreases are real small drops, so de_direction may lose a little.
- **Counter-claim against L2:**
  - bcatenin and Gata4 share lab, stage, driver and WT reference. Library depth may be a property of how knockout embryos were processed, so Gata4's depth may transfer with the right sign and gene pattern.
  - The cross-knockout flip comes from a different stage and WT, and is weak evidence about bcatenin.

## Compliance log
- Barred-class items met: none.
- Note, not a barred class: two WebFetch calls on arXiv PDFs returned binary, and the harness saved copies under the account directory. I did not open them (AGENTS.md); I read the ar5iv HTML renderings instead.
- No hidden-target property, knockout gene-level result, protected-window data, leaderboard or partner material, scorer docstring number, or test-board file was opened, searched for or recorded.
## END
