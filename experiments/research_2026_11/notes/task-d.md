---
task_id: d
role: Computational statistician (survey calibration, entropy balancing, I-projection, count models)
status: complete
sources_found: 15
---
## Sources
[1] Data-based RNA-seq simulations by binomial thinning (Gerard, BMC Bioinformatics 21:206; doi 10.1186/s12859-020-3450-9) | https://pmc.ncbi.nlm.nih.gov/articles/PMC7245910/ | Source-Type: academic | Accessibility: public | As Of: 2020-05 | Authority: 8/10
[2] Data thinning for convolution-closed distributions (Neufeld, Dharamshi, Gao, Witten; arXiv 2301.07276, JMLR 25) | https://arxiv.org/abs/2301.07276 | Source-Type: academic | Accessibility: public | As Of: 2023-11 | Authority: 9/10
[3] Separating measurement and expression models clarifies confusion in single-cell RNA sequencing analysis (Sarkar, Stephens; Nat Genet; doi 10.1038/s41588-021-00873-4; full text read through PubMed/PMC) | https://pmc.ncbi.nlm.nih.gov/articles/PMC8370014 | Source-Type: academic | Accessibility: public | As Of: 2021-05 | Authority: 9/10
[4] Entropy Balancing for Causal Effects (Hainmueller; Political Analysis 20(1):25-46; abstract record) | https://econpapers.repec.org/RePEc:cup:polals:v:20:y:2012:i:01:p:25-46_01 | Source-Type: academic | Accessibility: semi-public | As Of: 2012-01 | Authority: 8/10
[5] ebal: Entropy reweighting to create balanced samples, reference manual v0.2.1 (Hainmueller) | https://cran.r-project.org/web/packages/ebal/refman/ebal.html | Source-Type: official | Accessibility: public | As Of: 2026-04 | Authority: 7/10
[6] I-divergence geometry of probability distributions and minimization problems (Csiszar; Ann. Probab. 3(1):146-158; bibliographic record, page did not render) | https://projecteuclid.org/euclid.aop/1176996454 | Source-Type: academic | Accessibility: semi-public | As Of: 1975-02 | Authority: 9/10
[7] Devaud and Tille (2019, TEST): optimisation treatment of Deville-Sarndal calibration, its distance functions and bounded weights (repository record) | https://libra.unine.ch/handle/20.500.14713/63825 | Source-Type: academic | Accessibility: semi-public | As Of: 2019-12 | Authority: 8/10
[8] calibev: calibration estimator and its variance (R package sampling; linear, raking, truncated and logit g-weights) | https://rdrr.io/cran/sampling/man/calibev.html | Source-Type: other | Accessibility: public | As Of: 2026-11 | Authority: 6/10
[9] Efficient balanced sampling: the cube method (Deville, Tille; Biometrika 91(4):893-912; abstract record) | https://econpapers.repec.org/RePEc:oup:biomet:v:91:y:2004:i:4:p:893-912 | Source-Type: academic | Accessibility: semi-public | As Of: 2004-12 | Authority: 9/10
[10] Computational Optimal Transport (Peyre, Cuturi; Remark 2.31, the Gaussian OT map and the Bures term) | https://arxiv.org/pdf/1803.00567 | Source-Type: academic | Accessibility: public | As Of: 2018-03 | Authority: 9/10
[11] The Schrodinger bridge between Gaussian measures has a closed form (Bunne, Hsieh, Cuturi, Krause) | https://arxiv.org/pdf/2202.05722v2.pdf | Source-Type: academic | Accessibility: public | As Of: 2022-02 | Authority: 7/10
[12] An interpolating distance between optimal transport and Fisher-Rao (Chizat, Schmitzer, Peyre, Vialard) | https://arxiv.org/pdf/1506.06430v2.pdf | Source-Type: academic | Accessibility: public | As Of: 2015-06 | Authority: 8/10
[13] A benchmark study of simulation methods for single-cell RNA sequencing data (Cao, Yang, Yang; Nat Commun) | https://www.nature.com/articles/s41467-021-27130-w.pdf | Source-Type: academic | Accessibility: public | As Of: 2021-11 | Authority: 8/10
[14] Genomic encoding of transcriptional burst kinetics (Larsson et al.; Nature 565:251-254; abstract record) | https://pubmed.ncbi.nlm.nih.gov/30602787/ | Source-Type: academic | Accessibility: semi-public | As Of: 2019-01 | Authority: 9/10
[15] Droplet scRNA-seq is not zero-inflated (Svensson; bioRxiv 582064, preprint of the 2020 Nat Biotechnol note) | https://www.biorxiv.org/content/10.1101/582064v1 | Source-Type: academic | Accessibility: public | As Of: 2019-03 | Authority: 7/10

## Findings
1. Binomial thinning is exact for Poisson counts (y ~ Poisson(a), y'|y ~ Bin(y, b) gives y' ~ Poisson(ab)), so thinning real counts adds a known log-scale signal while leaving the data's own unwanted variation and correlation untouched, but it realises positive effects only by thinning everything else to the row maximum, which spends depth [1].
2. Simulation by thinning real data reproduced the real data's scree plot and mean-variance trend where a fitted NB simulator added zeros and lost the small factors, and a 12-method benchmark found that several simulators struggled to keep real-data heterogeneity [1][13].
3. Every convolution-closed family splits exactly, Poisson by binomial or multinomial thinning and NB(r, p) by beta-binomial (Dirichlet-multinomial) thinning that scales the size r rather than the scale, which needs r known and leaves the parts dependent when r is wrong [2].
4. UMI counts fit a Poisson measurement model over an expression model (a Gamma often suffices, only 0.6-9% of genes favour a point-Gamma, measurement overdispersion is small), with zeros explained by low rates and Poisson sampling rather than a separate dropout process [3][15].
5. Cell-type-specific expression is shaped mainly by burst frequency (enhancers) rather than burst size (core promoters), so a real regulation can move a gene's zero pattern differently from the pure rate scaling that binomial thinning models [14][2].
6. Entropy balancing chooses maximum-entropy weights that meet prespecified moment conditions exactly (squared and interaction terms for variances and co-moments), and its software can trim the max-to-mean weight ratio under the same moments or flag the target as infeasible [4][5].
7. The I-projection of a distribution onto a convex set is unique and satisfies a Pythagorean relation that holds with equality for sets defined by linear moment constraints, the formal basis of minimum-KL reweighting [6].
8. Calibration estimators form one family indexed by a distance function (linear for GREG, KL for raking, truncated and logit for bounded weights), all consistent with GREG-like large-sample efficiency, though bounded versions can be infeasible [7][8].
9. The cube method draws a fixed-size sample without replacement, at given inclusion probabilities, whose Horvitz-Thompson totals approximately match known auxiliary totals, with stratification as a special case [9].
10. The OT map between Gaussians is affine (x -> m2 + A(x - m1), A from the Bures formula, W2^2 = |dm|^2 + a Bures covariance term), and unbalanced OT (Wasserstein-Fisher-Rao) interpolates between moving mass and re-weighting it [10][11][12].

## First-principles analysis

**0. Bottom line.** The realisation layer is not where the 15-25 missing points are.
- T3 already uses near-exact kernels: thinning below 1 and the kNN-smoothed top-up above 1.
- On T1 and the heart, the zero-preserving own-count top-up in `vec.core.apply_factors` does under-realise every increase. Fixing it is worth about +0.1 to +0.3 points on T1 at the size of Feng's direction, and its sign on the heart is unclear because it is mixed up with the closed depth lever.
- Composition changes realised by quotas cost almost nothing at the realisation step.
- NB resampling, latent Gaussian OT with nearest-cell decoding, and many-moment tilts are inferior or closed.

The earlier digest already showed that count-space edits beat dense or log-space edits. This note answers the next question: which count kernel to use, and why.

**1. What "costs nothing" means.** Measurement model (A1, assumed; [3]): c_ig | theta_i ~ Poisson(L_i theta_ig), where theta_i is the cell's latent expression vector under any expression model. The scorer sees x = log1p(a_i c), with a_i = 1e4/L_i (about 1 on T1, which has about 1e4 UMI per cell; 50-130 on the MERFISH boards).

Every ranked expression statistic is a functional of the law of x:
- pb_g = E x_g, which feeds de_score (ranks), de_direction (rank partial correlation) and severity (slope on the truth's DE genes);
- the variogram, through E|x_g - x_h|^0.5, mostly per-gene detection P(x > 0) and E[sqrt x | x > 0] (digest Part 1);
- mmd_u, through the joint law projected on the truth's 30 PCs.

So a realisation costs nothing exactly when the realised cells have the same law as an ideal population in which each kept cell's latent theta_i was moved by the intended change. Then every statistic equals its ideal value up to sampling noise. Any remaining loss belongs to the direction, not to the realisation.

**2. Count kernels for an expression change (factor f per type and gene).**

*Decrease (f < 1).* Binomial(c, f) given theta is Poisson(f L theta) for every expression model and jointly across genes [1]. Thinning is exact in law. Zero fraction becomes E[e^(-f L theta)], DE ranks and the severity slope are realised fully, and the covariance is diag(f) Cov diag(f), as the ideal change would give.

*Increase (f > 1), own-count top-up (T1 and T2 B0 now).* c' = c + Poisson((f-1)c). The mean is right, but the kernel is wrong in law:
- P(c' = 0) = P(c = 0), so detection never rises;
- Var(c' | theta) = L theta (f^2 + f - 1), against the exact f L theta;
- a gene off in a cell can never come on.

Fraction of the intended change it realises (own exact calculation, Poisson rates lambda = L theta):

| a (counts per cell) | lambda | mean log1p (pseudobulk) | E sqrt x (variogram building block) | detection |
|---|---|---|---|---|
| about 1 (T1-like) | 0.02-0.2 | 42-58% | 17-27% | 0% |
| about 1 (T1-like) | 1 | 64-65% | 36-42% | 0% |
| about 1 (T1-like) | 5 | 85-92% | 78-90% | 0% |
| about 100 (MERFISH-like) | 0.02-0.2 | 10-16% | 4-8% | 0% |
| about 100 (MERFISH-like) | 1 | 24-31% | 12-18% | 0% |
| about 100 (MERFISH-like) | 5 | 76-90% | 67-87% | 0% |

Factors 1.25-4 throughout. Decreases are realised at 100% and increases only partly, so dp is biased by sign. On MERFISH, B0's realised change is almost "decreases plus detection loss". This matches the ledger's note that B0's heart gain lies on the pure-thinning variogram-versus-variance line.

*Increase, exact kernel: posterior-predictive (empirical-Bayes) top-up.* Prior theta_ig ~ Gamma(alpha_g, rate alpha_g/mu_ig), where mu_ig is a kNN-neighbourhood mean proportion (cell itself excluded, shrunk to its type mean) and alpha_g is the gene's within-neighbourhood shape. Draw theta~ ~ Gamma(alpha_g + c, rate alpha_g/mu_ig + L_i), then set c' = c + Poisson((f-1) L_i theta~).
- If the prior is the true conditional law, (c, theta~) has the joint law of (c, theta), so c' ~ Poisson(f L theta) marginally. Up-moves become exact in law, as thinning is for down-moves.
- T3's frozen kernel, c + Poisson((f-1) n P_kNN), is the plug-in (posterior-mean, alpha -> infinity) version. Because e^(-(f-1)L theta) is convex and decreasing in theta, it slightly over-detects: it spreads new counts too evenly over cells with zero counts. The EB draw corrects this through the cell's own count.
- Cross-gene structure: the added part is drawn independently per gene given the neighbourhood, so the within-neighbourhood co-variation of the added counts is lost. This is second order when f - 1 is moderate.

*Burst-mode caveat (A2).* Thinning and the posterior top-up are exact under a scale-mode change (theta -> f theta, burst size). A frequency-mode change, which [14] says dominates, has exact kernels of its own [2]:
- decrease: beta-binomial thinning, Bin(c, B) with B ~ Beta(f r, (1-f) r);
- increase: an independent prior top-up, Poisson(L dtheta) with dtheta ~ Gamma((f-1) r, rate r/mu).
Both need r. In the toy (section 6) a frequency-mode truth realised with scale-mode kernels stays within about 1-2 toy points of the oracle, so the mode choice is second order.

*Closure (keep_depth, T1's rule).* At a fixed expected depth, the per-cell constant s_i = sum c / sum f c thins every gene of a cell with net increases. This is exact for a compositional (CP10k) change at fixed depth. The alternative, a floating depth, is the closed depth lever: a median library shift above 5% goes to the captain.

*Pseudocount consistency (new).* The directions are log((m_b + eps)/(m_a + eps)) at a time fraction, with eps = 0.5 CP10k. They therefore imply an absolute change of (f - 1)(m_a + eps), while multiplying the cell's own rate delivers only (f - 1) m_a. For genes with m_a << eps the delivered fraction is m_a/(m_a + eps), about 0 for genes that switch on. The consistent realisation adds Poisson((f - 1) eps L_i/1e4) to up-genes. That uses no external level, only the ratio and our own constant. Whether to do it is a direction decision, not a realisation one (lever 4).

*NB resampling.* c' ~ Poisson(L theta^), theta^ ~ Gamma(alpha, rate alpha/(f mu)). The marginals are right if the model is right, but the cell's own residual structure is replaced by gene-independent noise, so mmd_u and the variogram pay (toy: 3-9 toy points below EB).

*Calibrated (exact-total) variants.* Hypergeometric thinning, fixed-total multinomial top-up, or cube-method balancing of which cells receive the new counts [9] remove the edit's own Monte Carlo noise in per-gene totals. That noise is far below the scorer's 10% reference-copy noise (try T5, herding, found +/-0.03), so these are not worth a try.

**3. Reweighting realisations (composition): every real-cell method is a reweighting.**

*Which weights.* Among weights that meet moment targets E_w g(x) = t, the I-projection w_i proportional to exp(lambda' g(x_i)) is KL-closest to the reference and the least-informative choice [6]. Entropy balancing and raking are this choice with sample moments [4][7].

*Latent Gaussian OT with nearest-real-cell decoding is the same thing.* For x ~ N(m, S) pushed by x + d and decoded to the nearest reference cells, the decoded density ratio tends to N(x; m+d, S)/N(x; m, S) = exp(d' S^-1 (x-m) - d' S^-1 d/2). That is exactly the exponential tilt with lambda = S^-1 d, i.e. entropy balancing on the mean (own derivation; the affine map is from [10]). A Bures covariance change adds a quadratic term, i.e. balancing on squared terms. So OT decoding is a composition change in disguise:
- closed on T3 and T2;
- quota-only on T1;
- on T2 it would also swap expression away from positions, the failure that neighborhood_mmd punishes (digest Part 1).

*Price of a tilt.*
- Effective sample size: a tilt that shifts a Gaussian score by k SD has ESS/N = exp(-k^2), i.e. 0.78, 0.37 and 0.11 for k = 0.5, 1 and 1.5.
- Ceiling without duplicates: with inclusion probabilities at most 1 and sampling fraction p = n/N, the largest realisable within-stratum shift is phi(Phi^-1(1-p))/p SD, i.e. 0.64 SD at p = 0.6 and 1.16 SD at p = 0.3. This is consistent with the pruned within-type maturity tilt realising only 32-42% of its target at a draw fraction near 0.6.
- Many moments degenerate: a 30-moment exact balance in the toy had ESS 0.04-0.28 and toy scores 40-66, below even the own-count edit. The pruned stage-classifier tilt (ESS 0.12-0.28) is the same failure. A tilt can carry only the part of a change that lies along variation already present in the reference.

*Quota draws (T1's legal route).* Cells stay real, so per-gene marginals are exact mixtures, the variogram and de reflect the composition exactly (dp = sum_s dpi_s mu_s), and severity is set by the share change. Repeats when a quota exceeds its pool:
- leave pseudobulk and the variogram unbiased (they are per-cell averages);
- raise unbiased mmd_u by about q(1 - kbar), with q ~ (r-1)/(N-1) for r copies per cell among N submitted cells and kbar ~ 0.3-0.5 at the median-heuristic bandwidths;
- for r = 2 and N in the thousands that is of order 1e-4, against a T1 mmd_u floor-to-ceiling gap of order 1e-1 on the validation anchors, so negligible.
So composition choices should be judged on biology, not on worry about realisation. Quotas, not importance weights, also remove the dependence of the realised shares on n/N (the 2026-11-04 (f) reason).

*Hybrid.* The minimum-information way to reach a target that lies partly outside the reference's support is to re-weight for the in-support part and count-edit (thin or EB top-up) the residual: a Wasserstein-Fisher-Rao-type split [12]. It is legal on T1 only, through quotas.

**4. Board consequences (from the code: `vec.core.apply_factors` and `vec.t1_scope.external_shift` on T1 and in the heart's B0; `apply_gene_factor` with the kNN profile on T3).**
- **T3:** the kernel is already near-exact (toy: kNN within about 1 point of the oracle-theta kernel; EB the same). Composition is closed and the rules are frozen, so nothing changes. The caveat is the depth term, realised by thinning: tries E7 and H3 showed that thinned cells do not read like a shallower acquisition, which is a test risk of the frozen s1, not a lever.
- **T1:** Feng's increases, and the momentum before ruling R1, are under-realised (table above). The expected effect of a kernel swap scales like (log f)^2.
  - The skill map gives dS ~ 100 w s^2 dVS/d_f. At M1's variogram skill (0.43, chamber-cut backtest) that is about 20 x 0.18 x dVS/d_f.
  - The toy's dVS (1.7e-3 for a strong change with mean |log f| about 0.18, 1.7e-4 at 0.3x) scaled to a mean |log f| of 0.05-0.1, against a variogram d_f of order 5e-3 (validation-board anchors, docs/CHALLENGE.md), gives dVS/d_f of about 0.03-0.1, i.e. +0.1 to +0.4 T1 points if the increases are right. This agrees with the ledger's pruned "exact top-up on T1" estimate of 0.1-0.4.
  - Assumption: Feng's mean |log f| at 2/5 is unknown here; `external_shift`'s report (mean_abs_log_f per group) gives it target-free.
- **T2 heart:** at MERFISH depth the own-count top-up realises 10-30% of B0's increases. Its detection loss mimics a depth decrease, which fits the falling heart depth series (233 to 56 median counts). An exact kernel would remove that hidden depth-like effect, so the swap is a depth-lever question in disguise (pruned b0_smooth, 3-5% to pass).
- **T2 embryo:** no count edits. A 1:1 mixture of brackets at depths 62 and 212 has a Jensen gap: detection follows the arithmetic-mean depth (137), not the geometric one (115), and single-count values form two levels. Symmetric harmonisation would remove this: thin the E8.0 half and EB-top-up the E7.5 half to the geometric mean, keeping the mean log depth, which is E7's failure mode avoided. The depth line is closed.

**5. A backtest that holds the direction fixed (designed to count where a counting task exists).**
- **R1, thin-and-restore (target-free; a gate in the spirit of E10's gate; uses the data-thinning identity [1][2]).**
  - Input stage: one that is not a counting target. For T1 that is released E9.5 (chamber classes); for the heart, E9.5 MERFISH; for T3, WT E8.75 if ever needed.
  - Factor: the entry's actual factor matrix F (T1: Feng at 2/5 per group with the held-flat genes; heart: B0's momentum), up part only, F+ = max(F, 1).
  - Split the cells by a seeded, type-stratified permutation into A and B. Set B- = Bin(B, 1/F+). B- is an exact before-population and A an exact after-population for any expression model.
  - Realise F+ on B- with each kernel (own-count, kNN plug-in, EB), keep_depth off, at the build's count, each with its own RNG stream so the cell draws are identical.
  - Score against A with B- as reference, using the pinned veckit metrics, the 10% working copies and the split half (vec.evaluate), 10 seeds, paired.
  - Readouts: board points; skills per metric; realised fraction of the intended dpb by sparsity decile; detection and E sqrt x errors.
  - Because the direction is exact, the gain here is an upper bound on what any backtest or the test can give. Fixed in advance: register a counted try only if EB or kNN beats own-count by more than 0.3 points in at least 8 of 10 seeds; otherwise prune.
  - It reads no target and no test file and scores synthetic truths, so by the synthetic-toy precedent it is not a try.
- **R2, kernel swap on the counting backtests (counts under the adoption rule).**
  - Only the up-kernel changes inside the entry's own pipeline: T1 chamber-cut E9.5 -> E10.5 with B0's momentum; heart, both backtests with B0's momentum. Same seeds and draws; row registered before the run.
  - Fit-not-effect conditions written in advance: (i) the gain must sit on the variogram (T1) or on the variogram, mmd_u and de (heart); (ii) the gain must not track a detection or depth shift toward the target, so report both arms' mean detected genes per cell and median library size against the target's; (iii) T1 is reported per category (CM against non-CM).
- **R3, oracle-composition diagnostic (a counted try that can keep nothing).** R2 with both arms drawn at the target's category shares. It removes the CM-share confound (64% against 19%) that otherwise makes R2 on T1 judge B0's direction rather than the kernel.

**6. Toy evidence.** Synthetic NB-Poisson world: 4 types plus a within-type programme, own numpy re-implementations of the published formulas. "MERFISH-like" is 500 genes at about 100 counts per cell; "T1-like" is 2,500 genes (oracle rows) or 1,500 genes (estimated rows) at about 1,000 counts per cell, close to T1's counts per gene. Toy points use floor = the reference copy and ceiling = an independent truth draw; they are not board scores. Run times: see the Compliance log.

| Regime and direction (seeds) | own-count (T1 and B0 now) | kNN (T3 now) | EB | NB resample | oracle theta |
|---|---|---|---|---|---|
| MERFISH-like, oracle factor, scale truth (2) | 82.7, 85.0 | 99.2, 98.6 | 98.9, 97.9 | 93.6, 92.8 | 99.9, 99.4 |
| MERFISH-like, oracle factor, frequency truth (2) | 80.5, 86.1 | 99.4, 98.9 | 97.9, 98.7 | 91.2, 94.0 | 99.8, 99.7 |
| MERFISH-like, estimated factor, scale truth (2) | 67.7, 72.8 | 78.0, 83.7 | 76.8, 83.1 | 70.0, 79.7 | 77.5, 83.8 |
| MERFISH-like, estimated factor, frequency truth (2) | 65.7, 74.1 | 74.5, 82.3 | 73.8, 82.7 | 69.4, 76.2 | 74.4, 83.1 |
| T1-like, oracle factor (2) | 81.1, 89.7 | 98.0, 99.9 | 96.3, 100.0 | 87.4, 91.2 | 98.1, 100.0 |
| T1-like, estimated factor (1) | 71.2 | 74.4 | 74.7 | - | 74.9 |
| T1-like, change at 0.3x, estimated factor (1) | 66.1 | 72.0 | 72.2 | - | 72.7 |

- The estimated factor comes from 300 cells per type with the eps pseudocount.
- Where own-count loses: the variogram (T1-like skill 0.68-0.78 against 0.77-1.0 for exact kernels), mmd_u, and on MERFISH also severity. With the oracle factor, own-count has log beta -0.7 to -1.0 and r2 0.47-0.62, against about 0 and 0.97 or more for exact kernels; with the estimated factor, -1.2 to -1.6 and 0.35-0.64, against -0.6 to -1.0 and 0.69-0.88.
- Toy gaps overstate the test value: the toy's floor-to-ceiling gaps are set by its own change, and its direction error is far milder than Feng's capture and platform effects.

## Candidate levers
Ranked by expected value to a top-3 place. None is worth more than about +0.3 Overall in expectation.

1. **L1. Exact up-kernel for T1's direction (EB posterior-predictive top-up, or T3's kNN plug-in, inside `apply_factors`).**
   - Board and group: T1:test, covariance (variogram, 0.20) first, then DE recovery and direction.
   - Mechanism: replace c + Poisson((f-1)c) with c + Poisson((f-1) L theta~), keep_depth unchanged; s1 and s2 both carry it.
   - Expected effect: +0.1 to +0.3 T1 points if Feng's increases are right (scales like (log f)^2); expected value about +0.1, range -0.3 to +0.5.
   - Legality: an expression edit, which stays open. Depth is kept in expectation, so median library size stays within 5% (report it). No composition, count, scope or cross-modal change. The core changes, so the golden table in tests/test_core.py moves in the same commit, and the s1 change sends s2 back to the captain under the frozen s2 rules.
   - Validation: R1 gate on E9.5, then R2 on the chamber cut. R2 is confounded by the CM share, so a pass is unlikely (about 15%) and a fail is weak evidence. Realistically R1 evidence plus the captain's blind bet.
   - Effort: 1-1.5 days.
   - Main risk: own-count under-delivery currently acts as implicit shrinkage of a noisy direction. Feng's sparse-gene ratios carry zone-capture and platform artefacts, and delivering them fully can lower de_score and the variogram.
2. **L2. Exact up-kernel in the heart's B0 (revives the pruned b0_smooth).**
   - Board and group: T2:heart:test_extrap, expression change and cell state.
   - Mechanism: as L1, in B0's momentum.
   - Expected effect: 0 +/- 0.3 board points (+/-0.15 Overall).
   - Legality: an expression edit only; counts, sizes, scope and composition untouched; no cross-modal input.
   - Validation: R1 on E9.5 MERFISH, then R2 on both heart backtests with the depth-confound condition. Pass chance 5-10%.
   - Effort: 1 day.
   - Main risk: the own-count detection loss is a hidden partial depth correction in the direction of the falling depth trend, so an exact kernel may lose at test.
3. **L3. T1 within-category state shares by a bounded single-moment I-projection, realised by quotas.**
   - Board and group: T1:test, distribution (mmd_u) and direction.
   - Mechanism: within each Feng group, raking over E10.5 sub-states with w ~ exp(lambda s), where s is the cell's score along the group's Feng direction. lambda is set so the share-weighted mean score moves by the intended amount, capped at the no-duplicate bound phi(Phi^-1(1-p))/p. Shares go through `quotas` at the frozen count, and the residual of the direction is count-edited (the hybrid).
   - Expected effect: 0 +/- 0.7. The oracle per-type shift moved mmd_u by only 0.12-0.30 points (JL-OT entry), so state shifts within a type are worth little.
   - Legality: composition by the legal route of 2026-11-04 (d) and (f) only. Rule and constants committed in the ledger row before any run; E10.5 plus legal Feng data only; quotas, never weights; report the library shift.
   - Validation: none counted (Feng's E9.5 lies inside both T1 windows, and both counted tasks carry E9.5's mix), so the captain's blind bet only.
   - Effort: 2-3 days.
   - Main risk: Feng's E9.5 -> E14.5 capture change (from mostly LV zone to zone-balanced) means the score partly selects cells by zone, a composition artefact.
4. **L4. Pseudocount-consistent switch-on term for T1's increases.**
   - Board and group: T1:test, DE recovery.
   - Mechanism: add Poisson((f-1) eps L_i/1e4) to genes with f > 1, so the realised change matches the estimator's implied (f-1)(m + eps).
   - Expected effect: 0 +/- 0.5, sign unknown.
   - Legality: an expression edit using no external level (ratio and our own constant), but it changes what the ruled Feng direction delivers, so it needs the captain's ruling.
   - Validation: R1 cannot test it (R1 is multiplicative by construction); R2 with B0's eps faces the CM-share confound. Captain's blind bet only.
   - Effort: 0.5 day.
   - Main risk: it imports Feng's zone, capture and platform artefacts as genes switched on in every cell of a group.
5. **L5. Embryo symmetric depth harmonisation. CLOSED** (the depth line was closed by the ledger after E7 and H3, and the embryo is at "stop polishing").
   - Board and group: T2:embryo:test_interp, cell state.
   - Mechanism: thin the E8.0 half and EB-top-up the E7.5 half to the brackets' geometric-mean depth, keeping the mean log depth.
   - Expected effect: +0.05 to +0.3 board points (half in Overall), only if the captain reopens the line.
   - Legality: no counts, sizes or scope change, but the depth line is closed.
   - Validation: both embryo backtests, with E7.25's depth closeness to E6.75 as a known feature written in advance.
   - Effort: 1 day.
   - Main risk: thinned cells do not read like shallow acquisitions (E7, H3).
6. **L6. T3 kernel change. CLOSED, and nothing to gain:** the frozen kNN kernel is near-exact in the toy; changing it would be a frozen-rule change by the captain.
7. **L7. Calibrated exact-total realisation (hypergeometric thinning, fixed-total top-up, cube balancing).** Legal everywhere, but value is negligible next to the 10% reference-copy noise (T5 precedent). Not worth a try.
8. **L8. NB resampling. Rejected:** it is legal as an expression edit, but it loses the cell's residual structure (toy: 3-9 points below EB, mostly on mmd_u and the variogram).
9. **L9. Latent Gaussian OT with nearest-real-cell decoding. Closed on T2 and T3, rejected on T1:** it is an exponential tilt (a composition change). On T1 it reduces to L3 with worse ESS and duplicates; on T2 it also breaks the expression-position coupling.

## Deep Read Notes
- [1] Gerard 2020. A Poisson log-link model with Theta-tilde = Theta + Q; thinning probabilities 2^q with q <= 0 after subtracting each row's maximum ("thin to the max"); Lemma 1 is Poisson thinning, with NB and mixture extensions in the supplement. Unwanted variation is untouched, and the mean log fold change is right even for non-Poisson data, but the mean-variance shape may not be. No upward kernel exists: positive effects cost depth.
- [2] Neufeld et al. Folds are independent and sum to X. Poisson splits by binomial or multinomial; NB(r, p) by beta-binomial (epsilon r, (1-epsilon) r), scaling r; gamma by a Beta split; binomial by hypergeometric. The NB split requires a known r, and a wrong r leaves fold covariance proportional to 1 - (r+1)/(r~+1). This supplies the frequency-mode kernels in section 2.
- [3] Sarkar and Stephens. Poisson measurement with independent capture per molecule; NB and ZINB are Poisson measurement over Gamma or point-Gamma expression. Point-Gamma is favoured for 0.6-9% of genes, a unimodal nonparametric expression model for 20-69%, multimodal for 0.4-3%. Measurement dispersion is small, and estimates of true expression are generally not simple functions of the counts, which supports an EB posterior over plug-ins.
- [5] ebal manual. Control weights from a maximum-entropy scheme satisfy researcher-specified moments (squared terms for variances, interactions for co-moments). ebalance.trim lowers the max/mean weight ratio under the same moments and flags infeasible targets, the practical form of bounded calibration.

## Gaps
- No source quantifies how a zero-preserving own-count Poisson top-up distorts log-normalised metrics, and none names a posterior-predictive top-up for adding increases to real counts: [1] avoids increases by thinning everything else. The table and kernel in section 2 are my own derivations plus the toy.
- I could not read the full texts of [4], [6], [7] and [9] (paywalled or PDF only); I used their records and the standard results.
- I found no evidence on whether embryonic or knockout changes act on burst frequency or burst size gene by gene, and no study of MERFISH measurement as Poisson (spillover and segmentation would break it).
- Nothing explains why thinned cells read unlike a shallower acquisition (E7, H3).
- **Counter-claims against L1:**
  - The current under-delivery may be beneficial shrinkage of a noisy direction; the toy's estimated direction was far cleaner than Feng's.
  - The expected gain scales like (log f)^2 and may be invisible at 2/5.
  - E7 and H3 show that changing detection profiles can cost mmd_u even when depth matches.
  - On the heart, exactness removes a depth-like effect that may be helping.

## Compliance log
- `.venv/lib/python3.12/site-packages/common/core_metrics.py`, the module docstring and the `de_genes` docstring: example truth-DE counts that `docs/DATA_POLICY.md` 2026-11-05 (i) says may describe a T3 board. Met while reading the formulas; not copied or used. Class A (possible).
- Process notes (not barred items, nothing opened): the harness saved two WebFetch PDFs (the Sarkar-Stephens bioRxiv PDF and a Csiszar-Shields tutorial PDF) and one oversized `sed` output of `experiments/harness/TRIES.md` into the session store under `~/.claude-account`; none was opened. No leaderboard, board file, `data/`, `runs/` or unscreened research file was opened, and no target was read.
- Run-time slip: the synthetic toy, which reads no data file, ran with the repository venv's `python -I -B` from the scratch folder. Two runs went over the one-minute limit on the loaded local machine (load average about 85): the first was killed at 120 s by BLAS thread oversubscription, and one T1-like run took 88 s. All later runs used 4 BLAS threads and smaller sizes and finished in 14-36 s.
## END
