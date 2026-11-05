---
task_id: f
role: Applied mathematician (dynamical systems, optimal transport, growth processes, time warping)
status: complete
sources_found: 18
---
## Sources
[1] Why the damped trend works (Gardner & McKenzie; revised preprint of J Oper Res Soc 2011) | https://www.bauer.uh.edu/gardner/docs/pdf/Why-the-damped-trend-works.pdf | Source-Type: academic | Accessibility: public | As Of: 2009-10 | Authority: 8/10
[2] Unmasking the Theta method (Hyndman & Billah; Int J Forecasting 19(2) 287-290; working-paper record) | https://econpapers.repec.org/paper/mshebswps/2001-5.htm | Source-Type: academic | Accessibility: public | As Of: 2003-04 | Authority: 8/10
[3] forecast R package, theta_model / forecast.theta_model source (Hyndman et al.) | https://raw.githubusercontent.com/robjhyndman/forecast/master/R/theta.R | Source-Type: other | Accessibility: public | As Of: 2026-11 | Authority: 7/10
[4] From Snapshots to Trajectories: Learning Single-Cell Gene Expression Dynamics via Conditional Flow Matching (scFM; Pu et al.) | https://arxiv.org/html/2605.22340v1 | Source-Type: academic | Accessibility: public | As Of: 2026-05 | Authority: 5/10
[5] TIMO: Modeling Cell Developmental Trajectory using Multinomial Unbalanced Optimal Transport (Zhu, Zhang, Zhang, Kong) | https://arxiv.org/html/2501.03501v2 | Source-Type: academic | Accessibility: public | As Of: 2026-03 | Authority: 5/10
[6] WFR-FM: Simulation-free dynamic unbalanced optimal transport (Peng et al.) | https://arxiv.org/html/2601.06810 | Source-Type: academic | Accessibility: public | As Of: 2026-01 | Authority: 5/10
[7] Reconstructing growth and dynamic trajectories from single-cell transcriptomics data (TIGON; Sha, Qiu, Zhou, Nie; Nat Mach Intell) | https://www.nature.com/articles/s42256-023-00763-w.pdf | Source-Type: academic | Accessibility: public | As Of: 2024-01 | Authority: 7/10
[8] Second order models for optimal transport and cubic splines on the Wasserstein space (Benamou, Gallouet, Vialard) | https://pdf.arxiv.org/pdf/1801.04144 | Source-Type: academic | Accessibility: public | As Of: 2018-01 | Authority: 7/10
[9] Approximation of splines in Wasserstein spaces | https://arxiv.org/pdf/2302.10682 | Source-Type: academic | Accessibility: public | As Of: 2023-02 | Authority: 5/10
[10] Fast and smooth interpolation on Wasserstein space (Chewi et al., AISTATS 2021) | https://proceedings.mlr.press/v130/chewi21a | Source-Type: academic | Accessibility: public | As Of: 2021-04 | Authority: 7/10
[11] Gene-level alignment of single-cell trajectories (Genes2Genes; Sumanaweera et al., Nat Methods) | https://suochenqu.github.io/publication/sumanaweera2025gene.pdf | Source-Type: academic | Accessibility: public | As Of: 2024-10 | Authority: 7/10
[12] Real age prediction from the transcriptome with RAPToR (Bulteau & Francesconi, Nat Methods) | https://www.nature.com/articles/s41592-022-01540-0 | Source-Type: academic | Accessibility: public | As Of: 2022-07 | Authority: 8/10
[13] Species-specific segmentation clock periods are due to differential biochemical reaction speeds (Matsuda et al., Science) | https://repositori.upf.edu/handle/10230/52326 | Source-Type: academic | Accessibility: public | As Of: 2020-09 | Authority: 8/10
[14] Waves, patterns and bifurcations: a tutorial review on the vertebrate segmentation clock | https://arxiv.org/pdf/2403.00457 | Source-Type: academic | Accessibility: public | As Of: 2024-03 | Authority: 6/10
[15] A statistical analysis of memory CD8 T cell differentiation: a hierarchical state space model for a short time course | https://arxiv.org/pdf/0712.1124 | Source-Type: academic | Accessibility: public | As Of: 2007-12 | Authority: 5/10
[16] M4 competition multi-factor analysis (MDPI Data 5(2):41) | https://mdpi-res.com/d_attachment/data/data-05-00041/article_deploy/data-05-00041-v2.pdf | Source-Type: academic | Accessibility: public | As Of: 2020-04 | Authority: 5/10
[17] Oh SnapMMD! Forecasting stochastic dynamics beyond the Schrodinger bridge's end | https://arxiv.org/html/2505.16082v1 | Source-Type: academic | Accessibility: public | As Of: 2025-05 | Authority: 4/10
[18] Benchmarking zero-shot single-cell foundation model embeddings for cellular dynamics reconstruction | https://www.biorxiv.org/content/10.64898/2026.03.10.710748v1.full | Source-Type: academic | Accessibility: public | As Of: 2026-03 | Authority: 4/10

## Findings
1. Fitting the damped trend to the 3,003 M3 series selected a damped drift or trend in 84% of series (local initial values; 70% with global ones) and an undamped linear trend in 0.1% (7.9%), and its random-coefficient state-space reading makes the damping phi the per-step probability that the current trend persists, with h-step forecast l + (phi + ... + phi^h) b [1].
2. The Theta method, one of the M4 statistical benchmarks beside damped and Comb, is simple exponential smoothing plus a drift equal to half the slope of the fitted linear trend, i.e. a fixed 0.5 damping of the extrapolated trend [2][3][16].
3. In a May 2026 benchmark on dense series (10-18 time points, latent-space Sinkhorn W2), a copy-last "Naive" baseline beat VGFM in 7 of 8 extrapolation settings while scNODE and the authors' scFM beat it everywhere, so gains over copy-last are method- and regime-dependent and shown only on dense, single-regime series [4].
4. A March 2026 preprint found single-cell foundation-model embeddings underperform a plain HVG baseline when the last time point is held out for extrapolation [18].
5. Per-cell growth rates from unbalanced OT show inflated within-type variance, metacell-level semi-relaxed OT gives steadier rates, and growth cannot be uniquely separated from differentiation using snapshots [5].
6. Growth read from hallmark gene signatures (Waddington-OT, PRESCIENT) depends strongly on which gene-set database is chosen [7], while Wasserstein-Fisher-Rao flow matching learns a velocity and a time-varying growth field jointly from several snapshots [6].
7. Cubic splines on Wasserstein space are well defined and computable, but spline trajectories overshoot where piecewise geodesics do not [8][9][10].
8. Developmental tempo is set by cell-intrinsic biochemical rates (mouse segmentation clock about 2-3 h against 5-6 h in human), and that scaling may be specific to the segmentation clock rather than general [13][14].
9. Transcriptome age-staging (RAPToR) reads a sample's age off an interpolated reference time series and gene-level DTW alignment (Genes2Genes) matches trajectories observed over the same range, so both need reference data covering the age being placed [11][12].
10. Hierarchical empirical-Bayes state-space models fit four-point expression time courses by borrowing strength across thousands of genes [15].

## First-principles analysis
Notation: per cell type k and gene g, y_t = stage pseudobulk of log1p CP10k; D = Feng's within-dataset E9.5 -> E14.5 log ratio; m = y(E10.5) - y(E9.5). Assumptions are marked A1-A5. Repository facts come from committed files (cited by path); nothing below reads a target.

**0. What the scorer rewards in an extrapolated delta (veckit cc3859a, formulas only).** Board points per group = 100 w (skill - 0.5), skill = d_floor / (d_floor + d), d = distance to the split-half ceiling. de_score ranks genes by the prediction's delta against the reference working copy, takes the truth's n_up / n_dn counts, and subtracts an expression-level null; de_direction is a rank partial correlation of predicted and true deltas with the reference level partialled out. Both are invariant to a positive rescaling of the delta and saturate once the shift dominates the sampling noise of the 10% reference copy and of our draw. mmd_u is an unbiased multi-bandwidth RBF MMD^2 in a 30-PC basis fitted on the truth; for a within-type mean-shift error e small against the bandwidth it is roughly quadratic in e. variogram compares E|x_i - x_j|^(1/2) over random gene pairs and behaves mostly as a per-gene marginal and sparsity statistic. So in T1 the rank metrics (weight 0.50) want a delta of the right identity and sign and are flat in its size; the magnitude metrics (0.50) want the size f* = <T, D> / ||D||^2, the regression of the true change T on the external direction D. If T = f0 nu D_t with corr(D_t, D) = rho and ||D_t|| = ||D||, then f* = f0 rho nu: damping read as regression-dilution (James-Stein) shrinkage.
- Synthetic toy (4,000 genes, 10 seeds, no data): with rho = 0.7 and nu = 1 the de_direction-like score is 0.627 at f = 0.28, 0.638 at 0.40 and 0.640 at 0.50; the de_score-like 0.456, 0.465, 0.467; the mean-shift error (relative to a copy) 0.50, 0.59, 0.80. With rho = 0.9 the error is lowest near 0.40; with rho = 0.5 near 0.20.
- Committed evidence on size: B0 -> B0 + Feng at 1/5 moved every T1 skill up a little (de_score 0.450 -> 0.480, de_direction 0.608 -> 0.620, mmd_u 0.493 -> 0.497, variogram 0.423 -> 0.428; `experiments/t1_m1/results/diag_v1/README.md`), and an oracle per-type mean shift moved mmd_u only +0.12 and +0.30 points (TRIES, JL-OT item). The mean-shift share of mmd_u is small, so T1 is flat in the step size over about 0.28-0.5 of D to within roughly +-0.2 points. The step size is not a lever.

**1. Momentum and damping as shrinkage (Q1, Q6).** A1: y_t = mu_t + delta_t + e_t, with mu smooth development, delta a sampling-regime offset (embryo portion against dissected heart) and e a stage-specific idiosyncrasy (litter, sex ratio, run) independent across stages. Under the random-coefficient model of [1], v_{t+1} = A_t v_t with P(A_t = 1) = phi, so E[y_{t+h} - y_t | v_t] = v_t (phi + ... + phi^h). For pre-break increments Delta_i over Delta t_i:
- Var(Delta_i) = s_v^2 Delta t_i^2 + 2 s_e^2;
- Cov(Delta_i, Delta_{i+1}) = phi s_v^2 Delta t_i Delta t_{i+1} - s_e^2.
The -s_e^2 term is the MA(1) fingerprint of stage noise. It is target-free and estimable by moments across thousands of genes [15], and it gives the empirical-Bayes velocity factor S(Delta t) = s_v^2 Delta t^2 / (s_v^2 Delta t^2 + 2 s_e^2). Theta's fixed 1/2 [2][3] and B0's damping 0.5 are this factor frozen at one value.
- Committed in-sample reading (`experiments/t1_direction/results/E10.5/per_type.csv`): the E8.5 -> E9.5 momentum predicts the E9.5 -> E10.5 change with Spearman 0.49 for CM, 0.44 for endocardium, -0.53 for endothelium, -0.14 for ventricular CM and -0.06 for blood. A negative lag-1 association is what a large s_e (plus delta) produces, so S is about 0 for those types. M1 gives the Feng groups no momentum (Feng's direction instead), and R1's hold is the same empirical-Bayes answer for blood.
- Momentum in W2: for per-type translations, x -> x + s v stays an optimal map for every s. For a cell-level Brenier map T = grad phi, the McCann extension grad[(1-s)|x|^2/2 + s phi] stays the gradient of a convex function only while (1-s) I + s Hess(phi) is positive semidefinite, i.e. s <= 1/(1 - lambda_min) when lambda_min < 1. A population that converges along some direction (lambda_min near 0) cannot be carried past s of about 1 without crossing paths. A 2-day push from a 1-day step is s = 2, so cell-level Wasserstein momentum is ill-posed at the test horizon; per-type translations, which the entries use, are the stable subclass.

**2. The sampling break decides the state-space answer (Q6, Q8).** A2: delta switches on at E10.5 (dissected heart) and persists to E12.5 (also a dissection). In a local-linear-trend Kalman model with an intervention at the last observed stage, the post-break level is that observation (effective alpha = 1). The across-break increment m = mu-dot + delta + e_10.5 - e_9.5 is not a velocity, so the slope must come from pre-break increments (S small, above) or from a regime-free external source. That is M1 after R1 exactly: level = E10.5 cells, slope = Feng's within-dataset direction, horizon factor = bracket fraction. It also explains why removing the reference's own deviation is degenerate (heart leverage 0.80-0.91, TRIES).
- The heart differs. B0 at test moves each lineage by half of E9.5 -> E10.5, an across-break increment. Both counted heart backtests move by half of a pre-break increment (E8.75 -> E9.5; E8.5 -> E8.75; `configs/pseudo_tasks.yaml`) into a post-break target. So the backtests score a different estimator from the one used at test (a known feature under TRIES rule 2), and a placebo reading on them (H5) inherits the gap. On the 1.75-day backtest, B0's gain over a copy also sat on neighbourhoods and cell state, not on expression change (METHOD_SUMMARY, heart table).

**3. Geodesic to the bracket against constant velocity and non-autonomous fields (Q2, Q3).** Let r = m - D/5, the deviation of E10.5 from the linear schedule along Feng's bracket. Every candidate path has the form y_hat(12.5) = y(10.5) + (2/5) D + kappa r:
- kappa = 0: constant Feng velocity (M1);
- kappa = -1/2: the McCann chord from E10.5 to the Feng-anchored endpoint y(9.5) + D at s = (12.5-10.5)/(14.5-10.5) = 1/2 (geodesic to the legal bracket);
- kappa = +1/2: quadratic Lagrange through E9.5, E10.5 and E14.5 (weights -0.8, +1.5, +0.3), a constant-acceleration non-autonomous field.
If r is per-gene switch timing (a gene that switched by E10.5 has nothing left; in a uniform-switch model an unswitched gene has 2/4 of D still to come), kappa = -1/2 follows. If r is the persistent regime offset delta, which E12.5 shares, kappa = 0. If r is acceleration, kappa = +1/2. The evidence disfavours kappa > 0: mom+feng is worse than Feng alone in the readout (CM Spearman 0.65 against 0.69; ventricular CM 0.33 against 0.67; in-sample), and the digest's within-Feng late-stage check found momentum adds nothing. kappa < 0 needs r to be development rather than delta, which one post-break stage cannot show. So kappa = 0 is the minimax choice.
- Spline (non-autonomous) extrapolation [8][10] needs at least three same-regime points, T1 has one post-break point, and splines overshoot [9].
- Embryo: a cubic through E6.75, E7.25, E7.5 and E8.0, evaluated at E7.75, weights them +0.07, -0.67, +1.33 and +0.27. That is Runge-type leverage on stages whose median depths run 30, 35, 62 and 212 counts (TRIES). The current piecewise geodesic at 1/2 weights them 0, 0, 0.5 and 0.5 and is the stable choice.

**4. Time warps and the bracket fraction (Q5).** f = [tau(12.5) - tau(10.5)] / [tau(14.5) - tau(9.5)] (exact arithmetic):

| Clock tau | f |
|---|---|
| linear (chronological) | 0.400 |
| log age | 0.412 |
| speed proportional to age^-2 | 0.420 |
| speed proportional to age^-3 | 0.422 |
| accelerating, speed proportional to age | 0.383 |
| relaxation 1 - exp(-(t - 9.5)/T), T = 1 | 0.320 |
| relaxation, T = 2 | 0.418 |
| relaxation, T = 3 to 3.5 (the maximum) | 0.430 |
| relaxation, T = 5 | 0.427 |
| relaxation, T = 10 | 0.417 |
| power progress, p = 0.5 | 0.327 |
| power progress, p = 1.5 | 0.375 |
| power progress, p = 2 | 0.320 |
| one switch per gene at a uniform time | 0.400 (mean) |

Staging noise of +-0.5 day on Feng's bracket gives 2/(5 +- 0.5), which is 0.36-0.44. Any smooth, moderately curved clock moves f by at most about +-0.04 from 2/5, and section 0 makes T1 flat over that range. Linear 2/5 is also the Bayes value under a uniform switch-time prior. The time warp is not a T1 lever.
- Embryo: a legal clock can be read only outside (E7.25, E8.0). Reference staging [12] and DTW alignment [11] need in-window references, which are barred. Clock-period scaling [13][14] gives tempo, not position. Depth (62 against 212 counts) confounds any data-driven distance clock between the brackets.
- Toy (b), synthetic: brackets N(0,1) and N(D,1), target an intermediate state N(t* D, 1), the scorer's five bandwidth multipliers. The mmd_u-optimal mixture weight equals t* for D = 1 but is 0.31 at D = 2 and 0.28 at D = 4 for t* = 1/3, i.e. pushed away from 1/2. A known warp would have to be applied amplified, but under a symmetric prior on t* the Bayes weight stays 1/2. The pruned clock warp stays pruned.

**5. Growth via unbalanced OT (Q4).** The state-level reaction term (WFR [6], McKendrick) gives w_s(t+h) proportional to w_s(t) exp(h g_s), with g_s = b_s - d_s plus net transition flux. In steady asynchronous growth the cell-age density is 2 lambda exp(-lambda a) on [0, T_c], so the fraction of cells in the last tau of the cycle (S+G2/M) is c = exp(lambda tau) - 1, and lambda = ln(1 + c)/tau (A3: tau of about 0.3-0.4 day, generic mammalian; a G0 fraction lowers c and lambda together). Example: c = 0.3 against 0.5 at tau = 0.4 day gives a gap of about 0.36 per day and a relative share factor of about 2 over two days. That is a large move driven by an unvalidated rate. Death and transition fluxes are not identified from one post-break snapshot [5], and signature-based rates depend on the gene-set database [7].
- Per-cell weights are both forbidden (DATA_POLICY 2026-11-04 (f)) and statistically noisy [5]; state-level quotas are the stable unit, and they match the legal route's form.
- An own-growth quota between states, using no external data, is legal on the counted E9.5 -> E10.5 tasks. But there the E9.5 embryo-portion mix (cardiomyocytes about 64% of cells against about 19% in the chamber-cut truth) swamps any one-day growth signal, so rule 2 fails before the run.
- A maturation variant moves the cycling against non-cycling split within each category along a legal same-platform bracket instead (lever L1).

**6. Benchmarks (Q7, beyond the digest).**
- Dense-series forecasting can beat copy-last in latent W2 [4], but one learned method lost to it in 7 of 8 extrapolation settings, foundation-model embeddings lose to HVG baselines [18], and MMD-fitted forecasting beyond a bridge's end exists only on synthetic and T-cell data [17].
- In forecasting competitions, damped and half-damped trends are the benchmarks to beat [1][2][16].
- None of this covers three stages with a regime break, which is our setting. Within VEC the digest already records that the organisers' neural ODE does not beat a constant shift.

**7. Keeping extrapolation honest (Q8).**

| Change | Can a counted backtest judge it as used at test? |
|---|---|
| Feng's scale, its fraction, kappa, any JAX- or Feng-derived quota | No: Feng's (and JAX's) E9.5 lies inside every T1 window |
| T1 own-growth quota | Counts formally, but is confounded by E9.5's portion mix (rule 2) |
| Heart momentum calibration | Only in its pre-break form, which the test does not use |
| Embryo mixing fraction | No (TRIES, campaign of 2026-11-05) |

The cross-fitting recipe for any damping or shrinkage constant: estimate it by the moment identities of section 1 from increments strictly before the backtest's reference, rerun the same estimator at test, horizon-normalise it as a per-day rate times (phi + ... + phi^h), and pair it with a decile-matched permuted-gene placebo (as K2b).

**8. Bottom line for the top-3 objective.**
- The dynamical-systems view confirms M1 after R1 as the state-space optimum given the break, and closes the time warp, the step size and kappa.
- It flags one untested structural risk, the heart's across-break momentum.
- It leaves one blind bet with real upside, L1.
- Expected Overall gain from everything below: about +0.2 to +0.6, against a needed +15-25.

When far below a threshold, P(top-3) is driven by upside variance, so L1 is worth more as a decorrelating bet in one upload slot than in both (A4; a captain's call, because the frozen s2 rule builds T1 s2 from s1).

## Candidate levers
Ranked by expected value to the top-3 objective. Nothing here passes the adoption rule on a counting backtest in the form used at test.

**L1. Maturation quota: the cycling-exit split within chamber categories (rank 1; open, blind)**
- Board and group: T1:test, distribution (mmd_u, 0.30; 15.08 left), secondarily covariance (variogram, 0.20).
- Mechanism: within cardiomyocytes, endocardium, endothelium and epicardium (blood excluded, as R1 does), call S/G2M against non-cycling nuclei with a generic cell-cycle gene set. Shift each category's cycling log-odds from E10.5's own value by f x [logit c_X(late) - logit c_X(E9.5)], drawn by vec.core quotas at the frozen count.
  - X is JAX whole-embryo nuclei (E9.5 and at or after E13.75, the deeper embryos thinned to the shallower depth before scoring), or Feng computed within zone so that its capture shift cancels.
  - f = 2/(t_late - 9.5) (2/5 when the late stage is E14.5).
- Expected effect: +0.3 to +1.5 T1 points if the cycling fraction falls as the bracket implies; -0.3 to -1.5 if the sign or size is wrong; EV about +0.2 to +0.5 (A5: sign right with p of about 0.7; magnitude unknown). The cell cycle is a major PCA axis, so within-category state shares move mmd_u, whereas per-type mean shifts barely do (section 0).
- Legality:
  - Composition of states within categories by the legal route: rule and constants committed in the TRIES row before any run; quotas, never weights; inputs only unlocked E10.5 and legal external data. JAX rows follow the Q7 staging rule. An X002/X003 amendment and an EXTERNAL_DATA row for the gene set (PLAN 4.4 provenance; no Ctnnb1/Gata4 perturbation origin) come before use.
  - The cardiomyocyte and blood shares are unchanged.
  - The drawn median library size is reported against E10.5's; a shift above 5% goes to the captain.
  - T1 depth is untouched (cells keep their counts). No size, count or scope change, and no cross-modal input.
  - s2 follows s1 under its frozen rule, and a changed s1 sends s2 back to the captain.
- Validation route: captain's blind bet only, after he writes a prior (TRIES item 1); no counted T1 task can judge it. Target-free pre-checks: split-half stability of the external log-odds change, and the depth shift.
- Effort: 2-3 days.
- Main risk: cycling calls are depth-confounded (JAX's late runs are shallower; nuclei against cells), so the shift may be technical. If the dissected chamber tissue does not lose cycling share by E12.5, it loses.

**L2. Heart: no momentum across the sampling break (an R1 analogue), or cardiomyocyte-only momentum (rank 2; open, ruling)**
- Board and group: T2:heart:test_extrap; expression_change, cell_state and local_organisation (0.25 each).
- Mechanism: replace "half the E9.5 -> E10.5 change of each harmonised lineage" with a hold, or keep it only where pre-break and across-break changes agree (the pruned-items note says per-type momentum carries over only for cardiomyocytes).
- Expected effect: -0.5 to +0.5 heart-board points (about +-0.25 Overall), EV about +0.05. Its value is removing an untested test-only estimator (section 2).
- Legality: expression only; no size, count, scope or composition change; no cross-modal input. The frozen s2 rule rebuilds the heart s2 from s1, so it needs the captain's ruling.
- Validation route: no counting backtest (neither crosses the break). Captain's ruling, as R1 was. An information-only per-lineage concordance reading (pre-break against across-break change) would count as a try.
- Effort: 0.5-1 day.
- Main risk: it discards real cardiomyocyte maturation signal, and part of B0's backtest gain on cell state and neighbourhoods.

**L3. Empirical-Bayes / Kalman calibration of the heart's pre-break velocity (closed by rules 5 and 2)**
- Board and group: heart, expression_change.
- Mechanism: the moment identities of section 1 on E8.25_late -> E9.5, a per-day S(Delta t) and phi, horizon-normalised.
- Expected effect: under 0.3 board points on the backtests, and it is not the estimator used at test (the break).
- Legality: clean; it reads no target if cross-fitted.
- Validation route: counted heart backtests, but rule 2's known feature applies.
- Effort: 1 day.
- Main risk: a lucky fit. Prune.

**L4. James-Stein scale of Feng's step, f = (2/5) rho with rho of about 0.7 (closed by analysis)**
- Board and group: T1, distribution and covariance.
- Expected effect: within about +-0.2 points; the rank metrics lose roughly what the magnitude metrics gain (section 0).
- Legality: a change to Feng's direction.
- Validation route: blind.
- Effort: 0.5 day.
- Main risk: none worth taking. Not recommended.

**L5. Time-warped bracket fraction for T1 (closed by analysis)**
- Board and group: T1.
- Mechanism: f from log age or a relaxation clock, 0.41-0.43.
- Expected effect: below 0.05 points (section 4).
- Validation route: blind.
- Effort: 0.5 day. Do not pursue.

**L6. kappa different from 0: the geodesic chord to the Feng-anchored endpoint, or the quadratic (not recommended)**
- Board and group: T1, every group.
- Expected effect: kappa = +1/2 is disfavoured by the in-sample readout and the within-Feng check; kappa = -1/2 removes half of a persistent regime offset if r is mostly delta. EV at or below 0.
- Validation route: blind.
- Effort: 1 day.
- Main risk: sign error on half of every gene's residual.

**L7. Own-growth (proliferation) quotas between states (prune)**
- Board and group: T1, distribution.
- Mechanism: lambda = ln(1 + c)/tau per state.
- Expected effect: share factors of about 2 from unvalidated rates; sign unknown; EV about 0.
- Validation route: counted tasks only, confounded (section 5).
- Effort: 1-2 days. Rule 2 fails before the run.

**L8. Embryo clock warp or spline interpolation of the fraction (stays pruned)**
- Board and group: embryo.
- Expected effect: no legal clock places E7.75; the cubic weights in section 3 are unstable; the Bayes weight is 1/2.
- Validation route: none (the backtests cannot count for fraction changes).
- Effort: none. Closed.

**L9. Cell-level W2 momentum / McCann extension past s = 1 (closed)**
- Board and group: T1 and heart.
- Mechanism: ill-posed at s = 2 (section 1); the per-cell latent transport variant was pruned (TRIES, JL-OT).
- Effort: none. Closed.

## Deep Read Notes
[1] Gardner & McKenzie (full): the damped trend's recurrence and error-correction forms, with h-step forecast l + (phi + ... + phi^h) b.
It is optimal for a random-coefficient model in which the trend survives each step with probability phi, so phi measures trend persistence.
On M3 the fitted parameters picked special cases most of the time (damped trend itself 43% / 28%, SES with damped drift about 24%, Holt 10% / 2%); damping was needed in 84% / 70% of series.
Lesson: damping is a persistence probability that can be estimated, not a knob.

[5] TIMO (full): SEACells metacells, then semi-relaxed OT between consecutive stages with a KL penalty on source mass; growth = transported mass / proportion; no growth prior.
It argues per-cell WOT/moscot growth is sparse with inflated variance, and shows lower within-type growth variance.
Its theory assumes oracle metacells, and it concedes that growth and differentiation cannot be separated from snapshots.
Lesson: state-level quotas are the stable unit; one post-break snapshot cannot identify rates.

[4] scFM (evaluation sections in full): "Naive replicates the last observed timepoint"; ZB, DR and SC series with 10-18 time points; Sinkhorn W2 and mean pairwise l2.
scFM and scNODE beat Naive in every column; VGFM loses to Naive in 7 of 8 extrapolation columns; PRESCIENT loses once (SC t = 5).
The authors run their own benchmark.
Lesson: no transfer to three stages with a regime break.

## Gaps
- Searched and not found:
  - a forecasting benchmark with a sampling or regime break at the last observed time, or with three or fewer time points;
  - a legal clock that places E12.5 or E7.75 inside its window (reference staging and DTW need in-window references);
  - a validation of signature-based birth rates against lineage-traced truth in heart or embryo;
  - accuracy figures for Wasserstein-spline extrapolation (only overshoot in interpolation).
- The Hyndman & Billah full text was blocked (HTTP 403, then timeouts); the half-slope drift is confirmed from the authors' own package code [3].
- Blood (24% of s1 chamber cells) has no legal direction with skill: the committed readout's JAX and Feng blood rows show weak Spearman (0.05-0.17) and are in-sample. Not proposed.
- Counter-claims against L1, the top lever:
  - The truth-fitted 30-PC basis may not carry the cycle axis strongly at E12.5, which makes the gain about 0.
  - JAX's late heart nuclei are few and look batch-dominated (digest 4.3), so the log-odds change may be technical even after depth matching.
  - Feng's within-zone change comes from whole cells.
  - Cycling nuclei carry more RNA, so the quota can move the median depth past the 5% line.
  - It couples s1 and s2 unless the captain splits it.
  - Its EV is two orders of magnitude below the needed gap.

## Compliance log
- /mnt/md0/IITM/BackUp/Home/vijayavallabh/.claude-account/projects/-mnt-md0-IITM-BackUp-Home-vijayavallabh/memory/virtual-embryo-challenge.md: L (screened overall leaderboard totals appeared; nothing copied or used).
- https://arxiv.org/html/2501.03501v2: W (its results section reports developmental findings over E6.5-E7.5, which includes the embryo validation stage; nothing copied).
- Process note, no barred class: WebFetch automatically saved a binary copy of source [1] under the session's tool-results store in ~/.claude-account. It was not opened; the copy that was read was downloaded into this task's scratchpad folder.

## END
