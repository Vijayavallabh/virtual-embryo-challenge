---
task_id: g
role: Tissue physicist and spatial statistician
status: complete
sources_found: 16
---
## Sources
Months marked * are approximate (the search result gave the year only).

[1] Durand M., Large-scale simulations of biological cell sorting driven by differential adhesion follow diffusion-limited domain coalescence regime (PLoS Comput Biol, doi 10.1371/journal.pcbi.1008576) | https://pmc.ncbi.nlm.nih.gov/articles/PMC8389523 | Source-Type: academic | Accessibility: public | As Of: 2021-08 | Authority: 8/10
[2] Allen S., Ginsbourger D., Ziegel J., Evaluating forecasts for high-impact events using transformed kernel scores | https://arxiv.org/pdf/2202.12732 | Source-Type: academic | Accessibility: public | As Of: 2022-02 | Authority: 8/10
[3] Gretton A., Borgwardt K., Rasch M., Schölkopf B., Smola A., A Kernel Two-Sample Test (JMLR 13) | https://www.jmlr.org/beta/papers/v13/gretton12a.html | Source-Type: academic | Accessibility: public | As Of: 2012-03 | Authority: 10/10
[4] Score Matched Neural Exponential Families for Likelihood-Free Inference (defines the kernel score and ties it to squared MMD) | https://arxiv.org/pdf/2012.10903 | Source-Type: academic | Accessibility: public | As Of: 2020-12 | Authority: 6/10
[5] Weighted scoringRules: Emphasising Particular Outcomes when Evaluating Probabilistic Forecasts (names the Gaussian-kernel kernel score the MMD score) | https://arxiv.org/pdf/2305.07312 | Source-Type: academic | Accessibility: public | As Of: 2023-05 | Authority: 6/10
[6] Beysens D., Forgacs G., Glazier J., Cell sorting is analogous to phase ordering in fluids (PNAS; index page) | https://www.bohrium.com/paper-details/cell-sorting-is-analogous-to-phase-ordering-in-fluids/812600463739322369-11299 | Source-Type: academic | Accessibility: public | As Of: 2000-08* | Authority: 7/10
[7] Morpheus model M0021, Differential Adhesion: Cell Sorting in Two Dimensions (the Graner-Glazier cellular Potts model) | https://morpheus.gitlab.io/model/M0021 | Source-Type: other | Accessibility: public | As Of: 1992-10* | Authority: 6/10
[8] PASTE2, Partial alignment of multislice spatially resolved transcriptomics data (Genome Research) | https://genome.cshlp.org/highwire_display/entity_view/node/1102275/full | Source-Type: academic | Accessibility: public | As Of: 2023-01* | Authority: 8/10
[9] RAFT-UP: Robust Alignment for Spatial Transcriptomics with Explicit Control of Spatial Distortion (preprint) | https://arxiv.org/abs/2603.18249 | Source-Type: academic | Accessibility: public | As Of: 2026-03 | Authority: 4/10
[10] stVCR: spatiotemporal dynamics of single cells (bioRxiv 10.1101/2024.06.02.596937) | https://www.biorxiv.org/content/10.1101/2024.06.02.596937.full.pdf | Source-Type: academic | Accessibility: public | As Of: 2024-06 | Authority: 5/10
[11] SOCS, Accurate trajectory inference in time-series spatial transcriptomics with structurally constrained optimal transport (third-party review page) | https://www.themoonlight.io/review/biorxiv/accurate-trajectory-inference-in-time-series-spatial-transcriptomics-with-structurally-constrained-optimal-transport | Source-Type: secondary-industry | Accessibility: public | As Of: 2024-01* | Authority: 3/10
[12] Spatialign documentation, About the hyperparameters (choice of k versus density) | https://spatialign.readthedocs.io/en/latest/about%20hyperparameter.html | Source-Type: other | Accessibility: public | As Of: 2024-01* | Authority: 4/10
[13] On Estimating Non-uniform Density Distributions using N Nearest Neighbors | https://arxiv.org/pdf/1301.5346 | Source-Type: academic | Accessibility: public | As Of: 2013-01 | Authority: 6/10
[14] Mémoli F., Gromov-Wasserstein distances and the metric approach to object matching (Found. Comput. Math.) | https://www.math.ucdavis.edu/~saito/data/acha.read.w12/memoli-gromov-dist.pdf | Source-Type: academic | Accessibility: public | As Of: 2011-08* | Authority: 9/10
[15] Clock and wavefront model (Cooke and Zeeman 1976; encyclopedia entry) | https://en.wikipedia.org/wiki/Clock_and_wavefront_model | Source-Type: community | Accessibility: public | As Of: 2026-11 | Authority: 5/10
[16] Kallenberg O., Limits of compound and thinned point processes (J. Appl. Probab. 12, 269-278) | https://www.cambridge.org/core/journals/journal-of-applied-probability/article/abs/limits-of-compound-and-thinned-point-processes/FB41475A110BBD65B42B8AE9376CE18B | Source-Type: academic | Accessibility: semi-public | As Of: 1975-06* | Authority: 9/10

## Findings
1. MMD is the RKHS distance between kernel mean embeddings, MMD^2 = E k(x,x') - 2 E k(x,y) + E k(y,y'); its unbiased U-statistic drops the diagonal and can be negative, and the median-distance bandwidth is only a heuristic [3].
2. Kernel scores (the energy score, the CRPS and the variogram score among them) are proper scoring rules whose divergence is a squared MMD, so the expected score is minimised by issuing the true distribution; the variogram score is proper but not strictly proper [2][4][5].
3. In large 2-D cellular Potts simulations of differential-adhesion sorting, domain size grows only as t^(1/4) (diffusion-limited coalescence after a t^(1/3) transient, with 1/5 expected in 3-D), and earlier studies reported logarithmic or other exponents, so sorting kinetics are slow and model-dependent [1][6][7].
4. Gromov-Wasserstein distances compare metric-measure spaces only up to isometry, so a cloud and its mirror image are indistinguishable to any structure-only coupling term [14].
5. Current optimal-transport methods for spatial omics address partial overlap and distortion in slice alignment (PASTE2's partial FGW, RAFT-UP) or growth and migration between snapshots (stVCR, SOCS); none is validated on a neighbourhood-profile metric of an interpolated whole tissue [8][9][10][11].
6. A fixed-k neighbour average spans the volume out to the k-th neighbour, so its smoothing scale and bias grow as local sampling density falls [13].
7. Independent p-thinning of a point process tends to a Poisson process (Rényi; Kallenberg), so the scorer's 5% truth half is close to an inhomogeneous Poisson sample of the tissue [16].
8. Practical spatial-omics guidance treats the neighbour count k as density-dependent and recommends sensitivity checks; it offers no theory of k-sensitivity [12].
9. The clock-and-wavefront model treats maturation as a wavefront moving along the body axis, a generic reason to expect spatially graded developmental progress within one embryo [15].

## First-principles analysis
Notation. We submit rows x_i (log1p CP10k over the panel) at positions u_i. The truth is a random 5% of the target's cells (a 10% working copy, then one half). The floor is copy_last's 10% working copy. skill = d_f / (d_f + d), with d = value - ceiling and d_f = floor - ceiling for lower-is-better terms. Formulas are from veckit cc3859a and `src/vec/evaluate.py`; no scorer docstring numbers are used.

**1. Separation (exact).** Fix the submitted rows and positions and permute rows over positions. de_score, de_direction, mmd_u and variogram read X only. d2_shape, occupancy_dice and scale_log_ratio read coordinates only. So the permutation moves only `neighborhood_mmd` (25 points, plus the gate through it). Any assignment lever is therefore a single-metric lever with built-in attribution: if any other metric moves, the build is wrong.

**2. What the neighbourhood term measures.** Each profile is h_i = (1/15) sum over the 15 nearest neighbours j of x_j, with self included, computed on the whole submitted cloud. The term is an unbiased MMD^2 in a 30-PC basis and bandwidth fitted on the truth's profiles, over a fixed-seed subsample of at most 2,000 profiles per side (`experiments/density/REPORT.md`). Write x(u) = xbar_r(u) + eta, the local mean over the 15-NN span r plus local scatter. Then
 h_i ~= xbar_r(u_i) + eps_i, with Var eps ~= sigma_loc^2(r)/15, and Var h ~= Var x / 15 + (14/15) Var xbar_r.
A coherent profile shift delta (stage-wide or regional: state, composition, depth) is not reduced by the averaging, but the spread is. In the wide-kernel limit, MMD^2 ~= |P_30 delta|^2 / sigma_h^2 + O(sigma^-4), where sigma_h is the truth's profile scale (far below the per-cell scale on the noise part). So to first order the neighbourhood term is a normalised pseudobulk-shift metric, much more sensitive to content shifts than mmd_u is.
 - Evidence that content dominates arrangement (repository):
   - On three pseudo-boards the floor sits 31-139 ceiling units above the ceiling.
   - Jitter of 2 nearest-neighbour spacings moves the raw value by about one ceiling unit.
   - Shuffling all expression is about 18x worse (digest, Part 1). My inference: even that would stay below such floors.
 - Assumption: the test boards' floor-to-ceiling ratios are of the same order.
 - Consequence: an assignment lever can recover only the arrangement share of our distance. At our skill of 0.52-0.64 that share is likely a minority.

**3. Span (point-process view).** The truth is a 5% independent thinning of the target, so it is roughly Poisson [16]. Its 15-NN therefore covers about 15/0.05 = 300 target cells. Ours covers 15 x N_ref / n_pred reference cells, and the floor's covers 150.
 - Under self-similar growth (G-fold more cells per domain), the spans match when n_pred is about truth_cells. That is the count rule's region; the density study's U-shaped raw curve (minimum near 1-2x) is this effect plus an overlapping-profile bias of the unbiased estimator.
 - A local *expression* shuffle at radius R ~= sqrt(r_t^2 - r_p^2) would emulate the truth's span without moving positions. Unlike positional jitter, it would cost no Dice.
 - But the span effect is about 1% of the gap in skill: at 8x the raw value rises 0.5-0.9 ceiling units, against floors of 31-139 (`experiments/density/REPORT.md`). Counts are closed, and the remainder is below rule 5's 0.3 points.

**4. Sorting physics.** Differential-adhesion coarsening grows domain length as t^n, with n about 1/4 in 2-D and about 1/5 in 3-D [1]. Purity therefore changes slowly, and the exponent is model-dependent. Proliferation changes the cells per domain by G, which is the span effect of item 3 and is diluted in the same way. A Potts-type sorting extrapolation would need a map from Monte Carlo steps to days and an adhesion matrix. The heart ladder cannot identify either, because scope and depth change between stages (`configs/pseudo_tasks.yaml` notes). No lever.

**5. Scoring-rule view: the mixture principle and its limits.** mmd_u, variogram and the neighbourhood term are kernel-score divergences [2][3]: D(P,Q) = ||Phi(P) - Phi(Q)||^2, with Phi linear in P. Under uncertainty about Q, E D(P,Q) = ||Phi(P) - E Phi(Q)||^2 + Var Phi(Q). For a fixed kernel this is minimised by the predictive mixture P* = E Q, exactly. Two scorer features break exactness:
 - the PCA basis and the bandwidth are refitted on each truth;
 - skill = d_f/(d_f + d) is convex in d, and d_f changes with the scenario.

I tested this in a synthetic toy that uses no competition data. It is a 2-D lattice tissue with 6 domain types and Poisson counts; the bracket depths 62 and 212 are the released medians recorded in TRIES.md, and the neighbourhood formula is re-implemented. Scripts: `/tmp/claude-1001/-mnt-md0-IITM-BackUp-Home-vijayavallabh/07d440b9-3ae8-4421-8a08-a046b42f307d/scratchpad/toy_g/` (toy_embryo.py, toy_oracle.py, toy_density.py, toy_pair.py; each runs in under 60 s). Points below are neighbourhood points out of 25.
 - At matched density, exactly intermediate cells against the whole-cell bracket mixture score 20.41 against 20.25 (mean over truth progress 0.2-0.8). This is a wash: the mixture wins on the early side and loses on the late side.
 - Spatially patchy bracket identity loses 0.6-1.5 points at every clock offset. Its neighbourhood after-fraction SD is 0.29-0.40, against the binomial 0.13.
 - Anti-clumped (balanced) assignment stays within +-0.15 of random.
 - Spatially coherent progress among intermediate cells gives no gain over per-cell spread (one seed).
 - Verdict: there is no within-submission spread or arrangement lever on this term.

**6. Floor-normalisation asymmetry (the one positive result).** d_f is the copy_last floor's distance to the realised truth. It shrinks when the truth sits near copy_last, so skill falls fastest there for any prediction that is not a copy. A closed form, for a fixed kernel, whole-cell bracket mixtures and a negligible ceiling: Phi(mix_a) = Phi_0 + a(Phi_1 - Phi_0). For a prediction at fraction t, a truth at t', and copy_last the earlier bracket (t' = 0):
 s(t; t') = t'^2 / (t'^2 + (t - t')^2).
 - At t = 1/2 this gives s = 0.31 at t' = 0.2, against 0.88 at t' = 0.8.
 - Add a second upload at t2 = t/2 = 1/4. Under a uniform prior on t' in [0.2, 0.8], E[max] - E[s1] = +0.13 skill per mixture-sensitive group, or +3.3 points per 25-point group in this idealisation.
 - The toy (refitted kernel, a ceiling, intermediate-cell truths, 3 seeds) agrees. An s2 at 0.35 adds +1.2 to +1.7 neighbourhood points and +1.5 to +1.9 cell-MMD points; an s2 at 0.65 adds only +0.7 to +1.1. Shading s1 itself gains nothing on average, consistent with the pruned clock warp.
 - Real-data instance (TRIES headroom table): on E7.25, whose truth sits near its copy_last E6.75, fgw's cell state is below the floor (-2.22). On E7.5, whose truth sits past nominal and away from copy_last, it is +4.21. Both pseudo-tasks use the earlier bracket as reference (`configs/pseudo_tasks.yaml`).

The expression group is nearly invariant to the mixture fraction. If the test reference is the earlier bracket, a mixture's predicted logFC is t x (pb_after - pb_before) plus reference-copy noise. That is a positive rescaling. de_score is rank-based and de_direction is a correlation after partialling out level, so both are invariant to it to first order; only the signal-to-noise falls slowly as t shrinks. A clock hedge therefore trades only cell_state and the neighbourhood term (and the gate).

The current embryo pair does not cover this. s1 (fgw) and s2 (fgw_shape0) share one mixture and differ only in geometry (`docs/S2_RULES.md`, T2 embryo). So the steep scenario, E7.75 near E7.5, is unhedged on 50 of the board's 100 points.

Transfer check (not my boards; one line each): T3's s2 (no factor) already covers its wild-type side. The heart s1 is close to a copy. T1's truth is probably far from copy_last because of the scope change, which is the gentle regime. The embryo is the board whose pair leaves its steep side open.

**7. Hygiene.** The scored subsample of profiles is a fixed-seed draw over submitted rows, so row order must stay random. It does: `sample_cells` and `interpolate_clouds` both permute. Exploiting the order would be gaming the scorer and is not proposed. If a build repeats pool rows, the jittered copies (at most 0.25 median NN apart) sit inside each other's 15-NN and overweight self.

**Bottom line for the top-3 objective.** Each T2 board counts 0.5 in Overall. This lane's realistic ceiling is about +1 Overall, from L1. The 15-25 points the team needs must come mostly from T1 and T3. The floor-normalisation asymmetry is the most transferable result: a pair's s2 should cover the "truth near the baseline" scenario wherever s1 is bold.

## Candidate levers
Ranked by expected value to the top-3 objective. Only L1 is a positive bet; L2 is a cheap check; L3-L8 are pruned or closed, recorded here to save tries.

**L1. Embryo s2 copy-side clock hedge (an "E7.5-leaning" s2).** Captain only.
- Board and group: T2:embryo:test_interp. It moves cell_state (0.25) and local_organisation (0.25), and protects the gate. expression_change is unchanged to first order (section 6).
- Mechanism: keep s1 (fgw at t = 1/2). Build s2 as fgw_shape0 (E7.5 geometry, the same coupling, size rule and count) with its whole-cell bracket mixture at t2 = t/2 = 1/4. The rule is declared before any run and is untuned: halfway between copy_last and nominal. Under better-of-two, E[max(S1,S2)] = E[S1] + E[(S2 - S1)+]. s2 covers the scenario that the floor normalisation makes costliest for s1: E7.75 near copy_last.
- Expected effect:
  - If that scenario happens: about +3 to +6 embryo points (toy, plus the E7.25 instance). The prior is about 0.3-0.5, from 1 of 2 released transitions.
  - Otherwise s1 wins, and part of the current s2's geometry-only edge is lost (about +1.5 shape points on E7.5).
  - Net: +0.5 to +2.5 embryo points, central +1.2, so +0.25 to +1.25 Overall (central about +0.6). Break-even prior on the steep scenario is about 0.15-0.2.
  - Uncertainty is wide: the magnitudes come from a toy and one real instance.
- Legality check:
  - Data: no new data; it uses the test inputs (E7.5, E8.0) and the existing pipeline.
  - Closed levers: no change to cell numbers, sizes or scope; not T3 or T1 composition; not T1 depth; no cross-modal input.
  - Reads and inference: it reads no target or score and infers nothing from scores. Research rules A, B, W, L and G: none touched.
  - Rule 2: its E7.25 caution concerns a single submission's backtest gain. L1 claims no backtest gain and uses E7.25 only as evidence that the steep scenario occurs.
  - Frozen rule: it changes the frozen embryo s2 rule, so it needs the captain's ruling. If s1 changes, its s2 goes back to him anyway.
- Validation route: captain's blind bet only. TRIES records that no counting backtest can judge a mixing-fraction change. Target-free support before the ruling, which is not a try:
  - build pseudo-truths from the test pair only, as whole-cell mixtures of E7.5 and E8.0 at t' = 0.1-0.9 on fgw positions;
  - score s1 and the declared s2 against each in the emulator;
  - report the curve's asymmetry and E[max] - E[s1] under a prior declared beforehand.
  The table checks the direction only; it never tunes t2.
- Effort: about 1 day (a mixture-fraction argument for fgw_shape0 if it lacks one, the target-free table, a ruling memo).
- Main risk: E7.75 sits at or past nominal, as E7.5 did. s2 then adds nothing and the geometry hedge is diluted. If the test reference is E7.25 rather than E7.5, the expression group is no longer invariant to t2.

**L2. Neighbourhood hygiene audit (target-free, code-only).**
- Board and group: both T2 boards; local_organisation.
- Mechanism: three checks, none of which prints a cell count.
  - (a) Do any submitted rows repeat pool rows? If so, the copies sit inside each other's 15-NN. Report a yes/no only.
  - (b) Is the row order random? (It is today.)
  - (c) Embryo: compute each cell's after-bracket fraction f_i among its 15-NN, and compare SD(f) with sqrt(t(1-t)/15) and with a within-type permutation null. The merged cloud carries two section-plane families (transverse E7.5, sagittal E8.0; `experiments/coords/REPORT.md`), which could segregate the brackets.
- If (c) shows more than about 1.5-2x the null, re-randomise bracket identity locally within coarse type. Only the neighbourhood term moves, and the toy says excess segregation costs 0.6-1.5 points.
- Expected effect: 0 to +0.3 board points if clean; up to about +1 embryo point if (c) finds segregation.
- Legality: no data beyond the build inputs and no target is read, so the audit is not a try. A fix would be a single-metric try on the two embryo backtests, with high power because every other metric is bit-identical.
- Validation route: target-free diagnostic. Any fix goes to the counted embryo backtests under the adoption rule.
- Effort: 0.5 day.
- Main risk: it finds nothing, the likely outcome.

**L3. Heart within-submission momentum mixture (pruned, rule 5).**
- Board and group: heart; cell_state and local_organisation.
- Mechanism: half the cells at no momentum and half at twice B0's, instead of all at B0's. This is the mixture principle of section 5.
- Expected effect: the toy shows a wash at matched density, and the heart's momentum direction is weak. About 0 +- 0.5, so below 0.3.
- Legality: open in form, needs no new data.
- Validation route: the two heart backtests.
- Effort: 1 day.
- Main risk: extra spread along a direction that carries no signal.

**L4. Embryo anti-clumped (balanced) bracket assignment (pruned).**
- Board and group: embryo; local_organisation only.
- Expected effect: within +-0.15 points in the toy. Rule 5 prunes it.
- Validation route: the embryo backtests, if it were ever run.

**L5. Embryo clock-field patches (do not).** Spatially coherent bracket identity costs 0.6-1.5 points at every offset in the toy. Its two reasons, the mixture principle and spatially graded maturation [15], fail against the refitted kernel and the floor normalisation.

**L6. Span-matched local expression shuffle (pruned).** Heart and embryo; local_organisation. Worth about 1% of the gap (density study). The count it would emulate is closed.

**L7. Differential-adhesion (Potts) sorting extrapolation on the heart (pruned; a blind bet at best).** Kinetics are slow and model-dependent [1]. Its time map and adhesion matrix cannot be identified on the ladder.

**L8. Choosing which rows the scorer's fixed profile subsample reads (closed, integrity).** Not proposed. Keep the row order random.

Closed lines not reopened: depth, blend, herding and coupling. The partial, unbalanced and entropic OT literature [8][9][10][11] addresses overlap and growth. The embryo's measured coupling defect is mirror symmetry, which isometry invariance makes invisible to structure-only terms [14], and E11 already tested a position-aware fix. Clock warp for s1 also stays closed: shading s1 gains nothing on average (section 6).

## Deep Read Notes
[1] Durand 2021 (PMC full text). A non-fragmenting cellular Potts algorithm simulates up to 320k cells. The mismatched-contact length decays as t^(-1/4) whatever the mixture ratio, with domain size growing as t^(1/4) after a t^(1/3) transient of entangled clusters. Diffusion-limited coalescence gives n = 1/(2 - alpha) with alpha = -1, so 1/4 in 2-D and 1/5 in 3-D. Earlier work found logarithmic decay (Graner-Glazier) or 1/3 (Nakajima), and experiments 0.49-0.74. Implication: sorting-driven purity change is slow and not calibratable over 2 days.
[2] Allen, Ginsbourger and Ziegel. Read in the arXiv HTML rendering at https://ar5iv.labs.arxiv.org/html/2202.12732 (about the first 100k of 157k characters). The kernel score is S(P,y) = E rho(X,y) - (1/2) E rho(X,X') - (1/2) rho(y,y) for a conditionally negative definite rho, and its divergence is a squared MMD. It is proper, and strictly proper with a characteristic kernel. The energy score, the CRPS and the variogram score are kernel scores; the variogram score is not strictly proper.
[3] Gretton et al. Read in the HTML rendering of the precursor report, https://ar5iv.labs.arxiv.org/html/0805.2368 (about the first 100k of 156k characters). Lemma 4: MMD = ||mu_p - mu_q||_H. Lemma 5: the expansion and the unbiased MMD^2_u (i != j, m(m-1)), which may be negative. Theorem 3: with a universal kernel, MMD = 0 if and only if p = q. Section 8: the median-distance bandwidth is a heuristic.

## Gaps
Searched for and not found:
- a theory of how kNN-profile MMD responds to sampling density or thinning (only practice guidance [12] and generic kNN smoothing [13]);
- any treatment of skill scores whose denominator is a baseline's scenario-dependent distance, or of better-of-two submission design, in the forecast-verification literature;
- a published whole-tissue interpolation validated on neighbourhood metrics;
- the full text of the Beysens-Forgacs-Glazier exponent [6].
In-vivo coarsening data for the heart were not sought: they would fall in the protected window.

Counter-claims against L1:
1. The asymmetry rests on a toy, a closed form and one real instance (E7.25). If E7.75 behaves like E7.5 did (past nominal), s2 adds nothing and its current geometry-only edge is diluted.
2. If the test board's reference is E7.25 rather than E7.5, d_f is larger in every scenario and the expression group gains a t2-dependent component, so the trade-off changes.
3. The real curvature depends on how E7.75's depth sits between 62 and 212. MERFISH depth behaved unlike binomial thinning (E7, H3), so the toy's magnitudes may overstate the effect.
4. E7.25's steep outcome is a known feature of one backtest. Treating it as a prior rather than a fit is a judgement for the captain.
5. The convex skill map concentrates s2's value in tail scenarios. If the steep-scenario prior is below about 0.15-0.2, the expected gain is near zero or negative against the current s2.

## Compliance log
none

Note (no barred class): two PDF fetches (https://arxiv.org/pdf/2202.12732 and the CMU-hosted copy of [3]) were saved automatically by the fetch tool under the harness's tool-results folder. Neither was opened; both texts were read from the HTML renderings instead.
## END
