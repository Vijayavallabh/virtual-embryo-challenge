---
task_id: b
role: Statistical physicist and systems biologist (fluctuation-response relations, causal inference from natural variation)
status: complete
sources_found: 15
---
## Sources
[1] Fluctuation structure predicts genome-wide perturbation outcomes (CIPHER; preprint deposited in PMC, PMID 40831507) | https://pmc.ncbi.nlm.nih.gov/articles/PMC12363937/ | Source-Type: academic | Accessibility: public | As Of: 2025-06 | Authority: 6/10
[2] Co-fluctuations of genes in unperturbed single cells predict transcriptome-wide outcomes to genetic perturbations (APS March Meeting abstract, same group as [1]) | https://meetings-archive.aps.org/smt/2026/mar-f60/4/ | Source-Type: academic | Accessibility: public | As Of: 2026-03 | Authority: 4/10
[3] Kaneko, Relationship among phenotypic plasticity, genetic and epigenetic fluctuations, robustness, and evolvability (arXiv 0807.3064; re-derives the Sato et al. relation) | https://ar5iv.labs.arxiv.org/html/0807.3064 | Source-Type: academic | Accessibility: public | As Of: 2008-07 | Authority: 7/10
[4] Sato, Ito, Yomo, Kaneko, On the relation between fluctuation and response in biological systems (PNAS 2003) | https://pmc.ncbi.nlm.nih.gov/articles/PMC283550 | Source-Type: academic | Accessibility: public | As Of: 2003-11 | Authority: 8/10
[5] Lehner and Kaneko, Fluctuation and response in biology (Cell Mol Life Sci 2011) | https://pmc.ncbi.nlm.nih.gov/articles/PMC3043249 | Source-Type: academic | Accessibility: public | As Of: 2011-03 | Authority: 7/10
[6] Furusawa and Kaneko, Formation of dominant mode by evolution in biological systems (arXiv 1704.01751; Phys Rev E 2018) | https://arxiv.org/pdf/1704.01751 | Source-Type: academic | Accessibility: public | As Of: 2017-04 | Authority: 7/10
[7] Gupta et al., Inferring gene regulation from stochastic transcriptional variation across single cells at steady state (PNAS 2022, doi 10.1073/pnas.2207392119; Europe PMC full text) | https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9407670/fullTextXML | Source-Type: academic | Accessibility: public | As Of: 2022-08 | Authority: 8/10
[8] Maathuis, Colombo, Kalisch, Buhlmann, Predicting causal effects in large-scale systems from observational data (Nat Methods 2010, doi 10.1038/nmeth0410-247) | https://stat.ethz.ch/Manuscripts/buhlmann/maathuisetal2010.pdf | Source-Type: academic | Accessibility: public | As Of: 2010-04 | Authority: 8/10
[9] Chevalley et al., CausalBench: a large-scale benchmark for network inference from single-cell perturbation data | https://arxiv.org/abs/2210.17283 | Source-Type: academic | Accessibility: public | As Of: 2022-10 | Authority: 7/10
[10] Aslyamov, Ptaszynski, Esposito, Macroscopic fluctuation-response theory and its use for gene regulatory networks (arXiv 2510.17587; Phys Rev Lett 136, 067102) | https://arxiv.org/abs/2510.17587 | Source-Type: academic | Accessibility: public | As Of: 2025-10 | Authority: 7/10
[11] Lehner, Genes confer similar robustness to environmental, stochastic, and genetic perturbations in yeast (PLoS One 2010, doi 10.1371/journal.pone.0009035) | https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2815791/ | Source-Type: academic | Accessibility: public | As Of: 2010-02 | Authority: 7/10
[12] Mitchel et al., Impact of segmentation errors in analysis of spatial transcriptomics data (bioRxiv; Nat Genet 2026, doi 10.1038/s41588-025-02497-4) | https://www.biorxiv.org/content/10.1101/2025.01.02.631135v1.full | Source-Type: academic | Accessibility: public | As Of: 2025-01 | Authority: 8/10
[13] Systema: a framework for evaluating genetic perturbation response prediction (Nat Biotech 2025; Brbic lab project page) | https://brbiclab.epfl.ch/projects/systema | Source-Type: academic | Accessibility: public | As Of: 2025-09 | Authority: 7/10
[14] Hilfinger, Norman, Paulsson, Exploiting natural fluctuations to identify kinetic mechanisms in sparsely characterized systems (Cell Systems 2016; lab publication list, not read) | https://paulsson.med.harvard.edu/papers | Source-Type: academic | Accessibility: public | As Of: 2016-04 | Authority: 6/10
[15] Evaluating parameter uncertainty in the Poisson lognormal model with corrected variational estimators (arXiv 2411.08524) | https://arxiv.org/html/2411.08524v1 | Source-Type: academic | Accessibility: public | As Of: 2024-11 | Authority: 5/10

## Findings
1. The Sato-Kaneko relation, (<x>_{a+da} - <x>_a)/da = C <dx^2>, is derived for a near-Gaussian distribution in which the parameter enters only through a linear coupling C(a - a0)(x - X0) with C and the variance held constant, and its authors note that C may itself change when the parameter change is large [3][4].
2. Furusawa-Kaneko theory, with E. coli laboratory-evolution data, constrains high-dimensional expression changes after evolution to one major axis that correlates with growth rate, so responses to diverse perturbations share a dominant low-dimensional mode [6].
3. CIPHER writes the multivariate relation as an exponential tilt of the unperturbed distribution (dX = Sigma u) and reports that control-cell covariance explains CRISPRi/a responses across ten Perturb-seq datasets far better than shuffled or variance-only nulls, but it fits the perturbation strength u to each observed response, computes covariance from raw counts with no normalisation described, and finds responses propagating along about three eigengenes (mean participation ratio 3.07) tied mostly to translation and homeostasis [1][2].
4. Even with perfect detection, a transcription factor's bursts produce target co-fluctuations that peak near Spearman 0.2 (Hill coefficient 2; 0.13-0.33 across Hill 1-4), and at 10% capture the time-shifted covariation drops to about 0.04-0.08, so roughly 50,000 cells per time point would be needed, while stable cell states add a time-invariant covariation that swamps the regulatory signal [7].
5. In K562 cells, the genes correlated with GATA1 at all six sampled time points were 3.8-fold enriched among genes changed by GATA1 knockdown, the clearest single-cell case found here of steady-state co-expression marking knockdown-responsive genes [7].
6. Observational-only causal-effect estimation (IDA, 63 wild-type yeast profiles over 5,361 genes) predicted the strongest single-gene deletion effects better than Lasso or elastic-net as a proof of concept, yet on genome-wide Perturb-seq data network-inference methods given interventional data did not beat observational-only ones, and causal methods did not beat non-causal ones [8][9].
7. Systema shows that common perturbation-prediction metrics reward systematic variation shared across perturbations, that simple baselines match published models, and that unseen perturbations are much harder to predict than those metrics suggest [13].
8. For nonequilibrium steady states, exact fluctuation-response relations tie the linear response to the power spectral density (the time correlations) of stationary fluctuations and allow reconstruction of the linearised kernel and diffusion matrix, whereas the classical equal-time form needs detailed balance [10].
9. Across 4,656 yeast gene deletions, genes that buffer environmental or stochastic variation also buffer genetic change with similar efficacy (Spearman 0.39), the main empirical support for coupling natural variability to perturbation sensitivity [11][5].
10. In imaging spatial transcriptomics including MERFISH, transcripts mis-assigned from neighbouring cells create spurious co-expression and can dominate top differential-expression and neighbour-influence hits (for example, 11 of 36 top up-regulated genes in one fibroblast comparison were neighbouring malignant-cell markers) [12].

## First-principles analysis
Notation: x is one cell's 500-gene panel expression, t a harmonised label, mu_t and Sigma_t the within-label mean and covariance in the WT reference, a the activity of the deleted regulator, and dp, dt the scorer's predicted and true pseudobulk deltas against WT. Assumptions are marked [As1]-[As6].

**1. What a wild-type covariance can predict (theory).**
- Exponential tilt, the multivariate Sato-Kaneko/CIPHER form [1][3][4]. If a perturbation reweights the WT ensemble, p_u(x) proportional to p_0(x) exp(u.x), then mu(u) = grad log Z(u) and d mu = Sigma u + (1/2) k3[u,u] + ... A knockout of regulator k is a tilt along e_k. Fixing u by the regulator's own loss (d mu_k = -mu_k) gives d mu_j = -mu_k Sigma_jk / Sigma_kk = beta_{j|k} d mu_k: the response is the regression on the regulator, proportional to its covariance column. [As1: the knockout reaches other genes only through a, within each label. As3: the WT file is a stationary ensemble of the same system.]
- Large perturbation. As u goes to minus infinity the exact tilt concentrates on the cells with the lowest x_k, so the non-linear tilt prediction becomes E_0[x_j | x_k minimal], the "natural knockout" cells. That is conditioning, not intervention. E[x_j | do(x_k = 0)] equals E[x_j | x_k = 0] only if nothing else in the stratum drives both genes and x_j does not feed back on x_k (the back-door condition with an empty or measured adjustment set). Sato-Kaneko make the same demand in one variable (a linear coupling with constant C and unchanged variance) and say C may change for large shifts [3]. A knockout of a lineage regulator is the large case; near a bifurcation it removes an attractor, which shows up as composition, not as a linear within-state shift (CIPHER's synthetic error also peaks near transitions [1]).
- Out of equilibrium. Linearise dx/dt = J(x - mu) + h + xi with <xi xi^T> = 2D. The static response is d mu = -J^{-1} h, and the equal-time covariance solves J Sigma + Sigma J^T + 2D = 0. Only when J is symmetric and D = d I (detailed balance, isotropic noise) is Sigma = -d J^{-1}, so that d mu = Sigma h / d (the FDT form). In general -J^{-1} = [integral over tau of C(tau)] Sigma^{-1} with C(tau) = e^{J tau} Sigma: the static susceptibility needs time-lagged correlations [10][14], which a snapshot lacks. Gene networks are driven and non-symmetric, so this is the general case.
- Timescales. In a cascade with regulator decay g_k, target decay g_j and gain b, the steady-state gain is b/g_j, but the snapshot slope is Cov/Var = b/(g_j + g_k) (from the Lyapunov equation). The snapshot under-reads the gain by g_j/(g_j + g_k), most when the regulator fluctuates faster than the target relaxes; Gupta et al. measure this lag and its small amplitude [7]. A co-activator regulated at protein level, as the test gene is, adds a further layer the mRNA snapshot cannot see.
- Where the relation does work. CIPHER fits u per perturbation, and its responses run along about three global eigengenes [1]. That is the systematic part of responses which Systema warns inflates metrics [13], the dominant mode of Kaneko's theory [6], and the variability-sensitivity coupling in Lehner's yeast data [11][5]. So a WT covariance predicts the subspace and the gene-wise susceptibility of a response, and its direction once one scalar per perturbation is fitted. Zero-shot sign and size for one specific knockout need exogenous variation of the regulator within the stratum, which a snapshot does not supply.

**2. Within-type against across-type covariation, and what MERFISH adds.**
- Law of total covariance: Sigma_total = sum_t w_t Sigma_t + sum_t w_t (mu_t - mu_bar)(mu_t - mu_bar)^T. Fluctuation-response on Sigma_total adds the across-label term, which predicts the pseudobulk change of a composition shift towards labels where the regulator is low. That is right for the wrong reason when a knockout shrinks its own compartment, and it is a closed lever on T3, so only within-label Sigma_t may be used.
- [As2] Within a label, Sigma_t = Sigma_reg + Sigma_state + Sigma_volume + Sigma_spill + diag(Poisson). Given the rates, Poisson sampling is independent across genes, so in size-factor-normalised counts it only inflates the diagonal. Off-diagonal covariances are then nearly unbiased, and a slope is attenuated only through the regulator's own variance. A method-of-moments correction repairs that: Var_signal(x_k) = Var(x_k) - E[x_k / s]. A Poisson-lognormal fit is the heavier, full-matrix alternative [15].
- Normalisation does not remove the rest. Sigma_volume is the captured cell volume and slice geometry of thin sections. Sigma_spill is segmentation admixture, which can dominate top hits in MERFISH [12]. Sigma_state (sub-states, position along gradients) is the dominant time-invariant covariation [7]. A snapshot cannot separate Sigma_state from regulation without over-control: conditioning on sub-state also blocks a knockout effect that acts through sub-state.
- Synthetic check: numpy only, no competition data; scripts in /tmp/claude-1001/-mnt-md0-IITM-BackUp-Home-vijayavallabh/07d440b9-3ae8-4421-8a08-a046b42f307d/scratchpad/frr_sim/ (sim_attenuation.py, sim_deatten_confound.py). [As5] The setup:
  - about 100 other counts per cell (the WT heart series' medians are 130 at E8.75 and 75 at E9.5, TRIES row H3);
  - a regulator at 0.2, 1 and 3 mean counts per cell (detected in 16%, 53% and 82% of cells; the released training knockout's deleted gene is detected in 13-23% of WT cells, digest 2.4);
  - a true relative knockout effect of 0.6.
- Results:
  - The naive within-label slope recovers about 7%, 26% and 49% of the effect at the three regulator levels.
  - The moment-corrected slope recovers it (0.57) at 1 count per cell or more, and is unstable at 0.2 (reliability 0.07).
  - A shared state axis of log-SD 0.4 inflates the corrected slope to 0.92-0.94 (1.13 at the sparsest level). With no causal effect at all, the same axis gives a spurious relative effect of about 0.2.
- At MERFISH depth the estimator is therefore either attenuated or confounded, and the confounded part is as large as a plausible real effect.

**3. Where the T3 points are (committed results; no new run).**
- Decomposition of the released knockouts' true delta, as share of dt (`experiments/t3_decomp/results/*/*_decomposition.csv`):

  | Component | Gata4 rep1 | Gata4 rep2 | Mab21l2 |
  |---|---|---|---|
  | Shared non-Mesp1 factor | 0.46 | 0.38 | 0.60 |
  | Depth | 0.34 | 0.40 | 0.07 |
  | Mesp1-lineage composition | 0.16 | 0.11 | 0.07 |
  | Stage | 0.07 | 0.12 | 0.15 |
  | Coverage | 0.01 | 0.14 | 0.01 |
  | Within-type residual | -0.04 | -0.15 | 0.10 (Spearman 0.18, DE sign agreement 0.55, chance level) |

- A WT-covariance predictor targets only the residual: the knockout-specific, type-specific within-label change. Its oracle ceiling, on top of an oracle shared factor at WT depth (`oracle:ko_within_wtdepth` minus `oracle:factor_all`, in-sample):
  - Gata4 rep1: 79.99 - 75.52, about +4.5 T3 points;
  - Gata4 rep2: 74.19 - 69.81, about +4.4;
  - Mab21l2: 89.22 - 85.39, about +3.8.
- Cross-knockout reality:
  - K2a (Gata4 factor on the Mab21l2 fold) scored 50.08 against 49.90 for its permuted control.
  - D11b's 59.40 / 62.15 is mostly the WT E9.5-against-E8.75 axis (`experiments/gate/results/t3_swap_review_2026-11-05.md`).
  - s2's stand-in sits at the severity floor in 14 of 20 seeds (the same review).
- Scorer mechanics for a small added component (pinned `common/core_metrics.py`; digest Part 1):
  - de_score and de_direction are rank-based, so magnitude does not matter to them.
  - Severity regresses dp on dt through the origin over truth-DE genes, with R^2 = 1 - ||dp - beta dt||^2 / ||dp - mean(dp)||^2. It sits at the floor when beta <= 0 or R^2 is under the small gate the digest quotes.
  - So a common offset in dp across DE genes lowers R^2: the residual keeps the offset, while the centred denominator drops it. A component must therefore add no global shift.
  - A truth down-gene needs WT pseudobulk of at least 0.25, so predicted decreases below that waste down slots.
- Synthetic scorer check (`.../frr_sim/sim_scorer_gain.py`; 500 genes, truth-DE = the top 20% by |dt|). Setup: a knockout-specific module carries 5-15% of dt on 10% of genes and is predicted with correlation 0.2-0.4. Results:
  - Raw top-n overlap rises by +0.001 to +0.010, and Spearman by +0.002 to +0.015.
  - When the base prediction has no relation to the truth, the module lifts severity off the floor in some draws (+0.02 to +0.5 mean log-beta units).
  - When the base tracks the shared part (correlation 0.6), the severity gain is only +0.002 to +0.05.
  - [As6] With local skill slopes near the floor, that is +0.05 to +0.7 T3 points, centre about +0.2. That is under the adoption rule's +0.3 on every fold except in the optimistic corner. A module with no real signal scores about 0 or slightly negative.

**4. A zero-shot predictor for our T3 setting (design).**
- Inputs: the fold's WT reference (WT E9.5 for Mab21l2; WT E8.75 for Gata4 and for the test), deduplicated, with harmonised labels and the Mesp1-lineage flag (`configs/lineage_map.yaml`; the rules.yaml non-Mesp1 and mixed lists). [As4: deletion is complete in Mesp1-lineage labels and absent elsewhere.] No parameter is read from any knockout file.
- Steps:
  1. Recover counts c from log1p CP10k; size factor s_c = N_c / median N.
  2. Within each label, regress out a fixed list of nuisances: log N_c; the label composition of the cell's 15 spatial nearest neighbours (spillover); z (section plane).
  3. Regulator activity a_c. On the folds it is the deleted gene's normalised counts averaged over the cell's 15 same-label spatial neighbours: this raises reliability, and the positional confounding it adds is controlled below. On the test it is the motif proxy of section 6.
  4. Slope per gene and label: beta_{g,t} = Cov_t(x_g, a) / Var_signal,t(a) on nuisance residuals.
     - Shrink it by James-Stein: kappa_t = max(0, 1 - Var_noise / Var_g(beta_hat)), with Var_noise = (1/4) Var(beta_A - beta_B) from two spatial halves of the label.
     - Both quantities come from the WT file only.
  5. Magnitude, fixed a priori: complete loss, d a_t = -E_t[a], in Mesp1-lineage labels; d a_t = 0 in strict non-Mesp1 and mixed labels.
  6. Factor phi_{g,t} = clip(1 + kappa_t beta_{g,t} d a_t / mu_{g,t}, 0.25, 4), using the rules.yaml clip frozen on 2026-10-18. phi = 1 for:
     - labels with fewer than 200 WT cells;
     - predicted decreases of genes whose WT pseudobulk is below 0.25.
     The deleted gene itself is treated identically in both arms.
  7. Realisation: the T3 script's count-space primitive (binomial thinning below 1, Poisson top-up from kNN-smoothed profiles above 1, depth kept, renormalised to CP10k). The renormalisation keeps the component free of a global offset.
- One non-linear variant, to be registered instead of the linear one, not as well: the tilt-limit contrast. It estimates beta_{g,t} as (E_t[x_g | a in its lowest decile] - E_t[x_g]) / (E_t[a | lowest decile] - E_t[a]) and extrapolates to complete loss as in step 5.
- Why every constant is a priori: complete loss is the genotype; the moment correction and the shrinkage are statistics of the WT file; the clip and the 0.25 floor are an existing rule and the scorer's threshold. Nothing is swept.

**5. What would count as evidence.**
- Folds: F1, the Mab21l2 knockout against WT E9.5; F2 and F3, the Gata4 replicates against WT E8.75. With nothing read from a knockout file, each is a zero-shot cross-knockout fold for this component.
- Arms: base plus component against base. Two bases, fixed in the TRIES row:
  - a WT copy (fac:none);
  - a cross-knockout factor base: K2a on F1, transfer:mab21l2_all on F2 and F3. This measures the increment where the test's s1 would sit.
- Control: the same component with beta permuted across genes within WT-expression deciles (one permutation, seed fixed in the row).
- Predicted metrics: de_score and de_direction first, severity second. A severity rescue counts only if the arm leaves the floor with beta > 0.
- F1 is the weak fold: its deleted gene is sparse and its knockout has no heart defect at its stage. A fail on F1 alone is weak evidence.

**6. The off-panel regulator.**
- Ctnnb1 is off-panel, and its activity is set post-translationally, so its mRNA would be a poor regressor anyway (digest 2.4).
- Proxy: a_c = the mean, over a motif-selected subset M of panel genes, of their within-label z-scored Pearson residuals, smoothed over same-label neighbours. Code computes M as the panel genes whose TSS-centred sequence window (fixed width) scores above a fixed PWM threshold for the TCF/LEF family.
- Resources and provenance (no lists):
  - JASPAR CORE vertebrate profiles (HT-SELEX and ChIP-seq binding data);
  - HOCOMOCO v12 (ChIP-seq and HT-SELEX binding data);
  - the mouse reference genome, with Ensembl or GENCODE TSS annotation;
  - optionally, the cisTarget motif-ranking databases built from TSS-centred sequence only.
- Stay sequence-only. Do not use epigenomic region sets that pool protected or validation stages (for example ENCODE's mouse embryonic heart biosamples).
- Barred or avoid:
  - MSigDB's Wnt/beta-catenin hallmark (built from beta-catenin activation experiments);
  - PROGENy footprints (perturbation-derived);
  - DoRothEA or CollecTRI regulons for TCF/LEF or CTNNB1 (curated from literature that includes perturbation results; the digest's screened note says such edges restate held-out-genotype results);
  - GO annotations with IMP, IGI or IEP evidence codes;
  - any expression signature of Wnt inhibition.
- KEGG and Reactome membership lists pathway components, not readouts, and adds little.
- Each resource needs an `EXTERNAL_DATA.md` row and `docs/PLAN.md` 4.4 provenance before use (`docs/DATA_POLICY.md` 2026-11-05 (b)(iii)). That policy line also names A15. If the captain reads it as gating any pathway or motif set until the organisers answer, L1's test stage and Stage 2 are closed until then; the on-panel folds (Stage 1) need no external set.
- The sign is generic biology: losing a co-activator lowers its targets' transcription. Gene directions come out of code as beta x d a, never from a list.
- Analogue for validation: a GATA-family PWM proxy on F2 and F3, with Gata4 and Gata6 excluded from M, compared with the direct-regulator version.
- Caveats:
  - TCF/LEF and GATA motifs are short and common, and promoter windows miss distal enhancers.
  - M is circular: its genes are predicted down because they define a. The proxy predicts only that the motif module goes down, plus its within-label co-variates.

**7. Failure modes, and what success on Mab21l2 and Gata4 would and would not imply.**
- Failure modes:
  - (a) state or positional confounding;
  - (b) admixture, where the regulator's module includes neighbouring labels' markers [12];
  - (c) a captured-volume mode;
  - (d) attenuation at sparse detection, or overshoot after correction;
  - (e) compensation and feedback: a paralog or feedback gene can co-vary positively in WT and rise after loss;
  - (f) the knockout acting through composition (a closed lever) or at protein level (adhesion), which mRNA covariance cannot see;
  - (g) the shared run factor plus depth carrying most of dt (0.78-0.80 of the decomposed delta on Gata4), which no WT covariance predicts;
  - (h) Cre timing and protein perdurance making the loss partial at the scored stage;
  - (i) the proxy tracking another positional programme.
- What success on F1-F3 would show: WT within-label covariance with an on-panel regulator, after depth and spillover control and at an a-priori magnitude, carries out-of-sample signal for that regulator's knockout in this assay and under this scorer. It would also give A12's go/no-go (r >= 0.2 in Mesp1 labels; de_direction +0.03) from MERFISH alone, without cross-modal input.
- What it would not show:
  - that the motif proxy tracks beta-catenin activity;
  - that beta-catenin's response is mainly transcriptional within surviving labels;
  - that the magnitude rule transfers from a TF (Gata4) or a sparse non-TF (Mab21l2) to a co-activator;
  - that the gain survives the test file's own run effect.
- A Gata4 pass could also come from the GATA module coinciding with myocardial sub-state structure.

**8. Co-expression modules and removed compartments (question 2, part 2).**
- No systematic estimate was found (Gaps).
- From first principles, a lineage-restricted regulator's across-label covariance column is its compartment's identity module. Whenever the knockout shrinks that compartment, across-label fluctuation-response "predicts" the composition change.
- In the released data:
  - On Gata4, the counted composition shifts sit in cardiomyocyte states, with both signs (`gata4_calibration.csv`); by generic biology that is the deleted gene's main WT domain. Composition carried 0.11-0.16 of dt.
  - On Mab21l2, composition carried 0.07.
- Composition is closed on T3, so this matters only as a confounder of the within-label estimate.

**9. Verdict.**
- A knockout's expression response is predictable from WT covariation only in its systematic, low-dimensional part, and only with one fitted scale.
- Zero-shot sign and size need exogenous regulator variation, which a MERFISH snapshot does not provide.
- On T3, the part a WT-covariance predictor can reach is 10% of dt or less on the released knockouts, worth about 4.5 T3 points at oracle accuracy and about +0.2 expected at realistic accuracy.
- It is not a top-3 lever. At most it is an information try (A12's MERFISH-only readout), after a target-free audit.

## Candidate levers
Ranked by risk-adjusted expected value to a top-3 place. Both open levers are small: together at most about +0.3 Overall expected, against the 15-25 points needed.

**L1. T3 fluctuation-response proxy with a zero-shot validation ladder (open only through the captain).**
- Board and scored group: T3:bcatenin. Response identity (de_score, 0.30) and response direction (de_direction, 0.25) first; response magnitude (severity, 0.25) second.
- Mechanism: sections 4-6. Each panel gene is regressed, within label, on the deleted regulator's activity (folds) or the TCF/LEF motif proxy (test). The regression is depth- and spillover-adjusted, moment-corrected and shrunk, then extrapolated to complete loss in Mesp1-lineage labels and realised as per-label, per-gene count-space factors on WT cells.
- Expected effect:
  - Oracle ceiling about +4.5 T3 points (in-sample); realistic +0.05 to +0.7 T3, centre about +0.2 (synthetic).
  - P(passes F1-F3 under the adoption rule) about 0.10 (0.05-0.2). The step from regulator to proxy loses perhaps half again, so EV is about +0.02 to +0.05 Overall.
  - It is worth more on an entry at the severity floor (s2), where the payoff is convex, but s2 is frozen.
- Legality:
  - No hidden target is read. The inputs are WT E8.75 and E9.5 (released), Mab21l2 (released) and the Gata4 replicates (unlocked, A20), the knockouts for the folds only.
  - Class B: no gene-level knockout findings; motif binding data only; perturbation-derived regulons and footprints excluded.
  - Class G: no lists; directions come from code.
  - Class W: sequence-only resources; no in-window or validation-stage epigenomes.
  - No composition, scope, size or depth change, and no cross-modal input.
  - rules.yaml stays untouched; T3's s1 changes only by the captain's ruling after the folds pass.
  - Every fold run is a registered try; resources are logged in `EXTERNAL_DATA.md` first; the factors run through the existing core primitive.
  - Open question for the captain: the policy line requiring an `EXTERNAL_DATA.md` row for any pathway set also names A15 (section 6). If that gates motif sets, only Stage 0 and Stage 1 stay open.
- Validation route:
  - Stage 0 (target-free, WT files only, not a try). The split-half (spatial halves) Spearman of beta columns must be at least 0.5 in at least half the Mesp1 labels with 200 or more cells. The Spearman of beta with against without nuisance covariates must be at least 0.7. Also report the share of beta variance that label microclusters (state only) explain. A fail on either of the first two prunes the lever.
  - Stage 1 (tries): F1, F2 and F3 zero-shot under the adoption rule, with the two bases and the permuted control of section 5.
  - Stage 2 (tries): the GATA-motif proxy on F2 and F3.
  - Test use of the TCF/LEF proxy: the captain's blind bet only, and only after Stages 1-2 pass.
- Effort: 3.5-4 days (Stage 0: 0.5-1; Stage 1: 1.5-2 on the DGX; Stage 2: 1).
- Main risk: confounded covariance (state, position, admixture) gives wrong-sign module moves. The likely outcome is a "no" that costs 3-4 T3 tries.

**L2. T1 susceptibility-weighted shrinkage of the external direction (open; the captain's blind bet).**
- Board and scored group: T1:test. DE recovery (de_score, 0.25) and direction (de_direction, 0.25).
- Mechanism: the variability-sensitivity coupling [3][5][11] and the time-dependent linear response: over a partial step, slow high-variance modes carry most of the remaining change.
  - Weight: within each Feng group, the E10.5 excess dispersion v_{g,t}, i.e. analytic Pearson-residual variance minus 1, with a fixed overdispersion and a fixed mean floor.
  - Shrinkage: f' = f x v / (v + median v); genes below the floor keep Feng's factor.
- Expected effect: -0.5 to +1.5 T1 points, sign uncertain; EV about 0 to +0.2 Overall. The truth's DE genes on E9.5 -> E10.5 were mostly highly expressed genes falling (`experiments/harness/results/t1_de_v1/`), and those are often low-dispersion, so the shrinkage can cut correct moves.
- Legality:
  - Data are E10.5 (unlocked, A20) and Feng's E9.5 and E14.5 (X002, legal for T1).
  - The direction lever is open; there is no composition, depth or cross-modal change.
  - It changes Feng's direction, which is the captain's ruling, as Feng's own adoption was.
- Validation route:
  - (a) Target-free legal check on Feng's own data (not a try): within each group, the partial Spearman of E9.5 excess dispersion with |E9.5 -> E14.5 log ratio|, given mean expression, must be at least 0.2, else prune.
  - (b) Mechanism tries on the two counted T1 backtests, with the same shrinkage applied to B0's released-data momentum. The adoption rule applies, and the chamber-mix confound is named in advance.
  - (c) Test use on Feng: the captain's blind bet only, since no backtest can judge Feng's direction.
- Effort: 1.5-2.5 days.
- Main risk: it reduces to expression-level shrinkage, which the de_score null and the de_direction partial already neutralise, or it removes real large changes of canalised genes.

**L3. Natural-knockout reweighting, the exact tilt limit (CLOSED).**
- Board and scored group: T3; mutant state and the DE groups.
- Mechanism: resample each label's WT cells towards low regulator activity.
- Closed: it is a within-label composition change by importance weights, and T3's composition is closed (importance weights are no legal composition route on any board). Its expression content stays available through L1's factor form.
- No validation route, effort or risk applies while closed.

**L4. WT stage-axis (dominant-mode) term on T3 (CLOSED).**
- Mechanism: move WT cells along the WT E8.75 <-> E9.5 axis, the dominant mode of embryo-to-embryo variation.
- Closed:
  - Both folds already showed their lean in-sample (fac:stage; stage weights 0.22 / 0.32 towards E9.5 on Gata4), so neither can count.
  - s1's Gata4 factor already carries the E8.75-referenced offset (the swap review).
  - s2 is frozen.

**L5. Fluctuation-response programme from T1 scRNA off-panel readouts (CLOSED).**
- This was the digest's P3 action and A12's original form.
- Closed by A15: no new cross-modal input until the organisers answer.
- Even if opened, Gupta's power analysis [7] says that at scRNA capture the regulatory co-fluctuation is weak, while the state covariation is strong but confounded.

**L6. Lineage-aware deleted-gene term (no-op at test).**
- Ctnnb1 is off-panel, so the term changes nothing on the test board.
- It is used only as a shared control inside L1's folds, so that the deleted gene's own drop is not credited to the component.

## Deep Read Notes
- [1] CIPHER: dX = Sigma u from an exponential tilt; covariance from raw control-cell counts; u fitted per perturbation by least squares (their Eq. 61). Ten Perturb-seq datasets: real Sigma beats shuffled and mean-field (variance-only) nulls (KS 0.88 and 0.65); same-cell-type cross-study transfer R^2 0.73, cross-type poor. Responses ride about three eigengenes (participation ratio 3.07, which correlates negatively with R^2). Stated limits: linearity, timescale separation, perturbation size. Not addressed: no zero-shot magnitude test, no tissue heterogeneity.
- [3] Kaneko 2008 (re-deriving [4]): P(x;a) = N exp(-(x - X0)^2 / 2 alpha + C(a - a0)(x - X0)) gives d<x>/da = C alpha, with Gaussian form, constant C and constant alpha assumed. Evidence: the GFP-fusion selection in E. coli (evolution speed against fluorescence variance) and network simulations. Per-gene proportionality was then shown only in silico ("direct experimental support is not available yet"). The evolutionary form V_g proportional to V_ip has its own stability conditions.
- [7] Gupta 2022: TF burst to target covariation peaks hours later at about Spearman 0.2 (Hill coefficient 2), and stable states give flat, symmetric covariation. At 10% capture about 50,000 cells per time point are needed. K562 scNT-seq (13,679 cells): GATA1-correlated genes were 3.8-fold enriched among knockdown DE genes. Conclusion: feasible in principle, underpowered now.
- [12] Mitchel et al.: segmentation admixture (same plane and adjacent z-planes) confounds DE, neighbour influence, ligand-receptor and multicellular-programme analyses in MERFISH, CosMx and Xenium. Examples: 11 of 36 top up-regulated genes; about 25-30% of top programme hits. Correction: subcellular-neighbourhood NMF plus a CRF; factor labelling is manual.

## Gaps
- Searched for and not found:
  - a zero-shot test (no fitted magnitude) of a fluctuation-response or covariance predictor on knockout data in tissue: CIPHER fits u, and IDA ranks effects;
  - a systematic estimate of how often a deleted regulator's co-expression module marks the compartment its knockout removes;
  - a MERFISH-specific measurement of how detectable within-type regulator-target covariance is.
- Not read in full:
  - [4] (PMC verification page; PNAS returned 403), used through [3];
  - [8] (the fetch tool could not parse the PDF);
  - [14], [9] and [13], abstract or summary level only.
  - The fetch tool saved unparsed PDF copies in the harness's tool-results area; they were not opened (AGENTS.md).
- Counter-claims against the top lever (L1):
  - (i) CIPHER's strong fits suggest covariance prediction may work better than my simulations imply. However, its magnitude is fitted, its cells are homogeneous lines, and its signal sits in about three global modes: here those modes correspond to the shared factor and depth, which s1 already takes from Gata4.
  - (ii) [removed: A] The allowed broad phenotype, though, is a growth and allocation defect (composition, closed) plus protein-level adhesion roles, neither visible to mRNA covariance.
  - (iii) The decomposition may understate within-type biology, because the non-Mesp1 factor can absorb non-autonomous biological effects.
  - (iv) A Gata4 pass could come from myocardial sub-state structure and would not generalise.

## Compliance log
none

## END
