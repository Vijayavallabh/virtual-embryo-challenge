# Task 1 methods: forecasting the E12.5 single-cell distribution (Human Team track)

Research report for the NeurIPS 2026 Virtual Embryo Challenge, Task 1 (E8.5 / E9.5 / E10.5 to E12.5).
Compiled 2026-10-18. Test phase opens 2026-10-20; final submissions are due 2026-12-02.

Confidence labels: **[High]** = verified in a primary source, in code, or by my own computation in this session. **[Medium]** = verified indirectly, or a strong inference from verified facts. **[Low]** = plausible but unverified, or an opinion. When something could not be checked, the report says so.

---

## 0. Bottom line

1. **The literature has no tested recipe for this task.** All forecasting methods with published held-out-future tests (scNODE, PRESCIENT, MIOFlow, PI-SDE, scIMF; benchmarks scNODE-2024, scTimeBench-2026) were compared only with each other. None was compared with a copy-last or constant-velocity baseline. Most flow-matching, Schrödinger-bridge and OT methods were validated only on *interpolation* of held-out intermediate time points, in 5 to 100 dimensions. [High]
2. **In VEC itself, simple predictions beat dynamics models.** The organisers' autonomous neural ODE does not beat a constant shift. `pseudobulk_shift` scores 48.58, below the copy_last floor of 50. [High, live leaderboard and site]
3. **Change the problem from extrapolation to interpolation.** Legal external data exist on both sides of E12.5: anything at or before E9.5, and anything after E13.5. Two good sources:
   - **Feng et al. 2022 (GSE193346):** 10x 3' v3, chamber-dissected hearts, E9.5 to P9. The E9.5 and E14.5+ samples are legal.
   - **The JAX time-lapse atlas:** 11.4M nuclei, E8 to P0. Stages E13.75+ and ≤E9.5 are legal.

   Waddington-OT's own validation found that OT interpolation of a held-out middle stage is "almost as accurate as the batch-to-batch baseline" for 1-day gaps, and degrades for 2-day gaps. [High] Competitors already do this: leaderboard model names include "Feng2022GSE193346_E95E145", "[Feng2022 WT E9.5+E14.5; MOSTA E14.5]" and "JAX hep+mes". [High that the names exist; Low on what those models actually do]
4. **What the metrics reward**, from reading the scorer and running my own probe [High]:
   - de_score and de_direction depend only on the prediction's **pseudobulk** (mean log1p per gene).
   - mmd_u depends on **composition and cell-state positions** in the top 30 PCs fitted on the truth.
   - The variogram on 32K genes is in practice a test of **per-gene marginals and sparsity**. In my probe, 72% of random gene pairs reduce to a single-gene term. Shuffling genes independently across cells left the variogram at the ceiling value (0.00263 vs 0.00252), while it more than doubled mmd_u (0.212 vs 0.091).
5. **Never add a dense delta to sparse cells.** Even an *oracle* mean shift added densely got the ceiling de_score but a variogram 3.2x *worse than copy_last* (0.0081 vs 0.0025). This is why `pseudobulk_shift` falls below the floor (variogram skill 38.9). Count-space edits (binomial thinning or Poisson top-up) and real cells keep sparsity. [High]
6. **Counts can be recovered exactly from the released `.X`.** The formula is `counts = expm1(X) * N / 1e4` with `N = 1e4 / expm1(min nonzero per cell)`. On the veckit sample files the result is integral (mean absolute deviation from integers 0.0000), with median depth of about 11K UMIs per cell. Count models (scVI NB, CFGen) can therefore be trained and sampled at the real depth. [High on the sample files; check on the full release]
7. **Recommended plan (Section 4):**
   - Quick win: a count-space pseudobulk shift (days).
   - Main bet: composition plus per-type count-space shifts, with deltas taken from the external bracketing data (1-2 weeks).
   - Second, diverse bet: moscot OT bracketing from E10.5 to Feng E14.5 in a shared scVI latent, with count-space decoding (2-3 weeks).
   - Whole-ladder non-autonomous flow matching with growth (VGFM-style) only if time allows.
   - With only 2 scored test submissions and best-of-two scoring, submit two diverse bets, not two near-copies.
8. **Foundation models (FMs) add little and several are illegal here** [High/Medium per model]:
   - CELLxGENE has carried in-window mouse embryo data since 2021-2022: LungMAP lung at TS21 (about E12.5-13) and developing thalamus at E12.5. The JAX atlas (11,441,407 cells, Theiler stages 16-22 included) was added on 2024-03-11.
   - Disallowed: UCE (its training list includes six in-window mouse datasets), TranscriptFormer TF-Exemplar/TF-Metazoa, scPRINT, Nicheformer, GeneCompass, Mouse-Geneformer, C2S-Scale, and the Navigo/Chreode checkpoints.
   - Allowed as far as the mouse window goes: human-only models (scGPT whole-human, Geneformer, scFoundation, State SE, TF-Sapiens, and Tabula, which is likely human-only). They are of little use for 32,285 mouse genes, of which 8,522 are Gm* and 2,524 are *Rik.

---

## 1. What the Task 1 metrics actually reward

### 1.1 Scorer facts relevant to method choice
These come from `veckit/common/core_metrics.py` and `T1/metrics_v2.py` (MIT, github.com/aristoteleo/veckit). [High]

| Metric (weight) | Depends on | Implication for the generator |
|---|---|---|
| de_score (0.25) | Only `dp = pb(pred) - pb(ref)`, ranked. Truth DE genes are those with \|Δmean log1p\| ≥ 0.25 against E10.5; rank-based, sign-split, minus an expression-level null. | Get the *ranking* of per-gene mean changes right. Scale is irrelevant; noise in the top genes hurts. |
| de_direction (0.25) | Rank partial correlation of dp vs dt over **all 32,285 genes**, controlling for rank pb(ref). | Every gene's rank counts equally, including thousands of near-zero genes. Avoid injecting noise into genes that should not change. Keep exact zeros. Submit many cells so pseudobulk noise stays low. |
| mmd_u (0.30) | Multi-scale RBF MMD² on ≤2,000 cells in a 30-PC space fitted **on the truth only**. | Composition and type-specific profiles of the high-variance genes (markers: globins, CM, hepatocyte, collagen and endothelial genes). Mode collapse and over-smoothing are both penalised. |
| variogram (0.20) | 20,000 random gene pairs on ≤1,500 cells; compares E\|x_i - x_j\|^0.5. | With about 88% zeros and 47% of genes all-zero in small samples, most pairs reduce to a single-gene moment E\|x\|^0.5. Per-gene detection rate and nonzero level dominate. |

- The skill map is hyperbolic: floor maps to 0.5, ceiling to 1.0. Per-metric skill is clipped at 1.
- The floor is copy_last. For the test board this is E10.5 verbatim, which is also the DE reference.
- T1 validation anchors: de_score ceiling 0.846; de_direction ceiling 0.790; mmd_u floor 0.0836 / ceiling 0.0041; variogram floor 0.00522 / ceiling 0.000158 (`index.json`).

### 1.2 Live leaderboard calibration
Source: `https://api.virtualembryo.ai/challenge/leaderboard?task=T1`, updated 2026-10-18 15:41 UTC, phase P2. [High]

- 200 rows (151 human, 49 agent). Every row is above 50.
- Best headline is 76.00 (human): de_score 63.8, de_direction 74.1, mmd_u 84.7, variogram 80.5.
- Per-metric maxima across rows: de_score 63.8, de_direction 78.5, mmd_u 84.7, variogram 84.5. Medians: 50.1 / 60.8 / 62.4 / 60.8.
- Reference row `pseudobulk_shift`: 48.58 (de_score 46.2, de_direction 53.9, mmd_u 52.5, **variogram 38.9**).
- Model names hint at approaches (names only, methods unpublished) [Low]:
  - "PopTransport-v4 … Endo0.24-Epi0.095-Thin0.85"
  - "t1_comp_thin05_trend0.5"
  - "e95_typeshift_sparse_f1.5"
  - "t1mean_Feng2022GSE193346_E95E145…"
  - "[Feng2022 WT E9.5+E14.5; MOSTA E14.5]"
  - "r228 m6 JAX hep+mes s5"
  - "r1_20260909_t1_bures_subspace"
  - One agent row is self-labelled "NONCOMPLIANT E10.5 prior".
- Validation is E9.5 to E10.5 (1 day). The test is E10.5 to E12.5 (2 days), so changes and errors will both be larger.

### 1.3 My metric probe (veckit's own code, the 150-cell sample files)
- Script: `scratchpad/scripts/t1_probe.py`, run on the DGX in `~/v/vec-t1-scratch` (uv env `envs/t1`).
- Setup: truth = random 75-cell half of `sample_9.5`; ref = `sample_8.5`; mean over 3 split seeds.
- Small samples, so the numbers are directional only. [High as a computation; Medium on how far it transfers to full size]

| Prediction | de_score | de_direction | mmd_u | variogram | frac. zeros |
|---|---|---|---|---|---|
| Ceiling analog: other 75-cell half of E9.5 | 0.585 | 0.377 | 0.0910 | 0.00252 | 0.873 |
| Floor analog: copy E8.5 | 0 | 0 | 0.1109 | 0.00251 | 0.878 |
| **Oracle dense shift**: clip(E8.5 + Δpb_true, 0) | 0.585 | 0.371 | 0.1142 | **0.00811** | 0.657 |
| Oracle **count-space** shift of E8.5 (thin/top-up; crude target) | 0.311 | 0.325 | 0.1276 | 0.00334 | 0.899 |
| kNN-mean (k=10) of E9.5 half (dense "decoder mean") | 0.515 | 0.311 | 0.1474 | **0.01890** | 0.656 |
| Same, then Poisson-resampled at each cell's own depth | 0.496 | 0.258 | 0.1391 | 0.00459 | 0.882 |
| Cell-type means (dense) | 0.585 | 0.377 | 0.1184 | 0.01062 | 0.733 |
| Type means, then Poisson-resampled | 0.525 | 0.303 | 0.1099 | 0.00356 | 0.888 |
| E9.5 half with each gene shuffled independently (marginals kept) | 0.585 | 0.377 | **0.2120** | **0.00263** | 0.873 |

Other facts from the same files [High, sample files only]:
- `expm1(X)` row sums are exactly 1e4. Recovered counts are integral.
- Median depth is 10,771 (E8.5) and 11,466 (E9.5) UMIs; 10th-90th percentile about 6.2K-22.4K.
- Zeros are 87.5-87.8% of entries; 47% of genes are zero in all 150 cells.
- Median genes detected is 3,862-4,016 per cell.
- Median mitochondrial fraction is 5.9-6.4% and Malat1 is 0.33%. That profile fits **whole dissociated cells (10x 3')**, not nuclei. The JAX atlas is nuclei, so this is a platform shift. [Medium]
- In seed 0, 22.6% of variogram pairs had both genes all-zero in both sets, and 72.1% had at least one such gene. 14,183 genes had dt exactly 0.
- The 32,285-gene panel order (Xkr4, Gm1992, Gm19938, Gm37381, Rp1, Sox17, …, AC149090.1) appears to be Cell Ranger's mm10-2020-A feature list (Ensembl 98). [Medium, from memory, not re-verified]

What this means for design:
- The three metric families can be engineered somewhat separately: pseudobulk ranks for DE, composition for MMD, per-gene sparsity for the variogram.
- Representation matters more than the dynamics model. Dense or decoded means are heavily penalised.

---

## 2. Method survey (Q1)

Column meanings:
- "Extrap.": was prediction past the last observed time point evaluated?
- "New states": can the method produce cell states absent from the source?
- "32K genes": how it scales to the full panel.

| Method (venue) | Mechanism | Extrap. | New states | Growth | Code / licence / last push | 32K genes | Published held-out evidence | T1 verdict |
|---|---|---|---|---|---|---|---|---|
| **Waddington-OT** (Schiebinger et al., Cell 2019) | Entropic unbalanced static OT between consecutive snapshots | **No** | Only if present at a later observed stage | Yes (proliferation/apoptosis prior, refined) | `broadinstitute/wot`, BSD-3, push 2022-09; PyPI 1.0.8.post2 (2020) | Local 30-PC | Held-out midpoints ≈ batch-to-batch noise at 1-day spacing; worse at 2 days [High] | Use only for **bracketing interpolation** |
| **moscot** TemporalProblem (Klein et al., Nature 2025) | Same idea, OTT-JAX/GPU, low-rank OT; push/pull, cell_transition, interpolation validators | **No** | Same as WOT | Yes (`score_genes_for_marginals`, tau_a/tau_b) | `theislab/moscot`, BSD-3, v0.5.2 (2026-08-25), active; Python ≥3.11 | Any embedding (PCA/scVI); 1.7M cells × 20 time points shown | WOT-style `compute_interpolated_distance` [High] | **Best tool for bracketing.** Its interpolation is a dense linear blend `(1-t)x_s + t x_t`, so replace it with count-space sampling (Section 5) |
| **TrajectoryNet** (Tong et al., ICML 2020) | CNF with dynamic-OT regulariser | No (leave-one-out interpolation) | Weak | Optional (WOT rates) | `KrishnaswamyLab/TrajectoryNet`, **Yale non-commercial; derivatives must be posted on GitHub** | ≤ ~5-10 dims | Interpolation only [High] | Skip (licence risk, low dims) |
| **MIOFlow** (Huguet et al., NeurIPS 2022) | Neural ODE in geodesic-AE latent, manifold OT | Interpolation in paper; extrapolation tested by others | Weak | No | `KrishnaswamyLab/MIOFlow`, **Yale non-commercial**, push 2026-10 | Latent | Lost to scNODE on SC/DR extrapolation; degrades on whole transcriptome (scTimeBench) [High] | Skip |
| **PRESCIENT** (Yeo et al., Nat Commun 2021) | Potential-gradient SDE in 30-50 PCs | Held-out day 4 only (interpolation) | No (assumes final-stage steady state) | Yes (signature or clonal) | `gifford-lab/prescient`, MIT, push 2025-01 | PCA | Beat endpoints and WOT linear interpolation on Weinreb day 4; lost to scNODE on extrapolation [High] | Skip, except as a growth idea |
| **scNODE** (Zhang et al., Bioinformatics 2024, ISMB) | VAE + latent neural ODE + dynamic regulariser; Wasserstein loss | **Yes** (last-k held out) | Only via latent drift | No | `rsinghlab/scNODE`, MIT, push 2024-10 | 2,000 HVGs (VAE) | Better than PRESCIENT/MIOFlow, especially under large shift (SC t=15: W 132.9 vs 162.1/150.5); error grows with horizon; **no naive baseline** [High] | Candidate latent-dynamics module only |
| **PI-SDE** (Jiang & Wan, Bioinformatics 2024, ISMB) | HJ-regularised neural SDE (PRESCIENT extension) | Held-out time points | No | No | `QiJiang-QJ/PI-SDE`, **no licence**, push 2025-05 | PCA | Mid-pack in scTimeBench; OOM on whole transcriptome [High] | Skip |
| **scIMF** (Jiang et al., arXiv 2505.16492) | Transformer-enhanced mean-field neural SDE (cell-cell interaction) | Yes (scNODE-style tasks) | No | No | `QiJiang-QJ/scIMF`, **no licence** | Latent | Best overall forecaster in scTimeBench; one reported medium-extrapolation task W1 21.14 vs scNODE 26.51 [High; dataset not checked] | Reimplement only if a latent-dynamics component is needed |
| **TIGON** (Sha et al., Nat Mach Intell 2024) | WFR dynamic unbalanced OT (velocity + growth) | No | No | **Yes** | `yutongo/TIGON`, MIT | Low-dim (AE/PCA) | Reconstruction at unmeasured times [High] | Concept only |
| **DeepRUOT** (Zhang, Li, Zhou, ICLR 2025 oral) | Regularised unbalanced OT ≈ unbalanced Schrödinger bridge; velocity + growth + score nets; simulation-based | No | Weak | **Yes**, learned | `zhenyiizhang/DeepRUOT`, MIT stated in setup.py only (no LICENSE file), push 2025-12 | 2-10 D shown; ~90 min on EB 50-D | Best on low-dim synthetic and hematopoiesis among SB/OT baselines [High] | Too slow and low-dim; concept only |
| **VGFM** (Wang et al., NeurIPS 2025; cited in VEC proposal [25]) | Joint velocity + growth flow matching from semi-relaxed OT; plus W1 fitting loss | **No** (hold-one-out interpolation) | Weak | **Yes** | `DongyiWang-66/VGFM`, **no licence**, push 2026-01 | ≤1000-D synthetic; EB/CITE 50-D; Pancreas 2,000 HVGs | Best or near-best on hold-one-out (EB-5D W1 0.676; CITE-50D lost to OT-MFM); 13 min vs DeepRUOT 90 min [High] | Best "growth-aware" dynamics; **reimplement** on torchcfm/POT |
| **Schrödinger bridges**: DMSB (NeurIPS 2023), Light SB (ICLR 2024), SF2M (AISTATS 2024), UDSB (2023), GSBM (ICLR 2024), CytoBridge (NeurIPS 2025 toolkit) | Entropic dynamic OT; multi-marginal / momentum / unbalanced variants | No (intermediate-day recovery) | Weak | UDSB / CytoBridge yes | DMSB MIT (2023); LightSB MIT; UDSB MIT; **GSBM CC BY-NC, archived**; CytoBridge **GPL-3.0** | 50-1000 D (LightSB 1000-D MSCI) | Light SB on MSCI days 2/3/4/7: ED 2.31 / 2.16 / 1.27 at 50 / 100 / 1000 D [Medium, from a search snippet, PDF not opened] | Only as interpolators between brackets |
| **Flow matching**: torchcfm OT-CFM (TMLR 2024), MMFM spline multi-marginal (ICLR 2025), MMSFM (ICML 2025), Metric FM (NeurIPS 2024), Curly-FM (NeurIPS 2025) | Simulation-free regression of velocity along (OT-)paired interpolants; splines across ≥3 marginals | No; MMFM imputes missing times | Weak | No (VGFM adds it) | torchcfm MIT, v1.0.7, 2.6K stars; MMFM Apache-2.0; MMSFM **no licence**; Metric FM MIT; Curly-FM MIT | Latent / PCA | Leave-one-time-point-out on EB/CITE/MULTI; MMFM beats baselines at imputing missing times [High] | **MMFM-style multi-marginal splines through E8.5…E10.5 + legal ≥E13.75 marginals** make E12.5 an *imputation* problem |
| **Navigo** (Fan et al., bioRxiv 2026; likely the `dynode_flow` reference) | Flow matching + reflow + RNA kinetics on (spliced, unspliced); metacells | Interpolation demos | No | No | `aristoteleo/Navigo-release`, BSD-2, push 2026-06 | 4,000 HVG atlas | Demo: held-out intermediate EMD 5.37 vs 5.97 start-anchor [High] | Code OK; **checkpoints trained on the JAX atlas incl. E9.75-E13.5 are illegal**; needs unspliced counts we do not have |
| **Action matching** (Neklyudov et al., ICML 2023) | Gradient-field dynamics from marginals; unbalanced variant | No | No | uAM | `necludov/action-matching`, **no licence** | Low-dim | uAM loses to DeepRUOT in its tables [High] | Skip |
| **scDiffusion** (Luo et al., Bioinformatics 2024) | Latent diffusion + foundation-model autoencoder; classifier guidance; "Gradient Interpolation" | No | Out-of-training combinations of conditions | No | `EperLuo/scDiffusion`, MIT | Via AE (which FM was not checked) | Generated a continuous mouse-embryo trajectory (qualitative) [High] | Skip (decoder means; FM dependency) |
| **Squidiff** (He et al., Nat Methods 2026) | Diffusion with semantic latent; linear arithmetic between condition latents | Not quantitatively (interpolated days 1-2 between 0 and 3) | Via latent arithmetic | No | `siyuh/Squidiff`, BSD-3 text (GitHub shows NOASSERTION because of merge-conflict markers) | Not stated | scTimeBench: best Gaussian MMD, worst on other metrics [High] | Skip |
| **CFGen** (Palma et al., ICLR 2025) | Latent flow matching + discrete count likelihoods; compositional guidance | No | Attribute combinations | No | `theislab/CFGen`, MIT (PyPI "cfgen" is an unrelated 2017 package) | Whole genome | Realistic counts, rare-type augmentation [High] | Optional **count decoder** for emergent states |
| **scVI** (Lopez et al., Nat Methods 2018) | NB/ZINB VAE with batch covariates | No (covariate extrapolation unreliable) | No | No | `scvi-tools` 1.5.1, BSD-3, Python ≥3.12 | All genes | Standard; posterior predictive gives counts [High] | **Shared latent + count decoder** for bracketing |
| **CellFlow** (Klein et al., bioRxiv 2025) | Conditional flow matching with condition embeddings | Unseen *conditions* | Via conditions | No | `theislab/CellFlow`, MIT, active | Latent | Generalises to unseen perturbation conditions [Medium] | Possible time-as-condition generator; unvalidated for time |
| **GeneSys** (bioRxiv 2025.08.20.671385) | Bi-LSTM + TD-VAE over a user-defined **lineage blueprint**; zero vectors for late types | No clean test | **Yes, by design (user specifies lineage)** | No | Repo not checked | 24,552 genes, ~100K cells | Recreation accuracy needs ~6-7 input bins (mouse) [High] | Idea worth borrowing: an explicit lineage map for emergent types |
| **dynamo** (Qiu et al., Cell 2022) / **CellRank 2** | Vector fields from RNA velocity or labelling; CellRank RealTimeKernel uses moscot couplings | Vector-field extrapolation possible | No | No | dynamo-release 1.5.3, BSD-3; cellrank 2.3.3, BSD-3 | Velocity genes | Fate mapping, not distribution forecasting [High] | Velocity unavailable (log1p `.X` only). CellRank RealTimeKernel can estimate which E10.5 states expand |
| **Latent ODEs** (Rubanova et al., NeurIPS 2019) | ODE-RNN encoders | Generic | No | No | various | Latent | Not single-cell specific | Same caveats as scNODE |
| **Chreode** (arXiv 2605.28111) | Cell "world model" (scVI encoder + DiT) pretrained on a 2.4M-cell mouse embryo atlas (7 datasets) | Not shown | Not shown | No | No repo or weights found | scVI | Fine-tuning gains vs PI-SDE/PRESCIENT (no numbers in abstract) | Weights presumed contaminated; skip |

Overall:
- No surveyed method has evidence of *creating* new cell types when forecasting past the last observed stage. Models trained only on E8.5-E10.5 can only drift existing states. [High]
- Emergent E12.5 states (e.g. epicardium-derived fibroblasts, smooth muscle and pericytes, coronary endothelium, definitive erythrocytes, immune cells) must come from **external data after E13.5** (Feng 2022 annotates Fb-like, SMC/pericyte, vas_EC, immune and RBC clusters) or from an explicit lineage prior. [High for Feng's annotations; Medium for which states appear by E12.5 in this dissection]

---

## 3. Benchmarks of temporal forecasting (Q2)

**scNODE (2024)** [High]
- Three datasets: ZB (zebrafish, 38,731 cells, 12 time points), DR (Drosophila, 27,386 cells, 11 time points), SC (Schiebinger reprogramming, 236,285 cells, 19 time points). All on 2,000 HVGs.
- Tasks: easy (interpolation), medium (last time points held out) and hard (both).
- scNODE beat PRESCIENT and MIOFlow in most cases. Its advantage correlated with distribution shift (Spearman up to 0.93), and accuracy fell with horizon.
- **Only deep baselines; no copy-last or constant-velocity baseline.**

**scIMF (2025)** [Medium]
- Same task design as scNODE. Reports W1 21.14 vs scNODE 26.51 on a medium (extrapolation) task. I did not check which dataset.

**scTimeBench** (Osakwe et al., Bioinformatics 2026, ISMB; MIT) [High]
- 9 methods (7 forecasting + 2 OT; the abstract says 10) on 8 datasets across 4 species. Metrics: Wasserstein, Gaussian MMD, energy MMD, Hausdorff (geomloss).
- scIMF best overall; scNODE, MIOFlow and PI-SDE similar; Squidiff best only on Gaussian MMD. PI-SDE ran out of memory on whole-transcriptome data. MIOFlow degrades on whole-transcriptome and cross-batch data.
- **No forecasting naive baseline.** A Spearman-correlation baseline was competitive for lineage fidelity. Adding pseudotime helped lineage recovery, but dataset-dependently.

**Interpolation evidence** [High]
- WOT: OT interpolation ≈ batch-to-batch noise at 1-day spacing, degrading at 2 days.
- PRESCIENT: held-out day 4 better than the endpoints and WOT linear interpolation.
- VGFM: hold-one-out wins on EB/CITE at 5-D.
- Navigo demo: prediction EMD 5.37 vs start-anchor 5.97 on held-out intermediate stages.

**VEC organisers** [High]
- The neural ODE (autonomous, flow-matched on global entropic OT) does not beat the constant shift.
- Constant-velocity extrapolation fitted off the evaluation stages lands below copy_last on every setting.
- Their stated need: "supervision across the whole stage ladder, a non-autonomous field, and a treatment of populations appearing".

**Adjacent field: perturbation prediction** [High]
- Five foundation models and two deep models did not beat deliberately simple linear baselines (Ahlmann-Eltze et al., Nat Methods 2025).
- Zero-shot scGPT and Geneformer were often beaten by simpler methods (Kedzierska et al., Genome Biol 2025).

**Two recent items I found but did not review in detail**
- A 2026 zero-shot FM-embedding benchmark for dynamics (bioRxiv 10.64898/2026.03.10.710748).
- "From Snapshots to Trajectories … Conditional Flow Matching" (arXiv 2605.22340).

**Do simple baselines win? What generalises?**
- The forecasting literature cannot tell us whether deep models beat simple baselines, because none tested one. In VEC, simple baselines win among the published references.
- What generalises [Medium]:
  1. **Interpolation over extrapolation.** Bracket the target whenever legal data exist.
  2. **Real cells and count-faithful outputs**, given the variogram and MMD findings above.
  3. **Shrinkage toward the last observed state.** copy_last is strong; momentum should be damped.
  4. **Explicit composition and growth modelling, as noisy priors.** PRESCIENT's gene-signature growth rates correlated only r = 0.207 with clonal growth, yet still helped fate prediction.

---

## 4. Candidate pipelines for E12.5 (Q3)

### 4.1 Ground rules

**Legal external windows for T1 extrapolation** [High, rules text]
- Allowed: data at or before E9.5, and anything strictly after E13.5. E13.5 itself is excluded. Disclose every source and its stages.
- Resources that span the window: remove in-window cells first, then describe the removal. Do this before normalisation, HVG selection or any model fitting.

**Bracketing sources** (stage coverage verified as stated)
- **Feng et al. 2022** (Nat Commun 13:7960; GSE193346) [High]
  - 10x 3' v3, the same chemistry family as T1.
  - CD1 (18 stages) and C57BL/6 (17 stages), E9.5 to P9.
  - LA/RA/LV/RV chamber dissections; 29,001 + 25,605 singlets; about 400 cells per sample.
  - Use **E9.5 (also usable to calibrate platform and dissection against our own E9.5)** and **E14.5, E15.5, …**. Drop E10.5-E13.5.
- **JAX time-lapse** (Qiu et al., Nature 2024) [High]
  - 11,441,407 nuclei from 74 embryos; 24,552 genes; median 2,545 UMIs per nucleus (shallow); whole embryo.
  - Somite-staged to E10, then one embryo per 6-hour bin.
  - Use ≤E9.5 and ≥E13.75 (bins E13.75, E14.0, E14.25 are listed). Drop E9.75-E13.5.
  - On CELLxGENE as "A single-cell transcriptional timelapse of mouse embryonic development".
- **MOSTA** Stereo-seq E14.5+ is used by at least one competitor. [Medium]
- The dataset agent's report should be the authority on further sources.

**Backtest protocol (available from 2026-10-20, when E10.5 answers are released)**
- Rebuild each pipeline to predict E10.5 from E8.5/E9.5 plus legal data, keeping the same t-logic. Bracket from our E9.5 to Feng E14.5 is t = 0.2; to JAX E13.75 it is t = 0.235.
- Score with veckit over many 10%-then-half resamples of E10.5. Use this to pick hyperparameters.
- The val leaderboard (8 submissions/day) stays open until then. Using scores to *choose* is fine; solving for target properties from scores is prohibited.

**Two test submissions, best-of-two counts**
- Pick two predictions that are each good, and good in *different* scenarios (composition, time-warp, aggressiveness). That maximises the expected maximum.
- Useful identity [High, from the formulas]: squared MMD and the variogram are convex in mixture weights. A population mixture of two predictions is never worse than the weighted average of their values on these two metrics. Mixtures are therefore a safe hedge, but they do not raise the maximum of two diverse bets.

### 4.2 The pipelines (ranked by expected score per unit effort)

#### P1. Count-space pseudobulk shift (quick floor-beater). Effort 1-2 days. Rank 1 on value per effort.
1. Recover counts.
2. Harmonise labels across E8.5/E9.5/E10.5 (scANVI or marker mapping; e.g. Endocardium → Endocardium-1/-2).
3. Per type k: Δ_k = λ·(mean_E10.5 − mean_E9.5) in log1p, with λ in about [0.3, 1] tuned on the backtest.
4. Apply Δ_k in **count space**: binomial thinning for decreases; Poisson/NB top-up proportional to existing counts for increases; renormalise to 1e4 and take log1p.
5. Output about 5K real E10.5 cells (test max_cells will be at least 5,000) at the E10.5 composition, with mild composition momentum.

Metric reasoning:
- variogram: stays near the floor or better, unlike dense pseudobulk_shift (38.9).
- de_direction: shares the sub-floor-to-floor gains of pseudobulk_shift (53.9 at val).
- de_score: limited, because momentum misses dissection-driven and emergent genes.
- mmd_u: about the floor.
- Expected: modestly above 50. [Medium] Also the base layer for P2.

#### P2. Composition + per-type count-space shifts with **external bracketing deltas** (main bet). Effort 1-2 weeks. Rank 2.

1. **Per-type gene trends across the gap.** For each harmonised type with an external counterpart (CM subtypes, endocardium, vascular EC, epicardium/mesothelium, fibroblast-like, SMC/pericyte, RBC, macrophage; hepatocyte and foregut from JAX):
   - Fit log-mean expression against time on the *legal* stages separately in Feng and JAX: linear or monotone spline in time, or in log time.
   - Take Δ_k^ext = f_k(12.5) − f_k(10.5). This is an **interpolation**: both ends are bracketed by data.
   - Combine with internal momentum: Δ_k = w_F·Δ_k^Feng + w_J·Δ_k^JAX + w_I·Δ_k^internal, with empirical-Bayes shrinkage per gene (precision from cell counts).
   - Tune the weights on the E10.5 backtest, where the analogous target is f(10.5) − f(9.5).
   - Use **within-dataset fold changes**, never absolute external levels. Cells vs nuclei and depth differ: 11K vs 2.5K UMIs.
2. **Composition π′(E12.5)** is the largest uncertainty, because the dissection changes by stage.
   - Start from E10.5 proportions (known from 2026-10-20).
   - Apply damped log-odds trends from E8.5 → E9.5 → E10.5 on harmonised types.
   - Add the biology prior that cardiac compartments expand and new lineages appear. Feng's compositions at E14.5 cover cardiac types only.
   - Keep every E10.5 type at nonzero weight.
3. **Emergent states** (absent at E10.5, present in Feng E14.5): seed a small mass, about 2-8% in total, set on the backtest. Build them from (a) the nearest E10.5 progenitor cells (e.g. epicardium for EPDC fibroblasts and SMCs; endocardium or endothelium for vascular EC; primitive erythroid for definitive), pushed in count space toward the external type's fold-change profile, or (b) the scVI posterior predictive with batch = ours (see P3).
4. **Generate.** Resample E10.5 cells per type to π′; apply count-space edits at each cell's own depth; renormalise; take log1p.
5. **Calibrate.** Optionally adjust the *ranks* of the prediction's pseudobulk to the composed target by rescaling nonzero entries per gene. This never turns zeros into small values.

Metric reasoning:
- de_score / de_direction gain directly from Δ^ext. External time trends carry the E10.5 → E12.5 biology that momentum cannot, e.g. the primitive-to-definitive globin switch and CM maturation.
- mmd_u gains if π′ is right.
- variogram is protected by count-space edits.

Main risks:
- Composition error.
- Gene symbol and annotation mismatches (map by Ensembl ID).
- Feng lacks non-cardiac tissue present in the T1 dissection.

#### P3. **OT bracketing** E10.5 → post-E13.5 in a shared latent, with count-space decoding (second, diverse bet). Effort 2-3 weeks. Rank 3.
1. scVI (NB) on our E8.5/E9.5/E10.5 plus Feng legal stages (plus optional JAX ≤E9.5 and ≥E13.75 heart-relevant subclusters), with dataset or technology as the batch covariate. Check that Feng E9.5 and our E9.5 co-embed.
2. moscot `TemporalProblem` from our E10.5 to Feng E14.5 (or a pooled E14.5-E15.5 anchor):
   - Proliferation/apoptosis marginals from mouse gene sets; tau_a, tau_b < 1 (unbalanced); low-rank if needed.
   - Validate with `compute_interpolated_distance` on the backtest (our E9.5 → Feng E14.5, evaluated at E10.5).
3. Sample pairs (i, j) from the coupling, accounting for growth, at t = 0.5. Calibrate the effective t on the backtest.
4. Produce cells **without dense blending**:
   - Either draw counts per gene from the source cell's counts with probability (1 − t) and from the target cell's *our-platform* counterfactual with probability t, using binomial mixing.
   - Or decode the latent interpolant with scVI posterior predictive sampling (batch = ours, library size from the E10.5 depth distribution).
5. Non-cardiac compartments (foregut, hepatocyte/ST, pharyngeal mesoderm, pericardium) are missing from Feng. Carry them from P2's composition model.

Pros:
- Emergent states and within-lineage changes come straight from real post-E13.5 cells.
- It is interpolation, which has the best evidence (WOT).

Risks:
- The 4-day bracket is wider than WOT's validated 1-2 days.
- Batch correction across dissections and platforms.
- Composition must still be controlled explicitly.

#### P4. Whole-ladder non-autonomous flow matching with growth. Effort 3+ weeks. Rank 4.
- Train v(z, t) and g(z, t) (VGFM-style: semi-relaxed OT couplings, growth regression, plus a W1 fitting loss) or MMFM-style multi-marginal splines in the shared scVI latent.
- Marginals: ours E8.5/E9.5/E10.5 plus Feng E9.5, E14.5, E15.5 plus JAX legal stages, after mapping into our batch.
- E12.5 then sits *between* supervised marginals.
- Integrate the E10.5 cells to t = 12.5 with growth weights, resample, and decode with NB sampling.
- This matches the organisers' hint (whole ladder, non-autonomous, appearing populations) but has the least evidence. VGFM and MMSFM code is unlicensed, so reimplement on torchcfm (MIT) + POT (MIT) or ott-jax (Apache-2.0).
- Use as a module (growth and composition estimates) unless P2/P3 plateau.

#### P5. Time-conditioned decoder generation (scVI / CFGen / CellFlow with continuous time + dataset covariates). Effort 2 weeks. Rank 5.
- Sample at t = 12.5 in our batch.
- Risks:
  - Covariate interpolation inside a decoder blends programmes.
  - Generated cells tend to be over-smoothed: see the probe's decoded-mean rows (variogram 0.0189 dense; 0.0046 even after resampling vs floor 0.0025).
- Best used only for emergent states inside P2/P3.

| Rank | Pipeline | Effort | de_score | de_direction | mmd_u | variogram | Main risk |
|---|---|---|---|---|---|---|---|
| 1 | P1 count-space shift | days | small + | + | ≈ floor | ≈ floor / + | Momentum misses emergent biology |
| 2 | P2 composition + external-delta shifts | 1-2 wk | ++ | ++ | + / ++ (composition) | + | Dissection composition |
| 3 | P3 OT bracketing + count decoding | 2-3 wk | ++ | ++ | ++ if integration is good | + / ≈ | Batch effects, wide gap |
| 4 | P4 whole-ladder FM + growth | 3+ wk | + / ++ | + / ++ | + / ++ | depends on decoder | Effort, evidence |
| 5 | P5 conditional decoder | 2 wk | + | + | ± | − / ≈ | Over-smoothing |

The +/− entries are qualitative expectations. [Low/Medium] Do not treat them as score predictions.

### 4.3 Practical defaults
- Bracketing time parameter: E10.5 → E14.5 gives t = 0.5; → E13.75 gives t = 0.615. Development is not linear in days, so fit a monotone warp on the backtest (t_true = 0.2 or 0.235) and carry it over.
- Submit the maximum allowed cells (≥5,000). Pseudobulk uses every submitted cell; MMD and the variogram subsample 2,000 and 1,500.
- Check with veckit before submitting:
  - `library_size_ratio` ≈ 1 (expm1 row sums exactly 1e4).
  - `variance_ratio` not ≪ 1.
  - No negatives or NaNs.
  - Exact gene order from `T1__val.genes.txt`; the test panel file follows at P3.

### 4.4 Licence, compliance and compute notes

**Code licences** [High on licence texts; Low on how organisers would treat prize money]
- The rules allow published code only "if licences permit".
- Unlicensed repositories grant no permission: VGFM, PI-SDE, scIMF, MMSFM, action-matching. Reimplement their *methods* from the papers.
- Yale non-commercial licences (TrajectoryNet, MIOFlow) require derivatives to be published on GitHub. That clashes with the no-code-sharing rule, and prize money may count as "monetary compensation".
- GSBM is CC BY-NC; Tabula is non-commercial; CytoBridge is GPL-3.0 (copyleft).
- Safe building blocks: moscot, scvi-tools, cellrank, dynamo (BSD-3); torchcfm, POT, geomloss, LightSB, DMSB, TIGON, scNODE, PRESCIENT, CFGen, CellFlow (MIT); ott-jax, MMFM (Apache-2.0); Navigo (BSD-2).

**External data compliance**
- Drop E9.75-E13.5 cells (including E13.5 exactly and Theiler stages 16-22) from JAX, Feng, MOSTA and any atlas *before* any processing. Log the exact filters and disclose sources and stages.
- Do not use the challenge's own Virtual Embryo atlas at target stages.

**Compute** [Medium]
- All pipelines fit on one H100. Low-rank moscot on tens of thousands of cells in a 30-50-D latent takes minutes. scVI on about 100-200K cells takes tens of minutes. Flow matching in a 50-D latent takes minutes to hours.
- Workspace: DGX `~/v/vec-t1-scratch`; currently a minimal env only (no anndata, scanpy, scvi-tools or moscot yet).
- Caveat: the backtest gap (1 day, t ≈ 0.2) is shorter than the test gap (2 days, t ≈ 0.5-0.6). Prefer hyperparameters with a physical meaning (time fractions, shrinkage toward E10.5) over free knobs that may not transfer.

---

## 5. Emitting realistic single cells in 32K-gene log1p space (Q4)

**Facts** [High on the sample files]
- `.X = log1p(1e4·c/N)` exactly. Counts and per-cell depth N can be recovered losslessly.
- Real cells are about 88% zeros, with about 3.9-4.0K genes detected per cell.
- The smallest nonzero value is log1p(1e4/N), about 0.38-0.9 for the observed depths. Values between 0 and that level never occur in real data. Dense decoders produce exactly those values, which inflate E|x|^0.5 terms.

**Recommendations**
1. **Prefer real cells plus count-space edits** for every type present at E10.5: resampling, reweighting, binomial thinning, Poisson/NB top-up. In the probe they keep sparsity and the variogram close to real data.
2. **If you decode, sample counts.** Use scVI `posterior_predictive_sample` (or CFGen's discrete likelihood) at a depth distribution drawn from E10.5, then renormalise to 1e4 and log1p. Never submit `get_normalized_expression` means.
   - In the probe, dense means made the variogram 7.5x worse than the floor and worsened MMD.
   - Poisson resampling recovered most of the variogram loss. NB with learned dispersion should do better than Poisson, because real data are overdispersed.
   - Use NB rather than ZINB for UMI data: droplet UMI data are not zero-inflated (Svensson, Nat Biotechnol 2020).
3. **Hybrid "retrieve-and-edit" for new states**: find the nearest real cell (ours, or external mapped to our batch) and edit its counts toward the target mean, rather than decoding a mean.
4. **Calibrate pseudobulk ranks** by rescaling nonzero entries per gene; do not add small dense increments.
5. **Watch the gene sets.**
   - 8,522 Gm* and 2,524 *Rik genes, plus genes absent from JAX (24,552 genes), have no external signal. Hold them at E10.5 levels (dp ≈ 0, exact zeros kept). That keeps de_direction from absorbing noise.
   - Map genes across datasets by Ensembl ID.

---

## 6. Single-cell foundation models: utility and legality (Q5)

**Rule** [High]: "No measured data from a held-out stage … by any route: not directly, **not through a model pre-trained on it**, and not through a public dataset that contains it." For T1 the excluded window is after E9.5 through E13.5. Cells from a resource that spans the window must be removed *before training*. You cannot remove them from someone else's pretrained weights.

**What CELLxGENE contains** [High]
- From the CZ CELLxGENE curation API dataset listing, fetched 2026-10-18 (2,243 datasets; 499 mouse).
- 49 mouse datasets carry Theiler stages 16-22 (about E10-E13.5). Among them:
  - JAX time-lapse, **11,441,407 primary cells**, published 2024-03-11, with per-cluster subsets including "Cardiomyocytes" (37,205 cells).
  - LungMAP murine lung (TS21), 2021.
  - Developing thalamus (TS20 / TS22), 2022.
  - Mouse islet atlas (TS20-22), 2023.
  - Mouse gonad (TS17-22), 2023.
  - Cerebellum (TS18-22), 2023.
  - Limb (TS16-22) and craniofacial (TS17-22).
- So any model pretrained on CELLxGENE *mouse* data contains in-window cells, from any census since 2023, and since 2024-03 very likely the JAX embryos as well.

| Model | Pretraining corpus (verified) | In-window mouse data? | Verdict for T1 | Utility here |
|---|---|---|---|---|
| **scGPT** whole-human (Cui et al., Nat Methods 2024) | 33M normal **human** cells (CELLxGENE) | No mouse | **Allowed** (mouse window) [High]. Human fetal data are included; the rules do not address cross-species stage equivalence (ask organisers if used) | Low: human vocabulary; ortholog mapping drops Gm/Rik/lncRNA genes |
| **Geneformer** V1 / V2 (Theodoris et al., Nature 2023) | ~30M (V1) / ~104M (V2) **human** cells | No | **Allowed** [High] | Low: rank encoding, no expression decoder |
| **scFoundation** (Hao et al., Nat Methods 2024) | >50M **human** cells | No | **Allowed** [High] | Low |
| **UCE** (Rosen et al., bioRxiv 2023) | 36M cells, 8 species; mouse 23 datasets / 4.14M cells (supplementary list) | **Yes**: gonad male/female/germ (TS17-22), islet (TS20-22), developing thalamus (TS20/22), LungMAP lung (TS21); 726K cells in those datasets [High] | **Disallowed** | n/a |
| **TranscriptFormer** (Pearce et al., bioRxiv 2025) | TF-Sapiens 57M human; TF-Exemplar 110M (human + mouse, zebrafish, fly, worm); TF-Metazoa 112M (12 species); "27 million mouse cells (88 datasets)", mostly via the CELLxGENE Discover API; released download code filters primary data only, no stage filter | TF-Exemplar / Metazoa: **yes** (old datasets are in-window regardless). Including JAX is very likely: non-JAX primary mouse RNA on CELLxGENE was about 6.4M at early 2025, too little to reach 27M [Medium-High] | TF-Sapiens **allowed**; TF-Exemplar / Metazoa **disallowed** | Low |
| **Tabula** (Ding … Qiu, arXiv 2607.19400; organisers' lab) | Not read (PDF too large to fetch). Released vocabulary: 60,697 tokens, all human gene nomenclature (HGNC symbols, human ORFs, RP11-clones) | Likely none [Medium] | **Likely allowed**; confirm with organisers. Non-commercial licence | Low |
| **State** SE (Adduri et al., bioRxiv 2025, Arc) | SE "trained on 167 million **human** cells" | No | **Allowed** [High] | Low (built for perturbation in cell lines) |
| **scPRINT** (Kalfon et al., Nat Commun 2025) | 54M primary **human and mouse** cells from CELLxGENE "mid-2023 release" | Yes (LungMAP, thalamus, islet, gonad were all public by then) [Medium-High] | **Disallowed** | n/a |
| **Nicheformer** (Nat Methods 2025) | Census v2023-07-15 (human + mouse, 33.47M) + GEO + sfaira + HCA, + 53M spatial | Yes (same reasoning) [Medium] | **Disallowed** | n/a |
| **GeneCompass** (Cell Res 2024) | 48.2M mouse + 53.6M human cells (GEO/ArrayExpress/CNCB/CELLxGENE) | Very likely (Table S1 not checked) [Medium] | **Disallowed** | n/a |
| **Mouse-Geneformer** (PLOS Genet 2025) | 20.6M mouse cells incl. embryonic "E" stages | Very likely [Medium] | **Disallowed** | n/a |
| **C2S-Scale** (bioRxiv 2025) | >50M human + mouse cells (CELLxGENE, HCA) | Likely [Medium] | **Disallowed** | n/a |
| **Navigo** checkpoints / **Chreode** weights | JAX atlas, 43 time points / 2.4M-cell mouse embryo atlas | Yes / presumed | **Disallowed** (Navigo code BSD-2 is fine to retrain on legal stages) | n/a |
| Text-derived gene embeddings (e.g. GenePT) | Gene descriptions, no expression | No | **Allowed** [Medium] | Small: gene-similarity prior for unobserved genes |

Overall:
- The legal FMs are human-only, and the evidence that they help transcriptome-level prediction is weak (see the two perturbation papers).
- The task needs mouse heart- and embryo-specific temporal trends, which the legal bracketing data give directly.
- **Recommendation:** do not use FMs. If one is used, restrict to human-only checkpoints, disclose it, and ask the organisers about cross-species stage equivalence. [Medium]

---

## 7. Open questions

1. E12.5 dissection composition: what fraction is cardiac vs foregut, liver/ST or blood. This dominates mmd_u and de_score, and no legal data observes it.
2. Do organisers count in-window mouse cells of *other tissues* inside a pretrained model as a violation? A strict reading says yes. Ask on Slack.
3. Is human fetal data inside human-only FMs acceptable (cross-species staging)? The rules are silent.
4. Confirm usage of Feng 2022 E9.5 and E14.5+ (chamber dissections) and JAX E13.75+. Is "after E13.5" satisfied by E13.75?
5. Are the veckit sample files true subsamples of the released stages? My count-recovery and cell-vs-nucleus findings come from them.
6. E12.5 truth size. The 360 MB file vs 647 MB for E10.5 (about 25K cells implied) suggests fewer cells. Test max_cells comes with P3.
7. Best time-warp for bracketing: linear in days, or a developmental-clock warp.

---

## 8. References

**Challenge**
- Challenge site and rules: https://virtualembryo.ai/challenge/rules ; data: https://virtualembryo.ai/challenge/data ; T1 leaderboard API: https://api.virtualembryo.ai/challenge/leaderboard?task=T1
- veckit scorer (MIT): https://github.com/aristoteleo/veckit

**OT, dynamics and forecasting methods**
- Schiebinger G. et al. Waddington-OT. *Cell* 176:928 (2019). https://doi.org/10.1016/j.cell.2019.01.006 ; code https://github.com/broadinstitute/wot
- Klein D. et al. moscot. *Nature* 638:1065 (2025). https://doi.org/10.1038/s41586-024-08453-2 ; code https://github.com/theislab/moscot
- Tong A. et al. TrajectoryNet. ICML 2020. https://arxiv.org/abs/2002.04461 ; code https://github.com/KrishnaswamyLab/TrajectoryNet
- Huguet G. et al. MIOFlow. NeurIPS 2022 (arXiv 2206.14928, from memory, not re-verified). Code https://github.com/KrishnaswamyLab/MIOFlow
- Yeo G.H.T. et al. PRESCIENT. *Nat Commun* 12:3222 (2021). https://doi.org/10.1038/s41467-021-23518-w ; code https://github.com/gifford-lab/prescient
- Zhang J. et al. scNODE. *Bioinformatics* 40(Suppl 2):ii146 (2024). https://doi.org/10.1093/bioinformatics/btae393 ; code https://github.com/rsinghlab/scNODE
- Jiang Q., Wan L. PI-SDE. *Bioinformatics* 40(Suppl 2):ii120 (2024). https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11373338/ ; code https://github.com/QiJiang-QJ/PI-SDE
- Jiang Q. et al. scIMF. arXiv 2505.16492. https://arxiv.org/abs/2505.16492 ; code https://github.com/QiJiang-QJ/scIMF
- Sha Y. et al. TIGON. *Nat Mach Intell* (2023/24). https://doi.org/10.1038/s42256-023-00763-w ; code https://github.com/yutongo/TIGON
- Zhang Z., Li T., Zhou P. DeepRUOT. ICLR 2025. https://arxiv.org/abs/2410.00844 ; code https://github.com/zhenyiizhang/DeepRUOT
- Wang D. et al. VGFM. NeurIPS 2025. https://arxiv.org/abs/2505.13413 ; code https://github.com/DongyiWang-66/VGFM
- Tong A. et al. OT-CFM / torchcfm. TMLR 2024. https://arxiv.org/abs/2302.00482 ; code https://github.com/atong01/conditional-flow-matching
- Rohbeck M. et al. MMFM. ICLR 2025. https://proceedings.iclr.cc/paper_files/paper/2025/file/206b764b00195b8c896c5af8bde80726-Paper-Conference.pdf ; code https://github.com/Genentech/MMFM
- Lee J. et al. MMSFM. ICML 2025. https://proceedings.mlr.press/v267/lee25x.html
- Petrović K. et al. Curly-FM. NeurIPS 2025. https://arxiv.org/abs/2510.26645
- Korotin A. et al. Light SB. ICLR 2024. https://arxiv.org/abs/2310.01174
- Chen T. et al. DMSB. NeurIPS 2023. https://arxiv.org/abs/2303.01751
- Neklyudov K. et al. Action Matching. ICML 2023. https://proceedings.mlr.press/v202/neklyudov23a.html
- CytoBridge: https://github.com/zhenyiizhang/CytoBridge
- Fan Y. et al. Navigo. bioRxiv 2026. https://www.biorxiv.org/content/10.64898/2026.06.18.733286v1 ; code https://github.com/aristoteleo/Navigo-release
- GeneSys. bioRxiv 2025. https://www.biorxiv.org/content/10.1101/2025.08.20.671385
- Chreode. arXiv 2605.28111. https://arxiv.org/abs/2605.28111

**Generative models and decoders**
- Luo E. et al. scDiffusion. *Bioinformatics* 40:btae518 (2024). https://doi.org/10.1093/bioinformatics/btae518
- He S. et al. Squidiff. *Nat Methods* 23:65 (2026). https://doi.org/10.1038/s41592-025-02877-y
- Palma A. et al. CFGen. ICLR 2025. https://arxiv.org/abs/2407.11734 ; code https://github.com/theislab/CFGen
- Klein D. et al. CellFlow. bioRxiv 2025. https://doi.org/10.1101/2025.04.11.648220
- Lopez R. et al. scVI. *Nat Methods* 15:1053 (2018). https://doi.org/10.1038/s41592-018-0229-2
- Svensson V. Droplet scRNA-seq is not zero-inflated. *Nat Biotechnol* 38:147 (2020). https://doi.org/10.1038/s41587-019-0379-5

**Vector fields and fate mapping**
- Qiu X. et al. dynamo. *Cell* 185:690 (2022). https://doi.org/10.1016/j.cell.2021.12.045
- Lange M. et al. CellRank. *Nat Methods* 19:159 (2022). https://doi.org/10.1038/s41592-021-01346-6

**Benchmarks**
- Osakwe A. et al. scTimeBench. *Bioinformatics* (2026). https://academic.oup.com/bioinformatics/article/42/Supplement_2/btag414/8767264 ; preprint https://www.biorxiv.org/content/10.64898/2026.03.16.712069v1 ; code https://github.com/li-lab-mcgill/scTimeBench
- Ahlmann-Eltze C., Huber W., Anders S. *Nat Methods* 22:1657 (2025). https://doi.org/10.1038/s41592-025-02772-6
- Kedzierska et al. *Genome Biol* (2025). https://doi.org/10.1186/s13059-025-03574-x

**External data**
- Qiu C. et al. JAX time-lapse. *Nature* 626:1084 (2024). https://doi.org/10.1038/s41586-024-07069-w
- Feng W. et al. Mouse heart E9.5-P9. *Nat Commun* 13:7960 (2022). https://doi.org/10.1038/s41467-022-35691-7 (GEO GSE193346)
- CZ CELLxGENE curation API (dataset listing used for the contamination analysis): https://api.cellxgene.cziscience.com/curation/v1/datasets

**Foundation models**
- Cui H. et al. scGPT. *Nat Methods* 21:1470 (2024). https://doi.org/10.1038/s41592-024-02201-0
- Theodoris C.V. et al. Geneformer. *Nature* 618:616 (2023). https://doi.org/10.1038/s41586-023-06139-9 ; model card https://huggingface.co/ctheodoris/Geneformer
- Hao M. et al. scFoundation. *Nat Methods* 21:1481 (2024). https://doi.org/10.1038/s41592-024-02305-7
- Rosen Y. et al. UCE. bioRxiv 2023. https://www.biorxiv.org/content/10.1101/2023.11.28.568918 (supplementary dataset list media-3.xlsx)
- Pearce et al. TranscriptFormer. bioRxiv 2025. https://www.biorxiv.org/content/10.1101/2025.04.25.650731v2 ; code https://github.com/czi-ai/transcriptformer
- Adduri et al. State. bioRxiv 2025. https://www.biorxiv.org/content/10.1101/2025.06.26.661135v2
- Kalfon J. et al. scPRINT. *Nat Commun* (2025). https://doi.org/10.1038/s41467-025-58699-1
- Nicheformer. *Nat Methods* (2025). https://doi.org/10.1038/s41592-025-02814-z
- GeneCompass. *Cell Res* (2024). https://doi.org/10.1038/s41422-024-01034-y
- Mouse-Geneformer. *PLOS Genet* (2025). https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11964219/
- C2S-Scale. bioRxiv 2025. https://www.biorxiv.org/content/10.1101/2025.04.14.648850
- Tabula: https://github.com/aristoteleo/tabula ; arXiv https://arxiv.org/abs/2607.19400

**Local artefacts from this session**
- Probe script: `SCRATCH/scripts/t1_probe.py`
- Nuclear-fraction check: `.../scratchpad/scripts/t1_nuclear_check.py`
- CELLxGENE listing: `.../scratchpad/dl_cxg/datasets.json`
- UCE dataset list: `.../scratchpad/dl_uce/media3.txt`
- DGX workspace: `~/v/vec-t1-scratch` (uv env `envs/t1`: numpy 2.4.6, scipy 1.17.1, scikit-learn 1.9.1, h5py 3.16.0)
