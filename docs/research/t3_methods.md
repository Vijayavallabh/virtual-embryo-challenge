# Task 3 methods: predicting the unseen β-catenin knockout (Mesp1-Cre; Ctnnb1 F/F, E8.75)

Report date: 2026-10-18 (P2; P3 opens 2026-10-20; final deadline 2026-12-02). Scope: Human Team track only.
Dimension: methods for predicting an unseen knockout response from wild-type development plus 1-2 observed knockouts.

Confidence labels: **[High]** verified in primary source or official code/data; **[Medium]** well-supported inference or single
reliable source; **[Low]** plausible hypothesis, not verified; **[Unverified]** could not check. Synthetic results are labelled
**[Synthetic]**. Public material written by other competitors is labelled **[Untrusted]**.

Sibling reports in the same directory cover the biology and the rule text in more depth
(`t3_biology.md`, `landscape.md`, `architecture.md`); this report is about methods, metric mechanics, and validation.

---

## 0. Bottom line

1. **What the T3 metric rewards is the population-level, pseudobulk change versus one specific wild-type sample, measured on
   sparse 500-gene MERFISH.** On this data almost every pseudobulk change is a change in detection rate, so composition
   (which cell types, how many), sequencing/imaging depth and section coverage dominate the differentially expressed (DE)
   gene set. Subtle within-type programme changes matter mostly for `de_direction`. **[Medium]**
2. **Why `shift_transfer` fails is largely mechanical, not just biological.** Adding a pseudobulk delta to every log-normalised
   cell and clipping at zero realises only ~4% (median gene) of an intended *decrease* but 100% of an *increase*, because ~93%
   of MERFISH entries are zero; the realised down-ranking then correlates 0.99 with wild-type expression level, which is
   exactly the de_score null and exactly what de_direction partials out. Up-shifts fill zeros in every cell, which is why its
   variogram is ~5x worse than doing nothing. On top of that, the source (Mab21l2, E9.5, different wild-type sample, different
   tissues lost) shares none of the target's composition or reference-sample idiosyncrasies. **[High]** for the
   clipping arithmetic (measured on the 150-cell veckit sample and from the leaderboard); **[Medium]** for the rest.
3. **The public board already shows what works on Gata4 (P2):** 179 of 200 teams beat the floor; best human 80.27; published
   baselines 42.5. A public competitor log shows that mixing *real knockout cells* (even Mab21l2 cells from E9.5) into the
   wild type scores ~64, while adding the same knockout's pseudobulk delta scores ~41-42 (organisers' `shift_transfer`
   42.48; the competitor's ×1.5 variant 40.98). **[High]** for the leaderboard numbers; **[Untrusted]** for the competitor log.
4. **Literature consensus that simple baselines rival deep models is strong** (Ahlmann-Eltze et al. 2025; Csendes et al. 2025;
   Kernfeld et al. 2023; Systema 2025), and the Arc Virtual Cell Challenge 2025 was won by hybrid statistical/deep systems that
   optimised the metrics directly. With one or two observed knockouts, GEARS/CPA/scGen/CellOT/biolord/PerturbNet/State cannot
   be trained meaningfully; foundation-model in-silico deletion of *Ctnnb1* mRNA is mis-specified because β-catenin is
   regulated post-translationally. **[High]** for the literature; **[Medium]** for the T3-specific verdicts.
5. **Recommended core pipeline (P1 + P3 + P2, combined by P4):** (a) carry the knockout-series "systematic" component over from
   the released Gata4 knockout (same stage, same wild-type reference, same Cre driver) using real cells, not additive shifts;
   (b) model Mesp1-lineage composition changes explicitly, with sign hedged across the two test submissions; (c) add a
   within-type Wnt-responsive programme learned only from wild-type T1 scRNA (Wnt *target* activity, not *Ctnnb1* mRNA);
   (d) weight the parts by leave-one-knockout-out validation with the official scorer. **[Medium]**
6. **The single biggest binary risk is `severity_slope`** (25%): if the predicted change on the truth's DE genes has slope ≤ 0
   or R² ≤ 0.01 it falls to the floor (skill 0.5 instead of ≥0.93), costing ~12 T3 points. Getting the sign of the DE-gene
   response right matters more than its fine detail. **[High]** (from scorer code).
7. **Day-1 (2026-10-20) diagnostics decide the strategy:** depth and detection per sample, harmonised composition per sample,
   non-Mesp1 compartments in Gata4 replicate 1 vs 2 vs wild type, whether knockout files contain littermate wild-type cells,
   and an "oracle decomposition" of the Gata4 truth into composition, depth/coverage and within-type parts, each scored with
   veckit. **[High]** that these are cheap and decisive.
8. **Legality:** everything in P1-P4 uses only released data plus generic, provenance-filtered resources. Published Ctnnb1
   loss-of-function work in cardiac mesoderm at E7.5-E10.5 ([redacted 2026-10-19: B]
   [redacted 2026-10-19: B]) and LLM priors that may reproduce it are a real rule risk under "another allele of the same gene at a
   comparable stage" and "not through a model pre-trained on it". Ask the organisers before using them (questions in §7).
9. **Two test submissions only:** pre-register the candidates, make submission 2 differ from submission 1 in the most uncertain
   high-impact choice (Mesp1-lineage composition sign/magnitude), and submit near the board's maximum cell count to minimise
   pseudobulk noise.

---

## 1. What the T3 scorer rewards, and what the public board already tells us

### 1.1 Metric mechanics and the design rules they imply

Source: `ext/veckit/common/core_metrics.py`, `T3/metrics_v2.py`, proposal Table 1, `index.json`. **[High]**

| Metric (weight) | Exact computation (T3) | Design implication |
|---|---|---|
| `de_score` (0.30) | Truth DE genes = abs(pseudobulk logFC) ≥ 0.25 vs matched WT (Mann-Whitney + BH is inert; the docstring reports all [removed 2026-11-05: a possible property of a hidden target; docs/DATA_POLICY.md 2026-11-05 (h)] genes identical to the abs(lfc) rule on a T3 board). Prediction contributes its own top-n_up and bottom-n_dn genes by dp = pb(pred) − pb(WT_ref); signed overlap; chance = best of ranking by +pb(WT) or −pb(WT). | Only gene *ranks* of dp matter. You must beat an expression-level null: predicting "highly expressed genes go up/down" scores 0. Down-gene predictions that merely track WT expression level are worthless. |
| `de_direction` (0.25) | Rank partial correlation of dp and dt over all 500 genes, controlling for rank of pb(WT_ref). | Small, non-DE changes count. Anything collinear with expression level is removed. |
| `severity_slope` (0.25) | log of zero-intercept OLS slope β of dp on dt over truth DE genes; worst value log(1e-3) = −6.91 if β ≤ 0 or R² ≤ 0.01; scored on abs(value). | Lenient on magnitude (abs(log β) = 0.5 still gives skill ≈ 0.93), brutal on sign. Rescaling dp does not change de metrics, so set the magnitude for severity. |
| `mmd_u` (0.12) | Unbiased multi-bandwidth RBF MMD² on 2,000 cells in a 30-PC space fitted on the truth. | Needs a realistic cell population (composition + per-type states + sparsity). |
| `variogram` (0.08) | 20,000 random gene pairs, 1,500 cells, E abs(x_i − x_j)^0.5 compared. | Punishes destroyed sparsity / co-expression (e.g., additive shifts). |
| Skill map | skill = d_floor / (d_floor + d(value)); floor = wt_identity (0.5), ceiling = half-vs-half of the truth replicate (1.0); NaN → 0. | Floor and ceiling are measured per board. The ceiling contains no embryo-to-embryo variation, so it is optimistic for any cross-embryo prediction. |

Scorer protocol: every stage is subsampled to 10% before scoring, then the knockout target is split in half. On the Gata4
board: `truth_cells` 2,320, `ref_cells` 2,483, cells allowed 1,000-7,449 (`index.json`). So the scored Gata4 replicate holds
~46,400 cells and the WT E8.75 reference ~24,830 cells; **knockout samples contain ~1.9x more cells than the WT reference**
(consistent with different section coverage; the cause is unknown). **[High]** for the numbers, **[Low]** for the cause.

Practical rules that follow:
- **Generate cells, do not add deltas.** Work in count space (thinning/augmenting counts, then `log1p(1e4·c/Σc)`) or resample
  real cells. The MERFISH files carry `layers["counts"]` (integers) and `.X` is exactly `log1p(1e4·c/Σc)` over the 500
  panel genes (expm1 row sums = 10,000 in all three 150-cell samples). **[High]**
- **Submit close to the maximum cell count** (7,449 on Gata4; test limits not yet published) and build from the full WT pool
  (all ~24.8K cells), never a 10% subsample. A public competitor log **[Untrusted]** reports de_score skill ≈ 39 (below the
  floor) for re-uploaded WT cells, which is what you expect when a submission's own sampling noise is scored as a
  "response" against the expression-level null. **[Medium]** for the explanation.
- Per-cell normalisation over only 500 genes couples genes: if a few highly expressed genes fall in a cell, every other gene's
  normalised value rises. Simulate in counts and re-normalise.

### 1.2 The Gata4 validation board (P2), fetched 2026-10-18 15:43 UTC

Source: public endpoint `https://api.virtualembryo.ai/challenge/leaderboard?task=T3` (the site's leaderboard data). Raw metric
values below were back-computed from published skills with the published anchors (this is about the scale of public scores,
not about the held-out truth). **[High]**

| Row | Headline | de_score skill | de_dir skill | severity skill | mmd skill | variogram skill | Implied raw values |
|---|---|---|---|---|---|---|---|
| ceiling | 100 | 100 | 100 | 100 | 100 | 100 | de 0.870, dir 0.925, sev −0.009 |
| wt_identity (floor) | 50 | 50 | 50 | 50 | 50 | 50 | mmd 0.0314, vario 0.0322 |
| shift_transfer | 42.48 | 40.4 | 48.7 | 50.0 | 37.1 | 15.5 | de ≈ −0.41, dir ≈ −0.05, sev at floor, mmd ≈ 0.053, vario ≈ 0.17 |
| gene_ko | 42.50 | 40.4 | 48.7 | 50.0 | 37.2 | 15.5 | as above (zeroing one gene changes little) |
| best agent row | 80.51 | 70.2 | 70.2 | 100 | 79.9 | 91.4 | de ≈ 0.50, dir ≈ 0.53 |
| best human row (team "Bat") | 80.27 | 74.1 | 72.7 | 99.5 | 71.0 | 80.8 | de ≈ 0.57, dir ≈ 0.58, abs(sev) ≈ 0.04 |

Distribution of 200 teams (148 human, 52 agent): median 65.7, 90th percentile 75.6, 179 above 50.

How to read this:
- P2 allowed 8 scored submissions per day for weeks on the Gata4 board, so top P2 scores include heavy leaderboard
  selection on the very target they are scored against. **The test board allows 2 submissions in total**, so P2 ranks are a
  poor guide to test performance. **[High]** (rules §11) / **[Medium]** (interpretation).
- Model names on the board hint at what people found: `t3cmp_fewerCM0.85`, `union_fewerCM_thin0.7`, `t3_lineage_P1null`,
  `Composition r1_..._control_cova`, `R4 STRING11.5-GATv2/PyG`. Composition (fewer cardiomyocytes) and lineage restriction
  recur. **[Medium]**
- A public, self-described methods journal by one human team (`wqty123/virtual-embryo-challenge-methods`,
  `experiments_t3.en.md`, last commit 2026-09-23) reports on the Gata4 board **[Untrusted]**:
  global Δ transfer 40.98 (variogram skill 49.8 → 12.2); restricting Δ to cardiac cells 44.4-44.9; a "Gata4 state programme"
  ≈ floor; **70% real Mab21l2 knockout cells + 30% WT E8.75 cells 64.6** (de 53.3, dir 59.8, severity 92.3); a 40% "early
  carrier" mix 50.1 (dir 58.4 but severity at the floor); cell-type-flattened resampling 65.3; any per-gene transform family
  collapsed the variogram (10.8-27.1). They also report that de/dir/severity are deterministic across scorer seeds while
  mmd_u varies by ±5.75 and variogram by ±2.2 skill points.
  The contrast "same knockout, real cells ≈ 64 vs additive delta ≈ 41" is the single most useful public fact for T3.

### 1.3 Facts from the official 150-cell samples (illustrative only)

Measured with h5py on `ext/veckit/data/` (scripts in `scratchpad/scripts/`). Tiny samples, possibly non-random; use only as
orientation. **[Medium]**

| Quantity | WT sample (heart vocabulary) | WT heart E9.5 sample | Mab21l2 KO sample (E9.5) |
|---|---|---|---|
| Counts per cell, median (p10-p90) | 75 (38-153) | 69 (39-155) | **102 (43-289)** |
| Genes detected per cell, median | 34 | 35 | 43 |
| Fraction of zeros in `.X` | 0.926 | 0.925 | 0.903 |
| Mab21l2 detected | 13% of cells | 23% | **0%** |
| Cell-type vocabulary | V-CM, A-CM, OFT-CM, SV-CM, Chamber-Endo, Peri, Dorsal-PAM, Lateral-PAM, aPHM, pPHM, Proepi, ST, BEC, ... (22 types) | 20 types, same vocabulary | vCM, OFT/RV-CM, SV-CM, Endocardium, Pericardium, DPM, LPM, Paraxial Mesoderm, Branchial Arch, Hepatocyte, ... (14 types) |
| Extra obs columns | celltype, cm_celltype | same | genotype, perturb_batch (0), sample_id, section, slide (M125), time, z_coord, leiden |

Consequences: (i) knockout and WT files use different annotation vocabularies, so any cell-type-resolved model needs a label
harmonisation step (the scorer's own T3 code notes this). (ii) The knockout sample is deeper per cell than the WT samples;
if that holds in the full files it is a first-order, predictable nuisance (§4.1). (iii) No Mab21l2 detection in the knockout
sample suggests the probed transcript is removed in all sampled cells, i.e. a germline or very broad deletion; the data page
only states that the Gata4 and β-catenin knockouts are Mesp1-Cre. Verify on the full files.

Clipping arithmetic on the WT sample: per-gene detection median 4.0% (p90 22.7%); non-zero log values 3.61-8.46 (median
5.07); "add −Δ then clip at 0" realises a median 4.0% of the intended decrease (p90 22.7%), "add +Δ" realises 100%; the
realised decrease magnitudes correlate r = 0.992 with WT pseudobulk level. **[High]** (`scripts/clip_attenuation.py`).

---

## 2. Q1: survey of perturbation-response methods for an unseen gene with one or two observed knockouts

### 2.1 Methods table

"Unseen-gene ready" = can produce a prediction for a gene never perturbed in training. "T3 fit" judges usefulness here,
where there are 1-2 observed knockouts, a 500-gene panel, and a gene (*Ctnnb1*) that is not on the panel.

| Method (ref) | Core idea | Needs | Unseen-gene ready? | T3 fit and verdict |
|---|---|---|---|---|
| GEARS (Roohani et al., Nat Biotechnol 2024) | GNN over a GO gene-gene graph; predicts post-perturbation mean | Perturb-seq with many perturbations; perturbed gene must be in the graph | Yes, via GO neighbours | **Low.** Cannot be trained on 2 knockouts; the GO-graph idea survives only as a prior. `cell-gears` 0.1.2 (2023-12). |
| CPA (Lotfollahi et al., Mol Syst Biol 2023) | Autoencoder with additive latent perturbation/covariate embeddings | Many perturbations; chemical embeddings for unseen drugs | Not for unseen genes (stated by Ahlmann-Eltze 2025 and Systema) | **Low.** `cpa-tools` 0.8.8. |
| scGen (Lotfollahi et al., Nat Methods 2019) | VAE latent arithmetic: transfer a perturbation vector across cell types | The *same* perturbation observed in another context | No | **Low** for the test gene. Conceptually `shift_transfer` in latent space; useful only to transfer *shared* (systematic) components. `scgen` 2.1.0; also in `pertpy` 1.4.0. |
| CellOT (Bunne et al., Nat Methods 2023) | Neural OT map control → perturbed per perturbation | Observed perturbed cells for that perturbation | No | **Low** for biology; an OT map WT → knockout *series* learned on non-Mesp1 cells could transfer technical shifts. |
| biolord (Piran et al., Nat Biotechnol 2024) | Disentangles known attributes; shifts cells across attribute values | Attributes with several observed values | Partly (needs attribute embedding) | **Low-medium.** Could disentangle stage vs genotype vs batch across 5 WT stages + 2 KOs; with 2 KOs the genotype attribute is barely identified. `biolord` 0.0.3. |
| PerturbNet (Yu et al., Mol Syst Biol 2025) | Conditional generative model mapping perturbation representations (GO/text/sequence) to cell-state distributions | Many perturbations for the mapping | Yes | **Low** (no training set of gene knockouts in this system). |
| State (Adduri et al., bioRxiv 2025) | Set-level transformer over cell populations; trained on >100M perturbed cells | Large perturbation corpora (human cell lines) | Within its training distribution | **Low.** Human cell-line perturbations; no mouse embryo knockouts. |
| TxPert (Wenkel et al., arXiv 2505.14919, 2025) | Several knowledge graphs (STRING, GO, phenomics/transcriptomics maps) for OOD perturbations | Large perturbation screens | Yes | **Low-medium** as a *prior source* (graph proximity to *Ctnnb1*). |
| Squidiff (He et al., Nat Methods 2025/2026) | Conditional diffusion model for development and perturbation | Training conditions | Interpolation between seen conditions | **Low** for the unseen gene; possible generator for realistic cells. |
| scGPT perturbation fine-tuning (Cui et al., Nat Methods 2024) | Transformer pre-trained on 33M human cells; fine-tuned on Perturb-seq | Perturb-seq fine-tuning data; human genes | Claimed | **Low.** Human-only pretraining; benchmarks show no gain over simple baselines. `scgpt` 0.2.4. |
| Geneformer / Mouse-Geneformer in-silico deletion (Theodoris et al., Nature 2023; Ito et al., PLoS Genet 2025: 21M mouse cells) | Delete a gene token, measure embedding shift | Only WT cells | Yes | **Low.** Deleting the *Ctnnb1* token models co-expression with *Ctnnb1* mRNA, which mostly reflects adhesion/epithelial state, not Wnt output (β-catenin is regulated post-translationally). Also an unknown pretraining corpus (rule risk if it contains Ctnnb1 loss-of-function data at comparable stages). |
| CellOracle (Kamimoto et al., Nature 2023) | Cluster-wise ridge GRN on a motif/ATAC base GRN; in-silico TF KO propagated 3 steps (`n_propagation=3`) | WT scRNA + base GRN (mouse: mm9 sci-ATAC atlas base GRN, mm10/mm9 promoter base GRNs) | Only for TFs in the base GRN | **Medium as a module.** Source check: `simulate_shift` raises "not TF or TF motif information is not available" for non-TF genes, so *Ctnnb1* cannot be perturbed directly; use TCF/LEF ([redacted 2026-10-19: B]) as entry points, or better a Wnt-activity latent variable (§5, P2). Directly testable on Gata4/Gata6 (TFs). `celloracle` 0.20.0 (2024-11). |
| dynamo (Qiu et al., Cell 2022) | Vector field from RNA velocity/labelling; perturbation response ΔY = J·ΔX via the analytical Jacobian; PCA basis only | Velocity (spliced/unspliced or labelling) | Yes | **Low-medium.** T1 h5ad has no spliced/unspliced layers (checked on the T1 sample: `layers` empty), so a velocity field would have to come from temporal OT/flow matching instead. `dynamo-release` 1.5.3 (2026-03). |
| Navigo (Qiu lab, GitHub `aristoteleo/Navigo-release`, BSD-2) | Flow matching + RNA kinetics on mouse embryogenesis atlases; zero-shot KO | Their checkpoints/atlas | Yes (zero-shot) | **Low-medium.** Likely the organisers' `dynode_flow` reference; that row is absent from the T3 board. Check checkpoint training data for rule issues before use. |
| scTenifoldKnk (Osorio et al., Patterns 2022) | Virtual KO: build GRN from WT scRNA, delete a gene's edges, manifold alignment → differentially regulated genes | WT scRNA only | Yes | **Low-medium.** Same mis-specification as Geneformer for *Ctnnb1* mRNA; output is a ranking without a reliable sign. Cheap cross-check. `scTenifoldpy` 0.5.1. |
| SCENIC+ (Bravo González-Blas et al., Nat Methods 2023) / pySCENIC | Enhancer-driven eGRN from multiome; regulon activity; perturbation simulation | Multiome (RNA+ATAC) or cisTarget databases | TFs only | **Medium (pySCENIC regulons on T1 RNA)** to score TCF/LEF regulon activity; full SCENIC+ needs ATAC, which was not released. `pyscenic` 0.12.1. |
| LINGER (Yuan & Duren, Nat Biotechnol 2024) | GRN from single-cell multiome with atlas-scale external bulk data | Multiome; external bulk compendium | TFs | **Low.** Needs ATAC; built around human external data (mouse support not verified). |
| GenePT (Chen & Zou, bioRxiv 2023) | GPT text embeddings of NCBI gene summaries as gene features | None | Yes (as features) | **Low** alone; a weak, unsigned similarity prior. |
| LangPert (Märtens et al., ICLR 2025 MLGenX workshop) | LLM chooses informative training perturbations; kNN aggregates their measured responses | Training perturbations to aggregate | Yes | **Low** (only 2 training knockouts to aggregate); the LLM-as-prior idea carries rule risk (§7). |
| Knowledge-graph / pathway priors: CollecTRI (Müller-Dott et al., NAR 2023), PROGENy (Schubert et al., Nat Commun 2018), NicheNet (Browaeys et al., Nat Methods 2020), STRING | Curated TF→target, pathway-responsive genes, ligand→target potential | None | Yes | **Medium** as priors, after provenance filtering (§7). CollecTRI lists 40 panel genes as CTNNB1 targets; some edges cite the held-out-genotype literature (those edges must be dropped; gene names redacted 2026-10-18), and signs are noisy ([redacted 2026-10-19: B]). `decoupler` 2.2.0 `op.collectri/op.progeny(organism=...)`. |

All method metadata from PubMed/journal pages or the cited repositories; package versions from PyPI JSON on 2026-10-18. **[High]**

### 2.2 Evidence that linear or mean baselines rival deep models

- **Ahlmann-Eltze, Huber, Anders, Nat Methods 22:1657 (2025).** Five foundation models (scGPT, scFoundation, scBERT,
  Geneformer, UCE) plus GEARS and CPA versus deliberately simple baselines. For double perturbations (Norman K562), all
  models had higher error than the *additive* baseline; none beat "no change" at finding genetic interactions. For unseen
  single perturbations (Replogle K562/RPE1, Adamson), no deep model consistently beat the *training mean* or a linear model
  `Y ≈ G W Pᵀ + b` with PCA embeddings. The best model was the same linear model with perturbation embeddings pre-trained on
  *another cell line's perturbation data*. Pre-training on atlas data gave little; pre-training on perturbation data helped.
  **[High]**
- **Csendes et al., BMC Genomics 26:393 (2025).** The mean of training examples outperformed scGPT and scFoundation; simple
  ML with biologically meaningful features beat scGPT by a large margin; Perturb-seq benchmarks have low
  perturbation-specific variance. **[High]**
- **Kernfeld et al., bioRxiv 10.1101/2023.07.28.551039.** 11 datasets; on simple metrics such as MSE it is uncommon for
  expression-forecasting methods to beat simple baselines; conclusions depend strongly on the metric. **[High]** (preprint).
- **Systema (Viñas Torné et al., Nat Biotechnol 44:1050, 2025).** Shows that high scores of "perturbed mean" baselines and
  of published models come from *systematic variation*: consistent differences between perturbed and control cells caused by
  panel selection, confounders or generic biology (for example, 46% vs 25% of cells in G1 in perturbed vs control RPE1).
  Reference-based metrics (Pearson of deltas against control) are inflated by it; with the perturbed centroid as reference,
  all methods score low. **[High]**

**What this means for T3.** All T3 response metrics are reference-based (WT as reference), so the T3 analogue of Systema's
"systematic variation" (any shift shared by knockout samples relative to the WT reference: depth, section coverage, batch,
staging, generic mutant effects) is *rewarded* by the scorer. The T3 analogue of the "perturbed mean" baseline is the
released Gata4 knockout population itself (same stage, same WT reference, same Cre driver). It is the baseline to beat, and
it should be part of any serious submission. Deep perturbation models add little when only 1-2 perturbations exist.

### 2.3 Lessons from the Arc Virtual Cell Challenge 2025 (VCC)

Sources: Roohani et al., Cell 188:3370 (2025); Arc wrap-up (arcinstitute.org/news/virtual-cell-challenge-2025-wrap-up);
a participant post-mortem (yraj.in). **[High]** for the facts below unless marked.
- Data: H1 hESC, 300 CRISPRi perturbations, ~300,000 cells (~1,000 per perturbation), 10x Flex. Metrics: differential
  expression score (DES), perturbation discrimination score (PDS, L1-based), MAE.
- Winners: 1st BM_xTVC (BioMap; improved scFoundation encoder + protein-model embeddings + training-set DEG frequency/mean
  expression features; loss fused PDS + DES + small MAE; trained on pseudobulk; *repeated outputs to satisfy DES*); 2nd
  XLearning Lab (fully connected net on pseudobulk with ESM-2 embeddings; residual delta learning; public Perturb-seq only);
  3rd Team Outlier (TransPert: purely statistical, summary-level, similarity-weighted aggregation from reference cell lines,
  *global linear scaling tuned for PDS*); Generalist Prize Altos Labs (flow matching in gene space, pretrained on ~7M cells).
- Metric pitfalls the organisers named: almost all models were worse than baseline on MAE, which then stopped shaping
  optimisation; PDS rises with prediction magnitude (participants exploited it with scaling); a Generalist Prize across seven
  metrics was added mid-competition rather than changing the rules; "perturbation prediction models are not yet consistently
  outperforming naive baselines across all metrics".
- Transfer to T3: (i) hybrids of simple statistics and modest learned parts win; (ii) loss design aligned with the metric
  matters as much as architecture, so optimise and validate with veckit's exact functions; (iii) VEC already closed several
  VCC-style holes (rank-based de_score with an expression-level null, sign-split matching, partial correlation, slope-based
  severity, adversarial control rows such as `ctrl_one_cell` and `ctrl_random_dir`), so repetition/scaling tricks do not
  help here; the open channel is systematic variation; (iv) do not switch approach days before the deadline (the
  post-mortem author lost 68 places that way; anecdotal).

---

## 3. Q2: constraints and what they imply

| Constraint | Consequence | What to do |
|---|---|---|
| *Ctnnb1* is not on the 500-gene panel; classic direct targets ([redacted 2026-10-19: B]) are absent. Present Wnt-related genes: [redacted 2026-10-19: B]. **[High]** (panel check) | `gene_ko` is a no-op for the test (nothing to zero). No direct readout of deletion or Wnt activity in MERFISH. | Infer Wnt activity in WT from T1 scRNA (which has all of these genes) and transfer per-cell-type programmes to the panel via label mapping. Predict panel genes indirectly. |
| Deletion restricted to the Mesp1 lineage. The widely used Mesp1-Cre (Saga) is the knock-in `Mesp1<tm2(cre)Ysa>` (RIKEN RBRC01145, listed as MGI:2176467; homozygous lethal), so Cre carriers are Mesp1 heterozygous; the challenge does not name the allele, so this is presumed. Saga et al. 1999 traced Mesp1 cells to amnion mesoderm first, then mainly heart-tube myocardium (not endocardium in that study) and head mesenchyme. **[High]** for the allele/lineage text, **[Medium]** that this allele is the one used; whether the matched WT carries Cre is not stated. | Non-Mesp1 compartments (foregut/hepatic endoderm, surface ectoderm, neural tube, NCC, much of the endothelium/blood) should not change cell-autonomously. Their changes in a knockout sample estimate nuisance (depth, coverage, batch, staging) plus non-autonomous effects. | Treat non-Mesp1 cells as internal controls (RUV-style "negative-control compartment") to estimate the knockout-series nuisance. Restrict biological changes to harmonised Mesp1 types. |
| Mab21l2 KO is at E9.5 with its own WT E9.5 reference; its allele is not stated (sample suggests broad deletion). Mab21l2 loss affects myocardium, septum transversum, proepicardium and hepatoblasts (Saito et al., PLoS One 2012). **[High]** for the paper; **[Medium]** for the allele inference. | Its delta mixes different tissues (including endoderm), a different stage and a different reference sample. Little transfers gene-by-gene; it is still useful for hyperparameters (typical magnitudes) and as a leave-one-out fold. | Use for the LOKO fold Gata4 → Mab21l2 and to estimate knockout-series nuisance at E9.5; never as a direct delta. |
| Gata4 KO (Mesp1-Cre; Gata4 F/F; Gata6 F/+), E8.75, 2 replicates (374 + 409 MB), answers released 2026-10-20. Gata4/Gata6 double loss blocks cardiomyocyte differentiation ([redacted 2026-10-19: B]). **[High]** | Same stage, same WT reference file, same driver as the test: the most informative training knockout, and the only one that can carry the shared knockout-series component. Its Mesp1-lineage biology (likely fewer CMs; P2 model names agree) is Gata-specific. | Split Gata4's delta into shared (non-Mesp1, depth, coverage, reference idiosyncrasy) and Gata-specific (Mesp1 composition/state) parts; transfer the shared part with real cells. Use replicate 1 vs 2 to measure how reproducible each part is. |
| β-catenin is broadly expressed; the organisers call its effect "diffuse" and the hardest. Its two replicates are 300 + 381 MB; never released. **[High]** | Expect many Mesp1 types to shift a little, plus possible composition changes, rather than one programme. Diffuse within-type changes rarely reach abs(lfc) ≥ 0.25 in pseudobulk, so DE genes will again be dominated by composition and nuisance. | Spend most effort on composition + nuisance (de_score, severity), add within-type programmes for de_direction. |
| T1 scRNA E8.5 (16,787 cells) and E9.5 (17,057 cells), whole transcriptome, log-normalised `.X` only (no counts or spliced/unspliced layers in the sample file); E8.5 includes aSHF 4.9%, pSHF 13.1%, JCF 6.7%, CM subtypes ~23%, pericardium 4.4%, endothelium 6.8%, EXEM 5.1%. **[High]** | WT Wnt activity can be measured per cell; E8.75 sits between the two scRNA stages. Scale differs between scRNA (32K genes) and MERFISH (500-gene normalisation). | Wnt-activity score from target genes; per-type regressions; interpolate E8.5/E9.5 to E8.75; calibrate fold changes into MERFISH count space. |
| Scoring: de_score on abs(lfc) ≥ 0.25 genes; de_direction over 500 genes partialling WT level; abs(severity) on DE genes; mmd_u/variogram on the mutant state. **[High]** | See §1.1. | Optimise directly with veckit functions under the scorer's subsampling protocol. |
| P3 data release (2026-10-20): validation answers for every task become training data: Gata4 KO, WT MERFISH E8.5 (heart val_interp) and E10.5 (val_extrap), T1 E10.5, T2 embryo E7.5. **[High]** (site text) | A dense WT MERFISH time series around E8.75 (E8.25, E8.5, E8.75, E9.5, E10.5) becomes available to estimate a "developmental delay" axis and the WT E8.75 sample's idiosyncrasy. | Use them for T3 components D (delay) and R (reference idiosyncrasy), §5. |

---

## 4. Q3: mechanisms — why `shift_transfer` gets the direction wrong, and how composition drives pseudobulk DE

### 4.1 Pseudobulk algebra on sparse MERFISH

Let π_k be the fraction of cells of type k, m_gk the mean log-normalised expression of gene g in type k (zeros included).
Then pb_g = Σ_k π_k m_gk and the knockout change splits exactly into

    dt_g = Σ_k (π'_k − π_k) m_gk          (composition)
         + Σ_k π'_k (m'_gk − m_gk)        (within-type change: biology, depth, batch)

On this MERFISH data a detected count maps to a log value of about 3.6-8.5 (median ~5), so m_gk ≈ 5 × d_gk, where d_gk is the
detection rate of gene g in type k. Rough numbers **[Medium]**:
- A marker detected in 50% of its type's cells has m ≈ 2.5. A 10-percentage-point drop in that type's share moves the
  pseudobulk by ≈ 0.25 = the DE threshold. Composition shifts of ≈ 10 points in a type with strong markers therefore produce
  DE calls for its markers (down if the type is lost, and up for markers of the types whose share rises by closure).
- A within-type detection change of 0.1 in a type making up 10% of cells moves the pseudobulk by ≈ 0.05: invisible to
  de_score, visible to de_direction.
- A depth change shifts detection of nearly every gene at once (Poisson: d = 1 − e^(−λ)). For a broadly expressed gene
  detected in 30% of *all* cells, 36% more counts per cell raises detection to ≈ 38%, while the log value per detected
  count falls from log1p(1e4/75) ≈ 4.9 to log1p(1e4/102) ≈ 4.6; the net pseudobulk change is ≈ +0.3, i.e. above the DE
  threshold. (For a gene confined to a type holding 10% of cells the same arithmetic gives ≈ +0.03.) So a knockout series
  that is ~36% deeper than the WT sample (as in the 150-cell samples) would by itself generate "up" DE calls for the broadly
  expressed genes, most of them correlated with expression level, which raises the de_score null and is partly partialled
  out of de_direction. Whether the full files differ in depth is the first thing to check.
- Section coverage: the Gata4 replicate has ~1.9x the cells of the WT reference. If knockout sections capture more of some
  non-heart tissue, that tissue's markers become "up" DE genes in every knockout at this stage.

Two further sources of shift that a predictor can partly anticipate:
- **Reference idiosyncrasy.** dt = g + (b_KO − b_WT) + noise, where b are embryo/sample effects. The WT E8.75 sample's own
  deviation b_WT from a typical E8.75 embryo enters every knockout's dt with a minus sign. It can be estimated from the WT
  MERFISH time series (E8.25, E8.5, E9.5) once E8.5 is released, and tested on Gata4 (does dt_Gata4 correlate with −b_WT?).
- **Shared knockout-series nuisance** E[b_KO]: if the Gata4 and β-catenin embryos were sectioned, imaged and annotated in
  the same series (same stage, same driver, similar file sizes; unverified), their depth/coverage/batch shift relative to the
  WT atlas sample is shared. **[Low]** until checked.

### 4.2 Why `shift_transfer` (and `gene_ko`) get the direction wrong

`shift_transfer` adds Δ = pb(Mab21l2 KO, E9.5) − pb(WT, E9.5) to every WT E8.75 cell and clips at 0. Five mechanisms, in
rough order of how much they matter:

1. **Clipping asymmetry on 93%-sparse data (mechanical).** Decreases are realised only in cells that express the gene
   (median 4% of the intended change); increases are realised in every cell (100%). The realised decrease ranking is ≈ the WT
   expression ranking (r = 0.99 on the sample), which (a) is the de_score null, so the predicted down set earns nothing, and
   (b) is removed by de_direction's partial correlation. Increases turn zeros into constant positive values in all cells,
   destroying sparsity and co-expression: published variogram skill 15.5 (raw ≈ 0.17 vs floor 0.032) and mmd skill 37.
   **[High]**
2. **Different tissues lost, [redacted 2026-10-19: B].** Each knockout's composition term is signed by *which* types
   change. Mab21l2 loss hits septum transversum/proepicardium/hepatoblasts and myocardium growth at E9.5 (Saito et al.
   2012); Gata4/Gata6 loss in the Mesp1 lineage reduces cardiomyocyte differentiation. [redacted 2026-10-19: B]
   [redacted 2026-10-19: B]
   **[Medium]**
3. **Stage mismatch.** The E9.5 delta includes genes and states that barely exist at E8.75 (more mature chamber
   cardiomyocyte genes, liver, proepicardium) and different baseline composition; added to E8.75 cells it moves the wrong
   genes. The organisers list "ignoring the stage difference" as a pitfall. **[Medium]**
4. **Different reference samples.** The source delta contains b_KO(Mab21l2) − b_WT(E9.5); the target contains
   b_KO(Gata4) − b_WT(E8.75). Unless the knockout series share nuisance, these are independent; in particular the target's
   −b_WT(E8.75) term is absent from the source. This is Systema's point in reverse: a transferred delta carries the source's
   systematic variation and none of the target's. **[Medium]**
5. **The knocked-out gene itself carries almost no signal.** `gene_ko` ≈ `shift_transfer` on the board (42.50 vs 42.48); for
   the test it is a pure no-op because *Ctnnb1* is not on the panel. **[High]**

Public contrast that supports (1) and (4) **[Untrusted]**: mixing real Mab21l2 KO cells (70%) into WT (30%) scored 64.6
with de_direction skill 59.8 (raw ≈ +0.30) and severity 92.3, whereas the additive delta scored ≈ 41-42. The mix keeps
MERFISH sparsity and carries the knockout series' real cell-level properties (and the E8.75 → E9.5 developmental shift);
the additive version keeps neither.

### 4.3 Synthetic demonstration with the official metric code **[Synthetic]**

`scratchpad/scripts/t3_mechanism_sim.py` builds MERFISH-like populations (500 genes, 10 types of which 5 are "Mesp1",
lognormal depth around 75 counts/cell, Poisson counts, `log1p(1e4·c/Σc)`), a WT reference (2,483 cells), a hidden truth
knockout (2 × 2,320 cells), an observed same-stage knockout, and a different-stage knockout with its own WT. By
construction the two same-stage knockouts share a knockout-series nuisance (depth 95 vs 75; extra neural-tube/ectoderm
coverage) but differ in Mesp1 biology. Scored with veckit's own `de_score`, `de_direction`, `severity_slope`, `mmd_unbiased`,
`variogram_score` and `skill`, T3 weights, floor = WT reference, ceiling = truth half vs half. Truth here: 35 up / 5 down DE
genes, expression-level chance 0.35.

| Prediction | T3 | de | dir | sev | mmd | vario |
|---|---|---|---|---|---|---|
| WT cells resampled from the full pool | 66.7 | 0.19 | 0.54 | −1.93 | 0.0051 | 0.0053 |
| shift_transfer of the same-stage KO delta (additive + clip) | 81.8 | 0.65 | 0.54 | −0.26 | 0.0035 | **0.0286** |
| shift_transfer of the different-stage KO delta | 56.5 | 0.00 | 0.62 | **−6.91** | 0.0039 | 0.0195 |
| 70% different-stage KO cells + 30% WT | 60.2 | −0.23 | 0.42 | −1.61 | 0.0067 | 0.0055 |
| Same-stage KO cells as-is ("perturbed mean") | 83.0 | 0.54 | 0.54 | −0.30 | 0.0024 | 0.0023 |
| Systematic only: same-stage KO non-Mesp1 cells + WT-state Mesp1 cells at KO depth | 93.6 | 0.58 | 0.81 | −0.20 | 0.0005 | 0.0008 |
| Systematic + correct Mesp1 composition (no programme) | 98.1 | 0.62 | 0.89 | −0.11 | 0.0000 | 0.0009 |
| Composition only (no depth/coverage) | 73.3 | 0.27 | 0.67 | −1.25 | 0.0036 | 0.0045 |
| Depth only | 62.9 | 0.15 | 0.59 | −6.91 | 0.0026 | 0.0027 |
| Systematic + **wrong-sign** Mesp1 composition | 54.8 | 0.12 | **0.18** | **−6.91** | 0.0050 | 0.0044 |
| Oracle (independent draw of the truth model) | 99.6 | 0.65 | 0.89 | −0.03 | 0.0003 | 0.0007 |

What it shows, and what it does not:
- Additive shifts damage the variogram even when the transferred delta is informative (0.029 vs floor 0.0066), matching the
  direction of the real board (≈ 0.17 vs 0.032). Real-cell transplant keeps the distribution metrics.
- A wrong-sign composition hypothesis collapses de_direction and sends severity to its floor even when the systematic part
  is right. This is the risk to hedge with the second submission.
- It does **not** show that the real knockout series share nuisance: that was built in. The simulation is also optimistic
  for cross-stage transfer and for "WT resampled" (its de_direction of 0.54 comes from noise shared with the scorer's
  reference subsample, which is small in real data: a public WT-copy upload scored de_direction skill 50.6, about +0.01 raw).
  Re-run the same harness on the real files on day 1.

---

## 5. Q4: concrete pipelines

Notation for components of the predicted change: **N** = shared knockout-series/systematic component; **R** = WT-reference
idiosyncrasy correction; **D** = developmental-delay axis; **C** = Mesp1-lineage composition change; **B** = within-type
Wnt-responsive programme. All pipelines generate a cell population (counts → `log1p(1e4·c/Σc)`, panel order from
`T3__gata4.genes.txt`-style file for the test board, `obsm["spatial_3D"]` copied from source cells), never additive log
shifts.

### P1. Systematic-component transfer with real cells ("same-stage perturbed mean, done right")

Idea: carry over from Gata4 everything that is not Gata-specific, using real Gata4 cells where possible.

Steps:
1. **Harmonise labels** across WT E8.25/E8.5/E8.75/E9.5/E10.5 MERFISH, Mab21l2 KO and both Gata4 replicates: train a
   classifier on the WT heart vocabulary (logistic regression on log-normalised panel, or scANVI on counts with `sample` as
   batch; `scvi-tools` 1.5.1) and assign KO cells; keep a crosswalk to the KO files' own labels. Make obs names unique (a
   public pitfalls list reports duplicates in released T2/T3 files **[Untrusted]**).
2. **Describe each sample**: counts/cell and genes/cell distributions; harmonised composition; per-type detection profiles;
   `obs` batch fields (`perturb_batch`, `slide`, `sample_id`); and whether `obs["condition"]` reveals littermate WT cells
   inside knockout files (the data page lists this column as "tbc"). Littermate WT cells, if present, give a batch-matched
   reference that separates batch from genotype directly.
3. **Estimate N** from compartments that should not respond: non-Mesp1 types (foregut/hepatic endoderm, surface ectoderm,
   neural tube, NCC, blood) in Gata4 rep1 and rep2 versus WT E8.75. Check concordance rep1 vs rep2 and against Mab21l2 vs
   WT E9.5 (shared series effects across stages?).
4. **Build the N-carrier population**: Gata4 non-Mesp1 cells in Gata4 proportions (they carry depth, coverage and batch),
   plus WT E8.75 Mesp1 cells re-expressed "as if" in the knockout series: either scVI counterfactual decoding with
   `transform_batch` set to the knockout series (supported in `SCVI.get_normalized_expression`) followed by count sampling,
   or empirical depth matching (binomial thinning if the KO series is shallower; Poisson count augmentation from a
   kNN-smoothed profile if deeper).
5. Optional **R**: residual of WT E8.75 from a smooth WT time trajectory (E8.25, E8.5, E9.5 per harmonised type and
   composition); add −R only if dt_Gata4 (and dt_Mab21l2 with the E9.5 analogue) correlate with −R.

Expected value: if the day-1 decomposition shows N is a large share of Gata4's dt, P1 alone should beat the cross-stage
real-cell mixes reported on the P2 board (~64-66 on Gata4) because it uses same-stage, same-reference cells and avoids
importing another knockout's Mesp1 biology; it mainly secures severity and part of de_score. This is a qualitative
projection from the P2 board and the simulation, not a measurement, and I give no number. **[Low]**
Risk: Gata-specific non-autonomous changes leak in (e.g., endoderm responding to the heart defect); mitigate by checking
replicate concordance and down-weighting non-Mesp1 types with large replicate disagreement.
Compute: CPU-minutes to hours; scVI training on ~150K cells × 500 genes is minutes on one H100.
Legality: released data only. **Allowed.**

### P2. WT-derived Wnt-programme in-silico knockout (GRN propagation done on pathway output, not on *Ctnnb1* mRNA)

Idea: learn, in WT cells, how panel genes co-vary with Wnt/β-catenin *output* within each Mesp1 cell type, then set that
output to its "off" level in Mesp1 cells.

Steps:
1. T1 scRNA E8.5 + E9.5; keep Mesp1-like types (aSHF, pSHF, JCF, CM subtypes, pericardium/proepicardium, ST,
   endocardium/endothelium flagged as partial, pharyngeal mesoderm, EXEM) and non-Mesp1 controls.
2. **Wnt activity W_i per cell** from direct targets that are *not* on the panel ([redacted 2026-10-19: B]
   [redacted 2026-10-19: B]; all present in T1, checked), e.g. AUCell or a z-score; cross-check with pySCENIC
   TCF/LEF regulon activity. [redacted 2026-10-19: B] Do **not** use *Ctnnb1* mRNA as the
   activity proxy.
3. **Per-type responsiveness**: for each panel gene g and type k, regress expression on W_i with covariates (library size,
   stage, cell-cycle score, sub-cluster), with empirical-Bayes shrinkage of β_gk; also test W against proliferation genes on
   the panel (Ccna2, Cdk1, Cenpf, Top2a, Prc1, Diaph3) to get a "Wnt-linked proliferation" index per type.
4. **In-silico KO**: Δlog-mean_gk = −β_gk·(W̄_k − W_off), with W_off a low quantile of W in that lineage; optional
   CellOracle-style propagation (ridge GRN fitted on T1 with a mouse base GRN, TCF/LEF as entry points, `n_propagation=3`)
   for second-order targets; scTenifoldKnk only as a ranking cross-check.
5. **Map to MERFISH**: harmonise T1 types to MERFISH types via a classifier on the shared 500 genes; interpolate E8.5/E9.5 to
   E8.75; convert fold changes into detection changes per gene and type (d' = 1 − exp(−λ·FC)) and apply in count space to the
   Mesp1 cells of the P1 carrier; re-normalise.
6. **Composition hook for P3**: per-type Wnt-dependence index (mean W, Wnt-linked proliferation).

Validation unique to P2: run the same machinery for **Gata4 + half Gata6** (CellOracle can perturb TFs directly) and for
Mab21l2 (co-variation programme) and score against the released truths' *within-type* changes (computed per harmonised type,
not pseudobulk), to learn whether WT-learned programmes predict knockout programmes in this system at all.
Expected value: mostly de_direction and tie-breaking within the de_score selection; small effect on DE genes. **[Medium]**
Risks: correlation ≠ causation (W also tracks progenitor vs differentiated state; cell-type covariates reduce but do not remove
this); scRNA-to-MERFISH calibration; E8.75 interpolation.
Compute: CPU; pySCENIC/GRNBoost2 on ~34K cells is a few hours on many cores (DGX).
Legality: released T1 data; public motif databases and base GRNs built from WT chromatin. **Allowed** (disclose).

### P3. Cell-type-resolved composition model with explicit hedging

Idea: the composition term C drives most DE genes; predict it explicitly and hedge its sign.

Steps:
1. Quantify composition changes in the two observed knockouts per harmonised Mesp1 type (scCODA for credible changes,
   `pertpy` 1.4.0 `Sccoda`/`Milo` for neighbourhood-level differential abundance; Büttner et al. 2021; Dann et al. 2022).
   Derive empirical-Bayes hyperparameters: typical |log fold change| of type shares in a Mesp1-Cre knockout at E8.75, and
   how much "loss" propagates to derived types.
2. Prior mean for β-catenin: log(π'_k/π_k) = −γ·WntDep_k (from P2), propagated to lineage-derived types with a lag
   weight; σ from step 1.
3. Scenario set (all from WT-derived reasoning, no knockout literature): S-prolif: Wnt-high progenitor types reduced and
   their derivatives modestly reduced; S-diff: progenitors reduced, early cardiomyocytes increased (premature
   differentiation); S-mild: small changes only. Choose γ grids per scenario.
4. Generate populations by resampling cells per type with new shares (no duplicates where avoidable) on top of the P1
   carrier with the P2 programme.
5. Choose the two test submissions by expected best-of-two score under scenario probabilities you set from day-1 analyses
   (e.g., whether Wnt activity is concentrated in progenitors, how Gata4's composition change looks).

Expected value: high if the sign is right; the wrong sign is the main way to fall to the floor on severity (simulation).
**[Medium]**
Legality: priors from WT data and released knockouts: **allowed**. Priors from Ctnnb1 cKO papers: **risk** (§7).

### P4. Empirical-Bayes / stacking combination fitted by leave-one-knockout-out (LOKO)

Model: predicted population = generator(N, R, D, C, B; weights w). For each training knockout build its *own* analogues:
N from the *other* knockout (or from its non-Mesp1 compartment), D from the WT time series, B from its own gene's WT
co-variation programme (Gata4/Gata6 via CellOracle; Mab21l2 via co-expression), C from its gene's WT-derived dependence.
Fit w by maximising the official T3 score (veckit functions, scorer subsampling protocol, several seeds), with strong priors
(two folds only): w ~ N(1, 0.5²) on N, half-normal on D/R/B/C. Report the posterior spread and prefer flat, robust optima.
Magnitude: calibrate the overall response so that the slope β on DE genes would be ≈ 1 for the LOKO folds (severity), then
check mmd/variogram did not degrade.
Legality: method only. **Allowed.**

### P5. LLM-elicited gene-level priors (optional; legal clearance first)

Protocol if cleared: for each Mesp1 type × 500 genes, ask an LLM for P(up), P(down) and a magnitude class under
"Mesp1-Cre; Ctnnb1 F/F at E8.75", with ~5 samples per prompt and majority aggregation; calibrate its weight by eliciting the
same for Gata4 (before looking at the Gata4 truth) and for Mab21l2, then scoring those elicitations against the truths.
Cost is negligible (batch 500 genes per prompt).
Why it is risky: the rules exclude data from the held-out genotype or "another allele of the same gene at a comparable
stage", "not through a model pre-trained on it". General LLMs have read the Ctnnb1 cardiac loss-of-function papers,
including one with exactly this genotype (§7). Treat as **prohibited until the organisers say otherwise**; if allowed,
disclose model, prompts and outputs.

### 5.6 Recommended assembly and the two test submissions

- **Submission 1 (robust):** P1 carrier (N, plus R if supported) + P2 programme (B) + P3 central scenario (C), weights from
  P4, magnitude calibrated for severity, ~maximum allowed cells.
- **Submission 2 (hedge):** same N, but the alternative Mesp1 composition scenario (or a stronger/weaker B), chosen to
  maximise expected best-of-two. Using submission 1's returned score to choose between pre-registered candidates is
  allowed ("selection"); solving for target properties from scores is not.
- Timing: submit 1 around 2026-11-20 to 11-23 and 2 by 2026-11-30, leaving buffer before 2026-12-02.

---

## 6. Validation protocol

### 6.1 Day-1 diagnostics (2026-10-20 to 10-22), all CPU, on the DGX

1. Depth, genes/cell, zero fraction per sample (WT E8.25/E8.5/E8.75/E9.5/E10.5, Mab21l2, Gata4 rep1/rep2).
2. Harmonised composition per sample; Gata4 rep1 vs rep2 vs WT for non-Mesp1 and Mesp1 compartments.
3. `obs["condition"]`/batch fields; littermate WT present or not.
4. **Oracle decomposition of Gata4** with veckit (scorer protocol, 5-10 subsample seeds), each candidate scored against both
   Gata4 replicates:
   - WT resampled to Gata4 composition (composition only);
   - WT at Gata4 depth (depth only);
   - Gata4 non-Mesp1 + WT Mesp1 (systematic only, P1);
   - WT non-Mesp1 + Gata4 Mesp1 (biology only);
   - Gata4 replicate A predicting replicate B (realistic upper bound across embryos).
   The same on Mab21l2 vs WT E9.5. This tells you where T3 points come from in *this* data.
5. Correlate dt_Gata4 with D (WT E8.5 − E8.75), with −R, and with dt_Mab21l2 (stage-matched as far as possible).

### 6.2 Leave-one-knockout-out folds

| Fold | Train | Predict | Score with | Notes |
|---|---|---|---|---|
| A (now, P2 board until 10-20; offline after) | Mab21l2 + WT | Gata4 rep1, rep2 | veckit vs WT E8.75 | Cross-stage: pessimistic for N transfer |
| B (offline) | Gata4 + WT | Mab21l2 | veckit vs WT E9.5 | Cross-stage the other way |
| C (offline) | Gata4 rep1 (+ Mab21l2) | Gata4 rep2 | veckit | Same-stage proxy for "how much transfers between knockout embryos" (upper bound for N) |
| Test | Mab21l2 + Gata4 + WT | β-catenin | official board, 2 submissions | — |

Report every metric with intervals over subsample seeds; treat mmd/variogram differences below ~3-6 skill points as noise
(public report **[Untrusted]**); de/dir/severity are deterministic given the files.

### 6.3 Decision rules

- Keep a component only if it improves the LOKO T3 score in both directions (A and B) or in fold C, and never if it flips
  severity's sign in any fold.
- Prefer components whose sign agrees across Gata4 replicates.
- Freeze the pipeline by ~2026-11-16; after that only choose among pre-registered variants.

---

## 7. Legality matrix and questions for the organisers

Rules §10 (verbatim points): no measured data from a held-out genotype "by any route: not directly, not through a model
pre-trained on it, and not through a public dataset that contains it"; "Data from the held-out knockouts, from another allele
of the same gene at a comparable stage, or from a perturbation that phenocopies them, is treated as the held-out condition";
general resources spanning held-out data may be used only after removing the in-window part and describing the removal;
"If a source sits near the boundary and you are unsure, ask before submitting."

| Source / technique | Status | Note |
|---|---|---|
| All released competition data (WT MERFISH stages, T1 scRNA, Mab21l2 KO, Gata4 KO after 10-20) | Allowed | Intended training material |
| Generic motif databases, base GRNs from WT chromatin atlases (CellOracle), pySCENIC cisTarget | Allowed (disclose) | WT information |
| CollecTRI / PROGENy / NicheNet / STRING priors | Allowed after provenance filtering (disclose) | Drop edges whose evidence comes from Ctnnb1 or Wnt-pathway loss-of-function in embryonic cardiac mesoderm at E7.5-E10.5 (e.g., [redacted 2026-10-19: B]); avoid STRING's text-mining channel |
| Foundation models pre-trained on large mouse corpora (Mouse-Geneformer 21M cells; Navigo checkpoints) | Gray | Must confirm the corpus contains no Ctnnb1 loss-of-function data at comparable stages; low expected benefit, so skip |
| Published Ctnnb1 cKO/GOF in cardiac progenitors: [redacted 2026-10-19: B] | **Treat as prohibited** unless cleared | Their abstracts name panel genes and phenotypes; using them, even qualitatively, is plausibly "data from the held-out genotype / another allele at a comparable stage" |
| Wnt-pathway perturbations that could phenocopy (e.g., mesodermal Wnt secretion/receptor loss, Dkk1 overexpression in embryos) | Treat as prohibited unless cleared | "perturbation that phenocopies them" |
| LLM-elicited gene priors | **Ask first** | "not through a model pre-trained on it" |
| Other teams' public code | Avoid | Rules prohibit code sharing across teams; reading public write-ups for ideas is fine |
| Using returned scores to *compute* predictions (inverting the skill map, differencing scores) | Prohibited | Selection among pre-registered candidates is allowed |

Questions to send now (organisers answer in advance, and an advance answer protects you):
1. May we use qualitative findings from published Ctnnb1 loss-of-function studies in cardiac mesoderm (any Cre, E7.5-E10.5),
   including [redacted 2026-10-19: B]? Our reading is no.
2. May we use gene-level priors elicited from general-purpose LLMs, given they were pre-trained on that literature?
3. May we use CollecTRI/PROGENy/NicheNet/STRING after removing interactions supported only by such studies?
4. Is it acceptable to use the released Gata4 knockout cells directly as part of the predicted β-catenin population
   (shared non-genotype component)? (We read this as allowed: it is released training data.)
5. Do the knockout files contain littermate WT cells (`obs["condition"]`), and may they be used as a batch-matched
   reference?

---

## 8. Six-week plan (2-4 H100s available; most work is CPU)

| Dates | Work | Output |
|---|---|---|
| Oct 18-19 | Register; set up uv envs on the DGX: a Python 3.12 env (scanpy 1.12.4, anndata 0.13.4, scvi-tools 1.5.1, pertpy 1.4.0, decoupler 2.2.0 [needs ≥3.11], veckit from GitHub HEAD) and a separate older-Python env for celloracle 0.20.0 / pyscenic 0.12.1 (2022 release; compatibility with 3.12 not verified); send §7 questions; offline harness reproducing the scorer protocol | Harness + emails |
| Oct 20-26 | Download (T1 E8.5/E9.5 571+590 MB; WT MERFISH 84/164/[redacted 2026-10-19: A]/225/[redacted 2026-10-19: A] MB; Mab21l2 458 MB; Gata4 374+409 MB); label harmonisation; §6.1 diagnostics and oracle decomposition | Decision memo: how large are N, C, B |
| Oct 27-Nov 2 | P1 carrier (scVI counterfactual or empirical depth matching); P2 Wnt score + per-type responsiveness; scRNA→MERFISH calibration | First LOKO numbers |
| Nov 3-9 | P3 composition scenarios; P2 validation on Gata4/Mab21l2 programmes; P4 weight fitting | Ranked candidate list |
| Nov 10-16 | Robustness (seeds, replicate swaps), magnitude calibration for severity, legality-dependent switches (P5 or literature only if cleared); freeze | Pre-registered S1/S2 |
| Nov 17-23 | Test submission 1; format checks unlimited | Score 1 |
| Nov 24-30 | Choose submission 2 among pre-registered variants; submit; start method report (due within 14 days of winning) | Score 2 |
| Dec 1-2 | Buffer | — |

---

## 9. Open questions and things I could not verify

- Whether the β-catenin and Gata4 knockout embryos share a processing series (depth, coverage, batch). Everything in P1 that
  goes beyond "use the same-stage, same-reference knockout" depends on it. **[Unverified]**
- Whether the full knockout files are deeper per cell than WT (seen only in 150-cell samples). **[Unverified]**
- Whether the matched WT E8.75 (the T2 heart-series file) is a littermate or Cre-carrying control. **[Unverified]**
- Whether knockout files contain WT cells (`obs["condition"]` marked "tbc"). **[Unverified]**
- The test board's cell limits and which β-catenin replicate the scorer uses (fixed by the scorer, not published).
- Allele and Cre status of the Mab21l2 knockout in this dataset. **[Unverified]**
- Mouse support in LINGER; exact `organism` values accepted by decoupler 2.2 (use `op.show_organisms()`). **[Unverified]**
- The competitor-reported mmd/variogram seed noise and P2 experiment scores. **[Untrusted]**
- Kernfeld et al. is cited as the 2023 preprint; I did not find a journal version.

---

## 10. References

Bibliographic metadata for journal articles was retrieved from PubMed (DOIs below); preprints from bioRxiv/arXiv.

Challenge and tooling
- Virtual Embryo Challenge site (tasks, data, rules, baselines, evaluation): https://virtualembryo.ai/challenge ;
  T3 task page https://virtualembryo.ai/challenge/tasks/perturbation ; baselines https://virtualembryo.ai/challenge/baselines ;
  rules https://virtualembryo.ai/challenge/rules
- Public leaderboard data (T3, fetched 2026-10-18): https://api.virtualembryo.ai/challenge/leaderboard?task=T3
- veckit official scorer (MIT): https://github.com/aristoteleo/veckit ; PyPI https://pypi.org/project/veckit/ (0.1.2)
- Cao N, Lu Y, Qiu X. Towards predictive virtual embryos with genomics and AI. Nat Methods (2026).
  https://doi.org/10.1038/s41592-026-03055-4
- Navigo: https://github.com/aristoteleo/Navigo-release
- CellOracle source (`simulate_shift`, base GRNs): https://github.com/morris-lab/CellOracle
- dynamo in-silico perturbation docs: https://dynamo-release.readthedocs.io/en/latest/introduction/index_silico.html
- scvi-tools SCVI API (`transform_batch`): https://docs.scvi-tools.org/en/stable/api/reference/scvi.model.SCVI.html
- decoupler OmniPath resources: https://decoupler.scverse.org/en/latest/api/op.html
- pertpy tools: https://pertpy.scverse.org/en/latest/api/tools_index.html
- OmniPath CollecTRI query used here:
  [redacted 2026-10-19: B]
- Competitor methods journal (untrusted, ideas only): https://github.com/wqty123/virtual-embryo-challenge-methods

Perturbation-prediction methods
- Roohani Y, Huang K, Leskovec J. GEARS. Nat Biotechnol 42:927-935 (2024). https://doi.org/10.1038/s41587-023-01905-6
- Lotfollahi M et al. CPA. Mol Syst Biol 19:e11517 (2023). https://doi.org/10.15252/msb.202211517
- Lotfollahi M, Wolf FA, Theis FJ. scGen. Nat Methods 16:715-721 (2019). https://doi.org/10.1038/s41592-019-0494-8
- Bunne C et al. CellOT. Nat Methods 20:1759-1768 (2023). https://doi.org/10.1038/s41592-023-01969-x
- Piran Z et al. biolord. Nat Biotechnol 42:1678-1683 (2024). https://doi.org/10.1038/s41587-023-02079-x
- Yu H et al. PerturbNet. Mol Syst Biol 21:960-982 (2025). https://doi.org/10.1038/s44320-025-00131-3
- Adduri A et al. State (bioRxiv 2025). https://doi.org/10.1101/2025.06.26.661135
- Wenkel F et al. TxPert (arXiv 2025). https://arxiv.org/abs/2505.14919
- He S et al. Squidiff. Nat Methods 23:65-77 (2026). https://doi.org/10.1038/s41592-025-02877-y
- Cui H et al. scGPT. Nat Methods 21:1470-1480 (2024). https://doi.org/10.1038/s41592-024-02201-0
- Theodoris CV et al. Geneformer. Nature 618:616-624 (2023). https://doi.org/10.1038/s41586-023-06139-9
- Ito K et al. Mouse-Geneformer. PLoS Genet 21:e1011420 (2025). https://doi.org/10.1371/journal.pgen.1011420
- Kamimoto K et al. CellOracle. Nature 614:742-751 (2023). https://doi.org/10.1038/s41586-022-05688-9
- Qiu X et al. dynamo. Cell 185:690-711 (2022). https://doi.org/10.1016/j.cell.2021.12.045
- Osorio D et al. scTenifoldKnk. Patterns 3:100434 (2022). https://doi.org/10.1016/j.patter.2022.100434
- Bravo González-Blas C et al. SCENIC+. Nat Methods 20:1355-1367 (2023). https://doi.org/10.1038/s41592-023-01938-4
- Yuan Q, Duren Z. LINGER. Nat Biotechnol (2024). https://doi.org/10.1038/s41587-024-02182-7
- Chen Y, Zou J. GenePT (bioRxiv). https://doi.org/10.1101/2023.10.16.562533
- Märtens K et al. LangPert (ICLR 2025 MLGenX workshop). https://mlanthology.org/iclrw/2025/martens2025iclrw-langpert
- Müller-Dott S et al. CollecTRI. Nucleic Acids Res 51:10934-10949 (2023). https://doi.org/10.1093/nar/gkad841
- Schubert M et al. PROGENy. Nat Commun 9:20 (2018). https://doi.org/10.1038/s41467-017-02391-6
- Browaeys R et al. NicheNet. Nat Methods 17:159-162 (2020). https://doi.org/10.1038/s41592-019-0667-5
- Büttner M et al. scCODA. Nat Commun 12:6876 (2021). https://doi.org/10.1038/s41467-021-27150-6
- Dann E et al. Milo. Nat Biotechnol 40:245-253 (2022). https://doi.org/10.1038/s41587-021-01033-z

Benchmarks and challenges
- Ahlmann-Eltze C, Huber W, Anders S. Nat Methods 22:1657-1661 (2025). https://doi.org/10.1038/s41592-025-02772-6
- Csendes G et al. BMC Genomics 26:393 (2025). https://doi.org/10.1186/s12864-025-11600-2
- Kernfeld EM et al. (bioRxiv 2023). https://doi.org/10.1101/2023.07.28.551039
- Viñas Torné R et al. Systema. Nat Biotechnol 44:1050-1059 (2025). https://doi.org/10.1038/s41587-025-02777-8
- Roohani YH et al. Virtual Cell Challenge. Cell 188:3370-3374 (2025). https://doi.org/10.1016/j.cell.2025.06.008
- Arc Institute. Virtual Cell Challenge 2025 wrap-up. https://arcinstitute.org/news/virtual-cell-challenge-2025-wrap-up
- Participant post-mortem (anecdotal): https://yraj.in/ai/virtual-cell/perturbation-modeling/single-cell/2025/12/17/virtual-cell-challenge.html
- Huang X et al. Single-cell, whole-embryo phenotyping of mammalian developmental disorders. Nature 623:772-781 (2023).
  https://doi.org/10.1038/s41586-023-06548-w (22 mutants at E13.5; shared changes among mutants; WT strain differences;
  evidence only, not proposed as training data)

Biology used for constraints
- Saga Y et al. Mesp1 lineage. Development 126:3437-3447 (1999). https://doi.org/10.1242/dev.126.15.3437
- RIKEN BRC RBRC01145, Mesp1<tm2(cre)Ysa>: https://brc.riken.jp/mus/RBRC01145
- Saito Y, Kojima T, Takahashi N. Mab21l2. PLoS One 7:e32991 (2012). https://doi.org/10.1371/journal.pone.0032991
- [redacted 2026-10-19: B]
- [redacted 2026-10-19: B]

Held-out-genotype literature (listed so the team can avoid it until cleared)
- [redacted 2026-10-19: B]
- [redacted 2026-10-19: B]
- [redacted 2026-10-19: B]
- [redacted 2026-10-19: B]
- [redacted 2026-10-19: B]

Scripts written for this report (scratchpad, tiny local runs)
- `scratchpad/scripts/inspect_h5ad.py`, `libsize_check.py`, `counts_check.py`, `gene_det.py`, `clip_attenuation.py`,
  `t1_inspect.py`, `lb_t3_table.py`, `collectri_panel.py`, `t3_mechanism_sim.py`
