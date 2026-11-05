---
task_id: h
role: Measurement scientist for MERFISH and single-cell RNA-seq (technical noise, batch and imaging effects)
status: complete
sources_found: 15
---
## Sources
[1] Multiplexed imaging of high-density libraries of RNAs with MERFISH and expansion microscopy (Wang, Moffitt, Zhuang; Sci Rep; read in full through the Europe PMC full-text mirror of this record) | https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5859009/ | Source-Type: academic | Accessibility: public | As Of: 2018-03 | Authority: 8/10
[2] Gene count normalization in single-cell imaging-based spatially resolved transcriptomics (Atta et al.; Genome Biol; DOI 10.1186/s13059-024-03303-w; read in full) | https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11167774/ | Source-Type: academic | Accessibility: public | As Of: 2024-06 | Authority: 8/10
[3] Library size confounds biology in spatial transcriptomics data (Bhuva et al.; Genome Biol; DOI 10.1186/s13059-024-03241-7; read in full through Europe PMC, PMC11025268) | https://GenomeBiology.biomedcentral.com/articles/10.1186/s13059-024-03241-7 | Source-Type: academic | Accessibility: public | As Of: 2024-04 | Authority: 8/10
[4] The effect of background noise and its removal on the analysis of single-cell expression data (Janssen et al.; Genome Biol; DOI 10.1186/s13059-023-02978-x; read through Europe PMC, PMC10278251, up to the clustering section) | https://consensus.app/papers/details/bccaab7df4745253935ddac54ac111e9/ | Source-Type: academic | Accessibility: public | As Of: 2023-06 | Authority: 8/10
[5] Impact and correction of segmentation errors in spatial transcriptomics (Mitchel et al.; Nat Genet; DOI 10.1038/s41588-025-02497-4; abstract and PubMed metadata only) | https://consensus.app/papers/details/3a0296b31207580ebac6905b7705b50a/ | Source-Type: academic | Accessibility: semi-public | As Of: 2026-01 | Authority: 8/10
[6] Cell segmentation in imaging-based spatial transcriptomics (Baysor; Petukhov et al.; Nat Biotechnol; DOI 10.1038/s41587-021-01044-w; abstract only) | https://consensus.app/papers/details/a67a407890eb5dc085c3988fc67f15c7/ | Source-Type: academic | Accessibility: semi-public | As Of: 2021-10 | Authority: 8/10
[7] Spatially resolved, highly multiplexed RNA profiling in single cells (MERFISH; Chen et al.; Science; full text blocked by an anti-bot page, search summary only) | https://pmc.ncbi.nlm.nih.gov/articles/PMC4662681/ | Source-Type: academic | Accessibility: public | As Of: 2015-04 | Authority: 8/10
[8] RNA imaging with MERFISH: design of oligonucleotide probes (protocol; search summary only) | https://www.protocols.io/view/rna-imaging-with-merfish-design-of-oligonucleotide-menc3de.html | Source-Type: community | Accessibility: public | As Of: 2018-03 | Authority: 6/10
[9] SoupX removes ambient RNA contamination from droplet-based single-cell RNA sequencing data (Young, Behjati; GigaScience; DOI 10.1093/gigascience/giaa151; abstract) | https://consensus.app/papers/details/5c45823b193e5917bebaef289bf2b16c/ | Source-Type: academic | Accessibility: public | As Of: 2020-12 | Authority: 7/10
[10] Decontamination of ambient RNA in single-cell RNA-seq with DecontX (Yang et al.; Genome Biol; abstract) | https://consensus.app/papers/details/8194c3b5a7695f9c8961df6e6e3186d6/ | Source-Type: academic | Accessibility: public | As Of: 2020-03 | Authority: 7/10
[11] Dissection of artifactual and confounding glial signatures by single-cell sequencing of mouse and human brain (Marsh et al.; Nat Neurosci; DOI 10.1038/s41593-022-01022-8; abstract) | https://consensus.app/papers/details/d1b61a3c7de75c1cba23f49bd3e87272/ | Source-Type: academic | Accessibility: semi-public | As Of: 2022-03 | Authority: 8/10
[12] Dissociation protocols used for sarcoma tissues bias the transcriptome observed in single-cell and single-nucleus RNA sequencing (Truong et al.; BMC Cancer; DOI 10.1186/s12885-023-10977-1; abstract) | https://consensus.app/papers/details/da815397966f513e9288a7dcdf5c75d5/ | Source-Type: academic | Accessibility: public | As Of: 2023-05 | Authority: 6/10
[13] Single-nucleus and single-cell transcriptomes compared in matched cortical cell types (Bakken et al.; PLoS One; DOI 10.1371/journal.pone.0209648; abstract) | https://consensus.app/papers/details/0b6e97913daf527abd32ee04160c978a/ | Source-Type: academic | Accessibility: public | As Of: 2018-12 | Authority: 7/10
[14] SpaNorm: spatially-aware normalisation for spatial transcriptomics data (bioRxiv preprint; search summary only) | https://www.biorxiv.org/content/10.1101/2024.05.31.596908.full.pdf | Source-Type: academic | Accessibility: public | As Of: 2024-05 | Authority: 6/10
[15] bioRxiv 2021.11.01.466826, MERFISH of human and mouse cortex with 250- and 4,000-gene panels (community mirror of the preprint; search summary only) | https://storage.prod.researchhub.com/uploads/papers/2025/11/22/2021.11.01.466826v1.full.pdf | Source-Type: academic | Accessibility: public | As Of: 2021-11 | Authority: 5/10

## Findings
1. In unexpanded samples at high total RNA density, MERFISH decoded only about a fifth of the copies smFISH counted (mean 21%, median 16%, 12 genes) against near-complete detection after 2.3-fold linear expansion, and simulation put a ~50% detection loss at about 3-4 molecules per square micrometre per 1-µm optical section [1].
2. Per-gene copy numbers from two independent MERFISH preparations correlated at r = 0.94 (expanded against unexpanded r = 0.88), and a 4,000-gene panel detected about 57% of what a 250-gene panel detected, so sizeable per-gene offsets between runs coexist with high replicate correlations [1][15].
3. Count-based (library-size) normalization of a targeted panel skewed toward one region's or cell type's genes deflates normalized expression where the panel is enriched and inflates it elsewhere, giving up to 13% DE errors and up to 19% fold-change sign switches on MERFISH, while cell-volume normalization avoided panel effects but depends on segmentation and partial z-capture [2].
4. Total detections per cell track tissue structure on four platforms even after adjusting for cell count, and region-specific library-size effects reach from about a quarter to nearly all genes depending on the dataset, so per-cell depth is part biology and part technique [3][14].
5. Background counts, mostly ambient RNA, averaged 3-35% of UMIs per cell, varied strongly between replicates, were several-fold higher in single-nucleus than in whole-cell libraries, and in nuclei did not mirror the sample's cell-type mix; soup composition and magnitude are experiment-specific [4][9][10].
6. Segmentation errors (molecules assigned to the wrong cell) confound most downstream cell-state analyses in imaging-based spatial transcriptomics, differential expression included, and often dominate the result [5].
7. The choice of segmentation method alone can nearly double the number of cells called from the same MERFISH measurements [6].
8. MERFISH reads each barcode over many rounds with an error-robust code (minimum Hamming distance 4 corrects one bit), so a gene's calling rate depends on the run's per-bit error rates, and unassigned "blank" barcodes estimate the false-positive floor (5-10% of codewords are left blank by design) [7][8].
9. Warm enzymatic dissociation induces an ex vivo stress signature in whole-cell preparations that can confound analyses, while single-nucleus preparations show a gene-length bias instead [11][12].
10. The nuclear share of a cell's mRNA ranged from about 20% to over 50% across neuron types, so whole-cell and nuclear abundances of the same gene differ with cell size and transcript localization [13].

## First-principles analysis

### 1. Observation model: where each technical component enters
For sample s (one embryo stack, imaged on run or slide r(s)), cell i of type t, panel gene g:

    E[c_sig] = kappa_r * V_si * eta_rg * psi_r(D_si) * rho_tg  +  sum_j w_ij,s * c_sjg  +  beta_r * A_g

- Per-gene detection efficiency: eta_rg = pi_g * prod_{b in B(g)} q_rb * sigma_g(SNR_r). pi_g is the probe set (probe count, which grows with usable transcript length until capped; GC and melting temperature), fixed within a series by the probe library. q_rb is the run's readout efficiency for bit b; with one-bit correction a weight-4 codeword is called with probability about q^4 + 4 q^3 (1 - q) when the bits behave alike (assumption: a weight-4, distance-4 code, the classic design [7][8]). sigma_g is an SNR threshold: genes with few probes (dim spots) drop out first when a run's signal-to-background falls.
- Optical crowding: psi_r(D) falls with local spot density D [1]; it is gene-specific only through abundance and shared bits.
- Cell volume, section depth and z-sampling: V_si is the captured volume (cell size times the fraction inside the section; partial capture at stack ends [2]); the segmentation granularity sets how much volume each object gets [6].
- Segmentation spillover: w_ij,s, molecules of neighbour j assigned to i, in xy and by z-projection within a section [5].
- Ambient signal: beta_r * A_g, off-target binding, autofluorescence-driven misdecodes and extracellular RNA; blanks measure its floor [7][8].
- Imaging run and batch: every r-indexed term (kappa, q, sigma, psi, beta, segmentation settings), plus the release pipeline (file stream, filters, synthetic records).
Assumption: serial sections with one 2D segmentation per section; the files store z as a section index times a per-file pitch (`docs/DATA_NOTES.md`), which supports this.

### 2. How these terms reach the scored quantities
The released X is log1p(1e4 c / N) with N = sum_g c. Proportions p = c/N cancel kappa and V in expectation, but sparsity does not. With lambda = E[c] much below 1 (the regime of ~90% zeros on the 500-gene panel), E[x] is about lambda * ln(1 + 1e4/N) = p * h(N), with h(N) = N ln(1 + 1e4/N). So the scorer's pseudobulk is PB_g = p_g * h(N), and:
- (a) Depth. A depth ratio r between two samples moves every sparse gene by Delta_g = [h(rN)/h(N) - 1] * PB_g. At N = 130 (WT E8.75's median), r = 0.55 gives about -0.37 PB and r = 1.27 about +0.20 PB. A synthetic count simulation (my scratchpad `synth_h/depth_pb.py`, arbitrary log-normal expression, no data) reproduced -0.39 to -0.33 PB up to PB 1.5, rolling off to about half that for the densest genes. Genes with PB above about 0.7 then cross the scorer's |lfc| >= 0.25 on depth alone.
- (b) A per-gene efficiency ratio phi_g moves PB_g by about PB_g * ln(phi_g) in the sparse regime, after the panel renormalization, which is what the T3 factor models.
- (c) Contamination adds (neighbour or ambient level) / (own level), so it is type-dependent and largest where a gene is foreign to the cell. A heart-skewed 500-gene panel also makes CP10k couple genes within a cell: losing cardiac transcripts inflates every other panel gene in those cells [2].

How each T3 metric responds (formulas in `common/core_metrics.py` at cc3859a):
- de_direction is the rank partial correlation of dp and dt with rank(PB_ref) held fixed. The depth artefact is monotone in PB_ref over most genes, so it is mostly partialled out. de_direction is robust to depth and sees (b) and (c).
- de_score is (overlap - null)/(1 - null), with the null the better of ranking genes by +PB_ref or by -PB_ref. A pure depth DE set is largely "explained" by the null; a depth term of the right sign adds overlap above it and one of the wrong sign removes overlap.
- severity is beta = (dt.dp)/(dt.dt) over the truth's DE genes, with r^2 = 1 - SS_res/SS_tot and SS_tot centred on dp; the floor applies when beta <= 0 or r^2 <= 0.01. When the truth carries a depth artefact, the artefact dominates its DE genes. A predicted depth term of the wrong sign, or one where the truth has none, is a near-uniform shift on those genes and puts severity at the floor.
- mmd_u and variogram: depth changes every cell's sparsity, and admixture shifts every cell.
Synthetic illustration (stylised vectors with arbitrary scales, no data: `synth_h/t3_pair.py`). With a reference component shared by dt and dp, de_direction is about 0.44 whatever the depth. "Factor plus Gata4-like depth" put severity at the floor in 78% of draws when the truth's depth ratio was 1 and in 100% when it was about 1.2, and in none when the depth matched. "Factor without depth" never hit the floor. The numbers illustrate the formulas only and estimate nothing at test.

### 3. How technical offsets behave across samples within one imaging series
Variance components, per gene g, for a sample on run r and type t:

    ln m_stg = mu_tg(stage) + a_s + e_rg + c_stg + b_stg + eps,     e_rg = sum_k L_gk beta_rk + xi_rg

- Within a series the loadings L (probe counts, bit assignment, abundance) are fixed by the library. The run coefficients beta_r (per-bit efficiencies, SNR, crowding level, segmentation granularity) are drawn per run, and samples on one slide share them.
- Size: with an ln-expression spread of SD about 1.5 over the panel (assumption), replicate correlations of 0.94-0.99 [1][15] leave a per-gene residual SD of 1.5 * sqrt(1 - r^2), about 0.2-0.5 ln units. In scorer units that is (0.2-0.5) * PB_g, enough to push many mid and high genes past 0.25. Run-to-run technical factors are therefore expected to be as large as the knockout factor itself (D10's WT-series factor had sd ln phi 0.456, `gata4_decisions.json`).
- Reference anchoring. Every contrast "X against W" has the factor F_X = e_X - e_W + Delta_stage(X,W) + b_X. Two contrasts against the same W share -e_W and W's own staging deviation. With independent runs of similar variance, Corr(F_G, Y_B) is about sigma_W^2 / (sigma_W^2 + sigma^2), about 0.5. If B shares G's slide it approaches the two replicates' agreement; if B shares W's run it falls to about 0 (only biology is left).
- Two models give different predictions. Model A has independent runs, so only the reference anchor is shared. Model B adds an offset common to all knockout runs. A predicts D11b > 0, vs95 about 0 and K2a about 0; B predicts all three > 0. Observed (committed): D11b 59.40 / 62.15 against WT cells 52.39 / 50.61; vs95 52.04 / 46.24; K2a +0.18 over its permuted control. That supports A.
- A does not explain R's failure. The WT-series interpolation scored 50.49 / 48.17 although it should share -e_W875. There are two readings: (i) R is noise- or development-dominated (three stages, per-label log-linear fits); (ii) the part that D11b shares with Gata4 is staging, with W875 younger than both knockout litters, which R's time interpolation removes by construction. The oracle NNLS stage weights fitted on Gata4 (E9.5 0.274, E8.25_late 0.010; `gata4_decisions.json` stage_alphas_mean) lean toward (ii), or toward a W875 offset aligned with the developmental axis.

### 4. What the committed files already say (repository facts)
- Two file streams (`docs/DATA_NOTES.md`). Every knockout file (Mab21l2, both Gata4 replicates) is dense, has a normalized_counts layer, is uncentred and carries slide, section, embryo, sample_id and perturb_batch columns. The WT files are sparse and centred, and E8.0, E8.25_late and E8.75 carry synthetic records (copies and averaged half-integer cells). The beta-catenin file is a knockout, so it almost surely comes through the knockout stream; W875 is a WT-stream file.
- The two Gata4 replicates sit on one slide (M164). Their agreement (D11a Spearman 0.724; D2 cross-replicate 68.05 / 64.09, with depth 74.50 / 72.15) therefore cannot separate run offsets from biology, and the ledger's rule that the replicate fold never counts is physically right.
- Depth (median counts per object) and objects per file: E8.25_late 233 (58,716); E8.5 166 (23,710); E8.75 WT 130 (24,826); E9.5 WT 75 (53,742); E10.5 56 (167,864); Mab21l2 95 (50,294); Gata4 82 and 68 (42,257 and 46,396). W875 sits on the WT stage trend (log-linear between E8.5 and E9.5 gives about 136 at E8.75).
- At equal stage the product of median depth and object count, a crude proxy for transcripts per stack, is nearly conserved: E8.75 WT 3.2M against Gata4 3.5M and 3.2M; E9.5 WT 4.0M against Mab21l2 4.5M (deduplicated). Gata4 has 46 and 44 levels at pitch 15 against W875's 35 at pitch 14, 920-1,050 objects per level against about 710, and about 0.8x the counts per level. Gata4's 0.55x is therefore "the same transcripts over more objects": either finer segmentation [6] or more tissue at lower capture per area. Mab21l2 shows the opposite: objects per level within 2% of WT E9.5's and about 1.12x the stack total.
- Gata4's per-label depth ratios are about 0.4-0.8 in non-Mesp1 labels (neural tube, neural crest, surface ectoderm, foregut) as well as in Mesp1 labels (`gata4_decisions.json` depth_ratio_pooled), so the deficit is a run property, not lineage biology.
- W875-specific artefacts: pitch 14 (most heart files use 15), 1,449 synthetic records (5.8%), and about 2x more vertical similarity than xy in E8.0-E8.75 (`experiments/coords/REPORT.md`). The scorer loads the reference as released, synthetic records included, while our predictions drop them; that small anchor is already shared by dp and dt.
- The knockout pipeline filters cells at 20 counts and 20 genes (TRIES, QC20). This one technical component transfers to beta-catenin with near certainty, but its lever is closed (a depth increase and a composition change).

### 5. What this means for the T3 pair
The test outcome hinges on three sample-level unknowns of the beta-catenin file that no released data can reveal:
- U1: was it imaged with W875? If so, the -e_W875 anchor is absent. My judgement is P(separate run) of about 0.85: different streams, and D11b shows Gata4's run differs.
- U2: its depth ratio to W875. Of the two knockout runs released, one is 0.55x and one 1.27x of its own stage-matched WT, and the per-object mechanisms differ. I put P(not Gata4-shallow, ratio above about 0.85) at about 0.35-0.45.
- U3: its litter's staging relative to W875. If the shared component is staging, its sign for beta-catenin is a fresh draw unless W875 itself is off-nominal.
s1 (Gata4 factor, Gata4 depth, composition) bets on U1 and U2 and carries the M164 run offset as noise. s2 (WT cells with composition) bets on nothing, and the synthetic check and the committed stand-in (severity at the floor in 4 and 10 of 10 seeds) say it will probably sit at the severity floor. The barred-branch factor (Mab21l2 against W875, no depth) shares U1 and U3 with s1 but not U2. So the two bets are coupled on the factor and decoupled on depth; the adversarial review of 2026-11-05 priced only the coupling.

Rough payoffs in T3 points, my judgement anchored on the committed Gata4 readings:
- s1: about 66 if U1 and U2 hold; about 50 if U1 holds and U2 fails (severity at the floor, mutant state worse); about 46 if U1 fails.
- s2: about 51 throughout.
- Barred-branch factor (s2'): about 58 if U1 holds (cf. D11b's 59.4 / 62.2 on targets with an unmodelled 0.55x depth); about 47 if U1 fails.
- E[max(s1, s2)] is about 58.7 and E[max(s1, s2')] about 60.6, before the U3 risk.

### 6. Target-free diagnostics that would tell which part of the s1 factor transfers (T3)
- D1, run fingerprint audit (metadata only). Per file: slide, perturb_batch, embryo, sample_id per level, pitch, levels, gaps, synthetic records, centring, matrix type, layers, obs schema, minimum counts and genes per cell, identical rows. Readout: whether Gata4 and Mab21l2 share a perturb_batch. If they do, the null cross-knockout transfers (vs95, K2a) show that batch membership does not carry per-gene offsets, so only the reference anchor transfers (Model A). If they do not, beta-catenin's sharing of the M164 offset stays possible but unknowable. The audit cannot see the hidden file; it only fixes the prior.
- D2, areal capture audit (structure only). Per section: tissue area (alpha shape of xy), objects per area, counts per area, object nearest-neighbour distance, and the counts-against-object-area slope where an area or volume column exists. Compare W875 with Gata4 and WT E9.5 with Mab21l2. Pre-registered decision: if Gata4 shows about 1.3x+ objects per area at about equal counts per area, the deficit is segmentation granularity; if counts per area are about 0.6x at equal objects per area, it is capture. Granularity that is a knockout-stream setting raises P(U2 holds); slide-specific capture leaves it near 0.6.
- D3, within-stack section structure (WT files only, W875 included). Per-level median depth, change points along z, and an SVD of per-level per-gene offsets within harmonised labels after a smooth z trend. Readout: the within-run section-to-section per-gene SD and its low-rank structure. If all sections sit on one slide, this is a lower bound on run-to-run variance and shows which gene covariates carry it.
- D4, factor anatomy against WT-only covariates. Regress the committed s1 factor (ln phi, 500 genes) on: ln PB_W875 (crowding and normalization); a spillover class from W875's own geometry (non-Mesp1 own level against the 15-NN Mesp1 neighbours' level); the WT developmental axis d_g; and a probe-count proxy from generic transcript annotation (needs an EXTERNAL_DATA row). The variance shares say how much is run-type technical (whose sign is unknown for the beta-catenin run), spillover (Gata4 lineage biology leaking into non-Mesp1 cells, which transfers only by phenocopy) or staging.
- D5, localising the shared component. Split F_G (s1's factor) and F_M|W875 (D11b's factor) into the part along d (WT E8.5 to E9.5 in non-Mesp1 labels, harmonised, share-weighted like the factor) and the orthogonal rest. Corr(rest_G, rest_M) estimates the technical W875 anchor; a_G * a_M * Var(d) estimates the staging part. Pre-registered decision: if the orthogonal part carries most of the covariance, U3 is minor and the straddle (lever 1) is positive; if the along-d part does, the straddle's value is about 0 or negative.
- D6, W875 staging placement (WT ladder only: E8.25_late, E8.5, E8.75, E9.5). Place W875 by a pre-registered estimator after a known-lag check, because the earlier timing estimator failed such a check. If W875 is young for its nominal stage, every other E8.75 embryo is older in expectation, the staging part regresses in a known direction, and U3 becomes a technical-like anchor.
Every diagnostic computes no score. D3 and D6 read WT files only. D1, D2, D4 and D5 read structure, or the factor vectors already written to `runs/t3_decomp/.../factors/` (whose sha256 values are in `gata4_decisions.json`; I did not open runs/). Whether that counts as a target read is the main agent's call.

### 7. T1 single-nucleus RNA: cross-stage library and batch effects
Model for library L, nucleus i of type t: E[c_Lig] = kappa_L * n_i * eta_nuc_g * rho_tg + alpha_Li * S_Lg, with ambient profile S_L and per-nucleus contamination alpha. The absolute contamination is roughly constant per droplet, so its fraction is higher in small nuclei [4]. X = log1p(1e4 c/N) over 32,285 genes. N runs to thousands, so more genes sit in the dense regime than on MERFISH, but most of the 32k remain sparse.
- What differs between the E10.5 base library and the E12.5 truth library: N and complexity (sparsity and the h(N) shift; closed as T1 depth); ambient alpha and S (nuclei background 17-35% against 3-11% in whole cells, varying between replicates, and in nuclei not composition-shaped [4], so not predictable from composition); isolation quality (cytoplasmic carryover); the embryos' sex mix (already held flat); doublets. All are fresh draws at E12.5.
- The reference-anchored part exists on T1 too. dt contains -e_E10.5. Our base copies E10.5, so it misses the "regression to the mean" of E10.5's own idiosyncrasy. That part cannot be identified for most genes, because E8.5 and E9.5 are embryo portions (tissue confound). Only tissue-independent modules (sex mix, preparation stress) could regress, and they are few genes, so this is worth about nothing.
- Feng's whole-cell direction carries three library differences that do not belong in a nuclei E10.5 to E12.5 change: a dissociation-stress difference between the E9.5 and E14.5 preparations [11][12]; an ambient difference that in whole cells follows each sample's composition [4][9]; and whole-cell against nuclear abundance for stable, cytoplasm-dominant mRNAs [13]. Large technical changes rank at the top of |lfc|, so they take de_score's predicted slots.
- T1 diagnostics:
  - T1-a, module audit of Feng's direction (Feng E9.5 and E14.5 plus released T1 E8.5 and E9.5): the share of each Feng group's top-100, top-300 and top-1000 |log ratio| genes in a published dissociation-stress signature, in mitochondrially encoded and ribosomal-protein annotation classes, and in genes that are ambient-dominated in Feng's own libraries (within-type mean below a fixed fraction of the gene's maximum-type mean). Compare with the shares among T1 E8.5 to E9.5 within-type changes (same nuclei assay) as the reference.
  - T1-b, nuclear detectability ratio at E9.5 (T1 E9.5 nuclei against Feng E9.5 cells in matched types): flags genes whose whole-cell changes nuclei cannot show.
  - T1-c, an E10.5 library fingerprint: it reads a counted T1 truth, so it is a try and not recommended.

## Candidate levers
Ranked by expected value to a top-3 overall place. None closes a 15-25 point gap: the measurement view mainly re-prices T3's risk structure. The realistic total is about +1.5 Overall from T3 and +0.3 from T1.

1. **T3 depth straddle: s2 becomes the barred branch's s1 (the Mab21l2 factor against WT E8.75, exponent 1, no depth, s2's current D7m composition), gated by D5/D6 and D2**
   - Board and scored group: T3:bcatenin. Response identity (de_score 0.30), response direction (de_direction 0.25) and response magnitude (severity 0.25); mutant state unchanged in expectation.
   - Mechanism: s1 and this factor share the W875 anchor (U1, U3); only s1 carries the depth sign bet (U2). The pair then straddles depth instead of the factor (section 5).
   - Expected effect: about +1.5 Overall, range -1 to +3. The gain is roughly P(U1) * P(U2 fails) * (+6 to +8) - P(U1 fails or U3 flips) * 3. The adversarial review of 2026-11-05 put the swap at -0.2 to +0.1 because it treated the bets as fully coupled.
   - Legality check:
     - The T3 rules are frozen, so this changes a frozen rule outcome (notes 3, 6, 8) and only the captain's ruling can make it; `rules.yaml` is not edited.
     - No closed lever is touched: the composition is unchanged, there are no cell-number or size changes, no T1 lever is involved and no cross-modal input is added.
     - It is Gata4-free, so it keeps the organiser-bar insurance.
     - It reads no score.
     - Research rules: data are the released Mab21l2 and W875 files only, with no gene lists and no barred literature.
   - Validation route: captain's blind bet only. Fold B (Mab21l2 to Gata4 rep2) is D11b's own reading (it missed its bar of 60 on rep1 by 0.6) and is selected-on for this candidate. Fold A (lever 4) tests the mechanism in mirror.
   - Effort: 0.5 day. The barred branch is already built by the build script.
   - Main risk: U3. The shared part may be litter staging (oracle stage weight 0.27 on E9.5 cells). If the beta-catenin litter is younger than W875, the factor's sign flips in both uploads, and today's s2 is the better insurance. The Mab21l2 factor also carries 0.75 day of development, which is pure bias at E8.75.

2. **T3 diagnostics D5 and D6 (shared-component localisation and W875 staging placement)**
   - Board and scored group: T3. Decision value only, for lever 1 and for T3's band in the 209 estimate; no scored group moves.
   - Mechanism: run offsets do not follow the developmental axis, but staging does.
   - Expected effect: no points by itself. It sets lever 1's sign: about +1.5 to +2.5 if the shared part is orthogonal to the axis, about 0 or negative if it lies along it.
   - Legality check:
     - Reads WT files plus the factor files already written; no score is computed and no test file is touched.
     - Reading the factor vectors is the main agent's call: they derive from the unlocked Gata4 files and are already inside s1.
     - No closed lever or frozen rule is touched.
   - Validation route: target-free analysis with a decision rule written before it runs.
   - Effort: 1 day.
   - Main risk: stage placement on this ladder failed a known-lag check before, and the developmental axis carries the WT runs' own offsets.

3. **T3 diagnostics D2 and D1 (areal capture audit and run fingerprint audit)**
   - Board and scored group: T3. Decision value: P(U2 holds), so whether s1's depth bet or lever 1 deserves the weight.
   - Mechanism: separates segmentation granularity (objects per area) from capture (counts per area) in Gata4's 0.55x and Mab21l2's 1.27x.
   - Expected effect: no points. It moves P(U2 fails) between about 0.15 (a granularity setting shared by every knockout slide) and about 0.45 (slide-specific capture).
   - Legality check:
     - Structure-only reads (coordinates, per-object totals, level ids, slide and perturb_batch metadata), the `docs/DATA_NOTES.md` intake standard.
     - No expression profile, no score and no test file are involved; whether the Gata4 and Mab21l2 structure reads count as target reads is the main agent's call.
   - Validation route: target-free analysis.
   - Effort: 0.5 day.
   - Main risk: the hidden file's slide cannot be known; the audit only classifies mechanisms.

4. **T3 Fold A mirror of D11b (try K3a)**
   - The candidate: WT E9.5 cells times the Gata4-against-WT-E9.5 non-Mesp1 factor (exponent 1, no depth, no composition), paired with K2a (anchored on W875) and our WT copy (50.74).
   - Board and scored group: T3 fold Gata4 to Mab21l2; the response groups.
   - Mechanism: Model A says a factor anchored on the target's own reference carries -e_ref, while K2a, anchored elsewhere, carried nothing.
   - Expected effect: K3a beats K2a; against the copy, +1 to +6 with P(pass) about 0.4-0.6. A pass is the first counted support in the mirror direction for lever 1's mechanism. A fail is weak evidence against it, since 0.75 day of development is in the factor and W95 may be less idiosyncratic.
   - Legality check:
     - A counted try, registered before it runs (T3 goes from 3 to 4).
     - Its inputs are inside the fold's allowed set (Gata4 files, W95; W875 for the control only).
     - No composition or depth; it changes no upload by itself; it reads no test score.
   - Validation route: it is the counted reading on Fold A. Adoption also needs Fold B, where the matching candidate is D11b's earlier reading, so whether that counts is the main agent's call.
   - Effort: 0.5-1 day, using the K2a machinery with the factor source swapped.
   - Main risk: it spends Fold A on a mechanism test that the development confound may sink.

5. **T1 technical-module hygiene of Feng's direction (T1-a, then hold flat any module whose share clears a pre-set margin)**
   - Board and scored group: T1:test. Mainly DE recovery (de_score 0.25); direction (de_direction 0.25) and distribution (mmd_u 0.30) are minor.
   - Mechanism: stress, ambient and whole-cell-only changes in Feng's two libraries do not belong in nuclei; removing them frees de_score's predicted slots.
   - Expected effect: +0.1 to +0.5 Overall. The sign is uncertain where modules carry real biology.
   - Legality check:
     - Generic module sets need EXTERNAL_DATA rows with provenance. None is built from a Ctnnb1 or Gata4 perturbation, and the rows record the published source rather than hand-typing genes.
     - Feng's E9.5 and E14.5 are outside the window.
     - No depth or composition is touched and no cross-modal input is added.
     - Changing Feng's direction is a T1 entry change for the captain, since no counted check exists (Feng's E9.5 sits inside every T1 backtest window).
   - Validation route: captain's blind bet, with T1-a as its target-free basis.
   - Effort: 1-1.5 days.
   - Main risk: the E12.5 library's own stress or ambient may by chance move like Feng's, and immediate-early genes can be developmental.

6. **T1 ambient-aware Feng direction (the direction recomputed on decontaminated Feng counts; CellBender if Feng's raw droplet matrices exist, otherwise DecontX)**
   - Board and scored group: T1:test, DE recovery and distribution.
   - Mechanism: in whole cells, ambient follows each sample's composition, so part of Feng's within-type change is a change in the soup [4][9][10].
   - Expected effect: 0 to +0.5 Overall.
   - Legality check:
     - Legal external data, with an EXTERNAL_DATA row for the tool.
     - Same captain route as lever 5; no closed lever is touched.
   - Validation route: captain's blind bet.
   - Effort: 2-3 days.
   - Main risk: GEO-processed data may lack empty droplets, and correction can distort fine structure [4].

Closed or declined:
- T3 spillover-aware factor (pruned 2026-11-04).
- T3 QC20 cell filter (closed: depth and composition).
- An absolute knockout-depth rule (pruned as a lucky fit).
- Any change to s1's depth or factor exponent (frozen D4/W1; captain only; not recommended).
- T1 depth matching (closed).
- T1 ambient matching (no legal anchor).
- T2: the depth line is closed, and there is nothing new from this angle.

## Deep Read Notes
- [1] Wang 2018. A high-abundance library (about 14x earlier totals) gave MERFISH/smFISH 0.21 +- 0.04 (median 0.16, r = 0.75) unexpanded and 1.05 +- 0.09 (r = 0.96) after 2.3x expansion; replicate preparations r = 0.94, expanded against unexpanded r = 0.88. The loss comes from overlapping spots of neighbouring molecules. Simulation puts about a 50% loss at 3-4 molecules per square micrometre per 1-µm section. Detection efficiency is a run-times-density property, not a constant.
- [2] Atta 2024. 100-gene panels skewed toward one region, drawn from a 483-gene MERFISH panel, raised the library-size scaling factors there (deflated values) and lowered them elsewhere: DE errors up to 13% and fold-change sign switches up to 19% on MERFISH (up to 30% and 60% on other platforms). DESeq2 and TMM were milder. Volume normalization was panel-independent but rests on segmentation and partial z-capture (simulated 7 planes x 1.5 µm).
- [3] Bhuva 2024. 25 samples on Xenium, CosMx, Visium and STOmics (no MERSCOPE). Detections per cell differ by region after a Poisson GLM adjustment for cell count. Permeability is proposed for the sequencing platforms; cell size, segmentation, thickness and DE are untested. Library-size normalization hurt domain finding (sctransform worst); the authors advise against it by default.
- [4] Janssen 2023. A two-subspecies design measured true background: median about 11% and under 3% in two whole-cell replicates, about 35% and 17% in two nuclei replicates. It is mostly ambient (profile rho 0.73-0.85 with empty droplets; barcode swapping under 10%), and in nuclei the ambient mix is shifted from the cells' mix. CellBender was the most precise; correction raised marker fold changes and changed clustering little.

## Gaps
- Searched for and not found: any study of run-to-run per-bit or per-round efficiency drift and its per-gene structure in MERFISH. Only single-run error rates exist [7][8]. No public codebook exists for this unpublished 500-gene panel, so a bit-structured decomposition of the factor is not identifiable from released data.
- Searched for and not found: a gene-level variance-components study of section-to-section or slide-to-slide technical variation in serial-section MERFISH. Only replicate correlations of per-gene means exist [1][15].
- Searched for and not found: a study of how segmentation granularity and section pitch set per-object depth in serial-section 3D MERFISH stacks. Only Atta's partial-volume simulation [2] and Baysor's object-count result [6] bear on it.
- Not read in full: Chen 2015 (an anti-bot page blocked PMC) and Mitchel 2026 (paywalled; abstract only).
- Counter-claims against lever 1:
  - The review of 2026-11-05 (wf_6b41d823-d5a): both bets ride the WT E9.5-against-E8.75 offset that s1's factor shares, 83% of the swap's stand-in gain was severity un-flooring, the downside was floored at K2a, and three of four ruling switches would break. I answer only the coupling point, through depth; the ruling-switch objection stands untested.
  - If the shared part is staging (oracle weight 0.27 on E9.5 cells) and the beta-catenin litter is delayed, both factor uploads lose and today's s2 is the right insurance.
  - R's failure may mean that the technical W875 anchor is small and that D11b's gain is a Mab21l2-specific coincidence.
  - If every knockout slide at E8.75 was segmented with Gata4's settings, P(U2 fails) is near 0.15, and lever 1 adds about nothing.

## Compliance log
- None. No barred-class item was met: every search was technical and no result concerned beta-catenin, Gata4, Wnt loss of function, in-window development, the leaderboard or another team.
- Transparency notes:
  - Two tool outputs were saved by the harness under the account folder, a PubMed full-text response and a PDF; I did not open either (`AGENTS.md`).
  - I did not open data/, runs/ or boards/, or any docs/research file except `digest_2026-10-19.md` (Parts 0-1 and sections 2.4, 2.5 and 4.2).
## END
