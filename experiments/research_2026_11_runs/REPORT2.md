# Depth, run profiles and two uploads: what the literature says about our measurements

2026-11-12 · Mode: Lightweight · AS_OF 2026-11-12 · Approved sources: 90 · Official sources: 16 of 90 (17.8%, below the 20% gate)

Final report (P6–P7), edited by Claude (AI) from the P5 draft, the registry, two counter-reviews and three screened notes; no web page opened, no code run on competition data.
PRIVATE: contains our measurements from released competition data. Keep it out of every outbound request.

## Summary

1. **Check one fact first.** Read from the experiment record what the render's gains on the interpolation board (+2.6 to +7.4 on three practice targets) were measured against, and whether its factor was fitted without the target. If the baseline was not slot one's pool, the case for a new second upload has no number.
2. **Then score one frozen table (recipe 1).** Slot one, the current slot two, the render, each bracket's own cells and two ablations, scored as complements to slot one under a rule signed off beforehand. The render's gross value is +1.8 to +2.4 board points under equal thirds (derived; +1.4 to +3.1 across weightings; Low confidence), before subtracting what the current slot two earns.
3. **Ask the captain one question (recipe 2).** Depth matching gained +9.7 and +22.4 in our zero-gap test on the scRNA board, the largest effect we measured, but that lever is closed. A shallow specialist's ceiling over the pooled slot one is 15.25 points. Build nothing before a ruling.
4. **Do not build** denoiser depth lifting or a predictor of the unseen imaging run's profile. The notes found no validated method for either, and both need closed or external inputs.

Firmness: High that depth mismatch moves the board score under idealised thinning (one controlled test). The cause of the MERFISH thinning failure is not established. That two different specialists beat two variants of one hedge is an unpublished derivation (Low).

## Reading key

[n] cites the registry; numbers beside a citation come from the screened notes. "Re-read" marks a counter-reviewer's own reading of a source; "we measured", our findings as the brief gave them to the draft's author (I could not open the brief); [unverified], what the notes do not support. A *pool* mixes cell sets in one upload; a *specialist* targets one scenario; a *complement* is the best second upload given the first. *The render* is the existing builder: the shallower bracket's cells times a per-gene factor, topped up to the deeper bracket's depth from kNN-smoothed profiles. The *seed spread* is how far a score moves when the scorer redraws its random gene pairs.

## 1. Depth in droplet libraries

**Log-normalised values keep a depth signature.** Halving the UMIs of half of one monocyte population produced more than 2,000 differentially expressed genes on log-normalised data, against 11 on Pearson residuals [7]; the 500 deepest against the 500 shallowest cells of one type gave 1,490 false positives [9]. Cell-to-cell distance tracked depth difference at 0.54 and 0.71 [9], clustering agreement fell from about 0.60 to 0.35–0.38 when half the cells were thinned to 25% depth [26], and depth variation alone gave an independent gene pair a correlation of 0.14 [10].

**Thinning fits the UMI noise model; its match to a real shallow library is untested.** Spike-in overdispersion was "effectively zero in six of the datasets and less than 0.015 in the other two" [19], and a thinned negative-binomial count keeps its dispersion [17][18]. But background RNA was 3–35% of a cell's UMIs [11], molecular recovery at fixed saturation ranged from 39% to 70% [15], and read-matched libraries differed by 11–22% in median UMIs per cell [16]. Rendering deeper has still less support (recipe 7).

**Against our findings.** Matching the base to a thinned truth gained +9.7 at 0.68 of full depth and +22.4 at 0.35, over an unmatched base (derived: the pooled base already scores 97.20 and 84.75 there). The literature predicts an effect on all three metric families [7][9][26][10]; we show one, variogram skill of 0.51 and 0.24 for a full-depth base at those truths, which alone could account for the gain if the families weigh about equally [unverified]. Two libraries of one kind scored 66 (median 4,799 counts per cardiomyocyte) and 80 (9,441) as targets, but that pair is uncontrolled [unverified].

The notes found no source measuring an RBF-kernel MMD against depth; the nearest proxy, an energy distance, "increases as the UMI count per cell increases" [2] (design [unverified]). The variogram score was not searched. Our five variogram points are two bases, not one curve (applicability review's reading): a base at 0.80 of full depth, presumably the pool, scores 0.84, 0.99 and 0.40 against truths at 1.0, 0.68 and 0.35; a full-depth base scores 0.51 and 0.24 against the thinned truths. No single-depth base shallower than its truth was measured, so a render's tolerance on that side is unknown.

**Confidence: High** that depth mismatch moves the board score under idealised thinning (one controlled test); **Medium** that it moves each metric family; **Low** on the size of any gain on a hidden target, where the only figure is a ceiling (recipe 2). Opposing reading: Key controversy 1.

## 2. Run-to-run detection in imaging data, and the unseen run

**Total counts vary between runs, and more between protocols.** Across 641 mouse-brain MERSCOPE sections, transcript density had an SD of about 0.41 of its mean and differed about 2-fold between sections of one mouse [41]; an academic atlas saw about 30% variation between experiments [44]. In cultured cells an older and an optimised protocol differed about 14-fold in RNA copies per cell [42]. The largest effect in [41] is section-wide and "presumably multiplicative", which uniform thinning would mimic.

**Per-gene effects have mechanisms, with limits.** A barcode-derived "readout position score" explained much of the gene-to-gene difference in detection between protocols in cultured cells but only "partially explained the variation" in tissue, and covering the reagents "almost completely removed the apparent drop" in later rounds [42]. A dropped readout image removes the 60–104 genes sharing it [41]. Sensitivity "is correlated with the number of probe sets used" and falls where transcripts crowd [54]. Part of the error is additive background [42]. One study models the effect as cell-type-specific [43], another assumes it shared and names that a limitation [50]; the shared part is unmeasured.

**No per-gene size, and no correction without the target.** Replicate agreement is reported only as log-scale correlations, 0.95 to 0.99 for sections and cultures [44][45][46][48] and at least 0.7 for separate tissue cores [47], which the spread between genes dominates. Routine normalisation removes one global factor [44][53][49]. Every correction in the sources found needs the target run's cells [50][43][61], and the notes found no method that predicts a profile before the run is seen.

**Against our findings.** (a) Uniform thinning of deeper MERFISH cells lost points even at the right depth. Gene-set losses and additive background [41][42] are compatible with that, but the cause is not established: the largest published effect is section-wide, which thinning would mimic. (b) The render gained +2.6 to +7.4 on the three targets sharing the profile it renders towards, but it also tops cells up and smooths them, and it gained +1.3 on the target sharing neither profile: the three parts are not separated. (c) Successive quarter-day changes anti-correlate at −0.44, near the −0.5 of independent per-sample deviations (derived): the change is mostly sample-specific (run, embryo, staging, section or estimator noise), not a smooth trend [unverified].

Our three large steps have a log-factor SD of 0.79 to 0.94 (2.2- to 2.6-fold in natural logs, 1.7- to 1.9-fold in base 2; the base is not stated), the other three 0.31 to 0.57. Three readings stay open: a protocol or batch change as in [42], biology, or estimator noise, which at 30 to 62 counts per cell and a hundred cells gives a per-cell-type log ratio a noise SD of about 0.3 to 0.4 (applicability review's arithmetic). One sample's factor is nearly one vector against three later samples (Spearman 0.95 and 0.82; which comparisons is not stated), as development would also produce, yet it is unrelated to neighbouring factors (−0.23 to +0.05) [unverified].

**So we hedge:** we found no way to predict an unseen run, and correcting it needs its cells. **Confidence: Medium** for the mechanisms (full-text sources [41][42][43]); **Low** that profiles are run-specific in our data. Opposing reading: Key controversy 2.

## 3. Two uploads where the better counts

**Jointly trained predictors specialise under a best-of-several loss:** under squared error each becomes the mean of its own region of outcomes, while a single predictor sits between the modes [68][66][67]. The measured benefit is mixed: 0.822 against 0.757 for two specialists against two independent models [68], but 59.7% to 60.0% for a second hypothesis elsewhere [66].

**Pooling has exact algebra, under conditions.** For kernel scores a pool scores the weighted average of its components minus a disagreement term that does not depend on the outcome [62]; a kernel score is half a squared MMD [63]. The note derives, beyond the sources, that a pool with weight w on the wrong component pays w² of the divergence between the two [62][63][65]. When no component is the truth, as with our imperfect renders, a pool can beat a single model (by 10–16% in one study [63]) or fall between its components (1.19 against 1.10 and 1.47, lower being better [62]).

**For two slots there is only the note's derivation,** unpublished and awaiting an independent check: the best pair splits the scenarios into two groups and submits each group's probability-weighted pool, so two pure specialists for two scenarios and, for three, a specialist plus a pool of the pair with the smallest probability-weighted divergence [66][67][62]. It assumes exact scenario models and a squared distance [63][64][77]; whether our scorer squares its distance was not established.

**Against our findings.** Both slot ones are pools (both brackets' real cells; three depths), right for one slot under expected squared distance [62][63]; the scRNA pool scores 84.75, 97.20 and 91.45 at depths 0.35, 0.68 and 1.0. The render behaves like a specialist: +2.6 to +7.4 on one side, −2.3 and −6.6 on the other. Both pairs put a variant of the same pool in slot two (the earlier geometry; a higher share of one cell type). The scRNA pair has no depth diversity, although the pooled base varies by 12.45 points across depths (derived); the interpolation pair has none on the axis where the render moves scores. The brief gives no figure for geometry or cell-type share, which does not make them small [unverified].

By the claims review's re-read, [63] lists the variogram score among kernel scores, so the pooling identity is sourced for its population form; its fit to our scorer, with any unbiased estimator and the transform skill = d_floor / (d_floor + d), is untested (recipe 3). That transform is convex in distance, which favours specialists (derived) [unverified], but that is no advice to change slot one: on the composite score the pool averages 91.13 over the three depths against at most 89.3 for a full-depth base (derived), and two pure specialists would drop the only upload that covers "neither".

**Confidence: Medium** for the pooling algebra (exact mathematics in two preprints); **Low** for the pair rule and its transfer to our composite score. Opposing readings: Key controversies 3 and 4.

## 4. Key controversies

1. **Does depth matter here, and is thinning a fair model of it?** Claim: log-normalised values carry depth into differential expression, distances and gene-pair correlation [7][9][26][10], and thinning matches the UMI noise model [17][18][19]. Counter-evidence: one benchmark found the shifted log with PCA kept a cell's neighbours as well as depth-aware transforms [1], and leading components survive shallow sequencing [30]; real libraries also differ in background and recovery [11][15], though clustering was "fairly robust towards background noise" [11] and fully concordant between two sequencing runs [16]. Standing: [1] and [30] measure neighbour retention, not whether two cell sets are one distribution, and in [1] the size factor "remained a strong variance component". If the hidden library differs in background or chemistry, a depth render chases the wrong variable.
2. **Is the imaging run effect one global scalar?** Claim: nearly; replicate correlations are 0.95 to 0.99 [44][45][46][48], same-panel experiments show "a strong similarity in gene-specific detection efficiency" [49], and the largest effect in [41] is section-wide. Counter-evidence: gene-set losses [41], round-dependent and additive background [42], crowding [54], cell-type-specific examples [43], and segmentation errors "frequently dominating the results" of cell-state analyses [59]. Standing: on the scalar reading our factor is mostly biology, segmentation or noise, and a sample between the brackets wants a partial render; the −0.44 argues against that without excluding it.
3. **Diversify the two uploads, or take the two best validated?** Claim: best-of-several objectives reward specialists [66][68], and one practitioner report picks from different model classes [79]. Counter-evidence: a second hypothesis added almost nothing [66]; such specialists are overconfident [69]; a meta-analysis of over one hundred competitions found little overfitting to public holdouts [80]; with five submissions allowed, breadth of sampling rather than ranking skill appeared to drive results [85]; rank-based reward can favour neither the honest pool nor the honest specialist [86][87]. Standing: practitioner advice is thin and split (a course hedges across two selection criteria [78]); the notes found no source on a best-of-two rule.
4. **Should the render take slot two on the interpolation board?** Claim (the draft): slot two is spent on the smaller uncertainty, and the render adds about +2 points gross. Counter-evidence (applicability review, from our data): slot one already holds the deeper bracket's real cells; the practice gains come from targets known afterwards to sit on that side; the gain may be a depth effect; and over-fitting a noisy selection criterion often costs about as much as the differences between learning algorithms [82]. Standing: recipe 1 decides.

## 5. Recipes, ranked

Ranked by what each can change in an upload decision for the work it costs; the order follows the applicability review. Gains are board points estimated on practice data, never measured on a hidden target. Recipes 3 to 5 have no measured gain, so their order is Low confidence; 6 and 7 are not worth doing. The other MERFISH board is not assessed.

### Recipe 1. Baseline check, then one frozen complement table (interpolation board) [draft 1]

- **Do:** (0) Read from the experiment record what the render's gains were measured against, whether its factor was fitted without the target, and how each practice target got its scenario label. (1) Freeze the rows before scoring, each written as an operation (whose cells, rendered towards which sample, at what depth ratio, in which geometry): slot one; the current slot two; the render in slot one's geometry; each bracket's own cells; top-up without factor; factor without top-up. The half-strength render is dropped: the brief lists no scenario for it. (2) With fixed seeds compute per target max(slot one, candidate) − max(slot one, current slot two). Weight the three scenarios equally, show each scenario's row, and stress with half the weight on "earlier". Report the seed spread, each leave-one-target-out result, the factor's split-half reliability and the top-up's share of a rendered cell.
- **Moves:** gross gain of the render as complement, if the baseline was slot one's pool: +2.6 to +7.4 on three "later" targets (possibly one batch), +1.3 on the one "neither" target, zero on the two "earlier" targets, where it lost 2.3 and 6.6. That averages +1.8 to +2.4 under equal thirds and +1.4 (half on "earlier") to +3.1 (the practice mix) across weightings (derived); each range reflects one unstated middle gain. It is gross of what the current slot two earns and an upper bound until the factor is refitted without the target. Low confidence.
- **Kill:** the baseline was not slot one's pool; or the render's advantage over the current slot two as complement is below the seed spread; or it changes sign when any one target is left out.
- **Legal:** no external input. The render raises counts per cell on a MERFISH board, closed if the "sizes" rule covers per-cell depth there: ask. Single-bracket rows may fall under the closed "cell numbers"; then they are diagnostic only. Slot one is taken as fixed; opening it, or replacing slot two, needs the captain's ruling. The decision rule is signed off before the table is scored, by whoever approves decision rules here; a choice that turns on the weighting is the captain's.

### Recipe 2. One question to the captain: reopen depth on the scRNA board? [draft 5]

- **Do:** Ask one neutral question, giving the ceiling and assumptions below; build nothing before the ruling. If reopened: thin the pooled base [35][36] towards the shallow library's totals per cell type, pairing cells with target totals by rank, and aim at UMI totals, not reads [34][16].
- **Moves:** mainly variogram skill. A ceiling, not an estimate: 15.25 points if the truth sits at the shallow depth (derived), nothing where slot one wins. Take the odds from the depth spread of the released libraries, not from our three test depths.
- **Kill:** the captain does not reopen depth; or, on the real pair (previous stage to the 4,799-count library, which scored 66), the thinned base recovers under half of the 14-point gap. That tests thinning on one real shallow library, not whether the hidden library is shallow.
- **Legal:** needs the closed depth lever on the scRNA board and the captain's ruling on a second upload, which would displace the cell-type-share hedge.

### Recipe 3. Pooling-identity diagnostic [draft 2]

- **Do:** After recipe 1's step 0. For MMD and the variogram score compute d1, d2 (each bracket against the truth) and d12 (bracket against bracket) as squared distances (square the scorer's output if it reports the root), with one bandwidth, one PC space and one library-size target. Predict the pool's distance as w1·d1 + w2·d2 − w1·w2·d12 and compare with the scored pool [62][63]. First rule out what breaks the identity for other reasons: a library-size target or bandwidth set from the scored sample, and the unbiased estimator.
- **Kill:** with those ruled out, predicted and scored distances differ on two targets by more than a tolerance fixed beforehand as a fraction of d12 (default 10%) and by more than the seed spread.
- **Legal:** the analysis needs nothing. Acting on the best-weight formula [62][63][65] would tune a weight on one or two practice targets and change the bracket split, probably closed under "cell numbers": ask first.

### Recipe 4. Sharper render, tuned on a round trip [draft 4]

Only if recipe 1 puts a render in a slot.

- **Do:** (a) Measure the ceiling by a round trip on the deeper bracket: thin its cells to the shallow depth, render them back, and score them against a held-out half with the board's metrics. The shortfall from the split-half score bounds what a sharper render can recover and picks the neighbourhood size without a practice target; molecular cross-validation [20] is not used, because it optimises the per-cell mean, not the cell-to-cell spread the variogram and MMD measure (applicability review). (b) Draw the top-up from a Gamma-Poisson with per-gene dispersion, since deeper counts should be sampled, not read off a mean [3]. (c) Read from the builder whether the factor is one shared vector or per cell type; a per-cell-type factor waits for its split-half reliability, with any threshold in molecules per gene [unverified]. The ceiling covers smoothing and top-up only.
- **Kill:** the round-trip shortfall is within the seed spread; or, for (b), dispersions come out near zero, likely because sparse counts look Poisson [4].
- **Legal:** as recipe 1.

### Recipe 5. Earlier-side specialist by gene-wise thinning [draft 3]

Build only if recipe 1 shows the pair is short on the "earlier" side and that bracket's own cells do not cover it.

- **Do:** thin the deeper bracket's cells gene by gene towards the shallower bracket's profile and depth. Genes that would need a keep probability above 1 (roughly 3% to 10% at the largest depth step, applicability review's arithmetic) get a Poisson top-up, so it imputes for those genes only. The gain is unknown: thinning keeps the count family [17][18], but uniform thinning already lost points at the right depth and thinning cannot create additive background [42].
- **Kill:** it does not beat both slot one and the current slot two on the two "earlier" targets by more than the seed spread.
- **Legal:** as recipe 1; it also needs slot one or displaces the other specialist.

### Recipe 6. Predicting or correcting the unseen run's profile: not worth doing [draft 7]

- **Do:** regress the per-gene factor on barcode round or probe number [42][54], or apply RCTD- or Crescendo-style correction [50][43].
- **Kill:** not tested. Rejected because every correction in the sources found needs the target's cells [50][43][57][61] and codebook tables are an external input. Our measurement (neighbouring factors unrelated, −0.23 to +0.05) shows one sample's factor does not predict its neighbour's; it does not test a covariate. Whether the large factor tracks in-data abundance is unchecked [unverified].
- **Legal:** codebook and probe tables are an external input needing disclosure; calibrating MERFISH against external sequencing data is cross-modal and closed.

### Recipe 7. Depth lifting with denoisers: not worth doing [draft 6]

- **Do:** denoise a shallow library and sample it deeper.
- **Kill:** not run. The route is untested on any distributional metric and on our scorer. Against a deep reference five denoisers "recovered nearly all zero expression values across all gene groups, regardless of sequencing depth" [21]; three produced 20–80% false positives on permuted data [23]; a supervised lifter raised correlation with the deep profile only from 0.776 to 0.800, and its authors say transfer to real low-depth data "has not yet been directly established" [29]. Those are denoised means; sampled renders have only FISH checks on 11–18 genes [3][19][27]. Rejected because the lever is closed and nothing validates the route.
- **Legal:** a depth change on the scRNA board (closed); a pretrained model would be an external input.

**Work order.** Recipe 1's step 0 and recipe 2's question now; then recipe 1's table; recipe 3 after step 0; recipes 4 and 5 only on their conditions. Freezing both pairs stays a live option (Key controversy 4).

## 6. Verification

- **Citations.** 131 instances over 54 sources in sections 1 to 5, all in the registry's approved list (numbers 1 to 90). The one dropped source has no number and is not cited. The most-cited source, number 42, has 9 instances (6.9%), below the 30% limit.
- **Gates.** Approved sources 90 (at least 6: pass). Hostnames 23 (at least 3: pass). Official sources 16 of 90, or 17.8% (at least 20%: **fail**; 6 of the 54 cited here). No new source may be introduced at this stage, so the shortfall stands.
- **Spot-checks against the screened notes.** I checked each literature figure above against the notes. Eight are listed; four failed in the draft and are fixed.
  1. Sources 7, 9: more than 2,000 against 11 genes; 1,490; 0.54 and 0.71 (task-a, findings 6, 7). Pass.
  2. Sources 11, 15, 16: 3–35%; 39% to 70%; 11–22% (task-a, findings 1, 2, 4). Pass.
  3. Sources 41, 44: 641 sections, SD about 0.41 of the mean, 60–104 genes; about 30% (task-b, findings 1, 4). Pass.
  4. Sources 68, 66, 62, 63: 0.822 against 0.757; 59.7% to 60.0%; 1.19 against 1.10 and 1.47; 10–16% (task-c, findings 1, 6 and deep-read note). Pass.
  5. Source 19, "six of nine": the note's quotation covers eight of nine datasets. Fail; now quoted as in the note.
  6. Source 43, "needs at least 100 cells": the note says "fitted ... with at least 100". Fail; deleted.
  7. Correlations "0.95 to 0.99": the note also has at least 0.7 for tissue cores (source 47). Fail; restated.
  8. Recipe gains "0 to +1" and "+0.5 to +2": no trace in any note. Fail; deleted.
- **Not traceable.** Our own measurements come from the brief, which I could not open; they are carried from the draft, and every derived figure is recomputed. The render's description and the "applicability review's arithmetic" rest on that review's reading of the brief; the "re-read" statement, on the claims review's reading of source 63. The claims review could not re-open sources 2 and 41.
- **Data rules.** No external measurement from the barred embryonic-day windows and no knockout finding is used; sources 3 and 20 are cited for method only.
- **Length.** 4,453 words by whitespace count, the reference and dropped lists excluded as in the registry's count of the draft; the whole file is 6,475 words.

## References

The registry's approved list under its own numbers. "Not cited" marks an approved source not used above; other flags show how a source was read.

[1] Comparison of transformations for single-cell RNA-seq data (Ahlmann-Eltze and Huber, Nat Methods). https://pmc.ncbi.nlm.nih.gov/articles/PMC10172138 (academic; as of 2023-04; authority 9/10)

[2] scPerturb: harmonized single-cell perturbation data (Peidli et al., Nat Methods). https://www.nature.com/articles/s41592-023-02144-y (academic; as of 2024-01; authority 9/10)

[3] SAVER: gene expression recovery for single-cell RNA sequencing (Huang et al., Nat Methods). https://doi.org/10.1038/s41592-018-0033-z (academic; as of 2018-06; authority 9/10)

[4] Comparison and evaluation of statistical error models for scRNA-seq (Choudhary and Satija, Genome Biol). https://doi.org/10.1186/s13059-021-02584-9 (academic; as of 2022-01; authority 9/10)

[5] Determining sequencing depth in a single-cell RNA-seq experiment (Zhang, Ntranos and Tse, Nat Commun). https://doi.org/10.1038/s41467-020-14482-y (academic; as of 2020-02; authority 8/10; not cited)

[6] Feature selection and dimension reduction for single-cell RNA-Seq based on a multinomial model (Townes et al., Genome Biol). https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6927135/ (academic; as of 2019-12; authority 9/10; not cited)

[7] Normalization and variance stabilization of single-cell RNA-seq data using regularized negative binomial regression (Hafemeister and Satija, Genome Biol). https://pmc.ncbi.nlm.nih.gov/articles/PMC6927181 (academic; as of 2019-12; authority 9/10)

[8] Overcoming systematic errors caused by log-transformation of normalized single-cell RNA sequencing data (Lun, bioRxiv preprint). https://www.biorxiv.org/content/10.1101/404962v1 (academic; as of 2018-08; authority 7/10; preprint; not cited)

[9] Depth normalization for single-cell genomics count data (Booeshaghi et al., bioRxiv preprint). https://www.biorxiv.org/content/10.1101/2022.05.06.490859v1 (academic; as of 2022-05; authority 6/10; preprint)

[10] Cell-type-specific co-expression inference from single cell RNA-sequencing data (Su et al., Nat Commun). https://doi.org/10.1038/s41467-023-40503-7 (academic; as of 2023-08; authority 8/10)

[11] The effect of background noise and its removal on the analysis of single-cell expression data (Janssen et al., Genome Biol). https://doi.org/10.1186/s13059-023-02978-x (academic; as of 2023-06; authority 9/10)

[12] Massively parallel digital transcriptional profiling of single cells (Zheng et al., Nat Commun). https://doi.org/10.1038/ncomms14049 (academic; as of 2017-01; authority 9/10; not cited)

[13] EmptyDrops: distinguishing cells from empty droplets in droplet-based single-cell RNA sequencing data (Lun et al., Genome Biol). https://doi.org/10.1186/s13059-019-1662-y (academic; as of 2019-03; authority 9/10; not cited)

[14] Biology-inspired data-driven quality control for scientific discovery in single-cell transcriptomics (Subramanian et al., Genome Biol). https://doi.org/10.1186/s13059-022-02820-w (academic; as of 2022-12; authority 8/10; not cited)

[15] Sequencing saturation does not uniquely determine molecular recovery in UMI transcriptomics (Wilson et al., Bioinformatics). https://doi.org/10.1093/bioinformatics/btag593 (academic; as of 2026-08; authority 7/10)

[16] Comparative performance of the BGI and Illumina sequencing technology for single-cell RNA-sequencing (Senabouth et al., NAR Genom Bioinform). https://doi.org/10.1093/nargab/lqaa034 (academic; as of 2020-05; authority 8/10)

[17] Negative binomial count splitting for single-cell RNA sequencing data (Neufeld et al., arXiv preprint). https://arxiv.org/abs/2307.12985 (academic; as of 2023-07; authority 7/10; preprint)

[18] Separating measurement and expression models clarifies confusion in single-cell RNA sequencing analysis (Sarkar and Stephens, bioRxiv preprint of the Nat Genet 2021 paper). https://www.biorxiv.org/content/10.1101/2020.04.07.030007v1.full (academic; as of 2020-04; authority 8/10; preprint)

[19] Gene expression distribution deconvolution in single-cell RNA sequencing (Wang et al., PNAS). https://doi.org/10.1073/pnas.1721085115 (academic; as of 2018-06; authority 9/10)

[20] Molecular Cross-Validation for Single-Cell RNA-seq (Batson, Royer and Webber, bioRxiv preprint). https://www.biorxiv.org/content/10.1101/786269v1 (academic; as of 2019-09; authority 7/10; preprint)

[21] Evaluating imputation methods for accurate estimation of cell population fractions in single-cell RNA sequencing (Valyaeva et al., NAR Genom Bioinform). https://doi.org/10.1093/nargab/lqaf204 (academic; as of 2026-01; authority 7/10)

[22] A systematic evaluation of single-cell RNA-sequencing imputation methods (Hou et al., Genome Biol). https://doi.org/10.1186/s13059-020-02132-x (academic; as of 2020-08; authority 9/10; not cited)

[23] False signals induced by single-cell imputation (Andrews and Hemberg, F1000Research, version 2). https://doi.org/10.12688/f1000research.16613.2 (academic; as of 2019-03; authority 7/10)

[24] Deep generative modeling for single-cell transcriptomics (Lopez et al., Nat Methods). https://doi.org/10.1038/s41592-018-0229-2 (academic; as of 2018-11; authority 9/10; not cited)

[25] Quantifying the tradeoff between sequencing depth and cell number in single-cell RNA-seq (Svensson, da Veiga Beltrame and Pachter, bioRxiv preprint). https://www.biorxiv.org/content/10.1101/762773v1 (academic; as of 2019-09; authority 6/10; preprint; not cited)

[26] Normalization by distributional resampling of high throughput single-cell RNA-sequencing data (Brown et al., Bioinformatics). https://doi.org/10.1093/bioinformatics/btab450 (academic; as of 2021-11; authority 7/10)

[27] bayNorm: Bayesian gene expression recovery, imputation and normalization for single-cell RNA-sequencing data (Tang et al., Bioinformatics). https://doi.org/10.1093/bioinformatics/btz726 (academic; as of 2020-02; authority 7/10)

[28] scDesign2: a transparent simulator that generates high-fidelity single-cell gene expression count data with gene correlations captured (Sun et al., Genome Biol). https://doi.org/10.1186/s13059-021-02367-2 (academic; as of 2021-05; authority 8/10; not cited)

[29] DepthDiff: Restoring Low-Depth Single-Cell RNA-Seq Signals via Diffusion Denoising (Hou et al., Biology). https://doi.org/10.3390/biology15151223 (academic; as of 2026-07; authority 5/10)

[30] Low dimensionality in gene expression data enables the accurate extraction of transcriptional programs from shallow sequencing (Heimberg et al., Cell Syst). https://doi.org/10.1016/j.cels.2016.04.001 (academic; as of 2016-04; authority 8/10)

[31] Evaluating Single-Cell Perturbation Response Models Is Far from Straightforward (Heidari et al., bioRxiv preprint). https://www.biorxiv.org/content/10.64898/2026.02.14.705879v1 (academic; as of 2026-02; authority 5/10; preprint; not cited)

[32] Advanced Single-Cell Analysis with Bioconductor, chapter "Normalization, redux". https://bioconductor.org/books/release/OSCA.advanced/more-norm.html (official; as of 2026; authority 8/10; not cited)

[33] Advanced Single-Cell Analysis with Bioconductor, chapter "Droplet processing". https://bioconductor.org/books/release/OSCA.advanced/droplet-processing.html (official; as of 2026; authority 8/10; not cited)

[34] DropletUtils downsampleReads reference page (rdrr.io mirror of the Bioconductor package documentation). https://rdrr.io/bioc/DropletUtils/man/downsampleReads.html (official; undated; authority 7/10; third-party mirror)

[35] scuttle downsampleMatrix reference page (rdrr.io mirror of the Bioconductor package documentation). https://rdrr.io/bioc/scuttle/man/downsampleMatrix.html (official; undated; authority 7/10; third-party mirror)

[36] scanpy.pp.downsample_counts API reference. https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.downsample_counts.html (official; undated; authority 7/10)

[37] scvi-tools API reference: scvi.criticism.PosteriorPredictiveCheck. https://docs.scvi-tools.org/en/stable/api/reference/scvi.criticism.PosteriorPredictiveCheck.html (official; undated; authority 7/10; not cited)

[38] scvi-tools API reference: scvi.model.SCVI. https://docs.scvi-tools.org/en/stable/api/reference/scvi.model.SCVI.html (official; undated; authority 7/10; not cited)

[39] Seurat vignette: Using sctransform in Seurat. https://satijalab.org/seurat/articles/sctransform_vignette (official; as of 2026-11; authority 8/10; not cited)

[40] Illumina: Index Hopping, intro and how to minimize it. https://www.illumina.com/techniques/sequencing/ngs-library-prep/multiplexing/index-hopping.html (official; undated; authority 8/10; not cited)

[41] MerQuaCo: a computational tool for quality control in image-based spatial transcriptomics (eLife reviewed preprint; bioRxiv version read in full). https://elifesciences.org/reviewed-preprints/105149 (academic; as of 2025-07; authority 8/10; reviewed preprint)

[42] Protocol optimization improves the performance of multiplexed RNA imaging (Scientific Reports). https://pmc.ncbi.nlm.nih.gov/articles/PMC12399773/ (academic; as of 2025-08; authority 8/10)

[43] Batch correcting single-cell spatial transcriptomics count data with Crescendo improves visualization and detection of spatial gene patterns (Genome Biology). https://pmc.ncbi.nlm.nih.gov/articles/PMC11863647/ (academic; as of 2025-02; authority 8/10)

[44] Spatially resolved cell atlas of the mouse primary motor cortex by MERFISH (Nature). https://pmc.ncbi.nlm.nih.gov/articles/PMC8494645/ (academic; as of 2021-10; authority 9/10)

[45] High-throughput single-cell gene-expression profiling with multiplexed error-robust fluorescence in situ hybridization (PNAS). https://pmc.ncbi.nlm.nih.gov/articles/PMC5047202/ (academic; as of 2016-09; authority 9/10)

[46] Concordance of MERFISH spatial transcriptomics with bulk and single-cell RNA sequencing (Life Science Alliance). https://pmc.ncbi.nlm.nih.gov/articles/PMC9760489/ (academic; as of 2022-12; authority 7/10)

[47] Systematic benchmarking of imaging spatial transcriptomics platforms in FFPE tissues (preprint version read; journal version Nature Communications 2025). https://pmc.ncbi.nlm.nih.gov/articles/PMC10723440 (academic; as of 2023-12; authority 7/10; preprint)

[48] Standardized metrics for assessment and reproducibility of imaging-based spatial transcriptomics datasets (Nature Biotechnology). https://pmc.ncbi.nlm.nih.gov/articles/PMC13368583/ (academic; as of 2025-12; authority 9/10)

[49] Optimizing Xenium In Situ data utility by quality assessment and best-practice analysis workflows (Nature Methods). https://pmc.ncbi.nlm.nih.gov/articles/PMC11978515/ (academic; as of 2025-03; authority 9/10)

[50] Robust decomposition of cell type mixtures in spatial transcriptomics (Nature Biotechnology). https://pmc.ncbi.nlm.nih.gov/articles/PMC8606190/ (academic; as of 2021-02; authority 9/10)

[51] Gene count normalization in single-cell imaging-based spatially resolved transcriptomics (Genome Biology). https://pmc.ncbi.nlm.nih.gov/articles/PMC11167774/ (academic; as of 2024-06; authority 8/10; not cited)

[52] FOV QC from single-cell gene expression in spatial dataset (Bruker Spatial Biology, CosMx Analysis Scratch Space). https://nanostring-biostats.github.io/CosMx-Analysis-Scratch-Space/posts/fov-qc (official; as of 2024-05; authority 6/10; not cited)

[53] QC and normalization of RNA data (Bruker Spatial Biology, CosMx Analysis Scratch Space). https://nanostring-biostats.github.io/CosMx-Analysis-Scratch-Space/posts/normalization (official; as of 2024-06; authority 6/10)

[54] Xenium Prime add-on panel design report for the Human 5K Pan Tissue and Pathways Panel (10x Genomics). https://cf.10xgenomics.com/supp/xenium/analysis-workshop/Xenium_Custom_Panel_Design_demo_summary.html (official; undated; authority 7/10)

[55] Xenium panel design concepts and terms (10x Genomics support documentation; search extract only). https://www.10xgenomics.com/support/software/xenium-panel-designer/latest/analysis/introduction/xenium-panel-design-concepts-terms (official; undated; authority 6/10; search extract only; not cited)

[56] Xenium Prime 5K data highlights technical note CG000775 Rev B (10x Genomics; search extract only). https://cdn.10xgenomics.com/image/upload/v1754601291/support-documents/CG000775_Prime5K_DataHighlightsTN_RevB.pdf (official; as of 2025-07; authority 6/10; search extract only; not cited)

[57] Biomarker Quantification in Breast Cancer using Xenium In Situ (bioRxiv preprint by vendor authors; abstract only). https://www.biorxiv.org/content/10.64898/2025.12.08.692193v1 (academic; as of 2025-12; authority 5/10; abstract only; preprint)

[58] Multiplexed imaging of high-density libraries of RNAs with MERFISH and expansion microscopy (Scientific Reports; abstract only). https://doi.org/10.1038/s41598-018-22297-7 (academic; as of 2018-03; authority 8/10; abstract only; not cited)

[59] Impact and correction of segmentation errors in spatial transcriptomics (Nature Genetics; abstract only). https://consensus.app/papers/details/w:082f8b81-8e24-4619-b538-33f954681326/?utm_source=claude_code (academic; as of 2026; authority 8/10; abstract only; index record page)

[60] A technical comparison of spatial transcriptomics platforms across six cancer types (Genome Biology; abstract only). https://consensus.app/papers/details/w:819eb4bb-eaa3-4e55-96cf-3c2c7f8a61e5/?utm_source=claude_code (academic; as of 2026; authority 7/10; abstract only; index record page; not cited)

[61] ResolVI: addressing noise and bias in spatial transcriptomics (Nature Methods; abstract only). https://consensus.app/papers/details/w:db5fb85f-1b07-4309-81b1-fddde8433957/?utm_source=claude_code (academic; as of 2026; authority 8/10; abstract only; index record page)

[62] A Kernel Score Perspective on Forecast Disagreement and the Linear Pool. https://arxiv.org/abs/2412.09430 (academic; as of 2026-05; authority 8/10; preprint)

[63] Efficient pooling of predictions via kernel embeddings. https://arxiv.org/html/2411.16246v1 (academic; as of 2024-11; authority 8/10; preprint)

[64] scipy.stats.energy_distance (SciPy v1.18.0 Manual). https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.energy_distance.html (official; as of 2026; authority 8/10)

[65] Robust Forecasting. https://arxiv.org/pdf/2011.03153 (academic; as of 2020-12; authority 8/10; arXiv record)

[66] Learning in an Uncertain World: Representing Ambiguity Through Multiple Hypotheses. https://arxiv.org/abs/1612.00197 (academic; as of 2017-08; authority 9/10; arXiv record)

[67] Winner-takes-all learners are geometry-aware conditional density estimators. https://arxiv.org/abs/2406.04706 (academic; as of 2024-06; authority 9/10; arXiv record)

[68] Stochastic Multiple Choice Learning for Training Diverse Deep Ensembles. https://arxiv.org/pdf/1606.07839 (academic; as of 2016; authority 9/10; arXiv record)

[69] Confident Multiple Choice Learning. https://arxiv.org/abs/1706.03475 (academic; as of 2017-09; authority 8/10; arXiv record)

[70] Diverse M-Best Solutions in Markov Random Fields. https://consensus.app/papers/details/w:4db65158-d869-4f48-a1c8-a399e6c0e36d/?utm_source=claude_code (academic; as of 2012; authority 8/10; index record page; not cited)

[71] Analyzing the Variety Loss in the Context of Probabilistic Trajectory Prediction. https://consensus.app/papers/details/w:c2ab871e-fda7-499b-8e9a-c992e5b35b14/?utm_source=claude_code (academic; as of 2019; authority 8/10; index record page; not cited)

[72] Principal Points and Self-Consistent Points of Elliptical Distributions. https://consensus.app/papers/details/w:bd85b678-082f-48b3-8838-44fb91ca680c/?utm_source=claude_code (academic; as of 1995; authority 9/10; index record page; not cited)

[73] Min-max-min robust combinatorial optimization. https://consensus.app/papers/details/w:1155f390-fdf6-4adc-83f8-29961cca6e34/?utm_source=claude_code (academic; as of 2016; authority 9/10; index record page; not cited)

[74] An analysis of approximations for maximizing submodular set functions-II. https://consensus.app/papers/details/w:ff962a6e-740f-4b29-99d9-10c0296e4cbd/?utm_source=claude_code (academic; as of 1978; authority 9/10; index record page; not cited)

[75] Clustering to minimize the maximum intercluster distance. https://consensus.app/papers/details/w:f6a31d0a-430c-4041-8d49-ab8606d8b356/?utm_source=claude_code (academic; as of 1985; authority 9/10; index record page; not cited)

[76] apricot documentation: Facility Location. https://apricot-select.readthedocs.io/en/latest/functions/facilityLocation.html (official; undated; authority 6/10; not cited)

[77] On the Consistency of Top-k Surrogate Losses. https://consensus.app/papers/details/w:a9577c95-f18a-4e3a-9ea0-a5dce3296aed/?utm_source=claude_code (academic; as of 2019; authority 8/10; index record page)

[78] Winning a Kaggle Competition in Python (course lesson, Modeling chapter). https://campus.datacamp.com/courses/winning-a-kaggle-competition-in-python/modeling?ex=12 (secondary-industry; undated; authority 6/10)

[79] How I Made It To The Kaggle Leaderboard. https://wandb.ai/lavanyashukla/kaggle-feature-encoding/reports/How-I-Made-It-To-The-Kaggle-Leaderboard--Vmlldzo2NzQyNQ (secondary-industry; as of 2020-04; authority 5/10)

[80] A Meta-Analysis of Overfitting in Machine Learning. https://papers.nips.cc/paper/2019/hash/ee39e503b6bedf0c98c388b7e8589aca-Abstract.html (academic; as of 2019; authority 9/10; abstract and published reviews only)

[81] Reviews: A Meta-Analysis of Overfitting in Machine Learning (NeurIPS 2019 peer reviews). https://papers.nips.cc/paper_files/paper/2019/file/ee39e503b6bedf0c98c388b7e8589aca-Reviews.html (academic; as of 2019; authority 6/10; peer-review page; not cited)

[82] On Over-fitting in Model Selection and Subsequent Selection Bias in Performance Evaluation. https://consensus.app/papers/details/w:5e54ec70-38b9-4b9d-a0ce-15446f54929d/?utm_source=claude_code (academic; as of 2010; authority 9/10; index record page)

[83] To combine or not to combine: selecting among forecasts and their combinations. https://consensus.app/papers/details/w:bc45bde4-2425-4d9c-a259-b927c9dd4354/?utm_source=claude_code (academic; as of 2004; authority 8/10; index record page; not cited)

[84] Assessment of Protein Complex Predictions in CASP16: Are We Making Progress? https://consensus.app/papers/details/w:45a7a882-22d3-414a-bbfa-68858b76ef0f/?utm_source=claude_code (academic; as of 2025; authority 8/10; index record page; not cited)

[85] AlphaFold3 at CASP16. https://consensus.app/papers/details/w:c5cd6572-963b-4223-a1be-7ce2e57bd31c/?utm_source=claude_code (academic; as of 2025; authority 7/10; index record page)

[86] Probability Elicitation, Scoring Rules, and Competition Among Forecasters. https://consensus.app/papers/details/w:0d3b9e20-5067-403a-b2fa-041393f183b1/?utm_source=claude_code (academic; as of 2007; authority 9/10; index record page)

[87] Hedging and Approximate Truthfulness in Traditional Forecasting Competitions. https://consensus.app/papers/details/w:b18f450b-160d-4817-b2ff-ec1c99acd087/?utm_source=claude_code (academic; as of 2025; authority 8/10; index record page)

[88] Results and Retrospective Analysis of the CODS 2025 AssetOpsBench Challenge. https://consensus.app/papers/details/w:503bd57a-f1e1-424b-9a07-9e2029ccbbb9/?utm_source=claude_code (academic; as of 2026; authority 6/10; index record page; not cited)

[89] The Ladder: A Reliable Leaderboard for Machine Learning Competitions. https://consensus.app/papers/details/w:84f3f555-e20e-4d51-b319-65840da6b022/?utm_source=claude_code (academic; as of 2015; authority 9/10; index record page; not cited)

[90] Approximation Guarantees for Min-Max-Min Robust Optimization and k-Adaptability Under Objective Uncertainty. https://consensus.app/papers/details/w:4c5361b0-5db2-4135-bd40-b286c9b3722f/?utm_source=claude_code (academic; as of 2024; authority 9/10; index record page; not cited)

## Dropped

- Open Problems in Single-Cell Analysis: denoising benchmark results. https://openproblems.bio/results/denoising/ (task-a, original number 41; community; as of 2026; authority 6/10). Withheld as a precaution: whether the platform has any tie to the challenge or its organisers could not be checked. No finding cites it.

## Counter-review

Two counter-reviews were applied: claims (C: 19 issues, 4 high and 15 medium, plus 10 small points) and applicability (A: 16 issues, 4 high, 9 medium, 3 low). Neither raised a critical issue. Every high and medium issue is applied, some only in part.

- **High.** C1, A1, A4: the headline gain is gross, by scenario and weighting, against slot one, at Low confidence, after a baseline check. C2: the cause of the MERFISH thinning failure is "not established". C3: the pair rule is an unpublished derivation. C4: "all three metric families" is removed. A2, A3: recipe 1 is one frozen complement table with two ablations, a new kill line and a sign-off.
- **Medium.** C5, C15: search results read "the notes found no". C6 to C8: the 100-cell figure is deleted; cell types and source 42 are scoped. C9, C10, C16: depth evidence is re-cited. C11, C12: pooling and practitioner claims are reworded. C13, C14, A6: "sample-specific", both log bases, estimator noise. C17, A9, A10: variogram points labelled by base; "+5 at equal odds" deleted. C18, C19, A7, A11, A12: recipes reworked; unsupported gain ranges deleted. A5, A8, A13: table rows. The ranking follows A.
- **In part.** Low and small points are applied; items needing the brief, the code, the experiment record or new measurements (C14, C17, A1, A6, A10, A12, A16, and which comparisons the two Spearman values are) are stated as open or written as steps; none was run. C and A disagree on draft recipe 7's kill line; C is followed. The registry's "most cited" line omits source 63 and is left unedited.
- **Editor's own.** Confidence in the size of a hidden-target gain is lowered to Low. One draft passage is cut and source 20 cited for method only, after the claims review's compliance flag. Eleven draft sources are no longer cited (6, 8, 22, 24, 28, 32, 39, 73, 74, 76, 90) and two are new (30, 80). Opposing interpretations moved into Key controversies.

## Removal log (classes only)

Added by the compliance auditor (Claude, AI) on 2026-11-12; this section and the next are not part of the edited report. I read this report only: not the research brief, the notes, the raw material, the experiment record or any source.

- Wording of a draft passage about an outside source, cut upstream on a compliance flag but still quoted in the counter-review record: 1, reworded. The record that a cut was made stays. Closest class: outside measurements under the barred-window rule; the report does not let me confirm that.
- Property of a hidden target (test-board limits, anchors, truth counts, any other): 0.
- Per-board leaderboard figures: 0.
- Gene-level knockout findings from outside the released data: 0.
- Other outside measurements from the barred embryonic-day windows: 0.
- Material on the challenge's hidden targets, leaderboard or other teams: 0.
- Recipe steps that send the reader to such material: 0.

Total: 1 passage reworded, none deleted. Everything outside that one line of the counter-review record is unchanged.

## Screening notes

Kept after doubt, with the reason.

1. Scores and gains on targets (the render's gains and losses, the zero-gap scores, the two-library pair and its gap). Kept as our own results. The report says its gains are estimated on practice data and never measured on a hidden target, and it proposes re-scoring the two-library pair, which only a released pair allows. I could not check any of this against the research brief or the experiment record.
2. Counts per cell of the two libraries in that pair, and of the imaging samples. Kept as our own measurements of released data, on the same reasoning.
3. The "ceiling" for a shallow specialist, and the other upper bounds. Kept: each is arithmetic on our own test scores and is stated as conditional. None is a limit read from a test board, and none says where the hidden library sits.
4. The three test depths. Kept as settings of our own thinning test. The report tells the reader not to take them as odds for the hidden library.
5. The scorer's transform, written with a symbolic floor. Kept: it describes the scorer, and no floor value for any board is given.
6. Statistics on quarter-day changes in run profile. Kept as our own findings from released samples. No stage is named, and the barred-window rule covers outside measurements only.
7. Scenario labels and weights, and mentions of the hidden target, the hidden library and the unseen run. Kept: the labels describe practice targets and hypotheses, the mentions are generic or conditional, and nothing states which scenario or which property holds for a hidden target.
8. "Legal" lines naming closed levers, and the descriptions of what sits in each slot. Kept: they state team rules and our own builds, with no limit, count or target value.
9. Literature figures. Kept: none is attributed to mouse embryos. Where the report names the material, it is blood cells, spike-ins, brain sections, cultured cells or tissue cores. Where it does not, I relied on the report's own "Data rules" statement, which I could not check without opening the sources.
10. Source 3 beside the number of genes in the FISH checks (recipe 7). Kept: it gives the size of a validation, not a measurement at any stage. From memory, not checked here, those checks were made on a cultured cell line. It does sit oddly with the report's statement that source 3 is cited for method only.
11. Steps that say to read the experiment record or the builder. Kept: they point to the team's own records and code, not to hidden-target, leaderboard or other-team material.
12. Recipe 6's "Do" line. Kept: the report rejects the recipe and already marks its tables as an external input needing disclosure; they are not in a removed class.
13. The "Dropped" entry and its link. Kept: it reports nothing from the page and is not a recipe step. Its own caution stands: the page's tie to the challenge is unchecked.
14. Other competitions cited in section 4 and named in the references. Kept: none is an embryo-forecasting challenge, and nothing about this challenge's leaderboard or other teams appears.

The Verification section's length figures describe the report before these two sections were added.
