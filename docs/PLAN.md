# Plan

Version 1.2, 2026-10-19, one day before the final phase opens. Update it weekly and log decisions at the bottom.
Evidence behind every claim is in `docs/research/` (one report per topic; `digest_2026-10-19.md` is the screened
summary of the verified research sweep and its gap-fills) and `experiments/*/REPORT.md` (experiments on released data,
each with an independent review); risks are tracked in `ASSUMPTIONS.md`.

## 1. What we are playing for

**Scoring.** Overall standing = T1 + mean(T2 embryo, T2 heart) + T3. A point on T1 or T3 is worth two points on a T2
board. Every metric becomes a skill `1/(1+x)`, where x is the fraction of the floor-to-ceiling gap still left: floor
0.5, halfway 0.667, a quarter left 0.8, ceiling 1.0. Board scores are 100 × the weighted sum of skills:

| Task | Formula (skills in [0, 1]) |
|---|---|
| T1 | 0.25 de_score + 0.25 de_direction + 0.30 mmd_u + 0.20 variogram |
| T2 (per board) | 0.25 mean(de_score, de_direction) + 0.25 (0.6 mmd_u + 0.4 variogram) + 0.25 mean(d2_shape, occupancy_dice, abs scale_log_ratio) × min(1, 2·skill(neighborhood_mmd)) + 0.25 neighborhood_mmd |
| T3 | 0.30 de_score + 0.25 de_direction + 0.25 abs(severity_slope) + 0.20 (0.6 mmd_u + 0.4 variogram) |

All three were re-derived from 1,000+ public leaderboard rows to within 1e-4, including the T2 structure gate. The
local emulator (`vec eval`) implements exactly this.

**The field.** 313 Human teams. On the validation boards the Human prize line (6th) is 220.8 and the leader 224.7.
Those numbers are inflated: teams had 8 scored submissions a day for two months, and the size and severity columns are
saturated by search. The test allows 2 blind submissions per board on harder targets (E12.5 two days out, the
"hardest" knockout), so test totals will drop.

**Target.** About 190+ on the test boards (T1 60-65, T2 62-68 averaged over its two boards, T3 58-63), which should be
competitive for a placed prize. The differentiators are the high-variance blind calls: new cell populations in T1 and
the direction of the T3 response. T2 embryo interpolation is a board everyone can do reasonably well.

## 2. Non-negotiables

Full detail in `docs/DATA_POLICY.md`. In short:

- Nothing measured in a held-out window, by any route. T1 and T2 heart: no external data after E9.5 up to and including
  E13.5; after E13.5 is fine with disclosure. T2 embryo: nothing strictly between E7.25 and E8.0. T2 heart
  interpolation: nothing strictly between E8.25 and E8.75. T3: no beta-catenin or Gata4 perturbation data, no other
  alleles at comparable stages, no phenocopies (GSE165300, Mesp1-Cre; Wls, is excluded).
- **No gene-level knockout findings from papers for T3** (organiser ruling, 2026-10-15). Broad phenotype background is
  allowed: beta-catenin loss in cardiac mesoderm impairs second-heart-field expansion and outflow-tract/right-ventricle
  growth.
- No contaminated pretrained models. Mouse single-cell foundation models trained on CELLxGENE or GEO corpora, and the
  Navigo checkpoints, contain in-window cells. We use none.
- **The leaked E12.5 heart size never enters a model.** The heart size comes from the pre-registered procedure in
  section 4.3.
- Scores choose between candidates; they never compute one.
- Every external source is logged in `EXTERNAL_DATA.md` before use, with the removal code committed.
- No scored upload without the user's explicit go-ahead.

## 3. Design principles for every board

These come from the scorer code, controlled probes and the public record of what worked and failed on the
validation boards.

1. **Real cells, count-space edits.** The variogram on 32,285 genes is mostly a per-gene sparsity test. Decoded means
   make it 7x worse; adding a dense delta to every cell doubles it (that is why `pseudobulk_shift` and `shift_transfer`
   sit below the floor). Recover integer counts (`c = expm1(x)·s/1e4`, `s = 1e4/expm1(min nonzero)`), edit by
   binomial thinning or Poisson top-up, renormalise to 1e4, then log1p. Sample counts if a decoder is ever used.
2. **Decouple the pseudobulk from the population.** de_score, de_direction and severity see only the pseudobulk. mmd_u
   and the variogram see the cells. Choose the cells (composition, real states) for the distributional metrics, then
   graft the intended pseudobulk change onto them without destroying sparsity.
3. **Pseudobulk hygiene.** Build deltas against the full released reference stage, not a subsample. Set the change to
   exactly 0 for genes we know nothing about and add no noise. Keep exact zeros for undetected genes. Submit at least
   2,000 cells (the MMD draw) and never fewer than 1,500 (the variogram's variance penalty). Predict decreases only for
   genes that can be called down (reference pseudobulk >= 0.25). de_direction rewards sharing the server's 10%
   reference-subsample noise, so a plain copy already scores above zero: never read de_direction > 0 as evidence that a
   direction is right, and make the local harness imitate the 10%/split-half protocol. Hold sex-linked (Xist, Tsix,
   Y-linked) and strain-polymorphic genes flat: the pooled T1 stages differ in sex ratio. Score T3 against each
   knockout replicate separately (the board ranks one) and optimise the worse or the mean.
4. **Composition is the main lever** on all four boards: which cell types, and how many.
5. **Interpolate whenever the rules allow.** E12.5 can be bracketed between our E10.5 and legal external hearts after
   E13.5, for both T1 and T2 heart. E7.75 is bracketed by E7.5 and E8.0 by design.
6. **Spatial hygiene** (`experiments/density/REPORT.md`, `experiments/coords/REPORT.md`).
   - **Cell count: never below truth_cells; otherwise as many as allowed, up to max_cells.** Below 1x every seed loses
     (0.5x costs 1-5 points, 0.25x 2-10), because occupancy_dice thins the truth to the submitted count. Above 1x a
     near-perfect prediction is flat, but imperfect ones keep gaining (intermediate quality +1.6 points from 1x to 8x;
     +1.0 at 3x on an 8k-cell truth half). Keying the count to max_cells also avoids choosing it from truth_cells.
     One community observation on the live heart-extrapolation board points the other way (10k and 17.6k cells lost
     0.5-0.7 points against 5k), so re-test truth_cells, 2x and max_cells on the E10.5 and E7.5 validation targets
     before fixing the test counts, with this rule fixed in advance: keep max_cells unless, on a validation board, a
     lower max_cells level beats it for the intermediate-quality condition by more than 0.3 points (paired mean) in at
     least 8 of 10 seeds. The evidence for max_cells reaches 8x truth_cells (3x on an 8k-cell truth half). The test
     boards may publish only their cell limits.
   - **Random rows only.** Rows are stored section by section. A contiguous block costs a near-perfect prediction 16-28
     points (mid-file) or 21-48 (first rows) at 1x, and still 12-21 at 8x.
   - **Jitter budget.** Up to 0.25 x the median nearest-neighbour distance is free (use it to break duplicate
     coordinates); 2 NN costs 1.2-1.8 points at 1x. One NN is about 20 units, 6-9% of the RMS radius, at 5% sampling.
   - **Proper rotations only.** All released heart frames share one handedness (Pitx2 landmarks; E8.0 unresolved). Under
     the pinned (fixed) scorer a mirrored copy of a correctly shaped cloud loses up to 0.12 Dice, about the whole
     floor-to-ceiling gap. Check E7.5 against E8.0 by proper vs improper shape alignment before building the embryo
     board.
   - **Scorer build.** The live board behaves like the pre-fix veckit build: a mirror twin scores the same as the
     original, and Dice depends on a parity coin set by the SVD sign pattern (row order), worth up to ~3.5 heart-board
     points at good shape quality. The emulator models both builds (`experiments/scorer_build/`): `build="main"`
     (the default) compares candidates on shape quality without the parity coin, and `build="prefix"` reproduces live
     scores (calibration and the anchor check). A "parity hedge"
     across the two submissions would exploit this defect: per the user's decision (2026-10-19) we ask the organisers
     first and use it only if they say it is acceptable or confirm the build.
   - **Native coordinates.** xy are slide micrometres and the stored z step is most likely the section pitch; never
     rescale z across stages. A 1.5x z mismatch costs ~12.6 heart-board points on shape for an otherwise perfect
     prediction; neighbourhoods barely notice.
   - **Duplicated records.** E8.0, E8.25_late, E8.75 and the Mab21l2 KO hold identical cells in adjacent sections (KO
     5.6%). Drop them (`vec.prep.load_stage(drop_duplicates=True)`) before neighbourhood, pseudobulk or composition
     work on prediction inputs; the emulator loads targets as released, since the server's truth includes them.
   - Expression moves with its cell.
7. **Pre-register.** Size models, composition scenarios and both submission candidates per board are fixed and hashed in
   git before the first test submission.

## 4. Board plans

Effort split: T1 ~35%, T3 ~25%, T2 heart ~15%, T2 embryo ~10-15%, shared core and harness ~10-15%.

### 4.1 T1: E12.5 single-cell population (weight 1)

Facts that shape it (`docs/research/digest_2026-10-19.md`, gap-fills 1 and 3):
- E12.5 is "only the heart chamber" (organisers' dissection note; a few oesophagus-junction cells), while E8.5-E10.5
  include surrounding tissue. The scorer counts every cell.
- **Scope is the biggest lever.** In a pseudo-task through the pinned scorer, re-weighting real reference-stage cells to
  a chamber-only composition, with no state change, moved T1 from 51 to 75 (+12 to +26 over 11 truth variants). Drop
  the surrounding tissue and also the second-heart-field, pericardial and pharyngeal types: keeping them at the
  reference share cost 13 points, dropping them when the truth keeps 15% cost 1-3. Blood about 1.5-2% (over-predicting
  by 10 points costs 7-9). The cardiomyocyte share of heart cells carries ~75% of the composition risk, and legal E14.5
  data disagree on it (40%, 42% and 77% across datasets and strains). Order the work: scope (~+20), per-type state
  shifts (~15), births (~1.5-5).
- **Direction.** Within-platform fold changes from Feng 2022 chamber hearts (legal stages only) beat internal momentum in
  a test-like geometry (sign agreement 0.92 vs 0.73, Spearman 0.80 vs 0.40). JAX nuclei are run- and depth-confounded
  and are dropped as a direction source. The pre-registered readout in `experiments/t1_direction/` decides on the
  E9.5 -> E10.5 transition whether M2 is built: the external direction must beat cross-fitted momentum by >= 0.05 on
  both DE metrics, with cardiomyocyte sign agreement >= 0.70.
- **Legal external data** (each registered in `EXTERNAL_DATA.md` with committed filters before use): Feng 2022
  (GSE193346; E9.5, CD1 only, and E14.5 and later, filtered on the per-cell stage column), Feulner 2026 whole hearts
  (E14.5-P7), JAX/Qiu 2024 embryos 1-22 and 54-74 only (23-24 and 50-53 come from boundary or in-window litters;
  1-12 are T1-only), ENCODE heart bulk from E14.5. Use within-dataset fold changes, never absolute external levels; map
  genes by Ensembl ID.
- Most T1 genes (8,522 Gm*, 2,524 *Rik) have no external signal: hold them at E10.5 levels.
- **Emergent populations, contingent.** Which classes appear by E12.5 can rest on qualitative developmental background
  inside the window only if the organisers allow it (asked 2026-10-18, unanswered). Until then use the legal route:
  classes present in the legal E14.5+ chamber data but absent at E10.5, with E12.5 shares interpolated between the
  two. Births are worth only ~1.5-5 points.

| Step | What | Effort | Kept if |
|---|---|---|---|
| B0 floor-plus | E10.5 real cells re-weighted to a chamber-only composition (surrounding, second-heart-field, pericardial and pharyngeal types dropped; blood ~2%) plus a per-type count-space momentum shift | 2 days | beats copy_last on the E10.5 backtest |
| M1 main | B0's composition with the cardiomyocyte share as the hedged parameter, per-type count-space shifts with direction from the Feng bracket (E10.5 to E14.5+, t ≈ 0.5) if the readout passes, damped momentum otherwise, and cheap proxy births by the legal route | 1-2 weeks | beats B0 on the backtest |
| M2 diverse bet | moscot OT bracketing from our E10.5 to external E14.5 in a shared scVI latent, decoded by count sampling; built only if the readout passes | 2-3 weeks | beats M1 on at least 2 of 3 pseudo-tasks |
| Stretch | Whole-ladder non-autonomous flow with growth (VGFM-style, reimplemented on torchcfm/POT) | only if time | beats M1 |

Validation: predict E10.5 from E8.5 and E9.5 (plus legal external data), scored by the emulator over 10 seeds. Bracket
fraction t ≈ 0.2-0.235 there versus ≈ 0.5 at test, so prefer parameters with a physical meaning (time fractions,
shrinkage) over free knobs. After the release, re-run the composition pseudo-task with E10.5 as the reference and
commit the s1/s2 rule before any upload.

Submissions: s1 = M1 with cardiomyocytes at about 58% of heart cells; s2 = identical except about 70% (the hedge
belongs on the cardiomyocyte share; expected best-of-two gain ~1.1-1.4 points). Both shares are refit on E10.5 and the
legal E14.5+ data and committed before any upload; M2 replaces s2 only if it wins the backtests.

### 4.2 T2 embryo: E7.75 interpolation (weight 0.5)

Pipeline:
1. Canonicalise E7.5 (released at P3) and E8.0 (centre, PCA rotation with det +1, divide by RMS).
2. Couple them with fused Gromov-Wasserstein (moscot `SpatioTemporalProblem` or POT; alpha 0.3-0.6, mild unbalancing).
3. Proper-rotation Procrustes under the coupling, then McCann interpolation of positions at t ≈ 0.5. Calibrate a
   developmental clock on the E7.5 pseudo-target (nominal t = 1/3) and carry the warp over.
4. Expression as a mixture of coupled real cells (each sampled cell keeps a real profile; count-space only).
5. Size from the pre-registered embryo procedure (`experiments/coords/size_procedure.py`): PCHIP on raw 3D log RMS
   through E6.75, E7.25, E7.5 and E8.0, with leave-one-out errors at E7.25 and E7.5 reported next to the prediction.
   The embryo curve is convex (two-point interpolation of E7.25 overshoots by 0.21 in log), so E7.5 matters. If E7.5's
   z step or sectioning orientation differs from E8.0's, record the z-convention risk with the submission (up to about
   0.25 in log if E7.75 were cut along its long axis like E6.75 and E7.25).
6. Cell count: max_cells (never below truth_cells); optional blue-noise relaxation if pseudo-targets show no loss on
   d2 or neighborhood_mmd.

Pseudo-targets: E7.5 from (E7.25, E8.0) and E7.25 from (E6.75, E7.5). Target emulator score of at least 70.

External (optional, legal only outside (E7.25, E8.0)): Pijuan-Sala/extended gastrulation atlas E6.5-E7.25 and
E8.0+; Sampath Kumar 2023 3D Slide-seq at E8.0+. Never the Chi-lab GSE176306 E7.5/E7.75 samples.

Submissions: s1 central clock; s2 alternative clock or the nearest-stage variant (E8.0 tissue shrunk to the predicted
size, mixture expression).

### 4.3 T2 heart: E12.5 extrapolation (weight 0.5)

Facts (`experiments/coords/REPORT.md`): xy are slide units (rigid copies of the slide coordinates, scale 1.0000, in
29 of 35 KO sections); z is the section position times a per-file step (10, 14 or 15), consistent with xy units in
E9.5 and the KO but not clearly so at E8.0-E8.75. Every file has 30-40 sections, so the z extent is set by acquisition,
and the stack cuts the heart in E9.5 and the KO. Raw RMS radii are not comparable across stages (E8.25_late includes
yolk sac and non-cardiac cells; E8.75 is mostly heart). Label vocabularies differ by file: only CM transfers cleanly
between E8.75 and E9.5 (`experiments/harmonise/REPORT.md`). Every geometry deformation tried publicly on extrapolation
lost points; rescaling helped.

**Pre-registered size procedure (fixed 2026-10-19, before the E10.5 release):** `experiments/coords/size_procedure.py`,
run by `experiments/coords/run_p3.sh`.
1. Core = the CM and Endocardium groups of `configs/lineage_map.yaml` at the commit (cushion endocardium in;
   epicardium, proepicardium and septum transversum out). New E10.5 labels are mapped by the rule in that file's header
   and committed before the run; the procedure refuses unmapped labels. Duplicated records are dropped first.
2. Growth uses only rate intervals that start at or after E9.0, which leaves E9.5 to E10.5. Rates measured over
   E8.0-E8.75 cannot predict later growth: the core RMS is flat from E8.0 to E8.75 and then nearly doubles at E9.5, so
   every rule scored like copy_last on the E9.5 pseudo-target. The planned pseudo-targets (E9.5 and E10.5 from stages up
   to E8.75) are dropped for the same reason; no released pseudo-target can rank rules in the E10.5-E12.5 regime.
3. E12.5 core log growth = the log-space average of three rules on the E9.5-E10.5 interval: latest rate, half rate and
   a power law in t, i.e. 1.58 x the observed log growth. The justification uses no assumed rate: if the mean rate over
   E10.5-E12.5 is between 0.5 and 1 times the E9.5-E10.5 rate, the average has the smallest worst-case bias of the
   candidates (0.58 x the observed log growth), at the cost of more noise amplification than half rate (3.2 vs 2.4
   per-embryo SDs). The table is in the procedure's output.
4. z: grow the in-plane extent and hold z at E10.5's values when the core's thinnest axis is within 30 degrees of z in
   both E9.5 and E10.5 (the stack, not the tissue, sets z); otherwise grow in 3D. All four z modes are logged next to
   the pre-registered choice.
5. Scope: three hypotheses for the full-cloud/core ratio (H1 trend, H0 hold, H2 half); s1 uses H1 and s2 uses H2, each
   rescaled about its centroid to its own full-cloud RMS.
6. No external size bracket (units, dissection scope and stack depth do not transfer). The JSON output, with the git
   hash and input checksums, is committed before any test submission. Nobody on this step opens veckit 0.1.2 or PR #6.

Pipeline:
1. Start from real E10.5 tissue (released at P3), duplicates dropped, sampled at random to max_cells (never below
   truth_cells).
2. Damped composition extrapolation on the merged lineage level (CM, Endothelial, CranialMesoderm, NCC, Endoderm,
   Ectoderm; labels harmonised with `vec.harmonise`), centred log-ratios over E8.75-E10.5, damping chosen on the
   1.75-day pseudo-target (expected 0.3-0.6). CardiacProgenitor and Pericardium have no stable cross-stage meaning and
   stay out of the trends.
3. Emergent populations (contingent, as in section 4.1): the class list comes from the legal route unless the
   organisers allow qualitative in-window background, and so does placement: each new cell goes next to the E10.5 cells
   it most resembles in expression, not by in-window anatomy. Expression comes from E10.5 precursors edited in count
   space, or from external post-E13.5 cells mapped to MERFISH (only if cross-task data use is confirmed).
4. Expression direction per type from the external bracket where available, otherwise damped momentum.
5. Grow coordinates with the pre-registered size and keep E10.5's shape. A shape operation (axis stretch, radial trim,
   type-stratified subsampling) is adopted only if it improves both released transitions (E9.5 -> E10.5 and
   E8.25_late -> E8.5) by >= 0.15 d2 skill without lowering Dice and has a physical rationale; then it is applied at
   half magnitude in at most one submission. The leaders' validation shape gains look transition-specific.
6. Scope factors (radius of a scope-reduced cloud relative to the whole) are re-measured on E10.5; on E9.5, dropping
   all shell lineages gives x0.95 and the cardiac core alone x0.83, and scope moves shape as well as size.

Submissions: two scope hypotheses, H1 continuation (shell fraction follows the E9.5-E10.5 trend) and H2 heart-focused
(shell halved or less), each with its own pre-registered size.

### 4.4 T3: beta-catenin knockout at E8.75 (weight 1)

Facts: Ctnnb1 and the classic Wnt targets are not on the 500-gene panel, so zeroing the gene does nothing. Deletion is
limited to the Mesp1 lineage. The released Gata4 knockout (from P3) has the same stage, driver and wild-type reference
as the test. Every published reference model gets the response direction wrong, partly for mechanical reasons: an
added delta on 93%-sparse data realises only about 4% of an intended decrease.

Day-1 diagnostics (CPU, DGX):
- depth, genes per cell, zero fraction and duplicated records per sample (`vec data qc`);
- z-coverage profiles of each knockout against its wild type before any composition comparison;
- harmonised composition per sample (`vec.harmonise`, knockout mapped into the matched wild type's vocabulary), quoted
  as a top-call/EM bracket with its z-trim range;
- Gata4 replicate 1 vs 2 vs wild type: the replicate difference is the noise bar every composition effect must clear;
- whether knockout files carry littermate wild-type cells;
- an oracle decomposition of the Gata4 truth into composition, depth/coverage and within-type parts, each scored with
  the emulator.

What the Mab21l2 knockout already shows (`experiments/harmonise/REPORT.md`, with its review): group-level differences
from the E9.5 wild type are file-against-file differences, not genotype effects, because the two stacks cover different
anatomy (trimming the stack ends removes or flips most of them). The robust difference is inside the heart: among
cardiomyocytes, V-CM rises from 0.39 to 0.61 and SV-CM falls from 0.16 to 0.04, under every trim, classifier and batch
correction tried. There is also a gene-level technical shift between the knockout and wild-type files (rank correlation
0.32-0.74 across matched, unrelated cell types) that component N must carry.

What the Mab21l2 decomposition shows (gap-fill 2 in `docs/research/digest_2026-10-19.md`; emulator scores, 10 seeds):
the response is mostly **one per-gene factor shared by every cell type**. Estimated from the non-Mesp1 compartments
only (endoderm, surface ectoderm, neural tube, neural crest) and applied to all WT cells, it scored 79.3, against 59.5
for the knockout's full harmonised composition and 66.1 for "real non-lineage knockout cells plus WT lineage cells".
Coverage adds nothing once the factor is present, lineage composition adds +1.5 on top, depth alone drives severity to
its floor and helps only together with the factor, and the type-specific remainder (6.5% of the delta) has chance-level
sign. A blind rule built from the broad phenotype failed (de_direction 0.07 against 0.06 for plain WT resampling;
severity at the floor in 8-9 of 10 seeds). About half of the factor is knockout-specific and half specific to the WT
reference: the WT half transfers to the test board for certain, the knockout half only if the series shares a batch.
Because the factor is measured in compartments the knockout cannot touch, much of it may be sample processing rather
than biology, so its value on the test rests on the beta-catenin knockout sharing processing with the released ones.
The first real check is cross-knockout transfer on release day (rule D11b: the Mab21l2 factor scored on Gata4); if it
fails, the two submissions must not both use the Gata4 factor (rules.yaml note 6), and composition (C) and the pathway
programme (B) stay as hedges rather than afterthoughts.

Components (all legal under the 2026-10-15 ruling):
- **F, per-gene factor (lead component):** estimated from the non-Mesp1 compartments of a released knockout against
  its matched WT and applied in count space to all WT E8.75 cells, with per-type depth only alongside it. The two
  submissions differ in how F is sourced: the full Gata4 factor plus Gata4 depth, or a reference-only factor built from
  the WT series.
- **C, composition (secondary):** rides on top of F in both submissions. The states come from the broad phenotype
  (`configs/t3_scenarios.yaml`; second-heart-field, outflow-tract and right-ventricle-lineage states depleted, first
  heart field and left ventricle spared); the magnitude comes from within-heart Gata4 effects (harmonised, bracketed,
  z-trimmed) and counts only where it exceeds the Gata4 replicate difference. Its sign is hedged only if the
  pre-registered rules say so.
- **B, pathway programme (low priority):** a WT-derived activity score built in the 32k-gene scRNA at cluster level from
  a generic pathway database set whose provenance excludes any perturbation of the held-out genes, regressed onto panel
  genes within type. Needs the cross-task ruling (A15); kept only if the analogous check passes on Gata4.
- **Weights:** leave-one-knockout-out (Mab21l2 to Gata4, Gata4 to Mab21l2, Gata4 rep1 to rep2) with the official
  scorer; calibrate the magnitude so the severity slope on likely DE genes is near 1; log the slope's R² on every fold
  and drop any component that adds a shared offset (it can push severity to the floor).
- **Pre-registration:** `experiments/t3_decomp/` with its `rules.yaml` is committed before the Gata4 download and run
  unchanged on both replicates on release day; its decision file chooses the components.
- **Pending confirmation:** using the released Gata4 files inside the beta-catenin prediction (follow-up question 3 to
  the organisers). If refused, both submissions use the reference-only factor.

Submissions: s1 = F (Gata4-sourced) + C; s2 = F (reference-sourced) + C. Both are built before the first upload. Any
rule that picks or changes s2 after s1's score is written beforehand in terms of component skills and never reads a
target property.

**Provenance rule for T3 gene directions (pre-registered 2026-10-18).**
1. No per-gene direction or magnitude is written by hand. None comes from knockout literature, from curated target
   lists built from beta-catenin or Gata4 knockout or perturbation experiments, or from LLM elicitation.
2. Every per-gene change in a T3 prediction is the output of committed code applied to released data only:
   - the Gata4 and Mab21l2 knockouts against their matched wild types (the per-gene factor from their non-Mesp1
     compartments, and composition effect sizes);
   - wild-type cell-state expression profiles combined with a composition scenario (the composition term);
   - wild-type co-variation with Wnt activity (T1 scRNA, only if cross-task use is confirmed).
3. The composition scenario names only cell states and lineage groups. They come from the broad phenotype
   (second-heart-field, outflow-tract and right-ventricle lineages depleted; first heart field and left ventricle
   relatively spared), with magnitudes measured in the released knockouts. The scenario lives in
   `configs/t3_scenarios.yaml` and is committed before any test submission.
4. Each test submission records the git commit, config hash and data checksums in the submission ledger.
   `METHOD_SUMMARY.md` describes this procedure.
5. The gene-level findings removed from the research notes on 2026-10-18 are never consulted (see `docs/DATA_POLICY.md`).

## 5. Shared architecture (Generality Award option)

One "transport-then-adjust" operator underlies every board:
1. Couple adjacent populations with unbalanced OT (fused Gromov-Wasserstein when positions exist).
2. Move cells (displacement interpolation, or a damped learned velocity with a growth term).
3. Reweight composition.
4. Decode by retrieval from real cells.
5. Add births from the best legal source.
6. Calibrate a few shrinkage parameters by leave-one-transition-out.

Task heads hold what differs: the 32K-gene vs 500-gene encoders, the coordinate head (T2) and the genotype head (T3).
copy_last, pseudobulk_shift, bracket interpolation and wild type plus response are all settings of this operator, so it
should cost nothing per board.

Decision point D3 (2026-11-12): if the shared core is within 1 point of the best specialised variant on every board's
honest validation, claim the award in the method summary; otherwise use specialised code where it wins.

## 6. Infrastructure

Done:
- `vec eval`, the emulator, with the gate formula confirmed against the live board.
- `vec check`, the free format checks; `vec baseline`, the reference baselines; `vec data`, the safe downloader.
- `vec.prep` loaders (exact count recovery, lineage map with a merged cross-stage level), `vec.spatial` (z snapping,
  duplicated-record detection, `vec data qc`), `vec.harmonise` (label transfer across files and stages).
- `vec.boards.update_boards`: board refresh that strips test-board anchors in memory.
- Experiments with independent reviews: density and sampling, coordinates and the pre-registered size procedure,
  label harmonisation (`experiments/*/REPORT.md`).
- DGX environment (torch 2.14 cu130, scvi-tools, moscot), assumption ledger, data policy, external-data register,
  release-day runbook (`docs/P3_RUNBOOK.md`).

Next:
1. Calibrate the emulator against the published reference rows on the validation boards (A2), with the scorer-build
   switch (pre-fix vs main) and the validation anchor check.
2. Pseudo-task harness with leave-one-transition-out splits and paired comparisons over 10 seeds.
3. External-data ingestion with window filters as committed code (Feng by per-cell stage, JAX by embryo ID: 1-22 and
   54-74), registered in `EXTERNAL_DATA.md` before use.
4. Run registry; submission client (feat-016) that defaults to format-only, counts distinct file hashes per board whose
   status is not rejected, treats the server's quota endpoint as the source of truth, keeps the filename and form
   fields identical on retries, and sends an idempotency key on scored uploads.
5. Set `POT_BACKEND_DISABLE_JAX=1` and `XLA_PYTHON_CLIENT_PREALLOCATE=false` before importing POT or JAX on the DGX.

## 7. Timeline and decision points

| Dates | Work | Decision point |
|---|---|---|
| Oct 18-20 | Token and download; organiser questions out; leak report sent; experiments on released data; size procedure and release-day scripts committed; redactions and history rewrite | D0: data on the DGX, else escalate |
| Oct 20-26 | Release-day runbook (`docs/P3_RUNBOOK.md`): board intake with anchors stripped, data, label map, pre-registered runs (size, T3 decomposition, T1 direction readout, harmonisation, density, anchor check); emulator calibration; T3 day-1 diagnostics; B0 on all four boards; external data ingested | D1 (Oct 22): test-board reference stages fixed |
| Oct 27-Nov 2 | T2 embryo FGW pipeline; T1 M1 v1; T3 F+C v1; T2 heart growth v1 | D2 (Oct 31): every board beats its floor on internal validation |
| Nov 3-9 | T1 external directions and births; T3 leave-one-knockout-out weights (pathway programme only if its check passes); T2 heart emergent populations and scope hypotheses; T1 M2 if the readout passed | - |
| Nov 10-16 | Ablations, seed robustness, freeze s1 candidates (Nov 16) | D3 (Nov 12): shared core vs specialised |
| Nov 17-23 | s1 on all four boards by Nov 22 (after user approval); s2 rules pre-registered and hashed | D4: read scores, apply the written rules |
| Nov 24-30 | s2 by Nov 30; method summary and disclosures; final nomination | - |
| Dec 1-2 | Buffer only. Archive every submitted file with its code hash for the 14-day winners' report | - |

## 8. Submission policy

- Format checks are free and unlimited. Scored uploads need the user's explicit go-ahead, and the client defaults to
  format-only.
- s1 = the best honestly validated candidate. s2 = a pre-registered alternative chosen by a written rule. Never submit a
  plain copy, and never in the last hours (the portal has had outages).
- Each submission carries its external-data disclosure text.

## 9. Risks

`ASSUMPTIONS.md` holds 19 scored assumptions with go/no-go tests. The biggest:
- **A11**, the T3 response sign (risk 5; early proxy: a blind Gata4 test);
- **A6**, T1 E12.5 composition;
- **A8**, heart scope and size;
- **A5**, the value of the external bracket;
- **A14**, a metric change from veckit PR #6 (mitigation: aim for slope ≈ 1 on T3 and realistic cell spacing on T2);
- **A15**, the cross-task data ruling;
- **A19**, whether released knockout-vs-WT composition shifts measure genotype effects (Gata4 replicates decide).

## 10. Open questions for the organisers

`docs/ORGANIZER_QUESTIONS.md`. Sent on 2026-10-18 (unanswered): cross-task use of released data, qualitative in-window
background, `truth_cells` for the submitted count, Generality scope, scorer version, deadline time zone. Follow-up,
sent 2026-10-19: scorer build and canonicalisation (and the parity defect), the wrong-sign severity rule, released Gata4
files inside the beta-catenin prediction, duplicated records, the z coordinate, `truth_cells` and tissue size, stage
conventions for external embryos, and an addendum to the leak report.

## Decision log

- 2026-10-18: Human Team track only (user). Python 3.12 uv environment. veckit pinned to git main cc3859a (the live
  board's metric set, not PyPI 0.1.2). No mouse foundation models. T3 built without gene-level knockout literature.
  Leaked heart size quarantined; heart size from the pre-registered procedure only.
- 2026-10-19: three pre-release experiments on released data, each independently reviewed (`experiments/*/REPORT.md`).
  Decisions: (1) cell count = max_cells, never below truth_cells (A13 revised); random rows only; jitter at most 0.25
  NN. (2) Heart size procedure fixed before the E10.5 release: core = CM + endocardium (epicardium is shell), growth
  from E9.5-E10.5 only, fixed average of three rules justified without an assumed rate, z held at E10.5 when the
  sectioning plane matches, three scope hypotheses, no external bracket. (3) Label transfer with `vec.harmonise`;
  cross-file compositions only through it, quoted as a top-call/EM bracket; Mab21l2 group-level shifts are not used as
  genotype effects, the within-CM shift is. (4) Lineage map: merged level for cross-stage work, KO-only overrides for
  LPM and DPM, a pre-registered rule for mapping new labels. (5) Duplicated records dropped from prediction inputs.

- 2026-10-19 (b): integrated the screened research digest (a separate agent screened the raw results; the lead read
  only the digest). T3 now leads with the per-gene factor measured on released knockouts; the blind broad-phenotype
  composition rule is demoted to a secondary term. T1 works in the order scope, then state shifts, then births; the
  Feng bracket is the direction source if its pre-registered readout passes; JAX is dropped as a direction source and
  restricted to embryos 1-22 and 54-74. Emergent classes and their placement use the legal route until the organisers
  answer. Live scorer modelled as the pre-fix build; no parity hedge without the organisers' approval (user). Test-board
  anchors and the public board file are quarantined. History to be rewritten again (user).
- 2026-10-19 (c): release-day scripts committed after screening and independent review: T3 decomposition with
  pre-registered rules (lead's notes: the Mab21l2 factor source is transfer:mab21l2_all; R includes E8.5; if the
  organisers bar the Gata4 files entirely, both submissions use the Mab21l2 factor and no Gata4 score chooses anything;
  replicate independence checked first), T1 direction readout (design signed off: Feng E14.5 alone decides, matched
  coverage, two-fold cross-fitting, cardiomyocyte sign agreement on the full target, no momentum+Feng option; guard
  against non-default configs and label maps), and the scorer-build switch (emulator default `main`, `prefix` for live
  reproduction). Feng 2022, JAX/Qiu 2024 and the vendored veckit 0.1.1 function registered in `EXTERNAL_DATA.md`.
