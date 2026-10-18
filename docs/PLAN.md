# Plan

Version 1.1, 2026-10-19, one day before the final phase opens. Update it weekly and log decisions at the bottom.
Evidence behind every claim is in `docs/research/` (one report per topic) and `experiments/*/REPORT.md` (experiments on
released data, each with an independent review); risks are tracked in `ASSUMPTIONS.md`.

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
   2,000 cells (the MMD draw) and never fewer than 1,500 (the variogram's variance penalty).
4. **Composition is the main lever** on all four boards: which cell types, and how many.
5. **Interpolate whenever the rules allow.** E12.5 can be bracketed between our E10.5 and legal external hearts after
   E13.5, for both T1 and T2 heart. E7.75 is bracketed by E7.5 and E8.0 by design.
6. **Spatial hygiene** (`experiments/density/REPORT.md`, `experiments/coords/REPORT.md`).
   - **Cell count: never below truth_cells; otherwise as many as allowed, up to max_cells.** Below 1x every seed loses
     (0.5x costs 1-5 points, 0.25x 2-10), because occupancy_dice thins the truth to the submitted count. Above 1x a
     near-perfect prediction is flat, but imperfect ones keep gaining (intermediate quality +1.6 points from 1x to 8x;
     +1.0 at 3x on an 8k-cell truth half). Keying the count to max_cells also avoids choosing it from truth_cells.
   - **Random rows only.** Rows are stored section by section. A contiguous block costs a near-perfect prediction 16-28
     points (mid-file) or 21-48 (first rows) at 1x, and still 12-21 at 8x.
   - **Jitter budget.** Up to 0.25 x the median nearest-neighbour distance is free (use it to break duplicate
     coordinates); 2 NN costs 1.2-1.8 points at 1x. One NN is about 20 units, 6-9% of the RMS radius, at 5% sampling.
   - **Proper rotations only.** All released heart frames share one handedness (Pitx2 landmarks; E8.0 unresolved). A
     mirrored copy of a correctly shaped cloud loses up to 0.12 Dice, about the whole floor-to-ceiling gap. Check E7.5
     against E8.0 by proper vs improper shape alignment before building the embryo board.
   - **Duplicated records.** E8.0, E8.25_late, E8.75 and the Mab21l2 KO hold identical cells in adjacent sections (KO
     5.6%). Drop them (`vec.prep.load_stage(drop_duplicates=True)`) before neighbourhood, pseudobulk or composition
     work on prediction inputs; the emulator loads targets as released, since the server's truth includes them.
   - Expression moves with its cell.
7. **Pre-register.** Size models, composition scenarios and both submission candidates per board are fixed and hashed in
   git before the first test submission.

## 4. Board plans

Effort split: T1 ~35%, T3 ~25%, T2 heart ~15%, T2 embryo ~10-15%, shared core and harness ~10-15%.

### 4.1 T1: E12.5 single-cell population (weight 1)

Facts that shape it:
- E12.5 is "only the heart chamber" (organisers' dissection note), while E8.5-E10.5 include surrounding tissue.
  Contamination is minimal at E10.5 and less at E12.5, and the scorer counts every cell.
- Legal brackets after E13.5: Feng 2022 chamber-dissected hearts (GSE193346, E9.5 and E14.5-P9, 10x 3' v3), Feulner
  2026 whole hearts (E14.5-P7, CELLxGENE), JAX time-lapse embryos 1-24 and 50-74 (sci-RNA-seq3 nuclei), and ENCODE heart
  bulk at E14.5+. Our T1 data looks like whole cells (10x), so use within-dataset fold changes, never absolute
  external levels, and map genes by Ensembl ID.
- Most T1 genes (8,522 Gm*, 2,524 *Rik) have no external signal: hold them at E10.5 levels.

| Step | What | Effort | Kept if |
|---|---|---|---|
| B0 floor-plus | E10.5 real cells, per-type count-space momentum shift (damped E9.5 to E10.5 change), E10.5 composition with non-cardiac types reduced | 2 days | beats copy_last on the E10.5 backtest |
| M1 main | Composition model for a heart-chamber-only E12.5 (damped log-odds trends plus dissection rule) + per-type count-space shifts whose direction comes from the external bracket (E10.5 to E14.5+, t ≈ 0.5, warp calibrated on the backtest) + births of emergent states (epicardium-derived fibroblasts, coronary endothelium, smooth muscle, valve mesenchyme, definitive erythrocytes) built from the nearest E10.5 progenitors and edited toward external profiles | 1-2 weeks | beats B0 on the backtest |
| M2 diverse bet | moscot OT bracketing from our E10.5 to external E14.5 in a shared scVI latent, decoded by count sampling | 2-3 weeks | beats M1 on at least 2 of 3 pseudo-tasks |
| Stretch | Whole-ladder non-autonomous flow with growth (VGFM-style, reimplemented on torchcfm/POT) | only if time | beats M1 |

Validation: predict E10.5 from E8.5 and E9.5 (plus legal external data), scored by the emulator over 10 seeds. Bracket
fraction t ≈ 0.2-0.235 there versus ≈ 0.5 at test, so prefer parameters with a physical meaning (time fractions,
shrinkage) over free knobs.

Submissions: s1 = M1 central. s2 = the pre-registered alternative: a different non-cardiac fraction, or M2 if it wins
the backtests.

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
3. Emergent populations placed in anatomical niches on the E10.5 geometry: epicardium on the outer myocardial shell,
   epicardium-derived cells just beneath, coronary endothelium on the dorsal surface and AV groove, valve mesenchyme in
   the cushions, smooth muscle in the outflow-tract walls. Their expression comes from E10.5 precursors edited in count
   space, or from external post-E13.5 cells mapped to MERFISH (only if cross-task data use is confirmed).
4. Expression direction per type from the external bracket where available, otherwise damped momentum.
5. Grow coordinates with the pre-registered size; keep the shape (no deformation unless pseudo-targets show a gain).

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

Components (all legal under the 2026-10-15 ruling):
- **N, systematic transfer:** carry Gata4's non-Mesp1 compartments and depth/coverage shift with real cells.
- **C, composition:** deplete second-heart-field, outflow-tract and right-ventricle-lineage states, sparing first heart
  field and left ventricle (broad phenotype background). Magnitudes come from within-heart composition effects in the
  Gata4 and Mab21l2 knockouts (harmonised, bracketed, z-trimmed), counted only where they exceed the Gata4 replicate
  difference.
- **B, Wnt programme:** learn how panel genes co-vary with Wnt target activity ([redacted 2026-10-19: B], off
  the panel) inside Mesp1 cell types of the wild-type T1 scRNA, and switch it off in Mesp1 cells. This needs the
  organisers to confirm cross-task use; otherwise use MERFISH-only co-variation.
- **Weights:** set by leave-one-knockout-out (Mab21l2 to Gata4, Gata4 to Mab21l2, Gata4 rep1 to rep2) with the official
  scorer. Calibrate the response magnitude so the slope on likely DE genes is near 1.

Submissions: s1 = central combination. s2 = the pre-registered alternative composition scenario. Any rule that picks
s2 from s1's score is written down beforehand, and the organisers are asked to confirm it counts as selection.

**Provenance rule for T3 gene directions (pre-registered 2026-10-18).**
1. No per-gene direction or magnitude is written by hand. None comes from knockout literature, from curated target
   lists built from beta-catenin or Gata4 knockout or perturbation experiments, or from LLM elicitation.
2. Every per-gene change in a T3 prediction is the output of committed code applied to released data only:
   - the Gata4 and Mab21l2 knockouts against their matched wild types (systematic component and composition effect
     sizes);
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
- `vec check`, the free format checks.
- `vec baseline`, the reference baselines.
- `vec data`, the safe downloader.
- DGX environment (torch 2.14 cu130, scvi-tools, moscot), assumption ledger, data policy, external-data register.

Next:
1. Calibrate the emulator against the published reference rows on the validation boards (A2).
2. Pseudo-task harness with leave-one-transition-out splits and paired comparisons over 10 seeds.
3. Data utilities: count recovery, CSR conversion with explicit zeros removed, label harmonisation across stages.
4. External-data ingestion with window filters as committed code (JAX by embryo ID, Feng by stage).
5. Run registry; submission client that defaults to format-only and sends an idempotency key on scored uploads.

## 7. Timeline and decision points

| Dates | Work | Decision point |
|---|---|---|
| Oct 18-20 | Token and download; organiser questions out; leak report sent | D0: data on the DGX, else escalate |
| Oct 20-26 | P3 intake (test boards in index.json, Data Protocol, scorer version); emulator calibration; T3 day-1 diagnostics; density experiment; size procedure; B0 on all four boards; external data ingested | D1 (Oct 22): test-board reference stages fixed |
| Oct 27-Nov 2 | T2 embryo FGW pipeline; T1 M1 v1; T3 N+C v1; T2 heart growth v1 | D2 (Oct 31): every board beats its floor on internal validation |
| Nov 3-9 | T1 external directions and births; T3 Wnt programme and leave-one-knockout-out; T2 heart niches and scope hypotheses; T1 M2 | - |
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

`ASSUMPTIONS.md` holds 18 scored assumptions with go/no-go tests. The biggest:
- **A11**, the T3 response sign (risk 5; early proxy: a blind Gata4 test);
- **A6**, T1 E12.5 composition;
- **A8**, heart scope and size;
- **A5**, the value of the external bracket;
- **A14**, a metric change from veckit PR #6 (mitigation: aim for slope ≈ 1 on T3 and realistic cell spacing on T2);
- **A15**, the cross-task data ruling.

## 10. Open questions for the organisers

`docs/ORGANIZER_QUESTIONS.md`: cross-task use of released data in the final phase, and the test-board reference
stages, cell bounds, scorer version and deadline time zone. The leak report is a separate private message.

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
