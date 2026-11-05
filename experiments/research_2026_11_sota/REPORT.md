# The 2025-2026 state of the art against our top-3 chance: a screened sweep

> Date 2026-11-06 | AS_OF 2026-11-06 | Standard mode.
> Sources: 163 registry entries (159 distinct works) from 41 domains (`notes/registry.md`):
> - by type: 114 academic, 28 official, 21 other;
> - official and academic together: 87%;
> - largest domain: biorxiv.org at 23%, under the 25% gate;
> - the "official" label alone is 17%, below the skill's 30% gate, because a methods question is answered by papers;
> - 125 entries (77%) are dated 2025-2026.
> Seven research notes (tasks a, a2 and b-f), each screened and audited. Five board designers, each attacked by an
> adversarial critic. One target-free aggregation (`experiments/gate/path_to_target.py`). Four counter-reviewers and a
> cold devil's advocate (section 9).

Written by the AI assistant for the captain's request of 2026-11-06: research the state of the art, then synthesise
solutions that maximise our chance of a top-3 place, with his target of a total above 209.

**Compliance.**
- No agent opened a test file, a board file or the Data page, and none saw a test-board limit.
- Task e's agent met organiser statements that may describe hidden targets; they were removed (class A, below).
- No try was run.

**What the screeners and auditors removed:**
- 6 hidden-target items (class A), all from the organisers' public task pages;
- 24 in-window passages (class W);
- 1 per-board leaderboard figure (class L);
- 5 residual items found in audit;
- 2 sources withdrawn whole;
- 1 residue of a class-A removal in the committed notes, removed by the main session after the counter-review.

**Met but not recorded:**
- class B twice, in search-result summaries read by task c's agent while it searched the Mesp1-Cre phenotype;
- classes L and Y in search listings read by task e.

**Slips.** Three agents (task a's and task e's screeners, task b's auditor) opened memory files under
`~/.claude-account`, against AGENTS.md; nothing was copied. All process notes are in `docs/DATA_POLICY.md` 2026-11-06
(a).

## Executive summary

**In the captain's terms first.** Every figure here is on reading A of our estimate, the upper reading (mean 191.5).
- **The chance of 209:** P(total >= 209) is 0.42% today and 0.59% with this research's package. With everything the
  record supports, it is about 1% at most.
- **Even odds:** a sure gain of about +16 Overall.
- **On the committed outlook** (184, range 170-198) the gap is about 25.
- **Under the measurement-based T3 payoffs** of 2026-11-05 (s1's best case about 66), 209 is out of reach.
- **What 209 is:** the leading Human total on the test boards when the captain set his gate. It is reported for
  information and chooses nothing here.
- **Our own 95th percentile (205.5):**
  - the chance today is 5% by construction;
  - every draw that reaches it takes T3 from stand-in A, s1's Gata4 factor scored on the other Gata4 replicate, a fold
    the adoption rule never counts. Without that stand-in, no draw reaches even our 90th percentile (section 2).
- **The top-3 chance** stays "a few percent at most", as the 2026-11-05 research found.

**Bottom line: the 2025-2026 state of the art holds no method, dataset or tactic that moves our total materially under
our rules.**
- **The surviving items** are worth +0.01 to +0.3 Overall each in expectation.
- **The largest effect on our upper tail** is a ruling already on record and still open: T1's strong-null down-tail
  bet (the 2026-11-05 report's sheet C). At P(strong-null) = 0.2-0.3 it adds +0.5 to +1.1 points of chance at our
  95th percentile.
- **Everything together** reaches about 7-8% at our 95th percentile, against 5% today.
- **Consistent with every earlier search:** four of them found +1 to +3 at most, and this one, with 159 works mostly
  from 2025-2026, does not change that.

Why, one line per board group:
- **T1 forecasting.**
  - Only one 2025-2026 paper compares a forecaster with copy-last at 2,000 or more genes; it reports 27-36% lower W2
    [2].
  - Task a2's derivation, not the authors', is that this margin matches what denoising alone buys in high-dimensional
    W2.
  - Our scorer punishes mean-like cells: kNN smoothing multiplied the variogram distance 7.4x, and two-cell blends lost
    on both embryo backtests (`docs/research/digest_2026-10-19.md` 2.1; TRIES E8).
  - Benchmarks report no naive baseline [1]. Every growth method sets birth and death only between observed snapshots
    and is validated only on interpolation [7][8][9][10].
- **T3 perturbation.** It supports the bet s1 already makes and offers no new lever.
  - Within one screen (one run, one control), the component all perturbations share carries 23-65% of response
    energy. It is predicted for a new perturbation at r = 0.62-0.75, and the perturbation-specific rest only at 0.20-0.26
    [55].
  - A correction estimated against the scoring control is inflated 1.5-2.4x when scored against that control, and it
    loses against an independent one [53].
  - In vivo, separate-control comparisons are dominated by sample-level shifts [65][66], and MERFISH run effects are
    large and gene-structured [71][72].
  - Our two released knockout runs' factors correlate only 0.15 (TRIES N8), so s1's upper tail needs the
    beta-catenin file to share Gata4's run (U2).
- **T2 spatial.** No 2025-2026 spatial generator is compared with a copy-nearest-stage baseline on a held-out stage;
  the one lookup comparison is on held-out slices [86]. Legal 3D heart morphology exists, but our released heart files
  are mostly surrounding tissue (`experiments/coords/REPORT.md`).
- **The organisers' reference rows** sit below the floor on T1, heart extrapolation and the Gata4 board; pseudobulk_shift
  reaches 55.22 on embryo interpolation (`docs/research/digest_2026-10-19.md`). Navigo needs spliced and unspliced
  counts that our files lack [105].
- **Competition post-mortems.** Winners were hybrids tuned to the metrics, close to simple pseudobulk statistics
  [136][141][144]. Picks that chased the public score fell, and robust validation rose [158][141].

What to do (section 8):
1. **Operational, with this commit.** `experiments/gate/path_to_target.py`, its result and this folder are committed.
   Left untracked, they would have made the upload client's clean-tree guard refuse every official upload.
2. **Operational, the captain's message.** A Slack DM to Weize Xu, the private route `docs/ORGANIZER_QUESTIONS.md`
   names. It repeats the fifth email's floored-severity question (A24) word for word and asks for a public answer,
   with an optional one-line question on the embryo stage-label convention. About +0.01 to +0.05.
3. **Operational, a dated upload runbook that keeps every ruling window open** (section 8). Do not front-load the second
   round: that would forfeit the frozen fallbacks for a late A20 or Gata4 ruling.
4. **One bundled T1 ruling for the captain, by about 20 Nov** (rebuild and re-hash by 24 Nov). Each item changes T1's
   s1, so T1's s2 goes back to him:
   - **sheet C, the strong-null down-tail bet:** +0.1 to +0.3 in expectation and +0.5 to +1.1 points of chance at
     our 95th percentile at P = 0.2-0.3; tail break-even near P = 0.12;
   - an A15 reading for T1-only external gene sets, for a Feng-direction hygiene step: about +0.10;
   - a cycling-state quota within T1 categories: about +0.05.
5. **Offered, default no:**
   - *The embryo s2 at mixture fraction 1/4.* About -0.3 in expectation at the residual our committed E7.5 row implies;
     worth a ruling only if the captain holds, from released-stage evidence alone, that the embryo labels are
     plug-dated with a staging spread of 0.15 day or more.
   - *Holding the heart's momentum at test.* About +0.03; sheet E's default stays decline.
6. **Not recommended.** Re-reading the organisers' test-knockout acquisition text to decide T3's depth step: about
   +0.05, close to the 2026-11-04 exposures, and a misread costs -6 to -8 T3.

## Contents
1. Question and method
2. Where the chance lives, and what it rests on
3. T1: what the forecasting literature shows
4. T3: what perturbation research says about s1's bet
5. T2: spatial generation, organ shape and the gastrulation brackets
6. The organisers' and the community's methods
7. What separated winners in single-cell competitions
8. Solutions: candidates, verdicts and recommendations
9. Counter-review
10. Key findings
11. Limitations and future directions
12. References

---

## 1. Question and method

**The question.** Which 2025-2026 papers, reports and resources could raise our expected total, or its upside, on the
four test boards, within the hard constraints of `notes/brief.md` section 4? Those constraints are:
- the data windows and the genotype rule;
- the 2026-10-15 knockout ruling;
- the post-exposure block and A15;
- the frozen rules and the adoption rule;
- "scores choose, never compute".

**The pipeline** (the deep-research skill's P0-P7 with the team's screened route):
- *Research.* Six specialists, plus a second independent pass on task a: a2, created when the workflow restarted task a
  (section 11). The specialists covered:
  - (a, a2) single-cell temporal forecasting;
  - (b) legal external heart data and ways to create counted evidence;
  - (c) perturbation-response prediction and in vivo knockout structure;
  - (d) 3D spatial generation and organ shape;
  - (e) the organisers' work and the challenge's public record;
  - (f) post-mortems of single-cell competitions.
- *Screen and audit.* Each raw note went to a screener, who removed the barred classes (A, B, W, L, Y, Q) under the
  standing list in `docs/DATA_POLICY.md` (Scores). An auditor then:
  - removed residues;
  - spot-checked 4 claims per note against their sources: 28 in all, 17 verified, 9 corrected, 2 not verifiable;
  - corrected 32 legality statements;
  - made 19 duplicate flags (on 16 levers and one mapping entry) against the pruned list.
  The main agent read the audited notes only (`notes/task-*.md`).
- *Design and critique.* Five designers (T1, T3, heart, embryo, cross-board) read the audited notes and the
  repository, with no web, and proposed 28 candidates. A critic per designer tried to kill each one
  (`notes/design/`).
- *Counter-review.* Four reviewers (claims, compliance, value, completeness) and a devil's advocate attacked the draft
  (`notes/review/`; section 9).

**Confidence:** High for the process; medium for coverage. ICLR 2027 submissions could not be read (an OpenReview browser
check), and the NeurIPS 2026 accepted-paper list was not scanned (section 11).

## 2. Where the chance lives, and what it rests on

`experiments/gate/path_to_target.py` draws reading A of our committed estimate 200,000 times with the estimator's own
loaders and seed. It is aggregation, not a try; 20 self-checks pass, including the 191.51 reproduction
(`experiments/gate/results/path_to_target_2026-11-06.json`).

**Our own distribution.**
- The mean is 191.5, with each board at the better of its two uploads.
- Our 90th, 95th and 99th percentiles are 203.95, 205.54 and 208.07.
- Reading A sits 7.5 above the committed outlook of 184; 6.6 of the difference is on T1, where reading A credits the
  oracle-base grid in full (`experiments/gate/README.md`).

**Which boards carry the upper tail.** In the draws that reach our 95th percentile:

| Board | Mean overall | Mean when total >= our 95th percentile |
|---|---|---|
| T1 | 69.6 | 72.4 |
| Heart | 53.9 | 58.5 |
| Embryo | 60.9 | 62.3 |
| T3 | 64.6 | 74.3 |

**What the tail rests on.** Every draw at or above our 95th percentile has T3 from stand-in A. That is s1's Gata4
factor and depth step, scored on the other Gata4 replicate, a fold the adoption rule never counts. The same draws with
stand-in A's weight lowered (the record's weight is 0.5; B and K2a share the rest evenly):

| Stand-in A's weight | T3 mean | Overall mean | P(>= our 90th) | P(>= our 95th) |
|---|---|---|---|---|
| 0.5 (the record) | 64.6 | 191.5 | 10.0% | 5.0% |
| 0.375 | 62.4 | 189.3 | 7.5% | 3.7% |
| 0.25 | 60.2 | 187.1 | 5.0% | 2.5% |
| 0 | 55.8 | 182.8 | 0.0% | 0.0% |

The value reviewer adds two working answers the figures assume (`notes/review/value.md`, issue 3; not recomputed
here):
- *Per-board final marks* (working answer 0.55). One mark per bundle would cost about 4.1 Overall, T3 included, and
  cut P(>= our 95th) to about 2.1%.
- *Skill 0.5 for a floored severity* (0.65). A zero floor costs about -2.8 T3 in expectation.
Each condition is larger than this research's whole package.

**How much a sure gain is worth.** A sure Overall gain g, any mix of boards, moves the chance of reaching our own
percentiles of today as follows. A half-weight board counts half, so +10 heart points is +5 Overall. An expected gain
is applied here as a sure shift.

| Expected gain applied as a sure shift | P(>= our 90th, 203.95) | P(>= our 95th, 205.54) |
|---|---|---|
| none (today) | 10% | 5% |
| +0.35 (the cross-board items) | 11.3% | 5.9% |
| +0.8 (with R1 and the T1 sheet-C bet) | 13.2% | 7.3% |
| +1 | 14.0% | 7.9% |
| +3 | 23.2% | 15.8% |
| +5 | 31.8% | 25.1% |
| +10 | 47.3% | 43.0% |

The curve is convex here. +1 buys 2.9 points of chance at our 95th percentile, and -1 costs 2.1 (value review). So a
right-skewed bet can raise the chance even with a mean near zero; section 8 prices the rulings by both measures.

**For information only, the captain's 209.** It was the leading Human total on the test boards when he set his gate
(`docs/PLAN.md` section 8). It chooses nothing here.

| Sure Overall gain | P(total >= 209) |
|---|---|
| none (today) | 0.42% |
| +0.35 | 0.59% |
| +0.8 | 0.89% |
| +5 | 9.8% |
| +10 | 31.7% |
| +16 | 49.7% |

**Confidence:** High for the arithmetic, which reuses the estimator's loaders and self-checks. Medium for its meaning:
the estimate inherits its stand-ins and working answers (above).

**Opposing reading.** Under task h's measurement-based T3 payoffs (s1's best case about 66;
`experiments/research_2026_11/REPORT.md` section 1), the chance of reaching even our 95th percentile is under 1%.

## 3. T1: what the forecasting literature shows

**The one copy-last comparison is mostly denoising.**
- scFM (arXiv, May 2026; venue unconfirmed) is the only 2025-2026 population forecaster that reports copy-last at
  2,000 or more genes, outside in-window evaluations. Its Sinkhorn W2 on extrapolated points is 27-36% below
  copy-last, 1-4 steps past the last training snapshot, and scNODE is within 1-10% of it [2][21].
- Task a2's derivation, not the authors': a mean-like decoded output costs about sigma*sqrt(d) in W2 against real
  cells, and a real-cell copy about sqrt(2)*sigma*sqrt(d). That predicts a ratio near 0.71; the observed ratios are
  0.64-0.73, with a mean of 0.69 [2].
- Our scorer punishes mean-like cells: kNN smoothing multiplied the variogram distance 7.4x
  (`docs/research/digest_2026-10-19.md` 2.1), and two-cell blends lost on both embryo backtests (TRIES E8). The
  organisers' pseudobulk_shift, an additive log shift that fills zeros, loses about 2.2 points to the variogram.

**Benchmarks and newer methods show no naive-baseline wins either.**
- scTimeBench (ISMB 2026) runs ten methods on eight datasets, one step ahead, with no copy-last or nearest-stage
  baseline; rankings flip with the metric [1].
- CellTempo beats a delta extrapolation only within overlapping spreads (pseudotime PCC 0.74 against 0.63) [3].
- PACE (NeurIPS 2026) reports 23.7% over learned competitors on low-dimensional held-out snapshots [15]. LGP-OT
  reports 1-4% over the best learned baseline in task a's reading and at most 2% in task a2's [16]. Neither runs a
  naive baseline.

**Growth is identified only between observed snapshots.**
- WFR-FM (ICLR 2026) regresses growth from a semi-coupling of adjacent snapshots; its held-out interpolation gains are
  1-4% W1, with one loss [7].
- MUST-FM adds lineage masks [8]. ARTEMIS fits kill rates to cell-number changes [9], and scMTG weights cells by
  proliferation gene sets, with no ablation [10]. USB and SUDO report no held-out test in their abstracts [24][27].
- None sets a rate past the last observed stage, and none creates a population with no source.

**Count realism is unvalidated at held-out times.** count-FM and scLDM report unconditional fits only [18][19].

**Pretrained weights.** Chreode and CellTempo plausibly contain mouse cells inside the window [13][3]. TabPFN-2.6 is
pretrained on synthetic tables only, under a licence that bars commercial or production use of its outputs [23].

**What transfers, and when.**
- *Within one study only.* A temporal signature from one study has 1.52 weeks error within it and 8.03 weeks outside
  it; an elastic net pooled over 11 studies holds 2.65 weeks [6].
- *Only for consistent genes.* Trends carry to unseen cell types only for genes whose direction agrees across types:
  MASE 0.78 against carry-forward above a consistency of 0.8, and no model beats carry-forward below 0.2 [5].
- *Within one protocol.* Composition is forecastable inside one protocol: TabPFN reaches Spearman 0.808 and R2 0.704 at
  a held-out zebrafish stage, against a mean baseline's 0.352. No copy-last baseline was run [4].

**New legal resources, each pending its register row and staging checks:**
- *Leonard et al. (2025).* E14.5 C57BL/6 hearts by 10x Multiome nuclei, T1's chemistry: 12,507 nuclei pooled over four
  genotype groups (WT and three mutant libraries), about 29% cardiomyocytes; GSE272154 [50].
  - Motivation only: never an input to a share prior.
  - Legal by its E-day label, pending the Q7 Theiler-stage check of its sample metadata and an `EXTERNAL_DATA.md` row.
- *mdCLA (2026).* A 10x v3 whole-cell atlas of 351 C57BL/6 embryos at 37 time points: 10 late points from E14.0 to
  E18.5 half a day apart, birth, and 18 points at or before E9.5 [33][34].
  - E14.0 is subject to our staging rule, which keeps JAX only from E14.25.
  - Its E7.5, E7.75 and E8.5 points sit in the embryo and heart-interpolation windows, so it is T1-only, or those
    points are removed by our code for any shared-core step.
- *Feulner et al. (2026, Version of Record).* Whole hearts by Drop-seq from E14.5 to P7; its E14.5 sample stays out under
  the staging rule [31].

**Why none of them gives counted evidence.** The T1 critic's decisive reason: one lab's single stage (Leonard) leaves
no within-dataset time contrast to cancel its lab and dissection offset, and its pass on our E9.5 -> E10.5 backtest
would read that backtest's own confound, embryo portion to dissected heart (`notes/design/t1_critique.md`). No study
measures how far such offsets move a per-type temporal direction (task b, Gaps). Older off-organ comparisons show
protocol moving cell-type shares by tens of points [45][46][47], [46] read at snippet level. They bear on share priors,
not directly on directions.

**Blood.** Legal-stage erythroid anchors:
- no primitive enucleation at or before E9.5, and more than 95% enucleated by E15.5 [40] (not verifiable in audit);
- definitive cells above 94% of circulating cells at E16.5 [39] (not checked in audit);
- the nucleated fraction falls from about 0.34 to 0.01 across the E14.5 litter harvest [38] (verified).
No counted backtest can show a blood term: blood is 2.6% of chamber cells in the E9.5 inputs both counted tasks draw
(TRIES row R1), against about 24% at E10.5 (`docs/PLAN.md` 4.1).

**Confidence:** High that no 2025-2026 T1 method has shown a naive-baseline win that would carry to our metrics.
Medium on coverage of very recent preprints.

**Opposing reading.** If scFM samples with noise, its margin is real skill; a split-half ceiling would settle it, and
the paper reports none [2]. Even then, its gains come from dense ladders of 6-13 snapshots (hours in zebrafish and
Drosophila, 24-48 h in reprogramming), against our three stages. A second reading (task a2, Gaps): mmd_u's
truth-fitted 30-PC space projects much per-cell noise out, so denoising may cost little on T1's largest-weight metric.
Only the variogram penalty is certain.

## 4. T3: what perturbation research says about s1's bet

**Within one run, shared components dominate and transfer; across contexts only a conserved part does.**
- *Four CRISPRi screens.* Split-half-corrected response variance divides into a shared cell-line template (27.8%), a
  conserved perturbation part (29.4%), an interaction (23.5%) and noise (19.3%). The template is 23-65% of response
  energy and lies 92-98% in the top control PCs [55].
- *Within a line.* The template is predicted for a new perturbation at r = 0.62-0.75, and template-removed responses
  only at 0.20-0.26 [55].
- *Across contexts.* Zero-shot transfer recovers only the conserved part: r = 0.25-0.39 with an MLP, 0.10-0.13 with
  ridge [55]. State, with 30% of perturbations seen in the new context, still loses to the context mean by 12% [60].
  PerturbMap's transport cuts MSE by 4.1% [61].
- *A common additive term.* In a 15-line by 88-drug atlas, one term common to all drugs of a context carries a median
  48% of response variance; removing it from prediction and truth collapses a calibration-mean score from 0.590 to
  -0.021 [53].
- *Shared-control inflation.* Scoring against a shared control inflates methods that read that control 1.5-2.4x [53].
- *Few-shot adaptation.* Gains stay tiny even with many measured responses: oracle-weighted adaptation rises only from
  0.0128 (k = 2) to 0.0179 (k = 40), and the best fixed rule, 0.3 shrinkage, gives 0.0148 at k = 5. Against
  independent controls, the calibration mean lost (0.340-0.419 against 0.447-0.463) [53].

**In vivo, the separate control is the shared component.**
- *Drosophila midgut (scPT-seq).* Edited and unedited cells of perturbed animals changed alike against control animals,
  while within-animal genotype contrasts were nearly orthogonal [65].
- *Mouse liver (Perturb-Multimodal).* With same-mouse control guides, 109 of 406 targeting guides against 0 of 50
  control guides changed transcriptional state [66].
- *Zebrafish TF mutants (MIC-Drop-seq).* The median mutant changed the abundance of one major cell type but the
  transcription of 12 [63].
- *ZSCAPE.* DE sets are power-limited [64].

**MERFISH run effects are large and gene-structured.**
- Across 641 MERSCOPE brain sections, detection differs about 2-fold between sections of one mouse and falls with depth
  in the section. In 31% of sections, dropped images removed groups of 60-104 genes sharing one readout probe [72].
- Per-gene detection depends on where a barcode's bits fall in the imaging rounds [71], and per-gene ratios between
  chemistries are not uniform [73].
- Segmentation admixture and gene-count normalisation create false within-type DE [70][74].

**What this means for s1.** s1 is WT cells moved by the Gata4 factor, measured where the knockout cannot act, plus
Gata4's depth step and a lineage depletion. It bets that the beta-catenin file shares Gata4's run-level component. The
literature makes such components common and large, but run-specific:
- our two released knockout runs' factors correlate only 0.15 (TRIES N8);
- with one or two observed runs, nothing estimates the shared part better than what s1 already uses.
So the research confirms the bet's logic and its fragility (A23, U2) and adds no lever. The T3 designer's four
candidates were cut to about +0.05 to +0.15 together (section 8).

**Phenotype only** (allowed by the 2026-10-15 ruling). After Mesp1-Cre-driven loss of canonical Wnt signalling:
- the crescent and linear heart tube form, but looping and right-ventricle formation fail;
- the result is a primitive left ventricle with a shortened outflow tract;
- atria and inflow form, and no general growth retardation is reported [78].
How this bears on s1's depletion:
- it agrees with sparing first-heart-field and left-ventricle states and with depleting OFT/RV cardiomyocytes;
- it gives the weakest support to depleting the progenitor states aSHF, pSHF, aPHM and pPHM, because pharyngeal
  mesoderm is reported grossly normal (task c, mapping 11).
It rests on one 2007 paper, read with a no-genes prompt; related second-heart-field Cre losses were left unopened to
limit class-B exposure. T3's composition is closed, so this is risk information only.

**Metrics.**
- DEG overlap has limited dynamic range, and DEG counts track statistical power [57].
- Control-referenced deltas reward collapse to the mean [58][59].
- A shared high-minus-low perturbation-strength pattern reproduces across perturbations [56].
- Of T3's three DE metrics only the severity slope is scale-sensitive; the mutant-state term (mmd_u, variogram) is a
  cell-level distance. s1's depth step is the part that serves the slope.

**Confidence:** High on the structure of shared components. Low on how it maps to the beta-catenin file, which no
allowed source describes.

**Opposing reading.** Two embryo knockout runs imaged apart from the WT (N6) may share less than a screen's template.
That is the U2 risk already priced in the 2026-11-06 ruling.

## 5. T2: spatial generation, organ shape and the gastrulation brackets

**No held-out-stage generator beats a copy:**
- *ContextFlow (2026):* about 10% weighted W2 over the best flow-matching baseline on a held-out middle stage; it
  generates expression only and runs no copy, identity or linear baseline [79].
- *NicheFlow (NeurIPS 2025):* tested in distribution only [80].
- *stVCR (Nature Methods 2026):* checked by leave-one-time-point-out interpolation against TIGON, SF2M and static OT,
  with no copy-nearest baseline [81].
- *TP-DATE (September 2026):* beats FGW-CFM, a learned static FGW coupling with displacement paths. No paper runs the
  closed-form FGW-plus-displacement version on real cells, which is our s1's class [82].
- *Mimyr (ICLR 2026):* the one generator with a lookup baseline; soft Spearman rises from 0.39-0.48 to 0.45-0.58 on
  held-out brain slices, not held-out stages [86].

**Isotropic 3D reconstruction would hurt us.**
- SINTER3D, isoST and UniST resample serial sections into isotropic volumes [84][85][88].
- Our files keep z as a section index times a per-file step, and a 1.5x z-factor mismatch costs 12.6-12.8 heart points
  on an otherwise perfect prediction (`docs/research/digest_2026-10-19.md` 4.4).

**Organ shape cannot be borrowed.**
- Legal late 3D heart morphology exists:
  - DMDD HREM, 215 wild-type C57BL/6N E14.5 embryos at 3 um voxels [94];
  - an E15.5 micro-CT average [96];
  - eMouseAtlas MRI models [97].
- But our released heart files are not the organ (`experiments/coords/REPORT.md`; `docs/research/digest_2026-10-19.md`
  4.4):
  - at E9.5, 63% of cells are non-cardiac shell;
  - every released file has 30-40 sections (the test stacks' counts are unknown and not sought);
  - a core crop alone moves Dice to 0.79.
- At most, whole-organ morphology could set a core sub-cloud's shape, and size and scope are closed.

**Gastrulation brackets.**
- Legal spatial atlases near E7.75 are thin: SEU-3D only at E8.0 [98], and Harland et al. only at or before E7.25 [99]
  (metadata only; its stage list unresolved).
- Every morphodynamic model of the interval is built on embryos inside the embryo window [91] (cited for method only).
- Flysta3D-v2 offers a window-free Drosophila serial-section series [103]. The embryo critic priced an asymmetry
  check on it at most +0.02 in decision value.

**Confidence:** High that no spatial SOTA method has a naive-baseline win at a held-out stage.

**Opposing reading.** A learned spatial model trained across our released heart ladder could in principle learn growth
and shape. But every released stage changes its imaged field non-monotonically, so it would learn the field, not the
heart (`notes/design/heart_critique.md`).

## 6. The organisers' and the community's methods

**What the organisers ask for.** Their Task 1 page names "supervision across the whole stage ladder, a non-autonomous
field, and a treatment of populations appearing rather than only existing ones moving" [109].

**Their own references miss that bar.**
- neural_ode "does not beat the one-line constant shift" [109][112].
- Their published practice-board rows (`docs/research/digest_2026-10-19.md`):
  - pseudobulk_shift scores 48.58 on T1 and 47.76 on heart extrapolation, below the floor, and 55.22 on embryo
    interpolation, above it;
  - shift_transfer and gene_ko score 42.48 and 42.50 on the Gata4 board.
- Navigo, the probable dynode_flow, is a non-autonomous rectified flow trained on 12.4M cells at 43 time points [105].
  It needs spliced and unspliced counts, has no proliferation, death or unbalanced transport, and is validated by
  interpolation only [105][122].

**The T3 page.** It lists two pitfalls, "ignoring the stage difference between the training knockout and the evaluation
ones" and copying the training delta, and asks models to "separate the wild-type developmental trajectory from the
genetic effect" [111].
- s1 measures its factor in compartments the knockout cannot touch, at the test stage's WT.
- The barred s1, which a "zero" answer would make T3's s2, applies the E9.5 Mab21l2 factor at E8.75. Its Gata4 credit
  is mostly that stage offset, the pitfall the page names.
- Task e's lever 1 prices this risk. The committed readings do not show the organisers' 7.5-point loss, and N11
  already found the swap's gain sensitive to it (+1.38 falling to +0.42).

**The scorer and the rules have not changed.**
- veckit main is still cc3859a of 2026-10-02, the commit we pin [114].
- An organiser's open PR says the server checks only that X is finite and non-negative [114].
- No public change after 2026-10-28 was found on the Rules, FAQ or Timeline [115][116][117].

**The proposal** scores T3's distribution term on the mutant state directly, and plans paired bootstrap comparisons of
the leading entries [108].

**Community material** is mostly tooling. Its method ideas rest on target-aware toys or on P2 validation-score
selection [128][129][130]. Several repositories carry no licence, and some size steps rest on board anchors; those were
skipped [130][133].

**Withheld by our screen.** The Task 1 and Task 2 pages carry organiser statements about the test targets; the screen
removed them (class A). Whether public organiser guidance may inform a legal-route prior is the captain's question (section 11).

**Confidence:** High. Nothing in the organisers' public record points to a method class we have not tried or pruned.

**Opposing reading.** The organisers' framing (whole ladder, births) may describe what wins at test even though no
public model shows it. Our evidence cannot test that, and a leaderboard read may not be used to.

## 7. What separated winners in single-cell competitions

**Virtual Cell Challenge 2025** (H1 cells, 300 CRISPRi perturbations, more than 1,200 teams):
- *The winners were hybrids tuned to the metrics.* The winner fused PDS, DES and a small MAE loss and repeated outputs
  to meet DES; second place used residual delta learning aimed at PDS; third place transferred summary statistics
  with a CV-tuned global scale [136][160].
- *The scoring shaped the field.* The aggregate clipped below-baseline metrics to zero, so MAE dropped out; without
  the clipping, every top-50 entry would have scored below zero [136][140]. DES plateaued near 40%, against about 39%
  for a technical replicate [140][142].
- *Simple methods sat close to the top.*
  - A 1970 regression finished 15th of 337 teams, a rank read from a search summary [141].
  - A pseudobulk kernel ridge finished 38th [140].
  - On the released test split, a deterministic centred transfer scored 0.3413 against the winning entry's 0.3401,
    after four frozen versions ran in sequence on that split, so not a single blind test [144].
- *What carries over to us.* T3's DE metrics are rank-based except the slope, so output-repetition and scale tricks
  have no analogue here. The tactic that transfers is matching the metric's own reference, which is what s1's
  factor-plus-depth already does.

**Kaggle and NeurIPS single-cell competitions:**
- *2022 multimodal:* winners beat a kNN baseline by only +0.06 and +0.03 Pearson. The public-test score was the best
  predictor of private rank (Spearman 0.944), yet one entry fell from 2nd to 133rd [148].
- *2021:* simple kernel ridge and linear methods won key subtasks, and re-running the overall winner moved it by up to
  8 ranks [149].
- *2023 perturbations:* the top four used target encodings, mixed losses, pseudo-labelling and median ensembles [151],
  read from an unchecked mirror.
- *Broad obesity challenge:* the MMD runner-up resampled real residuals around a predicted centroid, after Gaussian
  noise was "rejected by MMD" [155][156]. Our real-cell, count-space construction is the same idea.

**Final picks: score-chasing fails; that boldness pays is not shown.**
- M5 teams (a 2020 contest) that tuned multipliers to the public board fell, while teams relying on robust
  cross-validation rose [158].
- A BirdCLEF note kept a decorrelated second pick [159]; it was read as a search summary.
- Turbine led VCC's live board and finished 15th [141].

**Confidence:** Medium-high. Several sources are participant blogs or summaries [140][141][142][151][159]. The VCC 2026
results (a six-metric average over three held-out cell lines) are due mid-to-late November 2026 [146].

**Opposing reading.** Final picks at Kaggle scale say little about two uploads per board. The contest-economics case
for adding variance (`experiments/research_2026_11/REPORT.md` section 1) still holds. The bold options still open are
small and rest on the captain's priors: T1's sheet C (about +0.3 at a prior of 0.3) and the embryo s2 at 1/4
(about -0.3 at the committed residual).

## 8. Solutions: candidates, verdicts and recommendations

The designers proposed 28 candidates. The critics' verdicts and corrected gains (Overall units) are below; the
counter-review revised two of them (section 9). Files: `notes/design/`.

**T1** (`t1.md`, `t1_critique.md`):

| Candidate | Verdict | Corrected gain |
|---|---|---|
| A. Second direction from Leonard's E14.5 nuclei, pooled with Feng's | KILL | about 0 (-0.3 to +0.2) |
| B. Leonard nuclei drawn as a fraction t of each category quota | KILL | about -0.2 (-1.0 to +0.3) |
| C. Legal blood and cardiomyocyte anchors priced on the grid: keep today's pair | KEEP | 0 (status quo) |
| D. Cycling-state quota within categories | RULING-ONLY | about +0.05 (-1.0 to +0.4) |
| E. Larger fibroblast-like share (births revisited) | KILL | -0.2 if built (-2 to +1) |

**T3** (`t3.md`, `t3_critique.md`):

| Candidate | Verdict | Corrected gain |
|---|---|---|
| C1. Re-read the test-knockout acquisition text to decide the depth step | RULING-ONLY | about +0.05 |
| C2. Repeat the floored-severity question by Slack DM | KEEP-REVISED | +0.01 to +0.05 |
| C3. Section-artefact audit and a section-robust factor | KILL | about 0 |
| C4. Lineage-distance source factor | KILL | about 0 |

C1 would:
- reverse the 2026-11-01 (b) decision;
- reopen the closed T3 depth lever (TRIES rule 6, the QC20 entry);
- depart from the frozen T3 rules (D4);
- sit close to the 2026-11-04 exposures, by the team's own reasoning in `docs/DATA_POLICY.md` 2026-11-04 (f)(3).
C3 and C4 fail rule 5 and would taint the Gata4-free hedge.

**Heart** (`heart.md`, `heart_critique.md`):

| Candidate | Verdict | Corrected gain |
|---|---|---|
| H-A. Lineage-contrast momentum | KILL | 0 |
| H-B. Hold the momentum at test (the full hold only) | RULING-ONLY | +0.03 (-0.25 to +0.25) |
| H-C. Readout-structured run model | KILL | 0 |
| H-D. Run-noise shrinkage | KILL | 0 |
| H-E. Core shape toward a legal late anchor | KILL | 0 |
| H-F. Organiser question on external spatial data | KILL | 0 |

The designer's premise is that the lineage-shared part of B0's momentum carries the reference run's offset, with the
wrong sign at test. It is not new: the 2026-11-05 dynamics note and sheet E priced the same run term. The critic finds
it unidentifiable here, because the shared part also carries admixture, the depth drift, shared biology and, at test,
the dissection change. The cardiomyocyte-only variant rests on an uncommitted, possibly scRNA-derived line that A15
keeps out.

**Embryo** (`embryo.md`, `embryo_critique.md`):

| Candidate | Verdict | Corrected gain |
|---|---|---|
| C1. Lean-early s2 (fraction 1/4 in E7.5's geometry) | RULING-ONLY | about -0.3 at the committed residual (-0.7 to 0) |
| C2. Organiser question on the stage-label convention | KEEP-REVISED | +0.01 |
| C3. Emergence-first s2 | KILL | -0.15 |
| C4. Asymmetry check on real intermediates | KILL | 0 |
| C5. Retire the flat-geometry bet in s1 | KILL | -0.1 if adopted |

**Cross-board** (`cross.md`, `cross_critique.md`):

| Candidate | Verdict | Corrected gain |
|---|---|---|
| C1. Front-load the second round and the marks | KEEP-REVISED, as a dated runbook | +0.05 (about 0 as proposed) |
| C2. Embryo s2 at fraction 1/4 | RULING-ONLY | revised in section 9: about -0.3 |
| C3. Organiser push on A24, with one answer table | KEEP-REVISED, private route | +0.05 |
| C4a. A15 read as not covering T1-only external gene sets | RULING-ONLY | +0.10 |
| C4b. A15 answer for a cross-modal heart direction | KILL | 0 |
| C5. Captain-written legal T1 composition prior | KILL | 0 |
| C5-births. Nuclei-anchored births | KILL | 0 |
| C6. Late-fetal counting family revival | KILL | 0 |

### Recommendations

**Operational, no ruling.**
1. **The commit.** `experiments/gate/path_to_target.py`, its result and this folder are committed with this report.
   The clean-tree guard exempts only `results/`, `runs/` and `experiments/**/results/` (`src/vec/provenance.py`). After
   the sync, the DGX tree must pass the same rule (guard "dgx provenance"); its four old untracked harmonise files lie
   under `experiments/harmonise/results/` and are exempt.
2. **The DM.** The captain sends a Slack DM to Weize Xu, the private route `docs/ORGANIZER_QUESTIONS.md` names:
   - the fifth email's floored-severity question, word for word, asking as it did for a public answer in
     #virtual-embryo-challenge, and posting nothing in public himself;
   - optionally, one general line on whether the embryo-series stage labels are gestational (plug-dated) or
     morphological;
   - a second recipient only if the captain first adds that contact to `docs/ORGANIZER_QUESTIONS.md`.
   Any answer is recorded with its date in `docs/DATA_POLICY.md`. One that arrives after T3's s2 is up changes nothing
   (ruling 2026-11-06 (b)). Value: about +0.01 to +0.05, since no rules question has been answered since 2026-10-16.
3. **A dated upload runbook that keeps every frozen switch alive.** Each upload is run by the captain on his go-ahead
   at the time, and every T3 upload follows a screened check (`docs/S2_RULES.md`, T3).
   - 27 Nov 05:30 IST: rebuild at the freeze commit and checks. Then the first round, with T3's s2 chosen at the
     screened check of mail, Slack and the Rules page that ruling 2026-11-06 (b) requires (run no earlier than 26 Nov).
   - 29 Nov IST: T1, heart and embryo s2, each after its screened check.
   - 30 Nov IST: T3's s1, after its screened check.
   - 1 Dec IST: buffer, and marks as soon as both of a board's uploads are scored.

   Do not front-load: a late A20 bar or Gata4 restriction is absorbed only by the s2 and last-T3 slots. Also reconcile
   `docs/PLAN.md`, whose section 7 status note and section 8 deadline-exception bullet both say "by 2026-12-01 06:00
   UTC", with section 8's "by the end of 30 November IST" and `docs/S2_RULES.md`.

**Rulings, ranked by their effect on our upper tail** (ΔP is in points of chance at our 95th percentile):

| Ruling | Mean (Overall) | ΔP at our 95th | Deadline | Conditions |
|---|---|---|---|---|
| T1 sheet C: strong-null down-tail bet | +0.09 / +0.29 at P = 0.2 / 0.3 | +0.49 / +1.13 (break-even P about 0.12) | s1 freeze | T1's s2 back to the captain; 2026-11-05 report section 7 |
| T1 A15 reading for T1-only gene sets (Feng hygiene) | about +0.10 (-0.3 to +0.4) | about +0.3 (cross design) | about 20 Nov | see below |
| T1 cycling-state quota | about +0.05 (-1.0 to +0.4) | small | about 20 Nov | see below |
| Embryo s2 at fraction 1/4 | about -0.3 (-0.7 to 0) | -1.3 to +1.0 | the embryo's s2 | default no; see below |
| Heart: hold the momentum at test | about +0.03 (-0.25 to +0.25) | about 0 | the s1 freeze | sheet E's default: decline |

Sheet C is priced as a two-point bet: +1.75 T1 with probability P, else -0.33 (its stated sizes), read on the
survival curve of section 2 (value review; it matches `path_to_target.py`'s figures where they overlap). Cut by a third
(+1.17 if right), it is worth +0.05 to +0.5 points of chance.

The conditions:
- **T1 bundle (sheet C, A15 hygiene, cycling quota).** One ruling by about 20 Nov, so the rebuild and re-hash finish by
  24 Nov.
  - All three change T1's s1, so T1's s2 goes back to the captain. They enter by his ruling, not the adoption rule.
  - *A15 hygiene* needs an `EXTERNAL_DATA.md` row per source, with `docs/PLAN.md` 4.4-style provenance. Its gene
    classes are fixed in the ruling row before any build, none chosen from E10.5's DE genes, and R1 stays untouched.
  - *The cycling quota* is anchored on Leonard's WT nuclei only. It needs the Q7 staging check and a register row first,
    and its rule and constants go in a TRIES row before any run or any computation on E10.5's within-category shares,
    with the (f)(3) library-size report. Cell calls come from the released labels or a gene list under his A15 reading.
- **Embryo s2.** A tail bet on the staging spread, default today's s2.
  - The embryo critic's pricing stands; the cross critic's +0.15 used a loss below its own bound.
  - At the residual the committed E7.5 mixture row implies, it is about -0.3 in expectation: -0.7 at a staging spread
    of 0.05 day and 0 at 0.15.
  - It is worth a ruling only if the captain holds that the embryo labels are plug-dated with a spread of 0.15 day or
    more, about P(E7.75 is E7.5-like) of 0.37 or more.
  - His prior must rest only on committed released-stage evidence: the labels' staging spread, the 2026-11-03 depth
    diagnostic and the E7.5 and E7.25 backtests. Nothing seen on the Data page or in the board index on 2026-11-04 may
    enter it.
  - The ruling goes in its own commit with his words and the new sha256 (`docs/S2_RULES.md`, changes after the freeze),
    and it drops the closed depth reason.
- **Heart hold.** The full hold only. It reverses a component the adoption rule kept (+0.47 in 9 of 10 seeds, +0.90 in
  10 of 10), and it needs a new golden table and hash reference.

**Not recommended.** T3 C1, re-reading the acquisition text, for the reasons under the T3 table.

**Watch.**
- The VCC 2026 results (mid-to-late November), through a screening agent.
- Any organiser answer (A20, A24, A15, Q7).

**Confidence:** Medium. Each price comes from one designer and one critic, and the counter-review changed two of them.

## 9. Counter-review

**Who reviewed.** Four reviewers and a cold devil's advocate, with no web (`notes/review/`):
- claims: 1 high, 6 medium, 13 low issues;
- compliance: 3 high, 7 medium, 7 low;
- value: 2 high, 2 medium, 3 low;
- completeness: 2 high, 6 medium, 6 low.
Every issue was fixed in this version, or is answered below.

**What the review changed:**
- *Process and compliance.* The draft's process claims were wrong in three places:
  - it said class B was never met;
  - it pointed to a change note that did not yet exist;
  - it called the commit done.
  These are corrected, and the note `docs/DATA_POLICY.md` 2026-11-06 (a) is written. A class-A residue in
  `notes/task-e.md` and full session-store paths in four notes were removed. The design notes' passages that measure
  candidates against 209 are marked struck, and `path_to_target.py` no longer computes per-board 209 what-ifs.
- *Ranking.* Sheet C, the largest tail lever on record, was missing from the ruling list, which ranked by mean. The
  embryo s2 was overpriced: the embryo critic is right, and the residual in our committed numbers removes most of the
  asymmetry.
- *Hidden assumptions.* The chance figures rest on T3's in-sample stand-in A and on two working answers. Section 2 now
  says so, with the stand-in reweighting recomputed here.
- *Accuracy.* Twenty smaller fixes: citation forms, the heart critic's actual verdict, the axis of the T3 literature,
  the organisers' embryo row, the source of each figure, and weak-source caveats.

**Key controversies:**
1. **Mean or tail?** The devil's advocate puts P(MISTAKE) at 60%: a team about 2.6 SD below 209 should value spread,
   not mean. Its "bet harder on what is adopted" (for example, a larger T3 factor weight) is the same move as
   magnitude changes the frozen T3 rules close.
   - *What holds:* the ranking by ΔP at our 95th percentile, which now leads with sheet C.
   - *What we decline:* scaling the T3 factor. Three of T3's four metric groups (de_score, de_direction and the
     mutant-state distances) do not reward a larger dp if the direction is right, and the severity slope is already
     near its calibration on stand-in A.
   - *The critics' calibration:* the advocate's check, feeding a known adopted winner to a critic blind, was not run
     (section 11).
2. **Front-loading the second round.** The cross designer's +0.35 against its critic's about 0. The critic is right: the
   frozen rules absorb late rulings only in the s2 and last-T3 slots.
3. **The embryo s2 at 1/4.** The cross critic's +0.15 against the embryo critic's -0.05 and the value review's -0.3. At
   the committed residual it is a losing bet unless the staging spread is wide.
4. **Leonard's nuclei.** The T1 designer's +0.1 to +0.3 against its critic's 0 to -0.2: no within-dataset time contrast
   cancels the lab offset.

## 10. Key findings

1. **No 2025-2026 method has shown a naive-baseline win at a held-out stage that would carry to our metrics,** on any
   board. The one copy-last comparison at scale is mostly denoising, which our scorer punishes [2][1][79][86].
2. **New legal anchors exist** (Leonard 2025, mdCLA 2026) [50][33], each pending its register row and staging checks.
   Leonard is one lab's single stage, so no within-dataset contrast cancels its offsets, and no study measures how far
   such offsets move a per-type direction.
3. **T3's upside is s1's existing bet.** The literature makes shared, run-level components common and large
   [53][55][65][72], and our two knockout runs' factors correlate only 0.15. Our whole upper tail rests on T3's in-sample
   stand-in A.
4. **The cross-board items are worth about +0.35 Overall; the open T1 sheet-C bet is the largest tail lever** (+0.5 to
   +1.1 points of chance at our 95th percentile at P = 0.2-0.3). Everything together reaches about 7-8% at our 95th
   percentile, and about 1% at 209.
5. **The upload plan holds the largest preventable losses:** a clean tree on both machines on 27 Nov, and the ruling
   windows kept open.

## 11. Limitations and future directions

**Limitations.**
- *Coverage gaps:*
  - ICLR 2027 submissions could not be read (an OpenReview browser check);
  - the NeurIPS 2026 accepted-paper list was not scanned (only PACE was confirmed);
  - ICML 2026 proceedings were not checked;
  - one legal lead was not assessed: a 2026 whole-embryo endothelial time series at 0.5-day sampling (PubMed 41702400);
  - no 2025-2026 comparison of single-cell and single-nucleus data was found, and no mouse in vivo knockout
    single-cell study at E8-E9.5.
- *The staging spread that decides the embryo s2 is sourced nowhere.* A target-free diagnostic of per-embryo staging
  spread on released labels would inform it.
- *Read depth:* several sources were read only at abstract or summary level (marked in the notes).
- *The critics were not calibrated blind.* No adopted winner was fed back to a critic as a fresh candidate (the devil's
  advocate's check).
- *Prices come from designers and critics working on committed results, not from runs.*
- *Task a ran twice:* a workflow restart created a second agent, and the first agent, resumed by a coordination
  message, wrote a separate note. Both were screened and audited. They agree on every conclusion and differ on some
  figures (LGP-OT's margin, scFM's horizon, count-FM's comparators); section 3 uses the audit-corrected figure where one
  exists.
- *Organiser statements withheld.* The organisers' public Task 1 and Task 2 pages carry statements about the test
  targets, which our brief classes as hidden-target properties (A). The screen withheld them, and the main session saw
  none. Whether public organiser guidance may inform a legal-route prior (the E12.5 collection note is a precedent for
  using organiser text) is the captain's question. Not priced.

**Future directions, in deadline order:**
1. The DM (recommendation 2), and the dated runbook (recommendation 3) once the captain agrees its dates.
2. The bundled T1 ruling by about 20 Nov.
3. The embryo s2 only if the captain holds a wide staging spread on released evidence.
4. A screened read of the VCC 2026 results in late November.

## 12. References

The full registry, with source type, accessibility, date and authority, is `notes/registry.md`. Its local-to-global map
and barred-content tags are in the same file. The tags:
- (W): in-window results on the page, used for method or legal-stage points only;
- (B): gene-level knockout text on the page, read for phenotype only;
- (A): may describe a hidden target, read only through a screening prompt.

Cited above:

[1] Osakwe, Huang, Li. scTimeBench (Bioinformatics 42 Suppl 2, ISMB 2026). academic, 2026-08. https://www.biorxiv.org/content/10.64898/2026.03.16.712069v1.full
[2] Pu et al. scFM, conditional flow matching for single-cell dynamics (arXiv 2605.22340). academic, 2026-05. https://arxiv.org/html/2605.22340v1
[3] Luo et al. CellTempo (bioRxiv 2026.02.08.704720). academic, 2026-02. https://www.biorxiv.org/content/10.64898/2026.02.08.704720v1.full
[4] Palla, Hillsley, Kim, Royer. Tabular foundation models as perturbation predictors (bioRxiv 2026.06.28.735106). academic, 2026-07. https://www.biorxiv.org/content/10.64898/2026.06.28.735106v1.full
[5] Shi, Wu, Li, Duprey. Dynamic consistency and predictable genes (bioRxiv 2026.02.11.705300). academic, 2026-02. https://www.biorxiv.org/content/10.64898/2026.02.11.705300v1.full
[6] (W) Venkatesan, Werner, Li, Gillis. Cell type-agnostic maturation signatures (PLOS Biol 24(4) e3003757). academic, 2026-04. https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3003757
[7] Peng et al. WFR-FM (ICLR 2026; arXiv 2601.06810). academic, 2026-04. https://arxiv.org/html/2601.06810v2
[8] (W) Peng, Chen, Zhou. MUST-FM (arXiv 2605.16529). academic, 2026-10. https://arxiv.org/html/2605.16529v2
[9] Alatkar, Wang. ARTEMIS (Bioinformatics 41 Suppl 1). academic, 2025-07. https://pmc.ncbi.nlm.nih.gov/articles/PMC12261462/
[10] (W) Cui et al. scMTG (bioRxiv 2026.06.04.730241). academic, 2026-06. https://www.biorxiv.org/content/10.64898/2026.06.04.730241v1.full
[13] Qiu M et al. Chreode (arXiv 2605.28111). academic, 2026-05. https://arxiv.org/abs/2605.28111
[15] Yu et al. PACE (NeurIPS 2026; arXiv 2605.18587). academic, 2026-10. https://arxiv.org/abs/2605.18587
[16] Balik, Lähdesmäki. LGP-OT (arXiv 2605.20989). academic, 2026-05. https://arxiv.org/html/2605.20989v1
[18] Wei, Pearson. Flow matching for count data (arXiv 2605.07746). academic, 2026-10. https://arxiv.org/html/2605.07746v3
[19] Palla et al. scLDM (arXiv 2511.02986). academic, 2026-06. https://arxiv.org/html/2511.02986v2
[21] Zhang et al. scNODE (Bioinformatics 40 Suppl 2). academic, 2024-09. https://pmc.ncbi.nlm.nih.gov/articles/PMC11373355/
[23] Prior Labs. TabPFN-2.6 model card and licence. official, 2026-11. https://huggingface.co/Prior-Labs/tabpfn_2_6
[24] Ying et al. USB, discrete branching dynamics (arXiv 2605.00545). academic, 2026-06. https://arxiv.org/abs/2605.00545
[27] Ying et al. SUDO, unbalanced dynamic OT with general growth penalty (arXiv 2609.04710). academic, 2026-10. https://arxiv.org/abs/2609.04710
[31] Feulner et al. Transcriptional dynamics of the murine heart during perinatal development (Biol Open 15(9) bio062493). academic, 2026-09. https://journals.biologists.com/bio/article/15/9/bio062493/372856/Transcriptional-dynamics-of-the-murine-heart
[33] Cao S et al. Developmental chronology of mouse embryo from 2-cell stage through birth, mdCLA (Nat Cell Biol 28(7)). academic, 2026-06. https://pmc.ncbi.nlm.nih.gov/articles/PMC13364681/
[34] CNCB-NGDC OMIX014646, scRNA-seq and Xenium data of mdCLA. official, 2026-03. https://ngdc.cncb.ac.cn/omix/release/OMIX014646
[38] (W) Korporaal et al. Hemoglobin switching in mice carrying the Klf1 Nan variant (Haematologica 106(2)). academic, 2021-02. https://digital.haematologica.org/wp-content/uploads/flipbook/67/files/basic-html/page153.html
[39] (W) McGrath et al. A transient definitive erythroid lineage (Blood 117(17)). academic, 2011. https://pmc.ncbi.nlm.nih.gov/articles/PMC3099576/
[40] (W) Isern, Fraser, He, Baron. The fetal liver is a niche for maturation of primitive erythroid cells (PNAS 2008). academic, 2008. https://pmc.ncbi.nlm.nih.gov/articles/PMC2373346/
[45] Dewing et al. Cellular composition of the early fetal human heart (PLoS One 17(11)). academic, 2022-11. https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0259477
[46] Galow, Wolfien et al. Integrative cluster analysis of whole hearts (Cells 9(5)). academic, 2020-05. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7291011/
[47] Renaut et al. Paired single-cell and single-nucleus lung samples (PLoS Genet 2024). academic, 2024-06. https://pmc.ncbi.nlm.nih.gov/articles/PMC11166281/
[50] Leonard et al. Single-nuclei multiomics of a murine model of cardiac development (Nat Commun 16:6947; GSE272154). academic, 2025-07. https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12304165/fullTextXML
[53] Liang et al. A split-control criterion for few-shot adaptation (bioRxiv 2026.10.05.756906). academic, 2026-10. https://www.biorxiv.org/content/10.64898/2026.10.05.756906v1.full
[55] Molina, Zhang. Perturbation response decomposition (bioRxiv 2026.07.24.740459). academic, 2026-07. https://www.biorxiv.org/content/10.64898/2026.07.24.740459v1.full
[56] Sharma et al. Averaging hides perturbation-strength variation (bioRxiv 2026.10.20.760537). academic, 2026-10. https://www.biorxiv.org/content/10.64898/2026.10.20.760537v1.full
[57] Schäfer et al. Towards principled evaluation of perturbation prediction (bioRxiv 2026.07.23.740433). academic, 2026-07. https://www.biorxiv.org/content/10.64898/2026.07.23.740433v1.full
[58] Miller et al. Deep perturbation models on well-calibrated metrics (bioRxiv 2025.10.20.683304). academic, 2025-10. https://www.biorxiv.org/content/10.1101/2025.10.20.683304v1.full
[59] Mejia et al. Diversity by design (arXiv 2506.22641). academic, 2025-06. https://arxiv.org/abs/2506.22641
[60] Adduri et al. State (bioRxiv 2025.06.26.661135). academic, 2025-06. https://www.biorxiv.org/content/10.1101/2025.06.26.661135v2.full
[61] Cui, Liu, Sun. PerturbMap (arXiv 2607.28090). academic, 2026-07. https://arxiv.org/abs/2607.28090
[63] Carey et al. MIC-Drop-seq (Nat Commun 17:4738). academic, 2026-04. https://pmc.ncbi.nlm.nih.gov/articles/PMC13216532/
[64] Saunders et al. ZSCAPE (Nature 623). academic, 2023-11. https://pmc.ncbi.nlm.nih.gov/articles/PMC10665197/
[65] Hawkins et al. scPT-seq (bioRxiv 2025.12.23.696319). academic, 2025-12. https://www.biorxiv.org/content/10.64898/2025.12.23.696319v1.full
[66] Saunders RA et al. Perturb-Multimodal (Cell 188). academic, 2025-06. https://pmc.ncbi.nlm.nih.gov/articles/PMC12324982/
[70] Mitchel et al. Impact of segmentation errors in spatial transcriptomics (bioRxiv 2025.01.02.631135). academic, 2025-01. https://www.biorxiv.org/content/10.1101/2025.01.02.631135v1.full
[71] Luce et al. Protocol optimization for multiplexed RNA imaging (Sci Rep 15). academic, 2025-03. https://www.biorxiv.org/content/10.1101/2025.03.21.644527v1.full
[72] Martin et al. MerQuaCo (eLife reviewed preprint 105149). academic, 2025-10. https://elifesciences.org/reviewed-preprints/105149
[73] He et al. MERFISH 2.0 (bioRxiv 2026.03.06.710199, vendor-authored). other, 2026-03. https://www.biorxiv.org/content/10.64898/2026.03.06.710199v2.full
[74] Atta et al. Gene count normalization in imaging-based spatial transcriptomics (Genome Biol 25). academic, 2024-06. https://genomebiology.biomedcentral.com/articles/10.1186/s13059-024-03303-w
[78] (B) Klaus et al. Distinct roles of Wnt/beta-catenin and Bmp signaling during early cardiogenesis (PNAS 104). academic, 2007-11. https://pmc.ncbi.nlm.nih.gov/articles/PMC2141811
[79] (W) Rathod et al. ContextFlow (arXiv 2510.02952). academic, 2026-05. https://arxiv.org/html/2510.02952v3
[80] (W) Sakalyan et al. NicheFlow (NeurIPS 2025). academic, 2025-11. https://arxiv.org/html/2511.00977v2
[81] Peng, Zhou, Li. stVCR (Nature Methods 23). academic, 2026-03. https://www.biorxiv.org/content/10.1101/2024.06.02.596937.full
[82] (W) Ying et al. TP-DATE (arXiv 2609.20008). academic, 2026-09. https://arxiv.org/html/2609.20008v2
[84] Zhang et al. SINTER3D (Genome Biology 27). academic, 2026-06. https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13508360/fullTextXML
[85] (W) Li et al. isoST (bioRxiv 2025.08.15.670472). academic, 2025-08. https://www.biorxiv.org/content/10.1101/2025.08.15.670472v1.full
[86] Deshpande, Bei, Ma, Krieger. Mimyr (ICLR 2026). academic, 2025-11. https://www.biorxiv.org/content/10.1101/2025.11.24.690239v1.full
[88] Shui et al. UniST (bioRxiv 2026.03.12.711362). academic, 2026-03. https://www.biorxiv.org/content/10.64898/2026.03.12.711362v1.full
[91] Esteban et al. Pseudodynamic analysis of heart tube formation in the mouse (Nat Cardiovasc Res 1; abstract only). academic, 2022-05. https://doi.org/10.1038/s44161-022-00065-1
[94] Geyer et al. A staging system for mouse embryos harvested at E14.5 (J Anat 230). academic, 2017-02. https://pubmed.ncbi.nlm.nih.gov/28185240/
[96] Wong et al. A novel 3D mouse embryo atlas based on micro-CT (Development). academic, 2012-09. https://doi.org/10.1242/dev.082016
[97] eMouseAtlas EMA 3D model summary. official, 2026-11. https://emouseatlas.org/emap/ema/modelsummary/modelsummary.html
[98] (W) SEU-3D digital embryos (Cell 2025; not opened). academic, 2025-06. https://www.cell.com/cell/fulltext/S0092-8674(25)00622-1
[99] Harland et al. Spatiotemporal atlas of mouse gastrulation (Cell Reports 2025; metadata only). academic, 2025-07. https://www.biorxiv.org/content/10.1101/2024.12.12.628063v1
[103] Flysta3D-v2 Drosophila 3D spatiotemporal atlas (CNGB STOmics). academic, 2025-06. https://db.cngb.org/stomics/flysta3d-v2/resource
[105] (W) Fan et al. Generative modeling of mouse embryogenesis (Navigo; bioRxiv 2026.06.18.733286). academic, 2026-06. https://www.biorxiv.org/content/10.64898/2026.06.18.733286v1
[108] (A) Ding et al. The Virtual Embryo Challenge, NeurIPS 2026 competition proposal. official, 2026-08. https://virtualembryo.ai/challenge/paper/neurips-competition
[109] (A) VEC organisers. Task 1 (Temporal) page. official, 2026-11. https://virtualembryo.ai/challenge/tasks/temporal
[111] (A) VEC organisers. Task 3 (Perturbation) page. official, 2026-11. https://virtualembryo.ai/challenge/tasks/perturbation
[112] VEC organisers. Reference rows page. official, 2026-11. https://virtualembryo.ai/challenge/baselines
[114] aristoteleo/veckit (MIT): commits, issues and PRs. official, 2026-10. https://github.com/aristoteleo/veckit
[115] VEC Official Rules. official, 2026-11. https://virtualembryo.ai/challenge/rules
[116] VEC FAQ. official, 2026-11. https://virtualembryo.ai/challenge/faq
[117] VEC Timeline. official, 2026-11. https://virtualembryo.ai/challenge/timeline
[122] aristoteleo/Navigo-release (BSD-2-Clause). official, 2026-06. https://github.com/aristoteleo/Navigo-release
[128] i-habib/virtual-embryo-community (MIT). community, 2026-11. https://github.com/i-habib/virtual-embryo-community
[129] wqty123/virtual-embryo-challenge-methods (no licence). community, 2026-09. https://github.com/wqty123/virtual-embryo-challenge-methods
[130] tubager/vec-t2 (no licence). community, 2026-11. https://github.com/tubager/vec-t2
[133] GitHub repository search for "virtual embryo". other, 2026-11. https://api.github.com/search/repositories?q=virtual+embryo&sort=updated&order=desc&per_page=60
[136] Arc Institute. Virtual Cell Challenge 2025 wrap-up. official, 2025-12. https://arcinstitute.org/news/virtual-cell-challenge-2025-wrap-up
[140] LaPierre. The real value of the Virtual Cell Challenge (Substack). community, 2025-12. https://nathanlapierre731428.substack.com/p/the-real-value-of-the-virtual-cell
[141] Szalay. How did we get a regression model to the top of the Virtual Cell Challenge (Turbine AI blog). secondary-industry, 2025-12. https://blog.turbine.ai/p/how-did-we-get-a-regression-model
[142] Chu, Stevens. Virtual Cell Challenge shenanigans (Substack). community, 2025-12. https://claudiachu99.substack.com/p/virtual-cell-challenge-shenanigans
[144] Luan et al. Centered cross-context perturbation transfer on the VCC 2025 benchmark (bioRxiv 2026.10.20.760514). academic, 2026-10. https://www.biorxiv.org/content/10.64898/2026.10.20.760514v1.full
[146] Arc Institute. The 2026 Virtual Cell Challenge. official, 2026-08. https://arcinstitute.org/news/virtual-cell-challenge-2026
[148] Lance et al. Longitudinal modality prediction, insights from a single-cell competition (bioRxiv 2026.02.24.707614). academic, 2026-02. https://www.biorxiv.org/content/10.64898/2026.02.24.707614v1.full
[149] Lance et al. Multimodal single cell data integration challenge (PMLR 176). academic, 2022-04. https://www.biorxiv.org/content/10.1101/2022.04.11.487796v1.full
[151] Open Problems Single-Cell Perturbations, top-4 write-up summaries (kaggle-solutions mirror, unchecked). community, 2023-12. https://kaggle-solutions-quartz.pages.dev/competitions/Open-Problems-–-Single-Cell-Perturbations
[155] CrunchDAO. Broad Obesity ML Competition documentation. official, 2026-07. https://docs.crunchdao.com/competitions/competitions/broad-obesity
[156] Sudbrack Zimmermann. Competition write-ups (GitHub README). community, 2026-05. https://github.com/rafaelsudbrackzimmermann/-MachineLearning
[158] Kaggle Chronicles: 15 years of competitions (arXiv 2511.06304). other, 2025-11. https://arxiv.org/html/2511.06304
[159] BirdCLEF+ 2025 working note on diversity-aware final selection (CEUR-WS Vol-4283; search summary only). academic, 2026-01. https://ceur-ws.org/Vol-4283/paper369.pdf
[160] Wang group. Third place in the Virtual Cell Challenge. community, 2025-12. https://jingshuw.org/post/vcc-2025/
