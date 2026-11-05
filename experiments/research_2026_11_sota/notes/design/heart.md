# T2:heart:test_extrap: solution design from the audited SOTA notes of 2026-11-06

> Note 2026-11-06 (main session, after the counter-review): the passages of this file that measure candidates against the captain's 209 are struck (marked below). 209 gates uploads only and is never a threshold or evidence for a try, an adoption or a method choice (docs/DATA_POLICY.md 2026-11-05; docs/PLAN.md section 8).

Board: T2:heart:test_extrap. Target: E12.5 3D MERFISH, 500-gene panel plus coordinates. Baseline: copy_last, which is
E10.5 with its own coordinates. Weight: 0.5 in Overall. Team Mithril, Human track. AS_OF 2026-11-06.

Written by the solution-design agent of the 2026-11-06 workflow for the captain's request: SOTA research turned into
solutions that raise our chance of a top-3 place (his target: a total above 209).

**Citation key.**
- [task-x n] is source [n] of the audited note task-x (`final/task-a.md` to `final/task-f.md`).
- "task-x F n" is that note's finding n, and "task-x L n" its candidate lever n.
- Team records are cited by repository path.

---

## 1. Plain answer

**The audited notes open nothing large on the heart.** Facet by facet:

- **Expression change.** One small refinement of B0's momentum can be counted (H-A). It rests on MERFISH run-structure
  evidence and is worth about +0.1 board points on the backtests. It is registered only if a target-free check (D-H1)
  first shows that its mechanism is large. Unconditional value: about +0.02 Overall.
- **Cell state.** Nothing. What moves it (composition, depth, cell numbers) is closed. No 2025-2026 generator beats
  real-cell copies at a held-out stage: the one paper with a copy baseline needed an atlas of about 10M cells for
  +0.06 to +0.10 soft Spearman across a spatial gap [task-d 8].
- **Neighbourhoods.** Nothing beyond what a better pseudobulk direction carries. In effect neighborhood_mmd is a
  normalised pseudobulk-shift metric (`experiments/research_2026_11/REPORT.md` section 6), and the arrangement levers
  are spent: NicheFlow's gain comes from niche sampling, in-distribution only [task-d 2].
- **Shape, where it can be separated from scope.** Nothing legal today. The one SOTA route is a legal late organ anchor
  applied at fixed RMS radius (H-E). A15's plain words close it, it also needs a ruling on the size and scope block, and
  it is priced at 0 ± 0.5 Overall. The scored cloud's shape is set by field of view and stack: at E9.5, 63% of cells
  are non-cardiac shell, and a core crop alone moves Dice to 0.79 [task-d F12; task-d 29].
- **States that emerge, proportions.** Closed: the captain's composition ruling of 2026-11-01, plus A15 for any external
  heart state.

**Recommendation.** Spend at most two days on the heart:
1. Run D-H1 (half a day, target-free; not a try).
2. Register H-A only if D-H1 passes its rule.
3. Give the captain D-H1's numbers with REPORT sheet E (H-B).
4. Decline H-C to H-F this season unless D-H1 or an organiser answer changes the picture.

The heart cannot carry the path to 209 (section 6). [Struck 2026-11-06: 209 is not a threshold for choices; information only.]

---

## 2. Where the board stands (committed records only)

**The entries.**
- s1 = H1: B0 on real E10.5 tissue at the count rule, with the pre-registered size. B0 moves each harmonised lineage by
  half its E9.5 -> E10.5 change in count space. Its coordinates are E10.5's, with the core grown by 1.58 times the
  observed log growth, xy grown and z held.
- s2 = H2 as the half trim (frozen, `docs/S2_RULES.md`).

**The B0 stand-in on the two counted backtests** (board points, 10 seeds; `experiments/harness/results/headroom/headroom.txt`):

| Group (weight) | B0 vs floor: 1-day / 1.75-day | B0 vs our copy | Points left, board (mean) | Points left, Overall |
|---|---|---|---|---|
| expression_change (0.25) | +0.38 / -0.46 | +0.77 / -0.03 | 12.54 | 6.27 |
| cell_state (0.25) | -0.32 / +0.44 | -0.20 / +0.43 (mmd_u -0.42 / +0.18) | 12.44 | 6.22 |
| shape_growth (0.25) | -0.07 / +0.19 | 0 / 0 (the copy's coordinates) | 12.44 | 6.22 |
| local_organisation (0.25) | +0.44 / +1.28 | -0.10 / +0.50 | 11.64 | 5.82 |
| Board | 50.43 / 51.45 | +0.47 (9 of 10) / +0.90 (10 of 10) | 49.06 | 24.53 |

**The estimate** (`experiments/gate/results/estimate_2026-11-05.json`):
- reading A: 53.85 (10th-90th percentile 48.48-60.12);
- first round alone (s1): 49.20;
- outlook: 52;
- heart SD: about 4.5 board points (REPORT section 1).

**Closed levers.**
- The post-exposure block: cell numbers, sizes and scope.
- The captain's ruling of 2026-11-01: the composition trend.
- A15: cross-modal input, and so any external direction or heart births.
- The ledger: the depth line (H3: -3.92 and +0.45), herding (H4), momentum mixtures and Potts sorting (task g).
- The frozen s2 rule.
- TRIES item 4: "Heart: no try".
- Running count: heart tries = 3.

---

## 3. What the notes add: the reference run's offset enters B0's direction with the wrong sign

**The model.** Let a MERFISH file's per-gene run offset be u_f(g), shared by every cell type in the file.
- B0 at test moves cells by about dp(g) = 1/2 [b(E9.5->E10.5)(g) + u_E10.5(g) - u_E9.5(g)].
- The scorer's truth change is about dt(g) = b(E10.5->E12.5)(g) + u_E12.5(g) - u_E10.5(g), plus composition and scope
  terms.
- So Cov(dp, dt) carries -1/2 Var(u_E10.5) on every gene.
- On the cell-level metrics, B0's run offset relative to the truth is 1.5 u_E10.5 - 0.5 u_E9.5 - u_E12.5. With
  independent runs of equal variance that is a squared distance of about 3.5 sigma_u^2 per gene, against 2 sigma_u^2 for
  a copy.

**Why the term is real and shared.**
- **MERFISH run effects are multiplicative and gene-structured** [task-c F11, F12; task-c 19, 20, 21]. Detection differs
  about 2-fold between sections of one mouse and falls with depth in the section. Dropped readout images remove 60-104
  genes at once in 31% of sections. Per-gene detection follows readout position and chemistry.
- **In vivo, comparisons with a separate sample are dominated by sample-level shifts** [task-c 13].
- **Metrics scored against a shared reference credit or penalise what a prediction shares with that reference**, most of
  all when the truth sits near it [task-c 1, 2]. The heart, at the floor, is that case.
- **Team record: run offsets are one per-gene factor shared by all cell types.** The Mab21l2 factor, measured in
  non-Mesp1 compartments, scored 80.5 against 50.8 when applied to all WT cells (`METHOD_SUMMARY.md`). Factor plus depth
  carry 0.78-0.80 of the decomposed Gata4 delta (REPORT section 4).
- **Team record: run factors are idiosyncratic.** Two runs' factors correlate 0.15 (TRIES N8).
- **The 1-day backtest shows the signature of added run distance.** B0's mmd_u is 0.42 board points below a copy's
  (skill 0.466 against 0.494), and its neighborhood_mmd 0.10 below.

**Limits.**
- B0 still beats a copy on both backtests (+0.47, +0.90). On pre-break increments, biology outweighs the term, together
  with the Poisson top-up's detection loss (the b0_smooth note in TRIES).
- The lineage-shared part also holds shared biology and the monotone depth drift: median counts per cell fall 233, 166,
  130, 75, 56 from E8.25_late to E10.5.
- The size of sigma_u^2 against the biology is unmeasured.

**What follows.**
- The expression change of the E12.5 truth is mostly the E12.5 run's own offsets, plus composition and scope. None of
  that can be predicted legally.
- The one predictable technical term, -u_E10.5, cannot be estimated at the series endpoint (leverage 0.80-0.91, the
  pruned item).
- What can be done is to stop adding +1/2 u_E10.5 (H-A, H-B, H-C).

---

## 4. Candidates

### H-A. Run-aware momentum: drop the lineage-shared part of B0's direction, gated by a target-free persistence check

**Facet and groups.**
- Facet: expression, the per-gene change within each harmonised lineage. No state, share, count, size or scope changes.
- Targets expression_change (12.54 board / 6.27 Overall left).
- Secondarily cell_state (12.44 / 6.22) and local_organisation (11.64 / 5.82): both see any common run shift, and the
  neighbourhood metric behaves like a pseudobulk-shift metric.

**Mechanism.**
- B0's factor is log f_lin(g) = 1/2 [log(m_ref + 0.5) - log(m_start + 0.5)], over mean CP10k per merged lineage
  (`vec.t1_scope.momentum_shift`, called by `vec.harness._b0_t2`).
- H-A replaces each tracked lineage's log f by its deviation from the across-lineage median of the tracked lineages (a
  lineage contrast). This cancels any per-gene term every lineage shares, which is where u_ref - u_start sits.
- It keeps lineage-specific change, which is where the one carried-over signal was: "per-type momentum carries over
  only for cardiomyocytes" (TRIES, pruned list).
- Applied by `apply_factors` with keep_depth, at B0's damping 0.5. "Untracked" groups stay flat, as in B0.

**SOTA basis.**
- MERFISH run structure [task-c 19, 20, 21].
- Sample-level shifts in separate-sample comparisons [task-c 13].
- Shared-reference carriage [task-c 1, 2].
- A shared template recovers only the conserved part across contexts [task-c 3].
- [task-a 5] (DCI): in scRNA, only changes that agree across cell types forecast better than carry-forward (MASE 0.78
  against above 1 below DCI 0.2). In MERFISH the cross-lineage part is also where the run term sits, so its persistence
  must be measured, not assumed. D-H1 measures it.

**New against the pruned list and the 2026-11-05 report.**
- Nearest pruned item: "Heart: removing the reference file's own transient per-gene deviation", pruned because the
  trend estimator's leverage at the reference file is 0.80-0.91. H-A estimates no deviation: it cancels the shared term
  by construction, with no endpoint leverage.
- Distinct from "T2: a heart momentum mixture" (task g), which mixes shift sizes and was a wash.
- Distinct from REPORT sheet E's hold, which has no counted check. H-A is countable because the reference-run term sits
  in the backtests' momentum exactly as it does at test.
- Distinct from "Heart: two gap-anatomy diagnostics": D-H1 never reads E10.5.

**Legality.**
- Released heart MERFISH only: E8.25_late, E8.75 and E9.5 are training files; E8.5 and E10.5 are unlocked validation
  files used under A20. No external data, pretrained model or atlas.
- Genotype rules and the KO ruling are not engaged.
- Post-exposure block respected: expression only.
- A15 not engaged: same task, same modality.
- Frozen rules: `docs/S2_RULES.md` lets a change "in the heart pipeline that every heart entry runs" carry over to the s2
  trim. Its hedged parameters (scope and size) are untouched, but the captain is told.
- Before the run: a TRIES row, and the rule as a plain summary to the captain (`docs/PLAN.md` decision log
  2026-11-04 (d)).
- If the T2 output changes, the golden table in `tests/test_core.py` is updated in the same commit.
- No ruling is needed.

**Counted evidence.**
- Backtests: `T2h/<=E9.5->E10.5` (momentum E8.75 -> E9.5) and `T2h/<=E8.75->E10.5` (momentum E8.5 -> E8.75), 10 paired
  seeds against B0, also reported against floor and copy.
- Kept only under the adoption rule on both (paired mean above +0.3 board points, ahead in at least 8 of 10 seeds, no
  loss above 0.3), with the gain on expression_change and neighbourhoods as predicted, not on variogram alone.

Why the backtests' known features (TRIES rule 2) would not explain a gain:
- **The term is structural.** It is present in both backtests (u_E9.5 and u_E8.75 sit in their momentum and, with the
  opposite sign, in their dt) and at test (u_E10.5). It is a property of momentum that ends at the reference, not of
  one backtest.
- **The 1-day truth crosses the dissection change** (`configs/pseudo_tasks.yaml`). A removal of the part shared across
  lineages cannot imitate a scope or composition change.
- **Both backtests share the E10.5 truth,** so their verdicts are correlated (u_E10.5 is common noise). This is written
  into the row as a weakness, not hidden.
- **B0's 1.75-day gain over a copy lies on the top-up's detection-loss line.** H-A shrinks the factors, so a variogram
  loss there is expected and is not evidence against the direction. The row reports it separately.

**Pre-check D-H1** (target-free; not a try, as N1 and N2):
- **Reads:** E8.25_late, E8.5, E8.75 and E9.5 only, duplicates dropped. Never E10.5, no score, no cell count printed.
- **Transitions:** A = E8.25_late -> E8.5, B = E8.5 -> E8.75, C = E8.75 -> E9.5.
- **Per transition:**
  - per tracked lineage (at least 20 cells at both stages), Delta_lin(g) by B0's formula;
  - the shared part s(g), the median over lineages;
  - the shared share of energy, sum s^2 / mean_lin sum Delta_lin^2.
- **Persistence:** the Spearman correlation of s between consecutive transitions (A,B) and (B,C), partialled for the
  rank of mean expression (as de_direction partials pb(ref)).
- **Placebo:** genes permuted within expression deciles, seed fixed in the row.
- **Rule, fixed before the run.** Register H-A only if:
  - the shared share is at least 0.3 in each transition; and
  - both partial correlations are at most -0.2 and outside the placebo's 95% band, which means the shared part is
    run-dominated (the middle run's offset flips sign between consecutive transitions).
- **Otherwise** prune H-A (B0 stays).
- **Either way,** the numbers go to the captain with sheet E.

**Expected gain.**
- On the backtests: about +0.1 board points (range -0.5 to +0.8). The anchor is B0's own margins: its whole momentum is
  worth +0.77 and -0.03 board points on the DE group over a copy, and costs 0.42 on mmd_u on the 1-day backtest. H-A
  could at best recover that and keep part of the DE gain.
- The SOTA effect sizes are discounted for transfer: DCI's 22% error reduction comes from scRNA; the carriage sizes from
  cell lines.
- P(adoption): about 0.1, or about 0.15 if D-H1 passes. The "on every backtest" requirement and the 1.75-day backtest's
  top-up side effect make a double pass unlikely.
- If adopted: +0.2 to +0.5 board at test, i.e. +0.1 to +0.25 Overall.
- Unconditional: about +0.02 Overall.
- Upside: none measurable. It removes a noise term, and the heart's spread comes from size and scope.

**Effort and DGX need.**
- D-H1: half a day, minutes of DGX CPU.
- The try: about one day (a T2-head flag, the TRIES row, 2 backtests x 10 seeds on DGX CPU, the golden table).

**Risk.**
- The shared part may carry persistent biology or the depth drift; D-H1 is the guard.
- Variants inflate lucky-fit risk, so only this one form is registered.

**First step.** Commit D-H1's script and rule as a TRIES diagnostic row, send the captain the plain summary, and run it
on the DGX.

### H-B. Hold, or cardiomyocyte-only momentum, at test (the captain's ruling; REPORT sheet E revisited)

**Facet and groups.** Expression; the same groups as H-A.

**Mechanism.**
- At test, B0's increment is E9.5 -> E10.5. It crosses the dissection change (embryo portion to dissected heart), and
  its run term is the reference's own (+1/2 u_E10.5).
- The counted backtests scored only pre-break increments (task-f L2 in `experiments/research_2026_11/notes/`; REPORT
  section 6).
- H-B holds E10.5's expression (damping 0), or keeps momentum only for CM, the one lineage where per-type momentum
  carried over.

**SOTA basis.**
- The run-structure items of H-A [task-c 1, 2, 13, 20].
- [task-a 5]: lineage-specific changes do not beat carry-forward.
- [task-f 13]: a single held-out-time validation ranks models only moderately (Spearman 0.72-0.87 to private rank), so
  pre-break backtests say little about the across-break increment.

**New against the pruned list and the report.**
- This is sheet E and task-f L2, declined at +0.05.
- New: a mechanism that does not depend on the break (the sign of the reference run's term), and D-H1's number as its
  deciding input.

**Legality.**
- Expression only, no data issue.
- No counted check of the form used at test, so it enters only by the captain's ruling, as R1 did.
- The heart s2 rebuilds from s1, so the change carries over (the captain is told).

**Counted evidence.**
- None in the form used at test: a captain's bet.
- A hold on the backtests is a copy, which lost to B0 by 0.47 and 0.90; that measures pre-break momentum, not the test's.

**Expected gain.**
- The committed price is ±0.5 board points (±0.25 Overall), with an expected value of +0.1 board (+0.05 Overall).
- If D-H1 finds the shared part run-dominated: about +0.1 to +0.2 board. If it finds it persistent: about -0.1 board.
- Upside: a slightly narrower heart; no effect on our 90th-95th percentile.

**Effort.** Half a day plus the ruling.

**Risk.** It discards real CM maturation signal and B0's detection-loss side effect, which may help again at test,
since depth falls at every released stage.

**First step.** None beyond D-H1. The captain decides sheet E with D-H1's numbers before the 2026-11-27 freeze.

### H-C. Readout-structured run model of the momentum (codebook groups from co-dropout)

**Facet and groups.** Expression; the same groups as H-A.

**Mechanism.**
- Genes whose barcodes share a readout round lose counts together when that round's image is dim or dropped
  [task-c 19, 20].
- Infer gene-to-readout groups target-free, from co-dropout inside the released files: sections where a block of
  60-104 genes falls to near zero together.
- Fit the lineage-shared part of B0's Delta on that low-dimensional readout basis (16-24 coefficients over 500 genes),
  and remove only the fitted readout-structured part.
- The rest of the shared part (shared biology, depth drift) stays. It is a sharper form of H-A's full contrast.

**SOTA basis.**
- [task-c 19, 20, 21].
- [task-c 3]: 92-98% of a shared template lies in the top 100 control PCs.
- [task-d 6]: SINTER3D absorbs per-section nuisance terms with a section embedding and offset.

**New against the pruned list.**
- Nearest pruned items: "T3: a fluctuation-response (wild-type covariance) programme" (biological axes) and task-c L4
  (technical axes, for T3).
- On the heart, no pruned item uses the codebook structure.
- It competes with H-A for the same try.

**Legality.**
- As H-A.
- MerQuaCo code, if used, needs a licence check and an `EXTERNAL_DATA.md` row first (task-c audit).

**Counted evidence.**
- Both heart backtests, as H-A.
- Pre-check, rule fixed first: co-dropout blocks must exist in the released heart files, and the readout basis must
  explain the shared part of the pre-E10.5 transitions better than a permuted-membership control. Otherwise prune.

**Expected gain.**
- Backtests: about +0.05 board (range -0.3 to +0.5).
- P(adoption): at most 0.05.
- Unconditional: about +0.01 Overall.
- Upside: none.

**Effort.** 3-5 days on DGX CPU.

**Risk.** Our files may have no dropped-image sections. Runs of one protocol may differ along other axes (buffer age, z
decay, section depth) [task-c 20]. The codebook is unknown.

**First step.** Only if D-H1 passes: a target-free count of co-dropout blocks per released heart file.

### H-D. Run-noise shrinkage from the knockout files (needs the captain's reading of A15)

**Facet and groups.** Expression; the same groups as H-A.

**Mechanism.**
- The released knockout files are independent MERFISH runs on the same panel (TRIES N6: Gata4 on slide M164, Mab21l2 on
  M125).
- Their non-Mesp1 compartments against the matched WT are direct samples of between-run per-gene differences.
- Estimate each gene's run-noise variance (smoothed across genes by expression level), then shrink each gene's momentum
  by its reliability. The shrinkage changes ranks: run-noisy genes move toward zero.

**SOTA basis.**
- [task-c 13]: unperturbed cells of the same animal isolate the sample shift.
- [task-c 1]: release a correction only against an independent control.
- [task-c 10]: TRADE, estimation-noise-aware impact.

**New against the pruned list.**
- Nearest pruned items: "T1: a noise-aware, subtype-resolved Feng direction" (split halves found nothing to shrink) and
  "T3: a two-scenario factor exponent".
- New here: the noise is run-level, which split halves cannot see.

**Legality.**
- It uses T3 files for T2: same modality, other task.
- A15's working answer ("no new cross-modal input") does not bar it in its words, but the A15 row says "keep modalities
  separate per task". So it needs the captain's reading.
- Precedent: the T3 rules' WT series already uses heart E8.5 (component R).
- To keep the heart free of Gata4-derived numbers (the reason T3's s2 is Gata4-free), use Mab21l2 against WT E9.5 only.
  That leaves one run difference; with Gata4 there are three, two of them on one slide.
- No composition change, no external data.

**Counted evidence.** Both heart backtests, as H-A; the weights are target-free.

**Expected gain.**
- About 0 to +0.1 board on the backtests. One to three runs estimate per-gene variance poorly, and two run factors
  correlate 0.15 (N8).
- Unconditional: about 0 Overall.
- Upside: none.

**Effort.** 1-2 days after the ruling.

**Risk.** The noise estimate is too poor, and it adds a dependence on A20 and A15 readings to a fourth board.

**First step.** None now. It is an add-on to H-A only, after the captain's A15 reading.

### H-E (bold). Fixed-RMS core-shape continuation toward a legal late organ anchor

**Facet and groups.** Spatial organisation (shape). Targets shape_growth through d2_shape and occupancy_dice (12.44
board / 6.22 Overall); scale untouched.

**Mechanism** [task-d L2; task-d 13, 17, 18, 19, 20]:
- Register wild-type heart volumes from legal late anchors to the MERFISH core (CM plus endocardium): DMDD HREM embryos
  staged S23- or S23, the E15.5 micro-CT average, or the EMA TS23-TS24 MRI models.
- Interpolate the core's principal-axis ratios, or a deformation field, from E9.5 and E10.5 toward the anchor.
- Apply it to core cells only, at the pre-registered RMS radii and counts.

**New against the pruned list.**
- No pruned item repeats it.
- It overlaps the closed "heart's core-focused s2 (scope; a frozen rule)" and TRIES item 4, "Heart: no try".
- New: a legal late anchor for the core's shape.

**Legality: closed today.**
- **A15.** On its plain words, an imaging shape applied to MERFISH coordinates is a new cross-modal input (task-d audit;
  precedent: Feng's direction on the heart panel, pruned under A15).
- **Size and scope block.** A fixed-RMS shape change needs the captain's ruling that it lies outside the block.
- **If opened:**
  - wild-type embryos only, since DMDD also images 58 lethal mutant lines;
  - the staging rule's S23 filter;
  - an `EXTERNAL_DATA.md` row, with DMDD's data terms (still unchecked);
  - a step in s1 alone sends the heart s2 back to the captain.
- The windows hold: the anchors are E14.5 or later, and the E11.5-E14.0 OPT atlas is barred.

**Counted evidence.**
- task-d names both heart backtests, with the anchor at t = 1/5 and the matching fraction.
- Rule 2 bites. The scored cloud is set by field of view and stack: at E10.5 the core is 23% of cells, and a shell
  removal or core crop alone moves Dice to 0.86 or 0.79. A backtest gain would read E10.5's acquisition, not the organ.
- So in effect no clean counted check exists: a captain's bet.

**Expected gain.**
- 0 (range -1 to +1 board; -0.5 to +0.5 Overall) [task-d L2].
- The public record points the same way (selection-on-validation evidence only): every geometry deformation tried
  publicly on extrapolation lost points (`docs/PLAN.md` 4.3), and four morphology-based heart shape transfers failed or
  hurt [task-e 25].
- Upside: a little added variance in the shape group, worth at most +0.1 to +0.2 Overall at our 95th percentile.

**Effort.** 1.5-2 weeks (segmentation, cross-modal registration). Not feasible with its rulings before the 2026-11-27
freeze.

**Risk.** A wrong sign; scope dominates the scored shape.

**First step.** None recommended. If the captain wants it anyway, the first step is his reading of A15 (does
"cross-modal" cover non-omics imaging?), not code.

### H-F (bold, organiser question). Do external spatial and organ-imaging data count as cross-modal on T2?

**Facet and groups.** Access to expression and shape data; all four groups.

**The question.** One paragraph to the organisers:
- (a) Do other labs' imaging spatial transcriptomics (MERFISH, Xenium) or 3D organ imaging (HREM, micro-CT) count as
  cross-modal external data for the T2 heart board?
- (b) May released MERFISH files cross tasks within the modality, for example T3 knockouts for the T2 heart?

**Basis.**
- [task-d], Gaps: "if spatial-to-spatial is not cross-modal, late spatial heart data could open a heart-direction
  lever".
- [task-b] M6, M7: Hu 2026 Visium and Slide-seq at E14.5 or later, deconvolved against Feng [task-b 6]; postnatal
  Xenium [task-b 7].

**New.** A15's question (18 October email, question 1) covers released data across tasks, not external spatial or
imaging data.

**Legality.** An organiser answer is the legitimate route. It would open H-D's and H-E's data routes.

**Counted evidence.** None; it is a question.

**Expected value for the heart: about 0.**
- No legal mouse heart MERFISH exists at or after E14.0 ([task-d], Gaps).
- Visium and Slide-seq are sequencing-based, so cross-modal in any case, and they inherit Feng's reference [task-b M6].
- Xenium is postnatal only, with no early anchor on the same platform, and per-gene ratios between chemistries are not
  uniform [task-c 21].
- The target-free Feng-vs-MERFISH concordance was already null (cardiomyocyte Spearman 0.05-0.10; TRIES).
- The organisers have answered no team rule question since 2026-10-16.

**Effort and risk.** One hour; no risk.

**Recommendation.** Do not send it for the heart's sake. Append it to an email sent for another reason.

### Summary

| Candidate | Route | Counted check | E[gain], Overall | Range if done, Overall | Effort |
|---|---|---|---|---|---|
| H-A lineage-contrast momentum | D-H1, then one TRIES try | both heart backtests | about +0.02 (unconditional) | -0.25 to +0.4 | 1.5-2 days |
| H-B hold or CM-only at test | captain's ruling (sheet E) | none for the test form | about +0.05 | -0.25 to +0.25 | 0.5 day |
| H-C readout-structured removal | only if D-H1 passes; competes with H-A | both heart backtests | about +0.01 | -0.15 to +0.25 | 3-5 days |
| H-D run-noise shrinkage | captain's A15 reading, then a try | both heart backtests | about 0 | 0 to +0.05 | 1-2 days |
| H-E late organ shape (bold) | A15 plus a size/scope ruling | none clean | 0 | -0.5 to +0.5 | 1.5-2 weeks |
| H-F organiser question (bold) | organiser answer | none | about 0 for the heart | - | 1 hour |

---

## 5. Considered and not proposed

- **Learned spatiotemporal flows and niche generators** (ContextFlow, stVCR, TP-DATE, NicheFlow [task-d 1-4]).
  - No paper reports a gain over copy-nearest at a held-out stage.
  - Generated states change the heart's composition, which is closed.
  - GPU non-determinism breaks the freeze-morning hash check.
  - Expected at most 0 (task-d (a)); it repeats the latent-transport prunes.
- **Isotropic resampling and continuous 3D fields** (isoST, SINTER3D, UniST [task-d 6, 7, 10]). Harmful: a 1.5x z
  mismatch costs 12.6-12.8 board points [task-d F13]. Keep native z.
- **Heart composition, births and within-core shares.** Closed by the captain's ruling of 2026-11-01 and by A15. Feng's
  within-core shares also had the wrong sign (TRIES).
- **A blood-share prior** (the heart analogue of task-b L1). The heart MERFISH vocabulary has no erythroid category,
  only haemato-endothelial cells (`configs/lineage_map.yaml`; `experiments/labels/results/E10.5__vs__E9.5.csv`), so
  there is nothing to move.
- **An external expression direction on the panel** (Feng, JAX, Visium). A15 closes it, and the concordance check was
  null.
- **Depth, by thinning or as a normalisation trick** [task-c 22, 23]. The ledger closed the line (H3: -3.92 and +0.45),
  and nothing new reopens it. Any persistent depth signature would show up in D-H1's shared part.
- **Count generators** (count-FM, scLDM [task-a 18, 19]). No held-out count-realism evidence exists, and real cells
  already carry E10.5's sparsity.
- **Herding, momentum mixtures, Potts sorting, span-matched shuffles, tail decoders.** All pruned. N5 found 1-3 of the
  heart's bottom 40 genes infeasible.
- **A reseeded s1 in the s2 slot** [task-f L2]. The heart s2 is the scope hedge, so the post-exposure block closes it.
- **Reopening scope, size or cell numbers.** I do not propose reopening the block: it exists because test-board limits
  were seen.

---

## 6. Path to 209 for the heart

> Struck 2026-11-06 (main session, after the counter-review): this passage measures candidates against the captain's 209, which gates uploads only and is never a threshold or evidence for a try, an adoption or a method choice (docs/DATA_POLICY.md 2026-11-05; docs/PLAN.md section 8). Information only; nothing was chosen by it.

**What 209 takes.** Reading A totals 191.5. Reaching 209 needs +17.5 Overall, at or beyond reading A's 99th
percentile: its 90th-99th percentiles run about 204-208 (REPORT section 1).

**What the heart can supply.**
- The heart counts 0.5. Its 10th-90th band is 48.5-60.1 (SD about 4.5 board points).
- At its 90th percentile it adds about +3.1 Overall over its mean.

**The share I assume.**
- The heart at 58-60, i.e. +2 to +3 Overall over its mean.
- That needs the pre-registered size to land and the s2 trim to match a more heart-focused E12.5.
- On reading A's band this happens with probability about 0.2 (about 0.1 for 60).
- The other +14.5 to +15.5 must come from T1 and T3, landing near the tops of their bands (T1 about 73, T3 about 75)
  at the same time. The REPORT puts our top-3 chance at a few percent at most, and under 1% at our 95th percentile
  under the measurement-based T3 payoffs.

**What the candidates change.** Nothing above moves the heart's mean by more than about +0.05 Overall, or its tail by
more than about +0.2 (H-E, closed).

**Recommendation.** The heart's role is to be built and uploaded cleanly at the freeze, its rehearsed builds unchanged.
Its try budget should go elsewhere: D-H1, and H-A only if D-H1 passes, is the most it earns. Its coherence with T1's s2
is already in place: both s2s bet on a more chamber-dominated E12.5.

---

## 7. Compliance of this design

**Read.**
- The six audited notes and the brief.
- From the repository: `AGENTS.md`, `docs/CHALLENGE.md`, `METHOD_SUMMARY.md`, `experiments/harness/TRIES.md`,
  `docs/PLAN.md` (sections 1-3, 4.3 and 8, and the decision log of 2026-11-04 to 2026-11-06), `docs/DATA_POLICY.md`
  (windows, Scores, the notes of 2026-11-04), `docs/S2_RULES.md`, `ASSUMPTIONS.md`, `docs/ORGANIZER_QUESTIONS.md`,
  `experiments/research_2026_11/REPORT.md`, `results/upload_gate.json` and
  `experiments/gate/results/estimate_2026-11-05.json`.
- For mechanics and numbers:
  - `experiments/harness/results/headroom/headroom.txt` (heart rows) and `configs/pseudo_tasks.yaml` (heart tasks);
  - `configs/lineage_map.yaml` and `experiments/labels/results/E10.5__vs__E9.5.csv` (released and unlocked label
    counts, not a test target);
  - the heads of two `experiments/harmonise/results/` tables (committed released-data summaries);
  - `src/vec/harness.py` (`_b0_t2`) and `src/vec/t1_scope.py` (`momentum_shift`);
  - the heart lines of `experiments/research_2026_11/notes/task-f.md` and `task-g.md`;
  - the claim list of `experiments/coords/REPORT.md`.

**Not opened.** `boards/`, `data/`, `runs/`, `raw/`, `screened/`, `docs/research/`, anything under `~/.claude-account` or
`~/.claude`, and the web.

**Nothing else done.**
- No computation and no repository edit; this file is the only output.
- No hidden-target property, per-board leaderboard figure or in-window developmental claim is used.
- The relayed user request (sync the DGX and this machine, no unpushed commits) applies before the workflow. This design
  step made no git or DGX change.
