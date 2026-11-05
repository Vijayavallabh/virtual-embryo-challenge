# T2:embryo:test_interp (E7.75 MERFISH): solution design from the SOTA round of 2026-11-06

> Note 2026-11-06 (main session, after the counter-review): the passages of this file that measure candidates against the captain's 209 are struck (marked below). 209 gates uploads only and is never a threshold or evidence for a try, an adoption or a method choice (docs/DATA_POLICY.md 2026-11-05; docs/PLAN.md section 8).

Board: T2:embryo:test_interp, counts 0.5 in Overall. Published baseline copy_last (E7.5 with its coordinates, scores
50). Entries: s1 = fgw, s2 = fgw_shape0 (frozen, `docs/S2_RULES.md`). Written by the AI assistant as the board's
solution designer in the SOTA workflow. AS_OF 2026-11-06. Design only: nothing was run and no repository file was
changed.

## 1. Bottom line

- **The stop on s1 stands.** None of the six audited notes offers a counted check for any s1 change.
  - The 2025-2026 spatiotemporal generators compare only against other learned flows.
  - The one study with a copy-nearest baseline needed an atlas of about 10M cells to gain +0.06 to +0.10 soft
    Spearman, across a spatial gap, not a temporal one [task-d 8].
  - Everything else repeats a pruned item (section 5).
- **What this round changes is the one open embryo question: the s2 backup (REPORT sheet B).** Three derivations,
  none of which reads a target, sharpen it (section 3):
  - (a) The mixing fraction cannot hedge the expression group at all.
  - (b) A legal staging-noise prior puts P(E7.75 sits at t' < 0.375) at about 0.16-0.40, depending on how tightly
    the lab stages. That range runs from sheet B's break-even (0.15-0.2) to about twice it.
  - (c) Sheet B's second reason, the depth hedge, is weak. The usable released interior stages sat within 5% of
    log-linear depth, so the t = 1/2 mixture's log-mean depth (115) is the central forecast.
- **Recommended:**
  - C1: a target-free pseudo-truth table on held-out section halves, then the captain's ruling on sheet B before the
    embryo's first upload.
  - C2: one organiser question that would settle both the backup and the geometry hedge.
  - Optional: C4. The bold candidates C3 and C5 are not recommended as they stand.
- **Value:** about +0.15 Overall in expectation (range -0.5 to +1.0). In the world where E7.75 sits near E7.5
  (roughly one chance in four), it adds +2 to +4 Overall. That is a world in which teams that interpolate at 1/2 lose
  together.
- **Path to 209:** [struck 2026-11-06: information only; 209 chooses nothing] the embryo can supply about a tenth to a fifth of the gap, and only in its favourable world; in
  expectation it supplies under 1% (section 6).

## 2. Where the board stands (committed results only)

| Seed means over 10 seeds (`results/emb_fgw_v1`, `emb_shape_v1`) | E7.5 backtest (t = 1/3) | E7.25 backtest (t = 2/3) |
|---|---|---|
| s1 fgw, score (seed SD) | 61.12 (0.83) | 58.78 (2.20) |
| s2 fgw_shape0, score (seed SD) | 62.95 (0.76) | 53.67 (1.94) |
| s1 points left: expression / cell state / shape / neighbourhoods | 8.96 / 8.29 / 12.63 / 9.01 | 7.69 / 14.72 / 8.86 / 9.95 |
| s2 points left: shape / neighbourhoods | 11.16 / 8.64 | 12.99 / 10.94 |

- **Estimate** (`experiments/gate/estimate.py`): better-of-two 60.87 (10th-90th percentile 56.9-63.5), first round
  59.95. That is 30.4 Overall. The board has 20.02 Overall points left in all, the smallest gap of the four boards.
- **The two uploads share one mixture.**
  - Their expression and cell-state skills are therefore identical, and they differ only in geometry.
  - s1's moved cloud is flatter than either bracket in every test draw: its thinnest axis is 0.21-0.25 of the RMS
    radius, against 0.36 for E8.0 and 0.53 for E7.5 (N1, test panel).
  - s2 takes E7.5's geometry.
- **Closed or fixed:**
  - cell number, size (PCHIP) and scope, by the post-exposure block;
  - the depth, blend, herding and coupling lines, by the ledger;
  - the clock warp, pruned;
  - counted evidence for any mixing-fraction change: E7.25 and E7.5 cannot count for one (TRIES, DGX campaign entry
    "Embryo.");
  - s2, which is frozen and moves only by the captain's ruling;
  - cross-modal inputs, under A15.

## 3. Three findings behind the candidates (derivations; no target read)

1. **Rank invariance.**
   - For a mixture with a fraction t of E8.0 cells:
     dp = t·[pb(E8.0 draw) − pb(E7.5 draw)] + [pb(E7.5 draw) − pb(reference copy)].
   - de_score (the top n_up and bottom n_dn genes by dp) and de_direction (a rank partial correlation) are
     invariant to scaling dp. So every fraction t > 0 predicts the same gene ranking. Only its signal-to-noise
     against the scorer's 10% reference copy changes, and that halves at t = 1/4.
   - The expression group (25 of 100 board points) therefore cannot be hedged by the fraction. Cell state and
     neighbourhoods (50 points) can. This matches the REPORT's "unhedged on half the board".
   - It also adds a small cost of a lower fraction in every world, which the REPORT's toy does not price.
2. **Staging noise.** Sheet B's deciding prior is P(t' < 0.375), where t' is E7.75's position between its brackets.
   - Legal sources show that a single embryo can stray from its gestational label by a large fraction of a day:
     - E14.5-harvested embryos span six forelimb sub-stages, S21 to S23 [task-d 17];
     - whole-embryo atlases stage by somite count in 2-h bins [task-b 7], or by each embryo's morphology with a
       photo [task-b 3], because the gestational label is too coarse;
     - litter mis-staging can masquerade as mutant delay, and somite-matched individuals replace uterus-mates
       [task-c 16][task-c 17] (method level only; no in-window staging account is used).
   - Model: independent per-embryo staging noise σ (days) on E7.5, E7.75 and E8.0. Then t' < 0.375 exactly when
     a(7.75) − 0.625·a(7.5) − 0.375·a(8.0) < 0. That quantity is normal with mean 0.0625 day and SD 1.24σ.
   - Result: P = 0.16, 0.31, 0.37 and 0.40 at σ = 0.05, 0.10, 0.15 and 0.20 day. Equivalently, the SD of t' is
     about 2.45σ.
   - No backtest position enters this prior (REPORT section 6 bars them).
   - The quarter-day labels (E6.75, E7.25, E8.25_late) may mean morphological staging, which would put σ near the
     low end. Only the organisers can say (C2).
3. **Depth.**
   - From E7.25 on, the released embryo series rises steadily in log depth: medians 35, 62 and 212, or 2.29 and 2.46
     per day over E7.25 → E7.5 → E8.0.
   - The two usable interior stages sat within 5% of log-linear interpolation between their brackets: E7.5 at 62
     against 64, and heart E8.5 at 166 against 174 (TRIES rows E7 and H3). E7.25's 0.71x shortfall is part of the
     known feature that is barred as mixing-fraction evidence, so I set it aside.
   - So the t = 1/2 mixture's log-mean depth, 115, is the central forecast for E7.75. A t = 1/4 backup (log-mean 84)
     hedges the logged "E7.5-like depth" risk only in a tail, which I put at about 0.1.
   - The tail is not zero: one run can land off trend. MerQuaCo found 2-fold section-to-section detection
     differences and readout-group losses [task-c 20], and protocol and chemistry move counts [task-c 19][task-c 21].

## 4. Candidates

### C1. Lean-early backup: the embryo s2 at mixing fraction 1/4 (REPORT sheet B), re-priced. Recommended for a ruling once its table is in

- **Facet and targets.**
  - Facet: proportions (which bracket's real cells, and in what share).
  - Targets cell state (11.50 board points left on the mean; 8.29 and 14.72) and neighbourhoods (8.64 and 10.94 for
    s2).
  - Expression change (8.96 and 7.69 left) is out of reach (finding 1).
- **Mechanism.**
  - s2 keeps fgw's coupling, E7.5's geometry (t_pos = 0), the size rule and the count rule, with 3/4 E7.5 cells and
    1/4 E8.0 cells. `vec.t2_build.interpolate_clouds` already takes t; the builder needs one flag.
  - The skill map's floor distance shrinks when the truth sits near copy_last. A t = 1/2 prediction therefore loses
    most when E7.75 sits near E7.5: skill(t; t') = t'^2 / (t'^2 + (t − t')^2) (REPORT section 6). A lean-early
    backup covers that side, and t = 1/2 stays robust on the far side.
  - Today's s2 shares s1's mixture, so the near-E7.5 world is unhedged on the 50 mixture-sensitive points. Across
    two uploads, a two-point pair gains 0.08-0.11 skill under kernel scores. T1's frozen pair has that form; the
    embryo pair does not (REPORT section 2).
- **SOTA basis.**
  - The staging-noise prior of finding 2 [task-d 17][task-b 7][task-b 3][task-c 16][task-c 17].
  - Contest records: ranks are noisy under a time shift and only the top rank is bootstrap-robust [task-f 13]; a
    winner moved up to 8 ranks on re-runs [task-f 14]; a decorrelated second pick was kept [task-f 25]; picks tuned
    to the public board failed [task-f 24].
  - Copy-nearest is a strong predictor [task-d 8].
  - A laggard gains by cutting correlation with the field (REPORT section 1). If, as seems likely, most teams
    interpolate near 1/2, the near-E7.5 world hits them together.
- **What is new.**
  - Nearest pruned items: "a developmental-clock warp of the bracket fraction t" (nothing says which way E7.75
    sits) and the DGX campaign's "Embryo." entry (E7.25 and E7.5 cannot count).
  - C1 claims no direction. It hedges symmetric staging noise through the skill map's asymmetry.
  - Against REPORT sheet B:
    - (a) the prior comes from legal staging evidence, not from backtest positions;
    - (b) the expression group is shown to be unhedgeable;
    - (c) the depth hedge, sheet B's tiebreaker, is weak (finding 3);
    - (d) a table on held-out halves replaces the toy;
    - (e) the decision is timed before the embryo's first upload.
- **Legality.**
  - Inputs are E7.5 (unlocked, A20) and E8.0 (training), as for today's s2. No external data enters the build.
  - The prior uses no in-window account (class W) and no backtest position.
  - Cell numbers, size and scope are unchanged (post-exposure block). The embryo's composition is not among the
    block's closed items.
  - No cross-modal input (A15). No genotype or knockout-ruling question arises.
  - It amends the frozen embryo s2 rule after the fact, as R2 would have. That needs the captain's own ruling, in its
    own commit with his words and the file's new sha256 (`docs/S2_RULES.md`).
  - It is best ruled before the embryo's first upload (from 2026-11-27), so that no test score is public when the
    rule changes.
- **Counted evidence.** None exists and none can: E7.25 and E7.5 cannot count for a mixing-fraction change, and the
  s2 rule is frozen. It is priced as a captain's bet; the target-free table below decides whether it goes to him.
- **Expected gain.**
  - Sheet B's loss if wrong (about 1.5 board points) and break-even (0.15-0.2) imply a gain of about 6-8.5 board
    points when right.
  - Discounting for the unhedgeable expression group and for mixtures standing in for real intermediates:
    G ≈ 4-8.5 and L ≈ 1-2 board points.
  - With P = 0.16-0.37, E = P·G − (1 − P)·L ≈ +0.3 board points (range −1.0 to +2.0), or **+0.15 Overall (range
    −0.5 to +1.0)**.
  - In the near-E7.5 world (P ≈ 0.16-0.37) it adds +2 to +4 Overall over today's pair. Elsewhere it costs up to
    about 1 Overall.
  - Upside: its gain sits in a world independent of T1 and T3, so at our 90th-95th total percentile it adds about
    its mean (+0.1 to +0.3 Overall). Under the rank view, where the field shares the near-E7.5 loss, it is worth more
    than its mean.
- **Effort and risk.**
  - Effort: about 1 day of code (the builder flag and the table script) and 1-2 hours of DGX CPU at 16 threads. No
    GPU.
  - Risk 1: pseudo-truths built from bracket mixtures favour mixture predictions at a matching t. Held-out halves
    mitigate this, and C4 checks real intermediates.
  - Risk 2: it is an after-the-fact amendment, which the captain declined once, for R2. R2 rested on selected-on
    backtests; this rests on target-free arithmetic and legal priors.
- **First step.** Register an N row: a diagnostic, not a try, on N1's precedent. It reads no target and no score,
  and prints no cell count.
  1. Use the test pair only (E7.5, E8.0). Split each bracket by alternating sections into halves A and B.
  2. From half A, build pseudo-truths as mixtures at t' = 0.1, 0.2, ..., 0.9, in two geometries: E7.5's (t_pos = 0)
     and the aligned interpolation (t_pos = t').
  3. From half B, build s1 (fgw), today's s2 (fgw_shape0) and the backup (t = 1/4, t_pos = 0), with B's E7.5 as the
     floor.
  4. Score with the calibrated emulator over 10 seeds.
  5. Report each pair's expected better-of-two [task-f mapping 7] under t' ~ N(0.5, sd), for sd = 0.12, 0.25 and
     0.37 (σ = 0.05, 0.10 and 0.15 day), with a 50:50 geometry prior. Then judge the per-world payoffs by the method
     of `experiments/gate/options_tail.py` (P(Overall ≥ c) at our 75th-99th percentiles), the form in which the
     captain weighed R2 and D3.
  6. Decision rule, fixed in the row: put sheet B to the captain as "recommended" if the backup pair beats today's
     pair by at least +0.3 board points at sd = 0.25 and loses no more than 0.3 at sd = 0.12; otherwise as "not
     recommended".

### C2. One organiser question that would legitimately open both embryo hedges. Recommended (cheap)

- **Facet and targets.**
  - Facet: spatial organisation (geometry), and through it the choice of s2. It also sets C1's staging prior.
  - Targets shape and growth (s1 12.63 and 8.86 left; s2 11.16 and 12.99) and cell state.
- **Mechanism.** The test geometry and depth are set by acquisition, not by development, and the organisers know
  E7.75's protocol and staging convention.
  - E7.5 was cut in 34 sections at step 14, into a near-isotropic cloud: z range 462 against xy ranges 464 and 408
    (P3 file summary in `experiments/coords/results/p3/`).
  - z is the thin axis at E8.0 and the long axis at E6.75 and E7.25.
  - On the backtests the truth was once rounder than both brackets (E7.5: 0.53 against 0.38 and 0.36) and once
    thinner than both (E7.25), so no interpolation predicts thinness.
  - Two facts would flip decisions: whether E7.75 was cut and imaged like E7.5 or like E8.0, and whether the embryo
    labels are gestational ages or morphological stages.
- **SOTA basis.**
  - MERFISH detection and geometry are set by run and section: 2-fold detection differences between sections,
    attenuation with depth in the section, and readout-group losses in 31% of sections [task-c 20][task-c 19]
    [task-c 21].
  - The data lab assembles its stacks with non-rigid alignment and states no z rule [task-d 27].
  - Frames are per-embryo and unregistered [task-e 4].
  - Task d names a public acquisition protocol as what would change its conclusion (task-d, Gaps).
- **What is new.** It narrows follow-up Q5 (sent 2026-10-19, unanswered) and Q7, and ties them to the two embryo
  decisions they would settle.
- **Legality.**
  - Rules s.10 protects an answer given in advance.
  - Asking for a public answer means no team holds private information.
  - It asks for protocol metadata, not for a measured property, and quotes no limit, anchor or count.
  - Draft for the captain (a sixth email, sent privately, public answer requested): "Two short questions on the Task
    2 embryo series, whose answers apply to every team: (1) Was the E7.75 test stack acquired with the same protocol
    as the E7.5 validation stack (section orientation relative to the embryo, section pitch and number, MERFISH
    chemistry and imaging run), or as the E8.0 training stack? (2) Are the embryo-series stage labels (E6.75 to
    E8.0) gestational ages from the plug, or morphological stages of each embryo?"
- **Decision branches.** Each switch remains the captain's ruling.
  - Q1 "like E7.5": keep E7.5's geometry in s2. With depth near E7.5's, C1 is strongly favoured, and C5 turns
    positive.
  - Q1 "like E8.0": the aligned geometry (R2, about +0.15 Overall) becomes the better s2 geometry, and C1's depth
    tail vanishes.
  - Q2 "morphological": σ is low; keep today's s2 unless C1's table is strongly positive.
  - Q2 "gestational": σ is about 0.1-0.2, and C1 is favoured.
- **Counted evidence.** None. An organiser fact enters only through the captain's ruling.
- **Expected gain.**
  - P(an answer before 2026-11-26) ≈ 0.1: none of five emails has had a reply, and no public change is dated after
    2026-10-28 [task-e 11].
  - If answered: +1 to +3 board points from aiming the pair, i.e. +0.5 to +1.5 Overall.
  - Expected value: **+0.05 to +0.15 Overall**. It adds no variance.
- **Effort and risk.** An hour to draft; the captain sends it. The only risk is no answer.
- **First step.** Add the draft to `docs/ORGANIZER_QUESTIONS.md` as the sixth email.

### C3. Bold: an emergence-first backup (states that emerge, at E7.5 depth)

- **Facet and targets.**
  - Facet: which states emerge, and their proportions.
  - Targets expression change (8.96 and 7.69 left), the group C1 cannot reach, as well as cell state and
    neighbourhoods.
- **Mechanism.**
  - s2 keeps E7.5's real cells for every state E7.5 has.
  - It adds E8.0's real cells only for the states E7.5 lacks, at half their E8.0 share. They are placed by fgw's
    coupling in E7.5's geometry, at the size rule and count rule.
  - Its dp is the emergence direction alone: new-state markers rise and the rest is diluted, which de_direction's
    control for reference level absorbs. It does not carry the half-mixture's full E8.0 − E7.5 direction, so it is
    the only backup that decorrelates the expression group from s1.
  - Everywhere else it stays near copy_last, with depth mostly E7.5's.
- **SOTA basis.**
  - Carrying a type's expression forward is unbeaten for type-specific temporal changes; only changes shared across
    types can be forecast [task-a 5].
  - Copy-nearest is strong [task-d 8].
  - A graft keeps whatever property it copies [task-e 24].
  - The captain's principle asks which states emerge (`AGENTS.md`).
- **What is new.**
  - Nearest pruned items: "T1 births by the legal route", and "Embryo: exponential-growth state shares" (share change
    0.014-0.041 against the interpolation, about +0.05, and a window breach as proposed).
  - This version uses the test inputs only, serves as a backup, and moves shares far more.
- **Legality.** As for C1: inputs E7.5 and E8.0 only, labels harmonised within Task 2 by `vec.harmonise`, no
  external atlas. It is a new s2 method, so it needs the captain's ruling. Count, size and scope are unchanged.
- **Counted evidence.** None: E7.25 and E7.5 cannot count, and no other embryo pseudo-task exists. It is a captain's
  bet.
- **Expected gain.**
  - It competes with C1 for the one slot.
  - Board mean about +0.1 (range −1.4 to +2.0): **+0.05 Overall (range −0.7 to +1.0)**.
  - It carries more variance than C1. It wins where the near-E7.5 truth's DE is led by emergence and loses where it
    is led by maturation.
  - Upside: +2 to +4 Overall in its world.
- **Effort and risk.**
  - Effort: 1-2 days for a state-resolved draw in `interpolate_clouds` and its tests. The default path is unchanged,
    so the golden table stays.
  - Risk: no check of any kind, and a new code path three weeks before the freeze.
- **First step.** A target-free sizing row on the test pair.
  - It measures the share of E8.0 cells in states with under 1% share at E7.5 after harmonisation, and the Spearman
    correlation of this backup's dp with the half-mixture's dp over the panel genes.
  - Kill rule, fixed in the row: drop it if the share is under 5% (it is then a copy) or the Spearman exceeds 0.9 (it
    is then C1 with extra steps).
  - Otherwise, add it as a fourth column of C1's table.

### C4. A check of C1's asymmetry on real intermediates (evidence that does not count)

- **Facet.** None scored; this is evidence for C1's ruling.
- **Mechanism.** C1's table scores mixtures against mixtures. A series whose middle stage is a real intermediate
  tests whether a lean-early backup still pays when t' is small. Two sources:
  - Flysta3D-v2 (window-free): 43 Drosophila embryos at 0.5-2 h intervals, as 3D single-cell Stereo-seq from serial
    sections [task-d 26]. It gives many triples at exactly t = 1/2, bilaterally symmetric like the mouse at E7.75.
  - The heart triple E8.25_late → E8.5 → E8.75 (same platform): t = 1/2, with the earlier bracket as floor, which is
    the embryo test's structure. But it is a single triple, its scope changes (E8.25_late carries yolk sac and
    non-cardiac cells), and it reads a released truth, so its runs count as tries (rule 4).
- **What is new.**
  - Task d's lever 1 used Flysta3D for geometry: re-deciding E11, which was pruned as "s1 = fgw_aligned re-read on a
    test-like task", and R2, which was declined.
  - Here it answers only the mixing-fraction asymmetry, which no pruned item covers.
  - The heart triple repeats the pruned "third pseudo-task (heart E8.5)", in a role that does not count.
- **Legality.**
  - Flysta3D is window-free (not mouse) and raises no genotype question. Because it informs a method choice, it
    needs an `EXTERNAL_DATA.md` row first, and its CNGB terms must be read first (the page states no licence)
    [task-d 26].
  - Heart E8.5 is an unlocked validation file (A20), in the same task (T2) and the same modality.
- **Counted evidence.** None, by design. An external family is not an adoption reading, and the captain declined
  such a family for T1 this season.
- **Expected gain.** Decision value only, about +0.0 to +0.1 Overall.
- **Effort and risk.** 3-5 days and some DGX CPU hours for Flysta3D; half a day for the heart triple. Risks: transfer
  across species and platforms, and time.
- **First step.** Read CNGB's terms and draft the register row. Run the check only if C1's table lands close to its
  threshold.

### C5. Bold: retire the flat-geometry bet (conditional on C2)

- **Facet and targets.** Spatial organisation; targets shape and growth (s1 12.63 left on E7.5).
- **Mechanism.**
  - With C1 in the s2 slot, s1 would carry both the mid-t' bet and the geometry bet.
  - This candidate makes s1's geometry central instead of flat, by one of two routes:
    - the mirror-only correction of FGW's partners (keep FGW's moves and put each partner on its anchor's side; not
      tried);
    - the aligned moves of E11 (thinnest axis 0.43 at test).
  - At test, FGW's partners are mirrored in 10 of 10 draws (N1).
- **SOTA basis.**
  - No 2025-2026 GW method addresses the mirror ambiguity.
  - TP-DATE beats a learned static-FGW flow by 7-8% in spatial W2, but is never compared with real cells moved
    along a closed-form coupling [task-d 4].
  - Geometry is set by acquisition [task-c 20][task-d 27].
- **What is new.**
  - Nearest items: "A mirror-only correction of FGW's partners" (the coupling line is closed for s1 unless the
    captain rules otherwise), E11 (vetoed on E7.25 at −1.59) and R2 (declined).
  - It is new only as a package with C1, conditional on C2's answer.
- **Legality.**
  - Reopening the coupling line and overriding E11's veto are the captain's calls alone.
  - A changed s1 geometry sends the embryo s2 back to him (`docs/S2_RULES.md`).
  - Count, size and scope are unchanged; no external data.
- **Counted evidence.**
  - E11's own result: E7.5 +2.13 (ahead in 10 of 10 seeds), E7.25 −1.59 (a veto).
  - The mirror-only form would face the same veto in E7.25's five mirrored seeds.
  - There is no counted path.
- **Expected gain.**
  - As it stands, about −0.1 Overall (range −0.6 to +0.5).
  - It gains about +1.5 to +2 board points in mid- or far-t' worlds where the truth's geometry is round or
    intermediate.
  - It loses about 5 in a flat-truth world like E7.25 (fgw_shape0 scored 5.11 below fgw there). That world keeps a
    probability of about 0.15-0.25.
  - If C2 answers "cut like E7.5", the flat world becomes unlikely and C5 turns positive: +0.2 to +0.5 Overall.
  - It adds no upside beyond that branch.
- **Effort and risk.** Half a day (the aligned moves exist) to a day (code for the mirror-only form). Risk: a veto
  override, which the ledger's rules discourage.
- **First step.** None unless C2 is answered.

## 5. SOTA reviewed and set aside (the stop on s1 confirmed)

| Item (source) | Use on this board | Nearest pruned or closed item | Why it does not change the picture |
|---|---|---|---|
| Learned spatiotemporal flows: ContextFlow, stVCR, TP-DATE, NicheFlow [task-d 1][task-d 3][task-d 4][task-d 2] | generate states and positions at t = 1/2 | "displacement interpolation in a count-model latent space" (7-10 days, about +0.15); E8 | 5-20% margins over other learned flows; no copy-nearest baseline; NicheFlow tested in-distribution only. E8 showed that averaged cells form too narrow a population (mmd_u 0.671 → 0.638, ahead in 0 of 10). Expected ≤ 0 |
| Mimyr [task-d 8] | generation against lookup | - | Copy-nearest is strong: an atlas of about 10M cells buys +0.06 to +0.10. Supports real-cell copies |
| SOCS structure-GW term [task-d 5] | a coupling that keeps labelled structures together | fgw_tissue; partial or unbalanced FGW | Mirror-blind; the coupling line is closed. About 0 |
| SINTER3D, isoST, UniST [task-d 6][task-d 7][task-d 10] | continuous fields; isotropic resampling | - | The scorer reads native section geometry, and a 1.5x z mismatch costs 12.6-12.8 heart points. A guard only: keep native z |
| Gastrulation spatial atlases: SEU-3D, Harland, sc3D, ST-MAGIC [task-d 21][task-d 22][task-d 23][task-d 31]; Mittnenzweig [task-d 24] | bracket references, or a clock | clock warp; exponential-growth state shares | Legal only at E7.25 or earlier, or at E8.0 (sc3D's E8.5 is inside the heart window for a shared step). None places E7.75. In-window results are barred |
| Developmental clocks across datasets: Gene-Chronos [task-a 12]; a pooled maturation clock [task-a 6] | place the brackets on a time axis | clock warp | Errors of 1.52 days across human embryo datasets, and 8.03 weeks for one study's clock applied to other studies, against a 0.5-day bracket. Clocks trained on scRNA also meet A15 |
| Pretrained temporal and world models: Navigo [task-e 16][task-e 17], Chreode [task-a 13][task-e 23], CellTempo [task-a 3] | expression dynamics | foundation-model screen of 2026-10-19 | Contaminated corpora or unstated provenance. Navigo's atlas starts at E8.5 (no gastrulation) and needs unspliced counts |
| Count generators: count-FM, scLDM [task-a 18][task-a 19] | realistic counts | negative-binomial resampling; decoded sampling | Unconditional gains only |
| MERFISH QC and run structure: MerQuaCo, readout position, MERFISH 2.0 [task-c 20][task-c 19][task-c 21] | depth or section corrections | depth line (E7: −4.18 on E7.5); depth-shaped draw | Depth is closed, and a section filter changes scope (closed). Used here only as prior information (C1, C2) |
| The organisers' spatiotemporal_ode growth rate [task-e 4] | size | - | This is already our PCHIP size, and size is closed |
| Community same-origin pairing; spatial transplants [task-e 24][task-e 25] | coupling variants | coupling sweeps | The coupling line is closed |
| A reseeded s1 in the s2 slot [task-f lever 2] | harvest draw noise | - | Two iid seeds give E[max] ≈ mean + 0.56 seed SD (0.83 on E7.5, 2.20 on E7.25). That is less than today's geometry hedge: +0.92 on the backtests' better-of-two, and more at test, where s1 is flattened in every draw. It also sits near the integrity line on draw noise. Rejected |
| Embryo tail decoder (N5) | free the infeasible bottom-tail slots | pruned 2026-11-05 | 1 of the bottom 40 genes is infeasible; nothing new |

## 6. Path to 209 (embryo)

> Struck 2026-11-06 (main session, after the counter-review): this passage measures candidates against the captain's 209, which gates uploads only and is never a threshold or evidence for a try, an adoption or a method choice (docs/DATA_POLICY.md 2026-11-05; docs/PLAN.md section 8). Information only; nothing was chosen by it.

- **Where it stands.** The embryo board counts 0.5. On committed evidence it scores 60.87 (better-of-two;
  10th-90th percentile 56.9-63.5), which is 30.4 Overall.
- **The gap.** Reaching 209 takes +17.5 Overall over reading A (191.5), or +25 over the outlook (184).
- **The embryo's assumed share.**
  - About +0.15 Overall in expectation, from C1 after its table, with C2 sent.
  - +2 to +4 Overall only in the world where E7.75 sits near E7.5 and the captain has adopted C1. That world has
    probability 0.16-0.37, times the chance of his yes: roughly 0.1-0.3 in all.
  - So the embryo can supply about a tenth to a fifth of the gap in its best case (probability about 0.1-0.3), and
    under 1% of it in expectation.
- **No bigger lever exists.** Nothing on record moves the embryo's mean by more than about a board point. Even A7's
  old target of 70 would add only about +4.6 Overall.
- **Conclusion.** A total above 209 has to come from T3 and T1 landing near their best cases. The embryo's job is to:
  - keep its 61;
  - not lose the near-E7.5 world;
  - have its s2 ruled before the first embryo upload on 2026-11-27.

## 7. Process notes

- **Read.**
  - The brief and `final/task-a.md` to `task-f.md`.
  - `AGENTS.md`, `docs/CHALLENGE.md`, `METHOD_SUMMARY.md` and `experiments/harness/TRIES.md` (all).
  - `docs/PLAN.md`: sections 1-3, 4.2, 7 and 8, and the decision log from 2026-11-04 to 2026-11-06.
  - `docs/DATA_POLICY.md`: windows, Scores, the do-not-use list, and change notes 2026-11-04 to 2026-11-05 (k).
  - `docs/S2_RULES.md`, `ASSUMPTIONS.md`, `docs/ORGANIZER_QUESTIONS.md`.
  - `experiments/research_2026_11/REPORT.md`: the summary and sections 1-2 and 6-10.
  - `results/upload_gate.json` and `experiments/gate/results/estimate_2026-11-05.json`.
  - For figures:
    - the per-seed skills in the committed practice-task files of `emb_shape_v1`;
    - the embryo lines of `experiments/coords/REPORT.md` and the E7.5 row of its P3 file summary;
    - the signature of `vec.t2_build.interpolate_clouds` and the embryo options of
      `experiments/builds/build_test.py`.
- **Not opened.**
  - boards/, data/, runs/, docs/research/, and the workflow's raw/ and screened/ folders.
  - Anything under ~/.claude-account or ~/.claude.
  - The leaderboard and any Yukon material. No web.
  - The untracked `experiments/gate/path_to_target.py` and its result file (another agent's uncommitted work).
- **Seen incidentally.** Inside a committed practice-task file: the E7.5 validation board's spec and anchors, and that
  practice task's truth DE counts. Allowed since the unlock; none is used or quoted.
- **Computation.** Closed-form arithmetic only: finding 2's normal probabilities, the log-linear depths and sheet B's
  break-even. No data run, no toy and no try.
- **Not done here.** The DGX and local sync, and the check for unpushed commits, that the user asked for before this
  workflow belong to the main session. Seen read-only: the local repository is on main at 9ce2742 with no commits
  ahead of its upstream, and two files are untracked under experiments/gate/.
