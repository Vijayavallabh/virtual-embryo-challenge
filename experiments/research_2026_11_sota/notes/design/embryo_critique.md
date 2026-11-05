# Critique of the T2:embryo:test_interp design (E7.75 MERFISH, counts 0.5 in Overall)

Adversarial critic for team Mithril, SOTA workflow of 2026-11-06. AS_OF 2026-11-06. Object: `design/embryo.md`.
Read-only. Nothing in the repository was run or changed. No target, test-board limit, anchor, leaderboard figure or
barred item was read. One closed-form computation, in this agent's scratch with output to stdout only, extends the
skill-map toy with a non-mixture residual and the two-upload geometry pairing (section 3). It read no data.

## 1. Bottom line

| Candidate | Verdict | Designer's gain (Overall) | Corrected gain (Overall) | Main reason |
|---|---|---|---|---|
| C1: lean-early s2 (t = 1/4, E7.5 geometry; REPORT sheet B) | RULING-ONLY, with its pricing revised | +0.15 (-0.5 to +1.0) | about -0.05 (-1.0 to +0.8) | Legal and doable before the freeze. It is priced on sheet B's pure-mixture toy, though. Once the skill map carries the residual that committed skills show, break-even rises to P(t' < 0.375) of about 0.3-0.6, at or above the designer's own legal prior (0.16-0.40) |
| C2: organiser question | KEEP-REVISED (cut down) | +0.05 to +0.15 | about +0.01 (0 to +0.1) | Re-asks Q5 and Q7, both unanswered since 2026-10-19. Its E7.75-stack half asks for class-A territory on closed levers (z extent, so size; imaging run, so depth) |
| C3: emergence-first s2 | KILL | +0.05 (-0.7 to +1.0) | about -0.15 (-1.0 to +0.6) | Competes for C1's slot with a near-copy bet. No check is possible, its expression-decorrelation premise rests on finding 1 (which is overstated), and it needs a new code path |
| C4: real-intermediate check (Flysta3D; heart triple) | KILL | 0 to +0.1 | about 0 (at most +0.02) | The heart triple is the pruned third pseudo-task and counts as tries. Flysta3D transfers neither its residual nor its staging spread to mouse MERFISH, and no licence is stated. Committed mouse skills calibrate C1 better, at no cost |
| C5: retire the flat geometry | KILL | -0.1 (-0.6 to +0.5) | about -0.1 if adopted; about 0 as the conditional it is | It overrides a counted veto (E11 on E7.25) and a declined ruling (R2). Its trigger, C2's answer, has P about 0.03-0.05 by 2026-11-26 |

What this design is worth: about 0 Overall in expectation (C1 about -0.05, C2 about +0.01), not +0.15. The stop on
s1 is right. Nothing here changes the embryo's s1, its 60.87 or the path to 209.

## 2. The designer's three findings, checked

**Finding 1 (rank invariance) is partly wrong.**
- The scorer computes dp = pb(pred) - pb(ref) and dt = pb(true) - pb(ref) against the same 10% reference working copy.
  So dp and dt share a noise term, -pb(ref) noise.
- The committed results show this term is not small:
  - a fresh copy scores de_direction skill 0.526 and 0.628 on the two embryo backtests, above the floor's 0.5
    (`METHOD_SUMMARY.md`);
  - T5 found that the reference copy carries most of a prediction's dp noise.
- Lowering t gives the shared term more weight in dp. That helps where the truth's own change is small (the near-E7.5
  world) and costs where it is large. So the expression group is partly hedgeable, in the same direction as cell
  state.
- "A small cost of a lower fraction in every world" is wrong in the near world.
- Size: a few hundredths of de_direction skill. Treated as half-strength extra mixture weight, it moves C1 by +0.0 to
  +0.15 Overall. It changes no verdict, but it undercuts C3's claim to be "the only backup that decorrelates the
  expression group".

**Finding 2 (staging prior): the arithmetic is right and the input is judgement.**
- The numbers check out: P = 0.16 / 0.31 / 0.37 / 0.40 at sigma = 0.05 / 0.10 / 0.15 / 0.20 day, and the SD of t' is
  about 2.45 sigma.
- sigma itself is sourced nowhere:
  - [task-d 17] gives six forelimb sub-stages at E14.5 and no spread in days. Task d lists the per-embryo sub-stage
    spread as "read only partly".
  - [task-b 7] is a mis-citation. Source [7] of task b is a postnatal heart paper. The claim about somite counts in 2-h
    bins is task-b finding 7, which cites source [5] (Qiu 2024).
  - [task-c 16][task-c 17] are method notes.
- So "P = 0.16-0.40" restates sigma.
- The designer's own reading of the quarter-day labels and of `E8.25_late` (morphological staging) points to the low
  end: sigma about 0.05 and P about 0.16. There C1 loses in every model below.
- t' in the skill map is the metric position of a real, non-mixture embryo, not an age. Its spread also carries
  developmental nonlinearity and embryo or run idiosyncrasy, which legal data cannot price in either direction: the
  backtest positions are barred both ways (`docs/DATA_POLICY.md` 2026-11-05 (b)(v); REPORT section 6).
- Treat sigma as the captain's stated prior, not as a derived quantity.

**Finding 3 (depth) uses its evidence selectively.**
- Three released interior points exist. The designer keeps the two that sit within 5% of log-linear interpolation
  (E7.5 at 62 against 64; heart E8.5 at 166 against 174). It sets aside the third, E7.25 (35 against 49, 0.71x),
  because it is the barred mixing-fraction feature.
- But the same bar covers backtest positions both ways, and depth is set by acquisition: section and run detection
  differ about 2-fold [task-c 20], and the embryo series' per-day rate was 0.31 over E6.75 -> E7.25 and about 2.4
  after.
- One interval in three far off trend puts the "E7.5-like depth" tail at about 0.1-0.3, not 0.1.
- This cuts in C1's favour, by about +0.05 to +0.3 Overall from the depth hedge sheet B gave, and is included in
  C1's corrected range.

## 3. C1 re-priced: the residual and the two-upload pairing

**What is priced now.**
- Sheet B and the designer use skill(t; t') = t'^2 / (t'^2 + (t - t')^2): a world of pure bracket mixtures, where a
  matched mixture scores skill 1.0.
- The designer's G of 4-8.5 board points is backed out of sheet B's break-even of 0.15-0.2 and loss of 1.5.

**Two things the toy leaves out.**
1. *Residual (the transfer discount).* Real intermediates are not mixtures.
   - Committed cell-state skills at a nominal mixture are 0.671 (mmd_u) and 0.664 (variogram) on E7.5.
   - The better committed mixture in TRIES' `diag_emb_mix` table still sits at 0.0380 against a copy's 0.0858, a skill
     near 0.7.
   - Model: skill = (rho + t'^2) / (2 rho + t'^2 + (t - t')^2). A matched mixture then scores 0.8 at rho = 0.083 and
     0.7 at rho = 0.19.
   - This calibrates to committed skill levels, not to any backtest position.
   - The residual damps exactly the near-floor blow-up that creates the asymmetry. Per-metric skill of t = 1/4 minus
     t = 1/2:

   | t' | 0.1 | 0.2 | 0.3 | 0.5 | 0.7 |
   |---|---|---|---|---|---|
   | toy (rho = 0) | +0.25 | +0.63 | +0.28 | -0.20 | -0.22 |
   | rho = 0.083 | +0.19 | +0.17 | +0.09 | -0.11 | -0.16 |
   | rho = 0.19 | +0.12 | +0.09 | +0.04 | -0.06 | -0.11 |

2. *Pairing.*
   - With C1, t = 1/2 survives only in s1, whose moved cloud is flatter than either bracket in every test draw (N1).
   - So in the central world (t' near 1/2, a truth that is not flat), neither upload is right on both counts.
   - Pair gain = max(dMix + dS, 0) - max(dS, 0), where dS is s2's geometry edge over s1. The base case has:
     - +4 board points in a round world and +2 in an intermediate one. That is above the backtests' +1.83, because s1
       is flatter at test.
     - -5 in a flat world like E7.25, with prior 0.2.
   - dMix = 50 points x the per-metric difference (cell state and neighbourhoods).

**Result** (Overall; t' ~ N(0.5, sd), sd = 2.45 sigma; base geometry):

| Skill map | sigma 0.05 | sigma 0.10 | sigma 0.15 | Near-world gain G / far loss L (board pts) | Break-even P |
|---|---|---|---|---|---|
| toy, rho = 0 (what a same-embryo half-split table measures) | -0.23 | +1.39 | +1.50 | 10-14 / 2.3 | 0.14-0.19 (= sheet B) |
| rho = 0.083 | -0.72 | +0.03 | +0.33 | 2.8-5.7 / 2.2 | 0.28-0.44 |
| rho = 0.19 | -0.72 | -0.26 | -0.02 | 1.4-3.6 / 2.0 | 0.37-0.58 |

- Over geometry gaps from (+2, +1) to (+6, +3), the calibrated rows span -1.2 to +0.7 Overall.
- Adding the depth tail (+0.05 to +0.3) and finding 1's noise channel (+0.0 to +0.15) gives **about -0.05 Overall
  (range -1.0 to +0.8)**.
- In the near-E7.5 world it adds +0.7 to +2.9 Overall, not +2 to +4.
- The toy row reproduces sheet B's break-even. That supports this extension as a correction of sheet B's own model,
  not a different one.

**What follows for the designer's table (C1 first step).**
- Its pseudo-truths are pure mixtures of the same embryo's alternating-section halves, so the residual is about 0.
- Its "E7.5 geometry" world is E7.5's own shape, which today's s2 and C1's s2 share almost exactly.
- So it measures the toy row: positive at sd 0.25, with the sd 0.12 clause near or past its -0.3 line depending on
  the geometry gap.
- It cannot measure the two things the decision turns on: the residual and sigma.
- Read honestly, the designer's own rule is about as likely to say "not recommended" as "recommended". The headline
  +0.15 assumes the table says yes.

## 4. Candidates

### C1. Lean-early backup (embryo s2 at t = 1/4 in E7.5's geometry): RULING-ONLY

1. **Repeat?** Not a pruned item.
   - It is REPORT sheet B, already on the captain's list (decision log 2026-11-05 (f)).
   - It sits next to R2, declined on 2026-11-05: an after-the-fact amendment of a frozen rule, with a mean of +0.08
     Overall. By that precedent, an amendment needs a clearly positive price, which C1 does not have once calibrated.
   - What is new is arithmetic and judgement, not data.
2. **Rules.** Legal.
   - Inputs are E7.5 (A20) and E8.0 only.
   - Embryo composition is open; count, size rule and scope are unchanged.
   - No A15, window or genotype issue. No score is used.
   - Process gaps the design does not state:
     - (a) s2 is built "at the commit s1 is built from" (`docs/S2_RULES.md`). So the builder flag must be committed and
       tested by 2026-11-26, whatever the ruling date. `interpolate_clouds` already takes t, and `build_test.py` needs
       one option.
     - (b) The table's decision rule goes to the captain as a plain summary for sign-off before the run. Decision log
       2026-11-04 (d) asks this for counted runs; it applies here too, because this rule feeds a ruling on a frozen
       rule.
     - (c) The N row should say that E7.5 enters as the test input, not as the E7.5 backtest's truth.
     - (d) A pass/fail-only check, printing no count, that 3/4 of the count rule's N fits in E7.5's deduplicated pool.
       Otherwise `interpolate_clouds` repeats rows, and s2' differs from today's s2 in more than t.
3. **Evidence.** No counted evidence is possible. The table prices C1, and as built it is biased toward C1 (section 3).
4. **Gain.** About -0.05 Overall (-1.0 to +0.8). The tail value at our 90th-95th percentile is about +0.03-0.09
   (variance only), not +0.1-0.3. The rank-view argument ("teams interpolating at 1/2 lose together") needs lambda,
   which REPORT section 1 says cannot be read.
5. **Feasibility.** About 1 day of code and hours of DGX CPU, well before the freeze.

**Revision.**
- Keep t fixed at 1/4: no sweep, since each value would be a choice made on synthetic truths.
- Price it on the calibrated closed form, with rho from committed matched-mixture skill levels (0.7-0.8) and
  explicit round, intermediate and flat worlds.
- If the table runs at all, use it only to check the noise terms, against the rho = 0 row.
- Give the captain one line per sigma, with the break-even.
- Recommend adoption only if he states sigma >= 0.15 day (gestational labels) and accepts the low-residual end.
  Otherwise the default stays: today's s2.
- Rule before the freeze and before the first embryo upload, in its own commit with his words and the new sha256.

### C2. Organiser question: KEEP-REVISED

1. **Repeat?** Yes, of unanswered questions:
   - follow-up Q5 of 2026-10-19 already asks whether the test stacks follow the training files' protocol (section
     thickness, number of sections);
   - Q7 already asks for the staging convention.
   - There is no new reason to expect an answer: five emails, no reply, and no public change after 2026-10-28
     [task-e 11].
2. **Rules.**
   - "Section pitch and number" of the E7.75 stack is its z extent, so a size property. Sizes are closed by the
     post-exposure block, and "stack-limited embryo geometry (z extent)" is a closed scenario (TRIES, "After the first
     wave"). Task d itself calls the test stack "class A territory" (mapping (d)).
   - "Imaging run and chemistry" is asked in order to infer the test target's depth, and the depth line is closed.
   - Rules s.10's "answer given in advance" protects the use of borderline data. It does not cover asking for
     hidden-target metadata.
3. **Evidence.** None. An organiser answer enters only through the captain.
4. **Gain.** About +0.01 Overall.
   - P(an answer by 2026-11-26) is about 0.03-0.05.
   - The answer is worth little: C1 sits near break-even, and C5 cannot move a counted veto.
   - A public answer also helps the field.
5. **Feasibility.** An hour.

**Revision.** One general line, added to the next email the captain sends or as a reply on the fifth-email thread,
not a standalone sixth email: "Are the embryo-series stage labels (E6.75 to E8.0) gestational ages from the plug, or
morphological stages of each embryo?" Drop every E7.75-specific acquisition item. For those, Q5 already stands.

### C3. Emergence-first backup: KILL

1. **Repeat.** It is the pruned "exponential-growth state shares" bet (share changes against the interpolation, about
   +0.05) made bigger ("moves shares far more") with no new evidence. It competes for the same s2 slot as C1.
2. **Rules.** Legal on data, as a new s2 method that needs a ruling. "States E7.5 lacks" depends on forced
   cross-stage label calls. The births design found such calls unreliable: argmax calls forced, and its confident-call
   test caught only 28-63% of present cells (TRIES, T1 births).
3. **Evidence.** None possible. The kill rule only sizes it. Nothing target-free can tell an emergence-led near-E7.5
   truth from a maturation-led one.
4. **Gain.** About -0.15 Overall (-1.0 to +0.6).
   - It is effectively a near-copy backup: a skill near t = 0.1 on the mixture axis, with no maturation.
   - In mid and far worlds it loses more than C1, and the pair's only t = 1/2 entry there is the flattened s1.
   - Its SOTA basis is misread:
     - [task-a 5] is carry-forward for unseen types in a human immune extrapolation, not interpolation with a known
       later bracket;
     - [task-e 24] is a community repository of target-aware toys.
5. **Feasibility.** A new state-resolved code path 1-2 days long, three weeks before the freeze. That is a build risk
   for no check.

### C4. Real-intermediate check: KILL

1. **Repeat.**
   - The heart triple is the pruned "third pseudo-task (heart E8.5)". It reads a released truth, so its runs count as
     tries (rule 4), on a board the ledger told to stop polishing.
   - Task d's audit flags Flysta3D as a repeat of "s1 = fgw_aligned re-read on a test-like task" and of R2. The
     mixing-fraction role is new.
2. **Rules.**
   - Flysta3D states no licence (task-d audit). Rules s.10 needs licences that permit use, and it needs an
     `EXTERNAL_DATA.md` row first.
   - The heart files are legal.
3. **Evidence.** Non-counting by design. The two quantities C1 needs, the residual and the spread of the effective
   position, belong to the platform, species and collection:
   - Stereo-seq cell bins, 0.5-2 h Drosophila collections;
   - for the heart triple, a falling depth and E8.25_late's wider scope, n = 1.
   - None of them transfers. Committed mouse-MERFISH skills calibrate the residual at no cost (section 3).
4. **Gain.** About 0, at most +0.02 Overall of decision value.
5. **Feasibility.** 3-5 days plus a licence read, for a check the designer would run only "if C1's table lands close to
   its threshold".

### C5. Retire the flat-geometry bet: KILL

1. **Repeat.** E11 (vetoed on E7.25, -1.59), the pruned mirror-only correction (it meets the same veto) and R2
   (declined). The coupling line is closed for s1.
2. **Rules.** An organiser fact cannot turn a counted veto into a pass. Only a captain's override of the adoption rule
   could, which the ledger's rules discourage. A changed s1 also sends s2 back to the captain.
3. **Evidence.** The counted evidence is against it. E7.25 shows that a truth flatter than both brackets exists.
4. **Gain.** About -0.1 if adopted. As proposed, conditional on C2 (P about 0.03-0.05), it is about +0.01.
   - Packaged with C1, it removes the board's only cover in the flat world, since neither upload would then keep a
     flat cloud.
5. **Feasibility.** Only if C2 is answered before the freeze, which is unlikely.

## 5. The rest of the design

- **What holds.**
  - The stop on s1.
  - The set-aside table (section 5), including the rejected reseeded s1.
  - Every committed figure checked: the headroom rows, 60.87 (56.9-63.5), 59.95, the N1 geometry, E7.5's 34 sections
    at step 14, the depth medians and the staging arithmetic.
  - No leaderboard figure sets a threshold.
  - 209 frames section 6 and chooses nothing.
- **Corrections.**
  - The embryo's expected share of the gap to 209 is about 0, not +0.15.
  - The near-world upside is +0.7 to +2.9 Overall (calibrated), times P(ruled yes).
- **Citation and reading slips.**
  - [task-b 7] should be task-b finding 7, which cites source [5].
  - The contest records [task-f 13][task-f 14][task-f 24][task-f 25] support "score-chasing fails". Task f's own
    reading is that they show "no average payoff for boldness" (mapping 7), and [task-f 25] was read only as a search
    summary.
  - [task-d 8] (a copy-nearest baseline across a spatial gap) says nothing about where E7.75 sits.

## 6. Strongest candidates and the biggest risk

**Strongest.**
1. **C1, as a captain's ruling only, re-priced** on the residual-calibrated skill map with the geometry pairing:
   - the code goes into the freeze commit;
   - the ruling comes before the first embryo upload;
   - the default stays today's s2 unless the captain holds sigma >= 0.15 day.
2. **C2, cut to the one staging-convention line**, added to the next email.

**Biggest risk in the designer's reasoning.**
- C1's value is priced on sheet B's pure-mixture toy, which assumes a matched mixture scores skill 1.
- The skill-map asymmetry behind the lean-early backup is large only when the truth has no non-mixture residual.
- Committed skills (about 0.67-0.7 at the embryo backtests) show a substantial residual. With it:
  - the near-world gain shrinks 2-4x;
  - break-even rises from 0.15-0.2 to about 0.3-0.6, at or above the designer's own legal prior of 0.16-0.40;
  - the expected gain falls to about zero.
- The proposed same-embryo half-split table has no residual. It would reproduce the toy and could green-light an
  after-the-fact amendment of a frozen rule whose calibrated mean is about zero.

## 7. Process notes

- **Read.**
  - `design/embryo.md`, `brief.md`, `final/task-d.md` in full, and the cited items and findings of task-a (task-a.md),
    b, c, e and f. task-a2 was checked for numbering only.
  - In the repository:
    - `docs/CHALLENGE.md`, `METHOD_SUMMARY.md`, `experiments/harness/TRIES.md` (all), `docs/PLAN.md` 4.2, 7, 8 and
      the decision log 2026-11-04 to 2026-11-06;
    - `docs/DATA_POLICY.md` (data sections, Scores, do-not-use list, change notes 2026-11-04 to 2026-11-05 (b)),
      `docs/S2_RULES.md`, `ASSUMPTIONS.md`, `docs/ORGANIZER_QUESTIONS.md`;
    - `experiments/research_2026_11/REPORT.md` (all), `results/upload_gate.json`,
      `experiments/gate/results/estimate_2026-11-05.json` (embryo and total lines);
    - the header and E7.5 row of `experiments/coords/results/p3/file_summary.csv` (a validation stage);
    - the signature of `vec.t2_build.interpolate_clouds` and the embryo options of `experiments/builds/build_test.py`.
- **Not opened.**
  - boards/, data/, runs/, docs/research/, the workflow's raw/ and screened/ folders.
  - Anything under ~/.claude-account or ~/.claude. One oversized listing was saved there by the harness; it was not
    opened, and the listing was redone in smaller form.
  - The untracked `experiments/gate/path_to_target.py` and its result file.
  - No web, no leaderboard.
- **Not quoted.** No best mixing fraction derived for a counted embryo target (DATA_POLICY 2026-11-05 (b)(v)). The
  residual is calibrated to committed skill levels only.
- **Not done by this step.** The user asked for the DGX and this machine to be synced, and for unpushed commits to be
  checked, before the workflow. That belongs to the main session. Seen read-only here:
  - local `main` at 9ce2742, with nothing ahead of `origin/main` as last fetched;
  - two untracked files under `experiments/gate/`;
  - the DGX was not checked.
