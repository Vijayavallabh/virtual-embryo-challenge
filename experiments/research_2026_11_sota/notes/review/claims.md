# Claims and traceability review: `experiments/research_2026_11_sota/REPORT.md` (draft of 2026-11-06)

Reviewer role: claim validator and traceability checker (deep-research P6/P7 gates 3-5), team Mithril. Read-only
review; this file is the only write. Line numbers are those of the draft as read (654 lines). Evidence paths are
relative to the repository root unless they start with `notes/` (= `experiments/research_2026_11_sota/notes/`).

## Verdict

**Revise before commit.** Section 2's tables match `experiments/gate/results/path_to_target_2026-11-06.json`
exactly (28 of 28 table figures and every figure quoted in the text). Most literature and repository numbers trace to
an audited note, a design or critique file, or a committed file (about 75 checks pass, listed at the end). But:

- one citation is broken: the body's "[13 in task d, ...]" names a source the registry never approved, and the
  section 12 entry [13] is a different paper (Chreode) that no body sentence cites;
- the headline "+0.35" is the cross-board critic's sub-total, not the sum of the report's own surviving items
  (+0.39 to +0.43 on the report's own prices);
- two critics are misreported: the embryo critic's break-even is 0.3-0.6, not 0.2, and the disagreement is about the
  skill map's residual, not only the captain's prior; the heart critic calls the H-B premise "not new" and
  unidentifiable and offers the full hold as the one form;
- an audit-corrected overstatement comes back ("the organisers' own models sit at or below the floor");
- two anchors the report points to do not exist yet: `docs/DATA_POLICY.md` 2026-11-06 (a), and a committed
  `path_to_target.py` with its JSON (both untracked, as is the report's own folder).

## Issues, most severe first

### 1. [high] Section 5 "Gastrulation brackets", bullet 2, and section 12 entry [13]: broken citation

- **Sentences.** "Morphodynamic models of the E7.75-E8.25 interval are built on embryos inside the embryo window [13 in
  task d, cited for method only]." (REPORT.md:326-327). Entry "[13] Qiu M et al. Chreode (arXiv 2605.28111)"
  (REPORT.md:575).
- **What is wrong.**
  - "[13 in task d, ...]" is a non-registry citation form (gate 2).
  - task-d's local [13] is Raiola, Sendra, Domínguez, Torres, eLife 14:RP108559 (`notes/task-d.md`:23; finding 9 at
    :52). It has no global number: the map goes from `task-d [12] -> [90]` to `task-d [14] -> [91]`
    (`notes/registry.md`:283-284), and the registry's Stats line records "removed source lines skipped: 1". Its source
    line carries the screen's W redaction of the title.
  - Global [13] is Chreode (`notes/registry.md`:17, from task a), a different paper. No body sentence cites Chreode,
    so section 12's [13] has no body citation (gate 3), and the one body sentence that says "13" points outside the
    Approved list.
- **Fix.** Replace the sentence with one that [91] and task-d's mapping (e) support (`notes/task-d.md`:134, "Every
  morphodynamic model of the E7.25-E8.0 interval uses in-window embryos (Raiola, Esteban [13][14] ...)"):
  "Every morphodynamic model of the E7.25-E8.0 interval is built on embryos inside the embryo window [91] (cited for
  method only)." In section 12, delete the [13] Chreode line and add: "[91] Esteban et al. Pseudodynamic analysis of
  heart tube formation in the mouse (Nat Cardiovasc Res 1, abstract only). academic, 2022-05.
  https://doi.org/10.1038/s44161-022-00065-1". If Raiola itself is wanted, first add it to the registry's Approved
  list under a new number, with its screened title, and cite that number.

### 2. [medium] Executive summary, section 2 table, section 10 finding 4: "+0.35" is not the package's sum

- **Sentences.** "Every surviving item is worth +0.01 to +0.15 Overall. Together they come to about +0.35 Overall if
  every ruling goes as recommended." (REPORT.md:24-25); "| +0.35 (this research's package) |" (:137); "The surviving
  actions are worth about +0.35 Overall" (:537).
- **What is wrong.**
  - +0.35 is the cross-board critic's total of the cross-board items only: C2 +0.15, C4a +0.10, C1 revised +0.05, C3
    revised +0.05 (`notes/design/cross_critique.md`:228-230, "Cross-board total if every ruling goes as
    recommended: about +0.35").
  - The report's own surviving list (summary items 2-4, recommendations 2-4) also carries T1's cycling-state quota
    (+0.05; `notes/design/t1_critique.md`:14) and the heart hold (+0.03; `notes/design/heart_critique.md`:35). It
    prices the DM at +0.01 to +0.05 (`notes/design/t3_critique.md`:12), where the +0.35 uses +0.05.
  - On the report's own prices: 0.01-0.05 (DM) + 0.05 (runbook) + 0.15 (embryo s2) + 0.10 (A15) + 0.05 (cycling) +
    0.03 (heart hold) = **+0.39 to +0.43**.
  - The rows also apply expected values of bets (ranges down to -1.0) as sure shifts. `notes/design/cross.md`:459-464
    priced a similar set as two-point bets and got a different P(>= q95).
- **Fix.** Either:
  - (a) write "The cross-board items come to about +0.35 (`notes/design/cross_critique.md`); with T1's cycling-state
    quota (+0.05) and the heart hold (+0.03), the surviving set is about +0.4 (+0.39 to +0.43)". Then re-run
    `path_to_target.py` with 0.43 in `OVERALL_GAINS`, plus about 0.9 for the R1 + sheet C row (0.39-0.43 + R1's 0.2 +
    sheet C's 0.29; `notes/design/cross.md`:411-412). Replace 11.3%/5.9%/0.59% and 13.2%/7.3%/0.89% with its outputs;
    or
  - (b) keep +0.35 but call it "the cross-board items" at :24-25, :137 and :537.

  In both cases, head the section 2 column "Expected gain, applied as a sure shift".

### 3. [medium] Embryo s2: the threshold and the critics' disagreement are misreported

- **Sentences.** "the embryo s2 at mixture fraction 1/4: about +0.15, and only if he holds P(E7.75 is E7.5-like) at
  or above about 0.2" (REPORT.md:57-58); "worth it only if he himself puts P(E7.75 is E7.5-like) at about 0.2 or
  more" (:508); "The disagreement is about the captain's own prior, which only he holds." (:522).
- **What is wrong.**
  - 0.2 is the break-even of the pure-mixture toy, which the cross critic kept (`notes/design/cross_critique.md`:131,
    "Break-even stays near 0.15-0.2 (REPORT)"; :243-244).
  - The embryo critic calibrates the skill map to committed skill levels (matched-mixture skill about 0.67-0.7). That
    puts break-even at about 0.3-0.6 (`notes/design/embryo_critique.md`:12, table at :101-105, :272).
  - It also recommends adoption "only if he states sigma >= 0.15 day (gestational labels) and accepts the
    low-residual end" (:158).
  - So the critics differ on the skill map's non-mixture residual, a modelling input the committed files inform, and
    not only on the captain's prior.
- **Fix.**
  - At :57-58 and :508: "... only if he puts P(E7.75 is E7.5-like) at about 0.2 or more on the pure-mixture toy (cross
    critic), or at about 0.3-0.6 or more once the residual that committed skills show is allowed (embryo critic,
    section 3)".
  - At :522: "The critics differ on the skill map's non-mixture residual (embryo critic: break-even 0.3-0.6; cross
    critic: 0.15-0.2). Given that, the call rests on the captain's prior."

### 4. [medium] Section 8 "Heart" bullets, recommendation 4 and summary item 4: the heart critic is misreported

- **Sentences.** "The designer's new point stands: the part of B0's momentum shared across lineages carries the
  reference run's per-gene offset, which enters the test with the wrong sign. That is the reason for H-B."
  (REPORT.md:461-462); "Hold the heart's momentum at test (or keep it for cardiomyocytes only): about +0.03." (:515).
- **What is wrong.**
  - The critic says the mechanism is not new. The 2026-11-05 dynamics note already wrote the across-break increment
    with the run term inside e, and sheet E was priced on it (`notes/design/heart_critique.md`:147-150).
  - It also says "The premise under H-A to H-D cannot be identified in these files" (:22-23; section 5 at :280-301).
  - Its revision is to "Offer the full hold as the one form" (:186). The CM-only variant rests on an uncommitted
    carry-over line that may be scRNA-derived, which A15's working answer keeps out of uploads (:167-171).
- **Fix.**
  - Replace :461-462 with: "The designer's premise is not new. The premise: the lineage-shared part of B0's momentum
    carries the reference run's per-gene offset, with the wrong sign at test. The 2026-11-05 dynamics note and sheet E
    already priced it. The critic also finds it unidentifiable here: the shared part also carries admixture,
    depth-linked detection, shared biology and, at test, the E10.5 dissection change."
  - At :515, write: "Hold the heart's momentum at test (the full hold; sheet E's default stays no change): about
    +0.03."
  - Drop the CM-only option from the H-B row (:455).

### 5. [medium] Header: the process-notes anchor does not exist

- **Sentence.** "Process notes are in `docs/DATA_POLICY.md` 2026-11-06 (a)." (REPORT.md:19).
- **What is wrong.**
  - `docs/DATA_POLICY.md` has no 2026-11-06 entry. A grep for "11-06" finds nothing, and the last change note is
    2026-11-05 (k).
  - The header's counts are right. I recounted them from the notes' logs:
    - class A 6 (task e);
    - class W 24 (a2 1, b 6, c 6, d 5, e 6);
    - class L 1 (task e);
    - 5 audit removals (a 1, b 2, d 1, e 1);
    - 2 sources withdrawn (task b [22], task c [27]).
  - But the file they point to holds none of it. The notes also record slips that such an entry should carry:
    - the task-a and task-e screeners each opened a memory file under `~/.claude-account` (`notes/task-a.md`:232,
      `notes/task-e.md`:304);
    - the cross critic's grep displayed part of a practice-grid footer line (`notes/design/cross_critique.md`:274-276).
- **Fix.** Before committing the report, write the 2026-11-06 (a) change note in `docs/DATA_POLICY.md` with those
  counts and slips. Otherwise, change the sentence to: "Counts from the screening and audit logs of `notes/task-*.md`;
  the process notes go into `docs/DATA_POLICY.md` with this report's commit."

### 6. [medium] Summary item 1 and recommendation 1: "Done" is not yet true

- **Sentences.** "**Done, operational:** commit `experiments/gate/path_to_target.py` ..." (REPORT.md:50-51); "Commit
  `experiments/gate/path_to_target.py` and its output (done with this report)." (:491).
- **What is wrong.**
  - `GIT_OPTIONAL_LOCKS=0 git status --short` shows three untracked entries:
    - `?? experiments/gate/path_to_target.py`;
    - `?? experiments/gate/results/path_to_target_2026-11-06.json`;
    - `?? experiments/research_2026_11_sota/`.
  - The untracked research folder (REPORT.md and notes/) lies outside `results/`. It would trip the same clean-tree
    guard (`src/vec/submit.py`:8-10).
  - The cross critic asked to commit the script "with a review" (`notes/design/cross_critique.md`:31, :89-90). The T1
    critic called both files not committed evidence (`notes/design/t1_critique.md`:202-203). No review is named.
  - Section 2's "High for the arithmetic ... (19 checks passed)" rests on two things: the estimator's checks of its
    inputs, and the script's 191.51 reproduction (`path_to_target.py`:41-42). Neither checks the new outputs (the
    conditional means and the sure-gain rows).
- **Fix.**
  - Commit the script, its JSON and the whole `experiments/research_2026_11_sota/` folder in one commit. Until that
    commit exists, write "To do with this report's commit" instead of "Done".
  - Name a reviewer of the script's new outputs, or set section 2's confidence (:160) to "Medium-high: the 19 input
    checks and the 191.51 reproduction pass; the new outputs have no independent review."

### 7. [medium] Executive summary "The organisers' own models": an audit-corrected claim comes back

- **Sentence.** "**The organisers' own models** sit at or below the floor on the practice boards. Navigo needs
  spliced and unspliced counts that our files do not have [105][112]." (REPORT.md:44-45).
- **What is wrong.** task-e's audit corrected exactly this claim: "the embryo-interpolation row, 55.22, sits above
  the floor of 50, so 'at or below the floor' holds for T1, heart extrapolation and the Gata4 board only". The same
  audit found that "the reference-rows page [6] as rendered shows none of them" (`notes/task-e.md`:47, spot-check 1;
  task-e [6] is global [112]). The figures come from the screened digest, section 2.7
  (`docs/research/digest_2026-10-19.md`:529-530).
- **Fix.** "**The organisers' published reference rows** sit at or below the floor on T1, heart extrapolation and the
  Gata4 board; their embryo row (55.22) is above it (`docs/research/digest_2026-10-19.md` 2.7). Navigo needs spliced
  and unspliced counts that our files do not have [105]."

### 8. [low] Section 4: the few-shot range is tied to the wrong k

- **Sentence.** "With 2-5 measured responses in a new context, the best fixed rule is 0.3 shrinkage, and gains stay
  tiny (0.0128 to 0.0179) [53]." (REPORT.md:246-247).
- **What is wrong.** In [53], 0.0128 to 0.0179 is the oracle-weighted gain from k = 2 to k = 40. The 0.3 rule gives
  0.0148 at k = 5 (`notes/task-c.md`:39, finding 3; audit :318).
- **Fix.** "Gains stay tiny even with many measured responses: oracle-weighted adaptation rises only from 0.0128
  (k = 2) to 0.0179 (k = 40), and the best fixed rule, 0.3 shrinkage, gives 0.0148 at k = 5 [53]."

### 9. [low] Recommendation 3: the deadline conflict is located incompletely

- **Sentence.** "Also reconcile `docs/PLAN.md` section 7 ("by 2026-12-01 06:00 UTC") with section 8 and S2_RULES ("by
  the end of 30 Nov IST")." (REPORT.md:505-506).
- **What is wrong.** Section 8 holds both dates:
  - its deadline-exception bullet says "The s2 round follows by 2026-12-01 06:00 UTC" (`docs/PLAN.md`:630-631);
  - its later bullet says "Final uploads still finish by the end of 30 November IST" (:651), as `docs/S2_RULES.md`
    does (:49, :282).

  Fixing section 7 (:603-604) alone leaves section 8 inconsistent. The wording is inherited from
  `notes/design/cross_critique.md`:84-86.
- **Fix.** "Also reconcile `docs/PLAN.md` section 7's 2026-11-05 status note and section 8's deadline-exception bullet
  (both "by 2026-12-01 06:00 UTC") with section 8's "Final uploads still finish by the end of 30 November IST" and
  `docs/S2_RULES.md` ("by the end of 30 Nov IST")."

### 10. [low] Section 11 "they agree", and section 3's LGP-OT figure

- **Sentences.** "Both were screened and audited, and they agree (section 3)." (REPORT.md:552); "PACE (NeurIPS 2026)
  and LGP-OT report gains over learned models only (23.7% and 1-4%)" (:182).
- **What is wrong.** The two task-a notes agree on conclusions but differ on figures the report uses:
  - **LGP-OT.** task-a gives "1-4% W2 ... (31.88 against 32.80 at t=8)" (`notes/task-a.md`:46). task-a2 gives "at
    most 2% W2 (32.15, 32.52, 36.51 against 32.80, 32.88, 37.06)" (`notes/task-a2.md`:42). Neither was spot-checked,
    and the report takes 1-4% silently.
  - **scFM horizon.** task-a, as audit-corrected, gives "1-4 steps" (`notes/task-a.md`:35); task-a2's mapping A gives
    "1-3 steps" (`notes/task-a2.md`:85).
  - **count-FM comparators.** task-a gives "20.739 and 22.226" (`notes/task-a.md`:43); task-a2 gives "20.545-22.226"
    (`notes/task-a2.md`:52).
  - **USB.** task-a says "interpolation only" (`notes/task-a.md`:48); task-a2 says "no held-out evaluation in their
    abstracts" (`notes/task-a2.md`:50).
- **Fix.**
  - At :552: "Both were screened and audited. They agree on every conclusion and differ on some figures (LGP-OT's
    margin, scFM's horizon, count-FM's comparators); section 3 uses the audit-corrected figure where one exists."
  - At :182: "(23.7%; for LGP-OT 1-4% in task a and at most 2% in task a2)".

### 11. [low] Section 5: stVCR's baselines

- **Sentence.** "*stVCR (Nature Methods 2026).* Learns rigid transforms with unbalanced dynamical OT, checked against
  learned baselines only [81]." (REPORT.md:302-303).
- **What is wrong.** task-d finding 3 lists stVCR's baselines as "TIGON, SF2M and static OT" (`notes/task-d.md`:46).
  Static OT is not learned. What the finding supports is "no copy-nearest baseline".
- **Fix.** "... checked by leave-one-time-point-out interpolation against TIGON, SF2M and static OT, with no
  copy-nearest baseline [81]."

### 12. [low] Section 5: the Dice 0.79 source

- **Sentence.** "... and a core crop alone moves Dice to 0.79 (`experiments/coords/REPORT.md`)." (REPORT.md:320-321).
- **What is wrong.** `experiments/coords/REPORT.md` holds the 63% shell (0.630 at :17 and :45) and the 30-40 sections
  (:49), but no Dice of 0.79. The core-crop Dice is in `docs/research/digest_2026-10-19.md` 4.4 (:917-918), which
  task-d's [29] also names (`notes/task-d.md`:41).
- **Fix.** "(`experiments/coords/REPORT.md`; `docs/research/digest_2026-10-19.md` 4.4)".

### 13. [low] Section 3 "Blood": the 24% is not in row R1

- **Sentence.** "blood is 2.6% of chamber cells in the E9.5 inputs of both counted tasks, against 24% at E10.5 (TRIES
  row R1)." (REPORT.md:225-226).
- **What is wrong.** Row R1 gives the 2.6% only (`experiments/harness/TRIES.md`:260). The 24% is in `docs/PLAN.md`
  4.1 (:190, :216) and in brief section 3.
- **Fix.** "(TRIES row R1; `docs/PLAN.md` 4.1 for E10.5's 24%)".

### 14. [low] Executive summary, T1 bullet: the denoising reading is attributed to [2]

- **Sentence.** "Its margin matches what denoising alone buys in high-dimensional W2, and our variogram penalises
  denoising [2]." (REPORT.md:33-34).
- **What is wrong.** The denoising reading is task-a2's own derivation ("my derivation, not the authors'",
  `notes/task-a2.md`:41). The variogram cost comes from the screened digest, section 2.1
  (`docs/research/digest_2026-10-19.md`:120-121). [2] supplies only the observed ratios. Section 3 attributes this
  correctly (:173, "Task a2's reading").
- **Fix.** "... matches what denoising alone buys in high-dimensional W2 (task a2's reading of [2]'s figures), and
  our variogram penalises denoising (`docs/research/digest_2026-10-19.md` 2.1)."

### 15. [low] Section 3 "Blood": two anchors were not verified in audit

- **Sentences.** "no primitive enucleation at or before E9.5, and more than 95% enucleated by E15.5 [40]"; "definitive
  cells are more than 94% of circulating cells at E16.5 [39]" (REPORT.md:222-223).
- **What is wrong.** task-b's spot-check 4 on [40] is "not verifiable": no allowed route served the full text, and
  the abstract gives no stage-specific fraction. "The [9] part was not checked" (task-b [9] is global [39];
  `notes/task-b.md`:42). The report states both as fact.
- **Fix.** Append "(not verified in audit; task-b spot-check 4)" to both lines, or drop both figures and keep the
  verified [38] anchor.

### 16. [low] Section 1: the audit outcomes are miscounted in wording

- **Sentences.** "spot-checked 2-4 claims per note" (REPORT.md:91); "19 levers flagged as repeating a pruned item"
  (:106).
- **What is wrong.**
  - All seven audits ran 4 spot-checks each: 28 in all, and 17 verified + 9 corrected + 2 not verifiable = 28, which
    matches the report.
  - The 19 counts duplicate flags, not levers:
    - task-a's 4 marks fall on 2 levers (`notes/task-a.md`:258, "4 marks on 2 levers");
    - task-e's one flag is mapping 8 (Squidiff), not a lever (`notes/task-e.md`:359-360).

    By lever, the count is 16.
- **Fix.** "spot-checked 4 claims per note (28 in all)"; "19 duplicate flags (on 16 levers and one mapping entry)".

### 17. [low] Sections 4 and 7: "rank-based except the slope" overgeneralises its source

- **Sentences.** "de_score is rank-based, so of T3's four metrics only the severity slope rewards magnitude
  calibration." (REPORT.md:286); "our T3 metrics are rank-based except the slope" (:391).
- **What is wrong.** The source limits this to the DE family: "in VEC the only scale-sensitive DE-family metric is
  T3's severity_slope; de_score and de_direction are rank-based and immune" (`notes/task-f.md`:78). T3's fourth term,
  mutant state (0.20: mmd_u 0.6, variogram 0.4; brief section 2), is a cell-level distance. It is neither rank-based
  nor scale-free.
- **Fix.** At :286, "of T3's three DE metrics only the severity slope rewards magnitude calibration". At :391, "our T3
  DE metrics are rank-based except the slope".

### 18. [low] Section 2 "Opposing reading": stronger than its source

- **Sentence.** "Under task h's measurement-based T3 payoffs (s1's best case about 66; ... section 1), even our 95th
  percentile is out of reach." (REPORT.md:163-164).
- **What is wrong.** The source says the chance "is under 1%" (`experiments/research_2026_11/REPORT.md`:44-45; table
  at :108). It does not say the percentile is out of reach.
- **Fix.** "... the chance of reaching even our 95th percentile is under 1%."

### 19. [low] Section 2: reading A against the outlook

- **Sentence.** "Reading A is the upper reading: the committed outlook is 184, because reading A credits T1's
  oracle-base grid in full." (REPORT.md:117-118).
- **What is wrong.** The committed gap is +7.51, split across the boards (`experiments/gate/results/estimate_2026-11-05.json`,
  `A_minus_outlook`):

  | Board | Gap |
  |---|---|
  | T1 | +6.57 |
  | heart | +1.85 |
  | embryo | +0.87 |
  | T3 | -0.42 |

  The README says "(c) does not say how it reached 63" (`experiments/gate/README.md`, Result). The "because" explains
  only T1's share, and it is an inference.
- **Fix.** "Reading A sits 7.5 above the committed outlook of 184, 6.6 of it on T1, where reading A credits the
  oracle-base grid in full (`experiments/gate/README.md`, Result)."

### 20. [low] Executive summary, post-mortems: [148] is the wrong source

- **Sentence.** "Final picks that chased the public score fell, and robust validation rose [148][158]." (REPORT.md:47).
- **What is wrong.** [148]'s finding is that the public-test score was the best predictor of private rank (Spearman
  0.944), with one entry falling from 2nd to 133rd (`notes/task-f.md`, finding 8). Task f's score-chasing finding
  (finding 14) cites task-f [24][25][6], which are global [158][159][141], not [148].
- **Fix.** "... rose [158][141]". If wanted, add "(one 2022 entry fell from 2nd to 133rd [148])".

### 21. [low] Section 2 path, and the executive summary's source count

- "(`results/path_to_target_2026-11-06.json`)" (REPORT.md:114). No file exists at that path from the repository root.
  The file is `experiments/gate/results/path_to_target_2026-11-06.json`. Fix: give the full path.
- "this one, with 163 sources" (:28). The 163 registry entries include four duplicate works, so there are about 159
  distinct works:
  - [1]/[26] scTimeBench;
  - [13]/[127] Chreode;
  - [33]/[102] mdCLA;
  - [104]/[126] MERFISH+.

  Fix: "163 registry entries (about 159 distinct works)". The header's "163 registry entries" (:4) is accurate as
  written.

## Gate 4: section 2 against `path_to_target_2026-11-06.json`: no issue

| Report (section 2) | JSON key | JSON value | Match |
|---|---|---|---|
| mean 191.5 | mean_overall | 191.4969 | yes |
| 90th / 95th / 99th: 203.95 / 205.54 / 208.07 | levels.own_q90/95/99 | 203.9529 / 205.5388 / 208.0747 | yes |
| board means T1 69.6, heart 53.9, embryo 60.9, T3 64.6 | board_means | 69.5715, 53.8562, 60.8620, 64.5663 | yes |
| given >= q95: T1 72.4, heart 58.5, embryo 62.3, T3 74.3 | board_means_given_own_q95 | 72.3738, 58.4786, 62.3326, 74.2990 | yes |
| none: 10% / 5% | p_at_least.own_q90 / own_q95 | 0.1 / 0.05 | yes |
| +0.35: 11.3% / 5.9% | whatif "0.35" | 0.11305 / 0.058955 | yes |
| +0.8: 13.2% / 7.3% | whatif "0.8" | 0.13172 / 0.07259 | yes |
| +1: 14.0% / 7.9% | whatif "1" | 0.140055 / 0.079145 | yes |
| +3: 23.2% / 15.8% | whatif "3" | 0.232305 / 0.158235 | yes |
| +5: 31.8% / 25.1% | whatif "5" | 0.318405 / 0.251285 | yes |
| +10: 47.3% / 43.0% | whatif "10" | 0.472715 / 0.430175 | yes |
| P(>= 209): 0.42%, 0.59%, 0.89%, 9.8%, 31.7%, 49.7% | p_at_least.target_209; whatif 0.35/0.8/5/10/16 p_target_209 | 0.004225, 0.005925, 0.008935, 0.09832, 0.316635, 0.49724 | yes |
| "about +16 ... even bet" (:27, :154) | whatif "16"; whatif_p_target_after_sure_shift.T1.d_for_p_half | 0.49724; 16.105 | yes |
| 200,000 draws | estimate.py N_DRAWS | 200_000 | yes |
| "19 checks passed" | checks | 19 entries | yes |
| "+10 heart points is +5 Overall" | estimate.py OVERALL_WEIGHT heart | 0.5 | yes |
| "(T3 about 75, T1 about 73)" | board_means_given_target_209 | T3 75.06, T1 73.03 | yes |

The "near their own 90th percentiles" part of :158 holds against `estimate_2026-11-05.json` (readings A band 90th:
T1 73.01, T3 74.62). Cite that file for it.

## Gate 3: section 12 against the body

- The body cites 89 distinct registry numbers, and section 12 has 90 entries. Every body number has an entry.
- Entry [13] has no body citation, and one body citation is in non-registry form (issue 1).
- Every entry's URL equals its registry URL (a scripted comparison found no mismatch).

## Gate 2: [n] against the registry map

I checked these numbers against `notes/registry.md`'s map and against the finding or deep read that carries the
claim. Each is the right source for its sentence, except as noted:

- **Task a and a2.**
  - [1]: a[1].
  - [2]: a[2], a2[1].
  - [3]: a[3].
  - [4]: a[4], a2[16].
  - [5]: a[5].
  - [6]: a[6].
  - [7]: a[7], a2[5].
  - [8]: a[8], a2[6].
  - [9]: a[9].
  - [10]: a[10].
  - [15]: a[15].
  - [16]: a[16], a2[4].
  - [18], [19]: a2[13], a2[14].
  - [21], [22], [23], [24]: a[21]-[24].
  - [27]: a2[8].
- **Task b.**
  - [31]: b[1].
  - [33], [34]: b[3], b[4].
  - [38], [39], [40]: b[8]-[10].
  - [45], [46], [47]: b[15]-[17].
  - [50]: b[20].
- **Task c.**
  - [53]: c[1].
  - [55]-[61]: c[3]-c[9].
  - [63], [64], [65], [66]: c[11]-c[14].
  - [70]-[74]: c[18]-c[22].
  - [78]: c[26].
- **Task d.**
  - [79]-[82]: d[1]-d[4].
  - [84], [85], [86]: d[6]-d[8].
  - [88]: d[10].
  - [94]: d[17].
  - [96], [97], [98], [99]: d[19]-d[22].
  - [103]: d[26].
  - [105]: d[28], e[16].
- **Task e.**
  - [108], [109]: e[2], e[3].
  - [111], [112]: e[5], e[6].
  - [114]-[117]: e[8]-e[11].
  - [122]: e[17].
  - [128], [129], [130]: e[24]-e[26].
  - [133]: e[29].
- **Task f.**
  - [136]: f[1].
  - [140], [141], [142]: f[5]-f[7].
  - [144]: f[9].
  - [146]: f[11].
  - [148], [149]: f[13], f[14].
  - [151]: f[16].
  - [155], [156]: f[20], f[21].
  - [158], [159], [160]: f[24]-f[26].

  [160] backs only "third place"; the method itself is [136].

Wrong-source cases:

| Citation | Where | Issue |
|---|---|---|
| [13] | :327 | issue 1 |
| [112] | :45 | issue 7 |
| [2] | :34, for the denoising derivation | issue 14 |
| [148] | :47 | issue 20 |

## Gate 1: claims checked that pass (no issue)

| # | Claim (REPORT line) | Source | Verdict |
|---|---|---|---|
| 1 | 163 entries, 41 domains, 114/28/21, 87%, biorxiv 23%, official 17% (:4-6) | registry.md Stats :361-366 | no issue |
| 2 | A 6, W 24, L 1, 5 residual, 2 withdrawn, no class B (:13-19) | notes' screening and audit logs | counts right; anchor missing (issue 5) |
| 3 | 17 verified, 9 corrected, 2 not verifiable; 32 legality corrections (:105) | notes' audit logs (4+4+4+4+4+4+4; legality 4+5+4+4+4+7+4) | no issue |
| 4 | scFM 27-36% below copy-last, 1-4 steps, scNODE within 1-10% (:171-172) | task-a F1 (audit-corrected) | no issue |
| 5 | ratio near 0.71 against observed 0.64-0.73 (:173-174) | task-a2 F2, deep read [1] | no issue |
| 6 | pseudobulk_shift loses about 2.2 points to the variogram (:176) | digest 2.1, :120-121 | no issue |
| 7 | scTimeBench ten methods, eight datasets, one step, Squidiff first on MMD and last elsewhere (:179-180) | task-a F2 | no issue |
| 8 | CellTempo PCC 0.74 vs 0.63 (:181) | task-a F3 | no issue |
| 9 | PACE 23.7%, no naive baseline (:182) | task-a F12 | no issue (LGP-OT: issue 10) |
| 10 | WFR-FM d/dt ln m from a semi-coupling; 1-4% W1, one loss (:185-186) | task-a F7 and deep read [7]; task-a2 mapping :94 | no issue |
| 11 | ARTEMIS 4 of 6 by at most 0.5 W; scMTG no ablation; USB and SUDO (:187-189) | task-a2 F11 | no issue |
| 12 | count-FM 20.456 vs 20.545-22.226; scLDM 7.267 vs 11.658 (:193-194) | task-a2 F13 | no issue |
| 13 | 1.52 vs 8.03 weeks; 11 studies; 2.65 weeks (:197-198) | task-a F6 (audit-corrected) | no issue |
| 14 | MASE 0.78 above consistency 0.8; none beats carry-forward below 0.2 (:199-200) | task-a F5 | no issue |
| 15 | TabPFN Spearman 0.808, R2 0.704; mean baseline 0.352 (:201-202) | task-a F4 | no issue (0.352 is an R2) |
| 16 | TabPFN-2.6 licence bars commercial or production use (:202-203) | task-a2 F15 (audit-corrected) | no issue |
| 17 | Leonard 12,507 nuclei, four groups, about 29% CM, 26% fibroblasts, GSE272154 (:206-207) | task-b F1, spot-check 1 | no issue |
| 18 | mdCLA 351 embryos, 37 points, 10 late plus birth, 18 at or before E9.5 (:208-210) | task-b F2, spot-check 2 | no issue |
| 19 | Feulner Drop-seq E14.5-P7, six replicates, E14.5 out (:211-212) | task-b F4 | no issue |
| 20 | 75-80% vs "nearly half"; about 23% vs about 9% (:215-216) | task-b F14 | no issue |
| 21 | 0.84-0.86 vs 0.90-0.95 (:217-218) | task-b F15 | no issue |
| 22 | nucleated fraction about 0.34 to 0.01 across the E14.5 litter harvest (:224) | task-b F9, spot-check 3 | no issue |
| 23 | 27.8 / 29.4 / 23.5 / 19.3%; template 23-65%; 92-98% (:238-240) | task-c F1, spot-check 2 | no issue |
| 24 | r = 0.25-0.39; State 12%; PerturbMap 4.1% (:241-243) | task-c F5 | no issue |
| 25 | 15 x 88 atlas, median 48%, 0.590 to -0.021 (:244-245) | task-c F2, spot-check 1 | no issue |
| 26 | scPT-seq; 109 of 406 vs 0 of 50 (:250-253) | task-c F9 | no issue |
| 27 | MIC-Drop-seq: abundance of 1 type, transcription of 12, extrinsic more than half (:254-255) | task-c F7, spot-check 4 | no issue |
| 28 | ZSCAPE r = 0.62 (:256) | task-c F8 | no issue |
| 29 | 641 sections, about 2-fold, 31%, 60-104 genes; "massive" loss; chemistries (:259-262) | task-c F11, F12, spot-check 3 | no issue |
| 30 | admixture [70]; normalisation false DE [74] (:263-264) | task-c F13 | no issue |
| 31 | Klaus phenotype (:276-278) | task-c F15 and deep read [26] | no issue |
| 32 | split-half cosine 0.79-0.95; [57]; [58][59] (:282-285) | task-c F14 | no issue |
| 33 | ContextFlow about 10% weighted W2, expression only, no copy, identity or linear baseline (:299-300) | task-d F1 | no issue |
| 34 | TP-DATE beats FGW-CFM; no closed-form real-cell FGW baseline anywhere (:304-305) | task-d F4 and audit :261 | no issue |
| 35 | Mimyr 0.39-0.48 to 0.45-0.58 on held-out slices (:306-307) | task-d F5, spot-check 1 | no issue |
| 36 | 97-98% same-section neighbours; 1.5x z costs 12.6-12.8 (:311-313) | coords REPORT :19; digest :920-921; task-d F13 | no issue |
| 37 | DMDD 215 WT C57BL/6N E14.5, 3 um; E15.5 micro-CT; EMA TS23-TS24 (:317-319) | task-d F11, spot-check 4 | no issue |
| 38 | 63% shell at E9.5; 30-40 sections (:320-321) | coords REPORT :17, :45, :49 | no issue (Dice: issue 12) |
| 39 | SEU-3D only E8.0; Harland only at or before E7.25; Flysta3D-v2 (:325, :328) | task-d F14, F15 | no issue |
| 40 | embryo critic's at most +0.02 for the asymmetry check (:329-330) | embryo_critique.md :15, :223 | no issue |
| 41 | E8.25_late includes yolk sac; E10.5 77% surrounding tissue (:335-336) | heart_critique :92-95; PLAN.md :377; S2_RULES :185 | no issue |
| 42 | Task 1 page quote; neural_ode does not beat the constant shift (:340-343) | task-e F1, F2 | no issue |
| 43 | 48.58 / 47.76 / 55.22; 42.48 / 42.50 on the Gata4 board (:344-346) | digest :529-530; task-e F3 | no issue in section 6 (summary: issue 7) |
| 44 | Navigo 12.4M cells, 43 time points, 20-fold; spliced and unspliced; no proliferation; interpolation only (:349-352) | task-e F4, F5 | no issue |
| 45 | T3 page pitfalls (:354-356) | task-e F8, spot-check 2 | no issue |
| 46 | veckit cc3859a of 2026-10-02; PR says the server checks finite and non-negative; no change after 2026-10-28 (:360-362) | task-e F10, F11, spot-check 4 | no issue |
| 47 | proposal: distribution term on mutant state; paired bootstrap (:364-365) | task-e F12 | no issue |
| 48 | VCC 2025: H1, 300 perturbations, more than 1,200 teams; winners' methods (:381-384) | task-f F1 | no issue |
| 49 | clipping; top-50 below zero; DES about 40% vs 39% (:385-387) | task-f F2, F4 | no issue |
| 50 | 1970 regression 15th of 337; kernel ridge 38th; 0.3413 vs 0.3401 (:388-390) | task-f F5 | no issue |
| 51 | 2022: +0.06 / +0.03; Spearman 0.944; 2nd to 133rd (:396-397) | task-f F7, F8 | no issue |
| 52 | 2021: up to 8 ranks; 2023 top four; obesity MMD runner-up (:398-402) | task-f F9, F10, F12 | no issue |
| 53 | M5 1,000-5,000 places; BirdCLEF; Turbine 15th (:405-408) | task-f F14 | no issue |
| 54 | VCC 2026: six metrics, three lines, mid-to-late November (:411-412) | task-f F6 (scoring details unverified; the rest verified) | no issue |
| 55 | 28 candidates (:420) | 5 + 4 + 6 + 5 + 8 in the design files | no issue |
| 56 | every price in the five section 8 tables | t1, t3, heart, embryo and cross critiques' verdict tables; t1.md :393 for E's range | no issue |
| 57 | +0.8 = +0.35 plus R1 and sheet C as the design counts them (:25-26) | cross.md :411-412; cross_critique :230-231 | no issue |
| 58 | clean-tree guard: untracked files outside results/, runs/ and experiments/**/results/ (:492-493) | src/vec/submit.py :8-10 | no issue |
| 59 | DM route in ORGANIZER_QUESTIONS; no rules answer since 2026-10-16 (:494-498) | ORGANIZER_QUESTIONS.md :3; DATA_POLICY.md :269; t3_critique :146-149 | no issue |
| 60 | dated runbook (:500-503) | cross_critique :82-84 | no issue |
| 61 | C1 bytes-per-cell reasoning; reverses 2026-11-01 (b); misread costs -6 to -8 (:446-447, :516-517) | t3_critique :36-40, :75-79, :112-113 | no issue |
| 62 | 209 was the leading Human total when set (:144-145) | PLAN.md :618 | no issue |
| 63 | R1 kept 2026-11-05, not in the estimate (:25) | PLAN.md :642-643 | no issue |
| 64 | sheet C is the 2026-11-05 report's T1 strong-null bet (:26) | experiments/research_2026_11/REPORT.md :481-484 | no issue |
| 65 | four earlier searches found +1 to +3 (:28, :156-157) | brief.md :10-12 | no issue |
| 66 | ICLR 2027 unread because of an OpenReview browser check (:108-109, :545) | task-a :196; task-a2 :204 | no issue |

The following are process claims with no line in the evidence set, so I did not count them as issues: "a second
independent pass on task a (a2, created when the workflow restarted task a)" (:88-89, :551-552) and "The main agent
read the audited notes only" (:93). If the workflow log is to stand behind them, cite its id.

## Process and compliance of this review

- **Read.**
  - REPORT.md; notes/registry.md (all); notes/brief.md; all ten files in notes/design/ (the critiques in full, cross.md
    sections 2-5, heart.md section 3, t1.md for E's range).
  - From the notes: the Findings, Deep Read Notes, Screening logs and Audit logs of task-a, a2, b, c, d, e and f, by
    grep and line ranges.
  - In the repository:
    - `experiments/gate/path_to_target.py` and its JSON; `experiments/gate/results/estimate_2026-11-05.json`;
      `experiments/gate/README.md` (reading rules and Result); the `estimate.py` constants;
    - `docs/PLAN.md` sections 7-8 and decision log 2026-11-06; `docs/DATA_POLICY.md` change-note headings;
      `docs/S2_RULES.md` (grep); `docs/ORGANIZER_QUESTIONS.md` (grep); the TRIES R1 row (grep);
    - `experiments/coords/REPORT.md` (grep); `src/vec/submit.py` (docstring lines 8-10);
      `experiments/research_2026_11/REPORT.md` summary and section 1;
    - `docs/research/digest_2026-10-19.md`, the screened digest, as named lines only.
- **Git.** `GIT_OPTIONAL_LOCKS=0 git status --short` and `GIT_OPTIONAL_LOCKS=0 git log --oneline -8` only. No fetch,
  pull or push.
- **Computation.** Stdlib `python3 -I -B` one-liners that parsed the two gate JSONs, the registry and the report's
  citations. No repository code was run.
- **Not opened.** `boards/`, `data/`, `runs/`, `docs/research/` other than the digest, the workflow's raw/ and
  screened/ copies, anything under `~/.claude-account` or `~/.claude`, and the web. The harness put the auto-memory
  index into this session's context unasked; it was not opened or used.
- **Display note.** One `sed` of digest lines 525-534 displayed screened validation-phase overall totals of other
  teams, a kept class in the screened digest. They were not used and are not copied here. No class A, B, W, L, Y or
  Q item was met.
