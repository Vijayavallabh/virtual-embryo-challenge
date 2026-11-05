# T3:bcatenin: adversarial critique of the design of 2026-11-06

Critic for team Mithril, 2026-11-06. This reviews `design/t3.md` (candidates C1-C4) against the audited notes, the
brief and the committed records. It is design review only: no run, no score, no repository edit, no web. T3 counts 1 in
Overall, so T3 points and Overall units are the same here.

## Verdicts

| # | Candidate | Verdict | Designer's mean | Corrected mean (Overall) | Main reason |
|---|---|---|---|---|---|
| C1 | Screened release of the Data Protocol page's test-knockout acquisition text as s1's depth premise | RULING-ONLY | +0.2 (0 to +0.5) | about +0.05 (-0.05 to +0.15) | The price uses P(depth mismatch) where it needs P(the text triggers branch (b)). T3 depth is exposure-adjacent. The branch (b) action is wrong for one of its triggers |
| C2 | A second channel for the floored-severity answer (A24) | KEEP-REVISED | +0.05 to +0.15 | +0.01 to +0.05 | Nearly free and legal. The DM channel already went unanswered once |
| C3 | Section-artifact audit and a section-robust source factor | KILL | about +0.05 | about 0 | Rule 5: the expected gain on the folds is far below 0.3. As written it would break the Gata4-free hedge |
| C4 | Lineage-distance source factor | KILL | +0.03 to +0.05 | about 0 | It repeats a pruned item and fails rule 5. Its spatial selection is confounded and its control is ill-defined |

The corrected total is about +0.05 to +0.15 Overall, against the designer's +0.3 to +0.45. I agree with the
designer's main conclusion: nothing here moves T3's upper tail, and the open factor lever is close to empty.

## What I checked, and errors in the design text

These figures match the committed files:
- Reading A gives T3 64.58 with a band of 50.36-74.62, and the first round gives 50.84
  (`experiments/gate/results/estimate_2026-11-05.json`). The headroom rows match TRIES.
- `xrep:F` scored 68.05 / 64.09 and `xrep:F+depth` 74.50 / 72.15, so depth adds +6.45 / +8.06
  (`experiments/t3_decomp/results/gata4/gata4_summary.csv`).
- W1 gives 73.51 against 72.25 and 73.35. These are the means of `xrep:F+depth@w` on the replicate fold, which never
  counts.
- `vs95` scored 52.04 / 46.24 and `fac:none` 52.39 / 50.61.
- K2a's record says "change does not apply" on Mab21l2->Gata4rep2 (`experiments/harness/results/t3_k2_v1/`).
- `configs/pseudo_tasks.yaml` sets the folds' references.

Errors and overstatements:
1. **C2 misattributes a figure.** It says the swap's +1.9 to +2.5 points of top-3 chance hold "in the zero world (N10)".
   N10's figures are under the mixture of floor readings at lambda 1 (TRIES N10; REPORT section 3), not under the zero
   reading.
2. **C1 gives the wrong reason for the 2026-11-01 withholding.** It says the T3 text "stayed withheld on 2026-11-01,
   when T3 had no open decision it could inform". The record's reason is different: the text about the T2 stacks and
   the test knockout "stays withheld, because the T2 size and scope choices are pre-registered" (DATA_POLICY 2026-11-01
   (b)). T3's depth step is a pre-registered choice too (D4). So C1 asks the captain to reverse that reasoning for T3; it
   does not fill a gap.
3. **Section 2 says "At test the opposite holds".** The W1 evidence comes from the replicate fold, which shares both
   the reference and the knockout run (slide M164). At test only the reference is certainly shared, so the best
   exponent lies somewhere between what the cross folds and the replicate fold prefer. The control-arm logic does not
   need this claim.
4. **C1 relies on a builder flag that does not exist.** It says "one rebuild of s1 with the builder's depth flag off".
   `experiments/t3_decomp/build_submission.py` has no depth option. The depth step comes from the decision file, and only
   N4's diagnostic (`diagnostics.py t3-gate`) switches it off, through `spec["depth"]`. An upload file without depth needs
   a new build path in a hash-pinned builder, with tests, `vec check`, the server format check and the determinism check.
   Budget about a day, not hours. That still fits before 2026-11-27.
5. **C3 and C4 say they leave the Gata4-free hedge intact.** They do not (see C3).
6. **Section 5 judges the candidates against 209.** DATA_POLICY 2026-11-05 says the gate figure is never a threshold
   for a method choice, and the REPORT counter-review restated such distances against our own percentiles (reading A's
   90th-99th, about 204-208). No verdict in the design turns on it, so this is a wording fix.

## C1: organiser acquisition text as the depth premise (RULING-ONLY)

**(1) Repeat, closed lever or declined ruling.**
- C1 reopens the closed T3 depth line (TRIES rule 6; the QC20 entry) and departs from D4 (the frozen T3 rules). The
  designer says both.
- It also reverses the 2026-11-01 (b) decision to keep the test-knockout text withheld (error 2 above). The captain
  should see it as that reversal.
- The new information would come from the text itself, which is acceptable for a ruling request. The pattern "an
  explicit organiser statement switches a pre-registered T3 output" already has a precedent in the ruling of
  2026-11-06 (b).

**(2) Rules.**
- *Not class A in itself.* The text is a processing fact that the organisers published to every team, and their page
  "withholds every file property of a held-out target" (DATA_NOTES). T1's premise sets the precedent for protocol text.
- *T3 depth is exposure-adjacent, which the design misses.* The design says "the post-exposure block is untouched".
  Two exposures bear on this:
  - DATA_POLICY 2026-11-04 (b): the captain pasted the Data page into the main session, and it lists the two
    beta-catenin replicates' file sizes.
  - The first 2026-11-04 entry: the main session read the test boards' cell limits.

  (f)(3) closed T1 depth because "a cell count read with the ... file size gives bytes per cell, which bears on
  detection depth". The same reasoning reaches T3. N6 found every knockout file dense and every WT file sparse, so bytes
  per cell also show which pipeline wrote a file, and that is the fact C1 wants from the text. The block's text names
  T3 composition, not T3 depth, so C1 breaks no rule as written. But the captain has to rule on T3 depth after reading
  2026-11-04 (b) and (f)(3), and the person ruling can see the barred proxy for the same answer on the Data page.
- *Screen list gaps.* Two more items belong on the withheld list:
  - stage measurements of the test embryos (somite counts, morphological stage). The staging *method* may be released,
    but a measured stage would inform U3, another property of the hidden target;
  - slide and batch identifiers, which are file properties.

**(3) Evidence.** No fold can test C1, and the design prices it from committed records. I agree with that framing.

**(4) Expected gain.** The design computes 0.25 x 0.3 x (net +4), about +0.2. Four problems:
- **(a) The 0.3 is the wrong probability.** It is 1 - p_U2, the chance of a depth mismatch. The trigger is narrower:
  the text must say (i) the beta-catenin knockout was imaged with the WT controls, or (ii) it used a different chemistry,
  protocol or pipeline from Gata4.
  - "A separate run from Gata4" on the same chemistry falls in branch (c), which changes nothing. That is the commonest
    way U2 could fail: Gata4 and Mab21l2 already come from two independent runs (N6).
  - The released Data Protocol text already says, in general terms, that "Knockouts and their wild-type controls were
    processed in separate runs" (DATA_NOTES), so (i) is unlikely.
  - The page describes one pipeline for every file (Cellpose 2 on DAPI; log1p "checked on every released file").
  - So P((b) | the text speaks to the run relation) is about 0.1-0.2, not 0.3.
- **(b) A chemistry statement does not fix the direction of the depth mismatch.** An older protocol would leave the
  file shallower and a newer one deeper [task-c 19][task-c 21]. An independent knockout run on the same chemistry was
  still shallow (Mab21l2 at 0.73x WT E8.75 counts per object; N7). So P(the depth step still helps | (b)) is about
  0.25-0.35, not 0.2.
- **(c) The gain from a correct switch is optimistic.** The +6 to +9 rests on depth-free s1 passing the severity gate
  at test (bar 0.03 without depth, 0.22-0.66 with it). The only committed analogue is D11b's Mab21l2-against-WT-E8.75
  factor passing on Gata4 (severity skill 0.89). The 2026-11-05 review traced most of that carry to the WT
  E9.5-against-E8.75 offset, the part that couples the two bets. So it is one analogue, measured on shallow Gata4
  files, and not a measurement of the switch. Take +4 to +8.
- **(d) Branch (b) has one action for two different worlds.** If the beta-catenin file was imaged with the WT
  controls, the factor is noise too (task-c lever 1, risk line), and "drop depth, keep the factor" is the wrong action.

Recomputed: P(the text speaks to it) 0.15-0.25, times P((b)) 0.1-0.2, times a net of about +2 T3 (0.7 x +6 - 0.3 x 7;
range -1 to +4), times about 0.85. The last factor is there because the zero-reading world (0.35) keeps only about
half the gain: the s2 switch already puts a depth-free file in the pair there. That gives about +0.03 to +0.08 Overall,
central +0.05. A misread costs -6 to -8 T3 (D4's
replicate-fold figures), with a chance of about 1%. The upside at our 90th-95th percentile is 0, as the design says.

**(5) Feasibility.** The screening read, the ruling and a rebuild fit before 2026-11-27, given about a day for the build
path (error 4).

**Revisions, if it goes to the captain** (a one-paragraph plain sheet, before anything is read):
1. Present it as two things: reversing the 2026-11-01 (b) reasoning for T3, and reopening the T3 depth line after the
   2026-11-04 exposures. Quote (f)(3).
2. Rewrite the table:
   - (a) The same run, chemistry and pipeline as Gata4: keep s1.
   - (b) An explicit statement that the beta-catenin files used a different MERFISH chemistry or protocol version, or a
     different segmentation or decoding pipeline, than the Gata4 files: recommend s1 without the depth step.
   - (c) Anything else, including "a separate run" with no chemistry or pipeline difference: no change.
   - Drop the "imaged with the WT controls" trigger, which the released text already addresses. Or give it its own row
     (no factor and no depth), which would be a larger ruling.
3. The screening agent applies the table and relays only the branch letter and the triggering sentence, with numbers
   replaced. It withholds stage measurements and identifiers. The ruling records that the switch rests on the Data
   Protocol text alone.
4. Fallback (c)'s organiser question is optional and of low value. The follow-up email of 2026-10-19 (Q5) already asked
   whether the test stacks follow the training acquisition protocol, and no email has had a reply.

## C2: a second channel for the floored-severity answer (KEEP-REVISED)

1. **Repeat.** No. The question is the fifth email's; the channel is new.
2. **Rules.**
   - `docs/ORGANIZER_QUESTIONS.md` names a Slack DM to Weize Xu as a private channel, so the DM is in policy.
   - The optional public post departs from the same file's "not in the public channel" rule. It stays the captain's
     choice, and the default is no.
   - The question applies to every team and reveals no hidden value.
3. **Evidence.** None needed: an organiser answer switches the in-force rule of 2026-11-06 (b).
4. **Expected gain.**
   - An answer pays only if it says zero: 0.35 times the swap's gain under zero (+0.5 to +2.5, N11) is about +0.2 to
     +0.9 T3.
   - The design adds 0.1-0.2 to the chance an answer arrives in time. That is too high:
     - the private Slack note of 2026-10-18 to Weize Xu and Zehua Zeng is still unanswered (third email);
     - the organisers have answered no rules question since 2026-10-16.
   - Take +0.02 to +0.08 instead. That gives +0.01 to +0.05 Overall.
5. **Feasibility.** About fifteen minutes, well before the 2026-11-26 check.

Revision:
- Send the fifth email's question word for word as a DM to Weize Xu. Since it is a metrics question, optionally also
  send it to the organiser the proposal names as lead for tasks, metrics and baselines (Yifan Lu; task-e deep read of
  [2]).
- Ask for a public answer, and post nothing in public unless the captain overrides the ask-privately rule.
- Record any answer with its date in `docs/DATA_POLICY.md`. An answer that arrives after T3's s2 is up changes nothing.
- Correct the misattributed N10 figure (error 1).

## C3: section-artifact audit and a section-robust factor (KILL)

1. **Repeat.** It is not a pruned item. Its audit reads per-section, per-gene detection in both fold truths (Mab21l2,
   Gata4 rep2). That is acceptable as a structure-only read (the N6 and N7 precedent), provided no flag from a file is
   ever used on the fold where that file is the truth. The design says it is not.
2. **Rules.**
   - As a factor change it needs a harness candidate, a TRIES row first and the captain's ruling after a pass. The
     design has all three.
   - The claim that it "applies equally ... so it leaves the Gata4-free hedge intact" is wrong. On the barred file, the
     change would be gated by an audit that reads the Gata4 files, and adopted on the Mab21l2->Gata4rep2 fold, whose
     truth is a Gata4 file.
   - `docs/S2_RULES.md` treats a choice made on Gata4 scores as disqualifying a file from the full-bar case: "its 'no
     factor' was chosen on Gata4 scores".
   - It could therefore enter the allowed s1 only, as task-c's audit said of the related lever 3 (i). The barred file
     named in the 2026-11-06 (b) ruling would have to stay as built.
3. **Evidence.**
   - The shrinkage-matched control is the right design for any factor try on these folds.
   - Its keep rule ("on every fold where it changes the factor") counts a fold where the change runs but does nothing
     as not applicable. The approved reading of 2026-11-01 covers only a pseudo-task "a change cannot run on", so a
     one-fold pass needs a new reading from the captain.
   - The design's evidence from the released files is about the wrong thing. The 3.26x and 2.35x levels are cell-count
     anomalies, each a single labelled section (DATA_NOTES), not detection losses.
4. **Expected gain (rule 5).** A row must state a gain of at least 0.3 on its counting folds, or it is pruned before it
   runs.
   - A section-level dropout biases a gene's mean in the file by about the loss divided by the number of sections, per
     affected section. With 36-46 sections that is about 2-3% per section, and a few percent for a gene hit in two or
     three sections.
   - Ranks over 500 genes barely move, so de_score and de_direction move by hundredths of skill, and the severity gate
     does not depend on scale.
   - K2a's whole-factor gene identity was worth only +0.18 over its permuted control.
   - So the expected gain over the control is well below 0.1, and P(adoption | present) = 0.3 has no support. A pass at
     this power would be the lucky fit TRIES was set up to stop.
5. **Feasibility.** Fits before the freeze, but it would spend three more counted T3 runs.

Corrected gain: about 0 (at most +0.01).

## C4: lineage-distance source factor (KILL)

1. **Repeat.** It repeats the pruned "T3: a background- or spillover-aware factor", as the design says. The new parts
   are a fixed distance in place of a background threshold, and the control. Two of the prune's reasons still stand: it
   needs a new harness candidate and an amendment of rules note 5, which is hash-pinned.
2. **Rules.** It re-measures the source factor, so the same Gata4-free taint as C3 applies, and it would be allowed
   branch only (task-c audit of lever 3 (i)).
3. **Evidence.**
   - *The control is ill-defined here.* A re-estimate on fewer cells moves every gene, in both directions, and on
     average away from identity, because noisier estimates are larger in magnitude. A control that moves random genes
     toward identity does not absorb the folds' known preference.
   - *The spatial selection is confounded.* The far set must be built the same way in the knockout and in its WT, but a
     Mesp1-Cre knockout changes the very lineage the distance is measured from (Gata4 depletes heart states). The factor
     then picks up position-dependent biology, which is gene-structured and is not run processing.
   - *The label weights shift.* The endoderm compartment includes the pharyngeal endoderm next to the heart field, so the
     30% gate is likely to remove much of one of the factor's four compartments.
4. **Expected gain.** The design's own P(adoption) is about 0.1, so the expected gain on the folds is far below 0.3, and
   rule 5 prunes it.
5. **Feasibility.** Two to three DGX days and three more counted runs.

Corrected gain: about 0.

## Strongest candidates

1. **C2 (KEEP-REVISED).** It is nearly free, within policy and positive in expectation, and it feeds a ruling that is
   already in force.
2. **C1 (RULING-ONLY).** Its value is small (about +0.05), but it is the only item aimed at A23, the bet behind 43% of
   s1's dp. Take it to the captain as a short plain sheet with the revised table and the exposure note, not as a
   workflow.

## The biggest risk in the designer's reasoning

The design's headline is "the largest item is information": C1, priced at +0.2 Overall and called compliance-clean.
- **The price** uses P(depth mismatch) = 1 - p_U2 = 0.3 as the chance that the organisers' text will trigger a switch.
  - The released Data Protocol text already says knockouts and their WT controls were processed in separate runs.
  - "A separate run" on the same chemistry changes nothing under the design's own table.
  - An independent knockout run (Mab21l2) was still shallow.
  - Corrected, C1 is worth about a quarter of the design's figure.
- **The legality section** ("the post-exposure block is untouched") misses that T3's depth premise is the facet the
  2026-11-04 exposures bear on, by the team's own reasoning in DATA_POLICY 2026-11-04 (f)(3) (beta-catenin file sizes
  plus the T3 cell limits give bytes per cell).
- Acted on as written, C1 would spend a captain's ruling and a screened re-read of hidden-target text after an exposure
  to change a frozen T3 output. In exchange it offers roughly +0.05 in expectation and risks -6 to -8 T3 on a misread.

## Process and compliance

- **Read in full or in the requested parts:**
  - `design/t3.md`, `brief.md`, and `final/task-c.md`, `task-e.md` and `task-f.md` in full;
  - `final/task-a.md`, `task-a2.md`, `task-b.md` and `task-d.md`, searched for T3 content (none found, as the design
    says);
  - `docs/CHALLENGE.md`, `AGENTS.md`, `experiments/harness/TRIES.md` (rules, headroom, the K and N rows, "Pruned before
    running"), `docs/PLAN.md` (sections 3, 4.4 and 8, and the decision log from 2026-11-04 to 2026-11-06),
    `docs/DATA_POLICY.md` (the held-out rules, windows, Scores, rulings and change notes), `docs/S2_RULES.md` (preamble
    and T3), `ASSUMPTIONS.md`, `docs/ORGANIZER_QUESTIONS.md`, `experiments/research_2026_11/REPORT.md`,
    `results/upload_gate.json` and the estimate file.
- **For verification only:** `docs/DATA_NOTES.md`, `METHOD_SUMMARY.md` (T3 parts), `configs/pseudo_tasks.yaml` (T3
  folds), `experiments/t3_decomp/results/gata4/gata4_summary.csv`, the K2 decision files, and the option lists of
  `build_submission.py` and `diagnostics.py`. All are committed practice results or code.
- **Commands.** One read-only `git log --oneline -5`, run without `GIT_OPTIONAL_LOCKS=0`. Stdlib-only `python3 -I -B`
  one-liners parsed committed CSV and JSON files; no repository code ran.
- **Not opened:** `raw/` and `screened/` of this workflow, `boards/`, `data/`, `runs/`, and anything under
  `~/.claude-account` or `~/.claude`. No web. The only write is this file.
- **Barred material.**
  - Class Q: the harness put the auto-memory index into this agent's context. It was not opened or used.
  - Kept class: the decision-log read displayed screened overall leaderboard totals. They were neither used nor
    copied.
  - No class A, B, W, L or Y item was met. No hidden-target figure was computed, and no bytes-per-cell figure is
    written here.
- **Not done here.** The relayed user request (sync the DGX and this machine, and make sure neither has unpushed
  commits) is a step for the main session. This read-only critique ran no sync, push or DGX command.
