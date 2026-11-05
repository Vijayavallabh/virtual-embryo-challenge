# Cross-board design: adversarial critique

Written 2026-11-06 by the AI assistant (adversarial critic for the cross-board design in the state-of-the-art
workflow). Critique only: no try, no target read, no score computed, no repository edit, no web. Every number is
quoted from a committed or named file, or is hand arithmetic on such numbers (section "Process"). 209 appears only
where the design uses it; nothing here is chosen by it.

## Bottom line

- **The cross-board items are worth about +0.35 Overall, not +0.7 (+0.86 with C4a).** C1 falls from +0.35 to about
  +0.05 once revised (as proposed, its sign is uncertain), C2 from +0.31 to about +0.15 on a prior only the captain
  can hold, C4a from +0.15 to +0.10, and C3 stays at +0.05. C4b, C5 (prior and births) and C6 are killed. On the
  design's own method, P(>= our q95) moves from 5% to about 6% (not 8.4%), P(>= our q90) from 10% to about 11.5%, and
  P(>= 209) from 0.42% to about 0.6%.
- **C1's premise is false.** "No upload in that round depends ... on any pending answer except a Gata4 restriction"
  contradicts `docs/S2_RULES.md`: every s2 upload is where a late A20 bar is absorbed ("before a board's s2 upload"),
  and T3's last upload is where a Gata4 restriction or full bar is absorbed ("before the last T3 upload").
  Front-loading every second-round upload to 27-28 Nov turns a late adverse ruling into "no allowed entry" on the
  boards it touches. The design prices that failure mode at -6 to -10 Overall in R4 and cites it as C3's reason, but
  at -0.03 here. Revised: fixed dates after screened checks (29-30 Nov IST), 1 Dec kept as buffer.
- **C2 survives only as the captain's own prior bet.** Its facet is misstated (the expression-change group does not
  move with the mixture fraction); its recorded reason leans on the closed depth line (REPORT sheet B: "leans yes,
  because of the depth hedge"); its toy assumes E7.75 lies on the line between its brackets, which the committed E7.5
  row contradicts; the proposed pseudo-truth table is circular; and R2's ground for decline (amending the frozen
  embryo s2 after the fact) applies with less evidence than R2 had.
- **C3 is the most robust item**, revised to go private first (`docs/ORGANIZER_QUESTIONS.md`: "Send privately ... not
  in the public channel") and to drop its T3-s1 arithmetic, which rests on the replicate fold that never counts.
- **Operational blocker found:** `experiments/gate/path_to_target.py`, the source of the design's "committed basis",
  is untracked on both machines. The submit client's clean-tree guard refuses every official upload while an
  untracked file sits outside `results/`, `runs/` and `experiments/**/results/` (`docs/P3_RUNBOOK.md`, "Never
  today"; `src/vec/submit.py`). Commit it with a review, or delete it, before 27 Nov.

## Verdicts

| # | Candidate | Verdict | Design's dE | Corrected dE, Overall (range) | Decisive reason |
|---|---|---|---|---|---|
| C1 | Front-load the second round and marks | KEEP-REVISED | +0.35 (+0.1 to +0.8) | +0.05 (0 to +0.15) as revised; about 0 (-0.5 to +0.2) as proposed | Forfeits the frozen ruling switches; gain is an assumed outage reduction |
| C2 | Embryo s2 at mixture fraction 1/4 | RULING-ONLY | +0.31 (-0.2 to +1.5) | +0.15 (-0.3 to +0.8) | Captain's prior only; depth line, frozen rule, circular tool |
| C3 | Organiser push on A24 and an answer table | KEEP-REVISED | +0.05 (0 to +0.3) | +0.05 (0 to +0.2) | Cheap and robust; private route first; drop the s1 claim |
| C4a | A15 read as not covering T1 gene sets; hygiene | RULING-ONLY | +0.15 (-0.2 to +0.5) | +0.10 (-0.3 to +0.4) | No counted check; sign risk on maturation genes |
| C4b | A15 for a cross-modal heart direction | KILL | about 0 | 0 | Repeats the A15 prune; null concordance |
| C5 | Legal T1 composition prior | KILL | 0 (0 to +0.1) | 0 | No decision value; asymmetry over-read; disclosure risk |
| C5-b | Nuclei-anchored births | KILL | -0.4 | 0 (not done) | Pruned births line |
| C6 | XP2 revival | KILL | +0.10 (0 to +0.3) | 0 | Declined 2026-11-05 (c); core objections unanswered |

## C1. Front-load the second round: KEEP-REVISED

**(1) Repeat?** No pruned item, but most of it is already the working rule: marks as soon as a board's last upload
is scored and a recheck on 1 Dec IST (`docs/P3_RUNBOOK.md` section 5); rebuild, `vec check` and server format checks
before the uploads (decision log 2026-11-04 (c), "take the operational risk out of the last days"; the rehearsal of
2026-11-05). New are only the dates and the front-loading.

**(2) Rules.** It keeps the frozen order in letter but disables the frozen switches in practice:
- `docs/S2_RULES.md`: "If a ruling bars the validation answers (assumption A20) before a board's s2 upload, nothing
  built from E10.5, E7.5 or the Gata4 files is uploaded ... That board's s2 goes back to the captain."
- T3: "If they are restricted before the last T3 upload, nothing from the allowed branch is sent, and the last upload
  is the barred branch's s1"; "If such a bar arrives after the refused file is uploaded, the refused file is not
  marked final ... If no uploaded entry is allowed by the rulings received, the captain asks the organisers before
  marking." Decision log 2026-11-06 (b) adds the barred-s2 case, again keyed to "before the last upload".
- So the design's sentence that no second-round upload depends on a pending answer except a Gata4 restriction is
  wrong. With every slot used by 28 Nov, an A20 "no" or a Gata4 full bar arriving between 28 Nov and 2 Dec leaves the
  affected boards with two files built from barred inputs. Today's dates keep the s2 slot, and T3's last upload, as
  the one legal fallback.

**(3) Evidence.** None. The +0.35 is an assumed 1-3 percentage-point fall in miss probability times 21.3. Partial
misses are worth less than 21.3 (T3's s1 13.7, T1's s2 4.7, the heart's s2 2.3, the embryo's s2 0.5 Overall;
`experiments/gate/results/estimate_2026-11-05.json`, readings A and B).

**(4) Realistic gain.**
- Gain side: today's rule already finishes the second round by the end of 30 Nov IST and keeps 1 Dec for outages
  (`docs/PLAN.md` section 8). Two more days of buffer buy perhaps 0.3-1 point of miss probability: +0.05 to +0.2.
- Cost side, on the design's own arrival odds (about 1% for a ruling inside a three-day window; say 1-5% for the A20
  answer near the deadline): an A20 "no" (0.12 by the working answer) that voids uploads made before it (0.3-0.5)
  loses the rebuilt legal entries the s2 slots would otherwise carry, about 150 Overall: -0.05 to -0.45. A T3 full
  bar adds up to about -0.1. Net about 0 (-0.5 to +0.2).
- Unpriced: each upload is published when made (Rules s.11), so front-loading shows rivals our best boards about three
  days earlier.

**(5) Time.** Hours of runbook text.

**Revision.**
1. A dated runbook (`docs/P3_RUNBOOK.md` section 5): 27 Nov 05:30 IST rebuild at the freeze commit, checks and the
   first round; 29 Nov IST T1, heart and embryo s2, each after a screened mail, Slack and Rules check; 30 Nov IST
   (morning) T3's s1 after its screened check; 1 Dec IST buffer and marks recheck. Reconcile `docs/PLAN.md` section 7
   ("the s2 round by 2026-12-01 06:00 UTC") with section 8 and S2_RULES ("by the end of 30 Nov IST"); the design
   quotes the later date as today's plan without noting the conflict.
2. A pre-filled go-ahead sheet per upload: board key, file, sha256, build commit, disclosure ids.
3. Marks as soon as both of a board's uploads are scored (already the rule).
4. A clean tree on both machines before 27 Nov: commit `experiments/gate/path_to_target.py` with a review, or delete
   it (see Bottom line).
5. Front-load only if the organisers confirm A20 and the Gata4 files before 27 Nov; the option it forfeits is then
   worth nothing.

## C2. Embryo s2 at mixture fraction 1/4: RULING-ONLY

**(1) Repeat?** Not new in kind. It is the pruned clock warp moved into the s2 slot (TRIES, "Pruned before running":
"Nothing measured on released stages says which way E7.75 would sit"), on R2's procedure (amend the frozen embryo s2
after the fact; declined 2026-11-05 (c)) with less evidence than R2 had (R2: +0.154 board points on selected-on
backtests; C2: none). The design lists R2 but not the clock-warp prune, "Embryo: a depth-shaped draw of real bracket
cells" (pruned: "The depth line is closed") or TRIES item 3, "Embryo: stop polishing".

**(2) Rules.** Windows, genotype and A15 are not engaged, and count, size rule and scope are unchanged. But:
- it amends a frozen rule after the fact;
- its recorded reason is partly depth. A 1/4 mixture moves the s2's average depth from about 115 to about 84 counts
  per cell (geometric means of E7.5's 62 and E8.0's 212 at 1/2 and 1/4; TRIES, "Diagnostics behind these choices"),
  on a board whose depth line the ledger closed (E7, H3; `docs/DATA_POLICY.md` 2026-11-04 (f)(3)). The ruling must
  reopen that line explicitly or drop the depth reason, and REPORT sheet B's deciding prior ("near E7.5 or E7.5-like in
  depth") then has to be restated;
- provenance: the "E7.5-like depth" risk was logged in decision log 2026-11-04 (c), the day of the board-index
  exposure. Its basis must be the released-stage depth diagnostic of 2026-11-03, and the sheet should cite that alone.

**(3) Counted evidence.** None, as the design says. Two errors in what it cites:
- Expression change (about 8.3 left) does not move with the fraction. de_score and de_direction rank dp, and a
  real-cell mixture's dp is t times the bracket difference plus the shared reference noise; going from 1/2 to 1/4
  rescales the signal, which ranks ignore, and changes only its ratio to the shared noise. Only cell state and
  neighbourhoods are mixture-sensitive (REPORT section 6: "half the board", not "the expression and cell-state half").
- The cell-state gap it cites (11.50) "comes from E7.25" (TRIES headroom notes), rule 2's known feature.

**(4) Realistic gain.**
- The toy skill(t; t') = t'^2/(t'^2 + (t - t')^2) assumes E7.75 lies on the line between its brackets. The committed
  E7.5 row says real stages do not: at about midway, the best mixture still leaves 0.0380 of mmd_u against 0.083-0.086
  for a bracket copy (TRIES diagnostic table), far above the split-half ceiling (fgw's committed mmd_u skill of 0.671
  on that practice task puts it near 0.002). With that residual the floor distance does not vanish near E7.5, so the
  winning-world gain shrinks. What holds is the core: if E7.75 is E7.5-like, s1 sits below the floor on the
  distributional groups and a 1/4 s2 recovers part of it.
- The loss branch (-0.2 Overall) disagrees with REPORT sheet B ("loss if wrong: about 1.5 embryo points, 0.75
  Overall") without saying why. In expectation the loss is bounded by today's s2's hedge value (+0.92 embryo, 0.46
  Overall; `docs/S2_RULES.md`).
- Corrected: with the captain's own P(E7.75 E7.5-like) at 0.2-0.3, about +1 to +1.5 Overall there and -0.1 to -0.4
  elsewhere (bounded by the 0.46 above): +0.15 (-0.3 to +0.8); about +0.4-0.5 points of chance at q95 and +0.8 at q90,
  not +1.1 and +1.5. Break-even stays near 0.15-0.2 (REPORT).
- "The only s2 that does not hedge its board's main uncertainty" is wrong: the embryo s2 hedges the board's one
  measured test risk, the flattened s1 geometry (TRIES item 3; `METHOD_SUMMARY.md`).

**(5) Time.** A flag and a frozen-file amendment; the deadline is the embryo's second upload.

**The decision tool is circular.** Pseudo-truths built from the test pair at t' = 0.1-0.9 are mixtures of the same
bracket cells the predictions draw, so they have no off-line residual by construction. The table would show the
closed form's largest asymmetry and cannot inform the one input that decides the case. Skip the DGX day; give the
captain the closed form with a residual range as a stated assumption (not calibrated on E7.25 or E7.5, which cannot
count for a mixing-fraction change, `docs/DATA_POLICY.md` 2026-11-05 (b)(v)).

**Price for the sheet.** +0.15 if he holds the prior; the sheet states plainly: no evidence of direction, a depth
change on a closed line, a frozen rule amended after the fact (R2's ground), cell state and neighbourhoods only.

## C3. Organiser push and answer table: KEEP-REVISED

**(1) Repeat?** It repeats the fifth email; new are the channel and the bundle. No pruned item.

**(2) Rules.** The standing default is private (`docs/ORGANIZER_QUESTIONS.md`: "Send privately ... not in the public
channel: some questions reveal method choices"), and P3 methods stay out of Slack until 2026-12-11
(`docs/DATA_POLICY.md`). The severity question is generic and safe. The A20 re-ask must not say how we use the files
(the fourth email named E10.5 as the starting stage and the Gata4 replicates in the beta-catenin prediction). The
runbook already allows a public question about the Final column by 28 Nov IST. A public "zero" answer also informs
rivals who have not uploaded, while our s1 is frozen and our s2 switch is pre-wired, so its rank value is below the
design's framing ("information, not a loss").

**(3) Evidence.** "T3's s1 needs no answer-contingent switch" is argued from stand-in A's replicate-fold scores
(74.50/72.15 against 68.05/64.09), a fold that never counts, and s1 has no switch to make without a new ruling. Drop
it.

**(4) Gain.** +0.05 holds: a small lift in P(answer by 26 Nov) x 0.35 (zero) x about +1.2 T3.

**(5) Time.** Hours.

**Revision.** The severity question first by the Slack DM route `docs/ORGANIZER_QUESTIONS.md` names, a generic public
post only as the captain's explicit choice; the answer table built only from switches that S2_RULES and decision log
2026-11-06 (b) already define (it adds no rule). Optional, so that the A20 "no" row can actually be carried out:
pre-build and format-check A20-no fallback files (T1 and heart from E9.5, embryo from E7.25 and E8.0; the T3 barred
branch is built) under rules the captain signs.

## C4a. A15 reading and Feng direction hygiene: RULING-ONLY

**(1) Repeat?** Not pruned (REPORT sheet D is open). Neighbours: the pruned noise-aware Feng direction (another axis)
and task b's L5 consensus filter (flagged by its audit).

**(2) Rules.** A15 covers released data across tasks and new cross-modal input; a T1-only external gene set is
neither, so the captain's reading is the right route. Conditions: an `EXTERNAL_DATA.md` row per gene-set source with
`docs/PLAN.md` 4.4-style provenance; gene classes fixed in a ruling row before any build; no class chosen by reading
E10.5's DE genes (a counted truth); T1's s2 back to the captain (`docs/DATA_POLICY.md` 2026-11-04 (d)); R1 untouched.

**(3) Evidence.** None counts: Feng's E9.5 lies inside every T1 window (TRIES item 5). The support is cross-tissue
(human brain and immune cells, lung, iPSC cardiomyocytes).

**(4) Gain.** Real headroom (all of M1's gain is direction, +3.01 in-sample), but sign risk: cardiomyocyte maturation
after E10.5 includes genuine mitochondrial, oxidative and ribosomal-programme change, so holding detectability or
organelle modules flat can drop true DE genes from de_score's lists. A target-free size check first (the share of
Feng's top 400 up and down changes per type that fall in the modules; N2's precedent). +0.10 (-0.3 to +0.4), about
+0.3 points of chance at q95.

**(5) Time.** 1-2 days; ruled by about 24 Nov to rebuild and re-hash before the 27 Nov freeze.

## C4b, C5, C5-births, C6: KILL

- **C4b.** It repeats the A15 prune ("Heart: Feng's direction ... on the heart panel"): closed by A15's working
  answer, a null target-free concordance (cardiomyocyte Spearman 0.05-0.10), and no counted check. The design itself
  says nothing is new.
- **C5.**
  - No Overall gain by the design's own account. The work only moves an estimate that cannot cross 209 and gates
    nothing after 27 Nov.
  - The asymmetry is over-read. Under-shooting costs 11.1 against 2.0 only at the grid's 5%-blood point, a
    near-removal of the blood cluster. At its 40%-blood point, about as far above E10.5's own share as 5% is below it,
    the costs are about equal: own blood against the 40% truth 70.13 vs 74.65 matched (-4.5); 40% against the own
    truth 66.10 vs 69.92 (-3.8) (`experiments/s2_hedge/results/t1_grid/payoff.txt`). `docs/PLAN.md` 4.1's pseudo-task put over-prediction
    at 7-9 points per 10. The grid has no point near a legal prior's centre (about 12-17%), so "no legal prior that
    lowers E12.5's blood share should move s1" is not shown; the captain's ruling of 2026-11-01 is what keeps s1
    there.
  - Leaning on Leonard's pooled 29% (mutant libraries included) or the erythroid fractions for any s1 or s2 choice
    before register rows would recreate the disclosure gap of `docs/DATA_POLICY.md` 2026-11-05 (c).
- **C5-births.** The pruned births line (-0.2 to -0.4; the design rejects it too).
- **C6.**
  - It repeats the late-fetal family the captain declined on 2026-11-05 (c) and answers neither core objection
    (whole cells against nuclei; late-fetal against E10.5-E12.5 biology).
  - Approving a new family as an adoption-rule reading widens what counts under the adoption rule (the captain's
    call, under the campaign audit's conditions) for +0.1 at best, and every lever it could count is itself
    ruling-only.
  - 4-6 days plus register rows, an X003 amendment, the Q7 label check and an independent leakage review are unlikely
    to be finished and checked by 27 Nov. The design recommends declining.

## The rejected list (R1-R7)

Agreed, all seven. Two notes: R4's price for losing the Gata4-free insurance (-6 to -10 Overall) is the failure mode
C1's front-loading opens in the late window, so the two are priced inconsistently; R2 is correctly barred (a score
read as which world we are in computes).

## Corrected ranking and path

- By expected gain: C2 (+0.15, ruling) > C4a (+0.10, ruling) > C1 revised = C3 revised (+0.05) > the rest (0).
- By robustness (no ruling and no unmeasured prior): C3 revised, C1 revised.
- Cross-board total if every ruling goes as recommended: about +0.35 Overall. With R1 and sheet C as the design counts
  them (sheet C not checked here): about +0.8, so P(>= our q95) about 7.5%, not 8.4%. Under task h's T3 payoffs
  nothing changes: 209 stays out of reach.
- Presentation: take the P(>= 209) column out of the options table and report it once, apart. 209 is the leading
  test-board total (`docs/PLAN.md` section 8) and must not look like a column the options are ranked on.

## Strongest candidates

1. **C3, revised:** the severity question sent now (private route first) and one answer table built from existing
   switches. It is cheap and positive, and the only item that shrinks a disqualification-sized risk (A20).
2. **C1, revised:** fixed dates after screened checks with 1 Dec as buffer, pre-filled go-ahead sheets, marks as soon
   as scored, and the clean-tree fix. Operational de-risking that keeps every legal switch alive.
3. **C2, as the captain's own prior bet** with its rationale cleaned (no depth reason, cell state and neighbourhoods
   only, no pseudo-truth table): the largest corrected mean (+0.15), but only if he holds P(E7.75 E7.5-like) above
   about 0.2. C4a (+0.10) is close behind and procedurally cleaner.

## Biggest risk in the designer's reasoning

C1 is ranked first and needs no ruling, so it is the item most likely to be acted on unreviewed, and it rests on a
misreading of the frozen rules. Under `docs/S2_RULES.md`, the s2 uploads absorb a late A20 bar and T3's last upload
absorbs a Gata4 restriction or full bar. Front-loading them trades an assumed 1-3 point fall in outage risk for the
only legal fallback in late-ruling worlds, where the outcome is "no allowed entry", the disqualification-sized loss
the design itself prices at -6 to -10 Overall in R4 and gives as C3's reason. By the design's own arrival odds the
forfeited option is worth -0.05 to -0.5, so the headline +0.35 has an uncertain sign.

## Process and compliance

- **Read:**
  - the design file;
  - the brief;
  - the audited notes `final/task-a.md`, `task-a2.md`, `task-b.md` to `task-f.md`.
- **In the repository:**
  - `docs/CHALLENGE.md`, `docs/S2_RULES.md`, `ASSUMPTIONS.md` and `experiments/harness/TRIES.md` (all);
  - `docs/PLAN.md`: sections 1-4.4, 7 and 8, and the decision log of 2026-11-04 to 2026-11-06 (b);
  - `docs/DATA_POLICY.md`: lines 1-95 and 239-646 (windows, Scores, rulings, change notes);
  - `experiments/research_2026_11/REPORT.md`: sections 1-11, to the references;
  - `docs/ORGANIZER_QUESTIONS.md`: the header, the working answers, and the fourth and fifth emails;
  - `docs/P3_RUNBOOK.md`: section 5 and "Never today";
  - `METHOD_SUMMARY.md`: its headings and the passages a grep for "embryo" and "hedge" printed;
  - `results/upload_gate.json`;
  - `experiments/gate/`: `results/estimate_2026-11-05.json`, `results/options_tail_2026-11-05.txt`, and the untracked
    `path_to_target.py` with its result JSON;
  - `experiments/s2_hedge/results/t1_grid/payoff.txt`, lines 1-19 only;
  - the clean-tree guard lines of `src/vec/submit.py` (grep).
- **One display to record:** the grep that located `payoff.txt`'s footer printed the first 60 characters of line 20:
  a practice truth's mix and the start of its size. It is E10.5 validation material, not a test board, and is not
  used here.
- **Never opened:** `raw/`, `screened/`, `boards/`, `data/`, `runs/`, `docs/research/`, the leaderboard, any memory
  file, session log or harness-saved output. The harness put the auto-memory index into context; it was not opened or
  used.
- **No barred class met.** I derived no embryo target's best mixing fraction and quote none (`docs/DATA_POLICY.md`
  2026-11-05 (b)(v)). The E7.5 residual point uses the committed table and TRIES's own "about midway" only.
- **Computation:** hand arithmetic on committed numbers only. It includes one inversion of the skill map, on a
  committed practice result of released files (fgw's mmd_u skill on the E7.5 backtest, against the committed
  diagnostic table), to size a split-half ceiling. `docs/DATA_POLICY.md` (Scores) bars inversion for hidden targets
  only.
- **The relayed user request** (sync the DGX and this machine, and leave no commits unpushed) is the main session's
  step, so I did not fetch, pull, push or sync. A read-only check with `GIT_OPTIONAL_LOCKS=0` and no fetch:
  - local main matches its last-fetched origin/main at 9ce2742, with 0 commits ahead and 0 behind;
  - untracked: `experiments/gate/path_to_target.py` and its result JSON. The `.py` trips the upload guard (Bottom
    line);
  - the DGX was not checked here.
