# Cross-board design: evidence routes, rulings, organiser questions, upload strategy, path to 209

> Note 2026-11-06 (main session, after the counter-review): the passages of this file that measure candidates against the captain's 209 are struck (marked below). 209 gates uploads only and is never a threshold or evidence for a try, an adoption or a method choice (docs/DATA_POLICY.md 2026-11-05; docs/PLAN.md section 8).

Written 2026-11-06 by the AI assistant (cross-board designer in the state-of-the-art workflow). Design only: no try,
no target read, no score computed, no repository edit. Every number is either from a committed result (named) or
arithmetic on committed results (section 0 says how). 209 is the captain's target and upload gate. It is reported for
information only and is never used to choose anything here. Options are ranked by expected Overall gain and by the
chance of reaching our own 90th and 95th percentiles.

## Bottom line

- **No ruling or organiser answer on the table opens a lever worth more than about +0.3 Overall.** The three the
  captain asked about come out like this:
  - A captain-written T1 composition prior from legal data mostly confirms today's T1 pair. On the committed
    oracle-base grid, under-predicting blood costs 11 points and over-predicting costs 2. So no legal prior that
    lowers E12.5's blood share should move s1. The 60% cardiomyocyte s2 keeps its value unless the prior puts less
    than about 0.12 on a cardiomyocyte-rich E12.5.
  - The late-fetal counting family (XP2) is now buildable but can count only T1 expression post-processing:
    about +0.1 expected.
  - A15 for a cross-modal heart direction is worth about 0. The target-free concordance is null, and Feng's E9.5
    sits inside the heart backtests' windows, so no counted check exists.
- **The largest cross-board item needs no ruling: front-load the second round.** The second uploads are worth
  +21.3 Overall in expectation (reading A 191.5 against the first round's 170.2). T3's s1 alone is +13.7 of that
  and all of our upper tail, yet today's plan uploads it last, by 30 Nov IST. If every second-round upload goes in on
  27-28 November (the frozen order allows it), with the marks set at once and 29 Nov-1 Dec kept as buffer, the
  expected gain is about +0.35 Overall (+0.1 to +0.8).
- **The one open s2 decision with a positive tail value is the embryo backup at mixture fraction 1/4** (2026-11-05
  report, sheet B). It is the only s2 that does not hedge its board's main uncertainty. It is worth about +0.3
  Overall, +1.1 points of chance at our 95th percentile and +1.5 at our 90th. It is a captain's ruling, because it
  amends a frozen rule after the fact.
- **Path to 209.** [struck 2026-11-06: information only; 209 chooses nothing] Today's estimate (reading A) gives P(total >= 209) = 0.42%.
  - The draws that reach 209 have every board near its own 90th percentile at once: T3 about 75 (s1's Gata4 factor
    and depth both carrying), T1 about 73, heart about 60, embryo about 63.
  - The cross-board items plus the skewed T1 bet on record lift P(>= 209) to about 1.5%, P(>= our q95) from 5% to
    about 8.4%, and P(>= our q90) from 10% to about 14%.
  - The rest has to come from the entries' existing bets landing.
  - Under task h's measurement-based T3 payoffs (s1's best case about 66), 209 would need T1 and the T2 mean far above
    their 90th percentiles, so it is essentially out of reach.

## 0. How the options are priced

Committed basis: `experiments/gate/results/path_to_target_2026-11-06.json` and `options_tail_2026-11-05.txt`
(reading A: each board at the better of its two uploads, with K2a; 200,000 draws, boards independent).

| Level | Overall | P(Overall >= level) today |
|---|---|---|
| mean | 191.50 | - |
| our q50 / q75 | 192.89 / 200.57 | 50% / 25% |
| **our q90** | **203.95** | **10%** |
| **our q95** | **205.54** | **5%** |
| our q99 | 208.07 | 1% |
| the captain's 209 (information only) | 209.0 | 0.42% |

- Board 10th-90th bands (reading A): T1 64.8-73.0, heart 48.5-60.1, embryo 56.9-63.5, T3 50.4-74.6.
- Board means in the draws that reach our q95: T1 72.4, heart 58.5, embryo 62.3, T3 74.3. In the draws that reach
  209: T1 73.0, heart 60.4, embryo 63.0, T3 75.1.
- A sure +d on a weight-1 board gives P(>= 209) of 1.07% at d = 1, 2.2% at 2, 3.9% at 3, 9.8% at 5 and 23% at 8.
  It takes +16.1 for an even chance (`whatif_p_target_after_sure_shift`).
- **Method.** Each option is a sure shift, or a two-point bet (gain g with probability p, loss l otherwise),
  independent of the board outcomes. Its effect on P(>= c) comes from the committed survival curve above,
  interpolated linearly between the committed points. That is accurate to about +-0.3 points of chance.
- **Caveats.**
  - Reading A is an upper reading: T1 credits the oracle-base grid in full, and T3 puts half its weight on the
    in-sample replicate stand-in.
  - Under task h's measurement-based T3 payoffs (s1's best case about 66), every upper-tail probability here is
    smaller. The ranking of the options does not change.

Three adjustments to the estimate itself, if they materialise, move every probability down. A sure -1 Overall halves
P(>= q95) to about 3%:
- Per-set rather than per-board final marks: the s2s gain +5.54 instead of +7.63 (`docs/S2_RULES.md`, bundle_sim,
  even weights, T3 left out). That is about -2.1 Overall at a 0.45 chance.
- A floored severity scoring zero: in reading A's K2a world both T3 uploads floor, about -2.8 T3 at a 0.35 chance.
  Part of this is recovered only if the answer arrives and the ruled swap is made.
- A legal blood prior (C5): -0.9 to +0.4 on T1, depending on its mass on low blood.

## 1. Candidates

### C1. Front-load the second round and the final marks (no ruling needed)
- **Facet:** none directly. It protects the expected value of every facet that the second uploads carry: T3's
  per-gene response (s1), T1's proportions (the cardiomyocyte s2), the heart's scope and the embryo's geometry.
- **Groups and headroom:** all. Most of the value is T3's s1. On stand-in A it stands +23.3 board points over the
  floor, with response identity 11.89, direction 8.05, magnitude 1.03 and mutant state 5.70 left (TRIES headroom).
- **Mechanism:** the second-round uploads (T1 s2, heart s2, embryo s2 and T3's s1) are worth reading A minus reading
  B: +21.3 Overall in expectation.
  - By board: T1 +4.74, heart +2.33, embryo +0.46, T3 +13.74.
  - T3's s1 also carries the whole upper tail. In the draws that reach our q95, T3 averages 74.3, which only s1 can
    give.
  - Today's plan sends T3's s1 last, "by the end of 30 Nov IST", and the s2 round "by 2026-12-01 06:00 UTC".
  - No upload in that round depends on a score (the s2 rules read none) or on any pending answer except a Gata4
    restriction.
  - So: upload the first round on 27 November (T3's s2 per the ruling of 2026-11-06 (b)), and T3's s1 as soon as
    T3's s2 is scored that day or the next. Then upload T1, heart and embryo s2 on 28 November. Set each final mark
    as soon as both of a board's uploads are scored.
  - That leaves 29 November-1 December for outages, rejected files (a rejected file uses no scored attempt, Rules
    s.11), client guards (a stored format check blocks every POST until resolved) and the captain's availability.
- **SOTA basis:** the post-mortems show the top of a leaderboard is noise-sensitive. Only each task's top rank was
  bootstrap-robust in 2022 [task-f 13]. Four re-runs moved the 2021 winner by up to 8 ranks [task-f 14]. Top
  margins were 0.0013 composite on VCC TEST100 [task-f 9]. So losing an upload that carries the upper tail is the
  cheapest loss to avoid. The rest is our own committed arithmetic.
- **What is new:** a refinement of `docs/PLAN.md` section 8 ("never in the last hours") and the fallback dates of
  decision log 2026-11-04 (b).
  - New: pricing the second round as 21.3 Overall, two thirds of it T3's s1.
  - New: moving that s1 from "last, by 30 Nov" to the first or second day.
  - Nearest pruned item: none. Nearest standing rule: `docs/S2_RULES.md` T3 order ("s2 goes first, and s1 last,
    by the end of 30 Nov IST"), which this satisfies.
- **Legality:** no data, no window, no change to any file. The frozen order is kept. Each upload still needs the
  captain's go-ahead at the time.
  - Cost: if the organisers restrict the Gata4 files between 28 and 30 November, after s1 is up, the last-upload
    switch to the barred s1 is lost. Chance about 1% (no organiser reply since 2026-10-16). Loss about 2-4 T3 points
    (the barred s1 against the refused s1). That is about -0.03 in expectation.
- **Counted evidence:** not applicable (scheduling).
- **Expected gain:** +0.35 Overall (range +0.1 to +0.8).
  - Basis: the chance that some second-round upload misses the window falls by an assumed 1-3 percentage points,
    times 21.3.
  - Tail: +0.1 point of chance at q95 and +0.2 at q90. Small in absolute terms, because q95 is a 5% event, but it
    protects the only route to it.
- **Effort and risk:** about an hour of runbook text; no DGX. Risk low.
- **First step:** add a dated schedule and a pre-filled go-ahead sheet to `docs/P3_RUNBOOK.md` section 5 (main
  session):
  - per upload: board key, file, sha256 and build commit;
  - 27 November 05:30 IST: rebuild and checks;
  - then the first round, T3's s1 the same day or the next, T1, heart and embryo s2 on 28 November, marks at once;
  - 1 December IST: recheck the marks.

### C2. The s2 portfolio: keep the frozen pairs, and rule on the embryo backup at mixture fraction 1/4 before the freeze
- **Facet:**
  - embryo: proportions (the share of E7.5 and E8.0 cells) and expression;
  - the other boards: unchanged.
- **Groups and headroom:** embryo cell state 11.50, expression change about 8.3 and neighbourhoods about 9.5 board
  points left on the mean of the two backtests (Overall at half weight).
- **Mechanism.** An audit of every board's second slot, by E[max] and by its value in our upper tail:
  - T1's s2 (cardiomyocytes at 60%) hedges the board's largest uncertainty. Under any stated prior that lets the
    share move, it adds +4.74 to +7.90 to the better-of-two (S2_RULES), and C5 finds no legal prior that flips it.
  - The heart's s2 (the half trim) hedges scope, which the post-exposure block closes to any other choice.
  - T3's s2 is rule-risk insurance (Gata4-free), switched by the floored-severity answer (ruling 2026-11-06 (b)).
    Its upper-tail value is about 0 by design (N10: +0.01 / +0.00 points of chance at lambda 0). It must stay
    Gata4-free; see "Rejected", R4.
  - **The embryo's s2 shares s1's mixture and differs only in geometry.** So the case where E7.75 sits near E7.5
    (the skill map punishes a mid-way mixture most there) is unhedged on the board's expression and cell-state half.
    So is the logged risk that E7.75 was imaged at E7.5-like depth (decision log 2026-11-04 (c)).
  - Sheet B (REPORT section 6) changes only the embryo s2's mixture fraction to 1/4, keeping its geometry, size rule
    and count.
- **SOTA basis:**
  - On final picks the record shows score-chasing fails, not that boldness pays. Magic multipliers tuned to the
    public board ranked poorly, while robust-CV teams rose [task-f 24]. A public 2nd place fell to 133rd
    [task-f 13]. A decorrelated second pick was kept against a public-kernel monoculture (weak source) [task-f 25].
  - The E[max] arithmetic: a second slot is worth most when it is decorrelated along the largest uncertainty
    [task-f mapping 7].
  - For a trailing team, variance pays near its own upper percentiles (REPORT section 1, [66][67]).
- **What is new:** the audit applies the post-mortem evidence and finds that the embryo slot is the one exception.
  - Nearest pruned or declined items: R2 (embryo s2 with aligned moves, declined 2026-11-05 (c)), "Embryo: s1 =
    fgw_aligned re-read on a test-like task", and task f lever 1 (upload-pair audit; overlaps S2_RULES' bundle_sim).
  - Sheet B itself is open in REPORT section 7. This candidate adds its timing (rule it before the freeze, so the s2
    can go in C1's front-loaded round) and its tail pricing.
- **Legality:** no new data. It amends the frozen embryo s2 rule after the fact, which only the captain's ruling can
  do. The count, size rule and scope are unchanged (post-exposure block). It needs no window or A15 question.
- **Counted evidence: no counted check exists.** E7.25 and E7.5 cannot count for a change to the mixing fraction
  (`docs/DATA_POLICY.md` 2026-11-05 (b)(v); TRIES rule 2: both backtests prefer equal mixing). So it is priced as a
  captain's bet. Its target-free decision tool is REPORT's pseudo-truth table from the test pair at t' = 0.1-0.9.
- **Expected gain:** +0.6 embryo (+0.31 Overall), range -0.2 to +1.5 Overall.
  - Two-point model: +1.5 Overall with probability 0.3, -0.2 otherwise. The idealised toy's +1.2 embryo is
    discounted by about a third.
  - Tail: +1.1 points of chance at q95 and +1.5 at q90; P(>= 209) 0.42% -> 0.77%.
  - If the embryo designer prices it differently, use that price.
- **Effort and risk:** about one DGX day for the pseudo-truth table; a one-flag build that is already format-checkable.
  Risk: if E7.75 sits at or past mid-way, it dilutes today's geometry edge (about -0.2 Overall).
- **First step:** the target-free pseudo-truth table (registered as a diagnostic, not a try), then a plain-language
  ruling sheet for the captain before 2026-11-27.

### C3. Organiser answers: a public push on the floored-severity reading, and one answer-contingent table for every upload
- **Facet:**
  - T3 expression (the severity group, through T3's s2 choice);
  - all facets for the A20 and marks questions.
- **Groups and headroom:** T3 response magnitude, with 10.92 board points left on s2's stand-in and 1.03-2.68 on
  s1's. Reading A's T3 falls about 2.8 under the zero reading (section 0).
- **Mechanism.** The only open organiser answer that moves an upload is A24, the floored-severity reading: zero gives
  the barred s1 as T3's s2, and half or no answer keeps today's s2 (ruling 2026-11-06 (b)). Five private emails have
  had no reply, and the organisers have answered no team's rule question on Slack since 2026-10-16. The fifth email
  itself invited a public answer. So:
  - (a) the captain posts the same one-line scoring question in #virtual-embryo-challenge. It applies to every team
    and reveals no method.
  - (b) the same push re-asks three standing questions:
    - what counts when no final mark is set, and whether marks are per board (4th email Q2; working answer per
      board, 0.55);
    - the A20 confirmation (validation files as model inputs; a "no" would force rebuilds, so learning it before the
      freeze matters);
    - each test board's reference stage (A3).
  - (c) one answer-contingent table in the runbook (main session), fixing in advance what each answer does to each
    upload: zero / half / none x A20 yes / no x marks per board / per set. A late answer is then acted on
    mechanically, never by reading a score.
  - T3's s1 needs no answer-contingent switch, by my rough arithmetic. Under the zero reading, its depth step still
    pays in expectation: about +7 T3 in the matching-depth world (p_U2 0.7; replicate stand-in 74.50 / 72.15 with
    depth against 68.05 / 64.09 without) against about -3.4 where it floors (0.3 x 0.9 x 12.5). It also carries the
    whole T3 upper tail. The T3 designer should confirm this with N10's machinery.
- **SOTA basis:**
  - The proposal maps an undefined metric to skill 0, while the live scorer gives 0.5 [task-e 2].
  - The Timeline warns that "exact metric formulations may still be corrected" [task-e 11].
  - The organisers name "ignoring the stage difference between the training knockout and the evaluation ones" as a
    pitfall [task-e 5], which bears on the barred s1. Its carry to Gata4 is mostly the WT E9.5-against-E8.75 offset
    (t3_swap_review). Under the zero reading it is still the better s2, because today's s2 floors in 14 of 20 seeds
    (REPORT section 3).
- **What is new:** the fifth email exists, and ruling (b) covers s2. New:
  - the public channel, to raise the chance of an answer;
  - folding the marks, A20 and A3 questions into one push;
  - the explicit finding that s1 needs no contingency.
  - Nearest item: REPORT sheet A (ruled). Nothing pruned is repeated.
- **Legality:** asking is allowed and recommended (`docs/DATA_POLICY.md`: "an answer given in advance protects us").
  The questions concern scoring and rules, not a property of a hidden target. No method is revealed (P3 privacy
  holds). T3's rules move only by the organisers' answer through the ruled switches.
- **Counted evidence:** not applicable (rulings and answers). Pricing comes from TRIES N10 and N11.
- **Expected gain:** +0.05 Overall (range 0 to +0.3).
  - Derivation: the push raises P(an answer before 2026-11-26) by an assumed 0.05-0.15. Then x 0.35 (zero) x about
    +1.2 T3 (the swap under the zero reading, N11 range +0.5 to +2.5).
  - Tail: about 0 at lambda 0, because the swap helps in worlds where s1 has failed. It is up to +2 points of chance
    at lambda 1 (N10).
  - The A20 and marks answers carry no mean gain, but they remove a disqualification-sized tail risk and a possible
    -2 bias in our estimate.
- **Effort and risk:** 2-3 hours; no DGX. Risk low. A "zero" answer lowers our own estimate, which is information,
  not a loss.
- **First step:** draft for the captain one public Slack post (the floored-severity question, verbatim from the fifth
  email) and one short email with the three re-asked points. Then the answer table in `docs/P3_RUNBOOK.md` section 5.

### C4. A15, read two ways: the captain's reading on external gene sets for T1, and the organisers' answer for a cross-modal heart direction
- **Facet:** expression. (a) T1, the per-type direction; (b) the heart, per-lineage change.
- **Groups and headroom (Overall units):**
  - (a) T1 direction 9.49 and DE recovery 13.00;
  - (b) heart expression change 6.27, the heart's near-tied weakest group (TRIES).
- **Mechanism.**
  - **(a)** A15 asks whether released data may cross tasks. External, provenance-checked gene sets are not released
    data, yet the DGX campaign prune and REPORT sheet D hold three T1 bets behind "A15":
    - Feng direction hygiene: hold flat the dissociation-stress, ambient and nuclear-detectability modules;
    - a maturation quota on cycling cells;
    - a blood direction from legal globin anchors [task-b L3].
    The captain can rule that A15 does not cover them. Each then still needs an `EXTERNAL_DATA.md` row, `docs/PLAN.md`
    4.4-style provenance and a TRIES ruling row with the gene classes fixed before any run (R1's precedent). Hygiene
    is the cleanest. Task b's protocol-and-lab axis [task-b 17, 18, 23] and task a's cross-study agreement result
    (single-study error 8.03 against 2.65 weeks pooled [task-a 6]; carry-forward unbeatable for low-agreement genes
    [task-a 5]) support holding protocol-sensitive genes flat.
  - **(b)** An organiser "yes" on A15 would let Feng's scRNA direction enter the heart panel. But the target-free
    concordance is null (cardiomyocyte Spearman 0.05-0.10, TRIES prune "Heart: Feng's direction ... on the heart
    panel"). Per-type momentum carries over only for cardiomyocytes. By my reading, Feng's E9.5 would also lie inside
    the heart backtests' own external-data windows, as it does in T1's (8.5, 11.5]; the harness defines those
    windows, so this needs checking.
- **What is new:**
  - (a): the protocol-and-lab evidence [task-b 17, 18], the cross-study error sizes [task-a 5, 6], and the point that
    A15's own wording does not reach external gene sets. Nearest: REPORT sheet D (open, not pruned), "T1: a
    noise-aware, subtype-resolved Feng direction" (pruned; different axis).
  - (b): nothing; it repeats the A15 prune.
- **Legality:**
  - (a) Windows: gene sets carry no stage data, but any set built from in-window measurements stays out. The
    globin direction's start cannot be the screened-out early anchor (it must lie at or before E7.25, or at
    E8.0-E8.25). It amends R1 for blood, and a changed s1 sends T1's s2 back to the captain (DATA_POLICY 2026-11-04
    (d)).
  - (b) needs the organisers' A15 answer (working answer 0.75, never a ruling). A shape input from imaging is closed
    by A15's plain words [task-d lever 2].
- **Counted evidence:** no counted check exists for either. For (a), no T1 backtest can judge a change to Feng's
  direction, and blood is diluted about tenfold on the counted tasks. XP2 (C6) would be its only counting route. So
  each (a) item is priced as a captain's bet.
- **Expected gain:**
  - (a), one bet taken (hygiene): +0.15 Overall (range -0.2 to +0.5). Tail: +0.55 points of chance at q95, +0.75 at
    q90; P(>= 209) 0.52%.
  - (b): about 0 (+-0.25 Overall). An answer before the freeze is unlikely.
- **Effort and risk:**
  - (a) 1-2 days, DGX CPU (a target-free size diagnostic like N2, then build). Risks: holding genes flat can drop
    true DE genes from de_score's top lists; a wrong-sign maturation quota costs -0.3 to -1.5.
  - (b) 0 effort; the risk is spending a try on a null concordance.
- **First step:** (a) the captain's reading of A15's scope, then a ruling row for hygiene with its gene classes fixed
  in the row. (b) none: wait for the organisers.

### C5. A captain-written legal T1 composition prior from the new nuclei and erythroid anchors (the legal route)
- **Facet:** proportions: the blood share of chamber cells, the cardiomyocyte share of heart cells, and (sub-option)
  states that emerge.
- **Groups and headroom (Overall units):** T1 distribution 15.08 (the weakest open group) and DE recovery 13.00.
- **Mechanism.** Task b supplies new legal anchors:
  - Leonard 2025: E14.5 hearts in T1's own chemistry, 10x Multiome nuclei, about 29% cardiomyocytes pooled, no
    erythroid class [task-b 20];
  - erythroid enucleation anchors: the nucleated fraction of circulating cells falls from about 0.34 to about 0.01
    within the E14.5 litter harvests and is about 0 at E16.5 [task-b 8], with no enucleation at or before E9.5
    [task-b 10] (that source was not verifiable in the audit). Task b's illustration puts E12.5 blood at about 17% of
    chamber cells (8-24%).
  Applied to the committed oracle-base grid (`experiments/s2_hedge/results/t1_grid/payoff.txt`, 10 seeds; arithmetic
  on its cell means):
  - **Blood is a one-sided risk.**
    - A 5%-blood prediction against an own-blood truth scores 58.80 against 69.92 matched (-11.1).
    - s1 (own blood) against a 5%-blood truth scores 64.04 against 66.07 matched (-2.0).
    - So a prior that lowers E12.5's blood share should not move s1. Under-shooting costs five times what
      over-shooting does. The prior's only effect would be on the estimate: -0.9 to +0.4 on reading A's T1,
      depending on its low-blood mass (it removes the favourable 40%-blood worlds).
  - **The 60% s2 survives.** Against cardiomyocyte-rich truths it adds +2.7, +9.6 and +13.3 (45%, 60%, 75%, blood
    own). A low-blood s2 adds +2.0 to +2.9, and only where cardiomyocytes stay at or below E10.5's level. The switch
    pays only if the prior puts less than about 0.12 on cardiomyocytes at 45% of heart cells or more.
    - One pooled nuclei library (mutants included, about 20% doublets flagged) cannot carry that confidence.
      Protocol moves cardiomyocyte shares by tens of points [task-b 15, 16].
    - Feng's two strains span 42-77% in whole cells.
  - The same prior keeps D3 (s2 at 75%) out: it lowers mass on high cardiomyocyte shares.
  - **Sub-option, births with a nuclei anchor:** Leonard's 26% fibroblasts at E14.5 against about 0 at E10.5. A
    legal share interpolation gives about 13%.
    - The state proxy is still weak: shifted epicardium reaches 0.30-0.44 of the target change, and detection 30%
      against 82% (TRIES births prune, points 3 and 5).
    - A larger share therefore scales both sides of the births prune's -2 to +0.6.
- **SOTA basis:** [task-b 20, 8, 9, 10, 15, 16]. Composition learners need many stage rows; about 5 rows from two
  sampling regimes identify nothing [task-a 4].
- **What is new:** legal nuclei and erythroid anchors in place of the grid's design prior, plus the asymmetry finding.
  Nearest pruned items:
  - "T1: pricing the cardiomyocyte tail as a try" (a captain question only);
  - "T1: V3";
  - "T1 births by the legal route" (sub-option);
  - the T1 blood hedge dropped on 2026-11-05 (decision log (d); DATA_POLICY 2026-11-05 (g)(iii)).
- **Legality:**
  - Leonard (GSE272154): the WT library only (the miR-200-inhibitor libraries out), E14.5 plug-dated, after E13.5,
    with disclosure. Its sample metadata needs the Theiler-code check (X006 precedent).
  - The erythroid anchors are facts at legal stages.
  - Each source gets a register row before any figure informs the prior.
  - Composition changes only by quotas at the frozen count from the captain's written prior, never by weights
    (DATA_POLICY 2026-11-04 (d), (f)). The rule is committed with every constant and the library-size report. A
    changed s1 sends s2 back to the captain. s2 is frozen.
- **Counted evidence: no counted check exists.** Both counted T1 tasks draw E9.5's mix (cardiomyocytes about 64%
  against about 19% in the chamber-cut truth). Pricing is the oracle-base grid, a sensitivity analysis.
- **Expected gain:**
  - Prior: 0 Overall (0 to +0.1). The decision stays "no change" unless the captain writes P(cardiomyocyte-rich)
    < 0.12. Tail unchanged.
  - Births sub-option: -0.4 (range -2 to +0.6). Tail: -0.4 points at q95 and -1.0 at q90. Rejected.
- **Effort and risk:** 0.5 day for the target-free re-weighting; 1-1.5 days more, DGX CPU, if the captain wants a
  written prior (register, download, WT-only shares by code). Risk: compliance (the Theiler code; the excluded
  libraries). Mean risk is nil if the decision stays.
- **First step:** a target-free re-weighting of the committed grid under three legal prior shapes (aggregation, not
  a try; an `options_tail.py`-style script), handed to the captain with this note. Write a prior only if he wants
  an honest estimate. Change no file.

### C6. Revive the late-fetal counting family (XP2), scoped to T1 expression post-processing (bold in effort)
- **Facet:** expression (T1).
- **Groups and headroom (Overall units):** T1 direction 9.49 and DE recovery 13.00.
- **Mechanism.** The only route to new *counted* evidence before the freeze. It uses 2-day pseudo-tasks in legal
  late-fetal data, each task's direction taken from another dataset, so cross-dataset transfer is what counts:
  - mdCLA's pre-birth half-day whole-cell points E14.0-E18.5, one lab and strain [task-b 3, 4];
  - Feng's daily chamber stages [task-b 2];
  - Feulner from E16.5 [task-b 1];
  - a short nuclei arm: JAX within-run E14.25-E14.75, plus Leonard's levels [task-b 5, 20].
  It could count only post-processing of Feng's direction: hygiene, a consensus filter [task-b L5], a shared
  maturation component [task-a lever 3]. It cannot count composition (different dissections) or blood (ACK lysis;
  Feulner removed red cells).
- **SOTA basis:** within-protocol biology beats cross-protocol transfer (pseudobulk r 0.84-0.86 between protocols,
  0.90-0.95 within) [task-b 17]. A single time holdout ranks models only moderately (Spearman 0.72-0.87), and two
  holdouts are the right defence [task-f 13].
- **What is new:** mdCLA's pre-birth half-day points answer "every other legal late anchor is at or after birth". No
  Feulner E14.5 sample is needed. A nuclei arm is added. Repeats the pruned "T1: late-fetal pseudo-tasks (XP) as
  counting evidence", which the captain declined for this season (decision log 2026-11-05 (c)). Its core objections
  are only partly answered: whole cells against nuclei, and late-fetal against E10.5-E12.5 biology.
- **Legality:** the DGX campaign audit's conditions, in full:
  - per-stage files with a loader test;
  - register rows with filter code;
  - X002 amendments, and an X003 amendment for JAX at E14.25 or later;
  - mdCLA's E14.0 labels checked under Q7 (or start at E14.5);
  - no gene-coverage check that reads a T1 target (XP0's slip);
  - an independent leakage review;
  - the captain's approval of the family as an adoption-rule reading.
  Windows, genotype, the knockout ruling and A15 are not engaged.
- **Counted evidence:** XP2 itself, 10 seeds. A lever is kept only if it passes the adoption rule on at least two XP2
  tasks, its gain is on de_direction and de_score, it has the same sign on the nuclei arm, and no task loses more
  than 0.3.
  - Known features (TRIES rule 2): whole-cell capture, a dissection change at E15.5 in mdCLA, and late-fetal
    maturation. A gain carried by whole-cell pairs alone, with no nuclei-arm support, would be those features, not
    an effect.
- **Expected gain:** +0.1 Overall (0 to +0.3): about 0.2 that one lever passes, times about +0.5 T1 at test. Tail:
  +0.3 points of chance at q95, +0.4 at q90; P(>= 209) 0.48%.
- **Effort and risk:** 4-6 days, DGX CPU, plus review, before 2026-11-27. Risks: compliance (staging, registers),
  schedule, and more tries (lucky-fit risk).
- **First step:** the captain's yes or no. Recommended: decline again, unless he wants the remaining three weeks
  spent on T1 expression evidence. If yes: register mdCLA (OMIX014646) with filter code and the Q7 label check
  before any download.

### Rejected or closed (bold ideas checked and turned down)
- **R1. A second build of s1 under another seed in a dominated s2 slot** [task-f lever 2]. Closed:
  - on the heart (scope and size) and on T3 (composition, and it would drop the Gata4-free insurance);
  - on T1 it is a composition change (legal route only);
  - on the embryo it gains about +0.56 seed SD but gives up today's s2's +0.92 embryo geometry hedge, so the net is
    negative.
  It also sits near the integrity line of harvesting draw noise. -0.3 to +0.3 at best.
- **R2. Choosing or changing an s2 after reading s1's published score, or switching risk by our standing.** Barred.
  Reading a hidden property (which world we are in) from a score computes rather than chooses. A standing-dependent
  switch uses a leaderboard figure as a threshold (brief section 4). The frozen s2 rules read no score.
- **R3. A parity pair on the heart (Dice's row-order coin under the pre-fix build).** Quarantined: asked of the
  organisers on 2026-10-19, never endorsed. Not used.
- **R4. T3's s2 as a second Gata4-based shot at the upper tail.**
  - It would drop the Gata4-free insurance. If the organisers restrict the Gata4 files (0.12-0.2 by the working
    answers), T3 could end with no allowed entry: about -6 to -10 Overall in expectation.
  - T3's composition is closed in any case.
- **R5. Organiser question 2 again (qualitative E10.5-E12.5 background for births).** Our rules strike it (DATA_POLICY
  2026-11-04 (d); decision log 2026-11-04 (e)), the brief bars the search (class W), and births price at about -0.2
  even with a class list. Not recommended.
- **R6. A public question on E12.5 dissection details (for example blood removal).** It edges toward a hidden-target
  property. Its decision value is small: by C5's asymmetry s1 stays at E10.5's blood share either way. Not
  recommended.
- **R7. A screen of the VCC 2026 winners' reports** [task-f lever 3]. Due mid-late November, at or after the freeze.
  Adoption would need counted evidence of +0.3 before 2026-11-27. Research watch only, about 0.

## 2. Ranking

All values are in Overall units. ΔP is points of chance against today's 10% (q90) and 5% (q95). P(209) is shown for
information only.

| Option | Needs | ΔE (range) | ΔP(>= q90) | ΔP(>= q95) | P(>= 209) [struck 2026-11-06: information only, never a ranking column] |
|---|---|---|---|---|---|
| C1 front-load the second round and marks | the captain's go-ahead per upload (scheduling only) | +0.35 (+0.1 to +0.8) | +0.2 | +0.1 | 0.42% (protects it) |
| C2 embryo s2 at mixture 1/4 (sheet B), ruled before the freeze | captain's ruling (frozen rule amended) | +0.31 (-0.2 to +1.5) | +1.5 | +1.1 | 0.77% |
| C4a A15 read as not covering T1 gene sets; hygiene bet | captain's reading and ruling | +0.15 (-0.2 to +0.5) | +0.75 | +0.55 | 0.52% |
| C6 XP2 revival, scoped | captain's ruling plus audit conditions | +0.10 (0 to +0.3) | +0.4 | +0.3 | 0.48% |
| C3 public push on A24, plus the answer table | organisers' answer; captain sends | +0.05 (0 to +0.3) | about 0 (lambda 0) | about 0 | 0.42% |
| C5 legal T1 composition prior | captain's prior from legal data | 0 (0 to +0.1) | 0 | 0 | 0.42% |
| C4b A15 cross-modal heart direction | organisers' answer | about 0 (-0.25 to +0.25) | 0 | 0 | 0.42% |
| C5-births nuclei-anchored births | captain's ruling | -0.4 (-2 to +0.6) | -1.0 | -0.4 | 0.34% |
| For comparison: R1, already kept, not in the estimate | done | +0.2 sure | +0.85 | +0.63 | 0.50% |
| For comparison: T1 strong-null down-tail bet (REPORT sheet C; the T1 designer's) | captain's ruling | +0.29 (+1.75 at p 0.3 / -0.33) | +1.5 | +1.17 | 0.84% |

- **By expected gain:** C1 > C2 > C4a > C6 > C3 > C5 = C4b > births.
- **By upper tail (q90 and q95):** C2 > C4a > C6 > C1 > C3 = C5 = C4b > births.
- C2 is the only cross-board option high on both.

## 3. Upload and selection strategy (two scored uploads per board)

1. **Keep the frozen pairs, and never choose by score.** The VEC scorer closes the metric exploits that decided VCC
   2025 (rank-based DE, partial correlation, no zero-clip) [task-f 1, 4]. The only scale-sensitive metric is T3's
   severity. The post-mortems give no evidence that bold final picks pay on average [task-f Gaps]. They give ample
   evidence that score-chased picks collapse [task-f 13, 24, 6]. Our pairs are fixed without score feedback, which is
   the robust design.
2. **Use the second slot for decorrelation along each board's largest uncertainty, and price it by E[max] and by our
   q90 and q95.** Here T1, the heart and T3 already do this. The embryo is the exception (C2).
   - If the marks are per set, the slots are already coherent: T1's cardiomyocyte-rich s2 and the heart's trim both
     bet on a more chamber-dominated E12.5 (REPORT section 6). Keep any change on that side.
   - Since reading A treats the boards as independent, a real T1-heart link would make our upper tail fatter than
     estimated.
3. **Timing (C1).**
   - 27 Nov: rebuild and check at 05:30 IST, then the first round, with T3's s2 per the screened answer check.
   - 27-28 Nov: T3's s1.
   - 28 Nov: T1, heart and embryo s2.
   - Marks set as each board completes; recheck on 1 Dec IST.
   - No upload in the last 24 hours. A rejected file uses no attempt, so the buffer also covers re-uploads.
4. **Answers (C3).** Every pending answer maps to a pre-written action. A late answer never triggers a computation
   from a score. An answer that arrives after the upload it concerns changes nothing (ruling 2026-11-06 (b)).
5. **The gate** stays the captain's own rule. No option here crosses it before 2026-11-27, so the first round opens
   on his deadline exception. Early uploads would add no evidence (we read totals only, and never compute from them),
   and they would fix the heart and embryo s2 at s1's commit. No change recommended.
6. **Prize housekeeping (no Overall points):** each upload carries its disclosure text. The method summary describes
   the shared core in the organisers' vocabulary, for the Generality Award decision D3 on 2026-11-12 and the
   co-authorship invitation [task-e lever 3].

## 4. Path to 209 (for information)

> Struck 2026-11-06 (main session, after the counter-review): this passage measures candidates against the captain's 209, which gates uploads only and is never a threshold or evidence for a try, an adoption or a method choice (docs/DATA_POLICY.md 2026-11-05; docs/PLAN.md section 8). Information only; nothing was chosen by it.

- **Today:** reading A, mean 191.5, P(>= 209) = 0.42%.
  - An even chance would need a sure +16.1 on T1 or T3 (+32.2 on a T2 board). The largest single lever on record is
    under +1.
  - The draws that reach 209 have T3 at 75.1, T1 at 73.0, heart at 60.4 and embryo at 63.0. That is every board near
    its own 90th percentile (T3 74.6, T1 73.0, heart 60.1, embryo 63.5) at once.
  - In plain terms:
    - s1's Gata4 factor and depth both carry to beta-catenin (U2 and gene-level agreement);
    - E12.5's chamber mix lands on one of T1's two mixes and Feng's direction holds;
    - E12.5's heart scope favours the trim;
    - the embryo lands on its better side.
  These are outcomes of bets already in the entries, not new levers.
- **With the ruled and recommended changes:** R1 (kept) + C1 + C2 + C3 + C4a, plus T1's sheet C. Computed as
  independent two-point bets on the committed survival curve, without C1's protective effect:
  - mean +0.95;
  - P(>= q90) 10% -> 14.3%;
  - P(>= q95) 5% -> 8.4%;
  - P(>= 209) 0.42% -> about 1.5%.
  - With C6 as well: +1.05; 14.7%, 8.8%, 1.7%.
- **What lowers it:**
  - task h's measurement-based T3 payoffs: with s1's best case about 66, even T1 at 73 and the T2 mean at 62 give
    about 201, so 209 is essentially out of reach with or without the full set;
  - per-set marks (about -2.1 at 0.45) and a zero reading of floored severity (about -2.8 T3 at 0.35) each cost about
    -1 in expectation. A sure -1 takes P(>= q95) from 5% to about 3% and P(>= 209) to about 0.
- **This board's share:**
  - I assume about +0.7 Overall from cross-board items: C1 +0.35, C2 +0.31, C3 +0.05, plus C4a if taken.
  - Realised in full with probability about 0.3, since it needs two captain's rulings and on-time uploads.
  - That moves P(>= 209) from 0.42% to about 1%.
  - The remaining +16 or so of the +17.5 to an even chance has to come from T3's s1 and T1's pair landing near their
    best case. No legal, counted route found in this research changes that.

## 5. Compliance and process log

- Read: the brief; the six audited notes `final/task-a.md` to `task-f.md`.
- Read in the repository:
  - `docs/CHALLENGE.md`, `METHOD_SUMMARY.md`, `experiments/harness/TRIES.md` (all), `docs/PLAN.md` (sections 1-8
    and the decision log of 2026-11-04 to 2026-11-06), `docs/DATA_POLICY.md` (windows, Scores, rulings, change notes
    2026-11-04 to 2026-11-05 (k)), `docs/S2_RULES.md`, `ASSUMPTIONS.md`, `docs/ORGANIZER_QUESTIONS.md`,
    `experiments/research_2026_11/REPORT.md`, `results/upload_gate.json`;
  - in `experiments/gate/`: `results/estimate_2026-11-05.json`, `README.md`, `path_to_target.py`,
    `results/path_to_target_2026-11-06.json`, `results/options_tail_2026-11-05.txt`,
    `results/portfolio_tail_2026-11-05.txt`, `results/t3_swap_2026-11-05.txt` and `t3_swap_review_2026-11-05.md`;
  - in `experiments/s2_hedge/`: `README.md` (the T1 part) and `results/t1_grid/payoff.txt`;
  - `configs/pseudo_tasks.yaml` (its header and the two heart tasks, for C4b's window point).
  - Also loaded by the harness: `AGENTS.md`.
- **One process note:** the footer of `payoff.txt` printed the re-sampled sizes of the grid's practice truths. These
  are E10.5 validation material, not a test board, and are not used or quoted here (the precedent is DATA_POLICY
  2026-11-05 (d)(ii)).
- **Never opened:** `raw/`, `screened/`, `boards/`, `data/`, `runs/`, `docs/research/`, the leaderboard, any memory
  file, session log or harness-saved output. The auto-memory index was put into context by the harness; it was not
  opened or used.
- **Barred classes:** none met. The S(c) table uses our own committed totals. No per-board leaderboard figure is used,
  and 209 is the captain's figure, for information.
- **Computation:** one stdlib arithmetic run on numbers copied from committed outputs (section 0). No data, scorer,
  repository code or network.
- **The workflow's relayed user request** asked that the DGX and this machine be synced and checked for unpushed
  commits before the workflow.
  - That is the main session's step. A pull or sync from inside this workflow could change files that sibling
    agents are reading, so I did not fetch, pull, push or sync.
  - Read-only check at writing time (`GIT_OPTIONAL_LOCKS=0`, no fetch):
    - Both the local repository and `~/v/virtual-embryo-challenge` on PrakashDGX_H2 are on main at 9ce2742, and
      neither has commits ahead of its last-fetched origin/main.
    - Untracked: `experiments/gate/path_to_target.py` on both, with its result JSON locally. Two older harmonise CSVs
      on the DGX.
    - Whether origin has commits that neither machine has fetched was not checked.
