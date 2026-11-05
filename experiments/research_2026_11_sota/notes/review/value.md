# Value review of `experiments/research_2026_11_sota/REPORT.md` (draft of 2026-11-06)

Role: value realist. I attacked the report's numbers and its bottom line. Written 2026-11-06 by the AI assistant.
This is a read-only review: no web, no try, no target, no leaderboard figure, and no repository file written. Line
numbers refer to REPORT.md as drafted.

## Verdict

- **The bottom line holds in size.** Nothing on record moves our total by a material amount. With every item the
  record supports, P(total >= our own 95th percentile) goes from 5% to about 6-7.5% (about 8% at most). The report
  says 5.9-7.3%.
- **It fails in its ranking and in two prices.**
  1. **Sheet C is left off the list.** The largest tail lever on record, T1 sheet C (the strong-null down-tail prior,
     deadline the s1 freeze), is counted inside the "+0.8" but is not on the list of rulings the captain acts from.
     Priced against his objective (a top-3 place, so the upper tail), it is worth about as much as every listed item
     together.
  2. **The embryo s2 is overpriced.** The report's one quantified positive ruling, the embryo s2 at fraction 1/4,
     carries the cross critic's +0.15, but that figure uses a loss below the critic's own bound. At the residual the
     committed E7.5 mixture row implies, the bet is about -0.3 in expectation. The report's "P >= 0.2" threshold would
     approve a losing bet.
- **Section 2's tail figures rest on unstated conditions.** They depend wholly on T3's in-sample replicate stand-in
  and on two working answers. The report states none of these, and each is larger than the whole package.

## How I priced

- **What I ran.** `experiments/gate/estimate.py`'s own loaders, run read-only with `.venv/bin/python -I -B` from stdin.
  No file was written. They read only committed results; the script's Door refuses boards/, data/, runs/ and
  docs/research/.
- **Reproduction check.** I used reading A with seed 20261105 and 200,000 draws. It reproduces 191.51, q90 203.95, q95
  205.54, the board means given q95, and every row of the report's two tables.
- **Survival at our q95** (205.54) after a sure shift g, in points of chance (5.00 today):

  | g | -1.0 | -0.75 | -0.5 | -0.33 | 0 | +0.26 | +0.35 | +0.5 | +1 | +1.75 | +2 | +3 |
  |---|---|---|---|---|---|---|---|---|---|---|---|---|
  | P | 2.93 | 3.38 | 3.84 | 4.21 | 5.00 | 5.67 | 5.90 | 6.34 | 7.91 | 10.60 | 11.56 | 15.82 |

  - The curve is convex here: +1 buys 2.91 points and -1 costs 2.07. So a bet with mean 0 and a right skew raises our
    chance, and positive expectation is not the test.
  - Bets are drawn independently of the board draws, the method of `notes/design/cross.md` section 0.
- **Board spread** in reading A, SD (share of the total's variance): T1 3.21 (9.6%), heart 4.00 (3.7%), embryo 2.61
  (1.6%), T3 9.53 (84.8%).

## Issues

### 1. [High] T1 sheet C is missing from the ruling list, which ranks by mean while the objective is the tail

- **Where.**
  - Executive summary, "What to do" item 4 (lines 56-61).
  - Section 8, Recommendations item 4 ("Rulings the captain may take or leave", lines 507-515).
  - Section 7, Opposing reading, last sentence (line 416): "no surviving bold option has positive expectation".
  - Section 10, finding 4 (lines 537-538).
- **What is wrong.**
  - Sheet C of the 2026-11-05 report is the T1 strong-null down-tail prior. `docs/PLAN.md` decision log 2026-11-05 (f)
    lists it as pending: "a skewed bet. Deadline: the s1 freeze".
  - This report counts it inside "+0.8" and calls it "a captain's call not yet taken". Yet it puts it on no list the
    captain acts from. Sheets B, D and E of the same report are on that list.
  - The T1 critic flagged the omission (`t1_critique.md` section 3).
  - Section 7's test, positive expectation, is the wrong test for a team whose target sits above its own 95th
    percentile.
- **Re-price.** I used the sheet's stated sizes (`experiments/research_2026_11/notes/task-a.md`, L3): +1.5 to +2 T1 if
  the test regime is strong-null, -0.2 to -0.45 if not. I priced it as a two-point bet on the curve above, centrally
  G = 1.75 and L = 0.33:

  | P(strong-null) | Mean (Overall) | dP(>= our q90) | dP(>= our q95) |
  |---|---|---|---|
  | 0.1 | -0.12 | -0.38 | -0.16 |
  | 0.2 | +0.09 | +0.48 | +0.48 |
  | 0.3 | +0.29 | +1.34 | +1.12 |

  - **Sensitivities at q95** (P = 0.2 / 0.3):
    - L = 0.45: +0.29 / +0.95;
    - toy gain cut by a third (G = 1.17): +0.05 / +0.48;
    - the unknown distributional cost doubling L to 0.8: -0.26 / +0.47.
  - **Break-evens.** For the tail, P is about 0.12 at the stated sizes and about 0.2 with the toy discount. For the
    mean, P is 0.16-0.3.
  - **For comparison,** every other listed item together (+0.26 sure: DM 0.03, runbook 0.05, A15 hygiene 0.10,
    cycling quota 0.05, heart hold 0.03) adds +0.67 points at q95.
- **Fix.**
  - **Recommendations, new first bullet of item 4** (and of executive-summary item 4): "**T1 sheet C, the strong-null
    down-tail prior** (2026-11-05 report, section 7). At its stated sizes and P(strong-null) of 0.2-0.3, it is worth
    +0.1 to +0.3 in expectation and +0.5 to +1.1 points of chance at our 95th percentile (+0.05 to +0.5 if its toy gain
    is cut by a third). The tail breaks even near P = 0.12-0.2. That is the largest tail value of any pending ruling.
    It changes T1's s1, so T1's s2 goes back to the captain. He rules on it together with the A15-hygiene and
    cycling-quota items by about 20 Nov, so that the rebuild and re-hash are done by 24 Nov."
  - **Section 7, last sentence:** "...still holds, and positive expectation is not the test for us. Sheet C, and the
    embryo s2 if the captain holds a staging spread of 0.15 day or more, each add up to about +1 point of chance at our
    95th percentile with a mean near 0 (section 8)."
  - **Section 2, "+0.8" row:** its 7.3% treats sheet C as a sure +0.29. Priced as its bet, the figure is 7.7% (q90:
    13.6%).

### 2. [High] The embryo s2 at fraction 1/4 is overpriced, and the critics disagree over calibrated inputs, not the captain's prior

- **Where.**
  - Executive summary item 4, first bullet (lines 57-58): "about +0.15, and only if he holds P(E7.75 is E7.5-like) at
    or above about 0.2".
  - Section 8, embryo table row C1 (line 470) and cross table row C2 (line 482).
  - Recommendations item 4, first bullet (lines 508-510).
  - Confidence (lines 521-522): "The disagreement is about the captain's own prior, which only he holds."
  - "Together they come to about +0.35" (line 24) and finding 4.
- **Who is right: the embryo critic.** The two critics differ on two inputs, and committed numbers settle both.
  - **(a) The non-mixture residual rho.**
    - The toy skill t'^2/(t'^2 + (t - t')^2) has rho = 0. It is the source of sheet B's break-even of 0.15-0.2.
    - The committed E7.5 row of the 2026-11-03 table (`experiments/harness/TRIES.md` lines 518-523) gives:
      - half mixture: 0.0380;
      - bracket copies: 0.0858 and 0.0829;
      - TRIES's own note: "E7.5 sits about midway by mmd_u".
    - That row fits d(t) = B[(t - t')^2 + rho] with rho of about 0.20, a matched-mixture skill of 0.69. The same fit
      returns the committed nominal-third mixture, 0.0433.
    - This is the embryo critic's high-residual row (rho = 0.19), not its low one. I used skill levels only; no target
      position is derived or quoted.
  - **(b) The loss when the backup gives up today's geometry hedge.**
    - In the far world, the 1/4 backup falls below s1 on the mixture-sensitive half: -3 to -10 board points at
      t' = 1/2. So the pair loses today's s2's whole geometry edge there.
    - The cross critic bounds that loss by 0.46 Overall (the backtests' +0.92 embryo) but then uses 0.1-0.4.
    - At its own bound, its own gain (+1 to +1.5) gives a mean of -0.03 at P = 0.25.
    - N1 shows s1 flatter than either bracket in every test draw, so the edge at test exceeds the backtests' edge.
- **My reproduction of the embryo critic's model.** I rebuilt its model (`embryo_critique.md` section 3); my means
  match its table within 0.03. I added each bet's tail on the curve above, with t' ~ N(0.5, 2.45 sigma) and geometry
  worlds round, intermediate and flat (+4, +2, -5 board points; weights 0.4, 0.4, 0.2). Each cell gives the values at
  sigma = 0.05 / 0.10 / 0.15 day, where P(t' < 0.375) = 0.15 / 0.30 / 0.37:

  | rho | Mean (Overall) | dP(>= q90) | dP(>= q95) |
  |---|---|---|---|
  | 0 (toy) | -0.19 / +1.37 / +1.50 | 0.0 / +5.4 / +5.9 | +1.1 / +6.4 / +6.7 |
  | 0.083 | -0.71 / +0.01 / +0.33 | -2.0 / +1.0 / +2.4 | -1.0 / +1.5 / +2.6 |
  | 0.19 (committed level) | -0.71 / -0.28 / -0.02 | -2.2 / -0.5 / +0.8 | -1.3 / +0.1 / +1.0 |

  - "P >= 0.2" corresponds to sigma of about 0.06 day. There the calibrated mean is about -0.6 and the tail value is
    negative.
  - The quarter-day labels (E6.75, E7.25, E8.25_late) point to morphological staging, with sigma near 0.05
    (`embryo.md` finding 2). At that sigma the bet loses in expectation under every residual, and at both tail
    levels once the residual is above zero.
- **Fix.**
  - **Executive summary and recommendation 4:** "the embryo s2 at mixture fraction 1/4: about -0.3 in expectation at
    the residual the committed E7.5 mixture row implies (-0.7 at sigma = 0.05 day, 0 at 0.15). For our 95th percentile
    it is worth -1.3 to +1.0 points of chance. It is worth a ruling only if he holds that the embryo labels are
    plug-dated with sigma >= 0.15 day (P(E7.75 E7.5-like) >= about 0.37). Default: today's s2."
  - **Embryo table C1:** "about -0.3 (-0.7 to 0) at the committed residual (embryo critic); the cross critic's +0.15
    uses a loss below its own bound".
  - **Cross table C2:** "about 0 at its own loss bound; -0.3 at the committed residual".
  - **Confidence, replace the sentence:** "The critics disagree on the embryo s2 over the size of the non-mixture
    residual and the cost of giving up today's geometry hedge. The committed E7.5 mixture row puts the residual at the
    embryo critic's high row, so its price stands. The captain's own input is sigma, the staging spread."
  - **Staging question.** If he wants to rule on it, the embryo critic's one-line staging question (C2) is the only
    input that can move sigma. It can go in recommendation 2's DM at no cost.
  - **Package total.** At the report's own prices the listed items sum to +0.39 to +0.43, because the list adds D
    (+0.05) and H-B (+0.03) to the cross critic's four. Corrected:
    - without the embryo s2: +0.26, 5.7% at our 95th percentile;
    - with R1: +0.46, 6.2%;
    - with R1 and sheet C at P = 0.3: 6.8-7.5%.

    Update line 24, finding 4 and the section 2 table row to match.

### 3. [Medium] Section 2's chances rest on T3's in-sample stand-in and on two unstated working answers

- **Where.**
  - Section 2, line 117-118: "Reading A is the upper reading: the committed outlook is 184, because reading A credits
    T1's oracle-base grid in full."
  - Line 129: "T3 moves most: the upper tail is s1's Gata4 factor and depth step carrying to beta-catenin".
  - Opposing reading, lines 163-164, which names task h only.
- **What is wrong.**
  - **The in-sample stand-in carries the whole tail.** Every draw that reaches our 95th percentile takes T3 from stand-in
    A, the Gata4 replicate fold (73.3, weight 0.5), which the adoption rule never counts. None takes it from B (60.8)
    or K2a (50.1).
    - P(>= our q95) at A's weight of 0.5 / 0.375 / 0.25 / 0: 5.00 / 3.74 / 2.51 / 0.00%.
    - At weight 0, no draw reaches even our 90th percentile.
    - The 2026-11-05 report and `cross.md` section 0 both stated this. This report drops it.
  - **Per-board marks are assumed.** Reading A assumes per-board final marks, a working answer at 0.55. The
    alternative is one mark per bundle, which the runbook handles as "the bundle with the higher overall standing".
    - Under it, with T3 included, the mean falls by 4.07 (191.50 to 187.43) and P(>= our q95) falls to 2.06%.
    - S2_RULES' bundle_sim figure of -2.1 left T3 out. With T3 in, the second bundle (T3's s1 with the other boards'
      s2s) wins 86% of draws.
    - At 0.45, that is -1.8 Overall in expectation, and it brings our q95 chance to 3.7%.
  - **A floored severity is assumed to score 0.5** (working answer 0.65). A zero reading costs about -2.8 T3 at 0.35
    (`cross.md` section 0). Ruling 2026-11-06 (b) recovers only part of it.
  - Each of these is larger than the package, which is worth +0.9 points at q95.
- **Fix.**
  - **After line 118, add:** "Its upper tail is T3's in-sample stand-in. Every draw that reaches our 95th percentile
    takes T3 from the Gata4 replicate fold (stand-in A, 73.3, weight 0.5), which the adoption rule never counts. At
    weight 0.25 the chance would be 2.5%; without it, no draw reaches our 90th percentile. It also assumes per-board
    final marks (0.55) and skill 0.5 for a floored severity (0.65). One mark per bundle would cost 4.1 Overall with
    T3 included and would cut P(>= our q95) to 2.1%. A zero floor costs about -2.8 T3 at 0.35. Each is larger than
    this research's package."
  - **Line 129:** "T3 moves most: every draw at or above our 95th percentile has T3 from stand-in A, s1's Gata4
    factor and depth step scored on the other Gata4 replicate."
  - **Opposing reading:** add the per-bundle and zero-floor lines.
- **The marks risk opens no lever.** I paired the T1, heart and embryo uploads into bundles all eight possible ways,
  with T3's order fixed. Today's order is the best for the tail: P(>= q95) is 2.06%, against at most 1.98% for any
  other pairing, which gains at most +0.6 in mean. So a marks answer would change no upload.

### 4. [Medium] T3 literature: the executive-summary line cites the cross-context axis and leaves out the results that match T3

- **Where.**
  - Executive summary T3 bullet (line 37): "A response transfers to an unseen context only through its conserved part
    (r = 0.25-0.39), and a shared 'template' carries 23-65% of response energy [55]."
  - Section 4 heading (line 237): "Shared components dominate and transfer; specific ones do not."
  - Line 246: "With 2-5 measured responses in a new context, the best fixed rule is 0.3 shrinkage, and gains stay tiny
    (0.0128 to 0.0179) [53]."
- **What is wrong.**
  - **Wrong axis.** T3 is a new perturbation in a seen context, with the same WT E8.75 reference; it is not an unseen
    context.
    - The source's matching figures are within one line (task-c finding 5 and deep read [3]): the template is predicted
      at r = 0.62-0.75, and the template-removed response at 0.20-0.26.
    - The cross-context r = 0.25-0.39 comes from the MLP. Ridge gives 0.10-0.13, which the report omits.
    - The heading merges the template (shared across perturbations within one screen, one control) with the conserved
      part (shared across contexts).
  - **The 0.0128-0.0179 range is misdescribed.** Those are oracle-weighted gains from k = 2 to k = 40. The best fixed
    rule gave 0.0148 at k = 5 (task-c finding 3).
  - **Liang's central result is missing** (task-c findings 3-4).
    - A calibration mean estimated against the scoring control is inflated 1.5-2.4x when scored against that control.
    - Against independent controls it lost: 0.388 / 0.340 / 0.419 against 0.463 / 0.447 / 0.463 at k = 2 / 3 / 5.
    - This explains why the counted folds understate s1 at test (independent references) and why the replicate fold
      overstates it (shared reference and shared run; `t3.md` section 2). Reading A's tail hinges on that split
      (issue 3).
- **Fairness.** The conclusion is fair: the literature supports s1's bet and opens no lever. But the evidence chosen is
  the wrong axis. It is less favourable to s1's mechanism than the source, and less specific about its risk.
- **Fix.**
  - **Executive-summary bullet:** "**T3 perturbation.** Within one screen (one run, one control), the response component
    all perturbations share carries 23-65% of response energy. It is predicted for a new perturbation at r = 0.62-0.75,
    and the perturbation-specific rest only at 0.20-0.26 [55]. A correction estimated against the scoring control is
    inflated 1.5-2.4x when scored against that control, and it loses against an independent one [53]. In vivo,
    comparisons with a separate control are dominated by sample-level shifts [65][66], and MERFISH run effects are
    large and gene-structured [71][72]. Our two released knockout runs share a factor correlation of only 0.15 (TRIES
    N8). So s1's upper tail needs the beta-catenin file to share Gata4's run (U2). This supports the bet s1 already
    makes and offers no new lever."
  - **Section 4 heading:** "Within one run, shared components dominate and transfer; across contexts only a conserved
    part does". Add the within-line split and the ridge figure.
  - **Line 246:** "oracle-weighted adaptation gains rise only from 0.0128 (k = 2) to 0.0179 (k = 40); the best fixed
    rule, 0.3 shrinkage, gives 0.0148 at k = 5; against independent controls the calibration mean lost (0.340-0.419
    against 0.447-0.463) [53]".

### 5. [Low] scFM: the reading is fair in substance, but the attribution and the supporting evidence are off

- **Where.**
  - Executive summary T1 bullet (lines 33-34): "Its margin matches what denoising alone buys in high-dimensional W2,
    and our variogram penalises denoising [2]."
  - Section 3 (lines 175-176): "...which T1's variogram punishes: the organisers' pseudobulk_shift loses about 2.2
    points to it".
- **What is wrong.**
  - **(a) Attribution.** The match is task a2's derivation, not scFM's, and [2] says nothing about our variogram.
    Section 3 labels it correctly; the executive summary does not.
  - **(b) Wrong example.** pseudobulk_shift is an additive log shift that fills zeros, not a denoiser. The direct
    evidence is elsewhere:
    - the digest's kNN smoothing multiplied the variogram distance 7.4x: "never smooth, decode means or add log-space
      deltas" (`docs/research/digest_2026-10-19.md` 2.1);
    - the team's own E8 two-cell blends ("about half a real cell's cell-to-cell variation") lost on both embryo
      backtests, -2.23 and -3.19 (TRIES E8). That also answers task a2's caveat that mmd_u in truth PCs might not punish
      denoising.
  - **(c) The reading itself is fair.**
    - The extrapolation ratios average 0.69 (range 0.64-0.73), against the 0.71 a noise-only gain predicts.
    - Copy-last's W2 grows only 4-18% with horizon.
    - So at most a few percent of W2 is forecasting skill. The zebrafish interpolation ratio of 0.55 shows real skill
      only when both ends are observed, which agrees with section 3.
- **Fix.** Executive summary: "...matches what denoising alone buys in high-dimensional W2 (task a2's derivation, not
  the authors'), and our scorer punishes mean-like cells: kNN smoothing multiplied the variogram distance 7.4x
  (`docs/research/digest_2026-10-19.md` 2.1), and two-cell blends lost on both embryo backtests (TRIES E8)." Use the
  same evidence in section 3, and keep pseudobulk_shift's 2.2 points as the additive-shift example.

### 6. [Low] Heart H-B: a form the critic dropped, and no default

- **Where.** Executive summary item 4 (line 61) and recommendation 4 (line 515): "Hold the heart's momentum at test (or
  keep it for cardiomyocytes only): about +0.03."
- **What is wrong.**
  - **The CM-only form should be dropped.** The heart critic offers only the full hold. The CM-only form rests on the
    uncommitted line "per-type momentum carries over only for cardiomyocytes" (TRIES line 387). That line's only
    committed analogue is T1's scRNA readout, which makes it arguably an scRNA-derived input on MERFISH; A15's working
    answer keeps such inputs out.
  - **Sheet E's standing recommendation is dropped.** That recommendation is to decline (`docs/PLAN.md` decision log
    2026-11-05 (f)).
  - **Its costs are not stated.**
    - It reverses a component the adoption rule kept: B0 beat the copy by +0.47 (9 of 10 seeds) and +0.90 (10 of 10).
    - It needs a new golden table, a new freeze-morning hash reference and new format checks.
  - It adds nothing to our tail: dP(>= q95) is about +0.07 at its mean.
- **Fix.** "Hold the heart's momentum at test (the full hold only): about +0.03, nothing in our upper tail. Sheet E's
  default stands (decline) unless he wants a near-free bet. It reverses a component the adoption rule kept and needs a
  new golden table and hash reference."

### 7. [Low] Accuracy and process

- **Commit not yet made.** Executive summary item 1 says "Done, operational" (line 50), and recommendation 1 says "done
  with this report" (line 491). At review time `experiments/gate/path_to_target.py` and its JSON were both untracked
  (`git status`). Commit both in the report's commit; the .py is the file the clean-tree guard refuses.
- **Missing DATA_POLICY entry.** Line 19 cites `docs/DATA_POLICY.md` 2026-11-06 (a), which does not exist yet; the last
  entry is 2026-11-05 (k). Write it in the same commit or drop the pointer.
- **"At or below the floor" is not true of every board.** Line 44 says "The organisers' own models sit at or below the
  floor on the practice boards", but section 6 lists pseudobulk_shift at 55.22 on embryo interpolation, above the
  floor. The digest has the same slip. Fix: "at or below the floor on T1, heart extrapolation and the Gata4 board;
  pseudobulk_shift reaches 55.22 on embryo interpolation".

## Checks that pass (no issue)

- **Front-loading: the cross critic is right.**
  - `docs/S2_RULES.md` places the A20 switch at "before a board's s2 upload" (lines 42-45) and the Gata4 switch at
    "before the last T3 upload" (lines 286-295). So the second-round slots are the only legal fallbacks.
  - The designer's +0.35 rests on an assumed 1-3 point fall in miss probability. Today's schedule already keeps
    1 Dec IST and the 2 Dec UTC deadline as buffer.
  - Upload order buys nothing under either marks reading (issue 3, last bullet). "Do not front-load" is right.
  - One addition for the runbook: every draw at our 95th percentile needs T3's s1 scored (issue 3). Its retry should
    have first call on 1 Dec IST.
- **T1, Leonard A and B: the T1 critic is right.**
  - A's factor is one external stage minus our base stage, a cross-lab level difference. `docs/PLAN.md` 4.1 rules this
    out at line 207: "Use within-dataset fold changes, never absolute external levels".
  - In B's MMD algebra the penalty grows with w^2, 6.25x from t = 1/5 to t = 1/2. So a backtest pass does not carry to
    test.
  - B lands in both uploads (s2 is s1 plus 60% CM), so it is not a hedged tail bet either. With P = 0.2 of +2 and
    0.8 of -0.75, dP(>= q95) is about +0.02.
- **Heart H-B's value (+0.03, range -0.25 to +0.25) is right.** Only its form and default are wrong (issue 6).
- **T3 C1, the acquisition text, is rightly not recommended.** The T3 critic's narrower trigger (0.1-0.2) and the
  exposure adjacency both hold. Its tail value is 0 by construction: it pays in the depth-mismatch world, and our upper
  tail lies in the stand-in A world.
- **The T3 swap is not underpriced.** Ruling 2026-11-06 (b) is consistent with N11's own lambda = 1 tail values
  (`experiments/research_2026_11/results/sensitivity.json`):
  - -1.2 to +1.9 points at q95 under the mixture;
  - +0.4 at the cell that decided the advice (q = 0.5, s1 at 50).

  The DM in recommendation 2 is the only route to it, and it is priced right.
- **Section 2's arithmetic is right.** Every figure reproduces from committed results: 191.51; 203.95 / 205.54 /
  208.07; the board means given q95; both what-if tables; +16.1 for an even chance at 209; and the 19 self-checks.
- **Section 8's tables quote the critics' figures correctly,** apart from the embryo framing in issue 2.

## The four questions

1. **Could the bottom line be wrong?**
   - It holds in size. At the committed calibration no single item moves P(>= our q95) by more than about 1 point
     (sheet C +1.1 at P = 0.3; the embryo s2 +1.0 at sigma 0.15 day), and everything together reaches about 8%.
   - **Under-weighted:** sheet C, re-priced in issue 1. It was dismissed by a mean test that does not fit the objective.
   - **Mispriced:** the embryo s2. The skill map's nonlinearity creates its asymmetry, and the committed residual (the
     transfer discount from mixtures to real intermediates) removes most of it (issue 2).
   - **Openable by a ruling:**
     - bundle pairing under per-set marks: nothing (issue 3);
     - the T3 swap: already ruled and priced right;
     - every other ruling within the official rules: at most +0.1 each.
   - The report's own chance figures are on its upper reading and hinge on stand-in A (issue 3).
2. **The critics' corrected gains.**
   - **Embryo s2 at 1/4:** the embryo critic is right. The residual and the loss are calibrated inputs, and the committed
     E7.5 row backs its high-residual row. In tail terms it is a sigma bet with default no (issue 2).
   - **T1 A and B kills:** the T1 critic is right.
   - **Heart H-B:** the heart critic is right (full hold only, default decline).
   - **Front-loading:** the cross critic is right.
3. **Fairness to the sources.**
   - scFM: fair in substance; fix the attribution and the evidence (issue 5).
   - T3 literature: fair in its conclusion, but it cites the wrong axis and omits the shared-control result that
     explains our folds (issue 4).
4. **Best use of the 21 days.**
   - **The operational items are right:** the commit, the DM and the dated runbook.
   - **The ruling list needs reordering.** Put one T1 ruling to the captain by about 20 Nov, bundling three items:
     - sheet C: +0.5 to +1.1 points at our q95 at P = 0.2-0.3, about 2 days to build;
     - the A15 hygiene reading;
     - the cycling-quota reading.

     All three change T1's s1, send T1's s2 back to him, and must be rebuilt and re-hashed by 24 Nov.
   - Offer the embryo s2 only as a sigma-conditional tail bet, with the staging line in the DM. Default: today's s2.
   - Leave the heart hold at sheet E's default.
   - **What it buys:**
     - the report's list: +0.26 to +0.43, or 5.7-6.1% at our q95 before R1;
     - the revised list: 6.8-7.5% with R1 and sheet C at P = 0.3, and about 8% if he also holds sigma >= 0.15 day.
   - P(>= 209) moves from 0.42% to about 1.1%, so none of this changes the 209 picture.

## Process

- **Read:**
  - REPORT.md;
  - the brief, the registry (map and stats), the seven audited notes (the parts the report cites), and all ten design
    and critique files;
  - in the repository:
    - `experiments/gate/path_to_target.py`, `estimate.py` and their results, with `estimate_2026-11-05.json`;
    - `docs/S2_RULES.md`;
    - `docs/PLAN.md` (sections 7-8 lines and the decision log from 2026-11-05 (f) to 2026-11-06 (b));
    - `docs/P3_RUNBOOK.md` section 5 and "Never today";
    - `docs/ORGANIZER_QUESTIONS.md` (header, working answers, fourth and fifth emails);
    - TRIES (rules, the embryo diagnostic table, rows R1, E8, N10 and N11);
    - `experiments/research_2026_11/REPORT.md` (sections 1, 5-7) and its `notes/task-a.md` L3;
    - `experiments/research_2026_11/results/sensitivity.json`;
    - `docs/research/digest_2026-10-19.md` (the variogram and reference-row lines only).
- **Not opened:** boards/, data/, runs/, any other file under docs/research/, the leaderboard, the web, anything under
  ~/.claude-account or ~/.claude, and the other reviewers' files in this folder.
- **Git:** reads only, with `GIT_OPTIONAL_LOCKS=0`.
- **Computation:**
  - The estimator's loaders, run in memory as `python -I -B` from stdin; no file written.
  - The embryo critic's closed-form model, reproduced and extended to tail values.
  - No embryo target position is derived or quoted.
- **Barred classes:** none met.
- **Written:** this file only.
