# T2:heart:test_extrap: adversarial critique of the 2026-11-06 design

Board: T2:heart:test_extrap (E12.5 3D MERFISH, 500 genes plus coordinates; baseline copy_last = E10.5 with its
coordinates; counts 0.5 in Overall). Team Mithril, Human track. AS_OF 2026-11-06.

Written by the critic agent of the 2026-11-06 workflow, on `design/heart.md`. Citation keys follow the design: [task-x n]
is source n of the audited note `final/task-x.md`, and "task-x F n" its finding n. Team records are cited by repository
path. "The 2026-11-05 dynamics note" means `experiments/research_2026_11/notes/task-f.md` in the repository (audited and
committed), not this round's `final/task-f.md`.

---

## 1. Bottom line

- **No candidate survives as a try.** H-B survives only as a ruling that the captain already has (REPORT sheet E). It
  is worth about +0.03 Overall (range -0.25 to +0.25). The design's own summary adds up to about +0.08 Overall; the
  corrected figure is about +0.03.
- **H-A fails TRIES rule 5 on its own numbers.** It expects +0.1 board points on the backtests, against the 0.3 bar.
- **H-A also repeats a lever that the 2026-11-05 research closed by rules 5 and 2, and the design does not cite it.**
  That lever is L3 of the 2026-11-05 dynamics note: the same MA(1) moment calibration of the heart's momentum, with the
  same decile-matched permuted-gene placebo and the same two backtests.
- **The premise under H-A to H-D cannot be identified in these files.** The premise is that the part of B0's momentum
  shared across lineages is the reference run's offset. Neither D-H1 nor the backtests can identify it (section 5).
- **Recommendation: no heart try before the freeze.**
  - Leave sheet E at its default (no change) unless the captain wants a small bet that costs almost nothing.
  - If he does want it, the full hold is the cleaner form.
  - D-H1 can go with it as information only. It takes half a day, the captain signs its row first, and it has no
    "register H-A" branch.

## 2. Verdicts

| Candidate | Verdict | Corrected E[gain], Overall | Why, in one line |
|---|---|---|---|
| H-A lineage-contrast momentum, gated by D-H1 | KILL | 0 | Fails rule 5 by its own +0.1; repeats the closed L3 without answering rules 5 and 2; its premise is not identifiable; H-B beats it whatever D-H1 shows |
| H-B hold, or CM-only momentum, at test | RULING-ONLY | +0.03 (-0.25 to +0.25) | This is sheet E, still pending, with a decline recommendation; no counted check of the test form exists; the full hold is the clean form |
| H-C readout-structured removal | KILL | 0 | Fails rule 5 (+0.05); the codebook is unknown and co-dropout is unverified in our files; inherits H-A's prune |
| H-D run-noise shrinkage from the knockout files | KILL | 0 | Fails rule 5; needs an A15 reading; one run difference cannot estimate per-gene variance, and expression-level smoothing makes the weights rank-neutral |
| H-E fixed-RMS late organ shape | KILL | 0 (-0.5 to +0.5 if it could be built) | Needs A15 and block rulings, plus PLAN 4.3 step 5's shape gate and native z; no clean check; 1.5-2 weeks, so not buildable by 2026-11-27 |
| H-F organiser question on external spatial data | KILL (for this board) | 0 for the heart | No heart value; it would dilute the pending T3 severity question |

---

## 3. Candidate by candidate

The five checks: (1) repeat of a pruned or closed item, (2) rules, (3) counted evidence, (4) realistic gain,
(5) feasibility before the s1 freeze (2026-11-27 00:00 UTC).

### H-A. Lineage-contrast momentum, gated by D-H1: KILL (corrected gain 0)

**(1) A repeat, with the nearest item missing.** The nearest prior item is not one of the three the design names.
- It is "L3. Empirical-Bayes / Kalman calibration of the heart's pre-break velocity (closed by rules 5 and 2)", in the
  2026-11-05 dynamics note (lines 144-151). The design's read list includes "the heart lines of ...
  `experiments/research_2026_11/notes/task-f.md`".
- L3 and H-A share:
  - the fingerprint. "Cov(Delta_i, Delta_{i+1}) = phi s_v^2 Dt_i Dt_{i+1} - s_e^2" is "the MA(1) fingerprint of stage
    noise", and e there includes the run (lines 46-49).
  - the method: a target-free moment estimate from pre-break increments, paired with "a decile-matched
    permuted-gene placebo (as K2b)" (line 106). D-H1 is that recipe.
  - the two heart backtests.
- L3 was closed for two reasons. It gains under 0.3 on the backtests. And the backtests score the pre-break estimator,
  not the across-break one used at test (line 54).
- H-A adds one new element: the noise is assumed to sit in the part shared across lineages. That answers neither
  closing reason.

**(2) Rules.**
- **TRIES rule 5.** It says that "a change expected to gain less than the adoption rule's 0.3 points on its counting
  backtests is pruned before it runs (rule 3), since it could pass only by luck".
  - The design expects +0.1 board points on the backtests (range -0.5 to +0.8).
  - It puts P(adoption) at 0.1, or 0.15 if D-H1 passes.
  - D-H1 does not lift the expectation to 0.3, and the design never says it would.
  - So H-A is pruned before its run, whatever D-H1 shows.
- **Process.** D-H1's rule decides whether a counted try is registered. Under decision log 2026-11-04 (d), such a rule
  goes to the captain as a plain summary before the run. The design's first step ("send the captain the plain summary,
  and run it") has to wait for his sign-off.
- **No data rule is broken.**
  - Released heart MERFISH only, with E8.5 and E10.5 used under A20.
  - No A15, genotype, KO-ruling or post-exposure issue.
  - The frozen s2 follows s1 (`docs/S2_RULES.md`, "as s1").

**(3) Counted evidence. Even a pass would not show the mechanism.**
- **The backtests carry the run term but not the break term.**
  - At test the momentum is E9.5 -> E10.5. That increment crosses the dissection change (`configs/pseudo_tasks.yaml`,
    notes), so its shared part carries half of that sampling change.
  - Both backtests use pre-break momenta (E8.75 -> E9.5 and E8.5 -> E8.75).
  - "Exactly as it does at test" therefore holds for u_ref only.
- **The rule-2 defence fails for MERFISH.**
  - The design says that "a removal of the part shared across lineages cannot imitate a scope or composition change".
  - That ignores segmentation admixture. Neighbour-type molecules enter each cell in proportion to its neighbourhood
    mix, so neighbour markers dominate within-type differences when neighbourhoods differ [task-c F13; task-c 18].
  - A change of scope or composition therefore shifts every lineage's mean in the same per-gene direction.
  - The imaged field changes non-monotonically along the ladder:
    - E8.25_late includes yolk sac;
    - E8.75 is mostly heart;
    - E9.5 is 63% shell;
    - E10.5 is 77% surrounding tissue, and its core share fell from 37% to 23% (`docs/PLAN.md` 4.3; `docs/S2_RULES.md`;
      task-d F12).
  - That is a known feature of each backtest, and a shared-part removal can track it (TRIES rule 2, third condition).
  - The same alternation also produces D-H1's negative autocorrelation without any run effect.
- **Admixture also breaks "it keeps lineage-specific change".**
  - For abundant cardiomyocyte markers, spillover into the other lineages follows the cardiomyocytes' own change.
  - The across-lineage median then follows that change too.
  - So the contrast removes part of the one signal the design means to keep.
- **The depth drift persists.**
  - Counts per cell fall monotonically: 233, 166, 130, 75, 56 (TRIES H3).
  - Any depth-linked part of the per-gene run offsets is therefore a persistent shared trend. B0 extrapolates it with
    the right sign, and H-A removes it.
  - D-H1 sees this only if it dominates.
- **The design's own SOTA citation points the other way.**
  - DCI [task-a 5; task-a F5]: changes that agree across cell types are the ones that can be forecast, and changes
    specific to one type do not beat carry-forward.
  - H-A keeps exactly the type-specific part, for every tracked lineage.
  - The only heart evidence that a lineage-specific part carries over is the line "per-type momentum carries over only
    for cardiomyocytes" (`experiments/harness/TRIES.md` line 387). It has no committed source.
  - Its only committed analogue is T1's scRNA readout: CM 0.49, endocardium 0.44 and endothelium -0.53 (the
    2026-11-05 dynamics note, line 50). That is cross-modal for the heart.
- **The "signature of added run distance" is not specific.**
  - On the same 1-day backtest, B0's DE group gained +0.77 over the copy: de_score skill 0.510 against 0.464, and
    de_direction 0.521 against 0.505 (`experiments/harness/results/headroom/headroom.txt` lines 159-160).
  - That gain is hard to square with a dominant wrong-signed shared term.
  - B0's cell-level differences from the copy sit on the top-up's detection-loss line (TRIES, b0_smooth).
- **No released pair checks the test form.** No released stage gives a post-break increment followed by a post-break
  target. So the form H-A would use at test cannot be counted at all.

**(4) Realistic gain.**
- The heart sits at the floor. There, skill = d_f/(d_f + d) gives the least skill per unit of distance: a quarter of
  what it gives near the ceiling.
- A true +0.1 board effect passes the double bar (+0.3 and 8 of 10 seeds on both backtests) mostly by seed luck. Both
  backtests also share the E10.5 truth.
- **H-B beats H-A whatever D-H1 shows.**
  - If D-H1 finds the shared part transient, the same evidence says to drop the non-CM lineage-specific parts too: DCI,
    and the lineages that do not carry over. That is H-B.
  - If D-H1 finds it persistent, B0 stays.
- H-A's own form is therefore never the best option.
- As a ruling, it would fare the same way against H-B.

**(5) Feasibility.** It takes 1.5-2 days and fits before the freeze. Feasibility is not the problem.

**What survives:** D-H1 alone, as information for H-B (see H-B).

### H-B. Hold, or CM-only momentum, at test: RULING-ONLY (about +0.03 Overall, range -0.25 to +0.25)

**(1) A repeat.**
- This is REPORT sheet E (the heart's momentum hold, about +0.05, "Recommended: decline", default no change). It is
  also L2 of the 2026-11-05 dynamics note.
- The design calls it "declined at +0.05". It is in fact pending: no ruling on sheet E is on record (`docs/PLAN.md`
  decision log 2026-11-05 (f)).
- Its "new mechanism" is not new. The 2026-11-05 note already wrote the across-break increment as
  mu-dot + delta + e_10.5 - e_9.5, with the run inside e (lines 46 and 53), and sheet E was priced on that.
- The only new input would be D-H1's numbers. They describe how persistent pre-break changes are, not the across-break
  increment.

**(2) Rules.**
- It changes expression only and is legal. It can enter only by the captain's ruling, as R1 did.
- Its plain summary must say that it reverses a component the adoption rule kept: B0 beat the copy by that rule
  (+0.47, ahead in 9 of 10 seeds; +0.90, 10 of 10), although on the pre-break form.
- Under a hold, the heart's s1 predicts no expression change beyond the baseline. `METHOD_SUMMARY.md` would have to say
  so.
- **The s2 follows s1** ("as s1"; heart step 2 of `docs/S2_RULES.md`), so no s2 goes back to the captain.
  - The s2's payoffs were built with B0's momentum. That changed basis should be stated, but it does not trigger the
    rule.
- **Work it would need before the freeze:**
  - an updated golden table;
  - a new freeze-morning hash reference (today `runs/archive/builds_88b057f`);
  - fresh format checks.

**(3) Evidence.**
- Nothing tests the form used at test.
- The CM-only variant rests on the uncommitted carry-over line above.
  - If that line comes from T1's scRNA readout, the CM-only choice is arguably scRNA-derived and applied to MERFISH.
    A15's working answer keeps that kind of input out of uploads.
  - The full hold has neither problem.

**(4) Gain.**
- **What a hold removes at test:** both the run term and the break term.
- **What it also gives up:**
  - the top-up's detection-loss side effect. That effect carried B0's cell-state and neighbourhood gains on the
    1.75-day backtest, and depth falls at every released stage, so it may help again at E12.5.
  - any real cardiomyocyte maturation signal.
- **Net:** about +0.0 to +0.05 Overall, and nothing for our upper tail.

**(5) Feasibility.**
- Half a day, plus the rebuild and the checks.
- The ruling must come before the s1 freeze.

**The revision**, if the captain wants the bet at all:
- **Form.** Offer the full hold as the one form.
- **D-H1.** Attach it as information only:
  - drop its "register H-A" branch;
  - before the run, fix in its row a second gene set of broadly expressed genes, to blunt admixture. For example: genes
    whose highest lineage mean is within 3 times their median lineage mean at both stages of each transition, with the
    constant written in the row;
  - report the result as "transient or persistent shared change on pre-break transitions", never as "run offset
    measured";
  - the captain signs the row before it runs.
- **Cost and default.** It takes half a day and is not a try. The default stays as it is.

### H-C. Readout-structured run model: KILL (0)

- **(1) A repeat.**
  - Its technical-axis idea is task-c's lever 4. The task-c audit found that lever 4 partly repeats "a mean of runs"
    and the fluctuation-response prune, and that it competes for T3's one try.
  - On the heart, H-C is a sharper H-A, so it inherits H-A's prune: rule 5 and the closed L3.
- **(2) Rules.** MerQuaCo code needs a licence check and an `EXTERNAL_DATA.md` row. Otherwise it is legal.
- **(3) Evidence.**
  - The codebook is unknown.
  - Co-dropout blocks are unverified in our files. MerQuaCo's 31% of sections comes from MERSCOPE brain over 641
    sections [task-c F11], not from the released heart stacks.
  - Five heart runs cannot separate a readout basis from shared biology without the codebook.
  - Its pre-check tests whether the basis fits the shared part, not whether that part persists.
- **(4) Gain.** The design's own expectation is +0.05 board points on the backtests, under rule 5, with P(adoption) of
  at most 0.05.
- **(5) Feasibility.** It takes 3-5 days and would fit before the freeze, but the work buys nothing.

### H-D. Run-noise shrinkage from the knockout files: KILL (0)

- **(1) A repeat.**
  - Its shrinkage repeats "T3: a two-scenario factor exponent" and the T1 noise-aware direction (split halves found
    nothing to shrink).
  - Only the run-level variance is new.
- **(2) Rules.**
  - It needs the captain's reading of A15, whose row says "keep modalities separate per task".
  - It would add T3 files to the heart's dependence on A20.
- **(3) Evidence.**
  - One run difference (Mab21l2 against WT E9.5) gives one draw per gene.
  - Smoothing by expression level turns the weights into a function of expression level. de_direction partials out
    that level, and de_score's null is built on it, so the rank metrics barely move.
  - The knockout factors correlate only 0.15 (TRIES N8), so run factors do not repeat, and one run's variance does not
    describe the next run.
  - The matched WT E9.5 appears to be the heart's own E9.5 file: `docs/CHALLENGE.md` lists both at 164 MB. If so, the
    difference contains u_E9.5, which is the start of the test momentum itself.
  - The "three run differences" with Gata4 are two runs: TRIES N6 found the Gata4 replicates on one slide, M164.
- **(4) Gain.** 0 to +0.1 board points on the backtests, under rule 5.
- **(5) Feasibility.** 1-2 days, but only after a ruling that nobody has asked the captain for.

### H-E. Fixed-RMS core shape toward a legal late anchor: KILL (0; -0.5 to +0.5 if it could be built)

- **(1) An overlap the design misses.**
  - It overlaps the closed "heart's core-focused s2" and "Heart: no try", as the design says.
  - The design does not mention `docs/PLAN.md` 4.3, step 5. H-E is a shape operation under that step, which adopts one
    only if it improves both released transitions (E9.5 -> E10.5 and E8.25_late -> E8.5) by at least 0.15 d2 skill
    without lowering Dice. It is then applied at half magnitude in at most one upload.
  - The s2 is frozen as the trim, so that one upload would be s1. A step in s1 alone sends the s2 back to the captain.
- **(2) Rules.**
  - A15's working answer closes it on its plain words (task-d audit).
  - It needs a ruling that a fixed-RMS deformation lies outside the size and scope block.
  - z must stay native: "never rescale z across stages"; a 1.5x z mismatch costs about 12.6 board points
    (`docs/PLAN.md` section 3, item 6). So a 3-D organ deformation would have to stay in-plane.
  - It must use wild-type embryos only.
  - DMDD's data terms are unchecked, and an `EXTERNAL_DATA.md` row comes first.
- **(3) Evidence.** No clean counted check exists. The scored cloud's shape is set by the field of view and the stack
  (rule 2).
- **(4) Gain.** 0. The public record on shape changes (selected on validation scores) points down.
- **(5) Feasibility.** 1.5-2 weeks, plus two rulings and a register row. It cannot be built and checked by
  2026-11-27.

### H-F. Organiser question on external spatial and imaging data: KILL for this board (0)

- Its value to the heart runs through H-D and H-E, and both are dead.
- Visium and Slide-seq are sequencing-based, so they are cross-modal in any case [task-b M6].
- No legal mouse heart MERFISH exists at or after E14.0 (task-d, Gaps).
- The organiser question that matters now is the fifth email. Its answer on floored severity decides T3's s2, worth
  about +0.5 to +2.5 T3 under a "zero" answer (`docs/PLAN.md` decision log 2026-11-06 (b)).
- A further question on the same thread should not dilute that one. The organisers have answered no team rule question
  since 2026-10-16.
- If this question is ever asked, it should go after the T3 answer and not for the heart's sake.

---

## 4. The strongest candidates

1. **H-B as the captain's sheet-E ruling, in the full-hold form.** About +0.03 Overall (-0.25 to +0.25). The default
   stays "no change".
2. **D-H1, taken from H-A and stripped down to information for that ruling.** Half a day, not a try, with the row
   signed first. Its decision value is about +0.01.

Nothing else on this board survives.

## 5. The biggest risk in the designer's reasoning

The design's load-bearing premise is that **the part of B0's momentum shared across lineages is the reference run's
per-gene offset**. If so, removing it removes a noise term with the wrong sign. In these files, though, the shared
part also carries four other things:
- **Admixture of neighbouring types** [task-c F13; task-c 18]. It moves with every change of scope or composition, and
  the imaged field changes non-monotonically along the ladder:
  - E8.25_late includes yolk sac;
  - E8.75 is mostly heart;
  - E9.5 is 63% shell;
  - E10.5 is 77% surrounding tissue.
- **Depth-linked detection.** Depth falls monotonically, from 233 to 56 counts per cell, so this part persists.
- **Shared biology and the spillover of abundant markers.** DCI says this is the part that can be forecast
  [task-a 5].
- **The E10.5 dissection change, at test only.**

Three consequences follow:
- D-H1's negative autocorrelation can come from admixture alone, with no run effect.
- The backtests do not contain the dissection change.
- The design's rule-2 defence ("cannot imitate a scope or composition change") is false once admixture is allowed for.

So neither a D-H1 pass nor a backtest pass would show the mechanism. The same premise carries H-B's "new mechanism",
H-C and H-D. It also explains the process error: the design does not cite L3 of the 2026-11-05 dynamics note, which
already priced this calibration below rule 5.

## 6. Smaller corrections to the design's text

- **"Declined at +0.05" (sheet E).** Sheet E is pending, with a recommendation to decline; its default is no change.
- **"The heart, at the floor, is that case" [task-c 1, 2].**
  - The carriage bias grows when the truth sits near the reference, as with weak responses.
  - Our entry sitting near the floor says nothing about E12.5, which is two days, and possibly a change of scope, away
    from E10.5.
  - The carriage there is mostly the sampling noise of the scorer's 10% working copy.
- **"At E10.5 the core is 23% of cells, and a shell removal or core crop alone moves Dice to 0.86 or 0.79".** task-d
  F12 gives these Dice figures next to E9.5's 63% shell. It does not attribute them to E10.5.
- **[task-f 13]** (one held-out-time validation ranks models only moderately) says nothing about the sampling break.
  The break argument rests on the 2026-11-05 dynamics note, not on that source.
- **"H-A shrinks the factors".** This is true on average, because the contrast's energy is at most the factors' when
  the median is near the mean. It is not true for every lineage and gene. A predicted variogram loss written into the
  row in advance does not stop that loss from counting against the board total under the adoption rule.

## 7. Compliance of this critique

**Read.**
- `design/heart.md` and `brief.md`.
- The audited notes:
  - `final/task-c.md` and `final/task-d.md` in full;
  - `final/task-a.md`, `task-b.md`, `task-e.md` and `task-f.md` at their outlines and the sources the design cites;
  - `final/task-a2.md` at its header only.
- From the repository:
  - `AGENTS.md`, `docs/CHALLENGE.md`, `METHOD_SUMMARY.md` and `experiments/harness/TRIES.md` (all of it);
  - `docs/PLAN.md`: sections 1-3, 4.3 and 8, and the decision log from 2026-11-03 to 2026-11-06;
  - `docs/DATA_POLICY.md`: windows, Scores and the change notes to 2026-11-05 (k);
  - `docs/S2_RULES.md`: the general rules and the heart rule;
  - `ASSUMPTIONS.md`, `docs/ORGANIZER_QUESTIONS.md` and `experiments/research_2026_11/REPORT.md` (sections 1-7);
  - `results/upload_gate.json` and `experiments/gate/results/estimate_2026-11-05.json`.
- For checks:
  - `experiments/harness/results/headroom/headroom.txt` (the heart rows);
  - `configs/pseudo_tasks.yaml` (the heart tasks) and `configs/lineage_map.yaml` (the groups);
  - `src/vec/t1_scope.py` (`momentum_shift`) and `src/vec/harness.py` (`_b0_t2`);
  - the 2026-11-05 dynamics note (sections 1-2, 7 and levers L2-L3) and the heart lines of `task-g.md` in the same
    folder;
  - the Dice lines of `experiments/coords/REPORT.md`;
  - the first lines of `experiments/labels/results/E10.5__vs__E9.5.csv` (released and unlocked labels; no count is
    quoted).

**Not opened.**
- `boards/`, `data/`, `runs/`, `docs/research/`, `raw/` and `screened/`;
- anything under `~/.claude-account` or `~/.claude`;
- the web, any leaderboard page and any Yukon material.

Searches excluded `boards/` and `docs/research/`.

**Nothing else done.**
- No hidden-target property, per-board leaderboard figure or in-window finding is used.
- Nothing was computed beyond reading committed files. One read of a committed JSON used the system `python3 -I -B`;
  no repository code ran.
- No repository edit was made; this file is the only output.
- No git or DGX operation was run. The relayed user request (sync the DGX and this machine; no unpushed commits)
  belongs to the main session before the workflow, and this critique step touched neither.
