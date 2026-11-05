# T1:test (E12.5 scRNA): adversarial critique of the design of 2026-11-06

Written on 2026-11-06 by the AI assistant as the T1 critic of the state-of-the-art workflow. It critiques
`design/t1.md`. This is design review only: no run, no try, no repository edit, no web. T1 counts 1 in Overall, so
board points and Overall units are the same here. [task-x n] cites the audited notes in `final/`.

## 0. Bottom line

| Candidate | Verdict | Designer's Overall | Corrected Overall |
|---|---|---|---|
| A. Leonard direction, counted, then pooled with Feng | KILL | +0.1 (-0.3 to +0.5) | about 0 (-0.3 to +0.2) |
| B. Leonard nuclei drawn into category quotas (bold) | KILL | +0.3 (0 to +0.6) | about -0.2 (-1.0 to +0.3) |
| C. Legal blood and CM anchors priced on the grid: keep today's pair | KEEP | 0 | 0 (status quo) |
| D. Cycling-state quota anchored on Leonard | RULING-ONLY | +0.1 to +0.3 | about +0.05 (-1.0 to +0.4), by ruling only |
| E. Larger fibroblast-like share (births revisited) | KILL | -0.2 | -0.2 if built; 0 once killed |

The design's "+0.4 Overall, and +2 to +3 on T1's 95th percentile through B" does not survive. The realistic
contribution of this set is about +0.05 Overall in expectation, with no measurable lift of T1's upper tail.

The design's arithmetic is right; its premise is not. A and B treat "same chemistry" as "same scale". They then
count on E9.5 -> E10.5 to validate an anchor that shares that backtest's dominant known feature.

## 1. What I re-checked, and what holds

I re-weighted the committed per-seed grid files (`experiments/s2_hedge/results/per_seed/t1-grid/`, 15 truths x 15
predictions x 10 seeds), in memory with `python3 -I -B` and the standard library. Only score and skill means were
printed. All of the designer's grid figures reproduce:
- today's pair (bown_cown, bown_c0.6) gives 68.28 under the uniform prior;
- the best prediction for each truth gives 71.20, so exact knowledge of the mix is worth +2.92;
- that gain is +2.19 if blood only stays or falls, +1.58 under the nuclei-leaning CM prior, and +1.02 with CM known
  at E10.5's share;
- the best pair over all 105 pairs is today's under the first three priors;
- a low-blood s2 wins (+1.03) only when the CM share is known at its own value;
- the asymmetry is about 11 against 2: b0.05_cown scores 58.80 against the own-blood truth, where the matched
  prediction scores 69.92, while bown_cown scores 64.04 against the 5%-blood truth, where the matched prediction
  scores 66.07.

At a matched mix the points left are 8.4-9.0 (de_score), 6.2-7.0 (de_direction), 5.6-7.2 (mmd_u) and 6.7-7.2
(variogram) in the blood-own slice. The designer's ranges are a little narrow; this does not matter. The headroom
row, estimate figures, windows ((8.5, 11.5] and (7.5, 10.5]) and bracket fractions (1/5, 1/6, 1/2) are as stated.

## 2. Candidates

### A. A second direction from Leonard's nuclei: KILL

1. **It is not a direction. It is a cross-lab level difference.** The factor is
   log f = t x [log(m_Leonard + 0.5) - log(m_base + 0.5)], where m_base is our own base stage's category mean.
   - Leonard has one stage, so no within-Leonard time contrast cancels the lab's offset. Every difference enters
     log f in full:
     - CellBender ambient removal [task-b Findings 1];
     - pipeline and annotation;
     - C57BL/6 against the Chi lab's mice;
     - "Hearts were dissected from embryos", not a stated chamber cut [task-b audit, spot-check 1];
     - doublet handling.
   - At test, t = 1/2 applies half of that offset. Pooled with Feng it is still a quarter.
   - `docs/PLAN.md` 4.1 sets the standing rule this breaks: "Use within-dataset fold changes, never absolute external
     levels." The design does not mention it.
   - With apply_factors' keep-depth constant, a few large offset genes rescale every other gene. N2 showed that the
     constant alone moves thousands of genes past 10% (TRIES N2).
2. **Forced category calls.** The classifier is trained on the base stage's labels: E9.5 on the backtest, E10.5 at
   test. It must place every E14.5 class into one of those categories, including fibroblasts at about 26% pooled.
   The births scout met the same problem: forced argmax calls, mean probability 0.18 (TRIES, births prune).
   Category means at E14.5 then carry states the base category does not have, and A moves E10.5 cells halfway toward
   them.
3. **The counted evidence measures the wrong comparison, and it is confounded.**
   - Try 1 is Leonard against B0's momentum. The change that would be submitted is Feng pooled with Leonard against
     Feng alone. TRIES item 5 ("T1 direction: no other try") still governs that change in full, because no check
     exists for it.
   - The backtest's main known feature is embryo portion -> dissected heart. Any dissected-heart, ambient-cleaned
     anchor reproduces it.
   - One concrete form of it: E9.5's library holds non-heart tissue (hepatocytes 7%; `docs/CHALLENGE.md`). If our
     nuclei carry ambient RNA from that tissue (unknown), a CellBender-cleaned Leonard minus E9.5 predicts those
     genes falling. The E10.5 truth shows the same fall, for a reason absent at test.
   - So rule 2's third condition fails before the run, as it did for births ("No score pattern could tell such a gain
     from the scope and composition gaps the counted tasks already carry").
   - The pruned "ladder-consistency gate" was dropped for reading this same sampling change.
4. **The guard cannot do its job.** On E8.5 -> E9.5 at t = 1/6, Leonard's factor mixes four things:
   - development from E8.5 to E14.5;
   - dissection;
   - the lab offset;
   - that task's own scope changes (neural tube out, hepatocytes in).

   A loss above 0.3 cannot be pinned on dissection, and no loss clears nothing.
5. **Pooling adds little even if Leonard were clean.**
   - The screened digest, as cited and audited in `final/task-a2.md` (not read by the designer), already has an
     independent within-dataset long-range bracket that agrees with Feng: JAX, sign 0.82-0.91, Spearman 0.68-0.72.
     A second source mostly adds its own error.
   - [task-a 6]'s pooling gain came from 11 studies, each with its own time contrast, in a learned elastic net. Two
     sources, one of them with no time contrast, do not reproduce it.
6. **Cost.** It spends two counted tries on a transition read many times.
   - Decision log 2026-11-05 (b) spends no try on a pruned idea, "which would only raise the lucky-fit risk". A's
     counted reading is the one the ladder-consistency prune refused.
   - TRIES item 6 requires "a new reason that holds at test". The lab offset is the part that does not hold at test.

Corrected gain: about 0 (-0.3 to +0.2).

What survives is only the target-free intake (GSE272154 Theiler-code check, the register row, WT-only shares by
code), and only if C's record or D needs it. It carries no points by itself.

There is one legal direction-side use of Leonard: a hygiene mask from Feng-vs-Leonard E14.5 level ratios. That
mask would flag highly expressed cytoplasmic transcripts first, the genes that carry T1's down-tail on the released
transition (`METHOD_SUMMARY.md`, "What the T1 scorer rewards"). Its sign is unknown, so it is not recommended.

### B. Leonard's nuclei drawn into the category quotas: KILL

1. **The cited evidence is all same-dataset.** The embryo mixture draws two brackets from one lab and one pipeline.
   The obesity runner-up resampled residuals within its own experiment [task-f Findings 12; 21]. Mimyr's lookup copies from the
   same atlas [task-d 8]. None of them supports scoring another lab's nuclei in a 30-PC space fitted on the truth.
2. **The evidence logic runs backwards.** Write the prediction as (1 - w)Q + wR, where Q is the base cells, R is
   Leonard's cells and P is the truth. Then
   MMD^2 = ||a - w b||^2, with a = mu_P - mu_Q and b = mu_R - mu_Q.
   - The gain term, -2w<a,b>, is linear in w. On the backtest it is inflated by the shared portion -> heart component
     of a and b. At test, where both ends are dissected hearts, that component is gone.
   - The penalty term, w^2 ||b||^2, carries the whole lab offset. It grows (0.5 / 0.2)^2, about 6x, from the
     backtest's 1/5 to the test's 1/2.
   - So a pass is weak evidence for B and a fail is strong evidence against it, the reverse of the K2a-style reading
     the design claims.
3. **Depth.** The design itself expects the (f)(3) depth report to fire.
   - T1 depth is closed because the team has seen a depth proxy of the target (the file size with the count;
     DATA_POLICY 2026-11-04 (f)(3)). Asking the captain to admit a depth shift after that exposure loosens a
     safeguard (`AGENTS.md` rule 7).
   - Matching depth by thinning is the closed depth line: thinning lowered mmd_u in 3 of 4 cases (TRIES H3).
4. **Sparsity.** CellBender and another lab's depth change per-gene detection. At test that affects half of all
   non-blood cells, and the variogram on 32,285 genes is mostly a per-gene sparsity test (`docs/PLAN.md` 3.1).
5. **No hedge.** An s1 change reaches the s2, which is s1 plus 60% CM, so both uploads carry B.
   - Decision log 2026-11-04 (e) puts bolder bets in the backup slot, and that slot is frozen.
   - A bold bet that lands in both slots lowers the chance of a high total.
6. **Data terms.** The uploads would contain another lab's nuclei verbatim, while the GEO terms are unstated.

Corrected gain: about -0.2 (-1.0 to +0.3). P(pass) is perhaps 0.2-0.25, and much of any pass would be the confound.
At t = 1/2 the expected test effect is negative, plausibly -4 to +2. The claimed +2 to +3 on T1's 95th percentile
has no support.

### C. Keep today's pair, priced on the grid: KEEP

The arithmetic is verified (section 1). It costs hours, closes task-b L1 and L2 for this season, and protects the
60% s2's upside. Two edits are needed:
1. The "nuclei-leaning" CM prior leans on Leonard's pooled 29%. That figure includes the PMIS libraries and comes
   from a source not yet registered, so it is motivation only [task-b audit, L2 (4)]. The decision-log line should
   say so and name the WT-only share as the one a prior could use, computed by code after registration.
2. Do not run the optional 10%/15% blood extension unless the captain asks. It reads the counted E10.5 truth again
   (rule 4), and it can change nothing unless he is near-certain of the fall. Moving s1 to 5% blood breaks even at
   about 85% (11.1 against 2.0).

Corrected gain: 0.

### D. A cycling-state quota anchored on Leonard: RULING-ONLY

- **It cannot be counted.** At t = 1/5 the shift is a few points of cycling share within CMs. Paired seed noise has a
  median of 0.95 points (`ASSUMPTIONS.md` A17), so the effect cannot clear 0.3 points and 8 of 10 seeds. Rule 5
  prunes it as a try, so it enters only by ruling, as R1 did. The design's "it can be counted" claim gives no
  backtest estimate.
- **The sign it measures is depth-sensitive.** A cycling call on another lab's nuclei is a cross-lab level
  comparison too.
- **Revision.** Score both datasets with one call rule, after thinning both to a common depth. This is a diagnostic
  only, and no thinned cell enters a build.
- **Commit order.** E10.5 is both the test base and the counted truth, so the rule is committed before anything is
  computed on it (DATA_POLICY 2026-11-04 (f)(1)).
- **Gates, in order:**
  1. a check of the released label vocabulary, names only;
  2. if no cycling labels exist, the captain's A15 reading for T1 gene sets, an `EXTERNAL_DATA.md` row and PLAN 4.4
     provenance;
  3. the (f)(3) library-size report;
  4. the s2 back to the captain, since an s1 change reaches it.

Price: about +0.05 (-1.0 to +0.4). REPORT sheet D gives -0.3 to -1.5 if the sign is wrong, and Leonard's
measurement only partly settles the sign.

### E. A larger fibroblast-like share: KILL

- The births line is closed for this season, and the design admits that the prune's reasons (1) and (2) stand.
- The crosswalk it needs is not identified. Early VIC is a state inside the endocardial category at E10.5, and
  Leonard's "fibroblast" is a different class definition. A log-linear interpolation between them has no meaning.
- The 26% is pooled over PMIS libraries, so it is motivation only.
- The grid's under-/over-prediction asymmetry was measured on blood and CM, and it does not carry over to a new state.
- Like B, it lands in both uploads.

Corrected gain: -0.2 if built (the designer's own figure); 0 once killed.

## 3. The set-aside list and omissions

- I agree with every set-aside item.
- **Omitted.** REPORT sheet C, the strong-null down-tail prior, is the one skewed T1 ruling on record:
  - +1.5 to +2 T1 if right;
  - break-even at a prior of 0.2-0.3;
  - deadline: the s1 freeze.

  The design neither prices it nor checks how it interacts with its own candidates; A's pooling would reshuffle the
  same down-tail genes. Every T1 s1 change competes with it for the captain's s1 decisions and sends the s2 back to
  him.
- **Unread.** `final/task-a2.md` is a second, audited task-a note whose digest figures (A, reason 5) weaken A. Its
  lever 4, a counted JAX twin, is the design that does cancel lab offsets. It too only informs; it changes no entry.

## 4. Section 4, the path to 209

- **Compliance.** It uses the gate figure as a planning threshold: P(Overall >= 209), "17.5 points between reading
  A's mean and 209", and the gains T1 would need.
  - The REPORT's counter-review struck exactly this framing and restated distances against our own percentiles.
  - DATA_POLICY 2026-11-05 says the gate figure gates uploads only and is never a threshold for a try or a method
    choice.
  - Section 4 should be restated against our own percentiles.
- **Provenance.** It cites an untracked, unreviewed script (`experiments/gate/path_to_target.py`) and its output as
  evidence. Neither is committed evidence, and neither should be carried into a committed file.
- **Substance.** With A, B and E killed and D small, the T1 side adds about +0.05 Overall and no tail lift.

## 5. Strongest candidates and the biggest risk

**Strongest:**
1. C: keep today's pair and record the legal anchors as motivation. It is worth 0 Overall but costs hours, closes L1
   and L2, and keeps the 60% s2's upside.
2. D: a small positive, by the captain's ruling only, after the label check (or the A15 reading) and a depth-matched
   measurement.

**Biggest risk in the designer's reasoning.** It treats same chemistry as same scale:
- A's "direction" is our base stage subtracted from another lab's single stage;
- B's "real cells" are another lab's nuclei scored in the truth's PC space.

Both are to be validated on E9.5 -> E10.5. That transition's dominant known feature, embryo portion -> dissected
heart and the ambient change it brings, is exactly what any dissected-heart, CellBender-cleaned anchor reproduces.
A pass would be that confound. At test (heart -> heart) only the uncancelled lab offset remains, and for B its MMD
penalty grows about 6x.

## 6. Compliance and process (this critique)

- **Read:**
  - the brief, `design/t1.md` and all of `final/` (task-a to task-f, plus task-a2);
  - `docs/CHALLENGE.md`, `METHOD_SUMMARY.md` and `experiments/harness/TRIES.md` in full;
  - `docs/PLAN.md` sections 3, 4.1, 7 and 8 and the decision log of 2026-11-04 to 2026-11-06;
  - `docs/DATA_POLICY.md` (sections up to the do-not-use list, organiser rulings, change notes 2026-11-04 to
    2026-11-05 (e));
  - `docs/S2_RULES.md` (general part and T1), `ASSUMPTIONS.md`, `docs/ORGANIZER_QUESTIONS.md`;
  - `experiments/research_2026_11/REPORT.md` (summary, sections 1-2, 5, 7-9);
  - `results/upload_gate.json`, `experiments/gate/results/estimate_2026-11-05.json`;
  - `configs/pseudo_tasks.yaml` (the T1 tasks) and `experiments/s2_hedge/README.md`.
- **Computation:** in-memory re-weighting of the committed per-seed grid files (score and skill means only; no cell
  count printed). A helper script was briefly written to the scratchpad and deleted within a minute. Nothing else
  was written but this file.
- **Not opened:** `raw/`, `screened/`, `boards/`, `data/`, `runs/`, `docs/research/`, or anything under
  `~/.claude-account` or `~/.claude`. No web. The harness put the auto-memory index into this session's context; it
  was not opened or used.
- **Barred classes:** none met. No leaderboard figure is used here except to name the gate figure that section 4
  misuses.
- **The relayed user request** (sync the DGX and this machine; no unpushed commits on either) belongs to the
  workflow's owner, so this read-only step did not fetch, pull, push or touch the DGX. A read-only look at the local
  clone found main tracking origin/main, 0 ahead and 0 behind as of its last fetch (none was run), with two
  untracked files (`experiments/gate/path_to_target.py` and its result). The DGX was not checked.
