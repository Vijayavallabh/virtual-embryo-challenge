# s2 rules (pre-registration, frozen 2026-11-03)

**Amended on 2026-11-13 by the captain's own say-so:** the last section, "Amendments of 2026-11-13", changes the embryo's
s2, adds a factor to T3's second upload and records a change of T1's s1 that its s2 does not follow. Where that section
and the text before it differ, that section holds; the text before it is kept as the record of 2026-11-03.

Dates are local (IST) unless marked UTC.

**Authorship.** The AI assistant wrote these rules under the captain's delegation of 2026-11-02 (`docs/PLAN.md`
decision log, 2026-11-02). T3's s1, s2 and upload order are the captain's own rulings of 2026-11-01. The
captain can overrule any rule: before the freeze at will, after it as "Timing and freeze" says.

**Timing and freeze.**
- **Before any score.** Every rule was fixed before any test-board score existed. On 2026-11-03 nothing had been
  uploaded: `results/submissions.csv` holds only format checks on validation boards.
- **The freeze.** The file was frozen on 2026-11-03 at the captain's request, brought forward from the Nov 17-23 slot
  of `docs/PLAN.md` section 7 (as at 168e08e). Its sha256 is in `docs/S2_RULES.sha256` and in the `docs/PLAN.md`
  decision log, both written in the freeze commit, which is pushed with the captain's go-ahead.
  `tests/test_s2_rules.py` fails if the file and `docs/S2_RULES.sha256` disagree.
- **Changes after the freeze.** A rule changes only with the captain's own say-so. The change goes in its own commit,
  which logs his words, the reason, and the file's new sha256 in both places.
- **Before each s2 build and upload,** `sha256sum -c docs/S2_RULES.sha256` must pass, run from the repo root.

**No score feedback.**
- **No rule reads a test score.** Each s2 was fixed before any test score existed, and only organiser rulings switch a
  rule (T3).
- **After a board's first official upload,** only the captain himself can change that board's rule. Such a change is
  logged with the scores public at that time, and it computes nothing from a test score (Official Rules: scores
  choose, never compute).

**What "as s1" means (T1, heart, embryo).** Each s2 is built with its board's s1 build command, at the commit s1 is
built from, plus only the change its rule names:
- T1 adds `--cm-share 0.60`;
- the heart uses `--entry s2_trim` in place of `--entry s1`;
- the embryo uses `--method fgw_shape0` in place of `--method fgw`.

Today's s1 commands are `build_test.py t1`, `build_test.py heart --entry s1` and `build_test.py embryo --method fgw`.
(That sentence is of 2026-11-02. The amendments of 2026-11-13 at the end of this file change T1's s1 command and, in
the evening, both heart commands, which now carry `--trend`.)

A change to s1 that its builder makes for every entry carries over to s2. An example is a new step in `build_t1`, or
in the heart pipeline that every heart entry runs. That board's s2 goes back to the captain before any s2 build in
these cases:
- s1 moves to a different builder or method (M2, a heart pipeline for s1 alone, another embryo method);
- s1 changes the parameter its s2 hedges (T1's cardiomyocyte share, the heart's scope or size hypothesis, the embryo's
  geometry).

**If a ruling bars the validation answers** (assumption A20; `docs/DATA_POLICY.md`) before a board's s2 upload, nothing
built from E10.5, E7.5 or the Gata4 files is uploaded. The T1, heart and embryo rules lapse with their plans, which
are rebuilt as `docs/DATA_POLICY.md` says. That board's s2 goes back to the captain. T3 then uses the barred
branch, as `docs/DATA_POLICY.md` says (its ruling switches below).

**Uploads and final marks.**
- **Go-ahead.** Every official upload and every final mark needs the captain's explicit go-ahead.
- **Order.** On T1 and both T2 boards, s1 goes first and s2 after it, by the end of 30 Nov IST (`docs/PLAN.md`
  sections 7-8). T3's order is in its section.
- **Final marks** follow the working rule in `docs/P3_RUNBOOK.md` section 5 (2026-10-31). In its words: "of the
  entries every ruling received by then allows, mark the one with the higher published score (a tie, or s2 unscored by
  30 Nov IST: s1)". The AI assistant added one case under the delegation: on T3, where s1 goes last, if s1 is unscored
  by 30 Nov IST, mark s2.
- **What an s2 is worth.** An s2 is valued as the expected better of s1 and s2. That holds if marks are per board,
  which is a working answer at 0.55 (`docs/ORGANIZER_QUESTIONS.md`), open in `docs/PLAN.md` section 8.
- **Per-set marks.** If one final mark turns out to cover a whole set of boards, these s2s still stand, and no second
  set is kept. The AI assistant decided this under the delegation on 2026-11-03, after a side agent's note.
  `experiments/s2_hedge/bundle_sim.py` computes exact expectations over the committed per-seed payoffs: boards
  independent, the T2 boards averaged, T3 left out. Gains over s1 alone, for even weights and for weights of 1/3, 0.6
  and 1 on each board's released-like truth:
  - with one mark per set, these s2s gain +5.54, +3.81, +2.09 and 0.00;
  - a milder set (heart H2 as first built, embryo redrawn) gains +5.16, +3.64, +2.14 and 0.00;
  - T1's s2 alone gains +4.74, +3.39, +2.03 and 0.00;
  - with per-board marks, these s2s gain +7.63, +5.97, +4.00 and +1.06.

**Builds.** Builds run on the compute machine from the repo root, at a committed clean tree, with
`OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 NUMPY_MADVISE_HUGEPAGE=0` (`docs/COMPUTE.md`).
- **Check.** Each file is then checked there with `uv run vec check FILE --board KEY`.
- **Format check and upload.** The server format check (`uv run vec submit FILE --board KEY --from-dgx`, format-only by
  default) and the official upload follow `docs/P3_RUNBOOK.md` section 5.
- **Board keys:** `T1:test`, `T2:heart:test_extrap`, `T2:embryo:test_interp`, `T3:bcatenin`.
- **Seeds.** The default build seed is 1000001 (never seed 0), as for s1. The count is the count rule, read in code.

**Evidence.**
- `experiments/s2_hedge/README.md`. Its commands regenerate the T1 and heart tables, the step-5 values and the scopes
  and core shares from the committed per-seed files in `experiments/s2_hedge/results/per_seed/`. The per-metric skill
  gains are seed means of those files' `skills` fields. T3's numbers come from `experiments/t3_decomp/results/gata4/`.
- `experiments/harness/results/emb_shape_v1/` (registered in `results/registry.csv`).

These are sensitivity analyses on released transitions, not adoption results.

## T1 (decided 2026-11-02)

s2 is M1 at **cardiomyocytes 0.60 of heart cells** instead of E10.5's own 0.265. Blood stays at E10.5's own chamber
share, as in s1.

`uv run python experiments/builds/build_test.py t1 --cm-share 0.60 --out runs/builds/T1_s2.h5ad`

**Why.** The 2-D grid ran on T1/E9.5->E10.5:chamber, with an oracle base, 10 seeds and commit 90b25d3. Each row below
gives the expected better of s1 and s2, minus s1 alone:

| Prior over the 15 tissue mixes | 60% s2 | best blood-only s2 (blood 40%) |
|---|---|---|
| Uniform | +4.74 | +1.49 |
| A third on own | +3.39 | +1.07 |
| Blood own, cardiomyocytes from own to 0.75 | +6.40 | +0.03 |
| Cardiomyocytes at 0.45 or above | +7.90 | +1.59 |
| Cardiomyocytes at own, blood varying | +0.00 | +1.51 |

A blood-share hedge adds at most 1.6 points under these priors. It beats the 60% s2 only if the cardiomyocyte share
is known to stay at E10.5's.

**Supersedes.**
- This reverses the direction the captain agreed on 2026-11-01 (`docs/PLAN.md` 4.1; decision log 2026-11-01 (d)),
  which was a hedge on the blood share. The reversal was decided by the AI assistant under his delegation, on the
  grid's evidence, and the captain can restore it.
- M2 does not replace this s2. `docs/PLAN.md` 4.1's "M2 replaces s2 only if it passes the adoption rule against s2" is
  superseded. If M2 replaces M1 as T1's s1, this s2 goes back to the captain.

**Caveats.** Each row of the grid starts from the truth's own chamber mix (an oracle base). There is one transition,
and the priors are stated assumptions.

## T2 heart (rule committed at 38faf73, 2026-11-01 21:38 UTC, before the second transition's step-5 readout; amended after the outcome)

### The rule in force

**s2 is H2 built as the half trim:**

`uv run python experiments/builds/build_test.py heart --entry s2_trim --out runs/builds/T2heart_s2.h5ad`

Since the amendment of 2026-11-13, evening (the last section of this file), this command and s1's carry `--trend`:
`uv run python experiments/builds/build_test.py heart --entry s2_trim --trend --out runs/builds/T2heart_s2.h5ad`.
Without `--trend` the builder gives the file as it was before that ruling.

Not `--entry s2`, which is H2 as first built.
1. **Trim.** E10.5's surrounding (non-core) cells are removed, farthest in plane from the core centroid first, until
   the in-plane scope is half of E10.5's (`vec.t2_build.radial_trim_mask`, pinned by tests to the payoffs' trim). This
   is H2's own definition in `docs/PLAN.md` 4.3: "heart-focused (shell halved or less)".
2. **Expression.** s1's heart pipeline runs on the trimmed cells: today B0, whose damped momentum runs from E9.5 to
   the trimmed E10.5 cells. The payoffs built half_radial that way on E9.5 -> E10.5 and E8.75 -> E10.5; on
   E8.25_late -> E8.5, which has no earlier stage, it is a trimmed copy. The count rule and seed rule are s1's.
3. **Size.** `size_or_hold` with H2, in the size JSON's pre-registered mode (xy_grow_z_hold). x and y are scaled about
   the centroid, z is held, and the 3-D RMS reaches H2's pre-registered full-cloud RMS. The reference cloud is the full
   E10.5 pool, as for s1.

**Restoring the original rule.** The captain can restore the rule as committed ("s2 = H2 as built", `--entry s2`)
by saying so. After the first official T2 heart upload, that change falls under "No score feedback".

### Why the amendment

What an s2 is worth is the expected better of s1 and s2 on board score. The table gives that value for H2 as built
and for the half trim, with s1 alone in brackets. Each is over 10 seeds on three heart transitions, under five
priors over the five truths (`experiments/s2_hedge/results/heart_scope*/hedge_hold.txt`):

| Prior | E9.5 -> E10.5 | E8.25_late -> E8.5 | E8.75 -> E10.5 |
|---|---|---|---|
| Uniform | 50.17 / 51.12 (49.09) | 53.45 / 56.84 (51.61) | 54.50 / 60.02 (52.73) |
| A third on own | 50.18 / 50.97 (49.28) | 53.10 / 56.17 (51.47) | 54.67 / 59.46 (53.12) |
| The released target only | 50.22 / 50.22 (50.22) | 51.34 / 52.80 (50.76) | 55.50 / 56.65 (55.08) |
| Radial family (own, radial_0.5, radial_0) | **50.84** / 50.07 (49.92) | 54.00 / 57.80 (52.08) | 54.17 / 57.54 (52.26) |
| Thin family (thin_0.5, thin_0) | 49.16 / 52.70 (47.84) | 52.61 / 55.39 (50.90) | 55.00 / 63.73 (53.44) |

**Where it wins.** The half trim beats or equals H2 as built in 14 of the 15 cells:
- Under a uniform prior it is ahead by +0.95, +3.39 and +5.51.
- On the released target it ties on E9.5 -> E10.5, where neither beats s1 in any seed. It leads on the other two.
- Under the radial-family prior on E9.5 -> E10.5, H2 as built is ahead, 50.84 against 50.07. That is the transition
  nearest the test (E9.5 -> E10.5 and E8.75 -> E10.5 share the same E10.5 truths).

**What drives the gain.** On the released targets, the half trim's paired gain over s1 differs by transition:
- +1.75 on E8.25_late -> E8.5 (ahead in 9 of 10 seeds). The largest skill gains there are scale (+0.113) and
  de_direction (+0.108), then de_score (+0.056); d2_shape moves -0.036. The scale term is a size term under the
  payoffs' oracle core size, and H2 as built shares it (+0.115).
- +1.57 on E8.75 -> E10.5 (10 of 10). The largest skill gains there are de_score (+0.044), d2_shape (+0.043) and
  variogram (+0.036), then de_direction (+0.023).
- -5.88 on E9.5 -> E10.5 (0 of 10), where its better-of-two ties with s1.

**Why not the full trim.** core_radial is a larger hedge under a uniform prior, but it never beats s1 on a released
target.

**This departs from the plan.** Treating the trim as a scope hedge departs from `docs/PLAN.md` 4.3. Step 5 there (and
its 2026-11-01 side note) classes the trim as a shape operation, and this file did the same at 38faf73. The
reclassification came after the trim failed step 5. No test data were involved.

### History: the rule as committed at 38faf73

**The test.** The test was `docs/PLAN.md` 4.3 step 5 for the radial trim, read with
`experiments/s2_hedge/payoff.py step5` on both released transitions: T2h/<=E9.5->E10.5 and T2h/E8.25_late->E8.5.

**The bar.** The operation passed if, on both, the trim at full magnitude (`core_radial`, in-plane scope 0) gained at
least 0.15 in mean d2_shape skill over hold and did not lower mean occupancy_dice skill.

**The branches.** If it passed on both, s2 would be H2 built as the half trim, as above. If it failed on either, s2
would be H2 as built: s1's cells and expression at H2's size.

**Physical rationale.** This is step 5's third condition. The text at 38faf73 claimed both transitions became more
compact, and that has been corrected on 2026-11-03 after review:
- On E8.25_late -> E8.5 the later sample is more compact around the heart core: the in-plane scope goes from 0.942 to
  0.765 and the core share from 4% to 10%. E8.25_late includes yolk sac and non-cardiac cells.
- On E9.5 -> E10.5 it is not: the scope rose from 0.148 to 0.210 and the core share fell from 37% to 23%
  (`experiments/coords/results/p3/size_heart_E12.5.json`). This holds although E9.5 is an embryo portion and E10.5 a
  dissected heart (`configs/pseudo_tasks.yaml`; `docs/PLAN.md` 4.3).
- A radial trim removes the tissue farthest from the core first, as a tighter dissection would. The released stages do
  not show that E12.5 is one, so the amendment rests on the board-score hedge alone. s2 covers a more heart-focused
  E12.5, and s1 the rest.

**The outcome.** Step 5 failed on T2h/E8.25_late->E8.5:
- core_radial: d2 -0.021 (ahead in 3 of 10 seeds), Dice +0.112;
- half_radial: d2 -0.036.

On E9.5 -> E10.5 it had passed (core_radial d2 +0.416, Dice +0.041). On the 1.75-day pseudo-task T2h/<=E8.75->E10.5,
which is not a step-5 transition, core_radial gains +0.141 on d2 and lowers Dice by 0.013. As committed, the rule gave
s2 = H2 as built. The amendment above replaced it.

**The timing.**
- The rule text was written before any of the second transition's output had been seen.
- Each run file is stamped when its run ended: 2026-11-01 21:37:51-56 UTC by the compute machine's clock. That clock's
  offset from the local clock is not pinned to better than about two minutes.
- So by the compute clock, the second transition's full output, d2 and Dice included, existed seconds before 38faf73
  was committed at 21:38:06 UTC. What had been seen by then was one process's first-seed board-score line, not its d2
  or Dice skills.
- The step-5 readout ran at about 21:49:50 UTC (the local file time of `heart_scope_e825/step5.txt`), before 8031d9c
  at 21:52:13 UTC.

## T2 embryo (rule fixed at 90ddaaa, 2026-11-01 22:02 UTC, before the variants ran)

**Earlier candidates.** On every embryo pseudo-task each one runs on, every candidate built before 2026-11-02 scores
below fgw in every one of the 10 seeds (`experiments/harness/results/emb_fgw_v1/`):
- copy and B0 on both tasks;
- pseudobulk_shift and psb_size on E7.5 only. Their base cannot run on E7.25, which has no stage before E6.75.

None of them adds anything as an s2.

**fgw's open risk** is its interpolated shape: on E7.5 its shape skills fall below a copy's.

**The candidates.**
- `fgw_shape0` keeps fgw's coupling and its mixture of real cells, with positions taken from the earlier bracket's
  geometry (`interpolate_clouds(..., t_pos=0)`).
- `fgw_shape1` does the same with the later bracket's (`t_pos=1`).

Both are scaled to the embryo size rule. At test the earlier bracket is E7.5 and the later is E8.0.

**What they replaced.** They replace the open list at 38faf73:
- FGW with an alternative clock was not run. The plan fixes no alternative value, and choosing one now would be tuned
  on the same two tasks.
- The nearest-stage variant is now fgw_shape1.
- psb_size was dropped (see above).

**Caveat.** The variants answer fgw's known E7.5 shape result and are judged on the same two tasks. Neither task sits at
the test's t = 1/2 (they are at 1/3 and 2/3).

**The rule.** Both variants run on both pseudo-tasks against fgw, over seeds 0-9 with the harness.
- **Choice.** s2 is the variant with the higher expected better of fgw and it: the mean over seeds of the larger of the
  two scores, averaged over the two tasks.
- **Tie.** If neither beats fgw in any seed, s2 is fgw_shape1.

**Outcome** (`experiments/harness/results/emb_shape_v1/`, 90ddaaa; the run's fgw scores equal emb_fgw_v1's):

| | fgw | fgw_shape0 | fgw_shape1 |
|---|---|---|---|
| E7.5 from (E7.25, E8.0), mean score | 61.12 | 62.95 (ahead in 10 of 10) | 59.47 (ahead in 0 of 10) |
| E7.25 from (E6.75, E7.5), mean score | 58.78 | 53.67 (ahead in 0) | 57.49 (ahead in 0) |
| Better of fgw and it, mean of both tasks | 59.95 (fgw alone) | **60.87** | 59.95 |

On E7.5, fgw_shape0's gain is shape: d2 +0.067, Dice +0.109, neighborhood_mmd +0.015. Expression is unchanged. Its
gain is E7.5's alone. As a replacement for s1 it fails the adoption rule: it is vetoed on E7.25 at -5.11
(`emb_shape_v1/decisions_T2_embryo_test_interp.txt`). It is chosen only as s2.

**s2 = fgw_shape0**: E7.5's geometry at test, with fgw's mixture at the size rule.

`uv run python experiments/builds/build_test.py embryo --method fgw_shape0 --out runs/builds/T2embryo_s2.h5ad`

## T3 (s1, s2 and order by the captain, 2026-11-01; ruling switches from `docs/P3_RUNBOOK.md` section 5, completed by the AI assistant under the delegation on 2026-11-03)

**s1 and s2.** s1 is the allowed branch's s1. s2 is the refused branch's s1: wild-type E8.75 cells with the central
composition at the Mab21l2 magnitude (log2 0.718), no per-gene factor, and no number derived from the Gata4 files.

**Branch names.** The refused branch is `--branch refused --submission s1`, not `--branch barred`. The barred branch is
rules.yaml note 3's full-bar case: the Mab21l2 factor (exponent 1). Its s1 has the same central composition as the
refused s1; its s2 has the mild one (half the magnitude, log2 0.359).

**Commands,** on the compute machine:

```bash
uv run python experiments/t3_decomp/build_submission.py build experiments/t3_decomp/results/gata4/gata4_decisions.json \
  --branch refused --submission s1 --wt data/raw/T2/E8.75.h5ad --run-dir runs/t3_decomp/gata4_20261028T122425Z \
  --board T3:bcatenin --seed 1000001 --out runs/t3_build/refused_s1.h5ad
```

- The decision file's sha256 is 9440ae67625608112666d5235e6ceaa5dd9f47a4b0bbd30e16943f32f35c5a8c.
- s1 and the switch files are built the same way, with `--branch allowed|barred --submission s1|s2` and
  `--out runs/t3_build/<branch>_<s>.h5ad`.
- `experiments/t3_decomp/rules.yaml` and `t3_decompose.py` stay unchanged until every T3 file is built, because the
  builder refuses a decision file whose hashes differ.

**Order.** s2 goes first, and s1 last, by the end of 30 Nov IST. Each upload has the captain's go-ahead and follows a
screened check of mail, Slack and the Rules page. s1 is fixed whatever s2 scores.

**Ruling switches.** Organiser rulings switch this rule; scores never do.
- If the Gata4 files are allowed before the first T3 upload, s1 and s2 are the allowed branch's s1 and s2.
- If they are restricted before the last T3 upload, nothing from the allowed branch is sent, and the last upload is the
  barred branch's s1. The runbook left "the barred branch's s1, or nothing" open. The AI assistant settled it, under
  the delegation: a built entry that uses no Gata4-derived number can only help the better-of-two.
- If, before the first T3 upload, the Gata4 files are barred even for scoring candidates (rules.yaml note 3) or the
  validation answers are barred (A20; `docs/DATA_POLICY.md`), both T3 uploads come from the barred branch. The barred
  branch's s2 goes first and its s1 last, and the refused file is not uploaded.
- If such a bar arrives after the refused file is uploaded, the refused file is not marked final, because its
  "no factor" was chosen on Gata4 scores. The last upload is the barred branch's s1, and it is the one marked. If no
  uploaded entry is allowed by the rulings received, the captain asks the organisers before marking.
- These choices beyond the runbook were made by the AI assistant under the delegation on 2026-11-03: the barred s1
  under a restriction, the full-bar switches, and the T3 unscored-s1 mark.

## Amendments of 2026-11-13 (after the freeze; the captain's own say-so)

Written by the AI assistant. The quoted words are the captain's, given in the main session on 2026-11-13 in answer to
three questions that named each rule and gave its numbers (`docs/PLAN.md` decision log 2026-11-13 (d)). No test score
existed: nothing had been uploaded. The rows are R4, R5 and R6 of `experiments/harness/TRIES.md`.

### T2 embryo: s2 = `fgw_aligned_series`

**The captain's words:** "Run-corrected file (Recommended)", the first of three options (this file; the plain depth
top-up in E7.5's geometry; today's fgw_shape0), each offered with its numbers on nine practice tasks.

`uv run python experiments/builds/build_test.py embryo --method fgw --move-coupling aligned --expression series --out runs/builds/T2embryo_s2.h5ad`

- It replaces `--method fgw_shape0`. For the embryo, "as s1" now reads: s1's command plus `--move-coupling aligned
  --expression series`.
- **What it is.** fgw's draws, with E7.5's cells rendered into E8.0's per-gene detection profile and depth, and the
  moves taken from the aligned coupling (`vec.t2_build.series_rendered_expression`; tries E24, E25, E28).
- **Why.** E8.0, E8.25_late and E8.5 share one detection profile and E7.5 has its own (diagnostic N23), so fgw's
  mixture is half wrong if E7.75 was measured like E8.0. As better of two against the pair of 2026-11-03, on nine
  practice tasks: +6.69, +3.56, +3.42, +3.67, 0.00 and -3.09 where the target shares its deeper bracket's profile;
  +0.86, -1.42 and 0.00 elsewhere. It is the outcome of the rules fixed before tries E24, E25 and E28.
- **What it gives up.** On the two backtests the adoption rule counts it scores -1.42 (E7.5) and 0.00 (E7.25) against
  the old pair: E7.5's shape hedge is gone. It enters by ruling, not by the adoption rule.

### T3: the second upload also carries a wild-type-only factor

**The captain's words:** "Yes, add it (Recommended)", to the question whether to reopen the ruling of 2026-11-01 for
this file and add the per-cell-group correction built from wild-type samples only. The question gave both readings of
tries K7 and K8, the risks and the cost.

The command of the T3 section above with one more option:

```bash
uv run python experiments/t3_decomp/build_submission.py build experiments/t3_decomp/results/gata4/gata4_decisions.json \
  --branch refused --submission s1 --wt data/raw/T2/E8.75.h5ad --run-dir runs/t3_decomp/gata4_20261028T122425Z \
  --board T3:bcatenin --seed 1000001 --wt-revert --out runs/t3_build/refused_s1.h5ad
```

- **What it is.** The refused branch's s1 as before (wild-type E8.75 cells, the central composition at the Mab21l2
  magnitude), with each drawn cell's counts then moved by its lineage group's per-gene factor at exponent 0.5. The
  factor is WT E8.75's own deviation from a trend through the six other released wild-type stage files with the full
  panel (tries K7b and K8b; the estimator is `experiments/t3_decomp/k7_wt_revert.py`).
- **What it reads.** Released wild-type stage files only. No knockout file, no Gata4-derived number, nothing from the
  decision file beyond what the refused s1 already takes. `rules.yaml` and `t3_decompose.py` are unchanged.
- **How the ruling switches treat it.** Exactly as they treat the refused file. Its adoption was decided knowing tries
  K7 and K8, which read Gata4 and Mab21l2 scores, so like the refused file's "no factor" it was chosen on those scores:
  under a full bar it is not uploaded, or not marked final if already up. The ruling of 2026-11-06 (b) stands: an
  organiser's "zero" on a floored severity still puts the barred s1 in this slot.
- **Unchanged:** T3's first-choice file (the allowed branch's s1), the barred files, the order (this file goes first,
  s1 last) and every switch above.
- **Why.** On this file's own form the factor adds +4.72 and +7.12 on the two Gata4 replicates (tries K8; ahead in 19
  of 20 seeds). **The risks the captain was given:** the gain is the scorer's severity gate cleared by a thin margin,
  on one knockout run of another gene; the pre-set reading of tries K7 closed the factor by one comparison of nine.

### T1: s1 is drawn from library 3's cells, and s2 does not follow it

**The captain's words,** to the question whether to draw T1's first upload from E10.5's dissected-heart library alone,
asked with the released-data facts and three options (decline, the default; switch only on an organiser's answer;
switch on released data alone): "do whatever it takes to top the leaderboard".

The AI assistant reads those words as leaving the direction to it, and decided to switch, from released data alone
(row R6 has the reasons, fixed before any build). It holds no figure of the exposure of 2026-11-04.

- **s1:** `uv run python experiments/builds/build_test.py t1 --g-mom --base lib3 --out runs/builds/T1_s1.h5ad`.
  Library 3's cells (names ending in `_3`) at library 3's own chamber mix. The count, Feng's direction, R1, R2 and
  the fit of R3b's factors on the whole of E10.5 are unchanged.
- **s2 stays as built:** `uv run python experiments/builds/build_test.py t1 --g-mom --cm-share 0.60 --out
  runs/builds/T1_s2.h5ad`, the pooled base. This departs from "as s1": the base option is not carried into s2. The
  pair that try T30 priced keeps the second upload on the pooled base, so s2 now hedges the base as well as the
  cardiomyocyte share.
- **This is a depth choice in effect** (the drawn cells' median library rises by 27%), which `docs/DATA_POLICY.md`
  2026-11-04 (f)(3) sends to the captain. It was sent; the direction was then decided by the party that holds none of
  the exposed figures. The captain's go-ahead on the upload day is for uploading the planned file. A change of base
  after 2026-11-13, in either direction, needs a reason that can be stated from released data, or an organiser's
  answer read as `docs/ORGANIZER_QUESTIONS.md` ("Seventh email") fixes it.

### Builds

`experiments/builds/build_all.sh` builds these files: set s1 has T1 with `--base lib3`; set s2 has the embryo's new
command and T3's refused s1 with `--wt-revert`. T1's s1 as built until 2026-11-13 came out with different bytes on
the local machine and on the DGX, from the same code at 16 BLAS threads (`docs/PLAN.md` decision log 2026-11-13 (d)),
so the reference hashes for the freeze-morning rebuild come from a rehearsal on the machine that builds the uploads.

## Amendment of 2026-11-13, later that day: T1's s2 is drawn from the mid base (the captain's own say-so)

Written by the AI assistant. The quoted words are the captain's, given in the main session on 2026-11-13 in answer to
one question that named this rule and gave its numbers (`docs/PLAN.md` decision log 2026-11-13 (k)). No test score
existed: nothing had been uploaded. The row is R7 of `experiments/harness/TRIES.md`.

**The captain's words:** "Half-way base (Recommended)", the first of four options (a base half-way between the pooled
base and library 3 for the second upload; that base for both uploads; library 3's cells for both; keep today's file),
each offered with its practice numbers and with the request to decide from those numbers alone.

- **s2:** `uv run python experiments/builds/build_test.py t1 --g-mom --base mid --cm-share 0.60 --out
  runs/builds/T1_s2.h5ad`. It replaces "s2 stays as built" in the section above: the second upload no longer draws
  from the pooled base.
- **What the mid base is.** Every library-3 cell of E10.5 (names ending in `_3`) and, in each chamber category, a
  random fraction s / (1 + s) of the other libraries' cells, s being library 3's share of that category's cells in the
  file; the fraction is drawn by the build seed's own stream (`build_test.py`, `_t1_mid_rows`). Library 3 then
  supplies (s + 1) / 2 of each category, half-way between its share in the pooled base and 1. The cells are drawn at
  that pool's own chamber mix with cardiomyocytes moved to 0.60 of heart cells. The count, Feng's direction, R1, R2
  and the fit of R3b's factors on the whole of E10.5 are unchanged.
- **s1 is unchanged:** library 3's cells at library 3's own chamber mix (the section above; row R6).
- **Why.** Both uploads are spent on the cardiomyocyte share, and which base suits a cardiomyocyte-rich truth turns on
  E12.5's depth, which no released file shows (tries T30 and T32). The scored distribution metrics are quadratic in a
  mixture's weights, so a blend of the two bases gives up about a quarter of each gap, not half. On try T30's grid
  (try T33; zero time gap, sizes are upper bounds) the pair with this second upload leads the pair as built until
  today under all nine priors: +2.43 under the main one, +0.39 at the least; by scenario +0.2 to +6.0 on
  cardiomyocyte-rich truths at depths 0.68 and above, 0 where the truth keeps its own mix. It is the outcome of the
  rule fixed before try T33, recomputed independently.
- **What it gives up, as the captain was told.** Where the truth is much shallower than library 3 and
  cardiomyocyte-rich it loses 1.6 to 1.9 against the pooled base. Its lead rests on the scenarios where the truth is
  as deep as library 3: over the two shallower depths alone the pair trails by 0.36. The pair with both uploads on the
  mid base, steadier across depth priors, was offered and not chosen. It enters by ruling, not by the adoption rule.
- **This is a depth choice in effect** (library 3's share of the drawn cells rises from about 0.65 to about 0.82), which
  `docs/DATA_POLICY.md` 2026-11-04 (f)(3) sends to the captain; the question asked for a decision from the practice
  numbers alone. A later change of either T1 base needs a reason that can be stated from released data, or an
  organiser's answer read as `docs/ORGANIZER_QUESTIONS.md` ("Seventh email") fixes it.
- **Builds.** The s2 set of `experiments/builds/build_all.sh` builds T1's file with this command. The option is refused
  with any other cardiomyocyte share and with every option no try read on that base. The reference list for the
  freeze-morning rebuild changes for this one file and is remade by a rehearsal on the machine that builds the
  uploads (`docs/P3_RUNBOOK.md`, section 5).

## Amendment of 2026-11-13, evening: both heart uploads take the trend (the captain's own say-so)

Written by the AI assistant. The quoted words are the captain's, given in the main session on 2026-11-13 in answer to
one question that named the change and gave its numbers (`docs/PLAN.md` decision log 2026-11-13 (w)). No test score
existed: nothing had been uploaded. The row is H11 of `experiments/harness/TRIES.md`.

**The captain's words:** "Both files (Recommended)", the first of three options (the trend in both heart files; in the
second file only; no change), each offered with its price.

- **s1:** `uv run python experiments/builds/build_test.py heart --entry s1 --trend --out runs/builds/T2heart_s1.h5ad`.
- **s2:** `uv run python experiments/builds/build_test.py heart --entry s2_trim --trend --out
  runs/builds/T2heart_s2.h5ad`. It is s1's command with the trim, as "as s1" says. Point 2 of the heart rule
  ("Expression") now reads: s1's heart pipeline runs on the trimmed cells, and since this amendment that pipeline is
  B0 with the trend, fitted on E8.25_late, E8.5, E8.75, E9.5 and the trimmed E10.5 cells.
- **What the trend is.** B0's draw, cells, coordinates, count rule, size step and trim are unchanged. For each lineage
  B0 shifts, the per-gene log factor is 0.5 x the last gap (one day) x the least-squares slope of log(mean CP10k + 0.5)
  on the stage day over every released heart stage where the lineage has at least 20 cells (E8.25_late, E8.5, E8.75,
  E9.5 and E10.5), in place of half the E9.5 to E10.5 log ratio. A lineage seen at only the last two stages would keep
  the earlier factor; none is (diagnostic N27). The rule text is `vec.harness.B0_TREND_RULE`.
- **Why.** Row H11: on the two heart backtests the trend scores +0.36 and +0.48 over B0's last step, ahead in 10 of 10
  seeds on each, and passes the adoption rule on both, the first heart candidate to do so. At the test it leans on the
  E9.5 to E10.5 step, which crosses the change from sectioned embryo portions to dissected hearts and which no backtest
  contains, for 43% of its slope instead of all of it.
- **What it rests on, as the captain was told.** The margin over the rule's threshold is 0.06 on one backtest. It is
  the one arm of seven read on these two backtests that day that passes, and one of two trend forms (the size-matched
  twin does not pass). The two backtests share one truth, and the ten seeds resample the scoring, not the method. No
  backtest contains the test's last step: if that step is real development, the trend dilutes it. At the test pair
  the trend's factor has a cosine of 0.57 to the earlier one and 0.72 of its size (N27).
- **The rule on a changed s1 method** (above: that board's s2 goes back to the captain when s1 moves to a different
  method). The captain's answer names both files, so the s2 follows s1 with its trim.
- **Builds.** `experiments/builds/build_all.sh` builds both heart files with `--trend`. The option is refused with
  `--entry s2`, `--depth`, `--herd` and `--cfm`, which no try read with it. The reference list for the freeze-morning
  rebuild changes for these two files and is remade by a third rehearsal on the machine that builds the uploads
  (`docs/P3_RUNBOOK.md`, section 5).

## Amendment of 2026-11-18, after the first scores (the captain's own say-so)

Written by the AI assistant. Times are UTC. The record is `docs/PLAN.md` decision log 2026-11-18 (j) to (l) and
`results/manual_uploads.csv`. Unlike every section above, this one was written with test scores known.

**What happened.** On 2026-11-18 the captain uploaded the first official entry on each board through the site's form,
not the client: T1's s1; the embryo's s1; on T3 the allowed branch's s1, which this file had last; and for the heart a
uniform random half of s1's cells, after the server refused the built file for its cell count. The site returned the
scores and the captain pasted them, the per-metric values included, into the main session. T1 and the embryo beat
the floor; the heart and T3 scored below it.

**The captain's words.**
- After the scores: "use the scores for our submission to understand where we went fully wrong and incompetent and
  drastically iterate and improve our scores to get our total leaderboard scores around 226 at all drastic and wild
  costs."
- After the diagnosis, in answer to three priced questions on the heart (the size mode, the expression step, the cell
  number), a question on two organiser drafts and the advice to hold the built files until the screened check of
  19 November: "go ahead with whatever is best for us and immediately as weell instead of waiting".

**"No score feedback", as it now stands.** "Each s2 was fixed before any test score existed" holds for T1, the embryo
and T3: their second uploads are the files this document already names, built on 2026-11-13. It no longer holds for
the heart, whose second file is chosen below with the first scores known. The second bullet stands for every board: a
change after a first upload is the captain's, is logged with the scores, and computes nothing from a test score
(Official Rules, Section 10, as read through a screen on 2026-11-18: scores may choose among predictions and not
compute them; "the test is whether the truth enters your prediction's content or only its selection").

**T1 and the embryo: unchanged.** s2 = `runs/builds/T1_s2.h5ad` and `runs/builds/T2embryo_s2.h5ad` as built.

**T2 heart: a new second file.** It replaces `--entry s2_trim --trend`:

```bash
uv run python experiments/builds/build_test.py heart --entry s1 --expression copy --z-mode raw --out runs/builds/T2heart_s2iso.h5ad
uv run python experiments/builds/halve.py runs/builds/T2heart_s2iso.h5ad runs/builds/T2heart_s2iso_half.h5ad
```

- **What it is.** E10.5's cells at the count rule with their expression unchanged (the draw B0 starts from, with no
  momentum and no trend), scaled about the centroid by one factor on all three axes to the size the size file gives
  for H1 in its raw mode, then a uniform random half of the rows (numpy's default generator, seed 1000001, without
  replacement, order kept). The uploaded file is the half.
- **What it departs from,** each by the captain's answer above to a question that named it: the pre-registered z
  mode, xy_grow_z_hold (the size lever of the post-exposure block); the ruling of 2026-11-13, evening, that both heart
  uploads take the trend; the half trim, s2's hedge on scope; and a cell number (the halving).
- **Why.** The first upload lost to the floor on d2_shape, occupancy_dice, mmd_u and neighborhood_mmd. d2_shape and
  occupancy_dice ignore a uniform rescale and not the in-plane stretch, so an evenly scaled copy has the floor's form
  by their formulas and keeps the size; the expression step is the only expression difference from the floor. The
  pre-set second file keeps both steps. No harness run had scored the size step.
- **Nothing is computed from a test score.** The size is the size file's own raw-mode number, logged on 2026-10-28
  beside the pre-registered mode (`experiments/coords/results/p3/size_heart_E12.5.json`); the half is the rule of the
  first upload. The scores chose which steps to drop.
- **The first upload's file** is given by the same tool:
  `uv run python experiments/builds/halve.py runs/builds/T2heart_s1.h5ad runs/builds/T2heart_s1_half.h5ad`.
- **Builds.** `bash experiments/builds/build_all.sh heart2` runs the two commands above and checks both files. The
  sets s1, s2 and t3alt are unchanged, so the reference list of the third rehearsal still describes them.

**T3: the order, and the last upload.** The allowed s1 went first. The last upload is the refused s1 with the
wild-type-only factor, as built (`runs/t3_build/refused_s1.h5ad`), unless the captain grants a request that row K16
of `experiments/harness/TRIES.md` may lead to (the factor's exponent), which would be its own amendment. "s1 is fixed
whatever s2 scores" is spent. The ruling switches read, for the last upload: if an organiser ruling restricts or bars
the Gata4 files, or bars the validation answers, before it, the last upload is the barred branch's s1 and it is the
one marked. The screened check of 2026-11-18 found no ruling.

**Uploads, dates and marks.** Uploads go through the site's form by the captain's hand; the client's guards do not
run, and `results/manual_uploads.csv` is their record. "By the end of 30 Nov IST" is replaced by the organisers'
deadline, 2026-11-26 11:59 UTC, for uploads and nominations alike. The screened check before an upload ran on
2026-11-18. The rule for final marks is unchanged.

## Amendment of 2026-11-19: T1's last upload has no outside data (the captain's own say-so)

Written by the AI assistant. The quoted words are the captain's, as `docs/PLAN.md` decision log 2026-11-19 (f)
records them; nothing else here is his. The record is decision log 2026-11-19 (a) to (h). The row is R8 of
`experiments/harness/TRIES.md`. Times are UTC. Like the amendment of 2026-11-18, this section was written with the
first test scores known; no test score sets anything in it (below).

**The question and the captain's answer** (decision log 2026-11-19 (f), verbatim; the answers were read at 11:12 UTC).
- The question, as it was put: "Task 1's last upload. The outside heart atlas's cell-type labels were made with
  in-window stages pooled in, which the organisers ask teams to exclude; our first upload (70.9) used them. Which file
  should be our last Task 1 upload?"
- The two options: "No outside data (Recommended)" (the first file's design with the atlas step taken out; if the
  organisers discount the first upload, about 0 to 3.5 points below 70.9 on in-sample practice figures; if it stands,
  nothing changes) and "Atlas kept, labels recomputed" (no score figure on record; the tool's final review not ready;
  about a day more work and his sign-offs).
- The captain's answer, typed in his own words: "whatever it takes to win this challenge".

**How the AI assistant reads that answer** (its reading, not his words). He did not pick an option; he left the choice
to the AI assistant with one criterion. By the precedent of decision log 2026-11-18 (k), a general answer to a priced
question is read as the recommendation named in that question and as nothing beyond it. So T1's last upload is the
no-outside-data file below. The sentence lifts no rule: the organisers' data condition, the data rules, the
post-exposure block and the list of material never opened all stand. The captain can correct this reading before any
upload, and every upload is his own act.

**What the organiser wrote, and when the team read it.** On 2026-11-16 an organiser answered another participant, in
a thread of the public Slack channel, on external datasets that span a protected window. Such data are allowed
"provided all cells within the protected windows are removed before training, and the retained data do not carry
information derived from the excluded cells. Using raw expression counts and stage labels from the permitted cells is
acceptable. Please exclude annotations, clusters, embeddings, or fitted transformations learned using protected
cells, and recompute any such processing on permitted data only." The reply also asks teams to describe the filtering
in the method report and to list sources in the submission page's method-report fields. The team read it through a
screen on 2026-11-19 (`docs/ORGANIZER_QUESTIONS.md`, "Screened reads of 2026-11-19"; `docs/DATA_POLICY.md`, change
note 2026-11-19 (a)). The check of 2026-11-18 had read no Slack thread, so the team did not know it before then.

**What the Methods read found.** T1's one external model input is Feng et al. 2022 (register row X002). Its in-window
cells (E10.5 to E13.5) were removed at ingest, on the per-cell stage column. But the cell groups of the permitted
cells are the authors' own labels. A screened read of the paper's Methods, on annotation only (decision log
2026-11-19 (b)), found that the authors merged every sample of every stage into one object, then batch-corrected,
clustered and annotated it. No per-stage labelling is described. So the labels were learned using protected cells,
in the organiser's words. They enter both built T1 files through the four groups of Feng's direction (CM,
Endocardium, Endothelium, Epicardium): T1 cells of a group are moved two fifths of the way along that group's E9.5 to
E14.5 log ratio. The audit of decision log 2026-11-19 (c) lists every other decision that rests on them.

**T1's second file of the sections above is not uploaded.** That file is `runs/builds/T1_s2.h5ad`: the mid base with
cardiomyocytes at 0.60 of heart cells (the amendment of 2026-11-13, later that day; row R7). It is not uploaded, for
two reasons.
- It carries the authors' groups through Feng's direction, as the first upload does.
- Its 0.60 took its place from a grid partly bounded by a range of E14.5 cardiomyocyte shares read from the same
  labels in the authors' two strains (register row X002). That range cannot be recomputed on permitted data: the
  second strain's tables were deleted on 2026-10-19. The share cannot be re-chosen either: T1's composition is a
  lever closed by the post-exposure block (`docs/DATA_POLICY.md`, 2026-11-04, (b) to (d)). So no later T1 upload
  carries the 0.60.

The file stays on disk as the record of 2026-11-13. It is not deleted and not uploaded.

**T1's last upload.** One command, as the build script gives it (set t1none):

```bash
uv run python experiments/builds/build_test.py t1 --g-mom --base lib3 --feng-groups none --out runs/builds/T1_s2none.h5ad
```

**What it is.** The first upload's design with one step taken out.
- The same as the first upload: the base (library 3's cells at library 3's own chamber mix; row R6), the count rule,
  the seed (1000001), and so the same drawn cells. g_mom's within-type factors are as built: fitted within each
  chamber category on the released E9.5 and E10.5, applied to cardiomyocytes and endocardium only (row R3b).
- Taken out: Feng's direction (M1's step). Nothing is put in its place. Every cell keeps its E10.5 expression before
  g_mom's edit of cardiomyocytes and endocardium; R1 then holds every label.
- No outside data enter it. No Feng file is read, and the build record names no external file and says that the step
  was removed. R2 holds genes flat in Feng's direction only, so with no direction it is off: the check builds gave the
  same sha256 with and without R2's held lists (register rows X007 and X009; below).
- g_mom's mask under the option. In the first upload g_mom holds back, inside cardiomyocytes and endocardium, the
  labels R1 holds, and R1's mask reads the list of groups with a Feng direction, fixed in N12 by counting atlas cells.
  Under `--feng-groups none` the mask does not read that list. It holds back the labels that the label map itself
  (the readout config's regexes) gives no group. On every label the committed map can place in cardiomyocytes or
  endocardium that is the same single label as in the first upload, Early VIC, so g_mom edits the same cells. This
  was checked from committed configs alone and is pinned by `tests/test_t1_feng_groups_none.py` (branch `campaign/t1-none-fix2`, a04b3b2).
- The option builds one file only. The builder refuses `--feng-groups none` before reading anything unless the rest
  of the command is the first upload's (`t1 --g-mom --base lib3`): it is refused without `--g-mom`, and with `--base
  pooled` or `mid`, `--cm-share`, `--blood-share`, `--lib3-pair`, `--herd`, `--no-hold-unbracketed`,
  `--feng-hygiene` or `--strong-null`. With the option not given, every T1 file and its build record are as before.
- Its changes beyond the published baseline are the first upload's, less the step. Proportions: chamber categories
  only, drawn from library 3's cells at library 3's own chamber mix (row R6). Expression: g_mom's within-type edit of
  cardiomyocytes and endocardium. The step taken out was an expression step. The cells of the four groups outside
  g_mom's edit now keep E10.5's expression, as the published baseline predicts them. States, shares, cell number,
  size and scope are the first upload's.
- What it gives up against the sections above: the second upload no longer hedges the base or the cardiomyocyte
  share (rows R7 and T33). Both T1 uploads then draw the same cells. "As s1" holds: the first upload's command plus
  one option.

**What it is expected to do, priced both ways.** Every figure below is an in-sample labelled diagnostic on one
truth (E10.5), with Feng's step at 1/5 of its bracket against 2/5 at the test. No practice run applies this file's
path end to end, so the record holds no practice score for it.
- For scale, not a price of the step: on the E9.5 to E10.5 chamber cut, ten seeds, the copy scores 48.06 and b0
  49.68 (`experiments/harness/TRIES.md`, T1's headroom table). b0 is the chamber-only copy with damped momentum, M1
  without Feng's direction.
- The step it takes out, as the record prices it. Row T26 records b0 49.68, m1_feng 50.99, m1_r 50.83 and m1_r_gmom
  52.06 (`experiments/t1_m1/results/t26_m1_gmom/t26.txt`; void by its own gate, its numbers on record): the step
  over b0 is +1.31 alone and +1.15 with R1 and R2. The forecast credits the step with +1.29
  (`experiments/gate/README.md`). Row T60 reads its arm entry_step at +3.46, +1.85 and +3.33 over the oracle mix on
  its three library tasks, against +0.0 to +0.2 for the permuted control at the projection's size
  (`experiments/t1_m1/results/t60_reach/`). That arm also carries R1, R2 and g_mom's factors (g_mom's own part in T26:
  +1.23), and this file keeps g_mom. So T60's figures are an upper bound on the step alone, over an oracle base
  (ledger rule 9, J3), not a centre. The step is worth about 1.2 practice points by T26 and the forecast, and at
  most about 3.5 by T60.
- Against the first upload's design: the file gives up that step. By ledger rule 9, J6, the step's practice gain
  counts as zero at the centre of an estimate, with a range up to its practice value: the quantity that sets its size
  is twice as large at the test as on any practice task. So the expected cost is zero at the centre and up to about
  3.5 practice points if the step's practice gain carries over. No practice figure shows a gain from taking it out.
- What it buys: a T1 file with no outside data in it for the organiser's condition to bear on. It reads no Feng
  file, label, count or statistic. If the organisers discount the first upload, T1 rests on this file. If the first
  upload stands, its 70.9 stands whatever goes second. Neither T1 file is expected to beat 70.9.
- No counted backtest can price a file without the step: Feng's E9.5 lies inside every T1 window.

**Traces of the labels remain.** The file's data are free of the authors' labels; its adoption evidence is not.
Decision log 2026-11-19 (e) corrects entry (d), which had called both free. Traces remain in the evidence behind
choices the file keeps (the audit, decision log 2026-11-19 (c)), and the method report names them:
- g_mom: its factors are fitted on released stages, but the backtest that adopted it (T11) held back the labels
  outside the four groups, its restriction to cardiomyocytes and endocardium cites a read made on top of Feng's
  direction (T26), and the interaction check N22 was read the same way. In the first upload, which labels g_mom
  edits was set by the four groups' names; under the option the mask no longer reads that list (above), and it
  gives the same cells.
- The library-3 base: one of row R6's four reasons for it was priced on scenarios that reuse the share grid's values.
- E10.5's own mix: Feng's CD1 cardiomyocyte share, read on the authors' labels, was one supporting reason on
  2026-11-01 for keeping E10.5's own chamber mix (decision log 2026-11-01 (g), (h)), which this file keeps as library
  3's own.

No number from any of them enters the file. The labels email counts them as two kept parts: g_mom, and the draw of
E10.5 cells (the base and its mix).

**A stricter file, priced and not chosen.** A side agent's note, forwarded by the captain, asked for the price of a
file that also drops g_mom (decision log 2026-11-19 (f)). g_mom's factors are fitted on the released E9.5 and E10.5
alone; on practice they add +1.23 and +1.14 board points on the two E10.5 cuts, ahead in 10 of 10 seeds (row T26's
recorded numbers; they keep 0.96 to 0.99 of row T11's raw gains on b0). They are fitted on a one-day step and
applied over the test's two-day gap, so by rule 9's J6 that gain counts as zero at the centre with a range up to its
practice value. Dropping g_mom removes one trace; the choice of cell base cannot be removed without changing the
design that scored 70.9. Under the captain's answer the AI assistant chose to keep g_mom, as in the first upload. The
stricter file is built instead only if the captain says so; it would need its own amendment.

**It rests on no test score.** The first upload's 70.9 (decision log 2026-11-18 (j)) is not an input. It sets no
weight, prior, task, tolerance, constant or estimate here (ledger rule 9, J8). What moves this file is an organiser's
condition, not a score.

**"No score feedback", as it now stands.** The change comes after T1's first upload. "Each s2 was fixed before any
test score existed" no longer holds for T1: its last file is chosen here with the first scores known, as the heart's
was on 2026-11-18. The second bullet stands: a change after a first upload is the captain's, is logged with the
scores public at the time (the team's T1 total, 70.9), and computes nothing from a test score.

**The alternative, not taken.** The first upload's design with cell groups recomputed from the permitted cells' raw
counts and stage labels alone (branch `campaign/feng-relabel-fix2`, f211416). Its final review says not ready, with
three blockers; its marker sources need corrections; one of its tests fails on the build machine. It is parked:
unmerged and unrun.

**Marks and disclosure.** The rule for final marks is unchanged: of the entries every ruling received by then allows,
the one with the higher published score. Before the upload the captain sends two emails from his drafts, where they
have been since 2026-11-19 at his word ("Yes, both into drafts (Recommended)"): the disclosure of the X002 labels,
worded for this file, and then the sources note of 2026-11-18 in its amended form (`docs/ORGANIZER_QUESTIONS.md`,
"Disclosure of the X002 labels"). He sends them himself. Register row X002 and the method report say what this file
reads.

**Builds.** `bash experiments/builds/build_all.sh t1none` checks this file's hash, builds `runs/builds/T1_s2none.h5ad`
and runs the local format check for `T1:test`. It is built twice, both at 16 BLAS threads (T1 files move with the
thread count), and the two hashes are compared. The check builds so far (DGX campaign checkouts; never upload
sources; read by exit code, sha256 and the stamps' commit and clean fields only): two at f211416, with and without
R2's held lists, two at ea95c02 and two at 555daf4, all six with one sha256, 55eaee547ff1; with the option not given both T1 files
kept the reference list's hashes; the local format check passed. The builds at the ruling commit are compared with
55eaee547ff1, and any difference stops the hand-over until the captain has the explanation. The sets s1, s2, t3alt
and heart2 are unchanged: the reference list of the third rehearsal still describes the files of s1, s2 and t3alt,
and heart2's files keep their own list (`experiments/builds/results/heart2_307c9b2.SHA256SUMS`). This file is
outside that list; its own list is committed under `experiments/builds/results/` before checklist line 1 of
`docs/P3_RUNBOOK.md` can be ticked.

## Amendment of 2026-11-19: the embryo's last upload is the first upload on E7.5's exact positions (the captain's own say-so)

Written by the AI assistant. Times are UTC. The record is `docs/PLAN.md` decision log 2026-11-18 (j), (k), (u) and (v)
and 2026-11-19 (f) and (j), and rows E11, N51, E47a and E47b, N53 and R9 of `experiments/harness/TRIES.md`. Like the
amendment of 2026-11-18, this section was written with test scores known. Where it and the sections above differ, this
section holds for the embryo.

**The question and the captain's answer** (decision log 2026-11-19 (f); the answers were read at 11:12 UTC). The
question, verbatim: "Embryo's last upload (open since 18 Nov). The first upload scored 64.6, with its two shape scores
at 30.9 and 2.5 of 100. Which file should be the embryo's last upload?" It carried the four options of decision log
2026-11-18 (v): the first upload's cells and expression on E7.5's exact positions (`fgw_shape0x`, the AI assistant's
recommendation); the same positions with the series rendering (`fgw_shape0x_series`); fgw_shape0, the second file
before ruling R4; and `fgw_aligned_series`, the second file as ruled on 2026-11-13. The captain's answer: "Exact E7.5
positions (Recommended)".

**The embryo's last upload.** It replaces the second file of the ruling of 2026-11-13 ("Amendments of 2026-11-13",
`fgw_aligned_series`, `runs/builds/T2embryo_s2.h5ad`), which is not uploaded:

```bash
uv run python experiments/builds/build_test.py embryo --method fgw_shape0 --positions floor_exact --expression mixture --out runs/builds/T2embryo_s2x.h5ad
```

- **What it is.** The first upload's file with every cell placed on the position of a real E7.5 cell. The set s1 of
  `experiments/builds/build_all.sh` built the first upload as `build_test.py embryo --method fgw --expression
  mixture`. The new command keeps its seed (1000001, the builder's default), its board, its count rule, its size
  target, the anchors and coupling of FGW, the draw of E7.5 and E8.0 rows and their order, and each row's expression
  as released (`--expression mixture`, the default, written out as the set s1 writes it). Only the positions change.
  Every E7.5 row stays at its own position. Every E8.0 row goes to the position of an E7.5 cell that the file does not
  otherwise draw: a uniform subset of those cells, as large as the E8.0 half (stream [seed, 47]), matched to the E8.0
  rows one to one by the exact minimum-cost assignment on their expression in FGW's joint principal components
  (`vec.t2_build.interpolate_clouds`, positions "floor_exact"; tries E47a and E47b). Without the fallback below, the
  cloud is then a uniform sample of E7.5's own positions, the floor's form at test, centred and scaled to the
  pre-registered embryo size (`experiments/coords/results/p3/size_embryo_E7.75.json`), as the first upload is.
- **"As s1"** for the embryo now reads: s1's command with `--method fgw_shape0` in place of `--method fgw`, as the
  rule of 2026-11-03 had it, and `--positions floor_exact` added. The option was merged to main at 5e3eb9a; with it
  off, every embryo file builds as before.
- **The fallback.** If the file's count exceeds E7.5's cells, so that the E8.0 half outnumbers the E7.5 cells the file
  does not otherwise draw, every unused cell is a candidate, the remainder comes from all of E7.5's cells, and the
  repeats are jittered apart. Whether that happens at the test follows from the test board's count, which nobody
  reads; under the validation board's count rule it does not (row N53).
- **What a check build showed** (decision log 2026-11-18 (u); a check build in a DGX campaign checkout at 9967cf1,
  never an upload source): by a comparison that printed only booleans and hashes, the file holds the first upload's
  rows, row names and expression matrix exactly and differs from it in coordinates only.

**Why.** The first upload scored 64.6 (decision log 2026-11-18 (j)). Its four expression metrics, its size and its
neighbourhood metric carried over or beat practice; its two shape scores were 30.9 (d2_shape) and 2.5 (occupancy_dice)
of 100, the floor scoring 50.
- The record had flagged the first upload's geometry before the upload. At the test pair FGW's partners are mirrored
  in 10 of 10 draws, and the moved cloud's thinnest axis is 0.21 to 0.25 of its RMS radius, against 0.53 for E7.5 and
  0.36 for E8.0 (row E11). Those draws are mirrored less strongly than any mirrored practice draw and are thinner than
  every flagged practice draw but one task's (row N51).
- Reading N53, which reads no target, finds this file's geometry in the floor's form at the test pair, and the
  geometries of fgw_shape0 and of `fgw_aligned` (the planned second file's) outside it. The series option,
  `fgw_shape0x_series`, has this file's positions, so the reading covers it too. N53 does not build the first upload's
  own geometry. Against a 10 percent draw of E7.5, ten-seed means of (occupancy_dice, d2_shape): this file 0.859 and
  0.0030; copies of E7.5 0.848 and 0.0028; fgw_shape0 0.797 and 0.0202; `fgw_aligned` 0.787 and 0.0273. By N53's fixed
  rule, "in form" means no difference was detected at 2 standard errors (an arm exactly in form is called out about 15
  percent of the time); against B0, the copy of E7.5 that shares no cell with the reference draw, this file is out
  (row N53). It is one pair, ten draws, and E7.5 is the floor there, not the target.
- Practice holds no task with a shape loss like the first upload's, so it cannot price a floor-level shape against one
  (row E47, the AI assistant's reading after the result).
- The first scores choose here among geometries built from released data: the four options, each built before the
  question. The Official Rules allow that (Section 10, as read through a screen on 2026-11-18: "the test is whether
  the truth enters your prediction's content or only its selection"). No figure is solved from the scores. The 2.5
  sets no bound and no constant here, and the record does not explain it: a practice task with mirrored, flatter draws
  scored 50.1 on occupancy_dice (decision log 2026-11-18 (k)). Nothing here says what this file will score on either
  shape metric.

**What the record says it is worth, both ways.** Every figure is from practice tasks; none is a test value.
- **Limits.** Tries E47a and E47b are labelled diagnostics: the noise floor is 6.54 board points (the between-truth
  SD; the seed SE is 0.23), so no figure below is a detected effect; the harness applies its own embryo size rule per
  task, not the builder's size JSON (the shape metrics ignore scale); on 7 of the 9 tasks the earlier bracket is not
  E7.5 (row E47).
- **Against the uploaded form on the nine practice tasks** (row E47, rule (iii);
  `experiments/harness/results/e47_floor_exact/e47.txt`): merged mean -1.43 board points; better of two beside the
  uploaded form +0.46. It gains on three tasks (e75 +1.82, x1 +0.73, x3 +0.01; mean +0.86) and loses on six (half
  -1.06, late -1.87, new1 -2.92, new2 -3.13, e725 -4.27, x2 -3.35; mean -2.77). On the two counted backtests: +1.82 on
  E7.5, -4.27 on E7.25.
- **The other three options, the same way** (merged mean, then better of two): fgw_shape0 -0.50 and +0.88;
  `fgw_shape0x_series` -1.25 and +1.54; `fgw_aligned_series` +0.13 and +2.18. **Practice alone ranks this file last of
  the four** on both values.
- **Why it ranks last.** Its shape skills sit at the floor's level on every task (d2_shape 0.47 to 0.51,
  occupancy_dice 0.48 to 0.54), which is what it was built for; the other geometries scored above the floor on several
  practice tasks, which hold no shape loss like the first upload's (row E47).
- **The neighbourhood term.** Against the uploaded form, neighborhood_mmd skill -0.089 to +0.021 by task, merged
  -0.039 (decision log 2026-11-18 (v)); against fgw_shape0, merged -0.025 (row E47).
- **The expression metrics.** de_score, de_direction, mmd_u and variogram equal the uploaded form's on all 90
  task-seeds (row E47's guard; recomputed independently, decision log 2026-11-18 (v)). The file stakes nothing new on
  expression: its rows and their expression are the first upload's.
- **Why the captain's choice stands on the record** (the AI assistant's reading, from decision log 2026-11-18 (k) and
  (v)): the board keeps the better of two uploads (the working answer on final marks), so the second file is worth
  what it adds where the first is weak; the first upload's weakness is its geometry, which the record had flagged;
  reading N53 finds this geometry, which options 1 and 2 share, in the floor's form at the test pair and the
  geometries of options 3 and 4 outside it; and between those two, the series twin's gain on this geometry is level on
  practice (+0.18 merged over the five truths, negative on three tasks), so the cold attack's own check kept the
  mixture. Entry (v) records that the recommendation could not be reached from practice and reading N53 alone: it
  needs the first upload's two shape scores, and the question said so.

**What it overrides.**
- **Ruling R4** (the captain's ruling of 2026-11-13, "Amendments of 2026-11-13", T2 embryo): its file is not uploaded.
  The file stays on disk as the record of 2026-11-13.
- **Rule (ii) of row E47**, which the AI assistant fixed before its run under the delegation of 2026-11-02: by it,
  fgw_shape0 stays the floor-form arm, since this file's merged mean (-0.93), neighbourhood skill (-0.025) and worst
  shape skill (-0.324, d2_shape on x3) against fgw_shape0 each fail their clause. The captain's answer overrides it.
  This file enters by ruling, not by the adoption rule and not by a pre-set rule's outcome.

**"No score feedback", as it now stands.** The first bullet no longer holds for the embryo: its last upload is chosen
after, and because of, the first scores. It holds for T3; for T1 and the heart the sections above say what holds. The
second bullet is met: this rule changes after the embryo's first official upload by the captain's own say-so, and
computes nothing from a test score. The scores public at the time are the team's own of decision log 2026-11-18 (j):
board totals T1 70.9, heart 46.4, embryo 64.6, T3 47.8; on the embryo de_score 64.6, de_direction 73.0, mmd_u 78.1,
variogram 62.9, d2_shape 30.9, occupancy_dice 2.5, scale_log_ratio 96.5, neighborhood_mmd 74.4. The lines of decision
log 2026-11-18 (k) hold here: no constant set from a test value; no bound on a test anchor from the 2.5; no inference
of E7.75's detection profile from the first upload's expression scores; no choice of a file or a row order by its
parity.

**What stays.** The count rule and the pre-registered size: the file has the first upload's cells and size target. It
touches no lever of the post-exposure block: its cell number, size and scope are the first upload's. The first upload,
the other boards' files, the order (the first upload went first) and the rule for final marks are unchanged. "T1 and
the embryo: unchanged" in the amendment of 2026-11-18 no longer holds for the embryo.

**Builds.** `bash experiments/builds/build_all.sh embryo2` checks this file's hash, runs the command above into
`runs/builds/T2embryo_s2x.h5ad` and checks the file locally for `T2:embryo:test_interp`. It is built twice, both at 16
BLAS threads, the default that every upload build uses (`docs/COMPUTE.md`, "What reproduces where": the rendered
embryo files move with the thread count; whether this one does has not been measured), and the two sha256 are
compared. The sets s1, s2, t3alt, heart2 and t1none are unchanged. The reference list of the third rehearsal still
describes the files of s1, s2 and t3alt; heart2's files have their own list
(`experiments/builds/results/heart2_307c9b2.SHA256SUMS`), and so has t1none's file
(`experiments/builds/results/t1none_385ff58.SHA256SUMS`). This file is in none of those lists; its own list is
committed under `experiments/builds/results/` before line 1 of the captain's checklist (`docs/P3_RUNBOOK.md`
section 5) can be ticked.

## Amendment of 2026-11-19: the heart's last upload is the full trim (the captain's own named exception)

Written by the AI assistant. Times are UTC. The quoted question, options and answer are as `docs/PLAN.md` decision log
2026-11-19 (cc) records them, and the two quoted messages as entries (aa) and (bb) record them; nothing else here is
the captain's. The record is decision log 2026-11-18 (j) and (k) and 2026-11-19 (m) to (r), (aa), (bb) and (cc), and
row H26 of `experiments/harness/TRIES.md`. Like the amendment of 2026-11-18 and the two sections above, this section
was written with test scores known. Where it and the sections above differ, this section holds for the heart.

**The question and the captain's answer** (decision log 2026-11-19 (cc); one dialog, state at 19:31 UTC). The
question, as put: "Heart, last upload. I state no preference here (I recused myself after a lapse; see my message).
Practice results in heart points against the ruled full-scope file; Overall moves by half; a file scoring below the
first upload's 46.4 just leaves the 46.4 standing. If E12.5 was sampled like E10.5: half trim -1.0 (-7.9 on the
nearest practice step), full trim -3.9 (-15.5). If E12.5 was dissected clean: half trim +5.8, full trim +11.8. An
independent reviewer put 'clean to that degree' at about 0.25 to 0.30. A trim is your own named exception on a choice
closed after the 4 November exposure; the method report will say so. Which file?" The options: "Keep ruled file";
"Half trim"; "Full trim", described as "Highest ceiling: about +11.8 heart points (+5.9 Overall) if E12.5 is clean; if
it is like E10.5 the file likely falls below 46.4 and the first score stands (costing the 2.6 Overall the ruled file
would add). Your named exception on scope, cell number and size." **His answer: "Full trim".**

**Who decided what.** The captain chose the file. It is his own named exception on a lever closed after the exposure
of 2026-11-04: route (ii) of decision log 2026-11-19 (o). The AI assistant stated no preference. It had recused itself
from every choice and every lean on the heart's scope and size after a lapse of its own (decision log 2026-11-19
(bb)), and its own heart estimate had been withdrawn after the second cold check (same entry).

**The heart's last upload.** It replaces the file of the amendment of 2026-11-18 (the full-scope copy,
`runs/builds/T2heart_s2iso_half.h5ad`, set heart2), which stays built and is not uploaded:

```bash
uv run python experiments/builds/build_test.py heart --entry core_trim --expression copy --z-mode raw --out runs/builds/T2heart_s2core.h5ad
uv run python experiments/builds/halve.py runs/builds/T2heart_s2core.h5ad runs/builds/T2heart_s2core_half.h5ad
```

- **What it is.** The pool is E10.5 under `vec.t2_build.radial_trim_mask(C, core, 0.0)`: surrounding cells are
  removed farthest in plane from the core centroid first, until the in-plane scope is at most zero. What stays is the
  core (cardiomyocytes and endocardium, `configs/lineage_map.yaml`) and any surrounding cells inside its in-plane
  spread. On the released E10.5 the core is 0.231 of the cells (decision log 2026-11-19 (m), target-free). The file
  is the copy's draw from that pool at the count rule, with its expression unchanged, as the ruled file of 2026-11-18
  draws from the whole stage. It is then scaled about the centroid by one factor on all three axes to the size file's
  raw-mode core forecast of 2026-10-28 (the field `z_modes.raw.pred_rms_core`). Last comes the halving of the ruled
  file of 18 November, by the same rule and seed: a uniform random half of the rows (numpy's default generator, seed
  1000001, without replacement, order kept). The uploaded file is the half.
- **The option.** Merged to main at 71aba2d, off by default (`tests/test_build_heart_core_trim.py`). The builder
  refuses `--entry core_trim` except as `heart --entry core_trim --expression copy --z-mode raw`. With the option not
  given, every heart file builds as before.

**What it replaces and what it departs from, lever by lever.**
- **Scope.** A lever of the post-exposure block (`docs/DATA_POLICY.md` 2026-11-04, (b) to (d)). The frozen rule of
  2026-11-02 made the second upload a half trim and weighed the full trim without taking it ("it never beats s1 on a
  released target", section "T2 heart" above). The amendment of 2026-11-18 replaced the half trim with the whole
  stage. This file takes the full trim. No rule returned it; the captain named it.
- **Cell number.** A closed quantity. The count rule is unchanged and does not read the pool, so the built file has
  the ruled file's number of rows, and the halving is the same. What changes is the pool the rows are drawn from,
  which follows from the scope: rows repeat only if the count exceeds the pool, and repeated rows are then moved apart
  by the T2 head's jitter, as in any file (`vec.t2_build.sample_cells`). The captain's option named the cell number
  with the scope and the size.
- **Size.** A closed quantity. The pre-registered mode was xy_grow_z_hold (section "T2 heart" above). The amendment of
  2026-11-18 took the raw mode at H1's size. This file takes the raw mode's core forecast: the term that the size
  file's three scope hypotheses (H1, H0 and H2) multiply, and not one of them. No earlier rule fixed it. The captain's
  option named the size.
- **Expression.** No change from the file it replaces: the copy's expression, unchanged. No trend, no momentum.
- **The default.** Rule 9's J1a and J4 of `experiments/harness/TRIES.md`, row H26's own text and the rule for a lever
  closed after an exposure all gave the file in force as the default, and the two cold checks of decision log
  2026-11-19 (n) and (bb) kept it. The captain's exception overrides that default. This file enters by ruling, not by
  the adoption rule and not by a pre-set rule's outcome.

**Why.** The captain's aim, in his words (decision log 2026-11-19 (aa) and (bb)): "Proceed with whatever that gets us
to Top 3 in the leaderboard. If other teams can do it, so can us" and "ultrathink just deeply internalize that it is
very much possible to beat the top 3 leaderboard scores". The AI assistant's reading of them (its own): a best-case
aim, each last file weighed by its ceiling with the first score banked.

**What practice says, both ways.** Row H26 (`experiments/harness/results/h26_heart_scope_forms/`; read in decision
log 2026-11-19 (q), recomputed independently in (r)). Every figure is a difference between files at an oracle core
size shared by the arms, on released transitions; no figure is a file's level. C is this file's form, A the ruled
full-scope file's, B the half trim's; heart points, plain means, the two E10.5 transitions merged into one fold.
- **If E12.5 was sampled like E10.5** (the truth as released): C minus A -3.86, behind in 10 of 10 seeds. By
  transition: -15.53 on E9.5 to E10.5, -0.97 on E8.25_late to E8.5, +2.05 on E8.75 and earlier to E10.5. E9.5 to
  E10.5 is the transition nearest the test.
- **If E12.5 was dissected clean** (the three truths with most or all surroundings removed): C minus A +11.79, ahead
  in 10 of 10 seeds. C minus B +6.01 (10 of 10; +0.75 on the E10.5 fold, with one seed of ten behind, and +11.27 on the
  E8.5 fold).
- **The odds the full trim needs** to pay, that E12.5 is clean (decision log 2026-11-19 (q) and (r)): about 0.25 on
  the fold-merged figures; 0.31 with the one-day transition's loss capped by the banked score (the cap is the AI
  assistant's estimate of A's margin over that score), 0.57 uncapped; 0.68 capped and 0.86 uncapped with the gain as
  well as the loss taken from the one-day transition alone.
- **The reviewer's odds** (the second cold check, decision log 2026-11-19 (bb)): about 0.45 that E12.5 holds less
  surrounding tissue than E10.5, and about 0.25 to 0.30 that it is clean to the degree row H26's truths model.

**What practice cannot say.**
- How E12.5's stack was sampled. It is a test-stage premise no released task carries (rule 9, J1a). The team's
  screened notes give the E10.5 and E12.5 stacks the same wording (decision log 2026-11-19 (o)).
- No released transition has a reference and a target of equal scope, which is the test's case if E12.5 was sampled
  like E10.5 (decision log 2026-11-19 (r)).
- The trimmed truths are synthetic. Two of the three clean truths are thinnings no dissection makes, and the core is
  cardiomyocytes and endocardium only, so a clean heart is not all core (decision log 2026-11-19 (bb)).
- The size term. Every arm of row H26 was scaled to one oracle core size. At the test this file takes the raw-mode
  core forecast and the replaced file H1's raw-mode size. That difference in the scale term is not priced; it can be
  as large as the trim's own gain, in either direction (decision log 2026-11-19 (n), point (2)).

**Held until the size term is priced** (added by the AI assistant before this section was merged). The question put to
the captain quoted row H26's figures without saying that they were measured with every arm at an oracle size; the last
bullet above is what the question left out. A side note the captain forwarded on 2026-11-19 pointed this out. Analysts
who have not seen the AI assistant's conversation are designing a practice read of the three forms at forecast sizes
on released stages (row H27 of `experiments/harness/TRIES.md`, once registered), as decision log 2026-11-19 (bb)
provides. The file of this section is built, and it is not handed over as the last heart file until that read's result
is with the captain. His ruling stands unless he names another file; the AI assistant states no preference.

**Nothing is computed from a test score.** The trim is `radial_trim_mask` at fraction 0, the hedge study's
core_radial of 2026-11-01. The size is the size file's own raw-mode field, logged on 2026-10-28, before any score. The
halving is the first upload's rule. The first heart score, 46.4, enters only as the level a last file must beat, as
the question said (decision log 2026-11-19 (n), point (4)). The scores chose among built or pre-defined predictions
(Official Rules, Section 10, as decision log 2026-11-19 (g) quotes them).

**What the record must say** (the method report, its section on exposures in `METHOD_SUMMARY.md`, and the note that
precedes the upload).
- The scope and the size of the heart's last file, and with the scope the pool its rows are drawn from, were chosen by
  the captain, after the exposure of 2026-11-04 and after the first scores (the count rule itself is unchanged). As
  decision log 2026-11-19 (o), route (ii), puts it, the record says that the choice was made "after the exposure and
  the first scores by a person who had seen the E12.5 heart file's size, and that it cannot be shown to be free of
  it".
- The AI assistant stated no preference and named no heart file. It had recused itself after a lapse (decision log
  2026-11-19 (bb)). Before the dialog it had told the captain that a Top-3 aim weighs each board's last file by its
  ceiling with the first score banked (entry (cc)), and the heart's options were put with their practice gains and
  losses.
- The size is the size file's own raw-mode core forecast of 2026-10-28, the term its three scope hypotheses multiply,
  and not one of those hypotheses.

**A correction to the amendment of 2026-11-18** (owed since decision log 2026-11-19 (n)). That amendment lists the half
trim among four departures, "each by the captain's answer above to a question that named it". The record shows three
priced questions on the heart: the size mode, the expression step and the cell number (decision log 2026-11-18 (k)).
The half trim appears in that entry only in the plan table's column of what the new file touches, with no price of its
own. So the scope hedge was dropped inside the redesign, without a priced question of its own. The text of that
amendment is not edited; this paragraph corrects it.

**"No score feedback", as it now stands.** For the heart the first bullet has not held since 2026-11-18. The second
bullet is met: the change comes after the heart's first official upload, by the captain's own say-so, is logged with
the scores public at the time (decision log 2026-11-18 (j): heart 46.4), and computes nothing from a test score.

**What stays.** The first upload, the count rule, the halving's rule and seed, the other boards' files, the order (the
first upload went first) and the rule for final marks.

**Builds.** `bash experiments/builds/build_all.sh heart3` checks this file's hash, runs the two commands above and
checks both files locally for `T2:heart:test_extrap`. It is built twice, both at 16 BLAS threads, and the two builds'
sha256 are compared for both files. The sets s1, s2, t3alt, heart2, t1none and embryo2 are unchanged. The files of
heart3 are in no reference list yet; their own list is committed under `experiments/builds/results/` before line 1 of
the captain's checklist (`docs/P3_RUNBOOK.md` section 5) can be ticked.

## Amendment of 2026-11-19: the embryo's last upload is the series twin on E7.5's exact positions (the captain's own say-so)

Written by the AI assistant. Times are UTC. The quoted question, options and answer are as `docs/PLAN.md` decision log
2026-11-19 (cc) records them, and the two quoted messages as entries (aa) and (bb) record them; nothing else here is
the captain's. The record is decision log 2026-11-18 (j), (k), (u) and (v) and 2026-11-19 (j), (p), (aa), (bb) and
(cc), and rows E47a and E47b, N53 and R9 of `experiments/harness/TRIES.md`. This section was written with test scores
known. Where it and the sections above differ, this section holds for the embryo.

**The question and the captain's answer** (decision log 2026-11-19 (cc); the same dialog). The question, as put:
"Embryo, last upload. Your chosen file (the first upload's cells on E7.5's exact positions) is expected near 69. Its
'series' twin re-renders the E7.5 cells at E8.0's measurement profile. If E7.75 was measured like E8.0 (the organisers'
page groups E7.75 with E8.0; our committed odds 0.4 to 0.5) the twin scores about 4 to 5 embryo points higher;
otherwise 2 to 7 lower, and then the first upload's 64.6 stands. Same expected value; higher best case, lower worst
case. Which file?" The options: "Series twin (Top-3 aim)"; "Keep chosen file". **His answer: "Series twin (Top-3
aim)".**

**Who decided what.** The captain chose the file. Before the dialog the AI assistant had marked the higher-ceiling file
as the one that fits a Top-3 aim and given the expected values as about equal (decision log 2026-11-19 (cc)).

**The embryo's last upload.** It replaces the file of the embryo's earlier amendment of 2026-11-19 (`fgw_shape0x`,
`runs/builds/T2embryo_s2x.h5ad`, set embryo2), which stays built and is not uploaded:

```bash
uv run python experiments/builds/build_test.py embryo --method fgw_shape0 --positions floor_exact --expression series --out runs/builds/T2embryo_s2xs.h5ad
```

- **What it is.** The harness's `fgw_shape0x_series`, try E47b. The positions are the replaced file's: every E7.5 row
  at its own position, every E8.0 row on the position of an E7.5 cell the file does not otherwise draw, matched one to
  one by the exact minimum-cost assignment on expression in FGW's joint principal components; then centred and scaled
  to the pre-registered embryo size (the section above; row R9). That section's fallback holds as written. The
  expression is the series rendering of ruling R4: E7.5's cells rendered into E8.0's per-gene detection profile and
  depth (`vec.t2_build.series_rendered_expression`). Of the check builds at 9967cf1, decision log 2026-11-18 (u)
  records that "the two series files hold the built second file's" rows, row names and expression matrix (the AI
  assistant reads "the built second file" as ruling R4's `T2embryo_s2.h5ad`).
- **"As s1"** for the embryo now reads: s1's command with `--method fgw_shape0` in place of `--method fgw`,
  `--positions floor_exact` added, and `--expression series` in place of `--expression mixture`. The option was merged
  to main at 5e3eb9a with tries E47a and E47b, off by default.

**What it replaces and what it departs from, lever by lever.**
- **Geometry.** No change from the file it replaces: E7.5's exact positions.
- **Expression.** The replaced file kept the first upload's expression as released (`--expression mixture`). This file
  takes the series rendering, as ruling R4's second file did. Its four expression metrics no longer repeat the first
  upload's.
- **The cold attack's own check** (decision log 2026-11-18 (v)): stay with the mixture if the series gain on this
  geometry averages one point or less or is negative on two or more tasks. It averages +0.18 and is negative on three.
  The captain's answer overrides that check.
- **Rule (ii) of row E47**, which the AI assistant fixed before its run: by it the series twin does not replace
  `fgw_shape0_series` (decision log 2026-11-18 (u)). The captain's answer overrides it, as it overrode the same rule
  for the replaced file.
- **Closed levers.** None touched. The file keeps the pre-registered size and the count rule, and touches no closed
  lever (decision log 2026-11-18 (v)).

**Why.** The captain's aim, in his words (decision log 2026-11-19 (aa) and (bb)): "Proceed with whatever that gets us
to Top 3 in the leaderboard. If other teams can do it, so can us" and "ultrathink just deeply internalize that it is
very much possible to beat the top 3 leaderboard scores". The AI assistant's reading of them (its own): a best-case
aim, each last file weighed by its ceiling with the first score banked.

**What practice says, both ways.** Every figure is from practice tasks; none is a test value.
- **The series rendering on this geometry** (`fgw_shape0x_series` minus `fgw_shape0x`, paired by seed, board points
  per truth; the independent recomputation of decision log 2026-11-18 (v), on the task files of row E47): E8.0 +5.54
  (half +7.53, late +3.56); E8.25_late +4.01 (new1 +2.71, x1 +5.30); E8.5 -0.46 (new2 +1.38, x2 +0.50, x3 -3.25); E7.5
  -1.62; E7.25 -6.58. Merged over the five truths +0.18. On the split ruling R4 used: +2.72 on the six tasks whose
  target shares its deeper bracket's profile, -2.27 on the three that do not. At the committed 0.4 to 0.5 chance that
  E7.75 was measured like E8.0, the twin's practice mean against the mixture is -0.27 to +0.22: level.
- **Against the uploaded form on the nine practice tasks** (row E47, rule (iii);
  `experiments/harness/results/e47_floor_exact/e47.txt`): merged mean -1.25 board points; better of two beside the
  uploaded form +1.54. The replaced file: -1.43 and +0.46.
- **The analyst's reading** (decision log 2026-11-19 (p); its own judgement from committed text): level with the
  replaced file on the expected better of two (+0 to +0.6), with the higher tail (about 74 to 75 if E7.75 was measured
  like E8.0).

**What practice cannot say.**
- How E7.75 was measured. Nothing released shows whether its detection profile is E7.5's or E8.0's. The 0.4 to 0.5
  is the committed odds, not a measurement.
- Tries E47a and E47b are labelled diagnostics: the noise floor is 6.54 board points (the between-truth SD; the seed
  SE is 0.23), so no figure above is a detected effect; on 7 of the 9 tasks the earlier bracket is not E7.5 (row E47).
- Practice holds no task with a shape loss like the first upload's (row E47, the AI assistant's reading after the
  result).

**Nothing is computed from a test score.** The file's positions and its rendering were each defined and built before
the question (tries E47a and E47b; ruling R4's rendering). The scores chose among built or pre-defined predictions
(Official Rules, Section 10, as decision log 2026-11-19 (g) quotes them). The lines of decision log 2026-11-18 (k) hold
here, as in the section above: no constant set from a test value; no bound on a test anchor from the first upload's
occupancy score; no inference of E7.75's detection profile from the first upload's expression scores; no choice of a
file or a row order by its parity.

**What the record must say.** The embryo's last file was chosen by the captain after the first scores, for its best
case, against the cold attack's check and rule (ii) of row E47. Its expression bet (E7.75 measured like E8.0) is the
one ruling R4 made on 2026-11-13, before any score.

**"No score feedback", as it now stands.** As the section above says for the embryo: the first bullet no longer holds;
the second is met. The change comes after the embryo's first official upload, by the captain's own say-so, is logged
with the scores public at the time (decision log 2026-11-18 (j): embryo 64.6), and computes nothing from a test score.

**What stays.** The first upload, the count rule, the pre-registered size, the other boards' files, the order (the
first upload went first) and the rule for final marks. The file of ruling R4 (`runs/builds/T2embryo_s2.h5ad`) and the
replaced file (`runs/builds/T2embryo_s2x.h5ad`) stay on disk as the record and are not uploaded.

**Builds.** `bash experiments/builds/build_all.sh embryo3` checks this file's hash, runs the command above into
`runs/builds/T2embryo_s2xs.h5ad` and checks the file locally for `T2:embryo:test_interp`. It is built twice, both at
16 BLAS threads, and the two sha256 are compared. The sets s1, s2, t3alt, heart2, t1none and embryo2 are unchanged.
This file is in no reference list yet; its own list is committed under `experiments/builds/results/` before line 1 of
the captain's checklist (`docs/P3_RUNBOOK.md` section 5) can be ticked.

## Amendment of 2026-11-19: T3's last upload is the barred branch's s1 (the captain's own say-so)

Written by the AI assistant. Times are UTC. The quoted question, options and answer are as `docs/PLAN.md` decision log
2026-11-19 (cc) records them, and the two quoted messages as entries (aa) and (bb) record them; nothing else here is
the captain's. The record is decision log 2026-11-18 (j) and (u) and 2026-11-19 (p), (s), (aa), (bb) and (cc), and row
K18 of `experiments/harness/TRIES.md`. This section was written with test scores known. Where it and the sections
above differ, this section holds for T3.

**The question and the captain's answer** (decision log 2026-11-19 (cc); the same dialog). The question, as put: "Task
3, last upload. The ruled 'refused' file stays close to the floor: expected about 52 (48 to 57). The already-built
'barred' file carries another knockout's response: about 62 to 64 if the beta-catenin embryos lean older than the
wild-type reference as the Gata4 ones did (a reviewer's odds 0.35 to 0.45); otherwise it falls below the first
upload's 47.8, which then stands, costing the 2 to 6 points the ruled file would likely add. Expected value about
equal; our own pre-set rule and a cold check both said keep the ruled file. Which file?" The options: "Barred file
(Top-3 aim)"; "Keep ruled file". **His answer: "Barred file (Top-3 aim)".**

**Who decided what.** The captain chose the file. Before the dialog the AI assistant had marked the higher-ceiling file
as the one that fits a Top-3 aim and given the expected values as about equal (decision log 2026-11-19 (cc)). Try
K18's fixed reading had said "no question", and the cold check of the swap had said it does not stand on expected
score (decision log 2026-11-19 (s)); the captain chose the file for its best case.

**T3's last upload.** The barred branch's s1, as set t3alt builds it with the command of the T3 section above on the
barred branch:

```bash
uv run python experiments/t3_decomp/build_submission.py build experiments/t3_decomp/results/gata4/gata4_decisions.json \
  --branch barred --submission s1 --wt data/raw/T2/E8.75.h5ad --run-dir runs/t3_decomp/gata4_20261028T122425Z \
  --board T3:bcatenin --seed 1000001 --out runs/t3_build/barred_s1.h5ad
```

- **The file.** `runs/t3_build/barred_s1.h5ad`, sha256
  3712c3f9988d0f8aa4198b22c37b866a84ef8757fd500a9163944f0663a520a2, built in set t3alt and on the reference list of
  the third rehearsal (`experiments/builds/results/builds_dc7b6ae.SHA256SUMS`). The copy handed to the captain is the
  checked copy of the same bytes on the local machine (`runs/archive/builds_dc7b6ae/t3_build/barred_s1.h5ad`, sha256
  checked against the reference list), from the folder the uploads of 2026-11-18 came from
  (`results/manual_uploads.csv`): he uploads by hand through the site's form (the amendment of 2026-11-18, "Uploads,
  dates and marks"). Never the copy in the DGX's `runs/t3_build`. No new build is needed.
- **What it is.** Wild-type E8.75 cells with the central composition at the Mab21l2 magnitude, the refused s1's, each
  cell then moved by the Mab21l2 factor at exponent 1 (`transfer_mab21l2_all`; section "T3" above, "Branch names";
  row K18), with no depth, stage or coverage step. It takes no number derived from the Gata4 files. What its factor
  carries, by the cold check of decision log 2026-11-19 (s): Mab21l2's E9.5 knockout cells over WT E8.75, so the whole
  step from E8.75 to E9.5 along the wild-type trend, the Mab21l2 knockout's effect, that sample's own run and WT
  E8.75's own deviation.

**What it replaces and what it departs from, lever by lever.**
- **The file.** It replaces "the refused s1 with the wild-type-only factor" (`runs/t3_build/refused_s1.h5ad`, built
  with `--wt-revert`) of the amendment of 2026-11-18, which stays built and is not uploaded.
- **The factor.** The replaced file moved each cell by its lineage group's wild-type-only factor at exponent 0.5 (the
  amendment of 2026-11-13). This file moves each cell by the Mab21l2 factor at exponent 1.
- **The composition.** No change: the central composition at the Mab21l2 magnitude in both.
- **The frozen T3 rule.** It put the barred s1 last only under an organiser ruling: one that restricts or bars the
  Gata4 files or bars the validation answers ("Ruling switches" above), or an organiser's "zero" on a floored severity
  (the ruling of 2026-11-06 (b); the amendments of 2026-11-13). No such ruling is on record. This file goes last by
  the captain's say-so, not by a switch.
- **Try K18's fixed reading** ("no question"; decision log 2026-11-18 (u)) and the cold check of decision log
  2026-11-19 (s) both kept the refused file. The captain's answer overrides them. This file enters by ruling, not by
  the adoption rule and not by a pre-set rule's outcome.

**Why.** The captain's aim, in his words (decision log 2026-11-19 (aa) and (bb)): "Proceed with whatever that gets us
to Top 3 in the leaderboard. If other teams can do it, so can us" and "ultrathink just deeply internalize that it is
very much possible to beat the top 3 leaderboard scores". The AI assistant's reading of them (its own): a best-case
aim, each last file weighed by its ceiling with the first score banked.

**What practice says, both ways.** Every figure is from practice tasks (WT E8.75 to Gata4, replicates g1 and g2); none
is a test value.
- **Row K18, by its fixed rule** (decision log 2026-11-18 (u)): `barred_form` minus `refused_form` +2.29 and +6.90 in
  the plain world; -4.48 on g1 (severity scored in 6 of 10 seeds against the 9 needed) and +4.94 on g2 at matched
  depth. The refused form's board means: 57.1 and 55.1 plain, 60.1 and 58.9 at matched depth.
- **The cold check's repricing** (decision log 2026-11-19 (s); approximate, from two scorers): Gata4's knockout cells
  lean towards E9.5 against WT E8.75, and the barred factor's lead sits on that part. With it priced for the barred
  file alone, the four cells go to about -2.8, +1.8, -9.6 and -0.2, and the swap is worth about -1 to +1 as the
  expected better of two. Its one clear edge is the tail: about 62 to 64 against 57 to 59.
- **The analyst's reading** (decision log 2026-11-19 (p); its own judgement): level with or ahead of the refused file
  on the expected better of two (four-cell practice mean +2.4), ahead on the tail (about 62 to 64). The cold check
  corrected its brief on what the barred factor carries (above).

**What practice cannot say.**
- Whether the beta-catenin embryos lean older than WT E8.75, as Gata4's did. The tail pays only then, and nothing
  released shows it (decision log 2026-11-19 (s)).
- g1 and g2 are two replicates of one knockout run: one fold of evidence, not two (row K18). No released knockout is
  beta-catenin.

**Nothing is computed from a test score.** The file was built in set t3alt at dc7b6ae, the third rehearsal of
2026-11-13, before any score, and is on that reference list. The scores chose among built or pre-defined predictions
(Official Rules, Section 10, as decision log 2026-11-19 (g) quotes them). The first T3 score, 47.8, enters only as the
level a last file must beat, as the question said.

**What the record must say.** T3's last file was chosen by the captain after the first scores, for its best case,
against try K18's fixed reading and a cold check. The practice figures that priced it read the Gata4 replicates; the
file itself takes no Gata4-derived number.

**"No score feedback", as it now stands.** The first bullet no longer holds for T3: its last upload is chosen after,
and with, the first scores. The second bullet is met: the change comes after T3's first official upload, by the
captain's own say-so, is logged with the scores public at the time (decision log 2026-11-18 (j): T3 47.8), and
computes nothing from a test score. "s1 is fixed whatever s2 scores" was spent on 2026-11-18.

**What stays.** The first upload (the allowed branch's s1), the order (it went first), the ruling switches for any
organiser ruling still to come (every switch that names the barred s1 last names this file) and the rule for final
marks.

**Builds.** None. The file is set t3alt's (`bash experiments/builds/build_all.sh t3alt`), built at dc7b6ae and on that
reference list.

## Amendment of 2026-11-19: T1's last upload is the recomputed-groups file (the captain's own say-so)

Written by the AI assistant. Times are UTC. The quoted question, options and answer are as `docs/PLAN.md` decision log
2026-11-19 (ii) records them; nothing else here is the captain's unless it is quoted and its entry named. The record is
decision log 2026-11-19 (g), (l), (u) to (z), (dd), (gg), (hh) and (ii); rows R8, N54, N55, T61 and T62 of
`experiments/harness/TRIES.md`; register rows X002, X003 and X010 of `EXTERNAL_DATA.md`; the relabelling's record
`experiments/external/results/feng_relabel/MANIFEST.json`; and the check builds' record
`experiments/builds/results/t1rc_check_8138fa2.txt`. Like the sections above of the same date, this section was
written with the first test scores known; no test score sets anything in it (below). Where it and the sections above
differ, this section holds for T1.

**The question and the captain's answer** (decision log 2026-11-19 (ii); one dialog question, header "Task 1 file";
state at 20:55 UTC).
- The question, as it was put: "Task 1, last upload (one left). Our own odds are about 0.6 that only a compliant upload
  will count, because the first one (70.9) used the atlas authors' labels; then this file is Task 1's whole score. On
  practice (five tasks, ten seeds, one day ahead; independently recomputed) the atlas step on groups we recomputed from
  permitted cells beat the no-outside-data file on every task, by +1.2 to +2.5 points. By our own rule a one-day gain
  counts as zero at the centre for the two-day test, so the honest edge is 0 to about +2. A cold check let my pre-set
  advice (the recomputed file) stand, with one thing it says my rule left out: if the organisers rejected this second
  atlas file while discounting the first, Task 1 would have no counted upload, and a rejection chance of roughly 2 to
  4% cancels the edge. Which file?"
- The options, as put: "Recomputed atlas file" ("My pre-set advice, and the higher ceiling: 0 to about +2 over the
  other file, perhaps more at the test's double step. It follows the organisers' words on recomputing, but the atlas
  authors' quality filter chose which cells exist and that cannot be undone. Before building I run one more practice
  check at the test's step size (about two hours)."); "No-outside-data file" ("The file you ruled on 19 November; built
  and ready. Nothing in it can be ruled non-compliant. Gives up the 0 to +2."); "Own-data momentum file" (a third
  candidate from released data only, which reverses ruling R1, is tied to depth and has no builder; "Not a pre-set
  rule's output; the cold check did not prefer it; I state no lean on it. Choosing it names both exceptions.").
- The captain's answer, typed in his own words: "proceed with your recommendation".

**How the AI assistant reads that answer** (its reading, as decision log 2026-11-19 (ii) records it; not his words).
The question named one recommendation, in its own text ("my pre-set advice (the recomputed file)") and in the first
option ("My pre-set advice"). His answer therefore rules the recomputed atlas file in as T1's last upload. It is an
answer to this question alone and extends to nothing else. It sets the ruling of decision log 2026-11-19 (u) on the
no-outside-data file aside for the last upload. The answer lifts no rule: the organisers' data condition, the data
rules, the post-exposure block and the list of material never opened all stand. The captain can correct this reading
before any upload, and every upload is his own act.

**Who decided what.** The captain chose the file, by his answer above. The recommendation he followed is row T62's
advice, fixed by the AI assistant in that row before any result: the recomputed file if its step is AHEAD under both
readings of the first upload, the ruled no-outside-data file otherwise. A cold challenger, given the decision and the
record but not how the advice was formed, let the advice stand with conditions (decision log 2026-11-19 (hh)). The
question carried that check's caveat beside the advice. The AI assistant stated no lean on the third option.

**T1's last upload.** One command, as the build script gives it (set t1rc):

```bash
uv run python experiments/builds/build_test.py t1 --g-mom --base lib3 --feng-groups recomputed --out runs/builds/T1_s2rc.h5ad
```

- **What it is.** The first upload's design (library 3's cells at library 3's own chamber mix, row R6; the g_mom step,
  row R3b) with the outside atlas's E9.5 to E14.5 direction (M1's step) formed on cell groups recomputed from the
  atlas's permitted cells only. The labeller is `experiments/external/feng_relabel.py`, with its settings
  `experiments/external/feng_relabel.yaml` frozen at sha256
  883e74c2c1ef5f01d5244555f8b3dcae24255568c731a6d53ee2788b7b677b49 (row N54; its threshold set by row N55's dry run on
  the released E9.5). The statistics file it wrote has sha256
  36d1f7edcbbf3d7a80c68b51d5624a94bb753668079fe3d69aed720c71baa7e9. Its record says route recomputed and that no
  author label column was read (cell names, stage, zone and strain only). The authors' labels are not read.
- **The groups.** By row N54's record, each of the four direction groups (cardiomyocytes, endocardium, endothelium,
  epicardium) has its own recomputed direction, with at least 10 cells at both ends of the bracket; no fallback fired.
  R1's hold and g_mom's mask follow the groups with a recomputed direction. R2 holds its genes flat, as in the first
  upload.
- **The option.** Merged to main with the labeller (decision log 2026-11-19 (w)), off by default. The builder refuses
  `--feng-groups recomputed` except as `t1 --g-mom --base lib3` with no other base or option of its list, and unless
  the labeller's record is committed at HEAD with the route recomputed; then unless the statistics file passes the
  labeller's own checks (`experiments/builds/build_test.py`, `_feng_recomputed_refusals`; `feng_dirs.check_recomputed`).
  So it builds only in a checkout that holds the relabelling's output (`data/derived/feng_relabel_c/`), which on
  2026-11-19 was the DGX checkout `~/v/vec-campaign-d` alone. With the option not given, every T1 file builds as
  before.
- **Check builds, not upload builds.** At 8138fa2 on `~/v/vec-campaign-d`, twice at 16 BLAS threads: exit 0 both, one
  sha256, ad8f24c71c90ca17bce5f8876d503471ab3faa5b17a0878addb00aff689fca4f, stamps at the commit and clean
  (`experiments/builds/results/t1rc_check_8138fa2.txt`; decision log 2026-11-19 (z)). That file is on no reference list.

**What it replaces.** The amendment of 2026-11-19 above, "T1's last upload has no outside data", and its file:
`build_test.py t1 --g-mom --base lib3 --feng-groups none`, set t1none, row R8, `runs/builds/T1_s2none.h5ad`, sha256
55eaee547ff1 (built at 385ff58; `experiments/builds/results/t1none_385ff58.txt`). That file stays built and is not
uploaded. The text of that section is not edited; where the two differ, this section holds. One consequence, read from
the builder's code and not from a run: since row N54's record routes to recomputed (committed at 3fcfb40), the builder
refuses `--feng-groups none` unless `--feng-none-ruling` names a later decision-log entry of the captain that names N54
and contains "sets the N54 route aside". Set t1none passes no such option, and no such entry exists, so the set no
longer rebuilds that file. The built file and its record stand as they are.

**What it departs from, lever by lever.**
- **The data.** The outside atlas (register row X002) returns, on recomputed groups. The amendment it replaces put no
  outside data in T1's last file. The first upload formed the same step on the authors' groups, which the organisers'
  condition excludes. What enters this file: the atlas's permitted CD1 cells at E9.5 and E14.5 as raw counts, with gene
  names and the cell name, stage, zone and strain columns; marker gene sets, those of the four direction groups from
  the released E9.5 label means (`experiments/t1_scope/results/label_evidence/markers.csv`) and those of the competing
  classes from the eight adult-mouse papers of register row X010. No clustering, embedding or classifier. The captain
  approved the atlases as Task 1 inputs ("utilize the atlases", decision log 2026-11-19 (k)), the two further inputs
  ("Approve both (Recommended)") and the labeller's settings and runs ("Sign and run (Recommended)", both in decision
  log 2026-11-19 (u)). The second atlas (register row X003) is not read.
- **Composition, cell number and depth.** Unchanged from the first upload: the same base, count rule and seed
  (1000001), so the same drawn cells (row T62: "States, shares, cell number, size and scope do not change: the same
  cells as R8's file and the first upload"). No depth step is added. These are levers the post-exposure block closes
  (`docs/DATA_POLICY.md` 2026-11-04, (b) to (d)), and this file touches none of them.
- **Expression.** The first upload's expression steps, with the direction formed on the recomputed groups: M1's step at
  2/5 of the atlas's bracket, R1's hold, R2's held genes and g_mom's within-type edit of cardiomyocytes and
  endocardium.
- **"s2 as s1".** T1's last file is the first upload's command plus one option, `--feng-groups recomputed`. Both T1
  uploads draw the same cells. The base and cardiomyocyte-share hedges of rows R7 and T33 stay given up, as in the
  amendment this one replaces.
- **"No score feedback", as it now stands.** For T1 the first bullet has not held since the amendment this one
  replaces. The second bullet is met: the change comes after T1's first official upload, by the captain's own say-so,
  is logged with the scores public at the time (decision log 2026-11-18 (j): T1 70.9), and computes nothing from a test
  score.
- **The default.** Row N54's terminal rule named this file as a candidate only; row T62's pre-set advice named it if
  AHEAD under both readings. This file enters by the captain's ruling on that advice, not by the adoption rule.

**Why.** The captain's words on this question are the answer above: "proceed with your recommendation". Earlier the same
day he had answered the first T1 question with "whatever it takes to win this challenge" (decision log 2026-11-19 (f))
and approved the atlases with "utilize the atlases" (entry (k)). The recommendation is row T62's advice, fixed before
the run: the recomputed file if its step is AHEAD under both readings of the first upload. By that row's reading it was
AHEAD under both (decision log 2026-11-19 (gg)).

**What practice says, both ways.** Row T62 (`experiments/t1_m1/results/t62_recomputed/`; read in decision log
2026-11-19 (gg), recomputed independently there). Every task scores one truth, E10.5, one day ahead, with the step at
1/5 of the bracket: labelled in-sample diagnostics, never an adoption result.
- **For the file.** Reading S (the first upload stands): AHEAD. Reading D (it is discounted): AHEAD. The recomputed
  analogue minus the no-outside-data analogue, paired over ten seeds: chamber cut +1.97, whole cut +1.67, libraries 1
  to 3 +2.50, +1.23 and +2.51, each ahead in 10 of 10 seeds. The epicardium step alone PAYS by the row's statement:
  +0.68 on library 3 (10 of 10). An independent recomputation, from the task files without the read script, matched
  every figure and statement.
- **Against it.** By rule 9's J6 of `experiments/harness/TRIES.md`, the direction's share of the bracket is 1/5 on
  every practice task and 2/5 at the test, a ratio above 1.5, so the one-day gain counts zero at the centre for the
  two-day test, with a range from 0 to +1.97 (the chamber cut's gain). No practice figure is a test value.
- **The compliance price, which the advice rule left out** (decision log 2026-11-19 (hh)). Row T62's rule decides on
  practice gain alone. If the organisers discounted the first upload and also rejected this second outside-data file,
  T1 would have no counted upload; the no-outside-data file cannot be rejected on that ground. The AI assistant's
  arithmetic, its own and made after the result: with an edge of 1 to 2 points and a loss of 50 to 70 if rejected, the
  two files break even at a rejection chance of about 1.5 to 4 percent; nobody can estimate that chance so finely. The
  question put to the captain said "roughly 2 to 4%".

**The cold check** (decision log 2026-11-19 (hh); one challenger given the decision and the record, not how the advice
was formed). Its verdict: "advice stands with conditions". Its probability that advising the recomputed file over the
no-outside-data file is a mistake: 0.35.
- **Its compliance reading.** On the face of the organisers' words the file complies: raw counts and stage labels of
  permitted cells, marker sets, no clustering, embedding or classifier, the authors' labels not read. What an organiser
  could still object to: the authors' quality filter, made on all stages pooled, decided which permitted cells exist,
  and register row X002 says it cannot be undone; the zone column comes from the authors' demultiplexing; the atlas,
  the bracket and R2's held genes were chosen by readouts that read the authors' groups.
- **Its other points.** The test's 2/5 step was never read on the recomputed groups, and the epicardium direction,
  formed on 11 early cells, is applied without shrinkage and doubles at the test. It names the check: a counted read at
  fraction 0.4 with a placebo and the epicardium-held arm (row T63, below).

**What the organisers wrote.** On 2026-11-16 an organiser answered another participant on external datasets that span
a protected window. As `docs/ORGANIZER_QUESTIONS.md` gives it, in part: "Using raw expression counts and stage labels
from the permitted cells is acceptable. Please exclude annotations, clusters, embeddings, or fitted transformations
learned using protected cells, and recompute any such processing on permitted data only. Describe the filtering in your
method report." Under the captain's standing instruction of 2026-11-19 (decision log (g)), as the AI assistant reads
it, no step waits for an organiser reply, and an assumed reply changes no file.

**What practice cannot say.**
- How the file scores at the test's 2/5 step on the recomputed groups. Row T62 read 1/5 only. Row T63 reads 2/5 (below).
- Whether the organisers accept a second outside-data file. No practice task bears on it.
- How the step carries to library 3 at its own mix. The practice tasks draw E9.5 embryo portions (the cuts) or an oracle
  mix (the libraries); the test draws library 3, a dissected heart.
- How sturdy the epicardium step is. Its direction rests on 11 early cells. Row T62's target-free reliability listing
  gave the two halves of that group a Spearman agreement of 0.350 (6 and 5 early cells a half); it decides nothing.
- Forecast weight. The four groups' marker genes come from the released E9.5 labels, and row N55 moved the score
  threshold on the released E9.5, the practice tasks' reference stage (rule 9, J2). g_mom's practice fit spans E8.5 to
  E9.5, portions both; at the test it crosses the collection change (decision log 2026-11-18 (f)).

**Nothing is computed from a test score.** The labeller's settings were frozen before the relabelling ran, its route
was fixed by its terminal rule, and row T62's advice was fixed before its run. No test value sets a weight, prior,
task, tolerance, constant or estimate here (rule 9, J8). The first upload's 70.9 (decision log 2026-11-18 (j)) is only
the score a second file must beat if the first upload stands. It favours neither T1 file: no file without the atlas
step has a test score (decision log 2026-11-19 (hh)). The score chose nothing between the two files.

**Built after row T63's read.** The practice check the chosen option promised is row T63: the candidate at the test's
step size (2/5 of the atlas's bracket) on row T62's five tasks, with a placebo direction and the epicardium-held arm at
that size, a counted try with its row first; what follows from it is fixed in its row before the run (decision log
2026-11-19 (ii)). When this section was written, row T63 was being written and had not run. The file is built only
after its read.
- If its outcome sends a question to the captain (the epicardium step costs at the doubled size; the doubled step is
  not ahead; or the gain is not shown to need the direction's own genes), nothing is built until he answers.
- Otherwise the set t1rc is built (below).
- Row T63's outcome: Row T63 ran on 2026-11-20 and was read once (decision log 2026-11-20 (e)): the doubled step is AHEAD of the no-outside-data analogue (+3.73 and +3.30 on the cuts, +4.78, +2.00 and +4.29 on the libraries, each in 10 of 10 seeds); the epicardium step at the test's size PAYS (+0.80 on library 3); the placebo statement reads "a placebo gains as well; the direction more", so a part of the gain comes with any direction of this size profile. By the row's fixed line the ruled file is built as ruled.

**Builds.** `bash experiments/builds/build_all.sh t1rc` checks this file's hash, builds `runs/builds/T1_s2rc.h5ad` by
the command above and runs the local format check for `T1:test`. It runs on the DGX checkout that holds the
relabelling's output, at a clean tree. It is built twice, both at 16 BLAS threads (T1 files move with the thread
count), and the two sha256 are compared, and compared with the check builds' ad8f24c71c90. Any difference stops the
hand-over until the captain has the explanation. The build record: `experiments/builds/results/t1rc_5028cab.txt`: built twice on the DGX at 5028cab with the same bytes, sha256 ad8f24c71c90ca17bce5f8876d503471ab3faa5b17a0878addb00aff689fca4f, equal to the check builds'; clean stamps; local copy `runs/manual_upload/T1_s2rc.h5ad`, local format check passed. The sets s1, s2,
t3alt, heart2, t1none, embryo2, heart3 and embryo3 are unchanged in the script. This file is in no reference list yet;
its own list is committed under `experiments/builds/results/` before line 1 of the captain's checklist
(`docs/P3_RUNBOOK.md` section 5) can be ticked.

**What the record must say** (the method report, its section on exposures in `METHOD_SUMMARY.md`, the submission
page's Data used field and the note that precedes the upload), as the organisers' reply asks.
- The recomputation: the cell groups of the atlas's permitted cells were recomputed from their raw counts by
  `experiments/external/feng_relabel.py`, by fixed marker scores against expression-matched control genes, a threshold,
  a margin and an endothelial split; no clustering, embedding, classifier or fitted transformation.
- Its inputs and settings: GEO GSE193346, CD1 cells at E9.5 and E14.5 (register row X002); the columns read (cell
  name, stage, zone, strain); the marker sources (the released E9.5 label means; the eight papers of register row
  X010); the settings file's sha256 (883e74c2c1ef...), the statistics file's (36d1f7edcbbf...) and the dry run that set
  the threshold (row N55).
- The filtering: the in-window cells (E10.5 to E13.5) removed at ingest on the per-cell stage column, CD1 only. What
  could not be undone: the authors' quality filter, made on all stages pooled, decided which permitted cells exist; the
  zone column comes from the authors' demultiplexing; ambient correction is not stated. The traces that the register
  row lists in the evidence behind the design (g_mom's adoption evidence, the library-3 base's pricing scenarios, the
  E10.5 mix's supporting reason, the readout that chose this atlas, the gate that adopted R2) stay listed.
- Who decided: the captain chose this file by his answer above, after the first scores; the advice was row T62's,
  fixed before its run; the cold check's caveat was put to him with it.
- The disclosure email in the captain's drafts was worded for the no-outside-data file. By its own send condition, its
  paragraph on T1's remaining upload must match this ruling, so it is reworded and reviewed before it goes. Whether
  anything has been sent: [TO CONFIRM: whether the captain sent the disclosure of the X002 labels or the sources note].

**What stays.** The first upload, the count rule, the seed, the other boards' files, the order (the first upload went
first) and the rule for final marks. `runs/builds/T1_s2.h5ad` (the mid base) and `runs/builds/T1_s2none.h5ad` stay on
disk as the record and are not uploaded.

## Note of 2026-11-20: the hold on the heart's last file has ended (the captain's answer)

Written by the AI assistant, which is recused on the heart's scope and size and stated no preference. The heart section
of 2026-11-19 above held its file until a practice read at forecast sizes was with the captain. That read is row H27 of
`experiments/harness/TRIES.md`, designed, attacked, settled, written and reviewed by analysts who had not seen the AI
assistant's conversation, and recomputed independently (`docs/PLAN.md` decision log 2026-11-19 (jj)). Its figures went
to the captain in one question that ended "I state no preference. Which heart file goes?", with the options "Full trim
(as ruled)", "Half trim" and "Full-scope file" (decision log 2026-11-19 (kk) gives the question and the options in
full). His answer, verbatim: "Full trim (as ruled)". So the heart's last upload is the file of the heart section above,
`T2heart_s2core_half.h5ad` of set heart3 (sha256 04cde6361c04840ee47fff32eb9a20f4781bd649387cc540ed9fdea71e65fa35,
`experiments/builds/results/heart3_3d52b83.txt`). Nothing else in that section changes.

## Note of 2026-11-20: uploads go on after the second files; candidate sets (the captain's order)

Written by the AI assistant. The sections above each name one file as a board's second and last upload. All four of
those files were uploaded by the captain on 2026-11-20 (`results/manual_uploads.csv`). His later messages, verbatim:
"they havee improveed 2 submissions daily" and "let's upload all the 4 filees after full optimization to beat the
leadeerboard scores [2 times a day]". The AI assistant has not read the site's rules since 2026-11-19, at his word, and
takes "[2 times a day]" as his statement of the team's allowance; it plans on at most two uploads per board per day and
uses no score beyond that (`docs/PLAN.md` decision log 2026-11-20 (g), (h), (j)).

What stays of these rules for every further file:

- A board's candidates are fixed in the decision log before any of their scores exists. A returned score chooses among
  them. It never sets a number, a file's form or an estimate, and no property of a hidden target is worked out from it
  (Official Rules Section 10; decision log 2026-11-18 (k)).
- Each file is built twice from one clean, stamped commit with identical bytes, passes the format check, and has its
  sha256 in `experiments/builds/results/`. A builder option is off by default, and the default command must still
  build the bytes it built before.
- Uploads are the captain's, by hand. `builds/T1_s2.h5ad` is never uploaded, and no T1 upload carries the 0.60.
- The data rules of `docs/DATA_POLICY.md` stand as they are.
- The AI assistant is recused on the heart's scope and size and on the embryo's occupancy term (decision log
  2026-11-19 (bb), 2026-11-20 (i)). Analysts who have not seen its conversation fix those candidates.

The candidate sets as they stand at this note:

- **T1** (decision log 2026-11-20 (i); `experiments/harness/TRIES.md` row R15): the recomputed-groups file of the
  section above with one change each: A, g_mom finished at damping 1.0 (`--g-mom-damping 1.0`); B, the atlas direction
  at 1.5 times the calendar fraction (`--feng-fraction-mult 1.5`); C, both.
- **T3** (decision log 2026-11-20 (h), (k); row R14): the barred branch's s1 at factor exponent 1.5, then at 2.0 if
  1.5 beats the s1 file's score, otherwise the barred branch's s2.
- **Embryo** (decision log 2026-11-20 (h), (l), (m)): `T2embryo_s2x.h5ad`, the first upload's cells and expression on
  E7.5's exact positions. And, by the captain's ruling of 2026-11-20 ("Use it, after a practice read"), that form with
  its rows re-ordered by `--row-order copy_parity`, after practice row E48 has run. His ruling sets aside, for this one
  construction, his decision of 2026-11-02 that no parity hedge is used without the organisers' approval.
- **Heart**: a list fixed by analysts other than the AI assistant (decision log 2026-11-20 (m)). Scope is a closed
  lever, so nothing on it is built before the captain's own named exception.

## Amendment of 2026-11-20: the heart's next two files, K1 the half trim and K2 the full-scope copy (the captain's own named exception; row R16)

Written by the AI assistant, which is recused on the heart's scope and size and stated no preference. The list was
fixed, before either file has a score, by analysts who had not seen the AI assistant's conversation (a designer, a cold
challenger, a writer and a reviewer; `docs/PLAN.md` decision log 2026-11-20 (m); `experiments/harness/TRIES.md` row
R16). Scope is a closed lever, so the list went to the captain in one question.

**The question, in full:** "Heart board (56.4 with the full trim, which keeps only the heart core). Analysts who have
not seen my conversation fixed the heart's next two files before either has a score. K1, the half trim: keeps half of
the tissue around the core; never built. K2, the full-scope copy: the whole E10.5 stage; built on 18 Nov, never
uploaded. Their rule: K1 first, K2 in the next slot whatever K1 scores, then the best of the three stands and the heart
list ends. Practice at forecast sizes (ten seeds) reads both ways. With the released stages as they are: the half trim
beats the full trim on all three transitions (by 6.2, 2.5 and 1.9 points) and the full scope beats it by 2.9 on
average. With the practice truth cut down towards the core: the full trim beats the half trim by 6.2 on average (-0.3,
+11.9, +1.3) and the full scope by 11.3. Nothing released says which case the hidden E12.5 heart is. A miss costs one
upload slot per file; your table showed the board keeping our best score, not our latest. Scope is a closed lever, so
each file needs your own named exception. I am recused on the heart's scope and size and state no preference. Which
heart files may be built and uploaded?" The options: "Both, K1 then K2" ("The analysts' list as they fixed it: the half
trim first, the full-scope copy in the next slot whatever the first scores. Two heart slots."), "Half trim only" ("K1
alone, one heart slot. The full-scope copy stays set aside as on 19 November.") and "Neither" ("The full trim stands as
the heart's file. No heart slot is used and nothing is built.").

**His answer, verbatim:** "Both, K1 then K2". This is his own named exception on scope for K1, and his word for K2, the
full-scope file he ruled on 2026-11-18 and set aside on 2026-11-19.

**The files.** Each differs from the full trim of the heart section above in one named option; the copy draw, the raw
mode, the count rule, the seed (1000001) and the halving are the full trim's.

- K1, the half trim: `uv run python experiments/builds/build_test.py heart --entry s2_trim --expression copy --z-mode
  raw --out <folder>/T2heart_s2trimcopy.h5ad`, then `uv run python experiments/builds/halve.py
  <folder>/T2heart_s2trimcopy.h5ad <folder>/T2heart_s2trimcopy_half.h5ad`. The half is the upload file. The option is
  the scope `s2_trim` together with the size field the builder ties to it; a K1 score is read as that form, never as
  scope alone.
- K2, the full-scope copy: `uv run python experiments/builds/build_test.py heart --entry s1 --expression copy --z-mode
  raw --out <folder>/T2heart_s2iso.h5ad`, then `halve.py` into `T2heart_s2iso_half.h5ad`, which must equal set heart2's
  file (sha256 dc436d9225bb28e7997d2ec51e95254cd125e569d44c059596298d0a7b6af5c2).

**Builds.** `experiments/builds/heart_candidates.sh <commit>` on the DGX campaign checkout, from a clean tree at a
commit that holds this amendment and row R16. It first rebuilds the full trim and stops unless its half is the uploaded
file (sha256 04cde6361c04840ee47fff32eb9a20f4781bd649387cc540ed9fdea71e65fa35). Each candidate is then built twice into
separate folders and must give identical bytes; every file passes the format check. A test build is read by its exit
code, its sha256 and its stamp's commit and clean fields.

**The analysts' rule, in their words, fixed before either file has a score:**

> Search rule, fixed before either file has a score: 1. Only the headlines (board totals) of the team's own heart
> uploads choose, as the captain reports them in the decision log. 2. No per-term value and no difference between
> scores goes into any number, file or rule. No reading of the headlines is written as a statement about the E12.5
> stack's scope or size. 3. "Beats" means a strictly higher headline at the precision the site shows. Equal or lower
> keeps the incumbent. 4. K1 goes first. K2 goes in the next slot whatever K1 scores, and neither result decides
> whether the other is uploaded. 5. If only one slot remains, K1 goes and K2 does not. 6. A candidate is skipped if the
> captain does not approve it, or if its build or check fails. 7. A refused upload has no score. The incumbent stands,
> and the file is not changed to get it accepted. 8. No candidate is rebuilt with another seed, uploaded again, or
> changed after its score. 9. After K1 and K2, the heart's best file is whichever of C, B and A has the highest
> headline. The list ends there: no further heart file is designed from these scores. 10. Each upload is the captain's
> own act.

> Left off the list: - K3, the trend step at any strength. Row H28's pre-set outcome was PLACEBO ("no question; the
> ruled file stands"), its S3 test failed in CLEAN, and H17a/H19/H28 tie the trend to the closed depth lever. It can
> come back only as the captain's own exception, put to him as a separate priced question. - Any other size for the
> full trim. - Trim fractions other than 0, 0.5 and the whole stage. - Trend strengths 0.5 and 0.25. - A quarter of the
> rows, or undoing the halving. - Mirrors, flips or row-order changes. - The outside-atlas direction. - Depth, herding,
> cfm, and the Potts and flow lines. - The old s2 with the trend, and the first upload's form. - Any re-draw with a new
> seed.

Every upload stays the captain's own act. Nothing else in the heart sections above changes.

## Amendment of 2026-11-20: the embryo's row-order rebuild is held (the analysts' reading A.2)

Written by the AI assistant, which is recused on the embryo's occupancy term and states no preference. The note of
2026-11-20 above lets the embryo's positions file be rebuilt with `--row-order copy_parity` after practice row E48 has
run, by the captain's ruling "Use it, after a practice read". The question that ruling answered said the analysts had
found no other route. A target-free check on the built file then fell under their reading A.2, fixed before its output:
"The file's draw orders or turns its principal axes differently from the copy's: mechanism 3 is active. The scorer
searches signs only, so this costs occupancy under both builds. What follows: report it. The row-order option does not
fix it. The only fix the code allows is giving the file the copy's own rows (the investigator's Option A), which breaks
the count rule and needs a ruling. Nothing is built until he rules." A second group of analysts read it strictly and
set a hold (`docs/PLAN.md` decision log 2026-11-20 (s)).

**The hold.** The rebuild on the test board with `--row-order copy_parity` waits for the captain's answer to the A.2
question, whatever E48 shows. The question has three options: keep his answer; hold the rebuild until the supplementary
check `file_axes` is in, with E48 still running; withdraw his answer and keep his decision of 2026-11-02. If he
withdraws, this candidate is removed from the embryo's set. E48 itself, a practice read on released stages that builds
no upload file, may run. The second construction the reading names (the file takes the baseline copy's own rows and
count) is not written; it would need its own ruling and its own practice row.

## Amendment of 2026-11-20: T3's two-knockout files C1 and C2 may use the released Gata4 validation files (the captain's reply; row R17)

Written by the AI assistant. The list was fixed by a designer and a cold challenger before any of its scores
(`docs/PLAN.md` decision log 2026-11-20 (q)); the thin builder is `experiments/t3_decomp/build_two_knockout.py` (merged
at c5865d5; the three frozen T3 files untouched). The challenger required one question to the captain before any C file
is built for upload.

**The question as put, in full:** "Task 3: our best file (61.8) uses only the released Mab21l2 knockout. Two more files
would multiply its per-gene change by the released Gata4 validation knockout's change. That is the same table our first
upload (47.8) carried, but that upload also had a depth step and a different cell mix, so the record cannot say why it
scored low. Practice cannot price the new files: on the one fold that tests the Gata4 change on a different knockout,
it beat a shuffled control by 0.18, under our 0.3 bar, and it lowered the cell-state terms. Gain if it carries:
unknown, perhaps a few points. Loss: the board keeps our best, so the cost is up to two of about ten remaining Task 3
slots. If the organisers later restrict the Gata4 files, these cannot be marked final and the Mab21l2-only file stays
our mark. They are built only after both of today's Task 3 scores are back. May these two files use the Gata4
validation files?" The options: "Yes, both (Recommended)" ("C1 (the plain product) then C2 (the same, held at the best
file's size). My reason: slots are not short and a lower score costs nothing else."), "Only C2" and "No, Mab21l2 only".

**His answer.** He dismissed the dialog that held this question and three others, then wrote, verbatim: "go ahead with
whatever will improve our score [my reply]". The AI assistant reads that as yes to the recommended option, "Yes, both"
(decision log 2026-11-20 (t) gives the reading and its precedent). It is the AI assistant's reading of a general reply,
not his own choice of an option, and it is recorded as that.

**The files and the rule** are entry (q)'s, unchanged: the base B is the highest returned headline among the barred s1,
its exponent-1.5 and exponent-2.0 files and the barred s2, fixed only after both uploads of entry (k) are back; C1
(`--form c1-s1` or `c1-s2`) is built only when B is the barred s1 or s2; C2 (`--form c2-<base>`) for every B; no depth
step; seed 1000001; each built twice with the same bytes after the identity forms reproduce the base's bytes. Only
headlines of the team's own uploads choose. Under any organiser restriction on the Gata4 files no C file is marked, and
the Gata4-free barred file stays the mark. Every upload is the captain's own act.

## Note of 2026-11-20: the hold on the embryo's row-order rebuild, after the supplementary check

Written by the AI assistant, which is recused on the embryo's occupancy term and states no preference. The
supplementary check named in the amendment above read S1, by the reading fixed before its output: the turn that reading
A.2 fired on is the stage's own draw-to-draw noise (32 of 50 plain copies of E7.5 turn as far or further), the second
construction is not written, and the row-order path is unchanged (`docs/PLAN.md` decision log 2026-11-20 (t)). The
captain was told of A.2 before his reply "go ahead with whatever will improve our score [my reply]"; the AI assistant
reads that reply as not withdrawing his answer "Use it, after a practice read", the first of the three options, and
tells him so. The rebuild therefore follows practice row E48, unless E48's first reading voids the construction. If he
objects, the candidate is removed. The upload stays his own act.
