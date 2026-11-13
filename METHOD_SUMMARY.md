# Method summary (running draft)

Required for prize eligibility: winners must deliver a written report within 14 days of notification (approach, design
decisions that mattered, what was trained on what, enough detail to understand why it performs as it does). The
Generality Award is judged on this summary. Keep it current as methods change; every claim should point to a run in
`results/registry.csv` or an experiment report in `experiments/*/REPORT.md`.

Status 2026-11-10: the entry for each board is chosen (below). The test-board specs were taken in on 2026-11-04 by
`vec boards update --from-file` (634b5fd; nobody opened the saved file, but the board index had been displayed and
pasted once before it was saved: Data provenance, below). The current ten test-board files (s1, s2 and the Task 3
switch files) are the set in `runs/archive/builds_f337dc6/SHA256SUMS`: T1's s1 and s2 were rebuilt at f337dc6 with R3b,
and the other eight are byte-identical to the builds of ae6d6b3. All ten pass `vec check`. The eight unchanged files
passed the server's format check at ae6d6b3; T1's two rebuilt files have none on record (`results/submissions.csv`),
since the server format checks wait for the freeze morning, on the captain's watch. No scored submission has been made. The final phase opened on 2026-10-27, and the team unlocked the P2 validation answers on 2026-10-28 (Data
provenance, below). A change entered an entry in one of three ways: it passed the adoption rule on
pseudo-tasks built from released stages (Shared machinery), a pre-registered rule decided it, or the captain ruled it
in where no counted run could judge it as built (T1's tissue mix, Feng's direction on the pre-registered readout, R1,
R2 and R3b; each is named below).
Since 2026-11-04 the summary is framed by the team captain's guiding principle: each entry is described by the change
it predicts beyond the published baseline and by the held-out evidence for that change ("What each entry predicts
beyond the published baselines", below).

Status 2026-11-13: no entry has changed since 2026-11-10, and no scored submission has been made. The campaign the
captain asked for on 2026-11-12 ran four waves of tries on released data (`docs/PLAN.md` decision log 2026-11-12 (a)
to 2026-11-13 (c)) and kept nothing. What it found is under "Tries ledger" below. What it put to the captain (the
embryo's second upload, the base of T1's first upload, and two options that their pre-set rules closed) changes
nothing until the captain rules.

Entries changed on 2026-11-13 (later that day; `docs/S2_RULES.md`, "Amendments of 2026-11-13"; rows R4 to R6 of
`experiments/harness/TRIES.md`; `docs/PLAN.md` decision log 2026-11-13 (d)). None of the three passed the adoption
rule; each entered by the captain's ruling or delegation, after the captain was asked about each one with its
numbers. No scored submission has been made. The ten test-board files as they now stand are the set in
`runs/archive/builds_fd2683a/SHA256SUMS` on the build machine (rebuilt there twice on 2026-11-13, byte-identical,
all passing `vec check`; `experiments/builds/results/rehearsal_fd2683a.txt`): three new files (T1's s1, the embryo's
s2, T3's second upload) and seven byte-identical to the set of f337dc6 named above. The three new files and T1's s2
have no server format check on record yet.
- **Task 2 embryo, s2:** `fgw_aligned_series` in place of fgw_shape0 (the captain's choice among three options). E7.5's
  cells are rendered into E8.0's per-gene detection profile and depth, and the moves come from the aligned coupling.
  It pays if E7.75 was measured like E8.0 (+3.4 to +6.7 as better of two on the practice tasks of that kind) and gives
  up E7.5's shape hedge (-1.42 on the counted E7.5 backtest).
- **Task 3, second upload:** the refused branch's s1 also carries a per-gene factor fitted on released wild-type
  samples alone (WT E8.75's own deviation from a trend through the other wild-type stages), applied per lineage
  group at exponent 0.5 (the captain's ruling, which reopens the ruling of 2026-11-01 for this file). On this file's
  form it adds +4.72 and +7.12 on the two Gata4 replicates; the gain is the severity gate cleared by a thin margin.
- **Task 1, s1:** the cells are drawn from E10.5's dissected-heart library alone (library 3) at its own chamber mix,
  instead of from all three E10.5 libraries; every other step is unchanged, and s2 keeps the pooled base at 60%
  cardiomyocytes. The captain left the direction to the AI assistant ("do whatever it takes to top the leaderboard"),
  which decided from released data alone: the hidden sample is a dissected heart, as library 3 is; five of the seven
  released libraries are at least as deep as library 3; and across one real day library 3's base was ahead of the
  pooled one in 10 of 10 seeds. It is a depth choice in effect (+27% median library), and it loses if E12.5 is
  shallower than about 0.7 of library 3.

Status 2026-11-13, after the campaign's fifth round (`docs/PLAN.md` decision log 2026-11-13 (g) to (j)): no entry has
changed and no scored submission has been made. The tries stand at T1 30, embryo 36, heart 4, T3 15. Two pricings on
the zero-time-gap grid of try T30 looked at Task 1's second upload (tries T32 and T33; both recomputed
independently). Drawn from library 3's cells alone it gains where the hidden sample is at least as deep as library 3
and loses where it is shallower, so its rule recommended nothing. Drawn from a base half-way between the pooled
base and library 3 it leads today's pair under all nine priors of the grid (+2.43 under the main one), and its rule
puts it to the captain; its lead comes from the deeper scenarios, and nothing changes without the captain's ruling.
On the heart, the within-type step that Task 1 carries (R3b) was written for the heart head and stopped by a
target-free gate before any scored run (diagnostic N26). A forecast of the total as a distribution
(`experiments/gate/forecast.py`) is on record beside the gate's figure, which is unchanged.

## Team and track

Team Mithril, Track 1 (Human Team). A person reviews results and makes every submission decision; coding agents do
the work.

## What each entry predicts beyond the published baselines

The team captain's guiding principle (2026-11-04) treats the embryo as a process of becoming: an entry predicts which
cell states emerge, in what proportions, how their expression changes and, on Task 2, how cells organise in space.
The published baselines predict no change in any of them. copy_last submits the reference stage unchanged (expected: E10.5 on Task 1 and the heart, E7.5, the earlier
bracket, on the embryo; the test boards publish no reference stage, A3), wt_identity the matched wild type (Task 3), and
each scores 50. Below, per board and facet,
are the change an entry predicts beyond its baseline and the held-out evidence for it: backtests on released stages
or knockouts the change was not tuned on, over 10 seeds. Harness rows are scored with the calibrated emulator (bake
protocol), in board points. The Task 3 rows come from the pre-registered decomposition's own scorer (2026-10-28: the
emulator's metrics with its default draws and computed anchors, before the 2026-10-31 calibration), the Task 1 readout
is raw metric deltas, and the T1 and heart s2 rest on sensitivity analyses with an oracle base or core size. Where no
backtest exists, the cell names the nearest check and what it is: in-sample, selected-on, an oracle base or one knockout. The
per-group numbers are in `experiments/harness/TRIES.md` ("Where the entries stand against the published
baselines"); T1's lines there were built before R1, R2 and R3b and have not been rerun for them (`docs/PLAN.md` decision
log 2026-11-10 (q)). In this section, D-numbers are the pre-registered Task 3 rules (`experiments/t3_decomp/rules.yaml`),
C-numbers claims in `experiments/coords/REPORT.md`, A-numbers rows of `ASSUMPTIONS.md`, and E11, K2a, R1, R2 and R3b rows
of `TRIES.md` (R2 here is T1's hold in Feng's direction, not the embryo request of the same name declined on 2026-11-05).
B0 is the harness's baseline candidate (`vec.harness`, feat-009): on Task 1 a chamber-only copy of the reference with
damped momentum, on Task 2 a copy with damped momentum or bracket interpolation, both in count space.

**Task 1, E12.5, beyond E10.5 verbatim.** s1 = M1 with R1 (cells without a Feng group keep their E10.5 expression; the
captain's ruling of 2026-11-05), with R2's two gene modules (the captain's ruling of 2026-11-06) held flat in Feng's direction and, by the captain's ruling of
2026-11-10 (R3b), within-type factors on cardiomyocytes and endocardium; s2 = the same with cardiomyocytes at 60% of
heart cells. (Since 2026-11-13, row R6: s1 draws these cells from library 3 alone at its own chamber mix; s2 is
unchanged, on the pooled base.)

| Facet | Change beyond the baseline | Held-out evidence |
|---|---|---|
| States that emerge | None. Births by the legal route (states present in legal post-E14.5 chamber data and absent at E10.5) were designed and pruned on 2026-11-04 (`experiments/harness/TRIES.md`) | - |
| Proportions | Chamber categories only (the surrounding, cardiac-progenitor (first and second heart field, JCF, DMP), pericardial and pharyngeal categories dropped; `configs/t1_chamber_classes.yaml`), at E10.5's own chamber mix; s2 moves cardiomyocytes to 60% of heart cells | A premise rather than a measurement (the organisers' E12.5 collection note). On the chamber-cut backtest (E10.5 from E9.5) B0 beats a copy (+1.62, 9 of 10 seeds) but not the floor (-0.32, vetoed), and that backtest draws E9.5's mix, so E10.5's own mix is never tested (`experiments/harness/results/d2_v1/`). The s2 rests on a sensitivity grid on an oracle base, not on a backtest (`experiments/s2_hedge/`) |
| Expression | Count-space shifts: Feng 2022's E9.5 -> E14.5 direction at 2/5 for cardiomyocytes, endocardium, endothelium and epicardium; the other cells (mostly blood) keep their E10.5 expression (R1, the captain's ruling of 2026-11-05; until then half the E9.5 -> E10.5 per-type change, a step with no measured sign skill on blood); 23 genes held flat (13 Y-linked, Xist and Tsix, 2 named strain genes and 6 by the pool rule alone; `experiments/t1_direction/results/held_flat_genes.csv`); R2 holds 40 more flat in Feng's direction only (32 dissociation-stress and 8 haemoglobin genes; `experiments/t1_hygiene/results/n12_hygiene_size/`); R3b adds a per-metacell, category-mean-free within-type change on cardiomyocytes and endocardium, fitted from E9.5 to E10.5 and not scaled to the two-day step (below) | The pre-registered readout on E10.5 from E9.5: de_score +0.175 and de_direction +0.188 over cross-fitted momentum, cardiomyocyte sign agreement 0.912 (`experiments/t1_direction/results/E10.5/`). A labelled diagnostic at 1/5, in-sample for that readout: +0.97 over the floor and +1.29 over B0, 10 of 10 seeds each (measured before R1, R2 and R3b); over the floor all of the gain is on de_direction (skill 0.620, against 0.536 for a copy), and de_score (0.480) stays below the floor (`experiments/t1_m1/results/diag_v1/`). No harness task can judge it, and the 2/5 used at test is untested. R2 has a target-free size check only (N12: 0.7-1.1% of each group's summed absolute log factor) and no run of its own; R1 and R2 together read -0.16 and -0.14 in T26's void, in-sample diagnostic. R3b's evidence is under Task 1 below: one fit, cuts it was selected on |
| Spatial organisation | Not scored on Task 1 | - |

**Task 2 heart, E12.5, beyond E10.5 with its coordinates.** s1 = H1 (B0 with the pre-registered size); s2 = H2 as
the half trim.

| Facet | Change beyond the baseline | Held-out evidence |
|---|---|---|
| States that emerge, proportions | None: E10.5's composition is held (the trend was dropped on 2026-11-01). s2 removes surrounding cells, farthest from the core first, to half the in-plane scope | s2 against s1 on the released targets: -5.88 on E10.5 from E9.5 (ahead in 0 of 10 seeds), +1.75 on E8.5 from E8.25_late (9 of 10), +1.57 on E10.5 from <= E8.75 (10 of 10), every prediction scaled to an oracle core size (`docs/S2_RULES.md`, `experiments/s2_hedge/`) |
| Expression | Half the E9.5 -> E10.5 change of each harmonised lineage, in count space | E10.5 from E9.5 and from <= E8.75: +0.47 and +0.90 over a copy (9 and 10 of 10 seeds), +0.43 and +1.45 over the floor (7 and 10 of 10); every metric within 1.3 points of the floor (`experiments/harness/results/d2_v1/`). The expression-change group itself moves +0.77 and -0.03 against a copy (+0.38 and -0.46 against the floor); on the 1.75-day backtest the gain over a copy sits on neighbourhoods (+0.50) and cell state (+0.43), not on the group the momentum's reason targets (`TRIES.md` rule 2; `experiments/harness/results/headroom/`) |
| Spatial organisation | E10.5's tissue and shape; the core grown by 1.58 x the E9.5 -> E10.5 log growth (a fixed three-rule average), xy grown and z held | None is possible. No released pseudo-target reaches the E10.5-E12.5 regime, and on E9.5 every growth rule scores like copy_last (`experiments/coords/REPORT.md`, C7; A8). The one released scope check (n = 1) ranked s1's scope, which equals H0 as built, last of three (C10) |

**Task 2 embryo, E7.75, beyond E7.5 with its coordinates.** s1 = fgw; s2 = fgw_shape0 until 2026-11-13 and
`fgw_aligned_series` since (row R4, the captain's ruling). On the two backtests
that chose it as the entry (selected on, not tuned on), fgw gains +4.15 against the organisers' pseudobulk_shift on E7.5
and +4.92 against B0 on E7.25, ahead in 10 of 10 seeds each. After the organisers' protocol text showed E7.75 collected like
E8.0, tries E16 and E17 (2026-11-11) tested whether s2 should take E8.0's geometry instead: on a backtest that copies that
sampling structure, neither geometry beat fgw in any seed, so s2 is unchanged.

| Facet | Change beyond the baseline | Held-out evidence |
|---|---|---|
| States that emerge, proportions | Real cells, half from E7.5 and half from E8.0: E8.0's states enter at half their E8.0 share, and states seen only at E7.5 halve | Mixed. On E7.5 (t = 1/3) mmd_u skill 0.671 and variogram 0.664; on E7.25 (t = 2/3) 0.497 and 0.283, both below the floor of 0.5 and the variogram well below it. E7.25 sits much nearer E6.75 than its day says, and a mixture at its nominal fraction does worse there than a copy of E6.75 (`experiments/harness/results/diag_emb_mix/`). Nothing measured places E7.75 between its brackets (the clock warp was pruned) |
| Expression | The pseudobulk moves between the brackets; the cells are not edited | de_score 0.581 and 0.683, de_direction 0.703 and 0.702, against 0.450 and 0.522, and 0.526 and 0.628, for a copy (`experiments/harness/results/emb_fgw_v1/`) |
| Spatial organisation | s1: cells moved half-way along their FGW partners after a proper rotation; s2: E7.5's geometry. Both at the size given by PCHIP of log RMS radius through the released embryo stages, with the 15-NN neighbourhoods of real cells | Size: scale skill 0.727 and 0.624. Shape, s1: d2 0.372 and 0.717, Dice 0.385 and 0.595 (below a copy on E7.5); at test the moved cloud is flatter than either bracket in every draw (mirrored partners; `experiments/harness/results/diag_emb_coupling/`), and the fix, try E11, gained on E7.5 and was vetoed on E7.25. Shape, s2: d2 0.439 and 0.383, Dice 0.494 and 0.434. Neighbourhoods: 0.640 and 0.602 (s1), 0.654 and 0.563 (s2). Better of the two: 60.87, against 59.95 for fgw alone (`experiments/harness/results/emb_shape_v1/`) |

**Task 3, beta-catenin knockout at E8.75, beyond the unchanged wild type.** s1 = the allowed branch's s1; s2 = the
refused branch's s1, uploaded first. If the organisers say before that upload that a floored severity scores zero, s2 is
instead the barred branch's s1: the Mab21l2 factor with no depth step, Gata4-free (the captain's ruling of 2026-11-06,
`docs/PLAN.md` decision log 2026-11-06 (b)).

| Facet | Change beyond the baseline | Held-out evidence |
|---|---|---|
| States that emerge | None: the scenario only depletes | - |
| Proportions | Both: OFT/RV-CM, aSHF, pSHF, aPHM and pPHM depleted, first-heart-field and left-ventricle states spared; to 0.38x in s1 (the Gata4 replicates' median shift, D7), 0.61x in s2 (Mab21l2's, D7m) | None for beta-catenin. On Gata4 the term adds +0.22 and -0.43 on top of the borrowed factor (D8 did not fire). The same states at 0.5x on WT cells score 53.62 and 48.13, against 52.39 and 50.61 for WT cells, with severity at its floor in 4 and 10 of 10 seeds. The ledger's proxy for this composition (A11: de_direction > 0.2 on both Gata4 replicates, applied blind) failed, 0.113 and 0.022, so A11 is killed (2026-11-05) |
| Expression (the per-gene response) | s1: the Gata4 per-gene factor, measured in the compartments the knockout cannot touch (non-Mesp1) and applied in count space to every WT E8.75 cell, with Gata4's per-type depth. s2: none | Within one knockout only: the other replicate's factor scores 68.05 and 64.09, with depth 74.50 and 72.15 (rules D2 and D4; Gata4's shallow depth is a property of its imaging run, and whether the beta-catenin files share it is open and cannot be settled from released data, A23; the depth step carries 43% of s1's dp energy). One unseen knockout has tested this factor since (try K2a, 2026-11-04, on the Gata4 -> Mab21l2 fold): no support, +0.18 over its gene-permuted control (ahead in 9 of 10 seeds, short of the adoption rule; the miss rests on one seed's severity gate) and -0.65 against our WT copy (ahead in 4 of 10 seeds, a veto; on the nine seeds other than seed 3 still -0.11, ahead in 4 of 9); that fold is harsher than the test (`experiments/harness/TRIES.md`). The one factor carried across knockouts, Mab21l2's on Gata4 (D11b), scored 59.40 and 62.15 against 52.39 and 50.61 for WT cells, under its pre-registered bar of 60 on rep1; its gain is on severity (skill 0.893 and 0.892), and its de_score stays just below the floor (skill 0.490 and 0.498; `experiments/t3_decomp/results/gata4/`). Against Mab21l2's own WT E9.5 the same knockout's factor scores 52.04 and 46.24, below plain WT cells (52.39 and 50.61), and the review of 2026-11-05 read most of its carry to Gata4 as the WT E9.5-against-E8.75 offset, which s1's factor shares (`experiments/gate/results/t3_swap_review_2026-11-05.md`); with its component along a WT developmental axis removed it scores 55.67 over the 20 rows, against 60.77 with it and 49.14 for the E9.5-referenced factor (U3; no rule changed). On a WT pool and reference thinned to each replicate's depth, s1's stand-in (the Gata4 factor and depth) scores 53.89 with severity at its floor in 20 of 20 rows, and the Mab21l2 factor 59.39 with 5 of 20 at the floor (U3d, `experiments/t3_decomp/results/orth_diag/`). A10, that the factor carries across knockouts, was killed on 2026-11-05 |
| Spatial organisation | Not scored on Task 3 | - |

s2 is floor insurance by design: it carries no number derived from the Gata4 files (though its having no factor was
itself chosen on Gata4 scores; `docs/S2_RULES.md`, T3) and stays near wt_identity, so the better of the two stays near
the floor if s1's factor does not carry to beta-catenin or has the wrong sign (A10, killed on 2026-11-05). Its
composition is the same bet as s1's (A11, killed on 2026-11-05), at a smaller magnitude.

**What no entry predicts yet.** No Task 1, heart or Task 3 entry predicts a state that no released stage has: Task 1
births by the legal route were pruned on 2026-11-04, the heart holds E10.5's composition, and the Task 3 scenario only
depletes. On the embryo, states that E8.0 has and E7.5 lacks enter at half their E8.0 share, as E8.0's real cells, by
interpolation only; no backtest checks that at t = 1/2. The heart entry's checked changes stay within 1.3 points of
copy_last on every metric, so its gain at test rests on the pre-registered size, which no released stage can test.
Task 3's s2 is the wild type plus a five-state depletion, by design.

## Data provenance and compliance (applies to every task)

- **Data credit.** The challenge data are unpublished (Dr. Qingquan Zhang and Dr. Neil Chi, UCSD), released for research
  and education within the Virtual Embryo Challenge, which this summary acknowledges. Any publication or public
  presentation of this work needs the prior written approval of Dr. Neil Chi's group (Terms of Use clause 14).
- **Released competition data only, plus disclosed legal external data** (`EXTERNAL_DATA.md`). Nothing measured inside a
  protected window is used, by any route; windows and the do-not-use list are in `docs/DATA_POLICY.md`.
- **Validation answers (since 2026-10-28).** The P2 validation answers (T1 E10.5, T2 heart E8.5 and E10.5, T2 embryo
  E7.5, the two Gata4 knockout replicates) are opt-in in the final phase; the team captain unlocked them on the Data
  page on 2026-10-28. Rules s.10 makes them an exception to the held-out-data restrictions solely for preparing P3
  predictions; the team reads that as including their use as model inputs at prediction time (for example E10.5 as
  the starting stage for E12.5) and asked the organisers to confirm (21 and 28 October), with no reply by the last screened sweep on record (2026-11-05, `progress.md`).
  The team proceeds on that reading and on the predicted organiser answer (`docs/ORGANIZER_QUESTIONS.md`, working
  answers), not on a ruling; the Task 3 use of the Gata4 replicates rests on the same reading. External data from
  those stages and genotypes stays out, and no test inputs are distributed.
- **No board anchors and no leaderboard scores set any number in a prediction or a method choice.** Test-board anchors
  are stripped in memory when the board specs are refreshed (`vec.boards.update_boards`; `update_boards_from_file` does
  the same for a board index the captain saved in the browser, the route the captain cleared on 2026-11-01 and used on
  2026-11-04, 634b5fd). Since the final phase opened the public board file has been opened once on record: on 2026-11-04
  the captain displayed the board index in a browser and pasted it, and then the Data page, into the AI assistant's
  session. The paste held no test anchors and no truth counts but did show the test boards' cell limits, and the Data
  page lists hidden targets' file sizes. No value from either is written in the repository, its notes or its commits,
  and none is used (the paste stayed in the AI assistant's local session log, outside every repository): every use of a
  test board's limits is code fixed before that day, and since then no try may change a cell number, a size or a scope
  on any board, or T3's composition, and T1's composition only by the legal route, unless the captain rules otherwise
  (`docs/DATA_POLICY.md` 2026-11-04 and (b) to (d)). The organisers were told the same day
  (`docs/ORGANIZER_QUESTIONS.md`). The captain's upload gate takes its figure from an overall leaderboard total and
  decides only when to upload (`docs/PLAN.md` section 8).
- **Task 3 gene directions come only from released data** through committed code (organiser ruling of 2026-10-15): no
  gene-level knockout findings from the literature, no curated gene sets built from perturbing the held-out genes.
- **The leaked held-out heart size** (veckit PyPI 0.1.2 docstring) was reported to the organisers and never used; the
  heart size comes from a procedure fixed before the release.
- **Released files contain synthetic records** (copies of cells in adjacent sections and averaged, half-integer cells);
  prediction inputs drop them (`vec.prep.load_stage(drop_duplicates=True)`), the local scorer loads targets as released.
- Every exposure to a held-out value during research is recorded, with what was done about it, in `docs/DATA_POLICY.md`.

## Shared machinery

- **Local board emulator** (`vec.evaluate`): the official veckit metric code (git `cc3859a`), the published skill map
  and weights, the T2 structure gate, 10% working copies and split-half ceilings, several seeds. Two scorer builds:
  `main` for comparing candidates, `prefix` (veckit 0.1.1 canonicalisation) to reproduce live scores.
- **Real cells, count-space edits:** predictions are built from real released cells; expression changes are applied by
  binomial thinning or Poisson top-up on recovered integer counts, then renormalised to log1p CP10k, or not at all (the
  embryo entry places unedited real cells). No decoded means. The captain allowed additive log-space edits on the
  embryo board (2026-11-01), where the organisers' pseudobulk_shift beat our count-space B0; the chosen embryo entry
  does not use them.
- **Label harmonisation** (`vec.harmonise`): logistic regression on 50 PCs of the shared panel, group probabilities
  summed from fine labels, depth matching; cross-file compositions quoted as a top-call/EM bracket
  (`experiments/harmonise/REPORT.md`).
- **Pseudo-task harness** (`vec.harness`): paired backtests on released stages in each test board's geometry, scored
  with the calibrated emulator (the server's own draws reproduced: all 66 validation anchors exact and four of the five organiser baseline rows to their printed precision; the heart-extrapolation pseudobulk_shift row is not reproduced, its cell-type matching being unpublished, A2), the floor checked
  to score 50 on every seed, candidates built from the allowed stages only (no file paths, a read guard against files
  derived from targets) on shared draws, and a change the harness can judge kept only under the adoption rule (approved rule 7's test per
  pseudo-task); what it cannot judge entered by a pre-registered rule or the captain's ruling (Status, above).
- **Pre-registration:** the size procedure, the T3 decomposition rules, the T1 direction readout, the cell-count rule
  and the label-mapping rule were committed before the data they act on. The other steps were not: B0's momentum and fgw
  were chosen on backtests of released stages, T1's tissue mix is a premise, and R1, R2 and R3b came by the captain's
  ruling after E10.5 was unlocked on 2026-10-28, each with its row committed before its build. The s2 rules
  (`docs/S2_RULES.md`) were written by the AI assistant under the captain's delegation of 2026-11-02 from sensitivity
  analyses on released data, before any test score, and frozen on 2026-11-03. The heart rule was amended after its
  released-data outcome, and the amendment is logged as such.
- **Tries ledger** (since 2026-11-03, `experiments/harness/TRIES.md`): every change tried for an entry is written
  down before its run with its physical reason, the metrics where its gain should show and what would make a gain a
  fit to one backtest, and every try is counted (running totals on 2026-11-10: T1 23, embryo 11, heart 3,
  T3 6; T26 is counted although its reproduction gate voided it). One counted try had no row before it: T1's seventh, a
  campaign designer's unregistered in-sample comparison against E10.5, counted under rule 4. Tried and rejected on 2026-11-03: matching detection depth by
  binomial thinning (released stages differ several-fold in counts per cell: the embryo series rises from 30 to 212,
  the heart series falls from 233 to 56). It lost on one backtest of each board (-4.18 on E7.5, -3.92 on the 1-day heart
  backtest, both ahead in 0 of 10 seeds), lowering mmd_u and neighbourhood MMD there, and gained on the other (+0.74 on
  E7.25, +0.45 on the 1.75-day heart backtest); the adoption rule vetoed it, so it is not used. Blending coupled pairs of cells also lost on both embryo backtests:
  averaged cells are too narrow a population. A target-free diagnostic found that the FGW coupling cannot tell the
  embryo's left side from its right (expression is left-right symmetric at these stages), so its partners can be
  mirrored and the moved cloud flattens; it happens in every draw of the test pair. Taking the moves from an exact
  transport on expression and position in the rotated common frame fixed every draw and gained +2.13 on the E7.5
  backtest (shape only) but lost -1.59 on E7.25, where its shorter moves cost shape, so it is not used; the s2
  (E7.5's geometry) hedges the flattened s1 geometry. Since 2026-11-04 every row also names its board's published
  baseline, the change it predicts beyond the baseline (facet and metrics), the scored group it targets and the
  held-out evidence that would show its effect, and tries go first to the weakest group with an open lever (rules 5-8
  and the headroom tables, built from committed results by `experiments/harness/headroom.py`). The captain's campaign of
  2026-11-09 lifted the ledger's limits on research effort (the prune of candidates expected to gain under 0.3 points,
  "no further try" and the lines the ledger had closed), not the data rules, the row before each run, the try count, the
  adoption rule or the captain's ruling for any change to an upload (`AGENTS.md`, rule 8). Its outcome by 2026-11-10: T11
  passed the adoption rule and its controls and entered T1 as R3b by ruling; E12 and E15 were run and not kept; T8 to
  T10 stopped the within-type GPU line; U3 and U3d changed no T3 rule; the late-proxy tournament (T13 to T15) failed its
  gate; T1's strong-null bet was measured (T12a and T12b) and left out by the captain; H6, E13, E14, K4 and K4c were
  stopped by their target-free gates and not counted (`experiments/harness/TRIES.md`; `docs/PLAN.md` decision log
  2026-11-10).
  The campaign of 2026-11-12 and 2026-11-13 (running totals on 2026-11-13: T1 28, embryo 36, heart 4, T3 15; one
  reading, N25, ran before its row and is counted). Its finding: released samples differ more by how they were
  measured than by stage. On the embryo ladder some neighbouring samples share a per-gene detection profile and others
  do not, and the seven released T1 libraries span 4,030 to 14,755 counts per chamber cell in two kinds of sample.
  Tried and not kept on the embryo: every cell moved to one depth (E19), a count-model latent interpolation (E21,
  scored this time) and a wider top-up neighbourhood (E27). Priced as second uploads on nine practice tasks and put to
  the captain (request R-A): the shallower bracket's cells rendered into the deeper bracket's run profile (E24, E25,
  E28). On T1, a base drawn from the dissected-heart library alone was priced at zero time gap and across one real day
  (T27 to T31); it is a depth choice in effect and is not adopted. On T3, a per-gene factor estimated from wild-type
  samples alone, the reference's own deviation, gained on both knockout runs that can be read (applied per lineage
  group +4.7 to +7.1 against WT E8.75 and +20.0 against WT E9.5; pooled +1.2 to +5.1 and +20.0) but failed its
  pre-set reading by one comparison of nine, so it is not used (K7, K8).
- **Provenance:** every prediction file carries a build stamp (its sha256, the commit that built it, and whether the
  code was unchanged from the build's start to the write); an official upload is accepted only with a matching clean stamp
  (a format check records the stamp and warns when it is missing), and the submission ledger records the building commit
  of every file sent. T1's two files built with R3b depend on the BLAS thread count: at 16 threads, the default of
  `experiments/builds/build_all.sh`, a repeat is byte-identical; at 8 threads 0.41 of s1's rows differ, in 0.043 of its
  entries (`docs/PLAN.md` decision log 2026-11-10 (q); the cause is expected to be floating-point order in the fit's
  PCA, k-means or Sinkhorn steps and was not traced further). The other eight files were byte-identical at 1, 4, 8 and
  16 threads when last checked (`experiments/builds/results/determinism_6964b7d.txt`) and have not changed since. The
  two candidate second uploads of request R-A, which are no entries, reproduce at 16 threads on the local machine
  (`experiments/builds/results/determinism_1f8956b.txt`). On 2026-11-13 the ten files as amended were built three
  times on the build machine at 16 threads with identical bytes (`experiments/builds/results/rehearsal_fd2683a.txt`).
  Between two machines with different processors seven of the ten are byte-identical; the three T3 files that apply
  a per-gene factor hold the same cells in the same order, with the expression of a few cells drawn differently (3,
  84 and 14 cells), because a nearest-neighbour smoothing step follows the linear-algebra library's processor-specific
  kernel. The hashes on record are the build machine's.

## Task 1: temporal (E12.5) — M1

- **Scope.** E12.5 is a heart-chamber-only dissection, so the entry draws real E10.5 cells to E10.5's own composition
  restricted to the chamber categories (`configs/t1_chamber_classes.yaml`): blood 24% of chamber cells, cardiomyocytes
  26.5% of heart cells (the captain's choices of 2026-11-01). E10.5 is the same assay (nuclei, by the organisers' Data-page text; our own check on E9.5, N13, was unclear) on dissected hearts; the
  earlier stages were embryo portions, and their compositions do not transfer (blood 2.6% of chamber cells at E9.5; cardiomyocytes 66% of heart cells). Legal external data do not settle the cardiomyocyte share: Feng's own share barely moves from E9.5 to
  E14.5, but the E14.5 shares of Feng's two strains are 41.5% and 76.6% of heart cells (X002), which is why the s2
  hedges it. The s2 hedge: cardiomyocytes at 60% of heart cells, decided from the T1 grid
  (`docs/S2_RULES.md`).
- **Direction.** Cells of the groups Feng 2022 profiles (cardiomyocytes, endocardium, endothelium, epicardium) move in
  count space along Feng's within-dataset E9.5 -> E14.5 log ratio of mean CP10k, scaled by the step's share of that
  bracket (E10.5 -> E12.5: 2/5); the other cells (blood, mostly) keep their E10.5 expression (R1, the captain's ruling
  of 2026-11-05: the half step of E9.5 -> E10.5 momentum they took before had no measured sign skill on blood and spans
  T1's sampling change). A pre-registered readout on the E9.5 -> E10.5 transition chose Feng over cross-fitted momentum
  (de_score +0.175, de_direction +0.188, cardiomyocyte sign agreement 0.91; `experiments/t1_direction/results/E10.5/`),
  and the captain ruled that this check governs the direction, since no T1 pseudo-task can judge it (Feng's E9.5 cells
  lie inside each one's external window). A labelled end-to-end diagnostic in the readout's geometry, in-sample for that
  readout: +1.29 and +1.09 board points over B0, 10 of 10 seeds, measured before R1, R2 and R3b
  (`experiments/t1_m1/results/diag_v1/`). Sex-linked, Y-linked
  and strain-variable genes are held flat (a list of 23 genes fixed from E8.5 and E9.5 before E10.5 was released). R2 (the
  captain's ruling of 2026-11-06, built at 02abf39; not the embryo request also called R2, declined on 2026-11-05) holds
  40 more genes flat in Feng's direction only, for the cells M1 moves along it: module D, the 32 dissociation-stress
  genes of van den Brink 2017 Supplementary Table 1 (X007), and module A1, the 8 haemoglobin genes matching ^Hb[ab]-
  (chosen on the statement in X009). R2's premise, which no run tests, is that both are technical features of Feng's
  whole-cell data; Feng's bracket also spans the embryonic-to-adult globin switch. Two further modules (mitochondrially
  encoded genes and nucleus-retained ncRNAs) were dropped before R2's size check, because a target-free check (N13) could
  not confirm that T1 is nuclei. No backtest can judge R2 (Feng's E9.5 lies inside every T1 window), so it entered by
  ruling on a target-free size gate (N12, `experiments/t1_hygiene/results/n12_hygiene_size/`): the two modules hold
  0.7-1.1% of each group's summed absolute log factor, and the direction with and without them correlates at Spearman
  0.998. No run isolates R2. The one score that contains it, T26's labelled in-sample diagnostic (void by its reproduction
  gate, at 1/5 of Feng's step, with R1 in the same arm), has M1 with R1 and R2 at -0.16 on the chamber cut and -0.14 on the
  whole cut against M1 without them, ahead in 0 of 10 seeds on each (`experiments/t1_m1/results/t26_m1_gmom/`).
- **Within-type states (R3b, the captain's ruling of 2026-11-10).** On cardiomyocytes and endocardium only (R1's held
  labels, Early VIC inside endocardium, stay held), after Feng's direction, each cell also takes a per-metacell, category-mean-free change fitted by optimal transport within its
  category from E9.5 to E10.5 (damped by 0.5, zeroed below log 1.05, clipped at log 3, expected depth kept; shares, cell
  numbers and scope untouched). Evidence: on B0's base, T11 passed the adoption rule on two cuts of one E10.5
  library (+1.13 and +1.23, ahead in 10 of 10 seeds on each), and its controls T20 to T25 passed
  (`experiments/harness/TRIES.md`; `experiments/harness/results/t11_g_mom/`,
  `experiments/t1_scope/results/t20_t24_controls/`, `t22_composition/` and `t25_mixfree/`). All of it is one E8.5 -> E9.5
  fit read on cuts T11 was selected on: the ten seeds are draw and scoring noise, n = 1 for the biology, and none of it
  is adoption evidence for R3b. A size-matched control with nothing developing keeps about 40% of the gain (T23: +0.44
  and +0.46 over B0), and with the mix matched and each type's pseudobulk held the gain is +0.32, 6 of 10 seeds at the
  +0.3 bar (T22, reported, deciding nothing). On M1's base the factors kept 0.96 to 0.99 of their raw gain (T26, a
  labelled in-sample diagnostic, counted, void by its reproduction gate; the captain accepted its numbers as no
  interference).
  At the test pair the fit crosses the change from an E9.5 embryo portion to a dissected E10.5 heart, which no
  backtest's fit does; the per-gene sizes carry (Spearman 0.79 and 0.68) but the leading axis does not (absolute cosine
  0.20 and 0.29 against ceilings of 0.91 and 0.84; N21), and the one-day step is applied unscaled to the two-day
  horizon, which is untested.
- **What the T1 scorer rewards.** On the released E9.5 -> E10.5 transition the changed genes are mostly highly
  expressed genes going down, which de_score's expression-magnitude null names 39-41% of; a copy names 11%
  (`experiments/harness/results/t1_de_v1/`). Better direction raises de_score but has not cleared its floor on this transition: in the labelled diagnostic (B0 with Feng's direction at 1/5) de_score skill is 0.480 against B0's 0.450, with raw overlap 0.347 against the null's 0.388 (`experiments/t1_m1/results/diag_v1/`).
- **Emergent populations:** none in the entry yet. The organisers have not answered on qualitative in-window
  background; births by the legal route (states present in legal post-E14.5 chamber data and absent at E10.5, their
  states and shares taken only from E10.5 and those data by a rule written before any run; the cardiomyocyte and blood
  shares stay as the s1 and s2 rules set them; `docs/DATA_POLICY.md`, 2026-11-04 (c)) were designed and pruned on
  2026-11-04 (`experiments/harness/TRIES.md`).
- Training data: released T1 E8.5, E9.5, E10.5; Feng 2022 E9.5 and E14.5 CD1 hearts (X002), reduced from GEO by
  committed code that drops E10.5-E13.5 cells before anything is written. The E14.5 cardiomyocyte shares of Feng's two
  strains (label metadata only, in-window rows dropped at read time; X002) bound where the s2's 60% sits. JAX/Qiu 2024
  for diagnostics only (X003; embryos 23-53 excluded, as in that row). R2's held genes are module D, the 32 symbols of
  Supplementary Table 1 of van den Brink et al. 2017 (X007), and module A1, the ^Hb[ab]- haemoglobin genes among Feng's
  symbols (X002), chosen on the statement of Young and Behjati 2020 (X009); no cell, count or expression enters from
  them, and X008 supplied only the dropped nuclear module and part of R2's premise. A whole-heart E14.5 cardiomyocyte
  share from a sample coded Theiler stage 22 (X006) informed the s2's rationale in October and was withdrawn on
  2026-11-05 under our staging rule; no rule changed, and no cell, count or expression from it enters any model
  (`docs/DATA_POLICY.md` 2026-11-05 (c), (e)). The open-source libraries the method runs on are listed in X005 (POT also
  runs R3b's within-category transport).

## Task 2: spatial-temporal (heart E12.5, embryo E7.75)

- **Heart E12.5 — B0 with the pre-registered size:** real E10.5 tissue (synthetic records dropped), sampled at random
  to the pre-registered count (the board's maximum: the release-day density rerun, by a rule fixed in advance, kept it on both test
  boards; `vec.t2_build.n_cells_from_decision`, decision dfee546), moved in count space by half the E9.5 -> E10.5 change of
  each harmonised lineage; coordinates grown by the pre-registered size procedure (core = cardiomyocytes +
  endocardium; growth from E9.5 -> E10.5 only; a fixed three-rule average justified without an assumed rate; z held at
  E10.5 when the sectioning plane matches); the two submissions are the two scope hypotheses (H1; H2 as the half trim, amended after its released-data outcome,
  `docs/S2_RULES.md`). E10.5's
  composition is held: E9.5 and earlier were embryo portions and E10.5 a dissected heart, so a composition trend through
  them would carry that change of sampling on to E12.5 (dropped by the captain, 2026-11-01). On the pseudo-tasks B0
  beats a copy by the adoption rule (+0.47, 9 of 10 seeds; +0.90, 10 of 10) and is ahead of the floor (+0.43 and
  +1.45), short of the rule on the one-day task (7 of 10 seeds); the size cannot be backtested.
- **Embryo E7.75 — displacement interpolation (FGW):** a balanced fused Gromov-Wasserstein coupling of 2,000 cells per
  bracket (E7.5, E8.0) on 30 expression principal components and within-cloud distances (structure weight 0.5), the
  proper rotation that maps one bracket onto the other under the coupling, every cell moved by its nearest anchors'
  displacements t (or 1 - t) of the way, at the nominal time fraction t = 1/2; the population is a mixture of
  unedited real cells of both brackets (half each), scaled to the pre-registered embryo size (PCHIP of log RMS radius
  through the released embryo stages). On the two pseudo-tasks that chose it as the entry (selected on, not tuned on; its settings were fixed before the run): 61.12 on E7.5 against 56.97 for the organisers'
  pseudobulk_shift and 58.83 for it with the size rule, and 58.78 on E7.25 against 53.86 for B0, every gap ahead in
  10 of 10 seeds (`experiments/harness/results/emb_fgw_v1/`). It loses on shape on E7.5, where E8.0's sampling differs;
  the s2 is fgw_shape0, fgw's mixture in E7.5's geometry, chosen by a rule fixed before its run (`docs/S2_RULES.md`).
  Since the captain's ruling of 2026-11-13 the s2 is `fgw_aligned_series` instead (above; row R4).
- **Spatial hygiene:** random rows only, jitter of at most 0.25 nearest-neighbour distances to break duplicate
  coordinates, proper rotations only, native coordinates (`experiments/density/REPORT.md`, `experiments/coords/REPORT.md`).

## Task 3: perturbation (beta-catenin KO, E8.75) — pre-registered

- **Per-gene factor first.** On the released Mab21l2 knockout, most of the measured difference from its wild type is one
  per-gene factor shared by every cell type: estimated from that file's non-Mesp1 compartments and applied to all
  wild-type cells it scored 80.5 against 50.8 for wild-type resampling, in-sample on that file, while rules built from
  the broad phenotype alone scored about 51 (`experiments/t3_decomp/results/mab21l2/`). About half of the factor is
  specific to the WT reference (`docs/PLAN.md` 4.4), and the cross-knockout checks that followed did not support it
  (A10, above).
- **Composition second.** Second-heart-field, outflow-tract and right-ventricle-lineage states depleted at a magnitude
  calibrated on the released Gata4 replicates, counted only where it exceeds the replicate difference.
- **Decided by pre-registered rules** (`experiments/t3_decomp/rules.yaml`) run once on the released Gata4 files; the two
  submissions differ in how the factor is sourced. Because the factor is measured in compartments the knockout cannot
  touch, it may partly reflect sample processing: unless the Mab21l2 factor transfers to Gata4, the two submissions do
  not both rely on the Gata4 factor, and if no check supports it, the rules' s1 is unchanged wild-type cells
  (with the composition term) and only the rules' s2 carries the factor, at half strength (a fallback that did not arise:
  D2 passed, and the entries' s1 and s2 are those under Outcome, below). The released Gata4 files enter the entry's s1 on
  the team's reading of Rules s.10 and the predicted organiser answer, not on a ruling: the organisers were asked on 19
  and 28 October and reminded on 21 October, with no answer on record (Data provenance, above), and `docs/S2_RULES.md`
  says what switches if they restrict the files before a T3 upload. The allowed, refused and barred branches' files are
  all built. The rules were
  drafted with an AI assistant and reviewed and approved by the team before any Gata4 file was downloaded; while drafting
  them the author saw the T3:gata4 validation-board anchors quoted in `docs/CHALLENGE.md`, and no threshold, grid or
  magnitude was derived from them (`experiments/t3_decomp/rules.yaml`, notes).
- **Outcome (2026-10-28, run once):** the Gata4 replicates transfer to each other, the Mab21l2 factor does not transfer
  to Gata4. s1 = the allowed branch's first submission; s2 = the refused branch's first (wild-type cells, the central
  composition at the Mab21l2 magnitude, no number from the Gata4 files), uploaded first (approved 2026-11-01); since
  the captain's ruling of 2026-11-13 it also carries the wild-type-only factor per lineage group (above; row R5). The
  harness gives these entries no verdict (every fold is one they were selected or fitted on); they stand on the
  pre-registered rules and the captain's ruling of 2026-11-01, although the two assumptions behind those rules' bets
  have since been killed: A10 (D11b failed on 2026-10-28, K2a found no support and a review traced most of the factor's
  carry to a stage offset; killed 2026-11-05) and A11 (its proxy failed; killed 2026-11-05).

## Shared architecture (Generality Award)

**Request.** Team Mithril asks to be considered for the Generality Award (Rules s.14), for the tasks it enters. What is
common, what differs and where the boundary falls are set out below, including where the tasks share a design but no
code. Whether this structure counts as one shared architecture is for the organisers to judge.

- **Common core** (`src/vec/core.py`, feat-020). Every prediction is built from real released cells by one
  "transport-then-adjust" recipe, each head running the steps it needs:
  - allocate a target composition to whole cell counts;
  - draw real cells from the source pools;
  - edit their expression in count space (per-type, per-metacell or per-gene factors by binomial thinning or Poisson top-up, a
    depth change);
  - write log1p CP10k.

  - couple two populations with positions (fused Gromov-Wasserstein, a proper rotation under the coupling, barycentric
    partners, nearest anchors; added 2026-11-01 for the T2 embryo). The head turns these into the cells' moves.

  Births from a legal source and leave-one-transition-out calibration are not in the entries.

  Which steps run on which board. Task 1: allocation (`quotas`), draw (`draw_least_used`), count-space edit
  (`apply_factors`), log1p CP10k (`log1p_cp10k`). T2 heart: draw (`draw_balanced`), count-space edit (`apply_factors`,
  through the Task 1 head's `momentum_shift`), log1p CP10k. T2 embryo: draw (`draw_balanced`), coupling (`fgw_coupling`,
  `coupling_rotation`, `barycentric`, `nearest_anchors`) and the head's moves; the cells are unedited, so no edit or
  write step. Task 3: the pre-registered script's own `alloc`, `draw_rows`, `apply_gene_factor`, `rescale_depth` and
  `lognorm` (the refused-branch file, s2 as planned, edits no expression). No `vec.core` function runs on more than two
  of the four boards: what the boards share is the design (every submitted cell is a real released cell, edited only in
  count space or not at all), and each step is a separate function in `vec.core` or in the Task 3 script.
- **The core's steps and the four predictions** (the captain's guiding principle, 2026-11-04). Drawing real cells answers
  which states appear and in what proportions: to a target composition on Tasks 1 and 3, at the source pool's own
  composition on the heart, and as a two-bracket mixture on the embryo; births would add states no released stage
  has (not in the entries yet; on the embryo the later bracket's real cells bring them); count-space edits carry the
  expression change on Task 1 and the heart, and on Task 3 the genotype's per-gene response (the embryo's entries and
  Task 3's refused-branch file edit no counts); coupling and moves carry the spatial organisation on the embryo only,
  and the heart's comes from E10.5's tissue and the pre-registered size. The published baselines correspond to every
  step left at identity (our re-implementations, `vec.baselines`, are separate code outside the core), so an entry's
  change beyond its baseline is the list of steps it switches on, and the heads decide which.
- **Task heads** (what differs):
  - T1 (`vec.t1_scope`, `vec.t1_state_transport`, `experiments/t1_m1/feng_dirs.py`, with `vec.harness` computing the
    target shares and running the within-type fit, and `experiments/builds/build_test.py` choosing the cells that fit
    edits): the 32,285-gene output panel, chamber categories and target shares, and the expression factors (Feng's
    external direction with R2's held genes, `experiments/t1_hygiene/`; momentum in the backtests; and the within-type
    factors, fitted by entropic optimal transport inside each category and applied through the core).
  - T2 (`vec.t2_build`, with `vec.harness._b0_t2` running the heart's draw and per-lineage shift through
    `vec.t1_scope.momentum_shift`): each board's gene panel (about 500 genes) and the coordinate head (cell count, the
    random-row rule, jitter, the pre-registered size procedure, proper rotations, the embryo's displacement moves and
    bracket shares).
  - T3 (`experiments/t3_decomp`): the genotype head (knockout composition from harmonised labels, and a per-gene
    factor estimated in compartments the knockout cannot touch), the wild-type pool and its smoothed proportions, the depth
    ratios, each cell's position (its source cell's, with the script's own duplicate jitter) and the script's own log1p
    CP10k step.
- **Boundary:** heads supply the source pools, the target composition, the factors and the coordinates. On T1 and T2 the
  core's functions draw every submitted cell, make every count edit and, on the embryo, compute the coupling, rotation,
  barycentres and nearest anchors; each head orders those steps itself, and the T2 head makes the embryo's moves, its
  bracket counts and all coordinate changes. The T3 script was written first, then pre-registered and hash-pinned, so it
  builds its cells with its own functions and calls none of the core's. `tests/test_core.py` pins five of them to copies
  held in the core (`alloc`, `draw_rows`, `thin`, `rescale_depth`, `apply_gene_factor`: same body and signature, same
  outputs and generator states); the script's log1p CP10k step, `lognorm`, is not a copy and agrees with the core's to
  one float32 ulp on the test's synthetic counts. No build calls those five core functions. Golden tests on synthetic
  inputs and a comparison on released stages show the move into the core changed no T1 or T2 builder's output
  (`experiments/core/REPORT.md`); they predate and do not cover the steps added since (Feng's direction, the
  within-type factors, the embryo's coupling and moves).
- **Decision point D3 (read 2026-11-10, ahead of its date of 2026-11-12):** claim the award unless the best specialised
  variants beat the shared-core builds by more than 3 points of overall standing (T1 + mean of the T2 boards + T3) on
  honest validation, paired over the emulator seeds; in that case use specialised code where it wins and drop the claim.
  That comparison was not run. No specialised pipeline has a score on any board: the one coded, try E14 (a
  negative-binomial VAE for the embryo's expression), stopped at its target-free fidelity gate before any counted run
  (`experiments/harness/TRIES.md`), and the T1 latent-transport pipeline the plan named (M2) met its build condition on
  2026-10-28 and was not built (`experiments/t1_direction/results/E10.5/decision.json`; `feature_list.json` feat-015).
  The shared-core side can be scored as built only on the embryo: Task 1 and the heart are read through stand-ins, and no
  backtest runs the Task 3 entries (`experiments/harness/TRIES.md`, "Where the entries stand against the published
  baselines"). The claim therefore stands by the rule's default and not by a comparison passed; nothing here shows what
  specialised pipelines would have scored, in either direction. On T1 and both T2 boards, both submissions call the core
  to draw the submitted cells and for every count edit they make (the embryo entries make none), and the embryo for its
  coupling, rotation, barycentres and nearest anchors. The Task 3 files call no core function: the pre-registered script
  runs its own, five of them pinned identical to copies held in the core (Boundary, above).
