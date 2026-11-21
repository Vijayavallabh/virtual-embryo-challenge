# The embryo's copy-rows file: the specification in force (practice row E50)

Committed with row E50's placement (`docs/PLAN.md` decision log 2026-11-21 (l) and (m)). The texts below are the
analysts' own words, verbatim, from the rounds of 2026-11-21: a fresh analyst specified the construction from the
public scorer and builder code, a cold challenger amended it (the final version), and after a code review a further
fresh analyst settled the open points. The AI assistant (the main session) is recused on the embryo's occupancy term:
it placed these texts and wrote none of their content. Where an earlier text and a later decision differ, the later
decision holds; section 9 lists the differences. No number here comes from a test score.

## 1. The construction (the cold challenger's final version)

The copy-rows file (the investigator's Option A, 'the second construction'). Option `--copy-rows` (with
`--copy-rows-test-ok` on a test board), added to the handed-over form's command (`build_test.py embryo --method
fgw_shape0 --positions floor_exact --expression <that file's own>`).

1. ROWS, ORDER, COUNT, COORDINATES. These are exactly the floor copy `vec.evaluate.subsample(Stage(E7.5 as released,
read_stage, no dedup), 0.1, 0)`, as `build_test.floor_copy_coords('T2/E7.5.h5ad')` returns it. The coordinates are
float32, bit for bit, in the copy's sorted released order. The count is n_c = max(2, int(round(0.1 * N_E7.5
released))). There is no size rule, centring, translation, jitter, permutation or row-order step.

2. EXPRESSION. n_b = int(round(0.5 n_c)) slots keep the expression of their own E7.5 cell. n_a = n_c - n_b slots take
E8.0 cells:
- The E8.0 rows are `draw_balanced(len(E8.0 deduped), n_a, default_rng([seed, 1]))`.
- The E8.0 slots are every copy slot whose released row the duplicate filter drops
  (`load_stage(drop_duplicates=True).uns['qc']['dropped_rows']`). The rest is a uniform subset of the other slots,
  drawn by `default_rng([seed, 51])` (`COPY_ROWS_STREAM = 51`; streams 47, 48 and 50 are taken).
- If the dropped slots exceed n_a, the build refuses with no value shown.
- E8.0 rows go to slots by `linear_sum_assignment` on squared distance. The distance is between the E8.0 row and the
  slot cell's own released expression, both projected on `_joint_pc_loadings(Eb[ib], Ea[ia], 30)`, with anchors from
  `default_rng([seed, 4])`.
- Mixture: E7.5 slots carry Eb[the deduped index of the slot's row] and E8.0 slots carry Ea[row]. Series: the (stage,
  rows) go through `series_rendered_expression(before, after, res, genes, 0.5, [seed, 7])` unchanged.
- The --expression is the handed-over file's own. The main session reads it from that file's committed build record.

3. WRITE PATH. `build_submission(..., allow_duplicate_coords=True, allow_below_truth=True)`, in this path only.
- Before it, build_test checks [min_cells, max_cells] itself and exits with a message that holds no number.
- Any ValueError from build_submission in this path is replaced by a message that holds no value. The stock message
  prints n, min_cells and max_cells.

4. SELF-CHECK before writing, and again on the written file read back: the stored obsm['spatial_3D'] equals
floor_copy_coords exactly (shape, order, float32 bits). Otherwise the build exits.

5. HARNESS ARM `fgw_copyrows` (and `_series` if the handed-over file is the series form). The arm uses the very
`ref_cells` AnnData that run_task builds for the floor, after any within cut, and copy seed s = the scoring seed
(build seed minus BUILD_SEED_OFFSET). Only the arm that matches the handed-over file's expression counts.

## 2. The count and the size (the cold challenger)

The file takes the floor copy's count, a tenth of E7.5 as released, computed from released data alone. It departs from
the count rule (docs/PLAN.md 3.6, approved rule 7, the decision file's max_cells), may fall below the hidden truth
count (allow_below_truth is used, and says so), and drops the pre-registered embryo size (R9). The departure is
forced. Any other count, row, order or scale changes the rows that _match_n takes, or the stored bits, and the
equality with the floor is lost. Section 3.6's reason ('occupancy_dice thins the truth to the submitted count') does
not bite, because the floor is thinned the same way. The count can still cost on the expression and neighbourhood
terms, and practice row E50 prices that. No test-board count enters anywhere. `_t2_count` still runs for provenance,
but does not set the count.

## 3. What equals what under each scorer build (the first analyst)

IF ASSUMPTIONS (a) TO (c) HOLD

Exactly the floor's, under 'main' and under 'prefix' alike:
- The file's coordinate array is the floor's array (same count, rows, order and float32 bits). `_match_n` therefore
  hands both the same rows for every truth count and every metric seed.
- So d2_shape, occupancy_dice and scale_log_ratio, the three members of the shape_growth group, equal the floor's
  values exactly and sit at the floor's skill. The unranked sliced_wasserstein equals the floor's too.
- The raw shape_growth group skill is the floor's. With the emulator's assumed 'floor_normalised' gate, the gated
  group is the floor's group times min(1, the file's own neighbourhood skill / 0.5). It is exactly the floor's level
  when the file's neighbourhood skill is at least 0.5. The real gate is not published.

The file's own:
- de_score, de_direction and mmd_u read expression only. They depend on the n_c expression rows (count and slot
  choice).
- variogram and neighborhood_mmd read expression on positions. They depend on which expression sits on which copy
  position.

What the analysts' findings can still do:
- The handedness coin (prefix): the file's SVD parity on the matched rows is the floor's own parity, so the file
  always falls in the floor's mode. It shares the floor's luck whatever the hidden half's parity, and has no coin of
  its own.
- The turn of the second and third axes: the file's axes are the floor's axes, so there is no turn relative to the
  copy. Any turn between the copy and the hidden truth is the floor's own and is shared.
- Draw-to-draw frame noise (reading S1): the file has the floor's exact draw, so it carries the floor's noise and none
  of its own.
- Under 'main' there is no coin, and the same equalities hold.

IF THE SERVER'S COPY IS NOT DRAWN THE RELEASED WAY (another seed, another fraction, a re-ordered or cleaned E7.5, or
the submission subsampled or sorted):
- By the code, the file is then one more uniform 10% copy of E7.5, independent of the server's floor.
- Its position terms are distributed as a plain copy's against the truth. Each is centred near the floor's level, with
  the copy-to-copy spread the S1 check measured on released E7.5 (overlap with a copy at about the copy-against-copy
  level, sd 0.029). Under prefix its parity coin is independent of the floor's.
- On this term it is then in the same position as the handed-over exact-position file, which entry (l)'s check found
  distributed like a copy. It keeps none of the pin and still pays the count and size departures.
- If only the metric seed or LAPACK differs, the equality still holds.

## 4. The assumption about the server, and what if it fails (the first analyst)

Assumption:
- The live test board's floor anchor is subsample(E7.5 as released, 0.1, seed 0), the bake recovered on the validation
  boards.
- The submission is scored as held: full count, file order, no subsample, sort or dedup.
- The floor and the submission share one metric seed in the same metric code.
No metric-seed value and no LAPACK sign behaviour is assumed: equal arrays give equal values.

If it fails:
- The file is a fresh uniform 10% copy of E7.5, independent of the floor. Its three shape terms are distributed like a
  copy's around the floor's level, and under prefix its coin is independent of the floor's.
- On the shape group it is then on average like the handed-over file (also distributed like a copy) and keeps no pin.
  It still pays the count and size departures on the other terms.
- Nobody outside the organisers can confirm the assumption.
- The returned scores must not be used to test it: comparing this file's per-term values with another upload's to
  learn the server's draw, the floor anchor or the hidden half's mode is the route decision log 2026-11-18 (k) bars.

## 5. How likely the file lands at the baseline's level on these terms (the cold challenger)

About 0.7, given that the build passes the board's limits. I did not estimate that pass from the board, and must not.
The figure is roughly the product of three parts:
- About 0.85 that the test board's floor is subsample(E7.5 as released, 0.1, 0) in released order. The bake was
  recovered exactly on all three validation boards, but the test bake is not shown.
- About 0.95 that uploads are scored as held, with no subsample, sort or dedup, as in the pinned score_h5ad.
- About 0.9 that the stored floor anchor and the live scoring of uploads use the same build, metric seed and numeric
  library. Under main this matters little; under prefix, a LAPACK sign difference could split equal arrays.
Two further limits. Under the assumed floor-normalised gate, the shape group's board points can still fall below the
floor's if the file's own neighbourhood skill is under 0.5. And if the assumptions fail, the file is a plain 10% copy
on these terms.

## 6. The fix round's decisions (the further fresh analyst): the split

(b). Keep the code as built. The practice arm splits by the task's own time fraction t, n_b = int(round((1 - t) n)).
The row text and the read script's wording change to match. t2_build.COPY_ROWS_RULE already says round((1 - t) n), so
only the row text, the read script's docstring and the harness candidate's "what" text change.

Reason: the practice read exists to check the equality (reading (i)) and to price what the copy-rows file gives up
against the file it follows, the handed-over exact-position file (fgw_shape0x). At the test pair the two files differ
only in the forced departures: the floor copy's rows, order and count, the slot choice, and no size rule. Both use a
0.5 split there. On the practice tasks fgw_shape0x splits at the task's own t: _fgw passes frac = (target_day -
lo)/(hi - lo) to interpolate_clouds, which I read on the branch. Under (a), on the eight tasks where t is not 0.5, the
paired difference in (iii) and (iv) would carry a second change, a different E7.5/E8.0 share. That change does not
exist at the test, so it would make the price less faithful to the test comparison. Under (b) the arm differs from its
paired control only in what also differs at the test.

(b) also matches the specification's own wording: spec, "Kept: the 1:1 bracket share ... as interpolate_clouds splits
them", and fin test 3, "the n_b/n_a split equals interpolate_clouds'". The equality in reading (i) does not depend on
the split. The test builder keeps exactly 0.5 (build_test._build_copy_rows passes 0.5 to copy_rows_cloud and to the
series step), so the test file is the same under (a) and (b). Under (b), (iii) and (iv) still decide nothing.

## 6a. The cloud guard

DELETE it from the builder path, and delete it from the practice arm (harness._copy_rows) too, so the two are made the
same way, which was the engineer's aim.

Reason: the copy-rows coordinates are by construction a uniform 10% draw of the released E7.5 cloud. The guard
(t2_build.cloud_guard) compares a cloud with random same-count draws of that source, with a 3% RMS tolerance and block
and outside flags. On this path a flag could only be a chance event on a correct file, or on the floor's own array. It
can never catch a wrong file. The bit-for-bit self-check, before writing and on the file read back, is strictly
stronger. Reading (i) in the read script already checks the arm's coordinates against the floor copy bit for bit per
task-seed. Keeping the guard would only add a way for a correct file not to exist.

What the practice arm does: _copy_rows builds the copy, renders the expression, clips, then calls
build_submission(..., allow_duplicate_coords=True, allow_below_truth=True). There is no guard. Its coordinates are
checked by reading (i). Refusals are only those listed in (ii). Controls fgw_shape0x and fgw_shape0x_series keep their
own guard, which is unchanged.

Exact sentences that change:
- fin final_construction item 5: delete nothing there, since it names no guard. In fin final_builder_and_tests_spec,
  harness section, replace "Finish path: cloud_guard against the reference cells' coordinates, then
  build_submission(..., allow_duplicate_coords=True, allow_below_truth=True)." with "Finish path:
  build_submission(..., allow_duplicate_coords=True, allow_below_truth=True); no cloud guard on this path (the
  coordinates are the floor copy's by construction and reading (i) checks them bit for bit)."
- In the same section, "A refusal is recorded per task-seed" stays.
- final_buildability: no new condition. Its conditions stay 1 to 4 (amended statement below).
- Row text (ii): the refusal list does not gain a guard entry. The reviewer's lesser point 4 is moot.
- harness._copy_rows docstring: "Finish: the cloud guard against the reference cells' coordinates, then
  build_submission with duplicate positions and a count below the truth's allowed (this path only)." becomes "Finish:
  build_submission with duplicate positions and a count below the truth's allowed (this path only); no cloud guard,
  because the coordinates are the floor copy's by construction and reading (i) of row E50 checks them bit for bit."
- Engineer departure 4 is withdrawn.

## 6b. The reviewer's points on readings (i) and (ii)

Point 3: ACCEPT as worded. Reading (i) gains "and its coordinates equal the floor copy's bit for bit". The code
already enforces this (departure 8), and it is the construction's definition, so the row should say it.

Point 4: REJECT as moot. The cloud guard is removed from the arm and from the builder path (decision 2), so a guard
flag cannot be a refusal. The refusal list in (ii) gains instead "a copy of fewer than three rows". That follows from
accepting lesser point 2: REFUSALS then holds three categories (dropped_over_share, count_outside_limits,
copy_too_small), and the row lists the same three. Any other error from the arm raises, which makes the run void as a
missing arm, and the row's void list says so.

## 6c. Which arm counts, and the expression

Verified in the files.

experiments/builds/results/embryo2_2dccfee.txt records the set embryo2 at 2dccfee, built twice on the DGX. Its
command: "build_test.py embryo --method fgw_shape0 --positions floor_exact --expression mixture --out
<O>/builds/T2embryo_s2x.h5ad". Its sha256 is a03611a343faf1d37edf06ce0e2bccc86cce6f14d3580b823864e5faf3ba03f9. In
experiments/builds/build_all.sh, EMBRYO_EXPRESSION=mixture (line 62) and the set embryo2 passes --expression
"$EMBRYO_EXPRESSION" to T2embryo_s2x.h5ad (lines 115-119). The --expression series line (130) writes
T2embryo_s2xs.h5ad. In results/manual_uploads.csv the embryo row of 2026-11-20 is T2embryo_s2xs.h5ad, sha256
bc50d012..., outcome scored. No row names T2embryo_s2x.h5ad.

So the handed-over, unscored file is the mixture form (fgw_shape0x), and the series twin is the one already scored.

- The arm that counts: fgw_copyrows, the mixture arm. fgw_copyrows_series is not run and nothing of it counts. The
  read script runs only ARM, so this follows from the constant.
- HANDED_OVER_EXPRESSION = "mixture". ARM resolves to fgw_copyrows.
- Test-board command: --expression mixture.

What else changes because of it:
- The read script's docstring (lines 7-9) and the comment on line 77 claim series. Both are rewritten.
- The engineer's builder_usage text says "series for T2embryo_s2x". It is wrong and is replaced by mixture.
- The primary control PRIMARY = fgw_shape0x is now the handed-over form itself. fgw_shape0x_series stays only as the
  secondary table and in the gate. No code change is needed there.
- Row text: the candidate cell names fgw_copyrows alone, and the baseline cell names the handed-over file as
  fgw_shape0x, --expression mixture, T2embryo_s2x.h5ad.
- The decision-log draft and the rules-note amendment say "with the handed-over file's expression rule". Make it
  explicit: "(mixture, as the exact-position file: each cell carries its own measured panel expression)".
- The data-used text needs no change. It already says E8.0 is used "for the expression of the cells placed on half of
  those positions", and under mixture those cells carry their own measured expression with no series rendering.
- The method report line ("Its expression follows the exact-position file's rule") stays true.
- The fin note in attack (e), that the records disagree with the brief, is resolved: the records are right.
- The start condition "the default command reproduces file 1's bytes" uses the mixture command, with expected sha256
  a03611a3....

## 7. Buildability and the start conditions of the test build (as amended by the fix round)

The file can be built for the test board with no withheld quantity. It depends only on released E7.5 (positions,
order, count and the copy cells' own expression), released E8.0 (the expression of the cells on the other half of the
slots), frac 0.1, copy seed 0, the build seed, and --expression mixture, the handed-over file's own form. Equality
with the floor holds for every hidden truth count and every metric seed, because the file is the floor's array. A
bit-for-bit self-check enforces it before writing and on the file read back. There is no cloud guard on this path.

Conditions:
1. n_c must lie within the board's published limits. Code reads this, behind `--copy-rows-test-ok`. A refusal (or a
pass) shows one bit about the limits, so it is read by exit code only, and nothing about the limits or n_c is quoted.
2. `_hide_counts` takes its off-by-default branch: it checks len(X) == len(floor_copy_coords) and writes n_cells as
'floor_copy'.
3. The dropped slots must not exceed n_a, and the copy must hold at least three rows.
4. The stock build_submission limits error is never reached unwrapped on a test board, because it prints values.

Start conditions for the test build:
(a) Row E50 is placed, and its reading (i) holds.
(b) On the DGX campaign checkout, at the build commit, with the build_all default of 16 BLAS threads, the default
command `build_test.py embryo --method fgw_shape0 --positions floor_exact --expression mixture` reproduces sha256
a03611a343faf1d37edf06ce0e2bccc86cce6f14d3580b823864e5faf3ba03f9 (T2embryo_s2x.h5ad). Compare only with that DGX-built
hash.
(c) Two copy-rows builds from that one clean stamped commit give identical sha256.

If condition 1 or 3 refuses, or any start condition fails, there is no file and file 1 stands. The server assumptions
(floor drawn as subsample(E7.5 as released, 0.1, 0), uploads scored as held, anchors and uploads scored with the same
build and numeric library) decide whether the pin holds at the test. They do not decide whether the file can be built.

## 8. Data used (the cold challenger)

Released challenge data only. The E7.5 embryo stage as released, for the positions, order and count of a 10% copy
drawn the way the public scorer draws its baseline copy, and for its cells' expression. The E8.0 embryo stage, for the
expression of the cells placed on half of those positions. No outside data, no pretrained model, and no test score or
test-board count in any number.

## 9. Where an earlier text is superseded (listed by the main session; no content of its own)

- Section 1, item 2 and item 5, leave the expression to "the handed-over file's own": it is mixture (section 6c). The
  series arm is not run and a series file is not this file.
- Section 1, item 2, says n_b = int(round(0.5 n_c)): that is the test build. On the practice tasks the arm splits by
  the task's own time fraction (section 6).
- No cloud guard runs on this path, in the builder or in the practice arm (section 6a).
- The order of the embryo's files is the rule of decision log 2026-11-21 (l), which replaces this round's own order
  rule for two files.
- Sections 3 and 4 say that no LAPACK sign behaviour is assumed and that the equality still holds if only the metric
  seed or LAPACK differs. The cold challenger's version holds instead (sections 5 and 7): the stored floor anchor and
  the live scoring of an upload must also use the same scorer build, metric seed and numeric library; under the prefix
  build a LAPACK sign difference could split equal arrays.
- Section 3 cites "entry (l)'s check": that is decision log 2026-11-20 (l), the analysts' check of the occupancy term,
  not the order rule of 2026-11-21 (l).
