# Questions for the organisers

Send privately (email virtual.embryo.moonshot@gmail.com, or a Slack DM to Weize Xu), not in the public channel: some
questions reveal method choices. Rules s.10: "an answer given in advance is not a finding of misconduct afterwards".
Record each answer with its date in `docs/DATA_POLICY.md`.

Status: **sent** by email to virtual.embryo.moonshot@gmail.com from be23b041@smail.iitm.ac.in on 2026-10-18. Answers: _pending_ (record them in `docs/DATA_POLICY.md`). The leak report went separately by Slack DM the same day.

Dropped because they are already answered (see `docs/DATA_POLICY.md`):
- T3 literature: gene-level knockout findings are prohibited, broad phenotype background is allowed (2026-10-15).
- T2 heart extrapolation window equals T1's (2026-10-15).
- External data is allowed for T3 within the windows (GSE208162 ruling, 2026-08-18).
- Scoring conventions: no renormalisation; cell limits inclusive and enforced (2026-10-16).

Deferred to the 2026-10-20 release (ask only if `panels/index.json` leaves them unclear): each test board's reference
stage (`ref_cells`) and cell bounds.

Not asked because our plan already takes the strict reading: somite-to-day boundary cases (we drop ambiguous embryos),
Tier-B mutant datasets (unused), cross-species foundation models (unused), selection of the second submission among
pre-registered candidates (Rules s.10 already allows it).

---

**Subject:** Human Team: four rules questions before the final phase (cross-task data, developmental background, board metadata, Generality scope)

Hello organisers,

We are a Human Team (captain: Vijayavallabh J). Before the final phase we would like to confirm four points, so that
we follow the rules exactly.

1. **Released competition data across tasks in the final phase.** Once the validation answers are released, may a
   released file from one task be used as training input for another? Specifically: (a) Task 1 scRNA (E8.5, E9.5,
   E10.5) for the Task 2 heart E12.5 board and for Task 3, for example to learn wild-type gene programmes and to calibrate
   scRNA against MERFISH; (b) Task 2/3 MERFISH (heart E8.25 to E10.5, wild type E8.75 and E9.5) for the Task 1 E12.5
   board. Our reading is that released competition data is training data for every task in the final phase, and that
   the "after E9.5 up to and including E13.5" exclusion applies to external data only. Is that correct?

2. **General developmental-biology background inside a protected window.** For Task 1 and the Task 2 heart board we
   would like to use textbook-level qualitative knowledge of heart development between E10.5 and E12.5 (for example
   that the epicardium covers the heart, the coronary plexus begins to form, valve primordia appear and the outflow tract
   septates) to decide which cell populations should appear in our predictions. We would not use any measured values
   from that window (no proportions, sizes or expression data). By analogy with your 2026-10-15 answer about broad
   phenotype background for Task 3, we read this as allowed. Is it? (This overlaps with a question from another team on
   16 October in the public channel.)

3. **Published board metadata.** May we use the `truth_cells` value that `panels/index.json` publishes for a board to
   choose how many cells we submit? Neighbourhood statistics depend on spatial density, so the count matters for the
   Task 2 boards. We will not use the anchor values (floor and ceiling metric values) for anything other than reading
   scores.

4. **Generality Award scope.** If a team enters all three tasks, must the one shared architecture underlie all three, or
   is it enough that it answers at least two of them? And does a shared population-transport model count (unbalanced
   optimal-transport coupling, a conditional velocity/growth network, decoding by retrieval of real cells and one
   calibration protocol, with task-specific encoders, a coordinate head for Task 2 and a genotype head for Task 3)?

Also: (5) will final evaluation use veckit `main` at commit cc3859a (severity_slope, scale_log_ratio), or the metric set
of veckit PR #6? And (6) which time zone applies to the 2 December deadline?

Thank you for running the challenge.

Vijayavallabh J

---

## Follow-up (ready to send, 2026-10-19)

User's decision (2026-10-19): send the full follow-up. Not sent yet; the user sends it by email from their address.
Duplicate counts come from `experiments/qc/` (update them if the QC run changes them). Points whose answer would help
every team ask for a public answer, so we never hold private information about held-out data.

**Subject:** Human Team: follow-up on the scorer build, Task 3 data use and two data issues

Hello organisers,

Thank you for your earlier answers. A few follow-up points from checking the released files and the scorer. Where an
answer would help every team (points 1, 2, 5 and 6), please feel free to answer publicly instead of to us.

1. **Scorer build and canonicalisation.** Which veckit build will the test boards and the final evaluation use? The
   public validation scores look consistent with the build before the right-handed-frame fix that was merged on
   veckit `main` on 2 October. Under that build a mirror copy scores the same as the original, and `occupancy_dice`
   depends on the sign pattern of the SVD used for canonicalisation, which changes with row order. So two otherwise
   identical submissions can get different Dice values, and a team could deliberately submit opposite "parities". We
   have not done this and will not unless you tell us it is acceptable; we are reporting it so you can decide. If the
   build changes, will the published shape anchors be regenerated?
2. **Severity with the wrong sign.** When the predicted severity slope has the wrong sign (or R² ≤ 0.01), the current
   scorer returns log(1e-3), which equals the floor (skill 0.5). The Evaluation page and the proposal call this case
   "undefined". Which applies in the final evaluation: skill 0.5 or skill 0?
3. **Released Gata4 data in Task 3.** Once the Gata4 validation answers are released, may we use those released
   knockout files as training input for the beta-catenin test prediction? Rules s.10 still lists the Gata4 knockout at
   E8.75 among the held-out genotypes, while the Data page says validation answers become training material in the
   final phase. (We will keep excluding external beta-catenin and Gata4 perturbation data either way.)
4. **Duplicated and interpolated records.** Some files contain the same cell more than once: identical integer count
   rows in adjacent z sections, either at the same xy position (443 in E8.0, 1,164 in E8.25_late, 317 in E8.75) or
   shifted in xy (about 2,460, 560 and 520 more). In the Mab21l2 knockout file, sample_id 35 appears to be stored twice
   and sample_ids 18 and 19 look like one section placed at two z levels (2,811 records together). Separately, E8.0,
   E8.25_late and E8.75 contain 200, 2,558 and 614 records with half-integer counts, each equal to the average of two
   cells in the neighbouring sections. We found none of these in E6.75, E7.25 or E9.5. Are they artefacts of
   registration or section interpolation, and do the validation and test truth sets keep or remove them?
5. **z coordinate.** Is `spatial_3D` z the section position times the section thickness, in the same unit as xy, for
   every file? Were sections registered rigidly, or with expression-aware or non-rigid alignment (some knockout
   sections are not rigid copies of their slide coordinates)? If you can say so publicly: do the test stacks follow
   the same acquisition protocol (section thickness, similar number of sections) as the training files?
6. **truth_cells and tissue size.** Further to our earlier question about using `truth_cells` to choose the submitted
   cell count: may it also inform the predicted tissue size (for example a cell-count-to-volume relation fitted on the
   released files)? We will not use it for that unless you confirm.
7. **Stage conventions for external data.** Could you share your convention for converting somite counts, Theiler
   stages and morphometric stages to E-days, or the somite counts of your own E8.25_late, E8.75 and E9.5 embryos? Two
   concrete cases, which we currently exclude: atlas embryos staged E9.5 by somite count but harvested from a litter
   labelled "E9.75", and embryos staged E13.75-E14.0 morphometrically but harvested from gestational E12.5 or E13.5
   litters. Are they inside the excluded window?
8. **Leak report, addendum.** Further to our private Slack note of 18 October: the same held-out values also appear
   in the veckit 0.1.2 wheel on PyPI and in the description of veckit PR #6. We have not used them.

Thank you,

Vijayavallabh J
