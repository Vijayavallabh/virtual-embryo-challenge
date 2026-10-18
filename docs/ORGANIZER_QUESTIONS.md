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

## Follow-up (draft, 2026-10-19, not sent)

Waiting for the user's decision. Points 1 and 2 report data problems that affect every team's scores, so sending them
early also protects us; point 3 asks whether a published field may inform a prediction, before we would use it.

**Subject:** Human Team: data questions on the MERFISH coordinates and duplicated records

Hello organisers,

Three follow-up points from checking the released Task 2/3 files:

1. **Duplicated records.** Some files contain the same cell twice: identical integer count rows at the same xy
   position in adjacent z sections. We count about 440 such records in E8.0, 1,164 in E8.25_late and 317 in E8.75. In
   the Mab21l2 knockout file, sample_id 35 appears to be stored twice (about 1,400 records at near-identical xy), and
   sample_ids 18 and 19 look like one section placed at two z levels (about 1,400 records). We found none in E6.75,
   E7.25 or E9.5. Are these artefacts? Will the validation and test truth sets keep or remove such records?
2. **z coordinate.** Is `spatial_3D` z the section position times the section thickness, in the same unit as xy, for
   every file? Were the sections registered rigidly, or with expression-aware or non-rigid alignment? Some knockout
   sections are not rigid copies of their slide coordinates. Will the test stacks follow the same acquisition protocol
   (section thickness, roughly the same number of sections) as the training files?
3. **truth_cells and tissue size.** Further to our earlier question about using `truth_cells` to choose the submitted
   cell count: may it also inform the predicted tissue size (for example, a cell-count-to-volume relation fitted on the
   released files)? We will not use it for that unless you confirm.

Thank you,

Vijayavallabh J
