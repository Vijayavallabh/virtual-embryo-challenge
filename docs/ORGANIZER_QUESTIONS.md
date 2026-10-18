# Questions for the organisers

Send by email (virtual.embryo.moonshot@gmail.com) or in Slack `#virtual-embryo-challenge` **before** relying on any
of these. Rules s.10: "an answer given in advance is not a finding of misconduct afterwards". Record each answer with
its date in `docs/DATA_POLICY.md`.

Status: **draft, not sent** (2026-10-18).

Already answered in Slack, so dropped from this draft (see `docs/DATA_POLICY.md`):
- T3 literature priors: gene-level knockout findings are prohibited, broad phenotype background is allowed (2026-10-15).
- The T2 heart extrapolation window equals T1's (2026-10-15).
- External data for T3 is allowed within the windows (GSE208162 ruling, 2026-08-18).

---

**Subject:** Human Team: two s.10 questions (cross-task data use, test-board definitions)

We are a Human Team. Before building, we would like to confirm two points.

1. **Cross-task use of released competition data.** Task 1's E8.5 scRNA is measured at E8.5, which is a held-out
   validation stage of the Task 2 heart setting; likewise the Task 2/3 E8.75 and E9.5 MERFISH stages for Task 1.
   (a) In the final phase, once all validation answers are released, may any released file from any task be used for
   any task (e.g. the Task 1 scRNA as a reference for the Task 2 heart E12.5 and Task 3 models)? (b) Does the s.10
   held-out list stop applying to validation stages once their answers are released?

2. **Test-board definitions.** For the final-phase boards: (a) is the differential-expression reference for Task 1 E12.5
   and Task 2 heart E12.5 the E10.5 stage, and for Task 2 embryo E7.75 the E7.5 stage? (b) Is copy_last on those boards
   that same stage? (c) Will `panels/index.json` list the test boards with their cell bounds and truth/reference sizes?
   (d) Will final evaluation use veckit `main` at cc3859a (severity_slope, scale_log_ratio), or the 0.2.0 metric set
   from veckit PR #6 (severity_ratio, size_fidelity)? (e) What time zone applies to the 2 December deadline?

Thank you.
