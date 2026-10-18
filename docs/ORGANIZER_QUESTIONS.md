# Questions for the organisers

Send by email (virtual.embryo.moonshot@gmail.com) or Slack `#virtual-embryo-challenge` **before** relying on any of
these. Rules s.10: "an answer given in advance is not a finding of misconduct afterwards". Record each answer with its
date in `docs/DATA_POLICY.md`.

Status: **draft, not sent** (2026-10-18).

---

**Subject:** Human Team: rules s.10 questions on Task 3 priors, external data for Tasks 2/3, and cross-task data use

We are a Human Team. Before building, we would like to confirm four points under s.10.

1. **Published phenotype findings for the beta-catenin board.** The test condition is Mesp1-Cre; Ctnnb1 F/F at E8.75.
   The literature contains in situ, immunostaining and qPCR results from closely related genotypes (e.g. Klaus et al.,
   PNAS 2007, doi:10.1073/pnas.0703113104, the same driver at 6-10 somites; Lin et al. 2007, Ai et al. 2007 and
   Cohen et al. 2007 with other cardiac Cre drivers at E8.0-E9.5). We read s.10 ("another allele of the same gene at a
   comparable stage") as prohibiting their use, even as a qualitative gene-to-direction prior, and will proceed on
   that basis unless you say otherwise. Does the prohibition also cover (i) curated Wnt target-gene lists where some
   entries came from those experiments, and (ii) review articles restating them? We assume (iii) sequence resources
   such as TCF/LEF motif databases, and generic pathway membership not derived from these genotypes, are fine.

2. **External data for Tasks 2 and 3.** The Rules page says external public data may be used (with disclosure and the
   stage windows). The Data page says in places that "Task 1 is the only task where training on external public
   single-cell data is allowed". Which wording controls? Specifically, may we train on external public mouse embryo
   single-cell or spatial data for Task 2 and Task 3, subject to the held-out windows and disclosure?

3. **Cross-task use of released competition data.** Task 1's E8.5 scRNA is measured at E8.5, which is a held-out
   validation stage of the Task 2 heart setting; likewise the Task 2/3 E8.75 and E9.5 MERFISH stages for Task 1.
   (a) During the final phase, once all validation answers are released, may any released file from any task be used
   for any task (e.g. Task 1 scRNA as a reference for the Task 2 heart E12.5 and Task 3 models)? (b) Does the s.10
   held-out list stop applying to validation stages once their answers are released?

4. **Test-board definitions.** For the final-phase boards: (a) is the differential-expression reference for Task 1
   E12.5 and Task 2 heart E12.5 the E10.5 stage, and for Task 2 embryo E7.75 the E7.5 stage? (b) Is copy_last on those
   boards that same stage? (c) Will `index.json` list the test boards with their cell bounds and truth/reference sizes?
   (d) Will the scorer used for final evaluation be veckit `main` at cc3859a (severity_slope, scale_log_ratio), or the
   0.2.0 metric set from veckit PR #6 (severity_ratio, size_fidelity)?

Thank you. We would rather have the answers in advance than argue them afterwards.
