# Release-day rules: plain summary for sign-off

Drafted by the AI assistant (Claude) on 2026-10-20 from the files committed at `eb946ac`. Not yet approved.

AI agents wrote these four rules; please approve or change each before the 2026-10-27 release (Task 3 is already
approved). Steps are from `docs/P3_RUNBOOK.md`, PLAN is `docs/PLAN.md`, s1 and s2 are a board's two scored
submissions, * marks a file new on release day, and "(inferred)" marks what I worked out rather than read.

All four rules use released validation answers, and their use for test predictions is still unconfirmed (assumption
A20). Without them the heart size cannot run validly (inferred).

## 1. Heart and embryo size

- **Decides:** the size of the E12.5 heart and E7.75 embryo clouds as RMS radius (typical distance of a cell from the
  centre; the scorer's measure). Submissions are rescaled to it.
- **Reads, and when:** step 4.1 (`run_p3.sh`, on the DGX), after the label map is committed. Heart: E8.25_late, E8.5*,
  E8.75, E9.5, E10.5*. Embryo: E6.75, E7.25, E7.5*, E8.0.
- **Settings:**
  - Heart growth uses core cells (cardiomyocytes and endocardium) and only intervals starting at or after E9.0
    (`heart_min_rate_start_t: 9.0`), which leaves E9.5-E10.5.
  - Forecast: the log-scale average of three growth rules (the latest rate, half that rate, and a rate that slows over
    time), which comes to 1.581 x the observed log growth. The rule-picking step (tolerance 0.15) needs backtests these
    stages cannot give, so all three are always averaged.
  - If the core's thinnest axis is within 30 degrees of z (the axis across sections) at E9.5 and E10.5, only the
    in-plane size grows and z keeps E10.5's spread. Otherwise growth is 3D.
  - Scope (whole cloud relative to core): H1 continues the trend, H0 holds it, H2 halves it. s1 = H1, s2 = H2.
  - Embryo: a monotone curve (PCHIP) through the four stages' log 3D sizes, with no extrapolation.
- **Outcomes:**
  - `valid` true: the JSON is committed and sizes both submissions.
  - `valid` false (no interval from E9.0 on, e.g. E10.5 unusable): the builder refuses it; no fallback is written
    down.
  - The preflight stops on a missing file, an unmapped label or uncommitted code.
- **Blind spots:**
  - Growth after E10.5 is untested; the average is safest only if later growth runs at 0.5-1x the observed rate.
  - Unknown test stacks: about 0.08-0.14 log error for the heart, up to 0.25 for the embryo.
  - Embryo-to-embryo noise, amplified about 3.2x.
  - Compliance: the agent that set these choices had seen anchor-derived validation sizes (now redacted, logged). No
    number entered the code, but the E9.0 cut and 30-degree switch followed the E9.5 result.
- **Files:** `experiments/coords/size_procedure.py`, `run_p3.sh`, `p3_preflight.py`, `REPORT.md`;
  `src/vec/t2_build.py` (`apply_size_procedure`); PLAN 4.2-4.3.

## 2. T1 direction readout

- **Decides:** whether M2 (interpolating between our E10.5 and external E14.5 hearts) is built. Per PLAN 4.1 it also
  picks the gene-change direction of M1, the main T1 pipeline: Feng's if it passes, damped momentum if not.
- **Reads, and when:** step 4.3 (DGX): T1 E10.5*, E8.5 and E9.5, plus the filtered Feng 2022 CD1 heart summary
  (`EXTERNAL_DATA.md` entry X002).
- **Settings:** "build" needs all three checks:
  1. Feng's E9.5-to-E14.5 direction beats cross-fitted momentum by at least 0.05 in de_score.
  2. It also does so in de_direction.
  3. In cardiomyocytes, Feng gets the direction right for at least 0.70 of the genes that change by at least 0.25 on
     the scorer's log expression scale ("sign agreement").

  Momentum is the E8.5-to-E9.5 change, built and scored on separate halves of the E9.5 cells. de_score and
  de_direction, the scorer's "did genes move the right way" metrics, are averaged over 2 folds on groups both
  predictors cover (each needs 20 T1 and 10 external cells). Labels are grouped by fixed name patterns; an override
  must be a committed file, chosen from the printed label list alone.
- **Outcomes:**
  - All pass: "build". Any fail, or a missing reading: "do not build".
  - Refused: uncommitted files, an unregistered Feng source, a changed list of held-flat genes (sex-linked and strain
    genes kept unchanged), or a rerun without `--overwrite`.
- **Blind spots:** whether a one-day direction, with surrounding tissue, carries over to a two-day jump into a
  chamber-only dissection. It checks which way genes move, not how far, and Feng is another lab's data.
- **Files:** `experiments/t1_direction/config.yaml`, `readout.py`, `dirlib.py`, `t1_stats.py`; PLAN 4.1.

## 3. Cell count and the density rerun

- **Decides:** how many cells each T2 test file contains.
- **Reads, and when:** the test board's cell limits, at build time. The rerun (step 4.5, about 20 minutes) uses heart
  E10.5* (reference E9.5) and embryo E7.5* (reference E7.25), with the validation boards' max_cells and seeds 0-9.
- **Settings:**
  - Count = max_cells (the board's upper cell limit), never below truth_cells (the scorer's truth-half size; min_cells
    if that is not published).
  - The rerun tries 0.25, 0.5, 0.75 and 1 x max_cells on "intermediate quality" predictions: random real cells with
    half their expression rows shuffled.
  - Keep max_cells unless a lower level beats it by more than 0.3 points (paired mean, the average within-seed
    difference) in at least 8 of 10 seeds (random repeats).
- **Outcomes:**
  - No level passes: max_cells. One passes: that fraction, on the matching test board (inferred).
  - A count outside the allowed range, or a board without max_cells, stops the build (never clipped).
  - If the rerun cannot run: max_cells (inferred).
- **Blind spots:** our real prediction quality, since the rerun samples the answer itself. The evidence reaches 8x
  truth_cells (3x on an 8,000-cell truth half); counts above the source pool repeat rows, untested.
- **Files:** `src/vec/t2_build.py` (`n_cells`); PLAN section 3 item 6; `experiments/density/REPORT.md`,
  `run_density.py`, `summarise.py`.

## 4. Label mapping (`vec.harmonise` on release day)

- **Decides:** each new label's lineage group (coarse cell family), core status (core sets the heart size) and T1
  chamber category, and how cells are relabelled across files. It also measures a Gata4 replicate noise bar and
  cross-stage transfer, which no rule acts on automatically (inferred).
- **Reads, and when:**
  - Step 3, before anything touches a new file: every new file's labels.
  - Step 4.4: Gata4 replicates 1* and 2* against wild-type E8.75, plus E8.25_late, E8.5*, E8.75, E9.5 and E10.5*.
- **Settings:**
  - A label is core only if its name says cardiomyocytes or endocardium. Epicardium and its derived cells,
    fibroblasts, smooth muscle, pericytes, blood and vascular endothelium are shell.
  - Ambiguous names: check marker genes and classifier calls against the nearest released stage, record the evidence,
    and default to shell ("surrounding" for T1). Commit before running.
  - Classifier: logistic regression on 50 principal components (summary axes of the gene panel), or
    k-nearest-neighbours if that cross-validates better. Cells under 0.5 confidence are "Unassigned", labels under 10
    cells are left out, and the deeper file (more counts per cell) is thinned to match.
  - Cross-file compositions: a top-call/EM range (each cell's best label counted, against a statistical re-estimate).
  - Between stages, a group is clean if its worse F1 (an accuracy score) is at least 0.8, partial if at least 0.5,
    and poor otherwise.
  - Noise bar: the larger of the two estimators' replicate-to-replicate differences. A group clears it if both
    replicates' shifts from wild type, under both estimators, point the same way and all exceed it.
- **Outcomes:**
  - Unmapped labels stop the size procedure and T1 classification; harmonisation skips those files unless
    `--allow-unmapped`.
  - Missing files are skipped and logged in `release_status.json`; errors exit non-zero.
- **Blind spots:**
  - Whether labels mean the same across files ("accuracy" is agreement with the organisers' labels and our map).
  - Cell types absent from the reference go unflagged, and the noise bar rests on two embryos.
  - A core/shell mistake at E10.5 is amplified: in log terms the E12.5 core forecast moves about 2.6 times as much
    (inferred: 1 + 1.581), and the submitted whole-cloud size by less.
- **Files:** `configs/lineage_map.yaml` (header rule), `configs/t1_chamber_classes.yaml`, `src/vec/harmonise.py`,
  `experiments/harmonise/run_harmonise.py` and `REPORT.md`, `scripts/label_markers.py`.

## What you are approving

1. **Size:** the three-rule heart average, the 30-degree z switch, s1 = H1, s2 = H2, embryo PCHIP, valid results only.
2. **T1:** M2 only if Feng beats momentum by 0.05 on both metrics, with 0.70 cardiomyocyte sign agreement.
3. **Count:** max_cells, unless the 0.3-point, 8-of-10-seed test picks a lower level.
4. **Labels:** the core/shell naming rule (default shell) and the fixed classifier, range, cut-offs and noise bar.

## Open choices

Current choice first, then the alternative.

1. **Heart growth:** the three-rule average (worst-case bias 0.581 of the observed growth, noise amplification 3.187)
   or half rate alone (1 and 2.449).
2. **Heart scope:** H1 and H2, or H0 (hold) in one submission.
3. **Invalid heart result:** no fallback, or one written now (for example, keep E10.5's size).
4. **Epicardium:** shell or core.
5. **Embryo:** only record the z risk, or grow in-plane if E7.5 is cut differently from E8.0.
6. **T1:** all three checks or a looser test; should a fail also drop Feng from M1 (PLAN: yes; config: silent)?
7. **Count wording:** "mean gain above 0.3, with 8+ seeds ahead" (my reading) or "8+ seeds each ahead by 0.3"; also
   which level if several pass, and per board or both.
8. **Groups rated poor on new stage pairs:** no rule covers them. Keep PLAN 4.3's list, or drop them.

## Recommended answers (AI assistant, for your sign-off)

1. **Heart growth:** keep the three-rule average. Its worst-case bias (0.581 of the observed growth) is well under half
   rate's (1.0), for a little more noise (3.19 against 2.45).
2. **Heart scope:** keep s1 = H1 and s2 = H2. No new evidence favours H0, and both stay inside the shared core.
3. **Invalid heart result:** write a fallback now. If the procedure returns `valid` false while E10.5 is usable, keep
   E10.5's size, so the scale term sits at the floor instead of the build failing; s1 and s2 still differ by scope.
   If E10.5 is barred (A20 fails), the heart plan is rebuilt from E9.5 anyway and the size rule comes back to you then.
4. **Epicardium:** keep it in the shell. It is a thin surface layer; the core is the tissue whose growth the procedure
   measures.
5. **Embryo:** keep recording the z risk only. An in-plane rule would be a new free choice with no test behind it.
6. **T1:** keep all three checks, and write down what PLAN 4.1 already says: a fail also drops Feng from M1, which then
   uses damped momentum.
7. **Count:** a lower level replaces max_cells on a board only if, on that board's validation target, its paired mean
   gain over max_cells exceeds 0.3 points and it beats max_cells in at least 8 of the 10 seeds. If several pass, the
   largest mean gain wins (ties go to the higher count). Heart and embryo are decided separately, each for its own test
   board. If the rerun cannot run, keep max_cells. The rule goes into code so nobody reads it off a table.
8. **Poor groups:** keep groups rated poor between the new stage pairs out of the composition trends (their share holds
   at E10.5), as PLAN 4.3 already does for CardiacProgenitor and Pericardium. The harmonisation noise bar stays
   informational; Task 3 computes its own.
