# Method summary (running draft)

Required for prize eligibility: winners must deliver a written report within 14 days of notification (approach, design
decisions that mattered, what was trained on what, enough detail to understand why it performs as it does). The
Generality Award is judged on this summary. Keep it current as methods change; every claim should point to a run in
`results/registry.csv` or an experiment report in `experiments/*/REPORT.md`.

Status 2026-10-28: methods are pre-registered on released training data; no test submission exists yet. The final phase
opened on 2026-10-27 (postponed by the organisers from 2026-10-20), and the team unlocked the P2 validation answers on
2026-10-28 (Data provenance, below). Sections marked _planned_ describe the committed plan (`docs/PLAN.md` v1.2), not
results.

## Team and track

Team Mithril, Track 1 (Human Team). A person reviews results and makes every submission decision; coding agents do
the work.

## Data provenance and compliance (applies to every task)

- **Released competition data only, plus disclosed legal external data** (`EXTERNAL_DATA.md`). Nothing measured inside a
  protected window is used, by any route; windows and the do-not-use list are in `docs/DATA_POLICY.md`.
- **Validation answers (since 2026-10-28).** The P2 validation answers (T1 E10.5, T2 heart E8.5 and E10.5, T2 embryo
  E7.5, the two Gata4 knockout replicates) are opt-in in the final phase; the team captain unlocked them on the Data
  page on 2026-10-28. Rules s.10 makes them an exception to the held-out-data restrictions solely for preparing P3
  predictions; the team reads that as including their use as model inputs at prediction time (for example E10.5 as
  the starting stage for E12.5) and has asked the organisers to confirm. External data from those stages and genotypes
  stays out, and no test inputs are distributed.
- **No board anchors and no leaderboard scores set any number.** Test-board anchors are stripped in memory when the
  board specs are refreshed (`vec.boards.update_boards`, or `update_boards_from_file` for a board index the user
  saved in the browser); the public board file is never opened.
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
  binomial thinning or Poisson top-up on recovered integer counts, then renormalised to log1p CP10k. No decoded means,
  no additive log-space deltas (they break the sparsity-driven variogram).
- **Label harmonisation** (`vec.harmonise`): logistic regression on 50 PCs of the shared panel, group probabilities
  summed from fine labels, depth matching; cross-file compositions quoted as a top-call/EM bracket
  (`experiments/harmonise/REPORT.md`).
- **Pre-registration:** every rule that turns released data into a submission is committed before the data it acts on
  (size procedure, T3 decomposition rules, T1 direction readout, cell-count rule, label-mapping rule).
- **Provenance:** every prediction file carries a build stamp (its sha256, the commit that built it, and whether the
  code was unchanged from the build's start to the write); an upload is accepted only with a matching clean stamp,
  and the submission ledger records the building commit of every file sent.

## Task 1: temporal (E12.5) — _planned_

- **Scope first.** E12.5 is a heart-chamber-only dissection; re-weighting real E10.5 cells to a chamber-only composition
  is the largest single lever: in an E8.5 -> E9.5 pseudo-task on released data, re-weighting to the target's exact
  category shares gained about +31 points over copying the earlier stage (an upper bound, since the E12.5 shares must be
  chosen; `experiments/t1_scope/results/`). The cardiomyocyte share is the main unknown and the hedge between the two
  submissions; there, a share 10 points too high cost about 3 points and 10 points too low under 1.
- **Per-type count-space shifts.** Direction from within-dataset fold changes of legal external chamber data (Feng 2022:
  E9.5 and E14.5 onward) if a pre-registered readout on the E9.5 -> E10.5 transition shows it beats cross-fitted
  momentum; damped internal momentum otherwise. Sex-linked and strain-polymorphic genes are held flat.
- **Emergent populations** by the legal route (classes present in legal E14.5+ data but absent at E10.5), pending the
  organisers' answer about qualitative background.
- Training data: released T1 E8.5, E9.5, E10.5; Feng 2022 legal stages (X002). JAX/Qiu 2024 for diagnostics only (X003).

## Task 2: spatial-temporal (heart E12.5, embryo E7.75) — _planned_

- **Heart E12.5:** real E10.5 tissue (synthetic records dropped), sampled at random to the pre-registered count (the
  board's maximum cell count unless the release-day density rerun picks a lower level by a rule fixed in advance,
  never below the board's lower limit; `vec.t2_build.n_cells_from_decision`);
  damped composition trends on harmonised merged lineages; coordinates grown by the pre-registered size procedure
  (core = cardiomyocytes + endocardium; growth from E9.5 -> E10.5 only; a fixed three-rule average justified without an
  assumed rate; z held at E10.5 when the sectioning plane matches); two scope hypotheses as the two submissions.
- **Embryo E7.75:** interpolation between E7.5 and E8.0 (fused Gromov-Wasserstein coupling, proper-rotation alignment,
  displacement interpolation at a calibrated clock), size from a PCHIP curve through the released embryo stages.
- **Spatial hygiene:** random rows only, jitter of at most 0.25 nearest-neighbour distances to break duplicate
  coordinates, proper rotations only, native coordinates (`experiments/density/REPORT.md`, `experiments/coords/REPORT.md`).

## Task 3: perturbation (beta-catenin KO, E8.75) — _planned_

- **Per-gene factor first.** On the released Mab21l2 knockout, the response is mostly one per-gene factor shared by every
  cell type: estimated from non-Mesp1 compartments and applied to all wild-type cells it scored 80.5 against 50.8 for
  wild-type resampling, while rules built from the broad phenotype alone scored about 51
  (`experiments/t3_decomp/results/mab21l2/`).
- **Composition second.** Second-heart-field, outflow-tract and right-ventricle-lineage states depleted at a magnitude
  calibrated on the released Gata4 replicates, counted only where it exceeds the replicate difference.
- **Decided by pre-registered rules** (`experiments/t3_decomp/rules.yaml`) run once on the released Gata4 files; the two
  submissions differ in how the factor is sourced. Because the factor is measured in compartments the knockout cannot
  touch, it may partly reflect sample processing: unless the Mab21l2 factor transfers to Gata4, the two submissions do
  not both rely on the Gata4 factor, and if no check supports it, the first submission is unchanged wild-type cells
  (with the composition term) and only the second carries the factor, at half strength. Released Gata4 files are used
  only if the organisers confirm it (2026-10-28: the team reads the final-phase Rules s.10 as allowing it and has asked
  them to confirm; until they answer, both versions are built, with and without the Gata4 files). The rules were
  drafted with an AI assistant and reviewed and approved by the team before any Gata4 file was downloaded.

## Shared architecture (Generality Award)

_Planned claim_: the summary will ask to be considered for the award unless decision point D3 (below) says otherwise.

- **Common core** (`src/vec/core.py`, feat-020). Every prediction is built from real released cells by one
  "transport-then-adjust" operator:
  - allocate a target composition to whole cell counts;
  - draw real cells from the source pools;
  - edit their expression in count space (per-type or per-gene factors by binomial thinning or Poisson top-up, a
    depth change);
  - write log1p CP10k.

  The planned steps go into the same module: coupling and displacement interpolation between bracketing stages (T2
  embryo, the T1 external bracket), births from a legal source (T1), and leave-one-transition-out calibration of the
  few shrinkage parameters.
- **Task heads** (what differs):
  - T1 (`vec.t1_scope`): the 32,285-gene encoder, chamber categories and target shares, and the direction (momentum or
    external).
  - T2 (`vec.t2_build`): the 500-gene panel and the coordinate head (cell count, random rows, jitter, the
    pre-registered size procedure, proper rotations).
  - T3 (`experiments/t3_decomp`): the genotype head (knockout composition from harmonised labels, and a per-gene
    factor estimated in compartments the knockout cannot touch).
- **Boundary:** heads supply the source pools, the target composition, the factors and the coordinate rules; the core
  turns them into cells. The T3 script was written first, then pre-registered and hash-pinned, so it keeps its own
  primitives; the core holds identical copies, and `tests/test_core.py` proves them identical (same code, same
  outputs and random draws). Golden tests show the move into the core changed no T1 or T2 output
  (`experiments/core/REPORT.md`).
- **Decision point D3 (2026-11-12):** claim the award unless specialised variants beat the shared-core builds by more
  than 3 points of overall standing on honest validation. Both submissions on every board are built through the core.
