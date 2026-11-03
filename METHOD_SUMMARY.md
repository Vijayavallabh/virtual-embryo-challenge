# Method summary (running draft)

Required for prize eligibility: winners must deliver a written report within 14 days of notification (approach, design
decisions that mattered, what was trained on what, enough detail to understand why it performs as it does). The
Generality Award is judged on this summary. Keep it current as methods change; every claim should point to a run in
`results/registry.csv` or an experiment report in `experiments/*/REPORT.md`.

Status 2026-11-01: the entry for each board is chosen (below) and built end to end; no test submission exists yet (the
test boards' specs wait on the board intake, not before 2026-11-04). The final phase opened on 2026-10-27, and the team
unlocked the P2 validation answers on 2026-10-28 (Data provenance, below). Each change was kept only if it passed the
adoption rule on pseudo-tasks built from released stages (Shared machinery), unless a pre-registered rule decided it.

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
  the starting stage for E12.5) and asked the organisers to confirm (21 and 28 October), with no reply by 2026-10-31.
  The team proceeds on that reading and on the predicted organiser answer (`docs/ORGANIZER_QUESTIONS.md`, working
  answers), not on a ruling; the Task 3 use of the Gata4 replicates rests on the same reading. External data from
  those stages and genotypes stays out, and no test inputs are distributed.
- **No board anchors and no leaderboard scores set any number.** Test-board anchors are stripped in memory when the
  board specs are refreshed (`vec.boards.update_boards`; `update_boards_from_file` does the same for a board index
  saved in the browser, a route not cleared for use); the public board file is never opened.
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
  with the calibrated emulator (the server's own draws reproduced; all 66 validation anchors exact), the floor checked
  to score 50 on every seed, candidates built from the allowed stages only (no file paths, a read guard against files
  derived from targets) on shared draws, and a change kept only under the adoption rule (approved rule 7's test per
  pseudo-task).
- **Pre-registration:** every rule that turns released data into a submission is committed before the data it acts on
  (size procedure, T3 decomposition rules, T1 direction readout, cell-count rule, label-mapping rule). The s2 rules
  (`docs/S2_RULES.md`) were written by the AI assistant under the captain's delegation of 2026-11-02 from sensitivity
  analyses on released data, before any test score, and frozen on 2026-11-03. The heart rule was amended after its
  released-data outcome, and the amendment is logged as such.
- **Tries ledger** (since 2026-11-03, `experiments/harness/TRIES.md`): every change tried for an entry is written
  down before its run with its physical reason, the metrics where its gain should show and what would make a gain a
  fit to one backtest, and every try is counted. Tried and rejected on 2026-11-03: matching detection depth by
  binomial thinning (released stages differ several-fold in counts per cell: the embryo series rises from 30 to 212,
  the heart series falls from 233 to 56). It lost on the more test-like backtest of both boards (-4.18 embryo, -3.92
  heart), lowering mmd_u and neighbourhood MMD, so it is not used. Blending coupled pairs of cells also lost on both embryo backtests:
  averaged cells are too narrow a population. A target-free diagnostic found that the FGW coupling cannot tell the
  embryo's left side from its right (expression is left-right symmetric at these stages), so its partners can be
  mirrored and the moved cloud flattens; it happens in every draw of the test pair. Taking the moves from an exact
  transport on expression and position in the rotated common frame fixed every draw and gained +2.13 on the E7.5
  backtest (shape only) but lost -1.59 on E7.25, where its shorter moves cost shape, so it is not used; the s2
  (E7.5's geometry) hedges the flattened s1 geometry.
- **Provenance:** every prediction file carries a build stamp (its sha256, the commit that built it, and whether the
  code was unchanged from the build's start to the write); an upload is accepted only with a matching clean stamp,
  and the submission ledger records the building commit of every file sent.

## Task 1: temporal (E12.5) — M1

- **Scope.** E12.5 is a heart-chamber-only dissection, so the entry draws real E10.5 cells to E10.5's own composition
  restricted to the chamber categories (`configs/t1_chamber_classes.yaml`): blood 24% of chamber cells, cardiomyocytes
  26.5% of heart cells (the captain's choices of 2026-11-01). E10.5 is the same nuclei assay on dissected hearts; the
  earlier stages were embryo portions, and their compositions do not transfer (blood 1.25% at E9.5; cardiomyocytes 66%
  of heart cells). Legal external data agree on not raising the cardiomyocyte share: Feng's own share barely moves from
  E9.5 to E14.5. The s2 hedge: cardiomyocytes at 60% of heart cells, decided from the T1 grid
  (`docs/S2_RULES.md`).
- **Direction.** Cells of the groups Feng 2022 profiles (cardiomyocytes, endocardium, endothelium, epicardium) move in
  count space along Feng's within-dataset E9.5 -> E14.5 log ratio of mean CP10k, scaled by the step's share of that
  bracket (E10.5 -> E12.5: 2/5); the other cells (blood, mostly) move by half the E9.5 -> E10.5 per-type change. A
  pre-registered readout on the E9.5 -> E10.5 transition chose Feng over cross-fitted momentum (de_score +0.175,
  de_direction +0.188, cardiomyocyte sign agreement 0.91; `experiments/t1_direction/results/E10.5/`), and the captain
  ruled that this check governs the direction, since no T1 pseudo-task can judge it (Feng's E9.5 cells lie inside each
  one's external window). A labelled end-to-end diagnostic in the readout's geometry: +1.29 and +1.09 board points over
  B0, 10 of 10 seeds (`experiments/t1_m1/results/diag_v1/`). Sex-linked, Y-linked and strain-variable genes are held
  flat (a list fixed from E8.5 and E9.5 before E10.5 was released).
- **What the T1 scorer rewards.** On the released E9.5 -> E10.5 transition the changed genes are mostly highly
  expressed genes going down, which de_score's expression-magnitude null names 39-41% of; a copy names 11%
  (`experiments/harness/results/t1_de_v1/`). Direction quality is what moves de_score above its floor.
- **Emergent populations:** none in the entry (the organisers have not answered on qualitative background).
- Training data: released T1 E8.5, E9.5, E10.5; Feng 2022 E9.5 and E14.5 CD1 hearts (X002), reduced from GEO by
  committed code that drops E10.5-E13.5 cells before anything is written. JAX/Qiu 2024 for diagnostics only (X003).

## Task 2: spatial-temporal (heart E12.5, embryo E7.75)

- **Heart E12.5 — B0 with the pre-registered size:** real E10.5 tissue (synthetic records dropped), sampled at random
  to the pre-registered count (the board's maximum unless the release-day density rerun picks a lower level by a rule
  fixed in advance; `vec.t2_build.n_cells_from_decision`), moved in count space by half the E9.5 -> E10.5 change of
  each harmonised lineage; coordinates grown by the pre-registered size procedure (core = cardiomyocytes +
  endocardium; growth from E9.5 -> E10.5 only; a fixed three-rule average justified without an assumed rate; z held at
  E10.5 when the sectioning plane matches); the two submissions are the two scope hypotheses (H1; H2 as the half trim, amended after its released-data outcome,
  `docs/S2_RULES.md`). E10.5's
  composition is held: E9.5 and earlier were embryo portions and E10.5 a dissected heart, so a composition trend through
  them would carry that change of sampling on to E12.5 (dropped by the captain, 2026-11-01). On the pseudo-tasks B0
  beats a copy (+0.47, +0.90) and the floor (+0.43, +1.45); the size cannot be backtested.
- **Embryo E7.75 — displacement interpolation (FGW):** a balanced fused Gromov-Wasserstein coupling of 2,000 cells per
  bracket (E7.5, E8.0) on 30 expression principal components and within-cloud distances (structure weight 0.5), the
  proper rotation that maps one bracket onto the other under the coupling, every cell moved by its nearest anchors'
  displacements t (or 1 - t) of the way, at the nominal time fraction t = 1/2; the population is a mixture of
  unedited real cells of both brackets (half each), scaled to the pre-registered embryo size (PCHIP of log RMS radius
  through the released embryo stages). On the pseudo-tasks: 61.12 on E7.5 against 56.97 for the organisers'
  pseudobulk_shift and 58.83 for it with the size rule, and 58.78 on E7.25 against 53.86 for B0, every gap ahead in
  10 of 10 seeds (`experiments/harness/results/emb_fgw_v1/`). It loses on shape on E7.5, where E8.0's sampling differs;
  the s2 is fgw_shape0, fgw's mixture in E7.5's geometry, chosen by a rule fixed before its run (`docs/S2_RULES.md`).
- **Spatial hygiene:** random rows only, jitter of at most 0.25 nearest-neighbour distances to break duplicate
  coordinates, proper rotations only, native coordinates (`experiments/density/REPORT.md`, `experiments/coords/REPORT.md`).

## Task 3: perturbation (beta-catenin KO, E8.75) — pre-registered

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
- **Outcome (2026-10-28, run once):** the Gata4 replicates transfer to each other, the Mab21l2 factor does not transfer
  to Gata4. s1 = the allowed branch's first submission; s2 = the refused branch's first (wild-type cells, the central
  composition at the Mab21l2 magnitude, no number from the Gata4 files), uploaded first (approved 2026-11-01). The
  harness gives these entries no verdict (every fold is one they were selected or fitted on); they stand on the
  pre-registered rules.

## Shared architecture (Generality Award)

_Planned claim_: the summary will ask to be considered for the award unless decision point D3 (below) says otherwise.

- **Common core** (`src/vec/core.py`, feat-020). Every prediction is built from real released cells by one
  "transport-then-adjust" operator:
  - allocate a target composition to whole cell counts;
  - draw real cells from the source pools;
  - edit their expression in count space (per-type or per-gene factors by binomial thinning or Poisson top-up, a
    depth change);
  - write log1p CP10k.

  - couple two populations with positions and move cells between them (fused Gromov-Wasserstein, a proper rotation
    under the coupling, barycentric partners, nearest anchors; added 2026-11-01 for the T2 embryo).

  Births from a legal source and leave-one-transition-out calibration are not in the entries.
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
