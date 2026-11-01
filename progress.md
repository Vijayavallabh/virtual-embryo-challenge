# Session Progress Log

## Current State

**Last Updated:** 2026-11-01
**Active Feature:** feat-018 (s1 files): the s1 build path is reviewed and fixed (2026-11-01, 46089fe); feat-004's
step 1 waits for the board-spec save (not before 2026-11-04 UTC), then the s1 builds; the s2 hedge payoffs
(`experiments/s2_hedge/`) feed the s2 pre-registration (17-23 Nov). Earlier: feat-004 - Final-phase intake (`docs/P3_RUNBOOK.md`). Steps 2 (data) and 3 (labels) done on
2026-10-28 (feat-003 done), and step 4 (the pre-registered runs) ran the same day at the user's go-ahead, each output
committed before it was read. Step 1 is blocked by the 403 on the panel files (Cloudflare error code 1010, a
browser-signature block; no reply to any organiser email). Its saved-file intake (`vec boards update --from-file`,
cee1c92) is cleared by the captain (2026-11-01): he sends the drafted notice, then saves the file not before
2026-11-04 UTC (`docs/P3_RUNBOOK.md` step 1). On
2026-10-31, at the captain's instruction, the open organiser questions got predicted working answers
(`docs/ORGANIZER_QUESTIONS.md`; T3 selects the allowed branch), and the emulator calibration (feat-005) ran: the
scorer reproduces every published validation anchor exactly. The release-day rules are signed off and in code
(`docs/RELEASE_RULES_SUMMARY.md`).

## Status

### What's Done

- [x] Read the full challenge site, rules, data page, proposal and the veckit scorer; condensed in `docs/CHALLENGE.md`.
- [x] uv project (Python 3.12), harness, `vec` package (emulator, checks, baselines, safe downloader), DGX environment.
- [x] Research sweep (10 reports in `docs/research/`), PantheonOS Slack read in full, plan in `docs/PLAN.md`.
- [x] `ASSUMPTIONS.md` (19 rows, linted), data policy with organiser rulings and a do-not-use list.
- [x] Leaked E12.5 heart size quarantined: redacted from notes; private report to organisers sent on Slack.
- [x] Team API key stored at ~/.config/vec/token (mode 600, local only); 9 training files downloaded to the DGX and checksummed.
- [x] 2026-10-19: experiments on released data, each independently reviewed (`experiments/*/REPORT.md`): density and
      sampling (count rule max_cells, never below truth_cells; random rows only), coordinates and the pre-registered heart
      and embryo size procedure (committed 171695a before the release), label harmonisation (`vec.harmonise`), data QC
      (`vec.spatial`, `vec data qc`; copied and interpolated records dropped from prediction inputs).
- [x] 2026-10-19: research sweep screened by a separate agent; PLAN v1.2 built from the screened digest only
      (`docs/research/digest_2026-10-19.md`). Two meaning-level redaction passes on the notes; leaked-scorer copies,
      anchor-derived outputs and in-window external files deleted (user's decision).
- [x] 2026-10-19: test-board anchors quarantined at intake (`vec.boards.update_boards`, `vec boards update`);
      release-day runbook written. History rewritten twice and the GitHub repo recreated; organiser follow-up sent.
- [x] 2026-10-19: T3 rules notes 6 and 7: unless cross-knockout transfer (D11b) passes, s1 and s2 do not both use the
      Gata4 factor; note 7 keeps the D8 hedge in the replacement s2. The user read a plain summary of the T3 rules,
      approved them and chose note 8 (no check supports the Gata4 factor -> s1 = WT cells, s2 = the factor at half
      strength). All committed before any Gata4 file is downloaded; `build_submission.py plan` shows the result.
- [x] 2026-10-20: the organisers moved the final phase to 2026-10-27 (email). The two code paths that could read
      held-out validation material now wait for the final phase itself (the T1 composition at `vec boards update`; the
      anchor check); runbook, plan timeline, data policy, ledger (A20) and handoff updated; Yukon team material
      quarantined. Slack and the organiser email swept by screened agents (findings under Blockers).
- [x] 2026-10-20: the upload path tested against the live server while the validation phase runs a week longer: three
      format checks pulled from the DGX (T1:val, T2:heart:val_extrap, T3:gata4; practice files from
      `vec baseline --n-cells 1000`) each came back format_ok with no submission id; before and after, no server entry
      on any board and the daily scored and upload counts at 0 (feat-016 evidence; ledger in results/submissions.csv).
- [x] 2026-10-20 (organisers' email of 22:03 IST): a collaborative Yukon track (Eigen Labs; invite-only; separate
      $10K-$30K pool; collaborative phase ends before P3; one-way transfer; ineligible for the official prizes) and the
      deadline's time zone (through 2 December 2026, UTC). An agent read the official Rules, Terms, Privacy,
      partner-rules, Timeline and FAQ pages (nothing quarantined opened): ordinary-team rules are unchanged; official
      P3 scores are published when each submission is made; the P3 privacy rule (methods, prediction files and
      evaluation feedback private until 2026-12-11) is on the partner page only, and we follow it. Challenge
      reference, data policy, plan, organiser-questions log and AGENTS.md updated; the organiser draft cut to two
      questions.
- [x] 2026-10-20: prize strategy, approved by the user ("do whatever's best"; `docs/PLAN.md` decision log 2026-10-20
      (c)). Stay official; go for a Human-track place plus the Generality Award, which stacks on a place and is judged
      across both tracks among shared-architecture entries. D3 loosened to "claim unless specialised code wins by more
      than 3 points overall"; the shared core goes into code first (feat-020); ledger A16 revised, A21 added. The
      Generality question was added to the organiser draft (now three questions). Prize terms and a leaderboard census
      (323 Human teams, 159 above the floor, 6th place 223.5; Agent track 80) recorded in `docs/CHALLENGE.md` and PLAN
      section 1. Not pursued: Yukon; a switch to the Agent track.
- [x] 2026-10-21: feat-020, the shared transport-then-adjust core `src/vec/core.py` (Generality Award claim). The
      T1 and T2 heads' allocation, draws and count-space edits moved into it, source unchanged except names; the
      pre-registered T3 script keeps its own primitives and the core holds identical copies. Evidence
      (`experiments/core/REPORT.md`):
      - `tests/test_core.py`: 31 tests, including 12 golden cases recorded at eb946ac and T3 copies identical in
        body and on 300 random inputs per function;
      - released data on the DGX: old against new, all 64 hashed parts identical;
      - `./init.sh`: 320 passed;
      - an independent verifier: PASS. Its fixes were made: the module was renamed from `vec.operator`, which
        shadowed the standard library, and the wording about T3 was corrected.
      Release-day rules summary drafted for sign-off (`docs/RELEASE_RULES_SUMMARY.md`, with recommended answers).
- [x] 2026-10-21: the two approved answers that needed code (answers 3 and 7 of `docs/RELEASE_RULES_SUMMARY.md`),
      built by an agent and checked by three independent verifiers (fidelity, integrity, integration) over two fix
      rounds; every code finding was fixed:
      - `vec.t2_build.size_or_hold`, the heart size step: a valid JSON is applied exactly as `apply_size_procedure`
        applies it; `valid` false keeps E10.5's size; cells from another stage get no size (back to the user); every
        heart call's reference cloud must be the cloud the JSON measured as E10.5.
      - `count_override`, `read_count_decision` and `n_cells_from_decision` in `vec.t2_build`, and
        `experiments/density/summarise.py --decide`: the count decided per board from the release rerun, written once
        to `experiments/density/results/release/count_decision.json`; uncommitted rule code and a second decision
        without `--overwrite` are refused.
      - `docs/P3_RUNBOOK.md` (steps 4.1 and 4.5, section 5) and `run_p3.sh`'s closing note now name them.
        `./init.sh`: 355 passed (35 new tests).
- [x] 2026-10-21: the case approved answer 7 left open, in code as the user approved it ("the next passing count the
      board allows, or max_cells if none does"; read as the next count up, the reading the user left to the AI
      assistant and a decision panel chose, `experiments/density/floor_reading/README.md`):
      - `vec.t2_build.n_cells_from_decision` passes over a decided level below the test board's lower limit to the
        passing level with the fewest cells that the board allows, else max_cells, and records `count_floor`
        (`COUNT_FLOOR_RULE`);
      - `read_count_decision` now also checks the decision's record: the decided level heads its passing levels,
        every level's verdict agrees with its own numbers (both ways), and the rerun was complete;
      - runbook section 5, the rules summary, PLAN 3.6 and 4.2 step 6, the decision log, METHOD_SUMMARY and the
        ledger log updated; the third organiser email recorded as sent (`docs/ORGANIZER_QUESTIONS.md`).
      Independent verifiers in two rounds on the gain reading (fidelity, robustness and docs; then code and docs
      re-checks): every finding fixed; `./init.sh`: 357 passed (3 new tests, 1 replaced; d06727f). The reading of
      "next" was then put to the user, who left it to the AI assistant; a four-agent panel (simulation, data, an
      advocate for the gain reading, judge) chose the next count up, medium confidence. The switch to it: a code
      verifier (pass; 0 mismatches in 579,150 cases against an independent implementation) and a docs verifier, their
      findings fixed; `./init.sh`: 357 passed.
- [x] 2026-10-19: submission client and builders, each built by an agent and verified by an independent agent, then
      three rounds of independent review and fixes (every finding verified fixed):
      - `vec submit` (feat-016): format-only by default; official uploads behind --official, a hash confirmation and
        every gate in `src/vec/submit.py`; `vec submit status` is read-only.
      - `vec.t2_build`: count rule with the truth floor enforced, random rows, jitter, size procedure, cloud guard;
        smoke test on both validation boards (`experiments/t2_build/results/`).
      - `vec.t1_scope`: chamber categories for every released T1 label (`configs/t1_chamber_classes.yaml`),
        re-weighting to input shares, momentum shift; E8.5 -> E9.5 pseudo-task (`experiments/t1_scope/results/`).
      - `experiments/t3_decomp/build_submission.py`: s1/s2 from the decision file, no re-deciding.
- [x] 2026-10-27 to 2026-10-28 (details under Blockers / Risks): the final phase opened and runbook step 0 passed; the
      history of `main` was rewritten without the AI co-author trailer (the user's decision, `docs/COMMIT_MAP.md`);
      the user unlocked the validation answers, the data client reads them (e0e1d58) and step 2 is done (f423f6a;
      feat-003 done); step 3 (labels) done at the user's go-ahead (aebcb8c, d6b9566); A1 passed, and A20 passed on
      the user's reading (d47badf); the fourth organiser email, drafted by the AI assistant, was sent by the user
      (ed3b1bc). `./init.sh`: 366 passed.
- [x] 2026-10-28: runbook step 4 at the user's go-ahead ("build the intake and start step 4"), in order, on the DGX at
      725d19c with a clean tree, each output committed before it was read (the rules decide; nothing chosen by eye):
      - 4.1 size (f5be247): both JSONs valid. Heart E12.5 in the pre-registered xy_grow_z_hold mode: full RMS 614.7
        for s1 (H1) and 560.5 for s2 (H2), core 487.2; embryo E7.75 in raw mode: 260.3. Handedness: E8.5 and E10.5
        resolve right-handed by the Pitx2 landmarks.
      - 4.2 T3 Gata4: replicate overlap 0 shared rows of 42,257 and 46,396, so independent (95d8198); decision file
        (7e7c207), aggregates (3ff6bf0), plan (05ead40). Rule outcomes: D2, D2u, D4, D7, D7m and D11a pass; D1, D5,
        D6, D8, D9, D10 and D11b do not. Plan: allowed branch s1 = Gata4 factor at exponent 1 with the central
        composition (log2 1.383) and depth, s2 = WT cells with the same composition (note 6: D11b failed, D7
        passed); barred branch s1 = Mab21l2 factor with the central composition (log2 0.718), s2 = the same factor
        with the mild one (0.359); refused branch s1 = WT cells with the central composition, s2 = plain WT cells.
      - 4.3 T1 direction readout (3c5cc66, default config and label map; the label check mapped every E10.5 label):
        build M2. External direction minus cross-fitted momentum: de_score +0.175, de_direction +0.188 (bar 0.05);
        CM sign agreement 0.912 (bar 0.70).
      - 4.4 harmonisation release stage (0586bf6): Gata4 rep1 vs rep2 vs WT E8.75 (the replicate noise bar for T3
        composition effects, A19; e.g. IFT-CM falls by 0.03-0.04 in both replicates, 10x the replicate difference)
        and cross-stage transfer along E8.25_late-E10.5. The script's default Gata4 file names differ from the
        released ones, so the run passed the real paths.
      - 4.5 density rerun (dfee546 decision, ef1819f tables): 20 of 20 runs complete, none skipped; approved rule 7
        keeps max_cells on both test boards (no lower level gains more than 0.3 points in 8 of 10 seeds).
      - 4.6 scorer-build anchor check (12151e0, validation anchors only): the bake protocol is reproduced exactly on
        all three T2 validation boards (numpy Generator seed 0, rows sorted, veckit split_half seed 0, truth = first
        half, floor = the reference's working copy, metric seed 0; every d2_shape, Dice and scale_log_ratio anchor
        hit). The six Dice anchors are proper-mode (5) or tied (1), none mirror-mode; both builds agree on these
        draws, so the anchors do not say which build baked them. Feeds emulator calibration (feat-005).
- [x] 2026-10-28: saved-file board intake at the user's request (cee1c92): `vec boards update --from-file PATH`
      (`vec.boards.update_boards_from_file`) for a board index the user saves with the browser ("Save link as", never
      opened) outside every git work tree; the same core as the online refresh, all-or-nothing writes, typed test
      specs, field names dropped from test boards reported for D1 (names only). An independent review by three AI
      agents found ways a mistaken file in the saved folder could carry anchors into `boards/`; all refused now, each
      with a test that fails when its fix is undone (the online refresh had the composition gap too). 402 tests.
      Not cleared for use (2026-10-31): the browser save would itself work around the 403 (`docs/DATA_POLICY.md`).

- [x] 2026-10-31: emulator calibration (feat-005; `experiments/calibration/README.md`): vec.evaluate
      protocol="bake"; all 66 published validation anchors exact at seed 0; T1 pseudobulk_shift, T2 embryo
      pseudobulk_shift, T3 shift_transfer and gene_ko equal the published rows; the T2 heart-extrapolation
      pseudobulk_shift row is unexplained (14 constructions; cell-type matching across E8.75 and E9.5 is unpublished).
      A2 revised. ./init.sh: 412 passed.
- [x] 2026-10-31: working answers to the organiser questions (`docs/ORGANIZER_QUESTIONS.md`), predicted by seven AI
      agents at the captain's instruction after no reply; screened reads of Slack and the Data Protocol page
      (`docs/DATA_NOTES.md`). T3 selects the allowed branch; no new cross-modal input; the 403 is Cloudflare error
      code 1010 (browser-signature block); the D1 check is ready and never prints test-board counts.

- [x] 2026-11-01, the captain's "yes for all": (1) the adoption rule (PLAN section 3, item 8; `vec.adoption`,
      7 tests) replaces "improves >= 2 of 3 pseudo-tasks", which kept a change with no effect half the time;
      (2) the browser route for the test-board specs is cleared (notice first, save not before 2026-11-04 UTC);
      (3) T3 s2 = the refused branch's s1 (no Gata4-derived number), uploaded first. Commits pushed to GitHub.

- [x] 2026-11-01: pseudo-task harness (feat-008; `vec.harness`, `configs/pseudo_tasks.yaml`,
      `experiments/harness/run.py`, `vec.registry`): reviewed by four independent agents (32 issues fixed in
      57aa7ab), 13 of 13 safeguards mutation-tested, ./init.sh 461 passed; 10-seed baselines on all 10
      pseudo-tasks (`experiments/harness/results/baselines_v1/`). A copy at the count rule beats the floor on the
      embryo board only; on T1 a copy loses 2.0 points. Found: never build a test file at seed 0 (runbook).
      Awaiting the captain: how inapplicable, tuned-on and replicate pseudo-tasks count; the T1 tissue-mix
      geometry; whether the screen may release the qualitative E12.5 collection text.
      (Answered the same day, a second "yes for all": PLAN section 3, item 8.)

- [x] 2026-11-01: D2 check of B0 (feat-009; `experiments/harness/results/d2_v1/`, 10 seeds, clean trees). B0
      beats the floor on the T2 embryo board only (+1.21 and +3.86 points, 10 of 10 seeds). T2 heart: +0.43 on the
      one-day task but ahead in 7 of 10 seeds (8 needed), +1.45 on the 1.75-day task. T1: -0.32 on the chamber-cut
      truth (4 of 10 seeds; the 0.3-point veto), 44.1 on the whole E10.5. B0 beats a plain copy on all three
      boards. On T1 a fresh draw of reference cells starts about 2 points below the floor, all on de_score (by the
      metric's code: its null ranks genes by reference expression, and the floor scores exactly 0); and E9.5's
      chamber mix is far from E10.5's (blood about 2.5% against 24% of chamber cells), so B0's scope step gains
      nothing on mmd_u and loses on the variogram. T3 has no harness check (the pre-registered entries were selected
      or fitted on every fold). The first T1 run stopped when the read guard refused B0's own read of the held-flat
      list; the harness now hands the list over (20fb8b5). Awaiting the captain: the T1 blood share at test, the T3
      D2 reading, the push. (Answered the same day: next entry.)

- [x] 2026-11-01, the captain's third "yes for all": the T1 test entry keeps E10.5's own chamber share of blood
      (about 24% of chamber cells, not the plan's 1.5-2%; PLAN 4.1), and the s2 rule is to be proposed with the
      hedge on the blood share; the pre-registered T3 entries stand on the T3 rules, and new T3 changes face the
      harness on folds they were not tuned on (PLAN 4.4); the five D2 commits pushed (6c214e3..b2fc85c).

- [x] 2026-11-01: harness rows keep each candidate's metric panel (ab072da). T1 de_score audit
      (`experiments/harness/results/t1_de_v1/`, a rerun of the D2 T1 candidates that reproduces d2_v1 seed for seed):
      the copy handicap belongs to the E9.5 -> E10.5 change, not to fresh cells. Its changed genes are mostly highly
      expressed reference genes going down (426 of 558 on the chamber cut); the metric's expression null names 39-41%
      of them, a copy 11% (-2.8 to -3.0 points), B0 20-28% (-1.25 on the chamber cut). On E8.5 -> E9.5 the null names
      2% and a copy beats it.
- [x] 2026-11-01: embryo base: pseudobulk_shift replaces B0 by the adoption rule (+5.76 points, ahead in 10 of 10
      seeds on the E7.5 pseudo-task, the one where it runs; `experiments/harness/results/d2_v1/psb_vs_b0/`); vetoed on
      the heart board (-2.75 on the 1.75-day task). It edits log values, not counts (PLAN principle 1).
- [x] 2026-11-01: M1 direction step (401e3eb): `vec.t1_scope.external_shift` and the directions option of
      `build_submission` (golden outputs unchanged; 2 new tests), `experiments/t1_m1/feng_dirs.py` (the readout's
      primary predictor, manifest-checked). The harness cannot score Feng (its E9.5 cells lie inside every T1
      pseudo-task's window, and the readout chose Feng on the E9.5 -> E10.5 step both counting tasks use). Labelled diagnostic (`experiments/t1_m1/results/diag_v1/`, 4 registry rows): +1.29 over B0
      on the chamber-cut truth and +1.09 on the whole, 10 of 10 seeds each; +0.97 over the floor on the chamber cut,
      10 of 10. Found: E10.5's chamber mix puts cardiomyocytes at 26.5% of heart cells (E9.5 66%), against s1's
      planned 58%. Heart: the composition trend is proposed for removal (PLAN 4.3 step 2). Slip: the code was synced to
      the DGX while the T1 reruns ran (against the compute rule); the runs had loaded their code, and their scores
      reproduce d2_v1 exactly. Awaiting the captain: Feng's direction on the readout; the Feng re-cut from the raw GEO
      object and its licence; the heart composition step; the embryo's log-value edits; the push. (Answered the
      same day: next entry.)

- [x] 2026-11-01, the captain's fourth "yes for all" (PLAN decision log (f)): Feng's direction stands on the
      readout; the Feng re-cut went ahead and the licence is confirmed; the heart's composition trend is dropped; the
      embryo may use log-value edits; the four commits pushed (b2fc85c..8fca1cf). Then "put all the effort to win
      this challenge at all drastic costs", taken within the rules.
- [x] 2026-11-01: Feng clean reduction from GEO (feat-007; `experiments/external/feng2022_geo.py`, commit 6093147,
      DGX): the 5.35 GB Seurat object parsed in 118 s without R (rdata, in its own uv group `external`); 19,073 of
      27,445 cells kept; E9.5 and E14.5 cells and counts identical to the 2026-10-18 subsets, all 393 group totals
      equal to float rounding; raw files deleted (checked). M1 now reads these statistics; its directions match the
      readout's to 6e-7.
- [x] 2026-11-01: FGW embryo pipeline built (bd2e433): `vec.core` coupling step (fgw_coupling, coupling_rotation,
      barycentric, nearest_anchors), `vec.t2_build.interpolate_clouds` and `pchip_rms`, harness candidates `fgw`
      and `psb_size`; 152 tests in the touched files pass. 10-seed runs on both embryo pseudo-tasks started on the DGX
      (runs/harness/emb_fgw_v1).
- [x] 2026-11-01: embryo results (`experiments/harness/results/emb_fgw_v1/`, 10 registry rows): fgw +4.15 over
      pseudobulk_shift and +2.29 over psb_size on E7.5, +4.92 over B0 on E7.25, 10 of 10 seeds each; psb_size +1.86
      over pseudobulk_shift. Kept: psb_size over pseudobulk_shift (1 of 2), fgw over B0 (2 of 2); fgw over psb_size
      waits on the reading of a task the base cannot run on. Weak spots, left alone (any fix chosen now would be tuned
      on both tasks): FGW's shape on E7.5 and its variogram on E7.25.
- [x] 2026-11-01: test builds (`experiments/builds/build_test.py`, 55cbe3d): T1 M1, heart B0 s1/s2 with the
      pre-registered size, embryo pseudobulk_shift or fgw with the pre-registered size; dry runs against the
      validation boards passed their checks on the DGX (seed 1,000,001) and were deleted. Feng's E9.5 -> E14.5
      cardiomyocyte share: 89.5% -> 85.5% of heart cells (77% counting fibroblast-like cells).
- [x] 2026-11-01, the captain's fifth "yes for all": a pseudo-task the base cannot run on is not counted
      (`vec.harness.decide`, "base does not apply"; test with a mutation check), so fgw is kept over psb_size (+2.29)
      and pseudobulk_shift (+4.15) on E7.5 (1 of 2) and over B0 on both: **the embryo entry is fgw**. The s1
      cardiomyocyte share stays at E10.5's own (26.5% of heart cells). Pushed 8fca1cf..6a90a93.
- [x] 2026-11-01, the captain's sixth "yes for all": the four commits be612bf..77ab626 pushed.
- [x] 2026-11-01: independent review of the s1 build path (workflow: four reviewers, fidelity, maths, test-board
      robustness and compliance; each finding checked by one or two skeptics; 23 agents; 11 confirmed, 2 refuted).
      Fixed in 46089fe (./init.sh 504 passed; each critical fix mutation-checked):
      - critical: both board intakes compared the full sha256 of a gene list with the site's genes_sha256, which is
        the first 16 hex digits of sha256 over the newline-joined symbols (all five committed lists verify that way),
        so the 4 Nov intake would have refused every board. `vec.boards._genes_match` takes both forms.
      - critical: builds, `vec check`, `vec submit` and `vec boards` printed a test board's cell limits or a test file's
        count (max_cells under the count rule). `Board.counts_hidden`; the tools say "max_cells" or "hidden".
      - critical: the T3 builder defaulted to seed 0 and needed max_cells typed by hand: now seed 1,000,001 by default,
        refused below the offset, and `--n-cells max` (the default) reads the board in code; runbook updated. The T3
        decision file pins t3_decompose.py and rules.yaml only, so the pre-registered build is unchanged (its 22 local
        tests pass).
      - minor: null min_cells; obs_required/obsm_required loaded, checked and listed by the D1 check; T2 test counts
        through the committed count decision (both overrides are null, so max_cells; provenance recorded); the embryo
        fgw build no longer loads E7.25 (it lacks two panel genes); dry runs need a validation board and "dry" in the
        file name.
      - Compliance note: one verifier's grep printed a single line of docs/research/engineering.md (about the gene
        digest format; no board values, no sizes); its conclusion rests on its own computation over boards/. Reported
        to the captain.
- [x] 2026-11-01: assumption A5 revised (the test as written cannot run; the M1 diagnostic met the de_score bar,
      +0.110, not the de_direction bar, +0.026); feat-010 set to in progress (built and the entry; A7's 70 not met).
- [x] 2026-11-01: s2 hedge payoffs (`experiments/s2_hedge/`, sensitivity analysis, not a choice): T1 blood and
      cardiomyocyte axes on the chamber task (10 seeds, cf476aa; B0 at its own mix reproduces d2_v1 exactly). Tissue
      mix moves T1 by up to 12-19 points; both one-axis runs are confounded by E9.5's own mix, so a 2-D grid on an
      oracle base and a heart scope axis (prompted by a side note: heart s2 = H2 differs from s1 only in size, -9%,
      while E10.5's heart sample is 77% surrounding tissue) were started on the DGX at 46089fe.
- [x] 2026-11-01: no organiser reply on the fourth-email thread; the captain's notice went out 2026-10-31 18:59 UTC.
      Screened Slack sweep (a separate agent; #virtual-embryo-challenge and #general, 2026-10-28 00:00 to 2026-11-01
      13:40 UTC): no organiser post on the board specs, the 403 or browser saves (no objection, no fix), nothing
      addressed to our team, no rule, deadline, nomination or data-use change. One organiser reply (2026-10-29): the
      Data Protocol is on the site and the training files' SHA-256 checksums are on the Data page (A1 could add a
      checksum comparison; the protocol's held-out text stays withheld). Two threads skipped (other teams' standings;
      a team's final-phase method details).
- [x] 2026-11-01: the s2 hedge runs on the DGX stalled. All 27 processes spent 99% of their CPU time in the kernel:
      numpy's transparent-huge-page requests triggered memory compaction on a full, fragmented NUMA node 0, with about
      435,000 failed page migrations a second. In 6.9 hours a heart seed took 4.8 hours and no grid seed finished. The
      runs were stopped at 20:25 UTC and the failed migrations dropped to 0 at once. The logs are kept in
      `runs/s2_hedge/logs/stuck_2026-11-01/`; heart seeds 0 and 5 finished and serve as a check on the rerun.
      Relaunched with `NUMPY_MADVISE_HUGEPAGE=0`, and the grid now holds one prediction at a time (it held about
      20 GB per process; same scores, since every call is seeded). The fix and how to spot the stall are in
      `docs/COMPUTE.md`.

### What's Next

1. Until 2026-10-27: work that needs no release data (`docs/PLAN.md` section 7): external data through committed
   window filters (feat-007), v1 pipelines on released stages, B0 builds checked with `vec check` on the DGX.
   (2026-10-28: that week is over; these items are still open in `feature_list.json`.)
2. 2026-10-27, once `vec submit status` shows the final phase: follow `docs/P3_RUNBOOK.md` (board intake with anchors
   stripped, the quota shape, data, label map, pre-registered runs, commit each output), then emulator calibration
   (feat-005) and test-board intake (D1).
   (2026-10-28: step 0 passed on 2026-10-27, step 1 is blocked by the 403 (fallback built, not cleared for use), and
   steps 2, 3 and 4 are done (What's Done). Next is the runbook's section 5: emulator calibration against the
   validation reference rows (feat-005), then candidate files on the DGX; the test-board intake (D1) and every
   server format check wait for the test boards' specs.)
   (2026-10-31: feat-005 done; the D1 check script is ready (`experiments/intake/d1_check.py`); D1, candidate files
   and format checks wait for the captain's decision on the browser route.)
3. B0 floor-plus files on all four boards with the builders (feat-009); T3 per-gene factor on Gata4 (feat-012); T1
   scope shares by a written rule (feat-014).

## Blockers / Risks

- [ ] Awaiting organiser replies: leak report (Slack DM, sent 2026-10-18), the rules questions (email, 2026-10-18) and
      the follow-up (scorer build and parity defect, wrong-sign severity, released Gata4 in the T3 prediction,
      duplicated/interpolated records, z units, truth_cells and size, stage conventions; sent 2026-10-19). Screened
      sweeps of Slack and the organiser email on 2026-10-20 (to 22:43 IST): no organiser post on Slack since
      2026-10-16, no reply in any DM or by email (the sender's mail tracker reported both emails unopened 24 hours
      after sending); the platform emails are the postponement (07:41 IST) and the collaborative Yukon track (22:03
      IST), which gave the deadline's time zone (through 2 December 2026, UTC). Unanswered in public: qualitative
      in-window priors and (2026-10-19) whether released validation answers may be used for test predictions (A20).
      Final uploads still finish by the end of 30 November IST; keep the strict in-window reading. The third email
      (three questions: the test phase's start time, validation answers as training data (A20), the Generality scope
      (A21)) was sent by the user on 2026-10-21 at 01:24 IST (`docs/ORGANIZER_QUESTIONS.md`). Not checked by agents:
      the challenge account messages page (user only) and Gmail Spam.
      Screened Slack sweep on 2026-10-23 (from 2026-10-20 22:43 IST): still no organiser post or reply anywhere;
      two new unanswered participant questions (2026-10-21: textbook E9.5-E12.5 knowledge as a prior and whether to
      disclose it; how travel grants are chosen).
      The fourth email (the 403 on the panel files, what counts when the captain marks no final entry, confirmation
      of the A20 reading; the scorer build still open) was sent by the user on 2026-10-28 at 15:41 IST, signed "Team
      captain, Mithril" (`docs/ORGANIZER_QUESTIONS.md`). No reply to any of the four emails in that mailbox as of then.
- [ ] **Main blocker for submitting (2026-10-28):** the board specs. `panels/index.json` answers HTTP 403 to our
      client (not worked around), `boards/` holds only the P2 specs, and `vec submit` checks the board against
      `boards/index.json` first, so it refuses every upload to the four test boards, format checks included, until
      the specs are in (`docs/P3_RUNBOOK.md`, step 1). Step 4 does not need them; candidate files need the test
      boards' cell limits. Asked in the fourth email (sent 2026-10-28). The options are the user's call: wait for the
      reply, ask in #virtual-embryo-challenge, or lift the no-workaround rule so the user can save the file in the
      browser for the intake that strips the anchors in memory (built 2026-10-28, cee1c92, not cleared for use;
      `docs/P3_RUNBOOK.md` step 1). 2026-10-31: the online refresh still answered 403 (16:52 and 18:20 UTC; Cloudflare
      error code 1010, a browser-signature block), and none of the four organiser emails has a reply (Mailsuite
      listed the fourth as unopened on 2026-10-29). The predicted organiser answer is no objection to a browser save
      (0.8); the route waits for the captain's yes, then a dated notice and about 72 hours (runbook step 1).
- [ ] Heads-up from step 4.1 for the embryo build (no rule acts on these by itself): `docs/PLAN.md` asks to check E7.5
      against E8.0 by proper vs improper shape alignment before building the embryo board, and that check came out
      undecided (relative margin -0.019 +/- 0.023 over 5 seeds), so proper rotations stay. E8.5 against E8.75 reads
      improper although both heart stages resolve right-handed by landmarks. The z steps differ (E7.5 14, E8.0 15),
      so the z-convention risk goes on record with the embryo submission (`docs/PLAN.md` section 4.2, step 5).
- [ ] Release day (2026-10-27, checked to 17:43 UTC): the final phase is open (public phase endpoint: p3, split
      test, accepts submissions; official quota keys T1:test, T2:heart:test_extrap, T2:embryo:test_interp and
      T3:bcatenin, 0/2 used each; daily per task: scored 0/2, uploads 0/8). Runbook step 0 passed (357 tests, nothing
      unresolved). Three things stop the intake: (1) the data manifest was refused (HTTP 409) until the user accepted
      the current Challenge Rules and Terms of Use on the account page (done about 17:42 UTC; an API key cannot
      accept them); (2) it then lists only the 9 training files, with the note "Training inputs only. Validation and
      test targets are not distributed - a submission comes back as a score." (first recorded today; unknown whether
      it is new or a leftover from the validation phase, and it conflicts with the Data page's promise that the
      validation answers become training data); (3) `vec boards update` fails: panels/index.json answers HTTP 403
      from Cloudflare to our client, so nothing was written (not worked around). The screened read of the new Rules
      and Terms was stopped by the user before it reported: what changed in them is not yet known. Nothing was
      downloaded, built or uploaded. The user asked to try again on 2026-10-28.
      Retry 2026-10-28 04:37 UTC: phase, quota and manifest unchanged (9 training files, same note); `vec boards
      update` still 403, nothing written. A screened read of the site (04:40-04:49 UTC; quotes came through WebFetch,
      so check the live wording before relying on them) explains the missing files: the validation answers are now
      opt-in. Data page and Rules s.10: they are "available to eligible registered participants during P3, after
      explicit confirmation", and the opt-in files are "a specific exception to the held-out-data restrictions ...
      solely for preparing P3 predictions after access is granted" (no external E10.5/E12.5 data; no test truth).
      Rules s.11: obtaining access (even without downloading) or the first accepted P3 submission "permanently ends
      P2 submissions across all tasks for the entire ordinary team"; viewing the Data page or a format check does not.
      No test-input files will come ("Hidden test ground truth is never distributed"). P3 opened 2026-10-27 16:00
      UTC (answers question 1 of the third email). The Data page still links the panel files under /challenge/panels/
      and announces no outage, so the 403 is unexplained. No page says how the opt-in works, who is eligible, or
      whether a Human Team must nominate entries (the phase endpoint shows nominations_open). Waiting on the user:
      the opt-in (irreversible, team-wide, browser only) and the reading of "solely for preparing P3 predictions"
      for inputs at prediction time (A20 stays open until then).
      2026-10-28 09:32 UTC: the user unlocked the validation answers on the Data page (server state: granted
      09:32:18 UTC, lock_scope team, so P2 submissions are permanently disabled for the whole team; the page says it
      uses no P3 attempt). They have their own list and links (`/challenge/data/validation/` state, manifest and
      link, found in the site's public script); the unlock is a separate POST that no repo code sends
      (`tests/test_data.py` checks for it). The data client reads the validation list since e0e1d58, and all 6 files
      (2.09 GB) are on the DGX: runbook step 2 done (checksums, inventory and QC; `docs/DATA_NOTES.md`). The Data
      page's own words on use: "P2 validation answers can be used to prepare your final P3 models" and "After you
      explicitly unlock the official P2 validation release above, you may use those released validation files for
      P3 preparation. This does not permit access to hidden test data or other restricted external sources."
      Still open: the 403 on the panel files (step 1, the test-board specs) and the user's reading for A20 (inputs
      at prediction time), which steps 4.1-4.3 need. (Later on 2026-10-28 the user read Rules s.10 and the Data page
      as allowing it: A20 passed, `ASSUMPTIONS.md`; the 403 remains.)
- [x] 2026-10-28: runbook step 3 (labels) done, at the user's go-ahead. The 30 labels first seen in the validation
      answers are mapped in `configs/lineage_map.yaml` by its 2026-10-19 rule, and the 8 labels new to a T1 file in
      `configs/t1_chamber_classes.yaml` (6 new names plus OFT-CM and Gut Endoderm). Evidence: classifier calls
      against the nearest released stage (`experiments/labels/`, results committed; the T1 rule's own script in
      `experiments/t1_scope/results/label_evidence_E9.5_E10.5/`) and marker scores (`scripts/label_markers.py`, DGX).
      The AI assistant mapped them with three independent mappers (29 of 30 unanimous) and a judge (BW).
      - New core labels: Endocardium-1, Endocardium-2 and EndoMT Endocardium (T1), Valve-EndoMT-Endo (T2 E10.5, 2,621
        cells) and the Gata4 files' CM.
      - Calls a reader may question: Early VIC (lineage Unknown, so shell; T1 category endocardium, as cushion
        mesenchyme by rule 2); BW (pharyngeal mesoderm, by the classifier); NCC-2 (NCC by name, though 94% of
        replicate 1 is called surface ectoderm); Pericardium/Body Wall (pericardium by name; the calls are mostly
        JCF); Endocardium/Intra-Embryonic Endothelium (Endothelium, shell, 57-64% of calls).
      - Every label of all 15 files maps and the T1 config checks clean (DGX); 366 tests pass.
      - Names already in the map keep their mapping, as the rule says, though some mean other cells in the new
        files: T1 E10.5 BEC is 57% called E9.5 Endocardium (both chamber categories), T1 E10.5 Surface Ectoderm 80%
        Foregut, and T2 E10.5 Branch Arch 50% Dorsal-PAM (all shell or surrounding either way).
- [ ] Final-entry nomination (found 2026-10-28 in the site's public account-page script): "In the final phase the
      captain picks which submission counts" (a "Final" column in the captain's submission list while
      nominations_open); the timeline makes 2026-12-02 the last moment to nominate. Rules s.11 still says the board
      ranks a team by its best; what counts if the captain marks nothing is not stated. The user (captain) marks
      the final entries after the last official upload; ask the organisers what an unmarked board gets
      (`docs/P3_RUNBOOK.md`, section 5). Asked in the fourth email (sent 2026-10-28); no reply yet.
- [x] 2026-10-27: history of `main` rewritten (user's decision) to drop the AI co-author trailer from all 44 commit
      messages; trees, parents, authors and dates checked identical; old -> new IDs in `docs/COMMIT_MAP.md`;
      force-pushed, DGX .git refreshed. Older IDs quoted in files resolve through that table.
- [ ] The user's decision: the collaborative Yukon track (organisers' email of 2026-10-20 22:03 IST). The entry stays
      an official Human Team unless the user decides otherwise; nothing has been done toward joining. Moving is
      one-way and ends official-prize eligibility, a team left empty is archived, and the collaborative phase was
      to end before P3 opened on 2026-10-27. (2026-10-28: P3 is open and the entry is still an official Human
      Team, Mithril.)
- [x] 2026-10-21: the user signed off ("approve all") on the other agent-written rules that act on the
      release (heart and embryo size procedure, T1 direction readout, cell-count rule, label-mapping rule): all eight
      recommended answers in `docs/RELEASE_RULES_SUMMARY.md` (drafted by the AI assistant on 2026-10-20). The gaps it
      found are closed:
      - the T1 readout's second role is written down (`experiments/t1_direction/README.md`, PLAN 4.1);
      - poor harmonisation groups are held flat (PLAN 4.3);
      - the fallback for an invalid heart size and the count rule are now in code (`vec.t2_build` and
        `summarise.py --decide`; see What's Done).
- [x] 2026-10-21: the user's call on the case approved answer 7 left open, a winning count level below a test
      board's lower limit (truth_cells or min_cells; only 0.25 x max_cells can fall there): the user approved ("go
      ahead") the recommendation "the next passing count the board allows, or max_cells if none does". In code as
      `vec.t2_build.COUNT_FLOOR_RULE` (see What's Done).
- [x] 2026-10-21: how the code reads "next" in that rule. The user left it to the AI assistant ("Anything that would
      improve my chances of winning"). A decision panel chose the next count up over the next best gain: they differ
      only when 0.25 x max_cells wins below the limit while 0.5 x and 0.75 x both pass and 0.75 x measures at least as
      high, and there the next count up is ahead by about 0.2-0.4 points in most simulated models (a data analysis put
      it near zero; `experiments/density/floor_reading/README.md`). In code (see What's Done).
- [ ] The E12.5 emergent-class list waits on the organisers' answer about qualitative in-window background; the legal
      route is used meanwhile.
- [ ] Quota shape: `official_per_board` is unpublished. The client refuses official uploads until it can read it;
      check with `vec submit status` once the final phase opens and adapt the parser (with a review) if needed.
      2026-10-28 (as served since 2026-10-27): the final phase serves it (keys T1:test, T2:heart:test_extrap,
      T2:embryo:test_interp and T3:bcatenin, 0/2 used each), but the client reads a board's entry only for a board
      in `boards/index.json`, which still holds the P2 specs (step 1, 403), so the test boards' entries are not yet
      checked through it.
- [x] Upload logistics (decided by the user, 2026-10-19; built and verified the same day): `vec submit PATH --board KEY
      --from-dgx` pulls the built file into a private copy for that call only; the key never goes to the DGX.
- [x] Provenance across machines: "dgx provenance" alone checked only the code checked out at upload time, not the
      code that built the file (a side agent caught this; I had marked it closed too early). Fixed 2026-10-20 with
      build stamps: every prediction file gets `<file>.stamp.json` when it is built (sha256, the commit, whether the
      tree was clean from the build's start to the write), and an official upload needs a stamp that matches the
      bytes sent, says clean and names a commit in our history; the ledger's build_commit records it.
- [ ] T1 pseudo-task gain from exact category shares is +31 against about +24 in the screened digest; both are upper
      bounds (exact shares), the real gain depends on how well the E12.5 shares are chosen.

## Decisions Made

- **Human Team track only** (user, 2026-10-18).
- **Python 3.12; torch 2.14 cu130; veckit pinned to git main cc3859a.**
- **No mouse single-cell foundation models** (in-window contamination). **No gene-level T3 knockout literature.**
- **Heart size only from the pre-registered training-stage procedure** (docs/PLAN.md section 4.3).
- **Count = max_cells, never below truth_cells; random rows only** (A13 revised, 2026-10-19).
- **Heart size procedure fixed before the release** (core = CM + endocardium; E9.5-E10.5 only; fixed three-rule average;
  z held when the plane matches; H1/H2 scope). **T3 leads with the per-gene factor** from released knockouts.
- **No parity hedge without the organisers' approval** (user, 2026-10-19). **History rewritten again** (user).
- **2026-10-27 and 2026-10-28** (user): history rewritten once more to drop the AI co-author trailer; the validation
  answers unlocked (P2 submissions end for the whole team) and read as allowed as training data and model inputs
  (A20; organiser confirmation asked); the go-ahead for runbook step 3.
- **GPUs:** free VRAM only until the end of October; all eight H100s from November (user, 2026-10-19).
- **Uploads pull from the DGX** (user, 2026-10-19): `vec submit PATH --board KEY --from-dgx` streams the built file
  from the DGX into the client's private temporary copy at upload time; no lasting local copy, and the key never
  leaves this machine.

## Evidence of Completion

- [x] `./init.sh` on 2026-10-28: lock check, ruff clean, 402 tests pass. Recorded on 2026-10-19 (237 tests then):
      experiment tests 203 pass, 2 skipped locally (the T3 end-to-end tests, which pass on the DGX: 50 T3 tests in all).
- [x] DGX: torch 2.14.1+cu130 with CUDA, jax 0.11.2, moscot 0.5.2, scvi-tools 1.5.1 import (runs/logs/setup2.log).

## Notes for Next Session

- Slack sweeps: leave the group DM C0CA7QTEMDF out of date searches that reach back to 2026-10-18 (they surface our
  leak report). Six workspace accounts use eigenlabs.org addresses (none has posted); record their posts, or posts in a
  new Yukon channel, by place and date only.
- The Data Protocol (label table, units, sectioning) was promised before the test phase opens (an organiser,
  2026-10-16): look for it on the Data page this week (never the board file). (2026-10-28: the test phase opened on
  2026-10-27; whether the Data page now has the protocol is not recorded.)

- `neighborhood_mmd` builds 15-NN neighbourhoods over all submitted cells, so submitted density moves T2 scores.
- Ctnnb1 is not on the 500-gene MERFISH panel: gene_ko gives nothing for the beta-catenin test board.
- T1 E12.5 is a heart-chamber-only dissection; earlier stages include surrounding tissue.
- DGX research scratch removed on 2026-10-19 (user's OK); `~/v/vec-gapconc` keeps only `ext/` and `out/`, the
  legal inputs the committed filters in `experiments/external/` read (the JAX files are hard-linked into data/external).
- DGX: cap the BLAS threads for numpy and scikit-learn jobs (`OPENBLAS_NUM_THREADS` 64 or less; the runners set
  `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS` and `MKL_NUM_THREADS` themselves, e.g. 8 in `experiments/coords/run_p3.sh`
  and 16 in `experiments/t3_decomp/run_gata4.sh`). Without the cap a run crashed on 2026-10-28 ("OpenBLAS warning:
  precompiled NUM_THREADS exceeded", then "malloc(): corrupted top size").
- `experiments/harmonise_review/review.py` has 5 pre-existing ruff findings (outside the checked `src tests`).
- Release day (2026-10-27): re-check #virtual-embryo-challenge and the Rules/Data pages (never the board file) for
  the release time, submission rules, scorer version and any Yukon details the organisers post; an organiser said on
  2026-09-01 that T2 E8.5 is released when the final stage starts. Before the first upload, confirm that our account
  sign-in works (a participant saw "Sign-in is not available yet" on the account page on 2026-10-18).
  2026-10-28: P3 opened on 2026-10-27 at 16:00 UTC (Rules s.11); the Rules and Data pages were read, screened, on
  2026-10-28 (Blockers above); T2 E8.5 came as an opt-in validation answer and is on the DGX. Not re-checked yet:
  Slack (last swept 2026-10-23); the scorer build is runbook step 4.6. Sign-in works: the user signed in to accept
  the Rules (2026-10-27) and to unlock the validation answers (2026-10-28).
