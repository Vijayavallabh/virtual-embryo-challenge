# Session Handoff

## Current Objective

- Goal: a compliant Human-track entry on all four test boards that finishes as high as it can. The final phase runs
  through 2026-12-02 (UTC), with two scored uploads per board for the whole phase; none is used as of 2026-11-13.
- Current status (2026-11-13): the entries are built and their methods chosen. T1: M1 with R1, R2 and R3b, its first
  upload drawn from E10.5's dissected-heart library since 2026-11-13 and its second on the pooled base at 60%
  cardiomyocytes. Heart: B0 with the pre-registered size. Embryo: fgw, with `fgw_aligned_series` as second upload
  since 2026-11-13. T3: the frozen rules' files, the second upload now with a wild-type-only factor. Nothing is
  uploaded. The captain's upload gate holds every first upload until the team's own estimate passes 209 (it is 184,
  range 170 to 198) or the deadline exception opens at 2026-11-27 00:00 UTC; each upload then needs the captain's
  go-ahead at the time, and the captain runs the commands.
- Nothing waits for a ruling (`docs/PLAN.md` decision log 2026-11-13 (d); `docs/S2_RULES.md`, "Amendments of
  2026-11-13"). The ten files as amended were rehearsed on the DGX on 2026-11-13 (decision log (e)): two byte-identical
  sets, and the freeze-morning reference `runs/archive/builds_fd2683a/SHA256SUMS` there
  (`experiments/builds/results/rehearsal_fd2683a.txt`). What is owed:
  - on 19 and 23 Nov, the approved screened reads for an organiser answer to the seventh email
    (`docs/SCREENED_READ_SEVENTH_EMAIL.md`; a reply in the captain's inbox is not pasted into the session);
  - about 19 Nov, the Slack follow-up draft on floored T3 severity if no organiser has answered;
  - a new rehearsal and reference by 24 Nov only if a build changes after 2026-11-13.
- Branch / commit: `main` on GitHub (private). The DGX checkout (`PrakashDGX_H2:~/v/virtual-embryo-challenge`) is at
  fd2683a, synced on 2026-11-13 for the rehearsal; sync it again only when no build runs there.

## Completed This Session (2026-11-12 to 2026-11-13: the campaign the captain asked for on 2026-11-12)

- [x] First and second wave (2026-11-12; decision log 2026-11-12 (a) to (c)). The finding: released samples differ
      more by how they were measured than by stage (diagnostic N23; the T1 library census in `docs/DATA_NOTES.md`).
      Embryo tries E19 to E25 (depth and run renderings, pair pricing), T1 tries T27 to T30 (library kind and depth
      at zero time gap), T3 tries K6, a second research pass (`experiments/research_2026_11_runs/REPORT2.md`), and
      requests R-A, R-B and R-C.
- [x] Third wave (2026-11-13; decision log (a)). An independent recomputation of the record: every number held, the
      presentation was corrected. Tries E26 to E28 and T31, diagnostic N24. R-A restated with three options, priced
      on nine practice tasks.
- [x] Fourth wave (2026-11-13; decision log (b) and (c)). The R-C email, drafted after a three-lens review and sent
      by the captain. Reading N25 and T3 tries K7 and K8: a per-gene factor from wild-type samples alone, closed by
      its pre-set reading, with both readings given to the captain. Both R-A candidate files reproduce at 16 BLAS
      threads (`experiments/builds/results/determinism_1f8956b.txt`), and R-A's second option has a build command
      (`build_test.py embryo --method fgw_shape0 --expression render_hi`).
- [x] The captain's rulings of 2026-11-13 (decision log (d)): the embryo's second upload, T3's second upload and T1's
      first-upload base changed (rows R4, R5, R6), with `docs/S2_RULES.md` amended, two build options, two independent
      code reviews, and dry builds with one sha256 per new file. No upload. Tries: T1 28, embryo 36, heart 4, T3 15.
- [x] The DGX rehearsal of 2026-11-13 (decision log (e)), at the captain's instruction: the checkout synced to
      fd2683a, all ten files built twice from empty folders (byte-identical, clean stamps, `vec check` passed, seven
      equal to the archive of f337dc6), the set archived as the freeze-morning reference, and the freeze-morning
      command tried against it. T1's new s1 and the embryo's new s2 have the same sha256 on both machines; T3's
      refused s1 with the factor does not (the draws of 84 cells; a nearest-neighbour step follows the BLAS kernel),
      so only DGX hashes check DGX builds. No server format check, no upload.
- [x] A copy off the DGX (2026-11-13, after a side agent's note the captain forwarded): this machine holds a checked
      copy of the reference set (`runs/archive/builds_fd2683a/`) and of the T3 run folder
      (`runs/t3_decomp/gata4_20261028T122425Z/`), owner-only. The client's build-stamp guard and local format check
      pass on all ten copies (nothing sent), so an upload can go from the copy if the DGX is down, at the captain's
      word (`docs/P3_RUNBOOK.md`, "A copy off the DGX").

## Completed Earlier (2026-10-27 to 2026-10-28)

2026-11-13: kept as history; the sessions of 2026-10-29 to 2026-11-11 are logged in `progress.md`.

- [x] Release day (2026-10-27): the final phase opened and runbook step 0 passed; the user accepted the current Rules
      and Terms on the account page. Step 1 is blocked: `vec boards update` gets HTTP 403 on panels/index.json (not
      worked around), so `boards/` still holds the P2 specs.
- [x] History of `main` rewritten (the user's decision, 2026-10-27) to drop the AI co-author trailer; old -> new IDs in
      `docs/COMMIT_MAP.md`.
- [x] Validation answers: the user unlocked them on the Data page (2026-10-28 09:32 UTC; P2 submissions permanently
      disabled for the team); the data client reads their list and links (e0e1d58), and all six are on the DGX with
      checksums, inventory and QC (runbook step 2, f423f6a; feat-003 done). No test inputs are distributed.
- [x] Labels (runbook step 3, at the user's go-ahead; d6b9566): 30 labels in `configs/lineage_map.yaml` and 8 in
      `configs/t1_chamber_classes.yaml`, mapped by the AI assistant (three independent mappers, a judge for BW);
      evidence in `experiments/labels/` (aebcb8c) and `experiments/t1_scope/results/label_evidence_E9.5_E10.5/`.
- [x] A1 passed, and A20 passed on the user's reading (d47badf). The runbook now says the captain marks the final
      entries (ea4ae36). The fourth organiser email, drafted by the AI assistant, was sent by the user, signed
      "Team captain, Mithril" (ed3b1bc).
- [x] Saved-file board intake (cee1c92), at the user's request: `vec boards update --from-file` strips test anchors in
      memory like the online refresh; reviewed by three independent AI agents, every finding fixed with a test that
      fails when the fix is undone; 402 tests pass. Runbook step 1 gives the user's steps. Not cleared for use
      (2026-10-31): the browser save is itself a workaround for the 403, so it waits for the user to lift that rule.
- [x] Runbook step 4, in order on the DGX at 725d19c (`progress.md`, What's Done, has the numbers): 4.1 size (both
      JSONs valid; f5be247), 4.2 Gata4 (replicates independent; decision file, aggregates and plan: 95d8198, 7e7c207,
      3ff6bf0, 05ead40), 4.3 T1 readout (build M2; 3c5cc66), 4.4 harmonisation (0586bf6), 4.5 density (max_cells
      kept on both T2 test boards; dfee546, ef1819f), 4.6 anchor check (bake protocol reproduced, no mirror-mode
      Dice anchor; 12151e0).

## Completed Earlier (2026-10-18 to 2026-10-19)

2026-10-28: kept as history; the sessions of 2026-10-20 to 2026-10-23 are logged in `progress.md`.

- [x] Harness, uv environment, `vec` package (emulator, checks, baselines, safe downloader, loaders, QC, harmonise).
- [x] Training data on the DGX; research sweep, screened digest, PLAN v1.2, assumption ledger (19 rows).
- [x] Experiments with independent reviews: density, coordinates/size procedure, harmonisation, data QC.
- [x] Compliance: leaked heart size quarantined and reported; gene-level knockout findings and held-out-stage values
      redacted by meaning; leaked-scorer copies and in-window external files deleted; history rewritten and the GitHub
      repo recreated (first pass; a second pass is ready for the user to run).

## Verification Evidence

| Check | Command | Result | Notes |
|---|---|---|---|
| Lint and tests (2026-11-13) | `uv lock --check`; `uv run --frozen --no-sync ruff check src tests`; `uv run --frozen --no-sync pytest -q` | 1198 passed, 1 skipped (torch not installed) | local machine, started at 1435c22 on a clean tree (`runs/logs/checks_local_2026-11-13b.log`); the tests added later ran on their own; `./init.sh` itself not run |
| Reproducibility of R-A's candidate files | `experiments/builds/build_test.py embryo ...` at 16 BLAS threads, then `vec check` | one sha256 per file over four and three builds | `experiments/builds/results/determinism_1f8956b.txt`; build logs not read |
| DGX rehearsal of the amended entries (2026-11-13) | two chains of `BUILD_OUT=runs/rehearsal_fd2683a/<a\|b> bash experiments/builds/build_all.sh s1`, `s2`, `t3alt`; then `bash experiments/builds/determinism_check.sh runs/archive/builds_fd2683a/SHA256SUMS 16` | ten of ten files byte-identical across the two chains and the third build; seven equal to the archive of f337dc6 | `experiments/builds/results/rehearsal_fd2683a.txt`; exit codes, hashes and the stamps' commit and clean fields only |
| Dry builds of the amended entries | the commands of `docs/S2_RULES.md` (last section) at 16 BLAS threads, then `vec check` | one sha256 over three builds of each new file; unchanged files reproduce the DGX archive, T1's old s1 excepted | `experiments/builds/results/rulings_9cf5ca5.txt`; build logs not read |
| T3 tries K7 and K8 | `experiments/t3_decomp/k7_wt_revert.py` and `k8_s2_form.py` (`run` per task, then `read`) | `experiments/t3_decomp/results/k7_wt_revert/`, `k8_s2_form/` | 10 seeds, clean trees at 33204be and 9dc6af1 |
| Lint and tests | `./init.sh` | passes (2026-10-28) | `uv run pytest -q` 402 tests (90+ on 2026-10-19); experiment tests in `experiments/*/test_*.py` run separately |
| Pre-registered runs (step 4) | `docs/P3_RUNBOOK.md` step 4 commands (DGX, 2026-10-28) | all six ran, exit 0 | outputs committed before reading: f5be247 to 12151e0; the commits are listed in the runbook's step 4 status |
| Size procedure | `uv run pytest experiments/coords -q -p no:cacheprovider` | 33 passed (2026-10-19) | preflight requires committed procedure files |
| QC of released files | `bash experiments/qc/run_qc.sh <tag>` (DGX) | `experiments/qc/results/` | duplicate and interpolated counts in `docs/DATA_NOTES.md` |
| Validation answers (step 2) | `uv run vec data inventory`, `uv run vec data qc` (DGX, 2026-10-28) | 6 files, 2.09 GB, sizes as listed | SHA-256 of all 15 files in `data/raw/SHA256SUMS`; no duplicated or interpolated records; `docs/DATA_NOTES.md` (f423f6a) |
| Labels (step 3) | `uv run python experiments/labels/release_label_evidence.py` (DGX) | label evidence in `experiments/labels/results/` | every label of the 15 files maps (DGX check, 2026-10-28); the new entries pinned in `tests/test_prep.py` and `tests/test_t1_scope.py` (d6b9566) |

## Decisions Made

- Current: `docs/PLAN.md` decision log through 2026-11-13 (e), `docs/DATA_POLICY.md` change notes through
  2026-11-13 (e), the second-upload rules in `docs/S2_RULES.md` as amended on 2026-11-13, and the captain's upload gate in
  `results/upload_gate.json`. The campaign's method decisions are the AI assistant's, logged; every change to an
  upload is the captain's ruling.
- See `docs/PLAN.md` decision log (2026-10-18, 2026-10-19 a-d, 2026-10-20 a-c, 2026-10-21 a-c, 2026-10-28) and
  `docs/DATA_POLICY.md` records. The user's calls of 2026-10-27 and 2026-10-28 (the history rewrite, the unlock, the
  A20 reading, the go-ahead for step 3) are recorded in `progress.md` (Decisions Made, Blockers / Risks).

## Blockers / Risks

- 2026-11-13: the estimate (184) is about 25 under the gate, and no change was adopted in the campaign. What can
  still move the outcome is how the hidden samples were measured (the bets behind R-A, R-B and T3's first-choice
  file), which no released file shows.
- 2026-11-13: no organiser answer is recorded to any of the seven emails, the Slack DM or the public Slack
  question. The next screened check is on 26-27 Nov; one follow-up in the Slack thread may be drafted about 19 Nov.
- 2026-11-13: a ruling that changes a build is needed by about 24 Nov (a T3 builder change takes about a day, then
  a dry build and the rebuild rehearsal). A rendered embryo file reproduces only at the build's 16 BLAS threads.
- 2026-11-13: only DGX hashes check DGX builds. T3's refused s1 with the factor has other bytes on the local machine
  (as T1's s1 had until that day), so a change of the DGX's numerical libraries or processor before 27 Nov would show
  as a mismatch on the freeze morning; a file that differs is not uploaded until the captain has the explanation.
- 2026-11-13: the DGX's default folders `runs/builds` and `runs/t3_build` hold older files with clean stamps (T1's
  two date from 2026-11-04). Uploads are pulled only from the folder the freeze-morning command builds
  (`docs/P3_RUNBOOK.md`, "The freeze-morning rebuild").
- 2026-11-13: the account's usage limit stopped an agent workflow that morning; single reviewer agents ran later
  the same day. Larger agent workflows wait for the captain's word.
- Earlier items, kept as history:
- Organiser answers pending (leak report; rules email). Follow-up email sent 2026-10-19; a third (three questions)
  sent 2026-10-21 01:24 IST; a fourth (403, final-entry nomination, A20) sent 2026-10-28 15:41 IST
  (`docs/ORGANIZER_QUESTIONS.md`); no reply to any of the four as of 2026-10-28. Slack last swept 2026-10-23.
  Deadline: through 2026-12-02, UTC; the captain marks the final entries by then.
- Collaborative Yukon track (organisers' email of 2026-10-20 22:03 IST): joining is the user's decision alone; the
  entry stays an official Human Team meanwhile.
- Second history rewrite done 2026-10-19 (GitHub repo recreated; old-to-new commit map in `docs/DATA_POLICY.md`). A
  third, on 2026-10-27 (the user's decision), dropped the AI co-author trailer; its map is `docs/COMMIT_MAP.md`.
- Never open the public board file, veckit 0.1.2, PR #6, raw research outputs or Yukon team material (AGENTS.md,
  quarantined material).

## Next Session Startup

1. Read `AGENTS.md`.
2. Read `feature_list.json` and `progress.md` (Current State, What's Next).
3. Review this handoff, `docs/P3_RUNBOOK.md`, and `docs/PLAN.md` decision log 2026-11-13 (a) to (e).
4. Run `./init.sh` before editing. On the DGX do not run it as it is (`docs/COMPUTE.md`, `progress.md` notes).
5. Before choosing or changing work on an entry, read the guiding principle in `AGENTS.md` and the headroom tables in
   `experiments/harness/TRIES.md`.

## Recommended Next Step

- 2026-11-13 (latest): the rulings are given and built in, and the ten files are rehearsed on the DGX. Next, in
  order: the screened reads on 19 and 23 Nov and the Slack follow-up draft about 19 Nov; the method summary and the
  disclosure email (with the hashes of `runs/archive/builds_fd2683a/SHA256SUMS`) on 24-30 Nov; the screened check on
  26-27 Nov; on 27 Nov the DGX sync to the freeze commit, the freeze-morning command (`determinism_check.sh
  runs/archive/builds_fd2683a/SHA256SUMS 16`), the server format checks and the first round, each upload on the
  captain's go-ahead and pulled from the folder that command built; the second round on 29-30 Nov; the final marks on
  1 Dec. Before any reading that compares a predictor with a fold target, write its ledger row (row N25 records the
  one time this was skipped).

Earlier notes, kept as history (superseded where `progress.md` says otherwise):

- 2026-11-04: follow `progress.md`, What's Next: the T1 births decision (by 2026-11-10 at the latest), the
  T3 response-identity line through the adoption rule, and the headroom tables kept current. Official uploads only on
  the captain's go-ahead. Every change to an entry starts from its board's published baseline and a row in
  `experiments/harness/TRIES.md`, registered before its run (rules 5-8).

- 2026-11-01: the s1 build path was reviewed (four reviewers, skeptic-verified; 11 findings) and fixed in
  46089fe. The critical one: both board intakes would have refused every gene list on 4 Nov (the site publishes a
  16-hex digest of the newline-joined list); also test-board counts are no longer printed, and the T3 builder defaults
  to seed 1,000,001 and reads max_cells from the board. On or after 2026-11-04 UTC: screened Slack and mail check for
  an organiser objection, the captain saves the board index, then `uv run vec boards update --from-file
  ~/vec-panels/index.json --delete-source`, commit boards/, `uv run python experiments/intake/d1_check.py`, the s1
  builds (`experiments/builds/build_test.py t1|heart --entry s1|embryo --method fgw`, T3 per the runbook), `vec
  check` and format-only `vec submit ... --from-dgx`. s2: payoff analyses in `experiments/s2_hedge/` (T1 grid and
  heart scope ran at 46089fe) feed the s2 pre-registration (17-23 Nov), including the side note that heart s2 = H2
  differs from s1 only in size.

- Board specs (cleared by the captain on 2026-11-01): he sends the notice drafted in his Gmail (a reply on the
  fourth-email thread); not before 2026-11-04 UTC, if no objection or fix arrived, he checks the Data page link by
  hovering, saves it with "Save link as" into ~/vec-panels/ (mkdir -m 700), and the AI assistant runs `uv run vec
  boards update --from-file ~/vec-panels/index.json --delete-source`, commits `boards/`, and runs
  `experiments/intake/d1_check.py` (D1, A3). Test-board counts are read by code only, never quoted.
- T3 (approved 2026-11-01): s1 = the allowed branch's s1, s2 = the refused branch's s1, s2 uploaded first. Builds need
  the test board's max_cells (after the intake).
- Model changes are adopted only under the adoption rule (`vec.adoption.keep_change`; PLAN section 3, item 8),
  through the pseudo-task harness (feat-008 done 2026-11-01: `vec.harness`, `experiments/harness/run.py`;
  baselines in `experiments/harness/results/baselines_v1/`). D2 ran on 2026-11-01
  (`experiments/harness/results/d2_v1/`): B0 beats the floor on the T2 embryo board only, and a plain copy on T1,
  T2 heart and T2 embryo (`progress.md`). The captain answered the same day (a third "yes for all"): the T1 test
  entry keeps E10.5's own chamber share of blood (about 24%; the plan said about 2%), the pre-registered T3 entries
  stand on the T3 rules (the harness cannot check them), and the D2 commits are pushed. Done the same day: the T1 de_score audit
  (`experiments/harness/results/t1_de_v1/`: the handicap belongs to the E9.5 -> E10.5 change), the embryo base switch
  to pseudobulk_shift (adoption rule), and M1's direction step with a labelled diagnostic (+1.29 and +1.09 over B0;
  `experiments/t1_m1/`). The captain answered all five
  the same day (Feng on the readout; the Feng re-cut, done and identical; licence confirmed; the heart's composition
  trend dropped; the embryo's log-value edits allowed; pushed) and asked for all effort on winning. Built the same day:
  the FGW embryo pipeline (harness candidates `fgw` and `psb_size`; runs in runs/harness/emb_fgw_v1 on the DGX).
  FGW results: +4.15 over pseudobulk_shift (E7.5) and +4.92 over B0 (E7.25), 10 of 10 seeds. Test builds run end to
  end (`experiments/builds/build_test.py`; dry runs deleted). The captain approved the reading (a task the base cannot run
  on is not counted): the embryo entry is fgw. The s1 cardiomyocyte share stays at E10.5's own (26.5%). Open: the s2
  hedges (T1 blood or cardiomyocyte share; embryo shape), at the s2 pre-registration (17-23 Nov). Never build a test file at seed 0 (runbook section 5).
- Meanwhile (no specs needed): the v1 pipelines on the release data (`docs/PLAN.md` section 7, Nov 3-9): the T1
  pseudo-task with E10.5 as the truth (`experiments/t1_scope/pseudo_task.py`, honest momentum from E8.5), T2 embryo
  FGW, T2 heart growth, scored with the calibrated emulator (protocol="bake" reproduces the validation boards). On the
  DGX, cap the BLAS threads (`docs/COMPUTE.md`).
