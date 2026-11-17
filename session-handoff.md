# Session Handoff

## Current Objective

- **Read this first (2026-11-17).** The organisers moved the final deadline to **25 November 2026, 23:59 AoE (26
  November 11:59 UTC, 17:29 IST)**, for uploads and final-entry marks alike. Every later date in this file, in
  `docs/PLAN.md` section 7, `docs/P3_RUNBOOK.md`, `docs/S2_RULES.md` and `AGENTS.md` is superseded in fact but not
  yet in text: the dates are the captain's to move, and the captain has been asked (proposal: freeze and first round
  21 November from 00:00 UTC, second round 23 November, T3's last upload 24 November, marks 25 November by 12:00 UTC,
  any file change rehearsed by 20 November). **The upload gate's exception (27 November) is after the deadline: until
  the captain moves it the client refuses every official upload.** When the captain names the dates, one commit
  changes `results/upload_gate.json` (with the captain's words), its test in `tests/test_submit.py`, an amendment at
  the end of `docs/S2_RULES.md` with a new `docs/S2_RULES.sha256`, the runbook tables and `AGENTS.md`. Decision log
  2026-11-17 (a) has the whole list. The organisers also answered the exposure notice (accepted; they ask for the
  explanation and the precommitted procedure in the method report): a further reply to the thread's other questions
  is announced and must go through a screening agent, never be pasted.
- **What is live (2026-11-17, 14:00 UTC).** No candidate change is live. The step's request was read on both
  sides, attacked cold twice and is not put forward (decision log 2026-11-17 (b) to (e)); the AI assistant's advice
  is to keep 0.40 on both Task 1 uploads. Rows T48 to T54 are in the ledger with their outcomes (T1 93 tries).
  The DGX holds four idle
  campaign checkouts (`~/v/vec-campaign`, `-h7`, `-c`, `-d`; identical package sets; the reference checkout is
  untouched at dc7b6ae). The builder option for a step ruling stays on `campaign/t1-n` (2921b96), unused.
  Try T55 (`part_read.py`) runs in `~/v/vec-campaign-c` at e492c30 (exit file `runs/round6/t55.exit`; finish it
  with `lift6/finish_t1read.sh vec-campaign-c t55_m1_part t55 t55_logs e492c30`, commit unread, read once, outcome,
  totals then T1 100). The schedule amendments are drafted on `draft/schedule-nov21-fix` (490a51a, reviewed: ready,
  not merged; decision log 2026-11-17 (f)): fill its markers with the captain's words, recompute the rules hash,
  rebase and merge.
  **For the captain, in this order:** (1) the dates: the gate's exception (27 November) is after the organisers'
  deadline, so the client refuses every official upload until the captain moves it (proposal: freeze and first
  round 21 November 00:00 UTC, second round 23 November, T3's last upload 24 November, marks 25 November by 12:00
  UTC); (2) the six confirmations for the method report (next bullet); (3) for information, each with the AI
  assistant's advice to change nothing: the size of Task 3's factor (decision log 2026-11-17 (b)) and the size of
  Task 1's step (entry (e); if the captain wants it all the same, the least exposed form is the first upload only
  at 0.55, left in writing to the pre-set rule); (4) the older open decisions (Task 1's read depth: the rule gives
  declined; the embryo's moved-less family and the late-fetal practice family: advice no). Owed by the AI assistant
  after T52: the full suite on the DGX at the last commit, and, once the captain names the dates, the amendments
  (the gate's file with the captain's words and its test, an amendment at the end of `docs/S2_RULES.md` with a new
  hash, the runbook's tables, `AGENTS.md`), the freeze-morning rebuild and the disclosure email.
- **The method report's exposure section (2026-11-17).** `METHOD_SUMMARY.md` now holds "Exposure of hidden-target
  metadata on 2026-11-04 and the precommitted procedure", which the organisers asked for. It was drafted by an agent
  from the dated records, reviewed by three others (accuracy, framing, tone) and revised once; the AI assistant
  checked the commits it cites. Before the report leaves the repository the captain confirms what only the captain
  can: (1) what the captain saw and did on 2026-11-04 (the address opened in a browser, the index and then the
  signed-in Data page pasted, about 07:05 UTC) and that no copy exists outside that session's local log; (2) that the
  captain's choices of 2026-11-05, 2026-11-06, 2026-11-11 and 2026-11-13 (rows R4 and R7) did not draw on the figures
  seen; (3) that the email of 2026-11-04 went out as recorded and that the quoted answer is the one received; (4)
  whether the correction on Task 3's count default (set in code on 2026-11-01, after the release; the email said
  "before") should reach the organisers through this section, as it now does; (5) whether to answer the open question
  on Task 1's read depth first (the section says it is unanswered and that the rule gives "declined"); (6) whether
  older commits that still hold the Task 1 E12.5 file size are to be rewritten (a force-push; the captain's call).
  Not traced to a dated entry by the drafter and left as the record's own summary: who pruned the Task 1 additions
  of 2026-11-04, the date the chamber-only scope file was fixed, and the reason the embryo's depth is outside the
  block. Section 5 of the brief (what the organisers were told) quotes their answer in full.
  One bullet was added to the section after its review, on 2026-11-17: the step reads of that day (tries T44 to
  T54) and that no request followed. It has not been through the three reviewers.
- Goal: a compliant Human-track entry on all four test boards that finishes as high as it can. The final phase now
  ends on 2026-11-25 at 23:59 AoE, with two scored uploads per board for the whole phase; none is used as of
  2026-11-17.
- Current status (2026-11-14): the entries are built and their methods chosen. T1: M1 with R1, R2 and R3b, its first
  upload drawn from E10.5's dissected-heart library since 2026-11-13 and its second on the pooled base at 60%
  cardiomyocytes. Heart: B0 with the pre-registered size and, since the captain's ruling of the evening of
  2026-11-13 (row H11), each lineage's expression step taken from its trend over all released heart stages.
  Embryo: fgw, with `fgw_aligned_series` as second upload
  since 2026-11-13. T3: the frozen rules' files, the second upload now with a wild-type-only factor. Nothing is
  uploaded. The captain's upload gate holds every first upload until the team's own estimate passes 209 (it is 184,
  range 170 to 198) or the deadline exception opens at 2026-11-27 00:00 UTC; each upload then needs the captain's
  go-ahead at the time, and the captain runs the commands.
- No method ruling is open (`docs/PLAN.md` decision log 2026-11-13 (d), (k) and (w); `docs/S2_RULES.md`, its last
  three sections). After the captain's rulings of the evening of 2026-11-13 (T1's second upload on the mid base,
  row R7; both heart files with the trend, row H11) the ten files were rehearsed a third time on the DGX (decision
  log (bb)): three byte-identical builds, and **the freeze-morning reference is now
  `runs/archive/builds_dc7b6ae/SHA256SUMS` there** (`experiments/builds/results/rehearsal_dc7b6ae.txt`,
  `builds_dc7b6ae.SHA256SUMS`; `builds_38616fe` and `builds_fd2683a` are superseded and hold the replaced heart
  files). Every file of the set has passed the server's format check (six on 2026-11-13 at the captain's word,
  decision log (x)). The upload client refuses a scored upload of a file that is not on this list (guard
  "reference list", pointer `results/reference_list.json`; decision log (cc)).
  What is owed:
  - DONE on 2026-11-14 (decision log 2026-11-14 (c), item 1: the family closed with no request; kept here as the
    record of what was owed): the controls and checks of tries E41b. The embryo arm `fgw_partial` (branch `campaign/emb-f`, fcc1f6b)
    passed the adoption rule and the first-upload pair rule on 2026-11-14 (ledger row E41;
    `experiments/harness/results/n45_embryo_f_readings/`). Before any request to the captain: (1) build, review and
    read the three control arms the captain asked for (every anchor moved by the same share; the same number left
    unmoved at random; no movement), row first, on the nine practice tasks; the brief is drafted (session scratchpad,
    `lift6/briefs/embg.md`) and nothing of it is built; (2) rerun the cold adversarial check and the independent
    recomputation of the reading's tables, both of which stopped on the account's usage limit on 2026-11-14 without a
    result; (3) then a priced request, both ways, for the simplest arm that passes, or a recorded decline. A yes needs
    a builder port, a dry build, a fourth rehearsal and a new reference list by about 24 Nov. Then the seventh
    independent pass over decision log (ii), lines 13 to 18;
  - nothing to the organisers: on 2026-11-14 the captain had their logical replies assumed for the questions still
    open (decision log 2026-11-14 (a); `docs/ORGANIZER_QUESTIONS.md`, "Working answers (2026-11-14)"). The replies
    confirm the files as built; the items below about the eighth email, the Task 3 batch question and the screened
    reads stand only as the captain's options and as checks;
  - the captain's decision on a drafted question, if the captain wants it sent: an eighth email asking the organisers
    how final entries are marked (per board or one pick over several boards; what counts on an unmarked board; until
    when a pick can change). It is reviewed and staged in `docs/ORGANIZER_QUESTIONS.md` with its readings fixed
    beforehand, and a draft with placeholders for the name and the team is in the captain's mailbox. The per-board
    answer is a working answer rated 0.7 (the table of 2026-11-10); the estimate assumes it;
  - nothing on request R-D unless the captain wants it (`docs/PLAN.md` decision log 2026-11-13 (ff) and (gg); rows
    N37, N38 and N39 of `experiments/harness/TRIES.md`). The request: change the embryo's first upload from `fgw`
    to `fgw_trend_render` (the same cells and positions, each cell's expression moved to the profile a trend
    through the four released embryo samples gives at E7.75), second upload unchanged; tabled at +0.93 embryo
    points at the middle weighting by a pre-set pair rule on the nine practice tasks. **The AI assistant advises
    no**: a cold check and reading N39 showed that the gain follows the fit geometry and that the test's geometry
    is the losing kind where the truth is like its earlier bracket (fitted gain -0.60, uncovered). The request
    stays on record as the rule's outcome; if the captain rules for it all the same, it needs the builder port
    (prepared on the branch `port/emb-s1-trend`, off main, and dry-built once at 394f59f with exit 0, decision log
    (ii); a new dry build only if the port's commit changes), a fourth rehearsal, a new reference list
    and pointer, the server's format check of the new file at the captain's word, and the summary and disclosure
    lines;
  - nothing on the runner-up either unless the captain wants it: `fgw_aligned_trend` also passed the pair rule as a
    first upload (+0.59 embryo points at the middle weighting, reading N37). Reading N40 found its gain follows the
    same fit geometry (fitted -0.53 at the test's earlier-bracket geometry), and the AI assistant advises no; that
    lens was fixed after the scores existed, so it is advice and not a pre-set gate (decision log (ii));
  - on 19 and 23 Nov, the approved screened reads for an organiser answer to the seventh email
    (`docs/SCREENED_READ_SEVENTH_EMAIL.md`; a reply in the captain's inbox is not pasted into the session). One
    consequence was fixed on 2026-11-13 (decision log (hh)): if the answer is that the hidden embryo sample was
    hybridised and imaged with E7.5, and final marks are confirmed to be per board, an earlier-side first upload
    for the embryo goes to the captain as a priced request, after a target-free check of the test pair's factor;
    on any other answer, or none, nothing changes;
  - about 19 Nov, the Slack follow-up draft on floored T3 severity if no organiser has answered; with it, for the
    captain to keep or cut, the Task 3 batch question that was cut from the seventh email in review (was the
    beta-catenin knockout processed with the Gata4 replicates): an answer changes no file and would move the
    estimate by about 7 points either way (decision log 2026-11-13 (ee));
  - a new rehearsal and reference by 24 Nov only if a build changes after 2026-11-13 (and then the pointer
    `results/reference_list.json` is rewritten for the new list);
  The generator tests the captain ordered on 2026-11-13 are done (decision log (o) to (r); rows T34, T35, H8, H9,
  H10; results `experiments/harness/results/gen_t1_v1/` and `gen_heart_v1/`): five counted reads, nothing kept, no
  upload, build or rule changed. The code is on the branch `campaign/gen` (pushed; 459916f), off main. The heart
  tries that followed (rows H11, H11s, H12, H12s; branch `campaign/trend`, pushed): `b0_trend` passed the adoption
  rule and was ruled into both heart files; its size-matched twin and the two parts of its factor read alone were
  not kept. The sixth round that night (decision log (dd), (ee) and (ff)) kept nothing by the adoption rule: its
  first part (rows N33, N28, N29, N32, H13, K9, K10) ran three heart arms and two T3 readings and restated the
  forecast for the files as built (mean 194.4, the upper reading; the gate figure stays 184); its second part read
  every remaining allowed line (26 tries in 14 rows: H14 to H20, E29 and E31 to E35, T37, T38, K11) and produced
  request R-D above, which the AI assistant advises against. The options for taking our own estimate past the
  leading total were declared used up in decision log (ff) and (gg), too early: the round's audit ((ii)) then read the
  eight lines that independent passes over the record found open and the closures they found resting on judgement
  (tries T39 to T41, K12, K13 and E36 to E39; reading N40), and replaced every closure that rested on the AI
  assistant's own odds by a row-first read, a pre-set gate or arithmetic, or left it to the captain: nothing passes,
  nothing changes, and two critics found no open lever. Tries: T1 37, embryo 60, heart 25, T3 27.
  Settled on 2026-11-13 (decision log (k)): T1's second upload (the mid base; built in and rehearsed), the heart
  try (run, not kept), and the checked copy of the reference files on this machine (kept, by the captain's word).
- Branch / commit: `main` on GitHub (private). The DGX reference checkout
  (`PrakashDGX_H2:~/v/virtual-embryo-challenge`) is at dc7b6ae, the third rehearsal's commit, and stays untouched
  until the freeze-morning sync; campaign runs use `~/v/vec-campaign` and `~/v/vec-campaign-h7`, which sit at
  campaign-branch commits after the sixth round. The round's candidate code is on pushed branches named in decision
  log (ff), and the audit's on `campaign/t1-c`, `campaign/t1-d`, `campaign/t1-e`, `campaign/emb-c`,
  `campaign/emb-polar`, `campaign/emb-mirror`, `campaign/emb-d`, `campaign/t3-c` and `campaign/t3-d`; none of it is
  on main.

## Completed This Session (2026-11-12 to 2026-11-13: the campaign the captain asked for on 2026-11-12)

- [x] Sixth round, second part (2026-11-13, late night; decision log (ff)). Every remaining allowed line with a
      counted read: 26 tries in 14 rows on arms built and reviewed on campaign branches (heart H14 to H20, embryo
      E29 and E31 to E35, T1 T37 and T38, T3 K11), gate N34, readings N37, N38 and N39. Nothing passes the adoption
      rule; eight lines could not be built from the record; request R-D came out of reading N37 and the AI
      assistant advises against it after a cold check and reading N39 (decision log (gg)). A second option, an
      earlier-side file as first upload, was attacked cold before any advice and declined (decision log (hh)).
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
      copy of the reference set (`runs/archive/builds_fd2683a/`; since the second rehearsal
      `runs/archive/builds_38616fe/`) and of the T3 run folder
      (`runs/t3_decomp/gata4_20261028T122425Z/`), owner-only. The client's build-stamp guard and local format check
      pass on all ten copies (nothing sent), so an upload can go from the copy if the DGX is down, at the captain's
      word (`docs/P3_RUNBOOK.md`, "A copy off the DGX"). Every file also builds here now: seven of ten byte-identical
      to the DGX's, the three T3 factor files with a few cells drawn differently (3, 84 and 14 cells).
- [x] The harness files brought up to date a second time on 2026-11-13 (decision log (f)), at the captain's request:
      `AGENTS.md` (final-phase reading list, how a test build is read, where an upload's file comes from, the DGX
      sync rule, the repository map), `docs/COMPUTE.md` ("What reproduces where"), `README.md`, `init.sh`,
      `ASSUMPTIONS.md` (A28), `feature_list.json`, `METHOD_SUMMARY.md`, ledger rows R4 to R6.

- [x] The campaign's fifth round (2026-11-13; decision log (g) to (j)), after the captain renewed it through the
      ambition-lift skill. A fourth screened leaderboard read (sixth place 203.00, 124 Human entries); a cold
      review that stopped the first plan; a forecast of the total as a distribution
      (`experiments/gate/forecast.py`); tries T32 and T33 on T1's pair, both recomputed independently; the heart's
      target-free gate N26, which failed, so no heart try ran. Tries: T1 30, embryo 36, heart 4, T3 15. Campaign
      runs on the DGX now use a second checkout (`~/v/vec-campaign`); no GPU job ran.
- [x] The captain's answers of 2026-11-13 (decision log (k) to (m)) and what followed. T1's second upload moved to
      the mid base (row R7): builder option with tests and an independent code review, rules amended in the
      captain's words (new sha256 23ae5197...), `build_all.sh` changed, all ten files rehearsed on the DGX at
      38616fe (three byte-identical builds; nine equal to the first rehearsal's, T1_s2 new), the reference list and
      the checked local copy refreshed. The heart try H7 with its control was run and not kept (heart 4 to 6
      tries). The checked local copy is confirmed.

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

- Current: `docs/PLAN.md` decision log through 2026-11-13 (f), `docs/DATA_POLICY.md` change notes through
  2026-11-13 (e), the second-upload rules in `docs/S2_RULES.md` as amended on 2026-11-13, and the captain's upload gate in
  `results/upload_gate.json`. The campaign's method decisions are the AI assistant's, logged; every change to an
  upload is the captain's ruling.
- See `docs/PLAN.md` decision log (2026-10-18, 2026-10-19 a-d, 2026-10-20 a-c, 2026-10-21 a-c, 2026-10-28) and
  `docs/DATA_POLICY.md` records. The user's calls of 2026-10-27 and 2026-10-28 (the history rewrite, the unlock, the
  A20 reading, the go-ahead for step 3) are recorded in `progress.md` (Decisions Made, Blockers / Risks).

## Blockers / Risks

- 2026-11-13 (late night, after decision log (gg)): request R-D rests on nine practice tasks of which three share
  one target stage and two are the counted backtests that prompted the read, and reading N39 shows its gain does
  not carry to the test's geometry; the AI assistant advises against it. If it were ruled in all the same, the new
  first upload is a rendered file, which reproduces only at the build's 16 BLAS threads, and at the test one
  sample's file lacks two panel genes (they stay as `fgw` has them, a rule fixed before any build).
- 2026-11-13 (night, after the sixth round): in our own evidence a total past the leading 215.41 lies beyond the
  best case of every board at once (their best practice outcomes sum to 219.0), and the forecast's upper reading
  (194.4) leans on one judgement that the record itself rates lower: that the beta-catenin samples share the Gata4
  slide's profile (weight one half; at 0.15 the mean is 189.0). The committed 184 is the figure to quote.
- 2026-11-13: the estimate (184) is about 25 under the gate, and no change was adopted in the campaign. What can
  still move the outcome is how the hidden samples were measured (the bets behind R-A, R-B and T3's first-choice
  file), which no released file shows.
- 2026-11-13: no organiser answer is recorded to any of the seven emails, the Slack DM or the public Slack
  question. The next screened check is on 26-27 Nov; one follow-up in the Slack thread may be drafted about 19 Nov.
- 2026-11-13: a ruling that changes a build is needed by about 24 Nov (a T3 builder change takes about a day, then
  a dry build and the rebuild rehearsal). A rendered embryo file reproduces only at the build's 16 BLAS threads.
- 2026-11-13: only DGX hashes check DGX builds. The three T3 files that apply a per-gene factor have other bytes on
  the local machine (as T1's s1 had until that day), so a change of the DGX's numerical libraries or processor before 27 Nov would show
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
3. Review this handoff, `docs/P3_RUNBOOK.md`, and `docs/PLAN.md` decision log 2026-11-13 (a) to (f).
4. Run `./init.sh` before editing. On the DGX do not run it as it is (`docs/COMPUTE.md`, `progress.md` notes).
5. Before choosing or changing work on an entry, read the guiding principle in `AGENTS.md` and the headroom tables in
   `experiments/harness/TRIES.md`.

## Recommended Next Step

- 2026-11-14 (latest, after decision log (ii)): the round's audit is done. Every line that independent passes over the
  record found open is read, and no closure rests on the AI assistant's own odds; nothing passes and no file changes.
  Three things are the captain's, all optional: a ruling on request R-D or on its runner-up `fgw_aligned_trend` (the
  AI assistant advises no on both, and says the advice is an after-the-fact lens), whether to send the drafted eighth
  email on final marks, and whether the Task 3 batch question goes into the follow-up about 19 Nov. No run is queued,
  and the AI assistant knows of none whose result could change an upload under the standing rules. Lessons of the
  round, for the next brief: a check that reads practice scores gets its row first; a gate is tried on pure noise
  before it is trusted; a control is checked for whether it can differ from its arm; a brief names where a formula's
  translation goes; and a runner records an arm the builder refuses instead of stopping. The dated items follow as
  listed below. Report the expected total as 184 (170 to 198), with 194.4 named as the forecast's upper reading.
- 2026-11-13 (latest, after decision log (gg)): the sixth round is closed and the options for passing the leading
  total are used up; every allowed line with a counted read has been read. **No ruling is needed for the entry
  as built; one optional ruling is on record for the captain.** Request R-D (the embryo's first upload to
  `fgw_trend_render`) is on record as a pair rule's outcome
  with the AI assistant's advice against it (decision log (gg), reading N39); it changes nothing unless the captain
  rules for it, and then: port the branch `port/emb-s1-trend` to main after its review, set the embryo's
  first-upload expression in `experiments/builds/build_all.sh`, dry-build, rehearse a fourth time on the DGX, write
  the new reference list and the pointer `results/reference_list.json`, ask the captain before the one server
  format check of the new file, and update `METHOD_SUMMARY.md`. Before recommending any request built on practice
  tasks: name the property of a task the effect depends on, compute it for the test from target-free quantities,
  and run the cold check first (decision log (gg)). The dated items follow: the screened reads on 19 and 23 Nov,
  the follow-up draft about 19 Nov, the method summary and disclosure email on 24-30 Nov, and the freeze morning.
  Report the expected total as 184 (170 to 198), with 194.4 named as the forecast's upper reading.
- 2026-11-13 (latest, after decision log (ee)): the sixth campaign round is closed with nothing kept and nothing
  requested. No run is queued on any GPU, and the AI assistant knows of none whose result could change an upload
  under the standing rules. Next are the dated items: the screened reads on 19 and 23 Nov, the follow-up draft
  about 19 Nov (floored severity, and the Task 3 batch question if the captain wants it back), the method summary
  and disclosure email on 24-30 Nov, and the freeze morning. Report the expected total as 184 (170 to 198), with
  194.4 named as the forecast's upper reading. Where a counted backtest exists and a read costs minutes, write the
  row and run it rather than pruning on pre-run odds (`AGENTS.md` rule 8).
- 2026-11-13 (latest, after decision log (cc)): the heart ruling is built in and rehearsed (reference
  `runs/archive/builds_dc7b6ae/SHA256SUMS`; every file format-checked on the server; the upload client refuses a
  scored upload of a file off that list). Next: the dated items below, reading `builds_dc7b6ae` for the reference
  and skipping the server format checks, which are done. No run is queued on any GPU. On the upload mornings the
  12 characters to confirm come from the list: heart 3d6548382b24 (first upload) and 18c2cf552f98 (second).
- 2026-11-13 (after decision log (m)): the captain ruled, T1's second upload is on the mid base and the
  ten files are rehearsed again (reference `runs/archive/builds_38616fe/SHA256SUMS`, since superseded). The
  generator tests the captain ordered are done and changed nothing (decision log (r)).
- 2026-11-13 (earlier that day): the rulings are given and built in, and the ten files are rehearsed on the DGX. Next, in
  order: the screened reads on 19 and 23 Nov and the Slack follow-up draft about 19 Nov; the method summary and the
  disclosure email (with the hashes of `runs/archive/builds_38616fe/SHA256SUMS`, the list since the second
  rehearsal) on 24-30 Nov; the screened check on
  26-27 Nov; on 27 Nov the DGX sync to the freeze commit, the freeze-morning command (`determinism_check.sh
  runs/archive/builds_38616fe/SHA256SUMS 16`, the list since the second rehearsal), the server format checks and
  the first round, each upload on the
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
