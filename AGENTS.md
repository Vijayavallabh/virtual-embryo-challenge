# AGENTS.md

Team repository for the **Virtual Embryo Challenge** (NeurIPS 2026, https://virtualembryo.ai/challenge), entered in
**Track 1 (Human Team)**. We predict held-out embryo stages: T1 single-cell RNA (E12.5), T2 3D MERFISH expression plus
coordinates (heart E12.5, embryo E7.75), T3 a held-out knockout (beta-catenin, E8.75).

## Guiding principle (the captain, 2026-11-04)

> Treat the virtual embryo not as an image to perfect, but as a process of **becoming**. Predict which cell states
> emerge, in what proportions, how their expression changes, and, where scored, how they organize in space. A
> population is not an average cell; an embryo is not expression without structure; a perturbation is not an
> unchanged wild type.
>
> Begin with the published copy-last-stage and unchanged-wild-type baselines across the boards you intend to enter.
> Improve the weakest weighted component before polishing what already works. Let developmental biology suggest
> testable hypotheses, not excuse complexity. Retain changes that generalize to unseen stages or knockouts—not merely
> those that flatter observed data or validation scores.
>
> Ask of every component: “What meaningful change beyond the baseline does this predict, and what held-out evidence
> shows it?”

An agent applies it through the rules below, before and while changing an entry. They were written by the AI
assistant under the captain's delegation of 2026-11-02, and he can overrule them. Detail: `experiments/harness/TRIES.md`
(rules 5-8, the headroom tables) and `docs/PLAN.md` (section 3, item 9; decision log 2026-11-04).

1. **Baseline first.** Each board starts from its published baseline, which scores exactly 50: copy_last on T1:test
   and both T2 boards, wt_identity on T3:bcatenin (`docs/CHALLENGE.md`). An entry is measured by what it adds over that
   baseline, a component by what it adds over the entry without it. A copy we submit is a fresh draw and can score
   below 50 (T1: 48.06 on the chamber-cut backtest, `TRIES.md`), so beating our copy is not beating the baseline.
2. **Name the facet.** Say which prediction a step changes: which states emerge, their proportions, their expression
   or, on T2 only, their spatial organisation (T1 and T3 score no spatial term). `METHOD_SUMMARY.md` records each facet
   an entry leaves at the baseline.
3. **Answer the question before the run.** The `TRIES.md` row, committed before the run, names the predicted change
   beyond the baseline and the counting backtests or folds that would show it, none tuned or selected on (on T3 the
   cross-knockout folds). A test score is never that evidence. A reading that compares a predictor with a fold or
   backtest target is a try even when it runs from a scratch folder and never calls the scorer: its row comes first
   (`TRIES.md` rows N25 and N59 record the two breaches, 2026-11-13 and 2026-11-24).
4. **Weakest weighted component first.** Before choosing or changing work on an entry, read the headroom tables in
   `TRIES.md`. A board's next try targets its scored group with the most points left (as `headroom.py` prints it) and
   an open lever; a try elsewhere says why. Overall units, ties and closed levers: `TRIES.md` rule 6.
5. **Developmental reasons.** A biological reason makes a change testable by naming the metrics it should move. It
   never keeps a step the adoption rule rejects; an added parameter, step or data source needs its own held-out gain.
6. **Keep only what generalises.** A change is kept only if it passes the adoption rule on backtests it was not tuned
   on, with its gain where its reason said (`TRIES.md` rule 2). A tuned-on gain, one that a backtest's known feature
   explains, a public validation-board result or a better test score keeps nothing. This governs changes from
   2026-11-04 on; the entries as built stand on the rules and rulings that chose them (`docs/PLAN.md`, decision log
   2026-11-04).
7. **No safeguard loosens.** The principle opens no data, route or exception. The data rules and the post-exposure
   block (`docs/DATA_POLICY.md`, 2026-11-04, (b) to (d)), the adoption rule, the ledger's rules, the frozen s2 and T3
   rules, the captain's go-ahead for every official upload and his upload gate (209, below) all stand; a closed lever
   (`TRIES.md` rule 6) stays closed whatever its headroom.

8. **Campaign of 2026-11-09.** The captain asked for an ambitious strategy and high-risk, high-reward experiments
   (`docs/PLAN.md` decision log 2026-11-09 (a)). That lifts the ledger's limits on research effort: the prune of
   candidates expected to gain under 0.3 points, "no further try", and the lines the ledger closed after its own tries
   (blend, herding, coupling). It does not lift rule 7's data rules, the post-exposure block, A15, the frozen rules, row
   before run, try counting, the adoption rule or his ruling for any change to an upload, or his go-ahead for uploads.
   The captain renewed the campaign on 2026-11-12 with a wider brief (decision log 2026-11-12 (a), the AI assistant's
   reading of the captain's words). For that campaign heavy work may also run on this machine's otherwise-unused GPUs
   and cores, the ledger's own closure of the embryo's depth line is lifted for tries, and a candidate once pruned on
   an estimated gain may run as a counted try with its row first. Nothing else is lifted: a candidate that needs a
   standing safeguard lifted goes to the captain as a priced ruling request and changes nothing until the captain
   rules. The captain renewed it once more on 2026-11-13 (decision log 2026-11-13 (g)): nothing further is lifted, and
   campaign runs on the DGX go in a second checkout (`~/v/vec-campaign`, `docs/COMPUTE.md`), so the reference checkout
   stays as it was rehearsed. A further renewal that evening named the leading total as the target (decision log
   2026-11-13 (dd) and (ee)): again nothing is lifted, a candidate's code stays on a campaign branch until a ruling
   ports it, and where a counted backtest exists and a read costs minutes the row is written and the read is run
   rather than pruned on the AI assistant's own pre-run odds, which under-called both changes the rule has kept.

## Startup Workflow

Before writing code:

1. Confirm the working directory with `pwd`.
2. Read this file completely.
3. Read `docs/CHALLENGE.md` (tasks, metrics, weights, rules), `docs/PLAN.md` (strategy and milestones),
   `docs/DATA_POLICY.md`, and `docs/COMPUTE.md`. In the final phase also `docs/P3_RUNBOOK.md` (section 5: the builds,
   the freeze-morning rebuild, the uploads) and `docs/S2_RULES.md` (the second-upload rules; its hash is pinned in
   `docs/S2_RULES.sha256`).
4. Run `./init.sh` to verify the environment is healthy.
5. Read `feature_list.json` for the current feature state, `progress.md` for the last session and
   `session-handoff.md` for what waits; before choosing or changing work on an entry, the headroom tables in
   `experiments/harness/TRIES.md`.
6. Review recent commits with `git log --oneline -5`.

If baseline verification is failing, repair that first before adding new scope.

## Invariants (never violate)

- **Data-use rules.** Never train, tune or select on measured data inside a held-out window, by any route
  (directly, via a pretrained model, via a public atlas). The windows are listed in `docs/DATA_POLICY.md`, with two
  exceptions. The official validation files our team unlocked on 2026-10-28 may be used for preparing P3 predictions
  (Rules s.10; A20). By the captain's ruling of 2026-11-11, the cardiomyocyte or blood share in one organiser statement
  on the public Task 1 page may set T1's composition, by the rule of `docs/PLAN.md` decision log 2026-11-11 (c) and
  nothing else (`docs/DATA_POLICY.md` 2026-11-11 (a)); its one read found no such statement, so nothing is used
  (2026-11-11 (b)). Every external dataset, pretrained model or published code must
  be logged in `EXTERNAL_DATA.md` *before* it enters a model.
- **Quarantined material.** Never open veckit PyPI 0.1.2, veckit PR #6, the raw research-workflow outputs, the Yukon
  partner team's material (yukon.org), or the public board file / board endpoints (test-board anchors encode hidden
  targets). Refresh boards only with `vec boards update` (it strips test anchors in memory), and never work around its
  HTTP 403. Its `--from-file` intake reads a copy the captain saves in the browser, a route he cleared on 2026-11-01
  with conditions (notice first, not before 2026-11-04; `docs/DATA_POLICY.md`): nobody opens the saved file, and
  test-board counts are read by code only, never quoted. A test build is read by its exit code, its sha256 and its
  stamp's `commit` and `clean` fields: its log, its file size and its stamp's `size` field follow from a test board's
  count and are neither read nor quoted (`docs/DATA_POLICY.md` 2026-11-13 (c), (e)). Read research results through
  the screened digest; see
  `docs/DATA_POLICY.md`. Recursive searches skip `boards/` and `docs/research/` (`grep -r --exclude-dir=boards
  --exclude-dir=research`), and nothing searches the whole filesystem. Agents never open session logs, transcripts,
  harness-saved tool outputs or memory files (under `~/.claude-account` or `~/.claude`), and run the repository's
  Python as `python -I -B` (`docs/DATA_POLICY.md`, 2026-11-04 (g)). Read the pinned scorer's `common/core_metrics.py`
  for its formulas only, and skip its module docstring and its `de_genes` docstring (example counts that may describe
  a T3 board; `docs/DATA_POLICY.md` 2026-11-05 (i), (k)). Read `docs/reference/neurips-2026-competition-proposal.pdf`
  skipping its Data section and its page-7 evaluation on held-out knockout data (`docs/DATA_POLICY.md` 2026-11-06
  (a)). Revisions of `docs/PLAN.md`, `docs/DATA_POLICY.md`,
  `docs/research/digest_2026-10-19.md` and `docs/research/t3_methods.md` older than their cleanup of 2026-11-05
  (`docs/DATA_POLICY.md` 2026-11-05 (h), (i)) hold screened-out material (per-board leaderboard figures, a count that
  may describe a hidden target): never print them (`git show` of an older revision, `git log -p`, or a diff across the
  cleanup on those paths); a review of such a commit masks the removed side.
- **Scores choose, never compute.** Leaderboard returns may be used to pick between predictions, never to solve for a
  property of a hidden target (no inverting the skill map, no differencing scores). A screened leaderboard read hands
  the main agent overall totals and overall ranks only (`docs/DATA_POLICY.md`, Scores).
- **Competition data stays private.** Data lives under `data/` (gitignored) and is never committed, shared or shown
  publicly. The dataset is unpublished (Neil Chi lab, UCSD): no public presentation without written approval. P3
  methods, prediction files and evaluation feedback stay private from other participants until the 2026-12-11
  announcement. Any public release of the code is a fresh repository built from the final tree, never this repository or
  its history (the captain, 2026-11-05).
- **Official submissions need the user's explicit go-ahead.** The final phase allows only 2 scored submissions per
  board for the whole phase, and they cannot be withdrawn. The API treats an upload as **official unless it sends
  `format_only=true`**. The team API key (`~/.config/vec/token`, mode 600, local machine only, never on the DGX) can
  make official uploads: use it only through this repo's own tools, never with community "submit" tools or ad-hoc
  scripts. Validate files with `vec check` or format-only uploads; real uploads go through the submission client
  (feat-016), which defaults to format-only. An official upload sends only a file whose sha256 is in the current
  reference list (since the third rehearsal of 2026-11-13 `experiments/builds/results/builds_dc7b6ae.SHA256SUMS`,
  the DGX's `runs/archive/builds_dc7b6ae/SHA256SUMS`), taken from the folder the freeze-morning command builds and never from
  the DGX's `runs/builds` or `runs/t3_build`, or from an earlier archive folder (`builds_38616fe` and `builds_fd2683a`
  hold the two heart files the ruling of 2026-11-13 replaced), which hold older files whose clean stamps the client accepts
  (`docs/P3_RUNBOOK.md`, "The freeze-morning rebuild"); if the DGX cannot be reached, the checked copy of the same
  bytes on the local machine can go, at the captain's word. A file that differs from the list is not uploaded until the
  captain has the explanation. The client enforces the list since 2026-11-13: its guard "reference list" refuses an
  official upload unless the file's sha256 is in that list under one of the board's names, read through the pointer
  `results/reference_list.json`; a rehearsal that writes a new list rewrites the pointer (the list's path and
  sha256, the boards' names) in the same commit.
- **The captain's upload gate (2026-11-05).** In his words: "Until our estimated score in leaderboard crosses 209, don't
  submit our 1st submission." No first official upload on any board (T3's first is its s2) while our estimated
  leaderboard total is 209 or below. The estimate is the overall total the leaderboard would show for our entries, each
  board at the better of its two uploads, from our own committed evidence by the procedure in `docs/PLAN.md` section 8,
  recomputed by an independent reviewer before the gate counts as crossed; never from a test score or a leaderboard
  figure. The figure stays 209 unless the captain changes it, and only he lifts the gate. His deadline exception
  (2026-11-05): "I am fine with submitting our 1st submission at Nov 27 anyway": the s1 methods freeze at 2026-11-27
  00:00 UTC, and the first round may then go in whatever the estimate, each upload only on his go-ahead at the time; he
  confirmed this reading and the better-of-two estimate on 2026-11-05. The gate is about uploads only, never a try
  (`docs/DATA_POLICY.md` 2026-11-05). The record is `results/upload_gate.json`; the submission client refuses an
  official upload while it does not pass (guard "upload gate").
- **Heavy compute runs on the DGX** (`PrakashDGX_H2`, workspace `~/v/virtual-embryo-challenge`), never on the shared
  local machine. See `docs/COMPUTE.md`. One exception, by the captain's instruction of 2026-11-12: the campaign of
  that day (rule 8) may use this machine's otherwise-unused GPUs and cores, with `data/` and `runs/` here owner-only.
  Never sync, pull or check out on the DGX while a build runs there (the build stamps then say unclean). Never check a
  file built on one machine against a hash made on the other: three T3 files come out with a few cells drawn
  differently on the two machines (`docs/COMPUTE.md`, "What reproduces where").
- **uv only.** Use `uv add`, `uv sync`, `uv run`. No pip installs into the env, no conda.
- **Human Team track.** A person reviews results and makes the calls; there is no Agent-track evidence to keep.

## Working Rules

- **One feature at a time**: pick exactly one unfinished feature from `feature_list.json`.
- **Shared core**: the T1 and T2 heads build every prediction through `vec.core`, and the pre-registered T3 script
  runs identical copies of its primitives, pinned by `tests/test_core.py` (the Generality Award claim, `docs/PLAN.md`
  section 5). New steps (coupling, births, calibration) go into the core or a head, never into a separate per-task
  model. A change meant to alter a T1 or T2 output updates the golden table in `tests/test_core.py` in the same commit
  and says why.
- **Evidence over claims**: a model result counts only with a local-emulator score file (`runs/<run_id>/scores.json`)
  and a row in `results/registry.csv`.
- **Update artifacts**: before ending a session, update `progress.md` and `feature_list.json`.
- **Stay in scope**: don't modify files unrelated to the current feature.
- **Leave clean state**: the next session must be able to run `./init.sh` immediately.

## Repository Map

- `src/vec/`: package (`io.py` loading and submission writing, `boards.py` board specs, `evaluate.py` server-faithful
  local board emulator, `baselines.py` reference baselines, `cli.py` the `vec` command, `core.py` the shared core,
  `t1_scope.py` and `t2_build.py` the T1 and T2 heads).
- `configs/`: experiment configs. `runs/` (gitignored): per-run outputs; on the local machine also the checked copy
  of the reference upload files and of the T3 run folder (`docs/COMPUTE.md`). `results/registry.csv`: one row per
  scored run; `results/submissions.csv`: the upload ledger; `results/upload_gate.json`: the captain's gate;
  `results/reference_list.json`: the pointer to the current reference list, read by the client's guard.
- `experiments/`: one folder per study (scripts, reports, committed `results/`); `docs/P3_RUNBOOK.md` step 4 names the
  pre-registered ones, and `labels/` holds the classifier evidence for the validation answers (step 3, 2026-10-28).
- `experiments/harness/`: `TRIES.md` (the tries ledger and the headroom tables per test board), `headroom.py` (writes
  those tables' numbers to `results/headroom/` from committed results), `run.py` (harness runs and adoption decisions).
- `experiments/builds/`: the test-file builders (`build_test.py` for T1 and T2; T3's is
  `experiments/t3_decomp/build_submission.py`), `build_all.sh` (one submission set, built and checked),
  `determinism_check.sh` (rebuild into a fresh folder and compare with a reference list: the freeze-morning command)
  and `results/` (rehearsal records and the reference list, `builds_dc7b6ae.SHA256SUMS` since the third rehearsal
  of 2026-11-13).
- `boards/`: public board specs from the challenge site (`index.json`, panel gene lists, T1 composition). Since 634b5fd
  `index.json` also holds the test boards' cell limits: code reads it; nobody opens it or quotes a test board's limits
  (`docs/DATA_POLICY.md`).
- `data/` (gitignored): `raw/` competition files, `external/` disclosed external data.
- `docs/`: `CHALLENGE.md`, `PLAN.md` (strategy, timeline, decision log), `DATA_POLICY.md`, `COMPUTE.md`,
  `P3_RUNBOOK.md` (the final phase, step by step), `S2_RULES.md` (the second-upload rules, hash-pinned),
  `ORGANIZER_QUESTIONS.md` (every question sent and the reading of each answer fixed beforehand), `DATA_NOTES.md`,
  `research/` (background reports; screened digests only).
- `EXTERNAL_DATA.md`: disclosure register (required with every submission). `METHOD_SUMMARY.md`: running method summary
  (required for the winners' report and the Generality Award).

## Definition of Done

A feature is done only when ALL of the following are true:

- [ ] Target behavior is implemented.
- [ ] `./init.sh` passes (lock, lint, tests).
- [ ] Evidence recorded in `feature_list.json` or `progress.md` (commands run, score files, run IDs).
- [ ] Any new external source is logged in `EXTERNAL_DATA.md`.
- [ ] A change to an entry has its `experiments/harness/TRIES.md` row, committed before its run, and, if kept, a
      rerun of `experiments/harness/headroom.py` with its output and the copied tables committed, and updated
      `METHOD_SUMMARY.md` lines (its change beyond the baseline and the evidence).
- [ ] The repository remains restartable from the standard startup path.

## End of Session

1. Update `progress.md` with the current state and `session-handoff.md` with what waits.
2. Update `feature_list.json` with the new feature status and evidence.
3. Record unresolved risks or blockers.
4. Commit with a descriptive message once work is in a safe state.

## Verification Commands

```bash
# Full verification (recommended)
./init.sh
```

Required checks:
- `uv sync --locked`
- `uv run ruff check src tests`
- `uv run pytest -q`

`./init.sh` checks the environment, not the entries. It pins each task's published floor (`tests/test_boards.py`);
the headroom tables in `experiments/harness/TRIES.md` show where each entry stands against that floor, copied from a
rerun of `experiments/harness/headroom.py` after a kept change (`TRIES.md` rule 8).

## Escalation

- **Rule ambiguity** (external data near a window, cross-task data use): stop and ask the user; draft the question
  for the organisers (virtual.embryo.moonshot@gmail.com or Slack #virtual-embryo-challenge).
- **Architecture decisions**: consult `docs/PLAN.md`, otherwise ask the user.
- **Repeated test failures**: update `progress.md` and flag for human review.
- **Anything touching an official submission**: ask the user first.
