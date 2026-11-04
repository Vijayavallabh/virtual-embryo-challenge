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
   cross-knockout folds). A test score is never that evidence.
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
   rules and the captain's go-ahead for every official upload all stand; a closed lever (`TRIES.md` rule 6) stays
   closed whatever its headroom.

## Startup Workflow

Before writing code:

1. Confirm the working directory with `pwd`.
2. Read this file completely.
3. Read `docs/CHALLENGE.md` (tasks, metrics, weights, rules), `docs/PLAN.md` (strategy and milestones),
   `docs/DATA_POLICY.md`, and `docs/COMPUTE.md`.
4. Run `./init.sh` to verify the environment is healthy.
5. Read `feature_list.json` for the current feature state and `progress.md` for the last session; before choosing
   or changing work on an entry, the headroom tables in `experiments/harness/TRIES.md`.
6. Review recent commits with `git log --oneline -5`.

If baseline verification is failing, repair that first before adding new scope.

## Invariants (never violate)

- **Data-use rules.** Never train, tune or select on measured data inside a held-out window, by any route
  (directly, via a pretrained model, via a public atlas). The windows are listed in `docs/DATA_POLICY.md`, with the
  one exception: the official validation files our team unlocked on 2026-10-28 may be used for preparing P3
  predictions (Rules s.10; A20). Every external dataset, pretrained model or published code must be logged in
  `EXTERNAL_DATA.md` *before* it enters a model.
- **Quarantined material.** Never open veckit PyPI 0.1.2, veckit PR #6, the raw research-workflow outputs, the Yukon
  partner team's material (yukon.org), or the public board file / board endpoints (test-board anchors encode hidden
  targets). Refresh boards only with `vec boards update` (it strips test anchors in memory), and never work around
  its HTTP 403. Its `--from-file` intake reads a copy the captain saves in the browser, a route he cleared on
  2026-11-01 with conditions (notice first, not before 2026-11-04; `docs/DATA_POLICY.md`): nobody opens the saved
  file, and test-board counts are read by code only, never quoted. Read research results through the screened digest;
  see `docs/DATA_POLICY.md`. Recursive searches skip `boards/` and `docs/research/`
  (`grep -r --exclude-dir=boards --exclude-dir=research`), and nothing searches the whole filesystem. Agents never
  open session logs, transcripts or harness-saved tool outputs (under `~/.claude-account` or `~/.claude`), and run
  the repository's Python as `python -I -B` (`docs/DATA_POLICY.md`, 2026-11-04 (g)).
- **Scores choose, never compute.** Leaderboard returns may be used to pick between predictions, never to solve for a
  property of a hidden target (no inverting the skill map, no differencing scores).
- **Competition data stays private.** Data lives under `data/` (gitignored) and is never committed, shared or shown
  publicly. The dataset is unpublished (Neil Chi lab, UCSD): no public presentation without written approval. P3
  methods, prediction files and evaluation feedback stay private from other participants until the 2026-12-11
  announcement.
- **Official submissions need the user's explicit go-ahead.** The final phase allows only 2 scored submissions per
  board for the whole phase, and they cannot be withdrawn. The API treats an upload as **official unless it sends
  `format_only=true`**. The team API key (`~/.config/vec/token`, mode 600, local machine only, never on the DGX) can
  make official uploads: use it only through this repo's own tools, never with community "submit" tools or ad-hoc
  scripts. Validate files with `vec check` or format-only uploads; real uploads go through the submission client
  (feat-016), which defaults to format-only.
- **Heavy compute runs on the DGX** (`PrakashDGX_H2`, workspace `~/v/virtual-embryo-challenge`), never on the shared
  local machine. See `docs/COMPUTE.md`.
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
- `configs/`: experiment configs. `runs/` (gitignored): per-run outputs. `results/registry.csv`: one row per scored run.
- `experiments/`: one folder per study (scripts, reports, committed `results/`); `docs/P3_RUNBOOK.md` step 4 names the
  pre-registered ones, and `labels/` holds the classifier evidence for the validation answers (step 3, 2026-10-28).
- `experiments/harness/`: `TRIES.md` (the tries ledger and the headroom tables per test board), `headroom.py` (writes
  those tables' numbers to `results/headroom/` from committed results), `run.py` (harness runs and adoption decisions).
- `boards/`: public board specs from the challenge site (`index.json`, panel gene lists, T1 composition). Since 634b5fd
  `index.json` also holds the test boards' cell limits: code reads it; nobody opens it or quotes a test board's limits
  (`docs/DATA_POLICY.md`).
- `data/` (gitignored): `raw/` competition files, `external/` disclosed external data.
- `docs/`: `CHALLENGE.md`, `PLAN.md`, `DATA_POLICY.md`, `COMPUTE.md`, `research/` (background reports).
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

1. Update `progress.md` with the current state.
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
