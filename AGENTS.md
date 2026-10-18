# AGENTS.md

Team repository for the **Virtual Embryo Challenge** (NeurIPS 2026, https://virtualembryo.ai/challenge), entered in
**Track 1 (Human Team)**. We predict held-out embryo stages: T1 single-cell RNA (E12.5), T2 3D MERFISH expression plus
coordinates (heart E12.5, embryo E7.75), T3 a held-out knockout (beta-catenin, E8.75).

## Startup Workflow

Before writing code:

1. Confirm the working directory with `pwd`.
2. Read this file completely.
3. Read `docs/CHALLENGE.md` (tasks, metrics, weights, rules), `docs/PLAN.md` (strategy and milestones),
   `docs/DATA_POLICY.md`, and `docs/COMPUTE.md`.
4. Run `./init.sh` to verify the environment is healthy.
5. Read `feature_list.json` for the current feature state and `progress.md` for the last session.
6. Review recent commits with `git log --oneline -5`.

If baseline verification is failing, repair that first before adding new scope.

## Invariants (never violate)

- **Data-use rules.** Never train, tune or select on measured data inside a held-out window, by any route
  (directly, via a pretrained model, via a public atlas). The windows are listed in `docs/DATA_POLICY.md`. Every
  external dataset, pretrained model or published code must be logged in `EXTERNAL_DATA.md` *before* it enters a model.
- **Scores choose, never compute.** Leaderboard returns may be used to pick between predictions, never to solve for a
  property of a hidden target (no inverting the skill map, no differencing scores).
- **Competition data stays private.** Data lives under `data/` (gitignored) and is never committed, shared or shown
  publicly. The dataset is unpublished (Neil Chi lab, UCSD): no public presentation without written approval.
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
- **Evidence over claims**: a model result counts only with a local-emulator score file (`runs/<run_id>/scores.json`)
  and a row in `results/registry.csv`.
- **Update artifacts**: before ending a session, update `progress.md` and `feature_list.json`.
- **Stay in scope**: don't modify files unrelated to the current feature.
- **Leave clean state**: the next session must be able to run `./init.sh` immediately.

## Repository Map

- `src/vec/`: package (`io.py` loading and submission writing, `boards.py` board specs, `evaluate.py` server-faithful
  local board emulator, `baselines.py` reference baselines, `cli.py` the `vec` command).
- `configs/`: experiment configs. `runs/` (gitignored): per-run outputs. `results/registry.csv`: one row per scored run.
- `boards/`: public board specs from the challenge site (`index.json`, panel gene lists, T1 composition).
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

## Escalation

- **Rule ambiguity** (external data near a window, cross-task data use): stop and ask the user; draft the question
  for the organisers (virtual.embryo.moonshot@gmail.com or Slack #virtual-embryo-challenge).
- **Architecture decisions**: consult `docs/PLAN.md`, otherwise ask the user.
- **Repeated test failures**: update `progress.md` and flag for human review.
- **Anything touching an official submission**: ask the user first.
