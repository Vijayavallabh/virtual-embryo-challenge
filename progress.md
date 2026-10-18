# Session Progress Log

## Current State

**Last Updated:** 2026-10-18
**Active Feature:** feat-002 - Local board emulator (calibration pending real data)

## Status

### What's Done

- [x] Read the full challenge site, rules, data page, proposal and the veckit scorer code; condensed in `docs/CHALLENGE.md`.
- [x] uv project, harness (`AGENTS.md`, `feature_list.json`, `init.sh`, `progress.md`, `session-handoff.md`).
- [x] `vec` package: board specs, local board emulator on veckit's metric code, format checks, reference baselines.
- [x] Data policy and compute docs; external-data register; method summary skeleton.

### What's In Progress

- [ ] Research sweep across methods, external data, metrics, biology; plan in `docs/PLAN.md`.

### What's Next

1. Register the team and download data to the DGX (feat-003).
2. Reproduce the published reference rows locally to calibrate the emulator (feat-002).

## Blockers / Risks

- [ ] Data download needs the user's registration at virtualembryo.ai.
- [ ] The T2 shape-group gate formula is not published; the emulator assumes `floor_normalised` until calibrated.

## Decisions Made

- **Human Team track only** (user decision, 2026-10-18).
- **Python 3.11, uv only; CUDA torch from the cu128 index on Linux, synced only with `--group gpu`.**
- **veckit pinned at cc3859a** so local scores are reproducible; the emulator calls its metric functions directly.

## Evidence of Completion

- [x] `./init.sh`: lock check, ruff clean, 11 tests pass (2026-10-18).

## Notes for Next Session

- `neighborhood_mmd` builds 15-NN neighbourhoods over all submitted cells, so cell count (density) moves T2 scores
  even though count is not ranked. Tune the submitted cell count on the validation answers.
- Ctnnb1 is not on the 500-gene MERFISH panel: gene_ko gives nothing for the beta-catenin test board.
