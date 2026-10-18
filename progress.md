# Session Progress Log

## Current State

**Last Updated:** 2026-10-18
**Active Feature:** feat-003 - Data acquisition (blocked on the user's session token)

## Status

### What's Done

- [x] Read the full challenge site, rules, data page, proposal and the veckit scorer; condensed in `docs/CHALLENGE.md`.
- [x] uv project (Python 3.12), harness, `vec` package (emulator, checks, baselines, safe downloader), DGX environment.
- [x] Research sweep (10 reports in `docs/research/`), PantheonOS Slack read in full, plan in `docs/PLAN.md`.
- [x] `ASSUMPTIONS.md` (18 rows, linted), data policy with organiser rulings and a do-not-use list.
- [x] Leaked E12.5 heart size quarantined: redacted from notes; private report to organisers drafted in Slack.

### What's Next

1. User saves the session token to `~/.config/vec/token`; then `uv run vec data fetch` (feat-003).
2. 2026-10-20 final-phase intake: new index.json, Data Protocol, scorer version (feat-004), then emulator calibration
   against the published reference rows (feat-005).
3. T3 day-1 diagnostics and the T2 heart size procedure (feat-012, feat-011).

## Blockers / Risks

- [ ] Data download needs the user's session token.
- [ ] The leak report to the organisers is a draft until the user sends it (Slack group DM `C0CA7QTEMDF`).
- [ ] Organiser questions (cross-task data, test-board definitions) not sent yet.

## Decisions Made

- **Human Team track only** (user, 2026-10-18).
- **Python 3.12; torch 2.14 cu130; veckit pinned to git main cc3859a.**
- **No mouse single-cell foundation models** (in-window contamination). **No gene-level T3 knockout literature.**
- **Heart size only from the pre-registered training-stage procedure** (docs/PLAN.md section 4.3).

## Evidence of Completion

- [x] `./init.sh`: lock check, ruff clean, 18 tests pass.
- [x] DGX: torch 2.14.1+cu130 with CUDA, jax 0.11.2, moscot 0.5.2, scvi-tools 1.5.1 import (runs/logs/setup2.log).

## Notes for Next Session

- `neighborhood_mmd` builds 15-NN neighbourhoods over all submitted cells, so submitted density moves T2 scores.
- Ctnnb1 is not on the 500-gene MERFISH panel: gene_ko gives nothing for the beta-catenin test board.
- T1 E12.5 is a heart-chamber-only dissection; earlier stages include surrounding tissue.
- Agent scratch dirs on the DGX (`~/v/vec-t1-scratch`, `~/v/vec2026-env-probe`, ~7.5 GB) can be removed once no
  longer needed.
