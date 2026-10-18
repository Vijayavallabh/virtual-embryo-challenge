# Session Progress Log

## Current State

**Last Updated:** 2026-10-19
**Active Feature:** feat-004 - Final-phase intake on 2026-10-20 (`docs/P3_RUNBOOK.md`); feat-003 completes with the release files

## Status

### What's Done

- [x] Read the full challenge site, rules, data page, proposal and the veckit scorer; condensed in `docs/CHALLENGE.md`.
- [x] uv project (Python 3.12), harness, `vec` package (emulator, checks, baselines, safe downloader), DGX environment.
- [x] Research sweep (10 reports in `docs/research/`), PantheonOS Slack read in full, plan in `docs/PLAN.md`.
- [x] `ASSUMPTIONS.md` (18 rows, linted), data policy with organiser rulings and a do-not-use list.
- [x] Leaked E12.5 heart size quarantined: redacted from notes; private report to organisers drafted in Slack.
- [x] Team API key stored at ~/.config/vec/token (mode 600, local only); 9 training files downloaded to the DGX and checksummed.
- [x] 2026-10-19: experiments on released data, each independently reviewed (`experiments/*/REPORT.md`): density and
      sampling (count rule max_cells, never below truth_cells; random rows only), coordinates and the pre-registered heart
      and embryo size procedure (committed 171695a before the release), label harmonisation (`vec.harmonise`), data QC
      (`vec.spatial`, `vec data qc`; copied and interpolated records dropped from prediction inputs).
- [x] 2026-10-19: research sweep screened by a separate agent; PLAN v1.2 built from the screened digest only
      (`docs/research/digest_2026-10-19.md`). Two meaning-level redaction passes on the notes; leaked-scorer copies,
      anchor-derived outputs and in-window external files deleted (user's decision).
- [x] 2026-10-19: test-board anchors quarantined at intake (`vec.boards.update_boards`); release-day runbook written.

### What's Next

1. 2026-10-20: follow `docs/P3_RUNBOOK.md` (board intake with anchors stripped, data, label map, pre-registered runs,
   commit each output), then emulator calibration (feat-005) and test-board intake (D1).
2. User: run the history rewrite (`rewrite_history2.sh` in the session scratchpad; needs the gh `delete_repo` scope) and
   send the organiser follow-up in `docs/ORGANIZER_QUESTIONS.md`.
3. Then B0 floor-plus on all four boards (feat-009), T3 per-gene factor (feat-012), T1 scope re-weighting (feat-014).

## Blockers / Risks

- [ ] Awaiting organiser replies: leak report (Slack DM `C0CA7QTEMDF`, sent 2026-10-18) and the rules questions (email,
      sent 2026-10-18). The follow-up (scorer build and parity defect, wrong-sign severity, released Gata4 in the T3
      prediction, duplicated/interpolated records, z units, truth_cells and size, stage conventions) is ready, not sent.
- [ ] Old commits still hold redacted values until the user runs the history rewrite.
- [ ] The E12.5 emergent-class list waits on the organisers' answer about qualitative in-window background; the legal
      route is used meanwhile.

## Decisions Made

- **Human Team track only** (user, 2026-10-18).
- **Python 3.12; torch 2.14 cu130; veckit pinned to git main cc3859a.**
- **No mouse single-cell foundation models** (in-window contamination). **No gene-level T3 knockout literature.**
- **Heart size only from the pre-registered training-stage procedure** (docs/PLAN.md section 4.3).
- **Count = max_cells, never below truth_cells; random rows only** (A13 revised, 2026-10-19).
- **Heart size procedure fixed before the release** (core = CM + endocardium; E9.5-E10.5 only; fixed three-rule average;
  z held when the plane matches; H1/H2 scope). **T3 leads with the per-gene factor** from released knockouts.
- **No parity hedge without the organisers' approval** (user, 2026-10-19). **History rewritten again** (user).

## Evidence of Completion

- [x] `./init.sh`: lock check, ruff clean, 18 tests pass.
- [x] DGX: torch 2.14.1+cu130 with CUDA, jax 0.11.2, moscot 0.5.2, scvi-tools 1.5.1 import (runs/logs/setup2.log).

## Notes for Next Session

- `neighborhood_mmd` builds 15-NN neighbourhoods over all submitted cells, so submitted density moves T2 scores.
- Ctnnb1 is not on the 500-gene MERFISH panel: gene_ko gives nothing for the beta-catenin test board.
- T1 E12.5 is a heart-chamber-only dissection; earlier stages include surrounding tissue.
- Agent scratch dirs on the DGX (`~/v/vec-t1-scratch`, `~/v/vec2026-env-probe`, ~7.5 GB) can be removed once no
  longer needed.
