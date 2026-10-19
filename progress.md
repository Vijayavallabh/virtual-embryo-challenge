# Session Progress Log

## Current State

**Last Updated:** 2026-10-19
**Active Feature:** feat-004 - Final-phase intake on 2026-10-20 (`docs/P3_RUNBOOK.md`); feat-003 completes with the release files

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
- [x] 2026-10-19: submission client and builders, each built by an agent and verified by an independent agent, then
      three rounds of independent review and fixes (every finding verified fixed):
      - `vec submit` (feat-016): format-only by default; official uploads behind --official, a hash confirmation and
        every gate in `src/vec/submit.py`; `vec submit status` is read-only.
      - `vec.t2_build`: count rule with the truth floor enforced, random rows, jitter, size procedure, cloud guard;
        smoke test on both validation boards (`experiments/t2_build/results/`).
      - `vec.t1_scope`: chamber categories for every released T1 label (`configs/t1_chamber_classes.yaml`),
        re-weighting to input shares, momentum shift; E8.5 -> E9.5 pseudo-task (`experiments/t1_scope/results/`).
      - `experiments/t3_decomp/build_submission.py`: s1/s2 from the decision file, no re-deciding.

### What's Next

1. 2026-10-20: follow `docs/P3_RUNBOOK.md` (board intake with anchors stripped, `vec submit status` to read the quota
   shape, data, label map, pre-registered runs, commit each output), then emulator calibration (feat-005) and
   test-board intake (D1).
2. B0 floor-plus files on all four boards with the builders (feat-009); T3 per-gene factor on Gata4 (feat-012); T1
   scope shares by a written rule (feat-014).

## Blockers / Risks

- [ ] Awaiting organiser replies: leak report (Slack DM, sent 2026-10-18), the rules questions (email, 2026-10-18) and
      the follow-up (scorer build and parity defect, wrong-sign severity, released Gata4 in the T3 prediction,
      duplicated/interpolated records, z units, truth_cells and size, stage conventions; sent 2026-10-19).
- [ ] The E12.5 emergent-class list waits on the organisers' answer about qualitative in-window background; the legal
      route is used meanwhile.
- [ ] Quota shape: `official_per_board` is unpublished. The client refuses official uploads until it can read it;
      check with `vec submit status` once the final phase opens and adapt the parser (with a review) if needed.
- [ ] Upload logistics (user's decision, before the first official upload in November): the key is local-only and
      built files live on the DGX, so a file must reach this machine to be uploaded (server format checks included).
- [ ] Provenance across machines (before the first official upload): the client's clean-tree gate sees only the local
      clone, while files are built on the DGX. Planned: the DGX writes a stamp next to each built file (commit, clean
      flag, sha256), and an official upload is refused unless the stamp's commit is the local HEAD and the tree was
      clean.
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
- **GPUs:** free VRAM only until the end of October; all eight H100s from November (user, 2026-10-19).

## Evidence of Completion

- [x] `./init.sh`: lock check, ruff clean, 237 tests pass; experiment tests 203 pass, 2 skipped locally (the T3
      end-to-end tests, which pass on the DGX: 50 T3 tests in all).
- [x] DGX: torch 2.14.1+cu130 with CUDA, jax 0.11.2, moscot 0.5.2, scvi-tools 1.5.1 import (runs/logs/setup2.log).

## Notes for Next Session

- `neighborhood_mmd` builds 15-NN neighbourhoods over all submitted cells, so submitted density moves T2 scores.
- Ctnnb1 is not on the 500-gene MERFISH panel: gene_ko gives nothing for the beta-catenin test board.
- T1 E12.5 is a heart-chamber-only dissection; earlier stages include surrounding tissue.
- DGX research scratch removed on 2026-10-19 (user's OK); `~/v/vec-gapconc` keeps only `ext/` and `out/`, the
  legal inputs the committed filters in `experiments/external/` read (the JAX files are hard-linked into data/external).
- `experiments/harmonise_review/review.py` has 5 pre-existing ruff findings (outside the checked `src tests`).
