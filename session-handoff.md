# Session Handoff

## Current Objective

- Goal: a compliant Human-track entry on all four test boards; the final phase opens 2026-10-27, postponed from
  2026-10-20 by the organisers (2 scored uploads per board for the whole phase).
- Current status: every method choice that must precede the release is committed (size procedure, T3 decomposition
  rules, T1 direction readout, count rule, label-mapping rule, test-anchor quarantine). The user signed off the
  release-day rules on 2026-10-21 (`docs/RELEASE_RULES_SUMMARY.md`), and the same day the wording for the case the
  count rule left open (how the code reads "next" awaits the user's confirmation; Blockers); the heart-size fallback
  and the count rule run as code. Release-day steps are in `docs/P3_RUNBOOK.md`.
- Branch / commit: `main` on GitHub (private), mirrored to `PrakashDGX_H2:~/v/virtual-embryo-challenge`.

## Completed This Session (2026-10-18 to 2026-10-19)

- [x] Harness, uv environment, `vec` package (emulator, checks, baselines, safe downloader, loaders, QC, harmonise).
- [x] Training data on the DGX; research sweep, screened digest, PLAN v1.2, assumption ledger (19 rows).
- [x] Experiments with independent reviews: density, coordinates/size procedure, harmonisation, data QC.
- [x] Compliance: leaked heart size quarantined and reported; gene-level knockout findings and held-out-stage values
      redacted by meaning; leaked-scorer copies and in-window external files deleted; history rewritten and the GitHub
      repo recreated (first pass; a second pass is ready for the user to run).

## Verification Evidence

| Check | Command | Result | Notes |
|---|---|---|---|
| Lint and tests | `./init.sh` | passes | `uv run pytest -q` 90+ tests; experiment tests in `experiments/*/test_*.py` run separately |
| Size procedure | `uv run pytest experiments/coords -q -p no:cacheprovider` | 33 passed | preflight requires committed procedure files |
| QC of released files | `bash experiments/qc/run_qc.sh <tag>` (DGX) | `experiments/qc/results/` | duplicate and interpolated counts in `docs/DATA_NOTES.md` |

## Decisions Made

- See `docs/PLAN.md` decision log (2026-10-18, 2026-10-19 a-d, 2026-10-20 a-c, 2026-10-21 a-b) and
  `docs/DATA_POLICY.md` records.

## Blockers / Risks

- Organiser answers pending (leak report; rules email). Follow-up email sent 2026-10-19; a third (three questions)
  sent 2026-10-21 01:24 IST (`docs/ORGANIZER_QUESTIONS.md`). Deadline: through 2026-12-02, UTC.
- The user's confirmation of how the count floor rule reads "next" (rule 7's order, coded and recommended, or the
  next count up; `progress.md` Blockers), before 2026-10-27.
- Collaborative Yukon track (organisers' email of 2026-10-20 22:03 IST): joining is the user's decision alone; the
  entry stays an official Human Team meanwhile.
- Second history rewrite done 2026-10-19 (GitHub repo recreated; old-to-new commit map in `docs/DATA_POLICY.md`).
- Never open the public board file, veckit 0.1.2, PR #6, raw research outputs or Yukon team material (AGENTS.md,
  quarantined material).

## Next Session Startup

1. Read `AGENTS.md`.
2. Read `feature_list.json` and `progress.md`.
3. Review this handoff and `docs/P3_RUNBOOK.md`.
4. Run `./init.sh` before editing.

## Recommended Next Step

- On 2026-10-27, once `vec submit status` shows the final phase, run `docs/P3_RUNBOOK.md` top to bottom, committing each
  pre-registered output before anything else. Until then, only work that needs no release data (`docs/PLAN.md`
  section 7).
