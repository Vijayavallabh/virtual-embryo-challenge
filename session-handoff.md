# Session Handoff

## Current Objective

- Goal: a compliant Human-track entry on all four test boards; the final phase opened on 2026-10-27 at 16:00 UTC,
  postponed from 2026-10-20 by the organisers (2 scored uploads per board for the whole phase; none used as of
  2026-10-28).
- Current status: every method choice that must precede the release is committed (size procedure, T3 decomposition
  rules, T1 direction readout, count rule, label-mapping rule, test-anchor quarantine). The user signed off the
  release-day rules on 2026-10-21 (`docs/RELEASE_RULES_SUMMARY.md`), and the same day the wording for the case the
  count rule left open, whose reading of "next" the user left to the AI assistant (a decision panel chose the next
  passing count up); the heart-size fallback and the count rule run as code. Release-day steps are in
  `docs/P3_RUNBOOK.md`.
  The phase opened on 2026-10-27. On 2026-10-28 the user unlocked the validation answers on the Data page (P2
  submissions now permanently disabled for the team); all six are on the DGX, and runbook steps 2 (data) and 3
  (labels) are done. The user read the final-phase Rules s.10 and Data page as allowing the validation answers as
  model inputs (A20 passed; confirmation asked in the fourth organiser email, sent 2026-10-28). The board-spec file
  still answers 403 (step 1 blocked; `progress.md`, Blockers).
- Branch / commit: `main` on GitHub (private), mirrored to `PrakashDGX_H2:~/v/virtual-embryo-challenge`.

## Completed This Session (2026-10-27 to 2026-10-28)

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
| Lint and tests | `./init.sh` | passes (2026-10-28) | `uv run pytest -q` 366 tests (90+ on 2026-10-19); experiment tests in `experiments/*/test_*.py` run separately |
| Size procedure | `uv run pytest experiments/coords -q -p no:cacheprovider` | 33 passed (2026-10-19) | preflight requires committed procedure files |
| QC of released files | `bash experiments/qc/run_qc.sh <tag>` (DGX) | `experiments/qc/results/` | duplicate and interpolated counts in `docs/DATA_NOTES.md` |
| Validation answers (step 2) | `uv run vec data inventory`, `uv run vec data qc` (DGX, 2026-10-28) | 6 files, 2.09 GB, sizes as listed | SHA-256 of all 15 files in `data/raw/SHA256SUMS`; no duplicated or interpolated records; `docs/DATA_NOTES.md` (f423f6a) |
| Labels (step 3) | `uv run python experiments/labels/release_label_evidence.py` (DGX) | label evidence in `experiments/labels/results/` | every label of the 15 files maps (DGX check, 2026-10-28); the new entries pinned in `tests/test_prep.py` and `tests/test_t1_scope.py` (d6b9566) |

## Decisions Made

- See `docs/PLAN.md` decision log (2026-10-18, 2026-10-19 a-d, 2026-10-20 a-c, 2026-10-21 a-c, 2026-10-28) and
  `docs/DATA_POLICY.md` records. The user's calls of 2026-10-27 and 2026-10-28 (the history rewrite, the unlock, the
  A20 reading, the go-ahead for step 3) are recorded in `progress.md` (Decisions Made, Blockers / Risks).

## Blockers / Risks

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
2. Read `feature_list.json` and `progress.md`.
3. Review this handoff and `docs/P3_RUNBOOK.md`.
4. Run `./init.sh` before editing.

## Recommended Next Step

- Continue `docs/P3_RUNBOOK.md` at step 4 (the pre-registered runs, in order), committing each output before anything
  else; ask the user before starting. Step 1 waits for the panel files (403). On the DGX, cap the BLAS threads
  (`docs/COMPUTE.md`; `progress.md`, Notes for Next Session).
