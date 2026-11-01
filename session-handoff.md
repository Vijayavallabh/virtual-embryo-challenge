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
  still answers 403 (step 1 blocked; `progress.md`, Blockers; still so on 2026-10-31). At the user's request its
  saved-file intake is built (`vec boards update --from-file`, cee1c92) but not cleared for use: saving the file in a
  browser would work around the 403, which only the user can allow. Step 4, the pre-registered runs,
  ran on 2026-10-28 at the user's go-ahead, each output committed before it was read.
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
| Lint and tests | `./init.sh` | passes (2026-10-28) | `uv run pytest -q` 402 tests (90+ on 2026-10-19); experiment tests in `experiments/*/test_*.py` run separately |
| Pre-registered runs (step 4) | `docs/P3_RUNBOOK.md` step 4 commands (DGX, 2026-10-28) | all six ran, exit 0 | outputs committed before reading: f5be247 to 12151e0; the commits are listed in the runbook's step 4 status |
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

- 2026-11-01 (latest): the s1 build path was reviewed (four reviewers, skeptic-verified; 11 findings) and fixed in
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
