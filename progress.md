# Session Progress Log

## Current State

**Last Updated:** 2026-10-28
**Active Feature:** feat-004 - Final-phase intake (`docs/P3_RUNBOOK.md`). Step 2 (data) done on 2026-10-28 (feat-003
done); step 3 (labels) is next; step 1 is blocked by the 403 on the panel files; steps 4.1-4.3 wait for the user's
reading of A20. The release-day rules are signed off and in code (`docs/RELEASE_RULES_SUMMARY.md`).

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
- [x] 2026-10-20: the organisers moved the final phase to 2026-10-27 (email). The two code paths that could read
      held-out validation material now wait for the final phase itself (the T1 composition at `vec boards update`; the
      anchor check); runbook, plan timeline, data policy, ledger (A20) and handoff updated; Yukon team material
      quarantined. Slack and the organiser email swept by screened agents (findings under Blockers).
- [x] 2026-10-20: the upload path tested against the live server while the validation phase runs a week longer: three
      format checks pulled from the DGX (T1:val, T2:heart:val_extrap, T3:gata4; practice files from
      `vec baseline --n-cells 1000`) each came back format_ok with no submission id; before and after, no server entry
      on any board and the daily scored and upload counts at 0 (feat-016 evidence; ledger in results/submissions.csv).
- [x] 2026-10-20 (organisers' email of 22:03 IST): a collaborative Yukon track (Eigen Labs; invite-only; separate
      $10K-$30K pool; collaborative phase ends before P3; one-way transfer; ineligible for the official prizes) and the
      deadline's time zone (through 2 December 2026, UTC). An agent read the official Rules, Terms, Privacy,
      partner-rules, Timeline and FAQ pages (nothing quarantined opened): ordinary-team rules are unchanged; official
      P3 scores are published when each submission is made; the P3 privacy rule (methods, prediction files and
      evaluation feedback private until 2026-12-11) is on the partner page only, and we follow it. Challenge
      reference, data policy, plan, organiser-questions log and AGENTS.md updated; the organiser draft cut to two
      questions.
- [x] 2026-10-20: prize strategy, approved by the user ("do whatever's best"; `docs/PLAN.md` decision log 2026-10-20
      (c)). Stay official; go for a Human-track place plus the Generality Award, which stacks on a place and is judged
      across both tracks among shared-architecture entries. D3 loosened to "claim unless specialised code wins by more
      than 3 points overall"; the shared core goes into code first (feat-020); ledger A16 revised, A21 added. The
      Generality question was added to the organiser draft (now three questions). Prize terms and a leaderboard census
      (323 Human teams, 159 above the floor, 6th place 223.5; Agent track 80) recorded in `docs/CHALLENGE.md` and PLAN
      section 1. Not pursued: Yukon; a switch to the Agent track.
- [x] 2026-10-21: feat-020, the shared transport-then-adjust core `src/vec/core.py` (Generality Award claim). The
      T1 and T2 heads' allocation, draws and count-space edits moved into it, source unchanged except names; the
      pre-registered T3 script keeps its own primitives and the core holds identical copies. Evidence
      (`experiments/core/REPORT.md`):
      - `tests/test_core.py`: 31 tests, including 12 golden cases recorded at eb946ac and T3 copies identical in
        body and on 300 random inputs per function;
      - released data on the DGX: old against new, all 64 hashed parts identical;
      - `./init.sh`: 320 passed;
      - an independent verifier: PASS. Its fixes were made: the module was renamed from `vec.operator`, which
        shadowed the standard library, and the wording about T3 was corrected.
      Release-day rules summary drafted for sign-off (`docs/RELEASE_RULES_SUMMARY.md`, with recommended answers).
- [x] 2026-10-21: the two approved answers that needed code (answers 3 and 7 of `docs/RELEASE_RULES_SUMMARY.md`),
      built by an agent and checked by three independent verifiers (fidelity, integrity, integration) over two fix
      rounds; every code finding was fixed:
      - `vec.t2_build.size_or_hold`, the heart size step: a valid JSON is applied exactly as `apply_size_procedure`
        applies it; `valid` false keeps E10.5's size; cells from another stage get no size (back to the user); every
        heart call's reference cloud must be the cloud the JSON measured as E10.5.
      - `count_override`, `read_count_decision` and `n_cells_from_decision` in `vec.t2_build`, and
        `experiments/density/summarise.py --decide`: the count decided per board from the release rerun, written once
        to `experiments/density/results/release/count_decision.json`; uncommitted rule code and a second decision
        without `--overwrite` are refused.
      - `docs/P3_RUNBOOK.md` (steps 4.1 and 4.5, section 5) and `run_p3.sh`'s closing note now name them.
        `./init.sh`: 355 passed (35 new tests).
- [x] 2026-10-21: the case approved answer 7 left open, in code as the user approved it ("the next passing count the
      board allows, or max_cells if none does"; read as the next count up, the reading the user left to the AI
      assistant and a decision panel chose, `experiments/density/floor_reading/README.md`):
      - `vec.t2_build.n_cells_from_decision` passes over a decided level below the test board's lower limit to the
        passing level with the fewest cells that the board allows, else max_cells, and records `count_floor`
        (`COUNT_FLOOR_RULE`);
      - `read_count_decision` now also checks the decision's record: the decided level heads its passing levels,
        every level's verdict agrees with its own numbers (both ways), and the rerun was complete;
      - runbook section 5, the rules summary, PLAN 3.6 and 4.2 step 6, the decision log, METHOD_SUMMARY and the
        ledger log updated; the third organiser email recorded as sent (`docs/ORGANIZER_QUESTIONS.md`).
      Independent verifiers in two rounds on the gain reading (fidelity, robustness and docs; then code and docs
      re-checks): every finding fixed; `./init.sh`: 357 passed (3 new tests, 1 replaced; d06727f). The reading of
      "next" was then put to the user, who left it to the AI assistant; a four-agent panel (simulation, data, an
      advocate for the gain reading, judge) chose the next count up, medium confidence. The switch to it: a code
      verifier (pass; 0 mismatches in 579,150 cases against an independent implementation) and a docs verifier, their
      findings fixed; `./init.sh`: 357 passed.
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

1. Until 2026-10-27: work that needs no release data (`docs/PLAN.md` section 7): external data through committed
   window filters (feat-007), v1 pipelines on released stages, B0 builds checked with `vec check` on the DGX.
2. 2026-10-27, once `vec submit status` shows the final phase: follow `docs/P3_RUNBOOK.md` (board intake with anchors
   stripped, the quota shape, data, label map, pre-registered runs, commit each output), then emulator calibration
   (feat-005) and test-board intake (D1).
3. B0 floor-plus files on all four boards with the builders (feat-009); T3 per-gene factor on Gata4 (feat-012); T1
   scope shares by a written rule (feat-014).

## Blockers / Risks

- [ ] Awaiting organiser replies: leak report (Slack DM, sent 2026-10-18), the rules questions (email, 2026-10-18) and
      the follow-up (scorer build and parity defect, wrong-sign severity, released Gata4 in the T3 prediction,
      duplicated/interpolated records, z units, truth_cells and size, stage conventions; sent 2026-10-19). Screened
      sweeps of Slack and the organiser email on 2026-10-20 (to 22:43 IST): no organiser post on Slack since
      2026-10-16, no reply in any DM or by email (the sender's mail tracker reported both emails unopened 24 hours
      after sending); the platform emails are the postponement (07:41 IST) and the collaborative Yukon track (22:03
      IST), which gave the deadline's time zone (through 2 December 2026, UTC). Unanswered in public: qualitative
      in-window priors and (2026-10-19) whether released validation answers may be used for test predictions (A20).
      Final uploads still finish by the end of 30 November IST; keep the strict in-window reading. The third email
      (three questions: the test phase's start time, validation answers as training data (A20), the Generality scope
      (A21)) was sent by the user on 2026-10-21 at 01:24 IST (`docs/ORGANIZER_QUESTIONS.md`). Not checked by agents:
      the challenge account messages page (user only) and Gmail Spam.
      Screened Slack sweep on 2026-10-23 (from 2026-10-20 22:43 IST): still no organiser post or reply anywhere;
      two new unanswered participant questions (2026-10-21: textbook E9.5-E12.5 knowledge as a prior and whether to
      disclose it; how travel grants are chosen).
- [ ] Release day (2026-10-27, checked to 17:43 UTC): the final phase is open (public phase endpoint: p3, split
      test, accepts submissions; official quota keys T1:test, T2:heart:test_extrap, T2:embryo:test_interp and
      T3:bcatenin, 0/2 used each; daily per task: scored 0/2, uploads 0/8). Runbook step 0 passed (357 tests, nothing
      unresolved). Three things stop the intake: (1) the data manifest was refused (HTTP 409) until the user accepted
      the current Challenge Rules and Terms of Use on the account page (done about 17:42 UTC; an API key cannot
      accept them); (2) it then lists only the 9 training files, with the note "Training inputs only. Validation and
      test targets are not distributed - a submission comes back as a score." (first recorded today; unknown whether
      it is new or a leftover from the validation phase, and it conflicts with the Data page's promise that the
      validation answers become training data); (3) `vec boards update` fails: panels/index.json answers HTTP 403
      from Cloudflare to our client, so nothing was written (not worked around). The screened read of the new Rules
      and Terms was stopped by the user before it reported: what changed in them is not yet known. Nothing was
      downloaded, built or uploaded. The user asked to try again on 2026-10-28.
      Retry 2026-10-28 04:37 UTC: phase, quota and manifest unchanged (9 training files, same note); `vec boards
      update` still 403, nothing written. A screened read of the site (04:40-04:49 UTC; quotes came through WebFetch,
      so check the live wording before relying on them) explains the missing files: the validation answers are now
      opt-in. Data page and Rules s.10: they are "available to eligible registered participants during P3, after
      explicit confirmation", and the opt-in files are "a specific exception to the held-out-data restrictions ...
      solely for preparing P3 predictions after access is granted" (no external E10.5/E12.5 data; no test truth).
      Rules s.11: obtaining access (even without downloading) or the first accepted P3 submission "permanently ends
      P2 submissions across all tasks for the entire ordinary team"; viewing the Data page or a format check does not.
      No test-input files will come ("Hidden test ground truth is never distributed"). P3 opened 2026-10-27 16:00
      UTC (answers question 1 of the third email). The Data page still links the panel files under /challenge/panels/
      and announces no outage, so the 403 is unexplained. No page says how the opt-in works, who is eligible, or
      whether a Human Team must nominate entries (the phase endpoint shows nominations_open). Waiting on the user:
      the opt-in (irreversible, team-wide, browser only) and the reading of "solely for preparing P3 predictions"
      for inputs at prediction time (A20 stays open until then).
      2026-10-28 09:32 UTC: the user unlocked the validation answers on the Data page (server state: granted
      09:32:18 UTC, lock_scope team, so P2 submissions are permanently disabled for the whole team; the page says it
      uses no P3 attempt). They have their own list and links (`/challenge/data/validation/` state, manifest and
      link, found in the site's public script); the unlock is a separate POST that no repo code sends
      (`tests/test_data.py` checks for it). The data client reads the validation list since e0e1d58, and all 6 files
      (2.09 GB) are on the DGX: runbook step 2 done (checksums, inventory and QC; `docs/DATA_NOTES.md`). The Data
      page's own words on use: "P2 validation answers can be used to prepare your final P3 models" and "After you
      explicitly unlock the official P2 validation release above, you may use those released validation files for
      P3 preparation. This does not permit access to hidden test data or other restricted external sources."
      Still open: the 403 on the panel files (step 1, the test-board specs) and the user's reading for A20 (inputs
      at prediction time), which steps 4.1-4.3 need.
- [ ] Final-entry nomination (found 2026-10-28 in the site's public account-page script): "In the final phase the
      captain picks which submission counts" (a "Final" column in the captain's submission list while
      nominations_open); the timeline makes 2026-12-02 the last moment to nominate. Rules s.11 still says the board
      ranks a team by its best; what counts if the captain marks nothing is not stated. The user (captain) marks
      the final entries after the last official upload; ask the organisers what an unmarked board gets
      (`docs/P3_RUNBOOK.md`, section 5).
- [x] 2026-10-27: history of `main` rewritten (user's decision) to drop the AI co-author trailer from all 44 commit
      messages; trees, parents, authors and dates checked identical; old -> new IDs in `docs/COMMIT_MAP.md`;
      force-pushed, DGX .git refreshed. Older IDs quoted in files resolve through that table.
- [ ] The user's decision: the collaborative Yukon track (organisers' email of 2026-10-20 22:03 IST). The entry stays
      an official Human Team unless the user decides otherwise; nothing has been done toward joining. Moving is
      one-way and ends official-prize eligibility, a team left empty is archived, and the collaborative phase ends
      before 2026-10-27.
- [x] 2026-10-21: the user signed off ("approve all") on the other agent-written rules that act on the
      release (heart and embryo size procedure, T1 direction readout, cell-count rule, label-mapping rule): all eight
      recommended answers in `docs/RELEASE_RULES_SUMMARY.md` (drafted by the AI assistant on 2026-10-20). The gaps it
      found are closed:
      - the T1 readout's second role is written down (`experiments/t1_direction/README.md`, PLAN 4.1);
      - poor harmonisation groups are held flat (PLAN 4.3);
      - the fallback for an invalid heart size and the count rule are now in code (`vec.t2_build` and
        `summarise.py --decide`; see What's Done).
- [x] 2026-10-21: the user's call on the case approved answer 7 left open, a winning count level below a test
      board's lower limit (truth_cells or min_cells; only 0.25 x max_cells can fall there): the user approved ("go
      ahead") the recommendation "the next passing count the board allows, or max_cells if none does". In code as
      `vec.t2_build.COUNT_FLOOR_RULE` (see What's Done).
- [x] 2026-10-21: how the code reads "next" in that rule. The user left it to the AI assistant ("Anything that would
      improve my chances of winning"). A decision panel chose the next count up over the next best gain: they differ
      only when 0.25 x max_cells wins below the limit while 0.5 x and 0.75 x both pass and 0.75 x measures at least as
      high, and there the next count up is ahead by about 0.2-0.4 points in most simulated models (a data analysis put
      it near zero; `experiments/density/floor_reading/README.md`). In code (see What's Done).
- [ ] The E12.5 emergent-class list waits on the organisers' answer about qualitative in-window background; the legal
      route is used meanwhile.
- [ ] Quota shape: `official_per_board` is unpublished. The client refuses official uploads until it can read it;
      check with `vec submit status` once the final phase opens and adapt the parser (with a review) if needed.
- [x] Upload logistics (decided by the user, 2026-10-19; built and verified the same day): `vec submit PATH --board KEY
      --from-dgx` pulls the built file into a private copy for that call only; the key never goes to the DGX.
- [x] Provenance across machines: "dgx provenance" alone checked only the code checked out at upload time, not the
      code that built the file (a side agent caught this; I had marked it closed too early). Fixed 2026-10-20 with
      build stamps: every prediction file gets `<file>.stamp.json` when it is built (sha256, the commit, whether the
      tree was clean from the build's start to the write), and an official upload needs a stamp that matches the
      bytes sent, says clean and names a commit in our history; the ledger's build_commit records it.
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
- **Uploads pull from the DGX** (user, 2026-10-19): `vec submit PATH --board KEY --from-dgx` streams the built file
  from the DGX into the client's private temporary copy at upload time; no lasting local copy, and the key never
  leaves this machine.

## Evidence of Completion

- [x] `./init.sh`: lock check, ruff clean, 237 tests pass; experiment tests 203 pass, 2 skipped locally (the T3
      end-to-end tests, which pass on the DGX: 50 T3 tests in all).
- [x] DGX: torch 2.14.1+cu130 with CUDA, jax 0.11.2, moscot 0.5.2, scvi-tools 1.5.1 import (runs/logs/setup2.log).

## Notes for Next Session

- Slack sweeps: leave the group DM C0CA7QTEMDF out of date searches that reach back to 2026-10-18 (they surface our
  leak report). Six workspace accounts use eigenlabs.org addresses (none has posted); record their posts, or posts in a
  new Yukon channel, by place and date only.
- The Data Protocol (label table, units, sectioning) was promised before the test phase opens (an organiser,
  2026-10-16): look for it on the Data page this week (never the board file).

- `neighborhood_mmd` builds 15-NN neighbourhoods over all submitted cells, so submitted density moves T2 scores.
- Ctnnb1 is not on the 500-gene MERFISH panel: gene_ko gives nothing for the beta-catenin test board.
- T1 E12.5 is a heart-chamber-only dissection; earlier stages include surrounding tissue.
- DGX research scratch removed on 2026-10-19 (user's OK); `~/v/vec-gapconc` keeps only `ext/` and `out/`, the
  legal inputs the committed filters in `experiments/external/` read (the JAX files are hard-linked into data/external).
- `experiments/harmonise_review/review.py` has 5 pre-existing ruff findings (outside the checked `src tests`).
- Release day (2026-10-27): re-check #virtual-embryo-challenge and the Rules/Data pages (never the board file) for
  the release time, submission rules, scorer version and any Yukon details the organisers post; an organiser said on
  2026-09-01 that T2 E8.5 is released when the final stage starts. Before the first upload, confirm that our account
  sign-in works (a participant saw "Sign-in is not available yet" on the account page on 2026-10-18).
