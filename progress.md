# Session Progress Log

## Current State

**Last Updated:** 2026-10-21
**Active Feature:** feat-004 - Final-phase intake on 2026-10-27, postponed from 2026-10-20 by the organisers
(`docs/P3_RUNBOOK.md`); feat-003 completes with the release files. Before then: the release-day rules sign-off
(`docs/RELEASE_RULES_SUMMARY.md`) and the v1 pipelines on the shared core (`vec.core`, feat-020 done)

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
      (two questions) waits in the user's drafts. Not checked by agents: the challenge account messages page (user
      only) and Gmail Spam.
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
- [ ] The user's call, not covered by approved answer 7: a winning count level below a test board's lower limit
      (truth_cells or min_cells; only 0.25 x max_cells can fall there). The build stops and asks. Recommended: the
      next passing level the board allows, else max_cells.
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
