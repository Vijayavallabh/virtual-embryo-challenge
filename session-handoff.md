# Session Handoff

## Current Objective

- **State at 2026-11-26, 01:25 UTC (newest; it supersedes everything below where they differ).** Decision log through
  2026-11-26 (d); next (e). Ruling rows through R76; next R77 (heart), R78 (embryo), R79 (T3) for round 4. Tries: T3
  89. The last upload is 2026-11-26 11:59 UTC.
- **STANDING 204.865** (unchanged): T1 74.75 `T1_s2rcMSBm10.h5ad`; T3 63.94 `T3_c2_s2.h5ad`; heart 62.45
  `T2heart_organ_ts125_tk_half.h5ad`; embryo 69.9 `T2embryo_s2x_ov.h5ad`.
- **DUE NOW (23), all handed over** (`docs/UPLOAD_SHEET.md`): T1 eleven (R73's four, R67's MSBm075, six B2/B25
  search files); T3 five (R76: s085, s05, g15, g05, g2_s067); embryo five (R75: ovt712, ovt56, ovt512, s2xs_ovt23,
  ovd); heart OTK10 and OTK15 (R74). WAITING: OTK075 only if OTK10's headline is strictly above OTK125's, OTK175 only
  if OTK15's is, each read once at the first headline and at the latest 06:00 UTC; A20 (last on T1).
- **HEADLINES:** record with `lift8/record_headlines.py <repo> <label> "<YYYY-MM-DD HH:MM>" <spec.json>` (headline
  only; the rules' readings go in the spec in the caller's words; conditional files released through `due_add` and
  `sheet_edits`; dry-run on a throwaway copy passed). Fix by hand any sheet line it prints for a new best.
- **IN FLIGHT, round 4** (`lift8/round_1125_gap/TASK_ref4.md`): three blind designers write HEART/EMBRYO/T3_REF4_design.md
  (due about 02:20 UTC): next steps named before round 3's headlines, each sent at once or on a headline condition read
  once by 10:30 UTC. Then cold challengers, one short question per board (no recommendation where recused; T3 depth and
  the heart's scope and size need the captain's named words per file), ruling rows R77 to R79 placed before code,
  engineers in worktrees, an independent review with the blocker standard stated, merge with the guard tests, DGX gate
  dry-run, builds twice with one sha256 (target about 07:30 UTC), record, hand-over.
- **DGX state:** `~/v/scratch/t1line-wt` at ead4f34 (T1); `~/v/vec-campaign-c` at 726aaac (T3, R76);
  `~/v/vec-campaign` at ee22b75 (embryo, R75); `~/v/scratch/h29b-wt` at 4b3e7a9 (heart, R74; pyarrow installed).
  Moved aside, never deleted: `~/v/scratch/bf_moved/`, `~/v/scratch/oic_moved/`, `~/v/scratch/t3_moved/k31`. GPUs 3 to
  6 run another job's vLLM server; leave it. Local disk: / 27 GB free, /mnt/md0 175 GB free.
- **State at 2026-11-25, 21:55 UTC (newest; it supersedes everything below where they differ).** Decision log through
  2026-11-25 (hhhh); next (iiii), or (a) on 2026-11-26. Ruling rows through R73; next R74. Tries: T3 89 (N65 the
  round-8 designer's breach, 85; K31, 86 to 89). The last upload is 2026-11-26 11:59 UTC.
- **STANDING 204.865:** T1 74.75 `T1_s2rcMSBm10.h5ad`; T3 63.94 `T3_c2_s2.h5ad`; heart 62.45
  `T2heart_organ_ts125_tk_half.h5ad`; embryo 69.9 `T2embryo_s2x_ov.h5ad`. `docs/UPLOAD_SHEET.md`: 11 T1 files due (row
  R73's four, MSBm075 of R67, and six B2/B25 search files); nothing due on T3, heart or embryo.
- **CLOSED TODAY, nothing built:** heart WBH (H37, STOP at G0), embryo BF (E59, STOP), heart OIC125 (H38, STOP at G0:
  2 of 72 heart regions confirmed, both sources dropped), the T3 sign cube (K31, STOP; its first run INCOMPLETE on a
  start-file race, rerun after a code fix per K18's rule); round 8's MOON3 and MOON4 struck at challenge. Built and
  handed over: R67's three T1 files, R68's four hedges (all four scored below their boards' bests), R73's four T1
  files.
- **DGX state:** `~/v/scratch/t1line-wt` at ead4f34 (T1 builds); `~/v/vec-campaign-c` at b2e0c6d (T3; untracked K31
  reading files there equal main's); `~/v/vec-campaign` at b282384 (BF's reading files there equal main's);
  `~/v/scratch/h29b-wt` at 9864866 with pyarrow installed (its G0 records equal main's). Files moved aside, never
  deleted: `~/v/scratch/bf_moved/`, `~/v/scratch/oic_moved/`. The OIC125 fetch (about 57 GB) sits in
  `~/v/vec-ext/heart_oic` with its MANIFEST (committed as `experiments/heart_ext/results/fetch_MANIFEST.json`). GPUs 3
  to 6 run another job's vLLM server; leave it.
- **State at 2026-11-25, 18:00 UTC (newest; it supersedes everything below where they differ).** Decision log through
  2026-11-25 (qqq); next (rrr). Ruling rows through R71; next R72. Main at c15e09d plus this bullet. The last upload
  is 2026-11-26 11:59 UTC.
- **STANDING 204.865:** T1 74.75 `T1_s2rcMSBm10.h5ad`; T3 63.94 `T3_c2_s2.h5ad`; heart 62.45
  `T2heart_organ_ts125_tk_half.h5ad`; embryo 69.9 `T2embryo_s2x_ov.h5ad`. `docs/UPLOAD_SHEET.md` holds every due file.
- **IN FLIGHT (blind agents, each in its own worktree under `.cache/vec-tmp/`, branched from c15e09d):** R67, T1's
  damping refinement (`campaign/t1-r67`: MSBm125, MSBm15, MSBm075); R68, the four hedge files (`campaign/hedge-fixed`,
  `experiments/builds/hedge_fixed.py`); R69, the heart's world-blur file WBH (`campaign/heart-wbh`; its gate row H37
  is drafted to `lift8/round_1125_gap/H37_row.md` and goes into TRIES.md before any read); R71, the embryo's consensus
  form BF on O and SO only (`campaign/embryo-bf`; gate row drafted to `EBF_row.md`, placed before any read); R70, the
  heart's OIC125: two independent readers pin each outside sample's stage and licence (`OIC_readerA.md`,
  `OIC_readerB.md`), then register rows before any fetch; a cold challenge of its design (`HEART_MOON1_challenge.md`)
  runs beside them. Each branch then gets an independent review, a merge with the guard tests, a gate dry-run on the
  real records, DGX builds twice with one sha256, a record and a hand-over.
- **OPEN WITH THE CAPTAIN:** "depth line too" would add BF on O_hi and W (R71's words do not cover the depth line).
- **COMPUTE:** DGX GPUs 3 to 6 run another job's vLLM server (not ours; leave it alone); ours are 0 to 2 and 7. The
  DGX root disk is full, so TMPDIR stays under `~/v/scratch`. Local: / 28 GB free, /mnt/md0 176 GB free.
- **State at 2026-11-25, 13:55 UTC (newest; it supersedes everything below where they differ).** Decision log through 2026-11-25 (rr); next (ss). Ruling rows through R66; next R67. Gate rows done today: N64, T72, T73, H35, H36, E57, E58, K30. Tries: T1 480, heart 127, embryo 95, T3 84. The last upload is 2026-11-26 11:59 UTC.
- **STANDING 204.465:** T1 74.35 `T1_s2rcMS20m10.h5ad`; T3 63.94 `T3_c2_s2.h5ad`; heart 62.45 `T2heart_organ_ts125_tk_half.h5ad`; embryo 69.9 `T2embryo_s2x_ov.h5ad`. Uploads unlimited (R42); headlines only (R22); `docs/UPLOAD_SHEET.md` holds every due file (T1 only now: ten files, then A15 last; heart, T3 and embryo fully uploaded).
- **NEXT, T1 round 3 (row R66):** K1 `T1_s2rcMS20H10m10` and C3 `T1_s2rcMS20H20` hold. Code on branch `campaign/t1-r66` (1c30ad5), under independent review (result to `lift8/round_1125_gap/r66_review.json`); then merge with the guard tests (t2_build literals, n54 listing, s2_rules), deliver by bundle, build each with `experiments/builds/t1_r66_line.sh <COMMIT> <RUN> <FILE>` (gate dry-run passed at 1c30ad5 on `~/v/scratch/t1line-wt`; the second build on `~/v/vec-campaign-d`), record with `lift8/record_t1r66.py`, hand over.
- **NEXT, round 4 (big swings; `lift8/round_1125_gap/TASK_gap4.md`):** T3 none inside the rules; T3 under the captain's waiver of the team's T3 bar (log (pp)): `T3_BIO_design.md` (pathway databases pick genes, C2's tables the sign; 0.30 / best +4), challenger running; embryo `EMBRYO_SWING_design.md` (E8.0 twin; closed lever, pre-set rule declines without an exception), challenger running; heart `HEART_SWING_design.md` (MRM; about 3%; closed lever), challenger running. Each goes to the captain as one short question; heart and embryo with no view from the AI assistant (recusals), T3 depth with the pre-set rule only, and no view on which genes a T3 file moves.
- **COMPUTE:** all eight DGX H100s free (the captain killed every GPU process, about 13:30 UTC); DGX root disk full, so TMPDIR stays under `~/v/scratch`; local data disk 178 GB free after removing 329 clean, pushed agent worktrees (88 kept: unpushed, unclean or locked). The captain wants the DGX used more and the local machine less.
- **KNOWN:** `tests/test_t1_cal_read.py` fails one test (stream 48 is also `COPY_PARITY_STREAM` in `src/vec/t2_build.py`); not fixed.
- **State at 2026-11-24, 09:05 UTC (newest; it supersedes the bullets below where they differ).** Decision log through
  2026-11-24 (mm); next label (nn). Ruling rows through R41; next R42. Register through X020. Next ledger rows: E56,
  H35 (H33 is never used), K30, N61, T72. Tries: T1 463, embryo 92, heart 105, T3 80. NOTHING RUNS: no agent, no DGX
  run of ours. Main is at the commit of this bullet, pushed.
  - **STANDING** 201.38 (T1 72.29 `T1_s2rcB2.h5ad`; T3 63.94 `T3_c2_s2.h5ad`; heart 60.4 D; embryo 69.9 O). No score
    reported since log 2026-11-23 (gg).
  - **UPLOADS (row R34):** 20 over 24 and 25 Nov, no cap per board, on the captain's word; none reported used yet.
    `docs/UPLOAD_SHEET.md` (state 08:50 UTC) holds every built file and its condition. DUE NOW (8): T1
    `T1_s2rcB25.h5ad`, `T1_s2rcCm20.h5ad`; T3 `T3_k23_s2.h5ad`, `T3_c1_s2_g2.h5ad`; embryo `T2embryo_s2xs_ov.h5ad`,
    `T2embryo_s2.h5ad`, `T2embryo_s2x_ovhi.h5ad`; heart `T2heart_organ_ts05_half.h5ad`. WAITING ON A SCORE: B30 (B25 >
    72.29), Cm25 (Cm20 > C's 71.8), B2d15 (Cm20 > 72.29; row R35), B2d20 (B2d15 > the higher of the 2.0 file and
    Cm20), CK and the E = 4 file (each K23 > the barred s2's 62.30), OTS10 (OTS05 > D's 60.4), and the new T1 files
    below. T3 `T3_c2_s2_g2.h5ad`: no condition, after the T3 files queued now. LAST ON T1: MSB, A15, then A20 (A15 >
    A's 70.8).
  - **BUILT SINCE 01:50 UTC, each twice with one sha256 at a clean commit after the identity gates, recorded before
    any score:** `runs/manual_upload/T3_c2_s2_g2.h5ad` sha256
    014d2bf9d8ba7fdc56dcfba2a3b8aa77f0ea20bb806ebc2c1193e483d62712ac (log (ll); record
    `experiments/builds/results/t3c2g_fc57ca9_g1.txt`); `T1_s2rcCm30.h5ad`
    807fe097c88ee779430896002a34635aeffbc99b433c68ebdb7bd41d8af2cbbe (line 7: only if Cm25 > Cm20);
    `T1_s2rcB25d15.h5ad` 5bc2b32e4688439f29c54c1b1203325595d996c37df30ba2b9ff26f27a60f4c6 (line 6: only if Cm25 > both
    the 2.5 and the 2.0 file); `T1_s2rcB25d20.h5ad` e6520dbee0fdc3f1336bc3f9de1e99e8a501de2d440b3f9085c88b5c7ecfb3f6
    (only if B25d15 > Cm25); `T1_s2rc_pooled.h5ad` 91341a5f8c71dcc103a47f9c19a5d7794b1358c1c642c5c9251624eb6feb8bdc
    (SKIPPED by the captain's word, never uploaded) (log (mm); record
    `experiments/builds/results/t1l67pb_ce1fb91_l1.txt`). Data used: T1 files as R13's (released data and GEO
    GSE193346, row X002); the T3 file released data only (the table gata4_pooled at weight 2 inside C2's size-held
    product).
  - **READS SINCE 01:50 UTC (the pre-set statements only):** K28 "Try K28: STOP at G0: the database bars are not met;
    no fold read, no try, no file."; H34 "Try H34: STOPPED (a in CLEAN): OTS15 is not built."; E55 "STOP: no pass, no
    file. The handed-over embryo files stand."; K29 "Try K29: STRUCK (S3, S4): no file; C2 stays." Lapses recorded as
    rows N59 (embryo try 92) and N60 (T3 try 75).
  - **THE CAPTAIN'S WORDS SINCE 01:50 UTC (rows R31 to R41, log (m) to (kk)):** run CHIPQ-C and lift K27's stop for it
    (R31, R33; it stopped at G0); 20 uploads (R34); "skip TS10, TS15 and TS20" (third heart line ended); a blind round
    on the closed levers (heart nothing; embryo struck; T1 pooled base yes then skipped; T3 cD7 yes, bound to its
    read, struck); "Yes: T1_s2rcCm20's headline may serve line 5; line 5 stands." (R35); full approval for anything
    that could improve the scores (log (v): effort, compute and agents; not the organisers' rules, closed levers or
    named words); two extra T1 steps (R39, R40); "yes, T3_c2_s2_g2" (R41).
  - **WHEN HEADLINES COME:** record the headline per board only in `results/manual_uploads.csv` (stop reading at the
    headline), one log entry, apply each condition above by its letter, update the sheet, name the next file per
    board, then an independent rule check through the Agent tool (score-blind agents never through Workflow).
  - **SESSION NOTE.** The session restarted at about 08:00 UTC: the main session's working scripts and round files
    moved from the retired scratchpad to `/mnt/md0/IITM/BackUp/Home/vijayavallabh/.cache/vec-work/lift8` (owner only);
    agents of the earlier session cannot be resumed. All new agents run on Opus (the earlier model is rate limited
    until after the deadline). DGX: `~/v/vec-campaign-c` is detached at fc57ca9 (K29's untracked results copy moved,
    not deleted, to `~/v/scratch/k29_results_from_campaign_c`); `~/v/scratch/t1line-wt` at ce1fb91. The DGX's root
    disk has about 1 GB free: tests and temp go under `~/v/scratch`.

- **State at 2026-11-24, 01:50 UTC (superseded by the bullet above where they differ).** Decision log through
  2026-11-24 (l); next label (m). Ruling rows through R30; next R31. Register through X020. Next embryo row E56. Tries
  unchanged: T1 463, embryo 88, heart 97, T3 74 (row E55's three embryo tries, 89 to 91, count only when its run step
  starts). RUNNING: one reviewer blind to test scores on the embryo candidate's code, and a target-free pre-check on
  the DGX (below). Main is at a5d1823, pushed.
  - **STANDING** 201.38 (T1 72.29 `T1_s2rcB2.h5ad`; T3 63.94 `T3_c2_s2.h5ad`; heart 60.4 D; embryo 69.9 O). No score
    reported since log 2026-11-23 (gg).
  - **ORDERS NOW** (`docs/UPLOAD_SHEET.md`, state 01:42 UTC: 19 built files in the orders, 3 skipped). T1:
    `T1_s2rcB25.h5ad` (due), `T1_s2rcCm20.h5ad` (due); B30 only if B25 beats the 2.0 file's 72.29; Cm25 only if Cm20
    beats C's 71.8; `T1_s2rcB2d15.h5ad` only if Cm20 beats 72.29; `T1_s2rcB2d20.h5ad` only if B2d15 beats the higher
    of the 2.0 file's and Cm20's headline; then, kept last, `T1_s2rcMSB.h5ad`, `T1_s2rcA15.h5ad` (A20 only if A15
    beats A's 70.8). T3: `T3_k23_s2.h5ad`; CK and the E = 4 file each only if K23 beats the barred s2's 62.30; then
    `T3_c1_s2_g2.h5ad`. Embryo: `T2embryo_s2xs_ov.h5ad`, `T2embryo_s2.h5ad`, `T2embryo_s2x_ovhi.h5ad`. Heart:
    `T2heart_core_ts10_half.h5ad` (due), `T2heart_organ_ts05_half.h5ad` (due), OTS10 only if OTS05 beats D's 60.4. The
    main session's one recommendation on order: T1 'B25 next, then Cm20' (from ordinal facts; no heart or embryo
    lean).
  - **THE CAPTAIN'S WORDS SINCE 00:45 UTC.** (1) "keep MSB and A15 if they can improve our scores" = the three T1
    files skipped an hour earlier are kept after all and stand LAST in T1's order (log (i); rules note re-pinned; rows
    R26 and R22); 'MSB first' restores their old places. (2) "check on the general-purpose subagent that's been
    waiting since saturday" = answered: that agent (the K24 review) completed on 21 Nov; nothing was pending.
  - **T1 LINE 5 FILES ARE BUILT AND HANDED OVER (log (k), commit 131cef5).** `runs/manual_upload/T1_s2rcB2d15.h5ad`,
    sha256 fe4e7403d6257cc265a02d0a13d40cdb60a97e9fdc0647d6b5c7e318aedaed6c; `runs/manual_upload/T1_s2rcB2d20.h5ad`,
    sha256 bd12c66312ba32ab938b4ee71bd5c79e58a154d6fc56434012e6dad907e668dd. Build `experiments/builds/t1_b2d_line.sh
    7e7d9a6... b1` on the DGX worktree `~/v/scratch/t1line-wt`: every exit record 0, the eleven identity commands
    rebuilt the recorded bytes, each new file built twice with one sha256 (record
    `experiments/builds/results/t1b2dline_7e7d9a6_b1.txt`). Recorded before any score of them. Data used: as R13's
    file.
  - **EMBRYO CANDIDATE OSR (rows and entries: log (j), (l); ledger row E55).** The target-free gate passed (run g1 at
    d145d78; record committed under `experiments/coords/results/osr_gate_g1/`, gate.json sha256 7c0fc6b6...). Row E55
    was written by a rule writer blind to test scores, checked cold (13 corrections), completed by two rounds of
    rulings of the specification's author, and PLACED at a5d1823 after the J7 rerun of the gate (run g2) reproduced
    every g exactly (largest change 0.0 over 100 pair-seeds). Nothing has been read against a practice target. The
    code (form `overlay_stretch`, arm `fgw_shape0x_ov_sr`, placebos `_anti` and `_axes`, builder option
    `--shape-interp overlay_stretch`, read script `experiments/coords/e55_overlay_stretch_read.py`, launcher
    `experiments/coords/run_e55.sh`, tests `tests/test_e55_overlay_stretch.py`) is on branch `campaign/emb-osr` at
    05d9a21 (base 1f0c9ad; it merges `campaign/emb-osr-gate`), NOT merged: an independent reviewer blind to test
    scores, with the blocker standard in its brief, was started at 01:36 UTC through the Agent tool and writes
    `lift8/round_1124_t2/e55review.json`.
  - **DGX STATE FOR E55.** Worktree `~/v/scratch/e55-wt` (of `~/v/vec-campaign-c`) is detached at 05d9a21, clean. The
    launcher's refusal paths were traced there under a throwaway home folder (exits 2, 90, 91, 93, 99, each with its
    record); its tests step, run by hand, passes in `~/v/vec-campaign-c/.venv` (297 passed). A target-free pre-check
    of step 0 runs since 01:43 UTC: the gate at the branch commit in that worktree (run g3, scratch
    `~/v/scratch/osr_gate/g3`), to be compared with g1 (the read's step 0 voids the row if any g moves by more than
    1e-9, so the cause of any difference must be found BEFORE the real launch).
  - **NEXT FOR E55.** Reviewer READY (after any fixes by the engineer and a re-check): merge `campaign/emb-osr` into
    main (no assistant trailer), push, deliver the merge commit to `~/v/vec-campaign-c` as a ref, `git checkout
    --detach` it in `~/v/scratch/e55-wt`, then launch detached, with PATH and the uv variables set:
    `E55_VENV=/home/prachh/v/vec-campaign-c/.venv E55_DATA=/home/prachh/v/vec-campaign-c/data/raw bash
    experiments/coords/run_e55.sh <full commit> r1` (scratch `~/v/scratch/e55/r1`, final record `done.exit`; a
    background waiter, never a foreground wait). The three tries count when the run step starts. The reading lands in
    `experiments/coords/results/e55_overlay_stretch/` in the worktree: copy it to main, fill row E55's last cell with
    the pre-set statement only. PASS = one priced naming question to the captain, labelled 'diagnostic, selected on
    the nine tasks'; the file is built (twice, one clean commit) only if the captain names it. STOP or VOID = no file,
    nothing further under this row. The challenger's odds, committed before data: about 0.02 overall.
  - **UNANSWERED BY THE CAPTAIN (defaults stand).** Heart: 'yes to TS15 and TS20' and 'yes to OTS15' (default no; each
    needs a practice read with a stop rule first); 'run CHIPQ-C' (default no; about two in a hundred); the T3 file
    `T3_c1_s2_g2.h5ad` keeps its place after K23, CK and the E = 4 file unless the captain places it.

- **State at 2026-11-24, 00:45 UTC (superseded by the bullet above where they differ).** Decision log through
  2026-11-24 (h); next label (i). Ruling rows through R30; next R31. Register through X020. Tries unchanged: T1 463,
  embryo 88, heart 97, T3 74. RUNNING: the T1 line-5 build on the DGX (below) and one engineer blind to test scores
  writing the embryo candidate's target-free gate (branch `campaign/emb-osr-gate`).
  - **STANDING** 201.38 (T1 72.29 `T1_s2rcB2.h5ad`; T3 63.94 `T3_c2_s2.h5ad`; heart 60.4 D; embryo 69.9 O). No score
    reported since log 2026-11-23 (gg).
  - **ORDERS NOW** (`docs/UPLOAD_SHEET.md`: 14 built files in the orders, 6 skipped). T1: `T1_s2rcB25.h5ad` (due),
    `T1_s2rcCm20.h5ad` (due); B30 only if B25 beats the 2.0 file's 72.29; Cm25 only if Cm20 beats C's 71.8; then line
    5's files once built (`T1_s2rcB2d15.h5ad` only if Cm20 beats 72.29; `T1_s2rcB2d20.h5ad` only if B2d15 beats that
    line's best). T3: `T3_k23_s2.h5ad`; CK and the E = 4 file each only if K23 beats the barred s2's 62.30; then
    `T3_c1_s2_g2.h5ad` (BUILT, sha256 e01e24facb68..., handed to the captain). Embryo: `T2embryo_s2xs_ov.h5ad`,
    `T2embryo_s2.h5ad`, `T2embryo_s2x_ovhi.h5ad`. Heart: `T2heart_core_ts10_half.h5ad` (due),
    `T2heart_organ_ts05_half.h5ad` (due), OTS10 only if OTS05 beats D's 60.4.
  - **THE CAPTAIN'S WORDS SINCE 00:10 UTC.** (1) "[loosen it all the way]" = ruling row R30: the T3 line's 5% guard
    gates nothing; the file was coded (review READY), built twice (a3c511c; gates rebuilt the recorded bytes) and
    recorded (log (d), (f); rules notes re-pinned). (2) "[I want it to boost our scores as drastic as possible ...]" =
    a design round for the heart and the embryo (log (e)); it is NOT a named exception for any file. (3) "skip them if
    they don't improve" = `T1_s2rcMSB.h5ad`, `T1_s2rcA15.h5ad` and with it `T1_s2rcA20.h5ad` are skipped (log (h);
    rules note re-pinned; rows R26 and R22); the captain reverses with 'keep MSB' or 'keep A15'.
  - **T1 LINE 5 BUILD.** Code reviewed READY and merged at 7e7d9a6; delivered to the DGX worktree
    `~/v/scratch/t1line-wt` (clean at that commit); refusal path traced (exit 96, record written); the real run
    `experiments/builds/t1_b2d_line.sh 7e7d9a64a87aaf7f766ac79f9907e230bf6dacf6 b1` was launched at 00:35 UTC (scratch
    `~/v/scratch/t1b2dline/b1`, output `runs/t1b2dline_7e7d9a6_b1`; final record `done.exit`; waiter
    `lift8/wait_dgx_run.sh v/scratch/t1b2dline/b1 done.exit`). WHEN DONE: copy both files to `runs/manual_upload/`,
    check their sha256 against both DGX builds, write the build record, the rules note with the sha256 (re-pin), the
    sheet rows and a log entry, all BEFORE any score of them; hand over with the T1 Data-used text.
  - **THE HEART AND EMBRYO ROUND (log (e), (g)).** Heart: two line texts beyond the lines' caps are on record and NOT
    in force (TS15 and TS20 on the full trim; OTS15 on the organ cut); each needs the captain's own named yes ('yes to
    TS15 and TS20', 'yes to OTS15'), default no, and a practice read with a stop rule before any build; the new heart
    candidate was struck. Embryo: no search line (designer and cold challenger agree); new candidates: one designer
    found none; a second found one, OSR (the overlay file with one stretch to the brackets' log-midpoint elongation),
    which its cold challenger put at about 0.02 and sent to its target-free gate only
    (`lift8/round_1124_t2/embryo_osr_attack.json`). The test pair's blocks are 1|2,3 by the committed E49 reading, so
    stop clause S3 does not fire. WHEN THE GATE SCRIPT COMES BACK: deliver the branch commit to a clean DGX checkout
    with the embryo data (`~/v/vec-campaign-c` has data and an environment), run `experiments/coords/osr_gate.sh
    <COMMIT> <RUN>` detached, read `gate.txt`: STOP ends the candidate; PASS means a ledger row (rule writer and
    challenger, blind), then code, review, the practice read with the challenger's pass rule, and only after a pass
    the captain's naming. Not yet in the decision log: the embryo part of the round (write it with the gate's
    outcome).
  - **BLIND AGENTS:** agent tool only, tripwire in the brief, ledger only through `lift8/ledger_blind.py` (it now also
    withholds cells that speak of a test upload's per-metric result). The ledger's rule J8 and rows N53, E47a and E49a
    say in words that the first embryo upload met a shape loss; a designer told it was fully blind recused itself from
    the shape group on meeting them, so embryo designers are now told what J8 allows.
- **State at 2026-11-24, 00:16 UTC (superseded by the bullet above where they differ; everything not named
  here is as the bullet of 00:10 UTC).** Decision log through 2026-11-24 (d); next label (e). Ruling rows through R30;
  next R31. RUNNING: two engineers blind to test scores, each in its own worktree (T1 line 5 on `campaign/t1-line5`;
  T3's second line on `campaign/t3-gata4-weight-2`), and one rehearsal build on the DGX.
  - **THE CAPTAIN LIFTED THE 5% GUARD OF T3's SECOND LINE** (own words "[loosen it all the way]" under a side note;
    ruling row R30; log 2026-11-24 (d); rules amendment re-pinned). The line is alive again: its one file
    `T3_c1_s2_g2.h5ad` (C1's command with `--gata4-exponent 2.0 --captain-ruling R22`, phi_M * phi_G ** 2, nothing
    rescaled) is to be built; the guard's figures are recorded and stop nothing; the check that the factor applied is
    the factor computed stays. Its place: after row K23's file, CK and the E = 4 file on their conditions, unless the
    captain says otherwise.
  - **THE CHEAP CHECK FIRST, THIS TIME:** a rehearsal of the existing launcher runs on the DGX in `~/v/vec-campaign-c`
    at main 015cd3e (`build_two_knockout.sh <commit> s2`; scratch `~/v/scratch/t3two_015cd3e_s2`, output
    `runs/t3_two_knockout_015cd3e_s2`), launched 00:14 UTC, with a local background waiter. It must rebuild:
    identity-s1 3712c3f9988d, identity-s2 1b69521e3c28, c1-s2 8718330866429 (C1), c2-s2 47111ac1d07e (C2). If a hash
    differs, the checkout's environment has drifted and must be put right before the real run. That checkout's
    untracked K27 output was moved, not deleted, to `~/v/scratch/k27_results_from_campaign_c`.
  - **T3 LINE 2, WHAT REMAINS:** the engineer's hand-back; a reviewer blind to test scores (agent tool; blocker
    standard up front); merge; deliver the commit to `~/v/vec-campaign-c`, check it out detached and clean; trace a
    refusal path of the new launcher argument; `build_two_knockout.sh <commit> s2 g2` detached with a background
    waiter; copy the file to `runs/manual_upload/`; its sha256 into log entry 2026-11-23 (ff)'s line and the sheet
    BEFORE any score; hand over with the Data-used text of the line (released challenge data only: wild-type E8.75
    cells, the released Mab21l2 knockout and the released Gata4 validation knockout files; no outside data).
- **State at 2026-11-24, 00:10 UTC (superseded by the bullet above where they differ).** Decision log through
  2026-11-24 (c); next label (d). Ruling rows through R29; next R30 (row R22's first cell was amended twice,
  add-only). Register through X020; next X021. T1's next row is T73. Tries unchanged: T1 463, embryo 88, heart 97, T3
  74. RUNNING: one engineer, blind to test scores, on T1 line 5's builder pair and launcher (its own worktree, branch
  `campaign/t1-line5`). Nothing runs on the DGX. Main at the commit of this note.
  - **EIGHT HEADLINES WERE REPORTED ON 2026-11-23** (log (cc) and (gg); headline only). Best per board: T1
    `T1_s2rcB2.h5ad` 72.29; T3 `T3_c2_s2.h5ad` 63.94; heart `T2heart_organ_half.h5ad` 60.4; embryo
    `T2embryo_s2x_ov.h5ad` 69.9. Standing 201.38. An independent reader confirmed the rule application and every
    sha256 of the sheet after each set.
  - **ORDERS NOW** (`docs/UPLOAD_SHEET.md`: 16 built files in the orders, 3 skipped). T1: `T1_s2rcMSB.h5ad`, then at
    each slot line 1's due file (`T1_s2rcB25.h5ad`, due; B30 only if B25 is strictly above the 2.0 file's 72.29), line
    2's (`T1_s2rcA15.h5ad`, due; A20 only if A15 above A's 70.8), line 3's (`T1_s2rcCm20.h5ad`, due; Cm25 only if Cm20
    above C's 71.8), then line 5's (below; not built yet). T3: `T3_k23_s2.h5ad`; CK and the E = 4 file each only if
    K23 is strictly above the barred s2's 62.30. Embryo, no conditions: `T2embryo_s2xs_ov.h5ad`, `T2embryo_s2.h5ad`,
    `T2embryo_s2x_ovhi.h5ad`. Heart: `T2heart_core_ts10_half.h5ad` (due; the first line ends at it), then
    `T2heart_organ_ts05_half.h5ad` (due); OTS10 only if OTS05 above D's 60.4.
  - **TOLD TO THE CAPTAIN, NOT ANSWERED.** T1 order: the AI assistant's recommendation, from ordinal facts only, is
    `T1_s2rcB25.h5ad` next and then `T1_s2rcCm20.h5ad` (each one step from the file holding the best), ahead of MSB
    and A15. Heart order: no recommendation; the captain may put OTS05 before TS10. Without the captain's words the
    order on record stands.
  - **THE CAPTAIN'S WORDS OF LOG 2026-11-23 (hh)** ('yes for both if it would improve our scores'; '500 external
    datasets ... 1 TB ... at all drastic costs') AND WHAT CAME OF THEM. (1) T3's second line (the Gata4 table's weight
    at 2 on C1) was put in force and ENDED AT ITS OWN GUARD in an early trial on released data (38.2% against the 5%
    limit; log 2026-11-24 (a)); no file; the option stays named in R22's first cell and no builder takes it. (2) 'T1
    line 4' was not built and lapses. (3) T1 LINE 5 IS IN FORCE (log 2026-11-24 (c); rules amendment re-pinned; R22's
    first cell names `--g-mom-damping 1.5 --feng-fraction-mult 2.0` and `--g-mom-damping 2.0 --feng-fraction-mult
    2.0`): files `T1_s2rcB2d15.h5ad` then `T1_s2rcB2d20.h5ad`; its member is line 3's `T1_s2rcCm20.h5ad`; B2d15 goes
    only if Cm20's headline is strictly above the 2.0 file's 72.29, B2d20 only if B2d15 is above the line's best. The
    amendment rests on the AI assistant's reading of the captain's yes; the challenger's question was put to the
    captain unchanged (default: amend, line 5 last). (4) Outside data: a cold completeness check (log 2026-11-24 (b))
    agrees with the round of 23 Nov: nothing inside the rules that volume helps; one long shot, CHIPQ-C on T3, about
    two in a hundred, NOT STARTED; the captain orders it with 'run CHIPQ-C'. The T3 pathway tier stays out. (5) T3's
    design round of today wrote no line (a safety classifier stopped the designer's text; not pressed).
  - **T1 LINE 5, WHAT REMAINS** (all by 2026-11-25 12:00 UTC): the engineer's hand-back; a reviewer blind to test
    scores, started with the agent tool, blocker standard stated up front; merge; deliver the commit to the DGX
    (`~/v/vec-campaign-d`, detached worktree `~/v/scratch/t1line-wt`); trace the launcher's refusal path under a
    throwaway home; `experiments/builds/t1_b2d_line.sh <COMMIT> <RUN>` detached with a background waiter (eleven
    identity builds, then each new file twice); copy the two files to `runs/manual_upload/`; their sha256 into log
    entry (c) BEFORE any score; the upload sheet.
  - **STANDING RULES RESTATED.** Analysts, engineers and reviewers who must be blind to test scores are started ONLY
    with the agent tool (a workflow relays the captain's message, and the captain's latest messages hold scores), with
    the EXPOSED tripwire, the ledger only through `lift8/ledger_blind.py`, no git log, show or diff. A line's own
    guard or gate is run BEFORE its engineer (the T3 line's guard failed after its engineer had started). A text a
    safety classifier stops is not asked for again.
  - **WHEN THE CAPTAIN PASTES SCORES:** the headline per board only into `results/manual_uploads.csv`; a log entry
    (next label 2026-11-24 (d)); each rule by its letter; line 5's gate (Cm20 against the 2.0 file's headline); the
    sheet; the reply with the next file per board; then an independent rule check with the agent tool.
- **State at 2026-11-23, 23:30 UTC (superseded by the bullet above where they differ).** Decision log through
  2026-11-23 (ff); next label (gg). Ruling rows through R29; next R30. Register through X020; next X021. T1's next row
  is T73. Tries unchanged: T1 463, embryo 88, heart 97, T3 74. NOTHING RUNS, here or on the DGX; no workflow and no
  agent is open. Main at the commit of this note.
  - **HEADLINES AND ORDERS:** as the bullet of 23:10 UTC below, unchanged (standing 200.75; `docs/UPLOAD_SHEET.md`
    holds 20 files in the orders and 3 skipped; an independent reader confirmed the rule application and every
    sha256).
  - **RE-READING AGAINST THE NEW BESTS** (log (dd); scratch `lift8/record_expects_below_new_bests.json`): one file is
    now below its board's best by the record, the heart's `T2heart_core_ts05_half.h5ad` (its base returned 56.4 and
    the widest practice figure on record for the step is +1.38, under D's 60.4); not below: `T2embryo_s2x_ovhi.h5ad`,
    `T2embryo_s2.h5ad`; nine files without an expectation; the eight conditional files were not judged.
  - **FOUR ITEMS WAIT ON THE CAPTAIN'S OWN WORDS** (each was put in the reply; a side note forwarded without words of
    the captain's changed nothing). (1) 'Cm20 next': T1's third search first (`T1_s2rcCm20.h5ad`, then
    `T1_s2rcCm25.h5ad` if due); the main session recommended it, from ordinal facts only. (2) Heart: 'OTS05 next' or
    'skip TS05'; no recommendation (recused). (3) T1 line 4 (text in log (ff); scratch
    `lift8/round_1123b/T1_line_4_final.txt`): files `T1_s2rcCd15.h5ad` and `T1_s2rcCd20.h5ad`, g_mom's damping 1.5
    then 2.0 on C's command; it needs the captain's add-only amendment of row R22's first cell (`--g-mom-damping 1.5
    --feng-fraction-mult 1.5`, `--g-mom-damping 2.0 --feng-fraction-mult 1.5`) and, to get a slot, the captain's words
    naming the file Cd15 goes ahead of. On a yes: amend R22's first cell by record script, add the text to
    `docs/S2_RULES.md` and re-pin, then a blind engineer and reviewer, the nine identity rebuilds on
    `~/v/vec-campaign-d` (worktree `~/v/scratch/t1line-wt`), two builds with one sha256, and the sha256 into the log
    before any score; all by 2026-11-25 12:00 UTC. (4) T3 line 2 (text in log (ff); scratch
    `lift8/round_1123b/T3_line_1_final.txt`): one file `T3_c1_s2_g2.h5ad`, C1's command with a new `--gata4-exponent
    2.0 --captain-ruling R22`; VOID without the captain's yes, default no, and nothing of it is coded before the yes.
    On a yes: R22's first cell (add-only), the rules amendment, a blind engineer and reviewer, the 5% guard before any
    build, the identity gates (identity-s1, identity-s2, c1-s2, c2-s2), two builds with one sha256.
  - **PROCESS FAULT FOUND AND CORRECTED** (log (ee); `docs/DATA_POLICY.md` 2026-11-23 (e); rules note re-pinned).
    Every agent of a workflow receives the captain's triggering message verbatim ahead of its brief; an agent started
    with the plain agent tool does not. Twelve earlier workflows of analysts told they were blind had figures in that
    message, and today's first design round was voided. RULE FROM NOW: score-blind analysts are started ONLY with the
    agent tool, never as a workflow; their brief carries the EXPOSED tripwire; they read the ledger only through
    `lift8/ledger_blind.py` and run no git log, show or diff. The captain's latest message (a forwarded side note)
    holds headline values, so any workflow launched now relays them.
  - **WHEN THE CAPTAIN PASTES SCORES:** the headline per board only into `results/manual_uploads.csv`; a log entry
    (next label (gg)); each rule by its letter (T3's two conditions against the barred s2's headline, not the board's
    best); the sheet; the reply with the next file per board.
- **State at 2026-11-23, 23:10 UTC (superseded by the bullet above where they differ).** Decision log through
  2026-11-23 (cc); next label (dd). Ruling rows through R29; next R30. Register through X020; next X021. T1's next row
  is T73. Tries unchanged: T1 463, embryo 88, heart 97, T3 74. TWO WORKFLOWS ARE OPEN (below); nothing runs on the
  DGX. Main at the commit of this note.
  - **FOUR HEADLINES REPORTED** (log (cc); `results/manual_uploads.csv`, headline only; the paste's per-metric lines
    and ranks are not recorded and not used): T1 `T1_s2rcC.h5ad` 71.8 (new best), heart `T2heart_organ_half.h5ad` (D)
    60.4 (new best), T3 `T3_c1_s2.h5ad` 63.8 (new best), embryo `T2embryo_s2x_cr.h5ad` 63.99 (the board keeps 69.9).
    Standing 71.8 + 63.8 + (60.4 + 69.9) / 2 = 200.75 (197.04 before).
  - **ORDERS NOW** (`docs/UPLOAD_SHEET.md`: 20 files in the orders, 3 skipped). T1: `T1_s2rcB2.h5ad`,
    `T1_s2rcMSB.h5ad`, then at each slot line 1's due file (B25 only if the 2.0 file is strictly above B's 71.59, C's
    headline does not enter; B30 only if B25 above the 2.0 file), then line 2's (`T1_s2rcA15.h5ad`, due; A20 only if
    A15 above A's 70.8), then line 3's (`T1_s2rcCm20.h5ad`, DUE NOW since C is strictly above A; Cm25 only if Cm20
    above C's 71.8). T3: `T3_c2_s2.h5ad`, `T3_k23_s2.h5ad`; CK only after both are uploaded and each is strictly above
    the barred s2's 62.30 (the rule names base B's headline, NOT the board's best: the sheet's and the reply's earlier
    wording is corrected, log (cc)); the e4 file only if K23 is above 62.30. Embryo, no conditions:
    `T2embryo_s2x_crs.h5ad`, `T2embryo_s2xs_ov.h5ad`, `T2embryo_s2.h5ad`, `T2embryo_s2x_ovhi.h5ad`. Heart:
    `T2heart_core_ts05_half.h5ad` (TS05), then `T2heart_organ_ts05_half.h5ad` (OTS05, due since D's headline is
    reported); TS10 only if TS05 above the full trim's 56.4; OTS10 only if OTS05 above D's 60.4.
  - **PUT TO THE CAPTAIN, NOT ANSWERED:** Cm20 is the next step on the file that holds T1's best and stands last among
    T1's ready files by line 3's own order rule; the order is the captain's to change; the main session gave no lean.
  - **OPEN WORKFLOW 1, wf_61761ec3-e4f:** a re-reading of which waiting files the record expects below the NEW bests
    (four readers, four cold verifiers; the definitions of log (aa), with the rule that a pattern across test
    headlines is not evidence) plus one independent check of the rule application in the sheet and in entry (cc). When
    it returns: a log entry, a corrected sheet if the checker found a difference, and the report to the captain; no
    queue changes without the captain's word.
  - **OPEN WORKFLOW 2, wf_5a8901fb-2a7:** a score-blind design round for T1 and T3 under R22's 'Other lines' clause
    (designer, then cold challenger, per board; inputs `lift8/round_1123/inputs_T1.txt` and `inputs_T3.txt`, which
    hold ordinal facts only). If a line comes back 'place_as_final': place it in the log and the rules note before any
    build (row R22's first cell must name the option strings), then a blind engineer and reviewer, the identity builds
    and two builds with one sha256, all by 2026-11-25 12:00 UTC. The heart and the embryo got no round: a new file
    there needs the captain's named order.
  - **Blind analysts, widened today:** the ledger cells with first-day headlines are in rows R8, R9, R10 to R15, R17,
    K18 and T61 to T63; blind analysts read the ledger only through `lift8/ledger_blind.py` (ruling rows: first cell
    only), never open `progress.md`, `session-handoff.md`, `METHOD_SUMMARY.md`, `feature_list.json`, the top-level
    `results/`, `experiments/builds/results/t3two_fb01487.txt` or `tests/test_h28_heart_trend_strength.py`, and never
    run git log, show or diff.
- **State at 2026-11-23, 22:33 UTC (superseded by the bullet above where they differ).** Decision log through
  2026-11-23 (bb); next label (cc). Ruling rows through R29; next R30. Register through X020; next X021. T1's next row
  is T73. Tries unchanged: T1 463, embryo 88, heart 97, T3 74. NOTHING RUNS, here or on the DGX; no workflow is open.
  Main at the commit of this note.
  - **THE CAPTAIN SKIPPED THREE FILES** (own words "skip those three and update the upload sheet", log (bb); rules
    note re-pinned; rows R23 and R8): `T1_s2rc_MS.h5ad`, `T1_s2none.h5ad` (F), `T3_refused_s1.h5ad` are NOT uploaded
    unless the captain says otherwise; they stay on disk and stand in their own section of `docs/UPLOAD_SHEET.md`. The
    sheet now holds 24 files in the orders and 3 skipped (all checked against disk).
  - **ORDERS NOW.** T1: `T1_s2rcC.h5ad`, `T1_s2rcB2.h5ad` (the 2.0 file), `T1_s2rcMSB.h5ad`, `T1_s2rcA15.h5ad`; then
    on conditions B25 (only if the 2.0 file strictly beats B's headline), B30 (only if 2.5 beats 2.0), A20 (only if
    A15 strictly beats A's), Cm20 (only if C strictly beats A's), Cm25 (only if Cm20 beats C's). T3: `T3_c1_s2.h5ad`,
    `T3_c2_s2.h5ad`, `T3_k23_s2.h5ad`; then CK (only if C2 and K23 each strictly beat the board's best at that time),
    e4 (only if K23 does). Embryo: `T2embryo_s2x_cr.h5ad`, `T2embryo_s2x_crs.h5ad`, `T2embryo_s2xs_ov.h5ad`,
    `T2embryo_s2.h5ad`, `T2embryo_s2x_ovhi.h5ad`. Heart: D, TS05, TS10 on its condition, OTS05 after D's score, OTS10
    on its condition (rules in the bullet of 20:42 UTC below).
  - **THE CAPTAIN'S UPLOAD CONDITION OF LOG (z) STANDS AS STATED** (upload only where the expected score is closer to
    four per-board targets, which are NOT in the repository). The checked reading of log (aa)
    (`lift8/record_expects_below.json`): the record expects the three skipped files below their board's best; not
    below: `T2embryo_s2x_ovhi.h5ad`, `T2embryo_s2.h5ad`; no expectation for the other 13 that need no earlier score.
    When the captain pastes scores: headline per board only into `results/manual_uploads.csv`, one log entry, each
    rule applied mechanically, the sheet updated.
- **State at 2026-11-23, 20:42 UTC (superseded by the bullet above where they differ).** Decision log through
  2026-11-23 (y); next label (z), then (aa). Ruling rows through R29; next R30. Register through X020; next X021. T1's
  next row is T73. Tries unchanged: T1 463, embryo 88, heart 97, T3 74. NOTHING RUNS, here or on the DGX; no workflow
  is open. Main at the commit of this note.
  - **ADDED 22:23 UTC: A CHECKED LIST OF WHAT THE RECORD EXPECTS BELOW EACH BOARD'S BEST** (log (aa); next label (bb);
    scratch `lift8/record_expects_below.json`). A forwarded side note (no words of the captain's) suggested skipping
    only files the record expects below today's best; no rule or queue was changed. Four readers and four cold
    verifiers: BELOW: `T1_s2rc_MS.h5ad`, `T1_s2none.h5ad`, `T3_refused_s1.h5ad`. NOT BELOW: `T2embryo_s2x_ovhi.h5ad`,
    `T2embryo_s2.h5ad` (an old range, 60 to 74). NO EXPECTATION: the other 13 that need no earlier score. My reply of
    entry (z) was wrong on three points (the 2026-11-18 (k) plan table does hold expected board scores for two older
    files; the copy-rows embryo files are 'no expectation', not 'below'; F and the refused s1 are 'below'): corrected
    to the captain. LESSON: before saying 'the record holds none', search PLAN.md's plan tables as well as the ledger.
  - **ADDED 22:02 UTC: THE CAPTAIN'S UPLOAD CONDITION** (log (z); next label is now (aa)). The captain uploads alone
    and "only if the expected leaderboard scores are closer to these targets than our current scores", with four
    per-board target figures that are NOT in the repository (possible per-board leaderboard figures: never record
    them, use them for nothing but the subtraction asked). The comparison given: the record's central figure for every
    waiting file is no change; practice-sized upper ends are about a point or two (T1), under a point (T3 K23, heart
    trend), not forecastable (embryo, by the analysts' text), no estimate (D, the trend on D, most T3 files); nothing
    on the shelf comes near the targets. No estimate was raised; do not raise one if asked again (lead with the
    record). Uploads are the captain's acts.
  - **THE HEART'S SECOND SEARCH LINE IS BUILT AND HANDED OVER** (the captain's words "i want it too"; log (w) the
    order and my reading, (x) ruling row R29 and the challenger's final rule text in the rules note, (y) the builds;
    record `experiments/builds/results/heartlined_0d47a67_od1.txt`). `runs/manual_upload/T2heart_organ_ts05_half.h5ad`
    sha256 f918e76731aac18eda3b46c8a35c18f26ae00a1a6a1db36d33e7515704cc43ea (OTS05); `T2heart_organ_ts10_half.h5ad`
    sha256 ee9e0e68d62e8d1448f1dbbbce1146dc5fb4877e94d77fc289fa82744acb06a9 (OTS10). Built before any headline of D,
    TS05, TS10, OTS05 or OTS10 was reported (the rule's build clause). The first line's files were re-checked at the
    same commit (run ts2: same sha256 as run ts1).
  - **HEART RULES WHEN A SCORE COMES BACK (headline only; two independent lines; no headline of one line enters the
    other's comparisons).** First line (R28): TS05 goes after D; TS10 only if TS05's headline is strictly above the
    full trim C's (on record), whatever D scores; ends at TS10. Second line (R29): it waits until D is uploaded and
    D's headline reported; the line's best is D's FIRST headline; OTS05 is then its next file; OTS10 only if OTS05's
    headline is strictly above D's (a tie, a lower headline or no reported headline ends the line); if D is never
    uploaded, neither OTS file is uploaded. Slot order at each heart slot: D, then the first line's due file, then the
    second line's due file; at most two uploads a day until 2026-11-26 11:59 UTC; a file that finds no slot is not
    uploaded. Nothing is rebuilt, re-seeded, resized or changed after a score.
  - **UPLOAD SHEET:** `docs/UPLOAD_SHEET.md` lists 27 built files (checked against disk at 20:41 UTC). No score has
    been reported since 21 Nov. The other boards' rules are in the bullet of 20:10 UTC below.
  - **BRIEFING NOTE FOR SCORE-BLIND AGENTS (new):** ledger row H28's FIFTH cell describes in words (no values) on
    which metrics the first heart upload stood above or below the floor: give blind agents H28's cells 1 to 4 and 6 to
    9 only. The second line's challenger read it under my brief and reported it; it says nothing of it entered the
    rule text (log (x)).
  - **DGX:** `~/v/scratch/h29b-wt` at 0d47a67 (clean) with `runs/heartlined_0d47a67_od1` and
    `runs/heartline_0d47a67_ts2`; `~/v/scratch/t1line-wt` at 02562e4. Still off unless the captain names them: T1 line
    1's step back, the turn order of T1's third line, the T3 pathway tier, the closed T1 three-atlas line.
- **State at 2026-11-23, 20:18 UTC (superseded by the bullet above where they differ).** Decision log through
  2026-11-23 (w); next label (x). Ruling rows through R28; R29 is reserved for the second heart line (not placed yet).
  Everything else as the bullet below.
  - **THE CAPTAIN ORDERED THE TREND TRIED ON FILE D TOO** (own words "i want it too" under a side note; log (w): my
    reading = the named exception for D's command with the trend step at strength 0.5 and 1.0, two files, built before
    any score of D or of them; not an order that lets a score pick the base). The AI assistant is recused on the
    heart's scope and size and designs no line.
  - **RUNNING: workflow w3hq5inzx** (four agents, blind to test scores). Chain A: a rule writer
    (`lift8/heart_line_d/rule_text_draft.txt`) then a cold challenger (`lift8/heart_line_d/rule_text_final.txt`; a
    background waiter saves the challenger's hand-back to `lift8/heart_line_d/attack_final.json`). Chain B: an
    engineer on `campaign/heart-trend-line-d` (the option as `heart --entry organ_scope --trend-strength S --z-mode
    raw --captain-ruling R29`; files `T2heart_organ_ts05.h5ad`, `T2heart_organ_ts10.h5ad` and their halves; launcher
    `experiments/builds/heart_trend_line_d.sh <COMMIT> <RUN>` with identity steps for D and for the first line's TS05;
    scratch `~/v/scratch/heartlined/<RUN>`, output `runs/heartlined_<commit7>_<RUN>`) and a reviewer.
  - **NEXT, IN ORDER.** (1) When the challenger's final text is in (verdict place_as_final): `python3 -I -B
    lift8/place_r29.py <repo> <lift8> x <HH:MM>` places ruling row R29 (its first cell names the four strings the
    builder's check needs), the rules note's amendment with the final text verbatim, the pin and log (x); then the pin
    and ledger tests, commit, push. If the challenger put a question for the captain, pass it on unchanged with its
    default. If the verdict is do_not_place: no row; tell the captain why. (2) When the code review is ready: merge to
    main with the targeted tests, push, deliver to the DGX (`git push PrakashDGX_H2:v/vec-campaign
    <commit>:refs/heads/<name>`), check out in `~/v/scratch/h29b-wt` (idle, own environment), trace the launcher's
    refusal path under a throwaway HOME (forty zeros: expect a gate exit record), launch detached with PATH exported,
    wait with a background waiter. (3) Read exit files, the sha256 list and the stamps file only; copy the two halves
    to `runs/manual_upload/`, check sha256; record (build record, row R29's last cell, rules note and pin,
    `docs/UPLOAD_SHEET.md`, log) BEFORE any score; hand over with path, sha256, board and the Data-used text (released
    Task 2 heart stages and the unlocked E10.5 validation file only; no outside data).
- **State at 2026-11-23, 20:10 UTC (superseded by the bullet above where they differ).** Decision log through
  2026-11-23 (v); next label (w). Ruling rows through R28; next R29. Register through X020; next X021. T1's next row
  is T73 (the name T72 belongs to the MX2s row that was never placed). Tries unchanged: T1 463, embryo 88, heart 97,
  T3 74. NOTHING RUNS, here or on the DGX; no workflow is open.
  - **WAITING ON THE CAPTAIN'S UPLOADS.** `docs/UPLOAD_SHEET.md` lists 25 built files with sha256, condition and
    Data-used text (checked against disk at 19:55 UTC). No score has been reported since 21 Nov. When the captain
    pastes scores: the headline per board only into `results/manual_uploads.csv`, one log entry, each board's pre-set
    rule applied mechanically (T1 line 1: 2.5 only if the 2.0 file strictly beats B's; line 2: A20 only if A15
    strictly beats A's; line 3: Cm20 only if C strictly beats A's, Cm25 if Cm20 beats C's; T3: CK only if C2 and K23
    each strictly beat the current best, e4 only if K23 does; heart: TS10 only if TS05 strictly beats the full trim
    C's headline, whatever D scores; MSB, MS, D, embryo files: best-headline mark only), then the sheet updated.
  - **MX2s (T1 three-atlas line): CLOSED WITH NO TRY** (log (v); ruling row R27's last cell; register rows X019, X020;
    `docs/DATA_POLICY.md` note 2026-11-23 (d)). The plan step refused every sample of both new series on series-wide
    wording and two independent readers answered no on checklist points 2 and 5 for both. No data file was ever
    requested; row T72 never placed. Code on pushed, unmerged branches `campaign/t1-mx2s`, `-2`, `-3`, `-read`,
    `-read-2`. Scratch: `lift8/t72_row_cells_v3.json`, `lift8/t72_rule_v2/`, `lift8/t72_series_text_readers.json`,
    `lift8/t72_trial/`. A reopening needs the captain's own word and, for the series refused on a pooling word, a new
    named exception. DGX: `~/v/scratch/t1line-wt` at 02562e4 (clean), metadata in `~/v/scratch/t72_trial_plan` and
    `~/v/vec-ext/t1/meta` (kept).
  - **HEART: as the bullet below** (TS05 and TS10 handed over; DGX worktree `~/v/scratch/h29b-wt` at d41cd14 with
    `runs/heartline_d41cd14_ts1`).
  - **STILL OFF UNLESS THE CAPTAIN NAMES THEM:** T1 line 1's step back; the turn order of T1's third line; the T3
    pathway tier (stays out unless the captain says "build the T file"); MX2s (closed). Known, not yet fixed:
    `build_combination.sh`, `run_k22.sh`, `run_k23.sh` name the old T3 tables folder (fix before any new launch).
    About 390 agent worktrees sit under the repository (disks: / 37 GB, /mnt/md0 55 GB free; the background disk
    watcher still runs).
- **State at 2026-11-23, 20:03 UTC (superseded by the bullet above where they differ).** Decision log through
  2026-11-23 (u); next label (v). Ruling rows through R28; next R29. Nothing of ours runs on the DGX.
  - **HEART SEARCH LINE: BUILT AND HANDED OVER** (log (u); record
    `experiments/builds/results/heartline_d41cd14_ts1.txt`). `runs/manual_upload/T2heart_core_ts05_half.h5ad` sha256
    b6b8155be9b5469a2dbc8943a59df567f52c7cbe013bbf0b55e4da706fda27b7 (TS05: after file D);
    `T2heart_core_ts10_half.h5ad` sha256 cf53b0c0afe7d18b8f7c84d7595729e4fc46e15e6077e645b6dc0a2fbfed39af (TS10: only
    if TS05's headline is strictly above the full trim C's, whatever D scores; the line ends there).
    `docs/UPLOAD_SHEET.md` lists 25 files (checked against disk). WHEN A HEART SCORE COMES BACK: headline only into
    `results/manual_uploads.csv`; TS05 strictly above C's headline -> TS10 is next; otherwise the line ends and TS10
    is never uploaded. Declared limit: the builder does not refuse a hand-typed `--seed` with `--trend-strength` (the
    launcher passes none).
  - **MX2s (T1): THE PLAN STEP REFUSES BOTH NEW ATLASES ON SERIES-WIDE WORDING.** A trial of the plan step at the
    data-side branch's head (`campaign/t1-mx2s-3` 02562e4, DGX worktree `~/v/scratch/t1line-wt`, output in the
    throwaway folder `~/v/scratch/t72_trial_plan`, metadata only) ran clean on the real records (the two first-run
    faults are fixed) and marks every sample UNKNOWN: GSE230531 (X020) for a form of 'pool' in its series-wide
    extraction protocol (the exception of ruling row R27 says series text 'still refuses on pooling words': final
    under that wording); GSE244783 (X019) for the token 'mix' of the reviewers' WIDER word list (combine, mix,
    multiplex, superload, aggregate), which may be a reagent mix. Two independent screening readers (workflow
    whkom36i7; files in scratch `lift8/t72_trial/`, extracted by code, unread by the main session) classify the
    wording and answer checklist points 2 and 5 on the series-wide text. IF X019 fails either point with either
    reader: both atlases are out, row T72 is NOT placed, the line closes with no try: log it and tell the captain. IF
    X019 clears: a score-blind rule agent and a cold challenger decide whether hits of the wider list alone may be
    cleared by both readers pinning the field's sha256 (forms of pool and merge stay hard refusals); only then the
    engineering round (below).
  - **MX2s ENGINEERING STILL OWED IF THE LINE LIVES.** Data-side review (scratch `lift8/mx2s_data_review.json`): not
    ready, four blockers (the pooling stage condition over the record's whole series text; the imports test's scan
    restored with a literal allow-list; raw files after reduction: a README departure line; one test renamed and the
    real one added). Read-script review (`lift8/mx2s_read_build.json`): three blockers; two are being fixed in
    workflow w657mwmml (`campaign/t1-mx2s-read-2`), the third is the interface, closed at integration (one branch from
    both, an integration test, the launcher's data-step command lines). Row T72 v3 in scratch (sha256 f8302ec2...).
    Before placement also: run the launcher's test list on the DGX at the integration commit (every launcher step runs
    once: a tests-step failure would end the row).
- **State at 2026-11-23, 19:17 UTC (superseded by the bullet above where they differ).** Decision log through
  2026-11-23 (t); next label (u). Ruling rows through R28; next R29. Main at the commit of this note; nothing of ours
  runs on the DGX.
  - **CORRECTION to the heart route below:** the DGX worktree `~/v/scratch/h29b-wt` belongs to `~/v/vec-campaign` (not
    `-d`): deliver heart commits with `git push PrakashDGX_H2:v/vec-campaign <commit>:refs/heads/<name>`.
    `~/v/scratch/t1line-wt` belongs to `~/v/vec-campaign-d` (the T1 atlas line).
  - **ROW T72 IS NOW v3 (scratch `lift8/t72_row_cells_v3.json`, as one ledger line sha256
    f8302ec25e2fd06e08a87bef994ac4fbde794c6fc197a70269be0fb3da8dee14, also in `lift8/t72_row_v3.sha256`):** the rule
    agent settled the six open points (`lift8/t72_rule_v2/open_points_settled.md`, hand-back
    `lift8/t72_open_points.json`); its eight exact old-to-new wordings are applied by code, and the headroom check the
    fifth cell asks of the main session is done and written into that cell (the ledger's tables, rerun at c4501f5).
    This hash goes in the two readers' brief; after the screens start only the fills of cells 3, 4 and 10 may change.
    Code consequences still to hand to the engineers in a fix round: the pre-gate fails when X002's own CM count is
    below 200 at either end; the launcher's clean-tree check exempts only untracked files under
    `experiments/t1_m1/results/t72_mx2s/`; start records are committed with the task files after the read; cell 3
    names the screened files' sha256 as a placement fill.
  - **THE HEART SEARCH LINE IS IN FORCE** (the captain's words "start the heart search", ruling row R28, the rules
    note's amendment of 2026-11-23, log (t)). Two files on the full trim C's command with the trend step at strength
    0.5 (TS05) and 1.0 (TS10), nothing else changed. Rule: TS05 is uploaded after the queued file D; TS10 only if
    TS05's headline is strictly above C's (C's headline is on record in `results/manual_uploads.csv`: the line's
    base); the line ends at TS10. The AI assistant is recused on the heart's scope and size and states no preference.
    T1 line 1's step back and the third T1 line's turn order were NOT named and stay off.
  - **HEART, NEXT.** Workflow wai2reyjh codes it (`campaign/heart-trend-line`: builder option `--trend-strength
    {0.5,1.0}` only with `--entry core_trim --z-mode raw --captain-ruling R28`; launcher
    `experiments/builds/heart_trend_line.sh <COMMIT> <RUN>`; a reviewer with the blocker standard up front;
    specification in scratch `lift8/heart_line_final.txt`). On a ready review: merge to main with the targeted tests;
    deliver the commit to the DGX (`git push PrakashDGX_H2:v/vec-campaign-d <commit>:refs/heads/<name>`), check it out
    in a detached worktree with its own environment (file D was built in `~/v/scratch/h29b-wt`; nothing runs there),
    trace the launcher's refusal path under a throwaway HOME, launch detached with PATH exported, wait with a
    background waiter. Read exit files, sha256 and the stamps' commit and clean fields only. Then: build record under
    `experiments/builds/results/`, the two halves copied to `runs/manual_upload/` and checked, their sha256 into the
    log and the rules note BEFORE any score, `docs/UPLOAD_SHEET.md` updated (Data used: released challenge data only;
    no outside data), and the files handed to the captain with path, sha256, board and that text.
  - **MX2s (T1), RUNNING.** Workflow wai2reyjh also holds the data side (`campaign/t1-mx2s-3` cut from
    `campaign/t1-mx2s-2` 75356ec: the three pooling blockers of the re-review in scratch
    `lift8/mx2s_pooling_fix.json`, the two faults of the first real plan run (a 403 answered once by the repository:
    bounded retry; the parser's refusal of the standard field name with a slash), and engineer_notes.md section A)
    with a reviewer. Workflow w101ygojc holds the practice-read script (`campaign/t1-mx2s-read`) with a reviewer and
    the rule agent on the six open points (`lift8/t72_rule_v2/open_points_settled.md`). Then as the bullet below says:
    merge, the plan step on the DGX, the screens' files, two readers, pins; placement of row T72 only if at least one
    new atlas survives.
- **State at 2026-11-23, 19:11 UTC (superseded by the bullet above where they differ).** Decision log through
  2026-11-23 (s); next label (t). MX2s in progress; everything else as the bullets below say.
  - **MX2s STATE.** Row T72 v2 (NOT placed): scratch `lift8/t72_row_cells_v2.json` (row sha256 as a ledger line
    0e4ee2d1... before the placement fills), rule files `lift8/t72_rule_v2/` (points_settled.md, engineer_notes.md,
    readers_checklist.md; open_points_settled.md when the rule agent hands back), the attack
    `lift8/t72_rule/attack_faults.md`. Port: `campaign/t1-mx2s` 290da3f (hand-backs `lift8/mx2s_port_build.json`,
    `mx2s_port_review.json`). RUNNING: w7xf07u2a (pooling fix on `campaign/t1-mx2s-2`, then a check); w101ygojc (read
    script on `campaign/t1-mx2s-read` + reviewer; the open-points rule agent).
  - **MX2s NEXT.** (1) When the pooling fix is ready: launch the DATA-SIDE engineer + reviewer on a branch cut from
    `campaign/t1-mx2s-2`: the fetch retries a 403 or 5xx a bounded number of times with a pause (one host); the SOFT
    parser accepts the repository's standard field names that hold a slash (`Sample_contact_zip/postal_code` refused
    the whole of GSE230531); engineer_notes.md section A in full (per-field sha256 at the plan step for every sample
    record; stage stated by title AND characteristics; `fetch.py fetch-meta` into `~/v/scratch/t72_screen`; member
    tables without consensus variants; the labeller with the skm and derived classes; the guards; outputs
    `data/external/t1_ext_t72/`, `data/derived/t1_ext_t72/`). (2) Merge; DGX worktree `~/v/scratch/t1line-wt`
    (currently at 290da3f; the first plan left `~/v/vec-ext/t1` with one metadata file: the revised row wants the
    screens under `~/v/scratch/t72_screen`); plan; fetch-meta; copy the records and the per-field hashes to local
    scratch unread; TWO readers with `readers_checklist.md` and the row's sha256; pins through `pin-template`; freeze
    the yaml. If both atlases drop: do not place the row; log; tell the captain. (3) Otherwise: merge the read script;
    placement commit (row with its fills, pins, yaml, derived classes, X002/X019/X020 amendments naming row T72);
    launcher trace; run; OUTCOME line only.
- **State at 2026-11-23, 18:40 UTC (superseded by the bullet above where they differ).** Decision log through
  2026-11-23 (r); next label (s). Next ruling row R28; next register row X021; T1's next row T72; T3's next K28. Local
  main = origin/main at ab7444d or one notes commit after it. `docs/UPLOAD_SHEET.md` lists every waiting file with
  sha256, condition and Data-used text (no score in it; update it when a file is added or a score changes a
  condition).
  - **THE CAPTAIN: "You have my approval and word for all the above" (entry (r)).** Read as: (1) the word T1 candidate
    MX2s needed, with the stage-parser exception: ruling row R27, `docs/DATA_POLICY.md` note 2026-11-23 (c); practice
    read only, a file needs a further question after AHEAD. (2) The T3 pathway tier STAYS OUT (a general approval is
    not an order to build on an item under the organisers' genotype clause; the captain was told to say 'build the T
    file' if that was meant). Nothing else changes.
  - **MX2s IS RUNNING.** Register rows X019 (GEO GSE244783) and X020 (GEO GSE230531) are placed; X002 amended. Texts
    for agents in scratch: `lift8/drastic_spec_t1_1.txt` (governs), `drastic_t1_challenger_notes.txt`,
    `drastic_t1_designer_mx2.txt`. Workflow w7wcujlu0 (journal folder wf_3689f154-f29): chain 1 rule writer -> cold
    attack -> reviser (files under `lift8/t72_rule/` and `lift8/t72_rule_v2/`); chain 2 engineer (branch
    `campaign/t1-mx2s`: the parked `experiments/t1_ext` on main, strict form, parser exception, screens, the known
    blocker) -> reviewer. Waiter bktlrykz5 saves the port's review to `lift8/mx2s_port_review.json`.
  - **MX2s, NEXT STEPS IN ORDER.** (a) Port review ready: merge `campaign/t1-mx2s`, deliver to `v/vec-campaign-d`,
    make a detached worktree with its own env; run the plan step on the DGX (family SOFT metadata only); copy the
    plan's records to local scratch WITHOUT reading them; two independent readers apply
    `lift8/t72_rule_v2/readers_checklist.md` (or v1) and hand back accessions, record sha256 and yes/no per check;
    `pin-template`; commit pins; plan again; `groups.py derive` (commit its output and sha256); fetch barcodes and
    features first, the screens (a) to (d), then matrices; add file names and sha256 to X019/X020. (b) Launch the
    engineers for the consensus tables (the specification's steps 4 to 6, in `experiments/t1_ext`) and for the
    practice read (steps 7 and 8: a generalised `cmjax_read.py`, its launcher, tests), each with a reviewer, blocker
    standard up front, briefs pointing at the rule files. (c) Place row T72 from the revised cells (factual updates
    only at placement, declared in the log) BEFORE the reducer runs; trace the launcher; run in a DGX worktree with
    PATH exported; read the OUTCOME line only. (d) AHEAD: one priced question to the captain (the file, the
    X002/X019/X020 upload amendments); then the builder option, review, two builds, all by 25 Nov 12:00 UTC. NOT AHEAD
    or a failed pre-gate: the line closes.
- **State at 2026-11-23, 18:14 UTC (superseded by the bullet above where they differ).** Decision log through
  2026-11-23 (q); next label (r). Next ruling row R27; next register row X019; T3's next row K28. Local main =
  origin/main at 56ff9c4 or one notes commit after it. Nothing is running on the DGX and no workflow is running.
  Tries: T1 463, embryo 88, heart 97, T3 74 (unchanged today).
  - **T3 ROW K27 (CHIPQ-S) IS CLOSED: "Try K27: STOP at P1: the database bars are not met; no fold read, no try, no
    file."** (entry (q)): 5 kept experiments of the test's gene, a bound set of 1 panel gene against a bar of 15.
    Nothing counted; no table; listing empty, pins None; no second commit. A counts-only check on the 421 placebo
    antigens showed the code path works (bound sets 0 to 423, median 5). Not to be rerun under another row; the
    database is not read again on this fold. The build path was never started. DGX leftovers, all harmless:
    `~/v/vec-campaign-c` sits at 4b11323 with the untracked folder `experiments/t3_decomp/results/k27_chip_query/`
    (the two reading files, committed locally at 56ff9c4: move it aside before the next checkout there and compare);
    `~/v/vec-ext/t3/chipatlas` (354 MB, kept for the record); `~/v/vec-ext/t3/tables_k27` (empty); worktree
    `~/v/scratch/k27-wt`.
  - **THE OUTSIDE-DATA ROUND OF 23 NOV IS FINISHED FOR T3 AND T2.** What is left of it: T1's MX2s, which needs the
    captain's word before any code (default no). By today's measured pace (a candidate went from brief to gate in
    about three and a half hours with four review rounds) it would still fit if ordered by about 24 Nov 12:00 UTC;
    that is the AI assistant's estimate, not the challenger's (who said the evening of 23 Nov).
  - **READY FOR UPLOAD, ALL IN `runs/manual_upload/` (the captain's own acts; at most two a board a day; no score has
    been reported since 21 Nov).** T1: `T1_s2rcC.h5ad` (31477ef7...), `T1_s2rcB2.h5ad` (the 2.0 file, 52897fe3...),
    `T1_s2rcMSB.h5ad` (78bb17ee..., in the plain momentum file's place by the captain's word), `T1_s2rc_MS.h5ad`
    (6103d105...), then `T1_s2rcA15.h5ad` and the lines' next files on their conditions; `T1_s2none.h5ad` (F) cannot
    raise the mark. T3: `T3_c1_s2.h5ad`, `T3_c2_s2.h5ad`, `T3_k23_s2.h5ad`, then CK and the E = 4.0 file on their
    conditions, `T3_refused_s1.h5ad`. Embryo: `T2embryo_s2x_cr.h5ad` (C), `T2embryo_s2x_crs.h5ad`,
    `T2embryo_s2xs_ov.h5ad` (SO), `T2embryo_s2.h5ad`, `T2embryo_s2x_ovhi.h5ad` (O_hi). Heart:
    `T2heart_organ_half.h5ad` (D). WHEN SCORES ARE PASTED: the headline per board only, then each pre-set rule (see
    the earlier bullet 'WHEN SCORES ARE PASTED'); MSB and MS serve the best-headline mark only.
  - **OPEN WITH THE CAPTAIN (each default no change):** the T3 pathway tier (stays out; entries (e) to (g)); MX2s
    (stays parked).
- **State at 2026-11-23, 17:08 UTC (superseded by the bullet above where they differ).** Decision log through
  2026-11-23 (p); next label (q). Local main = origin/main at e69d0e1 or one notes commit after it.
  - **T1: THE COMBINED FILE IS BUILT AND HANDED OVER (entry (p)).** `runs/manual_upload/T1_s2rcMSB.h5ad`, sha256
    78bb17eea5721e590cfd2c9f707e9cd6f2f496e081acede4f697dc4035871da2 (set t1msc, run m1, a3e879a; record
    `experiments/builds/results/t1msc_a3e879a_m1.txt`). Place by the captain's "try it": the plain momentum file's
    place unless the captain names another; MS (`T1_s2rc_MS.h5ad`, 6103d105...) follows it. Neither has a score. When
    T1 scores are pasted: headline only; MSB and MS serve the best-headline mark only and enter no line.
  - **T3: ROW K27 IS PLACED, REGISTRATION COMMIT REG = 4b11323e8da674804e944d4d16e770a095284cd1 (entry (o)).**
    `~/v/vec-campaign-c` on the DGX is checked out at REG (clean; two untracked copies of K24's reading moved aside to
    `~/v/scratch/k24_results_dgx_copy`, identical to the committed ones); `~/v/vec-ext/t3/tables_k27` exists (empty);
    both launcher invocations traced under a throwaway HOME (exit 97, all.exit written). **WAITING ON** the
    target-file download (second run, in `~/v/scratch/k27-wt` at 66d2e21; waiter bbxct9705; exit record
    `~/v/scratch/k27_targets2.exit`; its two count lines name no antigen). ON exit 0: run `chipq_fetch.py verify`
    there (exit code only), then in `~/v/vec-campaign-c`: `export PATH=$HOME/.local/bin:$PATH UV_FROZEN=1
    UV_NO_SYNC=1; setsid nohup bash experiments/t3_decomp/run_k27.sh db 4b11323e8da674804e944d4d16e770a095284cd1`
    detached; waiter on `~/v/scratch/k27_4b11323/db/all.exit` (0 P2 passed; 28 P1 stop; 29 P2 stop; 22 refused). P1's
    two counts are in `runs/k27_4b11323/p1.out`. On 0: read `~/v/scratch/k27_4b11323/db/second_commit.txt` (sha256 and
    table names only), make the SECOND commit on a side branch cut at REG (`campaign/k27-run`): only
    `experiments/t3_ext/TABLES_K27.sha256` and the pin lines of `experiments/t3_decomp/k27_chip_query.py`; deliver;
    check out SECOND in `~/v/vec-campaign-c`; `run_k27.sh run <SECOND> 4b11323e8da674804e944d4d16e770a095284cd1`;
    waiter on `.../run/all.exit`; read only `grep '^OUTCOME'` of
    `experiments/t3_decomp/results/k27_chip_query/reading.txt`. If P1 passes, launch the build-path engineer and
    reviewer at once. Add the download's counts and the manifest's state to register row X018. Count the tries in the
    ledger's totals when P3 starts.
  - **OPEN WITH THE CAPTAIN (each default no change):** the T3 pathway tier (stays out); MX2s (stays parked).
- **State at 2026-11-23, 16:53 UTC (superseded by the bullet above where they differ).** Decision log through
  2026-11-23 (n); next label (o). Local main = origin/main at f2a0bdf or one notes commit after it. The combined MS
  build (run m1, waiter bcqv7ks8j) is still running; see the bullet below for what to do on done=0.
  - **CHIPQ-S, WHERE IT STANDS.** Database code merged (f3cc4d7); lists and rows fetched on the DGX; the screen done
    (two readers agree: 5 of 20 Ctnnb1 experiments kept, 7 of 32 Gata4); kept list pinned
    (`experiments/t3_ext/chipq_kept.yaml`, f2a0bdf); X018 amended. The `targets` download stopped at the fifth file:
    the archive answers a placebo antigen it has no file for with a 302 to another host. **RUNNING:** wrrab9haf
    (branch `campaign/t3-chipq-5`: the fetch records such placebo antigens in `screen/placebo_unfetched.json` and goes
    on; the pool leaves them out; reviewer follows); wgt8p6ct5 (branch `campaign/t3-k27-read-2`: the read script's
    tests made independent of the committed pin state, the one blocker of its review; reviewer follows). A safety
    classifier stopped a rule agent's answer on six small points (and very likely the first rule writer's): not
    pressed; they stand as the reviewed code has them (entry (n)).
  - **THE LAUNCH PLAN FOR ROW K27** (hand-backs in scratch: `lift8/k27_read_build.json`; row cells
    `lift8/k27_row_cells_v2.json`; rule files `lift8/k27_rule_v2/`). (1) Merge `campaign/t3-chipq-5` when ready;
    deliver; in `~/v/scratch/k27-wt` check out the merged commit and rerun `chipq_fetch.py targets` (it takes up from
    the four files), then `verify`. (2) Merge `campaign/t3-k27-read-2` when ready (new files only). (3) REGISTRATION
    COMMIT on main: place row K27 from the v2 cells (after T71's neighbours, as the last K row), commit the rule files
    for the record, the listing empty, the pins None. Departures of the read engineer to accept in the log:
    `run_k27.sh db REG` then `run_k27.sh run SECOND REG`; the listing empty at REG; listing lines `<sha256> <name>`
    without .csv; P2 loads the barred factor for the gene order only. (4) On the DGX: `mkdir -p
    ~/v/vec-ext/t3/tables_k27`; the run is in `~/v/vec-campaign-c` ITSELF (the row names it; K23's run folder
    `runs/k23_b706e82` is there; nothing else runs there; move aside untracked files that the new commit tracks, then
    check out REG); trace a refusal path under a throwaway HOME; `bash experiments/t3_decomp/run_k27.sh db <REG>`
    detached with PATH exported; waiter on `~/v/scratch/k27_<REG7>/db/all.exit` (0 ok; 28 P1 stop; 29 P2 stop). (5)
    SECOND COMMIT: ONLY the listing and the script's pin lines, from `~/v/scratch/k27_<REG7>/db/second_commit.txt`;
    make it on a SIDE BRANCH cut at REG (invocation 2 refuses any other difference between REG and SECOND, and main
    will have moved); deliver; `run_k27.sh run <SECOND> <REG>`; waiter on `.../run/all.exit`. (6) Read only the
    OUTCOME line (`grep '^OUTCOME'` of reading.txt); a PASS goes to the captain as one question (section 4.4's
    provenance rule, A15, share-alike) and the build path is then written; P1 passing is the trigger to start the
    build-path engineer.
- **State at 2026-11-23, 16:32 UTC (superseded by the bullet above where they differ).** Decision log through
  2026-11-23 (m); next label (n). Local main = origin/main at a3e879a or one notes commit after it.
  - **THE CAPTAIN: "[try it]" (entry (m)).** The combined file `T1_s2rcMSB` leaves the end of the T1 queue: unless the
    captain names another place it takes the plain momentum file's place (after C and the 2.0 file) and MS follows it.
    **ITS BUILD IS RUNNING:** run m1 at a3e879a in `~/v/scratch/t1line-wt` since 16:31 UTC (launcher
    `experiments/builds/t1_ms_combined.sh`; scratch `~/v/scratch/t1msc/m1`; output `runs/t1msc_a3e879a_m1`; steps
    gate, tests, D, A, B, C, X2, MS, MSBa, MSBb, done; waiter bcqv7ks8j). The builder was reviewed ready with no
    blocker; the refusal path was traced on the DGX (exit 96). ON done=0: read the sha256 list and stamps only; copy
    `MSBa/T1_s2rcMSB.h5ad` to `runs/manual_upload/`; write `experiments/builds/results/t1msc_a3e879a_m1.*`; add the
    sha256 to the log (the list's fifth condition: before its score), to row R26's last cell and a rules note
    (re-pin); hand over with path, sha256, board, Data used (as R13's file).
  - **CHIPQ-S.** The two lists are fetched on the DGX (`~/v/vec-ext/t3/chipatlas`, 16:21 UTC; hashes in scratch
    `lift8/chipq_lists_manifest.txt`, to go into register row X018) from the worktree `~/v/scratch/k27-wt`
    (vec-campaign-c's repo, at 145f5d3; run steps with `~/v/vec-campaign-c/.venv/bin/python -I -B`). The `rows` step
    waits for the final keyword filter: workflow w6mzzh0vc applies the last three pattern fixes on
    `campaign/t3-chipq-4` (cut from -3 = the rule agent's addendum on -2) and re-checks. Then: merge -4, deliver,
    `rows`, copy the two query rows files to scratch unread, two screening agents, kept list, `targets`, `verify`.
    Four one-line confirmations for a rule agent before P2: scratch `lift8/k27_rule_open_confirmations.txt`. The read
    script: workflow wg1onofkg (`campaign/t3-k27-read`). The corrected row: `lift8/k27_row_cells_v2.json`; settled
    points and engineer notes: `lift8/k27_rule_v2/`. The row names `~/v/vec-campaign-c` itself as the run checkout
    (K23's run folder is there).
  - **OPEN WITH THE CAPTAIN (each default no change):** the T3 pathway tier (stays out); MX2s (stays parked).
- **State at 2026-11-23, 16:11 UTC (superseded by the bullet above where they differ).** Decision log through
  2026-11-23 (k); next label (l). Next ruling row R27; next register row X019; T3's next row K27 (not placed). Local
  main = origin/main at 2646249 or one notes commit after it.
  - **THE CAPTAIN: "combine MS" (entries (j), (k)).** The combined MS list is in force since 16:10 UTC (row R26, rules
    note re-pinned), before any score of MS was reported: ONE file, `T1_s2rcMSB` = `build_test.py t1 --g-mom --base
    lib3 --feng-groups recomputed --feng-fraction-mult 1.5 --direction-momentum all --captain-ruling R26`, queued
    LAST, reading no score; the captain may move it (for example into MS's slot). Workflow wt7lfrgsp (journal folder
    wf_770e75f9-ca7) is at its engineer (branch `campaign/t1-ms-combined`, launcher
    `experiments/builds/t1_ms_combined.sh <COMMIT> <RUN>`, scratch `~/v/scratch/t1msc/<RUN>`) and reviewer. On a ready
    review: merge, tests, deliver to `v/vec-campaign-d`, check out in `~/v/scratch/t1line-wt`, dry-run trace under a
    throwaway HOME, launch with PATH exported, waiter; on done 0: copy the file to
    `runs/manual_upload/T1_s2rcMSB.h5ad`, record `experiments/builds/results/t1msc_<c7>_<RUN>.*`, add the sha256 to
    the log BEFORE its score, hand over. It drops if not built by 25 Nov 12:00 UTC.
  - **CHIPQ-S (T3).** Row K27 as written: scratch `lift8/k27_row_cells.json`; the cold attack
    (`lift8/k27_rule_attack.json`): place after 13 corrections, two-step placement (registration commit before P1;
    second commit after P2 fills the listing and pins; counted steps at the second). Workflow wj1jeljlg (reviser)
    writes the corrected row to `lift8/k27_rule_v2/`. Workflow wtz3mwdnj (journal wf_98fdb744-4cd): the fix engineer
    for the five database blockers (branch `campaign/t3-chipq-2`) then a reviewer; waiter b3f2wjjui saves the fix
    engineer's hand-back to `lift8/chipq_data_fix.json`. NEXT: when the reviser and the fix engineer are both in,
    launch the read engineer (branch cut from `campaign/t3-chipq-2`; files `experiments/t3_decomp/k27_chip_query.py`,
    `run_k27.sh` with TWO invocations, tests; chipq files untouched) and a reviewer; when the fix review is ready:
    merge, DGX fetch (lists, rows), two screening agents, kept list, targets. The registration commit needs the read
    script merged; then P1, P2, pins commit, P3 to the read. The build path is coded only if P1 passes.
  - **OPEN WITH THE CAPTAIN (each default no change):** the T3 pathway tier (stays out); MX2s (stays parked; had to
    start on the evening of 23 Nov); unless named: T1 line 1's step back, the heart search line's exception, line 3's
    turn order.
- **State at 2026-11-23, 15:45 UTC (superseded by the bullet above where they differ).** Decision log through
  2026-11-23 (i); next label (j). CHIPQ-S step 2 is running: workflow w18fza0f9 (journal folder wf_916677f6-312): rule
  writer 2 (files under scratch `lift8/k27_rule/`: points_settled.md, cell_01..11.txt, engineer_notes.md), data-path
  reviewer, then read engineer (branch `campaign/t3-chipq-read` cut from `campaign/t3-chipq` 4498ece) and read
  reviewer. Waiters: b8vke5ffw saves the data reviewer's hand-back to `lift8/chipq_data_review.json`; bakunpdcx saves
  the rule writer's to `lift8/k27_rule_handback.json`.
  - **ON THE DATA REVIEW (ready):** merge `campaign/t3-chipq` to main (tests: tests/test_chipq.py,
    test_k24_wider_tables.py, test_build_k24_form.py, the rules pin), deliver to `v/vec-campaign-c`, check out in
    `~/v/scratch/k26-wt` (detached, own env), run `python -I -B experiments/t3_ext/chipq_fetch.py lists` then `rows`
    (PATH exported; prints counts only), copy `~/v/vec-ext/t3/chipatlas/screen/Ctnnb1.rows.tsv` and `Gata4.rows.tsv`
    to local scratch `lift8/chipq_screen/` WITHOUT reading them, launch two independent screening agents with the
    specification's SCREEN paragraph verbatim (they hand back kept ids and drop counts per class only), intersect by
    code, commit `experiments/t3_ext/chipq_kept.yaml` + `CHIPQ_KEPT.sha256`, add the fetch date and hashes to register
    row X018, run `targets` and `verify`. Not ready: fix round on the branch first.
  - **ON THE RULE (eleven cells):** check the cells mechanically (11, no raw pipe, first cell opens 'K27 (one row, one
    script;'), wait for the read reviewer's attack on the rule before placing the row; place row K27 before P1; then
    P1 (`chipq.py p1`: two counts; exit 0 pass, 3 fail = line closed, no try). If P1 passes: launch the build-path
    engineer and reviewer; build tables, commit `TABLES_K27.sha256`; run `run_k27.sh` on the DGX with a dry-run trace
    first.
- **State at 2026-11-23, 15:40 UTC (superseded by the bullet above where they differ).** Decision log through
  2026-11-23 (h); next label (i). Local main = origin/main at 56a357c or one notes commit after it. Next ruling row
  R26; next register row X019; T3's next row K27.
  - **MS IS BUILT AND HANDED OVER (entry (h)).** `runs/manual_upload/T1_s2rc_MS.h5ad`, sha256
    6103d10535665634ee0cad3390590a5df132d1805682fb1d309bb19581fdaa5e (set t1ms, run r1, commit 6bc36d9; record
    `experiments/builds/results/t1ms_6bc36d9_r1.txt`). Place: after the files of entry 2026-11-20 (r)'s list (C, and
    the 2.0 file by the order of 2026-11-21 (y)), the captain may reorder; never stacked on A, B or C. No score of it
    exists. The DGX worktree `~/v/scratch/t1line-wt` now sits at 6bc36d9.
  - **THE T3 PATHWAY TIER IS WEIGHED (entries (e) to (g)); ONE QUESTION IS WITH THE CAPTAIN: should the tier stay out
    of the Task 3 file? Default: it stays out.** Verdict of three cold analysts: leave it out (P(mistake) if used
    0.80, if left out 0.12); the screened search found no organiser statement beyond the three on record. A no from
    the captain means a file from the T knockout GSE205117 only (0.35 odds of being counted held out; about 0.4 to be
    built in time; two new readers must pin it; a reviewed commit must cite the captain's ruling to move it off the
    held list). The AI assistant has told the captain it would not build on GSE165300 in any case (the plan has
    classed it as held out since October). Nothing of the tier is fetched or unheld without the captain's no.
  - **RUNNING.** wh5dbcldb: CHIPQ-S rule writer then data-path engineer (see the bullet below for what follows its
    hand-back). **OPEN WITH THE CAPTAIN:** the pathway-tier question; MX2s (reopen the parked T1 multi-atlas branch in
    its strict form with a stage-parser exception; default no; must start tonight); unless named: T1 line 1's step
    back, the heart search line's exception, line 3's turn order.
- **State at 2026-11-23, 15:16 UTC (superseded by the bullet above where they differ).** Decision log through
  2026-11-23 (e); next label (f). Next ruling row R26; next register row X019; T3's next row K27. The outside-data
  round handed back (entry (d)): T3 candidate CHIPQ-S kept (register row X018 placed), T1 candidate MX2s kept but
  waiting on the captain's word, nothing for the T2 boards. The captain asked for the held T3 pathway tier to be
  weighed (entry (e)).
  - **RUNNING (four workflows).** wdfh9z62z: MS's builder (see the bullet below). wh5dbcldb: CHIPQ-S rule writer, then
    the data-path engineer (branch `campaign/t3-chipq`). On its hand-back: save to lift8/, launch (a) a reviewer of
    the data path and of the rule's design (blocker standard up front) and (b) the read's engineer and reviewer
    (`k27_chip_query.py`, `run_k27.sh`, builder path; branch cut from `campaign/t3-chipq`); merge the data path when
    ready; deliver to `v/vec-campaign-c`, check out in `~/v/scratch/k26-wt`; fetch lists on the DGX (both hosts answer
    200 from there); copy the two query antigens' row files to local scratch WITHOUT reading them; two independent
    screening agents read them and hand back ids and counts only; commit `experiments/t3_ext/chipq_kept.yaml` and
    `CHIPQ_KEPT.sha256`; fetch targets; place row K27 (from the rule writer's cells) before P1; P1 (uncounted, two
    counts); tables; pins; the counted read. A pass goes to the captain as one question with the file (provenance rule
    of section 4.4, A15, share-alike). wdz9jjg9o: the cold weighing of the T3 pathway tier; on its hand-back put the
    challenger's ONE priced question to the captain with the weighing and the AI assistant's own compliance view;
    nothing of the tier is fetched before the captain's answer.
  - **OPEN WITH THE CAPTAIN.** (1) MX2s: reopen the parked T1 multi-atlas branch in its strict form, including a
    stage-parser exception that loosens a team guard (default no; about 20 to 24 h; odds about 0.12; a few tenths to a
    point on practice). (2) The T3 pathway tier, after the weighing. (3) Unless named: T1 line 1's step back, the
    heart search line's exception, line 3's turn order.
- **State at 2026-11-23, 14:48 UTC (superseded by the bullet above where they differ).** Decision log through
  2026-11-23 (c); the next label is 2026-11-23 (d). Local main = origin/main at ba05554 or one notes commit after it.
  **THE CAPTAIN: "You have my approval for all" (entry (c)).** Read as the yes to the four items the message before it
  put by name: MS (row R23), SO and O_hi (row R24; O_hi's naming lifts the embryo's depth line for that one file), the
  heart's file D for upload (row R25); rules amended and re-pinned. NOT read as covering T1 line 1's step back, the
  heart search line's exception or line 3's turn order: those stay as they are unless the captain names them (said to
  the captain).
  - **HANDED OVER FOR UPLOAD ON 23 NOV (all in `runs/manual_upload/`, sha256 checked on disk today).** Embryo: SO
    `T2embryo_s2xs_ov.h5ad` (77062082...), O_hi `T2embryo_s2x_ovhi.h5ad` (7e4d1eaa...), beside C
    `T2embryo_s2x_cr.h5ad`, the sized file `T2embryo_s2x_crs.h5ad` and `T2embryo_s2`. Heart: D
    `T2heart_organ_half.h5ad` (1b4342c1...). T1 (my advice for the 200 goal: C and the 2.0 file before F):
    `T1_s2rcC.h5ad`, `T1_s2rcB2.h5ad`, then `T1_s2rcA15.h5ad`. T3: `T3_c1_s2.h5ad`, `T3_c2_s2.h5ad`, then
    `T3_k23_s2.h5ad`.
  - **RUNNING.** wdfh9z62z: MS's builder (`campaign/t1-ms`; option `--direction-momentum all --captain-ruling R23`;
    launcher `experiments/builds/t1_ms.sh <COMMIT> <RUN>`, steps gate, tests, D, MSa, MSb, done). On a ready review:
    merge (keep both sides of the tests list), tests, commit, deliver to `v/vec-campaign-d`, check out in
    `~/v/scratch/t1line-wt`, dry-run trace, launch with PATH exported, waiter on every step; on done 0 copy
    `MSa/T1_s2rc_MS.h5ad` to `runs/manual_upload/`, record `experiments/builds/results/t1ms_<c7>_<RUN>.*`, fill row
    R23's last cell by a log entry, hand over (place: after entry 2026-11-20 (r)'s list; Data used as R13's file).
    wprn6jz7a: the drastic outside-data design round (see the bullet below).
  - **WHEN SCORES ARE PASTED:** headline per board only; then each pre-set rule (T1 line 1: 2.5 if the 2.0 file
    strictly beats 71.59; line 2: A20 if A15 strictly beats A's 70.8; line 3: Cm20 if C strictly beats A's 70.8, then
    Cm25; T3: CK if C2 and K23 each strictly beat 62.30, the E = 4.0 file if K23 strictly beats 62.30; embryo and
    heart: best-headline mark only).
- **State at 2026-11-23, 14:40 UTC (superseded by the bullet above where they differ).** Decision log through
  2026-11-23 (b); the next label is 2026-11-23 (c). **THE CAPTAIN'S INSTRUCTION OF 23 NOV (entry (b)):** "for the 2
  submissions of today, I will only submit it if it reaches the score more than 200. If you really believe that it
  doesn't, then improve it drastically and utilize external datasets extensively [Even if it takes 1TB of space in the
  dgx]". My answer, given plainly: I do not believe today's ready files reach 200 (standing 197.04; rough odds one in
  25 today, one in 13 by the deadline, heart left out); an upload cannot lower the standing and the search lines need
  uploads. It answers none of the six open questions.
  - **RUNNING: workflow wprn6jz7a (drastic-outside-data-round; hand-back to `lift8/`).** Six score-blind designers and
    three cold challengers; it returns final specifications for at most two T3, two T1 and one T2 candidate. For each
    kept candidate: engineer + reviewer workflow (blocker standard up front; worktree isolation; test rules), register
    rows first, ledger row before any read, DGX worktree with its own environment and PATH exported, practice read,
    then the test-board build twice; a file from it goes to the captain as one priced question (T2 files with an
    outside atlas need A15's question; closed levers need the captain's named exception). Disk cap for the round:
    about 600 GB under `/home/prachh/v` (check `df` before each download; the root disk is full). A built file must
    exist by 2026-11-25 12:00 UTC.
  - **An allowed class found unused (entry (b)):** generic pathway and gene-set databases for T3 (the organisers'
    ruling of 2026-10-15; `docs/DATA_POLICY.md`, 'Organiser rulings on record'), queried by code with the gene's name;
    never from a language model's memory or the literature; I name no pathway gene.
  - **Everything else as below:** the screened read of 23 Nov found no organiser answer; T1's third line and its two
    files are in force; six questions and the turn order are with the captain; nothing else runs.
- **State at 2026-11-23, 04:41 UTC (superseded by the bullet above where they differ).** Decision log through
  2026-11-23 (a); the next label is 2026-11-23 (b). The second screened read for the seventh email ran (the captain's
  approval of 13 Nov; reminder): NO ANSWER in the mail thread, the Slack channel or the Data Protocol page; by the
  reading fixed beforehand nothing changes; recorded in `docs/ORGANIZER_QUESTIONS.md` and `docs/DATA_POLICY.md`.
  Nothing runs on the DGX. Tries: T1 463, embryo 88, heart 97, T3 74. No score has been pasted since 20 Nov. Still
  with the captain: the six questions (SO, O_hi, MS, file D, T1 line 1's step back, the heart line's exception), the
  turn order of T1's third line, and the uploads of the handed-over files; all as in the bullets below.
- **State at 2026-11-22, 11:05 UTC (superseded by the bullet above where they differ).** Decision log through
  2026-11-22 (nn); the next label is 2026-11-22 (oo). Local main = origin/main at a6648c1 or one notes commit after
  it. Tries: T1 463, embryo 88, heart 97, T3 74. No score pasted today. Nothing runs on the DGX.
  - **T1'S THIRD SEARCH LINE IS IN FORCE (entry (nn); committed and pushed at 11:04 UTC on 22 Nov, eleven minutes
    after the captain's "yes" of 10:53 UTC, entry (mm)).** The atlas step on file C's command with the correction held
    at 1.0: `T1_s2rcCm20.h5ad` (multiplier 2.0), then `T1_s2rcCm25.h5ad` (2.5). It starts only if A and C are both
    uploaded and C's headline is strictly above A's; turn order line 1, then line 2, then this line, unless the
    captain says otherwise; it is void if C's headline was returned before the entry was committed (the captain was
    asked to hold C until told; the captain is told now). The captain was also told plainly that under this order the
    line may get no slot and that the order is the captain's to change.
  - **ITS TWO FILES ARE BUILT AND COUNT (entries (oo), (pp); state at 12:09 UTC):**
    `runs/manual_upload/T1_s2rcCm20.h5ad` (acce5f5f7097269bffebbfd583833e4c9c1c5ae77f1d5294dde6d13e197f863e) and
    `T1_s2rcCm25.h5ad` (c5fb96a64075b6fbfb2f3066604d622ea8b89ac3a5e4f8a99d0fea7d93011dfd); the builder change is
    merged after an independent review; both sha256 are in entry (nn). They are uploaded only by the line's rule (Cm20
    if C's headline is strictly above A's; Cm25 if Cm20 strictly beats the line's best), last in the turn order unless
    the captain says otherwise. The decision log runs through 2026-11-22 (pp); the next label is (qq). Nothing runs on
    the DGX.
  - **THE SIX OTHER QUESTIONS ARE STILL OPEN** (SO, O_hi, MS, file D, T1 line 1's step back, the heart line's
    exception): the captain's "yes" answered the combined-line question only.
- **State at 2026-11-22, 04:29 UTC (superseded by the bullet above where they differ).** Decision log through
  - **THE CAPTAIN ORDERED A THIRD T1 SEARCH LINE FROM FILE C ("yes", entry (mm)); C'S UPLOAD WAITS UNTIL THE LINE IS
    LOGGED.** Read as the answer to that one question only; the six others stay open. Steps: (1) a score-blind
    designer and a cold challenger write the line (workflow; hand-back to `lift8/`); (2) the main session writes the
    challenger's final text into the log and `docs/S2_RULES.md` (re-pinned) BEFORE C's score and tells the captain
    that C may go up; (3) an engineer and a reviewer make the builder take the line's files under `--captain-ruling
    R22` (row R22's first cell names the options before any build) with a launcher whose identity steps reproduce
    R13's and C's bytes (C: 31477ef7df5d5374a74a8249c570217020c593916ab4d3780049b7b3f87b0620); (4) two builds per file
    in `~/v/scratch/t1line-wt`; (5) record and hand over by the line's rule. The decision log runs through 2026-11-22
    (mm); the next label is (nn).
  2026-11-22 (kk); the next label is 2026-11-22 (ll). Local main = origin/main at a233119 or one notes commit after
  it. Tries: T1 463, embryo 88, heart 97, T3 74. No new score yet today. NOTHING RUNS on the DGX and no agent is out.
  Every read started in this session is recorded.
  - **SIX QUESTIONS ARE WITH THE CAPTAIN (plain text in my messages; act only on the captain's own words; I gave no
    lean on any; a blanket reply does not answer the recused or closed-lever ones: ask one short question per item).**
    (1) Embryo: name SO? (2) Embryo: name O_hi, which also lifts the embryo's depth line for that one file? [the
    analysts' text: `experiments/coords/results/e51_so_ohi/question_for_the_captain.md`; recused]. (3) T1: name MS, so
    that it is built (entry (bb))? (4) Heart: upload the organ-cut file D
    (`.../h29b_heart_organ/question_for_the_captain.md`; recused)? (5) T1 line 1: a step back (1.75, 1.25), to be
    fixed before the 2.0 file's score? (6) The heart search line's named exception (recused)? Defaults: no upload,
    nothing built, the rule stays.
  - **FILES BUILT TODAY AND WHERE THEY STAND (all in `runs/manual_upload/`; sha256 in
    `experiments/builds/results/`).** Handed over, the captain may upload: `T2embryo_s2x_crs.h5ad` (embryo, after C;
    323f2acc...) and `T1_s2rcA15.h5ad` (T1 line 2's next file; e2334137...). On the captain's naming only:
    `T2embryo_s2xs_ov.h5ad` (SO; 77062082d97555961901c811d8a1aaf9b5b2403417c19a6988c816ff553f1165),
    `T2embryo_s2x_ovhi.h5ad` (O_hi; 7e4d1eaafc40080d9d9c461d4b9396f0cee6353686538e79fb92ca31d1d2deeb),
    `T2heart_organ_half.h5ad` (D; 1b4342c1...). On their lines' rules only: `T1_s2rcB25.h5ad`, `T1_s2rcB30.h5ad`
    (after the 2.0 file's headline), `T1_s2rcA20.h5ad` (after A15's), `T3_k23_s2_e4.h5ad` (after the K23 file's;
    9934f1b1...). Queued from before, unchanged: T1 F, C, the 2.0 file; T3 C1, C2, K23, CK on its condition, the
    refused s1; embryo C, `T2embryo_s2`.
  - **ON A YES TO SO OR O_hi:** amend `docs/S2_RULES.md` (the file enters the embryo's order by the captain's naming;
    for O_hi the depth line is lifted for that one file), re-pin, log it, and hand over path, sha256, board and Data
    used (released challenge data only). **ON A YES TO MS:** engineer + reviewer for the builder option, traced dry
    run, two builds where the default command first rebuilds ad8f24c7; place after entry 2026-11-20 (r)'s list. **ON A
    YES TO D:** amend the rules (file D enters the heart's order), re-pin, hand over. **ON THE STEP BACK:** name 1.75
    and 1.25 files in the log before the 2.0 file's score; a builder change and review are needed. **ON THE HEART
    LINE'S EXCEPTION:** engineer + reviewer for the two trend files on C.
  - **WHEN THE CAPTAIN PASTES SCORES:** take the headline per board and stop reading; `results/manual_uploads.csv`,
    one log entry, then each board's pre-set rule mechanically (T1 line 1: 2.5 only if the 2.0 file strictly beats the
    line's best 71.59, then 3.0; T1 line 2: A20 only if A15 strictly beats A's; T3: CK only if C2 and K23 each
    strictly beat 62.30, the E = 4.0 file only if K23 strictly beats B's; embryo and heart: best-headline mark only).
  - **CLOSED TONIGHT, NO FILE:** K26, K25 (earlier), C2A, C1A (T71), T68, rc18 (T69), H30, H31, H32, E53, E54, T66,
    T67. **AHEAD OR PASSED, WITH THE CAPTAIN:** MS (T65), SO and O_hi (E51), and the built files D and the sized
    copy-rows file.
- **State at 2026-11-22, 04:09 UTC (superseded by the bullet above where they differ).** Decision log through
  2026-11-22 (ii); the next label is 2026-11-22 (jj). Local main = origin/main at c4465e9 or one notes commit after
  it. Tries: T1 463, embryo 88, heart 97, T3 74. No new score yet today. Nothing runs on the DGX. The four questions
  with the captain and the two handed-over files are as in the bullets below.
  - **ROW E51 IS READ (entry (ii)): SO passes and O_hi passes; by the read each question goes to the captain and the
    read picks neither.** I am recused on these terms and read the outcome line only. Workflow ww3hcck7v (hand-back to
    `lift8/`): (1) a fresh analyst and a cold checker write the captain's question from the reading, with one marked
    sha256 blank per file: pass the checker's final text on unchanged, the blanks filled once the builds agree; (2) an
    engineer and a reviewer write the build launcher `experiments/builds/embryo_so_ohi.sh` on branch
    `campaign/embryo-so-ohi-build` (no builder code is new). On a ready review: merge, commit, deliver to
    `v/vec-campaign`, dry-run trace, launch in a detached worktree with its own environment and `data/raw` linked,
    PATH exported; on done 0 copy both files to `runs/manual_upload/` (`T2embryo_s2xs_ov.h5ad`,
    `T2embryo_s2x_ovhi.h5ad`), record under `experiments/builds/results/`, fill the blanks. SO needs the captain's
    named order; O_hi needs the captain's named order covering the embryo's depth line. Their places by the order
    rule: SO, then `T2embryo_s2`, then O_hi. Data used: released challenge data only.
- **State at 2026-11-22, 03:42 UTC (superseded by the bullet above where they differ).** Decision log through
  2026-11-22 (hh); the next label is 2026-11-22 (ii). Local main = origin/main at 05e395c or one notes commit after
  it. Tries: T1 463, embryo 88, heart 97, T3 74. No new score yet today. The four questions of the bullet below are
  still with the captain, and the two files handed over today still wait for their uploads.
  - **CLOSED SINCE 03:25.** Row T69: STOP, rc18 ends on T1 (entry (gg)). Row T71 (C1A): its target-free pre-gate
    failed on the skeletal-muscle guard, no try (entry (hh)). With these every candidate of the outside-data round is
    closed.
  - **BUILT SINCE 03:25.** The T3 line's E = 4.0 file: `runs/manual_upload/T3_k23_s2_e4.h5ad`, sha256
    9934f1b104fafed0514cff7bb1dec3f5e362e6a5068ae32711cbfc67f62b90f1 (entry (ff)); uploaded only if the K23 file's
    headline is strictly above B's, after C1, C2, the K23 file, CK on its condition and the refused s1.
  - **STILL RUNNING ON THE DGX.** Only E51 (the embryo overlay reads; worktree `~/v/scratch/e51-wt` at d0f0d6e;
    `~/v/scratch/e51/chain.exit`; at 03:41 UTC its last task was on seed 7 of 10). On its read: copy
    `experiments/coords/results/e51_so_ohi/` back and check sha256; record the read's fixed words only; by the
    embryo's order rule (entry 2026-11-21 (y)) SO (`T2embryo_s2xs_ov`) and O_hi (`T2embryo_s2x_ovhi`) exist only if
    their row's practice read is not void and its pre-set pass rule holds, then each is built twice with one sha256; I
    am recused on their terms, so any figures for the captain are written by fresh analysts.
  - **NOTHING ELSE IS QUEUED WITH AN ENGINEER.** New work starts only from the captain's answers (MS on a yes: builder
    option, review, build; the heart line on a named exception; the step back if named before the 2.0 file's score) or
    from scores the captain pastes (headline only; apply each board's pre-set rule).
- **State at 2026-11-22, 03:22 UTC (superseded by the bullet above where they differ).** Decision log through
  2026-11-22 (ee); the next label is 2026-11-22 (ff). Local main = origin/main at 11b0e26 or one notes commit after
  it. Tries: T1 463 (483 if row T71's pre-gate passes), embryo 88, heart 97, T3 74. No new score yet today.
  - **FOUR QUESTIONS ARE WITH THE CAPTAIN (act only on the captain's own words; no lean from me on any).** (1) T1 line
    1, a step back (1.75, 1.25), to be fixed before the 2.0 file's score; default the rule stays. (2) The heart search
    line's named exception; recused. (3) File D (organ cut), `runs/manual_upload/T2heart_organ_half.h5ad`, sha256
    1b4342c151c19dfc827d4b7a310c6a0fa8908f1e369836ebffb22c27c93b4f00; default no; recused. (4) **Q-MS (entry (bb)):**
    row T65 read MS AHEAD on practice; does the captain name MS? On a yes ONLY: an engineer writes the builder option
    (R13's file with the released per-category momentum added on the direction-group cells before the atlas step, arm
    'all'), a reviewer passes it, a dry run is traced, the file is built twice in a checkout where the default command
    first rebuilds ad8f24c7; its place is after entry 2026-11-20 (r)'s list, before rc18. With no answer nothing is
    written or built.
  - **HANDED OVER TODAY, WAITING FOR THE CAPTAIN'S UPLOADS.** Embryo: `runs/manual_upload/T2embryo_s2x_crs.h5ad`
    (323f2acc32f731819b1a3cc0da32478a51d9e9a10f13dbd940773b113b5644fa), after C, with the analysts' note. T1 line 2's
    next file: `runs/manual_upload/T1_s2rcA15.h5ad`
    (e2334137b54d98b04e9dece90213987e61fe9aa15a3fb8d85e6c262d8fb5ba75), after the queued T1 files and after line 1's
    next file; `T1_s2rcA20.h5ad` (372713686782be273a52427140eab652517e35b679f54fc2a969484dee730218) only if A15
    strictly beats the line's best. T1 line 1's `T1_s2rcB25.h5ad` and `T1_s2rcB30.h5ad` wait on the 2.0 file's
    headline.
  - **CLOSED SINCE 03:00.** T68 NOT AHEAD (the blood line on T1 closes). E54's counted gate failed (embryo 88). H32's
    gate failed earlier. C2A stopped.
  - **RUNNING ON THE DGX (a waiter on each).** rc18's read chain (row T69; `~/v/scratch/t69chain2/chain.exit`;
    worktree `~/v/scratch/t69-wt` at 54e5a7f). Row T71 (C1A; `~/v/scratch/t71/chain.exit`; worktree
    `~/v/scratch/t71-wt` at aa93a69; 3 = pruned by its target-free pre-gate, no try; on a read: copy
    `experiments/t1_m1/results/t71_cmjax/` back, the read's fixed words, totals T1 483). The T3 line's E = 4.0 build,
    third launch (`~/v/scratch/t3outside_9fa509f_s2_e4.0/all.exit`; worktree `~/v/scratch/k26-wt` at 9fa509f; 47 = the
    guard failed and the line ends at 2.0; 50 to 52 = the K23 file did not reproduce): on 0 read `SHA256SUMS` and
    `stamps`, copy `runs/t3_outside_form_9fa509f_s2_e4.0/a/xall-s2.h5ad` to `runs/manual_upload/T3_k23_s2_e4.h5ad`,
    record as `experiments/builds/results/t3line_9fa509f_s2_e4.*`; never open the sidecar. E51 (embryo overlay reads;
    run since 01:28; on its reading the embryo's order rule decides SO and O_hi; I am recused on their terms: fixed
    words only, figures to fresh analysts).
  - **FOR EACH READ THAT COMES BACK AHEAD (T69, T71):** one priced question to the captain in the read's and the row's
    own statements, no lean; nothing is built before the captain's yes. T69's builder draft sits on the old branch
    `campaign/t1-rc18` and must be ported and reviewed afresh.
  - **KNOWN, NOT YET FIXED.** `experiments/t3_decomp/build_combination.sh`, `run_k22.sh` and `run_k23.sh` still name
    the old tables folder; each needs the reviewed fix of entry (dd) before any new launch. Four tests of
    `tests/test_t3_combination_build.py` fail when that file runs after certain launcher test files in one process and
    pass alone.
- **State at 2026-11-22, 02:58 UTC (superseded by the bullet above where they differ).** Decision log through
  2026-11-22 (y); the next label is 2026-11-22 (z), then (aa). Local main = origin/main at de7050b or one notes commit
  after it. Tries: T1 448 (463 once row T69's runs have started: check `~/v/scratch/t69chain2/chain.log` and update
  the totals line), embryo 87, heart 97, T3 74. No new score yet today.
  - **THE THREE QUESTIONS OF THE BULLET BELOW ARE STILL WITH THE CAPTAIN** (T1 step back; the heart search line's
    exception; file D yes or no). Nothing new was asked since.
  - **CLOSED SINCE 02:35.** Row H32's target-free gate said FAIL (entry (w)): no arm, no try. All three outside-data
    heart candidates of this round are closed without a counted read.
  - **BUILT AND WAITING.** T1 line 1: `runs/manual_upload/T1_s2rcB25.h5ad` (sha256
    ccb7c457c0e8ac8ac668bcea17bc61ed22f414f58ab5a63e946cb2f08c4f4e77) and `T1_s2rcB30.h5ad`
    (b639de29bd13d7005450ca53e7cd3dacb2e362a55beca8545937858e6582a024); uploaded only by the line's rule. Embryo:
    `runs/manual_upload/T2embryo_s2x_crs.h5ad` (the copy-rows file at the team's size, sha256
    323f2acc32f731819b1a3cc0da32478a51d9e9a10f13dbd940773b113b5644fa; the captain's order of 21 Nov; its place is
    after C); its price note is being written by fresh analysts (workflow w057t50fr; hand-back to `lift8/`), because I
    am recused: pass the checker's final text on unchanged with path, sha256, board and Data used (released challenge
    data only; no outside data).
  - **RUNNING ON THE DGX (a waiter on each).** T65, T68, E51, E54 as before. rc18's read chain
    (`~/v/scratch/t69chain2/chain.exit`; worktree `~/v/scratch/t69-wt` at 54e5a7f; 30 = a run failed): on 0 copy the
    five task files and the reading from `experiments/t1_m1/results/` there, check sha256, fill row T69's last cell in
    the read's fixed words, totals T1 463. T1 line 2's builds (`~/v/scratch/t1damp/r1/`; worktree
    `~/v/scratch/t1line-wt` at 9a10157): on done 0 copy `X15a/T1_s2rcA15.h5ad` and `X20a/T1_s2rcA20.h5ad` from
    `runs/t1damp_9a10157_r1/` to `runs/manual_upload/`, record as `experiments/builds/results/t1damp_9a10157_r1.*`;
    A15 is the line's next file (after line 1's next file). The T3 line's E = 4.0 builds, second launch
    (`~/v/scratch/t3outside_de7050b_s2_e4.0/all.exit`; worktree `~/v/scratch/k26-wt` at de7050b; 47 = the guard failed
    and the line ends at 2.0): on 0 read `SHA256SUMS` and `stamps` there, copy
    `runs/t3_outside_form_de7050b_s2_e4.0/a/xall-s2.h5ad` to `runs/manual_upload/T3_k23_s2_e4.h5ad`, record; never
    open the build's sidecar. It is uploaded only if the K23 file's headline is strictly above B's, after the queued
    T3 files.
  - **SEARCH LINES IN FORCE.** T1 line 1 (2.5 then 3.0), T1 line 2 (damping 1.5 then 2.0), the T3 line (exponent 4.0;
    entry (x); rules amended and pinned). Heart: only on the captain's named exception. Embryo: none.
  - **C1A (row T71; `campaign/t1-cm-nuclei`, 6e79b6c; reviewed ready).** Its screens are clear (entries (q), (r)).
    Before placement: a second literature screen (a background agent) looks for a published adult statement naming
    Mylpf and Myl1 as fast skeletal-muscle genes; with it all five guard genes have a quoted source and the row is
    placed with X010 amendments naming both statements (the first is Oe, Ojima and Muroya 2021, PMID 33459151);
    without it a fresh analyst decides between the specification's five genes and the labeller's source rule. Row text
    and needs: `lift8/search_lines_build.json`, key c1a. A file from it that reads embryos 13 to 22 carries the run_15
    note into the captain's question.
  - **LAUNCH RULES LEARNED TODAY.** Put `export PATH=$HOME/.local/bin:$PATH UV_FROZEN=1 UV_NO_SYNC=1` in every DGX
    launch command (an older launcher stopped at its tests step with exit 127 without it). A build whose gate wants
    the reading's code byte for byte is launched from a side branch cut at the reading's commit (done for the sized
    embryo file, `campaign/e52-build`). Record script, tests and commit go in one `&&` chain.
- **State at 2026-11-22, 02:34 UTC (superseded by the bullet above where they differ).** Decision log through
  2026-11-22 (t); the next label is 2026-11-22 (u). Local main = origin/main at b193249. Tries: T1 448, embryo 87,
  heart 97, T3 74. No new score yet today.
  - **THREE QUESTIONS ARE WITH THE CAPTAIN (plain text; act only on the captain's own words).** (1) T1's first search
    line, a step back (1.75, and 1.25) if the 2.0 file does not beat 1.5; it must be fixed before the 2.0 file's
    score; default: the rule stays. (2) The heart search line (trend 0.5, then 1.0, on the full trim C): needs the
    captain's own named exception; recused. (3) **File D, the organ cut:** built,
    `runs/manual_upload/T2heart_organ_half.h5ad`, sha256
    1b4342c151c19dfc827d4b7a310c6a0fa8908f1e369836ebffb22c27c93b4f00; the analysts' question is
    `experiments/harness/results/h29b_heart_organ/question_for_the_captain.md` (passed on unchanged; recused); it is
    uploaded only on the captain's yes to it by name; default no. On a yes: amend `docs/S2_RULES.md` (file D enters
    the heart's order), re-pin, hand over path, sha256, board and the Data-used text (released Task 2 heart stages and
    the unlocked E10.5 validation file only; no outside data).
  - **CLOSED SINCE 02:00.** Row K26 stopped at its stage 0 (a diagnostic on outside tables; no try, no file).
    Candidate C2A stopped on its own rule: the heart atlas's zones were assigned by a fitted demultiplexer on a pooled
    library that held in-window cells. **Said to the captain once (entry (q)):** the same step assigned that atlas's
    stage column, which every atlas-based T1 file reads; register row X002 now says so; my reading is that it is the
    residual on record since 19 Nov with its mechanism named; file F (no outside data) does not depend on it. The
    nuclei atlas's filter question is answered (entry (r)): per-experiment doublet models, no experiment mixes a
    permitted embryo with an embryo labelled inside the window; run_15 also held the boundary embryos 23 and 24.
  - **RUNNING ON THE DGX (a waiter on each).** T65 (libraries one after another; then commit the five task files
    UNOPENED, push, and in `~/v/scratch/t65-wt` move the untracked task files aside before `git checkout --detach`,
    then `t65_ms.sh read`). T68 (five runs since 02:15; ten arms counted). E51 (run since 01:28). E54 (fit since
    01:01; its gate counts when it runs). E52 second start (run since 01:54). T1 line builds
    (`~/v/scratch/t1line/r1/`). **rc18** (rows N58 and T69; worktree `~/v/scratch/t69-wt` at 3844836, own environment;
    chain `~/v/scratch/t69chain1/chain.exit`: label-tests, label; label exit 3 = stopped after preflight, not rerun).
    On 0: copy `experiments/external/results/feng_relabel_e18/MANIFEST.json` back, place row T69 (`lift8/row_t69.txt`)
    with the E18 statistics' sha256 in its marked blank and the count T1 449 to 463, commit (ROWC), push to
    `v/vec-campaign-d`, move the untracked record aside in the worktree, check out ROWC, then `rc18.sh tests ROWC`,
    the five runs (`cut chamber`, `cut whole`, `lib 1`, `lib 2`, `lib 3`; they may run in parallel), `read ROWC`; the
    launcher names its exit files with underscores (`label_tests.exit`). **H32** (worktree `~/v/scratch/h32-wt` at
    3844836, own environment with the outside-data group; chain `~/v/scratch/h32chain/chain.exit`: tests, tables,
    gate; 3 = gate FAIL, nothing runs or counts). On 0: commit `runs/h32/gate/gate.json`, `gate.txt` and
    `runs/h32/tables/tables_sha256.json` as `experiments/harness/results/h32_heart_gated_trend/gate/` and
    `.../tables_sha256.json` ON A SIDE BRANCH cut from 3844836 (the launcher wants RUN to differ from ROW in those
    files only), deliver, then `run_h32.sh run ROW RUN` and `read ROW RUN`; sixteen tries, heart 98 to 113, from the
    run step; merge the side branch afterwards.
  - **WORKFLOW RUNNING.** wn2g89uiu: T1 line 2's builder (`campaign/t1-damping-line`), the T3 line (re-based; write
    its text into the log and the rules BEFORE the K23 file's headline, add `--table-exponent 4.0` to row R22's first
    cell then), C1A (row T71, `campaign/t1-cm-nuclei`; its screens are clear; register the adult muscle statement,
    PMID 33459151, as an X010 amendment at placement; any file that reads embryos 13 to 22 carries the run_15 note
    into the captain's question). Do not merge a builder change before a build that must match an earlier reading's
    code.
  - **A slip to avoid (entry (r)).** Where a row sends a read's figures to fresh analysts, read the outcome line only
    (`grep '^OUTCOME'`), never the head of the reading: on H29b the first fourteen lines held a headline figure. I am
    recused on the heart's scope and size; I restated nothing and gave no lean.
- **State at 2026-11-22, 02:01 UTC (superseded by the bullet above where they differ).** Decision log through
  2026-11-22 (o); the next label is 2026-11-22 (p). Local main = origin/main at 95dd312. Tries: T1 438, embryo 87,
  heart 88, T3 74. No new score yet today.
  - **TWO QUESTIONS ARE WITH THE CAPTAIN (plain text; act only on the captain's own words).** (1) T1's first search
    line, a step back (1.75, and 1.25) if the 2.0 file does not beat 1.5? It must be fixed before the 2.0 file's
    score; default: the rule stays (climb only). I give no lean. (2) The heart search line (the trend step at 0.5,
    then 1.0, on the full trim C): needs the captain's own named exception for the two files; the challenger's
    question was passed on unchanged; I am recused. Also said to the captain: by entry 2026-11-20 (r) file F takes the
    first free T1 slot (so F, C, then the 2.0 file), or the captain may keep C then the 2.0 file today.
  - **CLOSED SINCE 01:32: gate N56 VOID (classifier balanced accuracy 0.746 below 0.8) -> rows H30 and H31 end, no try
    (entry (l)).**
  - **RUNNING ON THE DGX (a waiter on each).** T65 (libraries one after another, to about 03:50 UTC; exit
    `~/v/scratch/t65/chain2.exit`; then commit the five task files UNOPENED, push, checkout in the worktree,
    `t65_ms.sh read`). T68 (T66's rerun; worktree `~/v/scratch/t68-wt` at cd3c021; `~/v/scratch/t68/chain.exit`; 3 =
    pruned). E51 (run step since 01:28; `e51/chain.exit`). E54 (fits; `e54/chain.exit`). E52 second start (run step
    since 01:54; `e52b/chain.exit`; if reading (i) holds: `experiments/builds/embryo_copy_rows_sized.sh` twice, hand
    over after C). H29b, the organ cut's read (worktree `~/v/scratch/h29b-wt` at e3cc0e4, own env; `h29b/chain.exit`;
    gate 3 = STOP, 4 = REFUSED on an unmapped label: then fix `configs/t2_chamber_classes.yaml` in a new commit; on
    READ: copy `experiments/harness/results/h29b_heart_organ/` back, commit as B, have FRESH ANALYSTS write the
    captain's question with the figures, then `experiments/builds/heart_organ.sh B <RUN>`). T1 line builds (worktree
    `~/v/scratch/t1line-wt` at e3cc0e4, own env; `bash experiments/builds/t1_fraction_line.sh <commit> r1`; exits in
    `~/v/scratch/t1line/r1/`; on `done.exit` 0 copy `SHA256SUMS` and the stamps into `experiments/builds/results/`,
    copy the 2.5 and 3.0 files to `runs/manual_upload/T1_s2rcB25.h5ad` and `T1_s2rcB30.h5ad`, check sha256, record;
    they are uploaded only on the line's conditions).
  - **WORKFLOWS RUNNING.** wen60aw6i: C2A (row T70, `campaign/t1-cm-segments`) and H32
    (`campaign/heart-trend-sign-gate`) builds, plus four screens (zone provenance; embryo day and somite table;
    skeletal-muscle marker statement; the nuclei atlas's quality-filter question, which bears on every X003 use:
    record its answer in the register whichever way it falls and tell the captain). wn2g89uiu: T1 line 2's builder
    (`campaign/t1-damping-line`, cut from the fraction-line branch; row R22 already names its options), the T3 line
    (re-based onto barred s2 first; write its text into the log and the rules BEFORE the K23 file's headline, add
    `--table-exponent 4.0` to row R22's first cell then), C1A (row T71, `campaign/t1-cm-nuclei`). w2b0b518u: K26
    second fix round (`campaign/moon-t3-map-3`; the pre-gate counted; register rows X013, X014, X015, X017 must each
    name K26 before a run). wkbcxkeij: rc18 without its builder part (`campaign/t1-rc18-2`; rows N58 then T69; X002
    must hold the token `feng_relabel_e18`).
  - **SEARCH LINES IN FORCE (entries (j), (o); rules amended; hash in `docs/S2_RULES.sha256`).** T1 line 1: 2.5 if the
    2.0 file strictly beats the line's best, then 3.0; ends at the first file that does not. T1 line 2: A's headline
    is above R13's, so the next file is `T1_s2rcA15.h5ad` (damping 1.5), then `T1_s2rcA20.h5ad` only if A15 strictly
    beats the line's best (A's); line 1's next file goes first; drop date 25 Nov 12:00 UTC. Embryo: none. Heart: none
    without the captain's named exception. T3: pending its re-base.
- **State at 2026-11-22, 01:32 UTC (superseded by the bullet above where they differ).** Decision log through
  2026-11-22 (k); the next label is 2026-11-22 (l). Local main = origin/main at f8c46dd. Tries: T1 438, embryo 84,
  heart 88, T3 74. No new score yet today.
  - **THE CAPTAIN ANSWERED THE THREE OPEN QUESTIONS (entry (j), ruling row R22, rules amended, hash 93b6335f replaced
    by the one in `docs/S2_RULES.sha256`), verbatim: "yes to all three, go ahead".** (1) A SCORE-GUIDED SEARCH OF
    ROUND STEPS MAY RUN, inside a fixed framework: one open setting per line, round steps on its existing scale, every
    file named in the log before its own score, the next upload named by headline comparisons only (strictly better
    than the line's best, else the line ends), nothing fitted from score sizes, no per-metric value read, closed
    levers outside (sizes, scope, cell numbers; T1 composition, depth, advance draw), two uploads a board a day. THE
    FIRST LINE IS FIXED (before the 2.0 file's score): T1's atlas step, 2.5 if the 2.0 file strictly beats the line's
    best headline, then 3.0 if 2.5 does; it ends at the first file that does not, and at 3.0. Files `T1_s2rcB25.h5ad`
    and `T1_s2rcB30.h5ad` need the builder change (workflow w9wob19fo: branch `campaign/t1-fraction-line`; launcher
    `experiments/builds/t1_fraction_line.sh`), then two identical builds each on `~/v/vec-campaign-d` (nothing runs in
    that checkout itself; T65 runs in its own worktree). Other lines: the same workflow's designers and cold
    challengers per board; write each final rule text into the log BEFORE any of its files is scored; I design none.
    (2) T3'S REFUSED s1 IS RE-ADMITTED: `runs/manual_upload/T3_refused_s1.h5ad`, sha256
    0394575a11fcd4f915e0daf4972e5a783b6b89efa065ed4f85fe2445cbbb1545; first T3 slot free after C1, C2, K23 and CK (on
    its condition); Data used: released challenge data only, no Gata4 validation file, no outside data. (3) T1'S FILE
    F MAY GO: `runs/manual_upload/T1_s2none.h5ad`, sha256
    55eaee547ff18d7477da5b730c34e7de9ab74f85cc12fd01a9eff6450367611c; first free T1 slot by entry 2026-11-20 (r); Data
    used: released challenge data only, no outside data. T1 order from here by the rules on record: F, C, the 2.0
    file, then the line's 2.5 and 3.0 on their conditions; MS, T68 and rc18 only on their own paths.
  - **H30 SECOND START (entry (k)).** The reader is corrected by a judge of record's rule (unlisted matrix rows enter
    nothing), reviewed ready, merged; rows N56 and H30 now name e5ba8d9 on `campaign/moon-heart-chord-4`; ROW =
    f8c46dd06e20a65397190eb8a4c2f41b7d9a0e73; worktree `~/v/scratch/h30b-wt` with its own environment; chain `h30b`
    (tests, gate, tables), exit `~/v/scratch/h30b/chain.exit` (3 = gate FAIL or VOID, the row ends). On 0: copy
    `runs/h30/gate/gate.json`, `gate.txt` and `runs/h30/tables/tables_sha256.json` from the worktree into
    `experiments/harness/results/h30_heart_late_chord/` (gate files under `gate/`), commit those three files and
    NOTHING else as RUN, push to the DGX repository, `git checkout --detach RUN` in the worktree, then `bash
    experiments/harness/run_h30.sh run ROW RUN` and `read ROW RUN` (a second small chain). H31's fix round runs
    (workflow wf1u3dl54, branch `campaign/moon-heart-bridge-2`).
  - **STILL AS IN THE 01:20 BULLET:** T65, E51, E54, E52 (`e52b`) run on the DGX; workflows for T68 (wcust43nx), the
    organ cut (wseung8we), K26 and rc18 (w7y2lvts0), the outside-data design round (wfiigl10k).
- **State at 2026-11-22, 01:20 UTC (superseded by the bullet above where they differ).** Decision log through
  2026-11-22 (h); the next label is 2026-11-22 (i). Local main = origin/main at 2c3d568. Tries: T1 438, embryo 84,
  heart 88, T3 74. No new score. Upload queue unchanged (T1: C, then the 2x file; embryo: the copy-positions file,
  then `T2embryo_s2`; T3: C1, C2, K23, CK on barred s2; heart: nothing queued).
  - **THE CAPTAIN'S ORDER OF 22 NOV (entry (g), ruling row R21, rules amended and re-pinned), verbatim, on a forwarded
    side note about the heart: "yes and go ahead with all the external datasets as much as possible".** (1) The
    heart's organ-cut file D is coded and read: the captain's named exception on scope and size for that one file. The
    main session stays recused and gives no lean; D's specification is the score-blind round's
    (`lift8/job_heart_organ.txt`). The upload still needs: gate passes and read not void, two identical builds, and
    the captain's yes by name on a question with the read's figures WRITTEN BY ANALYSTS OTHER THAN THE MAIN SESSION. D
    drops if not built, read and reviewed by 2026-11-25 12:00 UTC. (2) Outside datasets as widely as possible: a
    renewal of the bulk approval; each source still gets its register row or amendment first; no organiser rule
    lifted; a T2 file with an outside atlas still needs its own priced question under A15.
  - **READ SINCE 00:40 (all closed, no file): T67 NOT AHEAD (T1 timing); H29a both heart trend candidates NOT
    QUALIFIED; K25 FAIL (T3 dispersion; T3 73 and 74 counted); E53 stopped by its target-free gate (clause g1; no
    try); T66 stopped by its reducer before any score (no try).**
  - **RUNNING ON THE DGX (a local waiter on each exit file).** T65 (T1 momentum; worktree `~/v/scratch/t65-wt`, own
    env; the chamber run took 38 min; whole cut and three libraries follow one after another, to about 03:50 UTC; exit
    `~/v/scratch/t65/chain2.exit`; then commit the five task files UNOPENED, push, checkout in the worktree,
    `t65_ms.sh read`). E51 (embryo overlay reads; `~/v/scratch/e51-wt` at d0f0d6e; at its reproduction gate; exit
    `~/v/scratch/e51/chain.exit`). E54 (embryo flow; `~/v/scratch/e54-wt` at a99f203; fits on GPUs 2 and 7; exit
    `~/v/scratch/e54/chain.exit`; its gate E54a counts once it runs). E52 second start (embryo copy-positions at own
    size; `~/v/scratch/e52b-wt` at 2c3d568; exit `~/v/scratch/e52b/chain.exit`, the launcher's own exits stay under
    `~/v/scratch/e52/`; if its reading (i) holds, build with `experiments/builds/embryo_copy_rows_sized.sh` twice and
    hand over at its place after C).
  - **WORKFLOWS RUNNING (hand-backs go to `lift8/`; merge only on a ready review; row and register first, then run).**
    (a) t1-blood-rerun-row (wcust43nx): row T68 = T66's rerun with the gene rule fixed by a fresh analyst (the atlas's
    count matrices carry 25,033 gene ids its gene table of 24,552 does not list); branch `campaign/moon-t1-blood-3`;
    run it in a worktree of `~/v/vec-campaign-d` with its own env (`uv sync --frozen --group external`), data linked
    as for T66; X003 needs a sentence naming T68. (b) heart-organ-cut (wseung8we): branches `campaign/heart-organ`,
    `-2`, maybe `-3`; then gate on E10.5, row H29b, read, builds. (c) heart-chord-reader-fix (wzsj9va23): H30's gate
    crashed in its reader on the same gene-table fact; branch `campaign/moon-heart-chord-4`; then rows N56 and H30 are
    updated to name the corrected code in a new commit, which becomes ROW, and the chain (tests, gate, tables; then
    commit gate.json, gate.txt and tables_sha256.json as RUN; then run, read) starts again in a NEW worktree with its
    own env (`uv sync --frozen --group external`), data linked from `~/v/vec-campaign/data`. (d) outside-data-resume
    (w7y2lvts0): K26 fix round (`campaign/moon-t3-map-2`) and rc18's first build (`campaign/t1-rc18`; rows N58 and
    T69). H31 (heart bridge) waits for (c)'s gene rule, then gets its own fix round. A further score-blind round on
    other outside data is still to launch.
  - **LESSONS OF THIS NIGHT, now in every brief:** (1) no test asserts today's state of the repository's own ledger or
    register (H30's test did; the tests step would have failed after placement); (2) no test pins the bytes of a file
    built on one machine (E52's did; DGX bytes differ); when a pin is needed, measure it on both machines with the old
    code and pin the pair; (3) tie every commit to its test run's result with `&&` (the second batch's commit went in
    beside a failing test); (4) a T1, T3 or heart read from a worktree needs the worktree's own environment (T65); (5)
    a DGX checkout with untracked result files cannot be moved to a commit that tracks them (`~/v/vec-campaign-c`
    holds K24's two reading files untracked): use a worktree.
  - **Open with the captain (plain text, default no each):** a score-guided search of round steps; T3's refused s1
    file; T1's no-outside-data file F. (The organ cut is answered: yes.) Told in the 00:40 report: X003's
    `experimental_batch` column, reported only.
- **State at 2026-11-22, 00:40 UTC (superseded by the bullet above where they differ).** Decision log through
  2026-11-22 (c); the next label is 2026-11-22 (d). Local main = origin/main at 9ff856b. Tries: T1 438, embryo 84,
  heart 88, T3 72. No new score; the upload queue is unchanged (T1: C, then the 2x file; embryo: the copy-positions
  file, then `T2embryo_s2`; T3: C1, C2, K23, CK on barred s2; heart: nothing queued).
  - **READ: T67 (T1 per-gene timing) is NOT AHEAD: the line is closed, nothing is built, no question.** Files in
    `experiments/t1_m1/results/t67_timing/`. `~/v/vec-campaign-d` is idle again at 6c01146 (its worktrees t65-wt and
    t66-wt are in use: do not remove them).
  - **RUNNING ON THE DGX (five chains; a local waiter polls each exit file; never wait in the foreground).** (1)
    **T65** (T1 momentum), second launch: worktree `~/v/scratch/t65-wt` at 34567906…, own `.venv`; chain
    `~/v/scratch/t65_chain2.sh`; exit `~/v/scratch/t65/chain2.exit`. On 0: copy the five task files from the
    worktree's `experiments/t1_m1/results/t65_ms/tasks/` back UNOPENED, commit, push that commit to the DGX repository
    (`git push PrakashDGX_H2:v/vec-campaign-d <commit>:refs/heads/<name>`), `git checkout --detach` it in the
    worktree, run `bash experiments/t1_m1/t65_ms.sh read 34567906a651f08316f5689003bcfa1c4e57e430` there (no
    UV_PROJECT_ENVIRONMENT: the worktree's own environment), read only `read.exit` and `read.out`, copy `t65.json` and
    `t65.txt` back, commit with the outcome in its fixed words. AHEAD = question Q-MS to the captain; else MS ends.
    (2) **T66** (T1 held blood along the outside nuclei atlas): worktree `~/v/scratch/t66-wt` at the row commit R1 =
    d0f0d6e54c5af45b3cd4038fc41354a73358e46c, own `.venv` (with the R-data reader), data linked; exit
    `~/v/scratch/t66/chain.exit` (0 read done: copy `experiments/t1_m1/results/t66_jaxblood/` back and commit; 3 its
    pre-gate pruned it, no try). (3) **E51** (embryo: overlay with the series expression; overlay at the deeper
    depth): worktree `~/v/scratch/e51-wt` at R1; exit `~/v/scratch/e51/chain.exit` (3 = its reproduction gate failed,
    nothing counted); reading in the worktree's `experiments/coords/results/e51_so_ohi/`. (4) **E53** (embryo geodesic
    overlay, GPUs 2 and 7): worktree `~/v/scratch/e53-wt` at R1; exit `~/v/scratch/e53/chain.exit` (4 = its
    target-free gate failed, nothing counted; 3 = a reproduction miss, void); reading in
    `experiments/coords/results/e53_geodesic/`. (5) **H29a** (heart trend step on forms A and B; eight tries counted,
    heart 81 to 88): in `~/v/vec-campaign` itself at R1 (do not move that checkout while it runs); exit
    `~/v/scratch/h29a/chain.exit`; reading in `experiments/harness/results/h29a_heart_trend_forms/`, outcome words in
    `~/v/scratch/h29a/read.out`. For each: copy the result files back, check them against the DGX sha256, fill the
    row's last cell with the read's fixed words, one log entry, update the totals line (E51 85 to 86, E53 87 to 88,
    T66 439 to 448, each only once its run step started). A pass leads to one priced question to the captain (and, for
    the embryo and heart candidates, two identical builds first as their rows say); anything else ends the line.
  - **WORKFLOW RUNNING: open-points-round (task wnrpzb0uz).** Fresh score-blind analysts decide the open design points
    of K25 (T3 dispersion), E54 (embryo flow) and H30 with its gate N56 (heart chord); an engineer and a reviewer
    close E52's one launcher blocker (branch `campaign/embryo-copy-rows-size-3`); a screen confirms the licence of the
    outside nuclei atlas's portal copy (H30's code needs the register sentence 'Licence confirmed for row H30' in
    X003, written only if confirmed). On hand-back: save to the scratch folder `lift8/`; merge any new branch (`-3`),
    rerun the touched tests on main; then ONE commit with rows K25 (text `experiments/t3_decomp/k25_row.txt`, its
    ROW_MUST_STATE phrases word for word; T3 72 to 74), N56 and H30 together (H30's code needs X002 and X003 each to
    hold a sentence starting 'T2-heart amendment for row H30'; the heart use of both atlases rests on the captain's
    bulk approval of outside data and a file would still need the A15 ruling), E54 (its tries after E53's: 89, then 90
    and 91 if its gate passes) and E52 if ready. Places: K25 in `~/v/vec-campaign-c` itself (sync by a git push to a
    branch ref and a checkout; nothing runs there); H30 in its own worktree with its own environment (`uv sync
    --frozen --group external`), because its script does not put its own source first and it moves HEAD from ROW to
    RUN; E54 and E52 in worktrees of `~/v/vec-campaign` (their scripts put their own source first).
  - **RULE LEARNT TODAY (T65's refused first launch).** An environment's package path names the source of the checkout
    it was made in. A script that does not insert its own `src` first (the T1 reads, K25, H29a, H30) therefore runs
    only in a real checkout at the run commit, or in a worktree with its OWN environment (`uv sync --frozen` there;
    0.8 GB; about a minute). The embryo read scripts insert their own source and may share `~/v/vec-campaign`'s
    environment. On the DGX never print `git worktree list` unfiltered (the repositories carry this machine's worktree
    records: thousands of lines).
  - **Open with the captain (plain text, default no each), unchanged:** a score-guided search of round steps; T3's
    refused s1 file; the heart's organ-cut file D; T1's no-outside-data file F. To tell the captain in the next
    report: X003's amendment for T66 names one column beyond the words of the approval of 19 Nov
    (`experimental_batch`, reported only).
- **State at 2026-11-21, 23:47 UTC (superseded by the bullet above where they differ).** Decision log
  2026-11-21 (c) to (bb); the next label is 2026-11-22 (a). Local main = origin/main. Standing by best headline per
  board: T1 71.59, T3 62.30, heart 56.4, embryo 69.9 = 197.04. Only headlines are recorded; the main session gives no
  design lean on any board (entry (t)).
  - **WITH THE CAPTAIN, IN ORDER (all in `runs/manual_upload/`):** T1: `T1_s2rcC.h5ad` (31477ef7df5d), then the 2x
    file `T1_s2rcB2.h5ad` (52897fe375b8a81afffa788592047c30361824148b00a8176514385060518837; row R20; built and handed
    over, entry (bb)). Embryo: `T2embryo_s2x_cr.h5ad` (f20730da8308), then `T2embryo_s2.h5ad` (41982d8afd78). T3:
    `T3_c1_s2.h5ad` (871833086642), `T3_c2_s2.h5ad` (47111ac1d07e), `T3_k23_s2.h5ad` (77ff0b947f7f), `T3_ck_s2.h5ad`
    (e470aaeeab28) only if C2 and K23 each beat 62.30. Heart: nothing queued. When a score comes: headline only into
    `results/manual_uploads.csv`, one log entry, the board's rule, the next file.
  - **OPEN WITH THE CAPTAIN (plain text, each default no):** (1) lift the team's own bars (no number set from a score;
    no per-metric reading) for a score-guided search of round steps, each fixed in the log before its own score: the
    Rules page as screened on 18 Nov permits trying a parameter value, seeing its score and trying another 'however
    many rounds it takes' and prohibits solving for a target property and building a submission from it (entry (aa));
    never lifted: solving hidden properties, the levers closed after the 4 Nov exposure, the main session's recusals.
    (2) T3's refused s1 (`T3_refused_s1.h5ad`, 0394575a11fc; built). (3) The heart's organ-cut file D. (4) T1's
    no-outside-data file F (`T1_s2none.h5ad`, 55eaee547ff1).
  - **RUNNING at 23:58 UTC (2026-11-21); the decision log is at 2026-11-21 (dd), next label 2026-11-22 (a).** (1) ROW
    T67's CHAIN on the DGX: `~/v/vec-campaign-d` at ROWC = 6c011467fba75ef842705ff6f518fee5cf6ba578 (do not sync or
    rsync that checkout while it runs); `~/v/scratch/t67_chain.sh` started 23:56 UTC; log
    `~/v/scratch/t67/chain3.log`, exit `~/v/scratch/t67/chain3.exit` (0 = read done; 4 = the idle pre-gate pruned the
    row, no try; else the failing step); a local waiter polls. On exit 0: copy `experiments/t1_m1/results/t67_timing/`
    back, commit with the outcome in the read's fixed words (t67.txt), fill row T67's last cell and the totals line
    (T1 428), log entry; AHEAD = one priced question to the captain, no file. Row N57 is run (the late labelling is
    usable; record committed). (2) Workflow `moonshots-fix-1` (task wee4h40co): fix and re-review of T66, K25, H30
    (gate row N56), E53, E54 on branches `campaign/moon-*-2`. (3) Workflow `blind-round-2-fix` (task wlgdeo580): fix
    and re-review of E51 (`campaign/embryo-so-ohi-2`; one specification point decided by a fresh analyst), H29a
    (`campaign/heart-trend-forms-2`), E52 the copy-rows file at the team's own size
    (`campaign/embryo-copy-rows-size-2`; the analysts' verdict is build as amended, placed after C). Hand-backs of the
    first engineering rounds: `lift8/blind_eng.json`, `lift8/moonshots_eng1.json`. Deferred: K26, H31.
  - **ROW T65 IS RUNNING TOO (00:08 UTC, 22 Nov; entry 2026-11-22 (a); next log label 2026-11-22 (b)).** Row placed at
    ROWC65 = 34567906a651f08316f5689003bcfa1c4e57e430; tries T1 429 to 438 once its first run starts (T67 took 419 to
    428: its pre-gate passed and its five runs started at 00:02 UTC; the totals line says T1 428). The commit was
    pushed into the DGX repository as branch ref `t65-run` (the running checkout untouched); worktree
    `~/v/scratch/t65-wt` with `data/raw`, `data/derived`, `data/external` linked to the checkout's; chain
    `~/v/scratch/t65_chain.sh` (tests, cut_chamber, cut_whole, lib_1, lib_2, lib_3, sequential), log
    `~/v/scratch/t65/chain.log`, exit `~/v/scratch/t65/chain.exit`; a local waiter polls. **When it exits 0:** copy
    the five task files from `~/v/scratch/t65-wt/experiments/t1_m1/results/t65_ms/tasks/` back UNOPENED, commit them,
    push that commit to the DGX ref, `git checkout --detach` it in the worktree, then `bash
    experiments/t1_m1/t65_ms.sh read 34567906…` (with UV_PROJECT_ENVIRONMENT set to the checkout's `.venv`); read only
    `read.exit` and `read.out`; copy `t65.json` and `t65.txt` back and commit with the outcome in its fixed words;
    fill row T65's last cell and the totals line (T1 438). AHEAD = question Q-MS to the captain, then the builder
    option; anything else ends MS.
  - **Local disk (22 Nov 00:06 UTC):** / 42 GB free, /mnt/md0 82 GB free. The repository's `.claude/worktrees` holds
    191 GB of agents' worktrees (several with their own environment of about 0.8 GB); nothing was deleted. The old
    disk watchdog ended (exit 144, not a low-disk trip); a new report-only one is running (wakes below 8 GB on / or 25
    GB on /mnt/md0). Deleting old worktrees needs the captain's yes and a look at each first (some hold unmerged
    work).
  - **DGX facts (21 Nov 23:43 UTC):** the outside nuclei atlas (`data/external/jax_qiu2024`, 43 embryo files) and the
    heart atlas's clean reduction are on `~/v/vec-campaign-d` and `~/v/vec-campaign`; torch 2.14, torchdiffeq and scvi
    are in the campaign-d environment; the `rdata` package is MISSING in both environments (T66, H30 and H31 need it:
    install with uv from the lock's external group while nothing runs there, and say so in the log).
    `~/v/vec-campaign-d` is at 3d18841 (holds `runs/t1x2_3d18841_r1` and its archive); `~/v/vec-campaign` is at
    59ad7fa; nothing of ours is running on the DGX.
  - **For every hand-back:** read the task's output file, save it under `lift8/`, never the transcript; on a ready
    review merge (keep both sides of the test list in `experiments/external/n54_step.sh`), run the related tests on
    main, amend the register first where a read uses an outside source in a new way, place the row and a log entry in
    one commit, sync the DGX checkout only when nothing runs in it, trace the launcher's refusal path, run, read in
    the row's fixed words; only a pass leads to a build, and several would first need the captain's word (R1 for the
    blood line; A15 for T3's map and the heart's outside data). **Tries:** T1 418, embryo 84, heart 80, T3 72.
    Register: next row X018.
- **State at 2026-11-21, 22:57 UTC (superseded by the bullet above where they differ).** Decision log 2026-11-21 (c) to (aa). Local main = origin/main. Nothing of ours is running on the DGX. **EIGHT SCORES ARE IN FOR
  THE DAY (entries (s) and (w)); only headlines are recorded; the pasted per-metric lines and ranks are recorded
  nowhere and used for nothing; the main session gives no design lean on any board (entry (t)).**
  - **Headlines (first, second upload of the day) and the next file by each pre-set rule:** T1: A 70.8, B 71.59 (the
    board's best; first upload 70.9; R13's file 70.60): both beat 70.60, so next is C
    `runs/manual_upload/T1_s2rcC.h5ad` (31477ef7df5d); the 'more atlases' line may now be weighed (B's sign against
    R13's file is positive). Heart: K1 54.1, K2 49.51; R16's list is spent; the best stays 56.4
    (`T2heart_s2core_half`). Embryo: O 69.9 (the board's best), P 67.11; next C `T2embryo_s2x_cr.h5ad` (f20730da8308).
    T3: exponent 1.5 61.8 (a tie with barred s1), barred s2 62.30 (the board's best): the base B is the BARRED S2;
    next C1 `T3_c1_s2.h5ad` (871833086642), then C2 `T3_c2_s2.h5ad` (47111ac1d07e) whatever C1 returned, then the K23
    file `T3_k23_s2.h5ad` (77ff0b947f7f), then CK `T3_ck_s2.h5ad` (e470aaeeab28) only if C2 and the K23 file each
    return strictly above 62.30. **Standing by best headline per board: 71.59 + 62.30 + (56.4 + 69.9)/2 = 197.04.**
  - **The blind candidate round handed back (scratchpad `lift8/blind_round2.json`; parts per board: design.ready,
    design.fresh, final).** T1: two files survive, both need code and a practice row that reads AHEAD first: MS
    (released per-type momentum on top of the atlas step; 6 to 8 h) and rc18 (the atlas direction to E18; about 1.5
    days; weak); F (`T1_s2none`, 55eaee547ff1, no outside data) is already listed in entry 2026-11-20 (r), gated on
    the captain's answer to QF. T3: one file, the refused s1 with the wild-type-only factor (built, reference-listed,
    `t3_build/refused_s1.h5ad`, 0394575a11fc; no code), to go after the queued list; it needs the captain's word
    because the amendment of 2026-11-19 says it is not uploaded. Heart: H-T1 and H-T1B (the ladder trend on K2's and
    K1's forms; each only if a whole-pool placebo read passes; no code for the files) and D (an organ cut; a day of
    code; a long shot). Embryo: SO (the overlay's positions with the series expression; a practice read first),
    T2embryo_s2 (R4's file, built, 41982d8afd78; needs a record check of rulings R9 and R11), O_hi (the overlay topped
    up to E8.0's depth; the same read). Nobody proposed the copy-rows file at the team's own size; the captain's yes
    of entry (v) covers it, so fresh analysts specify it. The challengers' odds that a list raises its board's best by
    a point: T1 about 0.05, T3 about 0.15, heart about 0.06, embryo about 0.2.
  - **Since the second scores (entries (x), (y)):** (1) THE CAPTAIN ORDERED A 2x T1 FILE, verbatim "I allow and they
    allow too, so go ahead" on a forwarded side note (entry (x); ruling row R20; rules amendment): R13's file with
    `--feng-fraction-mult 2.0 --captain-ruling R20`, to go after C; the main session's reading of the rule on scores
    is in (x) (trying a setting is fine; nothing is fitted from score sizes; about one in ten that the organisers
    would object). (2) The blind round's lists are fixed in entry (y). READY NOW AND HANDED OVER: the embryo's
    `runs/manual_upload/T2embryo_s2.h5ad` (41982d8afd78; R4's file; after C). BUILT BUT WAITING FOR THE CAPTAIN'S
    WORD: T3's `runs/manual_upload/T3_refused_s1.h5ad` (0394575a11fc; the rules note of 19 Nov says it is not
    uploaded). Questions put to the captain in plain text, each with a default of no: re-admit T3's refused s1 after
    the queued list; start the heart's organ-cut file D; upload T1's no-outside-data file F in a free slot.
  - **THE 2x T1 FILE IS BUILDING (entry 2026-11-21 (aa)).** Code reviewed READY and merged (0619ad2); DGX
    `~/v/vec-campaign-d` at 3d18841203032ed84f43cc941f4f73174c6e5e4f; `bash experiments/builds/t1_fraction2.sh <that
    commit> r1` started 23:26 UTC; exit files in `~/v/scratch/t1x2/r1/` (gate, tests, D must give R13's ad8f24c71c90,
    B must give c7240cb6a7f7, X2a, X2b, done); output `runs/t1x2_3d18841_r1/`; a local background waiter is polling.
    On done.exit 0: read stamps and SHA256SUMS, archive on both machines, copy `T1_s2rcB2.h5ad` to
    `runs/manual_upload/`, local format check by exit code, build record
    `experiments/builds/results/t1x2_3d18841.txt`, log entry, hand over (after C in T1's order; data-used text as
    B's). **THE RULES PAGE ON SCORES (screened quote in `docs/ORGANIZER_QUESTIONS.md`, 18 Nov): trying a parameter
    value, seeing how it scores and trying another is permitted however many rounds it takes; prohibited is solving
    for a property of the target and building a submission out of the solution.** Entry (t)'s refusal to tune from
    scores was stricter than that; the bars are the team's own. ONE QUESTION IS WITH THE CAPTAIN (default no): lift
    the team's own bars (no number set from a score; no per-metric reading) for a score-guided search of round steps
    on a few settings per board, each fixed in the log before its own score. Staying closed whatever the answer:
    solving for hidden properties; the levers closed after the 4 Nov exposure; the main session's recusals.
  - **The moonshot round handed back (entry 2026-11-21 (z); scratchpad `lift8/moonshots_round1.json`; per-job
    specifications `lift8/job_moon_*.txt`).** Eight kill tests are with engineers (workflow `moonshots-engineering-1`,
    task w82xiz87h; branches `campaign/moon-*`): T66 T1 blood along the outside nuclei atlas (needs X003 amended
    before the read; a file would amend ruling R1); T67 T1 per-gene timing from the heart atlas's later stages (X002
    amended first); K25 T3 within-label dispersion (released data only); K26 T3 gene-conditioned response map (a file
    needs the captain's ruling on question A15); H30 heart late-anchor chord, its target-free gate first (if the gate
    fails, H30 and H31 end at no cost); H31 heart latent bridge (GPU); E53 embryo geodesic overlay; E54 embryo
    multi-stage flow. Judges' odds: a kill test passes with about 0.05 to 0.2; three board points or more at the test
    about 0.01 to 0.02 per moonshot. Run order when several are ready, by the judges' ranks and board weight: T66,
    E53, H30's gate, T67, K25, then the rest.
  - **Running (four workflows):** `moonshots-round-1` (task w8ary0c7c; hand-back to `lift8/`; a kill test is a counted
    try, row first); `blind-round-2-engineering` (task wany2l7le): three practice reads, each an engineer and a
    reviewer on its own branch (`campaign/embryo-so-ohi`, row E51; `campaign/t1-ms-read`, row T65;
    `campaign/heart-trend-forms`, row H29a), and the copy-rows file at the team's own size (analyst, challenger,
    engineer, reviewer; `campaign/embryo-copy-rows-size`, row E52); `t1-fraction-2x` (task ws65954zt): engineer and
    reviewer on `campaign/t1-fraction-2x`, and a cold challenge of the rule reading. **When each hands back with a
    ready review:** merge (resolve the test-list and harness conflicts by keeping both sides), run the related tests
    on main, place the row and a log entry in one commit, sync a DGX checkout (T1 builds and reads need
    `~/v/vec-campaign-d`, which holds the relabelling's output; never sync a checkout while a run goes in it), trace
    the launcher's refusal path, run, read in the row's fixed words, and only on a pass build twice and hand over.
    **Tries:** T1 418, embryo 84, heart 80, T3 72. **Next labels:** decision log 2026-11-21 (z), then 2026-11-22 (a);
    register row X018.
- **State at 2026-11-21, 22:31 UTC (superseded by the bullet above where they differ).** Decision log 2026-11-21 (c) to (v). Local main = origin/main. Nothing of ours is running on the DGX. **FOUR SCORES CAME IN at
  22:31 UTC (entry (s)), one per board, reported by the captain as the first upload of the day; only the headlines are
  recorded (`results/manual_uploads.csv`); the paste's per-metric lines and ranks are recorded nowhere and are used
  for nothing.**
  - **Headlines and what each board's pre-set rule says next (the second upload of the day):** T1: A (`T1_s2rcA`)
    70.8, above R13's file (70.60), below the first upload (70.9); next `runs/manual_upload/T1_s2rcB.h5ad`
    (c7240cb6a7f7); C (`T1_s2rcC`, 31477ef7df5d) only if B also returns strictly above 70.60. Heart: K1 54.1 (best so
    far 56.4); next K2 `runs/manual_upload/T2heart_s2iso_half.h5ad` (dc436d9225bb) whatever K1 scored; then the
    heart's list is spent. Embryo: O (overlay) 69.9, the board's highest; next P
    `runs/manual_upload/T2embryo_s2x.h5ad` (a03611a343fa), then C `T2embryo_s2x_cr.h5ad` (f20730da8308), each whatever
    the others scored. T3: the exponent-1.5 file 61.8, EQUAL to barred s1 at the reported precision, so not above
    61.8: next the barred s2 `runs/archive/builds_dc7b6ae/t3_build/barred_s2.h5ad` (1b69521e3c28); the exponent-2.0
    file goes instead only if the captain reports a more precise headline for both files with the 1.5 file strictly
    higher. After the second T3 score the base B is fixed (highest headline among barred s1, the 1.5 file, the second
    file; ties to barred s1); then on B: C1 and C2 (B is barred s1 or s2), then the K23 file, then CK only if C2 and
    the K23 file each beat B. All of those files exist for every base in `runs/manual_upload/` and were re-checked on
    21 Nov.
  - **The files of the four scores are ASSUMED** (the first of each pair in the hand-over table); the captain was
    asked to say if any was another file. **Standing by the record's best headline per board:** T1 70.9, T3 61.8,
    heart 56.4, embryo 69.9: 195.85 (193.2 before).
  - **When the next scores come:** add rows to `results/manual_uploads.csv` (headline only), one decision-log entry,
    apply the rules above, name the next files with path, sha256, board and data-used text. Embryo headlines serve the
    final mark only (entry (l)); no per-term value is tabulated, compared or quoted on any board; nothing is
    re-planned from a score.
  - **THE CAPTAIN'S INSTRUCTION AFTER THE SCORES (entry 2026-11-21 (t)), verbatim: "With this feedback, do whatever it
    takes and to whatever extreme, go ahead with optimizing to ensure that our next upload is a lot better than this
    and we have a lot of chances to beat 226 score in the leaderboard".** The main session's reading: cost is lifted;
    the organisers' rule (scores choose among predictions, never compute them) is not; no file is shaped from the
    pasted per-metric lines; a lapse in thought on the embryo is disclosed there; the main session gives no design
    lean on any board. **Running:** workflow `blind-candidates-round-2` (task w4fjqduvg; twelve agents): per board two
    score-blind designers and a cold challenger return a ranked list of further candidates with an order rule. **When
    it hands back:** save the hand-back to the scratchpad `lift8/`; T1 and T3: place the lists in the log before any
    score (a practice row first where a read is asked), build with thin launchers in the manner of
    `embryo_overlay.sh`, review, hand over; embryo and heart: proposals only, nothing joins either list without the
    captain's own named word. **One question is open with the captain (entry (u)), default closed:** does the captain
    reopen, by name, the embryo's and the heart's lists for new files designed without scores (this includes a
    copy-rows file at the team's own size, raised by a side agent's note that the main session may not act on). The
    main session's honest odds for 226: about one in a hundred or less (said to the captain).
  - **THE CAPTAIN REOPENED BOTH LISTS AND ORDERED MOONSHOTS (entry 2026-11-21 (v)), verbatim: "yes, reopen both lists
    and go ahead with wilder moonshot approached to drastically increase our score at all costs"** (answering: reopen
    the embryo and heart lists for new files designed without scores, including a copy-positions file at our own
    size). So the question of the bullet above is ANSWERED: embryo and heart proposals of the blind round may be
    placed in the log and built; a copy-rows file at the team's own size may be specified by fresh analysts and built;
    the main session's recusals, T1's exposed levers and the organisers' rules stand. **Also running:** workflow
    `moonshots-round-1` (twelve agents; per board two score-blind designers and a cold judge; at most two moonshots
    per board, each with a kill test on released data whose row comes first and counts as a try). When it hands back:
    save to `lift8/`; for each surviving moonshot an engineer and a reviewer write the kill test; place its row and a
    log entry in one commit; run on the DGX (GPUs 2 to 7 are ours); read it in its fixed words; only a pass leads to a
    build path.
  - **Tries:** T1 418, embryo 84, heart 80, T3 72. **Next labels:** decision log 2026-11-21 (w) (or 2026-11-22 (a) after midnight UTC); register row X018.
- **State at 2026-11-21, 11:10 UTC (superseded by the bullet above where they differ).** Decision log 2026-11-21 (c) to (r). Local main = origin/main. No score line since 20 Nov 15:30 UTC. Nothing of ours is running on
  the DGX.
  - **THE EMBRYO HAS THREE FILES WITH THE CAPTAIN, IN THE ANALYSTS' ORDER (entry 2026-11-21 (l); the main session is
    recused, executes the rule and states no preference). Each embryo upload takes the first not yet uploaded:** (1)
    O, the overlay file, `runs/manual_upload/T2embryo_s2x_ov.h5ad`, sha256
    1cfc91be3a20bfdd1a4f7bcfebd7316e1143c6d4048e578c1f30bcfbe1e0e1ea; (2) P, the positions file,
    `runs/manual_upload/T2embryo_s2x.h5ad`, sha256 a03611a343faf1d37edf06ce0e2bccc86cce6f14d3580b823864e5faf3ba03f9;
    (3) C, the copy-positions file, `runs/manual_upload/T2embryo_s2x_cr.h5ad`, sha256
    f20730da83084702c50c43d4d4aaa98be8ee801ea443472792172804fea4228e (built twice at 59ad7fa, entry (q), record
    `experiments/builds/results/embcr_59ad7fa.txt`). Board T2:embryo:test_interp; data-used text for each: "Released
    challenge data only. No outside data." No file waits on another's headline. A headline goes into
    `results/manual_uploads.csv` and serves the final mark only (highest headline, earliest on a tie, with the
    captain's go-ahead): no differencing, no per-term value tabulated or compared with another upload's or with 50, no
    new file, re-seed or re-upload from a score. The list is closed at these three; a fourth embryo file needs the
    captain's own named order.
  - **The captain's question on C is ANSWERED (entry 2026-11-21 (r); 13:09 UTC), verbatim: "yes, build it that way"**
    (the file takes the copy's cell count, below the count rule, and E7.5's own size, not the pre-registered size).
    The file was already built that way; nothing changes; C stays third in the embryo's list. Nothing is open with the
    captain. C's practice figures, which decide nothing, were sent with the hand-over: row E50's reading (i) holds on
    all 90 task-seeds; against the positions form the board is -2.70 on practice, behind in every seed; the file aims
    at the baseline's level on the three shape terms, not above it.
  - **When the captain reports an embryo score:** add the row to `results/manual_uploads.csv` (headline only; per-term
    values are not tabulated or compared) and a decision-log entry; nothing is re-planned from it.
  - **The upload timeline (a side agent's note the captain forwarded at about 12:23 UTC, with no words of the
    captain's own; checked against entries 2026-11-20 (k) and (q) and the K23 and CK amendments: it is right).** Task
    3's chain is sequential and needs three upload days at two a day: day 1 exponent 1.5, then (on its score) exponent
    2.0 or the barred s2; day 2, on the base those two scores fix, two of the next files (C1 and C2 if the base is
    barred s1 or s2; C2 and the K23 file if it is an exponent file); day 3 the rest (the K23 file if not yet up; CK
    only if C2 and the K23 file each beat the base). Task 1 and the embryo need two days, the heart one. Starting by
    23 Nov leaves a spare day; 24 Nov is the last start that still reaches every file before 26 Nov 11:59 UTC. Every
    later Task 3 file exists for every possible base in `runs/manual_upload/` (T3_c1_*, T3_c2_*, T3_k23_*, T3_ck_* on
    s1, s1x1.5, s1x2.0, s2) and all fourteen were re-checked today against their recorded sha256 and pass the format
    check, so no step waits on a build: each needs only the pasted score of the step before.
  - **Unchanged queue on the other boards:** T1 `T1_s2rcA` then `T1_s2rcB` (`T1_s2rcC` only if A and B each beat
    70.60); heart K1 `T2heart_s2trimcopy_half` then K2 `T2heart_s2iso_half`; T3 `T3_barred_x1.5`, then
    `T3_barred_x2.0` if 1.5 beats 61.8, else the barred s2; then on base B the C files, the K23 file and CK by entry
    2026-11-20 (q). Every one of these files was re-checked today against its recorded sha256 and passes the format
    check. The heart check returned 'do not pursue' (entry (j)).
  - **DGX:** `~/v/vec-campaign` is at 59ad7fa (holds `runs/embov_3b813b5_r1`, `runs/embcr_59ad7fa_r1` and their
    archives); worktree `~/v/scratch/e50-wt` at 5fcf6b9 holds row E50's run folder `runs/e50`. Launchers on this
    machine may hold no removal of a variable path; the two embryo build launchers take a run name whose folders must
    not exist.
  - **Tries:** T1 418, embryo 84, heart 80, T3 72. **Next labels:** decision log 2026-11-21 (s); register row X018.
- **State at 2026-11-21, 09:20 UTC (superseded by the bullet above where they differ).** Decision log
  2026-11-21 (c) to (o). Local main = origin/main. No score line since 20 Nov 15:30 UTC.
  - **THE EMBRYO'S ORDER RULE IS IN FORCE (entry 2026-11-21 (l), the analysts' rule; the main session is recused,
    executes it and states no preference): each embryo upload takes the first not yet uploaded of O, P, C.** O, the
    overlay file: `runs/manual_upload/T2embryo_s2x_ov.h5ad`, sha256
    1cfc91be3a20bfdd1a4f7bcfebd7316e1143c6d4048e578c1f30bcfbe1e0e1ea (the captain's "Build it"; row R19; built and
    reviewed, entries (j), (k)). P, the positions file: `runs/manual_upload/T2embryo_s2x.h5ad`, sha256
    a03611a343faf1d37edf06ce0e2bccc86cce6f14d3580b823864e5faf3ba03f9. C, the copy-positions file: not built; it exists
    only if row E50's reading (i) holds and the test-board build passes. Board T2:embryo:test_interp; data-used text
    for each: "Released challenge data only. No outside data." No file waits on another's headline; a headline is used
    for the final mark (highest headline, earliest on a tie) and for nothing else: no differencing, no per-term value
    tabulated or compared, no new file from a score. The embryo's list is closed at these three; a fourth needs the
    captain's own named order.
  - **ROW E50 IS RUNNING ON THE DGX (the copy-positions file's practice read; one try, embryo 83 to 84 once its run
    step starts).** Run commit 5fcf6b9ad187e2c5377900426f6dd749c4a1f9d9 (the row's commit); worktree
    `~/v/scratch/e50-wt`; chain `~/v/scratch/e50_chain.sh` started 09:17 UTC (steps tests, gate, run, read; log
    `~/v/scratch/e50/chain.log`, exit `~/v/scratch/e50/chain.exit`; a local background waiter is polling). DO NOT sync
    or check out `~/v/vec-campaign` while it runs (the worktree shares its repository, environment and data). **When
    the chain ends:** exit 0: copy `experiments/coords/results/e50_copy_rows/reading.json` and `reading.txt` from the
    worktree, commit them with the row's outcome in the reading's fixed words, update the totals line (embryo 84). A
    gate exit of 3: stop, report, count nothing. Any other exit after the run step started: the try counts and the
    outcome is void, recorded by hand. If reading (i) HOLDS: sync `~/v/vec-campaign` to the new commit, trace the
    refusal path of `experiments/builds/embryo_copy_rows.sh`, run it (`bash experiments/builds/embryo_copy_rows.sh
    <commit> r1`; exit files in `~/v/scratch/embcr/r1/`; a non-zero CR step means no file and its reason is not read),
    copy the file to `runs/manual_upload/T2embryo_s2x_cr.h5ad`, write the build record and the file's own rules-note
    section. Report E50's price tables (iii) and (iv) to the captain (they decide nothing; the question itself was
    already put, see below). Hand over with path, sha256, board and the data-used text. If reading (i) FAILS or the
    run is void: no file; O and P stand; tell the captain plainly.
  - **The review of the main session's own record (entry 2026-11-21 (n)):** NOT READY on one blocker, the build
    launcher's gate; fixed at ef46975 (the reading must be a tracked file of the build commit, made at the run commit
    5fcf6b9ad187…; full-hash HEAD; builder code unchanged since the run commit; a build counts only with a stamp of
    this commit, clean). The same reviewer verified the fix: READY, no blocker (entry (o); the launcher's
    code-unchanged check was widened on its lesser point). **THE CAPTAIN'S QUESTION WAS PUT AT 09:27 UTC as plain
    text, in the analysts' words, default build** (count below the count rule, E7.5's own size instead of the
    pre-registered size): a 'no' means no copy-positions file; no answer means built by the default; record the answer
    verbatim, or 'no answer by <time>; built by the default', in the file's rules-note section and log entry. The
    earlier plan to hold the question until the price tables existed was the main session's own timing and is
    withdrawn; the tables are reported to the captain when they exist and decide nothing. In the next rules-note
    amendment also correct the sentence that says the question goes 'when the file is ready'.
  - **The heart check: do not pursue** (entry (j)). K1 then K2 stand; nothing is asked of the captain.
  - **Hand-backs in the session scratchpad `lift8/`:** `embryo_copy_rows.json` (first round),
    `embryo_copy_rows_fix.json` (fix round and order rule), `heart_frame_check.json`; scripts `e50_place.py`,
    `e50_chain.sh`.
  - **Launchers on this machine may hold no removal of a variable path** (a built-in safety check refuses the whole
    command): the two embryo build launchers delete nothing and take a run name whose folders must not exist.
  - **Tries:** T1 418, embryo 83 (84 once E50's run step starts), heart 80, T3 72. **Next labels:** decision log
    2026-11-21 (p); register row X018.
- **State at 2026-11-21, 08:42 UTC (superseded by the bullet above where they differ).** Decision log
  2026-11-21 (c) to (k). Local main = origin/main. No score line since 20 Nov 15:30 UTC.
  - **THE EMBRYO'S OVERLAY FILE IS BUILT (the captain's order of 21 Nov, verbatim "Build it", answering a one-question
    dialog on a forwarded side-agent note; entry 2026-11-21 (j), ruling row R19, the rules-note amendment of
    2026-11-21).** `runs/manual_upload/T2embryo_s2x_ov.h5ad`, sha256
    1cfc91be3a20bfdd1a4f7bcfebd7316e1143c6d4048e578c1f30bcfbe1e0e1ea, board T2:embryo:test_interp; data used: released
    challenge data only. Built twice with one sha256 at 3b813b5 on the DGX by `experiments/builds/embryo_overlay.sh`
    (run r1), after the default command gave the positions file's bytes; record
    `experiments/builds/results/embov_3b813b5.txt`; archives `runs/archive/embov_3b813b5/` on both machines. It is row
    E49's overlay form (arm E49c), which did NOT pass its pre-set rule (clause P5); it exists by the captain's order,
    not by the rule. The main session is recused on the embryo's shape forms and states no preference. The independent
    review of the record and the build path said READY, no blocker (entry 2026-11-21 (k), which also corrects three
    wordings; never open the `.build.json` sidecars of test builds). **Before any score of it:** fresh analysts fix
    the order of the embryo's unscored files (the positions file a03611a3, the copy-positions file of entry (h) when
    it exists, the overlay file) and what each score means, in a decision-log entry before any score of the new files.
  - **The heart check: do not pursue** (entry (j); hand-back in the session scratchpad
    `lift8/heart_frame_check.json`). K1 and K2 are fresh draws with zero-mean frame luck; a copy-rows heart file would
    score the baseline on every term but size and needs R16's closed list reopened. K1 then K2 stand; nothing is asked
    of the captain.
  - **The embryo's copy-positions file (the captain's order of entry (h); practice row E50, not yet placed).** State
    at 08:55 UTC. The first round (task w7at8dmdh) handed back: the cold challenger says build as amended (it puts
    about 0.7 on the file landing at the baseline's level on the three shape terms, if the build passes the board's
    limits); the code reviewer says NOT READY, one blocker (the practice arm splits slots by the task's time fraction,
    the row text says one half) and six lesser points. Hand-back saved in the session scratchpad
    `lift8/embryo_copy_rows.json` (parts spec, fin, build, review). The engineer's branch is
    `campaign/embryo-copy-rows` (5dcd7b7). **Running:** workflow `embryo-copy-rows-fix` (task wj4gyfldx): a fresh
    analyst decides the blocker and the open points and returns the final row text; an engineer applies them on
    `campaign/embryo-copy-rows-2` (main merged in); a reviewer checks; in parallel an analyst and a cold challenger
    fix the order of the three unscored embryo files (positions, overlay, copy-positions) and what each headline
    means. **When it hands back with a ready review:** merge the branch (no force); run the full related test list in
    the main checkout (include tests/test_e48_row_order_read.py, tests/test_harness.py, tests/test_hidden_counts.py);
    place row E50 (the analyst's final text, commit filled in; embryo 83 to 84 when its run step starts) with a
    decision-log entry in ONE commit; place the order rule in the log before any score of the new files; sync
    `~/v/vec-campaign` (nothing running there), make a detached worktree `~/v/scratch/e50-wt` at that commit, trace
    the launcher's refusal path, then run `experiments/coords/run_e50.sh` step by step (tests, gate, run, read; exit
    files in `~/v/scratch/e50/`; a gate exit of 3 means stop, report, count nothing); copy `reading.json` and
    `reading.txt` back and commit with the outcome in the row's fixed words. If reading (i) holds: build the test file
    twice with `--copy-rows --copy-rows-test-ok` (the handed-over file's expression is MIXTURE; read by exit code,
    sha256 and the stamp's two fields only; a refusal means no file and the positions file stands), copy under
    `runs/manual_upload/`, write the build record and the rules-note amendment (add there the sentence that the
    overlay dialog also gave 'about +1.2 Overall if it carried over'). **One short question goes to the captain when
    the file is ready, priced with E50's tables:** the file takes the copy's cell count (below the count rule) and
    E7.5's own size (not the pre-registered size); the analysts' default is build; the main session gives no
    preference.
  - **Launchers on this machine may hold no removal of a variable path:** a built-in safety check refuses any command
    whose text has `rm` on a shell variable. `embryo_overlay.sh` deletes nothing and takes a run name whose folders
    must not exist (a rerun takes a new name).
  - **Tries:** T1 418, embryo 83, heart 80, T3 72 (unchanged: R19 is a ruling, not a try). **Next labels:** decision
    log 2026-11-21 (l); register row X018.
- **State at 2026-11-21, 00:54 UTC (superseded by the bullet above where they differ).** Decision log
  2026-11-21 (c) to (i). Local main = origin/main (9100471 or later). No score line since 20 Nov 15:30 UTC.
  - **The captain, 21 Nov about 08:10 UTC, verbatim: "Proceed with decisions that help us place at the top of the
    leaderboard at all wild costs"** (entry 2026-11-21 (i), with the AI assistant's reading: score decisions without a
    question first; no organiser rule, recusal, exposed lever or captain-fixed list lifted). Running since 08:18 UTC
    besides the embryo round: workflow `heart-frame-check` (task wfxrciozn (relaunched after a schema slip of the main session); read-only: do the heart candidates face
    the same frame luck, are they already the copy's rows, does the construction apply). Its verdict may carry one
    short priced question for the captain (a new heart file needs a named exception to R16). Not built, available at
    one word: the K23 table at exponent 4 as an extra T3 candidate.
  - **THE CAPTAIN'S ORDER OF 21 NOV, ABOUT 08:00 UTC (entry 2026-11-21 (h)): the embryo's shelved second construction is
    reopened** (verbatim: "reopen it [my ordeer]", on a forwarded side-agent note). The construction: the file takes
    the baseline copy's own rows and count, so its positions equal the copy's row for row. The main session is recused
    on this term: fresh analysts decide, it executes and states no preference. Workflow `embryo-second-construction`
    (task w7at8dmdh; engineer branch `campaign/embryo-copy-rows`): specify, challenge, build, review. **When it hands
    back:** verdict `cannot_be_built` or `do_not_build` goes to the captain in plain words; on `build_as_amended` with
    a ready review: merge; if a practice read is specified, place its row and a log entry in one commit, sync a DGX
    checkout, run its launcher detached and read the result in its fixed words; then build the test-board file twice
    with the option's flag (read by exit code, sha256 and the stamp's commit and clean fields only), copy it under
    `runs/manual_upload/`, write the build record, the rules-note amendment and the data-used text the analysts
    drafted, and hand it over with path, sha256 and board. The order of the two unscored embryo files and what each
    score means is the analysts' rule, fixed in the log before any score. If the analysts name one short question for
    the captain, put only that question.
  - **Outside data on the two position boards: scoped, not pursued** (entry 2026-11-21 (g); hand-back in the session
    scratchpad `lift8/t2_outside_scoping.json`). Fresh analysts and challengers returned 'do not pursue' for the heart
    and for the embryo; no row is opened. With K24 this closes the outside-data push on every board except the two
    parked lines.
  - **K24 is read: NO FILE** (entry (f); row K24 filled; six tries, T3 72). The reading's statement: no wider table
    beats the named form at a like size; the K23 file stands alone. The three-study table lost to the two-study table
    at the same size in every cell; the bulk table (14 lines of the many-line set) did not beat its shuffled control
    in any cell; the joint was dropped by the row's window (correlation -0.03). By the row, these tables, their joint
    and other sizes are never read again on this fold. Register rows X015 and X017 say they enter no file. **The T3
    outside-table line is finished.** Nothing is asked of the captain.
  - **With the captain, unchanged** (uploads are the captain's acts; up to two per board a day): T1 `T1_s2rcA` then
    `T1_s2rcB` (C only if A and B each beat 70.60); heart K1 `T2heart_s2trimcopy_half` then K2 `T2heart_s2iso_half`;
    embryo `T2embryo_s2x`; T3 `T3_barred_x1.5`, then `T3_barred_x2.0` if 1.5 beats 61.8, else the barred s2; then on
    the base B that entry 2026-11-20 (q) fixes: C1 (s1 or s2 bases only), C2, the K23 file, and CK only if C2 and the
    K23 file each beat B. Hashes: `experiments/builds/results/t3two_fb01487.txt`, `t3k23_49f3fd5.txt`,
    `t3ck_20a0e1a.txt` and the earlier build records. When a score line comes: add it to `results/manual_uploads.csv`
    and the decision log and apply that board's pre-set rule.
  - **Nothing is running** on the DGX for this project. Checkout `~/v/vec-campaign-c` is at 88f381d (K24's commit;
    `runs/k24_88f381d` holds the task files); outside data under `~/v/vec-ext/t3` (tables = the final three-study and
    bulk tables; `tables_two_study`; `tables_first_bulk_run_keep` is the never-read first bulk run).
  - **Parked, one word reopens:** the GEO single-study bulk forms (closed in code); the T1 extra atlases in their
    strict form.
  - **Tries:** T1 418, embryo 83, heart 80, T3 72. **Next labels:** decision log 2026-11-21 (k); register row X018.

- **State at 2026-11-21, 00:42 UTC (superseded above where they differ).** Decision log
  2026-11-21 (c), (d), (e). Local main = origin/main (88f381d or later). Uploads and the pre-built T3 files are as in
  the bullets below; no score line since 20 Nov 15:30 UTC.
  - **Row K24 is PLACED and RUNNING** (row commit 88f381ddbcf29defcebd9bea7906a9122f7797b2; ledger row K24 after K23;
    entry 2026-11-21 (c); tables pinned in `experiments/t3_ext/TABLES_K24.sha256`: bulk 6833d8d7..., sc3 14aa9ad6...,
    two 9490f524...). Launched 00:41 UTC on the DGX checkout `~/v/vec-campaign-c` with `run_k24.sh <commit>`; waiter
    bk7qo7kud polls `~/v/scratch/k24_88f381d/all.exit`. **When it ends:** exit 0 is a reading (VOID included); 22 the
    pre-gate refused (nothing counted); 25 no try left; 26 the guard part failed (a VOID reading counts the failed
    guard arms only); 30 or 43 incomplete (rerun the launcher at the same commit); 41 a crash in the read. Copy back
    `experiments/t3_decomp/results/k24_wider_tables/` (reading.json, reading.txt) and `runs/k24_88f381d/pregate.json`
    and `guard.json` (into the results folder); commit; fill the row's last cell with the reading's fixed statements
    only and the count (T3 66 plus the tries that ran: 9, 6 or 3); update the totals line; a decision-log entry (label
    2026-11-21 (f)). **If the reading names an arm:** `bash experiments/t3_decomp/build_k24_form.sh <HEAD> <B>
    sc3|bulk|joint` for every base (s1, s1x1.5, s1x2.0, s2) on `~/v/vec-campaign-c`, each built twice with identical
    bytes; a rules-note amendment (`docs/S2_RULES.md` with its sha file); register rows already exist (X013 to X015,
    X017); a data-used text; one short question to the captain on the file's place with the default 'last, after CK'.
    **If none:** the K23 file stands alone; nothing changes in T3's list.
  - **What is NOT to be done:** no other size, weight or combination of these tables on this fold; no return to the
    first run's bulk table; no guard change for this row.
  - **Tries:** T1 418, embryo 83, heart 80, T3 66 (K24 adds 9, 6 or 3 when it runs its tries). **Next labels:**
    decision log 2026-11-21 (f); register row X018.

- **State at 2026-11-21, 00:06 UTC (superseded above where they differ).** Decision log
  2026-11-20 (hhh) to (jjj) and 2026-11-21 (d), (e). Local main = origin/main (9954da3 or later). Uploads and the
  pre-built T3 files are as in the bullets below; no score line since 20 Nov 15:30 UTC.
  - **E49 is read: no candidate file** (entry (jjj); row E49 filled; seven tries, embryo 83). Nothing is asked of the
    captain; the main session is recused on the embryo's shape forms and states no preference.
  - **The bulk guard G5 is changed for the bulk route only** (entry 2026-11-21 (d); merged at 8f59260): margin above
    zero. **The many-line set is rebuilt under it** (entry (e)): 14 lines, 63 knockout and 64 wild-type libraries;
    bulk table sha256 6833d8d7...; the first run's table (865a5530...) is never read. DGX checkout
    `~/v/vec-campaign-c` is at 103a1be.
  - **Row K24 is NOT placed.** Code: branch `campaign/t3-wider-tables-3` (36a0c53) has two review blockers closed;
    workflow `t3-k24-amend` (task wgd9buxh5; branch `campaign/t3-wider-tables-4`) implements the two clauses a cold
    challenger restated (entry (iii); specification in the session scratchpad `lift8/k24_cold_amendments.txt`) and a
    manifest field in `experiments/t3_ext/combine.py` (`lines_held_by_study`). STILL OPEN after it: the last delta
    review's blocker (in `pregate_k24`, check the exponent bound (s1) BEFORE forming T ** e_T, so that a nearly
    constant table drops instead of raising; the reviewer's replacement code is in its hand-back; reviewer agent
    af9b13d17501dbca2 can be resumed by SendMessage to check the fix). Then: merge; sync; run `all "E-MTAB-7322"` once
    more so the manifest gains the field (the bulk csv must stay 6833d8d7...); on the DGX run `k24_wider_tables.py
    listing` and `pins` with `VEC_T3_EXT_TABLES=~/v/vec-ext/t3/tables` and
    `VEC_T3_EXT_TABLES_TWO=~/v/vec-ext/t3/tables_two_study`; commit the listing and pins; place row K24 (draft
    `lift8/k24_row_final_draft.txt`, to be completed with the restated clauses, the try counts of the amendments file,
    both runs' counts from entry (e), the sentence that the first run's table is never read and that a pre-gate miss
    reads nothing, and the engineer's row-text changes) with decision-log entry 2026-11-21 (c) in ONE commit; sync;
    `run_k24.sh <commit>` detached; wait on `~/v/scratch/k24_<7>/all.exit`. The K24 builder cites entry 2026-11-21
    (c): keep that label for the placement.
  - **GEO bulk forms: parked** (entry (iii)); branch `campaign/t3-ext-geo-bulk` not merged.
  - **Tries:** T1 418, embryo 83, heart 80, T3 66. **Next labels:** decision log 2026-11-21 (c) (K24 placement,
    reserved), then (f); register row X018.

- **State at 2026-11-20, 23:32 UTC (superseded above where they differ).** Decision log
  2026-11-20 (eee) to (hhh). Local main = origin/main (1c75cd5 or later). Everything about uploads, the pre-built T3
  files, T64 and E49 is as in the next bullet.
  - **The many-line set E-MTAB-7322 is pinned, mapped and run once.** Pins (25 lines, 103 knockout and 114 wild-type
    embryos, two readers, register row X017) at 270fff8; the gene-id map `experiments/t3_ext/panel_ids.tsv` with
    `PANEL_IDS_SHA256` (reviewed: ready) at d51640a. First run on `~/v/vec-campaign-c` at 71bef8c: 12 of 25 lines pass
    every guard (57 knockout, 37 wild-type libraries); verdict copies in `~/v/scratch/dmdd_pins/first_run/` on the
    DGX. Combine refused because the other datasets' stored plans predate the dataset-file change, so the whole chain
    was started again for `GSE168768 GSE186934 GSE204906 GSE277839 GSE299068 E-MTAB-7322` (waiter btwyreynw; exit file
    `~/v/vec-ext/t3/logs/ALL.exit`). A copy of the tables folder before the bulk run is
    `~/v/vec-ext/t3/tables_before_bulk_keep`.
  - **The whole-embryo guard (G5) is with analysts** (entry (hhh); workflow `t3-ext-whole-embryo-guard`, task
    w87x9vftw; branch `campaign/t3-ext-g5` if they decide a change): at its design number it drops 46 of 111 genuine
    whole-embryo control libraries, cutting one tight cluster in two. They decide from guard-level figures only. If a
    change is decided, implemented and reviewed ready: merge, sync, run `all "E-MTAB-7322"` again (plans stay valid:
    the dataset file does not change), read verdicts only, write a decision-log entry with the new unit counts. If
    they say keep: the present tables stand. **No bulk table value has been read and none may be before row K24 is
    placed.**
  - **Row K24** (workflow `t3-wider-tables-row` finished: planner, challenger verdict run as amended, engineer branch
    `campaign/t3-wider-tables` 7304c05, reviewer NOT READY with one blocker). The package is in the session scratchpad
    `lift8/k24_package.json`; the row text draft in `lift8/k24_row.txt` (it holds raw pipe characters inside cells:
    write `abs` instead before placing). A fix round runs (workflow `t3-k24-fix`, task wm5m03cbo; branch
    `campaign/t3-wider-tables-2`). After a ready check: merge; freeze the tables (`~/v/vec-ext/t3/tables` for sc3 and
    bulk, `tables_two_study` for the named form's table); on the DGX run `k24_wider_tables.py listing` and `pins` with
    `VEC_T3_EXT_TABLES` and `VEC_T3_EXT_TABLES_TWO` set, commit the listing `experiments/t3_ext/TABLES_K24.sha256` and
    the pins; place row K24 and decision-log entry 2026-11-21 (c) in ONE commit that becomes the row's commit; sync;
    `run_k24.sh <commit>` detached; wait on `~/v/scratch/k24_<7>/all.exit`. Preconditions checked: `runs/k23_b706e82`
    is in the checkout; `tables_two_study/COMBINE.json` is at 719c8819.... The challenger's own odds that the row
    names a file are about one in ten: an arm must beat K23's named form head to head.
  - **GEO bulk forms** (workflow `t3-ext-geo-bulk-open`, task w7w5t1rq8; branch `campaign/t3-ext-geo-bulk`): closing
    the two blockers; merge its first commit after a ready verification and the flag commit only on `ready_to_open`.
    Then the Zmym2, Uhrf1 and Rreb1 studies go to the plan step (two readers each, register rows from X018). Any bulk
    table that changes after K24's listing is pinned needs its own row.
  - **Tries:** T1 418, embryo 76 (E49 running: 8 of 9 tasks done at 23:16 UTC), heart 80, T3 66. **Next labels:**
    decision log 2026-11-20 (iii) or 2026-11-21 (c); register row X018.

- **State at 2026-11-20, 22:51 UTC (superseded above where they differ).** Decision log
  2026-11-20 (yy) to (eee) and 2026-11-21 (a), (b). Local main = origin/main (936b196 or later).
  - **With the captain, no score reported yet** (uploads are the captain's acts): unchanged from the bullet below (T1
    A then B, C on its condition; heart K1 then K2; embryo `T2embryo_s2x.h5ad`; T3 `T3_barred_x1.5.h5ad`, then x2.0 or
    the barred s2 by entry (k)'s rule). When a score line comes: add it to `results/manual_uploads.csv` and the
    decision log and apply that board's pre-set rule.
  - **T3 files built ahead for every possible base B, now four per base:** C1 (s1 and s2 only), C2, the K23 file, and
    the combination file CK (C2's factor times the square of the two-study table; row R18, decision log 2026-11-21
    (b); uploaded only if C2 and the K23 file each beat B). Records `experiments/builds/results/t3two_fb01487.txt`,
    `t3k23_49f3fd5.txt`, `t3ck_20a0e1a.txt`; copies under `runs/manual_upload/`.
  - **T64 is read: candidate D is NOT AHEAD** (entry (zz)); the line ends; no question to the captain. Tries: T1 418.
  - **T3 outside knockouts.** Cmtr2 (GSE270969, X016) was fetched and then withdrawn before any table: its wild-type
    libraries are too shallow for the fixed cell call (entry (ddd)); no guard was changed. The macOS companion-file
    rule is merged. The Kmt2d pair makes a three-study single-cell table (`~/v/vec-ext/t3/tables`,
    `tables_three_study`), read by no row yet; `tables_two_study` is what K22, K23 and every built T3 file read.
  - **The bulk whole-embryo intake is merged** (e0a58a1; verified ready for the many-line form; the two GEO bulk forms
    stay closed in code, `t3x_common.GEO_BULK_OPEN = False`). **The many-line set E-MTAB-7322 is in the dataset list**
    (936b196; 29 withheld lines; no pins yet). Done on the DGX checkout `~/v/vec-campaign-c` (at 936b196): plan,
    header step, plan again: the line set is FINAL: 25 lines, 103 knockout and 114 wild-type embryos (plan table
    sha256 b13278ab...; header sha256 483aac11...). The hand check the verifier asked for passed on the real sample
    table: all 2,274 withheld rows are withheld by their line column, none by name or allele only, none linked. First
    reading: session scratchpad `lift8/dmdd_first_reading.txt`; the second reader (an independent agent) is reading
    the local metadata copies under `~/.cache/vec-ext/t3/meta/E-MTAB-7322/` on this machine. **Next, in order:**
    compare the two readings (any doubt on a line drops that line; a doubt on the set stops it); register row X017 in
    `EXTERNAL_DATA.md` BEFORE any count is fetched; pins with `fetch.py pin-template E-MTAB-7322 --dest D --series`
    and `--source-file` per arm with exactly the final line set, then `fetch.py datasets-sha`, tests, commit naming
    both readers; plan again with pins; the gene-id map (`make_panel_ids.py --dest ~/v/vec-ext/t3 --accession
    GSE168768 --gsm GSM5169142 --file GSM5169142_Tbx1-wt-E85_features.tsv.gz --refstats
    ~/v/vec-ext/t3/reference/refstats.npz --extra-file <the pinned lines>`; commit `panel_ids.tsv` and
    `PANEL_IDS_SHA256` in one reviewed commit); `run_dgx.sh "$PWD" all "E-MTAB-7322"` (the count matrix is about 380
    MB); read QC verdicts only. Expect units refused by the somite guard (the two arms' somite ranges do not overlap
    for seven lines) and by the sex guard (the map has no chromosome column).
  - **Row K24** (one practice read of the wider tables: the three-study single-cell table and the bulk table, each
    brought to the size of the form K23 named, with shuffled controls and size twins) is being planned, challenged,
    built and reviewed by workflow `t3-wider-tables-row` (task w3ofznofn; branch `campaign/t3-k24` when built). The
    planner did not know that several lines' knockouts are much younger than their controls (target-free somite facts
    in the first reading): give those facts to a fresh challenger before the row is placed. The row, its log entry and
    the new tables listing and pins go in before any run; K22's and K23's pins stay as they are.
  - **Running on the DGX:** E49 (embryo shape interpolation; `~/v/scratch/e49`, worktree `~/v/scratch/e49-wt`; 7 of 9
    run tasks were done at 22:30 UTC; when `chain.exit` is 0 copy `experiments/coords/results/e49_shape_interp/` back,
    commit with the outcome quoting only the fixed reading, apply the row's outcome-to-file rule; the main session is
    recused and states no preference).
  - **Tries:** T1 418, embryo 76 (E49 running), heart 80, T3 66. **Next labels:** decision log 2026-11-21 (c) (or
    2026-11-20 (fff) for an event of the 20th); register row X017.

- **State at 2026-11-20, 20:33 UTC (superseded above where they differ).** Decision log
  2026-11-20 (nn) to (xx). Local main = origin/main (4fb01c5 or later).
  - **With the captain, no score reported yet** (uploads are the captain's acts): T1 A then B (C on its condition); heart K1
    then K2; embryo `T2embryo_s2x.h5ad`; T3 `T3_barred_x1.5.h5ad`, then x2.0 or the barred s2 by entry (k)'s rule.
    When a score line comes: add it to `results/manual_uploads.csv` and the decision log and apply that board's
    pre-set rule.
  - **T3 files built ahead for every possible base B** (entry (q) fixes B from the returned headlines; only B's files
    go up; records `experiments/builds/results/t3two_fb01487.txt` and `t3k23_49f3fd5.txt`; copies under
    `runs/manual_upload/`). Order on B: C1 (`T3_c1_<B>`, only when B is s1 or s2), C2 (`T3_c2_<B>`), then row K23's
    file (`T3_k23_<B>`: the base's factor times the square of the two-study outside-knockout table; rules-note
    amendment committed at 3eac9ff; data used: released data plus GEO GSE204906 and GSE277839). Hashes are in the two
    records.
  - **K23 is read** (entries (ss), (uu)): valid, six tries (T3 66); it names the exponent-2 file; exponent 4 was ahead
    in three cells and behind in one, so it is not named, and no further exponent of this table is read on the fold. A
    combination candidate (C2 times the square of the table; upload only if C2 and the K23 file each beat B) is with
    agents: workflow `t3-combination-candidate` (task wn9uem8gt, branch `campaign/t3-combination`); its row, log entry
    and rules-note amendment must be placed before any score of C2 or of the K23 file exists.
  - **T3 outside knockouts, more data.** Register rows X011 to X015. The DGX tables folder `~/v/vec-ext/t3/tables`
    holds the TWO-study tables (what K22, K23 and the built files read); the three-study tables (with the Kmt2d pair,
    entry (qq)) are at `~/v/vec-ext/t3/tables_three_study` and are read by no row. Workflow `t3-ext-fix-and-bulk`
    (task wihmjv0k1; branches `campaign/t3-ext-4`, `campaign/t3-ext-bulk`) closes three review blockers and builds the
    reviewed bulk whole-embryo intake (entry (vv)), including the 73-line Mouse Genetics Project set (ArrayExpress
    E-MTAB-7322) through its public accession-keyed count matrix; metadata copies are in the session scratchpad
    `dmdd_meta/`. After a ready review: merge; write the gene-id map from a pinned sample's features file; add dataset
    entries (withheld lines: Trps1; B9d2, Cc2d2a, Ift140, Kif3b, Rpgrip1l, Anks6, Atxn10; and on the first reader's
    doubt Kifap3, Mad2l2, Tmem98); plan; two readers; register rows; pins; fetch, process, combine; a new listing;
    then ONE new practice row for the wider tables, placed before its run.
  - **Running on the DGX:** E49 (embryo shape interpolation; row E49, run commit d76ccd5; `~/v/scratch/e49`, worktree
    `~/v/scratch/e49-wt`; the run step started 19:30 UTC; when `chain.exit` is 0 copy
    `experiments/coords/results/e49_shape_interp/` back, commit with the outcome quoting only the fixed reading, apply
    the row's outcome-to-file rule; the main session is recused and states no preference). T64 (T1 candidate D's
    practice row; run commit 4fb01c5 on `~/v/vec-campaign-d`; chain `~/v/scratch/t64_chain.sh`, exits in
    `~/v/scratch/t64/`): after the five run steps exit 0, copy the five task files under
    `experiments/t1_m1/results/t64_heldmom/tasks/` back UNOPENED, commit them, sync the checkout, run `bash
    experiments/t1_m1/t64_heldmom.sh read 4fb01c5...` detached, copy `t64.json` and `t64.txt` back, commit with the
    outcome in the read's fixed words. AHEAD allows one question to the captain (reverse R1 for that file; default
    keep R1; no advice); anything else ends candidate D.
  - **T1 outside atlases: parked** (entries (rr), (tt)): the second cold check put a 35% chance of a mistake on the
    relaxed stage rules; under the strict rule the new atlases give little; the expected gain does not cover the risk
    of a barred sample. Branch `campaign/t1-multi-atlas-3` (6bdf01a) is kept with one known blocker. One word from the
    captain reopens it in the strict form. The atlas in use (X002) states no plug-day convention in its record or its
    papers' main text; its half-day labels are used as labelled (register row updated).
  - **Tries:** T1 413 (T64 adds five), embryo 76 (E49 running), heart 80, T3 66. **Scratch hand-backs** (session
    scratchpad `lift8/`): k23_package, t64_package, t3x_more, t3x_bulk_design, t1_atlas_fix3, kmt2d_first_reading.

- **State at 2026-11-20, 19:23 UTC (superseded above where they differ).** Decision log
  2026-11-20 (w) to (mm). Local main = origin/main (76a5601 or later).
  - **With the captain, no score reported yet** (uploads are his acts): the files of the next bullet, unchanged: T1 A
    then B (C on its condition), heart K1 then K2, embryo `T2embryo_s2x.h5ad`, T3 `T3_barred_x1.5.h5ad` then x2.0 or
    the barred s2 by rule. When a score line comes: add it to `results/manual_uploads.csv` and the decision log, apply
    that board's pre-set rule, and for T3 build C1 and C2 only after BOTH of entry (k)'s headlines
    (`build_two_knockout.sh <commit> <base>` on `~/v/vec-campaign-c`).
  - **No embryo row-order file exists or can exist** (entries (ff), (hh)): E48 ran (valid; the construction is not
    voided; no effect on practice), and the build for the test board stops because the option needs a count the board
    does not give. Nothing replaces it.
  - **Running on the DGX: E49, the embryo's shape-interpolation read** (row E49, run commit d76ccd5; chain script
    `~/v/scratch/e49_chain.sh`, exit and times in `~/v/scratch/e49/`, detached worktree `~/v/scratch/e49-wt`; steps
    tests, describe, gate, run, read; embryo count 76 to 83, or 81 if the target-free step drops F1). When
    `chain.exit` is 0: copy `experiments/coords/results/e49_shape_interp/` from the worktree, commit it with the
    outcome in E49's row quoting only the fixed reading, then apply the row's outcome-to-file rule (at most two files;
    a test-board build needs `--shape-form X --shape-form-test-ok` after the default reproduces a03611a3...). The main
    session is recused on the occupancy term: the analysts' rule decides, the main session states no preference.
  - **T3 outside knockouts.** The pipeline ran end to end on the DGX (`~/v/vec-campaign-c`, data under
    `~/v/vec-ext/t3`): register rows X011 to X014; the two-study table (Dnmt3a, Dnmt3b, Usb1) is listed in
    `experiments/t3_ext/TABLES.sha256` and pinned in `experiments/t3_decomp/outside_forms.py`. Practice read K22
    (entry (ll)): valid, three tries (T3 60), outcome "no file: nothing qualifies"; it missed the bar by 0.02 in one
    cell. Two workflows follow (entry (mm)): `t3-outside-table-exponents` (task wb840flet: plan, challenge, build and
    review row K23, the same table at larger exponents with shuffled controls; branch `campaign/t3-k23`) and
    `t3-outside-knockouts-more-sources` (task wpuefr322: archive blockers and an unfiltered-matrix intake on branch
    `campaign/t3-ext-3`, plus a hunt and a reviewed design for bulk whole-embryo knockout data). After a ready review:
    merge, add the datasets (two readers read every pinned record; register row first), rerun `run_dgx.sh "$PWD" all
    "<every pinned accession>"`, write the new listing and pins, place the next row before its run. Held on the first
    reader's doubt: GSE246181 and GSE317197. GSE299068 (Kmt2d) is a low-odds candidate left for the next batch.
  - **T1 outside atlases.** The second fix round was verified NOT READY (the stage parser again). Third round in
    flight (task wmbo8t1ua, branch `campaign/t1-multi-atlas-3`): pinned samples by the hash of the full verbatim
    record as the first line, the parser as a second line. When ready: merge, sync, plan run (metadata only,
    destination `~/v/vec-ext/t1`), two readers, pins and clearances, derive and freeze, register rows, fetch, load,
    tables; then builder options and practice rows. A file that carries the composition table needs one short question
    to the captain (closed lever). Open point: the atlas in use (X002) was never read under the stricter plug-day
    rule; read its series record at the plan run and amend the register row.
  - **Queued, not started:** T1's practice row D (needs its arm written; low committed odds; it runs rather than being
    pruned); the hedge file F and the advance file E only at the captain's specific word.
  - **Tries:** T1 413, embryo 76 (E49 running), heart 80, T3 60. **Scratch hand-backs** (session scratchpad `lift8/`):
    t3x_widen, t3x_stage2, t1_atlas_fix, embryo_shape, k22_row, e49_row.

- **State at 2026-11-20, 16:42 UTC (superseded above where they differ).** Decision log
  2026-11-20 (m) to (v); `docs/S2_RULES.md`, the notes and amendments of 2026-11-20.
  - **Handed to the captain, no score reported yet** (uploads are his acts; paths under the repository):
    T1 `runs/manual_upload/T1_s2rcA.h5ad` (76865826f5a3741c7585f81c4241048ca1a95b5791ebb3c8feb1f2b536678ce4) then
    `T1_s2rcB.h5ad` (c7240cb6a7f7ee93d037f437832efbaad004f1fbe608d80072f8efd5d20a0062); `T1_s2rcC.h5ad`
    (31477ef7df5d5374a74a8249c570217020c593916ab4d3780049b7b3f87b0620) only if A and B each beat R13's file (entry (n)).
    Heart K1 `runs/manual_upload/T2heart_s2trimcopy_half.h5ad`
    (47bb3fe2a082e5b221ee4ae8c003d47df77259867563a11f6c805e9bb82d91d5) then K2 `T2heart_s2iso_half.h5ad`
    (dc436d9225bb28e7997d2ec51e95254cd125e569d44c059596298d0a7b6af5c2); the heart list ends after them.
    Embryo `runs/manual_upload/T2embryo_s2x.h5ad`. T3 `runs/manual_upload/T3_barred_x1.5.h5ad`, then
    `T3_barred_x2.0.h5ad` if it beats the s1 file, else the barred s2.
  - **His reply after a dismissed four-question dialog, verbatim:** "go ahead with whatever will improve our score
    [my reply]". Read (entry (t), labelled as the AI assistant's reading): T3's two Gata4 files C1 and C2: yes
    (row R17); T1's hedge file: not covered, stays out; T1's advance draw: stays declined (closed levers need his
    specific word); the embryo's row-order answer: not withdrawn. No more multi-question dialogs: one short
    question with a default, only when the answer changes what is built.
  - **When T3 scores come back:** record them; after BOTH of entry (k)'s uploads have headlines, fix the base B
    (entry (q)) and build C1 and C2 on `~/v/vec-campaign-c` with `experiments/t3_decomp/build_two_knockout.sh
    <commit> <base>` (identity gates first), then hand over. K21's files, if any qualify, go after C1 and C2.
  - **Running on the DGX:** E48 (embryo row-order practice read; chain in `~/v/scratch/e48`, worktree
    `~/v/scratch/e48-wt` at 9cd4b6e; the run step started 16:24 UTC, embryo count 75 to 76). When it ends: copy
    `experiments/coords/results/e48_row_order/reading.txt` and `.json` from the worktree, commit them with the
    outcome in E48's row quoting only the fixed reading, treat any `files_identical_across_builds` false as void;
    unless reading (i) voids the construction, build the row-order file in `~/v/vec-campaign` itself by the
    commands in the stage-2 hand-back (default must reproduce a03611a3..., then `--row-order copy_parity
    --row-order-test-ok` twice), then hand over. K21 has ended and is recorded (next item).
  - **K21 is read** (decision log (w)): valid, six tries (T3 57), outcome "no file: nothing qualifies"; nothing to
    build from it. The reviewer's two lesser points on its build launcher stay unapplied (no K21 file exists).
  - **The outside-knockout line is parked** (decision log (x)): pipeline merged at 95c7f01 and verified; the
    metadata plan step ran (tables under `~/v/vec-ext/t3/meta` on the DGX, copies of three under
    `~/.cache/vec-ext/t3/meta` locally); nothing fetched. The knockout panel has no per-sample files and no
    in-series control; only six Tbx1 samples could be pinned today. Stage 2 (practice read K22 and its builder)
    is reviewed ready on branch `campaign/t3-outside-stage2` (1c5d964), unmerged. Resume only at the captain's
    word.
  - **Scratch hand-backs** (JSON, in the session scratchpad's `lift8/`): embryo_occ2, t1_candidates, t1_further,
    t3_two_ko, t3_forms, t3x_fix3, heart_candidates.

- **State at 2026-11-20, mid-afternoon (superseded above where they differ).** Decision log
  2026-11-20 (g) to (m); `docs/S2_RULES.md`, the note and the heart amendment of 2026-11-20.
  - **Scores.** The captain uploaded all four files of the next bullet on 20 November and reported their scores
    (`results/manual_uploads.csv`; per-term values in decision log (h), (i)). Heart and T3 rose; T1 and the embryo
    did not beat their first uploads. The board showed our best score, not our latest.
  - **Regime, in his words:** "let's upload all the 4 filees after full optimization to beat the leadeerboard scores
    [2 times a day]". At most two uploads per board per day until 2026-11-26 11:59 UTC, on his word for the
    allowance. A board's candidates are fixed in the decision log BEFORE their scores exist; a returned score chooses
    among them and computes nothing. He also said: no questions on the login or disclosure; do not read the site's
    pages; do not look for our position.
  - **Handed over, no score reported yet:** T3 `runs/manual_upload/T3_barred_x1.5.h5ad`
    (c060fcaa6c7fdc6d6b1a58c8f850a1895a8295abbc50f297de71a32c97f9c420), then by rule `T3_barred_x2.0.h5ad`
    (7393810fbe3be2cde67cf23ed54f45d22da738d7d3c9e124780a8f6281075661) if 1.5 beats the s1 file, else
    `runs/archive/builds_dc7b6ae/t3_build/barred_s2.h5ad` (1b69521e3c28bc84185410e838efc6fcbcff898434e73f9a7547208a4c82711f;
    it carries the same factor as s1, checked, entry (m)); embryo `runs/manual_upload/T2embryo_s2x.h5ad`
    (a03611a343faf1d37edf06ce0e2bccc86cce6f14d3580b823864e5faf3ba03f9).
  - **Building on the DGX.** T1 candidates A, B, C (row R15) on `~/v/vec-campaign-d` at 1420c70, scratch
    `~/v/scratch/t1cand_1420c70` (all.exit 0 means: default reproduced ad8f24c7, each candidate built twice,
    pairs agree). Heart K1 and K2 (row R16; his answer "Both, K1 then K2") on `~/v/vec-campaign` at 1ba7c67,
    scratch `~/v/scratch/heartcand` (done.exit 0 means: the full trim reproduced 04cde636, pairs agree, K2 equals
    dc436d92). After each: copy the upload files to `runs/manual_upload/`, check hashes, local `vec check` by exit
    code only, write `experiments/builds/results/` records, fill the rows, hand over path, sha256, board and a
    short "Data used" text. A test file's shape and size are never printed.
  - **With agents (read their hand-backs from the task output files, never their transcripts):**
    wf_666d097c-21a, embryo stage 2: the run script for practice row E48 (the row-order option; his ruling "Use
    it, after a practice read"), the earlier reviewer's points, and the analysts' call on the check's turned-axes
    reading (the main session is recused on the embryo's occupancy term and on the heart's scope and size).
    wf_e988fff9-6e1, the T3 outside-knockout screen fix; nothing is fetched before the screen passes its
    adversarial review, the register rows exist and the plan tables are read. wf_42122fd0-924, T3 candidates from
    the two released knockouts only (row R17 to come).
  - **Held under the organisers' genotype clause** ("a perturbation that phenocopies them counts as held out"):
    Wls (GSE165300), Hand2 (GSE126128), Mef2c (GSE280543). His ruling on the three tiers does not lift it.
  - **When he reports a score:** a row in `results/manual_uploads.csv`, a decision-log entry, then the board's
    pre-set rule (T3: entry (k); heart: the amendment's ten points; T1: the best of the first upload, R13's file,
    A, B and C stands; embryo: the analysts decide).

- **State at 2026-11-20, 10:17 UTC.** Decision log 2026-11-19 (jj), (kk) and 2026-11-20 (a) to (f). **All four
  last files exist, each as the captain named it, built twice with identical bytes (T3's was already on the
  reference list) and passing the local format check:** T1 `runs/manual_upload/T1_s2rc.h5ad`
  (ad8f24c71c90ca17bce5f8876d503471ab3faa5b17a0878addb00aff689fca4f; record
  `experiments/builds/results/t1rc_5028cab.txt`); heart `runs/manual_upload/T2heart_s2core_half.h5ad`
  (04cde6361c04840ee47fff32eb9a20f4781bd649387cc540ed9fdea71e65fa35; `heart3_3d52b83.txt`; the hold ended by his
  answer "Full trim (as ruled)" to row H27's figures, entry (kk)); embryo `runs/manual_upload/T2embryo_s2xs.h5ad`
  (bc50d012cd2c9481cd18c1f5bcf5b01093e22fd16ac65582a5c9f34321f13115; `embryo3_3d52b83.txt`); T3
  `runs/archive/builds_dc7b6ae/t3_build/barred_s1.h5ad`
  (3712c3f9988d0f8aa4198b22c37b866a84ef8757fd500a9163944f0663a520a2). They were handed to the captain on
  2026-11-20 with a short "Data used" text per board. **Uploads are his acts, by hand, by 2026-11-26 11:59 UTC;
  when he reports a score, add the row to `results/manual_uploads.csv` and a decision-log entry.**
  **His word of 2026-11-20 (entry (b)), verbatim: "don't worry about which login. and disclosure. Focus only on
  improving our score".** Do not raise the login question, the disclosure email or the notes to the organisers
  again; the note drafts stay in `docs/ORGANIZER_QUESTIONS.md`, unreviewed (their three-lens review was stopped).
  **Reads since the rulings:** T63 (T1 at the test's step size: AHEAD, epicardium PAYS, "a placebo gains as well;
  the direction more"; recomputed independently; build as ruled); H27 (heart at forecast sizes; recomputed; he kept
  the full trim); H28 (heart trend step at reduced strengths on the full trim: outcome PLACEBO, no question; not yet
  recomputed independently, a negative result that changes nothing); K20 (T3 factor exponent: written and
  reviewed on `campaign/t3-k20-exponent`, NOT registered and NOT run: the committed record already shows a weaker
  factor scores lower; entry (c)). Tries: T1 413, embryo 75, heart 80, T3 51. `docs/S2_RULES.md` sha256
  a5380cdfbc02 (in `docs/S2_RULES.sha256`).
  **The AI assistant stays recused on the heart's scope and size** (entry 2026-11-19 (bb)).
  **Open, none of it blocking:** the 23 November screened read of mail and the rules page (approved; every thread
  with replies; a reply is never pasted into the session); `tests/test_feng_relabel.py::test_the_dry_runs_refusals`
  reads the checkout's own data folder and fails on `~/v/vec-campaign-d`, where row N55's dry run left its start
  record (make it use a temporary folder); the builder's help text still calls `core_trim` "no ruling";
  METHOD_SUMMARY's `[TO CONFIRM]` and `[TO COMPLETE]` markers on what was sent wait for the captain's own word and
  are not to be asked about; run the complete suite with its temporary folder OUTSIDE any git work tree (the home
  directory is one; `tests/test_boards.py` fails inside it) and delete the folder after (the system disk filled on
  19 November from agents' leftovers).
- **State at 2026-11-19, 21:07 UTC.** Decision log 2026-11-19 (aa) to (ii). **The captain named the last
  file on all four boards** (each an answer to its own priced dialog question; entries (cc) and (ii)):
  heart = the full trim (his own named exception); embryo = the series twin on E7.5's exact positions; T3 = the
  barred branch's s1; T1 = the recomputed-groups atlas file ("proceed with your recommendation" to a question whose
  one recommendation was that file). `docs/S2_RULES.md` carries the first three (sha256 9456a7737140); T1's
  amendment is being written (wf_d0b288ec-143: rules, build set `t1rc`, row R13, report lines, a draft note).
  **Files:** embryo `runs/manual_upload/T2embryo_s2xs.h5ad` (bc50d012cd2c; built twice at 3d52b83, record
  `experiments/builds/results/embryo3_3d52b83.txt`); T3 `runs/archive/builds_dc7b6ae/t3_build/barred_s1.h5ad`
  (3712c3f9988d; on the reference list; no build); heart `runs/manual_upload/T2heart_s2core_half.h5ad`
  (04cde6361c04; `heart3_3d52b83.txt`), **built and HELD**; T1's recomputed file is not built for upload yet
  (check builds ad8f24c71c90 at 8138fa2). **Nothing is handed over.** Before any hand-over: the note that precedes
  each second upload (drafts in `docs/ORGANIZER_QUESTIONS.md`, last sections; into the captain's Gmail drafts only
  at his word, after a three-lens review; he sends them), and the question of which login made the uploads of
  18 November (our key's quota showed 0 of 2 used; he said not to search for the score for now).
  **The AI assistant is recused from every choice and lean on the heart's scope and size** (entry (bb): it began to
  invert the heart's scale score and stopped; no value recorded anywhere). A side note the captain forwarded was
  right that row H26's figures were at an oracle size: **row H27** (the three heart forms at forecast and perfect
  sizes on released stages; designed, attacked, settled, written and reviewed by agents that had not seen the
  session; registered at 4751471) **is running** on `~/v/vec-campaign` (launcher `~/v/scratch/h27/go_h27.sh`; exit
  files under `~/v/scratch/h27/` and `runs/h27/exit/`; then ONE read with `~/v/scratch/h27/go_h27_read.sh`; copy
  back reading.json, reading.txt and the 15 parts and commit them with the row's outcome cell in one commit; then
  an independent recomputation; then the figures go to the captain in the read's pre-set words beside row H26's,
  with no preference from the AI assistant; his ruling stands unless he names another file).
  **T1:** row T62 read and recomputed (the recomputed step AHEAD under both readings: +1.97, +1.67, +2.50, +1.23,
  +2.51; the epicardium step PAYS); the cold check let the advice stand with conditions (entry (hh): the
  compliance price; the test's 2/5 step unread). **Row T63** (the candidate at the test's step size with a placebo
  and the epicardium-held arm; wf_81156c91-e7c) is being attacked, written and reviewed: then the row on main
  first, the run on `~/v/vec-campaign-d`, one read, an independent recomputation. What follows is fixed in its
  row; the T1 file is built (set `t1rc`, twice, on `~/v/vec-campaign-d`) only after that read, or after the
  captain's answer if the read sends him a question.
  **Housekeeping:** the complete suite last passed at 1750ade (3207 passed); rerun it at HEAD once the tree
  settles (worktree `.claude/worktrees/t1-none`, GIT_* unset). `experiments/external/n54_step.sh` must list every
  new test file that touches the builder (its own test fails otherwise). The heart builder's golden-record test
  now tolerates last-bit float differences between machines (ce589a6). 22 junk commits by "t" stay in place.
- **State at 2026-11-19, 19:02 UTC.** Decision log 2026-11-19 (u) to (z). **The labeller ran** (the
  captain's sign-off, entry (u), is the ruling naming N54): rows N55 and N54 are read; settings frozen (sha256
  883e74c2c1ef, min_score 0.3 by N55's rule); route recomputed (unassigned 0.019 and 0.038; four groups with a
  direction; epicardium on 11 early and 32 late cells, minimum 10); the candidate
  `t1 --g-mom --base lib3 --feng-groups recomputed` check-builds twice to ad8f24c71c90 on `~/v/vec-campaign-d`
  (8138fa2; the only checkout with `data/derived/feng_relabel_c/`). T1's last upload is still row R8's file
  (55eaee547ff1) until the captain rules. **The listing of N54 is committed and must stay unopened** by the main
  session until that ruling (it compares with the authors'-group direction). **Since a RECORD routing to recomputed
  is at HEAD, `build_all.sh t1none` is refused** unless `--feng-none-ruling` names a later entry of the captain
  containing "sets the N54 route aside"; the agreement table needs one containing "rules on T1's last upload".
  **In hand:** row T62 (wf_bf253c90-a8c: `experiments/t1_m1/recomputed_read.py`, branch `campaign/t1-recomputed-read`):
  the candidate's practice analogue beside T61's none_gmom on T61's five tasks; on its return a fix round adds the
  two things of entry (z) (a split-half reliability listing per group; an arm with the epicardium group held), then
  the row on main first, the run on `~/v/vec-campaign-d` (launch with `PATH=$HOME/.local/bin:$PATH`; plain ssh
  shells lack `uv`), the read, and the priced question to the captain: the ruled file, the recomputed file, and, if
  the epicardium step is noise, the recomputed file with epicardium held (a builder option and an amendment owed).
  **Waiting for the captain (plain text, no dialog open):** the heart's route: keep the ruled full-scope file (in
  force; the default), his own exception (half or full trim), or leave the direction to the AI assistant (then only
  the ruled file or the half trim, after a second cold check); he said the organisers will not reply, so no
  organiser route. The advance draw of try T61: advised to decline; declined if he says nothing. **Read today:**
  T61 (T1 358), H26 (heart 50). **Closed without a run:** K19. **Stopped:** the late bracket. **Known and left:**
  22 empty "init" commits on main (entry (t)). Agents cannot outlast a 35-minute suite: the main session runs the
  complete suite itself (last: 3156 passed, 2 skipped at 4c57dff; targeted sets since).
- **State at 2026-11-19, 15:02 UTC (newest).** Decision log 2026-11-19 (o) to (r). **Waiting for the captain
  (asked in plain text; no dialog is open):** the heart's route: (a) keep the ruled full-scope file, the default;
  (b) keep it and try the clean route (a screened re-read of the protocol page's two heart rows, and the question
  to the organisers drafted, unreviewed, in `docs/ORGANIZER_QUESTIONS.md`, "Data Protocol question on the E12.5
  heart stack"; if he says yes: three-lens review, then his drafts; a reply is screened and audited, and row H26's
  conditional outcome applies only once he has signed it); (c) his own named exception. The heart upload stays
  held. **Try H26 is read and recomputed independently** (heart 41 to 50): the trims carry over to the copy's
  form; with the truth as released B -1.0 and C -3.9 against A (fold-merged; -7.9 and -15.5 on the one-day
  transition), over the clean truths B +5.8 and C +11.8; the odds a trim needs run from 0.15 to 0.86 by
  transition; the read decides nothing. **Coming to the captain next:** (1) T3's last file, the refused file by the
  rule against the built barred file by the analyst's expected-better-of-two reading, after the cold challenge in
  wf_dc1bec61-a73; (2) the labeller's sign-off: `campaign/feng-relabel-main-fix` at 9f1e047 (the complete suite
  passes: 3070, 1 skipped; settings file sha256 af3ac7224a65, status draft), final check running
  (wf_235124c8-92a); on "ready": merge to main, give him the fifteen-sentence summary
  (`lift8/labeller_port.json` in the session's scratch folder, key fix), then in order: his ruling entry, row N55
  and its dry run on the released E9.5, the freeze with row N54, the relabelling run, the listing, a check build
  of `t1 --g-mom --base lib3 --feng-groups recomputed`, and the priced question on T1's last upload (the built
  no-outside-data file against the recomputed one). **Running:** wf_dc1bec61-a73 (fix rounds and re-reviews of
  rows K19, `campaign/t3-variants-fix`, and T61, `campaign/t1-released-bets-fix`; the T3 challenge). On their
  return: rows on main first, then the reads on idle DGX campaign checkouts (`~/v/vec-campaign` is at b2e84dd,
  idle; `-d` is clean; `-h7` and `-c` hold leftovers). **Stopped:** the T1 late bracket (entry (p)). **Not
  reopened:** the embryo's series twin (level on expectation; the captain may ask for it). The upside wave's
  hand-back is `lift8/upside_wave.json`; the heart read's are `lift8/h26_prep.json`, `h26_fix2.json` and
  `h26_recompute.json`.
- **State at 2026-11-19, 13:31 UTC (newest).** Decision log 2026-11-19 (n). The heart upload stays held. The price
  table of entry (m) was recomputed independently and stands; a cold challenge rejected the advice (the half trim in
  today's form is a mix no rule returned): **the lean on the half trim is withdrawn**; A (the ruled full-scope file)
  is the file in force and the default; B and C go to the captain side by side with no lean. Owed first: a practice
  read of the actual forms, row first (see entry (n)), and the builder option from the wave. Both checks are saved
  in the session's scratch folder as `lift8/heart_scope_check.json`. The next amendment of `docs/S2_RULES.md`
  carries a correction: the half trim was dropped on 18 November without its own priced question.
- **State at 2026-11-19, 13:24 UTC (newest).** Decision log 2026-11-19 (m). **The heart upload is held**: do not hand
  over `runs/manual_upload/T2heart_s2iso_half.h5ad` as the last heart file until the captain has ruled on scope.
  The amendment of 18 November dropped the frozen rule's half trim; entry (m) prices three files as the better of
  two against the banked 46.4 from the committed hedge study (A the ruled full-scope file, B the same form on the
  half trim, C on the full trim). Advice by the exposed-lever rule: B, the pre-set outcome; on C the table and no
  lean; he decides from the released-data table alone or declines. Owed before the question: the builder option
  (`campaign/heart-iso-trim`, from wf_cd5ffbb3-985), its review, an independent recomputation of the table, check
  builds, the drafted `S2_RULES` amendment. The atlas inventory is back and saved in the session's scratch folder
  (`lift8/atlas_inventory.json`).
- **State at 2026-11-19, 13:09 UTC (newest).** Decision log 2026-11-19 (f) to (l). **Ruled today by the captain and
  built:** T1's last upload is the first file's design with no outside data (`runs/manual_upload/T1_s2none.h5ad`,
  55eaee547ff1; row R8; `experiments/builds/results/t1none_385ff58.txt`); the embryo's last upload is the first file
  on E7.5's exact positions (`runs/manual_upload/T2embryo_s2x.h5ad`, a03611a343fa; row R9;
  `experiments/builds/results/embryo2_2dccfee.txt`). Both pass the local format check; `docs/S2_RULES.md` carries
  both amendments (sha256 ce5a4584c9b2). The heart's and T3's last files are as built on 18 November. **In his
  Gmail drafts, not known sent:** the disclosure of the X002 labels and the amended sources note; the two short
  notes that precede T1's and the embryo's second uploads are drafted in `docs/ORGANIZER_QUESTIONS.md` and go into
  drafts at his word. **His standing instruction (19 November):** organisers are not expected to answer; assume the
  rational reply, as odds with the loss if wrong, and change no file on an assumed reply. **Uploads:** the Rules page
  and the server say two official uploads per board (not unlimited); the team key shows none used and lists no
  upload, and the overall board has no row under the team's name; the row "Star" is not his; he said not to search
  for the score now (entry (l)), so no screened read for it runs. The AI assistant's advice stands: no upload from
  any login until it is known which account made the uploads of 18 November; raise it when an upload is next
  discussed. **He approved the two outside atlases as T1 model inputs** ("utilize the atlases": permitted stages,
  raw counts and stage labels, groups recomputed; register rows X002, X003). **Running:** wf_cd5ffbb3-985 (the
  upside wave: T3 variants on `campaign/t3-variants`, T1 released-data bets on `campaign/t1-released-bets`, the
  heart's scope hypothesis on `campaign/heart-iso-trim`, a T1 late-bracket design with a red team, a pricing table
  of each board's last file by the expected better of two, reviewers); wf_5a6524a5-2ff (the labeller ported onto
  main, `campaign/feng-relabel-main`, with a review); wf_31229ec2-5eb (inventory of the atlases on the DGX,
  outputs under `~/v/scratch/atlas_inventory/`); the complete test run at 2dccfee in the worktree
  `.claude/worktrees/t1-none` (log `t1none/full4.log` in the session's scratch folder). **On each return:** read the
  hand-back (never a transcript); a fix round where a review says not ready; ledger rows on main BEFORE any run or
  reading; practice reads on an idle DGX campaign checkout (`~/v/vec-campaign` is at 2dccfee, clean, idle; rsync,
  then checkout); apply each row's pre-set statement; bring each candidate to the captain as a question priced both
  ways. Every change of an upload is his ruling with an `S2_RULES` amendment first. The reference checkout is
  untouched at dc7b6ae. Screened read approved for 23 November (every thread).
- **State at 2026-11-19, 09:15 UTC (newest).** Decision log 2026-11-19 (b) to (e). The Feng authors' labels were made
  with every stage pooled, so both T1 files rest on labels the organisers ask teams to leave out. `builds/T1_s2.h5ad`
  is never uploaded, and no later T1 upload carries its 0.60 (it cannot be recomputed on permitted data). **Waiting
  for the captain:** (1) T1's last upload: the AI assistant advises the first file's design with the atlas step
  taken out, `t1 --g-mom --base lib3 --feng-groups none` (no outside data; priced in entry (d), corrected in (e));
  the alternative, groups recomputed by a labeller, is parked on `campaign/feng-relabel-fix2` (final review not
  ready: three blockers, nine lesser points, a data-folder link problem on the DGX, corrections owed to its marker
  sources) and needs about a day and several sign-offs; (2) his word to put the two rewritten emails
  (`docs/ORGANIZER_QUESTIONS.md`: "Disclosure of the X002 labels" and the amended sources note of 18 November) into
  his drafts after he has ruled; the copy of the sources note now in his drafts is not to be sent as it stands;
  (3) the embryo's last file (2026-11-18 (v)). **Ready for his T1 ruling:** the option on `campaign/t1-none`
  (ea95c02; review ready; lesser points being applied on `campaign/t1-none-fix`, wf_79afa627-63b: on its return read
  the review, rebuild on the DGX campaign checkout and confirm the no-direction file is still 55eaee547ff1 and the
  two option-off files still have the reference hashes, run the complete suite); the ruling package as drafts on
  `campaign/t1-none-docs-fix` (1921e9f) with its ordered steps from the ruling to the hand-over, saved in the
  session's scratch folder as `lift8/t1none_package.json` (key package.fix: the draft amendment, row R8, the X002
  and X003 texts, the method-summary texts, the runbook lines, the DATA_POLICY note, the checklist). At the ruling:
  his words into a new decision-log entry; merge the option; squash the docs branch in uncommitted; take "DRAFT, not
  in force" out of the amendment and put his words in; recompute `docs/S2_RULES.sha256`; one commit with row R8 and
  the texts; build the set `t1none` twice at 16 threads on a DGX campaign checkout (never the reference checkout;
  set BUILD_OUT outside any linked folder); compare with 55eaee547ff1; copy the file here, check it, hand it over
  with the checklist and the note that precedes each second upload. The DGX campaign checkout `~/v/vec-campaign` is
  at ea95c02, clean; check builds (never upload sources) are in its `runs/builds_check_t1/` and `runs/t1none_ea95c02/`.
  The reference checkout is untouched at dc7b6ae. Screened read approved for 23 November (open every thread).
- **State at 2026-11-19, 04:46 UTC (newest).** Decision log 2026-11-19 (a). HOLD: `builds/T1_s2.h5ad` is not uploaded
  as built (an organiser's reply of 2026-11-16 asks teams to leave out annotations learned using protected cells;
  X002's groups are the authors' labels). Waiting for the captain: (1) what T1's last upload should be once labels
  are recomputed from permitted cells (the hedge rebuilt, or the first file's design rebuilt) and whether to ask the
  organisers about the uploaded first file; (2) the embryo's last file (2026-11-18 (v)). Running: wf_c9186108-9a8
  (provenance of the authors' labels, a relabel design, a challenger); on its return pick the option, have it
  implemented on a campaign branch with reviews, register its N row, run the target-free listing, and bring the
  captain the table. T3's severity question is answered (skill 0.5; no follow-up). Screened Slack checks open every
  thread with replies from now on. Sources are disclosed in the submission page's method-report fields with the
  report attached (the organiser's words): tell the captain at hand-over and add it to the checklist.
- **State at 2026-11-18, 22:09 UTC (newest).** The departures section of `METHOD_SUMMARY.md` is merged (a3a9911); its
  [TO COMPLETE] line and item 2's "not yet in the record" (whether and when the disclosure email and the three notes
  were sent) are filled when the captain says. If he picks embryo option 1, 2 or 3, the disclosure email draft in his
  Gmail no longer names the right embryo file and is corrected before it is sent (three-lens review again). Nothing
  runs; the bullets below hold the rest.
- **State at 2026-11-18, 21:46 UTC (newest).** Decision log 2026-11-18 (w). All of today's rows are read (H25, E47a,
  E47b, N53, K18, T60); nothing runs on the DGX. Waiting for the captain: the embryo's last file (entry (v); the
  bullet below lists what follows his answer). In review: the method summary's section on the departures
  (wf_e4e3dd12-c83; branch `campaign/method-report-departures-fix` when done): read its final check, apply the
  lesser points, merge, and complete its marked line when the captain has ruled on the embryo. Beside the rule's
  "no question" on T3, entry (w) gives the captain the four practice cells of try K18 (the barred form's mean is
  higher, 60.2 against 57.8, with one failing cell); the swap is his to order and the advice stays with the rule.
- **State at 2026-11-18, 20:56 UTC (newest).** Decision log 2026-11-18 (q) to (v). The captain holds the second
  uploads and asked for a week of research. Read so far: H25, K18, E47a, E47b, N53 (rows of those names; results
  under `experiments/harness/results/h25_heart_forms`, `e47_floor_exact`, `n53_floor_proximity` and
  `experiments/t3_decomp/results/k18_built_files`). Waiting for the captain: which file is the embryo's last upload
  (entry (v): `fgw_shape0x`, recommended; `fgw_shape0x_series`; fgw_shape0; or the ruled `fgw_aligned_series`). On his
  answer, if it is not the ruled file: amend `docs/S2_RULES.md` with his words and the new hash, add the command to
  `experiments/builds/build_all.sh` with its tests, build twice on the DGX at two thread counts and compare, write
  the reference list entry, copy the file here, check it locally, draft the note to the organisers (three-lens
  review) and hand it over. Running: try T60 on the DGX (`~/v/vec-campaign-h7` at 13fd89a; exit files
  `runs/logs/t60_lib{1,2,3}.exit`; then copy `experiments/t1_m1/results/t60_reach/tasks/lib{1,2,3}.json` back, commit
  them unread, run `reach_read.py read --out t60_reach`, record the outcome by row T60's statements). Check builds
  (never upload sources) sit in `~/v/vec-campaign/runs/builds_check_emb/` at 9967cf1. The DGX reference checkout is
  untouched at dc7b6ae. Not started, with reasons in entry (s) and the status messages: the learned models, the
  heart's tilted half, T3's cell-resolved factor. Owed: the method report's lines on every departure (now also: rule
  9 was shaped after the first scores; K18's reading machine; the research week's tries after the first scores).
- **State at 2026-11-18, 17:52 UTC (newest).** The captain holds the second uploads ("We don’t need to submit again
  right now") and asked for a week of research; decision log 2026-11-18 (q) to (s) has his words, the reading, the
  design round and the portfolio. The four files of the bullet below remain the default second uploads. Running:
  wf_56b120a7-a36 (embryo arm fgw_shape0x with tries E47a, E47b and reading N53; T3 try K18) and wf_d1838b6f-23d (heart
  try H25; the second-file value table; rule 9 and the pre-upload checklist as text). When each returns ready: merge
  its branch to main, add its ledger row on main, sync a DGX campaign checkout (nothing may run in it), start the run
  detached with a local background waiter, read by the row's fixed rule, and take the outcome to the captain as a
  priced question where the row says so. Nothing changes an upload without his ruling. Still owed from before: rows
  in `results/manual_uploads.csv` when he uploads; the method report's lines on every departure; his nomination of
  final entries before 26 November 11:59 UTC. Saved outside the repository: the design round's result and the wave's
  hand-backs, in the session's scratch folder `lift8/`.
- **State at 2026-11-18, 16:54 UTC (newest).** All four second-upload files are with the captain, who uploads by
  hand: `builds/T1_s2.h5ad` (f169ea492e8d), `builds/T2embryo_s2.h5ad` (41982d8afd78),
  `runs/manual_upload/T2heart_s2iso_half.h5ad` (dc436d9225bb) and `t3_build/refused_s1.h5ad` (0394575a11fc); decision
  log 2026-11-18 (m) to (p) has the amendment, the heart build, and tries K16 (void) and K17 (no request). Order for the
  captain: send the disclosure email, upload T1 and the embryo, send the heart note and upload the heart, send the T3
  note and upload T3; the scorer-build note is optional and goes before the embryo upload. When the captain pastes the
  new scores: add rows to `results/manual_uploads.csv`, a decision-log entry with the scores as pasted, and the method
  report's lines on the departures (the gate, T3's order, the heart's halving and its new second file, second uploads
  chosen with the first scores known, tries K16 and K17 after the first scores), reviewed like its exposure section.
  The DGX reference checkout is untouched at dc7b6ae; the campaign checkouts are at cccff68 (`vec-campaign`) and
  307c9b2 (`vec-campaign-h7`), idle. Nomination of final entries on the site is the captain's, due 26 November
  11:59 UTC.
- **State at 2026-11-18, 16:00 UTC (newer than the bullet below).** `docs/PLAN.md` decision log 2026-11-18 (k) and
  (l) hold the diagnosis, the plan and the captain's two answers. Second uploads: T1's and the embryo's built s2
  files were handed to the captain (hashes f169ea492e8d and 41982d8afd78); T3's last file waits for try K16 (the
  fallback is the built `t3_build/refused_s1.h5ad`, 0394575a11fc); the heart needs a new file (E10.5 unchanged,
  scaled evenly to the size file's raw-mode H1 size, a uniform random half), which the captain's go-ahead covers.
  Before that file is built: merge the reviewed builder options, append the drafted amendment to `docs/S2_RULES.md`,
  write its new hash to `docs/S2_RULES.sha256` and the decision log, add the set `heart2` to
  `experiments/builds/build_all.sh`, then build twice on the DGX's campaign checkout (two thread counts), check,
  halve, copy here, and give the captain the path and hash with a one-line note for the organisers. The captain
  ordered that organiser replies be assumed and nothing wait on them (entry (l) has the assumed replies, the odds
  and what each costs if wrong). Two drafts are in the captain's Gmail (the disclosure email; the scorer-build
  note). Owed after the uploads: the method report's lines on every departure, reviewed like the rest of its
  exposure section; the captain's nomination of final entries on the site before 26 November 11:59 UTC.
- **Read this before anything else (2026-11-18, 15:13 UTC): the first official uploads are in, by the captain's
  hand, and the second uploads are being redesigned at the captain's direction.** `docs/PLAN.md` decision log
  2026-11-18 (j) holds everything: the files and hashes, the per-metric scores as the captain pasted them, each
  departure from the written procedure, the captain's direction verbatim and the AI assistant's reading of it.
  Totals: T1 70.9, heart 46.4, embryo 64.6, T3 47.8, Overall 174.2. One upload is left on each board; the deadline
  is 26 November 11:59 UTC. What holds now: (1) the second uploads are no longer fixed by the pre-score rules of
  `docs/S2_RULES.md`; they may be chosen with the first scores known, by the captain's say-so of 18 November; (2)
  the Official Rules' line stays: scores choose among predictions and never compute one, and no property of a hidden
  target is solved from scores; (3) the data rules and the never-open list stay; (4) a candidate that touches the
  post-exposure block (cell numbers, sizes, scope; T3's composition; the heart's depth; T1's depth) goes to the
  captain as its own priced question; (5) every upload is the captain's, and the captain uploads through the site's
  form, so the client's guards and ledger do not see them (`results/manual_uploads.csv` is the record). Owed by the
  AI assistant: the disclosure email (due before the first upload; draft for the captain to send), the amendment to
  `docs/S2_RULES.md` with the captain's words and a new hash before any new second upload is built, and the method
  report's lines on the departures. Everything below this bullet that speaks of dates, the gate blocking uploads,
  or 'nothing runs' is of before the uploads.
- **Read this first (2026-11-17).** The organisers moved the final deadline to **25 November 2026, 23:59 AoE (26
  November 11:59 UTC, 17:29 IST)**, for uploads and final-entry marks alike. Every later date in this file, in
  `docs/PLAN.md` section 7, `docs/P3_RUNBOOK.md`, `docs/S2_RULES.md` and `AGENTS.md` is superseded in fact but not
  yet in text: the dates are the captain's to move, and the captain has been asked (proposal: freeze and first round
  21 November from 00:00 UTC, second round 23 November, T3's last upload 24 November, marks 25 November by 12:00 UTC,
  any file change rehearsed by 20 November). **The upload gate's exception (27 November) is after the deadline: until
  the captain moves it the client refuses every official upload.** When the captain names the dates, one commit
  changes `results/upload_gate.json` (with the captain's words), its test in `tests/test_submit.py`, an amendment at
  the end of `docs/S2_RULES.md` with a new `docs/S2_RULES.sha256`, the runbook tables and `AGENTS.md`. Decision log
  2026-11-17 (a) has the whole list. The organisers also answered the exposure notice (accepted; they ask for the
  explanation and the precommitted procedure in the method report): a further reply to the thread's other questions
  is announced and must go through a screening agent, never be pasted.
- **What is live (2026-11-18, 09:24 UTC; the paragraph below is of 07:38 UTC, and the captain's forwarded notes on g_mom at the test were examined after it: see item (4) of the captain's list).** No candidate change is live and nothing runs. **The captain cleared
  the session's goal on 18 November (no words): no new read or agent round is opened without the captain's word.**
  The ninth independent pass is complete (decision log 2026-11-18 (a) and (c): no result wrong, no open line). The
  reads it named were finished and recorded (entries (d) to (f)): K15 (nothing changes), T58 (the first ground
  against a larger step withdrawn in part; the depth ground stands), T59 (both pre-set conditions for a request on
  g_mom's damping met; a cold check then said stop, after two target-free listings showed the test pair's factors
  to be larger than the practice fit's and to hold a library component the practice fit lacks). H24 and E46 were
  not run: rows H14 and E29 had read them. The DGX holds four campaign checkouts, all idle and clean
  (`~/v/vec-campaign` at 6214e28, `-h7` at 50154a0, `-c` at 36c4d3d of `campaign/t3-f-fix2`, `-d`); the reference
  checkout is untouched at dc7b6ae. Built, reviewed and unused on campaign branches: `campaign/t1-n` (2921b96, a
  step fraction), `campaign/t1-w` (f51a278, the form `cn`), `campaign/herd-t2-fix2` (c0b4aa0, a herded heart
  entry). No builder option exists for g_mom's damping. The schedule amendments are drafted on
  `draft/schedule-nov21-fix` (490a51a, reviewed: ready, not merged; decision log 2026-11-17 (f)): fill its markers
  with the captain's words, recompute the rules hash, rebase and merge. Owed by the AI assistant when the captain
  next rules on anything: the amendments above, and one more dated sentence in the method report's exposure
  section on tries T58 and T59 (reviewed like the rest of that section).
  **For the captain, in this order:** (1) the dates: the gate's exception (27 November) is after the organisers'
  deadline, so the client refuses every official upload until the captain moves it (proposal: freeze and first
  round 21 November 00:00 UTC, second round 23 November, T3's last upload 24 November, marks 25 November by 12:00
  UTC); (2) the six confirmations for the method report (next bullet); (3) three priced questions, none closed by
  the AI assistant, each with its advice 'change nothing' labelled as its own and formed after the results: (a)
  Task 1's step, 0.60 in place of 0.40 on both uploads, the output of the rule fixed before the runs (0.69 without
  a cap that leaned on a statement its row barred): on Task 1's board, by the truth's depth against the cells the
  prediction is built from (0.3, 0.6, 0.9, 1.35, 1.8), +1.98, +1.21, +0.80, +0.25, -0.26 at a real share of 0.28
  and +2.16, +0.97, +0.30, -0.52, -0.96 at 0.22; try T58 adds that most of the gain beyond twice the step goes
  with Feng's own genes; the captain has seen figures that bear on that depth, so the captain decides from the
  table alone or declines, and declining is the default (decision log 2026-11-18 (a), (f)); (b) g_mom's damping,
  1.0 in place of 0.5: for, +1.47 and +1.36 on the two practice cuts, +0.71 to +0.90 beyond a twin without
  development and a permuted control, and +1.83, +1.65, +1.27, +0.68, +0.27 in five worlds from a shallower to a
  deeper truth; against, the covariance term loses in the deepest world, nearly half the gain comes with factor
  sets that carry no development, and at the test the factors are 1.09 to 1.48 times as large and hold a library component (0.30 of the cardiomyocytes' factors' squared size) that no
  practice run contains; a cold check says stop (entry (f)); (c) Task 3's factor at half size: +0.89 and +2.35 on
  the two cross-knockout folds with the factor alone, -0.55 and -1.67 on the two replicate folds in the file's own
  form (row K15); (4) g_mom's fit at the test (decision log 2026-11-18 (g) to (i)): its factors go with which library a cell state's cells came from (0.31 of the cardiomyocytes' factors' squared size); a uniform offset between libraries explains about 1 to 3% of the factors' size, the rest is a difference between states that released stages cannot split into development and identity; three options, none scorable: A the files as built, B the fit on the dissected-heart library alone on the first upload (builder option reviewed and ready on `campaign/t1-ac-fix3`, 86d19ff, off by default, full suite 2,672 passed on the DGX), C no g_mom step; advice, the AI assistant's own, A, as a default and not a tested result; (5) the older open decisions (Task 1's read depth: the
  rule gives declined; the embryo's moved-less family and the late-fetal practice family: advice no). Once the
  captain names the dates, the AI assistant owes the amendments (the gate's file with the captain's words and its
  test, an amendment at the end of `docs/S2_RULES.md` with a new hash, the runbook's tables, `AGENTS.md`), the
  freeze-morning rebuild and the disclosure email.
- **The method report's exposure section (2026-11-17).** `METHOD_SUMMARY.md` now holds "Exposure of hidden-target
  metadata on 2026-11-04 and the precommitted procedure", which the organisers asked for. It was drafted by an agent
  from the dated records, reviewed by three others (accuracy, framing, tone) and revised once; the AI assistant
  checked the commits it cites. Before the report leaves the repository the captain confirms what only the captain
  can: (1) what the captain saw and did on 2026-11-04 (the address opened in a browser, the index and then the
  signed-in Data page pasted, about 07:05 UTC) and that no copy exists outside that session's local log; (2) that the
  captain's choices of 2026-11-05, 2026-11-06, 2026-11-11 and 2026-11-13 (rows R4 and R7) did not draw on the figures
  seen; (3) that the email of 2026-11-04 went out as recorded and that the quoted answer is the one received; (4)
  whether the correction on Task 3's count default (set in code on 2026-11-01, after the release; the email said
  "before") should reach the organisers through this section, as it now does; (5) whether to answer the open question
  on Task 1's read depth first (the section says it is unanswered and that the rule gives "declined"); (6) whether
  older commits that still hold the Task 1 E12.5 file size are to be rewritten (a force-push; the captain's call).
  Not traced to a dated entry by the drafter and left as the record's own summary: who pruned the Task 1 additions
  of 2026-11-04, the date the chamber-only scope file was fixed, and the reason the embryo's depth is outside the
  block. Section 5 of the brief (what the organisers were told) quotes their answer in full.
  One bullet was added to the section after its review, on 2026-11-17: the step reads of that day (tries T44 to
  T54) and that no request followed. It has not been through the three reviewers.
- Goal: a compliant Human-track entry on all four test boards that finishes as high as it can. The final phase now
  ends on 2026-11-25 at 23:59 AoE, with two scored uploads per board for the whole phase; none is used as of
  2026-11-17.
- Current status (2026-11-14): the entries are built and their methods chosen. T1: M1 with R1, R2 and R3b, its first
  upload drawn from E10.5's dissected-heart library since 2026-11-13 and its second on the pooled base at 60%
  cardiomyocytes. Heart: B0 with the pre-registered size and, since the captain's ruling of the evening of
  2026-11-13 (row H11), each lineage's expression step taken from its trend over all released heart stages.
  Embryo: fgw, with `fgw_aligned_series` as second upload
  since 2026-11-13. T3: the frozen rules' files, the second upload now with a wild-type-only factor. Nothing is
  uploaded. The captain's upload gate holds every first upload until the team's own estimate passes 209 (it is 184,
  range 170 to 198) or the deadline exception opens at 2026-11-27 00:00 UTC; each upload then needs the captain's
  go-ahead at the time, and the captain runs the commands.
- No method ruling is open (`docs/PLAN.md` decision log 2026-11-13 (d), (k) and (w); `docs/S2_RULES.md`, its last
  three sections). After the captain's rulings of the evening of 2026-11-13 (T1's second upload on the mid base,
  row R7; both heart files with the trend, row H11) the ten files were rehearsed a third time on the DGX (decision
  log (bb)): three byte-identical builds, and **the freeze-morning reference is now
  `runs/archive/builds_dc7b6ae/SHA256SUMS` there** (`experiments/builds/results/rehearsal_dc7b6ae.txt`,
  `builds_dc7b6ae.SHA256SUMS`; `builds_38616fe` and `builds_fd2683a` are superseded and hold the replaced heart
  files). Every file of the set has passed the server's format check (six on 2026-11-13 at the captain's word,
  decision log (x)). The upload client refuses a scored upload of a file that is not on this list (guard
  "reference list", pointer `results/reference_list.json`; decision log (cc)).
  What is owed:
  - DONE on 2026-11-14 (decision log 2026-11-14 (c), item 1: the family closed with no request; kept here as the
    record of what was owed): the controls and checks of tries E41b. The embryo arm `fgw_partial` (branch `campaign/emb-f`, fcc1f6b)
    passed the adoption rule and the first-upload pair rule on 2026-11-14 (ledger row E41;
    `experiments/harness/results/n45_embryo_f_readings/`). Before any request to the captain: (1) build, review and
    read the three control arms the captain asked for (every anchor moved by the same share; the same number left
    unmoved at random; no movement), row first, on the nine practice tasks; the brief is drafted (session scratchpad,
    `lift6/briefs/embg.md`) and nothing of it is built; (2) rerun the cold adversarial check and the independent
    recomputation of the reading's tables, both of which stopped on the account's usage limit on 2026-11-14 without a
    result; (3) then a priced request, both ways, for the simplest arm that passes, or a recorded decline. A yes needs
    a builder port, a dry build, a fourth rehearsal and a new reference list by about 24 Nov. Then the seventh
    independent pass over decision log (ii), lines 13 to 18;
  - nothing to the organisers: on 2026-11-14 the captain had their logical replies assumed for the questions still
    open (decision log 2026-11-14 (a); `docs/ORGANIZER_QUESTIONS.md`, "Working answers (2026-11-14)"). The replies
    confirm the files as built; the items below about the eighth email, the Task 3 batch question and the screened
    reads stand only as the captain's options and as checks;
  - the captain's decision on a drafted question, if the captain wants it sent: an eighth email asking the organisers
    how final entries are marked (per board or one pick over several boards; what counts on an unmarked board; until
    when a pick can change). It is reviewed and staged in `docs/ORGANIZER_QUESTIONS.md` with its readings fixed
    beforehand, and a draft with placeholders for the name and the team is in the captain's mailbox. The per-board
    answer is a working answer rated 0.7 (the table of 2026-11-10); the estimate assumes it;
  - nothing on request R-D unless the captain wants it (`docs/PLAN.md` decision log 2026-11-13 (ff) and (gg); rows
    N37, N38 and N39 of `experiments/harness/TRIES.md`). The request: change the embryo's first upload from `fgw`
    to `fgw_trend_render` (the same cells and positions, each cell's expression moved to the profile a trend
    through the four released embryo samples gives at E7.75), second upload unchanged; tabled at +0.93 embryo
    points at the middle weighting by a pre-set pair rule on the nine practice tasks. **The AI assistant advises
    no**: a cold check and reading N39 showed that the gain follows the fit geometry and that the test's geometry
    is the losing kind where the truth is like its earlier bracket (fitted gain -0.60, uncovered). The request
    stays on record as the rule's outcome; if the captain rules for it all the same, it needs the builder port
    (prepared on the branch `port/emb-s1-trend`, off main, and dry-built once at 394f59f with exit 0, decision log
    (ii); a new dry build only if the port's commit changes), a fourth rehearsal, a new reference list
    and pointer, the server's format check of the new file at the captain's word, and the summary and disclosure
    lines;
  - nothing on the runner-up either unless the captain wants it: `fgw_aligned_trend` also passed the pair rule as a
    first upload (+0.59 embryo points at the middle weighting, reading N37). Reading N40 found its gain follows the
    same fit geometry (fitted -0.53 at the test's earlier-bracket geometry), and the AI assistant advises no; that
    lens was fixed after the scores existed, so it is advice and not a pre-set gate (decision log (ii));
  - on 19 and 23 Nov, the approved screened reads for an organiser answer to the seventh email
    (`docs/SCREENED_READ_SEVENTH_EMAIL.md`; a reply in the captain's inbox is not pasted into the session). One
    consequence was fixed on 2026-11-13 (decision log (hh)): if the answer is that the hidden embryo sample was
    hybridised and imaged with E7.5, and final marks are confirmed to be per board, an earlier-side first upload
    for the embryo goes to the captain as a priced request, after a target-free check of the test pair's factor;
    on any other answer, or none, nothing changes;
  - about 19 Nov, the Slack follow-up draft on floored T3 severity if no organiser has answered; with it, for the
    captain to keep or cut, the Task 3 batch question that was cut from the seventh email in review (was the
    beta-catenin knockout processed with the Gata4 replicates): an answer changes no file and would move the
    estimate by about 7 points either way (decision log 2026-11-13 (ee));
  - a new rehearsal and reference by 24 Nov only if a build changes after 2026-11-13 (and then the pointer
    `results/reference_list.json` is rewritten for the new list);
  The generator tests the captain ordered on 2026-11-13 are done (decision log (o) to (r); rows T34, T35, H8, H9,
  H10; results `experiments/harness/results/gen_t1_v1/` and `gen_heart_v1/`): five counted reads, nothing kept, no
  upload, build or rule changed. The code is on the branch `campaign/gen` (pushed; 459916f), off main. The heart
  tries that followed (rows H11, H11s, H12, H12s; branch `campaign/trend`, pushed): `b0_trend` passed the adoption
  rule and was ruled into both heart files; its size-matched twin and the two parts of its factor read alone were
  not kept. The sixth round that night (decision log (dd), (ee) and (ff)) kept nothing by the adoption rule: its
  first part (rows N33, N28, N29, N32, H13, K9, K10) ran three heart arms and two T3 readings and restated the
  forecast for the files as built (mean 194.4, the upper reading; the gate figure stays 184); its second part read
  every remaining allowed line (26 tries in 14 rows: H14 to H20, E29 and E31 to E35, T37, T38, K11) and produced
  request R-D above, which the AI assistant advises against. The options for taking our own estimate past the
  leading total were declared used up in decision log (ff) and (gg), too early: the round's audit ((ii)) then read the
  eight lines that independent passes over the record found open and the closures they found resting on judgement
  (tries T39 to T41, K12, K13 and E36 to E39; reading N40), and replaced every closure that rested on the AI
  assistant's own odds by a row-first read, a pre-set gate or arithmetic, or left it to the captain: nothing passes,
  nothing changes, and two critics found no open lever. Tries: T1 37, embryo 60, heart 25, T3 27.
  Settled on 2026-11-13 (decision log (k)): T1's second upload (the mid base; built in and rehearsed), the heart
  try (run, not kept), and the checked copy of the reference files on this machine (kept, by the captain's word).
- Branch / commit: `main` on GitHub (private). The DGX reference checkout
  (`PrakashDGX_H2:~/v/virtual-embryo-challenge`) is at dc7b6ae, the third rehearsal's commit, and stays untouched
  until the freeze-morning sync; campaign runs use `~/v/vec-campaign` and `~/v/vec-campaign-h7`, which sit at
  campaign-branch commits after the sixth round. The round's candidate code is on pushed branches named in decision
  log (ff), and the audit's on `campaign/t1-c`, `campaign/t1-d`, `campaign/t1-e`, `campaign/emb-c`,
  `campaign/emb-polar`, `campaign/emb-mirror`, `campaign/emb-d`, `campaign/t3-c` and `campaign/t3-d`; none of it is
  on main.

## Completed This Session (2026-11-12 to 2026-11-13: the campaign the captain asked for on 2026-11-12)

- [x] Sixth round, second part (2026-11-13, late night; decision log (ff)). Every remaining allowed line with a
      counted read: 26 tries in 14 rows on arms built and reviewed on campaign branches (heart H14 to H20, embryo
      E29 and E31 to E35, T1 T37 and T38, T3 K11), gate N34, readings N37, N38 and N39. Nothing passes the adoption
      rule; eight lines could not be built from the record; request R-D came out of reading N37 and the AI
      assistant advises against it after a cold check and reading N39 (decision log (gg)). A second option, an
      earlier-side file as first upload, was attacked cold before any advice and declined (decision log (hh)).
- [x] First and second wave (2026-11-12; decision log 2026-11-12 (a) to (c)). The finding: released samples differ
      more by how they were measured than by stage (diagnostic N23; the T1 library census in `docs/DATA_NOTES.md`).
      Embryo tries E19 to E25 (depth and run renderings, pair pricing), T1 tries T27 to T30 (library kind and depth
      at zero time gap), T3 tries K6, a second research pass (`experiments/research_2026_11_runs/REPORT2.md`), and
      requests R-A, R-B and R-C.
- [x] Third wave (2026-11-13; decision log (a)). An independent recomputation of the record: every number held, the
      presentation was corrected. Tries E26 to E28 and T31, diagnostic N24. R-A restated with three options, priced
      on nine practice tasks.
- [x] Fourth wave (2026-11-13; decision log (b) and (c)). The R-C email, drafted after a three-lens review and sent
      by the captain. Reading N25 and T3 tries K7 and K8: a per-gene factor from wild-type samples alone, closed by
      its pre-set reading, with both readings given to the captain. Both R-A candidate files reproduce at 16 BLAS
      threads (`experiments/builds/results/determinism_1f8956b.txt`), and R-A's second option has a build command
      (`build_test.py embryo --method fgw_shape0 --expression render_hi`).
- [x] The captain's rulings of 2026-11-13 (decision log (d)): the embryo's second upload, T3's second upload and T1's
      first-upload base changed (rows R4, R5, R6), with `docs/S2_RULES.md` amended, two build options, two independent
      code reviews, and dry builds with one sha256 per new file. No upload. Tries: T1 28, embryo 36, heart 4, T3 15.
- [x] The DGX rehearsal of 2026-11-13 (decision log (e)), at the captain's instruction: the checkout synced to
      fd2683a, all ten files built twice from empty folders (byte-identical, clean stamps, `vec check` passed, seven
      equal to the archive of f337dc6), the set archived as the freeze-morning reference, and the freeze-morning
      command tried against it. T1's new s1 and the embryo's new s2 have the same sha256 on both machines; T3's
      refused s1 with the factor does not (the draws of 84 cells; a nearest-neighbour step follows the BLAS kernel),
      so only DGX hashes check DGX builds. No server format check, no upload.
- [x] A copy off the DGX (2026-11-13, after a side agent's note the captain forwarded): this machine holds a checked
      copy of the reference set (`runs/archive/builds_fd2683a/`; since the second rehearsal
      `runs/archive/builds_38616fe/`) and of the T3 run folder
      (`runs/t3_decomp/gata4_20261028T122425Z/`), owner-only. The client's build-stamp guard and local format check
      pass on all ten copies (nothing sent), so an upload can go from the copy if the DGX is down, at the captain's
      word (`docs/P3_RUNBOOK.md`, "A copy off the DGX"). Every file also builds here now: seven of ten byte-identical
      to the DGX's, the three T3 factor files with a few cells drawn differently (3, 84 and 14 cells).
- [x] The harness files brought up to date a second time on 2026-11-13 (decision log (f)), at the captain's request:
      `AGENTS.md` (final-phase reading list, how a test build is read, where an upload's file comes from, the DGX
      sync rule, the repository map), `docs/COMPUTE.md` ("What reproduces where"), `README.md`, `init.sh`,
      `ASSUMPTIONS.md` (A28), `feature_list.json`, `METHOD_SUMMARY.md`, ledger rows R4 to R6.

- [x] The campaign's fifth round (2026-11-13; decision log (g) to (j)), after the captain renewed it through the
      ambition-lift skill. A fourth screened leaderboard read (sixth place 203.00, 124 Human entries); a cold
      review that stopped the first plan; a forecast of the total as a distribution
      (`experiments/gate/forecast.py`); tries T32 and T33 on T1's pair, both recomputed independently; the heart's
      target-free gate N26, which failed, so no heart try ran. Tries: T1 30, embryo 36, heart 4, T3 15. Campaign
      runs on the DGX now use a second checkout (`~/v/vec-campaign`); no GPU job ran.
- [x] The captain's answers of 2026-11-13 (decision log (k) to (m)) and what followed. T1's second upload moved to
      the mid base (row R7): builder option with tests and an independent code review, rules amended in the
      captain's words (new sha256 23ae5197...), `build_all.sh` changed, all ten files rehearsed on the DGX at
      38616fe (three byte-identical builds; nine equal to the first rehearsal's, T1_s2 new), the reference list and
      the checked local copy refreshed. The heart try H7 with its control was run and not kept (heart 4 to 6
      tries). The checked local copy is confirmed.

## Completed Earlier (2026-10-27 to 2026-10-28)

2026-11-13: kept as history; the sessions of 2026-10-29 to 2026-11-11 are logged in `progress.md`.

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
| Lint and tests (2026-11-13) | `uv lock --check`; `uv run --frozen --no-sync ruff check src tests`; `uv run --frozen --no-sync pytest -q` | 1198 passed, 1 skipped (torch not installed) | local machine, started at 1435c22 on a clean tree (`runs/logs/checks_local_2026-11-13b.log`); the tests added later ran on their own; `./init.sh` itself not run |
| Reproducibility of R-A's candidate files | `experiments/builds/build_test.py embryo ...` at 16 BLAS threads, then `vec check` | one sha256 per file over four and three builds | `experiments/builds/results/determinism_1f8956b.txt`; build logs not read |
| DGX rehearsal of the amended entries (2026-11-13) | two chains of `BUILD_OUT=runs/rehearsal_fd2683a/<a\|b> bash experiments/builds/build_all.sh s1`, `s2`, `t3alt`; then `bash experiments/builds/determinism_check.sh runs/archive/builds_fd2683a/SHA256SUMS 16` | ten of ten files byte-identical across the two chains and the third build; seven equal to the archive of f337dc6 | `experiments/builds/results/rehearsal_fd2683a.txt`; exit codes, hashes and the stamps' commit and clean fields only |
| Dry builds of the amended entries | the commands of `docs/S2_RULES.md` (last section) at 16 BLAS threads, then `vec check` | one sha256 over three builds of each new file; unchanged files reproduce the DGX archive, T1's old s1 excepted | `experiments/builds/results/rulings_9cf5ca5.txt`; build logs not read |
| T3 tries K7 and K8 | `experiments/t3_decomp/k7_wt_revert.py` and `k8_s2_form.py` (`run` per task, then `read`) | `experiments/t3_decomp/results/k7_wt_revert/`, `k8_s2_form/` | 10 seeds, clean trees at 33204be and 9dc6af1 |
| Lint and tests | `./init.sh` | passes (2026-10-28) | `uv run pytest -q` 402 tests (90+ on 2026-10-19); experiment tests in `experiments/*/test_*.py` run separately |
| Pre-registered runs (step 4) | `docs/P3_RUNBOOK.md` step 4 commands (DGX, 2026-10-28) | all six ran, exit 0 | outputs committed before reading: f5be247 to 12151e0; the commits are listed in the runbook's step 4 status |
| Size procedure | `uv run pytest experiments/coords -q -p no:cacheprovider` | 33 passed (2026-10-19) | preflight requires committed procedure files |
| QC of released files | `bash experiments/qc/run_qc.sh <tag>` (DGX) | `experiments/qc/results/` | duplicate and interpolated counts in `docs/DATA_NOTES.md` |
| Validation answers (step 2) | `uv run vec data inventory`, `uv run vec data qc` (DGX, 2026-10-28) | 6 files, 2.09 GB, sizes as listed | SHA-256 of all 15 files in `data/raw/SHA256SUMS`; no duplicated or interpolated records; `docs/DATA_NOTES.md` (f423f6a) |
| Labels (step 3) | `uv run python experiments/labels/release_label_evidence.py` (DGX) | label evidence in `experiments/labels/results/` | every label of the 15 files maps (DGX check, 2026-10-28); the new entries pinned in `tests/test_prep.py` and `tests/test_t1_scope.py` (d6b9566) |

## Decisions Made

- Current: `docs/PLAN.md` decision log through 2026-11-13 (f), `docs/DATA_POLICY.md` change notes through
  2026-11-13 (e), the second-upload rules in `docs/S2_RULES.md` as amended on 2026-11-13, and the captain's upload gate in
  `results/upload_gate.json`. The campaign's method decisions are the AI assistant's, logged; every change to an
  upload is the captain's ruling.
- See `docs/PLAN.md` decision log (2026-10-18, 2026-10-19 a-d, 2026-10-20 a-c, 2026-10-21 a-c, 2026-10-28) and
  `docs/DATA_POLICY.md` records. The user's calls of 2026-10-27 and 2026-10-28 (the history rewrite, the unlock, the
  A20 reading, the go-ahead for step 3) are recorded in `progress.md` (Decisions Made, Blockers / Risks).

## Blockers / Risks

- 2026-11-13 (late night, after decision log (gg)): request R-D rests on nine practice tasks of which three share
  one target stage and two are the counted backtests that prompted the read, and reading N39 shows its gain does
  not carry to the test's geometry; the AI assistant advises against it. If it were ruled in all the same, the new
  first upload is a rendered file, which reproduces only at the build's 16 BLAS threads, and at the test one
  sample's file lacks two panel genes (they stay as `fgw` has them, a rule fixed before any build).
- 2026-11-13 (night, after the sixth round): in our own evidence a total past the leading 215.41 lies beyond the
  best case of every board at once (their best practice outcomes sum to 219.0), and the forecast's upper reading
  (194.4) leans on one judgement that the record itself rates lower: that the beta-catenin samples share the Gata4
  slide's profile (weight one half; at 0.15 the mean is 189.0). The committed 184 is the figure to quote.
- 2026-11-13: the estimate (184) is about 25 under the gate, and no change was adopted in the campaign. What can
  still move the outcome is how the hidden samples were measured (the bets behind R-A, R-B and T3's first-choice
  file), which no released file shows.
- 2026-11-13: no organiser answer is recorded to any of the seven emails, the Slack DM or the public Slack
  question. The next screened check is on 26-27 Nov; one follow-up in the Slack thread may be drafted about 19 Nov.
- 2026-11-13: a ruling that changes a build is needed by about 24 Nov (a T3 builder change takes about a day, then
  a dry build and the rebuild rehearsal). A rendered embryo file reproduces only at the build's 16 BLAS threads.
- 2026-11-13: only DGX hashes check DGX builds. The three T3 files that apply a per-gene factor have other bytes on
  the local machine (as T1's s1 had until that day), so a change of the DGX's numerical libraries or processor before 27 Nov would show
  as a mismatch on the freeze morning; a file that differs is not uploaded until the captain has the explanation.
- 2026-11-13: the DGX's default folders `runs/builds` and `runs/t3_build` hold older files with clean stamps (T1's
  two date from 2026-11-04). Uploads are pulled only from the folder the freeze-morning command builds
  (`docs/P3_RUNBOOK.md`, "The freeze-morning rebuild").
- 2026-11-13: the account's usage limit stopped an agent workflow that morning; single reviewer agents ran later
  the same day. Larger agent workflows wait for the captain's word.
- Earlier items, kept as history:
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
2. Read `feature_list.json` and `progress.md` (Current State, What's Next).
3. Review this handoff, `docs/P3_RUNBOOK.md`, and `docs/PLAN.md` decision log 2026-11-13 (a) to (f).
4. Run `./init.sh` before editing. On the DGX do not run it as it is (`docs/COMPUTE.md`, `progress.md` notes).
5. Before choosing or changing work on an entry, read the guiding principle in `AGENTS.md` and the headroom tables in
   `experiments/harness/TRIES.md`.

## Recommended Next Step

- 2026-11-14 (latest, after decision log (ii)): the round's audit is done. Every line that independent passes over the
  record found open is read, and no closure rests on the AI assistant's own odds; nothing passes and no file changes.
  Three things are the captain's, all optional: a ruling on request R-D or on its runner-up `fgw_aligned_trend` (the
  AI assistant advises no on both, and says the advice is an after-the-fact lens), whether to send the drafted eighth
  email on final marks, and whether the Task 3 batch question goes into the follow-up about 19 Nov. No run is queued,
  and the AI assistant knows of none whose result could change an upload under the standing rules. Lessons of the
  round, for the next brief: a check that reads practice scores gets its row first; a gate is tried on pure noise
  before it is trusted; a control is checked for whether it can differ from its arm; a brief names where a formula's
  translation goes; and a runner records an arm the builder refuses instead of stopping. The dated items follow as
  listed below. Report the expected total as 184 (170 to 198), with 194.4 named as the forecast's upper reading.
- 2026-11-13 (latest, after decision log (gg)): the sixth round is closed and the options for passing the leading
  total are used up; every allowed line with a counted read has been read. **No ruling is needed for the entry
  as built; one optional ruling is on record for the captain.** Request R-D (the embryo's first upload to
  `fgw_trend_render`) is on record as a pair rule's outcome
  with the AI assistant's advice against it (decision log (gg), reading N39); it changes nothing unless the captain
  rules for it, and then: port the branch `port/emb-s1-trend` to main after its review, set the embryo's
  first-upload expression in `experiments/builds/build_all.sh`, dry-build, rehearse a fourth time on the DGX, write
  the new reference list and the pointer `results/reference_list.json`, ask the captain before the one server
  format check of the new file, and update `METHOD_SUMMARY.md`. Before recommending any request built on practice
  tasks: name the property of a task the effect depends on, compute it for the test from target-free quantities,
  and run the cold check first (decision log (gg)). The dated items follow: the screened reads on 19 and 23 Nov,
  the follow-up draft about 19 Nov, the method summary and disclosure email on 24-30 Nov, and the freeze morning.
  Report the expected total as 184 (170 to 198), with 194.4 named as the forecast's upper reading.
- 2026-11-13 (latest, after decision log (ee)): the sixth campaign round is closed with nothing kept and nothing
  requested. No run is queued on any GPU, and the AI assistant knows of none whose result could change an upload
  under the standing rules. Next are the dated items: the screened reads on 19 and 23 Nov, the follow-up draft
  about 19 Nov (floored severity, and the Task 3 batch question if the captain wants it back), the method summary
  and disclosure email on 24-30 Nov, and the freeze morning. Report the expected total as 184 (170 to 198), with
  194.4 named as the forecast's upper reading. Where a counted backtest exists and a read costs minutes, write the
  row and run it rather than pruning on pre-run odds (`AGENTS.md` rule 8).
- 2026-11-13 (latest, after decision log (cc)): the heart ruling is built in and rehearsed (reference
  `runs/archive/builds_dc7b6ae/SHA256SUMS`; every file format-checked on the server; the upload client refuses a
  scored upload of a file off that list). Next: the dated items below, reading `builds_dc7b6ae` for the reference
  and skipping the server format checks, which are done. No run is queued on any GPU. On the upload mornings the
  12 characters to confirm come from the list: heart 3d6548382b24 (first upload) and 18c2cf552f98 (second).
- 2026-11-13 (after decision log (m)): the captain ruled, T1's second upload is on the mid base and the
  ten files are rehearsed again (reference `runs/archive/builds_38616fe/SHA256SUMS`, since superseded). The
  generator tests the captain ordered are done and changed nothing (decision log (r)).
- 2026-11-13 (earlier that day): the rulings are given and built in, and the ten files are rehearsed on the DGX. Next, in
  order: the screened reads on 19 and 23 Nov and the Slack follow-up draft about 19 Nov; the method summary and the
  disclosure email (with the hashes of `runs/archive/builds_38616fe/SHA256SUMS`, the list since the second
  rehearsal) on 24-30 Nov; the screened check on
  26-27 Nov; on 27 Nov the DGX sync to the freeze commit, the freeze-morning command (`determinism_check.sh
  runs/archive/builds_38616fe/SHA256SUMS 16`, the list since the second rehearsal), the server format checks and
  the first round, each upload on the
  captain's go-ahead and pulled from the folder that command built; the second round on 29-30 Nov; the final marks on
  1 Dec. Before any reading that compares a predictor with a fold target, write its ledger row (row N25 records the
  one time this was skipped).

Earlier notes, kept as history (superseded where `progress.md` says otherwise):

- 2026-11-04: follow `progress.md`, What's Next: the T1 births decision (by 2026-11-10 at the latest), the
  T3 response-identity line through the adoption rule, and the headroom tables kept current. Official uploads only on
  the captain's go-ahead. Every change to an entry starts from its board's published baseline and a row in
  `experiments/harness/TRIES.md`, registered before its run (rules 5-8).

- 2026-11-01: the s1 build path was reviewed (four reviewers, skeptic-verified; 11 findings) and fixed in
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
