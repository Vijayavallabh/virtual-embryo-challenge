# Session Handoff

## Current Objective

- **State at 2026-11-20, 19:23 UTC (newest; it supersedes the bullets below where they differ).** Decision log
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
