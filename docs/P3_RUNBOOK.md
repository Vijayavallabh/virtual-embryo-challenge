# Release-day runbook (2026-10-27, start of the final phase)

The organisers moved the start from 2026-10-20 to 2026-10-27 (email of 2026-10-20). Files written before that email say
2026-10-20; they mean this release.

The validation answers (T1 E10.5, T2 heart E8.5 and E10.5, T2 embryo E7.5, T3 Gata4 replicates) are served only after
a team member unlocks them on the site's Data page, which permanently ends P2 submissions for the whole team (our
unlock: 2026-10-28 09:32 UTC). No test inputs are distributed. Everything below that decides a method was committed
before the release; run it as written and commit each output before looking at anything else. Steps marked **(user)**
need you.

Progress (2026-10-28): step 0 passed; step 1 is blocked (its saved-file fallback is built, cee1c92, and waits for the
user's call); steps 2 and 3 are done; step 4 ran on 2026-10-28 at the user's go-ahead. Each step's status is noted
under its heading.

## 0. Before anything

Status: passed on 2026-10-27 (phase p3, split test, accepting submissions; official quota keys T1:test,
T2:heart:test_extrap, T2:embryo:test_interp and T3:bcatenin, 2 scored submissions each for the whole phase, none used;
nothing unresolved).

- **(user)** Do not open `panels/index.json` on the challenge site or the public board endpoints, and skip Slack threads
  that discuss sizes or other values worked out from anchors. Test-board anchors, if published, encode the hidden
  targets (`docs/DATA_POLICY.md`, do-not-use list).
- `git pull`, `./init.sh` (lock, lint, tests).
- `uv run vec submit status` (local, read-only; uploads nothing): the phase, its split and, per board, the quota and
  the ledger and server view. If it still shows P2 (split `val`), the release has not happened: stop here and run
  nothing below. If the official quota (`official_per_board`) shows as unknown, the client keeps refusing
  official uploads until its parser is adapted to the real shape and reviewed; nothing else depends on it today.
- An `unresolved` line there means a format check came back as a submission (or unconfirmed): nothing can be sent to
  any board until it is resolved. Check the portal's submission list, then run the printed command for each board it
  names, e.g. `uv run vec submit resolve SHA12 --board KEY --as acknowledged --note "portal shows submission N"` (a
  stray is recorded under the board the server named and the board the command named, so expect up to two). An
  `ambiguous` line is an official upload with an unknown outcome; it counts until the server lists it, or until
  `uv run vec submit resolve SHA12 --board KEY --as not-stored --note "..."` after the portal shows no such upload.

## 1. Board specs (anchors stripped in memory)

Status (2026-11-04): done. `vec boards update` still answered 403 at 07:29 UTC; the screened Slack and email sweep
(wf_2489e0c9-951, two clean audits) found no objection, no fix and no rule change since the notice; the captain
saved the board index and the embryo test gene list, and `vec boards update --from-file` read them (634b5fd; the
saved index deleted, `~/vec-panels` removed). D1 (ae6d6b3): every test board publishes its cell limits, gene list
and keys, but no truth_cells or ref_cells, so A3 is unresolved. Exposures on the way (the file and the Data page
pasted into the AI assistant's session): `docs/DATA_POLICY.md`, 2026-11-04 and 2026-11-04 (b).

Earlier status (2026-10-31): blocked. Since 2026-10-27 `vec boards update` fails: `panels/index.json` answers HTTP 403
(Cloudflare) to our scripted client (still so on 2026-10-31 at 16:52 UTC). It is not worked around (no other user
agents or endpoints), so `boards/` still holds the P2 specs. The fourth organiser email (sent 2026-10-28; no reply by
2026-10-31) asks whether the block is intended and how to get the test boards' specs, and says we have not tried to
get around it; the test-board intake (section 5) waits for them. So does every upload: `vec submit` checks the board
against `boards/index.json` first (the board guard in `vec.submit.prepare`), so until the test boards are there it
refuses every upload to them, format checks included. This is the main blocker for submitting; the pre-registered
runs of step 4 do not need it. The saved-file intake below was built on 2026-10-28 at the user's request (cee1c92,
after an independent review). 2026-10-31: the 403 is Cloudflare error code 1010, a browser-signature block (the
message now names it), and the predicted organiser answer is no objection (`docs/ORGANIZER_QUESTIONS.md`, working
answers). 2026-11-01: the captain cleared the route ("yes for all"), with these conditions: he sends the dated notice
drafted in his Gmail first (a reply on the fourth-email thread; he sent it on 2026-10-31 at 18:59 UTC); the save
happens not before 2026-11-04 (UTC), and not if an objection or a fix arrives first. After the intake, run `uv run python experiments/intake/d1_check.py` (D1 and
A3; it never prints test-board counts).

```bash
uv run vec boards update      # same as vec.boards.update_boards(); test anchors stripped in memory, never printed
git add boards && git commit -m "Board specs for the final phase (test anchors stripped at intake)"
```

Fallback, **cleared by the captain on 2026-11-01** for this one route (`docs/DATA_POLICY.md`). Before step 1: the
notice is sent, it is 2026-11-04 (UTC) or later, no objection or fix has arrived (a plain `uv run vec boards update`
still gets the 403), and hovering over the Data page's board-index link (no click) shows
`https://virtualembryo.ai/challenge/panels/index.json`, the URL the tool uses. Create the folder with `mkdir -m 700
~/vec-panels`. On any refusal, delete the folder at once (`rm -r -- ~/vec-panels`) rather than keeping the file for a
re-run, and save it again only if the refusal names a fixable cause. The steps:

1. **(user)** On the site's Data page, right-click the link to the board index (`panels/index.json`) and choose "Save
   link as" into a new, empty folder outside every git work tree: `~/vec-panels/` on this machine (`~` is
   /home/sports), not next to the repository, whose parent folder is a git work tree too. If the browser runs on
   another computer, copy the file here without opening it (e.g. with `scp`) and delete the browser's copy. Never open
   the file, and never type its address into the browser (that displays it); if the page has no such link, stop and
   tell the AI assistant.
   The Data page itself lists every board's cell limits, says how the upper limit is set and gives hidden targets'
   file sizes: never copy its text to anyone, the AI assistant included (`docs/DATA_POLICY.md`, 2026-11-04 (b)).
2. Run the intake on it. Nobody opens the saved file; the intake reads it in memory and strips the test anchors as
   the online refresh does. The commit runs only if the intake succeeded, and stages only the files it writes:

   ```bash
   uv run vec boards update --from-file ~/vec-panels/index.json --delete-source && \
     git add boards/index.json boards/*.genes.txt boards/t1_composition.json && \
     git commit -m "Board specs for the final phase (saved index; test anchors stripped at intake)"
   rm -r -- ~/vec-panels   # afterwards, unless the intake kept the files for a re-run
   ```

   Each gene list comes from that folder or from a file already in `boards/` with the sha256 the index publishes. If
   the intake names a gene list it still needs, save that file from the Data page into the same folder and run it
   again (nothing is written until every list is found). A `t1_composition.json` saved there goes through the same
   released-stage rule; without one the old file stays. `--delete-source` deletes the saved index (and a saved
   composition file) after a successful intake; without it the intake prints the `rm` line. A refusal writes
   nothing and keeps the saved files: do not open them to find out why (the message says what was wrong). If the
   summary says the phase endpoint could not be read, the files are kept and the intake prints the command to run
   again. It also names other files in the folder that look like copies of the index (e.g. `index (1).json`).

The summary lists the boards, which ones were reduced to their submission spec, and whether the T1 composition file was
written (it is skipped if it describes any unreleased stage). Never print `boards/index.json` anchors for a reduced
board; they are not there. E10.5 counts as unreleased until the public phase endpoint reports the test split; the
summary says why the file was held back. If step 0 showed the final phase and the summary says the endpoint could not be
read, re-run `uv run vec boards update` before committing.

## 2. Data

Status: done on 2026-10-28 (f423f6a). After the user's unlock (09:32 UTC) the data client fetched the validation
answers (support added in e0e1d58): the six files (T1/E10.5_RNA, T2/E8.5, T2/E10.5, T2/E7.5, T3/E8.75_gata4_ko_1 and
_2; 2.09 GB) reached `data/raw/T*/` on the DGX at 09:51-09:53 UTC, sizes as listed. The site publishes no checksums,
so the SHA-256 sums of all 15 files are in `data/raw/SHA256SUMS` there; inventory and QC are in `docs/DATA_NOTES.md`.

```bash
uv run vec data manifest      # local: needs the team token; training files, the unlock state, the validation answers
uv run vec data fetch         # local: streams each signed link to curl on PrakashDGX_H2 (default --dest data/raw there)
ssh PrakashDGX_H2 'cd ~/v/virtual-embryo-challenge && export PATH=$HOME/.local/bin:$PATH && \
  uv run vec data inventory && uv run vec data qc data/raw/T*/*.h5ad --csv runs/qc.csv'
```

Check sizes and checksums against the manifest, record the new files and the QC table (z levels, steps, duplicated
records, row order) in `docs/DATA_NOTES.md`, and commit.

## 3. Labels (before any method touches a new file)

Status: done on 2026-10-28 (d6b9566), at the user's go-ahead. The AI assistant mapped the 30 labels first seen in the
validation answers in `configs/lineage_map.yaml` by its 2026-10-19 rule (three independent mappers, 29 of 30
unanimous, a judge for BW) and the 8 labels new to a T1 file in `configs/t1_chamber_classes.yaml`. Evidence:
`experiments/labels/` (script aebcb8c) and `experiments/t1_scope/results/label_evidence_E9.5_E10.5/`. New core labels:
Endocardium-1, Endocardium-2 and EndoMT Endocardium (T1), Valve-EndoMT-Endo (T2 E10.5) and the Gata4 files' CM. Every
label of all 15 files maps (DGX check).

Map every new label in `configs/lineage_map.yaml` with the rule in its header (core = cardiomyocytes and endocardium;
ambiguous labels default to shell after a marker check with `scripts/label_markers.py`). Commit the map. The size
procedure refuses unmapped labels.

## 4. Pre-registered runs, in order (commit each output with the git hash)

Status (2026-10-28): run at the user's go-ahead ("build the intake and start step 4"), in order, every run on the DGX
at 725d19c with a clean tree, each output committed before it was read. 4.1 f5be247 (both size JSONs valid); 4.2
95d8198 (replicate overlap: no shared rows, independent), 7e7c207 (decision file), 3ff6bf0 (aggregates), 05ead40
(plan); 4.3 3c5cc66 (verdict: build); 4.4 0586bf6; 4.5 dfee546 (decision: max_cells kept on both boards) and
ef1819f; 4.6 12151e0. On the DGX, a command here that no runner script wraps needs the BLAS thread caps in
`docs/COMPUTE.md`.

1. **Heart and embryo size:** `bash experiments/coords/run_p3.sh` (paths as arguments). Commit
   `experiments/coords/results/p3/`. The heart JSON sizes both heart submissions through `vec.t2_build.size_or_hold`
   (approved rule 3, `docs/RELEASE_RULES_SUMMARY.md`): as recorded when `valid` is true; when `valid` is false while
   E10.5 is usable, both keep E10.5's size. If E10.5 is barred (A20), the size rule goes back to the user
   (2026-10-28: A20 passed on the user's reading, so E10.5 is usable; organiser confirmation is asked in the fourth
   email). An embryo JSON with `valid` false is not used.
2. **T3 decomposition on Gata4.** First count identical count rows shared by the two replicate files (rules.yaml
   note 4: more than 1% shared means the borrowed-replicate rules are recorded as not independent):
   `uv run python experiments/t3_decomp/replicate_overlap.py REP1.h5ad REP2.h5ad --out
   runs/t3_decomp/replicate_overlap.json`; bring the JSON back and commit it with the decision file. Then:
   `tmux new -d -s vec-t3-gata4 'bash experiments/t3_decomp/run_gata4.sh REP1.h5ad REP2.h5ad data/raw/T2/E8.75.h5ad
   --e85 <released E8.5 file>'`. Check that `runs/logs/t3_decomp_gata4_*.log` exists (a wrong path exits 0 silently).
   Bring back the aggregates with the rsync line the run prints and commit `gata4_decisions.json` first. The rules in
   `experiments/t3_decomp/rules.yaml` decide, including the dated notes (approved by the user on 2026-10-19); nothing
   is chosen by eye. The decision file shows the tables before notes 6-8; `uv run python
   experiments/t3_decomp/build_submission.py plan <gata4_decisions.json> --out <same folder>/gata4_plan.json` prints
   what `build` will make for every branch with the notes applied. Commit the plan with the decision file and show
   it to the user.
3. **T1 direction readout:** `uv run python experiments/t1_direction/readout.py <released E10.5 file>` (DGX). It
   refuses uncommitted procedure files, a non-default config or label map that is not a committed repo file, and a
   second run without `--overwrite`. Commit `decision.json` (build the bracketing pipeline or not). Re-run the
   committed filters (`experiments/external/`) before ~/v/vec-gapconc is cleaned up.
4. **Harmonisation:** `experiments/harmonise/run_harmonise.py` release stages (Gata4 replicate 1 vs 2, cross-stage
   transfer with E8.5 and E10.5). The replicate difference is the noise bar for T3 composition effects (A19). The
   script's default Gata4 paths (`..._rep1.h5ad`, `..._rep2.h5ad`) are not the released names, so the run passed
   `--stages release --out runs/harmonise_release --gata4-rep1 data/raw/T3/E8.75_gata4_ko_1.h5ad --gata4-rep2
   data/raw/T3/E8.75_gata4_ko_2.h5ad`; tables in `experiments/harmonise/results/release/`.
5. **Density:** the rule code (`src/vec/t2_build.py`, `experiments/density/summarise.py`) must be committed and synced
   to the DGX first (the rsync line in `docs/COMPUTE.md`): the decision refuses uncommitted rule code. Then, on the
   DGX from `~/v/virtual-embryo-challenge` (the rerun takes about 20 minutes; use tmux if the session may drop):

   ```bash
   TARGETS="heart_E10.5 embryo_E7.5" SEEDS="0 1 2 3 4 5 6 7 8 9" RAW=runs/density/raw_release \
     TAG=density_release JOBS=16 THREADS=1 bash experiments/density/run_all.sh
   grep -l 'skipped:' runs/logs/density_release_*_s*.log  # must print nothing: the decision is written once
   uv run python experiments/density/summarise.py --raw runs/density/raw_release \
     --out experiments/density/results/release           # tables (CSV) and plots
   uv run python experiments/density/summarise.py --raw runs/density/raw_release \
     --out experiments/density/results/release --decide  # writes count_decision.json there
   ```

   A log the grep names found no release file (`run_all.sh` passes no `--target-file`). Rerun that target, then
   decide (finished conditions are skipped):

   ```bash
   uv run python experiments/density/run_density.py --target heart_E10.5 --seeds 0 1 2 3 4 5 6 7 8 9 \
     --out runs/density/raw_release --target-file <the released E10.5 file>  # or embryo_E7.5 and its file
   ```

   Decide mode applies approved rule 7 (`docs/PLAN.md` section 3.6) in code: a lower level replaces max_cells on a
   board only if, for the intermediate-quality condition, its paired mean gain over max_cells exceeds 0.3 points and
   it beats max_cells in at least 8 of the 10 seeds; the largest gain wins (ties go to the higher count); heart and
   embryo are decided separately. A missing seed or a level capped by the pool keeps max_cells. Use exactly that
   `--out`. A second decision needs `--overwrite`, which replaces `count_decision.json` (the file the T2 build reads)
   and keeps the first beside it as `count_decision.replaced-<UTC time>.json`. Then, from the local machine, bring the
   folder back and commit the decision before building any T2 file (untracked result files do not block an upload,
   so committing is this step's job):

   ```bash
   rsync -az -e ssh PrakashDGX_H2:v/virtual-embryo-challenge/experiments/density/results/release/ \
     experiments/density/results/release/
   git add experiments/density/results/release/count_decision.json
   git commit -m "Density count decision (approved rule 7) from the release rerun"
   ```
6. **Scorer build:** `uv run python experiments/scorer_build/anchor_check.py check` against the released validation
   answers (validation anchors only; it refuses test boards, the public board file and anything before 2026-10-27).
   Its verdict says which build the server baked the anchors with; it changes no prediction.

## 5. Then

- Emulator calibration against the published validation reference rows (feat-005, assumption A2). Done 2026-10-31
  (`experiments/calibration/README.md`): every published validation anchor and four of the five organiser baseline
  rows reproduce exactly under `protocol="bake"` and the prefix build; the heart-extrapolation pseudobulk_shift row
  does not (cell-type matching across E8.75 and E9.5 is unpublished). On the validation boards, a local score at
  seed 0 with the published anchors is the server's score.
- Test-board intake: reference stage, cell limits and what else each test board publishes (decision point D1).
- Pseudo-task harness (feat-008, `vec.harness`, `configs/pseudo_tasks.yaml`): every model change is compared with
  its base on its test board's pseudo-tasks and kept only under the adoption rule (`docs/PLAN.md` section 3, item 8):
  `uv run python experiments/harness/run.py run --board KEY --candidates floor,copy,... --out runs/harness/<id>`
  (DGX), then `... decide --runs runs/harness/<id> --board KEY --base floor --change NAME --out
  experiments/harness/results/<id> --register`. `decide` refuses results from different commits, a dirty tree or
  another config.
- **Build seed (2026-11-01): never build a test file at seed 0.** numpy's `default_rng([s, 0])`, which
  `vec.t2_build.sample_cells` draws rows with, is the same stream as `default_rng(s)`, and the server draws its
  reference working copy with seed 0: a T2 file built at seed 0 contains the server's own floor cells. Use the
  harness's offset (`vec.harness.BUILD_SEED_OFFSET`, 1,000,000) or any seed above it, and record it with the file.
- Candidate files, built on the DGX with the committed builders and checked there with `uv run vec check FILE --board
  KEY` (no token, nothing sent). A server format check runs from this machine with `uv run vec submit runs/<...>.h5ad
  --board <test board key> --from-dgx` (the file is pulled for that call only; nothing official without `--official`;
  its output goes to a file, never through `head`, which can cut the client off mid-call). Nothing is uploaded
  officially today. Every built file gets a build stamp (`<file>.stamp.json`) that an official upload checks: commit
  before building, and do not sync, pull or check out on the DGX while a build runs there (the stamp then is not clean;
  rebuild after the sync).
  - T3: `uv run python experiments/t3_decomp/build_submission.py build <gata4_decisions.json> --branch <b>
    --submission s1|s2 --wt data/raw/T2/E8.75.h5ad --run-dir runs/t3_decomp/gata4_<STAMP> --board <T3 board key>
    --seed 1000001 --out runs/t3_build/<b>_<s>.h5ad`, for `allowed`, `refused` and `barred` (s2 = `refused` s1; `docs/S2_RULES.md`). The cell count defaults to the
    board's max_cells, read in code and never typed or printed (s1 build review, 2026-11-01); the builder refuses a
    seed below `vec.harness.BUILD_SEED_OFFSET`. 2026-10-31: the
    allowed branch is selected on the predicted answer to follow-up question 3, at the captain's instruction
    (`docs/ORGANIZER_QUESTIONS.md`); the barred branch stays built as the switch if a ruling arrives before the T3
    uploads. Approved by the captain on 2026-11-01: s1 = `--branch allowed --submission s1`; s2 = `--branch refused
    --submission s1` (wild-type cells, the central composition at the Mab21l2 magnitude, no number derived from the
    Gata4 files). Upload s2 first and s1 last (by 30 Nov IST; since 2026-11-05 not before the captain's upload gate
    is crossed or its deadline exception opens on 27 Nov, `docs/PLAN.md` section 8), each with the captain's go-ahead
    after a screened check
    of mail, Slack and the Rules page. If the organisers allow the Gata4 files before the first T3 upload, s1 stays
    the allowed branch's s1 and s2 follows the 2026-11-06 switch below (decision log 2026-11-06 (b), Interactions: the
    switch takes the s2 slot even when the Gata4 files are allowed); if they restrict them before the last one, send nothing from the allowed branch
    (the barred branch's s1; `docs/S2_RULES.md`, which also covers a full bar). It re-decides nothing and refuses a decision file whose script or rules hash
    differs from the files on disk.
    2026-11-06, the captain's ruling ("go with your recommendation"; `docs/PLAN.md` decision log 2026-11-06 (b)):
    - The screened check before T3's first upload also looks for the organisers' answer to the fifth email, on how a
      floored severity scores (`docs/ORGANIZER_QUESTIONS.md`), and from 2026-11-11 in the DM and in the replies to the
      public Slack question once the captain has posted it. Only an organiser's answer counts; an answer other than
      zero or 0.5 (for example "dropped") goes to the captain before T3's s2.
    - If they have said zero, T3's s2 is `runs/t3_build/barred_s1.h5ad`, even if the Gata4 files are allowed.
    - If they have said half, or nothing by then, T3's s2 stays the refused s1.
    - If the barred s1 went first and the Gata4 files are restricted before the last upload, the last upload is
      `runs/t3_build/barred_s2.h5ad`, and the higher published score of the two Gata4-free uploads is marked.
    - An answer that arrives after T3's s2 is up changes nothing. The disclosure ids are unchanged (X001, X004, X005).
    2026-11-10 (`docs/PLAN.md` decision log 2026-11-10 (s)): the screened check before the first upload (26-27 Nov) also
    reads, from the Rules and Prizes pages, where and by when the method summary goes in, and how the source and stages
    of external data are to be stated with the submission: the client keeps the `--disclose` ids in the local ledger
    and sends none of them. It also looks for the organisers' answer to the sixth email (`ASSUMPTIONS.md` A21), and quotes the Rules' definition of a
    Human Team and anything they say about AI assistance (decision log 2026-11-10 (x)).
    2026-11-13 (`docs/PLAN.md` decision log 2026-11-13 (c)): the seventh organiser email, sent that day, asks which
    released samples the hidden E7.75 sections and E12.5 Multiome libraries were processed with. An answer may come by
    reply, in the Slack channel or as a line on the Data Protocol page, and each is read through a screen with the scope
    fixed in `docs/ORGANIZER_QUESTIONS.md` ("Seventh email"): procedure only, every number replaced, everything else
    withheld. Two extra screened reads are asked of the captain for about 19 and 23 Nov, because a build change has to
    be rehearsed by 24 Nov; the check of 26-27 Nov repeats the read. The captain approved both on 2026-11-13; the
    screening agent's brief is `docs/SCREENED_READ_SEVENTH_EMAIL.md`. A reply in the captain's inbox is not pasted into
    the session. What each answer changes is fixed in the same place, and nothing changes without the captain's ruling.
  - T2: `vec.t2_build`, run from the repo root (the paths are relative). Count by `n_cells_from_decision(board,
    pool.n_obs, "experiments/density/results/release/count_decision.json")` (step 4.5; null keeps max_cells; a
    decided level below the test board's lower limit gives way to the next passing count up, the passing level with
    the fewest cells the board allows, else max_cells, and the count record's `count_floor` says so; the user approved
    the words on 2026-10-21 and left the reading of "next" to the AI assistant, whose decision panel chose the count
    up: `experiments/density/floor_reading/README.md`). If step 4.5 wrote no decision because the rerun left no rule
    rows, approved rule 7 keeps max_cells: count with `n_cells(board, pool.n_obs)` and
    record why in `progress.md`. Then `sample_cells`. Heart size by `size_or_hold(coords,
    "experiments/coords/results/p3/size_heart_E12.5.json", "H1" (s1) or "H2" (s2; the heart s2 is H2 built as the half
    trim, `build_test.py heart --entry s2_trim` and, since the ruling of 2026-11-13, `--trend` on both heart
    entries; trim before the draw, `docs/S2_RULES.md`), source_stage="E10.5",
    reference_cloud=<the full E10.5 pool from load_stage(..., drop_duplicates=True), before any scope filter>,
    expect_board="heart")`; embryo size by `apply_size_procedure` with
    `experiments/coords/results/p3/size_embryo_E7.75.json`. Then `cloud_guard` and `build_submission`. The module
    docstring of `src/vec/t2_build.py` shows these calls in order; `experiments/t2_build/smoke.py` counts with plain
    `n_cells` and sizes with `rescale_to_rms`, so copy only its `sample_cells`, `cloud_guard` and `build_submission`
    calls.
  - T1: `vec.t1_scope.build_submission` once the shares are fixed by their written rule (docs/PLAN.md section 4.1);
    pass `ref_source` (the reference file's stem) so the per-file label remaps apply.
  - **The amendments of 2026-11-13** (`docs/S2_RULES.md`, last section; `docs/PLAN.md` decision log 2026-11-13 (d)).
    `build_all.sh` now builds: T1 s1 with `--base lib3` (library 3's cells; s2 unchanged on the pooled base); the
    embryo s2 as `build_test.py embryo --method fgw --move-coupling aligned --expression series`; T3's refused s1
    with `--wt-revert` (it reads the six other wild-type stage files under `data/raw/T2`). The other seven files are
    as archived. Dry builds on the local machine (`experiments/builds/results/rulings_9cf5ca5.txt`): each new file
    twice with one sha256, `vec check` passed, and the unchanged T1 s2, embryo s1 and T3 refused s1 without the
    option reproduce `runs/archive/builds_f337dc6/SHA256SUMS`.
  - **The amendment of 2026-11-13, later that day** (`docs/S2_RULES.md`, last section; row R7; `docs/PLAN.md`
    decision log 2026-11-13 (k) to (m)). By the captain's ruling T1's s2 is built as `build_test.py t1 --g-mom --base
    mid --cm-share 0.60` (the base half-way between the pooled base and library 3); `build_all.sh s2` builds it. The
    other nine files are unchanged. Everything below that names the reference was remade the same evening: **the
    reference is now `runs/archive/builds_38616fe/SHA256SUMS` on the DGX** (in the repository:
    `experiments/builds/results/builds_38616fe.SHA256SUMS`; record `experiments/builds/results/rehearsal_38616fe.txt`),
    and the freeze-morning command is `bash experiments/builds/determinism_check.sh
    runs/archive/builds_38616fe/SHA256SUMS 16`. Where the three paragraphs below say `builds_fd2683a`, read
    `builds_38616fe`; they are kept as the record of the first rehearsal.
  - **The amendment of 2026-11-13, evening** (`docs/S2_RULES.md`, last section; row H11; `docs/PLAN.md` decision
    log 2026-11-13 (w) to (bb)). By the captain's ruling both heart files are built with `--trend`
    (`build_test.py heart --entry s1 --trend` and `--entry s2_trim --trend`: each lineage's step is the
    least-squares trend over E8.25_late to E10.5 instead of the E9.5 to E10.5 step alone); `build_all.sh` builds
    them so. The other eight files are unchanged. The reference was remade by a third rehearsal: **the reference
    is now `runs/archive/builds_dc7b6ae/SHA256SUMS` on the DGX** (in the repository:
    `experiments/builds/results/builds_dc7b6ae.SHA256SUMS`; record
    `experiments/builds/results/rehearsal_dc7b6ae.txt`), and the freeze-morning command is `bash
    experiments/builds/determinism_check.sh runs/archive/builds_dc7b6ae/SHA256SUMS 16` (tried: exit 0, ten of ten
    OK). Wherever this section says `builds_38616fe` or `builds_fd2683a`, read `builds_dc7b6ae`. **Every file of
    the set has passed the server's format check** (six on 2026-11-13 at the captain's word, decision log (x); four
    on 2026-11-04), so the upload mornings need none unless a build changes. The earlier archive folders, on the
    DGX and on the local machine, still hold the two heart files this ruling replaced, with clean stamps: for the
    heart the 12 characters to confirm are 3d6548382b24 (first upload) and 18c2cf552f98 (second), never
    e73d03963648 or 5a309434d4a2. The checked local copy is `runs/archive/builds_dc7b6ae/` (ten of ten OK).
  - **The client checks the list (2026-11-13).** `vec submit upload ... --official` now carries the guard
    "reference list": it refuses unless the file's sha256 is in the list the pointer
    `results/reference_list.json` names (its own sha256 is pinned there), under one of the board's names. Tried on
    the three local copies: all ten files of `builds_dc7b6ae` pass; the two heart files of `builds_38616fe` and
    `builds_fd2683a`, and the older T1 second file of `builds_fd2683a`, are refused. A format check is sent as
    before and prints a warning when the file is off the list. The pointer and the list are tracked files, so
    an edit to either also fails the clean-tree guard. A rehearsal that writes a new list rewrites the pointer
    in the same commit (`python3 -c` or by hand: `list`, `list_sha256` = `sha256sum` of the list file, `boards`).
    The guard is a second check, not a replacement for reading the 12 confirm characters against the list.
  - **The reference hashes (2026-11-13).** The rehearsal on the DGX at fd2683a
    (`experiments/builds/results/rehearsal_fd2683a.txt`; `docs/PLAN.md` decision log 2026-11-13 (e)) rebuilt all ten
    files twice from empty folders. The two sets are byte-identical, every stamp is clean, every file passed `vec
    check`, and the seven files whose commands did not change equal the archive of f337dc6. The reference for the
    freeze morning is `runs/archive/builds_fd2683a/SHA256SUMS` on the DGX. A hash made on another machine cannot stand
    in for it: the three T3 files that apply a per-gene factor have other bytes on the local machine (the draws of a
    few cells, through the BLAS kernel's effect on a nearest-neighbour step), as T1's s1 had until 2026-11-13.
  - **The freeze-morning rebuild.** On the DGX, synced to the freeze commit with no build running there and a clean
    tree: `bash experiments/builds/determinism_check.sh runs/archive/builds_fd2683a/SHA256SUMS 16`. It builds the three
    sets into the fresh folder `runs/determinism_<HEAD>/t16/`, runs `vec check` on every file and compares every file
    with the reference; it ends with "all files byte-identical to the archived set" or names the files that differ
    (tried at fd2683a on 2026-11-13: ten of ten OK). The uploads are pulled from that folder (`vec submit
    runs/determinism_<HEAD>/t16/builds/T1_s1.h5ad --board T1:test --from-dgx ...`); where this runbook or
    `docs/S2_RULES.md` names `runs/builds/...` or `runs/t3_build/...`, read the same name under that folder. The
    DGX's own `runs/builds` and `runs/t3_build` hold older files (T1's two there date from 2026-11-04; the embryo s2
    and the refused s1 predate the amendments) and are never an upload source: a clean stamp at an ancestor commit
    passes the client's stamp guard, so only the confirmed hash tells such a file from the current one. `build_all.sh`
    at its default folder would also stop at the first T3 file, which the T3 builder will not write over. A file that
    differs from the reference on the morning is not uploaded until the difference is explained to the captain. The
    server format checks for the three new files and T1's two were to wait for the freeze morning; they were run on
    2026-11-13 at the captain's word, with the two new heart files' (the evening amendment above).
    With `--from-dgx` an official upload also needs the DGX checkout at the local HEAD with a clean tree (the client's
    guard "dgx provenance"), so sync the DGX after the last local commit and before the upload, with no build running.
  - **A copy off the DGX (2026-11-13).** The local machine holds a verified copy of the reference set
    (`runs/archive/builds_fd2683a/`, ten of ten OK against its SHA256SUMS; the list is also in the repository as
    `experiments/builds/results/builds_fd2683a.SHA256SUMS`) and of the T3 run folder
    (`runs/t3_decomp/gata4_20261028T122425Z/`, identical file by file). The client's build-stamp guard and its local
    format check pass on all ten copies (local logic only; nothing was sent). If the DGX cannot be reached on an
    upload day, the same bytes can go from the copy, at the captain's word: `uv run vec submit
    runs/archive/builds_fd2683a/builds/T1_s1.h5ad --board T1:test ...` without `--from-dgx`; the hash to confirm is
    the same, the ledger then names fd2683a as the build commit, and the client does not ask for the DGX's state. This
    skips the freeze-morning rebuild, so it is for an outage only. A rebuild on the local machine is a weaker
    stand-in: it gives seven of the ten reference files byte for byte, and the three T3 files that apply a per-gene
    factor (allowed s1, refused s1, barred s1) with a few cells drawn differently (3, 84 and 14 cells), so with
    hashes the disclosure email does not list. If a build changes after 2026-11-13, copy and check the new
    reference the same day (done for the second rehearsal: `runs/archive/builds_38616fe/` on the local machine,
    ten of ten OK, the client's stamp guard and local format check passing on all ten; the first rehearsal's copy
    is still there and superseded). The copy sets aside, for these files, the "no lasting local copy" part of the
    choice of 2026-10-19 (uploads pull from the DGX); it was made after a side agent's note the captain forwarded,
    and the captain confirmed it on 2026-11-13 ("Keep it (Recommended)", `docs/PLAN.md` decision log 2026-11-13
    (k)). Since the third rehearsal the current copy is `runs/archive/builds_dc7b6ae/` (ten of ten OK). To remove
    the copies: `rm -r runs/archive/builds_dc7b6ae runs/archive/builds_38616fe runs/archive/builds_fd2683a
    runs/t3_decomp/gata4_20261028T122425Z` on the local machine, at the captain's word.
- Update `ASSUMPTIONS.md` (A1, A2, A3, A13, A19) and `progress.md`. A1 is confirmed (2026-10-28).
- **Disclosure ids for each official upload** (prepared 2026-11-05 by the AI assistant; the captain confirms them on
  the morning). The client refuses an official upload without `--disclose` ids from `EXTERNAL_DATA.md`. Rule used:
  every register row whose Tasks column names the board's task, since Rules s.10 asks for every public dataset,
  model and published code used, and a figure that informed a choice counts; disclosing a row that only informed a
  choice is safe, leaving one out is not. T1:test (s1 and s2): `X001,X002,X003,X005,X006,X007,X008,X009` (X007-X009 since R2's build, 2026-11-06 (k)). T2:heart:test_extrap and
  T2:embryo:test_interp: `X001,X004,X005`. T3:bcatenin: `X001,X004,X005`. A row added later joins the boards its
  Tasks column names.
  2026-11-10, working answers (`docs/ORGANIZER_QUESTIONS.md`; decision log 2026-11-10 (w), (x)), unless the screened check finds
  otherwise: the ids stay in the local ledger, and `EXTERNAL_DATA.md` asks for every source "with each submission", so a private
  email to the organisers goes **before the first official upload is pressed**, copied by Slack DM to Weize Xu (drafted by the
  AI assistant, sent by the captain). For every planned upload it gives the board, the file's SHA-256, the register ids and
  each external source's stages, with the removal of E10.5-E13.5 cells and the all-stage object's brief holding (X002). It says
  T3's s1 uses the two released Gata4 replicates and T3's s2 no number from them, and it corrects the note of 4 November on one
  point: T3's cell-count default, the board's maximum, was fixed on 2026-11-01, after the release and before that exposure. An
  upload it does not list gets one line before it goes. Every format check and upload carries the same
  `--architecture mithril-recipe-v1` and `--model sources-disclosed-by-email`; nothing about the methods goes in those
  fields, whose visibility is unknown. If the screened read of `docs/PLAN.md` decision log 2026-11-11 (c) returned a
  cardiomyocyte or blood share, the email also says which T1 upload or uploads take the share stated on the organisers'
  Task 1 page, the date it was read, and that it is used without their confirmation (even if no file changed). The
  read of 2026-11-11 returned none (decision log 2026-11-11 (d)), so this does not apply.
- **(user) Before each second upload: the captain's checklist** (2026-11-18, written by the AI assistant after the
  first uploads skipped the email, the screened check and the reference list; `docs/PLAN.md` decision log 2026-11-18
  (j)). Uploads go by the captain's hand through the site's form, so the client and its guards do not run: each guard
  is a line here. Lines 1 and 4 to 8 need nothing beyond `sha256sum` and `grep`; line 2 needs the repository's `vec
  check` (local, nothing sent); line 3 uses the submission client, format-only, and only where the captain asks for it. One board at a time; a line that cannot be ticked stops
  that upload until the captain has the explanation. Deadline 2026-11-26 11:59 UTC, never in the last hours.
  1. [ ] **The hash.** On the machine the file goes from: `sha256sum <file>`. All 64 characters equal the board's line
     in the reference list: `experiments/builds/results/builds_dc7b6ae.SHA256SUMS` (T1, the embryo, T3; the list
     `results/reference_list.json` points to), or for the heart's half
     `experiments/builds/results/heart2_307c9b2.SHA256SUMS` (`builds/T2heart_s2iso_half.h5ad`). That second list is
     outside the pointer and was built in a campaign checkout (decision log 2026-11-18 (n)): the heart's line is
     ticked only once the captain has accepted that provenance or a rebuild in the reference checkout has matched it.
     A file in no committed list (a new embryo file, entry (s), item 1) needs its list committed, and the
     `docs/S2_RULES.md` amendment it requires, before this line can be ticked. By hand: `grep <the 64 characters>
     <list>` prints one line, under the name given above. A file whose 64 characters differ is not the listed
     file, wherever it was built. The same 64 characters go in the disclosure note (line 4).
  2. [ ] **Local format check of these bytes:** `uv run --frozen --no-sync vec check <file> --board <KEY>` on the same
     path as line 1 passes (no token; nothing is sent).
  3. [ ] **Server format check, only with the captain's word for this check, given each time.** A server format
     check can come back stored as a submission: it may take the board's last upload and blocks every board until it
     is resolved (section 0). Nobody pastes or quotes the log's server text into a session. Before it, `uv run --frozen
     --no-sync vec submit status > <log>` shows nothing `unresolved` or `ambiguous` (section 0). Then, format-only by
     default and never `--official`: `uv run --frozen --no-sync vec submit <local file> --board <KEY> > <log>`, or for
     a DGX path `uv run --frozen --no-sync vec submit <path on the DGX> --board <KEY> --from-dgx > <log>`; output to a
     file, never through `head`. After it, `vec submit status` again shows nothing unresolved and the portal's list
     shows no new scored entry for the board; if either does, stop and resolve it (`vec submit resolve ...`, section
     0) before anything else. A pass is not a promise: the heart's s1 passed on 2026-11-13 and the form refused it for
     its cell count on 2026-11-18.
  4. [ ] **The disclosure note went first.** The captain has sent the email or one line that names this board, this
     file's sha256 and its register ids (above), with the board's own note where one is drafted
     (`docs/ORGANIZER_QUESTIONS.md`: "Note before the heart's second upload", "Note before Task 3's second upload").
  5. [ ] **The screened check** of mail, Slack and the Rules page ran on the upload day, by a screening agent, audited,
     with the switches fixed beforehand (`docs/ORGANIZER_QUESTIONS.md`; for T3 `docs/S2_RULES.md`, last section: the
     barred s1 if a ruling restricts the Gata4 files or bars the validation answers). Nothing it found switches this
     file. No organiser reply is pasted into the session.
  6. [ ] **The warnings register** for the file, in the decision-log entry that hands it over: every committed
     reading that puts any step of this file below the floor on any metric (a practice row, a hedge run, a failed
     gate, a flagged geometry, a step no practice run applied), each with its answer: why it does not apply, what
     hedges it, or "unanswered". The captain has read it; an unanswered line has been put to him by name.
  7. [ ] **The go-ahead:** the captain's own word for this board and these 64 characters.
  8. [ ] **Afterwards, paste back** into the main session: the board, the time and file name the site lists, the first characters of the hash the site shows, the outcome (scored, or "rejected" and the kind of rejection in
     your own words, with no number from the message), the board total and every per-metric value. The main session
     records them in `results/manual_uploads.csv` (one row; the board total) and in a decision-log entry in
     `docs/PLAN.md` (the per-metric values). It checks the copy only (that the board total is the one the listed
     metrics give) and never differences values across uploads. A score names the banked file and chooses the final
     mark; it sets nothing else (`experiments/harness/TRIES.md` rule 9, J8; "Scores choose, never compute").
- **(user)** Final entries. The account page (seen in the site's public script, 2026-10-28): "The leaderboard shows
  your team's best score on each board, so there is nothing to nominate yet. In the final phase the captain picks
  which submission counts." While nominations are open, the captain sees a "Final" column in the team's submission
  list and marks a scored submission (or a bundle) as final. The timeline: 2026-12-02 is "the last moment a submission
  can be made or a final entry nominated. Nothing after this counts, in either track." Rules s.11 still says a board
  ranks a team by its best; what counts when no entry is marked is not stated (asked in the fourth organiser email,
  2026-10-28). Mark the final entries after the last official upload and before 2026-12-02 (the user is the captain).
  Working rule (2026-10-31, red team): mark each board as soon as its last official upload is scored; of the entries
  every ruling received by then allows, mark the one with the higher published score (a tie, or s2 unscored by
  30 Nov IST: s1); if marks are per bundle, the bundle with the higher overall standing. Look for the Final column
  after every official upload, ask on the email thread and in #virtual-embryo-challenge if it has not appeared by
  28 Nov IST, and recheck all four boards on 1 Dec IST.

## Never today

- No official upload on any board while the captain's upload gate stands (2026-11-05: "Until our estimated score in
  leaderboard crosses 209, don't submit our 1st submission"; `docs/PLAN.md` section 8), until his deadline exception
  opens it at 2026-11-27 00:00 UTC. The client refuses one unless `results/upload_gate.json` passes (guard "upload
  gate").
- No scored upload (format checks only; `format_only=true`). `vec submit` sends format checks unless given
  `--official` with the file's hash, which needs the user's go-ahead.
- Server format checks only on boards of the phase's split: in the final phase, the test boards that `vec boards
  update` adds (the client refuses validation boards then); `vec check` works locally on any board.
- Official uploads (later) need every stream's work committed: an untracked file outside results/, runs/ and
  experiments/**/results/ blocks them.
- No use of anything read from a test board beyond its submission spec.
- No external data that has not been registered in `EXTERNAL_DATA.md` with committed filter code.
