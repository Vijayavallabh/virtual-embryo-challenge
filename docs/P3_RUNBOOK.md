# Release-day runbook (2026-10-27, start of the final phase)

The organisers moved the start from 2026-10-20 to 2026-10-27 (email of 2026-10-20). Files written before that email say
2026-10-20; they mean this release.

The validation answers (T1 E10.5, T2 heart E8.5 and E10.5, T2 embryo E7.5, T3 Gata4 replicates) and the test inputs
arrive today. Everything below that decides a method was committed before the release; run it as written and commit
each output before looking at anything else. Steps marked **(user)** need you.

## 0. Before anything

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

```bash
uv run vec boards update      # same as vec.boards.update_boards(); test anchors stripped in memory, never printed
git add boards && git commit -m "Board specs for the final phase (test anchors stripped at intake)"
```

The summary lists the boards, which ones were reduced to their submission spec, and whether the T1 composition file was
written (it is skipped if it describes any unreleased stage). Never print `boards/index.json` anchors for a reduced
board; they are not there. E10.5 counts as unreleased until the public phase endpoint reports the test split; the
summary says why the file was held back. If step 0 showed the final phase and the summary says the endpoint could not be
read, re-run `uv run vec boards update` before committing.

## 2. Data

```bash
uv run vec data manifest      # local: needs the team token; lists the new files
uv run vec data fetch         # local: streams each signed link to curl on PrakashDGX_H2 (default --dest data/raw there)
ssh PrakashDGX_H2 'cd ~/v/virtual-embryo-challenge && export PATH=$HOME/.local/bin:$PATH && \
  uv run vec data inventory && uv run vec data qc data/raw/T*/*.h5ad --csv runs/qc.csv'
```

Check sizes and checksums against the manifest, record the new files and the QC table (z levels, steps, duplicated
records, row order) in `docs/DATA_NOTES.md`, and commit.

## 3. Labels (before any method touches a new file)

Map every new label in `configs/lineage_map.yaml` with the rule in its header (core = cardiomyocytes and endocardium;
ambiguous labels default to shell after a marker check with `scripts/label_markers.py`). Commit the map. The size
procedure refuses unmapped labels.

## 4. Pre-registered runs, in order (commit each output with the git hash)

1. **Heart and embryo size:** `bash experiments/coords/run_p3.sh` (paths as arguments). Commit
   `experiments/coords/results/p3/`. The heart JSON sizes both heart submissions through `vec.t2_build.size_or_hold`
   (approved rule 3, `docs/RELEASE_RULES_SUMMARY.md`): as recorded when `valid` is true; when `valid` is false while
   E10.5 is usable, both keep E10.5's size. If E10.5 is barred (A20), the size rule goes back to the user. An embryo
   JSON with `valid` false is not used.
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
   transfer with E8.5 and E10.5). The replicate difference is the noise bar for T3 composition effects (A19).
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

- Emulator calibration against the published validation reference rows (feat-005, assumption A2).
- Test-board intake: reference stage, cell limits and what else each test board publishes (decision point D1).
- Candidate files, built on the DGX with the committed builders and checked there with
  `uv run vec check FILE --board KEY` (no token, nothing sent). A server format check runs from this machine with
  `uv run vec submit runs/<...>.h5ad --board <test board key> --from-dgx` (the file is pulled for that call only;
  nothing official without `--official`). Nothing is uploaded officially today. Every built file gets a build stamp
  (`<file>.stamp.json`) that an official upload checks: commit before building, and do not sync, pull or check out
  on the DGX while a build runs there (the stamp then is not clean; rebuild after the sync).
  - T3: `uv run python experiments/t3_decomp/build_submission.py build <gata4_decisions.json> --branch <b>
    --submission s1|s2 --wt data/raw/T2/E8.75.h5ad --run-dir runs/t3_decomp/gata4_<STAMP> --board <T3 board key>
    --n-cells <the board's max_cells> --out runs/t3_build/<b>_<s>.h5ad`, for `allowed` and `barred` until the
    organisers answer follow-up question 3 (then only the branch their answer selects). It re-decides nothing and
    refuses a decision file whose script or rules hash differs from the files on disk.
  - T2: `vec.t2_build`, run from the repo root (the paths are relative). Count by `n_cells_from_decision(board,
    pool.n_obs, "experiments/density/results/release/count_decision.json")` (step 4.5; null keeps max_cells; a
    decided level below the test board's lower limit gives way to the next passing count up, the passing level with
    the fewest cells the board allows, else max_cells, and the count record's `count_floor` says so; the user approved
    the words on 2026-10-21 and left the reading of "next" to the AI assistant, whose decision panel chose the count
    up: `experiments/density/floor_reading/README.md`). If step 4.5 wrote no decision because the rerun left no rule
    rows, approved rule 7 keeps max_cells: count with `n_cells(board, pool.n_obs)` and
    record why in `progress.md`. Then `sample_cells`. Heart size by `size_or_hold(coords,
    "experiments/coords/results/p3/size_heart_E12.5.json", "H1" (s1) or "H2" (s2), source_stage="E10.5",
    reference_cloud=<the full E10.5 pool from load_stage(..., drop_duplicates=True), before any scope filter>,
    expect_board="heart")`; embryo size by `apply_size_procedure` with
    `experiments/coords/results/p3/size_embryo_E7.75.json`. Then `cloud_guard` and `build_submission`. The module
    docstring of `src/vec/t2_build.py` shows these calls in order; `experiments/t2_build/smoke.py` counts with plain
    `n_cells` and sizes with `rescale_to_rms`, so copy only its `sample_cells`, `cloud_guard` and `build_submission`
    calls.
  - T1: `vec.t1_scope.build_submission` once the shares are fixed by their written rule (docs/PLAN.md section 4.1);
    pass `ref_source` (the reference file's stem) so the per-file label remaps apply.
- Update `ASSUMPTIONS.md` (A1, A2, A3, A13, A19) and `progress.md`.

## Never today

- No scored upload (format checks only; `format_only=true`). `vec submit` sends format checks unless given
  `--official` with the file's hash, which needs the user's go-ahead.
- Server format checks only on boards of the phase's split: in the final phase, the test boards that `vec boards
  update` adds (the client refuses validation boards then); `vec check` works locally on any board.
- Official uploads (later) need every stream's work committed: an untracked file outside results/, runs/ and
  experiments/**/results/ blocks them.
- No use of anything read from a test board beyond its submission spec.
- No external data that has not been registered in `EXTERNAL_DATA.md` with committed filter code.
