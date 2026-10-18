# Release-day runbook (2026-10-20, start of the final phase)

The validation answers (T1 E10.5, T2 heart E8.5 and E10.5, T2 embryo E7.5, T3 Gata4 replicates) and the test inputs
arrive today. Everything below that decides a method was committed before the release; run it as written and commit
each output before looking at anything else. Steps marked **(user)** need you.

## 0. Before anything

- **(user)** Do not open `panels/index.json` on the challenge site or the public board endpoints, and skip Slack threads
  that discuss sizes or other values worked out from anchors. Test-board anchors, if published, encode the hidden
  targets (`docs/DATA_POLICY.md`, do-not-use list).
- `git pull`, `./init.sh` (lock, lint, tests).

## 1. Board specs (anchors stripped in memory)

```bash
uv run python -c "from vec.boards import update_boards; print(update_boards())"
git add boards && git commit -m "Board specs for the final phase (test anchors stripped at intake)"
```

The summary lists the boards, which ones were reduced to their submission spec, and whether the T1 composition file
was written (it is skipped if it describes any unreleased stage). Never print `boards/index.json` anchors for a reduced
board; they are not there.

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

1. **Heart and embryo size:** `bash experiments/coords/run_p3.sh` (paths as arguments). Use the heart output only if
   `valid` is true. Commit `experiments/coords/results/p3/`.
2. **T3 decomposition on Gata4.** First count identical count rows shared by the two replicate files (rules.yaml
   note 4: more than 1% shared means the borrowed-replicate rules are recorded as not independent). Then:
   `tmux new -d -s vec-t3-gata4 'bash experiments/t3_decomp/run_gata4.sh REP1.h5ad REP2.h5ad data/raw/T2/E8.75.h5ad
   --e85 <released E8.5 file>'`. Check that `runs/logs/t3_decomp_gata4_*.log` exists (a wrong path exits 0 silently).
   Bring back the aggregates with the rsync line the run prints and commit `gata4_decisions.json` first. The rules in
   `experiments/t3_decomp/rules.yaml` decide, including the lead's dated notes; nothing is chosen by eye.
3. **T1 direction readout:** `uv run python experiments/t1_direction/readout.py <released E10.5 file>` (DGX). It
   refuses uncommitted procedure files, a non-default config or label map that is not a committed repo file, and a
   second run without `--overwrite`. Commit `decision.json` (build the bracketing pipeline or not). Re-run the
   committed filters (`experiments/external/`) before ~/v/vec-gapconc is cleaned up.
4. **Harmonisation:** `experiments/harmonise/run_harmonise.py` release stages (Gata4 replicate 1 vs 2, cross-stage
   transfer with E8.5 and E10.5). The replicate difference is the noise bar for T3 composition effects (A19).
5. **Density:** `TARGETS="heart_E10.5 embryo_E7.5" SEEDS="0 1 2 3 4 5 6 7 8 9" RAW=runs/density/raw_release
   TAG=density_release JOBS=16 THREADS=1 bash experiments/density/run_all.sh` (about 20 minutes), then `summarise.py`.
   The pre-registered rule in `docs/PLAN.md` section 3.6 decides the count: keep max_cells unless a lower level beats it
   for the intermediate-quality condition by more than 0.3 points in at least 8 of 10 seeds.
6. **Scorer build:** `uv run python experiments/scorer_build/anchor_check.py check` against the released validation
   answers (validation anchors only; it refuses test boards, the public board file and anything before 2026-10-20).
   Its verdict says which build the server baked the anchors with; it changes no prediction.

## 5. Then

- Emulator calibration against the published validation reference rows (feat-005, assumption A2).
- Test-board intake: reference stage, cell limits and what else each test board publishes (decision point D1).
- Update `ASSUMPTIONS.md` (A1, A2, A3, A13, A19) and `progress.md`.

## Never today

- No scored upload (format checks only; `format_only=true`).
- No use of anything read from a test board beyond its submission spec.
- No external data that has not been registered in `EXTERNAL_DATA.md` with committed filter code.
