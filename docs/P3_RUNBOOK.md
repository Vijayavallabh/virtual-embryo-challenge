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
2. **T3 decomposition on Gata4:** the single command in `experiments/t3_decomp/` with both replicates and WT E8.75.
   Commit its decision file. The rules in `experiments/t3_decomp/rules.yaml` decide; nothing is chosen by eye.
3. **T1 direction readout:** the command in `experiments/t1_direction/` on E9.5 to E10.5. Commit its decision
   (build the bracketing pipeline or not).
4. **Harmonisation:** `experiments/harmonise/run_harmonise.py` release stages (Gata4 replicate 1 vs 2, cross-stage
   transfer with E8.5 and E10.5). The replicate difference is the noise bar for T3 composition effects (A19).
5. **Density:** rerun `experiments/density/run_all.sh` with the E10.5 and E7.5 targets; it decides the count check in
   `docs/PLAN.md` section 3.6.
6. **Scorer build:** `experiments/scorer_build/anchor_check.py` against the released validation answers (validation
   anchors only).

## 5. Then

- Emulator calibration against the published validation reference rows (feat-005, assumption A2).
- Test-board intake: reference stage, cell limits and what else each test board publishes (decision point D1).
- Update `ASSUMPTIONS.md` (A1, A2, A3, A13, A19) and `progress.md`.

## Never today

- No scored upload (format checks only; `format_only=true`).
- No use of anything read from a test board beyond its submission spec.
- No external data that has not been registered in `EXTERNAL_DATA.md` with committed filter code.
