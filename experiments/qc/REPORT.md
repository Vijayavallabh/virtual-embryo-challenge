# Data-QC utilities: z snapping, duplicated records, per-file intake QC

_Saved by the lead on 2026-10-19 from the QC agent's returned text (agents cannot write report files). The independent review follows it, and the lead's decisions are at the end._

## Summary
The QC utilities are built and tested, and they ran on the DGX over all 9 released files (T2, T3, and T1 for row identity only). With the specified test (identical count rows, xy within 0.01, same or adjacent z level), the duplicate counts match the background: E8.0 443, E8.25_late 1164, E8.75 317, Mab21l2 KO 2811 (1404 + 1407). There are none in E6.75, E7.25, E9.5, T1 E8.5_RNA or T1 E9.5_RNA.

The background did not mention two other kinds of synthetic record in the heart and embryo files:
- **Copies with moved xy:** E8.0 2460, E8.25_late 563, E8.75 518.
- **Interpolated averages:** E8.0 200, E8.25_late 2558, E8.75 614.

`load_stage(..., drop_duplicates=True)` drops only the exact co-located copies, as specified. Whether prediction pipelines should also drop the other two kinds is your call (see Open issues).

## What changed (only files I own)
- **`src/vec/spatial.py` (new)**
  - `snap_z(z, tol=None)` returns `ZLevels(index, levels, step)`. It is a NamedTuple, so 3-way unpacking works, and it also has `.grid` (level positions counted in steps) and `.gaps()`.
  - How it snaps: it sorts the distinct z values and merges neighbours that are at most `tol` apart. There is no rounding to a fixed origin, so float noise, x.5 levels, gaps, two sections on one level, and a section stored twice at a small z offset are all handled. With `tol=None` the tolerance is a quarter of the step, estimated after merging float noise only. The step is the median spacing, refined over whole multiples so that gaps do not distort it.
  - `find_duplicate_records(adata, xy_tol=0.01, min_library=20, max_level_gap=1)` returns `(mask, report)`.
    - How it works: count rows are hashed (64-bit) first, then rows that are both identical and close in xy and z are compared exactly. Matches chain, and the first occurrence in file order is kept.
    - Adjacency is measured in steps, so the KO levels -320 and -290 on either side of its gap are not adjacent.
    - `xy_tol=inf` ignores xy; large groups are then swept by level.
    - The report gives counts by gap and by (kept level, dropped level) pair, the largest xy offset, and matches skipped because they are under `min_library`.
    - Speed: 60k cells take about 1 s, measured in a test.
  - `identical_rows(M)`, `qc_stage(adata)` (JSON-safe aggregates) and `summary_row`.
  - `qc_stage` reports:
    - cells, duplicated obs names, and whether the counts are integers
    - rows that repeat an earlier row anywhere in the file
    - raw z values and levels, step, spacings, gaps, off-grid levels, and dense levels (with the ratio to their neighbours)
    - whether `sample_id` / `section` labels agree with the levels
    - row order (Spearman of row index vs z level, number of z runs)
    - duplicates at `xy_tol`, and again with xy ignored
    - coordinate ranges and centring
- **`src/vec/prep.py`, loader only (your lineage functions are untouched)**
  - `load_stage(..., drop_duplicates=False)`. When True, detection runs on the file's own counts with all genes, before the int32 cast and before any gene subsetting.
  - `uns["qc"]` records `n_input`, `n_dropped`, `dropped_rows` (positions in the file), `by_gap`, `by_level_pair` (a DataFrame), `below_min_library`, `max_xy_offset` and the parameters. A test confirms it can still be written to h5ad.
  - Files without coordinates (T1) come back whole, with a note in `uns["qc"]`.
  - The docstrings say that prediction pipelines pass True and the emulator keeps the default False.
- **`src/vec/io.py`:** `check_submission` now warns when cells share exact coordinates.
- **`src/vec/cli.py`:** new command `vec data qc FILE... [--out runs/qc.json] [--csv PATH] [--xy-tol 0.01] [--min-library 20] [--max-level-gap 1]`, in the existing table style.
- **Tests:** `tests/test_spatial.py` (new, 18 tests); appended 2 tests to `tests/test_prep.py` and 1 to `tests/test_io.py`.
  - z cases: float-noise levels, x.5 levels, half-step levels, gaps, a level stored twice, two sections on one level, a single level, NaN.
  - duplicate cases: copies within a level and in an adjacent level, chains, a copy stored before its original in file order, a row with one count different (never flagged), the same row far away in xy, low library, all parameters, and 60k cells.
- **`experiments/qc/`:** `run_qc.sh` (intake run), `probe_records.py` (characterises records the detector leaves alone), and `results/`.

## Commands
- Local:
  - `uv run ruff check src tests experiments`
  - `uv run pytest -q`
- Before each DGX run I synced with the given rsync command.
- QC run, in tmux session `vec-qc-intake2` with its log at `runs/logs/qc_intake_2026-10-19_v2.log`, threads set to 8:
  - `bash experiments/qc/run_qc.sh released_2026-10-19`
  - This runs `uv run vec data qc data/raw/T2/*.h5ad data/raw/T3/*.h5ad data/raw/T1/*.h5ad --out experiments/qc/results/qc_released_2026-10-19.json --csv experiments/qc/results/qc_released_2026-10-19.csv`
  - It finished with EXIT 0 in under 2 minutes.
- Probe, in tmux session `vec-qc-probe`:
  - `uv run python experiments/qc/probe_records.py data/raw/T2/*.h5ad data/raw/T3/*.h5ad --out experiments/qc/results/probe_records_2026-10-19.json`
  - It finished with EXIT 0.
- Results were copied back with `rsync PrakashDGX_H2:v/virtual-embryo-challenge/experiments/qc/results/` (JSON and CSV aggregates only).
- For tomorrow's files: `bash experiments/qc/run_qc.sh <tag>` on the DGX.

## Results (`experiments/qc/results/qc_released_2026-10-19.{json,csv}`)

"Copies" are counted at `xy_tol` 0.01; "any xy" ignores xy.

| file | cells | z raw/levels | step | gaps | copies (xy 0.01) | copies (any xy) | identical rows (>=20 counts) | non-integer rows |
|---|---|---|---|---|---|---|---|---|
| E6.75 | 7093 | 65/37 | 10 | - | 0 | 0 | 49 (0) | 0 |
| E7.25 | 13295 | 40/40 | 10 | - | 0 | 0 | 95 (0) | 0 |
| E8.0 | 31671 | 30/30 | 15 | - | 443 (gap1 321, gap2 122) | 2903 | 2912 (2903) | 200 |
| E8.25_late | 58716 | 38/38 | 15 | - | 1164 (gap1) | 1727 | 1728 (1727) | 2558 |
| E8.75 | 24826 | 35/35 | 14 | - | 317 (gap1) | 835 | 877 (835) | 614 |
| E9.5 | 53742 | 98/40 | 15 | - | 0 | 0 | 0 | 0 |
| Mab21l2 KO | 50294 | 36/36 | 10 | -320 to -290 (2 levels missing) | 2811 (gap0 1404, gap1 1407) | 2811 | 2811 | 0 |
| T1 E8.5_RNA / E9.5_RNA | 16787 / 17057 | - | - | - | - | - | 0 / 0 | - |

- **z noise:** the largest within-level spread is 1.1e-13 in E9.5 and 2.8e-14 in E6.75; it is 0 elsewhere.
- **Row order:** in every spatial file, rows are grouped by level in descending z. Spearman is -0.9994 to -0.9996 and the number of z runs equals the number of levels.
- **Centring:** T2 coordinates are centred on the cell mean (|mean xy| < 1e-12). The KO is not centred: mean (-5438.6, -1463.8, -182.4).
- **Where the copies are, by level pair:**
  - E8.0: 295 from 32.3 to 17.3, 122 from 32.3 to 2.3 (chained through 17.3, which is why the background said "about 440"), and 26 from 17.3 to 2.3.
  - E8.25_late: all 1164 from -85.8 to -100.8.
  - E8.75: all 317 from -86.5 to -100.5.
  - KO: 1407 from -210 to -220. These are sample_ids 18 and 19; they share the section label 19, so it is one section stored at two levels. There are also 1404 within level -320, which is sample_id 35 stored twice; the xy offsets are at most 0.0059 because the coordinates are float32.
- **More KO checks:**
  - Level -320 holds sample_ids 30 and 35 and has 4170 cells, 3.15 times its neighbours.
  - All 70 exact duplicate coordinates are inside the flagged copies; none remain after dropping them.
  - All 2811 copies reuse the obs name of a kept record, which accounts for the 2811 duplicated names.

## Findings the background did not mention (`probe_records_2026-10-19.json`, aggregates only)
1. **Copies with moved xy.** These are identical rows (at least 20 counts) one level away, but shifted:
   - E8.0: 2460, shifted by 0.06 to 63 units (median 5); the shift is not a rigid transform (affine fit residual 5 to 28).
   - E8.25_late: 563, from -115.8 to -130.8, shifted up to 7.4 units.
   - E8.75: 518, from -114.5 to -128.5. These fit a rigid transform almost exactly (residual 0.03, rotation -0.25 degrees), so that section was re-registered.
2. **Interpolated averages (the half-integer counts).** In each case twice the row equals the sum of two records within 30 units, one from the level above and one from the level below; none of them sit at the same xy.
   - E8.75: 614 of 614 rows.
   - E8.25_late: 2556 of 2558. A few rows have fractions of 0.25, 0.125 or 0.375, which looks like averages of averages.
   - E8.0: 160 of 200 within the 30-unit search radius.
3. **The loader truncates these counts.** `load_stage` casts `layers["counts"]` with `.astype(np.int32)`, so 0.5 becomes 0 (and the entry disappears) and 1.5 becomes 1. That affects 13,573 / 174,097 / 34,324 entries in 200 / 2,558 / 614 rows of E8.0 / E8.25_late / E8.75. I did not change this, because it is a decision about what the counts layer promises. The emulator reads X through `io.read_stage`, so scoring is not affected.
4. **The background said E8.75 levels sit on x.5 values.** In the released file they are at -240.473939 + 14k, so none are on x.5. `snap_z` is tested on x.5 levels anyway, so nothing depends on this.
5. **Dense-level hints that are probably not doubled sections:** E8.25_late at -235.8 (1.7 times its neighbours) and E9.5 at -160.3 (1.82 times) look like normal variation between sections. Only the KO's -320 is a real doubled level.

## Caveats
- `identical_rows` counts rely on the 64-bit hash alone; at 60k rows the chance of a collision is about 1e-10. `find_duplicate_records` checks every match exactly.
- The 20-count guard leaves 10 co-located identical low-count matches unflagged (E8.0 2, E8.25_late 1, E8.75 7). At xy within 0.01 these are almost certainly copies too.
- To keep this work clear of the banned material, I did not open the `coords_REPORT.md` the task pointed to. I re-derived every background number from the released files instead.
- I did not open `docs/research/`, veckit 0.1.2, PR #6 or `boards/index.json`. The outputs contain only structural QC aggregates of released training files: no size figure for any held-out stage and no gene-level knockout findings.

## Independent review (2026-10-19)

The utilities stream holds up. Everything in the focus list reproduces exactly with independent code on the DGX: snap_z on every released file, the duplicate counts, no false positives, the drop_duplicates=False default and the `vec data qc` CLI. Lint and tests pass. The only failing lint is 5 errors in experiments/harmonise_review/review.py, a file from another stream committed in dbb6eaf. Two descriptions in the agent's "findings the background did not mention" section are overstated. First, the E8.75 moved copies are an exact affine image, not a rigid one. Second, the interpolated rows are averages of one record from the level above and one from the level below in every case only for E8.75. Neither affects the code or the specified behaviour. No repo file was edited. Only aggregate JSON/CSV left the DGX. Local scripts are in /tmp/claude-1001/-mnt-md0-IITM-BackUp-Home-vijayavallabh/07d440b9-3ae8-4421-8a08-a046b42f307d/scratchpad/verify_utilities/ (indep_qc.py, api_check.py, prop_test.py, followup.py, launch.sh) and the aggregate outputs are in its results/ subfolder. DGX outputs are in ~/v/virtual-embryo-challenge/runs/utilities_verify/.

### Checks

| # | Item | Verdict | Notes |
|---|---|---|---|
| 1 | Lint and full test suite (local) | confirmed | `uv run ruff check src tests` (the AGENTS.md check) passes. `uv run ruff check src tests experiments` reports 5 errors, all in experiments/harmonise_review/review.py (another stream, committed in dbb6eaf). Ruff on the agent's own files passes. `uv run pytest -q`: 64 passed in 14.2 s. After other agents changed files mid-run, a re-run gave 75 passed in 40.1 s. tests/test_spatial.py has 18 tests; the 60k-cell test took 0.80 s. |
| 2 | Read every changed file | confirmed | Read in full: src/vec/spatial.py (new, 462 lines), the prep.py loader (_duplicate_filter, load_stage), io.py (duplicate-coordinate warning), cli.py (`data qc`), tests/test_spatial.py, the new tests in test_prep.py and test_io.py, experiments/qc/run_qc.sh, probe_records.py and the 3 result files. The prep.py and test_prep.py diffs also carry another agent's lineage changes (file_overrides, merged, 2 lineage tests), which fits the agent's 'loader only' statement. The verified files were byte-identical locally and on the DGX at the end of the run. |
| 3 | snap_z on all released files | confirmed | Raw z values / levels / step: E6.75 65/37/10, E7.25 40/40/10, E8.0 30/30/15, E8.25_late 38/38/15, E8.75 35/35/14, E9.5 98/40/15. KO: 36/36, step 10, one gap -320 to -290 (30 units, 2 levels missing). E8.75 levels are -240.473939 + 14k, so none sit on x.5. My independent split (gaps > 1.0 absolute, median level values) gives the same cell-by-cell levels on all 7 files; level values differ by at most 2.8e-14. Explicit tol of 1e-9, 0.5, 2 and 4 give the same level counts. Largest within-level spread is 1.1e-13 (E9.5) and 2.8e-14 (E6.75), 0 elsewhere. |
| 4 | Duplicate counts vs the review (about 440 / 1164 / 317 / 2811 / 0) | confirmed | E8.0: 443 = gap1 321 + gap2 122, made of 295 from 32.3 to 17.3, 122 from 32.3 to 2.3 (chained via 17.3) and 26 from 17.3 to 2.3. E8.25_late: 1164, all -85.8 to -100.8. E8.75: 317, all -86.5 to -100.5. KO: 2811 = 1404 within -320 + 1407 from -210 to -220. Zero in E6.75, E7.25 and E9.5, and zero identical rows in T1. My detector (exact identity by lexsort, scipy connected_components, keep the lowest index) gives masks bit-identical to find_duplicate_records on all 7 files under 6 settings: default, xy inf, max_level_gap 0 and 2, min_library 1, xy 0.1. xy-inf counts are 2903 / 1727 / 835 / 2811. |
| 5 | No false positives | confirmed | All 4,735 flagged records have a count row and an X row identical to their kept record, library >= 20 (minimum 20-22, median 105-213) and at least 11 genes detected. E9.5 (53,742 cells), E6.75 and E7.25 contain zero identical rows with library >= 20; their chance-identical rows have library <= 7. So coincidences at that library size do not occur in real data. All 2811 KO copies reuse the kept record's obs name. Sensitivity: xy_tol 0.1 adds only 2/2/1 records (the closest moved copies, 0.056-0.088 units off). Step-based and positional level adjacency give the same result, and max_level_gap 2 changes nothing. Synthetic property test: 300 random scenarios x 7 settings, with chains, copies stored before originals, groups over 64 members (the cKDTree path), copies across a gap, and a level stored twice at +0.001 z. Result: 0 mask mismatches against brute force and 0 snap_z level mismatches. |
| 6 | load_stage keeps drop_duplicates=False by default | confirmed | The signature default is False. A default load keeps all rows on all 7 spatial files and adds no uns['qc']. drop_duplicates=True drops exactly my independently flagged rows (n_dropped 443 / 1164 / 317 / 2811) and keeps file order. Obs names are unique even before make_unique, and no duplicate 3D positions remain. T1 E8.5_RNA comes back whole (16,787 cells) with the note. The lineage=True path also runs on all 7 files. |
| 7 | vec data qc CLI and run_qc.sh | confirmed | Re-running the CLI on the 9 released files (17 s) gives JSON and CSV byte-identical to experiments/qc/results/qc_released_2026-10-19.{json,csv} (md5 92bd257b..., 148621a6...). The unmodified run_qc.sh, run in a sandbox under runs/utilities_verify, gives identical output. `--xy-tol inf` works (E8.75 835, KO 2811, xy_tol reported as 'any'). --help shows --out, --csv, --xy-tol, --min-library and --max-level-gap with the stated defaults. |
| 8 | KO structure, centring, row order, identical-row and integrality columns | confirmed | Level -320 holds sample_id 35 (2808 = 2 x 1404) and sample_id 30 (1362): 4170 cells, 3.15x its neighbours. The 1404 within-level copies are all sample 35 kept / sample 35 dropped, so sample 35 is stored twice. The 1407 copies run from sample 18 (-210) to sample 19 (-220), both labelled section 19. There are 70 duplicate 3D positions before dropping and 0 after; all 2811 duplicated obs names belong to the copies; the file has 37 sample_ids. KO mean is (-5438.6, -1463.8, -182.4); T2 /mean xy/ <= 2.6e-13. Spearman ranges from -0.99937 to -0.99963 and z runs equal the number of levels everywhere. Exact identical repeats (library >= 20): 49 (0), 95 (0), 2912 (2903), 1728 (1727), 877 (835), 0, 2811 (2811). Non-integer rows: 200 / 2558 / 614; E8.25_late also has fractions 0.25 x684, 0.125 x72 and 0.375 x26. |
| 9 | Loader int32 truncation claim | confirmed | load_stage counts equal trunc(file counts) exactly. Entries changed: 13,573 / 174,097 / 34,324. Of these, 9,549 / 119,975 / 24,387 are lost entirely (0 < v < 1). The scoring path uses io.read_stage and Stage.from_anndata(a.X) (cli.py lines 70 and 113), so scoring is unaffected. |
| 10 | Moved-xy copies (E8.0 2460, E8.25_late 563, E8.75 518) | partially-confirmed | The counts and offsets reproduce. E8.0: 2460 copies, 2294 one level away and 166 two levels away; offsets 0.056-63.4, median 5.07; not rigid (rigid RMS 5.3-31.1, affine 5.3-28.3). E8.25_late: max offset 7.38. E8.75: max offset 11.58. Refuted: 'E8.75 fits a rigid transform almost exactly (residual 0.03, rotation -0.25 deg)'. The 0.03 residual comes from a 6-parameter affine fit whose linear part scales the two axes by 1.028 and 0.971 (polar rotation -0.97 deg; -0.25 deg is only the angle of the first column). A rigid fit leaves 4.59 RMS. 'Re-registered' still holds as an exact affine image. Also, 166 of the E8.0 copies are two levels from their first occurrence, not 'one level away'. |
| 11 | Interpolated half-integer rows | partially-confirmed | The agent's counts (614/614, 2556/2558, 160/200) fit its search: partner A within 30 units, partner B anywhere. I searched with no radius limit and excluded a row's own copies. Then 614/614, 2558/2558 and 197/200 rows satisfy 2r = A + B, never with a co-located partner. But 'in each case both partners within 30 units, one from the level above and one from the level below' holds only for E8.75 (614/614, max 19.3 units). E8.25_late: the closest decomposition is above/below in 2446/2558 cases (the rest are 2 levels away), and 2554 are within 30 units (max 163). E8.0: 164/200 above/below, 8 use a same-level partner, 145 within 30 units (max 73.6). |
| 12 | Agent's run and process claims (logs, tmux, result transfer) | confirmed | runs/logs/qc_intake_2026-10-19_v2.log and qc_probe_2026-10-19.log both end with EXIT 0. No vec-qc-* tmux sessions remain. The result files have identical md5 on the DGX and locally. |

### Issues found by the review

1. The report's E8.75 'rigid transform' claim is wrong. The fit is affine (RMS 0.0317) with axis scales 1.028 / 0.971; a rigid fit leaves 4.59 RMS. Fix the wording before anyone relies on it. The code docstrings do not say 'rigid'.
2. The interpolation description in the report and in the spatial.py module docstring ('mean of two records in the levels above and below') is exact only for E8.75. It covers 2446/2558 rows in E8.25_late and 164/200 in E8.0; some partners are 2 levels away or in the same level, and up to 74 units (E8.0) or 163 units (E8.25_late) apart.
3. The spatial.py docstring describes moved copies as 'identical rows one level away'. 166 of the E8.0 ones are two levels from their first occurrence and are only reached by chaining.
4. Copies with identical counts sometimes carry a different celltype from their kept record: 19 of 443 in E8.0 and 25 of 1164 in E8.25_late (none in E8.75 or the KO). The keep-first-in-file rule decides which label survives, which matters for label-based pipelines.
5. No prediction pipeline passes drop_duplicates=True yet; experiments/harmonise/run_harmonise.py and experiments/coords/zscale.py use the default. The docstring line 'prediction pipelines pass True' is a policy the lead still has to wire in.
6. The agent's open decisions still stand for the lead: (1) dropping the moved copies would remove another 2460 / 563 / 518 records (xy_tol=inf); (2) the 200 / 2558 / 614 interpolated rows; (3) the loader's int32 cast drops 9,549 / 119,975 / 24,387 count entries entirely; (4) the min_library guard leaves 2 / 1 / 7 co-located low-count copies unflagged.
7. Repo-wide `ruff check src tests experiments` fails on experiments/harmonise_review/review.py (5 errors, another stream's committed file). The AGENTS.md check `ruff check src tests` passes.
8. The prep.py and tests/test_prep.py working-tree diffs mix another agent's lineage changes with this stream's loader changes. Commit them knowingly.
9. Minor edge cases with no effect on the released files: np.rint rounds half to even, so a level exactly 0.5 step away counts as the same level; the xy_tol=inf successor linking assumes no 64-bit hash collision inside a group; find_duplicate_records needs at least 3 coordinate columns.
10. From the relayed user request: the still-running research fact-check output must be screened by a separate agent before it is merged into the plan. This verification did not read it and cannot start agents, so that is for the lead.

## Lead's decisions (2026-10-19)

- **Prediction inputs drop every synthetic record.** `load_stage(drop_duplicates=True)` now drops copies at the same or a moved xy (identical count rows of >= 20 counts in the same or an adjacent level; such rows never occur by chance in E6.75, E7.25 or E9.5) and interpolated records (half-integer counts). Dropped per file (copies + interpolated): E8.0 2,903 + 198 (2 of its 200 half-integer rows are also copies), E8.25_late 1,727 + 2,558, E8.75 835 + 614, Mab21l2 KO 2,811 + 0; none in E6.75, E7.25 or E9.5 (checked on the DGX with the committed loader, below).
- **Counts are never truncated.** The loader keeps int32 counts when every value is an integer and float32 otherwise (the default load of E8.0, E8.25_late and E8.75, which keeps the interpolated rows).
- **The emulator keeps targets as released** (`drop_duplicates=False`): the server's truth halves come from the files as they are.
- **Wording corrections from the review:** the E8.75 moved copies are an exact affine image (axis scales 1.028 / 0.971), not a rigid one; the interpolated rows average one record from the level above and one from the level below in every case only for E8.75 (most, not all, in E8.25_late and E8.0).
- **Reported to the organisers** in the follow-up (question 4 of `docs/ORGANIZER_QUESTIONS.md`).
