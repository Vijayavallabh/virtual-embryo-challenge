# feat-020: one shared core, `vec.core` (a pure refactor)

_Returned as text by the implementing agent on 2026-10-20; the harness does not let agents write report files. Base:
eb946ac, the commit before the move. The changes are uncommitted in the working tree; the lead commits them._

_Renamed after the checks: the module was written as `src/vec/operator.py` and renamed to `src/vec/core.py` (with
`tests/test_core.py` and `experiments/core/`) after an independent verifier found that the name shadowed the standard
library's `operator` when a file in `src/vec/` is run directly (`python src/vec/cli.py` failed). Paths below use the new
names; the sha256 values in the tables are of the files as checked under the old name. The rename and two docstring
corrections changed no function: the tests pass under the new name. The verifier (PASS) also corrected two statements
here: the core holds copies of the T3 script's functions (written there first), not the reverse; and the DGX run of the
end-to-end experiment tests on the old copy was not logged (locally the old code gives 58 passed, 2 skipped)._

## Summary
The primitives that the three task builders each carried now live in one module, `src/vec/core.py`. It is the
"transport-then-adjust" core of docs/PLAN.md section 5, which makes the Generality Award description (Rules s.14: one
shared architecture, with task-specific input and prediction heads) true in the code itself:
- The T1 and T2 heads call the core for allocation, draws and count-space edits.
- The pre-registered T3 script keeps its own primitives (written there first); the core holds identical copies,
  and the tests prove them identical.

No output changed:
- **Released data (DGX).** The T1 and T2 builders' X, obs and obsm are bit-identical before and after the move. So are
  every other part hashed: var, layers, reports, the sampled rows and coordinates, format checks and cloud guards.
- **Golden cases.** 12 synthetic cases, recorded before any source change, still pass on both machines: locally (AVX2)
  and on the DGX (AVX-512).

## What moved where
| Core function (`vec.core`) | Came from | Step | Called by |
|---|---|---|---|
| `quotas` | `vec.t1_scope.quotas` (moved) | allocation | T1 `reweight` |
| `alloc` | copy of `t3_decompose.alloc` | allocation | T3, as its frozen copy |
| `draw_least_used` | `vec.t1_scope._draw` (moved) | draw | T1 `reweight` |
| `draw_balanced` | `vec.t2_build._balanced_rows` (moved) | draw | T2 `sample_cells`, `cloud_guard` |
| `draw_then_repeat` | copy of `t3_decompose.draw_rows` | draw | T3, as its frozen copy |
| `counts_of` | `vec.t1_scope.counts_of` (moved) | count space | T1 `momentum_shift` |
| `mean_cp10k` | `vec.t1_scope.mean_cp10k` (moved) | count space | T1 `momentum_shift` |
| `apply_factors` | `vec.t1_scope.apply_factors` (moved) | count space | T1 `momentum_shift` |
| `apply_gene_factor` | copy of `t3_decompose.apply_gene_factor` | count space | T3, as its frozen copy |
| `rescale_depth`, `thin` | copies of `t3_decompose.rescale_depth`, `thin` | count space | T3, as its frozen copies |
| `log1p_cp10k` | `vec.prep.log1p_cp10k` (re-exported, not moved) | log1p CP10k | T1 `momentum_shift` |

- **Moved functions.** Their full source is unchanged except the name, docstrings included (`check_equivalence.py
  source`, below). The old names stay importable as the same objects: `t1_scope.quotas`, `_draw`, `counts_of`,
  `mean_cp10k`, `apply_factors` and `t2_build._balanced_rows`. A namespace diff against eb946ac shows that neither
  module lost a name: `t1_scope` gains `draw_least_used` and `t2_build` gains `draw_balanced`. `t1_scope` still has
  `recover_counts` and `TARGET_SUM`, which it no longer uses.
- **Head call sites.** `reweight`, `sample_cells` and `cloud_guard` call the core by its new names.
- **T3 copies.** The def line is identical except the rename of `draw_rows` to `draw_then_repeat`, and the body is
  identical. Only the docstrings differ: the core's docstrings name the frozen copy.
- **Imports.** `operator.py` imports numpy, scipy.sparse, anndata and `vec.prep` (`TARGET_SUM`, `log1p_cp10k`,
  `recover_counts`). It imports neither head and nothing under experiments/. A test checks this in the AST, and checks
  that importing the core alone loads neither head.
- **Module docstring.** It lists:
  - the steps implemented now;
  - the steps PLAN section 5 plans: coupling and displacement interpolation (feat-010, feat-015), births (feat-014),
    and leave-one-transition-out calibration (feat-008);
  - the task heads and the boundary between them and the core;
  - the frozen T3 copies.

## What stayed in each head
- **T1 (`vec.t1_scope`):**
  - the chamber config: `ChamberSpec`, `load_spec`, `check_spec`, `classify`;
  - the target shares: `load_shares`, `share_of`, `with_share`, `composition`;
  - the re-weighting plan: a pool per category, `within` and `fallback`, own pools before stand-ins, the final row
    permutation, `_take` and the report;
  - the momentum factors: `_type_rows`, the damped log ratios of mean CP10k with `eps`, clipping, held-flat genes, the
    stand-in skip and `ShiftReport`;
  - `build_submission` (board gene order, cell limits, report) and `write_checked`.
- **T2 (`vec.t2_build`):**
  - the count rule, `n_cells`;
  - `sample_cells`: seed streams, duplicate detection, the jitter (`_ball_shifts`) and the record;
  - the coordinate head: `centre`, `check_rotation`, `proper_rotation`, `rescale_to_rms`, `apply_size_procedure`;
  - `cloud_guard` (its reference draws call `draw_balanced`) and `build_submission`.
- **T3 (`experiments/t3_decomp/`):** all of it, untouched (it is hash-pinned), including the frozen copies and its own
  `lognorm`.

## Other changes
- The docstrings of `t1_scope` and `t2_build` each gain one paragraph naming the core and the head's part.
- In `experiments/t1_scope/pseudo_task.py`, `PROCEDURE_FILES` gains `src/vec/core.py`. The pseudo-task records the
  sha256 of its procedure files, and refuses to write results when they change during a run. Without this entry, the
  code that now draws and edits T1 cells would sit outside both. No builder output depends on this line; revert it if
  you prefer.
- New files: `src/vec/core.py`, `tests/test_core.py`, `experiments/core/check_equivalence.py`, and this
  report.
- Not touched by this feature:
  - `experiments/t3_decomp/`, `configs/`, `docs/`, `src/vec/prep.py`, `src/vec/spatial.py`;
  - `AGENTS.md`, `METHOD_SUMMARY.md`, `feature_list.json`, `progress.md`, `ASSUMPTIONS.md`.

  Some of these (`docs/`, `ASSUMPTIONS.md`, `feature_list.json`, `progress.md`) show as modified in the working tree
  because of the lead's concurrent edits. The DGX "new" copy, which was rsynced from the tree, carries those edits too.
  None of them is code.

## The optional single entry point: not added
I did not add a single core entry point for the T1 and T2 builders to call, for three reasons:
- **The draw policies have different contracts.** `draw_least_used` tracks use counts across a category's pools and
  reports repeats. `draw_balanced` returns positions in a pool. `draw_then_repeat` is frozen in the T3 script.
- **It would add a layer without changing what runs.** One entry point would only be a dispatcher over these
  functions.
- **T3 could not call it.** Its script is hash-pinned.

Instead, the heads call the core's functions by their core names, and a test checks this in the AST.

## Equivalence evidence
### Tests: tests/test_core.py (31 tests)
**1. Golden outputs:** `test_golden_outputs_unchanged` (12 cases) and `test_every_golden_case_has_recorded_hashes`.
- **How they were recorded.** At eb946ac, with `uv run python tests/test_core.py`, before any source change. The
  hashes are hard-coded in the test.
- **T1 cases:**
  - `reweight` in two configurations: restricted labels; and stand-ins sharing a pool with another category, with a
    sequence seed;
  - `build_submission` without a shift and with a shift;
  - `build_submission` with options: stand-ins at the reference end of the shift, `within`, `standin_momentum`,
    `keep_depth=False`, clipped factors, `eps` and `min_cells`.
- **T2 cases:**
  - `sample_cells`: below the pool; above it, with jitter; with duplicate positions in the pool; from an AnnData;
  - `cloud_guard`: a random sample; a contiguous block in its own frame; above the pool;
  - `build_submission`: an AnnData with reordered and extra genes; a matrix with duplicate coordinates allowed.
- **The moved primitives, called by their old names:**
  - `quotas`, with ties;
  - `_draw`, with the generator state afterwards;
  - `_balanced_rows`, `counts_of` and `mean_cp10k`;
  - `apply_factors`, with `keep_depth` both ways and restored cells.
- **How outputs are hashed.** Arrays (X, obs, var, obsm, layers, rows, coordinates, counts) are hashed bit for bit.
  Report dicts are hashed with their floats at 9 significant digits. The reason: numpy's log and exp loops differ in
  the last bit between the AVX2 machine here and the AVX-512 DGX, and the golden test has to pass on both.
- **Config.** The synthetic chamber config belongs to the test, so edits to `configs/` cannot move the hashes.
- **Results.** The tests pass with the code from before the move: locally, and on the DGX, where the final test file
  ran with eb946ac's `vec` (`-k golden`: 13 passed, 18 deselected). They also pass with the code after the move, both
  locally and on the DGX. The table printed after the move is byte-identical to the one printed before it.

**2. T3 equivalence:**
- `test_t3_copies_have_the_frozen_bodies` covers 5 pairs. For each pair, `inspect.getsource` without the def line and
  docstring is identical and the signatures are identical. Every global a body reads is the same object in both
  modules (numpy), or the proven copy `thin`. The body check holds for all five.
- `test_t3_copies_match_on_random_inputs` runs 300 random inputs per primitive. For each input it checks that arrays,
  dtypes and shapes are identical, that the arguments are identical afterwards, and that the generators end in the
  same state, so the same draws were consumed.
- `test_t3_copies_match_on_edge_cases`:
  - empty cells and no cells;
  - ratio = 1, per cell and as a scalar; phi = 1;
  - k > pool, k = 0, k < 0, and an empty pool;
  - ties and zero shares in the allocation, and shares that do not sum to 1.
- `test_t3_lognorm_is_the_core_log1p_cp10k_within_one_float32_ulp`.

**3. Old names and boundaries:**
- `test_old_names_are_the_core_objects` (6 cases);
- `test_heads_call_the_core_and_define_none_of_it`: `reweight` calls `quotas` and `draw_least_used`; `momentum_shift`
  calls `counts_of`, `mean_cp10k`, `apply_factors` and `log1p_cp10k`; `sample_cells` and `cloud_guard` call
  `draw_balanced`;
- `test_core_api_and_imports`.

**4. Literal guard:**
- `test_core_holds_no_stage_constants` checks numeric literals against the allow-list {0, 1, 5, 1e-9, 1e-6, 1e-2}:
  indices, the factor of no change, the number of genes an error message lists, and tolerances.
- `test_core_names_no_class_stage_or_class_list` checks two things:
  - no string outside the docstrings names a released label, a chamber category, a lineage group or a stage;
  - there is no module-level container of strings except `__all__`.

### Suites
| Where | Code | Command | Result |
|---|---|---|---|
| local | working tree | `./init.sh` | lock ok, ruff ok, 320 passed (31 of them new) |
| local | working tree | `uv run pytest experiments/t1_scope experiments/t2_build experiments/t3_decomp -q -p no:cacheprovider` | 58 passed, 2 skipped (the `RUN_T3_E2E` end-to-end tests) |
| DGX | new copy | `pytest tests` | 320 passed |
| DGX | new copy | `RUN_T3_E2E=1 pytest experiments/t1_scope experiments/t2_build experiments/t3_decomp` | 60 passed, 2 warnings |
| DGX | old copy | same | not logged (`tests.log` holds only the new copy's run); locally, eb946ac: 58 passed, 2 skipped |
| DGX | old `vec`, final test file | `pytest tests/test_core.py -k golden` | 13 passed, 18 deselected |

The 2 warnings are a multiprocessing `DeprecationWarning` (fork in a multi-threaded process), raised in the T3
end-to-end test before and after the change.

The source comparison `uv run python experiments/core/check_equivalence.py source --base eb946ac` reports ALL
IDENTICAL:
- for the six moved functions against their eb946ac source;
- for the five T3 primitives against `t3_decompose.py`.

### Released data on the DGX: old (eb946ac) vs new (working tree)
Each run used the mirror's interpreter with `PYTHONPATH=<copy>/src`. Each printed the `vec` it imported:
`~/v/tmp-operator-check/old/src/vec/__init__.py` (no `vec.core`) and
`~/v/tmp-operator-check/new/src/vec/__init__.py`.

The builder files' sha256 values in the run JSONs differ as expected:

| File | Old | New |
|---|---|---|
| `t1_scope.py` | 74596a62… | 8478211f… |
| `t2_build.py` | 5a11ab1b… | db7a517e… |
| `operator.py` | absent | 92bb0476… |

`prep.py`, `io.py`, `boards.py` and `spatial.py` are identical between the two copies. The old and new outputs, by the
first 16 hex digits of each sha256 (equal in both runs):

| Check | X | obs | obsm | Other parts compared | Old = new |
|---|---|---|---|---|---|
| T1 build_submission, shift, seed 0 | `db24d687a6b76799` | `fc856e1a19df268b` | `44136fa355b3678a` | var, layers, report | yes |
| T1 build_submission, shift, seed 1 | `3f13e9a55babd9d5` | `fc856e1a19df268b` | `44136fa355b3678a` | var, layers, report | yes |
| T1 build_submission, no shift, seed 0 | `88c88f70e0dd4740` | `fc856e1a19df268b` | `44136fa355b3678a` | var, layers, report | yes |
| T1 build_submission, no shift, seed 1 | `b4bcd01a8dd82870` | `fc856e1a19df268b` | `44136fa355b3678a` | var, layers, report | yes |
| T2 T2:heart:val_extrap (E9.5), seed 0 | `27a4b65ff9f593e7` | `3a09f1f46c83c797` | `2ae4dca375a8aa3a` | rows, coords, sample record, var, layers, check, guard | yes |
| T2 T2:heart:val_extrap (E9.5), seed 1 | `b018ac6322b9c543` | `3a09f1f46c83c797` | `002119e23582549c` | rows, coords, sample record, var, layers, check, guard | yes |
| T2 T2:heart:val_interp (E8.75), seed 0 | `80d90b54ab06bc8c` | `a0383947aed64e63` | `a481ed365c93c0c8` | rows, coords, sample record, var, layers, check, guard | yes |
| T2 T2:heart:val_interp (E8.75), seed 1 | `cf4978af2474fc20` | `a0383947aed64e63` | `b0606696a88c6b83` | rows, coords, sample record, var, layers, check, guard | yes |

`compare` reports ALL IDENTICAL over 64 hashed parts, the count records and the T1 reuse and shift summaries. The
reports and the sampling records are compared with exact floats. After a final re-sync (the only code change was
`check_equivalence.py`), the check ran again with the final script: ALL IDENTICAL again, and every hash equals the
first run's. The first run is kept in `first_run/`.

The T1 setup:
- reference E9.5_RNA (17,057 cells), start E8.5_RNA (16,787 cells);
- board T1:val (32,285 genes), 5,000 cells, damping 0.5, the file stems as `ref_source` and `start_source`, and the
  pre-registered held-flat genes;
- shares: the reference's own category shares with the chamber categories weighted x3 and epicardium x8. They are an
  arbitrary input for this check, not a prediction.
- In every run epicardium's quota exceeded its pool, so cells repeat; with the shift, six types moved.

The T2 setup:
- pools from `load_stage(..., drop_duplicates=True)` on the board's panel: E9.5 for T2:heart:val_extrap (53,742 cells,
  0 records dropped) and E8.75 for T2:heart:val_interp (23,377 cells, 1,449 dropped);
- n = max_cells: 25,179 and 17,616;
- the format check and the cloud guard were ok in every run.

Both counts are below the pool sizes, so only the golden tests cover drawing above the pool.

## Commands
```bash
# local, before any source change: record the golden table (then pasted into tests/test_core.py)
uv run python tests/test_core.py
# local, after the change
./init.sh
uv run pytest experiments/t1_scope experiments/t2_build experiments/t3_decomp -q -p no:cacheprovider
uv run python experiments/core/check_equivalence.py source --base eb946ac

# DGX scratch copies (the mirror ~/v/virtual-embryo-challenge is never modified or synced; its data/ is read only)
ssh -o BatchMode=yes PrakashDGX_H2 'mkdir -p ~/v/tmp-operator-check/old ~/v/tmp-operator-check/new'
git archive eb946acc1ee53109aef4ac3ec07c005ac34a2cd2 | ssh -o BatchMode=yes PrakashDGX_H2 'tar -x -C ~/v/tmp-operator-check/old'
git ls-files -co --exclude-standard -z > files.lst   # the working tree without ignored files (.venv, data, runs, caches)
rsync -a --from0 --files-from=files.lst --exclude .venv --exclude data --exclude runs ./ PrakashDGX_H2:v/tmp-operator-check/new/

# DGX (tmux): ~/v/tmp-operator-check/run_check.sh, which does
cd ~/v/tmp-operator-check
export OPENBLAS_NUM_THREADS=8 OMP_NUM_THREADS=8 MKL_NUM_THREADS=8
PY=~/v/virtual-embryo-challenge/.venv/bin/python
for v in old new; do
  PYTHONPATH=~/v/tmp-operator-check/$v/src $PY new/experiments/core/check_equivalence.py run \
      --data ~/v/virtual-embryo-challenge/data/raw --out $v.json
done
$PY new/experiments/core/check_equivalence.py compare old.json new.json
# and ~/v/tmp-operator-check/run_tests.sh: the three DGX pytest rows of the table above
```
The results are in `~/v/tmp-operator-check/` on the DGX (44 MB), left in place:
- `old.json` and `new.json` (the final re-run, full hashes), and `first_run/` (the first run, identical);
- `check.log`, `tests.log`, `run_check.sh`, `run_tests.sh`;
- the two copies, `old/` and `new/`.

## What is not identical, and risks
1. **File hashes.** `t1_scope.py` and `t2_build.py` change by construction (sha256 above). Records that pin these files
   by hash describe the old files. For example, both pseudo-task `run.json` files under `experiments/t1_scope/results/`
   record `t1_scope.py` as 74596a62…. The new code reproduces the old outputs bit for bit, but a check that compares a
   recorded hash with the file on disk will see a difference.
2. **`lognorm`.** The T3 script's `lognorm` is not a copy of the core's log1p CP10k. It computes `1e4·c / total`,
   while `vec.prep` computes `c · (1e4 / total)`, so the two agree within one float32 ulp (tested), not bit for bit. T3
   keeps its own, as pre-registered.
3. **T3 at run time.** The T3 head runs its own primitives rather than importing `vec.core`, because its script is
   hash-pinned. The method description for the Generality Award says so: the core holds identical copies of the T3
   script's primitives, proven identical by `tests/test_core.py`.
4. **Golden report floats.** The golden test compares report floats at 9 significant digits (see above). The DGX old
   vs new comparison is exact everywhere, reports included.
5. **Imports.** Importing `t2_build` now also loads `vec.prep`, through `vec.core`, and with it `vec.spatial`,
   pandas and yaml. Nothing `t2_build` computes changes.
6. **Module name (fixed).** The module was first named `vec.operator`, which shadowed the standard library's
   `operator` whenever a file in `src/vec/` ran as a script: `python src/vec/cli.py --help` failed (found by the
   verifier). It is now `vec.core`.
7. **Maintenance.** A later change that is meant to alter a T1 or T2 output must update the golden table in the same
   commit and say why.
