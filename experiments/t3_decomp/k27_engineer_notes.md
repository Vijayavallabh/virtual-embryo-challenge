> Record for ledger row K27 (T3 candidate CHIPQ-S): the rule reviser's notes for the read's engineers (step order, the
> launcher's two invocations, exit codes). Copied from the session's scratch at the row's registration commit; where
> the script or the launcher differ, the log entry of the placement says how.

# K27 (CHIPQ-S): what the read's engineer needs beyond the row

Revised after the cold attack (points_settled.md section D). The row is placed in two commits: the registration
commit, before P1, and the second commit, after P2 and before P3, which fills the listing and the pins. Every counted
step runs at the second commit.

## Files to write (nothing else)

- `experiments/t3_decomp/k27_chip_query.py`, cloned from `k24_wider_tables.py`; it imports `chipq.py` for the
  database part and replaces k18's ARMS, build_arm and fold_factors inside one process only, restored on exit.
- `experiments/t3_decomp/run_k27.sh`, run_k24.sh's steps and exit codes, split into the two invocations below. Trace
  it with `bash -x` on a harmless path before the first real launch.
- `tests/test_k27_chip_query.py` and `tests/test_k27_launcher.py`, synthetic data and hand-made listings only.
- `experiments/t3_ext/TABLES_K27.sha256`: named, and not yet filled, at the registration commit; its 11 lines
  (`<sha256>  <name>`, sorted by name), `K27_LISTING_SHA256` and the 11 table pins in the script are filled by the
  second commit only.
- Unchanged, and byte-identical (test it): `chipq.py`, `chipq_fetch.py`, `outside_forms.py`,
  `k22_outside_forms.py`, `k23_outside_exponents.py`, `k18_built_files.py`, `k21_factor_forms.py`,
  `build_outside_form.py`, `t3x_tables.py`, the frozen T3 files.

## Invocation 1, the database part: `run_k27.sh db <registration commit>` (nothing counted)

1. Gate: HEAD is the registration commit (97), the tracked tree is clean (98), the row naming the script is in
   TRIES.md (95); `TABLES_K27.sha256` is tracked with its 11 lines and the script's 11 pins and listing pin not yet
   filled (refused otherwise, 22); `chipq_kept.yaml`, `CHIPQ_KEPT.sha256` and K23's committed `pregate.json` and
   `reading.json` are tracked; `VEC_T3_EXT_TABLES_K27` is a real directory, not a symlink; `EXTERNAL_DATA.md` holds a
   row naming ChIP-Atlas (X018) (96).
2. The targeted tests, synthetic data only, 8 BLAS threads (10).
3. `chipq_fetch.py verify`: every manifest file re-hashed (2 -> launcher 22).
4. **P1** (`chipq.py p1`): prints two integers, writes nothing. Exit 0 pass, 3 fail (launcher finishes 28, the P1
   stop; nothing counted), 2 refused (launcher 22).
5. **P2** (`k27_chip_query.py p2 --out R`): the practice bars, the pool (at least 10), whether the analysis list holds
   an antigen for the barred factor's gene (recorded in `R/p2.json` and printed as 1 or 0; nothing else follows from
   it: no x2 arm is built either way), the 11 tables into `$VEC_T3_EXT_TABLES_K27`, `R/p2.json`. Exit 0 pass, 3 fail
   (launcher 29, the P2 stop; nothing counted), 2 refused (launcher 22), including a placebo table equal to the
   query's byte for byte or a sha256 that appears twice among the 11 listing lines. It never opens or writes Ctnnb1's
   table.
6. It ends by printing the 11 `<sha256>  <name>.csv` lines and the listing's own sha256, for the second commit.

Launcher exits: 0 (P2 passed, listing lines printed), 10, 22, 28, 29, 95, 96, 97, 98, and 99 when it stopped
unexpectedly. 28 and 29 stop the candidate: no try, no file, no second commit.

## Invocation 2, the counted part: `run_k27.sh run <second commit>`

1. Gate with filled pins: HEAD is the second commit (97), the registration commit is its ancestor, the tracked tree is
   clean (98), the row is in TRIES.md (95); the listing has its 11 lines and its sha256 and the 11 pins are filled in
   the script (22 otherwise); both tables folders (`VEC_T3_EXT_TABLES_K27`, `VEC_T3_EXT_TABLES_TWO`) are real
   directories, not symlinks; X018 is in `EXTERNAL_DATA.md` (96).
2. The targeted tests (10).
3. `chipq_fetch.py verify` (22).
4. The tables checked against the listing: the 11 tables' bytes match the listing, the pins, the manifest and the gene
   order (22 otherwise; nothing counted).
5. **P3** (`pregate --out R`): `R/pregate.json` and its lines, no score; it holds both correlations,
   corr(ln C_Gata4, ln phi_M) and corr(ln C_Gata4, ln two ** 2), each over genes != 1 in either. The eight tries
   count from the start of this step. Exit 0; 5 when it drops the arm on either correlation (launcher 25, the eight
   tries counted); 2 refused (launcher 22, nothing counted).
6. **Guard part** (`run --part guard`, both worlds, ten seeds, into `R/guard/`), then `guard --runs R` ->
   `R/guard.json`; exit 6 reads the run into the results folder as a VOID reading and finishes 26. A rerun after a
   guard failure finishes 26 again with no scoring part started, and the VOID reading is made once, only when the
   results folder holds no `reading.json`.
7. **Scoring parts**, five in parallel as K24: plain g1 seeds 0-4 and 5-9, plain g2 seeds 0-4 and 5-9, and the
   depth-matched world (8 workers by 2 threads). They refuse without a passing `guard.json` of this commit and
   machine. A part's exit code in `R/exit/<part>`, its log in `R/logs/<part>.log`, never printed; unexpected log
   lines counted, never printed (30 if a part did not exit 0).
8. **Read** (`read --runs R --out experiments/t3_decomp/results/k27_chip_query/`), output to /dev/null: 0 a reading
   (VOID included), 43 INCOMPLETE, 44 a reading is there already, 41 a crash in the read.

Launcher exits: 0, 10, 20, 22, 25, 26, 30, 41, 43, 44, 95, 96, 97, 98, and 99 when it stopped unexpectedly.

`NUMPY_MADVISE_HUGEPAGE=0`; `TMPDIR` under the scratch folder `S = ~/v/scratch/k27_<first 7>`; run folder
`R = runs/k27_<first 7>`, resumable and the same folder for both invocations; `--run-dir` to the reference checkout's
gata4 run folder; K23's run folder `runs/k23_<first 7 of K23's commit>` in the same checkout; run detached and wait
with a background poller, never a foreground sleep.
## Exit codes the rule implies

Script: 0 done; 1 a crash or a world whose chip_perm figures are not the pre-gate's (INCOMPLETE; not a reading, and
the eight tries stand counted from P3); 2 refused (no row, dirty tree, pins unfilled or filled at the wrong commit,
listing off its pin, a duplicated placebo table or a sha256 twice in the listing, kept list absent or uncommitted, a
folder variable unset, bytes or records disagreeing, a gene order mismatch, a scaled SD off S by more than a relative
1e-12, the frozen files off their pins, `two` off K23's figures, a scoring part without a passing guard.json); 3 a bar
not met (P1, P2); 4 the reading exists; 5 the pre-gate drops the arm on either correlation; 6 a guard failed.
Launcher: the two lists above.

## What each step prints

- P1: two integers and nothing else. No gene name, no antigen name, no score, no table.
- P2: counts only (the practice kept count and size of B, the pool size, genes changed per table, and 1 or 0 for
  whether the analysis list holds an antigen for the barred factor's gene) and `<sha256>  <name>.csv` lines for the 11
  tables plus the listing's sha256.
- P3: the pre-gate's lines (figures, e per size, a, k, both correlations), no score.
- The guard and scoring parts: progress lines only; no score in any log.
- The read: one of the row's fixed sentences, the reported figures, and the files' sha256.
- Nowhere: a placebo antigen's name, a binding set's gene names, a test score, a leaderboard figure.

## Things easy to get wrong

- The base of every paired difference is `x_e2`, not `barred_form`; `barred_form` is guard 1 and a reported arm.
- Clause (g) is a new check in the read, not in the arms: in every AHEAD cell each control of that size is also
  compared with `x_e2`, and a control AHEAD there disqualifies the size.
- P3 computes two correlations and either one at 0.5 or more, or undefined, drops the arm and its controls at both
  sizes; the second is against `ln two ** 2`, the table the base already carries.
- The twin's `a` and `k` are computed from the whole arm's product (x_e2 times the scaled table), one twin per size.
- `perm` is an index permutation within deciles of each world's basis, seed 20261113, 10 bins, and it permutes the
  SCALED table, so `perm(C ** e) = perm(C) ** e`; use K22's function unchanged.
- The placebo draw is sequential without replacement over seeds 0 to 9 on the name-sorted pool (the row's words); a
  per-seed independent draw is wrong here.
- P2 refuses a placebo table byte-identical to the query's, and any sha256 twice among the 11 listing lines.
- `finish` is called with a zero profile (no rescale). Do not pass a reference profile.
- Nothing is rescaled or clipped after any product.
- Ctnnb1's table is written by no step of this read; the test build after a ruling computes its own `e`.
- Gata4's exponent, the pool indices, `a`/`k` and both correlations go into `pregate.json` to full precision and are
  applied unchanged in both worlds.
- No x2 arm is built and no x2 try is counted, whatever P2 records about an antigen for the barred factor's gene: the
  launcher starts four cells.
- The listing's lines and the script's pins are filled only in the second commit; invocation 1 refuses if they are
  already filled, invocation 2 refuses if they are not.
