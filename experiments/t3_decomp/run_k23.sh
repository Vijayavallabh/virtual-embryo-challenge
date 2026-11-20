#!/usr/bin/env bash
# Row K23 (experiments/harness/TRIES.md; experiments/t3_decomp/k23_outside_exponents.py) on a DGX campaign checkout
# (~/v/vec-campaign-c, never the reference checkout), the machine rows K18, K21 and K22 ran on. Run from the checkout's
# root, detached (setsid nohup ... &), and wait on $S/all.exit with a background poller; never a foreground sleep.
#
#   bash experiments/t3_decomp/run_k23.sh COMMIT
#
# COMMIT: the full hash of the row's commit (the row is in TRIES.md at that commit; TABLES.sha256, K22's pins in
# outside_forms.py and K22's committed pregate.json and reading.json are there). The tables folder is fixed below,
# ~/v/vec-ext/t3/tables (combine.py's --out), and reaches the script only through the environment variable
# VEC_T3_EXT_TABLES. S = ~/v/scratch/k23_<first 7>, the run folder R = runs/k23_<first 7> (resumable: a rerun at the
# same commit skips the seeds already written).
# Steps, in order, stopping at the first failure:
#   1. HEAD is COMMIT (exit 97), the tracked tree is clean (98), the row naming the script is in TRIES.md (95),
#      experiments/t3_ext/TABLES.sha256 and K22's committed pregate.json and reading.json are tracked and the tables
#      folder is a real directory, no symlink (96);
#   2. the targeted tests, synthetic data only, 8 BLAS threads (exit 10);
#   3. the pre-gate, target-free, before any scoring: R/pregate.json and $S/pregate.log (no score). Exit 22 when the
#      script refuses (its exit 2: pins, listing, folder, bytes or records disagree, frozen files off their pins;
#      nothing counted), 25 when it drops both exponents (its exit 5: nothing runs, nothing counted), 20 any other
#      failure;
#   4. the five parts in parallel: plain g1 and g2, each as seeds 0-4 and 5-9 (8 BLAS threads a process), and the
#      depth-matched world (8 workers by 2 threads, 16 BLAS threads in its set-up, as rows U3d, K18, K21 and K22); each
#      part's exit code in R/exit/<part>, its log in R/logs/<part>.log (never printed); the log lines other than the
#      expected kinds are counted, never printed: R/exit/<part>.other_lines and the total in R/exit/other_lines (exit
#      30 if a part did not exit 0);
#   5. the read into experiments/t3_decomp/results/k23_outside_exponents/ (output to /dev/null): exit 0 on a reading
#      (whatever it says, VOID included), 43 INCOMPLETE (a seed file missing; rerun this launcher at the same commit),
#      44 a reading is there already (never overwritten), 41 a crash in the read (masked traceback; fix the reading
#      code only, amend the row, read again into the emptied folder).
# all.exit: 0, 10, 20, 22, 25, 30, 41, 43, 44, 95, 96, 97, 98 as above; 99 the launcher stopped unexpectedly (a failing
# command, a TERM, INT or HUP signal, or an exit that was not one of the finishes). NUMPY_MADVISE_HUGEPAGE=0; TMPDIR
# under $S.
set -u
COMMIT=${1:?usage: run_k23.sh COMMIT}
S=~/v/scratch/k23_${COMMIT:0:7}
mkdir -p "$S/tmp" || exit 99
finish() { trap - ERR TERM INT HUP EXIT; echo "$1" > "$S/all.exit"; exit "$1"; }
abort() { trap - ERR TERM INT HUP EXIT; echo 99 > "$S/all.exit"; exit 99; }
trap abort ERR TERM INT HUP EXIT

TABLES=$HOME/v/vec-ext/t3/tables  # combine.py's --out; fixed, no override
export NUMPY_MADVISE_HUGEPAGE=0 TMPDIR="$S/tmp" UV_PROJECT_ENVIRONMENT="$PWD/.venv" VEC_T3_EXT_TABLES="$TABLES"
PY="uv run --frozen --no-sync python -I -B"
threads() { export OMP_NUM_THREADS=$1 OPENBLAS_NUM_THREADS=$1 MKL_NUM_THREADS=$1; }
SCRIPT=experiments/t3_decomp/k23_outside_exponents.py
R=runs/k23_${COMMIT:0:7}
RUN_DIR=$HOME/v/virtual-embryo-challenge/runs/t3_decomp/gata4_20261028T122425Z  # read only; fixed, no override
OUT=experiments/t3_decomp/results/k23_outside_exponents

[ "$(git rev-parse HEAD)" = "$COMMIT" ] || { echo "refused: HEAD is not $COMMIT" >&2; finish 97; }
[ -z "$(git status --porcelain --untracked-files=no)" ] || { echo "refused: the tree is not clean" >&2; finish 98; }
grep -q '^| .*`experiments/t3_decomp/k23_outside_exponents.py`' experiments/harness/TRIES.md \
  || { echo "refused: no TRIES.md row names the script" >&2; finish 95; }
for f in experiments/t3_ext/TABLES.sha256 experiments/t3_decomp/results/k22_outside_forms/pregate.json \
    experiments/t3_decomp/results/k22_outside_forms/reading.json; do
  GIT_OPTIONAL_LOCKS=0 git ls-files --error-unmatch "$f" > /dev/null 2>&1 \
    || { echo "refused: $f is not tracked" >&2; finish 96; }
done
{ [ -d "$TABLES" ] && [ ! -L "$TABLES" ]; } || { echo "refused: $TABLES is not a real directory" >&2; finish 96; }
git rev-parse HEAD > "$S/head"
mkdir -p "$R/exit" "$R/logs"

threads 8
rc=0
$PY -m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp="$S/tmp/pytest" \
  tests/test_k23_outside_exponents.py tests/test_k22_outside_forms.py tests/test_t3_outside_forms.py \
  tests/test_k21_factor_forms.py tests/test_k18_built_files.py > "$S/tests.log" 2>&1 || rc=$?
echo $rc > "$S/tests.exit"
if [ $rc != 0 ]; then finish 10; fi

if [ ! -e "$R/pregate.json" ]; then
  rc=0
  $PY "$SCRIPT" pregate --out "$R" --run-dir "$RUN_DIR" > "$S/pregate.log" 2>&1 || rc=$?
  echo $rc > "$S/pregate.exit"
  if [ $rc = 2 ]; then finish 22; fi
  if [ $rc = 5 ]; then finish 25; fi
  if [ $rc != 0 ]; then finish 20; fi
fi

PARTS=("plain_g1_0-4" "plain_g2_0-4" "plain_g1_5-9" "plain_g2_5-9" "dm")
pids=()
for part in "${PARTS[@]}"; do
  (
    if [ "$part" = dm ]; then
      threads 16
      $PY "$SCRIPT" run --world dm --workers 8 --threads 2 --out "$R" --run-dir "$RUN_DIR" \
        > "$R/logs/$part.log" 2>&1
    else
      threads 8
      IFS=_ read -r _ task seeds <<< "$part"
      $PY "$SCRIPT" run --world plain --task "$task" --seeds "$seeds" --out "$R" --run-dir "$RUN_DIR" \
        > "$R/logs/$part.log" 2>&1
    fi
    echo $? > "$R/exit/$part"
  ) &
  pids+=($!)
done
for p in "${pids[@]}"; do wait "$p" || true; done
bad=0
for part in "${PARTS[@]}"; do
  [ "$(cat "$R/exit/$part" 2>/dev/null)" = 0 ] || bad=1
  # the expected kinds: progress and write lines of the script, K18's per-row timing line and its fixed self-check
  # line; anything else (imported code, warnings, the frozen script's set-up lines) is counted, never printed
  n=$(grep -cvE '^K23 plain g[12] seed [0-9]+: (scored|present, skipped)$|^wrote [^ ]+\.json|^\[rep[12] seed [0-9]+\] K18 depth-matched: [0-9]+ arms and the check, [0-9]+s$|^self-check passed ' \
      "$R/logs/$part.log" 2>/dev/null) || true
  echo "${n:-unreadable}" > "$R/exit/$part.other_lines"
done
cat "$R"/exit/*.other_lines | awk '{ if ($1 ~ /^[0-9]+$/) s += $1; else u = 1 } END { print (u ? "unreadable" : s) }' \
  > "$R/exit/other_lines"
if [ $bad != 0 ]; then finish 30; fi

rc=0
$PY "$SCRIPT" read --runs "$R" --out "$OUT" > /dev/null 2>&1 || rc=$?
echo $rc > "$S/read.exit"
if [ $rc = 3 ]; then finish 43; fi
if [ $rc = 4 ]; then finish 44; fi
if [ $rc != 0 ]; then finish 41; fi
rm -rf "$S/tmp"
finish 0
