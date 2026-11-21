#!/usr/bin/env bash
# Row E50 (experiments/harness/TRIES.md), one step per call, on the DGX. Run from the root of a worktree of the campaign
# checkout's repository (e.g. ~/v/scratch/e50-wt, made with `git -C ~/v/vec-campaign worktree add --detach`), so that
# the campaign checkout's HEAD never moves; detached, and wait on the exit file with a background poller, never a
# foreground sleep. The worktree uses the campaign checkout's venv and data (E50_VENV, E50_DATA below); the script
# imports this worktree's own src/ and refuses otherwise. Steps, in order; each runs only if the step before it wrote
# exit 0 for the same COMMIT (else exit 96, refused):
#   tests COMMIT   the unit tests of the construction, the builder option, the harness arm and the script (synthetic
#                  arrays only)
#   gate COMMIT    the floor at 50 and fgw_shape0x and fgw_shape0x_series under build main, all ten seeds, against
#                  E47's committed scores (0.1 board points per ten-seed mean, 0.5 per seed); nine processes, one per
#                  task; exit 0 all pass, 3 some task fails (stop and report; nothing is counted), 1 an error. A task
#                  whose gate file of this commit exists is read again, not re-run
#   run COMMIT     the arm and its controls under both builds, nine processes (one per task), seed by seed; a seed
#                  whose file exists is skipped, so a killed run is resumed by calling the step again
#   read COMMIT    the reading into experiments/coords/results/e50_copy_rows/ (reading.json, reading.txt); it refuses
#                  to overwrite either file
# COMMIT is the run commit (full hash): the row's commit. Each step refuses on an unclean tree (exit 98) or another
# HEAD (exit 99), runs with 8 BLAS threads a process and NUMPY_MADVISE_HUGEPAGE=0, and writes
# ~/v/scratch/e50/<STEP>.exit on every path (a call without its two arguments writes usage.exit, 2; an unknown step 97)
# and one exit file per task under runs/e50/exit/. The script prints only "<task> seed <s> done", "wrote <path>" and
# "gate <task> pass|fail"; any other log line is counted, never printed: per task in
# runs/e50/exit/<step>-<task>.other_lines.
set -u
S="$HOME/v/scratch/e50"
mkdir -p "$S/tmp"
if [ "$#" -lt 2 ] || [ -z "$1" ] || [ -z "$2" ]; then
  echo "usage: run_e50.sh tests|gate|run|read COMMIT" >&2
  echo 2 > "$S/usage.exit"
  exit 2
fi
STEP="$1"
COMMIT="$2"
finish() { echo "$1" > "$S/$STEP.exit"; exit "$1"; }
case "$STEP" in
  tests) PREV="" ;;
  gate) PREV="tests" ;;
  run) PREV="gate" ;;
  read) PREV="run" ;;
  *)
    echo "unknown step $STEP" >&2
    finish 97 ;;
esac
if [ -n "$PREV" ]; then
  if [ "$(cat "$S/$PREV.exit" 2>/dev/null)" != 0 ] || [ "$(cat "$S/$PREV.head" 2>/dev/null)" != "$COMMIT" ]; then
    echo "refused: step $PREV has not exited 0 at $COMMIT" >&2
    finish 96
  fi
fi
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 NUMPY_MADVISE_HUGEPAGE=0 TMPDIR="$S/tmp"
export UV_PROJECT_ENVIRONMENT="${E50_VENV:-$HOME/v/vec-campaign/.venv}"
DATA="${E50_DATA:-$HOME/v/vec-campaign/data/raw}"
PY="uv run --frozen --no-sync python -I -B"
if [ -n "$(git status --porcelain)" ]; then
  echo "refused: the tree is not clean" >&2
  finish 98
fi
if [ "$(git rev-parse HEAD)" != "$COMMIT" ]; then
  echo "refused: HEAD is not $COMMIT" >&2
  finish 99
fi
git rev-parse HEAD > "$S/$STEP.head"
SCRIPT=experiments/coords/e50_copy_rows_read.py
R=runs/e50
TASKS=(half late new1 new2 e75 e725 x1 x2 x3)
OK_LINE='^([a-z0-9]+ seed [0-9]+ done|wrote [^ ]+|gate [a-z0-9]+ (pass|fail))$'

parallel() {  # $1 = gate|run: one process per task, all nine at once
  mkdir -p "$R/logs" "$R/exit"
  local pids=() t
  for t in "${TASKS[@]}"; do
    (
      $PY "$SCRIPT" "$1" --task "$t" --commit "$COMMIT" --runs "$R" --data "$DATA" > "$R/logs/$1-$t.log" 2>&1
      echo $? > "$R/exit/$1-$t"
    ) &
    pids+=($!)
  done
  for p in "${pids[@]}"; do wait "$p"; done
  local bad=0 fail=0 e n
  for t in "${TASKS[@]}"; do
    e=$(cat "$R/exit/$1-$t" 2>/dev/null)
    if [ "$e" = 3 ]; then fail=1; elif [ "$e" != 0 ]; then bad=1; fi
    n=$(grep -cvE "$OK_LINE" "$R/logs/$1-$t.log" 2>/dev/null)
    echo "${n:-unreadable}" > "$R/exit/$1-$t.other_lines"
  done
  if [ $bad = 1 ]; then return 1; fi
  if [ $fail = 1 ]; then return 3; fi
  return 0
}

case "$STEP" in
  tests)
    $PY -m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp="$S/tmp/pytest" \
      tests/test_copy_rows.py tests/test_copy_rows_build.py tests/test_e50_copy_rows_read.py \
      tests/test_e47_floor_exact.py tests/test_t2_build.py tests/test_occupancy_frame.py \
      tests/test_shape_interp_build.py tests/test_harness.py tests/test_hidden_counts.py > "$S/tests.log" 2>&1
    finish $? ;;
  gate)
    parallel gate
    finish $? ;;
  run)
    parallel run
    finish $? ;;
  read)
    $PY "$SCRIPT" read --commit "$COMMIT" --runs "$R" --out experiments/coords/results/e50_copy_rows > /dev/null 2>&1
    finish $? ;;
esac
