#!/usr/bin/env bash
# Row E49 (experiments/harness/TRIES.md), one step per call, on the DGX. Run from the root of a worktree of the campaign
# checkout's repository (e.g. ~/v/scratch/e49-wt, made with `git -C ~/v/vec-campaign worktree add --detach`), so that
# the campaign checkout's HEAD never moves; detached, and wait on the exit file with a background poller, never a
# foreground sleep. The worktree uses the campaign checkout's venv and data (E49_VENV, E49_DATA below); the script
# imports this worktree's own src/ and refuses otherwise. Steps, in order:
#   tests COMMIT     the unit tests of the forms, the builder option and the script (synthetic arrays only)
#   describe COMMIT  step (0), target-free: E7.5 and E8.0 only, the validation board's panel and count rule, build
#                    seeds 1,000,001 to 1,000,010; the forms' shape facts and every task's seed-free F1 blocks (brackets
#                    only); fixes whether F1 and its control are dropped. Not a try. Exit 0 written, 1 an error
#   gate COMMIT      step (1): the floor and fgw_shape0x under build main, all ten seeds, against E47's committed scores
#                    (run at 9967cf1, committed in d5ee9b5) within 0.01 board points per seed and per ten-seed mean;
#                    nine processes, one per task; exit 0 all pass, 3 some task fails (stop and report; nothing is
#                    re-run selectively; nothing is counted), 1 an error
#   run COMMIT       step (2): the arms and controls under both builds, nine processes (one per task, at most 10 at
#                    once), seed by seed; a seed whose file exists is skipped, so a killed run is resumed by calling the
#                    step again; each process refuses before anything loads unless row E49 is in TRIES.md, the tree is
#                    clean, HEAD is COMMIT, describe.json is at COMMIT and all nine gate files pass at COMMIT
#   read COMMIT      the reading into experiments/coords/results/e49_shape_interp/ (reading.json, reading.txt); it
#                    refuses to overwrite either file
# COMMIT is the run commit (full hash), fixed before step (0). Each step refuses on an unclean tree (exit 98) or another
# HEAD (exit 99), runs with 8 BLAS threads a process and NUMPY_MADVISE_HUGEPAGE=0, and writes ~/v/scratch/e49/<STEP>.exit
# (and one exit file per task under runs/e49/exit/). The script prints only "<task> seed <s> done", "wrote <path>" and
# "gate <task> pass|fail"; any other log line (imported code, warnings) is counted, never printed: per task in
# runs/e49/exit/<step>-<task>.other_lines.
set -u
STEP=${1:?usage: run_e49.sh tests|describe|gate|run|read COMMIT}
COMMIT=${2:?the run commit (full hash)}
S=~/v/scratch/e49
mkdir -p "$S/tmp"
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 NUMPY_MADVISE_HUGEPAGE=0 TMPDIR="$S/tmp"
export UV_PROJECT_ENVIRONMENT=${E49_VENV:-$HOME/v/vec-campaign/.venv}
DATA=${E49_DATA:-$HOME/v/vec-campaign/data/raw}
PY="uv run --frozen --no-sync python -I -B"
finish() { echo "$1" > "$S/$STEP.exit"; exit "$1"; }
if [ -n "$(git status --porcelain)" ]; then
  echo "refused: the tree is not clean" >&2
  finish 98
fi
if [ "$(git rev-parse HEAD)" != "$COMMIT" ]; then
  echo "refused: HEAD is not $COMMIT" >&2
  finish 99
fi
git rev-parse HEAD > "$S/$STEP.head"
SCRIPT=experiments/coords/e49_shape_interp_read.py
R=runs/e49
TASKS=(half late new1 new2 e75 e725 x1 x2 x3)
OK_LINE='^([a-z0-9]+ seed [0-9]+ done|wrote [^ ]+|gate [a-z0-9]+ (pass|fail))$'

parallel() {  # $1 = gate|run: one process per task, all nine at once (at most 10)
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
      tests/test_shape_interp.py tests/test_shape_interp_build.py tests/test_e49_shape_interp_read.py \
      tests/test_t2_build.py tests/test_occupancy_frame.py > "$S/tests.log" 2>&1
    finish $? ;;
  describe)
    mkdir -p "$R/logs"
    $PY "$SCRIPT" describe --commit "$COMMIT" --runs "$R" --data "$DATA" > "$R/logs/describe.log" 2>&1
    e=$?
    grep -cvE "$OK_LINE" "$R/logs/describe.log" > "$R/describe.other_lines" 2>/dev/null
    finish $e ;;
  gate)
    parallel gate
    finish $? ;;
  run)
    parallel run
    finish $? ;;
  read)
    $PY "$SCRIPT" read --commit "$COMMIT" --runs "$R" --out experiments/coords/results/e49_shape_interp > /dev/null 2>&1
    finish $? ;;
  *)
    echo "unknown step $STEP" >&2
    finish 97 ;;
esac
