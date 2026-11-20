#!/usr/bin/env bash
# Row H28 (experiments/harness/TRIES.md), one step per call, on a DGX campaign checkout (~/v/vec-campaign, never the
# reference checkout). Run from the checkout's root, detached, and wait on the exit file with a background poller;
# never a foreground sleep. Steps, in order:
#   tests COMMIT   the unit tests of the script and of row H27's (synthetic arrays and hand-made tables only)
#   run COMMIT     the ten parts: per transition five processes of two seeds, the two transitions on E10.5 one after
#                  the other (E8.25_late -> E8.5 is not buildable: one input stage); each part refuses before anything
#                  loads unless row H28 is in TRIES.md, the tree is clean and its output does not exist. The script
#                  prints "<task> seed <s> (<sec> s) done" and "wrote <path>" only; any other log line is counted,
#                  never printed: the count per part is in runs/h28/exit/<short>-<seeds>.other_lines and the total in
#                  runs/h28/exit/other_lines.
#   read COMMIT    the one read: void conditions, gate and fixed statements, into experiments/harness/results/
#                  h28_heart_trend_strength (its two printed lines, the outcome and what follows, in read.out)
# COMMIT is row H28's commit (full hash). Each step refuses on an unclean tree (exit 98) or another HEAD (exit 99),
# runs with 8 BLAS threads a process and NUMPY_MADVISE_HUGEPAGE=0, and writes ~/v/scratch/h28/<STEP>.exit (and, for
# run, one exit file per part under runs/h28/exit/).
set -u
STEP=${1:?usage: run_h28.sh tests|run|read COMMIT}
COMMIT=${2:?the commit row H28 names (full hash)}
S=~/v/scratch/h28
mkdir -p "$S/tmp"
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 NUMPY_MADVISE_HUGEPAGE=0 TMPDIR="$S/tmp"
export UV_PROJECT_ENVIRONMENT="$PWD/.venv"
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
SCRIPT=experiments/harness/h28_heart_trend_strength.py
R=runs/h28
TASKS=("T2h/<=E9.5->E10.5:e95" "T2h/<=E8.75->E10.5:e875")
DONE_LINE='^T2h/[^ ]+ seed [0-9]+ \([0-9.]+ s\) done$'
WROTE_LINE='^wrote [^ ]+/h28\.json$'

case "$STEP" in
  tests)
    $PY -m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp="$S/tmp/pytest" \
      tests/test_h28_heart_trend_strength.py tests/test_h27_heart_forecast_size.py > "$S/tests.log" 2>&1
    finish $? ;;
  run)
    mkdir -p "$R/logs" "$R/exit"
    bad=0
    for entry in "${TASKS[@]}"; do
      key=${entry%:*}
      short=${entry##*:}
      pids=()
      for seeds in 0-1 2-3 4-5 6-7 8-9; do
        (
          $PY "$SCRIPT" run --task "$key" --seeds "$seeds" --out "$R/$short/$seeds" > "$R/logs/$short-$seeds.log" 2>&1
          echo $? > "$R/exit/$short-$seeds"
        ) &
        pids+=($!)
      done
      for p in "${pids[@]}"; do wait "$p"; done
      for seeds in 0-1 2-3 4-5 6-7 8-9; do
        [ "$(cat "$R/exit/$short-$seeds" 2>/dev/null)" = 0 ] || bad=1
        # imported code and Python warnings write to the same log: count the lines other than the two expected kinds
        # and keep only the count, in the exit folder; the lines are never printed
        n=$(grep -cvE "$DONE_LINE|$WROTE_LINE" "$R/logs/$short-$seeds.log" 2>/dev/null)
        echo "${n:-unreadable}" > "$R/exit/$short-$seeds.other_lines"
      done
    done
    cat "$R"/exit/*.other_lines | awk '{ if ($1 ~ /^[0-9]+$/) s += $1; else u = 1 } END { print (u ? "unreadable" : s) }' \
      > "$R/exit/other_lines"
    finish $bad ;;
  read)
    $PY "$SCRIPT" read --runs "$R" --commit "$COMMIT" --out experiments/harness/results/h28_heart_trend_strength \
      > "$S/read.out" 2>&1
    finish $? ;;
  *)
    echo "unknown step $STEP" >&2
    finish 97 ;;
esac
