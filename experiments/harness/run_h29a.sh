#!/usr/bin/env bash
# Row H29a (experiments/harness/TRIES.md), one step per call, on a DGX campaign checkout (~/v/vec-campaign, never the
# reference checkout). Run from the checkout's root, detached, and wait on the exit file with a background poller;
# never a foreground sleep. Steps, in order; each runs only if the step before it wrote exit 0 for the same COMMIT
# (else exit 96, refused):
#   tests COMMIT   the unit tests of the script and of rows H28's and H27's (synthetic arrays and hand-made tables only)
#   run COMMIT     the ten parts: per transition five processes of two seeds, the two transitions on E10.5 one after
#                  the other (E8.25_late -> E8.5 is not buildable: one input stage). A part whose h29a.json exists is
#                  skipped, so a killed run is resumed by calling the step again (the read voids a part of another
#                  commit). Each part refuses before anything loads unless row H29a is in TRIES.md, the tree is clean,
#                  HEAD is COMMIT and its output does not exist. The script prints "<task> seed <s> (<sec> s) done" and
#                  "wrote <path>" only; any other log line is counted, never printed: the count per part is in
#                  runs/h29a/exit/<short>-<seeds>.other_lines and the total in runs/h29a/exit/other_lines.
#   read COMMIT    the one read: void conditions, gates and fixed statements, into experiments/harness/results/
#                  h29a_heart_trend_forms (its printed lines, one outcome word per candidate, in read.out); refused
#                  (exit 95) if either reading file already exists, so a reading is never overwritten
# COMMIT is row H29a's run commit (full hash). Each step refuses on an unclean tree (exit 98) or another HEAD (exit
# 99), runs with 8 BLAS threads a process and NUMPY_MADVISE_HUGEPAGE=0, and writes ~/v/scratch/h29a/<STEP>.exit on
# every path (a call without its two arguments writes usage.exit, 2; an unknown step 97), and, for run, one exit file
# per part under runs/h29a/exit/. Nothing is removed.
set -u
S="$HOME/v/scratch/h29a"
mkdir -p "$S/tmp"
if [ "$#" -lt 2 ] || [ -z "$1" ] || [ -z "$2" ]; then
  echo "usage: run_h29a.sh tests|run|read COMMIT" >&2
  echo 2 > "$S/usage.exit"
  exit 2
fi
STEP="$1"
COMMIT="$2"
finish() { echo "$1" > "$S/$STEP.exit"; exit "$1"; }
case "$STEP" in
  tests) PREV="" ;;
  run) PREV="tests" ;;
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
export UV_PROJECT_ENVIRONMENT="$PWD/.venv"
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
SCRIPT=experiments/harness/h29a_heart_trend_forms.py
R=runs/h29a
OUT=experiments/harness/results/h29a_heart_trend_forms
TASKS=("T2h/<=E9.5->E10.5:e95" "T2h/<=E8.75->E10.5:e875")
DONE_LINE='^T2h/[^ ]+ seed [0-9]+ \([0-9.]+ s\) done$'
WROTE_LINE='^wrote [^ ]+/h29a\.json$'

case "$STEP" in
  tests)
    $PY -m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp="$S/tmp/pytest" \
      tests/test_h29a_heart_trend_forms.py tests/test_h28_heart_trend_strength.py \
      tests/test_h27_heart_forecast_size.py > "$S/tests.log" 2>&1
    finish $? ;;
  run)
    mkdir -p "$R/logs" "$R/exit"
    bad=0
    for entry in "${TASKS[@]}"; do
      key=${entry%:*}
      short=${entry##*:}
      pids=()
      for seeds in 0-1 2-3 4-5 6-7 8-9; do
        if [ -f "$R/$short/$seeds/h29a.json" ]; then
          continue  # written by an earlier call of this step; its exit file is that call's
        fi
        (
          $PY "$SCRIPT" run --task "$key" --seeds "$seeds" --commit "$COMMIT" --out "$R/$short/$seeds" \
            > "$R/logs/$short-$seeds.log" 2>&1
          echo $? > "$R/exit/$short-$seeds"
        ) &
        pids+=($!)
      done
      for p in ${pids[@]+"${pids[@]}"}; do wait "$p"; done
      for seeds in 0-1 2-3 4-5 6-7 8-9; do
        [ "$(cat "$R/exit/$short-$seeds" 2>/dev/null)" = 0 ] || bad=1
        [ -f "$R/$short/$seeds/h29a.json" ] || bad=1
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
    if [ -e "$OUT/reading.json" ] || [ -e "$OUT/reading.txt" ]; then
      echo "refused: $OUT already holds a reading" >&2
      finish 95
    fi
    $PY "$SCRIPT" read --runs "$R" --commit "$COMMIT" --out "$OUT" > "$S/read.out" 2>&1
    finish $? ;;
esac
