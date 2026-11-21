#!/usr/bin/env bash
# Row H30 (experiments/harness/TRIES.md), one step per call, on a DGX campaign checkout (~/v/vec-campaign, never the
# reference checkout). Run from the checkout's root, detached, and wait on the exit file with a background poller; never
# a foreground sleep. CPU only. Steps, in order; each runs only if the step before it wrote exit 0 for the same ROW
# (and, for read, the same RUN) commit (else exit 96, refused):
#   tests ROW        the unit tests of the script and of rows H28 and H27 (synthetic arrays and hand-made tables only)
#   gate ROW         STEP G, target-free (counts no try): runs/h30/gate/gate.json and gate.txt; exit 0 PASS, 3 FAIL or
#                    VOID (H30 ends: nothing further runs), anything else an error
#   tables ROW       STEP K's two chord tables and their sha256 (runs/h30/tables/); refused unless the gate of ROW says
#                    PASS. Then, by hand, before run: commit runs/h30/gate/gate.json and gate.txt as
#                    experiments/harness/results/h30_heart_late_chord/gate/ and runs/h30/tables/tables_sha256.json as
#                    experiments/harness/results/h30_heart_late_chord/tables_sha256.json, byte for byte, and nothing
#                    else: that commit is RUN.
#   run ROW RUN      the ten parts: per backtest five processes of two seeds (all three forms in each), the two backtests
#                    one after the other; refused unless HEAD is RUN, RUN descends from ROW and differs from it in
#                    exactly those three files (exit 95). The script prints "<task> seed <s> (<sec> s) done" and
#                    "wrote <path>" only; any other log line is counted, never printed (runs/h30/exit/*.other_lines).
#   read ROW RUN     the one read into experiments/harness/results/h30_heart_late_chord (reading.json, reading.txt; the
#                    script refuses to overwrite either); its two printed lines in read.out
# ROW is row H30's commit, RUN the run commit (full hashes). Each step refuses on another HEAD (exit 99) or an unclean
# tree (exit 98), runs with 8 BLAS threads a process and NUMPY_MADVISE_HUGEPAGE=0, and writes ~/v/scratch/h30/<STEP>.exit
# on every path (a call without its arguments writes usage.exit, 2; an unknown step 97). Nothing is ever removed.
set -u
S="$HOME/v/scratch/h30"
mkdir -p "$S/tmp"
if [ "$#" -lt 2 ] || [ -z "$1" ] || [ -z "$2" ]; then
  echo "usage: run_h30.sh tests|gate|tables ROW, or run|read ROW RUN" >&2
  echo 2 > "$S/usage.exit"
  exit 2
fi
STEP="$1"
ROW="$2"
RUN="${3:-}"
finish() { echo "$1" > "$S/$STEP.exit"; exit "$1"; }
case "$STEP" in
  tests) PREV="" ;;
  gate) PREV="tests" ;;
  tables) PREV="gate" ;;
  run) PREV="tables" ;;
  read) PREV="run" ;;
  *)
    echo "unknown step $STEP" >&2
    finish 97 ;;
esac
case "$STEP" in
  run|read)
    if [ -z "$RUN" ]; then
      echo "usage: run_h30.sh $STEP ROW RUN" >&2
      finish 2
    fi
    HEAD_WANT="$RUN" ;;
  *) HEAD_WANT="$ROW" ;;
esac
if [ -n "$PREV" ]; then
  if [ "$(cat "$S/$PREV.exit" 2>/dev/null)" != 0 ] || [ "$(cat "$S/$PREV.head" 2>/dev/null)" != "$ROW" ]; then
    echo "refused: step $PREV has not exited 0 at $ROW" >&2
    finish 96
  fi
  if [ "$STEP" = read ] && [ "$(cat "$S/run.run" 2>/dev/null)" != "$RUN" ]; then
    echo "refused: step run has not exited 0 at $RUN" >&2
    finish 96
  fi
fi
if [ "$(git rev-parse HEAD)" != "$HEAD_WANT" ]; then
  echo "refused: HEAD is not $HEAD_WANT" >&2
  finish 99
fi
if [ -n "$(git status --porcelain)" ]; then
  echo "refused: the tree is not clean" >&2
  finish 98
fi
H30R=experiments/harness/results/h30_heart_late_chord
if [ -n "$RUN" ]; then
  ALLOWED="$H30R/gate/gate.json $H30R/gate/gate.txt $H30R/tables_sha256.json"
  CHANGED="$(git diff --name-only "$ROW" "$RUN" | sort | tr '\n' ' ')"
  WANT="$(printf '%s\n' $ALLOWED | sort | tr '\n' ' ')"
  if [ "$RUN" = "$ROW" ] || ! git merge-base --is-ancestor "$ROW" "$RUN" || [ "$CHANGED" != "$WANT" ]; then
    echo "refused: RUN must descend from ROW and differ from it in the gate reading and the table checksums only" >&2
    finish 95
  fi
  if ! cmp -s runs/h30/gate/gate.json "$H30R/gate/gate.json" \
      || ! cmp -s runs/h30/tables/tables_sha256.json "$H30R/tables_sha256.json"; then
    echo "refused: the committed gate reading or table checksums are not the steps' own files" >&2
    finish 95
  fi
fi
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 NUMPY_MADVISE_HUGEPAGE=0 TMPDIR="$S/tmp"
export UV_PROJECT_ENVIRONMENT="$PWD/.venv"
PY="uv run --frozen --no-sync python -I -B"
echo "$ROW" > "$S/$STEP.head"
if [ -n "$RUN" ]; then echo "$RUN" > "$S/$STEP.run"; fi
SCRIPT=experiments/harness/h30_heart_late_chord.py
R=runs/h30
TASKS=("T2h/<=E9.5->E10.5:e95" "T2h/<=E8.75->E10.5:e875")
DONE_LINE='^T2h/[^ ]+ seed [0-9]+ \([0-9.]+ s\) done$'
WROTE_LINE='^wrote [^ ]+/h30\.json$'

case "$STEP" in
  tests)
    $PY -m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp="$S/tmp/pytest" \
      tests/test_h30_heart_late_chord.py tests/test_h28_heart_trend_strength.py \
      tests/test_h27_heart_forecast_size.py > "$S/tests.log" 2>&1
    finish $? ;;
  gate)
    mkdir -p "$R"
    $PY "$SCRIPT" gate --out "$R/gate" > "$S/gate.log" 2>&1
    finish $? ;;
  tables)
    $PY "$SCRIPT" tables --gate "$R/gate" --out "$R/tables" > "$S/tables.log" 2>&1
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
          $PY "$SCRIPT" run --task "$key" --seeds "$seeds" --tables "$R/tables" --out "$R/$short/$seeds" \
            > "$R/logs/$short-$seeds.log" 2>&1
          echo $? > "$R/exit/$short-$seeds"
        ) &
        pids+=($!)
      done
      for p in "${pids[@]}"; do wait "$p"; done
      for seeds in 0-1 2-3 4-5 6-7 8-9; do
        [ "$(cat "$R/exit/$short-$seeds" 2>/dev/null)" = 0 ] || bad=1
        n=$(grep -cvE "$DONE_LINE|$WROTE_LINE" "$R/logs/$short-$seeds.log" 2>/dev/null)
        echo "${n:-unreadable}" > "$R/exit/$short-$seeds.other_lines"
      done
    done
    finish $bad ;;
  read)
    $PY "$SCRIPT" read --runs "$R" --commit "$RUN" --out "$H30R" > "$S/read.out" 2>&1
    finish $? ;;
esac
