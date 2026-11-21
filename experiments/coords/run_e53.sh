#!/usr/bin/env bash
# Row E53 (experiments/harness/TRIES.md), the geodesic overlay's kill test, one step per call, on the DGX. Run from the
# root of a worktree of the campaign checkout's repository (e.g. ~/v/scratch/e53-wt, made with
# `git -C ~/v/vec-campaign worktree add --detach`), so that the campaign checkout's HEAD never moves; detached, and
# wait on the exit file with a background poller, never a foreground sleep. The worktree uses the campaign checkout's
# venv and data (E53_VENV, E53_DATA below); the venv must hold the gpu group (torch). The script imports this
# worktree's own src/ and refuses otherwise. Steps, in order; each runs only if the step before it wrote exit 0 for the
# same COMMIT (else exit 96, refused):
#   tests COMMIT     the unit tests of vec.lddmm, the geodesic forms, the harness arms and the script (synthetic arrays
#                    only; CPU)
#   describe COMMIT  the target-free gate, not a try: ten processes (the test pair and the nine tasks' brackets), the
#                    arm's fit at build seeds 1,000,001 to 1,000,010, one GPU each (E53_GPUS, round robin); then the
#                    pregate applies g1 to g3 to the test pair: exit 0 pass, 4 fail (stop: nothing runs, nothing
#                    counts), 1 an error
#   j7 COMMIT        J7's dry run of the 0.01 tolerance on E48's committed rescore, quoted in the reading; exit 0
#                    written, 1 an error
#   gate COMMIT      the reproduction gate: the floor and C0 against E47, E49c against E49's committed reading, build
#                    main, ten seeds; nine processes, one per task, CPU only; exit 0 all pass, 3 some task fails (stop
#                    and report; nothing is counted), 1 an error
#   run COMMIT       C0, E49c, the arm and the placebo under both builds, nine processes (one per task, one GPU each,
#                    round robin), seed by seed; a seed whose file exists is skipped, so a killed run is resumed by
#                    calling the step again
#   read COMMIT      the reading into experiments/coords/results/e53_geodesic/ (reading.json, reading.txt); refused if
#                    either exists (exit 95), and the script refuses to overwrite them too
# COMMIT is the run commit (full hash): the row's commit. Each step refuses on an unclean tree (exit 98) or another
# HEAD (exit 99), runs with 8 BLAS threads a process and NUMPY_MADVISE_HUGEPAGE=0, and writes
# ~/v/scratch/e53/<STEP>.exit on every path (a call without its two arguments writes usage.exit, 2; an unknown step
# 97) and one exit file per process under runs/e53/exit/. The script prints only step lines; any other log line is
# counted, never printed: per process in runs/e53/exit/<step>-<name>.other_lines. Nothing is removed.
set -u
S="$HOME/v/scratch/e53"
mkdir -p "$S/tmp"
if [ "$#" -lt 2 ] || [ -z "$1" ] || [ -z "$2" ]; then
  echo "usage: run_e53.sh tests|describe|j7|gate|run|read COMMIT" >&2
  echo 2 > "$S/usage.exit"
  exit 2
fi
STEP="$1"
COMMIT="$2"
finish() { echo "$1" > "$S/$STEP.exit"; exit "$1"; }
case "$STEP" in
  tests) PREV="" ;;
  describe) PREV="tests" ;;
  j7) PREV="describe" ;;
  gate) PREV="j7" ;;
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
export UV_PROJECT_ENVIRONMENT="${E53_VENV:-$HOME/v/vec-campaign/.venv}"
DATA="${E53_DATA:-$HOME/v/vec-campaign/data/raw}"
read -r -a GPUS <<< "${E53_GPUS:-2 3 4 5 6 7}"
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
SCRIPT=experiments/coords/e53_geodesic_read.py
R=runs/e53
OUT=experiments/coords/results/e53_geodesic
TASKS=(half late new1 new2 e75 e725 x1 x2 x3)
PAIRS=(test "${TASKS[@]}")
OK_LINE='^([a-z0-9]+ seed [0-9]+ done|wrote [^ ]+|gate [a-z0-9]+ (pass|fail)|describe [a-z0-9]+ done|pregate (pass|fail)|j7 (within|beyond))$'

parallel() {  # $1 = describe|gate|run, $2 = the option (--pair|--task), $3 = gpu|cpu, then the names
  local sub="$1" opt="$2" dev="$3" i=0 pids=() name g
  shift 3
  mkdir -p "$R/logs" "$R/exit"
  for name in "$@"; do
    if [ "$dev" = gpu ]; then g="${GPUS[$((i % ${#GPUS[@]}))]}"; else g=""; fi
    (
      CUDA_VISIBLE_DEVICES="$g" $PY "$SCRIPT" "$sub" "$opt" "$name" --commit "$COMMIT" --runs "$R" --data "$DATA" \
        > "$R/logs/$sub-$name.log" 2>&1
      echo $? > "$R/exit/$sub-$name"
    ) &
    pids+=($!)
    i=$((i + 1))
  done
  for p in "${pids[@]}"; do wait "$p"; done
  local bad=0 fail=0 e n
  for name in "$@"; do
    e=$(cat "$R/exit/$sub-$name" 2>/dev/null)
    if [ "$e" = 3 ]; then fail=1; elif [ "$e" != 0 ]; then bad=1; fi
    n=$(grep -cvE "$OK_LINE" "$R/logs/$sub-$name.log" 2>/dev/null)
    echo "${n:-unreadable}" > "$R/exit/$sub-$name.other_lines"
  done
  if [ $bad = 1 ]; then return 1; fi
  if [ $fail = 1 ]; then return 3; fi
  return 0
}

single() {  # $1 = the script's step: run it once, its log and exit under runs/e53
  mkdir -p "$R/logs" "$R/exit"
  CUDA_VISIBLE_DEVICES="" $PY "$SCRIPT" "$1" --commit "$COMMIT" --runs "$R" --data "$DATA" > "$R/logs/$1.log" 2>&1
  local e=$?
  echo "$e" > "$R/exit/$1"
  local n
  n=$(grep -cvE "$OK_LINE" "$R/logs/$1.log" 2>/dev/null)
  echo "${n:-unreadable}" > "$R/exit/$1.other_lines"
  return "$e"
}

case "$STEP" in
  tests)
    CUDA_VISIBLE_DEVICES="" $PY -m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp="$S/tmp/pytest" \
      tests/test_lddmm.py tests/test_geodesic_overlay.py tests/test_e53_geodesic_read.py \
      tests/test_e49_shape_interp_read.py tests/test_shape_interp.py tests/test_shape_interp_build.py \
      tests/test_e47_floor_exact.py tests/test_t2_build.py tests/test_harness.py > "$S/tests.log" 2>&1
    finish $? ;;
  describe)
    parallel describe --pair gpu "${PAIRS[@]}"
    e=$?
    if [ "$e" != 0 ]; then finish 1; fi
    single pregate
    e=$?
    if [ "$e" = 0 ]; then finish 0; elif [ "$e" = 4 ]; then finish 4; else finish 1; fi ;;
  j7)
    single j7
    finish $? ;;
  gate)
    parallel gate --task cpu "${TASKS[@]}"
    finish $? ;;
  run)
    parallel run --task gpu "${TASKS[@]}"
    finish $? ;;
  read)
    if [ -e "$OUT/reading.json" ] || [ -e "$OUT/reading.txt" ]; then
      echo "refused: the reading exists and is never overwritten" >&2
      finish 95
    fi
    $PY "$SCRIPT" read --commit "$COMMIT" --runs "$R" --out "$OUT" > /dev/null 2>&1
    finish $? ;;
esac
