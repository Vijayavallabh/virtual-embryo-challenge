#!/usr/bin/env bash
# Row E54 (experiments/harness/TRIES.md), the kill test of the multi-stage flow, one step per call, on the DGX. Run from
# the root of a worktree of the campaign checkout's repository (e.g. ~/v/scratch/e54-wt, made with
# `git -C ~/v/vec-campaign worktree add --detach`), so that the campaign checkout's HEAD never moves; detached, and wait
# on the exit file with a background poller, never a foreground sleep. The worktree uses the campaign checkout's venv
# (with the uv gpu group synced: torch, torchdiffeq) and data (E54_VENV, E54_DATA below); the script imports this
# worktree's own src/ and refuses otherwise. Steps, in order; each runs only if the step before it wrote exit 0 for the
# same COMMIT (else exit 96, refused):
#   tests COMMIT   the unit tests of vec.t2_mmflow, the harness arm, the script and this launcher (synthetic arrays only;
#                  exit 94 if torch or torchdiffeq cannot be imported, so the flow's tests cannot be skipped)
#   fit COMMIT     the eight fits (target-free: each loads only its training stages), one process per fold, spread over
#                  the GPUs listed in CUDA_VISIBLE_DEVICES (exit 95 when it is unset or empty); a fold whose done.json
#                  is of this commit is not refitted
#   repro COMMIT   the reproduction gate, not counted: the floor and C0 under build main, ten seeds, against E47 within
#                  0.01; nine processes; exit 0 all pass, 3 some task misses (stop and report), 1 an error
#   gate COMMIT    E54a, counted: the leave-one-stage-out gate; writes its reading (gate.json, gate.txt, never
#                  overwritten); exit 0 it passes, 3 it fails (E54b and E54c then do not run and are not counted)
#   run COMMIT     E54b and E54c with C0 and the floor, nine processes, seed by seed; a seed whose file exists is
#                  skipped, so a killed run is resumed by calling the step again
#   read COMMIT    the kill test's reading into experiments/coords/results/e54_mmflow/ (reading.json, reading.txt); it
#                  refuses to overwrite either file
# COMMIT is the run commit (full hash): the row's commit. Each step refuses on an unclean tree (exit 98) or another HEAD
# (exit 99), runs with 8 BLAS threads a process and NUMPY_MADVISE_HUGEPAGE=0, and writes ~/v/scratch/e54/<STEP>.exit
# on every path (a call without its two arguments writes usage.exit, 2; an unknown step 97) and one exit file per fold
# or task under runs/e54/exit/. The script prints only "fit <fold> done", "repro <task> pass|fail", "gate pass|fail",
# "<task> seed <s> done" and "wrote <path>"; any other log line is counted, never printed: per process in
# runs/e54/exit/<step>-<name>.other_lines. Nothing is removed by this launcher.
set -u
S="$HOME/v/scratch/e54"
mkdir -p "$S/tmp"
if [ "$#" -lt 2 ] || [ -z "$1" ] || [ -z "$2" ]; then
  echo "usage: run_e54.sh tests|fit|repro|gate|run|read COMMIT" >&2
  echo 2 > "$S/usage.exit"
  exit 2
fi
STEP="$1"
COMMIT="$2"
finish() { echo "$1" > "$S/$STEP.exit"; exit "$1"; }
case "$STEP" in
  tests) PREV="" ;;
  fit) PREV="tests" ;;
  repro) PREV="fit" ;;
  gate) PREV="repro" ;;
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
export CUBLAS_WORKSPACE_CONFIG=:4096:8
export UV_PROJECT_ENVIRONMENT="${E54_VENV:-$HOME/v/vec-campaign/.venv}"
DATA="${E54_DATA:-$HOME/v/vec-campaign/data/raw}"
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
SCRIPT=experiments/coords/e54_mmflow_read.py
R=runs/e54
OUT=experiments/coords/results/e54_mmflow
TASKS=(half late new1 new2 e75 e725 x1 x2 x3)
FOLDS=(x80_825_late x80 x80_825_late_85 x825_late_85 x85 x75 x725 x825_late)
OK_LINE='^(fit [a-z0-9_]+ done|repro [a-z0-9]+ (pass|fail)|gate (pass|fail)|[a-z0-9]+ seed [0-9]+ done|wrote [^ ]+)$'

parallel() {  # $1 = fit|repro|run; one process per fold (fit) or task, all at once
  mkdir -p "$R/logs" "$R/exit"
  local pids=() x i=0 flag items gpus=()
  if [ "$1" = fit ]; then
    flag=--fold
    items=("${FOLDS[@]}")
    IFS=, read -r -a gpus <<< "$CUDA_VISIBLE_DEVICES"
  else
    flag=--task
    items=("${TASKS[@]}")
  fi
  for x in "${items[@]}"; do
    (
      if [ "$1" = fit ]; then export CUDA_VISIBLE_DEVICES="${gpus[$((i % ${#gpus[@]}))]}"; fi
      $PY "$SCRIPT" "$1" "$flag" "$x" --commit "$COMMIT" --runs "$R" --data "$DATA" --out "$OUT" \
        > "$R/logs/$1-$x.log" 2>&1
      echo $? > "$R/exit/$1-$x"
    ) &
    pids+=($!)
    i=$((i + 1))
  done
  for p in "${pids[@]}"; do wait "$p"; done
  local bad=0 fail=0 e n
  for x in "${items[@]}"; do
    e=$(cat "$R/exit/$1-$x" 2>/dev/null)
    if [ "$e" = 3 ]; then fail=1; elif [ "$e" != 0 ]; then bad=1; fi
    n=$(grep -cvE "$OK_LINE" "$R/logs/$1-$x.log" 2>/dev/null)
    echo "${n:-unreadable}" > "$R/exit/$1-$x.other_lines"
  done
  if [ $bad = 1 ]; then return 1; fi
  if [ $fail = 1 ]; then return 3; fi
  return 0
}

case "$STEP" in
  tests)
    if ! $PY -c "import torch, torchdiffeq" > /dev/null 2>&1; then
      echo "refused: torch or torchdiffeq cannot be imported (sync the uv gpu group)" >&2
      finish 94
    fi
    $PY -m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp="$S/tmp/pytest" \
      tests/test_t2_mmflow.py tests/test_e54_mmflow_read.py tests/test_e54_launcher.py \
      tests/test_e49_shape_interp_read.py tests/test_e47_floor_exact.py tests/test_harness.py \
      tests/test_t2_build.py > "$S/tests.log" 2>&1
    finish $? ;;
  fit)
    if [ -z "${CUDA_VISIBLE_DEVICES:-}" ]; then
      echo "refused: CUDA_VISIBLE_DEVICES names no GPU for the fits" >&2
      finish 95
    fi
    parallel fit
    finish $? ;;
  repro)
    parallel repro
    finish $? ;;
  gate)
    mkdir -p "$R/logs" "$R/exit"
    $PY "$SCRIPT" gate --commit "$COMMIT" --runs "$R" --data "$DATA" --out "$OUT" > "$R/logs/gate.log" 2>&1
    e=$?
    echo "$e" > "$R/exit/gate"
    n=$(grep -cvE "$OK_LINE" "$R/logs/gate.log" 2>/dev/null)
    echo "${n:-unreadable}" > "$R/exit/gate.other_lines"
    finish "$e" ;;
  run)
    parallel run
    finish $? ;;
  read)
    $PY "$SCRIPT" read --commit "$COMMIT" --runs "$R" --data "$DATA" --out "$OUT" > /dev/null 2>&1
    finish $? ;;
esac
