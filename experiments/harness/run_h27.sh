#!/usr/bin/env bash
# Row H27 (experiments/harness/TRIES.md), one step per call, on a DGX campaign checkout (~/v/vec-campaign, never the
# reference checkout). Run from the checkout's root, detached, and wait on the exit file with a background poller;
# never a foreground sleep. Steps, in order:
#   tests COMMIT   the unit tests of the script and of row H26's (synthetic arrays and hand-made tables only)
#   run COMMIT     the fifteen parts: per transition five processes of two seeds, the transitions one after another;
#                  each part refuses before anything loads unless row H27 is in TRIES.md, the tree is clean and its
#                  output does not exist. Logs carry "<task> seed <s> (<sec> s) done" and "wrote <path>" only.
#   read COMMIT    the one read: void conditions and the fixed statements, into experiments/harness/results/
# COMMIT is row H27's commit (full hash). Each step refuses on an unclean tree (exit 98) or another HEAD (exit 99),
# runs with 8 BLAS threads a process and NUMPY_MADVISE_HUGEPAGE=0, and writes ~/v/scratch/h27/<STEP>.exit (and, for
# run, one exit file per part under runs/h27/exit/).
set -u
STEP=${1:?usage: run_h27.sh tests|run|read COMMIT}
COMMIT=${2:?the commit row H27 names (full hash)}
S=~/v/scratch/h27
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
SCRIPT=experiments/harness/h27_heart_forecast_size.py
R=runs/h27
TASKS=("T2h/<=E9.5->E10.5:e95" "T2h/E8.25_late->E8.5:e825" "T2h/<=E8.75->E10.5:e875")

case "$STEP" in
  tests)
    $PY -m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp="$S/tmp/pytest" \
      tests/test_h27_heart_forecast_size.py tests/test_h26_heart_scope_forms.py > "$S/tests.log" 2>&1
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
      done
    done
    finish $bad ;;
  read)
    $PY "$SCRIPT" read --runs "$R" --commit "$COMMIT" --out experiments/harness/results/h27_heart_forecast_size \
      > /dev/null 2>&1
    finish $? ;;
  *)
    echo "unknown step $STEP" >&2
    finish 97 ;;
esac
