#!/usr/bin/env bash
# Row T73 (experiments/harness/TRIES.md): gate G0 of row R58's file T1_s2rcMSBm20, experiments/t1_m1/m20_read.py, one
# step per call, run from the root of the checkout or detached worktree it is called in (it moves to that checkout's
# top level first), with that folder's own environment and data; the folder needs data/derived/feng_relabel_c/ (row
# N54's output) and the released T1 files. CPU only. Wait on the exit file with a background poller; never a
# foreground sleep. Every step takes ROWC, the full 40-character hash of the commit that holds row T73 and this code,
# and runs only at HEAD = ROWC. Steps, in order; each is refused (exit 94) unless the step before it wrote exit 0 for
# the same commit:
#   tests ROWC          the targeted tests (synthetic arrays only)
#   cut chamber ROWC    the chamber cut         (after tests; the fifteen tries count from the first task step)
#   cut whole ROWC      the whole cut           (after tests)
#   lib 1|2|3 ROWC      the three library tasks (after tests; the five task steps may run side by side)
#   read ROWC           once, after all five task steps exited 0 at ROWC, with the task files uncommitted in the
#                       row's task folder; the main session commits them and the reading together and quotes the
#                       outcome line ("Outcome: T73 PASS (m20's step is AHEAD)", "Outcome: T73 STOP (...)" or
#                       "Outcome: T73 VOID (...)") in row T73's outcome cell, which
#                       experiments/builds/t1_lm_m20_line.sh reads by code
# The gate's 0.02 tolerance is dry-run before the row is committed, on committed results only and with no row (J7):
#   uv run --frozen --no-sync python -I -B experiments/t1_m1/m20_read.py j7
# Refusals write their code to ~/v/scratch/t73/<step>.refused and never touch <step>.exit or <step>.head (a call
# without its arguments writes usage.refused, 2): 97 an unknown step; 98 an unclean tree (any change outside this row's
# result folder); 96 TRIES.md does not hold exactly one row T73 whose first cell starts 'T73 (' and names
# experiments/t1_m1/m20_read.py and experiments/t1_m1/t73_m20.sh, each as a whole path; 95 ROWC is not a full hash or
# HEAD is not ROWC; 94 the step before did not exit 0 at ROWC; 93 an input folder is not in place; 91 the step's output
# exists (nothing is ever overwritten); 88 the checkout's .venv does not name its own src. A step that runs writes
# "running" to <step>.exit and HEAD to <step>.head first, then its exit code to <step>.exit; it runs at 8 BLAS threads
# with NUMPY_MADVISE_HUGEPAGE=0, TMPDIR under ~/v/scratch/t73 and $HOME/.local/bin first on PATH, and writes its stdout
# to <step>.out (timings only for the runs; the outcome words for the read) and its stderr to <step>.err. Nothing here
# prints a log or removes a file.
set -u
export PATH=$HOME/.local/bin:$PATH
S="$HOME/v/scratch/t73"
mkdir -p "$S/tmp"
if [ "$#" -lt 2 ] || [ -z "$1" ]; then
  echo "usage: t73_m20.sh tests|cut chamber|cut whole|lib 1|lib 2|lib 3|read ROWC" >&2
  echo 2 > "$S/usage.refused"
  exit 2
fi
STEP="$1"
shift
ARG=""
case "$STEP" in
  cut|lib)
    if [ "$#" -lt 2 ] || [ -z "$1" ] || [ -z "$2" ]; then
      echo "usage: t73_m20.sh $STEP <which> ROWC" >&2
      echo 2 > "$S/usage.refused"
      exit 2
    fi
    ARG="$1"
    shift
    NAME="${STEP}_${ARG}" ;;
  *) NAME="$STEP" ;;
esac
ROWC="$1"
refuse() { echo "$1" > "$S/$NAME.refused"; exit "$1"; }
start() { echo running > "$S/$NAME.exit"; echo "$HEAD" > "$S/$NAME.head"; }
finish() { echo "$1" > "$S/$NAME.exit"; exit "$1"; }

case "$NAME" in
  tests) PREV="" ;;
  cut_chamber|cut_whole|lib_1|lib_2|lib_3) PREV="tests" ;;
  read) PREV="cut_chamber cut_whole lib_1 lib_2 lib_3" ;;
  *)
    echo "unknown step" >&2
    NAME=unknown
    refuse 97 ;;
esac
TOP=$(git rev-parse --show-toplevel 2>/dev/null) || refuse 98
cd "$TOP" || refuse 98
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 NUMPY_MADVISE_HUGEPAGE=0 TMPDIR="$S/tmp"
export UV_PROJECT_ENVIRONMENT="$PWD/.venv"
PY="uv run --frozen --no-sync python -I -B"
TRIES=experiments/harness/TRIES.md
OUT=t73_m20
RES=experiments/t1_m1/results/$OUT
TASKS=$RES/tasks
ROW_PATHS=('experiments/t1_m1/m20_read\.py' 'experiments/t1_m1/t73_m20\.sh')  # escaped for grep -E
names() { printf '%s\n' "$1" | grep -qE "(^|[^[:alnum:]_./-])$2([^[:alnum:]_./-]|\$)"; }
if git status --porcelain --untracked-files=all | grep -v "^?? $RES/" | grep -q .; then
  echo "refused: the tree is not clean" >&2
  refuse 98
fi
HEAD=$(git rev-parse HEAD)
[ "$(grep -c '^| T73 ' "$TRIES" 2>/dev/null)" = 1 ] || { echo "refused: TRIES.md does not hold exactly one row T73" >&2; refuse 96; }
first=$(grep '^| T73 ' "$TRIES" | cut -d'|' -f2)
printf '%s\n' "$first" | grep -q '^ T73 (' || { echo "refused: row T73's first cell does not start 'T73 ('" >&2; refuse 96; }
for p in "${ROW_PATHS[@]}"; do
  names "$first" "$p" || { echo "refused: row T73's first cell does not name $p as a whole path" >&2; refuse 96; }
done
case "$ROWC" in
  *[!0-9a-f]*) refuse 95 ;;
esac
[ "${#ROWC}" -eq 40 ] || refuse 95
[ "$HEAD" = "$ROWC" ] || { echo "refused: HEAD is not $ROWC" >&2; refuse 95; }
for P in $PREV; do
  if [ "$(cat "$S/$P.exit" 2>/dev/null)" != 0 ] || [ "$(cat "$S/$P.head" 2>/dev/null)" != "$HEAD" ]; then
    echo "refused: step $P has not exited 0 at $HEAD" >&2
    refuse 94
  fi
done
venv=0
for pth in .venv/lib/python3*/site-packages/vec.pth; do
  if [ -f "$pth" ] && [ "$(readlink -f "$(< "$pth")")" = "$(readlink -f src)" ]; then venv=1; fi
done
[ "$venv" = 1 ] || { echo "refused: this checkout has no .venv of its own whose vec.pth names its src" >&2; refuse 88; }
[ -f data/derived/feng_relabel_c/cd1_stats_relabel.npz ] \
  || { echo "refused: the relabelling's output is not in data/derived/feng_relabel_c" >&2; refuse 93; }

run_step() {  # a task run: its task file does not exist
  [ -e "$TASKS/$1" ] && { echo "refused: $TASKS/$1 exists" >&2; refuse 91; }
  shift
  start
  $PY experiments/t1_m1/m20_read.py "$@" --out "$OUT" --commit "$HEAD" > "$S/$NAME.out" 2> "$S/$NAME.err"
  finish $?
}

case "$NAME" in
  tests)
    start
    $PY -m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp="$S/tmp/pytest" \
      tests/test_t73_m20_read.py tests/test_t65_ms_read.py tests/test_t1_lm_m20_launchers.py \
      tests/test_t62_recomputed_read.py tests/test_t1_scope.py tests/test_core.py tests/test_g_mom_verify.py \
      tests/test_t1_state_transport.py > "$S/$NAME.out" 2> "$S/$NAME.err"
    finish $? ;;
  cut_chamber) run_step T1_E9.5_E10.5_chamber.json run-cut --cut chamber ;;
  cut_whole) run_step T1_E9.5_E10.5.json run-cut --cut whole ;;
  lib_1|lib_2|lib_3) run_step "lib$ARG.json" run-lib --library "$ARG" ;;
  read)
    { [ -e "$RES/t73.json" ] || [ -e "$RES/t73.txt" ]; } && { echo "refused: a reading exists" >&2; refuse 91; }
    start
    $PY experiments/t1_m1/m20_read.py read --out "$OUT" --commit "$HEAD" > "$S/$NAME.out" 2> "$S/$NAME.err"
    finish $? ;;
esac
