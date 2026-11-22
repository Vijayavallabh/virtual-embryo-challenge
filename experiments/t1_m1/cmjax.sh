#!/usr/bin/env bash
# Row T71 (experiments/harness/TRIES.md), one step per call, run from the root of the checkout or detached worktree it
# is called in (it moves to that checkout's top level first), with that folder's own environment and data. The folder
# needs data/derived/feng_relabel_c/ and data/external/jax_qiu2024/. Wait on the exit file with a background poller;
# never a foreground sleep. Every step takes ROWC, the full 40-character hash of the commit that holds row T71 and this
# code, and runs only at HEAD = ROWC. Steps, in order; each is refused (exit 94) unless the step before it wrote exit 0
# for the same commit:
#   tests ROWC          the targeted tests (synthetic arrays only)
#   reduce ROWC         experiments/external/jax_cm.py --mode practice -> data/derived/jax_cm_t71_practice/
#   pregate ROWC        target-free: exit 0 passes; exit 3 fails (nothing is scored, no try is counted; stop here)
#   cut chamber ROWC    the chamber cut         (after pregate)
#   cut whole ROWC      the whole cut           (after pregate)
#   lib 1|2|3 ROWC      the three library tasks (after pregate; the five task steps may run side by side)
#   read ROWC           once, after all five task steps exited 0 at ROWC
# Refusals write their code to ~/v/scratch/t71/<step>.refused and never touch <step>.exit or <step>.head (a call
# without its arguments writes usage.refused, 2): 97 an unknown step; 98 an unclean tree (any change outside this row's
# result folder); 96 no row T71 in the ledger; 95 ROWC is not a full hash or HEAD is not ROWC; 94 the step before did
# not exit 0 at ROWC; 93 an input folder is not in place; 91 the step's output exists (nothing is ever overwritten).
# A step that runs writes "running" to <step>.exit and HEAD to <step>.head first, then its exit code to <step>.exit;
# it runs at 8 BLAS threads with NUMPY_MADVISE_HUGEPAGE=0 and TMPDIR under ~/v/scratch/t71, and writes its stdout to
# <step>.out (timings only for the runs; the outcome words for the pre-gate and the read) and its stderr to
# <step>.err. Nothing here prints a log or removes a file.
set -u
S="$HOME/v/scratch/t71"
mkdir -p "$S/tmp"
if [ "$#" -lt 2 ] || [ -z "$1" ]; then
  echo "usage: cmjax.sh tests|reduce|pregate|cut chamber|cut whole|lib 1|lib 2|lib 3|read ROWC" >&2
  echo 2 > "$S/usage.refused"
  exit 2
fi
STEP="$1"
shift
ARG=""
case "$STEP" in
  cut|lib)
    if [ "$#" -lt 2 ] || [ -z "$1" ] || [ -z "$2" ]; then
      echo "usage: cmjax.sh $STEP <which> ROWC" >&2
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
  reduce) PREV="tests" ;;
  pregate) PREV="reduce" ;;
  cut_chamber|cut_whole|lib_1|lib_2|lib_3) PREV="pregate" ;;
  read) PREV="cut_chamber cut_whole lib_1 lib_2 lib_3" ;;
  *)
    echo "unknown step $NAME" >&2
    refuse 97 ;;
esac
TOP=$(git rev-parse --show-toplevel 2>/dev/null) || refuse 98
cd "$TOP" || refuse 98
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 NUMPY_MADVISE_HUGEPAGE=0 TMPDIR="$S/tmp"
PY="uv run --frozen --no-sync python -I -B"
OUT=t71_cmjax
RES=experiments/t1_m1/results/$OUT
TASKS=$RES/tasks
RED=data/derived/jax_cm_t71_practice
if git status --porcelain --untracked-files=all | grep -v "^?? $RES/" | grep -q .; then
  echo "refused: the tree is not clean" >&2
  refuse 98
fi
HEAD=$(git rev-parse HEAD)
grep -q '^| T71 ' experiments/harness/TRIES.md || refuse 96
case "$ROWC" in
  *[!0-9a-f]*) refuse 95 ;;
esac
[ "${#ROWC}" -eq 40 ] || refuse 95
[ "$HEAD" = "$ROWC" ] || refuse 95
for P in $PREV; do
  if [ "$(cat "$S/$P.exit" 2>/dev/null)" != 0 ] || [ "$(cat "$S/$P.head" 2>/dev/null)" != "$HEAD" ]; then
    echo "refused: step $P has not exited 0 at $HEAD" >&2
    refuse 94
  fi
done
[ -f data/derived/feng_relabel_c/cd1_stats_relabel.npz ] || refuse 93

run_step() {  # a task run: its task file does not exist
  [ -e "$TASKS/$1" ] && refuse 91
  shift
  start
  $PY experiments/t1_m1/cmjax_read.py "$@" --out "$OUT" --commit "$HEAD" > "$S/$NAME.out" 2> "$S/$NAME.err"
  finish $?
}

case "$NAME" in
  tests)
    start
    $PY -m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp="$S/tmp/pytest" \
      tests/test_cmjax_read.py tests/test_jax_cm.py tests/test_jax_blood.py tests/test_jaxblood_read.py \
      tests/test_t62_recomputed_read.py tests/test_t63_double_step_read.py tests/test_t1_bets.py \
      tests/test_t1_scope.py tests/test_core.py tests/test_g_mom_verify.py tests/test_t1_state_transport.py \
      > "$S/$NAME.out" 2> "$S/$NAME.err"
    finish $? ;;
  reduce)
    [ -f data/external/jax_qiu2024/MANIFEST.json ] || refuse 93
    [ -e "$RED" ] && refuse 91
    start
    $PY experiments/external/jax_cm.py --mode practice --commit "$HEAD" > "$S/$NAME.out" 2> "$S/$NAME.err"
    finish $? ;;
  pregate)
    [ -f "$RED/MANIFEST.json" ] || refuse 93
    [ -e "$TASKS/pregate.json" ] && refuse 91
    start
    $PY experiments/t1_m1/cmjax_read.py pregate --out "$OUT" --commit "$HEAD" > "$S/$NAME.out" 2> "$S/$NAME.err"
    finish $? ;;
  cut_chamber) run_step T1_E9.5_E10.5_chamber.json run-cut --cut chamber ;;
  cut_whole) run_step T1_E9.5_E10.5.json run-cut --cut whole ;;
  lib_1|lib_2|lib_3) run_step "lib$ARG.json" run-lib --library "$ARG" ;;
  read)
    { [ -e "$RES/t71.json" ] || [ -e "$RES/t71.txt" ]; } && refuse 91
    start
    $PY experiments/t1_m1/cmjax_read.py read --out "$OUT" --commit "$HEAD" > "$S/$NAME.out" 2> "$S/$NAME.err"
    finish $? ;;
esac
