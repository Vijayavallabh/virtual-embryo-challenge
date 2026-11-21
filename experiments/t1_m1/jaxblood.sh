#!/usr/bin/env bash
# Row T66 (experiments/harness/TRIES.md), one step per call, on the DGX campaign checkout ~/v/vec-campaign-d (the one
# checkout holding data/derived/feng_relabel_c/). Run from the checkout's root, detached, and wait on the exit file with
# a background poller; never a foreground sleep. Every step takes ROWC, the commit that holds row T66 and this code
# (full hash), and runs only at HEAD = ROWC. Steps, in order; each is refused (exit 94) unless the step before it wrote
# exit 0 for the same commit:
#   tests ROWC          the targeted tests (synthetic arrays only)
#   reduce ROWC         experiments/external/jax_blood.py: JAX embryos 1-12 and 59-63 -> data/derived/jax_blood/
#   pregate ROWC        target-free: exit 0 passes; exit 3 fails (nothing is scored, no try is counted; stop here)
#   cut chamber ROWC    the chamber cut         (after pregate)
#   cut whole ROWC      the whole cut           (after pregate)
#   lib 1|2|3 ROWC      the three library tasks (after pregate; the five task steps may run side by side)
#   read ROWC           once, after all five task steps exited 0 at ROWC
# Refusals, each with its exit code in ~/v/scratch/t66/<step>.exit (a call without its arguments writes usage.exit,
# 2): 97 an unknown step; 98 an unclean tree (any change outside this row's result folder); 96 no row T66 in the ledger;
# 95 HEAD is not ROWC; 94 the step before did not exit 0 at ROWC; 93 an input folder is not in place; 91 the step's
# output exists (nothing is ever overwritten: the reducer's folder, the pre-gate file, a task file, the reading).
# Each step records HEAD in <step>.head, runs at 8 BLAS threads with NUMPY_MADVISE_HUGEPAGE=0 and TMPDIR under
# ~/v/scratch/t66, and writes its stdout to <step>.out (timings only for the runs; the outcome words for the pre-gate
# and the read) and its stderr to <step>.err. Nothing here prints a log or removes a file.
set -u
S="$HOME/v/scratch/t66"
mkdir -p "$S/tmp"
if [ "$#" -lt 2 ] || [ -z "$1" ]; then
  echo "usage: jaxblood.sh tests|reduce|pregate|cut chamber|cut whole|lib 1|lib 2|lib 3|read ROWC" >&2
  echo 2 > "$S/usage.exit"
  exit 2
fi
STEP="$1"
shift
ARG=""
case "$STEP" in
  cut|lib)
    if [ "$#" -lt 2 ] || [ -z "$1" ] || [ -z "$2" ]; then
      echo "usage: jaxblood.sh $STEP <which> ROWC" >&2
      echo 2 > "$S/usage.exit"
      exit 2
    fi
    ARG="$1"
    shift
    NAME="${STEP}_${ARG}" ;;
  *) NAME="$STEP" ;;
esac
ROWC="$1"
finish() { echo "$1" > "$S/$NAME.exit"; exit "$1"; }

case "$NAME" in
  tests) PREV="" ;;
  reduce) PREV="tests" ;;
  pregate) PREV="reduce" ;;
  cut_chamber|cut_whole|lib_1|lib_2|lib_3) PREV="pregate" ;;
  read) PREV="cut_chamber cut_whole lib_1 lib_2 lib_3" ;;
  *)
    echo "unknown step $NAME" >&2
    finish 97 ;;
esac
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 NUMPY_MADVISE_HUGEPAGE=0 TMPDIR="$S/tmp"
PY="uv run --frozen --no-sync python -I -B"
OUT=t66_jaxblood
RES=experiments/t1_m1/results/$OUT
TASKS=$RES/tasks
if git status --porcelain --untracked-files=all | grep -v "^?? $RES/" | grep -q .; then
  echo "refused: the tree is not clean" >&2
  finish 98
fi
HEAD=$(git rev-parse HEAD)
echo "$HEAD" > "$S/$NAME.head"
grep -q '^| T66 ' experiments/harness/TRIES.md || finish 96
ROWFULL=$(git rev-parse --verify --quiet "$ROWC^{commit}") || finish 95
[ "$HEAD" = "$ROWFULL" ] || finish 95
for P in $PREV; do
  if [ "$(cat "$S/$P.exit" 2>/dev/null)" != 0 ] || [ "$(cat "$S/$P.head" 2>/dev/null)" != "$HEAD" ]; then
    echo "refused: step $P has not exited 0 at $HEAD" >&2
    finish 94
  fi
done
[ -f data/derived/feng_relabel_c/cd1_stats_relabel.npz ] || finish 93

run_step() {  # a task run: its task file does not exist
  [ -e "$TASKS/$1" ] && finish 91
  shift
  $PY experiments/t1_m1/jaxblood_read.py "$@" --out "$OUT" --commit "$HEAD" > "$S/$NAME.out" 2> "$S/$NAME.err"
  finish $?
}

case "$NAME" in
  tests)
    $PY -m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp="$S/tmp/pytest" \
      tests/test_jaxblood_read.py tests/test_jax_blood.py tests/test_t64_heldmom_read.py \
      tests/test_t62_recomputed_read.py tests/test_t1_bets.py tests/test_t1_scope.py tests/test_core.py \
      tests/test_g_mom_verify.py tests/test_t1_state_transport.py > "$S/$NAME.out" 2> "$S/$NAME.err"
    finish $? ;;
  reduce)
    [ -f data/external/jax_qiu2024/MANIFEST.json ] || finish 93
    [ -e data/derived/jax_blood ] && finish 91
    $PY experiments/external/jax_blood.py --commit "$HEAD" > "$S/$NAME.out" 2> "$S/$NAME.err"
    finish $? ;;
  pregate)
    [ -f data/derived/jax_blood/MANIFEST.json ] || finish 93
    [ -e "$TASKS/pregate.json" ] && finish 91
    $PY experiments/t1_m1/jaxblood_read.py pregate --out "$OUT" --commit "$HEAD" > "$S/$NAME.out" 2> "$S/$NAME.err"
    finish $? ;;
  cut_chamber) run_step T1_E9.5_E10.5_chamber.json run-cut --cut chamber ;;
  cut_whole) run_step T1_E9.5_E10.5.json run-cut --cut whole ;;
  lib_1|lib_2|lib_3) run_step "lib$ARG.json" run-lib --library "$ARG" ;;
  read)
    { [ -e "$RES/t66.json" ] || [ -e "$RES/t66.txt" ]; } && finish 91
    $PY experiments/t1_m1/jaxblood_read.py read --out "$OUT" --commit "$HEAD" > "$S/$NAME.out" 2> "$S/$NAME.err"
    finish $? ;;
esac
