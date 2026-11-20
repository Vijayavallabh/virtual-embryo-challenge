#!/usr/bin/env bash
# Row T64 (experiments/harness/TRIES.md), one step per call, on the DGX campaign checkout ~/v/vec-campaign-d (the one
# checkout holding data/derived/feng_relabel_c/). Run from the checkout's root, detached, and wait on the exit file with
# a background poller; never a foreground sleep. Steps, in order, each with the row's commit ROWC (the commit that
# placed row T64):
#   tests ROWC          the targeted tests (synthetic data only)
#   cut chamber ROWC    the chamber cut (refused unless tests.exit is 0 at this HEAD)
#   cut whole ROWC      the whole cut
#   lib 1|2|3 ROWC      the three library tasks
#   read ROWC           once, at the commit that holds the five task files (committed unopened), a descendant of ROWC
# Refusals, each with its exit code in ~/v/scratch/t64/<step>.exit: 98 an unclean tree (any change outside this row's
# task folder); 96 no row T64 in the ledger; 95 HEAD is not ROWC (tests and runs) or ROWC is not an ancestor (read);
# 93 the relabelling's output is not in place; 92 the tests did not pass at this HEAD; 91 the step's output exists
# (a run never overwrites a task file, the read never overwrites a reading); 90 a task file is missing or not
# committed (read); 94 a task file's stamp is not at ROWC (read; the check prints nothing). 97 an unknown step.
# Each step records HEAD, runs at 8 BLAS threads with NUMPY_MADVISE_HUGEPAGE=0 and TMPDIR under ~/v/scratch/t64, and
# writes its stdout (timings only for the runs; the outcome words for the read) to <step>.out and its stderr to
# <step>.err. Nothing here prints a log; read <step>.exit, and the read's <step>.out only.
set -u
STEP=${1:?usage: t64_heldmom.sh tests|cut chamber|cut whole|lib 1|lib 2|lib 3|read ROWC}
shift
case "$STEP" in
  cut|lib) ARG=${1:?the cut or library}; shift; NAME="${STEP}_${ARG}" ;;
  *) ARG=""; NAME="$STEP" ;;
esac
ROWC=${1:?the commit of the row}
S=~/v/scratch/t64
mkdir -p "$S/tmp"
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 NUMPY_MADVISE_HUGEPAGE=0 TMPDIR="$S/tmp"
PY="uv run --frozen --no-sync python -I -B"
OUT=t64_heldmom
TASKS=experiments/t1_m1/results/$OUT/tasks
finish() { echo "$1" > "$S/$NAME.exit"; exit "$1"; }

case "$NAME" in
  tests|cut_chamber|cut_whole|lib_1|lib_2|lib_3|read) ;;
  *) echo "unknown step $NAME" >&2; finish 97 ;;
esac
if git status --porcelain --untracked-files=all | grep -v "^?? $TASKS/" | grep -q .; then
  echo "refused: the tree is not clean" >&2
  finish 98
fi
HEAD=$(git rev-parse HEAD)
echo "$HEAD" > "$S/$NAME.head"
grep -q '^| T64 ' experiments/harness/TRIES.md || finish 96
ROWFULL=$(git rev-parse --verify --quiet "$ROWC^{commit}") || finish 95
[ -f data/derived/feng_relabel_c/cd1_stats_relabel.npz ] || finish 93
FILES="$TASKS/T1_E9.5_E10.5_chamber.json $TASKS/T1_E9.5_E10.5.json $TASKS/lib1.json $TASKS/lib2.json $TASKS/lib3.json"

run_step() {  # a run: HEAD is the row's commit, the tests passed at it, its task file does not exist
  [ "$HEAD" = "$ROWFULL" ] || finish 95
  [ "$(cat "$S/tests.exit" 2>/dev/null)" = 0 ] && [ "$(cat "$S/tests.head" 2>/dev/null)" = "$HEAD" ] || finish 92
  [ -e "$TASKS/$1" ] && finish 91
  shift
  $PY experiments/t1_m1/heldmom_read.py "$@" --out "$OUT" > "$S/$NAME.out" 2> "$S/$NAME.err"
  finish $?
}

case "$NAME" in
  tests)
    [ "$HEAD" = "$ROWFULL" ] || finish 95
    $PY -m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp="$S/tmp/pytest" \
      tests/test_t64_heldmom_read.py tests/test_t62_recomputed_read.py tests/test_t63_double_step_read.py \
      tests/test_t1_bets.py tests/test_t1_scope.py tests/test_core.py tests/test_g_mom_verify.py \
      tests/test_t1_state_transport.py > "$S/$NAME.out" 2> "$S/$NAME.err"
    finish $? ;;
  cut_chamber) run_step T1_E9.5_E10.5_chamber.json run-cut --cut chamber ;;
  cut_whole) run_step T1_E9.5_E10.5.json run-cut --cut whole ;;
  lib_1|lib_2|lib_3) run_step "lib$ARG.json" run-lib --library "$ARG" ;;
  read)
    git merge-base --is-ancestor "$ROWFULL" "$HEAD" || finish 95
    for f in $FILES; do
      git ls-files --error-unmatch "$f" > /dev/null 2>&1 || finish 90
    done
    [ -e "experiments/t1_m1/results/$OUT/t64.json" ] || [ -e "experiments/t1_m1/results/$OUT/t64.txt" ] && finish 91
    # each task file's stamp commit is the row's commit; only that field is read and nothing is printed
    $PY -c 'import json,sys; sys.exit(0 if all(json.load(open(f))["stamp"]["commit"] == sys.argv[1]
            for f in sys.argv[2:]) else 1)' "$ROWFULL" $FILES > /dev/null 2>&1 || finish 94
    $PY experiments/t1_m1/heldmom_read.py read --out "$OUT" > "$S/$NAME.out" 2> "$S/$NAME.err"
    finish $? ;;
esac
