#!/usr/bin/env bash
# Row T65 (experiments/harness/TRIES.md), one step per call, on the DGX campaign checkout ~/v/vec-campaign-d (the one
# checkout holding data/derived/feng_relabel_c/). Run from the checkout's root, detached, and wait on the exit file with
# a background poller; never a foreground sleep. Usage: t65_ms.sh STEP ROWC, with ROWC the full hash of the commit that
# placed row T65. Steps, in this order; each refuses (exit 96) unless the step before it wrote exit 0 with ROWC as its
# head:
#   tests        the targeted tests (synthetic data only)
#   cut_chamber  the chamber cut
#   cut_whole    the whole cut
#   lib_1, lib_2, lib_3   the three library tasks
#   read         once, at the commit that holds the five task files (committed unopened), a descendant of ROWC
# Exit codes, each also written to ~/v/scratch/t65/<STEP>.exit (a call without its two arguments writes usage.exit):
#   2 usage; 97 an unknown step; 96 the step before has not exited 0 at ROWC; 98 not a git checkout, or an unclean
#   tree (any change outside this row's task folder); 99 HEAD is not ROWC (tests and runs); 95 ROWC is not a commit or
#   not an ancestor of HEAD (read); 94 no row T65 in the ledger; 93 the relabelling's output is not in place; 91 the
#   step's output exists (a run never overwrites a task file, the read never overwrites a reading); 90 a task file is
#   missing or not committed (read); 89 a task file's stamp is not at ROWC (read; the check prints nothing); otherwise
#   the exit code of the step's Python.
# Each step records ROWC as its head once its checks pass, runs at 8 BLAS threads with NUMPY_MADVISE_HUGEPAGE=0 and
# TMPDIR under ~/v/scratch/t65, and writes its stdout (timings only for the runs; the outcome words for the read) to
# <STEP>.out and its stderr to <STEP>.err. Nothing here prints a log or removes a file; read <STEP>.exit, and the
# read's <STEP>.out only.
set -u
S="$HOME/v/scratch/t65"
mkdir -p "$S/tmp"
if [ "$#" -lt 2 ] || [ -z "$1" ] || [ -z "$2" ]; then
  echo "usage: t65_ms.sh tests|cut_chamber|cut_whole|lib_1|lib_2|lib_3|read ROWC" >&2
  echo 2 > "$S/usage.exit"
  exit 2
fi
STEP="$1"
ROWC="$2"
finish() { echo "$1" > "$S/$STEP.exit"; exit "$1"; }
case "$STEP" in
  tests) PREV="" ;;
  cut_chamber) PREV="tests" ;;
  cut_whole) PREV="cut_chamber" ;;
  lib_1) PREV="cut_whole" ;;
  lib_2) PREV="lib_1" ;;
  lib_3) PREV="lib_2" ;;
  read) PREV="lib_3" ;;
  *)
    echo "unknown step $STEP" >&2
    finish 97 ;;
esac
if [ -n "$PREV" ]; then
  if [ "$(cat "$S/$PREV.exit" 2>/dev/null)" != 0 ] || [ "$(cat "$S/$PREV.head" 2>/dev/null)" != "$ROWC" ]; then
    echo "refused: step $PREV has not exited 0 at $ROWC" >&2
    finish 96
  fi
fi
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 NUMPY_MADVISE_HUGEPAGE=0 TMPDIR="$S/tmp"
PY="uv run --frozen --no-sync python -I -B"
OUT=t65_ms
RES=experiments/t1_m1/results/$OUT
TASKS=$RES/tasks
if [ "$(git rev-parse --is-inside-work-tree 2>/dev/null)" != true ]; then
  echo "refused: not a git checkout" >&2
  finish 98
fi
if git status --porcelain --untracked-files=all | grep -v "^?? $TASKS/" | grep -q .; then
  echo "refused: the tree is not clean" >&2
  finish 98
fi
HEAD=$(git rev-parse HEAD)
if [ "$STEP" = read ]; then
  git cat-file -e "$ROWC^{commit}" 2>/dev/null || finish 95
  git merge-base --is-ancestor "$ROWC" "$HEAD" || finish 95
elif [ "$HEAD" != "$ROWC" ]; then
  echo "refused: HEAD is not $ROWC" >&2
  finish 99
fi
grep -q '^| T65 ' experiments/harness/TRIES.md 2>/dev/null || finish 94
[ -f data/derived/feng_relabel_c/cd1_stats_relabel.npz ] || finish 93
FILES="$TASKS/T1_E9.5_E10.5_chamber.json $TASKS/T1_E9.5_E10.5.json $TASKS/lib1.json $TASKS/lib2.json $TASKS/lib3.json"

run_step() {  # a run: its task file does not exist
  if [ -e "$TASKS/$1" ]; then finish 91; fi
  shift
  echo "$ROWC" > "$S/$STEP.head"
  $PY experiments/t1_m1/ms_read.py "$@" --out "$OUT" --commit "$ROWC" > "$S/$STEP.out" 2> "$S/$STEP.err"
  finish $?
}

case "$STEP" in
  tests)
    echo "$ROWC" > "$S/$STEP.head"
    $PY -m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp="$S/tmp/pytest" \
      tests/test_t65_ms_read.py tests/test_t64_heldmom_read.py tests/test_t62_recomputed_read.py \
      tests/test_t63_double_step_read.py tests/test_t1_bets.py tests/test_t1_scope.py tests/test_core.py \
      tests/test_g_mom_verify.py tests/test_t1_state_transport.py > "$S/$STEP.out" 2> "$S/$STEP.err"
    finish $? ;;
  cut_chamber) run_step T1_E9.5_E10.5_chamber.json run-cut --cut chamber ;;
  cut_whole) run_step T1_E9.5_E10.5.json run-cut --cut whole ;;
  lib_1) run_step lib1.json run-lib --library 1 ;;
  lib_2) run_step lib2.json run-lib --library 2 ;;
  lib_3) run_step lib3.json run-lib --library 3 ;;
  read)
    for f in $FILES; do
      git ls-files --error-unmatch "$f" > /dev/null 2>&1 || finish 90
    done
    if [ -e "$RES/t65.json" ] || [ -e "$RES/t65.txt" ]; then finish 91; fi
    # each task file's stamp commit is the row's commit; only that field is read and nothing is printed
    $PY -c 'import json,sys; sys.exit(0 if all(json.load(open(f))["stamp"]["commit"] == sys.argv[1]
            for f in sys.argv[2:]) else 1)' "$ROWC" $FILES > /dev/null 2>&1 || finish 89
    echo "$ROWC" > "$S/$STEP.head"
    $PY experiments/t1_m1/ms_read.py read --out "$OUT" --commit "$ROWC" > "$S/$STEP.out" 2> "$S/$STEP.err"
    finish $? ;;
esac
