#!/usr/bin/env bash
# Row T67 (experiments/harness/TRIES.md) and its late labelling, one step per call, on the DGX campaign checkout
# ~/v/vec-campaign-d (the one checkout holding data/derived/feng_relabel_c/). Run from the checkout's root, detached, and
# wait on the exit file with a background poller; never a foreground sleep. Two chains, each step refused unless the
# step before it exited 0 at the same commit:
#   late-tests NROWC NROW   the targeted tests at the late labelling's row commit NROWC (N row NROW; synthetic data)
#   late NROWC NROW         experiments/external/feng_relabel_late.py --row NROW (target-free, minutes of CPU); then the
#                           main session commits its RECORD, the X002 amendment and row T67 naming the late sha256
#   tests ROWC              the targeted tests at row T67's run commit ROWC (synthetic data only)
#   pregate ROWC            the target-free idle pre-gate: exit 0 not idle, 4 idle (the row is pruned, no try counted)
#   cut chamber ROWC        the chamber cut (needs pregate 0 at ROWC); cut whole ROWC the whole cut
#   lib 1|2|3 ROWC          the three library tasks (each needs pregate 0 at ROWC)
#   read ROWC               the reading, at ROWC with the five task files and the pre-gate uncommitted in the row's
#                           task folder (needs all five runs 0 at ROWC); the main session commits them with the reading
# Refusals, each with its exit code in ~/v/scratch/t67/<step>.exit: 2 usage (usage.exit); 97 an unknown step
# (unknown.exit); 98 an unclean tree (any change outside this row's results folder); 95 HEAD is not the commit given (full
# hash); 96 the row is not in the ledger; 92 the step before did not exit 0 at this HEAD; 93 an input folder is missing
# (or, for late, its output folder exists); 91 the step's output exists (nothing is ever overwritten); 90 a task file is
# missing (read). The Python steps exit 0 done, 1 an error or a refusal, 3 the late labelling stopped after preflight,
# 4 the pre-gate found the arm idle. Each step records HEAD, runs at 8 BLAS threads with NUMPY_MADVISE_HUGEPAGE=0 and
# TMPDIR under ~/v/scratch/t67, and writes its stdout (step lines only) to <step>.out and its stderr to <step>.err.
# Nothing here deletes anything.
set -u
S="$HOME/v/scratch/t67"
mkdir -p "$S/tmp"
STEP="${1:-}"
ARG=""
NROW=""
usage() {
  echo "usage: timing.sh late-tests|late NROWC NROW, or tests|pregate|read ROWC, or cut chamber|whole ROWC, or lib 1|2|3 ROWC (full hashes)" >&2
  echo 2 > "$S/usage.exit"
  exit 2
}
[ -n "$STEP" ] || usage
case "$STEP" in
  cut|lib) ARG="${2:-}"; ROWC="${3:-}"; NAME="${STEP}_${ARG}" ;;
  late-tests|late) ROWC="${2:-}"; NROW="${3:-}"; NAME="${STEP/-/_}" ;;
  *) ROWC="${2:-}"; NAME="$STEP" ;;
esac
case "$NAME" in
  late_tests|late|tests|pregate|cut_chamber|cut_whole|lib_1|lib_2|lib_3|read) ;;
  *)
    echo "unknown step" >&2
    echo 97 > "$S/unknown.exit"
    exit 97 ;;
esac
BAD=0
[[ "$ROWC" =~ ^[0-9a-f]{40}$ ]] || BAD=1
case "$NAME" in
  late_tests|late) [[ "$NROW" =~ ^N[0-9]+[a-z]?$ ]] || BAD=1 ;;
esac
[ "$BAD" = 0 ] || usage
finish() { echo "$1" > "$S/$NAME.exit"; exit "$1"; }
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 NUMPY_MADVISE_HUGEPAGE=0 TMPDIR="$S/tmp"
PY="uv run --frozen --no-sync python -I -B"
OUT=t67_timing
TASKS=experiments/t1_m1/results/$OUT/tasks
SCRIPT=experiments/t1_m1/timing_read.py
RELABEL=data/derived/feng_relabel_c/cd1_stats_relabel.npz
LATE_DIR=data/derived/feng_relabel_late
LATE=$LATE_DIR/cd1_stats_relabel_late.npz
FILES="$TASKS/T1_E9.5_E10.5_chamber.json $TASKS/T1_E9.5_E10.5.json $TASKS/lib1.json $TASKS/lib2.json $TASKS/lib3.json"

if git status --porcelain --untracked-files=all | grep -v "^?? experiments/t1_m1/results/$OUT/" | grep -q .; then
  echo "refused: the tree is not clean" >&2
  finish 98
fi
HEAD=$(git rev-parse HEAD)
if [ "$HEAD" != "$ROWC" ]; then
  echo "refused: HEAD is not the commit given" >&2
  finish 95
fi
echo "$HEAD" > "$S/$NAME.head"
if [ -n "$NROW" ]; then
  grep -q "^| $NROW " experiments/harness/TRIES.md || finish 96
else
  grep -q '^| T67 ' experiments/harness/TRIES.md || finish 96
fi

after() {  # every named step exited 0 at this HEAD
  local p
  for p in "$@"; do
    if [ "$(cat "$S/$p.exit" 2>/dev/null)" != 0 ] || [ "$(cat "$S/$p.head" 2>/dev/null)" != "$HEAD" ]; then
      echo "refused: step $p has not exited 0 at this HEAD" >&2
      finish 92
    fi
  done
}

# One test of tests/test_feng_relabel.py reads the checkout's own data folder and fails only where a real N55 run left
# its marker (the DGX checkout that holds the relabelling's output; the T1 candidates' launcher of 2026-11-20 deselected
# it for the same reason). It is deselected by name in both test steps; nothing else is.
KNOWN_DESELECT="--deselect tests/test_feng_relabel.py::test_the_dry_runs_refusals"
TESTS_T67="tests/test_t1_timing.py tests/test_t62_recomputed_read.py tests/test_t63_double_step_read.py
  tests/test_feng_relabel.py tests/test_t1_bets.py tests/test_t1_scope.py tests/test_core.py tests/test_g_mom_verify.py
  tests/test_t1_state_transport.py"

run_step() {  # a run: the pre-gate passed at this HEAD, both statistics are in place, its task file does not exist
  after pregate
  [ -f "$RELABEL" ] && [ -f "$LATE" ] || finish 93
  [ -e "$TASKS/$1" ] && finish 91
  shift
  $PY "$SCRIPT" "$@" --out "$OUT" > "$S/$NAME.out" 2> "$S/$NAME.err"
  finish $?
}

case "$NAME" in
  late_tests)
    $PY -m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp="$S/tmp/pytest" \
      tests/test_t1_timing.py tests/test_feng_relabel.py $KNOWN_DESELECT > "$S/$NAME.out" 2> "$S/$NAME.err"
    finish $? ;;
  late)
    after late_tests
    [ -f "$RELABEL" ] || finish 93
    [ -e "$LATE_DIR" ] && finish 93
    $PY experiments/external/feng_relabel_late.py --row "$NROW" > "$S/$NAME.out" 2> "$S/$NAME.err"
    finish $? ;;
  tests)
    $PY -m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp="$S/tmp/pytest" $TESTS_T67 $KNOWN_DESELECT \
      > "$S/$NAME.out" 2> "$S/$NAME.err"
    finish $? ;;
  pregate)
    after tests
    [ -f "$RELABEL" ] && [ -f "$LATE" ] || finish 93
    [ -e "$TASKS/pregate.json" ] || [ -e "$TASKS/pregate.txt" ] && finish 91
    $PY "$SCRIPT" pregate --out "$OUT" > "$S/$NAME.out" 2> "$S/$NAME.err"
    finish $? ;;
  cut_chamber) run_step T1_E9.5_E10.5_chamber.json run-cut --cut chamber ;;
  cut_whole) run_step T1_E9.5_E10.5.json run-cut --cut whole ;;
  lib_1|lib_2|lib_3) run_step "lib$ARG.json" run-lib --library "$ARG" ;;
  read)
    after cut_chamber cut_whole lib_1 lib_2 lib_3
    for f in $FILES; do
      [ -f "$f" ] || finish 90
    done
    [ -e "experiments/t1_m1/results/$OUT/t67.json" ] || [ -e "experiments/t1_m1/results/$OUT/t67.txt" ] && finish 91
    $PY "$SCRIPT" read --out "$OUT" > "$S/$NAME.out" 2> "$S/$NAME.err"
    finish $? ;;
esac
