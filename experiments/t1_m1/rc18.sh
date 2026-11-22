#!/usr/bin/env bash
# Rows N58 (the E18 labelling) and T69 (rc18's practice read) of experiments/harness/TRIES.md, one step per call, on a
# DGX checkout or detached worktree whose own data/ is a real directory holding data/derived/feng_relabel_c/ (row
# N54's output) and the clean reduction under data/external. Run from that folder's root, detached, and wait on the exit
# file with a background poller; never a foreground sleep. CPU only. Two chains; each step runs only if the step before
# it exited 0 at the same commit:
#   label-tests NROWC NROW   the targeted tests at row N58's commit NROWC (N row NROW; synthetic arrays only)
#   label NROWC NROW         experiments/external/feng_relabel_e18.py --row NROW (target-free, minutes of CPU); exit 3:
#                            it stopped after preflight (recorded, not rerun). Then the main session commits its
#                            RECORD and row T69 naming the E18 statistics' sha256: that commit is ROWC
#   tests ROWC               the targeted tests at row T69's run commit ROWC (synthetic data only); needs label exit 0
#                            at a commit ROWC descends from
#   cut chamber ROWC, cut whole ROWC, lib 1 ROWC, lib 2 ROWC, lib 3 ROWC
#                            the five task runs (each needs tests exit 0 at ROWC); the task files stay uncommitted in
#                            the row's task folder until the read
#   read ROWC                the reading (needs the five runs exit 0 at ROWC); the main session commits the task files
#                            and the reading together
# COMMITs are full hashes. A step that passes its checks first writes "running" to ~/v/scratch/t69/<step>.exit, then
# its .head file, runs, and writes its exit code there. A call refused before its step runs writes
# <step>.refused.exit only (usage.refused.exit without a step), so a refusal never touches a finished step's exit:
# 2 usage, 97 unknown step, 94 the step already exited 0 at this commit, 95 HEAD is not the commit given, 98 an unclean
# tree (any change outside the row's results folders), 96 a row is not in the ledger, 92 the step before did not exit 0
# at this commit, 93 an input folder is missing (or, for label, its output folder exists), 91 the step's output
# exists (nothing is ever overwritten), 90 a task file is missing (read). The Python steps exit 0 done, 1 an error or a
# refusal, 3 the labelling stopped after preflight. Each step runs at 8 BLAS threads with NUMPY_MADVISE_HUGEPAGE=0, this
# folder's own environment and TMPDIR under ~/v/scratch/t69, and writes its stdout (step lines only) to <step>.out; any
# other line goes to <step>.err and is counted, never printed (<step>.other_lines). Nothing here deletes anything.
set -u
S="$HOME/v/scratch/t69"
mkdir -p "$S/tmp" || exit 89
STEP="${1:-}"
ARG=""
ROWC=""
NROW=""
usage() {
  echo "usage: rc18.sh label-tests or label NROWC NROW, or tests or read ROWC, or cut chamber or whole ROWC, or lib 1 2 or 3 ROWC (full hashes)" >&2
  echo 2 > "$S/usage.refused.exit"
  exit 2
}
[ -n "$STEP" ] || usage
case "$STEP" in
  cut|lib) ARG="${2:-}"; ROWC="${3:-}"; NAME="${STEP}_${ARG}" ;;
  label-tests|label) ROWC="${2:-}"; NROW="${3:-}"; NAME="${STEP/-/_}" ;;
  *) ROWC="${2:-}"; NAME="$STEP" ;;
esac
case "$NAME" in
  label_tests|label|tests|cut_chamber|cut_whole|lib_1|lib_2|lib_3|read) ;;
  *)
    echo "unknown step" >&2
    echo 97 > "$S/unknown.refused.exit"
    exit 97 ;;
esac
BAD=0
[[ "$ROWC" =~ ^[0-9a-f]{40}$ ]] || BAD=1
case "$NAME" in
  label_tests|label) [[ "$NROW" =~ ^N[0-9]+[a-z]?$ ]] || BAD=1 ;;
esac
[ "$BAD" = 0 ] || usage
refuse() { echo "$1" > "$S/$NAME.refused.exit"; exit "$1"; }
finish() { echo "$1" > "$S/$NAME.exit"; exit "$1"; }

OUT=t69_rc18
RES=experiments/t1_m1/results/$OUT
TASKS=$RES/tasks
SCRIPT=experiments/t1_m1/rc18_read.py
LABELLER=experiments/external/feng_relabel_e18.py
RELABEL=data/derived/feng_relabel_c/cd1_stats_relabel.npz
E18_DIR=data/derived/feng_relabel_e18
E18=$E18_DIR/cd1_stats_relabel_e18.npz
FILES="$TASKS/T1_E9.5_E10.5_chamber.json $TASKS/T1_E9.5_E10.5.json $TASKS/lib1.json $TASKS/lib2.json $TASKS/lib3.json"

ok_at() {  # step $1 exited 0 at commit $2
  [ "$(cat "$S/$1.exit" 2>/dev/null)" = 0 ] && [ "$(cat "$S/$1.head" 2>/dev/null)" = "$2" ]
}
if ok_at "$NAME" "$ROWC"; then
  echo "refused: step $NAME has already exited 0 at this commit" >&2
  refuse 94
fi
if [ "$(git rev-parse HEAD 2>/dev/null)" != "$ROWC" ]; then
  echo "refused: HEAD is not the commit given" >&2
  refuse 95
fi
if git status --porcelain --untracked-files=all | grep -v "^?? $RES/" | grep -q .; then
  echo "refused: the tree is not clean" >&2
  refuse 98
fi
if [ -n "$NROW" ]; then
  grep -q "^| $NROW .*feng_relabel_e18\.py" experiments/harness/TRIES.md \
    || { echo "refused: no row $NROW names the labeller" >&2; refuse 96; }
else
  grep -q '^| T69 .*rc18_read\.py' experiments/harness/TRIES.md \
    || { echo "refused: no row T69 names the read" >&2; refuse 96; }
fi
after() {  # every named step exited 0 at this commit
  local p
  for p in "$@"; do
    if ! ok_at "$p" "$ROWC"; then
      echo "refused: step $p has not exited 0 at this commit" >&2
      refuse 92
    fi
  done
}

TESTS="tests/test_t1_rc18.py tests/test_t62_recomputed_read.py tests/test_t63_double_step_read.py
  tests/test_feng_relabel.py tests/test_t1_bets.py tests/test_t1_scope.py tests/test_core.py tests/test_g_mom_verify.py
  tests/test_t1_state_transport.py"
# One test of tests/test_feng_relabel.py reads the checkout's own data folder and fails only where a real N55 run left
# its marker (the DGX checkout that holds the relabelling's output; rows N57 and T67's launcher deselect it for the same
# reason). It is deselected by name in both test steps; nothing else is.
KNOWN_DESELECT="--deselect tests/test_feng_relabel.py::test_the_dry_runs_refusals"

case "$NAME" in
  label_tests) ;;
  label)
    after label_tests
    [ -f "$RELABEL" ] || refuse 93
    [ -e "$E18_DIR" ] && refuse 93 ;;
  tests)
    if [ "$(cat "$S/label.exit" 2>/dev/null)" != 0 ] \
        || ! git merge-base --is-ancestor "$(cat "$S/label.head" 2>/dev/null)" "$ROWC" 2>/dev/null; then
      echo "refused: the labelling has not exited 0 at a commit this one descends from" >&2
      refuse 92
    fi
    [ -f "$RELABEL" ] && [ -f "$E18" ] || refuse 93 ;;
  cut_chamber|cut_whole|lib_1|lib_2|lib_3)
    after tests
    [ -f "$RELABEL" ] && [ -f "$E18" ] || refuse 93
    case "$NAME" in
      cut_chamber) F=T1_E9.5_E10.5_chamber.json ;;
      cut_whole) F=T1_E9.5_E10.5.json ;;
      *) F="lib$ARG.json" ;;
    esac
    [ -e "$TASKS/$F" ] && refuse 91 ;;
  read)
    after cut_chamber cut_whole lib_1 lib_2 lib_3
    for f in $FILES; do
      [ -f "$f" ] || refuse 90
    done
    [ -e "$RES/t69.json" ] || [ -e "$RES/t69.txt" ] && refuse 91 ;;
esac

echo running > "$S/$NAME.exit"
echo "$ROWC" > "$S/$NAME.head"
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 NUMPY_MADVISE_HUGEPAGE=0 TMPDIR="$S/tmp"
export UV_PROJECT_ENVIRONMENT="$PWD/.venv"
PY="uv run --frozen --no-sync python -I -B"
OK_LINE='^(T1/E9\.5->E10\.5(:chamber)? seed [0-9]+|library [123] seed [0-9]+): [0-9]+ s$|^E18 labelling (written|stopped)|^Outcome: |^Parts: |^What follows: |^\(numbers in '
count_other() {  # the lines of the step's stdout that are not step lines: counted, never printed
  local n
  n=$(grep -cvE "$OK_LINE" "$S/$NAME.out" 2>/dev/null) || true
  echo "${n:-unreadable}" > "$S/$NAME.other_lines"
}
run_py() {
  $PY "$@" > "$S/$NAME.out" 2> "$S/$NAME.err"
  local rc=$?
  count_other
  finish $rc
}

case "$NAME" in
  label_tests)
    run_py -m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp="$S/tmp/pytest" \
      tests/test_t1_rc18.py tests/test_feng_relabel.py $KNOWN_DESELECT ;;
  label) run_py "$LABELLER" --row "$NROW" ;;
  tests) run_py -m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp="$S/tmp/pytest" $TESTS $KNOWN_DESELECT ;;
  cut_chamber) run_py "$SCRIPT" run-cut --cut chamber --out "$OUT" ;;
  cut_whole) run_py "$SCRIPT" run-cut --cut whole --out "$OUT" ;;
  lib_1|lib_2|lib_3) run_py "$SCRIPT" run-lib --library "$ARG" --out "$OUT" ;;
  read) run_py "$SCRIPT" read --out "$OUT" ;;
esac
