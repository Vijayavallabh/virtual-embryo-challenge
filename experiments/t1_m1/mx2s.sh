#!/usr/bin/env bash
# Row T72 (experiments/harness/TRIES.md), candidate MX2s, one step per call, run from the root of the checkout or
# detached worktree it is called in (it moves to that checkout's top level first), with that folder's own environment
# and data. The folder needs data/derived/feng_relabel_c/, the released T1 stages under data/raw/T1/, the frozen
# experiments/t1_ext/atlases.yaml and experiments/t1_ext/results/extra_competing.json. Wait on the exit file with a
# background poller; never a foreground sleep. Every step takes ROWC, the full 40-character hash of the commit that
# holds row T72 and this code, and runs only at HEAD = ROWC. Steps, in the row's order; each is refused (exit 94)
# unless the step before it wrote exit 0 for the same commit, and each runs once (its exit record refuses a second
# call, 91):
#   tests ROWC          the targeted tests (synthetic data only)
#   fetch ROWC          experiments/t1_ext/fetch.py: the pinned samples' matrices -> data/external/t1_ext_t72/
#   load ROWC           experiments/t1_ext/load.py
#   tables ROWC         experiments/t1_ext/tables.py -> data/derived/t1_ext_t72/tables/ (with its MANIFEST.json)
#   pregate ROWC        target-free: exit 0 passes; exit 3 fails (nothing is scored, no try is counted; the row ends)
#   cut chamber ROWC    the chamber cut         (after pregate)
#   cut whole ROWC      the whole cut           (after cut chamber)
#   lib 1|2|3 ROWC      the three library tasks (after cut whole, lib 1, lib 2)
#   read ROWC           once, after pregate exited 0 at ROWC and either lib 3 (or the runs entry) exited 0 at ROWC
#                       or the chain stopped: a task step's exit record at ROWC (cut_chamber, cut_whole, lib_1 to 3,
#                       runs_<step>) is neither 0 nor running, at least one start record exists in the tasks folder,
#                       and no step's exit record says running (the row's one VOID read after a stop; 94 otherwise)
#   runs ROWC           the one allowed restart, once in all for this row: after cut chamber refused before anything
#                       loaded (its exit is not 0), while no task file and no start record exists and the read has not
#                       run; it checks pregate.json and the arm tables (mx2s_read.py check-pregate), then runs cut
#                       chamber, cut whole and lib 1 to 3 in order, each only after the one before exited 0, writing
#                       runs_<step>.exit and runs_<step>.head
# Refusals write their code to ~/v/scratch/t72/<step>.refused and never touch <step>.exit or <step>.head (a call
# without its arguments writes usage.refused, 2): 97 an unknown step; 98 an unclean tree (any change, or any untracked
# path outside experiments/t1_m1/results/t72_mx2s/, whose untracked files are the run's own: pregate.json, the start
# records and the task files stay untracked there until the read, nothing deletes them, and they are committed with
# the read's files); 96 not exactly one row T72 in the ledger; 95 ROWC is not a full hash or HEAD is not ROWC; 94 the step
# before did not exit 0 at ROWC; 93 an input is not in place; 92 the restart's conditions do not hold; 91 the step's
# output exists or the step ran before (nothing is ever overwritten).
# A step that runs writes "running" to <step>.exit and HEAD to <step>.head first, then its exit code to <step>.exit;
# it runs with $HOME/.local/bin first on PATH, at 8 BLAS threads with NUMPY_MADVISE_HUGEPAGE=0 and TMPDIR under
# ~/v/scratch/t72, and writes its stdout to <step>.out (the fixed words for the pre-gate and the read, a task file's
# path for the runs) and its stderr to <step>.err. Nothing here prints a log or removes a file.
set -u
S="$HOME/v/scratch/t72"
mkdir -p "$S/tmp"
if [ "$#" -lt 2 ] || [ -z "$1" ]; then
  echo "usage: mx2s.sh tests|fetch|load|tables|pregate|cut chamber|cut whole|lib 1|lib 2|lib 3|read|runs ROWC" >&2
  echo 2 > "$S/usage.refused"
  exit 2
fi
STEP="$1"
shift
ARG=""
case "$STEP" in
  cut|lib)
    if [ "$#" -lt 2 ] || [ -z "$1" ] || [ -z "$2" ]; then
      echo "usage: mx2s.sh $STEP <which> ROWC" >&2
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
  fetch) PREV="tests" ;;
  load) PREV="fetch" ;;
  tables) PREV="load" ;;
  pregate) PREV="tables" ;;
  cut_chamber) PREV="pregate" ;;
  cut_whole) PREV="cut_chamber" ;;
  lib_1) PREV="cut_whole" ;;
  lib_2) PREV="lib_1" ;;
  lib_3) PREV="lib_2" ;;
  read) PREV="pregate" ;;
  runs) PREV="pregate" ;;
  *)
    echo "unknown step $NAME" >&2
    refuse 97 ;;
esac
TOP=$(git rev-parse --show-toplevel 2>/dev/null) || refuse 98
cd "$TOP" || refuse 98
export PATH="$HOME/.local/bin:$PATH"
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 NUMPY_MADVISE_HUGEPAGE=0 TMPDIR="$S/tmp"
PY="uv run --frozen --no-sync python -I -B"
OUT=t72_mx2s
RES=experiments/t1_m1/results/$OUT
TASKS=$RES/tasks
EXT=data/external/t1_ext_t72
DER=data/derived/t1_ext_t72
if git status --porcelain --untracked-files=all | grep -v "^?? $RES/" | grep -q .; then
  echo "refused: the tree is not clean" >&2
  refuse 98
fi
HEAD=$(git rev-parse HEAD)
[ "$(grep -c '^| T72 ' experiments/harness/TRIES.md)" = 1 ] || refuse 96
case "$ROWC" in
  *[!0-9a-f]*) refuse 95 ;;
esac
[ "${#ROWC}" -eq 40 ] || refuse 95
[ "$HEAD" = "$ROWC" ] || refuse 95
passed() {  # step $1 exited 0 at HEAD
  [ "$(cat "$S/$1.exit" 2>/dev/null)" = 0 ] && [ "$(cat "$S/$1.head" 2>/dev/null)" = "$HEAD" ]
}
for P in $PREV; do
  if ! passed "$P"; then
    echo "refused: step $P has not exited 0 at $HEAD" >&2
    refuse 94
  fi
done
[ -e "$S/$NAME.exit" ] && refuse 91
if [ "$NAME" = read ]; then  # after pregate: the chain ended (lib 3 or runs exited 0) or stopped
  RUNNING=no
  for F in "$S"/*.exit; do
    [ "$(cat "$F" 2>/dev/null)" = running ] && RUNNING=yes
  done
  STOPPED=no
  for T in cut_chamber cut_whole lib_1 lib_2 lib_3 runs_cut_chamber runs_cut_whole runs_lib_1 runs_lib_2 runs_lib_3; do
    E=$(cat "$S/$T.exit" 2>/dev/null)
    if [ -n "$E" ] && [ "$E" != 0 ] && [ "$E" != running ] && [ "$(cat "$S/$T.head" 2>/dev/null)" = "$HEAD" ]; then
      STOPPED=yes
    fi
  done
  BEGUN=no
  for F in "$TASKS"/*.started; do
    [ -e "$F" ] && BEGUN=yes
  done
  if [ "$RUNNING" = yes ]; then
    echo "refused: a step's exit record says running" >&2
    refuse 94
  fi
  if ! passed lib_3 && ! passed runs && ! { [ "$STOPPED" = yes ] && [ "$BEGUN" = yes ]; }; then
    echo "refused: lib 3 and runs have not exited 0 at $HEAD, and no task step stopped after a start record" >&2
    refuse 94
  fi
fi
for F in data/derived/feng_relabel_c/cd1_stats_relabel.npz data/raw/T1/E8.5_RNA.h5ad data/raw/T1/E9.5_RNA.h5ad \
    data/raw/T1/E10.5_RNA.h5ad experiments/t1_ext/atlases.yaml experiments/t1_ext/results/extra_competing.json; do
  [ -f "$F" ] || refuse 93
done
TASK_FILES="T1_E9.5_E10.5_chamber.json T1_E9.5_E10.5.json lib1.json lib2.json lib3.json"

run_step() {  # a task run: neither its task file nor its start record exists
  [ -e "$TASKS/$1" ] && refuse 91
  [ -e "$TASKS/${1%.json}.started" ] && refuse 91
  shift
  start
  $PY experiments/t1_m1/mx2s_read.py "$@" --out "$OUT" --commit "$HEAD" > "$S/$NAME.out" 2> "$S/$NAME.err"
  finish $?
}

case "$NAME" in
  tests)
    start
    $PY -m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp="$S/tmp/pytest" \
      tests/test_mx2s_read.py tests/test_t1_ext_load.py tests/test_t1_ext_stage.py tests/test_t1_ext_tables.py \
      tests/test_cmjax_read.py tests/test_jax_cm.py tests/test_t62_recomputed_read.py \
      tests/test_t63_double_step_read.py tests/test_t1_bets.py tests/test_t1_scope.py tests/test_core.py \
      tests/test_g_mom_verify.py tests/test_t1_state_transport.py \
      > "$S/$NAME.out" 2> "$S/$NAME.err"
    finish $? ;;
  fetch)
    { [ -e "$EXT" ] || [ -e "$DER" ]; } && refuse 91
    start
    $PY experiments/t1_ext/fetch.py fetch --row T72 --commit "$HEAD" > "$S/$NAME.out" 2> "$S/$NAME.err"
    finish $? ;;
  load)
    [ -d "$EXT" ] || refuse 93
    start
    $PY experiments/t1_ext/load.py --row T72 --commit "$HEAD" > "$S/$NAME.out" 2> "$S/$NAME.err"
    finish $? ;;
  tables)
    [ -e "$DER/tables" ] && refuse 91
    start
    $PY experiments/t1_ext/tables.py --row T72 --commit "$HEAD" > "$S/$NAME.out" 2> "$S/$NAME.err"
    finish $? ;;
  pregate)
    [ -f "$DER/tables/MANIFEST.json" ] || refuse 93
    { [ -e "$TASKS/pregate.json" ] || [ -e "$DER/consensus" ]; } && refuse 91
    start
    $PY experiments/t1_m1/mx2s_read.py pregate --out "$OUT" --commit "$HEAD" > "$S/$NAME.out" 2> "$S/$NAME.err"
    finish $? ;;
  cut_chamber) run_step T1_E9.5_E10.5_chamber.json run-cut --cut chamber ;;
  cut_whole) run_step T1_E9.5_E10.5.json run-cut --cut whole ;;
  lib_1|lib_2|lib_3) run_step "lib$ARG.json" run-lib --library "$ARG" ;;
  read)
    { [ -e "$RES/t72.json" ] || [ -e "$RES/t72.txt" ]; } && refuse 91
    start
    $PY experiments/t1_m1/mx2s_read.py read --out "$OUT" --commit "$HEAD" > "$S/$NAME.out" 2> "$S/$NAME.err"
    finish $? ;;
  runs)
    C=$(cat "$S/cut_chamber.exit" 2>/dev/null)
    { [ -z "$C" ] || [ "$C" = 0 ] || [ "$C" = running ]; } && refuse 92
    [ "$(cat "$S/cut_chamber.head" 2>/dev/null)" = "$HEAD" ] || refuse 92
    { [ -e "$S/read.exit" ] || [ -e "$RES/t72.json" ] || [ -e "$RES/t72.txt" ]; } && refuse 92
    for T in $TASK_FILES; do
      { [ -e "$TASKS/$T" ] || [ -e "$TASKS/${T%.json}.started" ]; } && refuse 92
    done
    for T in "$TASKS"/*.started; do
      [ -e "$T" ] && refuse 92
    done
    $PY experiments/t1_m1/mx2s_read.py check-pregate --out "$OUT" --commit "$HEAD" > "$S/runs_check.out" \
      2> "$S/runs_check.err" || refuse 92
    start
    rc=0
    for R in "cut_chamber run-cut --cut chamber" "cut_whole run-cut --cut whole" "lib_1 run-lib --library 1" \
        "lib_2 run-lib --library 2" "lib_3 run-lib --library 3"; do
      set -- $R
      N="runs_$1"
      shift
      echo running > "$S/$N.exit"
      echo "$HEAD" > "$S/$N.head"
      $PY experiments/t1_m1/mx2s_read.py "$@" --out "$OUT" --commit "$HEAD" > "$S/$N.out" 2> "$S/$N.err"
      rc=$?
      echo "$rc" > "$S/$N.exit"
      [ "$rc" = 0 ] || break
    done
    finish "$rc" ;;
esac
