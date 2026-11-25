#!/usr/bin/env bash
# Row H36 (experiments/harness/TRIES.md): the heart's XPC line (experiments/harness/h36_heart_run_profile.py), on a DGX
# campaign checkout or a detached worktree of it (never the reference checkout). Run from the root of a clean checkout
# or worktree at COMMIT, detached, and wait on done.exit with a background poller; never a foreground sleep. CPU only.
#   bash experiments/harness/run_h36.sh <COMMIT> <RUN> read           entry, profile, gate, run, read
#   bash experiments/harness/run_h36.sh <COMMIT> <RUN> test-profile   the test profile (T), at a later commit
# COMMIT is the full hash HEAD must be. For `read` it is the run commit: it holds row H36, and the entry, the practice
# profile, the gate, the run and the read are all made at it. `test-profile` runs only after the read's PASSED AT
# PRACTICE is committed (its reading tracked at HEAD, row H36's outcome cell holding the statement word for word):
# COMMIT is then that later commit (G1 before G0: E10.5's profile is read only after the practice read's outcome).
# RUN names this call's own scratch folder ~/v/scratch/h36/<RUN>, which may not exist: a run folder is never reused and
# nothing is ever deleted here.
# NO STEP TWICE. Each of entry, profile, gate, run, read and test-profile is started once, ever, at any commit and
# under any RUN: a step is started when this launcher's refusals and the script's own (asked first: `refusals`) have
# passed, and from then on its exit file ~/v/scratch/h36/<STEP>.exit exists (the word running, then the step's code).
# A call that finds that file refuses the step (88). A step that ends without a verdict of its own (a killed process,
# a full disk, a crash) is never started again: the outcome is then VOID, which the main session states.
#   check  HEAD is COMMIT, the tracked tree is clean, the script is here, row H36 is in TRIES.md and well formed
#          (experiments/harness/h36_row.py check: one row, eleven cells, its first cell the fixed wording), the phase's
#          output does not exist yet, the virtual environment and the T1 and T2 stages are where they are said to be
#   tests  the unit tests of the script, its table, its row module and this launcher (synthetic data only)
#   entry, profile, gate, run, read   (phase read) as the script's docstring says; the run step is five processes of
#          two seeds, each refused before anything loads unless the gate of COMMIT said PASS; the eleven tries count
#          from the run step's start
#   test-profile  (phase test-profile) the test profile; PASS prints the statement with the profile's sha256
# <results> is experiments/harness/results/h36_heart_run_profile: every record and the parts are written there, to be
# committed with the outcome; T's profile is experiments/heart_profile/results/xpc_profile.json.
# Environment: H36_VENV names the virtual environment (default $PWD/.venv; exported as UV_PROJECT_ENVIRONMENT; the
# script puts its own checkout's src first); H36_DATA names the data folder (default data/raw of this checkout).
# 8 BLAS threads a process, NUMPY_MADVISE_HUGEPAGE=0, TMPDIR under the run's scratch folder. Nothing here deletes or
# sends anything.
# Exit codes of the launcher's own refusals: 88 the step was started before, 89 the script refuses the step, 90 odd
# arguments, 91 the run folder exists, 92 parts exist that this launcher did not make, 93 row H36 missing or
# malformed, 94 not the checkout's root, 95 the phase's output exists, 96 HEAD is not COMMIT, 97 no T1 or T2 stages in
# the data folder, 98 tracked tree not clean, 99 no virtual environment. Otherwise the failed step's own code: tests
# (pytest's), entry and profile (3 STOP), gate (3 FAIL), run (1: a part failed or wrote no file), read, test-profile
# (3 STOP).
set -u
COMMIT=${1:?usage: run_h36.sh COMMIT RUN PHASE}
RUN=${2:?usage: run_h36.sh COMMIT RUN PHASE}
PHASE=${3:?usage: run_h36.sh COMMIT RUN PHASE}
[ "$#" = 3 ] || { echo "refused: COMMIT, RUN and PHASE only" >&2; exit 90; }
case "$RUN" in
  *[!A-Za-z0-9_-]*) echo "refused: RUN takes letters, digits, _ and - only" >&2; exit 90 ;;
esac
case "$COMMIT" in
  *[!0-9a-f]*) echo "refused: COMMIT must be a full hash" >&2; exit 90 ;;
esac
[ "${#COMMIT}" = 40 ] || { echo "refused: COMMIT must be a full hash" >&2; exit 90; }
case "$PHASE" in
  read) STEPS="entry profile gate run read" ;;
  test-profile) STEPS="test-profile" ;;
  *) echo "refused: PHASE is read or test-profile" >&2; exit 90 ;;
esac
export PATH=$HOME/.local/bin:$PATH
BASE=$HOME/v/scratch/h36
S=$BASE/$RUN
mkdir -p "$BASE"
mkdir "$S" 2> /dev/null || { echo "refused: the run folder $S exists (a run folder is never reused)" >&2; exit 91; }
mkdir "$S/tmp"
echo running > "$S/done.exit"
STEP=check
OWNED=none  # the step this call has started: only then is its exit file under $BASE this call's to write
finish() {
  echo "$1" > "$S/$STEP.exit"
  if [ "$OWNED" = "$STEP" ]; then echo "$1" > "$BASE/$STEP.exit"; fi
  echo "$1" > "$S/done.exit"
  exit "$1"
}
begin() { echo running > "$S/$STEP.exit"; }
step_ok() {
  echo 0 > "$S/$STEP.exit"
  if [ "$OWNED" = "$STEP" ]; then echo 0 > "$BASE/$STEP.exit"; fi
}
begin

SCRIPT=experiments/harness/h36_heart_run_profile.py
OUT=experiments/harness/results/h36_heart_run_profile
head=$(git rev-parse HEAD)
[ "$head" = "$COMMIT" ] || { echo "refused: HEAD is not $COMMIT" >&2; finish 96; }
if [ -n "$(git status --porcelain --untracked-files=no)" ]; then
  echo "refused: tracked files differ from the commit" >&2
  finish 98
fi
[ -f "$SCRIPT" ] || { echo "refused: run from the root of the checkout ($SCRIPT is not here)" >&2; finish 94; }
python3 -I -B experiments/harness/h36_row.py check experiments/harness/TRIES.md \
  || { echo "refused: row H36 is missing from experiments/harness/TRIES.md or malformed" >&2; finish 93; }
if [ "$PHASE" = read ] && { [ -e "$OUT/reading.json" ] || [ -e "$OUT/reading.txt" ]; }; then
  echo "refused: $OUT already holds a reading" >&2
  finish 95
fi
if [ "$PHASE" = test-profile ] && { [ -e "$OUT/test_profile" ] || [ -e experiments/heart_profile/results/xpc_profile.json ]; }; then
  echo "refused: the test profile exists" >&2
  finish 95
fi
VENV=${H36_VENV:-$PWD/.venv}
[ -x "$VENV/bin/python" ] \
  || { echo "refused: no virtual environment at $VENV (name one with H36_VENV)" >&2; finish 99; }
DATA=${H36_DATA:-data/raw}
{ [ -d "$DATA/T2" ] && [ -d "$DATA/T1" ]; } \
  || { echo "refused: no T1 and T2 stages under $DATA (name the data folder with H36_DATA)" >&2; finish 97; }
echo "$head $PHASE" > "$S/head"
step_ok

export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 NUMPY_MADVISE_HUGEPAGE=0 TMPDIR="$S/tmp"
export UV_PROJECT_ENVIRONMENT="$VENV"
PY="uv run --frozen --no-sync python -I -B"
DATA_ARGS=()
if [ -n "${H36_DATA:-}" ]; then DATA_ARGS=(--data "$H36_DATA"); fi
PARTS="0-1 2-3 4-5 6-7 8-9"

# Each step is started once, ever: refused (88) when its exit file under $BASE exists, whatever it holds. That file
# is made, and must be new, only when this launcher's refusals and the script's own for the step have passed.
start_once() {
  if [ -e "$BASE/$STEP.exit" ]; then
    echo "refused: the $STEP step was started before ($BASE/$STEP.exit exists); it is not started twice" >&2
    finish 88
  fi
  if [ "$STEP" = run ] && [ -e "$OUT/parts" ]; then
    echo "refused: $OUT/parts exists, made outside this launcher; the run step is not started" >&2
    finish 92
  fi
  $PY "$SCRIPT" refusals --step "$STEP" --commit "$COMMIT" > "$S/$STEP.refusals" 2>&1 \
    || { echo "refused: the script refuses the $STEP step (its words are in $S/$STEP.refusals)" >&2; finish 89; }
  ( set -o noclobber; echo running > "$BASE/$STEP.exit" ) 2> /dev/null \
    || { echo "refused: the $STEP step was started before ($BASE/$STEP.exit exists); it is not started twice" >&2
         finish 88; }
  OWNED=$STEP
  echo "$COMMIT $RUN" > "$BASE/$STEP.head"
}

STEP=tests
begin
$PY -m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp="$S/tmp/pytest" \
  tests/test_xpc_table.py tests/test_h36_heart_run_profile.py > "$S/tests.log" 2>&1 || finish $?
step_ok

for STEP in $STEPS; do
  begin
  start_once
  if [ "$STEP" = run ]; then
    mkdir -p "$S/logs" "$S/exit"
    pids=()
    for seeds in $PARTS; do
      (
        $PY "$SCRIPT" run --seeds "$seeds" --commit "$COMMIT" ${DATA_ARGS[@]+"${DATA_ARGS[@]}"} \
          > "$S/logs/$seeds.log" 2>&1
        echo $? > "$S/exit/$seeds"
      ) &
      pids+=($!)
    done
    for p in "${pids[@]}"; do wait "$p"; done
    bad=0
    for seeds in $PARTS; do
      [ "$(cat "$S/exit/$seeds" 2>/dev/null)" = 0 ] || bad=1
      [ -f "$OUT/parts/$seeds/h36.json" ] || bad=1
    done
    [ "$bad" = 0 ] || finish 1
  else
    $PY "$SCRIPT" "$STEP" --commit "$COMMIT" ${DATA_ARGS[@]+"${DATA_ARGS[@]}"} > "$S/$STEP.out" 2>&1 \
      || finish $?
  fi
  step_ok
done

echo 0 > "$S/done.exit"
exit 0
