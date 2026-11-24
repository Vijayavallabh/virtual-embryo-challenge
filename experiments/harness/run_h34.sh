#!/usr/bin/env bash
# Row H34 (experiments/harness/TRIES.md): the practice read the heart's fourth search line requires before any build
# (row H28's trend arms on row H29b's D_f form; experiments/harness/h34_heart_organ_trend.py), in one call, on a DGX
# campaign checkout or a detached worktree of it (never the reference checkout). Run from the root of a clean checkout
# or worktree at COMMIT, detached, and wait on done.exit with a background poller; never a foreground sleep. CPU only.
#   bash experiments/harness/run_h34.sh <COMMIT> <RUN>
# COMMIT is the run commit (the full hash): the commit at which the gate is started; it holds row H34. RUN names this
# call's own scratch folder ~/v/scratch/h34/<RUN>, which may not exist: a run folder is never reused and nothing is
# ever deleted here, so another call takes another name.
# NO STEP TWICE. The gate, the run and the read are each started once, ever, at any commit and under any RUN: a step is
# started when this launcher's refusals and the script's own (asked first: `refusals`) have passed, and from then on
# its exit file ~/v/scratch/h34/<STEP>.exit exists (the word running, then the step's exit code; <STEP>.head names the
# commit and the RUN). A later call never starts that step again: if the file holds 0 and the step was started at
# COMMIT, the step is done and the call goes on to the next one (so a call that was refused, or cut short between two
# steps, can be followed by another at the same commit); otherwise the call refuses the step (88) and stops. A gate,
# run or read that ends without a verdict of its own (a killed process, a full disk, a crash: its exit file keeps the
# word running, or holds a code other than 0 or the gate's 3) is therefore never started again: by row H34 the outcome
# is then VOID, which the main session states. The tests step reads no released data and decides nothing: it runs in
# every call, and a commit at which it does not exit 0 is not the run commit.
# In order, each step only after the one before it exited 0 at COMMIT, each first writing the word running to
# ~/v/scratch/h34/<RUN>/<STEP>.exit, then its exit code (or the word before: done at COMMIT in an earlier call). The
# call stops at the first step that does not exit 0. ~/v/scratch/h34/<RUN>/done.exit holds the word running from the
# start and the call's exit code at its end: 0 only when every step exited 0 here or before, otherwise the failed
# step's code.
#   check  HEAD is COMMIT, the tracked tree is clean, the script is here (the checkout's root), TRIES.md holds a row
#          H34, no reading exists yet, the virtual environment and the T2 stages are where they are said to be; then
#          it writes the head file
#   tests  the unit tests of the script and of this launcher, and rows H29b's and H27's (synthetic arrays, hand-made
#          tables, throwaway repositories and committed result files read for arithmetic), into tests.log
#   gate   first, counting nothing: the control D_f alone built and scored; G1, its record, its draw, its eighteen
#          means and its paired statistics against row H29b's committed reading.json, exactly; G2, the trend step's
#          mapping and damping at 1.5 against 1.5 times those at 1.0, exactly. One process, both transitions, ten
#          seeds; into <results>/gate (gate.json, gate.txt); exit 0 PASS, 3 FAIL (the gate states the VOID outcome
#          itself and nothing further runs), else it ended without a verdict. Its printed lines are in gate.out
#   run    the four counted arms and the control once more: per transition five processes of two seeds, the two
#          transitions on E10.5 one after the other, into <results>/parts/<short>/<seeds>/h34.json. Each part refuses
#          before anything loads unless row H34 is in TRIES.md, the tracked tree is clean, HEAD is COMMIT, the gate of
#          COMMIT said PASS and its output does not exist. The eight tries count from this step's start. The script
#          prints "<task> seed <s> (<sec> s) done" and "wrote <path>" only; any other log line is counted, never
#          printed (exit/*.other_lines and exit/other_lines in the run folder)
#   read   the one read (the void conditions first, then the stop rule) into <results> (reading.json, reading.txt;
#          never overwritten); its two printed lines, the outcome statement and what follows, in read.out
# <results> is experiments/harness/results/h34_heart_organ_trend: the gate's record, the parts and the read's tables
# are written there, to be committed with the outcome.
# Environment: H34_VENV names the virtual environment (default $PWD/.venv; exported as UV_PROJECT_ENVIRONMENT; the
# script puts its own checkout's src first, so an environment made in another checkout serves); H34_DATA names the
# data folder (default: the harness's own, data/raw of this checkout).
# 8 BLAS threads a process, NUMPY_MADVISE_HUGEPAGE=0, TMPDIR under the run's scratch folder. Nothing here deletes or
# sends anything.
# Exit codes of the launcher's own refusals: 88 the step was started before, 89 the script refuses the step (the step
# is then not started and counts nothing), 90 odd RUN, 91 the run folder exists, 92 parts exist that this launcher did
# not make, 93 no row H34 in TRIES.md, 94 not the checkout's root, 95 a reading exists, 96 HEAD is not COMMIT, 97 no
# T2 stages in the data folder, 98 tracked tree not clean, 99 no virtual environment. Otherwise the failed step's own
# code: tests (pytest's), gate (3 FAIL; anything else no verdict), run (1: a part failed or wrote no file), read.
set -u
COMMIT=${1:?usage: run_h34.sh COMMIT RUN}
RUN=${2:?usage: run_h34.sh COMMIT RUN}
case "$RUN" in
  *[!A-Za-z0-9_-]*) echo "refused: RUN takes letters, digits, _ and - only" >&2; exit 90 ;;
esac
export PATH=$HOME/.local/bin:$PATH
BASE=$HOME/v/scratch/h34
S=$BASE/$RUN
mkdir -p "$BASE"
mkdir "$S" 2> /dev/null || { echo "refused: the run folder $S exists (a run folder is never reused)" >&2; exit 91; }
mkdir "$S/tmp"
echo running > "$S/done.exit"
STEP=check
OWNED=none  # the gate, run or read this call has started: only then is its exit file under $BASE this call's to write
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

SCRIPT=experiments/harness/h34_heart_organ_trend.py
OUT=experiments/harness/results/h34_heart_organ_trend
head=$(git rev-parse HEAD)
[ "$head" = "$COMMIT" ] || { echo "refused: HEAD is not $COMMIT (give the full hash)" >&2; finish 96; }
if [ -n "$(git status --porcelain --untracked-files=no)" ]; then
  echo "refused: tracked files differ from the commit" >&2
  finish 98
fi
[ -f "$SCRIPT" ] || { echo "refused: run from the root of the checkout ($SCRIPT is not here)" >&2; finish 94; }
grep -q '^| H34 ' experiments/harness/TRIES.md \
  || { echo "refused: row H34 is not in experiments/harness/TRIES.md" >&2; finish 93; }
if [ -e "$OUT/reading.json" ] || [ -e "$OUT/reading.txt" ]; then
  echo "refused: $OUT already holds a reading" >&2
  finish 95
fi
VENV=${H34_VENV:-$PWD/.venv}
[ -x "$VENV/bin/python" ] \
  || { echo "refused: no virtual environment at $VENV (name one with H34_VENV)" >&2; finish 99; }
DATA=${H34_DATA:-data/raw}
[ -d "$DATA/T2" ] \
  || { echo "refused: no T2 stages under $DATA (name the data folder with H34_DATA)" >&2; finish 97; }
echo "$head" > "$S/head"
step_ok

export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 NUMPY_MADVISE_HUGEPAGE=0 TMPDIR="$S/tmp"
export UV_PROJECT_ENVIRONMENT="$VENV"
PY="uv run --frozen --no-sync python -I -B"
DATA_ARGS=()
if [ -n "${H34_DATA:-}" ]; then DATA_ARGS=(--data "$H34_DATA"); fi
TASKS=("T2h/<=E9.5->E10.5:e95" "T2h/<=E8.75->E10.5:e875")
PARTS="0-1 2-3 4-5 6-7 8-9"
DONE_LINE='^T2h/[^ ]+ seed [0-9]+ \([0-9.]+ s\) done$'
WROTE_LINE='^wrote [^ ]+/h34\.json$'
GATE_LINE='^GATE (PASS|FAIL)$'
VOID_LINE='^(Try H34: VOID: .*; nothing is built\.|A screen, not evidence of development\. .*)$'

# The gate, the run and the read: each started once, ever. start_once starts the step and returns 0, or returns 1
# when the step wrote exit 0 at COMMIT in an earlier call (it is then not started again: the word before goes into
# this call's exit file for it and the call goes on to the next step). A step that was started before and did not
# write exit 0, or was started at another commit, is refused (88). The step's exit file under $BASE is made, and must
# be new, only when this launcher's refusals and the script's own for the step have passed; "$@" are the step's
# arguments to `refusals`.
start_once() {
  if [ -e "$BASE/$STEP.exit" ]; then
    if [ "$(cat "$BASE/$STEP.exit")" = 0 ] && [ "$(cut -d" " -f1 "$BASE/$STEP.head" 2> /dev/null)" = "$COMMIT" ]; then
      echo before > "$S/$STEP.exit"
      return 1
    fi
    echo "refused: the $STEP step was started before ($BASE/$STEP.exit exists); it is not started twice" >&2
    finish 88
  fi
  if [ "$STEP" = run ] && [ -e "$OUT/parts" ]; then
    echo "refused: $OUT/parts exists, made outside this launcher; the run step is not started" >&2
    finish 92
  fi
  $PY "$SCRIPT" refusals --step "$STEP" --commit "$COMMIT" "$@" > "$S/$STEP.refusals" 2>&1 \
    || { echo "refused: the script refuses the $STEP step (its words are in $S/$STEP.refusals)" >&2; finish 89; }
  ( set -o noclobber; echo running > "$BASE/$STEP.exit" ) 2> /dev/null \
    || { echo "refused: the $STEP step was started before ($BASE/$STEP.exit exists); it is not started twice" >&2
         finish 88; }
  OWNED=$STEP
  echo "$COMMIT $RUN" > "$BASE/$STEP.head"
  return 0
}

STEP=tests
begin
$PY -m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp="$S/tmp/pytest" \
  tests/test_h34_heart_organ_trend.py tests/test_h29b_heart_organ.py tests/test_h27_heart_forecast_size.py \
  > "$S/tests.log" 2>&1 || finish $?
step_ok

STEP=gate
begin
if start_once --out "$OUT/gate"; then
  $PY "$SCRIPT" gate --out "$OUT/gate" --commit "$COMMIT" ${DATA_ARGS[@]+"${DATA_ARGS[@]}"} > "$S/gate.out" 2>&1
  code=$?
  # imported code and Python warnings write to the same file: the lines other than the expected kinds are counted
  n=$(grep -cvE "$DONE_LINE|$GATE_LINE|$VOID_LINE" "$S/gate.out" 2>/dev/null)
  echo "${n:-unreadable}" > "$S/gate.other_lines"
  [ "$code" = 0 ] || finish "$code"
  step_ok
fi

STEP=run
begin
if start_once --gate "$OUT/gate"; then
  mkdir -p "$S/logs" "$S/exit"
  bad=0
  for entry in "${TASKS[@]}"; do
    key=${entry%:*}
    short=${entry##*:}
    pids=()
    for seeds in $PARTS; do
      (
        $PY "$SCRIPT" run --task "$key" --seeds "$seeds" --gate "$OUT/gate" --commit "$COMMIT" \
          --out "$OUT/parts/$short/$seeds" ${DATA_ARGS[@]+"${DATA_ARGS[@]}"} > "$S/logs/$short-$seeds.log" 2>&1
        echo $? > "$S/exit/$short-$seeds"
      ) &
      pids+=($!)
    done
    for p in "${pids[@]}"; do wait "$p"; done
    for seeds in $PARTS; do
      [ "$(cat "$S/exit/$short-$seeds" 2>/dev/null)" = 0 ] || bad=1
      [ -f "$OUT/parts/$short/$seeds/h34.json" ] || bad=1
      # imported code and Python warnings write to the same log: count the lines other than the two expected kinds
      # and keep only the count; the lines are never printed
      n=$(grep -cvE "$DONE_LINE|$WROTE_LINE" "$S/logs/$short-$seeds.log" 2>/dev/null)
      echo "${n:-unreadable}" > "$S/exit/$short-$seeds.other_lines"
    done
    [ "$bad" = 0 ] || break
  done
  cat "$S"/exit/*.other_lines \
    | awk '{ if ($1 ~ /^[0-9]+$/) s += $1; else u = 1 } END { print (u ? "unreadable" : s) }' > "$S/exit/other_lines"
  [ "$bad" = 0 ] || finish 1
  step_ok
fi

STEP=read
begin
if start_once --out "$OUT"; then
  $PY "$SCRIPT" read --runs "$OUT/parts" --gate "$OUT/gate" --commit "$COMMIT" --out "$OUT" > "$S/read.out" 2>&1 \
    || finish $?
  step_ok
fi

echo 0 > "$S/done.exit"
exit 0
