#!/usr/bin/env bash
# Row H35 (experiments/harness/TRIES.md): the counted screen of row R55's four trend forms
# (experiments/harness/h35_heart_trend_forms.py), in one call, on a DGX campaign checkout or a detached worktree of it
# (never the reference checkout). Run from the root of a clean checkout or worktree at COMMIT, detached, and wait on
# done.exit with a background poller; never a foreground sleep. CPU only.
#   bash experiments/harness/run_h35.sh <COMMIT> <RUN>
# COMMIT is the run commit (the full hash): it holds row H35, with the screen's fixed first cell and no outcome yet.
# RUN names this call's own scratch folder ~/v/scratch/h35/<RUN>, which may not exist: a run folder is never reused and
# nothing is ever deleted here, so another call takes another name.
# NO STEP TWICE. G0, the gate, the run and the read are each started once, ever, at any commit and under any RUN: a
# step is started when this launcher's refusals and the script's own (asked first: `refusals`) have passed, and from
# then on its exit file ~/v/scratch/h35/<STEP>.exit exists (the word running, then the step's exit code; <STEP>.head
# names the commit and the RUN). A call that finds that file refuses the step (88) and stops, whatever the file holds
# and whatever the commit; a step that was refused left no such file and counts nothing. A step that ends without a
# verdict of its own (a killed process, a full disk, a crash) is therefore never started again: by row H35 the outcome
# is then VOID, which the main session states.
# The tests step reads no released data and decides nothing: it runs in every call, and a commit at which it does not
# exit 0 is not the run commit.
# In order, each step only after the one before it exited 0 in this call, each first writing the word running to
# ~/v/scratch/h35/<RUN>/<STEP>.exit, then its exit code. The call stops at the first step that does not exit 0.
# ~/v/scratch/h35/<RUN>/done.exit holds the word running from the start and the call's exit code at its end.
#   check  HEAD is COMMIT, the tracked tree is clean, the script is here (the checkout's root), TRIES.md holds exactly
#          one row H35 of eleven cells whose first cell is the screen's fixed text and whose outcome cell holds no
#          outcome yet, no reading exists, the virtual environment and the T2 stages are where they are said to be;
#          then it writes the head file
#   tests  the unit tests of the script, of the trend forms and of both launchers, and rows H34's, H29b's and H27's
#          (synthetic arrays, hand-made tables, throwaway repositories and committed result files of rows H27 and
#          H29b read for arithmetic), into tests.log
#   g0     target-free and not a try: on the test inputs as the heart builder loads them for row R55's files, each
#          form against OTS125's step; into <results>/g0 (g0.json, g0.txt: each form's word and the condition that
#          set it, no figure); exit 0 when a form passes, 3 when none does (STOP: G0 states the outcome itself and
#          nothing further runs), 5 VOID (a code fault)
#   gate   first on the practice tasks, counting nothing: the copy D_f built and scored, G1 against row H29b's committed
#          reading, and the identities in every seed's context; into <results>/gate; exit 0 PASS, 3 FAIL (the gate
#          states the VOID outcome itself and nothing further runs)
#   run    the control and the forms G0 passed, and the copy D_f once more: per transition five processes of two
#          seeds, the two transitions on E10.5 one after the other, into <results>/parts/<short>/<seeds>/h35.json. The
#          tries count from this step's start. Any log line other than the script's two kinds is counted, never
#          printed (exit/*.other_lines and exit/other_lines in the run folder)
#   read   the one read (the void conditions first, then each form's word) into <results> (reading.json,
#          reading.txt; never overwritten); its two printed lines, the outcome statement and what follows, in read.out
# <results> is experiments/harness/results/h35_heart_trend_forms: G0's and the gate's records, the parts and the
# read's tables are written there, to be committed with the outcome, which the main session quotes, as printed, in row
# H35's outcome cell (experiments/builds/heart_trend_line_h.sh reads it there by code).
# Environment: H35_VENV names the virtual environment (default $PWD/.venv; exported as UV_PROJECT_ENVIRONMENT; the
# script puts its own checkout's src first); H35_DATA names the data folder (default: data/raw of this checkout).
# 8 BLAS threads a process, NUMPY_MADVISE_HUGEPAGE=0, TMPDIR under the run's scratch folder. Nothing here deletes or
# sends anything.
# Exit codes of the launcher's own refusals: 88 the step was started before, 89 the script refuses the step (the step
# is then not started and counts nothing), 90 odd RUN or arguments, 91 the run folder exists, 92 parts exist that this
# launcher did not make, 93 no row H35, or one that is malformed (not one row of eleven cells, or its first cell not
# the fixed text), 94 not the checkout's root, 95 a reading exists or row H35's outcome cell already holds an outcome,
# 96 HEAD is not COMMIT, 97 no T2 stages in the data folder, 98 tracked tree not clean, 99 no virtual environment.
# Otherwise the failed step's own code: tests (pytest's), g0 (3 STOP, 5 VOID, anything else no verdict), gate (3 FAIL;
# anything else no verdict), run (1: a part failed or wrote no file), read.
set -u
set -o pipefail
COMMIT=${1:?usage: run_h35.sh COMMIT RUN}
RUN=${2:?usage: run_h35.sh COMMIT RUN}
[ "$#" = 2 ] || { echo "refused: COMMIT and RUN only" >&2; exit 90; }
case "$RUN" in
  *[!A-Za-z0-9_-]*) echo "refused: RUN takes letters, digits, _ and - only" >&2; exit 90 ;;
esac
export PATH=$HOME/.local/bin:$PATH
BASE=$HOME/v/scratch/h35
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

SCRIPT=experiments/harness/h35_heart_trend_forms.py
OUT=experiments/harness/results/h35_heart_trend_forms
head=$(git rev-parse HEAD)
[ "$head" = "$COMMIT" ] || { echo "refused: HEAD is not $COMMIT (give the full hash)" >&2; finish 96; }
if [ -n "$(git status --porcelain --untracked-files=no)" ]; then
  echo "refused: tracked files differ from the commit" >&2
  finish 98
fi
[ -f "$SCRIPT" ] || { echo "refused: run from the root of the checkout ($SCRIPT is not here)" >&2; finish 94; }
python3 -I -B -c '
import sys
FIRST = "H35 (the counted screen of row R55\x27s trend forms: ten tries; G0 first, target-free, not a try)"
rows = [ln for ln in open(sys.argv[1], encoding="utf-8").read().splitlines() if ln.startswith("| H35 ")]
line = rows[0].rstrip() if len(rows) == 1 else ""
cells = [c.strip() for c in line[1:-1].split("|")] if line.endswith("|") else []
if len(cells) != 11 or cells[0] != FIRST:
    sys.exit(3)
sys.exit(4 if "OUTCOME H35:" in cells[10] else 0)
' experiments/harness/TRIES.md
code=$?
case "$code" in
  0) ;;
  4) echo "refused: row H35's outcome cell already holds an outcome (the screen is read once)" >&2; finish 95 ;;
  *) echo "refused: TRIES.md must hold exactly one row H35 of eleven cells whose first cell is the fixed text" >&2
     finish 93 ;;
esac
if [ -e "$OUT/reading.json" ] || [ -e "$OUT/reading.txt" ]; then
  echo "refused: $OUT already holds a reading" >&2
  finish 95
fi
VENV=${H35_VENV:-$PWD/.venv}
[ -x "$VENV/bin/python" ] \
  || { echo "refused: no virtual environment at $VENV (name one with H35_VENV)" >&2; finish 99; }
DATA=${H35_DATA:-data/raw}
[ -d "$DATA/T2" ] \
  || { echo "refused: no T2 stages under $DATA (name the data folder with H35_DATA)" >&2; finish 97; }
echo "$head" > "$S/head"
step_ok

export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 NUMPY_MADVISE_HUGEPAGE=0 TMPDIR="$S/tmp"
export UV_PROJECT_ENVIRONMENT="$VENV"
PY="uv run --frozen --no-sync python -I -B"
DATA_ARGS=()
if [ -n "${H35_DATA:-}" ]; then DATA_ARGS=(--data "$H35_DATA"); fi
TASKS=("T2h/<=E9.5->E10.5:e95" "T2h/<=E8.75->E10.5:e875")
PARTS="0-1 2-3 4-5 6-7 8-9"
DONE_LINE='^T2h/[^ ]+ seed [0-9]+ \([0-9.]+ s\) done$'
WROTE_LINE='^wrote [^ ]+/h35\.json$'
G0_LINE='^G0 ((category|contrast|whole|portion) )?(PASS|STOP|VOID|DROP|IDLE)$'
GATE_LINE='^GATE (PASS|FAIL)$'
OUTCOME_LINE='^(OUTCOME H35: .*\.|A screen, not evidence of development\. .*)$'

# G0, the gate, the run and the read: each started once, ever. A step whose exit file under $BASE exists is refused
# (88), whatever the file holds and at any commit. That file is made, and must be new, only when this launcher's
# refusals and the script's own for the step have passed: only then is the step started. "$@" are the step's
# arguments to `refusals`.
start_once() {
  if [ -e "$BASE/$STEP.exit" ]; then
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
}
count_other() {  # file pattern: the lines of the file that are not of the expected kinds, counted, never printed
  local n
  n=$(grep -cvE "$2" "$1" 2> /dev/null)
  echo "${n:-unreadable}"
}

STEP=tests
begin
$PY -m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp="$S/tmp/pytest" \
  tests/test_h35_heart_trend_forms.py tests/test_run_h35_launcher.py tests/test_heart_trend_line_h.py \
  tests/test_heart_trend_line_h_launcher.py tests/test_h34_heart_organ_trend.py tests/test_h29b_heart_organ.py \
  tests/test_h27_heart_forecast_size.py tests/test_b0_trend.py \
  > "$S/tests.log" 2>&1 || finish $?
step_ok

STEP=g0
begin
start_once --out "$OUT/g0"
$PY "$SCRIPT" g0 --out "$OUT/g0" --commit "$COMMIT" ${DATA_ARGS[@]+"${DATA_ARGS[@]}"} > "$S/g0.out" 2>&1
code=$?
count_other "$S/g0.out" "$G0_LINE|$OUTCOME_LINE" > "$S/g0.other_lines"
[ "$code" = 0 ] || finish "$code"
step_ok

STEP=gate
begin
start_once --out "$OUT/gate" --g0 "$OUT/g0"
$PY "$SCRIPT" gate --out "$OUT/gate" --g0 "$OUT/g0" --commit "$COMMIT" ${DATA_ARGS[@]+"${DATA_ARGS[@]}"} \
  > "$S/gate.out" 2>&1
code=$?
count_other "$S/gate.out" "$DONE_LINE|$GATE_LINE|$OUTCOME_LINE" > "$S/gate.other_lines"
[ "$code" = 0 ] || finish "$code"
step_ok

STEP=run
begin
start_once --g0 "$OUT/g0" --gate "$OUT/gate"
mkdir -p "$S/logs" "$S/exit"
bad=0
for entry in "${TASKS[@]}"; do
  key=${entry%:*}
  short=${entry##*:}
  pids=()
  for seeds in $PARTS; do
    (
      $PY "$SCRIPT" run --task "$key" --seeds "$seeds" --g0 "$OUT/g0" --gate "$OUT/gate" --commit "$COMMIT" \
        --out "$OUT/parts/$short/$seeds" ${DATA_ARGS[@]+"${DATA_ARGS[@]}"} > "$S/logs/$short-$seeds.log" 2>&1
      echo $? > "$S/exit/$short-$seeds"
    ) &
    pids+=($!)
  done
  for p in "${pids[@]}"; do wait "$p"; done
  for seeds in $PARTS; do
    [ "$(cat "$S/exit/$short-$seeds" 2> /dev/null)" = 0 ] || bad=1
    [ -f "$OUT/parts/$short/$seeds/h35.json" ] || bad=1
    count_other "$S/logs/$short-$seeds.log" "$DONE_LINE|$WROTE_LINE" > "$S/exit/$short-$seeds.other_lines"
  done
  [ "$bad" = 0 ] || break
done
cat "$S"/exit/*.other_lines | awk '{ if ($1 ~ /^[0-9]+$/) s += $1; else u = 1 } END { print (u ? "unreadable" : s) }' \
  > "$S/exit/other_lines"
[ "$bad" = 0 ] || finish 1
step_ok

STEP=read
begin
start_once --out "$OUT"
$PY "$SCRIPT" read --runs "$OUT/parts" --g0 "$OUT/g0" --gate "$OUT/gate" --commit "$COMMIT" --out "$OUT" \
  > "$S/read.out" 2>&1 || finish $?
step_ok

echo 0 > "$S/done.exit"
exit 0
