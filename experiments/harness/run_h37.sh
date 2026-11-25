#!/usr/bin/env bash
# Row H37 (experiments/harness/TRIES.md): the gate of row R69's file WBH, the heart's world-blur hedge
# (experiments/harness/h37_heart_world_blur.py), in one call, on a DGX campaign checkout or a detached worktree of it
# (never the reference checkout). Run from the root of a clean checkout or worktree at COMMIT, detached, and wait on
# done.exit with a background poller; never a foreground sleep. CPU only.
#   bash experiments/harness/run_h37.sh <COMMIT> <RUN>
# COMMIT is the run commit (the full hash): it holds row H37, with the gate's fixed first cell and no outcome yet.
# RUN names this call's own scratch folder ~/v/scratch/h37/<RUN>, which may not exist: a run folder is never reused and
# nothing is ever deleted here, so another call takes another name.
# NO STEP TWICE. G0, the gate, the run and the read are each started once, ever, at any commit and under any RUN: a
# step is started when this launcher's refusals and the script's own (asked first: `refusals`) have passed, and from
# then on its exit file ~/v/scratch/h37/<STEP>.exit exists (the word running, then the step's exit code; <STEP>.head
# names the commit and the RUN). A call that finds that file refuses the step (88) and stops. A step that ends without
# a verdict of its own (a killed process, a full disk, a crash) is never started again: the outcome is then VOID, which
# the main session states.
# In order, each step only after the one before it exited 0 in this call, each first writing the word running to
# ~/v/scratch/h37/<RUN>/<STEP>.exit, then its exit code. ~/v/scratch/h37/<RUN>/done.exit holds the word running from
# the start and the call's exit code at its end.
#   check   HEAD is COMMIT, the tracked tree is clean, the script is here (the checkout's root), TRIES.md holds exactly
#           one row H37 of eleven cells whose first cell is the gate's fixed text and whose outcome cell holds no
#           outcome yet (experiments/harness/h37_row.py check), and one row R69 that is the captain's order naming the
#           file (h37_row.py ruling); no reading exists; the virtual environment and the T2 and T3 stages are there
#   tests   the unit tests of the blur, the script, its row module, both launchers and the file's builder, and rows
#           H35's and H34's (synthetic arrays, throwaway repositories), into tests.log; it reads no released data,
#           decides nothing and runs in every call
#   banked  the identity rebuild of the banked file (OTK125: build_test.py heart --entry organ_scope --trend-strength
#           1.25 --trend-form category --z-mode raw --captain-ruling R55, then halve.py at 1000001) into the run
#           folder; both files must give their recorded sha256 (experiments/builds/results/
#           hearth_8a8bc3a_h2.SHA256SUMS) and the build's stamp this commit and clean, or everything stops (20 the
#           file, 21 its half, 11 the build, 12 halve.py, 14 no stamp, 15 the stamp). Output to /dev/null: a test build
#           is read by its exit code, its sha256 and its stamp's commit and clean fields only
#   g0      target-free and not a try, into <results>/g0: exit 0 PASS, 3 STOP (G0 states the outcome itself)
#   gate    counting nothing, into <results>/gate: exit 0 PASS, 3 FAIL (VOID), 4 STOP (G0 on the practice outputs)
#   run     C, W, U and P per seed: per transition five processes of two seeds, the two transitions on E10.5 one after
#           the other, into <results>/parts/<short>/<seeds>/h37.json. The six tries count from this step's start. Any
#           log line other than the script's two kinds is counted, never printed
#   read    the one read into <results> (reading.json, reading.txt; never overwritten); its two printed lines, the
#           outcome statement and what follows, in read.out
# <results> is experiments/harness/results/h37_heart_world_blur: every record and the parts are written there, to be
# committed with the outcome, which the main session quotes, as printed, in row H37's outcome cell
# (experiments/builds/heart_wbh_build.sh and experiments/builds/world_blur.py read it there by code).
# Environment: H37_VENV names the virtual environment (default $PWD/.venv; exported as UV_PROJECT_ENVIRONMENT; the
# script puts its own checkout's src first); H37_DATA names the data folder (default: data/raw of this checkout).
# 8 BLAS threads a process (16 for the banked build, as the reference builds), NUMPY_MADVISE_HUGEPAGE=0, TMPDIR under
# the run's scratch folder. Nothing here deletes or sends anything.
# Exit codes of the launcher's own refusals: 88 the step was started before, 89 the script refuses the step, 90 odd
# RUN, COMMIT or arguments, 91 the run folder exists, 92 parts exist that this launcher did not make, 93 no row H37
# or one that is malformed, or no row R69 that orders the file, 94 not the checkout's root, 95 a reading exists or row
# H37's outcome cell already holds an outcome, 96 HEAD is not COMMIT, 97 no T2 or T3 stages in the data folder, 98
# tracked tree not clean, 99 no virtual environment. Otherwise the failed step's own code: tests (pytest's), banked (as
# above), g0 (3 STOP), gate (3 FAIL, 4 STOP), run (1: a part failed or wrote no file), read.
set -u
set -o pipefail
COMMIT=${1:?usage: run_h37.sh COMMIT RUN}
RUN=${2:?usage: run_h37.sh COMMIT RUN}
[ "$#" = 2 ] || { echo "refused: COMMIT and RUN only" >&2; exit 90; }
case "$RUN" in
  *[!A-Za-z0-9_-]*) echo "refused: RUN takes letters, digits, _ and - only" >&2; exit 90 ;;
esac
case "$COMMIT" in
  *[!0-9a-f]*) echo "refused: COMMIT must be a full hash" >&2; exit 90 ;;
esac
[ "${#COMMIT}" = 40 ] || { echo "refused: COMMIT must be a full hash" >&2; exit 90; }
export PATH=$HOME/.local/bin:$PATH
BASE=$HOME/v/scratch/h37
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

SCRIPT=experiments/harness/h37_heart_world_blur.py
ROWS=experiments/harness/h37_row.py
OUT=experiments/harness/results/h37_heart_world_blur
head=$(git rev-parse HEAD)
[ "$head" = "$COMMIT" ] || { echo "refused: HEAD is not $COMMIT (give the full hash)" >&2; finish 96; }
if [ -n "$(git status --porcelain --untracked-files=no)" ]; then
  echo "refused: tracked files differ from the commit" >&2
  finish 98
fi
{ [ -f "$SCRIPT" ] && [ -f "$ROWS" ]; } \
  || { echo "refused: run from the root of the checkout ($SCRIPT is not here)" >&2; finish 94; }
python3 -I -B "$ROWS" check experiments/harness/TRIES.md 2> "$S/check.row"
code=$?
case "$code" in
  0) ;;
  4) echo "refused: row H37's outcome cell already holds an outcome (the gate is read once)" >&2; finish 95 ;;
  *) echo "refused: TRIES.md must hold exactly one row H37 of eleven cells whose first cell is the fixed text" >&2
     finish 93 ;;
esac
python3 -I -B "$ROWS" ruling experiments/harness/TRIES.md \
  || { echo "refused: TRIES.md does not hold exactly one row R69 that is the captain's order naming file WBH" >&2
       finish 93; }
if [ -e "$OUT/reading.json" ] || [ -e "$OUT/reading.txt" ]; then
  echo "refused: $OUT already holds a reading" >&2
  finish 95
fi
VENV=${H37_VENV:-$PWD/.venv}
[ -x "$VENV/bin/python" ] \
  || { echo "refused: no virtual environment at $VENV (name one with H37_VENV)" >&2; finish 99; }
DATA=${H37_DATA:-data/raw}
{ [ -d "$DATA/T2" ] && [ -d "$DATA/T3" ]; } \
  || { echo "refused: no T2 and T3 stages under $DATA (name the data folder with H37_DATA)" >&2; finish 97; }
echo "$head" > "$S/head"
step_ok

export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 NUMPY_MADVISE_HUGEPAGE=0 TMPDIR="$S/tmp"
export UV_PROJECT_ENVIRONMENT="$VENV"
PY="uv run --frozen --no-sync python -I -B"
DATA_ARGS=()
if [ -n "${H37_DATA:-}" ]; then DATA_ARGS=(--data "$H37_DATA"); fi
TASKS=("T2h/<=E9.5->E10.5:e95" "T2h/<=E8.75->E10.5:e875")
PARTS="0-1 2-3 4-5 6-7 8-9"
DONE_LINE='^T2h/[^ ]+ seed [0-9]+ \([0-9.]+ s\) done$'
WROTE_LINE='^wrote [^ ]+/h37\.json$'
G0_LINE='^G0 (PASS|STOP)$'
GATE_LINE='^GATE (PASS|FAIL|STOP)$'
OUTCOME_LINE='^(OUTCOME H37: .*\.|A hedge.s gate on one truth, not evidence of development\. .*)$'
BANKED_SHA=d4d8d75f8bed48813bf143f5d7726828e5c30cdb27bd708b6fc2d9a455d21893       # hearth h2: OTK125's file
BANKED_HALF_SHA=9835fa5d47d7021760beb2056ce44c2c846ded6fc31dacf64bbb2d9f0f7d1b54  # hearth h2: OTK125's half
BANKED_ARGS="--entry organ_scope --trend-strength 1.25 --trend-form category --z-mode raw --captain-ruling R55"
BANKED=$S/banked/T2heart_organ_ts125_tk_half.h5ad

# G0, the gate, the run and the read: each started once, ever. "$@" are the step's arguments to `refusals`.
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
  tests/test_world_blur.py tests/test_h37_heart_world_blur.py tests/test_run_h37_launcher.py \
  tests/test_heart_wbh_build.py tests/test_h35_heart_trend_forms.py tests/test_h34_heart_organ_trend.py \
  tests/test_halve.py \
  > "$S/tests.log" 2>&1 || finish $?
step_ok

STEP=banked
begin
mkdir -p "$S/banked"
(
  export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16
  $PY experiments/builds/build_test.py heart $BANKED_ARGS --out "$S/banked/T2heart_organ_ts125_tk.h5ad" \
    > /dev/null 2>&1 || exit 11
  $PY experiments/builds/halve.py "$S/banked/T2heart_organ_ts125_tk.h5ad" "$BANKED" > /dev/null 2>&1 || exit 12
) || finish $?
got=$(sha256sum "$S/banked/T2heart_organ_ts125_tk.h5ad" | cut -d' ' -f1)
[ "$got" = "$BANKED_SHA" ] || { echo "stopped: OTK125's command does not give its recorded file" >&2; finish 20; }
got=$(sha256sum "$BANKED" | cut -d' ' -f1)
[ "$got" = "$BANKED_HALF_SHA" ] || { echo "stopped: OTK125's half is not its recorded half" >&2; finish 21; }
st=$($PY -c 'import json,sys; s=json.load(open(sys.argv[1])); print(s["commit"], s["clean"])' \
  "$S/banked/T2heart_organ_ts125_tk.h5ad.stamp.json" 2> /dev/null) || finish 14
echo "banked $st" > "$S/stamps"
[ "$st" = "$head True" ] || finish 15
step_ok

STEP=g0
begin
start_once --out "$OUT/g0" --banked "$BANKED"
$PY "$SCRIPT" g0 --out "$OUT/g0" --banked "$BANKED" --commit "$COMMIT" ${DATA_ARGS[@]+"${DATA_ARGS[@]}"} \
  > "$S/g0.out" 2>&1
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
    [ -f "$OUT/parts/$short/$seeds/h37.json" ] || bad=1
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
