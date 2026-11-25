#!/usr/bin/env bash
# Row H38 (experiments/harness/TRIES.md): the counted read of heart file OIC125
# (experiments/harness/h38_heart_insitu.py), on a DGX campaign checkout or a detached worktree of it (never the
# reference checkout). Run from the root of a clean checkout or worktree at COMMIT, detached, and wait on done.exit
# with a background poller; never a foreground sleep. CPU only.
#   bash experiments/harness/run_h38.sh <COMMIT> <RUN> <FETCHED>
# FETCHED is the fetch's folder: it must pass fetch_oic.py --check (MANIFEST.json against the pinned sizes and the
# files on disk), and the two practice chord files must name its files at the sha256 it records.
# COMMIT is the full hash HEAD must be: the run commit, which holds row H38 with G0's PASS statement in its outcome cell
# (and no read outcome) and the three committed chord files of the sha256 that statement records. RUN names this call's
# own scratch folder ~/v/scratch/h38/<RUN>, which may not exist: a run folder is never reused and nothing is deleted.
# NO STEP TWICE. Each of gate, run and read is started once, ever, at any commit and under any RUN: a step is started
# when this launcher's refusals and the script's own (asked first: `refusals`) have passed, and from then on its exit
# file ~/v/scratch/h38/<STEP>.exit exists (the word running, then the step's code). A call that finds that file refuses
# the step (88). A step that ends without a verdict of its own (a killed process, a full disk) is never started again:
# the outcome is then VOID, which the main session states. Steps, in order, each writing <S>/<STEP>.exit (running
# first, then its code), the call ending at the first that is not 0:
#   check  HEAD is COMMIT, the tracked tree is clean, the script is here, row H38's outcome cell holds G0's PASS
#          statement and no read outcome (experiments/harness/h38_row.py g0), no reading exists, the virtual
#          environment and the T2 stages are where they are said to be
#   tests  the script's, G0's, the in-situ form's and these launchers' tests (synthetic data only)
#   gate   counting nothing: C rebuilt to row H35's committed parts, the in-situ step with its chord out (FAIL: VOID)
#   run    ten processes, the two backtests by five pairs of seeds; the ten tries count from this step's start
#   read   the void conditions, the bars, the outcome line
# Results under experiments/harness/results/h38_heart_insitu (gate/, parts/<task>_<seeds>/h38.json, reading.json and
# reading.txt), to be committed with the outcome. H38_VENV names the virtual environment (default $PWD/.venv; exported
# as UV_PROJECT_ENVIRONMENT; the script puts its own checkout's src first); H38_DATA names the data folder (default
# data/raw). 8 BLAS threads a process, NUMPY_MADVISE_HUGEPAGE=0, TMPDIR under the run's scratch folder. Nothing here
# deletes or sends anything.
# Exit codes of the launcher's own refusals: 87 FETCHED or the chords fail the manifest check, 88 the step was
# started before, 89 the script refuses the step, 90 odd
# arguments, 91 the run folder exists, 92 parts exist that this launcher did not make, 93 row H38 lacks G0's PASS
# statement or holds a read outcome, 94 not the checkout's root, 95 a reading exists, 96 HEAD is not COMMIT, 97 no T2
# stages in the data folder, 98 tracked tree not clean, 99 no virtual environment. Otherwise the failed step's own code:
# tests (pytest's), gate (3 FAIL), run (1: a part failed or wrote no file), read.
set -u
COMMIT=${1:?usage: run_h38.sh COMMIT RUN FETCHED}
RUN=${2:?usage: run_h38.sh COMMIT RUN FETCHED}
FETCHED=${3:?usage: run_h38.sh COMMIT RUN FETCHED}
[ "$#" = 3 ] || { echo "refused: COMMIT, RUN and FETCHED only" >&2; exit 90; }
case "$RUN" in
  *[!A-Za-z0-9_-]*) echo "refused: RUN takes letters, digits, _ and - only" >&2; exit 90 ;;
esac
case "$COMMIT" in
  *[!0-9a-f]*) echo "refused: COMMIT must be a full hash" >&2; exit 90 ;;
esac
[ "${#COMMIT}" = 40 ] || { echo "refused: COMMIT must be a full hash" >&2; exit 90; }
export PATH=$HOME/.local/bin:$PATH
BASE=$HOME/v/scratch/h38
S=$BASE/$RUN
mkdir -p "$BASE"
mkdir "$S" 2> /dev/null || { echo "refused: the run folder $S exists (a run folder is never reused)" >&2; exit 91; }
mkdir "$S/tmp"
echo running > "$S/done.exit"
STEP=check
OWNED=none
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

SCRIPT=experiments/harness/h38_heart_insitu.py
OUT=experiments/harness/results/h38_heart_insitu
head=$(git rev-parse HEAD)
[ "$head" = "$COMMIT" ] || { echo "refused: HEAD is not $COMMIT" >&2; finish 96; }
if [ -n "$(git status --porcelain --untracked-files=no)" ]; then
  echo "refused: tracked files differ from the commit" >&2
  finish 98
fi
[ -f "$SCRIPT" ] || { echo "refused: run from the root of the checkout ($SCRIPT is not here)" >&2; finish 94; }
python3 -I -B experiments/harness/h38_row.py g0 experiments/harness/TRIES.md > "$S/chords" 2> "$S/row.err" \
  || { echo "refused: row H38 lacks G0's PASS statement or holds a read outcome (see $S/row.err)" >&2; finish 93; }
python3 -I -B experiments/heart_ext/fetch_oic.py --check "$FETCHED" --chords \
  experiments/heart_ext/results/chords_E9.5.json experiments/heart_ext/results/chords_E8.75.json \
  > "$S/manifest.out" 2>&1 \
  || { echo "refused: $FETCHED or the practice chords fail fetch_oic.py --check (see $S/manifest.out)" >&2; finish 87; }
if [ -e "$OUT/reading.json" ] || [ -e "$OUT/reading.txt" ]; then
  echo "refused: $OUT already holds a reading" >&2
  finish 95
fi
VENV=${H38_VENV:-$PWD/.venv}
[ -x "$VENV/bin/python" ] \
  || { echo "refused: no virtual environment at $VENV (name one with H38_VENV)" >&2; finish 99; }
DATA=${H38_DATA:-data/raw}
[ -d "$DATA/T2" ] || { echo "refused: no T2 stages under $DATA (name the data folder with H38_DATA)" >&2; finish 97; }
echo "$head" > "$S/head"
step_ok

export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 NUMPY_MADVISE_HUGEPAGE=0 TMPDIR="$S/tmp"
export UV_PROJECT_ENVIRONMENT="$VENV"
PY="uv run --frozen --no-sync python -I -B"
KEYS=("T2h/<=E9.5->E10.5" "T2h/<=E8.75->E10.5")
SLUGS=(e95 e875)
PARTS="0-1 2-3 4-5 6-7 8-9"

start_once() {
  if [ -e "$BASE/$STEP.exit" ]; then
    echo "refused: the $STEP step was started before ($BASE/$STEP.exit exists); it is not started twice" >&2
    finish 88
  fi
  if [ "$STEP" = run ] && [ -e "$OUT/parts" ]; then
    echo "refused: $OUT/parts exists, made outside this launcher; the run step is not started" >&2
    finish 92
  fi
  case "$STEP" in
    gate) ask=(--out "$OUT/gate") ;;
    run) ask=(--gate "$OUT/gate") ;;
    *) ask=(--out "$OUT") ;;
  esac
  $PY "$SCRIPT" refusals --step "$STEP" --commit "$COMMIT" "${ask[@]}" > "$S/$STEP.refusals" 2>&1 \
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
  tests/test_h38_heart_insitu.py tests/test_heart_oic_build.py tests/test_heart_ext_oic.py \
  tests/test_heart_ext_fetch_oic.py tests/test_oic_launchers.py > "$S/tests.log" 2>&1 || finish $?
step_ok

STEP=gate
begin
start_once
$PY "$SCRIPT" gate --out "$OUT/gate" --commit "$COMMIT" --data "$DATA" > "$S/gate.out" 2>&1 || finish $?
step_ok

STEP=run
begin
start_once
mkdir -p "$S/logs" "$S/exit"
pids=()
for i in 0 1; do
  for seeds in $PARTS; do
    (
      $PY "$SCRIPT" run --task "${KEYS[$i]}" --seeds "$seeds" --out "$OUT/parts/${SLUGS[$i]}_$seeds" \
        --gate "$OUT/gate" --commit "$COMMIT" --data "$DATA" > "$S/logs/${SLUGS[$i]}_$seeds.log" 2>&1
      echo $? > "$S/exit/${SLUGS[$i]}_$seeds"
    ) &
    pids+=($!)
  done
done
for p in "${pids[@]}"; do wait "$p"; done
bad=0
for i in 0 1; do
  for seeds in $PARTS; do
    [ "$(cat "$S/exit/${SLUGS[$i]}_$seeds" 2> /dev/null)" = 0 ] || bad=1
    [ -f "$OUT/parts/${SLUGS[$i]}_$seeds/h38.json" ] || bad=1
  done
done
[ "$bad" = 0 ] || finish 1
step_ok

STEP=read
begin
start_once
$PY "$SCRIPT" read --runs "$OUT/parts" --gate "$OUT/gate" --out "$OUT" --commit "$COMMIT" > "$S/read.out" 2>&1 \
  || finish $?
step_ok

echo 0 > "$S/done.exit"
exit 0
