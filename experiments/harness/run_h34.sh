#!/usr/bin/env bash
# Row H34 (experiments/harness/TRIES.md): the practice read the heart's fourth search line requires before any build
# (row H28's trend arms on row H29b's D_f form; experiments/harness/h34_heart_organ_trend.py), in one call, on a DGX
# campaign checkout or a detached worktree of it (never the reference checkout). Run from the root of a clean checkout
# or worktree at COMMIT, detached, and wait on done.exit with a background poller; never a foreground sleep. CPU only.
#   bash experiments/harness/run_h34.sh <COMMIT> <RUN>
# COMMIT is row H34's run commit (the full hash). RUN names this run's own scratch folder ~/v/scratch/h34/<RUN> and its
# output folder runs/h34_<commit7>_<RUN>. Neither may exist: a run folder is never reused and nothing is ever deleted
# here, so a second run takes a new name.
# In order, each step first writing the word running to ~/v/scratch/h34/<RUN>/<STEP>.exit, then its exit code. The call
# stops at the first step that does not exit 0. ~/v/scratch/h34/<RUN>/done.exit holds the word running from the start
# and the call's exit code at its end: 0 only when every step exited 0, otherwise the failed step's code.
#   check  HEAD is COMMIT, the tracked tree is clean, the script is here (the checkout's root), TRIES.md holds a row
#          H34, no reading exists yet, the virtual environment and the T2 stages are where they are said to be, the
#          output folder is new; then it writes the head file
#   tests  the unit tests of the script and of this launcher, and rows H29b's and H27's (synthetic arrays, hand-made
#          tables and throwaway repositories only), into tests.log
#   gate   first, counting nothing: the control D_f alone built and scored, its eighteen plain means held against row
#          H29b's committed reading exactly (G1), and the trend step's log factors and damping at 1.5 against 1.5 times
#          those at 1.0, exactly (G2); one process, both transitions, ten seeds; into <output>/gate (gate.json,
#          gate.txt); exit 0 PASS, 3 FAIL (the outcome is VOID: nothing further runs), else an error. Its printed lines
#          are in gate.out
#   run    the counted arms beside a rebuilt control: per transition five processes of two seeds, the two transitions
#          on E10.5 one after the other, into <output>/<short>/<seeds>/h34.json. Each part refuses before anything
#          loads unless row H34 is in TRIES.md, the tracked tree is clean, HEAD is COMMIT, the gate of COMMIT said PASS
#          and its output does not exist. The script prints "<task> seed <s> (<sec> s) done" and "wrote <path>" only;
#          any other log line is counted, never printed (<output>/exit/*.other_lines and other_lines)
#   read   the one read (the void conditions first, then the stop rule) into experiments/harness/results/
#          h34_heart_organ_trend (reading.json, reading.txt; never overwritten); its two printed lines, the outcome and
#          what follows, in read.out
# Environment: H34_VENV names the virtual environment (default $PWD/.venv; exported as UV_PROJECT_ENVIRONMENT; the
# script puts its own checkout's src first, so an environment made in another checkout serves); H34_DATA names the
# data folder (default: the harness's own, data/raw of this checkout).
# 8 BLAS threads a process, NUMPY_MADVISE_HUGEPAGE=0, TMPDIR under the run's scratch folder. Nothing here deletes or
# sends anything.
# Exit codes of the launcher's own refusals: 90 odd RUN, 91 the run folder exists, 92 the output folder exists, 93 no
# row H34 in TRIES.md, 94 not the checkout's root, 95 a reading exists, 96 HEAD is not COMMIT, 97 no T2 stages in the
# data folder, 98 tracked tree not clean, 99 no virtual environment. Otherwise the failed step's own code: tests
# (pytest's), gate (3 FAIL; anything else an error), run (1: a part failed or wrote no file), read.
set -u
COMMIT=${1:?usage: run_h34.sh COMMIT RUN}
RUN=${2:?usage: run_h34.sh COMMIT RUN}
case "$RUN" in
  *[!A-Za-z0-9_-]*) echo "refused: RUN takes letters, digits, _ and - only" >&2; exit 90 ;;
esac
export PATH=$HOME/.local/bin:$PATH
S=$HOME/v/scratch/h34/$RUN
mkdir -p "$HOME/v/scratch/h34"
mkdir "$S" 2> /dev/null || { echo "refused: the run folder $S exists (a run folder is never reused)" >&2; exit 91; }
mkdir "$S/tmp"
echo running > "$S/done.exit"
STEP=check
finish() { echo "$1" > "$S/$STEP.exit"; echo "$1" > "$S/done.exit"; exit "$1"; }
begin() { echo running > "$S/$STEP.exit"; }
step_ok() { echo 0 > "$S/$STEP.exit"; }
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
R=runs/h34_${head:0:7}_$RUN
mkdir -p runs
mkdir "$R" 2> /dev/null || { echo "refused: the output folder $R exists" >&2; finish 92; }
echo "$head" > "$S/head"
echo "$R" > "$S/output"
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

STEP=tests
begin
$PY -m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp="$S/tmp/pytest" \
  tests/test_h34_heart_organ_trend.py tests/test_h29b_heart_organ.py tests/test_h27_heart_forecast_size.py \
  > "$S/tests.log" 2>&1 || finish $?
step_ok

STEP=gate
begin
$PY "$SCRIPT" gate --out "$R/gate" --commit "$COMMIT" ${DATA_ARGS[@]+"${DATA_ARGS[@]}"} > "$S/gate.out" 2>&1
code=$?
# imported code and Python warnings write to the same file: the lines other than the expected kinds are counted
n=$(grep -cvE "$DONE_LINE|$GATE_LINE|$VOID_LINE" "$S/gate.out" 2>/dev/null)
echo "${n:-unreadable}" > "$S/gate.other_lines"
[ "$code" = 0 ] || finish "$code"
step_ok

STEP=run
begin
mkdir -p "$R/logs" "$R/exit"
bad=0
for entry in "${TASKS[@]}"; do
  key=${entry%:*}
  short=${entry##*:}
  pids=()
  for seeds in $PARTS; do
    (
      $PY "$SCRIPT" run --task "$key" --seeds "$seeds" --gate "$R/gate" --commit "$COMMIT" \
        --out "$R/$short/$seeds" ${DATA_ARGS[@]+"${DATA_ARGS[@]}"} > "$R/logs/$short-$seeds.log" 2>&1
      echo $? > "$R/exit/$short-$seeds"
    ) &
    pids+=($!)
  done
  for p in "${pids[@]}"; do wait "$p"; done
  for seeds in $PARTS; do
    [ "$(cat "$R/exit/$short-$seeds" 2>/dev/null)" = 0 ] || bad=1
    [ -f "$R/$short/$seeds/h34.json" ] || bad=1
    # imported code and Python warnings write to the same log: count the lines other than the two expected kinds and
    # keep only the count; the lines are never printed
    n=$(grep -cvE "$DONE_LINE|$WROTE_LINE" "$R/logs/$short-$seeds.log" 2>/dev/null)
    echo "${n:-unreadable}" > "$R/exit/$short-$seeds.other_lines"
  done
  [ "$bad" = 0 ] || break
done
cat "$R"/exit/*.other_lines | awk '{ if ($1 ~ /^[0-9]+$/) s += $1; else u = 1 } END { print (u ? "unreadable" : s) }' \
  > "$R/exit/other_lines"
[ "$bad" = 0 ] || finish 1
step_ok

STEP=read
begin
$PY "$SCRIPT" read --runs "$R" --gate "$R/gate" --commit "$COMMIT" --out "$OUT" > "$S/read.out" 2>&1 || finish $?
step_ok

echo 0 > "$S/done.exit"
exit 0
