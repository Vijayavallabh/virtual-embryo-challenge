#!/usr/bin/env bash
# Row T72 (experiments/harness/TRIES.md): the practice read that binds row R54's late-mix file,
# experiments/t1_scope/t72_late_mix.py, one step per call, run from the root of the checkout or detached worktree it is
# called in (it moves to that checkout's top level first), with that folder's own environment and data (the E10.5
# file data/raw/T1/E10.5_RNA.h5ad). CPU only. Wait on the exit file with a background poller; never a foreground sleep.
# Every step takes ROWC, the full 40-character hash of the commit that holds row T72 (after row N64's shares file and
# its `shares.json sha256 <hex>` are committed) and this code, and runs only at HEAD = ROWC. Steps, in order; each is
# refused (exit 94) unless the step before it wrote exit 0 for the same commit:
#   tests ROWC   the targeted tests (synthetic arrays only)
#   runs ROWC    the twelve scenarios (cm own 0.30 0.33 lm16 0.45 0.60 x blood own lm16), started together, one
#                process each (2 workers, 4 BLAS threads a process); each scenario's code in runs/<scenario>.exit and
#                its log in runs/<scenario>.log (never printed: lines other than the expected kinds are counted in
#                runs/<scenario>.other_lines); 30 if a scenario did not exit 0. THE TRY COUNTS FROM THE START OF THIS
#                STEP. The scenario files stay uncommitted in experiments/t1_scope/results/t72_late_mix/scenarios/
#   read ROWC    once, after runs exited 0 at ROWC: the bars, after every void condition, into
#                experiments/t1_scope/results/t72_late_mix/ (t72.json, t72.txt; the outcome words to read.out); the
#                main session commits the scenario files and the reading together and quotes the outcome line,
#                "Outcome: T72 PASS", "Outcome: T72 STOP (...)" or "Outcome: T72 VOID (...)", in row T72's outcome
#                cell, which experiments/builds/t1_lm_m20_line.sh reads by code
# Refusals write their code to ~/v/scratch/t72/<step>.refused and never touch <step>.exit or <step>.head (a call
# without its arguments writes usage.refused, 2): 97 an unknown step; 98 an unclean tree (any change outside the row's
# results folder); 96 TRIES.md does not hold exactly one row T72 whose first cell starts 'T72 (' and names
# experiments/t1_scope/t72_late_mix.py and experiments/t1_scope/t72_late_mix.sh, each as a whole path, and exactly
# one `shares.json sha256 <hex>`, the one row N64 names (row N64 naming exactly one); 95 ROWC is not a full hash or HEAD
# is not ROWC; 94 the step before has not exited 0 at ROWC; 93 an input is not in place (the E10.5 file, the shares
# file tracked with the named sha256, T30's committed scenario file); 91 the step's output exists (nothing is ever
# overwritten); 88 the checkout's .venv does not name its own src. A step that runs writes "running" to <step>.exit and
# HEAD to <step>.head first, then its exit code to <step>.exit; it runs with NUMPY_MADVISE_HUGEPAGE=0, TMPDIR under
# ~/v/scratch/t72 and $HOME/.local/bin first on PATH, and writes its stdout to <step>.out and its stderr to
# <step>.err. Nothing here deletes anything.
set -u
export PATH=$HOME/.local/bin:$PATH
S="$HOME/v/scratch/t72"
mkdir -p "$S/tmp"
if [ "$#" -lt 2 ] || [ -z "$1" ] || [ -z "$2" ]; then
  echo "usage: t72_late_mix.sh tests|runs|read ROWC" >&2
  echo 2 > "$S/usage.refused"
  exit 2
fi
NAME="$1"
ROWC="$2"
refuse() { echo "$1" > "$S/$NAME.refused"; exit "$1"; }
start() { echo running > "$S/$NAME.exit"; echo "$HEAD" > "$S/$NAME.head"; }
finish() { echo "$1" > "$S/$NAME.exit"; exit "$1"; }
case "$NAME" in
  tests) PREV="" ;;
  runs) PREV="tests" ;;
  read) PREV="runs" ;;
  *)
    echo "unknown step $NAME" >&2
    NAME=unknown
    refuse 97 ;;
esac
TOP=$(git rev-parse --show-toplevel 2>/dev/null) || refuse 98
cd "$TOP" || refuse 98
export NUMPY_MADVISE_HUGEPAGE=0 TMPDIR="$S/tmp"
export UV_PROJECT_ENVIRONMENT="$PWD/.venv"
PY="uv run --frozen --no-sync python -I -B"
TRIES=experiments/harness/TRIES.md
RES=experiments/t1_scope/results/t72_late_mix
SHARES=experiments/t1_scope/results/n64_late_mix/shares.json
T30=experiments/t1_scope/results/t30_pricing/scenarios/t30_d1.0_cmown.json
SCRIPT=experiments/t1_scope/t72_late_mix.py
ROW_PATHS=('experiments/t1_scope/t72_late_mix\.py' 'experiments/t1_scope/t72_late_mix\.sh')  # escaped for grep -E
SCEN=(own_own own_lm16 0.30_own 0.30_lm16 0.33_own 0.33_lm16 lm16_own lm16_lm16 0.45_own 0.45_lm16 0.60_own 0.60_lm16)
OK_LINE='^(own|0\.30|0\.33|lm16|0\.45|0\.60) (own|lm16) seed [0-9]+: [0-9]+ s$|^wrote t72_cm[a-z0-9.]+_blood[a-z0-9]+\.json$'
names() { printf '%s\n' "$1" | grep -qE "(^|[^[:alnum:]_./-])$2([^[:alnum:]_./-]|\$)"; }
if git status --porcelain --untracked-files=all | grep -v "^?? $RES/" | grep -q .; then
  echo "refused: the tree is not clean" >&2
  refuse 98
fi
HEAD=$(git rev-parse HEAD)
[ "$(grep -c '^| T72 ' "$TRIES" 2>/dev/null)" = 1 ] || { echo "refused: TRIES.md does not hold exactly one row T72" >&2; refuse 96; }
first=$(grep '^| T72 ' "$TRIES" | cut -d'|' -f2)
printf '%s\n' "$first" | grep -q '^ T72 (' || { echo "refused: row T72's first cell does not start 'T72 ('" >&2; refuse 96; }
for p in "${ROW_PATHS[@]}"; do
  names "$first" "$p" || { echo "refused: row T72's first cell does not name $p as a whole path" >&2; refuse 96; }
done
named=$(printf '%s\n' "$first" | grep -oE '`shares\.json sha256 [0-9a-f]{64}`' | sort -u)
n64=$(grep '^| N64 ' "$TRIES" | grep -oE '`shares\.json sha256 [0-9a-f]{64}`' | sort -u)
if [ "$(grep -c '^| N64 ' "$TRIES")" != 1 ] || [ "$(printf '%s\n' "$named" | grep -c .)" != 1 ] \
    || [ "$(printf '%s\n' "$n64" | grep -c .)" != 1 ] || [ "$named" != "$n64" ]; then
  echo "refused: rows N64 and T72 do not each name the one same \`shares.json sha256 <hex>\`" >&2
  refuse 96
fi
case "$ROWC" in
  *[!0-9a-f]*) refuse 95 ;;
esac
[ "${#ROWC}" -eq 40 ] || refuse 95
[ "$HEAD" = "$ROWC" ] || { echo "refused: HEAD is not $ROWC" >&2; refuse 95; }
for P in $PREV; do
  if [ "$(cat "$S/$P.exit" 2>/dev/null)" != 0 ] || [ "$(cat "$S/$P.head" 2>/dev/null)" != "$HEAD" ]; then
    echo "refused: step $P has not exited 0 at $HEAD" >&2
    refuse 94
  fi
done
venv=0
for pth in .venv/lib/python3*/site-packages/vec.pth; do
  if [ -f "$pth" ] && [ "$(readlink -f "$(< "$pth")")" = "$(readlink -f src)" ]; then venv=1; fi
done
[ "$venv" = 1 ] || { echo "refused: this checkout has no .venv of its own whose vec.pth names its src" >&2; refuse 88; }
inputs() {  # the E10.5 file, the shares file tracked with the named sha256, T30's committed scenario file
  [ -f data/raw/T1/E10.5_RNA.h5ad ] || { echo "refused: data/raw/T1/E10.5_RNA.h5ad is not in place" >&2; refuse 93; }
  git ls-files --error-unmatch "$SHARES" > /dev/null 2>&1 || { echo "refused: $SHARES is not tracked" >&2; refuse 93; }
  [ "$named" = "\`shares.json sha256 $(sha256sum "$SHARES" | cut -d' ' -f1)\`" ] \
    || { echo "refused: $SHARES is not the file rows N64 and T72 name" >&2; refuse 93; }
  git ls-files --error-unmatch "$T30" > /dev/null 2>&1 || { echo "refused: $T30 is not tracked" >&2; refuse 93; }
}

case "$NAME" in
  tests)
    start
    export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
    $PY -m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp="$S/tmp/pytest" \
      tests/test_t1_late_mix.py tests/test_t1_lm_m20_launchers.py tests/test_t1_scope.py \
      > "$S/$NAME.out" 2> "$S/$NAME.err"
    finish $? ;;
  runs)
    inputs
    [ -e "$RES/scenarios" ] && { echo "refused: $RES/scenarios exists" >&2; refuse 91; }
    { [ -e "$RES/t72.json" ] || [ -e "$RES/t72.txt" ]; } && { echo "refused: a reading exists" >&2; refuse 91; }
    [ -e "$S/runs" ] && { echo "refused: $S/runs exists" >&2; refuse 91; }
    start
    mkdir -p "$S/runs"
    export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
    pids=()
    for sc in "${SCEN[@]}"; do
      (
        $PY "$SCRIPT" run --cm "${sc%_*}" --blood "${sc#*_}" --commit "$HEAD" --jobs 2 > "$S/runs/$sc.log" 2>&1
        echo $? > "$S/runs/$sc.exit"
      ) &
      pids+=($!)
    done
    for p in "${pids[@]}"; do wait "$p" || true; done
    bad=0
    for sc in "${SCEN[@]}"; do
      [ "$(cat "$S/runs/$sc.exit" 2>/dev/null)" = 0 ] || bad=1
      n=$(grep -cvE "$OK_LINE" "$S/runs/$sc.log" 2>/dev/null) || true
      echo "${n:-unreadable}" > "$S/runs/$sc.other_lines"
    done
    [ "$bad" = 0 ] || finish 30
    finish 0 ;;
  read)
    inputs
    { [ -e "$RES/t72.json" ] || [ -e "$RES/t72.txt" ]; } && { echo "refused: a reading exists" >&2; refuse 91; }
    start
    export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
    $PY "$SCRIPT" read --commit "$HEAD" > "$S/$NAME.out" 2> "$S/$NAME.err"
    finish $? ;;
esac
