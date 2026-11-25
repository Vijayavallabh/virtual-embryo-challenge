#!/usr/bin/env bash
# Row N64 (experiments/harness/TRIES.md): step 0 of row R54's late-mix line, experiments/t1_scope/late_mix.py
# (target-free, not a try), one step per call, run from the root of the checkout or detached worktree it is called in
# (it moves to that checkout's top level first), with that folder's own environment and data (the E10.5 file
# data/raw/T1/E10.5_RNA.h5ad). Wait on the exit file with a background poller; never a foreground sleep. Every step
# takes NROWC, the full 40-character hash of the commit that holds row N64, the X002 amendment of 2026-11-25 and this
# code, and runs only at HEAD = NROWC. Steps, in order; the second is refused (exit 94) unless the first wrote exit 0
# for the same commit:
#   tests NROWC   the targeted tests (synthetic arrays only)
#   step0 NROWC   late_mix.py --commit NROWC: writes shares.json, n64.json and n64.txt into
#                 experiments/t1_scope/results/n64_late_mix/ and prints the shares file's sha256 to step0.out. The files
#                 stay uncommitted there; the main session commits them and writes `shares.json sha256 <hex>` into row
#                 N64 (the builder's --late-mix and rows T72's read take that string)
# Refusals write their code to ~/v/scratch/n64/<step>.refused and never touch <step>.exit or <step>.head (a call
# without its arguments writes usage.refused, 2): 97 an unknown step; 98 an unclean tree (any change outside the step's
# results folder); 96 TRIES.md does not hold exactly one row N64 whose first cell starts 'N64 (' and names
# experiments/t1_scope/late_mix.py and experiments/t1_scope/n64_step0.sh, each as a whole path; 95 NROWC is not a full
# hash or HEAD is not NROWC; 94 the tests have not exited 0 at NROWC; 93 the E10.5 file or a committed manifest is not
# in place; 92 EXTERNAL_DATA.md's row X002 does not hold the amendment of 2026-11-25 ("Amendment of 2026-11-25 (rows
# N64, T72 and R54)"); 91 an output exists (nothing is ever overwritten); 88 the checkout's .venv does not name its own
# src. A step that runs writes "running" to <step>.exit and HEAD to <step>.head first, then its exit code to
# <step>.exit; it runs at 8 BLAS threads with NUMPY_MADVISE_HUGEPAGE=0, TMPDIR under ~/v/scratch/n64 and
# $HOME/.local/bin first on PATH, and writes its stdout to <step>.out and its stderr to <step>.err. The Python step
# exits 0 written, 1 a refusal or a recomputation that differs (nothing written). Nothing here deletes anything.
set -u
export PATH=$HOME/.local/bin:$PATH
S="$HOME/v/scratch/n64"
mkdir -p "$S/tmp"
if [ "$#" -lt 2 ] || [ -z "$1" ] || [ -z "$2" ]; then
  echo "usage: n64_step0.sh tests|step0 NROWC" >&2
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
  step0) PREV="tests" ;;
  *)
    echo "unknown step $NAME" >&2
    NAME=unknown
    refuse 97 ;;
esac
TOP=$(git rev-parse --show-toplevel 2>/dev/null) || refuse 98
cd "$TOP" || refuse 98
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 NUMPY_MADVISE_HUGEPAGE=0 TMPDIR="$S/tmp"
export UV_PROJECT_ENVIRONMENT="$PWD/.venv"
PY="uv run --frozen --no-sync python -I -B"
TRIES=experiments/harness/TRIES.md
RES=experiments/t1_scope/results/n64_late_mix
ROW_PATHS=('experiments/t1_scope/late_mix\.py' 'experiments/t1_scope/n64_step0\.sh')  # escaped for grep -E
names() { printf '%s\n' "$1" | grep -qE "(^|[^[:alnum:]_./-])$2([^[:alnum:]_./-]|\$)"; }
if git status --porcelain --untracked-files=all | grep -v "^?? $RES/" | grep -q .; then
  echo "refused: the tree is not clean" >&2
  refuse 98
fi
HEAD=$(git rev-parse HEAD)
[ "$(grep -c '^| N64 ' "$TRIES" 2>/dev/null)" = 1 ] || { echo "refused: TRIES.md does not hold exactly one row N64" >&2; refuse 96; }
first=$(grep '^| N64 ' "$TRIES" | cut -d'|' -f2)
printf '%s\n' "$first" | grep -q '^ N64 (' || { echo "refused: row N64's first cell does not start 'N64 ('" >&2; refuse 96; }
for p in "${ROW_PATHS[@]}"; do
  names "$first" "$p" || { echo "refused: row N64's first cell does not name $p as a whole path" >&2; refuse 96; }
done
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

case "$NAME" in
  tests)
    start
    $PY -m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp="$S/tmp/pytest" \
      tests/test_t1_late_mix.py tests/test_t1_lm_m20_launchers.py tests/test_t1_scope.py \
      > "$S/$NAME.out" 2> "$S/$NAME.err"
    finish $? ;;
  step0)
    [ -f data/raw/T1/E10.5_RNA.h5ad ] || { echo "refused: data/raw/T1/E10.5_RNA.h5ad is not in place" >&2; refuse 93; }
    for f in experiments/external/results/feng_relabel_late/MANIFEST.json \
             experiments/external/results/feng_relabel/MANIFEST.json; do
      git ls-files --error-unmatch "$f" > /dev/null 2>&1 || { echo "refused: $f is not tracked" >&2; refuse 93; }
    done
    grep '^| X002 |' EXTERNAL_DATA.md | grep -qF 'Amendment of 2026-11-25 (rows N64, T72 and R54)' \
      || { echo "refused: EXTERNAL_DATA.md's row X002 does not hold the amendment of 2026-11-25" >&2; refuse 92; }
    for f in shares.json n64.json n64.txt; do
      [ -e "$RES/$f" ] && { echo "refused: $RES/$f exists" >&2; refuse 91; }
    done
    start
    $PY experiments/t1_scope/late_mix.py --commit "$HEAD" > "$S/$NAME.out" 2> "$S/$NAME.err"
    finish $? ;;
esac
