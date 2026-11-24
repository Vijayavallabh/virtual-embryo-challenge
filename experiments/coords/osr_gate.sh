#!/usr/bin/env bash
# The OSR gate (experiments/coords/osr_gate.py) on the DGX: target-free, scoreless. Run from the checkout's root:
#   bash experiments/coords/osr_gate.sh COMMIT RUN
# COMMIT is the full sha the checkout's HEAD must be (each pair command refuses otherwise, exit 2, and refuses an
# unclean tracked tree). RUN names the scratch folder ~/v/scratch/osr_gate/RUN: letters, digits, '.', '_' and '-' only,
# not starting with '.' (else exit 90); refused if it exists (exit 91). A call without both arguments exits 2.
# It starts the ten pair commands (test and the nine tasks) as ten background processes, 3 BLAS threads each, waits for
# them, runs read, and writes RUN/done.exit: "running", then 0 when every pair's process and exit file say 0 and read
# exited 0; 4 when a pair did not exit 0 (read still runs: stop clause S5); 5 when read did not exit 0. Data:
# OSR_DATA, default data/raw under the checkout. Each process's output goes to /dev/null; the results are RUN/<name>.json,
# RUN/<name>.exit, RUN/<name>.rc (the process status), RUN/gate.json and RUN/gate.txt. Nothing is deleted or uploaded.
set -u
if [ "$#" -ne 2 ] || [ -z "$1" ] || [ -z "$2" ]; then
  echo "usage: osr_gate.sh COMMIT RUN" >&2
  exit 2
fi
COMMIT="$1"
RUN="$2"
if ! [[ "$RUN" =~ ^[A-Za-z0-9_-][A-Za-z0-9._-]*$ ]]; then
  echo "refused: RUN may hold letters, digits, '.', '_' and '-' only, not starting with '.'" >&2
  exit 90
fi
S="$HOME/v/scratch/osr_gate/$RUN"
if [ -e "$S" ]; then
  echo "refused: $S exists" >&2
  exit 91
fi
mkdir -p "$S/tmp" || exit 92
echo running > "$S/done.exit"
export OMP_NUM_THREADS=3 OPENBLAS_NUM_THREADS=3 MKL_NUM_THREADS=3 NUMPY_MADVISE_HUGEPAGE=0 TMPDIR="$S/tmp"
export PATH="$HOME/.local/bin:$PATH"
DATA="${OSR_DATA:-data/raw}"
PY="uv run --frozen --no-sync python -I -B"
SCRIPT=experiments/coords/osr_gate.py
NAMES=(test half late new1 new2 e75 e725 x1 x2 x3)
pids=()
for x in "${NAMES[@]}"; do
  (
    $PY "$SCRIPT" pair --name "$x" --data "$DATA" --out "$S" --commit "$COMMIT" > /dev/null 2>&1
    echo $? > "$S/$x.rc"
  ) &
  pids+=($!)
done
for p in "${pids[@]}"; do wait "$p"; done
bad=0
for x in "${NAMES[@]}"; do
  if [ "$(cat "$S/$x.rc" 2>/dev/null)" != 0 ] || [ "$(cat "$S/$x.exit" 2>/dev/null)" != 0 ]; then bad=1; fi
done
$PY "$SCRIPT" read --out "$S" > /dev/null 2>&1
r=$?
if [ "$r" != 0 ]; then
  echo 5 > "$S/done.exit"
  exit 5
fi
if [ "$bad" != 0 ]; then
  echo 4 > "$S/done.exit"
  exit 4
fi
echo 0 > "$S/done.exit"
exit 0
