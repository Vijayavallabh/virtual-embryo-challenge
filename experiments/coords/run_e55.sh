#!/usr/bin/env bash
# Row E55 (experiments/harness/TRIES.md) on the DGX: the practice read of OSR (fgw_shape0x_ov_sr) and its two placebos,
# all steps in one call, from the root of a worktree of the campaign checkout's repository (HEAD = COMMIT, clean):
#   bash experiments/coords/run_e55.sh COMMIT RUN
# COMMIT is the run commit (full hash), the row's commit. RUN names the scratch folder ~/v/scratch/e55/RUN: letters,
# digits, '.', '_' and '-' only, not starting with '.' (else exit 90); refused if it exists (exit 91). A call without
# both arguments exits 2. Environment: E55_GATE_RECORD, the folder of the target-free gate's run at d145d78 (its ten
# pair files and gate.json; default experiments/coords/results/osr_gate_g1), E55_GATE_SHA256, the sha256 of its
# gate.json as the row quotes it (default 7c0fc6b6...7b4d67, the read's GATE_SHA256), E55_DATA (default
# ~/v/vec-campaign/data/raw), E55_VENV (default ~/v/vec-campaign/.venv). A record without gate.json exits 93 (preflight).
# Steps, in order; each writes RUN/<step>.exit, first the word running, then its exit code, and the call stops at the
# first step that does not exit 0; RUN/done.exit holds running, then the last step's code (0 when all exited 0):
#   preflight  the tracked tree clean (else 98), HEAD = COMMIT (else 99)
#   tests      the unit tests of the form, the builder option, the gate and the read (synthetic arrays only)
#   gate0      osr_gate.py pair for the ten pairs at COMMIT into RUN/runs/gate0, one process per pair, then
#              osr_gate.py read there (step 0: the read compares it with the d145d78 record, ruling 4)
#   gate       the floor, C0 and O under main against E47 and E49, one process per task (exit 3: some task fails)
#   run        O, the arm and the placebos under both builds, one process per task (refused without the gate record,
#              its recompute and the nine gate files; the script refuses without row E55)
#   read       the reading into experiments/coords/results/e55_overlay_stretch/ (never overwritten)
# BLAS threads 3 a process: at most ten at once, 30 threads. NUMPY_MADVISE_HUGEPAGE=0; TMPDIR under RUN. Every
# process's output goes to /dev/null; per-process exit codes are RUN/rc/<step>-<name>. Nothing is deleted or uploaded.
set -u
if [ "$#" -ne 2 ] || [ -z "$1" ] || [ -z "$2" ]; then
  echo "usage: run_e55.sh COMMIT RUN" >&2
  exit 2
fi
COMMIT="$1"
RUN="$2"
if ! [[ "$RUN" =~ ^[A-Za-z0-9_-][A-Za-z0-9._-]*$ ]]; then
  echo "refused: RUN may hold letters, digits, '.', '_' and '-' only, not starting with '.'" >&2
  exit 90
fi
S="$HOME/v/scratch/e55/$RUN"
if [ -e "$S" ]; then
  echo "refused: $S exists" >&2
  exit 91
fi
mkdir -p "$S/tmp" "$S/rc" "$S/runs" || exit 92
echo running > "$S/done.exit"
STEP=preflight
finish() { echo "$1" > "$S/$STEP.exit"; echo "$1" > "$S/done.exit"; exit "$1"; }
begin() { STEP="$1"; echo running > "$S/$STEP.exit"; }
begin preflight
GATE_RECORD="${E55_GATE_RECORD:-experiments/coords/results/osr_gate_g1}"
GATE_SHA="${E55_GATE_SHA256:-7c0fc6b6398092d54f6e4318bd8c79c1e092c4ec7e1edbbb6b86fcbede7b4d67}"
if [ ! -f "$GATE_RECORD/gate.json" ] || [ -z "$GATE_SHA" ]; then
  echo "refused: no gate record at $GATE_RECORD (its gate.json)" >&2
  finish 93
fi
export OMP_NUM_THREADS=3 OPENBLAS_NUM_THREADS=3 MKL_NUM_THREADS=3 NUMEXPR_NUM_THREADS=3 NUMPY_MADVISE_HUGEPAGE=0
export TMPDIR="$S/tmp"
export PATH="$HOME/.local/bin:$PATH"
export UV_PROJECT_ENVIRONMENT="${E55_VENV:-$HOME/v/vec-campaign/.venv}"
DATA="${E55_DATA:-$HOME/v/vec-campaign/data/raw}"
PY="uv run --frozen --no-sync python -I -B"
ST=$(git status --porcelain --untracked-files=no 2>/dev/null) || finish 98
if [ -n "$ST" ]; then
  echo "refused: the tracked tree is not clean" >&2
  finish 98
fi
HEAD=$(git rev-parse HEAD 2>/dev/null) || finish 99
if [ "$HEAD" != "$COMMIT" ]; then
  echo "refused: HEAD is not $COMMIT" >&2
  finish 99
fi
echo 0 > "$S/preflight.exit"
SCRIPT=experiments/coords/e55_overlay_stretch_read.py
GATE=experiments/coords/osr_gate.py
R="$S/runs"
OUT=experiments/coords/results/e55_overlay_stretch
TASKS=(half late new1 new2 e75 e725 x1 x2 x3)
PAIRS=(test half late new1 new2 e75 e725 x1 x2 x3)

parallel() {  # $1 = step: one process per task, all nine at once; 0 all exit 0, 3 a gate failure, 1 otherwise
  local pids=() t bad=0 fail=0 e
  for t in "${TASKS[@]}"; do
    (
      if [ "$1" = run ]; then
        $PY "$SCRIPT" run --task "$t" --commit "$COMMIT" --runs "$R" --data "$DATA" \
          --gate-record "$GATE_RECORD" --gate-sha256 "$GATE_SHA" > /dev/null 2>&1
      else
        $PY "$SCRIPT" "$1" --task "$t" --commit "$COMMIT" --runs "$R" --data "$DATA" > /dev/null 2>&1
      fi
      echo $? > "$S/rc/$1-$t"
    ) &
    pids+=($!)
  done
  for p in "${pids[@]}"; do wait "$p"; done
  for t in "${TASKS[@]}"; do
    e=$(cat "$S/rc/$1-$t" 2>/dev/null)
    if [ "$e" = 3 ]; then fail=1; elif [ "$e" != 0 ]; then bad=1; fi
  done
  if [ $bad = 1 ]; then return 1; fi
  if [ $fail = 1 ]; then return 3; fi
  return 0
}

begin tests
$PY -m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp="$S/tmp/pytest" \
  tests/test_e55_overlay_stretch.py tests/test_osr_gate.py tests/test_shape_interp.py \
  tests/test_shape_interp_build.py tests/test_e49_shape_interp_read.py tests/test_t2_build.py > /dev/null 2>&1
e=$?
[ "$e" = 0 ] || finish "$e"
echo 0 > "$S/tests.exit"

begin gate0
mkdir -p "$R/gate0"
pids=()
for x in "${PAIRS[@]}"; do
  (
    $PY "$GATE" pair --name "$x" --data "$DATA" --out "$R/gate0" --commit "$COMMIT" > /dev/null 2>&1
    echo $? > "$S/rc/gate0-$x"
  ) &
  pids+=($!)
done
for p in "${pids[@]}"; do wait "$p"; done
bad=0
for x in "${PAIRS[@]}"; do
  if [ "$(cat "$S/rc/gate0-$x" 2>/dev/null)" != 0 ]; then bad=1; fi
done
$PY "$GATE" read --out "$R/gate0" > /dev/null 2>&1 || finish 5
[ "$bad" = 0 ] || finish 4
echo 0 > "$S/gate0.exit"

begin gate
parallel gate
e=$?
[ "$e" = 0 ] || finish "$e"
echo 0 > "$S/gate.exit"

begin run
parallel run
e=$?
[ "$e" = 0 ] || finish "$e"
echo 0 > "$S/run.exit"

begin read
if [ -e "$OUT/reading.json" ] || [ -e "$OUT/reading.txt" ]; then
  echo "refused: a reading exists under $OUT and is never overwritten" >&2
  finish 95
fi
$PY "$SCRIPT" read --commit "$COMMIT" --runs "$R" --out "$OUT" --gate-record "$GATE_RECORD" \
  --gate-sha256 "$GATE_SHA" > /dev/null 2>&1
finish $?
