#!/usr/bin/env bash
# Row E58 (experiments/harness/TRIES.md) on the DGX: gate G1 of row R59's file T2embryo_s2x_ovt23.h5ad, the counted
# practice read of the arm fgw_shape0x_ov_t23 (O with its mix leaned a sixth toward the later bracket) against O, all
# steps in one call, from the root of a clean checkout or worktree of the campaign repository at the run commit:
#   bash experiments/coords/run_e58.sh COMMIT RUN
# Run it detached and wait on RUN/done.exit with a background poller; never a foreground sleep. COMMIT is the run commit
# (full hash), the commit of row E58. RUN names the scratch folder ~/v/scratch/e58/RUN: letters, digits, '.', '_' and
# '-' only, not starting with '.' (else exit 90); refused if it exists (91). A call without both arguments exits 2.
# Environment: E58_DATA (default ~/v/vec-campaign/data/raw), E58_VENV (default ~/v/vec-campaign/.venv).
# Steps, in order; each writes RUN/<step>.exit, first the word running, then its exit code; RUN/done.exit holds running,
# then the call's code (0 when all exited 0). The call stops at the first step that does not exit 0; when that step is
# gate or run, the read step is still called once (so the VOID reading, with its reasons, is written) and the call exits
# with the failing step's code:
#   preflight  HEAD = COMMIT (else 96), the tracked tree clean (else 98), and TRIES.md holds exactly one well-formed row
#              E58 (experiments/coords/e58_row.py row: eleven cells, its first cell opening with the gate's fixed words,
#              which say it moves the overlay's form; else 95); nothing else runs before these three checks
#   tests      the unit tests of the read, the row check, the builder option and the overlay (synthetic arrays only)
#   gate       the floor and O under main against E47 and E55 (the J7 guard), one process per task (3: a task fails)
#   run        the floor, O and the arm under main, one process per task (refused without the seven gate files)
#   read       the reading into experiments/coords/results/e58_mix_lean/ (never overwritten: 95)
# The seven tasks: half, late and new1 count (two truths); new2, x1, x2 and x3 are reported; e75 and e725 are never run
# (barred for any change to the mixing fraction). BLAS threads: 4 for the tests, 8 a process for gate, run and read
# (seven at once), as run_e55.sh's run, since the 0.01 guard compares with E55's run made at 8. NUMPY_MADVISE_HUGEPAGE=0;
# TMPDIR under RUN. Every process's output goes to /dev/null; per-process exit codes are RUN/rc/<step>-<task>. Nothing
# is deleted or uploaded.
set -u
if [ "$#" -ne 2 ] || [ -z "$1" ] || [ -z "$2" ]; then
  echo "usage: run_e58.sh COMMIT RUN" >&2
  exit 2
fi
COMMIT="$1"
RUN="$2"
if ! [[ "$RUN" =~ ^[A-Za-z0-9_-][A-Za-z0-9._-]*$ ]]; then
  echo "refused: RUN may hold letters, digits, '.', '_' and '-' only, not starting with '.'" >&2
  exit 90
fi
export PATH=$HOME/.local/bin:$PATH
S="$HOME/v/scratch/e58/$RUN"
if [ -e "$S" ]; then
  echo "refused: $S exists (a run folder is never reused)" >&2
  exit 91
fi
mkdir -p "$S/tmp" "$S/rc" "$S/runs" || exit 92
echo running > "$S/done.exit"
STEP=preflight
finish() { echo "$1" > "$S/$STEP.exit"; echo "$1" > "$S/done.exit"; exit "$1"; }
begin() { STEP="$1"; echo running > "$S/$STEP.exit"; }
begin preflight
HEAD=$(git rev-parse HEAD 2>/dev/null) || finish 96
if [ "$HEAD" != "$COMMIT" ]; then
  echo "refused: HEAD is not $COMMIT" >&2
  finish 96
fi
ST=$(git status --porcelain --untracked-files=no 2>/dev/null) || finish 98
if [ -n "$ST" ]; then
  echo "refused: the tracked tree is not clean" >&2
  finish 98
fi
if ! python3 -I -B experiments/coords/e58_row.py row experiments/harness/TRIES.md 2>/dev/null; then
  echo "refused: TRIES.md does not hold exactly one well-formed row E58 (experiments/coords/e58_row.py)" >&2
  finish 95
fi
echo 0 > "$S/preflight.exit"
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 NUMEXPR_NUM_THREADS=4 NUMPY_MADVISE_HUGEPAGE=0
export TMPDIR="$S/tmp"
export UV_PROJECT_ENVIRONMENT="${E58_VENV:-$HOME/v/vec-campaign/.venv}"
DATA="${E58_DATA:-$HOME/v/vec-campaign/data/raw}"
PY="uv run --frozen --no-sync python -I -B"
SCRIPT=experiments/coords/e58_mix_lean_read.py
R="$S/runs"
OUT=experiments/coords/results/e58_mix_lean
TASKS=(half late new1 new2 x1 x2 x3)

parallel() {  # $1 = step: one process per task, all seven at once; 0 all exit 0, 3 a gate failure, 1 otherwise
  local pids=() t bad=0 fail=0 e
  for t in "${TASKS[@]}"; do
    (
      $PY "$SCRIPT" "$1" --task "$t" --commit "$COMMIT" --runs "$R" --data "$DATA" > /dev/null 2>&1
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

call_read() {  # the read step: writes reading.json and reading.txt (never overwritten); its code
  begin read
  if [ -e "$OUT/reading.json" ] || [ -e "$OUT/reading.txt" ]; then
    echo "refused: a reading exists under $OUT and is never overwritten" >&2
    echo 95 > "$S/read.exit"
    return 95
  fi
  $PY "$SCRIPT" read --commit "$COMMIT" --runs "$R" --out "$OUT" > /dev/null 2>&1
  local e=$?
  echo "$e" > "$S/read.exit"
  return "$e"
}
fail_then_read() {  # $1 = the failing step's code: record it, call the read once, exit with that code
  local code="$1"
  echo "$code" > "$S/$STEP.exit"
  call_read
  echo "$code" > "$S/done.exit"
  exit "$code"
}

begin tests
$PY -m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp="$S/tmp/pytest" \
  tests/test_e58_mix_lean_read.py tests/test_mix_lean_build.py tests/test_shape_interp_build.py \
  tests/test_t2_build.py > /dev/null 2>&1
e=$?
[ "$e" = 0 ] || finish "$e"
echo 0 > "$S/tests.exit"

export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 NUMEXPR_NUM_THREADS=8
begin gate
parallel gate
e=$?
[ "$e" = 0 ] || fail_then_read "$e"
echo 0 > "$S/gate.exit"

begin run
parallel run
e=$?
[ "$e" = 0 ] || fail_then_read "$e"
echo 0 > "$S/run.exit"

call_read
finish $?
