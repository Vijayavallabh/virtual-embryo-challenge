#!/usr/bin/env bash
# Row E59 (experiments/harness/TRIES.md) on the DGX: the gate of row R71's files T2embryo_s2x_ovbf.h5ad (BF_O) and
# T2embryo_s2xs_ovbf.h5ad (BF_SO), all steps in one call, from the root of a clean checkout or worktree of the campaign
# repository at the run commit:
#   bash experiments/coords/run_bf.sh COMMIT RUN
# Run it detached and wait on RUN/done.exit with a background poller; never a foreground sleep. COMMIT is the run commit
# (full hash), the commit of row E59. RUN names the scratch folder ~/v/scratch/bf_gate/RUN: letters, digits, '.', '_'
# and '-' only, not starting with '.' (else exit 90); refused if it exists (91). A call without both arguments exits 2.
# Environment: BF_DATA (default ~/v/vec-campaign/data/raw), BF_VENV (default ~/v/vec-campaign/.venv; the read script
# puts this checkout's src/ first and refuses code from elsewhere).
# Steps, in order; each writes RUN/<step>.exit, first the word running, then its exit code; RUN/done.exit holds running,
# then the call's code (0 when all exited 0). The call stops at the first step that does not exit 0; when that step is
# g0read (4: G0 said STOP, nothing is counted) or run, the read step is still called once (so the G0 STOP or VOID
# reading, with its reasons, is written) and the call exits with the failing step's code:
#   preflight  the clock before 05:00 UTC on 2026-11-26 (else 84: the gate's budget is about 3 hours, 2.5 to 3.5, and the
#              challenge's fix 9 builds no BF file after 08:00 UTC, so a later start could spend the four tries on files
#              that cannot be built), HEAD = COMMIT (else 96), the tracked tree clean (else 98), TRIES.md holds
#              exactly one well-formed row E59 (experiments/coords/bf_row.py row: eleven cells, its first cell opening
#              with the gate's fixed words; else 95) and exactly one row R71 that is the captain's order naming BF on
#              O and on SO only and `--captain-ruling R71` (else 94); nothing else runs before these checks
#   tests      the unit tests of the form, the read, the row check, the builder option and the launchers (synthetic
#              arrays only)
#   g0         G0 at the ten build seeds, one process each (target-free: the test pair's inputs, never a target)
#   g0read     G0's verdict: 0 PASS, 4 STOP (the run never starts and nothing is counted)
#   run        the nine tasks, one process each (the four tries count from this step's start)
#   read       the reading into experiments/coords/results/bf_gate/ (never overwritten: 95); --run-started when the run
#              step began, so the reading counts the four tries from the run step's start (rule 4)
# BLAS threads: 4 for the tests, 8 a process for g0, run and read (nine at once), as run_e55.sh and run_e58.sh, since
# the 0.01 guard compares with E55's run made at 8. NUMPY_MADVISE_HUGEPAGE=0; TMPDIR under RUN. Every process's output
# goes to /dev/null; per-process exit codes are RUN/rc/<step>-<task or seed>. Nothing is deleted or uploaded.
set -u
if [ "$#" -ne 2 ] || [ -z "$1" ] || [ -z "$2" ]; then
  echo "usage: run_bf.sh COMMIT RUN" >&2
  exit 2
fi
COMMIT="$1"
RUN="$2"
if ! [[ "$RUN" =~ ^[A-Za-z0-9_-][A-Za-z0-9._-]*$ ]]; then
  echo "refused: RUN may hold letters, digits, '.', '_' and '-' only, not starting with '.'" >&2
  exit 90
fi
export PATH=$HOME/.local/bin:$PATH
S="$HOME/v/scratch/bf_gate/$RUN"
if [ -e "$S" ]; then
  echo "refused: $S exists (a run folder is never reused)" >&2
  exit 91
fi
mkdir -p "$S/tmp" "$S/rc" "$S/runs/g0" || exit 92
echo running > "$S/done.exit"
STEP=preflight
finish() { echo "$1" > "$S/$STEP.exit"; echo "$1" > "$S/done.exit"; exit "$1"; }
begin() { STEP="$1"; echo running > "$S/$STEP.exit"; }
begin preflight
LATEST_START=1795669200  # 2026-11-26T05:00:00Z: 08:00 UTC (fix 9's last build) minus the gate's 3-hour budget
if [ "$(date -u +%s)" -ge "$LATEST_START" ]; then
  echo "refused: after 05:00 UTC on 2026-11-26 the gate's files could not be built by 08:00 UTC (fix 9)" >&2
  finish 84
fi
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
if ! python3 -I -B experiments/coords/bf_row.py row experiments/harness/TRIES.md 2>/dev/null; then
  echo "refused: TRIES.md does not hold exactly one well-formed row E59 (experiments/coords/bf_row.py)" >&2
  finish 95
fi
if ! python3 -I -B -c 'import sys
rows = [ln for ln in open(sys.argv[1], encoding="utf-8").read().splitlines() if ln.startswith("| R71 ")]
first = rows[0][2:].split("|", 1)[0] if len(rows) == 1 else ""
ok = (first.startswith("R71 (the captain\x27s order") and "ruling, no run" in first
      and "BF on O and on SO only" in first and "`--captain-ruling R71`" in first)
sys.exit(0 if ok else 1)' experiments/harness/TRIES.md 2>/dev/null; then
  echo "refused: TRIES.md does not hold exactly one row R71, the captain's order naming the files and the ruling" >&2
  finish 94
fi
echo 0 > "$S/preflight.exit"
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 NUMEXPR_NUM_THREADS=4 NUMPY_MADVISE_HUGEPAGE=0
export TMPDIR="$S/tmp"
export UV_PROJECT_ENVIRONMENT="${BF_VENV:-$HOME/v/vec-campaign/.venv}"
DATA="${BF_DATA:-$HOME/v/vec-campaign/data/raw}"
PY="uv run --frozen --no-sync python -I -B"
SCRIPT=experiments/coords/bf_read.py
R="$S/runs"
OUT=experiments/coords/results/bf_gate
TASKS=(half late new1 new2 e75 e725 x1 x2 x3)
SEEDS=(1000001 1000002 1000003 1000004 1000005 1000006 1000007 1000008 1000009 1000010)

call_read() {  # the read step: writes reading.json and reading.txt (never overwritten); its code
  begin read
  if [ -e "$OUT/reading.json" ] || [ -e "$OUT/reading.txt" ]; then
    echo "refused: a reading exists under $OUT and is never overwritten" >&2
    echo 95 > "$S/read.exit"
    return 95
  fi
  local started=()
  if [ -e "$S/run.exit" ]; then started=(--run-started); fi
  $PY "$SCRIPT" read --commit "$COMMIT" --runs "$R" --out "$OUT" ${started[@]+"${started[@]}"} > /dev/null 2>&1
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
  tests/test_bf_form.py tests/test_bf_read.py tests/test_bf_build.py tests/test_bf_launchers.py \
  tests/test_shape_interp.py tests/test_shape_interp_build.py tests/test_t2_build.py > /dev/null 2>&1
e=$?
[ "$e" = 0 ] || finish "$e"
echo 0 > "$S/tests.exit"

export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 NUMEXPR_NUM_THREADS=8
begin g0
pids=()
for s in "${SEEDS[@]}"; do
  (
    $PY "$SCRIPT" g0 --seed "$s" --commit "$COMMIT" --out "$R/g0" --data "$DATA" > /dev/null 2>&1
    echo $? > "$S/rc/g0-$s"
  ) &
  pids+=($!)
done
for p in "${pids[@]}"; do wait "$p"; done
bad=0
for s in "${SEEDS[@]}"; do
  [ "$(cat "$S/rc/g0-$s" 2>/dev/null)" = 0 ] || bad=1
done
[ $bad = 0 ] || finish 1
echo 0 > "$S/g0.exit"

begin g0read
$PY "$SCRIPT" g0read --commit "$COMMIT" --out "$R/g0" > /dev/null 2>&1
e=$?
if [ "$e" = 4 ]; then fail_then_read 4; fi
[ "$e" = 0 ] || finish "$e"
echo 0 > "$S/g0read.exit"

begin run
pids=()
for t in "${TASKS[@]}"; do
  (
    $PY "$SCRIPT" run --task "$t" --commit "$COMMIT" --runs "$R" --data "$DATA" > /dev/null 2>&1
    echo $? > "$S/rc/run-$t"
  ) &
  pids+=($!)
done
for p in "${pids[@]}"; do wait "$p"; done
bad=0
for t in "${TASKS[@]}"; do
  [ "$(cat "$S/rc/run-$t" 2>/dev/null)" = 0 ] || bad=1
done
[ $bad = 0 ] || fail_then_read 1
echo 0 > "$S/run.exit"

call_read
finish $?
