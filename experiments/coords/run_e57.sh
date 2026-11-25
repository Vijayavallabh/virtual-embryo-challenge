#!/usr/bin/env bash
# Row E57 (experiments/harness/TRIES.md) on the DGX: step 0 (E57a), the reproduction gate, the run of RS (E57b,
# fgw_shape0x_ov_rs) and its placebo (E57c, fgw_shape0x_ov_rs_same), and the read, all in one call, from the root of a
# checkout or worktree of the campaign repository at the row's commit (HEAD = COMMIT, clean):
#   bash experiments/coords/run_e57.sh COMMIT RUN
# COMMIT is the run commit (full hash), the row's commit. RUN names the scratch folder ~/v/scratch/e57/RUN: letters,
# digits, '_' and '-' only (else exit 90); refused if it exists (exit 91). A call without both arguments exits 2.
# Environment: E57_DATA (default ~/v/vec-campaign/data/raw), E57_VENV (default ~/v/vec-campaign/.venv).
# Steps, in order; each writes RUN/<step>.exit, first the word running, then its exit code; RUN/done.exit holds
# running, then the call's code (0 when all exited 0). The call stops at the first step that does not exit 0:
#   preflight  HEAD = COMMIT (else 96), the tracked tree clean (98), docs/S2_RULES.sha256 checks (97), TRIES.md holds
#              exactly one row E57 whose first cell starts "E57 (" and names E57a, experiments/coords/e57_rs_read.py,
#              experiments/coords/run_e57.sh and the void clause's tolerance in the words "de_score and de_direction
#              more than 1e-4 off O's" (94; the read script's DE_TOL_TEXT), and neither results folder holds a reading
#              yet (93)
#   tests      the unit tests of the step, the builder option, the read and this launcher (synthetic arrays only)
#   step0      E57a, counted when it starts: e57_rs_read.py step0 for the ten pairs at COMMIT into RUN/runs/step0, one
#              process per pair, then step0-read into experiments/coords/results/e57a_step0/ (3 on STOP: E57b and
#              E57c never run, and the call ends there with 3)
#   gate       the floor and O under main against E55, one process per task (3: some task fails)
#   run        the floor, O, the arm and the placebo under both builds, one process per task (E57b and E57c are
#              counted when it starts; refused without step 0's PASS and the nine gate files)
#   read       the reading into experiments/coords/results/e57_rs/ (never overwritten: 95)
# When the gate or the run does not exit 0, the read is still called once (so the VOID reading, with its reasons and
# the sentence that a rerun needs a new row, is written) and the call exits with the failing step's code.
# BLAS threads: 3 a process for tests and step0 (ten at once), 8 for gate, run and read (nine at once), as run_e55.sh
# did, since the 0.01 guard compares with E55's run made at 8. NUMPY_MADVISE_HUGEPAGE=0; TMPDIR under RUN. Every
# process's output goes to /dev/null; per-process exit codes are RUN/rc/<step>-<name>. Nothing is deleted or uploaded.
set -u
if [ "$#" -ne 2 ] || [ -z "$1" ] || [ -z "$2" ]; then
  echo "usage: run_e57.sh COMMIT RUN" >&2
  exit 2
fi
COMMIT="$1"
RUN="$2"
case "$RUN" in
  *[!A-Za-z0-9_-]*) echo "refused: RUN takes letters, digits, _ and - only" >&2; exit 90 ;;
esac
case "$COMMIT" in
  *[!0-9a-f]*) echo "refused: COMMIT must be a full hash" >&2; exit 90 ;;
esac
[ "${#COMMIT}" = 40 ] || { echo "refused: COMMIT must be a full hash" >&2; exit 90; }
S="$HOME/v/scratch/e57/$RUN"
mkdir -p "$HOME/v/scratch/e57" || exit 92
mkdir "$S" 2> /dev/null || { echo "refused: the run folder $S exists (a run folder is never reused)" >&2; exit 91; }
mkdir -p "$S/tmp" "$S/rc" "$S/runs" || exit 92
echo running > "$S/done.exit"
STEP=preflight
finish() { echo "$1" > "$S/$STEP.exit"; echo "$1" > "$S/done.exit"; exit "$1"; }
begin() { STEP="$1"; echo running > "$S/$STEP.exit"; }
begin preflight
export OMP_NUM_THREADS=3 OPENBLAS_NUM_THREADS=3 MKL_NUM_THREADS=3 NUMEXPR_NUM_THREADS=3 NUMPY_MADVISE_HUGEPAGE=0
export TMPDIR="$S/tmp"
export PATH="$HOME/.local/bin:$PATH"
export UV_PROJECT_ENVIRONMENT="${E57_VENV:-$HOME/v/vec-campaign/.venv}"
DATA="${E57_DATA:-$HOME/v/vec-campaign/data/raw}"
PY="uv run --frozen --no-sync python -I -B"
HEAD=$(git rev-parse HEAD 2> /dev/null) || finish 96
if [ "$HEAD" != "$COMMIT" ]; then
  echo "refused: HEAD is not $COMMIT" >&2
  finish 96
fi
ST=$(git status --porcelain --untracked-files=no 2> /dev/null) || finish 98
if [ -n "$ST" ]; then
  echo "refused: the tracked tree is not clean" >&2
  finish 98
fi
sha256sum -c docs/S2_RULES.sha256 > /dev/null 2>&1 || { echo "refused: docs/S2_RULES.sha256 does not check" >&2;
                                                        finish 97; }
python3 -I -B -c 'import sys
rows = [ln for ln in open(sys.argv[1], encoding="utf-8").read().splitlines() if ln.startswith("| E57 ")]
first = rows[0][2:].split("|", 1)[0] if len(rows) == 1 else ""
ok = (first.startswith("E57 (") and "E57a" in first and "experiments/coords/e57_rs_read.py" in first
      and "experiments/coords/run_e57.sh" in first and "de_score and de_direction more than 1e-4 off O\x27s" in first)
sys.exit(0 if ok else 1)' experiments/harness/TRIES.md 2> /dev/null \
  || { echo "refused: TRIES.md does not hold exactly one row E57 naming E57a, the script, this launcher and the" \
            "de tolerance" >&2
       finish 94; }
SCRIPT=experiments/coords/e57_rs_read.py
STEP0_OUT=experiments/coords/results/e57a_step0
OUT=experiments/coords/results/e57_rs
for d in "$STEP0_OUT" "$OUT"; do
  if [ -e "$d/reading.json" ] || [ -e "$d/reading.txt" ]; then
    echo "refused: $d already holds a reading (a reading is never overwritten)" >&2
    finish 93
  fi
done
echo "$HEAD" > "$S/head"
echo 0 > "$S/preflight.exit"
R="$S/runs"
TASKS=(half late new1 new2 e75 e725 x1 x2 x3)
PAIRS=(test half late new1 new2 e75 e725 x1 x2 x3)

parallel() {  # $1 = step: one process per task, all nine at once; 0 all exit 0, 3 a gate failure, 1 otherwise
  local pids=() t bad=0 fail=0 e
  for t in "${TASKS[@]}"; do
    (
      if [ "$1" = run ]; then
        $PY "$SCRIPT" run --task "$t" --commit "$COMMIT" --runs "$R" --data "$DATA" --step0 "$STEP0_OUT" \
          > /dev/null 2>&1
      else
        $PY "$SCRIPT" "$1" --task "$t" --commit "$COMMIT" --runs "$R" --data "$DATA" > /dev/null 2>&1
      fi
      echo $? > "$S/rc/$1-$t"
    ) &
    pids+=($!)
  done
  for p in "${pids[@]}"; do wait "$p"; done
  for t in "${TASKS[@]}"; do
    e=$(cat "$S/rc/$1-$t" 2> /dev/null)
    if [ "$e" = 3 ]; then fail=1; elif [ "$e" != 0 ]; then bad=1; fi
  done
  if [ $bad = 1 ]; then return 1; fi
  if [ $fail = 1 ]; then return 3; fi
  return 0
}

call_read() {  # the read step: writes reading.json and reading.txt (never overwritten); its code
  begin read
  export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 NUMEXPR_NUM_THREADS=8
  if [ -e "$OUT/reading.json" ] || [ -e "$OUT/reading.txt" ]; then
    echo "refused: a reading exists under $OUT and is never overwritten" >&2
    echo 95 > "$S/read.exit"
    return 95
  fi
  $PY "$SCRIPT" read --commit "$COMMIT" --runs "$R" --out "$OUT" --step0 "$STEP0_OUT" > /dev/null 2>&1
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
  tests/test_run_sorted.py tests/test_e57_rs_read.py tests/test_embryo_rs_launcher.py \
  tests/test_series_back_build.py tests/test_e51_so_ohi_read.py tests/test_shape_interp.py \
  tests/test_shape_interp_build.py tests/test_t2_build.py > /dev/null 2>&1
e=$?
[ "$e" = 0 ] || finish "$e"
echo 0 > "$S/tests.exit"

begin step0
mkdir -p "$R/step0"
pids=()
for x in "${PAIRS[@]}"; do
  (
    $PY "$SCRIPT" step0 --name "$x" --commit "$COMMIT" --runs "$R" --data "$DATA" > /dev/null 2>&1
    echo $? > "$S/rc/step0-$x"
  ) &
  pids+=($!)
done
for p in "${pids[@]}"; do wait "$p"; done
$PY "$SCRIPT" step0-read --commit "$COMMIT" --runs "$R" --out "$STEP0_OUT" > /dev/null 2>&1
e=$?
echo "$e" > "$S/rc/step0-read"
[ "$e" = 0 ] || finish "$e"
echo 0 > "$S/step0.exit"

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
