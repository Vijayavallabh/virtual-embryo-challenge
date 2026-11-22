#!/usr/bin/env bash
# Row H29b (experiments/harness/TRIES.md), the practice read of the heart's organ-cut file D (row R21), one step per
# call, on a DGX campaign checkout (~/v/vec-campaign, never the reference checkout). Run from the root of the checkout
# or detached worktree it is called in, detached, and wait on the exit file with a background poller; never a
# foreground sleep. CPU only. Steps, in order; each runs only if the step before it wrote exit 0 at the same COMMIT:
#   tests COMMIT   the unit tests of the script, of its launchers, of row H27's script and of the organ cut (synthetic
#                  arrays, hand-made tables and throwaway repositories only)
#   gate COMMIT    the target-free gate, not a try: E10.5 only, no target and no score, into runs/h29b/COMMIT/gate
#                  (gate.json, gate.txt); exit 0 PASS, 3 STOP (D is not a new form: nothing further runs and nothing
#                  counts), 4 REFUSED (an unmapped label: the config is fixed first, in a new commit), else an error
#   run COMMIT     the parts: per transition five processes of two seeds, the three transitions one after another,
#                  into runs/h29b/COMMIT/<short>/<seeds>/h29b.json; each part refuses before anything loads unless row
#                  H29b is in TRIES.md, the tree is clean, HEAD is COMMIT, the gate of COMMIT said PASS and its output
#                  does not exist. The script prints "<task> seed <s> (<sec> s) done" and "wrote <path>" only; any
#                  other log line is counted, never printed (runs/h29b/COMMIT/exit/*.other_lines and other_lines)
#   read COMMIT    the one read into experiments/harness/results/h29b_heart_organ (reading.json, reading.txt; refused,
#                  exit 95, if either exists); its two printed lines in read.out
# COMMIT is row H29b's run commit (full hash). A call refused before its step runs (usage 2, unknown step 97, earlier
# step not done 96, the step already exited 0 at this COMMIT 94, no row H29b 93, unclean tree 98, another HEAD 99)
# writes ~/v/scratch/h29b/<STEP>.refused.exit only (usage.refused.exit without a step, unknown.refused.exit for an
# unknown step), so a refusal never touches a
# finished step's exit. A step that passes its checks first overwrites ~/v/scratch/h29b/<STEP>.exit with the word
# running, then writes <STEP>.head, then runs, then writes its exit code there, so a step killed part way leaves a
# non-zero exit file beside its new .head and the next step refuses (96). 8 BLAS threads a process,
# NUMPY_MADVISE_HUGEPAGE=0, TMPDIR under the scratch folder. The launcher deletes nothing.
set -u
S="$HOME/v/scratch/h29b"
mkdir -p "$S/tmp"
if [ "$#" -lt 2 ] || [ -z "$1" ] || [ -z "$2" ]; then
  echo "usage: run_h29b.sh tests|gate|run|read COMMIT" >&2
  echo 2 > "$S/usage.refused.exit"
  exit 2
fi
STEP="$1"
COMMIT="$2"
refuse() { echo "$1" > "$S/$STEP.refused.exit"; exit "$1"; }
finish() { echo "$1" > "$S/$STEP.exit"; exit "$1"; }
case "$STEP" in
  tests) PREV="" ;;
  gate) PREV="tests" ;;
  run) PREV="gate" ;;
  read) PREV="run" ;;
  *)
    echo "unknown step $STEP" >&2
    echo 97 > "$S/unknown.refused.exit"
    exit 97 ;;
esac
if [ -n "$PREV" ]; then
  if [ "$(cat "$S/$PREV.exit" 2>/dev/null)" != 0 ] || [ "$(cat "$S/$PREV.head" 2>/dev/null)" != "$COMMIT" ]; then
    echo "refused: step $PREV has not exited 0 at $COMMIT" >&2
    refuse 96
  fi
fi
if [ "$(cat "$S/$STEP.exit" 2>/dev/null)" = 0 ] && [ "$(cat "$S/$STEP.head" 2>/dev/null)" = "$COMMIT" ]; then
  echo "refused: step $STEP has already exited 0 at $COMMIT" >&2
  refuse 94
fi
if [ "$(git rev-parse HEAD)" != "$COMMIT" ]; then
  echo "refused: HEAD is not $COMMIT (give the full hash)" >&2
  refuse 99
fi
if [ -n "$(git status --porcelain)" ]; then
  echo "refused: the tree is not clean" >&2
  refuse 98
fi
if ! grep -q '^| H29b ' experiments/harness/TRIES.md; then
  echo "refused: row H29b is not in experiments/harness/TRIES.md" >&2
  refuse 93
fi
OUT=experiments/harness/results/h29b_heart_organ
if [ "$STEP" = read ] && { [ -e "$OUT/reading.json" ] || [ -e "$OUT/reading.txt" ]; }; then
  echo "refused: $OUT already holds a reading" >&2
  refuse 95
fi
echo running > "$S/$STEP.exit"  # past the checks, before anything else: a killed step never leaves an earlier 0 here
echo "$COMMIT" > "$S/$STEP.head"
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 NUMPY_MADVISE_HUGEPAGE=0 TMPDIR="$S/tmp"
export UV_PROJECT_ENVIRONMENT="$PWD/.venv"
PY="uv run --frozen --no-sync python -I -B"
SCRIPT=experiments/harness/h29b_heart_organ.py
R="runs/h29b/$COMMIT"  # one folder per commit: no part or gate of another commit is skipped or read
TASKS=("T2h/<=E9.5->E10.5:e95" "T2h/E8.25_late->E8.5:e825" "T2h/<=E8.75->E10.5:e875")
DONE_LINE='^T2h/[^ ]+ seed [0-9]+ \([0-9.]+ s\) done$'
WROTE_LINE='^wrote [^ ]+/h29b\.json$'

case "$STEP" in
  tests)
    $PY -m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp="$S/tmp/pytest" \
      tests/test_h29b_heart_organ.py tests/test_heart_organ_launchers.py tests/test_h27_heart_forecast_size.py \
      tests/test_organ_scope_mask.py tests/test_t2_chamber_classes.py > "$S/tests.log" 2>&1
    finish $? ;;
  gate)
    mkdir -p "$R"
    $PY "$SCRIPT" gate --out "$R/gate" --commit "$COMMIT" > "$S/gate.out" 2>&1
    finish $? ;;
  run)
    mkdir -p "$R/logs" "$R/exit"
    bad=0
    for entry in "${TASKS[@]}"; do
      key=${entry%:*}
      short=${entry##*:}
      pids=()
      for seeds in 0-1 2-3 4-5 6-7 8-9; do
        (
          $PY "$SCRIPT" run --task "$key" --seeds "$seeds" --gate "$R/gate" --commit "$COMMIT" \
            --out "$R/$short/$seeds" > "$R/logs/$short-$seeds.log" 2>&1
          echo $? > "$R/exit/$short-$seeds"
        ) &
        pids+=($!)
      done
      for p in "${pids[@]}"; do wait "$p"; done
      for seeds in 0-1 2-3 4-5 6-7 8-9; do
        [ "$(cat "$R/exit/$short-$seeds" 2>/dev/null)" = 0 ] || bad=1
        [ -f "$R/$short/$seeds/h29b.json" ] || bad=1
        # imported code and Python warnings write to the same log: count the lines other than the two expected kinds
        # and keep only the count; the lines are never printed
        n=$(grep -cvE "$DONE_LINE|$WROTE_LINE" "$R/logs/$short-$seeds.log" 2>/dev/null)
        echo "${n:-unreadable}" > "$R/exit/$short-$seeds.other_lines"
      done
    done
    cat "$R"/exit/*.other_lines | awk '{ if ($1 ~ /^[0-9]+$/) s += $1; else u = 1 } END { print (u ? "unreadable" : s) }' \
      > "$R/exit/other_lines"
    finish $bad ;;
  read)
    $PY "$SCRIPT" read --runs "$R" --gate "$R/gate" --commit "$COMMIT" --out "$OUT" > "$S/read.out" 2>&1
    finish $? ;;
esac
