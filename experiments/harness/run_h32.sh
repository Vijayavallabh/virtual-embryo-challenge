#!/usr/bin/env bash
# Row H32 (experiments/harness/TRIES.md; experiments/harness/h32_heart_gated_trend.py), one step per call, on a DGX
# campaign checkout or a detached worktree of one (never the reference checkout). It runs from the root of the checkout
# it is called in, with that checkout's own environment ($PWD/.venv after the cd), and its data linked there. Run it
# detached and wait on the exit file with a background poller; never a foreground sleep. CPU only. Steps, in order;
# each runs only if the step before it wrote exit 0 for the same ROW (and, for read, the same RUN) commit:
#   tests ROW        the unit tests of the script and of rows H30, H28 and H27 (synthetic arrays and hand-made tables)
#   tables ROW       row H30's two chord tables (h30.build_table, unchanged) and their sha256 into runs/h32/tables/; where
#                    row H30's tables step wrote its tables (runs/h30/tables) and their sha256 is committed and the reader
#                    is the same, those files are used instead of a rebuild (tables_origin.json says which)
#   gate ROW         G0, target-free (counts no try): runs/h32/gate/gate.json and gate.txt; exit 0 PASS, 3 FAIL (nothing
#                    further runs or counts), anything else an error. Then, by hand, before run: commit
#                    runs/h32/gate/gate.json and gate.txt as experiments/harness/results/h32_heart_gated_trend/gate/ and
#                    runs/h32/tables/tables_sha256.json as experiments/harness/results/h32_heart_gated_trend/
#                    tables_sha256.json, byte for byte, and nothing else: that commit is RUN.
#   run ROW RUN      the ten parts: per backtest five processes of two seeds, the two backtests one after the other;
#                    refused unless HEAD is RUN, RUN descends from ROW and differs from it in exactly those three files,
#                    each byte for byte the step's own file (exit 95). The script prints "<task> seed <s> (<sec> s) done"
#                    and "wrote <path>" only; any other log line is counted, never printed (runs/h32/exit/*.other_lines).
#   read ROW RUN     the one read into experiments/harness/results/h32_heart_gated_trend (reading.json, reading.txt; the
#                    script refuses to overwrite either); its two printed lines in read.out
# ROW is the commit of row H32 (with the X002 and X003 amendments naming it), RUN the run commit (full hashes). Every
# step also refuses while TRIES.md has no row H32 (95) or EXTERNAL_DATA.md's X002 and X003 lack the wording the script
# requires (93). A step that passes its checks first sets ~/v/scratch/h32/<STEP>.exit to "running", then records ROW
# in <STEP>.head, runs, and writes its exit code. A call refused before its step runs (usage 2, unknown step 97, earlier
# step not done 96, the step already exited 0 at this ROW 94, RUN not as above 95, unclean tree 98, another HEAD 99, no
# row 95, register wording 93) writes <STEP>.refused.exit only (usage.refused.exit without a step), so a refusal never
# overwrites a finished step's exit. 8 BLAS threads a process, NUMPY_MADVISE_HUGEPAGE=0, TMPDIR under the scratch
# folder. The launcher deletes nothing itself; pytest clears its own --basetemp folder (~/v/scratch/h32/tmp/pytest).
set -u
S="$HOME/v/scratch/h32"
mkdir -p "$S/tmp"
if [ "$#" -lt 2 ] || [ -z "$1" ] || [ -z "$2" ]; then
  echo "usage: run_h32.sh tests|tables|gate ROW, or run|read ROW RUN" >&2
  echo 2 > "$S/usage.refused.exit"
  exit 2
fi
STEP="$1"
ROW="$2"
RUN="${3:-}"
refuse() { echo "$1" > "$S/$STEP.refused.exit"; exit "$1"; }
finish() { echo "$1" > "$S/$STEP.exit"; exit "$1"; }
case "$STEP" in
  tests) PREV="" ;;
  tables) PREV="tests" ;;
  gate) PREV="tables" ;;
  run) PREV="gate" ;;
  read) PREV="run" ;;
  *)
    echo "unknown step $STEP" >&2
    refuse 97 ;;
esac
case "$STEP" in
  run|read)
    if [ -z "$RUN" ]; then
      echo "usage: run_h32.sh $STEP ROW RUN" >&2
      refuse 2
    fi
    HEAD_WANT="$RUN" ;;
  *) HEAD_WANT="$ROW" ;;
esac
TOP="$(git rev-parse --show-toplevel 2>/dev/null)"
if [ -z "$TOP" ] || ! cd "$TOP"; then
  echo "refused: not inside a checkout" >&2
  refuse 99
fi
if [ -n "$PREV" ]; then
  if [ "$(cat "$S/$PREV.exit" 2>/dev/null)" != 0 ] || [ "$(cat "$S/$PREV.head" 2>/dev/null)" != "$ROW" ]; then
    echo "refused: step $PREV has not exited 0 at $ROW" >&2
    refuse 96
  fi
  if [ "$STEP" = read ] && [ "$(cat "$S/run.run" 2>/dev/null)" != "$RUN" ]; then
    echo "refused: step run has not exited 0 at $RUN" >&2
    refuse 96
  fi
fi
if [ "$(cat "$S/$STEP.exit" 2>/dev/null)" = 0 ] && [ "$(cat "$S/$STEP.head" 2>/dev/null)" = "$ROW" ]; then
  echo "refused: step $STEP has already exited 0 at $ROW" >&2
  refuse 94
fi
if [ "$(git rev-parse HEAD)" != "$HEAD_WANT" ]; then
  echo "refused: HEAD is not $HEAD_WANT" >&2
  refuse 99
fi
if [ -n "$(git status --porcelain)" ]; then
  echo "refused: the tree is not clean" >&2
  refuse 98
fi
if ! grep -q '^| H32 ' experiments/harness/TRIES.md; then
  echo "refused: TRIES.md has no row H32" >&2
  refuse 95
fi
H32R=experiments/harness/results/h32_heart_gated_trend
if [ -n "$RUN" ]; then
  ALLOWED="$H32R/gate/gate.json $H32R/gate/gate.txt $H32R/tables_sha256.json"
  CHANGED="$(git diff --name-only "$ROW" "$RUN" | sort | tr '\n' ' ')"
  WANT="$(printf '%s\n' $ALLOWED | sort | tr '\n' ' ')"
  if [ "$RUN" = "$ROW" ] || ! git merge-base --is-ancestor "$ROW" "$RUN" || [ "$CHANGED" != "$WANT" ]; then
    echo "refused: RUN must descend from ROW and differ from it in the gate reading and the table checksums only" >&2
    refuse 95
  fi
  if ! cmp -s runs/h32/gate/gate.json "$H32R/gate/gate.json" \
      || ! cmp -s runs/h32/gate/gate.txt "$H32R/gate/gate.txt" \
      || ! cmp -s runs/h32/tables/tables_sha256.json "$H32R/tables_sha256.json"; then
    echo "refused: the committed gate reading or table checksums are not the steps' own files" >&2
    refuse 95
  fi
fi
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 NUMPY_MADVISE_HUGEPAGE=0 TMPDIR="$S/tmp"
export UV_PROJECT_ENVIRONMENT="$PWD/.venv"
PY="uv run --frozen --no-sync python -I -B"
SCRIPT=experiments/harness/h32_heart_gated_trend.py
if ! $PY "$SCRIPT" check > "$S/$STEP.check.log" 2>&1; then
  echo "refused: row H32 or the register wording the script requires is missing" >&2
  refuse 93
fi
echo running > "$S/$STEP.exit"
echo "$ROW" > "$S/$STEP.head"
if [ -n "$RUN" ]; then echo "$RUN" > "$S/$STEP.run"; fi
R=runs/h32
TASKS=("T2h/<=E9.5->E10.5:e95" "T2h/<=E8.75->E10.5:e875")
DONE_LINE='^T2h/[^ ]+ seed [0-9]+ \([0-9.]+ s\) done$'
WROTE_LINE='^wrote [^ ]+/h32\.json$'

case "$STEP" in
  tests)
    $PY -m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp="$S/tmp/pytest" \
      tests/test_h32_heart_gated_trend.py tests/test_h30_heart_late_chord.py tests/test_h28_heart_trend_strength.py \
      tests/test_h27_heart_forecast_size.py > "$S/tests.log" 2>&1
    finish $? ;;
  tables)
    mkdir -p "$R"
    $PY "$SCRIPT" tables --out "$R/tables" --h30 runs/h30/tables > "$S/tables.log" 2>&1
    finish $? ;;
  gate)
    $PY "$SCRIPT" gate --tables "$R/tables" --out "$R/gate" > "$S/gate.log" 2>&1
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
          $PY "$SCRIPT" run --task "$key" --seeds "$seeds" --tables "$R/tables" --out "$R/$short/$seeds" \
            > "$R/logs/$short-$seeds.log" 2>&1
          echo $? > "$R/exit/$short-$seeds"
        ) &
        pids+=($!)
      done
      for p in "${pids[@]}"; do wait "$p"; done
      for seeds in 0-1 2-3 4-5 6-7 8-9; do
        [ "$(cat "$R/exit/$short-$seeds" 2>/dev/null)" = 0 ] || bad=1
        n=$(grep -cvE "$DONE_LINE|$WROTE_LINE" "$R/logs/$short-$seeds.log" 2>/dev/null)
        echo "${n:-unreadable}" > "$R/exit/$short-$seeds.other_lines"
      done
    done
    finish $bad ;;
  read)
    $PY "$SCRIPT" read --runs "$R" --commit "$RUN" --out "$H32R" > "$S/read.out" 2>&1
    finish $? ;;
esac
