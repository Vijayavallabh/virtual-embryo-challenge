#!/usr/bin/env bash
# Row K25 (experiments/harness/TRIES.md; experiments/t3_decomp/k25_dispersion.py), one step per call, on a DGX campaign
# checkout (~/v/vec-campaign-c, never the reference checkout), the machine rows K18 and K21 to K24 ran on. Run from the
# checkout's root, detached (setsid nohup ... &), and wait on the step's exit file with a background poller; never a
# foreground sleep. CPU only, no GPU.
#
#   bash experiments/t3_decomp/run_k25.sh tests|pregate|run|read COMMIT
#
# COMMIT: the full hash of the run commit (the row's commit; the row naming the script is in TRIES.md there). Each step
# runs only if the step before it wrote exit 0 for the same COMMIT (else 96, refused), and refuses another HEAD (99),
# an unclean tracked tree (98) or a missing row (95). Every path writes ~/v/scratch/k25/<STEP>.exit (a call without
# its two arguments writes usage.exit, 2; an unknown step 97; a TERM, INT or HUP signal 89). Nothing is ever removed.
#   tests COMMIT    the unit tests (synthetic arrays only), 8 BLAS threads; exit 0 or 10
#   pregate COMMIT  the source-only pre-gate (no fold target, no try) into runs/k25/pregate.json: 0 passes; 25 fails
#                   (nothing runs, nothing counts; the script's exit 5); 22 the script refused (its exit 2); 20 any
#                   other failure. A pregate.json already there is read again (passes 0, fails 25), never remade
#   run COMMIT      seven parts in parallel: plain g1 and g2, each as seeds 0-4 and 5-9 (8 BLAS threads a process), the
#                   depth-matched world (8 workers by 2 threads, 16 BLAS threads in its set-up, as rows K18 and K21 to
#                   K24) and x2 as seeds 0-4 and 5-9 (8 threads each); each part's exit code in runs/k25/exit/<part>,
#                   its log in runs/k25/logs/<part>.log (never printed); log lines other than the expected kinds are
#                   counted, never printed; 0 when every part exits 0, else 30 (a seed already written is skipped, so
#                   calling the step again at the same commit resumes)
#   read COMMIT     the fixed reading into experiments/t3_decomp/results/k25_dispersion/ (output to /dev/null): 0 a
#                   reading (whatever it says, VOID included), 44 a reading is there already (never overwritten),
#                   43 not a reading (pregate.json missing), 41 a crash in the read
# NUMPY_MADVISE_HUGEPAGE=0; TMPDIR under ~/v/scratch/k25.
set -u
S="$HOME/v/scratch/k25"
mkdir -p "$S/tmp" || exit 90
if [ "$#" -lt 2 ] || [ -z "$1" ] || [ -z "$2" ]; then
  echo "usage: run_k25.sh tests|pregate|run|read COMMIT" >&2
  echo 2 > "$S/usage.exit"
  exit 2
fi
STEP="$1"
COMMIT="$2"
finish() { trap - TERM INT HUP; echo "$1" > "$S/$STEP.exit"; exit "$1"; }
case "$STEP" in
  tests) PREV="" ;;
  pregate) PREV="tests" ;;
  run) PREV="pregate" ;;
  read) PREV="run" ;;
  *)
    echo "unknown step $STEP" >&2
    finish 97 ;;
esac
trap 'finish 89' TERM INT HUP
if [ -n "$PREV" ]; then
  if [ "$(cat "$S/$PREV.exit" 2>/dev/null)" != 0 ] || [ "$(cat "$S/$PREV.head" 2>/dev/null)" != "$COMMIT" ]; then
    echo "refused: step $PREV has not exited 0 at $COMMIT" >&2
    finish 96
  fi
fi
if [ "$(git rev-parse HEAD 2>/dev/null)" != "$COMMIT" ]; then
  echo "refused: HEAD is not $COMMIT" >&2
  finish 99
fi
if [ -n "$(git status --porcelain --untracked-files=no 2>&1)" ]; then
  echo "refused: the tree is not clean" >&2
  finish 98
fi
grep -q '^| K25 .*`experiments/t3_decomp/k25_dispersion.py`' experiments/harness/TRIES.md \
  || { echo "refused: no TRIES.md row K25 names the script" >&2; finish 95; }
git rev-parse HEAD > "$S/$STEP.head"

export NUMPY_MADVISE_HUGEPAGE=0 TMPDIR="$S/tmp" UV_PROJECT_ENVIRONMENT="$PWD/.venv"
PY="uv run --frozen --no-sync python -I -B"
threads() { export OMP_NUM_THREADS=$1 OPENBLAS_NUM_THREADS=$1 MKL_NUM_THREADS=$1; }
threads 8
SCRIPT=experiments/t3_decomp/k25_dispersion.py
R=runs/k25
OUT=experiments/t3_decomp/results/k25_dispersion
RUN_DIR=$HOME/v/virtual-embryo-challenge/runs/t3_decomp/gata4_20261028T122425Z  # read only; fixed, no override
PASSES_CHECK='import json, sys; g = json.load(open(sys.argv[1])); print("pass" if g.get("passes") is True else "fail")'
OK_LINE='^K25 (plain g[12]|x2) seed [0-9]+: (scored|present, skipped)$|^K25 x2: every seed present, skipped$|^wrote [^ ]+\.json|^\[rep[12] seed [0-9]+\] K18 depth-matched: [0-9]+ arms and the check, [0-9]+s$|^self-check passed '

case "$STEP" in
  tests)
    rc=0
    $PY -m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp="$S/tmp/pytest" \
      tests/test_k25_dispersion.py experiments/t3_decomp/test_k9_residual.py tests/test_k18_built_files.py \
      tests/test_k21_factor_forms.py > "$S/tests.log" 2>&1 || rc=$?
    if [ $rc != 0 ]; then finish 10; fi
    finish 0 ;;
  pregate)
    if [ -e "$R/pregate.json" ]; then
      got=$($PY -c "$PASSES_CHECK" "$R/pregate.json" 2> /dev/null) || got=unreadable
      if [ "$got" = pass ]; then finish 0; fi
      if [ "$got" = fail ]; then finish 25; fi
      finish 20
    fi
    mkdir -p "$R"
    rc=0
    $PY "$SCRIPT" pregate --commit "$COMMIT" --out "$R" > "$S/pregate.log" 2>&1 || rc=$?
    if [ $rc = 0 ]; then finish 0; fi
    if [ $rc = 5 ]; then finish 25; fi
    if [ $rc = 2 ]; then finish 22; fi
    finish 20 ;;
  run)
    mkdir -p "$R/exit" "$R/logs"
    PARTS=("plain_g1_0-4" "plain_g1_5-9" "plain_g2_0-4" "plain_g2_5-9" "dm" "x2_0-4" "x2_5-9")
    pids=()
    for part in "${PARTS[@]}"; do
      (
        if [ "$part" = dm ]; then
          threads 16
          $PY "$SCRIPT" run --commit "$COMMIT" --world dm --workers 8 --threads 2 --out "$R" --run-dir "$RUN_DIR" \
            > "$R/logs/$part.log" 2>&1
        elif [ "${part%%_*}" = x2 ]; then
          threads 8
          $PY "$SCRIPT" run --commit "$COMMIT" --world x2 --seeds "${part#x2_}" --out "$R" > "$R/logs/$part.log" 2>&1
        else
          threads 8
          IFS=_ read -r _ task seeds <<< "$part"
          $PY "$SCRIPT" run --commit "$COMMIT" --world plain --task "$task" --seeds "$seeds" --out "$R" \
            --run-dir "$RUN_DIR" > "$R/logs/$part.log" 2>&1
        fi
        echo $? > "$R/exit/$part"
      ) &
      pids+=($!)
    done
    for p in "${pids[@]}"; do wait "$p" || true; done
    bad=0
    for part in "${PARTS[@]}"; do
      [ "$(cat "$R/exit/$part" 2>/dev/null)" = 0 ] || bad=1
      n=$(grep -cvE "$OK_LINE" "$R/logs/$part.log" 2>/dev/null) || true
      echo "${n:-unreadable}" > "$R/exit/$part.other_lines"
    done
    if [ $bad != 0 ]; then finish 30; fi
    finish 0 ;;
  read)
    if [ -e "$OUT/reading.json" ] || [ -e "$OUT/reading.txt" ]; then
      echo "refused: a reading is there already; it is never overwritten" >&2
      finish 44
    fi
    rc=0
    $PY "$SCRIPT" read --commit "$COMMIT" --runs "$R" --out "$OUT" > /dev/null 2>&1 || rc=$?
    if [ $rc = 0 ]; then finish 0; fi
    if [ $rc = 3 ]; then finish 43; fi
    if [ $rc = 4 ]; then finish 44; fi
    finish 41 ;;
esac
