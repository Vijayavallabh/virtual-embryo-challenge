#!/usr/bin/env bash
# Row K26 (experiments/harness/TRIES.md; experiments/t3_decomp/k26_query_map.py), one step per call, on a DGX campaign
# checkout (never the reference checkout). Run from the checkout's root, detached (setsid nohup ... &), and wait on the
# step's exit file with a background poller, never a foreground sleep. CPU only.
#
#   bash experiments/t3_decomp/run_k26.sh STEP COMMIT
#
# COMMIT: the full hash of the row's commit (the row is in TRIES.md there). Steps, in order; each runs only if the step
# before it wrote exit 0 for the same COMMIT (else 96):
#   tests     the new and the related tests, synthetic data only (exit 10 on a failure)
#   features  TARGET-FREE: T1 contexts, detection deciles, the random draws; R/features.json (22 refused, 24 exists,
#             20 other)
#   stage0    the COUNTED DIAGNOSTIC on outside perturbations only; R/stage0.json; 27 the gate fails (STOP: no fold
#             read; nothing after this step runs), 22 refused, 20 other
#   pregate   TARGET-FREE: the fold tables and their drops; R/pregate.json; 25 the arm is dropped (STOP: nothing runs),
#             22 refused, 24 exists, 20 other
#   guard     the guard part (base arms only) in seven parallel parts, then the reproduction gate; 26 the gate fails
#             (the VOID reading is written into OUT first; a rerun finishes 26 again and never reads twice), 30 a part
#             did not exit 0, 22 refused, 20 other
#   run       the scoring parts in seven parallel parts (plain g1 and g2 as seeds 0-4 and 5-9, dm with 8 workers by 2
#             threads, x2 as seeds 0-4 and 5-9); a seed whose file exists is skipped, so a killed run resumes; 30 a
#             part did not exit 0
#   read      the reading into experiments/t3_decomp/results/k26_query_map/ (output to /dev/null): 0 a reading (VOID
#             included), 43 INCOMPLETE, 44 a reading is there already (never overwritten), 41 a crash
# Every step first refuses: HEAD not COMMIT (99), a tracked file changed (98), no TRIES.md row naming the script (95),
# an input missing (94: the outside folders and files below are not real, TABLES_K24.sha256 untracked, or
# EXTERNAL_DATA.md rows X013, X014, X015 and X017 do not each name K26). A call without two arguments writes
# usage.exit (2); an unknown step 97. Each step writes $S/<STEP>.exit on every path and $S/<STEP>.head. Logs go to
# R/logs/ and are never printed; their lines other than the script's step lines are counted into
# R/exit/<part>.other_lines. Nothing is ever removed by this launcher. NUMPY_MADVISE_HUGEPAGE=0; TMPDIR under $S.
set -u
if [ "$#" -lt 2 ] || [ -z "$1" ] || [ -z "$2" ]; then
  mkdir -p "$HOME/v/scratch/k26" && echo 2 > "$HOME/v/scratch/k26/usage.exit"
  echo "usage: run_k26.sh tests|features|stage0|pregate|guard|run|read COMMIT" >&2
  exit 2
fi
STEP="$1"
COMMIT="$2"
S="$HOME/v/scratch/k26_${COMMIT:0:7}"
mkdir -p "$S/tmp" || exit 2
finish() { echo "$1" > "$S/$STEP.exit"; exit "$1"; }
case "$STEP" in
  tests) PREV="" ;;
  features) PREV="tests" ;;
  stage0) PREV="features" ;;
  pregate) PREV="stage0" ;;
  guard) PREV="pregate" ;;
  run) PREV="guard" ;;
  read) PREV="run" ;;
  *)
    echo "unknown step $STEP" >&2
    finish 97 ;;
esac
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
if [ -n "$(git status --porcelain --untracked-files=no 2>/dev/null)" ]; then
  echo "refused: the tree is not clean" >&2
  finish 98
fi
grep -q '^| K26 .*`experiments/t3_decomp/k26_query_map.py`' experiments/harness/TRIES.md 2>/dev/null \
  || { echo "refused: no TRIES.md row K26 names the script" >&2; finish 95; }

TABLES=$HOME/v/vec-ext/t3/tables  # sc3 and bulk; fixed, no override
TABLES_TWO=$HOME/v/vec-ext/t3/tables_two_study  # `two`; fixed, no override
PROCESSED=$HOME/v/vec-ext/t3/processed  # the per-gene unit files; fixed, no override
REFSTATS=$HOME/v/vec-ext/t3/reference/refstats.npz  # the units' reference statistics; fixed, no override
RUN_DIR=$HOME/v/virtual-embryo-challenge/runs/t3_decomp/gata4_20261028T122425Z  # read only; fixed, no override
for d in "$TABLES" "$TABLES_TWO" "$PROCESSED"; do
  { [ -d "$d" ] && [ ! -L "$d" ]; } || { echo "refused: $d is not a real directory" >&2; finish 94; }
done
{ [ -f "$REFSTATS" ] && [ ! -L "$REFSTATS" ]; } || { echo "refused: $REFSTATS is not a plain file" >&2; finish 94; }
GIT_OPTIONAL_LOCKS=0 git ls-files --error-unmatch experiments/t3_ext/TABLES_K24.sha256 > /dev/null 2>&1 \
  || { echo "refused: TABLES_K24.sha256 is not tracked" >&2; finish 94; }
for x in X013 X014 X015 X017; do
  grep "^| $x " EXTERNAL_DATA.md 2>/dev/null | grep -qF K26 \
    || { echo "refused: EXTERNAL_DATA.md row $x does not name K26" >&2; finish 94; }
done
git rev-parse HEAD > "$S/$STEP.head"

export NUMPY_MADVISE_HUGEPAGE=0 TMPDIR="$S/tmp" UV_PROJECT_ENVIRONMENT="$PWD/.venv" VEC_T3_EXT_TABLES="$TABLES" \
  VEC_T3_EXT_TABLES_TWO="$TABLES_TWO" VEC_T3_EXT_PROCESSED="$PROCESSED" VEC_T3_EXT_REFSTATS="$REFSTATS"
PY="uv run --frozen --no-sync python -I -B"
threads() { export OMP_NUM_THREADS=$1 OPENBLAS_NUM_THREADS=$1 MKL_NUM_THREADS=$1; }
threads 8
SCRIPT=experiments/t3_decomp/k26_query_map.py
R=runs/k26_${COMMIT:0:7}
OUT=experiments/t3_decomp/results/k26_query_map
OK_LINE='^(K26( guard)? (plain|dm|x2) [a-z0-9]+ seed [0-9]+: (scored|present, skipped)|wrote [^ ]+|\[rep[12] seed [0-9]+\] K18 depth-matched: [0-9]+ arms and the check, [0-9]+s|self-check passed .*)$'
PARTS=("plain_g1_0-4" "plain_g1_5-9" "plain_g2_0-4" "plain_g2_5-9" "dm" "x2_x2_0-4" "x2_x2_5-9")
mkdir -p "$R/exit" "$R/logs"

parts() {  # $1 = guard|score: the seven parts at once; returns 1 when a part did not exit 0
  local pids=() part p bad=0 n
  for part in "${PARTS[@]}"; do
    (
      if [ "$part" = dm ]; then
        threads 16
        $PY "$SCRIPT" run --part "$1" --world dm --workers 8 --threads 2 --out "$R" --run-dir "$RUN_DIR" \
          --commit "$COMMIT" > "$R/logs/$1_$part.log" 2>&1
      else
        IFS=_ read -r world task seeds <<< "$part"
        if [ "$world" = plain ]; then
          $PY "$SCRIPT" run --part "$1" --world plain --task "$task" --seeds "$seeds" --out "$R" --run-dir "$RUN_DIR" \
            --commit "$COMMIT" > "$R/logs/$1_$part.log" 2>&1
        else
          $PY "$SCRIPT" run --part "$1" --world x2 --seeds "$seeds" --out "$R" --run-dir "$RUN_DIR" \
            --commit "$COMMIT" > "$R/logs/$1_$part.log" 2>&1
        fi
      fi
      echo $? > "$R/exit/$1_$part"
    ) &
    pids+=($!)
  done
  for p in "${pids[@]}"; do wait "$p" || true; done
  for part in "${PARTS[@]}"; do
    [ "$(cat "$R/exit/$1_$part" 2>/dev/null)" = 0 ] || bad=1
    n=$(grep -cvE "$OK_LINE" "$R/logs/$1_$part.log" 2>/dev/null) || true
    echo "${n:-unreadable}" > "$R/exit/$1_$part.other_lines"
  done
  return $bad
}

void_read() {  # the VOID reading after a failed reproduction gate, made once
  local rc=0
  if [ -e "$OUT/reading.json" ]; then finish 26; fi
  $PY "$SCRIPT" read --runs "$R" --out "$OUT" --commit "$COMMIT" > /dev/null 2>&1 || rc=$?
  echo $rc > "$S/void_read.exit"
  if [ $rc = 3 ]; then finish 43; fi
  if [ $rc = 4 ]; then finish 44; fi
  if [ $rc != 0 ]; then finish 41; fi
  finish 26
}

rc=0
case "$STEP" in
  tests)
    $PY -m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp="$S/tmp/pytest" \
      tests/test_k26_query_map.py tests/test_k26_context_map.py tests/test_k26_launcher.py \
      tests/test_k24_wider_tables.py tests/test_k18_built_files.py tests/test_k22_outside_forms.py \
      > "$S/tests.log" 2>&1 || rc=$?
    if [ $rc != 0 ]; then finish 10; fi
    finish 0 ;;
  features)
    $PY "$SCRIPT" features --out "$R" --commit "$COMMIT" > "$R/logs/features.log" 2>&1 || rc=$?
    case $rc in 0) finish 0 ;; 2) finish 22 ;; 4) finish 24 ;; *) finish 20 ;; esac ;;
  stage0)
    $PY "$SCRIPT" stage0 --out "$R" --commit "$COMMIT" > "$R/logs/stage0.log" 2>&1 || rc=$?
    case $rc in 0) finish 0 ;; 7) finish 27 ;; 2) finish 22 ;; *) finish 20 ;; esac ;;
  pregate)
    $PY "$SCRIPT" pregate --out "$R" --run-dir "$RUN_DIR" --commit "$COMMIT" > "$R/logs/pregate.log" 2>&1 || rc=$?
    case $rc in 0) finish 0 ;; 5) finish 25 ;; 7) finish 27 ;; 2) finish 22 ;; 4) finish 24 ;; *) finish 20 ;; esac ;;
  guard)
    if [ ! -e "$R/guard.json" ]; then
      parts guard || finish 30
      $PY "$SCRIPT" guard --runs "$R" --commit "$COMMIT" > "$R/logs/guard.log" 2>&1 || rc=$?
      case $rc in 0) finish 0 ;; 6) void_read ;; 2) finish 22 ;; 3) finish 30 ;; *) finish 20 ;; esac
    fi
    gpass=$($PY -c 'import json, sys; g = json.load(open(sys.argv[1])); print("pass" if g.get("passes") is True else "fail")' \
      "$R/guard.json" 2> /dev/null) || gpass=unreadable
    if [ "$gpass" = pass ]; then finish 0; fi
    if [ "$gpass" = fail ]; then void_read; fi
    finish 20 ;;
  run)
    parts score || finish 30
    finish 0 ;;
  read)
    $PY "$SCRIPT" read --runs "$R" --out "$OUT" --commit "$COMMIT" > /dev/null 2>&1 || rc=$?
    case $rc in 0) finish 0 ;; 3) finish 43 ;; 4) finish 44 ;; *) finish 41 ;; esac ;;
esac
finish 20
