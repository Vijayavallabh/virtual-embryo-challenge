#!/usr/bin/env bash
# Row K26 (experiments/harness/TRIES.md; experiments/t3_decomp/k26_query_map.py), one step per call, on a DGX campaign
# checkout or a detached worktree (never the reference checkout), from that folder's root and with that folder's own
# environment ($PWD/.venv) and data ($PWD/data/raw, passed to the script). Run detached (setsid nohup ... &) and wait on
# the step's exit file with a background poller, never a foreground sleep. CPU only.
#
#   bash experiments/t3_decomp/run_k26.sh STEP COMMIT
#
# COMMIT: the full hash of the row's commit (the row is in TRIES.md there). Steps, in order; each runs only if the step
# before it wrote exit 0 for the same COMMIT:
#   tests     the new and the related tests, synthetic data only (exit 10 on a failure)
#   features  TARGET-FREE: T1 contexts, detection deciles, the random draws; R/features.json (22 refused, 24 exists,
#             20 other)
#   stage0    the COUNTED DIAGNOSTIC on outside perturbations only; R/stage0.json; 27 the gate fails (STOP: the
#             reading of the stop is written into OUT, no fold is read, nothing after this step runs), 22 refused,
#             20 other
#   pregate   TARGET-FREE (no fold's table is compared with a table built from that fold's own target file; x2 leaves
#             the barred factor out): the fold tables and their drops; R/pregate.json; 25 the arm is dropped (STOP:
#             the reading of the stop is written into OUT, nothing runs, no try counts), 22 refused, 24 exists,
#             20 other
#   guard     the guard part (base arms only) in seven parallel parts, then the reproduction gate; 26 the gate fails
#             (the VOID reading, which counts the three tries, is written into OUT first; a rerun finishes 26 again
#             and never reads twice), 30 a part did not exit 0, 22 refused, 20 other
#   run       the scoring parts in seven parallel parts (plain g1 and g2 as seeds 0-4 and 5-9, dm with 8 workers by 2
#             threads, x2 as seeds 0-4 and 5-9); a seed whose file exists is skipped, so a killed run resumes; 30 a
#             part did not exit 0
#   read      the reading into experiments/t3_decomp/results/k26_query_map/ (output to /dev/null): 0 a reading (VOID
#             included), 43 INCOMPLETE, 44 a reading is there already (never overwritten), 41 a crash
# A REFUSED call writes $S/<STEP>.refused.exit only (unknown.refused.exit for an unknown step) and never touches a
# step's .exit or .head: the step before has not exited 0 at COMMIT (96), HEAD is not COMMIT (99), a tracked file
# changed (98), no TRIES.md row K26 names the script (95), an input is missing (94: the outside folders and files below
# are not real, TABLES_K24.sha256 untracked, or EXTERNAL_DATA.md rows X013, X014, X015 and X017 do not each name K26),
# an unknown step (97). A call without two arguments writes $HOME/v/scratch/k26/usage.exit (2); when the scratch folder
# cannot be made the call exits 90 and writes nothing. A call that passes every check first overwrites $S/<STEP>.exit
# with the word running, then writes $S/<STEP>.head, runs, and writes its exit code to $S/<STEP>.exit (a TERM, INT or
# HUP signal writes 89), so a step killed part way never leaves an earlier call's 0. $S is
# $HOME/v/scratch/k26/<first 7 of COMMIT>. Logs go to R/logs/ and are never printed; their lines other than the
# script's step lines are counted into R/exit/<log name>.other_lines. A stop's reading is made once (a rerun finds
# OUT/reading.json and finishes with the same code; 43, 44 or 41 when that read is INCOMPLETE, finds a reading, or
# crashes). Nothing is ever removed. 8 BLAS threads a process (16 in the dm set-up), NUMPY_MADVISE_HUGEPAGE=0, TMPDIR
# under $S.
set -u
if [ "$#" -lt 2 ] || [ -z "$1" ] || [ -z "$2" ]; then
  mkdir -p "$HOME/v/scratch/k26" && echo 2 > "$HOME/v/scratch/k26/usage.exit"
  echo "usage: run_k26.sh tests|features|stage0|pregate|guard|run|read COMMIT" >&2
  exit 2
fi
STEP="$1"
COMMIT="$2"
S="$HOME/v/scratch/k26/${COMMIT:0:7}"
mkdir -p "$S/tmp" || exit 90
refused() { echo "$1" > "$S/$STEP.refused.exit"; exit "$1"; }
case "$STEP" in
  tests) PREV="" ;;
  features) PREV="tests" ;;
  stage0) PREV="features" ;;
  pregate) PREV="stage0" ;;
  guard) PREV="pregate" ;;
  run) PREV="guard" ;;
  read) PREV="run" ;;
  *)
    echo "unknown step" >&2
    echo 97 > "$S/unknown.refused.exit"
    exit 97 ;;
esac
if [ -n "$PREV" ]; then
  if [ "$(cat "$S/$PREV.exit" 2>/dev/null)" != 0 ] || [ "$(cat "$S/$PREV.head" 2>/dev/null)" != "$COMMIT" ]; then
    echo "refused: step $PREV has not exited 0 at $COMMIT" >&2
    refused 96
  fi
fi
if [ "$(git rev-parse HEAD 2>/dev/null)" != "$COMMIT" ]; then
  echo "refused: HEAD is not $COMMIT" >&2
  refused 99
fi
if [ -n "$(git status --porcelain --untracked-files=no 2>&1)" ]; then
  echo "refused: the tree is not clean" >&2
  refused 98
fi
grep -q '^| K26 .*`experiments/t3_decomp/k26_query_map.py`' experiments/harness/TRIES.md 2>/dev/null \
  || { echo "refused: no TRIES.md row K26 names the script" >&2; refused 95; }

TABLES=$HOME/v/vec-ext/t3/tables  # sc3 and bulk; fixed, no override
TABLES_TWO=$HOME/v/vec-ext/t3/tables_two_study  # `two`; fixed, no override
PROCESSED=$HOME/v/vec-ext/t3/processed  # the per-gene unit files; fixed, no override
REFSTATS=$HOME/v/vec-ext/t3/reference/refstats.npz  # the units' reference statistics; fixed, no override
RUN_DIR=$HOME/v/virtual-embryo-challenge/runs/t3_decomp/gata4_20261028T122425Z  # read only, as row K25; no override
for d in "$TABLES" "$TABLES_TWO" "$PROCESSED"; do
  { [ -d "$d" ] && [ ! -L "$d" ]; } || { echo "refused: $d is not a real directory" >&2; refused 94; }
done
{ [ -f "$REFSTATS" ] && [ ! -L "$REFSTATS" ]; } || { echo "refused: $REFSTATS is not a plain file" >&2; refused 94; }
GIT_OPTIONAL_LOCKS=0 git ls-files --error-unmatch experiments/t3_ext/TABLES_K24.sha256 > /dev/null 2>&1 \
  || { echo "refused: TABLES_K24.sha256 is not tracked" >&2; refused 94; }
for x in X013 X014 X015 X017; do
  grep "^| $x " EXTERNAL_DATA.md 2>/dev/null | grep -qF K26 \
    || { echo "refused: EXTERNAL_DATA.md row $x does not name K26" >&2; refused 94; }
done

# every check passed: from here each path writes $S/$STEP.exit
finish() { trap - TERM INT HUP; echo "$1" > "$S/$STEP.exit"; exit "$1"; }
echo running > "$S/$STEP.exit"
trap 'finish 89' TERM INT HUP
git rev-parse HEAD > "$S/$STEP.head"

export NUMPY_MADVISE_HUGEPAGE=0 TMPDIR="$S/tmp" UV_PROJECT_ENVIRONMENT="$PWD/.venv" VEC_T3_EXT_TABLES="$TABLES" \
  VEC_T3_EXT_TABLES_TWO="$TABLES_TWO" VEC_T3_EXT_PROCESSED="$PROCESSED" VEC_T3_EXT_REFSTATS="$REFSTATS"
PY="uv run --frozen --no-sync python -I -B"
threads() { export OMP_NUM_THREADS=$1 OPENBLAS_NUM_THREADS=$1 MKL_NUM_THREADS=$1; }
threads 8
SCRIPT=experiments/t3_decomp/k26_query_map.py
DATA="$PWD/data/raw"  # this folder's own data, never the installed package's
R=runs/k26_${COMMIT:0:7}
OUT=experiments/t3_decomp/results/k26_query_map
OK_LINE='^(K26( guard)? (plain|dm|x2) [a-z0-9]+ seed [0-9]+: (scored|present, skipped)|wrote [^ ]+|\[rep[12] seed [0-9]+\] K18 depth-matched: [0-9]+ arms and the check, [0-9]+s|self-check passed .*|features: (written|present)|stage 0: (passes|fails)( \(read again from stage0\.json\))?|pre-gate: .*|reproduction gate: (passes|FAILS))$'
PARTS=("plain_g1_0-4" "plain_g1_5-9" "plain_g2_0-4" "plain_g2_5-9" "dm" "x2_x2_0-4" "x2_x2_5-9")
mkdir -p "$R/exit" "$R/logs"

other_lines() {  # $1 = a log name under R/logs (no .log): its lines other than the step lines, counted, never printed
  local n
  n=$(grep -cvE "$OK_LINE" "$R/logs/$1.log" 2>/dev/null) || true
  echo "${n:-unreadable}" > "$R/exit/$1.other_lines"
}

parts() {  # $1 = guard|score: the seven parts at once; returns 1 when a part did not exit 0
  local pids=() part p bad=0
  for part in "${PARTS[@]}"; do
    (
      if [ "$part" = dm ]; then
        threads 16
        $PY "$SCRIPT" run --part "$1" --world dm --workers 8 --threads 2 --out "$R" --data "$DATA" \
          --run-dir "$RUN_DIR" --commit "$COMMIT" > "$R/logs/$1_$part.log" 2>&1
      else
        IFS=_ read -r world task seeds <<< "$part"
        if [ "$world" = plain ]; then
          $PY "$SCRIPT" run --part "$1" --world plain --task "$task" --seeds "$seeds" --out "$R" --data "$DATA" \
            --run-dir "$RUN_DIR" --commit "$COMMIT" > "$R/logs/$1_$part.log" 2>&1
        else
          $PY "$SCRIPT" run --part "$1" --world x2 --seeds "$seeds" --out "$R" --data "$DATA" \
            --run-dir "$RUN_DIR" --commit "$COMMIT" > "$R/logs/$1_$part.log" 2>&1
        fi
      fi
      echo $? > "$R/exit/$1_$part"
    ) &
    pids+=($!)
  done
  for p in "${pids[@]}"; do wait "$p" || true; done
  for part in "${PARTS[@]}"; do
    [ "$(cat "$R/exit/$1_$part" 2>/dev/null)" = 0 ] || bad=1
    other_lines "$1_$part"
  done
  return $bad
}

stop_read() {  # $1 = the step's finish (26, 25 or 27): the reading of a stop or a VOID gate, made once, then finish $1
  local rc=0
  if [ -e "$OUT/reading.json" ]; then finish "$1"; fi
  $PY "$SCRIPT" read --runs "$R" --out "$OUT" --commit "$COMMIT" > /dev/null 2>&1 || rc=$?
  echo $rc > "$S/stop_read.exit"
  if [ $rc = 3 ]; then finish 43; fi
  if [ $rc = 4 ]; then finish 44; fi
  if [ $rc != 0 ]; then finish 41; fi
  finish "$1"
}

rc=0
case "$STEP" in
  tests)
    $PY -m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp="$S/tmp/pytest" \
      tests/test_k26_query_map.py tests/test_k26_context_map.py tests/test_k26_launcher.py \
      tests/test_k24_wider_tables.py tests/test_build_k24_form.py tests/test_k23_outside_exponents.py \
      tests/test_k22_outside_forms.py tests/test_t3_outside_forms.py tests/test_k21_factor_forms.py \
      tests/test_k18_built_files.py tests/test_t3_ext_pipeline.py > "$S/tests.log" 2>&1 || rc=$?
    if [ $rc != 0 ]; then finish 10; fi
    finish 0 ;;
  features)
    $PY "$SCRIPT" features --out "$R" --data "$DATA" --commit "$COMMIT" > "$R/logs/features.log" 2>&1 || rc=$?
    other_lines features
    case $rc in 0) finish 0 ;; 2) finish 22 ;; 4) finish 24 ;; *) finish 20 ;; esac ;;
  stage0)
    $PY "$SCRIPT" stage0 --out "$R" --commit "$COMMIT" > "$R/logs/stage0.log" 2>&1 || rc=$?
    other_lines stage0
    case $rc in 0) finish 0 ;; 7) stop_read 27 ;; 2) finish 22 ;; *) finish 20 ;; esac ;;
  pregate)
    $PY "$SCRIPT" pregate --out "$R" --data "$DATA" --run-dir "$RUN_DIR" --commit "$COMMIT" \
      > "$R/logs/pregate.log" 2>&1 || rc=$?
    other_lines pregate
    case $rc in 0) finish 0 ;; 5) stop_read 25 ;; 2) finish 22 ;; 4) finish 24 ;; *) finish 20 ;; esac ;;
  guard)
    if [ ! -e "$R/guard.json" ]; then
      parts guard || finish 30
      $PY "$SCRIPT" guard --runs "$R" --commit "$COMMIT" > "$R/logs/guard.log" 2>&1 || rc=$?
      other_lines guard
      case $rc in 0) finish 0 ;; 6) stop_read 26 ;; 2) finish 22 ;; 3) finish 30 ;; *) finish 20 ;; esac
    fi
    gpass=$($PY -c 'import json, sys; g = json.load(open(sys.argv[1])); print("pass" if g.get("passes") is True else "fail")' \
      "$R/guard.json" 2> /dev/null) || gpass=unreadable
    if [ "$gpass" = pass ]; then finish 0; fi
    if [ "$gpass" = fail ]; then stop_read 26; fi
    finish 20 ;;
  run)
    parts score || finish 30
    finish 0 ;;
  read)
    $PY "$SCRIPT" read --runs "$R" --out "$OUT" --commit "$COMMIT" > /dev/null 2>&1 || rc=$?
    case $rc in 0) finish 0 ;; 3) finish 43 ;; 4) finish 44 ;; *) finish 41 ;; esac ;;
esac
finish 20
