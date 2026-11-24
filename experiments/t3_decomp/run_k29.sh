#!/usr/bin/env bash
# Row K29 (experiments/harness/TRIES.md; experiments/t3_decomp/k29_cd7_read.py): the practice read that T3's file cD7 on
# C2 is bound to (ruling row R36), on the DGX campaign checkout where row K18 ran (~/v/vec-campaign-c, never the
# reference checkout), all steps in one call from the checkout's root. Run it detached (setsid nohup ... &) and wait on
# $S/done.exit with a background poller; never a foreground sleep.
#
#   bash experiments/t3_decomp/run_k29.sh COMMIT RUN
#
# COMMIT: the full hash of the run commit, row K29's commit (the row naming the script is in TRIES.md there). RUN names
# the call's folders: the scratch folder S = ~/v/scratch/k29/RUN (exit records, the steps' logs, TMPDIR) and the run
# folder R = runs/k29_RUN in the checkout (gate.json, task files, guard.json); letters, digits, '.', '_' and '-' only,
# not starting with '.' (else 90). A RUN made for another commit is refused (91); the same RUN at the same commit
# resumes: the runs skip the seeds already written, and a gate.json or guard.json already there is read again, never
# remade. A call without both arguments exits 2. The gata4 run folder is fixed, no override: the reference checkout's
# runs/t3_decomp/gata4_20261028T122425Z, read only.
# Steps, in order; the call stops at the first that does not exit 0. Each writes S/<step>.exit (first the word running,
# then its code) and S/<step>.head (COMMIT), and runs only when the step before it wrote 0 at COMMIT (else 96).
# S/done.exit holds running, then the call's code. A failed gate or guard makes the VOID reading once (run_k24.sh's
# rule; its code in S/void_read.exit) and ends the call under its own code, no later step started. At or after the
# cutoff, 2026-11-25 06:00 UTC, the script refuses the gate, every part and the guard check (its exit 7): the call then
# makes the reading once (the VOID reading 'INCOMPLETE past its cutoff' for a folder still missing a file) and finishes
# 27. Once a reading is in experiments/t3_decomp/results/k29_cd7_read/, no part is started (44), whatever RUN.
#   preflight  HEAD is COMMIT (97), the tracked tree is clean (98), row K29 naming the script is in TRIES.md (95),
#              K18's three committed task files are tracked and the gata4 run folder is a directory (96)
#   tests      the targeted tests, synthetic data only, 8 BLAS threads (10)
#   gate       the target-free gate, R/gate.json (its counts and ratios in S/gate.log): 0 passes; 5 a check fails:
#              nothing is scored and no try counts, the VOID reading naming the check is made once and the call
#              finishes 25; 1 a crash: INCOMPLETE, no reading (23); 2 refused (22); other 20. A gate.json already there
#              is read again: passing at COMMIT goes on; failing finishes 25 (the VOID reading only when none is
#              there); of another commit 21; unreadable 20
#   guard      the guard part (44 first when a reading is there): oth_f only, both worlds, ten seeds, five parallel
#              parts as in the run step, into R/guard/ (30 if a part did not exit 0); then the guard check,
#              R/guard.json: 0 the guard holds; 6 it
#              fails: the VOID reading is made once (it counts one try, oth_f) and the call finishes 26, no scoring
#              part started; 2 refused (22: K18's files off their pins, files of another commit or machine); 3 a guard
#              file missing (30); other 20. A guard.json already there is read again: passing at COMMIT goes on;
#              failing finishes 26 (the VOID reading only when none is there); of another commit 21; unreadable 20
#   run        44 first when a reading is there (no arm is scored again). THE FIVE TRIES COUNT FROM THE START OF THIS
#              STEP (the script writes R/run_started.json). Five parallel
#              scoring parts, the six arms each: plain g1 and g2, each as seeds 0-4 and 5-9 (8 BLAS threads a
#              process), and the depth-matched world (8 workers by 2 threads, 16 BLAS threads in its set-up), as rows
#              K18 and K21 to K24 ran them; the script refuses without a passing gate.json and guard.json of this
#              commit and machine. Each part's code in R/exit/<part>, its log in R/logs/<part>.log (never printed);
#              log lines other than the expected kinds are counted, never printed (30 if a part did not exit 0)
#   read       the fixed reading into experiments/t3_decomp/results/k29_cd7_read/ (output to /dev/null): 0 a reading
#              (whatever it says, VOID included; at or after 2026-11-25 06:00 UTC a folder still missing a file gives
#              the VOID reading 'INCOMPLETE past its cutoff'), 43 INCOMPLETE, 44 a reading is there already (never
#              overwritten), 42 refused (no row, a dirty tree, another HEAD), 41 a crash (masked traceback)
# 99: the call stopped unexpectedly (a failing command, a TERM, INT or HUP signal, or an exit that was not a finish).
# NUMPY_MADVISE_HUGEPAGE=0; TMPDIR under S; no process above 16 BLAS threads. Nothing is deleted.
set -u
if [ "$#" -ne 2 ] || [ -z "$1" ] || [ -z "$2" ]; then
  echo "usage: run_k29.sh COMMIT RUN" >&2
  exit 2
fi
COMMIT="$1"
RUN="$2"
if ! [[ "$RUN" =~ ^[A-Za-z0-9_-][A-Za-z0-9._-]*$ ]]; then
  echo "refused: RUN holds letters, digits, '.', '_' and '-' only, not starting with '.'" >&2
  exit 90
fi
S="$HOME/v/scratch/k29/$RUN"
if [ -e "$S/commit" ] && [ "$(cat "$S/commit" 2> /dev/null)" != "$COMMIT" ]; then
  echo "refused: $S was made for another commit" >&2
  exit 91
fi
mkdir -p "$S/tmp" || exit 92
echo "$COMMIT" > "$S/commit" || exit 92
STEP=preflight
finish() { trap - ERR TERM INT HUP EXIT; echo "$1" > "$S/$STEP.exit"; echo "$1" > "$S/done.exit"; exit "$1"; }
abort() { trap - ERR TERM INT HUP EXIT; echo 99 > "$S/$STEP.exit"; echo 99 > "$S/done.exit"; exit 99; }
trap abort ERR TERM INT HUP EXIT
begin() {  # $1 the step, $2 the step before it (empty for the first)
  STEP="$1"
  echo running > "$S/$STEP.exit"
  echo "$COMMIT" > "$S/$STEP.head"
  if [ -n "$2" ]; then
    if [ "$(cat "$S/$2.exit" 2> /dev/null)" != 0 ] || [ "$(cat "$S/$2.head" 2> /dev/null)" != "$COMMIT" ]; then
      echo "refused: step $2 has not exited 0 at $COMMIT" >&2
      finish 96
    fi
  fi
}
passed() { echo 0 > "$S/$STEP.exit"; }
echo running > "$S/done.exit"

export NUMPY_MADVISE_HUGEPAGE=0 TMPDIR="$S/tmp" UV_PROJECT_ENVIRONMENT="$PWD/.venv"
PY="uv run --frozen --no-sync python -I -B"
threads() { export OMP_NUM_THREADS=$1 OPENBLAS_NUM_THREADS=$1 MKL_NUM_THREADS=$1; }
threads 8
SCRIPT=experiments/t3_decomp/k29_cd7_read.py
R="runs/k29_$RUN"
OUT=experiments/t3_decomp/results/k29_cd7_read
RUN_DIR=$HOME/v/virtual-embryo-challenge/runs/t3_decomp/gata4_20261028T122425Z  # read only; fixed, no override
K18_TASKS=experiments/t3_decomp/results/k18_built_files/tasks
PASSES_CHECK='import json, sys; g = json.load(open(sys.argv[1])); c = (g.get("provenance") or {}).get("commit"); print("foreign" if c != sys.argv[2] else "pass" if g.get("passes") is True else "fail")'
OK_LINE='^K29 (guard )?plain g[12] seed [0-9]+: (scored|present, skipped)$|^wrote [^ ]+\.json|^\[rep[12] seed [0-9]+\] K18 depth-matched: [0-9]+ arms and the check, [0-9]+s$|^self-check passed '
PART_SET=("plain_g1_0-4" "plain_g2_0-4" "plain_g1_5-9" "plain_g2_5-9" "dm")

# reading_there: a reading is in OUT already, so no part may score an arm again
reading_there() { [ -e "$OUT/reading.json" ] || [ -e "$OUT/reading.txt" ]; }

# late PREFIX: true when a part was refused past the cutoff (its exit 7)
late() {
  local part
  for part in "${PART_SET[@]}"; do
    if [ "$(cat "$R/exit/$1$part" 2> /dev/null)" = 7 ]; then return 0; fi
  done
  return 1
}

# void_read CODE: a failed gate (25) or guard (26), or a refusal past the cutoff (27): the reading once (VOID, or for a
# folder complete before the cutoff the reading it supports), then CODE; no later step starts
void_read() {
  rc=0
  $PY "$SCRIPT" read --commit "$COMMIT" --runs "$R" --out "$OUT" > /dev/null 2>&1 || rc=$?
  echo "$rc" > "$S/void_read.exit"
  if [ $rc = 3 ]; then finish 43; fi
  if [ $rc = 4 ]; then finish 44; fi
  if [ $rc = 2 ]; then finish 42; fi
  if [ $rc != 0 ]; then finish 41; fi
  finish "$1"
}

# again FILE CODE: a gate.json or guard.json already there, read again: passing at COMMIT goes on; failing finishes
# CODE (the VOID reading only when none is there); another commit 21; unreadable 20
again() {
  local got
  got=$($PY -c "$PASSES_CHECK" "$R/$1" "$COMMIT" 2> /dev/null) || got=unreadable
  if [ "$got" = pass ]; then return 0; fi
  if [ "$got" = fail ]; then
    if [ -e "$OUT/reading.json" ] || [ -e "$OUT/reading.txt" ]; then finish "$2"; fi
    void_read "$2"
  fi
  if [ "$got" = foreign ]; then finish 21; fi
  finish 20
}

# parts PREFIX [ARGS...]: the five parallel parts (plain g1 and g2 as seeds 0-4 and 5-9, and dm); each part's code in
# R/exit/PREFIX<part>, its log in R/logs/PREFIX<part>.log, its other lines counted; returns 1 when a part did not exit 0
parts() {
  local prefix="$1"
  shift
  local set=("${PART_SET[@]}")
  local pids=() part task seeds n bad=0
  for part in "${set[@]}"; do
    (
      if [ "$part" = dm ]; then
        threads 16
        $PY "$SCRIPT" run --commit "$COMMIT" "$@" --world dm --workers 8 --threads 2 --out "$R" --run-dir "$RUN_DIR" \
          > "$R/logs/$prefix$part.log" 2>&1
      else
        threads 8
        IFS=_ read -r _ task seeds <<< "$part"
        $PY "$SCRIPT" run --commit "$COMMIT" "$@" --world plain --task "$task" --seeds "$seeds" --out "$R" \
          --run-dir "$RUN_DIR" > "$R/logs/$prefix$part.log" 2>&1
      fi
      echo $? > "$R/exit/$prefix$part"
    ) &
    pids+=($!)
  done
  for p in "${pids[@]}"; do wait "$p" || true; done
  for part in "${set[@]}"; do
    [ "$(cat "$R/exit/$prefix$part" 2> /dev/null)" = 0 ] || bad=1
    n=$(grep -cvE "$OK_LINE" "$R/logs/$prefix$part.log" 2> /dev/null) || true
    echo "${n:-unreadable}" > "$R/exit/$prefix$part.other_lines"
  done
  return $bad
}

begin preflight ""
[ "$(git rev-parse HEAD 2> /dev/null)" = "$COMMIT" ] || { echo "refused: HEAD is not $COMMIT" >&2; finish 97; }
[ -z "$(git status --porcelain --untracked-files=no 2>&1)" ] || { echo "refused: the tree is not clean" >&2; finish 98; }
grep -q '^| K29 .*`experiments/t3_decomp/k29_cd7_read.py`' experiments/harness/TRIES.md 2> /dev/null \
  || { echo "refused: no TRIES.md row K29 names the script" >&2; finish 95; }
for f in plain_g1.json plain_g2.json dm.json; do
  GIT_OPTIONAL_LOCKS=0 git ls-files --error-unmatch "$K18_TASKS/$f" > /dev/null 2>&1 \
    || { echo "refused: $K18_TASKS/$f is not tracked" >&2; finish 96; }
done
[ -d "$RUN_DIR" ] || { echo "refused: the gata4 run folder $RUN_DIR is not a directory" >&2; finish 96; }
passed

begin tests preflight
rc=0
$PY -m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp="$S/tmp/pytest" tests/test_k29_cd7_read.py \
  tests/test_k18_built_files.py > "$S/tests.log" 2>&1 || rc=$?
if [ $rc != 0 ]; then finish 10; fi
passed

begin gate tests
mkdir -p "$R/exit" "$R/logs"
if [ ! -e "$R/gate.json" ]; then
  rc=0
  $PY "$SCRIPT" gate --commit "$COMMIT" --out "$R" --threads 8 > "$S/gate.log" 2>&1 || rc=$?
  if [ $rc = 5 ]; then void_read 25; fi
  if [ $rc = 7 ]; then void_read 27; fi
  if [ $rc = 1 ]; then finish 23; fi
  if [ $rc = 2 ]; then finish 22; fi
  if [ $rc != 0 ]; then finish 20; fi
fi
again gate.json 25
passed

begin guard gate
if [ ! -e "$R/guard.json" ]; then
  if reading_there; then echo "refused: a reading is there already; no arm is scored again" >&2; finish 44; fi
  rc=0
  parts guard_ --part guard || rc=$?
  if [ $rc != 0 ]; then
    if late guard_; then void_read 27; fi
    finish 30
  fi
  rc=0
  $PY "$SCRIPT" guard --commit "$COMMIT" --runs "$R" > "$S/guard.log" 2>&1 || rc=$?
  if [ $rc = 6 ]; then void_read 26; fi
  if [ $rc = 7 ]; then void_read 27; fi
  if [ $rc = 2 ]; then finish 22; fi
  if [ $rc = 3 ]; then finish 30; fi
  if [ $rc != 0 ]; then finish 20; fi
fi
again guard.json 26
passed

begin run guard
if reading_there; then echo "refused: a reading is there already; no arm is scored again" >&2; finish 44; fi
rc=0
parts "" --part score || rc=$?
cat "$R"/exit/*.other_lines | awk '{ if ($1 ~ /^[0-9]+$/) s += $1; else u = 1 } END { print (u ? "unreadable" : s) }' \
  > "$R/exit/other_lines"
if [ $rc != 0 ]; then
  if late ""; then void_read 27; fi
  finish 30
fi
passed

begin read run
if [ -e "$OUT/reading.json" ] || [ -e "$OUT/reading.txt" ]; then
  echo "refused: a reading is there already; it is never overwritten" >&2
  finish 44
fi
rc=0
$PY "$SCRIPT" read --commit "$COMMIT" --runs "$R" --out "$OUT" > /dev/null 2>&1 || rc=$?
if [ $rc = 3 ]; then finish 43; fi
if [ $rc = 4 ]; then finish 44; fi
if [ $rc = 2 ]; then finish 42; fi
if [ $rc != 0 ]; then finish 41; fi
finish 0
