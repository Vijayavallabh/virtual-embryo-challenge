#!/usr/bin/env bash
# Row K24 (experiments/harness/TRIES.md; experiments/t3_decomp/k24_wider_tables.py) on a DGX campaign checkout
# (~/v/vec-campaign-c, never the reference checkout), the machine rows K18 and K21 to K23 ran on. Run from the
# checkout's root, detached (setsid nohup ... &), and wait on $S/all.exit with a background poller; never a foreground
# sleep.
#
#   bash experiments/t3_decomp/run_k24.sh COMMIT
#
# COMMIT: the full hash of the row's commit (the row is in TRIES.md at that commit; TABLES_K24.sha256, the filled pins
# in k24_wider_tables.py and K23's committed pregate.json and reading.json are there). The two tables folders are fixed
# below, ~/v/vec-ext/t3/tables (sc3 and bulk) and ~/v/vec-ext/t3/tables_two_study (`two`), and reach the script only
# through VEC_T3_EXT_TABLES and VEC_T3_EXT_TABLES_TWO. K23's run folder runs/k23_<first 7 of K23's commit> must be in
# this checkout (guard 2 reads its per-seed x_e2 at the sha256 K23's reading records). S = ~/v/scratch/k24_<first 7>,
# the run folder R = runs/k24_<first 7> (resumable: a rerun at the same commit skips the seeds already written).
# Steps, in order, stopping at the first failure:
#   1. HEAD is COMMIT (exit 97), the tracked tree is clean (98), the row naming the script is in TRIES.md (95),
#      TABLES_K24.sha256 and K23's committed pregate.json and reading.json are tracked, both tables folders are real
#      directories (no symlink), and EXTERNAL_DATA.md holds a row naming E-MTAB-7322 and one naming GSE299068 (96);
#   2. the targeted tests, synthetic data only, 8 BLAS threads (exit 10);
#   3. the pre-gate, target-free, before any scoring: R/pregate.json and $S/pregate.log (no score). Exit 22 when the
#      script refuses (its exit 2: pins, listing, folders, bytes or records disagree, a missing or unpinned table,
#      frozen files off their pins, `two` off K23's figures, K23's per-seed files missing; nothing counted), 25 when it
#      drops both tables (its exit 5: nothing runs, nothing counted), 20 any other failure;
#   4. the GUARD PART: barred_form and x_e2 only, both worlds, ten seeds, in five parallel parts as step 6 (into
#      R/guard/; exit 30 if a part did not exit 0), then the guard check (R/guard.json): its exit 6 (a guard failed)
#      reads the run into OUT (a VOID reading that counts the failed guard arms only) and finishes 26; its exit 2
#      finishes 22, 3 finishes 30, anything else 20;
#   5. guard.json's key "passes", read by a python -I -B one-liner: not true finishes 26 again with no scoring part
#      started (a rerun after a guard failure; the VOID read is made only when OUT holds no reading.json, so nothing is
#      counted twice), unreadable finishes 20; the scoring parts also refuse without a passing guard.json of this
#      commit and machine;
#   6. the five scoring parts in parallel: plain g1 and g2, each as seeds 0-4 and 5-9 (8 BLAS threads a process), and
#      the depth-matched world (8 workers by 2 threads, 16 BLAS threads in its set-up, as rows U3d, K18 and K21 to
#      K23); each part's exit code in R/exit/<part>, its log in R/logs/<part>.log (never printed); the log lines other
#      than the expected kinds are counted, never printed (exit 30 if a part did not exit 0);
#   7. the read into experiments/t3_decomp/results/k24_wider_tables/ (output to /dev/null): exit 0 on a reading
#      (whatever it says, VOID included), 43 INCOMPLETE (a file missing; rerun this launcher at the same commit), 44 a
#      reading is there already (never overwritten), 41 a crash in the read (masked traceback; fix the reading code
#      only, amend the row, read again into the emptied folder).
# all.exit: 0, 10, 20, 22, 25, 26, 30, 41, 43, 44, 95, 96, 97, 98 as above; 99 the launcher stopped unexpectedly (a
# failing command, a TERM, INT or HUP signal, or an exit that was not one of the finishes). NUMPY_MADVISE_HUGEPAGE=0;
# TMPDIR under $S.
set -u
COMMIT=${1:?usage: run_k24.sh COMMIT}
S=~/v/scratch/k24_${COMMIT:0:7}
mkdir -p "$S/tmp" || exit 99
finish() { trap - ERR TERM INT HUP EXIT; echo "$1" > "$S/all.exit"; exit "$1"; }
abort() { trap - ERR TERM INT HUP EXIT; echo 99 > "$S/all.exit"; exit 99; }
trap abort ERR TERM INT HUP EXIT

TABLES=$HOME/v/vec-ext/t3/tables  # combine.py's --out (sc3, bulk); fixed, no override
TABLES_TWO=$HOME/v/vec-ext/t3/tables_two_study  # the two-study tables K22 and K23 read; fixed, no override
export NUMPY_MADVISE_HUGEPAGE=0 TMPDIR="$S/tmp" UV_PROJECT_ENVIRONMENT="$PWD/.venv" VEC_T3_EXT_TABLES="$TABLES" \
  VEC_T3_EXT_TABLES_TWO="$TABLES_TWO"
PY="uv run --frozen --no-sync python -I -B"
threads() { export OMP_NUM_THREADS=$1 OPENBLAS_NUM_THREADS=$1 MKL_NUM_THREADS=$1; }
SCRIPT=experiments/t3_decomp/k24_wider_tables.py
R=runs/k24_${COMMIT:0:7}
RUN_DIR=$HOME/v/virtual-embryo-challenge/runs/t3_decomp/gata4_20261028T122425Z  # read only; fixed, no override
OUT=experiments/t3_decomp/results/k24_wider_tables

[ "$(git rev-parse HEAD)" = "$COMMIT" ] || { echo "refused: HEAD is not $COMMIT" >&2; finish 97; }
[ -z "$(git status --porcelain --untracked-files=no)" ] || { echo "refused: the tree is not clean" >&2; finish 98; }
grep -q '^| .*`experiments/t3_decomp/k24_wider_tables.py`' experiments/harness/TRIES.md \
  || { echo "refused: no TRIES.md row names the script" >&2; finish 95; }
for f in experiments/t3_ext/TABLES_K24.sha256 experiments/t3_decomp/results/k23_outside_exponents/pregate.json \
    experiments/t3_decomp/results/k23_outside_exponents/reading.json; do
  GIT_OPTIONAL_LOCKS=0 git ls-files --error-unmatch "$f" > /dev/null 2>&1 \
    || { echo "refused: $f is not tracked" >&2; finish 96; }
done
for d in "$TABLES" "$TABLES_TWO"; do
  { [ -d "$d" ] && [ ! -L "$d" ]; } || { echo "refused: $d is not a real directory" >&2; finish 96; }
done
for acc in E-MTAB-7322 GSE299068; do
  grep '^| X[0-9]' EXTERNAL_DATA.md 2> /dev/null | grep -qF "$acc" \
    || { echo "refused: EXTERNAL_DATA.md holds no row naming $acc" >&2; finish 96; }
done
git rev-parse HEAD > "$S/head"
mkdir -p "$R/exit" "$R/logs"

threads 8
rc=0
$PY -m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp="$S/tmp/pytest" \
  tests/test_k24_wider_tables.py tests/test_build_k24_form.py tests/test_k23_outside_exponents.py \
  tests/test_k22_outside_forms.py tests/test_t3_outside_forms.py tests/test_k21_factor_forms.py \
  tests/test_k18_built_files.py > "$S/tests.log" 2>&1 || rc=$?
echo $rc > "$S/tests.exit"
if [ $rc != 0 ]; then finish 10; fi

if [ ! -e "$R/pregate.json" ]; then
  rc=0
  $PY "$SCRIPT" pregate --out "$R" --run-dir "$RUN_DIR" > "$S/pregate.log" 2>&1 || rc=$?
  echo $rc > "$S/pregate.exit"
  if [ $rc = 2 ]; then finish 22; fi
  if [ $rc = 5 ]; then finish 25; fi
  if [ $rc != 0 ]; then finish 20; fi
fi

# a failed guard part: the VOID reading (counts the failed guard arms only), then finish 26
void_read() {
  rc=0
  $PY "$SCRIPT" read --runs "$R" --out "$OUT" > /dev/null 2>&1 || rc=$?
  echo $rc > "$S/read.exit"
  if [ $rc = 3 ]; then finish 43; fi
  if [ $rc = 4 ]; then finish 44; fi
  if [ $rc != 0 ]; then finish 41; fi
  finish 26
}

# the guard part: barred_form and x_e2, both worlds, ten seeds, before any try arm
if [ ! -e "$R/guard.json" ]; then
  GSET=("plain_g1_0-4" "plain_g2_0-4" "plain_g1_5-9" "plain_g2_5-9" "dm")
  gpids=()
  for gpart in "${GSET[@]}"; do
    (
      if [ "$gpart" = dm ]; then
        threads 16
        $PY "$SCRIPT" run --part guard --world dm --workers 8 --threads 2 --out "$R" --run-dir "$RUN_DIR" \
          > "$R/logs/guard_$gpart.log" 2>&1
      else
        threads 8
        IFS=_ read -r _ task seeds <<< "$gpart"
        $PY "$SCRIPT" run --part guard --world plain --task "$task" --seeds "$seeds" --out "$R" --run-dir "$RUN_DIR" \
          > "$R/logs/guard_$gpart.log" 2>&1
      fi
      echo $? > "$R/exit/guard_$gpart"
    ) &
    gpids+=($!)
  done
  for p in "${gpids[@]}"; do wait "$p" || true; done
  gbad=0
  for gpart in "${GSET[@]}"; do
    [ "$(cat "$R/exit/guard_$gpart" 2>/dev/null)" = 0 ] || gbad=1
    n=$(grep -cvE '^K24 guard plain g[12] seed [0-9]+: (scored|present, skipped)$|^wrote [^ ]+\.json|^\[rep[12] seed [0-9]+\] K18 depth-matched: [0-9]+ arms and the check, [0-9]+s$|^self-check passed ' \
        "$R/logs/guard_$gpart.log" 2>/dev/null) || true
    echo "${n:-unreadable}" > "$R/exit/guard_$gpart.other_lines"
  done
  if [ $gbad != 0 ]; then finish 30; fi
  rc=0
  threads 8
  $PY "$SCRIPT" guard --runs "$R" > "$S/guard.log" 2>&1 || rc=$?
  echo $rc > "$S/guard.exit"
  if [ $rc = 6 ]; then void_read; fi
  if [ $rc = 2 ]; then finish 22; fi
  if [ $rc = 3 ]; then finish 30; fi
  if [ $rc != 0 ]; then finish 20; fi
fi
# step 5: a guard.json that does not pass (a rerun after a guard failure) finishes 26 again; no scoring part starts
gpass=$($PY -c 'import json, sys; g = json.load(open(sys.argv[1])); print("pass" if g.get("passes") is True else "fail")' \
  "$R/guard.json" 2> /dev/null) || gpass=unreadable
if [ "$gpass" = fail ]; then
  if [ -e "$OUT/reading.json" ]; then finish 26; fi
  void_read
fi
if [ "$gpass" != pass ]; then finish 20; fi

PARTS=("plain_g1_0-4" "plain_g2_0-4" "plain_g1_5-9" "plain_g2_5-9" "dm")
pids=()
for part in "${PARTS[@]}"; do
  (
    if [ "$part" = dm ]; then
      threads 16
      $PY "$SCRIPT" run --world dm --workers 8 --threads 2 --out "$R" --run-dir "$RUN_DIR" \
        > "$R/logs/$part.log" 2>&1
    else
      threads 8
      IFS=_ read -r _ task seeds <<< "$part"
      $PY "$SCRIPT" run --world plain --task "$task" --seeds "$seeds" --out "$R" --run-dir "$RUN_DIR" \
        > "$R/logs/$part.log" 2>&1
    fi
    echo $? > "$R/exit/$part"
  ) &
  pids+=($!)
done
for p in "${pids[@]}"; do wait "$p" || true; done
bad=0
for part in "${PARTS[@]}"; do
  [ "$(cat "$R/exit/$part" 2>/dev/null)" = 0 ] || bad=1
  # the expected kinds: progress and write lines of the script, K18's per-row timing line and its fixed self-check
  # line; anything else (imported code, warnings, the frozen script's set-up lines) is counted, never printed
  n=$(grep -cvE '^K24 plain g[12] seed [0-9]+: (scored|present, skipped)$|^wrote [^ ]+\.json|^\[rep[12] seed [0-9]+\] K18 depth-matched: [0-9]+ arms and the check, [0-9]+s$|^self-check passed ' \
      "$R/logs/$part.log" 2>/dev/null) || true
  echo "${n:-unreadable}" > "$R/exit/$part.other_lines"
done
cat "$R"/exit/*.other_lines | awk '{ if ($1 ~ /^[0-9]+$/) s += $1; else u = 1 } END { print (u ? "unreadable" : s) }' \
  > "$R/exit/other_lines"
if [ $bad != 0 ]; then finish 30; fi

rc=0
$PY "$SCRIPT" read --runs "$R" --out "$OUT" > /dev/null 2>&1 || rc=$?
echo $rc > "$S/read.exit"
if [ $rc = 3 ]; then finish 43; fi
if [ $rc = 4 ]; then finish 44; fi
if [ $rc != 0 ]; then finish 41; fi
rm -rf "$S/tmp"
finish 0
