#!/usr/bin/env bash
# Row K31 (experiments/harness/TRIES.md; experiments/t3_decomp/k31_cube_read.py) with row R72's G0
# (experiments/t3_decomp/cube_tables.py g0): the sign cube's target-free gate and its counted read, on the DGX campaign
# checkout where row K18 ran (~/v/vec-campaign-c, never the reference checkout), all steps in one call from the
# checkout's root. Run it detached (setsid nohup ... &) and wait on $S/done.exit with a background poller; never a
# foreground sleep.
#
#   bash experiments/t3_decomp/run_k31.sh COMMIT RUN
#
# COMMIT: the full hash of the run commit, row K31's commit (rows R72 and K31 are in TRIES.md there). RUN names the
# call's folders: the scratch folder S = ~/v/scratch/k31/RUN (exit records, logs, TMPDIR) and the run folder
# R = runs/k31_RUN in the checkout (the identity build, g0.json, dryrun.json, task files); letters, digits, '.', '_'
# and '-' only, not starting with '.' (else 90). A RUN made for another commit is refused (91); the same RUN at the
# same commit resumes: the parts skip the seeds already written, and a g0.json already there is read again, never
# remade. A call without both arguments exits 2. The gata4 run folder is fixed, no override: the reference checkout's
# runs/t3_decomp/gata4_20261028T122425Z, read only.
# Steps, in order; the call stops at the first that does not exit 0. Each writes S/<step>.exit (first the word running,
# then its code) and S/<step>.head (COMMIT), and runs only when the step before it wrote 0 at COMMIT (else 96).
# S/done.exit holds running, then the call's code. Once a reading is in experiments/t3_decomp/results/k31_cube_read/,
# no part is started (44), whatever RUN.
#   preflight  HEAD is COMMIT (97); the tracked tree is clean (98); TRIES.md holds exactly one row R72 whose first cell
#              opens "R72 (the captain's order", says "ruling, no run" and holds `--captain-ruling R72` and each of
#              the seven `--form c2-s2 --vertex <v>` as one whole backtick-delimited string (95); exactly one row K31
#              naming the read script and recording one cube tables manifest sha256, equal to the committed
#              MANIFEST.json's (94); K18's two committed plain task files are tracked and the gata4 run folder is a
#              directory (96)
#   tests      the targeted tests, synthetic data only, 8 BLAS threads (10)
#   identity   C2's identity rebuild: build_two_knockout.py --form c2-s2 (C2's own command) built once into R/gate/,
#              then `vec check --board T3:bcatenin`, its stamp (commit COMMIT, clean True) and its sha256 against C2's
#              recorded bytes (experiments/builds/results/t3two_fb01487.SHA256SUMS): 11 the build failed, 12 the
#              builder refused, 13 vec check, 14 no stamp, 15 stamp not this commit and clean, 16 not C2's bytes.
#              Its output goes to /dev/null: a test build is read by its exit code, sha256 and stamp only
#   g0         G0 (a to g), target-free, into R/g0.json and R/g0.txt (its log in S/g0.log): 0 CONTINUE goes on; 5
#              STOP: the STOP reading is made once (no try counts) and the call finishes 25; 2 refused (22); 1 a crash
#              (23); other 20. A g0.json already there is read again: CONTINUE at COMMIT goes on; STOP finishes 25 (the
#              reading only when none is there); of another commit 21; unreadable 20
#   dryrun     row K31's guard dry run on K21's committed reproduction (J7; committed results only, not a try), into
#              R/dryrun.json: 0 holds; 6 fails (26); other 20
#   run        44 first when a reading is there (no arm is scored again). THE TRIES COUNT FROM THE START OF THIS STEP.
#              The script's `start` writes R/run_started.json once, before any part starts (32 when it refuses; its
#              log in S/start.log); then four parallel parts, plain g1 and g2 as seeds 0-4 and 5-9, 8 BLAS threads a
#              process; each part's code in R/exit/<part>, its log in R/logs/<part>.log (never printed); log lines
#              other than the expected kinds are counted, never printed (30 if a part did not exit 0)
#   read       the fixed reading into experiments/t3_decomp/results/k31_cube_read/ (output to /dev/null): 0 a reading
#              (whatever it says, VOID included), 43 INCOMPLETE, 44 a reading is there already (never overwritten),
#              42 refused, 41 a crash (masked traceback)
# 99: the call stopped unexpectedly (a failing command, a TERM, INT or HUP signal, or an exit that was not a finish).
# NUMPY_MADVISE_HUGEPAGE=0; TMPDIR under S; no process above 16 BLAS threads. Nothing is deleted.
set -u
if [ "$#" -ne 2 ] || [ -z "$1" ] || [ -z "$2" ]; then
  echo "usage: run_k31.sh COMMIT RUN" >&2
  exit 2
fi
COMMIT="$1"
RUN="$2"
if ! [[ "$RUN" =~ ^[A-Za-z0-9_-][A-Za-z0-9._-]*$ ]]; then
  echo "refused: RUN holds letters, digits, '.', '_' and '-' only, not starting with '.'" >&2
  exit 90
fi
S="$HOME/v/scratch/k31/$RUN"
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
names() {  # cell string: whether the cell holds the string as one whole backtick-delimited string
  local rest=$1 tok
  while [[ $rest == *'`'*'`'* ]]; do
    rest=${rest#*'`'}
    tok=${rest%%'`'*}
    rest=${rest#*'`'}
    [ "$tok" = "$2" ] && return 0
  done
  return 1
}
echo running > "$S/done.exit"

export NUMPY_MADVISE_HUGEPAGE=0 TMPDIR="$S/tmp" UV_PROJECT_ENVIRONMENT="$PWD/.venv"
PY="uv run --frozen --no-sync python -I -B"
threads() { export OMP_NUM_THREADS=$1 OPENBLAS_NUM_THREADS=$1 MKL_NUM_THREADS=$1; }
threads 8
SCRIPT=experiments/t3_decomp/k31_cube_read.py
G0_SCRIPT=experiments/t3_decomp/cube_tables.py
MANIFEST=experiments/t3_decomp/results/t3_cube/MANIFEST.json
R="runs/k31_$RUN"
OUT=experiments/t3_decomp/results/k31_cube_read
RUN_DIR=$HOME/v/virtual-embryo-challenge/runs/t3_decomp/gata4_20261028T122425Z  # read only; fixed, no override
DEC=experiments/t3_decomp/results/gata4/gata4_decisions.json
WT=data/raw/T2/E8.75.h5ad
BOARD=T3:bcatenin
C2_SHA=47111ac1d07e73c0c26e990a474dff0de94ae22d6a4474a569bcb674ae2b85a7  # c2-s2, t3two_fb01487.SHA256SUMS
K18_TASKS=experiments/t3_decomp/results/k18_built_files/tasks
VERTICES="ppn pnp pnn npp npn nnp nnn"
G0_CHECK='import json, sys; g = json.load(open(sys.argv[1])); p = g.get("provenance") or {}; print("foreign" if p.get("commit") != sys.argv[2] or p.get("tree_clean") is not True else {"CONTINUE": "go", "STOP": "stop"}.get(g.get("outcome"), "unreadable"))'
OK_LINE='^K31 plain g[12] seed [0-9]+: (scored|present, skipped)$|^wrote [^ ]+\.json'
PART_SET=("plain_g1_0-4" "plain_g2_0-4" "plain_g1_5-9" "plain_g2_5-9")

reading_there() { [ -e "$OUT/reading.json" ] || [ -e "$OUT/reading.txt" ]; }

# stop_read: G0 said STOP: the STOP reading once (no try counts), then 25
stop_read() {
  local rc=0
  if ! reading_there; then
    $PY "$SCRIPT" read --runs "$R" --out "$OUT" > /dev/null 2>&1 || rc=$?
    echo "$rc" > "$S/stop_read.exit"
    if [ $rc != 0 ]; then finish 41; fi
  fi
  finish 25
}

begin preflight ""
[ "$(git rev-parse HEAD 2> /dev/null)" = "$COMMIT" ] || { echo "refused: HEAD is not $COMMIT" >&2; finish 97; }
[ -z "$(git status --porcelain --untracked-files=no 2>&1)" ] || { echo "refused: the tree is not clean" >&2; finish 98; }
[ "$(grep -c '^| R72 ' experiments/harness/TRIES.md 2> /dev/null)" = 1 ] \
  || { echo "refused: TRIES.md does not hold exactly one row R72" >&2; finish 95; }
first=$(grep '^| R72 ' experiments/harness/TRIES.md | cut -d'|' -f2)
[[ $first == " R72 (the captain's order"* ]] || { echo "refused: row R72 is not the captain's order" >&2; finish 95; }
[[ $first == *"ruling, no run"* ]] || { echo "refused: row R72 does not say: ruling, no run" >&2; finish 95; }
names "$first" "--captain-ruling R72" || { echo "refused: row R72 does not name --captain-ruling R72" >&2; finish 95; }
for v in $VERTICES; do
  names "$first" "--form c2-s2 --vertex $v" || { echo "refused: row R72 does not name vertex $v" >&2; finish 95; }
done
[ "$(grep -c '^| K31 .*`experiments/t3_decomp/k31_cube_read.py`' experiments/harness/TRIES.md 2> /dev/null)" = 1 ] \
  || { echo "refused: TRIES.md does not hold exactly one row K31 naming the script" >&2; finish 94; }
k31_sha=$($PY "$SCRIPT" row 2> /dev/null) || k31_sha=""
man_sha=$(sha256sum "$MANIFEST" 2> /dev/null | cut -d' ' -f1) || man_sha=""
[ -n "$k31_sha" ] && [ "$k31_sha" = "$man_sha" ] \
  || { echo "refused: row K31's manifest sha256 is not the committed manifest's" >&2; finish 94; }
GIT_OPTIONAL_LOCKS=0 git ls-files --error-unmatch "$MANIFEST" > /dev/null 2>&1 \
  || { echo "refused: the manifest is not tracked" >&2; finish 94; }
for f in plain_g1.json plain_g2.json; do
  GIT_OPTIONAL_LOCKS=0 git ls-files --error-unmatch "$K18_TASKS/$f" > /dev/null 2>&1 \
    || { echo "refused: $K18_TASKS/$f is not tracked" >&2; finish 96; }
done
[ -d "$RUN_DIR" ] || { echo "refused: the gata4 run folder $RUN_DIR is not a directory" >&2; finish 96; }
passed

begin tests preflight
rc=0
$PY -m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp="$S/tmp/pytest" tests/test_sign_cube.py \
  tests/test_cube_tables.py tests/test_k31_cube_read.py tests/test_k31_launcher.py tests/test_k18_built_files.py \
  > "$S/tests.log" 2>&1 || rc=$?
if [ $rc != 0 ]; then finish 10; fi
passed

begin identity tests
mkdir -p "$R/gate" "$R/exit" "$R/logs"
if [ ! -e "$R/gate/c2-s2.h5ad" ]; then
  threads 16
  rc=0
  $PY experiments/t3_decomp/build_two_knockout.py build "$DEC" --form c2-s2 --wt "$WT" --run-dir "$RUN_DIR" \
    --board "$BOARD" --seed 1000001 --out "$R/gate/c2-s2.h5ad" > /dev/null 2>&1 || rc=$?
  threads 8
  if [ $rc = 2 ]; then finish 12; fi
  if [ $rc != 0 ]; then finish 11; fi
fi
uv run --frozen --no-sync vec check "$R/gate/c2-s2.h5ad" --board "$BOARD" > /dev/null 2>&1 || finish 13
st=$($PY -c 'import json,sys; s=json.load(open(sys.argv[1])); print(s["commit"], s["clean"])' \
  "$R/gate/c2-s2.h5ad.stamp.json" 2> /dev/null) || finish 14
[ "$st" = "$COMMIT True" ] || { echo "refused: the identity build's stamp is not $COMMIT clean" >&2; finish 15; }
got=$(sha256sum "$R/gate/c2-s2.h5ad" | cut -d' ' -f1) || got=""
[ "$got" = "$C2_SHA" ] || { echo "stopped: the identity build is not C2's recorded bytes" >&2; finish 16; }
passed

begin g0 identity
if [ ! -e "$R/g0.json" ]; then
  threads 16
  rc=0
  $PY "$G0_SCRIPT" g0 --data data/raw --run-dir "$RUN_DIR" --decision "$DEC" --manifest "$MANIFEST" \
    --manifest-sha256 "$k31_sha" --board "$BOARD" --out "$R" --threads 16 > "$S/g0.log" 2>&1 || rc=$?
  threads 8
  if [ $rc = 5 ]; then stop_read; fi
  if [ $rc = 2 ]; then finish 22; fi
  if [ $rc = 1 ]; then finish 23; fi
  if [ $rc != 0 ]; then finish 20; fi
fi
got=$($PY -c "$G0_CHECK" "$R/g0.json" "$COMMIT" 2> /dev/null) || got=unreadable
if [ "$got" = stop ]; then stop_read; fi
if [ "$got" = foreign ]; then finish 21; fi
[ "$got" = go ] || finish 20
passed

begin dryrun g0
rc=0
$PY "$SCRIPT" dryrun --out "$R/dryrun.json" > "$S/dryrun.log" 2>&1 || rc=$?
if [ $rc = 6 ]; then finish 26; fi
if [ $rc != 0 ]; then finish 20; fi
passed

begin run dryrun
if reading_there; then echo "refused: a reading is there already; no arm is scored again" >&2; finish 44; fi
rc=0
$PY "$SCRIPT" start --out "$R" > "$S/start.log" 2>&1 || rc=$?
if [ $rc != 0 ]; then finish 32; fi
pids=()
for part in "${PART_SET[@]}"; do
  (
    IFS=_ read -r _ task seeds <<< "$part"
    $PY "$SCRIPT" run --task "$task" --seeds "$seeds" --out "$R" --data data/raw --run-dir "$RUN_DIR" \
      --decision "$DEC" --manifest "$MANIFEST" > "$R/logs/$part.log" 2>&1
    echo $? > "$R/exit/$part"
  ) &
  pids+=($!)
done
for p in "${pids[@]}"; do wait "$p" || true; done
bad=0
for part in "${PART_SET[@]}"; do
  [ "$(cat "$R/exit/$part" 2> /dev/null)" = 0 ] || bad=1
  n=$(grep -cvE "$OK_LINE" "$R/logs/$part.log" 2> /dev/null) || true
  echo "${n:-unreadable}" > "$R/exit/$part.other_lines"
done
if [ $bad != 0 ]; then finish 30; fi
passed

begin read run
if reading_there; then echo "refused: a reading is there already; it is never overwritten" >&2; finish 44; fi
rc=0
$PY "$SCRIPT" read --runs "$R" --out "$OUT" > /dev/null 2>&1 || rc=$?
if [ $rc = 3 ]; then finish 43; fi
if [ $rc = 4 ]; then finish 44; fi
if [ $rc = 2 ]; then finish 42; fi
if [ $rc != 0 ]; then finish 41; fi
finish 0
