#!/usr/bin/env bash
# Row R36's T3 file cD7 (experiments/harness/TRIES.md; the captain's named exception of 2026-11-24, decision log
# 2026-11-24 (aa) and (bb)): C2 (build_two_knockout.py --form c2-s2, recorded sha256 47111ac1d07e...) rebuilt by its
# own command with one change, `--composition-from allowed-s1 --captain-ruling R36`, which replaces its composition
# block by the decision file's submissions.allowed.s1.composition. Modelled on build_two_knockout.sh. Test-board builds
# on a DGX campaign checkout (never the reference checkout). Run from the checkout's root, detached
# (setsid nohup ... &), and wait on the exit files with a background poller; never a foreground sleep:
#   bash experiments/t3_decomp/t3_cd7_build.sh <COMMIT> <RUN>
# RUN names this run's scratch folder ~/v/scratch/t3cd7/<RUN> and its output folder runs/t3cd7_<commit7>_<RUN>.
# Neither may exist: a run folder is never reused and nothing is ever deleted here, so a second run takes a new name.
# In order, each step first writing the word running to ~/v/scratch/t3cd7/<RUN>/<STEP>.exit, then its exit code; the
# code that ends the run also goes to all.exit (running until then):
#   gate   HEAD is COMMIT (the full hash) and the tree is clean (git status --porcelain prints nothing); TRIES.md holds
#          exactly one row R36 whose first cell names `--composition-from allowed-s1`, `--form c2-s2` and
#          `--captain-ruling R36`, and exactly one row R22 whose first cell names `--gata4-exponent 2.0` and
#          `--captain-ruling R22` (the last identity gate's ruling); the output folder does not exist; then it writes
#          the head file
#   tests  the builder's and this launcher's targeted tests (synthetic data only), 8 BLAS threads
#   ids1 ids2 c1s2 c2s2 c1g2   the identity gates, each command built once into gate/, as a build step below, and its
#          sha256 the recorded one: identity-s1 and identity-s2 (the barred s1 and s2, builds_dc7b6ae.SHA256SUMS),
#          c1-s2 and c2-s2 (C1 and C2, t3two_fb01487.SHA256SUMS), and c1-s2 with --gata4-exponent 2.0
#          --captain-ruling R22 (T3_c1_s2_g2.h5ad, t3g2_a3c511c.SHA256SUMS)
#   g1     G1, before the file is built: `build_two_knockout.py g1`, the composition routine on WT E8.75 finds each of
#          the five scenario states at 0.3835131214774091 of its wild-type share before renormalising (within a
#          relative 1e-12); its five ratios go to g1.out
#   a b    the file T3_c2_s2_cd7.h5ad built twice, into folders a and b
#   pair   the two files byte for byte (one sha256)
#   g2     G2, before upload: `build_two_knockout.py g2`, the gate bar 0.01 + 0.99 cos^2(dp, 1) as row N4 computes it,
#          from each file's own build inputs, on each of N4's three gene sets: on every set the file's bar is at most
#          C2's plus 0.01; its six bars (C2's three, then the file's) go to g2.out. If it does not pass, the fixed line
#          "G2 failed: no upload; C2 stays." goes to g2.verdict and to the standard error, and the run ends with 32
#   done   one SHA256SUMS over every built file, copied into the output folder
# A build step is the builder's command, then `vec check --board T3:bcatenin`, then its stamp: it fails unless the
# stamp says COMMIT and clean True. Builds run at 16 BLAS threads, as the earlier T3 builds; NUMPY_MADVISE_HUGEPAGE=0.
# Every build and `vec check` output goes to /dev/null, and so does the standard error of g1 and g2: a test build is read
# by its exit code, its sha256 and its stamp's commit and clean fields only (AGENTS.md), never its log, size, cell count
# or sidecar record. The run folder of the factor tables is fixed below (read only, no override). Nothing here deletes,
# fetches or uploads anything.
# Exit codes: 90 odd RUN, 91 run folder exists, 96 HEAD is not COMMIT, 98 tree not clean, 95 row R36, 94 row R22, 93
# output folder exists, 10 tests, 11 build failed, 12 the builder refused (its exit 2), 13 vec check, 14 no stamp, 15
# stamp not this commit and clean, 20 21 22 23 24 an identity gate is not its recorded bytes (ids1 ids2 c1s2 c2s2
# c1g2), 30 G1 did not pass, 31 the pair differs, 32 G2 did not pass (its own exit in g2.rc), 33 hash list, 99 the
# launcher stopped unexpectedly (a failing command, a TERM, INT or HUP signal, or an exit that was not one of the
# finishes above).
set -u
COMMIT=${1:?usage: t3_cd7_build.sh COMMIT RUN}
RUN=${2:?usage: t3_cd7_build.sh COMMIT RUN}
case "$RUN" in
  *[!A-Za-z0-9_-]*) echo "refused: RUN takes letters, digits, _ and - only" >&2; exit 90 ;;
esac
export PATH=$HOME/.local/bin:$PATH
S=$HOME/v/scratch/t3cd7/$RUN
mkdir -p "$HOME/v/scratch/t3cd7"
mkdir "$S" 2> /dev/null || { echo "refused: the run folder $S exists (a run folder is never reused)" >&2; exit 91; }
mkdir "$S/tmp"
export NUMPY_MADVISE_HUGEPAGE=0 TMPDIR="$S/tmp"
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
PY="uv run --frozen --no-sync python -I -B"
BOARD=T3:bcatenin
DEC=experiments/t3_decomp/results/gata4/gata4_decisions.json
RUN_DIR=$HOME/v/virtual-embryo-challenge/runs/t3_decomp/gata4_20261028T122425Z  # read only; fixed, no override
WT=data/raw/T2/E8.75.h5ad
FILE=T3_c2_s2_cd7.h5ad
IDS1_SHA=3712c3f9988d0f8aa4198b22c37b866a84ef8757fd500a9163944f0663a520a2  # identity-s1, builds_dc7b6ae.SHA256SUMS
IDS2_SHA=1b69521e3c28bc84185410e838efc6fcbcff898434e73f9a7547208a4c82711f  # identity-s2, builds_dc7b6ae.SHA256SUMS
C1S2_SHA=8718330866429dc089eaffd2f96f66936fe54134eba2d50f8126266e6a3e54b8  # c1-s2, t3two_fb01487.SHA256SUMS
C2S2_SHA=47111ac1d07e73c0c26e990a474dff0de94ae22d6a4474a569bcb674ae2b85a7  # c2-s2 (C2), t3two_fb01487.SHA256SUMS
C1G2_SHA=e01e24facb6841a18a503655bcc5e9816304a993fe9a0dcfeca458c42ead1e16  # T3_c1_s2_g2, t3g2_a3c511c.SHA256SUMS
CD7_ARGS="--form c2-s2 --composition-from allowed-s1 --captain-ruling R36"
STEP=gate
finish() { trap - ERR TERM INT HUP EXIT; echo "$1" > "$S/$STEP.exit"; echo "$1" > "$S/all.exit"; exit "$1"; }
abort() { trap - ERR TERM INT HUP EXIT; echo 99 > "$S/$STEP.exit"; echo 99 > "$S/all.exit"; exit 99; }
begin() { echo running > "$S/$STEP.exit"; }
step_ok() { echo 0 > "$S/$STEP.exit"; }
echo running > "$S/all.exit"
trap abort ERR TERM INT HUP EXIT
begin

head=$(git rev-parse HEAD 2> /dev/null) || head=""
[ "$head" = "$COMMIT" ] || { echo "refused: HEAD is not $COMMIT (give the full hash)" >&2; finish 96; }
[ -z "$(git status --porcelain 2> /dev/null)" ] || { echo "refused: the tree is not clean" >&2; finish 98; }
[ "$(grep -c '^| R36 ' experiments/harness/TRIES.md)" = 1 ] \
  || { echo "refused: TRIES.md does not hold exactly one row R36" >&2; finish 95; }
first=$(grep '^| R36 ' experiments/harness/TRIES.md | cut -d'|' -f2)
for opt in '`--composition-from allowed-s1`' '`--form c2-s2`' '`--captain-ruling R36`'; do
  printf '%s\n' "$first" | grep -qF -- "$opt" \
    || { echo "refused: row R36's first cell does not name $opt" >&2; finish 95; }
done
[ "$(grep -c '^| R22 ' experiments/harness/TRIES.md)" = 1 ] \
  || { echo "refused: TRIES.md does not hold exactly one row R22" >&2; finish 94; }
first=$(grep '^| R22 ' experiments/harness/TRIES.md | cut -d'|' -f2)
for opt in '`--gata4-exponent 2.0`' '`--captain-ruling R22`'; do
  printf '%s\n' "$first" | grep -qF -- "$opt" \
    || { echo "refused: row R22's first cell does not name $opt" >&2; finish 94; }
done
O=runs/t3cd7_${head:0:7}_$RUN
mkdir -p runs
mkdir "$O" 2> /dev/null || { echo "refused: the output folder $O exists" >&2; finish 93; }
echo "$head" > "$S/head"
step_ok

STEP=tests
begin
$PY -m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp="$S/tmp/pytest" \
  tests/test_t3_cd7_build.py tests/test_t3_cd7_old_new.py tests/test_t3_cd7_launcher.py \
  tests/test_t3_two_knockout_build.py tests/test_t3_gata4_weight.py tests/test_t3_barred_exponent_build.py \
  tests/test_t3_combination_build.py tests/test_t3_outside_forms.py experiments/t3_decomp/test_build_submission.py \
  > "$S/tests.out" 2>&1 || finish 10
step_ok

export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16
stamp() {  # the stamp's commit and clean fields only
  $PY -c 'import json,sys; s=json.load(open(sys.argv[1])); print(s["commit"], s["clean"])' "$1.stamp.json"
}
build() {  # folder name args: the builder, then `vec check`; 11 build, 12 refused, 13 check, 14 no stamp, 15 stamp
  local d="$O/$1" name=$2 rc=0 st
  shift 2
  mkdir -p "$d"
  $PY experiments/t3_decomp/build_two_knockout.py build "$DEC" --wt "$WT" --run-dir "$RUN_DIR" --board "$BOARD" \
    --seed 1000001 "$@" --out "$d/$name" > /dev/null 2>&1 || rc=$?
  if [ $rc = 2 ]; then return 12; fi
  if [ $rc != 0 ]; then return 11; fi
  uv run --frozen --no-sync vec check "$d/$name" --board "$BOARD" > /dev/null 2>&1 || return 13
  st=$(stamp "$d/$name" 2> /dev/null) || return 14
  echo "$(basename "$d")/$name $st" >> "$S/stamps"
  [ "$st" = "$COMMIT True" ] || return 15
}
identity() {  # name recorded_sha256 code args: one identity gate, built once into gate/
  local name=$1 ref=$2 code=$3 got
  shift 3
  build gate "$name" "$@" || finish $?
  got=$(sha256sum "$O/gate/$name" | cut -d' ' -f1)
  [ "$got" = "$ref" ] || { echo "stopped: gate/$name is not its recorded bytes" >&2; finish "$code"; }
}
: > "$S/stamps"

STEP=ids1
begin
identity identity-s1.h5ad "$IDS1_SHA" 20 --form identity-s1
step_ok

STEP=ids2
begin
identity identity-s2.h5ad "$IDS2_SHA" 21 --form identity-s2
step_ok

STEP=c1s2
begin
identity c1-s2.h5ad "$C1S2_SHA" 22 --form c1-s2
step_ok

STEP=c2s2
begin
identity c2-s2.h5ad "$C2S2_SHA" 23 --form c2-s2
step_ok

STEP=c1g2
begin
identity T3_c1_s2_g2.h5ad "$C1G2_SHA" 24 --form c1-s2 --gata4-exponent 2.0 --captain-ruling R22
step_ok

STEP=g1
begin
rc=0
$PY experiments/t3_decomp/build_two_knockout.py g1 "$DEC" --wt "$WT" --composition-from allowed-s1 \
  > "$S/g1.out" 2> /dev/null || rc=$?
[ $rc = 0 ] || { echo "stopped: G1 did not pass (exit $rc; g1.out)" >&2; finish 30; }
step_ok

for STEP in a b; do begin; build "$STEP" "$FILE" $CD7_ARGS || finish $?; step_ok; done

STEP=pair
begin
cmp -s "$O/a/$FILE" "$O/b/$FILE" || { echo "stopped: the two builds of $FILE differ" >&2; finish 31; }
step_ok

STEP=g2
begin
rc=0
$PY experiments/t3_decomp/build_two_knockout.py g2 "$DEC" --wt "$WT" --run-dir "$RUN_DIR" --board "$BOARD" \
  --seed 1000001 --composition-from allowed-s1 > "$S/g2.out" 2> /dev/null || rc=$?
if [ $rc != 0 ]; then
  echo "$rc" > "$S/g2.rc"
  echo "G2 failed: no upload; C2 stays." > "$S/g2.verdict"
  echo "G2 failed: no upload; C2 stays." >&2
  finish 32
fi
step_ok

STEP=done
begin
(cd "$O" && sha256sum gate/*.h5ad a/*.h5ad b/*.h5ad) > "$S/SHA256SUMS" || finish 33
cp "$S/SHA256SUMS" "$O/SHA256SUMS" || finish 33
finish 0
