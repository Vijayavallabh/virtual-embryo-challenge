#!/usr/bin/env bash
# The three files of the options on C2's command (the T3 design of the design round of 2026-11-24, section 2; row R22's
# entry naming `--gata4-exponent inf` and `--form c2-s2 --size-scale 0.67`): C2's command (build_two_knockout.py
# --form c2-s2, recorded sha256 47111ac1d07e...) with `--captain-ruling R22` and
#   F1  T3_c2_s2_gonly.h5ad        --gata4-exponent inf                     g = inf: phi_G alone at C2's held size
#   F2  T3_c2_s2_s067.h5ad         --size-scale 0.67                        C2's ln factor times 0.67
#   F3  T3_c2_s2_gonly_s067.h5ad   --gata4-exponent inf --size-scale 0.67   F1's ln factor times 0.67
# The builder's two target-free guards (range on each file, distinctness on F1; exit 6, nothing written) and its factor
# check (exit 5) run inside each build. On the pattern of t3_c2g_build.sh. Test-board builds on a DGX campaign checkout
# (never the reference checkout). Run from the checkout's root, detached (setsid nohup ... &), and wait on the exit
# files with a background poller; never a foreground sleep:
#   bash experiments/t3_decomp/t3_gonly_build.sh <COMMIT> <RUN>
# RUN names this run's scratch folder ~/v/scratch/t3gonly/<RUN> and its output folder runs/t3gonly_<commit7>_<RUN>.
# Neither may exist: a run folder is never reused and nothing is ever deleted here, so a second run takes a new name.
# In order, each step first writing the word running to ~/v/scratch/t3gonly/<RUN>/<STEP>.exit, then its exit code; the
# code that ends the run also goes to all.exit (running until then):
#   gate   HEAD is COMMIT (the full hash) and the tree is clean (git status --porcelain prints nothing); TRIES.md holds
#          exactly one row R22 whose first cell names `--gata4-exponent inf` (F1 and F3), `--form c2-s2 --size-scale
#          0.67` (F2 and F3), `--form c2-s2 --gata4-exponent 2.0` and `--gata4-exponent 2.0` (lines 3 and 2: the last
#          two identity gates) and `--captain-ruling R22`; the output folder does not exist; then it writes the head
#          file
#   tests  the builder's and this launcher's targeted tests (synthetic data only), 8 BLAS threads
#   ids1 ids2 c1s2 c2s2 c1g2 c2g2   the identity gates, each command built once into gate/, as a build step below,
#          and its sha256 the recorded one: identity-s1 and identity-s2 (the barred s1 and s2,
#          builds_dc7b6ae.SHA256SUMS), c1-s2 and c2-s2 (C1 and C2, t3two_fb01487.SHA256SUMS), c1-s2 with
#          --gata4-exponent 2.0 --captain-ruling R22 (T3_c1_s2_g2.h5ad, t3g2_a3c511c.SHA256SUMS) and c2-s2 with
#          --gata4-exponent 2.0 --captain-ruling R22 (T3_c2_s2_g2.h5ad, t3c2g_fc57ca9_g1.SHA256SUMS)
#   f1a f1b f2a f2b f3a f3b   each file built twice, into folders a and b
#   pair   each file's two builds byte for byte (one sha256 per file)
#   done   one SHA256SUMS over every built file, copied into the output folder
# A build step is the builder's command, then `vec check --board T3:bcatenin`, then its stamp: it fails unless the
# stamp says COMMIT and clean True. Builds run at 16 BLAS threads, as the earlier T3 builds; NUMPY_MADVISE_HUGEPAGE=0.
# Every build and `vec check` output goes to /dev/null: a test build is read by its exit code, its sha256 and its
# stamp's commit and clean fields only (AGENTS.md), never its log, size, cell count or sidecar record. The run folder
# of the factor tables is fixed below (read only, no override). Nothing here deletes, fetches or uploads anything.
# Exit codes: 90 odd RUN, 91 run folder exists, 96 HEAD is not COMMIT, 98 tree not clean, 95 row R22, 93 output folder
# exists, 10 tests, 11 build failed, 12 the builder refused (its exit 2), 13 vec check, 14 no stamp, 15 stamp not this
# commit and clean, 16 a guard stopped the file (the builder's exit 6, nothing written), 17 the factor is not the one
# the guards describe (the builder's exit 5, nothing written), 20 21 22 23 24 25 an identity gate is not its recorded
# bytes (ids1 ids2 c1s2 c2s2 c1g2 c2g2), 31 32 33 a pair differs (F1 F2 F3), 39 hash list, 99 the launcher stopped
# unexpectedly (a failing command, a TERM, INT or HUP signal, or an exit that was not one of the finishes above).
set -u
COMMIT=${1:?usage: t3_gonly_build.sh COMMIT RUN}
RUN=${2:?usage: t3_gonly_build.sh COMMIT RUN}
case "$RUN" in
  *[!A-Za-z0-9_-]*) echo "refused: RUN takes letters, digits, _ and - only" >&2; exit 90 ;;
esac
export PATH=$HOME/.local/bin:$PATH
S=$HOME/v/scratch/t3gonly/$RUN
mkdir -p "$HOME/v/scratch/t3gonly"
mkdir "$S" 2> /dev/null || { echo "refused: the run folder $S exists (a run folder is never reused)" >&2; exit 91; }
mkdir "$S/tmp"
export NUMPY_MADVISE_HUGEPAGE=0 TMPDIR="$S/tmp"
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
PY="uv run --frozen --no-sync python -I -B"
BOARD=T3:bcatenin
DEC=experiments/t3_decomp/results/gata4/gata4_decisions.json
RUN_DIR=$HOME/v/virtual-embryo-challenge/runs/t3_decomp/gata4_20261028T122425Z  # read only; fixed, no override
WT=data/raw/T2/E8.75.h5ad
F1=T3_c2_s2_gonly.h5ad
F2=T3_c2_s2_s067.h5ad
F3=T3_c2_s2_gonly_s067.h5ad
IDS1_SHA=3712c3f9988d0f8aa4198b22c37b866a84ef8757fd500a9163944f0663a520a2  # identity-s1, builds_dc7b6ae.SHA256SUMS
IDS2_SHA=1b69521e3c28bc84185410e838efc6fcbcff898434e73f9a7547208a4c82711f  # identity-s2, builds_dc7b6ae.SHA256SUMS
C1S2_SHA=8718330866429dc089eaffd2f96f66936fe54134eba2d50f8126266e6a3e54b8  # c1-s2, t3two_fb01487.SHA256SUMS
C2S2_SHA=47111ac1d07e73c0c26e990a474dff0de94ae22d6a4474a569bcb674ae2b85a7  # c2-s2 (C2), t3two_fb01487.SHA256SUMS
C1G2_SHA=e01e24facb6841a18a503655bcc5e9816304a993fe9a0dcfeca458c42ead1e16  # T3_c1_s2_g2, t3g2_a3c511c.SHA256SUMS
C2G2_SHA=014d2bf9d8ba7fdc56dcfba2a3b8aa77f0ea20bb806ebc2c1193e483d62712ac  # T3_c2_s2_g2, t3c2g_fc57ca9_g1.SHA256SUMS
F1_ARGS="--form c2-s2 --gata4-exponent inf --captain-ruling R22"
F2_ARGS="--form c2-s2 --size-scale 0.67 --captain-ruling R22"
F3_ARGS="--form c2-s2 --gata4-exponent inf --size-scale 0.67 --captain-ruling R22"
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
[ "$(grep -c '^| R22 ' experiments/harness/TRIES.md)" = 1 ] \
  || { echo "refused: TRIES.md does not hold exactly one row R22" >&2; finish 95; }
first=$(grep '^| R22 ' experiments/harness/TRIES.md | cut -d'|' -f2)
for opt in '`--gata4-exponent inf`' '`--form c2-s2 --size-scale 0.67`' '`--form c2-s2 --gata4-exponent 2.0`' \
    '`--gata4-exponent 2.0`' '`--captain-ruling R22`'; do
  printf '%s\n' "$first" | grep -qF -- "$opt" \
    || { echo "refused: row R22's first cell does not name $opt" >&2; finish 95; }
done
O=runs/t3gonly_${head:0:7}_$RUN
mkdir -p runs
mkdir "$O" 2> /dev/null || { echo "refused: the output folder $O exists" >&2; finish 93; }
echo "$head" > "$S/head"
step_ok

STEP=tests
begin
$PY -m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp="$S/tmp/pytest" \
  tests/test_t3_gonly_build.py tests/test_t3_gonly_old_new.py tests/test_t3_gonly_launcher.py \
  tests/test_t3_c2g_build.py tests/test_t3_c2g_old_new.py tests/test_t3_c2g_launcher.py \
  tests/test_t3_cd7_build.py tests/test_t3_cd7_old_new.py tests/test_t3_cd7_launcher.py \
  tests/test_t3_two_knockout_build.py tests/test_t3_gata4_weight.py tests/test_t3_barred_exponent_build.py \
  tests/test_t3_combination_build.py tests/test_t3_outside_forms.py experiments/t3_decomp/test_build_submission.py \
  > "$S/tests.out" 2>&1 || finish 10
step_ok

export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16
stamp() {  # the stamp's commit and clean fields only
  $PY -c 'import json,sys; s=json.load(open(sys.argv[1])); print(s["commit"], s["clean"])' "$1.stamp.json"
}
build() {  # folder name args: the builder, then `vec check`; 11 build, 12 refused, 16 guard, 17 factor, 13 14 15
  local d="$O/$1" name=$2 rc=0 st
  shift 2
  mkdir -p "$d"
  $PY experiments/t3_decomp/build_two_knockout.py build "$DEC" --wt "$WT" --run-dir "$RUN_DIR" --board "$BOARD" \
    --seed 1000001 "$@" --out "$d/$name" > /dev/null 2>&1 || rc=$?
  if [ $rc = 2 ]; then return 12; fi
  if [ $rc = 6 ]; then return 16; fi
  if [ $rc = 5 ]; then return 17; fi
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

STEP=c2g2
begin
identity T3_c2_s2_g2.h5ad "$C2G2_SHA" 25 --form c2-s2 --gata4-exponent 2.0 --captain-ruling R22
step_ok

for STEP in f1a f1b; do begin; build "${STEP:2}" "$F1" $F1_ARGS || finish $?; step_ok; done
for STEP in f2a f2b; do begin; build "${STEP:2}" "$F2" $F2_ARGS || finish $?; step_ok; done
for STEP in f3a f3b; do begin; build "${STEP:2}" "$F3" $F3_ARGS || finish $?; step_ok; done

STEP=pair
begin
cmp -s "$O/a/$F1" "$O/b/$F1" || { echo "stopped: the two builds of $F1 differ" >&2; finish 31; }
cmp -s "$O/a/$F2" "$O/b/$F2" || { echo "stopped: the two builds of $F2 differ" >&2; finish 32; }
cmp -s "$O/a/$F3" "$O/b/$F3" || { echo "stopped: the two builds of $F3 differ" >&2; finish 33; }
step_ok

STEP=done
begin
(cd "$O" && sha256sum gate/*.h5ad a/*.h5ad b/*.h5ad) > "$S/SHA256SUMS" || finish 39
cp "$S/SHA256SUMS" "$O/SHA256SUMS" || finish 39
finish 0
