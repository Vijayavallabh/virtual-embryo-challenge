#!/usr/bin/env bash
# The depth line's files (row R53, the captain's order of 2026-11-24: his named exception on T3's depth for
# T3_c2_s2_dd, T3_c2_s2_ddh and T3_c2_s2_dg, bound to row K30's read; the T3 ERROR design of the step-change round,
# proposal P1, dD, as its cold challenge fixes it): C2's command (build_two_knockout.py --form c2-s2, recorded sha256
# 47111ac1d07e...) with `--captain-ruling R53` and
#   dd   T3_c2_s2_dd.h5ad    `--depth-from allowed-s1`                        the frozen per-label depth ratios
#   ddh  T3_c2_s2_ddh.h5ad   `--depth-from allowed-s1 --depth-exponent 0.5`   the ratios ** 0.5
#   dg   T3_c2_s2_dg.h5ad    `--depth-from allowed-s1 --depth-form global`    one ratio for every label
# ARM names one file, or all: the three in that order, after one set of identity gates. The builder's own checks run
# inside each build (the counts the depth step receives are C2's row by row, each cell's ratio its label's: exit 7,
# nothing written). On the pattern of t3_gonly_build.sh. Test-board builds on a DGX campaign checkout (never the
# reference checkout), on the machine that built C2 (its identity gate). Run from the checkout's root, detached
# (setsid nohup ... &), and wait on the exit files with a background poller; never a foreground sleep:
#   bash experiments/t3_decomp/t3_depth_build.sh <COMMIT> <RUN> <ROW> <ARM>
# ROW is the ruling row, R53 (the only row the builder takes for this line); ARM is dd, ddh, dg or all. RUN names this
# run's scratch folder ~/v/scratch/t3depth/<RUN> and its output folder runs/t3depth_<commit7>_<RUN>. Neither may exist:
# a run folder is never reused and nothing is ever deleted here, so a second run takes a new name.
# In order, each step first writing the word running to ~/v/scratch/t3depth/<RUN>/<STEP>.exit, then its exit code; the
# code that ends the run also goes to all.exit (running until then):
#   gate   HEAD is COMMIT (the full hash) and the tree is clean (git status --porcelain prints nothing); TRIES.md holds
#          exactly one row R22 whose first cell holds `--form c2-s2 --gata4-exponent 2.0`, `--gata4-exponent 2.0` and
#          `--captain-ruling R22`, each as one whole backtick-delimited string (the last two identity gates), and
#          exactly one row ROW whose first cell opens "ROW (the captain's order", says "ruling, no run" and holds each
#          named file's option string, `--form c2-s2` and `--captain-ruling ROW`, each as one whole backtick-delimited
#          string (the builder's own check); the output folder does not exist; then it writes the head and arm files
#   tests  the builder's and this launcher's targeted tests (synthetic data, and row R53's first cell read by code),
#          8 BLAS threads
#   ids1 ids2 c1s2 c2s2 c1g2 c2g2   the identity gates, each command built once into gate/, as a build step below,
#          and its sha256 the recorded one: identity-s1 and identity-s2 (the barred s1 and s2,
#          builds_dc7b6ae.SHA256SUMS), c1-s2 and c2-s2 (C1 and C2, t3two_fb01487.SHA256SUMS), c1-s2 with
#          --gata4-exponent 2.0 --captain-ruling R22 (T3_c1_s2_g2.h5ad, t3g2_a3c511c.SHA256SUMS) and c2-s2 with
#          --gata4-exponent 2.0 --captain-ruling R22 (T3_c2_s2_g2.h5ad, t3c2g_fc57ca9_g1.SHA256SUMS)
#   dda ddb, ddha ddhb, dga dgb   each named file built twice, into folders a and b
#   pair   each named file's two builds byte for byte (one sha256 per file)
#   done   one SHA256SUMS over every built file, copied into the output folder
# A build step is the builder's command, then `vec check --board T3:bcatenin`, then its stamp: it fails unless the
# stamp says COMMIT and clean True. Builds run at 16 BLAS threads, as the earlier T3 builds; NUMPY_MADVISE_HUGEPAGE=0.
# Every build and `vec check` output goes to /dev/null: a test build is read by its exit code, its sha256 and its
# stamp's commit and clean fields only (AGENTS.md), never its log, size, cell count or sidecar record. The run folder
# of the factor tables is fixed below (read only, no override). Nothing here deletes, fetches or uploads anything.
# Exit codes: 88 odd ARM, 89 a ROW other than R53 or not four arguments, 90 odd RUN, 91 run folder exists, 96 HEAD is
# not COMMIT, 98 tree not clean, 95 row R22, 94 row ROW, 93 output folder exists, 10 tests, 11 build failed, 12 the
# builder refused (its exit 2), 13 vec check, 14 no stamp, 15 stamp not this commit and clean, 16 line 3's guard stopped
# an identity gate (the builder's exit 6), 17 the factor is not the one the line's figures describe (the builder's exit
# 5), 18 the depth line's checks stopped the file (the builder's exit 7, nothing written), 20 21 22 23 24 25 an identity
# gate is not its recorded bytes (ids1 ids2 c1s2 c2s2 c1g2 c2g2), 31 32 33 a pair differs (dd ddh dg), 39 hash list, 99
# the launcher stopped unexpectedly (a failing command, a TERM, INT or HUP signal, or an exit that was not one of the
# finishes above).
set -u
set -o pipefail
COMMIT=${1:?usage: t3_depth_build.sh COMMIT RUN ROW ARM}
RUN=${2:?usage: t3_depth_build.sh COMMIT RUN ROW ARM}
ROW=${3:?usage: t3_depth_build.sh COMMIT RUN ROW ARM}
ARM=${4:?usage: t3_depth_build.sh COMMIT RUN ROW ARM}
[ $# -eq 4 ] || { echo "refused: four arguments (COMMIT RUN ROW ARM)" >&2; exit 89; }
case "$RUN" in
  *[!A-Za-z0-9_-]*) echo "refused: RUN takes letters, digits, _ and - only" >&2; exit 90 ;;
esac
[ "$ROW" = R53 ] || { echo "refused: ROW is R53, the one ruling row the builder takes for this line" >&2; exit 89; }
case "$ARM" in
  dd | ddh | dg) ARMS=$ARM ;;
  all) ARMS="dd ddh dg" ;;
  *) echo "refused: ARM takes dd, ddh, dg or all" >&2; exit 88 ;;
esac
export PATH=$HOME/.local/bin:$PATH
S=$HOME/v/scratch/t3depth/$RUN
mkdir -p "$HOME/v/scratch/t3depth"
mkdir "$S" 2> /dev/null || { echo "refused: the run folder $S exists (a run folder is never reused)" >&2; exit 91; }
mkdir "$S/tmp"
export NUMPY_MADVISE_HUGEPAGE=0 TMPDIR="$S/tmp"
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
PY="uv run --frozen --no-sync python -I -B"
BOARD=T3:bcatenin
DEC=experiments/t3_decomp/results/gata4/gata4_decisions.json
RUN_DIR=$HOME/v/virtual-embryo-challenge/runs/t3_decomp/gata4_20261028T122425Z  # read only; fixed, no override
WT=data/raw/T2/E8.75.h5ad
FORM_OPT="--form c2-s2"
DD_OPT="--depth-from allowed-s1"
DDH_OPT="--depth-from allowed-s1 --depth-exponent 0.5"
DG_OPT="--depth-from allowed-s1 --depth-form global"
DD=T3_c2_s2_dd.h5ad
DDH=T3_c2_s2_ddh.h5ad
DG=T3_c2_s2_dg.h5ad
IDS1_SHA=3712c3f9988d0f8aa4198b22c37b866a84ef8757fd500a9163944f0663a520a2  # identity-s1, builds_dc7b6ae.SHA256SUMS
IDS2_SHA=1b69521e3c28bc84185410e838efc6fcbcff898434e73f9a7547208a4c82711f  # identity-s2, builds_dc7b6ae.SHA256SUMS
C1S2_SHA=8718330866429dc089eaffd2f96f66936fe54134eba2d50f8126266e6a3e54b8  # c1-s2, t3two_fb01487.SHA256SUMS
C2S2_SHA=47111ac1d07e73c0c26e990a474dff0de94ae22d6a4474a569bcb674ae2b85a7  # c2-s2 (C2), t3two_fb01487.SHA256SUMS
C1G2_SHA=e01e24facb6841a18a503655bcc5e9816304a993fe9a0dcfeca458c42ead1e16  # T3_c1_s2_g2, t3g2_a3c511c.SHA256SUMS
C2G2_SHA=014d2bf9d8ba7fdc56dcfba2a3b8aa77f0ea20bb806ebc2c1193e483d62712ac  # T3_c2_s2_g2, t3c2g_fc57ca9_g1.SHA256SUMS
STEP=gate
finish() { trap - ERR TERM INT HUP EXIT; echo "$1" > "$S/$STEP.exit"; echo "$1" > "$S/all.exit"; exit "$1"; }
abort() { trap - ERR TERM INT HUP EXIT; echo 99 > "$S/$STEP.exit"; echo 99 > "$S/all.exit"; exit 99; }
begin() { echo running > "$S/$STEP.exit"; }
step_ok() { echo 0 > "$S/$STEP.exit"; }
names() {  # cell string: whether the cell holds the string as one whole backtick-delimited string (no pipe)
  local rest=$1 tok
  while [[ $rest == *'`'*'`'* ]]; do
    rest=${rest#*'`'}
    tok=${rest%%'`'*}
    rest=${rest#*'`'}
    [ "$tok" = "$2" ] && return 0
  done
  return 1
}
arm_opt() { case "$1" in dd) printf '%s' "$DD_OPT" ;; ddh) printf '%s' "$DDH_OPT" ;; dg) printf '%s' "$DG_OPT" ;; esac; }
arm_file() { case "$1" in dd) printf '%s' "$DD" ;; ddh) printf '%s' "$DDH" ;; dg) printf '%s' "$DG" ;; esac; }
arm_pair() { case "$1" in dd) echo 31 ;; ddh) echo 32 ;; dg) echo 33 ;; esac; }
echo running > "$S/all.exit"
trap abort ERR TERM INT HUP EXIT
begin

head=$(git rev-parse HEAD 2> /dev/null) || head=""
[ "$head" = "$COMMIT" ] || { echo "refused: HEAD is not $COMMIT (give the full hash)" >&2; finish 96; }
[ -z "$(git status --porcelain 2> /dev/null)" ] || { echo "refused: the tree is not clean" >&2; finish 98; }
[ "$(grep -c '^| R22 ' experiments/harness/TRIES.md)" = 1 ] \
  || { echo "refused: TRIES.md does not hold exactly one row R22" >&2; finish 95; }
first=$(grep '^| R22 ' experiments/harness/TRIES.md | cut -d'|' -f2)
for opt in '--form c2-s2 --gata4-exponent 2.0' '--gata4-exponent 2.0' '--captain-ruling R22'; do
  if ! names "$first" "$opt"; then
    echo "refused: row R22's first cell does not name \`$opt\` as one whole string" >&2
    finish 95
  fi
done
[ "$(grep -c "^| $ROW " experiments/harness/TRIES.md)" = 1 ] \
  || { echo "refused: TRIES.md does not hold exactly one row $ROW" >&2; finish 94; }
first=$(grep "^| $ROW " experiments/harness/TRIES.md | cut -d'|' -f2)
[[ $first == " $ROW (the captain's order"* ]] \
  || { echo "refused: row $ROW's first cell does not open: $ROW (the captain's order" >&2; finish 94; }
[[ $first == *"ruling, no run"* ]] || { echo "refused: row $ROW's first cell does not say: ruling, no run" >&2; finish 94; }
for a in $ARMS; do
  for opt in "$(arm_opt "$a")" "$FORM_OPT" "--captain-ruling $ROW"; do
    if ! names "$first" "$opt"; then
      echo "refused: row $ROW's first cell does not name \`$opt\` as one whole string" >&2
      finish 94
    fi
  done
done
O=runs/t3depth_${head:0:7}_$RUN
mkdir -p runs
mkdir "$O" 2> /dev/null || { echo "refused: the output folder $O exists" >&2; finish 93; }
echo "$head" > "$S/head"
echo "$ROW $ARM:$(for a in $ARMS; do printf ' %s' "$(arm_file "$a")"; done)" > "$S/arm"
step_ok

STEP=tests
begin
$PY -m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp="$S/tmp/pytest" \
  tests/test_t3_depth_build.py tests/test_t3_depth_old_new.py tests/test_t3_depth_launcher.py \
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
build() {  # folder name args: the builder, then `vec check`; 11 build, 12 refused, 16 guard, 17 factor, 18 depth
  local d="$O/$1" name=$2 rc=0 st
  shift 2
  mkdir -p "$d"
  $PY experiments/t3_decomp/build_two_knockout.py build "$DEC" --wt "$WT" --run-dir "$RUN_DIR" --board "$BOARD" \
    --seed 1000001 "$@" --out "$d/$name" > /dev/null 2>&1 || rc=$?
  if [ $rc = 2 ]; then return 12; fi
  if [ $rc = 6 ]; then return 16; fi
  if [ $rc = 5 ]; then return 17; fi
  if [ $rc = 7 ]; then return 18; fi
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
  got=$(sha256sum "$O/gate/$name" | cut -d' ' -f1) || got=""
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

for a in $ARMS; do
  for side in a b; do
    STEP=$a$side
    begin
    # shellcheck disable=SC2046  # the option string splits into its words, as the gonly launcher's _ARGS do
    build "$side" "$(arm_file "$a")" $FORM_OPT $(arm_opt "$a") --captain-ruling "$ROW" || finish $?
    step_ok
  done
done

STEP=pair
begin
for a in $ARMS; do
  f=$(arm_file "$a")
  cmp -s "$O/a/$f" "$O/b/$f" || { echo "stopped: the two builds of $f differ" >&2; finish "$(arm_pair "$a")"; }
done
step_ok

STEP=done
begin
(cd "$O" && sha256sum gate/*.h5ad a/*.h5ad b/*.h5ad) > "$S/SHA256SUMS" || finish 39
cp "$S/SHA256SUMS" "$O/SHA256SUMS" || finish 39
finish 0
