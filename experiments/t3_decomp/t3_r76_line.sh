#!/usr/bin/env bash
# Row R76's five files (the captain's order of 2026-11-26: T3's refinement round on C2's command, the T3 REF3 design as
# its cold challenge fixes it; his named exception on T3's depth for each file), one file per call, each C2's command
# (build_two_knockout.py --form c2-s2, recorded sha256 47111ac1d07e...) with one of row R76's option strings, only with
# `--captain-ruling R76`, and with `--den-ruling R76` for the three strings the row says "the 5% check gates nothing
# for it" of:
#   T3_c2_s2_s085      `--form c2-s2 --size-scale 0.85`                        --den-ruling R76
#   T3_c2_s2_s05       `--form c2-s2 --size-scale 0.5`                         --den-ruling R76
#   T3_c2_s2_g15       `--form c2-s2 --gata4-exponent 1.5`
#   T3_c2_s2_g05       `--form c2-s2 --gata4-exponent 0.5`
#   T3_c2_s2_g2_s067   `--form c2-s2 --gata4-exponent 2.0 --size-scale 0.67`   --den-ruling R76
# The builder's own checks run inside each build (the range guard and the 5% den check, exit 6, and the factor check,
# exit 5; nothing written). Test-board builds on a DGX campaign checkout with its own venv, made there by
# `uv sync --frozen` (never the reference checkout), on the machine that built the identity files' recorded bytes: a T3
# file built on one machine is never checked against a hash made on the other. Run from the root of that checkout,
# detached (setsid nohup ... &), and wait on the exit files with a background poller; never a foreground sleep:
#   bash experiments/t3_decomp/t3_r76_line.sh <COMMIT> <RUN> <FILE_NAME>
# FILE_NAME is one of the five in the literal case table below (tests/test_t3_r76_launcher.py holds it equal to the
# builder's R76_TABLE, the one table row R76's strings are checked against); any other name, or a fourth argument, is
# refused before any folder is made. RUN names this run's scratch folder ~/v/scratch/t3r76/<RUN> and its output folder
# runs/t3r76_<commit7>_<RUN>. Neither may exist: a run folder is never reused and nothing is ever deleted here, so a
# second file, or a second run of one, takes a new name. In order, each step first writing the word running to
# ~/v/scratch/t3r76/<RUN>/<STEP>.exit, then its exit code; the code that ends the run also goes to all.exit (running
# until then):
#   gate   HEAD is COMMIT (the full hash), the tracked tree is clean, TRIES.md holds exactly one row R76 whose first
#          cell starts "R76 (the captain's ", says "ruling, no run" and names FILE_NAME's option string and
#          `--captain-ruling R76` (and `--den-ruling R76` for the three lifted strings), each as a whole backticked
#          string, the option string followed there by " (FILE_NAME)" as the row pairs each string with its file; the
#          data (WT E8.75) and the decision file are files; the gata4 run folder holds the two factor tables at the
#          builder's pinned sha256; .venv's vec.pth names this checkout's src; then it writes the head file
#   tests  the builder's and this launcher's targeted tests (synthetic data only), 8 BLAS threads
#   c2 c2g2 s067 gonly   four identity rebuilds into gate/, each once, each of which must give its recorded bytes or
#          everything stops: C2 (`--form c2-s2`, t3two_fb01487.SHA256SUMS), T3_c2_s2_g2 (`--form c2-s2
#          --gata4-exponent 2.0 --captain-ruling R22`, t3c2g_fc57ca9_g1.SHA256SUMS), T3_c2_s2_s067 and T3_c2_s2_gonly
#          (`--form c2-s2 --size-scale 0.67` and `--form c2-s2 --gata4-exponent inf`, each with `--captain-ruling R22
#          --den-ruling R47`, t3gonly_466e94d_g1.SHA256SUMS): R76's code sits beside the g = 1 and g = inf branches
#          and its last file reuses line 3's weight
#   a b    FILE_NAME by its command, twice, into folders a and b
#   done   one sha256 list over every built file (in the scratch folder and the output folder), then the two builds of
#          FILE_NAME compared byte for byte
# A build step is the builder's command, then `vec check --board T3:bcatenin`, then its stamp: it fails unless the
# stamp says COMMIT and clean True. Builds run at 16 BLAS threads, as the earlier T3 builds; NUMPY_MADVISE_HUGEPAGE=0.
# Every builder and `vec check` output goes to /dev/null: a test build is read by its exit code, its sha256 and its
# stamp's commit and clean fields only (AGENTS.md), never its log, size, cell count or sidecar record. The run folder
# of the factor tables is fixed below (read only, no override). Nothing here deletes, fetches or uploads anything.
# Exit codes: 90 odd RUN, 89 FILE_NAME not one of row R76's five (or not three arguments), 91 run folder exists,
# 96 HEAD is not COMMIT, 98 tracked tree not clean, 95 row R76 not as above, 92 no WT E8.75 or decision file, 87 the
# gata4 run folder or its tables not at their pins, 88 .venv not this checkout's, 93 output folder exists, 10 tests,
# 11 build failed, 12 the builder refused (its exit 2), 13 vec check, 14 no stamp, 15 stamp not this commit and clean,
# 16 a check stopped the file (the builder's exit 6: the range guard or the 5% den check, nothing written; a stop, not
# a fault), 17 the factor is not the one the checks describe (the builder's exit 5, nothing written), 21 22 23 24 an
# identity rebuild is not its recorded bytes (c2 c2g2 s067 gonly), 40 hash list, 41 the two builds differ, 99 the
# launcher stopped unexpectedly (a failing command, a TERM, INT or HUP signal, or an exit that was not one of the
# finishes above).
set -u
set -o pipefail
COMMIT=${1:?usage: t3_r76_line.sh COMMIT RUN FILE_NAME}
RUN=${2:?usage: t3_r76_line.sh COMMIT RUN FILE_NAME}
FILE=${3:?usage: t3_r76_line.sh COMMIT RUN FILE_NAME}
[ "$#" = 3 ] || { echo "refused: give COMMIT RUN FILE_NAME and nothing else" >&2; exit 89; }
case "$RUN" in
  *[!A-Za-z0-9_-]*) echo "refused: RUN takes letters, digits, _ and - only" >&2; exit 90 ;;
esac
# FILE_NAME -> its option string, as row R76 pairs them (the builder's R76_TABLE, in the row's order), and whether the
# row lifts the 5% check for it (then the file is built with --den-ruling R76); a literal table: no name is read by a
# pattern
case "$FILE" in
  T3_c2_s2_s085) OPT="--form c2-s2 --size-scale 0.85"; LIFT=yes ;;
  T3_c2_s2_s05) OPT="--form c2-s2 --size-scale 0.5"; LIFT=yes ;;
  T3_c2_s2_g15) OPT="--form c2-s2 --gata4-exponent 1.5"; LIFT=no ;;
  T3_c2_s2_g05) OPT="--form c2-s2 --gata4-exponent 0.5"; LIFT=no ;;
  T3_c2_s2_g2_s067) OPT="--form c2-s2 --gata4-exponent 2.0 --size-scale 0.67"; LIFT=yes ;;
  *) echo "refused: FILE_NAME is not one of row R76's five files" >&2; exit 89 ;;
esac
ARGS="$OPT --captain-ruling R76"
if [ "$LIFT" = yes ]; then ARGS="$ARGS --den-ruling R76"; fi
export PATH=$HOME/.local/bin:$PATH
unset UV_PROJECT_ENVIRONMENT  # uv runs this checkout's own .venv, the one the gate checks
S=$HOME/v/scratch/t3r76/$RUN
mkdir -p "$HOME/v/scratch/t3r76"
mkdir "$S" 2> /dev/null || { echo "refused: the run folder $S exists (a run folder is never reused)" >&2; exit 91; }
mkdir "$S/tmp"
export NUMPY_MADVISE_HUGEPAGE=0 TMPDIR="$S/tmp"
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
PY="uv run --frozen --no-sync python -I -B"
BOARD=T3:bcatenin
DEC=experiments/t3_decomp/results/gata4/gata4_decisions.json
RUN_DIR=$HOME/v/virtual-embryo-challenge/runs/t3_decomp/gata4_20261028T122425Z  # read only; fixed, no override
WT=data/raw/T2/E8.75.h5ad
TABLE_M_SHA=64b1e0b8572c7193ed1e1f0b7339a26081e601f1e706e6b08193e1718dc5c7e8  # transfer_mab21l2_all.csv, the pin
TABLE_G_SHA=36f2d65a36d2d4fd5307cddfadb56e9724c362c982795f30428eda752b4734ef  # gata4_pooled.csv, the pin
C2S2_SHA=47111ac1d07e73c0c26e990a474dff0de94ae22d6a4474a569bcb674ae2b85a7  # c2-s2.h5ad (C2), t3two_fb01487.SHA256SUMS
C2G2_SHA=014d2bf9d8ba7fdc56dcfba2a3b8aa77f0ea20bb806ebc2c1193e483d62712ac  # T3_c2_s2_g2, t3c2g_fc57ca9_g1.SHA256SUMS
S067_SHA=808700e11c25065ce45a6bea2b50eb15c132e02d314e826c4e05fb37bb6331cd  # T3_c2_s2_s067, t3gonly_466e94d_g1
GONLY_SHA=e3c1a9e6abed8eaea4e3d047d4aac625779144c86e305005f90f76468a18b0ae  # T3_c2_s2_gonly, t3gonly_466e94d_g1
C2G2_ARGS="--form c2-s2 --gata4-exponent 2.0 --captain-ruling R22"
S067_ARGS="--form c2-s2 --size-scale 0.67 --captain-ruling R22 --den-ruling R47"
GONLY_ARGS="--form c2-s2 --gata4-exponent inf --captain-ruling R22 --den-ruling R47"
BT='`'
STEP=gate
finish() { trap - ERR TERM INT HUP EXIT; echo "$1" > "$S/$STEP.exit"; echo "$1" > "$S/all.exit"; exit "$1"; }
abort() { trap - ERR TERM INT HUP EXIT; echo 99 > "$S/$STEP.exit"; echo 99 > "$S/all.exit"; exit 99; }
begin() { echo running > "$S/$STEP.exit"; }
step_ok() { echo 0 > "$S/$STEP.exit"; }
names() {  # cell string: whether the cell holds the string as one whole backtick-delimited string (no pipe)
  local rest=$1 tok
  while [[ $rest == *"$BT"*"$BT"* ]]; do
    rest=${rest#*"$BT"}
    tok=${rest%%"$BT"*}
    rest=${rest#*"$BT"}
    [ "$tok" = "$2" ] && return 0
  done
  return 1
}
pinned() {  # PATH PIN: PATH is a file whose sha256 is PIN
  local got
  [ -f "$1" ] || return 1
  got=$(sha256sum < "$1" | cut -d' ' -f1) || return 1
  [ "$got" = "$2" ]
}
echo running > "$S/all.exit"
trap abort ERR TERM INT HUP EXIT
begin

head=$(git rev-parse HEAD 2> /dev/null) || head=""
[ "$head" = "$COMMIT" ] || { echo "refused: HEAD is not $COMMIT (give the full hash)" >&2; finish 96; }
tracked=$(git status --porcelain --untracked-files=no 2> /dev/null) || tracked="unknown"
if [ -n "$tracked" ]; then
  echo "refused: tracked files differ from the commit (the build stamps must say clean)" >&2
  finish 98
fi
n76=$(grep -c "^| R76 " experiments/harness/TRIES.md 2> /dev/null) || n76=0
[ "$n76" = 1 ] || { echo "refused: TRIES.md does not hold exactly one row R76" >&2; finish 95; }
first=$(grep "^| R76 " experiments/harness/TRIES.md | cut -d'|' -f2) || first=""
[[ $first == " R76 (the captain's "* ]] \
  || { echo "refused: row R76's first cell is not the captain's ruling" >&2; finish 95; }
[[ $first == *"ruling, no run"* ]] || { echo "refused: row R76's first cell does not say ruling, no run" >&2; finish 95; }
NAMED=("$OPT" "--captain-ruling R76")
if [ "$LIFT" = yes ]; then NAMED+=("--den-ruling R76"); fi
for named in "${NAMED[@]}"; do
  if ! names "$first" "$named"; then
    echo "refused: row R76's first cell does not name $BT$named$BT as a whole backticked string" >&2
    finish 95
  fi
done
[[ $first == *"$BT$OPT$BT ($FILE)"* ]] \
  || { echo "refused: row R76's first cell does not pair $BT$OPT$BT with ($FILE)" >&2; finish 95; }
[ -f "$WT" ] && [ -f "$DEC" ] \
  || { echo "refused: the data $WT or the decision file $DEC is not in this checkout" >&2; finish 92; }
[ -d "$RUN_DIR" ] && pinned "$RUN_DIR/factors/transfer_mab21l2_all.csv" "$TABLE_M_SHA" \
  && pinned "$RUN_DIR/factors/gata4_pooled.csv" "$TABLE_G_SHA" \
  || { echo "refused: the gata4 run folder $RUN_DIR does not hold the two tables at their pins" >&2; finish 87; }
venv=0
for pth in .venv/lib/python3*/site-packages/vec.pth; do
  if [ -f "$pth" ] && [ "$(readlink -f "$(< "$pth")")" = "$(readlink -f src)" ]; then venv=1; fi
done
[ "$venv" = 1 ] || {
  echo "refused: this checkout has no .venv of its own (uv sync --frozen here) whose vec.pth names its src" >&2
  finish 88
}
O=runs/t3r76_${head:0:7}_$RUN
mkdir -p runs
mkdir "$O" 2> /dev/null || { echo "refused: the output folder $O exists" >&2; finish 93; }
echo "$head" > "$S/head"
step_ok

STEP=tests
begin
$PY -m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp="$S/tmp/pytest" \
  tests/test_t3_r76_build.py tests/test_t3_r76_old_new.py tests/test_t3_r76_launcher.py \
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
build() {  # folder name args: the builder, then `vec check`; 11 build, 12 refused, 16 stop, 17 factor, 13 14 15
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
identity() {  # STEP name recorded_sha256 code args: one identity rebuild into gate/, built once
  STEP=$1
  local name=$2 ref=$3 code=$4 got
  shift 4
  begin
  build gate "$name" "$@" || finish $?
  got=$(sha256sum "$O/gate/$name" | cut -d' ' -f1) || got=""
  [ "$got" = "$ref" ] || { echo "stopped: gate/$name is not its recorded bytes" >&2; finish "$code"; }
  step_ok
}
: > "$S/stamps"

identity c2 c2-s2.h5ad "$C2S2_SHA" 21 --form c2-s2
identity c2g2 T3_c2_s2_g2.h5ad "$C2G2_SHA" 22 $C2G2_ARGS
identity s067 T3_c2_s2_s067.h5ad "$S067_SHA" 23 $S067_ARGS
identity gonly T3_c2_s2_gonly.h5ad "$GONLY_SHA" 24 $GONLY_ARGS

for STEP in a b; do begin; build "$STEP" "$FILE.h5ad" $ARGS || finish $?; step_ok; done

STEP=done
begin
(cd "$O" && sha256sum gate/*.h5ad a/*.h5ad b/*.h5ad) > "$S/SHA256SUMS" || finish 40
cp "$S/SHA256SUMS" "$O/SHA256SUMS" || finish 40
cmp -s "$O/a/$FILE.h5ad" "$O/b/$FILE.h5ad" || { echo "stopped: the two builds of $FILE.h5ad differ" >&2; finish 41; }
finish 0
