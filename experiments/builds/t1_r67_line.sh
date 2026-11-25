#!/usr/bin/env bash
# Row R67's three files (the captain's order of 2026-11-25), each T1_s2rcMSB's command (t1 --g-mom --base lib3
# --feng-groups recomputed --feng-fraction-mult 1.5 --direction-momentum all, row R26's) with one change, MS's own
# damping (--direction-momentum-damping) in place of harness.B0_DAMPING, each only with --captain-ruling R67: row R58's
# route (T1_s2rcMSBm20, the same command at 2.0) over three literal values. One file per call, given by its name:
#   T1_s2rcMSBm125  --direction-momentum-damping 1.25
#   T1_s2rcMSBm15   --direction-momentum-damping 1.5
#   T1_s2rcMSBm075  --direction-momentum-damping 0.75
# nothing else changed (default seed). Test-board builds on a DGX campaign worktree that holds the relabelling's output
# (data/derived/feng_relabel_c) and its own venv, made there by `uv sync --frozen` (never the reference checkout). Run
# from the root of that worktree, detached, and wait on the exit files with a background poller; never a foreground
# sleep:
#   bash experiments/builds/t1_r67_line.sh <COMMIT> <RUN> <FILE_NAME>
# FILE_NAME is one of the three in the table below (tests/test_t1_r67_launcher.py holds it equal to the builder's
# T1_R67_TABLE, the one table row R67's strings are checked against); any other name, or a fourth argument, is refused
# before any folder is made. RUN names this run's own scratch folder ~/v/scratch/t1r67/<RUN> and its output folder
# runs/t1r67_<commit7>_<RUN>. Neither may exist: a run folder is never reused and nothing is ever deleted here, so a
# second run, or a second file, takes a new name. STEM below is FILE_NAME without T1_s2rcMSB (m125, m15 or m075).
# In order, each step first writing the word running to ~/v/scratch/t1r67/<RUN>/<STEP>.exit, then its exit code:
#   gate   HEAD is COMMIT (the full hash), the tracked tree is clean, docs/S2_RULES.sha256 checks, TRIES.md holds exactly
#          one row R67 whose first cell starts "R67 (the captain's ", says "ruling, no run" and names FILE_NAME's option
#          string and `--captain-ruling R67`, each as a whole backticked string, the option string followed there by
#          " (FILE_NAME)" as the row pairs each string with its file, the relabelling's output folder
#          data/derived/feng_relabel_c exists, and .venv's vec.pth names this checkout's src; then it writes the head
#          file
#   tests  the builders', checkers' and launchers' targeted tests (synthetic stages only), 8 BLAS threads
#   MSB m10 m20  three identity rebuilds, each once, each of which must give its known bytes or everything stops:
#          MSB    T1_s2rcMSB (row R26's command, the base of row R67's files); MSB_SHA (t1msc_a3e879a_m1,
#                 t1r51_0cdcc6f_r1 and t1lmm20_89700ef_lm1.SHA256SUMS)
#          m10    T1_s2rcMSBm10 (row R51's file, MS's damping 1.0); M10_SHA (t1r51_0cdcc6f_r1.SHA256SUMS)
#          m20    T1_s2rcMSBm20 (row R58's file, MS's damping 2.0); M20_SHA (t1lmm20_89700ef_lm1.SHA256SUMS)
#   STEMbase  FILE_NAME's command once through experiments/builds/t1_r67_check.py base, on its own path with MS's
#          damping put back at B0's: it must give MSB_SHA, T1_s2rcMSB's bytes, or everything stops
#   STEMa STEMb  FILE_NAME by its command, twice, into separate folders
#   done   one sha256 list over every built file, the two builds of FILE_NAME compared byte for byte
# A build step also fails (15) unless its stamp says this commit and clean True.
# Builds run at 16 BLAS threads, as the reference builds (T1 files move with the thread count), seed and board as the
# reference commands (the builder's defaults: seed 1000001, board T1:test). Every builder, checker build and `vec check`
# output goes to /dev/null: a test build is read by its exit code, its sha256 and its stamp's commit and clean fields
# (AGENTS.md), never its log, size, record or cell count. Nothing here uploads anything.
# Exit codes: 90 odd RUN, 89 FILE_NAME not one of row R67's three (or a fourth argument), 91 run folder exists,
# 96 HEAD is not COMMIT, 98 tracked tree not clean, 97 S2_RULES.sha256, 95 row R67 not as above, 92 no relabelling
# output, 88 .venv not this checkout's, 93 output folder exists, 10 tests, 11 build, 12 checker build, 13 vec check,
# 14 no stamp, 15 stamp not this commit and clean, 23 MSB is not T1_s2rcMSB's bytes, 24 m10 is not T1_s2rcMSBm10's
# bytes, 25 m20 is not T1_s2rcMSBm20's bytes, 51 STEMbase is not MSB's bytes, 40 hash list, 41 STEMa and STEMb differ.
set -u
COMMIT=${1:?usage: t1_r67_line.sh COMMIT RUN FILE_NAME}
RUN=${2:?usage: t1_r67_line.sh COMMIT RUN FILE_NAME}
FILE=${3:?usage: t1_r67_line.sh COMMIT RUN FILE_NAME}
[ "$#" = 3 ] || { echo "refused: give COMMIT RUN FILE_NAME and nothing else" >&2; exit 89; }
case "$RUN" in
  *[!A-Za-z0-9_-]*) echo "refused: RUN takes letters, digits, _ and - only" >&2; exit 90 ;;
esac
MSB_OPT="--feng-fraction-mult 1.5 --direction-momentum all"
# FILE_NAME -> its option string, as row R67 pairs them (the builder's T1_R67_TABLE, in the row's order)
case "$FILE" in
  T1_s2rcMSBm125) OPT="$MSB_OPT --direction-momentum-damping 1.25" ;;
  T1_s2rcMSBm15) OPT="$MSB_OPT --direction-momentum-damping 1.5" ;;
  T1_s2rcMSBm075) OPT="$MSB_OPT --direction-momentum-damping 0.75" ;;
  *) echo "refused: FILE_NAME is not one of row R67's three files" >&2; exit 89 ;;
esac
STEM=${FILE#T1_s2rcMSB}
ARGS="--g-mom --base lib3 --feng-groups recomputed $OPT --captain-ruling R67"
export PATH=$HOME/.local/bin:$PATH
unset UV_PROJECT_ENVIRONMENT  # uv runs this checkout's own .venv, the one the gate checks
S=$HOME/v/scratch/t1r67/$RUN
mkdir -p "$HOME/v/scratch/t1r67"
mkdir "$S" 2> /dev/null || { echo "refused: the run folder $S exists (a run folder is never reused)" >&2; exit 91; }
mkdir "$S/tmp"
export NUMPY_MADVISE_HUGEPAGE=0 TMPDIR="$S/tmp"
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
PY="uv run --frozen --no-sync python -I -B"
BOARD=T1:test
MSB_SHA=78bb17eea5721e590cfd2c9f707e9cd6f2f496e081acede4f697dc4035871da2  # T1_s2rcMSB, t1msc_a3e879a_m1.SHA256SUMS
M10_SHA=a0c46efac87deb3756b3ca489a215694640707ea8a8053e254b8a86ee9d7d964  # T1_s2rcMSBm10, t1r51_0cdcc6f_r1
M20_SHA=d82c3db2e5f43a3abeefe32a04e4d58d2da41cefa8fa979eb642de884c435401  # T1_s2rcMSBm20, t1lmm20_89700ef_lm1
MSB_ARGS="--g-mom --base lib3 --feng-groups recomputed --feng-fraction-mult 1.5 --direction-momentum all --captain-ruling R26"
M10_ARGS="--g-mom --base lib3 --feng-groups recomputed --feng-fraction-mult 1.5 --direction-momentum all --direction-momentum-damping 1.0 --captain-ruling R51"
M20_ARGS="--g-mom --base lib3 --feng-groups recomputed --feng-fraction-mult 1.5 --direction-momentum all --direction-momentum-damping 2.0 --captain-ruling R58"
STEP=gate
finish() { echo "$1" > "$S/$STEP.exit"; exit "$1"; }
begin() { echo running > "$S/$STEP.exit"; }
step_ok() { echo 0 > "$S/$STEP.exit"; }
begin

head=$(git rev-parse HEAD)
[ "$head" = "$COMMIT" ] || { echo "refused: HEAD is not $COMMIT (give the full hash)" >&2; finish 96; }
if [ -n "$(git status --porcelain --untracked-files=no)" ]; then
  echo "refused: tracked files differ from the commit (the build stamps must say clean)" >&2
  finish 98
fi
sha256sum -c docs/S2_RULES.sha256 > /dev/null 2>&1 || finish 97
first_cell() {  # ID CODE: row ID's first cell, the row held exactly once and the cell the captain's ruling, else CODE
  local id=$1 code=$2 first
  [ "$(grep -c "^| $id " experiments/harness/TRIES.md)" = 1 ] \
    || { echo "refused: TRIES.md does not hold exactly one row $id" >&2; finish "$code"; }
  first=$(grep "^| $id " experiments/harness/TRIES.md | cut -d'|' -f2)
  printf '%s\n' "$first" | grep -q "^ $id (the captain's " \
    || { echo "refused: row $id's first cell is not the captain's ruling" >&2; finish "$code"; }
  printf '%s\n' "$first" | grep -qF -- "ruling, no run" \
    || { echo "refused: row $id's first cell does not say ruling, no run" >&2; finish "$code"; }
  printf '%s\n' "$first"
}
r67=$(first_cell R67 95) || exit $?
for opt in "$OPT" "--captain-ruling R67"; do
  printf '%s\n' "$r67" | grep -o '`[^`]*`' | grep -qxF -- "\`$opt\`" \
    || { echo "refused: row R67's first cell does not name \`$opt\` as a whole backticked string" >&2; finish 95; }
done
printf '%s\n' "$r67" | grep -qF -- "\`$OPT\` ($FILE)" \
  || { echo "refused: row R67's first cell does not pair \`$OPT\` with ($FILE)" >&2; finish 95; }
[ -d data/derived/feng_relabel_c ] \
  || { echo "refused: the relabelling's output data/derived/feng_relabel_c is not in this checkout" >&2; finish 92; }
venv=0
for pth in .venv/lib/python3*/site-packages/vec.pth; do
  if [ -f "$pth" ] && [ "$(readlink -f "$(< "$pth")")" = "$(readlink -f src)" ]; then venv=1; fi
done
[ "$venv" = 1 ] || {
  echo "refused: this checkout has no .venv of its own (uv sync --frozen here) whose vec.pth names its src" >&2
  finish 88
}
O=runs/t1r67_${head:0:7}_$RUN
mkdir -p runs
mkdir "$O" 2> /dev/null || { echo "refused: the output folder $O exists" >&2; finish 93; }
echo "$head" > "$S/head"
step_ok

STEP=tests
begin
$PY -m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp="$S/tmp/pytest" \
  tests/test_t1_r67_launcher.py tests/test_t1_r67_build.py tests/test_t1_lm_m20_build.py \
  tests/test_t1_lm_m20_launchers.py tests/test_t1_r51_build.py tests/test_t1_msb_build.py tests/test_t1_ms_build.py \
  tests/test_t1_candidates_build.py tests/test_hidden_counts.py tests/test_s2_rules.py \
  > "$S/tests.out" 2>&1 || finish 10
step_ok

export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16
stamp() {  # the stamp's commit and clean fields only
  $PY -c 'import json,sys; s=json.load(open(sys.argv[1])); print(s["commit"], s["clean"])' "$1.stamp.json"
}
checked() {  # folder name: `vec check`, then the stamp; exit 13 check, 14 no stamp, 15 stamp not this commit and clean
  local d=$1 name=$2 st
  uv run --frozen --no-sync vec check "$d/$name.h5ad" --board "$BOARD" > /dev/null 2>&1 || return 13
  st=$(stamp "$d/$name.h5ad" 2> /dev/null) || return 14
  echo "$(basename "$d")/$name $st" >> "$S/stamps"
  [ "$st" = "$head True" ] || return 15
}
BUILT=()  # every built file, folder/name.h5ad under $O, in build order
build() {  # folder name args: the build, then checked; exit 11 build
  local d="$O/$1" name=$2
  BUILT+=("$1/$2.h5ad")
  shift 2
  mkdir -p "$d"
  $PY experiments/builds/build_test.py t1 "$@" --out "$d/$name.h5ad" > /dev/null 2>&1 || return 11
  checked "$d" "$name"
}
build_base() {  # folder name args: the build on its own path with MS's damping at B0's; exit 12
  local d="$O/$1" name=$2
  BUILT+=("$1/$2.h5ad")
  shift 2
  mkdir -p "$d"
  $PY experiments/builds/t1_r67_check.py base t1 "$@" --out "$d/$name.h5ad" > /dev/null 2>&1 || return 12
  checked "$d" "$name"
}
same() {  # folder name sha code: stop with the code unless the built file has the known bytes
  local got
  got=$(sha256sum "$O/$1/$2.h5ad" | cut -d' ' -f1)
  [ "$got" = "$3" ] || { echo "stopped: step $1 does not give the known bytes of $2" >&2; finish "$4"; }
}
: > "$S/stamps"
identity() {  # STEP name sha code args: one identity rebuild
  STEP=$1
  local name=$2 sha=$3 code=$4
  shift 4
  begin
  build "$STEP" "$name" "$@" || finish $?
  same "$STEP" "$name" "$sha" "$code"
  step_ok
}
identity MSB T1_s2rcMSB "$MSB_SHA" 23 $MSB_ARGS
identity m10 T1_s2rcMSBm10 "$M10_SHA" 24 $M10_ARGS
identity m20 T1_s2rcMSBm20 "$M20_SHA" 25 $M20_ARGS

STEP=${STEM}base
begin
build_base "$STEP" "${FILE}_base" $ARGS || finish $?
same "$STEP" "${FILE}_base" "$MSB_SHA" 51
step_ok

for STEP in "${STEM}a" "${STEM}b"; do begin; build "$STEP" "$FILE" $ARGS || finish $?; step_ok; done

STEP=done
begin
(cd "$O" && sha256sum "${BUILT[@]}") > "$S/SHA256SUMS" || finish 40
cp "$S/SHA256SUMS" "$O/SHA256SUMS" || finish 40
cmp -s "$O/${STEM}a/$FILE.h5ad" "$O/${STEM}b/$FILE.h5ad" || finish 41
step_ok
