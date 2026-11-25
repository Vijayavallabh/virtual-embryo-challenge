#!/usr/bin/env bash
# Row R72's sign cube files (the captain's order of 2026-11-25; T3_MOON5_design.md with fixes F1 to F10): C2's
# command (build_two_knockout.py --form c2-s2, recorded sha256 47111ac1d07e...) with its factor replaced by a vertex of
# the cube, `experiments/t3_decomp/build_sign_cube.py build ... --form c2-s2 --vertex <v> --cube-tables MANIFEST
# --captain-ruling R72`, for exactly the vertices that G0 (a to g) and row K31's fixed outcome leave, as the committed
# reading in experiments/t3_decomp/results/k31_cube_read/ names them (build_sign_cube.py vertices). Test-board builds on
# a DGX campaign checkout (never the reference checkout), on the machine that built C2 (its identity gate). Run from
# the checkout's root, detached (setsid nohup ... &), and wait on the exit files with a background poller; never a
# foreground sleep:
#   bash experiments/t3_decomp/t3_cube_build.sh <COMMIT> <RUN>
# RUN names this run's scratch folder ~/v/scratch/t3cube/<RUN> and its output folder runs/t3_cube_<commit7>_<RUN>.
# Neither may exist: a run folder is never reused and nothing is ever deleted here, so a second run takes a new name.
# In order, each step first writing the word running to ~/v/scratch/t3cube/<RUN>/<STEP>.exit, then its exit code; the
# code that ends the run also goes to all.exit (running until then):
#   gate      HEAD is COMMIT (the full hash) and the tree is clean (git status --porcelain prints nothing); TRIES.md
#             holds exactly one row R72 whose first cell opens "R72 (the captain's order", says "ruling, no run" and
#             holds `--captain-ruling R72` and each of the seven `--form c2-s2 --vertex <v>` as one whole
#             backtick-delimited string; exactly one row K31 naming its script, whose recorded cube tables manifest
#             sha256 is the committed MANIFEST.json's; row K31's reading.json and G0's g0.json are committed there; the
#             output folder does not exist
#   tests     the builder's, the cube's and this launcher's targeted tests (synthetic data), 8 BLAS threads
#   vertices  the vertices the committed reading builds, checked against G0's reading (build_sign_cube.py vertices):
#             none (STOP or VOID) finishes 30 with nothing built; a reading that does not check finishes 31
#   ppp       the identity gate: --vertex ppp built once into gate/, as a build step below, and its sha256 C2's
#             recorded bytes (t3two_fb01487.SHA256SUMS); never a file to upload
#   <v>a <v>b each vertex built twice, into folders a and b
#   pair      each vertex's two builds byte for byte (one sha256 per file)
#   done      one SHA256SUMS over every built file, copied into the output folder
# A build step is the builder's command, then `vec check --board T3:bcatenin`, then its stamp: it fails unless the
# stamp says COMMIT and clean True. Builds run at 16 BLAS threads, as the earlier T3 builds; NUMPY_MADVISE_HUGEPAGE=0.
# Every build and `vec check` output goes to /dev/null: a test build is read by its exit code, its sha256 and its
# stamp's commit and clean fields only, never its log, size, cell count or sidecar record. The run folder of the
# factor tables is fixed below (read only, no override). Nothing here deletes, fetches or uploads anything.
# Exit codes: 89 not two arguments, 90 odd RUN, 91 run folder exists, 96 HEAD is not COMMIT, 98 tree not clean, 95 row
# R72, 94 row K31 or the manifest, 93 the readings are not committed, 92 output folder exists, 10 tests, 30 nothing to
# build, 31 the reading does not check against G0's, 11 build failed, 12 the builder refused (its exit 2), 13 vec
# check, 14 no stamp, 15 stamp not this commit and clean, 16 the vertex cannot be made (the builder's exit 6), 17 the
# factor is not the one checked (the builder's exit 5), 20 the identity gate is not C2's bytes, 21 a pair differs, 39
# hash list, 99 the launcher stopped unexpectedly (a failing command, a TERM, INT or HUP signal, or an exit that was
# not one of the finishes above).
set -u
set -o pipefail
[ $# -eq 2 ] && [ -n "$1" ] && [ -n "$2" ] || { echo "usage: t3_cube_build.sh COMMIT RUN" >&2; exit 89; }
COMMIT=$1
RUN=$2
case "$RUN" in
  *[!A-Za-z0-9_-]*) echo "refused: RUN takes letters, digits, _ and - only" >&2; exit 90 ;;
esac
export PATH=$HOME/.local/bin:$PATH
S=$HOME/v/scratch/t3cube/$RUN
mkdir -p "$HOME/v/scratch/t3cube"
mkdir "$S" 2> /dev/null || { echo "refused: the run folder $S exists (a run folder is never reused)" >&2; exit 91; }
mkdir "$S/tmp"
export NUMPY_MADVISE_HUGEPAGE=0 TMPDIR="$S/tmp"
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
PY="uv run --frozen --no-sync python -I -B"
BOARD=T3:bcatenin
DEC=experiments/t3_decomp/results/gata4/gata4_decisions.json
RUN_DIR=$HOME/v/virtual-embryo-challenge/runs/t3_decomp/gata4_20261028T122425Z  # read only; fixed, no override
WT=data/raw/T2/E8.75.h5ad
MANIFEST=experiments/t3_decomp/results/t3_cube/MANIFEST.json
READ_OUT=experiments/t3_decomp/results/k31_cube_read
BUILDER=experiments/t3_decomp/build_sign_cube.py
C2_SHA=47111ac1d07e73c0c26e990a474dff0de94ae22d6a4474a569bcb674ae2b85a7  # c2-s2 (C2), t3two_fb01487.SHA256SUMS
VERTICES="ppn pnp pnn npp npn nnp nnn"
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
echo running > "$S/all.exit"
trap abort ERR TERM INT HUP EXIT
begin

head=$(git rev-parse HEAD 2> /dev/null) || head=""
[ "$head" = "$COMMIT" ] || { echo "refused: HEAD is not $COMMIT (give the full hash)" >&2; finish 96; }
[ -z "$(git status --porcelain 2> /dev/null)" ] || { echo "refused: the tree is not clean" >&2; finish 98; }
[ "$(grep -c '^| R72 ' experiments/harness/TRIES.md)" = 1 ] \
  || { echo "refused: TRIES.md does not hold exactly one row R72" >&2; finish 95; }
first=$(grep '^| R72 ' experiments/harness/TRIES.md | cut -d'|' -f2)
[[ $first == " R72 (the captain's order"* ]] \
  || { echo "refused: row R72's first cell does not open: R72 (the captain's order" >&2; finish 95; }
[[ $first == *"ruling, no run"* ]] || { echo "refused: row R72's first cell does not say: ruling, no run" >&2; finish 95; }
names "$first" "--captain-ruling R72" || { echo "refused: row R72 does not name --captain-ruling R72" >&2; finish 95; }
for v in $VERTICES; do
  names "$first" "--form c2-s2 --vertex $v" || { echo "refused: row R72 does not name vertex $v" >&2; finish 95; }
done
[ "$(grep -c '^| K31 .*`experiments/t3_decomp/k31_cube_read.py`' experiments/harness/TRIES.md)" = 1 ] \
  || { echo "refused: TRIES.md does not hold exactly one row K31 naming its script" >&2; finish 94; }
k31_sha=$($PY experiments/t3_decomp/k31_cube_read.py row 2> /dev/null) || k31_sha=""
man_sha=$(sha256sum "$MANIFEST" 2> /dev/null | cut -d' ' -f1) || man_sha=""
[ -n "$k31_sha" ] && [ "$k31_sha" = "$man_sha" ] \
  || { echo "refused: row K31's manifest sha256 is not the committed manifest's" >&2; finish 94; }
for f in "$MANIFEST" "$READ_OUT/reading.json" "$READ_OUT/g0.json"; do
  GIT_OPTIONAL_LOCKS=0 git ls-files --error-unmatch "$f" > /dev/null 2>&1 \
    || { echo "refused: $f is not committed" >&2; finish 93; }
done
O=runs/t3_cube_${head:0:7}_$RUN
mkdir -p runs
mkdir "$O" 2> /dev/null || { echo "refused: the output folder $O exists" >&2; finish 92; }
echo "$head" > "$S/head"
step_ok

STEP=tests
begin
$PY -m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp="$S/tmp/pytest" \
  tests/test_t3_cube_build.py tests/test_t3_cube_build_launcher.py tests/test_sign_cube.py tests/test_cube_tables.py \
  tests/test_k31_cube_read.py tests/test_t3_two_knockout_build.py experiments/t3_decomp/test_build_submission.py \
  > "$S/tests.out" 2>&1 || finish 10
step_ok

STEP=vertices
begin
rc=0
BUILD=$($PY "$BUILDER" vertices --reading "$READ_OUT/reading.json" --g0 "$READ_OUT/g0.json" 2> "$S/vertices.err") \
  || rc=$?
if [ $rc = 3 ]; then echo "row K31's reading builds nothing" > "$S/vertices"; finish 30; fi
if [ $rc != 0 ]; then finish 31; fi
for v in $BUILD; do
  case " $VERTICES " in *" $v "*) ;; *) finish 31 ;; esac
done
echo "$BUILD" > "$S/vertices"
step_ok

export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16
stamp() {  # the stamp's commit and clean fields only
  $PY -c 'import json,sys; s=json.load(open(sys.argv[1])); print(s["commit"], s["clean"])' "$1.stamp.json"
}
build() {  # folder vertex: the builder, then `vec check`; 11 build, 12 refused, 16 not made, 17 factor
  local d="$O/$1" v=$2 name rc=0 st
  name=T3_cube_$v.h5ad
  mkdir -p "$d"
  $PY "$BUILDER" build "$DEC" --form c2-s2 --vertex "$v" --wt "$WT" --run-dir "$RUN_DIR" --board "$BOARD" \
    --seed 1000001 --cube-tables "$MANIFEST" --captain-ruling R72 --out "$d/$name" > /dev/null 2>&1 || rc=$?
  if [ $rc = 2 ]; then return 12; fi
  if [ $rc = 6 ]; then return 16; fi
  if [ $rc = 5 ]; then return 17; fi
  if [ $rc != 0 ]; then return 11; fi
  uv run --frozen --no-sync vec check "$d/$name" --board "$BOARD" > /dev/null 2>&1 || return 13
  st=$(stamp "$d/$name" 2> /dev/null) || return 14
  echo "$(basename "$d")/$name $st" >> "$S/stamps"
  [ "$st" = "$COMMIT True" ] || return 15
}
: > "$S/stamps"

STEP=ppp
begin
build gate ppp || finish $?
got=$(sha256sum "$O/gate/T3_cube_ppp.h5ad" | cut -d' ' -f1) || got=""
[ "$got" = "$C2_SHA" ] || { echo "stopped: the identity gate is not C2's recorded bytes" >&2; finish 20; }
step_ok

for v in $BUILD; do
  for side in a b; do
    STEP=$v$side
    begin
    build "$side" "$v" || finish $?
    step_ok
  done
done

STEP=pair
begin
for v in $BUILD; do
  cmp -s "$O/a/T3_cube_$v.h5ad" "$O/b/T3_cube_$v.h5ad" \
    || { echo "stopped: the two builds of T3_cube_$v.h5ad differ" >&2; finish 21; }
done
step_ok

STEP=done
begin
(cd "$O" && sha256sum gate/*.h5ad a/*.h5ad b/*.h5ad) > "$S/SHA256SUMS" || finish 39
cp "$S/SHA256SUMS" "$O/SHA256SUMS" || finish 39
finish 0
