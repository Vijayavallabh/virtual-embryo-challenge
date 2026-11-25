#!/usr/bin/env bash
# Row R66's thirty files (T1's third-round step-on slots, the captain's order of 2026-11-25), each on T1_s2rcMS20's
# command (t1 --g-mom --base lib3 --feng-groups recomputed --feng-fraction-mult 2.0 --direction-momentum all, row
# R51's file) with two or three of line 12's changes (slots K1, K2 and K3), g_mom's damping 2.5 (slot C1) or Blood's
# step at 2.0 (slot C3), each only with --captain-ruling R66. One file per call, given by its name, T1_s2rcMS20 plus
# the suffixes in d, H, m order:
#   d10 d15 d20 d25  --g-mom-damping 1.0 | 1.5 | 2.0 | 2.5 (g_mom's within-type damping; 2.5 only alone, C1)
#   H05 H10 H20      --held-momentum 0.5 | 1.0 | 2.0 (Blood's own released change after MS's step and before g_mom: R1
#                    reversed for Blood and T64's close lifted in each file with it, as in row R60; 2.0 only alone, C3)
#   m10 m20          --direction-momentum-damping 1.0 | 2.0 (MS's own damping)
# nothing else changed (default seed). Test-board builds on a DGX campaign worktree that holds the relabelling's output
# (data/derived/feng_relabel_c) and its own venv, made there by `uv sync --frozen` (never the reference checkout). Run
# from the root of that worktree, detached, and wait on the exit files with a background poller; never a foreground
# sleep:
#   bash experiments/builds/t1_r66_line.sh <COMMIT> <RUN> <FILE_NAME>
# FILE_NAME is one of the thirty in the table below (tests/test_t1_r66_launcher.py holds it equal to the builder's
# T1_R66_FILES, the one table row R66's strings are checked against); any other name, or a fourth argument, is refused
# before any folder is made. RUN names this run's own scratch folder ~/v/scratch/t1r66/<RUN> and its output folder
# runs/t1r66_<commit7>_<RUN>. Neither may exist: a run folder is never reused and nothing is ever deleted here, so a
# second run, or a second file, takes a new name. STEM below is FILE_NAME without T1_s2rcMS20 (e.g. d15H05m20).
# In order, each step first writing the word running to ~/v/scratch/t1r66/<RUN>/<STEP>.exit, then its exit code:
#   gate   HEAD is COMMIT (the full hash), the tracked tree is clean, docs/S2_RULES.sha256 checks, TRIES.md holds exactly
#          one row R66 whose first cell starts "R66 (the captain's ", says "ruling, no run" and names FILE_NAME's option
#          string and `--captain-ruling R66`, each as a whole backticked string, the option string followed there by
#          " (FILE_NAME)" as the row pairs each string with its file, the relabelling's output folder
#          data/derived/feng_relabel_c exists, and .venv's vec.pth names this checkout's src; then it writes the head
#          file
#   tests  the builders', checkers' and launchers' targeted tests (synthetic stages only), 8 BLAS threads
#   MS20   T1_s2rcMS20 (row R51's file, the base of row R66's), once: it must give MS20_SHA (t1r51_0cdcc6f_r1 and
#          t1l12_128dfa1_b1.SHA256SUMS), or everything stops
#   d10 d15 d20 H05 H10 m10 m20  the line-12 files (row R60's) whose step FILE_NAME uses, each once, in that order:
#          on each axis FILE_NAME changes, the line-12 file at FILE_NAME's value, and where line 12 has no file at that
#          value (C1's d25, C3's H20) every line-12 file on that axis; each must give the bytes
#          t1l12_128dfa1_b1.SHA256SUMS holds for it, or everything stops
#   STEMbase  FILE_NAME's command once through experiments/builds/t1_r66_check.py base, on its own path with every one
#          of its steps at T1_s2rcMS20's value (g_mom's damping at 0.5, Blood's factors at 1, MS's damping at the
#          default): it must give MS20_SHA, T1_s2rcMS20's bytes, or everything stops
#   STEMa STEMb  FILE_NAME by its command, twice, into separate folders
#   done   one sha256 list over every built file, the two builds of FILE_NAME compared byte for byte
# A build step also fails (15) unless its stamp says this commit and clean True.
# Builds run at 16 BLAS threads, as the reference builds (T1 files move with the thread count), seed and board as the
# reference commands (the builder's defaults: seed 1000001, board T1:test). Every builder, checker build and `vec check`
# output goes to /dev/null: a test build is read by its exit code, its sha256 and its stamp's commit and clean fields
# (AGENTS.md), never its log, size, record or cell count. Nothing here uploads anything.
# Exit codes: 90 odd RUN, 89 FILE_NAME not one of row R66's thirty (or a fourth argument), 91 run folder exists,
# 96 HEAD is not COMMIT, 98 tracked tree not clean, 97 S2_RULES.sha256, 95 row R66 not as above, 92 no relabelling
# output, 88 .venv not this checkout's, 93 output folder exists, 10 tests, 11 build, 12 checker build, 13 vec check,
# 14 no stamp, 15 stamp not this commit and clean, 39 MS20 is not T1_s2rcMS20's bytes, 71 d10 is not T1_s2rcMS20d10's
# bytes, 72 d15 is not T1_s2rcMS20d15's, 73 d20 is not T1_s2rcMS20d20's, 74 H05 is not T1_s2rcMS20H05's, 75 H10 is not
# T1_s2rcMS20H10's, 76 m10 is not T1_s2rcMS20m10's, 77 m20 is not T1_s2rcMS20m20's, 51 STEMbase is not MS20's bytes,
# 40 hash list, 41 STEMa and STEMb differ.
set -u
COMMIT=${1:?usage: t1_r66_line.sh COMMIT RUN FILE_NAME}
RUN=${2:?usage: t1_r66_line.sh COMMIT RUN FILE_NAME}
FILE=${3:?usage: t1_r66_line.sh COMMIT RUN FILE_NAME}
[ "$#" = 3 ] || { echo "refused: give COMMIT RUN FILE_NAME and nothing else" >&2; exit 89; }
case "$RUN" in
  *[!A-Za-z0-9_-]*) echo "refused: RUN takes letters, digits, _ and - only" >&2; exit 90 ;;
esac
L12="--g-mom --base lib3 --feng-groups recomputed --feng-fraction-mult 2.0 --direction-momentum all"
MS20_OPT="--feng-fraction-mult 2.0 --direction-momentum all"
# FILE_NAME -> its option string, as row R66 pairs them (the builder's T1_R66_TABLE, in the row's order)
case "$FILE" in
  T1_s2rcMS20d10H05) OPT="$MS20_OPT --g-mom-damping 1.0 --held-momentum 0.5" ;;
  T1_s2rcMS20d10H10) OPT="$MS20_OPT --g-mom-damping 1.0 --held-momentum 1.0" ;;
  T1_s2rcMS20d15H05) OPT="$MS20_OPT --g-mom-damping 1.5 --held-momentum 0.5" ;;
  T1_s2rcMS20d15H10) OPT="$MS20_OPT --g-mom-damping 1.5 --held-momentum 1.0" ;;
  T1_s2rcMS20d20H05) OPT="$MS20_OPT --g-mom-damping 2.0 --held-momentum 0.5" ;;
  T1_s2rcMS20d20H10) OPT="$MS20_OPT --g-mom-damping 2.0 --held-momentum 1.0" ;;
  T1_s2rcMS20d10m10) OPT="$MS20_OPT --g-mom-damping 1.0 --direction-momentum-damping 1.0" ;;
  T1_s2rcMS20d10m20) OPT="$MS20_OPT --g-mom-damping 1.0 --direction-momentum-damping 2.0" ;;
  T1_s2rcMS20d15m10) OPT="$MS20_OPT --g-mom-damping 1.5 --direction-momentum-damping 1.0" ;;
  T1_s2rcMS20d15m20) OPT="$MS20_OPT --g-mom-damping 1.5 --direction-momentum-damping 2.0" ;;
  T1_s2rcMS20d20m10) OPT="$MS20_OPT --g-mom-damping 2.0 --direction-momentum-damping 1.0" ;;
  T1_s2rcMS20d20m20) OPT="$MS20_OPT --g-mom-damping 2.0 --direction-momentum-damping 2.0" ;;
  T1_s2rcMS20H05m10) OPT="$MS20_OPT --held-momentum 0.5 --direction-momentum-damping 1.0" ;;
  T1_s2rcMS20H05m20) OPT="$MS20_OPT --held-momentum 0.5 --direction-momentum-damping 2.0" ;;
  T1_s2rcMS20H10m10) OPT="$MS20_OPT --held-momentum 1.0 --direction-momentum-damping 1.0" ;;
  T1_s2rcMS20H10m20) OPT="$MS20_OPT --held-momentum 1.0 --direction-momentum-damping 2.0" ;;
  T1_s2rcMS20d10H05m10) OPT="$MS20_OPT --g-mom-damping 1.0 --held-momentum 0.5 --direction-momentum-damping 1.0" ;;
  T1_s2rcMS20d10H05m20) OPT="$MS20_OPT --g-mom-damping 1.0 --held-momentum 0.5 --direction-momentum-damping 2.0" ;;
  T1_s2rcMS20d10H10m10) OPT="$MS20_OPT --g-mom-damping 1.0 --held-momentum 1.0 --direction-momentum-damping 1.0" ;;
  T1_s2rcMS20d10H10m20) OPT="$MS20_OPT --g-mom-damping 1.0 --held-momentum 1.0 --direction-momentum-damping 2.0" ;;
  T1_s2rcMS20d15H05m10) OPT="$MS20_OPT --g-mom-damping 1.5 --held-momentum 0.5 --direction-momentum-damping 1.0" ;;
  T1_s2rcMS20d15H05m20) OPT="$MS20_OPT --g-mom-damping 1.5 --held-momentum 0.5 --direction-momentum-damping 2.0" ;;
  T1_s2rcMS20d15H10m10) OPT="$MS20_OPT --g-mom-damping 1.5 --held-momentum 1.0 --direction-momentum-damping 1.0" ;;
  T1_s2rcMS20d15H10m20) OPT="$MS20_OPT --g-mom-damping 1.5 --held-momentum 1.0 --direction-momentum-damping 2.0" ;;
  T1_s2rcMS20d20H05m10) OPT="$MS20_OPT --g-mom-damping 2.0 --held-momentum 0.5 --direction-momentum-damping 1.0" ;;
  T1_s2rcMS20d20H05m20) OPT="$MS20_OPT --g-mom-damping 2.0 --held-momentum 0.5 --direction-momentum-damping 2.0" ;;
  T1_s2rcMS20d20H10m10) OPT="$MS20_OPT --g-mom-damping 2.0 --held-momentum 1.0 --direction-momentum-damping 1.0" ;;
  T1_s2rcMS20d20H10m20) OPT="$MS20_OPT --g-mom-damping 2.0 --held-momentum 1.0 --direction-momentum-damping 2.0" ;;
  T1_s2rcMS20d25) OPT="$MS20_OPT --g-mom-damping 2.5" ;;
  T1_s2rcMS20H20) OPT="$MS20_OPT --held-momentum 2.0" ;;
  *) echo "refused: FILE_NAME is not one of row R66's thirty files" >&2; exit 89 ;;
esac
STEM=${FILE#T1_s2rcMS20}
[[ $STEM =~ ^(d(10|15|20|25))?(H(05|10|20))?(m(10|20))?$ ]] \
  || { echo "refused: FILE_NAME is not one of row R66's thirty files" >&2; exit 89; }
D=${BASH_REMATCH[2]}
H=${BASH_REMATCH[4]}
M=${BASH_REMATCH[6]}
IDS=""  # the line-12 files whose step FILE_NAME uses, in d, H, m order (header: d10 ... m20)
case "$D" in 10 | 15 | 20) IDS="$IDS d$D" ;; 25) IDS="$IDS d10 d15 d20" ;; esac
case "$H" in 05 | 10) IDS="$IDS H$H" ;; 20) IDS="$IDS H05 H10" ;; esac
case "$M" in 10 | 20) IDS="$IDS m$M" ;; esac
ARGS="--g-mom --base lib3 --feng-groups recomputed $OPT --captain-ruling R66"
export PATH=$HOME/.local/bin:$PATH
unset UV_PROJECT_ENVIRONMENT  # uv runs this checkout's own .venv, the one the gate checks
S=$HOME/v/scratch/t1r66/$RUN
mkdir -p "$HOME/v/scratch/t1r66"
mkdir "$S" 2> /dev/null || { echo "refused: the run folder $S exists (a run folder is never reused)" >&2; exit 91; }
mkdir "$S/tmp"
export NUMPY_MADVISE_HUGEPAGE=0 TMPDIR="$S/tmp"
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
PY="uv run --frozen --no-sync python -I -B"
BOARD=T1:test
MS20_SHA=629d3a1cffa8151f0345a169ba97420e4ee3b43e379bc3bf7d33ab947dcb62ef  # T1_s2rcMS20, t1r51_0cdcc6f_r1.SHA256SUMS
D10_SHA=16f35132bb96919f1dcc89ace1dd12f5a3ced88749e889086904e189ab9fd2a8  # T1_s2rcMS20d10, t1l12_128dfa1_b1
D15_SHA=cfea3b689696bd2e08abe04ad607ab5a826b300ad6d841c3f2d8fc475f536f63  # T1_s2rcMS20d15, t1l12_128dfa1_b1
D20_SHA=c5340b1505e4e032e764465d69973910b7071cedbd266487866606c0dafb8cd9  # T1_s2rcMS20d20, t1l12_128dfa1_b1
H05_SHA=97589938a95ef957d2aed1117f84147cd027f2a0c3d1c48fbcaf3c67570a87cc  # T1_s2rcMS20H05, t1l12_128dfa1_b1
H10_SHA=5888b0f4da0124e7561d17a002006d79b26efe48c14e346cf3e45faf1222be3e  # T1_s2rcMS20H10, t1l12_128dfa1_b1
M10_SHA=316b67b44df2501db6f5d43ea099bec107ba9422498ed90f1861fc32d70cc1c5  # T1_s2rcMS20m10, t1l12_128dfa1_b1
M20_SHA=961df1c7b4b8cce40006729f8439701313cb5d2b41b159763e20d1dbd5097000  # T1_s2rcMS20m20, t1l12_128dfa1_b1
MS20_ARGS="--g-mom --base lib3 --feng-groups recomputed --feng-fraction-mult 2.0 --direction-momentum all --captain-ruling R51"
D10_ARGS="$L12 --g-mom-damping 1.0 --captain-ruling R60"
D15_ARGS="$L12 --g-mom-damping 1.5 --captain-ruling R60"
D20_ARGS="$L12 --g-mom-damping 2.0 --captain-ruling R60"
H05_ARGS="$L12 --held-momentum 0.5 --captain-ruling R60"
H10_ARGS="$L12 --held-momentum 1.0 --captain-ruling R60"
M10_ARGS="$L12 --direction-momentum-damping 1.0 --captain-ruling R60"
M20_ARGS="$L12 --direction-momentum-damping 2.0 --captain-ruling R60"
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
r66=$(first_cell R66 95) || exit $?
for opt in "$OPT" "--captain-ruling R66"; do
  printf '%s\n' "$r66" | grep -o '`[^`]*`' | grep -qxF -- "\`$opt\`" \
    || { echo "refused: row R66's first cell does not name \`$opt\` as a whole backticked string" >&2; finish 95; }
done
printf '%s\n' "$r66" | grep -qF -- "\`$OPT\` ($FILE)" \
  || { echo "refused: row R66's first cell does not pair \`$OPT\` with ($FILE)" >&2; finish 95; }
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
O=runs/t1r66_${head:0:7}_$RUN
mkdir -p runs
mkdir "$O" 2> /dev/null || { echo "refused: the output folder $O exists" >&2; finish 93; }
echo "$head" > "$S/head"
step_ok

STEP=tests
begin
$PY -m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp="$S/tmp/pytest" \
  tests/test_t1_r66_launcher.py tests/test_t1_r66_build.py tests/test_t1_l12_launcher.py tests/test_t1_l12_build.py \
  tests/test_t1_r51_launcher.py tests/test_t1_r51_build.py tests/test_t1_hm_launcher.py tests/test_t1_hm_build.py \
  tests/test_t1_l67pb_launcher.py tests/test_t1_l67pb_build.py tests/test_t1_b2d_line_launcher.py \
  tests/test_t1_line5_build.py tests/test_t1_line3_build.py tests/test_t1_c_line_launcher.py \
  tests/test_t1_damping_line_build.py tests/test_t1_damping_line_launcher.py tests/test_t1_fraction_line_build.py \
  tests/test_t1_fraction_line_launcher.py tests/test_t1_fraction2_build.py tests/test_t1_fraction2_launcher.py \
  tests/test_t1_ms_build.py tests/test_t1_msb_build.py tests/test_t1_candidates_build.py tests/test_t1_bets.py \
  tests/test_hidden_counts.py tests/test_s2_rules.py \
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
build_base() {  # folder name args: the build on its own path with every step at T1_s2rcMS20's value; exit 12
  local d="$O/$1" name=$2
  BUILT+=("$1/$2.h5ad")
  shift 2
  mkdir -p "$d"
  $PY experiments/builds/t1_r66_check.py base t1 "$@" --out "$d/$name.h5ad" > /dev/null 2>&1 || return 12
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
identity MS20 T1_s2rcMS20 "$MS20_SHA" 39 $MS20_ARGS
for id in $IDS; do
  case "$id" in
    d10) identity d10 T1_s2rcMS20d10 "$D10_SHA" 71 $D10_ARGS ;;
    d15) identity d15 T1_s2rcMS20d15 "$D15_SHA" 72 $D15_ARGS ;;
    d20) identity d20 T1_s2rcMS20d20 "$D20_SHA" 73 $D20_ARGS ;;
    H05) identity H05 T1_s2rcMS20H05 "$H05_SHA" 74 $H05_ARGS ;;
    H10) identity H10 T1_s2rcMS20H10 "$H10_SHA" 75 $H10_ARGS ;;
    m10) identity m10 T1_s2rcMS20m10 "$M10_SHA" 76 $M10_ARGS ;;
    m20) identity m20 T1_s2rcMS20m20 "$M20_SHA" 77 $M20_ARGS ;;
  esac
done

STEP=${STEM}base
begin
build_base "$STEP" "${FILE}_base" $ARGS || finish $?
same "$STEP" "${FILE}_base" "$MS20_SHA" 51
step_ok

for STEP in "${STEM}a" "${STEM}b"; do begin; build "$STEP" "$FILE" $ARGS || finish $?; step_ok; done

STEP=done
begin
(cd "$O" && sha256sum "${BUILT[@]}") > "$S/SHA256SUMS" || finish 40
cp "$S/SHA256SUMS" "$O/SHA256SUMS" || finish 40
cmp -s "$O/${STEM}a/$FILE.h5ad" "$O/${STEM}b/$FILE.h5ad" || finish 41
step_ok
