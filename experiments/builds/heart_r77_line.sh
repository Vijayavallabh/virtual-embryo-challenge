#!/usr/bin/env bash
# Row R77's two files (experiments/harness/TRIES.md row R77, the captain's order of 2026-11-26, ruling, no run; decision
# log 2026-11-26 (h); the captain's named exception on the heart's scope and size for each file): the heart's refinement
# round 4 on OTK125's command (row R55's category file: build_test.py heart --entry organ_scope --trend-strength 1.25
# --trend-form category --z-mode raw --captain-ruling R55, then halve.py at seed 1000001) with one change, the category
# form's strength, each only with --captain-ruling R77, on row R74's route in build_test.py with the strength taken from
# the builder's HEART_R77_TABLE. One file per call, given by its name:
#   OTK05   --trend-strength 0.5 --trend-form category    T2heart_organ_ts05_tk.h5ad, half T2heart_organ_ts05_tk_half.h5ad
#   OTK20   --trend-strength 2.0 --trend-form category    T2heart_organ_ts20_tk.h5ad and its half
# each as build_test.py heart --entry organ_scope <OPTION> --z-mode raw --captain-ruling R77 (the builder's default seed
# and board), then halved by halve.py --seed 1000001; the half goes up. Test-board builds on a DGX campaign checkout or
# detached worktree that holds the released stages (data/raw/T2) and its own venv, made there by `uv sync --frozen`
# (never the reference checkout). Run from its root, detached, and wait on the exit files with a background poller;
# never a foreground sleep. The two calls may run at the same time, each with its own RUN:
#   bash experiments/builds/heart_r77_line.sh <COMMIT> <RUN> <FILE_NAME>
# FILE_NAME is one of the two in the literal case table below (tests/test_heart_r77_line.py holds it equal to the
# builder's HEART_R77_TABLE, the one table row R77's strings are checked against); any other name (OTK125, OTK5, OTK2,
# OTK050, OTK200, otk05 and row R74's four among them), or a fourth argument, is refused before any folder is made. RUN
# names this run's own scratch folder ~/v/scratch/heartr77/<RUN> and its output folder runs/heartr77_<commit7>_<RUN>.
# Neither may exist: a run folder is never reused and nothing is ever deleted here, so a second run takes a new name.
# In order, each step first writing the word running to ~/v/scratch/heartr77/<RUN>/<STEP>.exit, then its exit code:
#   gate    HEAD is COMMIT (the full hash), the tracked tree is clean, docs/S2_RULES.sha256 checks, TRIES.md holds
#           exactly one row R77 whose first cell starts "R77 (the captain's ", says "ruling, no run" and names
#           `--entry organ_scope`, `--z-mode raw`, FILE_NAME's option string and `--captain-ruling R77`, each as a whole
#           backticked string, the option string followed there by " (FILE_NAME)" as the row pairs each string with its
#           file; the five released stages the builder reads (data/raw/T2: E10.5, E9.5, E8.75, E8.5 and E8.25_late) are
#           in this checkout, and .venv's vec.pth names this checkout's src; then it writes the head file
#   tests   the builder's and this launcher's tests, row R74's and the earlier heart lines' targeted tests (synthetic
#           stages only), 8 BLAS threads
#   OTK125  the identity step: OTK125 by its committed command, then its half by halve.py --seed 1000001; both must give
#           their recorded sha256 (experiments/builds/results/hearth_8a8bc3a_h2.SHA256SUMS), or everything stops
#   <FILE_NAME>a <FILE_NAME>b  the file by its command, twice, into separate folders, each then halved at seed 1000001
#   done    one sha256 list over every built file (written first, in the scratch and output folders); the two halves and
#           the two full files compared byte for byte; the file and its half must have neither OTK125's bytes (a build
#           at R55's strength 1.25 would), nor the recorded bytes of a file of row R74 (R74_SHAS below, its four lists
#           of experiments/builds/results), nor those of row R77's other file as any list of this line records them
#           (every run's scratch list ~/v/scratch/heartr77/*/SHA256SUMS, a run at the same time included, and any
#           committed experiments/builds/results/heartr77_*.SHA256SUMS)
# A build step also fails (15) unless its stamp says this commit and clean True; `vec check` runs on each half, the file
# that goes up. Builds run at 16 BLAS threads, as the reference builds. Every builder, halve.py and `vec check` output
# goes to /dev/null: a test build is read by its exit code, its sha256 and its stamp's commit and clean fields only
# (AGENTS.md), never its log, size, cell count or sidecar record. A build that exits non-zero means no file, and the
# reason is not read here. Nothing here uploads anything.
# Exit codes: 89 FILE_NAME not one of row R77's two (or a fourth argument), 90 an odd RUN or COMMIT, 91 the run folder
# exists, 96 HEAD is not COMMIT, 98 tracked tree not clean, 97 the rules pin, 95 row R77 not as above, 92 a released
# stage not in this checkout, 88 .venv not this checkout's, 93 the output folder exists, 10 the tests, 11 a build failed,
# 12 halve.py, 13 vec check, 14 no stamp, 15 stamp not this commit and clean, 20 OTK125's file is not its recorded
# bytes, 21 its half is not its recorded half, 30 the hash list, 31 the two builds differ, 32 the file or its half has
# OTK125's bytes, 33 the file or its half has the recorded bytes of a file of row R74, 34 the file or its half has the
# recorded bytes of row R77's other file (or a list of this line could not be read).
set -u
set -o pipefail
COMMIT=${1:?usage: heart_r77_line.sh COMMIT RUN FILE_NAME}
RUN=${2:?usage: heart_r77_line.sh COMMIT RUN FILE_NAME}
FILE=${3:?usage: heart_r77_line.sh COMMIT RUN FILE_NAME}
[ "$#" = 3 ] || { echo "refused: give COMMIT RUN FILE_NAME and nothing else" >&2; exit 89; }
case "$RUN" in
  *[!A-Za-z0-9_-]*) echo "refused: RUN takes letters, digits, _ and - only" >&2; exit 90 ;;
esac
case "$COMMIT" in
  *[!0-9a-f]*) echo "refused: COMMIT must be a full hash" >&2; exit 90 ;;
esac
[ "${#COMMIT}" = 40 ] || { echo "refused: COMMIT must be a full hash" >&2; exit 90; }
# FILE_NAME -> its option string and its file's stem, as row R77 pairs them (the builder's HEART_R77_TABLE and
# HEART_R77_STEMS, in the row's order); a literal table: no name is read by a pattern
case "$FILE" in
  OTK05) OPT="--trend-strength 0.5 --trend-form category"; STEM=T2heart_organ_ts05_tk ;;
  OTK20) OPT="--trend-strength 2.0 --trend-form category"; STEM=T2heart_organ_ts20_tk ;;
  *) echo "refused: FILE_NAME is not one of row R77's two files" >&2; exit 89 ;;
esac
R77_STEMS="T2heart_organ_ts05_tk T2heart_organ_ts20_tk"  # both files' stems, the table's order
ARGS="--entry organ_scope $OPT --z-mode raw --captain-ruling R77"
HALF_SEED=1000001
export PATH=$HOME/.local/bin:$PATH
unset UV_PROJECT_ENVIRONMENT  # uv runs this checkout's own .venv, the one the gate checks
S=$HOME/v/scratch/heartr77/$RUN
mkdir -p "$HOME/v/scratch/heartr77"
mkdir "$S" 2> /dev/null || { echo "refused: the run folder $S exists (a run folder is never reused)" >&2; exit 91; }
mkdir "$S/tmp"
export NUMPY_MADVISE_HUGEPAGE=0 TMPDIR="$S/tmp"
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
PY="uv run --frozen --no-sync python -I -B"
BOARD=T2:heart:test_extrap
OTK125_SHA=d4d8d75f8bed48813bf143f5d7726828e5c30cdb27bd708b6fc2d9a455d21893       # hearth_8a8bc3a_h2: OTK125's file
OTK125_HALF_SHA=9835fa5d47d7021760beb2056ce44c2c846ded6fc31dacf64bbb2d9f0f7d1b54  # hearth_8a8bc3a_h2: OTK125's half
OTK125_ARGS="--entry organ_scope --trend-strength 1.25 --trend-form category --z-mode raw --captain-ruling R55"
# row R74's four files and their halves as recorded (experiments/builds/results/heartr74_4b3e7a9_r74k10, _r74k15,
# _r74k075 and _r74k175.SHA256SUMS), file then half, in row R74's order
R74_SHAS="2692b14dc52098468621945a3f62d6ad5cde6bd71fed9217963938f134fdaae2
904d58ab182a69eaa58c2a89270f30753566cd1eb407314c30b17509ebf9d414
8971fec19eb4cb9e4b6466673b1429dd37c67dc935d9db89041a45882070fedb
27cca412a87f0fd27f006851f91d8cf59bd01aebce5d9de6d33e07da0a9f012e
b01cecc7918dd4116eec9bc9e33ed3beaff4a5b4f89e74839fabe56856d0d643
9dcfb21df2eb942a9f2876033e4e13a7b62b00315e6cbf49f7f82e06cceb3168
14591c2347b0ebefbd08ed9963da70f232649684fb1662139d950ef4fed71452
efaa72be2bbf6cb20d6d6c87b0b9e40f2ff25b5f1daba0831de2071064c49e0b"
STAGES="E10.5 E9.5 E8.75 E8.5 E8.25_late"
STEP=gate
finish() { echo "$1" > "$S/$STEP.exit"; exit "$1"; }
begin() { echo running > "$S/$STEP.exit"; }
step_ok() { echo 0 > "$S/$STEP.exit"; }
begin

head=$(git rev-parse HEAD)
[ "$head" = "$COMMIT" ] || { echo "refused: HEAD is not $COMMIT" >&2; finish 96; }
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
r77=$(first_cell R77 95) || exit $?
for opt in "--entry organ_scope" "--z-mode raw" "$OPT" "--captain-ruling R77"; do
  printf '%s\n' "$r77" | grep -o '`[^`]*`' | grep -qxF -- "\`$opt\`" \
    || { echo "refused: row R77's first cell does not name \`$opt\` as a whole backticked string" >&2; finish 95; }
done
printf '%s\n' "$r77" | grep -qF -- "\`$OPT\` ($FILE)" \
  || { echo "refused: row R77's first cell does not pair \`$OPT\` with ($FILE)" >&2; finish 95; }
for stage in $STAGES; do
  [ -f "data/raw/T2/$stage.h5ad" ] \
    || { echo "refused: the released stage data/raw/T2/$stage.h5ad is not in this checkout" >&2; finish 92; }
done
venv=0
for pth in .venv/lib/python3*/site-packages/vec.pth; do
  if [ -f "$pth" ] && [ "$(readlink -f "$(< "$pth")")" = "$(readlink -f src)" ]; then venv=1; fi
done
[ "$venv" = 1 ] || {
  echo "refused: this checkout has no .venv of its own (uv sync --frozen here) whose vec.pth names its src" >&2
  finish 88
}
O=runs/heartr77_${head:0:7}_$RUN
mkdir -p runs
mkdir "$O" 2> /dev/null || { echo "refused: the output folder $O exists" >&2; finish 93; }
echo "$head" > "$S/head"
step_ok

STEP=tests
begin
$PY -m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp="$S/tmp/pytest" \
  tests/test_heart_r77_line.py tests/test_heart_r77_build.py tests/test_heart_r74_line.py tests/test_heart_r74_build.py \
  tests/test_heart_trend_line_h.py tests/test_heart_wbh_build.py tests/test_heart_xpc_build.py \
  tests/test_heart_oic_build.py tests/test_heart_trend_line_g.py tests/test_heart_trend_line_f.py \
  tests/test_heart_trend_line_d.py tests/test_heart_trend_line.py tests/test_build_heart_organ_scope.py \
  tests/test_heart_trend_build.py tests/test_b0_trend.py tests/test_halve.py tests/test_hidden_counts.py \
  tests/test_s2_rules.py \
  > "$S/tests.out" 2>&1 || finish 10
step_ok

export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16
stamp() {  # the stamp's commit and clean fields only
  $PY -c 'import json,sys; s=json.load(open(sys.argv[1])); print(s["commit"], s["clean"])' "$1.stamp.json"
}
sha() { sha256sum "$1" | cut -d' ' -f1; }
BUILT=()  # every built file, folder/name.h5ad under $O, in build order
build() {  # folder stem args: the build, its half at HALF_SEED, `vec check` on the half, the stamp; exit 11 build,
           # 12 halve.py, 13 vec check, 14 no stamp, 15 stamp not this commit and clean
  local step=$1 name=$2 st
  local d="$O/$step"
  shift 2
  mkdir -p "$d"
  BUILT+=("$step/$name.h5ad" "$step/${name}_half.h5ad")
  $PY experiments/builds/build_test.py heart "$@" --out "$d/$name.h5ad" > /dev/null 2>&1 || return 11
  $PY experiments/builds/halve.py "$d/$name.h5ad" "$d/${name}_half.h5ad" --seed "$HALF_SEED" > /dev/null 2>&1 \
    || return 12
  uv run --frozen --no-sync vec check "$d/${name}_half.h5ad" --board "$BOARD" > /dev/null 2>&1 || return 13
  st=$(stamp "$d/$name.h5ad" 2> /dev/null) || return 14
  echo "$step/$name $st" >> "$S/stamps"
  [ "$st" = "$head True" ] || return 15
}
recorded_other() {  # the sha256 that the lists of this line record for row R77's other file and its half: every run's
                    # scratch list (a run at the same time included) and any committed heartr77 list
  local f s
  for f in "$HOME"/v/scratch/heartr77/*/SHA256SUMS experiments/builds/results/heartr77_*.SHA256SUMS; do
    [ -f "$f" ] || continue
    for s in $R77_STEMS; do
      [ "$s" = "$STEM" ] && continue
      awk -v a="$s.h5ad" -v b="${s}_half.h5ad" '{ n = $2; sub(/.*\//, "", n); if (n == a || n == b) print $1 }' "$f" \
        || return 34
    done
  done
}
: > "$S/stamps"

STEP=OTK125
begin
build OTK125 T2heart_organ_ts125_tk $OTK125_ARGS || finish $?
[ "$(sha "$O/OTK125/T2heart_organ_ts125_tk.h5ad")" = "$OTK125_SHA" ] \
  || { echo "stopped: OTK125's command does not give its recorded file" >&2; finish 20; }
[ "$(sha "$O/OTK125/T2heart_organ_ts125_tk_half.h5ad")" = "$OTK125_HALF_SHA" ] \
  || { echo "stopped: OTK125's half is not its recorded half" >&2; finish 21; }
step_ok

for STEP in "${FILE}a" "${FILE}b"; do begin; build "$STEP" "$STEM" $ARGS || finish $?; step_ok; done

STEP=done
begin
(cd "$O" && sha256sum "${BUILT[@]}") > "$S/SHA256SUMS" || finish 30
cp "$S/SHA256SUMS" "$O/SHA256SUMS" || finish 30
for f in "${STEM}_half" "$STEM"; do
  cmp -s "$O/${FILE}a/$f.h5ad" "$O/${FILE}b/$f.h5ad" || { echo "stopped: the two builds of $f differ" >&2; finish 31; }
done
file_sha=$(sha "$O/${FILE}a/$STEM.h5ad") || finish 30
half_sha=$(sha "$O/${FILE}a/${STEM}_half.h5ad") || finish 30
[ "$file_sha" != "$OTK125_SHA" ] && [ "$half_sha" != "$OTK125_HALF_SHA" ] \
  || { echo "stopped: $FILE or its half has OTK125's bytes (its strength was not taken)" >&2; finish 32; }
for x in $R74_SHAS; do
  [ "$file_sha" != "$x" ] && [ "$half_sha" != "$x" ] \
    || { echo "stopped: $FILE or its half has the recorded bytes of a file of row R74" >&2; finish 33; }
done
other=$(recorded_other) \
  || { echo "stopped: a sha256 list of this line could not be read" >&2; finish 34; }
for x in $other; do
  [ "$file_sha" != "$x" ] && [ "$half_sha" != "$x" ] \
    || { echo "stopped: $FILE or its half has the recorded bytes of row R77's other file" >&2; finish 34; }
done
step_ok
