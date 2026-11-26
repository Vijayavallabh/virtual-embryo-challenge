#!/usr/bin/env bash
# Row R75's five files (experiments/harness/TRIES.md, the captain's order of 2026-11-26, "ruling, no run": the embryo's
# refinement round on O, sent at once with no gate), each O's command (embryo --method fgw_shape0 --positions
# floor_exact --shape-interp overlay --shape-interp-test-ok) with one option string, each only with --captain-ruling
# R75: row R59's route plus the builder's one table EMBRYO_R75_TABLE. One file per call, given by its name:
#   T2embryo_s2x_ovt712.h5ad  --mix-lean later_half                 (t_mix 7/12)
#   T2embryo_s2x_ovt56.h5ad   --mix-lean later2                     (t_mix 5/6)
#   T2embryo_s2x_ovt512.h5ad  --mix-lean earlier_half               (t_mix 5/12)
#   T2embryo_s2xs_ovt23.h5ad  --expression series --mix-lean later  (SO's expression at row R59's t_mix 2/3)
#   T2embryo_s2x_ovd.h5ad     --expression depth                    (O's rows and positions, thinned to DEPTH_RULE)
# each after O's own options, then --captain-ruling R75, nothing else changed (default seed, test board). Test-board
# builds on a DGX campaign worktree that holds the released embryo stages (data/raw/T2/E7.5.h5ad and E8.0.h5ad; a link
# into the reference checkout's data is fine) and its own venv, made there by `uv sync --frozen` (never the reference
# checkout). Run from the root of that worktree, detached, and wait on the exit files with a background poller; never
# a foreground sleep:
#   bash experiments/builds/embryo_r75_line.sh <COMMIT> <RUN> <FILE_NAME>
# FILE_NAME is one of the five in the literal case table below (tests/test_embryo_r75_launcher.py holds it equal to the
# builder's EMBRYO_R75_TABLE, the one table row R75's strings are checked against); any other name, or a fourth
# argument, is refused before any folder is made. RUN names this run's own scratch folder ~/v/scratch/embryor75/<RUN>
# and its output folder runs/embryor75_<commit7>_<RUN>. Neither may exist: a run folder is never reused and nothing is
# ever deleted here, so a second run, or a second file, takes a new name. STEM is the file's short name in its steps
# (ovt712, ovt56, ovt512, sovt23 or ovd); BASES its base checks, each MODE:KNOWN, KNOWN the identity file whose bytes
# the file's own path must give with its one change put back by MODE (half: its t_mix at O's 1/2; mixture: its
# expression at O's mixture rows).
# In order, each step first writing the word running to ~/v/scratch/embryor75/<RUN>/<STEP>.exit, then its exit code:
#   gate   HEAD is COMMIT (the full hash), the tracked tree is clean, docs/S2_RULES.sha256 checks, TRIES.md holds exactly
#          one row R75 whose first cell starts "R75 (the captain's order", says "ruling, no run" and names FILE_NAME's
#          option string and `--captain-ruling R75`, each as a whole backticked string, the option string followed there
#          by " (FILE_NAME)" as the row pairs each string with its file, the released stages data/raw/T2/E7.5.h5ad and
#          data/raw/T2/E8.0.h5ad are in this checkout, .venv's vec.pth names this checkout's src, and each identity pin
#          below is the one its committed SHA256SUMS list holds (never a .txt record); then it writes the head file
#   tests  the builder's, checker's and launchers' targeted tests (synthetic stages only), 8 BLAS threads
#   O OT SO W  four identity rebuilds, each once, each of which must give its recorded bytes or everything stops:
#          O   T2embryo_s2x_ov     O_SHA   embov_3b813b5.SHA256SUMS OVa (experiments/builds/embryo_overlay.sh's OV_ARGS)
#          OT  T2embryo_s2x_ovt23  OT_SHA  embot23_045a92e_o2.SHA256SUMS OTa (row R59's file; embryo_ot23.sh's OT_ARGS)
#          SO  T2embryo_s2xs_ov    SO_SHA  embsohi_93f832e_r1.SHA256SUMS SOa (embryo_so_ohi.sh's SO_ARGS)
#          W   T2embryo_s2x_ovsb   W_SHA   embovsb_20b2381_w1.SHA256SUMS Wa (row R48's file; embryo_ovsb.sh's W_ARGS)
#   STEMhalf STEMmixture  FILE_NAME's command once per base check through experiments/builds/embryo_r75_check.py, on
#          the file's own path with its one change put back: half (a lean's t_mix at 1/2) must give O's bytes for the
#          three leans and SO's for T2embryo_s2xs_ovt23.h5ad; mixture (the expression at O's mixture rows) must give
#          row R59's Ot23's bytes for T2embryo_s2xs_ovt23.h5ad and O's for T2embryo_s2x_ovd.h5ad (O's coordinates bit
#          for bit), or everything stops
#   STEMa STEMb  FILE_NAME by its command, twice, into separate folders
#   done   one sha256 list over every built file, the two builds of FILE_NAME compared byte for byte
# A build step also fails (15) unless its stamp says this commit and clean True. Builds run at 16 BLAS threads, as the
# reference builds and the earlier embryo launchers' (the recorded files of O, Ot23, SO and W were built so), at the
# builder's default seed on the test board (no dry run: row R75's strings refuse one). Every builder, checker build and
# `vec check` output goes to /dev/null: a test build is read by its exit code, its sha256 and its stamp's commit and
# clean fields (AGENTS.md), never its log, size, record or cell count. Nothing here uploads anything.
# Exit codes: 90 odd RUN or COMMIT, 89 FILE_NAME not one of row R75's five (or a fourth argument), 91 run folder exists,
# 96 HEAD is not COMMIT, 98 tracked tree not clean, 97 S2_RULES.sha256, 95 row R75 not as above, 92 a released stage
# missing, 88 .venv not this checkout's, 85 an identity pin not its committed list's, 93 output folder exists,
# 10 tests, 11 build, 12 checker build, 13 vec check, 14 no stamp, 15 stamp not this commit and clean, 21 O is not O's
# bytes, 22 OT is not Ot23's, 23 SO is not SO's, 24 W is not W's, 51 a half base check is not its known bytes,
# 52 a mixture base check is not its known bytes, 40 hash list, 41 STEMa and STEMb differ.
set -u
COMMIT=${1:?usage: embryo_r75_line.sh COMMIT RUN FILE_NAME}
RUN=${2:?usage: embryo_r75_line.sh COMMIT RUN FILE_NAME}
FILE=${3:?usage: embryo_r75_line.sh COMMIT RUN FILE_NAME}
[ "$#" = 3 ] || { echo "refused: give COMMIT RUN FILE_NAME and nothing else" >&2; exit 89; }
case "$RUN" in
  *[!A-Za-z0-9_-]*) echo "refused: RUN takes letters, digits, _ and - only" >&2; exit 90 ;;
esac
case "$COMMIT" in
  *[!0-9a-f]*) echo "refused: COMMIT must be a full hash" >&2; exit 90 ;;
esac
[ "${#COMMIT}" = 40 ] || { echo "refused: COMMIT must be a full hash" >&2; exit 90; }
# FILE_NAME -> its option string, its short name and its base checks, as row R75 pairs them (the builder's
# EMBRYO_R75_TABLE, in the row's order); a literal table: no name is read by a pattern
case "$FILE" in
  T2embryo_s2x_ovt712.h5ad) OPT="--mix-lean later_half"; STEM=ovt712; BASES="half:O" ;;
  T2embryo_s2x_ovt56.h5ad) OPT="--mix-lean later2"; STEM=ovt56; BASES="half:O" ;;
  T2embryo_s2x_ovt512.h5ad) OPT="--mix-lean earlier_half"; STEM=ovt512; BASES="half:O" ;;
  T2embryo_s2xs_ovt23.h5ad) OPT="--expression series --mix-lean later"; STEM=sovt23; BASES="half:SO mixture:OT" ;;
  T2embryo_s2x_ovd.h5ad) OPT="--expression depth"; STEM=ovd; BASES="mixture:O" ;;
  *) echo "refused: FILE_NAME is not one of row R75's five files" >&2; exit 89 ;;
esac
NAME=${FILE%.h5ad}
ARGS="--method fgw_shape0 --positions floor_exact --shape-interp overlay --shape-interp-test-ok $OPT --captain-ruling R75"
export PATH=$HOME/.local/bin:$PATH
unset UV_PROJECT_ENVIRONMENT  # uv runs this checkout's own .venv, the one the gate checks
S=$HOME/v/scratch/embryor75/$RUN
mkdir -p "$HOME/v/scratch/embryor75"
mkdir "$S" 2> /dev/null || { echo "refused: the run folder $S exists (a run folder is never reused)" >&2; exit 91; }
mkdir "$S/tmp"
export NUMPY_MADVISE_HUGEPAGE=0 TMPDIR="$S/tmp"
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
PY="uv run --frozen --no-sync python -I -B"
BOARD=T2:embryo:test_interp
O_SHA=1cfc91be3a20bfdd1a4f7bcfebd7316e1143c6d4048e578c1f30bcfbe1e0e1ea  # embov_3b813b5, the overlay file O
OT_SHA=753900fab42738a8695d8048645cccd953ded07cd15539101964c26259aad22f  # embot23_045a92e_o2, row R59's Ot23
SO_SHA=77062082d97555961901c811d8a1aaf9b5b2403417c19a6988c816ff553f1165  # embsohi_93f832e_r1, SO
W_SHA=f576b64b8b93d79b1269b580df31bd5c6e58246cd4f3eda12cba7d77387087ab  # embovsb_20b2381_w1, row R48's W
O_SUMS=experiments/builds/results/embov_3b813b5.SHA256SUMS
O_ENTRY=OVa/T2embryo_s2x_ov.h5ad
OT_SUMS=experiments/builds/results/embot23_045a92e_o2.SHA256SUMS
OT_ENTRY=OTa/T2embryo_s2x_ovt23.h5ad
SO_SUMS=experiments/builds/results/embsohi_93f832e_r1.SHA256SUMS
SO_ENTRY=SOa/T2embryo_s2xs_ov.h5ad
W_SUMS=experiments/builds/results/embovsb_20b2381_w1.SHA256SUMS
W_ENTRY=Wa/T2embryo_s2x_ovsb.h5ad
O_ARGS="--method fgw_shape0 --positions floor_exact --expression mixture --shape-interp overlay --shape-interp-test-ok"
OT_ARGS="--method fgw_shape0 --positions floor_exact --expression mixture --shape-interp overlay --shape-interp-test-ok --mix-lean later --captain-ruling R59"
SO_ARGS="--method fgw_shape0 --positions floor_exact --expression series --shape-interp overlay --shape-interp-test-ok"
W_ARGS="--method fgw_shape0 --positions floor_exact --expression series_back --shape-interp overlay --shape-interp-test-ok --captain-ruling R48"
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
[ "$(grep -c "^| R75 " experiments/harness/TRIES.md)" = 1 ] \
  || { echo "refused: TRIES.md does not hold exactly one row R75" >&2; finish 95; }
r75=$(grep "^| R75 " experiments/harness/TRIES.md | cut -d'|' -f2)
printf '%s\n' "$r75" | grep -q "^ R75 (the captain's order" \
  || { echo "refused: row R75's first cell is not the captain's order" >&2; finish 95; }
printf '%s\n' "$r75" | grep -qF -- "ruling, no run" \
  || { echo "refused: row R75's first cell does not say ruling, no run" >&2; finish 95; }
for opt in "$OPT" "--captain-ruling R75"; do
  printf '%s\n' "$r75" | grep -o '`[^`]*`' | grep -qxF -- "\`$opt\`" \
    || { echo "refused: row R75's first cell does not name \`$opt\` as a whole backticked string" >&2; finish 95; }
done
printf '%s\n' "$r75" | grep -qF -- "\`$OPT\` ($FILE)" \
  || { echo "refused: row R75's first cell does not pair \`$OPT\` with ($FILE)" >&2; finish 95; }
for stage in E7.5 E8.0; do
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
for rec in "O_SHA $O_SUMS $O_ENTRY" "OT_SHA $OT_SUMS $OT_ENTRY" "SO_SHA $SO_SUMS $SO_ENTRY" "W_SHA $W_SUMS $W_ENTRY"; do
  read -r var sums entry <<< "$rec"
  got=$(awk -v e="$entry" '$2 == e {print $1}' "$sums" 2> /dev/null)
  [ "$got" = "${!var}" ] || { echo "refused: $sums does not give $var for $entry" >&2; finish 85; }
done
O=runs/embryor75_${head:0:7}_$RUN
mkdir -p runs
mkdir "$O" 2> /dev/null || { echo "refused: the output folder $O exists" >&2; finish 93; }
echo "$head" > "$S/head"
step_ok

STEP=tests
begin
$PY -m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp="$S/tmp/pytest" \
  tests/test_embryo_r75_launcher.py tests/test_embryo_r75_build.py tests/test_mix_lean_build.py \
  tests/test_embryo_ot23_launcher.py tests/test_series_back_build.py tests/test_embryo_ovsb_launcher.py \
  tests/test_embryo_so_ohi_launcher.py tests/test_bf_build.py tests/test_embryo_overlay.py tests/test_shape_interp.py \
  tests/test_shape_interp_build.py tests/test_e47_floor_exact.py tests/test_t2_build.py tests/test_hidden_counts.py \
  tests/test_s2_rules.py \
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
  $PY experiments/builds/build_test.py embryo "$@" --out "$d/$name.h5ad" > /dev/null 2>&1 || return 11
  checked "$d" "$name"
}
build_base() {  # folder name mode args: the build on the file's own path with its one change put back; exit 12
  local d="$O/$1" name=$2 mode=$3
  BUILT+=("$1/$2.h5ad")
  shift 3
  mkdir -p "$d"
  $PY experiments/builds/embryo_r75_check.py "$mode" embryo "$@" --out "$d/$name.h5ad" > /dev/null 2>&1 || return 12
  checked "$d" "$name"
}
same() {  # folder name sha code known: stop with the code unless the built file has the known file's bytes
  local got
  got=$(sha256sum "$O/$1/$2.h5ad" | cut -d' ' -f1)
  [ "$got" = "$3" ] || { echo "stopped: step $1 does not give the known bytes of $5" >&2; finish "$4"; }
}
: > "$S/stamps"
identity() {  # STEP name sha code args: one identity rebuild
  STEP=$1
  local name=$2 sha=$3 code=$4
  shift 4
  begin
  build "$STEP" "$name" "$@" || finish $?
  same "$STEP" "$name" "$sha" "$code" "$name"
  step_ok
}
identity O T2embryo_s2x_ov "$O_SHA" 21 $O_ARGS
identity OT T2embryo_s2x_ovt23 "$OT_SHA" 22 $OT_ARGS
identity SO T2embryo_s2xs_ov "$SO_SHA" 23 $SO_ARGS
identity W T2embryo_s2x_ovsb "$W_SHA" 24 $W_ARGS

for base in $BASES; do
  mode=${base%%:*}
  known=${base#*:}
  sha=${known}_SHA
  case "$mode" in
    half) code=51 ;;
    mixture) code=52 ;;
  esac
  STEP=$STEM$mode
  begin
  build_base "$STEP" "${NAME}_$mode" "$mode" $ARGS || finish $?
  same "$STEP" "${NAME}_$mode" "${!sha}" "$code" "$known"
  step_ok
done

for STEP in "${STEM}a" "${STEM}b"; do begin; build "$STEP" "$NAME" $ARGS || finish $?; step_ok; done

STEP=done
begin
(cd "$O" && sha256sum "${BUILT[@]}") > "$S/SHA256SUMS" || finish 40
cp "$S/SHA256SUMS" "$O/SHA256SUMS" || finish 40
cmp -s "$O/${STEM}a/$FILE" "$O/${STEM}b/$FILE" || finish 41
step_ok
