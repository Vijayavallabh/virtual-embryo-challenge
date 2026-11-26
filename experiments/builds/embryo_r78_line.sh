#!/usr/bin/env bash
# Row R78's two files (experiments/harness/TRIES.md, the captain's order of 2026-11-26, "ruling, no run": the embryo's
# fourth refinement round on O, both sent at once with no gate), each O's command (embryo --method fgw_shape0
# --positions floor_exact --shape-interp overlay --shape-interp-test-ok) with one option string, each only with
# --captain-ruling R78: row R75's route plus the builder's one table EMBRYO_R78_TABLE. One file per call, given by its
# name:
#   T2embryo_s2x_ovt1112.h5ad  --mix-lean later_11_12  (t_mix 11/12)
#   T2embryo_s2x_ovt13.h5ad    --mix-lean earlier      (t_mix 1/3)
# each after O's own options, then --captain-ruling R78, nothing else changed (default seed, test board). Test-board
# builds on a DGX campaign worktree that holds the released embryo stages (data/raw/T2/E7.5.h5ad and E8.0.h5ad; a link
# into the reference checkout's data is fine) and its own venv, made there by `uv sync --frozen` (never the reference
# checkout). Run from the root of that worktree, detached, and wait on the exit files with a background poller; never
# a foreground sleep:
#   bash experiments/builds/embryo_r78_line.sh <COMMIT> <RUN> <FILE_NAME>
# FILE_NAME is one of the two in the literal case table below (tests/test_embryo_r78_launcher.py holds it equal to the
# builder's EMBRYO_R78_TABLE, the one table row R78's strings are checked against); any other name, or a fourth
# argument, is refused before any folder is made, and so is a run started after 10:00 UTC on 2026-11-26 (NOT_AFTER,
# read from bash's own clock: every R78 file is sent at once, the last upload is at 11:59 UTC and a run takes about an
# hour). RUN names this run's own scratch folder ~/v/scratch/embryor78/<RUN> and its output folder
# runs/embryor78_<commit7>_<RUN>. Neither may exist: a run folder is never reused and nothing is ever deleted here, so
# a second run, or a second file, takes a new name. STEM is the file's short name in its steps (ovt1112 or ovt13);
# R75=yes marks the one run that also rebuilds row R75's Ot56 and Ot512 (T2embryo_s2x_ovt1112.h5ad's, the row's first
# file; the other run does not).
# In order, each step first writing the word running to ~/v/scratch/embryor78/<RUN>/<STEP>.exit, then its exit code:
#   gate   HEAD is COMMIT (the full hash), the tracked tree is clean, docs/S2_RULES.sha256 checks, TRIES.md holds exactly
#          one row R78 whose first cell starts "R78 (the captain's order", says "ruling, no run" and names FILE_NAME's
#          option string and `--captain-ruling R78`, each as a whole backticked string, the option string followed there
#          by " (FILE_NAME)" as the row pairs each string with its file, the released stages data/raw/T2/E7.5.h5ad and
#          data/raw/T2/E8.0.h5ad are in this checkout, .venv's vec.pth names this checkout's src, and each pin below is
#          the one its committed SHA256SUMS list holds (never a .txt record); then it writes the head file
#   tests  the builder's, checker's and launchers' targeted tests (synthetic stages only): the test lists of
#          embryo_r75_line.sh, embryo_ot23.sh, embryo_bf.sh and coords/run_e58.sh and this line's own two,
#          8 BLAS threads
#   O      O's identity rebuild (T2embryo_s2x_ov, O_SHA, embov_3b813b5.SHA256SUMS OVa; embryo_overlay.sh's OV_ARGS),
#          which must give its recorded bytes or everything stops
#   ovt56 ovt512  in the R75=yes run only: row R75's T2embryo_s2x_ovt56.h5ad and T2embryo_s2x_ovt512.h5ad by row R75's
#          own commands (embryo_r75_line.sh's ARGS), each of which must give its recorded bytes (OT56_SHA,
#          embryor75_ee22b75_r56.SHA256SUMS ovt56a; OT512_SHA, embryor75_ee22b75_r512.SHA256SUMS ovt512a) or everything
#          stops
#   STEMhalf  FILE_NAME's command once through experiments/builds/embryo_r78_check.py half, on the file's own path with
#          its t_mix put back at O's 1/2: it must give O's bytes, or everything stops
#   STEMa STEMb  FILE_NAME by its command, twice, into separate folders
#   done   one sha256 list over every built file, the two builds of FILE_NAME compared byte for byte, and FILE_NAME's
#          sha256 must differ from O's, row R59's Ot23's and each of row R75's five files' (their committed lists)
# A build step also fails (15) unless its stamp says this commit and clean True. Builds run at 16 BLAS threads, as the
# reference builds and the earlier embryo launchers' (the recorded files were built so), at the builder's default seed
# on the test board (no dry run: row R78's strings refuse one). Every builder, checker build and `vec check` output
# goes to /dev/null: a test build is read by its exit code, its sha256 and its stamp's commit and clean fields
# (AGENTS.md), never its log, size, record or cell count. Nothing here uploads anything.
# Exit codes: 90 odd RUN or COMMIT, 89 FILE_NAME not one of row R78's two (or a fourth argument), 87 a run started
# after 10:00 UTC on 2026-11-26, 91 run folder exists, 96 HEAD is not COMMIT, 98 tracked tree not clean,
# 97 S2_RULES.sha256, 95 row R78 not as above, 92 a released stage missing, 88 .venv not this checkout's, 85 a pin not
# its committed list's, 93 output folder exists, 10 tests, 11 build, 12 checker build, 13 vec check, 14 no stamp,
# 15 stamp not this commit and clean, 21 O is not O's bytes, 25 ovt56 is not Ot56's, 26 ovt512 is not Ot512's, 51 the
# half check is not O's bytes, 40 hash list, 41 STEMa and STEMb differ, 42 FILE_NAME has a known file's bytes.
set -u
COMMIT=${1:?usage: embryo_r78_line.sh COMMIT RUN FILE_NAME}
RUN=${2:?usage: embryo_r78_line.sh COMMIT RUN FILE_NAME}
FILE=${3:?usage: embryo_r78_line.sh COMMIT RUN FILE_NAME}
[ "$#" = 3 ] || { echo "refused: give COMMIT RUN FILE_NAME and nothing else" >&2; exit 89; }
case "$RUN" in
  *[!A-Za-z0-9_-]*) echo "refused: RUN takes letters, digits, _ and - only" >&2; exit 90 ;;
esac
case "$COMMIT" in
  *[!0-9a-f]*) echo "refused: COMMIT must be a full hash" >&2; exit 90 ;;
esac
[ "${#COMMIT}" = 40 ] || { echo "refused: COMMIT must be a full hash" >&2; exit 90; }
# FILE_NAME -> its option string, its short name and whether its run rebuilds row R75's Ot56 and Ot512, as row R78
# pairs them (the builder's EMBRYO_R78_TABLE, in the row's order); a literal table: no name is read by a pattern
case "$FILE" in
  T2embryo_s2x_ovt1112.h5ad) OPT="--mix-lean later_11_12"; STEM=ovt1112; R75=yes ;;
  T2embryo_s2x_ovt13.h5ad) OPT="--mix-lean earlier"; STEM=ovt13; R75=no ;;
  *) echo "refused: FILE_NAME is not one of row R78's two files" >&2; exit 89 ;;
esac
NOT_AFTER=1795687200  # 2026-11-26 10:00:00 UTC: no run of this launcher starts after it
printf -v now '%(%s)T' -1  # bash's own clock, in seconds since the epoch
[ "$now" -le "$NOT_AFTER" ] || { echo "refused: no run starts after 10:00 UTC on 2026-11-26" >&2; exit 87; }
NAME=${FILE%.h5ad}
ARGS="--method fgw_shape0 --positions floor_exact --shape-interp overlay --shape-interp-test-ok $OPT --captain-ruling R78"
export PATH=$HOME/.local/bin:$PATH
unset UV_PROJECT_ENVIRONMENT  # uv runs this checkout's own .venv, the one the gate checks
S=$HOME/v/scratch/embryor78/$RUN
mkdir -p "$HOME/v/scratch/embryor78"
mkdir "$S" 2> /dev/null || { echo "refused: the run folder $S exists (a run folder is never reused)" >&2; exit 91; }
mkdir "$S/tmp"
export NUMPY_MADVISE_HUGEPAGE=0 TMPDIR="$S/tmp"
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
PY="uv run --frozen --no-sync python -I -B"
BOARD=T2:embryo:test_interp
R=experiments/builds/results
O_SHA=1cfc91be3a20bfdd1a4f7bcfebd7316e1143c6d4048e578c1f30bcfbe1e0e1ea  # embov_3b813b5, the overlay file O
OT56_SHA=c58fe13bbca497fd8b7f11abab89094b6583b4fc6beebe0fa1a2c1db386d213c  # embryor75_ee22b75_r56, row R75's Ot56
OT512_SHA=598db8bb7cd125dabe45af22da7496a0d20b072f967a4ad4c87dc26bbeced7b5  # embryor75_ee22b75_r512, row R75's Ot512
OT712_SHA=ab14288faad3ad10c990073bb370c25dbb2cb3258d605cea28d5774d2dd35a42  # embryor75_ee22b75_r712, row R75's Ot712
SOT23_SHA=baaeff58340a0ab60fa921f11fdead7effd1e9df85bf61330caec673cb257d06  # embryor75_ee22b75_rs23, row R75's SOt23
OVD_SHA=96a13734d618eea3329c0a843c6e0e5b9ec2369d17a6c37aff2f7994733fbda8  # embryor75_ee22b75_rd, row R75's Od
OT23_SHA=753900fab42738a8695d8048645cccd953ded07cd15539101964c26259aad22f  # embot23_045a92e_o2, row R59's Ot23
O_SUMS=$R/embov_3b813b5.SHA256SUMS
O_ENTRY=OVa/T2embryo_s2x_ov.h5ad
OT56_SUMS=$R/embryor75_ee22b75_r56.SHA256SUMS
OT56_ENTRY=ovt56a/T2embryo_s2x_ovt56.h5ad
OT512_SUMS=$R/embryor75_ee22b75_r512.SHA256SUMS
OT512_ENTRY=ovt512a/T2embryo_s2x_ovt512.h5ad
OT712_SUMS=$R/embryor75_ee22b75_r712.SHA256SUMS
OT712_ENTRY=ovt712a/T2embryo_s2x_ovt712.h5ad
SOT23_SUMS=$R/embryor75_ee22b75_rs23.SHA256SUMS
SOT23_ENTRY=sovt23a/T2embryo_s2xs_ovt23.h5ad
OVD_SUMS=$R/embryor75_ee22b75_rd.SHA256SUMS
OVD_ENTRY=ovda/T2embryo_s2x_ovd.h5ad
OT23_SUMS=$R/embot23_045a92e_o2.SHA256SUMS
OT23_ENTRY=OTa/T2embryo_s2x_ovt23.h5ad
KNOWN="O_SHA OT23_SHA OT712_SHA OT56_SHA OT512_SHA SOT23_SHA OVD_SHA"  # FILE_NAME's sha256 must be none of these
O_ARGS="--method fgw_shape0 --positions floor_exact --expression mixture --shape-interp overlay --shape-interp-test-ok"
OT56_ARGS="--method fgw_shape0 --positions floor_exact --shape-interp overlay --shape-interp-test-ok --mix-lean later2 --captain-ruling R75"
OT512_ARGS="--method fgw_shape0 --positions floor_exact --shape-interp overlay --shape-interp-test-ok --mix-lean earlier_half --captain-ruling R75"
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
[ "$(grep -c "^| R78 " experiments/harness/TRIES.md)" = 1 ] \
  || { echo "refused: TRIES.md does not hold exactly one row R78" >&2; finish 95; }
r78=$(grep "^| R78 " experiments/harness/TRIES.md | cut -d'|' -f2)
printf '%s\n' "$r78" | grep -q "^ R78 (the captain's order" \
  || { echo "refused: row R78's first cell is not the captain's order" >&2; finish 95; }
printf '%s\n' "$r78" | grep -qF -- "ruling, no run" \
  || { echo "refused: row R78's first cell does not say ruling, no run" >&2; finish 95; }
for opt in "$OPT" "--captain-ruling R78"; do
  printf '%s\n' "$r78" | grep -o '`[^`]*`' | grep -qxF -- "\`$opt\`" \
    || { echo "refused: row R78's first cell does not name \`$opt\` as a whole backticked string" >&2; finish 95; }
done
printf '%s\n' "$r78" | grep -qF -- "\`$OPT\` ($FILE)" \
  || { echo "refused: row R78's first cell does not pair \`$OPT\` with ($FILE)" >&2; finish 95; }
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
for var in $KNOWN; do
  sums=${var%_SHA}_SUMS
  entry=${var%_SHA}_ENTRY
  got=$(awk -v e="${!entry}" '$2 == e {print $1}' "${!sums}" 2> /dev/null)
  [ "$got" = "${!var}" ] || { echo "refused: ${!sums} does not give $var for ${!entry}" >&2; finish 85; }
done
O=runs/embryor78_${head:0:7}_$RUN
mkdir -p runs
mkdir "$O" 2> /dev/null || { echo "refused: the output folder $O exists" >&2; finish 93; }
echo "$head" > "$S/head"
step_ok

STEP=tests
begin
$PY -m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp="$S/tmp/pytest" \
  tests/test_embryo_r78_launcher.py tests/test_embryo_r78_build.py tests/test_embryo_r75_launcher.py \
  tests/test_embryo_r75_build.py tests/test_mix_lean_build.py tests/test_embryo_ot23_launcher.py \
  tests/test_e58_mix_lean_read.py tests/test_series_back_build.py tests/test_embryo_ovsb_launcher.py \
  tests/test_embryo_so_ohi_launcher.py tests/test_bf_build.py tests/test_bf_launchers.py tests/test_bf_form.py \
  tests/test_bf_read.py tests/test_embryo_overlay.py tests/test_shape_interp.py tests/test_shape_interp_build.py \
  tests/test_e47_floor_exact.py tests/test_t2_build.py tests/test_hidden_counts.py tests/test_s2_rules.py \
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
build_base() {  # folder name mode args: the build on the file's own path with its t_mix put back; exit 12
  local d="$O/$1" name=$2 mode=$3
  BUILT+=("$1/$2.h5ad")
  shift 3
  mkdir -p "$d"
  $PY experiments/builds/embryo_r78_check.py "$mode" embryo "$@" --out "$d/$name.h5ad" > /dev/null 2>&1 || return 12
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
if [ "$R75" = yes ]; then
  identity ovt56 T2embryo_s2x_ovt56 "$OT56_SHA" 25 $OT56_ARGS
  identity ovt512 T2embryo_s2x_ovt512 "$OT512_SHA" 26 $OT512_ARGS
fi

STEP=${STEM}half
begin
build_base "$STEP" "${NAME}_half" half $ARGS || finish $?
same "$STEP" "${NAME}_half" "$O_SHA" 51 T2embryo_s2x_ov
step_ok

for STEP in "${STEM}a" "${STEM}b"; do begin; build "$STEP" "$NAME" $ARGS || finish $?; step_ok; done

STEP=done
begin
(cd "$O" && sha256sum "${BUILT[@]}") > "$S/SHA256SUMS" || finish 40
cp "$S/SHA256SUMS" "$O/SHA256SUMS" || finish 40
cmp -s "$O/${STEM}a/$FILE" "$O/${STEM}b/$FILE" || finish 41
new=$(sha256sum "$O/${STEM}a/$FILE" | cut -d' ' -f1)
for var in $KNOWN; do
  [ "$new" != "${!var}" ] || { echo "stopped: $FILE has the bytes of a known file ($var)" >&2; finish 42; }
done
step_ok
