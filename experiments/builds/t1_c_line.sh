#!/usr/bin/env bash
# The third T1 search line of row R22 (experiments/harness/TRIES.md; the captain's order of 2026-11-22; the line added
# in docs/PLAN.md decision log 2026-11-22): file C's command (R13's recomputed-groups file with g_mom's damping at 1.0)
# with the outside atlas's direction step at 2.0 and at 2.5 times the calendar fraction, the damping held at 1.0,
# nothing else changed (default seed). Test-board builds on the DGX campaign checkout that holds the relabelling's
# output (~/v/vec-campaign-d, never the reference checkout). Run from the root of the checkout or detached worktree it
# is called in, and wait on the exit files with a background poller; never a foreground sleep:
#   bash experiments/builds/t1_c_line.sh <COMMIT> <RUN>
# RUN names this run's own scratch folder ~/v/scratch/t1cline/<RUN> and its output folder runs/t1cline_<commit7>_<RUN>.
# Neither may exist: a run folder is never reused and nothing is ever deleted here, so a second run takes a new name.
# In order, each step first writing the word running to ~/v/scratch/t1cline/<RUN>/<STEP>.exit, then its exit code:
#   gate   HEAD is COMMIT (the full hash), the tracked tree is clean, docs/S2_RULES.sha256 checks, TRIES.md holds one
#          row R22 whose first cell names `--g-mom-damping 1.0 --feng-fraction-mult 2.0`,
#          `--g-mom-damping 1.0 --feng-fraction-mult 2.5` and `--captain-ruling R22`, and the relabelling's output
#          folder data/derived/feng_relabel_c exists; then it writes the head file
#   tests  the builder's and this launcher's targeted tests (synthetic stages only), 8 BLAS threads
#   D      the default command (R13's file, set t1rc) once; it must give R13_SHA or everything stops
#   A      the --g-mom-damping 1.0 command (candidate A) once; it must give A_SHA or everything stops
#   B      the --feng-fraction-mult 1.5 command (candidate B) once; it must give B_SHA or everything stops
#   C      both (candidate C, this line's base) once; it must give C_SHA or everything stops
#   X2     the 2.0 command with --captain-ruling R20 (row R20's file) once; it must give X2_SHA or everything stops
#   Cm20a Cm20b  C's command with the multiplier at 2.0 and --captain-ruling R22, twice, into separate folders
#   Cm25a Cm25b  C's command with the multiplier at 2.5 and --captain-ruling R22, twice, into separate folders
#   done   one sha256 list over every built file, each pair compared byte for byte
# A build step also fails (15) unless its stamp says this commit and clean True.
# Builds run at 16 BLAS threads, as the reference builds (T1 files move with the thread count), seed and board as the
# reference commands (the builder's defaults: seed 1000001, board T1:test). Every builder and `vec check` output goes
# to /dev/null: a test build is read by its exit code, its sha256 and its stamp's commit and clean fields only
# (AGENTS.md), never its log, size, cell count or sidecar record. Nothing here uploads anything.
# Exit codes: 90 odd RUN, 91 run folder exists, 96 HEAD is not COMMIT, 98 tracked tree not clean, 97 S2_RULES.sha256,
# 95 not exactly one row R22 in TRIES.md or its first cell does not name both combined options and the ruling, 92 no
# relabelling output, 93 output folder exists, 10 tests, 11 build, 13 vec check, 14 no stamp, 15 stamp not this commit
# and clean, 20 D is not R13's bytes, 21 A is not candidate A's bytes, 22 B is not candidate B's bytes, 23 C is not
# candidate C's bytes, 24 X2 is not the 2.0 file's bytes, 30 hash list, 31 Cm20a and Cm20b differ, 32 Cm25a and Cm25b
# differ.
set -u
COMMIT=${1:?usage: t1_c_line.sh COMMIT RUN}
RUN=${2:?usage: t1_c_line.sh COMMIT RUN}
case "$RUN" in
  *[!A-Za-z0-9_-]*) echo "refused: RUN takes letters, digits, _ and - only" >&2; exit 90 ;;
esac
export PATH=$HOME/.local/bin:$PATH
S=$HOME/v/scratch/t1cline/$RUN
mkdir -p "$HOME/v/scratch/t1cline"
mkdir "$S" 2> /dev/null || { echo "refused: the run folder $S exists (a run folder is never reused)" >&2; exit 91; }
mkdir "$S/tmp"
export NUMPY_MADVISE_HUGEPAGE=0 TMPDIR="$S/tmp"
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
PY="uv run --frozen --no-sync python -I -B"
BOARD=T1:test
R13_SHA=ad8f24c71c90ca17bce5f8876d503471ab3faa5b17a0878addb00aff689fca4f  # set t1rc, row R13 (the default command)
A_SHA=76865826f5a3741c7585f81c4241048ca1a95b5791ebb3c8feb1f2b536678ce4  # candidate A, set t1cand_1420c70
B_SHA=c7240cb6a7f7ee93d037f437832efbaad004f1fbe608d80072f8efd5d20a0062  # candidate B, set t1cand_1420c70
C_SHA=31477ef7df5d5374a74a8249c570217020c593916ab4d3780049b7b3f87b0620  # candidate C, set t1cand_1420c70
X2_SHA=52897fe375b8a81afffa788592047c30361824148b00a8176514385060518837  # row R20's file, t1x2_3d18841.SHA256SUMS
BASE_ARGS="--g-mom --base lib3 --feng-groups recomputed"
A_ARGS="--g-mom --base lib3 --feng-groups recomputed --g-mom-damping 1.0"
B_ARGS="--g-mom --base lib3 --feng-groups recomputed --feng-fraction-mult 1.5"
C_ARGS="--g-mom --base lib3 --feng-groups recomputed --g-mom-damping 1.0 --feng-fraction-mult 1.5"
X2_ARGS="--g-mom --base lib3 --feng-groups recomputed --feng-fraction-mult 2.0 --captain-ruling R20"
CM20_ARGS="--g-mom --base lib3 --feng-groups recomputed --g-mom-damping 1.0 --feng-fraction-mult 2.0 --captain-ruling R22"
CM25_ARGS="--g-mom --base lib3 --feng-groups recomputed --g-mom-damping 1.0 --feng-fraction-mult 2.5 --captain-ruling R22"
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
[ "$(grep -c '^| R22 ' experiments/harness/TRIES.md)" = 1 ] \
  || { echo "refused: TRIES.md does not hold exactly one row R22" >&2; finish 95; }
first=$(grep '^| R22 ' experiments/harness/TRIES.md | cut -d'|' -f2)
for opt in '`--g-mom-damping 1.0 --feng-fraction-mult 2.0`' '`--g-mom-damping 1.0 --feng-fraction-mult 2.5`' \
  '`--captain-ruling R22`'; do
  printf '%s\n' "$first" | grep -qF -- "$opt" \
    || { echo "refused: row R22's first cell does not name $opt" >&2; finish 95; }
done
[ -d data/derived/feng_relabel_c ] \
  || { echo "refused: the relabelling's output data/derived/feng_relabel_c is not in this checkout" >&2; finish 92; }
O=runs/t1cline_${head:0:7}_$RUN
mkdir -p runs
mkdir "$O" 2> /dev/null || { echo "refused: the output folder $O exists" >&2; finish 93; }
echo "$head" > "$S/head"
step_ok

STEP=tests
begin
$PY -m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp="$S/tmp/pytest" \
  tests/test_t1_c_line_launcher.py tests/test_t1_line3_build.py tests/test_t1_damping_line_build.py \
  tests/test_t1_damping_line_launcher.py tests/test_t1_fraction_line_build.py tests/test_t1_fraction_line_launcher.py \
  tests/test_t1_fraction2_build.py tests/test_t1_fraction2_launcher.py tests/test_t1_candidates_build.py \
  tests/test_t1_bets.py tests/test_hidden_counts.py tests/test_s2_rules.py > "$S/tests.out" 2>&1 || finish 10
step_ok

export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16
stamp() {  # the stamp's commit and clean fields only
  $PY -c 'import json,sys; s=json.load(open(sys.argv[1])); print(s["commit"], s["clean"])' "$1.stamp.json"
}
build() {  # folder name args: the build, then `vec check`; exit 11 build, 13 check, 14 no stamp, 15 stamp not clean
  local d="$O/$1" name=$2 st
  shift 2
  mkdir -p "$d"
  $PY experiments/builds/build_test.py t1 "$@" --out "$d/$name.h5ad" > /dev/null 2>&1 || return 11
  uv run --frozen --no-sync vec check "$d/$name.h5ad" --board "$BOARD" > /dev/null 2>&1 || return 13
  st=$(stamp "$d/$name.h5ad" 2> /dev/null) || return 14
  echo "$(basename "$d")/$name $st" >> "$S/stamps"
  [ "$st" = "$head True" ] || return 15
}
: > "$S/stamps"

STEP=D
begin
build D T1_s2rc $BASE_ARGS || finish $?
got=$(sha256sum "$O/D/T1_s2rc.h5ad" | cut -d' ' -f1)
[ "$got" = "$R13_SHA" ] || { echo "stopped: the default command does not give R13's file" >&2; finish 20; }
step_ok

STEP=A
begin
build A T1_s2rcA $A_ARGS || finish $?
got=$(sha256sum "$O/A/T1_s2rcA.h5ad" | cut -d' ' -f1)
[ "$got" = "$A_SHA" ] || { echo "stopped: the 1.0 command does not give candidate A's file" >&2; finish 21; }
step_ok

STEP=B
begin
build B T1_s2rcB $B_ARGS || finish $?
got=$(sha256sum "$O/B/T1_s2rcB.h5ad" | cut -d' ' -f1)
[ "$got" = "$B_SHA" ] || { echo "stopped: the 1.5 command does not give candidate B's file" >&2; finish 22; }
step_ok

STEP=C
begin
build C T1_s2rcC $C_ARGS || finish $?
got=$(sha256sum "$O/C/T1_s2rcC.h5ad" | cut -d' ' -f1)
[ "$got" = "$C_SHA" ] || { echo "stopped: the 1.0 and 1.5 command does not give candidate C's file" >&2; finish 23; }
step_ok

STEP=X2
begin
build X2 T1_s2rcB2 $X2_ARGS || finish $?
got=$(sha256sum "$O/X2/T1_s2rcB2.h5ad" | cut -d' ' -f1)
[ "$got" = "$X2_SHA" ] || { echo "stopped: the 2.0 command does not give row R20's file" >&2; finish 24; }
step_ok

for STEP in Cm20a Cm20b; do begin; build "$STEP" T1_s2rcCm20 $CM20_ARGS || finish $?; step_ok; done
for STEP in Cm25a Cm25b; do begin; build "$STEP" T1_s2rcCm25 $CM25_ARGS || finish $?; step_ok; done

STEP=done
begin
(cd "$O" && sha256sum D/*.h5ad A/*.h5ad B/*.h5ad C/*.h5ad X2/*.h5ad Cm20a/*.h5ad Cm20b/*.h5ad Cm25a/*.h5ad \
  Cm25b/*.h5ad) > "$S/SHA256SUMS" || finish 30
cp "$S/SHA256SUMS" "$O/SHA256SUMS" || finish 30
cmp -s "$O/Cm20a/T1_s2rcCm20.h5ad" "$O/Cm20b/T1_s2rcCm20.h5ad" || finish 31
cmp -s "$O/Cm25a/T1_s2rcCm25.h5ad" "$O/Cm25b/T1_s2rcCm25.h5ad" || finish 32
step_ok
