#!/usr/bin/env bash
# Blood's released momentum and MS on the banked command (the T1 design round of 2026-11-24, files F1, F7, F8 and F2 as
# its cold challenge fixes them), each one change to the banked file T1_s2rcB25 (t1 --g-mom --base lib3 --feng-groups
# recomputed --feng-fraction-mult 2.5), each only under a ruling row of the captain named on the command line:
#   T1_s2rcB25H05  the drawn Blood cells get Blood's own released E9.5 -> E10.5 change at damping 0.5 (fitted on
#                  Blood-labelled cells of E9.5 and E10.5's libraries 1 and 2), after the atlas step, before g_mom
#   T1_s2rcB25H10  the same at damping 1.0
#   T1_s2rcB25H20  the same at damping 2.0
#   T1_s2rcB25MS   candidate MS's step exactly as row R23 fixes it, on the banked command
# nothing else changed (default seed). Test-board builds on a DGX campaign worktree that holds the relabelling's output
# (data/derived/feng_relabel_c) and its own venv, made there by `uv sync --frozen` (never the reference checkout). Run
# from the root of that worktree, detached, and wait on the exit files with a background poller; never a foreground
# sleep:
#   bash experiments/builds/t1_hm_line.sh <COMMIT> <RUN> <ROW_H> <ROW_M>
# ROW_H is the row of experiments/harness/TRIES.md whose first cell is the captain's ruling naming `--held-momentum 0.5`,
# `--held-momentum 1.0`, `--held-momentum 2.0` and `--captain-ruling <ROW_H>`; ROW_M the row whose first cell is his
# ruling naming `--feng-fraction-mult 2.5 --direction-momentum all` and `--captain-ruling <ROW_M>` (the two may be one
# row). RUN names this run's own scratch folder ~/v/scratch/t1hm/<RUN> and its output folder runs/t1hm_<commit7>_<RUN>.
# Neither may exist: a run folder is never reused and nothing is ever deleted here, so a second run takes a new name.
# In order, each step first writing the word running to ~/v/scratch/t1hm/<RUN>/<STEP>.exit, then its exit code:
#   gate   HEAD is COMMIT (the full hash), the tracked tree is clean, docs/S2_RULES.sha256 checks, TRIES.md holds exactly
#          one row ROW_H whose first cell starts "<ROW_H> (the captain's ", says "ruling, no run" and names the three
#          held-momentum options and `--captain-ruling <ROW_H>`, and exactly one row ROW_M whose first cell starts
#          "<ROW_M> (the captain's ", says "ruling, no run" and names `--feng-fraction-mult 2.5 --direction-momentum
#          all` and `--captain-ruling <ROW_M>` (the builder's own check), the relabelling's output folder
#          data/derived/feng_relabel_c exists, and .venv's vec.pth names this checkout's src; then it writes the head
#          file
#   tests  the builders' and launchers' targeted tests (synthetic stages only), 8 BLAS threads
#   D ... B2d20  the fifteen identity rebuilds of t1_l67pb_line.sh, each once:
#   D      the default command (R13's file); it must give R13_SHA or everything stops
#   A      the --g-mom-damping 1.0 command (candidate A); A_SHA
#   B      the --feng-fraction-mult 1.5 command (candidate B); B_SHA
#   C      both (candidate C); C_SHA
#   X2     the 2.0 command with --captain-ruling R20 (row R20's file); X2_SHA
#   X25    the 2.5 command with --captain-ruling R22 (T1_s2rcB25, the banked file, the base of every new file); X25_SHA
#   A15    the --g-mom-damping 1.5 command with --captain-ruling R22 (T1_s2rcA15); A15_SHA
#   A20    the --g-mom-damping 2.0 command with --captain-ruling R22 (T1_s2rcA20); A20_SHA
#   Cm20   the --g-mom-damping 1.0 --feng-fraction-mult 2.0 command with --captain-ruling R22 (T1_s2rcCm20); CM20_SHA
#   MS     the --direction-momentum all command with --captain-ruling R23 (candidate MS); MS_SHA
#   MSB    the --feng-fraction-mult 1.5 --direction-momentum all command with --captain-ruling R26 (T1_s2rcMSB); MSB_SHA
#   X30    the 3.0 command with --captain-ruling R22 (T1_s2rcB30); X30_SHA
#   Cm25   the --g-mom-damping 1.0 --feng-fraction-mult 2.5 command with --captain-ruling R22 (T1_s2rcCm25); CM25_SHA
#   B2d15  the --g-mom-damping 1.5 --feng-fraction-mult 2.0 command with --captain-ruling R22 (T1_s2rcB2d15);
#          B2D15_SHA
#   B2d20  the --g-mom-damping 2.0 --feng-fraction-mult 2.0 command with --captain-ruling R22 (T1_s2rcB2d20);
#          B2D20_SHA (any of these fifteen giving other bytes stops everything)
#   H05one H10one H20one MS25one  each new command once through experiments/builds/t1_hm_check.py one, its momentum's
#          factors forced to 1 (every momentum_shift call's damping 0): each must give X25_SHA, the banked file's
#          bytes, or everything stops
#   H05a H05b    T1_s2rcB25H05 by its command, twice, into separate folders
#   H10a H10b    T1_s2rcB25H10, twice
#   H20a H20b    T1_s2rcB25H20, twice
#   MS25a MS25b  T1_s2rcB25MS, twice
#   rows   t1_hm_check.py rows on X25's file and H05a, H10a and H20a: the banked command's draw (the builder's own calls,
#          stopped before anything is written) puts the drawn Blood cells on rows the banked file holds at their drawn
#          counts, every other row of each H file is the banked file's byte for byte, and a Blood row moved; booleans
#          only, to rows.out
#   record t1_hm_check.py record on H05a, H10a and H20a: only each build record's held-momentum fields are read, by
#          code (the fit from E9.5 to E10.5's libraries 1 and 2 on Blood-labelled cells, library 3 out, Blood shifted,
#          its damping, stream 6); booleans only, to record.out
#   done   one sha256 list over every built file, each pair compared byte for byte
# A build step also fails (15) unless its stamp says this commit and clean True.
# Builds run at 16 BLAS threads, as the reference builds (T1 files move with the thread count), seed and board as the
# reference commands (the builder's defaults: seed 1000001, board T1:test). Every builder, checker build and `vec check`
# output goes to /dev/null: a test build is read by its exit code, its sha256 and its stamp's commit and clean fields
# (AGENTS.md), and by the two checks' booleans, never its log, size or cell count. Nothing here uploads anything.
# Exit codes: 90 odd RUN, 89 odd row ID, 91 run folder exists, 96 HEAD is not COMMIT, 98 tracked tree not clean,
# 97 S2_RULES.sha256, 95 row ROW_H not as above, 94 row ROW_M not as above, 92 no relabelling output, 88 .venv not this
# checkout's, 93 output folder exists, 10 tests, 11 build, 12 factors-at-1 build, 13 vec check, 14 no stamp, 15 stamp
# not this commit and clean, 20 D is not R13's bytes, 21 A is not candidate A's bytes, 22 B is not candidate B's bytes,
# 23 C is not candidate C's bytes, 24 X2 is not the 2.0 file's bytes, 25 X25 is not T1_s2rcB25's bytes, 26 A15 is not
# T1_s2rcA15's bytes, 27 A20 is not T1_s2rcA20's bytes, 28 Cm20 is not T1_s2rcCm20's bytes, 29 MS is not candidate
# MS's bytes, 30 MSB is not T1_s2rcMSB's bytes, 31 X30 is not T1_s2rcB30's bytes, 32 Cm25 is not T1_s2rcCm25's bytes,
# 33 B2d15 is not T1_s2rcB2d15's bytes, 34 B2d20 is not T1_s2rcB2d20's bytes, 35 H05one is not the banked file's
# bytes, 36 H10one is not, 37 H20one is not, 38 MS25one is not, 40 hash list, 41 H05a and H05b differ, 42 H10a and
# H10b differ, 43 H20a and H20b differ, 44 MS25a and MS25b differ, 45 rows check, 46 record check.
set -u
COMMIT=${1:?usage: t1_hm_line.sh COMMIT RUN ROW_H ROW_M}
RUN=${2:?usage: t1_hm_line.sh COMMIT RUN ROW_H ROW_M}
ROW_H=${3:?usage: t1_hm_line.sh COMMIT RUN ROW_H ROW_M}
ROW_M=${4:?usage: t1_hm_line.sh COMMIT RUN ROW_H ROW_M}
case "$RUN" in
  *[!A-Za-z0-9_-]*) echo "refused: RUN takes letters, digits, _ and - only" >&2; exit 90 ;;
esac
for id in "$ROW_H" "$ROW_M"; do
  [[ "$id" =~ ^R[0-9]+[a-z]?$ ]] || { echo "refused: a row ID is R, digits and an optional small letter" >&2; exit 89; }
done
export PATH=$HOME/.local/bin:$PATH
unset UV_PROJECT_ENVIRONMENT  # uv runs this checkout's own .venv, the one the gate checks
S=$HOME/v/scratch/t1hm/$RUN
mkdir -p "$HOME/v/scratch/t1hm"
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
X25_SHA=ccb7c457c0e8ac8ac668bcea17bc61ed22f414f58ab5a63e946cb2f08c4f4e77  # T1_s2rcB25, t1line_e3cc0e4_r1.SHA256SUMS
A15_SHA=e2334137b54d98b04e9dece90213987e61fe9aa15a3fb8d85e6c262d8fb5ba75  # T1_s2rcA15, t1damp_9a10157_r1.SHA256SUMS
A20_SHA=372713686782be273a52427140eab652517e35b679f54fc2a969484dee730218  # T1_s2rcA20, t1damp_9a10157_r1.SHA256SUMS
CM20_SHA=acce5f5f7097269bffebbfd583833e4c9c1c5ae77f1d5294dde6d13e197f863e  # T1_s2rcCm20, t1cline_2a1c248_c1.SHA256SUMS
MS_SHA=6103d10535665634ee0cad3390590a5df132d1805682fb1d309bb19581fdaa5e  # candidate MS, t1ms_6bc36d9_r1.SHA256SUMS
MSB_SHA=78bb17eea5721e590cfd2c9f707e9cd6f2f496e081acede4f697dc4035871da2  # T1_s2rcMSB, t1msc_a3e879a_m1.SHA256SUMS
X30_SHA=b639de29bd13d7005450ca53e7cd3dacb2e362a55beca8545937858e6582a024  # T1_s2rcB30, t1line_e3cc0e4_r1.SHA256SUMS
CM25_SHA=c5fb96a64075b6fbfb2f3066604d622ea8b89ac3a5e4f8a99d0fea7d93011dfd  # T1_s2rcCm25, t1cline_2a1c248_c1.SHA256SUMS
B2D15_SHA=fe4e7403d6257cc265a02d0a13d40cdb60a97e9fdc0647d6b5c7e318aedaed6c  # T1_s2rcB2d15, t1b2dline_7e7d9a6_b1
B2D20_SHA=bd12c66312ba32ab938b4ee71bd5c79e58a154d6fc56434012e6dad907e668dd  # T1_s2rcB2d20, t1b2dline_7e7d9a6_b1
BASE_ARGS="--g-mom --base lib3 --feng-groups recomputed"
A_ARGS="--g-mom --base lib3 --feng-groups recomputed --g-mom-damping 1.0"
B_ARGS="--g-mom --base lib3 --feng-groups recomputed --feng-fraction-mult 1.5"
C_ARGS="--g-mom --base lib3 --feng-groups recomputed --g-mom-damping 1.0 --feng-fraction-mult 1.5"
X2_ARGS="--g-mom --base lib3 --feng-groups recomputed --feng-fraction-mult 2.0 --captain-ruling R20"
X25_ARGS="--g-mom --base lib3 --feng-groups recomputed --feng-fraction-mult 2.5 --captain-ruling R22"
A15_ARGS="--g-mom --base lib3 --feng-groups recomputed --g-mom-damping 1.5 --captain-ruling R22"
A20_ARGS="--g-mom --base lib3 --feng-groups recomputed --g-mom-damping 2.0 --captain-ruling R22"
CM20_ARGS="--g-mom --base lib3 --feng-groups recomputed --g-mom-damping 1.0 --feng-fraction-mult 2.0 --captain-ruling R22"
MS_ARGS="--g-mom --base lib3 --feng-groups recomputed --direction-momentum all --captain-ruling R23"
MSB_ARGS="--g-mom --base lib3 --feng-groups recomputed --feng-fraction-mult 1.5 --direction-momentum all --captain-ruling R26"
X30_ARGS="--g-mom --base lib3 --feng-groups recomputed --feng-fraction-mult 3.0 --captain-ruling R22"
CM25_ARGS="--g-mom --base lib3 --feng-groups recomputed --g-mom-damping 1.0 --feng-fraction-mult 2.5 --captain-ruling R22"
B2D15_ARGS="--g-mom --base lib3 --feng-groups recomputed --g-mom-damping 1.5 --feng-fraction-mult 2.0 --captain-ruling R22"
B2D20_ARGS="--g-mom --base lib3 --feng-groups recomputed --g-mom-damping 2.0 --feng-fraction-mult 2.0 --captain-ruling R22"
H05_ARGS="--g-mom --base lib3 --feng-groups recomputed --feng-fraction-mult 2.5 --held-momentum 0.5 --captain-ruling $ROW_H"
H10_ARGS="--g-mom --base lib3 --feng-groups recomputed --feng-fraction-mult 2.5 --held-momentum 1.0 --captain-ruling $ROW_H"
H20_ARGS="--g-mom --base lib3 --feng-groups recomputed --feng-fraction-mult 2.5 --held-momentum 2.0 --captain-ruling $ROW_H"
MS25_ARGS="--g-mom --base lib3 --feng-groups recomputed --feng-fraction-mult 2.5 --direction-momentum all --captain-ruling $ROW_M"
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
ruling_row() {  # ID CODE OPTION...: one row ID, its first cell the captain's ruling naming each OPTION, else CODE
  local id=$1 code=$2 first opt
  shift 2
  [ "$(grep -c "^| $id " experiments/harness/TRIES.md)" = 1 ] \
    || { echo "refused: TRIES.md does not hold exactly one row $id" >&2; finish "$code"; }
  first=$(grep "^| $id " experiments/harness/TRIES.md | cut -d'|' -f2)
  printf '%s\n' "$first" | grep -q "^ $id (the captain's " \
    || { echo "refused: row $id's first cell is not the captain's ruling" >&2; finish "$code"; }
  printf '%s\n' "$first" | grep -qF -- "ruling, no run" \
    || { echo "refused: row $id's first cell does not say ruling, no run" >&2; finish "$code"; }
  for opt in "$@" "\`--captain-ruling $id\`"; do
    printf '%s\n' "$first" | grep -qF -- "$opt" \
      || { echo "refused: row $id's first cell does not name $opt" >&2; finish "$code"; }
  done
}
ruling_row "$ROW_H" 95 '`--held-momentum 0.5`' '`--held-momentum 1.0`' '`--held-momentum 2.0`'
ruling_row "$ROW_M" 94 '`--feng-fraction-mult 2.5 --direction-momentum all`'
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
O=runs/t1hm_${head:0:7}_$RUN
mkdir -p runs
mkdir "$O" 2> /dev/null || { echo "refused: the output folder $O exists" >&2; finish 93; }
echo "$head" > "$S/head"
step_ok

STEP=tests
begin
$PY -m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp="$S/tmp/pytest" \
  tests/test_t1_hm_launcher.py tests/test_t1_hm_build.py tests/test_t1_l67pb_launcher.py tests/test_t1_l67pb_build.py \
  tests/test_t1_b2d_line_launcher.py tests/test_t1_line5_build.py tests/test_t1_line3_build.py \
  tests/test_t1_c_line_launcher.py tests/test_t1_damping_line_build.py tests/test_t1_damping_line_launcher.py \
  tests/test_t1_fraction_line_build.py tests/test_t1_fraction_line_launcher.py tests/test_t1_fraction2_build.py \
  tests/test_t1_fraction2_launcher.py tests/test_t1_ms_build.py tests/test_t1_msb_build.py \
  tests/test_t1_candidates_build.py tests/test_t1_bets.py tests/test_hidden_counts.py tests/test_s2_rules.py \
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
build() {  # folder name args: the build, then checked; exit 11 build
  local d="$O/$1" name=$2
  shift 2
  mkdir -p "$d"
  $PY experiments/builds/build_test.py t1 "$@" --out "$d/$name.h5ad" > /dev/null 2>&1 || return 11
  checked "$d" "$name"
}
build_one() {  # folder name args: the build with the momentum's factors forced to 1, then checked; exit 12 build
  local d="$O/$1" name=$2
  shift 2
  mkdir -p "$d"
  $PY experiments/builds/t1_hm_check.py one t1 "$@" --out "$d/$name.h5ad" > /dev/null 2>&1 || return 12
  checked "$d" "$name"
}
same() {  # folder name sha code: stop with the code unless the built file has the known bytes
  local got
  got=$(sha256sum "$O/$1/$2.h5ad" | cut -d' ' -f1)
  [ "$got" = "$3" ] || { echo "stopped: step $1 does not give the known bytes of $2" >&2; finish "$4"; }
}
: > "$S/stamps"

STEP=D
begin
build D T1_s2rc $BASE_ARGS || finish $?
same D T1_s2rc "$R13_SHA" 20
step_ok

STEP=A
begin
build A T1_s2rcA $A_ARGS || finish $?
same A T1_s2rcA "$A_SHA" 21
step_ok

STEP=B
begin
build B T1_s2rcB $B_ARGS || finish $?
same B T1_s2rcB "$B_SHA" 22
step_ok

STEP=C
begin
build C T1_s2rcC $C_ARGS || finish $?
same C T1_s2rcC "$C_SHA" 23
step_ok

STEP=X2
begin
build X2 T1_s2rcB2 $X2_ARGS || finish $?
same X2 T1_s2rcB2 "$X2_SHA" 24
step_ok

STEP=X25
begin
build X25 T1_s2rcB25 $X25_ARGS || finish $?
same X25 T1_s2rcB25 "$X25_SHA" 25
step_ok

STEP=A15
begin
build A15 T1_s2rcA15 $A15_ARGS || finish $?
same A15 T1_s2rcA15 "$A15_SHA" 26
step_ok

STEP=A20
begin
build A20 T1_s2rcA20 $A20_ARGS || finish $?
same A20 T1_s2rcA20 "$A20_SHA" 27
step_ok

STEP=Cm20
begin
build Cm20 T1_s2rcCm20 $CM20_ARGS || finish $?
same Cm20 T1_s2rcCm20 "$CM20_SHA" 28
step_ok

STEP=MS
begin
build MS T1_s2rc_MS $MS_ARGS || finish $?
same MS T1_s2rc_MS "$MS_SHA" 29
step_ok

STEP=MSB
begin
build MSB T1_s2rcMSB $MSB_ARGS || finish $?
same MSB T1_s2rcMSB "$MSB_SHA" 30
step_ok

STEP=X30
begin
build X30 T1_s2rcB30 $X30_ARGS || finish $?
same X30 T1_s2rcB30 "$X30_SHA" 31
step_ok

STEP=Cm25
begin
build Cm25 T1_s2rcCm25 $CM25_ARGS || finish $?
same Cm25 T1_s2rcCm25 "$CM25_SHA" 32
step_ok

STEP=B2d15
begin
build B2d15 T1_s2rcB2d15 $B2D15_ARGS || finish $?
same B2d15 T1_s2rcB2d15 "$B2D15_SHA" 33
step_ok

STEP=B2d20
begin
build B2d20 T1_s2rcB2d20 $B2D20_ARGS || finish $?
same B2d20 T1_s2rcB2d20 "$B2D20_SHA" 34
step_ok

STEP=H05one
begin
build_one H05one T1_s2rcB25H05_one $H05_ARGS || finish $?
same H05one T1_s2rcB25H05_one "$X25_SHA" 35
step_ok

STEP=H10one
begin
build_one H10one T1_s2rcB25H10_one $H10_ARGS || finish $?
same H10one T1_s2rcB25H10_one "$X25_SHA" 36
step_ok

STEP=H20one
begin
build_one H20one T1_s2rcB25H20_one $H20_ARGS || finish $?
same H20one T1_s2rcB25H20_one "$X25_SHA" 37
step_ok

STEP=MS25one
begin
build_one MS25one T1_s2rcB25MS_one $MS25_ARGS || finish $?
same MS25one T1_s2rcB25MS_one "$X25_SHA" 38
step_ok

for STEP in H05a H05b; do begin; build "$STEP" T1_s2rcB25H05 $H05_ARGS || finish $?; step_ok; done
for STEP in H10a H10b; do begin; build "$STEP" T1_s2rcB25H10 $H10_ARGS || finish $?; step_ok; done
for STEP in H20a H20b; do begin; build "$STEP" T1_s2rcB25H20 $H20_ARGS || finish $?; step_ok; done
for STEP in MS25a MS25b; do begin; build "$STEP" T1_s2rcB25MS $MS25_ARGS || finish $?; step_ok; done

STEP=rows
begin
$PY experiments/builds/t1_hm_check.py rows "$O/X25/T1_s2rcB25.h5ad" "$O/H05a/T1_s2rcB25H05.h5ad" \
  "$O/H10a/T1_s2rcB25H10.h5ad" "$O/H20a/T1_s2rcB25H20.h5ad" > "$S/rows.out" 2> /dev/null || finish 45
step_ok

STEP=record
begin
$PY experiments/builds/t1_hm_check.py record "$O/H05a/T1_s2rcB25H05.h5ad" 0.5 "$O/H10a/T1_s2rcB25H10.h5ad" 1.0 \
  "$O/H20a/T1_s2rcB25H20.h5ad" 2.0 > "$S/record.out" 2> /dev/null || finish 46
step_ok

STEP=done
begin
(cd "$O" && sha256sum D/*.h5ad A/*.h5ad B/*.h5ad C/*.h5ad X2/*.h5ad X25/*.h5ad A15/*.h5ad A20/*.h5ad Cm20/*.h5ad \
  MS/*.h5ad MSB/*.h5ad X30/*.h5ad Cm25/*.h5ad B2d15/*.h5ad B2d20/*.h5ad H05one/*.h5ad H10one/*.h5ad H20one/*.h5ad \
  MS25one/*.h5ad H05a/*.h5ad H05b/*.h5ad H10a/*.h5ad H10b/*.h5ad H20a/*.h5ad H20b/*.h5ad MS25a/*.h5ad \
  MS25b/*.h5ad) > "$S/SHA256SUMS" || finish 40
cp "$S/SHA256SUMS" "$O/SHA256SUMS" || finish 40
cmp -s "$O/H05a/T1_s2rcB25H05.h5ad" "$O/H05b/T1_s2rcB25H05.h5ad" || finish 41
cmp -s "$O/H10a/T1_s2rcB25H10.h5ad" "$O/H10b/T1_s2rcB25H10.h5ad" || finish 42
cmp -s "$O/H20a/T1_s2rcB25H20.h5ad" "$O/H20b/T1_s2rcB25H20.h5ad" || finish 43
cmp -s "$O/MS25a/T1_s2rcB25MS.h5ad" "$O/MS25b/T1_s2rcB25MS.h5ad" || finish 44
step_ok
