#!/usr/bin/env bash
# Row R51's seven files (T1's line 10, the captain's order of 2026-11-24), each one change to T1_s2rcMSB's command
# (t1 --g-mom --base lib3 --feng-groups recomputed --feng-fraction-mult 1.5 --direction-momentum all, row R26's), each
# only with --captain-ruling R51:
#   T1_s2rcMS20    --feng-fraction-mult 2.0 --direction-momentum all (MS's step at the atlas multiplier 2.0)
#   T1_s2rcMS30    --feng-fraction-mult 3.0 --direction-momentum all (at 3.0)
#   T1_s2rcMSBd10  --feng-fraction-mult 1.5 --direction-momentum all --g-mom-damping 1.0 (g_mom's damping 1.0)
#   T1_s2rcMSBd15  --feng-fraction-mult 1.5 --direction-momentum all --g-mom-damping 1.5 (at 1.5)
#   T1_s2rcMSBH05  --feng-fraction-mult 1.5 --direction-momentum all --held-momentum 0.5 (Blood's own released change
#                  at damping 0.5, row R44's step, after MS's step and before g_mom)
#   T1_s2rcMSBH10  --feng-fraction-mult 1.5 --direction-momentum all --held-momentum 1.0 (at 1.0)
#   T1_s2rcMSBm10  --feng-fraction-mult 1.5 --direction-momentum all --direction-momentum-damping 1.0 (MS's own damping)
# nothing else changed (default seed). Test-board builds on a DGX campaign worktree that holds the relabelling's output
# (data/derived/feng_relabel_c) and its own venv, made there by `uv sync --frozen` (never the reference checkout). Run
# from the root of that worktree, detached, and wait on the exit files with a background poller; never a foreground
# sleep:
#   bash experiments/builds/t1_r51_line.sh <COMMIT> <RUN> R51
# The third argument is the ruling row's ID and must be R51, the one row the builder takes for these files. RUN names
# this run's own scratch folder ~/v/scratch/t1r51/<RUN> and its output folder runs/t1r51_<commit7>_<RUN>. Neither may
# exist: a run folder is never reused and nothing is ever deleted here, so a second run takes a new name.
# In order, each step first writing the word running to ~/v/scratch/t1r51/<RUN>/<STEP>.exit, then its exit code:
#   gate   HEAD is COMMIT (the full hash), the tracked tree is clean, docs/S2_RULES.sha256 checks, TRIES.md holds exactly
#          one row R51 whose first cell starts "R51 (the captain's ", says "ruling, no run" and names each of the seven
#          option strings and `--captain-ruling R51`, each as a whole backticked string (the builder's own check: a
#          backticked string holding an option inside a longer one does not name it), the relabelling's output folder
#          data/derived/feng_relabel_c exists, and .venv's vec.pth names this checkout's src; then it writes the head
#          file
#   tests  the builders', checkers' and launchers' targeted tests (synthetic stages only), 8 BLAS threads
#   D ... pooled  the nineteen identity rebuilds, each once: the fifteen of t1_hm_line.sh and t1_l67pb_line.sh and the
#          four files t1_l67pb_line.sh built twice with one sha256 (t1l67pb_ce1fb91_l1.SHA256SUMS); any of them giving
#          other bytes stops everything:
#   D      the default command (R13's file); R13_SHA
#   A      the --g-mom-damping 1.0 command (candidate A); A_SHA
#   B      the --feng-fraction-mult 1.5 command (candidate B); B_SHA
#   C      both (candidate C); C_SHA
#   X2     the 2.0 command with --captain-ruling R20 (row R20's file); X2_SHA
#   X25    the 2.5 command with --captain-ruling R22 (T1_s2rcB25); X25_SHA
#   A15    the --g-mom-damping 1.5 command with --captain-ruling R22 (T1_s2rcA15); A15_SHA
#   A20    the --g-mom-damping 2.0 command with --captain-ruling R22 (T1_s2rcA20); A20_SHA
#   Cm20   the --g-mom-damping 1.0 --feng-fraction-mult 2.0 command with --captain-ruling R22 (T1_s2rcCm20); CM20_SHA
#   MS     the --direction-momentum all command with --captain-ruling R23 (candidate MS); MS_SHA
#   MSB    the --feng-fraction-mult 1.5 --direction-momentum all command with --captain-ruling R26 (T1_s2rcMSB, the base
#          of every new file); MSB_SHA
#   X30    the 3.0 command with --captain-ruling R22 (T1_s2rcB30); X30_SHA
#   Cm25   the --g-mom-damping 1.0 --feng-fraction-mult 2.5 command with --captain-ruling R22 (T1_s2rcCm25); CM25_SHA
#   B2d15  the --g-mom-damping 1.5 --feng-fraction-mult 2.0 command with --captain-ruling R22 (T1_s2rcB2d15);
#          B2D15_SHA
#   B2d20  the --g-mom-damping 2.0 --feng-fraction-mult 2.0 command with --captain-ruling R22 (T1_s2rcB2d20);
#          B2D20_SHA
#   Cm30   the --g-mom-damping 1.0 --feng-fraction-mult 3.0 command with --captain-ruling R22 (T1_s2rcCm30); CM30_SHA
#   B25d15 the --g-mom-damping 1.5 --feng-fraction-mult 2.5 command with --captain-ruling R22 (T1_s2rcB25d15);
#          B25D15_SHA
#   B25d20 the --g-mom-damping 2.0 --feng-fraction-mult 2.5 command with --captain-ruling R22 (T1_s2rcB25d20);
#          B25D20_SHA
#   pooled the --base pooled command with --captain-ruling R37 (T1_s2rc_pooled); POOLED_SHA
#   MS20msb MS30msb d10msb d15msb H05msb H10msb m10msb  each new command once through experiments/builds/t1_r51_check.py
#          msb, on its own path with its one setting at T1_s2rcMSB's value (the multiplier at 1.5, g_mom's damping at
#          0.5, Blood's factors at 1, MS's damping at the default); each must give MSB_SHA, T1_s2rcMSB's bytes, or
#          everything stops
#   MS20a MS20b ... m10a m10b  each new file by its command, twice, into separate folders
#   rows   t1_r51_check.py rows on MSB's file and H05a and H10a: MSB's draw (the builder's own calls, stopped before
#          anything is written) puts the drawn Blood cells on rows MSB's file holds at their drawn counts, every other
#          row of each H file is MSB's byte for byte at the same seed, and a Blood row moved; booleans only, to rows.out
#   done   one sha256 list over every built file, each pair compared byte for byte
# A build step also fails (15) unless its stamp says this commit and clean True.
# Builds run at 16 BLAS threads, as the reference builds (T1 files move with the thread count), seed and board as the
# reference commands (the builder's defaults: seed 1000001, board T1:test). Every builder, checker build and `vec check`
# output goes to /dev/null: a test build is read by its exit code, its sha256 and its stamp's commit and clean fields
# (AGENTS.md), and by the rows check's booleans, never its log, size, record or cell count. Nothing here uploads
# anything.
# Exit codes: 90 odd RUN, 89 row ID not R51, 91 run folder exists, 96 HEAD is not COMMIT, 98 tracked tree not clean,
# 97 S2_RULES.sha256, 95 row R51 not as above, 92 no relabelling output, 88 .venv not this checkout's, 93 output folder
# exists, 10 tests, 11 build, 12 checker build, 13 vec check, 14 no stamp, 15 stamp not this commit and clean, 20 D is
# not R13's bytes, 21 A is not candidate A's bytes, 22 B is not candidate B's bytes, 23 C is not candidate C's bytes,
# 24 X2 is not the 2.0 file's bytes, 25 X25 is not T1_s2rcB25's bytes, 26 A15 is not T1_s2rcA15's bytes, 27 A20 is not
# T1_s2rcA20's bytes, 28 Cm20 is not T1_s2rcCm20's bytes, 29 MS is not candidate MS's bytes, 30 MSB is not
# T1_s2rcMSB's bytes, 31 X30 is not T1_s2rcB30's bytes, 32 Cm25 is not T1_s2rcCm25's bytes, 33 B2d15 is not
# T1_s2rcB2d15's bytes, 34 B2d20 is not T1_s2rcB2d20's bytes, 35 Cm30 is not T1_s2rcCm30's bytes, 36 B25d15 is not
# T1_s2rcB25d15's bytes, 37 B25d20 is not T1_s2rcB25d20's bytes, 38 pooled is not T1_s2rc_pooled's bytes, 51 MS20msb
# is not MSB's bytes, 52 MS30msb is not, 53 d10msb is not, 54 d15msb is not, 55 H05msb is not, 56 H10msb is not,
# 57 m10msb is not, 40 hash list, 41 MS20a and MS20b differ, 42 MS30a and MS30b differ, 43 d10a and d10b differ,
# 44 d15a and d15b differ, 45 H05a and H05b differ, 46 H10a and H10b differ, 47 m10a and m10b differ, 48 rows check.
set -u
COMMIT=${1:?usage: t1_r51_line.sh COMMIT RUN R51}
RUN=${2:?usage: t1_r51_line.sh COMMIT RUN R51}
ROW=${3:?usage: t1_r51_line.sh COMMIT RUN R51}
case "$RUN" in
  *[!A-Za-z0-9_-]*) echo "refused: RUN takes letters, digits, _ and - only" >&2; exit 90 ;;
esac
[ "$ROW" = R51 ] || { echo "refused: the ruling row of these files is R51, and only R51" >&2; exit 89; }
export PATH=$HOME/.local/bin:$PATH
unset UV_PROJECT_ENVIRONMENT  # uv runs this checkout's own .venv, the one the gate checks
S=$HOME/v/scratch/t1r51/$RUN
mkdir -p "$HOME/v/scratch/t1r51"
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
CM30_SHA=807fe097c88ee779430896002a34635aeffbc99b433c68ebdb7bd41d8af2cbbe  # T1_s2rcCm30, t1l67pb_ce1fb91_l1
B25D15_SHA=5bc2b32e4688439f29c54c1b1203325595d996c37df30ba2b9ff26f27a60f4c6  # T1_s2rcB25d15, t1l67pb_ce1fb91_l1
B25D20_SHA=e6520dbee0fdc3f1336bc3f9de1e99e8a501de2d440b3f9085c88b5c7ecfb3f6  # T1_s2rcB25d20, t1l67pb_ce1fb91_l1
POOLED_SHA=91341a5f8c71dcc103a47f9c19a5d7794b1358c1c642c5c9251624eb6feb8bdc  # T1_s2rc_pooled, t1l67pb_ce1fb91_l1
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
CM30_ARGS="--g-mom --base lib3 --feng-groups recomputed --g-mom-damping 1.0 --feng-fraction-mult 3.0 --captain-ruling R22"
B25D15_ARGS="--g-mom --base lib3 --feng-groups recomputed --g-mom-damping 1.5 --feng-fraction-mult 2.5 --captain-ruling R22"
B25D20_ARGS="--g-mom --base lib3 --feng-groups recomputed --g-mom-damping 2.0 --feng-fraction-mult 2.5 --captain-ruling R22"
POOLED_ARGS="--g-mom --base pooled --feng-groups recomputed --captain-ruling R37"
MS20_ARGS="--g-mom --base lib3 --feng-groups recomputed --feng-fraction-mult 2.0 --direction-momentum all --captain-ruling R51"
MS30_ARGS="--g-mom --base lib3 --feng-groups recomputed --feng-fraction-mult 3.0 --direction-momentum all --captain-ruling R51"
D10_ARGS="--g-mom --base lib3 --feng-groups recomputed --feng-fraction-mult 1.5 --direction-momentum all --g-mom-damping 1.0 --captain-ruling R51"
D15_ARGS="--g-mom --base lib3 --feng-groups recomputed --feng-fraction-mult 1.5 --direction-momentum all --g-mom-damping 1.5 --captain-ruling R51"
H05_ARGS="--g-mom --base lib3 --feng-groups recomputed --feng-fraction-mult 1.5 --direction-momentum all --held-momentum 0.5 --captain-ruling R51"
H10_ARGS="--g-mom --base lib3 --feng-groups recomputed --feng-fraction-mult 1.5 --direction-momentum all --held-momentum 1.0 --captain-ruling R51"
M10_ARGS="--g-mom --base lib3 --feng-groups recomputed --feng-fraction-mult 1.5 --direction-momentum all --direction-momentum-damping 1.0 --captain-ruling R51"
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
ruling_row() {  # ID CODE OPTION...: one row ID, its first cell the captain's ruling naming each OPTION whole, else CODE
  local id=$1 code=$2 first opt
  shift 2
  [ "$(grep -c "^| $id " experiments/harness/TRIES.md)" = 1 ] \
    || { echo "refused: TRIES.md does not hold exactly one row $id" >&2; finish "$code"; }
  first=$(grep "^| $id " experiments/harness/TRIES.md | cut -d'|' -f2)
  printf '%s\n' "$first" | grep -q "^ $id (the captain's " \
    || { echo "refused: row $id's first cell is not the captain's ruling" >&2; finish "$code"; }
  printf '%s\n' "$first" | grep -qF -- "ruling, no run" \
    || { echo "refused: row $id's first cell does not say ruling, no run" >&2; finish "$code"; }
  for opt in "$@" "--captain-ruling $id"; do
    printf '%s\n' "$first" | grep -o '`[^`]*`' | grep -qxF -- "\`$opt\`" \
      || { echo "refused: row $id's first cell does not name \`$opt\` as a whole backticked string" >&2
           finish "$code"; }
  done
}
ruling_row "$ROW" 95 "--feng-fraction-mult 2.0 --direction-momentum all" \
  "--feng-fraction-mult 3.0 --direction-momentum all" \
  "--feng-fraction-mult 1.5 --direction-momentum all --g-mom-damping 1.0" \
  "--feng-fraction-mult 1.5 --direction-momentum all --g-mom-damping 1.5" \
  "--feng-fraction-mult 1.5 --direction-momentum all --held-momentum 0.5" \
  "--feng-fraction-mult 1.5 --direction-momentum all --held-momentum 1.0" \
  "--feng-fraction-mult 1.5 --direction-momentum all --direction-momentum-damping 1.0"
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
O=runs/t1r51_${head:0:7}_$RUN
mkdir -p runs
mkdir "$O" 2> /dev/null || { echo "refused: the output folder $O exists" >&2; finish 93; }
echo "$head" > "$S/head"
step_ok

STEP=tests
begin
$PY -m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp="$S/tmp/pytest" \
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
build() {  # folder name args: the build, then checked; exit 11 build
  local d="$O/$1" name=$2
  shift 2
  mkdir -p "$d"
  $PY experiments/builds/build_test.py t1 "$@" --out "$d/$name.h5ad" > /dev/null 2>&1 || return 11
  checked "$d" "$name"
}
build_msb() {  # folder name args: the build on its own path with its one setting at MSB's value, then checked; exit 12
  local d="$O/$1" name=$2
  shift 2
  mkdir -p "$d"
  $PY experiments/builds/t1_r51_check.py msb t1 "$@" --out "$d/$name.h5ad" > /dev/null 2>&1 || return 12
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

STEP=Cm30
begin
build Cm30 T1_s2rcCm30 $CM30_ARGS || finish $?
same Cm30 T1_s2rcCm30 "$CM30_SHA" 35
step_ok

STEP=B25d15
begin
build B25d15 T1_s2rcB25d15 $B25D15_ARGS || finish $?
same B25d15 T1_s2rcB25d15 "$B25D15_SHA" 36
step_ok

STEP=B25d20
begin
build B25d20 T1_s2rcB25d20 $B25D20_ARGS || finish $?
same B25d20 T1_s2rcB25d20 "$B25D20_SHA" 37
step_ok

STEP=pooled
begin
build pooled T1_s2rc_pooled $POOLED_ARGS || finish $?
same pooled T1_s2rc_pooled "$POOLED_SHA" 38
step_ok

STEP=MS20msb
begin
build_msb MS20msb T1_s2rcMS20_msb $MS20_ARGS || finish $?
same MS20msb T1_s2rcMS20_msb "$MSB_SHA" 51
step_ok

STEP=MS30msb
begin
build_msb MS30msb T1_s2rcMS30_msb $MS30_ARGS || finish $?
same MS30msb T1_s2rcMS30_msb "$MSB_SHA" 52
step_ok

STEP=d10msb
begin
build_msb d10msb T1_s2rcMSBd10_msb $D10_ARGS || finish $?
same d10msb T1_s2rcMSBd10_msb "$MSB_SHA" 53
step_ok

STEP=d15msb
begin
build_msb d15msb T1_s2rcMSBd15_msb $D15_ARGS || finish $?
same d15msb T1_s2rcMSBd15_msb "$MSB_SHA" 54
step_ok

STEP=H05msb
begin
build_msb H05msb T1_s2rcMSBH05_msb $H05_ARGS || finish $?
same H05msb T1_s2rcMSBH05_msb "$MSB_SHA" 55
step_ok

STEP=H10msb
begin
build_msb H10msb T1_s2rcMSBH10_msb $H10_ARGS || finish $?
same H10msb T1_s2rcMSBH10_msb "$MSB_SHA" 56
step_ok

STEP=m10msb
begin
build_msb m10msb T1_s2rcMSBm10_msb $M10_ARGS || finish $?
same m10msb T1_s2rcMSBm10_msb "$MSB_SHA" 57
step_ok

for STEP in MS20a MS20b; do begin; build "$STEP" T1_s2rcMS20 $MS20_ARGS || finish $?; step_ok; done
for STEP in MS30a MS30b; do begin; build "$STEP" T1_s2rcMS30 $MS30_ARGS || finish $?; step_ok; done
for STEP in d10a d10b; do begin; build "$STEP" T1_s2rcMSBd10 $D10_ARGS || finish $?; step_ok; done
for STEP in d15a d15b; do begin; build "$STEP" T1_s2rcMSBd15 $D15_ARGS || finish $?; step_ok; done
for STEP in H05a H05b; do begin; build "$STEP" T1_s2rcMSBH05 $H05_ARGS || finish $?; step_ok; done
for STEP in H10a H10b; do begin; build "$STEP" T1_s2rcMSBH10 $H10_ARGS || finish $?; step_ok; done
for STEP in m10a m10b; do begin; build "$STEP" T1_s2rcMSBm10 $M10_ARGS || finish $?; step_ok; done

STEP=rows
begin
$PY experiments/builds/t1_r51_check.py rows "$O/MSB/T1_s2rcMSB.h5ad" "$O/H05a/T1_s2rcMSBH05.h5ad" \
  "$O/H10a/T1_s2rcMSBH10.h5ad" > "$S/rows.out" 2> /dev/null || finish 48
step_ok

STEP=done
begin
(cd "$O" && sha256sum D/*.h5ad A/*.h5ad B/*.h5ad C/*.h5ad X2/*.h5ad X25/*.h5ad A15/*.h5ad A20/*.h5ad Cm20/*.h5ad \
  MS/*.h5ad MSB/*.h5ad X30/*.h5ad Cm25/*.h5ad B2d15/*.h5ad B2d20/*.h5ad Cm30/*.h5ad B25d15/*.h5ad B25d20/*.h5ad \
  pooled/*.h5ad MS20msb/*.h5ad MS30msb/*.h5ad d10msb/*.h5ad d15msb/*.h5ad H05msb/*.h5ad H10msb/*.h5ad m10msb/*.h5ad \
  MS20a/*.h5ad MS20b/*.h5ad MS30a/*.h5ad MS30b/*.h5ad d10a/*.h5ad d10b/*.h5ad d15a/*.h5ad d15b/*.h5ad H05a/*.h5ad \
  H05b/*.h5ad H10a/*.h5ad H10b/*.h5ad m10a/*.h5ad m10b/*.h5ad) > "$S/SHA256SUMS" || finish 40
cp "$S/SHA256SUMS" "$O/SHA256SUMS" || finish 40
cmp -s "$O/MS20a/T1_s2rcMS20.h5ad" "$O/MS20b/T1_s2rcMS20.h5ad" || finish 41
cmp -s "$O/MS30a/T1_s2rcMS30.h5ad" "$O/MS30b/T1_s2rcMS30.h5ad" || finish 42
cmp -s "$O/d10a/T1_s2rcMSBd10.h5ad" "$O/d10b/T1_s2rcMSBd10.h5ad" || finish 43
cmp -s "$O/d15a/T1_s2rcMSBd15.h5ad" "$O/d15b/T1_s2rcMSBd15.h5ad" || finish 44
cmp -s "$O/H05a/T1_s2rcMSBH05.h5ad" "$O/H05b/T1_s2rcMSBH05.h5ad" || finish 45
cmp -s "$O/H10a/T1_s2rcMSBH10.h5ad" "$O/H10b/T1_s2rcMSBH10.h5ad" || finish 46
cmp -s "$O/m10a/T1_s2rcMSBm10.h5ad" "$O/m10b/T1_s2rcMSBm10.h5ad" || finish 47
step_ok
