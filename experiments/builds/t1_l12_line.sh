#!/usr/bin/env bash
# Row R60's seven files (T1's line 12, the captain's order of 2026-11-25), each one change to T1_s2rcMS20's command
# (t1 --g-mom --base lib3 --feng-groups recomputed --feng-fraction-mult 2.0 --direction-momentum all, row R51's file),
# each only with --captain-ruling R60, and row R46's pair (the captain's order of 2026-11-24, T1's lines 8 and 9) with
# --captain-ruling R46:
#   T1_s2rcMS20d10   ... --g-mom-damping 1.0 (g_mom's damping 1.0)
#   T1_s2rcMS20H05   ... --held-momentum 0.5 (Blood's own released change at damping 0.5 after MS's step and before
#                    g_mom: R1 reversed for Blood in this file and H10 only, as in row R51's MSBH05 and MSBH10)
#   T1_s2rcMS20d20   ... --g-mom-damping 2.0
#   T1_s2rcMS20m10   ... --direction-momentum-damping 1.0 (MS's own damping)
#   T1_s2rcMS20d15   ... --g-mom-damping 1.5
#   T1_s2rcMS20H10   ... --held-momentum 1.0
#   T1_s2rcMS20m20   ... --direction-momentum-damping 2.0
#   T1_s2rcB25MSH05  --feng-fraction-mult 2.5 --direction-momentum all --held-momentum 0.5 (row R44's Blood step on
#                    T1_s2rcB25MS's command, row R45's)
# nothing else changed (default seed). All eight are built here, the second wave included: the uploads are gated, not the
# build. Test-board builds on a DGX campaign worktree that holds the relabelling's output (data/derived/feng_relabel_c)
# and its own venv, made there by `uv sync --frozen` (never the reference checkout). Run from the root of that worktree,
# detached, and wait on the exit files with a background poller; never a foreground sleep:
#   bash experiments/builds/t1_l12_line.sh <COMMIT> <RUN>
# RUN names this run's own scratch folder ~/v/scratch/t1l12/<RUN> and its output folder runs/t1l12_<commit7>_<RUN>.
# Neither may exist: a run folder is never reused and nothing is ever deleted here, so a second run takes a new name.
# In order, each step first writing the word running to ~/v/scratch/t1l12/<RUN>/<STEP>.exit, then its exit code:
#   gate   HEAD is COMMIT (the full hash), the tracked tree is clean, docs/S2_RULES.sha256 checks, TRIES.md holds exactly
#          one row R60 whose first cell starts "R60 (the captain's ", says "ruling, no run" and names each of the seven
#          option strings and `--captain-ruling R60`, each as a whole backticked string, and exactly one row R46 whose
#          first cell starts "R46 (the captain's ", says "ruling, no run" and names T1_s2rcB25MSH05 as a whole word
#          (the row names no option string), the relabelling's output folder data/derived/feng_relabel_c exists, and
#          .venv's vec.pth names this checkout's src; then it writes the head file
#   tests  the builders', checkers' and launchers' targeted tests (synthetic stages only), 8 BLAS threads
#   D ... pooled MS20 B25MS  the twenty-one identity rebuilds, each once: t1_r51_line.sh's nineteen, then T1_s2rcMS20
#          (row R51's file, the base of row R60's; t1r51_0cdcc6f_r1.SHA256SUMS) and T1_s2rcB25MS (row R45's file, the
#          base of row R46's; t1hm_20d5129_h1.SHA256SUMS); any of them giving other bytes stops everything
#   d10base H05base d20base m10base d15base H10base m20base pairbase  each new command once through
#          experiments/builds/t1_l12_check.py base, on its own path with its one setting at its base file's value
#          (g_mom's damping at 0.5, Blood's factors at 1, MS's damping at the default); each of row R60's must give
#          MS20_SHA, T1_s2rcMS20's bytes, and row R46's B25MS_SHA, T1_s2rcB25MS's, or everything stops
#   d10a d10b ... pairb  each new file by its command, twice, into separate folders: the first wave (d10, H05, d20,
#          m10), the second (d15, H10, m20), then the pair
#   done   one sha256 list over every built file, each pair compared byte for byte
# A build step also fails (15) unless its stamp says this commit and clean True.
# Builds run at 16 BLAS threads, as the reference builds (T1 files move with the thread count), seed and board as the
# reference commands (the builder's defaults: seed 1000001, board T1:test). Every builder, checker build and `vec check`
# output goes to /dev/null: a test build is read by its exit code, its sha256 and its stamp's commit and clean fields
# (AGENTS.md), never its log, size, record or cell count. Nothing here uploads anything.
# Exit codes: 90 odd RUN, 91 run folder exists, 96 HEAD is not COMMIT, 98 tracked tree not clean, 97 S2_RULES.sha256,
# 95 row R60 not as above, 94 row R46 not as above, 92 no relabelling output, 88 .venv not this checkout's, 93 output
# folder exists, 10 tests, 11 build, 12 checker build, 13 vec check, 14 no stamp, 15 stamp not this commit and clean,
# as t1_r51_line.sh: 20 D is not R13's bytes, 21 A, 22 B, 23 C, 24 X2, 25 X25, 26 A15, 27 A20, 28 Cm20, 29 MS,
# 30 MSB, 31 X30, 32 Cm25, 33 B2d15, 34 B2d20, 35 Cm30, 36 B25d15, 37 B25d20, 38 pooled is not its known bytes;
# 39 MS20 is not T1_s2rcMS20's bytes, 50 B25MS is not T1_s2rcB25MS's bytes,
# 51 d10base is not MS20's bytes, 52 H05base is not, 53 d20base is not, 54 m10base is not, 55 d15base is not,
# 56 H10base is not, 57 m20base is not, 58 pairbase is not B25MS's bytes, 40 hash list, 41 d10a and d10b differ, 42 H05a
# and H05b differ, 43 d20a and d20b differ, 44 m10a and m10b differ, 45 d15a and d15b differ, 46 H10a and H10b differ,
# 47 m20a and m20b differ, 48 paira and pairb differ.
set -u
COMMIT=${1:?usage: t1_l12_line.sh COMMIT RUN}
RUN=${2:?usage: t1_l12_line.sh COMMIT RUN}
case "$RUN" in
  *[!A-Za-z0-9_-]*) echo "refused: RUN takes letters, digits, _ and - only" >&2; exit 90 ;;
esac
export PATH=$HOME/.local/bin:$PATH
unset UV_PROJECT_ENVIRONMENT  # uv runs this checkout's own .venv, the one the gate checks
S=$HOME/v/scratch/t1l12/$RUN
mkdir -p "$HOME/v/scratch/t1l12"
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
MS20_SHA=629d3a1cffa8151f0345a169ba97420e4ee3b43e379bc3bf7d33ab947dcb62ef  # T1_s2rcMS20, t1r51_0cdcc6f_r1.SHA256SUMS
B25MS_SHA=3127a2d9c92aa064e3f0fdf400447a074b7964ceacf2274da8a14ada7425f411  # T1_s2rcB25MS, t1hm_20d5129_h1.SHA256SUMS
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
B25MS_ARGS="--g-mom --base lib3 --feng-groups recomputed --feng-fraction-mult 2.5 --direction-momentum all --captain-ruling R45"
L12="--g-mom --base lib3 --feng-groups recomputed --feng-fraction-mult 2.0 --direction-momentum all"
D10_ARGS="$L12 --g-mom-damping 1.0 --captain-ruling R60"
H05_ARGS="$L12 --held-momentum 0.5 --captain-ruling R60"
D20_ARGS="$L12 --g-mom-damping 2.0 --captain-ruling R60"
M10_ARGS="$L12 --direction-momentum-damping 1.0 --captain-ruling R60"
D15_ARGS="$L12 --g-mom-damping 1.5 --captain-ruling R60"
H10_ARGS="$L12 --held-momentum 1.0 --captain-ruling R60"
M20_ARGS="$L12 --direction-momentum-damping 2.0 --captain-ruling R60"
PAIR_ARGS="--g-mom --base lib3 --feng-groups recomputed --feng-fraction-mult 2.5 --direction-momentum all --held-momentum 0.5 --captain-ruling R46"
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
r60=$(first_cell R60 95) || exit $?
for opt in "--feng-fraction-mult 2.0 --direction-momentum all --g-mom-damping 1.0" \
  "--feng-fraction-mult 2.0 --direction-momentum all --held-momentum 0.5" \
  "--feng-fraction-mult 2.0 --direction-momentum all --g-mom-damping 2.0" \
  "--feng-fraction-mult 2.0 --direction-momentum all --direction-momentum-damping 1.0" \
  "--feng-fraction-mult 2.0 --direction-momentum all --g-mom-damping 1.5" \
  "--feng-fraction-mult 2.0 --direction-momentum all --held-momentum 1.0" \
  "--feng-fraction-mult 2.0 --direction-momentum all --direction-momentum-damping 2.0" "--captain-ruling R60"; do
  printf '%s\n' "$r60" | grep -o '`[^`]*`' | grep -qxF -- "\`$opt\`" \
    || { echo "refused: row R60's first cell does not name \`$opt\` as a whole backticked string" >&2; finish 95; }
done
r46=$(first_cell R46 94) || exit $?
printf '%s\n' "$r46" | grep -qE '(^|[^A-Za-z0-9_])T1_s2rcB25MSH05([^A-Za-z0-9_]|$)' \
  || { echo "refused: row R46's first cell does not name T1_s2rcB25MSH05" >&2; finish 94; }
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
O=runs/t1l12_${head:0:7}_$RUN
mkdir -p runs
mkdir "$O" 2> /dev/null || { echo "refused: the output folder $O exists" >&2; finish 93; }
echo "$head" > "$S/head"
step_ok

STEP=tests
begin
$PY -m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp="$S/tmp/pytest" \
  tests/test_t1_l12_launcher.py tests/test_t1_l12_build.py tests/test_t1_r51_launcher.py tests/test_t1_r51_build.py \
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
build_base() {  # folder name args: the build on its own path with its one setting at its base's value; exit 12
  local d="$O/$1" name=$2
  shift 2
  mkdir -p "$d"
  $PY experiments/builds/t1_l12_check.py base t1 "$@" --out "$d/$name.h5ad" > /dev/null 2>&1 || return 12
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
identity D T1_s2rc "$R13_SHA" 20 $BASE_ARGS
identity A T1_s2rcA "$A_SHA" 21 $A_ARGS
identity B T1_s2rcB "$B_SHA" 22 $B_ARGS
identity C T1_s2rcC "$C_SHA" 23 $C_ARGS
identity X2 T1_s2rcB2 "$X2_SHA" 24 $X2_ARGS
identity X25 T1_s2rcB25 "$X25_SHA" 25 $X25_ARGS
identity A15 T1_s2rcA15 "$A15_SHA" 26 $A15_ARGS
identity A20 T1_s2rcA20 "$A20_SHA" 27 $A20_ARGS
identity Cm20 T1_s2rcCm20 "$CM20_SHA" 28 $CM20_ARGS
identity MS T1_s2rc_MS "$MS_SHA" 29 $MS_ARGS
identity MSB T1_s2rcMSB "$MSB_SHA" 30 $MSB_ARGS
identity X30 T1_s2rcB30 "$X30_SHA" 31 $X30_ARGS
identity Cm25 T1_s2rcCm25 "$CM25_SHA" 32 $CM25_ARGS
identity B2d15 T1_s2rcB2d15 "$B2D15_SHA" 33 $B2D15_ARGS
identity B2d20 T1_s2rcB2d20 "$B2D20_SHA" 34 $B2D20_ARGS
identity Cm30 T1_s2rcCm30 "$CM30_SHA" 35 $CM30_ARGS
identity B25d15 T1_s2rcB25d15 "$B25D15_SHA" 36 $B25D15_ARGS
identity B25d20 T1_s2rcB25d20 "$B25D20_SHA" 37 $B25D20_ARGS
identity pooled T1_s2rc_pooled "$POOLED_SHA" 38 $POOLED_ARGS
identity MS20 T1_s2rcMS20 "$MS20_SHA" 39 $MS20_ARGS
identity B25MS T1_s2rcB25MS "$B25MS_SHA" 50 $B25MS_ARGS

base_check() {  # STEP name sha code args: one new command with its one setting at its base's value
  STEP=$1
  local name=$2 sha=$3 code=$4
  shift 4
  begin
  build_base "$STEP" "$name" "$@" || finish $?
  same "$STEP" "$name" "$sha" "$code"
  step_ok
}
base_check d10base T1_s2rcMS20d10_base "$MS20_SHA" 51 $D10_ARGS
base_check H05base T1_s2rcMS20H05_base "$MS20_SHA" 52 $H05_ARGS
base_check d20base T1_s2rcMS20d20_base "$MS20_SHA" 53 $D20_ARGS
base_check m10base T1_s2rcMS20m10_base "$MS20_SHA" 54 $M10_ARGS
base_check d15base T1_s2rcMS20d15_base "$MS20_SHA" 55 $D15_ARGS
base_check H10base T1_s2rcMS20H10_base "$MS20_SHA" 56 $H10_ARGS
base_check m20base T1_s2rcMS20m20_base "$MS20_SHA" 57 $M20_ARGS
base_check pairbase T1_s2rcB25MSH05_base "$B25MS_SHA" 58 $PAIR_ARGS

for STEP in d10a d10b; do begin; build "$STEP" T1_s2rcMS20d10 $D10_ARGS || finish $?; step_ok; done
for STEP in H05a H05b; do begin; build "$STEP" T1_s2rcMS20H05 $H05_ARGS || finish $?; step_ok; done
for STEP in d20a d20b; do begin; build "$STEP" T1_s2rcMS20d20 $D20_ARGS || finish $?; step_ok; done
for STEP in m10a m10b; do begin; build "$STEP" T1_s2rcMS20m10 $M10_ARGS || finish $?; step_ok; done
for STEP in d15a d15b; do begin; build "$STEP" T1_s2rcMS20d15 $D15_ARGS || finish $?; step_ok; done
for STEP in H10a H10b; do begin; build "$STEP" T1_s2rcMS20H10 $H10_ARGS || finish $?; step_ok; done
for STEP in m20a m20b; do begin; build "$STEP" T1_s2rcMS20m20 $M20_ARGS || finish $?; step_ok; done
for STEP in paira pairb; do begin; build "$STEP" T1_s2rcB25MSH05 $PAIR_ARGS || finish $?; step_ok; done

STEP=done
begin
(cd "$O" && sha256sum D/*.h5ad A/*.h5ad B/*.h5ad C/*.h5ad X2/*.h5ad X25/*.h5ad A15/*.h5ad A20/*.h5ad Cm20/*.h5ad \
  MS/*.h5ad MSB/*.h5ad X30/*.h5ad Cm25/*.h5ad B2d15/*.h5ad B2d20/*.h5ad Cm30/*.h5ad B25d15/*.h5ad B25d20/*.h5ad \
  pooled/*.h5ad MS20/*.h5ad B25MS/*.h5ad d10base/*.h5ad H05base/*.h5ad d20base/*.h5ad m10base/*.h5ad d15base/*.h5ad \
  H10base/*.h5ad m20base/*.h5ad pairbase/*.h5ad d10a/*.h5ad d10b/*.h5ad H05a/*.h5ad H05b/*.h5ad d20a/*.h5ad \
  d20b/*.h5ad m10a/*.h5ad m10b/*.h5ad d15a/*.h5ad d15b/*.h5ad H10a/*.h5ad H10b/*.h5ad m20a/*.h5ad m20b/*.h5ad \
  paira/*.h5ad pairb/*.h5ad) > "$S/SHA256SUMS" || finish 40
cp "$S/SHA256SUMS" "$O/SHA256SUMS" || finish 40
cmp -s "$O/d10a/T1_s2rcMS20d10.h5ad" "$O/d10b/T1_s2rcMS20d10.h5ad" || finish 41
cmp -s "$O/H05a/T1_s2rcMS20H05.h5ad" "$O/H05b/T1_s2rcMS20H05.h5ad" || finish 42
cmp -s "$O/d20a/T1_s2rcMS20d20.h5ad" "$O/d20b/T1_s2rcMS20d20.h5ad" || finish 43
cmp -s "$O/m10a/T1_s2rcMS20m10.h5ad" "$O/m10b/T1_s2rcMS20m10.h5ad" || finish 44
cmp -s "$O/d15a/T1_s2rcMS20d15.h5ad" "$O/d15b/T1_s2rcMS20d15.h5ad" || finish 45
cmp -s "$O/H10a/T1_s2rcMS20H10.h5ad" "$O/H10b/T1_s2rcMS20H10.h5ad" || finish 46
cmp -s "$O/m20a/T1_s2rcMS20m20.h5ad" "$O/m20b/T1_s2rcMS20m20.h5ad" || finish 47
cmp -s "$O/paira/T1_s2rcB25MSH05.h5ad" "$O/pairb/T1_s2rcB25MSH05.h5ad" || finish 48
step_ok
