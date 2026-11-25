#!/usr/bin/env bash
# Rows R54's and R58's files (the step-change round's T1 lines), each one change to T1_s2rcMSB's command (t1 --g-mom
# --base lib3 --feng-groups recomputed --feng-fraction-mult 1.5 --direction-momentum all, row R26's):
#   R54  T1_s2rcMSB_LM16  --late-mix e15_17 --captain-ruling R54 (library 3's own chamber mix moved to row N64's
#        shares, experiments/t1_scope/results/n64_late_mix/shares.json); its gate is row T72
#   R58  T1_s2rcMSBm20    --direction-momentum-damping 2.0 --captain-ruling R58 (MS's own damping 2.0); its gate is
#        row T73
# nothing else changed (default seed). Test-board builds on a DGX campaign worktree that holds the relabelling's output
# (data/derived/feng_relabel_c) and its own venv, made there by `uv sync --frozen` (never the reference checkout). Run
# from the root of that worktree, detached, and wait on the exit files with a background poller; never a foreground
# sleep:
#   bash experiments/builds/t1_lm_m20_line.sh <COMMIT> <RUN> <ROWS...>
# ROWS: R54, R58 or both, each once: the files to build. A file is built only if its gate row's outcome cell (the
# row's last of eleven cells, read by code) holds "Outcome: <gate> PASS" as a whole word and neither "Outcome: <gate>
# STOP" nor "Outcome: <gate> VOID"; any requested row whose gate does not say PASS refuses the whole call (87). RUN names this
# run's own scratch folder ~/v/scratch/t1lmm20/<RUN> and its output folder runs/t1lmm20_<commit7>_<RUN>. Neither may
# exist: a run folder is never reused and nothing is ever deleted here, so a second run takes a new name.
# In order, each step first writing the word running to ~/v/scratch/t1lmm20/<RUN>/<STEP>.exit, then its exit code:
#   gate   HEAD is COMMIT (the full hash), the tracked tree is clean, docs/S2_RULES.sha256 checks; for each row asked
#          for: TRIES.md holds exactly one such row whose first cell starts "<ROW> (the captain's ", says "ruling, no
#          run" and names the file's option string and `--captain-ruling <ROW>`, each as a whole backticked string (the
#          builder's own check, _captain_order_names); its gate row's outcome cell says PASS as above; for R54 row N64
#          names exactly one `shares.json sha256 <hex>`, row T72's first cell names the same string, and the shares
#          file is tracked with that sha256; the relabelling's output folder data/derived/feng_relabel_c exists, and
#          .venv's vec.pth names this checkout's src; then it writes the head file
#   tests  the builders', checkers' and launchers' targeted tests (synthetic stages only), 8 BLAS threads
#   D B MS MSB  four identity rebuilds, each once: R13's file, candidate B's, candidate MS's (row R23) and T1_s2rcMSB
#          (row R26), the base of both new files; any of them giving other bytes stops everything:
#   D      the default command (R13's file); R13_SHA
#   B      the --feng-fraction-mult 1.5 command (candidate B); B_SHA
#   MS     the --direction-momentum all command with --captain-ruling R23 (candidate MS); MS_SHA
#   MSB    the --feng-fraction-mult 1.5 --direction-momentum all command with --captain-ruling R26 (T1_s2rcMSB); MSB_SHA
#   LM16own  (R54 asked) R54's command once through experiments/builds/t1_lm_m20_check.py lm16own, its moved shares
#          put back at the pool's own mix; it must give MSB_SHA, T1_s2rcMSB's bytes, or everything stops
#   m20msb (R58 asked) R58's command once through t1_lm_m20_check.py m20msb, MS's damping put back at B0's; it must
#          give MSB_SHA, or everything stops
#   LM16a LM16b  (R54 asked) T1_s2rcMSB_LM16 by its command, twice, into separate folders
#   m20a m20b    (R58 asked) T1_s2rcMSBm20 by its command, twice, into separate folders
#   done   one sha256 list over every built file, each pair compared byte for byte
# A build step also fails (15) unless its stamp says this commit and clean True.
# Builds run at 16 BLAS threads, as the reference builds (T1 files move with the thread count), seed and board as the
# reference commands (the builder's defaults: seed 1000001, board T1:test). Every builder, checker build and `vec check`
# output goes to /dev/null: a test build is read by its exit code, its sha256 and its stamp's commit and clean fields
# (AGENTS.md), never its log, size, record or cell count. Nothing here uploads anything.
# Exit codes: 90 odd RUN, 89 a row ID not R54 or R58 (or one given twice), 91 run folder exists, 96 HEAD is not
# COMMIT, 98 tracked tree not clean, 97 S2_RULES.sha256, 95 a ruling row not as above, 87 a gate row's outcome cell
# does not say PASS, 86 row N64's or T72's shares sha256 or the shares file not as above, 92 no relabelling output,
# 88 .venv not this checkout's, 93 output folder exists, 10 tests, 11 build, 12 checker build, 13 vec check, 14 no
# stamp, 15 stamp not this commit and clean, 20 D is not R13's bytes, 21 B is not candidate B's bytes, 22 MS is not
# candidate MS's bytes, 23 MSB is not T1_s2rcMSB's bytes, 51 LM16own is not MSB's bytes, 52 m20msb is not MSB's bytes,
# 40 hash list, 41 LM16a and LM16b differ, 42 m20a and m20b differ.
set -u
COMMIT=${1:?usage: t1_lm_m20_line.sh COMMIT RUN ROWS}
RUN=${2:?usage: t1_lm_m20_line.sh COMMIT RUN ROWS}
: "${3:?usage: t1_lm_m20_line.sh COMMIT RUN ROWS}"
shift 2
ROWS="$*"
case "$RUN" in
  *[!A-Za-z0-9_-]*) echo "refused: RUN takes letters, digits, _ and - only" >&2; exit 90 ;;
esac
DO54=0
DO58=0
for r in "$@"; do
  case "$r" in
    R54) [ "$DO54" = 0 ] || { echo "refused: row R54 is given twice" >&2; exit 89; }; DO54=1 ;;
    R58) [ "$DO58" = 0 ] || { echo "refused: row R58 is given twice" >&2; exit 89; }; DO58=1 ;;
    *) echo "refused: the rows of these files are R54 and R58, and only those" >&2; exit 89 ;;
  esac
done
export PATH=$HOME/.local/bin:$PATH
unset UV_PROJECT_ENVIRONMENT  # uv runs this checkout's own .venv, the one the gate checks
S=$HOME/v/scratch/t1lmm20/$RUN
mkdir -p "$HOME/v/scratch/t1lmm20"
mkdir "$S" 2> /dev/null || { echo "refused: the run folder $S exists (a run folder is never reused)" >&2; exit 91; }
mkdir "$S/tmp"
export NUMPY_MADVISE_HUGEPAGE=0 TMPDIR="$S/tmp"
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
PY="uv run --frozen --no-sync python -I -B"
BOARD=T1:test
TRIES=experiments/harness/TRIES.md
SHARES=experiments/t1_scope/results/n64_late_mix/shares.json
R13_SHA=ad8f24c71c90ca17bce5f8876d503471ab3faa5b17a0878addb00aff689fca4f  # set t1rc, row R13 (the default command)
B_SHA=c7240cb6a7f7ee93d037f437832efbaad004f1fbe608d80072f8efd5d20a0062  # candidate B, set t1cand_1420c70
MS_SHA=6103d10535665634ee0cad3390590a5df132d1805682fb1d309bb19581fdaa5e  # candidate MS, t1ms_6bc36d9_r1.SHA256SUMS
MSB_SHA=78bb17eea5721e590cfd2c9f707e9cd6f2f496e081acede4f697dc4035871da2  # T1_s2rcMSB, t1msc_a3e879a_m1.SHA256SUMS
BASE_ARGS="--g-mom --base lib3 --feng-groups recomputed"
B_ARGS="--g-mom --base lib3 --feng-groups recomputed --feng-fraction-mult 1.5"
MS_ARGS="--g-mom --base lib3 --feng-groups recomputed --direction-momentum all --captain-ruling R23"
MSB_ARGS="--g-mom --base lib3 --feng-groups recomputed --feng-fraction-mult 1.5 --direction-momentum all --captain-ruling R26"
LM16_ARGS="--g-mom --base lib3 --feng-groups recomputed --feng-fraction-mult 1.5 --direction-momentum all --late-mix e15_17 --captain-ruling R54"
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
ruling_row() {  # ID OPTION: one row ID, its first cell the captain's ruling naming OPTION whole, else 95
  local id=$1 first opt
  shift
  [ "$(grep -c "^| $id " "$TRIES")" = 1 ] || { echo "refused: TRIES.md does not hold exactly one row $id" >&2; finish 95; }
  first=$(grep "^| $id " "$TRIES" | cut -d'|' -f2)
  printf '%s\n' "$first" | grep -q "^ $id (the captain's " \
    || { echo "refused: row $id's first cell is not the captain's ruling" >&2; finish 95; }
  printf '%s\n' "$first" | grep -qF -- "ruling, no run" \
    || { echo "refused: row $id's first cell does not say ruling, no run" >&2; finish 95; }
  for opt in "$@" "--captain-ruling $id"; do
    printf '%s\n' "$first" | grep -o '`[^`]*`' | grep -qxF -- "\`$opt\`" \
      || { echo "refused: row $id's first cell does not name \`$opt\` as a whole backticked string" >&2; finish 95; }
  done
}
gate_passed() {  # GATE: exactly one row GATE of eleven cells whose last cell says PASS and neither STOP nor VOID, else 87
  local g=$1 cell
  [ "$(grep -c "^| $g " "$TRIES")" = 1 ] || { echo "refused: TRIES.md does not hold exactly one row $g" >&2; finish 87; }
  cell=$(grep "^| $g " "$TRIES" | awk -F'|' 'NF == 13 { print $12 }')
  [ -n "$cell" ] || { echo "refused: row $g does not have eleven cells" >&2; finish 87; }
  printf '%s\n' "$cell" | grep -qE -- "Outcome: $g PASS([^[:alnum:]_]|\$)" \
    || { echo "refused: row $g's outcome cell does not say Outcome: $g PASS" >&2; finish 87; }
  if printf '%s\n' "$cell" | grep -qE -- "Outcome: $g (STOP|VOID)([^[:alnum:]_]|\$)"; then
    echo "refused: row $g's outcome cell also says STOP or VOID" >&2
    finish 87
  fi
}
if [ "$DO54" = 1 ]; then
  ruling_row R54 "--late-mix e15_17"
  gate_passed T72
  [ "$(grep -c '^| N64 ' "$TRIES")" = 1 ] || { echo "refused: TRIES.md does not hold exactly one row N64" >&2; finish 86; }
  named=$(grep '^| N64 ' "$TRIES" | grep -oE '`shares\.json sha256 [0-9a-f]{64}`' | sort -u)
  [ "$(printf '%s\n' "$named" | grep -c .)" = 1 ] \
    || { echo "refused: row N64 does not name exactly one \`shares.json sha256 <hex>\`" >&2; finish 86; }
  grep '^| T72 ' "$TRIES" | cut -d'|' -f2 | grep -qF -- "$named" \
    || { echo "refused: row T72's first cell does not name row N64's shares sha256" >&2; finish 86; }
  git ls-files --error-unmatch "$SHARES" > /dev/null 2>&1 || { echo "refused: $SHARES is not tracked" >&2; finish 86; }
  got=$(sha256sum "$SHARES" | cut -d' ' -f1)
  [ "$named" = "\`shares.json sha256 $got\`" ] || { echo "refused: $SHARES is not the file row N64 names" >&2; finish 86; }
fi
if [ "$DO58" = 1 ]; then
  ruling_row R58 "--feng-fraction-mult 1.5 --direction-momentum all --direction-momentum-damping 2.0"
  gate_passed T73
fi
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
O=runs/t1lmm20_${head:0:7}_$RUN
mkdir -p runs
mkdir "$O" 2> /dev/null || { echo "refused: the output folder $O exists" >&2; finish 93; }
echo "$head" > "$S/head"
echo "$ROWS" > "$S/rows"
step_ok

STEP=tests
begin
$PY -m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp="$S/tmp/pytest" \
  tests/test_t1_lm_m20_build.py tests/test_t1_lm_m20_launchers.py tests/test_t1_r51_build.py \
  tests/test_t1_msb_build.py tests/test_t1_ms_build.py tests/test_t1_candidates_build.py tests/test_hidden_counts.py \
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
build() {  # folder name args: the build, then checked; exit 11 build
  local d="$O/$1" name=$2
  shift 2
  mkdir -p "$d"
  $PY experiments/builds/build_test.py t1 "$@" --out "$d/$name.h5ad" > /dev/null 2>&1 || return 11
  checked "$d" "$name"
}
build_check() {  # mode folder name args: the build through the checker's mode, then checked; exit 12
  local mode=$1 d="$O/$2" name=$3
  shift 3
  mkdir -p "$d"
  $PY experiments/builds/t1_lm_m20_check.py "$mode" t1 "$@" --out "$d/$name.h5ad" > /dev/null 2>&1 || return 12
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

STEP=B
begin
build B T1_s2rcB $B_ARGS || finish $?
same B T1_s2rcB "$B_SHA" 21
step_ok

STEP=MS
begin
build MS T1_s2rc_MS $MS_ARGS || finish $?
same MS T1_s2rc_MS "$MS_SHA" 22
step_ok

STEP=MSB
begin
build MSB T1_s2rcMSB $MSB_ARGS || finish $?
same MSB T1_s2rcMSB "$MSB_SHA" 23
step_ok

FOLDERS="D/*.h5ad B/*.h5ad MS/*.h5ad MSB/*.h5ad"
if [ "$DO54" = 1 ]; then
  STEP=LM16own
  begin
  build_check lm16own LM16own T1_s2rcMSB_LM16_own $LM16_ARGS || finish $?
  same LM16own T1_s2rcMSB_LM16_own "$MSB_SHA" 51
  step_ok
  FOLDERS="$FOLDERS LM16own/*.h5ad"
fi
if [ "$DO58" = 1 ]; then
  STEP=m20msb
  begin
  build_check m20msb m20msb T1_s2rcMSBm20_msb $M20_ARGS || finish $?
  same m20msb T1_s2rcMSBm20_msb "$MSB_SHA" 52
  step_ok
  FOLDERS="$FOLDERS m20msb/*.h5ad"
fi
if [ "$DO54" = 1 ]; then
  for STEP in LM16a LM16b; do begin; build "$STEP" T1_s2rcMSB_LM16 $LM16_ARGS || finish $?; step_ok; done
  FOLDERS="$FOLDERS LM16a/*.h5ad LM16b/*.h5ad"
fi
if [ "$DO58" = 1 ]; then
  for STEP in m20a m20b; do begin; build "$STEP" T1_s2rcMSBm20 $M20_ARGS || finish $?; step_ok; done
  FOLDERS="$FOLDERS m20a/*.h5ad m20b/*.h5ad"
fi

STEP=done
begin
# shellcheck disable=SC2086  # FOLDERS is a list of the built folders' globs, expanded in the output folder
(cd "$O" && sha256sum $FOLDERS) > "$S/SHA256SUMS" || finish 40
cp "$S/SHA256SUMS" "$O/SHA256SUMS" || finish 40
if [ "$DO54" = 1 ]; then
  cmp -s "$O/LM16a/T1_s2rcMSB_LM16.h5ad" "$O/LM16b/T1_s2rcMSB_LM16.h5ad" || finish 41
fi
if [ "$DO58" = 1 ]; then
  cmp -s "$O/m20a/T1_s2rcMSBm20.h5ad" "$O/m20b/T1_s2rcMSBm20.h5ad" || finish 42
fi
step_ok
