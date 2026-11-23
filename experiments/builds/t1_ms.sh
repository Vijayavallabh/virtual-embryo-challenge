#!/usr/bin/env bash
# Candidate MS of row R23 (experiments/harness/TRIES.md; the captain's order of 2026-11-23, docs/PLAN.md decision log
# 2026-11-23 (c)): R13's file (the default command, set t1rc) with B0's released per-category damped mean momentum added
# on the cells whose label's group has a recomputed direction, before the atlas step and g_mom, nothing else changed
# (default seed). Test-board builds on the DGX campaign checkout that holds the relabelling's output
# (~/v/vec-campaign-d, never the reference checkout). Run from the root of the checkout or detached worktree it is
# called in, and wait on the exit files with a background poller; never a foreground sleep:
#   bash experiments/builds/t1_ms.sh <COMMIT> <RUN>
# RUN names this run's own scratch folder ~/v/scratch/t1ms/<RUN> and its output folder runs/t1ms_<commit7>_<RUN>.
# Neither may exist: a run folder is never reused and nothing is ever deleted here, so a second run takes a new name.
# In order, each step first writing the word running to ~/v/scratch/t1ms/<RUN>/<STEP>.exit, then its exit code:
#   gate   HEAD is COMMIT (the full hash), the tracked tree is clean, docs/S2_RULES.sha256 checks, TRIES.md holds one
#          row R23 whose first cell names `--direction-momentum all` and `--captain-ruling R23`, and the relabelling's
#          output folder data/derived/feng_relabel_c exists; then it writes the head file
#   tests  the builder's and this launcher's targeted tests (synthetic stages only), 8 BLAS threads
#   D      the default command (R13's file, set t1rc) once; it must give R13_SHA or everything stops
#   MSa MSb  the default command with --direction-momentum all --captain-ruling R23, twice, into separate folders
#   done   one sha256 list over every built file, the pair compared byte for byte
# A build step also fails (15) unless its stamp says this commit and clean True.
# Builds run at 16 BLAS threads, as the reference builds (T1 files move with the thread count), seed and board as the
# reference commands (the builder's defaults: seed 1000001, board T1:test). Every builder and `vec check` output goes
# to /dev/null: a test build is read by its exit code, its sha256 and its stamp's commit and clean fields only
# (AGENTS.md), never its log, size, cell count or sidecar record. Nothing here uploads anything.
# Exit codes: 90 odd RUN, 91 run folder exists, 96 HEAD is not COMMIT, 98 tracked tree not clean, 97 S2_RULES.sha256,
# 95 not exactly one row R23 in TRIES.md or its first cell does not name the option and the ruling, 92 no relabelling
# output, 93 output folder exists, 10 tests, 11 build, 13 vec check, 14 no stamp, 15 stamp not this commit and clean,
# 20 D is not R13's bytes, 30 hash list, 31 MSa and MSb differ.
set -u
COMMIT=${1:?usage: t1_ms.sh COMMIT RUN}
RUN=${2:?usage: t1_ms.sh COMMIT RUN}
case "$RUN" in
  *[!A-Za-z0-9_-]*) echo "refused: RUN takes letters, digits, _ and - only" >&2; exit 90 ;;
esac
export PATH=$HOME/.local/bin:$PATH
S=$HOME/v/scratch/t1ms/$RUN
mkdir -p "$HOME/v/scratch/t1ms"
mkdir "$S" 2> /dev/null || { echo "refused: the run folder $S exists (a run folder is never reused)" >&2; exit 91; }
mkdir "$S/tmp"
export NUMPY_MADVISE_HUGEPAGE=0 TMPDIR="$S/tmp"
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
PY="uv run --frozen --no-sync python -I -B"
BOARD=T1:test
R13_SHA=ad8f24c71c90ca17bce5f8876d503471ab3faa5b17a0878addb00aff689fca4f  # set t1rc, row R13 (the default command)
BASE_ARGS="--g-mom --base lib3 --feng-groups recomputed"
MS_ARGS="--g-mom --base lib3 --feng-groups recomputed --direction-momentum all --captain-ruling R23"
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
[ "$(grep -c '^| R23 ' experiments/harness/TRIES.md)" = 1 ] \
  || { echo "refused: TRIES.md does not hold exactly one row R23" >&2; finish 95; }
first=$(grep '^| R23 ' experiments/harness/TRIES.md | cut -d'|' -f2)
for opt in '`--direction-momentum all`' '`--captain-ruling R23`'; do
  printf '%s\n' "$first" | grep -qF -- "$opt" \
    || { echo "refused: row R23's first cell does not name $opt" >&2; finish 95; }
done
[ -d data/derived/feng_relabel_c ] \
  || { echo "refused: the relabelling's output data/derived/feng_relabel_c is not in this checkout" >&2; finish 92; }
O=runs/t1ms_${head:0:7}_$RUN
mkdir -p runs
mkdir "$O" 2> /dev/null || { echo "refused: the output folder $O exists" >&2; finish 93; }
echo "$head" > "$S/head"
step_ok

STEP=tests
begin
$PY -m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp="$S/tmp/pytest" \
  tests/test_t1_ms_launcher.py tests/test_t1_ms_build.py tests/test_t65_ms_read.py tests/test_t1_line3_build.py \
  tests/test_t1_damping_line_build.py tests/test_t1_fraction_line_build.py tests/test_t1_fraction2_build.py \
  tests/test_t1_candidates_build.py tests/test_t1_bets.py tests/test_hidden_counts.py tests/test_s2_rules.py \
  > "$S/tests.out" 2>&1 || finish 10
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

for STEP in MSa MSb; do begin; build "$STEP" T1_s2rc_MS $MS_ARGS || finish $?; step_ok; done

STEP=done
begin
(cd "$O" && sha256sum D/*.h5ad MSa/*.h5ad MSb/*.h5ad) > "$S/SHA256SUMS" || finish 30
cp "$S/SHA256SUMS" "$O/SHA256SUMS" || finish 30
cmp -s "$O/MSa/T1_s2rc_MS.h5ad" "$O/MSb/T1_s2rc_MS.h5ad" || finish 31
step_ok
