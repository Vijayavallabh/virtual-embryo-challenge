#!/usr/bin/env bash
# Row R20's T1 file (experiments/harness/TRIES.md; the captain's order of 2026-11-21, docs/PLAN.md decision log
# 2026-11-21 (x); the last amendment of docs/S2_RULES.md): R13's recomputed-groups file with the outside atlas's
# direction step at 2.0 times the calendar fraction. Test-board builds on the DGX campaign checkout that holds the
# relabelling's output (~/v/vec-campaign-d, never the reference checkout). Run from the checkout's root, detached, and
# wait on the exit files with a background poller; never a foreground sleep:
#   bash experiments/builds/t1_fraction2.sh <COMMIT> <RUN>
# RUN names this run's own scratch folder ~/v/scratch/t1x2/<RUN> and its output folder runs/t1x2_<commit7>_<RUN>.
# Neither may exist: a run folder is never reused and nothing is ever deleted here, so a second run takes a new name.
# In order, each step writing ~/v/scratch/t1x2/<RUN>/<STEP>.exit:
#   gate   HEAD is COMMIT (the full hash), the tracked tree is clean, docs/S2_RULES.sha256 checks, row R20 is in
#          TRIES.md, docs/S2_RULES.md names row R20, and the relabelling's output folder data/derived/feng_relabel_c
#          exists
#   tests  the builder's targeted tests (synthetic stages only), 8 BLAS threads
#   D      the default command (R13's file, set t1rc) once; it must give R13_SHA or everything stops
#   B      the 1.5 command (candidate B) once; it must give B_SHA or everything stops
#   X2a X2b  the 2.0 command with --captain-ruling R20, twice, into separate folders
#   done   one sha256 list over every built file, X2a and X2b compared byte for byte
# A build step also fails (15) unless its stamp says this commit and clean True.
# Builds run at 16 BLAS threads, as the reference builds (T1 files move with the thread count), seed and board as the
# reference commands (the builder's defaults: seed 1000001, board T1:test). Every builder and `vec check` output goes
# to /dev/null: a test build is read by its exit code, its sha256 and its stamp's commit and clean fields only
# (AGENTS.md), never its log, size, cell count or sidecar record. Nothing here uploads anything.
# Exit codes: 90 odd RUN, 91 run folder exists, 96 HEAD is not COMMIT, 98 tracked tree not clean, 97 S2_RULES.sha256,
# 95 no row R20 in TRIES.md, 94 S2_RULES.md does not name R20, 92 no relabelling output, 93 output folder exists,
# 10 tests, 11 build, 13 vec check, 14 no stamp, 15 stamp not this commit and clean, 20 D is not R13's bytes, 21 B is
# not candidate B's bytes, 30 hash list, 31 X2a and X2b differ.
set -u
COMMIT=${1:?usage: t1_fraction2.sh COMMIT RUN}
RUN=${2:?usage: t1_fraction2.sh COMMIT RUN}
case "$RUN" in
  *[!A-Za-z0-9_-]*) echo "refused: RUN takes letters, digits, _ and - only" >&2; exit 90 ;;
esac
export PATH=$HOME/.local/bin:$PATH
S=$HOME/v/scratch/t1x2/$RUN
mkdir -p "$HOME/v/scratch/t1x2"
mkdir "$S" 2> /dev/null || { echo "refused: the run folder $S exists (a run folder is never reused)" >&2; exit 91; }
mkdir "$S/tmp"
export NUMPY_MADVISE_HUGEPAGE=0 TMPDIR="$S/tmp"
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
PY="uv run --frozen --no-sync python -I -B"
BOARD=T1:test
R13_SHA=ad8f24c71c90ca17bce5f8876d503471ab3faa5b17a0878addb00aff689fca4f  # set t1rc, row R13 (the default command)
B_SHA=c7240cb6a7f7ee93d037f437832efbaad004f1fbe608d80072f8efd5d20a0062  # candidate B, set t1cand_1420c70
BASE_ARGS="--g-mom --base lib3 --feng-groups recomputed"
B_ARGS="--g-mom --base lib3 --feng-groups recomputed --feng-fraction-mult 1.5"
X2_ARGS="--g-mom --base lib3 --feng-groups recomputed --feng-fraction-mult 2.0 --captain-ruling R20"
STEP=gate
finish() { echo "$1" > "$S/$STEP.exit"; exit "$1"; }
step_ok() { echo 0 > "$S/$STEP.exit"; }

head=$(git rev-parse HEAD)
[ "$head" = "$COMMIT" ] || { echo "refused: HEAD is not $COMMIT (give the full hash)" >&2; finish 96; }
if [ -n "$(git status --porcelain --untracked-files=no)" ]; then
  echo "refused: tracked files differ from the commit (the build stamps must say clean)" >&2
  finish 98
fi
sha256sum -c docs/S2_RULES.sha256 > /dev/null 2>&1 || finish 97
grep -q '^| R20 ' experiments/harness/TRIES.md || { echo "refused: row R20 is not in TRIES.md" >&2; finish 95; }
grep -q 'row R20' docs/S2_RULES.md || { echo "refused: docs/S2_RULES.md does not name row R20" >&2; finish 94; }
[ -d data/derived/feng_relabel_c ] \
  || { echo "refused: the relabelling's output data/derived/feng_relabel_c is not in this checkout" >&2; finish 92; }
O=runs/t1x2_${head:0:7}_$RUN
mkdir -p runs
mkdir "$O" 2> /dev/null || { echo "refused: the output folder $O exists" >&2; finish 93; }
echo "$head" > "$S/head"
step_ok

STEP=tests
$PY -m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp="$S/tmp/pytest" \
  tests/test_t1_fraction2_launcher.py tests/test_t1_fraction2_build.py tests/test_t1_candidates_build.py \
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
build D T1_s2rc $BASE_ARGS || finish $?
got=$(sha256sum "$O/D/T1_s2rc.h5ad" | cut -d' ' -f1)
[ "$got" = "$R13_SHA" ] || { echo "stopped: the default command does not give R13's file" >&2; finish 20; }
step_ok

STEP=B
build B T1_s2rcB $B_ARGS || finish $?
got=$(sha256sum "$O/B/T1_s2rcB.h5ad" | cut -d' ' -f1)
[ "$got" = "$B_SHA" ] || { echo "stopped: the 1.5 command does not give candidate B's file" >&2; finish 21; }
step_ok

for STEP in X2a X2b; do build "$STEP" T1_s2rcB2 $X2_ARGS || finish $?; step_ok; done

STEP=done
(cd "$O" && sha256sum D/*.h5ad B/*.h5ad X2a/*.h5ad X2b/*.h5ad) > "$S/SHA256SUMS" || finish 30
cp "$S/SHA256SUMS" "$O/SHA256SUMS" || finish 30
cmp -s "$O/X2a/T1_s2rcB2.h5ad" "$O/X2b/T1_s2rcB2.h5ad" || finish 31
step_ok
