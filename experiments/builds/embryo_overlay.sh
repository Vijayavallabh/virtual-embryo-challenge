#!/usr/bin/env bash
# Row R19's embryo overlay file (experiments/harness/TRIES.md; docs/S2_RULES.md, the amendment of 2026-11-21): check
# builds on a DGX campaign checkout (~/v/vec-campaign, never the reference checkout). Run from the checkout's root,
# detached, and wait on the exit files with a background poller; never a foreground sleep:
#   bash experiments/builds/embryo_overlay.sh <COMMIT> <RUN>
# RUN names this run's own scratch folder ~/v/scratch/embov/<RUN> and its output folder runs/embov_<commit7>_<RUN>.
# Neither may exist: a run folder is never reused and nothing is ever deleted here, so a second run takes a new name.
# In order, each step writing ~/v/scratch/embov/<RUN>/<STEP>.exit:
#   gate   HEAD is COMMIT, the tracked tree is clean, docs/S2_RULES.sha256 checks, row R19 is in TRIES.md and
#          docs/S2_RULES.md names R19 (the record comes before any build meant for upload)
#   tests  the embryo builder's targeted tests (synthetic stages only), 8 BLAS threads
#   D      the positions file by the default command; it must give POSITIONS_SHA or everything stops
#   OVa OVb  the overlay file, twice, into separate folders
#   done   one sha256 list over every built file, a and b compared byte for byte
# Builds run at 16 BLAS threads, as the reference builds. Every builder and `vec check` output goes to /dev/null: a
# test build is read by its exit code, its sha256 and its stamp's commit and clean fields only (AGENTS.md), never its
# log, size or cell count. The flag --shape-interp-test-ok is the captain's yes of 2026-11-21 (docs/PLAN.md decision
# log 2026-11-21 (j)). Nothing here uploads anything.
set -u
COMMIT=${1:?usage: embryo_overlay.sh COMMIT RUN}
RUN=${2:?usage: embryo_overlay.sh COMMIT RUN}
case "$RUN" in
  *[!A-Za-z0-9_-]*) echo "refused: RUN takes letters, digits, _ and - only" >&2; exit 90 ;;
esac
export PATH=$HOME/.local/bin:$PATH
S=$HOME/v/scratch/embov/$RUN
mkdir -p "$HOME/v/scratch/embov"
mkdir "$S" 2> /dev/null || { echo "refused: the run folder $S exists (a run folder is never reused)" >&2; exit 91; }
mkdir "$S/tmp"
export NUMPY_MADVISE_HUGEPAGE=0 TMPDIR="$S/tmp"
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
PY="uv run --frozen --no-sync python -I -B"
BOARD=T2:embryo:test_interp
POSITIONS_SHA=a03611a343faf1d37edf06ce0e2bccc86cce6f14d3580b823864e5faf3ba03f9  # set embryo2, the positions file
BASE_ARGS="--method fgw_shape0 --positions floor_exact --expression mixture"
OV_ARGS="--method fgw_shape0 --positions floor_exact --expression mixture --shape-interp overlay --shape-interp-test-ok"
STEP=gate
finish() { echo "$1" > "$S/$STEP.exit"; exit "$1"; }
step_ok() { echo 0 > "$S/$STEP.exit"; }

head=$(git rev-parse HEAD)
case "$head" in "$COMMIT"*) ;; *) echo "refused: HEAD is not $COMMIT" >&2; finish 96 ;; esac
if [ -n "$(git status --porcelain --untracked-files=no)" ]; then
  echo "refused: tracked files differ from the commit (the build stamps must say clean)" >&2
  finish 98
fi
sha256sum -c docs/S2_RULES.sha256 > /dev/null 2>&1 || finish 97
grep -q '^| R19 ' experiments/harness/TRIES.md || { echo "refused: row R19 is not in TRIES.md" >&2; finish 95; }
grep -q 'R19' docs/S2_RULES.md || { echo "refused: no docs/S2_RULES.md amendment names R19" >&2; finish 94; }
O=runs/embov_${head:0:7}_$RUN
mkdir -p runs
mkdir "$O" 2> /dev/null || { echo "refused: the output folder $O exists" >&2; finish 93; }
echo "$head" > "$S/head"
step_ok

STEP=tests
$PY -m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp="$S/tmp/pytest" \
  tests/test_embryo_overlay.py tests/test_shape_interp.py tests/test_shape_interp_build.py \
  tests/test_e47_floor_exact.py tests/test_hidden_counts.py tests/test_s2_rules.py > "$S/tests.out" 2>&1 || finish 10
step_ok

export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16
stamp() {  # the stamp's commit and clean fields only
  $PY -c 'import json,sys; s=json.load(open(sys.argv[1])); print(s["commit"], s["clean"])' "$1.stamp.json"
}
build() {  # folder name args: the build, then `vec check`; exit 11 build, 13 check
  local d="$O/$1" name=$2
  shift 2
  mkdir -p "$d"
  $PY experiments/builds/build_test.py embryo "$@" --out "$d/$name.h5ad" > /dev/null 2>&1 || return 11
  uv run --frozen --no-sync vec check "$d/$name.h5ad" --board "$BOARD" > /dev/null 2>&1 || return 13
  echo "$(basename "$d")/$name $(stamp "$d/$name.h5ad")" >> "$S/stamps"
}
: > "$S/stamps"

STEP=D
build D T2embryo_s2x $BASE_ARGS || finish $?
got=$(sha256sum "$O/D/T2embryo_s2x.h5ad" | cut -d' ' -f1)
[ "$got" = "$POSITIONS_SHA" ] || { echo "stopped: the default command does not give the positions file" >&2; finish 20; }
step_ok

for STEP in OVa OVb; do build "$STEP" T2embryo_s2x_ov $OV_ARGS || finish $?; step_ok; done

STEP=done
(cd "$O" && sha256sum D/*.h5ad OVa/*.h5ad OVb/*.h5ad) > "$S/SHA256SUMS"
cp "$S/SHA256SUMS" "$O/SHA256SUMS"
cmp -s "$O/OVa/T2embryo_s2x_ov.h5ad" "$O/OVb/T2embryo_s2x_ov.h5ad" || finish 31
step_ok
