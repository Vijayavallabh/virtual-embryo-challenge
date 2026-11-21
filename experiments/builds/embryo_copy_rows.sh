#!/usr/bin/env bash
# The embryo's copy-rows file (practice row E50 of experiments/harness/TRIES.md; docs/PLAN.md decision log 2026-11-21
# (l) and (m); the specification experiments/coords/E50_COPY_ROWS_SPEC.md): test-board builds on a DGX campaign
# checkout (~/v/vec-campaign, never the reference checkout). Run from the checkout's root, detached, and wait on the
# exit files with a background poller; never a foreground sleep:
#   bash experiments/builds/embryo_copy_rows.sh <COMMIT> <RUN>
# RUN names this run's own scratch folder ~/v/scratch/embcr/<RUN> and its output folder runs/embcr_<commit7>_<RUN>.
# Neither may exist: a run folder is never reused and nothing is ever deleted here, so a second run takes a new name.
# In order, each step writing ~/v/scratch/embcr/<RUN>/<STEP>.exit:
#   gate   HEAD is COMMIT (the full hash), the tracked tree is clean, docs/S2_RULES.sha256 checks, row E50 is in
#          TRIES.md, row E50's reading is a tracked file of this commit, made at the row's run commit, not void, and says
#          the test file is built (its reading (i) holds, mixture arm), and the code the builder runs (src/, the
#          builder, the count and size records, the lock file) is byte for byte the code that reading ran
#   tests  the construction's and the embryo builder's targeted tests (synthetic stages only), 8 BLAS threads
#   D      the positions file by the default command; it must give POSITIONS_SHA or everything stops
#   CRa CRb  the copy-rows file, twice, into separate folders
#   done   one sha256 list over every built file, a and b compared byte for byte
# A build step also fails (15) unless its stamp says this commit and clean True.
# Builds run at 16 BLAS threads, as the reference builds. Every builder and `vec check` output goes to /dev/null: a
# test build is read by its exit code, its sha256 and its stamp's commit and clean fields only (AGENTS.md), never its
# log, size, cell count or sidecar record. A copy-rows build that exits non-zero (step exit 11) means no file, and the
# reason is not read: a refusal at the board's limits, like a pass, shows one bit about them. The flag
# --copy-rows-test-ok is the captain's order of 2026-11-21 (docs/PLAN.md decision log 2026-11-21 (h)). Nothing here
# uploads anything.
set -u
COMMIT=${1:?usage: embryo_copy_rows.sh COMMIT RUN}
RUN=${2:?usage: embryo_copy_rows.sh COMMIT RUN}
case "$RUN" in
  *[!A-Za-z0-9_-]*) echo "refused: RUN takes letters, digits, _ and - only" >&2; exit 90 ;;
esac
export PATH=$HOME/.local/bin:$PATH
S=$HOME/v/scratch/embcr/$RUN
mkdir -p "$HOME/v/scratch/embcr"
mkdir "$S" 2> /dev/null || { echo "refused: the run folder $S exists (a run folder is never reused)" >&2; exit 91; }
mkdir "$S/tmp"
export NUMPY_MADVISE_HUGEPAGE=0 TMPDIR="$S/tmp"
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
PY="uv run --frozen --no-sync python -I -B"
BOARD=T2:embryo:test_interp
POSITIONS_SHA=a03611a343faf1d37edf06ce0e2bccc86cce6f14d3580b823864e5faf3ba03f9  # set embryo2, the positions file
READING=experiments/coords/results/e50_copy_rows/reading.json
E50_RUN_COMMIT=5fcf6b9ad187e2c5377900426f6dd749c4a1f9d9  # row E50's own commit: the code its reading ran
CODE_PATHS="src experiments/builds/build_test.py experiments/density/results/release/count_decision.json experiments/coords/results/p3 uv.lock pyproject.toml"
BASE_ARGS="--method fgw_shape0 --positions floor_exact --expression mixture"
CR_ARGS="--method fgw_shape0 --positions floor_exact --expression mixture --copy-rows --copy-rows-test-ok"
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
grep -q '^| E50 ' experiments/harness/TRIES.md || { echo "refused: row E50 is not in TRIES.md" >&2; finish 95; }
git ls-files --error-unmatch "$READING" > /dev/null 2>&1 \
  || { echo "refused: row E50 has no committed reading (the file must be tracked at this commit)" >&2; finish 94; }
python3 -I -B -c 'import json,sys; r=json.load(open(sys.argv[1])); sys.exit(0 if (r.get("void") is False and r.get("build_test_file") is True and r.get("i", {}).get("holds") is True and r.get("handed_over_expression") == "mixture" and r.get("commit") == sys.argv[2]) else 1)' "$READING" "$E50_RUN_COMMIT" \
  || { echo "refused: row E50's reading is void, is not the run commit's, or its reading (i) does not hold" >&2; finish 92; }
git diff --quiet "$E50_RUN_COMMIT" HEAD -- $CODE_PATHS \
  || { echo "refused: the builder's code or records differ from those row E50's reading ran" >&2; finish 89; }
O=runs/embcr_${head:0:7}_$RUN
mkdir -p runs
mkdir "$O" 2> /dev/null || { echo "refused: the output folder $O exists" >&2; finish 93; }
echo "$head" > "$S/head"
step_ok

STEP=tests
$PY -m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp="$S/tmp/pytest" \
  tests/test_embryo_copy_rows_launcher.py tests/test_copy_rows.py tests/test_copy_rows_build.py \
  tests/test_e47_floor_exact.py tests/test_hidden_counts.py tests/test_s2_rules.py > "$S/tests.out" 2>&1 || finish 10
step_ok

export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16
stamp() {  # the stamp's commit and clean fields only
  $PY -c 'import json,sys; s=json.load(open(sys.argv[1])); print(s["commit"], s["clean"])' "$1.stamp.json"
}
build() {  # folder name args: the build, then `vec check`; exit 11 build, 13 check, 14 no stamp, 15 stamp not clean
  local d="$O/$1" name=$2 st
  shift 2
  mkdir -p "$d"
  $PY experiments/builds/build_test.py embryo "$@" --out "$d/$name.h5ad" > /dev/null 2>&1 || return 11
  uv run --frozen --no-sync vec check "$d/$name.h5ad" --board "$BOARD" > /dev/null 2>&1 || return 13
  st=$(stamp "$d/$name.h5ad" 2> /dev/null) || return 14
  echo "$(basename "$d")/$name $st" >> "$S/stamps"
  [ "$st" = "$head True" ] || return 15
}
: > "$S/stamps"

STEP=D
build D T2embryo_s2x $BASE_ARGS || finish $?
got=$(sha256sum "$O/D/T2embryo_s2x.h5ad" | cut -d' ' -f1)
[ "$got" = "$POSITIONS_SHA" ] || { echo "stopped: the default command does not give the positions file" >&2; finish 20; }
step_ok

for STEP in CRa CRb; do build "$STEP" T2embryo_s2x_cr $CR_ARGS || finish $?; step_ok; done

STEP=done
(cd "$O" && sha256sum D/*.h5ad CRa/*.h5ad CRb/*.h5ad) > "$S/SHA256SUMS" || finish 30
cp "$S/SHA256SUMS" "$O/SHA256SUMS" || finish 30
cmp -s "$O/CRa/T2embryo_s2x_cr.h5ad" "$O/CRb/T2embryo_s2x_cr.h5ad" || finish 31
step_ok
