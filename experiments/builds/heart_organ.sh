#!/usr/bin/env bash
# The heart's organ-cut file D (row R21 of experiments/harness/TRIES.md, the captain's order of 2026-11-22; its practice
# read is row H29b; docs/S2_RULES.md, amendment of 2026-11-22): test-board builds on a DGX campaign checkout
# (~/v/vec-campaign, never the reference checkout). Run from the root of the checkout or detached worktree it is called
# in, detached, and wait on the exit files with a background poller; never a foreground sleep:
#   bash experiments/builds/heart_organ.sh <COMMIT> <RUN>
# RUN names this run's own scratch folder ~/v/scratch/heartorgan/<RUN> and its output folder
# runs/heartorgan_<commit7>_<RUN>. Neither may exist: a run folder is never reused and nothing is ever deleted here, so
# a second run takes a new name. In order, each step writing ~/v/scratch/heartorgan/<RUN>/<STEP>.exit (the word running
# first, once its checks pass, then its exit code):
#   gate   HEAD is COMMIT (the full hash), the tracked tree is clean, docs/S2_RULES.sha256 checks, rows R21 and H29b are
#          in TRIES.md, row H29b's reading is a tracked file of this commit that is not void, whose gate said PASS (not
#          STOP), made at a run commit this commit descends from, and the code the builder runs (src/, the builder,
#          halve.py and the records and configs they read, boards/, the lock file) is byte for byte the code that
#          reading ran
#   tests  the organ cut's and the heart builder's targeted tests (synthetic stages only), 8 BLAS threads
#   C      the identity step: the full trim C by its committed default command, then its half; the file must give
#          C_SHA and its half C_HALF_SHA (set heart3, experiments/builds/results/heart3_3d52b83.SHA256SUMS, the
#          DGX's builds), or everything stops
#   Da Db  file D, twice, into separate folders, each then halved by experiments/builds/halve.py (seed 1000001)
#   done   one sha256 list over every built file; Da and Db compared byte for byte, file and half
# A build step also fails (15) unless its stamp says this commit and clean True. Builds run at 16 BLAS threads, as the
# reference builds. Every builder, halve.py and `vec check` output goes to /dev/null: a test build is read by its exit
# code, its sha256 and its stamp's commit and clean fields only (AGENTS.md), never its log, size, cell count or
# sidecar record. A build that exits non-zero means no file, and the reason is not read here. Nothing here uploads
# anything; an upload is the captain's own act, on his yes to file D by name.
set -u
COMMIT=${1:?usage: heart_organ.sh COMMIT RUN}
RUN=${2:?usage: heart_organ.sh COMMIT RUN}
case "$RUN" in
  *[!A-Za-z0-9_-]*) echo "refused: RUN takes letters, digits, _ and - only" >&2; exit 90 ;;
esac
case "$COMMIT" in
  *[!0-9a-f]*) echo "refused: COMMIT must be a full hash" >&2; exit 90 ;;
esac
[ "${#COMMIT}" = 40 ] || { echo "refused: COMMIT must be a full hash" >&2; exit 90; }
export PATH=$HOME/.local/bin:$PATH
S=$HOME/v/scratch/heartorgan/$RUN
mkdir -p "$HOME/v/scratch/heartorgan"
mkdir "$S" 2> /dev/null || { echo "refused: the run folder $S exists (a run folder is never reused)" >&2; exit 91; }
mkdir "$S/tmp"
export NUMPY_MADVISE_HUGEPAGE=0 TMPDIR="$S/tmp"
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
export UV_PROJECT_ENVIRONMENT="$PWD/.venv"
PY="uv run --frozen --no-sync python -I -B"
BOARD=T2:heart:test_extrap
C_SHA=341058344dbd73cd4267f6512bce0f185431ab16b6691d87b28c379f8094c299       # set heart3: T2heart_s2core.h5ad
C_HALF_SHA=04cde6361c04840ee47fff32eb9a20f4781bd649387cc540ed9fdea71e65fa35  # set heart3: T2heart_s2core_half.h5ad
READING=experiments/harness/results/h29b_heart_organ/reading.json
CODE_PATHS="src experiments/builds/build_test.py experiments/builds/halve.py experiments/density/results/release/count_decision.json experiments/coords/results/p3 boards configs uv.lock pyproject.toml"
C_ARGS="--entry core_trim --expression copy --z-mode raw"
D_ARGS="--entry organ_scope --expression copy --z-mode raw"
STEP=gate
finish() { echo "$1" > "$S/$STEP.exit"; exit "$1"; }
step_ok() { echo 0 > "$S/$STEP.exit"; }

head=$(git rev-parse HEAD)
[ "$head" = "$COMMIT" ] || { echo "refused: HEAD is not $COMMIT" >&2; finish 96; }
if [ -n "$(git status --porcelain --untracked-files=no)" ]; then
  echo "refused: tracked files differ from the commit (the build stamps must say clean)" >&2
  finish 98
fi
sha256sum -c docs/S2_RULES.sha256 > /dev/null 2>&1 || finish 97
grep -q '^| R21 ' experiments/harness/TRIES.md || { echo "refused: row R21 is not in TRIES.md" >&2; finish 95; }
grep -q '^| H29b ' experiments/harness/TRIES.md || { echo "refused: row H29b is not in TRIES.md" >&2; finish 95; }
git ls-files --error-unmatch "$READING" > /dev/null 2>&1 \
  || { echo "refused: row H29b has no committed reading (the file must be tracked at this commit)" >&2; finish 94; }
RUNC=$(python3 -I -B -c 'import json,re,sys; r=json.load(open(sys.argv[1])); g=r.get("gate") or {}; c=str(r.get("commit")); sys.exit(1) if not (r.get("void") is False and r.get("outcome") == "READ" and g.get("verdict") == "PASS" and g.get("commit") == c and re.fullmatch("[0-9a-f]{40}", c)) else print(c)' "$READING") \
  || { echo "refused: row H29b's reading is void, its gate did not say PASS, or it names no run commit" >&2; finish 92; }
git merge-base --is-ancestor "$RUNC" HEAD 2> /dev/null \
  || { echo "refused: this commit does not descend from the reading's run commit" >&2; finish 88; }
git diff --quiet "$RUNC" HEAD -- $CODE_PATHS \
  || { echo "refused: the builder's code or records differ from those row H29b's reading ran" >&2; finish 89; }
O=runs/heartorgan_${head:0:7}_$RUN
mkdir -p runs
mkdir "$O" 2> /dev/null || { echo "refused: the output folder $O exists" >&2; finish 93; }
echo "$head" > "$S/head"
step_ok

STEP=tests
echo running > "$S/$STEP.exit"
$PY -m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp="$S/tmp/pytest" \
  tests/test_build_heart_organ_scope.py tests/test_organ_scope_mask.py tests/test_t2_chamber_classes.py \
  tests/test_build_heart_core_trim.py tests/test_build_heart_iso_trim.py tests/test_heart_candidates.py \
  tests/test_heart_organ_launchers.py tests/test_halve.py tests/test_hidden_counts.py tests/test_s2_rules.py \
  > "$S/tests.out" 2>&1 || finish 10
step_ok

export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16
stamp() {  # the stamp's commit and clean fields only
  $PY -c 'import json,sys; s=json.load(open(sys.argv[1])); print(s["commit"], s["clean"])' "$1.stamp.json"
}
build() {  # folder name args: the build, its half, `vec check` on both; exit 11 build, 12 halve, 13 check,
           # 14 no stamp, 15 stamp not this commit and clean
  local d="$O/$1" name=$2 st
  shift 2
  mkdir -p "$d"
  $PY experiments/builds/build_test.py heart "$@" --out "$d/$name.h5ad" > /dev/null 2>&1 || return 11
  $PY experiments/builds/halve.py "$d/$name.h5ad" "$d/${name}_half.h5ad" > /dev/null 2>&1 || return 12
  for f in "$d/$name.h5ad" "$d/${name}_half.h5ad"; do
    uv run --frozen --no-sync vec check "$f" --board "$BOARD" > /dev/null 2>&1 || return 13
  done
  st=$(stamp "$d/$name.h5ad" 2> /dev/null) || return 14
  echo "$(basename "$d")/$name $st" >> "$S/stamps"
  [ "$st" = "$head True" ] || return 15
}
: > "$S/stamps"

STEP=C
echo running > "$S/$STEP.exit"
build C T2heart_s2core $C_ARGS || finish $?
got=$(sha256sum "$O/C/T2heart_s2core.h5ad" | cut -d' ' -f1)
[ "$got" = "$C_SHA" ] || { echo "stopped: the full trim's default command does not give set heart3's file" >&2; finish 20; }
got=$(sha256sum "$O/C/T2heart_s2core_half.h5ad" | cut -d' ' -f1)
[ "$got" = "$C_HALF_SHA" ] || { echo "stopped: the full trim's half does not give set heart3's half" >&2; finish 21; }
step_ok

for STEP in Da Db; do
  echo running > "$S/$STEP.exit"
  build "$STEP" T2heart_organ $D_ARGS || finish $?
  step_ok
done

STEP=done
echo running > "$S/$STEP.exit"
(cd "$O" && sha256sum C/*.h5ad Da/*.h5ad Db/*.h5ad) > "$S/SHA256SUMS" || finish 30
cp "$S/SHA256SUMS" "$O/SHA256SUMS" || finish 30
for f in T2heart_organ T2heart_organ_half; do cmp -s "$O/Da/$f.h5ad" "$O/Db/$f.h5ad" || finish 31; done
step_ok
