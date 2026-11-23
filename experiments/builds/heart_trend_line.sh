#!/usr/bin/env bash
# The heart search line's two files, TS05 and TS10 (row R28 of experiments/harness/TRIES.md, the captain's order of
# 2026-11-23, "start the heart search"; docs/PLAN.md decision log 2026-11-23 (t)): the full trim C's command with row
# H28's trend step at strength 0.5 and 1.0, nothing else changed. Test-board builds on a DGX campaign checkout
# (~/v/vec-campaign, never the reference checkout). Run from the root of the checkout or detached worktree it is called
# in, detached, and wait on the exit files with a background poller; never a foreground sleep:
#   bash experiments/builds/heart_trend_line.sh <COMMIT> <RUN>
# RUN names this run's own scratch folder ~/v/scratch/heartline/<RUN> and its output folder
# runs/heartline_<commit7>_<RUN>. Neither may exist: a run folder is never reused and nothing is ever deleted here, so a
# second run takes a new name. In order, each step writing ~/v/scratch/heartline/<RUN>/<STEP>.exit (the word running
# first, once its checks pass, then its exit code):
#   gate   HEAD is COMMIT (the full hash), the tracked tree is clean, docs/S2_RULES.sha256 checks, and TRIES.md holds
#          exactly one row R28 whose first cell is the captain's order ("ruling, no run") naming
#          `--trend-strength 0.5`, `--trend-strength 1.0` and `--captain-ruling R28` (the builder's own check)
#   tests  the heart search line's and the heart builder's targeted tests (synthetic stages only), 8 BLAS threads
#   C      the identity step: the full trim C by its committed default command, then its half; the file must give
#          C_SHA and its half C_HALF_SHA (set heart3, experiments/builds/results/heart3_3d52b83.SHA256SUMS, the
#          DGX's builds, as experiments/builds/heart_organ.sh holds them), or everything stops
#   TS05a TS05b  TS05 (strength 0.5), twice, into separate folders, each then halved by experiments/builds/halve.py
#                (seed 1000001)
#   TS10a TS10b  TS10 (strength 1.0), likewise
#   done   one sha256 list over every built file; TS05a and TS05b, and TS10a and TS10b, compared byte for byte, file and
#          half
# A build step also fails (15) unless its stamp says this commit and clean True. Builds run at 16 BLAS threads, as the
# reference builds. Every builder, halve.py and `vec check` output goes to /dev/null: a test build is read by its exit
# code, its sha256 and its stamp's commit and clean fields only (AGENTS.md), never its log, size, cell count or
# sidecar record. A build that exits non-zero means no file, and the reason is not read here. Nothing here uploads
# anything; an upload goes by the line's rule and the captain's own act.
set -u
COMMIT=${1:?usage: heart_trend_line.sh COMMIT RUN}
RUN=${2:?usage: heart_trend_line.sh COMMIT RUN}
case "$RUN" in
  *[!A-Za-z0-9_-]*) echo "refused: RUN takes letters, digits, _ and - only" >&2; exit 90 ;;
esac
case "$COMMIT" in
  *[!0-9a-f]*) echo "refused: COMMIT must be a full hash" >&2; exit 90 ;;
esac
[ "${#COMMIT}" = 40 ] || { echo "refused: COMMIT must be a full hash" >&2; exit 90; }
export PATH=$HOME/.local/bin:$PATH
S=$HOME/v/scratch/heartline/$RUN
mkdir -p "$HOME/v/scratch/heartline"
mkdir "$S" 2> /dev/null || { echo "refused: the run folder $S exists (a run folder is never reused)" >&2; exit 91; }
mkdir "$S/tmp"
export NUMPY_MADVISE_HUGEPAGE=0 TMPDIR="$S/tmp"
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
export UV_PROJECT_ENVIRONMENT="$PWD/.venv"
PY="uv run --frozen --no-sync python -I -B"
BOARD=T2:heart:test_extrap
C_SHA=341058344dbd73cd4267f6512bce0f185431ab16b6691d87b28c379f8094c299       # set heart3: T2heart_s2core.h5ad
C_HALF_SHA=04cde6361c04840ee47fff32eb9a20f4781bd649387cc540ed9fdea71e65fa35  # set heart3: T2heart_s2core_half.h5ad
C_ARGS="--entry core_trim --expression copy --z-mode raw"
TS05_ARGS="--entry core_trim --trend-strength 0.5 --z-mode raw --captain-ruling R28"
TS10_ARGS="--entry core_trim --trend-strength 1.0 --z-mode raw --captain-ruling R28"
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
python3 -I -B -c 'import sys
rows = [ln for ln in open(sys.argv[1], encoding="utf-8").read().splitlines() if ln.startswith("| R28 ")]
first = rows[0][2:].split("|", 1)[0] if len(rows) == 1 else ""
need = ("`--trend-strength 0.5`", "`--trend-strength 1.0`", "`--captain-ruling R28`", "ruling, no run")
sys.exit(0 if first.startswith("R28 (the captain\x27s order") and all(s in first for s in need) else 1)' \
  experiments/harness/TRIES.md \
  || { echo "refused: TRIES.md does not hold exactly one row R28 that is the captain's order naming the line" >&2; finish 95; }
O=runs/heartline_${head:0:7}_$RUN
mkdir -p runs
mkdir "$O" 2> /dev/null || { echo "refused: the output folder $O exists" >&2; finish 93; }
echo "$head" > "$S/head"
step_ok

STEP=tests
echo running > "$S/$STEP.exit"
$PY -m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp="$S/tmp/pytest" \
  tests/test_heart_trend_line.py tests/test_heart_trend_line_launcher.py tests/test_build_heart_core_trim.py \
  tests/test_build_heart_iso_trim.py tests/test_heart_candidates.py tests/test_heart_trend_build.py \
  tests/test_h28_heart_trend_strength.py tests/test_halve.py tests/test_hidden_counts.py tests/test_s2_rules.py \
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

for STEP in TS05a TS05b; do
  echo running > "$S/$STEP.exit"
  build "$STEP" T2heart_core_ts05 $TS05_ARGS || finish $?
  step_ok
done
for STEP in TS10a TS10b; do
  echo running > "$S/$STEP.exit"
  build "$STEP" T2heart_core_ts10 $TS10_ARGS || finish $?
  step_ok
done

STEP=done
echo running > "$S/$STEP.exit"
(cd "$O" && sha256sum C/*.h5ad TS05a/*.h5ad TS05b/*.h5ad TS10a/*.h5ad TS10b/*.h5ad) > "$S/SHA256SUMS" || finish 30
cp "$S/SHA256SUMS" "$O/SHA256SUMS" || finish 30
for f in T2heart_core_ts05 T2heart_core_ts05_half; do cmp -s "$O/TS05a/$f.h5ad" "$O/TS05b/$f.h5ad" || finish 31; done
for f in T2heart_core_ts10 T2heart_core_ts10_half; do cmp -s "$O/TS10a/$f.h5ad" "$O/TS10b/$f.h5ad" || finish 32; done
step_ok
