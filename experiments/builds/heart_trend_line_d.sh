#!/usr/bin/env bash
# The heart search line's step on file D, files OTS05 and OTS10 (row R29 of experiments/harness/TRIES.md, the captain's
# order "i want it too"; docs/PLAN.md decision log 2026-11-23 (w)): the organ cut D's command with --expression copy
# replaced by row H28's trend step at strength 0.5 and 1.0, nothing else changed. Test-board builds on a DGX campaign
# checkout (~/v/vec-campaign, never the reference checkout). Run from the root of the checkout or detached worktree it
# is called in, detached, and wait on the exit files with a background poller; never a foreground sleep:
#   bash experiments/builds/heart_trend_line_d.sh <COMMIT> <RUN>
# RUN names this run's own scratch folder ~/v/scratch/heartlined/<RUN> and its output folder
# runs/heartlined_<commit7>_<RUN>. Neither may exist: a run folder is never reused and nothing is ever deleted here, so
# a second run takes a new name. In order, each step writing ~/v/scratch/heartlined/<RUN>/<STEP>.exit (the word running
# first, once its checks pass, then its exit code):
#   gate   HEAD is COMMIT (the full hash), the tracked tree is clean, docs/S2_RULES.sha256 checks, and TRIES.md holds
#          exactly one row R29 whose first cell is the captain's order ("ruling, no run") naming `--entry organ_scope`,
#          `--trend-strength 0.5`, `--trend-strength 1.0` and `--captain-ruling R29` (the builder's own check)
#   tests  the two heart search lines' and the heart builder's targeted tests (synthetic stages only), 8 BLAS threads
#   D      the first identity step: file D by its committed default command, then its half; the file must give D_SHA
#          and its half D_HALF_SHA (experiments/builds/results/heartorgan_72b9662_d1.SHA256SUMS, the DGX's builds),
#          or everything stops
#   T      the second identity step: row R28's TS05 by its committed command, then its half; the file must give
#          TS05_SHA and its half TS05_HALF_SHA (experiments/builds/results/heartline_d41cd14_ts1.SHA256SUMS, the
#          DGX's builds), or everything stops
#   OTS05a OTS05b  OTS05 (strength 0.5 on D), twice, into separate folders, each then halved by
#                  experiments/builds/halve.py (seed 1000001)
#   OTS10a OTS10b  OTS10 (strength 1.0 on D), likewise
#   done   one sha256 list over every built file; OTS05a and OTS05b, and OTS10a and OTS10b, compared byte for byte,
#          file and half
# A build step also fails (15) unless its stamp says this commit and clean True. Builds run at 16 BLAS threads, as the
# reference builds. Every builder, halve.py and `vec check` output goes to /dev/null: a test build is read by its exit
# code, its sha256 and its stamp's commit and clean fields only (AGENTS.md), never its log, size, cell count or
# sidecar record. A build that exits non-zero means no file, and the reason is not read here. Nothing here uploads
# anything; an upload goes by the line's rule and the captain's own act.
set -u
COMMIT=${1:?usage: heart_trend_line_d.sh COMMIT RUN}
RUN=${2:?usage: heart_trend_line_d.sh COMMIT RUN}
case "$RUN" in
  *[!A-Za-z0-9_-]*) echo "refused: RUN takes letters, digits, _ and - only" >&2; exit 90 ;;
esac
case "$COMMIT" in
  *[!0-9a-f]*) echo "refused: COMMIT must be a full hash" >&2; exit 90 ;;
esac
[ "${#COMMIT}" = 40 ] || { echo "refused: COMMIT must be a full hash" >&2; exit 90; }
export PATH=$HOME/.local/bin:$PATH
S=$HOME/v/scratch/heartlined/$RUN
mkdir -p "$HOME/v/scratch/heartlined"
mkdir "$S" 2> /dev/null || { echo "refused: the run folder $S exists (a run folder is never reused)" >&2; exit 91; }
mkdir "$S/tmp"
export NUMPY_MADVISE_HUGEPAGE=0 TMPDIR="$S/tmp"
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
export UV_PROJECT_ENVIRONMENT="$PWD/.venv"
PY="uv run --frozen --no-sync python -I -B"
BOARD=T2:heart:test_extrap
D_SHA=76e21cd2c754362c043ec32ca863a1e1b7f2cfde549e0366d6104354bfffe18b          # heartorgan d1: T2heart_organ.h5ad
D_HALF_SHA=1b4342c151c19dfc827d4b7a310c6a0fa8908f1e369836ebffb22c27c93b4f00     # heartorgan d1: T2heart_organ_half.h5ad
TS05_SHA=a6edffa9ab27bd8740db54db652832cbf56278f09154d7ce33d72b884db7d741       # heartline ts1: T2heart_core_ts05.h5ad
TS05_HALF_SHA=b6b8155be9b5469a2dbc8943a59df567f52c7cbe013bbf0b55e4da706fda27b7  # heartline ts1: its half
D_ARGS="--entry organ_scope --expression copy --z-mode raw"
TS05_ARGS="--entry core_trim --trend-strength 0.5 --z-mode raw --captain-ruling R28"
OTS05_ARGS="--entry organ_scope --trend-strength 0.5 --z-mode raw --captain-ruling R29"
OTS10_ARGS="--entry organ_scope --trend-strength 1.0 --z-mode raw --captain-ruling R29"
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
rows = [ln for ln in open(sys.argv[1], encoding="utf-8").read().splitlines() if ln.startswith("| R29 ")]
first = rows[0][2:].split("|", 1)[0] if len(rows) == 1 else ""
need = ("`--entry organ_scope`", "`--trend-strength 0.5`", "`--trend-strength 1.0`", "`--captain-ruling R29`",
        "ruling, no run")
sys.exit(0 if first.startswith("R29 (the captain\x27s order") and all(s in first for s in need) else 1)' \
  experiments/harness/TRIES.md \
  || { echo "refused: TRIES.md does not hold exactly one row R29 that is the captain's order naming the line" >&2; finish 95; }
O=runs/heartlined_${head:0:7}_$RUN
mkdir -p runs
mkdir "$O" 2> /dev/null || { echo "refused: the output folder $O exists" >&2; finish 93; }
echo "$head" > "$S/head"
step_ok

STEP=tests
echo running > "$S/$STEP.exit"
$PY -m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp="$S/tmp/pytest" \
  tests/test_heart_trend_line_d.py tests/test_heart_trend_line_d_launcher.py tests/test_heart_trend_line.py \
  tests/test_heart_trend_line_launcher.py tests/test_build_heart_organ_scope.py tests/test_heart_organ_launchers.py \
  tests/test_build_heart_core_trim.py tests/test_build_heart_iso_trim.py tests/test_heart_candidates.py \
  tests/test_heart_trend_build.py tests/test_h28_heart_trend_strength.py tests/test_halve.py \
  tests/test_hidden_counts.py tests/test_s2_rules.py \
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

STEP=D
echo running > "$S/$STEP.exit"
build D T2heart_organ $D_ARGS || finish $?
got=$(sha256sum "$O/D/T2heart_organ.h5ad" | cut -d' ' -f1)
[ "$got" = "$D_SHA" ] || { echo "stopped: the organ cut's default command does not give file D" >&2; finish 20; }
got=$(sha256sum "$O/D/T2heart_organ_half.h5ad" | cut -d' ' -f1)
[ "$got" = "$D_HALF_SHA" ] || { echo "stopped: the organ cut's half does not give file D's half" >&2; finish 21; }
step_ok

STEP=T
echo running > "$S/$STEP.exit"
build T T2heart_core_ts05 $TS05_ARGS || finish $?
got=$(sha256sum "$O/T/T2heart_core_ts05.h5ad" | cut -d' ' -f1)
[ "$got" = "$TS05_SHA" ] || { echo "stopped: row R28's TS05 command does not give its recorded file" >&2; finish 22; }
got=$(sha256sum "$O/T/T2heart_core_ts05_half.h5ad" | cut -d' ' -f1)
[ "$got" = "$TS05_HALF_SHA" ] || { echo "stopped: row R28's TS05 half does not give its recorded half" >&2; finish 23; }
step_ok

for STEP in OTS05a OTS05b; do
  echo running > "$S/$STEP.exit"
  build "$STEP" T2heart_organ_ts05 $OTS05_ARGS || finish $?
  step_ok
done
for STEP in OTS10a OTS10b; do
  echo running > "$S/$STEP.exit"
  build "$STEP" T2heart_organ_ts10 $OTS10_ARGS || finish $?
  step_ok
done

STEP=done
echo running > "$S/$STEP.exit"
(cd "$O" && sha256sum D/*.h5ad T/*.h5ad OTS05a/*.h5ad OTS05b/*.h5ad OTS10a/*.h5ad OTS10b/*.h5ad) > "$S/SHA256SUMS" \
  || finish 30
cp "$S/SHA256SUMS" "$O/SHA256SUMS" || finish 30
for f in T2heart_organ_ts05 T2heart_organ_ts05_half; do cmp -s "$O/OTS05a/$f.h5ad" "$O/OTS05b/$f.h5ad" || finish 31; done
for f in T2heart_organ_ts10 T2heart_organ_ts10_half; do cmp -s "$O/OTS10a/$f.h5ad" "$O/OTS10b/$f.h5ad" || finish 32; done
step_ok
