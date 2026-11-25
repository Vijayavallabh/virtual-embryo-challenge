#!/usr/bin/env bash
# File OXP125 (row R56 of experiments/harness/TRIES.md, the captain's order of 2026-11-24, ruling, no run: his named
# exception on the heart's scope and size for this file and his A15 exception, bound to row H36's read): OTS125's
# command with E10.5's run profile removed (build_test.py heart --entry organ_scope --trend-strength 1.25 --run-profile
# <sha256> --z-mode raw --captain-ruling R56, then halve.py at seed 1000001). Test-board builds on a DGX campaign
# checkout (~/v/vec-campaign, never the reference checkout). Run from the root of the checkout or detached worktree it
# is called in, detached, and wait on the exit files with a background poller; never a foreground sleep:
#   bash experiments/builds/heart_xpc_build.sh <COMMIT> <RUN>
# The row is fixed (R56; the builder pins it). The profile's sha256 is read by code from row H36's outcome cell (the
# test profile's PASS statement, experiments/harness/h36_row.py passed) and must be the sha256 of the committed
# experiments/heart_profile/results/xpc_profile.json. RUN names this run's own scratch folder ~/v/scratch/heartxpc/<RUN>
# and its output folder runs/heartxpc_<commit7>_<RUN>. Neither may exist: a run folder is never reused and nothing is
# ever deleted here. In order, each step writing ~/v/scratch/heartxpc/<RUN>/<STEP>.exit (the word running first, then
# its exit code):
#   gate    HEAD is COMMIT (the full hash), the tracked tree is clean, docs/S2_RULES.sha256 checks, TRIES.md holds
#           exactly one row R56 whose first cell is the captain's order ("ruling, no run") naming `--entry
#           organ_scope`, `--trend-strength 1.25`, `--run-profile` and `--captain-ruling R56` (the builder's own
#           check); row H36's outcome cell, read by code, holds the test profile's PASS statement with its sha256,
#           the practice read's PASSED AT PRACTICE before it and nothing that ends the line
#           (experiments/harness/h36_row.py passed); experiments/heart_profile/results/xpc_profile.json is committed
#           and is the file of that sha256; and the output folder is new. Nothing is built otherwise
#   tests   the run-profile builder's, OTS15's, file D's two further strengths', the two heart search lines' and the
#           heart builder's targeted tests (synthetic stages only), 8 BLAS threads
#   D       the first identity step: file D by its committed default command, then its half; the file must give D_SHA
#           and its half D_HALF_SHA (experiments/builds/results/heartorgan_72b9662_d1.SHA256SUMS), or everything stops
#   OTS125  row R49's OTS125 by its committed command, then its half: OTS125_SHA and OTS125_HALF_SHA
#           (experiments/builds/results/heartf_ee53e77_f1.SHA256SUMS, both of its builds), or everything stops
#           (the two identity steps come before any OXP125 build)
#   OXP125a OXP125b  OXP125, twice, into separate folders, each then halved by experiments/builds/halve.py (seed
#                  1000001)
#   done    one sha256 list over every built file; OXP125a and OXP125b compared byte for byte, file and half
# A build step also fails (15) unless its stamp says this commit and clean True. Builds run at 16 BLAS threads, as the
# reference builds. Every builder, halve.py and `vec check` output goes to /dev/null: a test build is read by its exit
# code, its sha256 and its stamp's commit and clean fields only (AGENTS.md), never its log, size, cell count or
# sidecar record. A build that exits non-zero means no file, and the reason is not read here. Nothing here uploads
# anything; an upload goes by the line's rule and the captain's own act.
# Exit codes of the gate: 90 odd arguments, 91 the run folder exists, 92 the profile file is not committed or not the
# file of row H36's sha256, 93 the output folder exists, 94 row H36's outcome cell does not say PASS, 95 no row R56
# that orders the file, 96 HEAD is not COMMIT, 97 docs/S2_RULES.sha256 fails, 98 tracked tree not clean.
set -u
COMMIT=${1:?usage: heart_xpc_build.sh COMMIT RUN}
RUN=${2:?usage: heart_xpc_build.sh COMMIT RUN}
case "$RUN" in
  *[!A-Za-z0-9_-]*) echo "refused: RUN takes letters, digits, _ and - only" >&2; exit 90 ;;
esac
case "$COMMIT" in
  *[!0-9a-f]*) echo "refused: COMMIT must be a full hash" >&2; exit 90 ;;
esac
[ "${#COMMIT}" = 40 ] || { echo "refused: COMMIT must be a full hash" >&2; exit 90; }
[ "$#" = 2 ] || { echo "refused: COMMIT and RUN only (the row, R56, is fixed)" >&2; exit 90; }
export PATH=$HOME/.local/bin:$PATH
S=$HOME/v/scratch/heartxpc/$RUN
mkdir -p "$HOME/v/scratch/heartxpc"
mkdir "$S" 2> /dev/null || { echo "refused: the run folder $S exists (a run folder is never reused)" >&2; exit 91; }
mkdir "$S/tmp"
export NUMPY_MADVISE_HUGEPAGE=0 TMPDIR="$S/tmp"
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
export UV_PROJECT_ENVIRONMENT="$PWD/.venv"
PY="uv run --frozen --no-sync python -I -B"
BOARD=T2:heart:test_extrap
PROFILE=experiments/heart_profile/results/xpc_profile.json
D_SHA=76e21cd2c754362c043ec32ca863a1e1b7f2cfde549e0366d6104354bfffe18b            # heartorgan d1: T2heart_organ.h5ad
D_HALF_SHA=1b4342c151c19dfc827d4b7a310c6a0fa8908f1e369836ebffb22c27c93b4f00       # heartorgan d1: its half
OTS125_SHA=79d0b813db758ce88bda22ed1f8c8b58a560214360bfbfc2795c7eee372419da       # heartf f1: OTS125's file
OTS125_HALF_SHA=adc3243665949f9802056216d86443e2ea0d45d07534013f99e26afb7e57ade0  # heartf f1: OTS125's half
D_ARGS="--entry organ_scope --expression copy --z-mode raw"
OTS125_ARGS="--entry organ_scope --trend-strength 1.25 --z-mode raw --captain-ruling R49"
STEP=gate
finish() { echo "$1" > "$S/$STEP.exit"; exit "$1"; }
step_ok() { echo 0 > "$S/$STEP.exit"; }
echo running > "$S/$STEP.exit"

head=$(git rev-parse HEAD)
[ "$head" = "$COMMIT" ] || { echo "refused: HEAD is not $COMMIT" >&2; finish 96; }
if [ -n "$(git status --porcelain --untracked-files=no)" ]; then
  echo "refused: tracked files differ from the commit (the build stamps must say clean)" >&2
  finish 98
fi
sha256sum -c docs/S2_RULES.sha256 > /dev/null 2>&1 || finish 97
python3 -I -B -c 'import sys
rows = [ln for ln in open(sys.argv[1], encoding="utf-8").read().splitlines() if ln.startswith("| R56 ")]
first = rows[0][2:].split("|", 1)[0] if len(rows) == 1 else ""
need = ("`--entry organ_scope`", "`--trend-strength 1.25`", "`--run-profile`", "`--captain-ruling R56`",
        "ruling, no run")
sys.exit(0 if first.startswith("R56 (the captain\x27s order") and all(s in first for s in need) else 1)' \
  experiments/harness/TRIES.md \
  || { echo "refused: TRIES.md does not hold exactly one row R56 that is the captain's order naming the file" >&2; finish 95; }
SHA=$(python3 -I -B experiments/harness/h36_row.py passed experiments/harness/TRIES.md 2> "$S/h36_row.err") \
  || { echo "refused: row H36's outcome cell does not say PASS (its reasons are in $S/h36_row.err)" >&2; finish 94; }
case "$SHA" in
  *[!0-9a-f]*|"") echo "refused: row H36's outcome cell gives no sha256" >&2; finish 94 ;;
esac
[ "${#SHA}" = 64 ] || { echo "refused: row H36's outcome cell gives no sha256" >&2; finish 94; }
[ -f "$PROFILE" ] || { echo "refused: $PROFILE is not here" >&2; finish 92; }
got=$(sha256sum "$PROFILE" | cut -d' ' -f1)
[ "$got" = "$SHA" ] || { echo "refused: $PROFILE is not the file of row H36's sha256" >&2; finish 92; }
git ls-files --error-unmatch "$PROFILE" > /dev/null 2>&1 \
  || { echo "refused: $PROFILE is not committed" >&2; finish 92; }
OXP125_ARGS="--entry organ_scope --trend-strength 1.25 --run-profile $SHA --z-mode raw --captain-ruling R56"
echo "$SHA" > "$S/profile_sha256"
O=runs/heartxpc_${head:0:7}_$RUN
mkdir -p runs
mkdir "$O" 2> /dev/null || { echo "refused: the output folder $O exists" >&2; finish 93; }
echo "$head $SHA" > "$S/head"
step_ok

STEP=tests
echo running > "$S/$STEP.exit"
$PY -m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp="$S/tmp/pytest" \
  tests/test_heart_xpc_build.py tests/test_heart_trend_line_g.py tests/test_heart_trend_line_f.py \
  tests/test_heart_trend_line_d.py tests/test_heart_trend_line.py tests/test_build_heart_organ_scope.py \
  tests/test_heart_trend_build.py tests/test_h28_heart_trend_strength.py tests/test_xpc_table.py \
  tests/test_halve.py tests/test_hidden_counts.py tests/test_s2_rules.py \
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

STEP=OTS125
echo running > "$S/$STEP.exit"
build OTS125 T2heart_organ_ts125 $OTS125_ARGS || finish $?
got=$(sha256sum "$O/OTS125/T2heart_organ_ts125.h5ad" | cut -d' ' -f1)
[ "$got" = "$OTS125_SHA" ] || { echo "stopped: row R49's OTS125 command does not give its recorded file" >&2; finish 28; }
got=$(sha256sum "$O/OTS125/T2heart_organ_ts125_half.h5ad" | cut -d' ' -f1)
[ "$got" = "$OTS125_HALF_SHA" ] || { echo "stopped: row R49's OTS125 half does not give its recorded half" >&2; finish 29; }
step_ok

for STEP in OXP125a OXP125b; do
  echo running > "$S/$STEP.exit"
  build "$STEP" T2heart_organ_ts125_xp $OXP125_ARGS || finish $?
  step_ok
done

STEP=done
echo running > "$S/$STEP.exit"
(cd "$O" && sha256sum D/*.h5ad OTS125/*.h5ad OXP125a/*.h5ad OXP125b/*.h5ad) > "$S/SHA256SUMS" || finish 30
cp "$S/SHA256SUMS" "$O/SHA256SUMS" || finish 30
for f in T2heart_organ_ts125_xp T2heart_organ_ts125_xp_half; do
  cmp -s "$O/OXP125a/$f.h5ad" "$O/OXP125b/$f.h5ad" || finish 31
done
step_ok
