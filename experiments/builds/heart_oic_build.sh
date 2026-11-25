#!/usr/bin/env bash
# File OIC125 (row R70 of experiments/harness/TRIES.md, the captain's order of 2026-11-25, ruling, no run: his A15
# exception for OIC125's outside sources and his named exception on the heart's scope and size for OIC125; read by row
# H38): OTK125's command with the in-situ chords added (build_test.py heart --entry organ_scope --trend-strength 1.25
# --trend-form insitu --insitu-chords <sha256> --z-mode raw --captain-ruling R70, then halve.py at seed 1000001).
# Test-board builds on a DGX campaign checkout (never the reference checkout), from the root of the checkout or
# detached worktree it is called in, detached; wait on the exit files with a background poller, never a foreground
# sleep:
#   bash experiments/builds/heart_oic_build.sh <COMMIT> <RUN> <FETCHED>
# FETCHED is the fetch's folder: it must pass fetch_oic.py --check, and chords_E10.5.json must name its files at the
# sha256 it records.
# The row is fixed (R70; the builder pins it). The chords' sha256 is read by code from row H38's outcome cell (G0's
# PASS statement, then the read's "OUTCOME H38: PASS.", experiments/harness/h38_row.py passed) and must be the sha256 of
# the committed experiments/heart_ext/results/chords_E10.5.json. RUN names this run's scratch folder
# ~/v/scratch/heartoic/<RUN> and its output folder runs/heartoic_<commit7>_<RUN>; neither may exist (never reused;
# nothing is deleted here). In order, each step writing ~/v/scratch/heartoic/<RUN>/<STEP>.exit (running first, then its
# exit code):
#   gate     HEAD is COMMIT (the full hash), the tracked tree is clean, docs/S2_RULES.sha256 checks, TRIES.md holds
#            exactly one row R70 whose first cell is the captain's order ("ruling, no run") naming OIC125 and
#            `--captain-ruling R70`; row H38's outcome cell holds G0's PASS statement and the read's PASS after it and
#            nothing that ends the line; chords_E10.5.json is committed and is the file of that sha256; the output
#            folder is new
#   tests    the in-situ form's, the builder option's, the read's and the launchers' tests and the heart builder's
#            targeted tests (synthetic stages only), 8 BLAS threads
#   OTK125   the identity rebuild of the banked heart file: OTK125 by its committed command, then its half; the file
#            must give OTK125_SHA and its half OTK125_HALF_SHA (experiments/builds/results/hearth_8a8bc3a_h2.SHA256SUMS,
#            both of its builds), or everything stops before any OIC125 build
#   OIC125a OIC125b  OIC125, twice, into separate folders, each then halved by experiments/builds/halve.py (seed
#            1000001)
#   done     one sha256 list over every built file; OIC125a and OIC125b compared byte for byte, file and half
# A build step also fails (15) unless its stamp says this commit and clean True. Builds run at 16 BLAS threads, as the
# reference builds. Every builder, halve.py and `vec check` output goes to /dev/null: a test build is read by its exit
# code, its sha256 and its stamp's commit and clean fields only (AGENTS.md), never its log, size, cell count or
# sidecar record. Nothing here uploads anything; an upload goes by the line's rule and the captain's own act.
# Exit codes of the gate: 89 FETCHED or the chords fail the manifest check, 90 odd arguments, 91 the run folder
# exists, 92 the chord file is not committed or not the file
# of row H38's sha256, 93 the output folder exists, 94 row H38's outcome cell does not say PASS, 95 no row R70 that
# orders the file, 96 HEAD is not COMMIT, 97 docs/S2_RULES.sha256 fails, 98 tracked tree not clean.
set -u
COMMIT=${1:?usage: heart_oic_build.sh COMMIT RUN FETCHED}
RUN=${2:?usage: heart_oic_build.sh COMMIT RUN FETCHED}
FETCHED=${3:?usage: heart_oic_build.sh COMMIT RUN FETCHED}
case "$RUN" in
  *[!A-Za-z0-9_-]*) echo "refused: RUN takes letters, digits, _ and - only" >&2; exit 90 ;;
esac
case "$COMMIT" in
  *[!0-9a-f]*) echo "refused: COMMIT must be a full hash" >&2; exit 90 ;;
esac
[ "${#COMMIT}" = 40 ] || { echo "refused: COMMIT must be a full hash" >&2; exit 90; }
[ "$#" = 3 ] || { echo "refused: COMMIT, RUN and FETCHED only (the row, R70, is fixed)" >&2; exit 90; }
export PATH=$HOME/.local/bin:$PATH
S=$HOME/v/scratch/heartoic/$RUN
mkdir -p "$HOME/v/scratch/heartoic"
mkdir "$S" 2> /dev/null || { echo "refused: the run folder $S exists (a run folder is never reused)" >&2; exit 91; }
mkdir "$S/tmp"
export NUMPY_MADVISE_HUGEPAGE=0 TMPDIR="$S/tmp"
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
export UV_PROJECT_ENVIRONMENT="$PWD/.venv"
PY="uv run --frozen --no-sync python -I -B"
BOARD=T2:heart:test_extrap
CHORDS=experiments/heart_ext/results/chords_E10.5.json
OTK125_SHA=d4d8d75f8bed48813bf143f5d7726828e5c30cdb27bd708b6fc2d9a455d21893       # hearth h2: OTK125's file
OTK125_HALF_SHA=9835fa5d47d7021760beb2056ce44c2c846ded6fc31dacf64bbb2d9f0f7d1b54  # hearth h2: OTK125's half
OTK125_ARGS="--entry organ_scope --trend-strength 1.25 --trend-form category --z-mode raw --captain-ruling R55"
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
rows = [ln for ln in open(sys.argv[1], encoding="utf-8").read().splitlines() if ln.startswith("| R70 ")]
first = rows[0][2:].split("|", 1)[0] if len(rows) == 1 else ""
need = ("OIC125", "`--captain-ruling R70`", "ruling, no run")
sys.exit(0 if first.startswith("R70 (the captain\x27s order") and all(s in first for s in need) else 1)' \
  experiments/harness/TRIES.md \
  || { echo "refused: TRIES.md does not hold exactly one row R70 that is the captain's order naming OIC125" >&2
       finish 95; }
SHA=$(python3 -I -B experiments/harness/h38_row.py passed experiments/harness/TRIES.md 2> "$S/h38_row.err") \
  || { echo "refused: row H38's outcome cell does not say PASS (its reasons are in $S/h38_row.err)" >&2; finish 94; }
case "$SHA" in
  *[!0-9a-f]*|"") echo "refused: row H38's outcome cell gives no sha256" >&2; finish 94 ;;
esac
[ "${#SHA}" = 64 ] || { echo "refused: row H38's outcome cell gives no sha256" >&2; finish 94; }
[ -f "$CHORDS" ] || { echo "refused: $CHORDS is not here" >&2; finish 92; }
got=$(sha256sum "$CHORDS" | cut -d' ' -f1)
[ "$got" = "$SHA" ] || { echo "refused: $CHORDS is not the file of row H38's sha256" >&2; finish 92; }
git ls-files --error-unmatch "$CHORDS" > /dev/null 2>&1 || { echo "refused: $CHORDS is not committed" >&2; finish 92; }
python3 -I -B experiments/heart_ext/fetch_oic.py --check "$FETCHED" --chords "$CHORDS" > "$S/manifest.out" 2>&1 \
  || { echo "refused: $FETCHED or $CHORDS fail fetch_oic.py --check (see $S/manifest.out)" >&2; finish 89; }
OIC125_ARGS="--entry organ_scope --trend-strength 1.25 --trend-form insitu --insitu-chords $SHA --z-mode raw"
OIC125_ARGS="$OIC125_ARGS --captain-ruling R70"
echo "$SHA" > "$S/chords_sha256"
O=runs/heartoic_${head:0:7}_$RUN
mkdir -p runs
mkdir "$O" 2> /dev/null || { echo "refused: the output folder $O exists" >&2; finish 93; }
echo "$head $SHA" > "$S/head"
step_ok

STEP=tests
echo running > "$S/$STEP.exit"
$PY -m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp="$S/tmp/pytest" \
  tests/test_heart_oic_build.py tests/test_h38_heart_insitu.py tests/test_oic_launchers.py \
  tests/test_heart_trend_line_h.py tests/test_build_heart_organ_scope.py tests/test_heart_trend_build.py \
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
  $PY experiments/builds/halve.py "$d/$name.h5ad" "$d/${name}_half.h5ad" --seed 1000001 > /dev/null 2>&1 || return 12
  for f in "$d/$name.h5ad" "$d/${name}_half.h5ad"; do
    uv run --frozen --no-sync vec check "$f" --board "$BOARD" > /dev/null 2>&1 || return 13
  done
  st=$(stamp "$d/$name.h5ad" 2> /dev/null) || return 14
  echo "$(basename "$d")/$name $st" >> "$S/stamps"
  [ "$st" = "$head True" ] || return 15
}
: > "$S/stamps"

STEP=OTK125
echo running > "$S/$STEP.exit"
build OTK125 T2heart_organ_ts125_tk $OTK125_ARGS || finish $?
got=$(sha256sum "$O/OTK125/T2heart_organ_ts125_tk.h5ad" | cut -d' ' -f1)
[ "$got" = "$OTK125_SHA" ] || { echo "stopped: OTK125's command does not give its recorded file" >&2; finish 20; }
got=$(sha256sum "$O/OTK125/T2heart_organ_ts125_tk_half.h5ad" | cut -d' ' -f1)
[ "$got" = "$OTK125_HALF_SHA" ] || { echo "stopped: OTK125's half does not give its recorded half" >&2; finish 21; }
step_ok

for STEP in OIC125a OIC125b; do
  echo running > "$S/$STEP.exit"
  build "$STEP" T2heart_organ_ts125_ti $OIC125_ARGS || finish $?
  step_ok
done

STEP=done
echo running > "$S/$STEP.exit"
(cd "$O" && sha256sum OTK125/*.h5ad OIC125a/*.h5ad OIC125b/*.h5ad) > "$S/SHA256SUMS" || finish 30
cp "$S/SHA256SUMS" "$O/SHA256SUMS" || finish 30
for f in T2heart_organ_ts125_ti T2heart_organ_ts125_ti_half; do
  cmp -s "$O/OIC125a/$f.h5ad" "$O/OIC125b/$f.h5ad" || finish 31
done
step_ok
