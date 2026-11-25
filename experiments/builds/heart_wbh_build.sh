#!/usr/bin/env bash
# File WBH (row R69 of experiments/harness/TRIES.md, the captain's order of 2026-11-25, ruling, no run: the heart's
# world-blur hedge on the banked half, built only after its gate row H37 reads PASS, only with --captain-ruling R69):
# OTK125's half (build_test.py heart --entry organ_scope --trend-strength 1.25 --trend-form category --z-mode raw
# --captain-ruling R55, then halve.py at 1000001) re-measured by experiments/builds/world_blur.py IN OUT --prior
# same_stage --blocks 32 --seed 1000001 --captain-ruling R69. Test-board builds on a DGX campaign checkout (never the
# reference checkout). Run from the root of the checkout or detached worktree it is called in, detached, and wait on
# the exit files with a background poller; never a foreground sleep:
#   bash experiments/builds/heart_wbh_build.sh <COMMIT> <RUN>
# The row is fixed (R69; world_blur.py pins it). RUN names this run's own scratch folder ~/v/scratch/heartwbh/<RUN> and
# its output folder runs/heartwbh_<commit7>_<RUN>. Neither may exist: a run folder is never reused and nothing is ever
# deleted here. In order, each step writing ~/v/scratch/heartwbh/<RUN>/<STEP>.exit (the word running first, then its
# exit code):
#   gate    HEAD is COMMIT (the full hash), the tracked tree is clean, docs/S2_RULES.sha256 checks; TRIES.md holds
#           exactly one row R69 that is the captain's order naming File WBH, row H37 and `--captain-ruling R69`, and
#           row H37's outcome cell, read by code, holds the read's PASS line (experiments/harness/h37_row.py ruling
#           and passed); the gate's records committed at COMMIT (experiments/harness/results/h37_heart_world_blur:
#           G0's g0.json and the read's reading.json) name one run commit, G0's from a clean tree, and the read says
#           PASS, not void; that run commit is in COMMIT's history and no file the gate read (GATE_FILES) differs
#           between it and COMMIT, so the file is built from the code the gate read; the output folder is new
#   tests   the blur's, the file's builder's, the gate's and both launchers' tests and the heart builder's targeted
#           tests (synthetic stages only), 8 BLAS threads
#   OTK125  the identity step: the banked file by its committed command, then its half; both must give their recorded
#           sha256 (experiments/builds/results/hearth_8a8bc3a_h2.SHA256SUMS), or everything stops
#   WBHa WBHb  world_blur.py on that half, twice, into separate folders
#   done    one sha256 list over every built file; WBHa and WBHb compared byte for byte
# A build step also fails (15) unless its stamp says this commit and clean True. Builds run at 16 BLAS threads, as the
# reference builds. Every builder, halve.py, world_blur.py and `vec check` output goes to /dev/null: a test build is
# read by its exit code, its sha256 and its stamp's commit and clean fields only (AGENTS.md), never its log, size, cell
# count or sidecar record. A build that exits non-zero means no file, and the reason is not read here. Nothing here
# uploads anything; the file goes in only on the captain's word at the time.
# Exit codes: 90 the arguments, 91 the run folder exists, 96 HEAD is not COMMIT, 98 tracked tree not clean, 97 the
# rules pin, 95 no row R69 that orders the file, 92 row H37 missing or malformed or its outcome cell without the
# read's PASS line, 87 the outcome is VOID, 86 the outcome is STOP, 84 the gate's records do not name one clean run
# commit in COMMIT's history that says PASS, 85 a file the gate read differs between that commit and COMMIT, 93 the
# output folder exists, 10 the tests, 11 a build failed, 12 halve.py, 13 vec check, 14 no stamp, 15 stamp not this
# commit and clean, 16 world_blur.py failed, 20 21 the identity step is not its recorded bytes (file, half), 30 the
# hash list, 31 the two builds differ.
set -u
set -o pipefail
COMMIT=${1:?usage: heart_wbh_build.sh COMMIT RUN}
RUN=${2:?usage: heart_wbh_build.sh COMMIT RUN}
case "$RUN" in
  *[!A-Za-z0-9_-]*) echo "refused: RUN takes letters, digits, _ and - only" >&2; exit 90 ;;
esac
case "$COMMIT" in
  *[!0-9a-f]*) echo "refused: COMMIT must be a full hash" >&2; exit 90 ;;
esac
[ "${#COMMIT}" = 40 ] || { echo "refused: COMMIT must be a full hash" >&2; exit 90; }
[ "$#" = 2 ] || { echo "refused: COMMIT and RUN only (the row, R69, is fixed)" >&2; exit 90; }
export PATH=$HOME/.local/bin:$PATH
S=$HOME/v/scratch/heartwbh/$RUN
mkdir -p "$HOME/v/scratch/heartwbh"
mkdir "$S" 2> /dev/null || { echo "refused: the run folder $S exists (a run folder is never reused)" >&2; exit 91; }
mkdir "$S/tmp"
export NUMPY_MADVISE_HUGEPAGE=0 TMPDIR="$S/tmp"
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
export UV_PROJECT_ENVIRONMENT="$PWD/.venv"
PY="uv run --frozen --no-sync python -I -B"
BOARD=T2:heart:test_extrap
ROWS=experiments/harness/h37_row.py
OTK125_SHA=d4d8d75f8bed48813bf143f5d7726828e5c30cdb27bd708b6fc2d9a455d21893       # hearth h2: OTK125's file
OTK125_HALF_SHA=9835fa5d47d7021760beb2056ce44c2c846ded6fc31dacf64bbb2d9f0f7d1b54  # hearth h2: OTK125's half
OTK125_ARGS="--entry organ_scope --trend-strength 1.25 --trend-form category --z-mode raw --captain-ruling R55"
WBH_ARGS="--prior same_stage --blocks 32 --seed 1000001 --captain-ruling R69"
GATE_RECORDS=experiments/harness/results/h37_heart_world_blur
# the files the gate read: the blur, the file's builder, the heart builder whose half it blurs, halve.py, every vec
# module and config, the size procedure and the gate's scripts as they import each other
GATE_FILES="src configs experiments/builds/build_test.py experiments/builds/halve.py experiments/builds/world_blur.py"
GATE_FILES="$GATE_FILES experiments/coords/size_procedure.py experiments/s2_hedge/payoff.py"
GATE_FILES="$GATE_FILES experiments/harness/h25_heart_forms.py experiments/harness/h26_heart_scope_forms.py"
GATE_FILES="$GATE_FILES experiments/harness/h27_heart_forecast_size.py"
GATE_FILES="$GATE_FILES experiments/harness/h28_heart_trend_strength.py experiments/harness/h29b_heart_organ.py"
GATE_FILES="$GATE_FILES experiments/harness/h34_heart_organ_trend.py experiments/harness/h35_heart_trend_forms.py"
GATE_FILES="$GATE_FILES experiments/harness/h37_heart_world_blur.py experiments/harness/h37_row.py"
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
python3 -I -B "$ROWS" ruling experiments/harness/TRIES.md \
  || { echo "refused: TRIES.md does not hold exactly one row R69 that is the captain's order naming file WBH" >&2
       finish 95; }
LABEL=$(python3 -I -B "$ROWS" passed experiments/harness/TRIES.md 2> "$S/h37_row.err")
code=$?
case "$code" in
  0) ;;
  5) echo "refused: row H37's outcome is VOID; nothing is built" >&2; finish 87 ;;
  6) echo "refused: row H37's outcome is STOP; nothing is built" >&2; finish 86 ;;
  *) echo "refused: row H37 is missing or malformed, or its outcome cell holds no PASS line of the read" >&2
     finish 92 ;;
esac
GATED=$(python3 -I -B -c '
import json
import re
import sys
def load(path):
    try:
        with open(path, encoding="utf-8") as f:
            rec = json.load(f)
    except (OSError, ValueError):
        sys.exit(1)
    if not isinstance(rec, dict):
        sys.exit(1)
    return rec
g0, read = load(sys.argv[1]), load(sys.argv[2])
commit = g0.get("commit")
prov = g0.get("provenance") if isinstance(g0.get("provenance"), dict) else {}
ok = (isinstance(commit, str) and re.fullmatch("[0-9a-f]{40}", commit) is not None and prov.get("commit") == commit
      and prov.get("tree_clean") is True and g0.get("verdict") == "PASS" and read.get("commit") == commit
      and read.get("void") is False and read.get("outcome") == "PASS" and read.get("label") == sys.argv[3])
if not ok:
    sys.exit(1)
print(commit)
' "$GATE_RECORDS/g0/g0.json" "$GATE_RECORDS/reading.json" "$LABEL") \
  || { echo "refused: the gate's records do not name one clean run commit whose read says PASS" >&2; finish 84; }
git merge-base --is-ancestor "$GATED" "$COMMIT" > /dev/null 2>&1 \
  || { echo "refused: the gate's run commit $GATED is not in the history of $COMMIT" >&2; finish 84; }
git diff --quiet "$GATED" "$COMMIT" -- $GATE_FILES \
  || { echo "refused: a file the gate read differs between its run commit $GATED and $COMMIT" >&2; finish 85; }
O=runs/heartwbh_${head:0:7}_$RUN
mkdir -p runs
mkdir "$O" 2> /dev/null || { echo "refused: the output folder $O exists" >&2; finish 93; }
echo "$head" > "$S/head"
echo "$GATED" > "$S/gated"
step_ok

STEP=tests
echo running > "$S/$STEP.exit"
$PY -m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp="$S/tmp/pytest" \
  tests/test_world_blur.py tests/test_heart_wbh_build.py tests/test_h37_heart_world_blur.py \
  tests/test_run_h37_launcher.py tests/test_heart_trend_line_h.py tests/test_build_heart_organ_scope.py \
  tests/test_heart_trend_build.py tests/test_halve.py tests/test_hidden_counts.py tests/test_s2_rules.py \
  > "$S/tests.out" 2>&1 || finish 10
step_ok

export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16
stamp() {  # the stamp's commit and clean fields only
  $PY -c 'import json,sys; s=json.load(open(sys.argv[1])); print(s["commit"], s["clean"])' "$1.stamp.json"
}
: > "$S/stamps"

STEP=OTK125
echo running > "$S/$STEP.exit"
mkdir -p "$O/OTK125"
$PY experiments/builds/build_test.py heart $OTK125_ARGS --out "$O/OTK125/T2heart_organ_ts125_tk.h5ad" \
  > /dev/null 2>&1 || finish 11
$PY experiments/builds/halve.py "$O/OTK125/T2heart_organ_ts125_tk.h5ad" "$O/OTK125/T2heart_organ_ts125_tk_half.h5ad" \
  > /dev/null 2>&1 || finish 12
st=$(stamp "$O/OTK125/T2heart_organ_ts125_tk.h5ad" 2> /dev/null) || finish 14
echo "OTK125/T2heart_organ_ts125_tk $st" >> "$S/stamps"
[ "$st" = "$head True" ] || finish 15
got=$(sha256sum "$O/OTK125/T2heart_organ_ts125_tk.h5ad" | cut -d' ' -f1)
[ "$got" = "$OTK125_SHA" ] || { echo "stopped: OTK125's command does not give its recorded file" >&2; finish 20; }
got=$(sha256sum "$O/OTK125/T2heart_organ_ts125_tk_half.h5ad" | cut -d' ' -f1)
[ "$got" = "$OTK125_HALF_SHA" ] || { echo "stopped: OTK125's half is not its recorded half" >&2; finish 21; }
step_ok

for STEP in WBHa WBHb; do
  echo running > "$S/$STEP.exit"
  mkdir -p "$O/$STEP"
  f="$O/$STEP/T2heart_organ_ts125_tk_wb_half.h5ad"
  $PY experiments/builds/world_blur.py "$O/OTK125/T2heart_organ_ts125_tk_half.h5ad" "$f" $WBH_ARGS \
    > /dev/null 2>&1 || finish 16
  uv run --frozen --no-sync vec check "$f" --board "$BOARD" > /dev/null 2>&1 || finish 13
  st=$(stamp "$f" 2> /dev/null) || finish 14
  echo "$STEP/T2heart_organ_ts125_tk_wb_half $st" >> "$S/stamps"
  [ "$st" = "$head True" ] || finish 15
  step_ok
done

STEP=done
echo running > "$S/$STEP.exit"
(cd "$O" && sha256sum OTK125/*.h5ad WBHa/*.h5ad WBHb/*.h5ad) > "$S/SHA256SUMS" || finish 30
cp "$S/SHA256SUMS" "$O/SHA256SUMS" || finish 30
cmp -s "$O/WBHa/T2heart_organ_ts125_tk_wb_half.h5ad" "$O/WBHb/T2heart_organ_ts125_tk_wb_half.h5ad" || finish 31
step_ok
