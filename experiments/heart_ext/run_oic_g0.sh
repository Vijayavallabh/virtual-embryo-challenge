#!/usr/bin/env bash
# Row H38's G0 for heart file OIC125 (experiments/heart_ext/oic_g0.py): target-free, not a try, CPU only. On a DGX
# campaign checkout or a detached worktree of it, from the root of a clean checkout at COMMIT, detached; wait on
# done.exit with a background poller, never a foreground sleep.
#   bash experiments/heart_ext/run_oic_g0.sh <COMMIT> <RUN> prep <FETCHED>
#   bash experiments/heart_ext/run_oic_g0.sh <COMMIT> <RUN> g0 <FETCHED> <PREP_RUN>
# FETCHED is the fetch's folder (fetch_oic.py --out; its MANIFEST.json is checked against the pinned sizes before any
# file is read, and every file's sha256 is recomputed as it is read). RUN names this call's scratch folder
# ~/v/scratch/oic_g0/<RUN>, which may not exist (never reused; nothing is deleted here).
# prep: check, tests, ranks (from the references alone, before any outside file is read), units (every fetched file
#   binned on the 12 um grid), regions (the plots for the reader and reader_sheet.json); all in <S>/work. prep may be
#   run again under a new RUN (it is deterministic and counts nothing).
# g0: check, tests, g0, at a later commit that holds the reader's committed verdicts
#   (experiments/heart_ext/oic_reader.json); it reads PREP_RUN's work folder and writes g0.json, g0.txt and the chord
#   files under experiments/heart_ext/results/ (to be committed with G0's statement in row H38's outcome cell). G0 is
#   started once, ever: its exit file ~/v/scratch/oic_g0/g0.exit (88 when it exists).
# Each step writes <S>/<STEP>.exit (running first, then its code); the call ends at the first step that is not 0.
#   check  HEAD is COMMIT, the tracked tree is clean, the scripts are here, row H38 is in TRIES.md, well formed and
#          its outcome cell empty (experiments/harness/h38_row.py check), FETCHED passes fetch_oic.py --check
#          (MANIFEST.json against the pinned sizes and the files on disk), the virtual
#          environment, the T2 stages and (for the Xenium table) pyarrow are there; for g0 also PREP_RUN's work folder
#          with its regions step at 0, the reader's file committed, and no G0 result yet
# OIC_VENV names the virtual environment (default $PWD/.venv, which `uv sync --frozen --inexact --group external`
# fills); OIC_DATA the competition data folder (default data/raw); OIC_JOBS the units step's processes (default 8).
# Cut-offs: neither phase starts after 2026-11-26 03:00 UTC (83); a source not complete by 2026-11-25 21:30 UTC
# (each file's own completion time) is dropped by the units step (fetch_oic.source_status; <S>/sources.json).
# Exit codes of the launcher's own refusals: 83 past the G0 cut-off, 84 the work folder of PREP_RUN is not a
# finished prep, 85 the reader's
# file is not committed, 86 a G0 result exists, 87 pyarrow is not in the environment, 88 G0 was started before, 90 odd
# arguments, 91 the run folder exists, 92 FETCHED fails fetch_oic.py --check, 93 row H38 missing, malformed or with a
# statement, 94 not the checkout's root, 96 HEAD is not COMMIT, 97 no T2 stages, 98 tracked tree not clean, 99 no
# virtual environment. Otherwise the failed step's own code (g0: 3 STOP, 5 VOID, 2 refused).
set -u
COMMIT=${1:?usage: run_oic_g0.sh COMMIT RUN PHASE FETCHED [PREP_RUN]}
RUN=${2:?usage: run_oic_g0.sh COMMIT RUN PHASE FETCHED [PREP_RUN]}
PHASE=${3:?usage: run_oic_g0.sh COMMIT RUN PHASE FETCHED [PREP_RUN]}
FETCHED=${4:?usage: run_oic_g0.sh COMMIT RUN PHASE FETCHED [PREP_RUN]}
case "$PHASE" in
  prep) [ "$#" = 4 ] || { echo "refused: prep takes COMMIT, RUN, prep and FETCHED only" >&2; exit 90; } ;;
  g0) [ "$#" = 5 ] || { echo "refused: g0 takes COMMIT, RUN, g0, FETCHED and PREP_RUN" >&2; exit 90; }
      PREP_RUN=$5 ;;
  *) echo "refused: PHASE is prep or g0" >&2; exit 90 ;;
esac
for v in "$RUN" "${PREP_RUN:-x}"; do
  case "$v" in
    *[!A-Za-z0-9_-]*|"") echo "refused: RUN and PREP_RUN take letters, digits, _ and - only" >&2; exit 90 ;;
  esac
done
case "$COMMIT" in
  *[!0-9a-f]*) echo "refused: COMMIT must be a full hash" >&2; exit 90 ;;
esac
[ "${#COMMIT}" = 40 ] || { echo "refused: COMMIT must be a full hash" >&2; exit 90; }
export PATH=$HOME/.local/bin:$PATH
BASE=$HOME/v/scratch/oic_g0
S=$BASE/$RUN
mkdir -p "$BASE"
mkdir "$S" 2> /dev/null || { echo "refused: the run folder $S exists (a run folder is never reused)" >&2; exit 91; }
mkdir "$S/tmp"
echo running > "$S/done.exit"
STEP=check
OWNED=none
finish() {
  echo "$1" > "$S/$STEP.exit"
  if [ "$OWNED" = "$STEP" ]; then echo "$1" > "$BASE/$STEP.exit"; fi
  echo "$1" > "$S/done.exit"
  exit "$1"
}
begin() { echo running > "$S/$STEP.exit"; }
step_ok() {
  echo 0 > "$S/$STEP.exit"
  if [ "$OWNED" = "$STEP" ]; then echo 0 > "$BASE/$STEP.exit"; fi
}
begin

SCRIPT=experiments/heart_ext/oic_g0.py
RESULTS=experiments/heart_ext/results
G0_CUTOFF=1795662000  # 2026-11-26 03:00 UTC: no G0 PASS by then is STOP, so neither phase starts after it
now=$(date -u +%s)
[ "$now" -lt "$G0_CUTOFF" ] \
  || { echo "refused: it is past 2026-11-26 03:00 UTC, the G0 cut-off (no G0 PASS by then is STOP)" >&2; finish 83; }
READER=experiments/heart_ext/oic_reader.json
head=$(git rev-parse HEAD)
[ "$head" = "$COMMIT" ] || { echo "refused: HEAD is not $COMMIT" >&2; finish 96; }
if [ -n "$(git status --porcelain --untracked-files=no)" ]; then
  echo "refused: tracked files differ from the commit" >&2
  finish 98
fi
[ -f "$SCRIPT" ] || { echo "refused: run from the root of the checkout ($SCRIPT is not here)" >&2; finish 94; }
python3 -I -B experiments/harness/h38_row.py check experiments/harness/TRIES.md 2> "$S/row.err" \
  || { echo "refused: row H38 is missing, malformed or already holds a statement (see $S/row.err)" >&2; finish 93; }
python3 -I -B experiments/heart_ext/fetch_oic.py --sources "$FETCHED" > "$S/sources.json" 2>&1
python3 -I -B experiments/heart_ext/fetch_oic.py --check "$FETCHED" > "$S/manifest.out" 2>&1 \
  || { echo "refused: $FETCHED is not as its MANIFEST.json and the pins say (see $S/manifest.out)" >&2; finish 92; }
VENV=${OIC_VENV:-$PWD/.venv}
[ -x "$VENV/bin/python" ] \
  || { echo "refused: no virtual environment at $VENV (name one with OIC_VENV)" >&2; finish 99; }
DATA=${OIC_DATA:-data/raw}
[ -d "$DATA/T2" ] || { echo "refused: no T2 stages under $DATA (name the data folder with OIC_DATA)" >&2; finish 97; }
if grep -q '"study": "X021"' "$FETCHED/MANIFEST.json" || grep -q 'transcripts.parquet' "$FETCHED/MANIFEST.json"; then
  "$VENV/bin/python" -I -B -c "import pyarrow.parquet" 2> /dev/null \
    || { echo "refused: pyarrow is not in $VENV (uv sync --frozen --inexact --group external)" >&2; finish 87; }
fi
if [ "$PHASE" = g0 ]; then
  W=$BASE/$PREP_RUN/work
  { [ -f "$W/ranks.json" ] && [ -f "$W/units.json" ] && [ -f "$W/regions.json" ] \
    && [ "$(cat "$BASE/$PREP_RUN/regions.exit" 2> /dev/null)" = 0 ]; } \
    || { echo "refused: $W is not a finished prep (ranks, units and regions at exit 0)" >&2; finish 84; }
  git ls-files --error-unmatch "$READER" > /dev/null 2>&1 \
    || { echo "refused: $READER is not committed" >&2; finish 85; }
  for f in g0.json g0.txt chords_E10.5.json chords_E9.5.json chords_E8.75.json; do
    [ -e "$RESULTS/$f" ] && { echo "refused: $RESULTS/$f exists" >&2; finish 86; }
  done
fi
echo "$head $PHASE" > "$S/head"
step_ok

export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 NUMPY_MADVISE_HUGEPAGE=0 TMPDIR="$S/tmp"
export UV_PROJECT_ENVIRONMENT="$VENV"
PY="uv run --frozen --no-sync python -I -B"

STEP=tests
begin
$PY -m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp="$S/tmp/pytest" \
  tests/test_heart_ext_oic.py tests/test_heart_ext_fetch_oic.py tests/test_oic_launchers.py > "$S/tests.log" 2>&1 \
  || finish $?
step_ok

if [ "$PHASE" = prep ]; then
  mkdir "$S/work"
  STEP=ranks
  begin
  $PY "$SCRIPT" ranks --out "$S/work" --data "$DATA" --commit "$COMMIT" > "$S/ranks.out" 2>&1 || finish $?
  step_ok
  STEP=units
  begin
  $PY "$SCRIPT" units --out "$S/work" --fetched "$FETCHED" --jobs "${OIC_JOBS:-8}" --commit "$COMMIT" \
    > "$S/units.out" 2>&1 \
    || finish $?
  step_ok
  STEP=regions
  begin
  $PY "$SCRIPT" regions --out "$S/work" --commit "$COMMIT" > "$S/regions.out" 2>&1 || finish $?
  step_ok
else
  STEP=g0
  begin
  ( set -o noclobber; echo running > "$BASE/g0.exit" ) 2> /dev/null \
    || { echo "refused: G0 was started before ($BASE/g0.exit exists); it is not started twice" >&2; finish 88; }
  OWNED=g0
  echo "$COMMIT $RUN $PREP_RUN" > "$BASE/g0.head"
  $PY "$SCRIPT" g0 --out "$W" --units "$W" --data "$DATA" --commit "$COMMIT" > "$S/g0.out" 2>&1 || finish $?
  step_ok
fi

echo 0 > "$S/done.exit"
exit 0
