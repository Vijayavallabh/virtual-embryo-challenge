#!/usr/bin/env bash
# File D's three further strengths, files OTS025, OTS075 and OTS125 (the design round of 2026-11-24, its FAST proposal
# for the heart; its ruling row, "R-F" there, is the ROW given here): the organ cut D's command with --expression copy
# replaced by row H28's trend step at strength 0.25, 0.75 and 1.25, nothing else changed (OTS05's draw, organ-cut pool,
# count rule, size procedure, seed and half). Test-board builds on a DGX campaign checkout (~/v/vec-campaign, never the
# reference checkout). Run from the root of the checkout or detached worktree it is called in, detached, and wait on
# the exit files with a background poller; never a foreground sleep:
#   bash experiments/builds/heart_trend_line_f.sh <COMMIT> <RUN> <ROW>
# ROW is the ruling row's id in experiments/harness/TRIES.md (R and its number, e.g. R44). RUN names this run's own
# scratch folder ~/v/scratch/heartf/<RUN> and its output folder runs/heartf_<commit7>_<RUN>. Neither may exist: a run
# folder is never reused and nothing is ever deleted here, so a second run takes a new name. In order, each step
# writing ~/v/scratch/heartf/<RUN>/<STEP>.exit (the word running first, once its checks pass, then its exit code):
#   gate   HEAD is COMMIT (the full hash), the tracked tree is clean, docs/S2_RULES.sha256 checks, and TRIES.md holds
#          exactly one row ROW whose first cell is the captain's order ("ruling, no run") naming
#          `--entry organ_scope`, `--trend-strength 0.25`, `--trend-strength 0.75`, `--trend-strength 1.25` and
#          `--captain-ruling ROW` (the builder's own check)
#   tests  file D's three further strengths', the two heart search lines' and the heart builder's targeted tests
#          (synthetic stages only), 8 BLAS threads
#   D      the first identity step: file D by its committed default command, then its half; the file must give D_SHA
#          and its half D_HALF_SHA (experiments/builds/results/heartorgan_72b9662_d1.SHA256SUMS and
#          heartlined_0d47a67_od1.SHA256SUMS, the DGX's builds), or everything stops
#   OTS05  the second identity step: row R29's OTS05 by its committed command, then its half; the file must give
#          OTS05_SHA and its half OTS05_HALF_SHA (experiments/builds/results/heartlined_0d47a67_od1.SHA256SUMS, the
#          DGX's builds), or everything stops
#   OTS10  the third identity step: row R29's OTS10 likewise, OTS10_SHA and OTS10_HALF_SHA (the same list)
#   OTS025a OTS025b  OTS025 (strength 0.25 on D), twice, into separate folders, each then halved by
#                    experiments/builds/halve.py (seed 1000001)
#   OTS075a OTS075b  OTS075 (strength 0.75 on D), likewise
#   OTS125a OTS125b  OTS125 (strength 1.25 on D), likewise
#   done   one sha256 list over every built file; OTS025a and OTS025b, OTS075a and OTS075b, and OTS125a and OTS125b,
#          compared byte for byte, file and half
# A build step also fails (15) unless its stamp says this commit and clean True. Builds run at 16 BLAS threads, as the
# reference builds. Every builder, halve.py and `vec check` output goes to /dev/null: a test build is read by its exit
# code, its sha256 and its stamp's commit and clean fields only (AGENTS.md), never its log, size, cell count or
# sidecar record. A build that exits non-zero means no file, and the reason is not read here. Nothing here uploads
# anything; an upload goes by the line's rule and the captain's own act.
set -u
COMMIT=${1:?usage: heart_trend_line_f.sh COMMIT RUN ROW}
RUN=${2:?usage: heart_trend_line_f.sh COMMIT RUN ROW}
ROW=${3:?usage: heart_trend_line_f.sh COMMIT RUN ROW}
case "$RUN" in
  *[!A-Za-z0-9_-]*) echo "refused: RUN takes letters, digits, _ and - only" >&2; exit 90 ;;
esac
case "$COMMIT" in
  *[!0-9a-f]*) echo "refused: COMMIT must be a full hash" >&2; exit 90 ;;
esac
[ "${#COMMIT}" = 40 ] || { echo "refused: COMMIT must be a full hash" >&2; exit 90; }
[[ "$ROW" =~ ^R[0-9]+[a-z]?$ ]] || { echo "refused: ROW is a ruling row id, R and its number (e.g. R44)" >&2; exit 90; }
export PATH=$HOME/.local/bin:$PATH
S=$HOME/v/scratch/heartf/$RUN
mkdir -p "$HOME/v/scratch/heartf"
mkdir "$S" 2> /dev/null || { echo "refused: the run folder $S exists (a run folder is never reused)" >&2; exit 91; }
mkdir "$S/tmp"
export NUMPY_MADVISE_HUGEPAGE=0 TMPDIR="$S/tmp"
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
export UV_PROJECT_ENVIRONMENT="$PWD/.venv"
PY="uv run --frozen --no-sync python -I -B"
BOARD=T2:heart:test_extrap
D_SHA=76e21cd2c754362c043ec32ca863a1e1b7f2cfde549e0366d6104354bfffe18b           # heartorgan d1: T2heart_organ.h5ad
D_HALF_SHA=1b4342c151c19dfc827d4b7a310c6a0fa8908f1e369836ebffb22c27c93b4f00      # heartorgan d1: its half
OTS05_SHA=c3b5fc9c7b47d0d0aafc4e5876bac7f06988cb2d628a0030b4b643806da8bceb       # heartlined od1: OTS05's file
OTS05_HALF_SHA=f918e76731aac18eda3b46c8a35c18f26ae00a1a6a1db36d33e7515704cc43ea  # heartlined od1: OTS05's half
OTS10_SHA=efaaabe0b744d82dc40135ed9269a6baf8938c9cdd83f85fb3917c0b102e7a2b       # heartlined od1: OTS10's file
OTS10_HALF_SHA=ee9e0e68d62e8d1448f1dbbbce1146dc5fb4877e94d77fc289fa82744acb06a9  # heartlined od1: OTS10's half
D_ARGS="--entry organ_scope --expression copy --z-mode raw"
OTS05_ARGS="--entry organ_scope --trend-strength 0.5 --z-mode raw --captain-ruling R29"
OTS10_ARGS="--entry organ_scope --trend-strength 1.0 --z-mode raw --captain-ruling R29"
OTS025_ARGS="--entry organ_scope --trend-strength 0.25 --z-mode raw --captain-ruling $ROW"
OTS075_ARGS="--entry organ_scope --trend-strength 0.75 --z-mode raw --captain-ruling $ROW"
OTS125_ARGS="--entry organ_scope --trend-strength 1.25 --z-mode raw --captain-ruling $ROW"
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
row = sys.argv[2]
rows = [ln for ln in open(sys.argv[1], encoding="utf-8").read().splitlines() if ln.startswith("| " + row + " ")]
first = rows[0][2:].split("|", 1)[0] if len(rows) == 1 else ""
need = ("`--entry organ_scope`", "`--trend-strength 0.25`", "`--trend-strength 0.75`", "`--trend-strength 1.25`",
        "`--captain-ruling " + row + "`", "ruling, no run")
sys.exit(0 if first.startswith(row + " (the captain\x27s order") and all(s in first for s in need) else 1)' \
  experiments/harness/TRIES.md "$ROW" \
  || { echo "refused: TRIES.md does not hold exactly one row $ROW that is the captain's order naming the line" >&2; finish 95; }
O=runs/heartf_${head:0:7}_$RUN
mkdir -p runs
mkdir "$O" 2> /dev/null || { echo "refused: the output folder $O exists" >&2; finish 93; }
echo "$head" > "$S/head"
echo "$ROW" > "$S/row"
step_ok

STEP=tests
echo running > "$S/$STEP.exit"
$PY -m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp="$S/tmp/pytest" \
  tests/test_heart_trend_line_f.py tests/test_heart_trend_line_f_launcher.py tests/test_heart_trend_line_d.py \
  tests/test_heart_trend_line_d_launcher.py tests/test_heart_trend_line.py tests/test_heart_trend_line_launcher.py \
  tests/test_build_heart_organ_scope.py tests/test_heart_organ_launchers.py tests/test_build_heart_core_trim.py \
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

STEP=D
echo running > "$S/$STEP.exit"
build D T2heart_organ $D_ARGS || finish $?
got=$(sha256sum "$O/D/T2heart_organ.h5ad" | cut -d' ' -f1)
[ "$got" = "$D_SHA" ] || { echo "stopped: the organ cut's default command does not give file D" >&2; finish 20; }
got=$(sha256sum "$O/D/T2heart_organ_half.h5ad" | cut -d' ' -f1)
[ "$got" = "$D_HALF_SHA" ] || { echo "stopped: the organ cut's half does not give file D's half" >&2; finish 21; }
step_ok

STEP=OTS05
echo running > "$S/$STEP.exit"
build OTS05 T2heart_organ_ts05 $OTS05_ARGS || finish $?
got=$(sha256sum "$O/OTS05/T2heart_organ_ts05.h5ad" | cut -d' ' -f1)
[ "$got" = "$OTS05_SHA" ] || { echo "stopped: row R29's OTS05 command does not give its recorded file" >&2; finish 22; }
got=$(sha256sum "$O/OTS05/T2heart_organ_ts05_half.h5ad" | cut -d' ' -f1)
[ "$got" = "$OTS05_HALF_SHA" ] || { echo "stopped: row R29's OTS05 half does not give its recorded half" >&2; finish 23; }
step_ok

STEP=OTS10
echo running > "$S/$STEP.exit"
build OTS10 T2heart_organ_ts10 $OTS10_ARGS || finish $?
got=$(sha256sum "$O/OTS10/T2heart_organ_ts10.h5ad" | cut -d' ' -f1)
[ "$got" = "$OTS10_SHA" ] || { echo "stopped: row R29's OTS10 command does not give its recorded file" >&2; finish 24; }
got=$(sha256sum "$O/OTS10/T2heart_organ_ts10_half.h5ad" | cut -d' ' -f1)
[ "$got" = "$OTS10_HALF_SHA" ] || { echo "stopped: row R29's OTS10 half does not give its recorded half" >&2; finish 25; }
step_ok

for STEP in OTS025a OTS025b; do
  echo running > "$S/$STEP.exit"
  build "$STEP" T2heart_organ_ts025 $OTS025_ARGS || finish $?
  step_ok
done
for STEP in OTS075a OTS075b; do
  echo running > "$S/$STEP.exit"
  build "$STEP" T2heart_organ_ts075 $OTS075_ARGS || finish $?
  step_ok
done
for STEP in OTS125a OTS125b; do
  echo running > "$S/$STEP.exit"
  build "$STEP" T2heart_organ_ts125 $OTS125_ARGS || finish $?
  step_ok
done

STEP=done
echo running > "$S/$STEP.exit"
(cd "$O" && sha256sum D/*.h5ad OTS05/*.h5ad OTS10/*.h5ad OTS025a/*.h5ad OTS025b/*.h5ad OTS075a/*.h5ad \
  OTS075b/*.h5ad OTS125a/*.h5ad OTS125b/*.h5ad) > "$S/SHA256SUMS" || finish 30
cp "$S/SHA256SUMS" "$O/SHA256SUMS" || finish 30
for f in T2heart_organ_ts025 T2heart_organ_ts025_half; do cmp -s "$O/OTS025a/$f.h5ad" "$O/OTS025b/$f.h5ad" || finish 31; done
for f in T2heart_organ_ts075 T2heart_organ_ts075_half; do cmp -s "$O/OTS075a/$f.h5ad" "$O/OTS075b/$f.h5ad" || finish 32; done
for f in T2heart_organ_ts125 T2heart_organ_ts125_half; do cmp -s "$O/OTS125a/$f.h5ad" "$O/OTS125b/$f.h5ad" || finish 33; done
step_ok
