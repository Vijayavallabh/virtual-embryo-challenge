#!/bin/bash
# Release-day run (2026-10-20) of the coordinate checks and the two pre-registered size procedures.
# Run on the DGX from the repo root, inside a tmux session (several minutes; the log goes to runs/logs/):
#   tmux new -d -s vec-coords-p3 'bash experiments/coords/run_p3.sh'
#
# Options:
#   --heart-e85 PATH     heart E8.5 file   (default data/raw/T2/E8.5.h5ad)
#   --heart-e105 PATH    heart E10.5 file  (default data/raw/T2/E10.5.h5ad)
#   --embryo-e75 PATH    embryo E7.5 file  (default data/raw/T2/E7.5.h5ad)
#   --t2-dir DIR         directory of the released T2 files (default data/raw/T2)
#   --timestamp TEXT     recorded in the size JSONs (default: the current UTC time)
#   --out DIR            full outputs (default runs/coords/p3)
#   --results DIR        aggregates copied here (default experiments/coords/results/p3)
#   --allow-uncommitted  run although the lineage map or the procedure code has uncommitted changes (the JSON logs
#                        them either way)
#   --rehearsal          code-path check before the release: released files stand in for the new stages, outputs go
#                        to runs/coords/p3_rehearsal, forecasts are redacted and nothing is copied to results
# The default paths follow the names of the released T2 files. The data manifest needs the team token, which stays on
# the local machine and was not used here, so check the real names after `vec data fetch` and pass them if they differ.
#
# Steps (stops at the first failure, with a message):
#   1. preflight (p3_preflight.py): every input exists; every label is in configs/lineage_map.yaml (map new labels by
#      the rule in its header and commit first); the lineage map and the procedure code are committed
#   2. vec data qc on the new files
#   3. characterise.py and zpitch_fit.py on the new files; handedness.py: landmark verdicts for E8.5 and E10.5, mirror
#      Dice for the new files, and proper-vs-improper shape alignment of E7.5 against E7.25 and E8.0, of E8.5 against
#      E8.25_late and E8.75, and of E10.5 against E9.5, with split-half calibration of every stage in those pairs
#      (a pair counts only when both stages are chiral: see usable_verdict)
#   4. heart size procedure (E8.25_late, E8.5, E8.75, E9.5, E10.5 -> E12.5) and embryo size procedure (E6.75, E7.25,
#      E7.5, E8.0 -> E7.75): one JSON each, with provenance, in --out
#   5. copy the aggregates (CSV, JSON, TXT; no per-cell rows) to --results. Then, on the local machine:
#        rsync -az -e ssh PrakashDGX_H2:v/virtual-embryo-challenge/experiments/coords/results/p3/ \
#          experiments/coords/results/p3/
#      Commit both JSONs before any test submission. The T2 build applies the heart JSON through
#      vec.t2_build.size_or_hold (approved rule 3, 2026-10-21): "valid" false keeps E10.5's size. The embryo JSON
#      goes through apply_size_procedure, which refuses "valid" false (docs/P3_RUNBOOK.md section 5).
set -euo pipefail
cd "$(dirname "$0")/../.."

T2=data/raw/T2
E85="" E105="" E75=""
TS=$(date -u +%Y-%m-%dT%H:%M:%SZ)
OUT="" RESULTS="" ALLOW="" REHEARSAL=0
while [ $# -gt 0 ]; do
  case "$1" in
    --heart-e85) E85=$2; shift 2 ;;
    --heart-e105) E105=$2; shift 2 ;;
    --embryo-e75) E75=$2; shift 2 ;;
    --t2-dir) T2=$2; shift 2 ;;
    --timestamp) TS=$2; shift 2 ;;
    --out) OUT=$2; shift 2 ;;
    --results) RESULTS=$2; shift 2 ;;
    --allow-uncommitted) ALLOW=--allow-uncommitted; shift ;;
    --rehearsal) REHEARSAL=1; shift ;;
    -h|--help) sed -n '2,36p' "$0"; exit 0 ;;
    *) echo "run_p3.sh: unknown option $1 (see --help)" >&2; exit 1 ;;
  esac
done

REDACT=""
if [ "$REHEARSAL" = 1 ]; then
  # stand-ins with the right roles (heart, heart, embryo); the outputs are code checks only
  E85=${E85:-$T2/E8.25_late.h5ad} E105=${E105:-$T2/E9.5.h5ad} E75=${E75:-$T2/E7.25.h5ad}
  OUT=${OUT:-runs/coords/p3_rehearsal} ALLOW=--allow-uncommitted REDACT=--redact-forecasts
else
  E85=${E85:-$T2/E8.5.h5ad} E105=${E105:-$T2/E10.5.h5ad} E75=${E75:-$T2/E7.5.h5ad}
  OUT=${OUT:-runs/coords/p3} RESULTS=${RESULTS:-experiments/coords/results/p3}
fi
E625=$T2/E6.75.h5ad E725=$T2/E7.25.h5ad E80=$T2/E8.0.h5ad
E825=$T2/E8.25_late.h5ad E875=$T2/E8.75.h5ad E95=$T2/E9.5.h5ad

mkdir -p "$OUT" runs/logs
SUFFIX=""
if [ "$REHEARSAL" = 1 ]; then SUFFIX=-rehearsal; fi
LOG=runs/logs/vec-coords-p3$SUFFIX-$(date -u +%Y%m%dT%H%M%SZ).log
exec > >(tee -a "$LOG") 2>&1
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-8} OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-8} \
  MKL_NUM_THREADS=${MKL_NUM_THREADS:-8} COLUMNS=${COLUMNS:-250}
echo "== run_p3.sh $(date -u +%FT%TZ) rehearsal=$REHEARSAL timestamp=$TS commit=$(git rev-parse HEAD) log=$LOG"
echo "   new files: E8.5=$E85 E10.5=$E105 E7.5=$E75"

step() { echo; echo "== $*"; }
fail() { echo "run_p3.sh STOPPED: $*" >&2; exit "${2:-1}"; }

step "1. preflight"
set +e
uv run python experiments/coords/p3_preflight.py $ALLOW "$E85" "$E105" "$E75" "$E625" "$E725" "$E80" "$E825" \
  "$E875" "$E95"
rc=$?
set -e
case $rc in
  0) ;;
  2) fail "missing input files (listed above): download them, or pass the real paths with --heart-e85, --heart-e105, --embryo-e75 or --t2-dir" 2 ;;
  3) fail "labels missing from configs/lineage_map.yaml (listed above): map them by the rule in the map's header, commit, sync to the DGX and re-run" 3 ;;
  4) fail "uncommitted changes to the procedure files (listed above): commit them first, or pass --allow-uncommitted" 4 ;;
  *) fail "preflight error (exit $rc)" "$rc" ;;
esac

step "2. intake QC (vec data qc) on the new files"
uv run vec data qc "$E85" "$E105" "$E75" --out "$OUT/qc_p3.json" --csv "$OUT/qc_p3.csv"

NEW=(--stage "E7.5=$E75" --stage "E8.5=$E85" --stage "E10.5=$E105")
step "3a. characterise"
uv run python experiments/coords/characterise.py --out "$OUT" "${NEW[@]}"
step "3b. z pitch (characterisation only; never used by a procedure)"
uv run python experiments/coords/zpitch_fit.py --out "$OUT" "${NEW[@]}"
step "3c. handedness: landmarks (E8.5, E10.5), mirror Dice, shape alignment against the neighbouring stages"
uv run python experiments/coords/handedness.py --out "$OUT" \
  --stage "E7.25=$E725" --stage "E7.5=$E75" --stage "E8.0=$E80" --stage "E8.25_late=$E825" --stage "E8.5=$E85" \
  --stage "E8.75=$E875" --stage "E9.5=$E95" --stage "E10.5=$E105" \
  --landmark-stages E8.5 E10.5 --mirror-stages E7.5 E8.5 E10.5 \
  --same-stages E7.25 E7.5 E8.0 E8.25_late E8.5 E8.75 E9.5 E10.5 \
  --pairs E7.25:E7.5 E7.5:E8.0 E8.25_late:E8.5 E8.5:E8.75 E9.5:E10.5

step "4a. heart size procedure -> E12.5"
uv run python experiments/coords/size_procedure.py --board heart --stage "$E825" 8.25 --stage "$E85" 8.5 \
  --stage "$E875" 8.75 --stage "$E95" 9.5 --stage "$E105" 10.5 --target 12.5 --timestamp "$TS" $REDACT \
  --out "$OUT/size_heart_E12.5.json"
step "4b. embryo size procedure -> E7.75"
uv run python experiments/coords/size_procedure.py --board embryo --stage "$E625" 6.75 --stage "$E725" 7.25 \
  --stage "$E75" 7.5 --stage "$E80" 8.0 --target 7.75 --timestamp "$TS" $REDACT \
  --out "$OUT/size_embryo_E7.75.json"

uv run python - "$OUT" > "$OUT/p3_summary.txt" <<'EOF'
import json
import sys
from pathlib import Path

out = Path(sys.argv[1])
for name in ("size_heart_E12.5.json", "size_embryo_E7.75.json"):
    r = json.loads((out / name).read_text())
    print(f"{name}: valid={r['valid']} timestamp={r['timestamp']} commit={r['provenance']['git_commit']} "
          f"uncommitted={r['provenance']['git_uncommitted']} pre-registered z mode={r['preregistered']['z_mode']}")
    for i in r["inputs"]:
        print(f"  input {i['name']} t={i['t']} records={i['n_input']} duplicates dropped={i['n_dropped_duplicates']}")
    for s in r["summary"]:
        print(f"  z mode {s['z_mode']}{' (pre-registered)' if s['preregistered'] else ''}: valid={s['valid']} "
              f"chosen={s.get('chosen')} core={s.get('pred_rms_core')} full={s.get('pred_rms_full')}")
for name in ("handedness_verdict.csv", "shape_alignment_summary.csv", "mirror_dice_summary.csv"):
    p = out / name
    if p.exists():
        print(f"--- {name}\n{p.read_text()}")
EOF
cat "$OUT/p3_summary.txt"

if [ "$REHEARSAL" = 1 ]; then
  step "rehearsal finished: outputs in $OUT (code checks only; nothing copied to results)"
  exit 0
fi
step "5. copy aggregates to $RESULTS"
mkdir -p "$RESULTS"
cp "$OUT"/*.csv "$OUT"/*.json "$OUT"/*.txt "$RESULTS"/
ls -la "$RESULTS"
echo
echo "Done. From the local machine: rsync -az -e ssh PrakashDGX_H2:v/virtual-embryo-challenge/$RESULTS/ $RESULTS/"
echo "Commit both size JSONs (with the git hash) before any test submission. The T2 build applies the heart"
echo "JSON through vec.t2_build.size_or_hold (approved rule 3: valid=false keeps E10.5's size); an invalid"
echo "embryo JSON is refused."
