#!/bin/bash
# Release-day run (2026-10-20) of the pre-registered T3 decomposition on the two Gata4 replicates against WT E8.75.
# The rules in experiments/t3_decomp/rules.yaml decide; nothing is chosen by eye. Run on the DGX from the repository
# root, in tmux (about 20-30 minutes; the log goes to runs/logs/):
#
#   tmux new -d -s vec-t3-gata4 'bash experiments/t3_decomp/run_gata4.sh REP1.h5ad REP2.h5ad data/raw/T2/E8.75.h5ad'
#
# Arguments: the two Gata4 replicate files and WT E8.75, in that order. Anything after them goes to
# `t3_decompose.py gata4` (e.g. --e85 PATH if WT E8.5 is not at data/raw/T2/E8.5.h5ad, --out DIR). If REP1, REP2 or
# WT E8.75 is missing, the script says so and exits with status 0 without running anything (before the release).
#
# Afterwards, on the local machine, bring back the aggregates (CSV/JSON/TXT; the per-gene factors/ stay on the DGX)
# and commit the decision file before anything else:
#   rsync -az --include='*.csv' --include='*.json' --include='*.txt' --exclude='*' \
#     PrakashDGX_H2:v/virtual-embryo-challenge/<the --out directory printed at the end>/ experiments/t3_decomp/results/gata4/
set -euo pipefail
cd "$(dirname "$0")/../.."

if [ $# -lt 3 ]; then
  echo "usage: run_gata4.sh REP1.h5ad REP2.h5ad WT_E8.75.h5ad [options for t3_decompose.py gata4]" >&2
  exit 2
fi
for f in "$1" "$2" "$3"; do
  if [ ! -f "$f" ]; then
    echo "run_gata4.sh: $f not found; nothing was run (expected before the 2026-10-20 release)" >&2
    exit 0
  fi
done

export PATH="$HOME/.local/bin:$PATH"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-16}" OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-16}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-16}"
export POT_BACKEND_DISABLE_JAX=1 XLA_PYTHON_CLIENT_PREALLOCATE=false
STAMP=$(date -u +%Y%m%dT%H%M%SZ)
mkdir -p runs/logs
uv run python experiments/t3_decomp/t3_decompose.py gata4 "$1" "$2" "$3" --out "runs/t3_decomp/gata4_$STAMP" \
  --workers 8 --threads 2 "${@:4}" 2>&1 | tee "runs/logs/t3_decomp_gata4_$STAMP.log"
