#!/bin/bash
# Intake QC of every released file with `vec data qc` (vec.spatial.qc_stage). Run on the DGX from the repo root:
#   bash experiments/qc/run_qc.sh [tag]
# Writes aggregates only (no per-cell rows, no expression values), safe to copy back:
#   experiments/qc/results/qc_<tag>.json (full report per file) and qc_<tag>.csv (one summary row per file).
# T1 has no coordinates, so it gets the row-identity part only.
set -euo pipefail
tag=${1:-released_$(date +%F)}
out=experiments/qc/results
mkdir -p "$out"
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-8} OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-8} \
    MKL_NUM_THREADS=${MKL_NUM_THREADS:-8} COLUMNS=${COLUMNS:-260}
shopt -s nullglob
files=(data/raw/T2/*.h5ad data/raw/T3/*.h5ad data/raw/T1/*.h5ad)
uv run vec data qc "${files[@]}" --out "$out/qc_${tag}.json" --csv "$out/qc_${tag}.csv"
