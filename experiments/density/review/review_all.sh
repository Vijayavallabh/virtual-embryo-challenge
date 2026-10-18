#!/bin/bash
# Reviewer's launcher: one CPU process per (target, seed), 24 at a time, 2 BLAS threads each.
# Copied from the independent review of 2026-10-19; only the script path changed. Run from the repo root.
set -u
export PATH=$HOME/.local/bin:$PATH
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
mkdir -p runs/logs runs/density_review/raw
for s in $(seq 0 19); do for t in heart_E9.5 heart_E8.75 embryo_E8.0 embryo_E7.25; do echo "$t $s"; done; done | \
  xargs -P 24 -n 2 sh -c 'uv run python experiments/density/review/review_run.py --target "$0" --seed "$1" --out runs/density_review/raw > "runs/logs/density_review_$0_s$1.log" 2>&1'
echo "ALL DONE $(date)" >> runs/logs/density_review_all.log
