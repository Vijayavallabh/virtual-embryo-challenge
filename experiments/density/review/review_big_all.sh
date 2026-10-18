#!/bin/bash
# Reviewer's large-truth launcher (score_board frac=0.3): one CPU process per (target, seed), 20 at a time.
# Copied from the independent review of 2026-10-19; only the script path changed. Run from the repo root.
set -u
export PATH=$HOME/.local/bin:$PATH
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
mkdir -p runs/logs runs/density_review/raw_big
for s in $(seq 0 9); do for t in heart_E9.5 heart_E8.75; do echo "$t $s"; done; done | \
  xargs -P 20 -n 2 sh -c 'uv run python experiments/density/review/review_big.py --target "$0" --seed "$1" --frac 0.3 --out runs/density_review/raw_big > "runs/logs/density_review_big_$0_s$1.log" 2>&1'
echo "ALL DONE $(date)" >> runs/logs/density_review_big_all.log
