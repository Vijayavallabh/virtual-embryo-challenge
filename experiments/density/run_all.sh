#!/bin/bash
# Launch the density experiment on the DGX: one CPU process per (target, seed), JOBS at a time, THREADS BLAS threads
# each (keep JOBS x THREADS at 16 or less). Run from ~/v/virtual-embryo-challenge, for example:
#   tmux new -d -s vec-density-v2 'bash experiments/density/run_all.sh'
#   tmux new -d -s vec-density-release 'TARGETS="heart_E10.5 embryo_E7.5" SEEDS="0 1 2 3 4 5 6 7 8 9" \
#     RAW=runs/density/raw_release TAG=density_release bash experiments/density/run_all.sh'
# Environment overrides [defaults]:
#   TARGETS   [heart_E9.5 heart_E8.75 embryo_E8.0 heart_E10.5 embryo_E7.5]; release targets without files are skipped
#   SEEDS     [0 1 2 3 4]
#   RAW       [runs/density/raw_v2]  output directory; finished conditions are skipped on rerun
#   PATTERNS  [all]  e.g. "block head refblock refhead"
#   TAG       [density]  logs go to runs/logs/${TAG}_<target>_s<seed>.log
#   JOBS      [15]   THREADS [1]
set -u
export PATH=$HOME/.local/bin:$PATH
export TARGETS=${TARGETS:-"heart_E9.5 heart_E8.75 embryo_E8.0 heart_E10.5 embryo_E7.5"}
export SEEDS=${SEEDS:-"0 1 2 3 4"}
export RAW=${RAW:-runs/density/raw_v2}
export PATTERNS=${PATTERNS:-}
export TAG=${TAG:-density}
JOBS=${JOBS:-15}
THREADS=${THREADS:-1}
export OMP_NUM_THREADS=$THREADS OPENBLAS_NUM_THREADS=$THREADS MKL_NUM_THREADS=$THREADS
mkdir -p runs/logs "$RAW"
echo "START $(date) targets=[$TARGETS] seeds=[$SEEDS] raw=$RAW patterns=[${PATTERNS:-all}]" >> "runs/logs/${TAG}_all.log"
for t in $TARGETS; do for s in $SEEDS; do echo "$t $s"; done; done | \
  xargs -P "$JOBS" -n 2 sh -c 'uv run python experiments/density/run_density.py --target "$0" --seeds "$1" \
    --out "$RAW" ${PATTERNS:+--patterns $PATTERNS} > "runs/logs/${TAG}_$0_s$1.log" 2>&1'
echo "ALL DONE $(date)" >> "runs/logs/${TAG}_all.log"
