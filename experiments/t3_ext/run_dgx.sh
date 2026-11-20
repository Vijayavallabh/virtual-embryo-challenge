#!/usr/bin/env bash
# Stage 1 of the T3 outside-knockout tables on the DGX: reference statistics, then per cleared dataset plan -> fetch ->
# process, then combine. Every output goes under ~/v/vec-ext/t3 (outside every checkout); each step writes
# logs/<step>.log and logs/<step>.exit; logs/ALL.exit is written last. One process at a time, 8 BLAS threads.
#
# Start it detached from a copied script file, from a campaign checkout at this branch's commit (docs/COMPUTE.md):
#   cd ~/v/vec-campaign-d
#   setsid nohup bash experiments/t3_ext/run_dgx.sh "$PWD" > ~/v/vec-ext/t3/launcher.log 2>&1 < /dev/null &
# and wait with a background poller for ~/v/vec-ext/t3/logs/ALL.exit.
set -u
CHECKOUT=${1:?usage: run_dgx.sh CHECKOUT}
D=$HOME/v/vec-ext/t3
TMP=$HOME/v/scratch/t3ext-tmp
mkdir -p "$D/logs" "$TMP"
cd "$CHECKOUT" || exit 3
export PATH=$HOME/.local/bin:$PATH
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 NUMPY_MADVISE_HUGEPAGE=0 TMPDIR=$TMP
WT=$CHECKOUT/data/raw/T2/E8.75.h5ad
REF=$D/reference/refstats.npz
ACCS="GSE168768 GSE186934 GSE205117 GSE204906 GSE137337"

step() {
  local name=$1; shift
  "$@" > "$D/logs/$name.log" 2>&1
  local rc=$?
  echo "$rc" > "$D/logs/$name.exit"
  return "$rc"
}

{ date -u; git -C "$CHECKOUT" rev-parse HEAD; git -C "$CHECKOUT" status --porcelain --untracked-files=no; df -h /home; } \
  > "$D/logs/start.txt" 2>&1
step reference uv run --frozen --no-sync python -I -B experiments/t3_ext/process.py reference --wt "$WT" --out "$REF" \
  || { echo 1 > "$D/logs/ALL.exit"; exit 1; }
for A in $ACCS; do
  step "plan_$A" uv run --frozen --no-sync python -I -B experiments/t3_ext/fetch.py plan "$A" --dest "$D" || continue
  step "fetch_$A" uv run --frozen --no-sync python -I -B experiments/t3_ext/fetch.py fetch "$A" --dest "$D" \
    --max-gb 30 || continue
  step "process_$A" uv run --frozen --no-sync python -I -B experiments/t3_ext/process.py run "$A" --dest "$D" \
    --wt "$WT" --refstats "$REF" --out "$D/processed" --threads 8 || true
done
step combine uv run --frozen --no-sync python -I -B experiments/t3_ext/combine.py --processed "$D/processed" \
  --refstats "$REF" --out "$D/tables"
echo $? > "$D/logs/ALL.exit"
