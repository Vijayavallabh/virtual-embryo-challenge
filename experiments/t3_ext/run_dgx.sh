#!/usr/bin/env bash
# Stage 1 of the T3 outside-knockout tables on the DGX. Every output goes under ~/v/vec-ext/t3 (outside every
# checkout); each step writes logs/<step>.log and logs/<step>.exit; logs/ALL.exit is written last, and also on a
# kill (99). One process at a time, 8 BLAS threads, CPU only.
#
#   run_dgx.sh CHECKOUT [MODE] [ACCESSIONS]
#   MODE plan (the default): only `fetch.py plan` per accession (GEO family SOFT metadata, the screen, the sample
#        table). Nothing else is downloaded and no count is read. A person reads every plan before MODE all.
#   MODE all: reference statistics, then fetch and process for each accession whose plan_<ACC>.exit is 0 (plans
#        are NOT rerun, so what is fetched is what was read), then combine.
#   ACCESSIONS: a quoted, space-separated subset (default: the five cleared datasets).
#   WT (environment, optional): the WT E8.75 file; default CHECKOUT/data/raw/T2/E8.75.h5ad.
#
# Start it detached from a scratch clone at this branch's commit (README, 'Running on the DGX'):
#   setsid nohup bash experiments/t3_ext/run_dgx.sh "$PWD" plan > ~/v/vec-ext/t3/launcher_plan.log 2>&1 < /dev/null &
# and wait with a background poller for ~/v/vec-ext/t3/logs/ALL.exit.
set -u
CHECKOUT=${1:?usage: run_dgx.sh CHECKOUT [plan|all] [ACCESSIONS]}
MODE=${2:-plan}
ACCS=${3:-"GSE168768 GSE186934 GSE205117 GSE204906 GSE137337"}
D=$HOME/v/vec-ext/t3
TMP=$HOME/v/scratch/t3ext-tmp
mkdir -p "$D/logs" "$TMP"
rm -f "$D/logs/ALL.exit"
trap 'echo 99 > "$D/logs/ALL.exit"; exit 99' TERM INT HUP
case "$MODE" in plan|all) ;; *) echo "MODE must be plan or all" > "$D/logs/start_$MODE.txt"; echo 3 > "$D/logs/ALL.exit"; exit 3;; esac
cd "$CHECKOUT" || { echo 3 > "$D/logs/ALL.exit"; exit 3; }
export PATH=$HOME/.local/bin:$PATH
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 NUMPY_MADVISE_HUGEPAGE=0 TMPDIR=$TMP
WT=${WT:-$CHECKOUT/data/raw/T2/E8.75.h5ad}
REF=$D/reference/refstats.npz
PY=(uv run --frozen --no-sync python -I -B)

step() {
  local name=$1; shift
  "$@" > "$D/logs/$name.log" 2>&1
  local rc=$?
  echo "$rc" > "$D/logs/$name.exit"
  return "$rc"
}

{ date -u; echo "mode $MODE; accessions $ACCS; WT $WT"; git -C "$CHECKOUT" rev-parse HEAD
  git -C "$CHECKOUT" status --porcelain --untracked-files=no; df -h /home
  "${PY[@]}" -c "import vec, sys; print('vec from', vec.__file__)"; } > "$D/logs/start_$MODE.txt" 2>&1
# the package must be imported from this checkout, never from another checkout's environment
if ! "${PY[@]}" -c "import vec, pathlib, sys; sys.exit(0 if pathlib.Path(vec.__file__).resolve().is_relative_to(pathlib.Path('$CHECKOUT').resolve()) else 1)"; then
  echo "vec is not imported from $CHECKOUT" >> "$D/logs/start_$MODE.txt"; echo 4 > "$D/logs/ALL.exit"; exit 4
fi

if [ "$MODE" = plan ]; then
  rc=0
  for A in $ACCS; do
    step "plan_$A" "${PY[@]}" experiments/t3_ext/fetch.py plan "$A" --dest "$D" || rc=2
  done
  echo "$rc" > "$D/logs/ALL.exit"
  exit 0
fi

[ -f "$WT" ] || { echo "WT file $WT missing" >> "$D/logs/start_$MODE.txt"; echo 5 > "$D/logs/ALL.exit"; exit 5; }
step reference "${PY[@]}" experiments/t3_ext/process.py reference --wt "$WT" --out "$REF" \
  || { echo 1 > "$D/logs/ALL.exit"; exit 1; }
for A in $ACCS; do
  if [ "$(cat "$D/logs/plan_$A.exit" 2>/dev/null)" != 0 ]; then
    echo "skipped: no passing plan read for $A" > "$D/logs/fetch_$A.log"; echo 2 > "$D/logs/fetch_$A.exit"; continue
  fi
  step "fetch_$A" "${PY[@]}" experiments/t3_ext/fetch.py fetch "$A" --dest "$D" --max-gb 30 || continue
  step "process_$A" "${PY[@]}" experiments/t3_ext/process.py run "$A" --dest "$D" \
    --wt "$WT" --refstats "$REF" --out "$D/processed" --threads 8 || true
done
step combine "${PY[@]}" experiments/t3_ext/combine.py --processed "$D/processed" \
  --refstats "$REF" --out "$D/tables"
echo $? > "$D/logs/ALL.exit"
