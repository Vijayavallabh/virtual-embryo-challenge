#!/usr/bin/env bash
# The T3 outside-knockout tables on the DGX. Every output goes under ~/v/vec-ext/t3 (outside every checkout); each
# step writes logs/<step>.log and logs/<step>.exit; logs/ALL.exit is written last, and also on a kill (99). One
# process at a time, 8 BLAS threads, CPU only.
#
#   run_dgx.sh CHECKOUT [MODE] [ACCESSIONS]
#   CHECKOUT: the campaign checkout ~/v/vec-campaign-c after the main session's usual sync of main into it. The DGX's
#        reference checkout ~/v/virtual-embryo-challenge is refused (nothing is fetched into it, nothing runs in it).
#   MODE plan (the default): only `fetch.py plan` per accession (GEO family SOFT metadata, the screen, the plan
#        tables). Nothing else is downloaded and no count is read. Two readers read every plan table and pin samples.
#   MODE all: first the PINNED gate (`fetch.py check-pins`): refused (ALL.exit 6) unless every accession is cleared,
#        has a non-empty pinned block, and datasets.yaml hashes to t3x_common.DATASETS_SHA256. Then, per accession, a
#        fresh plan that checks each pinned sample's live metadata against its pinned hash (replan_<ACC>), and only if
#        it passes: fetch (pinned files only), process; then the reference statistics first and combine last.
#   ACCESSIONS: a quoted, space-separated subset (default: the five cleared datasets).
#   WT (environment, optional): the WT E8.75 file, read only; default CHECKOUT/data/raw/T2/E8.75.h5ad.
#
# Start it detached and wait with a background poller for ~/v/vec-ext/t3/logs/ALL.exit:
#   setsid nohup bash experiments/t3_ext/run_dgx.sh "$PWD" plan > ~/v/vec-ext/t3/launcher_plan.log 2>&1 < /dev/null &
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
# never the reference checkout
REFCO=$(cd "$HOME/v/virtual-embryo-challenge" 2>/dev/null && pwd -P || echo /nonexistent)
CO=$(cd "$CHECKOUT" 2>/dev/null && pwd -P) || { echo 3 > "$D/logs/ALL.exit"; exit 3; }
case "$CO/" in "$REFCO"/*) echo "refused: $CO is the reference checkout" > "$D/logs/start_$MODE.txt"; echo 3 > "$D/logs/ALL.exit"; exit 3;; esac
cd "$CO" || { echo 3 > "$D/logs/ALL.exit"; exit 3; }
export PATH=$HOME/.local/bin:$PATH
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 NUMPY_MADVISE_HUGEPAGE=0 TMPDIR=$TMP
WT=${WT:-$CO/data/raw/T2/E8.75.h5ad}
REF=$D/reference/refstats.npz
PY=(uv run --frozen --no-sync python -I -B)

step() {
  local name=$1; shift
  "$@" > "$D/logs/$name.log" 2>&1
  local rc=$?
  echo "$rc" > "$D/logs/$name.exit"
  return "$rc"
}

{ date -u; echo "mode $MODE; accessions $ACCS; WT $WT"; git -C "$CO" rev-parse HEAD
  git -C "$CO" status --porcelain --untracked-files=no; df -h /home
  "${PY[@]}" -c "import vec, sys; print('vec from', vec.__file__)"; } > "$D/logs/start_$MODE.txt" 2>&1
# the package must be imported from this checkout, never from another checkout's environment
if ! "${PY[@]}" -c "import vec, pathlib, sys; sys.exit(0 if pathlib.Path(vec.__file__).resolve().is_relative_to(pathlib.Path('$CO').resolve()) else 1)"; then
  echo "vec is not imported from $CO" >> "$D/logs/start_$MODE.txt"; echo 4 > "$D/logs/ALL.exit"; exit 4
fi

if [ "$MODE" = plan ]; then
  rc=0
  for A in $ACCS; do
    step "plan_$A" "${PY[@]}" experiments/t3_ext/fetch.py plan "$A" --dest "$D" || rc=2
  done
  echo "$rc" > "$D/logs/ALL.exit"
  exit 0
fi

# the pinned gate: every accession pinned, and datasets.yaml exactly the committed, hash-pinned file
# shellcheck disable=SC2086
step pinned_gate "${PY[@]}" experiments/t3_ext/fetch.py check-pins $ACCS || { echo 6 > "$D/logs/ALL.exit"; exit 6; }
[ -f "$WT" ] || { echo "WT file $WT missing" >> "$D/logs/start_$MODE.txt"; echo 5 > "$D/logs/ALL.exit"; exit 5; }
step reference "${PY[@]}" experiments/t3_ext/process.py reference --wt "$WT" --out "$REF" \
  || { echo 1 > "$D/logs/ALL.exit"; exit 1; }
for A in $ACCS; do
  step "replan_$A" "${PY[@]}" experiments/t3_ext/fetch.py plan "$A" --dest "$D" || {
    echo "skipped: the plan with pins did not pass for $A" > "$D/logs/fetch_$A.log"; echo 2 > "$D/logs/fetch_$A.exit"; continue; }
  step "fetch_$A" "${PY[@]}" experiments/t3_ext/fetch.py fetch "$A" --dest "$D" --max-gb 30 || continue
  step "process_$A" "${PY[@]}" experiments/t3_ext/process.py run "$A" --dest "$D" \
    --wt "$WT" --refstats "$REF" --out "$D/processed" --threads 8 || true
done
step combine "${PY[@]}" experiments/t3_ext/combine.py --processed "$D/processed" \
  --refstats "$REF" --out "$D/tables"
echo $? > "$D/logs/ALL.exit"
