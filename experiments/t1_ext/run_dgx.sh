#!/usr/bin/env bash
# Row T72 (candidate MX2s), the steps BEFORE the row's commit, on the DGX: the plan, the screens' fetch of barcodes and
# features only, the screens' facts, and the derived classes. Every output of plan, fetch-meta and screens goes under
# ~/v/scratch/t72_screen (outside every checkout, never the row's data folders); each step writes
# logs/<step>.log and logs/<step>.exit there; logs/ALL.exit is written last, and also on a kill (99); a previous
# ALL.exit is moved to ALL.exit.prev first (this launcher deletes nothing). 8 BLAS threads, CPU only, TMPDIR under
# ~/v/scratch.
#
#   run_dgx.sh CHECKOUT MODE [ACCESSIONS]
#   CHECKOUT: a campaign checkout (its own environment; the DGX reference checkout ~/v/virtual-embryo-challenge is
#        refused); derive also needs the released T1 E8.5 file under its data/raw/T1.
#   MODE plan: fetch.py plan per accession (the series family SOFT only). Prints per accession the plan summary and
#        one line per protocol or processing field of every sample record of the series (field_sha256 GSE GSM field
#        sha256). Exit 2 per accession is expected before the pins are committed.
#   MODE fetch-meta: fetch.py fetch-meta per accession: the six named samples barcodes and features files only, then
#        screens.py on them. Needs PLAN_SHAS="ACC=<sha256 of meta/ACC/plan.json> ..." in the environment.
#   MODE screens: screens.py per accession again on the fetched files (re-hashed): the facts of checks (a) to (d) and
#        the yaml screen block. Exit 2 when a check says no (the readers decide; nothing is dropped by code here).
#   MODE derive: groups.py derive on the released E8.5 (a training stage; target-free), written to
#        experiments/t1_ext/results/extra_competing.json in CHECKOUT; the main session commits it with the frozen yaml.
#   The row's own data steps (fetch, load, tables) are NOT run here: experiments/t1_m1/mx2s.sh runs them at the row's
#   commit (README, "The row's steps").
#   ACCESSIONS: a quoted, space-separated list (default: GSE244783 GSE230531).
#
# Start it detached and wait with a background poller for ~/v/scratch/t72_screen/logs/ALL.exit:
#   setsid nohup bash experiments/t1_ext/run_dgx.sh "$PWD" plan > ~/v/scratch/t72_screen/launcher_plan.log 2>&1 \
#     < /dev/null &
set -u
CHECKOUT=${1:?usage run_dgx.sh CHECKOUT plan or fetch-meta or screens or derive [ACCESSIONS]}
MODE=${2:?usage run_dgx.sh CHECKOUT plan or fetch-meta or screens or derive [ACCESSIONS]}
D=$HOME/v/scratch/t72_screen
TMP=$HOME/v/scratch/t72-tmp
mkdir -p "$D/logs" "$TMP"
if [ -e "$D/logs/ALL.exit" ]; then mv -f "$D/logs/ALL.exit" "$D/logs/ALL.exit.prev"; fi
trap 'echo 99 > "$D/logs/ALL.exit"; exit 99' TERM INT HUP
ACCS=${3:-"GSE244783 GSE230531"}
case "$MODE" in
  plan|fetch-meta|screens|derive) ;;
  *) echo "MODE must be plan, fetch-meta, screens or derive (the row steps run from experiments/t1_m1/mx2s.sh)" \
       > "$D/logs/start_bad_mode.txt"
     echo 3 > "$D/logs/ALL.exit"; exit 3 ;;
esac
REFCO=$(cd "$HOME/v/virtual-embryo-challenge" 2>/dev/null && pwd -P || echo /nonexistent)
CO=$(cd "$CHECKOUT" 2>/dev/null && pwd -P) || { echo 3 > "$D/logs/ALL.exit"; exit 3; }
case "$CO/" in "$REFCO"/*) echo "refused: $CO is the reference checkout" > "$D/logs/start_$MODE.txt"
  echo 3 > "$D/logs/ALL.exit"; exit 3;; esac
cd "$CO" || { echo 3 > "$D/logs/ALL.exit"; exit 3; }
export PATH=$HOME/.local/bin:$PATH
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 NUMPY_MADVISE_HUGEPAGE=0 TMPDIR=$TMP
PY=(uv run --frozen --no-sync python -I -B)

step() {
  local name=$1; shift
  "$@" > "$D/logs/$name.log" 2>&1
  local rc=$?
  echo "$rc" > "$D/logs/$name.exit"
  return "$rc"
}

plan_sha() {  # the sha256 given for accession $1 in PLAN_SHAS, or nothing
  local kv out=""
  for kv in ${PLAN_SHAS:-}; do [ "${kv%%=*}" = "$1" ] && out=${kv#*=}; done
  printf '%s' "$out"
}

{ date -u; echo "mode $MODE; accessions $ACCS"; git -C "$CO" rev-parse HEAD
  git -C "$CO" status --porcelain --untracked-files=no; df -h "$HOME" "$TMP"
  "${PY[@]}" -c "import vec; print('vec from', vec.__file__)"; } > "$D/logs/start_$MODE.txt" 2>&1
if ! "${PY[@]}" -c "import vec, pathlib, sys; sys.exit(0 if pathlib.Path(vec.__file__).resolve().is_relative_to(pathlib.Path('$CO').resolve()) else 1)"; then
  echo "vec is not imported from $CO" >> "$D/logs/start_$MODE.txt"; echo 4 > "$D/logs/ALL.exit"; exit 4
fi
for A in $ACCS; do
  [[ "$A" =~ ^GSE[0-9]+$ ]] || { echo "refused: $A is not a GEO series accession" >> "$D/logs/start_$MODE.txt"
                                 echo 3 > "$D/logs/ALL.exit"; exit 3; }
done

rc=0
case "$MODE" in
  plan)
    for A in $ACCS; do step "plan_$A" "${PY[@]}" experiments/t1_ext/fetch.py plan "$A" --dest "$D" || rc=2; done ;;
  fetch-meta)
    for A in $ACCS; do
      SHA=$(plan_sha "$A")
      if [ -z "$SHA" ]; then echo "no plan sha256 given for $A" > "$D/logs/fetch-meta_$A.log"
        echo 2 > "$D/logs/fetch-meta_$A.exit"; rc=2; continue; fi
      step "fetch-meta_$A" "${PY[@]}" experiments/t1_ext/fetch.py fetch-meta "$A" --dest "$D" --plan-sha256 "$SHA" \
        || rc=2
    done ;;
  screens)
    for A in $ACCS; do step "screens_$A" "${PY[@]}" experiments/t1_ext/screens.py "$A" --dest "$D" || rc=2; done ;;
  derive)
    step derive "${PY[@]}" experiments/t1_ext/groups.py derive --released-e85 "$CO/data/raw/T1/E8.5_RNA.h5ad" \
      --out "$CO/experiments/t1_ext/results/extra_competing.json" || rc=1 ;;
esac
echo "$rc" > "$D/logs/ALL.exit"
exit 0
