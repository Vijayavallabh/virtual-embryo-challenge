#!/usr/bin/env bash
# The Task 1 outside-atlas tables on the DGX, strict form (candidate MX2s, ruling row R27). Every output goes under
# ~/v/vec-ext/t1 (outside every checkout); each step writes logs/<step>.log and logs/<step>.exit; logs/ALL.exit is
# written last, and also on a kill (99); a previous ALL.exit is moved to ALL.exit.prev first (this launcher deletes
# nothing). One process at a time, 8 BLAS threads, CPU only, TMPDIR under ~/v/scratch (the system disk was full before).
#
#   run_dgx.sh CHECKOUT [MODE] [ACCESSIONS]
#   CHECKOUT: a campaign checkout holding the permitted-only recomputed file of row N54 under data/derived (the
#        tables read the atlas in use through it) and the released T1 training stages under data/raw/T1. The DGX's
#        reference checkout ~/v/virtual-embryo-challenge is refused.
#   The steps, in order (README, "Steps"); the main session commits between them:
#   MODE plan (the default): fetch.py plan per accession (GEO family SOFT metadata only: plan.json, plan_table.tsv,
#        records.txt, records/<ACC>.txt and processing_text.txt under meta/ACC). Nothing else is downloaded and no
#        count is read. Two independent readers read every printed record line by line; for the named samples
#        `fetch.py pin-template ACC --dest ~/v/vec-ext/t1 --gsm G --gsm G ...` prints the pinned block and the
#        processing_cleared line (no hash typed by hand). The main session commits them, syncs and runs plan again; a
#        plan passes only when its pins match the records, cover every named sample, and its processing text is
#        cleared.
#   MODE derive: groups.py derive on the released E8.5 (a training stage; target-free): the extra competing classes,
#        written to ~/v/vec-ext/t1/derive/extra_competing.json. The main session commits it as
#        experiments/t1_ext/results/extra_competing.json, pins its sha256 in atlases.yaml and freezes the yaml.
#   MODE heads: fetch.py heads per accession: the barcode and feature files only, never a matrix. Needs PLAN_SHAS in
#        the environment, "ACC=<sha256 of meta/ACC/plan.json> ...", one per accession: a plan nobody read is not used.
#   MODE screens: screens.py per accession on those files (yes or no with counts; screens/ACC/screen.json). The main
#        session pastes each printed line into experiments/t1_ext/screens.yaml (kept or dropped) and commits it.
#   MODE fetch | load | tables | rest (fetch, load and tables in order): only after the screens record and the frozen
#        yaml are committed. fetch takes the matrices of kept members only (PLAN_SHAS as for heads); load deletes each
#        atlas's raw folder after its reduction (load.py --delete-raw); tables writes ~/v/vec-ext/t1/tables/$NAME
#        (NAME default consensus_v1) and reads each reduction only under the pins in force.
#   ACCESSIONS: a quoted, space-separated list (default: the two atlases in use, GSE244783 and GSE230531).
#
# Start it detached and wait with a background poller for ~/v/vec-ext/t1/logs/ALL.exit:
#   setsid nohup bash experiments/t1_ext/run_dgx.sh "$PWD" plan > ~/v/vec-ext/t1/launcher_plan.log 2>&1 < /dev/null &
set -u
CHECKOUT=${1:?usage: run_dgx.sh CHECKOUT [plan|derive|heads|screens|fetch|load|tables|rest] [ACCESSIONS]}
MODE=${2:-plan}
D=$HOME/v/vec-ext/t1
TMP=$HOME/v/scratch/t1ext-tmp
NAME=${NAME:-consensus_v1}
mkdir -p "$D/logs" "$TMP"
if [ -e "$D/logs/ALL.exit" ]; then mv -f "$D/logs/ALL.exit" "$D/logs/ALL.exit.prev"; fi
trap 'echo 99 > "$D/logs/ALL.exit"; exit 99' TERM INT HUP
case "$MODE" in
  plan|derive|heads|screens|fetch|load|tables|rest) ACCS=${3:-"GSE244783 GSE230531"} ;;
  *) echo "MODE must be plan, derive, heads, screens, fetch, load, tables or rest" > "$D/logs/start_$MODE.txt"
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

{ date -u; echo "mode $MODE; accessions $ACCS; name $NAME"; git -C "$CO" rev-parse HEAD
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
if [ "$MODE" = plan ]; then
  for A in $ACCS; do step "plan_$A" "${PY[@]}" experiments/t1_ext/fetch.py plan "$A" --dest "$D" || rc=2; done
  echo "$rc" > "$D/logs/ALL.exit"; exit 0
fi
if [ "$MODE" = derive ]; then
  step derive "${PY[@]}" experiments/t1_ext/groups.py derive --released-e85 "$CO/data/raw/T1/E8.5_RNA.h5ad" \
    --out "$D/derive/extra_competing.json" || rc=1
  echo "$rc" > "$D/logs/ALL.exit"; exit 0
fi
if [ "$MODE" = heads ] || [ "$MODE" = fetch ] || [ "$MODE" = rest ]; then
  SUB=fetch; [ "$MODE" = heads ] && SUB=heads
  for A in $ACCS; do
    SHA=$(plan_sha "$A")
    if [ -z "$SHA" ]; then echo "no plan sha256 given for $A" > "$D/logs/${SUB}_$A.log"; echo 2 > "$D/logs/${SUB}_$A.exit"
      rc=2; continue; fi
    step "${SUB}_$A" "${PY[@]}" experiments/t1_ext/fetch.py "$SUB" "$A" --dest "$D" --plan-sha256 "$SHA" || rc=2
  done
  [ "$MODE" = rest ] || { echo "$rc" > "$D/logs/ALL.exit"; exit 0; }
  [ "$rc" = 0 ] || { echo "$rc" > "$D/logs/ALL.exit"; exit 0; }
fi
if [ "$MODE" = screens ]; then
  for A in $ACCS; do step "screens_$A" "${PY[@]}" experiments/t1_ext/screens.py "$A" --dest "$D" || rc=2; done
  echo "$rc" > "$D/logs/ALL.exit"; exit 0
fi
if [ "$MODE" = load ] || [ "$MODE" = rest ]; then
  for A in $ACCS; do
    step "load_$A" "${PY[@]}" experiments/t1_ext/load.py "$A" --dest "$D" --delete-raw || rc=2
  done
  [ "$MODE" = load ] && { echo "$rc" > "$D/logs/ALL.exit"; exit 0; }
  [ "$rc" = 0 ] || { echo "$rc" > "$D/logs/ALL.exit"; exit 0; }
fi
# shellcheck disable=SC2086
step "tables_$NAME" "${PY[@]}" experiments/t1_ext/tables.py --dest "$D" --atlases $ACCS \
  --released-dir "$CO/data/raw/T1" --name "$NAME" || rc=2
echo "$rc" > "$D/logs/ALL.exit"
