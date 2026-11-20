#!/usr/bin/env bash
# The T3 outside-knockout tables on the DGX. Every output goes under ~/v/vec-ext/t3 (outside every checkout); each
# step writes logs/<step>.log and logs/<step>.exit; logs/ALL.exit is written last, and also on a kill (99). One
# process at a time, 8 BLAS threads, CPU only.
#
#   run_dgx.sh CHECKOUT [MODE] [ACCESSIONS]
#   CHECKOUT: the campaign checkout ~/v/vec-campaign-c after the main session's usual sync of main into it. The DGX's
#        reference checkout ~/v/virtual-embryo-challenge is refused (nothing is fetched into it, nothing runs in it).
#   MODE plan (the default): only `fetch.py plan` per accession (GEO family SOFT metadata of every screened series
#        and of a named control series, the screen, the plan tables). Nothing else is downloaded and no count is
#        read. Two readers read every plan table and pin samples, series files and the control series' samples.
#   MODE all: first the PINNED gate (`fetch.py check-pins`): refused (ALL.exit 6) unless every accession is cleared,
#        has a non-empty pinned block, and datasets.yaml hashes to t3x_common.DATASETS_SHA256. Then, per accession, a
#        fresh plan that checks each pinned sample's live metadata against its pinned hash (replan_<ACC>), and only if
#        it passes: fetch (pinned files only), process; then the reference statistics first and combine last.
#   Stage 6: a dataset with a named control series (control_series in datasets.yaml) plans, screens, shows and pins
#        that control series with its own; its re-plan checks the control series' live text too. fetch also
#        downloads the pinned series files (pinned_series_files, a series' own archives or loose matrices) into
#        raw/<ACC>/<GSE>/, never unpacked; process reads them member by member in memory and splits them by their
#        pinned rule. The --max-gb 30 cap covers series files too: a dataset over it refuses (raise it only in a
#        reviewed commit).
#   Stage 7: a dataset with a raw_matrices declaration (GSE270969 only) offers its samples' own unfiltered matrices in
#        its plan tables; once two readers have pinned them (README check 12), fetch downloads the pinned raw files
#        like any other and process reads each through the one fixed cell call of t3x_common (README "Unfiltered
#        per-sample matrices"), streamed in two passes, never unpacked. A raw archive may be several GiB: the
#        --max-gb 30 cap covers it too (raise it only in a reviewed commit). Name GSE270969 explicitly, e.g.
#        run_dgx.sh "$PWD" plan "GSE270969", while it is not in the default list below.
#   Stage 9 (bulk whole-embryo units, README "Bulk whole-embryo units"): MODE header runs only the header step
#        (`fetch.py header`) of a bulk count table after its plan passed: the first line of the table, streamed and
#        closed at the first newline, screened and tied; nothing else is downloaded. Name each item as
#        ACC:SERIES:FILE for a GEO series table, or ACC alone for the ArrayExpress resource (E-MTAB-7322), e.g.
#        run_dgx.sh "$PWD" header "GSE214161:GSE214161:GSE214161_counts.txt.gz E-MTAB-7322". Bulk datasets go
#        through the same plan, pinned gate, re-plan, fetch, process and combine as any other in MODE all; an
#        E-MTAB accession is accepted there besides GSE ones. A bulk dataset with id_type ensembl is refused at the
#        gate while t3x_common.PANEL_IDS_SHA256 is None.
#   ACCESSIONS: a quoted, space-separated subset (default: four cleared datasets; GSE205117 left the list in stage 8:
#        it is held now and on the team's held list, t3x_common.HELD_EXTRA).
#   WT (environment, optional): the WT E8.75 file, read only; default CHECKOUT/data/raw/T2/E8.75.h5ad.
#
# Start it detached and wait with a background poller for ~/v/vec-ext/t3/logs/ALL.exit:
#   setsid nohup bash experiments/t3_ext/run_dgx.sh "$PWD" plan > ~/v/vec-ext/t3/launcher_plan.log 2>&1 < /dev/null &
set -u
CHECKOUT=${1:?usage: run_dgx.sh CHECKOUT [plan|all] [ACCESSIONS]}
MODE=${2:-plan}
ACCS=${3:-"GSE168768 GSE186934 GSE204906 GSE137337"}
D=$HOME/v/vec-ext/t3
TMP=$HOME/v/scratch/t3ext-tmp
mkdir -p "$D/logs" "$TMP"
rm -f "$D/logs/ALL.exit"
trap 'echo 99 > "$D/logs/ALL.exit"; exit 99' TERM INT HUP
case "$MODE" in plan|all|header) ;; *) echo "MODE must be plan, header or all" > "$D/logs/start_$MODE.txt"; echo 3 > "$D/logs/ALL.exit"; exit 3;; esac
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

if [ "$MODE" = header ]; then
  # stage 9: the header step only (one first line per bulk count table); no count is read
  rc=0
  for ITEM in $ACCS; do
    IFS=: read -r A S F <<< "$ITEM"
    if [[ "$A" =~ ^E-MTAB-[0-9]+$ ]] && [ -z "$S" ]; then
      step "header_$A" "${PY[@]}" experiments/t3_ext/fetch.py header "$A" --dest "$D" || rc=2
    elif [[ "$A" =~ ^GSE[0-9]+$ ]] && [[ "$S" =~ ^GSE[0-9]+$ ]] && [[ "$F" =~ ^[A-Za-z0-9._+-]+$ ]]; then
      step "header_$A" "${PY[@]}" experiments/t3_ext/fetch.py header "$A" --dest "$D" --series "$S" --name "$F" || rc=2
    else
      echo "refused: $ITEM is not ACC:SERIES:FILE or an E-MTAB accession" >> "$D/logs/header_bad.log"; rc=2
    fi
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
  # no stale output: an accession's processed outputs go before its re-plan, so a refused re-plan, fetch or process
  # leaves nothing for combine to read (stage 5)
  if ! [[ "$A" =~ ^GSE[0-9]+$ || "$A" =~ ^E-MTAB-[0-9]+$ ]]; then echo "refused: not a GEO series or E-MTAB accession" >> "$D/logs/replan_bad.log"; continue; fi
  rm -rf "${D:?}/processed/$A"
  step "replan_$A" "${PY[@]}" experiments/t3_ext/fetch.py plan "$A" --dest "$D" || {
    echo "skipped: the plan with pins did not pass for $A" > "$D/logs/fetch_$A.log"; echo 2 > "$D/logs/fetch_$A.exit"; continue; }
  step "fetch_$A" "${PY[@]}" experiments/t3_ext/fetch.py fetch "$A" --dest "$D" --max-gb 30 || continue
  step "process_$A" "${PY[@]}" experiments/t3_ext/process.py run "$A" --dest "$D" \
    --wt "$WT" --refstats "$REF" --out "$D/processed" --threads 8 || true
done
step combine "${PY[@]}" experiments/t3_ext/combine.py --processed "$D/processed" \
  --refstats "$REF" --out "$D/tables" --dest "$D"
echo $? > "$D/logs/ALL.exit"
