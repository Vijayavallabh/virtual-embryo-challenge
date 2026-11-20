#!/usr/bin/env bash
# The T3 file of the one form row K24's committed reading names (experiments/t3_decomp/build_k24_form.py), built on a
# DGX campaign checkout (~/v/vec-campaign-c, never the reference checkout). Run from the checkout's root, detached
# (setsid nohup ... &), and wait on $S/all.exit with a background poller; never a foreground sleep.
#
#   bash experiments/t3_decomp/build_k24_form.sh HEAD_SHA BASE KIND
#
# HEAD_SHA: the full commit the checkout must be at (it holds the committed K24 reading, TABLES_K24.sha256 and the
# filled pins). BASE: s1, s1x1.5, s1x2.0 or s2, the base B of decision log 2026-11-20 (q), fixed there; nothing here
# chooses it. KIND: sc3, bulk or joint, the one arm K24's reading names; nothing here chooses it, and the builder
# refuses a kind the reading does not name. The tables folder is fixed below (~/v/vec-ext/t3/tables) and passed only
# through VEC_T3_EXT_TABLES. The scratch folder S (~/v/scratch/t3k24_<first 7 of HEAD_SHA>_<BASE>) is made first, so
# every finish below writes $S/all.exit. Steps, in order, stopping at the first failure (build_outside_form.sh's):
#   0. BASE one of the four and KIND one of the three (exit 2);
#   1. HEAD and a clean tree (exit 97, 98), the tables folder a real directory (96), and no output folder yet (95);
#   2. the targeted tests (synthetic data only), 8 BLAS threads (exit 10);
#   3. identity-s1 once through build_outside_form.py (unchanged): its sha256 must be the reference barred s1's (exit
#      20 build, 21 check, 22 hash);
#   4. identity of BASE, when BASE is not s1: its sha256 must be the base file's reference hash (exit 30, 31, 32);
#   5. the form k24KIND-BASE twice, into folders a and b, each format-checked (`vec check --board T3:bcatenin`), each
#      stamp's commit checked against HEAD_SHA and its clean field against True, the pair compared (exit 40 build a,
#      41 build b, 42 check, 44 a stamp not of HEAD_SHA or not clean, 43 the pair differs, 46 the builder refused the
#      run (its exit 2: pins, listing, folder, bytes, records or the reading)). A form the builder refuses on its own
#      (its exit 4) is not built and is listed in $S/refused: per the row, nothing replaces it (all.exit 45).
# all.exit: 0 every step passed; 45 the form was refused ($S/refused); 2, 10, 20-22, 30-32, 40-44, 46, 95-98 as above;
# 99 the launcher stopped unexpectedly. 16 BLAS threads for builds; NUMPY_MADVISE_HUGEPAGE=0. Build and check output
# goes to /dev/null: a test build is read by its exit code, its sha256 and its stamp's commit and clean fields only
# (AGENTS.md), never its log, size or cell count. Written: $S/<step>.exit, $S/all.exit, $S/head, $S/SHA256SUMS and
# $S/stamps.
set -u
HEAD_SHA=${1:?usage: build_k24_form.sh HEAD_SHA BASE KIND}
BASE=${2:?usage: build_k24_form.sh HEAD_SHA BASE KIND}
KIND=${3:?usage: build_k24_form.sh HEAD_SHA BASE KIND}
S=~/v/scratch/t3k24_${HEAD_SHA:0:7}_${BASE//[^A-Za-z0-9.]/_}
mkdir -p "$S/tmp" || exit 99
finish() { trap - ERR TERM INT HUP EXIT; echo "$1" > "$S/all.exit"; exit "$1"; }
abort() { trap - ERR TERM INT HUP EXIT; echo 99 > "$S/all.exit"; exit 99; }
trap abort ERR TERM INT HUP EXIT

case "$BASE" in
  s1|s1x1.5|s1x2.0|s2) ;;
  *) echo "unknown base $BASE (s1, s1x1.5, s1x2.0 or s2)" >&2; finish 2 ;;
esac
case "$KIND" in
  sc3|bulk|joint) ;;
  *) echo "unknown kind $KIND (sc3, bulk or joint)" >&2; finish 2 ;;
esac
FORM="k24$KIND-$BASE"
declare -A REF=(
  [s1]=3712c3f9988d0f8aa4198b22c37b866a84ef8757fd500a9163944f0663a520a2
  [s1x1.5]=c060fcaa6c7fdc6d6b1a58c8f850a1895a8295abbc50f297de71a32c97f9c420
  [s1x2.0]=7393810fbe3be2cde67cf23ed54f45d22da738d7d3c9e124780a8f6281075661
  [s2]=1b69521e3c28bc84185410e838efc6fcbcff898434e73f9a7547208a4c82711f
)
O=runs/t3_k24_form_${HEAD_SHA:0:7}_$BASE
DEC=experiments/t3_decomp/results/gata4/gata4_decisions.json
RUN_DIR=$HOME/v/virtual-embryo-challenge/runs/t3_decomp/gata4_20261028T122425Z  # read only; fixed, no override
TABLES=$HOME/v/vec-ext/t3/tables  # combine.py's --out; fixed, no override
WT=data/raw/T2/E8.75.h5ad
export NUMPY_MADVISE_HUGEPAGE=0 TMPDIR="$S/tmp" VEC_T3_EXT_TABLES="$TABLES"
PY="uv run --frozen --no-sync python -I -B"
threads() { export OMP_NUM_THREADS=$1 OPENBLAS_NUM_THREADS=$1 MKL_NUM_THREADS=$1; }
stamp() { $PY -c 'import json,sys; s=json.load(open(sys.argv[1])); print(s["commit"], s["clean"])' "$1.stamp.json"; }
one() {  # script form out -> 0, 1 (build failed), 2 (check failed), 3 (the form's own refusal: exit 4), 4 (refused: exit 2)
  local b=0
  $PY "experiments/t3_decomp/$1" build "$DEC" --form "$2" --wt "$WT" --run-dir "$RUN_DIR" \
    --board T3:bcatenin --seed 1000001 --out "$3" > /dev/null 2>&1 || b=$?
  if [ $b = 4 ]; then return 3; fi
  if [ $b = 2 ]; then return 4; fi
  if [ $b != 0 ]; then return 1; fi
  uv run --frozen --no-sync vec check "$3" --board T3:bcatenin > /dev/null 2>&1 || return 2
  return 0
}
sha() { sha256sum "$1" | cut -d' ' -f1; }

[ "$(git rev-parse HEAD)" = "$HEAD_SHA" ] || { echo "refused: HEAD is not $HEAD_SHA" >&2; finish 97; }
[ -z "$(git status --porcelain)" ] || { echo "refused: the tree is not clean" >&2; finish 98; }
{ [ -d "$TABLES" ] && [ ! -L "$TABLES" ]; } || { echo "refused: $TABLES is not a real directory" >&2; finish 96; }
git rev-parse HEAD > "$S/head"
if [ -e "$O" ]; then echo "refused: $O exists" >&2; finish 95; fi
mkdir -p "$O/gate" "$O/a" "$O/b"
: > "$S/SHA256SUMS"
: > "$S/stamps"

threads 8
rc=0
$PY -m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp="$S/tmp/pytest" \
  tests/test_build_k24_form.py tests/test_k24_wider_tables.py tests/test_t3_outside_forms.py \
  experiments/t3_decomp/test_build_submission.py > /dev/null 2>&1 || rc=$?
echo $rc > "$S/tests.exit"
if [ $rc != 0 ]; then finish 10; fi

threads 16
gate() {  # base code_offset
  local f="$O/gate/identity-$1.h5ad" rc=0 h
  one build_outside_form.py "identity-$1" "$f" || rc=$?
  echo $rc > "$S/identity-$1.exit"
  if [ $rc = 1 ] || [ $rc = 3 ] || [ $rc = 4 ]; then finish "$2"; fi
  if [ $rc = 2 ]; then finish $(($2 + 1)); fi
  h=$(sha "$f")
  echo "$h  gate/identity-$1.h5ad" >> "$S/SHA256SUMS"
  echo "gate/identity-$1 $(stamp "$f")" >> "$S/stamps"
  if [ "$h" != "${REF[$1]}" ]; then finish $(($2 + 2)); fi
}
gate s1 20
if [ "$BASE" != s1 ]; then gate "$BASE" 30; fi

for side in a b; do
  f="$O/$side/$FORM.h5ad"
  rc=0
  one build_k24_form.py "$FORM" "$f" || rc=$?
  echo $rc > "$S/$FORM.$side.exit"
  if [ $rc = 4 ]; then finish 46; fi
  if [ $rc = 3 ] && [ $side = a ]; then echo "$FORM" >> "$S/refused"; rm -rf "$S/tmp"; finish 45; fi
  if [ $rc = 1 ] || [ $rc = 3 ]; then if [ $side = a ]; then finish 40; else finish 41; fi; fi
  if [ $rc = 2 ]; then finish 42; fi
  echo "$(sha "$f")  $side/$FORM.h5ad" >> "$S/SHA256SUMS"
  st=$(stamp "$f") || st=""
  echo "$side/$FORM $st" >> "$S/stamps"
  read -r st_commit st_clean _ <<< "$st" || true
  if [ ${#st_commit} -lt 7 ] || [ "${HEAD_SHA#"$st_commit"}" = "$HEAD_SHA" ] || [ "$st_clean" != True ]; then
    finish 44
  fi
done
if [ "$(sha "$O/a/$FORM.h5ad")" != "$(sha "$O/b/$FORM.h5ad")" ]; then finish 43; fi
rm -rf "$S/tmp"
finish 0
