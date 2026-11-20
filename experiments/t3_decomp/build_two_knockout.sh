#!/usr/bin/env bash
# The two-knockout list's T3 files (experiments/t3_decomp/build_two_knockout.py; ledger row R17), built on a DGX
# campaign checkout (~/v/vec-campaign-c, never the reference checkout). Run from the checkout's root, detached
# (setsid nohup ... &), and wait on $S/all.exit with a background poller; never a foreground sleep.
#
#   bash experiments/t3_decomp/build_two_knockout.sh HEAD_SHA BASE
#
# HEAD_SHA: the full commit the checkout must be at. BASE: s1, s1x1.5, s1x2.0 or s2, the base B of the list, fixed by
# the decision log once both uploads of entry 2026-11-20 (k) have returned; nothing here chooses it. Steps, in order,
# stopping at the first failure:
#   1. HEAD and a clean tree (exit 97, 98);
#   2. the targeted tests (synthetic data only), 8 BLAS threads (exit 10);
#   3. identity-s1 once: its sha256 must be the reference barred s1's (exit 20 build, 21 check, 22 hash);
#   4. identity of BASE, when BASE is not s1: its sha256 must be the base file's reference hash (exit 30, 31, 32);
#   5. each candidate of BASE (c1 on s1 and s2 only; c2 on every base) twice, into folders a and b, each format-checked
#      as build_all.sh checks T3 files (`vec check --board T3:bcatenin`), the pair compared (exit 40 build a, 41 build
#      b, 42 check, 43 the pair differs).
# 16 BLAS threads for builds, as the earlier T3 builds; NUMPY_MADVISE_HUGEPAGE=0. Build and check output goes to
# /dev/null: a test build is read by its exit code, its sha256 and its stamp's commit and clean fields only (AGENTS.md),
# never its log, size or cell count. Written: $S/<step>.exit, $S/all.exit, $S/head, $S/SHA256SUMS and $S/stamps.
set -u
HEAD_SHA=${1:?usage: build_two_knockout.sh HEAD_SHA BASE}
BASE=${2:?usage: build_two_knockout.sh HEAD_SHA BASE}
case "$BASE" in
  s1|s2) CANDS="c1-$BASE c2-$BASE" ;;
  s1x1.5|s1x2.0) CANDS="c2-$BASE" ;;
  *) echo "unknown base $BASE (s1, s1x1.5, s1x2.0 or s2)" >&2; exit 2 ;;
esac
declare -A REF=(
  [s1]=3712c3f9988d0f8aa4198b22c37b866a84ef8757fd500a9163944f0663a520a2
  [s1x1.5]=c060fcaa6c7fdc6d6b1a58c8f850a1895a8295abbc50f297de71a32c97f9c420
  [s1x2.0]=7393810fbe3be2cde67cf23ed54f45d22da738d7d3c9e124780a8f6281075661
  [s2]=1b69521e3c28bc84185410e838efc6fcbcff898434e73f9a7547208a4c82711f
)
S=~/v/scratch/t3two_${HEAD_SHA:0:7}_$BASE
O=runs/t3_two_knockout_${HEAD_SHA:0:7}_$BASE
DEC=experiments/t3_decomp/results/gata4/gata4_decisions.json
RUN_DIR=${T3_RUN_DIR:-$HOME/v/virtual-embryo-challenge/runs/t3_decomp/gata4_20261028T122425Z}  # read only
WT=data/raw/T2/E8.75.h5ad
mkdir -p "$S/tmp"
export NUMPY_MADVISE_HUGEPAGE=0 TMPDIR="$S/tmp"
PY="uv run --frozen --no-sync python -I -B"
finish() { echo "$1" > "$S/all.exit"; exit "$1"; }
threads() { export OMP_NUM_THREADS=$1 OPENBLAS_NUM_THREADS=$1 MKL_NUM_THREADS=$1; }
stamp() { $PY -c 'import json,sys; s=json.load(open(sys.argv[1])); print(s["commit"], s["clean"])' "$1.stamp.json"; }
one() {  # form out -> 0 or 1 (build), 2 (check)
  $PY experiments/t3_decomp/build_two_knockout.py build "$DEC" --form "$1" --wt "$WT" --run-dir "$RUN_DIR" \
    --board T3:bcatenin --seed 1000001 --out "$2" > /dev/null 2>&1 || return 1
  uv run --frozen --no-sync vec check "$2" --board T3:bcatenin > /dev/null 2>&1 || return 2
  return 0
}
sha() { sha256sum "$1" | cut -d' ' -f1; }

[ "$(git rev-parse HEAD)" = "$HEAD_SHA" ] || { echo "refused: HEAD is not $HEAD_SHA" >&2; finish 97; }
[ -z "$(git status --porcelain)" ] || { echo "refused: the tree is not clean" >&2; finish 98; }
git rev-parse HEAD > "$S/head"
[ -e "$O" ] && { echo "refused: $O exists" >&2; finish 96; }
mkdir -p "$O/gate" "$O/a" "$O/b"
: > "$S/SHA256SUMS"
: > "$S/stamps"

threads 8
$PY -m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp="$S/tmp/pytest" \
  tests/test_t3_two_knockout_build.py tests/test_t3_barred_exponent_build.py \
  experiments/t3_decomp/test_build_submission.py > /dev/null 2>&1
rc=$?; echo $rc > "$S/tests.exit"; [ $rc = 0 ] || finish 10

threads 16
gate() {  # base code_offset
  local f="$O/gate/identity-$1.h5ad"
  one "identity-$1" "$f"; local rc=$?
  echo $rc > "$S/identity-$1.exit"
  [ $rc = 1 ] && finish "$2"
  [ $rc = 2 ] && finish $(($2 + 1))
  local h; h=$(sha "$f")
  echo "$h  gate/identity-$1.h5ad" >> "$S/SHA256SUMS"
  echo "gate/identity-$1 $(stamp "$f")" >> "$S/stamps"
  [ "$h" = "${REF[$1]}" ] || finish $(($2 + 2))
}
gate s1 20
[ "$BASE" = s1 ] || gate "$BASE" 30

for form in $CANDS; do
  for side in a b; do
    one "$form" "$O/$side/$form.h5ad"; rc=$?
    echo $rc > "$S/$form.$side.exit"
    [ $rc = 1 ] && { [ $side = a ] && finish 40 || finish 41; }
    [ $rc = 2 ] && finish 42
    echo "$(sha "$O/$side/$form.h5ad")  $side/$form.h5ad" >> "$S/SHA256SUMS"
    echo "$side/$form $(stamp "$O/$side/$form.h5ad")" >> "$S/stamps"
  done
  [ "$(sha "$O/a/$form.h5ad")" = "$(sha "$O/b/$form.h5ad")" ] || finish 43
done
rm -rf "$S/tmp"
finish 0
