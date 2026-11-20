#!/usr/bin/env bash
# The two-knockout list's T3 files (experiments/t3_decomp/build_two_knockout.py; ledger row R17), built on a DGX
# campaign checkout (~/v/vec-campaign-c, never the reference checkout). Run from the checkout's root, detached
# (setsid nohup ... &), and wait on $S/all.exit with a background poller; never a foreground sleep.
#
#   bash experiments/t3_decomp/build_two_knockout.sh HEAD_SHA BASE
#
# HEAD_SHA: the full commit the checkout must be at. BASE: s1, s1x1.5, s1x2.0 or s2, the base B of the list, fixed by
# the decision log once both uploads of entry 2026-11-20 (k) have returned; nothing here chooses it. The scratch folder
# S (~/v/scratch/t3two_<first 7 of HEAD_SHA>_<BASE>) is made first, so every finish below writes $S/all.exit. Steps,
# in order, stopping at the first failure:
#   0. BASE one of the four (exit 2);
#   1. HEAD and a clean tree (exit 97, 98), and no output folder yet (exit 96);
#   2. the targeted tests (synthetic data only), 8 BLAS threads (exit 10);
#   3. identity-s1 once: its sha256 must be the reference barred s1's (exit 20 build, 21 check, 22 hash);
#   4. identity of BASE, when BASE is not s1: its sha256 must be the base file's reference hash (exit 30, 31, 32);
#   5. each candidate of BASE (c1 on s1 and s2 only; c2 on every base) twice, into folders a and b, each format-checked
#      as build_all.sh checks T3 files (`vec check --board T3:bcatenin`), each stamp's commit checked against HEAD_SHA
#      and its clean field against True, the pair compared (exit 40 build a, 41 build b, 42 check, 44 a stamp not of
#      HEAD_SHA or not clean, 43 the pair differs).
# all.exit: 0 every step passed; 2, 10, 20-22, 30-32, 40-44, 96-98 as above; 99 the launcher stopped unexpectedly (a
# failing command, a TERM, INT or HUP signal, or an exit that was not one of the finishes above).
# 16 BLAS threads for builds, as the earlier T3 builds; NUMPY_MADVISE_HUGEPAGE=0. The run folder is fixed below (no
# override). Build and check output goes to /dev/null: a test build is read by its exit code, its sha256 and its
# stamp's commit and clean fields only (AGENTS.md), never its log, size or cell count. Written: $S/<step>.exit,
# $S/all.exit, $S/head, $S/SHA256SUMS and $S/stamps.
set -u
HEAD_SHA=${1:?usage: build_two_knockout.sh HEAD_SHA BASE}
BASE=${2:?usage: build_two_knockout.sh HEAD_SHA BASE}
S=~/v/scratch/t3two_${HEAD_SHA:0:7}_${BASE//[^A-Za-z0-9.]/_}
mkdir -p "$S/tmp" || exit 99
finish() { trap - ERR TERM INT HUP EXIT; echo "$1" > "$S/all.exit"; exit "$1"; }
abort() { trap - ERR TERM INT HUP EXIT; echo 99 > "$S/all.exit"; exit 99; }
trap abort ERR TERM INT HUP EXIT

case "$BASE" in
  s1|s2) CANDS="c1-$BASE c2-$BASE" ;;
  s1x1.5|s1x2.0) CANDS="c2-$BASE" ;;
  *) echo "unknown base $BASE (s1, s1x1.5, s1x2.0 or s2)" >&2; finish 2 ;;
esac
declare -A REF=(
  [s1]=3712c3f9988d0f8aa4198b22c37b866a84ef8757fd500a9163944f0663a520a2
  [s1x1.5]=c060fcaa6c7fdc6d6b1a58c8f850a1895a8295abbc50f297de71a32c97f9c420
  [s1x2.0]=7393810fbe3be2cde67cf23ed54f45d22da738d7d3c9e124780a8f6281075661
  [s2]=1b69521e3c28bc84185410e838efc6fcbcff898434e73f9a7547208a4c82711f
)
O=runs/t3_two_knockout_${HEAD_SHA:0:7}_$BASE
DEC=experiments/t3_decomp/results/gata4/gata4_decisions.json
RUN_DIR=$HOME/v/virtual-embryo-challenge/runs/t3_decomp/gata4_20261028T122425Z  # read only; fixed, no override
WT=data/raw/T2/E8.75.h5ad
export NUMPY_MADVISE_HUGEPAGE=0 TMPDIR="$S/tmp"
PY="uv run --frozen --no-sync python -I -B"
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
if [ -e "$O" ]; then echo "refused: $O exists" >&2; finish 96; fi
mkdir -p "$O/gate" "$O/a" "$O/b"
: > "$S/SHA256SUMS"
: > "$S/stamps"

threads 8
rc=0
$PY -m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp="$S/tmp/pytest" \
  tests/test_t3_two_knockout_build.py tests/test_t3_barred_exponent_build.py \
  experiments/t3_decomp/test_build_submission.py > /dev/null 2>&1 || rc=$?
echo $rc > "$S/tests.exit"
if [ $rc != 0 ]; then finish 10; fi

threads 16
gate() {  # base code_offset
  local f="$O/gate/identity-$1.h5ad" rc=0 h
  one "identity-$1" "$f" || rc=$?
  echo $rc > "$S/identity-$1.exit"
  if [ $rc = 1 ]; then finish "$2"; fi
  if [ $rc = 2 ]; then finish $(($2 + 1)); fi
  h=$(sha "$f")
  echo "$h  gate/identity-$1.h5ad" >> "$S/SHA256SUMS"
  echo "gate/identity-$1 $(stamp "$f")" >> "$S/stamps"
  if [ "$h" != "${REF[$1]}" ]; then finish $(($2 + 2)); fi
}
gate s1 20
if [ "$BASE" != s1 ]; then gate "$BASE" 30; fi

for form in $CANDS; do
  for side in a b; do
    f="$O/$side/$form.h5ad"
    rc=0
    one "$form" "$f" || rc=$?
    echo $rc > "$S/$form.$side.exit"
    if [ $rc = 1 ]; then if [ $side = a ]; then finish 40; else finish 41; fi; fi
    if [ $rc = 2 ]; then finish 42; fi
    echo "$(sha "$f")  $side/$form.h5ad" >> "$S/SHA256SUMS"
    st=$(stamp "$f") || st=""
    echo "$side/$form $st" >> "$S/stamps"
    read -r st_commit st_clean _ <<< "$st" || true
    # the stamp must be of this commit (its sha a prefix of HEAD_SHA, at least 7 characters) and say clean
    if [ ${#st_commit} -lt 7 ] || [ "${HEAD_SHA#"$st_commit"}" = "$HEAD_SHA" ] || [ "$st_clean" != True ]; then
      finish 44
    fi
  done
  if [ "$(sha "$O/a/$form.h5ad")" != "$(sha "$O/b/$form.h5ad")" ]; then finish 43; fi
done
rm -rf "$S/tmp"
finish 0
