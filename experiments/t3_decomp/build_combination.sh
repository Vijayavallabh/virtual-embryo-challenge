#!/usr/bin/env bash
# The combination file CK of row R18 (experiments/t3_decomp/build_combination.py), built on a DGX campaign checkout
# (~/v/vec-campaign-c, never the reference checkout). Run from the checkout's root, detached (setsid nohup ... &), and
# wait on $S/all.exit with a background poller; never a foreground sleep.
#
#   bash experiments/t3_decomp/build_combination.sh HEAD_SHA BASE
#
# HEAD_SHA: the full commit the checkout must be at (a main commit holding this builder, row R18 and its decision-log
# entry). BASE: s1, s1x1.5, s1x2.0 or s2; CK is built ahead for every base, and the base B of decision log 2026-11-20
# (q) is fixed there; nothing here chooses it. The tables folder is fixed below (~/v/vec-ext/t3/tables, the two-study
# copy the pins accept) and passed only through VEC_T3_EXT_TABLES. The scratch folder S (~/v/scratch/t3ck_<first 7 of
# HEAD_SHA>_<BASE>) is made first, so every finish below writes $S/all.exit. Steps, in order, stopping at the first
# failure:
#   0. BASE one of the four (exit 2);
#   1. HEAD and a clean tree (exit 97, 98), the tables folder a real directory (96), and no output folder yet (95);
#   2. the targeted tests (synthetic data only), 8 BLAS threads (exit 10);
#   3. identity-s1 once: its sha256 must be the reference barred s1's (exit 20 build, 21 check, 22 hash);
#   4. identity of BASE, when BASE is not s1: its sha256 must be the base file's reference hash (exit 30, 31, 32);
#   5. gate-c2-BASE (table exponent 0): its sha256 must be C2 on BASE as experiments/builds/results/
#      t3two_fb01487.SHA256SUMS records it (exit 50 build, 51 check, 52 hash);
#   6. gate-k23-BASE (no Gata4 table, table exponent 2): its sha256 must be the K23 file on BASE as experiments/builds/
#      results/t3k23_49f3fd5.SHA256SUMS records it (exit 60 build, 61 check, 62 hash);
#   7. ck-BASE twice, into folders a and b, each format-checked (`vec check --board T3:bcatenin`), each stamp's commit
#      checked against HEAD_SHA and its clean field against True, the pair compared (exit 40 build a, 41 build b, 42
#      check, 44 a stamp not of HEAD_SHA or not clean, 43 the pair differs, 46 the builder refused the run (its exit
#      2)). The form's own refusal (its exit 4) is listed in $S/refused, nothing is built in its place (exit 45).
# all.exit: 0 every step passed; 2, 10, 20-22, 30-32, 40-46, 50-52, 60-62, 95-98 as above; 99 the launcher stopped
# unexpectedly. 16 BLAS threads for builds, as the earlier T3 builds; NUMPY_MADVISE_HUGEPAGE=0. The run folder is fixed
# below (no override). Build and check output goes to /dev/null: a test build is read by its exit code, its sha256 and
# its stamp's commit and clean fields only (AGENTS.md), never its log, size or cell count. Written: $S/<step>.exit,
# $S/all.exit, $S/head, $S/SHA256SUMS and $S/stamps.
set -u
HEAD_SHA=${1:?usage: build_combination.sh HEAD_SHA BASE}
BASE=${2:?usage: build_combination.sh HEAD_SHA BASE}
S=~/v/scratch/t3ck_${HEAD_SHA:0:7}_${BASE//[^A-Za-z0-9.]/_}
mkdir -p "$S/tmp" || exit 99
finish() { trap - ERR TERM INT HUP EXIT; echo "$1" > "$S/all.exit"; exit "$1"; }
abort() { trap - ERR TERM INT HUP EXIT; echo 99 > "$S/all.exit"; exit 99; }
trap abort ERR TERM INT HUP EXIT

case "$BASE" in
  s1|s1x1.5|s1x2.0|s2) ;;
  *) echo "unknown base $BASE (s1, s1x1.5, s1x2.0 or s2)" >&2; finish 2 ;;
esac
declare -A REF=(
  [s1]=3712c3f9988d0f8aa4198b22c37b866a84ef8757fd500a9163944f0663a520a2
  [s1x1.5]=c060fcaa6c7fdc6d6b1a58c8f850a1895a8295abbc50f297de71a32c97f9c420
  [s1x2.0]=7393810fbe3be2cde67cf23ed54f45d22da738d7d3c9e124780a8f6281075661
  [s2]=1b69521e3c28bc84185410e838efc6fcbcff898434e73f9a7547208a4c82711f
)
declare -A C2REF=(
  [s1]=047506dbab19f007ce45820e3ab79c57b760fbbfc8da53058e4748b8c9930eac
  [s1x1.5]=12a0dffbfc4fa0cdbedaef992edf17778e2610792008cc114e12b02569f8c1c8
  [s1x2.0]=9559c9cd785b2eefc07847d4094d2d002daaac6a03bc2ba422eb00cd9e1d8d00
  [s2]=47111ac1d07e73c0c26e990a474dff0de94ae22d6a4474a569bcb674ae2b85a7
)
declare -A K23REF=(
  [s1]=f3af008c12e6a3f29e20a3fc770a43671012923d8dc25fbfada29138f1d8854f
  [s1x1.5]=0b31c03e301fca601957f06cbe1d31a204e425a287a70394e56b4a6ad9ee7934
  [s1x2.0]=7e2815bd6f1d7c43b116e22e9976d0629823b15718dda5e1b3a9cde5b5ef82f5
  [s2]=77ff0b947f7f409c9ce171bdc59a44de3d3c2ef2f3e5a7332011a459ce9ee9c3
)
O=runs/t3_combination_${HEAD_SHA:0:7}_$BASE
DEC=experiments/t3_decomp/results/gata4/gata4_decisions.json
RUN_DIR=$HOME/v/virtual-embryo-challenge/runs/t3_decomp/gata4_20261028T122425Z  # read only; fixed, no override
TABLES=$HOME/v/vec-ext/t3/tables  # combine.py's --out; fixed, no override
WT=data/raw/T2/E8.75.h5ad
export NUMPY_MADVISE_HUGEPAGE=0 TMPDIR="$S/tmp" VEC_T3_EXT_TABLES="$TABLES"
PY="uv run --frozen --no-sync python -I -B"
threads() { export OMP_NUM_THREADS=$1 OPENBLAS_NUM_THREADS=$1 MKL_NUM_THREADS=$1; }
stamp() { $PY -c 'import json,sys; s=json.load(open(sys.argv[1])); print(s["commit"], s["clean"])' "$1.stamp.json"; }
one() {  # form out -> 0, 1 (build failed), 2 (check failed), 3 (the form's own refusal: exit 4), 4 (refused: exit 2)
  local b=0
  $PY experiments/t3_decomp/build_combination.py build "$DEC" --form "$1" --wt "$WT" --run-dir "$RUN_DIR" \
    --board T3:bcatenin --seed 1000001 --out "$2" > /dev/null 2>&1 || b=$?
  if [ $b = 4 ]; then return 3; fi
  if [ $b = 2 ]; then return 4; fi
  if [ $b != 0 ]; then return 1; fi
  uv run --frozen --no-sync vec check "$2" --board T3:bcatenin > /dev/null 2>&1 || return 2
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
  tests/test_t3_combination_build.py tests/test_t3_two_knockout_build.py tests/test_t3_outside_forms.py \
  tests/test_k23_outside_exponents.py tests/test_k22_outside_forms.py tests/test_t3_barred_exponent_build.py \
  experiments/t3_decomp/test_build_submission.py > /dev/null 2>&1 || rc=$?
echo $rc > "$S/tests.exit"
if [ $rc != 0 ]; then finish 10; fi

threads 16
gate() {  # form reference_sha code_offset
  local f="$O/gate/$1.h5ad" rc=0 h
  one "$1" "$f" || rc=$?
  echo $rc > "$S/$1.exit"
  if [ $rc = 1 ] || [ $rc = 3 ] || [ $rc = 4 ]; then finish "$3"; fi
  if [ $rc = 2 ]; then finish $(($3 + 1)); fi
  h=$(sha "$f")
  echo "$h  gate/$1.h5ad" >> "$S/SHA256SUMS"
  echo "gate/$1 $(stamp "$f")" >> "$S/stamps"
  if [ "$h" != "$2" ]; then finish $(($3 + 2)); fi
}
gate identity-s1 "${REF[s1]}" 20
if [ "$BASE" != s1 ]; then gate "identity-$BASE" "${REF[$BASE]}" 30; fi
gate "gate-c2-$BASE" "${C2REF[$BASE]}" 50
gate "gate-k23-$BASE" "${K23REF[$BASE]}" 60

form=ck-$BASE
for side in a b; do
  f="$O/$side/$form.h5ad"
  rc=0
  one "$form" "$f" || rc=$?
  echo $rc > "$S/$form.$side.exit"
  if [ $rc = 4 ]; then finish 46; fi
  if [ $rc = 3 ] && [ $side = a ]; then echo "$form" >> "$S/refused"; rm -rf "$S/tmp"; finish 45; fi
  if [ $rc = 1 ] || [ $rc = 3 ]; then if [ $side = a ]; then finish 40; else finish 41; fi; fi
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
rm -rf "$S/tmp"
finish 0
