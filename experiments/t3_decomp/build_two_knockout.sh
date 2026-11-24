#!/usr/bin/env bash
# The two-knockout list's T3 files (experiments/t3_decomp/build_two_knockout.py; ledger row R17), built on a DGX
# campaign checkout (~/v/vec-campaign-c, never the reference checkout). Run from the checkout's root, detached
# (setsid nohup ... &), and wait on $S/all.exit with a background poller; never a foreground sleep.
#
#   bash experiments/t3_decomp/build_two_knockout.sh HEAD_SHA BASE
#   bash experiments/t3_decomp/build_two_knockout.sh HEAD_SHA s2 g2
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
#
# Third argument g2 (optional; its only value; BASE s2 only): row R22's second T3 search line, the Gata4-table weight
# on C1 raised to 2 (`--gata4-exponent 2.0 --captain-ruling R22`), as row R30 amends its guard (the figures are
# recorded and stop nothing; the builder's exit 5 is its factor integrity check only). Without it every step, folder
# and command above is as it was. With any third argument S ends in _<argument> (so _g2) and nothing is written in a
# run's folder without one. Refusals before the run (exit 2: a third argument other than g2, another base, a fourth
# argument; 97, 98, 96 as above) write $S/refused.exit only, never a finished run's all.exit; once they pass, all.exit
# holds the word running. The output folder is runs/t3_two_knockout_<first 7>_s2_g2. Steps 0-4 as above (the tests add
# tests/test_t3_gata4_weight.py; gates s1 then s2, exit 20-22 and 30-32), then the rule text's other identity gates,
# without the new option, into the gate folder: 4c. c1-s2, its sha256 must be C1's (exit 60 build, 61 check, 62 hash);
# 4d. c2-s2, its sha256 must be C2's (exit 70, 71, 72). Then 5g: T3_c1_s2_g2.h5ad twice, into folders a and b, with the
# format check, the stamp checks and the pair comparison of step 5 (exit 50 build a, 51 build b, 52 check, 53 the pair
# differs, 54 a stamp not of HEAD_SHA or not clean, 55 the builder's exit 5: the factor applied is not the one its
# figures describe, nothing written, 56 the builder refused (its exit 2: the ruling row, pins, ...)). The c1 and c2
# candidates of step 5 are not built in this mode.
#
# 16 BLAS threads for builds, as the earlier T3 builds; NUMPY_MADVISE_HUGEPAGE=0. The run folder is fixed below (no
# override). Build and check output goes to /dev/null: a test build is read by its exit code, its sha256 and its
# stamp's commit and clean fields only (AGENTS.md), never its log, size or cell count. Written: $S/<step>.exit,
# $S/all.exit, $S/head, $S/SHA256SUMS and $S/stamps. Nothing is deleted (the scratch folder's tmp/ is left in place).
set -u
HEAD_SHA=${1:?usage: build_two_knockout.sh HEAD_SHA BASE [g2]}
BASE=${2:?usage: build_two_knockout.sh HEAD_SHA BASE [g2]}
G2=${3:-}
S=~/v/scratch/t3two_${HEAD_SHA:0:7}_${BASE//[^A-Za-z0-9.]/_}
if [ -n "$G2" ]; then S=${S}_${G2//[^A-Za-z0-9.]/_}; fi
mkdir -p "$S/tmp" || exit 99
EXITF=all.exit  # with a third argument: refused.exit until every check passes
if [ -n "$G2" ]; then EXITF=refused.exit; fi
finish() { trap - ERR TERM INT HUP EXIT; echo "$1" > "$S/$EXITF"; exit "$1"; }
abort() { trap - ERR TERM INT HUP EXIT; echo 99 > "$S/$EXITF"; exit 99; }
trap abort ERR TERM INT HUP EXIT

case "$BASE" in
  s1|s2) CANDS="c1-$BASE c2-$BASE" ;;
  s1x1.5|s1x2.0) CANDS="c2-$BASE" ;;
  *) echo "unknown base $BASE (s1, s1x1.5, s1x2.0 or s2)" >&2; finish 2 ;;
esac
if [ -n "$G2" ]; then
  if [ "$G2" != g2 ]; then echo "unknown third argument $G2 (g2 only)" >&2; finish 2; fi
  if [ "$BASE" != s2 ]; then echo "g2 is for base s2 only" >&2; finish 2; fi
  if [ $# -gt 3 ]; then echo "g2 takes no further argument" >&2; finish 2; fi
fi
declare -A REF=(
  [s1]=3712c3f9988d0f8aa4198b22c37b866a84ef8757fd500a9163944f0663a520a2
  [s1x1.5]=c060fcaa6c7fdc6d6b1a58c8f850a1895a8295abbc50f297de71a32c97f9c420
  [s1x2.0]=7393810fbe3be2cde67cf23ed54f45d22da738d7d3c9e124780a8f6281075661
  [s2]=1b69521e3c28bc84185410e838efc6fcbcff898434e73f9a7547208a4c82711f
)
# the g2 line's other identity gates: C1 and C2 on s2 (experiments/builds/results/t3two_fb01487.SHA256SUMS; the rule text)
C1_S2=8718330866429dc089eaffd2f96f66936fe54134eba2d50f8126266e6a3e54b8
C2_S2=47111ac1d07e73c0c26e990a474dff0de94ae22d6a4474a569bcb674ae2b85a7
G2_FILE=T3_c1_s2_g2.h5ad
O=runs/t3_two_knockout_${HEAD_SHA:0:7}_$BASE
if [ -n "$G2" ]; then O=${O}_g2; fi
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
one_g2() {  # out -> 0, 1 (build failed), 2 (check failed), 4 (refused: exit 2), 5 (the factor check: exit 5)
  local b=0
  $PY experiments/t3_decomp/build_two_knockout.py build "$DEC" --form c1-s2 --gata4-exponent 2.0 \
    --captain-ruling R22 --wt "$WT" --run-dir "$RUN_DIR" --board T3:bcatenin --seed 1000001 --out "$1" \
    > /dev/null 2>&1 || b=$?
  if [ $b = 5 ]; then return 5; fi
  if [ $b = 2 ]; then return 4; fi
  if [ $b != 0 ]; then return 1; fi
  uv run --frozen --no-sync vec check "$1" --board T3:bcatenin > /dev/null 2>&1 || return 2
  return 0
}
sha() { sha256sum "$1" | cut -d' ' -f1; }

[ "$(git rev-parse HEAD)" = "$HEAD_SHA" ] || { echo "refused: HEAD is not $HEAD_SHA" >&2; finish 97; }
[ -z "$(git status --porcelain)" ] || { echo "refused: the tree is not clean" >&2; finish 98; }
if [ -n "$G2" ]; then
  if [ -e "$O" ]; then echo "refused: $O exists" >&2; finish 96; fi
  EXITF=all.exit
  echo running > "$S/all.exit"
fi
git rev-parse HEAD > "$S/head"
if [ -e "$O" ]; then echo "refused: $O exists" >&2; finish 96; fi
mkdir -p "$O/gate" "$O/a" "$O/b"
: > "$S/SHA256SUMS"
: > "$S/stamps"

threads 8
XTESTS=()
if [ -n "$G2" ]; then XTESTS=(tests/test_t3_gata4_weight.py); fi
rc=0
$PY -m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp="$S/tmp/pytest" \
  tests/test_t3_two_knockout_build.py tests/test_t3_barred_exponent_build.py \
  experiments/t3_decomp/test_build_submission.py ${XTESTS[@]+"${XTESTS[@]}"} > /dev/null 2>&1 || rc=$?
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

if [ -n "$G2" ]; then
  gate_form() {  # form reference_sha256 code_offset: a candidate's own command, without the new option
    local f="$O/gate/$1.h5ad" rc=0 h
    one "$1" "$f" || rc=$?
    echo $rc > "$S/$1.exit"
    if [ $rc = 1 ]; then finish "$3"; fi
    if [ $rc = 2 ]; then finish $(($3 + 1)); fi
    h=$(sha "$f")
    echo "$h  gate/$1.h5ad" >> "$S/SHA256SUMS"
    echo "gate/$1 $(stamp "$f")" >> "$S/stamps"
    if [ "$h" != "$2" ]; then finish $(($3 + 2)); fi
  }
  gate_form c1-s2 "$C1_S2" 60
  gate_form c2-s2 "$C2_S2" 70
  for side in a b; do
    f="$O/$side/$G2_FILE"
    rc=0
    one_g2 "$f" || rc=$?
    echo $rc > "$S/$G2_FILE.$side.exit"
    if [ $rc = 5 ]; then finish 55; fi
    if [ $rc = 4 ]; then finish 56; fi
    if [ $rc = 1 ]; then if [ $side = a ]; then finish 50; else finish 51; fi; fi
    if [ $rc = 2 ]; then finish 52; fi
    echo "$(sha "$f")  $side/$G2_FILE" >> "$S/SHA256SUMS"
    st=$(stamp "$f") || st=""
    echo "$side/$G2_FILE $st" >> "$S/stamps"
    read -r st_commit st_clean _ <<< "$st" || true
    # the stamp must be of this commit (its sha a prefix of HEAD_SHA, at least 7 characters) and say clean
    if [ ${#st_commit} -lt 7 ] || [ "${HEAD_SHA#"$st_commit"}" = "$HEAD_SHA" ] || [ "$st_clean" != True ]; then
      finish 54
    fi
  done
  if [ "$(sha "$O/a/$G2_FILE")" != "$(sha "$O/b/$G2_FILE")" ]; then finish 53; fi
  finish 0
fi

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
finish 0
