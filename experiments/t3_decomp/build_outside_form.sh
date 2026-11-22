#!/usr/bin/env bash
# The T3 files of the forms row K22's committed reading names (experiments/t3_decomp/build_outside_form.py), built on a
# DGX campaign checkout (~/v/vec-campaign-c, never the reference checkout). Run from the checkout's root, detached
# (setsid nohup ... &), and wait on $S/all.exit with a background poller; never a foreground sleep.
#
#   bash experiments/t3_decomp/build_outside_form.sh HEAD_SHA BASE KIND [KIND2]
#   bash experiments/t3_decomp/build_outside_form.sh --table-exponent E HEAD_SHA BASE xall
#   bash experiments/t3_decomp/build_outside_form.sh --table-exponent 4.0 --captain-ruling R22 HEAD_SHA s2 xall
#
# --table-exponent E (row K23; optional, first): E is 2.0 or 4.0, the exponent the committed K23 reading names, and the
# one form is xall on BASE (no KIND2); the builder refuses it unless that reading names x_eE. Without the option every
# step, folder and command below is as it was (the K22 forms at exponent 1). With it: S and the run folder end in _eE,
# the targeted tests add tests/test_k23_outside_exponents.py, the identity gates (steps 3 and 4) are unchanged and run
# without the option, and only step 5's builds pass --table-exponent E.
# --captain-ruling R22 (after --table-exponent 4.0; optional): row R22's T3 search line (the captain's order of
# 2026-11-22; the line's rule text in docs/PLAN.md decision log 2026-11-22): xall on the base in force, s2, at 4.0, with
# --captain-ruling R22 passed to the builder, nothing else changed. Only R22, only 4.0, only BASE s2 and xall alone
# (else exit 2). Its checks come before anything is written to all.exit: the call runs at the root of its checkout
# (93), HEAD and a clean tree (97, 98), the tables folder (96), exactly one row R22 in experiments/harness/TRIES.md
# whose first cell starts with "R22 (the captain's order", says "ruling, no run" and names `--table-exponent 4.0` and
# `--captain-ruling R22` (94), and no output folder yet (95); a refused call writes only $S/refused.exit and never
# touches a finished run's all.exit. Once they pass, all.exit holds the word running, then $S/head is written, then
# the steps run. After the identity gates (steps 3 and 4: s1 and s2) one more gate, 4b: the queued outside-knockout
# file (row K23's file on s2) by its own command (--table-exponent 2.0, no ruling) once, its sha256 must be K23_S2
# (exit 50 build or refusal, 51 check, 52 hash). Then step 5 builds xall-s2 with --table-exponent 4.0 --captain-ruling
# R22 twice; the builder's guard refusal (its exit 5: the den-weighted total moves beyond 5%, nothing written) stops
# the run with 47. The targeted tests add tests/test_t3_exponent_line_build.py and
# tests/test_t3_exponent_line_launcher.py. S and the run folder end in _s2_e4.0.
# HEAD_SHA: the full commit the checkout must be at (it holds the committed K22 reading, TABLES.sha256 and the filled
# pins). BASE: s1, s1x1.5, s1x2.0 or s2, the base B of decision log 2026-11-20 (q), fixed there; nothing here chooses
# it. KIND, KIND2: the one or two forms K22's reading names, in its order (xall or xunrel); nothing here chooses them,
# and the builder refuses a kind the reading does not name. The tables folder is fixed below
# (~/v/vec-ext/t3/tables_two_study: the two-study tables rows K22 and K23 read, whose bytes the builder pins through
# outside_forms.py's K22 pins; ~/v/vec-ext/t3/tables holds row K24's wider tables since 2026-11-21) and passed only
# through VEC_T3_EXT_TABLES. uv is put on the path here ($HOME/.local/bin). The scratch folder S
# (~/v/scratch/t3outside_<first 7 of HEAD_SHA>_<BASE>) is made first, so every finish below writes $S/all.exit. Steps, in order, stopping at the first
# failure:
#   0. BASE one of the four and each KIND one of the two, no KIND twice (exit 2);
#   1. HEAD and a clean tree (exit 97, 98), the tables folder a real directory (96), and no output folder yet (95);
#   2. the targeted tests (synthetic data only), 8 BLAS threads (exit 10);
#   3. identity-s1 once: its sha256 must be the reference barred s1's (exit 20 build, 21 check, 22 hash);
#   4. identity of BASE, when BASE is not s1: its sha256 must be the base file's reference hash (exit 30, 31, 32);
#   5. each KIND on BASE twice, into folders a and b, each format-checked (`vec check --board T3:bcatenin`), each
#      stamp's commit checked against HEAD_SHA and its clean field against True, the pair compared (exit 40 build a,
#      41 build b, 42 check, 44 a stamp not of HEAD_SHA or not clean, 43 the pair differs, 46 the builder refused the
#      run (its exit 2: pins, listing, folder, bytes, records or the reading)). A form the builder refuses on its own
#      (its exit 4: the product cannot be made or changes fewer than 10 genes) is not built, is listed in
#      $S/refused, and the next form still builds: per the row, nothing replaces it.
# all.exit: 0 every step passed; 45 every other step passed and a form was refused ($S/refused); 2, 10, 20-22, 30-32,
# 40-44, 46, 95-98 as above; 47, 50-52 and (in $S/refused.exit) 2, 93-98 with --captain-ruling R22; 99 the launcher
# stopped unexpectedly.
# 16 BLAS threads for builds, as the earlier T3 builds; NUMPY_MADVISE_HUGEPAGE=0. The run folder is fixed below (no
# override). Build and check output goes to /dev/null: a test build is read by its exit code, its sha256 and its
# stamp's commit and clean fields only (AGENTS.md), never its log, size or cell count. Written: $S/<step>.exit,
# $S/all.exit, $S/head, $S/SHA256SUMS and $S/stamps. Nothing is deleted.
set -u
EXP=""
RULING=""
if [ "${1:-}" = "--table-exponent" ]; then
  EXP=${2:?usage: build_outside_form.sh --table-exponent E HEAD_SHA BASE xall}
  shift 2
  if [ "${1:-}" = "--captain-ruling" ]; then
    RULING=${2:?usage: build_outside_form.sh --table-exponent 4.0 --captain-ruling R22 HEAD_SHA s2 xall}
    shift 2
  fi
fi
HEAD_SHA=${1:?usage: build_outside_form.sh HEAD_SHA BASE KIND [KIND2]}
BASE=${2:?usage: build_outside_form.sh HEAD_SHA BASE KIND [KIND2]}
KIND1=${3:?usage: build_outside_form.sh HEAD_SHA BASE KIND [KIND2]}
KIND2=${4:-}
S=~/v/scratch/t3outside_${HEAD_SHA:0:7}_${BASE//[^A-Za-z0-9.]/_}
if [ -n "$EXP" ]; then S=${S}_e${EXP//[^0-9.]/_}; fi
mkdir -p "$S/tmp" || exit 99
EXITF=all.exit  # with --captain-ruling: refused.exit until every check passes
if [ -n "$RULING" ]; then EXITF=refused.exit; fi
finish() { trap - ERR TERM INT HUP EXIT; echo "$1" > "$S/$EXITF"; exit "$1"; }
abort() { trap - ERR TERM INT HUP EXIT; echo 99 > "$S/$EXITF"; exit 99; }
trap abort ERR TERM INT HUP EXIT

case "$BASE" in
  s1|s1x1.5|s1x2.0|s2) ;;
  *) echo "unknown base $BASE (s1, s1x1.5, s1x2.0 or s2)" >&2; finish 2 ;;
esac
CANDS=""
for k in "$KIND1" $KIND2; do
  case "$k" in
    xall|xunrel) ;;
    *) echo "unknown form $k (xall or xunrel)" >&2; finish 2 ;;
  esac
  CANDS="$CANDS $k-$BASE"
done
if [ -n "$KIND2" ] && [ "$KIND1" = "$KIND2" ]; then echo "the same form twice" >&2; finish 2; fi
XFLAG=()
if [ -n "$EXP" ]; then
  case "$EXP" in
    2.0|4.0) ;;
    *) echo "unknown table exponent $EXP (2.0 or 4.0)" >&2; finish 2 ;;
  esac
  if [ "$KIND1" != xall ] || [ -n "$KIND2" ]; then echo "--table-exponent builds xall alone" >&2; finish 2; fi
  XFLAG=(--table-exponent "$EXP")
fi
if [ -n "$RULING" ]; then
  if [ "$RULING" != R22 ] || [ "$EXP" != 4.0 ] || [ "$BASE" != s2 ]; then
    echo "--captain-ruling takes R22 with --table-exponent 4.0 on base s2 only" >&2
    finish 2
  fi
  XFLAG=(--table-exponent "$EXP" --captain-ruling "$RULING")
fi
K23_S2=77ff0b947f7f409c9ce171bdc59a44de3d3c2ef2f3e5a7332011a459ce9ee9c3  # row K23's file on s2 (t3k23_49f3fd5)
r22_row() {  # exactly one row R22, its first cell the captain's order naming the option and the flag
  local n first
  n=$(grep -c '^| R22 ' experiments/harness/TRIES.md 2> /dev/null) || return 1
  [ "$n" = 1 ] || return 1
  first=$(grep '^| R22 ' experiments/harness/TRIES.md | cut -d'|' -f2)
  case "$first" in " R22 (the captain's order"*) ;; *) return 1 ;; esac
  case "$first" in *"ruling, no run"*) ;; *) return 1 ;; esac
  case "$first" in *'`--table-exponent 4.0`'*) ;; *) return 1 ;; esac
  case "$first" in *'`--captain-ruling R22`'*) ;; *) return 1 ;; esac
}
declare -A REF=(
  [s1]=3712c3f9988d0f8aa4198b22c37b866a84ef8757fd500a9163944f0663a520a2
  [s1x1.5]=c060fcaa6c7fdc6d6b1a58c8f850a1895a8295abbc50f297de71a32c97f9c420
  [s1x2.0]=7393810fbe3be2cde67cf23ed54f45d22da738d7d3c9e124780a8f6281075661
  [s2]=1b69521e3c28bc84185410e838efc6fcbcff898434e73f9a7547208a4c82711f
)
O=runs/t3_outside_form_${HEAD_SHA:0:7}_$BASE
if [ -n "$EXP" ]; then O=${O}_e$EXP; fi
DEC=experiments/t3_decomp/results/gata4/gata4_decisions.json
RUN_DIR=$HOME/v/virtual-embryo-challenge/runs/t3_decomp/gata4_20261028T122425Z  # read only; fixed, no override
TABLES=$HOME/v/vec-ext/t3/tables_two_study  # the two-study tables K22 and K23 read; fixed, no override
WT=data/raw/T2/E8.75.h5ad
export PATH=$HOME/.local/bin:$PATH  # uv, before its first use
export NUMPY_MADVISE_HUGEPAGE=0 TMPDIR="$S/tmp" VEC_T3_EXT_TABLES="$TABLES"
PY="uv run --frozen --no-sync python -I -B"
threads() { export OMP_NUM_THREADS=$1 OPENBLAS_NUM_THREADS=$1 MKL_NUM_THREADS=$1; }
stamp() { $PY -c 'import json,sys; s=json.load(open(sys.argv[1])); print(s["commit"], s["clean"])' "$1.stamp.json"; }
one() {  # form out [flags] -> 0, 1 (build failed), 2 (check failed), 3 (the form's own refusal: exit 4), 4 (refused: exit 2)
  local b=0
  $PY experiments/t3_decomp/build_outside_form.py build "$DEC" --form "$1" "${@:3}" --wt "$WT" --run-dir "$RUN_DIR" \
    --board T3:bcatenin --seed 1000001 --out "$2" > /dev/null 2>&1 || b=$?
  if [ $b = 4 ]; then return 3; fi
  if [ $b = 5 ]; then return 5; fi
  if [ $b = 2 ]; then return 4; fi
  if [ $b != 0 ]; then return 1; fi
  uv run --frozen --no-sync vec check "$2" --board T3:bcatenin > /dev/null 2>&1 || return 2
  return 0
}
sha() { sha256sum "$1" | cut -d' ' -f1; }

if [ -n "$RULING" ] && [ -n "$(git rev-parse --show-prefix)" ]; then
  echo "refused: run from the root of the checkout" >&2
  finish 93
fi
[ "$(git rev-parse HEAD)" = "$HEAD_SHA" ] || { echo "refused: HEAD is not $HEAD_SHA" >&2; finish 97; }
[ -z "$(git status --porcelain)" ] || { echo "refused: the tree is not clean" >&2; finish 98; }
{ [ -d "$TABLES" ] && [ ! -L "$TABLES" ]; } || { echo "refused: $TABLES is not a real directory" >&2; finish 96; }
if [ -n "$RULING" ]; then
  r22_row || { echo "refused: no single row R22 in TRIES.md is the captain's order naming the option" >&2; finish 94; }
  if [ -e "$O" ]; then echo "refused: $O exists" >&2; finish 95; fi
  EXITF=all.exit
  echo running > "$S/all.exit"
fi
git rev-parse HEAD > "$S/head"
if [ -e "$O" ]; then echo "refused: $O exists" >&2; finish 95; fi
mkdir -p "$O/gate" "$O/a" "$O/b"
: > "$S/SHA256SUMS"
: > "$S/stamps"

threads 8
XTESTS=()
if [ -n "$EXP" ]; then XTESTS=(tests/test_k23_outside_exponents.py); fi
if [ -n "$RULING" ]; then
  XTESTS+=(tests/test_t3_exponent_line_build.py tests/test_t3_exponent_line_launcher.py)
fi
rc=0
$PY -m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp="$S/tmp/pytest" \
  tests/test_t3_outside_forms.py tests/test_k22_outside_forms.py tests/test_t3_two_knockout_build.py \
  tests/test_t3_barred_exponent_build.py experiments/t3_decomp/test_build_submission.py ${XTESTS[@]+"${XTESTS[@]}"} \
  > /dev/null 2>&1 || rc=$?
echo $rc > "$S/tests.exit"
if [ $rc != 0 ]; then finish 10; fi

threads 16
gate() {  # base code_offset
  local f="$O/gate/identity-$1.h5ad" rc=0 h
  one "identity-$1" "$f" || rc=$?
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
if [ -n "$RULING" ]; then  # 4b: row K23's file on s2 by its own command
  f="$O/gate/xall-s2_e2.0.h5ad"
  rc=0
  one xall-s2 "$f" --table-exponent 2.0 || rc=$?
  echo $rc > "$S/xall-s2_e2.0.exit"
  if [ $rc = 2 ]; then finish 51; fi
  if [ $rc != 0 ]; then finish 50; fi
  h=$(sha "$f")
  echo "$h  gate/xall-s2_e2.0.h5ad" >> "$S/SHA256SUMS"
  echo "gate/xall-s2_e2.0 $(stamp "$f")" >> "$S/stamps"
  if [ "$h" != "$K23_S2" ]; then finish 52; fi
fi

refused=0
for form in $CANDS; do
  for side in a b; do
    f="$O/$side/$form.h5ad"
    rc=0
    one "$form" "$f" ${XFLAG[@]+"${XFLAG[@]}"} || rc=$?
    echo $rc > "$S/$form.$side.exit"
    if [ $rc = 4 ]; then finish 46; fi
    if [ $rc = 5 ]; then finish 47; fi
    if [ $rc = 3 ] && [ $side = a ]; then refused=1; echo "$form" >> "$S/refused"; continue 2; fi
    if [ $rc = 1 ] || [ $rc = 3 ]; then if [ $side = a ]; then finish 40; else finish 41; fi; fi
    if [ $rc = 2 ]; then finish 42; fi
    echo "$(sha "$f")  $side/$form.h5ad" >> "$S/SHA256SUMS"
    st=$(stamp "$f") || st=""
    echo "$side/$form $st" >> "$S/stamps"
    read -r st_commit st_clean _ <<< "$st" || true
    if [ ${#st_commit} -lt 7 ] || [ "${HEAD_SHA#"$st_commit"}" = "$HEAD_SHA" ] || [ "$st_clean" != True ]; then
      finish 44
    fi
  done
  if [ "$(sha "$O/a/$form.h5ad")" != "$(sha "$O/b/$form.h5ad")" ]; then finish 43; fi
done
if [ $refused = 1 ]; then finish 45; fi
finish 0
