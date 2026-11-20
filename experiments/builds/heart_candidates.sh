#!/usr/bin/env bash
# Row R16's heart candidates (experiments/harness/TRIES.md, as proposed): check builds on a DGX campaign checkout
# (~/v/vec-campaign, never the reference checkout). Run from the checkout's root, detached, and wait on the exit files
# with a background poller; never a foreground sleep:
#   bash experiments/builds/heart_candidates.sh <COMMIT>
# In order, each step writing ~/v/scratch/heartcand/<STEP>.exit:
#   gate   HEAD is COMMIT, the tracked tree is clean, docs/S2_RULES.sha256 checks, row R16 is in TRIES.md and
#          docs/S2_RULES.md names R16 (the record comes before any build meant for upload)
#   tests  the heart builder's targeted tests (synthetic stages only), 8 BLAS threads
#   C      the best file so far (row R10, set heart3) rebuilt once; its half must give HEART3_HALF or everything stops
#   K1a K1b  the half trim (form B), twice, into separate folders
#   K2a K2b  the full-scope copy (form A), twice; its half is also compared with set heart2's HEART2_HALF
#   done   one sha256 list over every built file, a and b compared byte for byte
# Builds run at 16 BLAS threads, as the reference builds. Every builder, halve.py and `vec check` output goes to
# /dev/null: a test build is read by its exit code, its sha256 and its stamp's commit and clean fields only (AGENTS.md),
# never its log, size or cell count. Nothing here uploads anything.
set -u
COMMIT=${1:?usage: heart_candidates.sh <COMMIT>}
export PATH=$HOME/.local/bin:$PATH
S=~/v/scratch/heartcand
mkdir -p "$S/tmp"
rm -f "$S"/*.exit "$S/head" "$S/stamps" "$S/SHA256SUMS" "$S/tests.out"  # no exit file of an earlier run is left
export NUMPY_MADVISE_HUGEPAGE=0 TMPDIR="$S/tmp"
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
PY="uv run --frozen --no-sync python -I -B"
BOARD=T2:heart:test_extrap
HEART3_HALF=04cde6361c04840ee47fff32eb9a20f4781bd649387cc540ed9fdea71e65fa35  # set heart3, row R10 (C's half)
HEART2_HALF=dc436d9225bb28e7997d2ec51e95254cd125e569d44c059596298d0a7b6af5c2  # set heart2 (A's half)
C_ARGS="--entry core_trim --expression copy --z-mode raw"
K1_ARGS="--entry s2_trim --expression copy --z-mode raw"
K2_ARGS="--entry s1 --expression copy --z-mode raw"
STEP=gate
finish() { echo "$1" > "$S/$STEP.exit"; exit "$1"; }
step_ok() { echo 0 > "$S/$STEP.exit"; }

head=$(git rev-parse HEAD)
case "$head" in "$COMMIT"*) ;; *) echo "refused: HEAD is not $COMMIT" >&2; finish 96 ;; esac
if [ -n "$(git status --porcelain --untracked-files=no)" ]; then
  echo "refused: tracked files differ from the commit (the build stamps must say clean)" >&2
  finish 98
fi
sha256sum -c docs/S2_RULES.sha256 > /dev/null 2>&1 || finish 97
grep -q '^| R16 ' experiments/harness/TRIES.md || { echo "refused: row R16 is not in TRIES.md" >&2; finish 95; }
grep -q 'R16' docs/S2_RULES.md || { echo "refused: no docs/S2_RULES.md amendment names R16" >&2; finish 94; }
O=runs/heartcand_${head:0:7}
rm -rf "$O" && mkdir -p "$O"
echo "$head" > "$S/head"
step_ok

STEP=tests
$PY -m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp="$S/tmp/pytest" \
  tests/test_heart_candidates.py tests/test_build_heart_core_trim.py tests/test_build_heart_iso_trim.py \
  tests/test_build_heart_floor_iso.py tests/test_heart_trend_build.py tests/test_halve.py tests/test_hidden_counts.py \
  tests/test_s2_rules.py > "$S/tests.out" 2>&1 || finish 10
step_ok

export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16
stamp() {  # the stamp's commit and clean fields only
  $PY -c 'import json,sys; s=json.load(open(sys.argv[1])); print(s["commit"], s["clean"])' "$1.stamp.json"
}
build() {  # folder name args: the build, the halving, `vec check` on both; exit 11 build, 12 halve, 13 check
  local d="$O/$1" name=$2
  shift 2
  mkdir -p "$d"
  $PY experiments/builds/build_test.py heart "$@" --out "$d/$name.h5ad" > /dev/null 2>&1 || return 11
  $PY experiments/builds/halve.py "$d/$name.h5ad" "$d/${name}_half.h5ad" > /dev/null 2>&1 || return 12
  for f in "$d/$name.h5ad" "$d/${name}_half.h5ad"; do
    uv run --frozen --no-sync vec check "$f" --board "$BOARD" > /dev/null 2>&1 || return 13
  done
  echo "$1/$name $(stamp "$d/$name.h5ad")" >> "$S/stamps"
}
: > "$S/stamps"

STEP=C
build C T2heart_s2core $C_ARGS || finish $?
got=$(sha256sum "$O/C/T2heart_s2core_half.h5ad" | cut -d' ' -f1)
[ "$got" = "$HEART3_HALF" ] || { echo "stopped: the best file's half does not reproduce set heart3" >&2; finish 20; }
step_ok

for STEP in K1a K1b; do build "$STEP" T2heart_s2trimcopy $K1_ARGS || finish $?; step_ok; done
for STEP in K2a K2b; do build "$STEP" T2heart_s2iso $K2_ARGS || finish $?; step_ok; done

STEP=done
(cd "$O" && sha256sum C/*.h5ad K1a/*.h5ad K1b/*.h5ad K2a/*.h5ad K2b/*.h5ad) > "$S/SHA256SUMS"
cp "$S/SHA256SUMS" "$O/SHA256SUMS"
for f in T2heart_s2trimcopy T2heart_s2trimcopy_half; do cmp -s "$O/K1a/$f.h5ad" "$O/K1b/$f.h5ad" || finish 31; done
for f in T2heart_s2iso T2heart_s2iso_half; do cmp -s "$O/K2a/$f.h5ad" "$O/K2b/$f.h5ad" || finish 32; done
got=$(sha256sum "$O/K2a/T2heart_s2iso_half.h5ad" | cut -d' ' -f1)
[ "$got" = "$HEART2_HALF" ] || { echo "the full-scope half differs from set heart2's" >&2; finish 33; }
step_ok
