#!/usr/bin/env bash
# Row R68's four hedge files (the captain's order of 2026-11-25: HEDGE_FIXED_design.md's four files, each row whole
# from one of two parents fixed by recipe), one file per call, each built by experiments/builds/hedge_fixed.py only
# with --captain-ruling R68. Test-board builds on a DGX checkout with its own venv (made there by `uv sync --frozen`),
# never the reference checkout. Run from the root of that checkout, detached, and wait on the exit files with a
# background poller; never a foreground sleep:
#   bash experiments/builds/hedge_fixed_line.sh <COMMIT> <RUN> <NAME> <ANCHOR> <FIRST> <SECOND>
# NAME is one of the four below (tests/test_hedge_fixed.py holds this table equal to the builder's TABLE); ANCHOR,
# FIRST and SECOND are the files in the DGX folders of their records, as paths from the checkout's root or absolute,
# each pinned by the sha256 its record's committed list holds (experiments/builds/results/<record>.SHA256SUMS):
#   T3_c2_s2_mxe2.h5ad        board T3:bcatenin            anchor C2  first C2   second xe2  (row parity)
#   T3_s2_mx_xe2_ck.h5ad      board T3:bcatenin            anchor C2  first xe2  second ck   (row parity)
#   T2embryo_s2x_ovmxso.h5ad  board T2:embryo:test_interp  anchor O   first O    second SO   (row parity)
#   T2embryo_s2x_ovbkso.h5ad  board T2:embryo:test_interp  anchor O   first O    second SO   (spatial blocks)
# with C2 <t3two_fb01487>/c2-s2.h5ad, xe2 <t3k23_49f3fd5>/xall-s2.h5ad, ck <t3ck_20a0e1a>/ck-s2.h5ad,
# O <embovsb_20b2381_w1>/O/T2embryo_s2x_ov.h5ad and SO <embovsb_20b2381_w1>/SO/T2embryo_s2xs_ov.h5ad.
# RUN names this run's scratch folder ~/v/scratch/hedgefix/<RUN> and its output folder runs/hedgefix_<commit7>_<RUN>.
# Neither may exist: a run folder is never reused and nothing is ever deleted here, so a second file, or a second run
# of one, takes a new name. In order, each step first writing the word running to ~/v/scratch/hedgefix/<RUN>/<STEP>.exit,
# then its exit code:
#   gate  HEAD is COMMIT (the full hash), the tracked tree is clean, TRIES.md holds exactly one row R68 whose first cell
#         starts "R68 (the captain's ", says "ruling, no run" and names NAME, experiments/builds/hedge_fixed.py and
#         --captain-ruling R68, each as a whole backticked string; ANCHOR, FIRST and SECOND are files with NAME's
#         pinned sha256 for them; .venv's vec.pth names this checkout's src; then it writes the head file
#   a b   NAME by hedge_fixed.py into runs/hedgefix_<commit7>_<RUN>/a, then into .../b, each then checked by `vec check`
#         on NAME's board, its stamp read for its commit and clean fields only
#   done  one sha256 list over the two built files (in the scratch folder and the output folder), the two compared
#         byte for byte
# Every builder and `vec check` output goes to /dev/null: a test build is read by its exit code, its sha256 and its
# stamp's commit and clean fields (AGENTS.md), never its log, size, sidecar or cell count. Nothing here uploads.
# Exit codes: 90 odd RUN, 89 NAME not one of the four (or not six arguments), 91 run folder exists, 96 HEAD is not
# COMMIT, 98 tracked tree not clean, 95 row R68 not as above, 81 ANCHOR not a file with its pin, 82 FIRST not a file
# with its pin, 83 SECOND not a file with its pin, 88 .venv not this checkout's, 93 output folder exists, 11 build,
# 12 the builder refused (its exit 2), 13 vec check, 14 no stamp, 15 stamp not this commit and clean, 40 hash list,
# 41 a and b differ.
set -u
COMMIT=${1:?usage: hedge_fixed_line.sh COMMIT RUN NAME ANCHOR FIRST SECOND}
RUN=${2:?usage: hedge_fixed_line.sh COMMIT RUN NAME ANCHOR FIRST SECOND}
NAME=${3:?usage: hedge_fixed_line.sh COMMIT RUN NAME ANCHOR FIRST SECOND}
ANCHOR=${4:?usage: hedge_fixed_line.sh COMMIT RUN NAME ANCHOR FIRST SECOND}
FIRST=${5:?usage: hedge_fixed_line.sh COMMIT RUN NAME ANCHOR FIRST SECOND}
SECOND=${6:?usage: hedge_fixed_line.sh COMMIT RUN NAME ANCHOR FIRST SECOND}
[ "$#" = 6 ] || { echo "refused: give COMMIT RUN NAME ANCHOR FIRST SECOND and nothing else" >&2; exit 89; }
case "$RUN" in
  *[!A-Za-z0-9_-]*) echo "refused: RUN takes letters, digits, _ and - only" >&2; exit 90 ;;
esac
C2_SHA=47111ac1d07e73c0c26e990a474dff0de94ae22d6a4474a569bcb674ae2b85a7  # c2-s2.h5ad, t3two_fb01487.SHA256SUMS
XE2_SHA=77ff0b947f7f409c9ce171bdc59a44de3d3c2ef2f3e5a7332011a459ce9ee9c3  # xall-s2.h5ad, t3k23_49f3fd5.SHA256SUMS
CK_SHA=e470aaeeab28918e157fb57cdc771122af35e46f0916cd1c6d63650c31faae26  # ck-s2.h5ad, t3ck_20a0e1a.SHA256SUMS
O_SHA=1cfc91be3a20bfdd1a4f7bcfebd7316e1143c6d4048e578c1f30bcfbe1e0e1ea  # O/T2embryo_s2x_ov.h5ad, embovsb_20b2381_w1
SO_SHA=77062082d97555961901c811d8a1aaf9b5b2403417c19a6988c816ff553f1165  # SO/T2embryo_s2xs_ov.h5ad, embovsb_20b2381_w1
# NAME -> its board and the pins of its anchor, its first and its second parent (the builder's TABLE)
case "$NAME" in
  T3_c2_s2_mxe2.h5ad) BOARD=T3:bcatenin ANCHOR_PIN=$C2_SHA FIRST_PIN=$C2_SHA SECOND_PIN=$XE2_SHA ;;
  T3_s2_mx_xe2_ck.h5ad) BOARD=T3:bcatenin ANCHOR_PIN=$C2_SHA FIRST_PIN=$XE2_SHA SECOND_PIN=$CK_SHA ;;
  T2embryo_s2x_ovmxso.h5ad) BOARD=T2:embryo:test_interp ANCHOR_PIN=$O_SHA FIRST_PIN=$O_SHA SECOND_PIN=$SO_SHA ;;
  T2embryo_s2x_ovbkso.h5ad) BOARD=T2:embryo:test_interp ANCHOR_PIN=$O_SHA FIRST_PIN=$O_SHA SECOND_PIN=$SO_SHA ;;
  *) echo "refused: NAME is not one of row R68's four files" >&2; exit 89 ;;
esac
export PATH=$HOME/.local/bin:$PATH
unset UV_PROJECT_ENVIRONMENT  # uv runs this checkout's own .venv, the one the gate checks
S=$HOME/v/scratch/hedgefix/$RUN
mkdir -p "$HOME/v/scratch/hedgefix"
mkdir "$S" 2> /dev/null || { echo "refused: the run folder $S exists (a run folder is never reused)" >&2; exit 91; }
mkdir "$S/tmp"
export NUMPY_MADVISE_HUGEPAGE=0 TMPDIR="$S/tmp"
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
PY="uv run --frozen --no-sync python -I -B"
STEP=gate
finish() { echo "$1" > "$S/$STEP.exit"; exit "$1"; }
begin() { echo running > "$S/$STEP.exit"; }
step_ok() { echo 0 > "$S/$STEP.exit"; }
begin

head=$(git rev-parse HEAD)
[ "$head" = "$COMMIT" ] || { echo "refused: HEAD is not $COMMIT (give the full hash)" >&2; finish 96; }
if [ -n "$(git status --porcelain --untracked-files=no)" ]; then
  echo "refused: tracked files differ from the commit (the build stamps must say clean)" >&2
  finish 98
fi
[ "$(grep -c "^| R68 " experiments/harness/TRIES.md)" = 1 ] \
  || { echo "refused: TRIES.md does not hold exactly one row R68" >&2; finish 95; }
first=$(grep "^| R68 " experiments/harness/TRIES.md | cut -d'|' -f2)
printf '%s\n' "$first" | grep -q "^ R68 (the captain's " \
  || { echo "refused: row R68's first cell is not the captain's ruling" >&2; finish 95; }
printf '%s\n' "$first" | grep -qF -- "ruling, no run" \
  || { echo "refused: row R68's first cell does not say ruling, no run" >&2; finish 95; }
for named in "$NAME" "experiments/builds/hedge_fixed.py" "--captain-ruling R68"; do
  printf '%s\n' "$first" | grep -o '`[^`]*`' | grep -qxF -- "\`$named\`" \
    || { echo "refused: row R68's first cell does not name \`$named\` as a whole backticked string" >&2; finish 95; }
done
pinned() {  # PATH PIN: PATH is a file whose sha256 is PIN
  [ -f "$1" ] && [ "$(sha256sum < "$1" | cut -d' ' -f1)" = "$2" ]
}
pinned "$ANCHOR" "$ANCHOR_PIN" || { echo "refused: ANCHOR is not a file with its pinned sha256" >&2; finish 81; }
pinned "$FIRST" "$FIRST_PIN" || { echo "refused: FIRST is not a file with its pinned sha256" >&2; finish 82; }
pinned "$SECOND" "$SECOND_PIN" || { echo "refused: SECOND is not a file with its pinned sha256" >&2; finish 83; }
venv=0
for pth in .venv/lib/python3*/site-packages/vec.pth; do
  if [ -f "$pth" ] && [ "$(readlink -f "$(< "$pth")")" = "$(readlink -f src)" ]; then venv=1; fi
done
[ "$venv" = 1 ] || {
  echo "refused: this checkout has no .venv of its own (uv sync --frozen here) whose vec.pth names its src" >&2
  finish 88
}
O=runs/hedgefix_${head:0:7}_$RUN
mkdir -p runs
mkdir "$O" 2> /dev/null || { echo "refused: the output folder $O exists" >&2; finish 93; }
echo "$head" > "$S/head"
step_ok

stamp() {  # the stamp's commit and clean fields only
  $PY -c 'import json,sys; s=json.load(open(sys.argv[1])); print(s["commit"], s["clean"])' "$1.stamp.json"
}
build() {  # folder: NAME by the builder into $O/folder, then `vec check` and the stamp; exit 11 build, 12 refused,
  local d="$O/$1" rc st  # 13 vec check, 14 no stamp, 15 stamp not this commit and clean
  $PY experiments/builds/hedge_fixed.py build --captain-ruling R68 --name "$NAME" --anchor "$ANCHOR" \
    --first "$FIRST" --second "$SECOND" --out "$d/$NAME" > /dev/null 2>&1
  rc=$?
  [ "$rc" = 0 ] || { [ "$rc" = 2 ] && return 12; return 11; }
  uv run --frozen --no-sync vec check "$d/$NAME" --board "$BOARD" > /dev/null 2>&1 || return 13
  st=$(stamp "$d/$NAME" 2> /dev/null) || return 14
  echo "$1/$NAME $st" >> "$S/stamps"
  [ "$st" = "$head True" ] || return 15
}
: > "$S/stamps"
for STEP in a b; do begin; build "$STEP" || finish $?; step_ok; done

STEP=done
begin
(cd "$O" && sha256sum "a/$NAME" "b/$NAME") > "$S/SHA256SUMS" || finish 40
cp "$S/SHA256SUMS" "$O/SHA256SUMS" || finish 40
cmp -s "$O/a/$NAME" "$O/b/$NAME" || finish 41
step_ok
