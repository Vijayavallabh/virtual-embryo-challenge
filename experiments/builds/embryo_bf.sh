#!/usr/bin/env bash
# Row R71's files (experiments/harness/TRIES.md, the captain's order of 2026-11-25: "My named exception on the embryo's
# overlay form and the z-extent frozen rule for each BF file, as R59, with the challenger's fixes."): the embryo's
# consensus-form files BF on O and on SO only, T2embryo_s2x_ovbf.h5ad (O's command with --shape-interp overlay_bayes
# --captain-ruling R71) and T2embryo_s2xs_ovbf.h5ad (SO's command with the same), each built twice at one clean commit,
# and only the files their gate row E59 (experiments/coords/bf_read.py) passed: row E59's outcome cell, read by code,
# holds one of the read's PASS texts, and the committed reading says the same PASS. Test-board builds on a DGX campaign
# checkout (~/v/vec-campaign, never the reference checkout). Run from the root of the checkout or detached worktree it is
# called in, detached, and wait on the exit files with a background poller; never a foreground sleep:
#   bash experiments/builds/embryo_bf.sh <COMMIT> <RUN>
# COMMIT is the full hash of the build commit. RUN names this run's own scratch folder ~/v/scratch/embbf/<RUN> and its
# output folder runs/embbf_<commit7>_<RUN>. Neither may exist: a run folder is never reused and nothing is ever deleted
# here, so a second run takes a new name. In order, each step writing ~/v/scratch/embbf/<RUN>/<STEP>.exit (the word
# running first, then its exit code):
#   gate     HEAD is COMMIT (else 96), the tracked tree is clean (98), docs/S2_RULES.sha256 checks (97), TRIES.md holds
#            exactly one row R71 whose first cell is the captain's order ("ruling, no run") naming BF on O and on SO
#            only and `--captain-ruling R71`, the builder's own check (95), exactly one well-formed row E59 (94) whose
#            outcome cell holds one PASS text and no STOP or VOID words (92; experiments/coords/bf_row.py passed, which
#            names the files), and the gate's reading is committed (89: READING_DIR's reading.txt and reading.json
#            tracked at this commit) and says the same PASS on a clean tree at its run commit (88), an ancestor of HEAD
#            (87) whose TRIES.md holds a well-formed row E59 (86). The identity pins of O and SO are read from their
#            committed SHA256SUMS lists (85; never from a .txt record) and must equal the ones written here
#   tests    the builder's, the gate's and this launcher's targeted tests (synthetic stages only), 8 BLAS threads
#   O        the identity file O by its committed command; it must give O_SHA or everything stops (21)
#   SO       the identity file SO by its committed command; it must give SO_SHA or everything stops (22)
#   BFOa BFOb    T2embryo_s2x_ovbf.h5ad by its command, twice, into separate folders (when the gate passed it)
#   BFSOa BFSOb  T2embryo_s2xs_ovbf.h5ad by its command, twice, into separate folders (when the gate passed it)
#   done     one sha256 list over every built file; each pair compared byte for byte (31)
# Where the identity pins and their commands come from:
#   O   1cfc91be3a20...: experiments/builds/results/embov_3b813b5.SHA256SUMS, OVa/T2embryo_s2x_ov.h5ad (built by
#       experiments/builds/embryo_overlay.sh at 3b813b5); its command is that launcher's OV_ARGS
#   SO  77062082d975...: experiments/builds/results/embsohi_93f832e_r1.SHA256SUMS, SOa/T2embryo_s2xs_ov.h5ad (built by
#       experiments/builds/embryo_so_ohi.sh, run r1); its command is that launcher's SO_ARGS
# A build step also fails (15) unless its stamp says this commit and clean True. Builds run at 16 BLAS threads, as the
# reference builds and embryo_ot23.sh's (the recorded files of O and SO were built so). Every builder and `vec check`
# output goes to /dev/null: a test build is read by its exit code, its sha256 and its stamp's commit and clean fields
# only (AGENTS.md), never its log, size, cell count or sidecar record. A build that exits non-zero means no file, and the
# reason is not read here. The flags --shape-interp-test-ok and --captain-ruling R71 are in each file's command. Nothing
# here uploads anything.
set -u
COMMIT=${1:?usage: embryo_bf.sh COMMIT RUN}
RUN=${2:?usage: embryo_bf.sh COMMIT RUN}
case "$RUN" in
  *[!A-Za-z0-9_-]*) echo "refused: RUN takes letters, digits, _ and - only" >&2; exit 90 ;;
esac
case "$COMMIT" in
  *[!0-9a-f]*) echo "refused: COMMIT must be a full hash" >&2; exit 90 ;;
esac
[ "${#COMMIT}" = 40 ] || { echo "refused: COMMIT must be a full hash" >&2; exit 90; }
export PATH=$HOME/.local/bin:$PATH
S=$HOME/v/scratch/embbf/$RUN
mkdir -p "$HOME/v/scratch/embbf"
mkdir "$S" 2> /dev/null || { echo "refused: the run folder $S exists (a run folder is never reused)" >&2; exit 91; }
mkdir "$S/tmp"
export NUMPY_MADVISE_HUGEPAGE=0 TMPDIR="$S/tmp"
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
export UV_PROJECT_ENVIRONMENT="$PWD/.venv"
PY="uv run --frozen --no-sync python -I -B"
BOARD=T2:embryo:test_interp
O_SHA=1cfc91be3a20bfdd1a4f7bcfebd7316e1143c6d4048e578c1f30bcfbe1e0e1ea  # embov_3b813b5, the overlay file O
SO_SHA=77062082d97555961901c811d8a1aaf9b5b2403417c19a6988c816ff553f1165  # embsohi_93f832e_r1, the file SO
O_SUMS=experiments/builds/results/embov_3b813b5.SHA256SUMS
O_ENTRY=OVa/T2embryo_s2x_ov.h5ad
SO_SUMS=experiments/builds/results/embsohi_93f832e_r1.SHA256SUMS
SO_ENTRY=SOa/T2embryo_s2xs_ov.h5ad
READING_DIR=experiments/coords/results/bf_gate
ROWCHECK=experiments/coords/bf_row.py
OV_ARGS="--method fgw_shape0 --positions floor_exact --expression mixture --shape-interp overlay --shape-interp-test-ok"
SO_ARGS="--method fgw_shape0 --positions floor_exact --expression series --shape-interp overlay --shape-interp-test-ok"
BFO_ARGS="--method fgw_shape0 --positions floor_exact --shape-interp overlay_bayes --shape-interp-test-ok --captain-ruling R71"
BFSO_ARGS="--method fgw_shape0 --positions floor_exact --expression series --shape-interp overlay_bayes --shape-interp-test-ok --captain-ruling R71"
BFO_FILE=T2embryo_s2x_ovbf
BFSO_FILE=T2embryo_s2xs_ovbf
STEP=gate
finish() { echo "$1" > "$S/$STEP.exit"; exit "$1"; }
step_ok() { echo 0 > "$S/$STEP.exit"; }
echo running > "$S/$STEP.exit"

head=$(git rev-parse HEAD)
[ "$head" = "$COMMIT" ] || { echo "refused: HEAD is not $COMMIT" >&2; finish 96; }
if [ -n "$(git status --porcelain --untracked-files=no)" ]; then
  echo "refused: tracked files differ from the commit (the build stamps must say clean)" >&2
  finish 98
fi
sha256sum -c docs/S2_RULES.sha256 > /dev/null 2>&1 || finish 97
python3 -I -B -c 'import sys
rows = [ln for ln in open(sys.argv[1], encoding="utf-8").read().splitlines() if ln.startswith("| R71 ")]
first = rows[0][2:].split("|", 1)[0] if len(rows) == 1 else ""
ok = (first.startswith("R71 (the captain\x27s order") and "ruling, no run" in first
      and "BF on O and on SO only" in first and "`--captain-ruling R71`" in first)
sys.exit(0 if ok else 1)' experiments/harness/TRIES.md 2> /dev/null \
  || { echo "refused: TRIES.md does not hold exactly one row R71 naming the files and the ruling" >&2; finish 95; }
python3 -I -B "$ROWCHECK" row experiments/harness/TRIES.md 2> /dev/null \
  || { echo "refused: TRIES.md does not hold exactly one well-formed row E59" >&2; finish 94; }
CELL_FILES=$(python3 -I -B "$ROWCHECK" passed experiments/harness/TRIES.md 2> /dev/null) \
  || { echo "refused: row E59's outcome cell does not say that the gate passed a file" >&2; finish 92; }
for f in reading.txt reading.json; do
  git ls-files --error-unmatch "$READING_DIR/$f" > /dev/null 2>&1 \
    || { echo "refused: the gate has no committed reading ($READING_DIR/$f must be tracked at this commit)" >&2
         finish 89; }
done
READ_OUT=$(python3 -I -B "$ROWCHECK" reading "$READING_DIR/reading.json" "$READING_DIR/reading.txt" 2> /dev/null) \
  || { echo "refused: the gate's committed reading does not say PASS on a clean tree at its run commit" >&2; finish 88; }
RUNC=$(printf '%s\n' "$READ_OUT" | head -n 1)
READ_FILES=$(printf '%s\n' "$READ_OUT" | tail -n +2)
[ "$READ_FILES" = "$CELL_FILES" ] \
  || { echo "refused: row E59's outcome cell and the committed reading pass different files" >&2; finish 88; }
git merge-base --is-ancestor "$RUNC" HEAD 2> /dev/null \
  || { echo "refused: the gate's run commit is not an ancestor of HEAD" >&2; finish 87; }
git show "$RUNC:experiments/harness/TRIES.md" > "$S/tries_at_run.md" 2> /dev/null \
  && python3 -I -B "$ROWCHECK" row "$S/tries_at_run.md" 2> /dev/null \
  || { echo "refused: row E59 is not in TRIES.md at the gate's run commit" >&2; finish 86; }
for rec in "O_SHA $O_SUMS $O_ENTRY" "SO_SHA $SO_SUMS $SO_ENTRY"; do
  read -r var sums entry <<< "$rec"
  got=$(awk -v e="$entry" '$2 == e {print $1}' "$sums" 2> /dev/null)
  [ "$got" = "${!var}" ] || { echo "refused: $sums does not give $var for $entry" >&2; finish 85; }
done
BUILD_O=0
BUILD_SO=0
while read -r f; do
  case "$f" in
    "$BFO_FILE.h5ad") BUILD_O=1 ;;
    "$BFSO_FILE.h5ad") BUILD_SO=1 ;;
    *) echo "refused: the gate names an unknown file" >&2; finish 92 ;;
  esac
done <<< "$CELL_FILES"
O=runs/embbf_${head:0:7}_$RUN
mkdir -p runs
mkdir "$O" 2> /dev/null || { echo "refused: the output folder $O exists" >&2; finish 93; }
echo "$head" > "$S/head"
echo "$RUNC" > "$S/gate_run_commit"
printf '%s\n' "$CELL_FILES" > "$S/files"
step_ok

STEP=tests
echo running > "$S/$STEP.exit"
$PY -m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp="$S/tmp/pytest" \
  tests/test_bf_build.py tests/test_bf_launchers.py tests/test_bf_form.py tests/test_bf_read.py \
  tests/test_embryo_ot23_launcher.py tests/test_mix_lean_build.py tests/test_embryo_so_ohi_launcher.py \
  tests/test_embryo_overlay.py tests/test_shape_interp.py tests/test_shape_interp_build.py tests/test_t2_build.py \
  tests/test_hidden_counts.py tests/test_s2_rules.py \
  > "$S/tests.out" 2>&1 || finish 10
step_ok

export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16
stamp() {  # the stamp's commit and clean fields only
  $PY -c 'import json,sys; s=json.load(open(sys.argv[1])); print(s["commit"], s["clean"])' "$1.stamp.json"
}
build() {  # folder name args: the build, then `vec check`; exit 11 build, 13 check, 14 no stamp, 15 stamp not this
           # commit and clean
  local d="$O/$1" name=$2 st
  shift 2
  mkdir -p "$d"
  $PY experiments/builds/build_test.py embryo "$@" --out "$d/$name.h5ad" > /dev/null 2>&1 || return 11
  uv run --frozen --no-sync vec check "$d/$name.h5ad" --board "$BOARD" > /dev/null 2>&1 || return 13
  st=$(stamp "$d/$name.h5ad" 2> /dev/null) || return 14
  echo "$(basename "$d")/$name $st" >> "$S/stamps"
  [ "$st" = "$head True" ] || return 15
}
: > "$S/stamps"

STEP=O
echo running > "$S/$STEP.exit"
build O T2embryo_s2x_ov $OV_ARGS || finish $?
got=$(sha256sum "$O/O/T2embryo_s2x_ov.h5ad" | cut -d' ' -f1)
[ "$got" = "$O_SHA" ] || { echo "stopped: the overlay command does not give O's committed file" >&2; finish 21; }
step_ok

STEP=SO
echo running > "$S/$STEP.exit"
build SO T2embryo_s2xs_ov $SO_ARGS || finish $?
got=$(sha256sum "$O/SO/T2embryo_s2xs_ov.h5ad" | cut -d' ' -f1)
[ "$got" = "$SO_SHA" ] || { echo "stopped: the SO command does not give SO's committed file" >&2; finish 22; }
step_ok

LIST="O/*.h5ad SO/*.h5ad"
if [ "$BUILD_O" = 1 ]; then
  for STEP in BFOa BFOb; do
    echo running > "$S/$STEP.exit"
    build "$STEP" "$BFO_FILE" $BFO_ARGS || finish $?
    step_ok
  done
  LIST="$LIST BFOa/*.h5ad BFOb/*.h5ad"
fi
if [ "$BUILD_SO" = 1 ]; then
  for STEP in BFSOa BFSOb; do
    echo running > "$S/$STEP.exit"
    build "$STEP" "$BFSO_FILE" $BFSO_ARGS || finish $?
    step_ok
  done
  LIST="$LIST BFSOa/*.h5ad BFSOb/*.h5ad"
fi

STEP=done
echo running > "$S/$STEP.exit"
(cd "$O" && sha256sum $LIST) > "$S/SHA256SUMS" || finish 30
cp "$S/SHA256SUMS" "$O/SHA256SUMS" || finish 30
if [ "$BUILD_O" = 1 ]; then
  cmp -s "$O/BFOa/$BFO_FILE.h5ad" "$O/BFOb/$BFO_FILE.h5ad" || finish 31
fi
if [ "$BUILD_SO" = 1 ]; then
  cmp -s "$O/BFSOa/$BFSO_FILE.h5ad" "$O/BFSOb/$BFSO_FILE.h5ad" || finish 31
fi
step_ok
