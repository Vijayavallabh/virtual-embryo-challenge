#!/usr/bin/env bash
# Row R59's file T2embryo_s2x_ovt23.h5ad (experiments/harness/TRIES.md, the captain's order of 2026-11-24: "Yes to
# T2embryo_s2x_ovt23.h5ad (Ot23): my named exception on the embryo's overlay form for this file, bound to its gate G1."):
# O's command with --mix-lean later --captain-ruling R59 (t_mix = 2/3), built twice at one clean commit, and only when
# its gate G1 (row E58, experiments/coords/e58_mix_lean_read.py) passed: row E58's outcome cell, read by code, holds the
# read's PASS words, and the committed reading says PASS. Test-board builds on a DGX campaign checkout (~/v/vec-campaign,
# never the reference checkout). Run from the root of the checkout or detached worktree it is called in, detached, and
# wait on the exit files with a background poller; never a foreground sleep:
#   bash experiments/builds/embryo_ot23.sh <COMMIT> <RUN>
# COMMIT is the full hash of the build commit. RUN names this run's own scratch folder ~/v/scratch/embot23/<RUN> and its
# output folder runs/embot23_<commit7>_<RUN>. Neither may exist: a run folder is never reused and nothing is ever deleted
# here, so a second run takes a new name. In order, each step writing ~/v/scratch/embot23/<RUN>/<STEP>.exit (the word
# running first, then its exit code):
#   gate     HEAD is COMMIT (else 96), the tracked tree is clean (98), docs/S2_RULES.sha256 checks (97), TRIES.md holds
#            exactly one row R59 whose first cell is the captain's order ("ruling, no run") naming `--mix-lean later`
#            and `--captain-ruling R59`, the builder's own check (95), exactly one well-formed row E58 (94) whose outcome
#            cell holds the read's PASS words and no STOP or VOID words (92; experiments/coords/e58_row.py passed), and
#            G1's reading is committed (89: READING_DIR's reading.txt and reading.json tracked at this commit) and says
#            PASS on a clean tree at its run commit (88), an ancestor of HEAD (87) whose TRIES.md holds a well-formed
#            row E58 (86). The identity pins of O and W are read from their committed SHA256SUMS lists (85; never from a
#            .txt record) and must equal the ones written here
#   tests    the builder's, the gate's and this launcher's targeted tests (synthetic stages only), 8 BLAS threads
#   O        the identity file O by its committed command; it must give O_SHA or everything stops (21)
#   W        the identity file W (row R48) by its committed command; it must give W_SHA or everything stops (22)
#   OTa OTb  T2embryo_s2x_ovt23.h5ad by the file's command, twice, into separate folders
#   done     one sha256 list over every built file; OTa and OTb compared byte for byte (31)
# Where the identity pins and their commands come from:
#   O  1cfc91be3a20...: experiments/builds/results/embov_3b813b5.SHA256SUMS, OVa/T2embryo_s2x_ov.h5ad (built by
#      experiments/builds/embryo_overlay.sh at 3b813b5); its command is that launcher's OV_ARGS
#   W  f576b64b8b93...: experiments/builds/results/embovsb_20b2381_w1.SHA256SUMS, Wa/T2embryo_s2x_ovsb.h5ad (built by
#      experiments/builds/embryo_ovsb.sh, run w1); its command is that launcher's W_ARGS
# A build step also fails (15) unless its stamp says this commit and clean True. Builds run at 16 BLAS threads, as the
# reference builds and embryo_ovsb.sh's (the recorded files of O and W were built so). Every builder and `vec check`
# output goes to /dev/null: a test build is read by its exit code, its sha256 and its stamp's commit and clean fields
# only (AGENTS.md), never its log, size, cell count or sidecar record. A build that exits non-zero means no file, and the
# reason is not read here. The flags --shape-interp-test-ok, --mix-lean later and --captain-ruling R59 are in the
# file's command. Nothing here uploads anything.
set -u
COMMIT=${1:?usage: embryo_ot23.sh COMMIT RUN}
RUN=${2:?usage: embryo_ot23.sh COMMIT RUN}
case "$RUN" in
  *[!A-Za-z0-9_-]*) echo "refused: RUN takes letters, digits, _ and - only" >&2; exit 90 ;;
esac
case "$COMMIT" in
  *[!0-9a-f]*) echo "refused: COMMIT must be a full hash" >&2; exit 90 ;;
esac
[ "${#COMMIT}" = 40 ] || { echo "refused: COMMIT must be a full hash" >&2; exit 90; }
export PATH=$HOME/.local/bin:$PATH
S=$HOME/v/scratch/embot23/$RUN
mkdir -p "$HOME/v/scratch/embot23"
mkdir "$S" 2> /dev/null || { echo "refused: the run folder $S exists (a run folder is never reused)" >&2; exit 91; }
mkdir "$S/tmp"
export NUMPY_MADVISE_HUGEPAGE=0 TMPDIR="$S/tmp"
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
export UV_PROJECT_ENVIRONMENT="$PWD/.venv"
PY="uv run --frozen --no-sync python -I -B"
BOARD=T2:embryo:test_interp
O_SHA=1cfc91be3a20bfdd1a4f7bcfebd7316e1143c6d4048e578c1f30bcfbe1e0e1ea  # embov_3b813b5, the overlay file O
W_SHA=f576b64b8b93d79b1269b580df31bd5c6e58246cd4f3eda12cba7d77387087ab  # embovsb_20b2381_w1, the file W (row R48)
O_SUMS=experiments/builds/results/embov_3b813b5.SHA256SUMS
O_ENTRY=OVa/T2embryo_s2x_ov.h5ad
W_SUMS=experiments/builds/results/embovsb_20b2381_w1.SHA256SUMS
W_ENTRY=Wa/T2embryo_s2x_ovsb.h5ad
READING_DIR=experiments/coords/results/e58_mix_lean
ROWCHECK=experiments/coords/e58_row.py
OV_ARGS="--method fgw_shape0 --positions floor_exact --expression mixture --shape-interp overlay --shape-interp-test-ok"
W_ARGS="--method fgw_shape0 --positions floor_exact --expression series_back --shape-interp overlay --shape-interp-test-ok --captain-ruling R48"
OT_ARGS="--method fgw_shape0 --positions floor_exact --expression mixture --shape-interp overlay --shape-interp-test-ok --mix-lean later --captain-ruling R59"
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
rows = [ln for ln in open(sys.argv[1], encoding="utf-8").read().splitlines() if ln.startswith("| R59 ")]
first = rows[0][2:].split("|", 1)[0] if len(rows) == 1 else ""
ok = (first.startswith("R59 (the captain\x27s order") and "ruling, no run" in first
      and "`--mix-lean later`" in first and "`--captain-ruling R59`" in first)
sys.exit(0 if ok else 1)' experiments/harness/TRIES.md 2> /dev/null \
  || { echo "refused: TRIES.md does not hold exactly one row R59 naming the option and the ruling" >&2; finish 95; }
python3 -I -B "$ROWCHECK" row experiments/harness/TRIES.md 2> /dev/null \
  || { echo "refused: TRIES.md does not hold exactly one well-formed row E58" >&2; finish 94; }
python3 -I -B "$ROWCHECK" passed experiments/harness/TRIES.md 2> /dev/null \
  || { echo "refused: row E58's outcome cell does not say that G1 passed" >&2; finish 92; }
for f in reading.txt reading.json; do
  git ls-files --error-unmatch "$READING_DIR/$f" > /dev/null 2>&1 \
    || { echo "refused: G1 has no committed reading ($READING_DIR/$f must be tracked at this commit)" >&2
         finish 89; }
done
RUNC=$(python3 -I -B "$ROWCHECK" reading "$READING_DIR/reading.json" "$READING_DIR/reading.txt" 2> /dev/null) \
  || { echo "refused: G1's committed reading does not say PASS on a clean tree at its run commit" >&2; finish 88; }
git merge-base --is-ancestor "$RUNC" HEAD 2> /dev/null \
  || { echo "refused: G1's run commit is not an ancestor of HEAD" >&2; finish 87; }
git show "$RUNC:experiments/harness/TRIES.md" > "$S/tries_at_run.md" 2> /dev/null \
  && python3 -I -B "$ROWCHECK" row "$S/tries_at_run.md" 2> /dev/null \
  || { echo "refused: row E58 is not in TRIES.md at G1's run commit" >&2; finish 86; }
for rec in "O_SHA $O_SUMS $O_ENTRY" "W_SHA $W_SUMS $W_ENTRY"; do
  read -r var sums entry <<< "$rec"
  got=$(awk -v e="$entry" '$2 == e {print $1}' "$sums" 2> /dev/null)
  [ "$got" = "${!var}" ] || { echo "refused: $sums does not give $var for $entry" >&2; finish 85; }
done
O=runs/embot23_${head:0:7}_$RUN
mkdir -p runs
mkdir "$O" 2> /dev/null || { echo "refused: the output folder $O exists" >&2; finish 93; }
echo "$head" > "$S/head"
echo "$RUNC" > "$S/g1_run_commit"
step_ok

STEP=tests
echo running > "$S/$STEP.exit"
$PY -m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp="$S/tmp/pytest" \
  tests/test_embryo_ot23_launcher.py tests/test_mix_lean_build.py tests/test_e58_mix_lean_read.py \
  tests/test_series_back_build.py tests/test_embryo_ovsb_launcher.py tests/test_embryo_overlay.py \
  tests/test_shape_interp.py tests/test_shape_interp_build.py tests/test_e47_floor_exact.py tests/test_t2_build.py \
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

STEP=W
echo running > "$S/$STEP.exit"
build W T2embryo_s2x_ovsb $W_ARGS || finish $?
got=$(sha256sum "$O/W/T2embryo_s2x_ovsb.h5ad" | cut -d' ' -f1)
[ "$got" = "$W_SHA" ] || { echo "stopped: the series_back command does not give W's committed file" >&2; finish 22; }
step_ok

for STEP in OTa OTb; do
  echo running > "$S/$STEP.exit"
  build "$STEP" T2embryo_s2x_ovt23 $OT_ARGS || finish $?
  step_ok
done

STEP=done
echo running > "$S/$STEP.exit"
(cd "$O" && sha256sum O/*.h5ad W/*.h5ad OTa/*.h5ad OTb/*.h5ad) > "$S/SHA256SUMS" || finish 30
cp "$S/SHA256SUMS" "$O/SHA256SUMS" || finish 30
cmp -s "$O/OTa/T2embryo_s2x_ovt23.h5ad" "$O/OTb/T2embryo_s2x_ovt23.h5ad" || finish 31
step_ok
