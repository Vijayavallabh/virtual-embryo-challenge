#!/usr/bin/env bash
# Row E57's file T2embryo_s2x_ovrs.h5ad ("RS: one run per side"; experiments/harness/TRIES.md): O's command with
# --neighbourhood run_sorted --captain-ruling <ROW>, an extra file beside O and never in its place, built twice at one
# clean commit, and only after row E57's read passed and the captain's ruling row ROW names the option. Test-board
# builds on a DGX campaign checkout (~/v/vec-campaign, never the reference checkout). Run from the root of the checkout
# or detached worktree it is called in, detached, and wait on the exit files with a background poller; never a
# foreground sleep:
#   bash experiments/builds/embryo_rs.sh <COMMIT> <RUN> <ROW>
# COMMIT is the full hash of the build commit; ROW the captain's ruling row (R, digits, an optional letter). RUN names
# this run's own scratch folder ~/v/scratch/embrs/<RUN> and its output folder runs/embrs_<commit7>_<RUN>. Neither may
# exist: a run folder is never reused and nothing is ever deleted here, so a second run takes a new name. In order,
# each step writing ~/v/scratch/embrs/<RUN>/<STEP>.exit (the word running first, then its exit code):
#   gate     HEAD is COMMIT (96), the tracked tree is clean (98), docs/S2_RULES.sha256 checks (97), TRIES.md holds
#            exactly one row ROW whose first cell is the captain's ruling ("ROW (the captain's ...", "ruling, no run")
#            naming `--neighbourhood run_sorted` and `--captain-ruling ROW` (95, the builder's own check), and row
#            E57's reading is committed: READING_DIR's reading.txt and reading.json are tracked files of this commit
#            (94), the text's first line is exactly "E57: PASS" and so is the JSON's outcome, not void, made on a clean
#            tree at its run commit (92), an ancestor of HEAD (88) whose TRIES.md holds row E57 (86). The committed
#            sha256 of O and of row R48's file are found in their build records (below) and must equal the ones written
#            here (85); the output folder must not exist (93)
#   tests    the builder's, the step's, the read's and this launcher's targeted tests (synthetic stages only)
#   O        the overlay file O by its committed command; it must give O_SHA or everything stops (21)
#   W        row R48's file by its committed command, the other file this builder takes with --captain-ruling; it
#            must give W_SHA or everything stops (22)
#   RSa RSb  T2embryo_s2x_ovrs.h5ad by the file's command, twice, into separate folders
#   identity booleans only, to identity.out: RS's obsm coordinates equal O's bit for bit, its obs_names and var_names
#            equal O's, its expression rows sorted equal O's sorted, and its X is not O's (O as built in step O)
#   done     one sha256 list over every built file; RSa and RSb compared byte for byte (31)
# Where the two committed files and their commands come from:
#   O  1cfc91be3a20...: experiments/builds/results/embov_3b813b5.SHA256SUMS, OVa/T2embryo_s2x_ov.h5ad (built by
#      experiments/builds/embryo_overlay.sh at 3b813b5); its command is that launcher's OV_ARGS
#   W  f576b64b8b93...: experiments/builds/results/embovsb_20b2381_w1.SHA256SUMS, Wa/T2embryo_s2x_ovsb.h5ad (built by
#      experiments/builds/embryo_ovsb.sh, run w1); its command is that launcher's W_ARGS
# A build step also fails (15) unless its stamp says this commit and clean True. Builds run at 16 BLAS threads, as the
# recorded files were built. Every builder and `vec check` output goes to /dev/null: a test build is read by its exit
# code, its sha256 and its stamp's commit and clean fields only (AGENTS.md), never its log, size, cell count or sidecar
# record. A build that exits non-zero means no file, and the reason is not read here. Nothing here uploads anything.
set -u
COMMIT=${1:?usage: embryo_rs.sh COMMIT RUN ROW}
RUN=${2:?usage: embryo_rs.sh COMMIT RUN ROW}
ROW=${3:?usage: embryo_rs.sh COMMIT RUN ROW}
case "$RUN" in
  *[!A-Za-z0-9_-]*) echo "refused: RUN takes letters, digits, _ and - only" >&2; exit 90 ;;
esac
case "$COMMIT" in
  *[!0-9a-f]*) echo "refused: COMMIT must be a full hash" >&2; exit 90 ;;
esac
[ "${#COMMIT}" = 40 ] || { echo "refused: COMMIT must be a full hash" >&2; exit 90; }
[[ "$ROW" =~ ^R[0-9]+[a-z]?$ ]] || { echo "refused: ROW must be a ruling row's ID (R, digits, a letter)" >&2; exit 90; }
export PATH=$HOME/.local/bin:$PATH
S=$HOME/v/scratch/embrs/$RUN
mkdir -p "$HOME/v/scratch/embrs"
mkdir "$S" 2> /dev/null || { echo "refused: the run folder $S exists (a run folder is never reused)" >&2; exit 91; }
mkdir "$S/tmp"
export NUMPY_MADVISE_HUGEPAGE=0 TMPDIR="$S/tmp"
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
export UV_PROJECT_ENVIRONMENT="$PWD/.venv"
PY="uv run --frozen --no-sync python -I -B"
BOARD=T2:embryo:test_interp
O_SHA=1cfc91be3a20bfdd1a4f7bcfebd7316e1143c6d4048e578c1f30bcfbe1e0e1ea  # embov_3b813b5, the overlay file O
W_SHA=f576b64b8b93d79b1269b580df31bd5c6e58246cd4f3eda12cba7d77387087ab  # embovsb_20b2381_w1, row R48's file
O_SUMS=experiments/builds/results/embov_3b813b5.SHA256SUMS
O_ENTRY=OVa/T2embryo_s2x_ov.h5ad
W_SUMS=experiments/builds/results/embovsb_20b2381_w1.SHA256SUMS
W_ENTRY=Wa/T2embryo_s2x_ovsb.h5ad
READING_DIR=experiments/coords/results/e57_rs
OV_ARGS="--method fgw_shape0 --positions floor_exact --expression mixture --shape-interp overlay --shape-interp-test-ok"
W_ARGS="--method fgw_shape0 --positions floor_exact --expression series_back --shape-interp overlay --shape-interp-test-ok --captain-ruling R48"
RS_ARGS="$OV_ARGS --neighbourhood run_sorted --captain-ruling $ROW"
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
row = sys.argv[2]
rows = [ln for ln in open(sys.argv[1], encoding="utf-8").read().splitlines() if ln.startswith(f"| {row} ")]
first = rows[0][2:].split("|", 1)[0] if len(rows) == 1 else ""
ok = (first.startswith(f"{row} (the captain\x27s ") and "ruling, no run" in first
      and "`--neighbourhood run_sorted`" in first and f"`--captain-ruling {row}`" in first)
sys.exit(0 if ok else 1)' experiments/harness/TRIES.md "$ROW" 2> /dev/null \
  || { echo "refused: TRIES.md does not hold exactly one row $ROW naming the option and the ruling" >&2; finish 95; }
for f in reading.txt reading.json; do
  git ls-files --error-unmatch "$READING_DIR/$f" > /dev/null 2>&1 \
    || { echo "refused: row E57 has no committed reading ($READING_DIR/$f must be tracked at this commit)" >&2
         finish 94; }
done
RUNC=$(python3 -I -B -c 'import json, re, sys
first = open(sys.argv[1], encoding="utf-8").read().split("\n", 1)[0]
r = json.load(open(sys.argv[2], encoding="utf-8"))
p = r.get("reading_code") if isinstance(r.get("reading_code"), dict) else {}
c = str(r.get("commit"))
if (first != "E57: PASS" or r.get("outcome") != "E57: PASS" or r.get("void") is not False
        or p.get("clean") is not True or p.get("commit") != c or not re.fullmatch("[0-9a-f]{40}", c)):
    sys.exit(1)
print(c)' "$READING_DIR/reading.txt" "$READING_DIR/reading.json" 2> /dev/null) \
  || { echo "refused: row E57's read did not pass (its committed reading is not E57: PASS on a clean tree)" >&2
       finish 92; }
git merge-base --is-ancestor "$RUNC" HEAD 2> /dev/null \
  || { echo "refused: row E57's run commit is not an ancestor of HEAD" >&2; finish 88; }
git show "$RUNC:experiments/harness/TRIES.md" 2> /dev/null | grep -q '^| E57 ' \
  || { echo "refused: row E57 is not in TRIES.md at its run commit" >&2; finish 86; }
for rec in "O_SHA $O_SUMS $O_ENTRY" "W_SHA $W_SUMS $W_ENTRY"; do
  read -r var sums entry <<< "$rec"
  got=$(awk -v e="$entry" '$2 == e {print $1}' "$sums" 2> /dev/null)
  [ "$got" = "${!var}" ] || { echo "refused: $sums does not give $var for $entry" >&2; finish 85; }
done
O=runs/embrs_${head:0:7}_$RUN
mkdir -p runs
mkdir "$O" 2> /dev/null || { echo "refused: the output folder $O exists" >&2; finish 93; }
echo "$head" > "$S/head"
echo "$RUNC" > "$S/e57_run_commit"
echo "$ROW" > "$S/ruling_row"
step_ok

STEP=tests
echo running > "$S/$STEP.exit"
$PY -m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp="$S/tmp/pytest" \
  tests/test_embryo_rs_launcher.py tests/test_run_sorted.py tests/test_e57_rs_read.py \
  tests/test_series_back_build.py tests/test_embryo_overlay.py tests/test_shape_interp.py \
  tests/test_shape_interp_build.py tests/test_e47_floor_exact.py tests/test_t2_build.py tests/test_hidden_counts.py \
  tests/test_s2_rules.py > "$S/tests.out" 2>&1 || finish 10
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
[ "$got" = "$W_SHA" ] || { echo "stopped: row R48's command does not give its committed file" >&2; finish 22; }
step_ok

for STEP in RSa RSb; do
  echo running > "$S/$STEP.exit"
  build "$STEP" T2embryo_s2x_ovrs $RS_ARGS || finish $?
  step_ok
done

STEP=identity
echo running > "$S/$STEP.exit"
$PY -c '
import sys
import anndata as ad
import numpy as np
o_path, w_path = sys.argv[1:3]
def x(a):
    return np.ascontiguousarray(a.X.toarray() if hasattr(a.X, "toarray") else np.asarray(a.X))
def same(p, q):
    p, q = np.ascontiguousarray(np.asarray(p)), np.ascontiguousarray(np.asarray(q))
    return p.dtype == q.dtype and p.shape == q.shape and p.tobytes() == q.tobytes()
def rows(m):
    return m[np.lexsort(m.T[::-1])] if m.size else m
o, w = ad.read_h5ad(o_path), ad.read_h5ad(w_path)
xo, xw = x(o), x(w)
res = {"rs_coords_equal_o": same(w.obsm["spatial_3D"], o.obsm["spatial_3D"]),
       "rs_obs_names_equal_o": list(w.obs_names) == list(o.obs_names),
       "rs_var_names_equal_o": list(w.var_names) == list(o.var_names),
       "rs_sorted_rows_equal_o": same(rows(xw), rows(xo)),
       "rs_x_not_o": not same(xw, xo)}
for k, v in res.items():
    print(k, bool(v))
sys.exit(0 if all(res.values()) else 1)
' "$O/O/T2embryo_s2x_ov.h5ad" "$O/RSa/T2embryo_s2x_ovrs.h5ad" > "$S/identity.out" 2> /dev/null || finish 40
step_ok

STEP=done
echo running > "$S/$STEP.exit"
(cd "$O" && sha256sum O/*.h5ad W/*.h5ad RSa/*.h5ad RSb/*.h5ad) > "$S/SHA256SUMS" || finish 30
cp "$S/SHA256SUMS" "$O/SHA256SUMS" || finish 30
cmp -s "$O/RSa/T2embryo_s2x_ovrs.h5ad" "$O/RSb/T2embryo_s2x_ovrs.h5ad" || finish 31
step_ok
