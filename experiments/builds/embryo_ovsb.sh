#!/usr/bin/env bash
# Row R48's file T2embryo_s2x_ovsb.h5ad (experiments/harness/TRIES.md, the captain's order of 2026-11-24, "Build
# T2embryo_s2x_ovsb.h5ad as named": an extra file beside O, never in its place): O's command with --expression
# series_back --captain-ruling R48, built twice at one clean commit, and only after its step 0 (row E56,
# experiments/coords/ovsb_step0.py) passed. Test-board builds on a DGX campaign checkout (~/v/vec-campaign, never the
# reference checkout). Run from the root of the checkout or detached worktree it is called in, detached, and wait on
# the exit files with a background poller; never a foreground sleep:
#   bash experiments/builds/embryo_ovsb.sh <COMMIT> <RUN>
# COMMIT is the full hash of the build commit. RUN names this run's own scratch folder ~/v/scratch/embovsb/<RUN> and its
# output folder runs/embovsb_<commit7>_<RUN>. Neither may exist: a run folder is never reused and nothing is ever
# deleted here, so a second run takes a new name. In order, each step writing ~/v/scratch/embovsb/<RUN>/<STEP>.exit
# (the word running first, then its exit code):
#   gate     HEAD is COMMIT, the tracked tree is clean, docs/S2_RULES.sha256 checks, TRIES.md holds exactly one row R48
#            whose first cell is the captain's order ("ruling, no run") naming `--expression series_back` and
#            `--captain-ruling R48` (the builder's own check), and step 0's reading is committed: READING_DIR's
#            reading.txt and reading.json are tracked files of this commit, the text's first line is exactly
#            "Step 0: PASS" and so is the JSON's outcome, made on a clean tree at a commit that is an ancestor of HEAD and
#            whose TRIES.md holds row E56. The committed sha256 of O, SO and O_hi are found in their build records
#            (below) and must equal the ones written here
#   tests    the builder's, step 0's and this launcher's targeted tests (synthetic stages only), 8 BLAS threads
#   O        the overlay file O by its committed command; it must give O_SHA or everything stops
#   SO       the series overlay file SO by its committed command; it must give SO_SHA or everything stops
#   OH       the topped-up overlay file O_hi by its committed command; it must give OH_SHA or everything stops
#   Wa Wb    T2embryo_s2x_ovsb.h5ad by the file's command, twice, into separate folders
#   identity booleans only, to identity.out: W's obsm coordinates equal O's bit for bit, its obs_names and var_names
#            equal O's, and its X is not O's (O as built in step O of this same run)
#   done     one sha256 list over every built file; Wa and Wb compared byte for byte
# Where the three committed files and their commands come from:
#   O     1cfc91be3a20...: experiments/builds/results/embov_3b813b5.SHA256SUMS, OVa/T2embryo_s2x_ov.h5ad (built by
#         experiments/builds/embryo_overlay.sh at 3b813b5); its command is that launcher's OV_ARGS
#   SO    77062082d975...: experiments/builds/results/embsohi_93f832e_r1.SHA256SUMS, SOa/T2embryo_s2xs_ov.h5ad (built
#         by experiments/builds/embryo_so_ohi.sh, run r1); its command is that launcher's SO_ARGS
#   O_hi  7e4d1eaafc40...: the same list, OHa/T2embryo_s2x_ovhi.h5ad; its command is that launcher's OH_ARGS
# A build step also fails (15) unless its stamp says this commit and clean True. Builds run at 16 BLAS threads, as the
# reference builds and embryo_so_ohi.sh's (the recorded files of O, SO and O_hi were built so). Every builder and
# `vec check` output goes to /dev/null: a test build is read by its exit code, its sha256 and its stamp's commit and
# clean fields only (AGENTS.md), never its log, size, cell count or sidecar record. A build that exits non-zero means no
# file, and the reason is not read here. The flags --shape-interp-test-ok and --captain-ruling R48 are in the file's
# command. Nothing here uploads anything.
set -u
COMMIT=${1:?usage: embryo_ovsb.sh COMMIT RUN}
RUN=${2:?usage: embryo_ovsb.sh COMMIT RUN}
case "$RUN" in
  *[!A-Za-z0-9_-]*) echo "refused: RUN takes letters, digits, _ and - only" >&2; exit 90 ;;
esac
case "$COMMIT" in
  *[!0-9a-f]*) echo "refused: COMMIT must be a full hash" >&2; exit 90 ;;
esac
[ "${#COMMIT}" = 40 ] || { echo "refused: COMMIT must be a full hash" >&2; exit 90; }
export PATH=$HOME/.local/bin:$PATH
S=$HOME/v/scratch/embovsb/$RUN
mkdir -p "$HOME/v/scratch/embovsb"
mkdir "$S" 2> /dev/null || { echo "refused: the run folder $S exists (a run folder is never reused)" >&2; exit 91; }
mkdir "$S/tmp"
export NUMPY_MADVISE_HUGEPAGE=0 TMPDIR="$S/tmp"
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
export UV_PROJECT_ENVIRONMENT="$PWD/.venv"
PY="uv run --frozen --no-sync python -I -B"
BOARD=T2:embryo:test_interp
O_SHA=1cfc91be3a20bfdd1a4f7bcfebd7316e1143c6d4048e578c1f30bcfbe1e0e1ea  # embov_3b813b5, the overlay file O
SO_SHA=77062082d97555961901c811d8a1aaf9b5b2403417c19a6988c816ff553f1165  # embsohi_93f832e_r1, SO
OH_SHA=7e4d1eaafc40080d9d9c461d4b9396f0cee6353686538e79fb92ca31d1d2deeb  # embsohi_93f832e_r1, O_hi
O_SUMS=experiments/builds/results/embov_3b813b5.SHA256SUMS
O_ENTRY=OVa/T2embryo_s2x_ov.h5ad
SO_SUMS=experiments/builds/results/embsohi_93f832e_r1.SHA256SUMS
SO_ENTRY=SOa/T2embryo_s2xs_ov.h5ad
OH_SUMS=experiments/builds/results/embsohi_93f832e_r1.SHA256SUMS
OH_ENTRY=OHa/T2embryo_s2x_ovhi.h5ad
READING_DIR=experiments/coords/results/ovsb_step0
OV_ARGS="--method fgw_shape0 --positions floor_exact --expression mixture --shape-interp overlay --shape-interp-test-ok"
SO_ARGS="--method fgw_shape0 --positions floor_exact --expression series --shape-interp overlay --shape-interp-test-ok"
OH_ARGS="--method fgw_shape0 --positions floor_exact --expression render_hi --shape-interp overlay --shape-interp-test-ok"
W_ARGS="--method fgw_shape0 --positions floor_exact --expression series_back --shape-interp overlay --shape-interp-test-ok --captain-ruling R48"
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
rows = [ln for ln in open(sys.argv[1], encoding="utf-8").read().splitlines() if ln.startswith("| R48 ")]
first = rows[0][2:].split("|", 1)[0] if len(rows) == 1 else ""
ok = (first.startswith("R48 (the captain\x27s order") and "ruling, no run" in first
      and "`--expression series_back`" in first and "`--captain-ruling R48`" in first)
sys.exit(0 if ok else 1)' experiments/harness/TRIES.md 2> /dev/null \
  || { echo "refused: TRIES.md does not hold exactly one row R48 naming the option and the ruling" >&2; finish 95; }
for f in reading.txt reading.json; do
  git ls-files --error-unmatch "$READING_DIR/$f" > /dev/null 2>&1 \
    || { echo "refused: step 0 has no committed reading ($READING_DIR/$f must be tracked at this commit)" >&2
         finish 94; }
done
RUNC=$(python3 -I -B -c 'import json, re, sys
first = open(sys.argv[1], encoding="utf-8").read().split("\n", 1)[0]
r = json.load(open(sys.argv[2], encoding="utf-8"))
p = r.get("provenance") if isinstance(r, dict) and isinstance(r.get("provenance"), dict) else {}
c = str(p.get("commit"))
if (first != "Step 0: PASS" or r.get("outcome") != "Step 0: PASS" or p.get("clean") is not True
        or not re.fullmatch("[0-9a-f]{40}", c)):
    sys.exit(1)
print(c)' "$READING_DIR/reading.txt" "$READING_DIR/reading.json" 2> /dev/null) \
  || { echo "refused: step 0 did not pass (its committed reading is not Step 0: PASS on a clean tree)" >&2; finish 92; }
git merge-base --is-ancestor "$RUNC" HEAD 2> /dev/null \
  || { echo "refused: step 0's commit is not an ancestor of HEAD" >&2; finish 88; }
git show "$RUNC:experiments/harness/TRIES.md" 2> /dev/null | grep -q '^| E56 ' \
  || { echo "refused: row E56 is not in TRIES.md at step 0's commit" >&2; finish 86; }
for rec in "O_SHA $O_SUMS $O_ENTRY" "SO_SHA $SO_SUMS $SO_ENTRY" "OH_SHA $OH_SUMS $OH_ENTRY"; do
  read -r var sums entry <<< "$rec"
  got=$(awk -v e="$entry" '$2 == e {print $1}' "$sums" 2> /dev/null)
  [ "$got" = "${!var}" ] || { echo "refused: $sums does not give $var for $entry" >&2; finish 85; }
done
O=runs/embovsb_${head:0:7}_$RUN
mkdir -p runs
mkdir "$O" 2> /dev/null || { echo "refused: the output folder $O exists" >&2; finish 93; }
echo "$head" > "$S/head"
echo "$RUNC" > "$S/step0_commit"
step_ok

STEP=tests
echo running > "$S/$STEP.exit"
$PY -m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp="$S/tmp/pytest" \
  tests/test_embryo_ovsb_launcher.py tests/test_series_back_build.py tests/test_ovsb_step0.py \
  tests/test_embryo_overlay.py tests/test_shape_interp.py tests/test_shape_interp_build.py \
  tests/test_e47_floor_exact.py tests/test_t2_build.py tests/test_hidden_counts.py tests/test_s2_rules.py \
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
[ "$got" = "$SO_SHA" ] || { echo "stopped: the series overlay command does not give SO's committed file" >&2; finish 22; }
step_ok

STEP=OH
echo running > "$S/$STEP.exit"
build OH T2embryo_s2x_ovhi $OH_ARGS || finish $?
got=$(sha256sum "$O/OH/T2embryo_s2x_ovhi.h5ad" | cut -d' ' -f1)
[ "$got" = "$OH_SHA" ] || { echo "stopped: the top-up overlay command does not give O_hi's committed file" >&2; finish 23; }
step_ok

for STEP in Wa Wb; do
  echo running > "$S/$STEP.exit"
  build "$STEP" T2embryo_s2x_ovsb $W_ARGS || finish $?
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
o, w = ad.read_h5ad(o_path), ad.read_h5ad(w_path)
res = {"w_coords_equal_o": same(w.obsm["spatial_3D"], o.obsm["spatial_3D"]),
       "w_obs_names_equal_o": list(w.obs_names) == list(o.obs_names),
       "w_var_names_equal_o": list(w.var_names) == list(o.var_names),
       "w_x_not_o": not same(x(w), x(o))}
for k, v in res.items():
    print(k, bool(v))
sys.exit(0 if all(res.values()) else 1)
' "$O/O/T2embryo_s2x_ov.h5ad" "$O/Wa/T2embryo_s2x_ovsb.h5ad" > "$S/identity.out" 2> /dev/null || finish 40
step_ok

STEP=done
echo running > "$S/$STEP.exit"
(cd "$O" && sha256sum O/*.h5ad SO/*.h5ad OH/*.h5ad Wa/*.h5ad Wb/*.h5ad) > "$S/SHA256SUMS" || finish 30
cp "$S/SHA256SUMS" "$O/SHA256SUMS" || finish 30
cmp -s "$O/Wa/T2embryo_s2x_ovsb.h5ad" "$O/Wb/T2embryo_s2x_ovsb.h5ad" || finish 31
step_ok
