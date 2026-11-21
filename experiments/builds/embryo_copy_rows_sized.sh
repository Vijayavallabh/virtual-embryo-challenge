#!/usr/bin/env bash
# The embryo's copy-rows file at the team's pre-registered size, C_s (practice row E52 of experiments/harness/TRIES.md;
# the captain's order of 2026-11-21, docs/PLAN.md decision log 2026-11-21 (v)): test-board builds on a DGX campaign
# checkout (~/v/vec-campaign, never the reference checkout). Cut from embryo_copy_rows.sh. Run from the checkout's root,
# detached, and wait on the exit files with a background poller; never a foreground sleep:
#   bash experiments/builds/embryo_copy_rows_sized.sh <COMMIT> <RUN> <E52_RUN_COMMIT>
# COMMIT is the build commit and E52_RUN_COMMIT the run commit of row E52's reading (both full hashes). RUN names this
# run's own scratch folder ~/v/scratch/embcrs/<RUN> and its output folder runs/embcrs_<commit7>_<RUN>. Neither may
# exist: a run folder is never reused and nothing is ever deleted here, so a second run takes a new name.
# In order, each step writing ~/v/scratch/embcrs/<RUN>/<STEP>.exit:
#   gate     HEAD is COMMIT, the tracked tree is clean, docs/S2_RULES.sha256 checks, row E52 is in TRIES.md here and at
#            E52_RUN_COMMIT (an ancestor of HEAD), row E52's reading is a tracked file of this commit, made at
#            E52_RUN_COMMIT, not void, of the arm fgw_copyrows_sized, and says the test file is built (its reading (i)
#            holds), and the code the builder runs (src/, the builder and the modules it imports, the count and size
#            records, boards/, configs/, the lock file) is byte for byte the code that reading ran
#   tests    the construction's and the embryo builder's targeted tests (synthetic stages only), 8 BLAS threads
#   D        the positions file P by its command; it must give POSITIONS_SHA or everything stops
#   CR       the copy-rows file C by its command; it must give C's committed sha256 (CR_SUMS) or everything stops
#   CSa CSb  the sized copy-rows file C_s, twice, into separate folders
#   identity booleans only: C_s's X, obs_names and var_names equal C's, and its coordinates equal
#            float32(rescale_to_rms(float64 of C's coordinates, target)) bit for bit, the target read by code from the
#            embryo size JSON (valid, raw mode, preregistered.pred_rms_full)
#   done     one sha256 list over every built file, CSa and CSb compared byte for byte
# A build step also fails (15) unless its stamp says this commit and clean True.
# Builds run at 16 BLAS threads, as the reference builds. Every builder and `vec check` output goes to /dev/null: a
# test build is read by its exit code, its sha256 and its stamp's commit and clean fields only (AGENTS.md), never its
# log, size, cell count or sidecar record. A build that exits non-zero (step exit 11) means no file, and the reason is
# not read. The flags --copy-rows-test-ok and --copy-rows-size-test-ok are the captain's orders of 2026-11-21 (decision
# log 2026-11-21 (h) and (v)). Nothing here uploads anything.
set -u
COMMIT=${1:?usage: embryo_copy_rows_sized.sh COMMIT RUN E52_RUN_COMMIT}
RUN=${2:?usage: embryo_copy_rows_sized.sh COMMIT RUN E52_RUN_COMMIT}
E52_RUN_COMMIT=${3:?usage: embryo_copy_rows_sized.sh COMMIT RUN E52_RUN_COMMIT}
case "$RUN" in
  *[!A-Za-z0-9_-]*) echo "refused: RUN takes letters, digits, _ and - only" >&2; exit 90 ;;
esac
export PATH=$HOME/.local/bin:$PATH
S=$HOME/v/scratch/embcrs/$RUN
mkdir -p "$HOME/v/scratch/embcrs"
mkdir "$S" 2> /dev/null || { echo "refused: the run folder $S exists (a run folder is never reused)" >&2; exit 91; }
mkdir "$S/tmp"
export NUMPY_MADVISE_HUGEPAGE=0 TMPDIR="$S/tmp"
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
PY="uv run --frozen --no-sync python -I -B"
BOARD=T2:embryo:test_interp
POSITIONS_SHA=a03611a343faf1d37edf06ce0e2bccc86cce6f14d3580b823864e5faf3ba03f9  # set embryo2, the positions file P
CR_SUMS=experiments/builds/results/embcr_59ad7fa.SHA256SUMS  # C's committed build record (row E50's file)
CR_ENTRY=CRa/T2embryo_s2x_cr.h5ad
READING=experiments/coords/results/e52_copy_rows_sized/reading.json
SIZE_JSON=experiments/coords/results/p3/size_embryo_E7.75.json  # read by code only, never shown
CODE_PATHS="src experiments/builds/build_test.py experiments/density/results/release/count_decision.json experiments/coords/results/p3 boards configs experiments/t1_m1 experiments/t1_hygiene experiments/t1_strong_null uv.lock pyproject.toml"
BASE_ARGS="--method fgw_shape0 --positions floor_exact --expression mixture"
CR_ARGS="--method fgw_shape0 --positions floor_exact --expression mixture --copy-rows --copy-rows-test-ok"
CS_ARGS="--method fgw_shape0 --positions floor_exact --expression mixture --copy-rows --copy-rows-size --copy-rows-test-ok --copy-rows-size-test-ok"
STEP=gate
finish() { echo "$1" > "$S/$STEP.exit"; exit "$1"; }
step_ok() { echo 0 > "$S/$STEP.exit"; }

head=$(git rev-parse HEAD)
[ "$head" = "$COMMIT" ] || { echo "refused: HEAD is not $COMMIT (give the full hash)" >&2; finish 96; }
if [ -n "$(git status --porcelain --untracked-files=no)" ]; then
  echo "refused: tracked files differ from the commit (the build stamps must say clean)" >&2
  finish 98
fi
sha256sum -c docs/S2_RULES.sha256 > /dev/null 2>&1 || finish 97
grep -q '^| E52 ' experiments/harness/TRIES.md || { echo "refused: row E52 is not in TRIES.md" >&2; finish 95; }
case "$E52_RUN_COMMIT" in
  *[!0-9a-f]*) echo "refused: E52_RUN_COMMIT is not a full hash" >&2; finish 88 ;;
esac
[ "${#E52_RUN_COMMIT}" = 40 ] || { echo "refused: E52_RUN_COMMIT is not a full hash" >&2; finish 88; }
git merge-base --is-ancestor "$E52_RUN_COMMIT" HEAD 2> /dev/null \
  || { echo "refused: E52_RUN_COMMIT is not an ancestor of HEAD" >&2; finish 87; }
git show "$E52_RUN_COMMIT:experiments/harness/TRIES.md" 2> /dev/null | grep -q '^| E52 ' \
  || { echo "refused: row E52 is not in TRIES.md at E52_RUN_COMMIT" >&2; finish 86; }
git ls-files --error-unmatch "$READING" > /dev/null 2>&1 \
  || { echo "refused: row E52 has no committed reading (the file must be tracked at this commit)" >&2; finish 94; }
python3 -I -B -c 'import json,sys; r=json.load(open(sys.argv[1])); sys.exit(0 if (r.get("void") is False and r.get("build_test_file") is True and r.get("i", {}).get("holds") is True and r.get("arm") == "fgw_copyrows_sized" and r.get("commit") == sys.argv[2]) else 1)' "$READING" "$E52_RUN_COMMIT" \
  || { echo "refused: row E52's reading is void, is not the run commit's, or its reading (i) does not hold" >&2; finish 92; }
git diff --quiet "$E52_RUN_COMMIT" HEAD -- $CODE_PATHS \
  || { echo "refused: the builder's code or records differ from those row E52's reading ran" >&2; finish 89; }
CR_SHA=$(awk -v e="$CR_ENTRY" '$2 == e {print $1}' "$CR_SUMS" 2> /dev/null)
case "$CR_SHA" in
  *[!0-9a-f]*|"") echo "refused: no committed sha256 of C in $CR_SUMS" >&2; finish 85 ;;
esac
[ "${#CR_SHA}" = 64 ] || { echo "refused: no committed sha256 of C in $CR_SUMS" >&2; finish 85; }
O=runs/embcrs_${head:0:7}_$RUN
mkdir -p runs
mkdir "$O" 2> /dev/null || { echo "refused: the output folder $O exists" >&2; finish 93; }
echo "$head" > "$S/head"
step_ok

STEP=tests
$PY -m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp="$S/tmp/pytest" \
  tests/test_embryo_copy_rows_sized_launcher.py tests/test_embryo_copy_rows_launcher.py \
  tests/test_copy_rows_sized.py tests/test_copy_rows.py tests/test_copy_rows_build.py \
  tests/test_e47_floor_exact.py tests/test_hidden_counts.py tests/test_s2_rules.py > "$S/tests.out" 2>&1 || finish 10
step_ok

export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16
stamp() {  # the stamp's commit and clean fields only
  $PY -c 'import json,sys; s=json.load(open(sys.argv[1])); print(s["commit"], s["clean"])' "$1.stamp.json"
}
build() {  # folder name args: the build, then `vec check`; exit 11 build, 13 check, 14 no stamp, 15 stamp not clean
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

STEP=D
build D T2embryo_s2x $BASE_ARGS || finish $?
got=$(sha256sum "$O/D/T2embryo_s2x.h5ad" | cut -d' ' -f1)
[ "$got" = "$POSITIONS_SHA" ] || { echo "stopped: the positions command does not give the positions file" >&2; finish 20; }
step_ok

STEP=CR
build CR T2embryo_s2x_cr $CR_ARGS || finish $?
got=$(sha256sum "$O/CR/T2embryo_s2x_cr.h5ad" | cut -d' ' -f1)
[ "$got" = "$CR_SHA" ] || { echo "stopped: the copy-rows command does not give C's committed file" >&2; finish 21; }
step_ok

for STEP in CSa CSb; do build "$STEP" T2embryo_s2x_crs $CS_ARGS || finish $?; step_ok; done

STEP=identity
$PY -c '
import json, sys
import anndata as ad
import numpy as np
from vec import t2_build
cr, cs, size_path = (ad.read_h5ad(sys.argv[1]), ad.read_h5ad(sys.argv[2]), sys.argv[3])
def dense(x):
    return np.ascontiguousarray(x.toarray() if hasattr(x, "toarray") else np.asarray(x), dtype=np.float32)
size = json.load(open(size_path))
ok_size = bool(size.get("valid")) and size["preregistered"]["z_mode"] == "raw"
target = float(size["preregistered"]["pred_rms_full"])
a, b = dense(cr.X), dense(cs.X)
C = np.asarray(cr.obsm["spatial_3D"], dtype=np.float32)[:, :3].astype(np.float64)
want = np.ascontiguousarray(t2_build.rescale_to_rms(C, target).astype(np.float32))
got = np.ascontiguousarray(np.asarray(cs.obsm["spatial_3D"])[:, :3], dtype=np.float32)
res = {"size_json_valid_raw": ok_size,
       "x_equal": a.shape == b.shape and a.tobytes() == b.tobytes(),
       "obs_names_equal": list(cr.obs_names) == list(cs.obs_names),
       "var_names_equal": list(cr.var_names) == list(cs.var_names),
       "coords_equal_rescaled_c": want.shape == got.shape and np.array_equal(want.view(np.uint32), got.view(np.uint32))}
for k, v in res.items():
    print(k, bool(v))
sys.exit(0 if all(res.values()) else 1)
' "$O/CR/T2embryo_s2x_cr.h5ad" "$O/CSa/T2embryo_s2x_crs.h5ad" "$SIZE_JSON" > "$S/identity.out" 2> /dev/null || finish 40
step_ok

STEP=done
(cd "$O" && sha256sum D/*.h5ad CR/*.h5ad CSa/*.h5ad CSb/*.h5ad) > "$S/SHA256SUMS" || finish 30
cp "$S/SHA256SUMS" "$O/SHA256SUMS" || finish 30
cmp -s "$O/CSa/T2embryo_s2x_crs.h5ad" "$O/CSb/T2embryo_s2x_crs.h5ad" || finish 31
step_ok
