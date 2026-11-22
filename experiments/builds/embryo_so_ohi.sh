#!/usr/bin/env bash
# The embryo's two candidates of row E51 (experiments/harness/TRIES.md; specified by a score-blind challenger, job
# job_embryo_so_ohi): SO, T2embryo_s2xs_ov (arm fgw_shape0x_series_ov: O's coordinates with s2xs's expression), and
# O_hi, T2embryo_s2x_ovhi (arm fgw_shape0x_ov_hi: O's coordinates, rows, count, seed and size, every drawn cell topped
# up to the deeper bracket's median depth, no per-gene factor). Test-board builds on a DGX campaign checkout
# (~/v/vec-campaign, never the reference checkout). Run from the root of the checkout or detached worktree it is called
# in, detached, and wait on the exit files with a background poller; never a foreground sleep:
#   bash experiments/builds/embryo_so_ohi.sh <COMMIT> <RUN>
# COMMIT is the full hash of the build commit. RUN names this run's own scratch folder ~/v/scratch/embsohi/<RUN> and its
# output folder runs/embsohi_<commit7>_<RUN>. Neither may exist: a run folder is never reused and nothing is ever
# deleted here, so a second run takes a new name. In order, each step writing ~/v/scratch/embsohi/<RUN>/<STEP>.exit
# (the word running first, then its exit code; a candidate the reading does not pass writes skipped for its steps):
#   gate     HEAD is COMMIT, the tracked tree is clean, docs/S2_RULES.sha256 checks, row E51 is in TRIES.md here and at
#            the reading's run commit (an ancestor of HEAD), row E51's reading is a tracked file of this commit and is
#            not void, and for each candidate its arm's own fields say it passes: arms.<arm>.void false,
#            arms.<arm>.passes true, arms.<arm>.clauses a, b and c all true, and its label and control those of the
#            row; a candidate whose pass the reading does not state so is not built, and with neither the gate refuses.
#            The committed sha256 of P, O and s2xs are found in their build records (below) and must equal the ones
#            written here
#   tests    the embryo builder's and this launcher's targeted tests (synthetic stages only), 8 BLAS threads
#   V        the specification's validation-board dry runs (T2:embryo:val_interp) of each candidate's command, into the
#            scratch folder, read by exit code only
#   D        the positions file P by the default command; it must give P_SHA or everything stops
#   O        the overlay file O by its committed command; it must give O_SHA or everything stops
#   S        the series file s2xs by its committed command; it must give S_SHA or everything stops
#   SOa SOb  T2embryo_s2xs_ov.h5ad by the specification's command, twice, into separate folders
#   OHa OHb  T2embryo_s2x_ovhi.h5ad by the specification's command, twice, into separate folders
#   identity booleans only, to identity.out: SO's obsm coordinates equal O's bit for bit and its X, obs_names and
#            var_names equal s2xs's; O_hi's coordinates, obs_names and var_names equal O's (O and s2xs as built in
#            steps O and S of this same run)
#   done     one sha256 list over every built file; SOa and SOb, OHa and OHb compared byte for byte
# Where the three committed files and their commands come from:
#   P     a03611a343fa...: experiments/builds/results/embryo2_2dccfee.SHA256SUMS, builds/T2embryo_s2x.h5ad (set embryo2
#         of experiments/builds/build_all.sh; the same bytes as D/T2embryo_s2x.h5ad of embov_3b813b5.SHA256SUMS); its
#         command is BASE_ARGS, as embryo_overlay.sh's step D
#   O     1cfc91be3a20...: experiments/builds/results/embov_3b813b5.SHA256SUMS, OVa/T2embryo_s2x_ov.h5ad (built by
#         experiments/builds/embryo_overlay.sh at 3b813b5); its command is that launcher's OV_ARGS
#   s2xs  bc50d012cd2c...: experiments/builds/results/embryo3_3d52b83.SHA256SUMS, builds/T2embryo_s2xs.h5ad (set embryo3
#         of experiments/builds/build_all.sh at 3d52b83); its command is that set's, S_ARGS
# A build step also fails (15) unless its stamp says this commit and clean True. Builds run at 16 BLAS threads, as the
# reference builds. Every builder and `vec check` output goes to /dev/null: a test build is read by its exit code, its
# sha256 and its stamp's commit and clean fields only (AGENTS.md), never its log, size, cell count or sidecar record. A
# build that exits non-zero means no file, and the reason is not read here. The flag --shape-interp-test-ok is in the
# specification's commands. Nothing here uploads anything.
set -u
COMMIT=${1:?usage: embryo_so_ohi.sh COMMIT RUN}
RUN=${2:?usage: embryo_so_ohi.sh COMMIT RUN}
case "$RUN" in
  *[!A-Za-z0-9_-]*) echo "refused: RUN takes letters, digits, _ and - only" >&2; exit 90 ;;
esac
case "$COMMIT" in
  *[!0-9a-f]*) echo "refused: COMMIT must be a full hash" >&2; exit 90 ;;
esac
[ "${#COMMIT}" = 40 ] || { echo "refused: COMMIT must be a full hash" >&2; exit 90; }
export PATH=$HOME/.local/bin:$PATH
S=$HOME/v/scratch/embsohi/$RUN
mkdir -p "$HOME/v/scratch/embsohi"
mkdir "$S" 2> /dev/null || { echo "refused: the run folder $S exists (a run folder is never reused)" >&2; exit 91; }
mkdir "$S/tmp"
export NUMPY_MADVISE_HUGEPAGE=0 TMPDIR="$S/tmp"
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
export UV_PROJECT_ENVIRONMENT="$PWD/.venv"
PY="uv run --frozen --no-sync python -I -B"
BOARD=T2:embryo:test_interp
VAL_BOARD=T2:embryo:val_interp
P_SHA=a03611a343faf1d37edf06ce0e2bccc86cce6f14d3580b823864e5faf3ba03f9  # set embryo2, the positions file P
O_SHA=1cfc91be3a20bfdd1a4f7bcfebd7316e1143c6d4048e578c1f30bcfbe1e0e1ea  # embov_3b813b5, the overlay file O
S_SHA=bc50d012cd2c9481cd18c1f5bcf5b01093e22fd16ac65582a5c9f34321f13115  # set embryo3, the series file s2xs
P_SUMS=experiments/builds/results/embryo2_2dccfee.SHA256SUMS
P_ENTRY=builds/T2embryo_s2x.h5ad
O_SUMS=experiments/builds/results/embov_3b813b5.SHA256SUMS
O_ENTRY=OVa/T2embryo_s2x_ov.h5ad
S_SUMS=experiments/builds/results/embryo3_3d52b83.SHA256SUMS
S_ENTRY=builds/T2embryo_s2xs.h5ad
READING=experiments/coords/results/e51_so_ohi/reading.json
BASE_ARGS="--method fgw_shape0 --positions floor_exact --expression mixture"
OV_ARGS="--method fgw_shape0 --positions floor_exact --expression mixture --shape-interp overlay --shape-interp-test-ok"
S_ARGS="--method fgw_shape0 --positions floor_exact --expression series"
SO_ARGS="--method fgw_shape0 --positions floor_exact --expression series --shape-interp overlay --shape-interp-test-ok"
OH_ARGS="--method fgw_shape0 --positions floor_exact --expression render_hi --shape-interp overlay --shape-interp-test-ok"
SO_DRY_ARGS="--method fgw_shape0 --positions floor_exact --expression series --shape-interp overlay"
OH_DRY_ARGS="--method fgw_shape0 --positions floor_exact --expression render_hi --shape-interp overlay"
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
grep -q '^| E51 ' experiments/harness/TRIES.md || { echo "refused: row E51 is not in TRIES.md" >&2; finish 95; }
git ls-files --error-unmatch "$READING" > /dev/null 2>&1 \
  || { echo "refused: row E51 has no committed reading (the file must be tracked at this commit)" >&2; finish 94; }
PASSES=$(python3 -I -B -c 'import json,re,sys
r = json.load(open(sys.argv[1]))
arms = r.get("arms") if isinstance(r.get("arms"), dict) else {}
row = {"fgw_shape0x_series_ov": ("SO", "fgw_shape0x_series"), "fgw_shape0x_ov_hi": ("O_hi", "fgw_shape0x_ov")}
def passes(arm):
    a = arms.get(arm)
    if not isinstance(a, dict):
        return False
    cl = a.get("clauses")
    return ((a.get("label"), a.get("control")) == row[arm] and a.get("void") is False and a.get("passes") is True
            and isinstance(cl, dict) and set(cl) == {"a", "b", "c"} and all(v is True for v in cl.values()))
c = str(r.get("commit"))
so, oh = passes("fgw_shape0x_series_ov"), passes("fgw_shape0x_ov_hi")
if r.get("void") is not False or not re.fullmatch("[0-9a-f]{40}", c) or not (so or oh):
    sys.exit(1)
print(c, int(so), int(oh))' "$READING" 2> /dev/null) \
  || { echo "refused: row E51's reading is void, names no run commit, or states no candidate's pass" >&2; finish 92; }
read -r RUNC SO_OK OH_OK <<< "$PASSES"
git merge-base --is-ancestor "$RUNC" HEAD 2> /dev/null \
  || { echo "refused: the reading's run commit is not an ancestor of HEAD" >&2; finish 88; }
git show "$RUNC:experiments/harness/TRIES.md" 2> /dev/null | grep -q '^| E51 ' \
  || { echo "refused: row E51 is not in TRIES.md at the reading's run commit" >&2; finish 86; }
for rec in "P_SHA $P_SUMS $P_ENTRY" "O_SHA $O_SUMS $O_ENTRY" "S_SHA $S_SUMS $S_ENTRY"; do
  read -r var sums entry <<< "$rec"
  got=$(awk -v e="$entry" '$2 == e {print $1}' "$sums" 2> /dev/null)
  [ "$got" = "${!var}" ] || { echo "refused: $sums does not give $var for $entry" >&2; finish 85; }
done
O=runs/embsohi_${head:0:7}_$RUN
mkdir -p runs
mkdir "$O" 2> /dev/null || { echo "refused: the output folder $O exists" >&2; finish 93; }
echo "$head" > "$S/head"
printf 'SO %s\nO_hi %s\n' "$SO_OK" "$OH_OK" > "$S/candidates"
step_ok

STEP=tests
echo running > "$S/$STEP.exit"
$PY -m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp="$S/tmp/pytest" \
  tests/test_embryo_so_ohi_launcher.py tests/test_embryo_overlay.py tests/test_shape_interp.py \
  tests/test_shape_interp_build.py tests/test_e47_floor_exact.py tests/test_e51_so_ohi_read.py \
  tests/test_hidden_counts.py tests/test_s2_rules.py > "$S/tests.out" 2>&1 || finish 10
step_ok

export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16
STEP=V
echo running > "$S/$STEP.exit"
mkdir "$S/V"
if [ "$SO_OK" = 1 ]; then
  $PY experiments/builds/build_test.py embryo --board "$VAL_BOARD" $SO_DRY_ARGS \
    --out "$S/V/T2embryo_s2xs_ov_dry.h5ad" > /dev/null 2>&1 || finish 16
fi
if [ "$OH_OK" = 1 ]; then
  $PY experiments/builds/build_test.py embryo --board "$VAL_BOARD" $OH_DRY_ARGS \
    --out "$S/V/T2embryo_s2x_ovhi_dry.h5ad" > /dev/null 2>&1 || finish 17
fi
step_ok

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

STEP=D
echo running > "$S/$STEP.exit"
build D T2embryo_s2x $BASE_ARGS || finish $?
got=$(sha256sum "$O/D/T2embryo_s2x.h5ad" | cut -d' ' -f1)
[ "$got" = "$P_SHA" ] || { echo "stopped: the default command does not give the positions file P" >&2; finish 20; }
step_ok

STEP=O
echo running > "$S/$STEP.exit"
build O T2embryo_s2x_ov $OV_ARGS || finish $?
got=$(sha256sum "$O/O/T2embryo_s2x_ov.h5ad" | cut -d' ' -f1)
[ "$got" = "$O_SHA" ] || { echo "stopped: the overlay command does not give O's committed file" >&2; finish 21; }
step_ok

STEP=S
echo running > "$S/$STEP.exit"
build S T2embryo_s2xs $S_ARGS || finish $?
got=$(sha256sum "$O/S/T2embryo_s2xs.h5ad" | cut -d' ' -f1)
[ "$got" = "$S_SHA" ] || { echo "stopped: the series command does not give s2xs's committed file" >&2; finish 22; }
step_ok

for STEP in SOa SOb; do
  if [ "$SO_OK" = 1 ]; then
    echo running > "$S/$STEP.exit"
    build "$STEP" T2embryo_s2xs_ov $SO_ARGS || finish $?
    step_ok
  else
    echo skipped > "$S/$STEP.exit"
  fi
done
for STEP in OHa OHb; do
  if [ "$OH_OK" = 1 ]; then
    echo running > "$S/$STEP.exit"
    build "$STEP" T2embryo_s2x_ovhi $OH_ARGS || finish $?
    step_ok
  else
    echo skipped > "$S/$STEP.exit"
  fi
done

STEP=identity
echo running > "$S/$STEP.exit"
SO_FILE=-
OH_FILE=-
[ "$SO_OK" = 1 ] && SO_FILE="$O/SOa/T2embryo_s2xs_ov.h5ad"
[ "$OH_OK" = 1 ] && OH_FILE="$O/OHa/T2embryo_s2x_ovhi.h5ad"
$PY -c '
import sys
import anndata as ad
import numpy as np
o_path, s_path, so_path, oh_path = sys.argv[1:5]
def x(a):
    return np.ascontiguousarray(a.X.toarray() if hasattr(a.X, "toarray") else np.asarray(a.X))
def same(p, q):
    p, q = np.ascontiguousarray(np.asarray(p)), np.ascontiguousarray(np.asarray(q))
    return p.dtype == q.dtype and p.shape == q.shape and p.tobytes() == q.tobytes()
def coords(a):
    return a.obsm["spatial_3D"]
o = ad.read_h5ad(o_path)
res = {}
if so_path != "-":
    s, so = ad.read_h5ad(s_path), ad.read_h5ad(so_path)
    res["so_coords_equal_o"] = same(coords(so), coords(o))
    res["so_x_equal_s2xs"] = same(x(so), x(s))
    res["so_obs_names_equal_s2xs"] = list(so.obs_names) == list(s.obs_names)
    res["so_var_names_equal_s2xs"] = list(so.var_names) == list(s.var_names)
if oh_path != "-":
    oh = ad.read_h5ad(oh_path)
    res["ohi_coords_equal_o"] = same(coords(oh), coords(o))
    res["ohi_obs_names_equal_o"] = list(oh.obs_names) == list(o.obs_names)
    res["ohi_var_names_equal_o"] = list(oh.var_names) == list(o.var_names)
for k, v in res.items():
    print(k, bool(v))
sys.exit(0 if res and all(res.values()) else 1)
' "$O/O/T2embryo_s2x_ov.h5ad" "$O/S/T2embryo_s2xs.h5ad" "$SO_FILE" "$OH_FILE" > "$S/identity.out" 2> /dev/null \
  || finish 40
step_ok

STEP=done
echo running > "$S/$STEP.exit"
(cd "$O" && sha256sum D/*.h5ad O/*.h5ad S/*.h5ad) > "$S/SHA256SUMS" || finish 30
for f in SOa SOb OHa OHb; do
  if [ -d "$O/$f" ]; then (cd "$O" && sha256sum "$f"/*.h5ad) >> "$S/SHA256SUMS" || finish 30; fi
done
(cd "$S/V" && sha256sum ./*.h5ad) > "$S/V.SHA256SUMS" 2> /dev/null || finish 30
cp "$S/SHA256SUMS" "$O/SHA256SUMS" || finish 30
if [ "$SO_OK" = 1 ]; then
  cmp -s "$O/SOa/T2embryo_s2xs_ov.h5ad" "$O/SOb/T2embryo_s2xs_ov.h5ad" || finish 31
fi
if [ "$OH_OK" = 1 ]; then
  cmp -s "$O/OHa/T2embryo_s2x_ovhi.h5ad" "$O/OHb/T2embryo_s2x_ovhi.h5ad" || finish 32
fi
step_ok
