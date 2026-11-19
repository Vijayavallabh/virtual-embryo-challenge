#!/usr/bin/env bash
# Rows N55 and N54 (experiments/harness/TRIES.md, as proposed), one step per call, on a DGX campaign checkout
# (~/v/vec-campaign, never the reference checkout). Run from the checkout's root, detached, and wait on the exit file
# with a background poller; never a foreground sleep. Steps, in order:
#   tests                       the relabeller's and the builder's unit tests (synthetic arrays only)
#   identity                    option off: T1_s1 and T1_s2 as the reference list has them, and main's none file
#   n55 RULING                  the dry run on the released E9.5 (before the freeze), its rule applied by code
#   relabel RULING              the permitted-only relabelling (exit 3: it stopped after preflight, route none)
#   listing RULING SHA          N54's listing, after the relabelling's RECORD is committed
#   build                       the check build of `t1 --g-mom --base lib3 --feng-groups recomputed`
#   agreement SHA               once, for the method report only, after the listing is committed
# Each step refuses on an unclean tree (exit 98, as the code's own gates do), records HEAD, runs with 8 BLAS threads
# (16 for the builds, as the reference builds) and NUMPY_MADVISE_HUGEPAGE=0, and writes ~/v/scratch/n54/<STEP>.exit.
# Build output goes to /dev/null: a test build is read by its exit code, its sha256 and its stamp's commit and clean
# fields only (AGENTS.md), never its log, size or cell count.
set -u
STEP=${1:?usage: n54_step.sh tests|identity|n55|relabel|listing|build|agreement [args]}
shift
S=~/v/scratch/n54
mkdir -p "$S/tmp"
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 NUMPY_MADVISE_HUGEPAGE=0 TMPDIR="$S/tmp"
PY="uv run --frozen --no-sync python -I -B"
finish() { echo "$1" > "$S/$STEP.exit"; exit "$1"; }
if [ -n "$(git status --porcelain)" ]; then
  echo "refused: the tree is not clean" >&2
  finish 98
fi
git rev-parse HEAD > "$S/$STEP.head"
stamp() {  # the stamp's commit and clean fields only
  $PY -c 'import json,sys; s=json.load(open(sys.argv[1])); print(s["commit"], s["clean"])' "$1.stamp.json"
}
threads16() { export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16; }

case "$STEP" in
  tests)
    $PY -m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp="$S/tmp/pytest" \
      tests/test_feng_relabel.py tests/test_t1_feng_groups_none.py tests/test_hidden_counts.py \
      tests/test_g_mom_verify.py tests/test_core.py tests/test_s2_rules.py tests/test_t1_late_proxy.py
    finish $? ;;
  identity)
    threads16
    O=runs/n54_identity
    rm -rf "$O" && mkdir -p "$O"
    $PY experiments/builds/build_test.py t1 --g-mom --base lib3 --out "$O/T1_s1.h5ad" > /dev/null 2>&1 || finish 11
    $PY experiments/builds/build_test.py t1 --g-mom --base mid --cm-share 0.60 --out "$O/T1_s2.h5ad" \
      > /dev/null 2>&1 || finish 12
    $PY experiments/builds/build_test.py t1 --g-mom --base lib3 --feng-groups none --out "$O/T1_s2none.h5ad" \
      > /dev/null 2>&1 || finish 13
    (cd "$O" && sha256sum T1_s1.h5ad T1_s2.h5ad T1_s2none.h5ad) > "$S/identity.sha"
    for f in T1_s1 T1_s2 T1_s2none; do echo "$f $(stamp "$O/$f.h5ad")"; done > "$S/identity.stamps"
    finish 0 ;;
  n55)
    $PY experiments/external/feng_relabel_dryrun.py --ruling "${1:?the decision-log entry naming N55}"
    finish $? ;;
  relabel)
    $PY experiments/external/feng_relabel.py --ruling "${1:?the decision-log entry naming N54}"
    rc=$?
    if [ -f data/derived/feng_relabel_c/cd1_stats_relabel.npz ] && [ "$rc" = 0 ]; then
      sha256sum data/derived/feng_relabel_c/cd1_stats_relabel.npz > "$S/relabel.sha"
    fi
    finish $rc ;;
  listing)
    $PY experiments/t1_m1/n54_relabel_listing.py --ruling "${1:?ruling}" --recomputed-sha256 "${2:?sha256}"
    finish $? ;;
  build)
    threads16
    O=runs/builds_rc
    mkdir -p "$O"
    $PY experiments/builds/build_test.py t1 --g-mom --base lib3 --feng-groups recomputed --out "$O/T1_s2rc.h5ad" \
      > /dev/null 2>&1 || finish 21
    (cd "$O" && sha256sum T1_s2rc.h5ad) > "$S/build.sha"
    echo "T1_s2rc $(stamp "$O/T1_s2rc.h5ad")" > "$S/build.stamps"
    finish 0 ;;
  agreement)
    $PY experiments/t1_m1/n54_authors_agreement.py --frozen-sha256 "${1:?sha256}" --method-report
    finish $? ;;
  *)
    echo "unknown step $STEP" >&2
    finish 97 ;;
esac
