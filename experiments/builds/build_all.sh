#!/bin/bash
# Build and check every test-board file of one submission set on the DGX, after the board-spec intake
# (docs/P3_RUNBOOK.md step 1 and section 5). Nothing is uploaded: each file is written with its build stamp and checked
# locally with `vec check`; server format checks and official uploads are separate calls from the local machine
# (`vec submit FILE --board KEY --from-dgx`), and an official upload needs the captain's go-ahead.
#
#   bash experiments/builds/build_all.sh s1     the s1 entries (T3: the allowed branch's s1, uploaded last)
#   bash experiments/builds/build_all.sh s2     the s2 entries of docs/S2_RULES.md (checked against its hash first)
#
# Since the amendments of 2026-11-13 (docs/S2_RULES.md, the sections of that day; rows R4 to R7 of
# experiments/harness/TRIES.md): T1's s1 is drawn from library 3's cells (--base lib3) and its s2 from the base
# half-way between the pooled base and library 3 (--base mid --cm-share 0.60); the embryo's s2 is the aligned series
# file; T3's refused s1 also carries the wild-type-only factor (--wt-revert). Both heart files are built with --trend
# (b0_trend: the captain's ruling of 2026-11-13, docs/PLAN.md decision log (w); row H11 of TRIES.md).
#   bash experiments/builds/build_all.sh t3alt  the T3 switch files of docs/S2_RULES.md (barred branch s1 and s2)
#   bash experiments/builds/build_all.sh heart2 the heart's second upload since the amendment of 2026-11-18
#                                                (docs/S2_RULES.md, that section): E10.5's cells with their expression
#                                                unchanged, scaled evenly on all three axes, then a uniform random half
#   bash experiments/builds/build_all.sh t1none T1's last upload with no outside data, by the amendment of
#                                                2026-11-19 (docs/S2_RULES.md, that section; row R8 of TRIES.md): the
#                                                first upload's command plus --feng-groups none (Feng's direction taken
#                                                out, nothing in its place; same base, seed and cells; g_mom as built).
#                                                A guard kept from the draft: the set refuses (exit 3) if the rules
#                                                ever hold the draft heading "DRAFT, not in force" again. T1 files move
#                                                with the BLAS thread count, so its two reference builds both use the
#                                                default 16.
#
# BUILD_OUT (default runs) and BUILD_THREADS (default 16) exist for the determinism check
# (experiments/builds/determinism_check.sh): the same builds into another folder at another BLAS thread count. Every
# build meant for an upload uses the defaults.
#
# Run from the repo root at a clean, committed tree, with no sync while it runs (the build stamps record the tree).
# The commands are the ones docs/P3_RUNBOOK.md and docs/S2_RULES.md give; this script only runs them in order and
# stops at the first failure. Test-board cell counts are never printed (the builders and `vec check` hide them).
set -euo pipefail
set_name=${1:?usage: build_all.sh s1|s2|heart2|t3alt|t1none}
T=${BUILD_THREADS:-16}
O=${BUILD_OUT:-runs}
export OMP_NUM_THREADS=$T OPENBLAS_NUM_THREADS=$T MKL_NUM_THREADS=$T NUMPY_MADVISE_HUGEPAGE=0
U=${UV:-uv}
SEED=1000001
EMBRYO_EXPRESSION=mixture  # the embryo s1's expression step (its s2 names its own since the ruling of 2026-11-13)
T3_DECISIONS=experiments/t3_decomp/results/gata4/gata4_decisions.json
T3_RUN_DIR=runs/t3_decomp/gata4_20261028T122425Z
mkdir -p "$O/builds" "$O/t3_build"

if [ -n "$(git status --porcelain --untracked-files=no)" ]; then
  echo "refused: tracked files differ from the commit; commit first (the build stamps must say clean)" >&2
  exit 2
fi
git log --oneline -1

t3() {  # branch submission out [more options]
  $U run python experiments/t3_decomp/build_submission.py build "$T3_DECISIONS" --branch "$1" --submission "$2" \
    --wt data/raw/T2/E8.75.h5ad --run-dir "$T3_RUN_DIR" --board T3:bcatenin --seed "$SEED" --out "$3" "${@:4}"
}

case "$set_name" in
  s1)
    $U run python experiments/builds/build_test.py t1 --g-mom --base lib3 --out $O/builds/T1_s1.h5ad
    $U run python experiments/builds/build_test.py heart --entry s1 --trend --out $O/builds/T2heart_s1.h5ad
    $U run python experiments/builds/build_test.py embryo --method fgw --expression "$EMBRYO_EXPRESSION" \
      --out $O/builds/T2embryo_s1.h5ad
    t3 allowed s1 $O/t3_build/allowed_s1.h5ad
    files="$O/builds/T1_s1.h5ad:T1:test $O/builds/T2heart_s1.h5ad:T2:heart:test_extrap
           $O/builds/T2embryo_s1.h5ad:T2:embryo:test_interp $O/t3_build/allowed_s1.h5ad:T3:bcatenin"
    ;;
  s2)
    sha256sum -c docs/S2_RULES.sha256
    $U run python experiments/builds/build_test.py t1 --g-mom --base mid --cm-share 0.60 --out $O/builds/T1_s2.h5ad
    $U run python experiments/builds/build_test.py heart --entry s2_trim --trend --out $O/builds/T2heart_s2.h5ad
    $U run python experiments/builds/build_test.py embryo --method fgw --move-coupling aligned --expression series \
      --out $O/builds/T2embryo_s2.h5ad
    t3 refused s1 $O/t3_build/refused_s1.h5ad --wt-revert
    files="$O/builds/T1_s2.h5ad:T1:test $O/builds/T2heart_s2.h5ad:T2:heart:test_extrap
           $O/builds/T2embryo_s2.h5ad:T2:embryo:test_interp $O/t3_build/refused_s1.h5ad:T3:bcatenin"
    ;;
  heart2)
    sha256sum -c docs/S2_RULES.sha256
    $U run python experiments/builds/build_test.py heart --entry s1 --expression copy --z-mode raw --out $O/builds/T2heart_s2iso.h5ad
    $U run python experiments/builds/halve.py $O/builds/T2heart_s2iso.h5ad $O/builds/T2heart_s2iso_half.h5ad
    files="$O/builds/T2heart_s2iso.h5ad:T2:heart:test_extrap $O/builds/T2heart_s2iso_half.h5ad:T2:heart:test_extrap"
    ;;
  t1none)
    sha256sum -c docs/S2_RULES.sha256
    if grep -q 'DRAFT, not in force' docs/S2_RULES.md; then echo "refused: the t1none amendment is a draft, not in force" >&2; exit 3; fi
    $U run python experiments/builds/build_test.py t1 --g-mom --base lib3 --feng-groups none --out $O/builds/T1_s2none.h5ad
    files="$O/builds/T1_s2none.h5ad:T1:test"
    ;;
  t3alt)
    sha256sum -c docs/S2_RULES.sha256
    t3 barred s1 $O/t3_build/barred_s1.h5ad
    t3 barred s2 $O/t3_build/barred_s2.h5ad
    files="$O/t3_build/barred_s1.h5ad:T3:bcatenin $O/t3_build/barred_s2.h5ad:T3:bcatenin"
    ;;
  *)
    echo "unknown set $set_name (s1, s2, heart2, t3alt or t1none)" >&2
    exit 2
    ;;
esac

for spec in $files; do
  f=${spec%%:*}
  board=${spec#*:}
  echo "== vec check $f ($board)"
  $U run vec check "$f" --board "$board"
done
echo "built and checked: $set_name"
