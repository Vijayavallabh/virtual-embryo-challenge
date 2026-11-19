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
#   bash experiments/builds/build_all.sh embryo2 the embryo's last upload since the amendment of 2026-11-19
#                                                (docs/S2_RULES.md, its section on the embryo; row R9 of TRIES.md): the
#                                                first upload's cells and expression, each cell on a real E7.5 position
#   bash experiments/builds/build_all.sh heart3 the heart's last upload since the amendment of 2026-11-19 for the
#                                                heart (docs/S2_RULES.md, that section; the captain's own named
#                                                exception): heart2's form on the full trim (E10.5 under
#                                                radial_trim_mask at fraction 0), sized to the size file's raw-mode
#                                                core forecast, then the same uniform random half; the half is uploaded
#   bash experiments/builds/build_all.sh embryo3 the embryo's last upload since the later amendment of 2026-11-19
#                                                for the embryo (docs/S2_RULES.md, that section): embryo2's positions
#                                                with the series rendering (the harness's fgw_shape0x_series, try E47b)
#   bash experiments/builds/build_all.sh t1rc   T1's last upload since the later amendment of 2026-11-19 for T1
#                                                (docs/S2_RULES.md, that section; the captain's own say-so): the first
#                                                upload's command plus --feng-groups recomputed (the atlas's direction
#                                                on groups recomputed from its permitted cells, row N54). The builder
#                                                refuses it unless the labeller's record at HEAD routes to recomputed
#                                                and the relabelling's output is in this checkout; built only after row
#                                                T63's read, twice at the default 16 threads. It replaces t1none's file.
#   T3 has no new set: its last upload since the amendment of 2026-11-19 for T3 is the barred branch's s1, which the
#   set t3alt builds.
#
# BUILD_OUT (default runs) and BUILD_THREADS (default 16) exist for the determinism check
# (experiments/builds/determinism_check.sh): the same builds into another folder at another BLAS thread count. Every
# build meant for an upload uses the defaults.
#
# Run from the repo root at a clean, committed tree, with no sync while it runs (the build stamps record the tree).
# The commands are the ones docs/P3_RUNBOOK.md and docs/S2_RULES.md give; this script only runs them in order and
# stops at the first failure. Test-board cell counts are never printed (the builders and `vec check` hide them).
set -euo pipefail
set_name=${1:?usage: build_all.sh s1|s2|heart2|t3alt|t1none|embryo2|heart3|embryo3|t1rc}
T=${BUILD_THREADS:-16}
O=${BUILD_OUT:-runs}
export OMP_NUM_THREADS=$T OPENBLAS_NUM_THREADS=$T MKL_NUM_THREADS=$T NUMPY_MADVISE_HUGEPAGE=0
U=${UV:-uv}
SEED=1000001
EMBRYO_EXPRESSION=mixture  # the embryo s1's expression step, which the set embryo2 shares (s2 names its own, R4)
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
  t1rc)
    sha256sum -c docs/S2_RULES.sha256
    $U run python experiments/builds/build_test.py t1 --g-mom --base lib3 --feng-groups recomputed --out $O/builds/T1_s2rc.h5ad
    files="$O/builds/T1_s2rc.h5ad:T1:test"
    ;;
  embryo2)
    sha256sum -c docs/S2_RULES.sha256
    $U run python experiments/builds/build_test.py embryo --method fgw_shape0 --positions floor_exact \
      --expression "$EMBRYO_EXPRESSION" --out $O/builds/T2embryo_s2x.h5ad
    files="$O/builds/T2embryo_s2x.h5ad:T2:embryo:test_interp"
    ;;
  heart3)
    sha256sum -c docs/S2_RULES.sha256
    $U run python experiments/builds/build_test.py heart --entry core_trim --expression copy --z-mode raw --out $O/builds/T2heart_s2core.h5ad
    $U run python experiments/builds/halve.py $O/builds/T2heart_s2core.h5ad $O/builds/T2heart_s2core_half.h5ad
    files="$O/builds/T2heart_s2core.h5ad:T2:heart:test_extrap $O/builds/T2heart_s2core_half.h5ad:T2:heart:test_extrap"
    ;;
  embryo3)
    sha256sum -c docs/S2_RULES.sha256
    $U run python experiments/builds/build_test.py embryo --method fgw_shape0 --positions floor_exact \
      --expression series --out $O/builds/T2embryo_s2xs.h5ad
    files="$O/builds/T2embryo_s2xs.h5ad:T2:embryo:test_interp"
    ;;
  t3alt)
    sha256sum -c docs/S2_RULES.sha256
    t3 barred s1 $O/t3_build/barred_s1.h5ad
    t3 barred s2 $O/t3_build/barred_s2.h5ad
    files="$O/t3_build/barred_s1.h5ad:T3:bcatenin $O/t3_build/barred_s2.h5ad:T3:bcatenin"
    ;;
  *)
    echo "unknown set $set_name (s1, s2, heart2, t3alt, t1none, embryo2, heart3, embryo3 or t1rc)" >&2
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
