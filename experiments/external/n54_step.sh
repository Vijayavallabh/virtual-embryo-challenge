#!/usr/bin/env bash
# Rows N55 and N54 (experiments/harness/TRIES.md, as proposed), one step per call, on a DGX campaign checkout
# (~/v/vec-campaign, never the reference checkout). Run from the checkout's root, detached, and wait on the exit file
# with a background poller; never a foreground sleep. Steps, in order:
#   tests                       the relabeller's and the builder's unit tests (synthetic arrays only): every test file
#                               that imports build_test, feng_dirs, feng_relabel or the harness's late proxy, with
#                               test_core, test_s2_rules, test_git_helpers and experiments/external/test_filters.py
#   identity                    option off: T1_s1 and T1_s2 as the reference list has them, and main's none file
#   n55 RULING                  the dry run on the released E9.5 (before the freeze), its rule applied by code (exit 3:
#                               it stopped, outcome stopped; any other non-zero exit after it started: --record-kill)
#   relabel RULING              the permitted-only relabelling (exit 3: it stopped after preflight, route none; any
#                               other non-zero exit: --record-kill writes the stop record if the run had started)
#   listing RULING SHA          N54's listing, after the relabelling's RECORD is committed
#   build                       the check build of `t1 --g-mom --base lib3 --feng-groups recomputed`
#   build2                      the same command again into its own file (T1_s2rc_b.h5ad, build2.sha, build2.stamps)
#   agreement SHA UPLOAD_RULING once, for the method report only, after the captain's ruling on T1's last upload
# Each step refuses on an unclean tree (exit 98, as the code's own gates do), records HEAD, runs with 8 BLAS threads
# (16 for the builds, as the reference builds) and NUMPY_MADVISE_HUGEPAGE=0, and writes ~/v/scratch/n54/<STEP>.exit.
# Build output goes to /dev/null: a test build is read by its exit code, its sha256 and its stamp's commit and clean
# fields only (AGENTS.md), never its log, size or cell count.
set -u
STEP=${1:?usage: n54_step.sh tests|identity|n55|relabel|listing|build|build2|agreement [args]}
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
      tests/test_g_mom_verify.py tests/test_core.py tests/test_s2_rules.py tests/test_t1_late_proxy.py \
      tests/test_b0_trend.py tests/test_build_heart_floor_iso.py tests/test_build_heart_iso_trim.py \
      tests/test_build_heart_core_trim.py tests/test_build_heart_organ_scope.py tests/test_e47_floor_exact.py \
      tests/test_h25_heart_forms.py tests/test_h6_cfm.py tests/test_heart_trend_build.py \
      tests/test_n12_hygiene_size.py tests/test_n2_momentum_size.py tests/test_n48_m1_direction_shrink.py \
      tests/test_t1_cal_read.py tests/test_t1_ctl_range_read.py tests/test_t1_eb_read.py tests/test_t1_eps_read.py \
      tests/test_t1_gmom_ctl_read.py tests/test_t1_hygiene.py tests/test_t1_part_depth_read.py \
      tests/test_t1_part_range_read.py tests/test_t1_reach_read.py tests/test_t1_scope.py tests/test_t1_share_read.py \
      tests/test_t1_split_read.py tests/test_t1_state_transport.py tests/test_t1_strong_null.py \
      tests/test_t1_zone_read.py tests/test_t2_latent.py tests/test_h26_heart_scope_forms.py tests/test_t1_bets.py \
      tests/test_t62_recomputed_read.py tests/test_h27_heart_forecast_size.py \
      tests/test_t63_double_step_read.py tests/test_h28_heart_trend_strength.py tests/test_h29a_heart_trend_forms.py \
      tests/test_heart_candidates.py tests/test_h29b_heart_organ.py tests/test_heart_organ_launchers.py \
      tests/test_t1_candidates_build.py tests/test_occupancy_frame.py tests/test_shape_interp_build.py \
      tests/test_t64_heldmom_read.py tests/test_t65_ms_read.py tests/test_copy_rows_build.py tests/test_copy_rows.py \
      tests/test_e50_copy_rows_read.py tests/test_embryo_overlay.py tests/test_embryo_copy_rows_launcher.py \
      tests/test_embryo_copy_rows_sized_launcher.py tests/test_e52_copy_rows_sized_read.py \
      tests/test_embryo_so_ohi_launcher.py \
      tests/test_series_back_build.py tests/test_ovsb_step0.py tests/test_embryo_ovsb_launcher.py \
      tests/test_t1_fraction2_build.py tests/test_t1_fraction2_launcher.py \
      tests/test_t1_fraction_line_build.py tests/test_t1_fraction_line_launcher.py tests/test_t1_timing.py \
      tests/test_t1_rc18.py tests/test_jaxblood_read.py tests/test_jax_blood.py tests/test_k25_dispersion.py \
      tests/test_t1_damping_line_build.py tests/test_t1_damping_line_launcher.py tests/test_cmjax_read.py \
      tests/test_jax_cm.py tests/test_t1_line3_build.py tests/test_t1_c_line_launcher.py \
      tests/test_t1_ms_build.py tests/test_t1_ms_launcher.py \
      tests/test_t1_msb_build.py tests/test_t1_ms_combined_launcher.py \
      tests/test_t1_line5_build.py tests/test_t1_b2d_line_launcher.py \
      tests/test_t1_l67pb_build.py tests/test_t1_l67pb_launcher.py \
      tests/test_t1_hm_build.py tests/test_t1_hm_launcher.py \
      tests/test_heart_trend_line.py tests/test_heart_trend_line_launcher.py \
      tests/test_heart_trend_line_d.py tests/test_heart_trend_line_d_launcher.py \
      tests/test_heart_trend_line_f.py tests/test_heart_trend_line_f_launcher.py \
      tests/test_e55_overlay_stretch.py tests/test_h34_heart_organ_trend.py \
      tests/test_git_helpers.py experiments/external/test_filters.py
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
    rc=$?
    if [ "$rc" != 0 ] && [ "$rc" != 3 ] && [ -f data/derived/feng_relabel_n55/STARTED.json ] \
        && [ ! -f experiments/external/results/feng_relabel_dryrun/n55.json ]; then
      # killed after it started (SIGKILL, OOM): record the stop (n55_rule (iv)); exit 3 when recorded
      $PY experiments/external/feng_relabel_dryrun.py --record-kill
      rc=$?
    fi
    finish $rc ;;
  relabel)
    $PY experiments/external/feng_relabel.py --ruling "${1:?the decision-log entry naming N54}"
    rc=$?
    if [ "$rc" != 0 ] && [ "$rc" != 3 ] && [ -f data/derived/feng_relabel_c/STARTED.json ] \
        && [ ! -f data/derived/feng_relabel_c/MANIFEST.json ]; then
      # killed after it started (SIGKILL, OOM): record the stop (terminal.on_stop, route none); exit 3 when recorded
      $PY experiments/external/feng_relabel.py --record-kill
      rc=$?
    fi
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
  build2)  # the same command into its own file, so the first build's file and sha stay as they were
    threads16
    O=runs/builds_rc
    mkdir -p "$O"
    $PY experiments/builds/build_test.py t1 --g-mom --base lib3 --feng-groups recomputed \
      --out "$O/T1_s2rc_b.h5ad" > /dev/null 2>&1 || finish 22
    (cd "$O" && sha256sum T1_s2rc_b.h5ad) > "$S/build2.sha"
    echo "T1_s2rc_b $(stamp "$O/T1_s2rc_b.h5ad")" > "$S/build2.stamps"
    finish 0 ;;
  agreement)
    $PY experiments/t1_m1/n54_authors_agreement.py --frozen-sha256 "${1:?sha256}" --method-report \
      --upload-ruling "${2:?the captain's decision-log entry on T1's last upload}"
    finish $? ;;
  *)
    echo "unknown step $STEP" >&2
    finish 97 ;;
esac
