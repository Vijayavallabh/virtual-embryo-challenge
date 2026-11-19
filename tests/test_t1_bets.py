"""Try T61 (experiments/harness/TRIES.md): `vec.t1_bets`, `experiments/t1_m1/bets_read.py` and the builder's
`--released-bet`, on synthetic stages and hand-made task files only. No released file, statistics file, board file or
score is read."""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pytest
from test_g_mom_verify import _small_inputs
from test_hidden_counts import T1_HELD, t1_env  # noqa: F401  (the synthetic build_t1 fixture)
from test_t1_feng_groups_none import TARGET, _sha, _target_env
from test_t1_scope import GENES, HELD
from test_t1_state_transport import _gmom_ctx

from vec import harness, t1_bets, t1_scope
from vec import t1_state_transport as ts
from vec.boards import Board
from vec.core import counts_of
from vec.t1_scope import build_submission, classify, load_spec, momentum_shift, reweight

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "experiments" / "t1_m1" / "bets_read.py"
BOARD = Board(key="pseudo:t61", task="T1", genes=list(GENES), needs_coords=False, min_cells=100, max_cells=600)
SHARES = {"cm": 0.6, "endocardium": 0.25, "blood": 0.15}


@pytest.fixture
def bt(monkeypatch):
    """The builder, with no permitted-only relabelling recorded, so a relabelling record committed later does not
    change these none builds (the final check of the fix round, item C); tests/test_feng_relabel.py covers the record
    and tests/test_t1_feng_groups_none.py pins the none path's refusal once one routes to recomputed."""
    import build_test

    monkeypatch.setattr(build_test, "_relabel_record", lambda: None)  # main's none path before any relabelling
    return build_test


@pytest.fixture
def br(monkeypatch, tmp_path):
    spec = importlib.util.spec_from_file_location("bets_read_module", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    monkeypatch.setattr(mod, "RESULTS", tmp_path / "results")
    return mod


def _dense(x):
    return np.asarray(x.todense()) if hasattr(x, "todense") else np.asarray(x)


def _classified():
    ref, start, _ = _small_inputs()
    classify(ref, source="E9.5_RNA")
    classify(start, source="E8.5_RNA")
    return ref, start


# ---------------------------------------------------------------------------------------------------------------------
# vec.t1_bets: share_trend


def test_trend_at_zero_is_the_own_mix_and_the_rule_is_own_times_the_ratio():
    own = {"cm": 0.5, "endocardium": 0.3, "blood": 0.2}
    later = {"cm": 600, "endocardium": 200, "blood": 200, "epicardium": 0}
    earlier = {"cm": 300, "endocardium": 300, "blood": 400, "epicardium": 0}
    assert t1_bets.trend_shares(own, later, earlier, 0.0)[0] == own
    got, rec = t1_bets.trend_shares(own, later, earlier, 1.0)
    ratio = {"cm": 0.6 / 0.3, "endocardium": 0.2 / 0.3, "blood": 0.2 / 0.4}
    raw = {c: own[c] * ratio[c] for c in own}
    assert got == pytest.approx({c: v / sum(raw.values()) for c, v in raw.items()})
    assert sum(got.values()) == pytest.approx(1.0) and rec["ratios"] == pytest.approx(ratio)
    half, _ = t1_bets.trend_shares(own, later, earlier, 0.5)
    raw = {c: own[c] * ratio[c] ** 0.5 for c in own}
    assert half == pytest.approx({c: v / sum(raw.values()) for c, v in raw.items()})


def test_trend_has_no_births_and_no_trend_below_min_cells():
    own = {"cm": 0.7, "endocardium": 0.3, "epicardium": 0.0}
    later = {"cm": 500, "endocardium": 19, "epicardium": 400}
    earlier = {"cm": 100, "endocardium": 400, "epicardium": 1}
    got, rec = t1_bets.trend_shares(own, later, earlier, 1.0)
    assert got["epicardium"] == 0.0  # a category the own mix lacks is never born
    assert rec["ratios"]["endocardium"] == 1.0 and rec["trended"] == ["cm"]  # 19 cells: no trend
    with pytest.raises(ValueError):
        t1_bets.trend_shares(own, later, earlier, -0.5)


def test_chamber_counts_counts_own_cells_and_leaves_the_stage_alone():
    ref, _ = _classified()
    before = list(ref.obs.columns)
    counts = t1_bets.chamber_counts(ref, "E9.5_RNA")
    assert list(ref.obs.columns) == before
    cats = ref.obs["category"].astype(str).to_numpy()
    assert counts["cm"] == int((cats == "cm").sum()) > 0 and set(counts) == set(load_spec().in_scope("chamber"))


# ---------------------------------------------------------------------------------------------------------------------
# vec.t1_bets: advance


def test_the_axes_are_momentum_shifts_factor_at_damping_one():
    ref, start = _classified()
    seen = []

    def keep(logf, counts, held):
        seen.append(np.array(logf))
        return logf, None

    pop, _ = reweight(ref, SHARES, 300, 1)
    _, rep = momentum_shift(start, ref, pop, 1.0, HELD, seed=[1, 1], logf_filter=keep,
                            start_fallback=load_spec().fallbacks("E8.5_RNA"),
                            ref_fallback=load_spec().fallbacks("E9.5_RNA"))
    shifted = [r["type"] for r in rep.rows if r["status"] == "shifted"]
    axes, rec = t1_bets.momentum_axes(start, ref, sorted(set(pop.obs["category"])), HELD,
                                      start_fallback=load_spec().fallbacks("E8.5_RNA"),
                                      later_fallback=load_spec().fallbacks("E9.5_RNA"))
    assert shifted and {t for t, v in axes.items() if v is not None} == set(shifted)
    for t, logf in zip(shifted, seen, strict=True):
        assert np.allclose(axes[t], logf, atol=1e-12)
        assert np.all(axes[t][np.isin(GENES, HELD)] == 0)


def test_advance_scores_are_standardised_within_each_category():
    ref, start = _classified()
    cats = sorted(set(ref.obs["category"]))
    axes, _ = t1_bets.momentum_axes(start, ref, cats, HELD)
    z = t1_bets.advance_scores(ref, axes)
    c = ref.obs["category"].astype(str).to_numpy()
    for t, axis in axes.items():
        m = c == t
        if axis is None:
            assert np.all(z[m] == 0)
        else:
            assert z[m].mean() == pytest.approx(0, abs=1e-9) and z[m].std() == pytest.approx(1.0)


def test_advance_at_zero_is_reweight_bit_for_bit_and_keeps_every_quota():
    ref, start = _classified()
    axes, _ = t1_bets.momentum_axes(start, ref, sorted(SHARES), HELD)
    z = t1_bets.advance_scores(ref, axes)
    a, _ = reweight(ref, SHARES, 300, 5)
    b, _ = t1_bets.advance_reweight(ref, SHARES, 300, 5, z, 0.0)
    assert list(a.obs["source_obs"]) == list(b.obs["source_obs"]) and (_dense(a.X) == _dense(b.X)).all()
    c, rep = t1_bets.advance_reweight(ref, SHARES, 300, 5, z, 1.0)
    assert list(c.obs["source_obs"]) != list(a.obs["source_obs"])
    qa = a.obs["category"].value_counts().to_dict()
    assert c.obs["category"].value_counts().to_dict() == qa  # the category shares never move
    # the draw leans towards the axis: higher mean z among drawn cm cells than among the pool's (0 by construction)
    rows = ref.obs_names.get_indexer(c.obs["source_obs"][c.obs["category"] == "cm"])
    assert z[rows].mean() > 0.2
    cm = next(r for r in rep.rows if r["category"] == "cm")["advance"]["bin_shares"]
    assert len(cm) == t1_bets.ADVANCE_BINS and list(cm.values()) == sorted(cm.values())  # quotas rise with z
    reused = next(r for r in rep.rows if r["category"] == "cm")["advance"]["bin_reused"]
    assert set(reused) == set(cm) and all(0.0 <= v <= 1.0 for v in reused.values())  # a fraction per bin, no count


def test_the_bins_are_quotas_at_the_count_never_weights():
    z = np.linspace(-2, 2, 100)
    bins = t1_bets.z_bins(z)
    assert np.array_equal(np.bincount(bins), [20] * 5) and bins[0] == 0 and bins[-1] == 4
    assert np.all(t1_bets.z_bins(np.ones(7)) == 0)
    flat = t1_bets.bin_quotas(z, bins, 0.0, 37)
    assert sum(flat.values()) == 37 and max(flat.values()) - min(flat.values()) <= 1  # tau 0: pool shares
    tilted = t1_bets.bin_quotas(z, bins, 1.0, 37)
    assert sum(tilted.values()) == 37 and tilted[4] > tilted[0]
    w = np.exp(z)
    expect = {b: w[bins == b].sum() / w.sum() * 370 for b in range(5)}
    assert t1_bets.bin_quotas(z, bins, 1.0, 370) == pytest.approx(expect, abs=1.0)  # whole cells at the count


def test_the_median_library_is_the_drawn_cells_own():
    ref, _ = _classified()
    pop, _ = reweight(ref, SHARES, 300, 2)
    assert t1_bets.median_library(pop) == float(np.median(np.asarray(counts_of(pop).sum(axis=1)).ravel())) > 0


def test_assemble_is_build_submissions_last_lines():
    ref, _ = _classified()
    hook_calls = []

    def hook(pop):
        hook_calls.append(pop.n_obs)
        return pop

    sub, _ = build_submission(ref, BOARD, SHARES, 300, 1000003, held_flat=HELD, ref_source="E9.5_RNA",
                              post_cells=hook)
    pop, _ = reweight(ref, SHARES, 300, 1000003, fallback=load_spec().fallbacks("E9.5_RNA"))
    other = t1_bets.assemble(pop, BOARD, hook)
    assert _dense(sub.X).tobytes() == _dense(other.X).tobytes() and hook_calls == [300, 300]


# ---------------------------------------------------------------------------------------------------------------------
# bets_read.py: the file's g_mom step and the arms


def _factors(ref, start):
    return harness.g_mom_pair_factors(ref, "E9.5_RNA", start, "E8.5_RNA")[0]


def test_the_hook_is_the_last_files_own(br, monkeypatch):
    """bets_read.file_hook equals build_test._t1_g_mom's hook under --feng-groups none, on the same factors."""
    ref, start = _classified()
    factors = _factors(ref, start)
    monkeypatch.setattr(br.bt, "load_stage", lambda *a, **k: start)
    monkeypatch.setattr(br.bt.harness, "g_mom_pair_factors", lambda *a, **k: (factors, None))
    monkeypatch.setattr(br.bt.t1_scope, "held_flat_genes", lambda *a, **k: list(HELD))
    groups = br.bt.feng_dirs.label_groups(ref.obs["celltype"].astype(str).unique())
    theirs, _ = br.bt._t1_g_mom(ref, sorted(factors.categories), 1000001, br.bt._g_mom_covered(groups, {}, "none"),
                                feng_direction=False)
    pop, _ = reweight(ref, SHARES, 300, 1000001)
    ours = br.file_hook(factors, 1000001, HELD)
    a, b = ours(pop), theirs(pop)
    assert (_dense(a.layers["counts"]) == _dense(b.layers["counts"])).all()
    edited = (_dense(a.layers["counts"]) != _dense(counts_of(pop))).any(axis=1)
    assert edited.any() and set(pop.obs["category"][edited]) <= {"cm", "endocardium"}


def test_the_cut_arms_reproduce_the_harness_and_read_no_outside_file(br, monkeypatch):
    ref, start, _ = _small_inputs()
    ctx = _gmom_ctx(ref, start, seed=1000002)
    fits = br.gmom_fits(ref, "E9.5_RNA", start, "E8.5_RNA")
    monkeypatch.setattr(br.bt.feng_dirs, "feng_directions", lambda *a, **k: pytest.fail("no Feng direction"))
    monkeypatch.setattr(ts, "r1_unbracketed_labels", lambda *a, **k: pytest.fail("R1's Feng list is not read"))
    arms, depth, reuse = br.cut_arms(ctx, fits)
    assert tuple(arms) == br.CUT_ARMS and set(depth) == set(br.COMPOSITION_ARMS) and set(reuse) == set(br.ADVANCE_ARMS)
    assert all(0.0 <= v <= 1.0 for r in reuse.values() for bins in r.values() for v in bins.values())
    assert all(np.isfinite(v) and v > 0 for v in depth.values())
    assert arms["copy"].X.tobytes() == harness._b0_t1(0.0)(ctx).X.tobytes()
    assert arms["b0"].X.tobytes() == harness._b0_t1(harness.B0_DAMPING)(ctx).X.tobytes()
    distinct = {a: arms[a].X.tobytes() for a in br.CUT_ARMS}
    assert len(set(distinct.values())) == len(br.CUT_ARMS)  # every arm is its own prediction
    assert all(arms[a].X.shape == arms["copy"].X.shape for a in br.CUT_ARMS)


def test_the_library_arms_on_t60s_base_and_on_the_own_mix(br):
    e95, e85 = _classified()
    spec = load_spec()
    fb = {"E8.5": spec.fallbacks("E8.5_RNA"), "E9.5": spec.fallbacks("E9.5_RNA")}
    within = {"surrounding": ["Foregut", "Epithelium"]}
    fits = br.gmom_fits(e95, "E9.5_RNA", e85, "E8.5_RNA")
    axes, _ = t1_bets.momentum_axes(e85, e95, sorted(SHARES), HELD, within=within, start_fallback=fb["E8.5"],
                                    later_fallback=fb["E9.5"])
    z = t1_bets.advance_scores(e95, axes)
    arms, depth, reuse = br.library_arms(e85, e95, SHARES, fits, HELD, 4, within, fb, z, n_cells=300)
    assert set(arms) == set(br.LIB_ARMS) - {"copy_last"} and set(depth) == set(br.COMPOSITION_ARMS)
    assert set(reuse) == set(br.ADVANCE_ARMS)
    assert (_dense(arms["b0_gmom"].X) != _dense(arms["b0"].X)).any()
    base, _ = reweight(e95, SHARES, 300, 4, within=within, fallback=fb["E9.5"])
    assert (_dense(arms["oracle_mix"].X) == _dense(base.X)).all()
    mom, _ = momentum_shift(e85, e95, base, 0.5, HELD, seed=[4, 1], within=within, start_fallback=fb["E8.5"],
                            ref_fallback=fb["E9.5"], min_cells=20, eps=0.5)
    assert (_dense(arms["b0"].X) == _dense(mom.X)).all()
    own = harness.t1_chamber_shares(e95, "E9.5_RNA")
    oc, _ = reweight(e95, own, 300, 4, within=within, fallback=fb["E9.5"])
    assert (_dense(arms["own_copy"].X) == _dense(oc.X)).all()
    assert not (_dense(arms["none_gmom"].X) == _dense(base.X)).all()


def test_the_counts_and_the_row_refusal(br, tmp_path, monkeypatch):
    # b0 is both the gate and a candidate on the cuts, so it counts there: 11 x 2 + 12 x 3 = 58 (T1 300 to 358)
    assert br.TRIES_COUNT == 58 and len(br.CUT_COUNTED) == 11 and len(br.LIB_COUNTED) == 12
    assert br.CUT_GATE in br.CUT_COUNTED and br.CUT_GATE in br.ARM_ORDER
    assert set(br.ARM_ORDER) <= set(br.CUT_COUNTED) and set(br.ARM_ORDER) - {"copy"} <= set(br.LIB_COUNTED)
    tries = tmp_path / "TRIES.md"
    tries.write_text("| T60 (one row) | T1 | x |\n")
    monkeypatch.setattr(br, "TRIES", tries)
    monkeypatch.setattr(br, "tries_rows", lambda path=tries: {"T60"})
    with pytest.raises(SystemExit, match="no row"):
        br.refuse(tmp_path / "out", "x.json")


# ---------------------------------------------------------------------------------------------------------------------
# bets_read.py: the reading


def _arm(score):
    return {"score": score, "groups": {g: score / 100 for g in ("de_recovery", "direction", "distribution",
                                                               "covariance")}}


def _cut(br, gains, gate_shift=0.0, seeds=range(10), per_seed=None, clean=True, commit="c" * 40):
    rows = []
    for s in seeds:
        arms = {a: _arm(50.0) for a in br.CUT_ARMS}
        arms["b0"] = _arm(49.0 + gate_shift)
        arms["copy"] = _arm(49.5)
        for a, g in gains.items():
            arms[a] = _arm(50.0 + (per_seed[a][s] if per_seed and a in per_seed else g))
        rows.append({"seed": s, "floor": 50.0, "arms": arms})
    return {"rows": rows, "stamp": {"commit": commit, "tree_clean": clean}}


def _lib(br, gains, commit="c" * 40):
    rows = []
    for s in range(10):
        arms = {a: _arm(60.0) for a in br.LIB_ARMS}
        arms["oracle_mix"], arms["own_copy"], arms["copy_last"] = _arm(59.0), _arm(59.0), _arm(40.0)
        for a, g in gains.items():
            arms[a] = _arm(60.0 + g)
        rows.append({"seed": s, "arms": arms, "inverted_anchors": []})
    return {"rows": rows, "stamp": {"commit": commit, "tree_clean": True}}


def _committed(br):
    cuts = {c: {s: {"b0": 49.0} for s in range(10)} for c in br.CUTS}
    libs = {lib: {s: {"copy_last": 40.0, "oracle_mix": 59.0} for s in range(10)} for lib in br.LIBRARIES}
    return cuts, libs


def test_an_arm_ahead_on_both_cuts_is_the_output_in_both_readings(br):
    cuts, libs = _committed(br)
    r = br.t61_reading({"chamber": _cut(br, {"gmom_s1": 1.2, "b0_gmom": 0.5}), "whole": _cut(br, {"gmom_s1": 0.8})},
                       {lib: _lib(br, {"gmom_s1": 0.4}) for lib in br.LIBRARIES}, cuts, libs)
    assert r["ahead"] == ["gmom_s1"]  # b0_gmom is behind on the whole cut (0 there): not ahead
    assert (r["reading_S"], r["reading_D"]) == ("gmom_s1", "gmom_s1")


def test_a_library_loss_keeps_the_arm_out_of_reading_d_only(br):
    cuts, libs = _committed(br)
    r = br.t61_reading({"chamber": _cut(br, {"share_trend_l1": 2.0}), "whole": _cut(br, {"share_trend_l1": 1.0})},
                       {"1": _lib(br, {"share_trend_l1": -1.0}), "2": _lib(br, {}), "3": _lib(br, {})}, cuts, libs)
    assert (r["reading_S"], r["reading_D"]) == ("share_trend_l1", "none") and r["composition_arm"]
    assert r["ranges"]["share_trend_l1"]["analogue_range"] == pytest.approx([-1.0, 2.0])
    assert r["ranges"]["share_trend_l1"]["j6_range"] == pytest.approx([0.0, 2.0])
    assert r["ranges"]["share_trend_l1"]["j6_range_d"] == pytest.approx([-1.0, 2.0])  # the library loss is shown


def test_a_depth_linked_output_names_the_best_arm_outside_those_sizes(br):
    cuts, libs = _committed(br)
    gains = {"mom_d1_gmom": 2.0, "b0_gmom": 0.8, "advance_t1": 0.6}
    chamber = _cut(br, gains)
    for row in chamber["rows"]:
        row["median_library_over_base"] = {"advance_t1": 1.08, "share_trend_l1": 1.0}
    r = br.t61_reading({"chamber": chamber, "whole": _cut(br, gains)}, {lib: _lib(br, {}) for lib in br.LIBRARIES},
                       cuts, libs)
    assert r["reading_S"] == "mom_d1_gmom" and r["depth_linked_output"] and r["depth_free_ahead"] == "b0_gmom"
    assert r["median_library_shift"] is None  # mom_d1_gmom draws the copy's cells: no record
    assert r["cuts"]["chamber"]["median_library_over_base"] == pytest.approx({"advance_t1": 1.08, "share_trend_l1": 1})
    lines = br.outcome_lines(r)
    assert any("released data alone or decline" in x and "b0_gmom" in x for x in lines)
    r = br.t61_reading({"chamber": chamber, "whole": _cut(br, {"advance_t1": 0.6})},
                       {lib: _lib(br, {}) for lib in br.LIBRARIES}, cuts, libs)
    assert r["reading_S"] == "advance_t1" and r["median_library_shift"] == {"ratio": pytest.approx(1.08),
                                                                            "beyond": True}
    assert any("beyond 5%" in x for x in br.outcome_lines(r))


R1_TEXT = ("The output arm reverses R1 (the captain's ruling of 2026-11-05): only his ruling naming R1 puts it in a "
           "file")


@pytest.mark.parametrize("arm", ["b0", "b0_gmom", "mom_d0.25_gmom", "mom_d1_gmom"])
def test_a_momentum_output_says_it_reverses_r1(br, arm):
    """B0's momentum on every cell reverses R1 under --feng-groups none: the read prints the line, word for word."""
    assert br.R1_LINE == R1_TEXT and set(br.R1_REVERSING) == {"b0", "b0_gmom", "mom_d0.25_gmom", "mom_d1_gmom"}
    cuts, libs = _committed(br)
    if arm == "b0":  # b0 sits at the gate's committed 49.0: put none_gmom and the rest below it
        gains = {**{a: -3.0 for a in br.CUT_ARMS if a != "b0"}, "none_gmom": -2.0}
        chamber, whole = _cut(br, gains), _cut(br, gains)
    else:
        chamber, whole = _cut(br, {arm: 1.5}), _cut(br, {arm: 1.0})
    r = br.t61_reading({"chamber": chamber, "whole": whole}, {lib: _lib(br, {}) for lib in br.LIBRARIES}, cuts, libs)
    assert r["reading_S"] == arm and r["r1_reversing_output"]
    assert R1_TEXT in br.outcome_lines(r)


@pytest.mark.parametrize("arm", ["copy", "gmom_s1", "share_trend_l0.5", "advance_t1"])
def test_an_output_outside_the_momentum_arms_has_no_r1_line(br, arm):
    cuts, libs = _committed(br)
    r = br.t61_reading({"chamber": _cut(br, {arm: 1.5}), "whole": _cut(br, {arm: 1.0})},
                       {lib: _lib(br, {}) for lib in br.LIBRARIES}, cuts, libs)
    assert r["reading_S"] == arm and not r["r1_reversing_output"] and R1_TEXT not in br.outcome_lines(r)
    r = br.t61_reading({"chamber": _cut(br, {}), "whole": _cut(br, {})}, {lib: _lib(br, {}) for lib in br.LIBRARIES},
                       cuts, libs)
    assert r["reading_S"] == "none" and R1_TEXT not in br.outcome_lines(r)


def test_too_few_seeds_ahead_or_a_small_mean_is_not_ahead(br):
    cuts, libs = _committed(br)
    seven = {"gmom_s1": [1.0] * 7 + [-0.1] * 3}
    r = br.t61_reading({"chamber": _cut(br, {"gmom_s1": 0}, per_seed=seven), "whole": _cut(br, {"gmom_s1": 1.0})},
                       {lib: _lib(br, {}) for lib in br.LIBRARIES}, cuts, libs)
    assert r["ahead"] == [] and r["reading_S"] == r["reading_D"] == "none"
    r = br.t61_reading({"chamber": _cut(br, {"gmom_s1": 0.3}), "whole": _cut(br, {"gmom_s1": 1.0})},
                       {lib: _lib(br, {}) for lib in br.LIBRARIES}, cuts, libs)
    assert r["ahead"] == []  # the bar is strict


def test_the_gate_voids_and_a_void_whole_cut_leaves_no_arm_ahead(br):
    cuts, libs = _committed(br)
    r = br.t61_reading({"chamber": _cut(br, {"gmom_s1": 1.0}, gate_shift=0.03), "whole": _cut(br, {})},
                       {lib: _lib(br, {}) for lib in br.LIBRARIES}, cuts, libs)
    assert r["outcome"] == "VOID"
    r = br.t61_reading({"chamber": _cut(br, {"gmom_s1": 1.0}), "whole": _cut(br, {}, gate_shift=0.03)},
                       {lib: _lib(br, {"gmom_s1": 1.0}) for lib in br.LIBRARIES}, cuts, libs)
    assert r["outcome"] != "VOID" and not r["whole_cut_valid"] and r["reading_S"] == "none"


def test_ties_break_by_the_fixed_order_and_a_late_commit_voids_a_library(br):
    cuts, libs = _committed(br)
    r = br.t61_reading({"chamber": _cut(br, {"advance_t1": 1.0, "copy": 1.0}),
                        "whole": _cut(br, {"advance_t1": 1.0, "copy": 1.0})},
                       {"1": _lib(br, {}, commit="d" * 40), "2": _lib(br, {}), "3": _lib(br, {})}, cuts, libs)
    assert r["reading_S"] == "copy" and r["valid_libraries"] == ["2", "3"]
    assert r["libraries"]["2"]["against_none_gmom"]["copy"]["base"] == "own_none_gmom"
    assert r["libraries"]["2"]["against_none_gmom"]["copy"]["change"] == "own_copy"
    assert r["reading_D"] == "none"  # on the libraries the copy (E9.5's own mix) is 1 point behind own_none_gmom


def test_the_read_prints_no_score(br, tmp_path, monkeypatch, capsys):
    cuts, libs = _committed(br)
    out = tmp_path / "results" / "t61"
    (out / "tasks").mkdir(parents=True)
    for f in br.CUT_FILES.values():
        (out / "tasks" / f).write_text(json.dumps(_cut(br, {"gmom_s1": 1.234567})))
    for lib in br.LIBRARIES:
        (out / "tasks" / f"lib{lib}.json").write_text(json.dumps(_lib(br, {"gmom_s1": 0.5})))
    monkeypatch.setattr(br, "committed_t11", lambda c: cuts[c])
    monkeypatch.setattr(br, "committed_t60", lambda lib: libs[lib])
    monkeypatch.setattr(br.gv, "stamp", lambda started, files: {"commit": "e" * 40, "tree_clean": True, "utc": "x"})
    monkeypatch.setattr(br.gv, "tree_state", lambda: {"commit": "e" * 40, "clean": True})
    res = br.read(out)
    printed = capsys.readouterr().out
    assert res["reading"]["reading_S"] == "gmom_s1"
    assert "1.23" not in printed and "51.2" not in printed and "Outcome" in printed
    text = (out / "t61.txt").read_text()
    assert "+1.23" in text and "over the floor: b0 -1.00" in text and "J6 under D" in text


# ---------------------------------------------------------------------------------------------------------------------
# The builder's --released-bet


def _ruling(bt, tmp_path, monkeypatch, arm="share_trend_l0.5", row="R99"):  # noqa: F811
    tries = tmp_path / "TRIES.md"
    tries.write_text(f"| {row} (the captain's ruling, test) | T1 | rules `{arm}` in | x |\n"
                     f"| R98 (a ruling request) | T1 | rules `{arm}` in |\n"
                     f"| R97 (the captain's ruling, test) | T1 | refuses `{arm}`; rules `{arm}` out |\n"
                     f"| R96 (awaits the captain's ruling) | T1 | rules `{arm}` in |\n"
                     f"| T96 (the captain's ruling, test) | T1 | rules `{arm}` in |\n")
    monkeypatch.setattr(bt, "TRIES_FILE", tries)


def test_the_ruling_guard_needs_a_ruling_row_and_the_dedicated_phrase(bt, tmp_path, monkeypatch):  # noqa: F811
    _ruling(bt, tmp_path, monkeypatch)
    arm = "share_trend_l0.5"
    assert bt._captain_ruling_names("R99", arm)
    assert not bt._captain_ruling_names("R98", arm)  # its first cell is no ruling
    assert not bt._captain_ruling_names("R97", arm)  # a negative mention of the arm
    assert not bt._captain_ruling_names("R96", arm)  # "awaits the captain's ruling" is no ruling
    assert not bt._captain_ruling_names("T96", arm)  # a try row, not a ruling row
    assert not bt._captain_ruling_names("R99", "share_trend_l1")  # another arm


def test_the_option_is_off_by_default_and_refused_without_its_ruling(bt, tmp_path, monkeypatch):  # noqa: F811
    ap = bt.build_parser()
    a = ap.parse_args(["t1", "--out", "x.h5ad"])
    assert a.released_bet is None and a.bet_strength is None and a.captain_ruling is None
    _ruling(bt, tmp_path, monkeypatch)
    monkeypatch.setattr(bt, "load_stage", lambda *a, **k: pytest.fail("nothing is read before the refusal"))
    none = [*TARGET, "--feng-groups", "none"]
    for argv, msg in (
            (["--released-bet", "share_trend", "--bet-strength", "0.5", "--captain-ruling", "R99"], "feng-groups none"),
            ([*none, "--released-bet", "share_trend", "--bet-strength", "0.7", "--captain-ruling", "R99"], "one of"),
            ([*none, "--released-bet", "share_trend", "--bet-strength", "0.5"], "post-exposure block"),
            ([*none, "--released-bet", "share_trend", "--bet-strength", "0.5", "--captain-ruling", "R98"], "ruling"),
            ([*none, "--released-bet", "share_trend", "--bet-strength", "1", "--captain-ruling", "R99"], "ruling"),
            ([*none, "--captain-ruling", "R99"], "with --released-bet only"),
            ([*none, "--released-bet", "share_trend", "--bet-strength", "0.5", "--captain-ruling", "R99", "--base",
              "pooled"], "builds one file")):
        with pytest.raises(SystemExit, match=msg):
            bt.main(["t1", *argv, "--out", str(tmp_path / "T1.h5ad")])


def test_share_trend_is_the_last_file_at_the_trended_mix(bt, t1_env, monkeypatch, tmp_path):  # noqa: F811
    seen = _target_env(bt, t1_env, monkeypatch)
    _ruling(bt, tmp_path, monkeypatch)
    none, bet = tmp_path / "T1_none.h5ad", tmp_path / "T1_trend.h5ad"
    bt.main(["t1", *TARGET, "--feng-groups", "none", "--out", str(none)])
    bt.main(["t1", *TARGET, "--feng-groups", "none", "--released-bet", "share_trend", "--bet-strength", "0.5",
             "--captain-ruling", "R99", "--out", str(bet)])
    rec_none = json.loads(none.with_suffix(".build.json").read_text())
    rec = json.loads(bet.with_suffix(".build.json").read_text())
    assert "released_bet" not in rec_none and rec["released_bet"]["arm"] == "share_trend_l0.5"
    assert rec["released_bet"]["captain_ruling"] == "R99" and rec["shares"] != rec_none["shares"]
    # the same bytes as build_submission called with the trended shares and the same hook
    ref = t1_env.ref.copy()
    lib = np.array([str(x).rsplit("_", 1)[-1] for x in ref.obs_names])
    pool = ref[np.flatnonzero(lib == "3")].copy()
    own = bt.harness.t1_chamber_shares(bt.harness._shallow(pool), "E10.5_RNA")
    trend, _ = t1_bets.trend_shares(own, t1_bets.chamber_counts(ref[np.flatnonzero(lib != "3")], "E10.5_RNA"),
                                    t1_bets.chamber_counts(t1_env.start, "E9.5_RNA"), 0.5)
    hook = seen[-1]["kwargs"]["post_cells"]
    assert seen[-1]["kwargs"]["post_cells"] is not None
    sub, _ = t1_env.build(pool, t1_env.board, t1_scope.load_shares(trend), bt.harness.count_for(t1_env.board,
                                                                                                ref.n_obs),
                          bt.harness.BUILD_SEED_OFFSET + 1, held_flat=list(T1_HELD), ref_source="E10.5_RNA",
                          directions={}, direction_fraction=0.0, direction_groups={}, herd=False, post_cells=hook)
    direct = tmp_path / "direct.h5ad"
    bt.t1_scope.write_checked(sub, direct, t1_env.board)
    assert _sha(bet) == _sha(direct) != _sha(none)
    # (f)(3): the drawn median library over the none file's draw at the same seed, a ratio and no count
    spec = load_spec()
    classify(pool, source="E10.5_RNA")
    seed = bt.harness.BUILD_SEED_OFFSET + 1
    n = bt.harness.count_for(t1_env.board, ref.n_obs)
    fb = spec.fallbacks("E10.5_RNA")
    want = (t1_bets.median_library(reweight(pool, t1_scope.load_shares(trend), n, seed, fallback=fb)[0])
            / t1_bets.median_library(reweight(pool, own, n, seed, fallback=fb)[0]))
    assert rec["released_bet"]["median_library_over_none_draw"] == pytest.approx(want) and want > 0


def test_advance_at_zero_is_the_last_file_and_at_its_strength_moves_it(bt, t1_env, monkeypatch, tmp_path):  # noqa: F811
    _target_env(bt, t1_env, monkeypatch)
    monkeypatch.setitem(bt.RELEASED_BETS, "advance", ("t", (0.0, *t1_bets.TAUS)))
    tries = tmp_path / "TRIES.md"
    tries.write_text("| R99 (the captain's ruling, test) | T1 | rules `advance_t0` in; rules `advance_t1` in |\n")
    monkeypatch.setattr(bt, "TRIES_FILE", tries)
    none, zero, one = tmp_path / "T1_none.h5ad", tmp_path / "T1_adv0.h5ad", tmp_path / "T1_adv1.h5ad"
    bt.main(["t1", *TARGET, "--feng-groups", "none", "--out", str(none)])
    for path, x in ((zero, "0"), (one, "1")):
        bt.main(["t1", *TARGET, "--feng-groups", "none", "--released-bet", "advance", "--bet-strength", x,
                 "--captain-ruling", "R99", "--out", str(path)])
    assert _sha(zero) == _sha(none) != _sha(one)
    rec = json.loads(one.with_suffix(".build.json").read_text())
    assert rec["released_bet"]["arm"] == "advance_t1" and rec["released_bet"]["tau"] == 1.0
    assert rec["n_cells"] == "max_cells" and "realised" not in rec
    assert np.isfinite(rec["released_bet"]["median_library_over_none_draw"])
    rec0 = json.loads(zero.with_suffix(".build.json").read_text())
    assert rec0["released_bet"]["median_library_over_none_draw"] == 1.0  # tau 0 draws the none file's cells
    assert json.loads(none.with_suffix(".build.json").read_text())["shares"] == rec["shares"]
