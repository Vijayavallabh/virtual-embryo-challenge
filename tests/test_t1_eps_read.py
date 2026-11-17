"""experiments/t1_m1/eps_read.py (try T52 of experiments/harness/TRIES.md) on synthetic data: the directions made with
another pseudocount, the arms, the placebos and the reading. No released file, real statistics file, board file or
score is read."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import re
import sys
import types

import numpy as np
import pytest
from test_t1_late_proxy import _manual, _stats
from test_t1_lib_step_read import GOLDEN_ARMS_5BC7, _arm
from test_t1_step_read import _arms_world

from vec import harness, t1_scope
from vec.boards import REPO_ROOT

SCRIPT = REPO_ROOT / "experiments" / "t1_m1" / "eps_read.py"
NEW = ("m1_r_gmom_e01", "m1_r_gmom_e01_x15", "m1_r_gmom_e01_x2", "m1_r_gmom_e002", "m1_r_gmom_e002_x15",
       "m1_r_gmom_e002_x2")
NEW_PERM = ("m1_r_gmom_e01_x15_perm", "m1_r_gmom_e002_x15_perm")
NEW_SPERM = ("m1_r_gmom_e01_x15_sperm", "m1_r_gmom_e002_x15_sperm")
SPERM = ("m1_r_gmom_sperm", *NEW_SPERM)
UNIFORM = ("m1_r_gmom_x15", "m1_r_gmom_x2", "m1_r_gmom_x3", "m1_r_gmom_x4")
CHAMBER = ("pooled", "lib1", "lib2", "lib3")
TRUTHS = (*CHAMBER, "whole")
B0, BASE = 50.0, 53.0
STAMP = {"commit": "f" * 40, "tree_clean": True, "utc": "x", "inputs": {}}


@pytest.fixture
def ep(monkeypatch, tmp_path):
    spec = importlib.util.spec_from_file_location("eps_read_module", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    monkeypatch.setattr(mod, "RESULTS", tmp_path / "results")
    monkeypatch.setattr(mod.gv, "stamp", lambda started, files: dict(STAMP))
    return mod


def _prints_a_count(txt: str, counts) -> bool:
    body = txt.rsplit("\nCommit ", 1)[0]
    return any(re.search(rf"(?<![\d.]){n}(?![\d.])", body) for n in counts)


# ---------------------------------------------------------------------------------------------------------------------
# The arms' names, fractions and pseudocounts


def test_the_arms_of_each_truth(ep):
    for k in CHAMBER:
        assert ep.arms_of(k) == ("b0", "m1_r_gmom", *UNIFORM, "m1_r_gmom_perm", *NEW, *NEW_PERM, *SPERM)
    assert ep.arms_of("whole") == ("b0", "m1_r_gmom", "m1_r_gmom_x15", "m1_r_gmom_x2", "m1_r_gmom_x3", *NEW)
    spec = ep.arm_spec("lib2")
    assert spec["m1_r_gmom"] == {"fraction": 0.2, "pseudocount": 0.5}
    assert [spec[a]["fraction"] for a in UNIFORM] == pytest.approx([0.3, 0.4, 0.6, 0.8])
    assert [spec[a]["fraction"] for a in NEW] == pytest.approx([0.2, 0.3, 0.4] * 2)
    assert [spec[a]["pseudocount"] for a in NEW] == [0.1] * 3 + [0.02] * 3
    assert spec["m1_r_gmom_e01_x15_perm"] == {"fraction": pytest.approx(0.3), "pseudocount": 0.1}
    assert spec["m1_r_gmom_e002_x15_perm"] == {"fraction": pytest.approx(0.3), "pseudocount": 0.02}
    assert spec["m1_r_gmom_perm"] == {"fraction": 0.2, "pseudocount": 0.5}
    assert spec["m1_r_gmom_sperm"] == {"fraction": 0.2, "pseudocount": 0.5}
    assert spec["m1_r_gmom_e01_x15_sperm"] == {"fraction": pytest.approx(0.3), "pseudocount": 0.1}
    assert spec["m1_r_gmom_e002_x15_sperm"] == {"fraction": pytest.approx(0.3), "pseudocount": 0.02}
    whole = ep.arm_spec("whole")
    assert not {"m1_r_gmom_x4", "m1_r_gmom_perm", *NEW_PERM, *SPERM} & set(whole)
    assert ep.PLACEBO_OF == {"m1_r_gmom_perm": "m1_r_gmom", "m1_r_gmom_sperm": "m1_r_gmom",
                             **dict(zip(NEW_PERM, ("m1_r_gmom_e01_x15", "m1_r_gmom_e002_x15"), strict=True)),
                             **dict(zip(NEW_SPERM, ("m1_r_gmom_e01_x15", "m1_r_gmom_e002_x15"), strict=True))}
    assert ep.SPERM_STREAM not in (ep.PLACEBO_STREAM, 0, 1, 2)


# ---------------------------------------------------------------------------------------------------------------------
# The directions on synthetic statistics


@pytest.fixture
def stats(tmp_path, monkeypatch, ep):
    monkeypatch.setattr(ep.feng_dirs, "REPO", tmp_path)
    ctx = _arms_world()[0]
    genes = list(ctx.inputs[ctx.task.reference_role].var_names)
    return _stats(tmp_path / "data" / "external" / "feng2022_gse193346_geo" / "cd1_stats.npz", genes), genes


def test_the_new_directions_are_lfc_ratio_at_their_pseudocount_with_the_entrys_groups_and_nans(ep, stats):
    path, genes = stats
    genes = [*genes, "not_in_feng"]  # an unmapped gene: NaN in every direction
    entry, rec = ep.feng_dirs.feng_directions(genes, stats_path=path)
    assert rec["pseudocount_cp10k"] == 0.5
    for eps in (0.1, 0.02):
        new, rec_new = ep.eps_directions(genes, eps, stats_path=path)
        assert rec_new["pseudocount_cp10k"] == eps and set(new) == set(entry) and rec_new["groups"] == rec["groups"]
        for g in entry:
            labels = tuple(rec["groups"][g]["feng_labels"])
            want = _manual(path, "E9.5", "E14.5", labels, eps=eps)
            np.testing.assert_allclose(new[g][:-1], want, rtol=1e-12)
            st = ep.dl.load_stats(path, fields=ep.dl.STAT_FIELDS)
            a, b = ep.dl.feng_sel(st, "E9.5", labels), ep.dl.feng_sel(st, "E14.5", labels)
            np.testing.assert_array_equal(new[g][:-1], ep.dl.lfc_ratio(a, b, eps))
            assert np.array_equal(np.isnan(new[g]), np.isnan(entry[g])) and np.isnan(new[g][-1])
            assert not np.allclose(new[g][:-1], entry[g][:-1])
        assert ep.directions_sha256(new).keys() == ep.directions_sha256(entry).keys()
        assert all(ep.directions_sha256(new)[g] != ep.directions_sha256(entry)[g] for g in entry)


def test_the_idle_pseudocount_gives_the_entrys_directions_bit_for_bit(ep, stats):
    path, genes = stats
    entry, _ = ep.feng_dirs.feng_directions(genes, stats_path=path)
    same, _ = ep.eps_directions(genes, 0.5, stats_path=path)
    assert ep.directions_sha256(same) == ep.directions_sha256(entry)


def test_no_change_reaches_the_config_file_or_the_entrys_build(ep, stats):
    """The config file and a loaded config object are unchanged after the new directions are made, and the entry's arm
    built after a new arm in the same process is the entry's arm built before it, bit for bit."""
    path, genes = stats
    before_file = hashlib.sha256(ep.feng_dirs.CONFIG.read_bytes()).hexdigest()
    cfg = ep.feng_dirs.load_config()
    snapshot = json.dumps(cfg, sort_keys=True)
    ctx, groups, _, held = _arms_world()
    entry, _ = ep.feng_dirs.feng_directions(genes, stats_path=path)
    first = ep.gc.m1_r_arms(ctx, entry, groups, 0.2, held, held)["m1_r_gmom"]
    copied = ep.eps_config(0.02, cfg)
    assert copied["metrics"]["pseudocount_cp10k"] == 0.02 and cfg["metrics"]["pseudocount_cp10k"] == 0.5
    for eps in (0.1, 0.02):
        new, _ = ep.eps_directions(genes, eps, stats_path=path)
        moved = ep.gc.m1_r_arms(ctx, new, groups, 0.2, held, held)["m1_r_gmom"]
        assert moved.X.tobytes() != first.X.tobytes()
    assert json.dumps(cfg, sort_keys=True) == snapshot
    assert ep.feng_dirs.load_config()["metrics"]["pseudocount_cp10k"] == 0.5
    assert hashlib.sha256(ep.feng_dirs.CONFIG.read_bytes()).hexdigest() == before_file
    again_dirs, _ = ep.feng_dirs.feng_directions(genes, stats_path=path)
    assert ep.directions_sha256(again_dirs) == ep.directions_sha256(entry)
    again = ep.gc.m1_r_arms(ctx, again_dirs, groups, 0.2, held, held)["m1_r_gmom"]
    assert again.X.tobytes() == first.X.tobytes()


def test_the_groups_and_held_genes_do_not_depend_on_the_pseudocount(ep, stats):
    """Groups come from the config's label rules and min_cells_ext; held genes are R2's committed list: neither reads
    the pseudocount."""
    path, genes = stats
    labels = ["V-CM", "Endocardium", "Blood", "Foregut", "Epicardium"]
    base = ep.feng_dirs.label_groups(labels)
    assert ep.feng_dirs.label_groups(labels, ep.eps_config(0.02)) == base
    a, cells_a, _ = ep.feng_dirs.directions_from_stats(ep.dl.load_stats(path), genes, ep.eps_config(0.5))
    b, cells_b, _ = ep.feng_dirs.directions_from_stats(ep.dl.load_stats(path), genes, ep.eps_config(0.02))
    assert set(a) == set(b) and cells_a == cells_b


def test_the_direction_stats_are_ratios_from_the_moved_genes(ep):
    genes = ["a", "b", "c", "d", "e"]
    entry = {"CM": np.array([1.0, 1.0, np.nan, 2.0, 1.0])}
    new = {"CM": np.array([3.0, 1.0, np.nan, 2.0, 5.0])}
    means = {"CM": np.array([0.1, 2.0, 0.1, 0.3, 4.0])}
    out = ep.direction_stats({"entry": entry, "e01": new}, means, genes, held=["e"])
    # moved: a, b, d (c NaN, e held); entry 1 + 1 + 4 = 6; new 9 + 1 + 4 = 14; low start mean (< 0.5): a, d
    assert out["entry"]["CM"] == {"squared_size_over_entry": 1.0, "share_on_low_start_genes": pytest.approx(5 / 6)}
    assert out["e01"]["CM"] == {"squared_size_over_entry": pytest.approx(14 / 6),
                                "share_on_low_start_genes": pytest.approx(13 / 14)}
    # with the E14.5 means: low at the end are b and c (c is not moved); either end: a, b, d
    end = {"CM": np.array([3.0, 0.2, 0.1, 1.0, 0.0])}
    out = ep.direction_stats({"entry": entry, "e01": new}, means, genes, held=["e"], end=end)
    assert out["e01"]["CM"]["share_on_low_end_genes"] == pytest.approx(1 / 14)
    assert out["e01"]["CM"]["share_on_low_either_end_genes"] == pytest.approx(14 / 14)
    assert out["entry"]["CM"]["share_on_low_end_genes"] == pytest.approx(1 / 6)


def test_the_end_means_are_the_groups_feng_e14_5_means(ep, stats):
    path, genes = stats
    _, rec = ep.feng_dirs.feng_directions(genes, stats_path=path)
    means = ep.end_means(genes, rec, stats_path=path)
    st = ep.dl.load_stats(path, fields=ep.dl.STAT_FIELDS)
    for g, c in rec["groups"].items():
        np.testing.assert_array_equal(means[g], ep.dl.feng_sel(st, "E14.5", c["feng_labels"])["cp"])


def test_the_half_split_agreement_is_target_free_and_per_pseudocount(ep, stats):
    """Per pseudocount and group: the correlation of Feng's half-A and half-B directions over the moved genes, and
    eb_shrink's per-decile factors, each from feng_dirs' own functions with the pseudocount's config copy."""
    path, genes = stats
    held = genes[:3]
    out = ep.half_agreement(genes, held, {"entry": 0.5, "e002": 0.02}, stats_path=path)
    st = ep.dl.load_stats(path, fields=ep.dl.STAT_FIELDS)
    for tag, eps in (("entry", 0.5), ("e002", 0.02)):
        halves = ep.feng_dirs.direction_halves(st, genes, ep.eps_config(eps))
        dirs, _, _ = ep.feng_dirs.directions_from_stats(st, genes, ep.eps_config(eps))
        assert set(out[tag]) == set(halves) & set(dirs)
        for g, rec in out[tag].items():
            h = halves[g]
            ok = np.isfinite(h["d_a"]) & np.isfinite(h["d_b"]) & ~np.isin(genes, held)
            assert rec["correlation"] == pytest.approx(np.corrcoef(h["d_a"][ok], h["d_b"][ok])[0, 1])
            assert rec["eb_shrink_factors"] == pytest.approx(
                ep.feng_dirs.eb_shrink(dirs[g], h["d_a"], h["d_b"], h["mean_cp"])[1])
            assert len(rec["eb_shrink_factors"]) == 10
    assert ep.feng_dirs.load_config()["metrics"]["pseudocount_cp10k"] == 0.5


def test_the_start_means_are_the_groups_feng_e9_5_means(ep, stats):
    path, genes = stats
    _, rec = ep.feng_dirs.feng_directions(genes, stats_path=path)
    means = ep.start_means(genes, rec, stats_path=path)
    st = ep.dl.load_stats(path, fields=ep.dl.STAT_FIELDS)
    for g, c in rec["groups"].items():
        np.testing.assert_array_equal(means[g], ep.dl.feng_sel(st, "E9.5", c["feng_labels"])["cp"])


def test_the_largest_log_factor_and_no_clip_on_the_route(ep, stats, monkeypatch):
    """The largest |log f| a new arm applies is fraction x max |log((b + eps) / (a + eps))|, at most fraction x
    log((max mean CP10k + eps) / eps); external_shift is called without max_abs_log_factor for every arm (no clip), and
    its only other bound, the restore of a cell left empty, binds no more often for the new arm than for the entry."""
    path, genes = stats
    st = ep.dl.load_stats(path, fields=ep.dl.STAT_FIELDS)
    top = {}
    for eps in (0.5, 0.1, 0.02):
        dirs, rec = ep.eps_directions(genes, eps, stats_path=path)
        top[eps] = max(float(np.nanmax(np.abs(v))) for v in dirs.values())
        cp_max = max(float(np.max(ep.dl.feng_sel(st, s, c["feng_labels"])["cp"])) for c in rec["groups"].values()
                     for s in ("E9.5", "E14.5"))
        assert top[eps] <= np.log((cp_max + eps) / eps)
    assert top[0.02] >= top[0.1] >= top[0.5]
    calls = []
    real = t1_scope.external_shift

    def spy(*args, **kw):
        out, rep = real(*args, **kw)
        calls.append({"clip": kw.get("max_abs_log_factor"), "restored": rep.restored_cells,
                      "max_abs": max(r.get("max_abs_log_f", 0.0) for r in rep.rows)})
        return out, rep

    monkeypatch.setattr(t1_scope, "external_shift", spy)
    ctx, groups, _, held = _arms_world()
    entry, _ = ep.feng_dirs.feng_directions(genes, stats_path=path)
    new, _ = ep.eps_directions(genes, 0.02, stats_path=path)
    ep.gc.m1_r_arms(ctx, entry, groups, 0.2, held, held)
    ep.gc.m1_r_arms(ctx, new, groups, 0.4, held, held)
    assert calls and all(c["clip"] is None for c in calls)
    assert max(c["max_abs"] for c in calls[len(calls) // 2:]) > max(c["max_abs"] for c in calls[:len(calls) // 2])
    half = len(calls) // 2  # the entry's two builds, then e002's at 0.4: the empty-cell restore binds no more often
    assert sum(c["restored"] for c in calls[half:]) <= sum(c["restored"] for c in calls[:half]) == 0


# ---------------------------------------------------------------------------------------------------------------------
# The builds


def _eps_world(ep, scale=(1.7, 3.1)):
    """_arms_world with two new direction sets: each group's direction rescaled per gene (no new NaN)."""
    ctx, groups, dirs, held = _arms_world()
    rng = np.random.default_rng(9)
    eps_dirs = {tag: {g: v * (1 + s * rng.random(len(v))) for g, v in dirs.items()}
                for tag, s in zip(ep.EPS, scale, strict=True)}
    return ctx, groups, dirs, eps_dirs, held


def test_a_new_arm_differs_from_the_entrys_build_in_the_directions_and_the_fraction_only(ep, monkeypatch):
    calls = []

    def spy(ctx, directions, groups, fraction, held_m1, held_factors):
        calls.append({"ctx": ctx, "directions": directions, "groups": groups, "fraction": fraction,
                      "held_m1": held_m1, "held_factors": held_factors})
        return {"m1_r": None, "m1_r_gmom": ("built", len(calls))}

    ctx, groups, dirs, eps_dirs, held = _eps_world(ep)
    monkeypatch.setattr(ep.gc, "m1_r_arms", spy)
    preds, _, _, _ = ep.build_seed(ctx, lambda c: "b0", "pooled", dirs, eps_dirs, groups, 0.2, held, held)
    assert list(preds) == list(ep.arms_of("pooled"))
    by_arm = dict(zip([a for a in ep.arms_of("pooled") if a != "b0"], calls, strict=True))
    base = by_arm["m1_r_gmom"]
    for a in ("m1_r_gmom", *UNIFORM, *NEW, "m1_r_gmom_perm", *NEW_PERM, *SPERM):
        c = by_arm[a]
        for key in ("ctx", "groups", "held_m1", "held_factors"):
            assert c[key] is base[key], f"{a} differs from the entry's build in {key}"
    for tag, m in (("e01", 1.0), ("e01_x15", 1.5), ("e01_x2", 2.0), ("e002", 1.0), ("e002_x15", 1.5),
                   ("e002_x2", 2.0)):
        c = by_arm[f"m1_r_gmom_{tag}"]
        assert c["directions"] is eps_dirs[tag.split("_")[0]]
        assert c["fraction"] == pytest.approx(0.2 * m, abs=1e-15)
        uniform = by_arm["m1_r_gmom" if m == 1.0 else {1.5: "m1_r_gmom_x15", 2.0: "m1_r_gmom_x2"}[m]]
        assert c["fraction"] == uniform["fraction"]  # the very float the uniform arm of that multiplier gets
    assert by_arm["m1_r_gmom_e01"]["fraction"] == base["fraction"]
    for a in UNIFORM:
        assert by_arm[a]["directions"] is dirs
    # the placebos: try T46's permutation of each arm's own directions, at the _x15 fraction
    for tag, p in (("e01", NEW_PERM[0]), ("e002", NEW_PERM[1])):
        c = by_arm[p]
        want = ep.lsr.placebo_directions(ctx, eps_dirs[tag], held)
        assert all(np.array_equal(c["directions"][g], want[g], equal_nan=True) for g in want)
        assert c["fraction"] == by_arm["m1_r_gmom_x15"]["fraction"]
    perm_entry = ep.lsr.placebo_directions(ctx, dirs, held)
    assert all(np.array_equal(by_arm["m1_r_gmom_perm"]["directions"][g], perm_entry[g], equal_nan=True) for g in dirs)
    assert by_arm["m1_r_gmom_perm"]["fraction"] == 0.2
    # the expression-matched placebos: each direction permuted within scored/other x reference-level deciles
    levels = ep.reference_levels(ctx, groups, dirs)
    x15 = by_arm["m1_r_gmom_x15"]["fraction"]
    for p, d, fr in (("m1_r_gmom_sperm", dirs, 0.2), (NEW_SPERM[0], eps_dirs["e01"], x15),
                     (NEW_SPERM[1], eps_dirs["e002"], x15)):
        want = ep.splacebo_directions(ctx, d, held, levels)
        assert all(np.array_equal(by_arm[p]["directions"][g], want[g], equal_nan=True) for g in want)
        assert by_arm[p]["fraction"] == fr


def test_the_whole_cut_builds_no_x4_and_no_placebo(ep, monkeypatch):
    calls = []
    monkeypatch.setattr(ep.gc, "m1_r_arms", lambda *a: calls.append(a[3]) or {"m1_r": None, "m1_r_gmom": a[3]})
    ctx, groups, dirs, eps_dirs, held = _eps_world(ep)
    preds, _, placebo, factors = ep.build_seed(ctx, lambda c: "b0", "whole", dirs, eps_dirs, groups, 0.2, held, held)
    assert list(preds) == list(ep.arms_of("whole")) and len(calls) == 4 + 6
    assert placebo == {"sha256": {}, "same_positions_as_m1_r_gmom_perm": {}, "same_positions_as_m1_r_gmom_sperm": {}}
    assert all(set(v) == set(ep.arms_of("whole")) - {"b0"} for v in factors.values())


def test_the_rescored_arms_keep_their_output_bit_for_bit(ep):
    """b0, m1_r_gmom, the four uniform arms and m1_r_gmom_perm as build_seed makes them now equal try T46's committed
    builds on the synthetic world (test_t1_lib_step_read's golden sha256)."""
    ctx, groups, dirs, eps_dirs, held = _eps_world(ep)
    preds, log_factors, placebo, factors = ep.build_seed(ctx, harness._b0_t1(harness.B0_DAMPING), "pooled", dirs,
                                                         eps_dirs, groups, 0.2, held, held)
    got = {a: ep.lsr.prediction_sha256(preds[a]) for a in ("b0", "m1_r_gmom", *UNIFORM, "m1_r_gmom_perm")}
    assert got == {a: GOLDEN_ARMS_5BC7[a] for a in got}
    for a in NEW:
        assert ep.lsr.prediction_sha256(preds[a]) != got["m1_r_gmom"]
    for p, real in zip(NEW_PERM, ("m1_r_gmom_e01_x15", "m1_r_gmom_e002_x15"), strict=True):
        assert ep.lsr.prediction_sha256(preds[p]) != ep.lsr.prediction_sha256(preds[real])
    assert set(log_factors) == set(dirs)
    assert all(set(v) == set(ep.arms_of("pooled")) - {"b0"} for v in log_factors.values())
    for v in log_factors.values():  # the board scores every gene here: the placebo keeps the arm's factor sizes
        assert v[NEW_PERM[0]] == pytest.approx(v["m1_r_gmom_e01_x15"])
        assert v["m1_r_gmom_perm"] == pytest.approx(v["m1_r_gmom"])
    assert placebo["same_positions_as_m1_r_gmom_perm"] == {p: True for p in NEW_PERM}
    assert placebo["same_positions_as_m1_r_gmom_sperm"] == {p: True for p in NEW_SPERM}
    assert set(placebo["sha256"]) == {"m1_r_gmom_perm", *NEW_PERM, *SPERM}
    for p, real in zip(NEW_SPERM, ("m1_r_gmom_e01_x15", "m1_r_gmom_e002_x15"), strict=True):
        assert ep.lsr.prediction_sha256(preds[p]) != ep.lsr.prediction_sha256(preds[real])
        assert ep.lsr.prediction_sha256(preds[p]) != ep.lsr.prediction_sha256(preds[p.replace("_sperm", "_perm")])
    assert set(factors) == set(dirs) and all(set(v) == set(ep.arms_of("pooled")) - {"b0"} for v in factors.values())
    pb = ep.pseudobulk_changes(preds)
    assert pb["m1_r_gmom"] == 0.0 and all(pb[a] > 0 for a in NEW)


def test_with_the_entrys_directions_the_new_arms_are_the_entrys_arms_bit_for_bit(ep):
    """The idle case: new direction sets equal to the entry's make each new arm the uniform arm of its fraction (and
    the new placebo the _x15 build along the entry's own permutation) bit for bit."""
    ctx, groups, dirs, held = _arms_world()
    idle = {tag: {g: v.copy() for g, v in dirs.items()} for tag in ep.EPS}
    preds, _, placebo, _ = ep.build_seed(ctx, harness._b0_t1(harness.B0_DAMPING), "pooled", dirs, idle, groups, 0.2,
                                         held, held)
    sha = {a: ep.lsr.prediction_sha256(p) for a, p in preds.items()}
    for tag in ep.EPS:
        assert sha[f"m1_r_gmom_{tag}"] == sha["m1_r_gmom"]
        assert sha[f"m1_r_gmom_{tag}_x15"] == sha["m1_r_gmom_x15"]
        assert sha[f"m1_r_gmom_{tag}_x2"] == sha["m1_r_gmom_x2"]
    perm = ep.lsr.placebo_directions(ctx, dirs, held)
    x15_perm = ep.gc.m1_r_arms(ctx, perm, groups, 0.2 * 1.5, held, held)["m1_r_gmom"]
    assert all(sha[p] == ep.lsr.prediction_sha256(x15_perm) for p in NEW_PERM)
    assert all(placebo["sha256"][p] == placebo["sha256"]["m1_r_gmom_perm"] for p in NEW_PERM)
    sperm = ep.splacebo_directions(ctx, dirs, held, ep.reference_levels(ctx, groups, dirs))
    x15_sperm = ep.gc.m1_r_arms(ctx, sperm, groups, 0.2 * 1.5, held, held)["m1_r_gmom"]
    assert all(sha[p] == ep.lsr.prediction_sha256(x15_sperm) for p in NEW_SPERM)
    assert all(placebo["sha256"][p] == placebo["sha256"]["m1_r_gmom_sperm"] for p in NEW_SPERM)


def test_the_permutation_of_positions_is_the_entrys_when_the_moved_genes_are_the_same(ep):
    ctx, groups, dirs, eps_dirs, held = _eps_world(ep)
    assert ep.same_positions(ctx, dirs, eps_dirs["e01"], held)
    other = {g: v.copy() for g, v in eps_dirs["e01"].items()}
    other["CM"][3] = np.nan  # one moved gene fewer: another permutation
    assert not ep.same_positions(ctx, dirs, other, held)
    genes = list(ctx.inputs[ctx.task.reference_role].var_names)
    a = ep.lsr.permuted_directions(ep.positions_of(dirs), genes, held, list(ctx.board.genes), ctx.seed)
    b = ep.lsr.permuted_directions(eps_dirs["e01"], genes, held, list(ctx.board.genes), ctx.seed)
    for g in dirs:  # the permuted new direction is the new direction read at the entry's permuted positions
        idx = a[g]
        ok = np.isfinite(idx)
        np.testing.assert_array_equal(b[g][ok], eps_dirs["e01"][g][idx[ok].astype(int)])


# ---------------------------------------------------------------------------------------------------------------------
# The expression-matched placebo and the count-space records


def test_the_reference_levels_are_each_groups_mean_cp10k_on_the_reference_cells(ep):
    ctx, groups, dirs, held = _arms_world()
    dirs = {**dirs, "Nowhere": dirs["CM"]}  # a group no reference cell maps to: the mean over all cells
    levels = ep.reference_levels(ctx, groups, dirs)
    ref = ctx.inputs[ctx.task.reference_role]
    C = np.asarray(ep.core.counts_of(ref).toarray(), dtype=np.float64)
    cp = C / np.maximum(C.sum(axis=1, keepdims=True), 1.0) * 1e4
    lab = ref.obs[t1_scope.LABEL_KEY].astype(str).to_numpy()
    cell_group = np.array([groups.get(x, "") for x in lab])
    for g in ("CM", "Endocardium"):
        assert (cell_group == g).any()
        np.testing.assert_allclose(levels[g], cp[cell_group == g].mean(axis=0), rtol=1e-10)
    np.testing.assert_allclose(levels["Nowhere"], cp.mean(axis=0), rtol=1e-10)


def _strata_world(n=3000, seed=0):
    """Synthetic genes: lognormal reference levels; an entry direction of small noise; a 'decompressed' direction that
    adds large factors on the genes weakly expressed in the reference (as a small pseudocount does)."""
    rng = np.random.default_rng(seed)
    genes = [f"g{i}" for i in range(n)]
    level = rng.lognormal(0.0, 2.0, n)
    entry = rng.normal(0.0, 0.3, n)
    weak = level < np.quantile(level, 0.4)
    new = entry + np.where(weak, rng.choice([-1.0, 1.0], n) * rng.uniform(2.0, 4.0, n), 0.0)
    entry[:5] = new[:5] = np.nan
    return genes, level, entry, new, genes[5:15], genes[::2]


def test_the_expression_matched_permutation_stays_within_its_strata_and_bins(ep):
    genes, level, entry, new, held, scored = _strata_world()
    out = ep.stratified_permuted_directions({"CM": new}, genes, held, scored, [7], {"CM": level})["CM"]
    held_m = np.isin(genes, held)
    moved = np.isfinite(new) & ~held_m
    np.testing.assert_array_equal(out[~moved], new[~moved])  # NaN and held positions stay
    on = np.isin(genes, scored)
    changed = 0
    for stratum in (moved & on, moved & ~on):
        idx = np.flatnonzero(stratum)
        for part in np.array_split(idx[np.argsort(level[idx], kind="stable")], ep.LEVEL_BINS):
            assert sorted(out[part]) == sorted(new[part])  # each bin keeps its own multiset of values
            changed += int((out[part] != new[part]).sum())
    assert changed > 0
    # the same moved genes are permuted the same way whatever the values; another stream from try T46's
    again = ep.stratified_permuted_directions({"CM": entry}, genes, held, scored, [7], {"CM": level})["CM"]
    pos = ep.positions_of({"CM": new})
    pa = ep.stratified_permuted_directions(pos, genes, held, scored, [7], {"CM": level})["CM"]
    ok = np.isfinite(pa)
    np.testing.assert_array_equal(again[ok], entry[pa[ok].astype(int)])
    t46 = ep.lsr.permuted_directions({"CM": new}, genes, held, scored, [7])["CM"]
    assert not np.array_equal(t46, out, equal_nan=True)


def test_the_expression_matched_placebo_moves_about_the_counts_the_real_direction_moves(ep):
    """The reviewer's design attack: with large factors on weakly expressed genes, try T46's permutation (scored/other
    only) moves them onto well-expressed genes, so its placebo moves several times the counts the real direction moves;
    the expression-matched permutation keeps the counts moved near the real direction's."""
    genes, level, _, new, held, scored = _strata_world()

    def moved(v):
        return ep.factor_stats(v, genes, held, scored, 0.3, level)["cp10k_moved_share"]

    real = moved(new)
    t46 = [moved(ep.lsr.permuted_directions({"CM": new}, genes, held, scored, [s])["CM"]) for s in range(5)]
    strat = [moved(ep.stratified_permuted_directions({"CM": new}, genes, held, scored, [s], {"CM": level})["CM"])
             for s in range(5)]
    assert min(t46) > 2.0 * real
    assert all(abs(x / real - 1.0) < 0.25 for x in strat)


def test_the_factor_records_are_external_shifts_own_and_count_large_factors(ep, stats, monkeypatch):
    """factor_stats' max_abs_log_f is the ShiftReport's per group (no clip on the route), scored_over_log10 counts the
    scored genes beyond log 10, and the moved shares are the reference levels times |f - 1| over their sum."""
    path, genes = stats
    rows = []
    real = t1_scope.external_shift

    def spy(*args, **kw):
        out, rep = real(*args, **kw)
        rows.extend((r["type"], r["max_abs_log_f"], kw.get("max_abs_log_factor")) for r in rep.rows
                    if "max_abs_log_f" in r)
        return out, rep

    monkeypatch.setattr(t1_scope, "external_shift", spy)
    ctx, groups, _, held = _arms_world()
    new, _ = ep.eps_directions(genes, 0.02, stats_path=path)
    ep.gc.m1_r_arms(ctx, new, groups, 0.4, held, held)
    levels = ep.reference_levels(ctx, groups, new)
    scored = list(ctx.board.genes)
    assert rows and all(clip is None for _, _, clip in rows)
    for g, top, _ in rows:
        assert ep.factor_stats(new[g], genes, held, scored, 0.4, levels[g])["max_abs_log_f"] == pytest.approx(top)
    v = np.array([1.0, -8.0, np.nan, 7.0, 0.5])
    names = ["a", "b", "c", "d", "e"]
    lev = np.array([10.0, 1.0, 5.0, 2.0, 2.0])
    fs = ep.factor_stats(v, names, ["d"], ["a", "b", "c"], 0.5, lev)
    logf = np.array([0.5, -4.0, 0.0, 0.0, 0.25])
    assert fs["max_abs_log_f"] == pytest.approx(4.0) and fs["scored_over_log10"] == 1
    assert fs["cp10k_moved_share"] == pytest.approx(float((lev * np.abs(np.expm1(logf))).sum() / lev.sum()))
    assert fs["cp10k_moved_share_scored"] == pytest.approx(float((lev * np.abs(np.expm1(logf)))[:3].sum() / lev.sum()))


def test_the_pseudobulk_change_is_against_m1_r_gmoms_prediction(ep):
    from vec.evaluate import Stage

    base = Stage(X=np.log1p(np.array([[1.0, 0.0], [3.0, 2.0]], dtype=np.float32)))
    arm = Stage(X=np.log1p(np.array([[1.0, 4.0], [3.0, 0.0]], dtype=np.float32)))
    out = ep.pseudobulk_changes({"m1_r_gmom": base, "a": arm})
    want = np.abs(np.log1p([2.0, 2.0]) - np.log1p([2.0, 1.0])).mean()
    assert out["m1_r_gmom"] == 0.0 and out["a"] == pytest.approx(want, rel=1e-6)


# ---------------------------------------------------------------------------------------------------------------------
# The reading


ROWS = {"pooled": 7919, "lib1": 1777, "lib2": 2711, "lib3": 3431, "whole": 9127}
CONTAM = {"pooled": 613, "lib1": 401, "lib2": 199, "lib3": 13, "whole": 0}
DIRS = {"entry": {"pseudocount": 0.5, "sha256": {"CM": "e-cm", "Endocardium": "e-en"}},
        "e01": {"pseudocount": 0.1, "sha256": {"CM": "a-cm", "Endocardium": "a-en"}},
        "e002": {"pseudocount": 0.02, "sha256": {"CM": "b-cm", "Endocardium": "b-en"}}}
# a per-cut read-back that matches the quoted library-3 sizes (m1_r_gmom at 0.58 of its nominal 0.2)
RULER = {"m1_r_gmom": 0.116, "m1_r_gmom_x15": 0.175, "m1_r_gmom_x2": 0.233, "m1_r_gmom_x3": 0.349,
         "m1_r_gmom_x4": 0.463, "m1_r_gmom_e01": 0.173, "m1_r_gmom_e01_x15": 0.2595, "m1_r_gmom_e01_x2": 0.347,
         "m1_r_gmom_e002": 0.207, "m1_r_gmom_e002_x15": 0.312, "m1_r_gmom_e002_x2": 0.416}


def _task(ep, truth, gains=None, per_seed=None, seeds=range(10), commit="c" * 40, clean=True, base_shift=0.0,
          identical=(), sha=None, fractions=None, pseudocounts=None, directions=None, drop=()):
    """A task file. gains: arm -> gain over m1_r_gmom on every seed; per_seed: arm -> ten gains. An arm with no gain
    given gets a distinct tiny offset (+d on seed 0, -d on seed 1); `identical` names arms that copy m1_r_gmom's
    record."""
    gains, per_seed = gains or {}, per_seed or {}
    arms_list = [a for a in ep.arms_of(truth) if a not in drop]
    rows = []
    for i, s in enumerate(seeds):
        arms = {"b0": _arm(B0), "m1_r_gmom": _arm(BASE + base_shift)}
        for j, a in enumerate(arms_list):
            if a in arms:
                continue
            if a in per_seed:
                g = per_seed[a][i]
            elif a in gains:
                g = gains[a]
            else:
                d = 1e-9 * (j + 1)
                g = d if i == 0 else -d if i == 1 else 0.0
            arms[a] = _arm(BASE + base_shift + g)
        for a in identical:
            arms[a] = dict(arms["m1_r_gmom"])
        lf = {g: {a: 0.1 for a in arms_list if a != "b0"} for g in ("CM", "Endocardium")}
        fs = {g: {a: {"max_abs_log_f": 0.5, "scored_over_log10": 0, "cp10k_moved_share": 0.2,
                      "cp10k_moved_share_scored": 0.1, "cp10k_moved_share_scored_weak": 0.03,
                      "cp10k_moved_share_scored_strong": 0.07} for a in arms_list if a != "b0"}
              for g in ("CM", "Endocardium")}
        rb = {g: {a: RULER[a] for a in ep.candidates(truth)} for g in ("CM", "Endocardium")}
        rows.append({"seed": s, "floor": 50.0, "arms": arms,
                     "sha256": {a: (sha or {}).get(a, f"{a}-{s}") for a in arms_list},
                     "sparsity": {a: {"zero_share": 0.25, "detected_per_cell": 6.0} for a in arms_list},
                     "scored_log_factors": lf, "factor_stats": fs, "ruler_readback": rb,
                     "pseudobulk_change": {a: 0.0 if a == "m1_r_gmom" else 0.05 for a in arms_list},
                     "placebo": {"same_positions_as_m1_r_gmom_perm": {p: True for p in NEW_PERM},
                                 "same_positions_as_m1_r_gmom_sperm": {p: True for p in NEW_SPERM}}})
    spec = ep.arm_spec(truth)
    return {"label": "synthetic", "task": ep.TASK_OF[truth], "truth": truth, "truth_rows": ROWS[truth],
            "truth_contamination": CONTAM[truth], "fraction": 0.2,
            "fractions": fractions or {a: v["fraction"] for a, v in spec.items()},
            "pseudocounts": pseudocounts or {a: v["pseudocount"] for a, v in spec.items()},
            "directions": json.loads(json.dumps(directions or DIRS)),
            "direction_stats": {t: {"CM": {"squared_size_over_entry": 1.0, "share_on_low_start_genes": 0.4,
                                           "share_on_low_end_genes": 0.3, "share_on_low_either_end_genes": 0.6}}
                                for t in DIRS},
            "half_agreement": {t: {"CM": {"correlation": 0.8, "eb_shrink_factors": [0.5] * 10}} for t in DIRS},
            "external": {t: {"pseudocount_cp10k": d["pseudocount"]} for t, d in DIRS.items()},
            "seeds": list(seeds), "rows": rows,
            "stamp": {"commit": commit, "tree_clean": clean, "utc": "x", "inputs": {}}}


def _write(ep, tmp_path, truths, committed_shift=None, committed_sha=None):
    """truths: truth -> task record; the committed files (T46's per chamber truth, T44's whole cut) hold the rescored
    arms' scores from the same rows unless shifted (truth -> arm -> shift), and on the chamber truths their sha256 from
    the same rows unless replaced (truth -> arm -> sha256 on every seed; T44's whole-cut file holds no sha256)."""
    committed_shift, committed_sha = committed_shift or {}, committed_sha or {}
    out = tmp_path / "out"
    (out / "tasks").mkdir(parents=True)
    for k, rec in truths.items():
        (out / "tasks" / f"{k}.json").write_text(json.dumps(rec))
    t46 = tmp_path / "t46"
    t46.mkdir()
    t44 = tmp_path / "t44_whole.json"
    for k in TRUTHS:
        repro = ep.REPRO_WHOLE if k == "whole" else ep.REPRO_CHAMBER
        src = truths.get(k) or _task(ep, k)
        rows = [{"seed": r["seed"], "arms": {a: {**r["arms"][a], "score": r["arms"][a]["score"]
                                                 + committed_shift.get(k, {}).get(a, 0.0)} for a in repro}}
                for r in src["rows"]]
        if k != "whole":
            for row, r in zip(rows, src["rows"], strict=True):
                row["sha256"] = {a: committed_sha.get(k, {}).get(a, r["sha256"][a]) for a in repro}
        (t44 if k == "whole" else t46 / f"{k}.json").write_text(json.dumps({"rows": rows}))
    return out, t46, t44


def _read(ep, tmp_path, by_truth=None, over=None, committed_shift=None, edit=None, committed_sha=None):
    """by_truth: truth -> gains (both cuts: give 'pooled' and 'whole'); over: truth -> _task overrides (None drops
    the truth); edit: truth -> a function changing the record after it is made (the committed files are written
    first, from the unedited rows)."""
    by_truth, over = by_truth or {}, over or {}
    truths = {k: _task(ep, k, gains=by_truth.get(k), **(over.get(k) or {})) for k in TRUTHS
              if not (k in over and over[k] is None)}
    out, t46, t44 = _write(ep, tmp_path, {k: v for k, v in truths.items()}, committed_shift, committed_sha)
    for k, f in (edit or {}).items():
        f(truths[k])
        (out / "tasks" / f"{k}.json").write_text(json.dumps(truths[k]))
    return ep.read(out, t46, t44)["reading"]


def _both(gains):
    return {"pooled": gains, "whole": gains}


SMALLER = "the smaller pseudocount scores better at the same delivered size"
HEAD2 = "Statement 2, decompression against uniform inflation at about the same delivered size: "


@pytest.mark.parametrize("gains, statement", [
    ({"m1_r_gmom_e01": 1.0, "m1_r_gmom_e01_x2": 1.0}, SMALLER),
    ({"m1_r_gmom_e01": -1.0, "m1_r_gmom_e01_x2": -1.0}, "uniform inflation scores better"),
    ({}, "no difference at the bar"),
    ({"m1_r_gmom_e01": 0.01, "m1_r_gmom_e01_x2": -0.2}, "no difference at the bar"),
    ({"m1_r_gmom_e01": 1.0, "m1_r_gmom_e01_x2": 0.01}, "mixed: passes at the smaller size, neither at the larger"),
    ({"m1_r_gmom_e01": 1.0, "m1_r_gmom_e01_x2": -1.0}, "mixed: passes at the smaller size, loses at the larger"),
    ({"m1_r_gmom_e01": -0.1, "m1_r_gmom_e01_x2": -1.0}, "mixed: neither at the smaller size, loses at the larger"),
    ({"m1_r_gmom_x15": 1.0, "m1_r_gmom_e01": 2.0, "m1_r_gmom_x3": 2.0, "m1_r_gmom_e01_x2": 1.0},
     "mixed: passes at the smaller size, loses at the larger"),
    ({"m1_r_gmom_x15": 1.0, "m1_r_gmom_e01": 2.0, "m1_r_gmom_x3": 2.0, "m1_r_gmom_e01_x2": 3.0}, SMALLER),
])
def test_statement_two(ep, tmp_path, gains, statement):
    r = _read(ep, tmp_path, _both(gains))
    assert r["valid"] and r["outcome"] == statement == r["statement_2"]["statement"]
    c1, c2 = r["statement_2"]["comparisons"]
    assert (c1["two_cuts"]["change"], c1["two_cuts"]["base"]) == ("m1_r_gmom_e01", "m1_r_gmom_x15")
    assert (c2["two_cuts"]["change"], c2["two_cuts"]["base"]) == ("m1_r_gmom_e01_x2", "m1_r_gmom_x3")
    if statement == "no difference at the bar":
        assert c1["two_cuts"]["label"] == c2["two_cuts"]["label"] == "neither"
    txt = (tmp_path / "out" / "t52.txt").read_text()
    assert HEAD2 + statement in txt
    assert "m1_r_gmom_e01 against m1_r_gmom_x15 on the two cuts: " in txt
    assert "m1_r_gmom_e01_x2 against m1_r_gmom_x3 on the two cuts: " in txt
    for k in ("pooled", "whole"):  # each pair's count-space sizes beside it, on each cut
        assert f"count-space sizes on {k} (ruler-matched is not count-matched): m1_r_gmom_e01: ruler size 0.173" in txt
    assert not _prints_a_count(txt, (*ROWS.values(), *[c for c in CONTAM.values() if c]))


def test_statement_two_on_the_library_truths(ep, tmp_path):
    r = _read(ep, tmp_path, {"lib1": {"m1_r_gmom_e01": 1.0}, "lib2": {"m1_r_gmom_x15": 1.0},
                             "lib3": {"m1_r_gmom_e01_x2": 0.2}})
    libs = r["statement_2"]["comparisons"][0]["libraries"]
    assert {k: v["label"] for k, v in libs.items()} == {"lib1": "ahead", "lib2": "behind", "lib3": "level"}
    assert r["statement_2"]["comparisons"][1]["libraries"]["lib3"]["label"] == "level"
    assert r["outcome"] == "no difference at the bar"


@pytest.mark.parametrize("gain, cuts, label", [
    (1.0, ("pooled", "whole"), "passes"), (-1.0, ("pooled", "whole"), "loses"), (1.0, ("pooled",), "neither"),
    (0.3, ("pooled", "whole"), "neither"), (0.31, ("pooled", "whole"), "passes")])
def test_statement_one(ep, tmp_path, gain, cuts, label):
    by = {k: {"m1_r_gmom_e002_x15": gain} for k in cuts}
    by.update({"lib1": {"m1_r_gmom_e002_x15": 1.0}, "lib2": {"m1_r_gmom_e002_x15": -1.0}})
    r = _read(ep, tmp_path, by)
    s = r["statement_1"]["m1_r_gmom_e002_x15"]
    assert s["two_cuts"]["label"] == label
    assert {k: v["label"] for k, v in s["libraries"].items()} == {"lib1": "ahead", "lib2": "behind", "lib3": "level"}
    assert all(r["statement_1"][a]["two_cuts"]["label"] == "neither" for a in NEW if a != "m1_r_gmom_e002_x15")
    txt = (tmp_path / "out" / "t52.txt").read_text()
    assert f"m1_r_gmom_e002_x15 against m1_r_gmom on the two cuts: {label}" in txt


def test_a_gain_of_exactly_point_three_does_not_pass(ep, tmp_path):
    r = _read(ep, tmp_path, _both({"m1_r_gmom_e01": 0.3}), over={k: {"base_shift": -BASE} for k in TRUTHS})
    s = r["statement_1"]["m1_r_gmom_e01"]
    assert r["valid"] and s["two_cuts"]["label"] == "neither"
    assert s["two_cuts"]["cuts"]["pooled"]["mean"] == pytest.approx(0.3) and s["two_cuts"]["cuts"]["pooled"][
        "seeds_ahead"] == 10
    r = _read(ep, tmp_path / "b", {"lib2": {"m1_r_gmom_e01": 0.3}}, over={k: {"base_shift": -BASE} for k in TRUTHS})
    assert r["statement_1"]["m1_r_gmom_e01"]["libraries"]["lib2"]["label"] == "level"


@pytest.mark.parametrize("gains, label", [
    ({"m1_r_gmom_e01_x15": 1.0, "m1_r_gmom_e01_x15_sperm": 0.0},
     "the gain needs gene identity (not separated from the entry's own direction-specific gain)"),
    ({"m1_r_gmom_e01_x15": 0.01, "m1_r_gmom_e01_x15_sperm": -2.0}, "no gain on this truth"),
    ({"m1_r_gmom_e01_x15": 1.0, "m1_r_gmom_e01_x15_sperm": 0.9}, "not separated"),
    ({"m1_r_gmom_e01_x15": 1.0, "m1_r_gmom_e01_x15_sperm": 1.5}, "not separated"),
    # try T46's placebo does far more damage (D with it is ahead) but the expression-matched one does not: not read
    ({"m1_r_gmom_e01_x15": 1.0, "m1_r_gmom_e01_x15_sperm": 0.9, "m1_r_gmom_e01_x15_perm": -3.0}, "not separated"),
])
def test_statement_three(ep, tmp_path, gains, label):
    r = _read(ep, tmp_path, {"lib3": gains})
    s = r["statement_3"]["lib3"]["e01"]
    assert s["label"] == label and s["same_positions"] is True
    d = s["did"]
    assert d["quantity"] == "(m1_r_gmom_e01_x15 - m1_r_gmom) - (m1_r_gmom_e01_x15_sperm - m1_r_gmom_sperm)"
    assert d["mean"] == pytest.approx(gains["m1_r_gmom_e01_x15"] - gains["m1_r_gmom_e01_x15_sperm"], abs=1e-7)
    t46 = s["did_t46_placebo_not_read"]
    assert t46["quantity"] == "(m1_r_gmom_e01_x15 - m1_r_gmom) - (m1_r_gmom_e01_x15_perm - m1_r_gmom_perm)"
    assert t46["mean"] == pytest.approx(gains["m1_r_gmom_e01_x15"] - gains.get("m1_r_gmom_e01_x15_perm", 0.0),
                                        abs=1e-7)
    txt = (tmp_path / "out" / "t52.txt").read_text()
    # every statement-3 label carries the caveat on its own line: D checks generic movement, not which genes
    assert f"on lib3, e01: {label} ({ep.S3_CAVEAT}); gate: " in txt and "reported, not read (try T46's placebo" in txt
    assert all(ep.S3_CAVEAT in line for line in txt.splitlines() if re.match(r"\s+on \w+, e0\w+: (?!\()", line)
               and "gate: " in line)
    assert r["statement_3"]["lib1"]["e01"]["label"] == "no gain on this truth"
    assert r["statement_3"]["lib3"]["e002"]["label"] == "no gain on this truth"
    assert "whole" not in r["statement_3"]


def test_d_uses_the_entrys_placebo_as_the_placebo_base(ep, tmp_path):
    r = _read(ep, tmp_path, {"lib2": {"m1_r_gmom_e002_x15": 1.0, "m1_r_gmom_e002_x15_sperm": 0.5,
                                      "m1_r_gmom_sperm": -0.5}})
    s = r["statement_3"]["lib2"]["e002"]
    assert s["did"]["mean"] == pytest.approx(1.0 - (0.5 + 0.5)) and s["label"] == "not separated"


def _sperm_differs(rec):
    rec["rows"][4]["placebo"]["same_positions_as_m1_r_gmom_sperm"]["m1_r_gmom_e002_x15_sperm"] = False


def test_statement_three_is_not_read_where_the_permutation_differs(ep, tmp_path):
    r = _read(ep, tmp_path, {"lib1": {"m1_r_gmom_e002_x15": 1.0, "m1_r_gmom_e002_x15_sperm": -1.0}},
              edit={"lib1": _sperm_differs})
    assert r["valid"]
    s = r["statement_3"]["lib1"]
    assert s["e002"]["label"] == "not read (permutation differs)" and s["e002"]["same_positions"] is False
    assert s["e01"]["label"] == "no gain on this truth"
    assert r["same_positions"]["lib1"]["m1_r_gmom_e002_x15_sperm"] is False
    assert r["same_positions"]["lib1"]["m1_r_gmom_e002_x15_perm"] is True
    txt = (tmp_path / "out" / "t52.txt").read_text()
    assert "on lib1, e002: not read (permutation differs)" in txt


def test_the_brackets_of_the_0_02_arms_and_the_reported_numbers(ep, tmp_path):
    by = {k: {"m1_r_gmom_e002_x2": 1.0, "m1_r_gmom_x4": 2.0, "m1_r_gmom_e01_x2": 0.5} for k in TRUTHS}
    r = _read(ep, tmp_path, by)
    b = r["brackets_0.02"]
    assert set(b) == {"m1_r_gmom_e002", "m1_r_gmom_e002_x15", "m1_r_gmom_e002_x2"}
    assert set(b["m1_r_gmom_e002"]) == {"m1_r_gmom_x15", "m1_r_gmom_x2"}
    assert b["m1_r_gmom_e002_x2"]["m1_r_gmom_x4"]["two_cuts"] is None  # x4 is not built on the whole cut
    assert b["m1_r_gmom_e002_x2"]["m1_r_gmom_x4"]["libraries"]["lib1"]["label"] == "behind"
    assert b["m1_r_gmom_e002_x2"]["m1_r_gmom_x3"]["two_cuts"]["label"] == "passes"
    assert r["highest_mean"] == {**{k: "m1_r_gmom_x4" for k in CHAMBER}, "whole": "m1_r_gmom_e002_x2"}
    nd = r["new_arm_differences"]
    assert nd["lib1"]["m1_r_gmom_e002_x2"]["nearest_uniform"] == "m1_r_gmom_x4"
    assert nd["whole"]["m1_r_gmom_e002_x2"]["nearest_uniform"] == "m1_r_gmom_x3"
    assert {a: nd["lib1"][a]["nearest_uniform"] for a in NEW} == {
        "m1_r_gmom_e01": "m1_r_gmom_x15", "m1_r_gmom_e01_x15": "m1_r_gmom_x2", "m1_r_gmom_e01_x2": "m1_r_gmom_x3",
        "m1_r_gmom_e002": "m1_r_gmom_x2", "m1_r_gmom_e002_x15": "m1_r_gmom_x3", "m1_r_gmom_e002_x2": "m1_r_gmom_x4"}
    c = nd["lib1"]["m1_r_gmom_e01_x2"]["against_m1_r_gmom"]
    assert c["mean"] == pytest.approx(0.5) and c["seeds_ahead"] == 10 and c["sd"] == pytest.approx(0.0, abs=1e-6)
    assert c["points"]["de"] + c["points"]["mmd_u_variogram"] == pytest.approx(0.5)
    assert set(r["sparsity"]["whole"]) == set(ep.arms_of("whole"))
    txt = (tmp_path / "out" / "t52.txt").read_text()
    assert ("m1_r_gmom_e002_x2 against m1_r_gmom_x4 on the two cuts: not made (m1_r_gmom_x4 is not built on "
            "whole)") in txt
    assert "truth whole: highest mean m1_r_gmom_e002_x2" in txt
    # the pairing that differs between the truths is named; the extrapolated size is marked; the count-space records
    assert r["nearest_uniform_differs"] == {"m1_r_gmom_e002_x2": {**{k: "m1_r_gmom_x4" for k in CHAMBER},
                                                                  "whole": "m1_r_gmom_x3"}}
    assert r["extrapolated_delivered"] == {"m1_r_gmom_e01_x15": pytest.approx(0.2595)}
    assert "differs between the truths for: m1_r_gmom_e002_x2 (pooled m1_r_gmom_x4" in txt
    assert "Extrapolated delivered sizes (not quoted read-backs): m1_r_gmom_e01_x15 0.2595" in txt
    assert "nearest delivered size m1_r_gmom_x2 (the arm's size is extrapolated)" in txt
    assert "Statement 1, each new arm against m1_r_gmom (does not separate decompression from step size):" in txt
    cs = r["count_space"]["lib2"]["m1_r_gmom_e002_x15_sperm"]
    assert cs["pseudobulk_change"] == pytest.approx(0.05) and cs["groups"]["CM"]["cp10k_moved_share"] == 0.2
    assert set(r["count_space"]["whole"]) == set(ep.arms_of("whole")) - {"b0"}
    assert "max 0.500 over-log10 0 moved 0.200/0.100" in txt and "Feng half-A against half-B directions" in txt
    for line in ep.LIMITS:
        assert f"Limit: {line}" in txt
    assert not _prints_a_count(txt, (*ROWS.values(), *[c for c in CONTAM.values() if c]))


def test_the_limits_are_in_the_docstring(ep):
    doc = " ".join(ep.__doc__.split())
    for line in ep.LIMITS:
        assert f"- {line}" in doc


def _same_dir(rec):
    rec["directions"]["e01"]["sha256"]["Endocardium"] = rec["directions"]["entry"]["sha256"]["Endocardium"]


def _other_groups(rec):
    rec["directions"]["e002"]["sha256"].pop("CM")


def _no_sha(rec):
    del rec["rows"][3]["sha256"]["m1_r_gmom_e01"]


def _no_sparsity(rec):
    del rec["rows"][2]["sparsity"]["m1_r_gmom_e002"]


def _anchor_failure(rec):
    rec["anchor_failures"] = [{"seed": 3, "reason": "inverted"}]


def _no_factor_stats(rec):
    del rec["rows"][1]["factor_stats"]["CM"]["m1_r_gmom_sperm"]


def _no_pseudobulk(rec):
    del rec["rows"][0]["pseudobulk_change"]["m1_r_gmom_e002"]


def _no_rows(rec):
    del rec["rows"]


def _ext_wrong(rec):  # e01 labelled 0.1 by the script's constants but made with 0.02 (the reviewer's mutant M14)
    rec["external"]["e01"]["pseudocount_cp10k"] = 0.02


def _ext_missing(rec):
    del rec["external"]


def _no_readback(rec):
    del rec["rows"][6]["ruler_readback"]["Endocardium"]["m1_r_gmom_e01_x2"]


@pytest.mark.parametrize("truth, arm", [("pooled", "m1_r_gmom"), ("lib2", "m1_r_gmom_perm"), ("lib3", "m1_r_gmom_x4"),
                                        ("lib1", "b0")])
def test_a_rescored_arm_with_another_sha256_than_t46s_is_void(ep, tmp_path, truth, arm):
    """A scorer-free check: on the chamber truths the rescored arms' prediction sha256 must equal T46's recorded ones
    (even with the scores equal); T44's whole-cut file holds none and is checked by score only."""
    r = _read(ep, tmp_path, _both({"m1_r_gmom_e01": 1.0}), committed_sha={truth: {arm: "other"}})
    assert not r["valid"] and r["outcome"] == "VOID"
    assert any(f"{truth}: the rescored arms' predictions do not have T46's recorded sha256" in x
               for x in r["validity"]["reasons"]), r["validity"]["reasons"]
    assert r["validity"]["reproduction"][truth]["scores_hold"] and not r["validity"]["reproduction"][truth]["holds"]
    assert r["validity"]["reproduction"]["whole"]["sha256_holds"]


@pytest.mark.parametrize("over, edit, shift, reason", [
    ({"lib2": None}, None, None, "missing truths ['lib2']"),
    ({"whole": None}, None, None, "missing truths ['whole']"),
    ({"lib1": {"seeds": range(9)}}, None, None, "lib1: seeds"),
    (None, None, {"pooled": {"m1_r_gmom_perm": 1e-6}}, "pooled: b0, m1_r_gmom"),
    (None, None, {"lib3": {"m1_r_gmom_x4": 1e-6}}, "lib3: b0, m1_r_gmom"),
    (None, None, {"whole": {"m1_r_gmom_x2": 1e-6}}, "whole: b0, m1_r_gmom"),
    (None, None, {"whole": {"b0": 1e-6}}, "whole: b0, m1_r_gmom"),
    ({"lib3": {"fractions": {"m1_r_gmom_e01_x15": 0.2}}}, None, None, "wrong fraction on"),
    ({"whole": {"pseudocounts": {"m1_r_gmom_e002": 0.1}}}, None, None, "wrong pseudocount on"),
    ({"lib1": {"commit": "d" * 40}}, None, None, "2 commits"),
    ({"lib2": {"clean": False}}, None, None, "tree not clean on ['lib2']"),
    ({"lib2": {"clean": None}}, None, None, "tree not clean on ['lib2']"),
    ({"pooled": {"identical": ("m1_r_gmom_e01_x2",)}}, None, None, "new arms score exactly as m1_r_gmom"),
    ({"lib1": {"sha": {"m1_r_gmom_e01_x15_perm": "s", "m1_r_gmom_e01_x15": "s"}}}, None, None,
     "the placebo m1_r_gmom_e01_x15_perm's predictions have m1_r_gmom_e01_x15's sha256"),
    ({"lib3": {"sha": {"m1_r_gmom_e002_x15_perm": "s", "m1_r_gmom_e002_x15": "s"}}}, None, None,
     "the placebo m1_r_gmom_e002_x15_perm's predictions have m1_r_gmom_e002_x15's sha256"),
    (None, {"lib2": _same_dir}, None, "lib2: the e01 directions do not differ from the entry's in every group"),
    (None, {"whole": _other_groups}, None, "whole: the e002 directions do not differ"),
    ({"lib1": {"drop": ("m1_r_gmom_e002_x2",)}}, None, None, "lib1: arms ['m1_r_gmom_e002_x2'] missing"),
    (None, {"lib3": _no_sha}, None, "lib3: sha256 not recorded for every arm"),
    (None, {"pooled": _no_sparsity}, None, "pooled: the sparsity or scored log-factor records are missing"),
    (None, {"lib1": _anchor_failure}, None, "anchors invert or go missing on seeds [3]"),
    (None, {"lib2": _no_factor_stats}, None, "lib2: the count-space factor or pseudobulk records are missing"),
    (None, {"whole": _no_pseudobulk}, None, "whole: the count-space factor or pseudobulk records are missing"),
    ({"lib3": {"sha": {"m1_r_gmom_sperm": "s", "m1_r_gmom": "s"}}}, None, None,
     "the placebo m1_r_gmom_sperm's predictions have m1_r_gmom's sha256"),
    ({"pooled": {"sha": {"m1_r_gmom_e01_x15_sperm": "s", "m1_r_gmom_e01_x15": "s"}}}, None, None,
     "the placebo m1_r_gmom_e01_x15_sperm's predictions have m1_r_gmom_e01_x15's sha256"),
    (None, {"lib2": _no_rows}, None, "lib2: seeds [], not"),
    (None, {"pooled": _ext_wrong}, None,
     "pooled: wrong pseudocount in the direction code's own record (external) for ['e01']"),
    (None, {"whole": _ext_missing}, None,
     "whole: wrong pseudocount in the direction code's own record (external) for ['entry', 'e01', 'e002']"),
    (None, {"lib3": _no_readback}, None, "lib3: the per-cut ruler read-back records are missing"),
])
def test_the_void_cases(ep, tmp_path, over, edit, shift, reason):
    r = _read(ep, tmp_path, _both({"m1_r_gmom_e01": 1.0}), over=over, edit=edit, committed_shift=shift)
    assert not r["valid"] and r["outcome"] == "VOID"
    assert any(reason in x for x in r["validity"]["reasons"]), r["validity"]["reasons"]
    txt = (tmp_path / "out" / "t52.txt").read_text()
    assert "Void: " in txt and not _prints_a_count(txt, (*ROWS.values(), *[c for c in CONTAM.values() if c]))


def test_a_valid_table_is_valid_and_a_shift_within_tolerance_too(ep, tmp_path):
    assert _read(ep, tmp_path, _both({"m1_r_gmom_e01": 1.0}))["valid"]
    r = _read(ep, tmp_path / "b", committed_shift={"lib1": {"m1_r_gmom_x3": 1e-12}})
    assert r["valid"]


def test_the_reading_refuses_to_overwrite(ep, tmp_path):
    truths = {k: _task(ep, k) for k in TRUTHS}
    out, t46, t44 = _write(ep, tmp_path, truths)
    ep.read(out, t46, t44)
    assert (out / "t52.json").exists() and (out / "t52.txt").exists()
    with pytest.raises(SystemExit, match="already holds"):
        ep.read(out, t46, t44)


# ---------------------------------------------------------------------------------------------------------------------
# The run's refusals


def test_the_run_needs_its_row(ep, monkeypatch):
    monkeypatch.setattr(ep.n12, "tries_rows", lambda: {"T46", "T51"})
    monkeypatch.setattr(ep.gc, "n22_passed", lambda path: True)
    with pytest.raises(SystemExit, match="T52 has no row"):
        ep.main(["run", "--truth", "lib3", "--out", "t52"])


def test_the_run_needs_n22(ep, monkeypatch):
    monkeypatch.setattr(ep.n12, "tries_rows", lambda: {"T52"})
    with pytest.raises(SystemExit, match="does not hold a passing N22"):
        ep.main(["run", "--truth", "pooled", "--out", "t52", "--n22", "n22x"])


def test_the_run_needs_the_entrys_pseudocount_in_the_config(ep, monkeypatch):
    monkeypatch.setattr(ep.n12, "tries_rows", lambda: {"T52"})
    monkeypatch.setattr(ep.gc, "n22_passed", lambda path: True)
    monkeypatch.setattr(ep.feng_dirs, "load_config", lambda: ep.eps_config(0.1, {"metrics": {}}))
    with pytest.raises(SystemExit, match="not the entry's 0.5"):
        ep.main(["run", "--truth", "pooled", "--out", "t52"])


@pytest.mark.parametrize("truth", TRUTHS)
def test_the_run_refuses_an_existing_task_file(ep, monkeypatch, truth):
    monkeypatch.setattr(ep.n12, "tries_rows", lambda: {"T52"})
    monkeypatch.setattr(ep.gc, "n22_passed", lambda path: True)
    tasks = ep.RESULTS / "t52" / "tasks"
    tasks.mkdir(parents=True)
    (tasks / f"{truth}.json").write_text("{}")
    with pytest.raises(SystemExit, match="already holds"):
        ep.main(["run", "--truth", truth, "--out", "t52"])


def test_the_run_refuses_an_unknown_truth(ep):
    with pytest.raises(SystemExit):
        ep.main(["run", "--truth", "lib4", "--out", "t52"])


# ---------------------------------------------------------------------------------------------------------------------
# Review fixes (campaign/t1-q-fix2)


def test_every_arm_of_a_seed_shares_one_cell_draw_and_the_limit_says_so(ep):
    """The pairs share the cell draw (cells, composition, count) and differ only in the count redraw: the limit that
    said every pair differed by a full redraw was false."""
    ctx, groups, dirs, eps_dirs, held = _eps_world(ep)
    preds, _, _, _ = ep.build_seed(ctx, harness._b0_t1(harness.B0_DAMPING), "pooled", dirs, eps_dirs, groups, 0.2,
                                   held, held)
    directional = [a for a in ep.arms_of("pooled") if a != "b0"]
    first = preds["m1_r_gmom"]
    for a in directional:
        assert preds[a].n == first.n and np.array_equal(preds[a].celltype, first.celltype), a
    assert not np.array_equal(preds["m1_r_gmom_e01"].X, first.X)  # the count redraw differs
    assert not any("full redraw" in x for x in ep.LIMITS)
    assert ("All arms of a seed share one cell draw (the same cells, composition and count); they differ in the count "
            "redraw of the direction step and of g_mom's factors (the same streams with other factors), so the "
            "cell-draw noise is shared within a pair and the redraw noise is not.") in ep.LIMITS


def test_the_expression_matched_placebos_flag_uses_the_expression_matched_permutation(ep, monkeypatch):
    """Mutant M6: build_seed must compute the _sperm placebos' same-position flag with the expression-matched
    permutation, and the _perm placebos' with try T46's."""
    calls = []
    real = ep.same_positions

    def spy(ctx, a, b, held, permute=None):
        calls.append(permute)
        return real(ctx, a, b, held, permute=permute)

    monkeypatch.setattr(ep, "same_positions", spy)
    monkeypatch.setattr(ep.gc, "m1_r_arms", lambda *a: {"m1_r": None, "m1_r_gmom": a[3]})
    ctx, groups, dirs, eps_dirs, held = _eps_world(ep)
    ep.build_seed(ctx, lambda c: "b0", "pooled", dirs, eps_dirs, groups, 0.2, held, held)
    assert len(calls) == 4 and calls[0] is None and calls[1] is None  # the two _perm placebos: T46's permutation
    levels = ep.reference_levels(ctx, groups, dirs)
    pos = ep.positions_of(dirs)
    want = ep.splacebo_directions(ctx, pos, held, levels)
    t46 = ep.lsr.placebo_directions(ctx, pos, held)
    for permute in calls[2:]:
        assert permute is not None
        got = permute(ctx, pos, held)
        assert all(np.array_equal(got[g], want[g], equal_nan=True) for g in want)
        assert not all(np.array_equal(got[g], t46[g], equal_nan=True) for g in t46)


def test_the_moved_shares_split_into_weak_and_strong_scored_genes(ep):
    v = np.array([1.0, -2.0, 0.5, 1.0, 3.0, np.nan])
    names = ["a", "b", "c", "d", "e", "f"]
    lev = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
    fs = ep.factor_stats(v, names, [], ["a", "b", "c", "d"], 0.5, lev)
    moved = lev * np.abs(np.expm1(0.5 * np.nan_to_num(v)))
    # scored levels 1, 2, 3, 4: median 2.5, weak a and b, strong c and d
    assert fs["cp10k_moved_share_scored_weak"] == pytest.approx(moved[:2].sum() / lev.sum())
    assert fs["cp10k_moved_share_scored_strong"] == pytest.approx(moved[2:4].sum() / lev.sum())
    assert fs["cp10k_moved_share_scored"] == pytest.approx(fs["cp10k_moved_share_scored_weak"]
                                                           + fs["cp10k_moved_share_scored_strong"])
    assert fs["cp10k_moved_share"] == pytest.approx(moved.sum() / lev.sum())


def test_the_ruler_readback_is_the_entrys_ruler_on_the_cuts_reference_cells(ep):
    """Target-free: each arm's direction step alone on the reference cells, one fixed stream; the slope of the
    realised log ratio (pseudocount 0.5) on the entry's direction. An idle direction reads 0, a larger fraction reads
    more, the stream is fixed, and the value is the one computed by hand from external_shift."""
    ctx, groups, dirs, eps_dirs, held = _eps_world(ep)
    zero = {g: np.zeros_like(v) for g, v in dirs.items()}
    arms = {"m1_r_gmom": (dirs, 0.2), "m1_r_gmom_x2": (dirs, 0.4), "idle": (zero, 0.2), "e01": (eps_dirs["e01"], 0.2)}
    rb = ep.ruler_readback(ctx, groups, dirs, arms, held)
    assert set(rb) == {"CM", "Endocardium"}
    for by in rb.values():
        assert by["idle"] == 0.0
        assert 0.0 < by["m1_r_gmom"] < by["m1_r_gmom_x2"]
    again = ep.ruler_readback(ctx, groups, dirs, {"m1_r_gmom": (dirs, 0.2)}, held)
    assert again["CM"]["m1_r_gmom"] == rb["CM"]["m1_r_gmom"]
    ref = ctx.inputs[ctx.task.reference_role]
    moved, _ = t1_scope.external_shift(ref, dirs, 0.2, groups, held, seed=[ctx.seed, ep.READBACK_STREAM])
    lab = ref.obs[t1_scope.LABEL_KEY].astype(str).to_numpy()
    rows = np.flatnonzero(np.array([groups.get(x) == "CM" for x in lab]))
    before = ep.core.mean_cp10k(ep.core.counts_of(ref)[rows])
    after = ep.core.mean_cp10k(ep.core.counts_of(moved)[rows])
    m = ep.moved_mask(dirs["CM"], list(ref.var_names), held)
    slope = np.polyfit(dirs["CM"][m], np.log((after + 0.5) / (before + 0.5))[m], 1)[0]
    assert rb["CM"]["m1_r_gmom"] == pytest.approx(slope, rel=1e-9)
    assert ep.READBACK_STREAM not in (ep.SPERM_STREAM, ep.PLACEBO_STREAM, 0, 1, 2, 4, 5, 6)


def test_the_ruler_arms_are_the_real_directional_arms_at_their_build_fractions(ep):
    dirs, eps = {"CM": np.ones(3)}, {"e01": {"CM": np.zeros(3)}, "e002": {"CM": np.full(3, 2.0)}}
    arms = ep.ruler_arms("pooled", dirs, eps, 0.2)
    assert set(arms) == set(ep.candidates("pooled"))
    assert set(ep.ruler_arms("whole", dirs, eps, 0.2)) == set(ep.candidates("whole"))
    fr = ep.sr.arm_fractions(0.2, ep.uniform_of("pooled"))
    assert all(arms[a][1] == f and arms[a][0] is dirs for a, f in fr.items())
    assert arms["m1_r_gmom_e01"] == (eps["e01"], 0.2) and arms["m1_r_gmom_e002_x2"][0] is eps["e002"]
    assert arms["m1_r_gmom_e01_x15"][1] == arms["m1_r_gmom_x15"][1]


def _off_ruler(rec):
    for row in rec["rows"]:
        row["ruler_readback"]["CM"]["m1_r_gmom_e01"] = 0.25  # x15 at 0.175, half its smaller gap 0.029


def test_a_pair_not_ruler_matched_on_a_cut_is_flagged_and_the_statement_says_so(ep, tmp_path):
    r = _read(ep, tmp_path, _both({"m1_r_gmom_e01": 1.0, "m1_r_gmom_e01_x2": 1.0}), edit={"pooled": _off_ruler})
    assert r["valid"]
    m = r["statement_2"]["comparisons"][0]["ruler_match"]
    assert m["pooled"]["flagged_groups"] == ["CM"] and not m["pooled"]["matched"] and m["whole"]["matched"]
    assert m["pooled"]["groups"]["CM"]["half_gap"] == pytest.approx(0.029)
    assert m["pooled"]["neighbours"] == ["m1_r_gmom", "m1_r_gmom_x2"]
    assert r["statement_2"]["comparisons"][1]["ruler_match"]["pooled"]["matched"]
    want = (SMALLER + " (not ruler-matched on the cuts' reference cells: m1_r_gmom_e01 against m1_r_gmom_x15 on pooled "
            "(CM))")
    assert r["outcome"] == want == r["statement_2"]["statement"]
    txt = (tmp_path / "out" / "t52.txt").read_text()
    assert HEAD2 + want in txt and "NOT ruler-matched" in txt and "FLAGGED" in txt
    # matched pairs: the statement carries no suffix and the read-back is printed beside the pair
    r = _read(ep, tmp_path / "b", _both({"m1_r_gmom_e01": 1.0, "m1_r_gmom_e01_x2": 1.0}))
    assert r["outcome"] == SMALLER and not r["statement_2"]["not_ruler_matched"]
    txt = (tmp_path / "b" / "out" / "t52.txt").read_text()
    assert "read-back on pooled's reference cells" in txt and "read-back on whole's reference cells" in txt
    assert "(scored weak 0.030, strong 0.070)" in txt
    # on the whole cut _x3 has one neighbour on the ladder (no _x4 there)
    m2 = r["statement_2"]["comparisons"][1]["ruler_match"]
    assert m2["whole"]["neighbours"] == ["m1_r_gmom_x2"]
    assert m2["pooled"]["neighbours"] == ["m1_r_gmom_x2", "m1_r_gmom_x4"]


def _cmd_run_world(ep, monkeypatch, stats):
    """cmd_run on the synthetic world: the released files, the boards, the scorer and the anchors replaced; the
    directions, the arms and the records made by the script's own code from synthetic Feng statistics."""
    path, genes = stats
    assert ep.feng_dirs.clean_stats_path() == path
    ctx, _, _, held = _arms_world()
    task = types.SimpleNamespace(key=ep.TASK_OF["whole"], panel_board="pb", reference="T1/E9.5_RNA.h5ad",
                                 target="T1/E10.5_RNA.h5ad", target_filter=None,
                                 inputs={"ref": "T1/E9.5_RNA.h5ad", "start": "T1/E8.5_RNA.h5ad"})
    monkeypatch.setattr(ep.n12, "tries_rows", lambda: {"T52"})
    monkeypatch.setattr(ep.gc, "n22_passed", lambda p: True)
    monkeypatch.setattr(ep.gv, "tree_state", lambda: {"commit": "f" * 40, "clean": True})
    monkeypatch.setattr(ep.harness, "load_tasks", lambda: {"T1:test": [task]})
    monkeypatch.setattr(ep, "load_boards", lambda: {"pb": types.SimpleNamespace(genes=list(ctx.board.genes))})
    monkeypatch.setattr(ep, "truth_and_contamination", lambda *a: (types.SimpleNamespace(n=5), None))
    monkeypatch.setattr(ep, "read_stage", lambda *a: None)
    monkeypatch.setattr(ep, "Stage", types.SimpleNamespace(from_anndata=lambda *a: types.SimpleNamespace(n=3)))
    monkeypatch.setattr(ep.harness, "task_board", lambda *a: "board")
    monkeypatch.setattr(ep.lsr, "all_anchors", lambda t, b, tr, rs, seeds: (
        {s: types.SimpleNamespace(score=50.0, anchors=None) for s in seeds}, []))
    by_name = {"E9.5_RNA.h5ad": ctx.inputs["ref"], "E8.5_RNA.h5ad": ctx.inputs["start"]}
    monkeypatch.setattr(ep, "load_stage", lambda p, lineage=False: by_name[p.name])
    monkeypatch.setattr(ep.t1_scope, "held_flat_genes", lambda: tuple(held))
    monkeypatch.setattr(ep.hygiene, "load_union", lambda: {"file": "x"})
    monkeypatch.setattr(ep.gc, "held_feng", lambda union, pre: list(held))
    monkeypatch.setattr(ep.harness.TaskView, "of", staticmethod(lambda t: ctx.task))
    monkeypatch.setattr(ep.harness, "_candidate_board", lambda b: ctx.board)
    monkeypatch.setattr(ep.er, "score_arms", lambda board, truth, ref, s, anchors, preds: {
        a: {"score": 50.0 + i * 1e-3, "skills": {}} for i, a in enumerate(preds)})
    captured = []
    real_build = ep.build_seed

    def spy(ctx_, b0, truth, directions, eps_dirs, *rest):
        captured.append(eps_dirs)
        return real_build(ctx_, b0, truth, directions, eps_dirs, *rest)

    monkeypatch.setattr(ep, "build_seed", spy)
    return path, genes, captured


def test_the_run_makes_each_direction_set_with_its_own_pseudocount_end_to_end(ep, monkeypatch, stats, tmp_path):
    """Mutant M14 (both new sets made with 0.02 inside cmd_run): the run's direction step on synthetic statistics
    gives build_seed lfc_ratio at 0.1 and 0.02, the file's external records say so, and a run with the mutant voids the
    reading through those records."""
    path, genes, captured = _cmd_run_world(ep, monkeypatch, stats)
    ep.main(["run", "--truth", "whole", "--out", "t52", "--seeds", "0", "--data", str(tmp_path)])
    rec = json.loads((ep.RESULTS / "t52" / "tasks" / "whole.json").read_text())
    assert {t: rec["external"][t]["pseudocount_cp10k"] for t in ("entry", "e01", "e002")} == {
        "entry": 0.5, "e01": 0.1, "e002": 0.02}
    (eps_dirs,) = captured
    for tag, eps in ep.EPS.items():
        want, _ = ep.eps_directions(genes, eps, stats_path=path)
        assert ep.directions_sha256(eps_dirs[tag]) == ep.directions_sha256(want)
        assert rec["directions"][tag]["sha256"] == ep.directions_sha256(want)
    row = rec["rows"][0]
    assert row["ruler_readback"] and all(set(ep.candidates("whole")) <= set(v)
                                         for v in row["ruler_readback"].values())
    v = ep.validity({"whole": rec}, {})
    assert not any("external" in x or "read-back" in x for x in v["reasons"]), v["reasons"]
    # the mutant: every new set made with 0.02 while the script's constants still say 0.1 for e01
    real = ep.eps_directions
    monkeypatch.setattr(ep, "eps_directions", lambda g, eps, stats_path=None: real(g, 0.02, stats_path))
    ep.main(["run", "--truth", "whole", "--out", "t52m", "--seeds", "0", "--data", str(tmp_path)])
    bad = json.loads((ep.RESULTS / "t52m" / "tasks" / "whole.json").read_text())
    assert bad["pseudocounts"]["m1_r_gmom_e01"] == 0.1  # the per-arm record cannot see it
    v = ep.validity({"whole": bad}, {})
    assert "whole: wrong pseudocount in the direction code's own record (external) for ['e01']" in v["reasons"]
