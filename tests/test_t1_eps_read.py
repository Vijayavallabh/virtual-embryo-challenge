"""experiments/t1_m1/eps_read.py (try T52 of experiments/harness/TRIES.md) on synthetic data: the directions made with
another pseudocount, the arms, the placebos and the reading. No released file, real statistics file, board file or
score is read."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import re
import sys

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
        assert ep.arms_of(k) == ("b0", "m1_r_gmom", *UNIFORM, "m1_r_gmom_perm", *NEW, *NEW_PERM)
    assert ep.arms_of("whole") == ("b0", "m1_r_gmom", "m1_r_gmom_x15", "m1_r_gmom_x2", "m1_r_gmom_x3", *NEW)
    spec = ep.arm_spec("lib2")
    assert spec["m1_r_gmom"] == {"fraction": 0.2, "pseudocount": 0.5}
    assert [spec[a]["fraction"] for a in UNIFORM] == pytest.approx([0.3, 0.4, 0.6, 0.8])
    assert [spec[a]["fraction"] for a in NEW] == pytest.approx([0.2, 0.3, 0.4] * 2)
    assert [spec[a]["pseudocount"] for a in NEW] == [0.1] * 3 + [0.02] * 3
    assert spec["m1_r_gmom_e01_x15_perm"] == {"fraction": pytest.approx(0.3), "pseudocount": 0.1}
    assert spec["m1_r_gmom_e002_x15_perm"] == {"fraction": pytest.approx(0.3), "pseudocount": 0.02}
    assert spec["m1_r_gmom_perm"] == {"fraction": 0.2, "pseudocount": 0.5}
    assert "m1_r_gmom_x4" not in ep.arm_spec("whole") and "m1_r_gmom_perm" not in ep.arm_spec("whole")


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
    preds, _, _ = ep.build_seed(ctx, lambda c: "b0", "pooled", dirs, eps_dirs, groups, 0.2, held, held)
    assert list(preds) == list(ep.arms_of("pooled"))
    by_arm = dict(zip([a for a in ep.arms_of("pooled") if a != "b0"], calls, strict=True))
    base = by_arm["m1_r_gmom"]
    for a in ("m1_r_gmom", *UNIFORM, *NEW, "m1_r_gmom_perm", *NEW_PERM):
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


def test_the_whole_cut_builds_no_x4_and_no_placebo(ep, monkeypatch):
    calls = []
    monkeypatch.setattr(ep.gc, "m1_r_arms", lambda *a: calls.append(a[3]) or {"m1_r": None, "m1_r_gmom": a[3]})
    ctx, groups, dirs, eps_dirs, held = _eps_world(ep)
    preds, _, placebo = ep.build_seed(ctx, lambda c: "b0", "whole", dirs, eps_dirs, groups, 0.2, held, held)
    assert list(preds) == list(ep.arms_of("whole")) and len(calls) == 4 + 6
    assert placebo == {"sha256": {}, "same_positions_as_m1_r_gmom_perm": {}}


def test_the_rescored_arms_keep_their_output_bit_for_bit(ep):
    """b0, m1_r_gmom, the four uniform arms and m1_r_gmom_perm as build_seed makes them now equal try T46's committed
    builds on the synthetic world (test_t1_lib_step_read's golden sha256)."""
    ctx, groups, dirs, eps_dirs, held = _eps_world(ep)
    preds, log_factors, placebo = ep.build_seed(ctx, harness._b0_t1(harness.B0_DAMPING), "pooled", dirs, eps_dirs,
                                                groups, 0.2, held, held)
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
    assert set(placebo["sha256"]) == {"m1_r_gmom_perm", *NEW_PERM}


def test_with_the_entrys_directions_the_new_arms_are_the_entrys_arms_bit_for_bit(ep):
    """The idle case: new direction sets equal to the entry's make each new arm the uniform arm of its fraction (and
    the new placebo the _x15 build along the entry's own permutation) bit for bit."""
    ctx, groups, dirs, held = _arms_world()
    idle = {tag: {g: v.copy() for g, v in dirs.items()} for tag in ep.EPS}
    preds, _, placebo = ep.build_seed(ctx, harness._b0_t1(harness.B0_DAMPING), "pooled", dirs, idle, groups, 0.2,
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
# The reading


ROWS = {"pooled": 7919, "lib1": 1777, "lib2": 2711, "lib3": 3431, "whole": 9127}
CONTAM = {"pooled": 613, "lib1": 401, "lib2": 199, "lib3": 13, "whole": 0}
DIRS = {"entry": {"pseudocount": 0.5, "sha256": {"CM": "e-cm", "Endocardium": "e-en"}},
        "e01": {"pseudocount": 0.1, "sha256": {"CM": "a-cm", "Endocardium": "a-en"}},
        "e002": {"pseudocount": 0.02, "sha256": {"CM": "b-cm", "Endocardium": "b-en"}}}


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
        rows.append({"seed": s, "floor": 50.0, "arms": arms,
                     "sha256": {a: (sha or {}).get(a, f"{a}-{s}") for a in arms_list},
                     "sparsity": {a: {"zero_share": 0.25, "detected_per_cell": 6.0} for a in arms_list},
                     "scored_log_factors": lf,
                     "placebo": {"same_positions_as_m1_r_gmom_perm": {p: True for p in NEW_PERM}}})
    spec = ep.arm_spec(truth)
    return {"label": "synthetic", "task": ep.TASK_OF[truth], "truth": truth, "truth_rows": ROWS[truth],
            "truth_contamination": CONTAM[truth], "fraction": 0.2,
            "fractions": fractions or {a: v["fraction"] for a, v in spec.items()},
            "pseudocounts": pseudocounts or {a: v["pseudocount"] for a, v in spec.items()},
            "directions": json.loads(json.dumps(directions or DIRS)),
            "direction_stats": {t: {"CM": {"squared_size_over_entry": 1.0, "share_on_low_start_genes": 0.4}}
                                for t in DIRS},
            "seeds": list(seeds), "rows": rows,
            "stamp": {"commit": commit, "tree_clean": clean, "utc": "x", "inputs": {}}}


def _write(ep, tmp_path, truths, committed_shift=None):
    """truths: truth -> task record; the committed files (T46's per chamber truth, T44's whole cut) hold the rescored
    arms' scores from the same rows unless shifted (truth -> arm -> shift)."""
    committed_shift = committed_shift or {}
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
        (t44 if k == "whole" else t46 / f"{k}.json").write_text(json.dumps({"rows": rows}))
    return out, t46, t44


def _read(ep, tmp_path, by_truth=None, over=None, committed_shift=None, edit=None):
    """by_truth: truth -> gains (both cuts: give 'pooled' and 'whole'); over: truth -> _task overrides (None drops
    the truth); edit: truth -> a function changing the record after it is made (the committed files are written
    first, from the unedited rows)."""
    by_truth, over = by_truth or {}, over or {}
    truths = {k: _task(ep, k, gains=by_truth.get(k), **(over.get(k) or {})) for k in TRUTHS
              if not (k in over and over[k] is None)}
    out, t46, t44 = _write(ep, tmp_path, {k: v for k, v in truths.items()}, committed_shift)
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
    ({"m1_r_gmom_e01": 1.0, "m1_r_gmom_e01_x2": 0.01}, "no difference at the bar"),
    ({"m1_r_gmom_e01": 1.0, "m1_r_gmom_e01_x2": -1.0}, "no difference at the bar"),
    ({"m1_r_gmom_x15": 1.0, "m1_r_gmom_e01": 2.0, "m1_r_gmom_x3": 2.0, "m1_r_gmom_e01_x2": 1.0},
     "no difference at the bar"),
    ({"m1_r_gmom_x15": 1.0, "m1_r_gmom_e01": 2.0, "m1_r_gmom_x3": 2.0, "m1_r_gmom_e01_x2": 3.0}, SMALLER),
])
def test_statement_two(ep, tmp_path, gains, statement):
    r = _read(ep, tmp_path, _both(gains))
    assert r["valid"] and r["outcome"] == statement == r["statement_2"]["statement"]
    c1, c2 = r["statement_2"]["comparisons"]
    assert (c1["two_cuts"]["change"], c1["two_cuts"]["base"]) == ("m1_r_gmom_e01", "m1_r_gmom_x15")
    assert (c2["two_cuts"]["change"], c2["two_cuts"]["base"]) == ("m1_r_gmom_e01_x2", "m1_r_gmom_x3")
    txt = (tmp_path / "out" / "t52.txt").read_text()
    assert HEAD2 + statement in txt
    assert "m1_r_gmom_e01 against m1_r_gmom_x15 on the two cuts: " in txt
    assert "m1_r_gmom_e01_x2 against m1_r_gmom_x3 on the two cuts: " in txt
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
    ({"m1_r_gmom_e01_x15": 1.0, "m1_r_gmom_e01_x15_perm": 0.0}, "the gain needs the direction's genes"),
    ({"m1_r_gmom_e01_x15": 0.01, "m1_r_gmom_e01_x15_perm": -2.0}, "no gain on this truth"),
    ({"m1_r_gmom_e01_x15": 1.0, "m1_r_gmom_e01_x15_perm": 0.9}, "not separated"),
    ({"m1_r_gmom_e01_x15": 1.0, "m1_r_gmom_e01_x15_perm": 1.5}, "not separated"),
])
def test_statement_three(ep, tmp_path, gains, label):
    r = _read(ep, tmp_path, {"lib3": gains})
    s = r["statement_3"]["lib3"]["e01"]
    assert s["label"] == label
    d = s["did"]
    assert d["mean"] == pytest.approx(gains["m1_r_gmom_e01_x15"] - gains["m1_r_gmom_e01_x15_perm"], abs=1e-7)
    assert r["statement_3"]["lib1"]["e01"]["label"] == "no gain on this truth"
    assert r["statement_3"]["lib3"]["e002"]["label"] == "no gain on this truth"
    assert "whole" not in r["statement_3"]


def test_d_uses_the_entrys_placebo_as_the_placebo_base(ep, tmp_path):
    r = _read(ep, tmp_path, {"lib2": {"m1_r_gmom_e002_x15": 1.0, "m1_r_gmom_e002_x15_perm": 0.5,
                                      "m1_r_gmom_perm": -0.5}})
    s = r["statement_3"]["lib2"]["e002"]
    assert s["did"]["mean"] == pytest.approx(1.0 - (0.5 + 0.5)) and s["label"] == "not separated"


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
