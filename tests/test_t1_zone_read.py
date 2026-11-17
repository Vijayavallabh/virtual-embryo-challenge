"""experiments/t1_m1/zone_read.py (try T53 of experiments/harness/TRIES.md) on synthetic data: the zone-fixed
directions, the new arms and the reading. No released file, statistics file, board file or score is read."""

from __future__ import annotations

import importlib.util
import json
import sys

import numpy as np
import pytest
from test_t1_step_read import _arms_world

from vec.boards import GROUPS, REPO_ROOT

SCRIPT = REPO_ROOT / "experiments" / "t1_m1" / "zone_read.py"
STAMP = {"commit": "f" * 40, "tree_clean": True, "utc": "x", "inputs": {}}
NEW = ("m1_r_gmom_zlv", "m1_r_gmom_zlv_x2", "m1_r_gmom_zbal", "m1_r_gmom_zbal_x2")
POOLED_OF = {"m1_r_gmom_zlv": "m1_r_gmom", "m1_r_gmom_zlv_x2": "m1_r_gmom_x2", "m1_r_gmom_zbal": "m1_r_gmom",
             "m1_r_gmom_zbal_x2": "m1_r_gmom_x2"}
TRUTHS = ("pooled", "lib1", "lib2", "lib3", "whole")


@pytest.fixture
def zr(monkeypatch, tmp_path):
    spec = importlib.util.spec_from_file_location("zone_read_module", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    monkeypatch.setattr(mod, "RESULTS", tmp_path / "results")
    monkeypatch.setattr(mod.gv, "stamp", lambda started, files: dict(STAMP))
    return mod


# ---------------------------------------------------------------------------------------------------------------------
# Synthetic Feng statistics: keys stage|zone|type|half


FENG_GENES = np.array([f"g{i}" for i in range(6)], dtype=object)
CFG = {"metrics": {"pseudocount_cp10k": 0.5, "min_cells_ext": 10},
       "external_labels": {"CM": [["ventricular_cm", "atrial_cm"]], "ventricularCM": [["ventricular_cm"]],
                           "Endocardium": [["endocardial_ec"]], "Blood": [["blood"]]},
       "cm_includes": ["ventricularCM"]}
T1_GENES = ["g0", "g1", "g2", "g3", "g4", "g5", "gX"]  # gX is not in Feng: NaN on the map


def _stats(cells):
    """cells: list of (stage, zone, type, half, n, cp profile): a statistics table whose CP10k sums are n x profile."""
    keys, n, S_cp = [], [], []
    for stage, zone, typ, half, k, prof in cells:
        keys.append(f"{stage}|{zone}|{typ}|{half}")
        n.append(float(k))
        S_cp.append(k * np.asarray(prof, dtype=np.float64))
    S_cp = np.vstack(S_cp)
    return {"genes": FENG_GENES, "keys": np.array(keys, dtype=object), "n": np.array(n), "S_cp": S_cp,
            "S_log": np.log1p(S_cp), "S_nz": np.ones_like(S_cp), "depth_sum": np.array(n)}


def _rng_world(seed=0):
    """Four zones, three groups, both halves, random profiles; CM has every zone, Endocardium two, Blood too few cells
    in RA at E14.5."""
    rng = np.random.default_rng(seed)
    cells = []
    for stage in ("E9.5", "E14.5"):
        for zone in ("LA", "LV", "RA", "RV"):
            for typ in ("ventricular_cm", "atrial_cm", "endocardial_ec", "blood"):
                if typ == "endocardial_ec" and zone in ("LA", "RA"):
                    continue
                for half in ("A", "B"):
                    k = 3 if (typ == "blood" and zone == "RA" and stage == "E14.5") else int(rng.integers(8, 40))
                    cells.append((stage, zone, typ, half, k, rng.gamma(2.0, 3.0, len(FENG_GENES))))
    return _stats(cells)


def _entry(zr, st):
    return zr.feng_dirs.directions_from_stats(st, T1_GENES, CFG)[0]


def test_feng_sel_pools_every_zone_and_both_halves(zr):
    st = _rng_world()
    pooled = zr.dl.feng_sel(st, "E9.5", ["blood"])
    parts = [zr.zone_sel(st, "E9.5", z, ["blood"]) for z in ("LA", "LV", "RA", "RV")]
    assert pooled["n"] == sum(p["n"] for p in parts)
    assert np.allclose(pooled["cp"] * pooled["n"], sum(p["cp"] * p["n"] for p in parts))


def test_the_within_zone_direction_is_lfc_ratio_on_that_zones_pooled_halves(zr):
    st = _rng_world()
    within = zr.within_zone_directions(st, T1_GENES, CFG, ["CM", "Endocardium", "Blood"])
    idx = zr.dl.sym_index(st["genes"], T1_GENES)
    labs = {"ventricular_cm", "atrial_cm"}
    for z in ("LA", "LV", "RA", "RV"):
        a = zr.dl.agg(st, lambda p, z=z: p[0] == "E9.5" and p[1] == z and p[2] in labs)
        b = zr.dl.agg(st, lambda p, z=z: p[0] == "E14.5" and p[1] == z and p[2] in labs)
        want = zr.dl.on_t1(zr.dl.lfc_ratio(a, b, 0.5), idx)
        assert np.array_equal(within["CM"][z], want, equal_nan=True)
        assert np.isnan(within["CM"][z][-1])  # the unmapped gene
    assert sorted(within["Endocardium"]) == ["LV", "RV"]


def test_a_zone_under_the_minimum_at_either_end_gives_no_direction(zr):
    st = _rng_world()
    within = zr.within_zone_directions(st, T1_GENES, CFG, ["Blood"])
    assert "RA" not in within["Blood"] and sorted(within["Blood"]) == ["LA", "LV", "RV"]
    # the same under the minimum at the start end only
    cells = [("E9.5", "LV", "blood", "A", 4, np.ones(6)), ("E9.5", "LV", "blood", "B", 4, np.ones(6)),
             ("E14.5", "LV", "blood", "A", 40, np.ones(6)), ("E9.5", "RV", "blood", "A", 20, np.ones(6)),
             ("E14.5", "RV", "blood", "B", 20, np.ones(6))]
    within = zr.within_zone_directions(_stats(cells), T1_GENES, CFG, ["Blood"])
    assert sorted(within["Blood"]) == ["RV"]


def test_a_key_without_four_parts_is_refused(zr):
    st = _rng_world()
    st["keys"] = np.array([*st["keys"][:-1], "E9.5|LV|blood"], dtype=object)
    with pytest.raises(ValueError, match="four parts"):
        zr.within_zone_directions(st, T1_GENES, CFG, ["Blood"])
    st["keys"] = np.array([*st["keys"][:-1], "E9.5|LV|blood|A|x"], dtype=object)
    with pytest.raises(ValueError, match="four parts"):
        zr.zones_of(st)


def test_zlv_changes_the_cardiomyocyte_group_only_and_is_its_lv_direction(zr):
    st = _rng_world()
    entry = _entry(zr, st)
    within = zr.within_zone_directions(st, T1_GENES, CFG, entry)
    dirs, used = zr.zlv_directions(entry, within)
    assert set(dirs) == set(entry)
    assert np.array_equal(dirs["CM"], within["CM"]["LV"], equal_nan=True)
    assert not np.array_equal(dirs["CM"], entry["CM"], equal_nan=True)
    for g in entry:
        if g != "CM":
            assert dirs[g] is entry[g]
    assert used == {"CM": ["LV"], "Endocardium": [], "Blood": []}


def test_zlv_is_void_without_an_lv_direction(zr):
    entry = {"CM": np.zeros(3), "Blood": np.ones(3)}
    with pytest.raises(ValueError, match="no LV direction"):
        zr.zlv_directions(entry, {"CM": {"RV": np.ones(3)}, "Blood": {}})


def test_zbal_is_the_unweighted_mean_over_defined_zones_per_gene_over_finite_values(zr):
    nan = np.nan
    entry = {"CM": np.array([9.0, 9.0, 9.0, 9.0]), "Blood": np.array([1.0, 2.0, 3.0, 4.0])}
    within = {"CM": {"LA": np.array([1.0, nan, 3.0, nan]), "LV": np.array([3.0, 2.0, nan, nan]),
                     "RV": np.array([5.0, 4.0, nan, nan])},
              "Blood": {}}
    dirs, used = zr.zbal_directions(entry, within)
    assert np.array_equal(dirs["CM"], np.array([3.0, 3.0, 3.0, nan]), equal_nan=True)
    assert used["CM"] == ["LA", "LV", "RV"]
    assert dirs["Blood"] is entry["Blood"] and used["Blood"] == []
    rec = zr.direction_record(entry, {"zlv": ({"CM": within["CM"]["LV"], "Blood": entry["Blood"]},
                                              {"CM": ["LV"], "Blood": []}), "zbal": (dirs, used)},
                              ["a", "b", "c", "d"], [], ["LA", "LV", "RV"])
    assert rec["zbal_fallback"] == ["Blood"] and rec["zbal"]["Blood"]["pooled"] is True
    assert rec["zbal"]["CM"]["zones"] == ["LA", "LV", "RV"] and rec["zbal"]["CM"]["pooled"] is False


def test_one_zone_holding_every_cell_gives_both_forms_the_entrys_direction(zr):
    rng = np.random.default_rng(4)
    cells = [(stage, "LV", typ, half, int(rng.integers(10, 30)), rng.gamma(2.0, 3.0, 6))
             for stage in ("E9.5", "E14.5") for typ in ("ventricular_cm", "atrial_cm", "blood") for half in ("A", "B")]
    st = _stats(cells)
    entry = _entry(zr, st)
    forms = zr.form_directions(entry, zr.within_zone_directions(st, T1_GENES, CFG, entry))
    for form, (dirs, _) in forms.items():
        for g in entry:
            assert dirs[g].tobytes() == entry[g].tobytes(), (form, g)
            assert zr.direction_sha256(dirs[g]) == zr.direction_sha256(entry[g])


def test_a_world_differing_only_in_zone_mix_gives_a_pooled_direction_and_zero_zone_fixed_ones(zr):
    """Per zone the profile is the same at both ends; only the share of LV against atrial cells changes."""
    lv, la = np.array([10.0, 1.0, 5.0, 0.2, 3.0, 7.0]), np.array([1.0, 10.0, 5.0, 4.0, 0.5, 7.0])
    cells = [("E9.5", "LV", "ventricular_cm", "A", 70, lv), ("E9.5", "LA", "atrial_cm", "A", 30, la),
             ("E14.5", "LV", "ventricular_cm", "A", 50, lv), ("E14.5", "LA", "atrial_cm", "A", 50, la)]
    st = _stats(cells)
    entry = _entry(zr, st)
    assert np.nanmax(np.abs(entry["CM"])) > 0.1
    forms = zr.form_directions(entry, zr.within_zone_directions(st, T1_GENES, CFG, entry))
    for form in ("zlv", "zbal"):
        d = forms[form][0]["CM"]
        assert np.all(d[np.isfinite(d)] == 0.0) and np.isnan(d[-1])


def test_the_bracket_cut_gives_the_entrys_directions_bit_for_bit(zr):
    st = _rng_world()
    other = _stats([("E12.5", "LV", "blood", "A", 30, np.ones(6))])
    full = {k: (np.concatenate([st[k], other[k]]) if k != "genes" else st[k]) for k in st}
    cut = zr.feng_dirs.bracket_stats(full, zr.BRACKET)
    a, b = _entry(zr, full), _entry(zr, cut)
    assert set(a) == set(b) and all(a[g].tobytes() == b[g].tobytes() for g in a)
    assert zr.zones_of(cut) == ["LA", "LV", "RA", "RV"]


def test_the_direction_record_holds_ratios_names_and_hashes(zr):
    entry = {"CM": np.array([1.0, 2.0, np.nan, 4.0, 1.0]), "Blood": np.array([1.0, 1.0, 1.0, 1.0, 1.0])}
    zlv = {"CM": np.array([2.0, 4.0, 5.0, 8.0, -9.0]), "Blood": entry["Blood"]}
    rec = zr.direction_record(entry, {"zlv": (zlv, {"CM": ["LV"], "Blood": []}),
                                      "zbal": (zlv, {"CM": ["LV"], "Blood": []})},
                              ["a", "b", "c", "d", "h"], ["h"], ["LV"])
    cm = rec["zlv"]["CM"]
    assert cm["slope"] == pytest.approx(2.0) and cm["pearson"] == pytest.approx(1.0)
    assert cm["size_ratio"] == pytest.approx(4.0)  # held gene h and the NaN gene c left out
    assert rec["zlv"]["Blood"]["slope"] == pytest.approx(1.0) and rec["zlv"]["Blood"]["size_ratio"] == 1.0
    assert rec["zlv"]["Blood"]["sha256"] == rec["entry"]["Blood"]["sha256"] != cm["sha256"]
    assert zr.differing_groups(rec, "zlv") == ["CM"]
    text = json.dumps(rec)
    assert "\"n\"" not in text and "cells" not in text


# ---------------------------------------------------------------------------------------------------------------------
# The arms


def test_a_new_arm_differs_from_its_pooled_arms_build_in_the_directions_only(zr, monkeypatch):
    calls = []

    def spy(ctx, directions, groups, fraction, held_m1, held_factors):
        calls.append({"ctx": ctx, "directions": directions, "groups": groups, "fraction": fraction,
                      "held_m1": held_m1, "held_factors": held_factors})
        return {"m1_r": ("m1_r", fraction), "m1_r_gmom": ("m1_r_gmom", id(directions), fraction)}

    monkeypatch.setattr(zr.gc, "m1_r_arms", spy)
    ctx, entry, groups, held_m1, held_f = object(), {"CM": np.zeros(3)}, {"V-CM": "CM"}, ["a", "b"], ["a"]
    forms = {"zlv": {"CM": np.ones(3)}, "zbal": {"CM": np.full(3, 2.0)}}
    preds = zr.build_seed(ctx, lambda c: "b0", entry, forms, groups, 0.2, held_m1, held_f)
    assert set(preds) == set(zr.ARMS)
    by_fraction = {c["fraction"]: c for c in calls if c["directions"] is entry}
    assert sorted(by_fraction) == pytest.approx([0.2, 0.3, 0.4])
    news = [c for c in calls if c["directions"] is not entry]
    assert len(news) == 4
    for c, arm in zip(news, NEW, strict=True):
        base = by_fraction[0.4 if arm.endswith("_x2") else 0.2]
        assert c["fraction"] == base["fraction"]
        assert c["directions"] is forms["zlv" if "zlv" in arm else "zbal"]
        for key in ("ctx", "groups", "held_m1", "held_factors"):
            assert c[key] is base[key], f"{arm} differs from its pooled arm's build in {key}"


def test_a_form_equal_to_the_entry_builds_the_pooled_arms_bit_for_bit(zr):
    """The idle case on a real build: a form whose directions are the entry's (copies) gives the pooled arm's
    prediction bit for bit at both fractions; a changed direction moves it."""
    ctx, groups, dirs, held = _arms_world()
    copy = {g: v.copy() for g, v in dirs.items()}
    moved = {g: (v * 0.5 if g == "CM" else v) for g, v in dirs.items()}
    preds = zr.build_seed(ctx, lambda c: None, dirs, {"zlv": copy, "zbal": moved}, groups, 0.2, held, held)
    assert preds["m1_r_gmom_zlv"].X.tobytes() == preds["m1_r_gmom"].X.tobytes()
    assert preds["m1_r_gmom_zlv_x2"].X.tobytes() == preds["m1_r_gmom_x2"].X.tobytes()
    assert preds["m1_r_gmom_zbal"].X.tobytes() != preds["m1_r_gmom"].X.tobytes()
    direct = zr.gc.m1_r_arms(ctx, dirs, groups, 0.2, held, held)["m1_r_gmom"]
    assert preds["m1_r_gmom"].X.tobytes() == direct.X.tobytes()


def test_the_fractions(zr):
    fr = zr.arm_fractions(0.2)
    assert fr["m1_r_gmom"] == fr["m1_r_gmom_zlv"] == fr["m1_r_gmom_zbal"] == 0.2
    assert fr["m1_r_gmom_x2"] == fr["m1_r_gmom_zlv_x2"] == fr["m1_r_gmom_zbal_x2"] == 0.2 * 2.0
    assert fr["m1_r_gmom_x15"] == 0.2 * 1.5


# ---------------------------------------------------------------------------------------------------------------------
# The reading on synthetic score tables


W = {m: gw * mw for _, gw, ms in GROUPS["T1"] for m, mw in ms}
B0, BASE = 50.0, 55.0
UNI = {"m1_r_gmom": 0.0, "m1_r_gmom_x15": 0.4, "m1_r_gmom_x2": 0.8}


def _arm(score):
    gain = score - BASE
    sk = {"de_score": 0.4 + gain / 4 / (100 * W["de_score"]),
          "de_direction": 0.4 + gain / 4 / (100 * W["de_direction"]),
          "mmd_u": 0.3 + gain * 0.3 / (100 * W["mmd_u"]), "variogram": 0.2 + gain * 0.2 / (100 * W["variogram"])}
    return {"score": score, "groups": {"distribution": score / 100},
            "metrics": {"mmd_u": 0.1, "variogram": 0.01, "de_score": 0.5 + gain / 100, "de_direction": 0.3},
            "skills": sk}


def _sha(entry_groups=("CM", "Endocardium"), zlv_diff=("CM",), zbal_diff=("CM", "Endocardium")):
    rec = {"zones": ["LA", "LV", "RA", "RV"], "entry": {g: {"sha256": f"e-{g}"} for g in entry_groups},
           "zbal_fallback": []}
    for form, diff in (("zlv", zlv_diff), ("zbal", zbal_diff)):
        rec[form] = {g: {"zones": ["LV"], "pooled": False, "slope": 0.75 if g in diff else 1.0, "pearson": 0.9,
                         "size_ratio": 0.8, "sha256": f"{form}-{g}" if g in diff else f"e-{g}"}
                     for g in entry_groups}
    return rec


def _uniform(seed):
    return {a: BASE + g + 0.01 * seed for a, g in UNI.items()}


def _task(truth, gains=None, per_seed=None, seeds=range(10), commit="c" * 40, clean=True, fractions=None,
          identical=(), directions=None, anchor_failures=None):
    """gains: new arm -> gain over its pooled arm on every seed; per_seed: new arm -> ten gains. A new arm with no gain
    gets +d on seed 0 and -d on seed 1, so that it never scores exactly as its pooled arm."""
    gains = gains or {}
    rows = []
    for i, s in enumerate(seeds):
        uni = _uniform(s)
        arms = {"b0": _arm(B0), **{a: _arm(x) for a, x in uni.items()}}
        for j, a in enumerate(NEW):
            if per_seed and a in per_seed:
                g = per_seed[a][i]
            elif a in gains:
                g = gains[a]
            else:
                d = 1e-9 * (j + 1)
                g = d if i == 0 else -d if i == 1 else 0.0
            arms[a] = _arm(uni[POOLED_OF[a]] + g)
        for a in identical:
            arms[a] = dict(arms[POOLED_OF[a]])
        rows.append({"seed": s, "floor": 50.0, "arms": arms, "sha256": {a: f"{a}-{s}" for a in arms}})
    fr = fractions or {"m1_r_gmom": 0.2, "m1_r_gmom_x15": 0.2 * 1.5, "m1_r_gmom_x2": 0.4,
                       **{a: (0.4 if a.endswith("_x2") else 0.2) for a in NEW}}
    task = "T1/E9.5->E10.5" if truth == "whole" else "T1/E9.5->E10.5:chamber"
    rec = {"label": "synthetic", "task": task, "truth": truth, "fraction": 0.2, "fractions": fr, "seeds": list(seeds),
           "rows": rows, "directions": directions if directions is not None else _sha(),
           "stamp": {"commit": commit, "tree_clean": clean, "utc": "x", "inputs": {}}}
    if anchor_failures is not None:
        rec["anchor_failures"] = anchor_failures
    return rec


def _read(zr, tmp_path, by_truth=None, over=None, committed_shift=None, drop_committed=(), **kw):
    by_truth, over, committed_shift = by_truth or {}, over or {}, committed_shift or {}
    out = tmp_path / "out"
    (out / "tasks").mkdir(parents=True)
    for k in TRUTHS:
        if k in over and over[k] is None:
            continue
        rec = _task(k, gains=by_truth.get(k), **{**kw, **over.get(k, {})})
        (out / "tasks" / f"{k}.json").write_text(json.dumps(rec))
    files = {}
    for k in TRUTHS:
        if k in drop_committed:
            files[k] = tmp_path / f"absent_{k}.json"
            continue
        rows = []
        for s in range(10):
            arms = {a: _arm(x + committed_shift.get((k, a), 0.0)) for a, x in {"b0": B0, **_uniform(s)}.items()}
            rows.append({"seed": s, "arms": arms})
        files[k] = tmp_path / f"committed_{k}.json"
        files[k].write_text(json.dumps({"rows": rows}))
    return zr.read(out, files)["reading"]


UP, DOWN = 1.0, -1.0


def _both(form, g):
    return {k: {f"m1_r_gmom_{form}": g, f"m1_r_gmom_{form}_x2": g} for k in TRUTHS}


@pytest.mark.parametrize("by_truth, zlv, zbal", [
    (_both("zlv", UP), "zone-fixed scores better", "no difference at the bar"),
    (_both("zlv", DOWN), "zone-fixed scores worse", "no difference at the bar"),
    (_both("zbal", UP), "no difference at the bar", "zone-fixed scores better"),
    ({k: {"m1_r_gmom_zlv": UP} for k in TRUTHS}, "no difference at the bar", "no difference at the bar"),
    ({"pooled": {"m1_r_gmom_zlv": UP, "m1_r_gmom_zlv_x2": UP}}, "no difference at the bar",
     "no difference at the bar"),
    ({**{k: {"m1_r_gmom_zlv": UP} for k in TRUTHS}, **{k: {"m1_r_gmom_zlv_x2": DOWN} for k in ("pooled", "whole")}},
     "no difference at the bar", "no difference at the bar"),
    ({}, "no difference at the bar", "no difference at the bar"),
])
def test_each_statement_is_reached(zr, tmp_path, by_truth, zlv, zbal):
    r = _read(zr, tmp_path, by_truth)
    assert r["valid"], r["validity"]["reasons"]
    assert r["outcome"] == {"zlv": zlv, "zbal": zbal}
    assert r["statements"]["zlv"]["statement"] == zlv
    txt = (tmp_path / "out" / "t53.txt").read_text()
    assert f"zlv: {zlv}" in txt and f"zbal: {zbal}" in txt


def test_a_gain_of_exactly_point_three_does_not_pass(zr, tmp_path):
    r = _read(zr, tmp_path, {k: {"m1_r_gmom_zlv": 0.3, "m1_r_gmom_zlv_x2": 0.3} for k in TRUTHS})
    c = r["statements"]["zlv"]["comparisons"]["m1_r_gmom_zlv"]
    assert c["test"] == "neither" and r["outcome"]["zlv"] == "no difference at the bar"
    assert r["statements"]["zlv"]["library_labels"]["lib1"]["m1_r_gmom_zlv"] == "level"


def test_the_library_labels_are_ahead_behind_and_level(zr, tmp_path):
    r = _read(zr, tmp_path, {"lib1": {"m1_r_gmom_zbal": UP}, "lib2": {"m1_r_gmom_zbal": DOWN}})
    labels = r["statements"]["zbal"]["library_labels"]
    assert labels["lib1"]["m1_r_gmom_zbal"] == "ahead"
    assert labels["lib2"]["m1_r_gmom_zbal"] == "behind"
    assert labels["lib3"]["m1_r_gmom_zbal"] == "level"
    assert r["outcome"]["zbal"] == "no difference at the bar"


def test_the_reported_numbers(zr, tmp_path):
    r = _read(zr, tmp_path, _both("zlv", UP))
    c = r["comparisons"]["whole"]["m1_r_gmom_zlv_x2"]
    assert c["base"] == "m1_r_gmom_x2" and c["mean"] == pytest.approx(1.0) and c["seeds_ahead"] == 10
    assert c["points"]["de"] == pytest.approx(0.5) and c["points"]["mmd_u_variogram"] == pytest.approx(0.5)
    assert set(c["skills"]) >= {"de_score", "mmd_u"}
    assert r["mean_scores"]["lib2"]["m1_r_gmom_x2"] == pytest.approx(BASE + 0.8 + 0.045)
    four = r["statements"]["zlv"]["comparisons"]
    assert set(four) == {"m1_r_gmom_zlv", "m1_r_gmom_zlv_x2"}
    assert all(set(c["cuts"]) == {"pooled", "whole"} for c in four.values())


def test_what_size_alone_would_give(zr, tmp_path):
    r = _read(zr, tmp_path)
    sa = r["size_alone"]["pooled"]
    # CM's zlv slope 0.75: at the base fraction outside 1 to 2 (not computed); at x2 the multiple 1.5, which the
    # uniform arms put 0.4 below x2
    assert sa["m1_r_gmom_zlv"]["CM"]["what_size_alone_would_give"] is None
    assert sa["m1_r_gmom_zlv_x2"]["CM"]["what_size_alone_would_give"] == pytest.approx(-0.4)
    assert sa["m1_r_gmom_zlv"]["Endocardium"]["what_size_alone_would_give"] == pytest.approx(0.0)
    means = {"m1_r_gmom": 0.0, "m1_r_gmom_x15": 0.4, "m1_r_gmom_x2": 0.6}
    assert zr.size_alone(means, 1.25, "m1_r_gmom") == pytest.approx(0.2)
    assert zr.size_alone(means, 1.75, "m1_r_gmom_x2") == pytest.approx(-0.1)
    assert zr.size_alone(means, 2.5, "m1_r_gmom") is None and zr.size_alone(means, None, "m1_r_gmom") is None


def test_the_limits_are_in_the_text_and_no_count_is_printed(zr, tmp_path):
    _read(zr, tmp_path, _both("zbal", UP))
    txt = (tmp_path / "out" / "t53.txt").read_text()
    for x in zr.LIMITS:
        assert x in txt
    assert "what size alone would give" in txt and "read by nothing" in txt


@pytest.mark.parametrize("over, kw, reason", [
    ({"lib2": None}, {}, "missing truths ['lib2']"),
    ({"whole": {"seeds": range(9)}}, {}, "whole: seeds"),
    ({"lib3": {"commit": "d" * 40}}, {}, "ran at 2 commits"),
    ({"lib1": {"clean": None}}, {}, "tree not clean on ['lib1']"),
    ({"pooled": {"fractions": {"m1_r_gmom": 0.2, "m1_r_gmom_x15": 0.3, "m1_r_gmom_x2": 0.4, "m1_r_gmom_zlv": 0.2,
                               "m1_r_gmom_zlv_x2": 0.2, "m1_r_gmom_zbal": 0.2, "m1_r_gmom_zbal_x2": 0.4}}}, {},
     "wrong fraction on ['m1_r_gmom_zlv_x2']"),
    ({"lib3": {"identical": ("m1_r_gmom_zbal_x2",)}}, {}, "score exactly as the pooled arm"),
    ({"whole": {"directions": _sha(zlv_diff=("CM", "Endocardium"))}}, {}, "not in CM alone"),
    ({"pooled": {"directions": _sha(zlv_diff=())}}, {}, "not in CM alone"),
    ({"lib1": {"directions": _sha(zbal_diff=())}}, {}, "zbal's directions are the entry's in every group"),
    ({"lib1": {"directions": {}}}, {}, "no direction record"),
    ({"lib2": {"anchor_failures": [{"seed": 3, "reason": "x"}]}}, {}, "anchors invert or go missing"),
])
def test_the_void_cases(zr, tmp_path, over, kw, reason):
    r = _read(zr, tmp_path, over=over, **kw)
    assert not r["valid"] and r["outcome"] == "VOID"
    assert any(reason in x for x in r["validity"]["reasons"]), r["validity"]["reasons"]
    assert "Outcome: VOID" in (tmp_path / "out" / "t53.txt").read_text()


@pytest.mark.parametrize("truth, arm", [("pooled", "b0"), ("lib2", "m1_r_gmom"), ("lib3", "m1_r_gmom_x15"),
                                        ("whole", "m1_r_gmom_x2")])
def test_a_rescored_arm_off_its_committed_score_by_1e_6_voids(zr, tmp_path, truth, arm):
    r = _read(zr, tmp_path, committed_shift={(truth, arm): 1e-6})
    assert not r["valid"] and any(x.startswith(f"{truth}: ") and "reproduce" in x for x in r["validity"]["reasons"])
    assert _read(zr, tmp_path / "ok", committed_shift={(truth, arm): 1e-10})["valid"]


def test_a_missing_committed_file_voids(zr, tmp_path):
    r = _read(zr, tmp_path, drop_committed=("whole",))
    assert not r["valid"] and any("whole: " in x and "reproduce" in x for x in r["validity"]["reasons"])


def test_the_reading_refuses_to_overwrite(zr, tmp_path):
    _read(zr, tmp_path)
    with pytest.raises(SystemExit, match="already holds"):
        zr.read(tmp_path / "out", {})


# ---------------------------------------------------------------------------------------------------------------------
# The run's refusals


def test_the_run_needs_its_row(zr, monkeypatch):
    monkeypatch.setattr(zr.n12, "tries_rows", lambda: {"T46", "T52"})
    monkeypatch.setattr(zr.gc, "n22_passed", lambda path: True)
    with pytest.raises(SystemExit, match="T53 has no row"):
        zr.main(["run", "--truth", "whole", "--out", "t53"])


def test_the_run_needs_n22(zr, monkeypatch):
    monkeypatch.setattr(zr.n12, "tries_rows", lambda: {"T53"})
    with pytest.raises(SystemExit, match="does not hold a passing N22"):
        zr.main(["run", "--truth", "pooled", "--out", "t53", "--n22", "n22x"])


@pytest.mark.parametrize("truth", TRUTHS)
def test_the_run_refuses_an_existing_task_file(zr, monkeypatch, truth):
    monkeypatch.setattr(zr.n12, "tries_rows", lambda: {"T53"})
    monkeypatch.setattr(zr.gc, "n22_passed", lambda path: True)
    tasks = zr.RESULTS / "t53" / "tasks"
    tasks.mkdir(parents=True)
    (tasks / f"{truth}.json").write_text("{}")
    with pytest.raises(SystemExit, match="already holds"):
        zr.main(["run", "--truth", truth, "--out", "t53"])


def test_the_run_refuses_an_unknown_truth(zr):
    with pytest.raises(SystemExit):
        zr.main(["run", "--truth", "lib4", "--out", "t53"])


def test_the_truths_go_to_their_cuts(zr):
    assert zr.TASK_OF == {"pooled": "T1/E9.5->E10.5:chamber", "lib1": "T1/E9.5->E10.5:chamber",
                          "lib2": "T1/E9.5->E10.5:chamber", "lib3": "T1/E9.5->E10.5:chamber",
                          "whole": "T1/E9.5->E10.5"}
    assert zr.COMMITTED["lib2"].parts[-3:] == ("t46_m1_lib_step", "tasks", "lib2.json")
    assert zr.COMMITTED["whole"].parts[-3:] == ("t44_m1_step", "tasks", "T1_E9.5_E10.5.json")
