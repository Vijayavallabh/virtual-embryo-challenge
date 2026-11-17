"""experiments/t1_m1/zone_read.py (try T53 of experiments/harness/TRIES.md, as amended after its first review) on
synthetic data: the zone-fixed directions (zlv, zw), their slope-matched versions and controls, the new arms and the
reading. No released file, statistics file, board file or score is read."""

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
RAW = ("m1_r_gmom_zlv", "m1_r_gmom_zlv_x2", "m1_r_gmom_zw", "m1_r_gmom_zw_x2")
SMA = ("m1_r_gmom_zlv_sm", "m1_r_gmom_zlv_sm_x2", "m1_r_gmom_zw_sm", "m1_r_gmom_zw_sm_x2")
CTRL = ("m1_r_gmom_zlv_sm_perm", "m1_r_gmom_zw_sm_perm")
NEW = (*RAW, *SMA, *CTRL)
POOLED_OF = {a: ("m1_r_gmom_x2" if a.endswith("_x2") else "m1_r_gmom") for a in NEW}
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
    with pytest.raises(ValueError, match="four parts"):
        zr.zw_directions(st, {"Blood": np.zeros(7)}, T1_GENES, CFG)
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
    assert zr.changed_groups(entry, dirs) == ["CM"]


def test_zlv_is_void_without_an_lv_direction(zr):
    entry = {"CM": np.zeros(3), "Blood": np.ones(3)}
    with pytest.raises(ValueError, match="no LV direction"):
        zr.zlv_directions(entry, {"CM": {"RV": np.ones(3)}, "Blood": {}})


# ---------------------------------------------------------------------------------------------------------------------
# zw


def _independent_zw(cells, labels, eps=0.5, mc=10):
    """zw by hand from the cell list: E9.5 zone shares over the kept zones, both ends' zone means at those same shares,
    on the six Feng genes."""
    labels = set(labels)
    n9, s9, n14, s14 = {}, {}, {}, {}
    for stage, zone, typ, _half, k, prof in cells:
        if typ not in labels:
            continue
        n, s = (n9, s9) if stage == "E9.5" else (n14, s14)
        n[zone] = n.get(zone, 0) + k
        s[zone] = s.get(zone, 0) + k * np.asarray(prof, float)
    total9 = sum(n9.values())
    keep = [z for z in n9 if n14.get(z, 0) >= mc]
    w = {z: n9[z] / total9 for z in keep}
    a = sum(w[z] * s9[z] / n9[z] for z in keep) / sum(w.values())
    b = sum(w[z] * s14[z] / n14[z] for z in keep) / sum(w.values())
    return np.log((b + eps) / (a + eps)), 1 - sum(w.values())


def test_zw_matches_an_independent_computation_on_unequal_zone_cell_numbers(zr):
    rng = np.random.default_rng(11)
    cells = []
    for stage, counts in (("E9.5", {"LA": 7, "LV": 61, "RA": 5, "RV": 23}),
                          ("E14.5", {"LA": 30, "LV": 17, "RA": 4, "RV": 52})):
        for zone, k in counts.items():
            typ = "atrial_cm" if zone in ("LA", "RA") else "ventricular_cm"
            for half, kk in (("A", k // 2 + 1), ("B", k - k // 2)):
                cells.append((stage, zone, typ, half, kk, rng.gamma(2.0, 3.0, 6)))
    st = _stats(cells)
    entry = _entry(zr, st)
    dirs, used, info = zr.zw_directions(st, entry, T1_GENES, CFG)
    want, lost = _independent_zw(cells, ["ventricular_cm", "atrial_cm"])
    assert np.allclose(dirs["CM"][:6], want, rtol=0, atol=1e-12) and np.isnan(dirs["CM"][6])
    assert used["CM"] == ["LA", "LV", "RV"]  # RA under the minimum at E14.5
    assert info["CM"]["dropped"] == ["RA"] and info["CM"]["dropped_share"] == pytest.approx(lost)
    assert info["CM"]["fallback"] is False
    assert not np.allclose(dirs["CM"][:6], entry["CM"][:6])


def test_zw_is_exactly_zero_in_a_world_that_differs_between_the_ends_only_in_zone_mix(zr):
    """Per zone the profile is the same at both ends; only the zone shares change (dyadic numbers, so the pooled E9.5
    mean and the reweighted E14.5 mean are exact)."""
    lv, la = np.array([10.0, 1.0, 5.0, 0.25, 3.0, 7.0]), np.array([1.0, 10.0, 5.0, 4.0, 0.5, 7.0])
    cells = [("E9.5", "LV", "ventricular_cm", "A", 48, lv), ("E9.5", "LA", "atrial_cm", "A", 16, la),
             ("E14.5", "LV", "ventricular_cm", "A", 32, lv), ("E14.5", "LA", "atrial_cm", "A", 32, la)]
    st = _stats(cells)
    entry = _entry(zr, st)
    assert np.nanmax(np.abs(entry["CM"])) > 0.1
    zw = zr.zw_directions(st, entry, T1_GENES, CFG)[0]["CM"]
    assert np.all(zw[:6] == 0.0) and np.isnan(zw[6])
    zlv = zr.zlv_directions(entry, zr.within_zone_directions(st, T1_GENES, CFG, entry))[0]["CM"]
    assert np.all(zlv[:6] == 0.0)


def test_zw_is_zero_to_rounding_in_a_random_zone_mix_world(zr):
    rng = np.random.default_rng(3)
    prof = {z: rng.gamma(2.0, 3.0, 6) for z in ("LA", "LV", "RA", "RV")}
    cells = [(s, z, "atrial_cm" if z[0] == "L" else "ventricular_cm", "A", int(rng.integers(10, 90)), prof[z])
             for s in ("E9.5", "E14.5") for z in prof]
    st = _stats(cells)
    entry = _entry(zr, st)
    zw = zr.zw_directions(st, entry, T1_GENES, CFG)[0]["CM"]
    assert np.nanmax(np.abs(entry["CM"])) > 1e-3
    assert np.nanmax(np.abs(zw)) < 1e-12


def test_zw_equals_the_entrys_direction_when_the_zone_mix_is_the_same_at_both_ends(zr):
    rng = np.random.default_rng(5)
    cells = []
    for stage, scale in (("E9.5", 1), ("E14.5", 3)):
        for zone, k in (("LA", 11), ("LV", 40), ("RV", 25)):
            for typ in ("ventricular_cm", "atrial_cm", "blood"):
                cells.append((stage, zone, typ, "A", k * scale, rng.gamma(2.0, 3.0, 6)))
    st = _stats(cells)
    entry = _entry(zr, st)
    dirs, used, info = zr.zw_directions(st, entry, T1_GENES, CFG)
    for g in entry:
        assert np.allclose(dirs[g], entry[g], rtol=0, atol=1e-12, equal_nan=True), g
        assert info[g]["dropped"] == [] and info[g]["dropped_share"] == 0.0
    assert np.nanmax(np.abs(entry["CM"])) > 0.1


def test_zw_keeps_the_entry_where_more_than_half_the_weight_is_dropped(zr):
    one = np.ones(6)
    # Blood: LV holds 3/4 of E9.5 and is under the minimum at E14.5 (dropped share 0.75, falls back);
    # CM: LV holds exactly half and is dropped (0.5 is not more than half: kept on RV alone)
    cells = [("E9.5", "LV", "blood", "A", 30, one), ("E9.5", "RV", "blood", "A", 10, one),
             ("E14.5", "LV", "blood", "A", 5, one * 2), ("E14.5", "RV", "blood", "A", 20, one * 3),
             ("E9.5", "LV", "ventricular_cm", "A", 20, one), ("E9.5", "RV", "ventricular_cm", "A", 20, one),
             ("E14.5", "LV", "ventricular_cm", "A", 4, one), ("E14.5", "RV", "ventricular_cm", "A", 20, one * 4)]
    st = _stats(cells)
    entry = _entry(zr, st)
    dirs, used, info = zr.zw_directions(st, entry, T1_GENES, CFG)
    assert dirs["Blood"] is entry["Blood"] and used["Blood"] == [] and info["Blood"]["fallback"] is True
    assert info["Blood"]["dropped"] == ["LV"] and info["Blood"]["dropped_share"] == 0.75
    assert info["CM"]["fallback"] is False and used["CM"] == ["RV"] and info["CM"]["dropped_share"] == 0.5
    assert np.allclose(dirs["CM"][:6], np.log(4.5 / 1.5))
    rec = zr.direction_record(entry, {"zw": (dirs, used)}, T1_GENES, [], ["LV", "RV"], info)
    assert rec["zw_fallback"] == ["Blood"] and rec["zw"]["Blood"]["dropped_share"] == 0.75
    assert "\"n\"" not in json.dumps(rec)


def test_one_zone_holding_every_cell_gives_every_form_the_entrys_direction(zr):
    rng = np.random.default_rng(4)
    cells = [(stage, "LV", typ, half, int(rng.integers(10, 30)), rng.gamma(2.0, 3.0, 6))
             for stage in ("E9.5", "E14.5") for typ in ("ventricular_cm", "atrial_cm", "blood") for half in ("A", "B")]
    st = _stats(cells)
    entry = _entry(zr, st)
    forms, rec = zr.all_directions(st, entry, T1_GENES, CFG, [])
    for form in ("zlv", "zw", "zlv_sm", "zw_sm"):
        for g in entry:
            assert forms[form][g].tobytes() == entry[g].tobytes(), (form, g)
        assert rec["changed"][form] == []


def test_the_bracket_cut_gives_the_entrys_directions_bit_for_bit(zr):
    st = _rng_world()
    other = _stats([("E12.5", "LV", "blood", "A", 30, np.ones(6))])
    full = {k: (np.concatenate([st[k], other[k]]) if k != "genes" else st[k]) for k in st}
    cut = zr.feng_dirs.bracket_stats(full, zr.BRACKET)
    a, b = _entry(zr, full), _entry(zr, cut)
    assert set(a) == set(b) and all(a[g].tobytes() == b[g].tobytes() for g in a)
    assert zr.zones_of(cut) == ["LA", "LV", "RA", "RV"]


# ---------------------------------------------------------------------------------------------------------------------
# Slope matching and the controls


def test_a_slope_matched_direction_has_slope_one_on_the_entrys(zr):
    # seed 1: in this world of independent random profiles per zone and stage, seed 0's zw Blood (RA dropped) has a
    # slope of 0.03 under the kept-zone E9.5 side and is refused (the refusal has its own tests)
    st = _rng_world(1)
    entry = _entry(zr, st)
    forms, rec = zr.all_directions(st, entry, T1_GENES, CFG, ["g5"])
    for form in ("zlv_sm", "zw_sm"):
        assert rec["changed"][form]
        for g in rec["changed"][form]:
            assert rec[form][g]["slope"] == pytest.approx(1.0, abs=1e-12)
            raw = rec[form.removesuffix("_sm")][g]["slope"]
            assert rec[form][g]["divisor"] == raw
            assert np.array_equal(forms[form][g], forms[form.removesuffix("_sm")][g] / raw, equal_nan=True)
        for g in set(entry) - set(rec["changed"][form]):
            assert forms[form][g] is entry[g]
    assert rec["changed"]["zlv_sm"] == ["CM"]


@pytest.mark.parametrize("scale", [0.2 - 1e-9, 0.1, -1.0, 0.0])
def test_the_run_refuses_a_changed_group_whose_slope_is_not_above_point_two(zr, scale):
    entry = {"CM": np.array([1.0, -2.0, 3.0, 0.5]), "Blood": np.ones(4)}
    form = {"CM": entry["CM"] * scale, "Blood": entry["Blood"]}
    with pytest.raises(ValueError, match="not above 0.2"):
        zr.slope_matched(entry, form, ["a", "b", "c", "d"], [])
    ok, div = zr.slope_matched(entry, {"CM": entry["CM"] * 0.21, "Blood": entry["Blood"]}, ["a", "b", "c", "d"], [])
    assert div == {"CM": pytest.approx(0.21)} and ok["Blood"] is entry["Blood"]


def test_all_directions_refuses_through_the_slope(zr, monkeypatch):
    st = _rng_world()
    entry = _entry(zr, st)
    monkeypatch.setattr(zr, "MIN_SLOPE", 1e9)
    with pytest.raises(ValueError, match="not above"):
        zr.all_directions(st, entry, T1_GENES, CFG, [])


def test_a_control_carries_the_same_departures_permuted_within_the_two_strata(zr):
    genes = [f"s{i}" for i in range(5)] + [f"o{i}" for i in range(6)] + ["h", "n"]
    scored = [f"s{i}" for i in range(5)]
    rng = np.random.default_rng(0)
    e = {"CM": np.append(rng.normal(size=12), np.nan), "Blood": np.append(rng.normal(size=12), np.nan),
         "Endo": np.append(rng.normal(size=12), np.nan)}
    sm = {"CM": e["CM"] + np.append(rng.normal(size=12), np.nan), "Blood": e["Blood"],
          "Endo": e["Endo"] * 1.0 + np.append(np.arange(12.0), np.nan)}
    sm["CM"][11] = 99.0  # the held gene h: not moved, so the control keeps the entry's value
    c = zr.control_directions(e, sm, genes, ["h"], scored, zr.control_rng([1000 + 3]))
    assert c["Blood"] is e["Blood"]
    for g in ("CM", "Endo"):
        r = sm[g] - e[g]
        d = c[g] - e[g]
        for sl in (slice(0, 5), slice(5, 11)):
            assert np.allclose(np.sort(d[sl]), np.sort(r[sl]), rtol=0, atol=1e-12)
        assert c[g][11] == e[g][11] and np.isnan(c[g][12])
        assert not np.allclose(d[:11], r[:11])
    again = zr.control_directions(e, sm, genes, ["h"], scored, zr.control_rng(1003))
    assert all(np.array_equal(c[g], again[g], equal_nan=True) for g in c)
    other = zr.control_directions(e, sm, genes, ["h"], scored, zr.control_rng([1004]))
    assert not np.array_equal(c["CM"], other["CM"], equal_nan=True)
    # the stream: [build seed, 53], groups in sorted order (CM before Endo), the scored stratum first
    g53 = np.random.default_rng([1003, 53])
    want = {}
    for g in ("CM", "Endo"):
        v = e[g].copy()
        r = sm[g] - e[g]
        for ix in (np.arange(5), np.arange(5, 11)):
            v[ix] = e[g][ix] + r[ix][g53.permutation(len(ix))]
        want[g] = v
    assert all(np.array_equal(c[g], want[g], equal_nan=True) for g in want)


# ---------------------------------------------------------------------------------------------------------------------
# The direction record


def test_the_direction_record_holds_ratios_names_and_hashes(zr):
    entry = {"CM": np.array([1.0, 2.0, np.nan, 4.0, 1.0]), "Blood": np.array([1.0, 1.0, 1.0, 1.0, 1.0])}
    zlv = {"CM": np.array([2.0, 4.0, 5.0, 8.0, -9.0]), "Blood": entry["Blood"]}
    rec = zr.direction_record(entry, {"zlv": (zlv, {"CM": ["LV"], "Blood": []})}, ["a", "b", "c", "d", "h"], ["h"],
                              ["LV"], divisors={"zlv": {"CM": 2.0}})
    cm = rec["zlv"]["CM"]
    assert cm["slope"] == pytest.approx(2.0) and cm["pearson"] == pytest.approx(1.0) and cm["divisor"] == 2.0
    assert cm["size_ratio"] == pytest.approx(4.0)  # held gene h and the NaN gene c left out
    assert rec["zlv"]["Blood"]["slope"] == pytest.approx(1.0) and rec["zlv"]["Blood"]["size_ratio"] == 1.0
    assert rec["zlv"]["Blood"]["sha256"] == rec["entry"]["Blood"]["sha256"] != cm["sha256"]
    assert zr.differing_groups(rec, "zlv") == ["CM"] == rec["changed"]["zlv"]
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

    ctx, entry, groups, held_m1, held_f = object(), {"CM": np.zeros(3)}, {"V-CM": "CM"}, ["a", "b"], ["a"]
    forms = {k: {"CM": np.full(3, float(i + 1))} for i, k in enumerate(("zlv", "zw", "zlv_sm", "zw_sm"))}
    controls = {"zlv_sm_perm": {"CM": np.full(3, 7.0)}, "zw_sm_perm": {"CM": np.full(3, 8.0)}}
    seen = []

    def fake_controls(c, e, f, h):
        seen.append((c, e, f, h))
        return controls

    monkeypatch.setattr(zr.gc, "m1_r_arms", spy)
    monkeypatch.setattr(zr, "seed_controls", fake_controls)
    preds = zr.build_seed(ctx, lambda c: "b0", entry, forms, groups, 0.2, held_m1, held_f)
    assert set(preds) == set(zr.ARMS) and len(zr.NEW_ARMS) == 10
    assert seen == [(ctx, entry, forms, held_m1)]
    by_fraction = {c["fraction"]: c for c in calls if c["directions"] is entry}
    assert sorted(by_fraction) == pytest.approx([0.2, 0.3, 0.4])
    news = [c for c in calls if c["directions"] is not entry]
    assert len(news) == 10
    want = {**forms, **controls}
    for c, (arm, (key, pooled)) in zip(news, zr.NEW_ARMS.items(), strict=True):
        base = by_fraction[0.4 if pooled == "m1_r_gmom_x2" else 0.2]
        assert pooled == POOLED_OF[arm]
        assert c["fraction"] == base["fraction"]
        assert c["directions"] is want[key]
        for k in ("ctx", "groups", "held_m1", "held_factors"):
            assert c[k] is base[k], f"{arm} differs from its pooled arm's build in {k}"
    assert {a: k for a, (k, _) in zr.NEW_ARMS.items()} == {
        "m1_r_gmom_zlv": "zlv", "m1_r_gmom_zlv_x2": "zlv", "m1_r_gmom_zw": "zw", "m1_r_gmom_zw_x2": "zw",
        "m1_r_gmom_zlv_sm": "zlv_sm", "m1_r_gmom_zlv_sm_x2": "zlv_sm", "m1_r_gmom_zw_sm": "zw_sm",
        "m1_r_gmom_zw_sm_x2": "zw_sm", "m1_r_gmom_zlv_sm_perm": "zlv_sm_perm", "m1_r_gmom_zw_sm_perm": "zw_sm_perm"}


def test_a_form_equal_to_the_entry_builds_the_pooled_arms_bit_for_bit(zr):
    """The idle case on a real build: forms whose directions are the entry's (copies) give the pooled arms' predictions
    bit for bit, controls included; a changed slope-matched direction moves its arm, and its control moves it
    differently."""
    ctx, groups, dirs, held = _arms_world()
    copy = {k: {g: v.copy() for g, v in dirs.items()} for k in ("zlv", "zw", "zlv_sm", "zw_sm")}
    preds = zr.build_seed(ctx, lambda c: None, dirs, copy, groups, 0.2, held, held)
    for arm, (_, pooled) in zr.NEW_ARMS.items():
        assert preds[arm].X.tobytes() == preds[pooled].X.tobytes(), arm
    direct = zr.gc.m1_r_arms(ctx, dirs, groups, 0.2, held, held)["m1_r_gmom"]
    assert preds["m1_r_gmom"].X.tobytes() == direct.X.tobytes()
    rng = np.random.default_rng(9)
    moved = {**copy, "zlv_sm": {g: (v + rng.normal(0, 0.3, len(v)) if g == "CM" else v) for g, v in dirs.items()}}
    preds = zr.build_seed(ctx, lambda c: None, dirs, moved, groups, 0.2, held, held)
    assert preds["m1_r_gmom_zlv_sm"].X.tobytes() != preds["m1_r_gmom"].X.tobytes()
    assert preds["m1_r_gmom_zlv_sm_perm"].X.tobytes() != preds["m1_r_gmom_zlv_sm"].X.tobytes()
    assert preds["m1_r_gmom_zw_sm_perm"].X.tobytes() == preds["m1_r_gmom"].X.tobytes()


def test_the_fractions(zr):
    fr = zr.arm_fractions(0.2)
    for a in NEW:
        assert fr[a] == (0.2 * 2.0 if a.endswith("_x2") else 0.2), a
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


def _sha(groups=("CM", "Endocardium"), zlv_diff=("CM",), zw_diff=("CM", "Endocardium"), zlv_sm_diff=None,
         zw_sm_diff=None, sm_slope=1.0, tag=""):
    rec = {"zones": ["LA", "LV", "RA", "RV"], "entry": {g: {"sha256": f"e-{g}{tag}"} for g in groups},
           "zw_fallback": []}
    diffs = {"zlv": zlv_diff, "zw": zw_diff, "zlv_sm": zlv_diff if zlv_sm_diff is None else zlv_sm_diff,
             "zw_sm": zw_diff if zw_sm_diff is None else zw_sm_diff}
    for form, diff in diffs.items():
        sm = form.endswith("_sm")
        rec[form] = {g: {"zones": ["LV"], "pooled": False,
                         "slope": (sm_slope if sm else 0.75) if g in diff else 1.0, "pearson": 0.9,
                         "size_ratio": 0.8, "sha256": f"{form}-{g}" if g in diff else f"e-{g}{tag}"}
                     for g in groups}
    rec["changed"] = {f: sorted(d) for f, d in diffs.items()}
    return rec


def _uniform(seed):
    return {a: BASE + g + 0.01 * seed for a, g in UNI.items()}


def _task(truth, gains=None, seeds=range(10), commit="c" * 40, clean=True, fractions=None, identical=(),
          directions=None, anchor_failures=None, same_control_sha=(), same_sha=(), control=None, extra=None):
    """gains: arm -> gain over its pooled arm (controls: over m1_r_gmom) on every seed, or ten gains. An arm with no
    gain gets +d on seed 0 and -d on seed 1, so that it never scores exactly as its pooled arm."""
    gains = gains or {}
    rows = []
    for i, s in enumerate(seeds):
        uni = _uniform(s)
        arms = {"b0": _arm(B0), **{a: _arm(x) for a, x in uni.items()}}
        for j, a in enumerate(NEW):
            if a in gains:
                g = gains[a][i] if isinstance(gains[a], list | tuple) else gains[a]
            else:
                d = 1e-9 * (j + 1)
                g = d if i == 0 else -d if i == 1 else 0.0
            arms[a] = _arm(uni[POOLED_OF[a]] + g)
        for a in identical:
            arms[a] = dict(arms[POOLED_OF[a]])
        sha = {a: f"{a}-{s}" for a in arms}
        for perm in same_control_sha:
            sha[perm] = sha[perm.removesuffix("_perm")]
        for a, b in same_sha:
            sha[a] = sha[b]
        row = {"seed": s, "floor": 50.0, "arms": arms, "sha256": sha}
        if control is not None:
            row["control_directions"] = control(s)
        rows.append(row)
    fr = fractions or {"m1_r_gmom": 0.2, "m1_r_gmom_x15": 0.2 * 1.5, "m1_r_gmom_x2": 0.4,
                       **{a: (0.4 if a.endswith("_x2") else 0.2) for a in NEW}}
    task = "T1/E9.5->E10.5" if truth == "whole" else "T1/E9.5->E10.5:chamber"
    rec = {"label": "synthetic", "task": task, "truth": truth, "fraction": 0.2, "fractions": fr, "seeds": list(seeds),
           "rows": rows, "directions": directions if directions is not None else _sha(),
           "stamp": {"commit": commit, "tree_clean": clean, "utc": "x", "inputs": {}}}
    if anchor_failures is not None:
        rec["anchor_failures"] = anchor_failures
    rec.update(extra or {})
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
BETTER = "the zone-fixed content scores better at equal size along the entry's direction"
WORSE = "the zone-fixed content scores worse at equal size along the entry's direction"
NODIFF = "no difference at the bar"


def _sm_both(form, g, extra=None):
    return {k: {f"m1_r_gmom_{form}_sm": g, f"m1_r_gmom_{form}_sm_x2": g, **(extra or {})} for k in TRUTHS}


@pytest.mark.parametrize("by_truth, zlv, zw", [
    (_sm_both("zlv", UP), BETTER, NODIFF),
    (_sm_both("zlv", DOWN), WORSE, NODIFF),
    (_sm_both("zw", UP), NODIFF, BETTER),
    (_sm_both("zw", DOWN), NODIFF, WORSE),
    ({k: {"m1_r_gmom_zlv_sm": UP} for k in TRUTHS}, NODIFF, NODIFF),
    ({"pooled": {"m1_r_gmom_zlv_sm": UP, "m1_r_gmom_zlv_sm_x2": UP}}, NODIFF, NODIFF),
    ({**{k: {"m1_r_gmom_zlv_sm": UP} for k in TRUTHS},
      **{k: {"m1_r_gmom_zlv_sm_x2": DOWN} for k in ("pooled", "whole")}}, NODIFF, NODIFF),
    ({k: {a: UP for a in RAW} for k in TRUTHS}, NODIFF, NODIFF),  # the raw arms carry no statement
    ({}, NODIFF, NODIFF),
])
def test_each_statement_is_reached(zr, tmp_path, by_truth, zlv, zw):
    r = _read(zr, tmp_path, by_truth)
    assert r["valid"], r["validity"]["reasons"]
    assert r["outcome"] == {"zlv": zlv, "zw": zw}
    txt = (tmp_path / "out" / "t53.txt").read_text()
    assert f"zlv (slope-matched): {zlv}" in txt and f"zw (slope-matched): {zw}" in txt
    for form, s in (("zlv", zlv), ("zw", zw)):
        if s != BETTER:
            assert r["statements"][form]["gene_statement"] is None


def test_the_gene_statement_takes_both_values(zr, tmp_path):
    # the control at the pooled arm: the slope-matched arm is ahead of it too
    r = _read(zr, tmp_path, _sm_both("zlv", UP))
    assert r["statements"]["zlv"]["gene_statement"] == "the gain needs the form's genes"
    assert r["statements"]["zlv"]["gene_comparison"]["test"] == "passes"
    # the control gains as much: not separated
    r = _read(zr, tmp_path / "b", _sm_both("zlv", UP, {"m1_r_gmom_zlv_sm_perm": UP}))
    assert r["outcome"]["zlv"] == BETTER
    assert r["statements"]["zlv"]["gene_statement"] == "a scrambled departure gains as well or is not separated"
    txt = (tmp_path / "b" / "out" / "t53.txt").read_text()
    assert "zlv's genes: a scrambled departure gains as well or is not separated" in txt
    # where the form does not read better the comparison is reported without a statement
    r = _read(zr, tmp_path / "c", {k: {"m1_r_gmom_zw_sm_perm": DOWN} for k in TRUTHS})
    assert r["statements"]["zw"]["gene_statement"] is None
    assert r["statements"]["zw"]["gene_comparison"]["test"] == "passes"
    assert "zw's genes: no statement" in (tmp_path / "c" / "out" / "t53.txt").read_text()


def test_a_gain_of_exactly_point_three_does_not_pass(zr, tmp_path):
    r = _read(zr, tmp_path, _sm_both("zlv", 0.3))
    c = r["statements"]["zlv"]["comparisons"]["m1_r_gmom_zlv_sm"]
    assert c["test"] == "neither" and r["outcome"]["zlv"] == NODIFF
    assert r["statements"]["zlv"]["library_labels"]["lib1"]["m1_r_gmom_zlv_sm"] == "level"


def test_the_library_labels_are_ahead_behind_and_level(zr, tmp_path):
    r = _read(zr, tmp_path, {"lib1": {"m1_r_gmom_zw_sm": UP}, "lib2": {"m1_r_gmom_zw_sm": DOWN}})
    labels = r["statements"]["zw"]["library_labels"]
    assert labels["lib1"]["m1_r_gmom_zw_sm"] == "ahead"
    assert labels["lib2"]["m1_r_gmom_zw_sm"] == "behind"
    assert labels["lib3"]["m1_r_gmom_zw_sm"] == "level"
    assert r["outcome"]["zw"] == NODIFF


def test_the_reported_numbers(zr, tmp_path):
    r = _read(zr, tmp_path, _sm_both("zlv", UP))
    c = r["comparisons"]["whole"]["m1_r_gmom_zlv_sm_x2"]
    assert c["base"] == "m1_r_gmom_x2" and c["mean"] == pytest.approx(1.0) and c["seeds_ahead"] == 10
    assert c["points"]["de"] == pytest.approx(0.5) and c["points"]["mmd_u_variogram"] == pytest.approx(0.5)
    assert set(c["skills"]) >= {"de_score", "mmd_u"}
    assert r["mean_scores"]["lib2"]["m1_r_gmom_x2"] == pytest.approx(BASE + 0.8 + 0.045)
    four = r["statements"]["zlv"]["comparisons"]
    assert set(four) == {"m1_r_gmom_zlv_sm", "m1_r_gmom_zlv_sm_x2"}
    assert all(set(c["cuts"]) == {"pooled", "whole"} for c in four.values())
    g = r["comparisons"]["pooled"]["m1_r_gmom_zlv_sm_vs_perm"]
    assert g["base"] == "m1_r_gmom_zlv_sm_perm" and g["change"] == "m1_r_gmom_zlv_sm" and g["mean"] == pytest.approx(1)
    assert set(r["raw_arms"]) == set(RAW)


def test_what_size_alone_would_give(zr, tmp_path):
    r = _read(zr, tmp_path)
    sa = r["raw_arms"]["m1_r_gmom_zlv"]["size_alone"]["pooled"]
    sa2 = r["raw_arms"]["m1_r_gmom_zlv_x2"]["size_alone"]["pooled"]
    # CM's zlv slope 0.75: at the base fraction outside 1 to 2 (with its reason); at x2 the multiple 1.5, which the
    # uniform arms put 0.4 below x2
    assert sa["CM"]["what_size_alone_would_give"] is None
    assert sa["CM"]["why_not"] == "slope below 1: outside the uniform arms 1 to 2"
    assert sa2["CM"]["what_size_alone_would_give"] == pytest.approx(-0.4) and sa2["CM"]["why_not"] is None
    assert sa["Endocardium"]["what_size_alone_would_give"] == pytest.approx(0.0)
    assert sa2["Endocardium"]["what_size_alone_would_give"] == pytest.approx(0.0)  # multiple 2: the x2 arm itself
    assert zr.size_alone_reason(1.2, "m1_r_gmom_x2") == "twice the slope above 2: outside the uniform arms 1 to 2"
    txt = (tmp_path / "out" / "t53.txt").read_text()
    assert "CM slope 0.750 -> slope below 1: outside the uniform arms 1 to 2" in txt
    means = {"m1_r_gmom": 0.0, "m1_r_gmom_x15": 0.4, "m1_r_gmom_x2": 0.6}
    assert zr.size_alone(means, 1.25, "m1_r_gmom") == pytest.approx(0.2)
    assert zr.size_alone(means, 1.75, "m1_r_gmom_x2") == pytest.approx(-0.1)
    assert zr.size_alone(means, 2.5, "m1_r_gmom") is None and zr.size_alone(means, None, "m1_r_gmom") is None
    assert zr.size_alone_reason(0.4, "m1_r_gmom_x2") == "twice the slope below 1: outside the uniform arms 1 to 2"
    assert zr.size_alone_reason(None, "m1_r_gmom") == "slope not defined"


def test_the_limits_and_choices_are_in_the_text_and_no_count_is_printed(zr, tmp_path):
    _read(zr, tmp_path, _sm_both("zw", UP))
    txt = (tmp_path / "out" / "t53.txt").read_text()
    assert txt.splitlines()[1] == "Choices the brief left open:"
    for x in (*zr.LIMITS, *zr.CHOICES):
        assert x in txt
    head = txt.split("Outcome")[0].split("zlv (slope-matched)")[0]
    for x in ("twice the slope for an _x2 arm", "per group", "ahead, behind", "four parts"):
        assert x in head
    assert "what size alone would give" in txt and "read by nothing" in txt


@pytest.mark.parametrize("over, kw, reason", [
    ({"lib2": None}, {}, "missing truths ['lib2']"),
    ({"whole": {"seeds": range(9)}}, {}, "whole: seeds"),
    ({"lib3": {"commit": "d" * 40}}, {}, "ran at 2 commits"),
    ({"lib1": {"clean": None}}, {}, "tree not clean on ['lib1']"),
    ({"pooled": {"fractions": {"m1_r_gmom": 0.2, "m1_r_gmom_x15": 0.3, "m1_r_gmom_x2": 0.4,
                               **{a: (0.4 if a.endswith("_x2") else 0.2) for a in NEW},
                               "m1_r_gmom_zw_sm_perm": 0.4}}}, {}, "wrong fraction on ['m1_r_gmom_zw_sm_perm']"),
    ({"lib3": {"identical": ("m1_r_gmom_zw_x2",)}}, {}, "score exactly as the pooled arm"),
    ({"lib3": {"identical": ("m1_r_gmom_zlv_sm_perm",)}}, {}, "score exactly as the pooled arm"),
    ({"whole": {"directions": _sha(zlv_diff=("CM", "Endocardium"))}}, {}, "zlv's directions differ"),
    ({"lib2": {"directions": _sha(zlv_sm_diff=("CM", "Endocardium"))}}, {}, "zlv_sm's directions differ"),
    ({"pooled": {"directions": _sha(zlv_diff=())}}, {}, "not in CM alone"),
    ({"lib1": {"directions": _sha(zw_diff=())}}, {}, "zw's directions are the entry's in every group"),
    ({"lib1": {"directions": {}}}, {}, "no direction record"),
    ({"lib3": {"directions": _sha(tag="x")}}, {}, "sha256 differ between the truths' files"),
    ({}, {"directions": _sha(sm_slope=1.0 + 2e-9)}, "recorded slope on the entry's is not within 1e-09 of 1"),
    ({"whole": {"same_control_sha": ("m1_r_gmom_zw_sm_perm",)}}, {}, "controls with their real arm's prediction"),
    ({"lib2": {"anchor_failures": [{"seed": 3, "reason": "x"}]}}, {}, "anchors invert or go missing"),
])
def test_the_void_cases(zr, tmp_path, over, kw, reason):
    r = _read(zr, tmp_path, over=over, **kw)
    assert not r["valid"] and r["outcome"] == "VOID"
    assert any(reason in x for x in r["validity"]["reasons"]), r["validity"]["reasons"]
    assert "Outcome: VOID" in (tmp_path / "out" / "t53.txt").read_text()


def test_a_slope_within_tolerance_of_one_is_valid(zr, tmp_path):
    r = _read(zr, tmp_path, directions=_sha(sm_slope=1.0 + 5e-10))
    assert r["valid"], r["validity"]["reasons"]


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
    assert zr.CONTROL_STREAM == 53


# ---------------------------------------------------------------------------------------------------------------------
# After the review of 636d1bb


def test_zw_is_exactly_zero_in_a_zone_mix_world_with_a_dropped_zone(zr):
    """The review's world: no change inside any zone, LA under the minimum at E14.5 holding a quarter of the E9.5
    weight. Both ends now carry the kept zones' weights, so zw is 0 (it was 1.263 in norm with the pooled E9.5 side)."""
    lv, la = np.array([10.0, 1.0, 5.0, 0.25, 3.0, 7.0]), np.array([1.0, 10.0, 5.0, 4.0, 0.5, 7.0])
    cells = [("E9.5", "LV", "ventricular_cm", "A", 48, lv), ("E9.5", "LA", "atrial_cm", "A", 16, la),
             ("E14.5", "LV", "ventricular_cm", "A", 32, lv), ("E14.5", "LA", "atrial_cm", "A", 9, la)]
    st = _stats(cells)
    entry = _entry(zr, st)
    assert np.nanmax(np.abs(entry["CM"])) > 0.05
    dirs, used, info = zr.zw_directions(st, entry, T1_GENES, CFG)
    assert info["CM"] == {"dropped": ["LA"], "dropped_share": 0.25, "fallback": False} and used["CM"] == ["LV"]
    assert np.all(dirs["CM"][:6] == 0.0) and np.isnan(dirs["CM"][6])
    # the first build's E9.5 side (the entry's pooled mean, LA's cells included) was not free of the zone mix
    labs = ["ventricular_cm", "atrial_cm"]
    old = zr.dl.lfc_ratio(zr.dl.feng_sel(st, "E9.5", labs), zr.zone_sel(st, "E14.5", "LV", labs), 0.5)
    assert np.max(np.abs(old)) > 0.5


def test_zw_is_zero_to_rounding_in_a_random_zone_mix_world_with_a_dropped_zone(zr):
    rng = np.random.default_rng(8)
    prof = {z: rng.gamma(2.0, 3.0, 6) for z in ("LA", "LV", "RA", "RV")}
    n9 = {"LA": 14, "LV": 70, "RA": 12, "RV": 30}
    n14 = {"LA": 40, "LV": 25, "RA": 3, "RV": 44}  # RA under the minimum at E14.5: dropped
    cells = [(s, z, "atrial_cm" if z[1] == "A" else "ventricular_cm", "A", n[z], prof[z])
             for s, n in (("E9.5", n9), ("E14.5", n14)) for z in prof]
    st = _stats(cells)
    entry = _entry(zr, st)
    dirs, used, info = zr.zw_directions(st, entry, T1_GENES, CFG)
    assert info["CM"]["dropped"] == ["RA"] and 0 < info["CM"]["dropped_share"] < 0.5
    assert used["CM"] == ["LA", "LV", "RV"]
    assert np.nanmax(np.abs(entry["CM"])) > 1e-3 and np.nanmax(np.abs(dirs["CM"])) < 1e-12


def test_slope_matching_equalises_the_size_along_the_entrys_direction_only(zr):
    """e = t + m with m orthogonal to t and |m|^2/|t|^2 = 0.5: a form F = t has s = 2/3, and F/s has slope 1 on e but
    squared size 1/s = 1.5 times the entry's: 'equal size' holds along the entry's direction only."""
    rng = np.random.default_rng(2)
    t = rng.normal(size=200)
    m = rng.normal(size=200)
    m -= (m @ t) / (t @ t) * t
    m *= np.sqrt(0.5 * (t @ t) / (m @ m))
    e, genes = t + m, [f"g{i}" for i in range(200)]
    s = zr.compare_directions(t, e, genes, [])["slope"]
    assert s == pytest.approx(2 / 3, abs=1e-12)
    sm = zr.compare_directions(t / s, e, genes, [])
    assert sm["slope"] == pytest.approx(1.0, abs=1e-12) and sm["size_ratio"] == pytest.approx(1 / s, abs=1e-12)
    assert "along the entry's direction" in zr.BETTER and "along the entry's direction" in zr.WORSE
    assert any("1/s" in x and "noise" in x for x in zr.LIMITS)


def test_the_statement_lines_print_one_over_s_and_the_size_alone_figure_at_it(zr, tmp_path):
    r = _read(zr, tmp_path, _sm_both("zlv", UP))
    assert r["outcome"]["zlv"] == BETTER
    g = r["statements"]["zlv"]["groups"]
    assert set(g) == {"CM"} and g["CM"]["inverse_slope"] == pytest.approx(4 / 3) and g["CM"]["sm_size_ratio"] == 0.8
    assert set(r["statements"]["zw"]["groups"]) == {"CM", "Endocardium"}
    sa = r["statements"]["zlv"]["size_alone"]
    # 1/s = 4/3 between the uniform arms 1 (gain 0) and 1.5 (gain 0.4); 2/s = 8/3 is above 2
    assert sa["m1_r_gmom_zlv_sm"]["pooled"]["CM"]["what_size_alone_would_give"] == pytest.approx(0.4 * (1 / 3) / 0.5)
    assert sa["m1_r_gmom_zlv_sm_x2"]["whole"]["CM"]["why_not"] == "2/s above 2: outside the uniform arms 1 to 2"
    assert sa["m1_r_gmom_zlv_sm_x2"]["whole"]["CM"]["what_size_alone_would_give"] is None
    txt = (tmp_path / "out" / "t53.txt").read_text()
    assert f"zlv (slope-matched): {BETTER}" in txt
    assert "CM 1/s 1.333, slope-matched squared size ratio 0.800" in txt
    assert "what size alone would give at the multiple 1/s" in txt and "CM x1.333 -> +0.27" in txt
    assert "what size alone would give at the multiple 2/s" in txt and "CM x2.667 -> 2/s above 2" in txt
    assert "m1_r_gmom_zlv_sm and m1_r_gmom_zlv_x2 are the same build" not in txt


def test_a_slope_matched_arm_that_is_the_raw_x2_build_is_reported_as_one_result(zr, tmp_path):
    r = _read(zr, tmp_path, same_sha=(("m1_r_gmom_zlv_sm", "m1_r_gmom_zlv_x2"),))
    assert r["valid"], r["validity"]["reasons"]
    assert all(r["statements"]["zlv"]["same_build_as_raw_x2"].values())
    assert not any(r["statements"]["zw"]["same_build_as_raw_x2"].values())
    txt = (tmp_path / "out" / "t53.txt").read_text()
    assert "m1_r_gmom_zlv_sm and m1_r_gmom_zlv_x2 are the same build" in txt and "one result, not two" in txt


def test_the_controls_share_one_generator_per_build_seed_and_record_their_slopes(zr):
    from types import SimpleNamespace
    genes = [f"s{i}" for i in range(5)] + [f"o{i}" for i in range(6)]
    scored = genes[:5]
    rng = np.random.default_rng(1)
    e = {"CM": rng.normal(size=11), "Endo": rng.normal(size=11)}
    cm_sm = e["CM"] + rng.normal(size=11)
    forms = {"zlv_sm": {"CM": cm_sm, "Endo": e["Endo"]},
             "zw_sm": {"CM": cm_sm.copy(), "Endo": e["Endo"] + rng.normal(size=11)}}
    ctx = SimpleNamespace(inputs={"ref": SimpleNamespace(var_names=genes)}, task=SimpleNamespace(reference_role="ref"),
                          board=SimpleNamespace(genes=scored), seed=1003)
    rec = {}
    c = zr.seed_controls(ctx, e, forms, [], rec)
    g53 = np.random.default_rng([1003, 53])
    want = {}
    for key, groups in (("zlv_sm", ("CM",)), ("zw_sm", ("CM", "Endo"))):
        want[key] = {}
        for g in groups:
            v, r = e[g].copy(), forms[key][g] - e[g]
            for ix in (np.arange(5), np.arange(5, 11)):
                v[ix] = e[g][ix] + r[ix][g53.permutation(len(ix))]
            want[key][g] = v
    assert np.array_equal(c["zlv_sm_perm"]["CM"], want["zlv_sm"]["CM"])
    assert np.array_equal(c["zw_sm_perm"]["CM"], want["zw_sm"]["CM"])
    assert np.array_equal(c["zw_sm_perm"]["Endo"], want["zw_sm"]["Endo"])
    assert c["zlv_sm_perm"]["Endo"] is e["Endo"]
    # the same CM residual, drawn twice from one stream: the two controls differ
    assert not np.array_equal(c["zlv_sm_perm"]["CM"], c["zw_sm_perm"]["CM"])
    assert set(rec) == {"zlv_sm_perm", "zw_sm_perm"} and set(rec["zlv_sm_perm"]) == {"CM"}
    assert set(rec["zw_sm_perm"]) == {"CM", "Endo"}
    for p, key in (("zlv_sm_perm", "zlv_sm"), ("zw_sm_perm", "zw_sm")):
        for g, x in rec[p].items():
            assert x == zr.compare_directions(want[key][g], e[g], genes, [])
    assert zr.seed_controls(ctx, e, forms, [])["zw_sm_perm"]["CM"].tobytes() == c["zw_sm_perm"]["CM"].tobytes()


def test_a_controls_slope_on_the_entry_is_not_one_where_the_means_are_not_zero(zr):
    """mean(e) 0.3, mean(r) -0.35: the control's slope on e is near 0.9; it is recorded, not assumed to be 1."""
    rng = np.random.default_rng(4)
    n = 4000
    genes = [f"g{i}" for i in range(n)]
    e = {"CM": rng.normal(0.3, 1.0, n)}
    sm = {"CM": e["CM"] + rng.normal(-0.35, 1.0, n)}
    c = zr.control_directions(e, sm, genes, [], genes[:200], zr.control_rng(7))
    s = zr.compare_directions(c["CM"], e["CM"], genes, [])["slope"]
    assert 0.85 < s < 0.95


def test_the_reading_reports_the_controls_slopes(zr, tmp_path):
    def ctl(s):
        return {"zlv_sm_perm": {"CM": {"slope": 0.9 + 0.001 * s}},
                "zw_sm_perm": {"CM": {"slope": 0.95}, "Endocardium": {"slope": 1.01}}}

    r = _read(zr, tmp_path, control=ctl)
    assert r["control_slopes"]["pooled"]["zlv_sm_perm"]["CM"] == pytest.approx([0.9 + 0.001 * s for s in range(10)])
    txt = (tmp_path / "out" / "t53.txt").read_text()
    assert "pooled zlv_sm_perm CM 0.9000 to 0.9090" in txt and "whole zw_sm_perm Endocardium 1.0100 to 1.0100" in txt
    _read(zr, tmp_path / "b")
    assert "Controls' slopes on the entry's direction: not recorded" in (tmp_path / "b" / "out" / "t53.txt").read_text()


def test_the_centred_slope_leaves_out_a_constant_part(zr):
    rng = np.random.default_rng(6)
    e = rng.normal(size=300)
    genes = [f"g{i}" for i in range(300)]
    x = zr.compare_directions(e + 2.0, e, genes, [])
    assert abs(x["slope"] - 1.0) > 0.01 and x["slope_centred"] == pytest.approx(1.0, abs=1e-12)
    rec = zr.direction_record({"CM": e}, {"zlv": ({"CM": e + 2.0}, {"CM": ["LV"]})}, genes, [], ["LV"])
    assert rec["zlv"]["CM"]["slope_centred"] == pytest.approx(1.0, abs=1e-12)


def test_slope_matching_names_every_refused_group(zr):
    entry = {"A": np.array([1.0, -2.0, 3.0]), "B": np.array([2.0, 1.0, -1.0]), "C": np.ones(3)}
    form = {"A": entry["A"] * 0.1, "B": -entry["B"], "C": entry["C"]}
    with pytest.raises(zr.SlopeRefusal) as e:
        zr.slope_matched(entry, form, ["a", "b", "c"], [], "zw")
    assert e.value.slopes == {"zw": {"A": pytest.approx(0.1), "B": pytest.approx(-1.0)}}
    assert "not read; the form is nearly orthogonal to the entry's direction" in str(e.value)


def test_a_slope_refusal_carries_both_forms_slopes_and_the_raw_record(zr, monkeypatch, tmp_path):
    st = _rng_world()
    entry = _entry(zr, st)
    monkeypatch.setattr(zr, "MIN_SLOPE", 1e9)
    with pytest.raises(zr.SlopeRefusal) as e:
        zr.all_directions(st, entry, T1_GENES, CFG, [])
    assert set(e.value.slopes) == {"zlv", "zw"} and set(e.value.slopes["zlv"]) == {"CM"}
    assert set(e.value.record) >= {"entry", "zlv", "zw"} and "zlv_sm" not in e.value.record
    res = {"label": "x", "task": "T1/E9.5->E10.5", "truth": "whole", "seeds": [0], **zr.slope_refusal_record(e.value),
           "rows": [], "stamp": dict(STAMP)}
    zr.write_task_file(tmp_path / "o", "whole.json", res)
    back = json.loads((tmp_path / "o" / "tasks" / "whole.json").read_text())
    assert back["slope_refusal"]["reason"] == zr.SLOPE_REFUSED and "\"n\"" not in json.dumps(back)
    with pytest.raises(SystemExit, match="already holds"):
        zr.gv.refuse_existing(tmp_path / "o" / "tasks", ("whole.json",))


def test_a_slope_refusal_voids_the_read_with_its_reason(zr, tmp_path):
    refusal = {"slope_refusal": {"reason": "x", "min_slope": 0.2, "slopes": {"zw": {"CM": 0.006}}}}
    r = _read(zr, tmp_path, over={"whole": {"extra": refusal}})
    assert not r["valid"]
    assert ("whole: not read; the form is nearly orthogonal to the entry's direction: zw in CM (slope 0.0060, not "
            "above 0.2)") in r["validity"]["reasons"]


def test_split_half_reliability_is_one_with_identical_halves_and_recorded(zr):
    rng = np.random.default_rng(12)
    cells = []
    for stage in ("E9.5", "E14.5"):
        for zone in ("LA", "LV", "RV"):
            for typ in ("ventricular_cm", "atrial_cm", "blood"):
                k, prof = int(rng.integers(10, 40)), rng.gamma(2.0, 3.0, 6)
                cells += [(stage, zone, typ, "A", k, prof), (stage, zone, typ, "B", k, prof)]
    st = _stats(cells)
    entry = _entry(zr, st)
    _, rec = zr.all_directions(st, entry, T1_GENES, CFG, [])
    assert rec["entry"]["CM"]["split_half_reliability"] == 1.0
    assert rec["zlv"]["CM"]["split_half_reliability"] == 1.0
    assert rec["changed"]["zw"] and all(rec["zw"][g]["split_half_reliability"] == 1.0 for g in rec["changed"]["zw"])
    assert "split_half_reliability" not in rec["zlv_sm"]["CM"]


def test_split_half_reliability_is_below_one_with_noisy_halves_and_none_without_a_half(zr):
    st = _rng_world()
    entry = _entry(zr, st)
    within = zr.within_zone_directions(st, T1_GENES, CFG, entry)
    raw = {"zlv": zr.zlv_directions(entry, within), "zw": zr.zw_directions(st, entry, T1_GENES, CFG)[:2]}
    rel = zr.split_half_reliability(st, entry, raw, T1_GENES, CFG, [])
    assert set(rel) == {"entry", "zlv", "zw"} and set(rel["zlv"]) == {"CM"}
    vals = [x for per in rel.values() for x in per.values()]
    assert all(x is None or 0.0 <= x <= 1.0 for x in vals) and any(x is not None and x < 1.0 for x in vals)
    keep = np.array([k.split("|")[3] == "A" for k in st["keys"]])
    only_a = {k: (v[keep] if k != "genes" else v) for k, v in st.items()}
    e_a = _entry(zr, only_a)
    raw_a = {"zlv": zr.zlv_directions(e_a, zr.within_zone_directions(only_a, T1_GENES, CFG, e_a)),
             "zw": zr.zw_directions(only_a, e_a, T1_GENES, CFG)[:2]}
    rel_a = zr.split_half_reliability(only_a, e_a, raw_a, T1_GENES, CFG, [])
    assert all(x is None for per in rel_a.values() for x in per.values())
