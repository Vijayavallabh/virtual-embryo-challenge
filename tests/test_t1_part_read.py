"""experiments/t1_m1/part_read.py (try T55 of experiments/harness/TRIES.md, brief AI) on synthetic data: the forms
a_g c + b s_g built from try T51's split, the per-group multiples, the refusals, the new arms and the reading. No
released file, statistics file, board file or score is read."""

from __future__ import annotations

import importlib.util
import json
import sys

import numpy as np
import pytest
from test_t1_step_read import _arms_world

from vec.boards import GROUPS, REPO_ROOT

SCRIPT = REPO_ROOT / "experiments" / "t1_m1" / "part_read.py"
STAMP = {"commit": "f" * 40, "tree_clean": True, "utc": "x", "inputs": {}}
NEW = ("m1_r_gmom_c1", "m1_r_gmom_s05", "m1_r_gmom_s15", "m1_r_gmom_cn", "m1_r_gmom_cn_x2", "m1_r_gmom_cp",
       "m1_r_gmom_cp_x2")
POOLED_OF = {a: ("m1_r_gmom_x2" if a.endswith("_x2") else "m1_r_gmom") for a in NEW}
TRUTHS = ("pooled", "lib1", "lib2", "lib3", "whole")


@pytest.fixture
def pr(monkeypatch, tmp_path):
    spec = importlib.util.spec_from_file_location("part_read_module", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    monkeypatch.setattr(mod, "RESULTS", tmp_path / "results")
    monkeypatch.setattr(mod.gv, "stamp", lambda started, files: dict(STAMP))
    return mod


# ---------------------------------------------------------------------------------------------------------------------
# The directions


GENES = [f"g{i}" for i in range(40)] + ["h", "n"]  # h held, n NaN in every group


def _entry(seed=0, groups=("CM", "Endocardium", "Endothelium", "Blood")):
    """Random directions sharing a common part (so sum(c d_g) > 0); h finite but held, n NaN everywhere."""
    rng = np.random.default_rng(seed)
    common = rng.normal(0, 1.0, 40)
    out = {}
    for g in groups:
        v = np.append(common + rng.normal(0, 0.8, 40), [rng.normal(), np.nan])
        out[g] = v
    return out


def test_one_one_rebuilds_the_entry_and_one_zero_is_c(pr):
    e = _entry()
    forms, rec = pr.all_directions(e, GENES, ["h"])
    sp = pr.split(e, GENES, ["h"])
    M = sp["M"]
    assert M.sum() == 40 and not M[40] and not M[41]
    c = np.mean([e[g][:40] for g in sorted(e)], axis=0)
    for g in e:
        one = pr.form_direction(e[g], sp["c"], sp["s"][g], M, 1.0, 1.0)
        assert np.max(np.abs(one[M] - e[g][M])) <= 1e-12
        assert rec["entry"][g]["rebuild_11_max_abs"] <= 1e-12
        assert np.allclose(forms["c1"][g][:40], c, rtol=0, atol=1e-12)
        assert np.array_equal(forms["c1"][g][:40], sp["c"][:40])
        # off M the entry's own values: the held gene and the NaN gene
        assert forms["c1"][g][40] == e[g][40] and np.isnan(forms["c1"][g][41])
        assert np.allclose(forms["s05"][g][M], c + 0.5 * (e[g][M] - c), rtol=0, atol=1e-12)
        assert np.allclose(forms["s15"][g][M], c + 1.5 * (e[g][M] - c), rtol=0, atol=1e-12)
    # c1 is the same vector in every group on M
    gs = sorted(e)
    assert all(np.array_equal(forms["c1"][g][M], forms["c1"][gs[0]][M]) for g in gs)


def test_a_n_gives_the_entrys_squared_size_and_a_p_its_projection(pr):
    e = _entry(3)
    forms, rec = pr.all_directions(e, GENES, ["h"])
    M = pr.split(e, GENES, ["h"])["M"]
    for g in e:
        d = e[g][M]
        cn, cp = forms["cn"][g][M], forms["cp"][g][M]
        assert float(cn @ cn) == pytest.approx(float(d @ d), rel=0, abs=1e-12 * float(d @ d))
        # projection on d_g: <cp, d> / |d| equals |d|, i.e. <cp, d> = <d, d>
        assert float(cp @ d) == pytest.approx(float(d @ d), rel=0, abs=1e-12 * float(d @ d))
        assert rec["cn"][g]["size_ratio"] == pytest.approx(1.0, abs=1e-12)
        assert rec["cp"][g]["slope"] == pytest.approx(1.0, abs=1e-12)
        assert rec["c1"][g]["slope"] == pytest.approx(1 / rec["entry"][g]["a_p"], abs=1e-12)
        assert rec["entry"][g]["a_p"] >= rec["entry"][g]["a_n"] - 1e-12  # Cauchy-Schwarz
        ep = rec["entry"][g]
        assert ep["e_c"] + ep["e_s"] + ep["cross"] == pytest.approx(1.0, abs=1e-12)
        assert ep["e_c"] == pytest.approx(1 / ep["a_n"] ** 2, abs=1e-12)
    text = json.dumps(rec)
    assert "\"n\"" not in text and "cells" not in text


def test_the_run_refuses_a_group_whose_moved_genes_are_not_m(pr):
    e = _entry()
    e["Blood"] = e["Blood"].copy()
    e["Blood"][5] = np.nan  # Blood moves one gene fewer than the others
    with pytest.raises(ValueError, match=r"moved genes of \['CM', 'Endocardium', 'Endothelium'\] are not M"):
        pr.all_directions(e, GENES, ["h"])
    # a gene held in one sense only cannot arise (the held list is common), but a gene finite in one group only can
    e = _entry()
    e["CM"] = e["CM"].copy()
    e["CM"][41] = 1.0
    with pytest.raises(ValueError, match=r"moved genes of \['CM'\] are not M"):
        pr.all_directions(e, GENES, ["h"])


def test_the_run_refuses_a_group_with_sum_c_d_not_positive(pr):
    rng = np.random.default_rng(1)
    c = rng.normal(size=40)
    e = {g: np.append(c + 0.1 * rng.normal(size=40), [0.0, np.nan]) for g in ("A", "B", "D")}
    # C opposite to the others: the mean is about c / 2, so sum(c d_C) < 0 while the others stay positive
    e["C"] = np.append(-c, [0.0, np.nan])
    with pytest.raises(ValueError, match=r"sum\(c d_g\) is not positive in \['C'\]"):
        pr.all_directions(e, GENES, ["h"])
    # exactly zero: the groups cancel on M, so c = 0 and every sum(c d_g) is 0
    z = {"A": np.append(c, [0.0, np.nan]), "B": np.append(-c, [0.0, np.nan])}
    with pytest.raises(ValueError, match="not positive in \\['A', 'B'\\]"):
        pr.all_directions(z, GENES, ["h"])


def test_the_split_is_try_t51s(pr):
    e = _entry(2)
    sp = pr.split(e, GENES, ["h"])
    held_set = {"h"}
    masks = {g: np.isfinite(e[g]) & np.array([x not in held_set for x in GENES]) for g in e}
    M = pr.split_read.common_mask(masks)
    sv = pr.split_read.split_vectors(e, M)
    assert np.array_equal(sp["M"], M)
    assert np.array_equal(sp["c"], sv["c"], equal_nan=True)
    assert all(np.array_equal(sp["s"][g], sv["s"][g], equal_nan=True) for g in e)


# ---------------------------------------------------------------------------------------------------------------------
# The arms


def test_a_new_arm_differs_from_its_pooled_arms_build_in_the_directions_only(pr, monkeypatch):
    calls = []

    def spy(ctx, directions, groups, fraction, held_m1, held_factors):
        calls.append({"ctx": ctx, "directions": directions, "groups": groups, "fraction": fraction,
                      "held_m1": held_m1, "held_factors": held_factors})
        return {"m1_r": ("m1_r", fraction), "m1_r_gmom": ("m1_r_gmom", id(directions), fraction)}

    ctx, entry, groups, held_m1, held_f = object(), {"CM": np.zeros(3)}, {"V-CM": "CM"}, ["a", "b"], ["a"]
    forms = {k: {"CM": np.full(3, float(i + 1))} for i, k in enumerate(pr.FORMS)}
    monkeypatch.setattr(pr.gc, "m1_r_arms", spy)
    preds = pr.build_seed(ctx, lambda c: "b0", entry, forms, groups, 0.2, held_m1, held_f)
    assert set(preds) == set(pr.ARMS) and len(pr.NEW_ARMS) == 7
    by_fraction = {c["fraction"]: c for c in calls if c["directions"] is entry}
    assert sorted(by_fraction) == pytest.approx([0.2, 0.3, 0.4])
    news = [c for c in calls if c["directions"] is not entry]
    assert len(news) == 7
    for c, (arm, (key, pooled)) in zip(news, pr.NEW_ARMS.items(), strict=True):
        base = by_fraction[0.4 if pooled == "m1_r_gmom_x2" else 0.2]
        assert pooled == POOLED_OF[arm]
        assert c["fraction"] == base["fraction"]
        assert c["directions"] is forms[key]
        for k in ("ctx", "groups", "held_m1", "held_factors"):
            assert c[k] is base[k], f"{arm} differs from its pooled arm's build in {k}"
    assert {a: k for a, (k, _) in pr.NEW_ARMS.items()} == {
        "m1_r_gmom_c1": "c1", "m1_r_gmom_s05": "s05", "m1_r_gmom_s15": "s15", "m1_r_gmom_cn": "cn",
        "m1_r_gmom_cn_x2": "cn", "m1_r_gmom_cp": "cp", "m1_r_gmom_cp_x2": "cp"}
    assert pr.FORMS == {"c1": (1.0, 0.0), "s05": (1.0, 0.5), "s15": (1.0, 1.5), "cn": ("a_n", 0.0),
                        "cp": ("a_p", 0.0)}


def test_a_form_equal_to_the_entry_builds_the_pooled_arms_bit_for_bit(pr):
    """The idle case on a real build: forms whose directions are the entry's (copies) give the pooled arms bit for bit;
    a changed form moves its arm."""
    ctx, groups, dirs, held = _arms_world()
    copy = {k: {g: v.copy() for g, v in dirs.items()} for k in pr.FORMS}
    preds = pr.build_seed(ctx, lambda c: None, dirs, copy, groups, 0.2, held, held)
    for arm, (_, pooled) in pr.NEW_ARMS.items():
        assert preds[arm].X.tobytes() == preds[pooled].X.tobytes(), arm
    direct = pr.gc.m1_r_arms(ctx, dirs, groups, 0.2, held, held)["m1_r_gmom"]
    assert preds["m1_r_gmom"].X.tobytes() == direct.X.tobytes()
    moved = {**copy, "c1": {g: v * 0.5 for g, v in dirs.items()}}
    preds = pr.build_seed(ctx, lambda c: None, dirs, moved, groups, 0.2, held, held)
    assert preds["m1_r_gmom_c1"].X.tobytes() != preds["m1_r_gmom"].X.tobytes()


def test_the_fractions(pr):
    fr = pr.arm_fractions(0.2)
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


def _dirs(groups=("CM", "Endocardium"), tag="", a_n=1.25, a_p=1.6, rebuild=1e-15, cp_size=1.44):
    rec = {"entry": {g: {"sha256": f"e-{g}{tag}", "a_n": a_n, "a_p": a_p, "e_c": 0.64, "e_s": 0.3, "cross": 0.06,
                         "rebuild_11_max_abs": rebuild} for g in groups},
           "forms": {}}
    for f in ("c1", "s05", "s15", "cn", "cp"):
        rec[f] = {g: {"slope": 0.6, "pearson": 0.8, "slope_centred": 0.6,
                      "size_ratio": cp_size if f == "cp" else 1.0 if f == "cn" else 0.64,
                      "sha256": f"{f}-{g}"} for g in groups}
    return rec


def _uniform(seed):
    return {a: BASE + g + 0.01 * seed for a, g in UNI.items()}


def _task(truth, gains=None, seeds=range(10), commit="c" * 40, clean=True, fractions=None, identical=(),
          directions=None, anchor_failures=None):
    """gains: arm -> gain over its pooled arm on every seed. An arm with no gain gets +d on seed 0 and -d on seed 1,
    so that it never scores exactly as its pooled arm."""
    gains = gains or {}
    rows = []
    for i, s in enumerate(seeds):
        uni = _uniform(s)
        arms = {"b0": _arm(B0), **{a: _arm(x) for a, x in uni.items()}}
        for j, a in enumerate(NEW):
            if a in gains:
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
           "rows": rows, "directions": directions if directions is not None else _dirs(),
           "stamp": {"commit": commit, "tree_clean": clean, "utc": "x", "inputs": {}}}
    if anchor_failures is not None:
        rec["anchor_failures"] = anchor_failures
    return rec


def _read(pr, tmp_path, by_truth=None, over=None, committed_shift=None, **kw):
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
        rows = []
        for s in range(10):
            arms = {a: _arm(x + committed_shift.get((k, a), 0.0)) for a, x in {"b0": B0, **_uniform(s)}.items()}
            rows.append({"seed": s, "arms": arms})
        files[k] = tmp_path / f"committed_{k}.json"
        files[k].write_text(json.dumps({"rows": rows}))
    return pr.read(out, files)["reading"]


UP, DOWN = 1.0, -1.0
NODIFF = "no difference at the bar"
DROP_BETTER, CARRY = "dropping them scores better", "they carry part of the score"
MOVE_BETTER = "the common part alone scores better at equal movement"
MOVE_WORSE = "the common part alone scores worse at equal movement"
PROJ_BETTER = "the common part alone scores better at equal projection"
PROJ_WORSE = "the common part alone scores worse at equal projection"


def _all(gains):
    return {k: dict(gains) for k in TRUTHS}


@pytest.mark.parametrize("gains, s1", [
    ({"m1_r_gmom_c1": UP}, {"m1_r_gmom_c1": DROP_BETTER, "m1_r_gmom_s05": NODIFF, "m1_r_gmom_s15": NODIFF}),
    ({"m1_r_gmom_c1": DOWN}, {"m1_r_gmom_c1": CARRY, "m1_r_gmom_s05": NODIFF, "m1_r_gmom_s15": NODIFF}),
    ({"m1_r_gmom_s05": UP, "m1_r_gmom_s15": DOWN},
     {"m1_r_gmom_c1": NODIFF, "m1_r_gmom_s05": DROP_BETTER, "m1_r_gmom_s15": CARRY}),
    ({}, {"m1_r_gmom_c1": NODIFF, "m1_r_gmom_s05": NODIFF, "m1_r_gmom_s15": NODIFF}),
])
def test_statement_one_takes_every_value(pr, tmp_path, gains, s1):
    r = _read(pr, tmp_path, _all(gains))
    assert r["valid"], r["validity"]["reasons"]
    assert r["outcome"]["statement_1"] == s1
    txt = (tmp_path / "out" / "t55.txt").read_text()
    for a, s in s1.items():
        assert f"  {a}: {s}" in txt


@pytest.mark.parametrize("gains, s2", [
    ({"m1_r_gmom_cn": UP, "m1_r_gmom_cn_x2": UP}, MOVE_BETTER),
    ({"m1_r_gmom_cn": DOWN, "m1_r_gmom_cn_x2": DOWN}, MOVE_WORSE),
    ({}, NODIFF),
    ({"m1_r_gmom_cn": UP}, "mixed (m1_r_gmom_cn: passes; m1_r_gmom_cn_x2: neither)"),
    ({"m1_r_gmom_cn": UP, "m1_r_gmom_cn_x2": DOWN}, "mixed (m1_r_gmom_cn: passes; m1_r_gmom_cn_x2: loses)"),
])
def test_statement_two_takes_every_value(pr, tmp_path, gains, s2):
    r = _read(pr, tmp_path, _all(gains))
    assert r["valid"], r["validity"]["reasons"]
    assert r["outcome"]["statement_2"] == s2
    assert f"Statement 2, content at equal movement: {s2}" in (tmp_path / "out" / "t55.txt").read_text()


@pytest.mark.parametrize("gains, s3", [
    ({"m1_r_gmom_cp": UP, "m1_r_gmom_cp_x2": UP}, PROJ_BETTER),
    ({"m1_r_gmom_cp": DOWN, "m1_r_gmom_cp_x2": DOWN}, PROJ_WORSE),
    ({}, NODIFF),
    ({"m1_r_gmom_cp_x2": DOWN}, "mixed (m1_r_gmom_cp: neither; m1_r_gmom_cp_x2: loses)"),
])
def test_statement_three_takes_every_value_and_is_read_by_nothing(pr, tmp_path, gains, s3):
    r = _read(pr, tmp_path, _all(gains))
    assert r["outcome"]["statement_3"] == s3 and r["statement_3"]["read_by_nothing"] is True
    assert "Statement 3, content at equal projection (read by nothing): " + s3 in (
        tmp_path / "out" / "t55.txt").read_text()


def test_statement_two_only_on_both_cuts(pr, tmp_path):
    r = _read(pr, tmp_path, {"pooled": {"m1_r_gmom_cn": UP, "m1_r_gmom_cn_x2": UP}})
    assert r["outcome"]["statement_2"] == NODIFF


def test_a_gain_of_exactly_point_three_does_not_pass(pr, tmp_path):
    r = _read(pr, tmp_path, _all({"m1_r_gmom_c1": 0.3, "m1_r_gmom_cn": 0.3, "m1_r_gmom_cn_x2": 0.3}))
    assert r["statement_1"]["m1_r_gmom_c1"]["comparison"]["test"] == "neither"
    assert r["outcome"]["statement_1"]["m1_r_gmom_c1"] == NODIFF and r["outcome"]["statement_2"] == NODIFF
    assert r["statement_1"]["m1_r_gmom_c1"]["library_labels"]["lib1"] == "level"


def test_the_library_labels_are_ahead_behind_and_level(pr, tmp_path):
    r = _read(pr, tmp_path, {"lib1": {"m1_r_gmom_c1": UP}, "lib2": {"m1_r_gmom_c1": DOWN}})
    labels = r["statement_1"]["m1_r_gmom_c1"]["library_labels"]
    assert labels == {"lib1": "ahead", "lib2": "behind", "lib3": "level"}
    assert r["outcome"]["statement_1"]["m1_r_gmom_c1"] == NODIFF


def test_the_reported_numbers(pr, tmp_path):
    r = _read(pr, tmp_path, _all({"m1_r_gmom_cn_x2": UP}))
    c = r["comparisons"]["whole"]["m1_r_gmom_cn_x2"]
    assert c["base"] == "m1_r_gmom_x2" and c["mean"] == pytest.approx(1.0) and c["seeds_ahead"] == 10
    assert c["points"]["de"] == pytest.approx(0.5) and c["points"]["mmd_u_variogram"] == pytest.approx(0.5)
    assert r["comparisons"]["pooled"]["m1_r_gmom_c1"]["base"] == "m1_r_gmom"
    assert r["mean_scores"]["lib2"]["m1_r_gmom_x2"] == pytest.approx(BASE + 0.8 + 0.045)
    assert set(r["comparisons"]["lib3"]) == set(NEW)
    txt = (tmp_path / "out" / "t55.txt").read_text()
    assert "truth lib3: mean scores b0 50.00" in txt and "de terms" in txt and "mmd_u and variogram" in txt
    assert "a_n 1.2500, a_p 1.6000" in txt


def test_statement_three_prints_the_uniform_arms_at_the_same_squared_size(pr, tmp_path):
    # cp's squared size ratio 1.44: the multiple 1.2 for m1_r_gmom_cp (between 1 and 1.5: 0.4 x 0.2 / 0.5 = 0.16),
    # 2.4 for m1_r_gmom_cp_x2 (outside the uniform arms)
    r = _read(pr, tmp_path)
    ss = r["statement_3"]["same_size"]
    v = ss["m1_r_gmom_cp"]["pooled"]["CM"]
    assert v["multiple"] == pytest.approx(1.2) and v["what_size_alone_would_give"] == pytest.approx(0.16)
    w = ss["m1_r_gmom_cp_x2"]["whole"]["Endocardium"]
    assert w["multiple"] == pytest.approx(2.4) and w["what_size_alone_would_give"] is None
    assert "outside the uniform arms" in w["why_not"]
    txt = (tmp_path / "out" / "t55.txt").read_text()
    assert "CM x1.200 -> +0.16" in txt and "outside the uniform arms" in txt
    assert pr.same_size({"m1_r_gmom": 0, "m1_r_gmom_x15": 0.4, "m1_r_gmom_x2": 0.8}, 0.5, "m1_r_gmom")["why_not"] \
        == "the same-size multiple below 1: outside the uniform arms 1 to 2"


def test_the_limits_and_choices_are_in_the_text(pr, tmp_path):
    _read(pr, tmp_path)
    txt = (tmp_path / "out" / "t55.txt").read_text()
    assert txt.splitlines()[1] == "Choices the brief left open:"
    for x in (*pr.LIMITS, *pr.CHOICES):
        assert x in txt
    assert "one result read twice" in txt and "equal squared size of the direction" in txt


@pytest.mark.parametrize("over, kw, reason", [
    ({"lib2": None}, {}, "missing truths ['lib2']"),
    ({"whole": {"seeds": range(9)}}, {}, "whole: seeds"),
    ({"lib3": {"commit": "d" * 40}}, {}, "ran at 2 commits"),
    ({"lib1": {"clean": False}}, {}, "tree not clean on ['lib1']"),
    ({"pooled": {"fractions": {"m1_r_gmom": 0.2, "m1_r_gmom_x15": 0.3, "m1_r_gmom_x2": 0.4,
                               **{a: (0.4 if a.endswith("_x2") else 0.2) for a in NEW}, "m1_r_gmom_cp_x2": 0.2}}},
     {}, "wrong fraction on ['m1_r_gmom_cp_x2']"),
    ({"lib3": {"identical": ("m1_r_gmom_cn_x2",)}}, {}, "score exactly as the pooled arm"),
    ({"lib3": {"directions": _dirs(tag="x")}}, {}, "differ between the truths' files"),
    ({"lib1": {"directions": _dirs(a_n=1.2500001)}}, {}, "differ between the truths' files"),
    ({"whole": {"directions": _dirs(a_p=1.61)}}, {}, "differ between the truths' files"),
    ({}, {"directions": _dirs(rebuild=2e-12)}, "the form (1, 1) is not the entry's direction within 1e-12"),
    ({"lib1": {"directions": {}}}, {}, "no full direction record"),
    ({"lib2": {"anchor_failures": [{"seed": 3, "reason": "x"}]}}, {}, "anchors invert or go missing"),
])
def test_the_void_cases(pr, tmp_path, over, kw, reason):
    r = _read(pr, tmp_path, over=over, **kw)
    assert not r["valid"] and r["outcome"] == "VOID"
    assert any(reason in x for x in r["validity"]["reasons"]), r["validity"]["reasons"]
    assert "Outcome: VOID" in (tmp_path / "out" / "t55.txt").read_text()


@pytest.mark.parametrize("truth, arm", [("pooled", "b0"), ("lib2", "m1_r_gmom"), ("lib3", "m1_r_gmom_x15"),
                                        ("whole", "m1_r_gmom_x2")])
def test_a_rescored_arm_off_its_committed_score_by_1e_6_voids(pr, tmp_path, truth, arm):
    r = _read(pr, tmp_path, committed_shift={(truth, arm): 1e-6})
    assert not r["valid"] and any(x.startswith(f"{truth}: ") and "reproduce" in x for x in r["validity"]["reasons"])
    assert _read(pr, tmp_path / "ok", committed_shift={(truth, arm): 1e-10})["valid"]


def test_the_reading_refuses_to_overwrite(pr, tmp_path):
    _read(pr, tmp_path)
    with pytest.raises(SystemExit, match="already holds"):
        pr.read(tmp_path / "out", {})


# ---------------------------------------------------------------------------------------------------------------------
# The run's refusals


def test_the_run_needs_its_row(pr, monkeypatch):
    monkeypatch.setattr(pr.n12, "tries_rows", lambda: {"T53", "T54"})
    monkeypatch.setattr(pr.gc, "n22_passed", lambda path: True)
    with pytest.raises(SystemExit, match="T55 has no row"):
        pr.main(["run", "--truth", "whole", "--out", "t55"])


def test_the_run_needs_n22(pr, monkeypatch):
    monkeypatch.setattr(pr.n12, "tries_rows", lambda: {"T55"})
    with pytest.raises(SystemExit, match="does not hold a passing N22"):
        pr.main(["run", "--truth", "pooled", "--out", "t55", "--n22", "n22x"])


@pytest.mark.parametrize("truth", TRUTHS)
def test_the_run_refuses_an_existing_task_file(pr, monkeypatch, truth):
    monkeypatch.setattr(pr.n12, "tries_rows", lambda: {"T55"})
    monkeypatch.setattr(pr.gc, "n22_passed", lambda path: True)
    tasks = pr.RESULTS / "t55" / "tasks"
    tasks.mkdir(parents=True)
    (tasks / f"{truth}.json").write_text("{}")
    with pytest.raises(SystemExit, match="already holds"):
        pr.main(["run", "--truth", truth, "--out", "t55"])


def test_the_truths_go_to_their_cuts(pr):
    assert pr.TASK_OF == {"pooled": "T1/E9.5->E10.5:chamber", "lib1": "T1/E9.5->E10.5:chamber",
                          "lib2": "T1/E9.5->E10.5:chamber", "lib3": "T1/E9.5->E10.5:chamber",
                          "whole": "T1/E9.5->E10.5"}
    assert pr.COMMITTED["lib2"].parts[-3:] == ("t46_m1_lib_step", "tasks", "lib2.json")
    assert pr.COMMITTED["whole"].parts[-3:] == ("t44_m1_step", "tasks", "T1_E9.5_E10.5.json")
