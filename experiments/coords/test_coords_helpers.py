"""Synthetic checks of the coords helpers (no competition data). Run: uv run pytest experiments/coords -q"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import coordlib as cl  # noqa: E402
import handedness as hd  # noqa: E402
import zpitch_fit as zp  # noqa: E402


def test_levels_are_paired_by_z_difference_not_index():
    levels = np.array([0.0, 10.0, 20.0, 50.0, 60.0])  # step 10, two levels missing between 20 and 50
    assert cl.level_step(levels) == 10.0
    assert cl.level_partners(levels, 1).tolist() == [1, 2, -1, 4, -1]
    assert cl.level_partners(levels, 3).tolist() == [-1, -1, 3, -1, -1]
    assert cl.level_partners(levels + 1e-9 * np.arange(5), 1).tolist() == [1, 2, -1, 4, -1]  # float noise


def test_vertical_offset_has_no_partner_across_a_gap():
    C = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 0.0, 10.0], [3.0, 0.0, 40.0]])
    lev, levels, _ = cl.snap_z(C[:, 2])
    vo = cl.vertical_offset(C, lev, levels, dk=1, step=10.0)
    assert vo[:2] == pytest.approx([0.5, 0.5])
    assert np.isnan(vo[2]) and np.isnan(vo[3])  # 10 -> 20 is empty; 40 -> 50 does not exist


def test_stage_rng_is_per_stage_and_reproducible():
    a1, a2 = cl.stage_rng("E9.5").random(3), cl.stage_rng("E9.5").random(3)
    np.testing.assert_array_equal(a1, a2)
    assert not np.allclose(a1, cl.stage_rng("E8.75").random(3))
    assert not np.allclose(a1, cl.stage_rng("E9.5", seed=1).random(3))


def test_multi_section_levels_and_section_ids():
    import anndata as ad

    a = ad.AnnData(obs=pd.DataFrame({"sample_id": ["30", "30", "35", "12"]}, index=list("abcd")))
    lev = np.array([0, 0, 0, 1])
    assert cl.multi_section_levels(a, lev) == [0]
    assert cl.section_ids(a, lev).tolist() == ["0:30", "0:30", "0:35", "1:12"]
    b = ad.AnnData(obs=pd.DataFrame(index=list("ab")))
    assert cl.multi_section_levels(b, np.array([0, 0])) == []


def _synthetic_pitch(p_true: float):
    dgrid = np.arange(5.0, 150.0, 2.0)
    gfit = 1 - np.exp(-dgrid / 20.0)
    dxy = np.arange(1.0, 40.0, 2.0)
    kk = np.repeat([1, 2, 3], dxy.size)
    dxy2 = np.tile(dxy ** 2, 3)
    obs = 1 - np.exp(-np.sqrt(dxy2 + (kk * p_true) ** 2) / 20.0)
    return dgrid, gfit, obs, np.ones_like(obs), dxy2, kk


def test_pitch_fit_recovers_an_identified_pitch_with_a_closed_range():
    f = zp.fit_pitch(*_synthetic_pitch(12.0))
    assert f["p_best"] == pytest.approx(12.0, abs=0.25)
    assert not f["p_lo_open"] and not f["p_hi_open"] and not f["p_best_in_clamped_zone"]
    assert f["range_p"].startswith("[")


def test_pitch_below_the_in_plane_curve_is_reported_open_ended():
    f = zp.fit_pitch(*_synthetic_pitch(3.2))
    assert f["p_clamp"] == pytest.approx(np.sqrt(5.0 ** 2 - 1.0), abs=1e-9)
    assert f["p_lo_open"] and f["range_p"].startswith("(open")


def tripod(n=1500, seed=0):
    """Three orthogonal arms of lengths 1, 2, 3 from the origin: chiral (arms are told apart by length)."""
    rng = np.random.default_rng(seed)
    arm = rng.integers(0, 3, n)
    pts = np.zeros((n, 3))
    pts[np.arange(n), arm] = rng.uniform(0, 1, n) * (arm + 1)
    return pts + rng.normal(scale=0.03, size=(n, 3))


def test_shape_alignment_tells_proper_from_improper():
    rng = np.random.default_rng(1)
    A = hd.normalise(tripod(seed=2))
    R = hd.random_rotations(1, rng)[1]
    B = hd.normalise(tripod(seed=3) @ R.T)
    dp, _ = hd.best_fit(A, B, False, np.random.default_rng(4), n_rot=300)
    di, _ = hd.best_fit(A, B, True, np.random.default_rng(4), n_rot=300)
    assert dp < di
    Bm = hd.normalise((tripod(seed=3) * hd.MIRROR) @ R.T)
    dp, _ = hd.best_fit(A, Bm, False, np.random.default_rng(5), n_rot=300)
    di, _ = hd.best_fit(A, Bm, True, np.random.default_rng(5), n_rot=300)
    assert di < dp


def test_handedness_sign_convention():
    # x = left, y = cranial, z = ventral is right-handed; a situs-solitus embryo then gives h = -1
    assert hd.handed_index(np.array([0.0, 1, 0]), np.array([0.0, 0, -1]), np.array([1.0, 0, 0])) == pytest.approx(-1)
    assert hd.handed_index(np.array([0.0, 1, 0]), np.array([0.0, 0, -1]), np.array([-1.0, 0, 0])) == pytest.approx(1)


def _stage(labels: dict[str, list[str]]):
    n = len(next(iter(labels.values())))
    return hd.StageData(name="s", C=np.zeros((n, 3)), core=np.zeros(n, bool), groups=np.array(["CM"] * n),
                        labels={k: np.asarray(v) for k, v in labels.items()})


def test_landmark_rule_picks_tiers_in_order():
    cells = ["OFT-CM"] * 60 + ["SV-CM"] * 60 + ["A-CM"] * 60 + ["Neural Tube"] * 60 + ["Forebrain"] * 60 \
        + ["Hindgut"] * 60 + ["Hepatocytes"] * 60
    lm = hd.rule_landmarks(_stage({"celltype": cells}))
    assert lm["heart_arterial"] == ("celltype", ["OFT-CM"])
    assert lm["heart_venous"] == ("celltype", ["SV-CM"])  # the SV tier wins over atrial labels
    assert lm["body_anterior"] == ("celltype", ["Forebrain"])
    assert lm["body_posterior"] == ("celltype", ["Hindgut"])
    assert lm["dorsal"] == ("celltype", ["Neural Tube"]) and lm["ventral"] == ("core", [True])
    early = ["D-FG"] * 60 + ["V-FG"] * 60 + ["Forebrain"] * 10
    lm = hd.rule_landmarks(_stage({"celltype": early, "cm_celltype": ["aSHF1"] * 70 + ["pSHF"] * 60}))
    assert lm["dorsal"] == ("celltype", ["D-FG"]) and lm["ventral"] == ("celltype", ["V-FG"])
    assert lm["heart_arterial"] == ("cm_celltype", ["aSHF1"]) and lm["heart_venous"] == ("cm_celltype", ["pSHF"])
    assert "body_anterior" not in lm  # 10 Forebrain cells are below the minimum


def test_rule_verdict_strict_and_secondary():
    rob = pd.DataFrame({"gene": ["Pitx2"] * 4, "h": [-0.6, -0.5, -0.4, 0.05]})
    assert hd.rule_verdict(rob)[0].startswith("unresolved")  # strict: one near-zero positive variant breaks it
    assert hd.rule_verdict(rob, hd.H_INFORMATIVE)[0].startswith("resolved: h < 0")  # secondary: it carries no sign
    rob = pd.DataFrame({"gene": ["Pitx2"] * 4, "h": [-0.6, 0.3, 0.05, -0.02]})
    assert hd.rule_verdict(rob, hd.H_INFORMATIVE)[0].startswith("unresolved")
    rob = pd.DataFrame({"gene": ["Pitx2"] * 4, "h": [-0.6, 0.05, 0.02, -0.01]})
    assert hd.rule_verdict(rob, hd.H_INFORMATIVE)[0].startswith("not enough informative")
    assert hd.combine_verdicts(["resolved: h < 0 x", "not testable y"]) == \
        "resolved: h < 0 (right-handed frame, situs solitus)"
    assert hd.combine_verdicts(["resolved: h < 0 x", "unresolved: y"]) == "unresolved"


def test_alignment_pairs_count_only_between_chiral_stages():
    def rows(a, b, kind, dp, di):
        return [{"A": a, "B": b, "kind": kind, "seed": k, "nn_dist_proper": dp, "nn_dist_improper": di,
                 "rel_margin": di / dp - 1, "dice_proper": 0.9, "dice_improper": 0.8} for k in range(3)]

    df = pd.DataFrame(rows("X", "X", "same_stage_halves", 0.07, 0.09) + rows("Y", "Y", "same_stage_halves", 0.07, 0.07)
                      + rows("Z", "Z", "same_stage_halves", 0.07, 0.08) + rows("X", "Y", "pair", 0.10, 0.09)
                      + rows("X", "Z", "pair", 0.10, 0.12) + rows("X", "W", "pair", 0.10, 0.12))
    out = hd.summarise_alignment(df).set_index(["A", "B"])
    assert out.loc[("X", "X"), "usable_verdict"] == "chiral"
    assert out.loc[("Y", "Y"), "usable_verdict"].startswith("not chiral")
    assert out.loc[("X", "Y"), "verdict"] == "improper" and out.loc[("X", "Y"), "usable_verdict"].startswith(
        "uninformative: Y")
    assert out.loc[("X", "Z"), "usable_verdict"] == "proper"
    assert out.loc[("X", "W"), "usable_verdict"].startswith("uninformative: no split-half calibration for W")
