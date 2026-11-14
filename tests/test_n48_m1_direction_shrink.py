"""Row N48 (experiments/t1_m1/n48_m1_direction_shrink.py): the target-free reliability shrinkage of M1's real
direction, on synthetic directions built here. No statistics file, no released file and no score is read."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments" / "t1_m1"))
sys.path.insert(0, str(ROOT / "experiments" / "t1_hygiene"))

import n48_m1_direction_shrink as n48  # noqa: E402


def halves(signal: np.ndarray, noise_sd: float, seed: int = 0) -> tuple[np.ndarray, dict]:
    """A direction and its two half-sample versions: each half is the signal plus independent noise of variance
    2 noise_sd^2, so their mean d has noise variance noise_sd^2."""
    rng = np.random.default_rng(seed)
    d_a = signal + rng.normal(0.0, noise_sd * np.sqrt(2.0), signal.size)
    d_b = signal + rng.normal(0.0, noise_sd * np.sqrt(2.0), signal.size)
    basis = np.linspace(0.1, 50.0, signal.size)
    return (d_a + d_b) / 2.0, {"d_a": d_a, "d_b": d_b, "mean_cp": basis}


def test_the_constants_are_t18s_stop():
    assert n48.ROW == "N48" and n48.NO_OP == 0.9 and n48.EXPRESSED_FROM == 3
    assert n48.STOP_GROUPS == ("CM", "Endocardium") and n48.REAL == ("E9.5", "E14.5")


def test_the_scalar_scale_is_signal_over_signal_plus_noise():
    rng = np.random.default_rng(1)
    signal = rng.normal(0.0, 1.0, 20000)
    d, h = halves(signal, 0.5)
    js = n48.js_scale(d, h["d_a"], h["d_b"])
    assert js["genes"] == 20000
    assert js["noise"] == pytest.approx(0.25, rel=0.05), "the mean d of two halves carries half a half's variance"
    assert js["scale"] == pytest.approx(1.0 / 1.25, abs=0.02)
    d, h = halves(signal, 0.01)
    assert n48.js_scale(d, h["d_a"], h["d_b"])["scale"] > 0.999
    d, h = halves(np.zeros(5000), 1.0)
    assert n48.js_scale(d, h["d_a"], h["d_b"])["scale"] < 0.1, "pure noise is shrunk to almost nothing"


def test_the_scalar_scale_by_hand_and_at_its_edges():
    d, d_a, d_b = np.array([1.0, -1.0, np.nan]), np.array([1.5, -0.5, 1.0]), np.array([0.5, -1.5, 1.0])
    js = n48.js_scale(d, d_a, d_b)
    assert js == {"genes": 2, "noise": 0.25, "mean_sq": 1.0, "scale": 0.75}
    assert n48.js_scale([np.nan], [1.0], [1.0]) == {"genes": 0, "noise": None, "mean_sq": None, "scale": None}
    assert n48.js_scale([0.0, 0.0], [0.0, 0.0], [0.0, 0.0])["scale"] == 1.0
    assert n48.js_scale([0.1, -0.1], [1.1, -1.1], [-0.9, 0.9])["scale"] == 0.0, "noise above the signal floors at 0"


def test_the_halves_correlation():
    a = np.array([1.0, 2.0, 3.0, np.nan])
    assert n48.half_correlation(a, np.array([2.0, 4.0, 6.0, 1.0])) == pytest.approx(1.0)
    assert n48.half_correlation(a, np.array([3.0, 2.0, 1.0, 1.0])) == pytest.approx(-1.0)
    assert n48.half_correlation(a, np.array([1.0, 1.0, 1.0, 1.0])) is None
    assert n48.half_correlation(a[:2], a[:2]) is None


def test_a_groups_record_holds_the_decile_factors_with_and_without_the_held_genes():
    rng = np.random.default_rng(2)
    signal = rng.normal(0.0, 1.0, 4000)
    d, h = halves(signal, 0.2)
    held = np.zeros(4000, dtype=bool)
    held[:400] = True  # the lowest-expression decile (the basis rises with the index)
    r = n48.group_shrink(d, h, held)
    assert r["genes_with_value"] == 4000 and len(r["eb_factors"]) == 10 and len(r["eb_factors_without_held"]) == 10
    assert all(0.9 < f <= 1.0 for f in r["eb_factors"]), "a clean direction is barely shrunk"
    assert r["js"]["scale"] == pytest.approx(1.0 / 1.04, abs=0.01)
    assert r["js_without_held"]["genes"] == 3600
    assert r["half_correlation"] > 0.9
    import feng_dirs
    _, direct = feng_dirs.eb_shrink(d, h["d_a"], h["d_b"], h["mean_cp"], 10)
    assert r["eb_factors"] == [float(f) for f in direct]


def test_the_stop_needs_both_groups_and_reads_deciles_four_to_ten():
    def rec(factors, scale):
        return {"eb_factors": factors, "js": {"scale": scale}}

    clean = [0.5, 0.6, 0.7, 0.95, 0.96, 0.97, 0.98, 0.99, 0.99, 0.99]  # only deciles 1-3 are shrunk
    noisy = [1.0, 1.0, 1.0, 0.95, 0.96, 0.85, 0.98, 0.99, 0.99, 0.99]  # decile 6 is under the bar
    s = n48.stop({"CM": rec(clean, 0.95), "Endocardium": rec(clean, 0.97), "Epicardium": rec(noisy, 0.5)})
    assert not s["eb_open"] and not s["js_open"], "the other groups do not open a form"
    assert s["eb_below_no_op"] == {"CM": False, "Endocardium": False} and s["deciles"] == "4-10"
    s = n48.stop({"CM": rec(noisy, 0.85), "Endocardium": rec(clean, 0.97)})
    assert not s["eb_open"] and not s["js_open"], "one group alone opens nothing"
    assert s["eb_below_no_op"] == {"CM": True, "Endocardium": False}
    assert s["js_below_no_op"] == {"CM": True, "Endocardium": False}
    s = n48.stop({"CM": rec(noisy, 0.85), "Endocardium": rec(noisy, 0.89)})
    assert s["eb_open"] and s["js_open"]
    s = n48.stop({"CM": rec(noisy, 0.9), "Endocardium": rec(noisy, 0.85)})
    assert s["eb_open"] and not s["js_open"], "a scale of exactly 0.9 is a no-op"
    s = n48.stop({"CM": rec(noisy, 0.5)})
    assert not s["eb_open"] and not s["js_open"] and s["eb_below_no_op"]["Endocardium"] is None
    s = n48.stop({"CM": rec(noisy, None), "Endocardium": rec(noisy, 0.5)})
    assert s["eb_open"] and not s["js_open"] and s["js_below_no_op"]["CM"] is None


def test_the_text_names_each_forms_outcome():
    def rec(factors, scale):
        js = {"genes": 10, "noise": 0.1, "mean_sq": 1.0, "scale": scale}
        return {"genes_with_value": 10, "eb_factors": factors, "js": js, "eb_factors_without_held": factors,
                "js_without_held": js, "half_correlation": 0.8}

    groups = {"CM": rec([0.95] * 10, 0.95), "Endocardium": rec([0.95] * 10, 0.97)}
    report = {"stop": n48.stop(groups), "groups": groups, "real_cells": {"CM": {"E9.5": 100, "E14.5": 200}},
              "stamp": {"commit": "abc", "tree_clean_at_end": True, "utc": "2026-11-14T00:00:00Z"}}
    txt = n48.text(report)
    assert txt.startswith("N48: the reliability shrinkage of M1's real direction")
    assert "decile form (T18's rule): no-op: closed" in txt and "scalar James-Stein scale: no-op: closed" in txt
    assert "  CM: 10; 0.9500 (0.9500); 0.800; " + " ".join(["0.950"] * 10) in txt
    groups = {"CM": rec([0.5] * 10, 0.5), "Endocardium": rec([0.5] * 10, 0.5)}
    txt = n48.text({**report, "stop": n48.stop(groups), "groups": groups})
    assert "decile form (T18's rule): OPEN" in txt and "scalar James-Stein scale: OPEN" in txt
