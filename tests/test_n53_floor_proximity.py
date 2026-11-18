"""experiments/harness/n53_floor_proximity.py (row N53 of experiments/harness/TRIES.md), on synthetic stages only: the
arms are the harness's own builders, the first upload's geometry is never built, the noise floor is a copy of E7.5 at
the arms' count independent of the reference draw, the rule reads as its row says, and nothing written or printed is a
cell count."""

from __future__ import annotations

import dataclasses
import importlib.util
import json
import re
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from vec import harness
from vec.boards import Board
from vec.io import COORD_KEY

REPO = Path(__file__).resolve().parents[1]
GENES = [f"g{i}" for i in range(18)]
SIZES = {"E7.5": 413, "E8.0": 527}  # the synthetic stages' cells; never printed or written
LIMITS = {"min_cells": 100, "max_cells": 300}  # stand for the validation board's limits; never printed or written


@pytest.fixture(scope="module")
def n53():
    spec = importlib.util.spec_from_file_location("harness_n53_floor_proximity",
                                                  REPO / "experiments" / "harness" / "n53_floor_proximity.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def _stage(name: str, seed: int) -> ad.AnnData:
    """A bracket as load_stage returns it: E7.5 a round cloud, E8.0 a long flat one, expression a function of position
    so the coupling has something to pair."""
    n = SIZES[name]
    rng = np.random.default_rng(seed)
    C = rng.normal(size=(n, 3)) * ([1.0, 1.0, 0.8] if name == "E7.5" else [3.0, 1.0, 0.3]) * 40.0
    centres = np.random.default_rng(0).normal(size=(len(GENES), 3)) * 40.0
    E = np.exp(-((C[:, None, :] / [1.0, 1.0, 1.0] - centres[None]) ** 2).sum(-1) / 4000.0)
    counts = rng.poisson(30.0 * E + 0.3)
    lib = counts.sum(1, keepdims=True).clip(min=1)
    a = ad.AnnData(X=np.log1p(1e4 * counts / lib).astype(np.float32),
                   obs=pd.DataFrame({"celltype": rng.choice(["a", "b"], n)}, index=[f"c{i}" for i in range(n)]),
                   var=pd.DataFrame(index=GENES))
    a.obsm[COORD_KEY] = C
    return a


def _board(n53):
    return Board(key=n53.BOARD, task="T2", genes=GENES, needs_coords=True, setting="embryo", split="val_interp",
                 **LIMITS)


def _inputs():
    return {"before": _stage("E7.5", 1), "after": _stage("E8.0", 2)}


def test_the_reference_draws_are_disjoint_and_seeded(n53):
    a, b = n53.reference_draws(413, 1_000_001)
    assert len(a) == len(b) == round(0.1 * 413) and not set(a) & set(b)
    a2, b2 = n53.reference_draws(413, 1_000_001)
    assert np.array_equal(a, a2) and np.array_equal(b, b2)
    assert not np.array_equal(n53.reference_draws(413, 1_000_002)[0], a)
    with pytest.raises(ValueError):
        n53.reference_draws(4, 1)
    assert n53.SEEDS == tuple(range(1_000_001, 1_000_011)) and n53.BOARD == "T2:embryo:val_interp"


def test_the_noise_draw_is_a_copy_at_the_arms_count_independent_of_draw_a(n53):
    """Draw B has the arms' count, no repeats while the pool suffices, its own stream, and shares about REF_FRAC of
    its cells with draw A, as a uniform sample of the pool does (fgw_shape0x's overlap without the fallback)."""
    b = n53.copy_draw(413, 300, 1_000_001)
    assert len(b) == 300 and len(np.unique(b)) == 300
    assert np.array_equal(b, n53.copy_draw(413, 300, 1_000_001))
    assert not np.array_equal(b, n53.copy_draw(413, 300, 1_000_002))
    assert len(n53.copy_draw(413, 500, 1)) == 500  # beyond the pool: every cell once, the rest repeat
    share = [np.isin(n53.copy_draw(4000, 3000, s), n53.reference_draws(4000, s)[0]).mean() for s in range(1, 41)]
    assert np.mean(share) == pytest.approx(n53.REF_FRAC, abs=0.01)
    assert n53.COPY_STREAM != n53.REF_STREAM


def test_the_first_uploads_geometry_is_not_an_arm(n53):
    assert n53.NOT_AN_ARM == "fgw" and "fgw" not in n53.ARMS
    assert n53.ARMS == ("fgw_shape0", "fgw_shape0x", "fgw_aligned")
    ctx = n53.context(_board(n53), _inputs(), 1_000_001)
    for arm in ("fgw", "fgw_series", "copy"):
        with pytest.raises(ValueError, match="not an arm"):
            n53.positions(ctx, arm)


def test_the_arms_are_the_harness_builders_and_can_differ(n53):
    board, inputs = _board(n53), _inputs()
    ctx = n53.context(board, inputs, 1_000_003)
    assert ctx.board.truth_cells is None and ctx.board.ref_cells is None  # as a candidate sees a board
    got = {arm: n53.positions(ctx, arm) for arm in n53.ARMS}
    for arm in n53.ARMS:
        np.testing.assert_array_equal(got[arm], harness.CANDIDATES[arm].build(ctx).coords)
    assert not np.allclose(got["fgw_shape0x"], got["fgw_shape0"])  # the control can differ from the new arm
    assert not np.allclose(got["fgw_aligned"], got["fgw_shape0"])


def test_the_shape_pair_is_the_scorers_values(n53):
    rng = np.random.default_rng(3)
    A = (rng.normal(size=(600, 3)) * [3.0, 2.0, 1.0]).astype(np.float32)
    same = n53.shape_pair(A, A)
    # the scorer samples pairs of cells for d2_shape, so a cloud against itself is near 0, not 0
    assert same["occupancy_dice"] == pytest.approx(1.0) and 0.0 <= same["d2_shape"] < 0.01
    flat = n53.shape_pair((A * [1.0, 1.0, 0.2]).astype(np.float32), A)
    assert flat["occupancy_dice"] < same["occupancy_dice"] and flat["d2_shape"] > same["d2_shape"]


def test_the_outside_draw_is_a_copy_at_the_arms_count_that_shares_nothing_with_draw_a(n53):
    """Draw B0 has the arms' count, no repeats, no cell of draw A, its own stream; beyond the cells outside A it is
    draw B and says so."""
    a, _ = n53.reference_draws(413, 1_000_001)
    b0, fell_back = n53.outside_draw(413, 300, a, 1_000_001)
    assert not fell_back and len(b0) == 300 and len(np.unique(b0)) == 300 and not np.isin(b0, a).any()
    assert np.array_equal(b0, n53.outside_draw(413, 300, a, 1_000_001)[0])
    assert not np.array_equal(b0, n53.outside_draw(413, 300, a, 1_000_002)[0])
    assert n53.OUTSIDE_STREAM not in (n53.REF_STREAM, n53.COPY_STREAM)
    over = 413 - len(a) + 1
    b0, fell_back = n53.outside_draw(413, over, a, 1_000_001)
    assert fell_back and np.array_equal(b0, n53.copy_draw(413, over, 1_000_001))
    assert not n53.outside_draw(413, over - 1, a, 1_000_001)[1]


def test_each_arm_is_read_against_the_floor_matched_to_its_overlap(n53):
    assert n53.FLOORS == {"fgw_shape0x": ("noise",), "fgw_shape0": ("noise", "noise_outside"),
                          "fgw_aligned": ("noise_outside",)}


def _rows(n53, rng, shift=None, overlap=0.0):
    """Ten seeds of synthetic values. `overlap` stands for the closeness that shared cells give: B and fgw_shape0x
    (about 10 percent shared) sit that much closer to A than B0 and fgw_aligned (none), fgw_shape0 half as much."""
    near = {"noise": 1.0, "noise_outside": 0.0, "noise_disjoint": 0.0,
            "fgw_shape0x": 1.0, "fgw_shape0": 0.5, "fgw_aligned": 0.0}
    rows = []
    for s in n53.SEEDS:
        r = {"seed": s}
        for key in ("noise", "noise_outside", "noise_disjoint", *n53.ARMS):
            d = (shift or {}).get(key, {})
            r[key] = {"occupancy_dice": rng.normal(0.6, 0.02) + overlap * near[key] + d.get("occupancy_dice", 0.0),
                      "d2_shape": rng.normal(0.05, 0.01) - 0.5 * overlap * near[key] + d.get("d2_shape", 0.0)}
        rows.append(r)
    return rows


def test_the_matched_floors_on_pure_noise_with_an_overlap_bias(n53):
    """Arms that are copies of E7.5 apart from their overlap with A: against its own floor each is called in form about
    85 percent of the time whatever the bias. Read against B alone, a large bias would push fgw_aligned and fgw_shape0
    out; the rule does not, and the comparison against B is reported beside it."""
    rng = np.random.default_rng(7)
    calls = {arm: [] for arm in n53.ARMS}
    beside_b = {"fgw_aligned": [], "fgw_shape0": []}
    for _ in range(300):
        s = n53.summary(_rows(n53, rng, overlap=0.04))
        for arm in n53.ARMS:
            calls[arm].append(s["in_form"][arm]["in_form"])
        for arm in beside_b:
            beside_b[arm].append(s["in_form"][arm]["beside"]["B"]["in_form"])
    for arm in n53.ARMS:
        assert 0.78 < np.mean(calls[arm]) < 0.92, arm
    for arm in beside_b:
        assert np.mean(beside_b[arm]) < 0.5, arm
    r = n53.summary(_rows(n53, rng))["in_form"]
    assert r["fgw_shape0x"]["floor"] == "B" and set(r["fgw_shape0x"]["beside"]) == {"B0", "mean of B and B0"}
    assert r["fgw_aligned"]["floor"] == "B0" and r["fgw_shape0"]["floor"] == "mean of B and B0"


def test_the_rule_on_pure_noise_and_on_shifted_arms(n53):
    """Pure noise: arms drawn like the noise floor are called in the floor's form unless the paired t of ten seeds
    passes 2 on a metric (about 7.6 percent per metric, so about 85 percent in form on both). An arm whose seeds sit
    1.4 single-draw noise SDs off on one metric, each seed beside its noise draw, is out (the first version's 2-SD bar
    would have called it in); one 4 SDs off is out on both bars."""
    rng = np.random.default_rng(4)
    inside = [[v["in_form"] for v in n53.summary(_rows(n53, rng))["in_form"].values()] for _ in range(300)]
    assert 0.78 < np.mean(inside) < 0.92
    sign = np.array([1.0, -1.0] * 5)
    rows = _rows(n53, rng)
    for r, x in zip(rows, sign, strict=True):
        r["noise"] = {"occupancy_dice": 0.6 + 0.02 * x, "d2_shape": 0.05 + 0.01 * x}
        r["noise_outside"] = dict(r["noise"])
        r["fgw_shape0"] = {"occupancy_dice": 0.6 + 0.02 * x - 0.03 + 0.001 * x * x, "d2_shape": 0.05 + 0.01 * x}
        r["fgw_aligned"] = {"occupancy_dice": 0.6 + 0.02 * x, "d2_shape": 0.05 + 0.01 * x + 0.042 + 0.001 * x}
    s = n53.summary(rows)
    r = s["in_form"]["fgw_shape0"]
    assert not r["in_form"] and not r["occupancy_dice"]["inside"] and r["occupancy_dice"]["inside_wide_bar"]
    assert r["d2_shape"]["inside"]
    a = s["in_form"]["fgw_aligned"]
    assert not a["in_form"] and not a["d2_shape"]["inside"] and not a["d2_shape"]["inside_wide_bar"]
    assert set(s) == {"noise", "noise_outside", "noise_disjoint", "arms", "in_form"}
    # the boundary: a paired mean just inside two standard errors is inside, just beyond is out
    d = np.array([0.01, -0.01] * 5)
    edge = 2.0 * d.std(ddof=1) / np.sqrt(10)

    def at(f):
        return [{"noise": {"occupancy_dice": 0.5, "d2_shape": 0.1},
                 "noise_outside": {"occupancy_dice": 0.4, "d2_shape": 0.2},
                 "fgw_shape0x": {"occupancy_dice": 0.5 + f * edge + x, "d2_shape": 0.1 - f * edge + x}} for x in d]

    assert n53.in_form(at(0.999), "fgw_shape0x")["in_form"]
    out = n53.in_form(at(1.001), "fgw_shape0x")
    assert not out["in_form"] and not out["occupancy_dice"]["inside"] and not out["d2_shape"]["inside"]


def _numbers(text: str) -> set[int]:
    return {int(x) for x in re.findall(r"(?<![\d.])\d+(?![\d.])", text.replace(",", ""))}


def test_a_run_writes_aggregates_only_and_never_builds_fgw(n53, monkeypatch, tmp_path, capsys):
    board = _board(n53)
    stages = _inputs()
    monkeypatch.setattr(n53, "load_boards", lambda: {n53.BOARD: board})
    monkeypatch.setattr(n53, "load_stage", lambda path, genes=None, drop_duplicates=False:
                        stages["before" if "E7.5" in str(path) else "after"][:, genes].copy())
    monkeypatch.setattr(n53.n12, "tries_rows", lambda path: {"N53"})
    monkeypatch.setattr(n53, "_state", lambda: ("f" * 40, True, "2026-11-18T00:00:00Z"))
    monkeypatch.setattr(n53, "SEEDS", n53.SEEDS[:3])

    def refuse(ctx):
        raise AssertionError("the first upload's geometry was built")

    monkeypatch.setitem(harness.CANDIDATES, "fgw", dataclasses.replace(harness.CANDIDATES["fgw"], build=refuse))
    out = tmp_path / "n53"
    n53.main(["run", "--out", str(out), "--data", str(tmp_path)])
    rep = json.loads((out / "n53.json").read_text())
    assert set(rep["arms"]) == set(n53.ARMS) and set(rep["in_form"]) == set(n53.ARMS)
    assert len(rep["rows"]) == 3 and all(set(r) == {"seed", "noise", "noise_outside", "noise_disjoint", *n53.ARMS}
                                         for r in rep["rows"])
    assert rep["floor_exact_fallback"] is False  # a boolean only (300 of 413 cells: the unused cells suffice)
    assert rep["outside_fallback"] is False  # a boolean only (300 of the cells outside draw A suffice)
    assert rep["floors"] == {"fgw_shape0": "mean of B and B0", "fgw_shape0x": "B", "fgw_aligned": "B0"}
    secret = {*SIZES.values(), *LIMITS.values(), round(0.1 * SIZES["E7.5"]), LIMITS["max_cells"] // 2}
    for written in ((out / "n53.json").read_text(), (out / "n53.txt").read_text(), capsys.readouterr().out):
        assert not _numbers(written) & secret
    with pytest.raises(SystemExit, match="refusing to overwrite"):
        n53.main(["run", "--out", str(out), "--data", str(tmp_path)])


def test_a_run_refuses_without_its_row_or_on_a_dirty_tree(n53, monkeypatch, tmp_path):
    monkeypatch.setattr(n53.n12, "tries_rows", lambda path: {"N51"})
    with pytest.raises(SystemExit, match="N53 has no row"):
        n53.main(["run", "--out", str(tmp_path / "a")])
    monkeypatch.setattr(n53.n12, "tries_rows", lambda path: {"N53"})
    monkeypatch.setattr(n53, "_state", lambda: ("f" * 40, False, "now"))
    with pytest.raises(SystemExit, match="not clean"):
        n53.main(["run", "--out", str(tmp_path / "b")])
    assert not (tmp_path / "a").exists() and not (tmp_path / "b").exists()
