"""experiments/harness/diag_emb_midpoint.py (N16 of experiments/harness/TRIES.md), on synthetic boards and stages only:
each pair reads its own panel and only its bracket stages, nothing written or printed is a cell count, the MMD effect
size separates a real tilt from a redraw, and the gate's rule holds at its boundaries."""

from __future__ import annotations

import importlib.util
import json
import re
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse
from scipy.spatial.distance import cdist

from vec.boards import Board
from vec.harness import BUILD_SEED_OFFSET
from vec.io import COORD_KEY

REPO = Path(__file__).resolve().parents[1]
ALL_GENES = [f"g{i}" for i in range(24)]
VAL_GENES = ALL_GENES[:20]
TEST_GENES = ALL_GENES[3:22]  # another panel, as the test board's is
SECRET = {"min_cells": 777, "max_cells": 4321, "truth_cells": 2468, "ref_cells": 1357}  # stand for a test board's
SIZES = {"E6.75": 229, "E7.25": 251, "E7.5": 263, "E8.0": 241}  # the synthetic stages' cells
COUNT_KEYS = {"cells", "n", "n_obs", "n_before", "n_after", "n_jittered", "n_cells", "n_anchor", "coupling_nonzero",
              *SECRET}
N_DRAW = 160


@pytest.fixture(scope="module")
def dm():
    spec = importlib.util.spec_from_file_location("harness_diag_emb_midpoint",
                                                  REPO / "experiments" / "harness" / "diag_emb_midpoint.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def _stage(name: str) -> ad.AnnData:
    """A stage as load_stage returns it (counts, X = log(1 + 1e4 c / N), three spatial coordinates), with four states
    placed in space so that a coupling can pair them."""
    n = SIZES[name]
    rng = np.random.default_rng(n)
    profiles = np.random.default_rng(1).gamma(2.0, 1.5, size=(4, len(ALL_GENES)))
    states = rng.integers(0, 4, n)
    depth = 2.0 + sorted(SIZES).index(name)  # later stages deeper
    C = rng.poisson(profiles[states] * depth)
    C[:, 3] += 1  # a count on a gene both panels hold: every cell has a panel gene
    N = C.sum(1, keepdims=True)
    a = ad.AnnData(X=sparse.csr_matrix(np.log1p(1e4 * C / N).astype(np.float32)),
                   obs=pd.DataFrame(index=[f"c{i}" for i in range(n)]), var=pd.DataFrame(index=ALL_GENES))
    a.layers["counts"] = sparse.csr_matrix(C.astype(np.int32))
    blocks = np.array([[-4.0, 0, 0], [4.0, 0, 0], [0, -4.0, 0], [0, 4.0, 0]])
    a.obsm[COORD_KEY] = blocks[states] + rng.normal(size=(n, 3))
    return a


def _patch(dm, monkeypatch, seeds=(0, 1)) -> list[tuple[str, list[str], bool]]:
    """Synthetic boards and stages in place of boards/ and data/; returns the load_stage calls."""
    boards = {"T2:embryo:val_interp": Board(key="T2:embryo:val_interp", task="T2", genes=VAL_GENES, needs_coords=True,
                                            split="val_interp"),
              "T2:embryo:test_interp": Board(key="T2:embryo:test_interp", task="T2", genes=TEST_GENES,
                                             needs_coords=True, split="test_interp", **SECRET)}
    calls = []

    def load_stage(path, genes=None, drop_duplicates=False):
        name = Path(path).name.removesuffix(".h5ad")
        calls.append((name, list(genes), drop_duplicates))
        return _stage(name)[:, genes].copy()

    monkeypatch.setattr(dm, "load_boards", lambda: boards)
    monkeypatch.setattr(dm, "load_stage", load_stage)
    monkeypatch.setattr(dm, "sha256_file", lambda path: "ab" * 32)
    monkeypatch.setattr(dm.provenance, "git_state", lambda repo, ignore=(): ("cafebabe", []))
    monkeypatch.setattr(dm, "SEEDS", seeds)
    monkeypatch.setattr(dm, "N_DRAW", N_DRAW)
    return calls


def _leaves(x, keys: set, ints: set, texts: list) -> None:
    if isinstance(x, dict):
        keys.update(x)
        for v in x.values():
            _leaves(v, keys, ints, texts)
    elif isinstance(x, list):
        for v in x:
            _leaves(v, keys, ints, texts)
    elif isinstance(x, int) and not isinstance(x, bool):
        ints.add(x)
    elif isinstance(x, str):
        texts.append(x)


def _whole_numbers(text: str) -> set[int]:
    """Every whole number in a text (thousands separators dropped), decimals and the stamp's UTC time left out."""
    text = re.sub(r"\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z", "", text)
    return {int(d.replace(",", "")) for d in re.findall(r"(?<![\d.,])\d+(?:,\d{3})*(?![.\d])", text)}


def test_the_gate_reads_each_pairs_panel_and_writes_ratios_only(dm, tmp_path, monkeypatch, capsys):
    calls = _patch(dm, monkeypatch)
    out = tmp_path / "gate"
    dm.main(["--out", str(out)])
    printed = capsys.readouterr().out
    # every stage comes from data/, deduplicated, each on its own pair's panel
    by_stage = {}
    for name, genes, dedup in calls:
        assert dedup
        by_stage.setdefault(name, []).append(genes)
    assert sorted(by_stage) == ["E6.75", "E7.25", "E7.5", "E8.0"]
    assert by_stage["E6.75"] == [VAL_GENES] and by_stage["E7.25"] == [VAL_GENES]
    assert TEST_GENES in by_stage["E7.5"] and VAL_GENES in by_stage["E7.5"]  # E7.5: an input of a backtest and the test
    rep = json.loads((out / "diag_emb_midpoint.json").read_text())
    txt = (out / "diag_emb_midpoint.txt").read_text()
    assert list(rep["pairs"]) == [k for k, *_ in dm.PAIRS] and rep["seeds"] == [0, 1]
    assert rep["verdict"] in ("pass", "fail") and printed.rstrip().endswith(f"verdict: {rep['verdict']}")
    assert rep["stamp"]["commit"] == "cafebabe" and rep["stamp"]["tree_clean"] is True and rep["stamp"]["utc"]
    assert set(rep["inputs"]["stages"]) == set(SIZES) and all("sha256" in v for v in rep["inputs"]["stages"].values())
    assert set(rep["inputs"]["panels"]) == {"T2:embryo:val_interp", "T2:embryo:test_interp"}
    test_pair = rep["pairs"]["test/E7.5+E8.0->E7.75"]
    assert test_pair["panel_board"] == "T2:embryo:test_interp" and test_pair["t"] == pytest.approx(0.5)
    assert [d["build_seed"] for d in test_pair["draws"]] == [BUILD_SEED_OFFSET, BUILD_SEED_OFFSET + 1]
    for p in rep["pairs"].values():
        for d in p["draws"]:
            assert 0.0 < d["ess_before"] <= 1.0 and 0.0 < d["ess_after"] <= 1.0
            assert d["mmd_uniform_uniform"] > 0.0 and d["mmd_ratio"] == pytest.approx(
                d["mmd_tilted_uniform"] / d["mmd_uniform_uniform"])
            assert d["library_max_shift"] == pytest.approx(max(abs(v - 1.0) for v in d["library_ratio"].values()))
            assert d["strata_exact"] is True  # the quotas hold the pool's stratum shares
            # the tilt moves each bracket's mean signed offset towards the midpoint (zero)
            for side in ("before", "after"):
                assert abs(d["s_mean"][side]["tilted"]) < abs(d["s_mean"][side]["pool"])
    # nothing written or printed is a cell count: no count field; the whole numbers are the seeds, the rule's and
    # streams' small numbers, and a date and times in the stamp
    keys, ints, texts = set(), set(), []
    _leaves(rep, keys, ints, texts)
    assert not keys & COUNT_KEYS
    assert ints <= {*range(13), *(BUILD_SEED_OFFSET + s for s in range(10))}
    forbidden = {*SECRET.values(), *SIZES.values(), N_DRAW, N_DRAW // 2, N_DRAW - N_DRAW // 3, N_DRAW // 3,
                 round(2 / 3 * N_DRAW), round(N_DRAW / 3), len(VAL_GENES), len(TEST_GENES)}
    for text in (*texts, txt, printed):
        assert not _whole_numbers(text) & forbidden, text


def test_a_partial_run_gives_no_verdict_and_outputs_are_never_overwritten(dm, tmp_path, monkeypatch, capsys):
    calls = _patch(dm, monkeypatch, seeds=tuple(range(10)))
    out = tmp_path / "gate"
    dm.main(["--out", str(out), "--seeds", "0"])
    rep = json.loads((out / "diag_emb_midpoint.json").read_text())
    assert rep["verdict"] == "partial run: no verdict" and rep["seed_span"] == f"{BUILD_SEED_OFFSET:,}"
    assert capsys.readouterr().out.rstrip().endswith("verdict: partial run: no verdict")
    n_calls = len(calls)
    with pytest.raises(SystemExit, match="already written"):
        dm.main(["--out", str(out), "--seeds", "0"])
    assert len(calls) == n_calls  # refused before reading any data
    with pytest.raises(SystemExit, match="at least one seed"):
        dm.main(["--out", str(tmp_path / "other"), "--seeds"])


def test_a_missing_board_stops_the_gate(dm, tmp_path, monkeypatch):
    _patch(dm, monkeypatch)
    monkeypatch.setattr(dm, "load_boards", lambda: {})
    with pytest.raises(SystemExit, match="not in boards"):
        dm.main(["--out", str(tmp_path / "gate")])


def _draw(ess=0.8, shift=0.01, ratio=4.0, flags=0, exact=True):
    return {"ess_before": ess, "ess_after": ess, "library_max_shift": shift, "mmd_ratio": ratio, "guard_flags": flags,
            "strata_exact": exact}


@pytest.mark.parametrize("change, passes", [
    ({}, True),
    ({"ess": 0.5}, True),  # exactly the floor
    ({"ess": 0.4999}, False),
    ({"shift": 0.05}, True),  # exactly the tolerance
    ({"shift": 0.0501}, False),
    ({"flags": 1}, False),  # any cloud-guard flag
    ({"exact": False}, False),  # stratum drift
])
def test_the_pair_rule_holds_at_its_boundaries(dm, change, passes):
    draws = [_draw() for _ in range(10)]
    draws[3] = _draw(**change)
    assert dm.pair_verdict(draws)["passes"] is passes


@pytest.mark.parametrize("at_floor, passes", [(10, True), (8, True), (7, False)])
def test_the_effect_size_needs_eight_of_ten_seeds_at_three_times_the_null(dm, at_floor, passes):
    draws = [_draw(ratio=3.0 if i < at_floor else 2.99) for i in range(10)]
    v = dm.pair_verdict(draws)
    assert v["passes"] is passes and v["mmd_seeds_at_floor"] == at_floor and v["mmd_ok"] is passes


def test_the_gate_passes_only_if_every_pair_passes(dm):
    ok = {"verdict": {"passes": True}}
    bad = {"verdict": {"passes": False}}
    assert dm.gate_verdict({"a": ok, "b": ok, "c": ok}, True) == "pass"
    assert dm.gate_verdict({"a": ok, "b": bad, "c": ok}, True) == "fail"
    assert dm.gate_verdict({"a": ok, "b": ok, "c": ok}, False) == "partial run: no verdict"


def test_effect_size_separates_a_tilted_draw_from_a_redraw_of_the_same_population(dm):
    rng = np.random.default_rng(50)
    pool = rng.normal(size=(20_000, 6))
    uni, uni2 = pool[rng.choice(20_000, 2000, replace=False)], pool[rng.choice(20_000, 2000, replace=False)]
    shifted = dm.effect_size(uni[:1000] + np.array([0.35, 0, 0, 0, 0, 0]), uni, uni2, seed=3)
    assert shifted["ratio"] > 3.0 and shifted["d_tilt"] > shifted["d_null"] > 0.0 and "ratio_flat" not in shifted
    redraw = pool[rng.choice(20_000, 2000, replace=False)]
    same = dm.effect_size(redraw, uni, uni2, seed=3, flat=redraw)
    assert same["ratio"] < 1.5  # a no-op looks like the two uniform draws
    assert same["ratio_flat"] == pytest.approx(same["ratio"])  # the no-op reference is computed like the numerator


def test_the_no_op_reference_sits_below_one_when_the_comparison_is_stratified(dm):
    """The reason the gate reports it: a stratified draw against a plain uniform draw has less composition noise than
    two plain uniform draws, so a no-op's MMD ratio is below 1 (the level N16's ratio should be read against)."""
    rng = np.random.default_rng(52)
    states = rng.integers(0, 8, 40_000)
    centres = 6.0 * rng.normal(size=(8, 5))
    pool = centres[states] + rng.normal(size=(40_000, 5))
    n = 4000
    uni = pool[rng.choice(40_000, n, replace=False)]
    uni2 = pool[rng.choice(40_000, n, replace=False)]
    parts = []
    for st in range(8):
        members = np.flatnonzero(states == st)
        parts.append(rng.choice(members, int(round(n * len(members) / 40_000)), replace=False))
    flat = pool[np.concatenate(parts)[:n]]
    eff = dm.effect_size(flat, uni, uni2, seed=5, flat=flat)
    assert eff["ratio_flat"] == pytest.approx(eff["ratio"]) and 0.4 < eff["ratio_flat"] < 1.0


def test_kernel_means_and_median_helpers(dm):
    rng = np.random.default_rng(51)
    A, B = rng.normal(size=(30, 4)), rng.normal(size=(50, 4))
    gamma = 0.7
    assert dm.kernel_mean(A, B, gamma) == pytest.approx(float(np.exp(-gamma * cdist(A, B, "sqeuclidean")).mean()))
    assert dm.mmd_v(0.5, 0.5, 0.5) == 0.0 and dm.mmd_v(1.0, 1.0, 0.5) == pytest.approx(1.0)
    assert dm.weighted_median(np.array([3.0, 1.0, 2.0]), np.ones(3)) == 2.0
    assert dm.weighted_median(np.array([1.0, 2.0, 3.0]), np.array([0.9, 0.05, 0.05])) == 1.0
    assert dm.weighted_median(np.array([1.0, 2.0]), np.array([0.0, 1.0])) == 2.0
    with pytest.raises(ValueError, match="all zero"):
        dm.weighted_median(np.array([1.0, 2.0]), np.zeros(2))
    assert dm.rbf_gamma(np.zeros((3, 2)), 1) > 0.0  # a degenerate sample still gives a finite rate
