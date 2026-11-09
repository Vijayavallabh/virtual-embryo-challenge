"""experiments/harness/n17_cfm_gate.py (row N17 of experiments/harness/TRIES.md, the target-free gate of heart try H6),
on synthetic boards and stages only: its pairs, its rule, and that nothing it writes or prints is a cell count."""

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

from vec.boards import Board

REPO = Path(__file__).resolve().parents[1]
BOARD = "T2:heart:test_extrap"
VAL = "T2:heart:val_extrap"
ALL_GENES = [f"g{i}" for i in range(220)]
VAL_GENES, TEST_GENES = ALL_GENES[:200], ALL_GENES[20:]
SIZES = {"CM": 1231, "Endothelial": 877, "Blood": 97, "Untracked": 353}  # the synthetic stages' cells per lineage
FILES = {"E8.5": 0.5, "E8.75": 1.0, "E9.5": 1.5, "E10.5": 2.0}  # stage -> the size factor of its cell numbers
SECRET = {"min_cells": 777, "max_cells": 4321, "truth_cells": 2468, "ref_cells": 1357}  # stand for a test board's


@pytest.fixture(scope="module")
def n17():
    spec = importlib.util.spec_from_file_location("harness_n17_cfm_gate",
                                                  REPO / "experiments" / "harness" / "n17_cfm_gate.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def _stage(stem: str, genes: list[str]) -> ad.AnnData:
    """A stage as load_stage returns it: CM with a coordinated programme, Endothelial without, Blood too small to
    filter, Untracked held flat."""
    seed = round(FILES[stem] * 100)
    rng, blocks, labels = np.random.default_rng(seed), [], []
    pi = np.random.default_rng(1).gamma(0.5, 1.0, len(ALL_GENES))
    pi /= pi.sum()
    a = np.zeros(len(ALL_GENES))
    a[:40] = 2.2 * np.random.default_rng(2).choice([-1.0, 1.0], 40)
    for lab, n0 in SIZES.items():
        n = max(30, int(n0 * FILES[stem]))
        s = rng.normal(size=n) if lab == "CM" else np.zeros(n)
        lam = pi * np.exp(np.outer(s, a) + (lab == "CM") * 0.8 * FILES[stem] * a / 2.2)
        depth = np.maximum(np.round(np.exp(rng.normal(np.log(70.0), 0.4, n))), 5)[:, None]
        C = rng.poisson(lam / lam.sum(1, keepdims=True) * depth).astype(np.int32)
        C[C.sum(1) == 0, 0] = 1
        blocks.append(C)
        labels += [lab] * n
    C = np.vstack(blocks)
    X = sparse.csr_matrix(np.log1p(1e4 * C / C.sum(1, keepdims=True)).astype(np.float32))
    a2 = ad.AnnData(X=X, obs=pd.DataFrame({"celltype": labels}, index=[f"{stem}_{i}" for i in range(len(C))]),
                    var=pd.DataFrame(index=ALL_GENES))
    a2.obs["merged"] = labels
    a2.layers["counts"] = sparse.csr_matrix(C)
    return a2[:, genes].copy()


def _patch(n17, monkeypatch, tmp_path) -> list[tuple]:
    boards = {BOARD: Board(key=BOARD, task="T2", genes=list(TEST_GENES), needs_coords=True, setting="heart",
                           split="test_extrap", **SECRET),
              VAL: Board(key=VAL, task="T2", genes=list(VAL_GENES), needs_coords=True, setting="heart",
                         split="val_extrap")}
    calls = []

    def load_stage(path, genes=None, drop_duplicates=False):
        stem = Path(path).name.removesuffix(".h5ad")
        calls.append((stem, tuple(genes), drop_duplicates))
        return _stage(stem, list(genes))

    monkeypatch.setattr(n17, "load_boards", lambda: boards)
    monkeypatch.setattr(n17, "load_stage", load_stage)
    monkeypatch.setattr(n17, "_file_sha256", lambda p: "0" * 64)
    monkeypatch.setattr(n17, "DATA", tmp_path)
    return calls


def _leaves(x, keys: set, ints: set, texts: list) -> None:
    """Every dict key, whole number (booleans aside) and string in a parsed JSON document."""
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


def test_pairs_are_the_test_builds_and_the_counting_backtests(n17):
    pairs = {p["key"]: p for p in n17.momentum_pairs()}
    assert pairs["test/E9.5->E10.5"] == {"key": "test/E9.5->E10.5", "start": "T2/E9.5.h5ad", "ref": "T2/E10.5.h5ad",
                                         "panel_board": BOARD}
    one, two = pairs["T2h/<=E9.5->E10.5"], pairs["T2h/<=E8.75->E10.5"]
    assert (one["start"], one["ref"], one["panel_board"]) == ("T2/E8.75.h5ad", "T2/E9.5.h5ad", VAL)
    assert (two["start"], two["ref"], two["panel_board"]) == ("T2/E8.5.h5ad", "T2/E8.75.h5ad", VAL)
    assert len(pairs) == 3  # no target stage among any pair's inputs but the test pair's E10.5 (its reference, A20)


def _row(name, eligible=True, k=3, rho=0.6):
    return {"lineage": name, "momentum": "shifted", "eligible": eligible,
            "filter": "filtered" if k else "kept", "directions": k, "retained_energy": 0.4 if k else None,
            "spearman": rho if k else None, "max_abs_ratio": 1.2 if k else None}


@pytest.mark.parametrize("rows, passes, outside", [
    ([_row("A")], True, []),
    ([_row("A", rho=0.30), _row("B", rho=0.95)], True, []),                           # the band is inclusive
    ([_row("A", rho=0.2999)], False, ["A"]),                                          # below: destroyed
    ([_row("A", rho=0.9501)], False, ["A"]),                                          # above: a no-op
    ([_row("A"), _row("B", rho=0.1)], False, ["B"]),                                  # every lineage, not CM alone
    ([_row("A", k=0)], False, []),                                                    # no direction: nothing to test
    ([_row("A"), _row("B", k=0)], True, []),                                          # K = 0 keeps B0's factor
    ([_row("A", eligible=False, k=0), _row("B", eligible=False, k=0)], False, []),    # no lineage eligible
    ([_row("A"), _row("B", eligible=False, k=0, rho=0.0)], True, []),                 # too few cells: not judged
    ([_row("A"), {**_row("B"), "spearman": None}], False, ["B"]),                     # not finite counts as outside
])
def test_the_rule_on_one_pair(n17, rows, passes, outside):
    got = n17.judge_pair(rows)
    assert got["passes"] is passes and got["outside_band"] == outside
    assert got["eligible"] == [r["lineage"] for r in rows if r["eligible"]]
    assert got["acting"] == [r["lineage"] for r in rows if r["eligible"] and r["directions"] >= 1]


def test_the_verdict_needs_all_three_pairs(n17):
    ok, bad = {"judged": {"passes": True}}, {"judged": {"passes": False}}
    assert n17.verdict({"a": ok, "b": ok, "c": ok}).startswith("pass")
    v = n17.verdict({"a": ok, "b": bad, "c": ok})
    assert v.startswith("prune") and "b" in v and "not run and not counted" in v


def test_a_rerun_refuses_to_overwrite_its_outputs_before_reading_anything(n17, tmp_path, monkeypatch):
    (tmp_path / "n17_cfm_gate.json").write_text("{}")

    def no_read(*a, **k):
        raise AssertionError("nothing is read before the refusal")

    monkeypatch.setattr(n17, "load_boards", no_read)
    monkeypatch.setattr(n17, "load_stage", no_read)
    with pytest.raises(SystemExit, match="new TRIES row"):
        n17.main(["--out", str(tmp_path)])
    assert (tmp_path / "n17_cfm_gate.json").read_text() == "{}"


def _numbers(text: str) -> set[int]:
    """The whole numbers in a report's text: decimals, hexadecimal hashes and ISO times are not counts of anything."""
    text = re.sub(r"\b[0-9a-f]{7,}\b", " ", text)
    text = re.sub(r"\d{4}-\d{2}-\d{2}T[\d:]+Z", " ", text)
    return {int(d) for d in re.findall(r"(?<![\d.])\d+(?![\d.])", text)}


def test_a_run_on_synthetic_stages_writes_no_cell_count_and_is_stamped(n17, tmp_path, monkeypatch, capsys):
    calls = _patch(n17, monkeypatch, tmp_path)
    out = tmp_path / "out"
    n17.main(["--out", str(out)])
    printed = capsys.readouterr().out
    # the pairs' two stages each, on the right panel, deduplicated, and the targets' truths nowhere read as such
    assert sorted((s, g[0]) for s, g, _ in calls) == sorted(
        [("E9.5", "g20"), ("E10.5", "g20"), ("E8.75", "g0"), ("E9.5", "g0"), ("E8.5", "g0"), ("E8.75", "g0")])
    assert all(dd for _, _, dd in calls)
    rep = json.loads((out / "n17_cfm_gate.json").read_text())
    txt = (out / "n17_cfm_gate.txt").read_text()
    assert txt.splitlines()[-1] == f"verdict: {rep['verdict']}" and printed.endswith(txt)
    assert rep["stamp"]["commit"] and isinstance(rep["stamp"]["tree_clean"], bool) and rep["stamp"]["utc"].endswith("Z")
    assert set(rep["inputs"]["files"]) == {"T2/E8.5.h5ad", "T2/E8.75.h5ad", "T2/E9.5.h5ad", "T2/E10.5.h5ad"}
    assert set(rep["inputs"]["panels"]) == {BOARD, VAL}
    assert rep["rule"]["band"] == [0.30, 0.95] and rep["rule"]["lineage_min_cells"] == 200
    # per lineage: a direction count and shares, the status of B0's momentum, no cell number anywhere
    judged = {key: p["judged"] for key, p in rep["pairs"].items()}
    assert set(judged) == {"test/E9.5->E10.5", "T2h/<=E9.5->E10.5", "T2h/<=E8.75->E10.5"}
    for p in rep["pairs"].values():
        by = {r["lineage"]: r for r in p["lineages"]}
        assert not by["Blood"]["eligible"] and by["Blood"]["directions"] == 0  # too few cells: B0's factor, reported
        assert by["Untracked@end"]["filter"] is None and by["Untracked@end"]["momentum"].startswith("skipped")
        assert by["CM"]["eligible"] and by["CM"]["directions"] >= 1  # the planted programme clears the edge
        assert 0.0 <= by["CM"]["retained_energy"] <= 1.0
        assert p["judged"]["passes"] == n17.judge_pair(p["lineages"])["passes"]
    keys, ints, texts = set(), set(), []
    _leaves(rep, keys, ints, texts)
    assert not keys & {"n", "n_obs", "n_cells", "n_ref", "n_start", "n_target", "cells", *SECRET}
    cell_numbers = {int(n * f) for n in SIZES.values() for f in FILES.values()}  # every lineage's cells at every stage
    assert not ints & cell_numbers, sorted(ints & cell_numbers)
    assert ints <= set(range(0, 300))  # the direction counts and the rule's smallest lineage
    for text in (*texts, txt, printed):
        assert not _numbers(text) & cell_numbers, text
        assert not _numbers(text) & set(SECRET.values()), text
