"""experiments/harness/diag_emb_coupling.py --test-panel (row N1 of experiments/harness/TRIES.md), on synthetic boards
and stages only: the option reads the test board's own panel and the test's two inputs, nothing it writes or prints is
a cell count, and its verdict reads the ranges the ledger recorded."""

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

from vec.boards import Board
from vec.harness import BUILD_SEED_OFFSET
from vec.io import COORD_KEY
from vec.t2_build import INTERP_N_ANCHOR

REPO = Path(__file__).resolve().parents[1]
ALL_GENES = [f"g{i}" for i in range(24)]
VAL_GENES = ALL_GENES[:20]
TEST_GENES = ALL_GENES[3:22]  # another panel, as the test board's is
SECRET = {"min_cells": 777, "max_cells": 4321, "truth_cells": 2468, "ref_cells": 1357}  # stand for a test board's
SIZES = {"E6.75": 229, "E7.25": 251, "E7.5": 263, "E8.0": 241}  # the synthetic stages' cells
COUNT_KEYS = {"cells", "n", "n_obs", "n_before", "n_after", "n_jittered", "n_cells", "n_anchor", "coupling_nonzero",
              *SECRET}
VERDICTS = ("hold", "risk note", "moved, within the note threshold")


@pytest.fixture(scope="module")
def dc():
    spec = importlib.util.spec_from_file_location("harness_diag_emb_coupling",
                                                  REPO / "experiments" / "harness" / "diag_emb_coupling.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def _stage(name: str) -> ad.AnnData:
    """A stage as load_stage returns it, x its thinnest axis and its expression left-right symmetric (|x| only)."""
    n = SIZES[name]
    C = np.random.default_rng(n).normal(size=(n, 3)) * [0.5, 1.5 + n / 200, 1.0]
    centres = np.random.default_rng(0).normal(size=(len(ALL_GENES), 3))
    E = np.exp(-((np.abs(C)[:, None, :] - np.abs(centres)[None]) ** 2).sum(-1))
    a = ad.AnnData(X=np.log1p(10 * E).astype(np.float32), obs=pd.DataFrame(index=[f"c{i}" for i in range(n)]),
                   var=pd.DataFrame(index=ALL_GENES))
    a.obsm[COORD_KEY] = C
    merged = np.where(C[:, 2] > 0, "Mesoderm", "Ectoderm").astype(object)
    merged[::9] = "Untracked"
    a.obs["merged"] = merged
    return a


def _patch(dc, monkeypatch) -> list[tuple[str, list[str], bool]]:
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

    monkeypatch.setattr(dc, "load_boards", lambda: boards)
    monkeypatch.setattr(dc, "load_stage", load_stage)
    monkeypatch.setattr(dc, "_git", lambda *args: "" if args[0] == "status" else "cafebabe")
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


def _numbers(text: str) -> set[int]:
    """Every whole number in a text, thousands separators dropped."""
    return {int(d) for d in re.findall(r"\d+", re.sub(r"(?<=\d),(?=\d{3})", "", text))}


def test_test_panel_reads_the_test_boards_panel_and_writes_no_cell_count(dc, tmp_path, monkeypatch, capsys):
    calls = _patch(dc, monkeypatch)
    monkeypatch.setattr(dc, "SEEDS", (0, 1))  # a full run in two draws, the s1 seed's included
    seen = []
    interpolate = dc.interpolate_clouds

    def spy(*args, **kw):
        seen.append((args[2].shape[1], kw.get("move_coupling", dc.MOVE_COUPLINGS[0])))
        return interpolate(*args, **kw)

    monkeypatch.setattr(dc, "interpolate_clouds", spy)
    dc.main(["--test-panel", "--out", str(tmp_path)])
    printed = capsys.readouterr().out
    # the test board's own panel, on the test's two inputs only, with the entry's moves
    assert sorted(name for name, _, _ in calls) == ["E7.5", "E8.0"]
    assert all(genes == TEST_GENES and dedup for _, genes, dedup in calls)
    assert seen == [(len(TEST_GENES), "fgw")] * 2
    rep = json.loads((tmp_path / "diag_emb_coupling_testpanel.json").read_text())
    txt = (tmp_path / "diag_emb_coupling_testpanel.txt").read_text()
    assert rep["genes_board"] == "T2:embryo:test_interp" and rep["n_genes"] == len(TEST_GENES)
    assert [d["build_seed"] for d in rep["draws"]] == [BUILD_SEED_OFFSET, BUILD_SEED_OFFSET + 1]
    assert rep["s1_draw"]["build_seed"] == BUILD_SEED_OFFSET + 1
    assert rep["values"]["s1_thin_ratio"] == [rep["s1_draw"]["fgw"]["thin_ratio"]]
    assert rep["verdict"] in VERDICTS and set(rep["departure"]) == set(dc.RECORDED)
    assert txt.splitlines()[-1] == f"verdict: {rep['verdict']}" and printed.endswith(txt)
    assert not (tmp_path / "diag_emb_coupling.json").exists()
    # nothing written or printed is a cell count: no count field, whole numbers only seeds and the panel's size
    keys, ints, texts = set(), set(), []
    _leaves(rep, keys, ints, texts)
    assert not keys & COUNT_KEYS
    assert ints <= {*range(10), *(BUILD_SEED_OFFSET + s for s in range(10)), len(TEST_GENES)}
    forbidden = {*SECRET.values(), *SIZES.values(), dc.N_DRAW, dc.N_DRAW // 2, dc.NN_SAMPLE, INTERP_N_ANCHOR}
    for text in (*texts, txt, printed):
        assert not _numbers(text) & forbidden, text


def test_a_partial_test_panel_run_gives_no_verdict(dc, tmp_path, monkeypatch, capsys):
    _patch(dc, monkeypatch)
    dc.main(["--test-panel", "--out", str(tmp_path), "--seeds", "0"])
    rep = json.loads((tmp_path / "diag_emb_coupling_testpanel.json").read_text())
    assert rep["verdict"] == "partial run: no verdict" and rep["departure"] is None and rep["s1_draw"] is None
    assert capsys.readouterr().out.rstrip().endswith("verdict: partial run: no verdict")


def test_without_the_option_every_pair_keeps_the_backtests_panel(dc, tmp_path, monkeypatch, capsys):
    calls = _patch(dc, monkeypatch)
    dc.main(["--out", str(tmp_path), "--seeds", "0"])
    assert sorted(name for name, _, _ in calls) == ["E6.75", "E7.25", "E7.5", "E8.0"]
    assert all(genes == VAL_GENES for _, genes, _ in calls)
    rep = json.loads((tmp_path / "diag_emb_coupling.json").read_text())
    assert rep["genes_board"] == "T2:embryo:val_interp" and len(rep["pairs"]) == 3
    assert rep["verdict"] == "partial run: no verdict"
    assert not (tmp_path / "diag_emb_coupling_testpanel.json").exists()


# The 2026-11-03 run's own extremes on the test pair (results/diag_emb_coupling/, 5baba46), which TRIES.md records
RUN_1103 = {"xcorr_min": [-0.5035, -0.2613], "thin": [0.2058, 0.2490], "s1_xcorr_min": [-0.3129],
            "s1_thin_ratio": [0.6736], "s2_pool_thin": [0.4697, 0.5062]}


@pytest.mark.parametrize("change, verdict", [
    ({}, "hold"),
    ({"xcorr_min": [-0.5349, -0.2613]}, "hold"),  # -0.53: 0.03 beyond its range
    ({"xcorr_min": [-0.5351, -0.2613]}, "moved, within the note threshold"),  # -0.54: 0.04
    ({"thin": [0.15, 0.25]}, "moved, within the note threshold"),  # far out, but the note reads two values only
    ({"s1_thin_ratio": [0.6249]}, "moved, within the note threshold"),  # 0.62: 0.05, not more
    ({"s1_thin_ratio": [0.6149]}, "risk note"),  # 0.61: 0.06
    ({"s2_pool_thin": [0.4697, 0.5551]}, "moved, within the note threshold"),  # 0.56: 0.05
    ({"s2_pool_thin": [0.4697, 0.5651]}, "risk note"),  # 0.57: 0.06
    ({"s2_pool_thin": [0.4149, 0.5062]}, "risk note"),  # 0.41: 0.06 below
])
def test_n1_verdict_reads_the_recorded_ranges_at_two_decimals(dc, change, verdict):
    got, dep = dc.n1_verdict({**RUN_1103, **change})
    assert got == verdict and set(dep) == set(dc.RECORDED) and min(dep.values()) >= 0
    if not change:
        assert max(dep.values()) == 0
    if "s1_thin_ratio" in change:
        assert dep["s1_thin_ratio"] == pytest.approx(0.67 - float(f"{change['s1_thin_ratio'][0]:.2f}"))


def test_recorded_values_are_the_ledgers(dc):
    text = (REPO / "experiments" / "harness" / "TRIES.md").read_text()
    flat = " ".join(text.split())
    diag = flat[flat.index("2026-11-03, the embryo coupling"):]
    cells = [c.strip() for c in re.search(r"\| E7\.5, E8\.0 -> E7\.75 \(1/2, test\) \|(.*)", text).group(1).split("|")]
    R, also = dc.RECORDED, dc.ALSO_RECORDED
    assert cells[0] == f"{R['xcorr_min'][1]:.2f} to {R['xcorr_min'][0]:.2f}"
    assert cells[2] == f"{R['thin'][0]:.2f}-{R['thin'][1]:.2f}"
    assert cells[5].startswith(f"{also['nn15_kept'][0]:.2f}-{also['nn15_kept'][1]:.2f} /")
    assert (f"cross-correlation {R['s1_xcorr_min'][0]:.2f}, thinnest axis {R['s1_thin_ratio'][0]:.2f} of the "
            "thinner bracket's") in diag
    assert f"thinnest axis is {R['s2_pool_thin'][0]:.2f}-{R['s2_pool_thin'][1]:.2f} with fgw's moves" in diag
    assert f"({also['thin_ratio'][0]:.2f}-{also['thin_ratio'][1]:.2f} of the thinner bracket's)" in flat
