"""Try H30 (experiments/harness/h30_heart_late_chord.py): synthetic arrays and hand-made tables only, and committed
result files of earlier rows read for arithmetic. Nothing here loads a competition file, an outside file or scores a
real cell; the launcher is traced on its refusal paths only, under a temporary home."""

from __future__ import annotations

import ast
import importlib.util
import inspect
import json
import os
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from vec import harness
from vec.adoption import SEEDS

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "experiments" / "harness" / "h30_heart_late_chord.py"
LAUNCHER = REPO_ROOT / "experiments" / "harness" / "run_h30.sh"


@pytest.fixture(scope="module")
def h30():
    spec = importlib.util.spec_from_file_location("h30_heart_late_chord_module", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


# ------------------------------------------------------------------------------------------------ synthetic counts
def _counts(n: int, profile: np.ndarray, rng, depth: float = 400.0) -> sparse.csr_matrix:
    p = np.exp(profile - profile.max())
    p /= p.sum()
    d = np.maximum(np.round(np.exp(rng.normal(np.log(depth), 0.3, n))), 40)[:, None]
    return sparse.csr_matrix(rng.poisson(p[None, :] * d))


def _gate_sets(rng, delta_m, delta_j, n=300, genes=250):
    base = rng.normal(0, 1.0, genes)
    bias = rng.normal(0, 0.5, genes)  # a platform's gene-by-gene detection, which the within-platform change cancels
    return (_counts(n, base, rng), _counts(n, base + delta_m, rng), _counts(n, base + bias, rng),
            _counts(n, base + bias + delta_j, rng))


def test_the_gate_passes_an_agreeing_change_and_fails_an_unrelated_one(h30):
    rng = np.random.default_rng(1)
    d = rng.normal(0, 1.0, 250)
    agree = h30.gate_lineage("CM", *_gate_sets(rng, d, d + rng.normal(0, 0.3, 250)))
    assert agree["evaluable"] and agree["passes"] and agree["rho"] >= 0.35 and agree["rho"] > agree["null_p95"]
    assert agree["genes"] >= 100 and -0.2 < agree["null_mean"] < 0.2
    assert agree["ceiling"]["merfish"] > 0.5 and agree["ceiling"]["jax"] > 0.5
    unrelated = h30.gate_lineage("CM", *_gate_sets(rng, d, rng.normal(0, 1.0, 250)))
    assert unrelated["evaluable"] and not unrelated["passes"]
    other = h30.gate_lineage("Endo", *_gate_sets(rng, d, d + rng.normal(0, 1.2, 250)))
    assert other["threshold"] == 0.25 and agree["threshold"] == 0.35


def test_the_gate_lines_at_their_edges(h30):
    rng = np.random.default_rng(2)
    d = rng.normal(0, 1.0, 250)
    sets = list(_gate_sets(rng, d, d))
    sets[3] = sets[3][:199]  # JAX E9.5 one cell short of 200
    r = h30.gate_lineage("CM", *sets)
    assert not r["evaluable"] and not r["passes"] and r["cells"]["jax_E9.5"] == 199
    sets = list(_gate_sets(rng, d, d))
    sets[0] = sparse.csr_matrix(sets[0].toarray() * (np.arange(250) < 90))  # only 90 genes can reach the floor
    r = h30.gate_lineage("CM", *sets)
    assert r["evaluable"] and r["too_few_genes"] and not r["passes"]
    # the outcome: VOID on the classifier or the gene floor, PASS on CM and one other, else FAIL
    p = {"passes": True, "too_few_genes": False}
    f = {"passes": False, "too_few_genes": False}
    assert h30.gate_outcome({"CM": p, "Endo": p, "Epi": f, "Hep": f}, 0.80) == ("PASS", [])
    assert h30.gate_outcome({"CM": p, "Endo": f, "Epi": f, "Hep": f}, 0.9)[0] == "FAIL"
    assert h30.gate_outcome({"CM": f, "Endo": p, "Epi": p, "Hep": p}, 0.9)[0] == "FAIL"
    assert h30.gate_outcome({"CM": p, "Endo": p}, 0.7999)[0] == "VOID"
    assert h30.gate_outcome({"CM": p, "Endo": p, "Hep": {"passes": False, "too_few_genes": True}}, 0.9)[0] == "VOID"
    assert h30.gate_outcome({"CM": p, "Endo": p, "Hep": {"passes": False, "evaluable": False}}, 0.9)[0] == "PASS"


def test_the_rho_line_needs_both_the_threshold_and_the_null(h30, monkeypatch):
    rng = np.random.default_rng(3)
    d = rng.normal(0, 1.0, 250)
    sets = _gate_sets(rng, d, d)
    monkeypatch.setattr(h30, "spearman", lambda a, b: 0.35)
    monkeypatch.setattr(h30, "binned_null", lambda *a, **k: np.full(1000, 0.0))
    assert h30.gate_lineage("CM", *sets)["passes"]  # exactly 0.35 and above the null
    monkeypatch.setattr(h30, "binned_null", lambda *a, **k: np.full(1000, 0.35))
    assert not h30.gate_lineage("CM", *sets)["passes"]  # not above the null's 95th percentile
    monkeypatch.setattr(h30, "spearman", lambda a, b: 0.3499)
    monkeypatch.setattr(h30, "binned_null", lambda *a, **k: np.full(1000, 0.0))
    assert not h30.gate_lineage("CM", *sets)["passes"]
    monkeypatch.setattr(h30, "spearman", lambda a, b: 0.25)
    assert h30.gate_lineage("Epi", *sets)["passes"]


def test_the_null_permutes_within_level_bins_and_is_fixed_by_its_seed(h30):
    rng = np.random.default_rng(4)
    x, y = rng.normal(size=300), rng.normal(size=300)
    basis = np.arange(300.0)  # sorted: bin b holds genes 30b to 30b + 29
    a = h30.binned_null(x, y, basis, "CM", n=50)
    np.testing.assert_array_equal(a, h30.binned_null(x, y, basis, "CM", n=50))
    assert not np.array_equal(a, h30.binned_null(x, y, basis, "Endo", n=50))
    # x equal to the bin index: any within-bin permutation leaves it, so every null rho is x's own
    xb = np.repeat(np.arange(10.0), 30)
    assert np.allclose(h30.binned_null(xb, y, basis, "CM", n=20), h30.spearman(xb, y))


def test_the_classifier_and_its_reject_option(h30, monkeypatch):
    rng = np.random.default_rng(5)
    groups = ["CM", "Endocardium", "Endothelium", "Neural"]
    profiles = {g: rng.normal(0, 1.5, 60) for g in groups}
    C = sparse.vstack([_counts(150, profiles[g], rng) for g in groups]).tocsr()
    labels = np.repeat(groups, 150)
    model, rec = h30.fit_classifier(C, labels)
    assert rec["balanced_accuracy"] > 0.9 and rec["n_held_out"] == 120 and rec["n_train"] == 480
    assert rec["classes"] == sorted(groups)
    lin, keep = h30.classify(model, C)
    assert set(lin[keep]) <= {"CM", "Endo", ""} and (lin[labels == "Neural"] == "").all()
    monkeypatch.setattr(h30, "CHUNK", 7)  # the chunking changes nothing
    lin7, keep7 = h30.classify(model, C)
    np.testing.assert_array_equal(lin7, lin)
    np.testing.assert_array_equal(keep7, keep)
    np.testing.assert_array_equal(h30.lineage_of(["Endocardium", "Endothelium", "EpicardiumST", "Hepatocyte", "LPM"]),
                                  ["Endo", "Endo", "Epi", "Hep", ""])

    class Fixed:
        classes_ = np.array(["CM", "Neural"])

        def predict_proba(self, F):
            return np.array([[0.7, 0.3], [0.6999, 0.3001], [0.9, 0.1], [0.2, 0.8]])

    X = sparse.csr_matrix(np.array([[30, 0], [40, 0], [29, 0], [50, 0]]))
    lin, keep = h30.classify(Fixed(), X)
    assert keep.tolist() == [True, False, False, True]  # 0.7 kept, 0.6999 not; a library of 29 not, 30 kept
    assert lin.tolist() == ["CM", "CM", "CM", ""]


def test_the_chord_within_each_source_its_floor_average_and_clip(h30):
    e = sparse.csr_matrix(np.array([[100, 100, 0, 9800, 0], [100, 100, 0, 9800, 0]]))
    late = sparse.csr_matrix(np.array([[400, 25, 0, 9575, 0], [400, 25, 0, 9575, 0]]))
    d = h30.source_chord(e, late, np.array([True, True, True, True, False]))
    assert d[0] == pytest.approx(np.log(400.5) - np.log(100.5))
    assert d[1] == pytest.approx(np.log(25.5) - np.log(100.5))
    assert d[2] == 0.0 and np.isnan(d[4])  # below 0.5 at both anchors: 0; not carried: NaN
    tiny_e = sparse.csr_matrix(np.array([[0, 4, 0, 99996, 0]]))  # gene 1 at 0.4 CP10k, later 0.3
    tiny_l = sparse.csr_matrix(np.array([[0, 3, 0, 99997, 0]]))
    assert h30.source_chord(tiny_e, tiny_l, np.ones(5, bool))[1] == 0.0
    a = np.array([0.2, np.nan, 2.0, np.nan, -3.0])
    b = np.array([0.4, 0.6, np.nan, np.nan, -0.5])
    # gene 4: the mean -1.75, then clipped at -log 3
    np.testing.assert_allclose(h30.combine([a, b], 5), [0.3, 0.6, np.log(3.0), 0.0, -np.log(3.0)])
    np.testing.assert_array_equal(h30.combine([], 3), np.zeros(3))


def _cells(n, genes=5):
    return sparse.csr_matrix(np.ones((n, genes), dtype=np.int64) * 10)


def test_a_source_contributes_only_with_200_kept_cells_at_each_of_its_anchors(h30):
    def lins(**n):
        return {k: _cells(n.get(k, 0)) for k in h30.LINEAGES}
    carried = {"jax": np.ones(5, bool), "feng": np.ones(5, bool)}
    # (i) Feng under 200 is not averaged; (v) Endo on JAX alone; Hep and Epi with no source at 200
    cells = {"jax": {"early": lins(CM=250, Endo=300, Hep=150, Epi=10), "late": lins(CM=300, Endo=300, Hep=500)},
             "feng": {"early": lins(CM=50, Hep=100, Epi=199), "late": lins(CM=100, Hep=100, Epi=300)}}
    D, rec = h30.chord_table(cells, carried, 5)
    assert rec["CM"]["moved"] and rec["CM"]["sources"] == ["jax"]
    assert rec["Endo"]["moved"] and rec["Endo"]["sources"] == ["jax"]
    assert not rec["Hep"]["moved"] and not rec["Epi"]["moved"] and set(D) == {"CM", "Endo"}
    assert rec["Hep"]["total"] == {"early": 250, "late": 600}  # summed over 200, yet no source holds 200 at each
    # (ii) JAX one cell short at its early anchor: Feng alone
    cells2 = {"jax": {"early": lins(CM=199), "late": lins(CM=300)},
              "feng": {"early": lins(CM=300), "late": lins(CM=300)}}
    r = h30.chord_table(cells2, carried, 5)[1]["CM"]
    assert r["moved"] and r["sources"] == ["feng"]
    # (iii) 150 and 150 in each source: the summed totals are 300, but no source contributes
    cells3 = {"jax": {"early": lins(CM=150), "late": lins(CM=150)},
              "feng": {"early": lins(CM=150), "late": lins(CM=150)}}
    D3, rec3 = h30.chord_table(cells3, carried, 5)
    assert rec3["CM"]["total"] == {"early": 300, "late": 300} and not rec3["CM"]["moved"] and "CM" not in D3
    assert rec3["CM"]["sources"] == []
    # (iv) the 1.75-day backtest: Feng has no early anchor and never contributes; its 900 late cells move nothing
    one = {"jax": {"early": lins(CM=250), "late": lins(CM=150)}, "feng": {"late": lins(CM=900)}}
    D4, rec4 = h30.chord_table(one, carried, 5)
    r = rec4["CM"]
    assert r["total"] == {"early": 250, "late": 1050} and r["sources"] == [] and not r["moved"] and "CM" not in D4


def test_the_code_guards_refuse_before_any_file_is_opened(h30, tmp_path):
    nowhere = tmp_path / "absent"
    with pytest.raises(PermissionError, match="forecast span"):
        h30.read_jax(nowhere, (13, 14, 22), ["g"], h30.REFUSED_JAX["T2h/<=E8.75->E10.5"])
    with pytest.raises(PermissionError, match="forecast span"):
        h30.read_jax(nowhere, (30,), ["g"], h30.REFUSED_JAX["T2h/<=E9.5->E10.5"])
    with pytest.raises(PermissionError, match="forecast span"):
        h30.read_feng(nowhere, ("E9.5", "E14.5"), ["g"], h30.REFUSED_FENG["T2h/<=E8.75->E10.5"])
    with pytest.raises(PermissionError, match="window"):
        h30.read_feng(nowhere, ("E11.5",), ["g"])
    with pytest.raises(PermissionError, match="E10.5"):
        h30.merfish_stage(nowhere, "T2/E10.5.h5ad", ["g"])
    assert h30.REFUSED_JAX["T2h/<=E8.75->E10.5"] >= set(range(15, 23))
    assert not ({22, 54, 55, 56} & h30.REFUSED_JAX["T2h/<=E9.5->E10.5"])
    assert h30.TRANSITIONS["T2h/<=E8.75->E10.5"]["jax_early"] == (13, 14)
    assert h30.TRANSITIONS["T2h/<=E8.75->E10.5"]["feng_early"] is None
    assert h30.TRANSITIONS["T2h/<=E9.5->E10.5"]["tau"] == 1 / 5
    assert h30.TRANSITIONS["T2h/<=E8.75->E10.5"]["tau"] == 1.75 / 5.75 and h30.TEST_TAU == 0.4


def test_the_jax_cell_table_gives_its_days_and_the_list_is_checked(h30, tmp_path):
    p = tmp_path / "df_cell.csv"
    rows = [("c1", "embryo_13", "E8.75", "author type A"), ("c2", "embryo_22", "E9.5", "B"),
            ("c3", "embryo_54", "E14.25", "C"), ("c4", "embryo_9", "E8.5", "D"), ("c5", "embryo_16", "E9.0", "E")]
    pd.DataFrame(rows, columns=["cell_id", "embryo_id", "day", "celltype"]).to_csv(p, index=False)
    days = h30.jax_days(p, (13, 16, 22, 54))
    assert days == {13: ["E8.75"], 16: ["E9.0"], 22: ["E9.5"], 54: ["E14.25"]}
    assert h30.check_jax_days(days) == {13: 8.75, 16: 9.0, 22: 9.5, 54: 14.25}
    for bad in ({13: ["E8.5"]}, {22: ["E9.5", "E9.75"]}, {16: ["E9.5"]}, {9: ["E8.5"]}, {55: []}):
        with pytest.raises(ValueError):
            h30.check_jax_days(bad)
    assert h30.check_jax_days({55: ["E14.333"], 56: ["E14.75"]}) == {55: 14.333, 56: 14.75}


def test_gene_matching_by_symbol(h30, tmp_path):
    p = tmp_path / "df_gene.csv"
    pd.DataFrame({"gene_ID": ["E1", "E2", "E3", "E4"], "gene_short_name": ["Myl7", "Tnnt2", "Myl7", "Nppa"],
                  "chr": ["1", "2", "3", "4"]}).to_csv(p, index=False)
    ids = h30.jax_gene_symbols(p)
    panel = ["Nppa", "Myl7", "Tnnt2", "Hand1"]
    syms, rec = h30.symbols_of_rows(["E1", "E2", "E3", "E4"], ids)
    assert rec == {"rows": 4, "listed": 4, "unlisted": 0}
    cols = h30.panel_columns(syms, panel, ids)
    assert cols == {"Nppa": 3, "Tnnt2": 1}  # Myl7 is held twice: not used
    # row names that are symbols are never matched as symbols: refused, since no row is a table gene ID
    with pytest.raises(ValueError, match="none of the count matrix's row names"):
        h30.symbols_of_rows(["Nppa", "Xyz"], ids)
    with pytest.raises(ValueError):
        h30.symbols_of_rows(["Q1", "Q2"], ids)
    C = sparse.csr_matrix(np.array([[1, 2, 3, 4], [5, 0, 0, 8]]))
    P, carried = h30.to_panel(C, cols, panel)
    assert carried.tolist() == [True, False, True, False]
    np.testing.assert_array_equal(P.toarray(), [[4, 0, 2, 0], [8, 0, 0, 0]])
    with pytest.raises(ValueError):
        h30.to_panel(sparse.csr_matrix(np.array([[1, 1, 1, 1.5]])), cols, panel)


# A gene table as JAX's is laid out ("", gene_id, gene_type, gene_short_name, chr): Myl7 held by two IDs.
_TABLE = {"": ["0", "1", "2", "3", "4"], "gene_id": ["ENSG1", "ENSG2", "ENSG3", "ENSG4", "ENSG5"],
          "gene_type": ["protein_coding"] * 5, "gene_short_name": ["Nppa", "Myl7", "Myl7", "Tnnt2", "Hand1"],
          "chr": ["1", "2", "3", "4", "5"]}
_PANEL = ["Nppa", "Myl7", "Tnnt2", "Hand1", "Gata4"]
# matrix rows: three table IDs (ENSG3, Myl7's second ID, absent), then three rows the table does not list: an ID
# outside it, a panel symbol given as a row name, and a version-suffixed table ID (no stripping)
_ROWS = ["ENSG1", "ENSG2", "ENSG4", "ENSGX9", "Gata4", "ENSG5.2"]
# genes x cells: cell 0 has 6 listed panel counts and 300 unlisted ones; cell 1 has 30 and 300
_X = np.array([[4, 20], [50, 7], [2, 10], [100, 100], [100, 100], [100, 100]])


def _fake_rds(h30, monkeypatch, matrices: dict[str, tuple[np.ndarray, list[str]]]):
    """Each rds path's name -> (genes x cells counts, row names); the manifest check and the rds parser replaced."""
    monkeypatch.setattr(h30.fio, "verify_manifest_output", lambda p: None)
    monkeypatch.setattr(h30, "_parse_rds", lambda p: Path(p).name)
    monkeypatch.setattr(h30.feng2022_geo, "r_class", lambda obj: "dgCMatrix")
    monkeypatch.setattr(h30.feng2022_geo, "dgc_matrix",
                        lambda obj: (sparse.csc_matrix(matrices[obj][0]), list(matrices[obj][1]), None))


def test_rows_the_gene_table_does_not_list(h30, tmp_path, monkeypatch):
    root = tmp_path / "jax"
    (root / "meta").mkdir(parents=True)
    pd.DataFrame(_TABLE).to_csv(root / "meta" / "df_gene.csv", index=False)
    pd.DataFrame({"cell_id": ["a", "b"], "embryo_id": ["embryo_13", "embryo_22"], "day": ["E8.75", "E9.5"]}).to_csv(
        root / "meta" / "df_cell.csv", index=False)
    ids = h30.jax_gene_symbols(root / "meta" / "df_gene.csv")

    # (1) exact ID matching; unlisted rows get no symbol and are never matched in their own right
    syms, rec = h30.symbols_of_rows(_ROWS, ids)
    assert syms == ["Nppa", "Myl7", "Tnnt2", None, None, None]
    assert rec == {"rows": 6, "listed": 3, "unlisted": 3}
    cols = h30.panel_columns(syms, _PANEL, ids)
    # Myl7: the table gives it two IDs, so dropped though the matrix holds only one; Hand1: its row is "ENSG5.2";
    # Gata4: a row name, not a table ID
    assert cols == {"Nppa": 0, "Tnnt2": 2}
    assert h30.symbols_of_rows(["ensg1", "ENSG4"], ids)[0] == [None, "Tnnt2"]  # no case folding
    with pytest.raises(ValueError, match="not unique"):
        h30.symbols_of_rows(["ENSG1", "ENSG1", "ENSG4"], ids)
    with pytest.raises(ValueError, match="none of"):
        h30.symbols_of_rows(["Nppa", "Tnnt2", "ENSGX9"], ids)
    dup = tmp_path / "dup.csv"
    pd.DataFrame({**_TABLE, "gene_id": ["ENSG1", "ENSG2", "ENSG3", "ENSG1", "ENSG5"]}).to_csv(dup, index=False)
    with pytest.raises(ValueError, match="more than once"):
        h30.jax_gene_symbols(dup)
    # Feng is unchanged: by symbol from its own list, a symbol held twice dropped
    assert h30.panel_columns(["Myl7", "Nppa", "Nppa", "Gata4"], _PANEL) == {"Myl7": 0, "Gata4": 3}

    # the reader end to end: the same counts with and without the unlisted rows give the same panel counts
    listed_only = (_X[:3], _ROWS[:3])
    _fake_rds(h30, monkeypatch, {"embryo_13_gene_count.rds": (_X, _ROWS),
                                 "embryo_22_gene_count.rds": listed_only})
    out, carried, jrec = h30.read_jax(root, (13, 22), _PANEL)
    assert carried.tolist() == [True, False, True, False, False]
    assert jrec["matrix_rows"] == {"13": {"rows": 6, "listed": 3, "unlisted": 3},
                                   "22": {"rows": 3, "listed": 3, "unlisted": 0}}
    np.testing.assert_array_equal(out[13].toarray(), [[4, 0, 2, 0, 0], [20, 0, 10, 0, 0]])
    np.testing.assert_array_equal(out[13].toarray(), out[22].toarray())
    gi = np.flatnonzero(carried)
    P = out[13][:, gi]

    # (2) the panel library: 6 (rejected) and 30 (kept); the unlisted rows' 300 counts add nothing
    class Model:
        classes_ = np.array(["CM"])

        def predict_proba(self, F):
            return np.ones((F.shape[0], 1))
    lin, keep = h30.classify(Model(), P)
    assert keep.tolist() == [False, True] and lin.tolist() == ["CM", "CM"]

    # (3) CP10k totals are the panel library over the carried genes
    np.testing.assert_allclose(h30.mean_cp10k(P), [(4 / 6 + 20 / 30) / 2 * 1e4, (2 / 6 + 10 / 30) / 2 * 1e4])
    np.testing.assert_allclose(h30.log1p_cp10k(P).toarray(),
                               np.log1p([[4 / 6 * 1e4, 2 / 6 * 1e4], [20 / 30 * 1e4, 10 / 30 * 1e4]]), rtol=1e-6)

    # (4) the chord: NaN where JAX carries no gene, then taken from the other source where it holds it
    late = sparse.csr_matrix(np.array([[8, 0, 2, 0, 0], [40, 0, 10, 0, 0]]))
    dj = h30.source_chord(out[13], late, carried)
    assert np.isnan(dj).tolist() == [False, True, False, True, True]
    feng = np.array([0.1, 0.2, 0.3, np.nan, 0.5])
    np.testing.assert_allclose(h30.combine([dj, feng], 5), [(dj[0] + 0.1) / 2, 0.2, (dj[2] + 0.3) / 2, 0.0, 0.5])

    # every embryo must carry the same panel genes
    _fake_rds(h30, monkeypatch, {"embryo_13_gene_count.rds": (_X, _ROWS),
                                 "embryo_22_gene_count.rds": (_X[:2], _ROWS[:2])})
    with pytest.raises(ValueError, match="different panel genes"):
        h30.read_jax(root, (13, 22), _PANEL)


# ------------------------------------------------------------------------------------------------ the steps' wiring
class _World:
    """Synthetic stages for the gate and the tables: every outside or released reader replaced by these arrays."""

    PANEL = [f"p{i}" for i in range(160)]
    GROUPS = ("CM", "Endocardium", "Neural", "EpicardiumST")

    def __init__(self, seed=11):
        rng = np.random.default_rng(seed)
        self.rng = rng
        self.base = {g: rng.normal(0, 1.0, 160) for g in self.GROUPS}
        self.delta = {g: rng.normal(0, 0.8, 160) for g in self.GROUPS}  # E8.75 -> E9.5, both platforms
        self.late = {g: rng.normal(0, 0.8, 160) for g in self.GROUPS}  # E9.5 -> E14.5
        self.bias = rng.normal(0, 0.4, 160)
        self.calls = []

    def stage(self, t, n=260, bias=False, groups=None):
        gs = groups or self.GROUPS
        b = self.bias if bias else 0.0
        C = sparse.vstack([_counts(n, self.base[g] + b + t[0] * self.delta[g] + t[1] * self.late[g], self.rng)
                           for g in gs]).tocsr()
        return C, np.repeat(gs, n)

    def merfish_stage(self, data_dir, f, genes):
        self.calls.append(("merfish", f))
        assert f != "T2/E10.5.h5ad"
        return self.stage((0.0, 0.0) if f == "T2/E8.75.h5ad" else (1.0, 0.0))

    def read_jax(self, root, embryos, panel, refused=frozenset()):
        assert not set(embryos) & set(refused)
        self.calls.append(("jax", tuple(sorted(embryos))))
        t = {13: (0, 0), 14: (0, 0), 22: (1, 0), 54: (1, 1), 55: (1, 1), 56: (1, 1), **{n: (0.5, 0) for n in
                                                                                       range(15, 22)}}
        carried = np.arange(160) >= 3  # JAX carries all but three panel genes
        out = {}
        for n in embryos:
            C = self.stage(t[n], n=230, bias=True)[0].tolil()
            C[:, :3] = 0
            out[n] = C.tocsr()
        days = {n: h for n, h in ((13, 8.75), (14, 8.75), (22, 9.5), (54, 14.25), (55, 14.333), (56, 14.75))
                if n in embryos}
        return out, carried, {"embryos": sorted(embryos), "days": {str(k): v for k, v in days.items()},
                              "matrix_rows": {str(n): {"rows": 170, "listed": 160, "unlisted": 10} for n in embryos}}

    def read_feng(self, root, stages, panel, refused=frozenset()):
        assert not set(stages) & set(refused)
        self.calls.append(("feng", tuple(stages)))
        t = {"E9.5": (1, 0), "E14.5": (1, 1)}
        return ({s: self.stage(t[s], n=300, groups=("CM", "Endocardium"))[0] for s in stages},
                np.ones(160, bool), {"stages": list(stages)})


ROW_C = "| H30 (x) | heart | ... named in this row before the read: C. Changes no file |"


def _wire(h30, monkeypatch, world):
    class Task:
        panel_board = "pb"

        def __init__(self, key):
            self.reference = h30.TRANSITIONS.get(key, {}).get("merfish", "T2/E9.5.h5ad")

    monkeypatch.setattr(h30, "tries_rows", lambda path=h30.TRIES: {"H30", "N56"})
    monkeypatch.setattr(h30, "row_text", lambda *a, **k: ROW_C)
    monkeypatch.setattr(h30, "register_names_row", lambda *a, **k: [])
    monkeypatch.setattr(h30.provenance, "git_state", lambda repo, ignore=(): ("c" * 40, []))
    monkeypatch.setattr(h30.payoff, "provenance", lambda: {"commit": "c" * 40, "tree_clean": True})
    monkeypatch.setattr(h30.payoff, "heart_task", Task)
    monkeypatch.setattr(h30, "need_rdata", lambda: None)
    monkeypatch.setattr(h30, "load_boards", lambda: {"pb": type("B", (), {"genes": world.PANEL})()})
    for f in ("merfish_stage", "read_jax", "read_feng"):
        monkeypatch.setattr(h30, f, getattr(world, f))


def test_the_gate_step_end_to_end_on_synthetic_stages(h30, monkeypatch, tmp_path):
    world = _World()
    _wire(h30, monkeypatch, world)
    res = h30.run_gate(tmp_path / "g", tmp_path, tmp_path)
    assert res["outcome"] == "PASS" and res["genes_matched"] == 157 and res["classifier"]["balanced_accuracy"] > 0.8
    assert res["lineages"]["CM"]["passes"] and res["lineages"]["Endo"]["passes"]
    assert not res["lineages"]["Hep"]["evaluable"] and res["lineages"]["Epi"]["evaluable"]
    assert res["lineages"]["CM"]["intermediate"]["rho"] is not None
    assert sorted(c for c in world.calls if c[0] == "merfish") == [("merfish", "T2/E8.75.h5ad"),
                                                                   ("merfish", "T2/E9.5.h5ad")]
    assert ("jax", tuple(range(13, 23))) in world.calls and not any(c[0] == "feng" for c in world.calls)
    assert h30.gate_passed(tmp_path / "g" / "gate.json", "c" * 40)
    assert (tmp_path / "g" / "gate.txt").read_text().startswith("Try H30, STEP G (commit cccc")
    with pytest.raises(SystemExit, match="exists"):  # a reading is never overwritten
        h30.run_gate(tmp_path / "g", tmp_path, tmp_path)
    # a JAX change unrelated to MERFISH's fails the gate
    w2 = _World()
    w2.delta = {g: w2.rng.normal(0, 0.8, 160) for g in w2.GROUPS}
    real = w2.read_jax
    shuffled = {g: w2.rng.normal(0, 0.8, 160) for g in w2.GROUPS}

    def jax_unrelated(*a, **k):
        keep, w2.delta = w2.delta, shuffled
        try:
            return real(*a, **k)
        finally:
            w2.delta = keep
    monkeypatch.setattr(h30, "read_jax", jax_unrelated)
    assert h30.run_gate(tmp_path / "g2", tmp_path, tmp_path)["outcome"] == "FAIL"


def test_the_tables_step_end_to_end_on_synthetic_stages(h30, monkeypatch, tmp_path):
    world = _World()
    _wire(h30, monkeypatch, world)
    (tmp_path / "g").mkdir()
    (tmp_path / "g" / "gate.json").write_text(json.dumps(
        {"axis": h30.AXIS + ":gate", "outcome": "PASS", "provenance": {"commit": "c" * 40, "tree_clean": True}}))
    sums = h30.run_tables(tmp_path / "t", tmp_path / "g", tmp_path, tmp_path)
    assert set(sums) == {"e95.json", "e875.json"}
    sha = tmp_path / "t" / "tables_sha256.json"
    assert json.loads(sha.read_text()) == sums
    for key in h30.TASK_KEYS:
        short = h30.TRANSITIONS[key]["short"]
        tb, s = h30.load_table(tmp_path / "t" / f"{short}.json", sha)
        assert tb["task"] == key and tb["stages_read"] == h30._expected_table(key)
        assert tb["tau_practice"] == h30.TRANSITIONS[key]["tau"] and tb["genes"] == world.PANEL
        assert set(tb["lineages"]) >= {"CM", "Endo"} and "Hep" not in tb["lineages"]
        d = np.asarray(tb["lineages"]["CM"])
        assert np.abs(d).max() <= np.log(3) + 1e-12
        assert all(c["trained_on"] == [h30.TRANSITIONS[key]["merfish"]] for c in tb["classifiers"].values())
        part = {"table_sha256": s, "table": tb}
        assert h30.table_checks(key, part, sums) == []
    e875 = json.loads((tmp_path / "t" / "e875.json").read_text())
    assert e875["record"]["CM"]["sources"] == ["jax"]  # no Feng early anchor on the 1.75-day backtest
    e95 = json.loads((tmp_path / "t" / "e95.json").read_text())
    assert e95["record"]["CM"]["sources"] == ["feng", "jax"]
    # what was read: the reference stage only, never E10.5; Feng E9.5 only on the one-day backtest
    assert [c for c in world.calls if c[0] == "merfish"] == [("merfish", "T2/E9.5.h5ad"),
                                                             ("merfish", "T2/E8.75.h5ad")]
    assert ("feng", ("E14.5",)) in world.calls and ("feng", ("E9.5", "E14.5")) in world.calls
    assert ("jax", (13, 14)) in world.calls and ("jax", (22, 54, 55, 56)) in world.calls
    # the tables are the same bytes when rebuilt
    world2 = _World()
    _wire(h30, monkeypatch, world2)
    again = h30.run_tables(tmp_path / "t2", tmp_path / "g", tmp_path, tmp_path)
    assert again == sums
    with pytest.raises(SystemExit, match="exists"):
        h30.run_tables(tmp_path / "t", tmp_path / "g", tmp_path, tmp_path)


# ------------------------------------------------------------------------------------------------ the arms
@pytest.fixture(scope="module")
def synth():
    spec = importlib.util.spec_from_file_location("b0_trend_tests_for_h30", REPO_ROOT / "tests" / "test_b0_trend.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    inputs = {"ref": mod._stage(1600, seed=2, shift=0.0, tiny=60),
              "start": mod._stage(1100, seed=3, shift=0.8, tiny=10),
              "earlier": mod._stage(900, seed=5, shift=1.6, tiny=10)}
    group = {"CM": "CM", "Endothelial": "Endothelium", "Untracked": "Unknown", "Tiny": "Hepatocyte"}
    for a in inputs.values():
        m = np.asarray(a.obs["merged"]).astype(str)
        a.obs["core"] = m == "CM"
        a.obs["group"] = [group[x] for x in m]
    return mod, inputs


def _size(h30, inputs):
    sp = h30.h27.sp
    geoms = [sp.geometry_from_arrays(n, t, np.asarray(inputs[r].obsm["spatial_3D"], float),
                                     np.asarray(inputs[r].obs["core"], bool))
             for r, n, t in (("earlier", "E8.5", 8.5), ("start", "E8.75", 8.75), ("ref", "E9.5", 9.5))]
    return h30.h27.analogue_size_result(geoms, 10.5)


def _table(genes, seed=7):
    rng = np.random.default_rng(seed)
    return {"genes": list(genes), "tau_practice": 0.2, "tau_test": 0.4,
            "lineages": {"CM": list(np.clip(rng.normal(0, 0.8, len(genes)), -np.log(3), np.log(3))),
                         "Endo": list(np.clip(rng.normal(0, 0.8, len(genes)), -np.log(3), np.log(3)))}}


def test_the_factors_of_every_arm(h30, synth):
    mod, inputs = synth
    genes = list(inputs["ref"].var_names)
    tb = _table(genes[::-1] + ["NotOnPanel"])  # the table's own gene order; a gene off the panel is ignored
    basis = h30.reference_basis(inputs["ref"], genes, ["CM", "Endo"])
    f, fixed = h30.arm_factors(tb, genes, basis)
    d_cm = np.asarray(tb["lineages"]["CM"])[: len(genes)][::-1]
    np.testing.assert_allclose(f["A"]["CM"], np.exp(0.2 * d_cm))
    np.testing.assert_allclose(f["A2"]["CM"], np.exp(0.4 * d_cm))
    assert all((v == 1.0).all() for v in f["chord0"].values())
    for k, p in enumerate(h30.PLACEBOS, start=1):
        assert fixed[p] == {"CM": 0, "Endo": 0}
        np.testing.assert_allclose(np.sort(f[p]["CM"]), np.sort(f["A2"]["CM"]))
        perm = h30.h28.derangement(basis["CM"], "CM", bins=10, seed=h30.SEED + k)
        np.testing.assert_array_equal(f[p]["CM"], f["A2"]["CM"][perm])
        where = np.empty(len(genes), int)
        for i, b in enumerate(np.array_split(np.argsort(basis["CM"], kind="stable"), 10)):
            where[b] = i
        assert np.array_equal(where[perm], where) and not np.any(perm == np.arange(len(genes)))
    assert not np.array_equal(f["P1"]["CM"], f["P2"]["CM"])


@pytest.mark.parametrize("form", ["full_iso_half", "core_iso_half"])
def test_every_arm_keeps_the_controls_cells_positions_size_and_depth(h30, synth, form):
    mod, inputs = synth
    PC = np.asarray(inputs["ref"].obsm["spatial_3D"], float)
    keep = h30.h26.pool_masks(PC, np.asarray(inputs["ref"].obs["core"], bool))[form]
    ins = {**inputs, "ref": inputs["ref"][keep].copy()}
    genes = list(inputs["ref"].var_names)
    basis = h30.reference_basis(inputs["ref"], genes, ["CM", "Endo"])
    factors, fixed = h30.arm_factors(_table(genes), genes, basis)
    ctx = mod._ctx(ins, seed=harness.BUILD_SEED_OFFSET + 3)
    out, extra = h30.build_arms(ctx, form, _size(h30, inputs), factors, fixed)
    c = extra["checks"]
    assert set(out) == set(h30.BUILT)
    for k in ("chord0_counts_equal_control", "same_draw", "pre_size_coords_equal", "post_size_coords_equal"):
        assert c[k] is True, k
    assert len(set(c["n_cells"].values())) == 1 and c["n_cells"]["control"] == mod.BOARD.max_cells // 2
    for arm in h30.BUILT:
        np.testing.assert_array_equal(out[arm].coords, out["control"].coords)
    for arm in (h30.CHORD0, *h30.COUNTED):
        rec = extra["records"][arm]["chord"]
        assert rec["unmoved_counts_equal"] and abs(rec["mean_total_change"]) < 0.01, arm
        assert rec["moved"] == ["CM", "Endo"]
    assert not np.array_equal(out["A2"].X, out["control"].X) and not np.array_equal(out["A2"].X, out["P1"].X)
    np.testing.assert_allclose(out["chord0"].X, out["control"].X, rtol=0, atol=1e-5)


def test_only_the_moved_lineages_counts_change(h30, synth):
    mod, inputs = synth
    genes = list(inputs["ref"].var_names)
    tb = _table(genes)
    del tb["lineages"]["Endo"]
    f = h30.factor_rows(tb, genes, 0.4)
    ctx = mod._ctx(inputs, seed=harness.BUILD_SEED_OFFSET + 5)
    _, new, rows, rec = h30.chord_expression(ctx, f)
    pool = inputs["ref"]
    drawn = h30.counts_of(pool[rows], genes)
    cm = np.asarray(pool.obs["group"])[rows] == "CM"
    assert (new[~cm] != drawn[~cm]).nnz == 0 and (new[cm] != drawn[cm]).nnz > 0
    assert rec["moved"] == ["CM"] and rec["moved_cells"]["CM"] == int(cm.sum())
    # a second call gives the same counts (seed [build seed, 29])
    assert (h30.chord_expression(ctx, f)[1] != new).nnz == 0


# ------------------------------------------------------------------------------------------------ hand-made parts
METRICS = ("de_score", "de_direction", "mmd_u", "variogram", "d2_shape", "occupancy_dice", "scale_log_ratio",
           "neighborhood_mmd")
GROUPS = ("expression_change", "cell_state", "shape_growth", "local_organisation")
COMMIT = "a" * 40


def _info(h30):
    t = {q: {"cells": 1000 + i, "core_share": 0.4, "xy_scope": 1.0} for i, q in enumerate(h30.TRUTHS)}
    return t, {f: {"pool_cells": 500 + i, "core_share": 0.9, "xy_scope": 0.2} for i, f in enumerate(h30.FORMS)}


def _table_rec(h30, key):
    want = h30._expected_table(key)
    tr = h30.TRANSITIONS[key]
    return {"task": key, "tau_practice": tr["tau"], "tau_test": 0.4, "record": {}, "anchors": {}, "jax_days": {},
            "carried": {}, "stages_read": want,
            "classifiers": {s: {"trained_on": [tr["merfish"]]} for s in ("jax", "feng")}}


SHA = {"e95.json": "s" * 64, "e875.json": "t" * 64}


def _part(h30, key, seeds, score=None, de=None, groups=None, commit=COMMIT, clean=True, checks=None):
    """score: name -> per-seed board points over 50 (or world -> array); de: name -> skill added to de_score and
    de_direction; groups: name -> {group: group skill added}."""
    score, de, groups = score or {}, de or {}, groups or {}
    t_info, p_info = _info(h30)
    rows = []
    for s in seeds:
        sc, sk, gr = {}, {}, {}
        for q in h30.TRUTHS:
            w = "own" if q == h30.OWN else "CLEAN"
            sc[q], sk[q], gr[q] = {}, {}, {}
            for nm in h30.ALL:
                v = score.get(nm, 0.0)
                v = v.get(w, 0.0) if isinstance(v, dict) else v
                sc[q][nm] = 50.0 + float(np.broadcast_to(v, (10,))[s])
                sk[q][nm] = {m: 0.5 + (de.get(nm, 0.0) if m in ("de_score", "de_direction") else 0.0)
                             for m in METRICS}
                gr[q][nm] = {g: 0.5 + groups.get(nm, {}).get(g, 0.0) for g in GROUPS}
        ck = {}
        for f in h30.FORMS:
            ck[h30.LETTER[f]] = {"chord0_counts_equal_control": True, "same_draw": True,
                                 "pre_size_coords_equal": True, "post_size_coords_equal": True,
                                 "chord0_x_max_abs_diff_control": 9.5e-7,
                                 "placebo_fixed_points": {p: {"CM": 0} for p in h30.PLACEBOS}}
            ck[h30.LETTER[f]].update((checks or {}).get(h30.LETTER[f], {}))
        recs = {nm: {"size": {}, "chord": {"mean_total_change": 0.001, "unmoved_counts_equal": True}}
                for nm in h30.ARM_NAMES}
        rows.append({"seed": s, "n_cells": {nm: 350 for nm in h30.ARM_NAMES}, "checks": ck, "records": recs,
                     "scores": sc, "skills": sk, "groups": gr})
    return {"axis": h30.AXIS, "task": key, "rows": rows, "truth_info": t_info, "pred_info": p_info,
            "forecast_targets": {f: 100.0 + i for i, f in enumerate(h30.FORMS)}, "v7": {"ok": True, "log_diff": 0.0},
            "table_sha256": SHA[f"{h30.TRANSITIONS[key]['short']}.json"], "table": _table_rec(h30, key),
            "window_guard": h30.WINDOW_BYPASS, "provenance": {"commit": commit, "tree_clean": clean}}


def _parts(h30, **kw):
    return {k: [_part(h30, k, range(0, 5), **kw), _part(h30, k, range(5, 10), **kw)] for k in h30.TASK_KEYS}


def _h27(h30, off=0.0):
    t_info, p_info = _info(h30)
    rows = [{"seed": s, "scores": {q: {h30.h27.arm_name(f, "f"): 50.0 + off for f in h30.FORMS} for q in h30.TRUTHS}}
            for s in SEEDS]
    return {k: [{"rows": rows, "truth_info": t_info, "pred_info": p_info,
                 "forecast_targets": {f: 100.0 + i for i, f in enumerate(h30.FORMS)}}] for k in h30.TASK_KEYS}


def _read(h30, parts=None, h27p=None, commit=COMMIT, form="core_iso_half", sha=None):
    return h30.reading(parts if parts is not None else _parts(h30), h27p if h27p is not None else _h27(h30),
                       sha if sha is not None else SHA, commit, form)


def _win(h30, form, gain=1.2, a2=1.0, plc=0.0, de=0.01, ec=0.02):
    """A on `form` ahead of the control by `gain`; A2 by `a2`; placebos by `plc`; A2's DE skills and group skills up."""
    n = lambda a: h30.name(form, a)  # noqa: E731
    score = {n("A"): np.full(10, gain), n("A2"): np.full(10, a2), **{n(p): np.full(10, plc) for p in h30.PLACEBOS}}
    return score, {n("A2"): de}, {n("A2"): {"expression_change": ec, "cell_state": ec}}


def _merge(*wins):
    sc, de, gr = {}, {}, {}
    for s, d, g in wins:
        sc.update(s)
        de.update(d)
        gr.update(g)
    return {"score": sc, "de": de, "groups": gr}


def test_a_clean_null_run_is_read_and_is_none(h30):
    res = _read(h30)
    assert not res["void"] and res["outcome"] == "NONE" and res["follows"] == h30.FOLLOWS["NONE"]
    assert res["guards"]["V4"][h30.TASK_KEYS[0]]["C"] == {"mean": 0.0, "seed": 0.0}
    txt = h30.text(res)
    assert txt.startswith(f"Try H30 (run commit {COMMIT}). A labelled in-sample diagnostic") and "OUTCOME: NONE" in txt
    assert all(ln.startswith(h30.TAG) for ln in txt.split("TABLES")[1].splitlines()[1:])
    assert "check_externals skipped by name" in txt


def test_pass_needs_the_form_in_force_and_one_other_form(h30):
    both = _merge(_win(h30, "core_iso_half"), _win(h30, "full_iso_half"))
    res = _read(h30, _parts(h30, **both))
    assert res["outcome"] == "PASS" and res["statements"]["C"]["all"] and res["statements"]["A"]["all"]
    assert "form C" not in res["follows"].split("and on")[1] and "form A" in res["follows"]
    assert "ONE priced question" in res["follows"] and "the file in force stands" in res["follows"]
    one = _merge(_win(h30, "core_iso_half"))
    assert _read(h30, _parts(h30, **one))["outcome"] == "NONE"  # in force alone: no pass
    other_only = _merge(_win(h30, "full_iso_half"), _win(h30, "half_iso_half"))
    assert _read(h30, _parts(h30, **other_only))["outcome"] == "NONE"  # not on the form in force
    assert _read(h30, _parts(h30, **both), form="half_iso_half")["outcome"] == "NONE"


def test_the_statements_at_their_edges(h30):
    f = "core_iso_half"
    n = lambda a: h30.name(f, a)  # noqa: E731

    def st(**kw):
        return _read(h30, _parts(h30, **_merge(_win(h30, f, **kw), _win(h30, "full_iso_half"))))["statements"]["C"]
    assert st(gain=1.0)["S1"]  # exactly +1.0 on the fold, +1.0 on each transition
    assert not st(gain=0.9999)["S1"]
    nine = np.array([1.2] * 9 + [-0.1])  # mean 1.07, 9 of 10 above 0
    eight = np.array([1.5] * 8 + [-0.1] * 2)
    base = _merge(_win(h30, f), _win(h30, "full_iso_half"))
    for arr, want in ((nine, True), (eight, False)):
        b = json.loads(json.dumps(base, default=lambda x: x.tolist()))
        b["score"][n("A")] = arr
        assert _read(h30, _parts(h30, **b))["statements"]["C"]["S1"] is want
    # +0.5 on each transition: the fold mean can pass while one transition is short
    p = _parts(h30, **base)
    for part in p["T2h/<=E8.75->E10.5"]:
        for r in part["rows"]:
            for q in h30.TRUTHS:
                r["scores"][q][n("A")] = 50.4999
    for part in p["T2h/<=E9.5->E10.5"]:
        for r in part["rows"]:
            for q in h30.TRUTHS:
                r["scores"][q][n("A")] = 52.0
    s = _read(h30, p)["statements"]["C"]
    assert not s["S1"]
    # S2: A2 minus each placebo at least +0.5
    assert st(a2=1.0, plc=0.5)["S2"] and not st(a2=1.0, plc=0.5001)["S2"]
    # S3: A2 not behind the control on either transition
    assert st(a2=0.0, plc=-0.5)["S3"] and not st(a2=-0.0001, plc=-0.6)["S3"]
    # S4: DE terms above zero, and half the margin over the placebos in expression_change plus cell_state
    assert not st(de=0.0)["S4"] and not st(de=-0.01)["S4"]
    # margin +1.0 board point; EC + CS each 0.25 x 100 x 0.01 = 0.25 points: exactly half
    assert st(ec=0.01)["S4"] and not st(ec=0.0099)["S4"]
    one_world = _merge(({n("A"): {"own": np.full(10, 1.2), "CLEAN": np.full(10, 0.9)}}, {}, {}))
    assert not _read(h30, _parts(h30, **one_world))["statements"]["C"]["S1"]


def test_each_outcome_and_its_follows_line(h30):
    f = "core_iso_half"
    placebo = _merge(_win(h30, f, plc=0.8))
    assert _read(h30, _parts(h30, **placebo))["outcome"] == "PLACEBO"
    off = _merge(_win(h30, f, de=-0.01))
    r = _read(h30, _parts(h30, **off))
    assert r["outcome"] == "OFF_DE" and r["follows"] == h30.FOLLOWS["OFF_DE"]
    assert _read(h30, _parts(h30, **_merge(_win(h30, f, gain=0.5))))["outcome"] == "NONE"
    void = _read(h30, h27p=_h27(h30, off=0.2))
    assert void["outcome"] == "VOID" and void["follows"].startswith("VOID: V4: ") and h30.text(void) == \
        void["follows"] + "\n"
    import itertools
    for bits in itertools.product([False, True], repeat=12):
        st = {}
        for i, fm in enumerate(h30.FORMS):
            s1, s2, s3, s4 = bits[4 * i:4 * i + 4]
            st[fm] = {"S1": s1, "S2": s2, "S3": s3, "S4": s4, "all": s1 and s2 and s3 and s4}
        o = h30.outcome(st, f, False)
        assert o in h30.FOLLOWS and o != "VOID"
        assert (o == "PASS") == (st[f]["all"] and any(st[x]["all"] for x in h30.FORMS if x != f))
        assert (o == "PLACEBO") == (st[f]["S1"] and not st[f]["S2"])
    assert h30.outcome(st, f, True) == "VOID"


def test_the_form_in_force_is_read_from_the_row(h30):
    assert h30.form_in_force("| H30 (x) | heart | ... named in this row before the read: C. Changes no file |") == \
        "core_iso_half"
    assert h30.form_in_force("named in this row before the read: full_iso_half.") == "full_iso_half"
    assert h30.form_in_force("named in this row before the read: <form>.") is None
    assert h30.form_in_force("named in this row before the read: A; named in this row before the read: B") is None
    assert h30.form_in_force(None) is None
    r = _read(h30, form=None)
    assert r["void"] and "names no form in force" in r["follows"]


def test_void_on_missing_repeated_mixed_or_dirty_parts(h30):
    p = _parts(h30)
    del p[h30.TASK_KEYS[1]]
    assert any(w.startswith("V1: missing transition") for w in _read(h30, p)["why"])
    p = _parts(h30)
    p[h30.TASK_KEYS[0]].append(_part(h30, h30.TASK_KEYS[0], [3]))
    assert any("more than one part" in w for w in _read(h30, p)["why"])
    p = _parts(h30)
    del p[h30.TASK_KEYS[0]][0]["rows"][0]["scores"]["thin_0"][h30.name("half_iso_half", "P3")]
    assert any("missing from scores" in w for w in _read(h30, p)["why"])
    p = _parts(h30)
    p[h30.TASK_KEYS[0]][1]["provenance"]["commit"] = "b" * 40
    assert any("different commits" in w for w in _read(h30, p)["why"])
    assert any("dirty tree" in w for w in _read(h30, _parts(h30, clean=False))["why"])
    assert any("run commit" in w for w in _read(h30, commit="a" * 12)["why"])  # a prefix is not the commit
    p = _parts(h30)
    p["T2h/E8.25_late->E8.5"] = [_part(h30, "T2h/<=E9.5->E10.5", SEEDS)]
    p["T2h/E8.25_late->E8.5"][0]["task"] = "T2h/E8.25_late->E8.5"
    assert any("other transitions" in w for w in _read(h30, p)["why"])


def test_the_reproduction_gate_per_form_at_its_edges(h30):
    assert not _read(h30, h27p=_h27(h30, off=0.0999))["void"]
    assert _read(h30, h27p=_h27(h30, off=0.1001))["void"]
    h = _h27(h30)
    h[h30.TASK_KEYS[0]][0]["rows"][0]["scores"]["own"][h30.h27.arm_name("half_iso_half", "f")] = 50.5001
    r = _read(h30, h27p=h)
    assert r["void"] and any(w.startswith("V4: B ") and "a seed is" in w for w in r["why"])
    for edit in (lambda p: p["truth_info"]["thin_0"].update(cells=1),
                 lambda p: p["pred_info"]["half_iso_half"].update(core_share=0.91),
                 lambda p: p["forecast_targets"].update({"full_iso_half": 100.0 * (1 + 1e-9)}),
                 lambda p: p["v7"].update(ok=False)):
        p = _parts(h30)
        edit(p[h30.TASK_KEYS[1]][0])
        r = _read(h30, p)
        assert r["void"] and any(w.startswith("V3:") for w in r["why"]), r["why"]


def test_void_on_every_identity_depth_and_table_failure(h30):
    for k in ("chord0_counts_equal_control", "same_draw", "pre_size_coords_equal", "post_size_coords_equal"):
        r = _read(h30, _parts(h30, checks={"B": {k: False}}))
        assert r["void"] and any(w.startswith("V5:") for w in r["why"]), k
    r = _read(h30, _parts(h30, checks={"A": {"placebo_fixed_points": {"P2": {"CM": 1}}}}))
    assert any("own factor" in w for w in r["why"])
    # CHORD0's X against the control's: void above 1e-5 (a normalisation offset), or when the record is missing
    for v, bad in ((1e-5, False), (-1e-5, False), (1.0001e-5, True), (float("nan"), True), (None, True),
                   (True, True)):
        r = _read(h30, _parts(h30, checks={"C": {"chord0_x_max_abs_diff_control": v}}))
        assert any("CHORD0's X differs" in w for w in r["why"]) is bad, v
    p = _parts(h30)
    p[h30.TASK_KEYS[0]][0]["rows"][0]["n_cells"][h30.name("core_iso_half", "A2")] = 349
    assert any("cell numbers differ" in w for w in _read(h30, p)["why"])
    p = _parts(h30)
    p[h30.TASK_KEYS[0]][0]["rows"][0]["skills"]["radial_0"][h30.name("full_iso_half", "P1")]["d2_shape"] += 1.1e-9
    assert any("d2_shape is not the control's" in w for w in _read(h30, p)["why"])
    for v, bad in ((0.0101, True), (-0.0101, True), (0.01, False)):
        p = _parts(h30)
        p[h30.TASK_KEYS[1]][1]["rows"][2]["records"][h30.name("core_iso_half", "A")]["chord"]["mean_total_change"] = v
        assert any("above 1%" in w for w in _read(h30, p)["why"]) is bad
    p = _parts(h30)
    p[h30.TASK_KEYS[1]][1]["rows"][2]["records"][h30.name("core_iso_half", "P1")]["chord"][
        "unmoved_counts_equal"] = False
    assert any("unmoved cell" in w for w in _read(h30, p)["why"])
    # the table: its checksum, the stages its fits read, the embryo list, the Feng stages, the taus
    assert any("not the committed one" in w for w in _read(h30, sha={"e95.json": "x", "e875.json": "t" * 64})["why"])
    edits = (("merfish", ["T2/E9.5.h5ad"]), ("merfish", ["T2/E10.5.h5ad"]), ("jax", [13, 14, 15, 54, 55, 56]),
             ("jax", [13, 14, 54, 55]), ("feng", ["E9.5", "E14.5"]))
    for field, val in edits:
        p = _parts(h30)
        p["T2h/<=E8.75->E10.5"][0]["table"]["stages_read"][field] = val
        assert any(w.startswith("V6:") for w in _read(h30, p)["why"]), (field, val)
    p = _parts(h30)
    p["T2h/<=E9.5->E10.5"][0]["table"]["classifiers"]["feng"]["trained_on"] = ["T2/E9.5.h5ad", "T2/E10.5.h5ad"]
    assert any("a fit read T2/E10.5.h5ad" in w for w in _read(h30, p)["why"])
    p = _parts(h30)
    p["T2h/<=E8.75->E10.5"][0]["table"]["classifiers"]["jax"]["trained_on"] = ["T2/E9.5.h5ad"]
    assert any("a fit read T2/E9.5.h5ad" in w for w in _read(h30, p)["why"])
    p = _parts(h30)
    p["T2h/<=E9.5->E10.5"][0]["table"]["tau_practice"] = 0.4
    assert any("taus" in w for w in _read(h30, p)["why"])
    p = _parts(h30)
    p["T2h/<=E9.5->E10.5"][0]["window_guard"] = ""
    assert any("bypass" in w for w in _read(h30, p)["why"])
    p = _parts(h30)
    del p[h30.TASK_KEYS[0]][0]["truth_info"]
    r = _read(h30, p)
    assert r["void"] and h30.text(r).startswith("VOID: ")


def test_folds_merge_the_two_backtests_seed_by_seed(h30):
    nm = h30.name("core_iso_half", "A")
    p = {h30.TASK_KEYS[0]: [_part(h30, h30.TASK_KEYS[0], SEEDS, score={nm: np.arange(10.0)})],
         h30.TASK_KEYS[1]: [_part(h30, h30.TASK_KEYS[1], SEEDS, score={nm: -np.arange(10.0)})]}
    res = _read(h30, p)
    lab = f"{nm} - {h30.name('core_iso_half', 'control')}"
    f = res["tables"]["own"]["fold E10.5"][lab]["score"]
    assert f["mean"] == 0.0 and f["lo"] == 0.0 and f["hi"] == 0.0
    assert res["tables"]["own"][h30.TASK_KEYS[0]][lab]["score"]["ahead"] == 9


# ------------------------------------------------------------------------------------------------ the script
def _prints(fn) -> list[str]:
    tree = ast.parse(inspect.getsource(fn).lstrip())
    return [ast.unparse(n) for n in ast.walk(tree) if isinstance(n, ast.Call) and getattr(n.func, "id", "") == "print"]


def test_every_step_prints_step_lines_only(h30):
    assert sorted(_prints(h30.run_transition)) == sorted([
        'print(f"{key} seed {s} ({row[\'seconds\']} s) done", flush=True)', "print(f'wrote {path}', flush=True)"])
    assert _prints(h30.run_gate) == ["print(f'gate {outcome}', flush=True)"]
    assert sorted(_prints(h30.run_tables)) == sorted([
        'print(f"tables {TRANSITIONS[key][\'short\']} written", flush=True)',
        'print(f"wrote {out / \'tables_sha256.json\'}", flush=True)'])
    assert sorted(_prints(h30.cmd_read)) == sorted(["print(f\"OUTCOME: {res['outcome']}\")",
                                                    "print(res['follows'])"])
    for fn in (h30.build_arms, h30.chord_expression, h30.build_table, h30.gate_lineage, h30.cmd_run, h30.cmd_gate,
               h30.cmd_tables, h30.read_jax, h30.read_feng):
        assert _prints(fn) == [], fn


def test_the_script_names_no_test_stage_or_board_file(h30):
    src = SCRIPT.read_text()
    for bad in ("results/p3", "size_heart_E12", "submissions", "leaderboard", "registry", "experiments/builds",
                "boards/", "upload", "12.5", "manual", "celltype", "cell_type", "lpm"):
        assert bad not in src, bad
    assert h30.JAX_CELL_COLUMNS == ("embryo_id", "day") and h30.FENG_OBS_COLUMNS == ("stage",)


def _refusal_env(h30, monkeypatch, tmp_path):
    def boom(*a, **k):
        raise AssertionError("something loaded")
    for mod, attr in ((h30.payoff, "heart_task"), (h30, "load_boards"), (h30, "read_jax"), (h30, "read_feng"),
                      (h30, "load_stage")):
        monkeypatch.setattr(mod, attr, boom)
    monkeypatch.setattr(h30, "tries_rows", lambda path=h30.TRIES: {"H27", "H28"})
    monkeypatch.setattr(h30, "register_names_row", lambda *a, **k: [])
    monkeypatch.setattr(h30.provenance, "git_state", lambda repo, ignore=(): ("c" * 40, []))


def test_each_step_refuses_before_anything_loads(h30, monkeypatch, tmp_path):
    _refusal_env(h30, monkeypatch, tmp_path)
    tables = tmp_path / "tables"
    calls = (lambda: h30.run_transition(h30.TASK_KEYS[0], (0,), tmp_path / "o", tmp_path, tables),
             lambda: h30.run_gate(tmp_path / "g", tmp_path, tmp_path),
             lambda: h30.run_tables(tmp_path / "t", tmp_path / "g", tmp_path, tmp_path))
    for call in calls:
        with pytest.raises(SystemExit, match="no row"):
            call()
    # the gate's row N56 is committed with row H30: every step refuses without it (blocker 2)
    monkeypatch.setattr(h30, "tries_rows", lambda path=h30.TRIES: {"H30"})
    for call in calls:
        with pytest.raises(SystemExit, match="gate's row N56"):
            call()
    monkeypatch.setattr(h30, "tries_rows", lambda path=h30.TRIES: {"N56", "H27"})
    for call in calls:
        with pytest.raises(SystemExit, match="no row"):
            call()
    # row H30 must name the form in force before anything runs (blocker 1: TRIES.md cannot change before the read)
    monkeypatch.setattr(h30, "tries_rows", lambda path=h30.TRIES: {"H30", "N56"})
    for txt in (None, "| H30 | ... named in this row before the read: <form>. |"):
        monkeypatch.setattr(h30, "row_text", lambda *a, t=txt, **k: t)
        for call in calls:
            with pytest.raises(SystemExit, match="form in force"):
                call()
    monkeypatch.setattr(h30, "row_text", lambda *a, **k: ROW_C)
    for state in (("c" * 40, ["src/x.py"]), ("unknown", None)):
        monkeypatch.setattr(h30.provenance, "git_state", lambda repo, ignore=(), s=state: s)
        for call in calls:
            with pytest.raises(SystemExit, match="not clean"):
                call()
    monkeypatch.setattr(h30.provenance, "git_state", lambda repo, ignore=(): ("c" * 40, []))
    monkeypatch.setattr(h30, "register_names_row", lambda *a, **k: ["X003"])
    for call in calls:
        with pytest.raises(SystemExit, match="required wording for row H30"):
            call()
    monkeypatch.setattr(h30, "register_names_row", lambda *a, **k: [])
    (tmp_path / "g").mkdir()
    (tmp_path / "g" / "gate.json").write_text("{}")
    with pytest.raises(SystemExit, match="exists"):
        calls[1]()
    with pytest.raises(SystemExit, match="does not say PASS"):  # the tables: the gate of this commit must say PASS
        calls[2]()
    (tmp_path / "g" / "gate.json").write_text(json.dumps(
        {"axis": h30.AXIS + ":gate", "outcome": "PASS", "provenance": {"commit": "d" * 40, "tree_clean": True}}))
    with pytest.raises(SystemExit, match="does not say PASS"):  # another commit's PASS is not enough
        calls[2]()
    # the run: the committed gate must say PASS and the table checksums must be committed
    monkeypatch.setattr(h30, "COMMITTED_GATE", tmp_path / "none.json")
    with pytest.raises(SystemExit, match="committed gate reading does not say PASS"):
        calls[0]()
    gate = tmp_path / "cg.json"
    gate.write_text(json.dumps({"axis": h30.AXIS + ":gate", "outcome": "FAIL"}))
    monkeypatch.setattr(h30, "COMMITTED_GATE", gate)
    with pytest.raises(SystemExit, match="committed gate reading does not say PASS"):
        calls[0]()
    # a committed PASS must come from a clean tree at a strict ancestor of this commit (lesser point 6)
    gate.write_text(json.dumps({"axis": h30.AXIS + ":gate", "outcome": "PASS"}))
    with pytest.raises(SystemExit, match="ancestor"):
        calls[0]()
    seen = []
    monkeypatch.setattr(h30, "is_ancestor", lambda old, new, repo=None: seen.append((old, new)) or old == "b" * 40)
    gate.write_text(json.dumps({"axis": h30.AXIS + ":gate", "outcome": "PASS",
                                "provenance": {"commit": "d" * 40, "tree_clean": True}}))
    with pytest.raises(SystemExit, match="ancestor"):
        calls[0]()
    assert seen[-1] == ("d" * 40, "c" * 40)
    gate.write_text(json.dumps({"axis": h30.AXIS + ":gate", "outcome": "PASS",
                                "provenance": {"commit": "b" * 40, "tree_clean": False}}))
    with pytest.raises(SystemExit, match="ancestor"):
        calls[0]()
    gate.write_text(json.dumps({"axis": h30.AXIS + ":gate", "outcome": "PASS",
                                "provenance": {"commit": "b" * 40, "tree_clean": True}}))
    monkeypatch.setattr(h30, "COMMITTED_SHA", tmp_path / "none_sha.json")
    with pytest.raises(SystemExit, match="sha256 is not committed"):
        calls[0]()
    sha = tmp_path / "sha.json"
    tables.mkdir()
    (tables / "e95.json").write_text('{"task": "T2h/<=E9.5->E10.5"}')
    sha.write_text(json.dumps({"e95.json": "0" * 64}))
    monkeypatch.setattr(h30, "COMMITTED_SHA", sha)
    with pytest.raises(SystemExit, match="not the committed one"):
        calls[0]()
    sha.write_text(json.dumps({"e95.json": h30.sha256_bytes((tables / "e95.json").read_bytes())}))
    with pytest.raises(AssertionError, match="something loaded"):  # past every refusal, the first load
        calls[0]()
    (tmp_path / "o").mkdir()
    (tmp_path / "o" / h30.PART).write_text("{}")
    with pytest.raises(SystemExit, match="exists"):
        calls[0]()


def test_the_read_refuses_without_the_row_a_pass_or_the_form_or_with_a_reading(h30, monkeypatch, tmp_path):
    _refusal_env(h30, monkeypatch, tmp_path)
    monkeypatch.setattr(h30, "load_parts", lambda *a: (_ for _ in ()).throw(AssertionError("something loaded")))
    argv = ["read", "--runs", str(tmp_path), "--out", str(tmp_path / "o"), "--commit", COMMIT]
    with pytest.raises(SystemExit, match="no row"):
        h30.main(argv)
    monkeypatch.setattr(h30, "tries_rows", lambda path=h30.TRIES: {"H30"})
    with pytest.raises(SystemExit, match="gate's row N56"):
        h30.main(argv)
    monkeypatch.setattr(h30, "tries_rows", lambda path=h30.TRIES: {"H30", "N56"})
    monkeypatch.setattr(h30, "COMMITTED_GATE", tmp_path / "none.json")
    with pytest.raises(SystemExit, match="does not say PASS"):
        h30.main(argv)
    gate = tmp_path / "cg.json"
    gate.write_text(json.dumps({"axis": h30.AXIS + ":gate", "outcome": "PASS"}))
    monkeypatch.setattr(h30, "COMMITTED_GATE", gate)
    with pytest.raises(SystemExit, match="ancestor of the run commit"):  # no provenance
        h30.main(argv)
    gate.write_text(json.dumps({"axis": h30.AXIS + ":gate", "outcome": "PASS",
                                "provenance": {"commit": "b" * 40, "tree_clean": True}}))
    monkeypatch.setattr(h30, "is_ancestor", lambda old, new, repo=None: (old, new) == ("b" * 40, COMMIT))
    monkeypatch.setattr(h30, "row_text", lambda *a, **k: "| H30 | ... named in this row before the read: <form>. |")
    with pytest.raises(SystemExit, match="form in force"):
        h30.main(argv)
    monkeypatch.setattr(h30, "row_text", lambda *a, **k: "| H30 | named in this row before the read: B. |")
    monkeypatch.setattr(h30, "COMMITTED_SHA", tmp_path / "sha.json")
    (tmp_path / "sha.json").write_text("{}")
    with pytest.raises(AssertionError, match="something loaded"):
        h30.main(argv)
    (tmp_path / "o").mkdir()
    (tmp_path / "o" / "reading.json").write_text("x")
    with pytest.raises(SystemExit, match="already holds"):
        h30.main(argv)


def test_the_register_must_carry_the_required_wording_in_both_outside_sources(h30, tmp_path):
    p = tmp_path / "EXTERNAL_DATA.md"
    amend, lic = "T2-heart amendment for row H30", "Licence confirmed for row H30"
    p.write_text(f"| X002 | Feng ... T1 only |\n| X003 | Qiu ... T1. {amend}: raw counts. {lic}: CC BY 4.0. |\n")
    assert h30.register_names_row(p) == ["X002"]
    p.write_text(f"| X002 | Feng. {amend}: raw counts |\n| X003 | Qiu. {amend}: raw counts |\n")
    assert h30.register_names_row(p) == ["X003"]  # X003 also needs the licence confirmed
    # merely naming the row, or the phrase inside wording that denies it, is not the amendment
    p.write_text(f"| X002 | Feng ... row H30 |\n| X003 | Qiu ... there is no {amend}. {lic} |\n")
    assert h30.register_names_row(p) == ["X002", "X003"]
    p.write_text(f"| X002 | Feng. {amend}0 |\n| X003 | Qiu. {amend}. Licence not confirmed for row H30 |\n")
    assert h30.register_names_row(p) == ["X002", "X003"]
    p.write_text(f"| X002 | Feng. {amend} (stages E9.5, E14.5) |\n| X003 | Qiu; {amend}; {lic} |\n")
    assert h30.register_names_row(p) == []
    # the committed register: both rows carry their wording since the amendments of 2026-11-22 (placed with rows
    # N56 and H30, before the gate); until then this line read ["X002", "X003"]
    assert h30.register_names_row(h30.REGISTER) == []


def test_the_gate_provenance_must_be_a_strict_ancestor(h30):
    def rev(r):
        return subprocess.run(["git", "rev-parse", r], cwd=REPO_ROOT, capture_output=True, text=True).stdout.strip()
    head, parent = rev("HEAD"), rev("HEAD^")
    assert h30.is_ancestor(parent, head) and not h30.is_ancestor(head, parent)
    assert not h30.is_ancestor(head, head) and not h30.is_ancestor("", head) and not h30.is_ancestor("0" * 40, head)


def test_a_table_is_read_only_with_its_committed_checksum(h30, tmp_path):
    t = tmp_path / "e875.json"
    b = h30.table_bytes({"task": "T2h/<=E8.75->E10.5", "genes": ["a"], "lineages": {"CM": [0.1]}})
    t.write_bytes(b)
    sha = tmp_path / "sha.json"
    sha.write_text(json.dumps({"e875.json": h30.sha256_bytes(b)}))
    got, s = h30.load_table(t, sha)
    assert got["lineages"]["CM"] == [0.1] and s == h30.sha256_bytes(b)
    t.write_bytes(b + b" ")
    with pytest.raises(SystemExit, match="not the committed one"):
        h30.load_table(t, sha)


# ------------------------------------------------------------------------------------------------ the record
def test_the_reproduction_gate_dry_runs_on_row_h27s_committed_parts(h30):
    """The gate's code on committed data: row H27's committed form at f against itself, every form, both backtests,
    every truth of this row (distance 0); a shifted copy fails."""
    h27p = h30.load_h27()
    for key in h30.TASK_KEYS:
        rows = h30.h26.merge_rows(h27p[key])[0]
        assert sorted(r["seed"] for r in rows) == list(SEEDS)
        for form in h30.FORMS:
            arm = h30.h27.arm_name(form, "f")
            mine = json.loads(json.dumps(rows))
            for r in mine:
                for q in h30.TRUTHS:
                    r["scores"][q][h30.name(form, "control")] = r["scores"][q][arm]
            why, big = h30.h28.gate(key, mine, rows, mine=h30.name(form, "control"), theirs=arm, truths=h30.TRUTHS)
            assert why == [] and big == {"mean": 0.0, "seed": 0.0}
            for r in mine:
                r["scores"]["own"][h30.name(form, "control")] += 0.2
            assert h30.h28.gate(key, mine, rows, mine=h30.name(form, "control"), theirs=arm, truths=h30.TRUTHS)[0]


# ------------------------------------------------------------------------------------------------ the launcher
def _run(tmp_path: Path, *argv: str) -> subprocess.CompletedProcess:
    env = {k: v for k, v in os.environ.items() if k not in ("GIT_DIR", "GIT_WORK_TREE", "GIT_INDEX_FILE")}
    env["HOME"] = str(tmp_path)
    return subprocess.run(["bash", str(LAUNCHER), *argv], cwd=REPO_ROOT, env=env, capture_output=True, text=True)


def test_the_launcher_steps_threads_and_no_removal():
    sh = LAUNCHER.read_text()
    code = [line.split("#")[0] for line in sh.splitlines()]
    assert not any(re.search(r"(^|[\s;&|(])(rm|rmdir|unlink|shred|truncate)\s", line) for line in code)
    for msg in re.findall(r"\$\{\d+:\?([^}]*)\}", sh):
        assert "'" not in msg and '"' not in msg, msg
    for need in ("OMP_NUM_THREADS=8", "OPENBLAS_NUM_THREADS=8", "MKL_NUM_THREADS=8", "NUMPY_MADVISE_HUGEPAGE=0",
                 'TMPDIR="$S/tmp"', 'S="$HOME/v/scratch/h30"', "status --porcelain", "python -I -B",
                 "merge-base --is-ancestor", "git diff --name-only", "0-1 2-3 4-5 6-7 8-9",
                 'cmp -s runs/h30/gate/gate.json "$H30R/gate/gate.json"',
                 'cmp -s runs/h30/gate/gate.txt "$H30R/gate/gate.txt"',
                 'cmp -s runs/h30/tables/tables_sha256.json "$H30R/tables_sha256.json"'):
        assert need in sh, need
    # a refusal writes <step>.refused.exit; only a step that ran writes <step>.exit (lesser point 1)
    body = sh.split("export OMP_NUM_THREADS")[0]
    assert "finish " not in body.split("finish() {")[1].split("\n", 1)[1]
    assert re.findall(r"refuse \d+", body) and "$S/$STEP.refused.exit" in sh
    order = dict(re.findall(r"^  (tests|gate|tables|run|read)\) PREV=\"?([a-z]*)\"? ;;$", sh, flags=re.M))
    assert order == {"tests": "", "gate": "tests", "tables": "gate", "run": "tables", "read": "run"}
    listed = re.findall(r"tests/test_[a-z0-9_]+\.py", sh)
    assert "tests/test_h30_heart_late_chord.py" in listed and all((REPO_ROOT / f).exists() for f in listed)
    assert sh.count("finish $") == 5 and sh.count("refuse 9") >= 7 and "$S/$STEP.exit" in sh
    assert "results/h30_heart_late_chord" in sh and "--commit \"$RUN\"" in sh


def test_the_launcher_refusal_paths_under_a_temporary_home(tmp_path):
    s = tmp_path / "v" / "scratch" / "h30"

    def refused(step):
        return (s / f"{step}.refused.exit").read_text().strip()
    out = _run(tmp_path)
    assert out.returncode == 2 and refused("usage") == "2"
    out = _run(tmp_path, "build", "0" * 40)
    assert out.returncode == 97 and refused("build") == "97"
    out = _run(tmp_path, "gate", "0" * 40)  # tests has not run
    assert out.returncode == 96 and refused("gate") == "96"
    out = _run(tmp_path, "run", "0" * 40)  # run needs RUN
    assert out.returncode == 2 and refused("run") == "2"
    out = _run(tmp_path, "read", "0" * 40, "1" * 40)
    assert out.returncode == 96 and refused("read") == "96"
    out = _run(tmp_path, "tests", "0" * 40)  # a wrong commit
    assert out.returncode == 99 and "HEAD is not" in out.stderr and refused("tests") == "99"
    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=REPO_ROOT, capture_output=True, text=True).stdout.strip()
    out = _run(tmp_path, "tests", head[:12])  # the right commit, shortened: never enough
    assert out.returncode == 99
    # a passed tests step at one commit does not open the gate at another
    (s / "tests.exit").write_text("0\n")
    (s / "tests.head").write_text("1" * 40 + "\n")
    assert _run(tmp_path, "gate", "0" * 40).returncode == 96
    (s / "tests.head").write_text("0" * 40 + "\n")
    out = _run(tmp_path, "gate", "0" * 40)
    assert out.returncode == 99 and refused("gate") == "99"
    assert not (s / "gate.head").exists()  # nothing recorded before the checks pass
    # a refused call never overwrites a finished step's exit (the reviewer's stray call: run.exit 0 became 99)
    (s / "gate.exit").write_text("0\n")
    (s / "gate.head").write_text("0" * 40 + "\n")
    (s / "tables.exit").write_text("0\n")
    (s / "tables.head").write_text("0" * 40 + "\n")
    (s / "run.exit").write_text("0\n")
    (s / "run.head").write_text("0" * 40 + "\n")
    (s / "run.run").write_text("1" * 40 + "\n")
    assert _run(tmp_path, "run", "0" * 40, "2" * 40).returncode == 94  # run already exited 0 at this ROW
    assert _run(tmp_path, "tests", "0" * 40).returncode == 94
    out = _run(tmp_path, "read", "0" * 40, "1" * 40)  # the right order, another HEAD
    assert out.returncode == 99 and refused("read") == "99"
    for step in ("tests", "gate", "tables", "run"):
        assert (s / f"{step}.exit").read_text().strip() == "0", step
    assert not (s / "read.exit").exists()
    assert sorted(p.name for p in s.iterdir()) == [
        "build.refused.exit", "gate.exit", "gate.head", "gate.refused.exit", "read.refused.exit", "run.exit",
        "run.head", "run.refused.exit", "run.run", "tables.exit", "tables.head", "tests.exit", "tests.head",
        "tests.refused.exit", "tmp", "usage.refused.exit"]


def test_a_start_at_a_new_commit_never_reads_an_earlier_starts_exit(tmp_path):
    """A scratch folder left by a start at another commit (tests exit 0, gate exit 1): while the new start's step runs,
    its exit file says "running", never the earlier start's code. Traced in a throwaway repository with a stand-in
    `uv` that only records what the scratch folder holds; nothing of this checkout runs."""
    repo = tmp_path / "repo"
    (repo / "experiments" / "harness").mkdir(parents=True)
    (repo / "experiments" / "harness" / "run_h30.sh").write_bytes(LAUNCHER.read_bytes())
    env = {k: v for k, v in os.environ.items() if k not in ("GIT_DIR", "GIT_WORK_TREE", "GIT_INDEX_FILE")}
    env.update(HOME=str(tmp_path), GIT_CEILING_DIRECTORIES=str(tmp_path))
    git = ["git", "-c", "user.name=t", "-c", "user.email=t@t", "-c", "commit.gpgsign=false"]
    for argv in (["init", "-q"], ["add", "-A"], ["commit", "-q", "-m", "x"]):
        subprocess.run([*git, *argv], cwd=repo, env=env, check=True, capture_output=True)
    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=repo, env=env, capture_output=True, text=True,
                          check=True).stdout.strip()
    s = tmp_path / "v" / "scratch" / "h30"
    s.mkdir(parents=True)
    old = "1" * 40
    for step, code in (("tests", "0"), ("gate", "1")):
        (s / f"{step}.exit").write_text(code + "\n")
        (s / f"{step}.head").write_text(old + "\n")
    bin_ = tmp_path / "bin"
    bin_.mkdir()
    (bin_ / "uv").write_text(f'#!/usr/bin/env bash\ncat "{s}/tests.exit" "{s}/tests.head" > "{tmp_path}/seen"\n')
    (bin_ / "uv").chmod(0o755)
    env["PATH"] = f"{bin_}{os.pathsep}{env['PATH']}"
    run = subprocess.run(["bash", "experiments/harness/run_h30.sh", "tests", head], cwd=repo, env=env,
                         capture_output=True, text=True)
    assert run.returncode == 0, run.stderr
    assert (tmp_path / "seen").read_text().split() == ["running", head]  # during the step: not the old 0
    assert (s / "tests.exit").read_text().strip() == "0" and (s / "tests.head").read_text().strip() == head
    assert (s / "gate.exit").read_text().strip() == "1"  # untouched until the gate itself starts
