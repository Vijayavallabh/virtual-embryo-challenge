"""experiments/t1_m1/n2_momentum_size.py (TRIES.md row N2): the factor-size arithmetic, the share computation and the
frozen s2 rule, on synthetic data, and a synthetic end-to-end run that writes no cell count."""

from __future__ import annotations

import importlib.util
import json
import re
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pytest
from scipy import sparse

from vec.core import apply_factors, mean_cp10k
from vec.prep import log1p_cp10k
from vec.t1_scope import with_share

ROOT = Path(__file__).resolve().parents[1]
GENES = [f"g{i}" for i in range(78)] + ["Xist", "Ddx3y"]  # the last two are on the held-flat list
# released labels (configs/t1_chamber_classes.yaml); every number of cells is above 100, so none can pass for a
# number of genes in the outputs
E95 = {"V-CM": 241, "Endocardium": 151, "Proepicardium": 107, "BEC": 113, "Blood": 139, "Hepatocyte": 181}
E105 = {"V-CM": 263, "Endocardium-1": 147, "EndoMT Endocardium": 119, "Early VIC": 131, "Epicardium": 103,
        "BEC": 109, "Blood": 173, "Macrophage": 127, "Hepatocyte": 157}
FENG = ("CM", "Endocardium", "Endothelium", "Epicardium")  # the groups with a Feng direction at 945de7d


@pytest.fixture(scope="module")
def n2():
    spec = importlib.util.spec_from_file_location("t1_m1_n2_momentum_size",
                                                  ROOT / "experiments" / "t1_m1" / "n2_momentum_size.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def _stage(path: Path, sizes: dict, seed: int):
    rng = np.random.default_rng(seed)
    blocks, labels = [], []
    for lab, n in sizes.items():
        mu = rng.gamma(0.4, 1.0, len(GENES))
        C = rng.poisson(1500.0 * mu / mu.sum() * rng.gamma(10.0, 0.1, (n, 1)))
        C[np.arange(n), rng.integers(0, len(GENES), n)] += 1
        blocks.append(C)
        labels += [lab] * n
    C = sparse.csr_matrix(np.vstack(blocks).astype(np.int32))
    a = ad.AnnData(X=log1p_cp10k(C), layers={"counts": C})
    a.var_names, a.obs_names = GENES, [f"c{i}" for i in range(C.shape[0])]
    a.obs["celltype"] = labels
    path.parent.mkdir(parents=True, exist_ok=True)
    a.write_h5ad(path)
    return C, np.array(labels)


@pytest.fixture(scope="module")
def stages(tmp_path_factory):
    data = tmp_path_factory.mktemp("data")
    return data, _stage(data / "T1" / "E9.5_RNA.h5ad", E95, 1), _stage(data / "T1" / "E10.5_RNA.h5ad", E105, 2)


def _run(n2, data, out, monkeypatch):
    monkeypatch.setattr(n2.feng_dirs, "feng_directions",
                        lambda genes, *a, **k: ({g: np.zeros(len(genes)) for g in FENG}, {"groups": {}}))
    n2.main(["--out", str(out), "--data", str(data)])
    return json.loads((out / "n2_momentum_size.json").read_text()), (out / "n2_momentum_size.txt").read_text()


def test_factor_size_counts_genes_by_threshold_among_expressed_genes(n2):
    f = np.array([1.05, 1.15, 0.85, 1.30, 0.70, 2.0, 1.5])
    expressed = np.array([1, 1, 1, 1, 1, 0, 1], dtype=bool)  # the 2.0 gene has no count: not in the denominator
    held = np.array([0, 0, 0, 0, 0, 0, 1], dtype=bool)  # held flat: f = 1 whatever its log factor
    cp10k = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 100.0, 10.0])
    d = n2.factor_size(np.log(f), expressed, held, cp10k)
    assert (d["denominator_genes"], d["held_flat_in_denominator"]) == (6, 1)
    assert {k: d["gt_0.10"][k] for k in ("genes", "up", "down")} == {"genes": 4, "up": 2, "down": 2}
    assert {k: d["gt_0.25"][k] for k in ("genes", "up", "down")} == {"genes": 2, "up": 1, "down": 1}
    assert d["gt_0.10"]["share_of_denominator"] == pytest.approx(4 / 6)
    assert (d["gt_0.10"]["cp10k_share"], d["gt_0.25"]["cp10k_share"]) == pytest.approx((14 / 25, 9 / 25))
    assert d["mean_abs_log_f"] == pytest.approx(np.abs(np.log(f[:5])).sum() / 6)
    s = n2.factor_size(np.log(f), expressed, held, cp10k, scale=1.1)  # f s: 1.155, 1.265, 0.935, 1.43, 0.77; held 1
    assert (s["gt_0.10"]["genes"], s["gt_0.25"]["genes"], s["gt_0.25"]["up"]) == (4, 2, 2)


def test_keep_depth_constant_is_apply_factors_own(n2):
    rng = np.random.default_rng(3)
    C = sparse.csr_matrix(rng.poisson(rng.gamma(0.5, 4.0, (1, 30)), (50, 30)))
    logf, held = rng.normal(0.0, 0.3, 30), np.arange(30) < 3
    logf[held] = 0.0
    _, scale, _ = apply_factors(C, np.exp(logf)[None], np.zeros(50, dtype=int), held, keep_depth=True)
    np.testing.assert_allclose(n2.keep_depth_constant(C, logf, held), scale, rtol=1e-12)


def test_cell_shares_from_the_mix_and_the_s2_rule(n2):
    s1 = {"blood": 0.24, "cm": 0.2, "endocardium": 0.36, "epicardium": 0.1, "vascular_endothelium": 0.1}
    s2 = with_share(s1, "cm", 0.60)  # the heart keeps 0.76; cm takes 0.60 of it, the rest in proportion
    assert s2["blood"] == pytest.approx(0.24) and s2["endocardium"] == pytest.approx(0.36 * 0.4 * 0.76 / 0.56)
    fraction = {"blood": 1.0, "endocardium": 0.25}
    assert n2.cell_shares(s1, fraction) == {"blood": {"category": 0.24, "moved": 0.24},
                                            "endocardium": {"category": 0.36, "moved": 0.09},
                                            "total": {"category": 0.6, "moved": 0.33}}
    assert n2.cell_shares(s2, fraction) == {"blood": {"category": 0.24, "moved": 0.24},
                                            "endocardium": {"category": 0.2, "moved": 0.05},
                                            "total": {"category": 0.44, "moved": 0.29}}


def test_s2_rule_is_read_from_the_frozen_file(n2, tmp_path):
    assert n2.s2_cm_share() == 0.60
    bad = tmp_path / "S2_RULES.md"
    bad.write_text(n2.S2_RULES.read_text().replace("--cm-share 0.60", "--cm-share 0.70"))
    with pytest.raises(SystemExit, match="frozen sha256"):
        n2.s2_cm_share(bad)


def test_synthetic_run_finds_the_momentum_cells_and_writes_no_cell_count(n2, stages, tmp_path, monkeypatch, capsys):
    data, (C95, lab95), (C105, lab105) = stages
    res, txt = _run(n2, data, tmp_path / "out", monkeypatch)
    printed = capsys.readouterr().out
    assert list(res["categories"]) == ["blood", "endocardium"]
    blood, endo = res["categories"]["blood"], res["categories"]["endocardium"]
    assert (blood["moved_labels"], blood["feng_labels"]) == (["Blood", "Macrophage"], [])
    assert (endo["moved_labels"], endo["feng_labels"]) == (["Early VIC"], ["EndoMT Endocardium", "Endocardium-1"])
    assert "factor_whole_category" in endo and "factor_whole_category" not in blood
    assert blood["label_mix"]["E9.5"] == {"Blood": 1.0}
    # shares from the synthetic E10.5 chamber mix and cm at 0.60 of heart cells (s2)
    chamber = sum(n for lab, n in E105.items() if lab != "Hepatocyte")
    blood_n = E105["Blood"] + E105["Macrophage"]
    endo_n = E105["Endocardium-1"] + E105["EndoMT Endocardium"] + E105["Early VIC"]
    heart = chamber - blood_n
    want = {"s1": (blood_n / chamber, endo_n / chamber),
            "s2": (blood_n / chamber, heart / chamber * 0.4 * endo_n / (heart - E105["V-CM"]))}
    for e, (b, en) in want.items():
        assert blood["share_of_cells"][e]["moved"] == pytest.approx(b, abs=0.005 + 1e-9)
        assert endo["share_of_cells"][e]["category"] == pytest.approx(en, abs=0.005 + 1e-9)
        assert endo["share_of_cells"][e]["moved"] == pytest.approx(en * E105["Early VIC"] / endo_n, abs=0.005 + 1e-9)
    # blood's factor from its category means: E10.5 Blood and Macrophage against E9.5 Blood
    b105, b95 = np.isin(lab105, ["Blood", "Macrophage"]), lab95 == "Blood"
    logf = 0.5 * (np.log(mean_cp10k(C105[b105]) + 0.5) - np.log(mean_cp10k(C95[b95]) + 0.5))
    expr, keep = np.asarray(C105[b105].sum(0)).ravel() > 0, np.arange(len(GENES)) < 78
    assert blood["factor"]["gt_0.10"]["genes"] == int((expr & keep & (np.abs(np.exp(logf) - 1) > 0.10)).sum())
    assert blood["factor"]["denominator_genes"] == int(expr.sum())
    assert printed == txt and res["settings"]["s2_cm_share"] == 0.60 and "min_cells" not in res["settings"]
    counts = {*E95.values(), *E105.values(), sum(E95.values()), sum(E105.values()), chamber, endo_n, blood_n, heart,
              blood_n + E105["Early VIC"], sum(n for lab, n in E95.items() if lab != "Hepatocyte")}
    text = (tmp_path / "out" / "n2_momentum_size.json").read_text() + txt
    assert not [n for n in counts if re.search(rf"(?<![\w.]){n}(?![\w.])", text)]


def test_run_stops_when_its_factor_is_not_momentum_shifts(n2, stages, tmp_path, monkeypatch):
    monkeypatch.setitem(n2.SHIFT, "eps", 0.25)  # momentum_shift keeps its default
    with pytest.raises(SystemExit, match="differs from momentum_shift's report"):
        _run(n2, stages[0], tmp_path / "out", monkeypatch)
    assert not (tmp_path / "out").exists()
