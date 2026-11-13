"""Tries K7a to K7c (experiments/harness/TRIES.md; experiments/t3_decomp/k7_wt_revert.py): the wild-type-only factor's
estimator, its tasks and its reading, on synthetic data only (no competition file, no score of a fold)."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from vec import harness
from vec.adoption import SEEDS

T3_DIR = Path(__file__).resolve().parents[1] / "experiments" / "t3_decomp"
if str(T3_DIR) not in sys.path:
    sys.path.insert(0, str(T3_DIR))
import k7_wt_revert as k7  # noqa: E402

G = 12


def _sample(groups: dict[str, np.ndarray], n: int, rng, depth: int = 400) -> ad.AnnData:
    """Cells of each group drawn as multinomial counts around the group's profile."""
    rows, labs = [], []
    for g, p in groups.items():
        rows.append(rng.multinomial(depth, p / p.sum(), size=n))
        labs += [g] * n
    C = np.vstack(rows).astype(np.int32)
    return ad.AnnData(X=sparse.csr_matrix(np.log1p(C)),
                      obs=pd.DataFrame({"group": labs}, index=[f"c{i}" for i in range(len(labs))]),
                      layers={"counts": sparse.csr_matrix(C)})


def test_the_registered_constants():
    assert (k7.EPS, k7.MIN_CELLS, k7.MIN_SAMPLES, k7.EXPONENT, k7.CLIP) == (1e-5, 200, 3, 0.5, (0.25, 4.0))
    assert (k7.PERM_SEED, k7.MIN_GAIN, k7.MIN_AHEAD, k7.STREAM) == (20261113, 0.3, 8, "t3_wt_revert")
    assert k7.ARMS == ("wt_draw", "wt_revert", "wt_revert_grp", "wt_revert_perm")
    assert k7.FULL_PANEL == ("E7.5", "E8.0", "E8.25_late", "E8.5", "E8.75", "E9.5", "E10.5")


def test_tasks_hold_no_knockout_input_and_only_full_panel_stages():
    assert set(k7.TASKS) == {"g1", "g2", "m"} and len(k7.KEYS) == 3
    for name, t in k7.TASKS.items():
        files = set(t.inputs.values())
        assert not files & harness.KNOCKOUTS and t.target in harness.KNOCKOUTS
        assert t.inputs["wt"] == t.reference and len(t.inputs) == len(k7.FULL_PANEL)
        assert not {"T2/E6.75.h5ad", "T2/E7.25.h5ad"} & files
        assert not t.adoption and not t.counted and t.test_board == "T3:bcatenin"
        assert harness.TaskView.of(t).reference_role == "wt"
        if name != "m":  # the sibling replicate is never an input
            assert set(t.forbidden) == {f for f in harness.KNOCKOUTS if "gata4" in f} - {t.target}
    assert k7.TASKS["g1"].reference == k7.TASKS["g2"].reference == "T2/E8.75.h5ad"
    assert k7.TASKS["m"].reference == "T2/E9.5.h5ad"
    assert all(c.applies(t) for c in k7.candidates()[1:] for t in k7.TASKS.values())
    assert not any(c.applies(t) for c in k7.candidates()[1:] for ts in harness.load_tasks().values() for t in ts)


def test_group_means_are_mean_within_cell_proportions(monkeypatch):
    monkeypatch.setattr(k7, "MIN_CELLS", 50)
    rng = np.random.default_rng(0)
    prof = {"A": rng.gamma(2.0, size=G), "B": rng.gamma(2.0, size=G), "Unknown": rng.gamma(2.0, size=G)}
    a, few = _sample(prof, 400, rng), _sample({"C": prof["A"]}, 20, rng)
    C = sparse.vstack([a.layers["counts"], few.layers["counts"]]).tocsr()
    a = ad.AnnData(X=C.astype(np.float32), obs=pd.DataFrame({"group": [*a.obs["group"], *few.obs["group"]]},
                                                           index=[f"c{i}" for i in range(C.shape[0])]),
                   layers={"counts": C})
    means, sizes = k7.group_means(a)
    assert set(means) == {"A", "B"} and sizes == {"A": 400, "B": 400}  # "Unknown" and the 20-cell group are left out
    for g in means:
        assert means[g].sum() == pytest.approx(1.0)
        assert np.allclose(means[g], prof[g] / prof[g].sum(), atol=0.01)


def test_the_trend_is_the_log_linear_fit_read_at_the_reference_day():
    rng = np.random.default_rng(1)
    base, slope = rng.gamma(2.0, size=G) / 20, rng.normal(0, 0.3, size=G)
    days = [7.5, 8.0, 8.5, 9.5, 10.5]
    means = [{"A": base * np.exp(slope * d), "B": base} for d in days]
    means[0].pop("B")
    means[1].pop("B")
    means[2].pop("B")  # B is held by two samples only: no trend for it
    tr = k7.trend_at(means, days, 8.75)
    assert set(tr) == {"A"}
    assert np.allclose(tr["A"], base * np.exp(slope * 8.75), rtol=2e-3, atol=k7.EPS)  # EPS bends the smallest values
    same = k7.trend_at([{"A": base}, {"A": 2 * base}, {"A": 4 * base}], [8.0, 8.0, 8.0], 9.0)
    assert np.allclose(same["A"], 2 * base, rtol=1e-3)  # one day: the geometric mean


def test_the_factor_undoes_a_planted_deviation_of_the_reference():
    rng = np.random.default_rng(2)
    typical = {g: rng.gamma(2.0, size=G) for g in ("A", "B", "C")}
    typical = {g: p / p.sum() for g, p in typical.items()}
    dev = np.exp(rng.normal(0, 0.4, size=G))  # the reference sample's own per-gene profile, shared by its groups
    ref = {g: p * dev / (p * dev).sum() for g, p in typical.items()}
    ref["D"] = typical["A"]  # a group the other samples do not hold
    sizes = {"A": 500, "B": 300, "C": 200, "D": 900}
    phi, per, rec = k7.revert_factors(ref, sizes, typical)
    assert rec["groups"] == ["A", "B", "C"] and set(per) == {"A", "B", "C"}
    assert np.corrcoef(np.log(phi), -np.log(dev))[0, 1] > 0.99
    assert phi.min() >= k7.CLIP[0] and phi.max() <= k7.CLIP[1]
    t3d = harness._t3_head().t3d
    w = np.array([500.0, 300.0, 200.0])
    want, _ = t3d.gene_factor(np.stack([typical[g] for g in "ABC"]), np.stack([ref[g] for g in "ABC"]), w / w.sum(),
                              np.ones(3, bool), k7.EPS, k7.CLIP)
    assert np.array_equal(phi, want)  # the T3 script's own function, weights the reference's group sizes
    for g in "ABC":
        assert np.corrcoef(np.log(per[g]), -np.log(dev))[0, 1] > 0.99
    with pytest.raises(ValueError):
        k7.revert_factors({"D": typical["A"]}, {"D": 900}, typical)


def _result(key: str, scores: dict[str, list[float]], clean: bool = True, commit: str = "abc") -> dict:
    rows = [{"seed": s, "scores": {a: v[i] for a, v in scores.items()},
             "skills": {a: {"de_score": 0.5} for a in scores},
             "t3": {a: {"severity": {"state": "scored" if a == k7.MAIN else "r2_gate"}} for a in scores}}
            for i, s in enumerate(SEEDS)]
    return {"task": {"key": key}, "rows": rows, "provenance": {"commit": commit, "tree_clean": clean},
            "factor_record": {"groups": ["A"], "sd_log_phi": 0.5, "genes_beyond_twofold": 3, "genes_at_clip": 0}}


def _scores(main: float, grouped: float = 50.0, control: float = 50.0) -> dict[str, list[float]]:
    n = len(SEEDS)
    return {harness.FLOOR: [50.0] * n, k7.BASE: [50.2] * n, k7.MAIN: [main] * n, k7.GROUPED: [grouped] * n,
            k7.CONTROL: [control] * n}


def test_the_reading_needs_every_comparison_on_every_task(tmp_path):
    def write(d, **per_task):
        d.mkdir()
        for name, sc in per_task.items():
            (d / f"{name}.json").write_text(json.dumps(_result(k7.TASKS[name].key, sc)))
        return k7.reading(d)

    yes = write(tmp_path / "yes", g1=_scores(55.0), g2=_scores(56.0), m=_scores(54.0))
    assert not yes["void"] and yes["k7a_supported"] and not yes["k7b_replaces_k7a"]
    assert yes["tasks"]["g1"]["paired"][f"{k7.MAIN} - {k7.BASE}"]["mean"] == pytest.approx(4.8)
    assert yes["tasks"]["m"]["severity_states"][k7.MAIN] == {"scored": len(SEEDS)}
    assert "supported" in k7.text(yes)
    one_short = write(tmp_path / "short", g1=_scores(55.0), g2=_scores(56.0), m=_scores(50.4))  # +0.2 over the base
    assert not one_short["k7a_supported"]
    control = write(tmp_path / "control", g1=_scores(55.0, control=54.9), g2=_scores(56.0), m=_scores(54.0))
    assert not control["k7a_supported"]  # the permuted factor does as well on g1
    grouped = write(tmp_path / "grouped", g1=_scores(55.0, 56.0), g2=_scores(56.0, 57.0), m=_scores(54.0, 55.0))
    assert grouped["k7a_supported"] and grouped["k7b_replaces_k7a"]
    missing = write(tmp_path / "missing", g1=_scores(55.0), g2=_scores(56.0))
    assert missing["void"] and any("missing task m" in w for w in missing["why"])
    d = tmp_path / "dirty"
    d.mkdir()
    for name in k7.TASKS:
        (d / f"{name}.json").write_text(json.dumps(_result(k7.TASKS[name].key, _scores(55.0), clean=name != "m")))
    assert k7.reading(d)["void"]
    d = tmp_path / "commits"
    d.mkdir()
    for i, name in enumerate(k7.TASKS):
        (d / f"{name}.json").write_text(json.dumps(_result(k7.TASKS[name].key, _scores(55.0), commit=f"c{i % 2}")))
    assert k7.reading(d)["void"]


def test_ahead_in_eight_seeds_is_required(tmp_path):
    n = len(SEEDS)
    sc = _scores(55.0)
    sc[k7.MAIN] = [80.0] * 7 + [40.0] * (n - 7)  # a large mean gain, ahead in 7 seeds only
    for name in k7.TASKS:
        (tmp_path / f"{name}.json").write_text(json.dumps(_result(k7.TASKS[name].key, sc)))
    res = k7.reading(tmp_path)
    assert not res["k7a_supported"] and res["tasks"]["g1"]["paired"][f"{k7.MAIN} - {k7.BASE}"]["ahead"] == 7
