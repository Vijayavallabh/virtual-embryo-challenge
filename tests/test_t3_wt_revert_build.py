"""The wild-type-only factor in the T3 builder (experiments/t3_decomp/build_submission.py --wt-revert; the captain's
ruling of 2026-11-13, experiments/harness/TRIES.md row R5), on synthetic data only: it is ruled for T3's second upload
alone, it is off unless asked for, its factors are the tried estimator's (k7_wt_revert), and each drawn cell is moved
by its own lineage group's factor through t3_decompose.apply_gene_factor."""

from __future__ import annotations

import inspect
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from vec.boards import Board

T3_DIR = Path(__file__).resolve().parents[1] / "experiments" / "t3_decomp"
if str(T3_DIR) not in sys.path:
    sys.path.insert(0, str(T3_DIR))
import build_submission as bs  # noqa: E402
import k7_wt_revert as k7  # noqa: E402

G = 14
GENES = [f"g{i}" for i in range(G)]


def test_the_option_is_off_by_default_and_ruled_for_the_second_upload_only(tmp_path):
    assert inspect.signature(bs.build).parameters["wt_revert"].default.default is False
    assert inspect.signature(bs.build_submission).parameters["wt_revert"].default is None
    assert bs.WT_REVERT_SLOT == ("refused", "s1")
    b = Board(key="T3:synthetic", task="T3", genes=GENES, needs_coords=True, max_cells=3000)
    for branch, sub in (("allowed", "s1"), ("allowed", "s2"), ("refused", "s2"), ("barred", "s1"), ("barred", "s2")):
        with pytest.raises(bs.Refused, match="second upload only"):  # before any file is read: none of these exists
            bs.build_submission(tmp_path / "d.json", branch, sub, tmp_path / "E8.75.h5ad", b, 100,
                                tmp_path / "out.h5ad", wt_revert={})


def test_the_series_is_every_other_full_panel_stage(tmp_path):
    want = sorted(n for n in k7.FULL_PANEL if n != "E8.75")
    files = {n: tmp_path / f"{n}.h5ad" for n in want}
    for f in files.values():
        f.write_bytes(b"x")
    assert bs.wt_revert_series(tmp_path / "E8.75.h5ad", files) == files
    with pytest.raises(bs.Refused, match="exactly the stages"):
        bs.wt_revert_series(tmp_path / "E8.75.h5ad", {n: f for n, f in files.items() if n != "E9.5"})
    with pytest.raises(bs.Refused, match="exactly the stages"):
        bs.wt_revert_series(tmp_path / "E8.75.h5ad", {**files, "E8.75": tmp_path / "E8.75.h5ad"})
    with pytest.raises(bs.Refused, match="not one of the stages"):
        bs.wt_revert_series(tmp_path / "E6.75.h5ad", files)
    files["E9.5"].unlink()
    with pytest.raises(bs.Refused, match="not found"):
        bs.wt_revert_series(tmp_path / "E8.75.h5ad", files)


def _stage(profiles: dict[str, np.ndarray], n: int, seed: int) -> ad.AnnData:
    rng = np.random.default_rng(seed)
    rows, labs = [], []
    for g, p in profiles.items():
        rows.append(rng.multinomial(300, p / p.sum(), size=n))
        labs += [g] * n
    C = sparse.csr_matrix(np.vstack(rows).astype(np.int32))
    obs = pd.DataFrame({"group": labs}, index=[f"c{i}" for i in range(len(labs))])
    return ad.AnnData(X=C.astype(np.float32), obs=obs, var=pd.DataFrame(index=GENES), layers={"counts": C})


def test_the_factors_are_the_tried_estimators(tmp_path, monkeypatch):
    monkeypatch.setattr(k7, "MIN_CELLS", 40)
    rng = np.random.default_rng(3)
    typical = {g: rng.gamma(2.0, size=G) for g in ("A", "B", "Unknown")}
    names = sorted(n for n in k7.FULL_PANEL if n != "E8.75")
    stages = {n: _stage(typical, 60, i) for i, n in enumerate(names)}
    dev = np.exp(rng.normal(0, 0.4, size=G))
    stages["E8.75"] = _stage({g: p * dev for g, p in typical.items()}, 60, 99)
    files = {n: tmp_path / f"{n}.h5ad" for n in stages}
    for f in files.values():
        f.write_bytes(f.name.encode())
    loads = []

    def load(path, genes=None, drop_duplicates=False, **kw):
        assert drop_duplicates and list(genes) == GENES
        loads.append(Path(path).stem)
        return stages[Path(path).stem]

    import vec.prep

    monkeypatch.setattr(vec.prep, "load_stage", load)
    series = {n: files[n] for n in names}
    phi, per, group, rec = bs.wt_revert_factors(files["E8.75"], series, GENES)
    assert loads == ["E8.75", *names]  # the samples in sorted name order, as the tries' roles were
    from vec.harness import STAGE_DAYS

    ref_means, ref_sizes = k7.group_means(stages["E8.75"])
    trend = k7.trend_at([k7.group_means(stages[n])[0] for n in names], [STAGE_DAYS[f"T2/{n}.h5ad"] for n in names],
                        STAGE_DAYS["T2/E8.75.h5ad"])
    want_phi, want_per, _ = k7.revert_factors(ref_means, ref_sizes, trend)
    assert np.array_equal(phi, want_phi) and set(per) == set(want_per) == {"A", "B"}
    assert all(np.array_equal(per[g], want_per[g]) for g in per)
    assert np.corrcoef(np.log(phi), -np.log(dev))[0, 1] > 0.9  # it undoes the reference's own deviation
    assert list(group) == list(stages["E8.75"].obs["group"])
    assert rec["exponent"] == k7.EXPONENT == 0.5 and rec["rule"] == k7.WT_REVERT_RULE
    assert set(rec["series"]) == set(names) and all(len(v["sha256"]) == 64 for v in rec["series"].values())
    assert "sizes" not in rec and not any("cells" in k for k in rec)  # a build record holds no cell number
    C = bs.t3d.dense_counts(stages["E8.75"])
    bs.wt_revert_factors(files["E8.75"], series, GENES, C)  # the pool it was fitted for: accepted
    with pytest.raises(bs.Refused, match="not the pool"):
        bs.wt_revert_factors(files["E8.75"], series, GENES, C[:-1])

    def lacking(path, **kw):
        raise KeyError(f"{path}: 1 panel genes missing")

    monkeypatch.setattr(vec.prep, "load_stage", lacking)
    with pytest.raises(bs.Refused, match="panel genes missing"):
        bs.wt_revert_factors(files["E8.75"], series, GENES)


def test_the_form_must_be_the_one_the_tries_read():
    import k8_s2_form as k8

    src = (T3_DIR / "build_submission.py").read_text()
    assert 'tried, comp = k8.spec()["composition"], spec.get("composition", {})' in src
    assert 'for k in ("kind", "depleted_states", "log2_magnitude")' in src and 'spec.get("coverage") or any(' in src
    assert set(k8.spec()["composition"]) == {"kind", "depleted_states", "log2_magnitude"}


def _reference(rng, n=90):
    t3d = bs.t3d
    C = rng.poisson(3.0, size=(n, G)).astype(np.int32)
    C[C.sum(1) == 0, 0] = 1
    lab = rng.integers(0, 3, n)
    P = (C / C.sum(1, keepdims=True)).astype(np.float32)
    pool = t3d.Pool(C, lab, P, 3)
    ref = bs.Reference(preset="gata4", rules={}, genes=GENES, V=None, pool=pool, coords=rng.normal(size=(n, 3)))
    return ref


def test_each_cell_is_moved_by_its_own_groups_factor():
    rng = np.random.default_rng(5)
    ref = _reference(rng)
    n = ref.pool.C.shape[0]
    group = np.array(["A", "B", "Unknown"])[rng.integers(0, 3, n)]
    row = rng.integers(0, n, 150)
    built = bs.Built(X=bs.t3d.lognorm(ref.pool.C[row]), lab=ref.pool.lab[row], source=np.full(150, "wt", dtype=object),
                     row=row)
    phi = np.exp(rng.normal(0, 0.5, G))
    per = {"A": np.exp(rng.normal(0, 0.5, G)), "B": np.exp(rng.normal(0, 0.5, G))}
    out, how = bs.apply_wt_revert(ref, built, phi, per, group, np.random.default_rng(11))
    assert how == {"factor_by_group": {"A": "its own factor", "B": "its own factor", "Unknown": "the pooled factor"}}
    assert np.array_equal(out.row, row) and np.array_equal(out.lab, built.lab) and np.all(out.source == "wt")
    want = ref.pool.C[row].copy()
    check = np.random.default_rng(11)  # the same generator, the groups in sorted order
    for g in ("A", "B", "Unknown"):
        rows = np.flatnonzero(group[row] == g)
        want[rows] = bs.t3d.apply_gene_factor(want[rows], per.get(g, phi) ** k7.EXPONENT, ref.pool.P[row][rows], check)
    assert np.array_equal(out.X, bs.t3d.lognorm(want)) and not np.array_equal(out.X, built.X)
    again, _ = bs.apply_wt_revert(ref, built, phi, per, group, np.random.default_rng(11))
    assert np.array_equal(again.X, out.X)  # a seed's file is reproducible
    mixed = bs.Built(X=built.X, lab=built.lab, source=np.array(["wt"] * 149 + ["stage:E9.5"], dtype=object), row=row)
    with pytest.raises(bs.Refused, match="reference pool"):
        bs.apply_wt_revert(ref, mixed, phi, per, group, np.random.default_rng(11))
    with pytest.raises(bs.Refused, match="line up"):
        bs.apply_wt_revert(ref, built, phi, per, group[:-1], np.random.default_rng(11))
