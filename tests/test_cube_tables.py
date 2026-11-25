"""Row R72's cube tables and G0 (experiments/t3_decomp/cube_tables.py) on synthetic data and hand-made files: the
stage table on a pinned table's den (the under-30-cells fallback keeps the den), F3's planted gene through the
tables and the axes, the manifest's and the tables' checks, G0's fixed outcome from its checks, and the committed
manifest's own bytes. No competition file is read and nothing is built."""

from __future__ import annotations

import json
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "experiments" / "t3_decomp"))
import cube_tables as ct  # noqa: E402
import sign_cube as sc  # noqa: E402

RULES = {"components": {"factor": {"eps": 1e-5, "clip": [0.25, 4.0], "min_cells": 20},
                        "stage": {"min_cells": 30}}}


def _ctx(seed: int = 0, L: int = 6, G: int = 300):
    rng = np.random.default_rng(seed)
    mu = rng.dirichlet(np.ones(G), size=L)
    mu[:, 5] *= 1e-3  # gene 5 holds a tiny share of every label (< 1e-4 of den)
    mu /= mu.sum(1, keepdims=True)
    p = np.full(L, 1.0 / L)
    return SimpleNamespace(rules=RULES, mu=mu, V=SimpleNamespace(p=p, L=L, labels=[f"t{i}" for i in range(L)]),
                           count=np.full(L, 500)), rng


def test_the_stage_table_is_gene_factor_on_the_pinned_tables_den_with_the_fallback():
    t3d = ct._bs().t3d
    ctx, rng = _ctx()
    mask = np.array([True, True, True, False, True, True])
    mu_x = rng.dirichlet(np.ones(ctx.mu.shape[1]), size=ctx.V.L)
    cnt = np.array([100, 12, 400, 3, 80, 31])  # label 1 under 30 cells: the reference's own mean
    phi, den, low = ct.stage_table(ctx, mu_x, cnt, mask)
    assert list(np.flatnonzero(low)) == [1]
    num = mu_x.copy()
    num[1] = ctx.mu[1]
    want, want_den = t3d.gene_factor(num, ctx.mu, ctx.V.p, mask, 1e-5, (0.25, 4.0))
    assert np.array_equal(phi, want) and np.array_equal(den, want_den)
    _, den_pinned = t3d.gene_factor(rng.dirichlet(np.ones(ctx.mu.shape[1]), size=ctx.V.L), ctx.mu, ctx.V.p, mask,
                                    1e-5, (0.25, 4.0))
    assert np.array_equal(den, den_pinned)  # the den is the pinned table's whatever the numerator
    same, _, _ = ct.stage_table(ctx, ctx.mu, cnt, mask)
    assert np.allclose(same, 1.0)  # WT E8.75 over itself: ln 0


def _planted(eps: float):
    ctx, rng = _ctx(1)
    ctx.rules = {"components": {**RULES["components"], "factor": {**RULES["components"]["factor"], "eps": eps}}}
    g = 5
    mask_m = np.array([True, True, True, True, False, True])
    mask_g = np.array([False, True, False, True, True, True])
    cnt = np.full(ctx.V.L, 200)
    mu95 = rng.dirichlet(np.ones(ctx.mu.shape[1]) * 5, size=ctx.V.L)
    mu80 = rng.dirichlet(np.ones(ctx.mu.shape[1]) * 5, size=ctx.V.L)
    for m in (mu95, mu80):
        m[:, g] = ctx.mu[:, g] * 1.1
    l_m, l_g = rng.normal(0, 0.3, ctx.mu.shape[1]), rng.normal(0, 0.3, ctx.mu.shape[1])
    c = sc.sd(l_m) / sc.sd(l_m + l_g)

    def tables(m95, m80):
        return {f"{p}_{k}": np.log(ct.stage_table(ctx, mm, cnt, mk)[0]) for p, mm in (("S", m95), ("Y", m80))
                for k, mk in (("M", mask_m), ("G", mask_g))}

    def axes(t):
        return sc.cube_axes(l_m, l_g, t["S_M"], t["Y_M"], t["S_G"], t["Y_G"], c)

    return SimpleNamespace(ctx=ctx, g=g, mask_m=mask_m, mask_g=mask_g, cnt=cnt, mu95=mu95, mu80=mu80, c=c,
                           tables=tables, axes=axes)


def test_f3_planted_gene_through_the_tables():
    """F3: on a gene holding < 1e-4 of den, before Q and the clip, to 1e-3: doubled in WT E9.5 only, X_S +c ln2/2,
    X_M -c ln2, X_G -c ln2/2, A +c ln2; halved in WT E8.0 only, X_S +c ln2/2, X_G +c ln2/2, A -c ln2. Every label
    of both masks holds >= 30 stage cells (check2: no under-30 label among the planted gene's labels). The bar needs
    eps far below the gene's den: with rules.yaml's eps (1e-5) a gene holding < 1e-4 of den is moved by
    ln((2 num + eps) / (num + eps)), not ln 2 (the next test), so this one runs at eps 1e-12."""
    w = _planted(1e-12)
    g, c = w.g, w.c
    _, den_m, _ = ct.stage_table(w.ctx, w.mu95, w.cnt, w.mask_m)
    assert den_m[g] / den_m.sum() < 1e-4
    t0 = w.tables(w.mu95, w.mu80)
    assert all(abs(v[g]) < np.log(4) - 0.1 for v in t0.values())  # off the clip
    a0 = w.axes(t0)
    up = w.mu95.copy()
    up[:, g] *= 2.0
    a1 = w.axes(w.tables(up, w.mu80))
    for k, v in {"S": c * np.log(2) / 2, "M": -c * np.log(2), "G": -c * np.log(2) / 2, "A": c * np.log(2)}.items():
        assert a1[k][g] - a0[k][g] == pytest.approx(v, abs=1e-3)
        assert np.abs(np.delete(a1[k] - a0[k], g)).max() <= 1e-3
    down = w.mu80.copy()
    down[:, g] *= 0.5
    a2 = w.axes(w.tables(w.mu95, down))
    for k, v in {"S": c * np.log(2) / 2, "M": 0.0, "G": c * np.log(2) / 2, "A": -c * np.log(2)}.items():
        assert a2[k][g] - a0[k][g] == pytest.approx(v, abs=1e-3)
        assert np.abs(np.delete(a2[k] - a0[k], g)).max() <= 1e-3


def test_the_planted_gene_at_rules_eps_moves_by_the_eps_exact_ratio():
    """With rules.yaml's eps the planted gene moves S_M by ln((2 n + eps) / (n + eps)) less the rescale's change, n
    its masked weighted numerator: the table is gene_factor on the pinned den, pole by pole."""
    w = _planted(1e-5)
    g = w.g
    up = w.mu95.copy()
    up[:, g] *= 2.0
    p = w.ctx.V.p * w.mask_m
    n0, n1 = float(p @ w.mu95[:, g]), float(p @ up[:, g])
    phi0, den, _ = ct.stage_table(w.ctx, w.mu95, w.cnt, w.mask_m)
    phi1, _, _ = ct.stage_table(w.ctx, up, w.cnt, w.mask_m)
    raw0 = (p @ w.mu95 + 1e-5) / (den + 1e-5)
    raw1 = (p @ up + 1e-5) / (den + 1e-5)
    r0, r1 = den @ raw0 / den.sum(), den @ raw1 / den.sum()
    want = np.log((n1 + 1e-5) / (n0 + 1e-5)) - np.log(r1 / r0)
    assert np.log(phi1[g] / phi0[g]) == pytest.approx(want, abs=1e-12)
    assert np.log(phi1[g] / phi0[g]) < np.log(2) - 1e-3  # eps is not small against such a gene's den


def _write_manifest(folder: Path, genes: list[str]) -> Path:
    rng = np.random.default_rng(2)
    recs = {}
    for t in ct.TABLES:
        phi = rng.uniform(0.5, 2, len(genes))
        recs[t] = {"file": f"{t}.csv", "sha256": ct.write_table(folder / f"{t}.csv", genes, phi)}
    import hashlib

    man = {"tables": recs, "genes_sha256": hashlib.sha256("\n".join(genes).encode()).hexdigest()}
    (folder / ct.MANIFEST).write_text(json.dumps(man))
    return folder / ct.MANIFEST


def test_the_manifest_and_its_tables_are_checked_byte_for_byte(tmp_path):
    genes = [f"g{i}" for i in range(20)]
    m = _write_manifest(tmp_path, genes)
    tabs, man, sha = ct.load_tables(m, genes, ct.sha256(m))
    assert sorted(tabs) == sorted(ct.TABLES) and sha == ct.sha256(m)
    with pytest.raises(ct.Refused, match="sha256"):
        ct.load_tables(m, genes, "0" * 64)
    with pytest.raises(ct.Refused, match="genes"):
        ct.load_tables(m, genes[::-1], ct.sha256(m))
    (tmp_path / "S_M.csv").write_text((tmp_path / "S_M.csv").read_text().replace("g3,", "g3x,"))
    with pytest.raises(ct.Refused, match="not the manifest's bytes"):
        ct.load_tables(m, genes, ct.sha256(m))
    with pytest.raises(ct.Refused, match="no manifest"):
        ct.load_manifest(tmp_path / "absent.json")


def test_a_table_off_its_genes_or_values_is_refused(tmp_path):
    genes = ["a", "b", "c"]
    ct.write_table(tmp_path / "t.csv", genes, [1.0, 2.0, 0.5])
    assert np.array_equal(ct.read_table(tmp_path / "t.csv", genes), [1.0, 2.0, 0.5])
    with pytest.raises(ct.Refused):
        ct.read_table(tmp_path / "t.csv", ["a", "c", "b"])
    ct.write_table(tmp_path / "u.csv", genes, [1.0, -2.0, 0.5])
    with pytest.raises(ct.Refused):
        ct.read_table(tmp_path / "u.csv", genes)


def _res(**over):
    ok = {"passes": True}
    res = {"G0a": dict(ok), "G0b": dict(ok), "G0g": dict(ok),
           "G0c": {"axes": {j: dict(ok) for j in sc.AXES}}, "G0e": {"axes": {j: dict(ok) for j in sc.AXES}},
           "G0f": {"axes": {j: dict(ok) for j in sc.AXES}},
           "G0d": {"vertices": {v: dict(ok) for v in sc.VERTICES}}}
    for k, v in over.items():
        res[k] = v
    return res


def test_g0_outcome_follows_its_checks():
    out, line, kept, leaves = ct.g0_outcome(_res())
    assert out == "CONTINUE" and kept == ["S", "M", "G"] and leaves == list(sc.NAMED)
    assert line.startswith("Cube G0: CONTINUE: axes left S, M, G; held at '+' none; vertices G0d refuses none;")
    out, line, _, _ = ct.g0_outcome(_res(G0a={"passes": False}))
    assert out == "STOP" and line == "Cube G0: STOP (G0a: the dens are not pinned): nothing is read or built."
    out, line, _, _ = ct.g0_outcome(_res(G0b={"passes": False}))
    assert out == "STOP" and "G0b" in line
    out, _, kept, leaves = ct.g0_outcome(_res(G0g={"passes": False}))
    assert out == "CONTINUE" and kept == ["M", "G"] and leaves == ["ppn", "pnp", "pnn"]
    r = _res()
    r["G0e"]["axes"]["M"]["passes"] = False
    r["G0g"]["passes"] = False
    out, line, kept, _ = ct.g0_outcome(r)
    assert out == "STOP" and kept == ["G"] and "neither S nor M" in line
    r = _res()
    r["G0d"]["vertices"]["pnn"]["passes"] = False
    r["G0c"]["axes"]["G"]["passes"] = False
    out, line, kept, leaves = ct.g0_outcome(r)
    assert kept == ["S", "M"] and leaves == ["pnp", "npp", "nnp"] and "refuses pnn" in line
    assert "G (G0c)" in line


def test_make_tables_and_g0_refuse_a_dirty_tree_and_never_overwrite(tmp_path, monkeypatch):
    monkeypatch.setattr(ct, "stamp", lambda: {"commit": "c" * 40, "tree_clean": False, "machine": "m", "utc": "t"})
    with pytest.raises(ct.Refused, match="not clean"):
        ct.make_tables(tmp_path, tmp_path / "out")
    a = SimpleNamespace(out=str(tmp_path / "g0"), data=".", run_dir=".", decision=".", manifest=".",
                        manifest_sha256="0" * 64, board="T3:bcatenin", threads=1)
    assert ct.cmd_g0(a) == ct.EXIT_REFUSED
    monkeypatch.setattr(ct, "stamp", lambda: {"commit": "c" * 40, "tree_clean": True, "machine": "m", "utc": "t"})
    (tmp_path / "out").mkdir()
    (tmp_path / "out" / ct.MANIFEST).write_text("{}")
    with pytest.raises(ct.Refused, match="never overwritten"):
        ct.make_tables(tmp_path, tmp_path / "out")
    (tmp_path / "g0").mkdir()
    (tmp_path / "g0" / "g0.json").write_text("{}")
    assert ct.cmd_g0(a) == ct.EXIT_REFUSED


def test_the_committed_manifest_holds_the_committed_tables():
    m = ct.DEFAULT_OUT / ct.MANIFEST
    man, sha = ct.load_manifest(m)
    assert man["script"] == ct.SCRIPT and man["provenance"]["tree_clean"] is True
    assert sorted(man["tables"]) == sorted(ct.TABLES) and man["n_genes"] > 0
    assert {t: r["den_of"] for t, r in man["tables"].items()} == {
        "S_M": "transfer_mab21l2_all", "Y_M": "transfer_mab21l2_all", "S_G": "gata4_pooled", "Y_G": "gata4_pooled"}
    assert {t: r["stage"] for t, r in man["tables"].items()} == {"S_M": "E9.5", "Y_M": "E8.0", "S_G": "E9.5",
                                                                 "Y_G": "E8.0"}
    assert set(man["inputs"]) == {"E8.75", "E9.5", "E8.0", "mab21l2", "rep1", "rep2"}
    assert sha == ct.sha256(m)
