"""Try K18 (experiments/harness/TRIES.md; experiments/t3_decomp/k18_built_files.py): T3's two built last-upload files
read against each other in two practice worlds. Synthetic data, hand-made score tables and committed records only (no
competition file, no score of a fold)."""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path
from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from vec import harness
from vec.adoption import SEEDS

REPO = Path(__file__).resolve().parents[1]
T3_DIR = REPO / "experiments" / "t3_decomp"
if str(T3_DIR) not in sys.path:
    sys.path.insert(0, str(T3_DIR))
import k7_wt_revert as k7  # noqa: E402
import k8_s2_form as k8  # noqa: E402
import k18_built_files as k18  # noqa: E402

G, N_POOL, N_DRAW = 14, 90, 150
GENES = [f"g{i}" for i in range(G)]


def test_arms_worlds_tasks_and_constants():
    assert k18.ARMS == ("refused_form", "barred_form", "resid_revert_k7", "resid_perm_k7", "resid_revert_up",
                        "resid_perm_up")
    assert k18.RULE_ARMS == (k18.REFUSED, k18.BARRED) and k18.WORLDS == ("plain", "dm")
    assert k18.STREAM == k8.STREAM == "t3_s2_form"
    assert k18.TASKS == k8.TASKS and set(k18.TASKS) == {"g1", "g2"}
    assert (k18.MIN_GAIN, k18.MIN_AHEAD, k18.MIN_SCORED, k18.GUARD_MEAN, k18.GUARD_SEED) == (0.3, 8, 9, 0.1, 0.5)
    cands = k18.candidates()
    assert [c.name for c in cands] == [harness.FLOOR, *k18.ARMS]
    assert all(c.applies(t) for c in cands[1:] for t in k18.TASKS.values())
    assert not any(c.applies(k7.TASKS["m"]) for c in cands[1:])  # the central scenario is WT E8.75's
    assert k18.DM_TARGETS == {"g1": "rep1", "g2": "rep2"}


def test_the_barred_arm_is_the_builders_barred_s1():
    """The spec comes from the committed decision file through the builder's own resolve_spec: the Mab21l2 factor at
    exponent 1, the upload's central composition, no extra step."""
    bs = harness._t3_head()
    spec, dec, rules = k18.barred_spec()
    want, _ = bs.resolve_spec(dec, rules, "barred", "s1", bs.check_hashes(dec, False))
    assert spec == want
    assert spec["factor"] == {"source": "mab21l2", "exponent": 1.0}
    assert {k: spec["composition"][k] for k in ("kind", "depleted_states", "log2_magnitude")} == {
        k: k8.spec()["composition"][k] for k in ("kind", "depleted_states", "log2_magnitude")}
    assert not any(spec.get(k) for k in bs.EXTRAS)
    assert bs.factor_name(spec["factor"]["source"], dec["preset"]) == "transfer_mab21l2_all"
    files = dec["submissions"]["parameters"]["factor_files"]
    assert "transfer_mab21l2_all" in files and len(files["transfer_mab21l2_all"]["sha256"]) == 64


def test_the_residual_is_orthogonal_to_the_pooled_factor_and_vanishes_on_a_power_of_it():
    rng = np.random.default_rng(1)
    pooled = np.exp(rng.normal(0, 0.2, G))
    per = {"A": np.exp(rng.normal(0, 0.2, G)), "B": pooled ** 0.7, "C": np.exp(rng.normal(0, 0.2, G))}
    res = k18.residual_factors(per, pooled)
    p = np.log(pooled)
    for g in ("A", "C"):
        assert abs(np.log(res[g]) @ p) < 1e-10
        assert not np.allclose(res[g], per[g])
    assert np.allclose(res["B"], 1.0)
    big = k18.residual_factors({"A": np.full(G, 100.0)}, np.exp(np.r_[1.0, np.zeros(G - 1)]))
    assert big["A"].max() == k7.CLIP[1]  # the clip holds
    assert np.allclose(k18.residual_factors(per, np.ones(G))["A"], per["A"])  # a zero pooled log leaves the factor


def test_gene_order():
    phi = np.arange(1.0, 5.0)
    assert np.array_equal(k18.in_gene_order(phi, ["a", "b", "c", "d"], ["c", "a"]), [3.0, 1.0])
    with pytest.raises(RuntimeError):
        k18.in_gene_order(phi, ["a", "b", "c", "d"], ["a", "z"])


# --- the builds, on a synthetic reference with the builder's draw stubbed ---


@pytest.fixture
def fold(monkeypatch):
    bs = harness._t3_head()
    rng = np.random.default_rng(5)
    C = rng.poisson(3.0, size=(N_POOL, G)).astype(np.int32)
    C[C.sum(1) == 0, 0] = 1
    lab = rng.integers(0, 3, N_POOL)
    P = (C / C.sum(1, keepdims=True)).astype(np.float32)
    V = SimpleNamespace(p=np.array([0.5, 0.3, 0.2]))
    ref = bs.Reference(preset="gata4", rules={}, genes=GENES, V=V, pool=bs.t3d.Pool(C, lab, P, 3),
                       coords=rng.normal(size=(N_POOL, 3)))
    pool = ad.AnnData(X=bs.t3d.lognorm(C).astype(np.float32), var=pd.DataFrame(index=GENES))
    group = np.array(["A", "B", "Unknown"])[rng.integers(0, 3, N_POOL)]
    phi = np.exp(rng.normal(0, 0.5, G))
    per = {"A": np.exp(rng.normal(0, 0.5, G)), "B": np.exp(rng.normal(0, 0.5, G))}
    prep = {"ref": ref, "phi": phi, "per": per, "group": group, "record": {}}
    f_s2 = np.array([0.6, 0.25, 0.15])
    file_genes = GENES[::-1] + ["extra"]  # the factor file's own gene order differs from the panel's
    barred_file = np.exp(rng.normal(0, 0.5, G + 1))
    spec, _, _ = k18.barred_spec()
    log = {"construct": [], "composition": [], "assemble": []}

    def construct(r, n, g, f=None, phi=None, w=1.0, ratio=None, alphas=None):
        assert r is ref and ratio is None and alphas is None
        log["construct"].append({"f": np.asarray(f).copy(), "phi": None if phi is None else np.asarray(phi).copy(),
                                 "w": w})
        rows = g.integers(0, N_POOL, n)
        Cr = C[rows]
        if phi is not None:
            Cr = bs.t3d.apply_gene_factor(Cr, np.asarray(phi) ** w, P[rows], g)
        return bs.Built(X=bs.t3d.lognorm(Cr), lab=lab[rows], source=np.full(n, "wt", dtype=object), row=rows)

    def composition(v, s, params):
        log["composition"].append(s)
        return f_s2, {"states_missing": []}

    def assemble(r, built, panel, g):
        log["assemble"].append({"row": built.row.copy(), "draw": g.random()})
        return built.X, None, {}

    monkeypatch.setattr(k7, "_prepared", lambda ctx: prep)
    monkeypatch.setattr(bs, "construct", construct)
    monkeypatch.setattr(bs, "composition", composition)
    monkeypatch.setattr(bs, "assemble", assemble)
    monkeypatch.setattr(harness, "count_for", lambda board, n: N_DRAW)
    monkeypatch.setattr(harness, "finish", lambda ctx, p, X, coords: X)
    monkeypatch.setattr(k18, "_BARRED", {"phi": barred_file, "genes": file_genes, "spec": spec, "record": {}})
    k18._FOLD.clear()

    def ctx(task: str, seed: int):
        return SimpleNamespace(task=SimpleNamespace(key=k18.TASKS[task].key, reference_role="wt"),
                               inputs={"wt": pool}, board=SimpleNamespace(genes=GENES, needs_coords=False), seed=seed)

    return SimpleNamespace(bs=bs, ctx=ctx, prep=prep, log=log, f_s2=f_s2, V=V,
                           phi_barred=k18.in_gene_order(barred_file, file_genes, GENES))


@pytest.mark.parametrize("seed", [0, 7])
def test_refused_form_is_k8b_value_for_value(fold, seed):
    for task in ("g1", "g2"):
        want = k8._build(fold.ctx(task, seed), k8.GROUPED)
        got = k18._build_plain(fold.ctx(task, seed), k18.REFUSED)
        assert np.array_equal(got, want)
        a, b = fold.log["assemble"][-2:]
        assert a["draw"] == b["draw"] and np.array_equal(a["row"], b["row"])


def test_the_barred_form_takes_the_builders_keywords_and_the_upload_arms_share_their_cells(fold):
    X = {}
    for arm in k18.ARMS:
        fold.log["construct"].clear()
        fold.log["assemble"].clear()
        X[arm] = k18._build_plain(fold.ctx("g2", 3), arm)
        (c,) = fold.log["construct"]
        form, step = k18.ARM_SPEC[arm]
        assert np.array_equal(c["f"], fold.f_s2 if form == "upload" else fold.V.p)
        if step == "barred":  # the Mab21l2 factor in the panel's order, at exponent 1, through prediction_inputs
            assert np.array_equal(c["phi"], fold.phi_barred) and c["w"] == 1.0
        else:
            assert c["phi"] is None
        fold.log[arm] = fold.log["assemble"][0]
    upload = [a for a in k18.ARMS if k18.ARM_SPEC[a][0] == "upload"]
    for a in upload:  # one stream: the same drawn cells and the same assemble stream in every arm of the upload form
        assert np.array_equal(fold.log[a]["row"], fold.log[k18.REFUSED]["row"])
        assert fold.log[a]["draw"] == fold.log[k18.REFUSED]["draw"]
    for a in k18.ARMS:  # every arm can differ from refused_form, and each control from its arm
        if a != k18.REFUSED:
            assert not np.array_equal(X[a], X[k18.REFUSED])
    for c, a in k18.CONTROL_OF.items():
        assert not np.array_equal(X[c], X[a])
    # the composition is asked for twice: the upload's (k8.spec) and the barred spec's, through prediction_inputs
    assert [s["factor"]["source"] for s in fold.log["composition"]] == ["none", "mab21l2"]
    assert fold.log["composition"][0] == k8.spec()


def test_a_depth_matched_row_draws_every_arm_from_the_matched_pool(monkeypatch):
    """dm_row on a synthetic context: the frozen script's draws and orth_diag's thinning run for real; the estimate
    (a depth ratio of 0.5 for every label), the scorer and de_genes are stubbed. Every arm and the check must be built
    from the pool thinned to the ratio, and the check on U3d's stream with the reference's composition."""
    import orth_diag as od
    from common import core_metrics

    from vec import evaluate

    bs = harness._t3_head()
    t3d = bs.t3d
    rng = np.random.default_rng(3)
    L, n_wt, n_ref, n_t = 3, 120, 80, 60
    C = rng.poisson(20.0, size=(n_wt, G)).astype(np.int32)
    lab = rng.integers(0, L, n_wt)
    P = (C / C.sum(1, keepdims=True)).astype(np.float32)
    V = SimpleNamespace(p=np.array([0.5, 0.3, 0.2]), L=L)
    rules = {"protocol": {"working_copy_fraction": 0.5, "split_seed": 0, "prediction_cells_multiplier": 2},
             "components": {"stage": {"min_cells": 5}}}
    T = SimpleNamespace(name="rep1", X_rel=rng.normal(size=(n_t, G)).astype(np.float32),
                        hash_rel=np.arange(n_t), hash_pool=np.arange(n_t + 20))
    ref_C = rng.poisson(20.0, size=(n_ref, G)).astype(np.int32)
    ctx = SimpleNamespace(targets=[T], rules=rules, V=V, ref_X=t3d.lognorm(ref_C), preset_name="gata4", genes=GENES,
                          wt=SimpleNamespace(pool=t3d.Pool(C, lab, P, L)), stages=[])
    group = np.array(["A", "B"])[rng.integers(0, 2, n_wt)]
    per = {"A": np.exp(rng.normal(0, 0.3, G)), "B": np.exp(rng.normal(0, 0.3, G))}
    phi = np.exp(rng.normal(0, 0.3, G))
    spec, _, _ = k18.barred_spec()
    monkeypatch.setattr(bs, "composition", lambda v, s, params: (np.array([0.6, 0.25, 0.15]), {"states_missing": []}))
    refw = bs.Reference(preset="gata4", rules=rules, genes=GENES, V=V, pool=ctx.wt.pool, coords=np.zeros((n_wt, 3)))
    phi_b = np.exp(rng.normal(0, 0.3, G))
    F = k18.fold_factors(bs, refw, per, phi, group, ctx.ref_X.mean(0), phi_b, spec)
    seen = []
    monkeypatch.setattr(t3d, "estimate", lambda c, t, rows, g: SimpleNamespace(ratio=np.full(L, 0.5)))
    monkeypatch.setattr(evaluate, "panel_metrics", lambda *a, **k: {})
    monkeypatch.setattr(core_metrics, "de_genes", lambda a, b: (np.arange(5), np.arange(5, 10), None))
    monkeypatch.setattr(k18, "_scored", lambda X, truth, ref, anchors, up, dn, seed: seen.append(
        np.expm1(np.asarray(X, dtype=np.float64))) or {"score": 50.0, "skills": {}, "severity": {"state": "scored",
                                                                                            "beta": 0.1, "r2": 0.1}})
    real_build = t3d.build
    streams = []

    def build(c, n, g, f, src=None, **kw):
        streams.append((json.dumps(g.bit_generator.state, default=int), np.asarray(f).copy(),
                        c.wt.pool if src is None else src))
        return real_build(c, n, g, f, src=src, **kw)

    monkeypatch.setattr(t3d, "build", build)
    st = {"ctx": ctx, "phi_b": phi_b, "ref_counts": ref_C, "ref_lab": rng.integers(0, L, n_ref),
          "ref_props": (ref_C / ref_C.sum(1, keepdims=True)), "F": F, "threads": 1}
    row = k18.dm_row(st, 0, 4)
    assert row["target"] == "rep1" and row["seed"] == 4 and set(row["scores"]) == set(k18.ARMS)
    assert 0.4 < row["depth_kept"]["pool"] < 0.6 and 0.4 < row["depth_kept"]["reference"] < 0.6
    check_state = json.dumps(t3d.stable_rng(4, od.B_STREAM).bit_generator.state, default=int)
    arm_state = json.dumps(t3d.stable_rng(4, k18.STREAM).bit_generator.state, default=int)
    assert streams[0][0] == check_state and np.array_equal(streams[0][1], V.p)
    assert all(s[0] == arm_state for s in streams[1:]) and len(streams) == 1 + len(k18.ARMS)
    assert {id(s[2]) for s in streams} == {id(streams[0][2])}  # the check and every arm: one pool
    assert streams[0][2] is not ctx.wt.pool and streams[0][2].C.sum() < 0.6 * C.sum()  # thinned, not as released
    assert len(seen) == 1 + len(k18.ARMS)


# --- the guard ---


def test_the_guard_on_the_committed_k8_files_and_the_distances_the_machines_showed():
    for n in k18.TASKS:
        rows = k18.committed_k8(n)
        assert [r["seed"] for r in sorted(rows, key=lambda r: r["seed"])] == list(SEEDS)
        d = k18.distance(sorted(rows, key=lambda r: r["seed"]), k8.GROUPED, rows, k8.GROUPED)
        assert d == {"mean_diff": 0.0, "max_abs_diff": 0.0, "seeds_off": []} and k18.guard_passes(d)
    # K16's base arm (K8b's arm, rerun at cccff68) against K8b, as K17 committed it, on both machines
    k17r = json.loads((T3_DIR / "results" / "k17_revert_dose_reread" / "reading.json").read_text())
    for machine in ("local", "dgx"):
        for n in ("g1", "g2"):
            d = k17r[machine]["base_distance"][n]
            assert k18.guard_passes(d)
            assert abs(d["mean_diff"]) < 0.002 and d["max_abs_diff"] < 0.02  # as row K18 quotes them
    assert k17r["local"]["base_distance"]["g1"]["max_abs_diff"] == 0.0


def test_the_guard_refuses_what_it_should():
    ref = [{"seed": s, "scores": {"x": 50.0 + s}} for s in SEEDS]

    def rows(shift):
        return [{"seed": s, "scores": {"x": 50.0 + s + shift(s)}} for s in SEEDS]

    assert k18.guard_passes(k18.distance(rows(lambda s: 0.09), "x", ref, "x"))
    assert not k18.guard_passes(k18.distance(rows(lambda s: 0.11), "x", ref, "x"))
    assert not k18.guard_passes(k18.distance(rows(lambda s: 0.6 if s == 3 else 0.0), "x", ref, "x"))
    assert k18.guard_passes(k18.distance(rows(lambda s: 0.45 if s == 3 else 0.0), "x", ref, "x"))
    assert not k18.guard_passes(k18.distance(rows(lambda s: 0.0)[:9], "x", ref[:8], "x"))  # a seed missing
    # pure noise: per-seed noise of sd 0.01 (a few cells drawn otherwise) passes; sd 1 (another draw) fails
    rng = np.random.default_rng(0)
    small = [k18.guard_passes(k18.distance(rows(lambda s: rng.normal(0, 0.01)), "x", ref, "x")) for _ in range(200)]
    large = [k18.guard_passes(k18.distance(rows(lambda s: rng.normal(0, 1.0)), "x", ref, "x")) for _ in range(200)]
    assert all(small) and not any(large)


# --- the rule ---


def _rows(scores: dict[str, float], gated: dict[str, int] | None = None, noise=None, seed_shift=0.1) -> list[dict]:
    gated = gated or {}
    out = []
    for s in SEEDS:
        sc = {a: v + seed_shift * s + (noise[a][s] if noise and a in noise else 0.0) for a, v in scores.items()}
        out.append({"seed": s, "scores": sc, "skills": {a: {"de_score": 0.5} for a in sc},
                    "t3": {a: {"severity": {"state": "r2_gate" if s < gated.get(a, 0) else "scored", "beta": 0.1,
                                            "r2": 0.05}} for a in sc}})
    return out


def _table(barred: float, refused: float = 55.0) -> dict[str, float]:
    return {harness.FLOOR: 50.0, k18.REFUSED: refused, k18.BARRED: barred, k18.RESID_K7: 54.0,
            k18.RESID_K7_PERM: 53.0, k18.RESID_UP: 55.5, k18.RESID_UP_PERM: 54.5}


def _per(barred: dict, gated: dict | None = None) -> dict:
    out = {}
    for w in k18.WORLDS:
        arms = (harness.FLOOR, *k18.ARMS) if w == "plain" else k18.ARMS
        out[w] = {n: k18.summarise(_rows({a: _table(barred[w, n])[a] for a in arms},
                                         (gated or {}).get((w, n))), arms) for n in k18.TASKS}
    return out


ALL = {(w, n): 56.0 for w in k18.WORLDS for n in k18.TASKS}
OK = {n: {"distance": None, "within_tolerance": True} for n in k18.TASKS}  # the depth-matched check, passed


def test_the_rule_needs_both_tasks_in_both_worlds_and_a_scored_severity():
    assert k18.decide(_per(ALL), OK)["question"]
    for k in ALL:  # one world and task short of the gain fails the whole
        assert not k18.decide(_per(ALL | {k: 55.2}), OK)["question"]
    two_gated = {("dm", "g2"): {k18.BARRED: 2}}
    res = k18.decide(_per(ALL, two_gated), OK)
    assert not res["question"] and res["dm_barred_severity_scored"] == {"g1": 10, "g2": 8}
    assert res["dm_refused_severity_scored"] == {"g1": 10, "g2": 10}
    assert k18.decide(_per(ALL, {("dm", "g2"): {k18.BARRED: 1}}), OK)["question"]  # 9 of 10 is enough
    assert k18.decide(_per(ALL, {("plain", "g2"): {k18.BARRED: 5}}), OK)["question"]  # the plain count is not read
    p = _per(ALL)
    p["dm"]["g1"]["paired"][f"{k18.BARRED} - {k18.REFUSED}"]["ahead"] = 7
    assert not k18.decide(p, OK)["question"]
    # the check can only block: outside tolerance on one replicate, missing, or not computable, no question
    for check in (OK | {"g2": {"distance": None, "within_tolerance": False}}, None, {}, {"g1": OK["g1"]}):
        res = k18.decide(_per(ALL), check)
        assert not res["question"] and not res["dm_check_ok"]
    # refused_form's scored count is reported and decides nothing (the rule's severity condition is barred's alone)
    res = k18.decide(_per(ALL, {("dm", "g1"): {k18.REFUSED: 6}}), OK)
    assert res["question"] and res["dm_refused_severity_scored"]["g1"] == 4


@pytest.mark.parametrize("shared", [False, True])
def test_the_rule_on_pure_noise_almost_never_asks(shared):
    """barred_form = refused_form + noise of sd 1 board point per seed, no true gain: the rule needs four (world, task)
    passes at once and asks in well under 1% of draws. With `shared`, the plain world's g1 and g2 carry the same noise
    (in that world they hold the same predicted cells: one pool, one stream, one seed), the worst case of that
    correlation."""
    rng = np.random.default_rng(11)
    asked = 0
    for _ in range(500):
        per = {}
        for w in k18.WORLDS:
            per[w] = {}
            common = rng.normal(0, 1.0, len(SEEDS))
            for n in k18.TASKS:
                noise = {k18.BARRED: common if shared and w == "plain" else rng.normal(0, 1.0, len(SEEDS))}
                per[w][n] = k18.summarise(_rows({k18.REFUSED: 55.0, k18.BARRED: 55.0}, noise=noise),
                                          (k18.REFUSED, k18.BARRED))
        asked += k18.decide(per, OK)["question"]
    assert asked <= 2


# --- the reading of a run folder ---


def _write_run(root: Path, barred=ALL, commit="abc", clean=True, machine="m1", guard_shift=0.0, dm_machine=None,
               dm_commit=None, check_shift=0.0, drop_floor=False):
    root.mkdir(parents=True)
    for n in k18.TASKS:
        arms = (harness.FLOOR, *k18.ARMS)
        rows = _rows({a: _table(barred["plain", n])[a] for a in arms})
        if drop_floor:  # the floor missing from some rows (not the first): must void, not raise
            for r in rows[3:]:
                del r["scores"][harness.FLOOR]
        for r, ref in zip(rows, sorted(k18.committed_k8(n), key=lambda x: x["seed"]), strict=True):
            r["scores"][k18.REFUSED] = ref["scores"][k8.GROUPED] + guard_shift
            r["scores"][k18.BARRED] = r["scores"][k18.REFUSED] + barred["plain", n] - 55.0
        (root / f"plain_{n}.json").write_text(json.dumps({
            "task": {"key": k18.TASKS[n].key}, "rows": rows,
            "provenance": {"commit": commit, "tree_clean": clean, "machine": machine}}))
    dm_rows = {n: [r | {"check": {k18.DM_CHECK: 60.0 + r["seed"] + check_shift}} for r in
                   _rows({a: _table(barred["dm", n])[a] for a in k18.ARMS})] for n in k18.TASKS}
    (root / "dm.json").write_text(json.dumps({
        "world": "dm", "rows": dm_rows,
        "provenance": {"commit": dm_commit or commit, "tree_clean": True, "machine": dm_machine or machine}}))
    return root


@pytest.fixture
def u3d(tmp_path):
    rows = [{"target": t, "seed": s, "world": w, "variant": v, "t3": 60.0 + s + (0.05 if t == "rep2" else 0.0)}
            for t in ("rep1", "rep2") for s in SEEDS for w, v in
            (("depth_matched", "dm:transfer:mab21l2_all"), ("depth_matched", "dm:xrep:F+depth"),
             ("as_released", "transfer:mab21l2_all"))]
    p = tmp_path / "orth_records.csv"
    pd.DataFrame(rows).to_csv(p, index=False)
    return p


def test_the_reading_of_a_run_folder(tmp_path, u3d):
    res = k18.reading(_write_run(tmp_path / "ask"), u3d=u3d)
    assert not res["void"] and res["question"] and res["machine"] == "m1"
    assert set(res["task_sha256"]) == {"plain_g1.json", "plain_g2.json", "dm.json"}
    assert all(k18.guard_passes(d) for d in res["guard"].values())
    assert res["dm_check"]["g1"]["within_tolerance"] and res["dm_check"]["g2"]["distance"]["mean_diff"] == \
        pytest.approx(-0.05)
    txt = k18.text(res)
    assert "WORLD plain" in txt and "WORLD dm" in txt and "a question to the captain" in txt
    assert "resid_perm_up - resid_revert_up" in txt and "refused_form - floor" in txt

    stands = k18.reading(_write_run(tmp_path / "stands", barred=ALL | {("dm", "g1"): 55.1}), u3d=u3d)
    assert not stands["void"] and not stands["question"] and "refused s1 stands" in stands["outcome"]

    off = k18.reading(_write_run(tmp_path / "check", check_shift=0.2), u3d=u3d)  # the check blocks, never voids
    assert not off["void"] and not off["question"] and not off["dm_check_ok"]
    assert not off["dm_check"]["g1"]["within_tolerance"] and "refused s1 stands" in off["outcome"]

    for name, kw in (("guard", {"guard_shift": 0.2}), ("dirty", {"clean": False}),
                     ("machines", {"dm_machine": "m2"}), ("commits", {"dm_commit": "def"}),
                     ("floor", {"drop_floor": True})):
        v = k18.reading(_write_run(tmp_path / name, **kw), u3d=u3d)
        assert v["void"] and "no question" in v["outcome"], name
        assert "VOID" in k18.text(v)
    assert any("different commits" in w for w in k18.reading(tmp_path / "commits", u3d=u3d)["why"])
    assert any("missing arms ['floor']" in w for w in k18.reading(tmp_path / "floor", u3d=u3d)["why"])


def test_a_missing_task_file_is_not_a_reading(tmp_path, u3d, capsys):
    """A run that crashed or was killed leaves no task file: the reading opens nothing, is not void, writes nothing
    and exits 3, so the world can be rerun blind."""
    root = _write_run(tmp_path / "run")
    (root / "dm.json").write_text("{ not JSON: the reading must not open this file")
    (root / "plain_g2.json").unlink()
    res = k18.reading(root, u3d=u3d)
    assert res == {"incomplete": True, "void": False, "missing": ["plain_g2.json"], "outcome": k18.INCOMPLETE}
    assert "missing task files: plain_g2.json" in k18.text(res)
    out = tmp_path / "reading"
    with pytest.raises(SystemExit) as e:
        k18.main(["read", "--runs", str(root), "--out", str(out)])
    assert e.value.code == 3 and not out.exists()
    assert "not a reading" in capsys.readouterr().out


def test_a_task_file_and_its_provenance(tmp_path):
    res = {"rows": [], "provenance": {"commit": "abc", "tree_clean": True, "machine": "m1", "utc": "x"}}
    prov = k18.write_task_file(tmp_path / "dm.json", res)
    assert prov["sha256"] == __import__("hashlib").sha256((tmp_path / "dm.json").read_bytes()).hexdigest()
    assert json.loads((tmp_path / "dm.json.prov.json").read_text()) == prov
    assert {k: prov[k] for k in ("commit", "tree_clean", "machine")} == {"commit": "abc", "tree_clean": True,
                                                                         "machine": "m1"}
    assert sorted(p.name for p in tmp_path.iterdir()) == ["dm.json", "dm.json.prov.json"]  # no temporary file left


# --- the reruns and the logs (row K18: no log holds a score; a code fix means a fresh folder) ---


def test_a_single_world_rerun_after_a_code_fix_cannot_be_read(tmp_path, u3d):
    """A folder whose dm.json was rerun at a later commit (a code fix) beside the plain files of the first commit is
    void: the one-commit check reads it as two runs. So after a code fix the three worlds go into a fresh folder."""
    root = _write_run(tmp_path / "run", commit="abc")
    dm = json.loads((root / "dm.json").read_text())
    dm["provenance"]["commit"] = "later"
    (root / "dm.json").write_text(json.dumps(dm))
    v = k18.reading(root, u3d=u3d)
    assert v["void"] and not v.get("question") and any("different commits" in w for w in v["why"])
    assert "rerun all three worlds at the new commit into a fresh folder" in k18.INCOMPLETE


def test_a_dm_row_without_the_check_voids_and_does_not_raise(tmp_path, u3d):
    root = _write_run(tmp_path / "run")
    dm = json.loads((root / "dm.json").read_text())
    del dm["rows"]["g2"][4]["check"]
    (root / "dm.json").write_text(json.dumps(dm))
    v = k18.reading(root, u3d=u3d)
    assert v["void"] and any("without the check arm" in w for w in v["why"])
    assert not v["dm_check"]["g2"]["within_tolerance"] and v["dm_check"]["g1"]["within_tolerance"]


@pytest.fixture
def fake_run(monkeypatch):
    """cmd_run with the barred factor's load, the harness and the stamp stubbed: the fake harness logs each seed's
    line in the harness's own format (every arm's score) and returns a minimal task result."""
    calls = {"log": None, "n": 0}
    lines = {}

    def run_task(t, cands, seeds, data_dir, log=print):
        calls["log"], calls["n"] = log, calls["n"] + 1
        for s in seeds:
            lines[s] = f"{t.key} seed {s}: floor 50.00, refused_form {55 + s / 7:.2f}, barred_form {60 - s / 9:.2f}"
            log(lines[s])
        return {"rows": [], "provenance": {"commit": "abc", "tree_clean": True}}

    monkeypatch.setattr(k18, "load_barred_for_plain", lambda *a: {})
    monkeypatch.setattr(harness, "run_task", run_task)
    monkeypatch.setattr(k7, "_PREP", {"fold": {"record": {}}})
    monkeypatch.setattr(k18, "stamp", lambda: {"commit": "abc", "tree_clean": True, "machine": "m1", "utc": "x"})
    return SimpleNamespace(calls=calls, lines=lines)


def test_the_plain_run_logs_progress_only(tmp_path, fake_run, capsys):
    out = tmp_path / "run"
    k18.main(["run", "--world", "plain", "--task", "g1", "--out", str(out)])
    assert fake_run.calls["log"] is not print and fake_run.calls["n"] == 1
    printed = capsys.readouterr()
    log_lines = [x for x in printed.out.splitlines() if " seed " in x]
    assert log_lines == [f"K18 plain g1 seed {s}: scored" for s in SEEDS]
    for x in log_lines:  # no number after the seed
        assert not re.search(r"\d", x.split(" seed ", 1)[1].split(":", 1)[1])
    for line in fake_run.lines.values():  # no score of the harness's line reaches the output
        for score in re.findall(r"\d+\.\d\d", line):
            assert score not in printed.out and score not in printed.err
    assert (out / "plain_g1.json").exists() and (out / "plain_g1.json.prov.json").exists()


def test_a_crash_prints_a_masked_traceback_and_writes_no_task_file(tmp_path, fake_run, monkeypatch, capsys):
    def boom(*a, **k):
        raise ValueError("T3w/WT->Gata4rep1 seed 0: shared anchors give 57.1234, fresh ones 57.1301; -1.5e-03 nan")

    monkeypatch.setattr(harness, "run_task", boom)
    out = tmp_path / "run"
    with pytest.raises(SystemExit) as e:
        k18.main(["run", "--world", "plain", "--task", "g1", "--out", str(out)])
    assert e.value.code == 1
    err = capsys.readouterr().err
    assert "ValueError" in err and "seed 0" in err and "Traceback" in err
    for x in ("57.1234", "57.1301", "1.5e-03", "nan"):
        assert x not in err
    assert not (out / "plain_g1.json").exists() and not (out / "plain_g1.json.prov.json").exists()


def test_masking_on_pure_noise():
    """Random scores in every format the harness or numpy prints: none survives the mask; integers stay."""
    rng = np.random.default_rng(4)
    for _ in range(300):
        v = rng.normal(50, 20, 6)
        txt = (f"seed 17: a {v[0]:.2f}, b {v[1]!r}, c {v[2]:.4e}, d {-abs(v[3]):.6f}, e {v[4]:g}, "
               f"f np.float64({v[5]})")
        m = k18.masked(txt)
        assert not re.search(r"\d\.\d|\d[eE][-+]?\d", m) and m.startswith("seed 17: ")


def test_a_run_refuses_a_folder_of_another_commit_or_machine(tmp_path, fake_run):
    out = tmp_path / "run"
    out.mkdir()
    for commit, machine, refused in (("old", "m1", True), ("abc", "m2", True), ("abc", "m1", False)):
        (out / "dm.json.prov.json").write_text(json.dumps({"commit": commit, "machine": machine, "file": "dm.json"}))
        if refused:
            with pytest.raises(SystemExit) as e:
                k18.main(["run", "--world", "plain", "--task", "g2", "--out", str(out)])
            assert "fresh folder" in str(e.value.code) and fake_run.calls["n"] == 0
            assert not (out / "plain_g2.json").exists()
        else:  # the same commit and machine: a rerun of a missing world, at the row's commit
            k18.main(["run", "--world", "plain", "--task", "g2", "--out", str(out)])
            assert fake_run.calls["n"] == 1 and (out / "plain_g2.json").exists()


def test_a_kill_between_the_two_writes_leaves_no_task_file(tmp_path, monkeypatch):
    real, n = k18._write_by_rename, []

    def once(path, raw):
        n.append(path.name)
        if len(n) == 2:
            raise KeyboardInterrupt  # the kill, after the provenance and before the task file
        real(path, raw)

    monkeypatch.setattr(k18, "_write_by_rename", once)
    res = {"rows": [], "provenance": {"commit": "abc", "tree_clean": True, "machine": "m1"}}
    with pytest.raises(KeyboardInterrupt):
        k18.write_task_file(tmp_path / "dm.json", res)
    assert n == ["dm.json.prov.json", "dm.json"] and not (tmp_path / "dm.json").exists()


# --- the third fix round: the dm self-check's score kept out of the log; a dirty tree is INCOMPLETE, never VOID ---


def test_the_dm_self_check_prints_no_score(capsys):
    """The frozen selfcheck's own line ends with a board score; quiet_selfcheck drops it and prints a fixed line.
    Pure noise: random scores printed through the frozen script's own log() never reach the output. A failing
    self-check still raises (cmd_run then masks the traceback). run_dm calls quiet_selfcheck, never selfcheck."""
    import inspect

    t3d = harness._t3_head().t3d
    rng = np.random.default_rng(8)
    for _ in range(200):
        v = float(rng.normal(55, 10))
        fake = SimpleNamespace(selfcheck=lambda ctx, seed, v=v: t3d.log(
            f"self-check passed: factored scoring equals vec.evaluate.score_board ({v:.4f})"))
        k18.quiet_selfcheck(fake, None)
        out = capsys.readouterr()
        assert out.out == k18.SELFCHECK_LINE + "\n" and out.err == ""
        assert not re.search(r"\d", out.out) and f"{v:.4f}" not in out.out

    def failing(ctx, seed):
        raise RuntimeError("self-check failed: factored scoring 57.1 vs vec.evaluate.score_board 57.2")

    with pytest.raises(RuntimeError):
        k18.quiet_selfcheck(SimpleNamespace(selfcheck=failing), None)
    assert k18.SELFCHECK_LINE not in capsys.readouterr().out
    src = inspect.getsource(k18.run_dm)
    assert "quiet_selfcheck(t3d, ctx)" in src and ".selfcheck(" not in src


@pytest.mark.parametrize("world", ["plain", "dm"])
def test_a_dirty_tree_is_refused_before_anything_loads(tmp_path, fake_run, monkeypatch, world):
    """A stamp whose tree is not clean: the run exits before loading the barred factor, the harness or the dm world,
    makes no folder and writes no task file, so the case is INCOMPLETE (read exits 3), never VOID."""
    loaded = []
    monkeypatch.setattr(k18, "stamp", lambda: {"commit": "abc", "tree_clean": False, "machine": "m1", "utc": "x"})
    monkeypatch.setattr(k18, "load_barred_for_plain", lambda *a: loaded.append("barred"))
    monkeypatch.setattr(k18, "run_dm", lambda *a: loaded.append("dm"))
    out = tmp_path / "run"
    argv = ["run", "--world", world, "--out", str(out)] + (["--task", "g1"] if world == "plain" else [])
    with pytest.raises(SystemExit) as e:
        k18.main(argv)
    assert e.value.code == k18.DIRTY and "INCOMPLETE, not VOID" in k18.DIRTY
    assert loaded == [] and fake_run.calls["n"] == 0 and not out.exists()
    res = k18.reading(out)
    assert res["incomplete"] and not res["void"] and len(res["missing"]) == 3


def test_the_stamp_rule_is_the_harness_rule(monkeypatch):
    """stamp()'s clean flag is vec.harness's: `git status --porcelain --untracked-files=no` prints nothing. A tracked
    change (any porcelain line) or a failing git ('unknown') is unclean."""
    seen = []
    status = {"now": ""}

    def git(*args):
        seen.append(args)
        return "abc" if args[0] == "rev-parse" else status["now"]

    monkeypatch.setattr(harness, "_git", git)
    for line, clean in (("", True), (" M progress.md", False), ("M  docs/x.md", False), (" D a.py", False),
                        ("UU b.py", False), ("unknown", False)):
        status["now"] = line
        assert k18.stamp()["tree_clean"] is clean, line
    assert ("status", "--porcelain", "--untracked-files=no") in seen


@pytest.mark.parametrize("prov", [{"commit": "abc", "tree_clean": False}, {"commit": "later", "tree_clean": True}])
def test_a_stamp_that_turned_dirty_or_moved_writes_no_task_file(tmp_path, fake_run, monkeypatch, prov):
    """The plain world's own stamp is the harness's, taken a little after the run's start: if it says unclean or
    another commit, no task file is written (INCOMPLETE), so the reading never meets a file it would void."""

    def run_task(t, cands, seeds, data_dir, log=print):
        return {"rows": [], "provenance": prov}

    monkeypatch.setattr(harness, "run_task", run_task)
    out = tmp_path / "run"
    with pytest.raises(SystemExit) as e:
        k18.main(["run", "--world", "plain", "--task", "g1", "--out", str(out)])
    assert "INCOMPLETE" in str(e.value.code)
    assert not (out / "plain_g1.json").exists() and not (out / "plain_g1.json.prov.json").exists()


def test_a_crash_in_read_is_masked_and_writes_nothing(tmp_path, monkeypatch, capsys):
    def boom(runs):
        raise KeyError("barred_form mean 57.1234 against 55.4321")

    monkeypatch.setattr(k18, "reading", boom)
    out = tmp_path / "reading"
    with pytest.raises(SystemExit) as e:
        k18.main(["read", "--runs", str(tmp_path), "--out", str(out)])
    err = capsys.readouterr().err
    assert e.value.code == 1 and "KeyError" in err and "57.1234" not in err and "55.4321" not in err
    assert not out.exists()
