"""Try K16 (experiments/harness/TRIES.md; experiments/t3_decomp/k16_revert_dose.py): the exponent of the wild-type-only
factor. Synthetic data and hand-made score tables only (no competition file, no score of a fold)."""

from __future__ import annotations

import json
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
import k16_revert_dose as k16  # noqa: E402

G, N_POOL, N_DRAW = 14, 90, 150
GENES = [f"g{i}" for i in range(G)]


def test_arms_tasks_and_constants():
    assert k16.BASE == "grp_e050" and k16.BASE_EXPONENT == k7.EXPONENT == 0.5
    assert k16.ARMS == ("grp_e050", "grp_e075", "grp_e100", "grp_e150", "perm_e050", "perm_e100")
    assert k16.TRIED == (0.75, 1.0, 1.5)
    assert k16.TASKS == k7.TASKS and set(k16.TASKS) == {"g1", "g2", "m"}
    assert k16.S2_KEYS == {k7.TASKS["g1"].key, k7.TASKS["g2"].key}
    assert (k16.MIN_GAIN, k16.MIN_AHEAD, k16.MIN_SCORED, k16.M_FLOOR) == (0.3, 8, 9, -0.3)
    cands = k16.candidates()
    assert [c.name for c in cands] == [harness.FLOOR, *k16.ARMS]
    assert all(c.applies(t) for c in cands[1:] for t in k16.TASKS.values())


# --- the builds, on a synthetic reference with every builder call the scripts share stubbed the same way ---


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
    basis = harness.t3_mean_expression(pool)
    prep = {"ref": ref, "phi": phi, "per": per, "group": group, "record": {},
            "phi_perm": harness.t3_decile_permutation(phi, basis, k7.PERM_SEED)}
    f_s2 = np.array([0.6, 0.25, 0.15])
    log = {"construct": [], "composition": [], "assemble": []}

    def construct(r, n, g, f=None, phi=None, w=1.0, ratio=None, alphas=None):
        assert r is ref and phi is None and ratio is None and alphas is None
        log["construct"].append(np.asarray(f).copy())
        rows = g.integers(0, N_POOL, n)
        return bs.Built(X=bs.t3d.lognorm(C[rows]), lab=lab[rows], source=np.full(n, "wt", dtype=object), row=rows)

    def composition(v, spec, params):
        log["composition"].append(spec)
        return f_s2, {"states_missing": []}

    def assemble(r, built, panel, g):
        log["assemble"].append({"X": built.X.copy(), "row": built.row.copy(), "lab": built.lab.copy(),
                                "draw": g.random()})
        return built.X, None, {}

    monkeypatch.setattr(k7, "_prepared", lambda ctx: prep)
    monkeypatch.setattr(bs, "construct", construct)
    monkeypatch.setattr(bs, "composition", composition)
    monkeypatch.setattr(bs, "assemble", assemble)
    monkeypatch.setattr(harness, "count_for", lambda board, n: N_DRAW)
    monkeypatch.setattr(harness, "finish", lambda ctx, p, X, coords: X)
    k16._PERM.clear()

    def ctx(task: str, seed: int):
        return SimpleNamespace(task=SimpleNamespace(key=k16.TASKS[task].key, reference_role="wt"),
                               inputs={"wt": pool}, board=SimpleNamespace(genes=GENES, needs_coords=False), seed=seed)

    return SimpleNamespace(bs=bs, ctx=ctx, prep=prep, basis=basis, log=log, f_s2=f_s2, V=V)


@pytest.mark.parametrize("seed", [0, 7])
def test_the_base_arm_is_k8b_on_the_replicates_and_k7b_on_m(fold, seed):
    for task, other, other_arm in (("g1", k8, k8.GROUPED), ("g2", k8, k8.GROUPED), ("m", k7, k7.GROUPED)):
        fold.log["assemble"].clear()
        fold.log["construct"].clear()
        want = other._build(fold.ctx(task, seed), other_arm)
        got = k16._build(fold.ctx(task, seed), k16.BASE)
        assert np.array_equal(got, want)
        a, b = fold.log["assemble"]
        assert a["draw"] == b["draw"] and np.array_equal(a["row"], b["row"])  # the same assemble stream
        assert len(fold.log["construct"]) == 2 and np.array_equal(*fold.log["construct"])
        assert np.array_equal(fold.log["construct"][0], fold.V.p if task == "m" else fold.f_s2)
        fold.log["construct"].clear()
    assert len(fold.log["composition"]) == 4 and all(s == k8.spec() for s in fold.log["composition"])
    perm = k16._build(fold.ctx("g1", seed), k16.arm(0.5, True))
    assert not np.array_equal(perm, k16._build(fold.ctx("g1", seed), k16.BASE))


def test_the_exponent_reaches_the_factor_as_a_power_and_nothing_else_changes(fold, monkeypatch):
    calls: dict[str, list] = {}
    real = fold.bs.t3d.apply_gene_factor

    def spy(C, phi, P, rng):
        calls[current].append({"C": C.copy(), "phi": np.asarray(phi).copy(), "P": P.copy(),
                               "state": json.dumps(rng.bit_generator.state, default=int)})
        return real(C, phi, P, rng)

    monkeypatch.setattr(fold.bs.t3d, "apply_gene_factor", spy)
    for task in ("g2", "m"):
        cells = {}
        for name in k16.ARMS:
            current = name
            calls[name] = []
            fold.log["assemble"].clear()
            k16._build(fold.ctx(task, 3), name)
            cells[name] = fold.log["assemble"][0]
        base = calls[k16.BASE]
        for name in k16.ARMS:
            e, perm = k16.ARM_SPEC[name]
            # every arm holds the same cells: rows, labels and the assemble stream
            assert np.array_equal(cells[name]["row"], cells[k16.BASE]["row"])
            assert np.array_equal(cells[name]["lab"], cells[k16.BASE]["lab"])
            assert cells[name]["draw"] == cells[k16.BASE]["draw"]
            assert len(calls[name]) == len(base) == 3  # A, B and the fallback group, in sorted order
            assert calls[name][0]["state"] == base[0]["state"]  # the per-group step starts from the same generator
            for c, b, g in zip(calls[name], base, ("A", "B", "Unknown"), strict=True):
                assert np.array_equal(c["C"], b["C"]) and np.array_equal(c["P"], b["P"])
                src = fold.prep["per"].get(g, fold.prep["phi"])
                if perm:  # K7c's permutation of the same group's factor, the same gene permutation for every group
                    src = harness.t3_decile_permutation(src, fold.basis, k7.PERM_SEED)
                assert np.array_equal(c["phi"], src ** e)
        # the fallback group's permuted factor is K7c's own control vector
        assert np.array_equal(calls[k16.arm(1.0, True)][2]["phi"], fold.prep["phi_perm"] ** 1.0)


# --- the reading, on hand-made score tables ---


def _result(key: str, scores: dict[str, float], gated: frozenset = frozenset(), n_gated: int = 0) -> dict:
    """Ten seeds; `scores` per arm (a constant offset per seed keeps every difference the same in every seed); the arms
    in `gated` are gated in their first `n_gated` seeds."""
    arms = {harness.FLOOR: 50.0, **scores}
    rows = []
    for s in SEEDS:
        rows.append({"seed": s, "scores": {a: v + 0.1 * s for a, v in arms.items()},
                     "skills": {a: {"de_score": 0.5, "severity": 0.6} for a in arms},
                     "t3": {a: {"severity": {"state": "r2_gate" if a in gated and s < n_gated else "scored",
                                             "beta": 0.1 * (1 + s), "r2": 0.05}} for a in arms}})
    return {"task": {"key": key}, "rows": rows, "provenance": {"commit": "abc", "tree_clean": True}}


def _table(e075, e100, e150, p050=50.0, p100=50.0, base=50.0) -> dict[str, float]:
    return {"grp_e050": base, "grp_e075": e075, "grp_e100": e100, "grp_e150": e150, "perm_e050": p050,
            "perm_e100": p100}


def _read(tmp_path, name, g, m=None, **kw):
    d = tmp_path / name
    d.mkdir()
    for n in ("g1", "g2"):
        (d / f"{n}.json").write_text(json.dumps(_result(k16.TASKS[n].key, g, **kw)))
    (d / "m.json").write_text(json.dumps(_result(k16.TASKS["m"].key, m or _table(50, 50, 50))))
    return k16.reading(d)


def test_every_exponent_carries_and_rises(tmp_path):
    r = _read(tmp_path, "rise", _table(51.0, 52.0, 53.0))
    assert not r["void"] and r["carries"] == {"grp_e075": True, "grp_e100": True, "grp_e150": True}
    assert r["factor_specific"] and r["named_exponent"] == 1.5 and r["request"]
    assert r["tasks"]["g1"]["paired"]["grp_e100 - grp_e050"]["mean"] == pytest.approx(2.0)
    assert r["tasks"]["g1"]["arms"]["grp_e100"]["severity_scored"] == 10
    assert r["tasks"]["g1"]["arms"]["grp_e100"]["median_beta"] == pytest.approx(0.55)
    assert "priced question" in k16.text(r)


def test_not_factor_specific_means_no_request(tmp_path):
    r = _read(tmp_path, "perm", _table(51.0, 52.0, 53.0, p050=50.0, p100=51.8))  # control gains 1.8 of the 2.0
    assert all(r["carries"].values()) and r["named_exponent"] == 1.5
    assert r["specific_gain_over_control"]["g1"] == pytest.approx(0.2)
    assert not r["factor_specific"] and not r["request"] and "stays at exponent 0.5" in r["outcome"]


def test_the_monotone_step(tmp_path):
    flat_top = _read(tmp_path, "flat_top", _table(51.0, 52.0, 52.2))  # 1.5 carries but is not 0.3 above 1.0
    assert flat_top["carries"]["grp_e150"] and not flat_top["step"]["grp_e150"]
    assert flat_top["named_exponent"] == 1.0 and flat_top["request"]
    plateau = _read(tmp_path, "plateau", _table(51.0, 51.2, 53.0))  # 1.0 is not above 0.75: no step past 0.75
    assert not plateau["step"]["grp_e100"] and plateau["step"]["grp_e150"]
    assert plateau["named_exponent"] == 0.75
    # 0.75 does not carry, so 1.0 needs no step over it; at the factor-specific check the control is level
    waived = _read(tmp_path, "waived", _table(50.1, 51.0, 51.1))
    assert not waived["carries"]["grp_e075"] and waived["step"]["grp_e100"] and waived["named_exponent"] == 1.0
    none = _read(tmp_path, "none", _table(49.0, 48.0, 47.0))
    assert not any(none["carries"].values()) and none["named_exponent"] is None and not none["request"]


def test_carrying_needs_severity_scored_and_m_level(tmp_path):
    gated = _read(tmp_path, "gated", _table(51.0, 52.0, 53.0), gated=frozenset({"grp_e150"}), n_gated=2)
    assert gated["tasks"]["g2"]["arms"]["grp_e150"]["severity_scored"] == 8
    assert not gated["carries"]["grp_e150"] and gated["named_exponent"] == 1.0
    m_low = _read(tmp_path, "m_low", _table(51.0, 52.0, 53.0), m=_table(50.0, 49.6, 49.8))
    assert m_low["carries"] == {"grp_e075": True, "grp_e100": False, "grp_e150": True}
    assert m_low["named_exponent"] == 1.5  # 1.0 did not carry, yet 1.5 carries and beats it by more than 0.3


def test_a_missing_or_dirty_task_voids_the_reading(tmp_path):
    d = tmp_path / "missing"
    d.mkdir()
    (d / "g1.json").write_text(json.dumps(_result(k16.TASKS["g1"].key, _table(51, 52, 53))))
    assert k16.reading(d)["void"]
    d = tmp_path / "dirty"
    d.mkdir()
    for n in k16.TASKS:
        r = _result(k16.TASKS[n].key, _table(51, 52, 53))
        r["provenance"]["tree_clean"] = n != "m"
        (d / f"{n}.json").write_text(json.dumps(r))
    assert k16.reading(d)["void"]
