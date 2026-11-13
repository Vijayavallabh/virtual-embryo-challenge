"""Tests for tries K12 (experiments/t3_decomp/k12_r8_arm.py), on synthetic data only: no competition file is opened,
nothing is scored, and no prediction is compared with any target."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pytest

from vec import harness
from vec.adoption import SEEDS

T3 = Path(__file__).resolve().parents[1] / "experiments" / "t3_decomp"
if str(T3) not in sys.path:
    sys.path.insert(0, str(T3))
import eb_factor as eb  # noqa: E402
import k7_wt_revert as k7  # noqa: E402
import k9_residual as k9  # noqa: E402
import k12_r8_arm as k12  # noqa: E402
import test_k9_residual as t9  # noqa: E402  (K9's synthetic world and reference, reused)

CLIP = t9.CLIP
world = t9.world  # the module fixture: synthetic files in task x1's roles


@pytest.fixture(scope="module")
def bs():
    return harness._t3_head()


# ----------------------------------------------------------------------------------------------- constants and tasks
def test_the_fixed_constants_and_arms():
    assert k12.PERM_SEED == k9.PERM_SEED == 20261113 and k12.STREAM == k9.STREAM == harness.K4_STREAM
    assert k12.GATE_PERMS == 200 and k12.GATE_REFERENCE == "E8.75"
    assert (k12.MIN_GAIN, k12.MIN_AHEAD, k12.X2_FLOOR) == (0.3, 8, -0.3)
    assert k12.ARMS == ("raw", "raw_r8", "raw_r8_perm") and k12.CONTROLS == {"raw_r8": "raw_r8_perm"}
    assert k12.TASKS is k9.TASKS
    assert [c.name for c in k12.candidates()] == [harness.FLOOR, *k12.ARMS]
    assert all(c.applies(t) for c in k12.candidates()[1:] for t in k12.TASKS.values())
    assert not any(c.applies(t) for c in k12.candidates()[1:] for ts in harness.load_tasks().values() for t in ts)


def test_the_gate_loads_only_the_reference_and_the_wild_type_series():
    for stem, name in k12.GATE_TASKS.items():
        t = k12.TASKS[name]
        files = k12.gate_files(name)
        assert t.reference == f"T2/{stem}.h5ad" and files[harness.base_role(t)] == t.reference
        assert set(files.values()) == {f"T2/{s}.h5ad" for s in k7.FULL_PANEL}
        assert t.target not in files.values() and not any(f.startswith("T3/") for f in files.values())
    assert k12.TASKS["r12"].reference == k12.TASKS["x1"].reference  # r12 shares x1's reference, so its gate


# ----------------------------------------------------------------------------------------------- the contrasts
def test_the_interpolation_weights_are_r2s():
    d = harness.STAGE_DAYS
    assert k12.interior_weight(d["T2/E8.5.h5ad"], d["T2/E8.75.h5ad"], d["T2/E9.5.h5ad"]) == pytest.approx(0.75)
    assert k12.interior_weight(d["T2/E8.75.h5ad"], d["T2/E9.5.h5ad"], d["T2/E10.5.h5ad"]) == pytest.approx(4 / 7)
    assert round(4 / 7, 3) == 0.571 and round(3 / 7, 3) == 0.429
    rng = np.random.default_rng(0)
    m = [{"A": rng.dirichlet(np.ones(8)), "B": rng.dirichlet(np.ones(8))} for _ in range(4)]
    days = [8.5, 8.75, 9.5, 10.5]
    M, labels = k12.run_contrasts(m, days, ["A", "B"])
    assert labels == [(1, "A"), (1, "B"), (2, "A"), (2, "B")]
    ly = [{g: np.log(x[g] + k7.EPS) for g in x} for x in m]
    assert np.allclose(M[0], ly[1]["A"] - (0.75 * ly[0]["A"] + 0.25 * ly[2]["A"]))
    assert np.allclose(M[3], ly[2]["B"] - (4 / 7 * ly[1]["B"] + 3 / 7 * ly[3]["B"]))
    # a group missing at any of the three stages is left out of that stage's contrast
    del m[3]["B"]
    M2, labels2 = k12.run_contrasts(m, days, ["A", "B"])
    assert labels2 == [(1, "A"), (1, "B"), (2, "A")] and np.array_equal(M2, M[:3])
    # a straight line in days leaves no residual
    line = [{"A": np.exp(0.3 * x + np.linspace(-3, -1, 8)) - k7.EPS} for x in days]
    assert np.allclose(k12.run_contrasts(line, days, ["A"])[0], 0.0)


def test_the_cleaning_step_removes_a_planted_monotone_part():
    rng = np.random.default_rng(1)
    n = 300
    level = rng.gamma(2.0, size=n)
    mono = np.tanh(level - 2.0) * 1.5
    v, direction = k12.level_monotone_fit(mono, level)
    assert direction == "up" and np.allclose(v, mono)
    M, rec = k12.clean_contrasts(np.vstack([mono, -mono]), level)
    assert np.allclose(M, 0.0) and [r["direction"] for r in rec] == ["up", "down"]
    assert rec[0]["removed_energy_share"] == pytest.approx(1.0)
    # a planted monotone part plus a signal that does not follow level: what is left is close to the signal
    sig = np.zeros(n)
    sig[rng.choice(n, 30, replace=False)] = rng.choice([-1.0, 1.0], 30)
    (left,), _ = k12.clean_contrasts((mono + sig)[None, :], level)
    assert np.corrcoef(left, sig)[0, 1] > 0.9
    assert abs(np.corrcoef(left, mono)[0, 1]) < 0.3
    # genes of equal level share one fitted value
    lv = np.repeat(np.arange(5.0), 4)
    fit, _ = k12.level_monotone_fit(rng.normal(size=20), lv)
    assert all(np.unique(fit[lv == x]).size == 1 for x in range(5))


def test_a_planted_shared_direction_is_recovered_and_the_gate_says_inside():
    rng = np.random.default_rng(2)
    n_genes, n_c = 200, 12
    level = rng.gamma(2.0, size=n_genes)
    d = rng.normal(size=n_genes)
    d /= np.linalg.norm(d)
    M = rng.normal(1.0, 0.3, n_c)[:, None] * d[None, :] * np.sign(rng.normal(size=n_c))[:, None] \
        + 0.02 * rng.normal(size=(n_c, n_genes))
    u = k12.rank1_basis(M)
    assert abs(float(u @ d)) > 0.99 and np.linalg.norm(u) == pytest.approx(1.0)
    g = k12.gate(M, level)
    assert g["verdict"] == "inside" and g["G"] > g["null_p95"] and g["n_contrasts"] == n_c
    assert g["n_permutations"] == 200 and g["null_median"] <= g["null_p95"]


def test_independent_noise_contrasts_give_not_inside():
    rng = np.random.default_rng(3)
    level = rng.gamma(2.0, size=200)
    g = k12.gate(rng.normal(size=(12, 200)), level)
    assert g["verdict"] == "not inside" and g["G"] <= g["null_p95"]
    assert k12.gate(rng.normal(size=(1, 200)), level)["verdict"] == "not inside"  # too few contrasts
    assert k12.gate(np.zeros((0, 200)), level)["G"] is None


def test_the_gate_statistic_is_the_mean_left_out_cos2():
    rng = np.random.default_rng(4)
    M = rng.normal(size=(5, 30))
    want = []
    for j in range(5):
        _, _, Vt = np.linalg.svd(np.delete(M, j, axis=0), full_matrices=False)
        want.append((M[j] @ Vt[0]) ** 2 / (M[j] @ M[j]))
    assert k12.gate_statistic(M) == pytest.approx(float(np.mean(want)))
    # the null's draws use new seeds derived from PERM_SEED, never PERM_SEED itself
    a = harness.t3_decile_permutation(np.arange(30.0), rng.gamma(2.0, size=30), [k12.PERM_SEED, 0])
    b = harness.t3_decile_permutation(np.arange(30.0), rng.gamma(2.0, size=30), k12.PERM_SEED)
    assert sorted(a) == sorted(b)


# ----------------------------------------------------------------------------------------------- the factors
def test_raw_r8_removes_exactly_the_component_along_u():
    rng = np.random.default_rng(5)
    y = rng.normal(0, 0.7, 50)
    u = rng.normal(size=50)
    u /= np.linalg.norm(u)
    y2 = k12.remove_basis(y, u)
    assert abs(float(y2 @ u)) < 1e-12
    assert np.allclose(y - y2, (y @ u) * u)
    den = rng.gamma(2.0, size=50)
    phi_raw = eb.finalize(y, den, CLIP)
    assert np.array_equal(k12.r8_factor(y, den, CLIP, u, phi_raw), eb.finalize(y2, den, CLIP))
    assert k12.removed_share(y, u) == pytest.approx((y @ u) ** 2 / (y @ y))


def test_a_zero_basis_gives_raw_bit_for_bit(bs):
    rng = np.random.default_rng(6)
    y = rng.normal(0, 0.7, t9.G)
    den = rng.gamma(2.0, size=t9.G)
    phi_raw = np.clip(np.exp(rng.normal(0, 0.6, t9.G)), *CLIP)
    for u in (None, np.zeros(t9.G)):
        assert np.array_equal(k12.r8_factor(y, den, CLIP, u, phi_raw), phi_raw)
    assert k12.rank1_basis(np.zeros((3, t9.G))) is None and k12.rank1_basis(np.zeros((0, 0))) is None
    ref = t9._reference(bs, rng)
    prep = {"phi_raw": phi_raw, "phi_r8": k12.r8_factor(y, den, CLIP, None, phi_raw),
            "phi_r8_perm": k12.r8_factor(y, den, CLIP, k12.permuted_basis(None, den), phi_raw), "clip": CLIP}
    raw = k12.arm_built(bs, ref, 200, 7, ref.V.p, "raw", prep)
    k9raw = k9.arm_built(bs, ref, 200, 7, ref.V.p, "raw", phi_raw, None, None, CLIP)
    assert np.array_equal(raw.X, k9raw.X) and np.array_equal(raw.row, k9raw.row)
    for arm in ("raw_r8", "raw_r8_perm"):
        b = k12.arm_built(bs, ref, 200, 7, ref.V.p, arm, prep)
        assert np.array_equal(b.X, raw.X) and np.array_equal(b.row, raw.row) and np.array_equal(b.lab, raw.lab)


def test_the_control_uses_a_permuted_u_of_the_same_norm():
    rng = np.random.default_rng(7)
    n = 60
    basis = rng.gamma(2.0, size=n)
    u = rng.normal(size=n)
    u /= np.linalg.norm(u)
    up = k12.permuted_basis(u, basis)
    assert np.linalg.norm(up) == pytest.approx(1.0) and not np.array_equal(up, u)
    assert np.array_equal(up, harness.t3_decile_permutation(u, basis, k9.PERM_SEED))
    for genes in np.array_split(np.argsort(basis, kind="stable"), 10):
        assert np.array_equal(np.sort(up[genes]), np.sort(u[genes]))


def test_the_arms_draw_the_same_cells(bs):
    rng = np.random.default_rng(8)
    ref = t9._reference(bs, rng)
    y = rng.normal(0, 0.7, t9.G)
    den = rng.gamma(2.0, size=t9.G)
    phi_raw = eb.finalize(y, den, CLIP)
    u = rng.normal(size=t9.G)
    u /= np.linalg.norm(u)
    prep = {"phi_raw": phi_raw, "phi_r8": k12.r8_factor(y, den, CLIP, u, phi_raw),
            "phi_r8_perm": k12.r8_factor(y, den, CLIP, k12.permuted_basis(u, den), phi_raw), "clip": CLIP}
    built = {a: k12.arm_built(bs, ref, 200, 9, ref.V.p, a, prep) for a in k12.ARMS}
    assert all(np.array_equal(b.row, built["raw"].row) for b in built.values())
    assert not np.array_equal(built["raw_r8"].X, built["raw"].X)
    assert not np.array_equal(built["raw_r8"].X, built["raw_r8_perm"].X)


# ----------------------------------------------------------------------------------------------- end to end, synthetic
def test_end_to_end_raw_is_k9s_raw_and_raw_r8_moves_it(world, monkeypatch):
    monkeypatch.setattr(k7, "MIN_CELLS", 20)
    t = k12.TASKS["x1"]
    k12._PREP.clear()
    k9._PREP.clear()
    ctx = t9._ctx(t, world)
    with harness.guarded_reads() as g:
        out = {arm: k12._build(ctx, arm) for arm in k12.ARMS}
        k9raw = k9._build(ctx, k9.BASE)
    assert not g["refused"]
    prep = k12._PREP["fold"]
    assert np.array_equal(prep["phi_raw"], k9._PREP["fold"]["phi_raw"])
    assert np.array_equal(out["raw"].X, k9raw.X) and np.array_equal(out["raw"].coords, k9raw.coords)
    rec = prep["record"]
    assert json.loads(json.dumps(rec)) == rec
    assert rec["basis"] and rec["n_contrasts"] >= 2 and 0.0 <= rec["removed_share_y_r8"] <= 1.0
    assert rec["series_days"] == sorted(rec["series_days"]) and len(rec["series_days"]) == len(k7.FULL_PANEL)
    assert not np.array_equal(out["raw_r8"].X, out["raw"].X)
    assert np.linalg.norm(prep["u_perm"]) == pytest.approx(np.linalg.norm(prep["u"]))
    # the gate on the same files gives the same contrasts as the run
    inputs = {r: world[r] for r in k12.gate_files("x1")}
    gr = k12.gate_reference(inputs, t)
    assert gr["n_contrasts"] == rec["n_contrasts"] and gr["verdict"] in ("inside", "not inside")
    # no basis: every arm is raw bit for bit
    k12._PREP.clear()
    real = k12.series_contrasts
    monkeypatch.setattr(k12, "series_contrasts", lambda *a, **kw: real(*a, **kw) | {"u": None})
    with harness.guarded_reads():
        for a in ("raw_r8", "raw_r8_perm"):
            off = k12._build(ctx, a)
            assert np.array_equal(off.X, out["raw"].X) and np.array_equal(off.coords, out["raw"].coords), a
    k12._PREP.clear()
    k9._PREP.clear()


# ----------------------------------------------------------------------------------------------- the reading
def _result(key, scores, clean=True, commit="abc"):
    rows = [{"seed": s, "scores": {a: v[i] for a, v in scores.items()},
             "skills": {a: {"de_score": 0.5, "severity": 0.6} for a in scores},
             "t3": {a: {"severity": {"state": "scored"}} for a in scores}} for i, s in enumerate(SEEDS)]
    return {"task": {"key": key}, "rows": rows, "provenance": {"commit": commit, "tree_clean": clean},
            "factor_record": {"n_contrasts": 9, "removed_share_y_r8": 0.1, "removed_share_y_r8_perm": 0.01},
            "gate": {"verdict": "inside"}}


def _scores(r8=54.0, perm=50.0, raw=52.0):
    n = len(SEEDS)
    return {harness.FLOOR: [50.0] * n, "raw": [raw] * n, "raw_r8": [r8] * n, "raw_r8_perm": [perm] * n}


def _read(d: Path, k9_dir=None, **per_task):
    d.mkdir()
    for name, sc in per_task.items():
        (d / f"{name}.json").write_text(json.dumps(_result(k12.TASKS[name].key, sc)))
    return k12.reading(d, k9_dir if k9_dir is not None else d / "no_k9")


def test_the_reading_supported(tmp_path):
    yes = _read(tmp_path / "yes", x1=_scores(), r12=_scores(), x2=_scores(r8=51.8))
    assert not yes["void"] and yes["supported"] and yes["verdict"].startswith("SUPPORTED")
    txt = k12.text(yes)
    assert "NOT the adoption rule" in txt and "share of y's squared norm removed" in txt


def test_the_reading_closed_by_x2(tmp_path):
    res = _read(tmp_path / "x2", x1=_scores(), r12=_scores(), x2=_scores(r8=51.6))
    assert not res["void"] and not res["supported"] and res["verdict"].startswith("CLOSED")
    assert not res["checks"]["x2: raw_r8 - raw not below -0.3"]


@pytest.mark.parametrize("fold", ["x1", "r12"])
def test_the_reading_closed_by_the_control(tmp_path, fold):
    sc = {n: _scores() for n in k12.TASKS}
    sc[fold] = _scores(perm=53.9)
    res = _read(tmp_path / fold, **sc)
    assert not res["void"] and not res["supported"] and not res["checks"][f"{fold}: raw_r8 - raw_r8_perm"]
    assert res["checks"][f"{fold}: raw_r8 - raw"]


def test_the_reading_needs_eight_seeds_ahead(tmp_path):
    sc = _scores()
    sc["raw_r8"] = [70.0] * 7 + [40.0] * (len(SEEDS) - 7)
    res = _read(tmp_path / "a", x1=sc, r12=_scores(), x2=_scores())
    assert not res["supported"]


def test_the_reading_is_void_on_a_raw_mismatch(tmp_path):
    k9_dir = tmp_path / "k9"
    k9_dir.mkdir()
    for name in k12.TASKS:
        (k9_dir / f"{name}.json").write_text(json.dumps(_result(k9.TASKS[name].key, {"raw": [52.0] * len(SEEDS)})))
    ok = _read(tmp_path / "ok", k9_dir, x1=_scores(), r12=_scores(), x2=_scores())
    assert not ok["void"] and ok["k9_reproduction"] == "reproduced"
    off = _read(tmp_path / "off", k9_dir, x1=_scores(raw=52.01), r12=_scores(), x2=_scores())
    assert off["void"] and any("does not reproduce" in w for w in off["why"])
    assert k12.text(off).startswith("VOID")


def test_the_reading_is_void_on_missing_arm_seed_or_dirty_tree(tmp_path):
    assert _read(tmp_path / "m", x1=_scores(), r12=_scores())["void"]
    sc = _scores()
    del sc["raw_r8_perm"]
    assert _read(tmp_path / "arm", x1=sc, r12=_scores(), x2=_scores())["void"]
    d = tmp_path / "seed"
    d.mkdir()
    for name in k12.TASKS:
        r = _result(k12.TASKS[name].key, _scores())
        if name == "x1":
            r["rows"] = r["rows"][:-1]
        (d / f"{name}.json").write_text(json.dumps(r))
    assert k12.reading(d, tmp_path / "none")["void"]
    d = tmp_path / "dirty"
    d.mkdir()
    for name in k12.TASKS:
        (d / f"{name}.json").write_text(json.dumps(_result(k12.TASKS[name].key, _scores(), clean=name != "x2")))
    assert k12.reading(d, tmp_path / "none")["void"]


# ----------------------------------------------------------------------------------------------- the commands
def test_run_refuses_without_an_inside_gate(tmp_path, monkeypatch):
    called = []
    monkeypatch.setattr(harness, "run_task", lambda *a, **kw: called.append(a) or (_ for _ in ()).throw(
        RuntimeError("run_task reached")))
    runs = tmp_path / "runs"
    args = ["run", "--task", "x1", "--out", str(runs), "--data", str(tmp_path / "nowhere")]
    with pytest.raises(SystemExit, match="does not say"):
        k12.main(args)  # no gate.json at all
    gate = {"references": {"E8.75": {"verdict": "not inside"}, "E9.5": {"verdict": "inside"}}, "decides": "E8.75"}
    (tmp_path / "gate.json").write_text(json.dumps(gate))
    with pytest.raises(SystemExit, match="does not say"):
        k12.main(args)
    other = tmp_path / "elsewhere.json"
    other.write_text(json.dumps({"references": {"E8.75": {"verdict": "inside"}}}))
    with pytest.raises(SystemExit, match="does not say"):
        k12.main(args)  # the parent's gate still says "not inside"
    assert not called and not (runs / "x1.json").exists()
    with pytest.raises(RuntimeError, match="run_task reached"):
        k12.main([*args, "--gate", str(other)])
    assert len(called) == 1


def test_run_gate_and_read_refuse_to_overwrite(tmp_path):
    (tmp_path / "gate.json").write_text(json.dumps({"references": {"E8.75": {"verdict": "inside"}}}))
    runs = tmp_path / "runs"
    runs.mkdir()
    (runs / "x1.json").write_text("{}")
    with pytest.raises(SystemExit, match="exists"):
        k12.main(["run", "--task", "x1", "--out", str(runs), "--data", str(tmp_path / "nowhere")])
    with pytest.raises(SystemExit, match="exists"):
        k12.main(["gate", "--out", str(tmp_path), "--data", str(tmp_path / "nowhere")])
    out = tmp_path / "out"
    out.mkdir()
    (out / "reading.json").write_text("{}")
    with pytest.raises(SystemExit, match="exists"):
        k12.main(["read", "--runs", str(runs), "--out", str(out)])
