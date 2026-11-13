"""The embryo frame alignment with a reflection allowed (vec.core.coupling_rotation(allow_reflection=True),
vec.t2_build.interpolate_clouds(reflect=True), candidate fgw_reflect of try E30, and the target-free gate of diagnostic
N34, experiments/harness/n34_reflect_gate.py): synthetic clouds, boards and stages only."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pytest
from test_t2_build import _left_right_pair, _thinnest, _xcorr_min

from vec import core
from vec import t2_build as tb
from vec.harness import CANDIDATES, FGW_REFLECT_RULE, load_tasks

REPO = Path(__file__).resolve().parents[1]


def _rounded(x):
    """Floats at 8 decimals (a digest that a last-bit difference between BLAS builds does not move), recursively."""
    if isinstance(x, float):
        return round(x, 8)
    if isinstance(x, dict):
        return {k: _rounded(v) for k, v in x.items()}
    if isinstance(x, list | tuple):
        return [_rounded(v) for v in x]
    return x


def _h(a) -> str:
    a = np.asarray(a)
    a = np.ascontiguousarray(np.round(a, 8) if a.dtype.kind == "f" else a)
    return hashlib.sha256(a.dtype.str.encode() + str(a.shape).encode() + a.tobytes()).hexdigest()


def _digest(move_coupling: str, **kw) -> str:
    """One digest of everything interpolate_clouds returns on a fixed synthetic left-right symmetric pair."""
    A, B, Eb, Ea = _left_right_pair(np.random.default_rng(7), n=300)
    r = tb.interpolate_clouds(A, B, Eb, Ea, 0.4, 350, seed=3, target_rms=1.0, n_anchor=150, k=4,
                              move_coupling=move_coupling, with_partners=True, **kw)
    parts = [_h(r.coords), _h(r.pool_coords), _h(r.pool_weights), _h(r.rows), _h(r.stage), _h(r.partners),
             hashlib.sha256(json.dumps(_rounded(r.record), sort_keys=True).encode()).hexdigest()]
    return hashlib.sha256("".join(parts).encode()).hexdigest()


def _rotation_digests(**kw) -> list[str]:
    """coupling_rotation on two fixed cases: a random coupling of two random clouds, and a cloud and its mirror image
    across x coupled by identity (the case where Kabsch's determinant correction acts)."""
    rng = np.random.default_rng(11)
    A, B = rng.normal(size=(40, 3)), rng.normal(size=(50, 3))
    T = rng.random((40, 50))
    T /= T.sum()
    M = A * [-1.0, 1.0, 1.0]
    return [_h(core.coupling_rotation(A, B, T, **kw)), _h(core.coupling_rotation(A, M, np.eye(40) / 40, **kw))]


# Taken at 4bf2c70, before the reflect option existed (the same digests on the code of that commit; the interpolation
# digests equal N32's, taken at 81ab462 on the same case)
PRE_EDIT = {"fgw": "2010611cc184936436cd37b1882b0b05632c816eb5f8be883ccff11afba7faa9",
            "aligned": "3cfce3e9edc30f139af57ab96ba9e7071ac4c4dfa9ba32c22a355f9cec04b0fd"}
PRE_EDIT_ROTATION = ["d97b0929c93f348118d9a365336849d5957a6f8ca58881bb864d9c11131b291d",
                     "63c5c6c3782488e863c67a848ef0c453b87e4b816d09dddaeb07005afa04f52c"]


@pytest.mark.parametrize("move_coupling", sorted(PRE_EDIT))
def test_interpolation_without_reflect_is_unchanged(move_coupling):
    assert _digest(move_coupling) == PRE_EDIT[move_coupling]
    assert _digest(move_coupling, reflect=False) == PRE_EDIT[move_coupling]


def test_coupling_rotation_without_the_keyword_is_unchanged():
    assert _rotation_digests() == PRE_EDIT_ROTATION
    assert _rotation_digests(allow_reflection=False) == PRE_EDIT_ROTATION


def _body(seed: int, n: int = 200) -> np.ndarray:
    """An anisotropic cloud, centred: x thinnest, y longest."""
    A = np.random.default_rng(seed).normal(size=(n, 3)) * [0.5, 2.0, 1.0]
    return A - A.mean(0)


def test_a_mirrored_frame_needs_a_reflection_which_the_proper_rotation_cannot_give():
    A = _body(0)
    B = A * [-1.0, 1.0, 1.0]  # the same cells in a frame mirrored across x
    T = np.eye(len(A)) / len(A)  # each cell coupled to itself
    Rr = core.coupling_rotation(A, B, T, allow_reflection=True)
    Rp = core.coupling_rotation(A, B, T)
    assert np.linalg.det(Rr) == pytest.approx(-1.0) and np.linalg.det(Rp) == pytest.approx(1.0)
    np.testing.assert_allclose(Rr.T @ Rr, np.eye(3), atol=1e-12)
    np.testing.assert_allclose(B @ Rr, A, atol=1e-10)  # aligned exactly
    # the proper rotation keeps two axes and leaves one anti-aligned: whitened cross-correlations 1, 1 and -1
    P = B @ Rp
    X, Y = A - A.mean(0), P - P.mean(0)
    ev = np.linalg.eigvalsh(X.T @ X)
    W = np.linalg.eigh(X.T @ X)[1] @ np.diag(ev ** -0.5) @ np.linalg.eigh(X.T @ X)[1].T
    K = W @ (X.T @ Y) @ W
    np.testing.assert_allclose(np.sort(np.linalg.eigvalsh((K + K.T) / 2)), [-1.0, 1.0, 1.0], atol=1e-8)
    assert _xcorr_min(A, P) == pytest.approx(-1.0, abs=1e-8)
    assert not np.allclose(P, A, atol=1e-3)


@pytest.mark.parametrize("case", ["identity", "random"])
def test_an_unmirrored_pair_gets_the_same_matrix_either_way(case):
    rng = np.random.default_rng(4)
    A = _body(1)
    Q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    Q = Q if np.linalg.det(Q) > 0 else Q @ np.diag([1.0, 1.0, -1.0])
    if case == "identity":
        B, T = A @ Q.T + rng.normal(scale=0.02, size=A.shape), np.eye(len(A)) / len(A)
    else:
        B = A @ Q.T
        T = np.eye(len(A)) + rng.random((len(A), len(A))) * 0.002
        T /= T.sum()
    Rr = core.coupling_rotation(A, B, T, allow_reflection=True)
    Rp = core.coupling_rotation(A, B, T)
    assert np.linalg.det(Rr) == pytest.approx(1.0)
    np.testing.assert_allclose(Rr, Rp, atol=1e-12)


def _mirrored_pair(seed: int, n: int = 300):
    """An earlier cloud (x thinnest) and the same cells in a frame mirrored across x, rows shuffled, slight noise;
    each cell's expression depends on its own signed position, so features pair every cell with itself and the
    coupling is the identity up to the shuffle."""
    rng = np.random.default_rng(seed)
    A = rng.normal(size=(n, 3)) * [0.6, 2.0, 1.0]
    perm = rng.permutation(n)
    B = (A * [-1.0, 1.0, 1.0])[perm] + rng.normal(scale=0.02, size=(n, 3))
    centres = rng.normal(size=(12, 3)) * [0.6, 2.0, 1.0]

    def expr(C):
        return np.exp(-((C[:, None, :] - centres[None]) ** 2).sum(-1))

    Eb = expr(A)
    return A, B, Eb, Eb[perm]


def _run(A, B, Eb, Ea, reflect: bool, seed: int = 2):
    res = tb.interpolate_clouds(A, B, Eb, Ea, 0.5, 400, seed=seed, target_rms=1.0, n_anchor=len(A), k=4,
                                reflect=reflect, with_coupling=True)
    c = res.coupling
    partners = (c["moves"] @ c["after_frame"][c["anchors_after"]]) / c["moves"].sum(1, keepdims=True)
    return res, {"xcorr": _xcorr_min(c["before_frame"][c["anchors_before"]], partners),
                 "thin": _thinnest(res.pool_coords, res.pool_weights) / min(_thinnest(A), _thinnest(B))}


def test_fgw_reflect_keeps_the_moved_cloud_thick_where_fgw_flattens_it():
    A, B, Eb, Ea = _mirrored_pair(0)
    base, mb = _run(A, B, Eb, Ea, False)
    refl, mr = _run(A, B, Eb, Ea, True)
    # the same anchors, FGW coupling, before frame, draws and order: only the frame alignment differs
    for key in ("anchors_before", "anchors_after", "fgw", "moves", "before_frame"):
        assert np.array_equal(base.coupling[key], refl.coupling[key]), key
    np.testing.assert_array_equal(base.rows, refl.rows)
    np.testing.assert_array_equal(base.stage, refl.stage)
    assert np.linalg.det(base.coupling["rotation"]) == pytest.approx(1.0)
    assert np.linalg.det(refl.coupling["rotation"]) == pytest.approx(-1.0)
    assert refl.record["rotation_det"] == -1 and refl.record["reflect"] is True
    assert "rotation_det" not in base.record and "reflect" not in base.record
    # fgw: partners anti-aligned along x and the moved cloud flattened towards the mid-plane; fgw_reflect: neither
    assert mb["xcorr"] < 0 and mb["thin"] < 0.5
    assert mr["xcorr"] > 0.9 and mr["thin"] > 0.9
    assert np.mean(refl.record["mean_anchor_shift"]) < np.mean(base.record["mean_anchor_shift"])


def test_fgw_reflect_on_an_unmirrored_pair_is_fgw_up_to_rounding():
    A, B, Eb, Ea = _mirrored_pair(1)
    B = B * [-1.0, 1.0, 1.0]  # undo the mirror: the same frame, rows shuffled
    base, _ = _run(A, B, Eb, Ea, False)
    refl, _ = _run(A, B, Eb, Ea, True)
    assert refl.record["rotation_det"] == 1
    np.testing.assert_allclose(refl.coupling["rotation"], base.coupling["rotation"], atol=1e-12)
    np.testing.assert_allclose(refl.coords, base.coords, atol=1e-9)
    np.testing.assert_array_equal(refl.rows, base.rows)


# ---------------------------------------------------------------------------------------------------- the candidate
def test_fgw_reflect_is_registered_and_applies_where_fgw_aligned_does():
    c, al = CANDIDATES["fgw_reflect"], CANDIDATES["fgw_aligned"]
    assert c.params["rule"] == FGW_REFLECT_RULE and "reflection" in FGW_REFLECT_RULE
    for f in ("tasks", "needs", "requires", "tuned_on", "externals", "derived", "derived_optional"):
        assert getattr(c, f) == getattr(al, f), f
    tasks = [t for ts in load_tasks().values() for t in ts]
    assert any(al.applies(t) for t in tasks) and not all(al.applies(t) for t in tasks)
    for t in tasks:
        assert c.applies(t) == al.applies(t), t.key
    keys = {t.key: t for t in tasks}
    for key in ("T2e/E7.25+E8.0->E7.5", "T2e/E6.75+E7.5->E7.25"):
        assert c.applies(keys[key])
    assert not c.applies(keys["T2h/<=E9.5->E10.5"])


def test_fgw_reflect_builds_fgw_with_reflect_and_nothing_else(monkeypatch):
    from test_harness import INTERP, _pool, _task

    from vec import harness
    from vec.boards import Board

    calls = []
    real = tb.interpolate_clouds

    def small(*args, **kw):
        calls.append(dict(kw))
        return real(*args, **{**kw, "n_anchor": 150})

    monkeypatch.setattr(harness.t2_build, "interpolate_clouds", small)
    from conftest import GENES

    rng = np.random.default_rng(21)
    base = rng.gamma(2.0, 2.0, size=(3, len(GENES)))
    earlier = _pool(1500, rng, base * 0.8, [0.6, 0.3, 0.1], scale=0.7)
    before = _pool(2500, rng, base, [0.5, 0.3, 0.2])
    after = _pool(4000, rng, base * 1.5, [0.3, 0.4, 0.3], scale=1.8, stretch=(1.0, 1.4, 1.0))
    board = Board(key="pseudo:synthetic_embryo", task="T2", genes=list(GENES), needs_coords=True, setting="embryo",
                  min_cells=100, max_cells=1500)
    task = _task(**INTERP, key="syn_embryo_reflect",
                 inputs={"before": "T2/E7.25.h5ad", "after": "T2/E8.0.h5ad", "earlier": "T2/E6.75.h5ad"})
    ctx = harness.BuildContext(task=harness.TaskView.of(task), board=board,
                               inputs={"before": before, "after": after, "earlier": earlier}, seed=7)
    fgw = CANDIDATES["fgw"].build(ctx)
    refl = CANDIDATES["fgw_reflect"].build(ctx)
    assert len(calls) == 2
    assert calls[0].pop("reflect") is False and calls[1].pop("reflect") is True
    assert calls[0] == calls[1]
    np.testing.assert_array_equal(refl.X, fgw.X)  # the same cells and expression rows
    assert harness.t2_build.rms_radius(refl.coords) == pytest.approx(harness.t2_build.rms_radius(fgw.coords),
                                                                     rel=0.01)


# ---------------------------------------------------------------------------------------------------- the gate N34
def _gate():
    spec = importlib.util.spec_from_file_location("n34_reflect_gate_test",
                                                  REPO / "experiments" / "harness" / "n34_reflect_gate.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def test_gate_constants_are_the_fixed_ones():
    g = _gate()
    assert (g.XCORR_MIN, g.THIN_RATIO_MIN, g.DRAWS_MIN) == (0.6, 0.9, 8)
    assert g.dc.TEST_PAIR == "test/E7.5+E8.0->E7.75" and set(g.BOARD_OF_PAIR) == {k for k, *_ in g.dc.PAIRS}
    assert [g.dc.BUILD_SEED_OFFSET + s for s in g.dc.SEEDS] == list(range(1_000_000, 1_000_010))


def test_gate_draw_rule_at_its_edges():
    g = _gate()
    assert g.draw_gate({"xcorr_min": 0.6, "thin_ratio": 0.9})["passes"]
    assert not g.draw_gate({"xcorr_min": 0.599, "thin_ratio": 0.95})["passes"]
    r = g.draw_gate({"xcorr_min": 0.95, "thin_ratio": 0.899})
    assert not r["passes"] and r["checks"] == {"xcorr": True, "thin": False}
    assert g.det_sign(np.diag([-1.0, 1.0, 1.0])) == -1 and g.det_sign(np.eye(3)) == 1


def _m(xcorr, thin, det=-1, shift=(0.1, 0.1), label=0.8):
    return {"xcorr_min": xcorr, "xcorr_eigs": [xcorr, 0.9, 1.0], "thin_ratio": thin, "label_agreement": label,
            "mean_anchor_shift": list(shift), "rotation_deg": 170.0, "shape_at_t": [0.8, 0.5, 0.3],
            "nn15_kept": 0.5, "det": det}


def _table(g, n_pass: int, n: int = 10) -> dict:
    draws = []
    for s in range(n):
        refl = _m(0.9, 0.95) if s < n_pass else _m(0.9, 0.5)
        draws.append({"seed": s, "build_seed": 1_000_000 + s, "fgw": _m(-0.4, 0.6, det=1), "fgw_reflect": refl,
                      "gate": g.draw_gate(refl)})
    return {"genes_board": "pseudo", "n_pass": sum(d["gate"]["passes"] for d in draws), "draws": draws}


def test_gate_verdict_reads_the_test_pair_only():
    g = _gate()
    test, bt = g.dc.TEST_PAIR, [k for k, *_ in g.dc.PAIRS if k != g.dc.TEST_PAIR]
    pairs = {bt[0]: _table(g, 0), bt[1]: _table(g, 2), test: _table(g, 8)}
    assert pairs[test]["n_pass"] == 8 and g.verdict(pairs, True, test) == "open"  # backtests failing decide nothing
    pairs[test] = _table(g, 7)
    assert g.verdict(pairs, True, test) == "prune"
    pairs = {bt[0]: _table(g, 10), bt[1]: _table(g, 10), test: _table(g, 7)}
    assert g.verdict(pairs, True, test) == "prune"
    assert g.verdict(pairs, False, test) == "partial run: no verdict"
    res = {"verdict": "prune", "seeds_label": "1,000,000 to 1,000,009", "pairs": pairs, "inputs": {"x.h5ad": "abc"},
           "gate": {"XCORR_MIN": 0.6, "THIN_RATIO_MIN": 0.9, "DRAWS_MIN": 8, "TEST_PAIR": test},
           "stamp": {"commit": "abc", "tree_clean": True, "utc": "2026-11-13T00:00:00Z"}}
    txt = g.report_text(res)
    assert txt.endswith("verdict: prune\n") and f"{test} (pseudo's panel; decides): 7 of 10 draws pass" in txt
    assert "reported only" in txt and "reflections chosen in 10 of 10 draws" in txt


def test_gate_refuses_to_overwrite(tmp_path):
    g = _gate()
    (tmp_path / "n34_reflect_gate.txt").write_text("earlier\n")
    with pytest.raises(SystemExit, match="already written"):
        g.main(["--out", str(tmp_path)])
    assert (tmp_path / "n34_reflect_gate.txt").read_text() == "earlier\n"
