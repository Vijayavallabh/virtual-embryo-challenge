"""The embryo coupling started from the aligned coupling (move_coupling="fgw_warm", candidate fgw_warm, diagnostic N32
of experiments/harness/TRIES.md and its script experiments/harness/n32_fgw_warm_gate.py): synthetic clouds, boards and
stages only."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import re
import sys
from pathlib import Path

import numpy as np
import pytest

from vec import core
from vec import t2_build as tb
from vec.harness import BUILD_SEED_OFFSET, CANDIDATES, FGW_WARM_RULE, load_tasks

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


def _digest(move_coupling: str) -> str:
    """One digest of everything interpolate_clouds returns on a fixed synthetic left-right symmetric pair."""
    from test_t2_build import _left_right_pair

    A, B, Eb, Ea = _left_right_pair(np.random.default_rng(7), n=300)
    r = tb.interpolate_clouds(A, B, Eb, Ea, 0.4, 350, seed=3, target_rms=1.0, n_anchor=150, k=4,
                              move_coupling=move_coupling, with_partners=True)
    parts = [_h(r.coords), _h(r.pool_coords), _h(r.pool_weights), _h(r.rows), _h(r.stage), _h(r.partners),
             hashlib.sha256(json.dumps(_rounded(r.record), sort_keys=True).encode()).hexdigest()]
    return hashlib.sha256("".join(parts).encode()).hexdigest()


# Taken at 81ab462, before fgw_warm existed (the same _digest on the code of that commit)
PRE_EDIT = {"fgw": "2010611cc184936436cd37b1882b0b05632c816eb5f8be883ccff11afba7faa9",
            "aligned": "3cfce3e9edc30f139af57ab96ba9e7071ac4c4dfa9ba32c22a355f9cec04b0fd"}


@pytest.mark.parametrize("move_coupling", sorted(PRE_EDIT))
def test_fgw_and_aligned_outputs_are_unchanged(move_coupling):
    assert _digest(move_coupling) == PRE_EDIT[move_coupling]


def _mirror_body(seed: int, n: int = 60):
    """One side of a body and its exact mirror image across x, features the same on both sides; B is a noisy copy of
    A in the same frame, rows shuffled: the coupling to the copy and the one to its mirror image cost the same."""
    rng = np.random.default_rng(seed)
    side = rng.normal(size=(n, 3)) * [1.0, 2.0, 0.7]
    side[:, 0] = np.abs(side[:, 0]) + 0.2
    A = np.vstack([side, side * [-1.0, 1.0, 1.0]])
    F = np.c_[np.abs(A[:, 0]), A[:, 1], A[:, 2]]
    perm = rng.permutation(2 * n)
    return A, F, A[perm] + rng.normal(scale=0.05, size=(2 * n, 3)), F[perm], np.argsort(perm)


def test_fgw_coupling_without_a_start_is_the_call_before_g0_existed():
    import ot

    A, F, B, FB, _ = _mirror_body(3, n=20)
    M, CA, CB = ot.dist(F, FB), ot.dist(A, A, metric="euclidean"), ot.dist(B, B, metric="euclidean")
    for C in (M, CA, CB):
        C /= C.max()
    p = np.full(40, 1 / 40)
    old = ot.gromov.fused_gromov_wasserstein(M, CA, CB, p, p, loss_fun="square_loss", alpha=0.5, max_iter=200,
                                             numItermaxEmd=100_000)
    assert np.array_equal(core.fgw_coupling(F, FB, A, B, 0.5, 200, 100_000), np.asarray(old, dtype=np.float64))
    assert np.array_equal(core.fgw_coupling(F, FB, A, B, 0.5, 200, 100_000, G0=None), old)


def test_fgw_objective_is_the_solvers_and_a_warm_start_keeps_the_side_the_start_is_on():
    import ot
    from test_t2_build import _xcorr_min

    A, F, B, FB, own = _mirror_body(0)
    cold = core.fgw_coupling(F, FB, A, B, 0.5, 200, 100_000)
    start = core.aligned_coupling(F, FB, A, B, 0.5, 100_000)
    warm = core.fgw_coupling(F, FB, A, B, 0.5, 200, 100_000, G0=start)
    # POT's own value of its objective at its solution
    M, CA, CB = ot.dist(F, FB), ot.dist(A, A, metric="euclidean"), ot.dist(B, B, metric="euclidean")
    for C in (M, CA, CB):
        C /= C.max()
    p = np.full(len(A), 1 / len(A))
    _, log = ot.gromov.fused_gromov_wasserstein(M, CA, CB, p, p, loss_fun="square_loss", alpha=0.5, max_iter=200,
                                                numItermaxEmd=100_000, G0=start, log=True)
    assert core.fgw_objective(F, FB, A, B, 0.5, warm) == pytest.approx(float(log["fgw_dist"]), rel=1e-9, abs=1e-15)
    # from the independent start FGW lands on the mirror image; from the aligned start it stays on each cell's own copy
    assert (cold.argmax(1) == own).mean() == 0.0 and (warm.argmax(1) == own).all()
    obj = [core.fgw_objective(F, FB, A, B, 0.5, C) for C in (cold, start, warm)]
    assert obj[2] <= obj[1] + 1e-12 and obj[2] == pytest.approx(obj[0], rel=0.05)  # the two optima cost the same
    assert _xcorr_min(A, B[cold.argmax(1)]) < 0 < 0.9 < _xcorr_min(A, B[warm.argmax(1)])
    np.testing.assert_allclose(warm.sum(1), 1 / len(A), atol=1e-9)
    with pytest.raises(ValueError, match="coupling shape"):
        core.fgw_objective(F, FB, A, B, 0.5, warm[:-1])


def _mirror_pair(seed: int, n: int = 200):
    """Two clouds of one mirror-symmetric body (x thinnest, expression a bump in |x|), the after cloud the before one
    stretched along y plus noise, rows shuffled: FGW can couple each side to either side of the other cloud."""
    rng = np.random.default_rng(seed)
    side = rng.normal(size=(n, 3)) * [0.6, 2.0, 1.0]
    side[:, 0] = np.abs(side[:, 0]) + 0.3
    A = np.vstack([side, side * [-1.0, 1.0, 1.0]])
    B = (A * [1.0, 1.3, 1.0])[rng.permutation(2 * n)] + rng.normal(scale=0.05, size=(2 * n, 3))
    centres = rng.normal(size=(12, 3)) * [0.5, 1.0, 1.0]
    centres[:, 0] = np.abs(centres[:, 0]) + 0.3

    def expr(C, length):
        f = np.c_[np.abs(C[:, 0]), C[:, 1] / length, C[:, 2]]
        return np.exp(-((f[:, None, :] - centres[None]) ** 2).sum(-1))

    return A, B, expr(A, 2.0), expr(B, 2.6)


def _run(A, B, Eb, Ea, mc, seed, n_anchor=300):
    from test_t2_build import _thinnest, _xcorr_min

    res = tb.interpolate_clouds(A, B, Eb, Ea, 0.5, 400, seed=seed, target_rms=1.0, n_anchor=n_anchor, k=4,
                                move_coupling=mc, with_coupling=True)
    c = res.coupling
    partners = (c["moves"] @ c["after_frame"][c["anchors_after"]]) / c["moves"].sum(1, keepdims=True)
    return res, {"xcorr": _xcorr_min(c["before_frame"][c["anchors_before"]], partners),
                 "thin": _thinnest(res.pool_coords, res.pool_weights) / min(_thinnest(A), _thinnest(B)),
                 "shift": float(np.mean(res.record["mean_anchor_shift"]))}


def test_fgw_warm_shares_everything_but_the_moves_with_fgw():
    A, B, Eb, Ea = _mirror_pair(0)
    base, _ = _run(A, B, Eb, Ea, "fgw", 2)
    warm, _ = _run(A, B, Eb, Ea, "fgw_warm", 2)
    for key in ("anchors_before", "anchors_after", "fgw", "rotation", "before_frame", "after_frame"):
        assert np.array_equal(base.coupling[key], warm.coupling[key]), key
    assert np.array_equal(base.rows, warm.rows) and np.array_equal(base.stage, warm.stage)
    assert warm.record["rotation_deg"] == base.record["rotation_deg"]  # the first solution's rotation is kept
    assert warm.record["coupling_nonzero"] == base.record["coupling_nonzero"]
    assert warm.record["move_coupling"] == "fgw_warm" and "fgw_objective" not in base.record
    obj = warm.record["fgw_objective"]
    assert set(obj) == {"fgw", "aligned_start", "fgw_warm"} and obj["fgw_warm"] <= obj["aligned_start"] + 1e-12
    assert not np.array_equal(base.coupling["moves"], warm.coupling["moves"])
    np.testing.assert_allclose(warm.coupling["moves"].sum(1), 1 / 300, atol=1e-9)


def test_on_a_mirrored_draw_fgw_warm_keeps_the_sides_and_the_thickness():
    A, B, Eb, Ea = _mirror_pair(0)
    out = {mc: _run(A, B, Eb, Ea, mc, 2)[1] for mc in ("fgw", "aligned", "fgw_warm")}
    assert out["fgw"]["xcorr"] < -0.5 and out["fgw"]["thin"] < 0.5  # FGW mirrored this draw: the moved cloud flattens
    assert out["fgw_warm"]["xcorr"] >= 0.6 and out["fgw_warm"]["thin"] >= 0.9  # N32's first two checks hold
    # ... while its moves are shorter than the mirrored ones, which cross the body: N32's third check would fail here
    assert out["fgw_warm"]["shift"] < 0.9 * out["fgw"]["shift"]
    assert out["fgw_warm"]["shift"] > out["aligned"]["shift"]  # longer than the aligned coupling's own moves


def test_on_an_unmirrored_draw_fgw_warm_keeps_fgws_long_moves():
    A, B, Eb, Ea = _mirror_pair(0)
    out = {mc: _run(A, B, Eb, Ea, mc, 1)[1] for mc in ("fgw", "fgw_warm")}
    assert out["fgw"]["xcorr"] > 0.3  # FGW already kept the sides this draw
    assert out["fgw_warm"]["shift"] >= 0.9 * out["fgw"]["shift"] and out["fgw_warm"]["xcorr"] >= out["fgw"]["xcorr"]


def test_on_the_left_right_pair_with_random_sides_the_warm_start_does_not_hold():
    """tests/test_t2_build.py's left-right pair (two different clouds, x thin): FGW's partners are random across x, not
    mirrored, the objective is nearly flat along x, and the warm start drifts back to FGW-like partners. Its long moves
    stay, but its partners do not reach the aligned coupling's: what N32's first check reads on real data."""
    from test_t2_build import _left_right_pair

    A, B, Eb, Ea = _left_right_pair(np.random.default_rng(6))
    out = {mc: _run(A, B, Eb, Ea, mc, 2)[1] for mc in ("fgw", "aligned", "fgw_warm")}
    assert out["fgw"]["xcorr"] < 0.2
    assert out["fgw_warm"]["shift"] >= 0.9 * out["fgw"]["shift"]
    assert out["fgw_warm"]["xcorr"] < out["aligned"]["xcorr"] - 0.2


def test_fgw_warm_is_registered_beside_fgw_aligned_and_refuses_what_it_refuses():
    c, twin = CANDIDATES["fgw_warm"], CANDIDATES["fgw_aligned"]
    assert c.tasks == twin.tasks and c.needs == twin.needs and c.requires is twin.requires
    assert not c.tuned_on and not c.externals and not c.derived
    assert c.params["rule"] == FGW_WARM_RULE and c.params["alpha"] == tb.INTERP_ALPHA
    tasks = [t for ts in load_tasks().values() for t in ts]
    assert [c.applies(t) for t in tasks] == [twin.applies(t) for t in tasks]
    assert any(twin.applies(t) for t in tasks) and not all(twin.applies(t) for t in tasks)
    with pytest.raises(ValueError, match="move_coupling"):
        A, B, Eb, Ea = _mirror_pair(0, n=20)
        tb.interpolate_clouds(A, B, Eb, Ea, 0.5, 10, seed=0, target_rms=1.0, move_coupling="fgw_cold")


# ------------------------------------------------------------------------- the gate script, on synthetic stages only
@pytest.fixture(scope="module")
def n32():
    spec = importlib.util.spec_from_file_location("harness_n32_fgw_warm_gate",
                                                  REPO / "experiments" / "harness" / "n32_fgw_warm_gate.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def _patch(n32, monkeypatch):
    """Synthetic boards and stages (tests/test_diag_emb_coupling.py's) in place of boards/ and data/; no git, no file
    hash, no stamp. Returns the load_stage calls."""
    from test_diag_emb_coupling import ALL_GENES, TEST_GENES, VAL_GENES, _stage

    from vec import harness, provenance
    from vec.boards import Board

    boards = {"T2:embryo:val_interp": Board(key="T2:embryo:val_interp", task="T2", genes=VAL_GENES, needs_coords=True,
                                            split="val_interp"),
              "T2:embryo:test_interp": Board(key="T2:embryo:test_interp", task="T2", genes=TEST_GENES,
                                             needs_coords=True, split="test_interp", min_cells=777, max_cells=4321)}
    calls = []

    def load_stage(path, genes=None, drop_duplicates=False):
        name = Path(path).name.removesuffix(".h5ad")
        calls.append((name, list(genes), drop_duplicates))
        return _stage(name)[:, genes].copy()

    assert ALL_GENES
    monkeypatch.setattr(n32.dc, "load_boards", lambda: boards)
    monkeypatch.setattr(n32.dc, "load_stage", load_stage)
    monkeypatch.setattr(harness, "_sha256", lambda p: "0" * 64)
    monkeypatch.setattr(provenance, "git_state", lambda *a, **k: ("cafebabe", []))
    monkeypatch.setattr(provenance, "write_stamp", lambda p, *a, **k: p)
    return calls, VAL_GENES, TEST_GENES


def test_n32_reads_the_pairs_inputs_only_and_writes_summaries(n32, tmp_path, monkeypatch, capsys):
    calls, val, test = _patch(n32, monkeypatch)
    monkeypatch.setattr(n32.dc, "SEEDS", (0, 1))  # a full run in two draws
    out = tmp_path / "n32"
    n32.main(["--out", str(out), "--seeds", "0", "1"])
    printed = capsys.readouterr().out
    # each stage loaded only as an input of a pair it brackets, on that pair's own board's panel
    want = [("E6.75", val), ("E7.5", val), ("E7.25", val), ("E8.0", val), ("E7.5", test), ("E8.0", test)]
    assert [(name, genes) for name, genes, _ in calls] == want and all(d for _, _, d in calls)
    rep = json.loads((out / "n32_fgw_warm_gate.json").read_text())
    txt = (out / "n32_fgw_warm_gate.txt").read_text()
    assert rep["verdict"] in ("open", "prune") and txt.splitlines()[-1] == f"verdict: {rep['verdict']}"
    assert printed.endswith(txt) and set(rep["pairs"]) == set(n32.BOARD_OF_PAIR)
    for p in rep["pairs"].values():
        assert [d["build_seed"] for d in p["draws"]] == [BUILD_SEED_OFFSET, BUILD_SEED_OFFSET + 1]
        assert p["n_pass"] == sum(d["gate"]["passes"] for d in p["draws"])
        for d in p["draws"]:
            assert set(d["objective"]) == {"fgw", "aligned_start", "fgw_warm"}
            assert d["gate"] == n32.draw_gate(d["fgw"], d["fgw_warm"])
    assert rep["gate"]["XCORR_MIN"] == 0.6 and rep["gate"]["DRAWS_MIN"] == 8
    # no cell count, coupling count or board limit anywhere: whole numbers are seeds only
    ints = set()

    def walk(x):
        if isinstance(x, dict):
            assert not {"cells", "n_obs", "coupling_nonzero", "n_before", "n_after", "min_cells", "max_cells"} & set(x)
            for v in x.values():
                walk(v)
        elif isinstance(x, list):
            for v in x:
                walk(v)
        elif isinstance(x, int) and not isinstance(x, bool):
            ints.add(x)

    walk(rep)
    assert ints <= {*range(11), BUILD_SEED_OFFSET, BUILD_SEED_OFFSET + 1}  # seeds, draws passing, the gate's 8
    from test_diag_emb_coupling import SIZES

    forbidden = {777, 4321, *SIZES.values(), n32.dc.N_DRAW, n32.dc.N_DRAW // 2}
    numbers = {int(x) for x in re.findall(r"\d+", re.sub(r"(?<=\d),(?=\d{3})", "", txt + printed))}
    assert not numbers & forbidden
    # it runs once
    with pytest.raises(SystemExit, match="already"):
        n32.main(["--out", str(out), "--seeds", "0"])


def test_n32_gate_and_verdict(n32):
    def m(xcorr, thin, shift):
        return {"xcorr_min": xcorr, "thin_ratio": thin, "mean_anchor_shift": [shift, shift]}

    assert n32.draw_gate(m(-0.3, 0.6, 0.6), m(0.6, 0.9, 0.54))["passes"]
    assert not n32.draw_gate(m(-0.3, 0.6, 0.6), m(0.59, 0.95, 0.6))["checks"]["xcorr"]
    assert not n32.draw_gate(m(-0.3, 0.6, 0.6), m(0.8, 0.89, 0.6))["checks"]["thin"]
    g = n32.draw_gate(m(-0.3, 0.6, 0.6), m(0.8, 0.95, 0.53))
    assert not g["checks"]["shift"] and not g["passes"] and g["shift_ratio"] == pytest.approx(0.53 / 0.6)
    pairs = {k: {"n_pass": 8} for k in n32.BOARD_OF_PAIR}
    assert n32.verdict(pairs, True) == "open" and n32.verdict(pairs, False) == "partial run: no verdict"
    assert n32.verdict({**pairs, "test/E7.5+E8.0->E7.75": {"n_pass": 7}}, True) == "prune"
