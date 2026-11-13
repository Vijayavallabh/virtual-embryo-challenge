"""The polar move (campaign branch campaign/emb-polar; experiments/research_2026_11_runs/REPORT.md recipe R10):
vec.t2_build.polar_moves, interpolate_clouds(..., moves="polar"), the two harness arms and the runner
experiments/harness/n41_polar_arms.py. Synthetic clouds only; no data file, committed result or scorer is read.

The routine interpolate_clouds computes its own FGW coupling and cannot take a fixed one, so the tests that need a
fixed pairing (the rotated copy, the stretch, the mirror, the degenerate patches) call the move step itself,
`polar_moves(Z, anchors, shifts, near, tau)`, with every cell an anchor, `shifts` the identity pairing's displacements
(moved copy minus cloud) and `near` the cells' k nearest anchors (vec.core.nearest_anchors), as interpolate_clouds
calls it for each stage."""

from __future__ import annotations

import importlib.util
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from vec import harness
from vec import t2_build as tb
from vec.core import nearest_anchors

REPO = Path(__file__).resolve().parents[1]
HARNESS = REPO / "experiments" / "harness"
BASE = "51b0672"  # campaign/emb-trend2, the branch's base
K = tb.INTERP_K


def _rot_z(a):
    return np.array([[np.cos(a), -np.sin(a), 0.0], [np.sin(a), np.cos(a), 0.0], [0.0, 0.0, 1.0]])


def _moved(A: np.ndarray, M: np.ndarray, tau: float, k: int = K):
    """The move step on cloud A paired with its copy A M^T (identity pairing, every cell an anchor): the polar move,
    its record and the entry's move, computed as interpolate_clouds computes it."""
    D = A @ M.T - A
    near = nearest_anchors(A, A, k)[1]
    pos, rec = tb.polar_moves(A, A, D, near, tau)
    return pos, rec, A + tau * D[near].mean(1)


def _pair(seed: int = 3, n: int = 400):
    rng = np.random.default_rng(seed)
    A = rng.normal(size=(n, 3)) * [4.0, 2.0, 1.0] + 20.0
    E = rng.gamma(2.0, 1.0, size=(n, 12))
    perm = rng.permutation(n)
    B, EB = (1.5 * (A - A.mean(0)) @ _rot_z(0.9).T - 30.0)[perm], E[perm]
    return A, B, E, EB


KW = dict(seed=1, target_rms=10.0, n_anchor=200)


# ------------------------------------------------------------------------------------------ 1. moves unset
def _base_module(tmp_path: Path):
    """vec.t2_build as it stands on the base commit, loaded from git (read-only)."""
    env = {**os.environ, "GIT_OPTIONAL_LOCKS": "0"}
    try:
        src = subprocess.run(["git", "show", f"{BASE}:src/vec/t2_build.py"], cwd=REPO, env=env, check=True,
                             capture_output=True, text=True).stdout
    except (OSError, subprocess.CalledProcessError):
        pytest.skip(f"the base commit {BASE} is not reachable here")
    path = tmp_path / "t2_build_base.py"
    path.write_text(src)
    spec = importlib.util.spec_from_file_location("t2_build_base", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def _same(a, b) -> None:
    for f in ("stage", "rows", "coords", "pool_coords", "pool_weights"):
        assert np.array_equal(getattr(a, f), getattr(b, f)), f
    assert json.dumps(a.record, sort_keys=True, default=str) == json.dumps(b.record, sort_keys=True, default=str)


@pytest.mark.parametrize("move_coupling", tb.MOVE_COUPLINGS)
def test_moves_unset_is_the_entry_bit_for_bit(tmp_path, move_coupling):
    A, B, E, EB = _pair()
    base = _base_module(tmp_path)
    old = base.interpolate_clouds(A, B, E, EB, 0.4, 300, move_coupling=move_coupling, with_partners=True, **KW)
    new = tb.interpolate_clouds(A, B, E, EB, 0.4, 300, move_coupling=move_coupling, with_partners=True, **KW)
    mean = tb.interpolate_clouds(A, B, E, EB, 0.4, 300, move_coupling=move_coupling, with_partners=True,
                                 moves="mean", **KW)
    _same(old, new)
    _same(new, mean)
    assert np.array_equal(old.partners, new.partners)
    assert "polar" not in new.record and "moves" not in new.record
    polar = tb.interpolate_clouds(A, B, E, EB, 0.4, 300, move_coupling=move_coupling, moves="polar", **KW)
    assert np.array_equal(polar.rows, new.rows) and np.array_equal(polar.stage, new.stage)  # the same draws
    assert not np.array_equal(polar.coords, new.coords)  # the idea is not idle here
    with pytest.raises(ValueError, match="moves"):
        tb.interpolate_clouds(A, B, E, EB, 0.4, 30, moves="affine", **KW)


def test_with_the_identity_map_the_polar_move_is_the_entrys_move():
    rng = np.random.default_rng(0)
    A = rng.normal(size=(200, 3))
    shift = np.array([0.7, -1.2, 0.3])
    near = nearest_anchors(A, A, K)[1]
    D = np.broadcast_to(shift, A.shape).copy()  # A (P - p) = Q - q with A = I
    pos, rec = tb.polar_moves(A, A, D, near, 0.35)
    assert rec["share"] == 1.0 and rec["median_rotation_deg"] == pytest.approx(0.0, abs=1e-6)
    assert np.allclose(pos, A + 0.35 * D[near].mean(1), rtol=0, atol=1e-12)


# ------------------------------------------------------------------------------------------ 2. a rigid rotation
def test_a_rigid_rotation_at_half_way():
    """A cloud and its copy rotated by 1.6 rad about the z axis through its centre, identity pairing, t = 1/2.

    Many patches (400 cells): the entry's move is close to the linear interpolation and shrinks the RMS radius about the
    centre by about cos(angle / 2). The polar move shrinks it too, a little less: each patch's centroid still moves
    along the chord, c(1/2) = (p + R p) / 2; only the cell's offset from its patch centroid is rotated. The brief's
    expectation that the polar move keeps the whole cloud's radius holds only where one patch is the whole rotating
    part (below).

    One patch (8 cells, each one's 8 nearest anchors are all of them): the polar move lands on the half-way rotation
    and keeps the radius to 1e-9; the entry's mean displacement is zero, so the entry leaves the cells where they
    are."""
    rng = np.random.default_rng(0)
    angle = 1.6
    A = rng.normal(size=(400, 3)) * [3.0, 3.0, 1.0]
    A -= A.mean(0)
    pos, rec, entry = _moved(A, _rot_z(angle), 0.5)
    r0 = tb.rms_radius(A)
    assert tb.rms_radius(entry) / r0 == pytest.approx(np.cos(angle / 2), abs=0.03)
    assert tb.rms_radius(entry) / r0 < tb.rms_radius(pos) / r0 < 0.8
    assert rec["share"] == 1.0 and rec["median_rotation_deg"] == pytest.approx(np.degrees(angle), abs=1e-6)
    assert rec["median_stretch_ratio"] == pytest.approx(1.0, abs=1e-9)
    P = A[:K]
    P = P - P.mean(0)
    pos, rec, entry = _moved(P, _rot_z(angle), 0.5)
    assert tb.rms_radius(pos) / tb.rms_radius(P) == pytest.approx(1.0, abs=1e-9)
    assert np.abs(pos - P @ _rot_z(angle / 2).T).max() < 1e-9
    assert np.allclose(entry, P, rtol=0, atol=1e-12)
    assert rec["share"] == 1.0 and rec["median_rotation_deg"] == pytest.approx(np.degrees(angle), abs=1e-6)


# ------------------------------------------------------------------------------------------ 3. a pure stretch
def test_a_pure_stretch_is_the_linear_interpolation():
    rng = np.random.default_rng(2)
    A = rng.normal(size=(300, 3))
    S = np.diag([2.0, 1.0, 0.5])
    for tau in (0.3, 0.5, 1.0):
        pos, rec, entry = _moved(A, S, tau)
        assert np.allclose(pos, (1 - tau) * A + tau * A @ S.T, rtol=0, atol=1e-9)
        assert rec["share"] == 1.0 and rec["median_rotation_deg"] == pytest.approx(0.0, abs=1e-6)
        assert rec["median_stretch_ratio"] == pytest.approx(4.0, rel=1e-9)
        assert not np.allclose(pos, entry)  # the entry's local translation is not the stretch


# ------------------------------------------------------------------------------------------ 4. a mirror
def test_a_mirrored_patch_keeps_the_entrys_move_exactly_and_lowers_the_share():
    rng = np.random.default_rng(4)
    A = rng.normal(size=(300, 3))
    mirror = np.diag([1.0, 1.0, -1.0])
    pos, rec, entry = _moved(A, mirror, 0.5)
    assert np.array_equal(pos, entry) and rec["share"] == 0.0
    assert rec["median_rotation_deg"] is None and rec["median_stretch_ratio"] is None
    # two clusters far apart: one rotates, one is mirrored; the share falls to the rotating half
    B = rng.normal(size=(300, 3)) + 100.0
    C = np.vstack([A, B])
    D = np.vstack([A @ _rot_z(0.5).T - A, B @ mirror.T - B])
    near = nearest_anchors(C, C, K)[1]
    pos, rec = tb.polar_moves(C, C, D, near, 0.5)
    assert rec["share"] == pytest.approx(0.5)
    assert np.array_equal(pos[300:], C[300:] + 0.5 * D[near[300:]].mean(1))
    _, rot = tb.polar_moves(A, A, A @ _rot_z(0.5).T - A, nearest_anchors(A, A, K)[1], 0.5)
    assert rot["share"] == 1.0 > rec["share"]


# ------------------------------------------------------------------------------------------ 5. degenerate patches
def test_degenerate_patches_fall_back_without_nan():
    rng = np.random.default_rng(5)
    flat = np.c_[rng.normal(size=(200, 2)), np.zeros(200)]  # coplanar anchors
    pos, rec, entry = _moved(flat, _rot_z(0.7), 0.5)
    assert rec["share"] == 0.0 and np.isfinite(pos).all() and np.array_equal(pos, entry)
    rep = np.repeat(rng.normal(size=(3, 3)), 10, axis=0)  # repeated anchors: three places, ten anchors each
    cells = rng.normal(size=(50, 3))
    D = rng.normal(size=rep.shape)
    near = nearest_anchors(cells, rep, K)[1]
    pos, rec = tb.polar_moves(cells, rep, D, near, 0.5)
    assert rec["share"] == 0.0 and np.isfinite(pos).all()
    assert np.array_equal(pos, cells + 0.5 * D[near].mean(1))
    # a non-finite displacement: the cells whose patch holds it keep the entry's move, the others stay finite
    A = rng.normal(size=(200, 3))
    D = A @ _rot_z(0.4).T - A
    D[7] = np.inf
    near = nearest_anchors(A, A, K)[1]
    pos, rec = tb.polar_moves(A, A, D, near, 0.5)
    hit = (near == 7).any(1)
    assert np.array_equal(pos[hit], (A + 0.5 * D[near].mean(1))[hit])
    assert np.isfinite(pos[~hit]).all() and rec["share"] == pytest.approx(float((~hit).mean()))
    # fewer anchors than three dimensions need: every patch is degenerate
    pos, rec = tb.polar_moves(A, A[:3], D[:3] * 0 + 1.0, nearest_anchors(A, A[:3], K)[1], 0.5)
    assert rec["share"] == 0.0 and np.isfinite(pos).all()


# ------------------------------------------------------------------------------------------ 6. reproducible
def test_the_polar_routine_is_reproducible_and_records_its_shares():
    A, B, E, EB = _pair()
    a = tb.interpolate_clouds(A, B, E, EB, 0.5, 300, moves="polar", **KW)
    b = tb.interpolate_clouds(A, B, E, EB, 0.5, 300, moves="polar", **KW)
    _same(a, b)
    rec = a.record["polar"]
    assert a.record["moves"] == "polar" and set(rec) == {"share", "median_rotation_deg", "median_stretch_ratio"}
    assert all(len(v) == 2 for v in rec.values()) and all(0.0 <= s <= 1.0 for s in rec["share"])
    assert tb.rms_radius(a.coords) == pytest.approx(10.0)
    other = tb.interpolate_clouds(A, B, E, EB, 0.5, 300, moves="polar", **{**KW, "seed": 2})
    assert not np.array_equal(a.coords, other.coords)


# ------------------------------------------------------------------------------------------ the harness arms
def test_the_two_arms_are_registered_for_t2_and_call_fgw_with_the_polar_move(monkeypatch):
    seen = {}
    monkeypatch.setattr(harness, "_fgw", lambda ctx, **kw: seen.setdefault(ctx, kw))
    for name, kw in (("fgw_polar", {"moves": "polar"}),
                     ("fgw_aligned_polar_series", {"move_coupling": "aligned", "expression": "series",
                                                   "moves": "polar"})):
        c = harness.CANDIDATES[name]
        assert c.tasks == frozenset({"T2"}) and c.name == name and c.params["moves"] == "polar"
        assert c.build(name) == kw
    assert harness.CANDIDATES["fgw_polar"].requires is harness.CANDIDATES["fgw"].requires
    assert harness.CANDIDATES["fgw_aligned_polar_series"].requires is \
        harness.CANDIDATES["fgw_aligned_series"].requires


# ------------------------------------------------------------------------------------------ 7. the runner
@pytest.fixture(scope="module")
def n41():
    sys.path.insert(0, str(HARNESS))
    spec = importlib.util.spec_from_file_location("n41_polar_arms", HARNESS / "n41_polar_arms.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_the_runner_scores_the_floor_and_the_two_arms_on_the_nine_tasks(n41):
    assert n41.ARMS == (harness.FLOOR, "fgw_polar", "fgw_aligned_polar_series")
    assert not set(n41.ARMS) & {"fgw", "fgw_aligned_series"}  # never rescored
    tasks = n41.nt.all_tasks()
    assert set(tasks) == {"half", "late", "new1", "new2", "e75", "e725", "x1", "x2", "x3"}
    for name in n41.NEW:
        assert all(harness.CANDIDATES[name].applies(t) for t in tasks.values())


def test_the_runners_task_file_loads_through_nine_tasks(tmp_path, monkeypatch, n41):
    calls = []

    def fake_run_task(t, cands, seeds, data_dir):  # synthetic scores: no data, no scorer
        calls.append((t, [c.name for c in cands], list(seeds)))
        rows = [{"seed": s, "scores": {c.name: 50.0 + 0.1 * s + i for i, c in enumerate(cands)}} for s in seeds]
        return {"task": {"key": t.key, "inputs": dict(t.inputs)}, "rows": rows,
                "provenance": {"commit": "c" * 40, "tree_clean": True}}

    monkeypatch.setattr(n41.harness, "run_task", fake_run_task)
    out = tmp_path / "runs"
    n41.main(["run", "--task", "e75", "--out", str(out), "--data", str(tmp_path / "no_data")])
    t = n41.nt.all_tasks()["e75"]
    assert calls == [(t, list(n41.ARMS), list(n41.SEEDS))]
    (path,) = out.glob("T2e_*.json")
    assert json.loads(path.read_text())["task_source"] == "experiments/harness/n41_polar_arms.py"
    scores, why, meta = n41.nt._load([out])
    assert why == [] and meta == ["runs:ccccccc"]
    assert set(scores[t.key]) == set(n41.ARMS) and all(len(v) == 10 for v in scores[t.key].values())
    with pytest.raises(SystemExit, match="exists"):
        n41.main(["run", "--task", "e75", "--out", str(out), "--data", str(tmp_path / "no_data")])
    assert len(calls) == 1
    with pytest.raises(SystemExit, match="unknown task"):
        n41.main(["run", "--task", "nope", "--out", str(out)])


# ------------------------------------------------------------------------------------------ the arms build
@pytest.fixture(scope="module")
def emb_expr():
    """The synthetic embryo stages and contexts of tests/test_emb_expr.py (no data file)."""
    spec = importlib.util.spec_from_file_location("emb_expr_helpers", REPO / "tests" / "test_emb_expr.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.mark.parametrize("arm, base", [("fgw_polar", "fgw"), ("fgw_aligned_polar_series", "fgw_aligned_series")])
def test_the_arms_build_on_a_synthetic_embryo_and_move_only_the_positions(emb_expr, arm, base):
    s = emb_expr.stages()
    c = emb_expr.ctx(inputs=s, roles=("before", "after"))
    a, again, b = harness.CANDIDATES[arm].build(c), harness.CANDIDATES[arm].build(c), harness.CANDIDATES[base].build(c)
    assert emb_expr.digest(a) == emb_expr.digest(again)  # same seed, same output
    assert np.array_equal(a.X, b.X)  # the same draws and expression
    assert not np.array_equal(a.coords, b.coords) and np.isfinite(a.coords).all()
