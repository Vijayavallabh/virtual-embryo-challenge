"""The mirror-only correction of FGW's partners (campaign branch campaign/emb-mirror):
vec.t2_build.mirror_fixed_partners, interpolate_clouds(..., partners="mirror_fixed"), the two harness arms and the
runner experiments/harness/n42_mirrorfix_arms.py. Synthetic clouds only; no data file, committed result or scorer is
read.

interpolate_clouds computes its own FGW coupling and cannot take a fixed one, so the tests that need a fixed pairing
(the mirrored copy, the mix, the singular and non-finite cases) call the correction itself,
`mirror_fixed_partners(X, Y)`, with X the anchors and Y their partners (the identity pairing, mirrored where stated),
as interpolate_clouds calls it for each stage."""

from __future__ import annotations

import dataclasses
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pytest

from vec import harness
from vec import t2_build as tb

REPO = Path(__file__).resolve().parents[1]
HARNESS = REPO / "experiments" / "harness"
BASE = "b7d6fa0"  # campaign/emb-polar, the branch's base
ARMS = ("fgw_mirrorfix", "fgw_mirrorfix_series")


def _helpers(name: str):
    spec = importlib.util.spec_from_file_location(f"{name}_helpers", REPO / "tests" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def polar():
    """tests/test_emb_polar.py's synthetic pair, base loader and comparison (no data file)."""
    return _helpers("test_emb_polar")


def _rotation(seed: int) -> np.ndarray:
    q, r = np.linalg.qr(np.random.default_rng(seed).normal(size=(3, 3)))
    q = q * np.sign(np.diag(r))
    return q if np.linalg.det(q) > 0 else -q


def _symmetric_cloud(seed: int, n_half: int = 500) -> tuple[np.ndarray, np.ndarray]:
    """A left-right symmetric cloud: a half with x in [1, 2] and its mirror image across x = 0, thinnest along x,
    turned by a fixed rotation Q so that the symmetric direction is Q e_x. Returns the cloud and Q."""
    rng = np.random.default_rng(seed)
    H = np.c_[rng.uniform(1.0, 2.0, n_half), rng.normal(scale=3.0, size=n_half), rng.normal(scale=2.0, size=n_half)]
    S = np.vstack([H, H * [-1.0, 1.0, 1.0]])
    Q = _rotation(seed + 100)
    return S @ Q.T + [5.0, -2.0, 1.0], Q


def _mirror(C: np.ndarray, Q: np.ndarray, rows) -> np.ndarray:
    """C with the given rows mirrored across the cloud's symmetric plane (direction Q e_x through the centre)."""
    c = np.asarray([5.0, -2.0, 1.0])
    M = Q @ np.diag([-1.0, 1.0, 1.0]) @ Q.T
    out = C.copy()
    out[rows] = (C[rows] - c) @ M.T + c
    return out


def _thinnest_sd(C: np.ndarray) -> float:
    X = C - C.mean(0)
    return float(np.sqrt(np.clip(np.linalg.eigvalsh(X.T @ X / len(X)), 0.0, None)[0]))


def _side(C: np.ndarray, Q: np.ndarray) -> np.ndarray:
    return np.sign((C - [5.0, -2.0, 1.0]) @ Q[:, 0])


# ------------------------------------------------------------------------------------------ 1. partners unset
@pytest.mark.parametrize("move_coupling", tb.MOVE_COUPLINGS)
def test_partners_unset_is_the_entry_bit_for_bit(tmp_path, polar, move_coupling):
    A, B, E, EB = polar._pair()
    base = polar._base_module(tmp_path, BASE)
    old = base.interpolate_clouds(A, B, E, EB, 0.4, 300, move_coupling=move_coupling, with_partners=True,
                                  **polar.KW)
    new = tb.interpolate_clouds(A, B, E, EB, 0.4, 300, move_coupling=move_coupling, with_partners=True, **polar.KW)
    coupled = tb.interpolate_clouds(A, B, E, EB, 0.4, 300, move_coupling=move_coupling, with_partners=True,
                                    partners="coupled", **polar.KW)
    polar._same(old, new)
    polar._same(new, coupled)
    assert np.array_equal(old.partners, new.partners) and np.array_equal(new.partners, coupled.partners)
    assert "mirror" not in new.record and "partner_fix" not in new.record
    for moves in ("polar", "polar_frame"):  # the base's other moves stay as they were
        polar._same(base.interpolate_clouds(A, B, E, EB, 0.4, 300, move_coupling=move_coupling, moves=moves,
                                            **polar.KW),
                    tb.interpolate_clouds(A, B, E, EB, 0.4, 300, move_coupling=move_coupling, moves=moves,
                                          **polar.KW))
    fixed = tb.interpolate_clouds(A, B, E, EB, 0.4, 300, move_coupling=move_coupling, with_partners=True,
                                  partners="mirror_fixed", **polar.KW)
    assert np.array_equal(fixed.rows, new.rows) and np.array_equal(fixed.stage, new.stage)  # the same draws
    assert np.array_equal(fixed.partners, new.partners)  # the coupling and its partner draw are unchanged
    assert max(fixed.record["mirror"]["share_mirrored"]) > 0.0
    assert not np.array_equal(fixed.coords, new.coords)  # the idea is not idle here
    with pytest.raises(ValueError, match="partners"):
        tb.interpolate_clouds(A, B, E, EB, 0.4, 30, partners="aligned", **polar.KW)
    for moves in ("polar", "polar_frame"):
        with pytest.raises(ValueError, match="moves"):
            tb.interpolate_clouds(A, B, E, EB, 0.4, 30, partners="mirror_fixed", moves=moves, **polar.KW)


# ------------------------------------------------------------------------------------------ 2. a mirrored copy
def test_a_mirrored_copy_is_restored():
    """A left-right symmetric cloud paired with its copy mirrored across the symmetric direction (identity pairing,
    every partner mirrored). The smallest eigenvalue is -1 before and +1 after; every partner returns to its anchor's
    side (to its anchor, to 1e-9); at t = 1/2 the entry's move flattens the cloud along the symmetric direction (its
    thinnest principal SD falls to 1e-9 of the brackets') while the corrected move keeps it (relative tolerance 1e-9,
    stated here)."""
    X, Q = _symmetric_cloud(1)
    Y = _mirror(X, Q, slice(None))
    assert np.array_equal(_side(Y, Q), -_side(X, Q))
    fixed, rec = tb.mirror_fixed_partners(X, Y)
    assert rec["xcorr_min_before"] == pytest.approx(-1.0, abs=1e-9)
    assert rec["xcorr_min_after"] == pytest.approx(1.0, abs=1e-9)
    assert rec["share_mirrored"] == 1.0
    assert np.array_equal(_side(fixed, Q), _side(X, Q))
    assert np.abs(fixed - X).max() < 1e-9
    bracket = _thinnest_sd(X)
    entry_half, fixed_half = X + 0.5 * (Y - X), X + 0.5 * (fixed - X)
    assert _thinnest_sd(entry_half) < 1e-9 * bracket
    assert _thinnest_sd(fixed_half) == pytest.approx(bracket, rel=1e-9)
    assert _thinnest_sd(Y) == pytest.approx(bracket, rel=1e-9)


# ------------------------------------------------------------------------------------------ 3. right sides
def test_partners_on_their_anchors_sides_are_left_bit_for_bit():
    X, Q = _symmetric_cloud(2)
    Y = X + np.random.default_rng(2).normal(scale=0.01, size=X.shape)
    fixed, rec = tb.mirror_fixed_partners(X, Y)
    assert np.array_equal(fixed, Y) and rec["share_mirrored"] == 0.0
    assert rec["xcorr_min_before"] == rec["xcorr_min_after"] > 0.0
    fixed, rec = tb.mirror_fixed_partners(X, X)  # the identity: every a * b > 0 or a = b = 0
    assert np.array_equal(fixed, X) and rec["share_mirrored"] == 0.0


def test_with_no_partner_on_the_wrong_side_interpolate_clouds_is_the_entry(polar, monkeypatch):
    """interpolate_clouds with partners="mirror_fixed", when the correction mirrors no partner (here: the routine
    returning its input as it does then), gives the entry's output bit for bit; only the record gains the shares."""
    A, B, E, EB = polar._pair()
    seen = []

    def idle(X, Y):
        seen.append(len(X))
        return np.asarray(Y, dtype=np.float64), {"xcorr_min_before": 0.5, "xcorr_min_after": 0.5,
                                                 "share_mirrored": 0.0}

    monkeypatch.setattr(tb, "mirror_fixed_partners", idle)
    entry = tb.interpolate_clouds(A, B, E, EB, 0.4, 300, with_partners=True, **polar.KW)
    fixed = tb.interpolate_clouds(A, B, E, EB, 0.4, 300, with_partners=True, partners="mirror_fixed", **polar.KW)
    assert seen == [200, 200]  # each stage's anchors, once
    for f in ("stage", "rows", "coords", "pool_coords", "pool_weights", "partners"):
        assert np.array_equal(getattr(entry, f), getattr(fixed, f)), f
    rec = dict(fixed.record)
    assert rec.pop("partner_fix") == "mirror_fixed" and rec.pop("mirror")["share_mirrored"] == [0.0, 0.0]
    assert rec == entry.record


# ------------------------------------------------------------------------------------------ 4. a mix
def test_only_the_mirrored_partners_change():
    X, Q = _symmetric_cloud(3)
    wrong = np.sort(np.random.default_rng(3).choice(len(X), 300, replace=False))
    Y = _mirror(X, Q, wrong)
    fixed, rec = tb.mirror_fixed_partners(X, Y)
    changed = np.flatnonzero((fixed != Y).any(1))
    assert np.array_equal(changed, wrong) and rec["share_mirrored"] == pytest.approx(len(wrong) / len(X))
    keep = np.setdiff1d(np.arange(len(X)), wrong)
    assert np.array_equal(fixed[keep], Y[keep])  # the others are the entry's, bit for bit
    assert np.array_equal(_side(fixed, Q), _side(X, Q))  # every mirrored partner is back on its anchor's side
    assert rec["xcorr_min_after"] > rec["xcorr_min_before"]


# ------------------------------------------------------------------------------------------ 5. singular, non-finite
def test_a_singular_covariance_or_a_non_finite_value_keeps_the_entrys_partners():
    rng = np.random.default_rng(5)
    flat = np.c_[rng.normal(size=(200, 2)), np.zeros(200)]  # coplanar anchors
    Y = flat * [1.0, 1.0, 1.0] + rng.normal(size=(200, 3))
    fixed, rec = tb.mirror_fixed_partners(flat, Y * [-1.0, 1.0, 1.0])
    assert np.array_equal(fixed, Y * [-1.0, 1.0, 1.0])
    assert rec == {"xcorr_min_before": None, "xcorr_min_after": None, "share_mirrored": 0.0}
    X, Q = _symmetric_cloud(5)
    Ym = _mirror(X, Q, slice(None))
    fixed, rec = tb.mirror_fixed_partners(X, np.c_[Ym[:, :2], np.zeros(len(X))])  # coplanar partners
    assert np.array_equal(fixed, np.c_[Ym[:, :2], np.zeros(len(X))]) and rec["share_mirrored"] == 0.0
    for bad in (np.nan, np.inf):
        Yb = Ym.copy()
        Yb[7, 1] = bad
        fixed, rec = tb.mirror_fixed_partners(X, Yb)
        assert np.array_equal(fixed, Yb, equal_nan=True) and rec["xcorr_min_before"] is None
        Xb = X.copy()
        Xb[3, 0] = bad
        fixed, rec = tb.mirror_fixed_partners(Xb, Ym)
        assert np.array_equal(fixed, Ym) and rec["share_mirrored"] == 0.0
    fixed, rec = tb.mirror_fixed_partners(X[:3], Ym[:3])  # fewer pairs than the dimensions need
    assert np.array_equal(fixed, Ym[:3]) and rec["xcorr_min_before"] is None


# ------------------------------------------------------------------------------------------ 6. reproducible
def test_the_corrected_routine_is_reproducible_and_records_its_shares(polar):
    A, B, E, EB = polar._pair()
    a = tb.interpolate_clouds(A, B, E, EB, 0.5, 300, partners="mirror_fixed", **polar.KW)
    b = tb.interpolate_clouds(A, B, E, EB, 0.5, 300, partners="mirror_fixed", **polar.KW)
    polar._same(a, b)
    rec = a.record["mirror"]
    assert a.record["partner_fix"] == "mirror_fixed"
    assert set(rec) == {"xcorr_min_before", "xcorr_min_after", "share_mirrored"}
    assert all(len(v) == 2 for v in rec.values()) and all(0.0 <= s <= 1.0 for s in rec["share_mirrored"])
    assert tb.rms_radius(a.coords) == pytest.approx(10.0)
    other = tb.interpolate_clouds(A, B, E, EB, 0.5, 300, partners="mirror_fixed", **{**polar.KW, "seed": 2})
    assert not np.array_equal(a.coords, other.coords)
    aligned = tb.interpolate_clouds(A, B, E, EB, 0.5, 300, partners="mirror_fixed", move_coupling="aligned",
                                    **polar.KW)  # composes with move_coupling
    assert aligned.record["move_coupling"] == "aligned" and "mirror" in aligned.record


# ------------------------------------------------------------------------------------------ the harness arms
def test_the_two_arms_are_registered_for_t2_and_call_fgw_with_the_correction(monkeypatch):
    seen = {}
    monkeypatch.setattr(harness, "_fgw", lambda ctx, **kw: seen.setdefault(ctx, kw))
    for name, kw in (("fgw_mirrorfix", {"partners": "mirror_fixed"}),
                     ("fgw_mirrorfix_series", {"partners": "mirror_fixed", "expression": "series"})):
        c = harness.CANDIDATES[name]
        assert c.tasks == frozenset({"T2"}) and c.name == name and c.params["partners"] == "mirror_fixed"
        assert c.build(name) == kw
    assert harness.CANDIDATES["fgw_mirrorfix"].requires is harness.CANDIDATES["fgw"].requires
    assert harness.CANDIDATES["fgw_mirrorfix_series"].requires is harness.CANDIDATES["fgw_series"].requires


@pytest.fixture(scope="module")
def emb_expr():
    """The synthetic embryo stages and contexts of tests/test_emb_expr.py (no data file)."""
    return _helpers("test_emb_expr")


@pytest.mark.parametrize("arm, base", [("fgw_mirrorfix", "fgw"), ("fgw_mirrorfix_series", "fgw_series")])
def test_the_arms_build_on_a_synthetic_embryo_and_move_only_the_positions(emb_expr, arm, base):
    s = emb_expr.stages()
    c = emb_expr.ctx(inputs=s, roles=("before", "after"))
    a, again, b = harness.CANDIDATES[arm].build(c), harness.CANDIDATES[arm].build(c), harness.CANDIDATES[base].build(c)
    assert emb_expr.digest(a) == emb_expr.digest(again)  # same seed, same output
    assert np.array_equal(a.X, b.X)  # the same draws and expression
    assert not np.array_equal(a.coords, b.coords) and np.isfinite(a.coords).all()


# ------------------------------------------------------------------------------------------ 7. the runner
@pytest.fixture(scope="module")
def n42():
    sys.path.insert(0, str(HARNESS))
    spec = importlib.util.spec_from_file_location("n42_mirrorfix_arms", HARNESS / "n42_mirrorfix_arms.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_the_runner_scores_the_floor_and_the_two_arms_on_the_nine_tasks(n42):
    assert n42.NEW == ARMS and n42.ARMS == (harness.FLOOR, *ARMS)
    assert not set(n42.ARMS) & {"fgw", "fgw_aligned_series"}  # never rescored
    tasks = n42.nt.all_tasks()
    assert set(tasks) == {"half", "late", "new1", "new2", "e75", "e725", "x1", "x2", "x3"}
    for name in n42.NEW:
        assert all(harness.CANDIDATES[name].applies(t) for t in tasks.values())


def test_the_runners_task_file_loads_through_nine_tasks(tmp_path, monkeypatch, n42):
    calls = []

    def fake_run_task(t, cands, seeds, data_dir):  # synthetic scores: no data, no scorer
        calls.append((t, [c.name for c in cands], list(seeds)))
        rows = [{"seed": s, "scores": {c.name: 50.0 + 0.1 * s + i for i, c in enumerate(cands)}} for s in seeds]
        return {"task": {"key": t.key, "inputs": dict(t.inputs)}, "rows": rows,
                "provenance": {"commit": "c" * 40, "tree_clean": True}}

    monkeypatch.setattr(n42.harness, "run_task", fake_run_task)
    out = tmp_path / "runs"
    n42.main(["run", "--task", "e75", "--out", str(out), "--data", str(tmp_path / "no_data")])
    t = n42.nt.all_tasks()["e75"]
    assert calls == [(t, list(n42.ARMS), list(n42.SEEDS))]  # one call with the floor and the two arms
    (path,) = out.glob("T2e_*.json")
    data = json.loads(path.read_text())
    assert data["task_source"] == "experiments/harness/n42_mirrorfix_arms.py" and data["not_buildable"] == {}
    scores, why, meta = n42.nt._load([out])
    assert why == [] and meta == ["runs:ccccccc"]
    assert set(scores[t.key]) == set(n42.ARMS) and all(len(v) == 10 for v in scores[t.key].values())
    with pytest.raises(SystemExit, match="exists"):
        n42.main(["run", "--task", "e75", "--out", str(out), "--data", str(tmp_path / "no_data")])
    assert len(calls) == 1
    with pytest.raises(SystemExit, match="unknown task"):
        n42.main(["run", "--task", "nope", "--out", str(out)])


def _fake_runner(score: dict[str, float]):
    """run_task's shape without data or scorer: each candidate is built (its build may raise) and given a constant."""
    calls = []

    def runner(t, cands, seeds, data_dir):
        calls.append([c.name for c in cands])
        rows = []
        for s in seeds:
            row = {"seed": s, "scores": {"floor": 50.0}, "skills": {"floor": {}}, "metrics": {"floor": {}},
                   "n_cells": {"floor": 1}, "groups": {"floor": {}}, "reads": {}}
            for c in cands:
                if c.name == "floor":
                    continue
                c.build(dataclasses.make_dataclass("Ctx", ["seed"])(seed=harness.BUILD_SEED_OFFSET + s))
                row["scores"][c.name] = score[c.name]
                for k in ("skills", "metrics", "groups"):
                    row[k][c.name] = {}
                row["n_cells"][c.name], row["reads"][c.name] = 1, []
            rows.append(row)
        return {"task": {"key": t.key, "inputs": dict(t.inputs)}, "rows": rows,
                "candidates": [c.name for c in cands], "params": {c.name: {} for c in cands},
                "provenance": {"commit": "a" * 40, "tree_clean": True}}

    return runner, calls


def _patch_builds(monkeypatch, n42, fail_on: set[str], key: str, seed_fails: int = 0):
    """Every new arm's build replaced: the arms in fail_on raise the builder's cloud-guard refusal from seed_fails."""
    cands = harness.CANDIDATES

    def ok(ctx):
        return "built"

    def bad(ctx):
        if ctx.seed - harness.BUILD_SEED_OFFSET >= seed_fails:
            raise ValueError(f"{key}: cloud guard flags ['rms: RMS radius is 0.000 x the source, expected 0.436']")
        return "built"

    for x in n42.NEW:
        monkeypatch.setitem(cands, x, dataclasses.replace(cands[x], build=bad if x in fail_on else ok))


def test_an_arm_the_builder_refuses_is_recorded_and_the_other_arm_is_still_scored(tmp_path, monkeypatch, n42):
    t = n42.nt.all_tasks()["e725"]
    score = {"fgw_mirrorfix": 53.0, "fgw_mirrorfix_series": 0.0}
    _patch_builds(monkeypatch, n42, {"fgw_mirrorfix_series"}, t.key, seed_fails=2)
    runner, calls = _fake_runner(score)
    res = n42.score_task(t, tmp_path, runner=runner)
    assert calls == [["floor", *n42.NEW], *[["floor", n] for n in n42.NEW]]  # together first, then arm by arm
    assert set(res["not_buildable"]) == {"fgw_mirrorfix_series"}
    rec = res["not_buildable"]["fgw_mirrorfix_series"]
    assert rec["status"] == n42.n36.NOT_BUILDABLE and rec["seed"] == 2
    assert rec["message"] == ("ValueError: <task>: cloud guard flags ['rms: RMS radius is # x the source, "
                              "expected #']")  # the task key named as such, every number masked
    scored = {"floor", "fgw_mirrorfix"}
    assert all(set(r["scores"]) == scored for r in res["rows"]) and set(res["candidates"]) == scored
    assert [r["seed"] for r in res["rows"]] == list(n42.SEEDS)
    monkeypatch.setattr(n42.harness, "run_task", runner)
    out = tmp_path / "runs"
    n42.main(["run", "--task", "e725", "--out", str(out), "--data", str(tmp_path / "no_data")])
    scores, why, _ = n42.nt._load([out])
    assert why == [] and set(scores[t.key]) == scored
    (path,) = out.glob("T2e_*.json")
    assert set(json.loads(path.read_text())["not_buildable"]) == {"fgw_mirrorfix_series"}
    # both arms failing still leaves the floor scored
    _patch_builds(monkeypatch, n42, set(n42.NEW), t.key)
    runner, calls = _fake_runner(score)
    res = n42.score_task(t, tmp_path, runner=runner)
    assert set(res["not_buildable"]) == set(n42.NEW) and calls[-1] == ["floor"]
    assert all(set(r["scores"]) == {"floor"} for r in res["rows"])
    # nothing fails: one run, nothing recorded
    _patch_builds(monkeypatch, n42, set(), t.key)
    runner, calls = _fake_runner(score)
    res = n42.score_task(t, tmp_path, runner=runner)
    assert calls == [["floor", *n42.NEW]] and res["not_buildable"] == {}


def test_a_failure_outside_a_build_stops_the_task(tmp_path, n42):
    def runner(t, cands, seeds, data_dir):
        raise OSError("the data are not there")

    with pytest.raises(OSError, match="not there"):
        n42.score_task(n42.nt.all_tasks()["x2"], tmp_path, runner=runner)
