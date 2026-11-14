"""The embryo batch I (campaign branch campaign/emb-i): the reference-residual render (vec.t2_build.RESID_RENDER_RULE,
the arm fgw_resid) and the runner experiments/harness/n49_embryo_i_arms.py. Synthetic stages (tests/test_emb_expr.py's
builder, with a fourth stage added here) and synthetic count matrices only; no data file, committed result or scorer
is read, and no prediction is compared with any target."""

from __future__ import annotations

import dataclasses
import importlib.util
import json
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np
import pytest
from scipy import sparse

from vec import harness
from vec import t2_build as tb

REPO = Path(__file__).resolve().parents[1]
HARNESS = REPO / "experiments" / "harness"
ARM = "fgw_resid"


@pytest.fixture(scope="module")
def ex():
    """tests/test_emb_expr.py's synthetic stages, contexts and pins (a private copy), with a fourth stage "later" at
    E8.5 so that both brackets of the E7.25 + E8.0 task have an earlier and a later sample."""
    spec = importlib.util.spec_from_file_location("test_emb_expr_helpers_i", REPO / "tests" / "test_emb_expr.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    mod.STAGES["later"] = (200, 200.0, 3.0, 6.0)
    mod.DAYS["later"] = 8.5
    mod.FILES["later"] = "T2/E8.5.h5ad"
    return mod


def four_ctx(ex, inputs=None, seed=harness.BUILD_SEED_OFFSET + 4):
    """E7.25 + E8.0 -> E7.5 with trend samples E6.75 and E8.5: four samples, both brackets interior."""
    return ex.trend_ctx(inputs=ex.stages() if inputs is None else inputs, extra=("earlier", "later"), seed=seed)


def _resid(ex, c):
    res, *_ = ex._interp(c)
    t = c.task
    return res, tb.resid_rendered_expression(c.inputs["before"], c.inputs["after"], res, list(ex.BOARD.genes),
                                             [(a, t.role_days[r]) for r, a in c.inputs.items()], [c.seed, 7],
                                             sample_names=[t.role_stems[r] for r in c.inputs])


# ------------------------------------------------------------------------------------- existing candidates unchanged
@pytest.mark.parametrize("roles, name", [(("earlier", "before", "after"), "fgw"),
                                         (("earlier", "before", "after"), "fgw_trend_render"),
                                         (("before", "after"), "fgw")])
def test_fgw_and_the_trend_render_are_bit_for_bit_what_they_were(ex, roles, name):
    c = ex.ctx(inputs=ex.stages(), roles=roles)
    assert ex.digest(harness.CANDIDATES[name].build(c)) == ex.PINS_AT_E709655[roles][name]


def test_fgw_on_the_four_sample_task_is_fgw_on_the_two_bracket_task(ex):
    s = ex.stages()
    c = four_ctx(ex, s)
    assert ex.digest(harness.CANDIDATES["fgw"].build(c)) == ex.PINS_AT_E709655[("before", "after")]["fgw"]


# ------------------------------------------------------------------------------------------------ neighbours and L
def test_neighbours_w1_w2_and_l_by_hand():
    days = [7.25, 7.5, 8.0, 6.75]
    assert tb.resid_neighbours(days, 1) == (0, 2)
    assert tb.resid_neighbours(days, 0) == (3, 1)
    assert tb.resid_neighbours(days, 3) == (None, 0)  # the earliest sample has no earlier one
    assert tb.resid_neighbours(days, 2) == (1, None)  # the latest has no later one
    w1, w2, L = tb.resid_weights(7.5, 7.25, 8.0)
    assert Fraction(w1).limit_denominator(100) == Fraction(2, 3) and w1 + w2 == 1.0
    assert L == pytest.approx(9 / 14, rel=1e-15)
    w1, w2, L = tb.resid_weights(7.25, 6.75, 8.0)
    assert (w1, w2) == pytest.approx((0.6, 0.4)) and L == pytest.approx(1 / 1.52)
    with pytest.raises(ValueError, match="strictly between"):
        tb.resid_weights(8.0, 7.25, 8.0)
    with pytest.raises(ValueError, match="nearest day"):
        tb.resid_neighbours([7.0, 7.5, 7.0, 8.0], 1)


def _prof(M):
    M = M.toarray().astype(float)
    return (M / np.maximum(M.sum(1, keepdims=True), 1)).mean(0)


def test_the_factor_before_rescale_and_clip_is_exp_of_minus_l_r_from_profiles_by_hand():
    rng = np.random.default_rng(7)
    C = [sparse.csr_matrix(rng.poisson(lam, size=(150, 8)) + 1) for lam in (3.0, 5.0, 9.0)]
    names = [np.array(["EXE-A"] * 150)] * 3
    w1, w2, L = tb.resid_weights(7.5, 7.25, 8.0)
    raw, phi, den = tb.resid_factor(C[1], names[1], C[0], names[0], C[2], names[2], w1, w2, L, ["EXE-A"])
    eps = tb.SERIES_EPS
    p1, own, p2 = (_prof(M) for M in C)
    r = np.log(own + eps) - (w1 * np.log(p1 + eps) + w2 * np.log(p2 + eps))
    np.testing.assert_allclose(raw, np.exp(-L * r), rtol=1e-12)
    np.testing.assert_allclose(den, own, rtol=1e-12)
    want = np.clip(raw / (own @ raw / own.sum()), *tb.SERIES_CLIP)
    np.testing.assert_allclose(phi, want, rtol=1e-12)
    # equal profiles: the residual is 0 and the factor exactly 1, before and after the rescale
    raw1, phi1, _ = tb.resid_factor(C[0], names[0], C[0], names[0], C[0], names[0], w1, w2, L, ["EXE-A"])
    assert np.all(raw1 == 1.0) and np.all(phi1 == 1.0)
    with pytest.raises(ValueError, match="no extraembryonic label"):
        tb.resid_factor(C[1], names[1], C[0], names[0], C[2], names[2], w1, w2, L, [])


# ------------------------------------------------------------------------------------------------ the arm
def test_registration(ex):
    c, tr = harness.CANDIDATES[ARM], harness.CANDIDATES["fgw_trend_render"]
    assert c.name == ARM and c.tasks == frozenset({"T2"}) and not c.tuned_on and not c.externals and not c.derived
    assert c.requires is tr.requires and c.params["rule"] == tb.RESID_RENDER_RULE and c.params["scope"] == "whole"
    assert c.applies(ex.trend_task()) and not c.applies(ex.task(("before", "after")))


def test_both_interior_brackets_are_corrected_at_their_own_day_on_fgws_cells_and_positions(ex):
    s = ex.stages()
    c = four_ctx(ex, s)
    rec: dict = {}
    x = harness._fgw_resid(c, build_record=rec)
    f = harness.CANDIDATES["fgw"].build(c)
    np.testing.assert_array_equal(x.coords, f.coords)  # the same draws and positions
    assert x.n == f.n == ex.BOARD.max_cells and not np.array_equal(x.X, f.X)
    b = rec["resid"]["brackets"]
    assert b["before"]["corrected"] and b["after"]["corrected"]
    assert b["before"]["neighbours"] == ["E6.75", "E8.0"] and b["after"]["neighbours"] == ["E7.25", "E8.5"]
    assert (b["before"]["w1"], b["before"]["w2"]) == pytest.approx((0.6, 0.4))
    assert (b["after"]["w1"], b["after"]["w2"]) == pytest.approx((0.4, 0.6))
    assert b["before"]["L"] == pytest.approx(1 / 1.52) and b["after"]["L"] == pytest.approx(1 / 1.52)
    for v in b.values():
        assert 0.0 <= v["clipped_share"] <= 1.0 and v["mean_abs_log_factor"] > 0.0
    res, (X, rec2) = _resid(ex, c)
    assert rec2 == rec["resid"]
    np.testing.assert_array_equal(X, x.X)
    assert np.isfinite(x.X).all() and (x.X >= 0).all() and ((x.X > 0).sum(1) >= 1).all()


def test_a_bracket_at_the_edge_of_the_series_is_left_as_drawn(ex):
    """E6.75 + (E7.25, E8.0): the before bracket has both neighbours (E6.75, E8.0), the after bracket no later one."""
    s = ex.stages()
    c = ex.trend_ctx(inputs=s)
    rec: dict = {}
    x = harness._fgw_resid(c, build_record=rec)
    f = harness.CANDIDATES["fgw"].build(c)
    res, *_ = ex._interp(c)
    assert rec["resid"]["brackets"]["after"] == {"corrected": False}
    assert rec["resid"]["brackets"]["before"]["corrected"]
    np.testing.assert_array_equal(x.coords, f.coords)
    np.testing.assert_array_equal(x.X[res.stage == 1], f.X[res.stage == 1])
    assert not np.array_equal(x.X[res.stage == 0], f.X[res.stage == 0])


def test_with_only_the_two_brackets_as_samples_the_arm_is_fgw_bit_for_bit(ex):
    c = ex.ctx(inputs=ex.stages(), roles=("before", "after"))
    rec: dict = {}
    x = harness._fgw_resid(c, build_record=rec)
    assert ex.digest(x) == ex.PINS_AT_E709655[("before", "after")]["fgw"]
    assert rec["resid"]["brackets"] == {"before": {"corrected": False}, "after": {"corrected": False}}


def test_a_task_where_no_bracket_has_both_neighbours_reads_no_label(ex):
    """The idle arm reads no label: samples without a shared extraembryonic label do not stop it."""
    s = ex.stages()
    for a in s.values():
        a.obs["celltype"] = a.obs["celltype"].astype(str).str.replace("EXE-", "EMB-")
    c = ex.ctx(inputs=s, roles=("before", "after"))
    assert ex.digest(harness.CANDIDATES[ARM].build(c)) == ex.digest(harness.CANDIDATES["fgw"].build(c))
    with pytest.raises(ValueError, match="no extraembryonic label held by every sample"):
        harness.CANDIDATES[ARM].build(four_ctx(ex, s))


def _flat_exe(a, row):
    """a with every EXE-Endoderm cell's panel and outside counts set to `row` (X recomputed as released)."""
    a = a.copy()
    C = a.layers["counts"].toarray()
    m = np.asarray(a.obs["celltype"]) == "EXE-Endoderm"
    C[m] = row
    a.layers["counts"] = sparse.csr_matrix(C.astype(np.int32))
    a.X = sparse.csr_matrix(np.log1p(1e4 * C / C.sum(1, keepdims=True)).astype(np.float32))
    return a


def test_a_bracket_on_its_neighbours_line_is_unchanged_its_factor_1_and_no_draw(ex):
    """The three samples hold the same cells (E6.75 and E8.0 copies of E7.25): the before bracket's profile lies on the
    line through its neighbours with a residual of exactly 0, so its factor is 1 on every gene and nothing is drawn;
    the after bracket has no later sample. The file is fgw's bit for bit."""
    s = ex.stages()
    s = {**s, "earlier": s["before"].copy(), "after": s["before"].copy()}
    c = ex.trend_ctx(inputs=s)
    rec: dict = {}
    x = harness._fgw_resid(c, build_record=rec)
    b = rec["resid"]["brackets"]
    assert b["before"]["corrected"] and b["before"]["mean_abs_log_factor"] == 0.0
    assert b["before"]["clipped_share"] == 0.0 and b["after"] == {"corrected": False}
    assert ex.digest(x) == ex.digest(harness.CANDIDATES["fgw"].build(c))


def test_a_profile_on_the_line_up_to_rounding_gets_a_factor_of_1_up_to_rounding(ex):
    """Every sample's EXE-Endoderm cells hold one row, the stages different cell numbers: the three profiles agree up
    to the rounding of their means, so the factor is 1 to within that rounding (the rule has no tolerance: such a
    factor is applied, and the cells come back as released up to the draw)."""
    s = ex.stages()
    row = np.asarray(s["before"].layers["counts"][0].toarray()).ravel() + 1
    s = {k: _flat_exe(a, row) for k, a in s.items()}
    c = ex.trend_ctx(inputs=s)
    rec: dict = {}
    harness._fgw_resid(c, build_record=rec)
    b = rec["resid"]["brackets"]["before"]
    assert b["corrected"] and b["mean_abs_log_factor"] < 1e-12 and b["clipped_share"] == 0.0


def test_same_seed_same_output_and_another_seed_another(ex):
    s = ex.stages()
    a = harness.CANDIDATES[ARM].build(four_ctx(ex, s))
    b = harness.CANDIDATES[ARM].build(four_ctx(ex, s))
    d = harness.CANDIDATES[ARM].build(four_ctx(ex, s, seed=harness.BUILD_SEED_OFFSET + 5))
    assert ex.digest(a) == ex.digest(b) != ex.digest(d)


def test_the_render_refuses_brackets_outside_its_samples(ex):
    c = ex.trend_ctx()
    res, *_ = ex._interp(c)
    with pytest.raises(ValueError, match="among the residual render's samples"):
        tb.resid_rendered_expression(c.inputs["before"], c.inputs["after"], res, list(ex.BOARD.genes),
                                     [(c.inputs["before"], 7.25)], 0)


# ================================================================================ the runner
@pytest.fixture(scope="module")
def n49():
    sys.path.insert(0, str(HARNESS))
    spec = importlib.util.spec_from_file_location("n49_embryo_i_arms", HARNESS / "n49_embryo_i_arms.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# per practice task, from its definition only (n36's sample rule and vec.harness.STAGE_DAYS): the corrected bracket,
# its neighbours and L; the other bracket sits at an end of its samples
CORRECTED = {"half": ("before", ["E7.25", "E8.5"]), "late": ("before", ["E7.25", "E8.25_late"]),
             "new1": ("before", ["E7.25", "E8.5"]), "new2": ("before", ["E7.25", "E8.75"]),
             "x1": ("before", ["E7.5", "E8.75"]), "x2": ("before", ["E7.5", "E8.75"]),
             "x3": ("before", ["E8.0", "E8.75"]), "e75": ("before", ["E6.75", "E8.0"]),
             "e725": ("after", ["E6.75", "E8.0"])}


def test_the_runner_scores_the_floor_fgw_the_trend_render_and_the_arm_on_the_nine_extended_tasks(n49):
    assert n49.NEW == (ARM,) and n49.ARMS == (harness.FLOOR, "fgw", "fgw_trend_render", ARM)
    tasks = n49.n36.all_tasks()
    assert set(tasks) == {"half", "late", "new1", "new2", "e75", "e725", "x1", "x2", "x3"}
    for n, t in tasks.items():
        assert all(harness.CANDIDATES[a].applies(t) for a in n49.ARMS[1:]), n
        tv = harness.TaskView.of(t)
        days = [tv.role_days[r] for r in tv.roles]
        stems = [tv.role_stems[r] for r in tv.roles]
        which, nbrs = CORRECTED[n]
        for side in ("before", "after"):
            e, lt = tb.resid_neighbours(days, tv.roles.index(side))
            assert (e is not None and lt is not None) == (side == which), (n, side)
            if side == which:
                assert [stems[e], stems[lt]] == nbrs, n


def test_the_test_pair_corrects_its_earlier_bracket_with_l_9_14(n49):
    days = {"T2/E7.75.h5ad": 7.75}
    files = n49.n36.trend_sample_files("T2/E7.5.h5ad", "T2/E8.0.h5ad", "T2/E7.75.h5ad", days=days)
    d = [n49.n36.stage_day(f, days) for f in files]
    i, j = files.index("T2/E7.5.h5ad"), files.index("T2/E8.0.h5ad")
    e, lt = tb.resid_neighbours(d, i)
    assert (files[e], files[lt]) == ("T2/E7.25.h5ad", "T2/E8.0.h5ad")
    assert tb.resid_weights(d[i], d[e], d[lt])[2] == pytest.approx(9 / 14)
    assert tb.resid_neighbours(d, j)[1] is None


def _patch_builds(monkeypatch, n49, fail_on: set[str], key: str, seed_fails: int = 0):
    cands = harness.CANDIDATES

    def ok(ctx, build_record=None):
        if build_record is not None:
            build_record["seen"] = 1
        return "built"

    def bad(ctx, build_record=None):
        if ctx.seed - harness.BUILD_SEED_OFFSET >= seed_fails:
            raise ValueError(f"{key}: no extraembryonic label held by every sample: the residual factor is not "
                             "defined 3")
        if build_record is not None:
            build_record["seen"] = 1
        return "built"

    for x in ("fgw_trend_render", ARM):
        monkeypatch.setitem(cands, x, dataclasses.replace(cands[x], build=bad if x in fail_on else ok))
    monkeypatch.setitem(cands, "fgw", dataclasses.replace(cands["fgw"], build=lambda ctx: "built"))


def _fake_runner(score: dict[str, float]):
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


def test_the_runners_task_file_loads_through_nine_tasks_with_the_rescored_arms_under_their_names(tmp_path,
                                                                                                monkeypatch, n49):
    score = {"fgw": 51.0, "fgw_trend_render": 52.0, ARM: 53.0}
    t = n49.n36.all_tasks()["half"]
    _patch_builds(monkeypatch, n49, set(), t.key)
    runner, calls = _fake_runner(score)
    monkeypatch.setattr(n49.harness, "run_task", runner)
    out = tmp_path / "runs"
    n49.main(["run", "--task", "half", "--out", str(out), "--data", str(tmp_path / "no_data")])
    assert calls == [list(n49.ARMS)]
    (path,) = out.glob("T2e_*.json")
    data = json.loads(path.read_text())
    assert data["task_source"] == "experiments/harness/n49_embryo_i_arms.py" and data["not_buildable"] == {}
    assert data["build_records"] == {ARM: {str(s): {"seen": 1} for s in n49.SEEDS}}
    assert data["task"]["inputs"] == dict(t.inputs)  # the extended task: the sample rule's inputs
    assert data["sample_rule"] == {"series": list(n49.n36.SERIES),
                                   "samples": n49.n36.trend_sample_files(t.inputs["before"], t.inputs["after"],
                                                                         t.target)}
    scores, why, meta = n49.nt._load([out])
    assert why == [] and meta == ["runs:aaaaaaa"]
    assert set(scores[t.key]) == set(n49.ARMS)
    assert all(len(v) == 10 for v in scores[t.key].values())
    with pytest.raises(SystemExit, match="exists"):
        n49.main(["run", "--task", "half", "--out", str(out), "--data", str(tmp_path / "no_data")])
    assert len(calls) == 1
    with pytest.raises(SystemExit, match="unknown task"):
        n49.main(["run", "--task", "nope", "--out", str(out)])


def test_an_arm_the_builder_refuses_is_recorded_and_the_others_are_still_scored(tmp_path, monkeypatch, n49):
    t = n49.n36.all_tasks()["e725"]
    score = {"fgw": 51.0, "fgw_trend_render": 52.0, ARM: 0.0}
    _patch_builds(monkeypatch, n49, {ARM}, t.key, seed_fails=2)
    runner, calls = _fake_runner(score)
    res = n49.score_task(t, tmp_path, runner=runner)
    assert calls == [list(n49.ARMS), ["floor", "fgw"], ["floor", "fgw_trend_render"], ["floor", ARM]]
    rec = res["not_buildable"][ARM]
    assert set(res["not_buildable"]) == {ARM} and rec["status"] == n49.n36.NOT_BUILDABLE and rec["seed"] == 2
    assert rec["message"] == ("ValueError: <task>: no extraembryonic label held by every sample: the residual factor "
                              "is not defined #")
    scored = {"floor", "fgw", "fgw_trend_render"}
    assert all(set(r["scores"]) == scored for r in res["rows"]) and set(res["candidates"]) == scored
    assert res["build_records"] == {}
    _patch_builds(monkeypatch, n49, {"fgw_trend_render", ARM}, t.key)
    runner, calls = _fake_runner(score)
    res = n49.score_task(t, tmp_path, runner=runner)
    assert set(res["not_buildable"]) == {"fgw_trend_render", ARM}
    assert all(set(r["scores"]) == {"floor", "fgw"} for r in res["rows"])


def test_a_failure_outside_a_guarded_build_stops_the_task(tmp_path, monkeypatch, n49):
    def runner(t, cands, seeds, data_dir):
        raise OSError("the data are not there")

    with pytest.raises(OSError, match="not there"):
        n49.score_task(n49.n36.all_tasks()["x2"], tmp_path, runner=runner)
    t = n49.n36.all_tasks()["x1"]
    _patch_builds(monkeypatch, n49, set(), t.key)

    def bad_fgw(ctx):
        raise ValueError("fgw fails")

    monkeypatch.setitem(harness.CANDIDATES, "fgw", dataclasses.replace(harness.CANDIDATES["fgw"], build=bad_fgw))
    runner, _ = _fake_runner({"fgw": 1.0, "fgw_trend_render": 1.0, ARM: 1.0})
    with pytest.raises(ValueError, match="fgw fails"):
        n49.score_task(t, tmp_path, runner=runner)
