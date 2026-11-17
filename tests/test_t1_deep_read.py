"""experiments/t1_m1/deep_read.py (try T54 of experiments/harness/TRIES.md) on synthetic data: the thinned reference,
where the route reads it, the arms and the reading. No released file, statistics file, board file or score is read."""

from __future__ import annotations

import importlib.util
import inspect
import json
import re
import sys
from types import SimpleNamespace

import anndata as ad
import numpy as np
import pytest
from scipy import sparse
from test_t1_lib_step_read import GOLDEN_ARMS_FD18
from test_t1_step_read import _arms_world

from vec import core, harness, t1_scope
from vec import t1_state_transport as ts
from vec.boards import REPO_ROOT
from vec.evaluate import Stage

SCRIPT = REPO_ROOT / "experiments" / "t1_m1" / "deep_read.py"
MULT = {"m1_r_gmom_x15": 1.5, "m1_r_gmom_x2": 2.0, "m1_r_gmom_x3": 3.0, "m1_r_gmom_x4": 4.0, "m1_r_gmom_x5": 5.0,
        "m1_r_gmom_x6": 6.0}
SCALED = tuple(MULT)
ALL_ARMS = ("b0", "m1_r_gmom", *SCALED)
THIN = {"ref67": (2 / 3, [54, 1]), "ref50": (1 / 2, [54, 2]), "ref67b": (2 / 3, [54, 3])}
WORLDS = ("ref100", "ref67", "ref67b", "ref50")
STAMP = {"commit": "f" * 40, "tree_clean": True, "utc": "x", "inputs": {}}
X2, X3, X4, X5, X6 = "m1_r_gmom_x2", "m1_r_gmom_x3", "m1_r_gmom_x4", "m1_r_gmom_x5", "m1_r_gmom_x6"


@pytest.fixture
def dr(monkeypatch, tmp_path):
    spec = importlib.util.spec_from_file_location("deep_read_module", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    monkeypatch.setattr(mod, "RESULTS", tmp_path / "results")
    monkeypatch.setattr(mod.gv, "stamp", lambda started, files: dict(STAMP))
    return mod


# ---------------------------------------------------------------------------------------------------------------------
# The thinned reference


def _world():
    """The synthetic route's inputs (an E9.5-like pool with its counts layer and an E8.5-like start), the scorer's
    reference stage on the board's panel, and the full reference file as read (its own .X, no counts layer)."""
    ctx, groups, dirs, held = _arms_world()
    ref = ctx.inputs["ref"]
    full = ad.AnnData(X=ref.X.copy(), obs=ref.obs[["celltype"]].copy(), var=ref.var.copy())
    genes = list(ctx.board.genes)
    ref_stage = Stage.from_anndata(full[:, genes].copy(), genes)
    return ctx, groups, dirs, held, full, ref_stage, genes


def test_the_thinning_is_t50s_binomial_on_the_worlds_streams_and_no_seed_enters(dr):
    rng = np.random.default_rng(9)
    C = sparse.csr_matrix(rng.poisson(4.0, (40, 30)))
    Cs = sparse.csr_matrix(C, dtype=np.int64)
    Cs.sort_indices()
    for world, (p, stream) in THIN.items():
        got = dr.thin_counts(C, world)
        want = np.random.default_rng(stream).binomial(Cs.data, p)
        assert (got != sparse.csr_matrix((want, Cs.indices, Cs.indptr), shape=C.shape)).nnz == 0
        assert (dr.thin_counts(C, world) != got).nnz == 0  # one thinned reference per world
    assert (dr.thin_counts(C, "ref67") != dr.thin_counts(C, "ref67b")).nnz > 0
    assert dr.THINNED == {k: {"p": p, "stream": s} for k, (p, s) in THIN.items()}
    assert dr.DRAW_NULL == ("ref67b", "ref67")
    assert set(dr.ob.THINNED) == {"lib3_t33", "lib3_t67", "lib3_t33b"}  # T50's table is left as it was
    streams = [tuple(v["stream"]) for v in dr.THINNED.values()] + [tuple(v["stream"]) for v in dr.ob.THINNED.values()]
    assert len(set(streams)) == len(streams)
    for f in (dr.thin_counts, dr.thinned_world, dr.world_reference):
        assert "seed" not in inspect.signature(f).parameters
    assert "thin_counts" in inspect.getsource(dr.thin_counts) and "ob.thin_counts" in inspect.getsource(dr.thin_counts)


@pytest.mark.parametrize("world", ["ref67", "ref67b", "ref50"])
def test_the_thinned_reference_keeps_the_cells_and_their_order_and_reaches_both_forms(dr, world):
    ctx, _, _, _, full, ref_stage, genes = _world()
    inputs = dict(ctx.inputs)
    released = core.counts_of(inputs["ref"])
    got, stage, rec = dr.thinned_world(world, inputs, "ref", full, ref_stage, genes)
    pool = got["ref"]
    assert got["start"] is inputs["start"]  # the start stage is the released one
    assert pool is not inputs["ref"] and inputs["ref"].layers["counts"] is ctx.inputs["ref"].layers["counts"]
    assert list(pool.obs_names) == list(inputs["ref"].obs_names) and list(pool.var_names) == list(full.var_names)
    assert np.array_equal(pool.obs["celltype"].to_numpy(), inputs["ref"].obs["celltype"].to_numpy())
    thinned = core.counts_of(pool)
    assert (thinned != dr.thin_counts(released, world)).nnz == 0
    assert (thinned.toarray() <= released.toarray()).all() and (thinned != released).nnz > 0
    assert pool.layers["counts"].dtype == inputs["ref"].layers["counts"].dtype
    assert dr._largest_difference(pool.X, dr.prep.log1p_cp10k(thinned)) == 0.0
    assert stage.n == ref_stage.n and np.array_equal(stage.celltype, ref_stage.celltype)
    direct = dr.ob.as_scored(thinned, full.obs[["celltype"]], list(full.var_names), genes)
    assert stage.X.tobytes() == direct.X.tobytes()
    assert not np.array_equal(stage.X, ref_stage.X)
    p, stream = THIN[world]
    assert rec["p"] == p and rec["stream"] == stream and rec["world"] == world
    assert abs(rec["depth_ratio"] - p) < 0.01 and rec["depth_ratio"] == pytest.approx(thinned.sum() / released.sum())
    assert rec["rebuild_largest_difference"] <= dr.REBUILD_TOL and rec["any_cell_emptied"] is False
    assert not re.search(r"\b(n|cells|count)\b", " ".join(rec))  # ratios only; no count in the record
    again, stage2, rec2 = dr.thinned_world(world, dict(ctx.inputs), "ref", full, ref_stage, genes)
    assert (core.counts_of(again["ref"]) != thinned).nnz == 0 and stage2.X.tobytes() == stage.X.tobytes()
    assert rec2 == rec


def test_the_two_draws_at_two_thirds_differ(dr):
    ctx, _, _, _, full, ref_stage, genes = _world()
    a, sa, _ = dr.thinned_world("ref67", dict(ctx.inputs), "ref", full, ref_stage, genes)
    b, sb, _ = dr.thinned_world("ref67b", dict(ctx.inputs), "ref", full, ref_stage, genes)
    assert (core.counts_of(a["ref"]) != core.counts_of(b["ref"])).nnz > 0 and not np.array_equal(sa.X, sb.X)


def test_ref100_thins_nothing_and_reads_no_file(dr):
    ctx, _, _, _, _, ref_stage, genes = _world()

    def reader():
        raise AssertionError("ref100 must not read the reference file again")

    inputs = dict(ctx.inputs)
    got, stage, rec = dr.world_reference("ref100", inputs, "ref", ref_stage, genes, reader)
    assert got is inputs and stage is ref_stage and rec is None
    with pytest.raises(ValueError, match="unknown world"):
        dr.world_reference("ref33", inputs, "ref", ref_stage, genes, reader)


def test_values_that_do_not_rebuild_are_refused(dr):
    ctx, _, _, _, full, ref_stage, genes = _world()
    scaled = full.copy()
    scaled.X = sparse.csr_matrix(scaled.X.toarray() * 1.01)
    off = Stage.from_anndata(scaled[:, genes].copy(), genes)
    with pytest.raises(ValueError, match="do not give back the matrices the route reads"):
        dr.thinned_world("ref67", dict(ctx.inputs), "ref", full, off, genes)
    other = full[::-1].copy()
    with pytest.raises(ValueError, match="cell types are not the input pool's"):
        dr.thinned_world("ref67", dict(ctx.inputs), "ref", other, ref_stage, genes)
    with pytest.raises(ValueError, match="cells or genes are not the input pool's"):
        dr.thinned_world("ref67", dict(ctx.inputs), "ref", full[:-1].copy(), ref_stage, genes)


# ---------------------------------------------------------------------------------------------------------------------
# Where the route reads the reference (spies)


def _route(dr, monkeypatch, world):
    """One seed of the world's route on the synthetic inputs, with spies on every place the reference is read: the
    floor's anchors, every arm's build (t1_scope.build_submission), g_mom's fit and every arm's scoring."""
    ctx, groups, dirs, held, full, ref_stage, genes = _world()
    harness._G_MOM_PREP.clear()
    released, start = core.counts_of(ctx.inputs["ref"]), core.counts_of(ctx.inputs["start"])
    inputs, stage, _ = dr.world_reference(world, dict(ctx.inputs), "ref", ref_stage, genes, lambda: full)
    seen = {"anchors": [], "score": [], "build_ref": [], "build_start": [], "fit_ref": [], "fit_start": []}
    real_build, real_fit = t1_scope.build_submission, ts.fit_state_factors

    def floor_spy(t, board, truth, ref, seed):
        seen["anchors"].append(ref)
        return SimpleNamespace(anchors={}, score=50.0)

    def build_spy(ref, *a, **kw):
        seen["build_ref"].append(core.counts_of(ref))
        if kw.get("start") is not None:
            seen["build_start"].append(core.counts_of(kw["start"]))
        return real_build(ref, *a, **kw)

    def fit_spy(early, late, *a, **kw):
        seen["fit_start"].append(core.counts_of(early))
        seen["fit_ref"].append(core.counts_of(late))
        return real_fit(early, late, *a, **kw)

    def score_spy(board, truth, ref, seed, anchors, preds):
        seen["score"].append(ref)
        return {k: {"score": 0.0} for k in preds}

    monkeypatch.setattr(dr.gv, "floor_anchors", floor_spy)
    monkeypatch.setattr(t1_scope, "build_submission", build_spy)
    monkeypatch.setattr(ts, "fit_state_factors", fit_spy)
    monkeypatch.setattr(dr.er, "score_arms", score_spy)
    truth = Stage(X=ref_stage.X[:50].copy(), celltype=ref_stage.celltype[:50])
    floors, failures = dr.ls.all_anchors(None, ctx.board, truth, stage, [0])
    rows = dr.seed_rows(ctx.task, ctx.board, ctx.board, truth, stage, inputs, floors, [0], dirs, groups, 0.2, held,
                        held)
    return seen, stage, ref_stage, core.counts_of(inputs["ref"]), released, start, rows, failures


def test_every_place_the_route_reads_the_reference_gets_the_thinned_one(dr, monkeypatch):
    seen, stage, ref_stage, thinned, released, start, rows, failures = _route(dr, monkeypatch, "ref67")
    assert not failures and len(rows) == 1 and set(rows[0]["arms"]) == set(ALL_ARMS)
    assert (thinned != released).nnz > 0  # a spy that saw the released reference would fail below
    assert seen["anchors"] and all(r is stage for r in seen["anchors"])
    assert len(seen["score"]) == 1 and seen["score"][0] is stage and stage is not ref_stage
    assert len(seen["build_ref"]) == 1 + 2 * len(MULT) + 2  # b0, then m1_r and m1_r_gmom for each of the seven arms
    assert all((c != thinned).nnz == 0 for c in seen["build_ref"])
    assert seen["build_start"] and all((c != start).nnz == 0 for c in seen["build_start"])  # B0's start: released
    assert len(seen["fit_ref"]) == 1  # g_mom's fit: once, on the thinned pool and the released start
    assert (seen["fit_ref"][0] != thinned).nnz == 0 and (seen["fit_start"][0] != start).nnz == 0
    assert set(rows[0]["sparsity"]) == set(ALL_ARMS) and set(rows[0]["sha256"]) == set(ALL_ARMS)


def test_in_ref100_every_place_sees_the_released_reference(dr, monkeypatch):
    seen, stage, ref_stage, counts, released, start, _, _ = _route(dr, monkeypatch, "ref100")
    assert stage is ref_stage and (counts != released).nnz == 0
    assert all(r is ref_stage for r in seen["anchors"] + seen["score"])
    assert all((c != released).nnz == 0 for c in seen["build_ref"] + seen["fit_ref"])


def test_g_moms_cache_refits_on_the_thinned_pool(dr):
    """harness._g_mom_factors caches its fit by the inputs' identity: a fit on the released pool is not reused for the
    thinned one in the same process."""
    ctx, groups, dirs, held, full, ref_stage, genes = _world()
    harness._G_MOM_PREP.clear()
    dr.build_seed(ctx, harness._b0_t1(harness.B0_DAMPING), dirs, groups, 0.2, held, held)
    first = harness._G_MOM_PREP["fold"]["factors"]
    inputs, _, _ = dr.world_reference("ref50", dict(ctx.inputs), "ref", ref_stage, genes, lambda: full)
    ctx2 = harness.BuildContext(task=ctx.task, board=ctx.board, inputs=inputs, seed=ctx.seed, derived=ctx.derived)
    dr.build_seed(ctx2, harness._b0_t1(harness.B0_DAMPING), dirs, groups, 0.2, held, held)
    assert harness._G_MOM_PREP["fold"]["factors"] is not first
    assert harness._G_MOM_PREP["fold"]["inputs"][0] is inputs["ref"]


# ---------------------------------------------------------------------------------------------------------------------
# The arms


def test_the_arms_are_t50s_up_to_six_times_the_step(dr):
    assert dr.MULTIPLIERS == MULT and dr.ARMS == ALL_ARMS
    fr = dr.arm_fractions(0.2)
    assert fr["m1_r_gmom"] == 0.2
    assert [fr[a] for a in SCALED] == pytest.approx([0.3, 0.4, 0.6, 0.8, 1.0, 1.2])


def test_ref100s_arms_are_t50s_bit_for_bit(dr):
    """b0, m1_r_gmom and the six scaled arms, built as try T50 builds them on the released inputs: the same bytes as
    over_read.build_seed's (and T46's pinned table for the arms it pins)."""
    ctx, groups, dirs, held = _arms_world()
    b0 = harness._b0_t1(harness.B0_DAMPING)
    mine = dr.build_seed(ctx, b0, dirs, groups, 0.2, held, held)
    t50 = dr.ob.build_seed(ctx, b0, dirs, groups, 0.2, held, held)
    assert list(mine) == list(ALL_ARMS)
    assert {a: dr.ls.prediction_sha256(p) for a, p in mine.items()} == \
        {a: dr.ls.prediction_sha256(t50[a]) for a in ALL_ARMS}
    assert {a: dr.ls.prediction_sha256(mine[a]) for a in GOLDEN_ARMS_FD18} == GOLDEN_ARMS_FD18


# ---------------------------------------------------------------------------------------------------------------------
# The reading


BASE = 53.0
W = {"de_score": 0.25, "de_direction": 0.25, "mmd_u": 0.30, "variogram": 0.20}
ROWS, CONTAM = 3431, 13


def _arm(score, de_part=0.5):
    gain = score - BASE
    sk = {"de_score": 0.4 + de_part * gain / 2 / (100 * W["de_score"]),
          "de_direction": 0.4 + de_part * gain / 2 / (100 * W["de_direction"]),
          "mmd_u": 0.3 + (1 - de_part) * gain * 0.6 / (100 * W["mmd_u"]),
          "variogram": 0.2 + (1 - de_part) * gain * 0.4 / (100 * W["variogram"])}
    return {"score": score, "groups": {"distribution": score / 100},
            "metrics": {"mmd_u": 0.1, "variogram": 0.01, "de_score": 0.5 + gain / 100, "de_direction": 0.3},
            "skills": sk}


def _task(world, gains=None, per_seed=None, seeds=range(10), commit="c" * 40, clean=True, fractions=None,
          thinning=None, anchor_failures=None, truth="lib3", task="T1/E9.5->E10.5:chamber", name=None, same=False):
    """A task file. gains: arm -> gain over m1_r_gmom on every seed; per_seed: arm -> ten gains. Each world's arms are
    offset by a tiny world-specific amount on seed 0 (none in ref100), unless `same`."""
    gains = gains or {}
    off = 0.0 if same else 1e-9 * WORLDS.index(world)
    out = []
    for i, s in enumerate(seeds):
        d = off if i == 0 else 0.0
        arms = {"b0": _arm(50.0 + d), "m1_r_gmom": _arm(BASE + d)}
        for a in SCALED:
            g = per_seed[a][i] if per_seed and a in per_seed else gains.get(a, 0.0)
            arms[a] = _arm(BASE + g + d)
        out.append({"seed": s, "floor": 50.0, "arms": arms, "sha256": {a: f"{world}-{a}-{s}" for a in arms},
                    "sparsity": {a: {"zero_share": 0.5 + 0.01 * j, "detected_per_cell": 100.0 - j}
                                 for j, a in enumerate(arms)}})
    rec = {"label": "synthetic", "task": task, "world": world if name is None else name, "truth": truth,
           "truth_rows": ROWS, "truth_contamination": CONTAM, "fraction": 0.2,
           "fractions": fractions or {"m1_r_gmom": 0.2, **{a: 0.2 * m for a, m in MULT.items()}},
           "multipliers": dict(MULT), "seeds": list(seeds), "rows": out, "thinning": None,
           "stamp": {"commit": commit, "tree_clean": clean, "utc": "x", "inputs": {}}}
    if world in THIN:
        p, stream = THIN[world]
        rec["thinning"] = {"world": world, "p": p, "stream": list(stream), "depth_ratio": p,
                           "rebuild_largest_difference": 0.0, "any_cell_emptied": False, **(thinning or {})}
    elif thinning is not None:
        rec["thinning"] = thinning
    if anchor_failures is not None:
        rec["anchor_failures"] = anchor_failures
    return rec


def _write(tmp_path, worlds, committed_shift=None):
    """worlds: world -> task record; try T50's committed lib3 file beside them, ref100's scores (and an x8 arm)."""
    out = tmp_path / "out"
    (out / "tasks").mkdir(parents=True)
    for k, rec in worlds.items():
        (out / "tasks" / f"{k}.json").write_text(json.dumps(rec))
    src = worlds.get("ref100") or _task("ref100")
    shift = committed_shift or {}
    rows = [{"seed": r["seed"], "arms": {**{a: {**v, "score": v["score"] + shift.get(a, 0.0)}
                                            for a, v in r["arms"].items()}, "m1_r_gmom_x8": _arm(BASE)}}
            for r in src["rows"]]
    t50 = tmp_path / "t50_lib3.json"
    t50.write_text(json.dumps({"truth": "lib3", "rows": rows}))
    return out, t50


def _read(dr, tmp_path, by_world=None, over=None, committed_shift=None):
    """by_world: world -> gains or {"per_seed": ...}; over: world -> task keyword overrides (None drops the world)."""
    by_world, over = by_world or {}, over or {}
    worlds = {}
    for k in WORLDS:
        if k in over and over[k] is None:
            continue
        g = by_world.get(k, {})
        args = {"per_seed": g["per_seed"]} if "per_seed" in g else {"gains": g}
        worlds[k] = _task(k, **args, **over.get(k, {}))
    out, t50 = _write(tmp_path, worlds, committed_shift)
    return dr.read(out, t50)["reading"]


def _prints_a_count(txt: str) -> bool:
    body = txt.rsplit("\nCommit ", 1)[0]
    return any(re.search(rf"(?<![\d.]){n}(?![\d.])", body) for n in (ROWS, CONTAM))


def _curve(top):
    """Mean gains over m1_r_gmom rising with the multiplier up to `top` and falling after it."""
    return {a: (m if m <= top else 2 * top - m) for a, m in MULT.items()}


@pytest.mark.parametrize("gains, against_x3, x3x2, x2x1, best, statement", [
    (_curve(6.0), {X4: "ahead", X5: "ahead", X6: "ahead"}, "ahead", "ahead", 6.0, "still rising at 6"),
    (_curve(3.0), {X4: "behind", X5: "behind", X6: "behind"}, "ahead", "ahead", 3.0, "turns at 3"),
    ({**_curve(4.0), X4: 3.2, X5: 3.1}, {X4: "level", X5: "level", X6: "behind"}, "ahead", "ahead", 4.0,
     "turns at 4"),
    ({X2: 0.1, X3: 0.2, X4: 0.25, X5: 0.2, X6: -0.2, "m1_r_gmom_x15": 0.05},
     {X4: "level", X5: "level", X6: "behind"}, "level", "level", 4.0, "turns at 4"),
    ({X2: -1.0, X3: -2.0, X4: -3.0, X5: -4.0, X6: -5.0, "m1_r_gmom_x15": -0.5},
     {X4: "behind", X5: "behind", X6: "behind"}, "behind", "behind", 1.0, "turns at 1"),
])
def test_statement_one_each_label_and_the_argmax(dr, tmp_path, gains, against_x3, x3x2, x2x1, best, statement):
    r = _read(dr, tmp_path, {k: gains for k in WORLDS})
    assert r["valid"], r["validity"]["reasons"]
    for k in WORLDS:
        s1 = r["statement_1"][k]
        assert s1["against_x3"] == against_x3 and s1["x3_against_x2"] == x3x2 and s1["x2_against_x1"] == x2x1
        assert s1["argmax"] == best and s1["statement"] == statement
    txt = (tmp_path / "out" / "t54.txt").read_text()
    assert f"in ref67: against x3 x4 {against_x3[X4]}" in txt and f"argmax {best:g}" in txt
    assert not _prints_a_count(txt)


def test_the_tie_rule_picks_the_smallest_multiplier(dr):
    means = {"m1_r_gmom": 53.0, "m1_r_gmom_x15": 54.0, X2: 55.0, X3: 56.0, X4: 57.0, X5: 57.0, X6: 56.5}
    assert dr.turn(means)["best"] == 4.0
    assert dr.turn({**means, X6: 57.0, X4: 56.0, X5: 56.0})["statement"] == "still rising at 6"


def test_a_gain_of_exactly_point_three_is_level_and_eight_seeds_are_needed(dr, tmp_path):
    exact = {X3: [0.0] * 10, X4: [0.25] * 8 + [0.5] * 2}  # x4 - x3: mean exactly 0.3, ten seeds ahead
    r = _read(dr, tmp_path, {"ref67": {"per_seed": exact}})
    c = r["statement_2"]["ref67"][f"{X4} against {X3}"]
    assert r["valid"] and c["mean"] == 0.3 and c["seeds_ahead"] == 10 and c["label"] == "level"
    r = _read(dr, tmp_path / "b", {"ref67": {"per_seed": {**exact, X4: [0.25] * 8 + [0.75] * 2}}})
    assert r["statement_2"]["ref67"][f"{X4} against {X3}"]["label"] == "ahead"
    r = _read(dr, tmp_path / "c", {"ref50": {"per_seed": {X3: [0.0] * 10, X4: [1.0] * 7 + [-0.1] * 3}}})
    c = r["statement_2"]["ref50"][f"{X4} against {X3}"]
    assert c["mean"] > 0.3 and c["seeds_ahead"] == 7 and c["label"] == "level"


def test_statement_two_reports_each_quoted_case_with_its_label_and_parts(dr, tmp_path):
    g = {X2: 1.0, X3: 2.0, X4: 2.2, X5: 1.5, X6: 1.0}
    r = _read(dr, tmp_path, {k: g for k in WORLDS})
    s2 = r["statement_2"]["ref50"]
    assert list(s2) == [f"{X4} against {X3}", f"{X5} against {X3}", f"{X4} against {X2}", f"{X5} against {X2}"]
    want = {f"{X4} against {X3}": (0.2, "level"), f"{X5} against {X3}": (-0.5, "behind"),
            f"{X4} against {X2}": (1.2, "ahead"), f"{X5} against {X2}": (0.5, "ahead")}
    for name, (mean, lab) in want.items():
        c = s2[name]
        assert c["mean"] == pytest.approx(mean, abs=1e-8) and c["label"] == lab
        assert c["seeds_ahead"] == (0 if mean < 0 else 10) and len(c["range"]) == 2
        assert c["points"]["de"] == pytest.approx(mean / 2, abs=1e-8)
        assert c["points"]["mmd_u_variogram"] == pytest.approx(mean / 2, abs=1e-8)
    assert set(r["statement_2"]) == set(WORLDS)
    txt = (tmp_path / "out" / "t54.txt").read_text()
    assert f"in ref50: {X5} against {X3}: -0.50 points, ahead 0 of 10" in txt
    assert "test: behind; de terms -0.25, mmd_u and variogram -0.25; de_score skill" in txt


def _depth(d67, d50, null=0.0, base=1.0):
    """(x4 - x2) per world: ref100's `base`, ref67's base + d67, ref50's base + d50, ref67b's ref67's + null."""
    def g(d):
        return {X2: 1.0, X3: 1.5, X4: 1.0 + d, X5: 0.5 + d, X6: 0.4 + d}
    return {"ref100": g(base), "ref67": g(base + d67), "ref50": g(base + d50), "ref67b": g(base + d67 + null)}


@pytest.mark.parametrize("d67, d50, statement", [
    (-1.0, -1.0, "The larger steps gain less when the truth is deeper than the reference"),
    (1.0, 1.0, "The larger steps gain more when the truth is deeper than the reference"),
    (0.1, -0.2, "no difference at the bar"),
    (-1.0, 0.1, "mixed (D67 behind, D50 level)"),
    (1.0, -1.0, "mixed (D67 ahead, D50 behind)"),
])
def test_statement_three_each_value(dr, tmp_path, d67, d50, statement):
    r = _read(dr, tmp_path, _depth(d67, d50))
    s3 = r["statement_3"]
    assert r["valid"] and s3["valid"] and s3["draw_null"]["holds"] and s3["draw_null"]["label"] == "level"
    assert s3["statement"] == statement and s3["label_if_valid"] == statement
    assert s3["D67"]["mean"] == pytest.approx(d67, abs=1e-8) and s3["D50"]["mean"] == pytest.approx(d50, abs=1e-8)
    assert s3["D67"]["quantity"] == f"({X4} - {X2}) on ref67 - ({X4} - {X2}) on ref100"
    assert r["outcome"].endswith(f"; depth: {statement}")
    assert set(s3["shapes"]) == set(WORLDS) and s3["shapes"]["ref50"] == r["statement_1"]["ref50"]
    txt = (tmp_path / "out" / "t54.txt").read_text()
    assert f"Statement 3, depth (deep side): {statement}" in txt and "statement 3 stands" in txt
    assert "D67 on the raw metrics (read by nothing): de_score" in txt


def test_d67_and_d50_are_paired_by_seed(dr, tmp_path):
    per = {a: [0.0] * 10 for a in SCALED}
    base = {**per, X2: [1.0] * 10, X4: [float(i) for i in range(10)]}
    r67 = {**per, X2: [1.0] * 10, X4: [float(i) - 0.5 for i in range(10)]}
    r50 = {**per, X2: [1.0] * 10, X4: [float(i) - 0.1 for i in range(10)]}
    r = _read(dr, tmp_path, {"ref100": {"per_seed": base}, "ref67": {"per_seed": r67}, "ref50": {"per_seed": r50},
                             "ref67b": {"per_seed": r67}})
    s3 = r["statement_3"]
    assert s3["D67"]["mean"] == pytest.approx(-0.5) and s3["D67"]["seeds_ahead"] == 0 and s3["D67"]["sd"] < 1e-8
    assert s3["D67"]["label"] == "behind" and s3["D50"]["label"] == "level"
    assert s3["statement"] == "mixed (D67 behind, D50 level)"


@pytest.mark.parametrize("null, label", [(0.5, "ahead"), (-0.5, "behind")])
def test_statement_three_is_void_when_the_draw_null_is_not_level(dr, tmp_path, null, label):
    r = _read(dr, tmp_path, _depth(-1.0, -1.0, null=null))
    assert r["valid"] and "statement_1" in r and "statement_2" in r
    s3 = r["statement_3"]
    assert not s3["valid"] and not s3["draw_null"]["holds"] and s3["draw_null"]["label"] == label
    assert s3["draw_null"]["mean"] == pytest.approx(null, abs=1e-8)
    assert s3["statement"] == dr.DEPTH_VOID and s3["label_if_valid"] == dr.LESS
    assert s3["D67"]["label"] == "behind"  # still reported
    txt = (tmp_path / "out" / "t54.txt").read_text()
    assert "statement 3 is void and only reported" in txt and dr.DEPTH_VOID in txt


def test_the_reported_quantities(dr, tmp_path):
    anchors = {m: {"floor": 0.1, "ceiling": 0.6} for m in W}
    worlds = {k: _task(k, gains=_curve(6.0)) for k in WORLDS}
    for row in worlds["ref67"]["rows"]:
        row["anchors"] = anchors
    out, t50 = _write(tmp_path, worlds)
    r = dr.read(out, t50)["reading"]
    assert r["valid"]
    n = r["neighbours"]["ref67"]
    assert n["m1_r_gmom_x15"]["base"] == "m1_r_gmom" and n[X6]["base"] == X5
    assert n["m1_r_gmom_x15"]["per_tenth"] == pytest.approx(1.5) and n[X5]["per_tenth"] == pytest.approx(0.5)
    assert r["mean_scores"]["ref50"][X6] == pytest.approx(BASE + 6, abs=1e-6)
    assert r["sparsity"]["ref67"]["b0"] == pytest.approx({"zero_share": 0.5, "detected_per_cell": 100.0})
    assert r["statement_3"]["anchor_gap"]["ref67"] == pytest.approx({m: 0.5 for m in W})
    assert r["statement_3"]["anchor_gap"]["ref100"] is None
    assert r["thinning"]["ref50"]["depth_ratio"] == 0.5
    txt = (out / "t54.txt").read_text()
    assert "world ref50: realised depth ratio 0.5000 (p 0.5000, stream [54, 2])" in txt
    assert f"{X6} against {X5} +1.00 (ahead), +0.50 per 0.1 of fraction" in txt
    assert "scored zero share and detected scored genes per cell: b0 0.5000 / 100.00" in txt
    assert "ref67: mean ceiling minus floor per metric (read by nothing): de_score +0.5" in txt
    for line in dr.LIMITS:
        assert f"Limit: {line}" in txt
    assert not _prints_a_count(txt)


def test_the_limits_are_the_briefs(dr):
    joined = " ".join(dr.LIMITS)
    for want in ("a thinned reference is not a deeper truth", "changes the floor (the copy of the reference)",
                 "changes g_mom's fit", "board points are not comparable across worlds",
                 "statement 3 is read with the raw metrics and the shapes", "the chamber cut and one library only",
                 "one thinning draw per world, with a single second draw as self-check",
                 "nothing here says how deep any hidden truth is",
                 "no outcome is evidence on task 1's read depth, which stays closed"):
        assert want in joined.lower(), want


FAILURE = [{"seed": 3, "reason": "T1/E9.5->E10.5:chamber seed 3: inverted"}]
FR_X5 = {"m1_r_gmom": 0.2, **{a: 0.2 * m for a, m in MULT.items()}, X5: 0.8}


@pytest.mark.parametrize("over, shift, reason", [
    ({"ref50": None}, None, "missing worlds ['ref50']"),
    ({"ref67b": None}, None, "missing worlds ['ref67b']"),
    ({"ref100": None}, None, "missing worlds ['ref100']"),
    ({"ref67": {"seeds": range(9)}}, None, "ref67: seeds"),
    ({}, {"m1_r_gmom_x3": 1e-6}, "do not reproduce try T50's committed lib3 scores"),
    ({}, {"b0": 1e-6}, "do not reproduce try T50's committed lib3 scores"),
    ({"ref50": {"fractions": FR_X5}}, None, "wrong fraction on ['m1_r_gmom_x5']"),
    ({"ref67": {"commit": "d" * 40}}, None, "2 commits"),
    ({"ref67b": {"clean": False}}, None, "tree not clean on ['ref67b']"),
    ({"ref100": {"clean": None}}, None, "tree not clean"),
    ({"ref67": {"thinning": {"depth_ratio": 2 / 3 + 0.02}}}, None, "ref67: realised depth ratio"),
    ({"ref50": {"thinning": {"depth_ratio": 0.5 - 0.02}}}, None, "ref50: realised depth ratio"),
    ({"ref67b": {"thinning": {"any_cell_emptied": True}}}, None, "ref67b: the thinning left a reference cell"),
    ({"ref50": {"thinning": {"any_cell_emptied": None}}}, None, "ref50: the thinning left a reference cell"),
    ({"ref67": {"thinning": {"rebuild_largest_difference": 1e-4}}}, None, "ref67: the rebuild check's largest"),
    ({"ref67": {"thinning": {"stream": [54, 3]}}}, None, "ref67: the file records p"),
    ({"ref50": {"thinning": {"p": 2 / 3}}}, None, "ref50: the file records p"),
    ({"ref67": {"same": True}}, None, "ref67: arms score exactly as in ref100 on every seed"),
    ({"ref50": {"same": True}}, None, "ref50: arms score exactly as in ref100 on every seed"),
    ({"ref67": {"anchor_failures": FAILURE, "seeds": ()}}, None, "anchors invert or go missing on seeds [3]"),
    ({"ref100": {"thinning": {"p": 1.0}}}, None, "ref100: the file records a thinning"),
    ({"ref67": {"name": "ref50"}}, None, "ref67: the file holds world 'ref50'"),
    ({"ref67": {"truth": "lib3_t67"}}, None, "ref67: the file holds truth 'lib3_t67'"),
])
def test_the_void_cases(dr, tmp_path, over, shift, reason):
    r = _read(dr, tmp_path, {k: _curve(4.0) for k in WORLDS}, over=over, committed_shift=shift)
    assert not r["valid"] and r["outcome"] == "VOID"
    assert any(reason in x for x in r["validity"]["reasons"]), r["validity"]["reasons"]
    txt = (tmp_path / "out" / "t54.txt").read_text()
    assert "Void: " in txt and not _prints_a_count(txt)


def test_a_missing_ref100_seed_and_a_missing_committed_file_are_void(dr, tmp_path):
    worlds = {k: _task(k) for k in WORLDS}
    worlds["ref100"]["rows"] = worlds["ref100"]["rows"][:-1]
    out, t50 = _write(tmp_path, worlds)
    r = dr.read(out, t50)["reading"]
    assert not r["valid"] and any("ref100: seeds" in x for x in r["validity"]["reasons"])
    worlds = {k: _task(k) for k in WORLDS}
    out, t50 = _write(tmp_path / "b", worlds)
    r = dr.read(out, tmp_path / "b" / "absent.json")["reading"]
    assert not r["valid"] and any("do not reproduce try T50's" in x for x in r["validity"]["reasons"])


def test_within_tolerance_stays_valid(dr, tmp_path):
    r = _read(dr, tmp_path, over={"ref67": {"thinning": {"depth_ratio": 2 / 3 + 0.009}}},
              committed_shift={X4: 1e-12})
    assert r["valid"]


def test_the_reading_refuses_to_overwrite(dr, tmp_path):
    out, t50 = _write(tmp_path, {k: _task(k) for k in WORLDS})
    dr.read(out, t50)
    assert (out / "t54.json").exists() and (out / "t54.txt").exists()
    with pytest.raises(SystemExit, match="already holds"):
        dr.read(out, t50)


# ---------------------------------------------------------------------------------------------------------------------
# The run's refusals


def test_the_run_needs_its_row(dr, monkeypatch):
    monkeypatch.setattr(dr.n12, "tries_rows", lambda: {"T50", "T53"})
    monkeypatch.setattr(dr.gc, "n22_passed", lambda path: True)
    with pytest.raises(SystemExit, match="T54 has no row"):
        dr.main(["run", "--world", "ref67", "--out", "t54"])


def test_the_run_needs_n22(dr, monkeypatch):
    monkeypatch.setattr(dr.n12, "tries_rows", lambda: {"T54"})
    with pytest.raises(SystemExit, match="does not hold a passing N22"):
        dr.main(["run", "--world", "ref100", "--out", "t54", "--n22", "n22x"])


@pytest.mark.parametrize("world", WORLDS)
def test_the_run_refuses_an_existing_task_file(dr, monkeypatch, world):
    monkeypatch.setattr(dr.n12, "tries_rows", lambda: {"T54"})
    monkeypatch.setattr(dr.gc, "n22_passed", lambda path: True)
    tasks = dr.RESULTS / "t54" / "tasks"
    tasks.mkdir(parents=True)
    (tasks / f"{world}.json").write_text("{}")
    with pytest.raises(SystemExit, match="already holds"):
        dr.main(["run", "--world", world, "--out", "t54"])


def test_the_run_thins_before_anything_else_is_made(dr):
    src = inspect.getsource(dr.cmd_run)
    assert src.index("world_reference(") < src.index("all_anchors(") < src.index("seed_rows(")
    assert src.index("world_reference(") < src.index("feng_directions(")
    assert "ob.truth_and_contamination(t, data, genes, HEART)" in src
