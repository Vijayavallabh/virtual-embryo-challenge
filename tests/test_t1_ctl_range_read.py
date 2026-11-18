"""experiments/t1_m1/ctl_range_read.py (try T58 of experiments/harness/TRIES.md, brief AM) on synthetic data: the arms
(spies on the route and small synthetic builds), the worlds' routes, the reading on synthetic per-seed score tables and
the run's refusals. No released file, statistics file, board file or score is read."""

from __future__ import annotations

import importlib.util
import inspect
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from test_t1_deep_read import _world
from test_t1_lib_step_read import GOLDEN_ARMS_FD18
from test_t1_part_depth_read import _arm
from test_t1_step_read import _arms_world

from vec import core, harness
from vec.boards import REPO_ROOT
from vec.evaluate import Stage

SCRIPT = REPO_ROOT / "experiments" / "t1_m1" / "ctl_range_read.py"
STAMP = {"commit": "f" * 40, "tree_clean": True, "utc": "x", "inputs": {}}
WORLDS = ("pooled", "whole", "lib3", "lib2", "ref67", "ref50")
CHAMBER, WHOLE_CUT = "T1/E9.5->E10.5:chamber", "T1/E9.5->E10.5"
THIN = {"ref67": (2 / 3, [54, 1]), "ref50": (1 / 2, [54, 2])}
TRUTH = {"pooled": "pooled", "whole": "whole", "lib3": "lib3", "lib2": "lib2", "ref67": "lib3", "ref50": "lib3"}
BASE = "m1_r_gmom"
X = {k: f"{BASE}_x{k}" for k in (2, 3, 4, 5)}
P = {k: f"{BASE}_p{k}" for k in (2, 3, 4, 5)}
XB = {k: f"{BASE}_x{k}b2" for k in (3, 4, 5)}
ARMS = ("b0", BASE, *X.values(), *P.values(), *XB.values())
PL = {"CM": "pl-CM", "Endocardium": "pl-Endocardium"}  # a placebo record (group -> sha256)
FENG = {"CM": "feng-CM", "Endocardium": "feng-Endocardium"}  # Feng's record


@pytest.fixture
def cr(monkeypatch, tmp_path):
    spec = importlib.util.spec_from_file_location("ctl_range_read_module", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    monkeypatch.setattr(mod, "RESULTS", tmp_path / "results")
    monkeypatch.setattr(mod.gv, "stamp", lambda started, files: dict(STAMP))
    return mod


# ---------------------------------------------------------------------------------------------------------------------
# The arms and their fractions


def test_the_worlds_arms_and_fractions(cr):
    assert cr.WORLDS == WORLDS and cr.TRUTH_OF == TRUTH and cr.ARMS == ARMS
    assert cr.TASK_OF == {k: (WHOLE_CUT if k == "whole" else CHAMBER) for k in WORLDS}
    assert cr.PLACEBO_AT[2] == cr.pm.P2 == P[2]  # try T45's name
    assert cr.STEP_MULTIPLIERS == {X[k]: float(k) for k in (2, 3, 4, 5)}
    assert all(cr.STEP_MULTIPLIERS[a] == cr.ob.MULTIPLIERS[a] for a in cr.STEP_MULTIPLIERS)  # try T50's values
    fr = cr.arm_fractions(0.2)
    assert set(fr) == set(ARMS) - {"b0"}
    step = cr.sr.arm_fractions(0.2, {X[k]: float(k) for k in (2, 3, 4, 5)})
    for k in (2, 3, 4, 5):
        assert fr[X[k]] == step[X[k]] and fr[P[k]] == fr[X[k]]  # the same float as the step arm of its multiple
    assert fr[P[2]] == cr.pm.arm_fractions(0.2)[cr.pm.P2]  # try T45's p2 fraction, bit for bit
    for k, want in ((3, [0.6, -0.2]), (4, [0.8, -0.4]), (5, [1.0, -0.6])):
        first, second = fr[XB[k]]
        assert first == fr[X[k]]  # the first edit's fraction is the x_k arm's, bit for bit
        assert second == (2 - k) * 0.2  # back to twice the step
        assert first + second == pytest.approx(fr[X[2]], abs=1e-15)
        assert [first, second] == pytest.approx(want, abs=1e-15)
    assert cr.TWIN == {**{P[k]: X[k] for k in P}, **{XB[k]: X[k] for k in XB}}


def test_the_committed_arms_each_world_reproduces_and_the_arms_no_file_holds(cr):
    assert {s: cr.repro_arms(s) for s in ("t44", "t45", "t46", "t50", "t54")} == {
        "t44": ("b0", BASE, X[2]), "t45": ("b0", BASE, X[2], P[2]), "t46": ("b0", BASE, X[2], X[3], X[4]),
        "t50": ("b0", BASE, X[2], X[3], X[4], X[5]), "t54": ("b0", BASE, X[2], X[3], X[4], X[5])}
    assert cr.SOURCES == {"pooled": ("t45", "t46", "t50"), "whole": ("t45", "t44"), "lib3": ("t46", "t50"),
                          "lib2": ("t46", "t50"), "ref67": ("t54",), "ref50": ("t54",)}
    unheld = cr.unheld_arms()
    assert {k: len(v) for k, v in unheld.items()} == {"pooled": 6, "whole": 9, "lib3": 7, "lib2": 7, "ref67": 7,
                                                      "ref50": 7}
    assert unheld["whole"] == (X[3], X[4], X[5], P[3], P[4], P[5], *XB.values())
    f = cr.committed_files()
    assert {k: p.parts[-3:] for k, p in f["t45"].items()} == {
        "pooled": ("t45_m1_perm", "tasks", "T1_E9.5_E10.5_chamber.json"),
        "whole": ("t45_m1_perm", "tasks", "T1_E9.5_E10.5.json")}
    assert {k: p.parts[-3:] for k, p in f["t44"].items()} == {"whole": ("t44_m1_step", "tasks", "T1_E9.5_E10.5.json")}
    assert {k: p.parts[-1] for k, p in f["t46"].items()} == {"pooled": "pooled.json", "lib3": "lib3.json",
                                                             "lib2": "lib2.json"}
    assert {k: p.parts[-3:] for k, p in f["t50"].items()} == {
        k: ("t50_m1_over_step", "tasks", f"{k}.json") for k in ("pooled", "lib3", "lib2")}
    assert {k: p.parts[-3:] for k, p in f["t54"].items()} == {
        k: ("t54_m1_deep_step", "tasks", f"{k}.json") for k in ("ref67", "ref50")}


class _Ctx:
    def __init__(self, seed):
        self.seed = seed


def _spy(cr, monkeypatch):
    """Spies on gmom_combo.m1_r_arms (every build) and perm_read.second_edit (the there-and-back's second edit)."""
    calls, edits = [], []
    real_hook = cr.gc.gmom_hook

    def build(ctx, directions, groups, fraction, held_m1, held_factors):
        calls.append({"ctx": ctx, "directions": directions, "groups": groups, "fraction": fraction,
                      "held_m1": held_m1, "held_factors": held_factors, "hooked": cr.gc.gmom_hook is not real_hook})
        return {"m1_r": None, BASE: ("built", len(calls))}

    real_edit = cr.pm.second_edit

    def edit(directions, groups, fraction, held, seed):
        edits.append({"directions": directions, "groups": groups, "fraction": fraction, "held": held, "seed": seed,
                      "at_call": len(calls)})
        return real_edit(directions, groups, fraction, held, seed)

    monkeypatch.setattr(cr.gc, "m1_r_arms", build)
    monkeypatch.setattr(cr.pm, "second_edit", edit)
    return calls, edits


def test_each_new_arm_is_the_earlier_trys_build_with_only_its_fractions_changed(cr, monkeypatch):
    calls, edits = _spy(cr, monkeypatch)
    groups, held_m1, held_f = {"V-CM": "CM"}, ["a", "b"], ["a"]
    dirs, placebo = {"CM": np.zeros(3)}, {"CM": np.ones(3)}
    ctx = _Ctx(1000004)
    preds = cr.build_seed(ctx, lambda c: ("b0", c.seed), dirs, placebo, groups, 0.2, held_m1, held_f)
    assert list(preds) == list(ARMS)
    order = [BASE, *X.values(), *P.values(), *XB.values()]
    assert len(calls) == len(order)
    by = dict(zip(order, calls, strict=True))
    for c in calls:  # every build shares the context, groups and held genes: only directions and fraction move
        for key in ("ctx", "groups", "held_m1", "held_factors"):
            assert c[key] is calls[0][key], key
    fr = cr.arm_fractions(0.2)
    for a in (BASE, *X.values(), *XB.values()):
        assert by[a]["directions"] is dirs
    for a in P.values():
        assert by[a]["directions"] is placebo and by[a]["fraction"] == fr[a] and not by[a]["hooked"]
    for k in (2, 3, 4, 5):
        assert by[X[k]]["fraction"] == 0.2 * k and not by[X[k]]["hooked"]  # try T50's build
    # the there-and-back arms: try T45's code, the first edit at the x_k arm's fraction, the second back to x2's
    assert len(edits) == 3
    for (k, a), e in zip(XB.items(), edits, strict=True):
        assert by[a]["fraction"] == by[X[k]]["fraction"] == fr[a][0]  # the same float, bit for bit
        assert by[a]["hooked"]  # the second edit sits in the g_mom hook during the build
        assert e["directions"] is dirs and e["groups"] is groups and e["held"] is held_m1 and e["seed"] == ctx.seed
        assert e["fraction"] == fr[a][1] == (2 - k) * 0.2
        assert e["at_call"] == order.index(a)  # set up just before its own build
    assert preds[P[3]] == ("built", order.index(P[3]) + 1)
    assert inspect.getsource(cr.build_seed).count("pm.there_and_back(") == 1


def test_at_multiplier_two_the_placebo_build_is_try_t45s_p2_call(cr, monkeypatch):
    calls, _ = _spy(cr, monkeypatch)
    groups, held_m1, held_f = {"V-CM": "CM"}, ["a", "b"], ["a"]
    dirs, placebo, placebo2 = {"CM": np.zeros(3)}, {"CM": np.ones(3)}, {"CM": np.full(3, 2.0)}
    ctx = _Ctx(1000007)
    cr.pm.perm_arms(ctx, dirs, placebo, groups, 0.2, held_m1, held_f, placebo2=placebo2)
    t45 = dict(zip(cr.pm.M1_ARMS, calls, strict=True))[cr.pm.P2]
    calls.clear()
    cr.build_seed(ctx, lambda c: None, dirs, placebo, groups, 0.2, held_m1, held_f)
    mine = calls[1 + 4]  # after m1_r_gmom and x2 to x5
    assert mine["directions"] is placebo
    assert {k: v for k, v in mine.items() if k != "hooked"} == {k: v for k, v in t45.items() if k != "hooked"}
    assert all(mine[k] is t45[k] for k in ("ctx", "directions", "groups", "held_m1", "held_factors"))
    assert mine["fraction"] == t45["fraction"]


def _placebo_of(dirs):
    return {g: v[::-1].copy() for g, v in dirs.items()}


def test_the_arms_on_a_small_synthetic_build(cr):
    """The step arms are try T50's builds bit for bit (and the pinned golden ones), p2 is try T45's p2 bit for bit, the
    there-and-back arm's first edit is the x_k arm's, and every arm moves its prediction."""
    ctx, groups, dirs, held = _arms_world()
    placebo = _placebo_of(dirs)
    b0 = harness._b0_t1(harness.B0_DAMPING)
    mine = cr.build_seed(ctx, b0, dirs, placebo, groups, 0.2, held, held)
    sha = {a: cr.ls.prediction_sha256(p) for a, p in mine.items()}
    t50 = cr.ob.build_seed(ctx, b0, dirs, groups, 0.2, held, held)
    step = ("b0", BASE, *X.values())
    assert {a: sha[a] for a in step} == {a: cr.ls.prediction_sha256(t50[a]) for a in step}
    golden = ("b0", BASE, X[2], X[3], X[4])
    assert {a: sha[a] for a in golden} == {a: GOLDEN_ARMS_FD18[a] for a in golden}
    t45 = cr.pm.perm_arms(ctx, dirs, placebo, groups, 0.2, held, held, placebo2=_placebo_of(placebo))["arms"]
    assert sha[P[2]] == cr.ls.prediction_sha256(t45[cr.pm.P2])
    assert sha[X[2]] == cr.ls.prediction_sha256(t45[cr.pm.X2])
    assert len(set(sha.values())) == len(ARMS)


def test_the_idle_cases_are_their_twins_bit_for_bit(cr):
    """A placebo equal to Feng's directions gives the x_k arm of its fraction bit for bit; try T45's there-and-back
    with no step back (multiplier 1 from twice the step) is the x2 arm bit for bit."""
    ctx, groups, dirs, held = _arms_world()
    same = {g: v.copy() for g, v in dirs.items()}
    preds = cr.build_seed(ctx, lambda c: None, dirs, same, groups, 0.2, held, held)
    for k in (2, 3, 4, 5):
        assert preds[P[k]].X.tobytes() == preds[X[k]].X.tobytes(), k
    twice = cr.arm_fractions(0.2)[X[2]]
    idle = cr.pm.there_and_back(ctx, dirs, groups, twice, held, held, multiplier=1.0)
    assert idle.X.tobytes() == preds[X[2]].X.tobytes()


# ---------------------------------------------------------------------------------------------------------------------
# The worlds' routes


def _route(cr, monkeypatch, world, target_filter=True):
    ctx, groups, dirs, held, full, ref_stage, genes = _world()
    released = core.counts_of(ctx.inputs["ref"])
    truths = {"pooled": Stage(X=ref_stage.X[:70].copy(), celltype=ref_stage.celltype[:70]),
              "lib3": Stage(X=ref_stage.X[:50].copy(), celltype=ref_stage.celltype[:50]),
              "lib2": Stage(X=ref_stage.X[50:90].copy(), celltype=ref_stage.celltype[50:90])}
    whole = Stage(X=ref_stage.X[10:80].copy(), celltype=ref_stage.celltype[10:80])
    asked = []

    def t46_spy(t, data, genes_, name):
        asked.append(("t46", name))
        return truths[name], 3

    def filter_spy(t, data, genes_):
        asked.append(("filter", t.key))
        return whole, {"contamination": 7}

    def read_spy(path, genes_):
        asked.append(("read", str(path)))
        return full[:, list(genes_)].copy()

    monkeypatch.setattr(cr.ob.ls, "truth_and_contamination", t46_spy)
    monkeypatch.setattr(cr.harness, "filter_target", filter_spy)
    monkeypatch.setattr(cr, "read_stage", read_spy)
    t = SimpleNamespace(key=cr.TASK_OF[world], target="target.h5ad",
                        target_filter={"keep": ["x"]} if target_filter else None)
    original = dict(ctx.inputs)
    inputs, stage, truth, contamination, thinning, name = cr.world_route(world, t, Path("data"), genes, original,
                                                                         "ref", ref_stage, lambda: full)
    return SimpleNamespace(inputs=inputs, stage=stage, truth=truth, contamination=contamination, thinning=thinning,
                           name=name, asked=asked, truths=truths, whole=whole, original=original,
                           ref_stage=ref_stage, released=released, full=full, genes=genes)


@pytest.mark.parametrize("world", ["pooled", "lib3", "lib2"])
def test_the_chamber_worlds_take_try_t46s_truth_and_thin_nothing(cr, monkeypatch, world):
    r = _route(cr, monkeypatch, world)
    assert r.asked == [("t46", world)] and r.truth is r.truths[world] and r.name == world and r.contamination == 3
    assert r.thinning is None and r.inputs is r.original and r.stage is r.ref_stage


@pytest.mark.parametrize("world", ["ref67", "ref50"])
def test_the_thinned_worlds_are_try_t54s_with_library_3s_truth(cr, monkeypatch, world):
    r = _route(cr, monkeypatch, world)
    assert r.asked == [("t46", "lib3")] and r.truth is r.truths["lib3"] and r.name == "lib3"
    p, stream = THIN[world]
    want = cr.dr.thin_counts(r.released, world)  # try T54's own thinning under the world's own name
    assert (core.counts_of(r.inputs["ref"]) != want).nnz == 0 and (want != r.released).nnz > 0
    assert r.thinning["p"] == p and r.thinning["stream"] == stream and r.thinning["world"] == world
    assert r.inputs["start"] is r.original["start"] and r.stage is not r.ref_stage


def test_the_whole_world_loads_its_truth_as_tries_t44_and_t45(cr, monkeypatch):
    r = _route(cr, monkeypatch, "whole")
    assert r.asked == [("filter", WHOLE_CUT)] and r.truth is r.whole and r.contamination == 7 and r.name == "whole"
    assert r.thinning is None and r.inputs is r.original and r.stage is r.ref_stage
    r = _route(cr, monkeypatch, "whole", target_filter=False)
    assert r.asked == [("read", str(Path("data") / "target.h5ad"))] and r.contamination is None
    assert np.array_equal(np.asarray(r.truth.X), np.asarray(r.ref_stage.X))
    src = inspect.getsource(cr.cmd_run)
    assert src.index("world_route(") < src.index("all_anchors(") < src.index("world_placebo(") \
        < src.index("seed_rows(")
    with pytest.raises(ValueError, match="unknown world"):
        cr.world_route("ref67b", None, None, [], {}, "ref", None, lambda: None)


def test_the_placebo_is_try_t45s_on_the_directions_as_loaded(cr, monkeypatch):
    rng = np.random.default_rng(1)
    dirs = {g: rng.normal(0, 1, 40) for g in ("CM", "Endocardium")}
    halves = {g: {"mean_cp": rng.gamma(1.0, 5.0, 40)} for g in dirs}
    monkeypatch.setattr(cr.feng_dirs, "feng_directions", lambda genes: (dirs, {"sha256": "s"}))
    monkeypatch.setattr(cr.er, "feng_halves", lambda genes, d: (halves, "s"))
    directions, placebo, ext, stats = cr.world_placebo([f"g{i}" for i in range(40)])
    want = cr.pm.placebo_directions(dirs, halves)
    assert directions is dirs and stats == "s" and all(np.array_equal(placebo[g], want[g]) for g in dirs)
    monkeypatch.setattr(cr.er, "feng_halves", lambda genes, d: (halves, "other"))
    with pytest.raises(ValueError, match="statistics the halves came from"):
        cr.world_placebo([f"g{i}" for i in range(40)])


# ---------------------------------------------------------------------------------------------------------------------
# The reading on synthetic per-seed score tables


SCORE = {"b0": 50.0, BASE: 53.0, X[2]: 53.6, X[3]: 54.4, X[4]: 55.0, X[5]: 55.3,
         P[2]: 53.0, P[3]: 53.1, P[4]: 53.2, P[5]: 53.2, XB[3]: 53.7, XB[4]: 53.8, XB[5]: 53.9}


def _task(world, scores=None, per_seed=None, seeds=range(10), commit="c" * 40, clean=True, fractions=None,
          thinning=None, placebo=None, identical=(), anchor_failures=None, truth=None, drop=(), idle=()):
    """scores: arm -> mean score (over SCORE); per_seed: arm -> ten values added seed by seed; identical: arms that
    score exactly as their twin on every seed."""
    sc = {**SCORE, **(scores or {})}
    rows = []
    for i, s in enumerate(seeds):
        arms = {}
        for a, v in sc.items():
            if a in drop:
                continue
            arms[a] = _arm(v + 0.01 * s + (per_seed or {}).get(a, [0.0] * 10)[i])
        for a in identical:
            arms[a] = dict(arms[{**{P[k]: X[k] for k in P}, **{XB[k]: X[k] for k in XB}}[a]])
        for a in idle:  # a placebo arm left at p2's fraction
            arms[a] = dict(arms[P[2]])
        rows.append({"seed": s, "floor": 50.0, "arms": arms})
    fr = {BASE: 0.2, **{X[k]: 0.2 * k for k in X}, **{P[k]: 0.2 * k for k in P},
          **{XB[k]: [0.2 * k, (2 - k) * 0.2] for k in XB}}
    rec = {"label": "synthetic", "task": WHOLE_CUT if world == "whole" else CHAMBER, "world": world,
           "truth": truth or TRUTH[world], "truth_rows": 3431, "truth_contamination": 13, "seeds": list(seeds),
           "fraction": 0.2, "fractions": fractions or fr, "rows": rows, "thinning": None,
           "placebo_sha256": placebo if placebo is not None else {a: dict(PL) for a in P.values()},
           "feng_sha256": dict(FENG),
           "applied_size": {a: {"fraction": fr[a], "feng": {"CM": {"mean_abs_logf": 0.1, "sum_sq_logf": 2.0}},
                                "placebo": {"CM": {"mean_abs_logf": 0.1, "sum_sq_logf": 2.1}}} for a in P.values()},
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


REPRO = {"t44": ("b0", BASE, X[2]), "t45": ("b0", BASE, X[2], P[2]), "t46": ("b0", BASE, X[2], X[3], X[4]),
         "t50": ("b0", BASE, X[2], X[3], X[4], X[5]), "t54": ("b0", BASE, X[2], X[3], X[4], X[5])}
FOLDERS = {"t44": "t44_m1_step", "t45": "t45_m1_perm", "t46": "t46_m1_lib_step", "t50": "t50_m1_over_step",
           "t54": "t54_m1_deep_step"}


def _write(cr, tmp_path, worlds, shift=None, t45_placebo=None, skip=()):
    """worlds: world -> task record; the committed files of tries T44, T45, T46, T50 and T54 beside them (the world's
    own scores, so that they reproduce); shift: (try, world, arm) -> added to that committed score."""
    shift = shift or {}
    out = tmp_path / "out"
    (out / "tasks").mkdir(parents=True)
    for k, rec in worlds.items():
        (out / "tasks" / f"{k}.json").write_text(json.dumps(rec))
    dirs = {s: tmp_path / FOLDERS[s] / "tasks" for s in FOLDERS}
    files = cr.committed_files(*(dirs[s] for s in ("t44", "t45", "t46", "t50", "t54")))
    for src, by in files.items():
        dirs[src].mkdir(parents=True, exist_ok=True)
        for k, path in by.items():
            if (src, k) in skip:
                continue
            rows = [{"seed": r["seed"], "arms": {a: {**r["arms"][a], "score": r["arms"][a]["score"]
                                                     + shift.get((src, k, a), 0.0)} for a in REPRO[src]}}
                    for r in (worlds.get(k) or _task(k))["rows"]]
            rec = {"rows": rows}
            if src == "t45":
                pl = t45_placebo or PL
                rec["placebo_sha256"] = {"m1_r_gmom_p1": dict(pl), "m1_r_gmom_p2": dict(pl),
                                         "m1_r_gmom_q1": {"CM": "q"}, "m1_r_gmom_q2": {"CM": "q"}}
            path.write_text(json.dumps(rec))
    return out, dirs


def _read(cr, tmp_path, over=None, shift=None, t45_placebo=None, skip=()):
    over = over or {}
    worlds = {k: _task(k, **over.get(k, {})) for k in WORLDS if not (k in over and over[k] is None)}
    out, dirs = _write(cr, tmp_path, worlds, shift, t45_placebo, skip)
    return cr.read(out, *(dirs[s] for s in ("t44", "t45", "t46", "t50", "t54")))["reading"]


def _text(tmp_path):
    return (tmp_path / "out" / "t58.txt").read_text()


NEEDS = "the gain beyond twice the step needs Feng's own genes"
ONE = "one placebo permutation, try T45's first"
NEEDS_OUT = f"{NEEDS} ({ONE})"
PLACEBO_TOO = "a placebo gains as well; Feng's direction more"
NO_GAIN = "no gain beyond twice the step to explain"


def test_the_gain_needs_fengs_own_genes_and_the_quantities(cr, tmp_path):
    r = _read(cr, tmp_path)
    assert r["valid"], r["validity"]["reasons"]
    assert r["outcome"] == NEEDS_OUT and r["statement_3"]["statement"] == NEEDS_OUT
    assert r["statement_3"]["kind"] == NEEDS and r["statement_3"]["r_not_ahead"] == []
    assert r["statement_3"]["p_ahead"] == []
    assert r["statement_3"]["readings"] == {f"{k} S_{m}": "ahead" for k in ("pooled", "whole", "lib3") for m in (3, 4)}
    q = r["quantities"]["ref50"]
    for k, (rk, pk, vk) in {3: (0.8, 0.1, 0.1), 4: (1.4, 0.2, 0.2), 5: (1.7, 0.2, 0.3)}.items():
        x = q[str(k)]
        assert x["R"]["mean"] == pytest.approx(rk) and x["R"]["change"] == X[k] and x["R"]["base"] == X[2]
        assert x["P"]["mean"] == pytest.approx(pk) and x["P"]["change"] == P[k] and x["P"]["base"] == P[2]
        assert x["S"]["mean"] == pytest.approx(rk - pk) and x["S"]["label"] == "ahead"
        assert x["V"]["mean"] == pytest.approx(vk) and x["V"]["change"] == XB[k] and x["V"]["base"] == X[2]
        assert x["N"]["mean"] == pytest.approx(rk - vk) and x["N"]["change"] == X[k] and x["N"]["base"] == XB[k]
        assert x["S"]["seeds_ahead"] == 10 and x["S"]["range"] == pytest.approx([rk - pk, rk - pk])
    # the parts on the de terms and on mmd_u with the variogram add up to S_k's mean (try T46's points by part)
    s = q["3"]["S"]
    assert s["points"]["de"] + s["points"]["mmd_u_variogram"] == pytest.approx(s["mean"], abs=1e-9)
    assert r["mean_scores"]["whole"][X[5]] == pytest.approx(55.3 + 0.045)
    txt = _text(tmp_path)
    assert f"Outcome: {NEEDS}" in txt and "lib2, S_5:" in txt and "whole, N_4:" in txt and "3431" not in txt


@pytest.mark.parametrize("world,m,per_seed,label", [
    ("pooled", 3, None, "level"),  # the placebo gains as much as the step at x3
    ("whole", 4, None, "level"),
    ("lib3", 4, "behind", "behind"),  # the placebo gains more than the step
])
def test_not_shown_names_each_of_the_six_that_is_level_or_behind(cr, tmp_path, world, m, per_seed, label):
    gain = SCORE[X[m]] - SCORE[X[2]] + (0.7 if per_seed else 0.0)
    over = {world: {"scores": {P[m]: SCORE[P[2]] + gain}}}
    r = _read(cr, tmp_path, over=over)
    assert r["valid"], r["validity"]["reasons"]
    assert r["outcome"] == f"not shown ({ONE}; {world} S_{m} {label})"
    assert r["statement_3"]["not_ahead"] == {f"{world} S_{m}": label}
    assert sum(v == "ahead" for v in r["statement_3"]["readings"].values()) == 5


def test_not_shown_names_every_reading_not_ahead(cr, tmp_path):
    level = {"scores": {P[3]: SCORE[P[2]] + 0.8, P[4]: SCORE[P[2]] + 1.4}}
    r = _read(cr, tmp_path, over={"pooled": level, "lib3": level})
    assert r["outcome"] == f"not shown ({ONE}; pooled S_3 level, pooled S_4 level, lib3 S_3 level, lib3 S_4 level)"


def test_seven_seeds_ahead_is_level(cr, tmp_path):
    per = [0.0] * 7 + [1.0] * 3  # the placebo's increment overtakes x3's on three seeds; R_3 stays ahead
    r = _read(cr, tmp_path, over={"whole": {"per_seed": {P[3]: per}}})
    assert r["quantities"]["whole"]["3"]["S"]["seeds_ahead"] == 7
    assert r["quantities"]["whole"]["3"]["S"]["mean"] == pytest.approx(0.4)
    assert r["quantities"]["whole"]["3"]["R"]["label"] == "ahead"
    assert r["outcome"] == f"not shown ({ONE}; whole S_3 level)"


def test_statement_3_reads_nothing_of_s_5_the_other_worlds_or_v_and_n(cr, tmp_path):
    flat = {"scores": {P[3]: SCORE[P[2]] + 0.8, P[4]: SCORE[P[2]] + 1.4, P[5]: SCORE[P[2]] + 1.7}}
    over = {"lib2": flat, "ref67": flat, "ref50": flat,
            "pooled": {"scores": {P[5]: SCORE[P[2]] + 1.7, XB[3]: SCORE[X[3]] + 0.5}}}
    r = _read(cr, tmp_path, over=over)
    assert r["valid"], r["validity"]["reasons"]
    assert r["outcome"] == NEEDS_OUT
    assert r["quantities"]["pooled"]["5"]["S"]["label"] == "level"
    assert r["quantities"]["lib2"]["3"]["S"]["label"] == "level"
    assert r["quantities"]["pooled"]["3"]["N"]["label"] == "behind"


def _incs(r3, p3, r4, p4):
    """Scores giving R_3, P_3, R_4, P_4 exactly (over x2 and p2) on every seed."""
    return {"scores": {X[3]: SCORE[X[2]] + r3, X[4]: SCORE[X[2]] + r4, P[3]: SCORE[P[2]] + p3,
                       P[4]: SCORE[P[2]] + p4}}


@pytest.mark.parametrize("r3,p3,r4,p4,r_label", [
    (0.0, -0.4, 0.0, -0.8, "level"),  # the reviewer's case 1: the step level, the placebo loses
    (-0.6, -1.5, -1.2, -2.5, "behind"),  # case 2: the step has turned, the placebo falls faster
])
def test_s_ahead_without_a_step_gain_is_no_gain_to_explain(cr, tmp_path, r3, p3, r4, p4, r_label):
    """S_k ahead with R_k level or behind is the expected dull result (a random direction of this size loses); it is
    never read as a gain that needs Feng's genes."""
    r = _read(cr, tmp_path, over={k: _incs(r3, p3, r4, p4) for k in ("pooled", "whole", "lib3")})
    assert r["valid"], r["validity"]["reasons"]
    s3 = r["statement_3"]
    assert all(v == "ahead" for v in s3["readings"].values())  # S_3 and S_4 ahead in all six, as the brief's bar reads
    assert s3["kind"] == NO_GAIN and NEEDS not in r["outcome"]
    assert s3["r_not_ahead"] == [f"{k} R_{m}" for k in ("pooled", "whole", "lib3") for m in (3, 4)]
    assert r["outcome"].startswith(f"{NO_GAIN} (pooled R_3 {r_label} {r3:+.2f} beside P_3 behind {p3:+.2f}, S_3 ahead: "
                                   "the placebo loses more; ")
    assert r["outcome"].count("the placebo loses more") == 6
    assert "Where R_k is negative (the step has turned) a positive S_k only says the placebo falls faster." \
        in _text(tmp_path)


def test_no_gain_names_only_the_readings_where_the_step_does_not_gain(cr, tmp_path):
    over = {"pooled": _incs(0.8, 0.1, 1.4, 0.2), "whole": _incs(0.8, 0.1, 0.1, -0.5), "lib3": _incs(0.8, 0.1, 1.4, 0.2)}
    r = _read(cr, tmp_path, over=over)
    assert r["outcome"] == (f"{NO_GAIN} (whole R_4 level +0.10 beside P_4 behind -0.50, S_4 ahead: "
                            "the placebo loses more)")
    over["whole"] = _incs(0.8, 0.1, 0.2, 0.1)  # R_4 level and S_4 level: no 'loses more'
    r = _read(cr, tmp_path / "b", over=over)
    assert r["outcome"] == f"{NO_GAIN} (whole R_4 level +0.20 beside P_4 level +0.10)"


def test_a_placebo_that_gains_as_well_is_named_with_its_means(cr, tmp_path):
    """The reviewer's case: R_3 +1.4 and P_3 +1.0, R_4 +2.0 and P_4 +1.5: S ahead, but most of the gain also comes
    with the placebo, so the outcome is not 'needs Feng's own genes'."""
    r = _read(cr, tmp_path, over={k: _incs(1.4, 1.0, 2.0, 1.5) for k in ("pooled", "whole", "lib3")})
    s3 = r["statement_3"]
    assert all(v == "ahead" for v in s3["readings"].values())
    assert s3["kind"] == PLACEBO_TOO and NEEDS not in r["outcome"]
    assert r["outcome"].startswith(f"{PLACEBO_TOO} ({ONE}; placebo ahead: pooled R_3 ahead +1.40 beside P_3 ahead "
                                   "+1.00; pooled R_4 ahead +2.00 beside P_4 ahead +1.50; whole R_3")
    assert len(s3["p_ahead"]) == 6
    # one P_k ahead in one world is enough to split
    r = _read(cr, tmp_path / "b", over={"lib3": _incs(1.4, 0.5, 1.4, 0.2)})
    assert r["outcome"] == f"{PLACEBO_TOO} ({ONE}; placebo ahead: lib3 R_3 ahead +1.40 beside P_3 ahead +0.50)"


def test_needs_fengs_genes_with_a_placebo_that_loses_and_a_step_that_gains(cr, tmp_path):
    r = _read(cr, tmp_path, over={k: _incs(0.8, -0.5, 1.4, -0.9) for k in ("pooled", "whole", "lib3")})
    assert r["outcome"] == NEEDS_OUT and r["statement_3"]["kind"] == NEEDS
    assert r["quantities"]["whole"]["4"]["P"]["label"] == "behind"


def test_statement_3_outcomes_directly(cr):
    """statement_3 on label tables alone: each of the four outcomes, with the precedence gate, split, not shown."""
    def q(r, p, s, rm=1.0, pm_=0.0):
        return {"R": {"label": r, "mean": rm}, "P": {"label": p, "mean": pm_}, "S": {"label": s, "mean": rm - pm_}}

    def table(**over):
        return {k: {m: over.get(f"{k}{m}", q("ahead", "level", "ahead")) for m in (3, 4, 5)}
                for k in ("pooled", "whole", "lib3", "lib2", "ref67", "ref50")}

    assert cr.statement_3(table())["kind"] == NEEDS
    assert cr.statement_3(table(lib25=q("behind", "ahead", "behind"), pooled5=q("level", "level", "level")))["kind"] \
        == NEEDS  # S_5, R_5 and the other worlds read nothing
    assert cr.statement_3(table(whole3=q("ahead", "ahead", "ahead")))["kind"] == PLACEBO_TOO
    assert cr.statement_3(table(whole3=q("ahead", "ahead", "level")))["kind"] == "not shown"
    assert cr.statement_3(table(whole3=q("ahead", "level", "level")))["statement"] == \
        f"not shown ({ONE}; whole S_3 level)"
    # the gate comes first: a level R_k gives no gain to explain whatever S_k and P_k read
    assert cr.statement_3(table(lib34=q("level", "behind", "ahead", 0.0, -1.0), whole3=q("ahead", "ahead", "level"))
                          )["kind"] == NO_GAIN


def test_the_limits_choices_statements_and_reported_numbers_are_in_the_text(cr, tmp_path):
    _read(cr, tmp_path)
    txt = _text(tmp_path)
    for x in (*cr.LIMITS, *cr.CHOICES, *cr.STATEMENTS):
        assert x in txt
    for x in ("One released day, in-sample.", "One placebo permutation (try T45's first).",
              "A thinned reference is one draw and is not a deeper truth.",
              "Nothing here says how deep any hidden truth is, and no statement may be used with a figure for it."):
        assert f"Limit: {x}" in txt
    assert "Arms no committed file holds, per world: pooled 6" in txt and "whole 9" in txt
    assert "the placebo's applied size beside Feng's, m1_r_gmom_p4 at fraction" in txt
    assert "ratio 1.000, 1.050" in txt
    assert "read by nothing" in txt and "world ref67: realised depth ratio" in txt
    assert "3431" not in txt and "truth_rows" not in txt


VOIDS = [
    ({"lib2": None}, None, "missing worlds ['lib2']"),
    ({"whole": {"seeds": range(9)}}, None, "whole: seeds"),
    ({"ref50": {"drop": (XB[5],)}}, None, "ref50: arms ['m1_r_gmom_x5b2'] missing"),
    ({"pooled": {"commit": "d" * 40}}, None, "the worlds ran at 2 commits"),
    ({"lib3": {"clean": False}}, None, "tree not clean on ['lib3']"),
    ({"ref67": {"clean": None}}, None, "tree not clean on ['ref67']"),
    ({"lib2": {"anchor_failures": [{"seed": 3}]}}, None, "lib2: the floor's anchors invert"),
    ({"lib2": {"truth": "lib3"}}, None, "lib2: the file holds truth 'lib3', not lib2"),
    ({"lib3": {"thinning": {"p": 0.5}}}, None, "lib3: the file records a thinning"),
    ({"ref67": {"thinning": {"depth_ratio": 0.6}}}, None, "ref67: realised depth ratio"),
    ({"ref50": {"thinning": {"stream": [54, 9]}}}, None, "ref50: the file records p"),
    ({"ref67": {"thinning": {"any_cell_emptied": True}}}, None, "ref67: the thinning left a cell"),
    ({"ref50": {"thinning": {"rebuild_largest_difference": 1e-3}}}, None, "ref50: the rebuild check's"),
    ({"whole": {"fractions": {BASE: 0.2, **{X[k]: 0.2 * k for k in X}, **{P[k]: 0.2 * k for k in P},
                              **{XB[k]: [0.2 * k, -0.2] for k in XB}}}}, None, "wrong fraction on ['m1_r_gmom_x4b2'"),
    # the placebo record
    ({"lib3": {"placebo": {**{a: dict(PL) for a in P.values()}, P[4]: {"CM": "other"}}}}, None,
     "lib3: the placebo arms"),
    ({"ref67": {"placebo": {a: {**PL, "CM": "x"} for a in P.values()}}}, None,
     "the placebo's recorded sha256 differ between the worlds' files"),
    # a placebo or there-and-back arm idle in a world
    ({"pooled": {"identical": (P[3],)}}, None,
     "pooled: placebo or there-and-back arms score exactly as their twin on every seed: ['m1_r_gmom_p3']"),
    ({"ref50": {"identical": (XB[4],)}}, None,
     "ref50: placebo or there-and-back arms score exactly as their twin on every seed: ['m1_r_gmom_x4b2']"),
    # a placebo arm left at p2's fraction
    ({"whole": {"idle": (P[4],)}}, None,
     "whole: placebo arms score exactly as m1_r_gmom_p2 on every seed (an idle fraction): ['m1_r_gmom_p4']"),
    ({"lib2": {"idle": (P[3], P[5])}}, None, "(an idle fraction): ['m1_r_gmom_p3', 'm1_r_gmom_p5']"),
    # the placebo record as each arm's build received it
    ({"pooled": {"placebo": {**{a: dict(PL) for a in P.values()}, P[5]: {"differs between seeds": ["a", "b"]}}}},
     None, "pooled: placebo arms ['m1_r_gmom_p5'] were built with different directions on different seeds"),
    ({"lib3": {"placebo": {a: dict(FENG) for a in P.values()}}}, None,
     "lib3: the placebo arms record Feng's direction in groups ['CM', 'Endocardium']"),
    ({"whole": {"placebo": {a: {**PL, "CM": FENG["CM"]} for a in P.values()}}}, None,
     "whole: the placebo arms record Feng's direction in groups ['CM']"),
    ({"ref67": {"placebo": {a: {"CM": "pl-CM"} for a in P.values()}}}, None,
     "ref67: the placebo's groups ['CM'] are not Feng's"),
]


@pytest.mark.parametrize("over,shift,reason", VOIDS)
def test_the_void_cases(cr, tmp_path, over, shift, reason):
    r = _read(cr, tmp_path, over=over, shift=shift)
    assert not r["valid"] and r["outcome"] == "VOID"
    assert any(reason in x for x in r["validity"]["reasons"]), r["validity"]["reasons"]
    assert "statement_3" not in r
    assert "3431" not in _text(tmp_path)


@pytest.mark.parametrize("src,world,arm", [
    ("t44", "whole", X[2]), ("t44", "whole", "b0"), ("t45", "pooled", P[2]), ("t45", "whole", BASE),
    ("t46", "pooled", X[4]), ("t46", "lib3", X[3]), ("t46", "lib2", "b0"), ("t50", "pooled", X[5]),
    ("t50", "lib3", X[2]), ("t50", "lib2", X[5]), ("t54", "ref67", X[5]), ("t54", "ref50", BASE)])
def test_a_rescored_arm_off_by_one_e_minus_six_in_each_committed_folder_voids(cr, tmp_path, src, world, arm):
    r = _read(cr, tmp_path, shift={(src, world, arm): 1e-6})
    assert not r["valid"]
    assert any(x.startswith(f"{world}: ") and f"try {src.upper()}'s committed scores" in x
               for x in r["validity"]["reasons"]), r["validity"]["reasons"]


def test_a_rescored_arm_within_tolerance_and_a_control_identical_on_nine_seeds_stay_valid(cr, tmp_path):
    per = [0.0] * 9 + [0.3]
    r = _read(cr, tmp_path, over={"lib2": {"per_seed": {P[3]: per}}},
              shift={("t50", "lib2", X[4]): 1e-12, ("t45", "whole", P[2]): 1e-10})
    assert r["valid"], r["validity"]["reasons"]


@pytest.mark.parametrize("arm,twin", [(P[3], X[3]), (XB[5], X[5]), (P[4], P[2])])
def test_a_twin_equal_on_nine_seeds_stays_valid_and_on_all_ten_voids(cr, tmp_path, arm, twin):
    """The arm equals its twin exactly on nine seeds and differs on the tenth: valid; on all ten: void."""
    nine = {"scores": {arm: SCORE[twin]}, "per_seed": {arm: [0.0] * 9 + [0.5]}}
    r = _read(cr, tmp_path, over={"whole": nine})
    rows = [row for row in json.loads((tmp_path / "out" / "tasks" / "whole.json").read_text())["rows"]]
    assert sum(row["arms"][arm]["score"] == row["arms"][twin]["score"] for row in rows) == 9
    assert r["valid"], r["validity"]["reasons"]
    r = _read(cr, tmp_path / "b", over={"whole": {"scores": {arm: SCORE[twin]}}})
    assert not r["valid"]
    assert any(x.startswith("whole: placebo") and arm in x for x in r["validity"]["reasons"]), r["validity"]["reasons"]


def test_the_placebo_record_is_what_each_placebo_build_received_per_seed(cr, monkeypatch):
    """seed_rows hashes, per seed, the directions object each placebo arm's build received (spied on
    gmom_combo.m1_r_arms), and perm_read.placebo_record folds them: the placebo's sha256, never Feng's."""
    calls, _ = _spy(cr, monkeypatch)
    monkeypatch.setattr(cr.harness, "_b0_t1", lambda damping: (lambda ctx: None))
    monkeypatch.setattr(cr.harness, "BuildContext", lambda **kw: _Ctx(kw["seed"]))
    monkeypatch.setattr(cr.er, "score_arms", lambda board, truth, ref, s, anchors, preds: {
        a: {"score": 0.0} for a in preds})
    monkeypatch.setattr(cr.ob, "anchor_record", lambda anchors: None)
    monkeypatch.setattr(cr.ls, "prediction_sha256", lambda p: "x")
    monkeypatch.setattr(cr.ls, "sparsity", lambda p: 0.0)
    rng = np.random.default_rng(3)
    dirs = {"CM": rng.normal(size=3)}
    placebo = {"CM": dirs["CM"][::-1].copy()}
    floors = {s: SimpleNamespace(score=50.0, anchors=None) for s in (0, 1, 2)}
    built = {}
    cr.seed_rows(None, None, None, None, None, {}, floors, (0, 1, 2), dirs, placebo, {"V-CM": "CM"}, 0.2, ["a"], ["a"],
                 built_sha=built)
    assert set(built) == set(P.values())
    got = cr.pm.placebo_record(built)
    seen = {json.dumps(cr.pm.direction_sha(c["directions"]), sort_keys=True) for c in calls
            if c["directions"] is not dirs}
    assert all(got[a] == cr.pm.direction_sha(placebo) != cr.pm.direction_sha(dirs) for a in P.values())
    assert seen == {json.dumps(got[P[3]], sort_keys=True)}
    # a build that received different directions on one seed is recorded as differing, which the read voids
    used = {}
    cr.build_seed(_Ctx(1), lambda c: None, dirs, placebo, {"V-CM": "CM"}, 0.2, ["a"], ["a"], used=used)
    assert all(used[a] is placebo for a in P.values())
    built[P[4]].add(json.dumps(cr.pm.direction_sha(dirs), sort_keys=True))
    assert cr.pm.DIFFERS in cr.pm.placebo_record(built)[P[4]]
    assert "built_sha=built_sha" in inspect.getsource(cr.cmd_run)
    assert "pm.placebo_record(built_sha)" in inspect.getsource(cr.cmd_run)


def test_a_placebo_record_that_differs_from_try_t45s_voids(cr, tmp_path):
    r = _read(cr, tmp_path, t45_placebo={**PL, "Endocardium": "other"})
    assert not r["valid"]
    assert any("differ from try T45's committed record (T1_E9.5_E10.5_chamber.json) in ['lib2', 'lib3', 'pooled', "
               "'ref50', 'ref67', 'whole']" in x for x in r["validity"]["reasons"]), r["validity"]["reasons"]


def test_a_missing_committed_file_voids(cr, tmp_path):
    r = _read(cr, tmp_path, skip={("t45", "whole")})
    assert not r["valid"]
    reasons = r["validity"]["reasons"]
    assert any("try T45's committed file T1_E9.5_E10.5.json is missing" in x for x in reasons)
    assert any(x.startswith("whole: ") and "try T45's committed scores" in x for x in reasons)
    r2 = _read(cr, tmp_path / "b", skip={("t54", "ref50")})
    assert any(x.startswith("ref50: ") and "try T54's committed scores" in x for x in r2["validity"]["reasons"])


def test_the_reading_refuses_to_overwrite(cr, tmp_path):
    out, dirs = _write(cr, tmp_path, {k: _task(k) for k in WORLDS})
    args = [dirs[s] for s in ("t44", "t45", "t46", "t50", "t54")]
    cr.read(out, *args)
    assert (out / "t58.json").exists() and (out / "t58.txt").exists()
    with pytest.raises(SystemExit, match="already holds"):
        cr.read(out, *args)


# ---------------------------------------------------------------------------------------------------------------------
# The run's refusals


def test_the_run_needs_its_row(cr, monkeypatch):
    monkeypatch.setattr(cr.n12, "tries_rows", lambda: {"T57"})
    monkeypatch.setattr(cr.gc, "n22_passed", lambda path: True)
    with pytest.raises(SystemExit, match="T58 has no row"):
        cr.main(["run", "--world", "pooled", "--out", "t58"])


def test_the_run_needs_n22(cr, monkeypatch):
    monkeypatch.setattr(cr.n12, "tries_rows", lambda: {"T58"})
    with pytest.raises(SystemExit, match="does not hold a passing N22"):
        cr.main(["run", "--world", "whole", "--out", "t58", "--n22", "n22x"])


@pytest.mark.parametrize("world", WORLDS)
def test_the_run_refuses_an_existing_task_file(cr, monkeypatch, world):
    monkeypatch.setattr(cr.n12, "tries_rows", lambda: {"T58"})
    monkeypatch.setattr(cr.gc, "n22_passed", lambda path: True)
    tasks = cr.RESULTS / "t58" / "tasks"
    tasks.mkdir(parents=True)
    (tasks / f"{world}.json").write_text("{}")
    with pytest.raises(SystemExit, match="already holds"):
        cr.main(["run", "--world", world, "--out", "t58"])


def test_the_run_refuses_an_unknown_world(cr):
    with pytest.raises(SystemExit):
        cr.main(["run", "--world", "ref67b", "--out", "t58"])
