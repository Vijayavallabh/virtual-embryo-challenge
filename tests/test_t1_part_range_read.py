"""experiments/t1_m1/part_range_read.py (try T57 of experiments/harness/TRIES.md, brief AL) on synthetic data: the dc
form on synthetic directions, the worlds' routes (spies), the arms and the reading on synthetic per-seed score tables.
No released file, statistics file, board file or score is read."""

from __future__ import annotations

import copy
import importlib.util
import inspect
import json
import sys
from types import SimpleNamespace

import numpy as np
import pytest
from test_t1_deep_read import _world
from test_t1_lib_step_read import GOLDEN_ARMS_FD18
from test_t1_part_depth_read import _arm, _dirs, _independent_thinning
from test_t1_step_read import _arms_world

from vec import core, harness, t1_scope
from vec import t1_state_transport as ts
from vec.boards import REPO_ROOT
from vec.evaluate import Stage

SCRIPT = REPO_ROOT / "experiments" / "t1_m1" / "part_range_read.py"
STAMP = {"commit": "f" * 40, "tree_clean": True, "utc": "x", "inputs": {}}
WORLDS = ("lib3", "lib2", "ref67", "ref50", "ref67_lib2", "ref50_lib2")
THIN = {"ref67": (2 / 3, [54, 1]), "ref50": (1 / 2, [54, 2]), "ref67_lib2": (2 / 3, [54, 1]),
        "ref50_lib2": (1 / 2, [54, 2])}
TRUTH = {"lib3": "lib3", "lib2": "lib2", "ref67": "lib3", "ref50": "lib3", "ref67_lib2": "lib2", "ref50_lib2": "lib2"}
BASE, X2, X3, X4 = "m1_r_gmom", "m1_r_gmom_x2", "m1_r_gmom_x3", "m1_r_gmom_x4"
ENTRY = {1: BASE, 2: X2, 3: X3, 4: X4}
CN = {m: f"{BASE}_cn" + ("" if m == 1 else f"_x{m}") for m in (1, 2, 3, 4)}
DC = {m: f"{BASE}_dc" + ("" if m == 1 else f"_x{m}") for m in (1, 2, 3, 4)}
ENTRY_OF = {**{CN[m]: ENTRY[m] for m in ENTRY}, **{DC[m]: ENTRY[m] for m in ENTRY}}
ARMS = ("b0", BASE, X2, X3, X4, *CN.values(), *DC.values())
FRACTIONS = {BASE: 0.2, X2: 0.4, X3: 0.6, X4: 0.8, **{a: {BASE: 0.2, X2: 0.4, X3: 0.6, X4: 0.8}[p]
                                                      for a, p in ENTRY_OF.items()}}
AHEAD_ALL = "ahead of the entry's direction in every cell read"


@pytest.fixture
def pr_(monkeypatch, tmp_path):
    spec = importlib.util.spec_from_file_location("part_range_read_module", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    monkeypatch.setattr(mod, "RESULTS", tmp_path / "results")
    monkeypatch.setattr(mod.gv, "stamp", lambda started, files: dict(STAMP))
    return mod


# ---------------------------------------------------------------------------------------------------------------------
# The form dc on synthetic directions


def _synthetic_entry(seed=3, n=30):
    rng = np.random.default_rng(seed)
    common = rng.normal(0.3, 1, n)
    entry = {g: common + rng.normal(0, 0.7, n) for g in ("CM", "Endocardium", "Endothelium")}
    for v in entry.values():
        v[5] = np.nan  # not moved in any group
    return entry, [f"g{i}" for i in range(n)], ["g0"]


def test_dc_is_a_n_c_plus_the_entrys_orthogonal_part(pr_):
    entry, genes, held = _synthetic_entry()
    sp = pr_.pr.split(entry, genes, held)
    M, c = sp["M"], sp["c"]
    mult = pr_.pr.multiples(entry, sp)
    forms, rec = pr_.dc_forms(entry, genes, held)
    cM = c[M]
    for g, d in entry.items():
        dM, f = d[M], forms[g]
        p = float(cM @ dM) / float(cM @ cM)
        a_n = mult[g]["a_n"]
        o = dM - p * cM
        assert np.allclose(f[M], a_n * cM + o, rtol=0, atol=1e-12)
        coef = float(cM @ f[M]) / float(cM @ cM)
        assert coef == pytest.approx(a_n, rel=1e-12)  # the coefficient on c is a_n
        assert np.allclose(f[M] - coef * cM, o, rtol=0, atol=1e-12)  # the entry's orthogonal part, unchanged
        assert abs(float(o @ cM)) < 1e-9
        assert np.array_equal(f[~M], d[~M], equal_nan=True)  # off M the entry's own values
        assert a_n > p  # Cauchy-Schwarz: a_n >= p_g always
        assert rec[g]["p_g"] == p and rec[g]["a_n"] == a_n and rec[g]["sha256"] == pr_.zr.direction_sha256(f)
        assert rec[g]["size_ratio"] == pytest.approx(1 + float(o @ o) / float(dM @ dM), rel=1e-12)
    # cn, by try T55's function: a_n c on M, no orthogonal part
    t55, _ = pr_.pr.all_directions(entry, genes, held)
    for g in entry:
        assert np.allclose(t55["cn"][g][M], mult[g]["a_n"] * cM, rtol=0, atol=1e-12)


def test_dc_equals_the_entrys_direction_bit_for_bit_where_a_n_equals_p_g(pr_):
    entry, genes, held = _synthetic_entry()
    sp = pr_.pr.split(entry, genes, held)
    for d in entry.values():
        _, p = pr_.dc_direction(d, sp["c"], sp["M"], 1.0)
        f, p2 = pr_.dc_direction(d, sp["c"], sp["M"], p)
        assert p2 == p and f.tobytes() == np.asarray(d, np.float64).tobytes()


def test_the_run_refuses_where_a_n_is_below_p_g(pr_, monkeypatch):
    entry, genes, held = _synthetic_entry()
    real = pr_.pr.multiples

    def low(entry_, sp):
        m = real(entry_, sp)
        m["CM"] = {**m["CM"], "a_n": 0.01}
        return m

    monkeypatch.setattr(pr_.pr, "multiples", low)
    with pytest.raises(ValueError, match=r"T57 is refused: a_n is below p_g in \['CM'\]"):
        pr_.dc_forms(entry, genes, held)
    src = inspect.getsource(pr_.cmd_run)
    assert "world_directions(" in src and "except ValueError as e:\n        sys.exit(str(e))" in src


def test_world_directions_builds_cn_by_try_t55_and_dc_on_the_same_entry(pr_, monkeypatch):
    entry, genes, held = _synthetic_entry()
    calls = []
    monkeypatch.setattr(pr_.feng_dirs, "feng_directions",
                        lambda g: (calls.append("t50"), ({k: v.copy() for k, v in entry.items()}, {"f": 1}))[1])
    monkeypatch.setattr(pr_.feng_dirs, "load_config", lambda: {"cfg": 1})
    monkeypatch.setattr(pr_.pr.share_read, "feng_inputs", lambda: ({"st": 1}, {"file": "s", "sha256": "h"}))
    monkeypatch.setattr(pr_.feng_dirs, "directions_from_stats",
                        lambda st, g, cfg, bracket: (calls.append(("t55", bracket)),
                                                     ({k: v.copy() for k, v in entry.items()}, {}, None))[1])
    directions, forms, rec, dc_rec, ext = pr_.world_directions(genes, held)
    assert calls == ["t50", ("t55", ("E9.5", "E14.5"))]
    t55_forms, t55_rec = pr_.pr.all_directions(entry, genes, held)
    assert json.dumps(rec, sort_keys=True) == json.dumps(t55_rec, sort_keys=True)
    assert set(forms) == {"cn", "dc"}
    assert all(forms["cn"][g].tobytes() == t55_forms["cn"][g].tobytes() for g in entry)
    want, want_rec = pr_.dc_forms(entry, genes, held)
    assert all(forms["dc"][g].tobytes() == want[g].tobytes() for g in entry) and dc_rec == want_rec
    assert ext["entry_agrees_with_t55s"] == dict.fromkeys(sorted(entry), True)


# ---------------------------------------------------------------------------------------------------------------------
# The worlds' routes (spies)


def _route(pr_, monkeypatch, world):
    """One seed of the world's route on the synthetic inputs, with spies on try T46's truth maker and on every place
    the reference is read: the floor's anchors, every arm's build, g_mom's fit and every arm's scoring."""
    ctx, groups, dirs, held, full, ref_stage, genes = _world()
    harness._G_MOM_PREP.clear()
    released, start = core.counts_of(ctx.inputs["ref"]), core.counts_of(ctx.inputs["start"])
    truths = {"lib3": Stage(X=ref_stage.X[:50].copy(), celltype=ref_stage.celltype[:50]),
              "lib2": Stage(X=ref_stage.X[50:90].copy(), celltype=ref_stage.celltype[50:90])}
    asked = []

    def t46_spy(t, data, genes_, name):
        asked.append(name)
        return truths[name], 0

    monkeypatch.setattr(pr_.ob.ls, "truth_and_contamination", t46_spy)
    original = dict(ctx.inputs)
    inputs, stage, truth, _, thinning = pr_.world_route(world, None, None, genes, original, "ref", ref_stage,
                                                        lambda: full)
    forms = _forms(pr_, dirs, list(inputs["ref"].var_names), held)
    seen = {"anchors": [], "score": [], "truth": [], "build_ref": [], "fit_ref": [], "fit_start": []}
    real_build, real_fit = t1_scope.build_submission, ts.fit_state_factors

    def floor_spy(t, board, truth_, ref, seed):
        seen["anchors"].append((truth_, ref))
        return SimpleNamespace(anchors={}, score=50.0)

    def build_spy(ref, *a, **kw):
        seen["build_ref"].append(core.counts_of(ref))
        return real_build(ref, *a, **kw)

    def fit_spy(early, late, *a, **kw):
        seen["fit_start"].append(core.counts_of(early))
        seen["fit_ref"].append(core.counts_of(late))
        return real_fit(early, late, *a, **kw)

    def score_spy(board, truth_, ref, seed, anchors, preds):
        seen["score"].append(ref)
        seen["truth"].append(truth_)
        return {k: {"score": 0.0} for k in preds}

    monkeypatch.setattr(pr_.gv, "floor_anchors", floor_spy)
    monkeypatch.setattr(t1_scope, "build_submission", build_spy)
    monkeypatch.setattr(ts, "fit_state_factors", fit_spy)
    monkeypatch.setattr(pr_.er, "score_arms", score_spy)
    floors, failures = pr_.ls.all_anchors(None, ctx.board, truth, stage, [0])
    rows = pr_.seed_rows(ctx.task, ctx.board, ctx.board, truth, stage, inputs, floors, [0], dirs, forms, groups, 0.2,
                         held, held)
    return SimpleNamespace(seen=seen, stage=stage, ref_stage=ref_stage, truth=truth, truths=truths, asked=asked,
                           inputs=inputs, original=original, counts=core.counts_of(inputs["ref"]), released=released,
                           start=start, rows=rows, failures=failures, thinning=thinning)


@pytest.mark.parametrize("world,twin", [("ref67_lib2", "ref67"), ("ref50_lib2", "ref50")])
def test_a_library_2_world_has_try_t54s_thinned_reference_and_try_t46s_library_2(pr_, monkeypatch, world, twin):
    before = (copy.deepcopy(pr_.ob.THINNED), copy.deepcopy(pr_.dr.THINNED))
    r = _route(pr_, monkeypatch, world)
    assert r.asked == ["lib2"] and r.truth is r.truths["lib2"]  # try T46's library 2
    p, stream = THIN[world]
    want = pr_.dr.thin_counts(r.released, twin)  # try T54's own thinning of the twin's name
    assert (r.counts != want).nnz == 0 and (r.counts != _independent_thinning(r.released, p, stream)).nnz == 0
    assert (r.counts != r.released).nnz > 0
    assert r.thinning["p"] == p and r.thinning["stream"] == stream and r.thinning["world"] == twin
    assert r.inputs["start"] is r.original["start"] and r.stage is not r.ref_stage
    assert r.seen["anchors"] and all(t is r.truth and s is r.stage for t, s in r.seen["anchors"])
    assert r.seen["score"] == [r.stage] and r.seen["truth"] == [r.truth]
    assert len(r.seen["build_ref"]) == 1 + 2 * 12  # b0, then m1_r and m1_r_gmom for each of the twelve M1 arms
    assert all((c != r.counts).nnz == 0 for c in r.seen["build_ref"])
    assert len(r.seen["fit_ref"]) == 1 and (r.seen["fit_ref"][0] != r.counts).nnz == 0
    assert (r.seen["fit_start"][0] != r.start).nnz == 0
    assert not r.failures and set(r.rows[0]["arms"]) == set(ARMS)
    assert (copy.deepcopy(pr_.ob.THINNED), copy.deepcopy(pr_.dr.THINNED)) == before


def test_the_library_2_world_and_its_twin_have_the_same_thinned_reference(pr_, monkeypatch):
    a = _route(pr_, monkeypatch, "ref67_lib2")
    b = _route(pr_, monkeypatch, "ref67")
    assert b.asked == ["lib3"] and b.truth is b.truths["lib3"]
    assert (a.counts != b.counts).nnz == 0 and a.thinning == b.thinning
    assert np.array_equal(np.asarray(a.stage.X), np.asarray(b.stage.X))


@pytest.mark.parametrize("world", ["lib3", "lib2"])
def test_the_library_worlds_take_try_t46s_truth_and_thin_nothing(pr_, monkeypatch, world):
    r = _route(pr_, monkeypatch, world)
    assert r.asked == [world] and r.truth is r.truths[world] and r.thinning is None
    assert r.inputs is r.original and r.stage is r.ref_stage
    assert all((c != r.released).nnz == 0 for c in r.seen["build_ref"] + r.seen["fit_ref"])


def test_the_worlds_and_their_routes(pr_):
    assert pr_.WORLDS == WORLDS and pr_.TRUTH_OF == TRUTH
    assert {k: (v["p"], list(v["stream"])) for k, v in pr_.THINNING_OF.items()} == THIN
    assert pr_.THINNING_OF["ref67_lib2"] is pr_.dr.THINNED["ref67"]
    assert pr_.THINNING_OF["ref50_lib2"] is pr_.dr.THINNED["ref50"]
    assert pr_.ARMS == ARMS and pr_.arm_fractions(0.2) == pytest.approx(FRACTIONS)
    assert {a: p for a, (_, p) in pr_.FORM_ARMS.items()} == ENTRY_OF
    src = inspect.getsource(pr_.world_route)
    assert "dr.world_reference(REFERENCE_OF[world]" in src and "ob.truth_and_contamination(" in src
    with pytest.raises(ValueError, match="unknown world"):
        pr_.world_route("ref67b", None, None, [], {}, "ref", None, lambda: None)
    src = inspect.getsource(pr_.cmd_run)
    assert src.index("world_route(") < src.index("all_anchors(") < src.index("world_directions(") \
        < src.index("seed_rows(")


# ---------------------------------------------------------------------------------------------------------------------
# The arms


def _forms(pr_, dirs, genes, held):
    t55, _ = pr_.pr.all_directions(dirs, genes, held)
    return {"cn": t55["cn"], "dc": pr_.dc_forms(dirs, genes, held)[0]}


def test_the_entry_arms_are_try_t50s_and_cn_try_t55s_bit_for_bit(pr_):
    ctx, groups, dirs, held = _arms_world()
    genes = list(ctx.inputs["ref"].var_names)
    forms = _forms(pr_, dirs, genes, held)
    b0 = harness._b0_t1(harness.B0_DAMPING)
    mine = pr_.build_seed(ctx, b0, dirs, forms, groups, 0.2, held, held)
    assert list(mine) == list(ARMS)
    sha = {a: pr_.ls.prediction_sha256(p) for a, p in mine.items()}
    entry = ("b0", BASE, X2, X3, X4)
    t50 = pr_.ob.build_seed(ctx, b0, dirs, groups, 0.2, held, held)
    assert {a: sha[a] for a in entry} == {a: pr_.ls.prediction_sha256(t50[a]) for a in entry}
    assert {a: sha[a] for a in entry} == {a: GOLDEN_ARMS_FD18[a] for a in entry}
    t55_forms, _ = pr_.pr.all_directions(dirs, genes, held)
    t55 = pr_.pr.build_seed(ctx, b0, dirs, t55_forms, groups, 0.2, held, held)
    assert {a: sha[a] for a in (CN[1], CN[2])} == {a: pr_.ls.prediction_sha256(t55[a]) for a in (CN[1], CN[2])}
    assert len(set(sha.values())) == len(ARMS)  # every form moves its arm


def test_a_form_equal_to_the_entry_builds_the_entry_arm_bit_for_bit(pr_):
    """The idle case: dc at a_n = p_g in every group is the entry's direction, and the arms built on it (and on cn
    replaced by the entry) are the entry's arms of their fractions bit for bit."""
    ctx, groups, dirs, held = _arms_world()
    genes = list(ctx.inputs["ref"].var_names)
    sp = pr_.pr.split(dirs, genes, held)
    idle = {}
    for g, d in dirs.items():
        _, p = pr_.dc_direction(d, sp["c"], sp["M"], 1.0)
        idle[g] = pr_.dc_direction(d, sp["c"], sp["M"], p)[0]
    forms = {"cn": {g: v.copy() for g, v in dirs.items()}, "dc": idle}
    preds = pr_.build_seed(ctx, lambda c: None, dirs, forms, groups, 0.2, held, held)
    for arm, base in ENTRY_OF.items():
        assert preds[arm].X.tobytes() == preds[base].X.tobytes(), arm


def test_a_form_arm_differs_from_the_entry_arm_of_its_fraction_in_the_directions_only(pr_, monkeypatch):
    ctx, groups, dirs, held = _arms_world()
    genes = list(ctx.inputs["ref"].var_names)
    forms = _forms(pr_, dirs, genes, held)
    calls = []
    real = t1_scope.build_submission

    def spy(ref, *a, **kw):
        calls.append((core.counts_of(ref), a, kw))
        return real(ref, *a, **kw)

    monkeypatch.setattr(t1_scope, "build_submission", spy)
    b0 = harness._b0_t1(harness.B0_DAMPING)
    pr_.build_seed(ctx, b0, dirs, forms, groups, 0.2, held, held)
    m1 = calls[1:]
    order = [BASE, X2, X3, X4, *pr_.FORM_ARMS]
    assert len(m1) == 2 * len(order)
    by_arm = {a: m1[2 * i: 2 * i + 2] for i, a in enumerate(order)}
    for arm, (form, base) in pr_.FORM_ARMS.items():
        for (c1, a1, k1), (c0, a0, k0) in zip(by_arm[arm], by_arm[base], strict=True):
            assert (c1 != c0).nnz == 0 and len(a1) == len(a0) and set(k1) == set(k0)
            assert a1[0] is a0[0] and a1[1:] == a0[1:]  # board, shares, n, seed
            assert k1["directions"] is forms[form] and k0["directions"] is dirs
            assert k1["direction_fraction"] == k0["direction_fraction"]  # the same float, bit for bit
            assert k0["direction_fraction"] == pytest.approx(FRACTIONS[base], abs=1e-15)
            for key in set(k1) - {"directions", "post_cells"}:
                assert k1[key] is k0[key] or k1[key] == k0[key], key
            assert (k1["post_cells"] is None) == (k0["post_cells"] is None)


# ---------------------------------------------------------------------------------------------------------------------
# The reading on synthetic per-seed score tables


SCORE = {BASE: 53.0, X2: 53.6, X3: 54.0, X4: 54.2}
GROUPS = ("CM", "Endocardium")


def _dc_rec(tag=""):
    return {g: {"sha256": f"dc-{g}{tag}", "p_g": 0.9, "a_n": 1.25, "size_ratio": 1.3, "slope": 1.1} for g in GROUPS}


def _gain(g, i):
    return g[i] if isinstance(g, list | tuple) else g


def _task(world, gains=None, seeds=range(10), commit="c" * 40, clean=True, fractions=None, thinning=None,
          directions=None, dc=None, identical=(), anchor_failures=None, agree=None, truth=None):
    """gains: arm -> gain over its entry arm (a number or ten per-seed values); default every cn arm +0.5, every dc
    arm +0.6. identical: form arms that score exactly as their entry arm on every seed."""
    gains = {**dict.fromkeys(CN.values(), 0.5), **dict.fromkeys(DC.values(), 0.6), **(gains or {})}
    rows = []
    for i, s in enumerate(seeds):
        entry = {a: v + 0.01 * s for a, v in SCORE.items()}
        arms = {"b0": _arm(50.0), **{a: _arm(entry[a], SCORE[a]) for a in SCORE}}
        for a, p in ENTRY_OF.items():
            arms[a] = dict(arms[p]) if a in identical else _arm(entry[p] + _gain(gains[a], i), SCORE[p])
        rows.append({"seed": s, "floor": 50.0, "arms": arms, "sha256": {a: f"{world}-{a}-{s}" for a in arms}})
    rec = {"label": "synthetic", "task": "T1/E9.5->E10.5:chamber", "world": world, "truth": truth or TRUTH[world],
           "truth_rows": 3431, "truth_contamination": 13, "seeds": list(seeds), "fraction": 0.2,
           "fractions": fractions or dict(FRACTIONS), "rows": rows,
           "external": {"entry_agrees_with_t55s": dict.fromkeys(GROUPS, True) if agree is None else agree},
           "directions": directions if directions is not None else _dirs(),
           "dc_directions": dc if dc is not None else _dc_rec(), "thinning": None,
           "stamp": {"commit": commit, "tree_clean": clean, "utc": "x", "inputs": {}}}
    if world in THIN:
        p, stream = THIN[world]
        rec["thinning"] = {"p": p, "stream": list(stream), "depth_ratio": p, "rebuild_largest_difference": 0.0,
                           "any_cell_emptied": False, **(thinning or {})}
    elif thinning is not None:
        rec["thinning"] = thinning
    if anchor_failures is not None:
        rec["anchor_failures"] = anchor_failures
    return rec


COMMITTED = {"t50": {"lib3": ("b0", BASE, X2, X3, X4), "lib2": ("b0", BASE, X2, X3, X4)},
             "t54": {"ref67": ("b0", BASE, X2, X3, X4), "ref50": ("b0", BASE, X2, X3, X4)},
             "t56": {"lib3": (CN[1], CN[2]), "ref67": (CN[1], CN[2]), "ref50": (CN[1], CN[2])},
             "t55": {"lib2": (CN[1], CN[2])}}


def _write(tmp_path, worlds, shift=None):
    """worlds: world -> task record; the committed files of tries T50, T54, T55 and T56 beside them (the world's own
    scores, so that they reproduce); shift: (folder, world, arm) -> added to that committed score."""
    shift = shift or {}
    out = tmp_path / "out"
    (out / "tasks").mkdir(parents=True)
    for k, rec in worlds.items():
        (out / "tasks" / f"{k}.json").write_text(json.dumps(rec))
    dirs = {f: tmp_path / f / "tasks" for f in COMMITTED}
    for folder, by in COMMITTED.items():
        dirs[folder].mkdir(parents=True)
        for k, arms in by.items():
            src = worlds.get(k) or _task(k)
            rows = [{"seed": r["seed"], "arms": {a: {**r["arms"][a], "score": r["arms"][a]["score"]
                                                     + shift.get((folder, k, a), 0.0)} for a in arms}}
                    for r in src["rows"]]
            (dirs[folder] / f"{k}.json").write_text(json.dumps({"world": k, "rows": rows}))
    (dirs["t55"] / "lib3.json").write_text(json.dumps({"truth": "lib3", "rows": [], "directions": _dirs()}))
    return out, dirs


def _read(pr_, tmp_path, by_world=None, over=None, shift=None, t55_directions=None):
    by_world, over = by_world or {}, over or {}
    worlds = {k: _task(k, gains=by_world.get(k), **over.get(k, {})) for k in WORLDS
              if not (k in over and over[k] is None)}
    out, dirs = _write(tmp_path, worlds, shift)
    if t55_directions is not None:
        (dirs["t55"] / "lib3.json").write_text(json.dumps({"truth": "lib3", "rows": [],
                                                           "directions": t55_directions}))
    return pr_.read(out, dirs["t50"], dirs["t54"], dirs["t55"], dirs["t56"])["reading"]


def _text(tmp_path):
    return (tmp_path / "out" / "t57.txt").read_text()


def _seeds(pos, n=10, up=0.5, down=-0.1):
    return [up if i < pos else down for i in range(n)]


def test_both_forms_ahead_in_every_cell(pr_, tmp_path):
    r = _read(pr_, tmp_path)
    assert r["valid"], r["validity"]["reasons"]
    assert r["outcome"]["statement_1"] == {"cn": AHEAD_ALL, "dc": AHEAD_ALL}
    assert set(r["differences"]["cn"]) == {f"{k} x{m}" for k in WORLDS for m in (1, 2, 3, 4)}
    c = r["differences"]["dc"]["ref50_lib2 x3"]
    assert c["change"] == DC[3] and c["base"] == X3 and c["mean"] == pytest.approx(0.6) and c["label"] == "ahead"
    txt = _text(tmp_path)
    assert f"Statement 1, cn: {AHEAD_ALL}" in txt and f"Statement 1, dc: {AHEAD_ALL}" in txt


def test_a_positive_mean_below_the_bar_with_eight_seeds_ahead_is_still_ahead(pr_, tmp_path):
    g = _seeds(8, up=0.1, down=-0.05)
    r = _read(pr_, tmp_path, by_world={k: {CN[m]: g for m in CN} for k in WORLDS})
    assert r["outcome"]["statement_1"]["cn"] == AHEAD_ALL
    assert all(c["label"] == "level" for c in r["differences"]["cn"].values())


def test_behind_somewhere_names_the_cells(pr_, tmp_path):
    r = _read(pr_, tmp_path, by_world={"ref50_lib2": {DC[3]: -0.5}, "lib2": {DC[4]: -0.4, DC[1]: _seeds(7)}})
    assert r["outcome"]["statement_1"]["dc"] == "behind somewhere (lib2 x4, ref50_lib2 x3)"
    assert r["statement_1"]["dc"]["neither"] == ["lib2 x1"]  # still recorded
    assert r["outcome"]["statement_1"]["cn"] == AHEAD_ALL
    txt = _text(tmp_path)
    assert "Statement 1, dc: behind somewhere (lib2 x4, ref50_lib2 x3)" in txt
    assert "also meeting neither condition: lib2 x1" in txt


def test_not_separated_names_every_cell_meeting_neither(pr_, tmp_path):
    r = _read(pr_, tmp_path, by_world={"lib2": {CN[2]: _seeds(7)}, "ref67": {CN[4]: _seeds(2, up=0.01)},
                                       "ref67_lib2": {CN[1]: _seeds(6, up=0.2, down=-0.1)}})
    assert r["differences"]["cn"]["ref67 x4"]["mean"] < 0 and r["differences"]["cn"]["ref67 x4"]["label"] == "level"
    assert r["outcome"]["statement_1"]["cn"] == "not separated (lib2 x2, ref67 x4, ref67_lib2 x1)"
    assert "Statement 1, cn: not separated (lib2 x2, ref67 x4, ref67_lib2 x1)" in _text(tmp_path)


def test_statement_2_counts_and_names_the_cells(pr_, tmp_path):
    r = _read(pr_, tmp_path, by_world={"lib3": {DC[1]: 1.0}, "lib2": {DC[3]: 1.2},
                                       "ref67": {DC[4]: 0.01}, "ref50": {DC[2]: _seeds(5, up=1.0, down=0.4)}})
    s2 = r["statement_2"]
    assert s2["counts"] == {"ahead": 2, "behind": 1, "level": 21}
    assert s2["cells"]["ahead"] == ["lib3 x1", "lib2 x3"] and s2["cells"]["behind"] == ["ref67 x4"]
    assert "ref50 x2" in s2["cells"]["level"]
    assert r["dc_minus_cn"]["lib3 x1"]["mean"] == pytest.approx(0.5)
    assert r["outcome"]["statement_2"].startswith("ahead in 2 of 24 cells (lib3 x1, lib2 x3); behind in 1 of 24 cells "
                                                  "(ref67 x4); level in 21 of 24 cells (")
    assert "Statement 2, dc minus cn per cell: ahead in 2 of 24 cells (lib3 x1, lib2 x3)" in _text(tmp_path)


def test_statement_3_above_and_not_above_the_best_plain_step(pr_, tmp_path):
    # in lib2 cn's best is x1 (53.5), below the entry's best x4 (54.2): not above; dc's best x4 (54.8): above
    r = _read(pr_, tmp_path, by_world={"lib2": {CN[2]: -2.0, CN[3]: -2.0, CN[4]: -2.0},
                                       "ref50": {DC[4]: -0.5, DC[3]: 0.3}})
    s = r["statement_3"]["lib2"]
    assert s["cn"]["form_multiple"] == 1 and s["cn"]["entry_multiple"] == 4
    assert s["cn"]["comparison"]["mean"] == pytest.approx(53.5 - 54.2) and s["cn"]["statement"] == "not above"
    assert s["dc"]["form_multiple"] == 4 and s["dc"]["statement"] == "above the best plain step read"
    t = r["statement_3"]["ref50"]["dc"]
    assert t["form_multiple"] == 3 and t["comparison"]["mean"] == pytest.approx(54.3 - 54.2)
    assert t["statement"] == "above the best plain step read"  # positive, every seed ahead, below the 0.3 bar
    assert r["outcome"]["statement_3"]["lib2"] == {"cn": "not above", "dc": "above the best plain step read"}
    txt = _text(tmp_path)
    assert "lib2, cn at x1 against the entry at x4: -0.70 points, ahead 0 of 10: not above" in txt


def test_statement_3_takes_the_smallest_multiple_on_a_tie(pr_):
    means = {BASE: 1.0, X2: 2.0, X3: 2.0, X4: 1.5}
    assert pr_.best_multiple(means, {1: BASE, 2: X2, 3: X3, 4: X4}.__getitem__) == 2


def test_the_reported_quantities_are_read_by_nothing_and_no_count_is_printed(pr_, tmp_path):
    r = _read(pr_, tmp_path)
    for k in WORLDS:
        assert r["mean_scores"][k][CN[3]] == pytest.approx(r["mean_scores"][k][X3] + 0.5)
    d = r["differences"]["cn"]["lib3 x1"]
    assert d["points"]["de"] + d["points"]["mmd_u_variogram"] == pytest.approx(0.5)
    txt = _text(tmp_path)
    assert "cn, lib3 x1: m1_r_gmom_cn minus m1_r_gmom: +0.50 points, ahead 10 of 10" in txt
    assert "de terms +0.25, mmd_u and variogram +0.25 (read by nothing)" in txt
    assert "raw metrics de_score +0.005" in txt and "world ref50_lib2: mean scores" in txt
    assert "world ref67_lib2: the reference of try T54's ref67, realised depth ratio 0.6667 (p 0.6667, stream " \
           "[54, 1])" in txt
    assert "3431" not in txt and "truth_rows" not in txt and "13 " not in txt


def test_the_limits_and_choices_are_in_the_text(pr_, tmp_path):
    _read(pr_, tmp_path)
    txt = _text(tmp_path)
    for x in pr_.LIMITS:
        assert f"Limit: {x}" in txt
    for x in pr_.CHOICES:
        assert f"  - {x}" in txt
    for want in ("One released day, in-sample, the chamber cut, two libraries",
                 "Each difference is taken inside one world, in that world's own points",
                 "one draw for all seeds and also changes the floor, the de terms' reference, B0's shift and "
                 "g_mom's fit",
                 "covers 1 to 4 times the practice step and the six worlds and nothing else",
                 "statement 3 is the reading that a larger plain step cannot give",
                 "what those parts carry over the two days of the test no released day shows",
                 "no statement may be used with a figure for it"):
        assert want in txt, want


VOIDS = [
    ({"lib2": None}, None, "missing worlds ['lib2']"),
    ({"ref67_lib2": None}, None, "missing worlds ['ref67_lib2']"),
    ({"ref50_lib2": {"seeds": range(9)}}, None, "ref50_lib2: seeds"),
    (None, {("t50", "lib3", X3): 1e-6}, "lib3: b0, m1_r_gmom, m1_r_gmom_x2, m1_r_gmom_x3, m1_r_gmom_x4 do not "
                                        "reproduce try T50's"),
    (None, {("t50", "lib2", X4): 1e-6}, "lib2: b0, m1_r_gmom, m1_r_gmom_x2, m1_r_gmom_x3, m1_r_gmom_x4 do not "
                                        "reproduce try T50's"),
    (None, {("t54", "ref67", BASE): 1e-6}, "ref67: b0, m1_r_gmom, m1_r_gmom_x2, m1_r_gmom_x3, m1_r_gmom_x4 do not "
                                           "reproduce try T54's"),
    (None, {("t54", "ref50", "b0"): 1e-6}, "ref50: b0, m1_r_gmom"),
    (None, {("t56", "lib3", CN[1]): 1e-6}, "lib3: m1_r_gmom_cn, m1_r_gmom_cn_x2 do not reproduce try T56's"),
    (None, {("t56", "ref67", CN[2]): 1e-6}, "ref67: m1_r_gmom_cn, m1_r_gmom_cn_x2 do not reproduce try T56's"),
    (None, {("t56", "ref50", CN[1]): 1e-6}, "ref50: m1_r_gmom_cn, m1_r_gmom_cn_x2 do not reproduce try T56's"),
    (None, {("t55", "lib2", CN[2]): 1e-6}, "lib2: m1_r_gmom_cn, m1_r_gmom_cn_x2 do not reproduce try T55's"),
    ({"ref67": {"fractions": {**FRACTIONS, DC[3]: 0.8}}}, None, "wrong fraction on ['m1_r_gmom_dc_x3']"),
    ({"ref50_lib2": {"commit": "d" * 40}}, None, "the worlds ran at 2 commits"),
    ({"lib2": {"clean": False}}, None, "tree not clean on ['lib2']"),
    ({"ref67_lib2": {"thinning": {"stream": [54, 3]}}}, None,
     "ref67_lib2: the thinning record's p and stream (0.6666666666666666, [54, 3]) are not ref67's"),
    ({"ref50_lib2": {"thinning": {"p": 2 / 3, "stream": [54, 1]}}}, None, "are not ref50's"),
    ({"ref50": {"thinning": {"depth_ratio": 0.48}}}, None, "ref50: realised depth ratio"),
    ({"ref67_lib2": {"thinning": {"any_cell_emptied": True}}}, None, "ref67_lib2: the thinning left a cell"),
    ({"ref50_lib2": {"thinning": {"rebuild_largest_difference": 1e-3}}}, None, "ref50_lib2: the rebuild check's"),
    ({"lib2": {"thinning": {"p": 0.5}}}, None, "lib2: the file records a thinning"),
    ({"lib2": {"truth": "lib3"}}, None, "lib2: the file holds truth 'lib3', not lib2"),
    ({"ref67_lib2": {"directions": _dirs("x")}}, None, "differ between the worlds' files"),
    ({"ref50": {"dc": _dc_rec("x")}}, None, "dc's recorded sha256, p_g or a_n differ between the worlds' files"),
    ({"lib3": {"dc": {}}}, None, "lib3: no full record of the form dc"),
    ({"lib2": {"identical": (CN[3],)}}, None, "lib2: form arms score exactly as the entry's arm"),
    ({"ref50_lib2": {"identical": (DC[1],)}}, None, "ref50_lib2: form arms score exactly as the entry's arm"),
    ({"ref67_lib2": {"anchor_failures": [{"seed": 3}]}}, None, "ref67_lib2: the floor's anchors invert"),
    ({"lib2": {"agree": {"CM": True, "Endocardium": False}}}, None,
     "lib2: the entry's directions as tries T50 and T54 load them and as try T55 loads them do not have the same"),
]


@pytest.mark.parametrize("over,shift,reason", VOIDS)
def test_the_void_cases(pr_, tmp_path, over, shift, reason):
    r = _read(pr_, tmp_path, over=over, shift=shift)
    assert not r["valid"] and r["outcome"] == "VOID"
    assert any(reason in x for x in r["validity"]["reasons"]), r["validity"]["reasons"]
    assert "statement_1" not in r
    assert "3431" not in _text(tmp_path)


def test_directions_that_differ_from_try_t55s_void(pr_, tmp_path):
    r = _read(pr_, tmp_path, t55_directions=_dirs("old"))
    assert not r["valid"]
    assert any("differ from try T55's lib3 file in ['lib2', 'lib3', 'ref50', 'ref50_lib2', 'ref67', 'ref67_lib2']"
               in x for x in r["validity"]["reasons"])


def test_a_rescored_arm_within_tolerance_and_a_form_identical_on_nine_seeds_stay_valid(pr_, tmp_path):
    g = [0.0] + [0.5] * 9
    r = _read(pr_, tmp_path, by_world={"ref67_lib2": {CN[4]: g}}, shift={("t50", "lib2", X4): 1e-12})
    assert r["valid"], r["validity"]["reasons"]


def test_the_reading_refuses_to_overwrite(pr_, tmp_path):
    out, dirs = _write(tmp_path, {k: _task(k) for k in WORLDS})
    pr_.read(out, dirs["t50"], dirs["t54"], dirs["t55"], dirs["t56"])
    assert (out / "t57.json").exists() and (out / "t57.txt").exists()
    with pytest.raises(SystemExit, match="already holds"):
        pr_.read(out, dirs["t50"], dirs["t54"], dirs["t55"], dirs["t56"])


def test_the_committed_files_are_the_tries_own(pr_):
    f = pr_.committed_files()
    assert {k: p.parts[-3:] for k, p in f["entry"].items()} == {
        "lib3": ("t50_m1_over_step", "tasks", "lib3.json"), "lib2": ("t50_m1_over_step", "tasks", "lib2.json"),
        "ref67": ("t54_m1_deep_step", "tasks", "ref67.json"), "ref50": ("t54_m1_deep_step", "tasks", "ref50.json")}
    assert {k: p.parts[-3:] for k, p in f["cn"].items()} == {
        "lib3": ("t56_cn_depth", "tasks", "lib3.json"), "ref67": ("t56_cn_depth", "tasks", "ref67.json"),
        "ref50": ("t56_cn_depth", "tasks", "ref50.json"), "lib2": ("t55_m1_part", "tasks", "lib2.json")}
    assert f["t55_lib3"]["lib3"].parts[-3:] == ("t55_m1_part", "tasks", "lib3.json")


# ---------------------------------------------------------------------------------------------------------------------
# The run's refusals


def test_the_run_needs_its_row(pr_, monkeypatch):
    monkeypatch.setattr(pr_.n12, "tries_rows", lambda: {"T55", "T56"})
    monkeypatch.setattr(pr_.gc, "n22_passed", lambda path: True)
    with pytest.raises(SystemExit, match="T57 has no row"):
        pr_.main(["run", "--world", "ref67_lib2", "--out", "t57"])


def test_the_run_needs_n22(pr_, monkeypatch):
    monkeypatch.setattr(pr_.n12, "tries_rows", lambda: {"T57"})
    with pytest.raises(SystemExit, match="does not hold a passing N22"):
        pr_.main(["run", "--world", "lib2", "--out", "t57", "--n22", "n22x"])


@pytest.mark.parametrize("world", WORLDS)
def test_the_run_refuses_an_existing_task_file(pr_, monkeypatch, world):
    monkeypatch.setattr(pr_.n12, "tries_rows", lambda: {"T57"})
    monkeypatch.setattr(pr_.gc, "n22_passed", lambda path: True)
    tasks = pr_.RESULTS / "t57" / "tasks"
    tasks.mkdir(parents=True)
    (tasks / f"{world}.json").write_text("{}")
    with pytest.raises(SystemExit, match="already holds"):
        pr_.main(["run", "--world", world, "--out", "t57"])
