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
X15, X5, X6, X8 = "m1_r_gmom_x15", "m1_r_gmom_x5", "m1_r_gmom_x6", "m1_r_gmom_x8"
ENTRY = {1: BASE, 2: X2, 3: X3, 4: X4}
# the entry's eight arms and fractions, in try T50's build order
PLAIN = {BASE: 0.2, X15: 0.3, X2: 0.4, X3: 0.6, X4: 0.8, X5: 1.0, X6: 1.2, X8: 1.6}
CN = {m: f"{BASE}_cn" + ("" if m == 1 else f"_x{m}") for m in (1, 2, 3, 4)}
DC = {m: f"{BASE}_dc" + ("" if m == 1 else f"_x{m}") for m in (1, 2, 3, 4)}
GS = {m: f"{BASE}_gs" + ("" if m == 1 else f"_x{m}") for m in (1, 2, 3, 4)}
ENTRY_OF = {**{CN[m]: ENTRY[m] for m in ENTRY}, **{DC[m]: ENTRY[m] for m in ENTRY}, **{GS[m]: ENTRY[m] for m in ENTRY}}
ARMS = ("b0", *PLAIN, *CN.values(), *DC.values(), *GS.values())
FRACTIONS = {**PLAIN, **{a: PLAIN[p] for a, p in ENTRY_OF.items()}}
AHEAD_ALL = "ahead of the entry's direction in every cell read"
AHEAD_ALL_DC = "the entry plus a further step along c is ahead of the entry in every cell read"
ABOVE = "above the best plain and group-scaled step read"


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


def test_gs_moves_along_c_as_the_forms_and_scales_the_orthogonal_part(pr_):
    """Must-fix after review 2: the group-scaled control has coefficient a_n on c (as both forms) and orthogonal part
    (a_n/p_g) o_g; it is the entry's direction where a_n equals p_g, bit for bit."""
    entry, genes, held = _synthetic_entry()
    sp = pr_.pr.split(entry, genes, held)
    M, c = sp["M"], sp["c"]
    cM = c[M]
    mult = pr_.pr.multiples(entry, sp)
    forms, rec = pr_.gs_forms(entry, genes, held)
    dc, dc_rec = pr_.dc_forms(entry, genes, held)
    ratios = []
    for g, d in entry.items():
        dM, f = d[M], forms[g]
        p = float(cM @ dM) / float(cM @ cM)
        a_n = mult[g]["a_n"]
        assert np.allclose(f[M], dM * (a_n / p), rtol=0, atol=1e-12)
        coef = float(cM @ f[M]) / float(cM @ cM)
        assert coef == pytest.approx(a_n, rel=1e-12)  # the same movement along c as dc and cn
        assert coef == pytest.approx(float(cM @ dc[g][M]) / float(cM @ cM), rel=1e-12)
        assert np.allclose(f[M] - coef * cM, (a_n / p) * (dM - p * cM), rtol=0, atol=1e-12)
        assert np.array_equal(f[~M], d[~M], equal_nan=True)
        assert rec[g]["ratio"] == a_n / p and rec[g]["p_g"] == p and rec[g]["a_n"] == a_n
        assert rec[g]["sha256"] == pr_.zr.direction_sha256(f)
        ratios.append(rec[g]["ratio"])
        _, p2 = pr_.dc_direction(d, c, M, 1.0)
        idle, ratio, _ = pr_.gs_direction(d, c, M, p2)  # a_n = p_g: the entry's direction bit for bit
        assert ratio == 1.0 and idle.tobytes() == np.asarray(d, np.float64).tobytes()
    assert len(set(ratios)) == len(ratios)  # a different ratio in each group: a per-group step size
    assert pr_.gs_record_off(rec, dc_rec) == []
    with pytest.raises(ValueError, match="p_g is not positive"):
        pr_.gs_direction(-entry["CM"], c, M, 1.0)


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
    directions, forms, rec, dc_rec, gs_rec, ext = pr_.world_directions(genes, held)
    assert calls == ["t50", ("t55", ("E9.5", "E14.5"))]
    t55_forms, t55_rec = pr_.pr.all_directions(entry, genes, held)
    assert json.dumps(rec, sort_keys=True) == json.dumps(t55_rec, sort_keys=True)
    assert set(forms) == {"cn", "dc", "gs"}
    assert all(forms["cn"][g].tobytes() == t55_forms["cn"][g].tobytes() for g in entry)
    want, want_rec = pr_.dc_forms(entry, genes, held)
    assert all(forms["dc"][g].tobytes() == want[g].tobytes() for g in entry) and dc_rec == want_rec
    want_gs, want_gs_rec = pr_.gs_forms(entry, genes, held)
    assert all(forms["gs"][g].tobytes() == want_gs[g].tobytes() for g in entry) and gs_rec == want_gs_rec
    assert pr_.gs_record_off(gs_rec, dc_rec) == []  # gs's ratio is dc's a_n / p_g bit for bit
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
    inputs, stage, truth, _, thinning, truth_name = pr_.world_route(world, None, None, genes, original, "ref",
                                                                    ref_stage, lambda: full)
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
                           start=start, rows=rows, failures=failures, thinning=thinning, truth_name=truth_name)


@pytest.mark.parametrize("world,twin", [("ref67_lib2", "ref67"), ("ref50_lib2", "ref50")])
def test_a_library_2_world_has_try_t54s_thinned_reference_and_try_t46s_library_2(pr_, monkeypatch, world, twin):
    before = (copy.deepcopy(pr_.ob.THINNED), copy.deepcopy(pr_.dr.THINNED))
    r = _route(pr_, monkeypatch, world)
    assert r.asked == ["lib2"] and r.truth is r.truths["lib2"] and r.truth_name == "lib2"  # try T46's library 2
    p, stream = THIN[world]
    want = pr_.dr.thin_counts(r.released, twin)  # try T54's own thinning of the twin's name
    assert (r.counts != want).nnz == 0 and (r.counts != _independent_thinning(r.released, p, stream)).nnz == 0
    assert (r.counts != r.released).nnz > 0
    assert r.thinning["p"] == p and r.thinning["stream"] == stream and r.thinning["world"] == twin
    assert r.inputs["start"] is r.original["start"] and r.stage is not r.ref_stage
    assert r.seen["anchors"] and all(t is r.truth and s is r.stage for t, s in r.seen["anchors"])
    assert r.seen["score"] == [r.stage] and r.seen["truth"] == [r.truth]
    assert len(r.seen["build_ref"]) == 1 + 2 * 20  # b0, then m1_r and m1_r_gmom for each of the twenty M1 arms
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
    assert r.asked == [world] and r.truth is r.truths[world] and r.thinning is None and r.truth_name == world
    assert r.inputs is r.original and r.stage is r.ref_stage
    assert all((c != r.released).nnz == 0 for c in r.seen["build_ref"] + r.seen["fit_ref"])


def test_the_worlds_and_their_routes(pr_):
    assert pr_.WORLDS == WORLDS and pr_.TRUTH_OF == TRUTH
    assert {k: (v["p"], list(v["stream"])) for k, v in pr_.THINNING_OF.items()} == THIN
    assert pr_.THINNING_OF["ref67_lib2"] is pr_.dr.THINNED["ref67"]
    assert pr_.THINNING_OF["ref50_lib2"] is pr_.dr.THINNED["ref50"]
    assert pr_.ARMS == ARMS and pr_.arm_fractions(0.2) == pytest.approx(FRACTIONS)
    assert pr_.ENTRY_MULTIPLIERS == pr_.ob.MULTIPLIERS  # try T50's, x1.5 included (should-fix after review 2)
    assert list(pr_.ENTRY_MULTIPLIERS) == list(pr_.ob.MULTIPLIERS)  # in try T50's build order
    assert pr_.ENTRY_RANGE == {1: BASE, 1.5: X15, 2: X2, 3: X3, 4: X4, 5: X5, 6: X6, 8: X8}
    assert pr_.REPRO_ENTRY == {"t50": ("b0", *PLAIN), "t54": ("b0", BASE, X15, X2, X3, X4, X5, X6)}
    assert set(pr_.REPRO_ENTRY["t54"]) - {"b0", BASE} == set(pr_.dr.MULTIPLIERS)  # try T54's own arms
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
    return {"cn": t55["cn"], "dc": pr_.dc_forms(dirs, genes, held)[0], "gs": pr_.gs_forms(dirs, genes, held)[0]}


def test_the_entry_arms_are_try_t50s_and_cn_try_t55s_bit_for_bit(pr_):
    ctx, groups, dirs, held = _arms_world()
    genes = list(ctx.inputs["ref"].var_names)
    forms = _forms(pr_, dirs, genes, held)
    b0 = harness._b0_t1(harness.B0_DAMPING)
    mine = pr_.build_seed(ctx, b0, dirs, forms, groups, 0.2, held, held)
    assert list(mine) == list(ARMS)
    sha = {a: pr_.ls.prediction_sha256(p) for a, p in mine.items()}
    entry = ("b0", *PLAIN)
    t50 = pr_.ob.build_seed(ctx, b0, dirs, groups, 0.2, held, held)
    assert {a: sha[a] for a in entry} == {a: pr_.ls.prediction_sha256(t50[a]) for a in entry}  # x5, x6, x8 included
    golden = ("b0", BASE, X2, X3, X4)
    assert {a: sha[a] for a in golden} == {a: GOLDEN_ARMS_FD18[a] for a in golden}
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
    idle_gs = {}
    for g, d in dirs.items():
        _, p = pr_.dc_direction(d, sp["c"], sp["M"], 1.0)
        idle_gs[g] = pr_.gs_direction(d, sp["c"], sp["M"], p)[0]  # ratio 1: the entry's direction
    forms = {"cn": {g: v.copy() for g, v in dirs.items()}, "dc": idle, "gs": idle_gs}
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
    order = [*PLAIN, *pr_.FORM_ARMS]
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


SCORE = {BASE: 53.0, X15: 53.3, X2: 53.6, X3: 54.0, X4: 54.2, X5: 54.1, X6: 53.9, X8: 53.0}  # turns at x4
RISING = {BASE: 53.0, X15: 53.15, X2: 53.3, X3: 53.6, X4: 53.8, X5: 54.0, X6: 54.1, X8: 54.21}  # still rising at x8
GROUPS = ("CM", "Endocardium")


def _dc_rec(tag="", p_g=0.9, a_n=1.25, size_ratio=None):
    size = 2 - (p_g / a_n) ** 2 if size_ratio is None else size_ratio
    return {g: {"sha256": f"dc-{g}{tag}", "p_g": p_g, "a_n": a_n, "size_ratio": size, "slope": 1.1} for g in GROUPS}


def _gs_rec(dc, tag=""):
    """gs's record matching a dc record: per group the ratio a_n/p_g."""
    return {g: {"sha256": f"gs-{g}{tag}", "ratio": v["a_n"] / v["p_g"], "p_g": v["p_g"], "a_n": v["a_n"]}
            for g, v in dc.items()}


def _gain(g, i):
    return g[i] if isinstance(g, list | tuple) else g


def _task(world, gains=None, seeds=range(10), commit="c" * 40, clean=True, fractions=None, thinning=None,
          directions=None, dc=None, identical=(), anchor_failures=None, agree=None, truth=None, scores=None,
          truth_sha256="default", gs=None):
    """gains: arm -> gain over its entry arm (a number or ten per-seed values); default every cn arm +0.5, every dc
    arm +0.6, the gs arms -0.2 at x1 to x3 and -0.5 at x4 (so that gs is behind the plain steps and its best, x3, is
    inside its range unless a test says otherwise). identical: form arms that score exactly as their entry arm on
    every seed."""
    gains = {**dict.fromkeys(CN.values(), 0.5), **dict.fromkeys(DC.values(), 0.6), **dict.fromkeys(GS.values(), -0.2),
             GS[4]: -0.5, **(gains or {})}
    plain = scores or SCORE
    rows = []
    for i, s in enumerate(seeds):
        entry = {a: v + 0.01 * s for a, v in plain.items()}
        arms = {"b0": _arm(50.0), **{a: _arm(entry[a], plain[a]) for a in plain}}
        for a, p in ENTRY_OF.items():
            arms[a] = dict(arms[p]) if a in identical else _arm(entry[p] + _gain(gains[a], i), plain[p])
        rows.append({"seed": s, "floor": 50.0, "arms": arms, "sha256": {a: f"{world}-{a}-{s}" for a in arms}})
    rec = {"label": "synthetic", "task": "T1/E9.5->E10.5:chamber", "world": world, "truth": truth or TRUTH[world],
           "truth_rows": 3431, "truth_contamination": 13, "seeds": list(seeds), "fraction": 0.2,
           "fractions": fractions or dict(FRACTIONS), "rows": rows,
           "external": {"entry_agrees_with_t55s": dict.fromkeys(GROUPS, True) if agree is None else agree},
           "directions": directions if directions is not None else _dirs(),
           "dc_directions": dc if dc is not None else _dc_rec(),
           "gs_directions": gs if gs is not None else _gs_rec(dc if dc else _dc_rec()), "thinning": None,
           "truth_sha256": f"truth-{truth or TRUTH[world]}" if truth_sha256 == "default" else truth_sha256,
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


COMMITTED = {"t50": {"lib3": ("b0", *PLAIN), "lib2": ("b0", *PLAIN)},
             "t54": {"ref67": ("b0", BASE, X15, X2, X3, X4, X5, X6), "ref50": ("b0", BASE, X15, X2, X3, X4, X5, X6)},
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
    size = "ahead by try T46's test in 24 of 24 cells, cell means"
    assert r["outcome"]["statement_1"] == {"cn": f"{AHEAD_ALL}; {size} +0.50 to +0.50",
                                           "dc": f"{AHEAD_ALL_DC}; {size} +0.60 to +0.60"}
    assert set(r["differences"]["cn"]) == {f"{k} x{m}" for k in WORLDS for m in (1, 2, 3, 4)}
    c = r["differences"]["dc"]["ref50_lib2 x3"]
    assert c["change"] == DC[3] and c["base"] == X3 and c["mean"] == pytest.approx(0.6) and c["label"] == "ahead"
    txt = _text(tmp_path)
    assert f"{AHEAD_ALL}; {size} +0.50 to +0.50" in txt and f"{AHEAD_ALL_DC}; {size} +0.60 to +0.60" in txt
    assert "Statement 1, dc (dc minus the entry's arm of its fraction: exactly the entry plus a further step along c " \
           "of m x fraction x (a_n - p_g) per group, a reading of depth along c, not of the form): " in txt


def test_a_positive_mean_below_the_bar_with_eight_seeds_ahead_is_still_ahead(pr_, tmp_path):
    g = _seeds(8, up=0.1, down=-0.05)
    r = _read(pr_, tmp_path, by_world={k: {CN[m]: g for m in CN} for k in WORLDS})
    assert r["statement_1"]["cn"]["said"] == AHEAD_ALL
    assert all(c["label"] == "level" for c in r["differences"]["cn"].values())
    # the phrase cannot be quoted without its size (fix after review): no cell is ahead by try T46's test
    want = f"{AHEAD_ALL}; ahead by try T46's test in 0 of 24 cells, cell means +0.07 to +0.07"
    assert r["outcome"]["statement_1"]["cn"] == want and r["statement_1"]["cn"]["ahead_by_test"] == []
    assert f"Statement 1, cn (cn minus the entry's arm of its fraction: the deleted orthogonal parts and the raised " \
           f"coefficient on c together): {want}" in _text(tmp_path)


def test_behind_somewhere_names_the_cells(pr_, tmp_path):
    r = _read(pr_, tmp_path, by_world={"ref50_lib2": {DC[3]: -0.5}, "lib2": {DC[4]: -0.4, DC[1]: _seeds(7)}})
    assert r["statement_1"]["dc"]["said"] == "behind somewhere (lib2 x4, ref50_lib2 x3)"
    assert r["outcome"]["statement_1"]["dc"].startswith("behind somewhere (lib2 x4, ref50_lib2 x3); ahead by try "
                                                        "T46's test in 21 of 24 cells, cell means -0.50 to +0.60")
    assert r["statement_1"]["dc"]["neither"] == ["lib2 x1"]  # still recorded
    assert r["statement_1"]["cn"]["said"] == AHEAD_ALL
    txt = _text(tmp_path)
    assert "not of the form): behind somewhere (lib2 x4, ref50_lib2 x3); ahead by try T46's test in 21 of 24" in txt
    assert "also meeting neither condition: lib2 x1" in txt


def test_not_separated_names_every_cell_meeting_neither(pr_, tmp_path):
    r = _read(pr_, tmp_path, by_world={"lib2": {CN[2]: _seeds(7)}, "ref67": {CN[4]: _seeds(2, up=0.01)},
                                       "ref67_lib2": {CN[1]: _seeds(6, up=0.2, down=-0.1)}})
    assert r["differences"]["cn"]["ref67 x4"]["mean"] < 0 and r["differences"]["cn"]["ref67 x4"]["label"] == "level"
    assert r["statement_1"]["cn"]["said"] == "not separated (lib2 x2, ref67 x4, ref67_lib2 x1)"
    assert r["outcome"]["statement_1"]["cn"].startswith("not separated (lib2 x2, ref67 x4, ref67_lib2 x1); ahead by "
                                                        "try T46's test in 21 of 24 cells")
    assert "together): not separated (lib2 x2, ref67 x4, ref67_lib2 x1); ahead by" in _text(tmp_path)


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
    # its quantity goes with it (design attack after review 2): o_g at the forms' coefficient a_n on c
    assert r["outcome"]["statement_2_quantity"] == s2["quantity"] == pr_.QUANTITY_2
    assert "at the forms' coefficient a_n on c (not at the entry's p_g)" in s2["quantity"]
    assert "Statement 2, dc minus cn per cell (dc minus cn: exactly the orthogonal parts o_g at the forms' " \
           "coefficient a_n on c (not at the entry's p_g), so it does not say what o_g is worth at the entry's own " \
           "step): ahead in 2 of 24 cells (lib3 x1, lib2 x3)" in _text(tmp_path)


def _both(plain, scaled):
    """Statement 3's two comparisons as the statement gives them: (arm, mean, seeds ahead, label) each."""
    return "; ".join(f"against the {w} {a}: {m:+.2f} points, ahead {n} of 10, try T46's test: {lab}"
                     for w, (a, m, n, lab) in (("plain best", plain), ("group-scaled best", scaled)))


def test_statement_3_above_and_not_above_the_best_plain_step(pr_, tmp_path):
    # in lib2 cn's best is x1 (53.4), below the entry's best x4 (54.2) and gs's best x3 (53.8): not above; dc's best x4
    # (54.8) is above both
    r = _read(pr_, tmp_path, by_world={"lib2": {CN[1]: 0.4, CN[2]: -2.0, CN[3]: -2.0, CN[4]: -2.0},
                                       "ref50": {DC[4]: -0.5, DC[3]: 0.3}})
    s = r["statement_3"]["lib2"]
    assert s["cn"]["form_multiple"] == 1 and s["cn"]["entry_multiple"] == 4 and s["cn"]["gs_multiple"] == 3
    assert "control" not in s["cn"] and "comparison" not in s["cn"]  # no single control (must-fix after review 3)
    assert s["cn"]["against_plain"]["mean"] == pytest.approx(53.4 - 54.2) and s["cn"]["said"] == "not above"
    assert s["cn"]["against_group_scaled"]["mean"] == pytest.approx(53.4 - 53.8)
    assert s["dc"]["form_multiple"] == 4 and s["dc"]["said"] == ABOVE
    t = r["statement_3"]["ref50"]["dc"]
    assert t["form_multiple"] == 3 and t["against_plain"]["mean"] == pytest.approx(54.3 - 54.2)
    assert t["said"] == ABOVE  # positive, every seed ahead, below the 0.3 bar ...
    # ... and the statement gives both comparisons, each with its size and its label by try T46's test
    assert t["statement"] == f"{ABOVE}; " + _both((X4, 0.1, 10, "level"), (GS[3], 0.5, 10, "ahead"))
    assert r["outcome"]["statement_3"]["lib2"] == {
        "cn": "not above; " + _both((X4, -0.8, 0, "behind"), (GS[3], -0.4, 0, "behind")),
        "dc": f"{ABOVE}; " + _both((X4, 0.6, 10, "ahead"), (GS[3], 1.0, 10, "ahead"))}
    txt = _text(tmp_path)
    assert "    lib2, cn at x1 (plain best x4, group-scaled best x3, c-matched plain multiple x1.39): not above; " \
           + _both((X4, -0.8, 0, "behind"), (GS[3], -0.4, 0, "behind")) in txt


def test_statement_3_above_needs_both_controls_the_reviewers_case(pr_, tmp_path):
    """Must-fix after review 3, the reviewer's case: gs's best has a higher mean than the plain best, and dc holds
    against gs (10 of 10) but not against the plain best (7 of 10). The old code tested the higher-mean control alone
    and said 'above'; both must hold, so it is 'not above', and the statement gives both comparisons."""
    gs3 = [0.7] * 7 + [-0.3] * 3  # gs_x3 - x4: +0.5 on 7 seeds, -0.5 on 3 (x3 is 0.2 below x4)
    dc4 = [1.0] * 7 + [-0.2] * 3  # dc_x4 - x4: +1.0 on 7 seeds, -0.2 on 3
    r = _read(pr_, tmp_path, by_world={"lib3": {GS[3]: gs3, DC[4]: dc4}})
    assert r["valid"], r["validity"]["reasons"]
    means = r["mean_scores"]["lib3"]
    assert means[GS[3]] > means[X4]  # the old code's control was gs
    s = r["statement_3"]["lib3"]["dc"]
    assert s["form_multiple"] == 4 and s["entry_multiple"] == 4 and s["gs_multiple"] == 3
    g, p = s["against_group_scaled"], s["against_plain"]
    assert g["mean"] == pytest.approx(0.44) and g["seeds_ahead"] == 10 and s["holds_group_scaled"]
    assert p["mean"] == pytest.approx(0.64) and p["seeds_ahead"] == 7 and not s["holds_plain"]
    assert s["said"] == "not above"
    assert s["statement"] == "not above; " + _both((X4, 0.64, 7, "level"), (GS[3], 0.44, 10, "ahead"))
    assert r["outcome"]["statement_3"]["lib3"]["dc"] == s["statement"]
    # held against both elsewhere: above
    assert r["statement_3"]["lib2"]["dc"]["said"] == ABOVE


def test_statement_3_a_per_group_step_beats_dc_and_is_not_above(pr_, tmp_path):
    """Must-fix after review 2 (the reviewer's case in score form): a truth along each group's own entry direction,
    deeper in one group than the others, puts dc above every plain step at one multiple, while the group-scaled step
    (the entry's direction times a_n/p_g per group) beats it. Read against the plain steps alone this would be 'above';
    against both it is 'not above'."""
    gains = {**{GS[m]: g for m, g in zip((1, 2, 3, 4), (0.3, 0.6, 1.0, 0.7), strict=True)},
             **dict.fromkeys(CN.values(), -1.0)}
    r = _read(pr_, tmp_path, by_world={"lib3": gains})
    assert r["valid"], r["validity"]["reasons"]
    s = r["statement_3"]["lib3"]["dc"]
    assert s["form_multiple"] == 4 and s["entry_multiple"] == 4 and s["gs_multiple"] == 3
    assert s["against_plain"]["mean"] == pytest.approx(0.6) and s["against_plain"]["seeds_ahead"] == 10
    assert s["holds_plain"] and holds_plain(s["against_plain"])  # the plain steps alone would have read 'above'
    assert s["against_group_scaled"]["mean"] == pytest.approx(-0.2) and not s["holds_group_scaled"]
    assert s["said"] == "not above"
    assert s["statement"] == "not above; " + _both((X4, 0.6, 10, "ahead"), (GS[3], -0.2, 0, "level"))
    # where gs is behind the plain steps (the other worlds), dc holds against both and stays above
    assert r["statement_3"]["lib2"]["dc"]["said"] == ABOVE
    # gs is no form: statements 1 and 2 never read it
    assert set(r["differences"]) == {"cn", "dc"} and set(r["statement_1"]) == {"cn", "dc"}
    assert "lib3, dc at x4 (plain best x4, group-scaled best x3, c-matched plain multiple x5.56): not above; " \
           + _both((X4, 0.6, 10, "ahead"), (GS[3], -0.2, 0, "level")) in _text(tmp_path)


def holds_plain(c):
    return c["mean"] > 0 and c["seeds_ahead"] >= 8


def test_statement_3_takes_the_smallest_multiple_on_a_tie(pr_):
    means = {BASE: 1.0, X2: 2.0, X3: 2.0, X4: 1.5}
    assert pr_.best_multiple(means, {1: BASE, 2: X2, 3: X3, 4: X4}.__getitem__) == 2
    means = {BASE: 1.0, X15: 1.2, X2: 2.0, X3: 2.0, X4: 1.5, X5: 3.0, X6: 3.0, X8: 1.0}
    assert pr_.best_multiple(means, pr_.ENTRY_RANGE.__getitem__, tuple(pr_.ENTRY_RANGE)) == 5


def test_statement_3_reads_the_plain_curve_at_x1_5(pr_, tmp_path):
    """Should-fix after review 2: a plain curve that peaks between x1 and x2 is read at x1.5 (try T50's arm), so the
    entry's best is not understated."""
    peak = {**SCORE, X15: 54.75}
    r = _read(pr_, tmp_path, over={"lib2": {"scores": peak}})
    assert r["valid"], r["validity"]["reasons"]
    s = r["statement_3"]["lib2"]
    assert s["dc"]["entry_multiple"] == 1.5 and s["dc"]["against_plain"]["base"] == X15
    assert s["dc"]["against_plain"]["mean"] == pytest.approx(54.8 - 54.75) and s["dc"]["said"] == ABOVE
    # without x1.5 the plain best would be x4 (54.2) and cn (54.7 at x4) would read 'above'
    assert s["cn"]["said"] == "not above" and s["cn"]["against_plain"]["mean"] == pytest.approx(54.7 - 54.75)
    assert r["mean_scores"]["lib2"][X15] == pytest.approx(54.795)
    assert "lib2, cn at x4 (plain best x1.5, group-scaled best x3" in _text(tmp_path)


def test_statement_3_with_a_rising_entry_curve_is_not_above(pr_, tmp_path):
    """Fix after review: in a world where the plain curve still rises past x4, a form that only moves further along c
    is ahead of the entry in every cell of x1 to x4 and above the entry's best of x1 to x4, but not above the entry's
    best of x1 to x8."""
    gains = {**dict.fromkeys(DC.values(), 0.19), **dict.fromkeys(CN.values(), 0.11)}
    r = _read(pr_, tmp_path, by_world={"lib2": gains}, over={"lib2": {"scores": RISING}})
    assert r["valid"], r["validity"]["reasons"]
    assert all(r["differences"]["dc"][f"lib2 x{m}"]["seeds_ahead"] == 10 for m in (1, 2, 3, 4))
    s = r["statement_3"]["lib2"]
    assert s["dc"]["form_multiple"] == 4 and s["dc"]["entry_multiple"] == 8 and s["dc"]["gs_multiple"] == 3
    assert s["dc"]["against_plain"]["mean"] == pytest.approx(53.99 - 54.21)
    assert s["dc"]["c_matched_multiple"] == pytest.approx(4 * 1.25 / 0.9)  # within the range read: readable
    assert s["dc"]["said"] == "not above" and s["cn"]["said"] == "not above"
    # the brief's range (the entry's best over x1 to x4) would have read "above" here
    four = pr_.best_multiple(r["mean_scores"]["lib2"], pr_.ENTRY_AT.__getitem__)
    assert four == 4 and r["mean_scores"]["lib2"][DC[4]] > r["mean_scores"]["lib2"][X4]
    assert "lib2, dc at x4 (plain best x8, group-scaled best x3, c-matched plain multiple x5.56): not above; " \
           + _both((X8, -0.22, 0, "level"), (GS[3], 0.59, 10, "ahead")) in _text(tmp_path)


def test_statement_3_is_not_readable_at_the_edge_beyond_the_c_matched_multiple(pr_, tmp_path):
    """The entry's best is x8 (the edge) and a form at x4 moves along c as far as the entry at 4 x a_n/p_g = 16.67."""
    dc = _dc_rec(p_g=0.3)
    over = {k: {"dc": dc} for k in WORLDS}
    over["lib2"]["scores"] = RISING
    r = _read(pr_, tmp_path, over=over)
    assert r["valid"], r["validity"]["reasons"]
    s = r["statement_3"]["lib2"]["dc"]
    assert s["entry_multiple"] == 8 and s["form_multiple"] == 4
    want = ("not readable: the entry's best is at the edge of the range read (x8) and the form at x4 moves along c as "
            "far as the entry at x16.67")
    assert s["said"] == want
    assert s["statement"] == f"{want}; " + _both((X8, 0.19, 10, "level"), (GS[3], 1.0, 10, "ahead"))
    assert r["outcome"]["statement_3"]["lib2"]["dc"] == s["statement"] and s["statement"] in _text(tmp_path)
    # where the plain curve turned inside the range, the same ratio stays readable
    assert r["statement_3"]["lib3"]["dc"]["entry_multiple"] == 4
    assert r["statement_3"]["lib3"]["dc"]["said"] == ABOVE


def test_statement_3_is_not_readable_where_the_group_scaled_best_is_at_its_edge(pr_, tmp_path):
    """Should-fix after review 3: gs's best at x4, the edge of its range, could be beaten by a larger group-scaled
    step; statement 3 is not readable there, though the form holds against both controls read."""
    r = _read(pr_, tmp_path, by_world={"ref67": {GS[4]: -0.1}})
    assert r["valid"], r["validity"]["reasons"]
    s = r["statement_3"]["ref67"]
    assert s["dc"]["gs_multiple"] == 4 and s["dc"]["holds_plain"] and s["dc"]["holds_group_scaled"]
    want = "not readable: the group-scaled best is at the edge of its range read (x4)"
    assert s["dc"]["said"] == s["cn"]["said"] == want
    assert s["dc"]["statement"] == f"{want}; " + _both((X4, 0.6, 10, "ahead"), (GS[4], 0.7, 10, "ahead"))
    assert r["outcome"]["statement_3"]["ref67"]["cn"].startswith(want + "; against the plain best")
    assert s["dc"]["statement"] in _text(tmp_path)
    # gs's best inside its range (the default x3): readable
    assert r["statement_3"]["lib3"]["dc"]["gs_multiple"] == 3 and r["statement_3"]["lib3"]["dc"]["said"] == ABOVE
    assert pr_.GS_EDGE == 4


def test_statement_3_names_both_edges(pr_, tmp_path):
    dc = _dc_rec(p_g=0.3)
    over = {k: {"dc": dc} for k in WORLDS}
    over["lib2"]["scores"] = RISING
    r = _read(pr_, tmp_path, over=over, by_world={"lib2": {GS[4]: 0.01}})
    assert r["valid"], r["validity"]["reasons"]
    s = r["statement_3"]["lib2"]["dc"]
    assert s["entry_multiple"] == 8 and s["gs_multiple"] == 4
    assert s["said"] == ("not readable: the entry's best is at the edge of the range read (x8) and the form at x4 "
                         "moves along c as far as the entry at x16.67; the group-scaled best is at the edge of its "
                         "range read (x4)")


def test_c_ratio_is_the_largest_a_n_over_p_g_and_infinite_where_p_g_is_not_positive(pr_):
    rec = {"CM": {"a_n": 0.98, "p_g": 0.67}, "Endocardium": {"a_n": 1.0, "p_g": 0.9}}
    assert pr_.c_ratio(rec) == pytest.approx(0.98 / 0.67)
    assert pr_.c_ratio({**rec, "Blood": {"a_n": 1.0, "p_g": 0.0}}) == float("inf")


def test_the_dc_refusal_is_unreachable_by_cauchy_schwarz(pr_):
    """Should-fix after review: a_n >= |p_g| on every synthetic entry, so dc_forms never refuses unpatched."""
    for seed in range(40):
        rng = np.random.default_rng(seed)
        n = 25
        common = rng.normal(0.5, 1, n)
        entry = {g: common * rng.uniform(0.2, 2) + rng.normal(0, rng.uniform(0.01, 3), n) for g in ("A", "B", "C")}
        genes = [f"g{i}" for i in range(n)]
        try:
            _, rec = pr_.dc_forms(entry, genes, [])
        except ValueError as e:  # only try T55's own refusal (sum(c d_g) not positive) may stop it
            assert "T55 is refused" in str(e)
            continue
        for v in rec.values():
            assert v["a_n"] >= abs(v["p_g"]) and v["size_ratio"] >= 1 - 1e-12
            assert v["size_ratio"] == pytest.approx(2 - (v["p_g"] / v["a_n"]) ** 2, rel=1e-9)
    assert any("cannot fire" in x and "Cauchy-Schwarz" in x for x in pr_.CHOICES)


def test_the_run_records_the_routes_truth_name_and_the_truth_hash(pr_, monkeypatch):
    """Should-fix after review: a route that took the wrong library is caught by the read. The name recorded is the
    one passed to try T46's truth maker, and the truth matrix's sha256 is recorded."""
    monkeypatch.setitem(pr_.TRUTH_OF, "ref67_lib2", "lib3")  # the reviewer's mutant
    r = _route(pr_, monkeypatch, "ref67_lib2")
    assert r.asked == ["lib3"] and r.truth_name == "lib3"
    src = inspect.getsource(pr_.cmd_run)
    assert '"truth": truth_name' in src and '"truth_sha256": ls.prediction_sha256(truth)' in src
    assert "TRUTH_OF[a.world]" not in src


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
                 "by a different ratio a_n/p_g in each group",
                 "can put a form ahead with no direction effect at all",
                 "The c-match is on c alone: the group-scaled step also moves o_g by m x a_n/p_g, further than dc (m) "
                 "and cn (0)",
                 "so the other groups are over-matched there",
                 "it does not rule out a per-group depth at any other ratio",
                 "the form's four against the controls' twelve: eight plain, four group-scaled",
                 "its mean and its label by try T46's test are printed with it and go with it",
                 "the forms' a_n and not the entry's p_g",
                 "dc - cn = o_g exactly", "dc - entry = (a_n - p_g) x c exactly",
                 "a consistent tiny shift passes it in every cell", "is evidence of step size",
                 "no world here varies how much of them the truth carries",
                 "The library-2 thinned worlds have no draw control",
                 "8 of the 24 cells (lib3 x1, lib3 x2, lib2 x1, lib2 x2, ref67 x1, ref67 x2, ref50 x1, ref50 x2) "
                 "repeat for cn the committed scores of tries T55 and T56",
                 "a rescored world (lib3, lib2, ref67 and ref50) is no new evidence for its committed arms",
                 "and every arm of ref67_lib2 and ref50_lib2",
                 "which side that selection favours is not shown",
                 "its 'above' must hold against both",
                 "and where the group-scaled best is at 4, the edge of its own range",
                 "what those parts carry over the two days of the test no released day shows",
                 "no statement may be used with a figure for it"):
        assert want in txt, want
    for x in pr_.STATEMENTS:  # the statements as the code reads them, for the row to quote (should-fix, review 3)
        assert f"  - {x}" in txt
    assert json.loads((tmp_path / "out" / "t57.json").read_text())["statements"] == list(pr_.STATEMENTS)
    assert "where against both the mean is positive with at least 8 of 10 seeds ahead" in pr_.STATEMENTS[2]


def test_the_limit_counts_the_repeated_cells_from_the_reads_gates(pr_, monkeypatch):
    """Must-fix after review 3: the cells repeating committed results are counted from the committed-score gates the
    read applies (8 of 24: lib2's cn and cn_x2 are try T55's), and no rescored world is called new evidence."""
    want = tuple(f"{k} x{m}" for k in WORLDS for m, a in CN.items()
                 if any(a in COMMITTED[f].get(k, ()) for f in ("t55", "t56")))
    assert pr_.repeated_cn_cells() == want and len(want) == 8
    limit = pr_.LIMITS[-1]
    assert limit == pr_.repeated_limit() and limit.startswith(f"8 of the 24 cells ({', '.join(want)}) repeat for cn")
    assert "new evidence is" not in limit and "the library-2 worlds" not in limit
    assert "the entry's m1_r_gmom_x8 in ref67 and ref50" in limit  # try T54 has no x8: not a rescored arm
    # the count follows the gates, not a figure written by hand
    monkeypatch.delitem(pr_.CN_SOURCE, "lib2")
    assert len(pr_.repeated_cn_cells()) == 6 and pr_.repeated_limit().startswith("6 of the 24 cells (lib3 x1, ")
    monkeypatch.setattr(pr_, "REPRO_CN", (CN[1],))
    assert pr_.repeated_limit().startswith("3 of the 24 cells (lib3 x1, ref67 x1, ref50 x1)")


def test_no_hard_coded_ratio_and_the_recorded_ratios_are_printed(pr_, tmp_path):
    """Should-fix after review 3: the per-group ratios a_n/p_g come from the record, not from figures in the text."""
    for x in (*pr_.LIMITS, *pr_.CHOICES, pr_.__doc__, *pr_.STATEMENTS):
        assert "1.46" not in x and "1.16" not in x and "1.20" not in x
    _read(pr_, tmp_path)
    assert "CM: a_n 1.2500, p_g 0.9000, dc's squared size over the entry's 1.4816, gs's ratio a_n/p_g 1.3889" \
        in _text(tmp_path)


VOIDS = [
    ({"lib2": None}, None, "missing worlds ['lib2']"),
    ({"ref67_lib2": None}, None, "missing worlds ['ref67_lib2']"),
    ({"ref50_lib2": {"seeds": range(9)}}, None, "ref50_lib2: seeds"),
    (None, {("t50", "lib3", X3): 1e-6}, "lib3: b0, m1_r_gmom, m1_r_gmom_x15, m1_r_gmom_x2, m1_r_gmom_x3, "
                                        "m1_r_gmom_x4, m1_r_gmom_x5, m1_r_gmom_x6, m1_r_gmom_x8 do not reproduce try "
                                        "T50's"),
    (None, {("t50", "lib2", X4): 1e-6}, "lib2: b0, m1_r_gmom, m1_r_gmom_x15, m1_r_gmom_x2, m1_r_gmom_x3, "
                                        "m1_r_gmom_x4, m1_r_gmom_x5, m1_r_gmom_x6, m1_r_gmom_x8 do not reproduce try "
                                        "T50's"),
    (None, {("t50", "lib2", X8): 1e-6}, "lib2: b0, m1_r_gmom, m1_r_gmom_x15, m1_r_gmom_x2"),
    (None, {("t50", "lib3", X5): 1e-6}, "lib3: b0, m1_r_gmom, m1_r_gmom_x15, m1_r_gmom_x2"),
    # x1.5 reproduced in each committed folder (should-fix after review 2)
    (None, {("t50", "lib3", X15): 1e-6}, "lib3: b0, m1_r_gmom, m1_r_gmom_x15"),
    (None, {("t50", "lib2", X15): 1e-6}, "lib2: b0, m1_r_gmom, m1_r_gmom_x15"),
    (None, {("t54", "ref67", X15): 1e-6}, "ref67: b0, m1_r_gmom, m1_r_gmom_x15"),
    (None, {("t54", "ref50", X15): 1e-6}, "ref50: b0, m1_r_gmom, m1_r_gmom_x15"),
    ({"lib2": {"fractions": {**FRACTIONS, X15: 0.4}}}, None, "wrong fraction on ['m1_r_gmom_x15']"),
    (None, {("t54", "ref67", BASE): 1e-6}, "ref67: b0, m1_r_gmom, m1_r_gmom_x15, m1_r_gmom_x2, m1_r_gmom_x3, "
                                           "m1_r_gmom_x4, m1_r_gmom_x5, m1_r_gmom_x6 do not reproduce try T54's"),
    (None, {("t54", "ref50", X6): 1e-6}, "ref50: b0, m1_r_gmom"),
    ({"ref67": {"fractions": {**FRACTIONS, X8: 0.8}}}, None, "wrong fraction on ['m1_r_gmom_x8']"),
    ({"ref50_lib2": {"dc": _dc_rec(a_n=1.2, size_ratio=1.4)}}, None,
     "ref50_lib2: dc's a_n in CM (1.2) is not the entry's recorded a_n (1.25)"),
    ({"lib3": {"dc": _dc_rec(p_g=1.3, size_ratio=1.0)}}, None, "lib3: dc's p_g in CM (1.3) is above its a_n (1.25)"),
    ({"lib2": {"dc": _dc_rec(size_ratio=0.9)}}, None,
     "lib2: dc's squared size over the entry's in CM (0.9) is not at least 1"),
    ({"ref67": {"truth_sha256": None}}, None, "ref67: the truth matrix's sha256 is not recorded"),
    ({"ref67_lib2": {"truth_sha256": "truth-lib3"}}, None,
     "the truth matrix's sha256 differs between the worlds of lib2 ['lib2', 'ref50_lib2', 'ref67_lib2']"),
    ({"ref50": {"truth_sha256": "truth-other"}}, None,
     "the truth matrix's sha256 differs between the worlds of lib3 ['lib3', 'ref50', 'ref67']"),
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
    # the whole thinning record against the twin's, field for field (should-fix after review 3): each passes try
    # T54's own checks and has the twin's p and stream
    ({"ref67_lib2": {"thinning": {"depth_ratio": 0.6667}}}, None,
     "ref67_lib2: the thinning record is not ref67's field for field (it differs in ['depth_ratio'])"),
    ({"ref50_lib2": {"thinning": {"rebuild_largest_difference": 1e-7}}}, None,
     "ref50_lib2: the thinning record is not ref50's field for field (it differs in ['rebuild_largest_difference'])"),
    ({"ref67_lib2": {"thinning": {"rebuild_pool": "other"}}}, None,
     "ref67_lib2: the thinning record is not ref67's field for field (it differs in ['rebuild_pool'])"),
    ({"ref50": {"thinning": {"depth_ratio": 0.48}}}, None, "ref50: realised depth ratio"),
    ({"ref67_lib2": {"thinning": {"any_cell_emptied": True}}}, None, "ref67_lib2: the thinning left a cell"),
    ({"ref50_lib2": {"thinning": {"rebuild_largest_difference": 1e-3}}}, None, "ref50_lib2: the rebuild check's"),
    ({"lib2": {"thinning": {"p": 0.5}}}, None, "lib2: the file records a thinning"),
    ({"lib2": {"truth": "lib3"}}, None, "lib2: the file holds truth 'lib3', not lib2"),
    ({"ref67_lib2": {"directions": _dirs("x")}}, None, "differ between the worlds' files"),
    ({"ref50": {"dc": _dc_rec("x")}}, None, "dc's recorded sha256, p_g or a_n differ between the worlds' files"),
    ({"lib3": {"dc": {}}}, None, "lib3: no full record of the form dc"),
    ({"lib2": {"identical": (CN[3],)}}, None, "lib2: form (or gs) arms score exactly as the entry's arm"),
    ({"ref50_lib2": {"identical": (DC[1],)}}, None, "ref50_lib2: form (or gs) arms score exactly as the entry's arm"),
    # added after review 2: the group-scaled control and the two libraries' truths
    ({"ref67": {"identical": (GS[2],)}}, None,
     "ref67: form (or gs) arms score exactly as the entry's arm of their fraction on every seed: ['m1_r_gmom_gs_x2']"),
    ({"ref50": {"fractions": {**FRACTIONS, GS[4]: 0.6}}}, None, "wrong fraction on ['m1_r_gmom_gs_x4']"),
    ({"lib3": {"gs": {}}}, None, "lib3: no full record of the group-scaled control gs"),
    ({"ref67_lib2": {"gs": _gs_rec(_dc_rec(), "x")}}, None,
     "gs's recorded sha256 or ratio differ between the worlds' files"),
    ({"lib2": {"gs": {g: {**v, "ratio": 1.3} for g, v in _gs_rec(_dc_rec()).items()}}}, None,
     "lib2: gs's ratio in CM (1.3) is not dc's a_n/p_g"),
    ({"ref50": {"gs": {"CM": _gs_rec(_dc_rec())["CM"]}}}, None,
     "ref50: gs's groups ['CM'] are not dc's ['CM', 'Endocardium']"),
    ({k: {"truth_sha256": "truth-lib3"} for k in ("lib2", "ref67_lib2", "ref50_lib2")}, None,
     "a world of lib3 and a world of lib2 record the same truth matrix sha256 ['truth-lib3']"),
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
