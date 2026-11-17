"""experiments/t1_m1/part_depth_read.py (try T56 of experiments/harness/TRIES.md, brief AJ) on synthetic data: the
worlds' routes (spies), the arms, the directions and the reading on synthetic per-seed score tables. No released file,
statistics file, board file or score is read."""

from __future__ import annotations

import copy
import importlib.util
import inspect
import json
import sys
from types import SimpleNamespace

import numpy as np
import pytest
from scipy import sparse
from test_t1_deep_read import _world
from test_t1_lib_step_read import GOLDEN_ARMS_FD18
from test_t1_step_read import _arms_world

from vec import core, harness, t1_scope
from vec import t1_state_transport as ts
from vec.boards import REPO_ROOT
from vec.evaluate import Stage

SCRIPT = REPO_ROOT / "experiments" / "t1_m1" / "part_depth_read.py"
STAMP = {"commit": "f" * 40, "tree_clean": True, "utc": "x", "inputs": {}}
WORLDS = ("lib3", "lib3_t67", "lib3_t33", "lib3_t33b", "ref67", "ref67b", "ref50")
T50_THIN = {"lib3_t67": (2 / 3, [50, 2]), "lib3_t33": (1 / 3, [50, 1]), "lib3_t33b": (1 / 3, [50, 3])}
T54_THIN = {"ref67": (2 / 3, [54, 1]), "ref67b": (2 / 3, [54, 3]), "ref50": (1 / 2, [54, 2])}
THIN = {**T50_THIN, **T54_THIN}
SHALLOWER, DEEPER = "truth shallower than the reference", "truth deeper than the reference"
BASE, X2 = "m1_r_gmom", "m1_r_gmom_x2"
C1, CN, CN_X2 = "m1_r_gmom_c1", "m1_r_gmom_cn", "m1_r_gmom_cn_x2"
ARMS = ("b0", BASE, X2, C1, CN, CN_X2)
ENTRY_OF = {C1: BASE, CN: BASE, CN_X2: X2}
FRACTIONS = {BASE: 0.2, X2: 0.4, C1: 0.2, CN: 0.2, CN_X2: 0.4}
HOLDS = "the common part's gain at equal movement stays positive at every depth read"
TURNS = "it turns"
D_CN, D_CN_X2 = f"{CN} minus {BASE}", f"{CN_X2} minus {X2}"


@pytest.fixture
def pd_(monkeypatch, tmp_path):
    spec = importlib.util.spec_from_file_location("part_depth_read_module", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    monkeypatch.setattr(mod, "RESULTS", tmp_path / "results")
    monkeypatch.setattr(mod.gv, "stamp", lambda started, files: dict(STAMP))
    return mod


# ---------------------------------------------------------------------------------------------------------------------
# The worlds: each is its try's own, made by its try's own code


def test_the_worlds_and_their_thinnings_are_the_tries_own(pd_):
    assert pd_.WORLDS == WORLDS
    assert {k: (v["p"], list(v["stream"])) for k, v in pd_.THINNING_OF.items()} == THIN
    assert pd_.THINNING_OF["lib3_t67"] is pd_.ob.THINNED["lib3_t67"]
    assert pd_.THINNING_OF["ref67"] is pd_.dr.THINNED["ref67"]
    assert pd_.ARMS == ARMS and pd_.arm_fractions(0.2) == pytest.approx(FRACTIONS)
    assert {a: k for a, (k, _) in pd_.FORM_ARMS.items()} == {C1: "c1", CN: "cn", CN_X2: "cn"}
    assert {a: p for a, (_, p) in pd_.FORM_ARMS.items()} == ENTRY_OF
    src = inspect.getsource(pd_.world_route)
    assert "ob.truth_and_contamination(t, data, genes, world)" in src
    assert "dr.world_reference(world" in src and "ob.truth_and_contamination(t, data, genes, HEART)" in src
    assert pd_.THINNING_OF["lib3_t33b"] is pd_.ob.THINNED["lib3_t33b"]
    assert pd_.THINNING_OF["ref67b"] is pd_.dr.THINNED["ref67b"]
    assert pd_.DRAWS == (("lib3_t33", "lib3_t33b"), ("ref67", "ref67b"))
    assert pd_.SIDE == {**dict.fromkeys(WORLDS[:4], SHALLOWER), **dict.fromkeys(WORLDS[4:], DEEPER)}
    with pytest.raises(ValueError, match="unknown world"):
        pd_.world_route("ref100", None, None, [], {}, "ref", None, lambda: None)


def _truth(ref_stage):
    return Stage(X=ref_stage.X[:50].copy(), celltype=ref_stage.celltype[:50])


def _route(pd_, monkeypatch, world):
    """One seed of the world's route on the synthetic inputs, with spies on the truth's maker (try T50's
    thinned_truth and try T46's truth_and_contamination, behind over_read.truth_and_contamination) and on every place
    the reference is read: the floor's anchors, every arm's build, g_mom's fit and every arm's scoring."""
    ctx, groups, dirs, held, full, ref_stage, genes = _world()
    harness._G_MOM_PREP.clear()
    released, start = core.counts_of(ctx.inputs["ref"]), core.counts_of(ctx.inputs["start"])
    truths = {"thinned": [], "t46": []}
    lib3, thin_truth = _truth(ref_stage), Stage(X=ref_stage.X[50:100].copy(), celltype=ref_stage.celltype[50:100])

    def thinned_spy(t, data, genes_, name):
        truths["thinned"].append(name)
        return thin_truth, 0, {"p": pd_.ob.THINNED[name]["p"]}

    def t46_spy(t, data, genes_, name):
        truths["t46"].append(name)
        return lib3, 0

    monkeypatch.setattr(pd_.ob, "thinned_truth", thinned_spy)
    monkeypatch.setattr(pd_.ob.ls, "truth_and_contamination", t46_spy)
    original = dict(ctx.inputs)
    inputs, stage, truth, _, thinning = pd_.world_route(world, None, None, genes, original, "ref", ref_stage,
                                                        lambda: full)
    forms, _ = pd_.pr.all_directions(dirs, list(inputs["ref"].var_names), held)
    seen = {"anchors": [], "score": [], "truth": [], "build_ref": [], "build_dirs": [], "fit_ref": [],
            "fit_start": []}
    real_build, real_fit = t1_scope.build_submission, ts.fit_state_factors

    def floor_spy(t, board, truth_, ref, seed):
        seen["anchors"].append((truth_, ref))
        return SimpleNamespace(anchors={}, score=50.0)

    def build_spy(ref, *a, **kw):
        seen["build_ref"].append(core.counts_of(ref))
        seen["build_dirs"].append(kw.get("directions"))
        return real_build(ref, *a, **kw)

    def fit_spy(early, late, *a, **kw):
        seen["fit_start"].append(core.counts_of(early))
        seen["fit_ref"].append(core.counts_of(late))
        return real_fit(early, late, *a, **kw)

    def score_spy(board, truth_, ref, seed, anchors, preds):
        seen["score"].append(ref)
        seen["truth"].append(truth_)
        return {k: {"score": 0.0} for k in preds}

    monkeypatch.setattr(pd_.gv, "floor_anchors", floor_spy)
    monkeypatch.setattr(t1_scope, "build_submission", build_spy)
    monkeypatch.setattr(ts, "fit_state_factors", fit_spy)
    monkeypatch.setattr(pd_.er, "score_arms", score_spy)
    floors, failures = pd_.ls.all_anchors(None, ctx.board, truth, stage, [0])
    rows = pd_.seed_rows(ctx.task, ctx.board, ctx.board, truth, stage, inputs, floors, [0], dirs, forms, groups, 0.2,
                         held, held)
    return SimpleNamespace(seen=seen, stage=stage, ref_stage=ref_stage, truth=truth, lib3=lib3, thin_truth=thin_truth,
                           inputs=inputs, original=original, counts=core.counts_of(inputs["ref"]), released=released,
                           start=start, rows=rows, failures=failures, truths=truths, thinning=thinning, dirs=dirs,
                           forms=forms)


def test_in_lib3_t67_the_truth_is_try_t50s_and_the_predictions_are_lib3s(pd_, monkeypatch):
    r = _route(pd_, monkeypatch, "lib3_t67")
    assert r.truths == {"thinned": ["lib3_t67"], "t46": []}  # try T50's thinned_truth under its own name
    assert r.truth is r.thin_truth and r.thinning == {"p": 2 / 3}
    assert r.inputs is r.original and r.stage is r.ref_stage  # the inputs and the reference as loaded
    assert not r.failures and len(r.rows) == 1 and set(r.rows[0]["arms"]) == set(ARMS)
    assert all(t is r.thin_truth for t, _ in r.seen["anchors"]) and r.seen["truth"] == [r.thin_truth]
    assert all(s is r.ref_stage for _, s in r.seen["anchors"]) and r.seen["score"] == [r.ref_stage]
    assert all((c != r.released).nnz == 0 for c in r.seen["build_ref"] + r.seen["fit_ref"])
    assert (r.seen["fit_start"][0] != r.start).nnz == 0


def test_in_lib3_the_truth_is_try_t46s_lib3_and_nothing_is_thinned(pd_, monkeypatch):
    r = _route(pd_, monkeypatch, "lib3")
    assert r.truths == {"thinned": [], "t46": ["lib3"]} and r.truth is r.lib3 and r.thinning is None
    assert r.inputs is r.original and r.stage is r.ref_stage
    assert all((c != r.released).nnz == 0 for c in r.seen["build_ref"] + r.seen["fit_ref"])


def test_in_ref67_every_place_the_route_reads_the_reference_gets_try_t54s_thinned_one(pd_, monkeypatch):
    r = _route(pd_, monkeypatch, "ref67")
    assert r.truths == {"thinned": [], "t46": ["lib3"]} and r.truth is r.lib3  # the truth lib3 as released
    assert r.thinning["p"] == 2 / 3 and r.thinning["stream"] == [54, 1] and r.thinning["world"] == "ref67"
    want = pd_.dr.thin_counts(r.released, "ref67")
    assert (r.counts != want).nnz == 0 and (r.counts != r.released).nnz > 0
    assert r.inputs["start"] is r.original["start"]
    assert r.seen["anchors"] and all(s is r.stage for _, s in r.seen["anchors"]) and r.stage is not r.ref_stage
    assert r.seen["score"] == [r.stage] and r.seen["truth"] == [r.lib3]
    assert len(r.seen["build_ref"]) == 1 + 2 * 5  # b0, then m1_r and m1_r_gmom for each of the five M1 arms
    assert all((c != r.counts).nnz == 0 for c in r.seen["build_ref"])
    assert len(r.seen["fit_ref"]) == 1  # g_mom's fit: once, on the thinned pool and the released start
    assert (r.seen["fit_ref"][0] != r.counts).nnz == 0 and (r.seen["fit_start"][0] != r.start).nnz == 0
    # the forms' directions are the entry world's: built from the gene names and the held list only
    _, _, _, held, _, _, _ = _world()
    released_forms, released_rec = pd_.pr.all_directions(r.dirs, list(r.original["ref"].var_names), held)
    for f in ("c1", "cn"):
        assert {g: pd_.zr.direction_sha256(v) for g, v in r.forms[f].items()} == \
            {g: pd_.zr.direction_sha256(v) for g, v in released_forms[f].items()}
    used = {id(d) for d in r.seen["build_dirs"] if d is not None}
    assert used == {id(r.dirs), id(r.forms["c1"]), id(r.forms["cn"])}


def test_ref50_is_try_t54s_half(pd_, monkeypatch):
    r = _route(pd_, monkeypatch, "ref50")
    assert (r.counts != pd_.dr.thin_counts(r.released, "ref50")).nnz == 0
    assert r.thinning["p"] == 1 / 2 and r.thinning["stream"] == [54, 2]


def _independent_thinning(C, p, stream):
    """The thinning written out with the earlier try's p and stream as literals: every count a Binomial(c, p) draw on
    the generator seeded by the stream, cells in row order and within a cell the columns' order."""
    Cs = sparse.csr_matrix(C, dtype=np.int64, copy=True)
    Cs.sort_indices()
    want = sparse.csr_matrix((np.random.default_rng(stream).binomial(Cs.data, p), Cs.indices, Cs.indptr),
                             shape=Cs.shape)
    want.eliminate_zeros()
    return want


@pytest.mark.parametrize("world", list(T50_THIN))
def test_each_thinned_truth_is_drawn_at_try_t50s_own_p_and_stream(pd_, monkeypatch, world):
    """The real route (over_read.thinned_truth and thin_counts, not a spy) on synthetic lib3 counts: the truth equals
    an independent binomial draw at try T50's literal p and stream; the inputs and the reference are as loaded; the
    tries' tables are left as they were."""
    before = (copy.deepcopy(pd_.ob.THINNED), copy.deepcopy(pd_.dr.THINNED))
    ctx, _, _, _, full, ref_stage, genes = _world()
    ref = ctx.inputs["ref"]
    C, obs, var_names = core.counts_of(ref), ref.obs[["celltype"]].copy(), list(ref.var_names)
    lib3 = pd_.ob.as_scored(C, obs, var_names, genes)
    monkeypatch.setattr(pd_.ob, "lib3_counts", lambda t, data, genes_: (lib3, 7, C, obs, var_names))
    inputs = dict(ctx.inputs)
    got_inputs, stage, truth, contamination, rec = pd_.world_route(world, None, None, genes, inputs, "ref", ref_stage,
                                                                   lambda: full)
    p, stream = T50_THIN[world]
    want = pd_.ob.as_scored(_independent_thinning(C, p, stream), obs, var_names, genes)
    assert np.array_equal(np.asarray(truth.X), np.asarray(want.X)) and not np.array_equal(truth.X, lib3.X)
    assert rec["p"] == p and rec["stream"] == stream and contamination == 7
    assert got_inputs is inputs and stage is ref_stage
    assert (copy.deepcopy(pd_.ob.THINNED), copy.deepcopy(pd_.dr.THINNED)) == before


@pytest.mark.parametrize("world", list(T54_THIN))
def test_each_thinned_reference_is_drawn_at_try_t54s_own_p_and_stream(pd_, monkeypatch, world):
    """The real route (deep_read.world_reference) on synthetic inputs: the pool's counts equal an independent binomial
    draw at try T54's literal p and stream; the truth is lib3 as released; try T50's table holds only its own three
    truths afterwards and try T54's is unchanged."""
    before = (copy.deepcopy(pd_.ob.THINNED), copy.deepcopy(pd_.dr.THINNED))
    ctx, _, _, _, full, ref_stage, genes = _world()
    lib3 = _truth(ref_stage)
    monkeypatch.setattr(pd_.ob.ls, "truth_and_contamination", lambda t, data, genes_, name: (lib3, 0))
    released = core.counts_of(ctx.inputs["ref"])
    inputs, stage, truth, _, rec = pd_.world_route(world, None, None, genes, dict(ctx.inputs), "ref", ref_stage,
                                                   lambda: full)
    p, stream = T54_THIN[world]
    want = _independent_thinning(released, p, stream)
    assert (core.counts_of(inputs["ref"]) != want).nnz == 0 and (want != released).nnz > 0
    assert rec["p"] == p and rec["stream"] == stream and truth is lib3 and stage is not ref_stage
    assert (copy.deepcopy(pd_.ob.THINNED), copy.deepcopy(pd_.dr.THINNED)) == before
    assert set(pd_.ob.THINNED) == set(T50_THIN)


def test_the_directions_depend_on_the_gene_names_and_the_held_list_only(pd_, monkeypatch):
    """world_directions loads the entry's directions as tries T50 and T54 do and as try T55 does, and builds try T55's
    forms with try T55's own function; the same gene names give the same record in every world."""
    rng = np.random.default_rng(3)
    common = rng.normal(0, 1, 30)
    stat = {g: common + rng.normal(0, 0.5, 30) for g in ("CM", "Endocardium", "Endothelium")}
    calls = []
    monkeypatch.setattr(pd_.feng_dirs, "feng_directions",
                        lambda genes: (calls.append("t50"), ({g: v.copy() for g, v in stat.items()}, {"f": 1}))[1])
    monkeypatch.setattr(pd_.feng_dirs, "load_config", lambda: {"cfg": 1})
    monkeypatch.setattr(pd_.pr.share_read, "feng_inputs", lambda: ({"st": 1}, {"file": "s", "sha256": "h"}))
    monkeypatch.setattr(pd_.feng_dirs, "directions_from_stats",
                        lambda st, genes, cfg, bracket: (calls.append(("t55", bracket)),
                                                         ({g: v.copy() for g, v in stat.items()}, {}, None))[1])
    genes = [f"g{i}" for i in range(30)]
    d1, f1, rec1, ext1 = pd_.world_directions(genes, ["g0"])
    d2, f2, rec2, ext2 = pd_.world_directions(list(genes), ["g0"])
    assert calls == ["t50", ("t55", ("E9.5", "E14.5"))] * 2
    assert json.dumps(rec1, sort_keys=True) == json.dumps(rec2, sort_keys=True)
    assert ext1["entry_agrees_with_t55s"] == dict.fromkeys(sorted(stat), True) and ext1["t55_source"]["file"] == "s"
    _, want = pd_.pr.all_directions(stat, genes, ["g0"])
    assert json.dumps(rec1, sort_keys=True) == json.dumps(want, sort_keys=True)
    assert set(inspect.signature(pd_.world_directions).parameters) == {"pool_genes", "held_m1"}
    # the two routes disagree on a group, and one gives a group the other lacks: both recorded as not agreeing
    other = {**{g: v.copy() for g, v in stat.items()}, "Epicardium": common.copy()}
    other["CM"] = other["CM"] + 1e-6
    monkeypatch.setattr(pd_.feng_dirs, "feng_directions", lambda genes: (other, {"f": 1}))
    _, _, _, ext3 = pd_.world_directions(genes, ["g0"])
    assert ext3["entry_agrees_with_t55s"] == {"CM": False, "Endocardium": True, "Endothelium": True,
                                              "Epicardium": False}


def test_the_run_makes_the_world_before_anything_else_and_refuses_where_t55_refuses(pd_):
    src = inspect.getsource(pd_.cmd_run)
    assert src.index("world_route(") < src.index("all_anchors(") < src.index("world_directions(") \
        < src.index("seed_rows(")
    assert "except ValueError as e:\n        sys.exit(str(e))" in src


# ---------------------------------------------------------------------------------------------------------------------
# The arms


def test_the_entry_arms_are_try_t50s_and_the_form_arms_try_t55s_bit_for_bit(pd_):
    ctx, groups, dirs, held = _arms_world()
    forms, _ = pd_.pr.all_directions(dirs, list(ctx.inputs["ref"].var_names), held)
    b0 = harness._b0_t1(harness.B0_DAMPING)
    mine = pd_.build_seed(ctx, b0, dirs, forms, groups, 0.2, held, held)
    assert list(mine) == list(ARMS)
    sha = {a: pd_.ls.prediction_sha256(p) for a, p in mine.items()}
    t50 = pd_.ob.build_seed(ctx, b0, dirs, groups, 0.2, held, held)
    assert {a: sha[a] for a in ("b0", BASE, X2)} == {a: pd_.ls.prediction_sha256(t50[a]) for a in ("b0", BASE, X2)}
    t55 = pd_.pr.build_seed(ctx, b0, dirs, forms, groups, 0.2, held, held)
    assert {a: sha[a] for a in (C1, CN, CN_X2)} == {a: pd_.ls.prediction_sha256(t55[a]) for a in (C1, CN, CN_X2)}
    assert {a: sha[a] for a in GOLDEN_ARMS_FD18 if a in sha} == {a: GOLDEN_ARMS_FD18[a] for a in ("b0", BASE, X2)}
    assert len({sha[a] for a in ARMS}) == len(ARMS)  # every form moves its arm


def test_a_form_equal_to_the_entry_builds_the_entry_arm_bit_for_bit(pd_):
    """The idle case: forms whose directions are the entry's give the entry's arm of their fraction bit for bit."""
    ctx, groups, dirs, held = _arms_world()
    copy = {k: {g: v.copy() for g, v in dirs.items()} for k in pd_.pr.FORMS}
    preds = pd_.build_seed(ctx, lambda c: None, dirs, copy, groups, 0.2, held, held)
    for arm, base in ENTRY_OF.items():
        assert preds[arm].X.tobytes() == preds[base].X.tobytes(), arm


# ---------------------------------------------------------------------------------------------------------------------
# The reading on synthetic per-seed score tables


W = {"de_score": 0.25, "de_direction": 0.25, "mmd_u": 0.30, "variogram": 0.20}
B, BX2 = 53.0, 53.6
GROUPS = ("CM", "Endocardium")


def _arm(score, base=B):
    gain = score - base
    sk = {"de_score": 0.4 + gain / 4 / (100 * W["de_score"]),
          "de_direction": 0.4 + gain / 4 / (100 * W["de_direction"]),
          "mmd_u": 0.3 + gain * 0.3 / (100 * W["mmd_u"]), "variogram": 0.2 + gain * 0.2 / (100 * W["variogram"])}
    return {"score": score, "groups": {"distribution": score / 100},
            "metrics": {"mmd_u": 0.1 - gain / 1000, "variogram": 0.01, "de_score": 0.5 + gain / 100,
                        "de_direction": 0.3},
            "skills": sk}


def _dirs(tag=""):
    rec = {"entry": {g: {"sha256": f"e-{g}{tag}", "a_n": 1.25, "a_p": 1.6} for g in GROUPS}}
    for f in ("c1", "s05", "s15", "cn", "cp"):
        rec[f] = {g: {"sha256": f"{f}-{g}", "slope": 0.6, "size_ratio": 1.0} for g in GROUPS}
    return rec


def _gains(g, i):
    return g[i] if isinstance(g, list | tuple) else g


def _task(world, gains=None, seeds=range(10), commit="c" * 40, clean=True, fractions=None, thinning=None,
          directions=None, identical=(), anchor_failures=None, sha_tag=None, agree=None, truth_rows=3431,
          truth_contamination=13):
    """gains: arm -> gain over its entry arm (a number, or ten per-seed values); default cn and cn_x2 +0.5 on every
    seed, c1 -0.2. identical: form arms that score exactly as their entry arm on every seed. Try T50's worlds record
    lib3's prediction hashes (they change the truth only) unless sha_tag says otherwise; agree: the entry directions'
    agreement record (default every group agrees)."""
    gains = {CN: 0.5, CN_X2: 0.5, C1: -0.2, **(gains or {})}
    tag = sha_tag or ("lib3" if world.startswith("lib3") else world)
    rows = []
    for i, s in enumerate(seeds):
        entry = {BASE: B + 0.01 * s, X2: BX2 + 0.01 * s}
        arms = {"b0": _arm(50.0), BASE: _arm(entry[BASE]), X2: _arm(entry[X2], BX2)}
        for a, p in ENTRY_OF.items():
            arms[a] = (dict(arms[p]) if a in identical
                       else _arm(entry[p] + _gains(gains[a], i), B if p == BASE else BX2))
        rows.append({"seed": s, "floor": 50.0, "arms": arms, "sha256": {a: f"{tag}-{a}-{s}" for a in arms}})
    rec = {"label": "synthetic", "task": "T1/E9.5->E10.5:chamber", "world": world,
           "truth": world if world.startswith("lib3") else "lib3", "truth_rows": truth_rows,
           "truth_contamination": truth_contamination, "seeds": list(seeds), "fraction": 0.2,
           "fractions": fractions or dict(FRACTIONS), "rows": rows,
           "external": {"entry_agrees_with_t55s": dict.fromkeys(GROUPS, True) if agree is None else agree},
           "directions": directions if directions is not None else _dirs(), "thinning": None,
           "stamp": {"commit": commit, "tree_clean": clean, "utc": "x", "inputs": {}}}
    if world in THIN:
        p, stream = THIN[world]
        rec["thinning"] = {"p": p, "stream": list(stream), "depth_ratio": p, "rebuild_largest_difference": 0.0,
                           "any_cell_emptied": False, **(thinning or {})}
    if anchor_failures is not None:
        rec["anchor_failures"] = anchor_failures
    return rec


def _write(tmp_path, worlds, shift=None, libraries=True):
    """worlds: world -> task record; the committed try T50 and T54 files beside them (the entry arms' default scores)
    and try T55's lib3 file (lib3's record's scores) with its lib1 and lib2 files (cn +0.4, cn_x2 +0.45 on every seed
    unless libraries is False); shift: (folder, world, arm) -> added to that committed score."""
    shift = shift or {}
    out = tmp_path / "out"
    (out / "tasks").mkdir(parents=True)
    for k, rec in worlds.items():
        (out / "tasks" / f"{k}.json").write_text(json.dumps(rec))
    t50, t54, t55 = tmp_path / "t50" / "tasks", tmp_path / "t54" / "tasks", tmp_path / "t55" / "tasks"
    for d in (t50, t54, t55):
        d.mkdir(parents=True)
    for k in WORLDS:
        src = _task(k)
        folder, d = ("t50", t50) if k.startswith("lib3") else ("t54", t54)
        rows = [{"seed": r["seed"], "arms": {a: {**r["arms"][a], "score": r["arms"][a]["score"]
                                                 + shift.get((folder, k, a), 0.0)} for a in ("b0", BASE, X2)}}
                for r in src["rows"]]
        (d / f"{k}.json").write_text(json.dumps({"world": k, "rows": rows}))
    src = worlds.get("lib3") or _task("lib3")
    rows = [{"seed": r["seed"], "arms": {a: {**r["arms"][a], "score": r["arms"][a]["score"]
                                             + shift.get(("t55", "lib3", a), 0.0)} for a in ARMS}}
            for r in src["rows"]]
    (t55 / "lib3.json").write_text(json.dumps({"truth": "lib3", "rows": rows, "directions": _dirs()}))
    if libraries:
        for k in ("lib1", "lib2"):
            lib = _task(k, gains={CN: 0.4, CN_X2: 0.45})
            (t55 / f"{k}.json").write_text(json.dumps({"truth": k, "rows": lib["rows"], "directions": _dirs()}))
    return out, t50, t54, t55


def _read(pd_, tmp_path, by_world=None, over=None, shift=None, t55_directions=None, libraries=True):
    """by_world: world -> gains; over: world -> _task keyword overrides (None drops the world)."""
    by_world, over = by_world or {}, over or {}
    worlds = {k: _task(k, gains=by_world.get(k), **over.get(k, {})) for k in WORLDS
              if not (k in over and over[k] is None)}
    out, t50, t54, t55 = _write(tmp_path, worlds, shift, libraries)
    if t55_directions is not None:
        rec = json.loads((t55 / "lib3.json").read_text())
        (t55 / "lib3.json").write_text(json.dumps({**rec, "directions": t55_directions}))
    return pd_.read(out, t50, t54, t55)["reading"]


def _text(tmp_path):
    return (tmp_path / "out" / "t56.txt").read_text()


def _seeds(pos, n=10, up=0.5, down=-0.1):
    return [up if i < pos else down for i in range(n)]


def test_the_statement_stays_positive(pd_, tmp_path):
    r = _read(pd_, tmp_path)
    assert r["valid"], r["validity"]["reasons"]
    assert r["outcome"] == HOLDS and r["statement"]["turned"] == [] and r["statement"]["neither"] == []
    assert set(r["differences"]) == set(WORLDS)
    for k in WORLDS:
        d = r["differences"][k]
        assert set(d) == {D_CN, D_CN_X2}
        for c in d.values():
            assert c["mean"] == pytest.approx(0.5) and c["seeds_ahead"] == 10 and c["label"] == "ahead"
            assert c["range"] == pytest.approx([0.5, 0.5])
    assert f"The statement: {HOLDS}" in _text(tmp_path)


def test_a_positive_mean_below_the_bar_with_eight_seeds_ahead_still_stays_positive(pd_, tmp_path):
    g = _seeds(8, up=0.1, down=-0.05)
    r = _read(pd_, tmp_path, by_world={k: {CN: g, CN_X2: g} for k in WORLDS})
    assert r["outcome"] == HOLDS
    assert all(c["label"] == "level" for k in WORLDS for c in r["differences"][k].values())


def test_a_positive_mean_on_seven_of_ten_seeds_does_not_stay_positive(pd_, tmp_path):
    g = _seeds(7)
    r = _read(pd_, tmp_path, by_world={"ref50": {CN_X2: g}})
    c = r["differences"]["ref50"][D_CN_X2]
    assert c["mean"] > 0 and c["seeds_ahead"] == 7
    assert r["outcome"] == f"not separated (ref50 ({DEEPER}): {D_CN_X2})"
    assert r["statement"]["neither"] == [{"world": "ref50", "difference": D_CN_X2, "side": DEEPER}]


@pytest.mark.parametrize("world,side", [("lib3", SHALLOWER), ("lib3_t67", SHALLOWER), ("ref50", DEEPER)])
def test_behind_in_a_world_with_one_draw_turns(pd_, tmp_path, world, side):
    r = _read(pd_, tmp_path, by_world={world: {CN: -0.5}, "ref67": {CN: _seeds(5)}})
    assert r["differences"][world][D_CN]["label"] == "behind"
    assert r["outcome"] == f"{TURNS} ({world} ({side}): {D_CN})"
    assert r["statement"]["turned_in"] == [world] and r["statement"]["neither_in"] == ["ref67"]
    txt = _text(tmp_path)
    assert f"The statement: {TURNS} ({world} ({side}): {D_CN})" in txt
    assert f"also meeting neither condition: ref67 ({DEEPER}): {D_CN}" in txt


@pytest.mark.parametrize("pair,side", [(("lib3_t33", "lib3_t33b"), SHALLOWER), (("ref67", "ref67b"), DEEPER)])
def test_behind_in_one_draw_only_is_not_a_turn(pd_, tmp_path, pair, side):
    a, b = pair
    r = _read(pd_, tmp_path, by_world={a: {CN_X2: -0.5}})
    assert r["differences"][a][D_CN_X2]["label"] == "behind" and r["differences"][b][D_CN_X2]["label"] == "ahead"
    assert r["statement"]["turned"] == []
    assert r["outcome"] == f"not separated ({a} ({side}): {D_CN_X2})"
    # behind in the other draw only, with the first level: both named, neither a turn
    r = _read(pd_, tmp_path / "b", by_world={a: {CN_X2: _seeds(5)}, b: {CN_X2: -0.5}})
    assert r["outcome"] == f"not separated ({a} ({side}): {D_CN_X2}; {b} ({side}): {D_CN_X2})"
    assert f"The statement: not separated ({a} ({side}): {D_CN_X2}; {b} ({side}): {D_CN_X2})" in _text(tmp_path / "b")


@pytest.mark.parametrize("pair,side", [(("lib3_t33", "lib3_t33b"), SHALLOWER), (("ref67", "ref67b"), DEEPER)])
def test_behind_in_both_draws_is_a_turn_and_names_both_and_their_side(pd_, tmp_path, pair, side):
    a, b = pair
    r = _read(pd_, tmp_path, by_world={a: {CN: -0.5}, b: {CN: -0.4}})
    assert r["outcome"] == f"{TURNS} ({a} ({side}): {D_CN}; {b} ({side}): {D_CN})"
    assert r["statement"]["turned_in"] == [a, b] and r["statement"]["neither"] == []
    assert f"The statement: {TURNS} ({a} ({side}): {D_CN}; {b} ({side}): {D_CN})" in _text(tmp_path)


def test_behind_in_both_draws_on_different_differences_is_not_a_turn(pd_, tmp_path):
    r = _read(pd_, tmp_path, by_world={"ref67": {CN: -0.5}, "ref67b": {CN_X2: -0.5}})
    assert r["statement"]["turned"] == []
    assert r["outcome"] == f"not separated (ref67 ({DEEPER}): {D_CN}; ref67b ({DEEPER}): {D_CN_X2})"


def test_a_negative_mean_with_two_seeds_ahead_but_level_by_the_test_is_not_a_turn(pd_, tmp_path):
    """The brief's second turn condition is replaced: only try T46's 'behind' turns."""
    g = _seeds(2, up=0.01, down=-0.1)
    r = _read(pd_, tmp_path, by_world={"lib3_t67": {CN: g}})
    c = r["differences"]["lib3_t67"][D_CN]
    assert c["mean"] < 0 and c["seeds_ahead"] == 2 and c["label"] == "level"
    assert r["outcome"] == f"not separated (lib3_t67 ({SHALLOWER}): {D_CN})"


def test_not_separated_names_every_world_and_difference_with_its_side(pd_, tmp_path):
    r = _read(pd_, tmp_path, by_world={"lib3": {CN: _seeds(3, up=0.01, down=-0.1)},
                                       "ref67": {CN_X2: _seeds(6)}, "ref50": {CN: _seeds(6), CN_X2: _seeds(4)}})
    assert r["outcome"] == (f"not separated (lib3 ({SHALLOWER}): {D_CN}; ref67 ({DEEPER}): {D_CN_X2}; "
                            f"ref50 ({DEEPER}): {D_CN}; ref50 ({DEEPER}): {D_CN_X2})")
    assert r["statement"]["neither_in"] == ["lib3", "ref67", "ref50"]


def test_the_text_prints_the_means_the_two_draws_and_the_libraries_beside_the_statement(pd_, tmp_path):
    r = _read(pd_, tmp_path, by_world={"ref67": {CN: 0.4}, "ref67b": {CN: 0.6}})
    assert r["means"]["ref67"] == pytest.approx({D_CN: 0.4, D_CN_X2: 0.5})
    assert set(r["draws"]) == {"lib3_t33 and lib3_t33b", "ref67 and ref67b"}
    assert r["draws"]["ref67 and ref67b"][D_CN]["ref67b"]["mean"] == pytest.approx(0.6)
    libs = r["libraries_1_and_2"]
    assert set(libs) == {"lib1", "lib2"}
    assert libs["lib1"][D_CN]["mean"] == pytest.approx(0.4) and libs["lib2"][D_CN_X2]["mean"] == pytest.approx(0.45)
    txt = _text(tmp_path)
    assert (f"each world's two means: lib3 ({SHALLOWER}): {D_CN} +0.50, {D_CN_X2} +0.50;") in txt
    assert f"ref67 ({DEEPER}): {D_CN} +0.40, {D_CN_X2} +0.50;" in txt
    assert (f"ref67 and ref67b, {D_CN}: ref67 +0.40, ahead 10 of 10, per seed +0.40 to +0.40, test: ahead | "
            f"ref67b +0.60, ahead 10 of 10, per seed +0.60 to +0.60, test: ahead") in txt
    assert f"lib1 as released (try T55's committed file): {D_CN}: +0.40 points" in txt
    assert f"lib2 as released (try T55's committed file): {D_CN_X2}: +0.45 points" in txt
    assert txt.index("The statement:") < txt.index("each world's two means") < txt.index("The two draws side by side")


def test_missing_library_files_are_printed_as_not_available_and_void_nothing(pd_, tmp_path):
    r = _read(pd_, tmp_path, libraries=False)
    assert r["valid"] and r["libraries_1_and_2"] == {"lib1": None, "lib2": None}
    assert "lib1 as released (try T55's committed file): not available" in _text(tmp_path)


def test_the_reported_quantities_are_read_by_nothing(pd_, tmp_path):
    r = _read(pd_, tmp_path, by_world={k: {C1: -3.0} for k in WORLDS})
    assert r["outcome"] == HOLDS  # c1 far behind changes nothing
    for k in WORLDS:
        c = r["reported"][k][f"{C1} minus {BASE}"]
        assert c["mean"] == pytest.approx(-3.0) and c["label"] == "behind"
        assert r["mean_scores"][k][CN] == pytest.approx(r["mean_scores"][k][BASE] + 0.5)
        d = r["differences"][k][D_CN]
        assert d["points"]["de"] + d["points"]["mmd_u_variogram"] == pytest.approx(0.5)
        assert d["metrics"]["de_score"] == pytest.approx(0.005)
    txt = _text(tmp_path)
    for k in WORLDS:
        assert f"in {k}: {D_CN}: +0.50 points, ahead 10 of 10, per seed +0.50 to +0.50, test: ahead" in txt
        assert f"world {k}: mean scores" in txt
    assert "de terms +0.25, mmd_u and variogram +0.25 (read by nothing)" in txt
    assert "raw metrics de_score +0.005, de_direction +0, mmd_u -0.0005, variogram +0 (read by nothing)" in txt
    assert "world ref50: realised depth ratio 0.5000 (p 0.5000, stream [54, 2])" in txt
    assert "world lib3_t33b: realised depth ratio 0.3333 (p 0.3333, stream [50, 3])" in txt
    assert "world ref67b: realised depth ratio 0.6667 (p 0.6667, stream [54, 3])" in txt
    assert "3431" not in txt and "truth_rows" not in txt and "13 " not in txt  # no count is printed


def test_the_limits_and_choices_are_in_the_text(pd_, tmp_path):
    _read(pd_, tmp_path)
    txt = _text(tmp_path)
    for x in pd_.LIMITS:
        assert f"Limit: {x}" in txt
    for x in pd_.CHOICES:
        assert f"  - {x}" in txt
    assert "0.3 to 1.8 times as deep as the reference and no other" in txt
    assert "equal squared size of the direction, not equal counts moved" in txt
    assert "not compared across worlds" in txt and "one draw for all seeds" in txt
    assert "the floor, the de terms' reference, B0's shift and g_mom's fit" in txt
    assert "no statement may be used with a figure for it" in txt
    assert "lib3 reproduces try T55's committed scores by rule and adds no evidence" in txt
    assert "the statement rests on the six thinned worlds" in txt
    assert "'Stays positive' says nothing about size" in txt
    assert "with no raw metric changing sign, so it is read beside the raw-metric differences" in txt
    assert "two manipulations, not one axis" in txt


VOIDS = [
    ({"ref50": None}, None, "missing worlds ['ref50']"),
    ({"lib3_t33b": None}, None, "missing worlds ['lib3_t33b']"),
    ({"ref67b": None}, None, "missing worlds ['ref67b']"),
    ({"lib3_t33": {"seeds": range(9)}}, None, "lib3_t33: seeds"),
    (None, {("t50", "lib3_t67", BASE): 1e-6}, "lib3_t67: b0, m1_r_gmom, m1_r_gmom_x2 do not reproduce try T50's"),
    (None, {("t50", "lib3_t33b", "b0"): 1e-6}, "lib3_t33b: b0, m1_r_gmom, m1_r_gmom_x2 do not reproduce try T50's"),
    (None, {("t54", "ref67", X2): 1e-6}, "ref67: b0, m1_r_gmom, m1_r_gmom_x2 do not reproduce try T54's"),
    (None, {("t54", "ref67b", BASE): 1e-6}, "ref67b: b0, m1_r_gmom, m1_r_gmom_x2 do not reproduce try T54's"),
    (None, {("t55", "lib3", CN): 1e-6}, "do not reproduce try T55's committed lib3 scores"),
    ({"ref67": {"fractions": {**FRACTIONS, CN_X2: 0.2}}}, None, "wrong fraction on ['m1_r_gmom_cn_x2']"),
    ({"ref50": {"commit": "d" * 40}}, None, "the worlds ran at 2 commits"),
    ({"lib3_t67": {"clean": False}}, None, "tree not clean on ['lib3_t67']"),
    ({"ref67": {"directions": _dirs("x")}}, None, "differ between the worlds' files"),
    ({"lib3_t33": {"thinning": {"depth_ratio": 1 / 3 + 0.02}}}, None, "lib3_t33: realised depth ratio"),
    ({"lib3_t33b": {"thinning": {"depth_ratio": 0.3}}}, None, "lib3_t33b: realised depth ratio"),
    ({"ref50": {"thinning": {"any_cell_emptied": True}}}, None, "ref50: the thinning left a cell with no count"),
    ({"ref67b": {"thinning": {"any_cell_emptied": True}}}, None, "ref67b: the thinning left a cell with no count"),
    ({"ref67": {"thinning": {"rebuild_largest_difference": 1e-3}}}, None, "ref67: the rebuild check's largest"),
    ({"lib3_t67": {"thinning": {"stream": [50, 9]}}}, None, "lib3_t67: the file records p"),
    ({"lib3_t33b": {"thinning": {"stream": [50, 1]}}}, None, "lib3_t33b: the file records p"),
    ({"ref67b": {"thinning": {"stream": [54, 1]}}}, None, "ref67b: the file records p"),
    ({"ref67": {"identical": (CN,)}}, None, "ref67: form arms score exactly as the entry's arm"),
    ({"lib3": {"identical": (CN_X2,)}}, None, "lib3: form arms score exactly as the entry's arm"),
    ({"lib3_t33": {"anchor_failures": [{"seed": 3}]}}, None, "lib3_t33: the floor's anchors invert"),
    # the amendment's three checks
    ({"ref67b": {"agree": {"CM": True, "Endocardium": False}}}, None,
     "ref67b: the entry's directions as tries T50 and T54 load them and as try T55 loads them do not have the same"),
    ({"lib3": {"agree": {}}}, None, "lib3: the entry's directions as tries T50 and T54 load them"),
    ({"ref50": {"agree": "missing"}}, None, "ref50: the entry's directions as tries T50 and T54 load them"),
    ({"lib3_t67": {"sha_tag": "lib3_t67"}}, None, "lib3_t67: the predictions' sha256 do not equal this run's lib3"),
    ({"lib3_t33b": {"sha_tag": "other"}}, None, "lib3_t33b: the predictions' sha256 do not equal this run's lib3"),
    ({"lib3_t33": {"truth_rows": 3430}}, None, "the truth's cells: lib3_t33 does not equal lib3"),
    ({"lib3_t33b": {"truth_contamination": 12}}, None, "the truth's contamination cells: lib3_t33b does not equal"),
]


@pytest.mark.parametrize("over,shift,reason", VOIDS)
def test_the_void_cases(pd_, tmp_path, over, shift, reason):
    r = _read(pd_, tmp_path, over=over, shift=shift)
    assert not r["valid"] and r["outcome"] == "VOID"
    assert any(reason in x for x in r["validity"]["reasons"]), r["validity"]["reasons"]
    assert "statement" not in r
    txt = _text(tmp_path)
    assert "3431" not in txt and "3430" not in txt  # a count check names the world, never the count


def test_a_sha_check_is_made_in_try_t50s_thinned_worlds_only(pd_, tmp_path):
    """try T54's worlds change the inputs, so their predictions differ from lib3's and void nothing."""
    r = _read(pd_, tmp_path, over={k: {"sha_tag": k} for k in ("ref67", "ref67b", "ref50")})
    assert r["valid"], r["validity"]["reasons"]


def test_directions_that_differ_from_try_t55s_void(pd_, tmp_path):
    r = _read(pd_, tmp_path, t55_directions=_dirs("old"))
    assert not r["valid"]
    assert any("differ from try T55's lib3 file in ['lib3', 'lib3_t33', 'lib3_t33b', 'lib3_t67', 'ref50', 'ref67', "
               "'ref67b']" in x for x in r["validity"]["reasons"])
    a_n = _dirs()
    a_n["entry"]["CM"]["a_n"] = 1.2500001
    r = _read(pd_, tmp_path / "b", t55_directions=a_n)
    assert not r["valid"]


def test_a_rescored_arm_within_tolerance_and_a_form_identical_on_nine_seeds_stay_valid(pd_, tmp_path):
    g = [0.0] + [0.5] * 9
    r = _read(pd_, tmp_path, by_world={"ref67": {CN: g}}, shift={("t50", "lib3", "b0"): 1e-12})
    assert r["valid"], r["validity"]["reasons"]
    assert r["differences"]["ref67"][D_CN]["seeds_ahead"] == 9


def test_the_reading_refuses_to_overwrite(pd_, tmp_path):
    out, t50, t54, t55 = _write(tmp_path, {k: _task(k) for k in WORLDS})
    pd_.read(out, t50, t54, t55)
    assert (out / "t56.json").exists() and (out / "t56.txt").exists()
    with pytest.raises(SystemExit, match="already holds"):
        pd_.read(out, t50, t54, t55)


def test_the_committed_files_are_the_tries_own(pd_):
    for k in ("lib3", "lib3_t67", "lib3_t33", "lib3_t33b"):
        assert pd_.committed_file(k).parts[-3:] == ("t50_m1_over_step", "tasks", f"{k}.json")
    for k in ("ref67", "ref67b", "ref50"):
        assert pd_.committed_file(k).parts[-3:] == ("t54_m1_deep_step", "tasks", f"{k}.json")
    assert pd_.T55_LIB3.parts[-3:] == ("t55_m1_part", "tasks", "lib3.json")
    assert pd_.T55_TASKS.parts[-2:] == ("t55_m1_part", "tasks")


# ---------------------------------------------------------------------------------------------------------------------
# The run's refusals


def test_the_run_needs_its_row(pd_, monkeypatch):
    monkeypatch.setattr(pd_.n12, "tries_rows", lambda: {"T54", "T55"})
    monkeypatch.setattr(pd_.gc, "n22_passed", lambda path: True)
    with pytest.raises(SystemExit, match="T56 has no row"):
        pd_.main(["run", "--world", "ref67", "--out", "t56"])


def test_the_run_needs_n22(pd_, monkeypatch):
    monkeypatch.setattr(pd_.n12, "tries_rows", lambda: {"T56"})
    with pytest.raises(SystemExit, match="does not hold a passing N22"):
        pd_.main(["run", "--world", "lib3", "--out", "t56", "--n22", "n22x"])


@pytest.mark.parametrize("world", WORLDS)
def test_the_run_refuses_an_existing_task_file(pd_, monkeypatch, world):
    monkeypatch.setattr(pd_.n12, "tries_rows", lambda: {"T56"})
    monkeypatch.setattr(pd_.gc, "n22_passed", lambda path: True)
    tasks = pd_.RESULTS / "t56" / "tasks"
    tasks.mkdir(parents=True)
    (tasks / f"{world}.json").write_text("{}")
    with pytest.raises(SystemExit, match="already holds"):
        pd_.main(["run", "--world", world, "--out", "t56"])
