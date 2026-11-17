"""experiments/t1_m1/over_read.py (try T50 of experiments/harness/TRIES.md) on synthetic data: the thinned truths, the
arms and the reading. No released file, statistics file, board file or score is read."""

from __future__ import annotations

import importlib.util
import inspect
import json
import re
import sys

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse
from test_t1_lib_step_read import GOLDEN_ARMS_FD18
from test_t1_step_read import _arms_world

from vec import harness
from vec.boards import REPO_ROOT, Board
from vec.evaluate import Stage
from vec.prep import log1p_cp10k

SCRIPT = REPO_ROOT / "experiments" / "t1_m1" / "over_read.py"
OLD = {"m1_r_gmom_x15": 1.5, "m1_r_gmom_x2": 2.0, "m1_r_gmom_x3": 3.0, "m1_r_gmom_x4": 4.0}
NEW = {"m1_r_gmom_x5": 5.0, "m1_r_gmom_x6": 6.0, "m1_r_gmom_x8": 8.0}
MULT = {**OLD, **NEW}
SCALED = tuple(MULT)
ALL_ARMS = ("b0", "m1_r_gmom", *SCALED)
RESCORED = ("b0", "m1_r_gmom", *OLD)
T46_TRUTHS = ("pooled", "lib1", "lib2", "lib3")
THIN = {"lib3_t33": (1 / 3, [50, 1]), "lib3_t67": (2 / 3, [50, 2]), "lib3_t33b": (1 / 3, [50, 3])}
TRUTHS = (*T46_TRUTHS, *THIN)
STAMP = {"commit": "f" * 40, "tree_clean": True, "utc": "x", "inputs": {}}


@pytest.fixture
def orr(monkeypatch, tmp_path):
    spec = importlib.util.spec_from_file_location("over_read_module", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    monkeypatch.setattr(mod, "RESULTS", tmp_path / "results")
    monkeypatch.setattr(mod.gv, "stamp", lambda started, files: dict(STAMP))
    return mod


def _prints_a_count(txt: str, counts) -> bool:
    body = txt.rsplit("\nCommit ", 1)[0]
    return any(re.search(rf"(?<![\d.]){n}(?![\d.])", body) for n in counts)


# ---------------------------------------------------------------------------------------------------------------------
# The thinned truths


PANEL = [f"g{i}" for i in range(8)]
OTHER_GENES = [f"h{i}" for i in range(12)]
LABELS = ["V-CM", "Endocardium", "Blood", "Foregut", "Neural Tube", "aPHM", "Epicardium"]


def _counts_target(n=300, seed=3, tags=("1", "2", "3")):
    """A small E10.5-like target as the released T1 files hold it: .X the log1p CP10k of integer counts over every
    gene of the file (the panel and other genes), no counts layer; names ending in a library tag."""
    rng = np.random.default_rng(seed)
    genes = [*PANEL, *OTHER_GENES]
    C = rng.poisson(3.0, (n, len(genes)))
    C[:, 0] = 1  # no empty cell, and each cell's smallest count is 1 (as recover_counts needs)
    lib = rng.choice(list(tags), n)
    names = [f"AC{i:04d}-x_{lib[i]}" for i in range(n)]
    obs = pd.DataFrame({"celltype": rng.choice(LABELS, n)}, index=names)
    return ad.AnnData(X=log1p_cp10k(sparse.csr_matrix(C)), obs=obs, var=pd.DataFrame(index=genes)), C


def _chamber(monkeypatch, orr, a):
    monkeypatch.setattr(harness, "read_stage", lambda path, genes: a[:, list(genes)].copy())
    monkeypatch.setattr(orr, "read_stage", lambda path, genes=None: a.copy() if genes is None
                        else a[:, list(genes)].copy())
    return {x.key: x for x in harness.load_tasks()["T1:test"]}["T1/E9.5->E10.5:chamber"]


@pytest.fixture
def world(monkeypatch, orr):
    a, C = _counts_target()
    return _chamber(monkeypatch, orr, a), a, C


def test_lib3_is_t46s_and_its_counts_are_the_files_rows_in_that_order(orr, world, tmp_path):
    t, a, C = world
    lib3, contamination = orr.ls.truth_and_contamination(t, tmp_path, PANEL, "lib3")
    got, c, counts, obs, var = orr.lib3_counts(t, tmp_path, PANEL)
    assert got.X.tobytes() == lib3.X.tobytes() and c == contamination
    _, names, _ = orr.ls.filtered_target(t, tmp_path, PANEL)
    rows = [a.obs_names.get_loc(n) for n in np.asarray(names)[orr.ls.library_rows(names, "lib3")]]
    assert rows == sorted(rows)
    assert np.array_equal(counts.toarray(), C[rows]) and var == list(a.var_names)
    assert list(obs["celltype"]) == list(a.obs["celltype"].iloc[rows])


@pytest.mark.parametrize("name", ["lib3_t33", "lib3_t67", "lib3_t33b"])
def test_thinning_keeps_the_rows_their_order_and_count_and_gives_a_ratio_near_p(orr, world, tmp_path, name):
    t, a, C = world
    lib3, c3 = orr.ls.truth_and_contamination(t, tmp_path, PANEL, "lib3")
    th, c, rec = orr.truth_and_contamination(t, tmp_path, PANEL, name)
    assert th.n == lib3.n and th.X.shape == lib3.X.shape and c == c3
    assert np.array_equal(th.celltype, lib3.celltype)  # the same rows in the same order
    assert th.X.dtype == lib3.X.dtype and not np.array_equal(th.X, lib3.X)
    p, stream = THIN[name]
    assert rec["p"] == p and rec["stream"] == stream and rec["of"] == "lib3"
    assert abs(rec["depth_ratio"] - p) < 0.03
    assert rec["rebuild_largest_difference"] <= orr.REBUILD_TOL and rec["any_cell_emptied"] is False
    _, _, counts, _, _ = orr.lib3_counts(t, tmp_path, PANEL)
    thinned = orr.thin_counts(counts, name)
    assert rec["depth_ratio"] == pytest.approx(thinned.sum() / counts.sum())
    assert (thinned.toarray() <= counts.toarray()).all() and (thinned.toarray() >= 0).all()


def test_thinning_is_binomial_on_its_own_stream_the_same_on_every_call_and_differs_between_names(orr):
    rng = np.random.default_rng(9)
    C = sparse.csr_matrix(rng.poisson(4.0, (40, 30)))
    a, b = orr.thin_counts(C, "lib3_t33"), orr.thin_counts(C, "lib3_t33")
    assert (a != b).nnz == 0  # one thinned truth per name, whatever the build seed: no seed enters
    assert "seed" not in inspect.signature(orr.thin_counts).parameters
    assert "seed" not in inspect.signature(orr.thinned_truth).parameters
    Cs = sparse.csr_matrix(C, dtype=np.int64)
    Cs.sort_indices()
    want = np.random.default_rng([50, 1]).binomial(Cs.data, 1 / 3)
    got = sparse.csr_matrix((want, Cs.indices, Cs.indptr), shape=C.shape)
    assert (a != got).nnz == 0
    c67 = orr.thin_counts(C, "lib3_t67")
    assert (c67 != a).nnz > 0
    want67 = np.random.default_rng([50, 2]).binomial(Cs.data, 2 / 3)
    assert (c67 != sparse.csr_matrix((want67, Cs.indices, Cs.indptr), shape=C.shape)).nnz == 0
    c33b = orr.thin_counts(C, "lib3_t33b")
    assert (c33b != a).nnz > 0  # the second draw at p = 1/3 is a different draw
    want33b = np.random.default_rng([50, 3]).binomial(Cs.data, 1 / 3)
    assert (c33b != sparse.csr_matrix((want33b, Cs.indices, Cs.indptr), shape=C.shape)).nnz == 0
    assert orr.THINNED == {"lib3_t33": {"p": 1 / 3, "stream": [50, 1]}, "lib3_t67": {"p": 2 / 3, "stream": [50, 2]},
                           "lib3_t33b": {"p": 1 / 3, "stream": [50, 3]}}
    assert orr.DRAW_NULL == ("lib3_t33b", "lib3_t33")
    streams = [tuple(v["stream"]) for v in orr.THINNED.values()]
    assert len(set(streams)) == len(streams)


def test_the_matrix_the_scorer_gets_is_made_by_the_repositorys_own_routines(orr, world, tmp_path, monkeypatch):
    t, a, C = world
    calls = {"counts_of": 0, "log1p_cp10k": 0, "from_anndata": 0}
    real_counts, real_log, real_from = orr.core.counts_of, orr.prep.log1p_cp10k, Stage.from_anndata.__func__

    def counts_spy(x, *k, **kw):
        calls["counts_of"] += 1
        return real_counts(x, *k, **kw)

    def log_spy(x):
        calls["log1p_cp10k"] += 1
        return real_log(x)

    def from_spy(cls, x, genes, need_coords=False):
        calls["from_anndata"] += 1
        return real_from(cls, x, genes, need_coords=need_coords)

    monkeypatch.setattr(orr.core, "counts_of", counts_spy)
    monkeypatch.setattr(orr.prep, "log1p_cp10k", log_spy)
    monkeypatch.setattr(Stage, "from_anndata", classmethod(from_spy))
    th, _, _ = orr.thinned_truth(t, tmp_path, PANEL, "lib3_t33")
    assert calls["counts_of"] == 1 and calls["log1p_cp10k"] == 2  # the rebuild check and the thinned truth
    assert calls["from_anndata"] >= 2
    _, _, counts, obs, var = orr.lib3_counts(t, tmp_path, PANEL)
    direct = ad.AnnData(X=real_log(orr.thin_counts(counts, "lib3_t33")), obs=obs, var=pd.DataFrame(index=var))
    assert th.X.tobytes() == real_from(Stage, direct, PANEL).X.tobytes()


def test_the_rebuild_check_refuses_values_the_counts_do_not_give_back(orr, monkeypatch, tmp_path):
    """Normalised over more genes than the file holds (the file's other genes dropped after normalising): the counts
    recovered from .X do not give back the matrix the scorer reads, and T50 does not thin it."""
    a, _ = _counts_target()
    cut = a[:, [*PANEL, *OTHER_GENES[:4]]].copy()
    t = _chamber(monkeypatch, orr, cut)
    with pytest.raises(ValueError, match="do not give back lib3's matrix as the scorer reads it"):
        orr.thinned_truth(t, tmp_path, PANEL, "lib3_t67")


def test_values_that_are_not_counts_are_refused(orr, monkeypatch, tmp_path):
    a, _ = _counts_target()
    a.X = sparse.csr_matrix(np.asarray(a.X.toarray()) * 1.37)  # normalised values with no counts behind them
    t = _chamber(monkeypatch, orr, a)
    with pytest.raises(ValueError):
        orr.thinned_truth(t, tmp_path, PANEL, "lib3_t33")


def test_the_run_refuses_a_thinned_truth_whose_counts_do_not_rebuild(orr, monkeypatch):
    a, _ = _counts_target()
    t = _chamber(monkeypatch, orr, a[:, [*PANEL, *OTHER_GENES[:4]]].copy())
    monkeypatch.setattr(orr.n12, "tries_rows", lambda: {"T50"})
    monkeypatch.setattr(orr.gc, "n22_passed", lambda path: True)
    monkeypatch.setattr(orr, "load_boards", lambda: {t.panel_board: Board(key=t.panel_board, task="T1", genes=PANEL,
                                                                          needs_coords=False, split="test")})
    with pytest.raises(SystemExit, match="T50 is refused on lib3_t33: .*does not thin normalised values"):
        orr.main(["run", "--truth", "lib3_t33", "--out", "t50"])
    assert not (orr.RESULTS / "t50" / "tasks" / "lib3_t33.json").exists()


# ---------------------------------------------------------------------------------------------------------------------
# The arms


def test_the_fractions_are_t46s_and_one_one_point_two_and_one_point_six(orr):
    assert orr.MULTIPLIERS == MULT
    fr = orr.arm_fractions(0.2)
    assert fr["m1_r_gmom"] == 0.2
    assert [fr[a] for a in SCALED] == pytest.approx([0.3, 0.4, 0.6, 0.8, 1.0, 1.2, 1.6])


def test_the_new_arms_differ_from_x4s_build_in_the_fraction_only(orr, monkeypatch):
    calls = {}

    def spy(ctx, directions, groups, fraction, held_m1, held_factors):
        calls[round(fraction, 12)] = {"ctx": ctx, "directions": directions, "groups": groups, "fraction": fraction,
                                      "held_m1": held_m1, "held_factors": held_factors}
        return {"m1_r": ("m1_r", fraction), "m1_r_gmom": ("m1_r_gmom", fraction)}

    monkeypatch.setattr(orr.sr.gc, "m1_r_arms", spy)
    ctx, dirs, groups, held_m1, held_f = object(), {"CM": np.zeros(3)}, {"V-CM": "CM"}, ["a", "b"], ["a"]
    arms = orr.scaled_arms(ctx, dirs, groups, 0.2, held_m1, held_f)
    assert set(arms) == {"m1_r_gmom", *SCALED} and len(calls) == 8
    x4 = calls[0.8]
    for f in (1.0, 1.2, 1.6):
        c = calls[f]
        for key in ("ctx", "directions", "groups", "held_m1", "held_factors"):
            assert c[key] is x4[key], f"a new arm differs from x4's build in {key}"
        assert c["fraction"] / x4["fraction"] == pytest.approx({1.0: 5 / 4, 1.2: 6 / 4, 1.6: 8 / 4}[f])
    assert arms["m1_r_gmom_x8"] == ("m1_r_gmom", 0.2 * 8)


def test_t46s_arms_keep_their_output_and_the_new_arms_build_past_a_fraction_of_one(orr, monkeypatch):
    """b0, m1_r_gmom and T46's four scaled arms are fd18e0f8da's bit for bit (T46's pinned table); the route takes
    fractions 1.0, 1.2 and 1.6 with no refusal and moves the cells; the idle multiplier 1 is the base bit for bit."""
    ctx, groups, dirs, held = _arms_world()
    preds = orr.build_seed(ctx, harness._b0_t1(harness.B0_DAMPING), dirs, groups, 0.2, held, held)
    assert list(preds) == ["b0", "m1_r_gmom", *SCALED]
    assert {a: orr.ls.prediction_sha256(preds[a]) for a in GOLDEN_ARMS_FD18} == GOLDEN_ARMS_FD18
    shas = {orr.ls.prediction_sha256(p) for p in preds.values()}
    assert len(shas) == len(preds)
    for a in NEW:
        assert preds[a].X.shape == preds["m1_r_gmom_x4"].X.shape and np.isfinite(preds[a].X).all()
        assert not np.array_equal(preds[a].X, preds["m1_r_gmom_x4"].X)
    monkeypatch.setattr(orr, "MULTIPLIERS", {"idle": 1.0})
    idle = orr.scaled_arms(ctx, dirs, groups, 0.2, held, held)
    assert idle["idle"].X.tobytes() == preds["m1_r_gmom"].X.tobytes()


# ---------------------------------------------------------------------------------------------------------------------
# The reading


B0, BASE = 50.0, 53.0
W = {"de_score": 0.25, "de_direction": 0.25, "mmd_u": 0.30, "variogram": 0.20}
ROWS = {"pooled": 7919, "lib1": 1777, "lib2": 2711, "lib3": 3431, "lib3_t33": 3431, "lib3_t67": 3431,
        "lib3_t33b": 3431}
CONTAM = {"pooled": 613, "lib1": 401, "lib2": 199, "lib3": 13, "lib3_t33": 13, "lib3_t67": 13, "lib3_t33b": 13}


def _arm(score, de_part=0.5):
    gain = score - BASE
    sk = {"de_score": 0.4 + de_part * gain / 2 / (100 * W["de_score"]),
          "de_direction": 0.4 + de_part * gain / 2 / (100 * W["de_direction"]),
          "mmd_u": 0.3 + (1 - de_part) * gain * 0.6 / (100 * W["mmd_u"]),
          "variogram": 0.2 + (1 - de_part) * gain * 0.4 / (100 * W["variogram"])}
    return {"score": score, "groups": {"distribution": score / 100},
            "metrics": {"mmd_u": 0.1, "variogram": 0.01, "de_score": 0.5 + gain / 100, "de_direction": 0.3},
            "skills": sk}


def _default_gain(j: int, i: int) -> float:
    """A distinct tiny offset per arm (+d on seed 0, -d on seed 1), so that no two arms score alike on every seed."""
    d = 1e-9 * (j + 1)
    return d if i == 0 else -d if i == 1 else 0.0


def _task(truth, gains=None, per_seed=None, seeds=range(10), commit="c" * 40, clean=True, fractions=None,
          rows=None, contamination=None, sha=None, thinning=None, anchor_failures=None, identical=()):
    """A task file. gains: arm -> gain over m1_r_gmom on every seed; per_seed: arm -> ten gains. sha: a function
    (seed, arm) -> sha256 string (default: the same on every truth). identical: new arms that copy x4's record."""
    gains = gains or {}
    out = []
    for i, s in enumerate(seeds):
        arms = {"b0": _arm(B0), "m1_r_gmom": _arm(BASE)}
        for j, a in enumerate(SCALED):
            g = per_seed[a][i] if per_seed and a in per_seed else gains.get(a, _default_gain(j, i))
            arms[a] = _arm(BASE + g)
        for a in identical:
            arms[a] = dict(arms["m1_r_gmom_x4"])
        out.append({"seed": s, "floor": 50.0, "arms": arms,
                    "sha256": {a: (sha or (lambda s_, a_: f"{a_}-{s_}"))(s, a) for a in ALL_ARMS}})
    rec = {"label": "synthetic", "task": "T1/E9.5->E10.5:chamber", "truth": truth,
           "truth_rows": ROWS[truth] if rows is None else rows,
           "truth_contamination": CONTAM[truth] if contamination is None else contamination,
           "fraction": 0.2, "fractions": fractions or {"m1_r_gmom": 0.2, **{a: 0.2 * m for a, m in MULT.items()}},
           "multipliers": dict(MULT), "seeds": list(seeds), "rows": out,
           "stamp": {"commit": commit, "tree_clean": clean, "utc": "x", "inputs": {}}}
    if truth in THIN:
        p, stream = THIN[truth]
        rec["thinning"] = {"p": p, "stream": list(stream), "of": "lib3", "depth_ratio": p,
                           "rebuild_largest_difference": 0.0, "any_cell_emptied": False, **(thinning or {})}
    if anchor_failures is not None:
        rec["anchor_failures"] = anchor_failures
    return rec


def _write(tmp_path, truths, committed_shift=None):
    """truths: truth -> task record; T46's committed files for its four truths beside them, with the rescored arms'
    scores as in the run's files (committed_shift: (truth, arm) -> a shift of the committed score)."""
    committed_shift = committed_shift or {}
    out = tmp_path / "out"
    (out / "tasks").mkdir(parents=True)
    for k, rec in truths.items():
        (out / "tasks" / f"{k}.json").write_text(json.dumps(rec))
    t46 = tmp_path / "t46"
    t46.mkdir()
    for k in T46_TRUTHS:
        src = truths.get(k) or _task(k)
        rows = [{"seed": r["seed"], "arms": {a: {**r["arms"][a], "score": r["arms"][a]["score"]
                                                + committed_shift.get((k, a), 0.0)}
                                            for a in RESCORED if a in r["arms"]}} for r in src["rows"]]
        (t46 / f"{k}.json").write_text(json.dumps({"rows": rows}))
    return out, t46


def _read(orr, tmp_path, by_truth=None, over=None, committed_shift=None):
    """by_truth: truth -> gains or {"per_seed": ...}; over: truth -> task keyword overrides (None drops the truth)."""
    by_truth, over = by_truth or {}, over or {}
    truths = {}
    for k in TRUTHS:
        if k in over and over[k] is None:
            continue
        g = by_truth.get(k, {})
        args = {"per_seed": g["per_seed"]} if "per_seed" in g else {"gains": g}
        truths[k] = _task(k, **args, **over.get(k, {}))
    out, t46 = _write(tmp_path, truths, committed_shift)
    return orr.read(out, t46)["reading"]


def _rising(top=8.0):
    """Mean gains over m1_r_gmom rising with the multiplier up to `top` and falling after it."""
    return {a: (m if m <= top else 2 * top - m) for a, m in MULT.items()}


@pytest.mark.parametrize("gains, labels, best, statement", [
    (_rising(8.0), {"m1_r_gmom_x5": "ahead", "m1_r_gmom_x6": "ahead", "m1_r_gmom_x8": "ahead"}, 8.0,
     "still rising at 8"),
    (_rising(4.0), {"m1_r_gmom_x5": "behind", "m1_r_gmom_x6": "behind", "m1_r_gmom_x8": "behind"}, 4.0, "turns at 4"),
    ({**_rising(5.0), "m1_r_gmom_x6": 4.1}, {"m1_r_gmom_x5": "ahead", "m1_r_gmom_x6": "level",
                                              "m1_r_gmom_x8": "behind"}, 5.0, "turns at 5"),
    ({**_rising(4.0), "m1_r_gmom_x5": 4.2, "m1_r_gmom_x6": 4.05, "m1_r_gmom_x8": 3.8},
     {"m1_r_gmom_x5": "level", "m1_r_gmom_x6": "level", "m1_r_gmom_x8": "level"}, 5.0, "turns at 5"),
    ({**_rising(3.0), "m1_r_gmom_x4": 2.9, "m1_r_gmom_x5": 2.95, "m1_r_gmom_x6": 1.0, "m1_r_gmom_x8": -1.0},
     {"m1_r_gmom_x5": "level", "m1_r_gmom_x6": "behind", "m1_r_gmom_x8": "behind"}, 3.0, "turns at 3"),
])
def test_statement_one_each_label_and_the_argmax(orr, tmp_path, gains, labels, best, statement):
    r = _read(orr, tmp_path, {k: gains for k in TRUTHS})
    assert r["valid"]
    for k in TRUTHS:
        s1 = r["statement_1"][k]
        assert s1["against_x4"] == labels and s1["argmax"] == best and s1["statement"] == statement
    txt = (tmp_path / "out" / "t50.txt").read_text()
    assert f"on lib3: m1_r_gmom_x5 {labels['m1_r_gmom_x5']}" in txt and f"argmax {statement}" in txt
    assert not _prints_a_count(txt, (*ROWS.values(), *CONTAM.values()))


def test_a_gain_of_exactly_point_three_is_level_and_eight_seeds_are_needed(orr, tmp_path):
    # x5 - x4 per seed: 0.25 on eight seeds and 0.5 on two, exact in floating point; the mean is 0.3 exactly
    exact = {"m1_r_gmom_x4": [0.0] * 10, "m1_r_gmom_x5": [0.25] * 8 + [0.5] * 2}
    r = _read(orr, tmp_path, {"lib3": {"per_seed": exact}})
    c = r["over_comparisons"]["lib3"]["m1_r_gmom_x5"]
    assert r["valid"] and c["mean"] == 0.3 and c["seeds_ahead"] == 10 and c["label"] == "level"
    r = _read(orr, tmp_path / "b", {"lib3": {"per_seed": {**exact, "m1_r_gmom_x5": [0.25] * 8 + [0.75] * 2}}})
    assert r["over_comparisons"]["lib3"]["m1_r_gmom_x5"]["label"] == "ahead"
    seven = {"per_seed": {"m1_r_gmom_x5": [1.0] * 7 + [-0.1] * 3}}
    r = _read(orr, tmp_path / "c", {"lib1": seven})
    c = r["over_comparisons"]["lib1"]["m1_r_gmom_x5"]
    assert c["mean"] > 0.3 and c["seeds_ahead"] == 7 and c["label"] == "level"


def test_statement_two_reports_each_worst_case_with_its_label_and_parts(orr, tmp_path):
    g = {"m1_r_gmom_x3": 2.0, "m1_r_gmom_x4": 2.2, "m1_r_gmom_x5": 1.5, "m1_r_gmom_x6": 1.0, "m1_r_gmom_x8": 0.0}
    r = _read(orr, tmp_path, {"lib3": g})
    s2 = r["statement_2"]["lib3"]
    assert list(s2) == ["m1_r_gmom_x4 against m1_r_gmom_x3", "m1_r_gmom_x5 against m1_r_gmom_x3",
                        "m1_r_gmom_x6 against m1_r_gmom_x3", "m1_r_gmom_x6 against m1_r_gmom_x4",
                        "m1_r_gmom_x8 against m1_r_gmom_x4", "m1_r_gmom_x8 against m1_r_gmom_x5"]
    want = {"m1_r_gmom_x4 against m1_r_gmom_x3": (0.2, "level"),
            "m1_r_gmom_x5 against m1_r_gmom_x3": (-0.5, "behind"),
            "m1_r_gmom_x6 against m1_r_gmom_x3": (-1.0, "behind"),
            "m1_r_gmom_x6 against m1_r_gmom_x4": (-1.2, "behind"),
            "m1_r_gmom_x8 against m1_r_gmom_x4": (-2.2, "behind"),
            "m1_r_gmom_x8 against m1_r_gmom_x5": (-1.5, "behind")}
    for name, (mean, lab) in want.items():
        c = s2[name]
        assert c["mean"] == pytest.approx(mean, abs=1e-8) and c["label"] == lab
        assert c["seeds_ahead"] == (0 if mean < 0 else 10) and len(c["range"]) == 2
        assert c["points"]["de"] == pytest.approx(mean / 2, abs=1e-8)
        assert c["points"]["mmd_u_variogram"] == pytest.approx(mean / 2, abs=1e-8)
        assert set(c["skills"]) == set(W) and "de_score" in c["metrics"]
    assert set(r["statement_2"]) == set(TRUTHS)
    up = _read(orr, tmp_path / "b", {"lib3": {**g, "m1_r_gmom_x4": 2.5}})
    assert up["statement_2"]["lib3"]["m1_r_gmom_x4 against m1_r_gmom_x3"]["label"] == "ahead"
    txt = (tmp_path / "out" / "t50.txt").read_text()
    assert "on lib3: m1_r_gmom_x8 against m1_r_gmom_x5: -1.50 points, ahead 0 of 10" in txt
    assert "test: behind; de terms -0.75, mmd_u and variogram -0.75; de_score skill" in txt


def _depth(x4_minus_x2_thin, x4_minus_x2_lib3=1.0, x8_extra_thin=0.37, null=0.0):
    """Gains with (x4 - x2) set per truth: lib3's and each thinned truth's; x8 moved by x8_extra_thin on the thinned
    truths only (so x8 - x4 differs from x6 - x4 between truths); lib3_t33b's (x4 - x2) is lib3_t33's plus null."""
    def g(d, e=0.0):
        return {**{a: 0.0 for a in SCALED}, "m1_r_gmom_x2": 1.0, "m1_r_gmom_x4": 1.0 + d,
                "m1_r_gmom_x5": 0.5 + d, "m1_r_gmom_x6": 0.4 + d, "m1_r_gmom_x8": 0.3 + d + e}
    return {"lib3": g(x4_minus_x2_lib3), "lib3_t33": g(x4_minus_x2_thin, x8_extra_thin),
            "lib3_t67": g(x4_minus_x2_thin, x8_extra_thin), "lib3_t33b": g(x4_minus_x2_thin + null, x8_extra_thin)}


@pytest.mark.parametrize("thin, statement, label", [
    (2.0, "The larger steps gain more against the shallower truth", "ahead"),
    (0.0, "The larger steps gain less against the shallower truth", "behind"),
    (1.2, "no difference at the bar", "level"),
])
def test_statement_three_each_depth_label(orr, tmp_path, thin, statement, label):
    r = _read(orr, tmp_path, _depth(thin))
    s3 = r["statement_3"]
    assert s3["valid"] and s3["draw_null"]["holds"] and s3["draw_null"]["x4_minus_x2"]["label"] == "level"
    assert s3["D33"]["statement"] == statement and s3["D33"]["label"] == label
    assert s3["D67"]["statement"] == statement and s3["D67"]["label"] == label
    assert s3["D33"]["mean"] == pytest.approx(thin - 1.0, abs=1e-8)
    for k in THIN:
        d = s3["x8_minus_x4"][k]
        # x8 - x4 is -0.7 on lib3 and -0.7 + 0.37 on the thinned truths; x6 - x4 is -0.6 on every truth
        assert d["mean"] == pytest.approx(0.37, abs=1e-8)
    assert s3["draw_null"]["x8_minus_x4"]["mean"] == pytest.approx(0.0, abs=1e-8)
    assert f"; depth: {statement}" in r["outcome"]
    assert set(s3["reported_not_read"]) == {"lib1", "lib2"}
    assert s3["reported_not_read"]["lib1"]["mean"] == pytest.approx(-1.0, abs=1e-8)  # lib1's 0 minus lib3's 1
    txt = (tmp_path / "out" / "t50.txt").read_text()
    assert "Statement 3, depth: D33, " in txt and statement in txt
    assert "read by nothing (different truths, not paired in kind)" in txt


def test_d33_and_d67_are_read_separately_and_paired_by_seed(orr, tmp_path):
    per = {a: [0.0] * 10 for a in SCALED}
    lib3 = {**per, "m1_r_gmom_x2": [1.0] * 10, "m1_r_gmom_x4": [float(i) for i in range(10)]}
    t33 = {**per, "m1_r_gmom_x2": [1.0] * 10, "m1_r_gmom_x4": [float(i) + 0.5 for i in range(10)]}
    t67 = {**per, "m1_r_gmom_x2": [1.0] * 10, "m1_r_gmom_x4": [float(i) + 0.1 for i in range(10)]}
    for d in (lib3, t33, t67):
        d["m1_r_gmom_x5"] = [x + 1e-3 for x in d["m1_r_gmom_x4"]]
        d["m1_r_gmom_x6"] = [x + 2e-3 for x in d["m1_r_gmom_x4"]]
        d["m1_r_gmom_x8"] = [x + 3e-3 for x in d["m1_r_gmom_x4"]]
    r = _read(orr, tmp_path, {"lib3": {"per_seed": lib3}, "lib3_t33": {"per_seed": t33},
                              "lib3_t67": {"per_seed": t67}, "lib3_t33b": {"per_seed": t33}})
    s3 = r["statement_3"]
    assert s3["D33"]["mean"] == pytest.approx(0.5) and s3["D33"]["seeds_ahead"] == 10 and s3["D33"]["sd"] < 1e-9
    assert s3["D33"]["statement"] == "The larger steps gain more against the shallower truth"
    assert s3["D67"]["mean"] == pytest.approx(0.1) and s3["D67"]["statement"] == "no difference at the bar"


@pytest.mark.parametrize("null, label", [(0.5, "ahead"), (-0.5, "behind")])
def test_statement_three_is_void_when_the_second_draw_is_not_level_against_the_first(orr, tmp_path, null, label):
    """The draw null (lib3_t33b against lib3_t33, the same predictions) past the bar: one draw alone can pass it, so
    statement 3 is void and only reported; statements 1 and 2 stand."""
    r = _read(orr, tmp_path, _depth(2.0, null=null))
    assert r["valid"] and "statement_1" in r and "statement_2" in r
    s3 = r["statement_3"]
    assert not s3["valid"] and not s3["draw_null"]["holds"]
    assert s3["draw_null"]["x4_minus_x2"]["label"] == label
    assert s3["draw_null"]["x4_minus_x2"]["mean"] == pytest.approx(null, abs=1e-8)
    assert s3["D33"]["statement"] == orr.DEPTH_VOID and s3["D67"]["statement"] == orr.DEPTH_VOID
    assert s3["D33"]["label"] == "ahead"  # still reported
    assert s3["D33"]["label_if_valid"] == "The larger steps gain more against the shallower truth"
    assert r["outcome"].endswith(f"; depth: {orr.DEPTH_VOID}")
    txt = (tmp_path / "out" / "t50.txt").read_text()
    assert "statement 3 is void and only reported" in txt and orr.DEPTH_VOID in txt


def test_statement_three_stands_with_a_level_null_printed_beside_d33(orr, tmp_path):
    per = {a: [0.0] * 10 for a in SCALED}
    base = {**per, "m1_r_gmom_x2": [1.0] * 10, "m1_r_gmom_x4": [2.0] * 10}
    t33 = {**base, "m1_r_gmom_x4": [3.0] * 10}
    t33b = {**base, "m1_r_gmom_x4": [3.0 + (0.2 if i % 2 else -0.1) for i in range(10)]}  # null: mean +0.05, level
    for d in (base, t33, t33b):
        for j, a in enumerate(("m1_r_gmom_x5", "m1_r_gmom_x6", "m1_r_gmom_x8")):
            d[a] = [x + 1e-3 * (j + 1) for x in d["m1_r_gmom_x4"]]
    r = _read(orr, tmp_path, {"lib3": {"per_seed": base}, "lib3_t33": {"per_seed": t33},
                              "lib3_t67": {"per_seed": t33}, "lib3_t33b": {"per_seed": t33b}})
    s3 = r["statement_3"]
    assert s3["valid"] and s3["draw_null"]["x4_minus_x2"]["label"] == "level"
    assert s3["draw_null"]["x4_minus_x2"]["mean"] == pytest.approx(0.05)
    assert s3["D33"]["statement"] == "The larger steps gain more against the shallower truth"
    txt = (tmp_path / "out" / "t50.txt").read_text()
    assert "statement 3 stands" in txt
    assert "beside it the draw null's mean +0.05, per seed -0.10 to +0.20" in txt


def test_the_tie_rule_picks_the_smallest_multiplier(orr):
    means = {"m1_r_gmom": 53.0, "m1_r_gmom_x15": 54.0, "m1_r_gmom_x2": 55.0, "m1_r_gmom_x3": 56.0,
             "m1_r_gmom_x4": 57.0, "m1_r_gmom_x5": 57.0, "m1_r_gmom_x6": 56.5, "m1_r_gmom_x8": 57.0}
    t = orr.turn(means)
    assert t["best"] == 4.0 and t["statement"] == "turns at 4"
    t = orr.turn({**means, "m1_r_gmom_x4": 56.0, "m1_r_gmom_x5": 56.0, "m1_r_gmom_x8": 57.0})
    assert t["best"] == 8.0 and t["statement"] == "still rising at 8"
    t = orr.turn({**means, "m1_r_gmom": 57.0})
    assert t["best"] == 1.0


def test_the_outcome_line_gives_the_labels_first_and_marks_an_unseparated_argmax(orr, tmp_path):
    flat = {**_rising(4.0), "m1_r_gmom_x5": 4.2, "m1_r_gmom_x6": 4.05, "m1_r_gmom_x8": 3.8}  # all level, argmax 5
    r = _read(orr, tmp_path, {**{k: _rising(8.0) for k in TRUTHS}, "lib3": flat})
    assert "lib3 x5 level, x6 level, x8 level; argmax 5 (not separated)" in r["outcome"]
    assert "lib1 x5 ahead, x6 ahead, x8 ahead; argmax 8 (still rising at 8)" in r["outcome"]
    txt = (tmp_path / "out" / "t50.txt").read_text()
    assert "on lib3: m1_r_gmom_x5 level, m1_r_gmom_x6 level, m1_r_gmom_x8 level; argmax turns at 5 (not separated)" \
        in txt
    r = _read(orr, tmp_path / "b", {k: _rising(4.0) for k in TRUTHS})
    assert "lib3 x5 behind, x6 behind, x8 behind; argmax 4;" in r["outcome"]
    assert "(not separated)" not in r["outcome"]


def test_statement_three_reports_the_raw_metric_differences_the_anchor_gaps_and_the_shapes(orr, tmp_path):
    anchors = {m: {"floor": 0.1, "ceiling": 0.6} for m in W}
    thin_anchors = {m: {"floor": 0.2, "ceiling": 0.45} for m in W}
    truths = {k: _task(k, gains=g) for k, g in {**{k: {} for k in TRUTHS}, **_depth(2.0)}.items()}
    for k, rec in truths.items():
        for row in rec["rows"]:
            if k == "lib3":
                row["anchors"] = anchors
            elif k in THIN:
                row["anchors"] = thin_anchors
    out, t46 = _write(tmp_path, truths)
    r = orr.read(out, t46)["reading"]
    s3 = r["statement_3"]
    raw = s3["D33"]["raw_by_metric"]
    assert raw["de_score"]["mean"] == pytest.approx(1.0 / 100)  # D33 of +1 point is +0.01 of raw de_score here
    assert raw["mmd_u"]["mean"] == pytest.approx(0.0) and set(raw) == set(W)
    assert s3["anchor_gap"]["lib3"] == pytest.approx({m: 0.5 for m in W})
    assert s3["anchor_gap"]["lib3_t33"] == pytest.approx({m: 0.25 for m in W})
    assert set(s3["shapes"]) == {"lib3", *THIN}
    assert s3["shapes"]["lib3_t33"] == r["statement_1"]["lib3_t33"]
    txt = (tmp_path / "out" / "t50.txt").read_text()
    assert "D33 on the raw metrics (free of each truth's anchors; read by nothing): de_score +0.01" in txt
    assert "lib3: mean ceiling minus floor per metric (read by nothing): de_score +0.5" in txt
    assert "the curve's shape side by side" in txt
    r = _read(orr, tmp_path / "b", _depth(2.0))
    assert r["statement_3"]["anchor_gap"]["lib3"] is None
    assert "lib3: mean ceiling minus floor per metric (read by nothing): not recorded" in \
        (tmp_path / "b" / "out" / "t50.txt").read_text()


def test_the_run_records_each_seeds_anchors(orr):
    rec = orr.anchor_record({"mmd_u": {"floor": np.float32(0.25), "ceiling": 0.5, "other": 1}})
    assert rec == {"mmd_u": {"floor": 0.25, "ceiling": 0.5}} and type(rec["mmd_u"]["floor"]) is float
    assert '"anchors": anchor_record(f.anchors)' in inspect.getsource(orr.cmd_run)


def test_the_new_limits_are_in_the_text(orr):
    for want in ("No outcome of T50 is evidence on T1 read depth, a lever closed after an exposure that this read does "
                 "not reopen.",
                 "Board points are not comparable across truths: each truth, the thinned ones included, has its own "
                 "floor and ceiling",
                 "the argmax alone is not a finding.",
                 "'level' is not 'no loss'.",
                 "differ from it in contamination share and capture",
                 "statement 3 stands only when a second draw at p = 1/3 (lib3_t33b) reads level against the first.",
                 "it shows no biological change and no change in Feng's share",
                 "The depth-ratio check catches a wrong p or a wrong stream only"):
        assert any(want in x for x in orr.LIMITS), want


def test_the_thinning_record_is_printed(orr, tmp_path):
    _read(orr, tmp_path, over={"lib3_t33b": {"thinning": {"rebuild_largest_difference": 2.5e-7}}})
    txt = (tmp_path / "out" / "t50.txt").read_text()
    assert ("truth lib3_t33b: realised depth ratio 0.3333 (p 0.3333, stream [50, 3]), rebuild check's largest "
            "difference 2.5e-07, any cell left with no count: False") in txt


def test_a_thinning_that_empties_a_cell_is_recorded(orr, world, tmp_path, monkeypatch):
    t, _, _ = world

    def empty_first(C, name):
        out = sparse.csr_matrix(C, dtype=np.int64, copy=True)
        out[0] = 0
        out.eliminate_zeros()
        return out

    monkeypatch.setattr(orr, "thin_counts", empty_first)
    th, _, rec = orr.thinned_truth(t, tmp_path, PANEL, "lib3_t33")
    assert rec["any_cell_emptied"] is True
    assert th.X[0].sum() == 0  # the all-zero row the scorer would get


def test_a_file_with_no_rows_key_is_void_not_an_error(orr, tmp_path):
    truths = {k: _task(k) for k in TRUTHS}
    out, t46 = _write(tmp_path, truths)
    (out / "tasks" / "lib1.json").write_text(json.dumps({k: v for k, v in truths["lib1"].items() if k != "rows"}))
    r = orr.read(out, t46)["reading"]
    assert not r["valid"] and any("lib1: seeds []" in x for x in r["validity"]["reasons"])


def test_the_rebuild_refuses_rows_whose_cell_types_differ_from_lib3s(orr, monkeypatch, tmp_path):
    """The full file read by name gives rows whose cell types are not lib3's: refused."""
    a, _ = _counts_target()
    t = _chamber(monkeypatch, orr, a)
    shuffled = a.copy()
    shuffled.obs["celltype"] = np.random.default_rng(1).permutation(shuffled.obs["celltype"].to_numpy())
    monkeypatch.setattr(orr, "read_stage", lambda path, genes=None: shuffled.copy() if genes is None
                        else shuffled[:, list(genes)].copy())
    with pytest.raises(ValueError, match="the rebuilt rows' cell types differ from lib3's"):
        orr.thinned_truth(t, tmp_path, PANEL, "lib3_t33")


def test_the_reported_quantities(orr, tmp_path):
    r = _read(orr, tmp_path, {k: _rising(8.0) for k in TRUTHS},
              over={"lib3_t33": {"thinning": {"depth_ratio": 0.335}}})
    assert r["valid"] and r["depth_ratio"] == pytest.approx({"lib3_t33": 0.335, "lib3_t67": 2 / 3,
                                                             "lib3_t33b": 1 / 3})
    n = r["neighbours"]["lib2"]
    assert n["m1_r_gmom_x15"]["base"] == "m1_r_gmom" and n["m1_r_gmom_x8"]["base"] == "m1_r_gmom_x6"
    assert n["m1_r_gmom_x15"]["per_tenth"] == pytest.approx(1.5)  # +1.5 over 0.1 of fraction
    assert n["m1_r_gmom_x5"]["per_tenth"] == pytest.approx(0.5)  # +1 over 0.2
    assert n["m1_r_gmom_x8"]["per_tenth"] == pytest.approx(0.5)  # +2 over 0.4
    assert r["mean_scores"]["lib1"]["m1_r_gmom_x8"] == pytest.approx(BASE + 8)
    assert r["sha256_across_truths"]["all_agree"]
    txt = (tmp_path / "out" / "t50.txt").read_text()
    assert "truth lib3_t33: realised depth ratio 0.3350 (p 0.3333, stream [50, 1])" in txt
    assert "m1_r_gmom_x8 against m1_r_gmom_x6 +2.00 (ahead), +0.50 per 0.1 of fraction" in txt
    assert "predictions' sha256 across the truths: the same on every seed and arm" in txt
    for line in orr.LIMITS:
        assert f"Limit: {line}" in txt
    assert not _prints_a_count(txt, (*ROWS.values(), *CONTAM.values()))
    r = _read(orr, tmp_path / "b", over={"lib1": {"sha": lambda s, a: f"lib1-{a}-{s}"}})
    assert r["valid"] and not r["sha256_across_truths"]["all_agree"]


FAILURE = [{"seed": 3, "reason": "T1/E9.5->E10.5:chamber seed 3: inverted"}]
FR_X6 = {"m1_r_gmom": 0.2, **{a: 0.2 * m for a, m in MULT.items()}, "m1_r_gmom_x6": 1.0}


@pytest.mark.parametrize("over, shift, reason", [
    ({"lib2": None}, None, "missing truths ['lib2']"),
    ({"lib3_t67": None}, None, "missing truths ['lib3_t67']"),
    ({"lib3_t33b": None}, None, "missing truths ['lib3_t33b']"),
    ({"lib3_t33b": {"thinning": {"stream": [50, 1]}}}, None, "lib3_t33b: the file records p"),
    ({"lib3_t33b": {"sha": lambda s, a: "other" if (s, a) == (0, "b0") else f"{a}-{s}"}}, None,
     "lib3_t33b: the predictions' sha256 do not equal lib3's"),
    ({"lib3_t33": {"thinning": {"any_cell_emptied": True}}}, None, "lib3_t33: the thinning left a truth cell"),
    ({"lib3_t67": {"thinning": {"any_cell_emptied": None}}}, None, "lib3_t67: the thinning left a truth cell"),
    ({"lib3_t33b": {"thinning": {"rebuild_largest_difference": 1e-4}}}, None,
     "lib3_t33b: the rebuild check's largest difference"),
    ({"lib1": {"seeds": range(9)}}, None, "lib1: seeds"),
    ({"lib3_t33": {"seeds": range(9)}}, None, "lib3_t33: seeds"),
    ({}, {("lib2", "m1_r_gmom_x3"): 1e-6}, "do not reproduce T46's committed scores"),
    ({}, {("pooled", "b0"): 1e-6}, "do not reproduce T46's committed scores"),
    ({"lib3_t33": {"fractions": FR_X6}}, None, "wrong fraction on ['m1_r_gmom_x6']"),
    ({"lib1": {"commit": "d" * 40}}, None, "2 commits"),
    ({"lib3_t67": {"clean": False}}, None, "tree not clean on ['lib3_t67']"),
    ({"lib2": {"clean": None}}, None, "tree not clean"),
    ({"lib3_t33": {"sha": lambda s, a: "other" if (s, a) == (4, "m1_r_gmom_x5") else f"{a}-{s}"}}, None,
     "lib3_t33: the predictions' sha256 do not equal lib3's seed for seed (seed:arm ['4:m1_r_gmom_x5'])"),
    ({"lib3_t67": {"thinning": {"depth_ratio": 2 / 3 + 0.02}}}, None, "lib3_t67: realised depth ratio"),
    ({"lib3_t33": {"thinning": {"depth_ratio": 1 / 3 - 0.02}}}, None, "lib3_t33: realised depth ratio"),
    ({"lib3_t33": {"thinning": {"p": 0.5}}}, None, "lib3_t33: the file records p"),
    ({"lib3_t67": {"thinning": {"stream": [50, 1]}}}, None, "lib3_t67: the file records p"),
    ({"lib1": {"identical": ("m1_r_gmom_x6",)}}, None, "new arms score exactly as m1_r_gmom_x4"),
    ({"lib3_t67": {"identical": ("m1_r_gmom_x8",)}}, None, "new arms score exactly as m1_r_gmom_x4"),
    ({"lib2": {"anchor_failures": FAILURE, "seeds": ()}}, None, "anchors invert or go missing on seeds [3]"),
    ({"lib3_t33": {"rows": 3430}}, None, "lib3_t33 does not equal lib3"),
    ({"lib3": {"rows": 3430}}, None, "do not add up"),
    ({"lib2": {"contamination": 198}}, None, "contamination cells: the library truths do not add up"),
])
def test_the_void_cases(orr, tmp_path, over, shift, reason):
    r = _read(orr, tmp_path, {"lib3": _rising(4.0)}, over=over, committed_shift=shift)
    assert not r["valid"] and r["outcome"] == "VOID"
    assert any(reason in x for x in r["validity"]["reasons"]), r["validity"]["reasons"]
    txt = (tmp_path / "out" / "t50.txt").read_text()
    assert "Void: " in txt and not _prints_a_count(txt, (*ROWS.values(), *CONTAM.values(), 3430, 198))


def test_a_rescored_arm_within_tolerance_and_a_ratio_within_point_zero_one_stay_valid(orr, tmp_path):
    r = _read(orr, tmp_path, over={"lib3_t33": {"thinning": {"depth_ratio": 1 / 3 + 0.009}}},
              committed_shift={("lib1", "m1_r_gmom_x4"): 1e-12})
    assert r["valid"]


def test_the_reading_refuses_to_overwrite(orr, tmp_path):
    truths = {k: _task(k) for k in TRUTHS}
    out, t46 = _write(tmp_path, truths)
    orr.read(out, t46)
    assert (out / "t50.json").exists() and (out / "t50.txt").exists()
    with pytest.raises(SystemExit, match="already holds"):
        orr.read(out, t46)


# ---------------------------------------------------------------------------------------------------------------------
# The run's refusals


def test_the_run_needs_its_row(orr, monkeypatch):
    monkeypatch.setattr(orr.n12, "tries_rows", lambda: {"T46", "T49"})
    monkeypatch.setattr(orr.gc, "n22_passed", lambda path: True)
    with pytest.raises(SystemExit, match="T50 has no row"):
        orr.main(["run", "--truth", "lib3_t33", "--out", "t50"])


def test_the_run_needs_n22(orr, monkeypatch):
    monkeypatch.setattr(orr.n12, "tries_rows", lambda: {"T50"})
    with pytest.raises(SystemExit, match="does not hold a passing N22"):
        orr.main(["run", "--truth", "pooled", "--out", "t50", "--n22", "n22x"])


@pytest.mark.parametrize("truth", TRUTHS)
def test_the_run_refuses_an_existing_task_file(orr, monkeypatch, truth):
    monkeypatch.setattr(orr.n12, "tries_rows", lambda: {"T50"})
    monkeypatch.setattr(orr.gc, "n22_passed", lambda path: True)
    tasks = orr.RESULTS / "t50" / "tasks"
    tasks.mkdir(parents=True)
    (tasks / f"{truth}.json").write_text("{}")
    with pytest.raises(SystemExit, match="already holds"):
        orr.main(["run", "--truth", truth, "--out", "t50"])
