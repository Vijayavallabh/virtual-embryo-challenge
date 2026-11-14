"""experiments/t1_m1/perm_read.py (try T45 of experiments/harness/TRIES.md) on synthetic data: the control arms and the
reading. No released file, statistics file or score is read."""

from __future__ import annotations

import copy
import importlib.util
import json
import sys

import anndata as ad
import numpy as np
import pytest
from scipy import sparse
from test_t1_step_read import _arms_world

from vec import harness, t1_scope
from vec import t1_state_transport as ts
from vec.boards import REPO_ROOT
from vec.core import counts_of, log1p_cp10k

SCRIPT = REPO_ROOT / "experiments" / "t1_m1" / "perm_read.py"


@pytest.fixture
def pr(monkeypatch, tmp_path):
    spec = importlib.util.spec_from_file_location("perm_read_module", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    monkeypatch.setattr(mod, "RESULTS", tmp_path / "results")
    return mod


# ---------------------------------------------------------------------------------------------------------------------
# The placebo


def _dirs_and_halves(seed=0, n=203):
    rng = np.random.default_rng(seed)
    dirs, halves = {}, {}
    for g in ("CM", "Endocardium", "Endothelium", "Epicardium"):
        d = rng.normal(0, 0.5, n)
        d[rng.random(n) < 0.1] = np.nan
        mean_cp = rng.gamma(1.0, 5.0, n)
        dirs[g], halves[g] = d, {"mean_cp": mean_cp, "d_a": d + 0.1, "d_b": d - 0.1}
    return dirs, halves


def test_the_placebo_is_try_t15s_permuted_control_on_m1s_directions(pr):
    dirs, halves = _dirs_and_halves()
    pl = pr.placebo_directions(dirs, halves)
    assert set(pl) == set(dirs)
    for g, d in dirs.items():
        ok = np.isfinite(d)
        want = d.copy()
        want[ok] = harness.t3_decile_permutation(d[ok], halves[g]["mean_cp"][ok],
                                                 [harness.LP_PERM_SEED, harness.LP_GROUPS.index(g)], harness.LP_BINS)
        assert np.array_equal(pl[g], want, equal_nan=True)
        assert np.array_equal(np.isnan(pl[g]), np.isnan(d)), "NaN where Feng has NaN, and only there"
        idx = np.flatnonzero(ok)
        for run in np.array_split(idx[np.argsort(halves[g]["mean_cp"][idx], kind="stable")], harness.LP_BINS):
            assert np.array_equal(np.sort(pl[g][run]), np.sort(d[run])), "each bin keeps its multiset of values"
        assert not np.array_equal(pl[g], d, equal_nan=True)
    assert pr.direction_sha(pl) != pr.direction_sha(dirs)
    assert all(pr.direction_sha(pl)[g] != pr.direction_sha(dirs)[g] for g in dirs)
    # the inputs are untouched and a second call is the same permutation
    assert np.array_equal(pr.placebo_directions(dirs, halves)["CM"], pl["CM"], equal_nan=True)


def test_the_placebo_refuses_a_group_without_a_permutation_seed_or_halves(pr):
    dirs, halves = _dirs_and_halves()
    with pytest.raises(ValueError, match="not one of try T15's groups"):
        pr.placebo_directions({"Blood": dirs["CM"]}, {"Blood": halves["CM"]})
    with pytest.raises(ValueError, match="no mean CP10k"):
        pr.placebo_directions({"CM": dirs["CM"]}, {})


# ---------------------------------------------------------------------------------------------------------------------
# The arms (spy)


def test_the_arms_differ_from_the_base_build_in_the_fraction_and_the_directions_only(pr, monkeypatch):
    calls = []

    def spy(ctx, directions, groups, fraction, held_m1, held_factors):
        calls.append({"ctx": ctx, "directions": directions, "groups": groups, "fraction": fraction,
                      "held_m1": held_m1, "held_factors": held_factors})
        return {"m1_r": ("m1_r", fraction), "m1_r_gmom": ("m1_r_gmom", len(calls))}

    monkeypatch.setattr(pr.gc, "m1_r_arms", spy)
    groups, held_m1, held_f = {"V-CM": "CM"}, ["a", "b"], ["a"]
    dirs, placebo = {"CM": np.zeros(3)}, {"CM": np.ones(3)}

    class Ctx:
        def __init__(self, seed):
            self.seed = seed

    seen = {}
    for seed in (1000000, 1000003):
        calls.clear()
        ctx = Ctx(seed)
        arms = pr.perm_arms(ctx, dirs, placebo, groups, 0.2, held_m1, held_f)
        assert set(arms) == {"m1_r_gmom", *pr.CONTROLS, "m1_r_gmom_x2"}
        assert [c["fraction"] for c in calls] == [0.2, 0.4, 0.0, 0.2, 0.4, 0.4]
        assert [c["directions"] for c in calls] == [dirs, dirs, dirs, placebo, placebo, dirs]
        assert all(c["directions"] is (placebo if i in (3, 4) else dirs) for i, c in enumerate(calls))
        for c in calls:
            for key in ("ctx", "groups", "held_m1", "held_factors"):
                assert c[key] is calls[0][key], f"an arm differs from the base build in {key}"
        assert arms["m1_r_gmom_p1"] == ("m1_r_gmom", 4) and arms["m1_r_gmom_p2"] == ("m1_r_gmom", 5)
        assert arms["m1_r_gmom_x2b"] == ("m1_r_gmom", 6)
        seen[seed] = [c["directions"] for c in calls[3:5]]
    # the same placebo object at both fractions and for two build seeds
    assert all(d is placebo for ds in seen.values() for d in ds)
    assert pr.arm_fractions(0.2) == {"m1_r_gmom": 0.2, "m1_r_gmom_x2": 0.4, "m1_r_gmom_x0": 0.0,
                                     "m1_r_gmom_p1": 0.2, "m1_r_gmom_p2": 0.4, "m1_r_gmom_x2b": [0.4, -0.2]}


# ---------------------------------------------------------------------------------------------------------------------
# The arms (small synthetic builds)


def _spy_shift(monkeypatch, record):
    real = t1_scope.external_shift

    def spy(target_like, directions, fraction, group_of, held_flat=(), **kw):
        out, rep = real(target_like, directions, fraction, group_of, held_flat, **kw)
        record.append({"before": counts_of(target_like).toarray(), "after": counts_of(out).toarray(),
                       "fraction": fraction, "seed": list(np.atleast_1d(kw.get("seed"))), "held": list(held_flat),
                       "labels": target_like.obs[t1_scope.LABEL_KEY].astype(str).to_numpy(), "report": rep})
        return out, rep

    monkeypatch.setattr(t1_scope, "external_shift", spy)


def test_fraction_zero_redraws_no_count(pr, monkeypatch):
    """T45a: external_shift accepts a fraction of 0; every factor is exactly 1 and the covered cells' counts before
    g_mom equal the drawn cells'."""
    ctx, groups, dirs, held = _arms_world()
    record, factors = [], []
    _spy_shift(monkeypatch, record)
    real_apply = t1_scope.apply_factors

    def spy_apply(C, F, tidx, held_mask, **kw):
        factors.append(np.asarray(F).copy())
        return real_apply(C, F, tidx, held_mask, **kw)

    monkeypatch.setattr(t1_scope, "apply_factors", spy_apply)
    arms = pr.sr.step_arms(ctx, dirs, groups, 0.2, held, held, multipliers={pr.X0: 0.0})
    zero = [r for r in record if r["fraction"] == 0.0]
    assert zero and all(np.array_equal(r["before"], r["after"]) for r in zero)
    assert any(np.array_equal(f, np.ones_like(f)) for f in factors)
    assert all(r["report"].rows[i]["max_abs_log_f"] == 0.0 for r in zero for i in range(len(r["report"].rows) - 1))
    assert not np.array_equal(arms[pr.X0].X, arms["m1_r_gmom"].X)


def test_there_and_back_second_edit_covered_cells_held_genes_new_stream_first_edit_is_x2s(pr, monkeypatch):
    ctx, groups, dirs, held = _arms_world()
    held_m1 = list(held)  # R2's held set: here the same genes, a separate list
    record, gmom_in, gmom_streams, moved = [], [], [], []
    _spy_shift(monkeypatch, record)
    real_hook, real_moved, real_asf = pr.gc.gmom_hook, pr.gc.moved_by_feng, ts.apply_state_factors

    def spy_hook(*a, **k):
        h = real_hook(*a, **k)

        def f(pop):
            gmom_in.append(counts_of(pop).toarray())
            return h(pop)
        return f

    def spy_moved(pop, directions, fraction, groups_, held_, seed):
        out = real_moved(pop, directions, fraction, groups_, held_, seed)
        moved.append({"before": counts_of(pop).toarray(), "after": counts_of(out).toarray(), "seed": list(seed),
                      "fraction": fraction, "held": list(held_), "directions": directions,
                      "labels": pop.obs[t1_scope.LABEL_KEY].astype(str).to_numpy()})
        return out

    def spy_asf(pop, factors, held_, seed, **kw):
        gmom_streams.append(list(seed))
        return real_asf(pop, factors, held_, seed, **kw)

    monkeypatch.setattr(pr.gc, "gmom_hook", spy_hook)
    monkeypatch.setattr(pr.gc, "moved_by_feng", spy_moved)
    monkeypatch.setattr(ts, "apply_state_factors", spy_asf)

    x2 = pr.sr.step_arms(ctx, dirs, groups, 0.2, held_m1, held, multipliers={pr.X2: 2.0})[pr.X2]
    x2_first = [r for r in record if r["fraction"] == 0.4][-1]
    n_before = len(record)
    gmom_in.clear()
    gmom_streams.clear()
    x2b = pr.there_and_back(ctx, dirs, groups, 0.2, held_m1, held)
    assert pr.gc.gmom_hook is spy_hook, "the hook is restored"
    new = record[n_before:]
    first = [r for r in new if r["fraction"] == 0.4]
    second = [r for r in new if r["fraction"] == -0.2]
    assert first and len(second) == 1 and len(moved) == 1
    # the first edit is bit for bit the x2 arm's
    assert all(np.array_equal(f["after"], x2_first["after"]) and f["seed"] == x2_first["seed"] for f in first)
    # the second edit: Feng's directions, R2's held genes, covered cells only, on a new stream
    m = moved[0]
    assert m["directions"] is dirs and m["held"] == list(held_m1) and m["fraction"] == -0.2
    assert m["seed"] == [ctx.seed, pr.SECOND_STREAM] == second[0]["seed"]
    assert pr.SECOND_STREAM not in {1, 2, 4, 5}
    covered = np.array([groups.get(lab) in dirs for lab in m["labels"]])
    assert covered.any() and (~covered).any()
    assert np.array_equal(m["after"][~covered], m["before"][~covered]), "uncovered cells as they are"
    genes = list(ctx.inputs["ref"].var_names)
    hcols = np.array([g in set(held_m1) for g in genes])
    assert hcols.any() and np.array_equal(m["after"][:, hcols], m["before"][:, hcols]), "held genes held"
    assert not np.array_equal(m["after"][covered], m["before"][covered])
    # g_mom's factors act on the second edit's output, last
    assert np.array_equal(gmom_in[-1], m["after"])
    streams = [s[-1] for s in (*[r["seed"] for r in new], *gmom_streams)]
    assert sorted(set(streams)) == sorted({2, pr.SECOND_STREAM, 5})
    assert x2b.X.shape == x2.X.shape and not np.array_equal(x2b.X, x2.X)


def _toy(n_cells=12, n_genes=24, seed=3):
    rng = np.random.default_rng(seed)
    C = rng.integers(40, 400, size=(n_cells, n_genes))
    labels = np.array(["V-CM"] * (n_cells - 3) + ["Foregut"] * 3)
    a = ad.AnnData(X=log1p_cp10k(sparse.csr_matrix(C)), obs={t1_scope.LABEL_KEY: labels})
    a.var_names = [f"g{i}" for i in range(n_genes)]
    a.obs_names = [f"c{i}" for i in range(n_cells)]
    a.layers["counts"] = sparse.csr_matrix(C.astype(np.int32))
    d = rng.normal(0, 0.6, n_genes)
    d[0] = np.nan
    return a, {"CM": d}, {"V-CM": "CM", "Foregut": "other:Foregut"}, ["g1", "g2"]


def test_there_and_backs_expected_counts_are_the_bases(pr):
    """Over many draws the two edits' mean counts equal one edit at the net fraction's, to the depth constants."""
    a, dirs, groups, held = _toy()
    draws = 400
    base = np.zeros(a.shape)
    back = np.zeros(a.shape)
    for i in range(draws):
        base += counts_of(pr.gc.moved_by_feng(a, dirs, 0.2, groups, held, [i, 2])).toarray()
        first = pr.gc.moved_by_feng(a, dirs, 0.4, groups, held, [i, 2])
        back += counts_of(pr.second_edit(dirs, groups, -0.2, held, i)(first)).toarray()
    base, back = base / draws, back / draws
    C = counts_of(a).toarray()
    cov = a.obs[t1_scope.LABEL_KEY].to_numpy() == "V-CM"
    assert np.array_equal(back[~cov], C[~cov]) and np.array_equal(base[~cov], C[~cov])
    assert np.array_equal(back[:, [1, 2]], C[:, [1, 2]])
    rel = np.abs(back[cov] - base[cov]) / base[cov]
    assert rel.mean() < 0.01 and rel.max() < 0.05, (rel.mean(), rel.max())
    # and not one edit at 0.4: the larger step's expectation differs
    once = np.mean([counts_of(pr.gc.moved_by_feng(a, dirs, 0.4, groups, held, [i, 2])).toarray()
                    for i in range(50)], axis=0)
    assert (np.abs(once[cov] - base[cov]) / base[cov]).mean() > 3 * rel.mean()


# ---------------------------------------------------------------------------------------------------------------------
# The reading


def _arm(score, mmd_u=0.1):
    return {"score": score, "groups": {"distribution": score / 100},
            "metrics": {"mmd_u": mmd_u, "variogram": 0.01},
            "skills": {"de_score": score / 200, "de_direction": 0.4, "mmd_u": 0.3, "variogram": 0.2}}


B0, BASE = 50.0, 53.0
CONTROLS = ("m1_r_gmom_x0", "m1_r_gmom_p1", "m1_r_gmom_p2", "m1_r_gmom_x2b")
FENG = {g: f"{i}" * 64 for i, g in enumerate(("CM", "Endocardium", "Endothelium", "Epicardium"))}
PERM = {g: f"{i + 4}" * 64 for i, g in enumerate(FENG)}


def _task(gains, seeds=range(10), commit="c" * 40, clean=True, task=None, fraction=0.2, fractions=None,
          identical=(), per_seed=None, base=BASE, placebo=None):
    """A task file. gains: arm -> gain over m1_r_gmom (every seed); per_seed: arm -> ten gains. A control's raw metric
    differs from the base's by 1e-9 unless named in `identical`."""
    rows = []
    for i, s in enumerate(seeds):
        arms = {"b0": _arm(B0), "m1_r_gmom": _arm(base)}
        for a in ("m1_r_gmom_x2", *CONTROLS):
            g = per_seed[a][i] if per_seed and a in per_seed else gains.get(a, 0.0)
            arms[a] = _arm(base + g, mmd_u=0.1 if a in identical else 0.1 + 1e-9)
        rows.append({"seed": s, "floor": 50.0, "arms": arms})
    fr = fractions or {"m1_r_gmom": fraction, "m1_r_gmom_x2": fraction * 2, "m1_r_gmom_x0": 0.0,
                       "m1_r_gmom_p1": fraction, "m1_r_gmom_p2": fraction * 2,
                       "m1_r_gmom_x2b": [fraction * 2, -fraction]}
    pl = placebo or {"m1_r_gmom_p1": dict(PERM), "m1_r_gmom_p2": dict(PERM)}
    return {"label": "synthetic", "task": task, "fraction": fraction, "fractions": fr, "seeds": list(seeds),
            "feng_sha256": dict(FENG), "placebo_sha256": pl, "rows": rows,
            "stamp": {"commit": commit, "tree_clean": clean, "utc": "x", "inputs": {}}}


def _t43_rows(redraw=3.8):
    return {"task": {"key": "x"}, "rows": [
        {"seed": s, "scores": {"b0": B0, "b0_mn": B0 + redraw},
         "skills": {"b0": {"de_score": 0.1, "mmd_u": 0.2}, "b0_mn": {"de_score": 0.15, "mmd_u": 0.25}}}
        for s in range(10)]}


def _write(pr, tmp_path, cuts, committed_shift=None):
    """cuts: cut key -> task record. T42's (b0, m1_r_gmom) and T44's (m1_r_gmom_x2) committed files copied from the
    cut's own scores, then shifted by committed_shift: (source, arm, amount) on the first seed of every cut; T43's
    beside them."""
    out = tmp_path / "out"
    (out / "tasks").mkdir(parents=True)
    dirs = {name: tmp_path / name for name in ("t42", "t44", "t43")}
    for d in dirs.values():
        d.mkdir()
    for k, f in pr.CUT_FILES.items():
        rec = cuts.get(k)
        if rec is not None:
            rec = {**rec, "task": k if rec.get("task") is None else rec["task"]}
            (out / "tasks" / f).write_text(json.dumps(rec))
        src = cuts.get(k) or _task({})
        for name, arms in (("t42", ("b0", "m1_r_gmom")), ("t44", ("b0", "m1_r_gmom", "m1_r_gmom_x2"))):
            rows = [{"seed": r["seed"], "arms": {a: copy.deepcopy(r["arms"][a]) for a in arms}} for r in src["rows"]]
            if committed_shift and committed_shift[0] == name:
                rows[0]["arms"][committed_shift[1]]["score"] += committed_shift[2]
            (dirs[name] / f).write_text(json.dumps({"rows": rows}))
        (dirs["t43"] / f).write_text(json.dumps(_t43_rows()))
    return out, dirs


def _read(pr, tmp_path, gains=None, per_seed=None, committed_shift=None, **over):
    gains = gains or {}
    cuts = {k: _task(gains, per_seed=per_seed) for k in pr.CUT_FILES}
    for k, kw in over.items():
        key = getattr(pr.gc, k)
        cuts[key] = None if kw is None else _task(gains, per_seed=per_seed, **kw)
    out, dirs = _write(pr, tmp_path, {k: v for k, v in cuts.items() if v is not None}, committed_shift)
    return pr.read(out, dirs["t42"], dirs["t44"], dirs["t43"])["reading"]


@pytest.mark.parametrize("gains, genes", [
    ({"m1_r_gmom_x2": 1.0}, "needs the direction's genes"),
    ({"m1_r_gmom_x2": 1.0, "m1_r_gmom_p2": 1.0}, "a placebo gains as well"),
    ({"m1_r_gmom_x2": 2.0, "m1_r_gmom_p2": 1.0}, "a placebo gains as well"),
    ({"m1_r_gmom_x2": 0.1}, "not separated"),
    ({}, "not separated"),
])
def test_each_gene_identity_statement_is_reached_by_a_table_built_for_it(pr, tmp_path, gains, genes):
    r = _read(pr, tmp_path, gains)
    assert r["valid"] and r["statements"]["gene identity"] == genes
    a = r["comparisons"]["a"]
    assert a["base"] == "m1_r_gmom_p2" and a["change"] == "m1_r_gmom_x2"
    want = gains.get("m1_r_gmom_x2", 0.0) - gains.get("m1_r_gmom_p2", 0.0)
    cut = a["cuts"][pr.gc.TASK_CHAMBER]
    assert cut["mean"] == pytest.approx(want) and cut["range"] == pytest.approx([want, want])
    assert set(cut["skills"]) == {"de_score", "de_direction", "mmd_u", "variogram"}
    assert cut["skills"]["de_score"] == pytest.approx(want / 200)


@pytest.mark.parametrize("gains, noise", [
    ({"m1_r_gmom_x2": 1.0, "m1_r_gmom_x2b": 1.0}, "is count noise"),
    ({"m1_r_gmom_x2": 1.0}, "is not count noise"),
    ({"m1_r_gmom_x2": 2.0, "m1_r_gmom_x2b": 1.0}, "not separated"),
    ({"m1_r_gmom_x2": 0.1, "m1_r_gmom_x2b": 0.05}, "not separated"),
])
def test_each_count_noise_statement_is_reached_by_a_table_built_for_it(pr, tmp_path, gains, noise):
    r = _read(pr, tmp_path, gains)
    assert r["valid"] and r["statements"]["count noise"] == noise
    both = r["comparisons"]["c"]["passes"] and r["comparisons"]["d"]["passes"]
    assert ("count_noise_note" in r) == both
    if both:
        assert "+1.00" in r["count_noise_note"]
    c, d = r["comparisons"]["c"], r["comparisons"]["d"]
    assert (c["base"], c["change"]) == ("m1_r_gmom_x2b", "m1_r_gmom_x2")
    assert (d["base"], d["change"]) == ("m1_r_gmom", "m1_r_gmom_x2b")
    # T43c's committed difference sits beside (b) and (d)
    assert r["t43c"]["cuts"][pr.gc.TASK_WHOLE]["mean"] == pytest.approx(3.8)
    assert r["t43c"]["base"] == "b0" and r["t43c"]["change"] == "b0_mn"


def test_the_seven_comparisons_and_the_reported_numbers(pr, tmp_path):
    gains = {"m1_r_gmom_x2": 1.5, "m1_r_gmom_x0": -0.7, "m1_r_gmom_p1": -0.4, "m1_r_gmom_p2": 0.2,
             "m1_r_gmom_x2b": 0.6}
    r = _read(pr, tmp_path, gains)
    want = {"a": 1.3, "b": 0.6, "c": 0.9, "d": 0.6, "e": 0.4, "f": 0.3, "g": 0.7}
    assert set(r["comparisons"]) == set(want)
    for c, w in want.items():
        assert r["comparisons"][c]["cuts"][pr.gc.TASK_WHOLE]["mean"] == pytest.approx(w)
    assert r["mean_scores"][pr.gc.TASK_CHAMBER]["m1_r_gmom_x0"] == pytest.approx(BASE - 0.7)
    out = tmp_path / "out" / "t45.txt"
    text = out.read_text()
    assert "(e) reported, decides nothing" in text and "try T43c's full count redraw" in text
    assert text.count("try T43c's full count redraw") == 2


def test_a_gain_of_exactly_point_three_does_not_pass(pr, tmp_path):
    # 0.5 - 0.2 is 0.3 exactly in floating point
    r = _read(pr, tmp_path, {"m1_r_gmom_x2": 0.3}, TASK_CHAMBER={"base": 0.2}, TASK_WHOLE={"base": 0.2})
    assert r["valid"]
    assert [row["arms"]["m1_r_gmom_x2"]["score"] - row["arms"]["m1_r_gmom_p2"]["score"]
            for row in [_task({"m1_r_gmom_x2": 0.3}, base=0.2)["rows"][0]]] == [0.3]
    assert not r["comparisons"]["a"]["passes"] and r["statements"]["gene identity"] == "not separated"
    r = _read(pr, tmp_path / "more", {"m1_r_gmom_x2": 0.31})
    assert r["comparisons"]["a"]["passes"] and r["statements"]["gene identity"] == "needs the direction's genes"


def test_a_loss_beyond_point_three_on_one_cut_vetoes(pr, tmp_path):
    cuts = {pr.gc.TASK_CHAMBER: _task({"m1_r_gmom_x2": 1.0}), pr.gc.TASK_WHOLE: _task({"m1_r_gmom_x2": -0.4})}
    out, dirs = _write(pr, tmp_path, cuts)
    r = pr.read(out, dirs["t42"], dirs["t44"], dirs["t43"])["reading"]
    assert r["valid"] and not r["comparisons"]["a"]["passes"]


@pytest.mark.parametrize("over, reason", [
    ({"TASK_WHOLE": None}, "missing cuts"),
    ({"TASK_CHAMBER": {"seeds": range(9)}}, "seeds"),
    ({"TASK_WHOLE": {"fraction": 0.4}}, "fraction 0.4"),
    ({"TASK_WHOLE": {"fractions": {"m1_r_gmom": 0.2, "m1_r_gmom_x2": 0.4, "m1_r_gmom_x0": 0.0, "m1_r_gmom_p1": 0.2,
                                   "m1_r_gmom_p2": 0.6, "m1_r_gmom_x2b": [0.4, -0.2]}}}, "wrong fraction"),
    ({"TASK_WHOLE": {"fractions": {"m1_r_gmom": 0.2, "m1_r_gmom_x2": 0.4, "m1_r_gmom_x0": 0.0, "m1_r_gmom_p1": 0.2,
                                   "m1_r_gmom_p2": 0.4, "m1_r_gmom_x2b": 0.2}}}, "wrong fraction"),
    ({"TASK_CHAMBER": {"commit": "d" * 40}}, "2 commits"),
    ({"TASK_WHOLE": {"clean": False}}, "not clean"),
    ({"TASK_CHAMBER": {"task": "T1/E9.5->E10.5"}}, "the file holds task"),
    ({"TASK_CHAMBER": {"identical": ("m1_r_gmom_x0",)}, "TASK_WHOLE": {"identical": ("m1_r_gmom_x0",)}},
     "controls score exactly as m1_r_gmom"),
    ({"TASK_WHOLE": {"placebo": {"m1_r_gmom_p1": {**PERM, "CM": FENG["CM"]},
                                 "m1_r_gmom_p2": {**PERM, "CM": FENG["CM"]}}}}, "equals Feng's direction"),
    ({"TASK_WHOLE": {"placebo": {"m1_r_gmom_p1": dict(PERM), "m1_r_gmom_p2": {**PERM, "CM": "f" * 64}}}},
     "do not record the same permuted directions"),
    ({"TASK_WHOLE": {"placebo": {"m1_r_gmom_p1": dict(PERM)}}}, "do not record the same permuted directions"),
])
def test_the_void_cases(pr, tmp_path, over, reason):
    r = _read(pr, tmp_path, {"m1_r_gmom_x2": 1.0, "m1_r_gmom_p1": 0.1}, **over)
    assert not r["valid"] and r["outcome"] == "VOID"
    assert any(reason in x for x in r["validity"]["reasons"]), r["validity"]["reasons"]


@pytest.mark.parametrize("source, arm", [("t42", "b0"), ("t42", "m1_r_gmom"), ("t44", "m1_r_gmom_x2")])
def test_a_rescored_arm_off_its_committed_score_by_one_e_minus_six_voids(pr, tmp_path, source, arm):
    r = _read(pr, tmp_path, {"m1_r_gmom_x2": 1.0}, committed_shift=(source, arm, 1e-6))
    assert not r["valid"] and any("do not reproduce" in x for x in r["validity"]["reasons"])
    assert r["validity"]["reproduction"][pr.gc.TASK_CHAMBER][source]["largest_difference"][arm] > 1e-9
    r = _read(pr, tmp_path / "ok", {"m1_r_gmom_x2": 1.0}, committed_shift=(source, arm, 1e-12))
    assert r["valid"]


def test_a_missing_t44_file_voids(pr, tmp_path):
    out, dirs = _write(pr, tmp_path, {k: _task({"m1_r_gmom_x2": 1.0}) for k in pr.CUT_FILES})
    (dirs["t44"] / pr.CUT_FILES[pr.gc.TASK_WHOLE]).unlink()
    r = pr.read(out, dirs["t42"], dirs["t44"], dirs["t43"])["reading"]
    assert not r["valid"] and any("T44" in x for x in r["validity"]["reasons"])


def test_a_control_identical_on_one_cut_only_stays_valid(pr, tmp_path):
    r = _read(pr, tmp_path, {"m1_r_gmom_x2": 1.0}, TASK_WHOLE={"identical": CONTROLS})
    assert r["valid"]


def test_the_reading_refuses_to_overwrite(pr, tmp_path):
    out, dirs = _write(pr, tmp_path, {k: _task({"m1_r_gmom_x2": 1.0}) for k in pr.CUT_FILES})
    pr.read(out, dirs["t42"], dirs["t44"], dirs["t43"])
    assert (out / "t45.json").exists() and (out / "t45.txt").exists()
    with pytest.raises(SystemExit, match="already holds"):
        pr.read(out, dirs["t42"], dirs["t44"], dirs["t43"])


def test_the_run_needs_its_row(pr, monkeypatch):
    monkeypatch.setattr(pr.n12, "tries_rows", lambda: {"T44", "T42"})
    monkeypatch.setattr(pr.gc, "n22_passed", lambda path: True)
    with pytest.raises(SystemExit, match="T45 has no row"):
        pr.main(["run", "--task", pr.gc.TASK_CHAMBER, "--out", "t45"])


def test_the_run_needs_n22(pr, monkeypatch):
    monkeypatch.setattr(pr.n12, "tries_rows", lambda: {"T45"})
    with pytest.raises(SystemExit, match="does not hold a passing N22"):
        pr.main(["run", "--task", pr.gc.TASK_CHAMBER, "--out", "t45", "--n22", "n22x"])


def test_the_run_refuses_an_existing_task_file(pr, monkeypatch):
    monkeypatch.setattr(pr.n12, "tries_rows", lambda: {"T45"})
    monkeypatch.setattr(pr.gc, "n22_passed", lambda path: True)
    tasks = pr.RESULTS / "t45" / "tasks"
    tasks.mkdir(parents=True)
    (tasks / pr.CUT_FILES[pr.gc.TASK_CHAMBER]).write_text("{}")
    with pytest.raises(SystemExit, match="already holds"):
        pr.main(["run", "--task", pr.gc.TASK_CHAMBER, "--out", "t45"])
