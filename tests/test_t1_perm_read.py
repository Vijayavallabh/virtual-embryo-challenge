"""experiments/t1_m1/perm_read.py (try T45 of experiments/harness/TRIES.md) on synthetic data: the control arms and the
reading. No released file, statistics file or score is read."""

from __future__ import annotations

import copy
import importlib.util
import json
import os
import re
import subprocess
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
PREVIOUS = "bfdebbd56c108eedddb1407d7df3ad5e6b37d76c"  # the commit before T45e and the second placebo pair


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


def test_the_second_placebo_is_t15s_permutation_at_the_next_seed(pr):
    dirs, halves = _dirs_and_halves()
    q = pr.placebo_directions(dirs, halves, pr.PLACEBO2_OFFSET)
    p = pr.placebo_directions(dirs, halves)
    assert pr.PLACEBO2_OFFSET == 1
    for g, d in dirs.items():
        ok = np.isfinite(d)
        want = d.copy()
        want[ok] = harness.t3_decile_permutation(d[ok], halves[g]["mean_cp"][ok],
                                                 [harness.LP_PERM_SEED + 1, harness.LP_GROUPS.index(g)],
                                                 harness.LP_BINS)
        assert np.array_equal(q[g], want, equal_nan=True)
    assert all(pr.direction_sha(q)[g] != pr.direction_sha(p)[g] != pr.direction_sha(dirs)[g] for g in dirs)
    # offset 0 is the first placebo, bit for bit
    assert pr.direction_sha(pr.placebo_directions(dirs, halves, 0)) == pr.direction_sha(p)


def test_the_between_bin_share_is_the_variance_of_t15s_bin_means_over_the_directions(pr):
    dirs, halves = _dirs_and_halves()
    share = pr.between_bin_share(dirs, halves)
    assert set(share) == set(dirs)
    for g, d in dirs.items():
        ok = np.isfinite(d)
        v, basis = d[ok], halves[g]["mean_cp"][ok]
        m = np.empty_like(v)
        for run in np.array_split(np.argsort(basis, kind="stable"), harness.LP_BINS):
            m[run] = v[run].mean()
        assert share[g] == pytest.approx(m.var() / v.var())
        assert 0 <= share[g] <= 1
        # a permutation within the bins keeps the bin means, so the placebo's share is the direction's
        pl = pr.placebo_directions({g: d}, {g: halves[g]})
        assert pr.between_bin_share(pl, {g: halves[g]})[g] == pytest.approx(share[g])
    # a direction that is its bin mean has share 1; a constant one has no variance
    basis = np.arange(40.0)
    step = np.repeat(np.arange(10.0), 4)
    assert pr.between_bin_share({"CM": step}, {"CM": {"mean_cp": basis}})["CM"] == pytest.approx(1.0)
    assert np.isnan(pr.between_bin_share({"CM": np.ones(40)}, {"CM": {"mean_cp": basis}})["CM"])


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
    dirs, placebo, placebo2 = {"CM": np.zeros(3)}, {"CM": np.ones(3)}, {"CM": np.full(3, 2.0)}

    class Ctx:
        def __init__(self, seed):
            self.seed = seed

    seen = {}
    for seed in (1000000, 1000003):
        calls.clear()
        ctx = Ctx(seed)
        built = pr.perm_arms(ctx, dirs, placebo, groups, 0.2, held_m1, held_f, placebo2=placebo2)
        arms = built["arms"]
        assert set(arms) == {"m1_r_gmom", *pr.CONTROLS, "m1_r_gmom_x2"}
        assert set(pr.CONTROLS) == {"m1_r_gmom_x0", "m1_r_gmom_p1", "m1_r_gmom_p2", "m1_r_gmom_x2b", "m1_r_gmom_x15b",
                                    "m1_r_gmom_q1", "m1_r_gmom_q2"}
        # the directions each arm was built with are the objects the build received (the placebo hashes come from these)
        assert all(built["directions"][a] is calls[i]["directions"] for i, a in enumerate(pr.M1_ARMS))
        assert built["directions"]["m1_r_gmom_p1"] is placebo and built["directions"]["m1_r_gmom_p2"] is placebo
        assert built["directions"]["m1_r_gmom_q1"] is placebo2 and built["directions"]["m1_r_gmom_q2"] is placebo2
        # x15b's first edit is at T44's x15 fraction, the same float step_read builds
        assert [c["fraction"] for c in calls] == [0.2, 0.4, 0.0, 0.2, 0.4, 0.4, 0.2 * 1.5, 0.2, 0.4]
        assert calls[6]["fraction"] == pr.sr.arm_fractions(0.2)["m1_r_gmom_x15"]
        want = {3: placebo, 4: placebo, 7: placebo2, 8: placebo2}
        assert all(c["directions"] is want.get(i, dirs) for i, c in enumerate(calls))
        for c in calls:
            for key in ("ctx", "groups", "held_m1", "held_factors"):
                assert c[key] is calls[0][key], f"an arm differs from the base build in {key}"
        assert arms["m1_r_gmom_p1"] == ("m1_r_gmom", 4) and arms["m1_r_gmom_p2"] == ("m1_r_gmom", 5)
        assert arms["m1_r_gmom_x2b"] == ("m1_r_gmom", 6) and arms["m1_r_gmom_x15b"] == ("m1_r_gmom", 7)
        assert arms["m1_r_gmom_q1"] == ("m1_r_gmom", 8) and arms["m1_r_gmom_q2"] == ("m1_r_gmom", 9)
        seen[seed] = [c["directions"] for c in calls[3:5]]
    # the same placebo object at both fractions and for two build seeds
    assert all(d is placebo for ds in seen.values() for d in ds)
    # the run's record hashes what each arm was built with: one record when every seed agrees, a voiding one if not
    sha = {a: {json.dumps(pr.direction_sha(placebo), sort_keys=True)} for a in pr.PLACEBO_ARMS}
    assert pr.placebo_record(sha) == {a: pr.direction_sha(placebo) for a in pr.PLACEBO_ARMS}
    sha["m1_r_gmom_p2"].add(json.dumps(pr.direction_sha(dirs), sort_keys=True))
    assert pr.DIFFERS in pr.placebo_record(sha)["m1_r_gmom_p2"]
    assert pr.arm_fractions(0.2) == {"m1_r_gmom": 0.2, "m1_r_gmom_x2": 0.4, "m1_r_gmom_x0": 0.0,
                                     "m1_r_gmom_p1": 0.2, "m1_r_gmom_p2": 0.4, "m1_r_gmom_x2b": [0.4, -0.2],
                                     "m1_r_gmom_x15b": [0.2 * 1.5, -0.1], "m1_r_gmom_q1": 0.2, "m1_r_gmom_q2": 0.4}
    assert pr.arm_fractions(0.2)["m1_r_gmom_x15b"] == [0.30000000000000004, -0.1]


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


def test_x15b_first_edit_is_t44s_x15_build_and_its_second_edit_goes_back_on_the_same_stream(pr, monkeypatch):
    """T45e: the first edit is bit for bit T44's m1_r_gmom_x15 build's (same fraction, input and stream), the second at
    (1 - 1.5) x 0.2 = -0.1 on the covered cells with R2's held genes on stream [seed, SECOND_STREAM], before g_mom."""
    ctx, groups, dirs, held = _arms_world()
    held_m1 = list(held)
    record, moved = [], []
    _spy_shift(monkeypatch, record)
    real_moved = pr.gc.moved_by_feng

    def spy_moved(pop, directions, fraction, groups_, held_, seed):
        out = real_moved(pop, directions, fraction, groups_, held_, seed)
        moved.append({"before": counts_of(pop).toarray(), "after": counts_of(out).toarray(), "seed": list(seed),
                      "fraction": fraction, "held": list(held_), "directions": directions})
        return out

    monkeypatch.setattr(pr.gc, "moved_by_feng", spy_moved)
    x15 = pr.sr.step_arms(ctx, dirs, groups, 0.2, held_m1, held)["m1_r_gmom_x15"]
    t44_first = [r for r in record if r["fraction"] == 0.2 * 1.5][-1]
    n_before = len(record)
    x15b = pr.there_and_back(ctx, dirs, groups, 0.2, held_m1, held, multiplier=pr.X15_MULT)
    new = record[n_before:]
    first = [r for r in new if r["fraction"] == 0.2 * 1.5]
    second = [r for r in new if r["fraction"] == -0.1]
    assert first and len(second) == 1 and len(moved) == 1
    assert all(np.array_equal(f["before"], t44_first["before"]) and np.array_equal(f["after"], t44_first["after"])
               and f["seed"] == t44_first["seed"] and f["seed"][-1] == 2 for f in first)
    m = moved[0]
    assert m["directions"] is dirs and m["held"] == held_m1 and m["fraction"] == -0.1
    assert m["seed"] == [ctx.seed, pr.SECOND_STREAM] == second[0]["seed"]
    assert pr.arm_fractions(0.2)["m1_r_gmom_x15b"] == [t44_first["fraction"], m["fraction"]]
    assert x15b.X.shape == x15.X.shape and not np.array_equal(x15b.X, x15.X)
    # the default multiplier is x2b's, bit for bit as before
    a = pr.there_and_back(ctx, dirs, groups, 0.2, held_m1, held)
    b = pr.there_and_back(ctx, dirs, groups, 0.2, held_m1, held, multiplier=2.0)
    assert a.X.tobytes() == b.X.tobytes()


_SKIP_DIRS = {"boards", "data", "runs", "notes", "research", "reference", ".venv", ".venv-gpu", ".git", ".claude",
              "results"}


def test_no_existing_step_uses_the_second_edits_stream(pr):
    """Every literal stream index written as 'seed, <n>]' in src/ and experiments/ (*.py, perm_read.py left out): the
    second edits' stream is none of them."""
    used = set()
    files = 0
    for top in ("src", "experiments"):
        for path in (REPO_ROOT / top).rglob("*.py"):
            rel = path.relative_to(REPO_ROOT)
            if path.name == "perm_read.py" or _SKIP_DIRS & set(rel.parts[:-1]):
                continue
            files += 1
            used |= {int(n) for n in re.findall(r"seed, ?(\d+)\]", path.read_text(errors="replace"))}
    assert files > 10 and {1, 2, 4, 5} <= used, (files, sorted(used))
    assert pr.SECOND_STREAM not in used, sorted(used)


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


def _sparse_toy(n_cells=40, n_genes=60, seed=11):
    """Poisson-gamma counts with mean about 0.4 (mostly zeros), as in a shallow single-cell matrix."""
    rng = np.random.default_rng(seed)
    C = rng.poisson(rng.gamma(0.5, 0.8, size=(n_cells, n_genes)))
    C[:, 3] += 1  # no empty cell
    labels = np.array(["V-CM"] * (n_cells - 6) + ["Foregut"] * 6)
    a = ad.AnnData(X=log1p_cp10k(sparse.csr_matrix(C)), obs={t1_scope.LABEL_KEY: labels})
    a.var_names = [f"g{i}" for i in range(n_genes)]
    a.obs_names = [f"c{i}" for i in range(n_cells)]
    a.layers["counts"] = sparse.csr_matrix(C.astype(np.int32))
    d = rng.normal(0, 0.6, n_genes)
    d[0] = np.nan
    return a, {"CM": d}, {"V-CM": "CM", "Foregut": "other:Foregut"}, ["g1", "g2"]


def test_there_and_backs_expected_counts_are_the_bases_on_sparse_counts(pr):
    """The same claim on mostly-zero counts, where zeros stay zero and the depth constants vary most: the mean over
    the draws of x2b's counts minus the base's is within four standard errors of 0, and the larger step's is not."""
    a, dirs, groups, held = _sparse_toy()
    C = counts_of(a).toarray()
    assert (C == 0).mean() > 0.5 and C.mean() < 0.7
    draws = 600
    base, back, once = [], [], []
    for i in range(draws):
        base.append(counts_of(pr.gc.moved_by_feng(a, dirs, 0.2, groups, held, [i, 2])).toarray())
        first = pr.gc.moved_by_feng(a, dirs, 0.4, groups, held, [i, 2])
        once.append(counts_of(first).toarray())
        back.append(counts_of(pr.second_edit(dirs, groups, -0.2, held, i)(first)).toarray())
    base, back, once = np.array(base, dtype=float), np.array(back, dtype=float), np.array(once, dtype=float)
    cov = a.obs[t1_scope.LABEL_KEY].to_numpy() == "V-CM"
    free = np.ones(C.shape[1], dtype=bool)
    free[[1, 2]] = False
    m = np.zeros(C.shape, dtype=bool)
    m[np.ix_(cov, free)] = True
    m &= C > 0
    assert np.array_equal(back[:, ~cov], np.broadcast_to(C[~cov], back[:, ~cov].shape))
    assert np.array_equal(back[:, :, [1, 2]], np.broadcast_to(C[:, [1, 2]], back[:, :, [1, 2]].shape))
    assert (back[:, C == 0] == 0).all(), "zeros stay zero"

    # the mean over entries of the draws' mean difference from the base, in standard errors of that mean
    diff = back[:, m].mean(axis=0) - base[:, m].mean(axis=0)
    se = np.sqrt((back[:, m].var(axis=0, ddof=1) + base[:, m].var(axis=0, ddof=1)).sum() / draws) / m.sum()
    assert abs(diff.mean() / se) < 4, diff.mean() / se
    # per entry, x2b's mean is the base's to the sampling error, and one edit at 0.4 is not
    dev_b = (back[:, m].mean(axis=0) - base[:, m].mean(axis=0))
    dev_o = (once[:, m].mean(axis=0) - base[:, m].mean(axis=0))
    se_e = np.sqrt((back[:, m].var(axis=0, ddof=1) + base[:, m].var(axis=0, ddof=1)) / draws)
    assert np.mean((dev_b / se_e) ** 2) < 1.5, "x2b's per-entry means match the base's within the sampling error"
    assert np.mean((dev_o / se_e) ** 2) > 3 * np.mean((dev_b / se_e) ** 2)
    # and x2b carries more count variance than the base: the excess noise the count-noise statement names
    assert back[:, m].var(axis=0).mean() > base[:, m].var(axis=0).mean()


def test_x15b_has_the_bases_expected_counts_and_about_x2s_zeroing_on_sparse_counts(pr):
    """T45e on mostly-zero counts: its mean counts are the base's within four standard errors, its share of non-zero
    counts set to zero is within 25% of x2's, and x2b's is above x2's by more than x15b's differs from x2's. (The
    reviewer's 1.5 times x2's for x2b holds on the reviewer's world, not on this toy: 0.095 against 0.076 here.)"""
    a, dirs, groups, held = _sparse_toy()
    C = counts_of(a).toarray()
    draws = 600
    base, x15b, x2, x2b = [], [], [], []
    for i in range(draws):
        base.append(counts_of(pr.gc.moved_by_feng(a, dirs, 0.2, groups, held, [i, 2])).toarray())
        f15 = pr.gc.moved_by_feng(a, dirs, 0.2 * 1.5, groups, held, [i, 2])
        x15b.append(counts_of(pr.second_edit(dirs, groups, (1.0 - 1.5) * 0.2, held, i)(f15)).toarray())
        f2 = pr.gc.moved_by_feng(a, dirs, 0.4, groups, held, [i, 2])
        x2.append(counts_of(f2).toarray())
        x2b.append(counts_of(pr.second_edit(dirs, groups, -0.2, held, i)(f2)).toarray())
    base, x15b, x2, x2b = (np.array(v, dtype=float) for v in (base, x15b, x2, x2b))
    cov = a.obs[t1_scope.LABEL_KEY].to_numpy() == "V-CM"
    free = np.ones(C.shape[1], dtype=bool)
    free[[1, 2]] = False
    m = np.zeros(C.shape, dtype=bool)
    m[np.ix_(cov, free)] = True
    m &= C > 0
    diff = x15b[:, m].mean(axis=0) - base[:, m].mean(axis=0)
    se = np.sqrt((x15b[:, m].var(axis=0, ddof=1) + base[:, m].var(axis=0, ddof=1)).sum() / draws) / m.sum()
    assert abs(diff.mean() / se) < 4, diff.mean() / se

    def zeroed(x):
        return float((x[:, m] == 0).mean())

    assert abs(zeroed(x15b) - zeroed(x2)) <= 0.25 * zeroed(x2), (zeroed(x15b), zeroed(x2))
    assert zeroed(x2b) - zeroed(x2) > abs(zeroed(x15b) - zeroed(x2)), (zeroed(x2b), zeroed(x15b), zeroed(x2))


def test_the_redraw_probe_passes_every_arm_through_bit_for_bit_and_measures_each_arms_redraw(pr):
    """perm_arms with its redraw probe builds every arm bit for bit as the same call without it, restores what it
    patched, and records target-free redraw measures: none at fraction 0, more at 0.4 than at 0.2."""
    ctx, groups, dirs, held = _arms_world()
    rng = np.random.default_rng(9)
    placebo = {g: rng.permutation(d) for g, d in dirs.items()}
    placebo2 = {g: rng.permutation(d) for g, d in dirs.items()}
    shift, hook = t1_scope.external_shift, pr.gc.gmom_hook
    built = pr.perm_arms(ctx, dirs, placebo, groups, 0.2, held, held, placebo2=placebo2)
    assert t1_scope.external_shift is shift and pr.gc.gmom_hook is hook
    plain = {"m1_r_gmom": (dirs, 0.2), "m1_r_gmom_x2": (dirs, 0.4), "m1_r_gmom_x0": (dirs, 0.0),
             "m1_r_gmom_p1": (placebo, 0.2), "m1_r_gmom_p2": (placebo, 0.4), "m1_r_gmom_q1": (placebo2, 0.2),
             "m1_r_gmom_q2": (placebo2, 0.4)}
    for a, (d, f) in plain.items():
        want = pr.gc.m1_r_arms(ctx, d, groups, f, held, held)["m1_r_gmom"]
        assert built["arms"][a].X.tobytes() == want.X.tobytes(), a
    back = pr.there_and_back(ctx, dirs, groups, 0.2, held, held)
    assert built["arms"]["m1_r_gmom_x2b"].X.tobytes() == back.X.tobytes()
    back15 = pr.there_and_back(ctx, dirs, groups, 0.2, held, held, multiplier=1.5)
    assert built["arms"]["m1_r_gmom_x15b"].X.tobytes() == back15.X.tobytes()
    rd = built["redraw"]
    assert set(rd) == set(pr.M1_ARMS) and all(set(v) == set(pr.REDRAW_KEYS) for v in rd.values())
    x0 = rd["m1_r_gmom_x0"]
    assert x0["share_nonzero_changed"] == 0.0 and x0["share_nonzero_zeroed"] == 0.0
    assert x0["detected_per_cell"] == x0["detected_per_cell_drawn"]
    assert len({v["detected_per_cell_drawn"] for v in rd.values()}) == 1, "every arm starts from the same drawn cells"
    base, x2, x2b = rd["m1_r_gmom"], rd["m1_r_gmom_x2"], rd["m1_r_gmom_x2b"]
    assert 0 < base["share_nonzero_changed"] < x2["share_nonzero_changed"]
    assert x2b["share_nonzero_changed"] > base["share_nonzero_changed"]
    assert x2b["share_nonzero_zeroed"] > base["share_nonzero_zeroed"]
    x15b = rd["m1_r_gmom_x15b"]
    assert x15b["share_nonzero_changed"] > base["share_nonzero_changed"]


def _previous_module(tmp_path, monkeypatch):
    """perm_read.py as the previous commit (PREVIOUS) has it, loaded beside the current one."""
    env = {**os.environ, "GIT_OPTIONAL_LOCKS": "0"}
    got = subprocess.run(["git", "-C", str(REPO_ROOT), "show", f"{PREVIOUS}:experiments/t1_m1/perm_read.py"],
                         capture_output=True, text=True, env=env)
    if got.returncode != 0:
        pytest.skip(f"the previous commit {PREVIOUS[:7]} is not in this checkout: {got.stderr.strip()}")
    path = tmp_path / "perm_read_previous.py"
    path.write_text(got.stdout)
    monkeypatch.syspath_prepend(str(SCRIPT.parent))
    spec = importlib.util.spec_from_file_location("perm_read_previous", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def test_every_existing_arm_is_bit_for_bit_the_previous_commits_build(pr, tmp_path, monkeypatch):
    """On the synthetic world: base, x2, x0, p1, p2 and x2b, and the rescored b0, as the previous commit builds them;
    the first placebo and the applied size unchanged as well."""
    old = _previous_module(tmp_path, monkeypatch)
    ctx, groups, dirs, held = _arms_world()
    rng = np.random.default_rng(9)
    placebo = {g: rng.permutation(d) for g, d in dirs.items()}
    placebo2 = {g: rng.permutation(d) for g, d in dirs.items()}
    was = old.perm_arms(ctx, dirs, placebo, groups, 0.2, held, held)
    now = pr.perm_arms(ctx, dirs, placebo, groups, 0.2, held, held, placebo2=placebo2)
    assert set(was["arms"]) == {"m1_r_gmom", "m1_r_gmom_x2", "m1_r_gmom_x0", "m1_r_gmom_p1", "m1_r_gmom_p2",
                                "m1_r_gmom_x2b"}
    for a, st in was["arms"].items():
        assert now["arms"][a].X.tobytes() == st.X.tobytes(), a
        assert now["redraw"][a] == was["redraw"][a], a
    b0 = harness._b0_t1(harness.B0_DAMPING)
    assert b0(ctx).X.tobytes() == b0(ctx).X.tobytes()  # b0 is harness code, which neither commit changes
    d, h = _dirs_and_halves()
    assert pr.direction_sha(pr.placebo_directions(d, h)) == old.direction_sha(old.placebo_directions(d, h))
    assert old.arm_fractions(0.2).items() <= pr.arm_fractions(0.2).items()


def test_redraw_measures_and_applied_size_on_a_hand_example(pr):
    B = np.array([[2, 0, 1, 3], [0, 4, 1, 0]])
    A = np.array([[2, 0, 0, 5], [0, 4, 1, 0]])
    m = pr.redraw_measures(sparse.csr_matrix(B), A)
    assert m == {"share_nonzero_changed": 2 / 5, "share_nonzero_zeroed": 1 / 5, "detected_per_cell": 2.0,
                 "detected_per_cell_drawn": 2.5}
    with pytest.raises(ValueError, match="not the drawn ones"):
        pr.redraw_measures(B, A[:1])
    d = {"CM": np.array([1.0, np.nan, -2.0, 4.0])}
    s = pr.applied_size(d, ["g3"], ["g0", "g1", "g2", "g3"], 0.5)
    assert s == {"CM": {"mean_abs_logf": pytest.approx((0.5 + 0 + 1.0) / 3), "sum_sq_logf": pytest.approx(1.25)}}
    sizes = {"fraction": 0.2, "feng": {"CM": {"mean_abs_logf": 1.0, "sum_sq_logf": 2.0}},
             "placebo": {"CM": {"mean_abs_logf": 1.05, "sum_sq_logf": 2.9}}}
    flags = pr.size_flags(sizes)
    assert len(flags) == 1 and "sum_sq_logf" in flags[0]
    assert pr.size_flags({**sizes, "placebo": {"CM": {"mean_abs_logf": 1.09, "sum_sq_logf": 1.82}}}) == []


# ---------------------------------------------------------------------------------------------------------------------
# The reading


def _arm(score, mmd_u=0.1):
    return {"score": score, "groups": {"distribution": score / 100},
            "metrics": {"mmd_u": mmd_u, "variogram": 0.01},
            "skills": {"de_score": score / 200, "de_direction": 0.4, "mmd_u": 0.3, "variogram": 0.2}}


B0, BASE = 50.0, 53.0
CONTROLS = ("m1_r_gmom_x0", "m1_r_gmom_p1", "m1_r_gmom_p2", "m1_r_gmom_x2b", "m1_r_gmom_x15b", "m1_r_gmom_q1",
            "m1_r_gmom_q2")
FENG = {g: f"{i}" * 64 for i, g in enumerate(("CM", "Endocardium", "Endothelium", "Epicardium"))}
PERM = {g: f"{i + 4}" * 64 for i, g in enumerate(FENG)}
PERM2 = {g: "ab"[i % 2] * 63 + f"{i}" for i, g in enumerate(FENG)}
REDRAW = {"m1_r_gmom": 0.3, "m1_r_gmom_x2": 0.4, "m1_r_gmom_x0": 0.0, "m1_r_gmom_p1": 0.3, "m1_r_gmom_p2": 0.4,
          "m1_r_gmom_x2b": 0.5, "m1_r_gmom_x15b": 0.42, "m1_r_gmom_q1": 0.3, "m1_r_gmom_q2": 0.4}
SHARE = {g: 0.25 + 0.1 * i for i, g in enumerate(FENG)}
SIZES = {"fraction": 0.2, "feng": {g: {"mean_abs_logf": 0.1, "sum_sq_logf": 2.0} for g in FENG},
         "placebo": {g: {"mean_abs_logf": 0.1, "sum_sq_logf": 2.9 if g == "CM" else 2.0} for g in FENG}}


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
        redraw = {a: {"share_nonzero_changed": REDRAW[a], "share_nonzero_zeroed": REDRAW[a] / 4,
                      "detected_per_cell": 128.2 - 30 * REDRAW[a], "detected_per_cell_drawn": 128.2} for a in REDRAW}
        rows.append({"seed": s, "floor": 50.0, "arms": arms, "redraw": redraw})
    fr = fractions or {"m1_r_gmom": fraction, "m1_r_gmom_x2": fraction * 2, "m1_r_gmom_x0": 0.0,
                       "m1_r_gmom_p1": fraction, "m1_r_gmom_p2": fraction * 2,
                       "m1_r_gmom_x2b": [fraction * 2, -fraction], "m1_r_gmom_x15b": [fraction * 1.5, -fraction / 2],
                       "m1_r_gmom_q1": fraction, "m1_r_gmom_q2": fraction * 2}
    pl = placebo or {"m1_r_gmom_p1": dict(PERM), "m1_r_gmom_p2": dict(PERM), "m1_r_gmom_q1": dict(PERM2),
                     "m1_r_gmom_q2": dict(PERM2)}
    return {"label": "synthetic", "task": task, "fraction": fraction, "fractions": fr, "seeds": list(seeds),
            "feng_sha256": dict(FENG), "placebo_sha256": pl, "rows": rows, "applied_size": copy.deepcopy(SIZES),
            "between_bin_share": dict(SHARE),
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
    ({"m1_r_gmom_x2": 2.0, "m1_r_gmom_p2": 1.0}, "a placebo gains as well; the direction more"),
    ({"m1_r_gmom_x2": 0.5, "m1_r_gmom_p2": 0.25}, "not separated"),
    ({"m1_r_gmom_x2": 1.0, "m1_r_gmom_p1": 0.5, "m1_r_gmom_p2": 1.4}, "a placebo gains as well"),
    ({"m1_r_gmom_x2": 0.1}, "no larger-step gain to explain"),
    ({}, "no larger-step gain to explain"),
])
def test_each_gene_identity_statement_is_reached_by_a_table_built_for_it(pr, tmp_path, gains, genes):
    r = _read(pr, tmp_path, gains)
    assert r["valid"] and r["statements"]["gene identity"] == genes
    a = r["comparisons"]["a"]
    assert a["base"] == "m1_r_gmom_p2 - m1_r_gmom_p1" and a["change"] == "m1_r_gmom_x2 - m1_r_gmom"
    x2, p1, p2 = (gains.get(k, 0.0) for k in ("m1_r_gmom_x2", "m1_r_gmom_p1", "m1_r_gmom_p2"))
    want = x2 - (p2 - p1)  # the larger step's increment against the placebo's, paired by seed
    cut = a["cuts"][pr.gc.TASK_CHAMBER]
    assert cut["mean"] == pytest.approx(want) and cut["range"] == pytest.approx([want, want])
    assert set(cut["skills"]) == {"de_score", "de_direction", "mmd_u", "variogram"}
    assert cut["skills"]["de_score"] == pytest.approx(want / 200)
    assert r["comparisons"]["a0"]["cuts"][pr.gc.TASK_CHAMBER]["mean"] == pytest.approx(x2 - p2)
    assert ("gene_identity_note" in r) == (genes == "a placebo gains as well; the direction more")


def test_a_is_the_increment_not_the_level(pr, tmp_path):
    """The reviewer's case: the entry's own step gains 1.0 over the placebo's at 0.2 on every seed (e), and the larger
    step's increment is the same for Feng's direction and the placebo (+0.5). The level x2 over p2 gains 1.0 and would
    have read 'needs the direction's genes'; the increment comparison does not pass."""
    inc_p = [1.0] * 6 + [0.5] + [-0.5] * 3  # p2 - p1: mean 0.5, 7 of 10 seeds ahead
    per_seed = {"m1_r_gmom_p1": [-1.0] * 10, "m1_r_gmom_p2": [-1.0 + x for x in inc_p],
                "m1_r_gmom_x2": [0.5] * 10}  # x2 - base = +0.5 every seed
    r = _read(pr, tmp_path, per_seed=per_seed)
    assert r["valid"]
    c = r["comparisons"]
    assert c["e"]["cuts"][pr.gc.TASK_WHOLE]["mean"] == pytest.approx(1.0)
    assert c["a0"]["passes"] and c["a0"]["cuts"][pr.gc.TASK_WHOLE]["mean"] == pytest.approx(1.0)
    assert c["a"]["cuts"][pr.gc.TASK_WHOLE]["mean"] == pytest.approx(0.0) and not c["a"]["passes"]
    assert not c["b"]["passes"] and c["h"]["passes"]
    assert r["statements"]["gene identity"] == "not separated"


def test_no_larger_step_gain_means_neither_statement_is_read(pr, tmp_path):
    """(a) and (c) pass but (h), x2 over m1_r_gmom, does not: both statements read 'no larger-step gain to explain'."""
    r = _read(pr, tmp_path, {"m1_r_gmom_x2": 0.2, "m1_r_gmom_p2": -0.5, "m1_r_gmom_x2b": -0.5})
    c = r["comparisons"]
    assert r["valid"] and c["a"]["passes"] and c["c"]["passes"] and not c["h"]["passes"]
    assert r["statements"] == {"gene identity": "no larger-step gain to explain",
                               "count noise": "no larger-step gain to explain"}
    assert (c["h"]["base"], c["h"]["change"]) == ("m1_r_gmom", "m1_r_gmom_x2")
    text = (tmp_path / "out" / "t45.txt").read_text()
    assert "the two statements are not read" in text and "(h) gates the two statements" in text


@pytest.mark.parametrize("gains, noise", [
    ({"m1_r_gmom_x2": 1.0, "m1_r_gmom_x15b": 1.0}, "count noise of x2b's size gains"),
    ({"m1_r_gmom_x2": 1.0}, "is not count noise"),
    ({"m1_r_gmom_x2": 2.0, "m1_r_gmom_x15b": 1.0}, "part of each"),
    ({"m1_r_gmom_x2": 0.5, "m1_r_gmom_x15b": 0.25}, "not separated"),
    ({"m1_r_gmom_x2": 0.1, "m1_r_gmom_x15b": 0.05}, "no larger-step gain to explain"),
    # (c) and (d), x2b's, decide nothing: their passing or not leaves the value to (c') and (d')
    ({"m1_r_gmom_x2": 1.0, "m1_r_gmom_x2b": 1.0}, "is not count noise"),
    ({"m1_r_gmom_x2": 1.0, "m1_r_gmom_x2b": 1.0, "m1_r_gmom_x15b": 1.0}, "count noise of x2b's size gains"),
    ({"m1_r_gmom_x2": 2.0, "m1_r_gmom_x2b": 0.0, "m1_r_gmom_x15b": 1.0}, "part of each"),
])
def test_each_count_noise_statement_is_reached_by_a_table_built_for_it(pr, tmp_path, gains, noise):
    r = _read(pr, tmp_path, gains)
    assert r["valid"] and r["statements"]["count noise"] == noise
    both = r["comparisons"]["c'"]["passes"] and r["comparisons"]["d'"]["passes"]
    assert ("count_noise_note" in r) == both
    if both:
        assert "+1.00" in r["count_noise_note"] and "(c')" in r["count_noise_note"]
    c, d = r["comparisons"]["c"], r["comparisons"]["d"]
    assert (c["base"], c["change"]) == ("m1_r_gmom_x2b", "m1_r_gmom_x2")
    assert (d["base"], d["change"]) == ("m1_r_gmom", "m1_r_gmom_x2b")
    c1, d1 = r["comparisons"]["c'"], r["comparisons"]["d'"]
    assert (c1["base"], c1["change"]) == ("m1_r_gmom_x15b", "m1_r_gmom_x2")
    assert (d1["base"], d1["change"]) == ("m1_r_gmom", "m1_r_gmom_x15b")
    # T43c's committed difference sits beside (b) and (d)
    assert r["t43c"]["cuts"][pr.gc.TASK_WHOLE]["mean"] == pytest.approx(3.8)
    assert r["t43c"]["base"] == "b0" and r["t43c"]["change"] == "b0_mn"


def test_the_comparisons_and_the_reported_numbers(pr, tmp_path):
    gains = {"m1_r_gmom_x2": 1.5, "m1_r_gmom_x0": -0.7, "m1_r_gmom_p1": -0.4, "m1_r_gmom_p2": 0.2,
             "m1_r_gmom_x2b": 0.6, "m1_r_gmom_x15b": 0.2, "m1_r_gmom_q1": -0.1, "m1_r_gmom_q2": 0.8}
    r = _read(pr, tmp_path, gains)
    want = {"a": 0.9, "a0": 1.3, "a2": 0.6, "b": 0.6, "b2": 0.9, "c": 0.9, "c'": 1.3, "d": 0.6, "d'": 0.2, "e": 0.4,
            "f": 0.3, "g": 0.7, "h": 1.5}
    assert set(r["comparisons"]) == set(want) == set(pr.COMPARISONS)
    assert len(pr.COMPARISONS) == 13  # the docstring's "used thirteen times on these files"
    assert "used thirteen times on these files" in " ".join(pr.__doc__.split())
    for c in pr.COMPARISONS:
        assert f"\n  ({c})" in pr.__doc__, f"the docstring lists ({c})"
    assert r["comparisons"]["a2"]["base"] == "m1_r_gmom_q2 - m1_r_gmom_q1"
    assert r["comparisons"]["a2"]["change"] == "m1_r_gmom_x2 - m1_r_gmom"
    for c, w in want.items():
        assert r["comparisons"][c]["cuts"][pr.gc.TASK_WHOLE]["mean"] == pytest.approx(w)
    assert r["mean_scores"][pr.gc.TASK_CHAMBER]["m1_r_gmom_x0"] == pytest.approx(BASE - 0.7)
    out = tmp_path / "out" / "t45.txt"
    text = out.read_text()
    for c in ("a0", "a2", "b2", "c", "d", "e", "f", "g"):
        assert f"({c}) reported, decides nothing" in text
    for c in ("a", "b", "c'", "d'"):
        assert f"({c}): " in text, c
    assert text.count("the upper bound of the noise reading") == 2
    assert text.count("this second permutation (LP_PERM_SEED + 1) shows the spread across permutations") == 2
    assert "(a): m1_r_gmom_x2 - m1_r_gmom over m1_r_gmom_p2 - m1_r_gmom_p1" in text
    assert text.count("try T43c's full count redraw") == 2
    # the redraw measures beside (c) and (d), the placebo's applied size (and its >10% flag) beside (a) and (b)
    assert text.count("redraw of the covered cells before g_mom") == 8
    assert "m1_r_gmom_x2b changed 0.5000, zeroed 0.1250" in text
    assert "m1_r_gmom_x15b changed 0.4200, zeroed 0.1050" in text
    assert text.count("differs from Feng's by more than 10%: CM sum_sq_logf") == 4
    assert text.count("the threshold decides nothing") == 4
    # every group's ratio placebo over Feng for both measures, and the between-bin share, beside (a) and (b)
    assert text.count("the placebo's applied size over Feng's") == 4
    assert "CM 1.000, 1.450; Endocardium 1.000, 1.000" in text
    assert r["size_ratios"][pr.gc.TASK_WHOLE]["CM"] == {"mean_abs_logf": pytest.approx(1.0),
                                                        "sum_sq_logf": pytest.approx(1.45)}
    assert text.count("between-bin share of the direction's variance") == 4
    assert "CM 0.250; Endocardium 0.350; Endothelium 0.450; Epicardium 0.550" in text
    assert r["between_bin_share"][pr.gc.TASK_WHOLE] == SHARE
    assert r["size_flags"][pr.gc.TASK_WHOLE] == ["CM sum_sq_logf: Feng 2, placebo 2.9"]
    assert r["redraw"][pr.gc.TASK_WHOLE]["m1_r_gmom_x0"]["share_nonzero_changed"] == 0.0


def test_a_gain_of_exactly_point_three_does_not_pass(pr, tmp_path):
    # 0.5 - 0.2 is 0.3 exactly in floating point
    r = _read(pr, tmp_path, {"m1_r_gmom_x2": 0.3}, TASK_CHAMBER={"base": 0.2}, TASK_WHOLE={"base": 0.2})
    assert r["valid"]
    assert [row["arms"]["m1_r_gmom_x2"]["score"] - row["arms"]["m1_r_gmom_p2"]["score"]
            for row in [_task({"m1_r_gmom_x2": 0.3}, base=0.2)["rows"][0]]] == [0.3]
    assert not r["comparisons"]["a"]["passes"] and not r["comparisons"]["h"]["passes"]
    assert r["statements"]["gene identity"] == "no larger-step gain to explain"
    # (a) at exactly 0.3 with (h) passing: x2 gains 0.5 over the base, the placebo's increment 0.2
    r = _read(pr, tmp_path / "inc", {"m1_r_gmom_x2": 0.5, "m1_r_gmom_p2": 0.2}, TASK_CHAMBER={"base": 0.0},
              TASK_WHOLE={"base": 0.0})
    assert [0.5 - 0.0 - (0.2 - 0.0)] == [0.3]
    assert r["comparisons"]["h"]["passes"] and not r["comparisons"]["a"]["passes"]
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
    ({"TASK_WHOLE": {"placebo": {"m1_r_gmom_p1": dict(PERM), "m1_r_gmom_p2": {"differs between seeds": ["x", "y"]}}}},
     "different directions on different seeds"),
    ({"TASK_WHOLE": {"fractions": {"m1_r_gmom": 0.2, "m1_r_gmom_x2": 0.4, "m1_r_gmom_x0": 0.0, "m1_r_gmom_p1": 0.2,
                                   "m1_r_gmom_p2": 0.4, "m1_r_gmom_x2b": [0.4, -0.2], "m1_r_gmom_x15b": [0.3, -0.2],
                                   "m1_r_gmom_q1": 0.2, "m1_r_gmom_q2": 0.4}}}, "wrong fraction on ['m1_r_gmom_x15b']"),
    ({"TASK_WHOLE": {"fractions": {"m1_r_gmom": 0.2, "m1_r_gmom_x2": 0.4, "m1_r_gmom_x0": 0.0, "m1_r_gmom_p1": 0.2,
                                   "m1_r_gmom_p2": 0.4, "m1_r_gmom_x2b": [0.4, -0.2], "m1_r_gmom_x15b": [0.3, -0.1],
                                   "m1_r_gmom_q1": 0.2}}}, "wrong fraction on ['m1_r_gmom_q2']"),
    ({"TASK_CHAMBER": {"identical": ("m1_r_gmom_x15b",)}, "TASK_WHOLE": {"identical": ("m1_r_gmom_x15b",)}},
     "controls score exactly as m1_r_gmom on every seed of both cuts: ['m1_r_gmom_x15b']"),
    ({"TASK_CHAMBER": {"identical": ("m1_r_gmom_q1",)}, "TASK_WHOLE": {"identical": ("m1_r_gmom_q1",)}},
     "controls score exactly as m1_r_gmom on every seed of both cuts: ['m1_r_gmom_q1']"),
    ({"TASK_WHOLE": {"placebo": {"m1_r_gmom_p1": dict(PERM), "m1_r_gmom_p2": dict(PERM),
                                 "m1_r_gmom_q1": {**PERM2, "Endothelium": PERM["Endothelium"]},
                                 "m1_r_gmom_q2": {**PERM2, "Endothelium": PERM["Endothelium"]}}}},
     "the second placebo equals the first placebo in groups ['Endothelium']"),
    ({"TASK_WHOLE": {"placebo": {"m1_r_gmom_p1": dict(PERM), "m1_r_gmom_p2": dict(PERM),
                                 "m1_r_gmom_q1": {**PERM2, "CM": FENG["CM"]},
                                 "m1_r_gmom_q2": {**PERM2, "CM": FENG["CM"]}}}},
     "the second placebo equals Feng's direction in groups ['CM']"),
    ({"TASK_WHOLE": {"placebo": {"m1_r_gmom_p1": dict(PERM), "m1_r_gmom_p2": dict(PERM), "m1_r_gmom_q1": dict(PERM2),
                                 "m1_r_gmom_q2": {**PERM2, "CM": "f" * 64}}}},
     "m1_r_gmom_q1 and m1_r_gmom_q2 do not record the same permuted directions"),
    ({"TASK_WHOLE": {"placebo": {"m1_r_gmom_p1": dict(PERM), "m1_r_gmom_p2": dict(PERM), "m1_r_gmom_q1": dict(PERM2),
                                 "m1_r_gmom_q2": {"differs between seeds": ["x", "y"]}}}},
     "a placebo arm of the second placebo was built with different directions on different seeds"),
    ({"TASK_WHOLE": {"placebo": {"m1_r_gmom_p1": dict(PERM), "m1_r_gmom_p2": dict(PERM)}}},
     "m1_r_gmom_q1 and m1_r_gmom_q2 do not record the same permuted directions"),
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


def _t44_reading(tmp_path, valid=True, passing=("m1_r_gmom_x2",)):
    """A committed T44 reading in step_read.t44_reading's form: comparisons maps each scaled arm to its record."""
    path = tmp_path / "t44.json"
    comps = {a: {"base": "m1_r_gmom", "change": a, "passes": a in passing, "cuts": {}}
             for a in ("m1_r_gmom_x15", "m1_r_gmom_x2")}
    reading = {"valid": valid, "outcome": "x", **({"comparisons": comps} if valid else {})}
    path.write_text(json.dumps({"label": "synthetic", "reading": reading}))
    return path


@pytest.mark.parametrize("make, ok", [
    (lambda tp: _t44_reading(tp), True),
    (lambda tp: _t44_reading(tp, passing=("m1_r_gmom_x2", "m1_r_gmom_x15")), True),
    (lambda tp: _t44_reading(tp, passing=("m1_r_gmom_x15",)), False),
    (lambda tp: _t44_reading(tp, passing=()), False),
    (lambda tp: _t44_reading(tp, valid=False), False),
    (lambda tp: tp / "absent.json", False),
])
def test_the_run_needs_t44s_reading_to_name_x2_passing(pr, monkeypatch, tmp_path, make, ok):
    monkeypatch.setattr(pr.n12, "tries_rows", lambda: {"T45"})
    monkeypatch.setattr(pr.gc, "n22_passed", lambda path: True)
    path = make(tmp_path)
    monkeypatch.setattr(pr, "T44_READING", path)
    assert pr.t44_names_x2(path) is ok
    tasks = pr.RESULTS / "t45" / "tasks"
    tasks.mkdir(parents=True)
    (tasks / pr.CUT_FILES[pr.gc.TASK_CHAMBER]).write_text("{}")
    # a passing reading lets the run reach its next refusal (the existing task file); otherwise it stops here
    with pytest.raises(SystemExit, match="already holds" if ok else "naming m1_r_gmom_x2 among the arms that pass"):
        pr.main(["run", "--task", pr.gc.TASK_CHAMBER, "--out", "t45"])


def test_the_t44_reading_path_is_t44s_committed_reading(pr):
    assert pr.T44_READING == REPO_ROOT / "experiments" / "t1_m1" / "results" / "t44_m1_step" / "t44.json"


def test_the_run_refuses_an_existing_task_file(pr, monkeypatch, tmp_path):
    monkeypatch.setattr(pr.n12, "tries_rows", lambda: {"T45"})
    monkeypatch.setattr(pr.gc, "n22_passed", lambda path: True)
    monkeypatch.setattr(pr, "T44_READING", _t44_reading(tmp_path))
    tasks = pr.RESULTS / "t45" / "tasks"
    tasks.mkdir(parents=True)
    (tasks / pr.CUT_FILES[pr.gc.TASK_CHAMBER]).write_text("{}")
    with pytest.raises(SystemExit, match="already holds"):
        pr.main(["run", "--task", pr.gc.TASK_CHAMBER, "--out", "t45"])
