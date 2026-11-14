"""T1 batch t1-e, on synthetic data only: candidate b0_tail_filter, the tail decoder in the heart's form
(campaign/heart-b, try H17b) on B0's momentum factors (`vec.t1_scope.tail_infeasible_mask`,
`vec.t1_scope.tail_logf_filter` through `momentum_shift`'s logf_filter hook). With no negative
factor on a down-infeasible gene it is B0 bit for bit."""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse
from test_t1_batch_b import TASK, _build, _pair, _sha

from vec import harness
from vec import t1_scope as ts
from vec.prep import log1p_cp10k

B0_SHA = "23fab0fea7d206ed9d5f35050274d70e9c23ef5941888e869e8faf2a17a82508"  # b0 on _pair (test_t1_batch_b)
SEED = harness.BUILD_SEED_OFFSET + 1


# ---------------------------------------------------------------------------------------------------------------------
# 1. The mask and the filter

def test_the_mask_is_panel_genes_below_the_level():
    genes = [f"g{i}" for i in range(8)]
    panel = ["g5", "g0", "g1", "g2", "g3"]  # g4, g6, g7 are outside the panel
    pb_ref = np.array([0.01, 0.1, ts.TAIL_MIN_PB, 0.3, np.nextafter(ts.TAIL_MIN_PB, 0.0)])
    mask = ts.tail_infeasible_mask(pb_ref, panel, genes)
    assert mask.dtype == bool and mask.shape == (8,)
    assert np.flatnonzero(mask).tolist() == [0, 3, 5]  # g0 (0.1), g3 (just below), g5 (0.01)
    assert not mask[1] and not mask[2]  # at the level, above it
    assert not mask[[4, 6, 7]].any()  # outside the panel
    with pytest.raises(ValueError, match="shape"):
        ts.tail_infeasible_mask(pb_ref[:3], panel, genes)


def test_the_filter_sets_exactly_the_negative_infeasible_factors_to_zero():
    rng = np.random.default_rng(0)
    g = 200
    logf = rng.normal(0.0, 0.5, g)
    logf[:5] = 0.0  # exact zeros stay
    infeasible = rng.random(g) < 0.4
    held = np.zeros(g, dtype=bool)
    held[[7, 8, 9]] = True
    logf[held] = 0.0  # as momentum_shift hands them over
    hits: list = []
    out, rec = ts.tail_logf_filter(infeasible, hits)(logf.copy(), None, held)
    want = infeasible & (logf < 0)
    assert np.array_equal(out == 0.0, (logf == 0.0) | want)
    assert np.array_equal(out[~want], logf[~want])  # every other factor bit for bit, positives on infeasible included
    assert (out[infeasible & (logf > 0)] > 0).all()
    assert rec["n_compressed"] == int(want.sum()) > 0 and len(hits) == 1 and np.array_equal(hits[0], want)
    held_neg = logf.copy()
    held_neg[7] = -1.0
    infeasible[7] = True
    out2, _ = ts.tail_logf_filter(infeasible)(held_neg, None, held)
    assert out2[7] == -1.0  # a held gene is never the filter's (momentum_shift zeroes it anyway)
    none_hit, rec0 = ts.tail_logf_filter(np.zeros(g, dtype=bool))(logf.copy(), None, held)
    assert np.array_equal(none_hit, logf) and rec0["n_compressed"] == 0
    with pytest.raises(ValueError, match="shape"):
        ts.tail_logf_filter(infeasible)(logf[:10], None, held[:10])


# ---------------------------------------------------------------------------------------------------------------------
# 2. Through momentum_shift: only the intended factors change, held genes stay, the compressed fall is gone

GENES = [f"g{i}" for i in range(30)]


def _stage(n_per_type: int, mu_a: np.ndarray, mu_b: np.ndarray, seed: int, draw: bool = True) -> ad.AnnData:
    """Two types A and B; counts Poisson(mu) per cell, or exactly mu in every cell (draw False)."""
    rng = np.random.default_rng(seed)
    mu = np.vstack([np.tile(mu_a, (n_per_type, 1)), np.tile(mu_b, (n_per_type, 1))])
    C = sparse.csr_matrix((rng.poisson(mu) if draw else np.rint(mu)).astype(np.int64))
    types = np.r_[np.full(n_per_type, "A"), np.full(n_per_type, "B")]
    a = ad.AnnData(X=log1p_cp10k(C), layers={"counts": C.astype(np.int32)},
                   obs=pd.DataFrame({ts.CATEGORY_KEY: types, ts.LABEL_KEY: types},
                                    index=[f"s{seed}_{i}" for i in range(2 * n_per_type)]))
    a.var_names = GENES
    return a


def _momentum_case():
    """Type A: g0-g6 fall threefold from start to ref (so every other gene of A rises a little in CP10k); type B: g0
    rises threefold and g7 falls by the same counts, so B's other genes have exactly B0's factor 1 (start and ref are
    noise-free). Panel: every gene but g3. pb_ref (given, not read from the cells): g0-g3 and g5 below the level, g6
    exactly at it, the rest above. g5 is held. The cells to edit have g0-g6 well expressed, so the realised change is
    measured precisely."""
    base = np.full(30, 30.0)
    start_a, start_b = base.copy(), base.copy()
    start_a[:7] = 90.0
    start_b[0], start_b[7] = 10.0, 50.0
    start = _stage(300, start_a, start_b, 1, draw=False)
    ref = _stage(300, base, base, 2, draw=False)
    cells = _stage(1500, base, base, 3)
    panel = [g for g in GENES if g != "g3"]
    pb_ref = np.ones(len(panel))
    for g, v in (("g0", 0.1), ("g1", 0.1), ("g2", 0.2), ("g5", 0.05), ("g6", ts.TAIL_MIN_PB)):
        pb_ref[panel.index(g)] = v
    return start, ref, cells, panel, pb_ref


def _shift(monkeypatch, start, ref, cells, held, logf_filter=None):
    seen = []
    real = ts.apply_factors

    def spy(C, F, tidx, held_mask, **kw):
        res = real(C, F, tidx, held_mask, **kw)
        seen.append((F.copy(), tidx.copy(), res[1].copy()))
        return res

    monkeypatch.setattr(ts, "apply_factors", spy)
    out, rep = ts.momentum_shift(start, ref, cells, 0.5, held, min_cells=20, seed=[7, 1], logf_filter=logf_filter)
    monkeypatch.setattr(ts, "apply_factors", real)
    (F, tidx, scale), = seen
    return out, rep, F, tidx, scale


def test_through_momentum_shift_only_the_infeasible_falls_are_compressed(monkeypatch):
    start, ref, cells, panel, pb_ref = _momentum_case()
    held = ["g5"]
    mask = ts.tail_infeasible_mask(pb_ref, panel, GENES)
    assert np.flatnonzero(mask).tolist() == [0, 1, 2, 5]
    hits: list = []
    out0, rep0, F0, _, s0 = _shift(monkeypatch, start, ref, cells, held)
    out1, rep1, F1, tidx, s1 = _shift(monkeypatch, start, ref, cells, held, ts.tail_logf_filter(mask, hits))
    a, b = 0, 1  # sorted types
    assert (F0[a, :7][[0, 1, 2, 3, 4, 6]] < 0.8).all() and F0[a, 5] == 1.0 and F0[b, 0] > 1.2
    changed = F0 != F1
    assert np.flatnonzero(changed[a]).tolist() == [0, 1, 2] and not changed[b].any()
    assert np.all(F1[a, [0, 1, 2]] == 1.0)
    assert F1[b, 0] == F0[b, 0] > 1.0  # a rise on an infeasible gene is left
    assert F1[a, 3] == F0[a, 3] < 1.0  # outside the panel: untouched
    assert F1[a, 6] == F0[a, 6] < 1.0  # at the level: untouched
    assert F1[a, 4] == F0[a, 4] < 1.0  # feasible: untouched
    assert F1[a, 5] == F0[a, 5] == 1.0  # held
    assert [int(h.sum()) for h in hits] == [3, 0]
    assert [r["filter"]["n_compressed"] for r in rep1.rows] == [3, 0]
    C0, Cn = ts.counts_of(cells), ts.counts_of(out1)
    assert (Cn[:, 5] != C0[:, 5]).nnz == 0 and (ts.counts_of(out0)[:, 5] != C0[:, 5]).nnz == 0  # held untouched
    # The realised change of a compressed gene is the depth constant's alone. Its count in cell i is Binomial(c, s_i)
    # (s_i <= 1) or c + Poisson((s_i - 1) c), so its expected total over type A is sum s_i c_i, with variance at most
    # sum |1 - s_i| c_i. Tolerance: 4 such standard deviations below that expectation. Without the filter the same
    # gene's total sits near sum s_i f c_i, f about 0.58, far beyond the tolerance.
    in_a = tidx == a
    for j in (0, 1, 2):
        c = C0[:, j].toarray().ravel()[in_a].astype(float)
        exp1 = float((s1[in_a] * c).sum())
        tol = 4.0 * np.sqrt(float((np.abs(1.0 - s1[in_a]) * c).sum()) + 1.0)
        got1 = float(Cn[:, j].toarray().ravel()[in_a].sum())
        assert got1 >= exp1 - tol and abs(got1 - exp1) <= tol
        got0 = float(ts.counts_of(out0)[:, j].toarray().ravel()[in_a].sum())
        assert got0 < exp1 - 10 * tol  # b0's fall, which the filter removed
    assert np.array_equal(s0 == 1.0, s1 == 1.0)  # the same cells keep the constant at 1 (type B: nothing filtered)
    assert np.array_equal(s0[~in_a], s1[~in_a])


# ---------------------------------------------------------------------------------------------------------------------
# 3. The harness arm: B0 bit for bit with nothing to compress; the same cells, streams and held genes otherwise

def _check_lines(capsys) -> list[str]:
    return [x for x in capsys.readouterr().out.splitlines() if x.startswith("b0_tail_filter check:")]


def test_idle_arm_is_b0_bit_for_bit(monkeypatch, capsys):
    ref, start = _pair()
    a = _build("b0", ref, start)
    capsys.readouterr()
    b = _build("b0_tail_filter", ref, start)  # synthetic genes are all well expressed: nothing infeasible
    assert _sha(np.asarray(b.X, np.float32)) == _sha(np.asarray(a.X, np.float32)) == B0_SHA
    assert _check_lines(capsys) == [f"b0_tail_filter check: task={TASK.key} seed={SEED} infeasible=0 compressed=0 "
                                    "genes_compressed=0"]
    # every gene infeasible, but no negative factor: the stages are equal, so every log factor is exactly 0
    monkeypatch.setattr(ts, "TAIL_MIN_PB", 10.0)
    flat_a = _build("b0", ref, ref.copy())
    capsys.readouterr()
    flat_b = _build("b0_tail_filter", ref, ref.copy())
    assert _sha(np.asarray(flat_b.X, np.float32)) == _sha(np.asarray(flat_a.X, np.float32))
    assert _check_lines(capsys) == [f"b0_tail_filter check: task={TASK.key} seed={SEED} infeasible=50 compressed=0 "
                                    "genes_compressed=0"]


def _spy_shift(monkeypatch) -> list:
    seen = []
    real_shift, real_apply = ts.momentum_shift, ts.apply_factors

    def apply_spy(C, F, tidx, held_mask, **kw):
        seen[-1]["F"] = F.copy()
        return real_apply(C, F, tidx, held_mask, **kw)

    def shift_spy(*args, **kw):
        seen.append({"kw": dict(kw), "cells_in": args[2].copy()})
        out, rep = real_shift(*args, **kw)
        seen[-1]["out"] = out
        return out, rep

    monkeypatch.setattr(ts, "momentum_shift", shift_spy)
    monkeypatch.setattr(ts, "apply_factors", apply_spy)
    return seen


def test_active_arm_changes_only_the_filtered_factors(monkeypatch, capsys):
    ref, start = _pair()
    monkeypatch.setattr(ts, "TAIL_MIN_PB", 10.0)  # synthetic genes are all well expressed: make every gene infeasible
    seen = _spy_shift(monkeypatch)
    a = _build("b0", ref, start)
    capsys.readouterr()
    b = _build("b0_tail_filter", ref, start)
    lines = _check_lines(capsys)
    sa, sb = seen
    assert "logf_filter" not in sa["kw"] and "logf_filter" in sb["kw"]
    assert sa["kw"]["seed"] == sb["kw"]["seed"] == [SEED, 1]  # B0's stream, no new one
    assert {k: v for k, v in sb["kw"].items() if k != "logf_filter"} == sa["kw"]
    pd.testing.assert_frame_equal(sa["cells_in"].obs, sb["cells_in"].obs)  # the same drawn cells
    assert (ts.counts_of(sa["cells_in"]) != ts.counts_of(sb["cells_in"])).nnz == 0
    Fa, Fb = sa["F"], sb["F"]
    held = GENES_PAIR.index("g49")
    want = Fa < 1.0
    want[:, held] = False
    assert np.all(Fb[want] == 1.0) and np.array_equal(Fb[~want], Fa[~want])
    n_comp, n_genes = int(want.sum()), int(want.any(0).sum())
    assert n_comp > 0
    assert lines == [f"b0_tail_filter check: task={TASK.key} seed={SEED} infeasible=50 compressed={n_comp} "
                     f"genes_compressed={n_genes}"]
    assert (np.asarray(a.X) != np.asarray(b.X)).any()
    ca, cb = ts.counts_of(sa["out"]), ts.counts_of(sb["out"])
    assert (ca[:, held] != cb[:, held]).nnz == 0  # held: untouched in both
    lib_in = np.asarray(ts.counts_of(sb["cells_in"]).sum(1)).ravel()
    assert abs(np.asarray(cb.sum(1)).ravel().mean() / lib_in.mean() - 1) < 0.01  # depth kept, as in B0


GENES_PAIR = [f"g{i}" for i in range(50)]


def test_same_seed_same_output(monkeypatch, capsys):
    ref, start = _pair()
    monkeypatch.setattr(ts, "TAIL_MIN_PB", 10.0)
    b1 = _build("b0_tail_filter", ref, start)
    l1 = _check_lines(capsys)
    b2 = _build("b0_tail_filter", ref, start)
    l2 = _check_lines(capsys)
    assert _sha(np.asarray(b1.X, np.float32)) == _sha(np.asarray(b2.X, np.float32)) and l1 == l2 and len(l1) == 1


# ---------------------------------------------------------------------------------------------------------------------
# 4. The check line; registration; the other arms unchanged

def test_the_check_line_holds_only_the_allowed_fields():
    line = harness.tail_filter_check_line("T1/E9.5->E10.5", 1001,
                                          {"n_infeasible": 31, "n_compressed": 12, "n_genes_compressed": 7})
    assert line == "b0_tail_filter check: task=T1/E9.5->E10.5 seed=1001 infeasible=31 compressed=12 genes_compressed=7"
    keys = [t.split("=")[0] for t in line.split()[2:]]
    assert keys == ["task", "seed", "infeasible", "compressed", "genes_compressed"]


def test_registered_like_b0_tail():
    tasks = harness.load_tasks()["T1:test"]
    c, ref = harness.CANDIDATES["b0_tail_filter"], harness.CANDIDATES["b0_tail"]
    assert c.tasks == ref.tasks == frozenset({"T1"}) and c.derived == ref.derived == ("t1_held_flat",)
    assert c.externals == () and c.tuned_on == frozenset()
    assert [t.key for t in tasks if c.applies(t)] == [t.key for t in tasks if ref.applies(t)]
    for k in ("damping", "scope", "min_pb"):
        assert c.params[k] == ref.params[k]
    assert c.params["min_pb"] == ts.TAIL_MIN_PB


@pytest.mark.parametrize("name", ["b0", "b0_tail", "b0_tail_profile"])
def test_the_other_arms_are_unchanged(name, capsys):
    ref, start = _pair()
    assert _sha(np.asarray(_build(name, ref, start).X, np.float32)) == B0_SHA
