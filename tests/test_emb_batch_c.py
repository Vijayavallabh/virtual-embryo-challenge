"""Embryo batch emb-c (campaign branch): candidates fgw_tail (the embryo tail decoder) and fgw_geoshares
(exponential-growth state shares). Synthetic data only: arrays and AnnData objects built here."""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from vec import harness
from vec import t2_build as tb
from vec.boards import Board
from vec.harness import CANDIDATES, PseudoTask
from vec.io import COORD_KEY

N_PANEL = 40
N_EXTRA = 20
KW = dict(target_rms=1.0, n_anchor=150, k=4)
STREAM, GSTREAM = harness.EMB_TAIL_STREAM, harness.GEOSHARES_STREAM


CENTRES = {"A": (0.0, 0.0, 0.0), "epi": (0.0, 0.0, 0.0), "meso": (6.0, 0.0, 0.0), "exe": (0.0, 6.0, 0.0)}


def _counted(rng, n, depth, profile, labels=None, scale=1.0):
    """An AnnData as the released files store it: Poisson counts on the panel plus genes outside it, X = log(1 + 1e4 c
    / N) with N over every gene, counts in layers["counts"], 3D coordinates and obs["group"] / obs["celltype"]."""
    prof = np.asarray(profile, float) / np.sum(profile)
    C = rng.poisson(prof * depth, size=(n, len(prof)))
    C[C[:, :N_PANEL].sum(1) == 0, 0] += 1  # every cell has a panel count, as every released cell does
    N = C.sum(1)
    genes = [f"g{i}" for i in range(len(prof))]
    lab = np.array(["A"] * n) if labels is None else np.asarray(labels)
    a = ad.AnnData(X=sparse.csr_matrix(np.log1p(1e4 * C / N[:, None]).astype(np.float32)),
                   obs=pd.DataFrame({"group": lab, "celltype": lab}, index=[f"c{i}" for i in range(n)]),
                   var=pd.DataFrame(index=genes))
    a.layers["counts"] = sparse.csr_matrix(C.astype(np.int32))
    centres = np.array([CENTRES[k] for k in lab]) * scale  # clustered by label, as test_harness's synthetic embryo
    a.obsm[COORD_KEY] = (centres + rng.normal(size=(n, 3)) * 1.5 * scale).astype(np.float32)
    return a, genes[:N_PANEL]


def _profile(rng, rare=(), zero=()):
    p = rng.gamma(2.0, 1.0, size=N_PANEL + N_EXTRA)
    for g in rare:
        p[g] = 0.0
    for g in zero:
        p[g] = 0.0
    p = p / p.sum()
    for g in rare:  # about 0.03 counts per cell at depth 120: detected in about 3% of cells, mean X well below 0.25
        p[g] = 0.03 / 120.0
    return p


def _pair(seed=0, n=1200, n_draw=1000, t=0.5, rare=(), zero=()):
    rng = np.random.default_rng(seed)
    prof = _profile(rng, rare, zero)
    before, panel = _counted(rng, n, 120.0, prof)
    after, _ = _counted(rng, n, 120.0, prof)
    Eb, Ea = (np.asarray(x[:, panel].X.toarray(), np.float32) for x in (before, after))
    res = tb.interpolate_clouds(before.obsm[COORD_KEY], after.obsm[COORD_KEY], Eb, Ea, t, n_draw, seed=3, **KW)
    X = np.empty((len(res.rows), len(panel)), np.float32)
    X[res.stage == 0], X[res.stage == 1] = Eb[res.rows[res.stage == 0]], Ea[res.rows[res.stage == 1]]
    return before, after, panel, res, X


# ------------------------------------------------------------------------------------------------ fgw_tail
def test_tail_selection_is_every_gene_below_the_level_with_a_negative_change():
    before, after, panel, res, X = _pair(seed=1, rare=(1, 2, 3, 4, 5, 6))
    m = X.astype(np.float64).mean(0)
    assert (m[1:7] < 0.2).all() and (m[7:9] > 0.25).all()
    pb_ref = np.full(len(panel), 1.0)
    # genes 1-4: below the level, change negative -> selected; 5-6: below the level, change positive or zero -> not;
    # 7-8: change negative but the level at or above 0.25 -> not (no rank cut: all qualifying genes, not a bottom few)
    pb_ref[1:5] = m[1:5] + np.array([0.001, 0.005, 0.01, 0.02])
    pb_ref[5] = m[5] / 2
    pb_ref[6] = m[6]
    pb_ref[7], pb_ref[8] = max(m[7] + 0.01, 0.25), m[8] + 0.5
    expect = (pb_ref < 0.25) & (m - pb_ref < 0)
    _, rep = tb.tail_decoded_expression(before, after, res, panel, X, pb_ref, seed=[0, STREAM], report_bottom=80)
    assert rep["n_selected"] == int(expect.sum())
    assert expect[1:5].all() and not expect[5:9].any()
    assert rep["n_lifted"] + rep["n_skipped"] == rep["n_selected"]


def test_a_liftable_gene_lands_near_flat_and_out_of_reach_genes_are_skipped_and_counted():
    """Gene 3 (rare) starts at dp = -0.01, needs a factor below 2 and is lifted. Tolerance: over eight decoder seeds its
    realised dp averages within 0.003 of flat and each lies within 0.008 of it (the first-order solve's bias and the
    Poisson top-up's noise on about 1,000 drawn cells). Gene 4 (rare) needs far more than 2 and gene 5 has no counts in
    any cell: both are skipped and counted."""
    before, after, panel, res, X = _pair(seed=2, rare=(3, 4), zero=(5,))
    m = X.astype(np.float64).mean(0)
    assert m[5] == 0.0 and 0.0 < m[3] < 0.2 and 0.0 < m[4] < 0.2
    pb_ref = np.maximum(m + 0.5, 0.3)  # everything else at or above the level: not selected
    pb_ref[3] = m[3] + 0.01
    pb_ref[4] = min(m[4] + 0.15, 0.2499)
    pb_ref[5] = 0.1
    lifted_dp = []
    for s in range(8):
        out, rep = tb.tail_decoded_expression(before, after, res, panel, X, pb_ref, seed=[s, STREAM], report_bottom=80)
        assert rep["n_selected"] == 3 and rep["n_lifted"] == 1 and rep["n_skipped"] == 2
        assert 1.0 < rep["max_factor"] <= tb.EMB_TAIL_F_MAX and rep["median_factor"] == rep["max_factor"]
        lifted_dp.append(out[:, 3].astype(np.float64).mean() - pb_ref[3])
    assert abs(np.mean(lifted_dp)) < 0.003 and np.abs(lifted_dp).max() < 0.008
    dp_after = out.astype(np.float64).mean(0) - pb_ref
    assert dp_after[4] == pytest.approx(m[4] - pb_ref[4], abs=0.01)  # skipped: moves only by the depth constant
    assert (out[:, 5] == 0).all()
    assert out.dtype == np.float32 and out.shape == X.shape and np.isfinite(out).all()


def test_tail_keeps_each_cells_expected_panel_total():
    before, after, panel, res, X = _pair(seed=3, rare=(3, 6, 9))
    m = X.astype(np.float64).mean(0)
    pb_ref = np.maximum(m + 0.5, 0.3)
    pb_ref[[3, 6, 9]] = m[[3, 6, 9]] + 0.01
    _, rep = tb.tail_decoded_expression(before, after, res, panel, X, pb_ref, seed=[0, STREAM], report_bottom=80)
    assert rep["n_lifted"] == 3
    assert rep["expected_panel_max_error"] < 1e-9  # each cell's expected panel total is exactly its own
    assert rep["depth_ratio_mean"] == pytest.approx(1.0, abs=0.01)  # realised: Poisson and binomial noise only


def test_tail_with_nothing_selected_or_nothing_liftable_returns_the_expression_as_given():
    before, after, panel, res, X = _pair(seed=4, zero=(5,))
    m = X.astype(np.float64).mean(0)
    out, rep = tb.tail_decoded_expression(before, after, res, panel, X, np.maximum(m + 0.5, 0.3), seed=0,
                                          report_bottom=80)
    assert out is X and rep["n_selected"] == 0 and rep["n_lifted"] == 0 and np.isnan(rep["max_factor"])
    pb_ref = np.maximum(m + 0.5, 0.3)
    pb_ref[5] = 0.1  # selected, but no counts: nothing liftable
    out, rep = tb.tail_decoded_expression(before, after, res, panel, X, pb_ref, seed=0, report_bottom=80)
    assert out is X and rep["n_selected"] == 1 and rep["n_skipped"] == 1


def test_tail_is_deterministic_for_a_seed():
    before, after, panel, res, X = _pair(seed=5, rare=(3,))
    m = X.astype(np.float64).mean(0)
    pb_ref = np.maximum(m + 0.5, 0.3)
    pb_ref[3] = m[3] + 0.01
    a, ra = tb.tail_decoded_expression(before, after, res, panel, X, pb_ref, seed=[7, STREAM], report_bottom=80)
    b, rb = tb.tail_decoded_expression(before, after, res, panel, X, pb_ref, seed=[7, STREAM], report_bottom=80)
    np.testing.assert_array_equal(a, b)
    assert ra == rb or all(ra[k] == rb[k] or (ra[k] != ra[k] and rb[k] != rb[k]) for k in ra)


def test_tail_reads_the_drawn_cells_counts_as_the_renders_do():
    """With nothing lifted the re-rendering is not run; when it is, a cell whose counts are unchanged comes back as
    log(1 + 1e4 c / N) of its released counts and library, which is the released X up to float32 rounding."""
    before, after, panel, res, X = _pair(seed=6)
    C, lib, other = tb._drawn_counts(before, after, res, panel)
    np.testing.assert_allclose(tb._log_normalise(C, lib), X, atol=1e-5)
    np.testing.assert_array_equal(np.asarray(C.sum(1)).ravel() + other, lib)


# ------------------------------------------------------------------------------------------------ harness level
def _task(key):
    return PseudoTask(key=key, test_board="T2:heart:test_extrap", task="T2", setting="embryo", kind="interp",
                      target="T2/E7.5.h5ad", reference="T2/E7.25.h5ad",
                      inputs={"before": "T2/E7.25.h5ad", "after": "T2/E8.0.h5ad"}, twin_board="T2:embryo:val_interp")


def _ctx(seed=0, rare=False, shares=((0.5, 0.3, 0.2), (0.2, 0.5, 0.3))):
    rng = np.random.default_rng(100 + seed)
    pools = {}
    prof = _profile(rng)
    for role, depth, p, sc in (("before", 120.0, shares[0], 1.0), ("after", 150.0, shares[1], 1.8)):
        lab = rng.choice(["epi", "meso", "exe"], size=900, p=p)
        pr = prof.copy()
        if rare:  # genes 2-9 at a low level in both brackets: the draw puts some below the reference, by chance
            pr[2:10] = 0.03 / depth
        a, panel = _counted(rng, 900, depth, pr, labels=lab, scale=sc)
        pools[role] = a
    board = Board(key="pseudo:synthetic_embryo", task="T2", genes=panel, needs_coords=True, setting="embryo",
                  min_cells=100, max_cells=700)
    task = _task(f"syn_emb_c_{seed}")
    return harness.BuildContext(task=harness.TaskView.of(task), board=board, inputs=pools, seed=seed), task


def test_fgw_tail_idle_equals_fgw_bit_for_bit_and_prints_its_check_line(capsys):
    ctx, task = _ctx(seed=1)
    assert CANDIDATES["fgw_tail"].applies(task) and CANDIDATES["fgw_geoshares"].applies(task)
    pb = harness.reference_pb(ctx)
    assert (pb >= 0.25).all()  # every gene well expressed: nothing can be selected
    fgw = CANDIDATES["fgw"].build(ctx)
    capsys.readouterr()
    v = CANDIDATES["fgw_tail"].build(ctx)
    np.testing.assert_array_equal(v.X, fgw.X)
    np.testing.assert_array_equal(v.coords, fgw.coords)
    line = capsys.readouterr().out.strip()
    assert line == (f"fgw_tail check: task={task.key} seed=1 infeasible=0 lifted=0 skipped=0 in_bottom80=0 "
                    "max_factor=nan median_factor=nan")


def test_fgw_tail_lifts_on_fgws_cells_and_positions(capsys):
    ctx, task = _ctx(seed=2, rare=True)
    fgw = CANDIDATES["fgw"].build(ctx)
    capsys.readouterr()
    v = CANDIDATES["fgw_tail"].build(ctx)
    line = capsys.readouterr().out.strip()
    fields = dict(kv.split("=") for kv in line.split(": ", 1)[1].split())
    assert set(fields) == {"task", "seed", "infeasible", "lifted", "skipped", "in_bottom80", "max_factor",
                           "median_factor"}
    assert fields["task"] == task.key and fields["seed"] == "2"
    assert int(fields["infeasible"]) >= 2 and int(fields["lifted"]) >= 1
    assert int(fields["lifted"]) + int(fields["skipped"]) == int(fields["infeasible"])
    assert int(fields["in_bottom80"]) == int(fields["infeasible"])  # 40 panel genes: every gene is in the bottom 80
    np.testing.assert_array_equal(v.coords, fgw.coords)  # the same cells at the same positions
    assert not np.array_equal(v.X, fgw.X)
    again = CANDIDATES["fgw_tail"].build(ctx)
    np.testing.assert_array_equal(v.X, again.X)


# ------------------------------------------------------------------------------------------------ fgw_geoshares
def test_geometric_shares_keep_one_bracket_labels_linear_and_sum_to_one():
    pb = np.array([0.5, 0.3, 0.2, 0.0])
    pa = np.array([0.2, 0.5, 0.0, 0.3])
    t = 0.4
    q = tb.geometric_shares(pb, pa, t)
    lin = (1 - t) * pb + t * pa
    assert q.sum() == pytest.approx(1.0)
    assert q[2] == pytest.approx(lin[2]) and q[3] == pytest.approx(lin[3])  # present in one bracket: linear share
    g = pb[:2] ** (1 - t) * pa[:2] ** t
    np.testing.assert_allclose(q[:2], g / g.sum() * (1 - lin[2] - lin[3]))
    np.testing.assert_allclose(tb.geometric_shares(pb, pb, t), pb)  # identical brackets: the linear shares


def test_ipf_table_margins_and_rounded_counts():
    pb = np.array([0.5, 0.3, 0.2, 0.0])
    pa = np.array([0.2, 0.5, 0.0, 0.3])
    t, n = 0.4, 997
    n_b = int(round((1 - t) * n))
    n_a = n - n_b
    q = tb.geometric_shares(pb, pa, t)
    table, rec = tb.ipf_table(np.vstack([(1 - t) * pb, t * pa]), np.array([n_b, n_a]) / n, q)
    np.testing.assert_allclose(table.sum(1), [n_b / n, n_a / n], atol=1e-12)
    np.testing.assert_allclose(table.sum(0), q, atol=1e-12)
    assert table[0, 3] == 0.0 and table[1, 2] == 0.0  # zero cells stay zero
    assert rec["max_margin_error"] < 1e-12 and rec["sweeps"] >= 1
    Q = tb._geoshare_quotas(table, n_b, n_a)
    assert Q[0].sum() == n_b and Q[1].sum() == n_a  # both bracket counts exact
    assert np.abs(Q.sum(0) - n * q).max() < 2.0  # columns: q within the two rows' rounding
    assert Q[0, 3] == 0 and Q[1, 2] == 0


LAB_B = ("epi", "meso", "exe")


def _geo_pair(seed=0, t=0.4, n_pool=(600, 700), labels_b=None, labels_a=None):
    rng = np.random.default_rng(seed)
    Cb, Ca = rng.normal(size=(n_pool[0], 3)), rng.normal(size=(n_pool[1], 3)) * 1.3
    Eb, Ea = rng.gamma(2.0, 1.0, size=(n_pool[0], 12)), rng.gamma(2.0, 1.0, size=(n_pool[1], 12))
    lb = rng.choice(["epi", "meso", "exe"], size=n_pool[0], p=[0.5, 0.3, 0.2]) if labels_b is None else labels_b
    la = rng.choice(["epi", "meso", "late"], size=n_pool[1], p=[0.2, 0.5, 0.3]) if labels_a is None else labels_a
    return Cb, Ca, Eb, Ea, np.asarray(lb), np.asarray(la), t


def test_geoshares_realised_shares_match_q_bracket_counts_stay_and_edits_are_minimal():
    Cb, Ca, Eb, Ea, lb, la, t = _geo_pair(seed=11)
    n = 500
    base = tb.interpolate_clouds(Cb, Ca, Eb, Ea, t, n, seed=5, **KW)
    geo = tb.interpolate_clouds(Cb, Ca, Eb, Ea, t, n, seed=5, geo_labels=(lb, la), geo_stream=GSTREAM, **KW)
    rec = geo.record["geoshares"]
    assert not rec["idle"] and rec["tv_q_linear"] > 0
    # the bracket counts equal the entry's, place by place
    np.testing.assert_array_equal(geo.stage, base.stage)
    assert geo.record["n_before"] == base.record["n_before"] and geo.record["n_after"] == base.record["n_after"]
    # realised label shares of the drawn cells match q within the rounding
    names = np.unique(np.r_[lb, la])
    pb = np.array([(lb == k).mean() for k in names])
    pa = np.array([(la == k).mean() for k in names])
    q = tb.geometric_shares(pb, pa, t)
    drawn = tb.drawn_labels(geo, lb, la)
    got = np.array([(drawn == k).sum() for k in names])
    assert np.abs(got - n * q).max() < 2.0
    assert rec["tv_q_linear"] == pytest.approx(0.5 * np.abs(q - ((1 - t) * pb + t * pa)).sum())
    # rows not edited are the entry's own rows at their own places; the edited share is what the record says
    changed = geo.rows != base.rows
    assert changed.mean() == pytest.approx(rec["edited_share"]) and 0 < changed.mean() < 0.5
    # each bracket's edits drop over-quota labels and add under-quota ones, so the labels of the unchanged rows agree
    base_lab = tb.drawn_labels(base, lb, la)
    np.testing.assert_array_equal(drawn[~changed], base_lab[~changed])
    # the guard's pool weights are the edited draw's chances and sum to 1
    assert geo.pool_weights.sum() == pytest.approx(1.0)
    assert geo.coords.shape == base.coords.shape and np.isfinite(geo.coords).all()


def test_geoshares_adds_unused_rows_before_repeating():
    Cb, Ca, Eb, Ea, lb, la, t = _geo_pair(seed=12)
    base = tb.interpolate_clouds(Cb, Ca, Eb, Ea, t, 500, seed=5, **KW)
    geo = tb.interpolate_clouds(Cb, Ca, Eb, Ea, t, 500, seed=5, geo_labels=(lb, la), geo_stream=GSTREAM, **KW)
    for s in (0, 1):
        m = geo.stage == s
        # the pools are larger than the draw: no row is repeated when an unused row of its label exists
        assert len(np.unique(geo.rows[m])) == int(m.sum()) == len(np.unique(base.rows[m]))


@pytest.mark.parametrize("case", ["identical", "t0", "t1"])
def test_geoshares_with_q_equal_to_the_linear_shares_is_fgw_bit_for_bit(case):
    Cb, Ca, Eb, Ea, lb, la, t = _geo_pair(seed=13)
    if case == "identical":  # the same label shares in both brackets (different cells, different pool sizes)
        lb = np.repeat(["epi", "meso", "exe"], [300, 180, 120])
        la = np.repeat(["epi", "meso", "exe"], [350, 210, 140])
    else:
        t = 0.0 if case == "t0" else 1.0
    base = tb.interpolate_clouds(Cb, Ca, Eb, Ea, t, 500, seed=5, **KW)
    geo = tb.interpolate_clouds(Cb, Ca, Eb, Ea, t, 500, seed=5, geo_labels=(lb, la), geo_stream=GSTREAM, **KW)
    for f in ("stage", "rows", "coords", "pool_coords", "pool_weights"):
        np.testing.assert_array_equal(getattr(geo, f), getattr(base, f))
    assert geo.record["geoshares"]["idle"] and geo.record["geoshares"]["edited_share"] == 0.0
    assert {k: v for k, v in geo.record.items() if k != "geoshares"} == base.record


def test_geoshares_is_deterministic_and_refuses_other_draws():
    Cb, Ca, Eb, Ea, lb, la, t = _geo_pair(seed=14)
    a = tb.interpolate_clouds(Cb, Ca, Eb, Ea, t, 500, seed=5, geo_labels=(lb, la), geo_stream=GSTREAM, **KW)
    b = tb.interpolate_clouds(Cb, Ca, Eb, Ea, t, 500, seed=5, geo_labels=(lb, la), geo_stream=GSTREAM, **KW)
    for f in ("stage", "rows", "coords", "pool_coords", "pool_weights"):
        np.testing.assert_array_equal(getattr(a, f), getattr(b, f))
    with pytest.raises(ValueError, match="uniform draw only"):
        tb.interpolate_clouds(Cb, Ca, Eb, Ea, t, 500, seed=5, geo_labels=(lb, la), geo_stream=GSTREAM, herd=True, **KW)
    with pytest.raises(ValueError, match="one label per cell"):
        tb.interpolate_clouds(Cb, Ca, Eb, Ea, t, 500, seed=5, geo_labels=(lb[:-1], la), geo_stream=GSTREAM, **KW)


def test_fgw_geoshares_candidate_keeps_fgws_count_and_prints_shares_only(capsys):
    ctx, task = _ctx(seed=3)
    fgw = CANDIDATES["fgw"].build(ctx)
    capsys.readouterr()
    v = CANDIDATES["fgw_geoshares"].build(ctx)
    line = capsys.readouterr().out.strip()
    assert line.startswith(f"fgw_geoshares check: task={task.key} seed=3 idle=False tv_q_linear=")
    fields = dict(kv.split("=") for kv in line.split(": ", 1)[1].split())
    assert set(fields) == {"task", "seed", "idle", "tv_q_linear", "edited_share"}
    assert 0 < float(fields["tv_q_linear"]) < 1 and 0 < float(fields["edited_share"]) < 1
    assert v.X.shape == fgw.X.shape and v.coords.shape == fgw.coords.shape
    # identical label shares in both brackets: the candidate is fgw, bit for bit
    ctx2, _ = _ctx(seed=4, shares=((0.5, 0.3, 0.2), (0.5, 0.3, 0.2)))
    for a in ctx2.inputs.values():  # the same exact shares: 450 / 270 / 180 of 900 in each bracket
        a.obs["group"] = np.repeat(["epi", "meso", "exe"], [450, 270, 180])
    np.testing.assert_array_equal(CANDIDATES["fgw_geoshares"].build(ctx2).X, CANDIDATES["fgw"].build(ctx2).X)
    np.testing.assert_array_equal(CANDIDATES["fgw_geoshares"].build(ctx2).coords, CANDIDATES["fgw"].build(ctx2).coords)
