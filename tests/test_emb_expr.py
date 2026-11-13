"""The embryo expression batch (campaign branch campaign/emb-expr): `fgw_trend_render`, `fgw_aligned_hi` and
`fgw_early_specialist`, on synthetic labelled stages built here. Nothing is read from data/, nothing is scored, and
no prediction is compared with any target.

Each arm differs from its base in one idea only: the same drawn cells, positions and cell number as its base (fgw, or
fgw_aligned for fgw_aligned_hi), and with the idea's strength at zero, or on an input where it has nothing to do, the
arm's file is its base's bit for bit. The existing candidates the arms reuse are pinned by sha256 of their built
expression and coordinates on the same synthetic context, computed with main at f63f964 before this batch's edit."""

from __future__ import annotations

import hashlib

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from vec import harness
from vec import t2_build as tb
from vec.boards import Board
from vec.harness import CANDIDATES, PseudoTask, load_tasks

PANEL = [f"g{i}" for i in range(24)]
OUTSIDE = [f"x{i}" for i in range(3)]
# label -> (lineage group, centre of its cells): one extraembryonic label both brackets hold, and an embryonic one
LABELS = {"EXE-Endoderm": ("Endoderm", (0.0, 0.0, 0.0)), "SOM": ("ParaxialMesoderm", (6.0, 0.0, 0.0)),
          "PS": ("Mesoderm", (0.0, 6.0, 0.0))}
# stage -> (cells per label, median depth, run factor on panel genes 0-5, development switch on genes 12-17 in SOM)
STAGES = {"earlier": (140, 30.0, 0.8, 1.0), "before": (160, 40.0, 1.0, 2.0), "after": (220, 160.0, 2.2, 5.0)}
DAYS = {"earlier": 6.75, "before": 7.25, "after": 8.0}
FILES = {"earlier": "T2/E6.75.h5ad", "before": "T2/E7.25.h5ad", "after": "T2/E8.0.h5ad"}

# sha256 of (X float32, coords float32) built by each existing candidate on `ctx()`, computed with main at f63f964
# before this batch's edit (OPENBLAS/OMP/MKL threads 4)
PRE_EDIT_SHA256 = {
    "fgw": "6aff89bee3a37e2e0ffdebec29584091733a115596a1e8445d90599d2ed470d6",
    "fgw_aligned": "555b101e9078b37646ab79534b668bddb811575d3da4fe4ed3849e14aec4bea6",
    "fgw_render_hi": "a7266e10d066ee44b1b035c8eac332c87427042ee598c3e4f3334c09d206016a",
    "fgw_series": "a0bcb83d4b73d44e31734c535e7fe07d4cba04e1a569b011bfb2a1e87d7975e2",
    "fgw_aligned_series": "747804b40dfc25b95aa12437c3a9cd7a926fe18ac6372c13bf79839ea0ac1f9b",
}


def _stage(name: str, seed: int, run: float | None = None) -> ad.AnnData:
    """A released MERFISH stage as vec.prep.load_stage returns it: panel and outside genes, counts layer, X = log(1 +
    1e4 c / N) with N over every gene, obs celltype and group, 3D coordinates. `run` overrides the stage's run
    factor."""
    per, depth, run0, dev = STAGES[name]
    run = run0 if run is None else run
    rng = np.random.default_rng(seed)
    base = np.random.default_rng(5).gamma(2.0, 1.0, len(PANEL) + len(OUTSIDE))
    parts, obs, xyz = [], [], []
    for lab, (grp, centre) in LABELS.items():
        prof = base.copy()
        prof[:6] *= run
        if lab == "SOM":
            prof[12:18] *= dev
        prof /= prof.sum()
        d = np.maximum(np.round(np.exp(rng.normal(np.log(depth), 0.3, per))), 5)[:, None]
        C = rng.poisson(prof * d)
        C[C[:, :len(PANEL)].sum(1) == 0, 0] += 1
        parts.append(C)
        obs += [(lab, grp)] * per
        xyz.append(np.asarray(centre) + rng.normal(size=(per, 3)) * 1.5)
    C = np.vstack(parts)
    N = C.sum(1, keepdims=True)
    a = ad.AnnData(X=sparse.csr_matrix(np.log1p(1e4 * C / N).astype(np.float32)),
                   obs=pd.DataFrame(obs, columns=["celltype", "group"], index=[f"{name}{i}" for i in range(len(C))]),
                   var=pd.DataFrame(index=PANEL + OUTSIDE))
    a.layers["counts"] = sparse.csr_matrix(C.astype(np.int32))
    a.obsm["spatial_3D"] = np.vstack(xyz).astype(np.float32)
    return a


def stages(**runs) -> dict[str, ad.AnnData]:
    return {name: _stage(name, seed, runs.get(name)) for seed, name in enumerate(STAGES, start=11)}


BOARD = Board(key="pseudo:synthetic_embryo", task="T2", genes=list(PANEL), needs_coords=True, setting="embryo",
              min_cells=100, max_cells=400)


def task(roles=("earlier", "before", "after"), key="syn_emb_expr") -> PseudoTask:
    return PseudoTask(key=key, test_board="T2:embryo:test_interp", task="T2", setting="embryo", kind="interp",
                      target="T2/E7.5.h5ad", reference="T2/E7.25.h5ad", inputs={r: FILES[r] for r in roles},
                      twin_board="T2:embryo:val_interp")


def ctx(inputs=None, roles=("earlier", "before", "after"), seed=harness.BUILD_SEED_OFFSET + 4):
    inputs = stages() if inputs is None else inputs
    return harness.BuildContext(task=harness.TaskView.of(task(roles)), board=BOARD,
                                inputs={r: inputs[r] for r in roles}, seed=seed)


def digest(st) -> str:
    h = hashlib.sha256(np.ascontiguousarray(st.X, dtype=np.float32).tobytes())
    h.update(np.ascontiguousarray(st.coords, dtype=np.float32).tobytes())
    return h.hexdigest()


def _interp(c, move_coupling="fgw"):
    """fgw's own interpolation on context c (the draws, positions and brackets every arm keeps)."""
    t, g = c.task, list(BOARD.genes)
    before, after = c.inputs["before"], c.inputs["after"]
    Eb, Ea = (harness.dense(a[:, g].X).astype(np.float32) for a in (before, after))
    frac = (t.target_day - t.role_days["before"]) / (t.role_days["after"] - t.role_days["before"])
    return harness.t2_build.interpolate_clouds(
        np.asarray(before.obsm["spatial_3D"], float)[:, :3], np.asarray(after.obsm["spatial_3D"], float)[:, :3],
        Eb, Ea, frac, harness.count_for(BOARD, before.n_obs), c.seed, target_rms=harness.target_rms(c),
        move_coupling=move_coupling), Eb, Ea, frac


def _depths(monkeypatch, before_lib: float, after_lib: float, c) -> None:
    """Every cell of the before bracket at library size before_lib and of the after bracket at after_lib, as
    `_stage_counts` reports them (expression and counts as released): the depth step's strength set by hand."""
    real = tb._stage_counts

    def fixed(a, genes):
        E, C, lib, other = real(a, genes)
        return E, C, np.full_like(lib, before_lib if a is c.inputs["before"] else after_lib), other

    monkeypatch.setattr(tb, "_stage_counts", fixed)


@pytest.fixture(scope="module")
def built():
    c = ctx()
    names = (*PRE_EDIT_SHA256, "fgw_trend_render", "fgw_aligned_hi", "fgw_early_specialist")
    return c, {n: CANDIDATES[n].build(c) for n in names}


# ------------------------------------------------------------------------------------- existing candidates unchanged
@pytest.mark.parametrize("name", sorted(PRE_EDIT_SHA256))
def test_existing_candidates_are_bit_for_bit_what_they_were(built, name):
    assert digest(built[1][name]) == PRE_EDIT_SHA256[name]


# ------------------------------------------------------------------------------------------------ registration
def test_registration_and_where_each_arm_applies():
    tasks = {t.key: t for ts in load_tasks().values() for t in ts}
    counted = (tasks["T2e/E7.25+E8.0->E7.5"], tasks["T2e/E6.75+E7.5->E7.25"])
    heart = tasks["T2h/<=E9.5->E10.5"]
    two = task(roles=("before", "after"))
    for name in ("fgw_trend_render", "fgw_aligned_hi", "fgw_early_specialist"):
        c = CANDIDATES[name]
        assert c.name == name and c.tasks == {"T2"} and not c.tuned_on and not c.externals and not c.derived
        assert all(c.applies(t) for t in counted) and not c.applies(heart)
        assert c.params["scope"] == "whole" and c.params["rule"]
    # the trend needs a third sample: on a two-bracket task it does not apply, the other two arms do
    assert not CANDIDATES["fgw_trend_render"].applies(two)
    assert CANDIDATES["fgw_aligned_hi"].applies(two) and CANDIDATES["fgw_early_specialist"].applies(two)
    assert CANDIDATES["fgw_trend_render"].params["rule"] == tb.TREND_RENDER_RULE
    assert CANDIDATES["fgw_early_specialist"].params["rule"] == tb.EARLY_RULE
    assert CANDIDATES["fgw_aligned_hi"].params["rule"] == tb.RENDER_RULE


# ------------------------------------------------------------------- same cells, positions and count as the base
@pytest.mark.parametrize("arm, base", [("fgw_trend_render", "fgw"), ("fgw_early_specialist", "fgw"),
                                       ("fgw_aligned_hi", "fgw_aligned")])
def test_each_arm_keeps_its_bases_cells_positions_and_cell_number(built, arm, base):
    _, b = built
    assert b[arm].n == b[base].n == BOARD.max_cells
    np.testing.assert_array_equal(b[arm].coords, b[base].coords)
    assert not np.array_equal(b[arm].X, b[base].X)  # the expression is the arm's one idea


def test_aligned_hi_is_render_hi_at_the_aligned_positions(built):
    """fgw_aligned_hi: fgw_aligned's positions, fgw_render_hi's expression, cell for cell (the same draws: the moves'
    coupling changes positions only)."""
    _, b = built
    np.testing.assert_array_equal(b["fgw_aligned_hi"].X, b["fgw_render_hi"].X)
    np.testing.assert_array_equal(b["fgw_aligned_hi"].coords, b["fgw_aligned"].coords)
    assert not np.array_equal(b["fgw_aligned_hi"].coords, b["fgw_render_hi"].coords)


def test_aligned_hi_with_equal_depths_is_fgw_aligned_bit_for_bit(monkeypatch):
    """With the top-up at strength zero (both brackets at one depth, so RENDER_RULE's ratio is 1 for both) the arm is
    its base: fgw_aligned."""
    c = ctx()
    _depths(monkeypatch, 100.0, 100.0, c)
    assert digest(CANDIDATES["fgw_aligned_hi"].build(c)) == PRE_EDIT_SHA256["fgw_aligned"]


# ----------------------------------------------------------------------------------------------- the trend render
def test_trend_factor_reads_the_line_through_every_sample_at_the_target_day():
    rng = np.random.default_rng(3)
    C = [sparse.csr_matrix(rng.poisson(5.0, size=(200, 6)) + 1) for _ in range(3)]
    names = [np.array(["EXE-A"] * 200)] * 3
    labels = ["EXE-A"]
    days = [6.75, 7.25, 8.0]
    phi = tb.trend_factor(C[1], names[1], list(zip(C, names, strict=True)), days, 7.5, labels)

    def prof(M):
        M = M.toarray().astype(float)
        return (M / np.maximum(M.sum(1, keepdims=True), 1)).mean(0)

    Y = np.log(np.vstack([prof(M) for M in C]) + tb.SERIES_EPS)
    fit = np.polyfit(days, Y, 1)
    den = prof(C[1])
    want = np.exp(fit[0] * 7.5 + fit[1]) / (den + tb.SERIES_EPS)
    want = np.clip(want / (den @ want / den.sum()), *tb.SERIES_CLIP)
    np.testing.assert_allclose(phi, want, rtol=1e-9)
    # one profile in every sample: the factor is 1 to rounding
    same = tb.trend_factor(C[0], names[0], [(C[0], names[0])] * 3, days, 7.5, labels)
    np.testing.assert_allclose(same, 1.0, rtol=1e-12)
    with pytest.raises(ValueError, match="no extraembryonic label"):
        tb.trend_factor(C[0], names[0], list(zip(C, names, strict=True)), days, 7.5, [])
    with pytest.raises(ValueError, match="distinct days"):
        tb.trend_factor(C[0], names[0], list(zip(C, names, strict=True)), [7.0, 7.0, 7.0], 7.5, labels)
    with pytest.raises(ValueError, match="two samples"):
        tb.trend_labels([names[0]], np.array(["Endoderm"] * 200))


def test_trend_render_moves_both_brackets_towards_the_trend_and_leaves_development_alone(built):
    """Run factors 0.8 (E6.75), 1.0 (E7.25) and 2.2 (E8.0) on genes 0-5: the trend at E7.5 lies above the before
    bracket's and below the after bracket's, so genes 0-5's share rises in the before cells and falls in the after
    cells. The embryonic label's development switch (genes 12-17 of SOM) is not in the factor."""
    c, b = built
    res, Eb, Ea, _ = _interp(c)
    X = b["fgw_trend_render"].X
    _, rec = tb.trend_rendered_expression(c.inputs["before"], c.inputs["after"], res, list(BOARD.genes),
                                          [(a, c.task.role_days[r]) for r, a in c.inputs.items()],
                                          c.task.target_day, [c.seed, 7])
    assert rec["labels"] == ["EXE-Endoderm"] and rec["moved"]["before"] > 0 and rec["moved"]["after"] > 0
    for s, E, sign in ((0, Eb, 1), (1, Ea, -1)):
        m = res.stage == s
        new, old = np.expm1(X[m].astype(float)), np.expm1(E[res.rows[m]].astype(float))
        share = (new[:, :6].sum(1) / new.sum(1)).mean() / (old[:, :6].sum(1) / old.sum(1)).mean()
        assert sign * (share - 1.0) > 0.05
        som = np.asarray(c.inputs["before" if s == 0 else "after"].obs["celltype"])[res.rows[m]] == "SOM"
        dev = ((new[som][:, 12:18].sum(1) / new[som].sum(1)).mean()
               / (old[som][:, 12:18].sum(1) / old[som].sum(1)).mean())
        assert abs(np.log(dev)) < 0.3
    assert np.isfinite(X).all() and (X >= 0).all() and ((X > 0).sum(1) >= 1).all()


def test_trend_render_at_strength_zero_is_fgw_bit_for_bit(monkeypatch):
    """A factor of exactly 1 on every gene (no trend to follow) leaves every row as released: the arm is fgw."""
    monkeypatch.setattr(tb, "trend_factor", lambda C, *a, **k: np.ones(C.shape[1]))
    c = ctx()
    assert digest(CANDIDATES["fgw_trend_render"].build(c)) == PRE_EDIT_SHA256["fgw"]


def test_trend_render_refuses_brackets_outside_its_samples(built):
    c, _ = built
    res, *_ = _interp(c)
    with pytest.raises(ValueError, match="among the trend's samples"):
        tb.trend_rendered_expression(c.inputs["before"], c.inputs["after"], res, list(BOARD.genes),
                                     [(c.inputs["earlier"], 6.75), (c.inputs["before"], 7.25)], 7.5, 0)


# ----------------------------------------------------------------------------------------- the earlier-side specialist
def test_early_specialist_thins_the_deeper_bracket_to_the_shallower_profile_and_depth(built):
    c, b = built
    res, Eb, Ea, frac = _interp(c)
    X = b["fgw_early_specialist"].X
    lo, hi = res.stage == 0, res.stage == 1  # before (40) is the shallower bracket, after (160) the deeper
    np.testing.assert_array_equal(X[lo], Eb[res.rows[lo]])  # the shallower bracket keeps its released rows
    _, rec = tb.early_specialist_expression(c.inputs["before"], c.inputs["after"], res, list(BOARD.genes), frac,
                                            [c.seed, 7])
    assert rec["moved"] == "after" and rec["labels"] == ["EXE-Endoderm"] and rec["ratio"] < 0.5
    assert rec["rendered_median"] == pytest.approx(rec["median_depth"]["before"], rel=0.2)
    new, old = np.expm1(X[hi].astype(float)), np.expm1(Ea[res.rows[hi]].astype(float))
    assert (X[hi] > 0).sum(1).mean() < (Ea[res.rows[hi]] > 0).sum(1).mean()  # fewer detected genes: thinned
    share = (new[:, :6].sum(1) / new.sum(1)).mean() / (old[:, :6].sum(1) / old.sum(1)).mean()
    assert share < 0.8  # the deeper bracket's run factor 2.2 is taken back towards the shallower one's 1.0
    assert np.isfinite(X).all() and (X >= 0).all() and ((X > 0).sum(1) >= 1).all()


def test_early_specialist_tops_up_only_the_genes_whose_keep_factor_exceeds_one(monkeypatch):
    c = ctx()
    res, Eb, Ea, frac = _interp(c)
    _depths(monkeypatch, 100.0, 200.0, c)
    phi = np.ones(len(BOARD.genes))
    phi[0] = 4.0  # q = 4 x 0.5 = 2 on gene 0; 0.5 on every other gene
    monkeypatch.setattr(tb, "series_factor", lambda *a, **k: phi)
    X, rec = tb.early_specialist_expression(c.inputs["before"], c.inputs["after"], res, list(BOARD.genes), frac, 1)
    assert rec["ratio"] == 0.5 and rec["genes_topped_up"] == 1
    hi = res.stage == 1
    C0 = c.inputs["after"][:, list(BOARD.genes)].layers["counts"].toarray()[res.rows[hi]]
    assert np.isfinite(X).all() and ((X > 0).sum(1) >= 1).all()
    assert (X[hi][:, 1:] > 0).sum() < (C0[:, 1:] > 0).sum()  # the other genes are thinned: fewer detections
    assert (X[hi][:, 0] > 0).sum() >= (C0[:, 0] > 0).sum() * 0.95  # gene 0 is topped up, not thinned


def test_early_specialist_at_strength_zero_is_fgw_bit_for_bit(monkeypatch):
    """Both brackets at one depth (ratio 1) and a factor of exactly 1: every keep factor is 1, so the arm is fgw."""
    c = ctx()
    _depths(monkeypatch, 100.0, 100.0, c)
    monkeypatch.setattr(tb, "series_factor", lambda C_from, *a, **k: np.ones(C_from.shape[1]))
    assert digest(CANDIDATES["fgw_early_specialist"].build(c)) == PRE_EDIT_SHA256["fgw"]


def test_both_factor_arms_refuse_samples_without_a_matched_label():
    s = stages()
    for a in s.values():
        a.obs["celltype"] = a.obs["celltype"].astype(str).str.replace("EXE-", "EMB-")
    c = ctx(inputs=s)
    with pytest.raises(ValueError, match="no matched extraembryonic label"):
        CANDIDATES["fgw_early_specialist"].build(c)
    with pytest.raises(ValueError, match="no extraembryonic label held by every sample"):
        CANDIDATES["fgw_trend_render"].build(c)
