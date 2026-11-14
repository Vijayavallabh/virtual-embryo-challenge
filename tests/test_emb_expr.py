"""The embryo expression batch (campaign branch campaign/emb-expr): `fgw_trend_render`, `fgw_aligned_hi` and
`fgw_early_specialist`, on synthetic labelled stages built here. Nothing is read from data/, nothing is scored, and
no prediction is compared with any target.

Each arm differs from its base in one idea only: the same drawn cells, positions and cell number as its base (fgw, or
fgw_aligned for fgw_aligned_hi), and with the idea's strength at zero, or on an input where it has nothing to do, the
arm's file is its base's bit for bit. The existing candidates the arms reuse are pinned by sha256 of their built
expression and coordinates on the same synthetic context, computed with main at f63f964 before this batch's edit.

Tries E35a and E35b (campaign/emb-trend2) add `fgw_aligned_trend` and trend-sample roles; every pre-existing candidate
that builds on the synthetic contexts is pinned at e709655, before that edit (PINS_AT_E709655)."""

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


# sha256 of the batch's three arms on `ctx()` at ede98cd, before the empty-cell rule was added to the early specialist
# (OPENBLAS/OMP/MKL threads 4): where no cell is left empty the arms are bit for bit what they were
ARMS_AT_EDE98CD = {
    "fgw_early_specialist": "50c27daf50514a078346b8175914731fc46343ba0a99aba9d85936c939281565",
    "fgw_aligned_hi": "4447d06a545ab542e6bc8bd22bff814eee05fbd03de3404fd68aa64ba5ec6f7b",
    "fgw_trend_render": "a380d451db01db4d3f8ab02a6bea78eb5ed803708d620d84cecd182a75dfb7b7",
}


@pytest.mark.parametrize("name", sorted(ARMS_AT_EDE98CD))
def test_the_empty_cell_rule_changes_no_arm_where_no_cell_is_left_empty(name):
    c = ctx()
    assert digest(CANDIDATES[name].build(c)) == ARMS_AT_EDE98CD[name]
    if name == "fgw_early_specialist":
        res, _, _, frac = _interp(c)
        _, rec = tb.early_specialist_expression(c.inputs["before"], c.inputs["after"], res, list(BOARD.genes), frac,
                                                [c.seed, 7])
        assert rec["kept_as_released_share"] == 0.0


def test_early_specialist_keeps_the_released_row_of_a_cell_no_retry_can_fill(monkeypatch):
    """A depth ratio so small that no drawn cell of the deeper bracket keeps a panel count in any retry: the build
    does not stop, every such cell keeps its released row, and the record says what share that was."""
    c = ctx()
    res, Eb, Ea, frac = _interp(c)
    _depths(monkeypatch, 1e-9, 1.0, c)  # ratio 1e-9: with the factor at 1, every keep probability is 1e-9
    monkeypatch.setattr(tb, "series_factor", lambda C_from, *a, **k: np.ones(C_from.shape[1]))
    X, rec = tb.early_specialist_expression(c.inputs["before"], c.inputs["after"], res, list(BOARD.genes), frac, 3)
    lo, hi = res.stage == 0, res.stage == 1
    assert rec["moved"] == "after" and rec["kept_as_released_share"] == 1.0
    np.testing.assert_array_equal(X[hi], Ea[res.rows[hi]])  # the deeper bracket's cells keep their released rows
    np.testing.assert_array_equal(X[lo], Eb[res.rows[lo]])
    assert np.isfinite(X).all() and ((X > 0).sum(1) >= 1).all()


def test_early_specialist_mixes_moved_and_kept_rows_cell_by_cell(monkeypatch):
    """Some cells keep a count and some never do: the kept ones are thinned as before, the empty ones are released
    rows, and the retry helper still refuses an empty row for every caller that does not ask for it."""
    c = ctx()
    res, Eb, Ea, frac = _interp(c)
    _depths(monkeypatch, 1.0, 40.0, c)  # ratio 0.025: a cell with few panel counts is often left empty
    monkeypatch.setattr(tb, "series_factor", lambda C_from, *a, **k: np.ones(C_from.shape[1]))
    monkeypatch.setattr(tb, "DEPTH_ZERO_RETRIES", 0)
    X, rec = tb.early_specialist_expression(c.inputs["before"], c.inputs["after"], res, list(BOARD.genes), frac, 3)
    hi = res.stage == 1
    same = (X[hi] == Ea[res.rows[hi]]).all(1)
    assert 0.0 < rec["kept_as_released_share"] < 1.0 and same.mean() == pytest.approx(rec["kept_as_released_share"])
    assert ((X[hi][~same] > 0).sum(1) <= (Ea[res.rows[hi]][~same] > 0).sum(1)).all()  # the others are thinned
    assert np.isfinite(X).all() and ((X > 0).sum(1) >= 1).all()

    def never(idx, g):
        return sparse.csr_matrix((len(idx), 3), dtype=np.int64), np.zeros(len(idx))

    with pytest.raises(ValueError, match="without panel counts after every retry"):
        tb._keep_some_counts(never, 4, np.random.default_rng(0))
    C, lib, _ = tb._keep_some_counts(never, 4, np.random.default_rng(0), allow_empty=True)
    assert C.shape == (4, 3) and C.nnz == 0


# ------------------------------------------------- tries E35a and E35b: trend samples and the aligned trend arm
# sha256 of (X float32, coords float32) of every pre-existing candidate that builds on the synthetic contexts, with
# campaign/emb-exprfix at e709655 before the edit of tries E35a and E35b (OPENBLAS/OMP/MKL threads 4). fgw_latent
# needs scvi-tools, which this environment lacks; b0 and b0_herd need the lineage map's `merged` column, which the
# synthetic stages lack (both refusals are pinned below).
PINS_AT_E709655 = {
    ("earlier", "before", "after"): {
        "copy": "d5d07c592f4190fa9726e6891c9fa84e4fb790ff97755b1e69fe85ab3792814e",
        "pseudobulk_shift": "6fec7ada2a2e0fc4ebc882c97378d41c8fc33cc633b47e57e12da4f69ae48414",
        "psb_size": "ced60d6ee26f4c671801693f47d5ef6bf3095eb825cca220d0d4903daf38d880",
        "fgw": "6aff89bee3a37e2e0ffdebec29584091733a115596a1e8445d90599d2ed470d6",
        "fgw_shape0": "9d9d06bc5eada31e28572da397205855be9945db2d6caa37fb6b5d650dbc669d",
        "fgw_shape1": "778957272eba5d1dbe5519b18fecad09908898fa4ea4ed7de267f03c4c932e3f",
        "fgw_depth": "2a3bc3b13c588893de098ceb077db38e2a06f638bb200b04ffd94986afa8f530",
        "fgw_series": "a0bcb83d4b73d44e31734c535e7fe07d4cba04e1a569b011bfb2a1e87d7975e2",
        "fgw_shape0_hi": "7d62c8d1074f67b3dc3d5d8f293089fed176bd0d94ea4856f3851a703a11a198",
        "fgw_shape0_series": "faf1aa70df05ca2a0e7ee88b1cf6e86cea66a8c84f62176043ec532fba2c449a",
        "fgw_aligned_series": "747804b40dfc25b95aa12437c3a9cd7a926fe18ac6372c13bf79839ea0ac1f9b",
        "fgw_aligned_series_wide": "281d16ebc9fd1bce50ccf8b8c63e85241687315d5b2d80ef3411ca3d07513614",
        "fgw_trend_render": "a380d451db01db4d3f8ab02a6bea78eb5ed803708d620d84cecd182a75dfb7b7",
        "fgw_aligned_hi": "4447d06a545ab542e6bc8bd22bff814eee05fbd03de3404fd68aa64ba5ec6f7b",
        "fgw_early_specialist": "50c27daf50514a078346b8175914731fc46343ba0a99aba9d85936c939281565",
        "fgw_shape1_series": "ecfe79d8dc18f0330c9782f4313b768ac3af3ccc05dbbf2f022790869d80d4ec",
        "fgw_late_shape1": "f2bc7e197b5c70e3891e97ce99109e8a385026a84751aad774c987bf0280393a",
        "fgw_late": "520c29a0c66f3db2f6c58be5e84fbb148ee6ab08384d6fab79878bbed91bde64",
        "fgw_render_mid": "16cda62b90a93302baba00f5c3a9ad39eb573aacb1c8a1ecbca466d8910f46e8",
        "fgw_render_lo": "59ad777cf35f866ca980af05efc2956371e2c10866cc32de6fea6a99a414479b",
        "fgw_render_hi": "a7266e10d066ee44b1b035c8eac332c87427042ee598c3e4f3334c09d206016a",
        "fgw_blend": "c138d23dc9556f0a23c61bb208718394964da3d856d780bd86239e5deafba9e0",
        "fgw_herd": "28727e19047a59c3f1cd4980cea8eabe9908a812c6df08b3fda6a595095e655d",
        "fgw_aligned": "555b101e9078b37646ab79534b668bddb811575d3da4fe4ed3849e14aec4bea6",
        "fgw_shape0_aligned": "00a0802a6d0f9468c12133013ca75f7397e4d8eac206729544af72e034f94621",
        "fgw_perm": "4c7867c5e5405b7a8d2ad97b2ffc66625ca623004f3eaef74c5a67161b1b737e",
        "fgw_midpoint": "929b3c8ee11e7129d528411fa863bc2ec17fb7525e36781149b65eccf83f980f",
    },
    ("before", "after"): {
        "copy": "d5d07c592f4190fa9726e6891c9fa84e4fb790ff97755b1e69fe85ab3792814e",
        "fgw": "4a53a67d0be4e5195b3ddf8e554b5ac2ac0b506bfb8808f747ab48cea1d8164f",
        "fgw_shape0": "56e1a92b9cc1d160f60e39ad22d54aa1df8700111c4cd57ce85536a21cc455db",
        "fgw_shape1": "30423881b2d0cbf81f22cba905e65df8071fd3754d1971856321ddae203ae210",
        "fgw_depth": "923cb6a98d679d0f2f1ca874f538505352df7b2fd77acad86a315307e0825000",
        "fgw_series": "2b9230fa43fc5979dfeab8ff542835d385ddf768f8dc6a4e15e2f5399286bf3c",
        "fgw_shape0_hi": "ea169e3a237c872b637ba9268033b85e787b1dbc0cfbc42fdb2c4e999ec03d78",
        "fgw_shape0_series": "c124ad333a5c6bf49ce5a6e8a201c7ebd37247c197c99f243af43a692476406d",
        "fgw_aligned_series": "d552e97cec14166c5f0e67199a8531509a2f6ae761d755732c2ac7409fa9b879",
        "fgw_aligned_series_wide": "df5bbf531b51660afea22e195e8d2e62fc63de422f19fe3f4ede7febf9a225e5",
        "fgw_aligned_hi": "ae2604fece0c6d63e15263e86a050886eb760bb060eb1d0307dffddd5402d4f1",
        "fgw_early_specialist": "9d4e8fd1247df9e028658e04a6db6dc37e1317d637cf8e5ca0ff814691940441",
        "fgw_shape1_series": "94563ff5cdf5fb42a0772f3401c2f79cbde0c44ca4165b6c2f5b90eba59b4b57",
        "fgw_late_shape1": "1c54058aaccc90061ab83afd0223f8241716b715f0607f1a53f8332e02ec2b90",
        "fgw_late": "a452e42bbdd4cd932d4b36df6ecc635603dbf16b1799836cc2ab696bfba4fd5e",
        "fgw_render_mid": "381dd143998fce19a858b42484c0cb3e049885ca1242e70463a599f60548477f",
        "fgw_render_lo": "1e7afa57c3187258190b30e18d4a00bba885c5d02249399b99274b5fc302e604",
        "fgw_render_hi": "b4b33e85d5804b6dcea60aa1290dbf03610a4b8e576c7937318d01f39baee362",
        "fgw_blend": "ca4211cdbf0aa5af168cabc94d19ef6f2afd4bed8d533e084289c524b600da33",
        "fgw_herd": "ecd9a690452084639e7a9a7b305781a4bf1ce429d27094ca92930fd78b34cd2c",
        "fgw_aligned": "b533295c10151cea7c6ddbcd4935bdf2de08584bf8dbf5854f7796a4ed803c10",
        "fgw_shape0_aligned": "d4a560a8d1399f5dc9ae76fb8ee39cdb40355c1075e10a86b543a08437a71ed9",
        "fgw_perm": "a580e416e2bc47bcfb50df2ae5dcfb083f6ccc6a57714ac1447efea03c9dcb0b",
        "fgw_midpoint": "58651884a5805d8d152dec5905f51c039397db409b2af086ab0809d0ada022d0",
    },
}
NEW_AT_E35 = ("fgw_aligned_trend",)  # registered by the edit of tries E35a and E35b
NEW_AT_POLAR = ("fgw_polar", "fgw_aligned_polar_series",  # campaign/emb-polar (tests/test_emb_polar.py)
                "fgw_polar_frame", "fgw_aligned_polar_frame_series")
NEW_AT_MIRROR = ("fgw_mirrorfix", "fgw_mirrorfix_series",  # campaign/emb-mirror (tests/test_emb_mirror.py)
                 "fgw_mirrorgate", "fgw_mirrorgate_series")
NEW_AT_EMB_D = ("fgw_tissue", "fgw_tissue_series",  # campaign/emb-d (tests/test_emb_batch_d.py)
                "fgw_aligned_series_xt")
NEW_AT_EMB_F = ("fgw_potts", "fgw_partial")  # campaign/emb-f (tests/test_emb_batch_f.py)
NEW_AT_EMB_G = ("fgw_damp", "fgw_part_rand", "fgw_still",  # campaign/emb-g (tests/test_emb_batch_g.py)
                "fgw_damp_m")


@pytest.fixture(scope="module")
def pinned_builds():
    s = stages()
    return {roles: {n: digest(CANDIDATES[n].build(ctx(inputs=s, roles=roles))) for n in pins}
            for roles, pins in PINS_AT_E709655.items()}


@pytest.mark.parametrize("roles, name", [(r, n) for r, pins in PINS_AT_E709655.items() for n in sorted(pins)])
def test_every_pre_existing_candidate_is_what_it_was_at_e709655(pinned_builds, roles, name):
    assert pinned_builds[roles][name] == PINS_AT_E709655[roles][name]


def test_the_pins_cover_every_pre_existing_candidate_that_builds_here():
    """Every candidate applicable to a synthetic context is pinned, refused as pinned, or new at E35 or with the
    polar arms."""
    for roles, pins in PINS_AT_E709655.items():
        t = task(roles)
        here = {n for n, c in CANDIDATES.items() if n != harness.FLOOR and c.applies(t)}
        assert here - set(pins) <= {"b0", "b0_herd", "fgw_latent", *NEW_AT_E35, *NEW_AT_POLAR, *NEW_AT_MIRROR,
                                    *NEW_AT_EMB_D, *NEW_AT_EMB_F, *NEW_AT_EMB_G}
    c = ctx()
    for name in ("b0", "b0_herd"):
        with pytest.raises(KeyError, match="merged"):
            CANDIDATES[name].build(c)


def trend_task(extra=("earlier",)) -> PseudoTask:
    """The two-bracket synthetic task with each stage in `extra` handed over as a trend sample (role trend_<stem>)."""
    inputs = {"before": FILES["before"], "after": FILES["after"]}
    inputs.update({harness.TREND_SAMPLE_PREFIX + harness.stem(FILES[r]): FILES[r] for r in extra})
    return PseudoTask(key="syn_emb_trend", test_board="T2:embryo:test_interp", task="T2", setting="embryo",
                      kind="interp", target="T2/E7.5.h5ad", reference="T2/E7.25.h5ad", inputs=inputs,
                      twin_board="T2:embryo:val_interp")


def trend_ctx(inputs=None, extra=("earlier",), seed=harness.BUILD_SEED_OFFSET + 4):
    inputs = stages() if inputs is None else inputs
    t = trend_task(extra)
    by_file = {FILES[r]: inputs[r] for r in FILES}
    return harness.BuildContext(task=harness.TaskView.of(harness.check_task(t)), board=BOARD,
                                inputs={r: by_file[f] for r, f in t.inputs.items()}, seed=seed)


@pytest.mark.parametrize("name", ["fgw", "fgw_aligned"])
def test_trend_samples_change_nothing_a_bracket_only_path_builds(name):
    """fgw and fgw_aligned on the task with a trend sample are their builds on the two-bracket task, bit for bit: the
    size rule does not read a trend-sample role."""
    s = stages()
    two = CANDIDATES[name].build(ctx(inputs=s, roles=("before", "after")))
    ext = trend_ctx(inputs=s)
    assert harness.target_rms(ext) == harness.target_rms(ctx(inputs=s, roles=("before", "after")))
    assert digest(CANDIDATES[name].build(ext)) == digest(two) == PINS_AT_E709655[("before", "after")][name]


def test_target_rms_still_reads_every_configured_role():
    """The configured third stage ("earlier") stays in the size rule, as on the two counted backtests."""
    s = stages()
    assert harness.target_rms(ctx(inputs=s)) != harness.target_rms(ctx(inputs=s, roles=("before", "after")))


def test_trend_render_on_trend_samples_keeps_fgws_two_bracket_positions():
    s = stages()
    two = CANDIDATES["fgw"].build(ctx(inputs=s, roles=("before", "after")))
    ext = trend_ctx(inputs=s)
    assert CANDIDATES["fgw_trend_render"].applies(trend_task()) and not CANDIDATES["fgw_trend_render"].applies(
        task(("before", "after")))
    tr = CANDIDATES["fgw_trend_render"].build(ext)
    np.testing.assert_array_equal(tr.coords, two.coords)
    assert tr.n == two.n and not np.array_equal(tr.X, two.X)


def test_aligned_trend_is_registered_beside_the_trend_render():
    names = list(CANDIDATES)
    assert names.index("fgw_aligned_trend") == names.index("fgw_trend_render") + 1
    c, tr = CANDIDATES["fgw_aligned_trend"], CANDIDATES["fgw_trend_render"]
    assert c.name == "fgw_aligned_trend" and c.tasks == {"T2"} and not c.tuned_on and not c.externals
    assert not c.derived and c.requires is tr.requires
    assert c.params["rule"] == tb.TREND_RENDER_RULE and c.params["alpha"] == tb.INTERP_ALPHA
    assert c.params["scope"] == "whole"
    tasks = {t.key: t for ts in load_tasks().values() for t in ts}
    assert c.applies(tasks["T2e/E7.25+E8.0->E7.5"]) and c.applies(tasks["T2e/E6.75+E7.5->E7.25"])
    assert not c.applies(tasks["T2h/<=E9.5->E10.5"]) and not c.applies(task(("before", "after")))
    assert c.applies(trend_task())


@pytest.mark.parametrize("make", [ctx, trend_ctx])
def test_aligned_trend_is_the_trend_render_at_the_aligned_positions(make):
    """fgw_aligned_trend: fgw_aligned's positions, fgw_trend_render's expression, cell for cell (the same draws: the
    moves' coupling changes positions only)."""
    c = make()
    b = {n: CANDIDATES[n].build(c) for n in ("fgw_aligned_trend", "fgw_trend_render", "fgw_aligned")}
    r0, r1 = _interp(c)[0], _interp(c, move_coupling="aligned")[0]
    np.testing.assert_array_equal(r0.rows, r1.rows)
    np.testing.assert_array_equal(r0.stage, r1.stage)
    np.testing.assert_array_equal(b["fgw_aligned_trend"].X, b["fgw_trend_render"].X)
    np.testing.assert_array_equal(b["fgw_aligned_trend"].coords, b["fgw_aligned"].coords)
    assert not np.array_equal(b["fgw_aligned_trend"].coords, b["fgw_trend_render"].coords)
    assert b["fgw_aligned_trend"].n == BOARD.max_cells


def test_aligned_trend_at_strength_zero_is_fgw_aligned_bit_for_bit(monkeypatch):
    monkeypatch.setattr(tb, "trend_factor", lambda C, *a, **k: np.ones(C.shape[1]))
    assert digest(CANDIDATES["fgw_aligned_trend"].build(ctx())) == PRE_EDIT_SHA256["fgw_aligned"]


def test_both_factor_arms_refuse_samples_without_a_matched_label():
    s = stages()
    for a in s.values():
        a.obs["celltype"] = a.obs["celltype"].astype(str).str.replace("EXE-", "EMB-")
    c = ctx(inputs=s)
    with pytest.raises(ValueError, match="no matched extraembryonic label"):
        CANDIDATES["fgw_early_specialist"].build(c)
    with pytest.raises(ValueError, match="no extraembryonic label held by every sample"):
        CANDIDATES["fgw_trend_render"].build(c)
    with pytest.raises(ValueError, match="no extraembryonic label held by every sample"):
        CANDIDATES["fgw_aligned_trend"].build(c)
