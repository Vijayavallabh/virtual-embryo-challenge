"""The embryo trend rendering `fgw_trend_render`, ported alone from campaign/emb-trend2 at 51b0672 (port/emb-s1-trend),
and its missing-gene rule (vec.t2_build.TREND_COMMON_RULE), on synthetic labelled stages built here. Nothing is read
from data/, nothing is scored, and no prediction is compared with any target.

The trend parts of 51b0672's tests/test_emb_expr.py, unchanged where they test the trend arm: the arm keeps fgw's drawn
cells, positions and cell number, and with its factor at 1 it is fgw bit for bit. Its built output on `ctx()` is pinned
at the checksum 51b0672's own code gives (TREND_AT_51B0672), and every pre-existing candidate that builds on the
synthetic contexts is pinned at main 4f39a11, before the port (PINS_AT_MAIN_4F39A11)."""

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

# sha256 of (X float32, coords float32) fgw_trend_render builds on `ctx()`, computed with the src of 51b0672 (git
# archive of that commit, OPENBLAS/OMP/MKL threads 4); 51b0672's own test pins the same value
TREND_AT_51B0672 = "a380d451db01db4d3f8ab02a6bea78eb5ed803708d620d84cecd182a75dfb7b7"
# sha256 of (X float32, coords float32) of every pre-existing candidate that builds on the synthetic contexts, computed
# with main at 4f39a11 before the port (OPENBLAS/OMP/MKL threads 4). fgw_latent needs scvi-tools, which this environment
# lacks; b0 and b0_herd need the lineage map's `merged` column, which the synthetic stages lack (refusals pinned below).
PINS_AT_MAIN_4F39A11 = {
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
NEW_ON_PORT = ("fgw_trend_render",)  # registered by the port


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


def digest_X(X, coords) -> str:
    h = hashlib.sha256(np.ascontiguousarray(np.clip(X, 0.0, None), dtype=np.float32).tobytes())
    h.update(np.ascontiguousarray(coords, dtype=np.float32).tobytes())
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


def _samples(c):
    return [(a, c.task.role_days[r]) for r, a in c.inputs.items()]


@pytest.fixture(scope="module")
def built():
    c = ctx()
    names = ("fgw", "fgw_trend_render")
    return c, {n: CANDIDATES[n].build(c) for n in names}


# ------------------------------------------------------------------------------------- existing candidates unchanged
@pytest.fixture(scope="module")
def pinned_builds():
    s = stages()
    return {roles: {n: digest(CANDIDATES[n].build(ctx(inputs=s, roles=roles))) for n in pins}
            for roles, pins in PINS_AT_MAIN_4F39A11.items()}


@pytest.mark.parametrize("roles, name", [(r, n) for r, pins in PINS_AT_MAIN_4F39A11.items() for n in sorted(pins)])
def test_every_pre_existing_candidate_is_what_it_was_on_main(pinned_builds, roles, name):
    assert pinned_builds[roles][name] == PINS_AT_MAIN_4F39A11[roles][name]


def test_the_pins_cover_every_pre_existing_candidate_that_builds_here():
    """Every candidate applicable to a synthetic context is pinned, refused as pinned, or new on the port."""
    for roles, pins in PINS_AT_MAIN_4F39A11.items():
        t = task(roles)
        here = {n for n, c in CANDIDATES.items() if n != harness.FLOOR and c.applies(t)}
        assert here - set(pins) <= {"b0", "b0_herd", "fgw_latent", *NEW_ON_PORT}
    c = ctx()
    for name in ("b0", "b0_herd"):
        with pytest.raises(KeyError, match="merged"):
            CANDIDATES[name].build(c)


def test_only_the_trend_arm_is_ported():
    for name in ("fgw_aligned_hi", "fgw_early_specialist", "fgw_aligned_trend"):
        assert name not in CANDIDATES
    assert not hasattr(harness, "TREND_SAMPLE_PREFIX") and not hasattr(tb, "EARLY_RULE")


# ------------------------------------------------------------------------------------------------ the ported arm
def test_trend_render_is_what_51b0672_built(built):
    assert digest(built[1]["fgw_trend_render"]) == TREND_AT_51B0672


def test_registration_and_where_the_arm_applies():
    tasks = {t.key: t for ts in load_tasks().values() for t in ts}
    counted = (tasks["T2e/E7.25+E8.0->E7.5"], tasks["T2e/E6.75+E7.5->E7.25"])
    heart = tasks["T2h/<=E9.5->E10.5"]
    c = CANDIDATES["fgw_trend_render"]
    assert c.name == "fgw_trend_render" and c.tasks == {"T2"} and not c.tuned_on and not c.externals and not c.derived
    assert all(c.applies(t) for t in counted) and not c.applies(heart)
    assert c.params["scope"] == "whole" and c.params["rule"] == tb.TREND_RENDER_RULE
    # the trend needs a third sample: on a two-bracket task it does not apply
    assert not c.applies(task(roles=("before", "after")))


def test_the_arm_keeps_fgws_cells_positions_and_cell_number(built):
    _, b = built
    assert b["fgw_trend_render"].n == b["fgw"].n == BOARD.max_cells
    np.testing.assert_array_equal(b["fgw_trend_render"].coords, b["fgw"].coords)
    assert not np.array_equal(b["fgw_trend_render"].X, b["fgw"].X)  # the expression is the arm's one idea


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
                                          _samples(c), c.task.target_day, [c.seed, 7])
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
    fgw = PINS_AT_MAIN_4F39A11[("earlier", "before", "after")]["fgw"]
    assert digest(CANDIDATES["fgw_trend_render"].build(c)) == fgw


def test_trend_render_refuses_brackets_outside_its_samples(built):
    c, _ = built
    res, *_ = _interp(c)
    with pytest.raises(ValueError, match="among the trend's samples"):
        tb.trend_rendered_expression(c.inputs["before"], c.inputs["after"], res, list(BOARD.genes),
                                     [(c.inputs["earlier"], 6.75), (c.inputs["before"], 7.25)], 7.5, 0)


def test_the_arm_refuses_samples_without_a_matched_label():
    s = stages()
    for a in s.values():
        a.obs["celltype"] = a.obs["celltype"].astype(str).str.replace("EXE-", "EMB-")
    with pytest.raises(ValueError, match="no extraembryonic label held by every sample"):
        CANDIDATES["fgw_trend_render"].build(ctx(inputs=s))


# ------------------------------------------------------------------------------- the missing-gene rule (main session)
def test_with_every_gene_common_the_rule_is_the_ported_code_bit_for_bit(built):
    """trend_rendered_expression_common on samples that hold the whole panel: the ported function's output exactly
    (and so the 51b0672 checksum through the candidate's route), with 0 genes left unrendered."""
    c, b = built
    res, *_ = _interp(c)
    args = (c.inputs["before"], c.inputs["after"], res, list(BOARD.genes), _samples(c), c.task.target_day, [c.seed, 7])
    X0, rec0 = tb.trend_rendered_expression(*args)
    X1, rec1 = tb.trend_rendered_expression_common(*args)
    np.testing.assert_array_equal(X0, X1)
    assert digest_X(X1, res.coords) == digest(b["fgw_trend_render"]) == TREND_AT_51B0672
    assert rec1["panel_genes_unrendered"] == 0 and rec1["rule"] == tb.TREND_COMMON_RULE
    assert {k: v for k, v in rec1.items() if k not in ("rule", "panel_genes_unrendered")} == \
        {k: v for k, v in rec0.items() if k != "rule"}
    assert tb.TREND_COMMON_RULE.startswith(tb.TREND_RENDER_RULE)
    assert tb.trend_common_genes(BOARD.genes, _samples(c)) == list(BOARD.genes)


def test_a_sample_lacking_two_panel_genes_leaves_them_as_fgw_has_them(built):
    """The earlier sample lacks g3 and g9: their columns are fgw's mixture rows exactly, every other column is the
    ported code on a panel made of the other 22 genes (its proportions, rescale and top-up over those genes), and the
    record says 2."""
    c, b = built
    res, Eb, Ea, _ = _interp(c)
    lacking = ["g3", "g9"]
    common = [g for g in PANEL if g not in lacking]
    earlier = c.inputs["earlier"][:, [g for g in c.inputs["earlier"].var_names if g not in lacking]].copy()
    samples = [(earlier, 6.75), (c.inputs["before"], 7.25), (c.inputs["after"], 8.0)]
    X, rec = tb.trend_rendered_expression_common(c.inputs["before"], c.inputs["after"], res, PANEL, samples, 7.5,
                                                 [c.seed, 7])
    assert rec["panel_genes_unrendered"] == 2 and tb.trend_common_genes(PANEL, samples) == common
    mix = np.empty_like(X)
    mix[res.stage == 0], mix[res.stage == 1] = Eb[res.rows[res.stage == 0]], Ea[res.rows[res.stage == 1]]
    lack = [PANEL.index(g) for g in lacking]
    np.testing.assert_array_equal(X[:, lack], mix[:, lack])
    sub, _ = tb.trend_rendered_expression(c.inputs["before"], c.inputs["after"], res, common, samples, 7.5,
                                          [c.seed, 7])
    keep = [PANEL.index(g) for g in common]
    np.testing.assert_array_equal(X[:, keep], sub)
    assert not np.array_equal(X[:, keep], mix[:, keep])
    # the factors are the existing code's on the panel without the two genes
    names = [np.asarray(a.obs["celltype"]).astype(str) for a, _ in samples]
    counts = [sparse.csr_matrix(a[:, common].layers["counts"], dtype=np.int64) for a, _ in samples]
    labels = tb.trend_labels([names[1], names[0], names[2]], np.asarray(c.inputs["before"].obs["group"]).astype(str))
    phi = tb.trend_factor(counts[1], names[1], list(zip(counts, names, strict=True)), [6.75, 7.25, 8.0], 7.5, labels)
    assert phi.shape == (22,) and np.any(phi != 1.0)
    assert "g3" not in str(rec) and "g9" not in str(rec)
    # a bracket must hold the whole panel
    with pytest.raises(ValueError, match="a bracket lacks panel genes"):
        tb.trend_rendered_expression_common(c.inputs["before"][:, common].copy(), c.inputs["after"], res, PANEL,
                                            samples, 7.5, 0)
