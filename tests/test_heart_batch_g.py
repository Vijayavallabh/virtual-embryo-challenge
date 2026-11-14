"""Heart batch G (campaign/heart-g; vec.trend_flow): `b0_trend` with a per-cell, per-gene factor D inside each moved
lineage, from H10's field (`b0_trend_flow`), its stage-permuted placebo (`b0_trend_flow_perm`), the field on latents
centred per lineage and stage (`b0_trend_flow_pool`) and that field's placebo (`b0_trend_flow_pool_perm`); and the merge
of campaign/gen into main that carries them. Synthetic data only. The count model is row E14's stand-in of
tests/test_latent_extrap.py; torch is not in the main environment, so the flow network (`fit_flow`) is replaced by a
deterministic stand-in that still reads its pairs."""

from __future__ import annotations

import hashlib

import numpy as np
import pytest
from scipy import sparse
from test_b0_trend import BOARD, FLOW_INPUTS, _heart_task, _stage
from test_latent_extrap import _standin

from vec import harness, latent_extrap, t1_scope, t2_latent
from vec import trend_flow as tf
from vec.harness import CANDIDATES

ARMS = {"flow": "b0_trend_flow", "flow_perm": "b0_trend_flow_perm", "flow_pool": "b0_trend_flow_pool",
        "flow_pool_perm": "b0_trend_flow_pool_perm"}
SEED = harness.BUILD_SEED_OFFSET + 1


@pytest.fixture(scope="module")
def three():
    """Three heart-like stages (days 8.5, 8.75, 9.5): CM and Endothelial moved, Untracked, Tiny too small at start."""
    return {"ref": _stage(1600, seed=2, shift=0.0, tiny=60), "start": _stage(1100, seed=3, shift=0.8, tiny=10),
            "earlier": _stage(900, seed=5, shift=1.6, tiny=10)}


def _ctx(inputs, seed: int = SEED) -> harness.BuildContext:
    return harness.BuildContext(task=harness.TaskView.of(_heart_task(inputs=FLOW_INPUTS)), board=BOARD, inputs=inputs,
                                seed=seed)


def _fake_fit_flow(calls: list):
    """A deterministic stand-in for the flow network that reads its pairs: per lineage, the omega-weighted mean pair
    velocity a, mean pair midpoint b and SD s of the pair velocities, v(z) = a - s tanh(z - b) (nonlinear, so the
    centred move depends on the pairs); its checksum is that of (a, b, s)."""

    def fit(pairs, lineages):
        calls.append(pairs)
        L = len(lineages)
        a, b, s = [], [], []
        for li in range(L):
            sel = pairs["lineage"] == li
            w = pairs["omega"][sel][:, None]
            a.append((w * pairs["u"][sel]).sum(0) / w.sum())
            b.append((w * 0.5 * (pairs["z_i"][sel] + pairs["z_j"][sel])).sum(0) / w.sum())
            s.append(np.sqrt((w * (pairs["u"][sel] - a[-1]) ** 2).sum(0) / w.sum()))
        a, b, s = np.array(a), np.array(b), np.array(s)
        sha = hashlib.sha256(a.tobytes() + b.tobytes() + s.tobytes()).hexdigest()
        return latent_extrap.FlowField(lineages=list(lineages),
                                       velocity=lambda z, li: a[li] - s[li] * np.tanh(z - b[li]),
                                       record={"params_sha256": sha, "stand_in": True})

    return fit


@pytest.fixture
def stubs(monkeypatch):
    calls: dict = {"train": [], "fit": []}
    monkeypatch.setattr(t2_latent, "train_nbvae", _standin(calls["train"]))
    monkeypatch.setattr(latent_extrap, "fit_device", lambda: "cpu")
    monkeypatch.setattr(latent_extrap, "fit_flow", _fake_fit_flow(calls["fit"]))
    monkeypatch.setattr(latent_extrap, "_FOLD", {})
    monkeypatch.setattr(latent_extrap, "_FLOW", {})
    monkeypatch.setattr(latent_extrap, "RECORDS", [])
    monkeypatch.setattr(tf, "_FIELDS", {})
    return calls


def _spy(monkeypatch, mod, name):
    seen = []
    real = getattr(mod, name)

    def wrapped(*a, **kw):
        out = real(*a, **kw)
        seen.append({"args": a, "kw": kw, "out": out})
        return out

    monkeypatch.setattr(mod, name, wrapped)
    return seen


def _drawn_types(ctx) -> np.ndarray:
    pool, rows, _, _ = harness.draw_copy(ctx)
    return np.asarray(pool.obs["merged"]).astype(str)[rows]


# ------------------------------------------------------------------------------------------------ the merge
# sha256 of (X float32 bytes, then coords float64 bytes) at build seeds BUILD_SEED_OFFSET + 1 on `three`, with the
# stand-in count model and, for b0_lat_flow, the stand-in field of `_digest_field` (fold_flow replaced; no torch).
# Computed before the merge on each candidate's own branch: b0_trend on main (7f4683d), b0_lat_ratio and b0_lat_flow
# on campaign/gen (459916f), by a scratch script that put that branch's src first on sys.path; the merged tree gives
# the same digests.
PRE_MERGE = {"b0": "1b2d4174129ff6f17192af0cd43eaab95e63070dcc0479d70cfdac7cea71d186",
             "b0_trend": "404f50c201ad7a01fc1afd461aab2a304372ca978e82b96b0b495d6811c1356f",
             "b0_lat_ratio": "8806c11d59d8fab5271b3d029b92241b8ddba96d3bba6768b50d8d5d2292d8fb",
             "b0_lat_flow": "f8541db7160671d5227c19207ee32ad003a0fa5f463e29d92a3844aef513dd89"}


def _digest_field(fold, stages, pool, start, lineages):
    mp = np.asarray(pool.obs["merged"]).astype(str)
    ms = np.asarray(start.obs["merged"]).astype(str)
    gap = stages[-1][1] - stages[-2][1]
    c = np.array([(fold.z_ref[mp == lin].mean(0) - fold.z_start[ms == lin].mean(0)) / gap for lin in lineages])
    zbar = fold.z_ref.mean(0)
    return latent_extrap.FlowField(lineages=list(lineages),
                                   velocity=lambda z, li: c[li] - 0.5 * (np.asarray(z) - zbar), record={"stub": True})


def test_the_merge_leaves_b0_trend_and_the_latent_arms_bit_for_bit(stubs, three, monkeypatch):
    monkeypatch.setattr(latent_extrap, "fold_flow", _digest_field)
    for name, want in PRE_MERGE.items():
        st = CANDIDATES[name].build(_ctx(three))
        h = hashlib.sha256(np.ascontiguousarray(st.X, dtype=np.float32).tobytes())
        h.update(np.ascontiguousarray(st.coords, dtype=np.float64).tobytes())
        assert h.hexdigest() == want, name


# ------------------------------------------------------------------------------------------------ registration
def test_registration_is_b0_lat_flows():
    tasks = {t.key: t for ts_ in harness.load_tasks().values() for t in ts_}
    ref = CANDIDATES["b0_lat_flow"]
    for arm, name in ARMS.items():
        c = CANDIDATES[name]
        assert {k for k, t in tasks.items() if c.applies(t)} == {"T2h/<=E9.5->E10.5", "T2h/<=E8.75->E10.5"}
        assert all(c.applies(t) == ref.applies(t) for t in tasks.values())
        assert c.derived == ref.derived == () and c.params["flow_seed"] == latent_extrap.FLOW_SEED
        assert c.params["scope"] == CANDIDATES["b0_trend"].params["scope"]
        assert ("perm_stream" in c.params) == arm.endswith("_perm")
        if arm.endswith("_perm"):
            assert c.params["perm_stream"] == tf.PERM_STREAM
    assert tf.PERM_STREAM not in (1, 2)  # flow_pairs' streams under FLOW_SEED


# ------------------------------------------------------------------------------------------------ D = 1 is b0_trend
def test_with_d_one_everywhere_the_arm_is_b0_trend_bit_for_bit(three):
    ctx = _ctx(three)
    base = CANDIDATES["b0_trend"].build(ctx)
    one = harness._b0_t2(ctx, trend=True, cell_factors=lambda c, pool, rows, like: np.ones(like.shape))
    np.testing.assert_array_equal(one.X, base.X)
    np.testing.assert_array_equal(one.coords, base.coords)


def test_the_hook_refuses_d_on_an_unshifted_type(three):
    ctx = _ctx(three)

    def bad(c, pool, rows, like):
        D = np.ones(like.shape)
        D[np.asarray(like.obs["b0_type"]).astype(str) == "Tiny", 0] = 2.0
        return D

    with pytest.raises(ValueError, match="not shifted"):
        harness._b0_t2(ctx, trend=True, cell_factors=bad)
    with pytest.raises(ValueError, match="trend=True"):
        harness._b0_t2(ctx, cell_factors=lambda *a: None)


@pytest.mark.parametrize("arm", list(ARMS))
def test_a_move_alike_for_every_cell_of_a_lineage_is_b0_trend_bit_for_bit(stubs, three, monkeypatch, arm):
    """dz the same for all drawn cells of each lineage (a different vector per lineage): du = 0, D = 1, b0_trend."""

    def alike(arm_, fold, stages, pool, start, types, moved, z_own, damping):
        dz = np.zeros_like(z_own)
        for k, lin in enumerate(sorted(moved)):
            dz[np.asarray(types).astype(str) == lin] = 0.3 * (k + 1)
        return dz, {"params_sha256": None, "tau_days": 0.375}

    monkeypatch.setattr(tf, "latent_moves", alike)
    ctx = _ctx(three)
    np.testing.assert_array_equal(CANDIDATES[ARMS[arm]].build(ctx).X, CANDIDATES["b0_trend"].build(ctx).X)
    assert latent_extrap.RECORDS[-1]["check"]["sd_log_d"] == 0.0


def test_a_constant_field_gives_d_one_up_to_rounding(stubs, three, monkeypatch):
    """Through the real RK4 a constant field moves every cell alike up to float rounding of z + h v: D is 1 within
    1e-12 (not bit for bit: the exact case is the test above)."""

    def fold_flow(fold, stages, pool, start, lineages):
        c = np.arange(1, len(lineages) + 1, dtype=np.float64)[:, None] * np.linspace(-1, 1, fold.z_ref.shape[1])
        return latent_extrap.FlowField(lineages=list(lineages), velocity=lambda z, li: c[li].copy(), record={})

    monkeypatch.setattr(latent_extrap, "fold_flow", fold_flow)
    seen = _spy(monkeypatch, tf, "recentre")
    CANDIDATES["b0_trend_flow"].build(_ctx(three))
    np.testing.assert_allclose(seen[0]["out"], 1.0, rtol=0, atol=1e-12)


# ------------------------------------------------------------------------------------------------ every arm
@pytest.mark.parametrize("arm", list(ARMS))
def test_each_arm_keeps_b0_trends_cells_and_recentres(stubs, three, monkeypatch, arm):
    draws = _spy(monkeypatch, harness, "draw_copy")
    rec_seen = _spy(monkeypatch, tf, "recentre_after_depth")
    applied = _spy(monkeypatch, t1_scope, "apply_factors")
    ctx = _ctx(three)
    out = CANDIDATES[ARMS[arm]].build(ctx)
    base = CANDIDATES["b0_trend"].build(ctx)
    # the same rows, coordinates and number
    np.testing.assert_array_equal(draws[0]["out"][1], draws[1]["out"][1])
    np.testing.assert_array_equal(out.coords, base.coords)
    assert out.n == base.n == BOARD.max_cells
    merged = _drawn_types(ctx)
    # lineages b0_trend does not move, and Untracked cells: D = 1 and their counts b0_trend's
    D = rec_seen[0]["out"][0]
    for kept in ("Untracked", "Tiny"):
        m = merged == kept
        assert m.any() and (D[m] == 1.0).all()
        np.testing.assert_array_equal(out.X[m], base.X[m])
    # D acts inside CM and leaves each gene's expected total over each moved lineage, after b0_trend's depth rule, as
    # it was: from the factors and the per-cell depth constants apply_factors itself used in the arm's call and in
    # b0_trend's
    pool, rows, _, _ = harness.draw_copy(ctx)
    C = pool.layers["counts"][rows].toarray().astype(np.float64)
    assert len(applied) == 2
    arm_call, base_call = applied[0], applied[1]
    f_arm = arm_call["args"][1][arm_call["args"][2]] * arm_call["out"][1][:, None]
    f_base = base_call["args"][1][base_call["args"][2]] * base_call["out"][1][:, None]
    for lin in ("CM", "Endothelial"):
        m = merged == lin
        tot = C[m].sum(0)
        has = tot > 0
        np.testing.assert_allclose((C[m] * f_arm[m]).sum(0)[has], (C[m] * f_base[m]).sum(0)[has], rtol=1e-9, atol=0)
        assert (D[m][:, ~has] == 1.0).all()
    assert np.abs(np.log(D[merged == "CM"])).max() > 0
    assert not np.array_equal(out.X[merged == "CM"], base.X[merged == "CM"])
    rec = latent_extrap.RECORDS[-1]
    assert rec["candidate"] == ARMS[arm] and rec["check"]["max_rel_departure"] < 1e-9
    # the one-step G (the count-weighted mean of r) departs after the depth rule; the iteration removes it
    assert rec["check"]["max_rel_departure_one_step"] > 1e-6 and rec["check"]["recentre_iterations"] > 0
    assert set(rec["check"]["within_share"]) == {"CM", "Endothelial"}
    assert all(0 < v <= 1 for v in rec["check"]["within_share"].values())
    assert rec["check"]["lineages_moved"] == ["CM", "Endothelial"] and rec["check"]["sd_log_d"] > 0
    assert rec["check"]["latent_sha256"] == latent_extrap._FOLD["fold"]["fold"].record["latent_sha256"]
    # the same seed, the same output
    np.testing.assert_array_equal(CANDIDATES[ARMS[arm]].build(ctx).X, out.X)


def test_one_count_model_and_one_field_per_arm_for_every_seed(stubs, three, capsys):
    checks = {}
    for seed in (SEED, SEED + 1):
        for name in ARMS.values():
            CANDIDATES[name].build(_ctx(three, seed=seed))
        CANDIDATES["b0_lat_flow"].build(_ctx(three, seed=seed))
    assert len(stubs["train"]) == 1  # one count model for the four arms and b0_lat_flow
    assert len(stubs["fit"]) == 4  # one field per arm; b0_trend_flow's is b0_lat_flow's
    for r in latent_extrap.RECORDS:
        if r["candidate"] in ARMS.values():
            checks.setdefault(r["candidate"], set()).add(r["check"]["params_sha256"])
    assert all(len(v) == 1 for v in checks.values()) and len({next(iter(v)) for v in checks.values()}) == 4
    lines = [x for x in capsys.readouterr().out.splitlines() if x.startswith("check ")]
    assert len(lines) == 8 and all("params_sha256=" in x and "latent_sha256=" in x for x in lines)


# ------------------------------------------------------------------------------------------------ the placebo
def test_the_permutation_keeps_each_stages_cells_per_lineage_and_moves_with_its_stream_only():
    rng = np.random.default_rng(0)
    labs = [rng.choice(["A", "B", "C"], size=n) for n in (300, 250, 400)]
    zs = [rng.normal(size=(len(x), 4)) + k for k, x in enumerate(labs)]
    a = tf.stage_permutation(labs, [latent_extrap.FLOW_SEED, tf.PERM_STREAM])
    np.testing.assert_array_equal(a, tf.stage_permutation(labs, [latent_extrap.FLOW_SEED, tf.PERM_STREAM]))
    assert not np.array_equal(a, tf.stage_permutation(labs, [latent_extrap.FLOW_SEED, tf.PERM_STREAM + 1]))
    lab = np.concatenate(labs)
    orig = np.concatenate([np.full(len(x), k) for k, x in enumerate(labs)])
    assert (a != orig).any()
    for lin in "ABC":
        for k in range(3):
            assert ((a == k) & (lab == lin)).sum() == ((orig == k) & (lab == lin)).sum()
    pz, pl = tf.permute_stages(zs, labs, [latent_extrap.FLOW_SEED, tf.PERM_STREAM])
    for k in range(3):
        assert sorted(pl[k]) == sorted(labs[k]) and len(pz[k]) == len(zs[k])


def test_the_placebos_field_ignores_the_build_seed_and_differs_from_the_real_one(stubs, three, monkeypatch):
    pairs_seen = _spy(monkeypatch, latent_extrap, "flow_pairs")
    for seed in (SEED, SEED + 7):
        CANDIDATES["b0_trend_flow_perm"].build(_ctx(three, seed=seed))
    tf._FIELDS.clear()
    CANDIDATES["b0_trend_flow_perm"].build(_ctx(three, seed=SEED + 7))
    CANDIDATES["b0_trend_flow"].build(_ctx(three))
    assert len(pairs_seen) == 3  # once per fit: the perm field (twice, after a cache clear) and the real one
    p1, p2, real = (s["args"] for s in pairs_seen)
    for k in range(len(p1[0])):
        np.testing.assert_array_equal(p1[0][k], p2[0][k])  # the build seed does not enter the permutation
    assert not np.array_equal(p1[0][0], real[0][0])  # the real field's earliest stage is not the permuted one
    # the permuted stages hold the same cells per lineage as the real ones
    for k in range(len(real[1])):
        assert sorted(p1[1][k]) == sorted(real[1][k])
    # and they are the real stages with each lineage's true stage means kept and only the residuals permuted
    want, wlab, _ = tf.permute_residuals(real[0], real[1], [latent_extrap.FLOW_SEED, tf.PERM_STREAM], add_means=True)
    for k in range(len(real[1])):
        np.testing.assert_array_equal(p1[0][k], want[k])
        np.testing.assert_array_equal(p1[1][k], wlab[k])


# ------------------------------------------------------------------------------------------------ the pooled form
@pytest.mark.parametrize("arm", ["flow_pool", "flow_pool_perm"])
def test_the_pooled_forms_u_is_centred_per_lineage_and_stage_and_its_pairs_come_from_u(stubs, three, monkeypatch,
                                                                                       arm):
    pairs_seen = _spy(monkeypatch, latent_extrap, "flow_pairs")
    moves = _spy(monkeypatch, latent_extrap, "integrate")
    ctx = _ctx(three)
    CANDIDATES[ARMS[arm]].build(ctx)
    zs_u, labs = pairs_seen[0]["args"][0], pairs_seen[0]["args"][1]
    if arm == "flow_pool":
        for u, lab in zip(zs_u, labs, strict=True):
            for lin in set(lab):
                np.testing.assert_allclose(u[lab == lin].mean(0), 0.0, atol=1e-12)
    # the u's are the stages' latents less their (lineage, stage) means, the stages permuted first for the placebo
    fold = latent_extrap._FOLD["fold"]["fold"]
    stages = [(r, d, ctx.inputs[r]) for r, d in latent_extrap.flow_stages(ctx.task)]
    zs, lb, _ = latent_extrap._stage_latents(fold, stages, ctx.inputs["ref"], ctx.inputs["start"])
    if arm.endswith("_perm"):
        # centred on the true stages first, then the residuals permuted (no re-centring)
        want, _, _ = tf.permute_residuals(zs, lb, [latent_extrap.FLOW_SEED, tf.PERM_STREAM], add_means=False)
        cu, _ = tf.centre_by_stage(zs, lb)
        new = tf.stage_permutation(lb, [latent_extrap.FLOW_SEED, tf.PERM_STREAM])
        for k, a in enumerate(want):
            np.testing.assert_array_equal(a, np.concatenate(cu)[new == k])
    else:
        want, _ = tf.centre_by_stage(zs, lb)
    for a, b in zip(zs_u, want, strict=True):
        np.testing.assert_array_equal(a, b)
    # the pairs' velocities are differences of u, and the drawn cells start from their own u
    p = pairs_seen[0]["out"][0]
    assert np.abs(p["z_i"].mean(0)).max() < 0.5 * np.abs(np.concatenate(zs)).mean()
    x0 = moves[0]["args"][1]
    pool, rows, _, _ = harness.draw_copy(ctx)
    merged = np.asarray(pool.obs["merged"]).astype(str)[rows]
    m = np.isin(merged, ["CM", "Endothelial"])
    ref_lab = lb[-1]  # the true reference stage (lb, zs are not permuted): its lineage means, for both arms
    means = {lin: zs[-1][ref_lab == lin].mean(0) for lin in ("CM", "Endothelial")}
    np.testing.assert_allclose(x0, fold.z_ref[rows][m] - np.array([means[t] for t in merged[m]]), atol=1e-12)


def test_centre_moves_and_recentre_on_their_own():
    types = np.array(["A"] * 4 + ["B"] * 3 + ["U"] * 2)
    dz = np.vstack([np.arange(4)[:, None] * [1.0, 2.0], np.tile([5.0, -1.0], (3, 1)), np.ones((2, 2))])
    du = tf.centre_moves(dz, types, ["A", "B"])
    np.testing.assert_allclose(du[:4].mean(0), 0.0, atol=1e-15)
    assert (du[4:] == 0.0).all()  # B moved alike, U not a moved lineage
    rng = np.random.default_rng(3)
    C = sparse.csr_matrix(rng.poisson(2.0, size=(9, 5)) * (rng.random((9, 5)) > 0.3))
    C = C.toarray()
    C[:4, 4] = 0  # a gene with no count in A
    r = np.exp(rng.normal(scale=0.3, size=(9, 5)))
    D = tf.recentre(r, C, types, ["A", "B"])
    for lin in "AB":
        m = types == lin
        tot = C[m].sum(0)
        np.testing.assert_allclose((C[m] * D[m]).sum(0)[tot > 0], tot[tot > 0], rtol=1e-12)
    assert (D[:4, 4] == 1.0).all() and (D[types == "U"] == 1.0).all()


def test_momentum_shift_without_cell_factors_is_unchanged():
    """The hook's default leaves momentum_shift's code path as it was (the pinned pre-change hash of
    tests/test_b0_trend.py covers its output); here: D = 1 through the hook equals no hook, bit for bit."""
    from test_b0_trend import _adata

    rng = np.random.default_rng(21)
    genes = [f"g{i}" for i in range(50)]
    mu = rng.gamma(2.0, 2.0, 50)
    start = _adata(rng.poisson(mu, (300, 50)), ["A"] * 200 + ["B"] * 100, genes, "s")
    ref = _adata(rng.poisson(mu * np.linspace(0.5, 2.0, 50), (400, 50)), ["A"] * 250 + ["B"] * 120 + ["C"] * 30,
                 genes, "r")
    like = ref[np.arange(0, 400, 3)].copy()
    a, _ = t1_scope.momentum_shift(start, ref, like, 0.5, key="b0_type", seed=[7, 3])
    b, _ = t1_scope.momentum_shift(start, ref, like, 0.5, key="b0_type", seed=[7, 3], cell_factors=np.ones(like.shape))
    np.testing.assert_array_equal(a.layers["counts"].toarray(), b.layers["counts"].toarray())


# ------------------------------------------------------------------------------------------------ review fixes
def test_the_placebo_keeps_each_lineages_mean_development_and_permutes_only_residuals():
    """flow_perm's stages: residual about the true (lineage, stage) mean, permuted, the destination stage's true mean
    added back. With no within-lineage spread the placebo's stages are the true stages: only within-lineage content is
    permuted (a plain stage permutation mixes the stages' means)."""
    rng = np.random.default_rng(4)
    labs = [rng.choice(["A", "B"], size=n) for n in (120, 90, 150)]
    mu = {(k, lin): rng.normal(size=3) * 3 for k in range(3) for lin in "AB"}
    flat = [np.array([mu[(k, lin)] for lin in lab]) for k, lab in enumerate(labs)]
    seed = [latent_extrap.FLOW_SEED, tf.PERM_STREAM]
    pz, pl, _ = tf.permute_residuals(flat, labs, seed, add_means=True)
    for k in range(3):
        for lin in "AB":
            m = pl[k] == lin
            assert m.sum() == (labs[k] == lin).sum()
            np.testing.assert_allclose(pz[k][m], np.tile(mu[(k, lin)], (m.sum(), 1)), atol=1e-12)
    old, _ = tf.permute_stages(flat, labs, seed)
    assert not np.allclose(np.concatenate(old), np.concatenate(pz))  # the stage permutation mixes the means
    # with spread: each cell's residual is carried to its new stage, the means the true stages'
    zs = [f + rng.normal(size=f.shape) for f in flat]
    pz, pl, means = tf.permute_residuals(zs, labs, seed, add_means=True)
    cu, _ = tf.centre_by_stage(zs, labs)
    new = tf.stage_permutation(labs, seed)
    U, lab = np.concatenate(cu), np.concatenate(labs)
    for k in range(3):
        want = U[new == k] + np.array([means[k][x] for x in lab[new == k]])
        np.testing.assert_array_equal(pz[k], want)
        np.testing.assert_array_equal(pl[k], lab[new == k])
    pu, _, _ = tf.permute_residuals(zs, labs, seed, add_means=False)
    for k in range(3):
        np.testing.assert_array_equal(pu[k], U[new == k])


def test_the_recentring_holds_after_apply_factors_depth_rule():
    """sum_i c F D s'_i = sum_i c F s_i per lineage and gene, s and s' the per-cell constants apply_factors returns;
    the one-step G (the count-weighted mean of r) misses it by far more than the tolerance."""
    from vec.core import apply_factors

    rng = np.random.default_rng(11)
    types = np.array(["A"] * 60 + ["B"] * 40 + ["U"] * 20)
    C = rng.poisson(rng.gamma(1.0, 3.0, 30), size=(120, 30)) * (rng.random((120, 30)) > 0.2)
    C[:60, 7] = 0  # a gene with no count in A
    held = np.zeros(30, dtype=bool)
    held[3] = True
    F = np.exp(rng.normal(scale=0.5, size=(3, 30)))
    F[:, held] = 1.0
    tidx = np.array([{"A": 0, "B": 1, "U": 2}[t] for t in types])
    Fc = F[tidx]
    r = np.exp(rng.normal(scale=0.4, size=C.shape))
    D, info = tf.recentre_after_depth(r, C, Fc, held, types, ["A", "B"])
    assert (D[types == "U"] == 1.0).all() and (D[:60, 7] == 1.0).all() and (D[:, held] == 1.0).all()
    _, s0, _ = apply_factors(C, F, tidx, held, seed=0)
    _, s1, _ = apply_factors(C, Fc * D, np.arange(len(C)), held, seed=0)
    one = np.where(held, 1.0, tf.recentre(r, C, types, ["A", "B"]))
    _, s2, _ = apply_factors(C, Fc * one, np.arange(len(C)), held, seed=0)
    worst = 0.0
    for lin in "AB":
        m = types == lin
        without = (C[m] * Fc[m] * s0[m, None]).sum(0)
        ok = (without > 0) & ~held
        np.testing.assert_allclose((C[m] * Fc[m] * D[m] * s1[m, None]).sum(0)[ok], without[ok], rtol=1e-9, atol=0)
        got = (C[m] * Fc[m] * one[m] * s2[m, None]).sum(0)[ok]
        worst = max(worst, float(np.abs(got / without[ok] - 1).max()))
    assert worst > 1e-3
    assert info["max_rel_departure"] <= tf.RECENTRE_TOL
    assert info["max_rel_departure_one_step"] == pytest.approx(worst, rel=1e-6)


def test_momentum_shift_calls_a_callable_hook_with_the_trend_factors_and_held_mask():
    from test_b0_trend import _adata

    rng = np.random.default_rng(21)
    genes = [f"g{i}" for i in range(50)]
    mu = rng.gamma(2.0, 2.0, 50)
    start = _adata(rng.poisson(mu, (300, 50)), ["A"] * 200 + ["B"] * 100, genes, "s")
    ref = _adata(rng.poisson(mu * np.linspace(0.5, 2.0, 50), (400, 50)), ["A"] * 250 + ["B"] * 120 + ["C"] * 30,
                 genes, "r")
    like = ref[np.arange(0, 400, 3)].copy()
    seen = []

    def hook(Fc, held):
        seen.append((Fc.copy(), held.copy()))
        return np.ones(Fc.shape)

    a, _ = t1_scope.momentum_shift(start, ref, like, 0.5, ("g4",), key="b0_type", seed=[7, 3])
    b, _ = t1_scope.momentum_shift(start, ref, like, 0.5, ("g4",), key="b0_type", seed=[7, 3], cell_factors=hook)
    np.testing.assert_array_equal(a.layers["counts"].toarray(), b.layers["counts"].toarray())
    Fc, held = seen[0]
    assert Fc.shape == like.shape and held.tolist() == [g == "g4" for g in genes] and (Fc[:, 4] == 1.0).all()
    t = np.asarray(like.obs["b0_type"]).astype(str)
    assert (Fc[t == "C"] == 1.0).all() and not (Fc[t == "A"] == 1.0).all()


@pytest.mark.parametrize("arm", list(ARMS))
def test_every_arm_moves_over_damping_times_the_last_gap(stubs, three, monkeypatch, arm):
    """tau = B0_DAMPING x (reference day - start day), in the record and in the integrator's call."""
    moves = _spy(monkeypatch, latent_extrap, "integrate")
    ctx = _ctx(three)
    CANDIDATES[ARMS[arm]].build(ctx)
    days = dict(latent_extrap.flow_stages(ctx.task))
    want = harness.B0_DAMPING * (days[ctx.task.reference_role] - days["start"])
    assert latent_extrap.RECORDS[-1]["check"]["tau_days"] == pytest.approx(want, rel=1e-15)
    assert moves and all(float(m["args"][3]) == pytest.approx(want, rel=1e-15) for m in moves)
