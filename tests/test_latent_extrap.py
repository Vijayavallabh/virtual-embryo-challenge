"""vec.latent_extrap and the candidates `b0_lat_dec` and `b0_lat_ratio` (B0 moved in a count model's latent space, the
extrapolation counterpart of embryo row E21), on synthetic T1 and heart stages only. A deterministic stand-in replaces
the scvi-tools training (torch and scvi-tools live in the `gpu` dependency group); the one test of the real training,
on the CPU with a few epochs, is skipped where they are not installed."""

from __future__ import annotations

import inspect
import json
import subprocess
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest
from scipy import sparse
from test_h7_gmom import BOARD as HEART_BOARD
from test_h7_gmom import _ctx as heart_ctx
from test_h7_gmom import _heart_task
from test_h7_gmom import _stage as heart_stage
from test_t1_scope import GENES as T1_GENES
from test_t1_scope import HELD
from test_t1_state_transport import _gmom_ctx as t1_ctx
from test_t1_state_transport import states_stage

from vec import harness, latent_extrap, t1_scope, t2_latent
from vec import t1_state_transport as ts
from vec.core import counts_of
from vec.harness import CANDIDATES
from vec.prep import log1p_cp10k

REPO = Path(__file__).resolve().parents[1]
ARMS = ("b0_lat_dec", "b0_lat_ratio")


def _standin(calls: list):
    """A stand-in for row E14's training: the latent is each cell's log gene proportions (smoothed), the decoder the
    softmax of the latent times the library, a large inverse dispersion. Deterministic; records its calls."""

    def train(counts, seed=t2_latent.TRAIN_SEED, *, device=None, history=None, **kw):
        C = sparse.csr_matrix(counts).toarray().astype(np.float64)
        calls.append({"shape": C.shape, "device": device, "seed": seed})
        lib = C.sum(1, keepdims=True)
        latent = np.log((C + 0.5) / (lib + 0.5 * C.shape[1]))

        def decode(z, library):
            z = np.asarray(z, dtype=np.float64)
            e = np.exp(z - z.max(1, keepdims=True))
            return np.asarray(library, dtype=np.float64)[:, None] * e / e.sum(1, keepdims=True), \
                np.full(C.shape[1], 200.0)

        def encode(counts):
            X = sparse.csr_matrix(counts).toarray().astype(np.float64)
            return np.log((X + 0.5) / (X.sum(1, keepdims=True) + 0.5 * X.shape[1]))

        if history is not None:
            history["train_loss"] = 1.0
        return t2_latent.NBVAE(latent=latent, decode=decode, record={"stand_in": True}, encode=encode)

    return train


@pytest.fixture
def standin(monkeypatch):
    calls: list = []
    monkeypatch.setattr(t2_latent, "train_nbvae", _standin(calls))
    monkeypatch.setattr(latent_extrap, "fit_device", lambda: "cpu")
    monkeypatch.setattr(latent_extrap, "_FOLD", {})
    monkeypatch.setattr(latent_extrap, "_FLOW", {})
    monkeypatch.setattr(latent_extrap, "RECORDS", [])
    return calls


@pytest.fixture(scope="module")
def heart_inputs():
    return {"ref": heart_stage(1600, seed=2, shift=0.0, tiny=60),
            "start": heart_stage(1100, seed=3, shift=0.8, tiny=10)}


@pytest.fixture(scope="module")
def t1_inputs():
    """E9.5-like and E8.5-like stages (the config's per-file rules apply) with a category B0 does not move for each of
    its two reasons: blood below B0's minimum of cells at the start stage, and epicardium, which E8.5 has only through
    its stand-in label JCF."""
    ref = states_stage(7, {"V-CM": 1200, "Endocardium": 350, "Blood": 100, "Proepicardium": 80, "Foregut": 300},
                       later=True)
    start = states_stage(8, {"V-CM": 700, "Endothelium": 250, "Blood": 10, "JCF": 80, "Foregut": 170})
    return ref, start


def _spy(monkeypatch, name):
    """Wrap latent_extrap.<name>, recording each call's arguments and result."""
    seen = []
    real = getattr(latent_extrap, name)

    def wrapped(*a, **kw):
        out = real(*a, **kw)
        seen.append({"args": a, "kw": kw, "out": out})
        return out

    monkeypatch.setattr(latent_extrap, name, wrapped)
    return seen


# ------------------------------------------------------------------------------------------------ constants and scope
def test_the_constants_are_b0s_t11s_and_row_e14s():
    assert latent_extrap.B0_MIN_CELLS == inspect.signature(t1_scope.momentum_shift).parameters["min_cells"].default
    assert latent_extrap.CLIP == ts.CLIP_AT
    assert len({latent_extrap.RATIO_STREAM, latent_extrap.DECODE_STREAM, t2_latent.EXPRESSION_STREAM, 1, 3, 5}) == 6
    # nothing of the count model is set here: train_nbvae is called with row E14's seed and its own defaults
    src = inspect.getsource(latent_extrap.fold_model)
    assert "t2_latent.TRAIN_SEED" in src and "max_epochs" not in src and "n_latent" not in src


def test_the_candidates_apply_to_extrapolation_tasks_with_a_start_and_to_no_other():
    tasks = {t.key: t for ts_ in harness.load_tasks().values() for t in ts_}
    want = {"T1/E9.5->E10.5", "T1/E9.5->E10.5:chamber", "T2h/<=E9.5->E10.5", "T2h/<=E8.75->E10.5"}
    for name in ARMS:
        c = CANDIDATES[name]
        assert {k for k, t in tasks.items() if c.applies(t)} == want
        assert not c.applies(tasks["T1/E8.5->E9.5"])  # no start stage
        assert not c.applies(_heart_task(inputs={"ref": "T2/E9.5.h5ad"}))
        assert c.derived == ("t1_held_flat",) and c.params["damping"] == harness.B0_DAMPING
        assert c.params["scope"] == CANDIDATES["b0"].params["scope"]
        # B0's own momentum arms on the two boards apply on the same tasks
        for t in tasks.values():
            if t.task == "T1":
                assert c.applies(t) == CANDIDATES["g_mom"].applies(t)
            elif t.kind == "extrap":
                assert c.applies(t) == CANDIDATES["b0_cfm"].applies(t)


def test_the_module_and_the_harness_import_without_torch_or_scvi():
    code = ("import sys; sys.path.insert(0, 'src'); import vec.latent_extrap, vec.harness; "
            "assert 'torch' not in sys.modules and 'scvi' not in sys.modules")
    subprocess.run([sys.executable, "-I", "-B", "-c", code], check=True, cwd=REPO, capture_output=True)


# ------------------------------------------------------------------------------------------------ (a) zero displacement
def test_heart_ratio_with_zero_displacement_is_b0_at_damping_zero(standin, heart_inputs, monkeypatch):
    ctx = heart_ctx(heart_inputs)
    monkeypatch.setattr(harness, "B0_DAMPING", 0.0)
    b0_flat = CANDIDATES["b0"].build(ctx)
    ratio = CANDIDATES["b0_lat_ratio"].build(ctx)
    np.testing.assert_array_equal(ratio.X, b0_flat.X)
    np.testing.assert_array_equal(ratio.coords, b0_flat.coords)
    # and that is the drawn real cells' own expression (the fixture's X up to its own float rounding)
    pool, rows, _, _ = harness.draw_copy(ctx)
    np.testing.assert_allclose(ratio.X, pool.X[rows].toarray(), atol=2e-6)
    _, _, _, rec = harness.b0_lat_population(ctx, "ratio")
    assert rec["latent"]["mean_abs_latent_step"] == 0.0 and rec["expression"]["n_cells_moved"] == 0


def test_t1_ratio_with_zero_displacement_is_b0_at_damping_zero_and_the_drawn_cells(standin, t1_inputs, monkeypatch):
    ref, start = t1_inputs
    ctx = t1_ctx(ref, start)
    monkeypatch.setattr(harness, "B0_DAMPING", 0.0)
    ratio = CANDIDATES["b0_lat_ratio"].build(ctx)
    np.testing.assert_array_equal(ratio.X, CANDIDATES["b0"].build(ctx).X)
    np.testing.assert_array_equal(ratio.X, CANDIDATES["b0_scope"].build(ctx).X)
    # the drawn cells' released rows, exactly (the fixture's X is log1p CP10k of its counts)
    pop = ts.b0_population(harness._shallow(ref), harness._shallow(start),
                           harness.t1_chamber_shares(harness._shallow(ref), "E9.5_RNA"), ctx.board.max_cells,
                           ctx.seed, HELD, ref_source="E9.5_RNA", start_source="E8.5_RNA", damping=0.0)
    rows = ref.obs_names.get_indexer(pop.obs["source_obs"].astype(str))
    np.testing.assert_array_equal(ratio.X, ref.X[rows].toarray())


# ------------------------------------------------------------------------------------------------ (b) the decoded arm
def test_heart_decoded_arm_is_b0s_rows_and_coordinates_at_the_real_libraries(standin, heart_inputs, monkeypatch):
    ctx = heart_ctx(heart_inputs)
    seen = _spy(monkeypatch, "decoded_counts")
    lat = _spy(monkeypatch, "moved_latent")
    b0 = CANDIDATES["b0"].build(ctx)
    dec = CANDIDATES["b0_lat_dec"].build(ctx)
    assert dec.X.shape == b0.X.shape == (HEART_BOARD.max_cells, len(HEART_BOARD.genes))
    assert np.isfinite(dec.X).all() and (dec.X >= 0).all() and (dec.X > 0).any(1).all()
    np.testing.assert_array_equal(dec.coords, b0.coords)
    assert not np.array_equal(dec.X, b0.X)
    # the rows are B0's draw, and each cell is decoded at its real cell's own library
    pool, rows, _, _ = harness.draw_copy(ctx)
    np.testing.assert_array_equal(lat[0]["args"][1], rows)
    real_lib = np.asarray(counts_of(pool)[rows].sum(1)).ravel()
    library = seen[0]["args"][2]
    np.testing.assert_array_equal(library, real_lib)
    sampled = np.asarray(seen[0]["out"][1].sum(1)).ravel()
    assert np.mean(sampled / library) == pytest.approx(1.0, abs=0.02)  # up to sampling
    # another seed draws other cells and other counts; the same seed the same file
    assert not np.array_equal(CANDIDATES["b0_lat_dec"].build(heart_ctx(heart_inputs, seed=ctx.seed + 1)).X, dec.X)
    np.testing.assert_array_equal(CANDIDATES["b0_lat_dec"].build(ctx).X, dec.X)


def test_t1_decoded_arm_is_b0s_rows_at_the_cells_total_counts(standin, t1_inputs, monkeypatch):
    ref, start = t1_inputs
    ctx = t1_ctx(ref, start)
    seen = _spy(monkeypatch, "decoded_counts")
    lat = _spy(monkeypatch, "moved_latent")
    dec = CANDIDATES["b0_lat_dec"].build(ctx)
    b0 = CANDIDATES["b0"].build(ctx)
    assert dec.X.shape == b0.X.shape and np.isfinite(dec.X).all() and (dec.X >= 0).all()
    pop = ts.b0_population(harness._shallow(ref), harness._shallow(start),
                           harness.t1_chamber_shares(harness._shallow(ref), "E9.5_RNA"), ctx.board.max_cells,
                           ctx.seed, HELD, ref_source="E9.5_RNA", start_source="E8.5_RNA", damping=0.0)
    rows = ref.obs_names.get_indexer(pop.obs["source_obs"].astype(str))
    np.testing.assert_array_equal(lat[0]["args"][1], rows)
    np.testing.assert_array_equal(seen[0]["args"][2], np.asarray(counts_of(ref)[rows].sum(1)).ravel())
    assert not seen[0]["args"][3].any()  # every gene is modelled: nothing outside
    # log1p CP10k of the sampled counts, as the released files are
    np.testing.assert_allclose(np.expm1(dec.X.astype(np.float64)).sum(1), 1e4, rtol=1e-4)


def test_decoded_counts_keeps_the_library_and_counts_outside_the_modelled_genes():
    rng = np.random.default_rng(3)
    p = np.array([0.5, 0.3, 0.2])
    model = t2_latent.NBVAE(latent=np.zeros((1, 1)), decode=lambda z, lib: (lib[:, None] * p, np.full(3, 1e6)))
    lib = np.r_[np.full(1500, 400.0), np.full(1500, 40.0)]
    other = np.r_[np.full(1500, 100.0), np.zeros(1500)]
    X, C, retried = latent_extrap.decoded_counts(model, np.zeros((3000, 1)), lib, other, rng, chunk=700)
    tot = np.asarray(C.sum(1)).ravel()
    assert tot[:1500].mean() == pytest.approx(400, rel=0.01) and tot[1500:].mean() == pytest.approx(40, rel=0.02)
    implied = np.expm1(X.astype(np.float64)).sum(1) / 1e4
    assert implied[:1500].mean() == pytest.approx(400 / 500, abs=0.01)
    np.testing.assert_allclose(implied[1500:], 1.0, rtol=1e-4)
    with pytest.raises(ValueError):
        latent_extrap.decoded_counts(model, np.zeros((2, 1)), np.array([1.0, 0.0]), np.zeros(2), rng)


# -------------------------------------------------------------------------------------------- (c) unmoved categories
def test_heart_moves_exactly_b0s_lineages(standin, heart_inputs, monkeypatch):
    ctx = heart_ctx(heart_inputs)
    lat = _spy(monkeypatch, "moved_latent")
    ratio = CANDIDATES["b0_lat_ratio"].build(ctx)
    b0 = CANDIDATES["b0"].build(ctx)
    # B0's own report says which lineages its momentum moves
    pool, rows, _, _ = harness.draw_copy(ctx)
    tp, ts_ = latent_extrap._typed(pool, "end"), latent_extrap._typed(ctx.inputs["start"], "start")
    _, rep = t1_scope.momentum_shift(ts_, tp, tp[rows], 0.5, (), key="b0_type", seed=[0, 3])
    shifted = sorted(r["type"] for r in rep.rows if r["status"] == "shifted")
    assert shifted == ["CM", "Endothelial"]
    moved = lat[0]["args"][3]
    assert sorted(moved) == shifted
    z_own, z_new, _ = lat[0]["out"]
    merged = np.asarray(pool.obs["merged"]).astype(str)[rows]
    for kept in ("Untracked", "Tiny"):
        m = merged == kept
        assert m.any()
        np.testing.assert_array_equal(z_new[m], z_own[m])
        np.testing.assert_array_equal(ratio.X[m], b0.X[m])  # B0 leaves them too: the same counts, the same rows
    cm = merged == "CM"
    assert np.abs(z_new[cm] - z_own[cm]).max() > 0 and not np.array_equal(ratio.X[cm], b0.X[cm])
    # the displacement is the lineage's mean latent change, damped
    fold = latent_extrap._FOLD["fold"]["fold"]
    s_rows, r_rows = moved["CM"]
    d = fold.z_ref[r_rows].mean(0) - fold.z_start[s_rows].mean(0)
    np.testing.assert_allclose(z_new[cm] - z_own[cm], np.tile(harness.B0_DAMPING * d, (cm.sum(), 1)), atol=1e-12)


def test_t1_moves_exactly_b0s_categories(standin, t1_inputs, monkeypatch):
    ref, start = t1_inputs
    ctx = t1_ctx(ref, start)
    lat = _spy(monkeypatch, "moved_latent")
    CANDIDATES["b0_lat_ratio"].build(ctx)
    r0 = harness._shallow(ref)
    shares = harness.t1_chamber_shares(r0, "E9.5_RNA")
    _, rep = t1_scope.build_submission(r0, ctx.board, shares, ctx.board.max_cells, ctx.seed,
                                       start=harness._shallow(start), damping=0.5, held_flat=HELD,
                                       ref_source="E9.5_RNA", start_source="E8.5_RNA")
    rows = rep["shift"]["rows"]
    shifted = sorted(r["type"] for r in rows if r["status"] == "shifted")
    assert shifted and len(shifted) < len(rows)  # the fixture has categories B0 does not move
    assert sorted(lat[0]["args"][3]) == shifted
    z_own, z_new, _ = lat[0]["out"]
    cats = np.asarray(lat[0]["args"][2]).astype(str)
    still = ~np.isin(cats, shifted)
    assert still.any()
    np.testing.assert_array_equal(z_new[still], z_own[still])


# ------------------------------------------------------------------------------------------------ (d) the ratio arm
def test_ratio_counts_hold_held_genes_keep_zeros_and_clip():
    rng = np.random.default_rng(5)
    G = 6
    C = sparse.csr_matrix(rng.poisson(20.0, size=(300, G)) * (rng.random((300, G)) > 0.2))
    w = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0])

    def decode(z, lib):  # mean proportional to exp(z): a step of 3 in z[0] is a factor e^3 > 3, so it is clipped
        e = np.exp(np.asarray(z)) * w
        return np.asarray(lib)[:, None] * e / e.sum(1, keepdims=True), np.full(G, 1e3)

    model = t2_latent.NBVAE(latent=np.zeros((1, G)), decode=decode)
    z_own = np.zeros((300, G))
    z_new = z_own.copy()
    z_new[:200, 0] += 3.0
    z_new[:200, 1] -= 0.2
    held = np.zeros(G, dtype=bool)
    held[5] = True
    new, rec = latent_extrap.ratio_counts(model, C, z_own, z_new, held, [1, 21], chunk=128)
    Cd, Nd = C.toarray(), new.toarray()
    np.testing.assert_array_equal(Nd[:, 5], Cd[:, 5])  # held genes exactly as they were
    np.testing.assert_array_equal(Nd[200:], Cd[200:])  # unmoved cells exactly as they were
    assert not Nd[Cd == 0].any()  # zeros stay zero
    assert rec["n_cells_moved"] == 200 and rec["n_at_clip"] > 0
    assert rec["p99_log_factor"] <= latent_extrap.CLIP + 1e-12
    # each gene moves by its clipped log ratio (the depth constant cancels between two genes of a cell)
    mu_new, _ = decode(z_new[:1], np.ones(1))
    mu_own, _ = decode(z_own[:1], np.ones(1))
    lf = np.clip(np.log(mu_new[0] / mu_own[0]), -latent_extrap.CLIP, latent_extrap.CLIP)
    assert lf[0] == pytest.approx(latent_extrap.CLIP)  # the step of gene 0 is past the clip
    f = np.exp(lf)
    f[5] = 1.0
    free = ~held
    s = (Cd[:200, free]).sum(1) / (Cd[:200, free] * f[free]).sum(1)  # apply_factors' keep-depth constant per cell
    expect = (Cd[:200] * f * np.where(free, 1.0, 0.0) * s[:, None]).sum(0)
    np.testing.assert_allclose(Nd[:200, :5].sum(0), expect[:5], rtol=0.06)
    # expected depth kept over the non-held genes
    assert Nd[:200, :5].sum() / Cd[:200, :5].sum() == pytest.approx(1.0, abs=0.02)
    # the same seed the same counts
    again, _ = latent_extrap.ratio_counts(model, C, z_own, z_new, held, [1, 21], chunk=128)
    assert (again != new).nnz == 0


def test_t1_ratio_arm_leaves_held_flat_genes_unchanged(standin, t1_inputs, monkeypatch):
    ref, start = t1_inputs
    ctx = t1_ctx(ref, start)
    seen = _spy(monkeypatch, "ratio_counts")
    out = CANDIDATES["b0_lat_ratio"].build(ctx)
    b0 = CANDIDATES["b0"].build(ctx)
    a = seen[0]
    C_in, held = a["args"][1], a["args"][4]
    genes = list(ref.var_names)
    assert [genes[i] for i in np.flatnonzero(held)] == HELD
    new, rec = a["out"]
    np.testing.assert_array_equal(new[:, held].toarray(), C_in[:, held].toarray())
    assert rec["n_cells_moved"] > 0 and (new != C_in).nnz > 0
    assert out.X.shape == b0.X.shape and not np.array_equal(out.X, b0.X)


# ------------------------------------------------------------------------------------------------ (e) once per fold
def test_the_model_is_fitted_once_per_fold(standin, heart_inputs, t1_inputs):
    for s in range(2):
        for name in ARMS:
            CANDIDATES[name].build(heart_ctx(heart_inputs, seed=harness.BUILD_SEED_OFFSET + s))
    assert len(standin) == 1
    n_ref, n_start = heart_inputs["ref"].n_obs, heart_inputs["start"].n_obs
    assert standin[0]["shape"] == (n_start + n_ref, len(HEART_BOARD.genes))
    assert standin[0]["seed"] == t2_latent.TRAIN_SEED and standin[0]["device"] == "cpu"
    fold = latent_extrap._FOLD["fold"]["fold"]
    assert fold.n_start == n_start and fold.record["final_train_loss"] == 1.0
    ref, start = t1_inputs
    for s in range(2):
        CANDIDATES["b0_lat_ratio"].build(t1_ctx(ref, start, seed=harness.BUILD_SEED_OFFSET + s))
    assert len(standin) == 2 and standin[1]["shape"] == (start.n_obs + ref.n_obs, len(T1_GENES))
    # new input objects are a new fold
    fresh = {"ref": heart_inputs["ref"].copy(), "start": heart_inputs["start"]}
    CANDIDATES["b0_lat_dec"].build(heart_ctx(fresh))
    assert len(standin) == 3


def test_interpolation_and_a_missing_start_are_refused(standin, heart_inputs):
    from test_h7_gmom import BOARD

    t = harness.TaskView.of(_heart_task(kind="interp", reference="T2/E8.75.h5ad", target="T2/E9.5.h5ad",
                                        inputs={"before": "T2/E8.75.h5ad", "after": "T2/E10.5.h5ad"}))
    ctx = harness.BuildContext(task=t, board=BOARD, inputs={"before": heart_inputs["start"],
                                                            "after": heart_inputs["ref"]}, seed=1)
    with pytest.raises(ValueError, match="start"):
        harness.b0_lat_population(ctx, "ratio")


# ------------------------------------------------------------------------------------------------ the real training
def test_both_arms_build_with_the_real_count_model_on_the_cpu(heart_inputs, monkeypatch):
    """Skipped in the project's main environment (no torch, no scvi-tools); runs with the gpu dependency group, on the
    CPU and with a few epochs (a test-only override of train_nbvae's max_epochs; the candidates never pass one)."""
    pytest.importorskip("torch")
    pytest.importorskip("scvi")
    real = t2_latent.train_nbvae

    def short(counts, seed=t2_latent.TRAIN_SEED, **kw):
        return real(counts, seed, **{**kw, "max_epochs": 3, "threads": 2})

    monkeypatch.setattr(t2_latent, "train_nbvae", short)
    monkeypatch.setattr(latent_extrap, "fit_device", lambda: "cpu")
    monkeypatch.setattr(latent_extrap, "_FOLD", {})
    ctx = heart_ctx(heart_inputs)
    b0 = CANDIDATES["b0"].build(ctx)
    for arm in ("decoded", "ratio"):
        pool, X, coords, rec = harness.b0_lat_population(ctx, arm)
        assert X.shape == (HEART_BOARD.max_cells, len(HEART_BOARD.genes)) and np.isfinite(X).all()
        np.testing.assert_array_equal(coords, b0.coords)
        assert rec["latent"]["categories_moved"] == ["CM", "Endothelial"]
        assert rec["fold"]["n_latent"] == t2_latent.N_LATENT and np.isfinite(rec["fold"]["final_train_loss"])
        assert rec["fold"]["device"] == "cpu"
    for name in ARMS:
        st = CANDIDATES[name].build(ctx)
        assert st.X.shape == b0.X.shape
    m = latent_extrap._FOLD["fold"]["fold"]
    assert m.model.latent.shape == (heart_inputs["start"].n_obs + heart_inputs["ref"].n_obs, t2_latent.N_LATENT)


def test_log1p_cp10k_is_what_b0_writes_for_unmoved_counts():
    """The ratio arm normalises as B0's momentum_shift does (vec.prep.log1p_cp10k of the edited counts)."""
    C = sparse.csr_matrix(np.array([[1, 0, 3], [0, 2, 2]]))
    np.testing.assert_array_equal(log1p_cp10k(C).toarray(),
                                  np.log1p(np.diag(1e4 / np.array([4.0, 4.0])) @ C.toarray()).astype(np.float32))


# ================================================================================================ b0_lat_flow (heart)
FLOW_INPUTS = {"ref": "T2/E9.5.h5ad", "start": "T2/E8.75.h5ad", "earlier": "T2/E8.5.h5ad"}


@pytest.fixture(scope="module")
def flow_inputs():
    """Three heart-like stages (days 8.5, 8.75, 9.5): CM and Endothelial lineages, Untracked cells, and a lineage too
    small at the start stage for B0 (Tiny)."""
    return {"ref": heart_stage(1600, seed=2, shift=0.0, tiny=60),
            "start": heart_stage(1100, seed=3, shift=0.8, tiny=10),
            "earlier": heart_stage(900, seed=5, shift=1.6, tiny=10)}


def flow_ctx(inputs, seed: int = harness.BUILD_SEED_OFFSET + 1, view=None) -> harness.BuildContext:
    view = view or harness.TaskView.of(_heart_task(inputs=FLOW_INPUTS))
    return harness.BuildContext(task=view, board=HEART_BOARD, inputs=inputs, seed=seed)


def _constant_field(monkeypatch):
    """fold_flow replaced by the constant per-lineage field (mean z_ref(l) - mean z_start(l)) / D_{K-1}."""

    def fold_flow(fold, stages, pool, start, lineages):
        mp = np.asarray(pool.obs["merged"]).astype(str)
        ms = np.asarray(start.obs["merged"]).astype(str)
        gap = stages[-1][1] - stages[-2][1]
        c = np.array([(fold.z_ref[mp == lin].mean(0) - fold.z_start[ms == lin].mean(0)) / gap for lin in lineages])
        return latent_extrap.FlowField(lineages=list(lineages), velocity=lambda z, li: c[li].copy(),
                                       record={"constant": True})

    monkeypatch.setattr(latent_extrap, "fold_flow", fold_flow)


def test_flow_constant_field_integrates_to_the_ratio_arms_z_new(standin, flow_inputs, monkeypatch):
    """Spec section 3 identity, part (a): with the constant field, RK4 over B0's time budget lands on b0_lat_ratio's
    z_new to 1e-12."""
    _constant_field(monkeypatch)
    seen = _spy(monkeypatch, "ratio_counts")
    ctx = flow_ctx(flow_inputs)
    harness.b0_lat_population(ctx, "ratio")
    harness.b0_lat_population(ctx, "flow")
    (r_own, r_new), (f_own, f_new) = ((s["args"][2], s["args"][3]) for s in seen)
    np.testing.assert_array_equal(f_own, r_own)
    np.testing.assert_allclose(f_new, r_new, rtol=0, atol=1e-12)
    assert np.abs(r_new - r_own).max() > 0


def test_flow_with_the_ratio_arms_z_new_is_the_ratio_arm_bit_for_bit(standin, flow_inputs, monkeypatch):
    """Spec section 3 identity, part (b), the vehicle: handed the very z_new of b0_lat_ratio (injected, not
    integrated), b0_lat_flow's file is b0_lat_ratio's to the last bit (one model, the same cells, ratio_counts call and
    stream)."""

    def inject(fold, stages, pool, start, types, moved, z_own, damping):
        z_new = z_own.copy()
        for lin, (s_rows, r_rows) in moved.items():
            d = fold.z_ref[r_rows].mean(0) - fold.z_start[s_rows].mean(0)
            m = types == lin
            z_new[m] = z_own[m] + float(damping) * d
        return z_new, {"injected": True}

    monkeypatch.setattr(latent_extrap, "flow_latent", inject)
    ctx = flow_ctx(flow_inputs)
    ratio, flow = CANDIDATES["b0_lat_ratio"].build(ctx), CANDIDATES["b0_lat_flow"].build(ctx)
    np.testing.assert_array_equal(flow.X, ratio.X)
    np.testing.assert_array_equal(flow.coords, ratio.coords)
    assert len(standin) == 1  # one count model for both arms
    shas = {r["fold"]["latent_sha256"] for r in latent_extrap.RECORDS}
    assert len(shas) == 1 and {r["candidate"] for r in latent_extrap.RECORDS} == {"b0_lat_ratio", "b0_lat_flow"}


def test_flow_weights_are_inverse_variance_and_average_one():
    np.testing.assert_allclose(latent_extrap.block_weights([0.2, 0.25, 0.75]), [0.0602, 0.0940, 0.8459], atol=5e-5)
    np.testing.assert_allclose(latent_extrap.block_weights([0.2, 0.25]), [0.3902, 0.6098], atol=5e-5)
    rng = np.random.default_rng(0)
    # lineage A at all three stages, B too small at the first: A has two blocks, B one
    labs = [np.array(["A"] * 60 + ["B"] * 5), np.array(["A"] * 50 + ["B"] * 40), np.array(["A"] * 70 + ["B"] * 30)]
    zs = [rng.normal(size=(len(x), 3)) + k for k, x in enumerate(labs)]
    pairs, rec = latent_extrap.flow_pairs(zs, labs, [8.5, 8.75, 9.5], ["A", "B"])
    assert rec["n_blocks"] == 3 and rec["n_pairs"] == 3 * latent_extrap.FLOW_PAIRS
    assert pairs["omega"].mean() == pytest.approx(1.0)
    w = {(b["lineage"], b["transition"]): b["weight"] for b in rec["blocks"]}
    assert w[("A", 1)] == pytest.approx(0.0625 / (0.0625 + 0.5625)) and w[("B", 2)] == 1.0
    assert all(b["converged"] for b in rec["blocks"])
    with pytest.raises(ValueError, match="last transition"):  # a moved lineage must have the last transition
        latent_extrap.flow_pairs(zs, [labs[0], labs[1], np.array(["A"] * 100)], [8.5, 8.75, 9.5], ["A", "B"])


def test_flow_stages_are_in_day_order_with_start_then_the_reference_last():
    tasks = {t.key: t for ts_ in harness.load_tasks().values() for t in ts_}
    st = latent_extrap.flow_stages(harness.TaskView.of(tasks["T2h/<=E9.5->E10.5"]))
    assert [r for r, _ in st] == ["earliest", "earlier", "start", "ref"]
    assert [d for _, d in st] == [8.3, 8.5, 8.75, 9.5]
    assert [r for r, _ in latent_extrap.flow_stages(harness.TaskView.of(tasks["T2h/<=E8.75->E10.5"]))] == [
        "earlier", "start", "ref"]
    with pytest.raises(ValueError, match="three"):
        latent_extrap.flow_stages(harness.TaskView.of(_heart_task()))
    v = harness.TaskView.of(_heart_task(inputs=FLOW_INPUTS))
    with pytest.raises(ValueError, match="start"):  # an input after the start stage that is not the reference
        latent_extrap.flow_stages(replace(v, role_days={**v.role_days, "earlier": 9.0}))
    with pytest.raises(ValueError, match="reference"):
        latent_extrap.flow_stages(replace(v, role_days={**v.role_days, "ref": 8.6}))


def _planted(rng, v, days, n=600, d=4):
    """Two lineages drifting along latent axis 0 at +v and -v per day through the given days."""
    zs, labs = [], []
    for day in days:
        z = rng.normal(size=(2 * n, d)) * 0.5
        z[:n, 0] += v * (day - days[0])
        z[n:, 0] -= v * (day - days[0])
        z[n:, 1] += 3.0
        zs.append(z)
        labs.append(np.array(["A"] * n + ["B"] * n))
    return zs, labs


def test_flow_recovers_a_planted_linear_flow(monkeypatch):
    pytest.importorskip("torch")
    monkeypatch.setattr(latent_extrap, "FLOW_STEPS", 1500)
    days, v = [8.5, 8.75, 9.5], 2.0
    zs, labs = _planted(np.random.default_rng(1), v, days)
    pairs, _ = latent_extrap.flow_pairs(zs, labs, days, ["A", "B"])
    ff = latent_extrap.fit_flow(pairs, ["A", "B"])
    tau = 0.5 * (days[-1] - days[-2])
    for li, sign in ((0, 1.0), (1, -1.0)):
        ref = zs[-1][labs[-1] == ["A", "B"][li]]
        vel = ff.velocity(ref, np.full(len(ref), li))[:, 0].mean()
        assert np.sign(vel) == sign and abs(vel) == pytest.approx(v, rel=0.25)
        z_new = latent_extrap.integrate(ff.velocity, ref, np.full(len(ref), li), tau)
        assert (z_new - ref)[:, 0].mean() == pytest.approx(sign * v * tau, rel=0.25)


def test_flow_fit_is_deterministic(monkeypatch):
    pytest.importorskip("torch")
    monkeypatch.setattr(latent_extrap, "FLOW_STEPS", 200)
    days = [8.5, 8.75, 9.5]
    zs, labs = _planted(np.random.default_rng(2), 1.0, days, n=2500)  # above FLOW_SIDE_MAX: subsampled sides
    p1, r1 = latent_extrap.flow_pairs(zs, labs, days, ["A", "B"])
    p2, _ = latent_extrap.flow_pairs(zs, labs, days, ["A", "B"])
    for k in p1:
        np.testing.assert_array_equal(p1[k], p2[k])
    assert all(b["n_earlier"] == latent_extrap.FLOW_SIDE_MAX for b in r1["blocks"])
    a, b = latent_extrap.fit_flow(p1, ["A", "B"]), latent_extrap.fit_flow(p2, ["A", "B"])
    assert a.record["params_sha256"] == b.record["params_sha256"] and a.record["loss_last"] == b.record["loss_last"]


def test_flow_never_reads_the_target(standin, flow_inputs, monkeypatch):
    pytest.importorskip("torch")
    monkeypatch.setattr(latent_extrap, "FLOW_STEPS", 100)
    view = harness.TaskView.of(_heart_task(inputs=FLOW_INPUTS))
    blind = replace(view, target_day=float("nan"))
    assert set(flow_inputs) == set(FLOW_INPUTS)  # no target among the inputs
    with harness.guarded_reads() as g:
        a = CANDIDATES["b0_lat_flow"].build(flow_ctx(flow_inputs, view=blind))
    assert not g["refused"] and np.isfinite(a.X).all()
    rec = latent_extrap.RECORDS[-1]
    assert rec["flow"]["stages"] == ["earlier", "start", "ref"] and rec["flow"]["tau_days"] == pytest.approx(0.375)
    assert rec["flow"]["n_blocks"] >= 2 and np.isfinite(rec["flow"]["loss_last"])
    b = CANDIDATES["b0_lat_flow"].build(flow_ctx(flow_inputs, view=view))  # the target's day changes nothing
    np.testing.assert_array_equal(a.X, b.X)
    assert len(standin) == 1  # the fold model is fitted once


def test_flow_leaves_unmoved_lineages_and_untracked_cells_at_their_own_latent(standin, flow_inputs, monkeypatch):
    _constant_field(monkeypatch)
    seen = _spy(monkeypatch, "ratio_counts")
    ctx = flow_ctx(flow_inputs)
    flow, ratio = CANDIDATES["b0_lat_flow"].build(ctx), CANDIDATES["b0_lat_ratio"].build(ctx)
    pool, rows, _, _ = harness.draw_copy(ctx)
    merged = np.asarray(pool.obs["merged"]).astype(str)[rows]
    z_own, z_new = seen[0]["args"][2], seen[0]["args"][3]
    for kept in ("Untracked", "Tiny"):
        m = merged == kept
        assert m.any()
        np.testing.assert_array_equal(z_new[m], z_own[m])
        np.testing.assert_array_equal(flow.X[m], ratio.X[m])  # factor 1: the cells' own counts, as in the ratio arm
    rec = [r for r in latent_extrap.RECORDS if r["candidate"] == "b0_lat_flow"][0]
    assert rec["latent"]["categories_moved"] == ["CM", "Endothelial"]
    per = rec["latent"]["per_lineage"]
    assert set(per["lineages"]) == {"CM", "Endothelial"}
    for lin in per["lineages"].values():  # the constant field is B0's displacement: cosine 1, norm ratio 1
        assert lin["cosine_to_b0"] == pytest.approx(1.0) and lin["norm_ratio_to_b0"] == pytest.approx(1.0)
        assert lin["within_share"] == pytest.approx(0.0, abs=1e-12)


def test_encode_reproduces_the_training_latents_of_the_standin():
    C = sparse.csr_matrix(np.random.default_rng(3).poisson(4.0, size=(50, 12)) + 1)
    m = _standin([])(C)
    np.testing.assert_allclose(m.encode(C), m.latent, atol=1e-12)


def test_encode_reproduces_the_training_latents_of_a_real_fit():
    pytest.importorskip("torch")
    pytest.importorskip("scvi")
    rng = np.random.default_rng(30)
    C = sparse.csr_matrix(rng.poisson(rng.gamma(2.0, 2.0, size=(1, 40)), size=(300, 40)) + 1)
    m = t2_latent.train_nbvae(C, 123, max_epochs=3, threads=2, device="cpu")
    np.testing.assert_allclose(m.encode(C), m.latent, atol=1e-5)


def test_flow_applies_to_the_heart_extrapolation_tasks_only():
    tasks = {t.key: t for ts_ in harness.load_tasks().values() for t in ts_}
    c = CANDIDATES["b0_lat_flow"]
    assert {k for k, t in tasks.items() if c.applies(t)} == {"T2h/<=E9.5->E10.5", "T2h/<=E8.75->E10.5"}
    test_like = _heart_task(key="heart_test_like", inputs={"ref": "T2/E9.5.h5ad", "start": "T2/E8.75.h5ad",
                                                           "earlier": "T2/E8.5.h5ad",
                                                           "earliest": "T2/E8.25_late.h5ad"})
    assert c.applies(test_like)
    for k in ("T1/E9.5->E10.5", "T1/E9.5->E10.5:chamber", "T1/E8.5->E9.5", "T2e/E7.25+E8.0->E7.5",
              "T2e/E6.75+E7.5->E7.25"):
        assert not c.applies(tasks[k])
    assert not c.applies(replace(test_like, setting="embryo"))
    assert c.derived == () and c.params["flow_seed"] == latent_extrap.FLOW_SEED
    assert c.params["scope"] == CANDIDATES["b0_lat_ratio"].params["scope"]


def test_flow_is_refused_on_t1(standin, t1_inputs):
    ref, start = t1_inputs
    with pytest.raises(ValueError, match="heart"):
        harness.b0_lat_population(t1_ctx(ref, start), "flow")


def test_records_reach_a_file_per_task(standin, heart_inputs, tmp_path):
    ctx = heart_ctx(heart_inputs)
    for name in ARMS:
        CANDIDATES[name].build(ctx)
    latent_extrap.RECORDS.append({"task": "another", "candidate": "x"})  # left by a failed task: dropped
    n = latent_extrap.write_records(tmp_path / "latent_records" / "t.json", ctx.task.key)
    recs = json.loads((tmp_path / "latent_records" / "t.json").read_text())
    assert n == 2 and [r["candidate"] for r in recs] == ["b0_lat_dec", "b0_lat_ratio"]
    assert len({r["fold"]["latent_sha256"] for r in recs}) == 1 and not latent_extrap.RECORDS
    assert latent_extrap.write_records(tmp_path / "none.json", "nothing") == 0
    assert not (tmp_path / "none.json").exists()
    src = (REPO / "experiments" / "harness" / "run.py").read_text()
    assert 'latent_extrap.write_records(out / "latent_records" / f"{safe(t.key)}.json", t.key)' in src
