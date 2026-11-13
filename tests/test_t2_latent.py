"""vec.t2_latent (embryo try E14, candidate fgw_latent), its target-free fidelity gate N18
(experiments/harness/n18_latent_fidelity.py) and its smoke abort (experiments/harness/e14_smoke.py), on synthetic data
only, with a stand-in count model in place of the scvi-tools training (torch and scvi-tools live in the `gpu` dependency
group, not in the project's main environment). The one test of the real training, `test_scvi_training_...`, is skipped
there and runs where the gpu group is installed."""

from __future__ import annotations

import importlib.util
import json
import re
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from conftest import GENES
from scipy import sparse

from vec import harness, t2_build, t2_latent
from vec.boards import Board
from vec.harness import CANDIDATES, load_tasks
from vec.io import COORD_KEY
from vec.provenance import stamp_path

REPO = Path(__file__).resolve().parents[1]
G = len(GENES)


# ------------------------------------------------------------------------------------------------ the latent formula
def test_interpolation_preserves_the_latent_variance_where_an_average_loses_it():
    rng = np.random.default_rng(1)
    n, rho = 200_000, 0.6
    zb = rng.standard_normal((n, 3))
    za = rho * zb + np.sqrt(1 - rho**2) * rng.standard_normal((n, 3))
    cov = np.eye(3)
    for t in (0.25, 0.5, 0.75):
        z = t2_latent.interpolate_latent(zb, za, t, rho, cov, np.random.default_rng(2))
        assert z.var(0) == pytest.approx(np.ones(3), abs=0.02)  # the noise term restores V
    mid = (zb + za) / 2  # a plain average, what the noise term corrects: variance (1 + rho) / 2 of V
    assert mid.var(0) == pytest.approx(np.full(3, (1 + rho) / 2), abs=0.02)
    # the endpoints and rho = 1 add no noise: the cells' own latents
    np.testing.assert_array_equal(t2_latent.interpolate_latent(zb, za, 0.0, rho, cov, rng), zb)
    np.testing.assert_array_equal(t2_latent.interpolate_latent(zb, za, 1.0, rho, cov, rng), za)
    np.testing.assert_array_equal(t2_latent.interpolate_latent(zb, za, 0.5, 1.0, cov, rng), (zb + za) / 2)
    with pytest.raises(ValueError):
        t2_latent.interpolate_latent(zb, za, 1.5, rho, cov, rng)
    with pytest.raises(ValueError):
        t2_latent.interpolate_latent(zb, za[:-1], 0.5, rho, cov, rng)


def test_coupling_correlation_reads_the_pairs_of_the_coupling():
    rng = np.random.default_rng(3)
    n = 400
    zb = rng.normal(size=(n, 5))
    za = 0.8 * zb + 0.6 * rng.normal(size=(n, 5)) + 3.0  # a shift does not matter
    identity = np.eye(n) / n
    assert t2_latent.coupling_correlation(zb, za, identity) == pytest.approx(0.8, abs=0.04)
    # a coupling of unrelated partners reads about 0, one fully crossed reads below 0
    perm = rng.permutation(n)
    assert abs(t2_latent.coupling_correlation(zb, za, identity[perm])) < 0.12
    assert t2_latent.coupling_correlation(zb, -zb, identity) == pytest.approx(-1.0)
    with pytest.raises(ValueError):
        t2_latent.coupling_correlation(zb, za, np.ones((n, n - 1)))
    with pytest.raises(ValueError):
        t2_latent.coupling_correlation(zb, za, np.zeros((n, n)))


def test_sample_nb_has_the_nb_mean_and_variance_and_rejects_bad_parameters():
    rng = np.random.default_rng(4)
    mu = np.array([[2.0, 40.0]] * 300_000)
    theta = np.array([3.0, 5.0])
    c = t2_latent.sample_nb(mu, theta, rng)
    assert c.shape == mu.shape and c.dtype.kind == "i"
    assert c.mean(0) == pytest.approx([2.0, 40.0], rel=0.01)
    assert c.var(0) == pytest.approx(mu[0] + mu[0] ** 2 / theta, rel=0.03)
    a, b = (t2_latent.sample_nb(mu[:100], theta, np.random.default_rng(5)) for _ in range(2))
    np.testing.assert_array_equal(a, b)
    assert not t2_latent.sample_nb(np.zeros((3, 2)), theta, rng).any()  # a zero mean samples zero
    for bad in (dict(mu=-mu[:1], theta=theta), dict(mu=mu[:1] * np.nan, theta=theta), dict(mu=mu[:1], theta=-theta),
                dict(mu=mu[:1], theta=theta * 0.0)):
        with pytest.raises(ValueError):
            t2_latent.sample_nb(rng=rng, **bad)


def test_batch_size_avoids_a_last_batch_of_one_or_two_cells():
    assert t2_latent._batch_size(1000, 128) == 128  # 1000 % 128 = 104
    assert t2_latent._batch_size(129, 128) == 129   # 129 % 128 = 1
    assert t2_latent._batch_size(130, 128) == 130   # 130 % 128 = 2, 130 % 129 = 1, 130 % 130 = 0
    for n in range(300, 2000, 97):
        assert n % t2_latent._batch_size(n, 128) not in (1, 2)


# ------------------------------------------------------------------------------------------------ the expression step
def _stage(n, rng, profiles, props, depth, extra=15.0) -> ad.AnnData:
    """A bracket as the harness loads it: panel counts of per-type profiles at log-normal depths, plus counts outside
    the panel (`extra`, so N exceeds the panel's sum), X = log(1 + 1e4 c / N), and coordinates."""
    types = rng.choice(len(props), size=n, p=props)
    lib = depth * rng.lognormal(0.0, 0.3, size=n)
    C = rng.poisson(lib[:, None] * profiles[types])
    C[C.sum(1) == 0, 0] = 1
    other = rng.poisson(extra, size=n)
    N = C.sum(1) + other
    a = ad.AnnData(X=sparse.csr_matrix(np.log1p(1e4 * C / N[:, None]).astype(np.float32)),
                   obs=pd.DataFrame({"celltype": [f"type{t}" for t in types]}, index=[f"c{i}" for i in range(n)]),
                   var=pd.DataFrame(index=GENES))
    a.layers["counts"] = sparse.csr_matrix(C.astype(np.int32))
    a.obsm[COORD_KEY] = (types[:, None] * 5.0 + rng.normal(size=(n, 3))).astype(np.float32)
    a.obs["true_type"] = types
    a.uns["other"] = other
    return a


def _profiles(rng, k=3):
    """k cell-type profiles on the panel, each dominated by its own block of genes."""
    p = rng.gamma(1.0, 1.0, size=(k, G)) + 0.05
    for i in range(k):
        p[i, i * 10:(i + 1) * 10] *= 12.0
    return p / p.sum(1, keepdims=True)


def _oracle_trainer(profiles, flat=False):
    """A stand-in for the training: the latent means are each cell's type (read from its counts), the decoder the type's
    profile at the given library (flat=True: a uniform profile, a bad generator)."""
    calls = []

    def train(counts, seed=t2_latent.TRAIN_SEED, *, threads=None, device=None, **kw):
        C = sparse.csr_matrix(counts).toarray().astype(np.float64)
        calls.append(C.shape)
        loglik = C @ np.log(profiles.T + 1e-12)
        types = loglik.argmax(1)
        latent = np.eye(len(profiles))[types] + 1e-3 * np.random.default_rng(0).normal(size=(C.shape[0], len(profiles)))

        def decode(z, lib):
            z, lib = np.asarray(z), np.asarray(lib)
            prof = np.full(G, 1.0 / G)[None, :] if flat else profiles[np.argmax(z, 1)]
            return lib[:, None] * prof, np.full(G, 1e5)

        return t2_latent.NBVAE(latent=latent, decode=decode, record={"stand_in": True})

    return train, calls


def test_decode_and_sample_uses_the_given_libraries_and_counts_outside_the_panel():
    rng = np.random.default_rng(6)
    p = np.array([0.5, 0.3, 0.2])
    z = np.zeros((4000, 2))
    library = np.full(4000, 600.0)
    library[2000:] = 60.0
    other = np.full(4000, 150.0)
    other[1000:2000] = 0.0
    model = t2_latent.NBVAE(latent=z, decode=lambda zz, lib: (lib[:, None] * p, np.full(3, 1e6)))
    X, retried = t2_latent.decode_and_sample(model, z, library, other, rng)
    assert X.dtype == np.float32 and X.shape == (4000, 3) and retried == 0
    implied = np.expm1(X.astype(np.float64)).sum(1) / 1e4  # the sampled panel's share of N
    assert implied[:1000].mean() == pytest.approx(600 / 750, abs=0.01)
    assert implied[1000:2000].mean() == pytest.approx(1.0, abs=0.01)  # no counts outside the panel: the panel is N
    assert implied[2000:].mean() == pytest.approx(60 / 210, abs=0.02)  # the low-library cells: their own library
    assert np.expm1(X).sum(0)[0] / np.expm1(X).sum() == pytest.approx(0.5, abs=0.02)
    with pytest.raises(ValueError):  # a cell without a panel library cannot be decoded
        t2_latent.decode_and_sample(model, z[:2], np.array([3.0, 0.0]), other[:2], rng)


def test_decode_and_sample_resamples_a_cell_left_without_panel_counts():
    rng = np.random.default_rng(7)
    z = np.zeros((3000, 1))
    model = t2_latent.NBVAE(latent=z, decode=lambda zz, lib: (lib[:, None] * np.array([1.0]), np.array([1e6])))
    lib = np.full(3000, 0.7)  # Poisson(0.7): about half the cells draw zero and are sampled again
    X, retried = t2_latent.decode_and_sample(model, z, lib, np.ones(3000), rng)
    assert retried > 500 and (X.sum(1) > 0).all()  # no cell is left empty


def _interp(n_b, n_a, stage, rows, partners, n_anchor=40):
    T = np.eye(n_anchor) / n_anchor
    return SimpleNamespace(stage=np.asarray(stage), rows=np.asarray(rows), partners=np.asarray(partners),
                           coupling={"anchors_before": np.arange(n_anchor), "anchors_after": np.arange(n_anchor),
                                     "moves": T})


def test_latent_midpoint_counts_is_the_formula_at_the_drawn_cells_library():
    rng = np.random.default_rng(8)
    profiles = _profiles(rng)
    before = _stage(300, rng, profiles, [0.5, 0.3, 0.2], 80.0)
    after = _stage(400, rng, profiles, [0.2, 0.3, 0.5], 200.0)
    genes = GENES
    Zb, Za = rng.normal(size=(300, 4)), rng.normal(size=(400, 4))
    Za[:40] = 0.9 * Zb[:40] + 0.3 * rng.normal(size=(40, 4)) + 1.0  # anchor pairs correlated
    seen = {}

    def decode(z, lib):
        seen["z"], seen["lib"] = np.array(z), np.array(lib)
        return lib[:, None] * profiles[0], np.full(G, 1e4)

    model = t2_latent.NBVAE(latent=np.vstack([Zb, Za]), decode=decode)
    stage = np.array([0, 0, 1, 1, 0])
    rows = np.array([5, 17, 7, 300, 2])
    partners = np.array([11, 300, 9, 22, 4])  # the other stage's rows
    interp = _interp(300, 400, stage, rows, partners)
    t = 1 / 3
    X, rec = t2_latent.latent_midpoint_counts(model, before, after, interp, genes, t, [99, 20])
    assert X.dtype == np.float32 and X.shape == (5, G) and (X > 0).any(1).all()
    assert rec["t"] == t and -1 <= rec["rho"] <= 1
    # the z the decoder got is the formula's, with the coupling's rho and the t-weighted within-bracket covariance
    rho = t2_latent.coupling_correlation(Zb[:40], Za[:40], np.eye(40) / 40)
    assert rec["rho"] == pytest.approx(rho)
    cov = (1 - t) * t2_latent.latent_covariance(Zb) + t * t2_latent.latent_covariance(Za)
    zb = Zb[[5, 17, 9, 22, 2]]
    za = Za[[11, 300, 7, 300, 4]]
    want = t2_latent.interpolate_latent(zb, za, t, rho, cov, t2_latent._streams([99, 20])[0])
    np.testing.assert_allclose(seen["z"], want)
    # each output cell's library is the drawn real cell's panel library (before cells at their rows, after cells too)
    _, _, lb, ob = t2_build._stage_counts(before, genes)
    _, _, la, oa = t2_build._stage_counts(after, genes)
    np.testing.assert_allclose(seen["lib"], [lb[5] - ob[5], lb[17] - ob[17], la[7] - oa[7], la[300] - oa[300],
                                             lb[2] - ob[2]])
    # the same call again gives the same cells; another seed other NB draws
    X2, _ = t2_latent.latent_midpoint_counts(model, before, after, interp, genes, t, [99, 20])
    np.testing.assert_array_equal(X, X2)
    X3, _ = t2_latent.latent_midpoint_counts(model, before, after, interp, genes, t, [100, 20])
    assert not np.array_equal(X, X3)


def test_latent_midpoint_counts_needs_partners_coupling_and_matching_latents():
    rng = np.random.default_rng(9)
    profiles = _profiles(rng)
    before = _stage(100, rng, profiles, [0.5, 0.3, 0.2], 80.0)
    after = _stage(120, rng, profiles, [0.2, 0.3, 0.5], 200.0)
    model = t2_latent.NBVAE(latent=np.zeros((220, 3)), decode=lambda z, lib: (lib[:, None] / G, np.full(G, 1e4)))
    i = _interp(100, 120, [0, 1], [3, 4], [5, 6], n_anchor=30)
    no_partners = SimpleNamespace(**{**i.__dict__, "partners": None})
    no_coupling = SimpleNamespace(**{**i.__dict__, "coupling": None})
    for bad in (no_partners, no_coupling):
        with pytest.raises(ValueError, match="with_partners"):
            t2_latent.latent_midpoint_counts(model, before, after, bad, GENES, 0.5, [1, 20])
    wrong = t2_latent.NBVAE(latent=np.zeros((219, 3)), decode=model.decode)
    with pytest.raises(ValueError, match="brackets"):
        t2_latent.latent_midpoint_counts(wrong, before, after, i, GENES, 0.5, [1, 20])


def test_task_model_trains_once_per_inputs_and_keys_on_them(monkeypatch):
    rng = np.random.default_rng(10)
    profiles = _profiles(rng)
    before = _stage(200, rng, profiles, [0.5, 0.3, 0.2], 80.0)
    after = _stage(210, rng, profiles, [0.2, 0.3, 0.5], 200.0)
    other = _stage(200, np.random.default_rng(11), profiles, [0.5, 0.3, 0.2], 80.0)
    train, calls = _oracle_trainer(profiles)
    monkeypatch.setattr(t2_latent, "train_nbvae", train)
    monkeypatch.setattr(t2_latent, "_MODELS", {})
    a = t2_latent.task_model(before, after, GENES)
    assert calls == [(410, G)] and a.latent.shape[0] == 410
    assert t2_latent.task_model(before, after, GENES) is a and len(calls) == 1  # the 10 seeds reuse one model
    b = t2_latent.task_model(other, after, GENES)
    assert len(calls) == 2 and b is not a  # other inputs train again
    assert t2_latent.input_key([before.layers["counts"]], GENES) != t2_latent.input_key([other.layers["counts"]], GENES)
    assert t2_latent.input_key([before.layers["counts"]], GENES) == t2_latent.input_key(
        [before.layers["counts"].tocsc()], GENES)  # the key depends on the counts, not on their storage
    monkeypatch.setattr(t2_latent, "_MODELS", {})

    def wrong(counts, *a, **kw):
        return t2_latent.NBVAE(latent=np.zeros((1, 2)), decode=lambda z, lib: None)

    monkeypatch.setattr(t2_latent, "train_nbvae", wrong)
    with pytest.raises(ValueError, match="latent"):
        t2_latent.task_model(before, after, GENES)


def test_the_module_imports_without_torch_or_scvi():
    code = ("import sys; sys.path.insert(0, 'src'); import vec.t2_latent; "
            "assert 'torch' not in sys.modules and 'scvi' not in sys.modules")
    subprocess.run([sys.executable, "-I", "-B", "-c", code], check=True, cwd=REPO, capture_output=True)


# ------------------------------------------------------------------------------------------------ the harness candidate
def test_fgw_latent_keeps_fgw_geometry_and_changes_only_the_expression(monkeypatch):
    """Try E14 (experiments/harness/TRIES.md): fgw's cells and positions, expression from the count model (the
    stand-in's NB samples here); the model is trained once for all seeds of a task."""
    rng = np.random.default_rng(12)
    profiles = _profiles(rng)
    before = _stage(2500, rng, profiles, [0.5, 0.3, 0.2], 80.0)
    after = _stage(4000, rng, profiles, [0.2, 0.3, 0.5], 200.0)
    board = Board(key="pseudo:synthetic_embryo", task="T2", genes=list(GENES), needs_coords=True, setting="embryo",
                  min_cells=100, max_cells=1500)
    INTERP = dict(task="T2", setting="embryo", kind="interp", target="T2/E7.5.h5ad", reference="T2/E7.25.h5ad",
                  twin_board="T2:embryo:val_interp", test_board="T2:embryo:test_interp")
    task = harness.PseudoTask(key="syn_embryo_latent",
                              inputs={"before": "T2/E7.25.h5ad", "after": "T2/E8.0.h5ad"}, **INTERP)
    tasks = {t.key: t for ts in load_tasks().values() for t in ts}
    c = CANDIDATES["fgw_latent"]
    assert c.params["scope"] == CANDIDATES["fgw"].params["scope"]
    assert c.applies(task) and c.applies(tasks["T2e/E6.75+E7.5->E7.25"]) and c.applies(tasks["T2e/E7.25+E8.0->E7.5"])
    assert not c.applies(tasks["T2h/<=E9.5->E10.5"]) and not c.applies(tasks["T1/E9.5->E10.5"])
    train, calls = _oracle_trainer(profiles)
    monkeypatch.setattr(t2_latent, "train_nbvae", train)
    monkeypatch.setattr(t2_latent, "_MODELS", {})
    ctx = harness.BuildContext(task=harness.TaskView.of(task), board=board, inputs={"before": before, "after": after},
                               seed=7)
    fgw = CANDIDATES["fgw"].build(ctx)
    v = c.build(ctx)
    assert calls == [(6500, G)]  # one model, on both brackets' cells
    np.testing.assert_array_equal(v.coords, fgw.coords)  # positions and draws are fgw's
    assert v.X.shape == fgw.X.shape and not np.array_equal(v.X, fgw.X)
    assert ((v.X > 0).sum(1) >= 1).all() and np.isfinite(v.X).all()
    ctx8 = harness.BuildContext(task=harness.TaskView.of(task), board=board, inputs={"before": before, "after": after},
                                seed=8)
    v8 = c.build(ctx8)
    assert calls == [(6500, G)]  # the next seed reuses the model: the seeds vary pairs and NB draws only
    assert not np.array_equal(v8.X, v.X)
    # a seed's cells are reproducible, and the guarded read path holds (no file is opened by the expression step)
    with harness.guarded_reads() as g:
        v2 = c.build(ctx)
    np.testing.assert_array_equal(v2.X, v.X)
    assert not g["refused"]


def test_the_expression_option_is_only_added_for_the_new_try():
    src = (REPO / "experiments" / "builds" / "build_test.py").read_text()
    assert '"latent"' in src and 'with_partners=a.expression in ("blend", "latent")' in src
    # the entry's default stays the released cells ("series" joined the choices with tries E24 and E25, "render_hi"
    # with request R-A's second option; neither is an entry)
    assert 'choices=("mixture", "depth", "blend", "latent", "series", "render_hi"),' in src
    assert src.count('default="mixture"') == 1
    sys.path.insert(0, str(REPO / "experiments" / "builds"))
    import build_test

    ap = build_test.build_parser()
    base = ["embryo", "--method", "fgw", "--out", "runs/x/T2embryo_s1.h5ad"]
    assert ap.parse_args(base).expression == "mixture" and ap.parse_args(base).seed >= harness.BUILD_SEED_OFFSET
    assert ap.parse_args([*base, "--expression", "latent"]).expression == "latent"
    assert ap.parse_args([*base, "--expression", "series"]).expression == "series"
    assert ap.parse_args([*base, "--expression", "render_hi"]).expression == "render_hi"
    assert 'depth_rendered_expression(before, after, res, genes, 0.5, [a.seed, 7], mode="hi")' in src


# ------------------------------------------------------------------------------------------------ the fidelity gate N18
@pytest.fixture(scope="module")
def n18():
    spec = importlib.util.spec_from_file_location("n18_latent_fidelity",
                                                  REPO / "experiments" / "harness" / "n18_latent_fidelity.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def test_the_rule_reads_the_row_at_its_boundaries(n18):
    v = n18.check_metric("variogram", 1.3, 1.0, 5.0)
    assert v["pass"] and v["ratio_to_control"] == pytest.approx(1.3)
    assert not n18.check_metric("variogram", 1.31, 1.0, 5.0)["pass"]
    assert not n18.check_metric("variogram", 1.0, 0.0, 5.0)["pass"]  # a control of zero: no ratio, no pass
    m = n18.check_metric("mmd_u", 0.11, 0.01, 1.01)  # excess 0.10 of the bracket gap 1.00
    assert m["pass"] and m["share_of_bracket_gap"] == pytest.approx(0.1)
    assert not n18.check_metric("mmd_u", 0.12, 0.01, 1.01)["pass"]
    assert not n18.check_metric("mmd_u", 0.01, 0.01, 0.0)["pass"]  # no real gap between the brackets: no pass
    with pytest.raises(ValueError):
        n18.check_metric("energy", 1.0, 1.0, 2.0)


def test_the_verdict_needs_all_four_checks(n18):
    good = {"generated": 1.0, "control": 1.0, "other_bracket": 2.0}
    values = {s: {m: dict(good) for m in n18.METRICS} for s in ("before", "after")}
    assert n18.verdict(values)[0] is True
    values["after"]["variogram"]["generated"] = 5.0
    passed, checks = n18.verdict(values)
    assert passed is False and checks["after"]["mmd_u"]["pass"] and not checks["after"]["variogram"]["pass"]


def test_the_holdout_split_and_the_repeat_rows(n18):
    train, held = n18.holdout_split(1003, [1, 2])
    assert len(held) == round(0.2 * 1003) and len(train) == 1003 - len(held)
    assert sorted(np.r_[train, held]) == list(range(1003))
    t2, h2 = n18.holdout_split(1003, [1, 2])
    np.testing.assert_array_equal(train, t2)
    np.testing.assert_array_equal(held, h2)
    assert not np.array_equal(held, n18.holdout_split(1003, [1, 3])[1])
    ha, hb, rows, others = n18.repeat_rows(np.random.default_rng(1), held, len(train), 500)
    assert len(ha) == len(hb) == len(rows) == len(others) == len(held) // 2
    assert not set(ha) & set(hb) and set(ha) | set(hb) <= set(held)
    assert rows.max() < len(train) and others.max() < 500
    with pytest.raises(ValueError):
        n18.repeat_rows(np.random.default_rng(1), np.arange(2), 10, 10)
    with pytest.raises(ValueError):
        n18.repeat_rows(np.random.default_rng(1), np.arange(100), 10, 200)


SIZES = {"E7.25": 2417, "E8.0": 2963}


def _gate_stages(profiles):
    rng = np.random.default_rng(21)
    return {"E7.25": _stage(SIZES["E7.25"], rng, profiles, [0.6, 0.3, 0.1], 60.0),
            "E8.0": _stage(SIZES["E8.0"], rng, profiles, [0.2, 0.3, 0.5], 240.0)}


def _patch_gate(n18, monkeypatch, profiles, flat=False):
    stages = _gate_stages(profiles)
    reads = []

    def load_stage(path, genes=None, drop_duplicates=False):
        name = Path(path).name.removesuffix(".h5ad")
        reads.append((name, list(genes), drop_duplicates))
        return stages[name][:, genes].copy()

    monkeypatch.setattr(n18, "load_stage", load_stage)
    monkeypatch.setattr(n18, "load_boards", lambda: {
        "T2:embryo:val_interp": Board(key="T2:embryo:val_interp", task="T2", genes=list(GENES), needs_coords=True)})
    monkeypatch.setattr(n18.harness, "_sha256", lambda p: "0" * 64)
    train, calls = _oracle_trainer(profiles, flat=flat)
    monkeypatch.setattr(n18.t2_latent, "train_nbvae", train)
    return reads, calls


def _numbers(text: str) -> set[int]:
    return {int(d) for d in re.findall(r"\d+", re.sub(r"(?<=\d),(?=\d{3})", "", text))}


def _gate_counts():
    """Every cell count the gate must not print or write: stage sizes, held-out and training sizes, half sizes."""
    counts = set(SIZES.values())
    for n in SIZES.values():
        k = int(round(0.2 * n))
        counts |= {k, n - k, k // 2, (n - k) // 2}
    return counts


def test_the_gate_passes_a_good_generator_on_the_two_brackets_only(n18, tmp_path, monkeypatch, capsys, own_repo):
    profiles = _profiles(np.random.default_rng(20))
    reads, calls = _patch_gate(n18, monkeypatch, profiles)
    out = tmp_path / "results"
    n18.main(["--out", str(out), "--data", str(tmp_path / "raw")])
    printed = capsys.readouterr().out
    # the backtest's two brackets, as the harness loads them, and never its target
    assert sorted(r[0] for r in reads) == ["E7.25", "E8.0"] and all(r[1] == list(GENES) and r[2] for r in reads)
    # one model, on the training cells of both brackets
    expected = sum(n - int(round(0.2 * n)) for n in SIZES.values())
    assert calls == [(expected, G)]
    res = json.loads((out / n18.NAMES["json"]).read_text())
    assert res["verdict"] == "pass" and printed.rstrip().endswith("verdict: pass")
    assert set(res["checks"]) == {"before", "after"} and set(res["checks"]["before"]) == set(n18.METRICS)
    for per in res["checks"].values():
        assert per["variogram"]["statistic"] == "ratio_to_control"
        assert per["mmd_u"]["statistic"] == "share_of_bracket_gap"
        assert per["variogram"]["ratio_to_control"] < 1.3 and per["mmd_u"]["share_of_bracket_gap"] < 0.1
    # provenance: a commit, a clean flag and input hashes; each output stamped
    prov = res["provenance"]
    assert prov["commit"] and prov["tree_clean"] is True and res["inputs"]["before"]["sha256"] == "0" * 64
    for name in n18.NAMES.values():
        assert stamp_path(out / name).exists()
    # no cell count anywhere, printed or written
    text = (out / n18.NAMES["text"]).read_text()
    bad = _gate_counts()
    for blob in (printed, text, (out / n18.NAMES["json"]).read_text()):
        assert not _numbers(blob) & bad


def test_the_gate_fails_a_generator_without_the_cell_state(n18, tmp_path, monkeypatch, capsys, own_repo):
    profiles = _profiles(np.random.default_rng(20))
    _patch_gate(n18, monkeypatch, profiles, flat=True)
    out = tmp_path / "results"
    n18.main(["--out", str(out), "--data", str(tmp_path / "raw")])
    res = json.loads((out / n18.NAMES["json"]).read_text())
    assert res["verdict"] == "fail" and capsys.readouterr().out.rstrip().endswith("verdict: fail")
    assert not all(ck["pass"] for per in res["checks"].values() for ck in per.values())


def test_the_gate_refuses_to_overwrite_its_outputs(n18, tmp_path, monkeypatch, own_repo):
    profiles = _profiles(np.random.default_rng(20))
    reads, calls = _patch_gate(n18, monkeypatch, profiles)
    out = tmp_path / "results"
    out.mkdir()
    (out / "n18_latent_fidelity.json").write_text("{}\n")
    with pytest.raises(SystemExit, match="rerun"):
        n18.main(["--out", str(out), "--data", str(tmp_path / "raw")])
    assert not reads and not calls  # nothing was read or trained


def test_pick_device_is_cpu_unless_a_gpu_is_free(n18, monkeypatch):
    assert n18.pick_device("cpu") == "cpu"
    assert n18.pick_device("cuda:3") == "cuda:3"
    fake_torch = SimpleNamespace(cuda=SimpleNamespace(is_available=lambda: True, device_count=lambda: 3,
                                                       mem_get_info=lambda i: (10, 100) if i < 2 else (95, 100)))
    monkeypatch.setitem(sys.modules, "torch", fake_torch)
    assert n18.pick_device("auto") == "cuda:2"
    fake_busy = SimpleNamespace(cuda=SimpleNamespace(is_available=lambda: True, device_count=lambda: 1,
                                                     mem_get_info=lambda i: (10, 100)))
    monkeypatch.setitem(sys.modules, "torch", fake_busy)
    assert n18.pick_device("auto") == "cpu"


# ------------------------------------------------------------------------------------------------ E14's smoke abort
@pytest.fixture(scope="module")
def smoke():
    spec = importlib.util.spec_from_file_location("e14_smoke", REPO / "experiments" / "harness" / "e14_smoke.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def _smoke_rows(mmd_diffs, variogram_diffs):
    """The run's rows with fgw at fixed skills and fgw_latent at those skills plus the given differences."""
    base = {"mmd_u": 0.67, "variogram": 0.66, "de_score": 0.5}
    return [{"seed": s,
             "skills": {"fgw": dict(base), "fgw_latent": {**base, "mmd_u": 0.67 + dm, "variogram": 0.66 + dv}}}
            for s, (dm, dv) in enumerate(zip(mmd_diffs, variogram_diffs, strict=True))]


def _write_smoke_run(folder, rows, **over):
    """A folder as run.py run writes it for the E7.5 backtest's smoke (the fields e14_smoke.py reads)."""
    run = {"task": {"key": "T2e/E7.25+E8.0->E7.5"}, "seeds": [0, 1], "candidates": ["floor", "fgw", "fgw_latent"],
           "rows": rows, "provenance": {"commit": "abc123", "tree_clean": True}}
    run.update(over)
    folder.mkdir(parents=True)
    (folder / "T2e_E7.25_E8.0_E7.5.json").write_text(json.dumps(run))
    return folder


def test_the_smoke_rule_aborts_only_past_the_rows_margin_on_either_metric(smoke):
    assert smoke.smoke_rule(_smoke_rows([0.01, 0.02], [0.0, -0.01]))["abort"] is False
    out = smoke.smoke_rule(_smoke_rows([-0.04, -0.05], [0.0, 0.0]))
    assert out["abort"] is False and out["checks"]["mmd_u"]["mean_difference"] == pytest.approx(-0.045)
    assert out["checks"]["mmd_u"]["difference_per_seed"] == pytest.approx([-0.04, -0.05])
    assert smoke.smoke_rule(_smoke_rows([-0.06, -0.05], [0.0, 0.0]))["abort"] is True  # mean -0.055 on mmd_u
    assert smoke.smoke_rule(_smoke_rows([0.0, 0.0], [-0.08, -0.03]))["abort"] is True  # mean -0.055 on the variogram
    with pytest.raises(ValueError, match="seeds"):
        smoke.smoke_rule(_smoke_rows([0.0], [0.0]))  # seed 1 is missing
    rows = _smoke_rows([0.0, 0.0], [0.0, 0.0])
    del rows[1]["skills"]["fgw_latent"]
    with pytest.raises(ValueError, match="fgw_latent"):
        smoke.smoke_rule(rows)
    rows = _smoke_rows([0.0, 0.0], [0.0, 0.0])
    rows[0]["skills"]["fgw_latent"]["variogram"] = float("nan")
    with pytest.raises(ValueError, match="non-finite"):
        smoke.smoke_rule(rows)


def test_the_smoke_script_decides_stamps_and_refuses(smoke, tmp_path, capsys, own_repo):
    ok = _write_smoke_run(tmp_path / "ok", _smoke_rows([0.01, -0.02], [0.03, 0.0]))
    out = tmp_path / "verdict"
    assert smoke.main(["--run", str(ok), "--out", str(out)]) == 0
    printed = capsys.readouterr().out
    assert printed.rstrip().endswith("verdict: continue")
    res = json.loads((out / smoke.NAMES["json"]).read_text())
    assert res["abort"] is False and res["run"]["commit"] == "abc123" and res["provenance"]["tree_clean"] is True
    assert res["run"]["sha256"] == harness._sha256(ok / smoke.RUN_FILE)
    for name in smoke.NAMES.values():
        assert stamp_path(out / name).exists()
    with pytest.raises(SystemExit, match="rerun"):  # it never overwrites
        smoke.main(["--run", str(ok), "--out", str(out)])
    bad = _write_smoke_run(tmp_path / "bad", _smoke_rows([-0.2, -0.2], [0.0, 0.0]))
    out2 = tmp_path / "verdict2"
    assert smoke.main(["--run", str(bad), "--out", str(out2)]) == smoke.ABORT_EXIT
    assert capsys.readouterr().out.rstrip().endswith("verdict: abort")
    assert json.loads((out2 / smoke.NAMES["json"]).read_text())["abort"] is True
    for over, why in (({"provenance": {"commit": "x", "tree_clean": False}}, "dirty"),
                      ({"seeds": [0, 1, 2]}, "seeds"),
                      ({"candidates": ["floor", "fgw"]}, "fgw_latent"),
                      ({"task": {"key": "T2e/E6.75+E7.5->E7.25"}}, "run")):
        folder = _write_smoke_run(tmp_path / f"refused_{why}", _smoke_rows([0.0, 0.0], [0.0, 0.0]), **over)
        with pytest.raises(SystemExit, match=why):
            smoke.main(["--run", str(folder), "--out", str(tmp_path / f"out_{why}")])


# ------------------------------------------------------------------------------------------------ the real training
def test_scvi_training_trains_decodes_and_is_deterministic_on_cpu():
    """Skipped in the project's main environment (no torch, no scvi-tools); runs with the gpu dependency group."""
    pytest.importorskip("torch")
    pytest.importorskip("scvi")
    rng = np.random.default_rng(30)
    profiles = _profiles(rng)
    a = _stage(400, rng, profiles, [0.5, 0.3, 0.2], 120.0)
    counts = sparse.csr_matrix(a.layers["counts"])
    m1 = t2_latent.train_nbvae(counts, 123, max_epochs=3, threads=2, device="cpu")
    m2 = t2_latent.train_nbvae(counts, 123, max_epochs=3, threads=2, device="cpu")
    assert m1.latent.shape == (400, t2_latent.N_LATENT)
    np.testing.assert_allclose(m1.latent, m2.latent, atol=1e-5)  # one seed, one model
    assert m1.record["gene_likelihood"] == "nb" and m1.record["n_latent"] == t2_latent.N_LATENT
    lib = np.asarray(counts.sum(1)).ravel().astype(np.float64)
    mu, theta = m1.decode(m1.latent, lib)
    assert mu.shape == (400, G) and theta.shape == (G,) and (theta > 0).all()
    assert mu.sum(1) == pytest.approx(lib, rel=1e-3)  # the NB means sum to the library they were decoded at
