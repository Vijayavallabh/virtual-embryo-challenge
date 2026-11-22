"""Try T65 (experiments/harness/TRIES.md): experiments/t1_m1/ms_read.py and its launcher on synthetic stages and
hand-made task files only. No released file, statistics file, board file or score is read; the committed chamber config
is read for its labels and the committed T62 task folder for its file names. It reuses row T62's and T64's fixtures.
The script reaches the builder (build_test) through recomputed_read and bets_read."""

from __future__ import annotations

import importlib.util
import inspect
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy import sparse
from test_g_mom_verify import _small_inputs
from test_t1_scope import GENES, HELD
from test_t1_state_transport import _gmom_ctx

from vec import harness, t1_scope
from vec.core import counts_of, mean_cp10k
from vec.io import make_submission
from vec.t1_scope import classify, load_spec, reweight

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "experiments" / "t1_m1" / "ms_read.py"
LAUNCHER = ROOT / "experiments" / "t1_m1" / "t65_ms.sh"
SHARES = {"cm": 0.6, "endocardium": 0.25, "blood": 0.15}
DIRECTED = ("CM", "Endocardium", "Endothelium")
HELD_M1 = sorted({*HELD, GENES[3], GENES[7]})
WITHIN = {"surrounding": ["Foregut", "Epithelium"]}
C = "c" * 40


@pytest.fixture
def ms(monkeypatch, tmp_path):
    spec = importlib.util.spec_from_file_location("ms_read_module", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    monkeypatch.setattr(mod, "RESULTS", tmp_path / "results")
    return mod


def _dense(x):
    return np.asarray(x.todense()) if hasattr(x, "todense") else np.asarray(x)


def _classified():
    ref, start, _ = _small_inputs()
    classify(ref, source="E9.5_RNA")
    classify(start, source="E8.5_RNA")
    return ref, start


def _directions(n_genes, groups=DIRECTED, size=0.4):
    rng = np.random.default_rng(11)
    return {g: rng.normal(0.0, size, n_genes) for g in groups}


def _fb():
    spec = load_spec()
    return {"E8.5": spec.fallbacks("E8.5_RNA"), "E9.5": spec.fallbacks("E9.5_RNA")}


def _without(stage, labels, keep=0):
    """The stage with all but `keep` of each of its `labels`' cells removed."""
    out = stage
    for label in labels:
        out = _without_one(out, label, keep)
    return out


def _without_one(stage, label, keep):
    lab = stage.obs["celltype"].astype(str).to_numpy()
    rows = np.flatnonzero(lab == label)
    keep_rows = np.sort(np.concatenate([np.flatnonzero(lab != label), rows[:keep]]))
    return stage[keep_rows].copy()


def _cut_draw(ctx):
    """The draw build_submission makes on a harness cut (the same reweight call), in row order."""
    ref = harness._shallow(ctx.inputs[ctx.task.reference_role])
    src = ctx.task.role_stems[ctx.task.reference_role]
    classify(ref, source=src)
    shares = harness.t1_chamber_shares(ref, src)
    pop, _ = reweight(ref, shares, harness.count_for(ctx.board, ref.n_obs), ctx.seed,
                      fallback=load_spec().fallbacks(src))
    return pop


def _lib(ms, e95, e85, seed=4, dirs=None, fb=None, fit_start=None):
    fit_start = e85 if fit_start is None else fit_start  # g_mom's fit; the momentum's start is e85
    factors = ms.br.gmom_fits(e95, "E9.5_RNA", fit_start, "E8.5_RNA")["file"]
    dirs = _directions(e95.n_vars) if dirs is None else dirs
    fb = _fb() if fb is None else fb
    arms, recs = ms.library_arms(e85, e95, SHARES, factors, HELD, HELD_M1, seed, WITHIN, fb, dirs, 0.2, n_cells=300)
    base0, _ = reweight(e95, SHARES, 300, seed, within=WITHIN, fallback=fb["E9.5"])
    return arms, recs, base0, factors, dirs


def _cut(ms, seed, start_filter=None):
    ref, start, _ = _small_inputs()
    fits = ms.br.gmom_fits(ref, "E9.5_RNA", start, "E8.5_RNA")  # g_mom's fit on the whole start
    if start_filter is not None:
        start = start_filter(start)
    ctx = _gmom_ctx(ref, start, seed=seed)
    dirs = _directions(ref.n_vars)
    arms, recs = ms.cut_arms(ctx, fits["file"], dirs, 0.2, HELD_M1)
    return arms, recs, ctx, fits, dirs


# ---------------------------------------------------------------------------------------------------------------------
# (a) only covered cells change, on both routes and in both arms


def test_on_the_cut_route_only_covered_cells_change(ms):
    arms, recs, ctx, _, _ = _cut(ms, 1000002)
    assert tuple(arms) == ms.ARMS
    labels = _cut_draw(ctx).obs["celltype"].astype(str).to_numpy()
    for arm in ms.COUNTED:
        rec = recs[arm]
        cov = np.isin(labels, rec["covered_set"])
        differs = (_dense(arms["m1rc_gmom"].X) != _dense(arms[arm].X)).any(axis=1)
        assert rec["covered_set"] == rec["builder_covered_set"] == ["Endocardium", "V-CM"]
        assert rec["non_covered_equal"] is True and differs.any() and not differs[~cov].any()
        assert int(differs.sum()) == rec["covered_cells_changed"] and rec["covered_cells"] == int(cov.sum())
        assert any(r["type"] == "cm" and r["status"] == "shifted" for r in rec["shift"])
    assert recs["m1rc_gmom_mom"]["perm"] is None and recs["m1rc_gmom_momperm"]["perm"]
    assert _dense(arms["m1rc_gmom_mom"].X).tobytes() != _dense(arms["m1rc_gmom_momperm"].X).tobytes()


def test_on_the_library_route_only_covered_cells_change(ms):
    e95, e85 = _classified()
    arms, recs, base0, _, _ = _lib(ms, e95, e85)
    labels = base0.obs["celltype"].astype(str).to_numpy()
    cov = np.isin(labels, ["V-CM", "Endocardium"])
    a = _dense(counts_of(arms["m1rc_gmom"]))
    for arm in ms.COUNTED:
        differs = (a != _dense(counts_of(arms[arm]))).any(axis=1)
        assert differs.any() and not differs[~cov].any() and recs[arm]["non_covered_equal"] is True
        assert int(differs.sum()) == recs[arm]["covered_cells_changed"]


def test_a_cell_outside_the_covered_set_that_differs_stops_the_run(ms):
    e95, _ = _classified()
    pop, _ = reweight(e95, SHARES, 300, 4, within=WITHIN)
    X = _dense(pop.X).astype(np.float32)
    Y = X.copy()
    Y[np.flatnonzero(pop.obs["celltype"].astype(str).to_numpy() == "Blood")[0], 0] += 1.0
    groups = ms.feng_dirs.label_groups(sorted(set(pop.obs["celltype"].astype(str))))
    dirs = _directions(e95.n_vars)
    assert ms.covered_record(pop, groups, dirs, X, X.copy(), "m1rc_gmom_mom", None)["non_covered_equal"] is True
    Z = X.copy()
    Z[np.flatnonzero(pop.obs["celltype"].astype(str).to_numpy() == "V-CM")[0], 0] += 1.0  # a covered cell may move
    assert ms.covered_record(pop, groups, dirs, X, Z, "m1rc_gmom_mom", None)["covered_cells_changed"] == 1
    with pytest.raises(SystemExit, match="outside the covered set"):
        ms.covered_record(pop, groups, dirs, X, Y, "m1rc_gmom_mom", None)


# ---------------------------------------------------------------------------------------------------------------------
# (b) the momentum by hand, before the atlas step


def _hand_momentum(Cm, cats, S, s_cats, R, r_cats, held_genes, seed, shifted):
    """B0's damped per-category momentum on counts Cm by hand: for each category in `shifted`, mean CP10k over its own
    cells at each stage, f = ((m_ref + 0.5) / (m_start + 0.5)) ** 0.5; f = 1 on held genes and for every other
    category; the keep_depth constant; from np.random.default_rng(seed) a Binomial over every entry with f < 1 in CSR
    order, then a Poisson top-up over every entry with f > 1; a cell left empty gets its counts back."""
    def mean_cp(M):
        M = np.asarray(M, dtype=np.float64)
        return (M / M.sum(axis=1, keepdims=True) * 1e4).mean(axis=0)

    F = {}
    for c in sorted(set(cats)):
        f = (((mean_cp(R[r_cats == c]) + 0.5) / (mean_cp(S[s_cats == c]) + 0.5)) ** 0.5 if c in shifted
             else np.ones(R.shape[1]))
        f[held_genes] = 1.0
        F[c] = f
    M = sparse.csr_matrix(np.asarray(Cm, dtype=np.int64))
    M.sort_indices()
    row = np.repeat(np.arange(M.shape[0]), np.diff(M.indptr))
    c = M.data.astype(np.float64)
    f = np.array([F[cats[i]][j] for i, j in zip(row, M.indices, strict=True)])
    free = ~held_genes[M.indices]
    num = np.bincount(row, weights=c * free, minlength=M.shape[0])
    den = np.bincount(row, weights=c * f * free, minlength=M.shape[0])
    scale = np.divide(num, den, out=np.ones(M.shape[0]), where=den > 0)
    f = np.where(free, f * scale[row], 1.0)
    rng = np.random.default_rng(seed)
    new = M.data.copy()
    down, up = f < 1.0, f > 1.0
    new[down] = rng.binomial(M.data[down], f[down])
    new[up] = M.data[up] + rng.poisson((f[up] - 1.0) * M.data[up])
    dead = np.bincount(row, weights=new, minlength=M.shape[0]) == 0
    new[dead[row]] = M.data[dead[row]]
    return _dense(sparse.csr_matrix((new, M.indices, M.indptr), shape=M.shape))


def test_the_covered_cells_momentum_is_the_hand_computation(ms, monkeypatch):
    e95, e85 = _classified()
    base0, _ = reweight(e95, SHARES, 300, 4, within=WITHIN, fallback=_fb()["E9.5"])
    dirs = _directions(e95.n_vars)
    groups = ms.feng_dirs.label_groups(sorted(set(e95.obs["celltype"].astype(str))))
    monkeypatch.setattr(t1_scope, "external_shift", lambda cells, *a, **k: (cells, None))  # the momentum alone
    shift_kw = {"within": WITHIN, "start_fallback": _fb()["E8.5"], "ref_fallback": _fb()["E9.5"], "min_cells": 20,
                "eps": 0.5}
    pop, rep = ms.ms_population(base0, e85, e95, dirs, 0.2, groups, HELD_M1, 4, shift_kw)
    labels = base0.obs["celltype"].astype(str).to_numpy()
    idx = np.flatnonzero(np.isin(labels, ["V-CM", "Endocardium"]))
    shifted = {r["type"] for r in rep.rows if r["status"] == "shifted"}
    assert shifted == {"cm", "endocardium"}  # the synthetic start's Endothelium cells are E8.5's endocardium
    held_genes = np.isin(np.asarray(e95.var_names), HELD_M1)
    want = _hand_momentum(_dense(counts_of(base0))[idx], base0.obs["category"].astype(str).to_numpy()[idx],
                          _dense(counts_of(e85)), e85.obs["category"].astype(str).to_numpy(),
                          _dense(counts_of(e95)), e95.obs["category"].astype(str).to_numpy(), held_genes, [4, 1],
                          shifted)
    got = _dense(counts_of(pop))[idx]
    assert np.array_equal(got, want)
    drawn = _dense(counts_of(base0))[idx]
    assert not np.array_equal(got, drawn) and np.array_equal(got[:, held_genes], drawn[:, held_genes])
    rest = np.setdiff1d(np.arange(base0.n_obs), idx)
    assert _dense(pop.X)[rest].tobytes() == _dense(base0.X)[rest].tobytes()


def test_the_atlas_step_takes_the_moved_cells_on_its_own_stream(ms, monkeypatch):
    real, calls = t1_scope.external_shift, []

    def spy(cells, directions, fraction, groups, held=(), **kw):
        calls.append({"cells": cells, "fraction": fraction, "held": list(held), "seed": kw.get("seed")})
        return real(cells, directions, fraction, groups, held, **kw)

    monkeypatch.setattr(t1_scope, "external_shift", spy)
    e95, e85 = _classified()
    base0, _ = reweight(e95, SHARES, 300, 4, within=WITHIN, fallback=_fb()["E9.5"])
    groups = ms.feng_dirs.label_groups(sorted(set(e95.obs["celltype"].astype(str))))
    shift_kw = {"within": WITHIN, "start_fallback": _fb()["E8.5"], "ref_fallback": _fb()["E9.5"]}
    ms.ms_population(base0, e85, e95, _directions(e95.n_vars), 0.2, groups, HELD_M1, 4, shift_kw)
    assert len(calls) == 1 and calls[0]["seed"] == [4, 2] and calls[0]["fraction"] == 0.2
    assert calls[0]["held"] == HELD_M1 and "counts" in calls[0]["cells"].layers  # momentum's output, not the draw
    idx = np.flatnonzero(np.isin(base0.obs["celltype"].astype(str).to_numpy(), ["V-CM", "Endocardium"]))
    assert list(calls[0]["cells"].obs_names) == list(base0.obs_names[idx])


# ---------------------------------------------------------------------------------------------------------------------
# (c) the base is T62's candidate as called; with the momentum idle the arms are the base


def test_the_base_is_t62s_candidate_as_called_on_both_routes(ms):
    arms, _, ctx, fits, dirs = _cut(ms, 1000003)
    t62 = ms.rc.cut_pair(ctx, fits["file"], dirs, 0.2, HELD_M1)
    assert arms["m1rc_gmom"].X.tobytes() == t62["m1rc_gmom"].X.tobytes()
    e95, e85 = _classified()
    arms, _, _, factors, dirs = _lib(ms, e95, e85, seed=5)
    t62 = ms.rc.library_pair(e95, SHARES, factors, HELD, HELD_M1, 5, WITHIN, _fb(), dirs, 0.2, n_cells=300)
    assert _dense(arms["m1rc_gmom"].X).tobytes() == _dense(t62["m1rc_gmom"].X).tobytes()


def test_with_every_covered_category_skipped_the_arms_are_the_base(ms):
    e95, e85 = _classified()
    arms, recs, _, _, _ = _lib(ms, e95, _without(e85, ("V-CM", "Endothelium"), keep=5), fit_start=e85)
    for arm in ms.COUNTED:
        assert ms.prediction_sha(arms["m1rc_gmom"].X) == ms.prediction_sha(arms[arm].X)
        assert recs[arm]["covered_cells_changed"] == 0
    arms, recs, _, _, _ = _cut(ms, 1000005, start_filter=lambda s: _without(s, ("V-CM", "Endothelium"), keep=5))
    for arm in ms.COUNTED:
        assert ms.prediction_sha(arms["m1rc_gmom"].X) == ms.prediction_sha(arms[arm].X)
        assert recs[arm]["covered_cells_changed"] == 0 and recs[arm]["covered_cells"] > 0


def test_the_cut_arm_is_the_same_population_as_the_library_route_on_the_cut_draw(ms):
    """On a cut the hook receives build_submission's draw; the arm's population before g_mom is ms_population on that
    draw with T61's b0_gmom momentum settings (the chamber config's fallbacks for both file stems)."""
    ref, start, _ = _small_inputs()
    ctx = _gmom_ctx(ref, start, seed=1000007)
    fits = ms.br.gmom_fits(ref, "E9.5_RNA", start, "E8.5_RNA")
    dirs = _directions(ref.n_vars)
    arms, _ = ms.cut_arms(ctx, fits["file"], dirs, 0.2, HELD_M1)
    draw = _cut_draw(ctx)
    cref, cstart = harness._shallow(ref), harness._shallow(start)
    classify(cref, source="E9.5_RNA")
    classify(cstart, source="E8.5_RNA")
    groups, covered, with_direction = ms.rc.direction_kit(cref.obs["celltype"].astype(str).unique(), dirs)
    spec = load_spec()
    pop, _ = ms.ms_population(draw, cstart, cref, dirs, 0.2, groups, HELD_M1, ctx.seed,
                              {"start_fallback": spec.fallbacks("E8.5_RNA"),
                               "ref_fallback": spec.fallbacks("E9.5_RNA")})
    moved = ms.rc.rc_hook(fits["file"], ctx.seed, HELD, with_direction, covered)(pop)
    cols = moved.var_names.get_indexer(list(ctx.board.genes))
    sub = make_submission(_dense(moved.X[:, cols]), list(ctx.board.genes))
    want = harness.finish(ctx, cref, ms.rc.dense(sub.X), None)
    assert arms["m1rc_gmom_mom"].X.tobytes() == want.X.tobytes()


# ---------------------------------------------------------------------------------------------------------------------
# (d) the placebo


def test_the_placebo_permutes_within_bins_of_non_held_genes_by_hand(ms):
    rng = np.random.default_rng(3)
    n = 60
    logf = rng.normal(0, 0.3, n)
    held = np.zeros(n, dtype=bool)
    held[[2, 17, 40]] = True
    logf[held] = 0.0
    counts = sparse.csr_matrix(rng.poisson(rng.gamma(1.0, 2.0, n), size=(25, n)))
    filt, records = ms.perm_filter()
    out, rec = filt(logf, counts, held)
    free = np.flatnonzero(~held)
    basis = mean_cp10k(counts)[free]
    hand = logf[free].copy()
    g = np.random.default_rng(20261110)
    for genes in np.array_split(np.argsort(basis, kind="stable"), 10):
        hand[genes] = logf[free][genes[g.permutation(len(genes))]]
    assert np.array_equal(out[free], hand) and np.all(out[held] == 0.0)
    assert rec["same_values"] is True and rec["seed"] == 20261110 and rec["bins"] == 10 and records == [rec]
    assert rec["sum_abs_before"] == pytest.approx(rec["sum_abs_after"]) and rec["moved_genes"] > 0
    out2, _ = ms.perm_filter()[0](logf, counts, held)
    assert np.array_equal(out, out2)  # one generator per category, the same seed each time
    assert ms.PERM_SEED == harness.LP_PERM_SEED == 20261110 and ms.PERM_BINS == harness.LP_BINS == 10


def test_the_placebo_arm_applies_the_permuted_factor(ms, monkeypatch):
    real, calls = t1_scope.momentum_shift, []

    def spy(start, ref, target_like, damping, held_flat=(), **kw):
        calls.append(kw.get("logf_filter"))
        return real(start, ref, target_like, damping, held_flat, **kw)

    monkeypatch.setattr(t1_scope, "momentum_shift", spy)
    e95, e85 = _classified()
    arms, recs, base0, _, _ = _lib(ms, e95, e85, seed=6)
    assert calls[0] is None and callable(calls[1]) and len(calls) == 2
    perm = recs["m1rc_gmom_momperm"]["perm"]
    assert perm and all(p["same_values"] is True for p in perm)
    assert len(perm) == sum(r["status"] == "shifted" for r in recs["m1rc_gmom_momperm"]["shift"])


# ---------------------------------------------------------------------------------------------------------------------
# (e) streams and fit inputs; (f) the covered set on the committed config


def test_the_momentums_stream_and_fit_inputs(ms, monkeypatch):
    real, calls = t1_scope.momentum_shift, []

    def spy(start, ref, target_like, damping, held_flat=(), **kw):
        calls.append({"start": start, "ref": ref, "damping": damping, "seed": kw.get("seed"), "held": list(held_flat),
                      "kw": {k: v for k, v in kw.items() if k != "logf_filter"},
                      "labels": set(target_like.obs["celltype"].astype(str))})
        return real(start, ref, target_like, damping, held_flat, **kw)

    monkeypatch.setattr(t1_scope, "momentum_shift", spy)
    e95, e85 = _classified()
    _lib(ms, e95, e85, seed=6)
    assert len(calls) == 2 and all(c["start"] is e85 and c["ref"] is e95 for c in calls)
    assert all(c["seed"] == [6, 1] and c["damping"] == 0.5 and c["held"] == HELD_M1 for c in calls)
    assert calls[0]["kw"]["min_cells"] == 20 and calls[0]["kw"]["eps"] == 0.5 and calls[0]["kw"]["within"] == WITHIN
    assert all(c["labels"] <= {"V-CM", "Endocardium"} for c in calls)  # the covered cells alone
    calls.clear()
    ref, start, target = _small_inputs()
    ctx = _gmom_ctx(ref, start, seed=1000006)
    fits = ms.br.gmom_fits(ref, "E9.5_RNA", start, "E8.5_RNA")
    ms.cut_arms(ctx, fits["file"], _directions(ref.n_vars), 0.2, HELD_M1)
    assert len(calls) == 2 and all(c["seed"] == [1000006, 1] and c["damping"] == 0.5 for c in calls)
    spec = load_spec()
    assert calls[0]["kw"]["start_fallback"] == spec.fallbacks("E8.5_RNA")
    assert calls[0]["kw"]["ref_fallback"] == spec.fallbacks("E9.5_RNA")
    assert "within" not in calls[0]["kw"] and "min_cells" not in calls[0]["kw"]  # build_submission's defaults: 20, 0.5
    assert list(calls[0]["start"].obs_names) == list(start.obs_names)
    assert list(calls[0]["ref"].obs_names) == list(ref.obs_names) != list(target.obs_names)
    for fn in (ms.cut_arms, ms.library_arms, ms.ms_population):  # no truth reaches an arm
        assert not {"truth", "target", "e105"} & set(inspect.signature(fn).parameters)
    assert ms.MOMENTUM_STREAM == 1 and ms.DIRECTION_STREAM == 2 and ms.DAMPING == 0.5 == harness.B0_DAMPING


E95_CHAMBER = ("AVC-CM", "IFT-CM", "OFT/RV-CM", "SV-CM", "V-CM", "Endocardium", "BEC", "Proepicardium", "Blood")
RECOMPUTED = ("CM", "Endocardium", "Endothelium", "Epicardium")


def test_the_covered_set_is_the_builders_and_on_e95s_chamber_labels_it_is_all_but_blood(ms):
    import yaml

    cfg = yaml.safe_load((ROOT / "configs" / "t1_chamber_classes.yaml").read_text())
    dirs = {g: np.zeros(3) for g in RECOMPUTED}
    for labs in (E95_CHAMBER, sorted(cfg["labels"]), ("V-CM", "Blood", "Mystery cells")):
        groups = ms.feng_dirs.label_groups(sorted(labs))
        cov = sorted(lab for lab in labs if ms.covered_mask([lab], groups, dirs)[0])
        assert cov == ms.builder_covered(labs, dirs)
    assert ms.builder_covered(E95_CHAMBER, dirs) == sorted(set(E95_CHAMBER) - {"Blood"})


# ---------------------------------------------------------------------------------------------------------------------
# The reading on hand-made task files


def _arm(score):
    return {"score": score, "groups": {g: score / 100 for g in ("de_recovery", "direction", "distribution",
                                                               "covariance")}}


def _relabel(ms, **change):
    return {"sha256": ms.rc.STATS_SHA256, "stats_sha256_hashed": ms.rc.STATS_SHA256, "yaml_sha256": ms.rc.YAML_SHA256,
            "route": "recomputed", "record_route": "recomputed", "authors_labels_read": False, **change}


def _covered(idle=False, perm=True):
    rec = {"covered_set": ["V-CM"], "builder_covered_set": ["V-CM"], "covered_cells": 30,
           "covered_cells_changed": 0 if idle else 30, "non_covered_equal": True,
           "shift": [{"type": "cm", "status": "shifted", "n_start": 300, "n_ref": 420}], "perm": None}
    if perm:
        rec["perm"] = [{"seed": 20261110, "bins": 10, "free_genes": 50, "moved_genes": 45, "sum_abs_before": 3.0,
                        "sum_abs_after": 3.0, "same_values": True}]
    return rec


def _task(ms, gain, lib=False, per_seed=None, perm_gain=0.0, base=50.0, commit=C, clean=True, seeds=range(10),
          inverted=False, relabel=None, idle=False, **top):
    rows = []
    for s in seeds:
        g = per_seed[s] if per_seed is not None else gain
        pg = perm_gain[s] if isinstance(perm_gain, list | tuple) else perm_gain
        arms = {"m1rc_gmom": _arm(base + s / 100), "m1rc_gmom_mom": _arm(base + s / 100 + g),
                "m1rc_gmom_momperm": _arm(base + s / 100 + pg)}
        sha = {"m1rc_gmom": "a" * 64, "m1rc_gmom_mom": "a" * 64 if idle else "b" * 64, "m1rc_gmom_momperm": "d" * 64}
        row = {"seed": s, "arms": arms, "sha256": sha,
               "covered": {"m1rc_gmom_mom": _covered(idle, perm=False), "m1rc_gmom_momperm": _covered()}}
        if lib:
            arms["copy_last"], sha["copy_last"] = _arm(40.0), "c" * 64
            row["inverted_anchors"] = ["x"] if inverted and s == 3 else []
        else:
            row["floor"] = 50.0
        rows.append(row)
    out = {"rows": rows, "stamp": {"commit": commit, "tree_clean": clean}, "fraction": 0.2, "damping": 0.5,
           "momentum_stream": 1, "direction_stream": 2, "perm_seed": 20261110, "perm_bins": 10,
           "relabel": relabel if relabel is not None else _relabel(ms)}
    out.update(top)
    return out


def _committed(base=50.0):
    return {k: {s: base + s / 100 for s in range(10)} for k in ("chamber", "whole", "1", "2", "3")}


def _loaded(ms, chamber=1.0, whole=0.5, libs=(0.0, 0.0, 0.0), chamber_perm=0.2, **kw):
    out = {"chamber": _task(ms, chamber, perm_gain=chamber_perm), "whole": _task(ms, whole)}
    out.update({lib: _task(ms, g, lib=True) for lib, g in zip(("1", "2", "3"), libs, strict=True)})
    out.update(kw)
    return out


OK = {"reasons": []}


def _p(mean, ahead=10):
    return {"mean": mean, "ahead": ahead}


def test_ahead(ms):
    r = ms.t65_reading(_loaded(ms), _committed(), OK, (), C)
    assert r["outcome"] == "MS is AHEAD" and r["ahead"] and r["follows"] == ms.FOLLOWS_AHEAD
    assert r["j6_range"] == pytest.approx([0.0, 1.0]) and r["valid_libraries"] == ["1", "2", "3"]
    assert r["versus"]["chamber"]["mean"] == pytest.approx(0.8)
    assert r["placebo"]["chamber"]["mean"] == pytest.approx(0.2)
    lines = ms.outcome_lines(r)
    assert lines[0] == "Outcome: MS is AHEAD" and all(b in lines for b in ms.BESIDE) and lines[-1] == ms.FOLLOWS_AHEAD


BEATS = _p(1.0)


@pytest.mark.parametrize("tasks,valid,beats,part", [
    ({"chamber": _p(0.3), "whole": _p(1.0)}, ["1", "2"], BEATS, "chamber cut's mean"),
    ({"chamber": _p(1.0, 7), "whole": _p(1.0)}, ["1", "2"], BEATS, "fewer than 8 of 10 seeds are above 0"),
    ({"chamber": _p(1.0), "whole": _p(0.0)}, ["1", "2"], BEATS, "whole cut's mean"),
    ({"chamber": _p(1.0), "whole": _p(-0.01)}, ["1", "2"], BEATS, "whole cut's mean"),
    ({"chamber": _p(1.0), "whole": _p(1.0), "2": _p(-0.31)}, ["1", "2"], BEATS, "below -0.3 (library 2)"),
    ({"chamber": _p(1.0), "whole": _p(1.0)}, ["3"], BEATS, "fewer than 2 library tasks"),
    ({"chamber": _p(1.0), "whole": _p(1.0)}, ["1", "2"], _p(1.0, 7), "beats m1rc_gmom_momperm on the chamber cut in "
                                                                     "fewer than 8"),
    ({"chamber": _p(1.0), "whole": _p(1.0)}, ["1", "2"], _p(0.0, 9), "minus m1rc_gmom_momperm is not positive"),
    ({"chamber": _p(1.0), "whole": _p(1.0)}, ["1", "2"], _p(-0.1, 8), "is not positive"),
])
def test_not_ahead(ms, tasks, valid, beats, part):
    tasks = {"1": _p(0.0), "2": _p(0.0), "3": _p(0.0), **tasks}
    ahead, failed = ms.statement(tasks, valid, beats)
    assert not ahead and any(part in f for f in failed)


def test_the_boundaries_that_pass(ms):
    tasks = {"chamber": _p(0.3000001, 8), "whole": _p(1e-9), "1": _p(-0.3), "2": _p(0.0), "3": _p(-5.0)}
    assert ms.statement(tasks, ["1", "2"], _p(1e-9, 8)) == (True, [])  # library 3 invalid: its -5 is ignored


def test_the_placebo_clause_end_to_end(ms):
    r = ms.t65_reading(_loaded(ms, chamber_perm=1.0), _committed(), OK, (), C)  # mom = momperm on every seed
    assert r["outcome"] == "MS is NOT AHEAD" and r["follows"] == ms.FOLLOWS_END
    assert any("fewer than 8" in f for f in r["failed"]) and any("not positive" in f for f in r["failed"])
    per = [0.5] * 7 + [1.5] * 3  # mom ahead of momperm on 7 seeds only, the mean positive
    loaded = _loaded(ms)
    loaded["chamber"] = _task(ms, 1.0, perm_gain=per)
    r = ms.t65_reading(loaded, _committed(), OK, (), C)
    assert r["outcome"] == "MS is NOT AHEAD" and r["versus"]["chamber"]["ahead"] == 7
    assert r["versus"]["chamber"]["mean"] > 0


def test_an_invalid_library_is_ignored_and_two_invalid_libraries_stop_it(ms):
    loaded = _loaded(ms, libs=(0.0, 0.0, -5.0))
    loaded["3"] = _task(ms, -5.0, lib=True, inverted=True)
    r = ms.t65_reading(loaded, _committed(), OK, (), C)
    assert r["outcome"] == "MS is AHEAD" and r["valid_libraries"] == ["1", "2"]
    loaded["2"] = _task(ms, 0.0, lib=True, inverted=True)
    r = ms.t65_reading(loaded, _committed(), OK, (), C)
    assert r["outcome"] == "MS is NOT AHEAD" and any("fewer than 2 library" in f for f in r["failed"])
    assert ms.outcome_lines(r)[-1] == ms.FOLLOWS_END


def test_a_negative_chamber_mean_gives_the_j6_range_mean_to_zero(ms):
    r = ms.t65_reading(_loaded(ms, chamber=-0.5, chamber_perm=-1.0), _committed(), OK, (), C)
    assert r["outcome"] == "MS is NOT AHEAD" and r["j6_range"] == pytest.approx([-0.5, 0.0])


def test_the_idle_case_is_named_and_read_as_it_is(ms):
    loaded = _loaded(ms)
    loaded["2"] = _task(ms, 0.0, lib=True, idle=True)
    r = ms.t65_reading(loaded, _committed(), OK, (), C)
    assert r["idle_tasks"] == {"m1rc_gmom_mom": ["2"], "m1rc_gmom_momperm": []} and r["outcome"] == "MS is AHEAD"
    assert "m1rc_gmom_mom moved no cell on 2: nothing of it is priced there" in ms.outcome_lines(r)


def test_no_line_or_text_carries_advice(ms):
    for loaded in (_loaded(ms), _loaded(ms, chamber=0.0)):
        r = ms.t65_reading(loaded, _committed(), OK, (), C)
        res = {"label": ms.LABEL, "reading": r, "tries": 10, "stamp": {"commit": "e", "tree_clean": True, "utc": "x"}}
        assert "advice" not in ms.t65_text(res).lower()
    assert "advice" not in SCRIPT.read_text().lower()


@pytest.mark.parametrize("case", ["missing_task", "missing_seed", "missing_arm", "missing_perm_arm", "missing_sha",
                                  "missing_copy_last", "two_commits", "not_row_commit", "unclean", "gate", "stats",
                                  "yaml", "route", "authors", "checks", "source", "fraction", "damping", "stream",
                                  "direction_stream", "perm_seed", "perm_bins", "perm_record", "non_covered",
                                  "covered_set"])
def test_every_void_reason_fires(ms, case):
    loaded, committed, checks, source, row_commit = _loaded(ms), _committed(), OK, [], C
    if case == "missing_task":
        loaded["1"] = None
    elif case == "missing_seed":
        loaded["whole"] = _task(ms, 0.5, seeds=range(9))
    elif case == "missing_arm":
        del loaded["chamber"]["rows"][2]["arms"]["m1rc_gmom_mom"]
    elif case == "missing_perm_arm":
        del loaded["3"]["rows"][2]["arms"]["m1rc_gmom_momperm"]
    elif case == "missing_sha":
        del loaded["whole"]["rows"][5]["sha256"]["m1rc_gmom"]
    elif case == "missing_copy_last":
        del loaded["3"]["rows"][0]["arms"]["copy_last"]
    elif case == "two_commits":
        loaded["1"] = _task(ms, 0.0, lib=True, commit="d" * 40)
    elif case == "not_row_commit":
        row_commit = "d" * 40
    elif case == "unclean":
        loaded["chamber"] = _task(ms, 1.0, perm_gain=0.2, clean=False)
    elif case == "gate":
        committed["2"][4] += 0.0201
    elif case == "stats":
        loaded["chamber"] = _task(ms, 1.0, perm_gain=0.2, relabel=_relabel(ms, stats_sha256_hashed="0" * 64))
    elif case == "yaml":
        loaded["whole"] = _task(ms, 0.5, relabel=_relabel(ms, yaml_sha256="0" * 64))
    elif case == "route":
        loaded["1"] = _task(ms, 0.0, lib=True, relabel=_relabel(ms, route="none"))
    elif case == "authors":
        loaded["2"] = _task(ms, 0.0, lib=True, relabel=_relabel(ms, authors_labels_read=True))
    elif case == "checks":
        checks = {"reasons": ["the settings file's sha256 is x"]}
    elif case == "source":
        source = ["a feng_directions call without the recomputed groups"]
    elif case == "fraction":
        loaded["chamber"] = _task(ms, 1.0, perm_gain=0.2, fraction=0.4)
    elif case == "damping":
        loaded["whole"] = _task(ms, 0.5, damping=1.0)
    elif case == "stream":
        loaded["3"] = _task(ms, 0.0, lib=True, momentum_stream=2)
    elif case == "direction_stream":
        loaded["3"] = _task(ms, 0.0, lib=True, direction_stream=1)
    elif case == "perm_seed":
        loaded["whole"] = _task(ms, 0.5, perm_seed=1)
    elif case == "perm_bins":
        loaded["whole"] = _task(ms, 0.5, perm_bins=5)
    elif case == "perm_record":
        loaded["2"]["rows"][3]["covered"]["m1rc_gmom_momperm"]["perm"][0]["same_values"] = False
    elif case == "non_covered":
        loaded["2"]["rows"][6]["covered"]["m1rc_gmom_momperm"]["non_covered_equal"] = False
    elif case == "covered_set":
        loaded["1"]["rows"][1]["covered"]["m1rc_gmom_mom"]["builder_covered_set"] = ["V-CM", "Blood"]
    r = ms.t65_reading(loaded, committed, checks, source, row_commit)
    assert r["outcome"] == "VOID" and r["reasons"] and r["follows"] == ms.FOLLOWS_END
    lines = ms.outcome_lines(r)
    assert lines[0].startswith("Outcome: VOID") and lines[-1] == ms.FOLLOWS_END
    assert not any(b in lines for b in ms.BESIDE)
    if case == "gate":
        assert "the gate failed" in lines[0] and "0.0201" in lines[0]


def test_the_gate_passes_at_a_difference_of_0_02(ms):
    loaded, committed = _loaded(ms), _committed()
    committed["2"][4] = 0.0
    for a in ms.ARMS:
        loaded["2"]["rows"][4]["arms"][a]["score"] = 0.02
    assert ms.t65_reading(loaded, committed, OK, (), C)["outcome"] != "VOID"
    loaded["2"]["rows"][4]["arms"]["m1rc_gmom"]["score"] = 0.0201
    assert ms.t65_reading(loaded, committed, OK, (), C)["outcome"] == "VOID"
    assert ms.REPRO_TOL == 0.02 == ms.br.REPRO_TOL


def test_committed_t62_reads_m1rc_gmom_per_seed(ms, tmp_path):
    (tmp_path / "lib1.json").write_text(json.dumps({"rows": [{"seed": s, "arms": {"m1rc_gmom": {"score": 60.0 + s}}}
                                                            for s in range(10)]}))
    assert ms.committed_t62("lib1.json", tmp_path) == {s: 60.0 + s for s in range(10)}
    assert ms.T62_TASKS == SCRIPT.parent / "results" / "t62_recomputed" / "tasks"
    assert all((ms.T62_TASKS / f).exists() for f in ms.task_files().values())  # the gate's committed files


def test_the_counts_and_constants(ms):
    assert ms.TRIES_COUNT == 10 and ms.TRY_RANGE == (418, 428)
    assert ms.COUNTED == ("m1rc_gmom_mom", "m1rc_gmom_momperm") and ms.GATED == ("m1rc_gmom",)
    assert ms.ARMS == ("m1rc_gmom", "m1rc_gmom_mom", "m1rc_gmom_momperm") and ms.LIB_ARMS[0] == "copy_last"
    assert ms.BASE == ms.rc.CANDIDATE and ms.ROW == "T65" and ms.OUT_FILES == ("t65.json", "t65.txt")
    assert ms.AHEAD_BAR == 0.3 and ms.MIN_AHEAD == 8 and ms.WHOLE_ABOVE == 0.0 and ms.LIB_FLOOR == -0.3
    assert ms.MIN_LIBRARIES == 2 and ms.MIN_BEATS == 8 and ms.PRACTICE_FRACTION == 0.2


# ---------------------------------------------------------------------------------------------------------------------
# Refusals


def _tries(tmp_path, monkeypatch, ms, row="T65"):
    tries = tmp_path / "TRIES.md"
    tries.write_text(f"| {row} (one row) | T1 | x |\n")
    monkeypatch.setattr(ms, "TRIES", tries)


def test_the_run_and_the_read_are_refused_without_the_row(ms, tmp_path, monkeypatch):
    _tries(tmp_path, monkeypatch, ms, row="T64")
    monkeypatch.setattr(ms.rc, "relabel_checks", lambda *a, **k: pytest.fail("nothing before the row's refusal"))
    monkeypatch.setattr(ms.gv, "tree_state", lambda: pytest.fail("nothing before the row's refusal"))
    with pytest.raises(SystemExit, match="no row"):
        ms.refuse(tmp_path / "out", "x.json", C)
    with pytest.raises(SystemExit, match="no row"):
        ms.read(tmp_path / "out", C)


def test_the_run_is_refused_off_the_run_commit_unclean_or_on_another_damping(ms, tmp_path, monkeypatch):
    _tries(tmp_path, monkeypatch, ms)
    monkeypatch.setattr(ms.gv, "tree_state", lambda: {"commit": C, "clean": True})
    monkeypatch.setattr(ms.rc, "relabel_checks", lambda *a, **k: {"passes": True, "reasons": []})
    assert ms.refuse(tmp_path / "out", "x.json", C)[0]["clean"] is True
    for commit in ("e" * 40, "", None):
        with pytest.raises(SystemExit, match="HEAD is not the run commit"):
            ms.refuse(tmp_path / "out", "x.json", commit)
    (tmp_path / "out").mkdir()
    (tmp_path / "out" / "x.json").write_text("{}")
    with pytest.raises(SystemExit, match="already holds"):
        ms.refuse(tmp_path / "out", "x.json", C)
    monkeypatch.setattr(ms, "DAMPING", 1.0)
    with pytest.raises(SystemExit, match="damping"):
        ms.refuse(tmp_path / "out", "y.json", C)
    monkeypatch.setattr(ms.gv, "tree_state", lambda: {"commit": C, "clean": False})
    with pytest.raises(SystemExit, match="not clean"):
        ms.refuse(tmp_path / "out", "y.json", C)
    monkeypatch.setattr(ms.gv, "tree_state", lambda: {"commit": C, "clean": True})
    monkeypatch.setattr(ms.rc, "relabel_checks", lambda *a, **k: {"passes": False, "reasons": ["no record"]})
    monkeypatch.setattr(ms, "DAMPING", 0.5)
    with pytest.raises(SystemExit, match="no record"):
        ms.refuse(tmp_path / "out", "y.json", C)


def test_the_source_guard(ms, tmp_path):
    assert ms.source_guard() == [] and ms.rc.source_guard(SCRIPT) == []
    for text in ('x = fd.feng_directions(list(g))\n', 'x = fd.feng_directions(g, groups="none")\n',
                 "import perm_read\n", "from heldmom_read import ARM\n", "import double_step_read\n",
                 "import n5" + "4_authors_agreement as a\n", "m = import" + "lib.import_module('x')\n",
                 "p = 'boards/index" + ".json'\n"):
        bad = tmp_path / "bad.py"
        bad.write_text(text)
        assert ms.source_guard(bad), text


def test_the_read_prints_no_score_and_never_overwrites(ms, tmp_path, monkeypatch, capsys):
    out = tmp_path / "results" / "t65"
    (out / "tasks").mkdir(parents=True)
    loaded = _loaded(ms, 1.234567, 0.5, (0.4, 0.4, 0.4))
    for k, f in ms.task_files().items():
        (out / "tasks" / f).write_text(json.dumps(loaded[k]))
    committed = _committed()
    _tries(tmp_path, monkeypatch, ms)
    monkeypatch.setattr(ms, "committed_t62", lambda f: committed[{v: k for k, v in ms.task_files().items()}[f]])
    monkeypatch.setattr(ms.rc, "relabel_checks", lambda **k: {"reasons": [], "passes": True})
    monkeypatch.setattr(ms.gv, "stamp", lambda started, files: {"commit": "e" * 40, "tree_clean": True, "utc": "x"})
    monkeypatch.setattr(ms.gv, "tree_state", lambda: {"commit": "e" * 40, "clean": True})
    with pytest.raises(SystemExit, match="row's commit"):
        ms.read(out, None)
    res = ms.read(out, C)
    printed = capsys.readouterr().out
    assert res["reading"]["outcome"] == "MS is AHEAD" and res["tries"] == 10
    assert "1.23" not in printed and "51." not in printed and "50." not in printed and "Outcome" in printed
    assert all(b in printed for b in ms.BESIDE) and "advice" not in printed.lower()
    text = (out / "t65.txt").read_text()
    assert "+1.23" in text and "J6 range" in text and "momentum cm: shifted" in text and "T1 418 to 428" in text
    assert "permutation: seed 20261110, 10 bins" in text and "m1rc_gmom_momperm - m1rc_gmom" in text
    with pytest.raises(SystemExit, match="already holds"):
        ms.read(out, C)  # no second read into the same folder


def test_the_read_is_refused_on_an_unclean_tree(ms, tmp_path, monkeypatch):
    _tries(tmp_path, monkeypatch, ms)
    monkeypatch.setattr(ms.gv, "tree_state", lambda: {"commit": "e" * 40, "clean": False})
    with pytest.raises(SystemExit, match="not clean"):
        ms.read(tmp_path / "out", C)


# ---------------------------------------------------------------------------------------------------------------------
# The launcher, traced under a temporary HOME


def _env(tmp_path):
    env = {k: v for k, v in os.environ.items() if not k.startswith("GIT_")}
    env.update(HOME=str(tmp_path), GIT_CEILING_DIRECTORIES=str(tmp_path.parent))
    return env


def _launch(tmp_path, *args, cwd=None, trace=False):
    cmd = ["bash", *(["-x"] if trace else []), str(LAUNCHER), *args]
    return subprocess.run(cmd, cwd=cwd or tmp_path, env=_env(tmp_path), capture_output=True, text=True, timeout=60)


def _git(repo, *args):
    return subprocess.run(["git", "-c", "user.name=t", "-c", "user.email=t@t", "-c", "commit.gpgsign=false", *args],
                          cwd=repo, env=_env(repo.parent), capture_output=True, text=True, timeout=60, check=True)


def test_the_launcher_refuses_out_of_order_and_writes_an_exit_file(tmp_path):
    S = tmp_path / "v" / "scratch" / "t65"
    r = _launch(tmp_path)
    assert r.returncode == 2 and "usage" in r.stderr and (S / "usage.exit").read_text().strip() == "2"
    r = _launch(tmp_path, "tests")
    assert r.returncode == 2
    r = _launch(tmp_path, "bogus", C)
    assert r.returncode == 97 and (S / "bogus.exit").read_text().strip() == "97"
    order = ("tests", "cut_chamber", "cut_whole", "lib_1", "lib_2", "lib_3", "read")
    for prev, step in zip(order, order[1:], strict=False):
        r = _launch(tmp_path, step, C)
        assert r.returncode == 96 and f"step {prev} has not exited 0" in r.stderr
        assert (S / f"{step}.exit").read_text().strip() == "96"
    (S / "tests.exit").write_text("0\n")
    (S / "tests.head").write_text("e" * 40 + "\n")  # passed, but at another commit
    r = _launch(tmp_path, "cut_chamber", C)
    assert r.returncode == 96 and (S / "cut_chamber.exit").read_text().strip() == "96"
    r = _launch(tmp_path, "tests", C)  # tests has no step before it; here not a git checkout
    assert r.returncode == 98 and (S / "tests.exit").read_text().strip() == "98"
    trace = _launch(tmp_path, "read", C, trace=True)  # the dry-run trace of a harmless path: it stops before Python
    assert trace.returncode == 96 and "finish 96" in trace.stderr and "+ uv run" not in trace.stderr


def test_the_launcher_refuses_an_unclean_tree_another_head_and_a_missing_row(tmp_path):
    S = tmp_path / "v" / "scratch" / "t65"
    repo = tmp_path / "repo"
    repo.mkdir()
    _git(repo, "init", "-q")
    _git(repo, "commit", "-q", "--allow-empty", "-m", "x")
    head = _git(repo, "rev-parse", "HEAD").stdout.strip()
    (repo / "stray.txt").write_text("x")
    r = _launch(tmp_path, "tests", head, cwd=repo)
    assert r.returncode == 98 and "not clean" in r.stderr and (S / "tests.exit").read_text().strip() == "98"
    _git(repo, "add", "stray.txt")
    _git(repo, "commit", "-q", "-m", "y")
    r = _launch(tmp_path, "tests", head, cwd=repo)
    assert r.returncode == 99 and (S / "tests.exit").read_text().strip() == "99"
    head = _git(repo, "rev-parse", "HEAD").stdout.strip()
    r = _launch(tmp_path, "tests", head, cwd=repo, trace=True)
    assert r.returncode == 94 and "+ uv run" not in r.stderr and (S / "tests.exit").read_text().strip() == "94"
    (repo / "experiments" / "harness").mkdir(parents=True)
    (repo / "experiments" / "harness" / "TRIES.md").write_text("| T65 (one row) | T1 | x |\n")
    _git(repo, "add", "experiments")
    _git(repo, "commit", "-q", "-m", "row")
    head = _git(repo, "rev-parse", "HEAD").stdout.strip()
    r = _launch(tmp_path, "tests", head, cwd=repo, trace=True)
    assert r.returncode == 93 and "+ uv run" not in r.stderr  # the relabelling's output is not in place
    (S / "lib_3.exit").write_text("0\n")
    (S / "lib_3.head").write_text(head + "\n")
    r = _launch(tmp_path, "read", "f" * 40, cwd=repo)
    assert r.returncode == 96  # lib_3 passed at another commit
    (S / "lib_3.head").write_text("f" * 40 + "\n")
    r = _launch(tmp_path, "read", "f" * 40, cwd=repo)
    assert r.returncode == 95 and (S / "read.exit").read_text().strip() == "95"  # not a commit of this checkout


def test_the_launcher_refuses_to_overwrite_and_needs_committed_task_files(tmp_path):
    S = tmp_path / "v" / "scratch" / "t65"
    repo = tmp_path / "repo"
    (repo / "experiments" / "harness").mkdir(parents=True)
    (repo / "experiments" / "harness" / "TRIES.md").write_text("| T65 (one row) | T1 | x |\n")
    (repo / "data" / "derived" / "feng_relabel_c").mkdir(parents=True)
    (repo / "data" / "derived" / "feng_relabel_c" / "cd1_stats_relabel.npz").write_text("x")
    (repo / ".gitignore").write_text("data/\n")
    tasks = repo / "experiments" / "t1_m1" / "results" / "t65_ms" / "tasks"
    tasks.mkdir(parents=True)
    _git(repo, "init", "-q")
    _git(repo, "add", ".")
    _git(repo, "commit", "-q", "-m", "row")
    rowc = _git(repo, "rev-parse", "HEAD").stdout.strip()
    (S).mkdir(parents=True, exist_ok=True)
    (S / "tests.exit").write_text("0\n")
    (S / "tests.head").write_text(rowc + "\n")
    (tasks / "T1_E9.5_E10.5_chamber.json").write_text("{}")  # an untracked task file: the tree still counts clean
    r = _launch(tmp_path, "cut_chamber", rowc, cwd=repo, trace=True)
    assert r.returncode == 91 and "+ uv run" not in r.stderr and (S / "cut_chamber.exit").read_text().strip() == "91"
    assert (tasks / "T1_E9.5_E10.5_chamber.json").read_text() == "{}"  # never overwritten
    (S / "lib_3.exit").write_text("0\n")
    (S / "lib_3.head").write_text(rowc + "\n")
    r = _launch(tmp_path, "read", rowc, cwd=repo, trace=True)
    assert r.returncode == 90 and "+ uv run" not in r.stderr  # the task files are not committed
    for f in ("T1_E9.5_E10.5.json", "lib1.json", "lib2.json", "lib3.json"):
        (tasks / f).write_text("{}")
    _git(repo, "add", ".")
    _git(repo, "commit", "-q", "-m", "tasks")
    (repo / "experiments" / "t1_m1" / "results" / "t65_ms" / "t65.txt").write_text("a reading")
    _git(repo, "add", ".")
    _git(repo, "commit", "-q", "-m", "reading")
    r = _launch(tmp_path, "read", rowc, cwd=repo, trace=True)
    assert r.returncode == 91 and "+ uv run" not in r.stderr and (S / "read.exit").read_text().strip() == "91"


def test_the_launcher_text(tmp_path):
    sh = LAUNCHER.read_text()
    assert ":?" not in sh and "rm " not in sh and "rm -" not in sh  # no unsafe message, no removal
    assert subprocess.run(["bash", "-n", str(LAUNCHER)], capture_output=True).returncode == 0
    tests = sh.split("  tests)\n")[1].split(";;")[0]
    for f in ("test_t65_ms_read.py", "test_t64_heldmom_read.py", "test_t62_recomputed_read.py",
              "test_t63_double_step_read.py", "test_t1_bets.py", "test_t1_scope.py", "test_core.py",
              "test_g_mom_verify.py", "test_t1_state_transport.py"):
        assert f"tests/{f}" in tests and (ROOT / "tests" / f).exists()
    assert "ms_read.py" in sh and "NUMPY_MADVISE_HUGEPAGE=0" in sh and "OUT=t65_ms" in sh
    assert "feng_relabel_c/cd1_stats_relabel.npz" in sh and '--commit "$ROWC"' in sh
    runs = [ln for ln in sh.splitlines() if "ms_read.py" in ln and "$PY" in ln] + [tests]
    assert len(runs) == 3 and all('> "$S/$STEP.out" 2> "$S/$STEP.err"' in r for r in runs)  # no log is printed
