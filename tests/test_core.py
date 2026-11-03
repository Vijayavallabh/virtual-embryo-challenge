"""vec.core, the shared transport-then-adjust core (feat-020, docs/PLAN.md section 5), on synthetic data only.

1. Golden outputs. The builders that go through the core (vec.t1_scope: reweight, build_submission with and without
   a momentum shift; vec.t2_build: sample_cells, cloud_guard, build_submission) and the moved primitives themselves,
   on small synthetic inputs with fixed seeds, called by their old names. The hashes were recorded at eb946ac, before
   the primitives moved into vec.core, and must not change: the move is a pure refactor. Arrays (X, obs, var,
   obsm, layers, rows, coordinates, counts) are hashed bit for bit; report dicts with their floats at 9 significant
   digits, because report floats that pass through log or exp (the momentum rows) can differ in the last bit between
   CPU code paths (numpy's AVX2 loops here, AVX-512 on the DGX). An output change that is meant must update the table
   in the same commit and say why (`uv run python tests/test_core.py` prints it).
2. The pre-registered T3 script (experiments/t3_decomp/t3_decompose.py, hash-pinned by its build_submission.py) keeps
   frozen copies of the core's T3 primitives: the same bodies, and the same outputs and generator states as the core's
   for the same inputs and seeds, edge cases included.
3. The old names in vec.t1_scope and vec.t2_build are the core's objects, the heads call the core and define none of
   its functions, and the core imports neither head.
4. No size, composition or class list of any stage is written into the core.

vec.core is imported by a fixture, not at the top, so that the golden tests also run against the code from before
the move (`pytest tests/test_core.py -k golden` with that code on the path).
"""

from __future__ import annotations

import ast
import copy
import hashlib
import importlib
import importlib.util
import inspect
import io
import json
import re
import subprocess
import sys
import textwrap
import tokenize
from dataclasses import asdict
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import yaml
from scipy import sparse

from vec import t1_scope as ts
from vec import t2_build as tb
from vec.boards import Board
from vec.prep import log1p_cp10k

REPO = Path(__file__).resolve().parents[1]
REPORT_DIGITS = ".9g"  # report floats (see the module docstring); arrays are hashed bit for bit

# ---------------------------------------------------------------------------------------------------------------------
# Hashing


def _sha(*parts) -> str:
    h = hashlib.sha256()
    for p in parts:
        h.update(p if isinstance(p, bytes) else str(p).encode())
        h.update(b"\x1e")
    return h.hexdigest()


def h_array(a) -> str:
    """dtype, shape and bytes of a numeric array; a CSR/CSC matrix by its format and its stored arrays as they are."""
    if sparse.issparse(a):
        assert a.format in ("csr", "csc"), a.format
        return _sha(a.format, a.dtype.str, a.shape, h_array(a.data), h_array(a.indices), h_array(a.indptr))
    a = np.asarray(a)
    assert a.dtype != object
    return _sha(a.dtype.str, a.shape, np.ascontiguousarray(a).tobytes())


def h_frame(df: pd.DataFrame) -> str:
    """Index, column names, dtypes and every value (repr) of a data frame."""
    cols = [[str(c), repr(df[c].dtype), [repr(v) for v in df[c].tolist()]] for c in df.columns]
    return _sha(json.dumps([[repr(i) for i in df.index], repr(df.index.dtype), cols]))


def h_adata(a: ad.AnnData) -> dict[str, str]:
    layers = {k: h_array(a.layers[k]) for k in sorted(k for k in a.layers.keys() if k is not None)}
    return {"X": h_array(a.X), "obs": h_frame(a.obs), "var": h_frame(a.var),
            "obsm": _sha(json.dumps({k: h_array(a.obsm[k]) for k in sorted(a.obsm.keys())})),
            "layers": _sha(json.dumps(layers))}


def canon(o):
    """A report as JSON-able values: numpy scalars and arrays as Python values, floats at REPORT_DIGITS."""
    if isinstance(o, dict):
        return {str(k): canon(v) for k, v in o.items()}
    if isinstance(o, list | tuple):
        return [canon(v) for v in o]
    if isinstance(o, np.ndarray):
        return canon(o.tolist())
    if isinstance(o, bool | np.bool_):
        return bool(o)
    if isinstance(o, int | np.integer):
        return int(o)
    if isinstance(o, float | np.floating):
        return format(float(o), REPORT_DIGITS)
    if o is None or isinstance(o, str):
        return o
    raise TypeError(f"no canonical form for {type(o).__name__}")


def h_record(o) -> str:
    return _sha(json.dumps(canon(o), sort_keys=True))


def prefixed(prefix: str, d: dict[str, str]) -> dict[str, str]:
    return {f"{prefix}_{k}": v for k, v in d.items()}


# ---------------------------------------------------------------------------------------------------------------------
# Synthetic inputs, built as tests/test_t1_scope.py and tests/test_t2_build.py build theirs. Every number is made up:
# no size, share or count of any stage.

T1_GENES = [f"g{i}" for i in range(58)] + ["Xist", "Ddx3y"]
T1_HELD = ["Xist", "Ddx3y"]
T1_SIZES = {"V-CM": 300, "Endocardium": 90, "Proepicardium": 40, "Blood": 25, "aSHF": 70, "JCF": 50, "Foregut": 120,
            "Neural Tube": 60}
T1_NO_EPI = {k: v for k, v in T1_SIZES.items() if k != "Proepicardium"}  # an E8.5-like file: JCF stands in
T1_SHARES = {"cm": 0.4, "endocardium": 0.2, "epicardium": 0.1, "blood": 0.08, "progenitor": 0.12, "surrounding": 0.1}
# The chamber config is the test's own, so that a change to configs/t1_chamber_classes.yaml cannot move these hashes.
T1_CONFIG = """\
categories:
  cm: {scope: chamber}
  endocardium: {scope: chamber}
  epicardium: {scope: chamber}
  blood: {scope: chamber}
  progenitor: {scope: adjacent}
  surrounding: {scope: surrounding}
heart: [cm, endocardium, epicardium]
labels: {V-CM: cm, Endocardium: endocardium, Proepicardium: epicardium, Blood: blood, aSHF: progenitor,
         JCF: progenitor, Foregut: surrounding, Neural Tube: surrounding}
file_fallbacks:
  E8.5_RNA: {epicardium: [JCF]}
"""
T2_GENES = [f"g{i}" for i in range(60)]


def chamber_spec(directory: Path) -> ts.ChamberSpec:
    path = directory / "chamber.yaml"
    path.write_text(T1_CONFIG)
    return ts.load_spec(path)


@pytest.fixture(scope="module")
def spec(tmp_path_factory):
    return chamber_spec(tmp_path_factory.mktemp("chamber"))


def t1_profiles(seed=0):
    rng = np.random.default_rng(seed)
    out = {}
    for lab in T1_SIZES:
        mu = rng.gamma(0.4, 1.0, size=len(T1_GENES))  # many low genes, so cells have exact zeros
        mu[T1_GENES.index("Xist")] = 3.0
        mu = 2000.0 * mu / mu.sum()
        mu[T1_GENES.index("g1")] = 20.0
        out[lab] = mu
    return out


def t1_stage(seed, sizes=T1_SIZES, boost=None):
    """Poisson counts per label with cell-depth variation; boost = (label, gene, factor) scales one gene's mean."""
    rng = np.random.default_rng(seed)
    profiles = t1_profiles()
    blocks, labels = [], []
    for lab, n in sizes.items():
        mu = profiles[lab].copy()
        if boost and boost[0] == lab:
            mu[T1_GENES.index(boost[1])] *= boost[2]
        C = rng.poisson(mu * rng.gamma(10.0, 0.1, size=(n, 1)))
        C[np.arange(n), rng.integers(0, len(T1_GENES), n)] += 1  # every cell has counts
        blocks.append(C)
        labels += [lab] * n
    C = sparse.csr_matrix(np.vstack(blocks).astype(np.int32))
    a = ad.AnnData(X=log1p_cp10k(C))
    a.var_names = T1_GENES
    a.obs_names = [f"s{seed}_{i}" for i in range(a.n_obs)]
    a.obs["celltype"] = labels
    a.layers["counts"] = C
    return a


def t1_board():
    return Board(key="T1:synthetic", task="T1", genes=T1_GENES[::-1], needs_coords=False, min_cells=100,
                 max_cells=2000)  # reversed gene order: the builder must reorder


def t2_stack(n, seed=0, extent=(120.0, 80.0, 60.0), step=15.0):
    """A uniformly filled box cut into z sections, rows stored section by section in descending z."""
    rng = np.random.default_rng(seed)
    C = rng.uniform(-1.0, 1.0, size=(n, 3)) * np.asarray(extent)
    C[:, 2] = np.round(C[:, 2] / step) * step
    return C[np.argsort(-C[:, 2], kind="stable")]


def t2_cells(n, seed=0):
    rng = np.random.default_rng(seed)
    counts = rng.poisson(rng.gamma(2.0, 1.0, size=(1, len(T2_GENES))), size=(n, len(T2_GENES)))
    lib = counts.sum(1, keepdims=True).clip(min=1)
    return np.log1p(1e4 * counts / lib).astype(np.float32)


# ---------------------------------------------------------------------------------------------------------------------
# Golden cases: each returns part -> sha256


def case_core_quotas(spec):
    cases = [(T1_SHARES, 0), (T1_SHARES, 7), (T1_SHARES, 1234), ({"b": 0.25, "a": 0.25, "d": 0.25, "c": 0.25}, 2),
             ({"x": 1 / 3, "y": 1 / 3, "z": 1 / 3}, 10), ({"x": 0.5, "y": 0.5, "z": 0.0}, 3), ({"p": 0.7, "q": 0.3}, 1)]
    return {"quotas": h_record([list(ts.quotas(s, n).items()) for s, n in cases])}


def case_core_draw_least_used(spec):
    rng = np.random.default_rng(11)
    used = np.zeros(40, dtype=np.int64)
    out = []
    for pool, k in [(np.arange(10), 25), (np.arange(5, 20), 12), (np.arange(0, 40, 3), 0), (np.arange(30, 40), 7),
                    (np.array([3, 1, 4, 15, 9, 2, 6]), 30), (np.arange(40), 41)]:
        picks, reused = ts._draw(pool, k, used, rng)
        out.append([h_array(picks), reused, h_array(used)])
    return {"draws": h_record(out), "rng_state": h_record(rng.bit_generator.state)}


def case_core_draw_balanced(spec):
    out = []
    for n_pool, n, seed in [(10, 3, 0), (10, 10, 1), (10, 25, 2), (7, 0, 3), (1, 5, 4), (100, 333, [5, 1])]:
        rng = np.random.default_rng(seed)
        out.append([h_array(tb._balanced_rows(n_pool, n, rng)), rng.bit_generator.state])
    return {"rows": h_record(out)}


def case_core_counts_and_means(spec):
    rng = np.random.default_rng(12)
    C = rng.poisson(rng.gamma(0.5, 2.0, size=(1, 30)), size=(80, 30)).astype(np.int32)
    C[np.arange(80), rng.integers(0, 30, 80)] += 1  # every cell has counts (recovery from X needs it)
    genes = [f"h{i}" for i in range(30)]
    a = ad.AnnData(X=log1p_cp10k(C))
    a.var_names = genes
    from_x = ts.counts_of(a)
    a.layers["counts"] = sparse.csr_matrix(C.astype(np.float32))  # whole numbers stored as floats
    from_layer = ts.counts_of(a, genes[::-1])
    E = C.copy()
    E[5] = 0  # an empty cell
    return {"counts_from_x": h_array(from_x), "counts_from_layer": h_array(from_layer),
            "mean_cp10k_sparse": h_array(ts.mean_cp10k(sparse.csr_matrix(E))),
            "mean_cp10k_dense": h_array(ts.mean_cp10k(E))}


def case_core_apply_factors(spec):
    rng = np.random.default_rng(13)
    n, g = 300, 25
    D = rng.poisson(rng.gamma(0.6, 3.0, size=(1, g)), size=(n, g))
    D[0] = 0  # an empty cell
    C = sparse.csr_matrix(D)
    F = rng.lognormal(0.0, 0.6, size=(3, g))
    F[0, :4] = 1.0
    F[2] = 0.02  # type 2 nearly erased: without the depth correction some cells lose every count and are restored
    tix = rng.integers(0, 3, n)
    held = np.zeros(g, dtype=bool)
    held[[3, 7]] = True
    out = []
    for keep_depth, seed in [(True, 0), (False, 1), (True, [2, 5])]:
        new, scale, dead = ts.apply_factors(C, F, tix, held, keep_depth=keep_depth, seed=seed)
        out.append([h_array(new), h_array(scale), h_array(dead), int(dead.sum())])
    return {"apply_factors": h_record(out)}


def case_t1_reweight(spec):
    a = t1_stage(1)
    ts.classify(a, spec=spec)
    out_a, rep_a = ts.reweight(a, T1_SHARES, 700, seed=3, within={"surrounding": ["Foregut"]}, spec=spec)
    b = t1_stage(2, sizes=T1_NO_EPI)  # epicardium from JCF stand-ins, a pool progenitor also draws from
    ts.classify(b, spec=spec)
    out_b, rep_b = ts.reweight(b, T1_SHARES, 1234, seed=[7, 2], fallback={"epicardium": ["JCF"]}, spec=spec)
    return {**prefixed("a", h_adata(out_a)), "a_report": h_record(rep_a.as_dict()),
            **prefixed("b", h_adata(out_b)), "b_report": h_record(rep_b.as_dict())}


def _t1_build(ref, n_cells, seed, spec, **kw):
    sub, rep = ts.build_submission(ref, t1_board(), T1_SHARES, n_cells, seed, spec=spec, **kw)
    return {**h_adata(sub), "report": h_record({**rep, "config": Path(rep["config"]).name})}


def case_t1_build_no_shift(spec):
    return _t1_build(t1_stage(22), 300, 4, spec, ref_source="")


def case_t1_build_shift(spec):
    start, ref = t1_stage(20), t1_stage(21, boost=("V-CM", "g1", 4.0))
    return _t1_build(ref, 300, 0, spec, start=start, damping=0.5, held_flat=T1_HELD, ref_source="", start_source="")


def case_t1_build_shift_options(spec):
    """Stand-ins in the draw and at the reference end of the shift (an E8.5-like reference), restricted labels,
    momentum for stand-in types, no depth correction, clipped factors."""
    start = t1_stage(23, boost=("Endocardium", "g2", 0.25))
    ref = t1_stage(24, sizes=T1_NO_EPI, boost=("V-CM", "g1", 3.0))
    return _t1_build(ref, 450, 5, spec, start=start, damping=1.0, held_flat=T1_HELD, ref_source="E8.5_RNA",
                     start_source="", within={"surrounding": ["Foregut"]}, keep_depth=False, max_abs_log_factor=0.8,
                     standin_momentum=True, eps=0.25, min_cells=15)


def case_t2_sample_cells(spec):
    pool = t2_stack(2000, seed=1)
    dup = t2_stack(1000, seed=2)
    dup[11] = dup[10]
    dup[500] = dup[10]
    dup[700] = dup[20] + np.array([0.0, 0.0, 1e-12])  # one position as the scorer reads it (float32)
    a = ad.AnnData(X=np.zeros((800, len(T2_GENES)), np.float32))
    a.obsm["spatial_3D"] = t2_stack(800, seed=3).astype(np.float32)
    out = {}
    for name, (stage, n, seed, frac) in {"below_pool": (pool, 500, 3, 0.25), "above_pool": (pool, 5300, 0, 0.25),
                                         "pool_duplicates": (dup, 1000, 0, 0.1), "anndata": (a, 700, 9, 0.2)}.items():
        s = tb.sample_cells(stage, n, seed=seed, jitter_frac=frac)
        out.update({f"{name}_rows": h_array(s.rows), f"{name}_coords": h_array(s.coords),
                    f"{name}_record": h_record(s.record)})
    return out


def case_t2_cloud_guard(spec):
    pool = t2_stack(6000, seed=11)
    sample = tb.sample_cells(pool, 1500, seed=1).coords
    over = tb.sample_cells(pool, 8000, seed=2).coords  # above the pool: the reference draws repeat rows evenly
    return {"random": h_record(tb.cloud_guard(sample, pool, seed=0).as_dict()),
            "block_own_frame": h_record(tb.cloud_guard(pool[:1500], pool, seed=3, frame="own", n_ref=3).as_dict()),
            "above_pool": h_record(tb.cloud_guard(over, pool, seed=4).as_dict())}


def case_t2_build_submission(spec):
    board = Board(key="T2:synthetic", task="T2", genes=T2_GENES, needs_coords=True, min_cells=10, max_cells=6000,
                  truth_cells=400)
    pool = t2_stack(2500, seed=31)
    s = tb.sample_cells(pool, 3000, seed=6)
    X = t2_cells(2500, seed=32)
    order = np.random.default_rng(33).permutation(len(T2_GENES))
    cells = ad.AnnData(X=np.hstack([X[:, order], np.ones((2500, 2), np.float32)])[s.rows])
    cells.var_names = [T2_GENES[i] for i in order] + ["Extra1", "Extra2"]
    sub_a, rep_a = tb.build_submission(board, cells, s.coords)
    plain = Board(key="T2:synthetic", task="T2", genes=T2_GENES, needs_coords=True, min_cells=10, max_cells=6000)
    coords = pool[:300].copy()
    coords[1] = coords[0]
    sub_b, rep_b = tb.build_submission(plain, X[:300], coords, genes=T2_GENES, allow_duplicate_coords=True)
    return {**prefixed("a", h_adata(sub_a)), "a_check": h_record(asdict(rep_a)),
            **prefixed("b", h_adata(sub_b)), "b_check": h_record(asdict(rep_b))}


CASES = {name.removeprefix("case_"): fn for name, fn in globals().items() if name.startswith("case_")}

# Recorded at eb946ac (before the move) with `uv run python tests/test_core.py`.
GOLDEN: dict[str, dict[str, str]] = {
    "core_quotas": {
        "quotas": "e57116ccf364b52da783c518012566bf40d413951024a285748dae982a8c86c3",
    },
    "core_draw_least_used": {
        "draws": "03eff7afbf42f402749bfc367250c4e48e780eff3563c1cbe22f148b085b1a0d",
        "rng_state": "1acb607916dc2889ae8a224907c7838b5d1769666e6d278c48d7c0e36dd2821f",
    },
    "core_draw_balanced": {
        "rows": "1175448dd8dd3df1d21578a737e137d137dee5ba37ca1fd1b7f13844df28d1ca",
    },
    "core_counts_and_means": {
        "counts_from_x": "a7a0b39ca9bc6a5120a601cfee7dd159e824466f4b7bba37801f670dab04f4e0",
        "counts_from_layer": "19020008b18b021dc5f481fe6ffd3f3652286193498cbd6f3e9c088ddcd77ced",
        "mean_cp10k_sparse": "f8bee5fcfe2fe2814e8f5cb67622db8fa7a1fdd1a6cb3315ed7c6ec0fb4025ab",
        "mean_cp10k_dense": "f8bee5fcfe2fe2814e8f5cb67622db8fa7a1fdd1a6cb3315ed7c6ec0fb4025ab",
    },
    "core_apply_factors": {
        "apply_factors": "562d4508973bf22aec27a19db2644f74a1d3824316baf8ba83abc9c4f4c7330d",
    },
    "t1_reweight": {
        "a_X": "543eaacfef6df469fc7d3d2011d0ed9ec79743a5f094f61fd182af9e80b15c00",
        "a_obs": "b34ef217e78492c68133e7696b6a84f7e056b9ec67e47a34c51d3d6c29165083",
        "a_var": "b8413d19e77df5383d58f6ccba66a6623cf93d9ea9a2c131a3adfbbb146e6c14",
        "a_obsm": "f0594ba16cd5803b42e7f66606fab0db7f1efcfc305afd5cc09b53086dfe8c6a",
        "a_layers": "629c9546e81f55445271c27f728df70dcfacf23fe563bb89de36ebd31d79c56e",
        "a_report": "0f5cde0c53fadaafd020778c1fba0a688cef2f9d382161aa5691eefcc20f3f45",
        "b_X": "7659a0b69de16e84479030631f0294bc5d98e7f305994ea68555bcd537d22215",
        "b_obs": "57ad317a1d992c3ef6f5932de2973af9542859f164b6f3d783f9b940ea3e6ea0",
        "b_var": "b8413d19e77df5383d58f6ccba66a6623cf93d9ea9a2c131a3adfbbb146e6c14",
        "b_obsm": "f0594ba16cd5803b42e7f66606fab0db7f1efcfc305afd5cc09b53086dfe8c6a",
        "b_layers": "15e0a3625ad37249770bc61260a4fbd38275c34d53b72cbc7b8abcd6a4af4f0e",
        "b_report": "527844a596df0424ebad0789634f1f4af407d0a7236f0de8a80128b9785d880f",
    },
    "t1_build_no_shift": {
        "X": "0812a292b5befc854df5e4e4693e3844d500e63b65029ab30b66ed3e5f6e4499",
        "obs": "bb4dd322eb56536c2a496bd86105bca0e846dfb4ecc7abde8f2ac7d44cd3ee47",
        "var": "1a415014433a4b3a9467e1e19058f267f83d8b77f5ab1e6fd06f64123e18acc4",
        "obsm": "f0594ba16cd5803b42e7f66606fab0db7f1efcfc305afd5cc09b53086dfe8c6a",
        "layers": "f0594ba16cd5803b42e7f66606fab0db7f1efcfc305afd5cc09b53086dfe8c6a",
        "report": "66e2918ee4abb340c3ca1e6233bbb456177cd4a45783cd754e52a9d138778669",
    },
    "t1_build_shift": {
        "X": "5b330590772c05776d3c03a5c134c494ec3fe37bd1a9419f09d255b15eedfa20",
        "obs": "bb4dd322eb56536c2a496bd86105bca0e846dfb4ecc7abde8f2ac7d44cd3ee47",
        "var": "1a415014433a4b3a9467e1e19058f267f83d8b77f5ab1e6fd06f64123e18acc4",
        "obsm": "f0594ba16cd5803b42e7f66606fab0db7f1efcfc305afd5cc09b53086dfe8c6a",
        "layers": "f0594ba16cd5803b42e7f66606fab0db7f1efcfc305afd5cc09b53086dfe8c6a",
        "report": "a5d5c35c3575ac35e540605adcafb6b58cffdba3f90062dbe067a3341100c5d5",
    },
    "t1_build_shift_options": {
        "X": "317a33677ea4febee9e3a19f9cdbcb1ae9e52096d82c32e946cb00f96d079fb1",
        "obs": "a8722d49946768121c3d39d756d502d02c74e5caa550d474cd38c0cb02359d7e",
        "var": "1a415014433a4b3a9467e1e19058f267f83d8b77f5ab1e6fd06f64123e18acc4",
        "obsm": "f0594ba16cd5803b42e7f66606fab0db7f1efcfc305afd5cc09b53086dfe8c6a",
        "layers": "f0594ba16cd5803b42e7f66606fab0db7f1efcfc305afd5cc09b53086dfe8c6a",
        "report": "f20a8714a853c8a77003eb7a2a9688f62f4cf4d86fc45ed8108c29f24db7bef1",
    },
    "t2_sample_cells": {
        "below_pool_rows": "3490a2e245c4b08a52c565d05f3469db71eea97bef1eb6c7da40cabb017a0d5b",
        "below_pool_coords": "df96f82d632eec0ab976b52be5518be2d313031c8065312a9c8ae0b4ddeb9d96",
        "below_pool_record": "396daf38cbee7edee310a3e674c696fb4270d4f1ae1c2b3b2d4bf0f7cf5d380b",
        "above_pool_rows": "9176d69b58da38ed3fb0b4bd0714cbb81b0e959e18a2c407e7f91e03cdfc7e02",
        "above_pool_coords": "cd88726462307b7f3b861b47e2a387f2b4a9a66b5271edffc542a90eaae33c23",
        "above_pool_record": "77a4a7fd99906a6cd88feb9c7302349bf144fc89e15f8df2a55e60908b0d4b90",
        "pool_duplicates_rows": "89f10f127864cd50c47f899cfa3540e2cc99382887a85c8251b87a4ad4eb3c3f",
        "pool_duplicates_coords": "4c057c763ce48f6f4286cfb1c6413860c24c4b2b0a124e21fc76f52d86816e4c",
        "pool_duplicates_record": "d6c1678a384566327c3111845b0a1d5db661d6dbc568c9640aa32df73dfb2135",
        "anndata_rows": "1766c384946cb85ee51db9dcd3e75be958c0e44bfcefeebf00df5470b9b3ae45",
        "anndata_coords": "cd7c72eb3b4143022ffe87233844b7b3c09d2b986d855d32d12f9c79028c40e7",
        "anndata_record": "f82bb29292ca240dbd96957e95c7f031b0aa1f56424da73831fcea08287f1113",
    },
    "t2_cloud_guard": {
        "random": "82938f2b26986adf749dc837ab40d31bcac516d568fa5bb51dbf447dca14a512",
        "block_own_frame": "8738697de442f0c43c53d76efb75d3c46bb8bd98f5b75f8c76bd29d2a9a76641",
        "above_pool": "4479b7ab1462a198f748a6e0451e8e6793bb15b0bb4863f62f94db4eceab048c",
    },
    "t2_build_submission": {
        "a_X": "f9511079625d8a3a9ab93aec09f3cb49e3f7320a182a75f0ed80187e578ae1a4",
        "a_obs": "2d443dc6f43519cd7de4403d10d494fcba34db67dd91d4554240cd0ddaba3c07",
        "a_var": "148073843defa287b4006776c72ced25be2cf79919a571738eb963e9cdd60dcc",
        "a_obsm": "46ab414c89e804dd2a152b0c70331ab232aa51c8df2a8adf6ecb3544848dd722",
        "a_layers": "f0594ba16cd5803b42e7f66606fab0db7f1efcfc305afd5cc09b53086dfe8c6a",
        "a_check": "6e2d06bbe57935488dc276230edd3a58a846573d0abdc5fc49987558e78bb4b0",
        "b_X": "71932b1db43fbdf9af1fd532f645fc146847e6acb528c32805a29f37e46f5ef5",
        "b_obs": "bb4dd322eb56536c2a496bd86105bca0e846dfb4ecc7abde8f2ac7d44cd3ee47",
        "b_var": "148073843defa287b4006776c72ced25be2cf79919a571738eb963e9cdd60dcc",
        "b_obsm": "f27e8306f4832e4a4634eae20049f0597ede1c131cadc939fef6d48224391702",
        "b_layers": "f0594ba16cd5803b42e7f66606fab0db7f1efcfc305afd5cc09b53086dfe8c6a",
        "b_check": "7a72e9096e941a6c5a59e77c3e9932e2a024a84c29c6e02744793ce0cf94fbaf",
    },
}


def test_every_golden_case_has_recorded_hashes():
    assert set(CASES) == set(GOLDEN)


@pytest.mark.parametrize("name", sorted(CASES))
def test_golden_outputs_unchanged(name, spec):
    got = CASES[name](spec)
    want = GOLDEN[name]
    assert got == want, {k: (got.get(k), want.get(k)) for k in sorted(set(got) | set(want))
                         if got.get(k) != want.get(k)}


# ---------------------------------------------------------------------------------------------------------------------
# 2. The frozen T3 copies

T3_SCRIPT = REPO / "experiments" / "t3_decomp" / "t3_decompose.py"
T3_PAIRS = {"alloc": "alloc", "draw_rows": "draw_then_repeat", "thin": "thin", "rescale_depth": "rescale_depth",
            "apply_gene_factor": "apply_gene_factor"}  # name in the pre-registered script -> name in the core
N_RANDOM = 300  # random inputs per primitive


@pytest.fixture(scope="module")
def op():
    return importlib.import_module("vec.core")


@pytest.fixture(scope="module")
def t3d():
    """experiments/t3_decomp/t3_decompose.py loaded by path, as experiments/t3_decomp/test_t3_decompose.py loads it
    (under a module name of its own, so that the two tests never share or replace each other's module object)."""
    spec_ = importlib.util.spec_from_file_location("t3_decompose_for_test_core", T3_SCRIPT)
    mod = importlib.util.module_from_spec(spec_)
    sys.modules[spec_.name] = mod
    spec_.loader.exec_module(mod)
    return mod


def body_source(fn) -> str:
    """A function's source without its def line and docstring."""
    src = textwrap.dedent(inspect.getsource(fn))
    body = ast.parse(src).body[0].body
    if isinstance(body[0], ast.Expr) and isinstance(body[0].value, ast.Constant) and isinstance(body[0].value.value,
                                                                                                 str):
        body = body[1:]
    return "\n".join(src.splitlines()[body[0].lineno - 1:body[-1].end_lineno])


@pytest.mark.parametrize(("frozen", "core"), T3_PAIRS.items())
def test_t3_copies_have_the_frozen_bodies(t3d, op, frozen, core):
    a, b = getattr(t3d, frozen), getattr(op, core)
    assert body_source(a) == body_source(b)
    assert inspect.signature(a) == inspect.signature(b)
    # the globals the body reads are the same objects in both modules (numpy), or copies proven here (thin)
    for name in set(a.__code__.co_names) & (set(a.__globals__) | set(b.__globals__)):
        if name in T3_PAIRS:
            assert T3_PAIRS[name] == name, name
        else:
            assert a.__globals__.get(name) is b.__globals__.get(name), name


def run_pair(fa, fb, args, seed=None):
    """fa and fb on deep copies of the same arguments (and identically seeded generators): the same result, dtype and
    shape, the same arguments afterwards, and the generators left in the same state (the same draws consumed)."""
    args_a, args_b = copy.deepcopy(args), copy.deepcopy(args)
    if seed is None:
        a, b = fa(*args_a), fb(*args_b)
    else:
        ra, rb = np.random.default_rng(seed), np.random.default_rng(seed)
        a, b = fa(*args_a, ra), fb(*args_b, rb)
        assert ra.bit_generator.state == rb.bit_generator.state
    assert type(a) is type(b) and a.dtype == b.dtype and a.shape == b.shape
    assert np.array_equal(a, b)
    for x, y in zip(args_a, args_b, strict=True):
        assert np.array_equal(np.asarray(x), np.asarray(y))
    return a


def random_counts(rng, n=None, g=None) -> np.ndarray:
    """Dense int32 counts, cells x genes, with some empty cells."""
    n = int(rng.integers(0, 40)) if n is None else n
    g = int(rng.integers(1, 15)) if g is None else g
    C = rng.poisson(rng.gamma(0.5, 2.0, size=(1, g)), size=(n, g)).astype(np.int32)
    C[rng.random(n) < 0.15] = 0
    return C


def random_props(rng, C) -> np.ndarray:
    P = rng.random(C.shape) + 1e-3
    return (P / P.sum(1, keepdims=True)).astype(np.float32)


def random_keep_prob(rng, n):
    """A keep probability as thin takes it: a scalar, per cell, or per cell as a column, outside [0, 1] at times."""
    kind = int(rng.integers(0, 3))
    if kind == 0:
        return float(rng.choice([0.0, 1.0, 0.3, rng.random()]))
    p = rng.uniform(-0.2, 1.2, n)
    p[rng.random(n) < 0.2] = 1.0
    p[rng.random(n) < 0.1] = 0.0
    return p if kind == 1 else p.reshape(-1, 1)


def random_ratio(rng, n):
    """A depth ratio per cell, exactly 1, below 1 (0 included) and above 1; or a scalar of at most 1."""
    if rng.random() < 0.1:
        return float(rng.choice([1.0, 0.5]))
    r = rng.uniform(0.0, 3.0, n)
    r[rng.random(n) < 0.25] = 1.0
    r[rng.random(n) < 0.05] = 0.0
    return r


def random_phi(rng, g):
    phi = rng.uniform(0.0, 3.0, g)
    phi[rng.random(g) < 0.3] = 1.0
    return phi


def random_shares(rng):
    m = int(rng.integers(1, 12))
    kind = int(rng.integers(0, 4))
    if kind == 0:
        return np.full(m, 1.0 / m)  # ties
    f = rng.dirichlet(np.ones(m))
    f[rng.random(m) < 0.2] = 0.0
    if kind == 1:
        return list(f / f.sum() if f.sum() > 0 else np.full(m, 1.0 / m))
    return f * rng.uniform(0.8, 1.2) if kind == 2 else f  # shares that do not sum to 1


def test_t3_copies_match_on_random_inputs(t3d, op):
    rng = np.random.default_rng(2026)
    for i in range(N_RANDOM):
        n = int(rng.integers(0, 3000))
        run_pair(t3d.alloc, op.alloc, (n, random_shares(rng)))
        m = int(rng.integers(0, 30))
        pool = np.sort(rng.choice(1000, m, replace=False))
        run_pair(t3d.draw_rows, op.draw_then_repeat, (pool if i % 4 else list(pool), int(rng.integers(-2, 80))), i)
        C = random_counts(rng)
        run_pair(t3d.thin, op.thin, (C, random_keep_prob(rng, C.shape[0])), i)
        run_pair(t3d.rescale_depth, op.rescale_depth, (C, random_ratio(rng, C.shape[0]), random_props(rng, C)), i)
        run_pair(t3d.apply_gene_factor, op.apply_gene_factor, (C, random_phi(rng, C.shape[1]), random_props(rng, C)),
                 [i, 1])


def test_t3_copies_match_on_edge_cases(t3d, op):
    rng = np.random.default_rng(7)
    C = random_counts(rng, n=30, g=8)
    C[[0, 3, 4]] = 0  # empty cells
    P = random_props(rng, C)
    empty = np.zeros((0, 8), dtype=np.int32)  # no cells at all
    for seed in range(20):
        # empty cells stay empty; no cells give no cells
        out = run_pair(t3d.thin, op.thin, (C, rng.uniform(0.0, 1.0, 30)), seed)
        assert not out[[0, 3, 4]].any() and (out.sum(1) > 0).sum() == (C.sum(1) > 0).sum()
        assert run_pair(t3d.thin, op.thin, (empty, 0.5), seed).shape == (0, 8)
        assert run_pair(t3d.rescale_depth, op.rescale_depth, (empty, np.ones(0), P[:0]), seed).shape == (0, 8)
        assert run_pair(t3d.apply_gene_factor, op.apply_gene_factor, (empty, random_phi(rng, 8), P[:0]),
                        seed).shape == (0, 8)
        # ratio = 1 and phi = 1 change nothing
        assert np.array_equal(run_pair(t3d.rescale_depth, op.rescale_depth, (C, np.ones(30), P), seed), C)
        assert np.array_equal(run_pair(t3d.rescale_depth, op.rescale_depth, (C, 1.0, P), seed), C)
        assert np.array_equal(run_pair(t3d.apply_gene_factor, op.apply_gene_factor, (C, np.ones(8), P), seed), C)
        # k > pool: the whole pool once, then draws from it with replacement; k = 0 (or below): nothing
        pool = np.array([5, 9, 2, 40])
        out = run_pair(t3d.draw_rows, op.draw_then_repeat, (pool, 11), seed)
        assert np.array_equal(out[:4], pool) and set(out[4:]) <= set(pool) and len(out) == 11
        assert np.array_equal(np.sort(run_pair(t3d.draw_rows, op.draw_then_repeat, (pool, 4), seed)), np.sort(pool))
        for k in (0, -1):
            assert run_pair(t3d.draw_rows, op.draw_then_repeat, (pool, k), seed).size == 0
        assert run_pair(t3d.draw_rows, op.draw_then_repeat, (np.array([], dtype=int), 5), seed).size == 0
    # ties in the allocation (equal remainders), zero shares, no cells, shares that do not sum to 1
    for n, f in [(10, [1 / 3] * 3), (2, [0.25] * 4), (5, [0.1] * 10), (1, [0.5, 0.5]), (7, [0.5, 0.0, 0.5]),
                 (0, [0.2, 0.8]), (100, [0.3, 0.3]), (100, [0.6, 0.6]), (3, [1.0])]:
        k = run_pair(t3d.alloc, op.alloc, (n, f))
        floors = int(np.floor(n * np.asarray(f)).sum())
        assert k.sum() == floors + min(max(n - floors, 0), len(f))  # at most one more cell per share


def test_t3_lognorm_is_the_core_log1p_cp10k_within_one_float32_ulp(t3d, op):
    """Not a copy (log1p of 1e4 * c / n in the script, of c * (1e4 / n) in vec.prep), so equal only to float32
    rounding; the core's log1p CP10k step is vec.prep's own function."""
    from vec import prep

    assert op.log1p_cp10k is prep.log1p_cp10k
    rng = np.random.default_rng(3)
    for _ in range(100):
        C = random_counts(rng)
        np.testing.assert_array_max_ulp(t3d.lognorm(C), op.log1p_cp10k(C).toarray(), maxulp=1)


# ---------------------------------------------------------------------------------------------------------------------
# 3. Old names, calls and imports

MOVED = {("vec.t1_scope", "quotas"): "quotas", ("vec.t1_scope", "_draw"): "draw_least_used",
         ("vec.t1_scope", "counts_of"): "counts_of", ("vec.t1_scope", "mean_cp10k"): "mean_cp10k",
         ("vec.t1_scope", "apply_factors"): "apply_factors", ("vec.t2_build", "_balanced_rows"): "draw_balanced"}
# head function -> the core functions it must call
CALLS = {("t1_scope", "reweight"): {"quotas", "draw_least_used"},
         ("t1_scope", "momentum_shift"): {"counts_of", "mean_cp10k", "apply_factors", "log1p_cp10k"},
         ("t2_build", "sample_cells"): {"draw_balanced"}, ("t2_build", "cloud_guard"): {"draw_balanced"},
         ("t2_build", "interpolate_clouds"): {"fgw_coupling", "coupling_rotation", "barycentric", "nearest_anchors",
                                              "draw_balanced", "aligned_coupling"}}


@pytest.mark.parametrize(("home", "core"), MOVED.items(), ids=[f"{m}.{n}" for m, n in MOVED])
def test_old_names_are_the_core_objects(op, home, core):
    module, old = home
    assert getattr(importlib.import_module(module), old) is getattr(op, core)
    assert getattr(op, core).__module__ == "vec.core"


def test_heads_call_the_core_and_define_none_of_it(op):
    for (head, fn), core in CALLS.items():
        module = importlib.import_module(f"vec.{head}")
        tree = ast.parse(Path(module.__file__).read_text())
        defined = {n.name for n in ast.walk(tree) if isinstance(n, ast.FunctionDef)}
        assert not defined & {*op.__all__, *(old for _, old in MOVED)}, head
        node = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == fn)
        called = {c.func.id for c in ast.walk(node) if isinstance(c, ast.Call) and isinstance(c.func, ast.Name)}
        assert core <= called, (head, fn, core - called)
        for name in core:
            assert getattr(module, name) is getattr(op, name), (head, name)


def test_core_api_and_imports(op):
    assert all(callable(getattr(op, name)) for name in op.__all__)
    tree = ast.parse(Path(op.__file__).read_text())
    imported = {n.module for n in ast.walk(tree) if isinstance(n, ast.ImportFrom)}
    imported |= {a.name for n in ast.walk(tree) if isinstance(n, ast.Import) for a in n.names}
    assert not {m for m in imported if m and ("t1_scope" in m or "t2_build" in m or "t3_decomp" in m)}, imported
    code = "import sys, vec.core; print(sorted(m for m in sys.modules if m.startswith(('vec.t1_', 'vec.t2_'))))"
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=True, cwd=REPO)
    assert out.stdout.strip() == "[]", out.stdout  # importing the core alone loads neither head


# ---------------------------------------------------------------------------------------------------------------------
# 4. No stage constants in the core
# Every numeric literal allowed in vec/core.py: 0 and 1 (indices, axes, the factor and ratio of no change), 5 (the
# genes an error message lists) and tolerances (1e-9 for floors, 1e-6 for whole numbers, 1e-2 for count recovery). A
# size, a count or a composition share of any stage fails.
OPERATOR_LITERALS = {0.0, 1.0, 5.0, 1e-9, 1e-6, 1e-2,
                     20.0, 256.0}  # 2026-11-03, herd_rows: proposals per drawn row and the batch size (computation)


def test_core_holds_no_stage_constants(op):
    src = Path(op.__file__).read_text()
    lits = {float(t.string) for t in tokenize.generate_tokens(io.StringIO(src).readline) if t.type == tokenize.NUMBER}
    assert lits <= OPERATOR_LITERALS, sorted(lits - OPERATOR_LITERALS)


def released_class_names() -> set[str]:
    """Every chamber category and released label (configs/t1_chamber_classes.yaml) and every lineage group and label
    (configs/lineage_map.yaml)."""
    chamber = yaml.safe_load((REPO / "configs" / "t1_chamber_classes.yaml").read_text())
    lineage = yaml.safe_load((REPO / "configs" / "lineage_map.yaml").read_text())
    names = set(chamber["categories"]) | set(chamber["labels"]) | set(lineage["groups"]) | set(lineage["labels"])
    for per_file in (chamber.get("file_overrides") or {}, lineage.get("file_overrides") or {}):
        for labels in per_file.values():
            names |= set(labels)
    return {str(x) for x in names}


def test_core_names_no_class_stage_or_class_list(op):
    tree = ast.parse(Path(op.__file__).read_text())
    docstrings = {id(n.body[0].value) for n in ast.walk(tree)
                  if isinstance(n, ast.Module | ast.FunctionDef | ast.ClassDef) and n.body
                  and isinstance(n.body[0], ast.Expr) and isinstance(n.body[0].value, ast.Constant)}
    strings = [n.value for n in ast.walk(tree)
               if isinstance(n, ast.Constant) and isinstance(n.value, str) and id(n) not in docstrings]
    names = released_class_names()
    assert len(names) > 50
    for s in strings:
        named = [c for c in names if re.search(rf"(?<![A-Za-z0-9_]){re.escape(c)}(?![A-Za-z0-9_])", s)]
        assert not named, (s, named)
        assert not re.search(r"E\d+(?:\.\d+)?", s), s  # a stage name
    for node in tree.body:  # no module-level list, tuple, set or dict of strings (a class list), except __all__
        if not isinstance(node, ast.Assign | ast.AnnAssign) or node.value is None:
            continue
        targets = node.targets if isinstance(node, ast.Assign) else [node.target]
        if [getattr(t, "id", None) for t in targets] == ["__all__"]:
            continue
        for x in ast.walk(node.value):
            if isinstance(x, ast.List | ast.Tuple | ast.Set | ast.Dict):
                assert not any(isinstance(y, ast.Constant) and isinstance(y.value, str) for y in ast.walk(x)), \
                    ast.unparse(node)


if __name__ == "__main__":  # prints the table, for an output change that is meant (see the module docstring)
    import tempfile

    with tempfile.TemporaryDirectory() as tmp:
        sp = chamber_spec(Path(tmp))
        print("GOLDEN: dict[str, dict[str, str]] = {")
        for case, fn in CASES.items():
            print(f'    "{case}": {{')
            for part, digest in fn(sp).items():
                print(f'        "{part}": "{digest}",')
            print("    },")
        print("}")


# ---------------------------------------------------------------------------------------------------------------------
# 5. Coupling and displacement (feat-010)
def _rotation(a: float, b: float) -> np.ndarray:
    rz = np.array([[np.cos(a), -np.sin(a), 0.0], [np.sin(a), np.cos(a), 0.0], [0.0, 0.0, 1.0]])
    rx = np.array([[1.0, 0.0, 0.0], [0.0, np.cos(b), -np.sin(b)], [0.0, np.sin(b), np.cos(b)]])
    return rz @ rx


def test_coupling_rotation_recovers_a_proper_rotation_and_never_returns_a_mirror(op):
    rng = np.random.default_rng(0)
    A = rng.normal(size=(200, 3)) * [3.0, 2.0, 1.0]
    A -= A.mean(0)
    R0 = _rotation(0.7, -0.4)
    T = np.eye(200) / 200
    R = op.coupling_rotation(A, A @ R0.T, T)
    np.testing.assert_allclose(R, R0, atol=1e-8)
    Rm = op.coupling_rotation(A, A * [1.0, 1.0, -1.0], T)  # a mirror image: the best proper rotation, det +1
    assert np.linalg.det(Rm) == pytest.approx(1.0)
    with pytest.raises(ValueError, match="shape"):
        op.coupling_rotation(A, A, np.eye(199))


def test_fgw_coupling_pairs_a_moved_permuted_copy_cell_for_cell(op):
    rng = np.random.default_rng(1)
    A = rng.normal(size=(40, 3)) * [3.0, 2.0, 1.0]
    F = rng.normal(size=(40, 5))
    perm = rng.permutation(40)
    B, FB = (A @ _rotation(1.1, 0.3).T + 7.0)[perm], F[perm]  # rotated, moved, rows shuffled
    T = op.fgw_coupling(F, FB, A, B, alpha=0.5, max_iter=200, emd_max_iter=100_000)
    assert T.shape == (40, 40) and T.sum() == pytest.approx(1.0)
    np.testing.assert_allclose(T.sum(1), 1 / 40, atol=1e-9)
    np.testing.assert_allclose(T.sum(0), 1 / 40, atol=1e-9)
    assert (T.argmax(1) == np.argsort(perm)).all()  # each cell goes to its own copy
    with pytest.raises(ValueError, match="alpha"):
        op.fgw_coupling(F, FB, A, B, alpha=1.5, max_iter=10, emd_max_iter=10)
    with pytest.raises(ValueError, match="feature row"):
        op.fgw_coupling(F[:-1], FB, A, B, alpha=0.5, max_iter=10, emd_max_iter=10)


def test_aligned_coupling_keeps_each_cell_on_its_own_side_and_reduces_to_each_cost(op):
    import ot

    rng = np.random.default_rng(2)
    side = rng.normal(size=(30, 3)) * [1.0, 2.0, 0.7]
    side[:, 0] = np.abs(side[:, 0]) + 0.2
    A = np.vstack([side, side * [-1.0, 1.0, 1.0]])  # one side of a body and its mirror image
    F = np.c_[np.abs(A[:, 0]), A[:, 1], A[:, 2]]  # expression the same on both sides: features cannot tell them apart
    perm = rng.permutation(60)
    B, FB = A[perm] + rng.normal(scale=0.02, size=(60, 3)), F[perm]  # the same cells in the same frame, rows shuffled
    T = op.aligned_coupling(F, FB, A, B, alpha=0.5, emd_max_iter=100_000)
    assert T.shape == (60, 60) and T.sum() == pytest.approx(1.0)
    np.testing.assert_allclose(T.sum(1), 1 / 60, atol=1e-12)
    np.testing.assert_allclose(T.sum(0), 1 / 60, atol=1e-12)
    assert ((T > 0).sum(1) == 1).all() and ((T > 0).sum(0) == 1).all()  # one to one
    assert (T.argmax(1) == np.argsort(perm)).all()  # each cell goes to its own copy, never to its mirror twin
    M, D = ot.dist(F, FB), ot.dist(A, B)
    M, D = M / M.max(), D / D.max()
    p = np.full(60, 1 / 60)
    for alpha, C in ((0.0, M), (1.0, D)):  # alpha 0: the features alone; alpha 1: the positions alone
        Ta = op.aligned_coupling(F, FB, A, B, alpha=alpha, emd_max_iter=100_000)
        assert (Ta * C).sum() == pytest.approx(ot.emd2(p, p, C, numItermax=100_000), abs=1e-12)
    with pytest.raises(ValueError, match="alpha"):
        op.aligned_coupling(F, FB, A, B, alpha=-0.1, emd_max_iter=10)
    with pytest.raises(ValueError, match="feature row"):
        op.aligned_coupling(F[:-1], FB, A, B, alpha=0.5, emd_max_iter=10)
    with pytest.raises(ValueError, match="dimensions"):
        op.aligned_coupling(F, FB, A, B[:, :2], alpha=0.5, emd_max_iter=10)


def test_barycentric_and_nearest_anchors(op):
    T = np.array([[0.25, 0.25], [0.0, 0.5]])
    np.testing.assert_allclose(op.barycentric(T, np.array([[0.0, 0.0], [2.0, 4.0]])), [[1.0, 2.0], [2.0, 4.0]])
    with pytest.raises(ValueError, match="no mass"):
        op.barycentric(np.array([[0.0, 0.0]]), np.zeros((2, 2)))
    d, i = op.nearest_anchors(np.array([[0.1, 0, 0], [5, 5, 5]]), np.array([[0, 0, 0], [5, 5, 4], [9, 9, 9]]), 2)
    assert i[:, 0].tolist() == [0, 1] and d.shape == (2, 2) and d[0, 0] == pytest.approx(0.1)




# ---------------------------------------------------------------------------------------------------------------------
# herd_rows (2026-11-03, T1 try T5 of experiments/harness/TRIES.md)
def test_herd_rows_brings_the_draw_mean_to_the_pool_mean_and_keeps_its_shape(op):
    rng = np.random.default_rng(40)
    X = sparse.csr_matrix(rng.poisson(rng.gamma(1.0, 1.0, size=60), size=(3000, 60)).astype(float))
    pool = np.arange(500, 3000)
    picked = rng.choice(pool, 400, replace=False)
    new, rec = op.herd_rows(X, pool, picked, np.random.default_rng(1))
    assert len(new) == 400 and len(np.unique(new)) == 400 and np.isin(new, pool).all()
    goal = np.asarray(X[pool].mean(0)).ravel()
    before = np.asarray(X[picked].mean(0)).ravel() - goal
    after = np.asarray(X[new].mean(0)).ravel() - goal
    assert rec["sq_distance_after"] == pytest.approx(float(after @ after))
    assert after @ after < 0.2 * (before @ before) and rec["accepted"] > 0
    kept = new == picked  # a swapped-in row takes the place of the row it replaces
    assert kept.any() and (~kept).any()
    same, rec2 = op.herd_rows(X, pool, np.r_[picked[:10], picked[:10]], np.random.default_rng(1))
    assert np.array_equal(same, np.r_[picked[:10], picked[:10]]) and "skipped" in rec2  # repeats: left alone
    whole, rec3 = op.herd_rows(X, pool, pool, np.random.default_rng(1))
    assert np.array_equal(whole, pool) and "skipped" in rec3  # nothing to swap in
