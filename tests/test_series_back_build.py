"""Row R48's file T2embryo_s2x_ovsb.h5ad: the builder's off-by-default `--expression series_back`
(experiments/builds/build_test.py embryo; vec.t2_build.series_back_expression, the mirror of the series rendering), on
synthetic stages and a synthetic test board whose limits are never shown. Nothing is read from data/ or boards/, a
synthetic TRIES.md stands in for the ledger (one test asks, by code, whether the real ledger holds row R48 as the
builder needs it), and nothing is scored.

Pinned: the factor runs from E8.0 to E7.5 on the series rendering's matched extraembryonic labels and is weighted by
E8.0's label shares; it is applied in count space as the series rendering applies its factor, then the E8.0 rows are
thinned by RENDER_RULE's lo ratio (E7.5's median counts per cell over E8.0's), step by step from the seed; the E7.5
rows stay as released; on O's command the file keeps O's rows, order, positions, count and names and changes only the
E8.0 rows' expression; every refusal comes before anything is read and shows no count; and every other embryo command
gives what the builder of BASE (main before the option) gives, run in this same process from git's copy of its two
changed files: the same h5ad bytes and record, or the same refusal. Built through experiments/builds/build_test.py."""

from __future__ import annotations

import hashlib
import json
import re
import sys
import types
from functools import partial
from pathlib import Path
from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from conftest import git
from scipy import sparse

from vec import core
from vec import t2_build as tb
from vec.boards import Board

REPO = Path(__file__).resolve().parents[1]
BASE = "9494483f1eea66f7b6f3984ed00ad043b2d94afa"  # main before --expression series_back existed
PANEL = [f"g{i}" for i in range(12)]
EXTRA = ["x0", "x1", "x2"]  # genes outside the panel: part of each cell's library, never of a file
SECRET_MAX, SECRET_TRUTH = 400, 120  # stand for the test board's max_cells and truth_cells: never shown
# label: (lineage group, cells at E7.25 and E7.5, cells at E8.0, spatial cluster)
WORLD = {"EXE-Endoderm": ("Endoderm", 300, 250, 0), "ExM": ("ExtraEmbryonicMeso", 300, 250, 1),
         "Epi": ("Epiblast", 300, 200, 2), "Meso": ("Mesoderm", 300, 200, 2), "Neural": ("Ectoderm", 0, 100, 1)}
BUMP = {"EXE-Endoderm": [0, 4], "ExM": [1, 5], "Epi": [8, 10], "Meso": [2, 9], "Neural": [3, 11], "EXE-A": [0, 1],
        "EXE-B": [6, 7]}
RUN = np.r_[np.full(4, 2.0), np.full(4, 0.5), np.ones(4), np.ones(3)]  # the run's per-gene factor, E7.5 -> E8.0
DEV = {"Meso": np.r_[np.ones(8), np.full(2, 4.0), np.ones(5)]}  # an embryonic label's own development at E8.0
DEPTH = {"E7.25": 50.0, "E7.5": 60.0, "E8.0": 200.0}
SHIFT = {"E7.25": -0.1, "E7.5": 0.0, "E8.0": 0.3}
FLOOR = ["--method", "fgw_shape0", "--positions", "floor_exact"]
OV = ["--shape-interp", "overlay", "--shape-interp-test-ok"]
O_ARGS = [*FLOOR, *OV]
W_ARGS = [*FLOOR, "--expression", "series_back", *OV, "--captain-ruling", "R48"]
R48_ROW = ("| R48 (the captain's order of 2026-11-24, ruling, no run: \"Build T2embryo_s2x_ovsb.h5ad as named\"; the "
           "embryo builder is to take `--expression series_back` on O's command only with `--captain-ruling R48`; not "
           "a try) | T2 embryo | synthetic |")


# ------------------------------------------------------------------------------------------------ the synthetic world
def _stage(name: str, seed: int, depth: float | None = None, world=None, run=None, bump=None,
           dev=None) -> ad.AnnData:
    """A bracket as a released file holds it after the duplicate filter: Poisson counts over the panel and three genes
    outside it, X = log(1 + 1e4 c / N) with N over every gene, obs celltype and group, three spatial clusters. Each
    label's profile has its own bumps (BUMP); at E8.0 every label carries the run's factor (`run`: one vector, or one
    per label) and its own development (DEV)."""
    rng = np.random.default_rng(seed)
    world = WORLD if world is None else world
    depth = DEPTH[name] if depth is None else depth
    run = RUN if run is None else run
    bump = BUMP if bump is None else bump
    dev = DEV if dev is None else dev
    centres = np.array([[-4.0, 0.0, 0.0], [0.0, 3.0, 0.0], [4.0, 0.0, 1.0]]) * (1.0 + SHIFT.get(name, 0.0))
    C, labs, groups, P = [], [], [], []
    for lab, (grp, n_early, n_late, cl) in world.items():
        n = n_late if name == "E8.0" else n_early
        if not n:
            continue
        prof = np.r_[np.ones(8), np.full(4, 0.3), np.ones(3)]
        prof[bump.get(lab, [])] += 2.0
        if name == "E8.0":
            prof = prof * (run[lab] if isinstance(run, dict) else run) * dev.get(lab, 1.0)
        c = rng.poisson(prof / prof.sum() * depth, size=(n, len(prof)))
        c[c[:, :len(PANEL)].sum(1) == 0, 0] += 1  # every cell has a panel count, as every released cell does
        C.append(c)
        labs += [lab] * n
        groups += [grp] * n
        P.append(centres[cl] + rng.normal(size=(n, 3)) * [1.5, 1.0, 0.5])
    order = rng.permutation(sum(len(c) for c in C))
    return _anndata(name, np.vstack(C)[order], np.asarray(labs)[order], np.asarray(groups)[order],
                    np.vstack(P)[order])


def _anndata(name, C, labs, groups, P) -> ad.AnnData:
    N = C.sum(1)
    a = ad.AnnData(X=sparse.csr_matrix(np.log1p(1e4 * C / N[:, None]).astype(np.float32)),
                   obs=pd.DataFrame({"celltype": labs, "group": groups}, index=[f"{name}_{i}" for i in range(len(C))]),
                   var=pd.DataFrame(index=PANEL + EXTRA))
    a.layers["counts"] = sparse.csr_matrix(C.astype(np.int32))
    a.obsm["spatial_3D"] = np.asarray(P, dtype=np.float32)
    return a


@pytest.fixture(scope="module")
def stages():
    return {"E7.25.h5ad": _stage("E7.25", 11), "E7.5.h5ad": _stage("E7.5", 12), "E8.0.h5ad": _stage("E8.0", 13)}


def _board(key="T2:embryo:test_interp", split="test_interp") -> Board:
    return Board(key=key, task="T2", genes=list(PANEL), needs_coords=True, split=split, setting="embryo",
                 min_cells=50, max_cells=SECRET_MAX, truth_cells=SECRET_TRUTH)


def install(mod, setattr_, setitem_, shared: Path, reads: list, stages: dict) -> None:
    """Point a build_test module (this checkout's or BASE's) at the synthetic stages, the synthetic test and validation
    boards, a synthetic size JSON and a synthetic ledger, both in `shared`; the interpolation runs with 150 anchors (the
    code path is the builder's)."""
    boards = {b.key: b for b in (_board(), _board("T2:embryo:val_interp", "val_interp"))}

    def load(path, genes=None, drop_duplicates=False):
        name = Path(path).name
        assert Path(path).parent.name == "T2" and genes == PANEL and drop_duplicates, path
        reads.append(name)
        return stages[name].copy()

    def read(path, genes=None):
        reads.append("released:" + Path(path).name)
        a = stages[Path(path).name]
        return (a if genes is None else a[:, list(genes)]).copy()

    shared.mkdir(parents=True, exist_ok=True)
    size = shared / "size_embryo.json"
    size.write_text(json.dumps({"board": "embryo", "valid": True,
                                "preregistered": {"z_mode": "raw", "valid": True, "pred_rms_full": 6.0}}))
    ledger = shared / "TRIES.md"
    ledger.write_text("| ID | what |\n|---|---|\n| E56 (step 0) | T2 embryo |\n" + R48_ROW + "\n")
    setattr_(mod, "load_stage", load)
    setattr_(mod, "read_stage", read)
    setattr_(mod, "load_boards", lambda: boards)
    setattr_(mod, "_t2_count", lambda b, n: {"decision": "synthetic"})
    setattr_(mod.t2_build, "interpolate_clouds", partial(mod.t2_build.interpolate_clouds, n_anchor=150))
    setitem_(mod.SIZES, "embryo", size)
    setattr_(mod, "TRIES_FILE", ledger)


@pytest.fixture(scope="module")
def bt():
    sys.path.insert(0, str(REPO / "experiments" / "builds"))
    import build_test

    return build_test


@pytest.fixture
def env(bt, monkeypatch, own_repo, tmp_path, stages):
    reads: list[str] = []
    install(bt, monkeypatch.setattr, monkeypatch.setitem, tmp_path / "shared", reads, stages)
    return reads


def _shows(text: str, *values: int) -> bool:
    """Whether any of the integers appears in the text as a number of its own (not as digits inside a float)."""
    return any(re.search(rf"(?<![\d.]){v}(?![\d.])", text) for v in values)


def _build(mod, tmp_path: Path, name: str, *argv: str) -> tuple[Path, dict, str]:
    out = tmp_path / name / "T2embryo.h5ad"
    mod.main(["embryo", *argv, "--out", str(out)])
    text = out.with_suffix(".build.json").read_text()
    return out, json.loads(text), text


def _file(path: Path) -> tuple[np.ndarray, np.ndarray, list, list]:
    a = ad.read_h5ad(path)
    X = np.asarray(a.X.toarray() if sparse.issparse(a.X) else a.X)
    return X, np.asarray(a.obsm["spatial_3D"]), list(a.obs_names), list(a.var_names)


def _mix(before, after, n_before=150, n_after=180, seed=5):
    """A stand-in for the interpolation: drawn rows of each bracket, interleaved."""
    rng = np.random.default_rng(seed)
    b = rng.choice(before.n_obs, n_before, replace=False)
    a = rng.choice(after.n_obs, n_after, replace=False)
    order = rng.permutation(n_before + n_after)
    return SimpleNamespace(rows=np.r_[b, a][order], stage=np.r_[np.zeros(n_before, int), np.ones(n_after, int)][order])


def _props(C, rows):
    sub = sparse.csr_matrix(C[rows], dtype=np.float64)
    tot = np.asarray(sub.sum(1)).ravel()
    return np.asarray((sparse.diags(1.0 / np.maximum(tot, 1.0)) @ sub).mean(0)).ravel()


# ------------------------------------------------------------------------------------------------ the expression step
def test_the_factor_runs_from_e80_to_e75_on_the_matched_extraembryonic_labels(stages):
    before, after = stages["E7.5.h5ad"], stages["E8.0.h5ad"]
    phi, labels = tb.series_back_factor(before, after, PANEL)
    # the series rendering's labels: both brackets hold them with 100 cells, extraembryonic by name or by group
    assert labels == ["EXE-Endoderm", "ExM"] == tb.series_labels(before.obs["celltype"], before.obs["group"],
                                                                 after.obs["celltype"])
    # the run's factor reversed: genes the run doubled are halved and the other way round
    assert np.all(phi[:4] < 0.7) and np.all(phi[4:8] > 1.4)
    so = tb.series_factor(tb._stage_counts(before, PANEL)[1], before.obs["celltype"],
                          tb._stage_counts(after, PANEL)[1], after.obs["celltype"], labels)
    assert np.corrcoef(np.log(phi), np.log(so))[0, 1] < -0.95  # the series rendering's factor, the other way
    # the embryonic label's own development (Meso's genes 8 and 9, fourfold) is not picked up
    assert np.all(np.abs(np.log(phi[8:10])) < 0.4)


def test_the_factor_is_weighted_by_e80s_label_shares():
    """Two extraembryonic labels with their own run factors, whose shares swap between the brackets: the factor is
    gene_factor with E8.0's label shares as weights, not E7.5's."""
    world = {"EXE-A": ("Endoderm", 400, 100, 0), "EXE-B": ("Endoderm", 100, 400, 1), "Epi": ("Epiblast", 200, 200, 2)}
    run = {"EXE-A": np.r_[np.full(4, 3.0), np.ones(11)], "EXE-B": np.r_[np.ones(4), np.full(4, 3.0), np.ones(7)],
           "Epi": np.ones(15)}
    before, after = _stage("E7.5", 21, world=world), _stage("E8.0", 22, world=world, run=run)
    phi, labels = tb.series_back_factor(before, after, PANEL)
    assert labels == ["EXE-A", "EXE-B"]
    C75, C80 = tb._stage_counts(before, PANEL)[1], tb._stage_counts(after, PANEL)[1]
    n75, n80 = (np.asarray(x.obs["celltype"]).astype(str) for x in (before, after))

    def factor(weights_from):
        w = np.array([(weights_from == lab).sum() for lab in labels], float)
        w /= w.sum()
        num = sum(wi * _props(C75, np.flatnonzero(n75 == lab)) for wi, lab in zip(w, labels, strict=True))
        den = sum(wi * _props(C80, np.flatnonzero(n80 == lab)) for wi, lab in zip(w, labels, strict=True))
        f = (num + tb.SERIES_EPS) / (den + tb.SERIES_EPS)
        return np.clip(f / (den @ f / den.sum()), *tb.SERIES_CLIP)

    np.testing.assert_allclose(phi, factor(n80), rtol=1e-12)
    assert not np.allclose(phi, factor(n75), rtol=1e-2)
    assert phi[4:8].mean() < phi[:4].mean()  # E8.0's dominant label (EXE-B, genes 4-7 up) is pulled back the most


def test_the_moved_rows_are_the_count_space_factor_then_the_lo_thinning_step_by_step(stages):
    before, after = stages["E7.5.h5ad"], stages["E8.0.h5ad"]
    mix = _mix(before, after)
    X, rec = tb.series_back_expression(before, after, mix, PANEL, 0.5, seed=[3, 7])
    # by hand, in the rule's order: the factor in count space as the series rendering applies it, then RENDER_RULE's
    # lo thinning with its retry, on one random stream
    rng = np.random.default_rng([3, 7])
    E75, C75, lib75, _ = tb._stage_counts(before, PANEL)
    E80, C80, lib80, oth80 = tb._stage_counts(after, PANEL)
    d75, d80 = float(np.median(lib75)), float(np.median(lib80))
    phi, _ = tb.series_back_factor(before, after, PANEL)
    a = mix.stage == 1
    rows = mix.rows[a]
    S = tb._smoothed_props(E80, C80, lib80, rows, tb.RENDER_K)
    C1 = sparse.csr_matrix(core.apply_gene_factor(np.asarray(C80[rows].toarray(), np.int32), phi,
                                                  S / S.sum(1, keepdims=True), rng).astype(np.int64))
    w = np.full(len(rows), d75 / d80)

    def make(idx, g):
        Cn, on = tb._binomial_rows(C1[idx], oth80[rows][idx], w[idx], g)
        return Cn, np.asarray(Cn.sum(1)).ravel() + on

    Cn, ln, _ = tb._keep_some_counts(make, len(rows), rng)
    np.testing.assert_array_equal(X[a], tb._log_normalise(Cn, ln))
    assert (Cn.toarray() <= C1.toarray()).all()  # the thinning only removes counts
    np.testing.assert_array_equal(X[~a], E75[mix.rows[~a]])  # the E7.5 rows as released, bit for bit
    assert rec["ratio"] == d75 / d80 < 1 and rec["moved"] == "after" and rec["labels"] == ["EXE-Endoderm", "ExM"]
    assert rec["median_depth"] == {"before": d75, "after": d80} and rec["rule"] == tb.SERIES_BACK_RULE
    assert rec["thinned"] == {"before": 0, "after": int(a.sum())} and rec["rendered_median"] == float(np.median(ln))
    assert rec["sd_log_factor"] == float(np.log(phi).std())
    assert rec["genes_above_2x"] == int((np.abs(np.log(phi)) > np.log(2.0)).sum())


def test_the_moved_rows_sit_at_e75s_depth_and_in_its_run_profile(stages):
    before, after = stages["E7.5.h5ad"], stages["E8.0.h5ad"]
    mix = _mix(before, after, 200, 600)
    X, rec = tb.series_back_expression(before, after, mix, PANEL, 0.5, seed=[1, 7])
    a = mix.stage == 1
    E80 = tb._stage_counts(after, PANEL)[0][mix.rows[a]]
    assert rec["rendered_median"] == pytest.approx(rec["median_depth"]["before"], rel=0.12)
    assert (X[a] > 0).sum(1).mean() < (E80 > 0).sum(1).mean()  # thinned: fewer detected genes than as released
    assert np.isfinite(X).all() and (X >= 0).all() and ((X > 0).sum(1) >= 1).all()
    c, c0 = np.expm1(X[a].astype(np.float64)), np.expm1(E80.astype(np.float64))
    down = (c[:, :4].sum(1) / c.sum(1)).mean() / (c0[:, :4].sum(1) / c0.sum(1)).mean()
    up = (c[:, 4:8].sum(1) / c.sum(1)).mean() / (c0[:, 4:8].sum(1) / c0.sum(1)).mean()
    assert down < 0.8 and up > 1.2  # the run's doubled genes fall back, its halved genes rise
    # reproducible from its seed; another seed changes the moved rows only
    X2, _ = tb.series_back_expression(before, after, mix, PANEL, 0.5, seed=[1, 7])
    np.testing.assert_array_equal(X, X2)
    X3, _ = tb.series_back_expression(before, after, mix, PANEL, 0.5, seed=[2, 7])
    np.testing.assert_array_equal(X3[~a], X[~a])
    assert not np.array_equal(X3[a], X[a])


def test_equal_depths_apply_the_factor_alone_and_the_record_keeps_cell_numbers_in_its_dicts(stages):
    """The later bracket: the earlier one's cells with panel genes 0 and 1 swapped with 8 and 9, so every cell keeps its
    library and the medians are equal: nothing is thinned, the factor alone moves the rows."""
    before = stages["E7.5.h5ad"]
    C = np.asarray(before.layers["counts"].toarray())
    swap = np.r_[8, 9, 2, 3, 4, 5, 6, 7, 0, 1, 10, 11, 12, 13, 14]
    after = _anndata("late", C[:, swap], before.obs["celltype"].to_numpy(), before.obs["group"].to_numpy(),
                     before.obsm["spatial_3D"])
    mix = _mix(before, after, 100, 120)
    X, rec = tb.series_back_expression(before, after, mix, PANEL, 0.5, seed=[4, 7])
    assert rec["ratio"] == 1.0 and rec["thinned"] == {"before": 0, "after": 0}
    rng = np.random.default_rng([4, 7])
    E, Ca, lib, oth = tb._stage_counts(after, PANEL)
    phi, _ = tb.series_back_factor(before, after, PANEL)
    assert np.abs(np.log(phi)).max() > np.log(1.5)  # a factor that moves something
    a = mix.stage == 1
    rows = mix.rows[a]
    S = tb._smoothed_props(E, Ca, lib, rows, tb.RENDER_K)
    C1 = core.apply_gene_factor(np.asarray(Ca[rows].toarray(), np.int32), phi, S / S.sum(1, keepdims=True), rng)
    C1 = sparse.csr_matrix(C1.astype(np.int64))
    np.testing.assert_array_equal(X[a], tb._log_normalise(C1, np.asarray(C1.sum(1)).ravel() + oth[rows]))
    flat = {k: v for k, v in rec.items() if not isinstance(v, (dict, list))}
    assert set(flat) == {"rule", "t", "moved", "ratio", "sd_log_factor", "genes_above_2x", "rendered_median"}


def test_the_option_is_refused_when_the_earlier_bracket_is_the_deeper_or_no_label_matches(stages):
    before, after = stages["E7.5.h5ad"], stages["E8.0.h5ad"]
    mix = _mix(before, after)
    with pytest.raises(ValueError, match="earlier bracket is the shallower"):
        tb.series_back_expression(after, before, mix, PANEL, 0.5, seed=1)
    world = {"Epi": ("Epiblast", 300, 300, 2)}
    lone, lone_late = _stage("E7.5", 41, world=world), _stage("E8.0", 42, world=world)
    with pytest.raises(ValueError, match="no matched extraembryonic label"):
        tb.series_back_expression(lone, lone_late, _mix(lone, lone_late, 50, 50), PANEL, 0.5, seed=1)
    # the series rendering itself is untouched by the option
    _, rec = tb.series_rendered_expression(before, after, mix, PANEL, 0.5, seed=[3, 7])
    assert rec["moved"] == "before" and rec["rule"] == tb.SERIES_RULE


# ------------------------------------------------------------------------------------------------ the builder
def _capture(bt, monkeypatch) -> list:
    """Record each interpolation the builder makes (its rows and brackets, in file order)."""
    seen = []
    inner = bt.t2_build.interpolate_clouds

    def spy(*a, **kw):
        res = inner(*a, **kw)
        seen.append(res)
        return res

    monkeypatch.setattr(bt.t2_build, "interpolate_clouds", spy)
    return seen


def test_on_os_command_the_file_is_os_but_for_the_e80_rows_expression(bt, env, tmp_path, monkeypatch, capsys):
    seen = _capture(bt, monkeypatch)
    o, rec_o, _ = _build(bt, tmp_path, "o", *O_ARGS)
    w, rec_w, text = _build(bt, tmp_path, "w", *W_ARGS)
    Xo, Co, no, vo = _file(o)
    Xw, Cw, nw, vw = _file(w)
    assert Cw.dtype == Co.dtype and Cw.tobytes() == Co.tobytes()  # O's positions, bit for bit, in O's order
    assert nw == no and vw == vo and len(Xw) == len(Xo)  # O's names, genes and count
    np.testing.assert_array_equal(seen[1].rows, seen[0].rows)
    np.testing.assert_array_equal(seen[1].stage, seen[0].stage)
    early = seen[0].stage == 0
    np.testing.assert_array_equal(Xw[early], Xo[early])  # the E7.5 rows as released, as in O
    assert (np.abs(Xw[~early] - Xo[~early]).max(1) > 0).mean() > 0.9  # the E8.0 rows moved
    # the record: O's with the expression step named, the series_back block, no cell number
    r = rec_w["record"]
    assert r["expression"] == "series_back" and r["series_back_moved"] == "after"
    assert r["series_back_rule"] == tb.SERIES_BACK_RULE and 0 < r["series_back_ratio"] < 1
    assert not any(k.startswith("series_back_") and k.endswith(("thinned", "retried", "labels")) for k in r)
    blk = rec_w["series_back"]
    assert blk["captain_ruling"] == "R48" and blk["file"] == "T2embryo_s2x_ovsb.h5ad"
    assert blk["command"] == bt.EMBRYO_OVSB_COMMAND == "embryo " + " ".join(W_ARGS)
    assert blk["o_command"] == "embryo " + " ".join(O_ARGS) and "E56" in blk["step0"]
    assert all(not isinstance(v, (dict, list)) for v in blk.values())
    assert rec_w["shape_interp"] == rec_o["shape_interp"] and rec_w["n_cells"] == rec_o["n_cells"] == "max_cells"
    assert {k: v for k, v in r.items() if not k.startswith("series_back_") and k != "expression"} == rec_o["record"]
    assert "series_back" not in json.dumps(rec_o)
    assert not _shows(text + capsys.readouterr().out, SECRET_MAX, SECRET_TRUTH)
    w2, _, _ = _build(bt, tmp_path, "w2", *W_ARGS)  # built again: the same bytes
    assert w2.read_bytes() == w.read_bytes()


def _refused(bt, env, tmp_path, argv, match, name="T2embryo.h5ad"):
    out = tmp_path / "r" / name
    with pytest.raises(SystemExit, match=match) as e:
        bt.main(["embryo", *argv, "--out", str(out)])
    assert not out.exists() and env == []  # refused before any stage was read
    assert not _shows(str(e.value), SECRET_MAX, SECRET_TRUTH)
    return str(e.value)


@pytest.mark.parametrize("argv", [
    [*FLOOR, "--expression", "series_back", *OV],  # no ruling
    [*FLOOR, "--expression", "series_back", *OV, "--captain-ruling", "R47"],
    [*FLOOR, "--expression", "series_back", *OV, "--captain-ruling", "R28"],
])
def test_the_ruling_must_be_r48(bt, env, tmp_path, argv):
    assert "series_back" in _refused(bt, env, tmp_path, argv, "R48")


@pytest.mark.parametrize("argv", [
    ["--method", "fgw", "--expression", "series_back", "--captain-ruling", "R48"],
    ["--method", "fgw_shape0", "--expression", "series_back", *OV, "--captain-ruling", "R48"],
    [*FLOOR, "--expression", "series_back", "--captain-ruling", "R48"],  # no form
    [*FLOOR, "--expression", "series_back", "--shape-interp", "overlay", "--captain-ruling", "R48"],  # no yes flag
    [*FLOOR, "--expression", "series_back", "--shape-interp", "stretch_ratio", "--shape-interp-test-ok",
     "--captain-ruling", "R48"],
    [*W_ARGS, "--row-order", "copy_parity", "--row-order-test-ok"],
    [*W_ARGS, "--herd"],
    [*W_ARGS, "--draw", "midpoint"],
    [*W_ARGS, "--move-coupling", "aligned"],
    [*W_ARGS, "--copy-rows", "--copy-rows-test-ok"],
    [*W_ARGS, "--seed", "1000002"],
])
def test_any_command_but_os_is_refused(bt, env, tmp_path, argv):
    _refused(bt, env, tmp_path, argv, "O's command")


def test_no_dry_run_either(bt, env, tmp_path):
    _refused(bt, env, tmp_path, [*W_ARGS, "--board", "T2:embryo:val_interp"], "O's command", "T2embryo_dry.h5ad")


def test_a_ruling_without_the_option_is_refused(bt, env, tmp_path):
    for argv in ([*O_ARGS, "--captain-ruling", "R48"], ["--method", "fgw", "--captain-ruling", "R48"],
                 [*FLOOR, "--expression", "series", *OV, "--captain-ruling", "R48"]):
        _refused(bt, env, tmp_path, argv, "--captain-ruling goes with --expression series_back only")


@pytest.mark.parametrize("ledger", [
    "| E56 (step 0) | T2 embryo |\n",  # no row R48
    R48_ROW + "\n" + R48_ROW + "\n",  # two rows R48
    R48_ROW.replace("`--expression series_back`", "--expression series_back") + "\n",  # not in backticks
    R48_ROW.replace("`--captain-ruling R48`", "--captain-ruling R48") + "\n",
    R48_ROW.replace("ruling, no run", "ruling") + "\n",
    R48_ROW.replace("(the captain's order", "(the AI assistant's order") + "\n",
    "| R48 (the captain's order of 2026-11-24, ruling, no run) | `--expression series_back` `--captain-ruling R48` |\n",
    R48_ROW.replace("| R48 ", "| R48a ") + "\n",
])
def test_the_ledger_must_hold_exactly_one_row_r48_naming_the_option_and_the_ruling(bt, env, tmp_path, monkeypatch,
                                                                                   ledger):
    p = tmp_path / "ledger" / "TRIES.md"
    p.parent.mkdir()
    p.write_text("| ID | what |\n|---|---|\n" + ledger)
    monkeypatch.setattr(bt, "TRIES_FILE", p)
    _refused(bt, env, tmp_path, W_ARGS, "exactly one row R48")


def test_the_real_ledger_holds_row_r48_as_the_builder_needs_it(bt):
    assert bt._series_back_ordered() is True


def test_the_option_is_off_by_default_and_named_in_the_parser(bt):
    a = bt.build_parser().parse_args(["embryo", *O_ARGS, "--out", "x.h5ad"])
    assert a.expression == "mixture" and a.captain_ruling is None
    src = (REPO / "experiments" / "builds" / "build_test.py").read_text()
    assert src.count('default="mixture"') == 1 and "series_back = _series_back_refusals(a)" in src
    assert bt.EMBRYO_O_COMMAND == "embryo " + " ".join(O_ARGS)


def test_the_launcher_list_names_this_file():
    sh = (REPO / "experiments" / "external" / "n54_step.sh").read_text()
    assert "tests/test_series_back_build.py" in sh


# ------------------------------------------------------------------------------------------------ every other option
def _module(name: str, rel: str) -> types.ModuleType:
    """A module made from BASE's text of `rel`, registered under `name` (dataclasses look their module up there)."""
    mod = types.ModuleType(name)
    mod.__file__ = str(REPO / rel)
    sys.modules[name] = mod
    exec(compile(git(REPO, "show", f"{BASE}:{rel}") + "\n", str(REPO / rel), "exec"), mod.__dict__)
    return mod


@pytest.fixture(scope="module")
def base_bt(bt):
    """BASE's builder with BASE's vec.t2_build in place of this checkout's: the option changes these two files only,
    and every other module is this checkout's in both builders, so a difference below is the option's."""
    t2_base = _module("_vec_t2_build_at_base", "src/vec/t2_build.py")
    mod = _module("_build_test_at_base", "experiments/builds/build_test.py")
    mod.t2_build = t2_base
    assert not hasattr(t2_base, "series_back_expression") and "series_back" not in mod.__doc__
    return mod


DRY = "T2:embryo:val_interp"
# name: (argv, file name, what the builder of BASE does with it: "built", "refused" or None for either)
COMMANDS = {
    "s1": (["--method", "fgw"], "T2embryo.h5ad", "built"),
    "s2": (["--method", "fgw_shape0"], "T2embryo.h5ad", "built"),
    "shape1": (["--method", "fgw_shape1"], "T2embryo.h5ad", "built"),
    "pseudobulk": (["--method", "pseudobulk_shift"], "T2embryo.h5ad", None),
    "P": (FLOOR, "T2embryo.h5ad", "built"),
    "s2xs": ([*FLOOR, "--expression", "series"], "T2embryo.h5ad", "built"),
    "P_hi": ([*FLOOR, "--expression", "render_hi"], "T2embryo.h5ad", "built"),
    "O": (O_ARGS, "T2embryo.h5ad", "built"),
    "SO": ([*FLOOR, "--expression", "series", *OV], "T2embryo.h5ad", "built"),
    "O_hi": ([*FLOOR, "--expression", "render_hi", *OV], "T2embryo.h5ad", "built"),
    "stretch": ([*FLOOR, "--shape-interp", "stretch_ratio", "--shape-interp-test-ok"], "T2embryo.h5ad", "built"),
    "overlay_stretch": ([*FLOOR, "--shape-interp", "overlay_stretch", "--shape-interp-test-ok"], "T2embryo.h5ad",
                        None),
    "depth": (["--method", "fgw", "--expression", "depth"], "T2embryo.h5ad", "built"),
    "blend": (["--method", "fgw", "--expression", "blend"], "T2embryo.h5ad", "built"),
    "series": (["--method", "fgw", "--expression", "series"], "T2embryo.h5ad", "built"),
    "render_hi": (["--method", "fgw", "--expression", "render_hi"], "T2embryo.h5ad", "built"),
    "herd": (["--method", "fgw", "--herd"], "T2embryo.h5ad", "built"),
    "aligned": (["--method", "fgw", "--move-coupling", "aligned"], "T2embryo.h5ad", "built"),
    "midpoint": (["--method", "fgw", "--draw", "midpoint"], "T2embryo.h5ad", "built"),
    "copy_parity": ([*FLOOR, "--row-order", "copy_parity", "--row-order-test-ok"], "T2embryo.h5ad", None),
    "copy_rows": ([*FLOOR, "--copy-rows", "--copy-rows-test-ok"], "T2embryo.h5ad", None),
    "copy_rows_series": ([*FLOOR, "--expression", "series", "--copy-rows", "--copy-rows-test-ok"], "T2embryo.h5ad",
                         None),
    "copy_rows_size": ([*FLOOR, "--copy-rows", "--copy-rows-test-ok", "--copy-rows-size", "--copy-rows-size-test-ok"],
                       "T2embryo.h5ad", None),
    "seed": ([*O_ARGS, "--seed", "1000003"], "T2embryo.h5ad", "built"),
    "dry_O": ([*FLOOR, "--shape-interp", "overlay", "--board", DRY], "T2embryo_dry.h5ad", "built"),
    # the refusals of the options as they were
    "positions": (["--method", "fgw", "--positions", "floor_exact"], "T2embryo.h5ad", "refused"),
    "form_needs_ok": ([*FLOOR, "--shape-interp", "overlay"], "T2embryo.h5ad", "refused"),
    "form_floor_only": (["--method", "fgw_shape0", *OV], "T2embryo.h5ad", "refused"),
    "ok_alone": ([*FLOOR, "--shape-interp-test-ok"], "T2embryo.h5ad", "refused"),
    "row_order_needs_ok": ([*FLOOR, "--row-order", "copy_parity"], "T2embryo.h5ad", "refused"),
    "row_order_ok_alone": ([*FLOOR, "--row-order-test-ok"], "T2embryo.h5ad", "refused"),
    "copy_rows_needs_ok": ([*FLOOR, "--copy-rows"], "T2embryo.h5ad", "refused"),
    "copy_rows_method": (["--method", "fgw", "--copy-rows", "--copy-rows-test-ok"], "T2embryo.h5ad", "refused"),
    "copy_rows_expression": ([*FLOOR, "--expression", "render_hi", "--copy-rows", "--copy-rows-test-ok"],
                             "T2embryo.h5ad", "refused"),
    "copy_rows_form": ([*O_ARGS, "--copy-rows", "--copy-rows-test-ok"], "T2embryo.h5ad", "refused"),
    "size_alone": ([*FLOOR, "--copy-rows-size"], "T2embryo.h5ad", "refused"),
    "size_expression": ([*FLOOR, "--expression", "series", "--copy-rows", "--copy-rows-test-ok", "--copy-rows-size",
                         "--copy-rows-size-test-ok"], "T2embryo.h5ad", "refused"),
    "seed_low": ([*O_ARGS, "--seed", "5"], "T2embryo.h5ad", "refused"),
    "dry_name": ([*FLOOR, "--board", DRY], "T2embryo.h5ad", "refused"),
    "test_name": (FLOOR, "T2embryo_dry.h5ad", "refused"),
    "unknown_expression": ([*FLOOR, "--expression", "nonsense"], "T2embryo.h5ad", "refused"),
}


def _outcome(mod, tmp: Path, stages, which: str, argv: list[str], name: str):
    """What one builder does with one command: the h5ad's sha256 and its record's text, or ("refused", its message)."""
    with pytest.MonkeyPatch.context() as mp:
        install(mod, mp.setattr, mp.setitem, tmp / "shared", [], stages)
        out = tmp / which / name
        try:
            mod.main(["embryo", *argv, "--out", str(out)])
        except SystemExit as e:
            return "refused", str(e.code).replace(str(out), "<out>")
        return hashlib.sha256(out.read_bytes()).hexdigest(), out.with_suffix(".build.json").read_text()


@pytest.mark.parametrize("key", list(COMMANDS))
def test_every_other_command_does_what_base_does(bt, base_bt, own_repo, tmp_path, stages, key, capsys):
    argv, name, kind = COMMANDS[key]
    want = _outcome(base_bt, tmp_path, stages, "base", argv, name)
    got = _outcome(bt, tmp_path, stages, "new", argv, name)
    base_did = "refused" if want[0] == "refused" else "built"
    assert kind is None or base_did == kind, (base_did, kind)
    assert got[0] == want[0] and got[1] == want[1]  # the same bytes and record, or the same refusal
