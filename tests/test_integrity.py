"""N14 (experiments/harness/TRIES.md), on synthetic data only: the integrity pass of the ten frozen upload files
(experiments/research_2026_11/integrity.py). A small synthetic archive of the four board kinds raises no flag when it
is correct, and each check catches its named failure on a damaged copy. No project data, board file or build is read.
"""

from __future__ import annotations

import hashlib
import json
import re
import shutil
import sys
from pathlib import Path
from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from vec import provenance, t2_build
from vec.io import COORD_KEY
from vec.prep import log1p_cp10k

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "experiments" / "research_2026_11"))
import integrity as ig  # noqa: E402

BOARDS = {"t1": "T1:test", "heart": "T2:heart:test_extrap", "embryo": "T2:embryo:test_interp", "t3": "T3:bcatenin"}
T1_GENES = [f"G{i:03d}" for i in range(140)]
PANEL = [f"P{i:03d}" for i in range(230)]
# Cell counts: distinctive numbers that must never reach an output of the pass
FILE_CELLS = {"t1": 1733, "heart": 1597, "embryo": 1481, "t3": 1327}
SOURCE_CELLS = {"t1": 2111, "heart": 2003, "e75": 1709, "e80": 1801, "t3": 1907}
FORBIDDEN_FIELDS = ("n_cells", "pred_cells", "truth_cells", "count_log_ratio")
UNION = "ab" * 32


# ---------------------------------------------------------------------------------------------------------------------
# The synthetic archive: released stages, four upload files, build records, SHA256SUMS
def _cloud(n: int, rng) -> np.ndarray:
    return np.c_[rng.uniform(-40, 40, n), rng.uniform(-40, 40, n), rng.choice(np.arange(8) * 2.5, n)].astype(np.float64)


def _field(n_genes: int, rng) -> SimpleNamespace:
    """A smooth spatial field per gene (first 60% of the genes) and a flat level for the rest."""
    return SimpleNamespace(freq=rng.uniform(0.05, 0.2, (n_genes, 2)), phase=rng.uniform(0, 2 * np.pi, n_genes),
                           base=np.exp(rng.normal(0.0, 0.7, n_genes)),
                           amp=np.where(np.arange(n_genes) < 0.6 * n_genes, 1.4, 0.0))


def _spatial_counts(coords: np.ndarray, f: SimpleNamespace, rng) -> np.ndarray:
    mu = f.base * np.exp(f.amp * np.cos(coords[:, :2] @ f.freq.T + f.phase))
    return rng.poisson(mu).astype(np.int32)


def write_stage(path: Path, counts: np.ndarray, coords: np.ndarray | None, genes: list[str]) -> None:
    """A released stage as vec.prep.load_stage expects it: log1p CP10k in X, counts in a layer, celltype in obs."""
    n = len(counts)
    obs = pd.DataFrame({"celltype": ["Unknown"] * n}, index=[f"c{i}" for i in range(n)])
    a = ad.AnnData(X=log1p_cp10k(counts), obs=obs, var=pd.DataFrame(index=list(genes)))
    a.layers["counts"] = sparse.csr_matrix(counts)
    if coords is not None:
        a.obsm[COORD_KEY] = np.asarray(coords, dtype=np.float64)
    path.parent.mkdir(parents=True, exist_ok=True)
    a.write_h5ad(path)


def write_file(path: Path, X, genes, coords=None) -> None:
    """An upload file as vec.io.make_submission writes it: dense float32 X, float32 coordinates in obsm."""
    a = ad.AnnData(X=np.asarray(X, dtype=np.float32))
    a.var_names = [str(g) for g in genes]
    a.obs_names = [f"cell_{i}" for i in range(a.n_obs)]
    if coords is not None:
        a.obsm[COORD_KEY] = np.asarray(coords, dtype=np.float32)
    path.parent.mkdir(parents=True, exist_ok=True)
    a.write_h5ad(path)


def sha(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_sums(archive: Path, jobs) -> None:
    (archive / "SHA256SUMS").write_text("".join(f"{sha(archive / j.file)}  {j.file}\n" for j in jobs))


def make_world(root: Path) -> SimpleNamespace:
    """A correct archive: one upload file per board kind, each built from its released stage the way the real builder
    builds it (T1 and the embryo keep exact released rows, the heart moves xy only, T3 keeps WT positions)."""
    rng = np.random.default_rng(20261109)
    data, archive = root / "data", root / "archive"
    jobs = [ig.Job("builds/T1_s1.h5ad", BOARDS["t1"], "t1", entry="M1"),
            ig.Job("builds/T2heart_s1.h5ad", BOARDS["heart"], "heart", entry="B0 s1 (H1)", hypothesis="H1"),
            ig.Job("builds/T2embryo_s1.h5ad", BOARDS["embryo"], "embryo", entry="fgw"),
            ig.Job("t3_build/allowed_s1.h5ad", BOARDS["t3"], "t3", branch="allowed", submission="s1")]
    panels = {BOARDS["t1"]: T1_GENES, BOARDS["heart"]: PANEL, BOARDS["embryo"]: PANEL, BOARDS["t3"]: PANEL}

    # T1: E10.5 cells re-drawn, a share moved in count space (log1p CP10k of the new counts), all library 1e4
    labels = rng.integers(0, 3, SOURCE_CELLS["t1"])
    means = np.exp(rng.normal(0.0, 1.0, (3, len(T1_GENES))))
    counts_t1 = rng.poisson(means[labels]).astype(np.int32)
    write_stage(data / "T1" / "E10.5_RNA.h5ad", counts_t1, None, T1_GENES)
    rows = rng.choice(SOURCE_CELLS["t1"], FILE_CELLS["t1"], replace=False)
    X_t1 = log1p_cp10k(counts_t1).toarray()[rows].astype(np.float32)
    moved = rng.random(FILE_CELLS["t1"]) < 0.3
    X_t1[moved] = log1p_cp10k(rng.binomial(counts_t1[rows[moved]], 0.75)).toarray()
    write_file(archive / "builds" / "T1_s1.h5ad", X_t1, T1_GENES)

    # Heart: E10.5 cells with a small momentum in count space, x and y grown about the centroid, z held
    field = _field(len(PANEL), rng)
    cloud = _cloud(SOURCE_CELLS["heart"], rng)
    counts_h = _spatial_counts(cloud, field, rng)
    write_stage(data / "T2" / "E10.5.h5ad", counts_h, cloud, PANEL)
    rows = rng.choice(SOURCE_CELLS["heart"], FILE_CELLS["heart"], replace=False)
    X_h = log1p_cp10k(counts_h[rows] + rng.poisson(0.3, (len(rows), len(PANEL)))).toarray()
    coords_h = cloud[rows].copy()
    mu = coords_h.mean(0)
    coords_h[:, :2] = mu[:2] + (coords_h[:, :2] - mu[:2]) * 1.8
    target_heart = t2_build.rms_radius(coords_h)
    write_file(archive / "builds" / "T2heart_s1.h5ad", X_h, PANEL, coords_h)

    # Embryo: a mixture of whole cells of E7.5 and E8.0 (exact rows), cells moved a little, then scaled to the target
    field = _field(len(PANEL), rng)
    c75, c80 = _cloud(SOURCE_CELLS["e75"], rng), _cloud(SOURCE_CELLS["e80"], rng)
    counts_75, counts_80 = _spatial_counts(c75, field, rng), _spatial_counts(c80, field, rng)
    write_stage(data / "T2" / "E7.5.h5ad", counts_75, c75, PANEL)
    write_stage(data / "T2" / "E8.0.h5ad", counts_80, c80, PANEL)
    n_b = FILE_CELLS["embryo"] // 2
    rows_b = rng.choice(SOURCE_CELLS["e75"], n_b, replace=False)
    rows_a = rng.choice(SOURCE_CELLS["e80"], FILE_CELLS["embryo"] - n_b, replace=False)
    X_e = np.vstack([log1p_cp10k(counts_75).toarray()[rows_b], log1p_cp10k(counts_80).toarray()[rows_a]])
    coords_e = np.vstack([c75[rows_b], c80[rows_a]]) * 0.9 + rng.normal(0, 0.2, (FILE_CELLS["embryo"], 3))
    perm = rng.permutation(FILE_CELLS["embryo"])
    X_e, coords_e = X_e[perm], coords_e[perm]
    target_embryo = t2_build.rms_radius(coords_e)
    write_file(archive / "builds" / "T2embryo_s1.h5ad", X_e, PANEL, coords_e)

    # T3: WT cells thinned in count space (a factor / depth step), each keeping its WT cell's exact position
    field = _field(len(PANEL), rng)
    cloud = _cloud(SOURCE_CELLS["t3"], rng)
    counts_t = _spatial_counts(cloud, field, rng)
    write_stage(data / "T2" / "E8.75.h5ad", counts_t, cloud, PANEL)
    rows = rng.choice(SOURCE_CELLS["t3"], FILE_CELLS["t3"], replace=False)
    X_t = log1p_cp10k(rng.binomial(counts_t[rows], 0.7)).toarray()
    write_file(archive / "t3_build" / "allowed_s1.h5ad", X_t, PANEL, cloud[rows])

    # The inputs the pass reads beside them
    sizes = {"heart": root / "size_heart.json", "embryo": root / "size_embryo.json"}
    sizes["heart"].write_text(json.dumps({
        "board": "heart", "valid": True,
        "preregistered": {"z_mode": "xy_grow_z_hold", "valid": True,
                          "scope": {"H1": {"pred_rms_full": target_heart},
                                    "H2": {"pred_rms_full": 0.8 * target_heart}}}}))
    sizes["embryo"].write_text(json.dumps({"board": "embryo", "valid": True,
                                           "preregistered": {"z_mode": "raw", "pred_rms_full": target_embryo}}))
    decisions = root / "decisions.json"
    decisions.write_text("{}\n")

    # Build records (fields by name, as the builders write them) and the archive's hashes
    def record(job: ig.Job, extra: dict) -> None:
        (archive / job.record).write_text(json.dumps({"entry": job.entry, "board": job.board, "dry_run": False,
                                                      "seed": ig.BUILD_SEED, "n_cells": "max_cells", **extra}))
    record(jobs[0], {"hold_unbracketed": True, "herd": False, "feng_hygiene": {"sha256": UNION}})
    record(jobs[1], {"size": {"hypothesis": "H1"}})
    record(jobs[2], {})
    (archive / jobs[3].record).write_text(json.dumps({
        "branch": "allowed", "submission": "s1",
        "inputs": {"board": {"key": BOARDS["t3"]}, "n_cells": "max_cells", "seed": ig.BUILD_SEED,
                   "wt": {"sha256": sha(data / "T2" / "E8.75.h5ad")}},
        "output": {"sha256": sha(archive / jobs[3].file), "n_cells": "max_cells", "check": {"ok": True}},
        "decision_file": {"sha256": sha(decisions)}, "hashes": {"mismatch": []}}))
    write_sums(archive, jobs)
    return SimpleNamespace(root=root, data=data, archive=archive, jobs=jobs, panels=panels, sizes=sizes,
                           decisions=decisions, target_heart=target_heart, target_embryo=target_embryo)


@pytest.fixture(scope="module")
def world(tmp_path_factory) -> SimpleNamespace:
    return make_world(tmp_path_factory.mktemp("world"))


@pytest.fixture
def damaged(world, tmp_path) -> SimpleNamespace:
    """A private copy of the archive that a test can damage (the released stages stay shared)."""
    archive = tmp_path / "archive"
    shutil.copytree(world.archive, archive)
    return SimpleNamespace(**{**vars(world), "archive": archive})


def run(w: SimpleNamespace) -> dict:
    return ig.run_pass(w.jobs, w.archive, w.data, w.panels, sizes=w.sizes, decisions=w.decisions, union_sha=UNION)


def flags_of(res: dict, rel: str, mark: str) -> list[str]:
    return [f for f in res["flags"] if f.startswith(f"{rel}: {mark}")]


def rewrite(w: SimpleNamespace, rel: str, change, genes=None) -> None:
    """Replace one upload file with change(X, coords) -> (X, coords) applied, and repair the hashes the archive and the
    T3 sidecar record, so that only the damage shows."""
    path = w.archive / rel
    a = ad.read_h5ad(path)
    coords = np.asarray(a.obsm[COORD_KEY], dtype=np.float64) if COORD_KEY in a.obsm else None
    X, coords = change(np.array(a.X, dtype=np.float32), coords)
    write_file(path, X, [str(g) for g in a.var_names] if genes is None else genes, coords)
    write_sums(w.archive, w.jobs)
    if rel.startswith("t3_build/"):
        side = path.with_suffix(".json")
        rec = json.loads(side.read_text())
        rec["output"]["sha256"] = sha(path)
        side.write_text(json.dumps(rec))


T1, HEART = "builds/T1_s1.h5ad", "builds/T2heart_s1.h5ad"
EMBRYO, T3 = "builds/T2embryo_s1.h5ad", "t3_build/allowed_s1.h5ad"


# ---------------------------------------------------------------------------------------------------------------------
# The whole pass
def test_a_correct_archive_raises_no_flag(world):
    res = run(world)
    assert res["verdict"] == "no flag", res["flags"]
    t1, heart, embryo, t3 = (res["files"][r] for r in (T1, HEART, EMBRYO, T3))
    assert t1["pb"]["pass"] and t1["x"]["pass"] and t1["library"]["pass"] and t1["provenance"]["pass"]
    assert 0.5 < t1["library"]["median_over_1e4"] < 1.5 and t1["provenance"]["exact_row_fraction"] > 0.5
    assert heart["moran"]["pass"] and heart["moran"]["agreement"] > ig.MORAN_FLOOR
    assert heart["size_as_preregistered"] is True
    assert heart["provenance"]["z_level_fraction"] == 1.0
    assert embryo["moran"]["pass"] and embryo["size_as_preregistered"] is True
    assert embryo["provenance"]["exact_row_fraction"] == 1.0
    assert 0.4 <= embryo["provenance"]["share_from_first_bracket"] <= 0.6
    assert t3["moran"]["pass"] and t3["provenance"]["exact_coord_fraction"] == 1.0
    assert res["source_rows_in_build_records"] is False
    for r in (T1, HEART, EMBRYO, T3):
        rec = res["files"][r]["record"]
        assert rec["present"] and rec["problems"] == [] and rec["missing"] == []


def test_the_output_holds_no_cell_count_and_no_field_that_can_carry_one(world):
    res = run(world)
    dumped = json.dumps(res) + ig.report_text(res, {"commit": "c" * 40, "tree_clean": True, "utc": "t"}, "archive")
    for name in FORBIDDEN_FIELDS:
        assert name not in dumped
    for n in (*FILE_CELLS.values(), *SOURCE_CELLS.values()):
        assert not re.search(rf"\b{n}\b", dumped), n


def test_a_column_order_slip_is_caught_by_the_pseudobulk_check(damaged):
    rewrite(damaged, T1, lambda X, c: (X[:, ::-1], c))
    res = run(damaged)
    assert flags_of(res, T1, "(a) pseudobulk")
    assert res["files"][T1]["pb"]["pearson"] < 0.9


def test_gene_names_out_of_order_are_flagged_even_with_the_right_columns(damaged):
    rewrite(damaged, T1, lambda X, c: (X, c), genes=list(reversed(T1_GENES)))
    res = run(damaged)
    assert flags_of(res, T1, "(a) the file's gene names")
    assert res["files"][T1]["panel_names_in_order"] is False


def test_x_and_coordinates_misaligned_are_caught_by_moran_not_by_pseudobulk(damaged):
    def shuffle(X, c):
        return X[np.random.default_rng(5).permutation(X.shape[0])], c
    rewrite(damaged, HEART, shuffle)
    res = run(damaged)
    assert not flags_of(res, HEART, "(a)")  # pseudobulk is invariant to the cell order
    assert flags_of(res, HEART, "(b) Moran")
    m = res["files"][HEART]["moran"]
    assert m["agreement"] < ig.MORAN_FLOOR and m["level_ratio"] < ig.MORAN_LEVEL_FLOOR


def test_a_partial_misalignment_is_caught_by_the_level_floor(damaged):
    """70% of the rows shuffled among themselves: every gene's I shrinks by about (1 - 0.7)^2, the Pearson of the
    profiles hardly moves; the level floor is what sees it."""
    def shuffle_part(X, c):
        rng = np.random.default_rng(6)
        idx = np.flatnonzero(rng.random(X.shape[0]) < 0.7)
        X = X.copy()
        X[idx] = X[rng.permutation(idx)]
        return X, c
    rewrite(damaged, T3, shuffle_part)
    res = run(damaged)
    m = res["files"][T3]["moran"]
    assert m["agreement"] > ig.MORAN_FLOOR and m["level_ratio"] < ig.MORAN_LEVEL_FLOOR
    assert flags_of(res, T3, "(b) Moran")


def test_moran_flag_when_a_gene_is_exactly_the_panel_but_positions_are_wrong_on_the_embryo(damaged):
    def wrong_positions(X, c):
        return X, np.random.default_rng(7).permutation(c)
    rewrite(damaged, EMBRYO, wrong_positions)
    res = run(damaged)
    assert flags_of(res, EMBRYO, "(b) Moran")


@pytest.mark.parametrize("bad, what", [(np.nan, "finite False"), (np.inf, "finite False"),
                                       (-1.0, "non-negative False")])
def test_non_finite_or_negative_x_is_flagged(damaged, bad, what):
    def poison(X, c):
        X = X.copy()
        X[3, 2] = bad
        return X, c
    rewrite(damaged, T1, poison)
    assert any(what in f for f in flags_of(run(damaged), T1, "(c) X"))


def test_a_t1_library_off_1e4_is_flagged(damaged):
    rewrite(damaged, T1, lambda X, c: (np.log1p(np.expm1(X) * 2.0).astype(np.float32), c))
    res = run(damaged)
    assert flags_of(res, T1, "(c) implied library")
    assert not res["files"][T1]["library"]["pass"]


def test_the_heart_size_is_asserted_and_its_value_never_written(damaged):
    rewrite(damaged, HEART, lambda X, c: (X, c * np.array([1.05, 1.05, 1.0])))
    res = run(damaged)
    assert res["files"][HEART]["size_as_preregistered"] is False
    assert flags_of(res, HEART, "(d)")
    dumped = json.dumps(res) + ig.report_text(res, {"commit": "c", "tree_clean": True, "utc": "t"}, "archive")
    for form in (repr(damaged.target_heart), f"{damaged.target_heart:.3f}", f"{damaged.target_heart:.1f}"):
        assert form not in dumped


def test_the_heart_size_check_needs_the_files_own_hypothesis(damaged):
    size = json.loads(damaged.sizes["heart"].read_text())
    size["preregistered"]["scope"]["H1"]["pred_rms_full"], size["preregistered"]["scope"]["H2"]["pred_rms_full"] = (
        size["preregistered"]["scope"]["H2"]["pred_rms_full"], size["preregistered"]["scope"]["H1"]["pred_rms_full"])
    path = damaged.root / "size_heart_swapped.json"
    path.write_text(json.dumps(size))
    damaged.sizes = {**damaged.sizes, "heart": path}
    assert flags_of(run(damaged), HEART, "(d)")


def test_the_embryo_size_is_asserted_on_its_raw_mode_forecast(damaged):
    rewrite(damaged, EMBRYO, lambda X, c: (X, c * 1.2))
    res = run(damaged)
    assert res["files"][EMBRYO]["size_as_preregistered"] is False
    assert flags_of(res, EMBRYO, "(d)")


def test_a_hash_that_differs_from_the_archive_is_flagged(damaged):
    sums = (damaged.archive / "SHA256SUMS").read_text()
    (damaged.archive / "SHA256SUMS").write_text(sums.replace(sha(damaged.archive / T1), "0" * 64))
    assert flags_of(run(damaged), T1, "(e) sha256")


def test_a_build_record_field_that_differs_is_flagged_by_name_only(damaged):
    path = damaged.archive / "builds" / "T1_s1.build.json"
    rec = json.loads(path.read_text())
    rec["seed"] = 123456789
    path.write_text(json.dumps(rec))
    res = run(damaged)
    f = flags_of(res, T1, "(e) build record fields")
    assert f and "seed" in f[0] and "123456789" not in json.dumps(res)


def test_a_build_record_that_keeps_source_rows_stops_the_pass(damaged):
    path = damaged.archive / "builds" / "T2embryo_s1.build.json"
    rec = json.loads(path.read_text())
    rec["source_rows"] = list(range(500))
    path.write_text(json.dumps(rec))
    res = run(damaged)
    assert res["source_rows_in_build_records"] is True
    assert flags_of(res, EMBRYO, "(e) a build record keeps source rows")


def test_row_level_provenance_flags_a_wrong_stage(damaged):
    def other_stage(X, c):  # another population: no row of it is a released E10.5 row
        counts = np.random.default_rng(9).poisson(1.5, X.shape).astype(np.int32) + 1
        return log1p_cp10k(counts).toarray().astype(np.float32), c
    rewrite(damaged, T1, other_stage)
    res = run(damaged)
    assert res["files"][T1]["provenance"]["exact_row_fraction"] == 0.0
    assert flags_of(res, T1, "(e) row-level provenance")


def test_the_embryo_mixture_must_be_about_half_each_bracket(world, damaged):
    """All rows from one bracket: every row is still an exact released row, but the share of the mixture is wrong."""
    stage = ad.read_h5ad(world.data / "T2" / "E7.5.h5ad")
    rows = np.random.default_rng(8).choice(stage.n_obs, FILE_CELLS["embryo"], replace=True)
    rewrite(damaged, EMBRYO, lambda X, c: (stage.X[rows].toarray().astype(np.float32), c))
    res = run(damaged)
    p = res["files"][EMBRYO]["provenance"]
    assert p["exact_row_fraction"] == 1.0 and p["share_from_first_bracket"] == 1.0
    assert flags_of(res, EMBRYO, "(e) row-level provenance")


def test_t3_coordinates_not_from_the_wt_stage_are_flagged(damaged):
    rewrite(damaged, T3, lambda X, c: (X, c + 0.5))
    res = run(damaged)
    assert res["files"][T3]["provenance"]["exact_coord_fraction"] == 0.0
    assert flags_of(res, T3, "(e) row-level provenance")


def test_heart_z_not_on_the_released_levels_is_flagged(damaged):
    rewrite(damaged, HEART, lambda X, c: (X, c + np.array([0.0, 0.0, 0.123])))
    res = run(damaged)
    assert res["files"][HEART]["provenance"]["z_level_fraction"] == 0.0
    assert flags_of(res, HEART, "(e) row-level provenance")


# ---------------------------------------------------------------------------------------------------------------------
# The run: stamps, refusals, the diagnostics.py entry
def test_main_writes_stamped_outputs_and_refuses_a_rerun(world, tmp_path, monkeypatch):
    monkeypatch.setattr(ig, "tree_state", lambda exempt: {"commit": "a" * 40, "clean": True})
    monkeypatch.setattr(ig.hygiene, "load_config", lambda *a, **k: {"union_sha256": UNION})
    out = tmp_path / "out"
    argv = ["--out", str(out), "--archive", str(world.archive), "--data", str(world.data)]
    ig.main(argv, jobs=world.jobs, panels=world.panels, sizes=world.sizes, decisions=world.decisions)
    body = json.loads((out / "integrity.json").read_text())
    text = (out / "integrity.txt").read_text()
    assert body["verdict"] == "no flag" and "VERDICT: no flag" in text
    assert body["provenance"]["commit"] == "a" * 40 and body["provenance"]["tree_clean"] is True
    assert body["provenance"]["utc"].endswith("Z") and body["inputs"]["files_sha256"][T1]
    assert body["inputs"]["sources_sha256"]["T2/E8.75.h5ad"] == sha(world.data / "T2" / "E8.75.h5ad")
    for name in ("integrity.json", "integrity.txt"):
        assert provenance.read_stamp(out / f"{name}.stamp.json")["sha256"] == sha(out / name)
    dumped = (out / "integrity.json").read_text() + text
    for n in (*FILE_CELLS.values(), *SOURCE_CELLS.values()):
        assert not re.search(rf"\b{n}\b", dumped), n
    for name in FORBIDDEN_FIELDS:
        assert name not in dumped
    with pytest.raises(SystemExit, match="already written"):
        ig.main(argv, jobs=world.jobs, panels=world.panels, sizes=world.sizes, decisions=world.decisions)


@pytest.mark.parametrize("state", [{"commit": "a" * 40, "clean": False}, {"commit": "unknown", "clean": None}])
def test_main_refuses_a_tree_that_is_not_clean(world, tmp_path, monkeypatch, state):
    monkeypatch.setattr(ig, "tree_state", lambda exempt: state)
    out = tmp_path / "out"
    with pytest.raises(SystemExit, match="clean tree"):
        ig.main(["--out", str(out), "--archive", str(world.archive), "--data", str(world.data)],
                jobs=world.jobs, panels=world.panels, sizes=world.sizes, decisions=world.decisions)
    assert not out.exists()


def test_diagnostics_integrity_runs_the_pass_and_touches_nothing_else(monkeypatch, tmp_path):
    import diagnostics as dg

    seen = []
    monkeypatch.setattr(ig, "main", lambda argv=None, **kw: seen.append(argv))
    out = tmp_path / "o"
    dg.main(["integrity", "--out", str(out)])
    assert seen == [["--out", str(out)]]
    assert not out.exists()  # the command itself writes nothing: the pass refuses and writes on its own


# ---------------------------------------------------------------------------------------------------------------------
# The parts
def test_the_ten_jobs_are_the_frozen_set_and_their_records():
    assert [j.file for j in ig.JOBS] == ["builds/T1_s1.h5ad", "builds/T1_s2.h5ad", "builds/T2heart_s1.h5ad",
                                         "builds/T2heart_s2.h5ad", "builds/T2embryo_s1.h5ad", "builds/T2embryo_s2.h5ad",
                                         "t3_build/allowed_s1.h5ad", "t3_build/refused_s1.h5ad",
                                         "t3_build/barred_s1.h5ad", "t3_build/barred_s2.h5ad"]
    assert ig.JOBS[0].record == "builds/T1_s1.build.json" and ig.JOBS[6].record == "t3_build/allowed_s1.json"
    assert ig.JOBS[1].cm_share == 0.60
    assert [j.hypothesis for j in ig.JOBS if j.kind == "heart"] == ["H1", "H2"]
    assert [j.entry for j in ig.JOBS if j.kind == "heart"] == ["B0 s1 (H1)", "B0 s2_trim (H2)"]
    assert [j.entry for j in ig.JOBS if j.kind == "embryo"] == ["fgw", "fgw_shape0"]
    assert sorted({j.kind for j in ig.JOBS}) == ["embryo", "heart", "t1", "t3"]


def test_the_record_expectations_name_the_builders_fields():
    t1s2 = ig.JOBS[1]
    exp = ig.record_expectations(t1s2, "f", {}, "d", "u")
    assert exp["entry"] == "M1" and exp["cm_share_of_heart"] == 0.60 and exp["feng_hygiene.sha256"] == "u"
    assert exp["seed"] == ig.BUILD_SEED == 1_000_001 and exp["hold_unbracketed"] is True and exp["herd"] is False
    heart = ig.record_expectations(ig.JOBS[3], "f", {}, "d", "u")
    assert heart["size.hypothesis"] == "H2" and "cm_share_of_heart" not in heart and "feng_hygiene.sha256" not in heart
    t3 = ig.record_expectations(ig.JOBS[9], "f", {"T2/E8.75.h5ad": "wt"}, "d", "u")
    assert t3["branch"] == "barred" and t3["submission"] == "s2" and t3["inputs.wt.sha256"] == "wt"
    assert t3["output.sha256"] == "f" and t3["decision_file.sha256"] == "d" and t3["hashes.mismatch"] == []
    assert t3["inputs.seed"] == ig.BUILD_SEED and t3["output.n_cells"] == t3["inputs.n_cells"] == "max_cells"


def test_pseudobulk_is_the_mean_of_x_for_sparse_and_dense():
    X = np.random.default_rng(0).poisson(1.5, (7001, 5)).astype(np.float32)
    np.testing.assert_allclose(ig.pseudobulk(X), X.mean(0))
    np.testing.assert_allclose(ig.pseudobulk(sparse.csr_matrix(X)), X.mean(0))


def test_pb_agreement_needs_both_correlations_and_fails_on_nan():
    a = np.arange(50, dtype=float)
    assert ig.pb_agreement(a, a * 2 + 1, 0.95)["pass"]
    assert not ig.pb_agreement(a, np.r_[a[:25], a[:25][::-1]], 0.95)["pass"]
    assert not ig.pb_agreement(a, np.full(50, 1.0), 0.95)["pass"]  # constant: undefined, a fail


def test_moran_profile_equals_the_formula_by_brute_force():
    """I_g = sum_i z_ig mean_(j in kNN(i)) z_jg / sum_i z_ig^2 (the recipe of veckit's morans_I_profile at the pinned
    git), evaluated cell by cell."""
    rng = np.random.default_rng(4)
    coords = rng.uniform(0, 10, (60, 3))
    X = np.c_[coords[:, 0] + rng.normal(0, 0.1, 60), rng.normal(0, 1, 60), np.full(60, 2.0)].astype(np.float32)
    idx = ig.knn_graph(coords, k=5)
    assert idx.shape == (60, 5) and not (idx == np.arange(60)[:, None]).any()
    Z = X.astype(np.float64) - X.astype(np.float64).mean(0)
    want = [sum(Z[i, g] * Z[idx[i], g].mean() for i in range(60)) / ((Z[:, g] ** 2).sum() + 1e-12) for g in range(3)]
    np.testing.assert_allclose(ig.moran_profile(X, idx, np.arange(3)), want)
    assert ig.moran_profile(X, idx, np.arange(3))[0] > 0.5  # a smooth gene is autocorrelated
    assert abs(ig.moran_profile(X, idx, np.arange(3))[2]) < 1e-6  # a constant gene is not


def test_top_variance_genes_are_chosen_on_the_source_as_the_scorer_does():
    X = np.random.default_rng(1).normal(0, 1, (200, 8)) * np.array([1, 5, 2, 9, 3, 8, 4, 7])
    assert list(ig.top_variance_genes(X, 3)) == [3, 5, 7]
    np.testing.assert_array_equal(ig.top_variance_genes(sparse.csr_matrix(X), 3), ig.top_variance_genes(X, 3))


def test_moran_check_reports_the_control_and_no_ratio_without_structure():
    src = np.linspace(0.1, 0.9, 40)
    r = ig.moran_check(src * 0.8, src, np.random.default_rng(2).normal(0, 0.01, 40))
    assert r["agreement"] > 0.99 and abs(r["shuffled_control_agreement"]) < 0.5
    assert r["level_ratio"] == pytest.approx(0.8)
    assert np.isnan(ig.moran_check(src, -src, src)["level_ratio"])  # a source without structure has no level to keep


def test_x_sanity_and_the_library_band():
    X = np.log1p(np.array([[1e4, 0.0], [9e3, 1e3]])).astype(np.float32)
    assert ig.x_sanity(X)["pass"] and ig.implied_library(X)["pass"]
    assert not ig.x_sanity(np.array([[np.nan, 1.0]]))["finite"]
    assert not ig.x_sanity(np.array([[-0.5, 1.0]]))["non_negative"]
    assert not ig.implied_library(X + 1.0)["pass"]
    assert not ig.implied_library(np.array([[np.nan, 1.0], [1.0, 1.0]], dtype=np.float32))["pass"]


def test_exact_rows_are_found_by_digest_and_never_stored():
    X = np.random.default_rng(3).random((20, 6)).astype(np.float32)
    digests = ig.row_digests(sparse.csr_matrix(X))
    assert digests == ig.row_digests(X) and all(len(d) == 16 for d in digests)
    other = X[:5] + 1e-3
    np.testing.assert_array_equal(ig.exact_row_mask(np.vstack([X[:3], other]), digests), [True] * 3 + [False] * 5)


def test_coordinate_keys_match_exactly():
    C = np.random.default_rng(4).random((30, 3)).astype(np.float32)
    keys = ig.coord_keys(C)
    assert ig.in_keys(C[:4], keys).all() and not ig.in_keys(C[:4] + 1e-3, keys).any()


def test_size_matches_is_assert_only_and_refuses_a_bad_json(tmp_path):
    coords = np.random.default_rng(0).random((50, 3)) * 10
    rms = t2_build.rms_radius(coords)
    heart = {"valid": True, "preregistered": {"z_mode": "xy_grow_z_hold", "scope": {"H2": {"pred_rms_full": rms}}}}
    embryo = {"valid": True, "preregistered": {"z_mode": "raw", "pred_rms_full": rms}}
    assert ig.size_matches(coords, "heart", heart, "H2", 1e-6) is True
    assert ig.size_matches(coords * 1.01, "heart", heart, "H2", 1e-6) is False
    assert ig.size_matches(coords, "heart", heart, "H1", 1e-6) is False           # another hypothesis: no forecast
    assert ig.size_matches(coords, "embryo", embryo, None, 1e-6) is True
    assert ig.size_matches(coords, "embryo", {**embryo, "valid": False}, None, 1e-6) is False
    assert ig.size_matches(coords, "embryo", heart, None, 1e-6) is False          # the wrong JSON: no raw forecast
    path = tmp_path / "s.json"
    path.write_text(json.dumps(embryo))
    assert ig.size_matches(coords, "embryo", path, None, 1e-6) is True


def test_source_row_fields_finds_per_cell_lists_only():
    rec = {"a": {"rows": list(range(12))}, "b": [{"partners": ["x"] * 10}], "rows": [{"n": 1}] * 20,
           "short": {"rows": [1, 2]}, "names": {"genes": ["a"] * 50}}
    assert sorted(ig.source_row_fields(rec)) == ["a.rows", "b[0].partners"]


def test_check_record_separates_missing_from_differing_and_keeps_values_out(tmp_path):
    job = ig.Job("x/f.h5ad", "B", "heart", entry="E", hypothesis="H1")
    (tmp_path / "x").mkdir()
    rec = {"entry": "other", "seed": 5, "size": {"hypothesis": "H1"}}
    (tmp_path / "x" / "f.build.json").write_text(json.dumps(rec))
    out = ig.check_record(job, tmp_path, {"entry": "E", "seed": 5, "n_cells": "max_cells", "size.hypothesis": "H1"})
    assert out["present"] and out["problems"] == ["entry"] and out["missing"] == ["n_cells"]
    assert "other" not in json.dumps(out)
    assert not ig.check_record(ig.Job("y/none.h5ad", "B", "heart"), tmp_path, {"entry": "E"})["present"]


def test_equal_compares_floats_within_tolerance_only():
    assert ig._equal(0.6000000000000001, 0.6) and not ig._equal(0.6001, 0.6)
    assert ig._equal("max_cells", "max_cells") and not ig._equal(True, False) and ig._equal([], [])


def test_check_stamp_reads_a_real_stamp(tmp_path):
    f = tmp_path / "file.bin"
    f.write_bytes(b"content")
    assert ig.check_stamp(f, sha(f)) == {"present": False, "sha256_matches": None, "commit": None, "clean": None}
    provenance.write_stamp(f)
    s = ig.check_stamp(f, sha(f))
    assert s["present"] and s["sha256_matches"] is True
    assert ig.check_stamp(f, "0" * 64)["sha256_matches"] is False


def test_read_sums_parses_both_sha256sum_formats(tmp_path):
    p = tmp_path / "SHA256SUMS"
    p.write_text(f"{'a' * 64}  builds/T1_s1.h5ad\n{'b' * 64} *t3_build/allowed_s1.h5ad\n\n")
    assert ig.read_sums(p) == {"builds/T1_s1.h5ad": "a" * 64, "t3_build/allowed_s1.h5ad": "b" * 64}
