"""Synthetic checks of the external-data window filters (no competition or external data).

Run: uv run pytest experiments/external -q
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parent))
import feng2022 as feng  # noqa: E402
import filterio as fio  # noqa: E402
import jax_qiu2024 as jax  # noqa: E402
import windows as w  # noqa: E402


# ------------------------------------------------------------------ windows and stage labels
@pytest.mark.parametrize("label,eday", [("E9.5", 9.5), ("E14.333", 14.333), ("E18", 18.0), ("P0", 19.5),
                                        ("p9", 28.5), (" E13.75 ", 13.75)])
def test_parse_stage(label, eday):
    assert w.parse_stage(label) == pytest.approx(eday)


@pytest.mark.parametrize("label", ["", "E", "9.5", "embryo_3", "E9.5.1", "Theiler 15", "nan"])
def test_parse_stage_rejects_unknown_labels(label):
    with pytest.raises(ValueError):
        w.parse_stage(label)


@pytest.mark.parametrize("task,eday,inside", [
    ("T1", 9.5, False), ("T1", 9.51, True), ("T1", 12.5, True), ("T1", 13.5, True), ("T1", 13.75, False),
    ("T1", 8.5, False), ("T2-heart", 8.25, False), ("T2-heart", 8.5, True), ("T2-heart", 8.75, False),
    ("T2-heart", 9.5, False), ("T2-heart", 10.5, True), ("T2-heart", 14.5, False), ("T2-embryo", 7.25, False),
    ("T2-embryo", 7.5, True), ("T2-embryo", 7.75, True), ("T2-embryo", 8.0, False),
])
def test_window_boundaries(task, eday, inside):
    assert w.in_window(eday, task) is inside


def test_union_of_windows_and_task_checks():
    assert w.allowed(8.5, ["T1"]) and not w.allowed(8.5, ["T1", "T2-heart"])
    with pytest.raises(NotImplementedError):
        w.allowed(9.0, ["T3"])
    with pytest.raises(ValueError):
        w.allowed(9.0, ["T4"])
    with pytest.raises(ValueError):
        w.allowed(9.0, [])


# ------------------------------------------------------------------ Feng 2022
@pytest.mark.parametrize("stage,keep", [("E9.5", True), ("E9.0", False), ("E10.5", False), ("E11.5", False),
                                        ("E12.5", False), ("E13.5", False), ("E14.0", False), ("E14.5", True),
                                        ("E18", True), ("P0", True), ("P9", True)])
def test_feng_stage_rule(stage, keep):
    assert feng.keep_stage(stage)[0] is keep
    assert feng.decide((stage, "CD1"))[0] is keep
    assert feng.decide((stage, "C57BL/6"))[0] is False


def test_feng_unknown_stage_aborts():
    with pytest.raises(ValueError):
        feng.decide(("adult", "CD1"))


def write_csv(path: Path, header: list[str], rows: list[list[str]]) -> None:
    with open(path, "w", newline="") as f:
        wr = csv.writer(f, lineterminator="\n")
        wr.writerow(header)
        wr.writerows(rows)


def read_csv(path: Path) -> list[list[str]]:
    with open(path, newline="") as f:
        return list(csv.reader(f))


FENG_HEADER = ["", "orig.ident", "strain", "stage", "zone", "cell_type_spec"]
FENG_ROWS = [["c0", "a", "CD1", "E9.5", "LV", "ventricular_cm"], ["c1", "a", "CD1", "E10.5", "LV", "ventricular_cm"],
             ["c2", "a", "CD1", "E14.5", "RA", "atrial_cm"], ["c3", "a", "CD1", "E13.5", "LA", "atrial_cm"],
             ["c4", "a", "B6", "E14.5", "RV", "ventricular_cm"], ["c5", "a", "CD1", "P0", "LV", "endocardial_ec"],
             ["c6", "a", "CD1", "E12.5", "LV", "epicardial, outer"]]


def test_feng_obs_and_counts_filters(tmp_path):
    src = tmp_path / "obs.csv"
    write_csv(src, FENG_HEADER, FENG_ROWS)
    rep = feng.filter_obs(src, tmp_path / "out" / "obs.csv")
    out = read_csv(tmp_path / "out" / "obs.csv")
    assert out[0] == FENG_HEADER  # header (including the unnamed index column) preserved
    assert [r[0] for r in out[1:]] == ["c0", "c2", "c5"]
    assert all(not w.in_window(w.parse_stage(r[3]), "T1") for r in out[1:])
    assert rep["rows_in"] == 7 and rep["rows_kept"] == 3
    assert rep["positions"].tolist() == [0, 2, 5]
    X = sparse.csr_matrix(np.arange(7 * 4).reshape(7, 4))
    sparse.save_npz(tmp_path / "counts.npz", X)
    feng.filter_counts(tmp_path / "counts.npz", tmp_path / "out" / "counts.npz", rep["positions"], 7)
    Y = sparse.load_npz(tmp_path / "out" / "counts.npz").toarray()
    assert np.array_equal(Y, X.toarray()[[0, 2, 5]])
    with pytest.raises(ValueError):  # rows must align with the obs table
        feng.filter_counts(tmp_path / "counts.npz", tmp_path / "out" / "c2.npz", rep["positions"], 8)


def test_failed_filter_leaves_no_output(tmp_path):
    src = tmp_path / "obs.csv"
    write_csv(src, FENG_HEADER, [*FENG_ROWS, ["c7", "a", "CD1", "late fetal", "LV", "x"]])
    with pytest.raises(ValueError):
        feng.filter_obs(src, tmp_path / "out.csv")
    assert not (tmp_path / "out.csv").exists() and not list(tmp_path.glob(".*partial*"))
    write_csv(src, FENG_HEADER, [*FENG_ROWS, ["c8", "a", "CD1"]])  # malformed row
    with pytest.raises(ValueError):
        feng.filter_obs(src, tmp_path / "out.csv")
    assert not (tmp_path / "out.csv").exists()


def make_stats(path: Path, keys: list[str], extra: dict | None = None) -> None:
    k = len(keys)
    np.savez_compressed(path, genes=np.array(["g1", "g2"], dtype=object), keys=np.array(keys, dtype=object),
                        n=np.arange(1, k + 1, dtype=float), S_log=np.ones((k, 2)), S_cp=np.ones((k, 2)),
                        S_nz=np.ones((k, 2)), depth_sum=np.ones(k), **(extra or {}))


def test_feng_stats_filter(tmp_path):
    keys = ["E9.5|LV|ventricular_cm|A", "E10.5|LV|ventricular_cm|A", "E13.5|RA|atrial_cm|B",
            "E14.5|RA|atrial_cm|B", "P9|LV|blood|A"]
    make_stats(tmp_path / "s.npz", keys)
    rep = feng.filter_stats(tmp_path / "s.npz", tmp_path / "f.npz")
    z = np.load(tmp_path / "f.npz")  # written without object arrays
    assert z["keys"].tolist() == ["E9.5|LV|ventricular_cm|A", "E14.5|RA|atrial_cm|B", "P9|LV|blood|A"]
    assert z["n"].tolist() == [1.0, 4.0, 5.0] and z["S_log"].shape == (3, 2)
    assert rep["cells_kept"] == 10 and rep["cells_dropped"] == 5
    make_stats(tmp_path / "bad.npz", keys, extra={"stage": np.array("mixed")})
    with pytest.raises(ValueError):
        feng.filter_stats(tmp_path / "bad.npz", tmp_path / "g.npz")


# ------------------------------------------------------------------ JAX / Qiu 2024
@pytest.mark.parametrize("n,keep_t1,keep_t2", [(1, True, False), (12, True, False), (13, True, True),
                                               (22, True, True), (23, False, False), (24, False, False),
                                               (25, False, False), (49, False, False), (50, False, False),
                                               (53, False, False), (54, True, True), (74, True, True),
                                               (75, False, False), (90, False, False)])
def test_jax_embryo_rule(n, keep_t1, keep_t2):
    assert jax.keep_embryo(n, ("T1",))[0] is keep_t1
    assert jax.keep_embryo(n, ("T1", "T2-heart"))[0] is keep_t2


def test_jax_rows_and_day_check(tmp_path):
    header = ["", "cell_id", "day", "embryo_id", "celltype_update"]
    rows = [["0", "a", "E8.5", "embryo_3", "First heart field"], ["1", "b", "E9.75", "embryo_23", "x"],
            ["2", "c", "E11.5", "embryo_37", "x"], ["3", "d", "E13.75", "embryo_50", "x"],
            ["4", "e", "E14.25", "embryo_54", "Ventricular cardiomyocytes"], ["5", "f", "E9.5", "embryo_22", "y"]]
    src = tmp_path / "df_cell.csv"
    write_csv(src, header, rows)
    rep = jax.filter_cells(src, tmp_path / "t1.csv")
    assert [r[1] for r in read_csv(tmp_path / "t1.csv")[1:]] == ["a", "e", "f"]
    assert rep["rows_dropped"] == 3
    jax.filter_cells(src, tmp_path / "t2.csv", tasks=("T2-heart",))
    assert [r[1] for r in read_csv(tmp_path / "t2.csv")[1:]] == ["e", "f"]  # embryo 3 is T1-only
    write_csv(src, header, [["0", "a", "E10.5", "embryo_22", "x"]])  # kept embryo with an in-window day label
    with pytest.raises(AssertionError):
        jax.filter_cells(src, tmp_path / "t3.csv")
    assert not (tmp_path / "t3.csv").exists()
    with pytest.raises(ValueError):
        jax.embryo_number("E8.5_embryo")


def test_jax_file_gate(tmp_path):
    src_dir, dst_dir = tmp_path / "src", tmp_path / "dst"
    src_dir.mkdir()
    for n in (1, 23, 30, 54):
        np.savez(src_dir / f"embryo_{n}.npz", embryo=np.array(f"embryo_{n}"))
    rep = jax.gate_files(src_dir / "embryo_{n}.npz", dst_dir / "stats" / "embryo_{n}.npz", check_npz=True)
    assert rep["kept"] == [1, 54] and rep["skipped_not_kept"] == [23, 30]
    assert sorted(p.name for p in (dst_dir / "stats").iterdir()) == ["embryo_1.npz", "embryo_54.npz"]
    assert rep["files"][1]["mode"] in ("hardlink", "copy")
    rep2 = jax.gate_files(src_dir / "embryo_{n}.npz", dst_dir / "t2" / "embryo_{n}.npz", tasks=("T1", "T2-heart"))
    assert rep2["kept"] == [54]
    np.savez(src_dir / "embryo_2.npz", embryo=np.array("embryo_9"))  # contents disagree with the name
    with pytest.raises(ValueError):
        jax.gate_files(src_dir / "embryo_{n}.npz", dst_dir / "x" / "embryo_{n}.npz", check_npz=True)


# ------------------------------------------------------------------ manifests
def test_verify_manifest_output(tmp_path):
    d = tmp_path / "ext"
    d.mkdir()
    f = d / "a.npz"
    np.savez(f, x=np.arange(3))
    fio.write_manifest(d, {"outputs": {"a.npz": {"sha256": fio.sha256_file(f)}}})
    assert fio.verify_manifest_output(f)["sha256"] == fio.sha256_file(f)
    np.savez(f, x=np.arange(4))
    with pytest.raises(ValueError):
        fio.verify_manifest_output(f)
    g = d / "b.npz"
    np.savez(g, x=np.arange(3))
    with pytest.raises(ValueError):
        fio.verify_manifest_output(g)
    json.loads((d / "MANIFEST.json").read_text())
    with pytest.raises(FileNotFoundError):
        fio.verify_manifest_output(tmp_path / "loose.npz")
