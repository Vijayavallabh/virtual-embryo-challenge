import json
import time

import anndata as ad
import numpy as np
import pytest
from scipy import sparse

from vec.spatial import find_duplicate_records, identical_rows, qc_stage, snap_z

STEP = 15.0


def _stage(per_level=(200,) * 6, n_genes=40, seed=0, origin=-37.5, lam=3.0):
    """Poisson counts, uniform xy, z on levels origin + STEP * k, rows in descending z like the released files."""
    rng = np.random.default_rng(seed)
    levels = origin + STEP * np.arange(len(per_level))
    z = np.repeat(levels[::-1], per_level[::-1])
    xy = rng.uniform(-300, 300, size=(z.size, 2))
    counts = rng.poisson(lam, size=(z.size, n_genes)).astype(np.float32)
    return counts, np.c_[xy, z]


def _adata(counts, coords, **obs):
    a = ad.AnnData(X=np.log1p(counts).astype(np.float32))
    a.layers["counts"] = sparse.csr_matrix(counts)
    a.obsm["spatial_3D"] = np.asarray(coords, dtype=np.float64)
    for k, v in obs.items():
        a.obs[k] = v
    return a


# --- snap_z -----------------------------------------------------------------------------------------------------------

def test_snap_z_absorbs_float_noise():
    rng = np.random.default_rng(0)
    levels = -280.256503 + STEP * np.arange(40)
    idx = rng.integers(0, 40, 5000)
    z = levels[idx] + rng.choice([0.0, 1.1e-13, -5.7e-14, 2.8e-14], size=idx.size)
    assert np.unique(z).size > 40
    zl = snap_z(z)
    assert zl.levels.size == 40
    assert zl.step == pytest.approx(STEP)
    np.testing.assert_array_equal(zl.index, idx)
    np.testing.assert_allclose(zl.levels, levels, atol=1e-12)
    np.testing.assert_allclose(zl.grid, np.arange(40), atol=1e-9)
    assert zl.gaps() == []


@pytest.mark.parametrize(("origin", "step"), [(0.5, 1.0), (7.5, 15.0), (-240.5, 14.0), (-0.5, 0.5)])
def test_snap_z_levels_on_half_values_stay_distinct(origin, step):
    levels = origin + step * np.arange(12)
    zl = snap_z(np.repeat(levels, 3)[::-1])
    assert zl.levels.size == 12
    np.testing.assert_array_equal(zl.levels, levels)  # single-valued levels come back exactly
    assert zl.step == pytest.approx(step)
    np.testing.assert_array_equal(zl.index, np.repeat(np.arange(12), 3)[::-1])


def test_snap_z_half_step_level_gets_a_half_grid_position():
    levels = np.sort(np.r_[STEP * np.arange(10), STEP * 4 + STEP / 2])  # one extra section between two levels
    zl = snap_z(np.repeat(levels, 5))
    assert zl.levels.size == 11
    assert zl.step == pytest.approx(STEP)
    assert np.isclose(zl.grid, 4.5).sum() == 1


def test_snap_z_gap_keeps_the_step_and_reports_missing_levels():
    levels = np.r_[-370 + 10.0 * np.arange(6), -290 + 10.0 * np.arange(30)]  # -320 then -290: two levels missing
    zl = snap_z(np.repeat(levels, 4))
    assert zl.levels.size == 36
    assert zl.step == pytest.approx(10.0)
    assert zl.gaps() == [{"from": -320.0, "to": -290.0, "steps": 3.0, "missing": 2}]
    assert zl.grid[6] - zl.grid[5] == pytest.approx(3.0)


def test_snap_z_merges_a_level_stored_twice_at_a_small_offset():
    levels = 10.0 * np.arange(20)
    z = np.r_[np.repeat(levels, 4), np.full(4, levels[7] + 0.001)]
    zl = snap_z(z)
    assert zl.levels.size == 20
    assert zl.step == pytest.approx(10.0, rel=1e-3)
    assert (zl.index[-4:] == 7).all()
    tight = snap_z(z, tol=1e-6)  # an explicit tolerance below the offset keeps the two copies apart
    assert tight.levels.size == 21
    assert tight.step == pytest.approx(10.0, rel=1e-3)


def test_snap_z_two_sections_on_one_level_and_degenerate_inputs():
    z = np.r_[np.repeat(STEP * np.arange(5), 10), np.full(10, STEP * 2)]  # level 2 holds two sections
    zl = snap_z(z)
    assert zl.levels.size == 5 and np.bincount(zl.index)[2] == 20
    one = snap_z(np.full(7, 3.25))
    assert one.levels.tolist() == [3.25] and np.isnan(one.step) and (one.index == 0).all()
    assert one.grid.tolist() == [0.0] and one.gaps() == []
    assert snap_z([]).levels.size == 0
    with pytest.raises(ValueError, match="NaN"):
        snap_z([0.0, np.nan])
    with pytest.raises(ValueError, match="tol"):
        snap_z([0.0, 1.0], tol=-1.0)


# --- duplicated records ---------------------------------------------------------------------------------------------

def _with_copies():
    """A 6-level stage plus copies covering every case; returns (adata, expected redundant rows, named rows)."""
    counts, coords = _stage()
    counts[7] = 0
    counts[7, :5] = 1  # library 5, below min_library
    extra_c, extra_p = [], []

    def add(src, dxy=(0.0, 0.0), dz=0.0, bump=None):
        row = counts[src].copy()
        if bump is not None:
            row[bump] += 1
        extra_c.append(row)
        extra_p.append(coords[src] + [dxy[0], dxy[1], dz])
        return counts.shape[0] + len(extra_c)  # final row index (one copy is prepended below)

    flagged = [
        add(0, dz=-STEP + 1e-13),          # adjacent level (with float noise in z)
        add(1, dxy=(0.001, -0.001)),       # same level, ~0.001 xy offset
        add(6, dz=-STEP),                  # copied into two neighbouring levels: both copies chain to the original
        add(6, dz=-2 * STEP),
    ]
    named = {
        "two_levels": add(2, dz=-2 * STEP),         # two levels away: flagged only with max_level_gap=2
        "one_count_off": add(3, dz=-STEP, bump=5),  # one count differs: never flagged
        "moved": add(4, dxy=(5.0, 0.0), dz=-STEP),  # identical row 5 units away: flagged only with a wide xy_tol
        "low_library": add(7, dz=-STEP),            # low library: reported, not flagged
    }
    first = (counts[10], coords[10] + [0, 0, -STEP])  # a copy stored before its original in file order
    C = np.vstack([first[0][None], counts, np.array(extra_c)])
    P = np.vstack([first[1][None], coords, np.array(extra_p)])
    return _adata(C, P), sorted([*flagged, 1 + 10]), named


def test_find_duplicate_records_flags_copies_and_keeps_the_first_in_file_order():
    a, expected, _ = _with_copies()
    mask, rep = find_duplicate_records(a)
    assert np.flatnonzero(mask).tolist() == expected
    assert not mask[0]  # the earlier copy is kept, its later original is flagged
    assert rep["n_redundant"] == 5 and rep["n_with_copies"] == 4
    assert rep["by_gap"] == {"0": 1, "1": 3, "2": 1}
    assert sum(p["n"] for p in rep["by_level_pair"]) == 5
    assert rep["below_min_library"] == 1
    assert rep["max_xy_offset"] == pytest.approx(0.001, abs=1e-9)
    assert rep["z_step"] == pytest.approx(STEP)
    json.dumps(rep, allow_nan=False)


def test_find_duplicate_records_parameters():
    a, expected, named = _with_copies()
    assert find_duplicate_records(a, max_level_gap=2)[1]["n_redundant"] == 6
    assert find_duplicate_records(a, max_level_gap=0)[1]["by_gap"] == {"0": 1}
    assert find_duplicate_records(a, min_library=1)[1]["n_redundant"] == 6
    assert find_duplicate_records(a, xy_tol=1e-4)[1]["n_redundant"] == 4  # the 0.001-offset copy is missed
    assert find_duplicate_records(a, xy_tol=10.0)[1]["n_redundant"] == 6  # the copy 5 units away now counts
    mask, rep = find_duplicate_records(a, xy_tol=np.inf)  # identity in the same or an adjacent level, xy ignored
    assert rep["n_redundant"] == 6 and rep["xy_tol"] == "any" and rep["max_xy_offset"] == pytest.approx(5.0)
    assert np.flatnonzero(mask).tolist() == sorted([*expected, named["moved"]])


def test_find_duplicate_records_no_false_positives_on_near_identical_rows():
    counts, coords = _stage(per_level=(300,) * 4, lam=0.5)
    n = counts.shape[0]
    near = counts[:200].copy()
    cols = np.random.default_rng(1).integers(0, counts.shape[1], 200)
    near[np.arange(200), cols] += 1  # every row differs from its source by a single count
    a = _adata(np.vstack([counts, near]), np.vstack([coords, coords[:200] - [0, 0, STEP]]))
    mask, rep = find_duplicate_records(a, min_library=1)
    assert mask.shape == (n + 200,) and not mask.any()
    assert rep["n_redundant"] == 0 and rep["by_level_pair"] == [] and rep["by_gap"] == {}


def test_find_duplicate_records_requires_counts_and_coordinates():
    counts, coords = _stage()
    a = _adata(counts, coords)
    del a.layers["counts"]
    with pytest.raises(KeyError, match="counts"):
        find_duplicate_records(a)
    b = _adata(counts, coords)
    del b.obsm["spatial_3D"]
    with pytest.raises(KeyError, match="spatial_3D"):
        find_duplicate_records(b)


def test_find_duplicate_records_scales_to_60k_cells():
    rng = np.random.default_rng(3)
    n, g, n_levels = 60_000, 100, 40
    counts = rng.poisson(0.4, size=(n, g)).astype(np.float32)
    z = np.repeat(STEP * np.arange(n_levels)[::-1], n // n_levels)
    coords = np.c_[rng.uniform(-500, 500, size=(n, 2)), z]
    src = rng.choice(np.flatnonzero((counts.sum(1) >= 20) & (z > 0)), 1000, replace=False)
    single = np.zeros((3000, g), np.float32)
    single[:, 0] = 1  # 3000 identical one-count rows: one big identical group
    single_pos = np.c_[rng.uniform(-500, 500, size=(3000, 2)), rng.choice(z, 3000)]
    single_pos[1:6] = single_pos[0] - [0, 0, 0]  # five of them at the same xy and level as the first
    C = np.vstack([counts, counts[src], single])
    P = np.vstack([coords, coords[src] - [0, 0, STEP], single_pos])
    t0 = time.perf_counter()
    mask, rep = find_duplicate_records(_adata(C, P))
    elapsed = time.perf_counter() - t0
    assert np.flatnonzero(mask).tolist() == list(range(n, n + 1000))
    assert rep["by_gap"] == {"1": 1000} and rep["below_min_library"] == 5
    assert elapsed < 60, elapsed  # generous bound for a loaded machine; typically about a second
    mask_any, rep_any = find_duplicate_records(_adata(C, P), xy_tol=np.inf)  # level sweep inside the big group
    assert np.array_equal(mask_any, mask) and rep_any["below_min_library"] == 2999


# --- identical rows and qc_stage --------------------------------------------------------------------------------------

def test_identical_rows_counts_repeats_anywhere():
    M = np.random.default_rng(4).poisson(2.0, size=(50, 30)).astype(np.float32)
    M[[10, 20, 30]] = M[5]
    M[40] = 0
    M[41, :] = 0
    M[41, 0] = 1
    M[42] = M[41]
    rep = identical_rows(sparse.csc_matrix(M), min_library=10)
    assert rep["repeats"] == 4 and rep["groups"] == 2 and rep["largest_group"] == 4
    assert rep["repeats_min_library"] == 3 and rep["empty_rows"] == 1


def test_qc_stage_reports_levels_sections_order_duplicates_and_centring():
    per_level = (150, 160, 320, 170, 150, 140)  # level 2 holds two sections
    counts, coords = _stage(per_level=per_level, seed=7)
    rng = np.random.default_rng(8)
    coords[:, 2] += rng.choice([0.0, 1e-13], size=coords.shape[0])
    lv = snap_z(coords[:, 2]).index
    sample = np.array([f"s{i}" for i in lv], dtype=object)
    at2 = np.flatnonzero(lv == 2)
    sample[at2[: at2.size // 2]] = "s2b"  # two sample_ids on the doubled level
    sample[lv == 5] = "s4"                # one sample_id on two levels
    keep = lv != 1                        # drop level 1: a one-level gap
    a = _adata(counts[keep], coords[keep], sample_id=sample[keep])
    dup = _adata(np.vstack([a.layers["counts"].toarray(), a.layers["counts"][:3].toarray()]),
                 np.vstack([a.obsm["spatial_3D"], a.obsm["spatial_3D"][:3] - [0, 0, STEP]]),
                 sample_id=np.r_[a.obs["sample_id"].to_numpy(), ["s4"] * 3])
    dup.obsm["spatial_3D"][:, :2] -= dup.obsm["spatial_3D"][:, :2].mean(0)  # centre xy only
    qc = qc_stage(dup)
    json.dumps(qc, allow_nan=False)
    assert qc["cells"] == dup.n_obs and qc["counts"]["non_integer_rows"] == 0
    z = qc["z"]
    assert z["levels"] == 5 and z["raw_values"] > 5 and z["step"] == pytest.approx(STEP)
    assert z["max_level_spread"] <= 1e-12
    assert z["gaps"] == [{"from": -37.5, "to": -7.5, "steps": 2.0, "missing": 1}]
    assert z["off_grid_levels"] == [] and z["dense_levels"] == {"-7.5": {"cells": 320, "neighbour_ratio": 2.0}}
    assert z["cells_per_level"] == [150, 320, 170, 153, 140]
    s = qc["sections"]["sample_id"]
    assert s["levels_with_several_labels"] == {"-7.5": ["s2", "s2b"]}
    assert s["labels_on_several_levels"] == {"s4": [22.5, 37.5]}
    assert not s["one_level_per_label"] and not s["one_label_per_level"]
    order = qc["row_order"]
    assert order["spearman_row_vs_z"] < -0.9 and order["direction"] == "descending"
    assert order["z_runs"] == 6 and not order["grouped_by_level"]  # the appended copies start a new run
    assert qc["duplicates"]["n_redundant"] == 3 and qc["duplicates"]["by_gap"] == {"1": 3}
    assert qc["duplicates_any_xy"]["n_redundant"] == 3 and qc["duplicates_any_xy"]["xy_tol"] == "any"
    assert qc["coords"]["centred_on"] != "mean"  # xy are centred, z is not
    assert max(abs(v) for v in qc["coords"]["mean"][:2]) < 1e-9 and abs(qc["coords"]["mean"][2]) > 1
    assert qc["coords"]["duplicate_positions"] == 0


def test_qc_stage_without_coordinates_reports_row_identity_only():
    M = np.random.default_rng(9).poisson(1.0, size=(40, 25)).astype(np.float32)
    M[[3, 4]] = M[2]
    a = ad.AnnData(X=sparse.csc_matrix(np.log1p(M)))
    qc = qc_stage(a)
    assert qc["matrix"] == "X" and "z" not in qc and "counts" not in qc
    assert qc["identical_rows"]["repeats"] == 2 and "repeats_min_library" not in qc["identical_rows"]


def test_qc_stage_flags_half_integer_counts_and_centring():
    counts, coords = _stage(per_level=(100,) * 4, seed=11)
    counts[[0, 5], 0] += 0.5
    coords -= coords.mean(0)
    qc = qc_stage(_adata(counts, coords))
    assert qc["counts"] == {"non_integer_entries": 2, "non_integer_rows": 2, "half_integer_entries": 2,
                            "negative_entries": 0}
    assert qc["coords"]["centred_on"] == "mean"
    assert qc["row_order"]["grouped_by_level"] and qc["row_order"]["z_runs"] == 4
