"""vec.t2_build.organ_scope_mask (the organ cut, file D; row R21) and its size helpers, on synthetic stages only.

A detached chamber-category island is dropped; a bridge of chamber-category cells joins; a non-chamber cell is never
kept, and never bridges; core cells are always kept; the mask is invariant to row order and deterministic; L is the
99th percentile of the core's nearest-core distances, and "closer than L" is strict; an unmapped label refuses; the
pair search scales to a few hundred thousand cells. Nothing is read from data/."""

from __future__ import annotations

import numpy as np
import pytest

from vec import t2_build

CORE_LAB, ENDO_LAB, VESSEL_LAB, EPI_LAB, OTHER_LAB = "V-CM", "Chamber-Endo", "BEC", "Epi", "Dorsal-PAM"


def _grid(n: int, origin, step: float = 1.0) -> np.ndarray:
    """n x n x 2 lattice points at spacing `step` from `origin`."""
    g = np.stack(np.meshgrid(np.arange(n), np.arange(n), np.arange(2), indexing="ij"), -1).reshape(-1, 3)
    return np.asarray(origin, dtype=float) + g * step


def _stage(rng=None):
    """A core block (lattice spacing 1, so L = 1), chamber-category cells touching it, a detached vessel island 10
    away, and non-chamber cells next to the core."""
    core = _grid(6, (0, 0, 0))
    touching = _grid(3, (5.6, 0, 0), 0.9)  # vessels next to the core, links of 0.6 and 0.9 < L
    island = _grid(3, (20, 0, 0), 0.9)  # vessels far away
    other = _grid(3, (-3, 0, 0), 0.9)  # pharyngeal/paraxial cells beside the core
    C = np.vstack([core, touching, island, other])
    labels = np.array([CORE_LAB] * (len(core) // 2) + [ENDO_LAB] * (len(core) - len(core) // 2)
                      + [VESSEL_LAB] * len(touching) + [EPI_LAB] * len(island) + [OTHER_LAB] * len(other))
    is_core = np.zeros(len(C), dtype=bool)
    is_core[:len(core)] = True
    parts = {"core": slice(0, len(core)), "touching": slice(len(core), len(core) + len(touching)),
             "island": slice(len(core) + len(touching), len(core) + len(touching) + len(island)),
             "other": slice(len(C) - len(other), len(C))}
    return C, labels, is_core, parts


def test_the_link_length_is_the_99th_percentile_of_core_nearest_core_distances():
    rng = np.random.default_rng(0)
    C = rng.normal(size=(500, 3))
    core = rng.random(500) < 0.4
    from scipy.spatial import cKDTree

    P = C[core]
    d = np.array([np.sort(np.linalg.norm(P - p, axis=1))[1] for p in P])
    assert t2_build.organ_link_length(C, core) == pytest.approx(np.percentile(d, 99.0), rel=1e-12)
    assert t2_build.ORGAN_LINK_PERCENTILE == 99.0
    assert cKDTree(P).query(P, k=2)[0][:, 1] == pytest.approx(d)
    with pytest.raises(ValueError, match="at least two core cells"):
        t2_build.organ_link_length(C, np.eye(1, 500, dtype=bool)[0])


def test_a_detached_island_is_dropped_and_touching_chamber_cells_kept():
    C, labels, core, parts = _stage()
    assert t2_build.organ_link_length(C, core) == pytest.approx(1.0)
    keep = t2_build.organ_scope_mask(C, labels, core)
    assert keep[parts["core"]].all() and keep[parts["touching"]].all()
    assert not keep[parts["island"]].any()
    assert not keep[parts["other"]].any()


def test_a_bridge_of_chamber_cells_joins_the_island():
    C, labels, core, parts = _stage()
    x = np.arange(8.2, 20.0, 0.8)  # from the touching block (x up to 7.4) to the island (x from 20) in steps < L
    bridge = np.column_stack([x, np.zeros(len(x)), np.zeros(len(x))])
    C2 = np.vstack([C, bridge])
    labels2 = np.concatenate([labels, [EPI_LAB] * len(bridge)])
    core2 = np.concatenate([core, np.zeros(len(bridge), dtype=bool)])
    keep = t2_build.organ_scope_mask(C2, labels2, core2)
    assert keep[parts["island"]].all() and keep[len(C):].all()
    # the same bridge of non-chamber cells joins nothing and is not kept
    labels3 = np.concatenate([labels, [OTHER_LAB] * len(bridge)])
    keep3 = t2_build.organ_scope_mask(C2, labels3, core2)
    assert not keep3[parts["island"]].any() and not keep3[len(C):].any()


def test_a_non_chamber_cell_is_never_kept_even_inside_the_core():
    C, labels, core, parts = _stage()
    inside = np.array([[2.5, 2.5, 0.5]])
    keep = t2_build.organ_scope_mask(np.vstack([C, inside]), np.append(labels, "Late-a-PHM"),
                                     np.append(core, False))
    assert not keep[-1] and keep[: len(C)][parts["core"]].all()


def test_closer_than_l_is_strict():
    """Two core cells at distance 1 (so L = 1) and a vessel at exactly 1 from the nearest core cell: not linked."""
    C = np.array([[0.0, 0, 0], [1.0, 0, 0], [2.0, 0, 0], [3.0, 0, 0]])
    labels = np.array([CORE_LAB, CORE_LAB, VESSEL_LAB, VESSEL_LAB])
    core = np.array([True, True, False, False])
    assert t2_build.organ_link_length(C, core) == 1.0
    assert list(t2_build.organ_scope_mask(C, labels, core)) == [True, True, False, False]
    C[2, 0] = 1.999
    # 0.999 from the core: linked; the last vessel is 1.001 from it, not closer than L, so it stays out
    assert list(t2_build.organ_scope_mask(C, labels, core)) == [True, True, True, False]


def test_core_cells_are_always_kept_and_must_be_chamber_labelled():
    C, labels, core, _ = _stage()
    far_core = np.array([[100.0, 100.0, 100.0]])
    keep = t2_build.organ_scope_mask(np.vstack([C, far_core]), np.append(labels, CORE_LAB), np.append(core, True))
    assert keep[-1]  # a core cell in its own component: kept
    with pytest.raises(ValueError, match="core cells with a label outside the chamber categories"):
        t2_build.organ_scope_mask(np.vstack([C, far_core]), np.append(labels, OTHER_LAB), np.append(core, True))


def test_the_mask_is_invariant_to_row_order_and_deterministic():
    rng = np.random.default_rng(3)
    n = 4000
    C = rng.normal(size=(n, 3)) * np.array([10.0, 10.0, 3.0])
    labels = rng.choice([CORE_LAB, ENDO_LAB, VESSEL_LAB, EPI_LAB, OTHER_LAB, "Unknown", "Blood"], size=n)
    core = np.isin(labels, [CORE_LAB, ENDO_LAB])
    keep, rec = t2_build.organ_scope(C, labels, core)
    assert 0 < keep.sum() < n and rec["components_total"] > rec["components_kept"] >= 1
    for seed in range(3):
        p = np.random.default_rng(seed).permutation(n)
        np.testing.assert_array_equal(t2_build.organ_scope_mask(C[p], labels[p], core[p]), keep[p])
    np.testing.assert_array_equal(t2_build.organ_scope_mask(C, labels, core), keep)
    assert t2_build.organ_scope(C, labels, core)[1] == rec


def test_the_record_names_the_rule_config_and_l():
    C, labels, core, _ = _stage()
    keep, rec = t2_build.organ_scope(C, labels, core)
    spec = t2_build.load_t2_chamber_spec()
    assert rec["label_config"] == "configs/t2_chamber_classes.yaml" and rec["label_config_sha256"] == spec.sha256
    assert rec["link_percentile"] == 99.0 and rec["link_length"] == pytest.approx(1.0)
    assert rec["kept_categories"] == ["cm", "endocardium", "epicardium", "vascular_endothelium", "blood"]
    assert rec["kept_share"] == pytest.approx(keep.mean())
    assert rec["core_share_of_kept"] == pytest.approx(core.sum() / keep.sum())
    assert rec["rule"] == t2_build.ORGAN_SCOPE_RULE


def test_an_unmapped_label_refuses_the_mask():
    C, labels, core, _ = _stage()
    labels = labels.copy()
    labels[-1] = "Mystery"
    with pytest.raises(t2_build.UnmappedT2LabelError, match="Mystery"):
        t2_build.organ_scope_mask(C, labels, core)


def test_a_few_hundred_thousand_cells_with_a_k_d_tree():
    """300,000 cells: the pair list of a k-d tree, no dense matrix (a dense float64 one would take 720 GB)."""
    rng = np.random.default_rng(7)
    n = 300_000
    C = rng.uniform(0, 100, size=(n, 3)) * np.array([1.0, 1.0, 0.2])
    labels = rng.choice(np.array([CORE_LAB, VESSEL_LAB, OTHER_LAB]), size=n, p=[0.3, 0.3, 0.4])
    core = labels == CORE_LAB
    keep = t2_build.organ_scope_mask(C, labels, core)
    assert keep[core].all() and not keep[labels == OTHER_LAB].any()


# ------------------------------------------------------------------------------------------------ the size
def test_the_organ_core_rms_ratio():
    C, labels, core, _ = _stage()
    keep = t2_build.organ_scope_mask(C, labels, core)
    r = t2_build.organ_core_rms_ratio(C, keep, core)
    assert r == pytest.approx(t2_build.rms_radius(C[keep]) / t2_build.rms_radius(C[core]))
    assert r > 1.0
    with pytest.raises(ValueError, match="every core cell"):
        t2_build.organ_core_rms_ratio(C, keep & ~core, core)


def test_the_size_field_is_new_and_the_old_ones_unchanged():
    assert t2_build.SIZE_ORGAN_SCOPE == "organ"
    assert t2_build.SIZE_ORGAN_SCOPE not in (*t2_build.SIZE_HYPOTHESES, t2_build.SIZE_CORE_SCOPE)
    assert t2_build.SIZE_HYPOTHESES == ("H1", "H0", "H2") and t2_build.SIZE_CORE_SCOPE == "core"
    assert t2_build.SIZE_CORE_SCOPE_MODES == ("raw",)
    assert "z_modes.raw.pred_rms_core" in t2_build.SIZE_ORGAN_SOURCE and "E10.5" in t2_build.SIZE_ORGAN_SOURCE


@pytest.mark.parametrize("ratio", [0.0, -1.0, float("nan"), float("inf"), True, "2"])
def test_the_size_refuses_a_bad_ratio(ratio):
    with pytest.raises(ValueError, match="ratio"):
        t2_build.size_organ_scope(np.ones((4, 3)), {}, ratio, source_stage="E10.5", reference_cloud=np.ones((4, 3)))
