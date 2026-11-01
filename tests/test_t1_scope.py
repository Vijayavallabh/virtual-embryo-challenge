import hashlib
import importlib.util
import inspect
import json
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pytest
import yaml
from scipy import sparse

from vec.boards import BOARDS_DIR, REPO_ROOT, Board
from vec.prep import log1p_cp10k
from vec.provenance import read_stamp, stamp_path
from vec.t1_scope import (
    UnmappedLabelError,
    apply_factors,
    build_submission,
    check_spec,
    classify,
    composition,
    counts_of,
    external_shift,
    load_shares,
    load_spec,
    mean_cp10k,
    momentum_shift,
    quotas,
    reweight,
    share_of,
    with_share,
    write_checked,
)

GENES = [f"g{i}" for i in range(78)] + ["Xist", "Ddx3y"]
HELD = ["Xist", "Ddx3y"]
# released labels -> chamber categories (configs/t1_chamber_classes.yaml): cm, endocardium, epicardium,
# vascular_endothelium, blood, progenitor, progenitor, surrounding, surrounding
SIZES = {"V-CM": 400, "Endocardium": 150, "Proepicardium": 60, "BEC": 40, "Blood": 30, "aSHF": 120, "JCF": 80,
         "Foregut": 200, "Neural Tube": 120}
SHARES = {"cm": 0.42, "endocardium": 0.25, "epicardium": 0.1, "vascular_endothelium": 0.06, "blood": 0.05,
          "surrounding": 0.12}  # synthetic test values; no E12.5 share is set before the release
# E8.5- and E9.5-like stages for the per-file rules: E8.5 has one mixed "Endothelium" label (endocardium by the E8.5
# override, and the vascular stand-in) and JCF (the epicardium stand-in), but no Endocardium, BEC or Proepicardium
E85_SIZES = {"V-CM": 400, "Endothelium": 150, "Blood": 30, "aSHF": 120, "JCF": 80, "Foregut": 200}
E95_SIZES = {"V-CM": 400, "Endocardium": 150, "Proepicardium": 60, "BEC": 40, "Blood": 30, "aSHF": 120, "Foregut": 200}


def _profiles(seed=0):
    rng = np.random.default_rng(seed)
    out = {}
    for lab in SIZES:
        mu = rng.gamma(0.4, 1.0, size=len(GENES))  # many low genes, so cells have exact zeros
        mu[GENES.index("Xist")] = 3.0
        mu = 2000.0 * mu / mu.sum()
        mu[GENES.index("g1")] = 20.0  # a moderately expressed gene for the direction test
        out[lab] = mu
    out["Endothelium"] = out["Endocardium"]  # mostly endocardial (the E8.5 override)
    return out


def make_stage(seed, sizes=SIZES, profiles=None, boost=None):
    """Poisson counts per label with cell-depth variation; boost = (label, gene, factor) scales one gene's mean."""
    rng = np.random.default_rng(seed)
    profiles = profiles or _profiles()
    blocks, labels = [], []
    for lab, n in sizes.items():
        mu = profiles[lab].copy()
        if boost and boost[0] == lab:
            mu[GENES.index(boost[1])] *= boost[2]
        C = rng.poisson(mu * rng.gamma(10.0, 0.1, size=(n, 1)))
        C[np.arange(n), rng.integers(0, len(GENES), n)] += 1  # every cell has counts
        blocks.append(C)
        labels += [lab] * n
    C = sparse.csr_matrix(np.vstack(blocks).astype(np.int32))
    a = ad.AnnData(X=log1p_cp10k(C))
    a.var_names = GENES
    a.obs_names = [f"s{seed}_{i}" for i in range(a.n_obs)]
    a.obs["celltype"] = labels
    a.layers["counts"] = C
    return a


def t1_board(genes, min_cells=10, max_cells=1000):
    return Board(key="T1:synthetic", task="T1", genes=list(genes), needs_coords=False, min_cells=min_cells,
                 max_cells=max_cells)


# ---------------------------------------------------------------------------------------------------------------------
# config and classification


def test_config_maps_every_released_t1_label_consistently():
    spec = load_spec()
    released = json.loads((BOARDS_DIR / "t1_composition.json").read_text())["stages"]
    labels = {lab for stage in ("E8.5", "E9.5") for lab in released[stage]["composition"]}
    assert labels <= set(spec.labels), sorted(labels - set(spec.labels))
    assert check_spec(spec) == []
    assert set(spec.heart) <= set(spec.in_scope("chamber"))
    assert spec.labels["ST"] == "surrounding" and "ST" in spec.exceptions


def test_config_maps_every_label_of_the_released_e10_5_file():
    """E10.5_RNA (validation answer, unlocked 2026-10-28) is not in boards/t1_composition.json yet (the panel files
    answer 403), so its 23 labels are listed here."""
    from test_prep import VALIDATION_LABELS

    spec = load_spec()
    e105 = VALIDATION_LABELS["T1/E10.5_RNA"]
    assert len(e105) == 23 and e105 <= set(spec.labels), sorted(e105 - set(spec.labels))
    assert check_spec(spec) == []
    assert spec.labels["Early VIC"] == "endocardium" and "Early VIC" in spec.exceptions
    assert spec.labels["Macrophage"] == "blood" and spec.labels["Epicardium"] == "epicardium"


def test_no_default_target_composition_or_size():
    for fn, params in ((load_shares, ["shares"]), (reweight, ["target_shares", "n_cells"]),
                       (build_submission, ["target_shares", "n_cells"])):
        sig = inspect.signature(fn)
        for p in params:
            assert sig.parameters[p].default is inspect.Parameter.empty, (fn.__name__, p)
    assert "shares:" not in Path(load_spec().path).read_text()


def test_unmapped_labels_fail_loudly():
    a = make_stage(0)
    a.obs["celltype"] = a.obs["celltype"].astype(str)
    a.obs.iloc[:3, a.obs.columns.get_loc("celltype")] = "Mystery Cell"
    with pytest.raises(UnmappedLabelError) as err:
        classify(a)
    assert err.value.labels == {"Mystery Cell": 3}
    assert "category" not in a.obs


def test_classify_and_file_overrides():
    a = ad.AnnData(X=np.zeros((4, 2), np.float32))
    a.obs["celltype"] = ["V-CM", "Endothelium", "ST", "Proepicardium"]
    assert list(classify(a)) == ["cm", "vascular_endothelium", "surrounding", "epicardium"]
    assert list(a.obs["category"]) == ["cm", "vascular_endothelium", "surrounding", "epicardium"]
    assert classify(a, source="E8.5_RNA", key_added=None)[1] == "endocardium"
    assert load_spec().fallbacks("E8.5_RNA")["epicardium"] == ("JCF",)
    assert load_spec().fallbacks("E9.5_RNA") == {}


def test_bad_config_rejected(tmp_path):
    p = tmp_path / "bad.yaml"
    p.write_text("categories: {cm: {scope: chamber}}\nheart: [cm]\nlabels: {V-CM: cm, BEC: vessels}\n")
    with pytest.raises(ValueError, match="undefined category 'vessels'"):
        load_spec(p)


def test_a_source_is_a_released_file_stem_a_config_key_or_empty(tmp_path):
    a = ad.AnnData(X=np.zeros((2, 2), np.float32))
    a.obs["celltype"] = ["Endothelium", "JCF"]
    for source, endothelium in ((None, "vascular_endothelium"), ("", "vascular_endothelium"),
                                ("E9.5_RNA", "vascular_endothelium"), ("E10.5_RNA", "vascular_endothelium"),
                                ("E8.5_RNA", "endocardium")):
        assert classify(a, source=source, key_added=None)[0] == endothelium
    for bad in ("E8.5", "E8.5_rna", "E8.5_RNA.h5ad", "E12.5_RNA"):  # a stage name, typos, a stage never released
        with pytest.raises(ValueError, match="unknown source"):  # each would silently select no per-file rules
            classify(a, source=bad)
        with pytest.raises(ValueError, match="unknown source"):
            load_spec().fallbacks(bad)
    assert "category" not in a.obs
    p = tmp_path / "own_keys.yaml"
    p.write_text("categories: {cm: {scope: chamber}}\nheart: [cm]\nlabels: {V-CM: cm}\n"
                 "file_fallbacks: {lab_copy: {cm: [V-CM]}}\n")
    assert load_spec(p).fallbacks("lab_copy") == {"cm": ("V-CM",)}


def test_check_spec_reports_per_file_keys_that_are_not_released_file_stems(tmp_path):
    """A per-file section keyed by a stage name fails the config test (check_spec(spec) == []): the file's stem
    would select none of its rules, and nothing else would notice."""
    raw = yaml.safe_load(Path(load_spec().path).read_text())
    for key in ("E10.5", "E10.5_rna", "E12.5_RNA"):  # a stage name, a typo, a stage never released
        p = tmp_path / f"{key}.yaml"
        sections = {"file_overrides": {**raw["file_overrides"], key: {"Endothelium": "endocardium"}},
                    "file_fallbacks": {**raw["file_fallbacks"], key: {"epicardium": ["JCF"]}}}
        p.write_text(yaml.safe_dump({**raw, **sections}, sort_keys=False))
        spec = load_spec(p)
        problems = check_spec(spec)
        assert [x.split(f" key {key!r} ")[0] for x in problems] == ["file_overrides", "file_fallbacks"], problems
        assert spec.fallbacks(key) == {"epicardium": ("JCF",)}  # the library still accepts it when passed
    assert load_spec(tmp_path / "E10.5.yaml").fallbacks("E10.5_RNA") == {}


# ---------------------------------------------------------------------------------------------------------------------
# target shares


def test_load_shares_validates(tmp_path):
    full = load_shares(SHARES)
    assert set(full) == set(load_spec().categories) and full["pericardium"] == 0.0
    with pytest.raises(ValueError, match="unknown categories"):
        load_shares({**SHARES, "cushion": 0.0})
    with pytest.raises(ValueError, match="sum to"):
        load_shares({**SHARES, "cm": 0.5})
    with pytest.raises(ValueError, match="lie in"):
        load_shares({"cm": 1.2, "blood": -0.2})
    (tmp_path / "a.yaml").write_text("shares: {cm: 0.6, endocardium: 0.4}\n")
    (tmp_path / "b.yaml").write_text("cm: 0.6\nendocardium: 0.4\n")
    assert load_shares(tmp_path / "a.yaml") == load_shares(tmp_path / "b.yaml")


def test_with_share_moves_only_the_cm_share_of_heart_cells():
    heart = load_spec().heart
    new = with_share(SHARES, "cm", 0.6)
    assert share_of(new, "cm", heart) == pytest.approx(0.6)
    assert sum(new[c] for c in heart) == pytest.approx(sum(SHARES[c] for c in heart))
    for c in ("blood", "surrounding"):
        assert new[c] == pytest.approx(SHARES[c])
    ratio = new["endocardium"] / new["epicardium"]
    assert ratio == pytest.approx(SHARES["endocardium"] / SHARES["epicardium"])
    assert sum(new.values()) == pytest.approx(1.0)


@pytest.mark.parametrize("n", [0, 1, 7, 333, 5000])
def test_quotas_sum_exactly(n):
    shares = load_shares(SHARES)
    q = quotas(shares, n)
    assert sum(q.values()) == n
    for c, s in shares.items():
        assert abs(q[c] - s * n) < 1.0
        if s == 0:
            assert q[c] == 0


# ---------------------------------------------------------------------------------------------------------------------
# re-weighting


def test_reweight_hits_quotas_exactly_with_real_cells():
    a = make_stage(1)
    classify(a)
    out, rep = reweight(a, SHARES, 500, seed=3)
    q = quotas(load_shares(SHARES), 500)
    got = composition(out)
    for c, n in q.items():
        assert round(got.get(c, 0.0) * 500) == n
    assert out.n_obs == 500 and out.obs_names.is_unique
    assert rep.reused == []  # every category had enough cells: no repeats
    assert out.obs["source_obs"].is_unique
    src = a[out.obs["source_obs"].to_numpy()]
    assert (out.X != src.X).count_nonzero() == 0  # rows are the released values of real cells
    assert (out.layers["counts"] != src.layers["counts"]).count_nonzero() == 0
    again, _ = reweight(a, SHARES, 500, seed=3)
    assert list(again.obs["source_obs"]) == list(out.obs["source_obs"])


def test_reweight_reuses_as_few_cells_as_possible_and_reports_it():
    a = make_stage(2)
    classify(a)
    shares = {"cm": 0.5, "blood": 0.5}  # 30 blood cells for a quota of 100
    out, rep = reweight(a, shares, 200, seed=0)
    row = next(r for r in rep.rows if r["category"] == "blood")
    assert rep.reused == ["blood"] and row["reused"] == 70 and row["unique"] == 30
    mult = out.obs.loc[out.obs["category"] == "blood", "source_obs"].value_counts()
    assert mult.max() - mult.min() <= 1 and mult.max() == 4  # ceil(100 / 30)


def test_within_restricts_labels_and_fallback_fills_missing_categories():
    sizes = {k: v for k, v in SIZES.items() if k not in ("Proepicardium", "BEC")}
    a = make_stage(3, sizes=sizes)
    classify(a)
    with pytest.raises(ValueError, match="no cells to draw from"):
        reweight(a, SHARES, 400, seed=0)
    out, rep = reweight(a, SHARES, 400, seed=0, within={"surrounding": ["Foregut"]},
                        fallback={"epicardium": ["JCF"], "vascular_endothelium": ["Endocardium"]})
    obs = out.obs
    assert set(obs.loc[obs["category"] == "surrounding", "celltype"]) == {"Foregut"}
    assert set(obs.loc[obs["category"] == "epicardium", "celltype"]) == {"JCF"}
    rows = {r["category"]: r for r in rep.rows}
    assert rows["epicardium"]["drawn_from"] == "stand-in" and rows["surrounding"]["drawn_from"] == "own, restricted"
    # vascular stand-ins share the endocardium pool, which has enough unused cells: no cell is drawn twice
    assert rep.reused == [] and obs["source_obs"].is_unique
    assert (obs.loc[obs["category"] == "vascular_endothelium", "source_category"] == "endocardium").all()


# ---------------------------------------------------------------------------------------------------------------------
# count-space shift


def _shift_inputs():
    start = make_stage(10)
    ref = make_stage(11, boost=("V-CM", "g1", 4.0))  # g1 rises 4x in cardiomyocytes between start and ref
    for x in (start, ref):
        classify(x)
    pop, _ = reweight(ref, SHARES, 600, seed=1)
    return start, ref, pop


def test_shift_keeps_integers_zeros_and_held_flat_genes():
    start, ref, pop = _shift_inputs()
    out, rep = momentum_shift(start, ref, pop, damping=0.5, held_flat=HELD, seed=0)
    old, new = counts_of(pop), counts_of(out)
    assert new.dtype.kind == "i" and (new.data > 0).all()
    nz_new, nz_old = (new != 0).astype(np.int8), (new != 0).multiply(old != 0).astype(np.int8)
    assert (nz_new - nz_old).count_nonzero() == 0  # no zero became non-zero
    held = [GENES.index(g) for g in HELD]
    assert (new[:, held] != old[:, held]).count_nonzero() == 0
    assert (new != old).count_nonzero() > 0
    X = out.X.toarray() if sparse.issparse(out.X) else out.X
    np.testing.assert_allclose(X, log1p_cp10k(new).toarray(), atol=1e-6)
    assert X.dtype == np.float32
    assert rep.restored_cells == 0 and rep.n_held_flat == 2


def test_shift_follows_the_damped_direction_and_keeps_depth():
    start, ref, pop = _shift_inputs()
    out, rep = momentum_shift(start, ref, pop, damping=0.5, held_flat=HELD, seed=0)
    cm = (pop.obs["category"] == "cm").to_numpy()
    old, new = counts_of(pop)[cm], counts_of(out)[cm]
    g1 = GENES.index("g1")
    m_s = mean_cp10k(counts_of(start)[(start.obs["category"] == "cm").to_numpy()])[g1]
    m_r = mean_cp10k(counts_of(ref)[(ref.obs["category"] == "cm").to_numpy()])[g1]
    f = np.sqrt((m_r + 0.5) / (m_s + 0.5))
    rest = [i for i in range(len(GENES)) if i != g1 and GENES[i] not in HELD]
    ratio_g1 = new[:, g1].sum() / old[:, g1].sum()
    ratio_rest = new[:, rest].sum() / old[:, rest].sum()
    assert f > 1.7
    assert ratio_g1 / ratio_rest == pytest.approx(f, rel=0.08)
    depth = np.asarray(counts_of(out).sum(1)).ravel() / np.asarray(counts_of(pop).sum(1)).ravel()
    assert depth.mean() == pytest.approx(1.0, abs=0.01)
    assert {r["type"]: r["status"] for r in rep.rows}["cm"] == "shifted"


def test_zero_damping_is_the_identity_and_small_types_are_skipped():
    start, ref, pop = _shift_inputs()
    out, _ = momentum_shift(start, ref, pop, damping=0.0, held_flat=HELD, seed=0)
    assert (counts_of(out) != counts_of(pop)).count_nonzero() == 0
    _, rep = momentum_shift(start, ref, pop, damping=0.5, min_cells=35, seed=0)
    assert rep.skipped == ["blood"]  # 30 blood cells per stage


def test_shift_uses_stand_ins_for_a_type_missing_at_the_start():
    start = make_stage(12, sizes={k: v for k, v in SIZES.items() if k != "Proepicardium"})
    ref = make_stage(13)
    for x in (start, ref):
        classify(x)
    pop, _ = reweight(ref, SHARES, 400, seed=2)
    _, rep = momentum_shift(start, ref, pop, damping=0.5, seed=0)
    assert "epicardium" in rep.skipped
    _, rep = momentum_shift(start, ref, pop, damping=0.5, start_fallback={"epicardium": ["JCF"]},
                            standin_momentum=True, seed=0)
    row = {r["type"]: r for r in rep.rows}["epicardium"]
    assert row["status"] == "shifted" and row["start_cells"] == "stand-in JCF" and row["n_start"] == 80


def test_apply_factors_expectation():
    rng = np.random.default_rng(0)
    C = sparse.csr_matrix(rng.poisson(5.0, size=(4000, 3)))
    F = np.array([[0.5, 1.0, 2.0]])
    new, scale, dead = apply_factors(C, F, np.zeros(4000, dtype=int), keep_depth=False, seed=1)
    ratio = np.asarray(new.sum(0)).ravel() / np.asarray(C.sum(0)).ravel()
    np.testing.assert_allclose(ratio, [0.5, 1.0, 2.0], rtol=0.03)
    assert ratio[1] == 1.0 and not dead.any() and (scale == 1).all()
    with pytest.raises(ValueError, match="positive"):
        apply_factors(C, np.array([[0.0, 1.0, 1.0]]), np.zeros(4000, dtype=int))


# ---------------------------------------------------------------------------------------------------------------------
# submission


def _direction(**genes):
    d = np.zeros(len(GENES))
    for g, v in genes.items():
        d[GENES.index(g)] = v
    return d


def test_external_shift_moves_only_mapped_cells_by_the_fraction_of_the_direction():
    ref = make_stage(30)
    classify(ref)
    pop, _ = reweight(ref, SHARES, 600, seed=1)
    d = _direction(g1=np.log(4.0), Xist=np.log(8.0))
    d[GENES.index("g2")] = np.nan  # no external value: left alone
    out, rep = external_shift(pop, {"CM": d, "Absent": d}, 0.5, {"V-CM": "CM"}, HELD, seed=0)
    old, new = counts_of(pop), counts_of(out)
    cm = (pop.obs["celltype"] == "V-CM").to_numpy()
    assert (new[~cm] != old[~cm]).count_nonzero() == 0  # cells without a direction keep their counts
    g1, g2, xist = GENES.index("g1"), GENES.index("g2"), GENES.index("Xist")
    assert (new[:, xist] != old[:, xist]).count_nonzero() == 0  # held flat
    rest = [i for i in range(len(GENES)) if i not in (g1, g2) and GENES[i] not in HELD]
    base = new[cm][:, rest].sum() / old[cm][:, rest].sum()
    assert (new[cm][:, g1].sum() / old[cm][:, g1].sum()) / base == pytest.approx(2.0, rel=0.08)  # exp(0.5 log 4)
    assert (new[cm][:, g2].sum() / old[cm][:, g2].sum()) / base == pytest.approx(1.0, rel=0.08)
    nz_new, nz_old = (new != 0).astype(np.int8), (new != 0).multiply(old != 0).astype(np.int8)
    assert (nz_new - nz_old).count_nonzero() == 0  # zeros stay zero
    rows = {r["type"]: r for r in rep.rows}
    assert set(rows) == {"CM", "(no direction)"} and rows["CM"]["n_target"] == int(cm.sum())
    assert rows["CM"]["n_no_value"] == 1 and rep.damping == 0.5
    same, _ = external_shift(pop, {"CM": d}, 0.0, {"V-CM": "CM"}, HELD, seed=0)
    assert (counts_of(same) != old).count_nonzero() == 0  # fraction 0 is the identity
    with pytest.raises(ValueError, match="shape"):
        external_shift(pop, {"CM": d[:-1]}, 0.5, {"V-CM": "CM"})
    with pytest.raises(ValueError, match="finite"):
        external_shift(pop, {"CM": d}, float("nan"), {"V-CM": "CM"})


def test_build_submission_moves_mapped_cells_by_direction_and_the_rest_by_momentum():
    start, ref = make_stage(40), make_stage(41, boost=("V-CM", "g1", 4.0))
    board = t1_board(GENES[::-1], min_cells=100, max_cells=600)
    dirs, groups = {"CM": _direction(g5=np.log(9.0))}, {"V-CM": "CM"}
    plain, rep0 = build_submission(ref, board, SHARES, 400, seed=3, ref_source="")
    sub, rep = build_submission(ref, board, SHARES, 400, seed=3, ref_source="", directions=dirs,
                                direction_fraction=0.5, direction_groups=groups, held_flat=HELD)
    assert rep["realised"] == rep0["realised"] and rep["shift"] is None and rep["external"]["fraction"] == 0.5
    src = {r.tobytes() for r in ref[:, board.genes].X.toarray().astype(np.float32)}
    labels = ref.obs["celltype"].to_numpy()
    cats = {r["type"]: r for r in rep["external"]["rows"]}
    assert cats["CM"]["n_target"] == round(400 * SHARES["cm"]) and cats["CM"]["labels"] == "V-CM"
    unchanged = np.array([r.tobytes() in src for r in sub.X])
    assert unchanged.sum() == 400 - cats["CM"]["n_target"]  # only the mapped cells moved (no momentum given)
    g5 = board.genes.index("g5")
    cm_ref = labels == "V-CM"
    m_ref = np.expm1(ref[:, board.genes].X.toarray()[cm_ref, g5]).mean()
    assert np.expm1(sub.X[~unchanged, g5].astype(np.float64)).mean() / m_ref == pytest.approx(3.0, rel=0.15)
    both, rep2 = build_submission(ref, board, SHARES, 400, seed=3, ref_source="", start=start, start_source="",
                                  damping=0.5, directions=dirs, direction_fraction=0.5, direction_groups=groups,
                                  held_flat=HELD)
    shifted = {r["type"] for r in rep2["shift"]["rows"]}
    assert "cm" not in shifted and "endocardium" in shifted  # momentum moved only the cells without a direction
    assert {r["type"] for r in rep2["external"]["rows"]} == {"CM", "(no direction)"}
    with pytest.raises(ValueError, match="direction_groups"):
        build_submission(ref, board, SHARES, 400, seed=3, ref_source="", directions=dirs, direction_fraction=0.5)


def test_build_submission_passes_the_format_check(tmp_path, own_repo):
    start, ref = make_stage(20), make_stage(21, boost=("V-CM", "g1", 4.0))
    board = t1_board(GENES[::-1], min_cells=100, max_cells=600)  # reversed order: the builder must reorder
    sub, rep = build_submission(ref, board, SHARES, 300, seed=0, start=start, damping=0.5, held_flat=HELD,
                                ref_source="", start_source="")
    assert list(sub.var_names) == board.genes and sub.n_obs == 300
    assert sub.X.dtype == np.float32 and sub.obs_names.is_unique and len(sub.obs.columns) == 0
    np.testing.assert_allclose(np.expm1(sub.X.astype(np.float64)).sum(1), 1e4, rtol=1e-4)
    assert rep["shift"]["damping"] == 0.5 and rep["reweight"]["n_cells"] == 300
    check = write_checked(sub, tmp_path / "pred.h5ad", board)
    assert check.ok and check.n_cells == 300 and check.warnings == []
    stamp = read_stamp(stamp_path(tmp_path / "pred.h5ad"))  # stamped by default, against the test's repository
    assert stamp["clean"] and stamp["sha256"] == hashlib.sha256((tmp_path / "pred.h5ad").read_bytes()).hexdigest()


def test_build_submission_without_shift_returns_released_rows():
    ref = make_stage(22)
    board = t1_board(GENES[::-1], min_cells=100, max_cells=600)
    sub, rep = build_submission(ref, board, SHARES, 250, seed=4, ref_source="")
    assert rep["shift"] is None
    src = ref[:, board.genes].X.toarray()
    rows = {r.tobytes() for r in src.astype(np.float32)}
    assert all(r.tobytes() in rows for r in sub.X)


def test_build_submission_refuses_bad_inputs():
    ref = make_stage(23)
    board = t1_board(GENES, min_cells=100, max_cells=600)
    with pytest.raises(ValueError, match="outside the board's limits"):
        build_submission(ref, board, SHARES, 601)
    with pytest.raises(ValueError, match="needs the earlier stage"):
        build_submission(ref, board, SHARES, 300, damping=0.5)
    with pytest.raises(ValueError, match="T1 predictions"):
        build_submission(ref, Board(key="T3:x", task="T3", genes=GENES, needs_coords=True), SHARES, 300)
    ref.obs["celltype"] = ref.obs["celltype"].astype(str)
    ref.obs.iloc[0, ref.obs.columns.get_loc("celltype")] = "Newcomer"
    with pytest.raises(ValueError, match="per-file overrides"):  # the file's per-file rules must be chosen
        build_submission(ref, board, SHARES, 300)
    with pytest.raises(UnmappedLabelError):
        build_submission(ref, board, SHARES, 300, ref_source="")


def test_build_submission_needs_a_known_source_for_each_classified_file():
    ref, start = make_stage(24), make_stage(25)
    board = t1_board(GENES, min_cells=100, max_cells=600)
    shift = {"start": start, "damping": 0.5, "held_flat": HELD}
    with pytest.raises(ValueError, match="per-file overrides"):  # the reference's rules must be chosen
        build_submission(ref, board, SHARES, 300, start_source="")
    with pytest.raises(ValueError, match="per-file overrides"):  # and, with a shift, the start's
        build_submission(ref, board, SHARES, 300, ref_source="", **shift)
    _, rep = build_submission(ref, board, SHARES, 300, ref_source="", start=start)  # no shift: start not classified
    assert rep["shift"] is None and "category" not in start.obs
    for bad in ({"ref_source": "E9.5"}, {"ref_source": "", "start_source": "E8.5"},
                {"ref_source": "E9.5_RNA", "start_source": "E8.5", **shift}):  # stage names in place of file stems
        with pytest.raises(ValueError, match="unknown source"):
            build_submission(ref, board, SHARES, 300, **bad)
    no_epi = make_stage(26, sizes={k: v for k, v in SIZES.items() if k != "Proepicardium"})
    _, rep = build_submission(ref, board, SHARES, 300, ref_source="E10.5_RNA", start_source="E8.5_RNA",
                              **{**shift, "start": no_epi})
    epi = next(r for r in rep["shift"]["rows"] if r["type"] == "epicardium")
    assert epi["start_cells"] == "stand-in JCF" and epi["status"].startswith("skipped: stand-in")  # E8.5's stand-in
    assert rep["shift"]["standin_momentum"] is False


def test_build_submission_classifies_each_file_by_its_own_rules():
    """Momentum from an E8.5-like start to an E9.5-like reference (the release-day geometry). By the start file's
    rules its Endothelium is endocardium's own cells and the vascular type's stand-in: only the latter is skipped."""
    start, ref = make_stage(40, sizes=E85_SIZES), make_stage(41, sizes=E95_SIZES)
    board = t1_board(GENES, min_cells=100, max_cells=600)
    _, rep = build_submission(ref, board, SHARES, 300, start=start, damping=0.5, held_flat=HELD, ref_source="E9.5_RNA",
                              start_source="E8.5_RNA")
    rows = {r["type"]: r for r in rep["shift"]["rows"]}
    endo, vasc, epi = rows["endocardium"], rows["vascular_endothelium"], rows["epicardium"]
    assert endo["start_cells"] == "own" and endo["n_start"] == 150 and endo["status"] == "shifted"
    assert vasc["start_cells"] == "stand-in Endothelium" and vasc["status"].startswith("skipped: stand-in")
    assert epi["start_cells"] == "stand-in JCF" and epi["status"].startswith("skipped: stand-in")


def test_build_submission_draws_and_shifts_with_the_reference_files_stand_ins():
    """An E8.5-like reference: the config's E8.5 stand-ins fill the categories it lacks and stand for those types at
    the reference end of the shift, where they get momentum only with standin_momentum=True."""
    start, ref = make_stage(42), make_stage(43, sizes=E85_SIZES)  # the start has every type's own cells
    board = t1_board(GENES, min_cells=100, max_cells=600)
    for standin in (False, True):
        _, rep = build_submission(ref, board, SHARES, 300, start=start, damping=0.5, held_flat=HELD,
                                  ref_source="E8.5_RNA", start_source="", standin_momentum=standin)
        assert rep["fallback"] == {c: list(labs) for c, labs in load_spec().fallbacks("E8.5_RNA").items()}
        drawn = {r["category"]: (r["drawn_from"], r["labels"]) for r in rep["reweight"]["rows"]}
        assert {c: drawn[c] for c in ("endocardium", "vascular_endothelium", "epicardium")} == {
            "endocardium": ("own", "Endothelium"), "vascular_endothelium": ("stand-in", "Endothelium"),
            "epicardium": ("stand-in", "JCF")}
        rows = {r["type"]: r for r in rep["shift"]["rows"]}
        assert rows["endocardium"]["ref_cells"] == "own" and rows["endocardium"]["status"] == "shifted"
        for t, labs in (("vascular_endothelium", "Endothelium"), ("epicardium", "JCF")):
            assert rows[t]["start_cells"] == "own" and rows[t]["ref_cells"] == f"stand-in {labs}"
            assert rows[t]["status"].startswith("shifted" if standin else "skipped: stand-in"), (t, rows[t]["status"])


def test_build_submission_passes_within_and_an_explicit_fallback_through():
    """`within` restricts a category to its labels in the draw and at both ends of the shift; an explicit `fallback`
    replaces the config's stand-ins in the draw and in the report."""
    start, ref = make_stage(44), make_stage(45)
    board = t1_board(GENES, min_cells=100, max_cells=600)
    _, rep = build_submission(ref, board, SHARES, 300, start=start, damping=0.5, held_flat=HELD, ref_source="",
                              start_source="", within={"surrounding": ["Foregut"]})
    drawn = next(r for r in rep["reweight"]["rows"] if r["category"] == "surrounding")
    assert (drawn["drawn_from"], drawn["labels"], drawn["pool"]) == ("own, restricted", "Foregut", 200)
    row = next(r for r in rep["shift"]["rows"] if r["type"] == "surrounding")
    assert row["start_cells"] == row["ref_cells"] == "own, restricted" and row["n_start"] == row["n_ref"] == 200
    ref = make_stage(46, sizes={k: v for k, v in SIZES.items() if k != "Proepicardium"})
    _, rep = build_submission(ref, board, SHARES, 300, ref_source="", fallback={"epicardium": ["JCF"]})
    assert rep["fallback"] == {"epicardium": ["JCF"]}
    drawn = next(r for r in rep["reweight"]["rows"] if r["category"] == "epicardium")
    assert (drawn["drawn_from"], drawn["labels"]) == ("stand-in", "JCF")


@pytest.mark.parametrize("stand_in_at", ["start", "ref", "both"])
def test_stand_in_types_get_no_momentum_unless_asked(stand_in_at):
    """A type drawn from stand-ins at either stage is skipped by default. The shipped geometries (the pseudo-task's
    shift arms, and momentum from E8.5 on release day) have the stand-in at the start only."""
    no_epi = {k: v for k, v in SIZES.items() if k != "Proepicardium"}  # without Proepicardium, JCF stands in for it
    start = make_stage(30, sizes=no_epi if stand_in_at in ("start", "both") else SIZES)
    ref = make_stage(31, sizes=no_epi if stand_in_at in ("ref", "both") else SIZES)
    for a in (start, ref):
        classify(a, source="")
    fb = {"epicardium": ["JCF"]}
    pop, _ = reweight(ref, SHARES, 400, seed=0, fallback=fb)
    _, rep = momentum_shift(start, ref, pop, 1.0, HELD, start_fallback=fb, ref_fallback=fb)
    row = next(r for r in rep.rows if r["type"] == "epicardium")
    stand_ins = [s for s in ("start", "ref") if row[f"{s}_cells"] == "stand-in JCF"]
    assert stand_ins == (["start", "ref"] if stand_in_at == "both" else [stand_in_at])
    assert row["status"].startswith("skipped: stand-in") and row["max_abs_log_f"] == 0.0
    assert rep.as_dict()["standin_momentum"] is False
    _, rep = momentum_shift(start, ref, pop, 1.0, HELD, start_fallback=fb, ref_fallback=fb, standin_momentum=True)
    row = next(r for r in rep.rows if r["type"] == "epicardium")
    assert row["status"] == "shifted" and row["max_abs_log_f"] > 0.0 and rep.as_dict()["standin_momentum"] is True


# ---------------------------------------------------------------------------------------------------------------------
# pseudo-task script (experiments/t1_scope/pseudo_task.py)


@pytest.fixture(scope="module")
def pt():
    """experiments/t1_scope/pseudo_task.py, loaded by path."""
    spec = importlib.util.spec_from_file_location("t1_scope_pseudo_task",
                                                  REPO_ROOT / "experiments" / "t1_scope" / "pseudo_task.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def test_pseudo_task_writes_nothing_when_the_procedure_changed(pt, tmp_path, monkeypatch):
    out = tmp_path / "pseudo"
    out.mkdir()
    (out / "scores.csv").write_text("run A\n")
    (out / "run.json").write_text('{"run": "A"}\n')
    at_start = {"sha256": {"src/vec/t1_scope.py": "a"}}
    files = {"scores.csv": "run B\n", "run.json": '{"run": "B"}\n'}
    monkeypatch.setattr(pt, "git_state", lambda: {"sha256": {"src/vec/t1_scope.py": "b"}})
    for d in (out, tmp_path / "first_run"):
        with pytest.raises(SystemExit, match="nothing was written"):
            pt.write_results(d, files, at_start)
    assert (out / "scores.csv").read_text() == "run A\n" and (out / "run.json").read_text() == '{"run": "A"}\n'
    assert not (tmp_path / "first_run").exists()
    monkeypatch.setattr(pt, "git_state", lambda: at_start)
    pt.write_results(out, files, at_start)
    assert (out / "scores.csv").read_text() == "run B\n" and (out / "run.json").read_text() == '{"run": "B"}\n'


def test_pseudo_task_refuses_too_few_contamination_cells_at_a_stage(pt):
    def stage(counts):
        a = ad.AnnData(X=np.zeros((sum(counts.values()), 1), np.float32))
        a.obs["celltype"] = [lab for lab, n in counts.items() for _ in range(n)]
        return a

    labels = ["Foregut", "Epithelium"]
    truth = {"truth": stage({"Epithelium": 1, "V-CM": 3})}  # its contamination is drawn from them: one cell will do
    endpoints = {"sibling": stage({"Foregut": 12, "Epithelium": 8}), "reference": stage({"Foregut": 20, "V-CM": 2})}
    assert pt.check_contamination(truth, endpoints, labels, 20) == {
        "truth": {"Foregut": 0, "Epithelium": 1}, "sibling": {"Foregut": 12, "Epithelium": 8},
        "reference": {"Foregut": 20, "Epithelium": 0}}
    endpoints["momentum start"] = stage({"Foregut": 19, "JCF": 5})  # surrounding would get no momentum, silently
    with pytest.raises(pt.typer.BadParameter, match=r"\] at momentum start \(19 < 20\);"):
        pt.check_contamination(truth, endpoints, labels, 20)
    with pytest.raises(pt.typer.BadParameter, match=r"\] at truth \(0 < 1\);"):
        pt.check_contamination({"truth": stage({"V-CM": 3})}, {}, labels, 20)
    with pytest.raises(pt.typer.BadParameter, match=r"\] at sibling \(0 < 1\);"):  # an endpoint always needs one
        pt.check_contamination(truth, {"sibling": stage({"V-CM": 3})}, labels, 0)


def test_pseudo_task_checks_every_momentum_endpoint_against_min_cells(pt, tmp_path, monkeypatch):
    """main() checks the contamination labels at load time, before the truth is drawn: every momentum endpoint is held
    to --min-cells and the truth library to one cell (synthetic files; nothing is scored or written)."""

    def libraries(*stages):  # obs-name suffixes are the T1 library tags (pseudo_task.library_of)
        a = ad.concat(stages)
        a.obs_names = [f"c{i}_{k + 1}" for k, s in enumerate(stages) for i in range(s.n_obs)]
        return a

    files = {"E8.5_RNA.h5ad": make_stage(60, sizes={**E85_SIZES, "Foregut": 19}),
             "E9.5_RNA.h5ad": make_stage(61, sizes={**E95_SIZES, "Foregut": 30}),
             "E10.5_RNA.h5ad": libraries(make_stage(62, sizes={**SIZES, "Foregut": 25}),
                                         make_stage(63, sizes={**SIZES, "Foregut": 5}))}

    class Reached(Exception):  # the run got past the check
        pass

    def make_truth(*args, **kwargs):
        raise Reached

    monkeypatch.setattr(pt, "load_stage", lambda path, lineage=True: files[Path(path).name].copy())
    monkeypatch.setattr(pt, "get_board", lambda key: t1_board(GENES))
    monkeypatch.setattr(pt, "git_state", lambda: {})
    monkeypatch.setattr(pt, "make_truth", make_truth)
    monkeypatch.setattr(pt, "RESULTS", tmp_path)
    run = pt.typer.main.get_command(pt.app).main
    args = ["--ref", "E9.5_RNA.h5ad", "--truth-stage", "E10.5_RNA.h5ad", "--momentum-start", "E8.5_RNA.h5ad",
            "--contamination-labels", "Foregut", "--n-cells", "100"]
    with pytest.raises(pt.typer.BadParameter, match=r"\] at momentum start E8\.5 \(19 < 20\);"):  # the default 20
        run(args, standalone_mode=False)
    with pytest.raises(Reached):  # the truth library's 5 cells are enough
        run([*args, "--min-cells", "19"], standalone_mode=False)
    with pytest.raises(pt.typer.BadParameter, match=r"\] at sibling E10\.5 libraries 1 \(25 < 31\), reference E9\.5 "
                                                    r"\(30 < 31\), momentum start E8\.5 \(19 < 31\);"):
        run([*args, "--min-cells", "31"], standalone_mode=False)
    with pytest.raises(pt.typer.BadParameter, match=r"\] at truth E10\.5 library 2 \(0 < 1\), sibling "):
        run([*args, "--contamination-labels", "Epithelium"], standalone_mode=False)
    plain = ["--ref", "E9.5_RNA.h5ad", "--truth-stage", "E10.5_RNA.h5ad", "--contamination-labels", "Foregut",
             "--n-cells", "100", "--min-cells", "26"]  # the default geometry, without a momentum start, is checked too
    with pytest.raises(pt.typer.BadParameter, match=r"\] at sibling E10\.5 libraries 1 \(25 < 26\);"):
        run(plain, standalone_mode=False)
