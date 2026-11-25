"""Synthetic tests of OIC125's G0 code (experiments/heart_ext/oic_core.py, oic_units.py, oic_g0.py): no data file, no
outside file, no target and no score. Planted: a per-gene platform efficiency (must cancel), a lineage-specific late
change (must be recovered), a typing that the reader confirms, and the gate's refusals."""

from __future__ import annotations

import gzip
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy import sparse

REPO = Path(__file__).resolve().parents[1]
EXT = REPO / "experiments" / "heart_ext"
sys.path.insert(0, str(EXT))
import fetch_oic  # noqa: E402
import oic_core as oc  # noqa: E402
import oic_units as ou  # noqa: E402


def _g0():
    spec = importlib.util.spec_from_file_location("oic_g0_t", EXT / "oic_g0.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


G = 200
GENES = [f"g{i:03d}" for i in range(G)]
MARK = {"cm": range(0, 20), "endocardium": range(20, 40), "epicardium": range(40, 60),
        "vascular_endothelium": range(60, 80), "blood": range(80, 90)}


def profiles(seed=0):
    rng = np.random.default_rng(seed)
    base = np.exp(rng.normal(0.0, 0.6, G)) * 5.0
    out = {}
    for k, idx in MARK.items():
        p = base.copy()
        p[list(idx)] *= 25.0
        p[90:] *= np.exp(rng.normal(0.0, 0.5, G - 90))  # the reference's own contrasts
        out[k] = p
    return out


def ref_cells(prof, rng, n=None):
    n = n or {"cm": 400, "endocardium": 250, "epicardium": 150, "vascular_endothelium": 250, "blood": 100}
    rows, cats = [], []
    for k, m in n.items():
        p = prof[k] / prof[k].sum()
        rows.append(rng.poisson(300 * p, (m, G)))
        cats += [k] * m
    return sparse.csr_matrix(np.vstack(rows)), np.array(cats)


class FakeRef:
    def __init__(self, name, C, cats, start_cells=None):
        self.name, self.genes, self.C, self.cats = name, list(GENES), sparse.csr_matrix(C), np.asarray(cats)
        self.ref_cells = {k: int((self.cats == k).sum()) for k in oc.CATS}
        self.start_cells = start_cells or {k: 100 for k in oc.CATS}
        self.record = {"source": "synthetic"}


def section(prof, late, eff, rng, size=100, heart=60, depth=300):
    """A section on a grid: a heart square (cm connected, the three other categories one unit in four, each next to cm
    units) amid background units; counts Poisson of each unit's late profile times the platform efficiency."""
    bx, by, cats = [], [], []
    off = (size - heart) // 2
    cyc = ["endocardium", "epicardium", "vascular_endothelium"]
    n = 0
    for i in range(size):
        for j in range(size):
            bx.append(j)
            by.append(i)
            hi, hj = i - off, j - off
            if 0 <= hi < heart and 0 <= hj < heart:
                if hi % 2 == 0 and hj % 2 == 0:
                    cats.append(cyc[n % 3])
                    n += 1
                else:
                    cats.append("cm")
            else:
                cats.append("bg")
    bg = np.exp(rng.normal(0, 0.3, G)) * 5.0
    rows = []
    for k in cats:
        p = (late[k] if k != "bg" else bg) * eff
        rows.append(rng.poisson(depth * p / p.sum()))
    M = sparse.csr_matrix(np.vstack(rows))
    return {"bx": np.array(bx), "by": np.array(by), "total": np.asarray(M.sum(1)).ravel() + 50, "panel": M,
            "genes": list(GENES), "held": list(GENES), "unit_um": 12.0, "cats": np.array(cats)}


def late_profiles(prof, d, shared):
    return {k: prof[k] * np.exp(d.get(k, 0.0) + shared) for k in prof}


@pytest.fixture(scope="module")
def world():
    return build_world()


def build_world(d_sigma: float = 0.6, seed: int = 7) -> dict:
    rng = np.random.default_rng(seed)
    prof = profiles()
    C, cats = ref_cells(prof, rng)
    ref = FakeRef("E10.5", C, cats)
    d = {k: np.r_[np.zeros(90), rng.normal(0, d_sigma, G - 90)] for k in oc.CATS}  # the lineage-specific change
    shared = rng.normal(0, 0.3, G)
    late = late_profiles(prof, d, shared)
    eff = {s: np.exp(rng.normal(0, 0.7, G)) for s in ("X021", "X022")}
    units = {"pup_P1": section(prof, late, eff["X021"], rng),
             "E15.5_E1S1": section(prof, late, eff["X022"], rng),
             "E16.5_E1S1": section(prof, late, eff["X022"], rng)}
    doc = {"sections": [{"section": "pup_P1", "study": "X021", "stage": "P1"},
                        {"section": "E15.5_E1S1", "study": "X022", "stage": "E15.5"},
                        {"section": "E16.5_E1S1", "study": "X022", "stage": "E16.5"}]}
    rk = oc.marker_ranks(C, cats, GENES)
    return {"ref": ref, "units": units, "doc": doc, "rk": rk, "d": d, "eff": eff}


def confirm_all(g0, w, ref="E10.5", verdict="heart"):
    out = {}
    secs = g0.study_sections(w["doc"])
    for s in g0.STUDIES:
        genes = g0.study_genes(w["units"], secs[s], w["ref"].genes)
        markers = oc.markers_for(w["rk"], genes)
        for r in secs[s]:
            u = w["units"][r["section"]]
            _, region = g0.region_of(u, genes, markers)
            out[(ref, r["section"])] = (verdict, g0.region_sha256(u, region))
    return out


# ---- core ---------------------------------------------------------------------------------------------------------
def test_auroc_and_markers_take_the_reference_alone(world):
    rk = world["rk"]
    assert set(rk["ranks"]) == set(oc.CATS)
    for k in oc.CATS:
        top = set(rk["ranks"][k][:20])
        assert top == {GENES[i] for i in MARK[k]}
    held = GENES[5:]
    m = oc.markers_for(rk, held)
    assert all(len(v) == oc.N_MARKERS for v in m.values()) and not ({"g000", "g004"} & set(m["cm"]))
    X = np.array([[0.0], [1.0], [2.0], [3.0]])
    assert oc.auroc_columns(X, np.array([False, False, True, True]))[0] == 1.0
    assert oc.auroc_columns(X, np.array([True, True, False, False]))[0] == 0.0


def test_unit_scores_cancel_a_per_gene_efficiency():
    rng = np.random.default_rng(3)
    C = rng.poisson(5, (50, 40))
    markers = {k: [f"h{i}" for i in range(j * 10, j * 10 + 10)] for j, k in enumerate(oc.CATS)}
    genes = [f"h{i}" for i in range(40)]
    a = oc.unit_scores(C, genes, markers, min_counts=1)
    eff = np.exp(rng.normal(0, 1, 40))
    b = oc.unit_scores(C * eff, genes, markers, min_counts=1)
    np.testing.assert_allclose(a["S"], b["S"], rtol=0.01)  # up to each unit's renormalisation


def test_heart_region_is_the_largest_component_dilated():
    bx, by = np.meshgrid(np.arange(30), np.arange(30))
    bx, by = bx.ravel(), by.ravel()
    seed = np.zeros(900, bool)
    big = (bx >= 5) & (bx < 15) & (by >= 5) & (by < 15)
    small = (bx >= 25) & (by >= 25)
    seed[big | small] = True
    reg = oc.heart_region(bx, by, seed, 12.0, 25.0)
    assert reg[big].all() and not reg[small].any()
    # 25 um at 12 um a unit: two units out along an axis, not (2, 1)
    assert reg[(bx == 3) & (by == 10)].all() and reg[(bx == 4) & (by == 4)].all()
    assert not reg[(bx == 2) & (by == 10)].any() and not reg[(bx == 3) & (by == 4)].any()
    assert not oc.heart_region(bx, by, np.zeros(900, bool), 12.0).any()


def test_contrasts_cancel_efficiency_and_centre_on_the_set():
    rng = np.random.default_rng(1)
    means = {k: np.exp(rng.normal(5, 1, 60)) for k in oc.CATS}
    eff = np.exp(rng.normal(0, 1, 60))
    a = oc.contrasts(means, list(oc.CATS), pseudo=0.0)
    b = oc.contrasts({k: v * eff for k, v in means.items()}, list(oc.CATS), pseudo=0.0)
    for k in oc.CATS:
        np.testing.assert_allclose(a[k], b[k], atol=1e-10)
        assert abs(np.median(a[k])) < 1e-12
    three = oc.contrasts(means, ["cm", "endocardium", "vascular_endothelium"])
    assert set(three) == {"cm", "endocardium", "vascular_endothelium"}


def test_residual_removes_the_reference_term_and_guards():
    rng = np.random.default_rng(2)
    cref = rng.normal(0, 1, 100)
    r0 = rng.normal(0, 0.3, 100)
    csrc = 0.2 + 0.6 * cref + r0
    guard = np.ones(100, bool)
    guard[:10] = False
    r, fit = oc.residual(csrc, cref, guard)
    assert np.isnan(r[:10]).all() and abs(fit["b"] - 0.6) < 0.15
    assert abs(np.corrcoef(r[10:], cref[10:])[0, 1]) < 1e-9


def test_consensus_needs_two_studies_one_sign_no_typing_gene_and_clips():
    a = np.array([0.5, 0.5, -0.2, 0.3, np.nan, 2.0, 0.4])
    b = np.array([0.3, -0.1, -0.4, 0.3, 0.2, 3.0, 0.4])
    ex = np.array([False, False, False, True, False, False, False])
    P, rec = oc.consensus({"X021": a, "X022": b}, ex)
    assert P[0] == pytest.approx(0.4) and np.isnan(P[1]) and P[2] == pytest.approx(-0.3)
    assert np.isnan(P[3]) and np.isnan(P[4]) and P[5] == pytest.approx(np.log(2)) and rec["clipped"] == 1


def test_spearman_test_passes_shared_residuals_and_fails_noise():
    rng = np.random.default_rng(5)
    x = rng.normal(0, 1, 120)
    level = rng.normal(0, 1, 120)
    good = oc.spearman_test(x + rng.normal(0, 0.5, 120), x + rng.normal(0, 0.5, 120), level, "cm", n_perm=200)
    bad = oc.spearman_test(rng.normal(0, 1, 120), rng.normal(0, 1, 120), level, "cm", n_perm=200)
    assert good["pass"] and not bad["pass"]


def test_derangement_has_no_fixed_point_and_stays_in_its_bin():
    level = np.arange(57, dtype=float)
    for seed in (20261126, 20261127, 20261128):
        p = oc.derange_within_bins(level, "cm", seed)
        assert sorted(p) == list(range(57)) and not (p == np.arange(57)).any()
        bins = np.array_split(np.argsort(level, kind="stable"), 10)
        where = {int(i): b for b, idx in enumerate(bins) for i in idx}
        assert all(where[int(i)] == where[int(p[i])] for i in range(57))
    assert not (oc.derange_within_bins(np.arange(3.0), "x", 1) == np.arange(3)).any()


# ---- the end to end on synthetic sections -------------------------------------------------------------------------
def test_chords_recover_the_planted_change_through_two_platforms(world):
    g0 = _g0()
    rec = g0.chords_at(world["ref"], world["units"], world["doc"], world["rk"], confirm_all(g0, world))
    assert rec["pass"], rec["why"]
    for s in g0.STUDIES:
        assert rec["sources"][s]["pass"] and rec["sources"][s]["agreement"]["share"] >= oc.AGREE_ALL
    for k in rec["P"]:
        P = np.array([np.nan if x is None else x for x in rec["P"][k]])
        held = np.isfinite(P)
        assert held.sum() > 20
        assert not held[:90][[i for i in range(90) if GENES[i] in oc.typing_genes(
            rec["sources"]["X021"]["markers"])]].any()
        d = world["d"][k][held]
        assert np.corrcoef(P[held], d - d.mean())[0, 1] > 0.6, k
    assert rec["n_big"] >= oc.N_P_MIN and rec["lineage_tests"]["cm"]["pass"]


def test_a_doubt_or_a_region_the_reader_did_not_see_drops_the_section(world):
    g0 = _g0()
    v = confirm_all(g0, world)
    v[("E10.5", "pup_P1")] = ("doubt", v[("E10.5", "pup_P1")][1])
    rec = g0.chords_at(world["ref"], world["units"], world["doc"], world["rk"], v)
    assert not rec["pass"] and any("fewer than two sources" in w for w in rec["why"])
    v = confirm_all(g0, world)
    v[("E10.5", "pup_P1")] = ("heart", "0" * 64)
    assert not g0.chords_at(world["ref"], world["units"], world["doc"], world["rk"], v)["pass"]


def test_a_category_with_too_few_start_cells_gets_no_chord(world):
    g0 = _g0()
    ref = FakeRef("E10.5", world["ref"].C, world["ref"].cats,
                  start_cells={"cm": 100, "endocardium": 100, "epicardium": 3, "vascular_endothelium": 100})
    rec = g0.chords_at(ref, world["units"], world["doc"], world["rk"], confirm_all(g0, world))
    assert "epicardium" not in rec["P"] and "cm" in rec["P"]


def test_idle_chords_stop(world):
    g0 = _g0()
    flat = {k: v.copy() for k, v in world["units"].items()}
    rng = np.random.default_rng(11)
    prof = profiles()
    for name in flat:  # late sources with no lineage-specific change at all: the residuals are noise
        eff = np.exp(rng.normal(0, 0.7, G))
        flat[name] = section(prof, late_profiles(prof, {}, 0.0), eff, rng, depth=3000)
    rec = g0.chords_at(world["ref"], flat, world["doc"], world["rk"],
                       confirm_all(g0, {**world, "units": flat}))
    assert not rec["pass"]
    assert any(w.startswith("(c)") or w.startswith("(d)") for w in rec["why"])


def test_chord_file_round_trip_and_sha(world, tmp_path):
    g0 = _g0()
    rec = g0.chords_at(world["ref"], world["units"], world["doc"], world["rk"], confirm_all(g0, world))
    doc = {"row": "H38", "file": "OIC125", "ref": "E10.5", "genes": rec["genes"], "P": rec["P"]}
    p = tmp_path / "chords_E10.5.json"
    p.write_text(json.dumps(doc))
    sha = hashlib.sha256(p.read_bytes()).hexdigest()
    other = list(reversed(GENES))
    got = g0.load_chords(p, sha, other)
    for k, v in rec["P"].items():
        a = np.array([np.nan if x is None else x for x in v])
        np.testing.assert_array_equal(got[k][::-1], a)
    with pytest.raises(ou.Refused):
        g0.load_chords(p, "0" * 64)


# ---- units --------------------------------------------------------------------------------------------------------
def _gem(path: Path, rows, header=True):
    with gzip.open(path, "wt") as fh:
        fh.write("#FileFormat=GEMv0.1\n#BinSize=1\n")
        if header:
            fh.write("geneID\tx\ty\tMIDCounts\tExonCount\n")
        for g, x, y, c in rows:
            fh.write(f"{g}\t{x}\t{y}\t{c}\t{c}\n")


def test_mosta_units_bin24_and_sha(tmp_path):
    rows = [("g001", 0, 0, 2), ("g001", 23, 23, 1), ("g002", 24, 0, 5), ("Other", 47, 1, 7), ("g001", 50, 30, 1)]
    p = tmp_path / "E15.5_E1S1_GEM_bin1.tsv.gz"
    _gem(p, rows)
    sha = hashlib.sha256(p.read_bytes()).hexdigest()
    u = ou.mosta_units(p, ["g001", "g002", "g003"], sha256=sha, chunk_rows=2)
    got = {(int(a), int(b)): int(t) for a, b, t in zip(u["bx"], u["by"], u["total"], strict=True)}
    assert got == {(0, 0): 3, (1, 0): 12, (2, 1): 1}
    M = u["panel"].toarray()
    i = {(int(a), int(b)): n for n, (a, b) in enumerate(zip(u["bx"], u["by"], strict=True))}
    assert M[i[(0, 0)]].tolist() == [3, 0, 0] and M[i[(1, 0)]].tolist() == [0, 5, 0]
    assert u["held"] == ["g001", "g002"] and u["unit_um"] == 12.0 and u["sha256"] == sha
    with pytest.raises(ou.Refused):
        ou.mosta_units(p, ["g001"], sha256="0" * 64)
    q = tmp_path / "bad.tsv.gz"
    with gzip.open(q, "wt") as fh:
        fh.write("geneID\tx\ty\tfoo\n")
    with pytest.raises(ou.Refused):
        ou.mosta_units(q, ["g001"])


def test_xenium_units_read_only_allowed_columns_and_filter():
    batch = {"feature_name": np.array([b"g001", b"g002", b"NegControlProbe_1", b"g001", b"g003"], dtype=object),
             "x_location": np.array([1.0, 13.0, 2.0, 5.0, 30.0]), "y_location": np.array([1.0, 1.0, 1.0, 11.9, 25.0]),
             "qv": np.array([30.0, 30.0, 30.0, 10.0, 25.0]), "is_gene": np.array([True, True, False, True, True])}
    u = ou.xenium_units([batch], ["g001", "g002"])
    got = {(int(a), int(b)): int(t) for a, b, t in zip(u["bx"], u["by"], u["total"], strict=True)}
    assert got == {(0, 0): 1, (1, 0): 1, (2, 2): 1}
    assert u["held"] == ["g001", "g002"]
    del batch["is_gene"]
    u2 = ou.xenium_units([batch], ["g001", "g002"])
    assert int(u2["total"].sum()) == 3  # the control name is dropped by its prefix
    pq = pytest.importorskip("pyarrow.parquet")
    import pyarrow as pa

    batch["cell_id"] = np.array(["a", "b", "c", "d", "e"], dtype=object)
    t = pa.table({k: (v.astype(str) if v.dtype == object else v) for k, v in batch.items()})
    path = Path(pytest.importorskip("tempfile").mkdtemp()) / "t.parquet"
    pq.write_table(t, path)
    seen = list(ou.xenium_batches(path))
    assert "cell_id" not in seen[0] and set(seen[0]) == set(ou.XENIUM_COLUMNS)


def test_manifest_is_checked_against_the_pins(tmp_path):
    e = next(x for x in fetch_oic.PINNED if x["kind"] == "file")
    f = tmp_path / e["name"]
    f.write_bytes(b"x")
    man = {"files": [{"name": e["name"], "url": e["url"], "size": e["bytes"], "sha256": "a" * 64}]}
    (tmp_path / "MANIFEST.json").write_text(json.dumps(man))
    with pytest.raises(ou.Refused):  # the file on disk is not of its recorded size
        ou.manifest_files(tmp_path)
    man["files"][0]["size"] = 1
    (tmp_path / "MANIFEST.json").write_text(json.dumps(man))
    with pytest.raises(ou.Refused):  # the record's size is not the pin
        ou.manifest_files(tmp_path)
    man["files"][0].update(size=e["bytes"], url="https://example.org/x")
    (tmp_path / "MANIFEST.json").write_text(json.dumps(man))
    with pytest.raises(ou.Refused):
        ou.manifest_files(tmp_path)


def test_units_round_trip(tmp_path):
    u = {"bx": np.array([0, 1]), "by": np.array([0, 0]), "total": np.array([3, 4]),
         "panel": sparse.csr_matrix(np.array([[1, 0], [0, 2]])), "genes": ["a", "b"], "held": ["a", "b"],
         "unit_um": 12.0, "sha256": "f" * 64, "source": "X022"}
    sha = ou.save_units(u, tmp_path / "s.npz")
    assert sha == hashlib.sha256((tmp_path / "s.npz").read_bytes()).hexdigest()
    v = ou.load_units(tmp_path / "s.npz")
    assert (v["panel"].toarray() == u["panel"].toarray()).all() and v["genes"] == ["a", "b"]
    with pytest.raises(FileExistsError):
        _g0().write_new(tmp_path / "s.json", "x") or _g0().write_new(tmp_path / "s.json", "y")


def test_plot_has_no_text_and_is_new(world, tmp_path):
    g0 = _g0()
    u = world["units"]["pup_P1"]
    region = np.zeros(len(u["bx"]), bool)
    region[:50] = True
    sha = g0.plot_region(u, region, tmp_path / "p.png")
    assert sha == hashlib.sha256((tmp_path / "p.png").read_bytes()).hexdigest()
    with pytest.raises(ou.Refused):
        g0.plot_region(u, region, tmp_path / "p.png")


def test_regions_and_g0_steps_end_to_end(tmp_path, monkeypatch):
    """The regions step's plots and reader sheet, then the g0 step with the references and the reader's file stood in:
    the chord files, g0.json, g0.txt and the statement row H38's module reads; nothing is overwritten."""
    sys.path.insert(0, str(REPO / "experiments" / "harness"))
    import h38_row

    g0 = _g0()
    world = build_world(d_sigma=1.5)  # a change large enough to clear (d) at the practice strengths too
    work = tmp_path / "work"
    (work / "units").mkdir(parents=True)
    stage = {"pup_P1": "P1", "E15.5_E1S1": "E15.5", "E16.5_E1S1": "E16.5"}
    secs = []
    for name, u in world["units"].items():
        sha = ou.save_units({**u, "sha256": "f" * 64, "source": "X"}, work / "units" / f"{name}.npz")
        secs.append({"section": name, "study": "X021" if name == "pup_P1" else "X022", "stage": stage[name],
                     "file": f"{name}.gz", "source_sha256": "f" * 64, "npz": f"{name}.npz", "npz_sha256": sha})
    (work / "units.json").write_text(json.dumps({"code_sha256": g0.code_sha256(), "sections": secs}))
    refs_doc = {"refs": {r: world["rk"] for r in g0.REFS}}
    (work / "ranks.json").write_text(json.dumps(refs_doc))
    assert g0.run_regions(work, work) == 0
    sheet = json.loads((work / "reader_sheet.json").read_text())
    assert len(sheet["verdicts"]) == 3 * len(world["units"]) and all(v["verdict"] is None for v in sheet["verdicts"])
    for v in sheet["verdicts"]:
        assert hashlib.sha256((work / v["plot"]).read_bytes()).hexdigest() == v["plot_sha256"]
    refs = {r: FakeRef(r, world["ref"].C, world["ref"].cats) for r in g0.REFS}
    monkeypatch.setattr(g0, "utc_now", lambda: g0.G0_CUTOFF_EPOCH - 3600.0)  # never the real clock
    monkeypatch.setattr(g0, "load_refs", lambda data_dir: refs)
    monkeypatch.setattr(g0, "ranks_doc", lambda _refs: json.loads(json.dumps(refs_doc)))
    verdicts = {(v["ref"], v["section"]): ("heart", v["region_sha256"]) for v in sheet["verdicts"]}
    monkeypatch.setattr(g0, "reader_verdicts", lambda path: verdicts)
    reader = tmp_path / "reader.json"
    reader.write_text("{}")
    results = tmp_path / "results"
    assert g0.run_g0(work, work, None, reader=reader, results=results) == 0
    line = (results / "g0.txt").read_text().strip()
    m = h38_row.G0_FORM.fullmatch(line)
    assert m
    for ref, sha in zip(g0.REFS, m.groups(), strict=True):
        assert hashlib.sha256((results / f"chords_{ref}.json").read_bytes()).hexdigest() == sha
        doc = json.loads((results / f"chords_{ref}.json").read_text())
        assert doc["ref"] == ref and doc["sources"] == {f"{n}.gz": "f" * 64 for n in world["units"]}
        assert set(doc["centre"]) == {"X021", "X022"} and "cm" in doc["lineages_passed"]  # FIX 5: the set reported
    assert json.loads((results / "g0.json").read_text())["verdict"] == "PASS"
    with pytest.raises(FileExistsError):
        g0.run_g0(work, work, None, reader=reader, results=results)
    # a practice reference whose reader doubted every section: VOID (exit 5), no chord file written
    v2 = {k: (("doubt" if k[0] == "E8.75" else "heart"), s) for k, (_, s) in verdicts.items()}
    monkeypatch.setattr(g0, "reader_verdicts", lambda path: v2)
    r2 = tmp_path / "results2"
    assert g0.run_g0(work, work, None, reader=reader, results=r2) == 5
    assert not list(r2.glob("chords_*.json"))
    assert (r2 / "g0.txt").read_text().startswith("G0 H38: VOID AT PRACTICE")
    # the same PASS, ended after 2026-11-26 03:00 UTC: STOP (exit 3), no chord file
    monkeypatch.setattr(g0, "reader_verdicts", lambda path: verdicts)
    monkeypatch.setattr(g0, "utc_now", lambda: g0.G0_CUTOFF_EPOCH + 1.0)
    r3 = tmp_path / "results3"
    assert g0.run_g0(work, work, None, reader=reader, results=r3) == 3
    assert not list(r3.glob("chords_*.json"))
    assert "cut-off" in (r3 / "g0.txt").read_text() and json.loads((r3 / "g0.json").read_text())["verdict"] == "STOP"



# ---- G0(c) against the review of 38eabe2 (blocker 1): no late change must fail, a shared late change must pass -------
G3 = 300
GENES3 = [f"h{i:03d}" for i in range(G3)]
MARK3 = {"cm": range(0, 20), "endocardium": range(20, 40), "epicardium": range(40, 60),
         "vascular_endothelium": range(60, 80)}


def no_change_world(seed: int, change: float = 0.0, ref_n=(2000, 400, 200, 400)):
    """The reviewer's control: a sampled reference (2000/400/200/400 cells) and two studies whose units are drawn,
    with true labels and each study's own per-gene efficiency, from the reference's own profiles; `change` adds one
    lineage-specific late change shared by both studies (0: no late change at all)."""
    rng = np.random.default_rng(seed)
    base = np.exp(rng.normal(0, 0.8, G3)) * 3
    prof = {}
    for k, idx in MARK3.items():
        p = base.copy()
        p[list(idx)] *= 20
        p[80:] *= np.exp(rng.normal(0, 0.5, G3 - 80))
        prof[k] = p
    rows, cats = [], []
    for k, n in zip(oc.CATS, ref_n, strict=True):
        rows.append(rng.poisson(200 * prof[k] / prof[k].sum(), (n, G3)))
        cats += [k] * n
    d = {k: np.r_[np.zeros(80), rng.normal(0, 0.5, G3 - 80)] * change for k in oc.CATS}

    def stage_means(eff, n=(3000, 800, 400, 800)):
        out = {}
        for k, m in zip(oc.CATS, n, strict=True):
            p = prof[k] * np.exp(d[k]) * eff
            out[k] = oc.cp10k(rng.poisson(150 * p / p.sum(), (m, G3))).mean(0)
        return out
    studies = [{"name": "X021", "genes": GENES3, "centre": list(oc.CATS),
                "stages": {"P1": stage_means(np.exp(rng.normal(0, 0.7, G3)))}},
               {"name": "X022", "genes": GENES3, "centre": list(oc.CATS),
                "stages": {"E15.5": stage_means(np.exp(rng.normal(0, 0.7, G3))),
                           "E16.5": stage_means(np.exp(rng.normal(0, 0.7, G3)))}}]
    return sparse.csr_matrix(np.vstack(rows)), np.array(cats), studies


TYPING3 = np.r_[np.ones(80, bool), np.zeros(G3 - 80, bool)]


@pytest.mark.parametrize("seed", range(6))
def test_g0c_fails_with_no_late_change(seed):
    refC, cats, studies = no_change_world(seed)
    t = oc.g0c_tests(refC, cats, GENES3, studies, TYPING3, "E10.5", n_perm=200)["tests"]
    assert not t["cm"]["pass"] and abs(t["cm"]["rho"]) < oc.SPEARMAN_MIN
    assert not any(v["pass"] for v in t.values())


def test_the_shared_reference_alone_made_the_old_test_pass():
    """The fault the review found, kept as a record: least squares of both studies on one estimated reference
    contrast agrees with no late change at all."""
    refC, cats, studies = no_change_world(0)
    code = np.array([oc.CATS.index(c) for c in cats])
    whole = oc.category_means(refC, code)
    cref = oc.contrasts(whole, list(oc.CATS))
    res = []
    for st in studies:
        per = []
        for m in st["stages"].values():
            cs = oc.contrasts(m, list(oc.CATS))
            per.append(oc.residual(cs["cm"], cref["cm"], (m["cm"] >= 5) & (whole["cm"] >= 5))[0])
        res.append(oc.study_mean(per))
    assert oc.spearman_test(res[0], res[1], np.log(whole["cm"] + 0.5), "cm", n_perm=200)["pass"]


@pytest.mark.parametrize("seed", range(3))
def test_g0c_passes_a_shared_late_change(seed):
    refC, cats, studies = no_change_world(seed, change=1.0)
    t = oc.g0c_tests(refC, cats, GENES3, studies, TYPING3, "E10.5", n_perm=200)["tests"]
    assert t["cm"]["pass"] and t["cm"]["rho"] >= oc.SPEARMAN_MIN


def test_g0c_keeps_the_typing_genes_out_and_its_halves_are_fixed():
    refC, cats, studies = no_change_world(1, change=1.0)
    a = oc.g0c_tests(refC, cats, GENES3, studies, TYPING3, "E10.5", n_perm=50)["tests"]["cm"]
    assert a["n"] <= G3 - 80
    h1, h2 = oc.ref_halves(cats, "E10.5"), oc.ref_halves(cats, "E10.5")
    assert np.array_equal(h1, h2) and not np.array_equal(h1, oc.ref_halves(cats, "E9.5"))
    for k in oc.CATS:
        m = cats == k
        assert abs(int(h1[m].sum()) - int((1 - h1[m]).sum())) <= 1
    r, fit = oc.iv_residual(np.arange(20.0), np.arange(20.0), -np.arange(20.0), np.ones(20, bool))
    assert fit["weak"] and np.isnan(r).all()  # halves that disagree: no residual, the lineage fails


def test_every_g0_step_refuses_without_row_h38_head_and_a_clean_tree(tmp_path):
    g0 = _g0()
    why = g0.refusals("0" * 40)
    assert any("HEAD is not" in w for w in why)
    text = (REPO / "experiments" / "harness" / "TRIES.md").read_text(encoding="utf-8")
    if not any(ln.startswith("| H38 ") for ln in text.splitlines()):
        assert any("no row H38" in w for w in why)
    assert g0.main(["ranks", "--out", str(tmp_path / "w"), "--commit", "0" * 40]) == 2
    assert not (tmp_path / "w").exists()


# ---- the GEM files' forms on the DGX (2026-11-25): headers repeated mid-file, comment lines, MIDCount or MIDCounts --
GEM_ROWS = [("g001", 0, 0, 2), ("g001", 23, 23, 1), ("g002", 24, 0, 5), ("Other", 47, 1, 7), ("g001", 50, 30, 1),
            ("g002", 51, 31, 3), ("g003", 100, 3, 4)]
COMMENTS = ["#FileFormat=GEMv0.1", "#SortedBy=None", "#BinType=Bin", "#BinSize=1", "#OffsetX=0", "#OffsetY=0"]


def _gem_lines(path: Path, lines: list[str]) -> str:
    with gzip.open(path, "wt") as fh:
        fh.write("\n".join(lines) + "\n")
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _data(rows=GEM_ROWS):
    return [f"{g}\t{x}\t{y}\t{c}" for g, x, y, c in rows]


def _units_of(u):
    M = u["panel"].toarray()
    return {(int(a), int(b)): (int(t), tuple(int(v) for v in M[i]))
            for i, (a, b, t) in enumerate(zip(u["bx"], u["by"], u["total"], strict=True))}


@pytest.mark.parametrize("block", [7, 64, 1 << 24])
@pytest.mark.parametrize("form", ["plain_midcounts", "comments_midcount", "plain_midcount", "repeated"])
def test_gem_forms_all_give_the_same_units(tmp_path, form, block):
    panel = ["g001", "g002", "g003"]
    ref_path = tmp_path / "ref.tsv.gz"
    _gem_lines(ref_path, ["geneID\tx\ty\tMIDCounts", *_data()])
    want = _units_of(ou.mosta_units(ref_path, panel, chunk_rows=3))
    d = _data()
    lines = {"plain_midcounts": ["geneID\tx\ty\tMIDCounts", *d],
             "comments_midcount": [*COMMENTS, "geneID\tx\ty\tMIDCount", *d],
             "plain_midcount": ["geneID\tx\ty\tMIDCount", *d],
             "repeated": ["geneID\tx\ty\tMIDCounts", *d[:3], "geneID\tx\ty\tMIDCounts", *d[3:5], "#a comment",
                          "geneID\tx\ty\tMIDCount", *d[5:]]}[form]
    p = tmp_path / f"{form}.tsv.gz"
    sha = _gem_lines(p, lines)
    u = ou.mosta_units(p, panel, sha256=sha, chunk_rows=2, block=block)
    assert _units_of(u) == want and u["held"] == ["g001", "g002", "g003"]
    n_comments = {"plain_midcounts": 0, "comments_midcount": 6, "plain_midcount": 0, "repeated": 1}[form]
    n_headers = 3 if form == "repeated" else 1
    assert u["comment_lines"] == n_comments and u["header_lines"] == n_headers


@pytest.mark.parametrize("lines,why", [
    (["geneID\tx\ty\tMIDCounts", "g001\t1\t2\t3", "g001\tabc\t2\t3"], "not data"),          # a non-numeric field
    (["geneID\tx\ty\tMIDCounts", "g001\t1\t2\t3", "g001\t1\t2"], "not data"),               # a missing field
    (["geneID\tx\ty\tMIDCounts", "g001\t1\t2\t3", "geneID\tx\ty\tExonCount"], "other columns"),  # another header
    (["g001\t1\t2\t3", "geneID\tx\ty\tMIDCounts"], "before the first header"),              # data first
    (["#only a comment"], "no header"),
    (["geneID\tx\ty\tUMI", "g001\t1\t2\t3"], "lack"),                                       # no count column
])
@pytest.mark.parametrize("block", [5, 1 << 24])
def test_gem_refuses_anything_else(tmp_path, lines, why, block):
    p = tmp_path / "bad.tsv.gz"
    _gem_lines(p, lines)
    with pytest.raises(ou.Refused, match=why):
        ou.mosta_units(p, ["g001"], chunk_rows=2, block=block)
