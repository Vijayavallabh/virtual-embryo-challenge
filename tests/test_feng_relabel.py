"""The permitted-only relabelling of Feng's cells (experiments/external/feng_relabel.py and feng_relabel.yaml), its
selector (feng_dirs groups, build_test.py t1 --feng-groups, the harness's late-proxy option) and row N54's listing, on
synthetic arrays and hand-made files only. No Feng file, released file, board or score is read; the one committed
table read here is the released-data marker evidence the yaml cites (experiments/t1_scope/results/label_evidence/
markers.csv, aggregates of the released E9.5 file)."""

from __future__ import annotations

import ast
import copy
import csv
import json
import re
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
import yaml
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
for _p in ("experiments/external", "experiments/t1_m1", "experiments/t1_direction", "experiments/t1_hygiene",
           "experiments/builds"):
    sys.path.insert(0, str(ROOT / _p))

import feng_dirs  # noqa: E402
import feng_relabel as fr  # noqa: E402
import filterio as fio  # noqa: E402
import n54_relabel_listing as n54  # noqa: E402
import test_t1_feng_groups_none  # noqa: E402, F401  (imported here, before any fixture: its REAL_FENG_DIRECTIONS)
from test_hidden_counts import UNION, _reference_file, t1_env  # noqa: E402, F401  (the synthetic build_t1 fixture)

from vec import harness  # noqa: E402

COMMITTED = fr.load_config()
N55_ = {"outcome": "kept", "setting": {"min_score": 0.5, "min_margin": 0.25}, "commit": "n55c", "ruling": None}
# what feng_relabel.preflight returns on a clean, committed, ruled tree (synthetic: the gate is never run here on data)
GATE_ = {"row": "N54", "ruling": "2026-11-20 (a)", "commit": "abc", "clean": True, "config_at_head": True,
         "n55": N55_}
GATE = dict(GATE_)
TARGET = ["--g-mom", "--base", "lib3"]  # the first upload's command, the one design either route builds


@pytest.fixture
def record_ok(monkeypatch):
    """feng_dirs.check_recomputed's last check (the manifest committed at HEAD as the labeller's RECORD) passes."""
    monkeypatch.setattr(feng_dirs, "record_matches", lambda p: True)

# ------------------------------------------------------------------ a synthetic configuration and data set

SETS = {
    "cm": {"role": "target", "label": "r_cm", "genes": {"Cm1": "syn", "Cm2": "syn", "Cm3": "syn"}},
    "ec": {"role": "target", "label": "r_ec", "genes": {"Ec1": "syn", "Ec2": "syn", "Ec3": "syn"}},
    "epi": {"role": "target", "label": "r_epicardial", "genes": {"Ep1": "syn", "Ep2": "syn", "Ep3": "syn"}},
    "endo_split": {"role": "split", "genes": {"En1": "syn", "En2": "syn"}},
    "vasc_split": {"role": "split", "genes": {"Va1": "syn", "Va2": "syn"}},
    "blood": {"role": "competing", "label": "r_blood", "genes": {"Bl1": "syn", "Bl2": "syn"}},
}
FILLER = [f"f{i:04d}" for i in range(1500)]
GENES = [g for s in SETS.values() for g in s["genes"]] + FILLER
MARKERS = {g for s in SETS.values() for g in s["genes"]}
BASE = np.where([g in MARKERS for g in GENES], 0.2, np.random.default_rng(99).gamma(0.5, 4.0, len(GENES)))
# cell type -> the marker genes it expresses highly
PROFILE = {"cm": ["Cm1", "Cm2", "Cm3"], "endo": ["Ec1", "Ec2", "Ec3", "En1", "En2"],
           "vasc": ["Ec1", "Ec2", "Ec3", "Va1", "Va2"], "ambig": ["Ec1", "Ec2", "Ec3", "En1", "En2", "Va1", "Va2"],
           "epi": ["Ep1", "Ep2", "Ep3"], "blood": ["Bl1", "Bl2"], "none": [],
           "doublet": ["Cm1", "Cm2", "Cm3", "Ep1", "Ep2", "Ep3"]}


def syn_cfg(**over) -> dict:
    cfg = copy.deepcopy({k: v for k, v in COMMITTED.items() if k not in ("sets", "sources")})
    cfg["sets"] = copy.deepcopy(SETS)
    cfg["sources"] = {"syn": "synthetic"}
    cfg.update(over)
    fr.check_config(cfg)
    return cfg


def counts(types: list[str], seed: int = 0, high: float = 30.0) -> sparse.csr_matrix:
    """Poisson counts: every gene at its base rate (spread like real genes, markers low), a type's markers high."""
    rng = np.random.default_rng(seed)
    rate = np.tile(BASE, (len(types), 1))
    pos = {g: i for i, g in enumerate(GENES)}
    for i, t in enumerate(types):
        for g in PROFILE[t]:
            rate[i, pos[g]] = high
    C = rng.poisson(rate).astype(np.int32)
    C[:, -1] += 1  # no empty cell
    return sparse.csr_matrix(C)


def cells(spec: dict[str, int]) -> list[str]:
    return [t for t, n in spec.items() for _ in range(n)]


SPEC = {"cm": 30, "endo": 15, "vasc": 15, "ambig": 6, "epi": 14, "blood": 8, "none": 6, "doublet": 4}


# ------------------------------------------------------------------ the committed yaml


def test_the_committed_yaml_is_well_formed_and_not_runnable_while_draft(tmp_path):
    cfg = COMMITTED
    assert cfg["method"] == "c" and cfg["stages"] == ["E9.5", "E14.5"]
    assert tuple(cfg["obs_columns"]) == fr.OBS_WHITELIST
    assert set(cfg["groups"]) == set(harness.LP_GROUPS)
    assert cfg["groups"]["Endocardium"] == ["r_endocardial_ec"]
    assert "r_ec_ambiguous" in cfg["groups"]["Endothelium"] and "r_ec_ambiguous" not in cfg["groups"]["Endocardium"]
    assert cfg["fallback"] == {"CM": "none", "Endocardium": "Endothelium", "Endothelium": "none", "Epicardium": "none"}
    assert cfg["unassigned_cap"] == 0.40 and cfg["assign"]["min_score"] == 0.5 and cfg["assign"]["min_margin"] == 0.25
    assert cfg["assign"]["split_margin"] == 0.1
    fr.check_min_cells(cfg, feng_dirs.load_config())
    # since 2026-11-19 no gene is UNSOURCED; the yaml stays draft until the captain's sign-off, and the code refuses it
    assert not fr.unsourced_genes(cfg)
    assert cfg["status"] == "draft"
    with pytest.raises(fr.Refused, match="status is 'draft'"):
        fr.preflight(fr.CONFIG, "2026-11-20 (a)", tmp_path / "out")  # the real gate, on the real committed yaml
    with pytest.raises(fr.Refused, match="status is 'draft'"):
        fr.run(cfg, fr.CONFIG, Path("/nonexistent"), Path("/nonexistent"), GATE)
    with pytest.raises(fr.Refused, match="no gate"):
        fr.run(cfg, fr.CONFIG, Path("/nonexistent"), Path("/nonexistent"), None)
    assert not (tmp_path / "out").exists()
    # decisions (1) to (3) of 2026-11-19: no PENDING choice is left; lib3, g_mom on under none, a stop takes none
    assert fr.pending_choices(cfg) == []
    assert cfg["terminal"] == {"requires_direction": ["CM", "Endocardium"], "otherwise": "none",
                               "none_with_g_mom": True, "base": "lib3", "on_stop": "none"}
    # decision (4): N55's rule, its first setting the one under assign
    assert cfg["n55_rule"]["order"] == [{"min_score": 0.5, "min_margin": 0.25}, {"min_score": 0.3, "min_margin": 0.25},
                                        {"min_score": 0.5, "min_margin": 0.15}]
    assert cfg["n55_rule"]["chamber_labels"]["CM"]["released"] == cfg["sets"]["cm"]["released_labels"]
    # decision (5): the output folder is inside data/derived, never data/external or data/raw
    assert cfg["output"]["dir"] == "data/derived/feng_relabel_c"
    assert "R2" in cfg["r2"] and "selects nothing" in cfg["authors_labels"]
    import pandas as pd

    labels = set(pd.read_csv(ROOT / "experiments/t1_scope/results/label_evidence/markers.csv").query(
        "stage == 'E9.5'")["label"])
    covered = {lab for s in cfg["sets"].values() if s["role"] != "split" for lab in s.get("released_labels", [])}
    assert labels == covered
    # no Feng author label appears as a recomputed label or group member
    for labs in cfg["groups"].values():
        assert all(lab.startswith("r_") for lab in labs)


def _panel() -> dict[str, list[str]]:
    src = (ROOT / "experiments/t1_scope/label_evidence.py").read_text()
    start = src.index("MARKERS = {")
    return ast.literal_eval(src[start + len("MARKERS = "):src.index("}\n", start) + 1])


def test_every_released_data_gene_passes_the_yamls_rule_and_every_panel_gene_left_out_fails_it():
    """sources.released_e95: a panel gene joins a class if its highest released E9.5 mean over the class's labels is at
    least 0.3 and at least twice its highest mean over the other labels (for a split set, the other split label)."""
    import pandas as pd

    d = pd.read_csv(ROOT / "experiments/t1_scope/results/label_evidence/markers.csv")
    d = d[d.stage == "E9.5"].set_index("label")
    panel = _panel()
    checked = 0
    used_by_groups = {p for s in COMMITTED["sets"].values() if s["role"] != "competing" for p in s.get("panel", [])}
    for name, s in COMMITTED["sets"].items():
        if "panel" not in s:  # the late classes rest on literature sources alone
            assert all(sid.startswith("lit_") for src in s["genes"].values() for sid in fr.gene_sources(src)), name
            continue
        if s["role"] == "competing":  # a competing class never draws on a target's or split's panel group
            assert not set(s["panel"]) & used_by_groups, name
        own = s["released_labels"]
        other = s.get("against") or [lab for lab in d.index if lab not in own]
        excluded = s.get("excluded") or {}
        for g in [g for p in s["panel"] for g in panel[p]]:
            a, b = d.loc[own, f"mean:{g}"].max(), d.loc[other, f"mean:{g}"].max()
            passes = a >= 0.3 and a >= 2 * b
            # a gene that passes is a member, or is left out by name with its reason
            assert passes == (s["genes"].get(g) == "released_e95" or g in excluded), (name, g, a, b)
            checked += 1
        assert all(g in [x for p in s["panel"] for x in panel[p]] for g, src in s["genes"].items()
                   if src == "released_e95"), name
    assert checked > 50
    assert set(COMMITTED["sets"]["vasc_split"]["excluded"]) == {"Lyve1"}


# ------------------------------------------------------------------ the score


def test_log1p_cp10k_uses_each_cells_own_library():
    C = sparse.csr_matrix(np.array([[1, 3, 0], [10, 30, 0]]))
    L = fr.log1p_cp10k(C).toarray()
    assert np.allclose(L[0], L[1]) and np.allclose(L[0], np.log1p([2500.0, 7500.0, 0.0]))
    with pytest.raises(fr.Refused):
        fr.log1p_cp10k(sparse.csr_matrix(np.array([[0, 0, 0]])))


def test_set_score_is_the_set_mean_minus_the_control_mean_cell_by_cell():
    L = sparse.csr_matrix(np.array([[2.0, 4.0, 1.0, 3.0], [0.0, 0.0, 1.0, 1.0]]))
    assert np.allclose(fr.set_score(L, [0, 1], [2, 3]), [3.0 - 2.0, 0.0 - 1.0])
    assert np.allclose(fr.set_score(L, [0, 1], []), [3.0, 0.0])


def test_a_cells_score_does_not_depend_on_the_other_cells_once_the_controls_are_fixed():
    cfg = syn_cfg()
    a = counts(cells(SPEC), seed=1)
    b = counts(cells({"cm": 200, "epi": 3}), seed=2)  # a very different population mix
    La, Lab = fr.log1p_cp10k(a), fr.log1p_cp10k(sparse.vstack([a, b]).tocsr())
    idx = [GENES.index(g) for g in SETS["epi"]["genes"]]
    ctrl = [GENES.index(g) for g in FILLER[:20]]
    assert np.allclose(fr.set_score(La, idx, ctrl), fr.set_score(Lab, idx, ctrl)[:a.shape[0]])
    assert cfg["score"]["control_seed"] == 20261119


def test_control_sets_are_seeded_expression_matched_and_exclude_every_set_gene():
    rng = np.random.default_rng(0)
    genes = [f"g{i}" for i in range(500)]
    means = np.r_[rng.gamma(1.0, 1.0, 495), np.zeros(5)]
    sets = {"a": [0, 1, 2], "b": [3, 4]}
    excl = {0, 1, 2, 3, 4, 10, 11}
    c1 = fr.control_sets(means, genes, sets, excl, 25, 10, 7)
    c2 = fr.control_sets(means, genes, sets, excl, 25, 10, 7)
    c3 = fr.control_sets(means, genes, sets, excl, 25, 10, 8)
    assert all(np.array_equal(c1[k], c2[k]) for k in sets)
    assert any(not np.array_equal(c1[k], c3[k]) for k in sets)
    assert len(c1["a"]) == 30 and len(c1["b"]) == 20
    for k in sets:
        assert not set(c1[k].tolist()) & excl
        assert np.all(means[c1[k]] > 0)
    pool = np.flatnonzero(means > 0)
    order = pool[np.lexsort((np.array(genes)[pool], means[pool]))]
    rank = {int(g): r for r, g in enumerate(order)}
    width = len(order) / 25
    for g in sets["a"]:  # controls come from the set gene's own bin
        assert any(abs(rank[int(c)] - rank[g]) < 2 * width for c in c1["a"])
    zero = fr.control_sets(means, genes, {"z": [497]}, set(), 25, 10, 7)
    assert len(zero["z"]) == 0  # a set gene with a zero mean has no bin and draws nothing


def test_gene_means_weigh_the_stages_equally():
    L = sparse.csr_matrix(np.array([[2.0], [2.0], [2.0], [0.0]]))
    stage = np.array(["E9.5", "E9.5", "E9.5", "E14.5"])
    assert fr.gene_means(L, stage, ["E9.5", "E14.5"])[0] == pytest.approx(1.0)
    with pytest.raises(fr.Refused):
        fr.gene_means(L, stage, ["E9.5", "E18"])


def test_duplicated_symbols_are_not_used():
    assert fr.unique_index(["A", "B", "A", "C"]) == {"B": 1, "C": 3}


# ------------------------------------------------------------------ thresholds, margins and the split


def _scores(**cols):
    base = {n: np.zeros(len(next(iter(cols.values())))) for n in SETS}
    base.update({k: np.asarray(v, dtype=float) for k, v in cols.items()})
    return base


def test_threshold_and_margin():
    cfg = syn_cfg()
    s = _scores(cm=[0.49, 0.5, 0.9, 0.9, 0.9], epi=[0.0, 0.0, 0.66, 0.65, 0.0], blood=[0, 0, 0, 0, 0.8])
    lab = fr.assign(s, cfg)
    # below min_score; at min_score with margin; margin 0.24 < 0.25; margin 0.25; margin 0.1
    assert lab.tolist() == ["r_unassigned", "r_cm", "r_unassigned", "r_cm", "r_unassigned"]
    assert fr.assign(s, cfg, min_score=0.3, min_margin=0.05)[0] == "r_cm"


def test_competing_class_takes_its_own_label_and_enters_no_group():
    cfg = syn_cfg()
    lab = fr.assign(_scores(blood=[2.0], cm=[0.5]), cfg)
    assert lab.tolist() == ["r_blood"]
    assert not any("r_blood" in labs for labs in cfg["groups"].values())
    bad = copy.deepcopy(cfg)
    bad["groups"]["CM"] = ["r_cm", "r_blood"]
    with pytest.raises(fr.Refused, match="competing"):
        fr.check_config(bad)


def test_the_ec_split_and_the_ambiguous_class():
    cfg = syn_cfg()
    s = _scores(ec=[1, 1, 1, 1, 1], endo_split=[1.0, 0.0, 0.5, 0.55, 0.0], vasc_split=[0.0, 1.0, 0.5, 0.5, 0.1])
    lab = fr.assign(s, cfg)
    assert lab.tolist() == ["r_endocardial_ec", "r_vascular_ec", "r_ec_ambiguous", "r_ec_ambiguous",
                            "r_vascular_ec"]
    groups = cfg["groups"]
    assert "r_ec_ambiguous" in groups["Endothelium"] and "r_ec_ambiguous" not in groups["Endocardium"]
    # the split never acts on a cell that is not called ec, whatever its split scores
    assert fr.assign(_scores(cm=[2.0], endo_split=[3.0]), cfg).tolist() == ["r_cm"]


def test_relabel_end_to_end_on_synthetic_counts():
    cfg = syn_cfg()
    types = cells(SPEC) + cells(SPEC)
    stage = np.array(["E9.5"] * sum(SPEC.values()) + ["E14.5"] * sum(SPEC.values()))
    res = fr.relabel(counts(types, seed=3), GENES, stage, cfg)
    want = {"cm": "r_cm", "endo": "r_endocardial_ec", "vasc": "r_vascular_ec", "epi": "r_epicardial",
            "blood": "r_blood", "none": "r_unassigned"}
    lab = res["labels"]
    for t, w in want.items():
        m = np.array(types) == t
        assert np.mean(lab[m] == w) >= 0.9, (t, w, np.unique(lab[m], return_counts=True))
    # A known limit, not a rule: a cm + epi doublet is not reliably dropped, since the cm set's controls come from a
    # higher expression bin than the epi set's (the more common type's markers have the higher pooled mean)
    assert set(lab[np.array(types) == "doublet"]) <= {"r_cm", "r_epicardial", "r_unassigned"}
    amb = lab[np.array(types) == "ambig"]  # both split sets high: always an endothelial class, often the ambiguous one
    assert set(amb) <= {"r_endocardial_ec", "r_vascular_ec", "r_ec_ambiguous"} and "r_ec_ambiguous" in set(amb)
    assert res["group_cells"]["Endothelium"]["E9.5"] >= res["group_cells"]["Endocardium"]["E9.5"]
    assert res["cap"]["passed"] and all(v < 0.40 for v in res["unassigned_share"].values())
    assert all(f["own"] for f in res["fallbacks"].values())
    assert set(res["sets"]["cm"]["found"]) == {"Cm1", "Cm2", "Cm3"}
    assert not set(res["sets"]["cm"]["controls"]) & {g for s in SETS.values() for g in s["genes"]}
    with pytest.raises(fr.Refused):  # a cell outside the yaml's stages
        fr.relabel(counts(["cm"]), GENES, np.array(["E18"]), cfg)


def test_a_set_with_too_few_genes_in_fengs_symbols_stops_the_run():
    cfg = syn_cfg()
    cfg["sets"]["epi"]["genes"] = {"Ep1": "syn", "Absent1": "syn", "Absent2": "syn"}
    with pytest.raises(fr.Refused, match="set epi"):
        fr.relabel(counts(cells(SPEC)), GENES, np.array(["E9.5"] * sum(SPEC.values())), cfg)


# ------------------------------------------------------------------ the fallback and the cap


def test_group_level_fallback():
    cells_ = {"CM": {"E9.5": 50, "E14.5": 80}, "Endocardium": {"E9.5": 9, "E14.5": 40},
              "Endothelium": {"E9.5": 30, "E14.5": 60}, "Epicardium": {"E9.5": 4, "E14.5": 30}}
    fb = fr.fallbacks(cells_, 10, COMMITTED["fallback"])
    assert fb["CM"] == {"own": True, "fired": False, "takes": "CM"}
    assert fb["Endocardium"] == {"own": False, "fired": True, "takes": "Endothelium"}
    assert fb["Epicardium"] == {"own": False, "fired": True, "takes": None}
    # never chained: Endocardium's fallback must have its own direction
    fb2 = fr.fallbacks(dict(cells_, Endothelium={"E9.5": 2, "E14.5": 60}), 10, COMMITTED["fallback"])
    assert fb2["Endocardium"]["takes"] is None and fb2["Endothelium"]["takes"] is None
    # feng_dirs applies the same rule to directions
    d = {"CM": np.ones(3), "Endothelium": np.full(3, 2.0)}
    out, rec = feng_dirs.apply_fallback(d, COMMITTED)
    assert np.array_equal(out["Endocardium"], d["Endothelium"]) and out["Endocardium"] is not d["Endothelium"]
    assert "Epicardium" not in out and rec["Epicardium"] == {"own": False, "fired": True, "takes": None}
    out2, _ = feng_dirs.apply_fallback({"CM": np.ones(3)}, COMMITTED)
    assert set(out2) == {"CM"}


def test_the_per_stage_unassigned_cap():
    lab = np.array(["r_unassigned"] * 4 + ["r_cm"] * 6 + ["r_unassigned"] * 5 + ["r_cm"] * 5)
    stage = np.array(["E9.5"] * 10 + ["E14.5"] * 10)
    sh = fr.unassigned_shares(lab, stage, ["E9.5", "E14.5"], "r_unassigned")
    assert sh == {"E9.5": 0.4, "E14.5": 0.5}
    c = fr.cap_check(sh, 0.40)
    assert c["stages_over"] == ["E14.5"] and not c["passed"]
    assert fr.cap_check({"E9.5": 0.4, "E14.5": 0.39}, 0.40)["passed"]
    assert not fr.cap_check({"E9.5": None, "E14.5": 0.1}, 0.40)["passed"]
    # a competing class is assigned, not unassigned
    assert fr.unassigned_shares(np.array(["r_blood"]), np.array(["E9.5"]), ["E9.5"], "r_unassigned")["E9.5"] == 0.0


# ------------------------------------------------------------------ the whitelist


HEADER = ["", "orig.ident", "strain", "stage", "zone", "cell_type_spec", "seurat_clusters", "UMAP_1"]


def _obs(path: Path, n: int, author: list[str]) -> None:
    with open(path, "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(HEADER)
        for i in range(n):
            w.writerow([f"c{i}", "s1", "CD1", "E9.5" if i % 2 == 0 else "E14.5", ["LA", "RA", "LV", "RV"][i % 4],
                        author[i], str(i % 7), f"{i * 0.1:.2f}"])


def test_read_obs_reads_the_whitelist_and_nothing_else(tmp_path):
    p = tmp_path / "obs.csv"
    _obs(p, 6, ["ventricular_cm"] * 6)
    got = fr.read_obs(p)
    assert set(got) == {"cell", "stage", "zone", "strain"}
    assert got["cell"].tolist() == [f"c{i}" for i in range(6)] and got["strain"].tolist() == ["CD1"] * 6
    for col in ("cell_type_spec", "seurat_clusters", "UMAP_1", "orig.ident"):
        with pytest.raises(fr.Refused):
            fr.read_obs(p, ["stage", col])
    for col in ("cell_type", "celltype_2", "clusters", "annotation", "X_pca", "harmony_1", "leiden", "label"):
        with pytest.raises(fr.Refused):
            fr.check_whitelist([col])
    bad = copy.deepcopy(COMMITTED)
    bad["obs_columns"] = ["stage", "cell_type_spec"]
    with pytest.raises(fr.Refused):
        fr.check_config(bad)


def test_no_selector_or_listing_module_names_the_authors_column():
    """Only the method-report script may name cell_type_spec. The labeller, feng_dirs, the builder, the harness and
    N54's listing never do, and none of them imports the agreement script."""
    for rel in ("experiments/external/feng_relabel.py", "experiments/external/feng_relabel_dryrun.py",
                "experiments/t1_m1/feng_dirs.py",
                "experiments/builds/build_test.py", "src/vec/harness.py", "experiments/t1_m1/n54_relabel_listing.py"):
        tree = ast.parse((ROOT / rel).read_text())
        docs = {id(n.body[0].value) for n in ast.walk(tree)
                if isinstance(n, (ast.Module, ast.FunctionDef, ast.ClassDef)) and ast.get_docstring(n)}
        code_strings = [n.value for n in ast.walk(tree)
                        if isinstance(n, ast.Constant) and isinstance(n.value, str) and id(n) not in docs]
        assert not [s for s in code_strings if "cell_type_spec" in s or "type_spec" == s], rel  # comments aside
        imported = {a.name for n in ast.walk(tree) if isinstance(n, ast.Import) for a in n.names}
        imported |= {n.module for n in ast.walk(tree) if isinstance(n, ast.ImportFrom) and n.module}
        assert not {m for m in imported if m.startswith("n54_authors_agreement")}, rel
    agreement = (ROOT / "experiments/t1_m1/n54_authors_agreement.py").read_text()
    assert "FOR THE METHOD REPORT ONLY" in agreement and 'AUTHORS_COLUMN = "cell_type_spec"' in agreement


def _src(tmp_path: Path, author: list[str], types: list[str]) -> Path:
    src = tmp_path / "src"
    src.mkdir(parents=True, exist_ok=True)
    n = len(types)
    _obs(src / "cd1_legal_obs.csv", n, author)
    C = counts(types, seed=5)
    sparse.save_npz(src / "cd1_legal_counts.npz", C)
    (src / "cd1_legal_genes.txt").write_text("\n".join(GENES) + "\n")
    fio.write_manifest(src, {"outputs": {f: {"sha256": fio.sha256_file(src / f)} for f in
                                         ("cd1_legal_obs.csv", "cd1_legal_counts.npz", "cd1_legal_genes.txt")}})
    return src


def _cfg_file(tmp_path: Path, cfg: dict) -> Path:
    p = tmp_path / "relabel.yaml"
    p.write_text(yaml.safe_dump(cfg, sort_keys=False))
    return p


def test_run_writes_outputs_that_do_not_depend_on_the_authors_column(tmp_path):
    types = [t for t in cells(SPEC) for _ in range(2)]
    cfg = syn_cfg(status="frozen")
    cfg_path = _cfg_file(tmp_path, cfg)
    outs = []
    for k, author in enumerate((["ventricular_cm"] * len(types),
                                [f"x{(i * 7) % 5}" for i in range(len(types))])):
        src = _src(tmp_path / f"r{k}", author, types)
        dst = tmp_path / f"r{k}" / "relabel_c"
        man = fr.run(cfg, cfg_path, src, dst, GATE)
        outs.append((dst, man))
    (d0, m0), (d1, m1) = outs
    for name in (cfg["output"]["stats"], cfg["output"]["cells"]):
        assert (d0 / name).read_bytes() == (d1 / name).read_bytes(), name
    assert m0["authors_labels_read"] is False and m0["obs_columns_read"] == ["<cell names>", "stage", "zone", "strain"]
    assert m0["config"]["sha256"] == fio.sha256_file(cfg_path) and m0["usable"] is True
    assert m0["gate"] == GATE and m0["route"] == "recomputed" and "allow_unsourced" not in m0
    assert set(m0["unassigned_kinds"]["E9.5"]) == {"no_class_reached_min_score", "margin_failed"}
    assert m0["outputs"][cfg["output"]["stats"]]["usable"] is True
    z = np.load(d0 / cfg["output"]["stats"])
    keys = [str(k) for k in z["keys"]]
    assert all(k.split("|")[2].startswith("r_") and len(k.split("|")) == 4 for k in keys)
    assert {k.split("|")[0] for k in keys} == {"E9.5", "E14.5"}
    assert fio.verify_manifest_output(d0 / cfg["output"]["stats"])["usable"] is True
    # a strain other than the expected one stops the run
    bad = copy.deepcopy(cfg)
    bad["expect_strain"] = "B6"
    with pytest.raises(fr.Refused, match="strain"):
        fr.run(bad, cfg_path, d0.parent / "src", tmp_path / "bad", GATE)


# ------------------------------------------------------------------ the selector


GROUPS4 = ("CM", "Endocardium", "Endothelium", "Epicardium")


def _relabel_stats(dirpath: Path, cfg_path: Path, labels: dict[str, int], usable: bool = True,
                   source: str = "feng_relabel", genes: list[str] | None = None, **man_over) -> Path:
    """A synthetic recomputed statistics file with the manifest feng_relabel writes."""
    genes = genes or [f"f{i}" for i in range(30)]
    rng = np.random.default_rng(4)
    keys, n, cp = [], [], []
    for lab, k in labels.items():
        for stage in ("E9.5", "E14.5"):
            prof = rng.gamma(0.6, 20.0, len(genes))
            for half in "AB":
                keys.append(f"{stage}|LV|{lab}|{half}")
                n.append(k)
                cp.append(k * prof)
    dirpath.mkdir(parents=True, exist_ok=True)
    path = dirpath / "cd1_stats_relabel.npz"
    n = np.array(n, dtype=float)
    np.savez(path, genes=np.array(genes), keys=np.array(keys), n=n, S_cp=np.array(cp),
             S_log=np.zeros((len(keys), len(genes))), S_nz=np.zeros((len(keys), len(genes))), depth_sum=n * 1e3)
    man = {"source": source, "config": {"sha256": fio.sha256_file(cfg_path), "status": "frozen"},
           "cap": {"passed": usable}, "unsourced": {}, "gate": dict(GATE),
           "provenance": {"git_commit": "abc", "git_uncommitted": []},
           "outputs": {path.name: {"sha256": fio.sha256_file(path), "usable": usable}}}
    for k, v in man_over.items():  # "config.status": "draft" style overrides
        head, _, tail = k.partition(".")
        if tail:
            man[head] = {**man[head], tail: v}
        else:
            man[head] = v
    fio.write_manifest(dirpath, man)
    return path


def test_recomputed_directions_use_the_yaml_groups_and_its_fallback(tmp_path, record_ok):
    genes = [f"f{i}" for i in range(30)]
    path = _relabel_stats(tmp_path / "rc", fr.CONFIG, {"r_cm": 30, "r_endocardial_ec": 4, "r_vascular_ec": 30,
                                                      "r_ec_ambiguous": 5, "r_epicardial": 12, "r_blood": 40})
    dirs, rec = feng_dirs.feng_directions(genes, stats_path=path, groups="recomputed")
    assert set(dirs) == set(GROUPS4) and "Blood" not in dirs
    assert rec["feng_groups"] == "recomputed" and rec["relabel_config_sha256"] == fio.sha256_file(fr.CONFIG)
    assert rec["fallbacks"]["Endocardium"] == {"own": False, "fired": True, "takes": "Endothelium"}
    assert np.array_equal(dirs["Endocardium"], dirs["Endothelium"])
    assert rec["groups"]["Endothelium"]["E9.5"] == 2 * (30 + 4 + 5)  # every endothelial class, both halves
    assert rec["sha256"] == fio.sha256_file(path)
    assert rec["route"] == "recomputed" and rec["relabel_gate"] == GATE
    assert rec["with_direction"] == sorted(GROUPS4)
    # Epicardium under the minimum: it has no direction (absent), and the route stays recomputed
    path2 = _relabel_stats(tmp_path / "rc2", fr.CONFIG, {"r_cm": 30, "r_endocardial_ec": 30, "r_vascular_ec": 30,
                                                        "r_epicardial": 3})
    dirs2, rec2 = feng_dirs.feng_directions(genes, stats_path=path2, groups="recomputed")
    assert "Epicardium" not in dirs2 and rec2["with_direction"] == ["CM", "Endocardium", "Endothelium"]
    assert rec2["route"] == "recomputed" and rec2["fallbacks"]["Epicardium"]["takes"] is None


def test_recomputed_directions_refuse_a_manifest_not_committed_as_the_record(tmp_path, monkeypatch):
    path = _relabel_stats(tmp_path / "rc", fr.CONFIG, {"r_cm": 30})
    man = path.parent / "MANIFEST.json"
    for head in (None, b"{}\n"):
        monkeypatch.setattr(fr, "committed_at_head", lambda p, h=head: h)
        with pytest.raises(ValueError, match="not committed at HEAD as the labeller's RECORD"):
            feng_dirs.feng_directions(["f1"], stats_path=path, groups="recomputed")
    monkeypatch.setattr(fr, "committed_at_head", lambda p: man.read_bytes() if Path(p) == fr.RECORD else None)
    feng_dirs.feng_directions(["f1"], stats_path=path, groups="recomputed")  # committed byte for byte: read


@pytest.mark.parametrize("kw, msg", [({"usable": False}, "unusable"), ({"source": "feng2022_geo"}, "not written by"),
                                     ({"labels": {"r_cm": 30, "ventricular_cm": 30}}, "not recomputed"),
                                     ({"config.status": "draft"}, "not frozen"),
                                     ({"unsourced": {"immune": ["Ptprc"]}}, "unsourced"),
                                     ({"allow_unsourced": True}, "unsourced"),
                                     ({"provenance.git_uncommitted": [" M experiments/external/feng_relabel.py"]},
                                      "uncommitted"),
                                     ({"gate": None}, "without the labeller's gate"),
                                     ({"gate": {**GATE_, "clean": False}}, "without the labeller's gate"),
                                     ({"gate": {**GATE_, "row": "N12"}}, "without the labeller's gate")])
def test_recomputed_directions_refuse(tmp_path, kw, msg, record_ok):
    labels = kw.pop("labels", {"r_cm": 30})
    path = _relabel_stats(tmp_path / "rc", fr.CONFIG, labels, **kw)
    with pytest.raises(ValueError, match=msg):
        feng_dirs.feng_directions(["f1", "f2"], stats_path=path, groups="recomputed")


def test_recomputed_directions_refuse_a_file_made_from_another_yaml(tmp_path, record_ok):
    other = tmp_path / "other.yaml"
    other.write_text(fr.CONFIG.read_text() + "\n# edited\n")
    path = _relabel_stats(tmp_path / "rc", other, {"r_cm": 30})
    with pytest.raises(ValueError, match="another feng_relabel.yaml"):
        feng_dirs.feng_directions(["f1"], stats_path=path, groups="recomputed")
    with pytest.raises(ValueError, match="groups must be"):
        feng_dirs.feng_directions(["f1"], groups="jax")


def test_the_authors_default_is_unchanged(tmp_path, monkeypatch):
    """feng_directions with no groups argument and with groups='authors' are the same call on the same file, and the
    recomputed path is never touched by the default."""
    monkeypatch.setattr(feng_dirs, "REPO", tmp_path)
    genes = [f"f{i}" for i in range(30)]
    stats = feng_dirs.clean_stats_path()
    rng = np.random.default_rng(3)
    keys, n, cp = [], [], []
    for lab in ("ventricular_cm", "atrial_cm", "endocardial_ec", "vascular_ec", "epicardial"):
        for stage in ("E9.5", "E14.5"):
            for half in "AB":
                keys.append(f"{stage}|LV|{lab}|{half}")
                n.append(20)
                cp.append(20 * rng.gamma(0.6, 20.0, len(genes)))
    stats.parent.mkdir(parents=True, exist_ok=True)
    n = np.array(n, dtype=float)
    np.savez(stats, genes=np.array(genes), keys=np.array(keys), n=n, S_cp=np.array(cp),
             S_log=np.zeros((len(keys), 30)), S_nz=np.zeros((len(keys), 30)), depth_sum=n)
    fio.write_manifest(stats.parent, {"outputs": {stats.name: {"sha256": fio.sha256_file(stats)}}})
    monkeypatch.setattr(feng_dirs, "recomputed_directions", lambda *a, **k: pytest.fail("not the default"))
    d0, r0 = feng_dirs.feng_directions(genes)
    d1, r1 = feng_dirs.feng_directions(genes, groups="authors")
    assert r0 == r1 and "feng_groups" not in r0 and set(d0) == set(GROUPS4)
    assert all(np.array_equal(d0[g], d1[g], equal_nan=True) for g in d0)


def test_groups_config_replaces_only_the_feng_lists():
    cfg = feng_dirs.load_config()
    g = feng_dirs.groups_config(cfg, COMMITTED)
    assert g["external_labels"]["Endothelium"][0] == ["r_endocardial_ec", "r_vascular_ec", "r_ec_ambiguous"]
    assert g["external_labels"]["Blood"][0] is None and g["external_labels"]["ventricularCM"][0] is None
    assert g["external_labels"]["CM"][1:] == cfg["external_labels"]["CM"][1:]
    assert cfg["external_labels"]["CM"][0] == ["ventricular_cm", "atrial_cm"]  # the readout's config is not edited
    assert g["harmonise"] == cfg["harmonise"] and g["metrics"] == cfg["metrics"]


@pytest.fixture
def bt():
    import build_test

    return build_test


def test_builder_default_is_authors_and_makes_the_old_call(bt, monkeypatch):
    ap = bt.build_parser()
    assert ap.parse_args(["t1", "--out", "x.h5ad"]).feng_groups == "authors"
    assert ap.parse_args(["t1", "--feng-groups", "recomputed", "--out", "x.h5ad"]).feng_groups == "recomputed"
    with pytest.raises(SystemExit):
        ap.parse_args(["t1", "--feng-groups", "jax", "--out", "x.h5ad"])
    calls = []
    monkeypatch.setattr(bt.feng_dirs, "feng_directions", lambda genes, *a, **k: calls.append((genes, a, k)) or
                        ({"CM": np.zeros(1)}, {"file": "f", "sha256": "s", "route": "recomputed"}))
    bt._feng_directions(["g"])
    bt._feng_group_names()
    bt._feng_directions(["g"], "recomputed")
    assert calls == [(["g"], (), {}), ([""], (), {}), (["g"], (), {"groups": "recomputed"})]


def test_builder_refuses_an_unusable_recomputed_file_before_anything_else(bt, monkeypatch):
    def refuse(genes, *a, **k):
        raise ValueError("is marked unusable")

    monkeypatch.setattr(bt.feng_dirs, "feng_directions", refuse)
    with pytest.raises(SystemExit, match="--feng-groups recomputed is refused"):
        bt._feng_directions(["g"], "recomputed")


RC_FILE = "data/derived/feng_relabel_c/cd1_stats_relabel.npz"


def _recomputed_env(bt, env, monkeypatch, groups=GROUPS4, route="recomputed", record_route="recomputed"):
    """The build_t1 fixture made ready for --feng-groups recomputed: library suffixes, the E9.5 load g_mom's fit
    needs, R2's union (synthetic), the labeller's record at HEAD with `record_route`, and an answer for the recomputed
    call (the fixture's synthetic directions on `groups` only, with the recomputed record and `route`). Returns the
    log of g_mom's loads."""
    from test_g_mom_verify import _g_mom_loads
    from test_t1_feng_groups_none import _lib3

    _lib3(env)
    synthetic = bt.feng_dirs.feng_directions

    def directions(genes, *a, **k):
        if k.get("groups") != "recomputed":
            return synthetic(genes, *a, **k)
        env.reads.append(("recomputed", list(genes)))
        d, _ = synthetic(genes)
        d = {g: v for g, v in d.items() if g in groups}
        return d, {"file": RC_FILE, "sha256": "ab" * 32, "relabel_config_sha256": "cd" * 32, "route": route,
                   "with_direction": sorted(d), "relabel_gate": GATE,
                   "fallbacks": {g: {"own": g in d, "fired": g not in d, "takes": g if g in d else None}
                                 for g in GROUPS4}}

    monkeypatch.setattr(bt.feng_dirs, "feng_directions", directions)
    head = None if record_route is None else json.dumps({"route": record_route}).encode()
    monkeypatch.setattr(fr, "committed_at_head", lambda p: head if Path(p) == fr.RECORD else None)
    monkeypatch.setattr(bt.hygiene, "load_union", lambda *a, **k: UNION)
    return _g_mom_loads(bt, env, monkeypatch)


def test_builder_records_the_recomputed_file_only_when_asked(bt, t1_env, monkeypatch, tmp_path):  # noqa: F811
    """End to end through build_t1 on synthetic stages (the fixture of tests/test_hidden_counts.py, every input
    monkeypatched): the default file is byte-identical to the reference built directly and its record has no
    feng_groups field; with recomputed the builder asks feng_dirs for the recomputed groups (also for R2's Blood guard,
    never for the authors' groups) and the record names the file, its sha256, the yaml's sha256 and the route."""
    today = _reference_file(bt, t1_env, tmp_path / "today.h5ad", True)
    out = tmp_path / "T1_s1.h5ad"
    bt.main(["t1", "--no-feng-hygiene", "--out", str(out)])
    assert out.read_bytes() == today.read_bytes()
    assert "feng_groups" not in json.loads(out.with_suffix(".build.json").read_text())
    t1_env.reads.clear()
    loads = _recomputed_env(bt, t1_env, monkeypatch)
    out2 = tmp_path / "T1_rc.h5ad"
    bt.main(["t1", *TARGET, "--feng-groups", "recomputed", "--out", str(out2)])
    assert [r[0] for r in t1_env.reads if isinstance(r, tuple)] == ["recomputed", "Feng", "recomputed", "Feng"]
    assert loads == ["E10.5_RNA.h5ad", "E9.5_RNA.h5ad"]
    rec = json.loads(out2.with_suffix(".build.json").read_text())
    fg = rec["feng_groups"]
    assert fg["groups"] == "recomputed" and fg["sha256"] == "ab" * 32 and fg["file"] == RC_FILE
    assert fg["relabel_config_sha256"] == "cd" * 32 and fg["route"] == "recomputed" and fg["relabel_gate"] == GATE
    assert fg["with_direction"] == sorted(GROUPS4) and rec["external"]["sha256"] == "ab" * 32
    assert rec["base"] == "lib3" and rec["feng_hygiene"]["sha256"] == UNION["sha256"]  # R2 on, as in the first upload
    assert "recomputed groups" in rec["g_mom"]["rule"] and "R1_FENG_GROUPS is not read" in rec["g_mom"]["mask"]
    assert rec["fraction"] == bt.feng_dirs.bracket_fraction(10.5, 12.5)


def test_harness_late_proxy_option_is_off_by_default():
    reg = [c for name, c in harness.CANDIDATES.items() if name.startswith("lp_")]
    assert reg and all(c.derived == ("t1_late_proxy",) and "feng_groups" not in c.params for c in reg)
    c = harness._lp("lp_rc", "synthetic", feng_groups="recomputed")
    assert c.derived == ("t1_late_proxy_recomputed",) and c.params["feng_groups"] == "recomputed"
    assert "lp_rc" not in harness.CANDIDATES
    assert harness.LP_STATS_RECOMPUTED == feng_dirs.recomputed_stats_path(COMMITTED)
    assert harness.DERIVED["t1_late_proxy"].path == harness.LP_STATS
    with pytest.raises(ValueError):
        harness._late_proxy(feng_groups="jax")


# ------------------------------------------------------------------ row N54's listing


def test_listing_compare_counts_sign_changes_and_large_ones():
    old = np.array([1.0, -1.0, 0.5, 0.01, np.nan, 2.0])
    new = np.array([1.0, 1.0, -0.5, -0.01, 3.0, np.nan])
    genes = ["a", "b", "c", "d", "e", "held"]
    moved = n54.moved_mask(old, genes, ["held"])
    assert moved.tolist() == [True, True, True, True, False, False]
    c = n54.compare(new, old, moved)
    assert c["genes_moved"] == 4 and c["sign_changes"] == 3
    # 0.4 x 0.01 is below log(1.1) in both versions; 0.4 x 1 and 0.4 x 0.5 are above
    assert c["sign_changes_large"] == 2
    assert c["cosine"] == pytest.approx(float(old[:4] @ new[:4] / np.linalg.norm(old[:4]) / np.linalg.norm(new[:4])))
    assert c["pearson"] == pytest.approx(float(np.corrcoef(old[:4], new[:4])[0, 1]))
    same = n54.compare(old, old, moved)
    assert same["sign_changes"] == 0 and same["cosine"] == pytest.approx(1.0)
    gap = n54.compare(np.array([np.nan] * 6), old, moved)
    assert gap["genes_without_recomputed_value"] == 4 and gap["cosine"] is None


def test_listing_group_rows_and_sensitivity():
    man = {"group_cells": {g: {"E9.5": 12, "E14.5": 40} for g in GROUPS4}}
    rec = {"fallbacks": {"Endocardium": {"own": False, "fired": True, "takes": "Endothelium"}}}
    rows = n54.group_rows(man, rec)
    assert rows["Endocardium"]["fallback_fired"] and rows["Endocardium"]["takes"] == "Endothelium"
    assert not rows["CM"]["fallback_fired"]
    cfg = syn_cfg()
    s = _scores(cm=[0.6, 0.4], ec=[0.0, 0.0])
    rows2 = n54.sensitivity(s, np.array(["E9.5", "E14.5"]), cfg)
    assert [r["setting"] for r in rows2] == cfg["sensitivity"]
    assert rows2[0]["group_cells"]["CM"] == {"E9.5": 1, "E14.5": 1}  # min_score 0.3
    assert rows2[1]["group_cells"]["CM"] == {"E9.5": 0, "E14.5": 0}  # min_score 0.7


def test_the_ruling_is_found_anywhere_in_its_decision_log_entry(tmp_path):
    """The real decision log: '- <id>, headline' at the margin, the body indented below it."""
    plan = tmp_path / "PLAN.md"
    plan.write_text("## Decision log\n\n"
                    "- 2026-11-20 (a), the captain's ruling on T1's second file (state at 09:00 UTC). Written by\n"
                    "  the AI assistant.\n"
                    "  - **The ruling.** Recomputed groups by feng_relabel.yaml; row N54 runs after it.\n"
                    "- 2026-11-20 (b), something else; N540 only\n"
                    "- 2026-11-20 (aa), not (a)\n"
                    "## Next\n  N54 outside any entry\n")
    assert n54.ruling_logged("2026-11-20 (a)", "N54", plan) and fr.ruling_logged("2026-11-20 (a)", "N54", plan)
    assert not n54.ruling_logged("2026-11-20 (b)", "N54", plan)  # N540 is another id
    assert not n54.ruling_logged("2026-11-21 (a)", "N54", plan)
    assert not n54.ruling_logged("2026-11-20 (aa)", "N54", plan)
    assert not n54.ruling_logged("2026-11-20", "N54", plan) and not n54.ruling_logged(None, "N54", plan)


def test_listing_reads_no_authors_label_and_decides_nothing():
    src = (ROOT / "experiments/t1_m1/n54_relabel_listing.py").read_text()
    assert not re.search(r"import\s+n54_authors_agreement|from\s+n54_authors_agreement", src)
    for rel in ("experiments/builds/build_test.py", "src/vec/harness.py", "experiments/t1_m1/feng_dirs.py",
                "experiments/external/feng_relabel.py"):  # nothing that builds or labels reads N54's outputs
        assert not re.search(r"import\s+n54|from\s+n54|n54_relabel_listing|results/n54", (ROOT / rel).read_text())


def test_agreement_is_aggregate_only():
    import n54_authors_agreement as ag

    a = {"CM": np.array([1, 1, 0, 0], bool)}
    r = {"CM": np.array([1, 0, 1, 0], bool)}
    out = ag.agreement(a, r)["CM"]
    assert out == {"authors": 2, "recomputed": 2, "both": 1, "recovered_share": 0.5, "elsewhere_share": 0.5,
                   "jaccard": pytest.approx(1 / 3)}
    json.dumps(out)


# ------------------------------------------------------------------ fixes after review (2026-11-19)


def test_controls_never_draw_a_panel_gene_or_a_gene_left_out_by_name():
    cfg = syn_cfg()
    cfg["sets"]["vasc_split"]["excluded"] = {"Kept1": "left out by name"}
    panel = sorted(fr.panel_genes())
    assert {"Myl7", "Lyve1", "Hapln1", "Tbx18"} <= set(panel)
    genes = GENES + panel + ["Kept1"]
    types = cells(SPEC)
    extra = np.random.default_rng(1).poisson(2.0, (len(types), len(panel) + 1))
    C = sparse.hstack([counts(types, seed=8), sparse.csr_matrix(extra)]).tocsr()
    res = fr.relabel(C, genes, np.array(["E9.5"] * len(types)), {**cfg, "stages": ["E9.5"]})
    drawn = {g for s in res["sets"].values() for g in s["controls"]}
    assert drawn and not drawn & (set(panel) | {"Kept1"})


def test_an_unsplit_ec_cell_is_refused():
    cfg = syn_cfg()
    s = _scores(ec=[1.0], endo_split=[np.nan], vasc_split=[0.0])  # NaN falls in no split band
    with pytest.raises(fr.Refused, match="unsplit"):
        fr.assign(s, cfg)


def test_check_config_keeps_unsourced_genes_in_competing_classes_and_needs_the_terminal_rule():
    bad = syn_cfg()
    bad["sets"]["cm"]["genes"]["CmX"] = fr.UNSOURCED
    with pytest.raises(fr.Refused, match="only a competing class"):
        fr.check_config(bad)
    for key, val, msg in (("otherwise", "recomputed", "otherwise"), ("requires_direction", ["Blood"], "requires"),
                          ("none_with_g_mom", "maybe", "none_with_g_mom"),
                          ("none_with_g_mom", False, "none_with_g_mom"),
                          ("base", "pooled", "base"), ("base", "mid", "base"), ("on_stop", "recomputed", "on_stop"),
                          ("on_stop", None, "on_stop")):
        bad = syn_cfg()
        bad["terminal"] = {**bad["terminal"], key: val}
        with pytest.raises(fr.Refused, match=msg):
            fr.check_config(bad)


def test_the_unsourced_rule_deletes_genes_then_emptied_classes():
    assert fr.without_unsourced(COMMITTED) == COMMITTED  # nothing left to delete since 2026-11-19
    cfg = copy.deepcopy(COMMITTED)  # the rule itself, on a copy with unsourced genes put back
    cfg["sets"]["blood"]["genes"]["Hba-a1"] = fr.UNSOURCED
    cfg["sets"]["lymphatic"]["genes"]["Flt4"] = fr.UNSOURCED
    fr.check_config(cfg)
    out = fr.without_unsourced(cfg)
    assert not fr.unsourced_genes(out)
    assert set(cfg["sets"]) - set(out["sets"]) == {"lymphatic"}  # one gene left: the class goes whole
    assert set(out["sets"]["blood"]["genes"]) == {"Hba-x", "Hbb-y", "Gata1", "Hbb-bs", "Alas2"}
    assert out["groups"] == COMMITTED["groups"] and out["unassigned_cap"] == COMMITTED["unassigned_cap"]


# --- the literature sources of 2026-11-19 and the three dropped genes

LIT = {  # gene -> its source ids, as the main session's decisions of 2026-11-19 set them
    "Hba-a1": ["lit_E"], "Hbb-bs": ["lit_E"], "Alas2": ["lit_E"], "Plp1": ["lit_A"],
    "Col1a2": ["lit_B"], "Col3a1": ["lit_D"], "Dcn": ["lit_C"], "Myh11": ["lit_D"], "Pdgfrb": ["lit_A"],
    "Cspg4": ["lit_A"], "Rgs5": ["lit_C"], "Kcnj8": ["lit_C"], "Ptprc": ["lit_D"], "C1qb": ["lit_G"],
    "Lyz2": ["lit_H"], "Cd68": ["lit_A"], "Prox1": ["lit_J"], "Flt4": ["lit_J"]}


def test_every_formerly_unsourced_gene_has_its_literature_source_and_sentence():
    cfg = COMMITTED
    got = {g: fr.gene_sources(src) for s in cfg["sets"].values() for g, src in s["genes"].items()
           if any(sid.startswith("lit_") for sid in fr.gene_sources(src))}
    assert got == LIT
    for g, sids in LIT.items():
        for sid in sids:
            ref = cfg["sources"][sid]
            assert all(ref[k] for k in fr.LITERATURE_KEYS) and ref["doi"].startswith("10.")
            assert g in ref["names"][g], (g, sid)  # the quoted sentence names the gene
    assert set(cfg["sources"]) == {"released_e95", "UNSOURCED", "lit_A", "lit_B", "lit_C", "lit_D", "lit_E", "lit_G",
                                   "lit_H", "lit_J"}
    # Pdgfrb keeps lit_B's caveat; Prox1 and Flt4 keep lit_J's
    assert "fibroblast" in cfg["sources"]["lit_B"]["caveats"]["Pdgfrb"]
    assert set(cfg["sources"]["lit_J"]["caveats"]) >= {"Prox1", "Flt4", "Pdpn"}


def test_the_quotation_checks_corrections_are_applied():
    """The check of 2026-11-19 against the eight papers: lit_C names the globins only in a quality-control filter, so it
    is no source for them; lit_E's class is erythroid progenitors; the caveats are the papers' own sentences, with the
    team's reading marked; the ages and sexes as published."""
    src = COMMITTED["sources"]
    assert not {"Hba-a1", "Hbb-bs"} & set(src["lit_C"]["names"])
    assert "contamination by immune cells and erythrocytes" in src["lit_C"]["caveats"]["Hba-a1"]
    assert "progenitor" in src["lit_E"]["class_named"] and "both sexes" in src["lit_E"]["samples"]
    assert "commonly used markers Pdgfrb" in src["lit_B"]["caveats"]["Pdgfrb"]
    assert "team's reading" in src["lit_B"]["caveats"]["Pdgfrb"]
    assert "PROX1 was mainly expressed by valvular endothelial cells" in src["lit_J"]["caveats"]["Prox1"]
    assert "PDPN was expressed by the epicardium and fibroblasts" in src["lit_J"]["caveats"]["Pdpn"]
    assert "male and female" in src["lit_J"]["samples"] and "8-week-old male" in src["lit_A"]["samples"]
    for sid in ("lit_C", "lit_G"):
        assert "both sexes" in src[sid]["samples"], sid
    assert "BDNF" in src["lit_D"]["samples"] and "one week after" in src["lit_H"]["samples"]
    assert "3 weeks" in src["lit_G"]["samples"]
    # the check counted 20 gene statements; lit_C's two globin statements are now caveats, so 18 remain, one per gene
    assert sum(len(src[s]["names"]) for s in src if s.startswith("lit_")) == 18
    assert "PDPN" in COMMITTED["dropped"]["Pdpn"] and "protein" in COMMITTED["dropped"]["Pdpn"]


def test_the_three_dropped_genes_are_gone_dated_and_kept_out_of_the_controls():
    cfg = COMMITTED
    assert set(cfg["dropped"]) == {"Phox2b", "Fcer1g", "Pdpn"}
    assert all(cfg["dropped"][g].startswith("2026-11-19") for g in cfg["dropped"])
    members = {g for s in cfg["sets"].values() for g in s["genes"]}
    assert not members & set(cfg["dropped"])
    assert set(cfg["dropped"]) <= fr.named_genes(cfg)  # never drawn as a control
    bad = copy.deepcopy(cfg)
    bad["sets"]["immune"]["genes"]["Fcer1g"] = "released_e95"
    with pytest.raises(fr.Refused, match="both dropped"):
        fr.check_config(bad)
    bad["sets"]["immune"]["genes"]["Fcer1g"] = "lit_D"
    del bad["dropped"]["Fcer1g"]
    with pytest.raises(fr.Refused, match="quotes no sentence"):  # lit_D does not name it
        fr.check_config(bad)


def test_the_changed_classes_keep_the_minimum_of_genes():
    """With the three genes gone: neural_crest 6 (Phox2b gone), immune 4 (Fcer1g gone), lymphatic 2 (Pdpn gone, at
    the minimum exactly), mural 5 (unchanged); a class whose genes are all found always has at least the minimum."""
    cfg = COMMITTED
    n = {name: len(s["genes"]) for name, s in cfg["sets"].items()}
    need = int(cfg["score"]["min_genes_found"])
    assert need == 2 and all(v >= need for v in n.values())
    got = (n["mural"], n["immune"], n["lymphatic"], n["neural_crest"], n["blood"], n["mesenchyme"])
    assert got == (5, 4, 2, 6, 6, 7)


@pytest.mark.parametrize("bad_src, msg", [("lit_Z", "not under"), (["lit_C", "nope"], "not under"),
                                          (["lit_C", "lit_A"], "lit_A quotes no sentence"), ([], "an id or"),
                                          ({"lit_A": 1}, "an id or"), ("UNSOURCED_x", "not under")])
def test_check_config_refuses_a_bad_gene_source(bad_src, msg):
    bad = copy.deepcopy(COMMITTED)
    bad["sets"]["mural"]["genes"]["Rgs5"] = bad_src
    with pytest.raises(fr.Refused, match=msg):
        fr.check_config(bad)
    bad = copy.deepcopy(COMMITTED)
    del bad["sources"]["lit_C"]["doi"]
    with pytest.raises(fr.Refused, match="lacks"):
        fr.check_config(bad)


def test_run_refuses_a_draft_yaml_or_a_gate_without_the_row_or_a_clean_tree():
    cfg = syn_cfg(status="frozen")
    for gate, msg in (({**GATE, "row": "N53"}, "row N54"), ({**GATE, "clean": False}, "clean tree"),
                      ({k: v for k, v in GATE.items() if k != "clean"}, "clean tree")):
        with pytest.raises(fr.Refused, match=msg):
            fr.run(cfg, fr.CONFIG, Path("/nonexistent"), Path("/nonexistent"), gate)
    with pytest.raises(fr.Refused, match="status is 'draft'"):
        fr.run(syn_cfg(status="draft"), fr.CONFIG, Path("/nonexistent"), Path("/nonexistent"), GATE)


@pytest.mark.parametrize("usable, cells_, want", [
    (True, {"CM": 50, "Endocardium": 40, "Endothelium": 60, "Epicardium": 30}, "recomputed"),
    (True, {"CM": 50, "Endocardium": 4, "Endothelium": 60, "Epicardium": 3}, "recomputed"),  # fallback; none
    (True, {"CM": 9, "Endocardium": 40, "Endothelium": 60, "Epicardium": 30}, "none"),  # CM has no direction
    (True, {"CM": 50, "Endocardium": 4, "Endothelium": 8, "Epicardium": 30}, "none"),  # Endocardium has none
    (False, {"CM": 50, "Endocardium": 40, "Endothelium": 60, "Epicardium": 30}, "none")])  # the cap failed
def test_the_terminal_route_is_fixed_by_the_rule(usable, cells_, want):
    fb = fr.fallbacks({g: {"E9.5": n, "E14.5": 100} for g, n in cells_.items()}, 10, COMMITTED["fallback"])
    assert fr.route(usable, fb, COMMITTED) == want


def test_the_build_command_per_route_is_the_first_uploads_command(bt):
    """Decisions (1) and (2): both routes build the first upload's command; route none is main's none file as built
    (build_all.sh's set t1none, row R8)."""
    assert fr.build_command("recomputed", COMMITTED) == ["t1", "--g-mom", "--base", "lib3", "--feng-groups",
                                                         "recomputed"]
    assert fr.build_command("none", COMMITTED) == ["t1", "--g-mom", "--base", "lib3", "--feng-groups", "none"]
    assert fr.FIRST_UPLOAD_COMMAND == bt.FENG_NONE_COMMAND == " ".join(fr.build_command("none", COMMITTED)[:4])
    with pytest.raises(fr.Refused, match="route must be"):
        fr.build_command("authors", COMMITTED)
    cfg = copy.deepcopy(COMMITTED)
    cfg["terminal"]["base"] = fr.PENDING
    with pytest.raises(fr.Refused, match="PENDING"):
        fr.build_command("none", cfg)
    build_all = (ROOT / "experiments/builds/build_all.sh").read_text()
    assert "build_test.py t1 --g-mom --base lib3 --feng-groups none --out" in build_all  # main's set t1none
    assert set(fr.BASE_ARGS) == {"lib3"}  # the mid base and its 0.60 are in no later file


# --- the gate before any run


@pytest.fixture
def gate_env(tmp_path, monkeypatch):
    """A frozen, sourced, ruled yaml at a stand-in CONFIG path, with TRIES, PLAN, HEAD and the tree faked."""
    from vec import provenance

    draft = copy.deepcopy(fr.without_unsourced(COMMITTED))  # the yaml as N55 read it (draft)
    draft_bytes = yaml.safe_dump(draft, sort_keys=False).encode()
    cfg = copy.deepcopy(draft)
    cfg["status"] = "frozen"
    path = tmp_path / "feng_relabel.yaml"
    path.write_text(yaml.safe_dump(cfg, sort_keys=False))
    tries = tmp_path / "TRIES.md"
    tries.write_text("| N53 | x |\n| N54 (diagnostic) | T1 |\n")
    plan = tmp_path / "PLAN.md"
    plan.write_text("- 2026-11-20 (a), the captain's ruling.\n  - Recomputed groups; N54 after it.\n")
    n55_out = tmp_path / "n55" / "n55.json"
    n55 = {"rule": {"outcome": "kept", "setting": {"min_score": 0.5, "min_margin": 0.25}},
           "inputs": {"config_sha256": fio.hashlib.sha256(draft_bytes).hexdigest()}, "stamp": {"commit": "n55c"}}
    head = {path.resolve(): path.read_bytes(), n55_out.resolve(): json.dumps(n55).encode()}
    at_commit = {"n55c": draft_bytes}
    monkeypatch.setattr(fr, "CONFIG", path)
    monkeypatch.setattr(fr, "TRIES", tries)
    monkeypatch.setattr(fr, "PLAN", plan)
    monkeypatch.setattr(fr, "RECORD", tmp_path / "results" / "feng_relabel" / "MANIFEST.json")
    monkeypatch.setattr(fr, "N55_OUT", n55_out)
    monkeypatch.setattr(fr, "committed_at_head", lambda p: head.get(Path(p).resolve()))
    monkeypatch.setattr(fr, "yaml_at_commit", lambda commit, path=None: at_commit.get(commit))
    monkeypatch.setattr(provenance, "git_state", lambda repo, ignore=(): ("abc", []))
    return SimpleNamespace(cfg=cfg, path=path, tries=tries, plan=plan, head=head, dst=tmp_path / "out",
                           provenance=provenance, monkeypatch=monkeypatch, n55=n55, n55_out=n55_out,
                           at_commit=at_commit, draft=draft)


def test_the_gate_passes_on_a_frozen_ruled_committed_clean_tree(gate_env):
    g = fr.preflight(gate_env.path, "2026-11-20 (a)", gate_env.dst)
    assert g == {"row": "N54", "ruling": "2026-11-20 (a)", "commit": "abc", "clean": True, "config_at_head": True,
                 "n55": N55_}


def test_the_gate_takes_the_setting_n55s_rule_gave(gate_env):
    """N55's rule changed the thresholds: the frozen yaml must carry the new setting, and only it."""
    env = gate_env
    env.n55["rule"] = {"outcome": "changed", "setting": {"min_score": 0.3, "min_margin": 0.25}}
    env.head[env.n55_out.resolve()] = json.dumps(env.n55).encode()
    with pytest.raises(fr.Refused, match="not the setting row N55's rule gave"):
        fr.preflight(env.path, "2026-11-20 (a)", env.dst)
    _rewrite(env, assign={**env.cfg["assign"], "min_score": 0.3})
    assert fr.preflight(env.path, "2026-11-20 (a)", env.dst)["n55"]["outcome"] == "changed"


def _rewrite(env, **over):
    cfg = copy.deepcopy(env.cfg)
    for k, v in over.items():
        if k == "terminal":
            cfg["terminal"].update(v)
        else:
            cfg[k] = v
    env.path.write_text(yaml.safe_dump(cfg, sort_keys=False))
    env.head[env.path.resolve()] = env.path.read_bytes()


def _output(env, folder: str) -> None:
    """The stand-in yaml with another output folder, in the version N55 read too (so check_n55 still passes)."""
    _rewrite(env, output={**env.cfg["output"], "dir": folder})
    draft = copy.deepcopy(env.draft)
    draft["output"] = {**draft["output"], "dir": folder}
    raw = yaml.safe_dump(draft, sort_keys=False).encode()
    env.at_commit["n55c"] = raw
    env.n55["inputs"]["config_sha256"] = fio.hashlib.sha256(raw).hexdigest()
    env.head[env.n55_out.resolve()] = json.dumps(env.n55).encode()


@pytest.mark.parametrize("case, msg", [
    ("other_config", "only the committed"), ("draft", "status"), ("unsourced", "UNSOURCED"),
    ("pending", "PENDING"), ("no_row", "N54 has no row"), ("no_ruling", "no decision-log entry"),
    ("wrong_ruling", "no decision-log entry"), ("not_at_head", "differs from its committed version"),
    ("dirty", "not clean"), ("git_unknown", "not clean"), ("record_in_tree", "already recorded"),
    ("record_at_head", "already recorded"), ("dst_exists", "already recorded"), ("started", "already recorded"),
    ("no_n55", "N55's dry run is not committed"), ("n55_other_value", "more than its status"),
    ("n55_yaml_missing", "cannot be found"), ("n55_yaml_other_sha", "cannot be found"),
    ("n55s_own_entry", "the entry that ordered row N55")])
def test_the_gate_refuses_before_anything_is_read(gate_env, tmp_path, case, msg):
    env, ruling, path = gate_env, "2026-11-20 (a)", gate_env.path
    if case == "other_config":
        path = tmp_path / "copy.yaml"
        path.write_bytes(env.path.read_bytes())
    elif case == "draft":
        _rewrite(env, status="draft")
    elif case == "unsourced":
        cfg = copy.deepcopy(env.cfg)
        cfg["sets"]["blood"]["genes"]["Alas2"] = fr.UNSOURCED
        env.path.write_text(yaml.safe_dump(cfg, sort_keys=False))
    elif case == "pending":
        _rewrite(env, terminal={"base": fr.PENDING})
    elif case == "no_row":
        env.tries.write_text("| N53 | x |\n")
    elif case == "no_ruling":
        ruling = None
    elif case == "wrong_ruling":
        ruling = "2026-11-21 (a)"
    elif case == "not_at_head":
        env.head[env.path.resolve()] = b"older"
    elif case == "dirty":
        env.monkeypatch.setattr(env.provenance, "git_state", lambda repo, ignore=(): ("abc", [" M x.py"]))
    elif case == "git_unknown":
        env.monkeypatch.setattr(env.provenance, "git_state", lambda repo, ignore=(): ("unknown", None))
    elif case == "record_in_tree":
        fr.RECORD.parent.mkdir(parents=True)
        fr.RECORD.write_text("{}")
    elif case == "record_at_head":
        env.head[fr.RECORD.resolve()] = b"{}"
    elif case == "dst_exists":
        env.dst.mkdir()
        (env.dst / "MANIFEST.json").write_text("{}")
    elif case == "started":  # a killed run left only its marker
        env.dst.mkdir()
        (env.dst / fr.STARTED).write_text("{}")
    elif case == "no_n55":
        del env.head[env.n55_out.resolve()]
    elif case == "n55_other_value":  # anything but status, min_score and min_margin moved after N55
        _rewrite(env, unassigned_cap=0.5)
    elif case == "n55_yaml_missing":
        env.at_commit.clear()
    elif case == "n55_yaml_other_sha":
        env.at_commit["n55c"] = b"other: yaml\n"
    elif case == "n55s_own_entry":  # the final check of the fix round, item D: N55's entry names N54 too
        env.n55["ruling"] = ruling
        env.head[env.n55_out.resolve()] = json.dumps(env.n55).encode()
    with pytest.raises(fr.Refused, match=msg):
        fr.preflight(path, ruling, env.dst)


def test_main_refuses_without_the_gate_and_outside_data(gate_env, monkeypatch, tmp_path):
    """The stand-in repository has its own real data/ directory, so the test holds whatever data/external and data/raw
    link to in the checkout running it (on the DGX a campaign checkout links them into the reference checkout)."""
    monkeypatch.setattr(fr, "run", lambda *a, **k: pytest.fail("nothing is read before the refusal"))
    repo = tmp_path / "repo"
    (repo / "data").mkdir(parents=True)
    monkeypatch.setattr(fr.fio, "REPO", repo)
    with pytest.raises(SystemExit, match="no decision-log entry"):
        fr.main(["--ruling", "2026-11-21 (a)"])
    for bad in ("../outside", "data/../elsewhere"):
        _output(gate_env, bad)
        with pytest.raises(SystemExit, match="must be inside"):
            fr.main(["--ruling", "2026-11-20 (a)"])
    with pytest.raises(SystemExit):  # --config is gone: no other yaml can be run
        fr.main(["--ruling", "2026-11-20 (a)", "--config", "x.yaml"])


def test_the_output_folder_must_resolve_inside_the_checkouts_own_data_directory(tmp_path):
    repo, elsewhere = tmp_path / "repo", tmp_path / "reference"
    (repo / "data").mkdir(parents=True)
    (elsewhere / "external").mkdir(parents=True)
    assert fr.output_dir_ok(repo / "data/derived/feng_relabel_c", repo) == (repo / "data/derived/feng_relabel_c"
                                                                            ).resolve()
    (repo / "data" / "external").symlink_to(elsewhere / "external")  # as a DGX campaign checkout's data/external
    with pytest.raises(fr.Refused, match="must be inside"):
        fr.output_dir_ok(repo / "data/external/relabel_c", repo)
    assert fr.output_dir_ok(repo / "data/derived/feng_relabel_c", repo)  # data/derived is real
    (repo / "data" / "derived").symlink_to(elsewhere)  # a linked data/derived is refused too
    with pytest.raises(fr.Refused, match="must be inside"):
        fr.output_dir_ok(repo / "data/derived/feng_relabel_c", repo)
    linked = tmp_path / "linked_repo"
    linked.mkdir()
    (linked / "data").symlink_to(repo / "data")  # a data/ that is itself a link is not the checkout's own
    with pytest.raises(fr.Refused, match="is a link"):
        fr.output_dir_ok(linked / "data/x", linked)


def test_main_writes_the_record_copy_after_the_run(gate_env, monkeypatch, tmp_path):
    seen = {}

    def fake_run(cfg, cfg_path, src, dst, gate):
        assert (dst / fr.STARTED).exists()  # the marker comes before any Feng file is opened
        (dst / "MANIFEST.json").write_text('{"gate": "x"}\n')
        seen["dst"], seen["gate"] = dst, gate
        return {"usable": True, "route": "recomputed", "unassigned_share": {}, "fallbacks": {}}

    monkeypatch.setattr(fr, "run", fake_run)
    monkeypatch.setattr(fr.fio, "REPO", tmp_path)  # data/ under the stand-in repository
    _output(gate_env, "data/derived/relabel_c")
    assert fr.main(["--ruling", "2026-11-20 (a)"]) == 0
    assert seen["gate"]["ruling"] == "2026-11-20 (a)"
    assert fr.RECORD.read_bytes() == (seen["dst"] / "MANIFEST.json").read_bytes()


@pytest.mark.parametrize("stop", ["refused", "error", "missing_prox1"])
def test_a_stop_after_preflight_takes_route_none_and_blocks_a_rerun(gate_env, monkeypatch, tmp_path, capsys, stop):
    """Blocker 1 of the final review, decision (3) of 2026-11-19: a run that stops after preflight writes the stop
    record (manifest and RECORD: reason, usable false, route none), prints the none command and exits STOP_EXIT; the
    gate then refuses a second run; no reader takes a file from it; the listing lists nothing."""
    def fake_run(cfg, cfg_path, src, dst, gate):
        if stop == "refused":
            raise fr.Refused("set lymphatic: 1 of its genes are in Feng's symbols (unique); the yaml needs 2")
        if stop == "error":
            raise MemoryError("synthetic")
        C, genes, stage, _ = _committed_extremes({"x": {}}, n_per=40, drop=("Prox1",))
        fr.relabel(C, genes, stage, cfg)  # the real refusal of a class at the minimum with a gene missing

    monkeypatch.setattr(fr, "run", fake_run)
    monkeypatch.setattr(fr.fio, "REPO", tmp_path)
    monkeypatch.setattr(fr.fio, "provenance", lambda code: {"git_commit": "abc", "git_uncommitted": []})
    _output(gate_env, "data/derived/relabel_c")
    assert fr.main(["--ruling", "2026-11-20 (a)"]) == fr.STOP_EXIT
    dst = tmp_path / "data/derived/relabel_c"
    man = json.loads((dst / "MANIFEST.json").read_text())
    assert man["stopped"] is True and man["usable"] is False and man["route"] == "none" and man["outputs"] == {}
    assert man["gate"]["ruling"] == "2026-11-20 (a)" and man["config"]["status"] == "frozen"
    assert ("lymphatic" in man["reason"]) == (stop != "error")
    assert fr.RECORD.read_bytes() == (dst / "MANIFEST.json").read_bytes()
    assert "route none" in capsys.readouterr().out and (dst / fr.STARTED).exists()
    with pytest.raises(SystemExit, match="already recorded"):  # nothing is rerun under N54
        fr.main(["--ruling", "2026-11-20 (a)"])
    assert "stopped after preflight" in n54.stopped(dst / "MANIFEST.json")
    (dst / "cd1_stats_relabel.npz").write_bytes(b"partial")  # even a partial file is not readable through it
    with pytest.raises(ValueError, match="is not listed"):
        feng_dirs.check_recomputed(dst / "cd1_stats_relabel.npz")


def test_an_interrupt_after_preflight_is_recorded_and_raised(gate_env, monkeypatch, tmp_path):
    def fake_run(*a, **k):
        raise KeyboardInterrupt

    monkeypatch.setattr(fr, "run", fake_run)
    monkeypatch.setattr(fr.fio, "REPO", tmp_path)
    monkeypatch.setattr(fr.fio, "provenance", lambda code: {})
    _output(gate_env, "data/derived/relabel_c")
    with pytest.raises(KeyboardInterrupt):
        fr.main(["--ruling", "2026-11-20 (a)"])
    assert json.loads(fr.RECORD.read_text())["route"] == "none"


# --- the builder's routes


def test_builder_refuses_recomputed_unless_the_committed_record_routes_to_it(bt, t1_env, monkeypatch, tmp_path):  # noqa: F811
    """Refused before any T1 file is read: no record at HEAD; a record whose route is none; then, with a record that
    routes to recomputed, a recomputed file whose own route is none (feng_dirs' answer)."""
    out = str(tmp_path / "x.h5ad")
    for record_route, msg in ((None, "no permitted-only relabelling is committed"),
                              ("none", "committed relabelling's route is 'none'")):
        _recomputed_env(bt, t1_env, monkeypatch, record_route=record_route)
        monkeypatch.setattr(bt, "load_stage", lambda *a, **k: pytest.fail("refused before any T1 file is read"))
        with pytest.raises(SystemExit, match=msg):
            bt.main(["t1", *TARGET, "--feng-groups", "recomputed", "--out", out])
    _recomputed_env(bt, t1_env, monkeypatch, groups={"Endocardium", "Endothelium", "Epicardium"}, route="none")
    monkeypatch.setattr(bt, "load_stage", lambda *a, **k: pytest.fail("refused before any T1 file is read"))
    with pytest.raises(SystemExit, match="routes this relabelling to 'none'"):
        bt.main(["t1", *TARGET, "--feng-groups", "recomputed", "--out", out])


@pytest.mark.parametrize("argv", [
    ["--base", "lib3"], [], ["--g-mom"], ["--g-mom", "--base", "pooled"],
    ["--g-mom", "--base", "mid", "--cm-share", "0.60"], [*TARGET, "--cm-share", "0.60"],
    [*TARGET, "--blood-share", "0.10"], [*TARGET, "--lib3-pair", "--cm-share", "0.60"], [*TARGET, "--herd"],
    [*TARGET, "--no-hold-unbracketed"], [*TARGET, "--no-feng-hygiene"], [*TARGET, "--strong-null"],
    [*TARGET, "--feng-none-ruling", "2026-11-20 (a)"]])
def test_recomputed_builds_only_the_first_uploads_command(bt, monkeypatch, tmp_path, argv):
    for obj, name in ((bt, "load_stage"), (bt, "load_boards"), (bt.hygiene, "load_union"),
                      (bt.feng_dirs, "feng_directions"), (fr, "committed_at_head")):
        monkeypatch.setattr(obj, name, lambda *a, **k: pytest.fail("nothing is read before the refusal"))
    with pytest.raises(SystemExit, match="recomputed builds one file: the first upload's command"):
        bt.main(["t1", *argv, "--feng-groups", "recomputed", "--out", str(tmp_path / "x.h5ad")])


def test_recomputed_with_epicardium_missing_moves_none_of_its_cells_and_g_mom_edits_the_entrys_labels(
        bt, t1_env, monkeypatch, tmp_path):  # noqa: F811
    """Epicardium without a recomputed direction (its fallback is none): its cells keep their released counts (R1
    holds them); g_mom edits cm and endocardium cells whose group has a direction (V-CM, Endocardium-1, EndoMT
    Endocardium) and leaves Early VIC alone, as in the first upload."""
    from test_t1_feng_groups_none import EDITED, _edited_labels, _released

    _recomputed_env(bt, t1_env, monkeypatch, groups={"CM", "Endocardium", "Endothelium"})
    seen = []
    real = bt.t1_scope.build_submission

    def capture(*a, **k):
        hook = k["post_cells"]
        entry = {"kwargs": k}

        def wrapped(pop):
            entry["before"] = (pop.copy(), np.asarray(pop.layers["counts"].todense()))
            out_ = hook(pop)
            entry["after"] = (out_.copy(), np.asarray(out_.layers["counts"].todense()))
            return out_

        seen.append(entry)
        return real(*a, **{**k, "post_cells": wrapped})

    monkeypatch.setattr(bt.t1_scope, "build_submission", capture)
    out = tmp_path / "T1_rc.h5ad"
    bt.main(["t1", *TARGET, "--feng-groups", "recomputed", "--out", str(out)])
    assert set(seen[0]["kwargs"]["directions"]) == {"CM", "Endocardium", "Endothelium"}
    assert _edited_labels(seen[0]) == EDITED
    pop, before = seen[0]["before"]
    labels = pop.obs["celltype"].astype(str).to_numpy()
    epi = labels == "Epicardium"
    assert epi.any() and np.array_equal(before[epi], _released(t1_env, pop)[epi])  # no direction: as drawn
    rec = json.loads(out.with_suffix(".build.json").read_text())
    assert rec["feng_groups"]["with_direction"] == ["CM", "Endocardium", "Endothelium"]
    assert "['CM', 'Endocardium', 'Endothelium']" in rec["g_mom"]["mask"]


def test_g_mom_under_recomputed_holds_back_the_labels_whose_group_has_no_direction(
        bt, t1_env, monkeypatch, tmp_path):  # noqa: F811
    """The mask follows the groups with a direction (here only CM, which the terminal rule would route to none, so
    this checks the mask's arithmetic alone): only V-CM is edited, and the covered-label check agrees."""
    from test_t1_feng_groups_none import _edited_labels

    _recomputed_env(bt, t1_env, monkeypatch, groups={"CM"})
    seen = []
    real = bt.t1_scope.build_submission

    def capture(*a, **k):
        hook = k["post_cells"]
        entry = {}

        def wrapped(pop):
            entry["before"] = (pop.copy(), np.asarray(pop.layers["counts"].todense()))
            out_ = hook(pop)
            entry["after"] = (out_.copy(), np.asarray(out_.layers["counts"].todense()))
            return out_

        seen.append(entry)
        return real(*a, **{**k, "post_cells": wrapped})

    monkeypatch.setattr(bt.t1_scope, "build_submission", capture)
    bt.main(["t1", *TARGET, "--feng-groups", "recomputed", "--out", str(tmp_path / "T1_cm.h5ad")])
    assert _edited_labels(seen[0]) == {"V-CM"}
    assert bt._labels_without_direction(["V-CM", "Endocardium-1", "BEC", "Early VIC"], ["CM", "Endothelium"]) == {
        "Endocardium-1", "Early VIC"}


def _none_env(bt, t1_env, monkeypatch, record, where):  # noqa: F811
    """main's none route on the fixture, with the labeller's RECORD in the tree as `record` (None: no relabelling)."""
    from test_g_mom_verify import _g_mom_loads
    from test_t1_feng_groups_none import _lib3, _real_none

    _lib3(t1_env)
    _real_none(bt, t1_env, monkeypatch)
    _g_mom_loads(bt, t1_env, monkeypatch)
    rec_path = Path(where) / "relabel_record" / "MANIFEST.json"
    monkeypatch.setattr(fr, "RECORD", rec_path)
    monkeypatch.setattr(fr, "committed_at_head", lambda p: None)
    if record is not None:
        rec_path.parent.mkdir(parents=True, exist_ok=True)
        rec_path.write_text(json.dumps(record))
    plan = Path(where) / "PLAN.md"
    plan.write_text("- 2026-11-20 (a), the captain's ruling for the run; row N54.\n"
                    "- 2026-11-22 (a), the captain sets the N54 route aside for T1.\n"
                    "- 2026-11-22 (b), the listing of row N54 is recorded.\n")
    monkeypatch.setattr(fr, "PLAN", plan)
    # the record was committed on 2026-11-21, when docs/PLAN.md held the run's ruling only
    monkeypatch.setattr(fr, "last_commit", lambda path: ("rc1", "2026-11-21"))
    monkeypatch.setattr(fr, "plan_at_commit",
                        lambda commit: "- 2026-11-20 (a), the captain's ruling for the run; row N54.\n")


def test_none_before_any_relabelling_builds_as_on_main(bt, t1_env, monkeypatch, tmp_path):  # noqa: F811
    _none_env(bt, t1_env, monkeypatch, None, tmp_path)
    out = tmp_path / "T1_none.h5ad"
    bt.main(["t1", *TARGET, "--feng-groups", "none", "--out", str(out)])
    rec = json.loads(out.with_suffix(".build.json").read_text())
    assert set(rec["feng_groups"]) == {"groups", "file", "sha256", "direction_step"}  # main's record, nothing added
    with pytest.raises(SystemExit, match="--feng-none-ruling is refused: no relabelling routes to recomputed"):
        bt.main(["t1", *TARGET, "--feng-groups", "none", "--feng-none-ruling", "2026-11-22 (a)", "--out", str(out)])
    with pytest.raises(SystemExit, match="--feng-none-ruling is refused without --feng-groups none"):
        bt.main(["t1", *TARGET, "--feng-none-ruling", "2026-11-22 (a)", "--out", str(out)])


def test_none_after_a_relabelling_that_routes_to_none_builds_as_on_main(bt, t1_env, monkeypatch, tmp_path):  # noqa: F811
    _none_env(bt, t1_env, monkeypatch, {"route": "none", "stopped": True}, tmp_path)
    out = tmp_path / "T1_none.h5ad"
    bt.main(["t1", *TARGET, "--feng-groups", "none", "--out", str(out)])
    assert "ruling_over_the_relabelling_route" not in json.loads(out.with_suffix(".build.json").read_text())[
        "feng_groups"]


def test_none_after_a_relabelling_that_routes_to_recomputed_needs_the_captains_ruling(
        bt, t1_env, monkeypatch, tmp_path):  # noqa: F811
    """The final review's lesser point 1, with main's none: refused without --feng-none-ruling, refused with an entry
    that does not exist or does not name N54, built with the captain's entry, and the record names it. The file's
    bytes do not depend on it."""
    _none_env(bt, t1_env, monkeypatch, None, tmp_path)
    plain = tmp_path / "T1_none_plain.h5ad"
    bt.main(["t1", *TARGET, "--feng-groups", "none", "--out", str(plain)])
    fr.RECORD.parent.mkdir(parents=True, exist_ok=True)
    fr.RECORD.write_text(json.dumps({"route": "recomputed"}))  # a relabelling now exists and routes to recomputed
    out = tmp_path / "T1_none.h5ad"
    for extra in ([], ["--feng-none-ruling", "2026-11-23 (a)"], ["--feng-none-ruling", "2026-11-22 (b)"]):
        with pytest.raises(SystemExit, match="routes T1's file to recomputed"):
            bt.main(["t1", *TARGET, "--feng-groups", "none", *extra, "--out", str(out)])
    with pytest.raises(SystemExit, match="does not contain the fixed phrase 'sets the N54 route aside'"):
        bt.main(["t1", *TARGET, "--feng-groups", "none", "--feng-none-ruling", "2026-11-22 (b)", "--out", str(out)])
    assert not out.exists()
    bt.main(["t1", *TARGET, "--feng-groups", "none", "--feng-none-ruling", "2026-11-22 (a)", "--out", str(out)])
    rec = json.loads(out.with_suffix(".build.json").read_text())
    assert rec["feng_groups"]["ruling_over_the_relabelling_route"] == "2026-11-22 (a)"
    assert out.read_bytes() == plain.read_bytes()


def test_feng_directions_none_reads_no_file():
    d, rec = feng_dirs.feng_directions(["g1"], groups="none")
    assert d == {} and rec["feng_groups"] == "none" and rec["file"] is None


def test_builder_catches_a_labeller_refusal(bt, monkeypatch):
    def refuse(genes, *a, **k):
        raise fr.Refused("malformed yaml")

    monkeypatch.setattr(bt.feng_dirs, "feng_directions", refuse)
    with pytest.raises(SystemExit, match="--feng-groups recomputed is refused: malformed yaml"):
        bt._feng_directions(["g"], "recomputed")


def test_builder_blood_guard_reads_the_recomputed_group_names(bt, monkeypatch):
    monkeypatch.setattr(bt.hygiene, "load_union", lambda *a, **k: {"genes": []})
    monkeypatch.setattr(bt.feng_dirs, "feng_directions", lambda genes, *a, **k: (
        {"CM": np.zeros(1), "Blood": np.zeros(1)}, {"route": "recomputed"}))
    with pytest.raises(SystemExit, match="Blood is among"):
        bt._t1_union(True, "recomputed")
    monkeypatch.setattr(bt.feng_dirs, "feng_directions", lambda genes, *a, **k: (
        {"CM": np.zeros(1)}, {"route": "recomputed"}))
    assert bt._t1_union(True, "recomputed") == {"genes": []}


def test_the_output_paths_agree():
    """The yaml's output folder, sources.yaml's feng_relabel entry and the harness's late-proxy path name one place."""
    src = yaml.safe_load((ROOT / "experiments/external/sources.yaml").read_text())["feng_relabel"]
    assert src["dst"] == COMMITTED["output"]["dir"] == "data/derived/feng_relabel_c"
    assert harness.LP_STATS_RECOMPUTED == ROOT / COMMITTED["output"]["dir"] / COMMITTED["output"]["stats"]


# --- the listing, the agreement and the dry run


def test_listing_needs_the_record_committed(tmp_path, monkeypatch):
    man = tmp_path / "MANIFEST.json"
    man.write_text('{"a": 1}\n')
    monkeypatch.setattr(fr, "committed_at_head", lambda p: None)
    assert not n54.record_committed(man)
    monkeypatch.setattr(fr, "committed_at_head", lambda p: b'{"a": 2}\n')
    assert not n54.record_committed(man)
    monkeypatch.setattr(fr, "committed_at_head", lambda p: man.read_bytes())
    assert n54.record_committed(man)


def test_listing_main_refuses_without_row_or_ruling(tmp_path, monkeypatch):
    import n12_hygiene_size as n12

    monkeypatch.setattr(n12, "tries_rows", lambda *a, **k: set())
    with pytest.raises(SystemExit, match="N54 has no row"):
        n54.main(["--ruling", "2026-11-20 (a)", "--recomputed-sha256", "x"])
    monkeypatch.setattr(n12, "tries_rows", lambda *a, **k: {"N54"})
    plan = tmp_path / "PLAN.md"
    plan.write_text("- 2026-11-20 (a), something else\n")
    monkeypatch.setattr(n54, "PLAN", plan)
    with pytest.raises(SystemExit, match="no decision-log entry"):
        n54.main(["--ruling", "2026-11-20 (a)", "--recomputed-sha256", "x"])


def test_agreement_main_refuses_before_the_listing_is_committed(tmp_path, monkeypatch):
    import n54_authors_agreement as ag

    listing = tmp_path / "n54.json"
    monkeypatch.setattr(ag, "OUT", tmp_path / "agree")
    monkeypatch.setattr(ag, "LISTING", listing)
    args = ["--frozen-sha256", "x", "--method-report", "--upload-ruling", "2026-11-23 (a)"]
    with pytest.raises(SystemExit, match="has not run"):
        ag.main(args)
    listing.write_text('{"inputs": {"recomputed": {"sha256": "x"}}}\n')
    monkeypatch.setattr(fr, "committed_at_head", lambda p: None)
    with pytest.raises(SystemExit, match="not committed at HEAD"):
        ag.main(args)
    with pytest.raises(SystemExit):  # the upload ruling is a required argument
        ag.main(["--frozen-sha256", "x", "--method-report"])


def test_the_released_dry_run_on_synthetic_counts():
    import feng_relabel_dryrun as dr

    cfg = syn_cfg()
    types = cells(SPEC)
    name = {"cm": "V-CM", "endo": "Endocardium", "vasc": "BEC", "ambig": "BEC", "epi": "Proepicardium",
            "blood": "Blood", "none": "ST", "doublet": "V-CM"}
    released = np.array([name[t] for t in types])
    rep = dr.dry_run(counts(types, seed=11), GENES, released, cfg)
    assert set(rep["by_released_label"]) == set(released.tolist())
    assert rep["by_released_label"]["V-CM"]["n"] == SPEC["cm"] + SPEC["doublet"]
    assert rep["by_released_label"]["Blood"]["shares"].get("r_blood", 0) >= 0.9
    assert set(rep["cap"]["shares"]) == {"E9.5"} and set(rep["under_min_cells"]) == set(cfg["groups"])
    assert len(rep["sensitivity"]) == len(cfg["sensitivity"])
    for s in rep["sensitivity"]:  # each setting: counts, shares per released label and the cap's verdict
        assert set(s) == {"setting", "group_cells", "unassigned_share", "cap", "by_released_label"}
        assert set(s["by_released_label"]) == set(released.tolist()) and "passed" in s["cap"]
    assert rep["rule"]["outcome"] == "kept" and rep["rule"]["setting"] == cfg["n55_rule"]["order"][0]
    assert set(rep["rule"]["settings"][0]["chambers"]) == {"CM", "Endocardium", "BEC", "Proepicardium"}
    json.dumps(rep, default=float)
    code = (ROOT / "experiments/external/feng_relabel_dryrun.py").read_text().split('"""', 2)[2]
    assert not re.search(r"feng2022|cd1_|cell_type_spec", code)  # the dry run's code names no Feng file


def test_n55s_rule_keeps_changes_or_records_a_failed_guard():
    """The rule written before N55 runs (feng_relabel.n55_outcome), on hand-made scores: (i) the first setting passes:
    kept; (ii) it fails and a later one passes: changed to the first that passes, in the yaml's order; (iii) none
    passes: the first is kept and the outcome says the guard failed."""
    cfg = syn_cfg()
    released = np.array(["V-CM", "V-CM", "V-CM", "Endocardium", "BEC", "Proepicardium"])
    zero = np.zeros(6)

    def sc(cm, epi=None, ec=None, endo=None, vasc=None):
        return {"cm": np.array(cm, float), "ec": np.array(ec if ec is not None else [0, 0, 0, 0.9, 0.9, 0], float),
                "epi": np.array(epi if epi is not None else [0, 0, 0, 0, 0, 0.9], float),
                "endo_split": np.array(endo if endo is not None else [0, 0, 0, 1, 0, 0], float),
                "vasc_split": np.array(vasc if vasc is not None else [0, 0, 0, 0, 1, 0], float), "blood": zero}

    kept = fr.n55_outcome(released, sc([0.9, 0.9, 0.9, 0, 0, 0]), cfg)
    assert kept["outcome"] == "kept" and kept["setting"] == {"min_score": 0.5, "min_margin": 0.25}
    # the CM cells score 0.4: under min_score 0.5 they fall to r_unassigned; at 0.3 they pass
    changed = fr.n55_outcome(released, sc([0.4, 0.4, 0.4, 0, 0, 0]), cfg)
    assert changed["outcome"] == "changed" and changed["setting"] == {"min_score": 0.3, "min_margin": 0.25}
    assert not changed["settings"][0]["chambers"]["CM"]["passed"] and changed["settings"][1]["passed"]
    # the margin: CM 0.9 against epi 0.7 (margin 0.2) passes only at min_margin 0.15
    margin = fr.n55_outcome(released, sc([0.9, 0.9, 0.9, 0, 0, 0], epi=[0.7, 0.7, 0.7, 0, 0, 0.9]), cfg)
    assert margin["outcome"] == "changed" and margin["setting"] == {"min_score": 0.5, "min_margin": 0.15}
    failed = fr.n55_outcome(released, sc([0.1, 0.1, 0.1, 0, 0, 0]), cfg)
    assert failed["outcome"] == "guard_failed" and failed["setting"] == {"min_score": 0.5, "min_margin": 0.25}
    # a chamber label absent from the released file fails, and the cap counts too
    c = fr.chamber_check(np.array(["V-CM"]), np.array(["r_cm"]), cfg["n55_rule"]["chamber_labels"])
    assert c["CM"]["passed"] and not c["Endocardium"]["passed"] and c["Endocardium"]["n"] == 0
    tie = fr.chamber_check(np.array(["BEC", "BEC"]), np.array(["r_vascular_ec", "r_unassigned"]),
                           cfg["n55_rule"]["chamber_labels"])
    assert not tie["BEC"]["passed"]  # half and half: not "most often"
    pooled = fr.chamber_check(np.array(["BEC"] * 3), np.array(["r_vascular_ec", "r_ec_ambiguous", "r_unassigned"]),
                              cfg["n55_rule"]["chamber_labels"])
    assert pooled["BEC"]["passed"] and pooled["BEC"]["own_share"] == pytest.approx(2 / 3)


def test_the_dry_runs_refusals(monkeypatch, tmp_path):
    import feng_relabel_dryrun as dr

    monkeypatch.setattr(fr, "row_in_tries", lambda row, path=None: row == "N55")
    monkeypatch.setattr(fr, "ruling_logged", lambda ruling, row="N54", plan=None: ruling == "2026-11-20 (b)")
    monkeypatch.setattr(fr, "N55_OUT", tmp_path / "n55.json")
    monkeypatch.setattr(fr, "committed_at_head", lambda p: None)
    cfg = copy.deepcopy(COMMITTED)
    assert dr.refusals(cfg, "2026-11-20 (b)") is None
    assert "naming N55" in dr.refusals(cfg, "2026-11-20 (c)")
    assert "before the freeze" in dr.refusals({**cfg, "status": "frozen"}, "2026-11-20 (b)")
    other = copy.deepcopy(cfg)
    other["assign"]["min_score"] = 0.3
    assert "first setting" in dr.refusals(other, "2026-11-20 (b)")
    monkeypatch.setattr(fr, "committed_at_head", lambda p: b"{}")  # deleting an untracked output allows no rerun
    assert "already written" in dr.refusals(cfg, "2026-11-20 (b)")
    monkeypatch.setattr(fr, "committed_at_head", lambda p: None)
    (tmp_path / "n55.json").write_text("{}")
    assert "already written" in dr.refusals(cfg, "2026-11-20 (b)")


# ------------------------------------------------------------------ the committed sets on synthetic extremes


def _committed_extremes(types_spec: dict[str, dict[str, float]], n_per: int = 20, seed: int = 0, drop: tuple = ()):
    """Poisson counts over the committed yaml's set genes, the three dropped genes and 2000 filler genes: every gene at
    its base rate (set genes low, filler spread like real genes), a type's listed genes at the listed rate."""
    cfg = COMMITTED
    members = [g for s in cfg["sets"].values() for g in s["genes"]]
    named = [g for g in members + sorted(cfg["dropped"]) if g not in drop]
    genes = named + [f"f{i:04d}" for i in range(2000)]
    pos = {g: i for i, g in enumerate(genes)}
    base = np.where(np.arange(len(genes)) < len(named), 0.2, np.random.default_rng(99).gamma(0.5, 4.0, len(genes)))
    types = [t for t in types_spec for _ in range(n_per)]
    rate = np.tile(base, (len(types), 1))
    for i, t in enumerate(types):
        for g, v in types_spec[t].items():
            if g in pos:
                rate[i, pos[g]] = v
    C = np.random.default_rng(seed).poisson(rate).astype(np.int32)
    C[:, -1] += 1
    stage = np.array(["E9.5" if i % 2 else "E14.5" for i in range(len(types))])
    return sparse.csr_matrix(C), genes, stage, np.array(types)


def test_the_committed_sets_on_synthetic_extremes():
    """The changed sets (literature genes in; Phox2b, Fcer1g and Pdpn out) under the unchanged score and assignment.
    Synthetic: what it shows is the rule's arithmetic, not how Feng's cells will fall."""
    cfg = COMMITTED
    G = {n: list(s["genes"]) for n, s in cfg["sets"].items()}
    H = 30.0

    def hi(gs):
        return {g: H for g in gs}

    spec = {
        "cm": hi(G["cm"]), "endo": hi(G["ec"] + G["endo_split"]), "vasc": hi(G["ec"] + G["vasc_split"]),
        "epi": hi(G["epi"]),
        "epi_pdpn": hi(G["epi"] + ["Pdpn"]),  # Pdpn is out: an epicardial cell expressing it stays epicardial
        "blood_definitive": hi(["Hba-a1", "Hbb-bs", "Alas2"]),  # E14.5-like red cells without embryonic globins
        "neural_crest": hi(["Sox10", "Foxd3", "Tfap2a", "Dlx2", "Dlx5"]),
        "fibroblast_pdgfrb": hi(["Col1a1", "Col1a2", "Col3a1", "Dcn", "Postn", "Pdgfrb"]),  # lit_B's caveat
        "pericyte": hi(["Pdgfrb", "Rgs5", "Kcnj8", "Cspg4"]), "smooth_muscle": hi(["Myh11", "Pdgfrb"]),
        "macrophage": hi(["Ptprc", "C1qb", "Lyz2", "Cd68"]),
        "lec": {**{g: 8.0 for g in G["ec"]}, "Prox1": H, "Flt4": H},
        "lec_full_ec": {**hi(G["ec"]), "Prox1": H, "Flt4": H},
        "endo_flt4": {**hi(G["ec"] + G["endo_split"]), "Flt4": H},  # lit_J's caveat on Flt4
        "valve_endo_prox1": {**hi(G["ec"] + G["endo_split"]), "Prox1": H},  # lit_J's caveat on Prox1
        "glia_plp1_only": hi(["Plp1"]),
        "none": {}}
    C, genes, stage, types = _committed_extremes(spec)
    res = fr.relabel(C, genes, stage, cfg)
    lab = res["labels"]
    want = {"cm": "r_cm", "endo": "r_endocardial_ec", "vasc": "r_vascular_ec", "epi": "r_epicardial",
            "epi_pdpn": "r_epicardial", "blood_definitive": "r_blood", "neural_crest": "r_neural_crest",
            "fibroblast_pdgfrb": "r_mesenchyme", "pericyte": "r_mural", "smooth_muscle": "r_mural",
            "macrophage": "r_immune", "lec": "r_lymphatic", "endo_flt4": "r_endocardial_ec",
            "valve_endo_prox1": "r_endocardial_ec"}
    for t, w in want.items():
        m = types == t
        assert np.mean(lab[m] == w) >= 0.9, (t, w, np.unique(lab[m], return_counts=True))
    in_groups = {lab_ for labs in cfg["groups"].values() for lab_ in labs}
    # a lymphatic cell with the pan-endothelial genes as high as its own two: never in a group (lymphatic or unassigned)
    assert not set(lab[types == "lec_full_ec"]) & in_groups
    # a cell carrying one late gene alone (Plp1) scores low everywhere: mostly neural crest or unassigned; on this seed
    # one cell in twenty falls into Epicardium by noise on the three-gene epi set (reported, not tuned)
    assert np.mean(np.isin(lab[types == "glia_plp1_only"], list(in_groups))) <= 0.1
    assert np.mean(np.isin(lab[types == "none"], list(in_groups))) <= 0.1
    # the dropped genes are never drawn as controls
    assert not {g for s in res["sets"].values() for g in s["controls"]} & set(cfg["dropped"])


def test_the_lymphatic_class_at_the_minimum_stops_the_run_if_a_gene_is_missing():
    """lymphatic holds exactly min_genes_found (2) genes since Pdpn was dropped: if Feng's symbols lack Prox1 or Flt4,
    the run stops before any label is written (no class is dropped after the data are read)."""
    for gone in ("Prox1", "Flt4"):
        C, genes, stage, _ = _committed_extremes({"x": {}}, n_per=40, drop=(gone,))
        with pytest.raises(fr.Refused, match="set lymphatic"):
            fr.relabel(C, genes, stage, COMMITTED)


# ------------------------------------------------------------------ the review of the port (second fixer, 2026-11-19)


def _dry_env(monkeypatch, tmp_path, cfg, C, genes, released):
    """feng_relabel_dryrun.main on a stand-in repository: the released file's reader returns synthetic counts."""
    import feng_relabel_dryrun as dr
    import pandas as pd

    from vec import core, prep, provenance

    (tmp_path / "data").mkdir(exist_ok=True)
    monkeypatch.setattr(fr.fio, "REPO", tmp_path)
    monkeypatch.setattr(dr, "RELEASED", tmp_path / "data/raw/T1/E9.5_RNA.h5ad")
    monkeypatch.setattr(fr, "N55_OUT", tmp_path / "experiments/external/results/feng_relabel_dryrun/n55.json")
    monkeypatch.setattr(fr, "N55_STARTED", tmp_path / "data/derived/feng_relabel_n55/STARTED.json")
    monkeypatch.setattr(fr, "load_config", lambda path=None: copy.deepcopy(cfg))
    monkeypatch.setattr(fr, "row_in_tries", lambda row, path=None: row == "N55")
    monkeypatch.setattr(fr, "ruling_logged", lambda ruling, row="N54", plan=None, text=None: ruling == "2026-11-20 (b)")
    monkeypatch.setattr(fr, "committed_at_head", lambda p: None)
    monkeypatch.setattr(provenance, "git_state", lambda repo, ignore=(): ("n55c", []))
    ref = SimpleNamespace(var_names=list(genes), obs=pd.DataFrame({"celltype": np.asarray(released)}))
    monkeypatch.setattr(prep, "load_stage", lambda *a, **k: ref)
    monkeypatch.setattr(core, "counts_of", lambda r: C)
    monkeypatch.setattr(fr.fio, "sha256_file", lambda path: "sha-" + Path(path).name)
    return dr


@pytest.mark.parametrize("stop", ["missing_prox1", "no_control"])
def test_a_dry_run_that_stops_is_recorded_with_the_first_setting_and_never_rerun(monkeypatch, tmp_path, capsys, stop):
    """The review of the port, its blocker: n55_rule (iv). A dry run that stops (a set short of two genes in the
    released symbols, a set that draws no control gene) writes n55.json with outcome stopped and the first setting,
    exits STOP_EXIT, and a second call is refused."""
    if stop == "missing_prox1":
        cfg = copy.deepcopy(COMMITTED)
        C, genes, _, _ = _committed_extremes({"x": {}}, n_per=40, drop=("Prox1",))
        want = "set lymphatic"
    else:
        cfg = syn_cfg()
        types = cells(SPEC)
        C = counts(types, seed=2).tolil()
        for g in ("Bl1", "Bl2"):  # every gene of the blood set has a zero mean: no bin, so no control gene
            C[:, GENES.index(g)] = 0
        C, genes = sparse.csr_matrix(C), GENES
        want = "drew no control gene"
    dr = _dry_env(monkeypatch, tmp_path, cfg, C, genes, ["V-CM"] * C.shape[0])
    assert dr.main(["--ruling", "2026-11-20 (b)"]) == fr.STOP_EXIT
    rep = json.loads(fr.N55_OUT.read_text())
    assert rep["stopped"] is True and want in rep["reason"]
    assert rep["rule"] == {"outcome": "stopped", "setting": cfg["n55_rule"]["order"][0]}
    assert rep["stamp"]["commit"] == "n55c" and rep["inputs"]["config_sha256"]
    assert "outcome stopped" in capsys.readouterr().out and fr.N55_STARTED.exists()
    with pytest.raises(SystemExit, match="already written"):
        dr.main(["--ruling", "2026-11-20 (b)"])
    fr.N55_OUT.unlink()  # even with the output gone, the marker blocks a silent rerun
    with pytest.raises(SystemExit, match="was started before"):
        dr.main(["--ruling", "2026-11-20 (b)"])


def test_an_interrupted_dry_run_is_recorded_and_raised_and_a_killed_one_recorded_afterwards(monkeypatch, tmp_path):
    from vec import core

    cfg = syn_cfg()
    dr = _dry_env(monkeypatch, tmp_path, cfg, counts(cells(SPEC)), GENES, ["V-CM"] * sum(SPEC.values()))

    def interrupt(r):
        raise KeyboardInterrupt

    monkeypatch.setattr(core, "counts_of", interrupt)
    with pytest.raises(KeyboardInterrupt):
        dr.main(["--ruling", "2026-11-20 (b)"])
    assert json.loads(fr.N55_OUT.read_text())["rule"]["outcome"] == "stopped"
    with pytest.raises(SystemExit, match="n55.json exists"):  # nothing to record: the stop is on file
        dr.main(["--record-kill"])
    # a kill leaves the marker alone: --record-kill writes the stop record from it
    fr.N55_OUT.unlink()
    assert dr.main(["--record-kill"]) == fr.STOP_EXIT
    rep = json.loads(fr.N55_OUT.read_text())
    assert rep["rule"]["outcome"] == "stopped" and rep["reason"].startswith("killed")
    assert rep["ruling"] == "2026-11-20 (b)" and rep["stamp"]["commit"] == "n55c"
    fr.N55_STARTED.unlink()
    fr.N55_OUT.unlink()
    with pytest.raises(SystemExit, match="no started dry run"):
        dr.main(["--record-kill"])


def test_the_gate_accepts_a_stopped_n55_with_the_first_setting_and_refuses_a_setting_its_outcome_cannot_give(
        gate_env):
    env, ruling = gate_env, "2026-11-20 (a)"

    def n55(outcome, setting, **extra):
        env.n55["rule"] = {"outcome": outcome, "setting": setting}
        env.n55.update(extra)
        env.head[env.n55_out.resolve()] = json.dumps(env.n55).encode()

    first, later = {"min_score": 0.5, "min_margin": 0.25}, {"min_score": 0.5, "min_margin": 0.15}
    n55("stopped", first, stopped=True, reason="Refused: set lymphatic: 1 of its genes")
    g = fr.preflight(env.path, ruling, env.dst)["n55"]
    assert g["outcome"] == "stopped" and g["setting"] == first and "lymphatic" in g["reason"]
    for outcome, setting, msg in (("stopped", later, "cannot give the setting"),
                                  ("guard_failed", later, "cannot give the setting"),
                                  ("kept", later, "cannot give the setting"),
                                  ("changed", first, "cannot give the setting"),
                                  ("changed", {"min_score": 0.2, "min_margin": 0.25}, "cannot give the setting"),
                                  ("rerun", first, "is not one of")):
        n55(outcome, setting)
        with pytest.raises(fr.Refused, match=msg):
            fr.preflight(env.path, ruling, env.dst)


def test_a_killed_relabelling_is_recorded_afterwards_as_route_none(gate_env, monkeypatch, tmp_path, capsys):
    """The review of the port, lesser point 3: a run killed outright left only STARTED.json; --record-kill writes the
    stop record from it (reason killed, route none), and refuses without a marker or with a manifest or RECORD."""
    monkeypatch.setattr(fr.fio, "REPO", tmp_path)
    monkeypatch.setattr(fr.fio, "provenance", lambda code: {"git_commit": "abc", "git_uncommitted": []})
    (tmp_path / "data").mkdir()
    _output(gate_env, "data/derived/relabel_c")
    with pytest.raises(SystemExit, match="no started relabelling"):
        fr.main(["--record-kill"])
    dst = tmp_path / "data/derived/relabel_c"
    dst.mkdir(parents=True)
    (dst / fr.STARTED).write_text(json.dumps({"gate": GATE, "config_sha256": "other"}))
    with pytest.raises(SystemExit, match="differs from the one the killed run started with"):
        fr.main(["--record-kill"])
    (dst / fr.STARTED).write_text(json.dumps({"gate": GATE, "config_sha256": fr.file_sha256(fr.CONFIG)}))
    assert fr.main(["--record-kill"]) == fr.STOP_EXIT
    man = json.loads(fr.RECORD.read_text())
    assert man["stopped"] is True and man["route"] == "none" and man["reason"].startswith("killed")
    assert man["gate"] == GATE and man["outputs"] == {} and "route none" in capsys.readouterr().out
    with pytest.raises(SystemExit, match="already has a manifest or a RECORD"):
        fr.main(["--record-kill"])
    with pytest.raises(SystemExit):  # --ruling and --record-kill are one or the other
        fr.main(["--ruling", "2026-11-20 (a)", "--record-kill"])


def test_later_ruling_needs_an_entry_after_the_commit_other_than_the_runs(tmp_path, monkeypatch):
    plan = tmp_path / "PLAN.md"
    plan.write_text("- 2026-11-20 (a), the run's ruling; N54.\n- 2026-11-21 (a), old; N54.\n"
                    "- 2026-11-22 (a), the captain's later ruling; N54.\n- 2026-11-22 (b), something else.\n")
    monkeypatch.setattr(fr, "PLAN", plan)
    monkeypatch.setattr(fr, "last_commit", lambda path: ("c1", "2026-11-21"))
    monkeypatch.setattr(fr, "plan_at_commit", lambda c: "- 2026-11-20 (a), the run's ruling; N54.\n"
                                                        "- 2026-11-21 (a), old; N54.\n")
    after = fr.RECORD
    assert fr.later_ruling("2026-11-22 (a)", after, excluded=["2026-11-20 (a)"]) is None
    assert "no decision-log entry" in fr.later_ruling("2026-11-22 (b)", after)
    assert "authorised the run" in fr.later_ruling("2026-11-20 (a)", after, excluded=["2026-11-20 (a)"])
    assert "dated before" in fr.later_ruling("2026-11-20 (a)", after)
    assert "already in docs/PLAN.md" in fr.later_ruling("2026-11-21 (a)", after)
    monkeypatch.setattr(fr, "last_commit", lambda path: None)
    assert "is not committed" in fr.later_ruling("2026-11-22 (a)", after)


def test_later_ruling_refuses_when_the_plan_cannot_be_read_at_the_commit(tmp_path, monkeypatch):
    """The final check of the fix round, item E: a docs/PLAN.md that cannot be read at the commit of `after` (a
    shallow or foreign checkout) shows no ruling to be later, so the check refuses rather than passes."""
    plan = tmp_path / "PLAN.md"
    plan.write_text("- 2026-11-22 (a), the captain sets the N54 route aside.\n")
    monkeypatch.setattr(fr, "PLAN", plan)
    monkeypatch.setattr(fr, "last_commit", lambda path: ("c1", "2026-11-21"))
    monkeypatch.setattr(fr, "plan_at_commit", lambda c: None)
    assert "cannot be read at the commit" in fr.later_ruling("2026-11-22 (a)", fr.RECORD)
    monkeypatch.setattr(fr, "plan_at_commit", lambda c: "")
    assert fr.later_ruling("2026-11-22 (a)", fr.RECORD) is None


def test_later_ruling_needs_the_fixed_phrase_when_one_is_given(tmp_path, monkeypatch):
    """The final check of the fix round, item F: an entry that merely names N54 (here one that records the listing)
    does not pass where a fixed phrase is asked; the phrase may run over a line break and use a typographic
    apostrophe."""
    plan = tmp_path / "PLAN.md"
    plan.write_text("- 2026-11-22 (a), the listing of row N54 is recorded.\n"
                    "- 2026-11-22 (b), the captain sets the N54\n  route aside for T1.\n"
                    "- 2026-11-22 (c), the captain rules on T1\u2019s last upload: the none file; N54.\n"
                    "- 2026-11-22 (d), the captain sets the route of N54 aside.\n")
    monkeypatch.setattr(fr, "PLAN", plan)
    monkeypatch.setattr(fr, "last_commit", lambda path: ("c1", "2026-11-21"))
    monkeypatch.setattr(fr, "plan_at_commit", lambda c: "")
    none_, upload = fr.NONE_RULING_PHRASE, fr.UPLOAD_RULING_PHRASE
    assert (none_, upload) == ("sets the N54 route aside", "rules on T1's last upload")
    assert fr.later_ruling("2026-11-22 (a)", fr.RECORD) is None  # no phrase asked: as before
    for ruling, phrase in (("2026-11-22 (a)", none_), ("2026-11-22 (a)", upload), ("2026-11-22 (d)", none_),
                           ("2026-11-22 (b)", upload), ("2026-11-22 (c)", none_)):
        assert "does not contain the fixed phrase" in fr.later_ruling(ruling, fr.RECORD, phrase=phrase)
    assert fr.later_ruling("2026-11-22 (b)", fr.RECORD, phrase=none_) is None
    assert fr.later_ruling("2026-11-22 (c)", fr.RECORD, phrase=upload) is None


def test_last_commit_and_plan_at_commit_read_a_real_history(tmp_path, monkeypatch):
    import subprocess

    from conftest import git_env

    def g(*args):  # conftest's environment: never GIT_DIR or its kin (docs/PLAN.md decision log 2026-11-19 (t))
        subprocess.run(["git", "-C", str(tmp_path), "-c", "user.name=t", "-c", "user.email=t@example.org", "-c",
                        "commit.gpgsign=false", *args], cwd=tmp_path, check=True, capture_output=True, env=git_env())

    g("init", "-q")
    (tmp_path / "docs").mkdir()
    (tmp_path / "docs/PLAN.md").write_text("- 2026-11-20 (a), first; N54.\n")
    (tmp_path / "rec.json").write_text("{}")
    g("add", "-A")
    g("commit", "-q", "-m", "one")
    monkeypatch.setattr(fr.fio, "REPO", tmp_path)
    monkeypatch.setattr(fr, "PLAN", tmp_path / "docs/PLAN.md")
    commit, day = fr.last_commit(tmp_path / "rec.json")
    assert len(commit) == 40 and len(day) == 10
    assert fr.ruling_logged("2026-11-20 (a)", "N54", text=fr.plan_at_commit(commit))
    assert fr.last_commit(tmp_path / "absent.json") is None and fr.plan_at_commit("0" * 40) is None


def test_none_after_a_relabelling_refuses_the_runs_own_ruling_and_an_earlier_entry(
        bt, t1_env, monkeypatch, tmp_path):  # noqa: F811
    """The review of the port, lesser point 2: the ruling that authorised the run always exists after a relabelling,
    so it cannot set the route aside; nor can an entry already in docs/PLAN.md when the record was committed."""
    _none_env(bt, t1_env, monkeypatch, {"route": "recomputed", "gate": {"ruling": "2026-11-20 (a)"}}, tmp_path)
    out = tmp_path / "T1_none.h5ad"
    with pytest.raises(SystemExit, match="authorised the run"):
        bt.main(["t1", *TARGET, "--feng-groups", "none", "--feng-none-ruling", "2026-11-20 (a)", "--out", str(out)])
    monkeypatch.setattr(fr, "plan_at_commit", lambda commit: fr.PLAN.read_text())  # both entries were there already
    with pytest.raises(SystemExit, match="already in docs/PLAN.md"):
        bt.main(["t1", *TARGET, "--feng-groups", "none", "--feng-none-ruling", "2026-11-22 (a)", "--out", str(out)])
    monkeypatch.setattr(fr, "last_commit", lambda path: ("rc1", "2026-11-23"))
    with pytest.raises(SystemExit, match="dated before"):
        bt.main(["t1", *TARGET, "--feng-groups", "none", "--feng-none-ruling", "2026-11-22 (a)", "--out", str(out)])
    assert not out.exists()


def test_the_agreement_waits_for_the_captains_ruling_on_the_upload(tmp_path, monkeypatch):
    """The review of the port, lesser point 1: the agreement table runs only after a ruling on T1's last upload that
    came after the listing's commit, and never on the run's own ruling."""
    import n54_authors_agreement as ag

    listing = tmp_path / "n54.json"
    listing.write_text('{"inputs": {"recomputed": {"sha256": "x"}}}\n')
    monkeypatch.setattr(ag, "OUT", tmp_path / "agree")
    monkeypatch.setattr(ag, "LISTING", listing)
    record = json.dumps({"route": "recomputed", "gate": {"ruling": "2026-11-20 (a)"}}).encode()
    monkeypatch.setattr(fr, "committed_at_head", lambda p: listing.read_bytes() if Path(p) == listing else record)
    plan = tmp_path / "PLAN.md"
    plan.write_text("- 2026-11-20 (a), the run's ruling; N54.\n- 2026-11-21 (a), the listing noted; N54.\n"
                    "- 2026-11-23 (a), the captain rules on T1's last upload; N54.\n"
                    "- 2026-11-23 (b), the listing's numbers are noted for N54.\n")
    monkeypatch.setattr(fr, "PLAN", plan)
    monkeypatch.setattr(fr, "last_commit", lambda path: ("lc1", "2026-11-21"))
    monkeypatch.setattr(fr, "plan_at_commit", lambda c: "- 2026-11-20 (a), the run's ruling; N54.\n"
                                                        "- 2026-11-21 (a), the listing noted; N54.\n")

    def reached(*a, **k):
        raise RuntimeError("reached the reads")

    monkeypatch.setattr(feng_dirs, "load_relabel_config", reached)
    for ruling, msg in (("2026-11-20 (a)", "authorised the run"), ("2026-11-21 (a)", "already in docs/PLAN.md"),
                        ("2026-11-24 (a)", "no decision-log entry"), ("2026-11-23 (b)", "fixed phrase")):
        with pytest.raises(SystemExit, match=msg):
            ag.main(["--frozen-sha256", "x", "--method-report", "--upload-ruling", ruling])
    with pytest.raises(RuntimeError, match="reached the reads"):
        ag.main(["--frozen-sha256", "x", "--method-report", "--upload-ruling", "2026-11-23 (a)"])


@pytest.mark.parametrize("spec, route", [
    ({"cm": 30, "endo": 4, "vasc": 15, "ambig": 6, "epi": 4, "blood": 8, "none": 4}, "recomputed"),
    ({"cm": 30, "endo": 15, "vasc": 15, "ambig": 6, "epi": 14, "blood": 8, "none": 4}, "recomputed"),
    ({"cm": 4, "endo": 15, "vasc": 15, "ambig": 6, "epi": 14, "blood": 8, "none": 4}, "none")])
def test_the_labellers_route_and_the_builders_route_agree(tmp_path, monkeypatch, record_ok, spec, route):
    """The review of the port, lesser point 12: on one synthetic relabelling, the route the labeller writes in its
    manifest (group_cells and fallbacks on the labels) equals the route feng_dirs.recomputed_directions recomputes from
    the statistics file (directions_from_stats and apply_fallback), and the groups with a direction are those whose
    fallback takes something."""
    monkeypatch.setattr(fr.fio, "provenance", lambda code: {"git_commit": "abc", "git_uncommitted": []})
    types = [t for t in cells(spec) for _ in range(2)]
    cfg = syn_cfg(status="frozen")
    cfg_path = _cfg_file(tmp_path, cfg)
    src = _src(tmp_path, ["x"] * len(types), types)
    dst = tmp_path / "relabel_c"
    man = fr.run(cfg, cfg_path, src, dst, GATE)
    assert man["route"] == route and man["usable"] is True
    _, rec = feng_dirs.recomputed_directions(GENES, stats_path=dst / cfg["output"]["stats"], relabel_config=cfg_path)
    assert rec["route"] == man["route"]
    assert rec["fallbacks"] == man["fallbacks"]
    assert rec["with_direction"] == sorted(g for g, f in man["fallbacks"].items() if f["takes"] is not None)
    for g, c in man["group_cells"].items():  # the same cells at both ends
        assert (g in rec["groups"]) == (min(c.values()) >= cfg["min_cells_ext"])
        if g in rec["groups"]:
            assert {s: rec["groups"][g][s] for s in ("E9.5", "E14.5")} == c


def test_the_step_script_runs_every_builder_test_and_two_separate_check_builds():
    """The review of the port, lesser points 3 to 5: the tests step lists every test file that imports build_test,
    feng_dirs, feng_relabel or the harness's late proxy; build2 writes its own file; a kill is recorded."""
    sh = (ROOT / "experiments/external/n54_step.sh").read_text()
    pat = re.compile(r"build_test|feng_dirs|feng_relabel|late_proxy|LP_STATS")
    want = {f"tests/{p.name}" for p in (ROOT / "tests").glob("test_*.py") if pat.search(p.read_text())}
    listed = set(re.findall(r"tests/test_[a-z0-9_]+\.py", sh.split("  tests)")[1].split(";;")[0]))
    assert want <= listed and {"tests/test_core.py", "tests/test_s2_rules.py", "tests/test_git_helpers.py"} <= listed
    assert {"tests/test_h26_heart_scope_forms.py", "tests/test_t1_bets.py"} <= want  # main's builder tests since T61
    assert all((ROOT / f).exists() for f in listed)  # no listed file is missing
    assert "experiments/external/test_filters.py" in sh
    assert "T1_s2rc_b.h5ad" in sh and "build2.sha" in sh and sh.count("--record-kill") >= 2
    assert '--upload-ruling "${2:?' in sh
