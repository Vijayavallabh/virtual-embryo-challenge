"""Try T62 (experiments/harness/TRIES.md): experiments/t1_m1/recomputed_read.py on synthetic stages and hand-made task
files only. No released file, statistics file, board file or score is read; the committed relabelling record and the
frozen settings file are read for their fields and hashes."""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pytest
from test_g_mom_verify import _small_inputs
from test_t1_scope import GENES, HELD
from test_t1_state_transport import _gmom_ctx

from vec import harness
from vec.boards import Board
from vec.core import counts_of
from vec.t1_scope import build_submission, classify, load_spec, reweight

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "experiments" / "t1_m1" / "recomputed_read.py"
BOARD = Board(key="pseudo:t62", task="T1", genes=list(GENES), needs_coords=False, min_cells=100, max_cells=600)
SHARES = {"cm": 0.6, "endocardium": 0.25, "blood": 0.15}
DIRECTED = ("CM", "Endocardium", "Endothelium")


@pytest.fixture
def rc(monkeypatch, tmp_path):
    spec = importlib.util.spec_from_file_location("recomputed_read_module", SCRIPT)
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


HELD_M1 = sorted({*HELD, GENES[3], GENES[7]})


# ---------------------------------------------------------------------------------------------------------------------
# The guards on the source and on the relabelling


def test_the_source_asks_for_the_recomputed_groups_only(rc, tmp_path):
    assert rc.source_guard() == []
    calls = rc._calls(SCRIPT.read_text())
    assert calls and all('groups="recomputed"' in c for c in calls)
    for text in ('x = fd.feng_directions(list(g))\n', 'x = fd.feng_directions(g, groups="auth' + 'ors")\n',
                 'x = fd.feng_directions(list(g), groups="none")\n', 'p = fd.clean_' + 'stats_path()\n',
                 'open("n54' + '.json")\n'):
        bad = tmp_path / "bad.py"
        bad.write_text(text)
        assert rc.source_guard(bad), text


def test_the_relabelling_checks_on_the_committed_record_and_settings(rc, tmp_path):
    ok = rc.relabel_checks(need_stats=False)
    assert ok["passes"] and ok["route"] == "recomputed" and ok["authors_labels_read"] is False
    assert ok["yaml_sha256"] == rc.YAML_SHA256 and ok["record_stats_sha256"] == rc.STATS_SHA256
    other = tmp_path / "stats.npz"
    other.write_bytes(b"not the file")
    r = rc.relabel_checks(stats_path=other)
    assert not r["passes"] and any("statistics file's sha256" in x for x in r["reasons"])
    r = rc.relabel_checks(stats_path=tmp_path / "missing.npz")
    assert not r["passes"]
    rec = json.loads(rc.RECORD.read_text())
    for change, msg in (({"route": "none"}, "route"), ({"authors_labels_read": True}, "authors_labels_read"),
                        ({"outputs": {}}, "statistics sha256")):
        bad = tmp_path / "MANIFEST.json"
        bad.write_text(json.dumps({**rec, **change}))
        r = rc.relabel_checks(record=bad, need_stats=False)
        assert not r["passes"] and any(msg in x for x in r["reasons"])
    yaml = tmp_path / "feng_relabel.yaml"
    yaml.write_text("status: draft\n")
    assert not rc.relabel_checks(yaml_path=yaml, need_stats=False)["passes"]


def test_recomputed_inputs_asks_feng_dirs_for_the_recomputed_groups_and_refuses_another_route(rc, monkeypatch):
    seen = []

    def fake(genes, **kw):
        seen.append(kw)
        return {"CM": np.zeros(len(genes))}, {"route": kw.get("groups")}

    monkeypatch.setattr(rc.feng_dirs, "feng_directions", fake)
    dirs, ext = rc.recomputed_inputs(["a", "b"])
    assert seen == [{"groups": "recomputed"}] and set(dirs) == {"CM"}
    monkeypatch.setattr(rc.feng_dirs, "feng_directions", lambda genes, **kw: ({}, {"route": "none"}))
    with pytest.raises(ValueError, match="routes to"):
        rc.recomputed_inputs(["a"])
    monkeypatch.setattr(rc.feng_dirs, "feng_directions",
                        lambda genes, **kw: ({"Blood": np.zeros(1)}, {"route": "recomputed"}))
    with pytest.raises(ValueError, match="Blood"):
        rc.recomputed_inputs(["a"])


def test_the_practice_fraction_is_one_fifth_and_the_count_is_ten(rc):
    assert rc.feng_dirs.bracket_fraction(9.5, 10.5) == pytest.approx(rc.PRACTICE_FRACTION)
    assert rc.feng_dirs.bracket_fraction(10.5, 12.5) == pytest.approx(2 * rc.PRACTICE_FRACTION)
    assert rc.TRIES_COUNT == 10 and rc.COUNTED == ("m1rc_gmom", "m1rc_gmom_epi_held")
    assert rc.ARMS == ("none_gmom", "m1rc_gmom", "m1rc_gmom_epi_held") and rc.LIB_ARMS[0] == "copy_last"
    assert rc.EPI_GROUP == "Epicardium" and rc.EPI_TASK == "3"


def test_the_run_is_refused_without_its_row(rc, tmp_path, monkeypatch):
    tries = tmp_path / "TRIES.md"
    tries.write_text("| T61 (one row) | T1 | x |\n")
    monkeypatch.setattr(rc, "TRIES", tries)
    monkeypatch.setattr(rc, "relabel_checks", lambda *a, **k: pytest.fail("nothing before the row's refusal"))
    with pytest.raises(SystemExit, match="no row"):
        rc.refuse(tmp_path / "out", "x.json")


# ---------------------------------------------------------------------------------------------------------------------
# The candidate's g_mom step: the builder's under --feng-groups recomputed


def test_the_candidates_hook_is_the_builders_under_recomputed(rc, monkeypatch):
    ref, start = _classified()
    factors = rc.br.gmom_fits(ref, "E9.5_RNA", start, "E8.5_RNA")["file"]
    monkeypatch.setattr(rc.bt, "load_stage", lambda *a, **k: start)
    monkeypatch.setattr(rc.bt.harness, "g_mom_pair_factors", lambda *a, **k: (factors, None))
    monkeypatch.setattr(rc.bt.t1_scope, "held_flat_genes", lambda *a, **k: list(HELD))
    for directed in (DIRECTED, ("CM",)):
        dirs = _directions(ref.n_vars, directed)
        groups, covered, with_direction = rc.direction_kit(ref.obs["celltype"].astype(str).unique(), dirs)
        theirs, _ = rc.bt._t1_g_mom(ref, sorted(factors.categories), 1000001, covered, with_direction=with_direction)
        pop, _ = reweight(ref, SHARES, 300, 1000001)
        a, b = rc.rc_hook(factors, 1000001, HELD, with_direction, covered)(pop), theirs(pop)
        assert (_dense(a.layers["counts"]) == _dense(b.layers["counts"])).all()
        edited = (_dense(a.layers["counts"]) != _dense(counts_of(pop))).any(axis=1)
        labels = pop.obs["celltype"].astype(str).to_numpy()
        assert edited.any() and set(pop.obs["category"][edited]) <= {"cm", "endocardium"}
        if directed == ("CM",):  # Endocardium has no direction here: its cells are held back
            assert not edited[labels == "Endocardium"].any()


def test_with_cm_endocardium_and_endothelium_directed_the_mask_is_none_gmoms(rc):
    """On every label the committed chamber config can place in cm or endocardium, the candidate's mask (groups with a
    recomputed direction) and none_gmom's (labels with a group) hold back the same labels, Early VIC alone."""
    import yaml

    cfg = yaml.safe_load((ROOT / "configs" / "t1_chamber_classes.yaml").read_text())
    edit = set(rc.bt.G_MOM_EDIT_CATEGORIES)
    labels = {str(lab) for lab, cat in cfg["labels"].items() if cat in edit}
    for per_file in cfg["file_overrides"].values():
        labels |= {str(lab) for lab, cat in per_file.items() if cat in edit}
    for per_file in cfg["file_fallbacks"].values():
        for cat, stand_ins in per_file.items():
            if cat in edit:
                labels |= set(map(str, stand_ins))
    for directed in (DIRECTED, (*DIRECTED, "Epicardium")):
        assert rc.bt._labels_without_direction(labels, directed) == rc.bt._labels_without_group(labels) == {
            "Early VIC"}
    assert "V-CM" in rc.bt._labels_without_direction(labels, ("Endocardium", "Endothelium"))


# ---------------------------------------------------------------------------------------------------------------------
# The arms


def test_the_cut_pair_none_gmom_is_t61s_and_the_candidate_is_the_builders_path(rc, monkeypatch):
    ref, start, _ = _small_inputs()
    ctx = _gmom_ctx(ref, start, seed=1000002)
    fits = rc.br.gmom_fits(ref, "E9.5_RNA", start, "E8.5_RNA")
    dirs = _directions(ref.n_vars)
    monkeypatch.setattr(rc.feng_dirs, "feng_directions", lambda *a, **k: pytest.fail("the directions are given"))
    ours = rc.cut_pair(ctx, fits["file"], dirs, 0.2, HELD_M1)
    t61, _, _ = rc.br.cut_arms(ctx, fits)
    assert tuple(ours) == rc.ARMS
    assert ours["none_gmom"].X.tobytes() == t61["none_gmom"].X.tobytes()
    assert ours["m1rc_gmom"].X.tobytes() != ours["none_gmom"].X.tobytes()
    # the candidate: build_submission with the recomputed direction, R2's genes and the builder's own hook
    monkeypatch.setattr(rc.bt, "load_stage", lambda *a, **k: start)
    monkeypatch.setattr(rc.bt.harness, "g_mom_pair_factors", lambda *a, **k: (fits["file"], None))
    monkeypatch.setattr(rc.bt.t1_scope, "held_flat_genes", lambda *a, **k: list(HELD))
    cref = harness._shallow(ref)
    classify(cref, source="E9.5_RNA")
    groups, covered, with_direction = rc.direction_kit(cref.obs["celltype"].astype(str).unique(), dirs)
    hook, _ = rc.bt._t1_g_mom(cref, sorted(fits["file"].categories), ctx.seed, covered, with_direction=with_direction)
    shares = harness.t1_chamber_shares(cref, "E9.5_RNA")
    sub, _ = build_submission(cref, ctx.board, shares, harness.count_for(ctx.board, cref.n_obs), ctx.seed,
                              held_flat=HELD_M1, ref_source="E9.5_RNA", directions=dirs, direction_fraction=0.2,
                              direction_groups=groups, post_cells=hook)
    want = harness.finish(ctx, cref, rc.dense(sub.X), None)
    assert ours["m1rc_gmom"].X.tobytes() == want.X.tobytes()
    # no label of this stage maps to Epicardium: the two candidate arms are the same file
    assert ours["m1rc_gmom_epi_held"].X.tobytes() == ours["m1rc_gmom"].X.tobytes()


def test_a_group_without_a_direction_keeps_its_drawn_counts_through_the_candidate(rc):
    """Endocardium without a recomputed direction (only CM directed): its cells are neither moved nor edited."""
    e95, e85 = _classified()
    spec = load_spec()
    fb = {"E8.5": spec.fallbacks("E8.5_RNA"), "E9.5": spec.fallbacks("E9.5_RNA")}
    within = {"surrounding": ["Foregut", "Epithelium"]}
    factors = rc.br.gmom_fits(e95, "E9.5_RNA", e85, "E8.5_RNA")["file"]
    arms = rc.library_pair(e95, SHARES, factors, HELD, HELD_M1, 4, within, fb, _directions(e95.n_vars, ("CM",)), 0.2,
                           n_cells=300)
    base, _ = reweight(e95, SHARES, 300, 4, within=within, fallback=fb["E9.5"])
    labels = base.obs["celltype"].astype(str).to_numpy()
    held = np.isin(labels, ["Endocardium", "Blood", "Foregut"])
    got = _dense(counts_of(arms["m1rc_gmom"]))
    assert np.array_equal(got[held], _dense(counts_of(base))[held])
    assert not np.array_equal(got[labels == "V-CM"], _dense(counts_of(base))[labels == "V-CM"])


def test_the_library_pair_none_gmom_is_t61s_and_the_route_is_the_build_paths(rc, monkeypatch):
    e95, e85 = _classified()
    spec = load_spec()
    fb = {"E8.5": spec.fallbacks("E8.5_RNA"), "E9.5": spec.fallbacks("E9.5_RNA")}
    within = {"surrounding": ["Foregut", "Epithelium"]}
    fits = rc.br.gmom_fits(e95, "E9.5_RNA", e85, "E8.5_RNA")
    dirs = _directions(e95.n_vars)
    monkeypatch.setattr(rc.feng_dirs, "feng_directions", lambda *a, **k: pytest.fail("the directions are given"))
    ours = rc.library_pair(e95, SHARES, fits["file"], HELD, HELD_M1, 4, within, fb, dirs, 0.2, n_cells=300)
    from vec import t1_bets

    axes, _ = t1_bets.momentum_axes(e85, e95, sorted(SHARES), HELD, within=within, start_fallback=fb["E8.5"],
                                    later_fallback=fb["E9.5"])
    t61, _, _ = rc.br.library_arms(e85, e95, SHARES, fits, HELD, 4, within, fb, t1_bets.advance_scores(e95, axes),
                                   n_cells=300)
    assert (_dense(ours["none_gmom"].X) == _dense(t61["none_gmom"].X)).all()
    assert (_dense(counts_of(ours["none_gmom"])) == _dense(counts_of(t61["none_gmom"]))).all()
    # moved_by_feng on the drawn base is build_submission's direction step on the same draw (stream [seed, 2])
    groups, covered, with_direction = rc.direction_kit(e95.obs["celltype"].astype(str).unique(), dirs)
    seen = []
    build_submission(e95, BOARD, SHARES, 300, 4, held_flat=HELD_M1, ref_source="E9.5_RNA", within=within,
                     fallback=fb["E9.5"], directions=dirs, direction_fraction=0.2, direction_groups=groups,
                     post_cells=lambda pop: seen.append(pop) or pop)
    hooked = rc.rc_hook(fits["file"], 4, HELD, with_direction, covered)(seen[0])
    assert (_dense(counts_of(hooked)) == _dense(counts_of(ours["m1rc_gmom"]))).all()
    assert not (_dense(ours["m1rc_gmom"].X) == _dense(ours["none_gmom"].X)).all()


# ---------------------------------------------------------------------------------------------------------------------
# The second counted arm: the Epicardium group's step taken out

EPI_SHARES = {"cm": 0.5, "endocardium": 0.2, "epicardium": 0.15, "blood": 0.05, "surrounding": 0.1}
WITH_EPI = (*DIRECTED, "Epicardium")


def _epi_world():
    """An E9.5-like stage with Proepicardium cells (the readout config maps them to Epicardium) and an E8.5-like
    start, classified; g_mom's file factors on them."""
    from test_t1_state_transport import classified_pair

    start, ref = classified_pair(seed=5, n_scale=0.5)
    return ref, start


def test_the_epi_held_arm_differs_only_in_epicardium_cells_on_the_library_route(rc):
    e95, e85 = _epi_world()
    spec = load_spec()
    fb = {"E8.5": spec.fallbacks("E8.5_RNA"), "E9.5": spec.fallbacks("E9.5_RNA")}
    within = {"surrounding": ["Foregut"]}
    factors = rc.br.gmom_fits(e95, "E9.5_RNA", e85, "E8.5_RNA")["file"]
    dirs = _directions(e95.n_vars, WITH_EPI)
    arms = rc.library_pair(e95, EPI_SHARES, factors, HELD, HELD_M1, 6, within, fb, dirs, 0.2, n_cells=300)
    base, _ = reweight(e95, EPI_SHARES, 300, 6, within=within, fallback=fb["E9.5"])
    labels = base.obs["celltype"].astype(str).to_numpy()
    groups = rc.feng_dirs.label_groups(sorted(set(labels)))
    epi = np.array([groups[lab] == "Epicardium" for lab in labels])
    assert epi.any() and (~epi).any() and set(labels[epi]) == {"Proepicardium"}
    full, held = (_dense(counts_of(arms[a])) for a in ("m1rc_gmom", "m1rc_gmom_epi_held"))
    drawn = _dense(counts_of(base))
    differs = (full != held).any(axis=1)
    assert differs.any() and not differs[~epi].any()  # only cells whose label maps to Epicardium differ
    assert np.array_equal(held[epi], drawn[epi])  # they keep their drawn counts: not moved, no momentum, no g_mom
    assert not np.array_equal(full[epi], drawn[epi])  # the candidate moved them
    # none_gmom leaves them as drawn too (g_mom edits cardiomyocytes and endocardium only)
    assert np.array_equal(_dense(counts_of(arms["none_gmom"]))[epi], drawn[epi])
    # the same arm with no Epicardium direction at all is the candidate itself
    same = rc.library_pair(e95, EPI_SHARES, factors, HELD, HELD_M1, 6, within, fb, _directions(e95.n_vars), 0.2,
                           n_cells=300)
    assert np.array_equal(_dense(counts_of(same["m1rc_gmom"])), _dense(counts_of(same["m1rc_gmom_epi_held"])))


def test_the_epi_held_arm_differs_only_in_epicardium_cells_on_the_cut_route(rc, monkeypatch):
    ref, start = _epi_world()
    ctx = _gmom_ctx(ref, start, seed=1000003)
    fits = rc.br.gmom_fits(ref, "E9.5_RNA", start, "E8.5_RNA")
    seen = []
    real = rc.hold_group

    def spy(moved, drawn, groups, group="Epicardium"):
        seen.append(np.array([groups.get(x) == group for x in drawn.obs["celltype"].astype(str)]))
        return real(moved, drawn, groups, group)

    monkeypatch.setattr(rc, "hold_group", spy)
    ours = rc.cut_pair(ctx, fits["file"], _directions(ref.n_vars, WITH_EPI), 0.2, HELD_M1)
    assert tuple(ours) == rc.ARMS and len(seen) == 1
    epi = seen[0]
    full, held = ours["m1rc_gmom"].X, ours["m1rc_gmom_epi_held"].X
    differs = (np.asarray(full) != np.asarray(held)).any(axis=1)
    assert epi.any() and differs.any() and not differs[~epi].any()
    assert (np.asarray(held)[epi] == np.asarray(ours["none_gmom"].X)[epi]).all()  # as drawn, as in none_gmom


def test_taking_the_direction_out_before_the_step_would_redraw_the_other_groups(rc):
    """Why the arm holds the Epicardium cells after the step and not before it: vec.core.apply_factors draws one stream
    over every moved cell, so a direction step without the Epicardium group gives other groups' cells new draws."""
    e95, _ = _epi_world()
    spec = load_spec()
    base, _ = reweight(e95, EPI_SHARES, 300, 6, within={"surrounding": ["Foregut"]},
                       fallback=spec.fallbacks("E9.5_RNA"))
    labels = base.obs["celltype"].astype(str).to_numpy()
    groups = rc.feng_dirs.label_groups(sorted(set(labels)))
    dirs = _directions(e95.n_vars, WITH_EPI)
    a = rc.gc.moved_by_feng(base, dirs, 0.2, groups, HELD_M1, [6, 2])
    b = rc.gc.moved_by_feng(base, rc.without_epicardium(dirs), 0.2, groups, HELD_M1, [6, 2])
    epi = np.array([groups[lab] == "Epicardium" for lab in labels])
    assert ((_dense(counts_of(a)) != _dense(counts_of(b))).any(axis=1) & ~epi).any()
    held = rc.hold_group(a, base, groups)
    assert np.array_equal(_dense(counts_of(held))[~epi], _dense(counts_of(a))[~epi])
    assert np.array_equal(_dense(counts_of(held))[epi], _dense(counts_of(base))[epi])
    assert rc.without_epicardium(dirs).keys() == set(DIRECTED)
    with pytest.raises(ValueError, match="same draw"):
        rc.hold_group(a, base[::-1].copy(), groups)


# ---------------------------------------------------------------------------------------------------------------------
# The reliability listing (target-free)


def _relabel_stats(genes, seed=2, small=None):
    """Synthetic recomputed statistics in the file's key layout ('stage|zone|label|half'), per group label a noisy
    profile with a stage trend; `small` gives that label too few cells at E14.5 for a direction of its own."""
    rng = np.random.default_rng(seed)
    base = rng.gamma(1.0, 5.0, len(genes))
    keys, n, cp, cells = [], [], [], {}
    for label in ("r_cm", "r_endocardial_ec", "r_vascular_ec", "r_ec_ambiguous", "r_epicardial"):
        for i, stage in enumerate(("E9.5", "E14.5")):
            trend = np.exp((-1.0) ** (np.arange(len(genes)) % 2) * 0.3 * i)
            for half in "AB":
                k = 2 if (label, stage) == small else int(rng.integers(20, 60))
                keys.append(f"{stage}|LV|{label}|{half}")
                n.append(k)
                cp.append(k * base * trend * rng.lognormal(0.0, 0.3, len(genes)))
                cells[(label, stage, half)] = k
    st = {"genes": np.array(genes, dtype=object), "keys": np.array(keys, dtype=object), "n": np.array(n, float),
          "S_cp": np.array(cp), "S_log": np.zeros((len(keys), len(genes))), "S_nz": np.zeros((len(keys), len(genes))),
          "depth_sum": np.array(n, float) * 1000.0}
    return st, cells


def test_the_reliability_listing_reads_each_groups_own_halves(rc):
    from scipy.stats import spearmanr

    genes = [f"G{i}" for i in range(120)]
    rcfg = rc.feng_dirs.load_relabel_config()
    cfg = rc.feng_dirs.groups_config(rc.feng_dirs.load_config(), rcfg)
    st, _ = _relabel_stats(genes, small=("r_epicardial", "E14.5"))
    own, _, _ = rc.feng_dirs.directions_from_stats(st, genes, cfg)
    dirs, _ = rc.feng_dirs.apply_fallback(own, rcfg)
    assert "Epicardium" not in rc.reliability_listing(st, genes, cfg, dirs)  # no direction (no fallback): not listed
    st2, cells2 = _relabel_stats(genes)
    own2, _, _ = rc.feng_dirs.directions_from_stats(st2, genes[:100] + ["absent"], cfg)
    out = rc.reliability_listing(st2, genes[:100] + ["absent"], cfg, own2)
    assert set(out) == {"CM", "Endocardium", "Endothelium", "Epicardium"}
    epi = out["Epicardium"]
    assert epi["own_direction"] and epi["labels"] == ["r_epicardial"]
    assert epi["cells"] == {s: {h: cells2[("r_epicardial", s, h)] for h in "AB"} for s in ("E9.5", "E14.5")}
    ec = out["Endothelium"]["cells"]["E9.5"]["A"]
    assert ec == sum(cells2[(lab, "E9.5", "A")] for lab in ("r_endocardial_ec", "r_vascular_ec", "r_ec_ambiguous"))
    h = rc.feng_dirs.direction_halves(st2, genes[:100] + ["absent"], cfg)["Epicardium"]
    assert epi["mapped_genes"] == 100  # the unmapped gene is left out
    ok = np.isfinite(h["d_a"])
    assert epi["spearman_all"] == pytest.approx(spearmanr(h["d_a"][ok], h["d_b"][ok]).statistic)
    _, factors = rc.feng_dirs.eb_shrink(own2["Epicardium"], h["d_a"], h["d_b"], h["mean_cp"], 10)
    assert epi["eb_factors"] == pytest.approx(factors) and epi["eb_factor_mean"] == pytest.approx(np.mean(factors))
    assert len(epi["deciles"]) == 10 and sum(x["genes"] for x in epi["deciles"]) == 100
    low = np.argsort(h["mean_cp"][ok], kind="stable")[:10]
    assert epi["deciles"][0]["spearman"] == pytest.approx(spearmanr(h["d_a"][ok][low], h["d_b"][ok][low]).statistic)
    # a group whose direction came from the yaml's fallback has no halves of its own
    st3, _ = _relabel_stats(genes, small=("r_endocardial_ec", "E14.5"))
    own3, _, _ = rc.feng_dirs.directions_from_stats(st3, genes, cfg)
    dirs3, fired3 = rc.feng_dirs.apply_fallback(own3, rcfg)
    out3 = rc.reliability_listing(st3, genes, cfg, dirs3, fired3)
    assert out3["Endocardium"] == {"own_direction": False, "fallback": fired3["Endocardium"]}
    assert fired3["Endocardium"]["takes"] == "Endothelium"
    text = rc.reliability_text({"label": rc.RELIABILITY_LABEL, "bracket": ["E9.5", "E14.5"], "genes_from": "T1:val",
                                "of_genes": len(genes), "with_direction": sorted(dirs3), "groups": out3,
                                "stamp": {"commit": "c" * 40, "tree_clean": True, "utc": "x"}})
    assert "Endocardium: no direction of its own" in text and "decile 10" in text


def test_the_reliability_step_is_refused_without_its_row_and_reads_no_t1_cell(rc, tmp_path, monkeypatch):
    tries = tmp_path / "TRIES.md"
    tries.write_text("| T61 (one row) | T1 | x |\n")
    monkeypatch.setattr(rc, "TRIES", tries)
    monkeypatch.setattr(rc, "relabel_checks", lambda *a, **k: pytest.fail("nothing before the row's refusal"))
    with pytest.raises(SystemExit, match="no row"):
        rc.main(["reliability", "--out", "x"])
    import inspect

    body = inspect.getsource(rc.cmd_reliability)
    for banned in ("load_stage", "read_stage", "load_tasks", "score_arms", "filter_target", "a.data", "DATA_DIR",
                   "committed_t61"):
        assert banned not in body, banned


# ---------------------------------------------------------------------------------------------------------------------
# The reading


def _arm(score):
    return {"score": score, "groups": {g: score / 100 for g in ("de_recovery", "direction", "distribution",
                                                               "covariance")}}


def _relabel(rc, **change):
    return {"sha256": rc.STATS_SHA256, "stats_sha256_hashed": rc.STATS_SHA256, "yaml_sha256": rc.YAML_SHA256,
            "route": "recomputed", "record_route": "recomputed", "authors_labels_read": False,
            "epicardium_directed": True, **change}


def _task(rc, gain, lib=False, per_seed=None, control=50.0, commit="c" * 40, clean=True, seeds=range(10),
          inverted=False, relabel=None, epi=0.0, epi_per_seed=None):
    """`epi`: m1rc_gmom - m1rc_gmom_epi_held per seed (the epicardium step's gain)."""
    rows = []
    for s in seeds:
        g = per_seed[s] if per_seed is not None else gain
        e = epi_per_seed[s] if epi_per_seed is not None else epi
        arms = {"none_gmom": _arm(control + s / 100), "m1rc_gmom": _arm(control + s / 100 + g),
                "m1rc_gmom_epi_held": _arm(control + s / 100 + g - e)}
        row = {"seed": s, "arms": arms}
        if lib:
            arms["copy_last"] = _arm(40.0)
            row["inverted_anchors"] = ["x"] if inverted and s == 3 else []
        else:
            row["floor"] = 50.0
        rows.append(row)
    return {"rows": rows, "stamp": {"commit": commit, "tree_clean": clean},
            "relabel": relabel if relabel is not None else _relabel(rc)}


def _committed(control=50.0):
    return {k: {s: control + s / 100 for s in range(10)} for k in ("chamber", "whole", "1", "2", "3")}


def _loaded(rc, chamber=1.0, whole=0.0, libs=(0.0, 0.0, 0.0), **kw):
    out = {"chamber": _task(rc, chamber), "whole": _task(rc, whole)}
    out.update({lib: _task(rc, g, lib=True) for lib, g in zip(("1", "2", "3"), libs, strict=True)})
    out.update(kw)
    return out


OK = {"reasons": []}


def test_ahead_on_both_readings(rc):
    r = rc.t62_reading(_loaded(rc, 1.0, -0.1, (0.2, -0.29, 1.0)), _committed(), OK)
    assert (r["reading_S"], r["reading_D"]) == ("AHEAD", "AHEAD")
    assert r["advice_fixed_before_the_run"] == "the recomputed file"
    assert r["j6_range"] == pytest.approx([0.0, 1.0]) and r["tasks"]["chamber"]["ahead"] == 10
    assert r["tasks"]["chamber"]["mean_scores"]["m1rc_gmom"] == pytest.approx(51.045)
    assert r["tasks"]["chamber"]["skill_differences"]["direction"] == pytest.approx(0.01)


@pytest.mark.parametrize("chamber,whole,per_seed", [(0.3, 0.0, None), (1.0, -0.31, None),
                                                    (None, 0.0, [1.0] * 7 + [-0.1] * 3)])
def test_not_ahead(rc, chamber, whole, per_seed):
    loaded = _loaded(rc, 0.0, whole)
    loaded["chamber"] = _task(rc, chamber, per_seed=per_seed)
    r = rc.t62_reading(loaded, _committed(), OK)
    assert (r["reading_S"], r["reading_D"]) == ("NOT AHEAD", "NOT AHEAD")
    assert r["advice_fixed_before_the_run"] == "the ruled file (row R8)"


def test_the_whole_cut_at_its_floor_is_ahead_and_a_library_below_its_floor_stops_d(rc):
    r = rc.t62_reading(_loaded(rc, 1.0, -0.3, (0.0, -0.31, 0.0)), _committed(), OK)
    assert (r["reading_S"], r["reading_D"]) == ("AHEAD", "NOT AHEAD") and r["libraries_below_floor"] == ["2"]
    assert r["advice_fixed_before_the_run"] == "the ruled file (row R8)"


def test_too_few_valid_libraries_stop_d(rc):
    loaded = _loaded(rc, 1.0, 0.0)
    for lib in ("1", "2"):
        loaded[lib] = _task(rc, 0.0, lib=True, inverted=True)
    r = rc.t62_reading(loaded, _committed(), OK)
    assert r["valid_libraries"] == ["3"] and (r["reading_S"], r["reading_D"]) == ("AHEAD", "NOT AHEAD")
    loaded["2"] = _task(rc, 0.0, lib=True)
    assert rc.t62_reading(loaded, _committed(), OK)["reading_D"] == "AHEAD"


def _with_lib3_epi(rc, loaded, epi=0.0, epi_per_seed=None, inverted=False):
    loaded["3"] = _task(rc, 0.0, lib=True, epi=epi, epi_per_seed=epi_per_seed, inverted=inverted)
    return loaded


@pytest.mark.parametrize("epi,per_seed,word", [(0.5, None, "PAYS"), (-0.5, None, "COSTS"), (0.29, None, "within noise"),
                                               (-0.29, None, "within noise"), (None, [1.0] * 7 + [-0.1] * 3,
                                                                              "within noise"),
                                               (None, [-1.0] * 8 + [0.1] * 2, "COSTS"),
                                               (None, [1.0] * 8 + [-0.1] * 2, "PAYS")])
def test_the_epicardium_statement_on_library_3(rc, epi, per_seed, word):
    loaded = _with_lib3_epi(rc, _loaded(rc, 1.0, 0.0), epi=epi or 0.0, epi_per_seed=per_seed)
    r = rc.t62_reading(loaded, _committed(), OK)
    assert r["epicardium_step"] == word and r["epicardium_task_valid"]
    assert set(r["epicardium"]) == {"chamber", "whole", "1", "2", "3"}  # reported on every task


def test_the_epicardium_statement_reads_library_3_alone(rc):
    loaded = _loaded(rc, 1.0, 0.0)
    loaded["chamber"] = _task(rc, 1.0, epi=-2.0)
    loaded["1"] = _task(rc, 0.0, lib=True, epi=-2.0)
    r = rc.t62_reading(loaded, _committed(), OK)
    assert r["epicardium_step"] == "within noise" and r["epicardium"]["chamber"]["mean"] == pytest.approx(-2.0)
    assert r["epicardium"]["1"]["behind"] == 10


def test_a_third_file_is_named_only_when_the_recomputed_file_is_advised_and_the_step_costs(rc):
    costs = rc.t62_reading(_with_lib3_epi(rc, _loaded(rc, 1.0, 0.0), epi=-0.5), _committed(), OK)
    assert costs["advice_fixed_before_the_run"] == "the recomputed file" and costs["third_file_named"]
    assert "third file (the recomputed groups with the epicardium group held): yes" in rc.outcome_lines(costs)[2]
    pays = rc.t62_reading(_with_lib3_epi(rc, _loaded(rc, 1.0, 0.0), epi=0.5), _committed(), OK)
    assert pays["advice_fixed_before_the_run"] == "the recomputed file" and not pays["third_file_named"]
    noise = rc.t62_reading(_with_lib3_epi(rc, _loaded(rc, 1.0, 0.0), epi=0.1), _committed(), OK)
    assert not noise["third_file_named"]
    ruled = rc.t62_reading(_with_lib3_epi(rc, _loaded(rc, 0.1, 0.0), epi=-0.5), _committed(), OK)
    assert ruled["epicardium_step"] == "COSTS" and ruled["advice_fixed_before_the_run"] == "the ruled file (row R8)"
    assert not ruled["third_file_named"]


def test_an_invalid_library_3_is_reported_beside_the_epicardium_statement(rc):
    r = rc.t62_reading(_with_lib3_epi(rc, _loaded(rc, 1.0, 0.0), epi=-0.5, inverted=True), _committed(), OK)
    assert r["epicardium_step"] == "COSTS" and not r["epicardium_task_valid"] and r["valid_libraries"] == ["1", "2"]
    assert "a valid library: no" in rc.outcome_lines(r)[1]


def test_a_task_without_an_epicardium_direction_is_named(rc):
    loaded = _loaded(rc, 1.0, 0.0)
    loaded["2"] = _task(rc, 0.0, lib=True, relabel=_relabel(rc, epicardium_directed=False))
    r = rc.t62_reading(loaded, _committed(), OK)
    assert r["epicardium_not_directed"] == ["2"]
    assert "no recomputed direction on: 2" in rc.outcome_lines(r)[-1]


@pytest.mark.parametrize("case", ["control_off", "missing_seed", "missing_arm", "missing_epi_arm", "missing_task",
                                  "two_commits", "unclean", "stats", "yaml", "route", "authors", "checks", "source"])
def test_every_guard_voids(rc, case):
    loaded, committed, checks, source = _loaded(rc), _committed(), OK, []
    if case == "control_off":
        committed["2"][4] += 2e-9
    elif case == "missing_seed":
        loaded["whole"] = _task(rc, 0.0, seeds=range(9))
    elif case == "missing_arm":
        del loaded["3"]["rows"][0]["arms"]["copy_last"]
    elif case == "missing_epi_arm":
        del loaded["chamber"]["rows"][7]["arms"]["m1rc_gmom_epi_held"]
    elif case == "missing_task":
        loaded["1"] = None
    elif case == "two_commits":
        loaded["1"] = _task(rc, 0.0, lib=True, commit="d" * 40)
    elif case == "unclean":
        loaded["chamber"] = _task(rc, 1.0, clean=False)
    elif case == "stats":
        loaded["chamber"] = _task(rc, 1.0, relabel=_relabel(rc, stats_sha256_hashed="0" * 64))
    elif case == "yaml":
        loaded["whole"] = _task(rc, 0.0, relabel=_relabel(rc, yaml_sha256="0" * 64))
    elif case == "route":
        loaded["1"] = _task(rc, 0.0, lib=True, relabel=_relabel(rc, route="none"))
    elif case == "authors":
        loaded["2"] = _task(rc, 0.0, lib=True, relabel=_relabel(rc, authors_labels_read=True))
    elif case == "checks":
        checks = {"reasons": ["the settings file's sha256 is x"]}
    elif case == "source":
        source = ["a feng_directions call without the recomputed groups"]
    r = rc.t62_reading(loaded, committed, checks, source)
    assert r["outcome"] == "VOID" and r["reasons"]
    assert rc.outcome_lines(r)[0].startswith("Outcome: VOID")


def test_none_gmom_within_the_tolerance_passes(rc):
    committed = _committed()
    committed["2"][4] += 5e-10
    assert rc.t62_reading(_loaded(rc), committed, OK)["outcome"] != "VOID"


def test_committed_t61_reads_none_gmom_per_seed(rc, tmp_path):
    (tmp_path / "lib1.json").write_text(json.dumps({"rows": [{"seed": s, "arms": {"none_gmom": {"score": 60.0 + s}}}
                                                            for s in range(10)]}))
    assert rc.committed_t61("lib1.json", tmp_path) == {s: 60.0 + s for s in range(10)}
    assert rc.task_files() == {"chamber": "T1_E9.5_E10.5_chamber.json", "whole": "T1_E9.5_E10.5.json",
                               "1": "lib1.json", "2": "lib2.json", "3": "lib3.json"}


def test_the_read_prints_no_score(rc, tmp_path, monkeypatch, capsys):
    out = tmp_path / "results" / "t62"
    (out / "tasks").mkdir(parents=True)
    loaded = _loaded(rc, 1.234567, 0.5, (0.4, 0.4, 0.4))
    for k, f in rc.task_files().items():
        (out / "tasks" / f).write_text(json.dumps(loaded[k]))
    committed = _committed()
    monkeypatch.setattr(rc, "committed_t61", lambda f: committed[{v: k for k, v in rc.task_files().items()}[f]])
    monkeypatch.setattr(rc, "relabel_checks", lambda **k: {"reasons": [], "passes": True})
    monkeypatch.setattr(rc.gv, "stamp", lambda started, files: {"commit": "e" * 40, "tree_clean": True, "utc": "x"})
    monkeypatch.setattr(rc.gv, "tree_state", lambda: {"commit": "e" * 40, "clean": True})
    res = rc.read(out)
    printed = capsys.readouterr().out
    assert res["reading"]["reading_S"] == "AHEAD" and res["tries"] == 10
    assert "1.23" not in printed and "51." not in printed and "Outcome" in printed
    text = (out / "t62.txt").read_text()
    assert "+1.23" in text and "J6 range" in text and "mean scores none_gmom" in text
    with pytest.raises(SystemExit):
        rc.read(out)  # no second read into the same folder
