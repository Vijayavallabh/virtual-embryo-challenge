"""Row H29b (experiments/harness/h29b_heart_organ.py): synthetic arrays and hand-made tables only, and committed result
files of earlier rows read for arithmetic (the guard's dry run). Nothing here loads a competition file or scores a real
cell; the launchers are tested in tests/test_heart_organ_launchers.py."""

from __future__ import annotations

import ast
import importlib.util
import inspect
import json
import re
import sys
from pathlib import Path

import numpy as np
import pytest

from vec import harness, t2_build
from vec.adoption import SEEDS

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "experiments" / "harness" / "h29b_heart_organ.py"
RESULTS = REPO_ROOT / "experiments" / "harness" / "results"
COMMIT = "c" * 40


@pytest.fixture(scope="module")
def h29b():
    spec = importlib.util.spec_from_file_location("h29b_heart_organ_module", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


LABELS = np.array(["V-CM", "BEC", "NCC", "Epi"])  # cm (core), vascular_endothelium, surrounding, epicardium


@pytest.fixture(scope="module")
def synth():
    """test_b0_trend's synthetic T2 stages, with T2 labels by lineage: CM (the core) V-CM, Endothelial BEC, Untracked
    NCC, the tiny lineage Epi."""
    spec = importlib.util.spec_from_file_location("b0_trend_tests_for_h29b", REPO_ROOT / "tests" / "test_b0_trend.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    inputs = {"ref": mod._stage(1600, seed=2, shift=0.0, tiny=60),
              "start": mod._stage(1100, seed=3, shift=0.8, tiny=10),
              "earlier": mod._stage(900, seed=5, shift=1.6, tiny=10)}
    for a in inputs.values():
        merged = np.asarray(a.obs["merged"]).astype(str)
        a.obs["core"] = merged == "CM"
        kinds = {m: i for i, m in enumerate(mod.LINEAGES)}
        a.obs["celltype"] = LABELS[[kinds[m] for m in merged]]
    return mod, inputs


def _arrays(a):
    return (np.asarray(a.obsm["spatial_3D"], float), a.obs["celltype"].to_numpy(), np.asarray(a.obs["core"], bool))


def _size(h29b, inputs):
    sp = h29b.h27.sp
    geoms = [sp.geometry_from_arrays(n, t, np.asarray(inputs[r].obsm["spatial_3D"], float),
                                     np.asarray(inputs[r].obs["core"], bool))
             for r, n, t in (("earlier", "E8.5", 8.5), ("start", "E8.75", 8.75), ("ref", "E9.5", 9.5))]
    return h29b.h27.analogue_size_result(geoms, 10.5)


# ------------------------------------------------------------------------------------------------ the gate
def test_the_gate_verdict_at_its_edges(h29b):
    assert h29b.GATE_STOP_JACCARD == 0.90
    assert h29b.gate_verdict(0.8999999, 0.8999999) == "PASS"
    assert h29b.gate_verdict(0.90, 0.0) == "STOP" and h29b.gate_verdict(0.0, 0.90) == "STOP"
    assert h29b.gate_verdict(1.0, 1.0) == "STOP"
    assert h29b.GATE_EXIT == {"PASS": 0, "STOP": 3, "REFUSED": 4}
    a = np.array([1, 1, 0, 0], bool)
    assert h29b.jaccard(a, a) == 1.0 and h29b.jaccard(a, ~a) == 0.0
    assert h29b.jaccard(a, np.array([1, 0, 1, 0], bool)) == pytest.approx(1 / 3)


def test_the_gate_figures_on_a_synthetic_stage(h29b, synth):
    _, inputs = synth
    C, labels, core = _arrays(inputs["ref"])
    g = h29b.gate_figures(C, labels, core)
    keep, rec = t2_build.organ_scope(C, labels, core)
    c_pool, b_pool = t2_build.radial_trim_mask(C, core, 0.0), t2_build.radial_trim_mask(C, core, 0.5)
    assert g["pool_share"] == float(keep.mean()) and g["core_share"] == float(core[keep].mean())
    assert g["link_length"] == rec["link_length"] == t2_build.organ_link_length(C, core)
    assert g["jaccard_C"] == h29b.jaccard(keep, c_pool) and g["jaccard_B"] == h29b.jaccard(keep, b_pool)
    assert g["verdict"] == h29b.gate_verdict(g["jaccard_C"], g["jaccard_B"])
    # NCC is never kept, wherever it sits; the label pool keeps every chamber-category cell
    assert not keep[labels == "NCC"].any()
    assert g["label_pool_share"] == float((labels != "NCC").mean())


def test_the_gate_stops_when_d_is_the_full_trim(h29b):
    """A stage whose only chamber-category cells are the core: D's pool is the core, which the full trim keeps with
    the surrounding cells inside the core's own scope; far surrounding cells only. Jaccard(D, C) = 1: STOP."""
    rng = np.random.default_rng(0)
    core_xyz = rng.normal(size=(400, 3))
    far = rng.normal(size=(300, 3)) + np.array([60.0, 0.0, 0.0])
    C = np.vstack([core_xyz, far])
    core = np.r_[np.ones(400, bool), np.zeros(300, bool)]
    labels = np.r_[np.full(400, "V-CM"), np.full(300, "NCC")]
    g = h29b.gate_figures(C, labels, core)
    assert g["jaccard_C"] >= 0.9 and g["verdict"] == "STOP" and g["core_share"] == 1.0


def _gate_world(h29b, monkeypatch, stage):
    monkeypatch.setattr(h29b, "refuse_state", lambda commit, what="run": None)
    monkeypatch.setattr(h29b, "load_stage", lambda path, **k: stage)


def test_the_gate_step_writes_its_reading_and_prints_one_line(h29b, synth, monkeypatch, tmp_path, capsys):
    _, inputs = synth
    _gate_world(h29b, monkeypatch, inputs["ref"])
    code = h29b.run_gate(tmp_path / "g", tmp_path, COMMIT)
    g = json.loads((tmp_path / "g" / "gate.json").read_text())
    assert code == h29b.GATE_EXIT[g["verdict"]] and g["commit"] == COMMIT and g["stage"] == "T2/E10.5.h5ad"
    assert capsys.readouterr().out == f"GATE {g['verdict']}\n"
    txt = (tmp_path / "g" / "gate.txt").read_text()
    assert "not a try" in txt and h29b.GATE_LINES[g["verdict"]] in txt
    with pytest.raises(SystemExit, match="already holds a gate reading"):
        h29b.run_gate(tmp_path / "g", tmp_path, COMMIT)
    assert h29b.gate_passed(tmp_path / "g", COMMIT) == (g["verdict"] == "PASS")
    assert not h29b.gate_passed(tmp_path / "g", "d" * 40) and not h29b.gate_passed(tmp_path / "none", COMMIT)


def test_the_gate_refuses_an_unmapped_label(h29b, synth, monkeypatch, tmp_path, capsys):
    _, inputs = synth
    st = inputs["ref"].copy()
    lab = st.obs["celltype"].to_numpy().copy()
    lab[-3:] = "Not-A-Label"
    st.obs["celltype"] = lab
    _gate_world(h29b, monkeypatch, st)
    assert h29b.run_gate(tmp_path / "g", tmp_path, COMMIT) == 4
    g = json.loads((tmp_path / "g" / "gate.json").read_text())
    assert g["verdict"] == "REFUSED" and g["unmapped"] == {"Not-A-Label": 3} and "pool_share" not in g
    assert capsys.readouterr().out == "GATE REFUSED\n"
    assert not h29b.gate_passed(tmp_path / "g", COMMIT)


# ------------------------------------------------------------------------------------------------ the arms
def test_the_controls_are_row_h27s_f_arms_bit_for_bit(h29b, synth):
    mod, inputs = synth
    whole = inputs["ref"]
    PC, _, pcore = _arrays(whole)
    size = _size(h29b, inputs)
    pm = h29b.h26.pool_masks(PC, pcore)
    for form in h29b.FORMS:
        pool = whole if h29b.h26.TRIM[form] is None else whole[pm[form]].copy()
        ctx = mod._ctx({**inputs, "ref": pool}, seed=harness.BUILD_SEED_OFFSET + 4)
        mine, rec = h29b.build_control(form, ctx, {"even": 1.0}, size, {form: 1.0})
        ref = h29b.h27.build_forms(form, mod._ctx({**inputs, "ref": pool}, seed=harness.BUILD_SEED_OFFSET + 4),
                                   {"even": 1.0}, size, {f: 1.0 for f in h29b.FORMS})["f"]
        np.testing.assert_array_equal(mine.X, ref[0].X)
        np.testing.assert_array_equal(mine.coords, ref[0].coords)
        assert rec == ref[1]


def test_the_organ_arms_sizes_draw_and_half(h29b, synth):
    mod, inputs = synth
    whole = inputs["ref"]
    PC, labels, pcore = _arrays(whole)
    size = _size(h29b, inputs)
    masks, _ = h29b.organ_pools(PC, labels, pcore)
    assert (masks["organ"] <= masks["label"]).all() and masks["label"].sum() == (labels != "NCC").sum()
    ratios = {p: t2_build.organ_core_rms_ratio(PC, masks[p], pcore) for p in h29b.POOLS}
    ft = h29b.h27.forecast_targets(size)
    pt = {"core_iso_half": 3.0}
    targets = h29b.d_targets(ft, pt, ratios)
    assert targets["organ:f"]["target"] == ft["core_iso_half"] * ratios["organ"]
    assert targets["organ:p"]["target"] == 3.0 * ratios["organ"]
    assert targets["label:f"]["target"] == ft["core_iso_half"] * ratios["label"]
    seed = 6
    st = {}
    for pool in h29b.POOLS:
        ctx = mod._ctx({**inputs, "ref": whole[masks[pool]].copy()}, seed=harness.BUILD_SEED_OFFSET + seed)
        stages, recs, draw = h29b.build_pool_arms(pool, ctx, targets)
        st.update(stages)
        for arm, rec in recs.items():
            assert abs(rec["log_rms_after_vs_target"]) <= 1e-9 and rec["target_rms_full"] == targets[arm]["target"]
        assert draw["pool_cells"] == int(masks[pool].sum()) and draw["drawn"] == mod.BOARD.max_cells
        assert draw["short"] == (mod.BOARD.max_cells > masks[pool].sum())
        # the draw before the size step and the half: harness._copy on the same context, then the same half rows
        raw = harness._copy(mod._ctx({**inputs, "ref": whole[masks[pool]].copy()}, seed=harness.BUILD_SEED_OFFSET + seed))
        rows = h29b.h25.kept_rows(raw.n, seed, 1)
        np.testing.assert_array_equal(stages[f"{pool}:f"].X, np.asarray(raw.X)[rows])
    assert set(st) == set(h29b.COUNTED)
    # D_p and D_f: one draw and one half, positions differing by one scale about the same centroid
    np.testing.assert_array_equal(st["organ:p"].X, st["organ:f"].X)
    f, p = st["organ:f"].coords, st["organ:p"].coords
    k = targets["organ:p"]["target"] / targets["organ:f"]["target"]
    full_f = f - f.mean(0)
    np.testing.assert_allclose(p - p.mean(0), full_f * k, rtol=1e-5, atol=1e-4)


def test_d_f_is_size_organ_scopes_arithmetic(h29b, synth, monkeypatch):
    """The builder's size step for D (vec.t2_build.size_organ_scope) is rescale_to_rms to the core forecast times the
    ratio once its E10.5 checks pass (stood in here by a record with the forecast); the practice arm's step is the
    same arithmetic on the same cells."""
    _, inputs = synth
    C, _, _ = _arrays(inputs["ref"])
    fc, ratio = 7.5, 1.37
    monkeypatch.setattr(t2_build, "size_or_hold", lambda *a, **k: (None, {"target_rms_full": fc,
                                                                          "target_source": "stand-in"}))
    built, rec = t2_build.size_organ_scope(C, {}, ratio, source_stage="E10.5", reference_cloud=None)
    mine, mrec = h29b.size_to(C, fc * ratio, "f", ratio, fc, h29b.F_SOURCE)
    np.testing.assert_array_equal(built, mine)
    assert rec["target_rms_full"] == mrec["target_rms_full"]


# ------------------------------------------------------------------------------------------------ hand-made tables
GROUPS = ("expression_change", "cell_state", "shape_growth", "local_organisation")
METRICS = ("de_score", "de_direction", "mmd_u", "variogram", "neighborhood_mmd", "d2_shape", "occupancy_dice",
           "scale_log_ratio")
BASE = {"full_iso_half:f": 50.0, "half_iso_half:f": 51.0, "core_iso_half:f": 52.0, "organ:f": 53.0,
        "organ:p": 54.0, "label:f": 52.5, "floor": 50.0}


def _row(h29b, key, s, score=None):
    score = score or {}
    sc = {q: {a: BASE[a] + score.get(a, 0.0) + 0.01 * s + 0.001 * h29b.TRUTHS.index(q) * (a == "organ:f")
              for a in h29b.ALL} for q in h29b.TRUTHS}
    sk = {q: {a: {m: 0.5 for m in METRICS} for a in h29b.ALL} for q in h29b.TRUTHS}
    gr = {q: {a: {g: 0.5 for g in GROUPS} for a in h29b.ALL} for q in h29b.TRUTHS}
    recs = {a: {"log_rms_after_vs_target": 0.0, "target_rms_full": 2.0} for a in h29b.BUILT}
    draws = {p: {"pool_cells": 900, "drawn": 800, "short": False, "repeated_rows_at_least": 0} for p in h29b.POOLS}
    return {"seed": s, "build_seed": harness.BUILD_SEED_OFFSET + s, "n_cells": {a: 400 for a in h29b.BUILT},
            "size_records": recs, "draws": draws, "scores": sc, "skills": sk, "groups": gr}


def _info(h29b):
    t = {q: {"cells": 1000, "core_share": 0.4, "xy_scope": 1.0} for q in h29b.TRUTHS}
    p = {f: {"pool_cells": 900, "core_share": 0.5, "xy_scope": 1.0} for f in (*h29b.FORMS, *h29b.POOLS)}
    return t, p


def _part(h29b, key, seeds, commit=COMMIT, clean=True, score=None):
    t, p = _info(h29b)
    tg = {a: {"target": 2.0, "base": 1.0, "ratio": 2.0, "source": "x"} for a in h29b.COUNTED}
    return {"axis": h29b.AXIS, "task": key, "unmapped": None, "truth_info": t, "pred_info": p,
            "forecast_targets": {f: 2.0 for f in h29b.FORMS}, "perfect_targets": {f: 2.0 for f in h29b.FORMS},
            "ratios": {"organ": 1.5, "label": 1.8}, "d_targets": tg, "v7": {"ok": True, "log_diff": 0.0},
            "organ_record": {"reference": {"link_length": 1.0}, "target": {"link_length": 1.0}},
            "rows": [_row(h29b, key, s, score) for s in seeds],
            "provenance": {"commit": commit, "tree_clean": clean, "utc": "x"}}


def _parts(h29b, **kw):
    return {k: [_part(h29b, k, (0, 1, 2, 3, 4), **kw), _part(h29b, k, (5, 6, 7, 8, 9), **kw)]
            for k in h29b.TASK_KEYS}


def _h27(h29b):
    """Row H27's committed parts and reading as the hand-made controls make them (the same scores, exactly)."""
    t, p = _info(h29b)
    parts = {}
    for k in h29b.TASK_KEYS:
        rows = [{"seed": s, "scores": {q: {c: _row(h29b, k, s)["scores"][q][c] for c in h29b.CONTROLS}
                                       for q in h29b.H27_TRUTHS}} for s in SEEDS]
        parts[k] = [{"truth_info": {q: t[q] for q in h29b.H27_TRUTHS}, "pred_info": {f: p[f] for f in h29b.FORMS},
                     "forecast_targets": {f: 2.0 for f in h29b.FORMS}, "rows": rows}]
    merged = {k: parts[k][0]["rows"] for k in h29b.TASK_KEYS}
    tables = {w: {u: {"mean": {"f": {h29b.LETTER[f]: float(h29b.h27._series(merged, fs, ts, f, "f").mean())
                                     for f in h29b.FORMS}}} for u, fs in h29b.UNITS.items()} for w, ts in
              h29b.h27.WORLDS}
    return parts, {"void": False, "tables": tables}


def _gate(verdict="PASS", commit=COMMIT):
    return {"verdict": verdict, "commit": commit, "pool_share": 0.3, "core_share": 0.7, "link_length": 1.0,
            "link_percentile": 99.0, "jaccard_C": 0.5, "jaccard_B": 0.4, "label_pool_share": 0.5,
            "C_pool_share": 0.4, "B_pool_share": 0.6, "components_kept": 1, "components_total": 3,
            "label_config_sha256": "0" * 64}


def _read(h29b, parts=None, h27=None, gate="default", commit=COMMIT):
    hp, hr = h27 or _h27(h29b)
    return h29b.reading(parts if parts is not None else _parts(h29b), hp, hr, _gate() if gate == "default" else gate,
                        commit)


def test_a_complete_clean_read_is_read_with_its_stated_differences(h29b):
    res = _read(h29b)
    assert res["why"] == [] and res["outcome"] == "READ" and res["void"] is False
    c = res["tables"]["own"]["mean of folds"]
    assert c["stated"]["organ:f - core_iso_half:f"]["mean"] == pytest.approx(1.0)
    assert c["stated"]["organ:f - half_iso_half:f"]["mean"] == pytest.approx(2.0)
    assert c["stated"]["organ:f - full_iso_half:f"]["mean"] == pytest.approx(3.0)
    assert c["stated"]["organ:f - core_iso_half:f"]["ahead"] == 10
    assert c["diagnostic"]["label:f - organ:f"]["mean"] == pytest.approx(-0.5)
    assert c["diagnostic"]["organ:p - organ:f"]["mean"] == pytest.approx(1.0)
    assert set(res["tables"]) == {w for w, _ in h29b.WORLDS} and "truth organ" in res["tables"]
    assert res["headline_world"] == "own" and h29b.STATED_PAIRS[0] == ("organ:f", "core_iso_half:f")
    txt = h29b.text(res)
    assert txt.startswith("OUTCOME H29b: READ\n" + h29b.FOLLOWS["READ"] + "\n")
    for need in ("tautology", "cannot test D's bet", "Rule 2 caveat", "J5", "HEADLINE (J4)"):
        assert need in txt, need
    own = [ln for ln in txt.splitlines() if "stages] own " in ln]
    assert len(own) == len(h29b.UNITS) and all(ln.startswith("HEADLINE [practice") for ln in own)
    assert not any(ln.startswith("HEADLINE [") for ln in txt.splitlines() if "stages] own " not in ln)


def test_paired_se_and_seeds_ahead(h29b):
    s = h29b._stats(np.array([1.0, 2.0, 3.0, -1.0]))
    assert s["mean"] == 1.25 and s["ahead"] == 3 and s["lo"] == -1.0 and s["hi"] == 3.0
    assert s["se"] == pytest.approx(np.std([1, 2, 3, -1], ddof=1) / 2)


def test_every_void_condition(h29b):
    def why(**kw):
        return _read(h29b, **kw)["why"]

    def bad(edit):
        p = _parts(h29b)
        edit(p)
        return why(parts=p)
    assert any(w.startswith("G: no gate") for w in why(gate=None))
    assert any("'STOP'" in w for w in why(gate=_gate("STOP")))
    assert any("'REFUSED'" in w for w in why(gate=_gate("REFUSED")))
    assert any("not of the run commit" in w for w in why(gate=_gate(commit="d" * 40)))
    assert any("missing transition" in w for w in bad(lambda p: p.pop(h29b.TASK_KEYS[1])))
    assert any("seeds" in w for w in bad(lambda p: p[h29b.TASK_KEYS[0]][1]["rows"].pop()))
    assert any("occur in more than one part" in w for w in bad(lambda p: p[h29b.TASK_KEYS[0]][1]["rows"].append(
        p[h29b.TASK_KEYS[0]][0]["rows"][0])))
    assert any("organ/" in w for w in bad(lambda p: p[h29b.TASK_KEYS[2]][0]["rows"][0]["scores"].pop("organ")))
    assert any("label:f" in w for w in bad(lambda p: p[h29b.TASK_KEYS[0]][0]["rows"][3]["skills"]["own"].pop(
        "label:f")))
    assert any("draws missing" in w for w in bad(lambda p: p[h29b.TASK_KEYS[0]][0]["rows"][0].pop("draws")))
    assert any("different commits" in w for w in bad(lambda p: p[h29b.TASK_KEYS[0]][0]["provenance"].update(
        commit="d" * 40)))
    assert any("dirty tree" in w for w in why(parts=_parts(h29b, clean=False)))
    assert any("run commit" in w for w in why(commit="d" * 40))
    assert any("run commit" in w for w in why(commit=COMMIT[:12]))

    def unmapped(p):
        p[h29b.TASK_KEYS[1]] = [{"axis": h29b.AXIS, "task": h29b.TASK_KEYS[1], "rows": [],
                                 "unmapped": {"where": "reference", "labels": {"New-Label": 4}},
                                 "provenance": {"commit": COMMIT, "tree_clean": True}}]
    assert any("without a T2 chamber category" in w and "New-Label" in w for w in bad(unmapped))
    assert any("v7" in w or "reference cloud" in w for w in bad(lambda p: p[h29b.TASK_KEYS[0]][0]["v7"].update(
        ok=False)))
    assert any("truth radial_0" in w for w in bad(lambda p: p[h29b.TASK_KEYS[0]][0]["truth_info"]["radial_0"].update(
        cells=999)))
    assert any("core_iso_half's pool" in w for w in bad(
        lambda p: p[h29b.TASK_KEYS[0]][0]["pred_info"]["core_iso_half"].update(core_share=0.6)))
    assert any("forecast size" in w for w in bad(
        lambda p: p[h29b.TASK_KEYS[0]][0]["forecast_targets"].update(core_iso_half=2.0 + 1e-9)))
    assert any("disagree on the ratios" in w for w in bad(lambda p: p[h29b.TASK_KEYS[0]][1]["ratios"].update(
        organ=1.6)))
    assert any("not at its target" in w for w in bad(
        lambda p: p[h29b.TASK_KEYS[0]][0]["rows"][0]["size_records"]["organ:f"].update(log_rms_after_vs_target=1e-6)))
    assert any("expression skills" in w for w in bad(
        lambda p: p[h29b.TASK_KEYS[0]][0]["rows"][0]["skills"]["own"]["organ:p"].update(de_score=0.6)))
    assert any("another cell number" in w for w in bad(
        lambda p: p[h29b.TASK_KEYS[0]][0]["rows"][0]["n_cells"].update({"organ:p": 401})))
    # a malformed part is void, never a crash
    assert any("cannot be read" in w for w in bad(lambda p: p[h29b.TASK_KEYS[0]][0].pop("v7")))


def test_the_guard_is_bit_for_bit(h29b):
    """J7: a control one part in 1e12 from row H27's committed score voids the read; so does a reading mean."""
    def nudge(p):
        r = p[h29b.TASK_KEYS[2]][1]["rows"][2]
        r["scores"]["thin_0"]["half_iso_half:f"] += 1e-12
    p = _parts(h29b)
    nudge(p)
    w = _read(h29b, parts=p)["why"]
    assert w and all(x.startswith("J7: ") for x in w) and "half_iso_half:f" in w[0]
    hp, hr = _h27(h29b)
    hr["tables"]["CLEAN"]["fold E8.5"]["mean"]["f"]["C"] += 1e-12
    w = _read(h29b, h27=(hp, hr))["why"]
    assert w == [w[0]] and w[0].startswith("J7: CLEAN fold E8.5 C_f mean")
    assert _read(h29b, h27=(hp, {"void": True}))["why"] == ["J7: row H27's committed reading is void or unreadable"]
    # the organ truth is not row H27's: the guard compares row H27's five truths only
    assert h29b.H27_TRUTHS == ("own", "radial_0.5", "radial_0", "thin_0.5", "thin_0")


def test_a_short_pool_is_reported_never_voided(h29b):
    p = _parts(h29b)
    for ps in p.values():
        for part in ps:
            for r in part["rows"]:
                r["draws"]["organ"] = {"pool_cells": 300, "drawn": 800, "short": True, "repeated_rows_at_least": 500}
    res = _read(h29b, parts=p)
    assert res["outcome"] == "READ"
    sh = res["pools_short"][h29b.TASK_KEYS[0]]["organ"]
    assert sh == {"pool_cells": 300, "drawn": [800], "short_in_seeds": 10, "repeated_rows_at_least": 500}
    assert "short in 10 of 10 seeds" in h29b.text(res)


def test_the_void_text_and_the_printed_lines_hold_no_figure(h29b):
    res = _read(h29b, gate=None)
    assert res["outcome"] == "VOID" and "tables" not in res
    assert h29b.text(res).startswith("OUTCOME H29b: VOID\n" + h29b.FOLLOWS["VOID"] + "\nVOID: G: no gate reading")
    for line in h29b.FOLLOWS.values():
        assert not re.search(r"\d", re.sub(r"H29b|DGX|sha256", "", line)), line


# ------------------------------------------------------------------------------------------------ the guard's dry run
def test_the_guard_dry_runs_on_committed_results(h29b):
    """Rule 9 J7, before any run: the bit-for-bit guard's code on committed data other than the pair it will read.
    Row H27's rebuilt oracle twins of A, B and C against row H26's committed A, B and C, every transition, seed and
    truth: equal exactly. And row H27's committed parts against its own committed reading.json: every plain mean at
    the forecast size equal exactly. A copy shifted by one part in 1e12 fails both."""
    h27p = h29b.load_h27()
    h26p = h29b.h26.load_parts(h29b.H26_PARTS)
    merged = {}
    for key in h29b.TASK_KEYS:
        rows = h29b.h26.merge_rows(h27p[key])[0]
        merged[key] = rows
        assert [r["seed"] for r in rows] == list(SEEDS)
        pairs = tuple((h29b.h27.arm_name(f, "o"), f) for f in h29b.FORMS)
        assert h29b.guard_scores(key, rows, h29b.h26.merge_rows(h26p[key])[0], pairs) == [], key
    assert h29b.guard_reading(merged, h29b.load_h27_reading()) == []
    # the guard can fail
    key = h29b.TASK_KEYS[0]
    rows = json.loads(json.dumps(merged[key]))
    rows[4]["scores"]["radial_0"]["core_iso_half:o"] *= 1 + 1e-12
    pairs = (("core_iso_half:o", "core_iso_half"),)
    w = h29b.guard_scores(key, rows, h29b.h26.merge_rows(h26p[key])[0], pairs)
    assert len(w) == 1 and "radial_0 seed 4" in w[0]
    shifted = {**merged, key: rows}
    rows[4]["scores"]["own"]["core_iso_half:f"] += 1e-9
    assert h29b.guard_reading(shifted, h29b.load_h27_reading())


# ------------------------------------------------------------------------------------------------ refusals and prints
def _prints(fn):
    tree = ast.parse(inspect.getsource(fn).lstrip())
    return [ast.unparse(n) for n in ast.walk(tree) if isinstance(n, ast.Call) and getattr(n.func, "id", "") == "print"]


def test_each_step_prints_only_its_step_lines(h29b):
    assert _prints(h29b.run_transition) == ['print(f"{key} seed {s} ({row[\'seconds\']} s) done", flush=True)']
    assert _prints(h29b._write) == ["print(f'wrote {path}', flush=True)"]
    for fn in (h29b.build_pool_arms, h29b.build_control, h29b.reading, h29b.text, h29b.void_reasons, h29b.cmd_run,
               h29b.gate_figures, h29b.organ_pools):
        assert _prints(fn) == [], fn
    assert _prints(h29b.run_gate) == ["print(GATE_LINES[g['verdict']].split(':')[0], flush=True)"]
    assert sorted(_prints(h29b.cmd_read)) == sorted(['print(f"OUTCOME {ROW}: {res[\'outcome\']}")',
                                                     "print(res['follows'])"])
    assert {v.split(":")[0] for v in h29b.GATE_LINES.values()} == {"GATE PASS", "GATE STOP", "GATE REFUSED"}


def test_the_script_names_no_test_stage_file_and_reads_no_test_score(h29b):
    src = SCRIPT.read_text()
    for bad in ("results/p3", "/p3", "size_heart_E12", "submissions", "leaderboard", "registry", "build_test",
                "boards/", "upload", "12.5", "manual_uploads", "index.json", "test_extrap", "size_or_hold("):
        assert bad not in src, bad


def _refusing(h29b, monkeypatch):
    def boom(*a, **k):
        raise AssertionError("something loaded")
    monkeypatch.setattr(h29b.payoff, "heart_task", boom)
    monkeypatch.setattr(h29b, "load_stage", boom)
    monkeypatch.setattr(h29b, "load_parts", boom)
    monkeypatch.setattr(h29b, "load_h27", boom)


def test_gate_run_and_read_refuse_before_anything_loads(h29b, monkeypatch, tmp_path):
    _refusing(h29b, monkeypatch)
    g = tmp_path / "gate"
    g.mkdir()
    (g / "gate.json").write_text(json.dumps(_gate()))
    key = h29b.TASK_KEYS[0]
    run = ["run", "--task", key, "--seeds", "0-1", "--out", str(tmp_path / "o"), "--gate", str(g), "--commit", COMMIT]
    read = ["read", "--runs", str(tmp_path), "--gate", str(g), "--out", str(tmp_path / "r"), "--commit", COMMIT]
    gate = ["gate", "--out", str(tmp_path / "go"), "--commit", COMMIT]
    monkeypatch.setattr(h29b, "tries_rows", lambda path=h29b.TRIES: {"H27", "H29a", "R21"})
    for args in (run, read, gate):
        with pytest.raises(SystemExit, match="no row"):
            h29b.main(args)
    monkeypatch.setattr(h29b, "tries_rows", lambda path=h29b.TRIES: {"H27", "H29b"})
    monkeypatch.setattr(h29b.provenance, "git_state", lambda repo, ignore=(): (COMMIT, ["src/x.py"]))
    for args in (run, read, gate):
        with pytest.raises(SystemExit, match="not clean"):
            h29b.main(args)
    monkeypatch.setattr(h29b.provenance, "git_state", lambda repo, ignore=(): ("unknown", None))
    with pytest.raises(SystemExit, match="not clean"):
        h29b.main(run)
    monkeypatch.setattr(h29b.provenance, "git_state", lambda repo, ignore=(): ("d" * 40, []))
    for args in (run, read, gate):
        with pytest.raises(SystemExit, match="HEAD is not the run commit"):
            h29b.main(args)
    monkeypatch.setattr(h29b.provenance, "git_state", lambda repo, ignore=(): (COMMIT, []))
    for other in (COMMIT[:12], ""):
        with pytest.raises(SystemExit, match="HEAD is not the run commit"):
            h29b.run_transition(key, (0,), tmp_path / "o", tmp_path, other, g)
    # the gate of the commit must have said PASS
    for verdict, commit in (("STOP", COMMIT), ("REFUSED", COMMIT), ("PASS", "d" * 40)):
        (g / "gate.json").write_text(json.dumps(_gate(verdict, commit)))
        with pytest.raises(SystemExit, match="has not said PASS"):
            h29b.main(run)
    (g / "gate.json").write_text(json.dumps(_gate()))
    (tmp_path / "o").mkdir()
    (tmp_path / "o" / h29b.PART).write_text("{}")
    with pytest.raises(SystemExit, match="exists"):
        h29b.main(run)
    with pytest.raises(SystemExit, match="must be one of"):
        h29b.main(["run", "--task", "T2h/other", "--seeds", "0-1", "--out", str(tmp_path / "o2"), "--gate", str(g),
                   "--commit", COMMIT])
    with pytest.raises(AssertionError, match="something loaded"):  # past every refusal, the first load
        h29b.main([*run[:6], str(tmp_path / "o3"), *run[7:]])
    (tmp_path / "r").mkdir()
    for f in ("reading.txt", "reading.json"):
        (tmp_path / "r" / f).write_text("x")
        with pytest.raises(SystemExit, match="already holds"):
            h29b.main(read)
        (tmp_path / "r" / f).rename(tmp_path / f"moved_{f}")
    with pytest.raises(AssertionError, match="something loaded"):
        h29b.main(read)
    with pytest.raises(AssertionError, match="something loaded"):
        h29b.main(gate)


def test_the_fixed_constants(h29b):
    assert h29b.TASK_KEYS == ("T2h/<=E9.5->E10.5", "T2h/E8.25_late->E8.5", "T2h/<=E8.75->E10.5")
    assert h29b.TRUTHS == ("own", "radial_0.5", "radial_0", "thin_0.5", "thin_0", "organ")
    assert h29b.CONTROLS == ("full_iso_half:f", "half_iso_half:f", "core_iso_half:f")
    assert h29b.COUNTED == ("organ:f", "organ:p", "label:f")
    assert t2_build.ORGAN_LINK_PERCENTILE == 99.0 and h29b.HEADLINE == "own"
    assert h29b.h27.TASK_KEYS is h29b.TASK_KEYS  # imported, not copied
