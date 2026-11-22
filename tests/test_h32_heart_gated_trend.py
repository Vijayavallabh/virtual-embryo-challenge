"""Try H32 (experiments/harness/h32_heart_gated_trend.py): synthetic arrays, hand-made tables and temporary files only.
Nothing here loads a competition or outside file, scores a real cell or reads the repository's ledger or register as
they stand."""

from __future__ import annotations

import importlib.util
import json
import os
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from vec import harness
from vec.adoption import SEEDS

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "experiments" / "harness" / "h32_heart_gated_trend.py"
LAUNCHER = REPO_ROOT / "experiments" / "harness" / "run_h32.sh"


@pytest.fixture(scope="module")
def h32():
    spec = importlib.util.spec_from_file_location("h32_heart_gated_trend_module", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def synth():
    """test_b0_trend's synthetic heart stages, with the released 'group' labels and a core flag."""
    spec = importlib.util.spec_from_file_location("b0_trend_tests_for_h32", REPO_ROOT / "tests" / "test_b0_trend.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    inputs = {"ref": mod._stage(1600, seed=2, shift=0.0, tiny=60),
              "start": mod._stage(1100, seed=3, shift=0.8, tiny=10),
              "earlier": mod._stage(900, seed=5, shift=1.6, tiny=10)}
    for a in inputs.values():
        m = np.asarray(a.obs["merged"]).astype(str)
        a.obs["core"] = m == "CM"
        a.obs["group"] = np.where(m == "Endothelial", "Endothelium", m)
    return mod, inputs


def _trimmed(h32, inputs):
    pool = inputs["ref"]
    PC = np.asarray(pool.obsm["spatial_3D"], dtype=float)
    keep = h32.h26.pool_masks(PC, np.asarray(pool.obs["core"], bool))[h32.FORM]
    return {**inputs, "ref": pool[keep].copy()}


def _fallback_inputs(inputs):
    """Endothelial usable at start and the reference only: the mapping leaves it out (B0's own factor)."""
    early = inputs["earlier"].copy()
    m = np.asarray(early.obs["merged"]).astype(str)
    early.obs["merged"] = np.where(m == "Endothelial", "Elsewhere", m)
    return {**inputs, "earlier": early}


def _table(genes, seed=0, lineages=("CM", "Endo")):
    rng = np.random.default_rng(seed)
    return {"task": "T2h/<=E9.5->E10.5", "genes": list(genes),
            "lineages": {k: list(np.clip(rng.normal(0, 0.6, len(genes)), -np.log(3), np.log(3))) for k in lineages}}


# ------------------------------------------------------------------------------------------------ the trend
@pytest.mark.parametrize("case", ["mapping", "fallback"])
def test_a_mask_of_ones_gives_u_bit_for_bit(h32, synth, case):
    mod, inputs = synth
    ins = _trimmed(h32, inputs if case == "mapping" else _fallback_inputs(inputs))
    ctx = mod._ctx(ins, seed=harness.BUILD_SEED_OFFSET + 4)
    parts = h32.trend_parts(ctx)
    b, cases = h32.base_factors(ctx, parts)
    assert cases["Endothelial"] == case and cases["CM"] == "mapping"
    st, _ = h32.express(ctx, parts, {k: v * np.ones(len(v)) for k, v in b.items()})
    ref, _, rows, _ = h32.h28.trend_expression(ctx, 1.0)
    np.testing.assert_array_equal(st.X, ref.X)
    np.testing.assert_array_equal(st.coords, ref.coords)
    np.testing.assert_array_equal(rows, parts["rows"])
    # a mask of zeros on CM is a different matrix (the hook is live)
    st0, _ = h32.express(ctx, parts, {"CM": np.zeros(len(b["CM"]))})
    assert not np.array_equal(st0.X, ref.X)


def test_a_type_skipped_or_absent_stays_out_of_every_arm(h32, synth):
    mod, inputs = synth
    ctx = mod._ctx(_trimmed(h32, inputs), seed=harness.BUILD_SEED_OFFSET + 1)
    parts = h32.trend_parts(ctx)
    b, cases = h32.base_factors(ctx, parts, types=("CM", "Tiny", "Nowhere"))
    assert cases == {"CM": "mapping", "Tiny": "skipped", "Nowhere": "absent"} and set(b) == {"CM"}


# ------------------------------------------------------------------------------------------------ the sets
def _vectors(n=400, seed=1):
    rng = np.random.default_rng(seed)
    b = rng.normal(0, 0.3, n)
    b[rng.random(n) < 0.1] = 0.0
    d = rng.normal(0, 0.5, n)
    basis = rng.normal(0, 1, n)
    return b, d, basis


def test_the_gate_sets(h32):
    b = np.array([0.2, -0.2, 0.2, 0.0, 0.3, -0.1])
    d = np.array([0.3, 0.3, 0.1, 0.9, -0.25, -0.25])
    s = h32.gate_sets(b, d)
    np.testing.assert_array_equal(s["E"], [True, True, False, False, True, True])
    np.testing.assert_array_equal(s["M"], [True, False, False, False, False, True])
    np.testing.assert_array_equal(s["Fl"], [False, True, False, False, True, False])


def test_r_matches_gs_count_per_bin_and_draws_only_from_e(h32):
    b, d, basis = _vectors()
    s = h32.gate_sets(b, d)
    draws = []
    for j in (1, 2, 3):
        m = h32.random_gate(s["E"], s["M"], basis, "CM", j)
        assert not np.any(m & ~s["E"])
        assert h32.bin_counts(m, s["E"], basis) == h32.bin_counts(s["M"], s["E"], basis)
        np.testing.assert_array_equal(m, h32.random_gate(s["E"], s["M"], basis, "CM", j))
        draws.append(m)
    assert not np.array_equal(draws[0], draws[1])
    assert not np.array_equal(draws[0], h32.random_gate(s["E"], s["M"], basis, "Endothelial", 1))


def test_f_is_gs_complement_in_e_and_p_has_no_fixed_point(h32):
    b, d, basis = _vectors()
    vecs, rec = h32.arm_vectors({"CM": b}, {"CM": d}, {"CM": basis})
    s = h32.gate_sets(b, d)
    r = rec["CM"]
    assert r["f_complement"] is True and r["eligible"] == s["E"].sum() == r["kept"] + r["flipped"]
    np.testing.assert_array_equal(vecs[h32.F]["CM"] != 0, s["Fl"])
    np.testing.assert_array_equal(vecs[h32.G]["CM"], b * s["M"])
    np.testing.assert_array_equal(vecs[h32.U]["CM"], b)
    for p in h32.PS:
        assert r["P_fixed_points"][p] == 0
        np.testing.assert_array_equal(np.sort(vecs[p]["CM"]), np.sort(vecs[h32.G]["CM"]))
    for nm in h32.RS:
        assert r["R"][nm]["inside_E"] and r["R"][nm]["bin_counts_match"] and r["R"][nm]["kept"] == r["kept"]
    perm = h32.h28.derangement(basis, "CM", bins=h32.N_BINS, seed=h32.P_SEED + 1)
    assert not np.any(perm == np.arange(len(perm)))
    # a type with no chord is in no arm
    vecs2, rec2 = h32.arm_vectors({"CM": b}, {"CM": None}, {"CM": basis})
    assert rec2 == {"CM": {"gated": False}} and all(not v for v in vecs2.values())


def test_the_chord_vector_aligns_genes_and_misses_an_absent_lineage(h32):
    table = {"genes": ["a", "b", "c"], "lineages": {"CM": [0.1, -0.5, 0.9]}}
    np.testing.assert_array_equal(h32.chord_vector(table, ["c", "z", "a"], "CM"), [0.9, 0.0, 0.1])
    assert h32.chord_vector(table, ["a"], "Endo") is None


# ------------------------------------------------------------------------------------------------ G0
def test_g0_passes_on_agreement_and_fails_on_noise_or_few_genes(h32):
    rng = np.random.default_rng(3)
    n = 600
    b = rng.normal(0, 0.3, n)
    basis = rng.normal(0, 1, n)
    d_agree = np.sign(b) * rng.uniform(0.3, 1.0, n) * np.where(rng.random(n) < 0.8, 1, -1)
    r = h32.g0_lineage(b, d_agree, basis, "CM", n_perm=200)
    assert r["n"] == n and r["a"] > 0.7 and h32.g0_outcome(r) == "PASS"
    d_noise = rng.choice([-1, 1], n) * rng.uniform(0.3, 1.0, n)
    assert h32.g0_outcome(h32.g0_lineage(b, d_noise, basis, "CM", n_perm=200)) == "FAIL"
    few = h32.g0_lineage(b[:40], d_agree[:40], basis[:40], "CM", n_perm=50)
    assert few["n"] == 40 and h32.g0_outcome(few) == "FAIL"
    none = h32.g0_lineage(b, np.zeros(n), basis, "CM", n_perm=10)
    assert none["n"] == 0 and h32.g0_outcome(none) == "FAIL"
    again = h32.g0_lineage(b, d_noise, basis, "CM", n_perm=200)
    assert again == h32.g0_lineage(b, d_noise, basis, "CM", n_perm=200)


# ------------------------------------------------------------------------------------------------ the arms
def test_every_arm_has_the_controls_cells_coordinates_and_the_checks_hold(h32, synth):
    mod, inputs = synth
    ins = _trimmed(h32, inputs)
    sp = h32.h27.sp
    geoms = [sp.geometry_from_arrays(n, t, np.asarray(inputs[r].obsm["spatial_3D"], float),
                                     np.asarray(inputs[r].obs["core"], bool))
             for r, n, t in (("earlier", "E8.5", 8.5), ("start", "E8.75", 8.75), ("ref", "E9.5", 9.5))]
    size = h32.h27.analogue_size_result(geoms, 10.5)
    ctx = mod._ctx(ins, seed=harness.BUILD_SEED_OFFSET + 3)
    out, extra = h32.build_arms(ctx, size, _table(mod.GENES), inputs["ref"])
    c = extra["checks"]
    assert set(out) == set(h32.BUILT)
    for k in ("u_equals_h28_trend_s1", "same_draw", "pre_size_coords_equal", "post_size_coords_equal"):
        assert c[k] is True, k
    assert len(set(c["n_cells"].values())) == 1
    assert c["cases"] == {"CM": "mapping", "Endothelial": "mapping"}
    assert c["chord_lineages"] == {"CM": "CM", "Endothelial": "Endo"}
    for k in ("CM", "Endothelial"):
        s = c["sets"][k]
        assert s["gated"] and s["f_complement"]
        assert all(r["inside_E"] and r["bin_counts_match"] for r in s["R"].values())
    for st in out.values():
        np.testing.assert_array_equal(st.coords, out[h32.CONTROL].coords)
    assert not np.array_equal(out[h32.G].X, out[h32.U].X)
    assert not np.array_equal(out[h32.G].X, out[h32.F].X)
    # an Endo chord missing from the table leaves Endothelial ungated in every arm
    _, extra2 = h32.build_arms(ctx, size, _table(mod.GENES, lineages=("CM",)), inputs["ref"])
    assert extra2["checks"]["sets"]["Endothelial"] == {"gated": False}
    assert extra2["checks"]["chord_lineages"]["Endothelial"] is None


# ------------------------------------------------------------------------------------------------ hand-made parts
METRICS = ("de_score", "de_direction", "mmd_u", "variogram", "d2_shape", "occupancy_dice", "scale_log_ratio",
           "neighborhood_mmd")
GROUPS = ("expression_change", "cell_state", "shape_growth", "local_organisation")


def _checks(h32, over=None):
    sets = {k: {"gated": True, "eligible": 100, "kept": 60, "flipped": 40, "f_complement": True,
                "R": {r: {"inside_E": True, "bin_counts_match": True, "kept": 60} for r in h32.RS},
                "P_fixed_points": {p: 0 for p in h32.PS}, "kept_bin_counts": [6] * 10} for k in h32.PAIR}
    c = {"u_equals_h28_trend_s1": True, "same_draw": True, "pre_size_coords_equal": True,
         "post_size_coords_equal": True, "cases": {k: "mapping" for k in h32.PAIR}, "sets": sets,
         "n_cells": {n: 350 for n in h32.BUILT}}
    c.update(over or {})
    return c


def _part(h32, key, seeds, score=None, de=None, group=None, commit="f" * 40, checks=None, sha="s1"):
    """score: arm -> board points added to 50 (scalar or per seed); de: arm -> added to de_score and de_direction;
    group: arm -> added to the expression_change and cell_state skills."""
    score, de, group = score or {}, de or {}, group or {}
    rows = []
    for s in seeds:
        sc, sk, gr = {}, {}, {}
        for q in h32.TRUTHS:
            sc[q], sk[q], gr[q] = {}, {}, {}
            for name in h32.ALL:
                sc[q][name] = 50.0 + float(np.broadcast_to(score.get(name, 0.0), (10,))[s])
                sk[q][name] = {m: 0.5 + (de.get(name, 0.0) if m in h32.DE_TERMS else 0.0) for m in METRICS}
                gr[q][name] = {g: 0.5 + (group.get(name, 0.0) if g in h32.DE_GROUPS else 0.0) for g in GROUPS}
        ck = _checks(h32, checks)
        rows.append({"seed": s, "n_cells": ck["n_cells"], "checks": ck, "scores": sc, "skills": sk, "groups": gr,
                     "records": {n: {} for n in h32.BUILT}})
    tr = h32.h30.TRANSITIONS[key]
    return {"axis": h32.AXIS, "task": key, "rows": rows, "table_sha256": sha,
            "table": {"task": key, "tau_practice": tr["tau"], "tau_test": h32.h30.TEST_TAU,
                      "classifiers": {"jax": {"trained_on": [tr["merfish"]]}, "feng": {"trained_on": [tr["merfish"]]}},
                      "stages_read": h32.h30._expected_table(key)},
            "window_guard": h32.WINDOW_BYPASS, "provenance": {"commit": commit, "tree_clean": True}}


def _refs(h32, u_shift=0.0):
    """Row H27's control and row H28's trend_s1 as committed parts would hold them (scores only)."""
    h27p, h28p = {}, {}
    for key in h32.TASK_KEYS:
        r27, r28 = [], []
        for s in SEEDS:
            r27.append({"seed": s, "scores": {q: {h32.H27_ARM: 50.0} for q in h32.TRUTHS}})
            r28.append({"seed": s, "scores": {q: {h32.h28.arm(1.0): 50.0 + 0.4 + u_shift} for q in h32.TRUTHS}})
        h27p[key] = [{"rows": r27}]
        h28p[key] = [{"rows": r28}]
    return h27p, h28p


def _parts(h32, **kw):
    return {k: [_part(h32, k, range(0, 5), **kw), _part(h32, k, range(5, 10), **kw)] for k in h32.TASK_KEYS}


SHA = {"e95.json": "s1", "e875.json": "s1"}
AHEAD_SCORES = {"U": 0.4, "G": 1.5, "R1": 0.9, "R2": 0.95, "R3": 1.0, "F": 0.2, "P1": 0.7, "P2": 0.7, "P3": 0.7}


def test_the_read_says_ahead_only_when_s1_to_s6_hold(h32):
    h27p, h28p = _refs(h32)
    res = h32.reading(_parts(h32, score=AHEAD_SCORES, de={"G": 0.01}, group={"G": 0.04}), h27p, h28p, SHA, "f" * 40,
                      True)
    assert res["void"] is False, res["why"]
    assert res["outcome"] == "AHEAD" and all(res["statements"]["both"].values())
    # G only 0.2 over mean R: S2 fails, NOT AHEAD
    sc = {**AHEAD_SCORES, "R1": 1.3, "R2": 1.3, "R3": 1.3}
    res = h32.reading(_parts(h32, score=sc, de={"G": 0.01}, group={"G": 0.04}), h27p, h28p, SHA, "f" * 40, True)
    assert res["outcome"] == "NOT AHEAD" and not res["statements"]["both"]["S2"]
    # the margin off the DE terms: S4 fails
    res = h32.reading(_parts(h32, score=AHEAD_SCORES), h27p, h28p, SHA, "f" * 40, True)
    assert res["outcome"] == "NOT AHEAD" and not res["statements"]["both"]["S4"]
    # G behind U: S5 fails; G not ahead of F: S6 fails
    res = h32.reading(_parts(h32, score={**AHEAD_SCORES, "U": 1.6, "F": 1.6}, de={"G": 0.01}, group={"G": 0.04}),
                      h27p, h28p, SHA, "f" * 40, True)
    st = res["statements"]["both"]
    assert res["outcome"] == "NOT AHEAD" and not st["S5"] and not st["S6"] and st["S1"]
    txt = h32.text(res)
    assert "OUTCOME: NOT AHEAD" in txt and "Row N56" in txt and "J6" in txt


@pytest.mark.parametrize("what,why", [
    ("seeds", "V1"), ("commit", "V2"), ("control", "V3"), ("u_bits", "V4"), ("u_score", "V4"), ("shape", "V5"),
    ("r_outside", "V6"), ("r_bins", "V6"), ("f", "V6"), ("sha", "V6"), ("e105", "V6"), ("bypass", "V6"),
    ("g0", "V6"), ("case", "V6")])
def test_each_void_condition(h32, what, why):
    h27p, h28p = _refs(h32)
    parts = _parts(h32, score=AHEAD_SCORES, de={"G": 0.01}, group={"G": 0.04})
    key = h32.TASK_KEYS[0]
    p, r = parts[key][0], parts[key][0]["rows"][0]
    commit, gate_ok, sha = "f" * 40, True, dict(SHA)
    if what == "seeds":
        parts[key][1]["rows"].pop()
    elif what == "commit":
        commit = "e" * 40
    elif what == "control":
        h27p = {k: [{"rows": [{"seed": s, "scores": {q: {h32.H27_ARM: 49.8} for q in h32.TRUTHS}}
                              for s in SEEDS]}] for k in h32.TASK_KEYS}
    elif what == "u_bits":
        r["checks"]["u_equals_h28_trend_s1"] = False
    elif what == "u_score":
        h28p = _refs(h32, u_shift=0.02)[1]
    elif what == "shape":
        r["skills"][h32.OWN][h32.G]["d2_shape"] += 1e-6
    elif what == "r_outside":
        r["checks"]["sets"]["CM"]["R"]["R2"]["inside_E"] = False
    elif what == "r_bins":
        r["checks"]["sets"]["Endothelial"]["R"]["R1"]["bin_counts_match"] = False
    elif what == "f":
        r["checks"]["sets"]["CM"]["f_complement"] = False
    elif what == "sha":
        sha["e95.json"] = "other"
    elif what == "e105":
        p["table"]["classifiers"]["jax"]["trained_on"] = ["T2/E10.5.h5ad"]
    elif what == "bypass":
        p["window_guard"] = "no"
    elif what == "g0":
        gate_ok = False
    elif what == "case":
        r["checks"]["cases"] = {"CM": "mapping"}
    res = h32.reading(parts, h27p, h28p, sha, commit, gate_ok)
    assert res["outcome"] == "VOID" and any(w.startswith(why) for w in res["why"]), res["why"]
    assert h32.text(res).startswith("VOID: ")


def test_the_fixed_numbers(h32):
    assert h32.T_D == 0.25 and h32.BAR == 0.3 and h32.MIN_AHEAD == 9 and h32.F_MIN_AHEAD == 8
    assert (h32.R_SEED, h32.P_SEED, h32.G0_SEED) == (20261130, 20261126, 20261122)
    assert h32.G0_PERM == 1000 and h32.G0_MIN_GENES == 50 and h32.G0_MIN_A == 0.60 and h32.N_BINS == 10
    assert h32.PAIR == {"CM": "CM", "Endothelial": "Endo"}
    assert h32.h30.LINEAGE_OF_GROUP["Endocardium"] == h32.h30.LINEAGE_OF_GROUP["Endothelium"] == "Endo"
    assert h32.FORM == "core_iso_half" and h32.U_SCORE_TOL == 0.01
    assert h32.COUNTED == ("G", "R1", "R2", "R3", "F", "P1", "P2", "P3")
    assert float(harness._momentum_default("eps")) == 0.5
    assert h32.TEST_TABLE_KEY == "T2h/<=E9.5->E10.5" and h32.h30.TRANSITIONS[h32.TEST_TABLE_KEY]["jax_early"] == (22,)
    assert h32.h30.symbols_of_rows(["x", "y"], {"x": "A"})[0] == ["A", None]  # the corrected reader is the one imported


# ------------------------------------------------------------------------------------------------ refusals and files
def test_the_register_and_row_refusals_on_temporary_files(h32, tmp_path):
    reg = tmp_path / "EXTERNAL_DATA.md"
    reg.write_text("| X002 | Feng. T2-heart amendment for row H32: practice read only. |\n| X003 | Qiu. |\n")
    assert h32.register_missing(reg) == ["X003"]
    reg.write_text("| X002 | Feng | T2-heart amendment for row H32 (practice). |\n"
                   "| X003 | Qiu. No T2-heart amendment for row H32 yet. T2-heart amendment for row H30. |\n")
    assert h32.register_missing(reg) == ["X003"]
    reg.write_text("| X002 | a. T2-heart amendment for row H32: x. |\n| X003 | b; T2-heart amendment for row H32. |\n")
    assert h32.register_missing(reg) == []
    tries = tmp_path / "TRIES.md"
    tries.write_text("| H30 | heart | x |\n")
    with pytest.raises(SystemExit, match="no row"):
        h32.refuse_rows_and_register(tries, reg)
    tries.write_text("| H32 (one row) | heart | x |\n")
    h32.refuse_rows_and_register(tries, reg)
    reg.write_text("| X002 | a. |\n| X003 | b. |\n")
    with pytest.raises(SystemExit, match="X002"):
        h32.refuse_rows_and_register(tries, reg)


def test_h30s_tables_are_reused_only_when_committed_and_the_reader_is_the_same(h32, tmp_path):
    folder = tmp_path / "h30tables"
    folder.mkdir()
    sums = {}
    for key in h32.TASK_KEYS:
        nm = f"{h32.h30.TRANSITIONS[key]['short']}.json"
        by = json.dumps({"task": key, "provenance": {"commit": "c" * 40}}).encode()
        (folder / nm).write_bytes(by)
        sums[nm] = h32.h30.sha256_bytes(by)
    committed = tmp_path / "sha.json"
    assert h32.h30_reusable(folder, committed, same=lambda c: True) is None  # nothing committed
    committed.write_text(json.dumps(sums))
    got = h32.h30_reusable(folder, committed, same=lambda c: True)
    assert got is not None and set(got) == set(sums)
    assert h32.h30_reusable(folder, committed, same=lambda c: False) is None  # another reader
    assert h32.h30_reusable(None, committed, same=lambda c: True) is None
    (folder / "e95.json").write_bytes(b"{}")
    assert h32.h30_reusable(folder, committed, same=lambda c: True) is None  # a sha256 off


def test_gate_readings_and_ancestry(h32, tmp_path):
    g = tmp_path / "gate.json"
    g.write_text(json.dumps({"axis": h32.AXIS + ":gate", "outcome": "FAIL",
                             "provenance": {"commit": "a" * 40, "tree_clean": True}}))
    assert not h32.gate_passed(g)
    g.write_text(json.dumps({"axis": h32.AXIS + ":gate", "outcome": "PASS",
                             "provenance": {"commit": "a" * 40, "tree_clean": True}}))
    assert h32.gate_passed(g) and h32.gate_passed(g, "a" * 40) and not h32.gate_passed(g, "b" * 40)
    assert not h32.committed_gate_ok("a" * 40, g)  # not a strict ancestor of itself
    assert h32.n56_outcome(tmp_path / "missing.json") == "not committed"


# ------------------------------------------------------------------------------------------------ the launcher
def _env(tmp_path):
    env = {k: v for k, v in os.environ.items() if k not in ("GIT_DIR", "GIT_WORK_TREE", "GIT_INDEX_FILE")}
    env.update(HOME=str(tmp_path), GIT_CEILING_DIRECTORIES=str(tmp_path))
    return env


def test_the_launcher_text(h32):
    sh = LAUNCHER.read_text()
    code = [line.split("#")[0] for line in sh.splitlines()]
    assert not any(re.search(r"(^|[\s;&|(])(rm|rmdir|unlink|shred|truncate)\s", line) for line in code)
    for msg in re.findall(r"\$\{[A-Za-z0-9_]+:\?([^}]*)\}", sh):
        assert "'" not in msg and '"' not in msg, msg
    for need in ("OMP_NUM_THREADS=8", "OPENBLAS_NUM_THREADS=8", "MKL_NUM_THREADS=8", "NUMPY_MADVISE_HUGEPAGE=0",
                 'TMPDIR="$S/tmp"', 'S="$HOME/v/scratch/h32"', 'UV_PROJECT_ENVIRONMENT="$PWD/.venv"',
                 "git rev-parse --show-toplevel", "status --porcelain", "python -I -B", "merge-base --is-ancestor",
                 "git diff --name-only", "0-1 2-3 4-5 6-7 8-9", "'^| H32 '", '"$SCRIPT" check',
                 'cmp -s runs/h32/gate/gate.json "$H32R/gate/gate.json"',
                 'cmp -s runs/h32/tables/tables_sha256.json "$H32R/tables_sha256.json"'):
        assert need in sh, need
    assert "/v/vec-campaign" not in sh.split("set -u", 1)[1]  # no hard-coded checkout
    order = dict(re.findall(r"^  (tests|tables|gate|run|read)\) PREV=\"?([a-z]*)\"? ;;$", sh, flags=re.M))
    assert order == {"tests": "", "tables": "tests", "gate": "tables", "run": "gate", "read": "run"}
    listed = re.findall(r"tests/test_[a-z0-9_]+\.py", sh)
    assert "tests/test_h32_heart_gated_trend.py" in listed and all((REPO_ROOT / f).exists() for f in listed)
    body = sh.split('echo running > "$S/$STEP.exit"')[0]
    assert "finish " not in body.split("finish() {")[1].split("\n", 1)[1]  # before a step runs: refusals only
    assert body.index("refuse 99") < body.index("refuse 98") < body.index("refuse 93")


def _repo(tmp_path, tries_row=True):
    repo = tmp_path / "repo"
    (repo / "experiments" / "harness").mkdir(parents=True)
    (repo / "experiments" / "harness" / "run_h32.sh").write_bytes(LAUNCHER.read_bytes())
    (repo / "experiments" / "harness" / "TRIES.md").write_text("| H32 (x) | heart |\n" if tries_row else "| H30 |\n")
    env = _env(tmp_path)
    git = ["git", "-c", "user.name=t", "-c", "user.email=t@t", "-c", "commit.gpgsign=false"]
    for argv in (["init", "-q"], ["add", "-A"], ["commit", "-q", "-m", "x"]):
        subprocess.run([*git, *argv], cwd=repo, env=env, check=True, capture_output=True)
    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=repo, env=env, capture_output=True, text=True,
                          check=True).stdout.strip()
    return repo, env, head


def _stub_uv(tmp_path, env, check_exit: int):
    """A stand-in `uv`: the check exits `check_exit`; anything else records the scratch folder's tests files."""
    s = tmp_path / "v" / "scratch" / "h32"
    bin_ = tmp_path / "bin"
    bin_.mkdir(exist_ok=True)
    (bin_ / "uv").write_text("#!/usr/bin/env bash\n"
                             f'for a in "$@"; do [ "$a" = check ] && exit {check_exit}; done\n'
                             f'cat "{s}/tests.exit" "{s}/tests.head" > "{tmp_path}/seen"\n'
                             f'echo "$UV_PROJECT_ENVIRONMENT" >> "{tmp_path}/seen"\n')
    (bin_ / "uv").chmod(0o755)
    env["PATH"] = f"{bin_}{os.pathsep}{env['PATH']}"


def test_the_launcher_runs_a_step_from_the_checkout_root_and_writes_running_first(tmp_path):
    repo, env, head = _repo(tmp_path)
    s = tmp_path / "v" / "scratch" / "h32"
    s.mkdir(parents=True)
    (s / "tests.exit").write_text("1\n")
    (s / "tests.head").write_text("1" * 40 + "\n")
    _stub_uv(tmp_path, env, 0)
    run = subprocess.run(["bash", str(repo / "experiments/harness/run_h32.sh"), "tests", head],
                         cwd=repo / "experiments", env=env, capture_output=True, text=True)
    assert run.returncode == 0, run.stderr
    seen = (tmp_path / "seen").read_text().split()
    assert seen == ["running", head, f"{repo}/.venv"]  # the checkout's own environment, from its root
    assert (s / "tests.exit").read_text().strip() == "0" and (s / "tests.head").read_text().strip() == head


@pytest.mark.parametrize("case,code", [("head", 99), ("dirty", 98), ("row", 95), ("register", 93)])
def test_the_launcher_refusals_never_touch_a_finished_exit(tmp_path, case, code):
    repo, env, head = _repo(tmp_path, tries_row=case != "row")
    s = tmp_path / "v" / "scratch" / "h32"
    s.mkdir(parents=True)
    _stub_uv(tmp_path, env, 1 if case == "register" else 0)
    for step in ("tests", "tables"):
        (s / f"{step}.exit").write_text("0\n")
        (s / f"{step}.head").write_text(head + "\n")
    want = "0" * 40 if case == "head" else head
    if case == "head":
        for step in ("tests", "tables"):
            (s / f"{step}.head").write_text(want + "\n")
    if case == "dirty":
        (repo / "experiments" / "harness" / "TRIES.md").write_text("| H32 (y) | heart |\n")
    run = subprocess.run(["bash", "experiments/harness/run_h32.sh", "gate", want], cwd=repo, env=env,
                         capture_output=True, text=True)
    assert run.returncode == code, run.stderr
    assert (s / "gate.refused.exit").read_text().strip() == str(code)
    assert not (s / "gate.exit").exists() and not (s / "gate.head").exists()
    assert (s / "tables.exit").read_text().strip() == "0"


def test_the_launcher_order_and_usage(tmp_path):
    repo, env, head = _repo(tmp_path)
    s = tmp_path / "v" / "scratch" / "h32"

    def call(*argv):
        return subprocess.run(["bash", "experiments/harness/run_h32.sh", *argv], cwd=repo, env=env,
                              capture_output=True, text=True)
    assert call().returncode == 2 and (s / "usage.refused.exit").read_text().strip() == "2"
    assert call("build", head).returncode == 97
    assert call("tables", head).returncode == 96  # tests has not exited 0
    assert call("run", head).returncode == 2  # run needs RUN
    assert call("read", head, "1" * 40).returncode == 96
    (s / "tests.exit").write_text("0\n")
    (s / "tests.head").write_text(head + "\n")
    assert call("tests", head).returncode == 94  # already exited 0 at this ROW
    assert (s / "tests.exit").read_text().strip() == "0"
