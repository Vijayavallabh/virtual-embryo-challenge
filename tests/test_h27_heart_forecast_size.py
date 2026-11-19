"""experiments/harness/h27_heart_forecast_size.py (try H27 of experiments/harness/TRIES.md) on synthetic arrays and
hand-made result tables. Nothing is read from data/ or runs/ and nothing is scored. Committed result files are read
for arithmetic only: the size procedure's pre-release results (stage_geometry.csv, size_heart_modes.csv,
size_heart_pseudo.csv, size_heart_scope.csv), row H25's committed reading and task files, and row H26's committed parts
(the guards' dry run and the noise floor). No forecast is formed for any released target here."""

from __future__ import annotations

import ast
import csv
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
SCRIPT = REPO_ROOT / "experiments" / "harness" / "h27_heart_forecast_size.py"
COORDS = REPO_ROOT / "experiments" / "coords" / "results"
RESULTS = REPO_ROOT / "experiments" / "harness" / "results"


@pytest.fixture(scope="module")
def h27():
    spec = importlib.util.spec_from_file_location("h27_heart_forecast_size_module", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def _geom(h27, name, t, core, full, n=1000):
    """A synthetic stage geometry: only the RMS radii the raw mode reads are meaningful."""
    return h27.sp.StageGeometry(name=name, t=float(t), n_cells=n, n_core=n // 4, n_levels=40, z_step=10.0,
                                z_range=390.0, rms_full=float(full), rms_core=float(core),
                                rms_xy_full=float(full) * 0.8,
                                sd_z_full=float(full) * 0.6, rms_xy_core=float(core) * 0.8, sd_z_core=float(core) * 0.6,
                                core_thin_axis_cos_z=0.5)


def _four(h27):
    """Four synthetic stages on which the pre-registered tie_tol 0.15 would keep half_rate alone (mean |error| 0.667,
    0.300 and 0.620 for latest_rate, half_rate and power_t), with every rate interval before E9.0."""
    logs = (5.0, 5.5, 5.6, 6.5)
    days = (8.25, 8.5, 8.75, 9.5)
    return [_geom(h27, f"S{i}", t, np.exp(v), np.exp(v) * (1.9 - 0.2 * i)) for i, (t, v) in enumerate(zip(days, logs,
                                                                                                      strict=True))]


def _cloud_with_core(n=3000, seed=0):
    rng = np.random.default_rng(seed)
    core = rng.random(n) < 0.25
    sd = np.where(core[:, None], [[60.0, 50.0, 30.0]], [[300.0, 220.0, 60.0]])
    C = rng.normal(size=(n, 3)) * sd + np.array([500.0, -200.0, 40.0])
    return C, core


# ------------------------------------------------------------------------------------------------ the forecast
def test_the_analogue_lifts_the_cut_off_and_averages_all_three_rules(h27):
    geoms = _four(h27)
    sp = h27.sp
    pre = sp.heart_size_from_geometry(geoms, 10.5, z_mode="raw")  # the pre-registered settings
    assert pre["valid"] is False  # every rate interval starts before E9.0
    picky = sp.heart_size_from_geometry(geoms, 10.5, z_mode="raw", min_start_t=-np.inf)  # tie_tol 0.15
    assert picky["growth"]["core"]["chosen"] == ["half_rate"]  # selection would drop two rules: not the test's outcome
    size = h27.analogue_size_result(list(reversed(geoms)), 10.5)
    r = size["z_modes"]["raw"]
    assert size["valid"] is True and r is size["preregistered"]
    assert r["growth"]["core"]["chosen"] == list(h27.HEART_RULES) == ["latest_rate", "half_rate", "power_t"]
    assert r["min_start_t"] == -np.inf and r["tie_tol"] == np.inf and r["z_mode"] == "raw"
    ts, ys = np.array([g.t for g in geoms]), np.log([g.rms_core for g in geoms])
    want = np.exp(np.mean([sp.RULES[k](ts, ys, 10.5) for k in h27.HEART_RULES]))
    assert r["pred_rms_core"] == pytest.approx(want, rel=1e-12)
    assert r["reference"] == "S3"
    assert size["provenance"]["analogue"] == h27.ANALOGUE_RATES
    rec = h27.forecast_record(size, 10.5)
    assert rec["interval"] == ["S2", "S3"] and rec["gap"] == pytest.approx(1.0)
    assert rec["h"] == pytest.approx(1.0 / 0.75)


def test_one_input_stage_gives_the_copy_last_analogue(h27):
    g = _geom(h27, "E8.25_late", 8.25, 150.0, 350.0)
    size = h27.analogue_size_result([g], 8.5)
    r = size["z_modes"]["raw"]
    s = np.log(350.0 / 150.0)
    assert r["pred_rms_core"] == pytest.approx(150.0, rel=1e-12)
    assert r["scope"]["H1"]["pred_log_scope"] == pytest.approx(s) == r["scope"]["H0"]["pred_log_scope"]
    assert r["scope"]["H2"]["pred_log_scope"] == pytest.approx(s / 2)
    assert r["scope"]["H1"]["pred_rms_full"] == pytest.approx(350.0, rel=1e-12)
    assert r["scope"]["H2"]["pred_rms_full"] == pytest.approx(150.0 * np.exp(s / 2), rel=1e-12)
    assert r["growth"]["core"]["chosen"] == ["copy_last"]
    assert size["provenance"]["analogue"] == h27.ANALOGUE_COPY
    rec = h27.forecast_record(size, 8.5)
    assert rec["interval"] is None and rec["h"] is None and rec["gap"] == pytest.approx(0.25)


def test_apply_size_procedure_accepts_both_analogue_dicts(h27):
    C, _ = _cloud_with_core(seed=3)
    for size in (h27.analogue_size_result(_four(h27), 10.5),
                 h27.analogue_size_result([_geom(h27, "E8.25_late", 8.25, 150.0, 350.0)], 8.5)):
        targets = h27.forecast_targets(size)
        for form in h27.FORMS:
            out, rec = t2_build.apply_size_procedure(C, size, h27.HYP[form], z_mode="raw", expect_board="heart")
            assert rec["z_mode_used"] == "raw" and rec["target_rms_full"] == targets[form]
            assert abs(rec["log_rms_after_vs_target"]) <= 1e-12
            assert not rec["warnings"]
        assert json.loads(json.dumps(size, default=float))["board"] == "heart"  # the part can hold it


def test_each_form_takes_the_builders_field(h27):
    """A to scope H1 (the copy s1 path, `"s1": "H1"` in build_heart), B to the field build_test names for the half
    trim, C to the field it names for the full trim; each arm lands on its target."""
    sys.path.insert(0, str(REPO_ROOT / "experiments" / "builds"))
    import build_test

    opts = build_test.HEART_SCOPE_OPTIONS
    assert h27.HYP["half_iso_half"] == opts["s2_trim"]["hypothesis"] == "H2"
    assert h27.HYP["core_iso_half"] == opts["core_trim"]["hypothesis"] == t2_build.SIZE_CORE_SCOPE
    assert opts["s2_trim"]["frac"] == h27.h26.TRIM["half_iso_half"] and opts["core_trim"]["frac"] == \
        h27.h26.TRIM["core_iso_half"]
    assert '"s1": "H1"' in inspect.getsource(build_test.build_heart) and h27.HYP["full_iso_half"] == "H1"
    size = h27.analogue_size_result(_four(h27), 10.5)
    targets = h27.forecast_targets(size)

    def field(path):
        v = size
        for k in path.split("."):
            v = v[k]
        return v

    assert targets["half_iso_half"] == field(opts["s2_trim"]["size_source"])
    assert targets["core_iso_half"] == field(opts["core_trim"]["size_source"])
    assert targets["full_iso_half"] == field("z_modes.raw.scope.H1.pred_rms_full")
    C, _ = _cloud_with_core(seed=4)
    X = np.zeros((len(C), 2))
    for form in h27.FORMS:
        _, _, rec = h27.size_arm(X, C, form, "f", {"even": 1.0, "in_plane": 1.0}, size, {}, 3)
        assert rec["rms_after"] == pytest.approx(targets[form], rel=1e-9)


def test_perfect_twin_targets_and_c_offset_is_reference_only(h27):
    TC, tcore = _cloud_with_core(seed=5)
    PC, pcore = _cloud_with_core(n=2500, seed=6)
    pt = h27.perfect_targets(TC, tcore)
    assert pt["full_iso_half"] == pt["half_iso_half"] == t2_build.rms_radius(TC)
    assert pt["core_iso_half"] == t2_build.rms_radius(TC[tcore])
    pm = h27.h26.pool_masks(PC, pcore)
    pool_rms = {f: t2_build.rms_radius(PC[pm[f]]) for f in h27.FORMS}
    ref_core = t2_build.rms_radius(PC[pcore])
    offs = []
    for scale in (1.0, 1.7):  # another truth: C's offset does not move, A's does
        T = TC * scale
        f = h27.h26.oracle_factors(T, tcore, PC, pcore)
        o = h27.convention_offsets(h27.perfect_targets(T, tcore), f["even"], pool_rms, ref_core)
        assert o["core_iso_half"] == pytest.approx(o["C_parts"]["log_rms_ref_core"] - o["C_parts"]["log_rms_c_pool"],
                                                   abs=1e-12)
        assert o["full_iso_half"] == pytest.approx(np.log(t2_build.rms_radius(T) / t2_build.rms_radius(T[tcore]))
                                                   - np.log(pool_rms["full_iso_half"] / ref_core), abs=1e-12)
        offs.append(o)
    assert offs[0]["core_iso_half"] == pytest.approx(offs[1]["core_iso_half"], abs=1e-12)
    C = np.asarray(PC, float)
    _, _, rec = h27.size_arm(np.zeros((len(C), 1)), C, "full_iso_half", "p", {"even": 1.0}, {}, pt, 2)
    assert abs(rec["log_rms_after_vs_target"]) <= 1e-12 and rec["target_rms_full"] == pt["full_iso_half"]


def test_the_size_files_days_are_run_p3s(h27):
    text = (REPO_ROOT / "experiments" / "coords" / "run_p3.sh").read_text()
    files = {}
    for var, f in re.findall(r"\b(E\d+)=(?:\$\{E\d+:-)?\$T2/([\w.]+\.h5ad)", text):
        files[var] = f"T2/{f}"  # the last assignment wins (the release-day branch follows the rehearsal's)
    line = next(ln for ln in text.splitlines() if "size_procedure.py --board heart" in ln)
    line += next(ln for ln in text.splitlines()[text.splitlines().index(line) + 1:])
    days = {files[v]: float(d) for v, d in re.findall(r'--stage "\$(E\d+)" ([\d.]+)', line)}
    assert days.pop("T2/E10.5.h5ad") == 10.5
    assert days == h27.SIZE_FILE_DAYS


def test_the_target_never_reaches_the_size_code(h27, monkeypatch, tmp_path):
    seen = []

    def spy(path, t, name=None, **kw):
        seen.append((Path(path).relative_to(tmp_path).as_posix(), t))
        return _geom(h27, Path(path).stem, t, 100.0 + t, 200.0 + t), {"path": str(path)}

    monkeypatch.setattr(h27.sp, "measure_stage", spy)
    for key in h27.TASK_KEYS:
        seen.clear()
        task = h27.payoff.heart_task(key)
        h27.input_geometries(task, tmp_path)
        assert task.target not in [p for p, _ in seen]
        assert sorted(p for p, _ in seen) == sorted(set(task.inputs.values()))
        assert all(t == h27.SIZE_FILE_DAYS[p] for p, t in seen)
        if key == "T2h/<=E8.75->E10.5":
            assert "T2/E9.5.h5ad" not in [p for p, _ in seen]
    src = inspect.getsource(h27.run_transition)
    assert src.index("input_geometries(t, data_dir)") < src.index("read_stage(data_dir / t.target")
    assert "t.target" not in inspect.getsource(h27.input_geometries).replace("task.target in files", "")


# ------------------------------------------------------------------------------------------------ end to end, synthetic
@pytest.fixture(scope="module")
def synth():
    spec = importlib.util.spec_from_file_location("b0_trend_tests_for_h27", REPO_ROOT / "tests" / "test_b0_trend.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    inputs = {"ref": mod._stage(1600, seed=2, shift=0.0, tiny=60),
              "start": mod._stage(1100, seed=3, shift=0.8, tiny=10),
              "earlier": mod._stage(900, seed=5, shift=1.6, tiny=10)}
    for a in inputs.values():
        a.obs["core"] = np.asarray(a.obs["merged"]).astype(str) == "CM"
    return mod, inputs


def test_o_p_f_share_rows_and_expression(h27, synth):
    """One draw, one half: o is row H26's arm exactly; p and f are the same cells and expression at another size."""
    mod, inputs = synth
    sp = h27.sp
    geoms = [sp.geometry_from_arrays(n, t, np.asarray(inputs[r].obsm["spatial_3D"], float),
                                     np.asarray(inputs[r].obs["core"], bool))
             for r, n, t in (("earlier", "E8.5", 8.5), ("start", "E8.75", 8.75), ("ref", "E9.5", 9.5))]
    size = h27.analogue_size_result(geoms, 10.5)
    ft = h27.forecast_targets(size)
    pool = inputs["ref"]
    PC = np.asarray(pool.obsm["spatial_3D"], dtype=float)
    pm = h27.h26.pool_masks(PC, np.asarray(pool.obs["core"], bool))
    factors = {"even": 1.3, "in_plane": 1.4}
    pt = {"full_iso_half": 9.0, "half_iso_half": 9.0, "core_iso_half": 4.0}
    seed = harness.BUILD_SEED_OFFSET + 3
    for form in h27.FORMS:
        p = pool if h27.h26.TRIM[form] is None else pool[pm[form]].copy()
        ctx = mod._ctx({**inputs, "ref": p}, seed=seed)
        built = h27.build_forms(form, ctx, factors, size, pt)
        theirs = h27.h26.build_arm(form, mod._ctx({**inputs, "ref": p}, seed=seed), factors, "copy")
        o, p_, f = (built[z][0] for z in h27.SIZINGS)
        np.testing.assert_array_equal(o.X, theirs.X)
        np.testing.assert_array_equal(o.coords, theirs.coords)
        for other in (p_, f):
            np.testing.assert_array_equal(other.X, o.X)
            assert other.n == o.n == mod.BOARD.max_cells // 2
            a = other.coords.astype(float) - other.coords.astype(float).mean(0)
            b = o.coords.astype(float) - o.coords.astype(float).mean(0)
            k = np.sqrt((a ** 2).sum() / (b ** 2).sum())
            np.testing.assert_allclose(a, b * k, rtol=0, atol=2e-3 * np.abs(a).max())
        assert built["f"][1]["target_rms_full"] == ft[form]
        assert abs(built["f"][1]["log_rms_after_vs_target"]) <= 1e-9
        assert built["p"][1]["target_rms_full"] == pt[form]
        assert abs(built["p"][1]["log_rms_after_vs_target"]) <= 1e-9
        assert built["o"][1]["factor"] == 1.3


def test_v7_passes_the_untrimmed_reference_and_fails_a_trimmed_pool(h27):
    C, core = _cloud_with_core(seed=11)
    g = h27.sp.geometry_from_arrays("E9.5", 9.5, C, core)
    size = h27.analogue_size_result([_geom(h27, "E8.75", 8.75, 100.0, 300.0), g], 10.5)
    assert h27.v7_crosscheck(C, size)["ok"]
    trimmed = C[t2_build.radial_trim_mask(C, core, 0.5)]
    v = h27.v7_crosscheck(trimmed, size)
    assert not v["ok"] and abs(v["log_diff"]) > t2_build.SIZE_CROSSCHECK_LOG_TOL
    assert not h27.v7_crosscheck(C * (1 + 1e-5), size)["ok"]  # beyond float noise


# ------------------------------------------------------------------------------------------------ hand-made tables
def _committed_geometry():
    with open(COORDS / "stage_geometry.csv", newline="") as f:
        return {r["name"]: r for r in csv.DictReader(f)}


def _stage_dict(h27, row: dict) -> dict:
    fields = h27.sp.StageGeometry.__dataclass_fields__
    d = {k: (row[k] if k == "name" else float(row[k])) for k in fields}
    d["log_scope"] = float(row["log_scope"])
    return d


FAKE_E85 = {"name": "E8.5", "t": 8.5, "n_cells": 1000, "n_core": 100, "n_levels": 40, "z_step": 10.0,
            "z_range": 390.0, "rms_full": 300.0, "rms_core": 150.0, "rms_xy_full": 250.0, "sd_z_full": 150.0,
            "rms_xy_core": 120.0, "sd_z_core": 80.0, "core_thin_axis_cos_z": 0.5, "log_scope": float(np.log(2.0))}
STAGES = {"T2h/<=E9.5->E10.5": ("E8.25_late", "E8.5", "E8.75", "E9.5"),
          "T2h/<=E8.75->E10.5": ("E8.25_late", "E8.5", "E8.75"),
          "T2h/E8.25_late->E8.5": ("E8.25_late",)}


def _stub_size(h27, key: str) -> dict:
    """A stub of the analogue result: committed pre-release geometries (and a made-up E8.5), the rules a valid run
    records and made-up forecast numbers. No procedure runs on them."""
    geo = _committed_geometry()
    stages = [FAKE_E85 if n == "E8.5" else _stage_dict(h27, geo[n]) for n in STAGES[key]]
    chosen = ["copy_last"] if key == h27.ONE_STAGE_KEY else list(h27.HEART_RULES)
    r = {"valid": True, "z_mode": "raw", "reference": stages[-1]["name"], "stages": stages,
         "growth": {"core": {"chosen": chosen}}, "pred_rms_core": 100.0,
         "scope": {"H0": {"pred_rms_full": 130.0}, "H1": {"pred_rms_full": 120.0}, "H2": {"pred_rms_full": 110.0}},
         "min_start_t": None, "tie_tol": None}
    return {"board": "heart", "valid": True, "preregistered": r, "z_modes": {"raw": r},
            "provenance": {"analogue": h27.ANALOGUE_COPY if len(stages) == 1 else h27.ANALOGUE_RATES}}


def _info(h26):
    return {"truth_info": {q: {"cells": 100 + i, "core_share": 0.2 + 0.1 * i, "xy_scope": 0.1}
                           for i, q in enumerate(h26.TRUTHS)},
            "pred_info": {"full_iso_half": {"pool_cells": 900, "core_share": 0.23, "xy_scope": 0.2},
                          "half_iso_half": {"pool_cells": 600, "core_share": 0.3, "xy_scope": 0.1},
                          "core_iso_half": {"pool_cells": 300, "core_share": 1.0, "xy_scope": 0.0},
                          "s1form_half": {"pool_cells": 900, "core_share": 0.23, "xy_scope": 0.2}},
            "oracle_factors": {"label": "oracle size", "even": 1.1, "in_plane": 1.2}}


def _skills(**over):
    s = {m: 0.5 for m in ("de_score", "de_direction", "mmd_u", "variogram", "d2_shape", "occupancy_dice",
                          "scale_log_ratio", "neighborhood_mmd")}
    return {**s, **over}


def _groups():
    return {g: 0.125 for g in ("expression_change", "cell_state", "shape_growth", "local_organisation")}


def _part(h27, key, seeds, score=None, skills=None, commit="abc123", clean=True):
    """A hand-made part. score: {truth: {arm: value or per-seed list}}; default 50. skills: {(truth, arm): {...}}."""
    score, skills = score or {}, skills or {}
    info = _info(h27.h26)
    size = _stub_size(h27, key)
    ft = h27.forecast_targets(size)
    pt = {"full_iso_half": 200.0, "half_iso_half": 200.0, "core_iso_half": 90.0}
    rows = []
    for s in seeds:
        sc = {q: {a: (lambda v: v[s] if isinstance(v, list) else v)(score.get(q, {}).get(a, 50.0)) for a in h27.ALL}
              for q in h27.TRUTHS}
        sk = {q: {a: _skills(**skills.get((q, a), {})) for a in h27.ALL} for q in h27.TRUTHS}
        gr = {q: {a: _groups() for a in h27.ALL} for q in h27.TRUTHS}
        recs = {}
        for f in h27.FORMS:
            recs[h27.arm_name(f, "o")] = {"sizing": "o", "factor": 1.1}
            recs[h27.arm_name(f, "p")] = {"sizing": "p", "target_rms_full": pt[f], "log_rms_after_vs_target": 0.0}
            recs[h27.arm_name(f, "f")] = {"sizing": "f", "target_rms_full": ft[f], "log_rms_after_vs_target": 0.0}
        rows.append({"seed": s, "n_cells": {a: 500 for a in h27.ARM_NAMES}, "size_records": recs, "scores": sc,
                     "skills": sk, "groups": gr})
    fc = {"inputs": list(STAGES[key]), "chosen": size["z_modes"]["raw"]["growth"]["core"]["chosen"],
          "interval": list(STAGES[key][-2:]) if len(STAGES[key]) > 1 else None, "gap": 1.0,
          "h": 1.33 if len(STAGES[key]) > 1 else None}
    return {"axis": h27.AXIS, "task": key, "rows": rows, "size_result": size, "forecast": fc, "forecast_targets": ft,
            "perfect_targets": pt, "truth_rms": {"full": 200.0, "core": 90.0},
            "offsets": {"full_iso_half": 0.1, "half_iso_half": 0.2, "core_iso_half": -0.3,
                        "C_parts": {"log_rms_ref_core": 5.0, "log_rms_c_pool": 5.3}},
            "v7": {"ok": True, "log_diff": 0.0}, "departures": h27.DEPARTURES,
            "provenance": {"commit": commit, "tree_clean": clean}, **info}


def _parts(h27, score=None, **kw):
    return {k: [_part(h27, k, (2 * j, 2 * j + 1), score=(score or {}).get(k), **kw) for j in range(5)]
            for k in h27.TASK_KEYS}


def _h26_parts(h27, score=None):
    """Row H26's committed parts, hand-made: its arms' scores (default 50) and the same records."""
    h26 = h27.h26
    out = {}
    for key in h27.TASK_KEYS:
        ps = []
        for j in range(5):
            rows = []
            for s in (2 * j, 2 * j + 1):
                sc = {q: {a: 50.0 for a in (*h26.ALL,)} for q in h26.TRUTHS}
                for (k, q, a), v in (score or {}).items():
                    if k == key:
                        sc[q][a] = v[s] if isinstance(v, list) else v
                rows.append({"seed": s, "scores": sc})
            ps.append({"axis": h26.AXIS, "task": key, "rows": rows, **_info(h26)})
        out[key] = ps
    return out


@pytest.fixture(scope="module")
def committed(h27):
    return h27.load_committed()


def _read(h27, committed, parts=None, h26p=None, commit="abc"):
    return h27.reading(parts if parts is not None else _parts(h27), h26p or _h26_parts(h27), committed, commit)


def test_a_clean_hand_made_run_is_read_and_prices_nothing_more(h27, committed):
    res = _read(h27, committed)
    assert not res["void"], res["why"]
    assert res["decides"].startswith("nothing") and "pick" not in {k.lower() for k in res}
    assert res["h26_quoted"] == {"own": {"B-A": -1.00, "C-A": -3.86, "C-B": -2.85},
                                 "CLEAN": {"B-A": 5.78, "C-A": 11.79, "C-B": 6.01}}
    assert res["guards"]["V6"]["geometry"] == 0.0 and res["guards"]["V6"]["procedure"] <= 1e-9
    assert set(res["labels"]) == {"f", "p", "o"}
    g = res["groups"][h27.TASK_KEYS[0]]["own"]
    assert set(g["arms"]) == set(h27.ARM_NAMES) and set(g["size_terms"]) == {"A", "B", "C"}


def test_void_on_missing_repeated_mixed_or_dirty_parts(h27, committed):
    p = _parts(h27)
    del p[h27.TASK_KEYS[1]]
    assert _read(h27, committed, p)["void"]
    p = _parts(h27)
    p[h27.TASK_KEYS[0]] = p[h27.TASK_KEYS[0]][:4]
    assert _read(h27, committed, p)["void"]
    p = _parts(h27)
    p[h27.TASK_KEYS[0]].append(_part(h27, h27.TASK_KEYS[0], (9,)))
    assert any("more than one part" in w for w in _read(h27, committed, p)["why"])
    p = _parts(h27)
    p[h27.TASK_KEYS[2]][0]["provenance"]["commit"] = "def456"
    assert any("different commits" in w for w in _read(h27, committed, p)["why"])
    p = _parts(h27)
    p[h27.TASK_KEYS[2]][1]["provenance"]["tree_clean"] = False
    assert any("dirty" in w for w in _read(h27, committed, p)["why"])
    assert any("row's commit" in w for w in _read(h27, committed, _parts(h27), commit="fff")["why"])
    assert _read(h27, committed, _parts(h27), commit="")["void"]
    p = _parts(h27)
    del p[h27.TASK_KEYS[0]][0]["rows"][0]["scores"]["thin_0"]["core_iso_half:f"]
    assert any(w.startswith("V1") for w in _read(h27, committed, p)["why"])
    p = _parts(h27)
    del p[h27.TASK_KEYS[0]][0]["rows"][0]["size_records"]["core_iso_half:p"]
    assert any(w.startswith("V1") for w in _read(h27, committed, p)["why"])
    p = _parts(h27)
    del p[h27.TASK_KEYS[0]][0]["v7"]  # a malformed part voids, it does not crash
    assert _read(h27, committed, p)["void"]


def test_void_when_a_truth_a_pool_or_a_factor_differs_from_h26(h27, committed):
    for change in (lambda p: p["truth_info"]["thin_0"].update(cells=99),
                   lambda p: p["truth_info"]["own"].update(core_share=0.21),
                   lambda p: p["pred_info"]["core_iso_half"].update(pool_cells=301),
                   lambda p: p["oracle_factors"].update(even=1.1 * (1 + 1e-11)),
                   lambda p: p["oracle_factors"].update(in_plane=1.3)):
        p = _parts(h27)
        change(p[h27.TASK_KEYS[1]][2])
        assert any(w.startswith("V3") for w in _read(h27, committed, p)["why"])
    p = _parts(h27)
    p[h27.TASK_KEYS[1]][2]["oracle_factors"]["even"] = 1.1 * (1 + 1e-13)  # within 1e-12
    assert not _read(h27, committed, p)["void"]


def test_void_when_a_twin_departs_from_h26(h27, committed):
    key, arm = h27.TASK_KEYS[2], "half_iso_half"
    for mean_off, seed_off, void in ((0.09, 0.0, False), (0.11, 0.0, True), (0.0, 0.49, False),
                                     (0.0, 0.51, True)):
        theirs = [50.0 + mean_off + (seed_off if s == 3 else 0.0) - (seed_off / 9 if s != 3 else 0.0)
                  for s in SEEDS]
        res = _read(h27, committed, h26p=_h26_parts(h27, {(key, "radial_0", arm): theirs}))
        assert res["void"] is void, (mean_off, seed_off, res["why"])
        if not void:
            assert res["guards"]["V4"][key]["seed"] <= 0.5
    t = h27.text(_read(h27, committed, h26p=_h26_parts(h27, {(key, "own", arm): 50.05})))
    assert "largest distance to row H26's committed scores 0.050 (ten-seed mean) and 0.050 (one seed)" in t


def test_void_when_the_rescale_identity_breaks(h27, committed):
    key = h27.TASK_KEYS[0]
    for arm, over, void in (("full_iso_half:f", {"mmd_u": 0.5 + 0.9e-9}, False),
                            ("full_iso_half:f", {"mmd_u": 0.5 + 2e-9}, True),
                            ("core_iso_half:p", {"de_score": 0.5 - 2e-9}, True),
                            ("half_iso_half:p", {"neighborhood_mmd": 0.5 + 0.9e-3}, False),
                            ("half_iso_half:f", {"neighborhood_mmd": 0.5 + 1.1e-3}, True),
                            ("half_iso_half:f", {"variogram": 0.9, "d2_shape": 0.1}, False)):
        p = _parts(h27)
        p[key][1] = _part(h27, key, (2, 3), skills={("thin_0.5", arm): over})
        assert _read(h27, committed, p)["void"] is void, (arm, over)
    p = _parts(h27)
    p[key][0]["rows"][1]["n_cells"]["core_iso_half:f"] = 499
    assert any(w.startswith("V5") for w in _read(h27, committed, p)["why"])


def test_void_when_rules_chosen_are_not_all_three_or_size_does_not_reproduce(h27, committed):
    def check(edit, prefix):
        p = _parts(h27)
        edit(p)
        why = _read(h27, committed, p)["why"]
        assert any(w.startswith(prefix) for w in why), why

    k0, k1, k2 = h27.TASK_KEYS
    check(lambda p: p[k0][0]["size_result"]["z_modes"]["raw"]["growth"]["core"].update(chosen=["half_rate"]), "V8")
    check(lambda p: p[k1][0]["size_result"]["z_modes"]["raw"]["growth"]["core"].update(
        chosen=["power_t", "latest_rate", "half_rate"]), "V8")
    check(lambda p: p[k2][0]["size_result"]["z_modes"]["raw"]["growth"]["core"].update(chosen=["latest_rate"]), "V8")
    check(lambda p: p[k0][0]["size_result"].update(valid=False), "V8")
    check(lambda p: p[k0][0]["forecast_targets"].update(core_iso_half=float("nan")), "V8")
    check(lambda p: p[k0][0]["perfect_targets"].update(half_iso_half=-1.0), "V8")
    check(lambda p: p[k0][0]["rows"][0]["size_records"]["full_iso_half:f"].update(log_rms_after_vs_target=2e-9), "V8")
    check(lambda p: p[k0][0]["rows"][0]["size_records"]["full_iso_half:p"].update(log_rms_after_vs_target=2e-9), "V8")
    check(lambda p: p[k0][0]["rows"][0]["size_records"]["half_iso_half:f"].update(target_rms_full=111.0), "V8")
    check(lambda p: p[k1][3].update(v7={"ok": False, "log_diff": 1e-3}), "V7")
    check(lambda p: p[k1][4]["forecast_targets"].update(full_iso_half=121.0), "V8")  # one forecast per transition
    # V6: a geometry off by more than 1e-9 relative, or the procedure's pseudo-target off the committed row
    check(lambda p: p[k0][0]["size_result"]["z_modes"]["raw"]["stages"][2].update(
        rms_core=154.4831472408458 * (1 + 2e-9)), "V6")
    off = {"geometry": committed["geometry"], "modes_raw": {**committed["modes_raw"],
                                                           "full_H1": committed["modes_raw"]["full_H1"] * 1.001}}
    assert any(w.startswith("V6") for w in h27.reading(_parts(h27), _h26_parts(h27), off, "abc")["why"])
    p = _parts(h27)
    p[k0][0]["size_result"]["z_modes"]["raw"]["stages"][2]["rms_core"] *= 1 + 0.2e-9  # inside 1e-9
    assert not _read(h27, committed, p)["void"]


def test_the_void_reading_prints_only_the_void_line(h27, committed):
    p = _parts(h27)
    del p[h27.TASK_KEYS[1]]
    res = _read(h27, committed, p)
    t = h27.text(res)
    assert t.count("\n") == 1 and t.startswith("VOID: V1: missing transition T2h/E8.25_late->E8.5")
    assert t.endswith(". No forecast-size figure exists under row H27 and there is no rerun under it; row H26's table "
                      "stands alone, at its core-ratio oracle size.\n")
    assert "tables" not in res and "ranges" not in res


# ------------------------------------------------------------------------------------------------ the fixed words
def _scored(h27):
    """own truth on every transition: A at f 50, B at f 52, C at f alternating 49/51 by seed; A at p 49."""
    alt = [49.0 if s % 2 == 0 else 51.0 for s in SEEDS]
    sc = {"own": {"half_iso_half:f": 52.0, "core_iso_half:f": alt, "full_iso_half:p": 49.0}}
    return _parts(h27, {k: sc for k in h27.TASK_KEYS})


def test_the_text_is_the_fixed_words_with_labels_on_every_table_line_and_ranges_first(h27, committed):
    res = _read(h27, committed, _scored(h27))
    assert not res["void"], res["why"]
    lines = h27.text(res).splitlines()
    assert lines[0] == h27.L0 and lines[1] == h27.L1.format(commit="abc123")
    assert lines[2] == (
        "T2h/<=E9.5->E10.5: inputs E8.25_late, E8.5, E8.75, E9.5; rules latest_rate, half_rate, power_t; rate interval "
        "E8.75 to E9.5 (test: E9.5 to E10.5); horizon 1.00 days, h 1.33 (test: 2 days, h 2). Forecast log error on "
        "this released target, in the procedure's units: core +0.105, full H1 -0.511, full H2 -0.598. Convention "
        "offset log(p/o): A +0.100, B +0.200, C -0.300 (C's from the reference alone: log rms(ref core) 5.000 minus "
        "log rms(C pool) 5.300).")
    assert lines[3].startswith("T2h/E8.25_late->E8.5: inputs E8.25_late; rules copy_last; rate interval none (one "
                               "input stage) (test: E9.5 to E10.5); horizon 1.00 days, h none (test: 2 days, h 2).")
    s0 = lines.index(h27.S0_HEAD)
    s1, s2, s3 = (lines.index(h) for h in (h27.S1_HEAD, h27.S2_HEAD, h27.S3_HEAD))
    assert lines[s0 - 1].startswith("ORACLE TWINS") and s0 < s1 < s2 < s3
    assert lines[s0 + 1] == ("own: over the three transitions, at the forecast size B-A from +2.00 to +2.00, C-A from "
                             "+0.00 to +0.00, C-B from -2.00 to -2.00; forecast error (f-p) A from +1.00 to +1.00, B "
                             "from +2.00 to +2.00, C from +0.00 to +0.00.")
    assert lines[s0 + 2].startswith("CLEAN: over the three transitions, at the forecast size B-A from +0.00 to +0.00")
    assert lines[s0 + 3] == h27.S0_TAIL
    n_units = len(h27.UNITS)
    for a, b in ((s1, s2), (s2, s3), (s3, len(lines) - 2)):
        body = lines[a + 1:b]
        assert len(body) == len(h27.WORLDS) * n_units
        assert all(ln.startswith(h27.TAG + " ") for ln in body)
        units = [ln.split(":")[0] for ln in body[:n_units]]
        assert units[-1].endswith("mean of folds (the two folds weighed equally)")
        assert units[-2].endswith("fold E8.5 (copy_last analogue, no growth step)")
    assert lines[s1 + 1] == ("[practice analogue, released stages] own T2h/<=E9.5->E10.5: A 50.00 B 52.00 C 50.00 | "
                             "B-A +2.00 [+2.00, +2.00] 10/10; C-A +0.00 [-1.00, +1.00] 5/10; C-B -2.00 [-3.00, -1.00] "
                             "0/10")
    assert lines[s2 + 1] == ("[practice analogue, released stages] own T2h/<=E9.5->E10.5: A f-p +1.00 [+1.00, +1.00] "
                             "10/10, p-o -1.00, f-o +0.00; B f-p +2.00 [+2.00, +2.00] 10/10, p-o +0.00, f-o +2.00; C "
                             "f-p +0.00 [-1.00, +1.00] 5/10, p-o +0.00, f-o +0.00")
    assert lines[s3 + 1] == ("[practice analogue, released stages] own T2h/<=E9.5->E10.5: f-o B-A +2.00, C-A +0.00, "
                             "C-B -2.00; f-p B-A +1.00, C-A -1.00, C-B -2.00")
    assert lines[-2] == h27.H26_LINE and lines[-1] == h27.CLOSING
    worlds = [ln.split(" ")[4] for ln in lines[s1 + 1:s2:n_units]]
    assert worlds == ["own", "CLEAN", "truth", "truth", "truth", "truth", "truth"]
    t = "\n".join(lines).lower()
    assert "recommend" not in t and "we should" not in t and "prefer " not in t


def test_folds_merge_the_two_e10_5_transitions_seed_by_seed_and_mean_of_folds_prints_last(h27, committed):
    k95, k85, k875 = h27.TASK_KEYS
    up = [50.0 + s for s in SEEDS]
    down = [50.0 - s for s in SEEDS]
    p = _parts(h27, {k95: {"own": {"half_iso_half:f": up}}, k875: {"own": {"half_iso_half:f": down}},
                     k85: {"own": {"half_iso_half:f": 54.0}}})
    res = _read(h27, committed, p)
    c = res["tables"]["own"]
    assert c["fold E10.5"]["diff"]["f"]["B-A"] == {"mean": 0.0, "lo": 0.0, "hi": 0.0, "ahead": 0, "behind": 0}
    assert c[k95]["diff"]["f"]["B-A"]["mean"] == pytest.approx(4.5)
    assert c["fold E8.5"]["diff"]["f"]["B-A"]["mean"] == pytest.approx(4.0)
    assert c["mean of folds"]["diff"]["f"]["B-A"]["mean"] == pytest.approx(2.0)
    assert list(h27.UNITS)[-1] == "mean of folds" and list(h27.UNITS)[:3] == list(h27.TASK_KEYS)
    assert res["ranges"]["own"]["diff_f"]["B-A"] == [pytest.approx(-4.5), pytest.approx(4.5)]


def test_the_run_logs_print_no_score(h27):
    """The run prints only "<task> seed <s> (<sec> s) done" and "wrote <path>"; every other output is a refusal."""
    for fn in (h27.run_transition, h27.cmd_run, h27.build_forms, h27.size_arm, h27.input_geometries,
               h27.analogue_size_result):
        tree = ast.parse(inspect.getsource(fn).lstrip())
        prints = [ast.unparse(n) for n in ast.walk(tree) if isinstance(n, ast.Call) and getattr(n.func, "id", "") ==
                  "print"]
        for p in prints:
            assert p in ('print(f"{key} seed {s} ({row[\'seconds\']} s) done", flush=True)',
                         "print(f'wrote {path}', flush=True)"), p
        assert len(prints) == (2 if fn is h27.run_transition else 0), fn


def test_the_script_names_no_test_stage_file_and_reads_no_test_score(h27):
    src = SCRIPT.read_text()
    for bad in ("results/p3", "/p3/", "size_heart_E12", "submissions", "leaderboard", "registry", "build_test",
                "boards/", "upload", "12.5"):
        assert bad not in src, bad
    assert "SIZES" not in src and "size_or_hold(" not in src


def test_each_run_refuses_before_anything_loads(h27, monkeypatch, tmp_path):
    def boom(*a, **k):
        raise AssertionError("something loaded")

    monkeypatch.setattr(h27.payoff, "heart_task", boom)
    monkeypatch.setattr(h27.sp, "measure_stage", boom)
    monkeypatch.setattr(h27, "tries_rows", lambda path=h27.TRIES: {"H26"})
    with pytest.raises(SystemExit, match="no row"):
        h27.run_transition(h27.TASK_KEYS[0], (0,), tmp_path, tmp_path)
    monkeypatch.setattr(h27, "tries_rows", lambda path=h27.TRIES: {"H26", "H27"})
    monkeypatch.setattr(h27.provenance, "git_state", lambda repo, ignore=(): ("c" * 40, ["src/x.py"]))
    with pytest.raises(SystemExit, match="not clean"):
        h27.run_transition(h27.TASK_KEYS[0], (0,), tmp_path, tmp_path)
    monkeypatch.setattr(h27.provenance, "git_state", lambda repo, ignore=(): ("unknown", None))
    with pytest.raises(SystemExit, match="not clean"):
        h27.run_transition(h27.TASK_KEYS[0], (0,), tmp_path, tmp_path)
    monkeypatch.setattr(h27.provenance, "git_state", lambda repo, ignore=(): ("c" * 40, []))
    (tmp_path / h27.PART).write_text("{}")
    with pytest.raises(SystemExit, match="exists"):
        h27.run_transition(h27.TASK_KEYS[0], (0,), tmp_path, tmp_path)
    (tmp_path / h27.PART).unlink()
    with pytest.raises(AssertionError, match="something loaded"):  # past every refusal, the first load
        h27.run_transition(h27.TASK_KEYS[0], (0,), tmp_path, tmp_path)


def test_the_row_id_reader_finds_ledger_rows(h27, tmp_path):
    f = tmp_path / "TRIES.md"
    f.write_text("| # | Board |\n|---|---|\n| H26 (one row) | heart |\n| H27 (one row; eighteen tries) | heart |\n"
                 "text | H28 |\n")
    assert h27.tries_rows(f) == {"H26", "H27"}
    assert "H26" in h27.tries_rows()


# ------------------------------------------------------------------------------------------------ the record
def test_the_guards_dry_run_on_committed_results(h27, committed):
    """Rule 9 J7: the guards on committed data other than the pair the row reads.
    V4's tolerance: row H25's copy against row H11's committed copy, distance 0.000 on both E10.5 backtests.
    V6: the procedure's committed pseudo-target rows reproduce its committed mode row (relative 0), H1 and H2 recomputed
    from stage_geometry.csv equal size_heart_scope.csv, and the procedure on stage_geometry.csv's E8.25_late and E8.75
    reproduces size_heart_modes.csv's raw no-cut-off row at the pre-registered tie_tol and at tie_tol = inf.
    V7: row H25's reference_rms_full (an independent load) equals stage_geometry.csv's rms_full within 1e-9 relative."""
    g25 = json.loads((RESULTS / "h25_heart_forms" / "reading.json").read_text())["guard"]
    assert all(v == {"mean": 0.0, "largest_seed": 0.0} for v in g25.values()) and len(g25) == 2
    sp = h27.sp
    with open(COORDS / "size_heart_pseudo.csv", newline="") as f:
        pseudo = [r for r in csv.DictReader(f) if r["variant"] == "raw" and r["series"] == "core"
                  and r["target"] == "E9.5" and r["train"] == "E8.25_late,E8.75" and r["rule"] in h27.HEART_RULES]
    assert len(pseudo) == 3
    core = np.exp(np.mean([float(r["pred_log"]) for r in pseudo]))
    assert core == pytest.approx(committed["modes_raw"]["pred_rms_core"], rel=1e-14)
    geo = _committed_geometry()
    s = [float(geo[n]["log_scope"]) for n in ("E8.25_late", "E8.75")]
    with open(COORDS / "size_heart_scope.csv", newline="") as f:
        scope = {r["hypothesis"]: float(r["pred_log_scope"]) for r in csv.DictReader(f)}
    for h in ("H1", "H2"):
        assert sp.scope_log_ratio(s, [8.25, 8.75], 9.5, h) == pytest.approx(scope[h], rel=1e-12)
    geoms = [h27._geometry(_stage_dict(h27, geo[n])) for n in ("E8.25_late", "E8.75")]
    for tie in (None, np.inf):
        r = sp.heart_size_from_geometry(geoms, 9.5, z_mode="raw", min_start_t=-np.inf, tie_tol=tie)
        assert r["pred_rms_core"] == pytest.approx(committed["modes_raw"]["pred_rms_core"], rel=1e-12)
        assert r["scope"]["H1"]["pred_rms_full"] == pytest.approx(committed["modes_raw"]["full_H1"], rel=1e-12)
    rel = {}
    for f, n in (("T2h_E9.5_E10.5.json", "E9.5"), ("T2h_E8.75_E10.5.json", "E8.75")):
        t = json.loads((RESULTS / "h25_heart_forms" / "tasks" / f).read_text())
        rel[n] = abs(t["oracle_size"]["reference_rms_full"] / float(geo[n]["rms_full"]) - 1)
    print(f"\nH27 V7 dry run: relative distance of row H25's reference RMS to stage_geometry.csv {rel}")
    assert max(rel.values()) <= 1e-9


def test_the_h26_parts_agree_on_the_records_v3_checks(h27):
    parts = h27.load_h26()
    assert set(parts) == set(h27.TASK_KEYS)
    for ps in parts.values():
        assert len(ps) == 5 and sorted(r["seed"] for p in ps for r in p["rows"]) == list(SEEDS)
        for p in ps[1:]:
            assert p["truth_info"] == ps[0]["truth_info"] and p["pred_info"] == ps[0]["pred_info"]
            assert p["oracle_factors"] == ps[0]["oracle_factors"]
            assert p["provenance"]["commit"].startswith("b2e84dd") and p["provenance"]["tree_clean"] is True


def test_the_noise_floor_of_rule_9_j5(h27):
    """J5 from row H26's committed parts: A minus the floor with the truth as released, seed SE per fold and the
    between-fold SD of the fold means."""
    parts = h27.load_h26()
    rows = {k: h27.h26.merge_rows(ps)[0] for k, ps in parts.items()}

    def d(key):
        return np.array([r["scores"]["own"]["full_iso_half"] - r["scores"]["own"]["floor"] for r in rows[key]])

    e105 = (d("T2h/<=E9.5->E10.5") + d("T2h/<=E8.75->E10.5")) / 2
    e85 = d("T2h/E8.25_late->E8.5")
    se = [x.std(ddof=1) / np.sqrt(len(x)) for x in (e105, e85)]
    between = np.std([e105.mean(), e85.mean()], ddof=1)
    print(f"\nH27 J5: seed SE {se[0]:.3f} (fold E10.5), {se[1]:.3f} (fold E8.5); between folds {between:.3f}")
    assert round(se[0], 3) == 0.167 and round(se[1], 3) == 0.512 and round(between, 2) == 0.94


def test_the_launcher_runs_the_fifteen_parts_and_the_tests_first(h27):
    sh = (REPO_ROOT / "experiments" / "harness" / "run_h27.sh").read_text()
    assert "tests/test_h27_heart_forecast_size.py tests/test_h26_heart_scope_forms.py" in sh
    assert all(f'"{k}:' in sh for k in h27.TASK_KEYS) and "0-1 2-3 4-5 6-7 8-9" in sh
    for need in ("OMP_NUM_THREADS=8", "OPENBLAS_NUM_THREADS=8", "MKL_NUM_THREADS=8", "NUMPY_MADVISE_HUGEPAGE=0",
                 'TMPDIR="$S/tmp"', "S=~/v/scratch/h27", "status --porcelain", "python -I -B", "--commit"):
        assert need in sh, need
