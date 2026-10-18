"""Synthetic checks of size_procedure.py (no competition data). Run: uv run pytest experiments/coords -q

Growth rates in these tests are drawn at random: the procedure's correctness must not depend on any rate value.
"""

import hashlib
import io
import json
import re
import subprocess
import sys
import tokenize
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import size_procedure as sp  # noqa: E402

LAYOUT = (9.5, 10.5, 12.5)  # the 2026-10-20 heart layout: previous stage, reference stage, target


def cloud(sd, n=4000, seed=0, core_frac=0.5, step=15.0):
    """Gaussian cloud with per-axis SDs `sd`, z snapped to sections of the given step; first half is core."""
    rng = np.random.default_rng(seed)
    C = rng.normal(size=(n, 3)) * np.asarray(sd, float)
    C[:, 2] = np.round(C[:, 2] / step) * step
    core = np.zeros(n, bool)
    core[: int(n * core_frac)] = True
    return C, core


def geom(name, t, scale, sd=(100.0, 60.0, 30.0), seed=0):
    C, core = cloud(np.asarray(sd) * scale, seed=seed)
    return sp.geometry_from_arrays(name, t, C, core)


# ---------------------------------------------------------------------------------------------- rules and table
@pytest.mark.parametrize("seed", range(5))
def test_candidate_rules_are_multipliers_on_the_observed_growth(seed):
    rng = np.random.default_rng(seed)
    t1 = rng.uniform(8.0, 11.0)
    t2 = t1 + rng.uniform(0.2, 2.0)
    t3 = t2 + rng.uniform(0.2, 3.0)
    y1, g = rng.normal(5.0, 1.0), rng.normal(0.0, 1.0)  # any log size and any observed log growth, of either sign
    ts, ys = np.array([t1, t2]), np.array([y1, y1 + g])
    for r in sp.TABLE_RULES:
        m = sp.rule_multiplier(r, t1, t2, t3)
        assert sp.RULES[r](ts, ys, t3) - ys[-1] == pytest.approx(m * g, abs=1e-9)
    assert sp.rule_multiplier("latest_rate", t1, t2, t3) == pytest.approx((t3 - t2) / (t2 - t1))
    assert sp.rule_multiplier("half_rate", t1, t2, t3) == pytest.approx(0.5 * (t3 - t2) / (t2 - t1))
    assert sp.rule_multiplier("power_t", t1, t2, t3) == pytest.approx(np.log(t3 / t2) / np.log(t2 / t1))
    assert sp.rule_multiplier("copy_last", t1, t2, t3) == 0.0


def test_rule_table_closed_forms_for_the_release_layout():
    tab = {r["rule"]: r for r in sp.rule_table(*LAYOUT)}
    h = (LAYOUT[2] - LAYOUT[1]) / (LAYOUT[1] - LAYOUT[0])
    assert set(tab) == {"copy_last", "half_rate", "power_t", "latest_rate", "average"}
    m_avg = np.mean([sp.rule_multiplier(r, *LAYOUT) for r in sp.PREREGISTERED["heart_rules"]])
    assert tab["average"]["multiplier_m"] == pytest.approx(m_avg)
    for r in tab.values():
        m = r["multiplier_m"]
        assert r["horizon_ratio"] == pytest.approx(h)
        worst = max(abs(m - h * rho) for rho in np.linspace(*sp.RHO_RANGE, 1001))
        assert r["worst_case_bias_coef"] == pytest.approx(worst, abs=1e-12)
        assert r["noise_amplification"] == pytest.approx(np.sqrt((1 + m) ** 2 + m ** 2 + 1))
    # the pre-registered average has the smallest worst-case bias of the rules in the table
    assert tab["average"]["worst_case_bias_coef"] == min(r["worst_case_bias_coef"] for r in tab.values())


def test_docstring_table_matches_rule_table():
    rows = {r["rule"]: r for r in sp.rule_table(*LAYOUT)}
    found = 0
    for line in sp.__doc__.splitlines():
        parts = line.split()
        if len(parts) == 4 and parts[0] in rows:
            r = rows[parts[0]]
            got = [float(x) for x in parts[1:]]
            want = [r["multiplier_m"], r["worst_case_bias_coef"], r["noise_amplification"]]
            assert got == pytest.approx(want, abs=5e-4), parts[0]
            found += 1
    assert found == len(rows)


def test_decaying_rate_decelerates_and_power_law_is_between():
    rng = np.random.default_rng(7)
    ts = np.array([9.0, 9.5, 10.5])
    rb = rng.uniform(0.05, 1.0)
    ra = rb * rng.uniform(1.2, 3.0)  # decelerating, positive
    ys = np.array([0.0, ra * 0.5, ra * 0.5 + rb * 1.0]) + rng.normal(5.0, 1.0)
    d = sp.rule_decaying_rate(ts, ys, 12.5)
    assert ys[-1] < d < sp.rule_latest_rate(ts, ys, 12.5)
    p = sp.rule_power_t(ts, ys, 12.5)
    assert sp.rule_half_rate(ts, ys, 12.5) < p < sp.rule_latest_rate(ts, ys, 12.5)


# ---------------------------------------------------------------------------------------------- selection and validity
def test_applicability_and_validity():
    early = [geom("a", 8.25, 1.0), geom("b", 8.75, 1.0, seed=1), geom("c", 9.5, 1.8, seed=2)]
    res = sp.heart_size_from_geometry(early, 12.5, z_mode="raw")
    assert res["valid"] is False  # the only rate interval starts at 8.75 < 9.0
    late = [geom("b", 8.75, 1.0, seed=1), geom("c", 9.5, 1.8, seed=2), geom("d", 10.5, 2.4, seed=3)]
    res = sp.heart_size_from_geometry(late, 12.5, z_mode="raw")
    assert res["valid"] is True
    g = res["growth"]["core"]
    assert set(g["applicable"]) == {"latest_rate", "half_rate", "power_t"}
    assert g["common_pseudo_targets"] == []  # no pseudo-target whose last interval starts at or after E9.0
    assert set(g["chosen"]) == set(g["applicable"])  # so all three are averaged
    expected = np.mean([g["all_rule_predictions"][r] for r in g["chosen"]])
    assert np.log(res["pred_rms_core"]) == pytest.approx(expected)
    # in raw mode the predicted log growth is the average multiplier times the observed log growth
    obs_g = np.log(late[2].rms_core / late[1].rms_core)
    m_avg = [r for r in sp.rule_table(*LAYOUT) if r["rule"] == "average"][0]["multiplier_m"]
    assert res["pred_log_growth_core"] == pytest.approx(m_avg * obs_g)


def test_select_rules_compares_on_a_common_pseudo_target_set():
    def p(target, rule, err, eligible=True):
        return {"target": target, "kind": "next", "rule": rule, "error": err, "eligible": eligible}

    pt = [p("A", "r1", 0.10), p("A", "r2", 0.30), p("B", "r1", 0.50), p("B", "r2", 0.0, eligible=False),
          p("C", "r1", 0.20), p("C", "r2", 0.0)]
    chosen, mae, common = sp.select_rules(pt, ["r1", "r2"], tie_tol=0.01)
    assert common == ["A:next", "C:next"]  # B is left out: r2 is not eligible there
    assert mae == pytest.approx({"r1": 0.15, "r2": 0.15})  # over A and C only (r1's 0.5 at B does not count)
    assert set(chosen) == {"r1", "r2"}
    pt[1]["error"] = 0.5  # r2 now worse on the common set
    chosen, mae, _ = sp.select_rules(pt, ["r1", "r2"], tie_tol=0.05)
    assert chosen == ["r1"] and mae["r2"] == pytest.approx(0.25)
    chosen, _, common = sp.select_rules([p("A", "r1", 0.1), p("A", "r2", 0.2, eligible=False)], ["r1", "r2"], 0.0)
    assert common == [] and chosen == ["r1", "r2"]  # nothing to compare on: keep all


def test_split_mode_keeps_one_validity_per_series():
    late = [geom("b", 8.75, 1.0, seed=1), geom("c", 9.5, 1.8, seed=2), geom("d", 10.5, 2.4, seed=3)]
    res = sp.heart_size_from_geometry(late, 12.5, z_mode="split")
    assert set(res["growth"]) == {"core_xy", "core_z"}
    assert res["valid"] == all(g["valid"] for g in res["growth"].values())
    early = late[:2]
    res = sp.heart_size_from_geometry(early, 12.5, z_mode="split")
    assert res["valid"] is False and not any(g["valid"] for g in res["growth"].values())


@pytest.mark.parametrize("seed", range(3))
def test_scope_hypotheses(seed):
    rng = np.random.default_rng(seed)
    s_prev = rng.uniform(0.2, 1.0)
    s_ref = s_prev * rng.uniform(0.2, 0.9)
    ts = [9.5, 10.5]
    assert sp.scope_log_ratio([s_prev, s_ref], ts, 12.5, "H0") == pytest.approx(s_ref)
    assert sp.scope_log_ratio([s_prev, s_ref], ts, 12.5, "H2") == pytest.approx(0.5 * s_ref)
    assert sp.scope_log_ratio([s_prev, s_ref], ts, 12.5, "H1") == pytest.approx(s_ref * (s_ref / s_prev) ** 2)
    assert sp.scope_log_ratio([s_ref, s_prev], ts, 12.5, "H1") == pytest.approx(s_prev)  # increasing: held


def test_auto_z_mode_uses_in_plane_growth_when_the_thin_axis_is_z():
    flat = [geom("c", 9.5, 1.0, sd=(100, 80, 30), seed=2), geom("d", 10.5, 1.3, sd=(100, 80, 30), seed=3)]
    assert sp.resolve_heart_z_mode(flat, "auto") == "xy_grow_z_hold"
    tall = [geom("c", 9.5, 1.0, sd=(30, 80, 100), seed=2), geom("d", 10.5, 1.3, sd=(30, 80, 100), seed=3)]
    assert sp.resolve_heart_z_mode(tall, "auto") == "raw"
    res = sp.heart_size_from_geometry([geom("b", 9.0, 0.8, sd=(100, 80, 30), seed=1)] + flat, 12.5)
    assert res["z_mode"] == "xy_grow_z_hold"
    # z is held at the reference: predicted core RMS^2 = predicted xy^2 + reference z SD^2
    ref = flat[-1]
    assert res["pred_rms_core"] == pytest.approx(np.hypot(res["pred_rms_xy_core"], ref.sd_z_core))


def test_heart_report_logs_every_z_mode_and_scope_hypothesis():
    late = [geom("b", 8.75, 1.0, seed=1), geom("c", 9.5, 1.8, seed=2), geom("d", 10.5, 2.4, seed=3)]
    rep = sp.heart_report(late, 12.5)
    assert set(rep["z_modes"]) == set(sp.Z_MODES)
    assert [s["z_mode"] for s in rep["summary"]] == list(sp.Z_MODES)
    assert sum(s["preregistered"] for s in rep["summary"]) == 1
    for s in rep["summary"]:
        assert set(s["pred_rms_full"]) == {"H1", "H0", "H2"} and s["valid"] is True
    assert [s["z_mode"] for s in rep["summary"] if s["preregistered"]] == [rep["preregistered"]["z_mode"]]
    assert {r["rule"] for r in rep["rule_table"]} == {"copy_last", "half_rate", "power_t", "latest_rate", "average"}


def test_embryo_pchip_and_loo():
    gs = [geom("e1", 6.75, 1.0), geom("e2", 7.25, 1.2, seed=1), geom("e3", 7.5, 1.5, seed=2),
          geom("e4", 8.0, 2.4, seed=3)]
    res = sp.embryo_size_from_geometry(gs, 7.75)
    assert gs[2].rms_full < res["pred_rms_full"] < gs[3].rms_full  # monotone data: stays inside the bracket
    assert res["bracket"] == ["e3", "e4"]
    assert [r["left_out"] for r in res["loo_interior"]] == ["e2", "e3"]
    with pytest.raises(ValueError):
        sp.embryo_size_from_geometry(gs, 8.5)  # extrapolation is refused
    hold = sp.embryo_size_from_geometry(gs, 7.75, z_mode="xy_grow_z_hold")
    xy = sp.embryo_size_from_geometry(gs, 7.75, z_mode="split")
    assert hold["pred_rms_full"] > 0 and xy["pred_rms_full"] > 0
    rep = sp.embryo_report(gs, 7.75)
    assert set(rep["z_modes"]) == set(sp.Z_MODES) and rep["preregistered"]["z_mode"] == "raw"


def test_rescale_to_rms_is_a_proper_scaling():
    C, _ = cloud((50, 40, 30), seed=5)
    D = sp.rescale_to_rms(C, 123.0)
    assert sp.rms_radius(D) == pytest.approx(123.0)
    np.testing.assert_allclose(D.mean(0), C.mean(0), atol=1e-9)
    f = sp.rms_radius(D) / sp.rms_radius(C)
    np.testing.assert_allclose(D - D.mean(0), f * (C - C.mean(0)), atol=1e-8)


# ---------------------------------------------------------------------------------------------- no size numbers
# Every numeric literal allowed in size_procedure.py, by kind. Anything else (a radius, a rate, a log size) fails.
ALLOWED_LITERALS = {
    # E-days: stage times and the rate-interval cut-off
    6.75, 7.25, 7.5, 7.75, 8.0, 8.25, 8.5, 8.75, 9.0, 9.5, 10.5, 12.5,
    # angle (degrees) of the z-mode switch
    30.0,
    # tolerances: z-level snapping, rule tie tolerance (log units)
    1e-3, 0.15,
    # rule forms and structure: "half" (rate, scope, rho range end), indices, stage counts, minimum core cells
    0.0, 0.5, 1.0, 2.0, 3.0, 4.0,
}
_NUM = re.compile(r"(?<![\w.#])\d+(?:\.\d+)?(?:[eE][-+]?\d+)?")
_TEXT_TOKENS = {tokenize.STRING, tokenize.COMMENT} | ({tokenize.FSTRING_MIDDLE} if hasattr(tokenize, "FSTRING_MIDDLE")
                                                     else set())


def _tokens():
    return list(tokenize.generate_tokens(io.StringIO((HERE / "size_procedure.py").read_text()).readline))


def test_numeric_literals_in_code_are_on_the_allow_list():
    lits = {float(t.string) for t in _tokens() if t.type == tokenize.NUMBER}
    assert lits <= ALLOWED_LITERALS, sorted(lits - ALLOWED_LITERALS)


def test_numbers_in_docstrings_and_comments_are_on_the_allow_list():
    table = {round(v, 3) for r in sp.rule_table(*LAYOUT) for v in (r["multiplier_m"], r["worst_case_bias_coef"],
                                                                   r["noise_amplification"], r["horizon_ratio"])}
    allowed = {round(v, 3) for v in ALLOWED_LITERALS} | table
    bad = []
    for t in _tokens():
        if t.type in _TEXT_TOKENS:
            text = re.sub(r"\d{4}-\d{2}-\d{2}", " ", t.string)  # ISO dates
            bad += [(t.start[0], x) for x in _NUM.findall(text) if round(float(x), 3) not in allowed]
    assert not bad, bad


# ---------------------------------------------------------------------------------------------- end to end
def _write_stage(path: Path, sd, seed: int, n: int = 600, step: float = 15.0, duplicate: bool = False,
                 labels=("V-CM", "Endocardium", "Neural Tube", "SE")):
    """A tiny MERFISH-like h5ad: integer counts, log1p CP10k X, spatial_3D on sections, mapped labels."""
    import anndata as ad
    from scipy import sparse

    rng = np.random.default_rng(seed)
    counts = rng.poisson(3.0, size=(n, 40)).astype(np.float64)
    counts[:, 0] += 30  # every row has a library well above the duplicate check's minimum
    C = rng.normal(size=(n, 3)) * np.asarray(sd, float)
    C[:, 2] = np.round(C[:, 2] / step) * step
    if duplicate:  # one record copied into the next section up, same xy and counts
        counts[1] = counts[0]
        C[1] = C[0] + np.array([0.0, 0.0, step])
    X = np.log1p(counts / counts.sum(1, keepdims=True) * 1e4)
    a = ad.AnnData(X=sparse.csr_matrix(X.astype(np.float32)))
    a.layers["counts"] = sparse.csr_matrix(counts)
    a.obsm["spatial_3D"] = C
    a.obs["celltype"] = np.asarray(labels)[np.arange(n) % len(labels)]
    a.var_names = [f"g{i}" for i in range(counts.shape[1])]
    a.write_h5ad(path)
    return path


def test_command_line_writes_provenance_every_mode_and_duplicates(tmp_path):
    stages = [(_write_stage(tmp_path / "s_a.h5ad", (100, 60, 30), 1), 8.75),
              (_write_stage(tmp_path / "s_b.h5ad", (140, 90, 30), 2, duplicate=True), 9.5),
              (_write_stage(tmp_path / "s_c.h5ad", (180, 120, 30), 3), 10.5)]
    out = tmp_path / "heart.json"
    cmd = [sys.executable, str(HERE / "size_procedure.py"), "--board", "heart", "--target", "12.5",
           "--timestamp", "2026-10-20T00:00:00Z", "--out", str(out)]
    for p, t in stages:
        cmd += ["--stage", str(p), str(t)]
    subprocess.run(cmd, check=True, capture_output=True, text=True, cwd=sp.REPO)
    res = json.loads(out.read_text())
    assert res["timestamp"] == "2026-10-20T00:00:00Z" and res["valid"] is True
    assert set(res["z_modes"]) == set(sp.Z_MODES)
    assert all(set(s["pred_rms_full"]) == {"H1", "H0", "H2"} for s in res["summary"])
    prov = res["provenance"]
    assert prov["size_procedure_sha256"] == hashlib.sha256((HERE / "size_procedure.py").read_bytes()).hexdigest()
    for p, _ in stages:
        assert prov["inputs_sha256"][str(p)] == hashlib.sha256(p.read_bytes()).hexdigest()
    assert len(prov["git_commit"]) >= 7 and "configs/lineage_map.yaml" in prov["code_sha256"]
    dropped = {i["name"]: i["n_dropped_duplicates"] for i in res["inputs"]}
    assert dropped == {"s_a": 0, "s_b": 1, "s_c": 0}
    # redaction leaves no forecast number anywhere
    red = sp.redact_forecasts(res)
    flat = json.dumps(red)
    assert sp.REDACTED in flat and "all_rule_predictions\": {" not in flat
    assert all(v == sp.REDACTED for s in red["summary"] for k, v in s.items() if k.startswith("pred_"))


def test_unmapped_labels_are_refused(tmp_path):
    p = _write_stage(tmp_path / "s_x.h5ad", (100, 60, 30), 4, labels=("V-CM", "Not-A-Mapped-Label"))
    with pytest.raises(ValueError, match="labels not in configs/lineage_map.yaml"):
        sp.measure_stage(p, 9.5)
