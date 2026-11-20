"""Try T63 (experiments/harness/TRIES.md): experiments/t1_m1/double_step_read.py on synthetic stages, synthetic
directions and hand-made task tables only. No released file, statistics file, board file or score is read; the
committed relabelling record and the frozen settings file are read for their fields and hashes."""

from __future__ import annotations

import importlib.util
import itertools
import json
import sys
from pathlib import Path

import numpy as np
import pytest
from test_g_mom_verify import _small_inputs
from test_t1_scope import HELD
from test_t1_state_transport import _gmom_ctx

from vec import harness
from vec.core import counts_of
from vec.t1_scope import classify, load_spec, reweight

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "experiments" / "t1_m1" / "double_step_read.py"
DIRECTED = ("CM", "Endocardium", "Endothelium")
WITH_EPI = (*DIRECTED, "Epicardium")
EPI_SHARES = {"cm": 0.5, "endocardium": 0.2, "epicardium": 0.15, "blood": 0.05, "surrounding": 0.1}


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture
def ds(monkeypatch, tmp_path):
    mod = _load("double_step_read_module", SCRIPT)
    monkeypatch.setattr(mod, "RESULTS", tmp_path / "results")
    return mod


def _dense(x):
    return np.asarray(x.todense()) if hasattr(x, "todense") else np.asarray(x)


def _directions(n_genes, groups=WITH_EPI, size=0.4, seed=11, nan_every=17):
    rng = np.random.default_rng(seed)
    out = {}
    for g in groups:
        d = rng.normal(0.0, size, n_genes)
        d[::nan_every] = np.nan
        out[g] = d
    return out


def _halves(directions, seed=3):
    rng = np.random.default_rng(seed)
    out = {}
    for g, d in directions.items():
        n = len(d)
        out[g] = {"mean_cp": rng.gamma(1.0, 5.0, n), "d_a": d + rng.normal(0, 0.3, n),
                  "d_b": d + rng.normal(0, 0.3, n)}
    return out


HELD_M1 = sorted({*HELD})


def _held_m1(genes):
    return sorted({*HELD, genes[3], genes[7]} & set(genes))


# ---------------------------------------------------------------------------------------------------------------------
# Constants, the source guard and the refusals


def test_the_fractions_the_arms_and_the_count(ds):
    assert ds.feng_dirs.bracket_fraction(9.5, 10.5) == pytest.approx(ds.PRACTICE_FRACTION, abs=1e-12)
    assert abs(ds.feng_dirs.bracket_fraction(10.5, 12.5) - ds.TEST_FRACTION) <= 1e-12
    assert ds.TEST_FRACTION == pytest.approx(2 * ds.PRACTICE_FRACTION, abs=1e-12)
    assert ds.ARMS[:2] == ("none_gmom", "m1rc_gmom") and ds.GATED == ("none_gmom", "m1rc_gmom")
    assert ds.COUNTED == ("m1rc_gmom_x2", "m1rc_gmom_x2_epi_held", "m1rc_gmom_x2_epi_eb", "placebo_x1", "placebo_x2",
                          "placebo_x2_q1", "placebo_x2_q2", "placebo_x2_q3", "placebo_x2_q4")
    assert ds.TRIES_COUNT == 45 and ds.TRY_RANGE == (368, 413)
    assert ds.LIB_ARMS[0] == "copy_last" and ds.PLACEBO_OFFSET["placebo_x1"] == ds.PLACEBO_OFFSET["placebo_x2"] == 0
    assert sorted(ds.PLACEBO_OFFSET[q] for q in ds.X2_PLACEBOS) == [0, 1, 2, 3, 4]
    assert ds.arm_fractions(0.2) == pytest.approx(ds.FRACTIONS)
    ds.check_fractions(ds.feng_dirs.bracket_fraction(9.5, 10.5), "x")
    with pytest.raises(SystemExit):
        ds.check_fractions(0.25, "x")


def test_the_source_guard(ds, tmp_path):
    assert ds.source_guard() == [] and ds.rc.source_guard() == []
    for text in ("import perm_read as pr\n", "import eb_read\n", "from step_read import step_arms\n",
                 "import deep_read as dr, json\n", "x = __imp" + "ort__('perm_read')\n",
                 'x = fd.feng_directions(list(g))\n', 'p = fd.clean_' + 'stats_path()\n',
                 'open("n54' + '_relabel_listing")\n', 'g = "auth' + 'ors"\n'):
        bad = tmp_path / "bad.py"
        bad.write_text(text)
        assert ds.source_guard(bad), text
    ok = tmp_path / "ok.py"
    ok.write_text("import json\nfrom vec import harness\n# perm_read is named in a comment only\n")
    assert ds.source_guard(ok) == []


def test_the_runs_and_the_read_are_refused_without_the_row(ds, tmp_path, monkeypatch):
    tries = tmp_path / "TRIES.md"
    tries.write_text("| T62 (one row) | T1 | x |\n")
    monkeypatch.setattr(ds, "TRIES", tries)
    monkeypatch.setattr(ds.rc, "relabel_checks", lambda *a, **k: pytest.fail("nothing before the row's refusal"))
    monkeypatch.setattr(ds.gv, "tree_state", lambda: pytest.fail("nothing before the row's refusal"))
    for argv in (["run-cut", "--cut", "chamber", "--out", "x"], ["run-lib", "--library", "3", "--out", "x"],
                 ["read", "--out", "x"]):
        with pytest.raises(SystemExit, match="no row"):
            ds.main(argv)


def test_a_run_is_refused_on_an_existing_file_an_unclean_tree_and_a_failed_check(ds, tmp_path, monkeypatch):
    tries = tmp_path / "TRIES.md"
    tries.write_text("| T63 (one row) | T1 | x |\n")
    monkeypatch.setattr(ds, "TRIES", tries)
    out = tmp_path / "out"
    out.mkdir()
    (out / "lib3.json").write_text("{}")
    with pytest.raises(SystemExit, match="already holds"):
        ds.refuse(out, "lib3.json")
    monkeypatch.setattr(ds.gv, "tree_state", lambda: {"commit": "c", "clean": False})
    with pytest.raises(SystemExit, match="not clean"):
        ds.refuse(out, "lib1.json")
    monkeypatch.setattr(ds.gv, "tree_state", lambda: {"commit": "c", "clean": True})
    monkeypatch.setattr(ds.rc, "relabel_checks", lambda *a, **k: {"passes": False, "reasons": ["the yaml"]})
    with pytest.raises(SystemExit, match="the yaml"):
        ds.refuse(out, "lib1.json")
    monkeypatch.setattr(ds, "source_guard", lambda: ["the source imports ['perm_read']"])
    with pytest.raises(SystemExit, match="perm_read"):
        ds.refuse(out, "lib1.json")


# ---------------------------------------------------------------------------------------------------------------------
# The placebo and the shrunk direction (target-free)


def test_the_placebo_is_row_t45s_byte_for_byte(ds):
    pr = _load("perm_read_for_t63_test", ROOT / "experiments" / "t1_m1" / "perm_read.py")
    dirs = _directions(400)
    halves = _halves(dirs)
    for offset in range(5):
        ours, theirs = ds.placebo_directions(dirs, halves, offset), pr.placebo_directions(dirs, halves, offset)
        assert set(ours) == set(theirs) == set(WITH_EPI)
        for g in ours:
            assert ours[g].tobytes() == theirs[g].tobytes()
    free = ds.free_genes(HELD_M1, [f"G{i}" for i in range(400)])
    assert ds.between_bin_share(dirs, halves) == pytest.approx(pr.between_bin_share(dirs, halves), nan_ok=True)
    assert ds.between_bin_share(dirs, halves, free) == pytest.approx(pr.between_bin_share(dirs, halves, free))
    with pytest.raises(ValueError, match="no halves"):
        ds.placebo_directions(dirs, {g: h for g, h in halves.items() if g != "Epicardium"})
    with pytest.raises(ValueError, match="permutation seed"):
        ds.placebo_directions({"Blood": dirs["CM"]}, {"Blood": halves["CM"]})


def test_the_placebo_keeps_the_size_profile_and_none_of_the_gene_order(ds):
    genes = [f"G{i}" for i in range(500)]
    dirs = _directions(500)
    halves = _halves(dirs)
    placebos = ds.placebo_set(dirs, halves)
    assert placebos["placebo_x1"] is placebos["placebo_x2"]  # one permutation at both sizes
    for arm in ds.X2_PLACEBOS:
        for g, d in dirs.items():
            p = placebos[arm][g]
            ok = np.isfinite(d)
            assert np.array_equal(np.isfinite(p), ok)  # NaN stays NaN
            basis = halves[g]["mean_cp"][ok]
            for b in np.array_split(np.argsort(basis, kind="stable"), ds.BINS):  # each bin keeps its values
                assert np.array_equal(np.sort(p[ok][b]), np.sort(d[ok][b]))
            assert np.mean(p[ok] == d[ok]) < 0.2  # not the direction's gene order
            assert abs(np.corrcoef(p[ok], d[ok])[0, 1]) < 0.2
    shas = {a: ds.direction_sha(placebos[a]) for a in ds.X2_PLACEBOS}
    for a, b in itertools.combinations(ds.X2_PLACEBOS, 2):
        assert all(shas[a][g] != shas[b][g] for g in dirs)
    ratios = ds.size_ratios(placebos["placebo_x2"], dirs, [], genes)
    assert all(abs(r["mean_abs_logf"] - 1) < 1e-12 and abs(r["rms_logf"] - 1) < 1e-12 for r in ratios.values())
    held = genes[:250]
    ratios = ds.size_ratios(placebos["placebo_x2"], dirs, held, genes)
    assert all(abs(r["mean_abs_logf"] - 1) < 0.3 for r in ratios.values())
    big = {g: 2.0 * v for g, v in dirs.items()}
    assert all(r["flagged"] for r in ds.size_ratios(big, dirs, [], genes).values())


def test_the_shrunk_epicardium_direction(ds):
    dirs = _directions(300)
    halves = _halves(dirs)
    shrunk, factors = ds.epi_shrunk(dirs, halves)
    want, f = ds.feng_dirs.eb_shrink(dirs["Epicardium"], halves["Epicardium"]["d_a"], halves["Epicardium"]["d_b"],
                                     halves["Epicardium"]["mean_cp"], 10)
    assert factors == pytest.approx(f) and shrunk["Epicardium"].tobytes() == want.tobytes()
    assert all(shrunk[g] is dirs[g] for g in DIRECTED)
    same, none = ds.epi_shrunk({g: dirs[g] for g in DIRECTED}, halves)
    assert none is None and set(same) == set(DIRECTED)
    rec = ds.direction_record(dirs, ds.placebo_set(dirs, halves), shrunk, factors, halves, [],
                              [f"G{i}" for i in range(300)])
    assert rec["epi_eb_factor_mean"] == pytest.approx(np.mean(factors)) and ds.placebo_reasons(rec) == []


def test_recomputed_halves_refuses_a_group_without_its_own(ds, monkeypatch):
    dirs = _directions(10, DIRECTED)
    monkeypatch.setattr(ds.feng_dirs, "recomputed_stats_path", lambda *a, **k: Path("x"))
    monkeypatch.setattr(ds.feng_dirs, "check_recomputed", lambda *a, **k: ({}, {}))
    monkeypatch.setattr(ds.feng_dirs.dl, "load_stats", lambda p: {})
    monkeypatch.setattr(ds.feng_dirs, "check_recomputed_keys", lambda st: None)
    monkeypatch.setattr(ds.feng_dirs, "load_config", lambda *a: {})
    monkeypatch.setattr(ds.feng_dirs, "load_relabel_config", lambda *a: {})
    monkeypatch.setattr(ds.feng_dirs, "groups_config", lambda c, r: {})
    seen = {}

    def halves(st, genes, cfg, bracket, groups=None):
        seen["bracket"], seen["groups"] = bracket, set(groups)
        return {g: {} for g in groups if g != "Endothelium"}

    monkeypatch.setattr(ds.feng_dirs, "direction_halves", halves)
    with pytest.raises(ValueError, match="Endothelium"):
        ds.recomputed_halves(["a"], dirs)
    assert seen == {"bracket": ("E9.5", "E14.5"), "groups": set(DIRECTED)}


# ---------------------------------------------------------------------------------------------------------------------
# The arms on synthetic stages


def _epi_world():
    from test_t1_state_transport import classified_pair

    start, ref = classified_pair(seed=5, n_scale=0.5)
    return ref, start


def _lib_setup(ds, ref, start):
    spec = load_spec()
    fb = {"E8.5": spec.fallbacks("E8.5_RNA"), "E9.5": spec.fallbacks("E9.5_RNA")}
    within = {"surrounding": ["Foregut"]}
    factors = ds.br.gmom_fits(ref, "E9.5_RNA", start, "E8.5_RNA")["file"]
    return fb, within, factors


def _kit(ds, n_vars, groups=WITH_EPI):
    dirs = _directions(n_vars, groups)
    halves = _halves(dirs)
    shrunk, _ = ds.epi_shrunk(dirs, halves)
    return dirs, ds.placebo_set(dirs, halves), shrunk


def test_the_library_arms(ds):
    e95, e85 = _epi_world()
    fb, within, factors = _lib_setup(ds, e95, e85)
    genes = list(e95.var_names)
    held_m1 = _held_m1(genes)
    dirs, placebos, shrunk = _kit(ds, e95.n_vars)
    arms = ds.library_arms(e95, EPI_SHARES, factors, HELD, held_m1, 6, within, fb, dirs, placebos, shrunk, 0.2,
                           n_cells=300)
    assert tuple(arms) == ds.ARMS
    c = {a: _dense(counts_of(v)) for a, v in arms.items()}
    # the gates: T62's none_gmom and m1rc_gmom, byte for byte
    t62 = ds.rc.library_pair(e95, EPI_SHARES, factors, HELD, held_m1, 6, within, fb, dirs, 0.2, n_cells=300)
    for a in ("none_gmom", "m1rc_gmom"):
        assert np.array_equal(c[a], _dense(counts_of(t62[a]))) and (_dense(arms[a].X) == _dense(t62[a].X)).all()
    # the doubled arm is T62's candidate at 2/5, and equals the single arm when the fraction is set back
    t62x2 = ds.rc.library_pair(e95, EPI_SHARES, factors, HELD, held_m1, 6, within, fb, dirs, 0.4, n_cells=300)
    assert np.array_equal(c["m1rc_gmom_x2"], _dense(counts_of(t62x2["m1rc_gmom"])))
    assert np.array_equal(c["m1rc_gmom_x2_epi_held"], _dense(counts_of(t62x2["m1rc_gmom_epi_held"])))
    back = ds.library_arms(e95, EPI_SHARES, factors, HELD, held_m1, 6, within, fb, dirs, placebos, shrunk, 0.1,
                           n_cells=300)
    assert np.array_equal(_dense(counts_of(back["m1rc_gmom_x2"])), c["m1rc_gmom"])
    assert not np.array_equal(c["m1rc_gmom_x2"], c["m1rc_gmom"])
    base, _ = reweight(e95, EPI_SHARES, 300, 6, within=within, fallback=fb["E9.5"])
    labels = base.obs["celltype"].astype(str).to_numpy()
    groups = ds.feng_dirs.label_groups(sorted(set(labels)))
    epi = np.array([groups.get(lab) == "Epicardium" for lab in labels])
    covered = np.array([groups.get(lab) in dirs for lab in labels])
    drawn = _dense(counts_of(base))
    assert epi.any() and (~covered).any()
    # the epi-held arm differs from the doubled arm only in Epicardium cells, which keep their drawn counts
    differs = (c["m1rc_gmom_x2"] != c["m1rc_gmom_x2_epi_held"]).any(axis=1)
    assert differs.any() and not differs[~epi].any() and np.array_equal(c["m1rc_gmom_x2_epi_held"][epi], drawn[epi])
    # the shrunk arm differs from the doubled arm only in Epicardium cells, and is the doubled arm with no shrinkage
    differs = (c["m1rc_gmom_x2"] != c["m1rc_gmom_x2_epi_eb"]).any(axis=1)
    assert differs.any() and not differs[~epi].any()
    same = ds.library_arms(e95, EPI_SHARES, factors, HELD, held_m1, 6, within, fb, dirs, placebos, dict(dirs), 0.2,
                           n_cells=300)
    assert np.array_equal(_dense(counts_of(same["m1rc_gmom_x2_epi_eb"])), c["m1rc_gmom_x2"])
    # the placebo moves the same cells: every arm leaves the cells of undirected groups as drawn
    for a in ds.ARMS:
        assert np.array_equal(c[a][~covered], drawn[~covered]), a
    for a in ("placebo_x1", "placebo_x2", "placebo_x2_q3"):
        moved = (c[a] != drawn).any(axis=1)
        assert moved[covered].mean() > 0.5 and not np.array_equal(c[a], c["m1rc_gmom_x2"])
    # the placebo arm is the candidate's build on the placebo directions
    pl = ds.library_arms(e95, EPI_SHARES, factors, HELD, held_m1, 6, within, fb, placebos["placebo_x2"], placebos,
                         shrunk, 0.2, n_cells=300)
    assert np.array_equal(_dense(counts_of(pl["m1rc_gmom_x2"])), c["placebo_x2"])
    assert np.array_equal(_dense(counts_of(pl["m1rc_gmom"])), c["placebo_x1"])


def test_the_cut_arms(ds):
    ref, start = _epi_world()
    ctx = _gmom_ctx(ref, start, seed=1000003)
    fits = ds.br.gmom_fits(ref, "E9.5_RNA", start, "E8.5_RNA")
    held_m1 = _held_m1(list(ref.var_names))
    dirs, placebos, shrunk = _kit(ds, ref.n_vars)
    arms = ds.cut_arms(ctx, fits["file"], dirs, placebos, shrunk, 0.2, held_m1)
    assert tuple(arms) == ds.ARMS
    X = {a: np.asarray(v.X) for a, v in arms.items()}
    t62 = ds.rc.cut_pair(ctx, fits["file"], dirs, 0.2, held_m1)
    assert X["none_gmom"].tobytes() == t62["none_gmom"].X.tobytes()
    assert X["m1rc_gmom"].tobytes() == t62["m1rc_gmom"].X.tobytes()
    t62x2 = ds.rc.cut_pair(ctx, fits["file"], dirs, 0.4, held_m1)
    assert X["m1rc_gmom_x2"].tobytes() == t62x2["m1rc_gmom"].X.tobytes()
    assert X["m1rc_gmom_x2_epi_held"].tobytes() == t62x2["m1rc_gmom_epi_held"].X.tobytes()
    back = ds.cut_arms(ctx, fits["file"], dirs, placebos, shrunk, 0.1, held_m1)
    assert np.asarray(back["m1rc_gmom_x2"].X).tobytes() == X["m1rc_gmom"].tobytes()
    assert X["m1rc_gmom_x2"].tobytes() != X["m1rc_gmom"].tobytes()
    seen = []
    real = ds.rc.hold_group

    def spy(moved, other, groups, group="Epicardium"):
        seen.append(np.array([groups.get(x) == group for x in other.obs["celltype"].astype(str)]))
        return real(moved, other, groups, group)

    ds.rc.hold_group = spy
    try:
        ds.cut_arms(ctx, fits["file"], dirs, placebos, shrunk, 0.2, held_m1)
        same = ds.cut_arms(ctx, fits["file"], dirs, placebos, dict(dirs), 0.2, held_m1)
    finally:
        ds.rc.hold_group = real
    epi = seen[0]
    for a in ("m1rc_gmom_x2_epi_held", "m1rc_gmom_x2_epi_eb"):
        differs = (X["m1rc_gmom_x2"] != X[a]).any(axis=1)
        assert epi.any() and differs.any() and not differs[~epi].any(), a
    # with no shrinkage the spliced cells are the build's own direction step: the arm is the doubled arm
    assert np.asarray(same["m1rc_gmom_x2_epi_eb"].X).tobytes() == X["m1rc_gmom_x2"].tobytes()
    pl = ds.cut_arms(ctx, fits["file"], placebos["placebo_x2_q2"], placebos, shrunk, 0.2, held_m1)
    assert np.asarray(pl["m1rc_gmom_x2"].X).tobytes() == X["placebo_x2_q2"].tobytes()
    assert X["placebo_x2"].tobytes() != X["m1rc_gmom_x2"].tobytes()


def test_without_an_epicardium_direction_the_epicardium_arms_are_the_doubled_arm(ds):
    ref, start, _ = _small_inputs()
    classify(ref, source="E9.5_RNA")
    ctx = _gmom_ctx(ref, start, seed=1000002)
    fits = ds.br.gmom_fits(ref, "E9.5_RNA", start, "E8.5_RNA")
    dirs, placebos, shrunk = _kit(ds, ref.n_vars, DIRECTED)
    arms = ds.cut_arms(ctx, fits["file"], dirs, placebos, shrunk, 0.2, _held_m1(list(ref.var_names)))
    for a in ("m1rc_gmom_x2_epi_held", "m1rc_gmom_x2_epi_eb"):
        assert arms[a].X.tobytes() == arms["m1rc_gmom_x2"].X.tobytes()
    assert harness.finish  # the harness's finishing step is the one the arms go through


# ---------------------------------------------------------------------------------------------------------------------
# The reading on hand-made tables


def _arm(score):
    return {"score": float(score), "groups": {g: score / 100 for g in ("de_recovery", "direction", "distribution",
                                                                      "covariance")}}


GAINS = {"m1rc_gmom": 1.0, "m1rc_gmom_x2": 1.5, "m1rc_gmom_x2_epi_held": 1.0, "m1rc_gmom_x2_epi_eb": 1.3,
         "placebo_x1": 0.1, "placebo_x2": 0.1, "placebo_x2_q1": 0.1, "placebo_x2_q2": 0.1, "placebo_x2_q3": 0.1,
         "placebo_x2_q4": 0.1}


def _seq(v, s):
    return v[s] if isinstance(v, (list, tuple)) else v


def _directions_record(ds):
    cand = {g: g * 4 for g in WITH_EPI}
    pl = {a: {g: f"{g}{o}" for g in WITH_EPI} for a, o in ds.PLACEBO_OFFSET.items()}
    return {"candidate_sha256": cand, "placebo_sha256": pl,
            "size_ratios": {a: {g: {"mean_abs_logf": 1.0, "rms_logf": 1.0, "flagged": False} for g in WITH_EPI}
                            for a in ds.PLACEBO_OFFSET}}


def _relabel(ds, **change):
    return {"sha256": ds.rc.STATS_SHA256, "stats_sha256_hashed": ds.rc.STATS_SHA256, "yaml_sha256": ds.rc.YAML_SHA256,
            "route": "recomputed", "record_route": "recomputed", "authors_labels_read": False,
            "epicardium_directed": True, **change}


def _task(ds, lib=False, seeds=range(10), commit="c" * 40, clean=True, inverted=False, relabel=None, control=50.0,
          **gains):
    g = {**GAINS, **gains}
    rows = []
    for s in seeds:
        arms = {"none_gmom": _arm(control + s / 100)}
        arms.update({a: _arm(control + s / 100 + _seq(v, s)) for a, v in g.items()})
        sha = {a: f"{a}-{s}" for a in arms}
        row = {"seed": s, "arms": arms, "sha256": sha}
        if lib:
            arms["copy_last"] = _arm(40.0)
            sha["copy_last"] = f"copy-{s}"
            row["inverted_anchors"] = ["x"] if inverted and s == 3 else []
        else:
            row["floor"] = 50.0
        rows.append(row)
    return {"rows": rows, "stamp": {"commit": commit, "tree_clean": clean}, "fraction": 0.2,
            "fractions": ds.arm_fractions(0.2), "relabel": relabel if relabel is not None else _relabel(ds),
            "directions": _directions_record(ds)}


def _committed(control=50.0, gain=1.0):
    return {k: {"none_gmom": {s: control + s / 100 for s in range(10)},
                "m1rc_gmom": {s: control + s / 100 + gain for s in range(10)}}
            for k in ("chamber", "whole", "1", "2", "3")}


def _loaded(ds, chamber=None, whole=None, libs=None, **tasks):
    out = {"chamber": _task(ds, **(chamber or {})), "whole": _task(ds, **(whole or {}))}
    for lib in ("1", "2", "3"):
        out[lib] = _task(ds, lib=True, **((libs or {}).get(lib) or {}))
    out.update(tasks)
    return out


OK = {"reasons": []}


def _read(ds, loaded, committed=None, source=()):
    return ds.t63_reading(loaded, _committed() if committed is None else committed, OK, source)


def test_the_default_table_reads_build_as_ruled(ds):
    r = _read(ds, _loaded(ds))
    assert r["outcome"] == "read" and r["S1"]["holds"] and r["S1"]["words"] == ds.S1_WORDS[True]
    assert r["S2"] == "GAINS" and r["S3"] == "PAYS" and r["S4_x2"]["words"] == ds.NEEDS
    assert r["S4_x1_reported"]["words"] == ds.NEEDS
    assert r["follows"]["key"] == ds.BUILD and r["j6_range"] == pytest.approx([0.0, 1.5])
    lines = ds.outcome_lines(r)
    assert lines[0].startswith('S1: "the doubled step is AHEAD of the no-outside-data analogue" (the test\'s absolute '
                               "step on one-day truths at the released depth ratio; not a check of the two-day horizon")
    assert "never evidence for any fraction" in lines[1] and lines[-1] == r["follows"]["line"]
    inc = r["reported"]["increment_vs_placebo_increment"]["chamber"]
    assert inc["mean"] == pytest.approx(0.5)  # (1.5 - 1.0) - (0.1 - 0.1)


@pytest.mark.parametrize("chamber,whole,holds", [
    ({"m1rc_gmom_x2": 0.3}, {}, False),  # the bar is strict: above +0.3
    ({"m1rc_gmom_x2": 0.31}, {}, True),
    ({"m1rc_gmom_x2": [1.0] * 8 + [-0.1] * 2}, {}, True),  # 8 of 10 above 0
    ({"m1rc_gmom_x2": [1.0] * 7 + [-0.1] * 3}, {}, False),  # 7 of 10
    ({}, {"m1rc_gmom_x2": -0.3}, True),  # the whole cut at its floor
    ({}, {"m1rc_gmom_x2": -0.31}, False),
])
def test_s1_on_the_cuts(ds, chamber, whole, holds):
    r = _read(ds, _loaded(ds, chamber=chamber, whole=whole))
    assert r["S1"]["holds"] is holds and r["S1"]["cuts_ahead"] is holds
    assert r["S1"]["words"] == ds.S1_WORDS[holds]


def test_s1_needs_the_library_check(ds):
    r = _read(ds, _loaded(ds, libs={"2": {"m1rc_gmom_x2": -0.31}}))
    assert not r["S1"]["holds"] and r["S1"]["cuts_ahead"] and r["S1"]["libraries_below_floor"] == ["2"]
    r = _read(ds, _loaded(ds, libs={"2": {"m1rc_gmom_x2": -0.3}}))
    assert r["S1"]["holds"]
    r = _read(ds, _loaded(ds, libs={"1": {"inverted": True}, "2": {"inverted": True}}))
    assert r["valid_libraries"] == ["3"] and not r["S1"]["holds"]
    r = _read(ds, _loaded(ds, libs={"1": {"inverted": True}, "2": {"m1rc_gmom_x2": -5.0, "inverted": True}}))
    assert not r["S1"]["holds"]
    r = _read(ds, _loaded(ds, libs={"2": {"m1rc_gmom_x2": -5.0, "inverted": True}}))
    assert r["S1"]["holds"]  # an invalid library's mean does not count


@pytest.mark.parametrize("x2,word", [(1.31, "GAINS"), (1.3, "level"), (0.71, "level"), (0.69, "LOSES"), (0.0, "LOSES"),
                                     ([2.0] * 8 + [0.9] * 2, "GAINS"), ([2.0] * 7 + [0.9] * 3, "level"),
                                     ([0.0] * 8 + [1.1] * 2, "LOSES"), ([0.0] * 7 + [1.1] * 3, "level")])
def test_s2(ds, x2, word):
    r = _read(ds, _loaded(ds, chamber={"m1rc_gmom_x2": x2}))
    assert r["S2"] == word


@pytest.mark.parametrize("held,word", [(1.19, "PAYS"), (1.2, "within noise"), (1.8, "within noise"), (1.81, "COSTS"),
                                       ([1.0] * 8 + [1.6] * 2, "PAYS"), ([1.0] * 7 + [1.6] * 3, "within noise"),
                                       ([2.0] * 8 + [1.4] * 2, "COSTS"), ([2.0] * 7 + [1.4] * 3, "within noise")])
def test_s3_on_library_3(ds, held, word):
    r = _read(ds, _loaded(ds, libs={"3": {"m1rc_gmom_x2_epi_held": held}}))
    assert r["S3"] == word and set(r["S3_tasks"]) == {"chamber", "whole", "1", "2", "3"}


def test_s3_reads_library_3_alone_and_is_unread_when_it_cannot(ds):
    r = _read(ds, _loaded(ds, chamber={"m1rc_gmom_x2_epi_held": 5.0}, libs={"1": {"m1rc_gmom_x2_epi_held": 5.0}}))
    assert r["S3"] == "PAYS" and r["S3_tasks"]["chamber"]["mean"] == pytest.approx(-3.5)
    r = _read(ds, _loaded(ds, libs={"3": {"m1rc_gmom_x2_epi_held": 5.0, "inverted": True}}))
    assert r["S3"] == "unread" and not r["S3_task_valid"]
    r = _read(ds, _loaded(ds, libs={"3": {"m1rc_gmom_x2_epi_held": 5.0,
                                          "relabel": _relabel(ds, epicardium_directed=False)}}))
    assert r["S3"] == "unread" and r["follows"]["key"] == ds.BUILD  # unread counts as not COSTS


def _s4(ds, cut_gains):
    return _read(ds, _loaded(ds, chamber=cut_gains, whole=cut_gains))["S4_x2"]


def test_s4_categories(ds):
    assert _s4(ds, {})["words"] == ds.NEEDS
    # a placebo passes against none_gmom (any one permutation) while the direction beats all five: direction more
    assert _s4(ds, {"placebo_x2_q3": 0.5})["words"] == ds.MORE
    # the direction fails to beat one permutation and that placebo passes: a placebo gains as well
    assert _s4(ds, {"placebo_x2_q1": 1.4})["words"] == ds.PLACEBO_GAINS
    # it fails to beat one permutation and no placebo passes: not separated
    s = _s4(ds, {"m1rc_gmom_x2": 0.35, "placebo_x2_q4": 0.29})
    assert s["words"] == ds.NOT_SEPARATED and s["worst"] == "placebo_x2_q4"
    # the 8-of-10 edges, on the direction against one permutation and on a placebo against none_gmom
    assert _s4(ds, {"placebo_x2_q2": [1.6] * 3 + [0.0] * 7})["words"] == ds.NOT_SEPARATED  # beats on 7 of 10
    assert _s4(ds, {"placebo_x2_q2": [1.6] * 2 + [0.0] * 8})["words"] == ds.NEEDS
    assert _s4(ds, {"placebo_x2": [0.5] * 8 + [-0.1] * 2})["words"] == ds.MORE
    assert _s4(ds, {"placebo_x2": [0.6] * 7 + [-0.1] * 3})["words"] == ds.NEEDS
    # the whole cut's floor applies to both comparisons
    r = _read(ds, _loaded(ds, whole={"placebo_x2_q1": 1.81}))
    assert r["S4_x2"]["words"] == ds.NOT_SEPARATED  # behind by 0.31 on the whole cut: not beaten; it does not pass
    r = _read(ds, _loaded(ds, whole={"placebo_x2_q1": 1.8}))
    assert r["S4_x2"]["words"] == ds.NEEDS


def test_s4_at_one_fifth_is_reported_and_decides_nothing(ds):
    r = _read(ds, _loaded(ds, chamber={"placebo_x1": 0.95}, whole={"placebo_x1": 0.95}))
    assert r["S4_x1_reported"]["words"] == ds.PLACEBO_GAINS and r["follows"]["key"] == ds.BUILD
    assert r["S4_x1_reported"]["permutations"] == ["placebo_x1"]


def test_every_combination_maps_to_exactly_one_line(ds):
    seen = set()
    for s1, s3, s4 in itertools.product((True, False), ds.S3_WORDS, ds.S4_WORDS):
        f = ds.follows(s1, s3, s4)
        assert f["key"] in ds.FOLLOW_KEYS and f["line"].startswith("What follows:")
        assert f["files"][0] == ds.RULED
        build = s1 and s3 != "COSTS" and s4 in (ds.NEEDS, ds.MORE)
        assert (f["key"] == ds.BUILD) is build
        assert (ds.THIRD in f["files"]) is (s3 == "COSTS")
        assert (ds.NONE_FILE in f["files"]) is ((not s1) or s4 in (ds.PLACEBO_GAINS, ds.NOT_SEPARATED))
        if not build:
            assert "his ruling (the recomputed file) stands unless he changes it" in f["line"].lower()
            assert "nothing is built until he answers" in f["line"]
        seen.add(f["key"])
    v = ds.follows(None, None, None, void=True)
    assert v["key"] == ds.VOID and "ruled file is not built until his word" in v["line"]
    seen.add(v["key"])
    assert seen == set(ds.FOLLOW_KEYS)
    with pytest.raises(ValueError):
        ds.follows(True, "maybe", ds.NEEDS)


def test_the_outcomes_drive_what_follows(ds):
    r = _read(ds, _loaded(ds, libs={"3": {"m1rc_gmom_x2_epi_held": 2.0}}))
    assert r["S3"] == "COSTS" and r["follows"]["key"] == "ask+third"
    r = _read(ds, _loaded(ds, chamber={"m1rc_gmom_x2": 0.2}))
    assert r["follows"]["key"] == "ask+none"
    r = _read(ds, _loaded(ds, chamber={"placebo_x2_q1": 1.4, "m1rc_gmom_x2": 1.5},
                          whole={"placebo_x2_q1": 1.4}, libs={"3": {"m1rc_gmom_x2_epi_held": 2.0}}))
    assert r["S4_x2"]["words"] == ds.PLACEBO_GAINS and r["follows"]["key"] == "ask+third+none"


# ---------------------------------------------------------------------------------------------------------------------
# The gates and the void conditions


@pytest.mark.parametrize("case", ["control_off", "candidate_off", "missing_seed", "missing_arm", "missing_sha",
                                  "missing_task", "two_commits", "unclean", "stats", "yaml", "route", "authors",
                                  "checks", "source", "fraction", "x2_fraction", "placebo_is_candidate",
                                  "placebo_pair", "permutations_shared", "placebo_missing", "x2_is_x1",
                                  "placebo_is_x2"])
def test_every_guard_voids(ds, case):
    loaded, committed, checks, source = _loaded(ds), _committed(), OK, []
    if case == "control_off":
        committed["2"]["none_gmom"][4] += 2e-9
    elif case == "candidate_off":
        committed["chamber"]["m1rc_gmom"][9] -= 2e-9
    elif case == "missing_seed":
        loaded["whole"] = _task(ds, seeds=range(9))
    elif case == "missing_arm":
        del loaded["3"]["rows"][0]["arms"]["copy_last"]
    elif case == "missing_sha":
        del loaded["chamber"]["rows"][2]["sha256"]["placebo_x2_q4"]
    elif case == "missing_task":
        loaded["1"] = None
    elif case == "two_commits":
        loaded["1"] = _task(ds, lib=True, commit="d" * 40)
    elif case == "unclean":
        loaded["chamber"] = _task(ds, clean=False)
    elif case == "stats":
        loaded["chamber"] = _task(ds, relabel=_relabel(ds, stats_sha256_hashed="0" * 64))
    elif case == "yaml":
        loaded["whole"] = _task(ds, relabel=_relabel(ds, yaml_sha256="0" * 64))
    elif case == "route":
        loaded["1"] = _task(ds, lib=True, relabel=_relabel(ds, route="none"))
    elif case == "authors":
        loaded["2"] = _task(ds, lib=True, relabel=_relabel(ds, authors_labels_read=True))
    elif case == "checks":
        checks = {"reasons": ["the settings file's sha256 is x"]}
    elif case == "source":
        source = ["the source imports ['perm_read']"]
    elif case == "fraction":
        loaded["chamber"]["fraction"] = 0.2 + 1e-9
    elif case == "x2_fraction":
        loaded["3"]["fractions"]["m1rc_gmom_x2"] = 0.2
    elif case == "placebo_is_candidate":
        loaded["whole"]["directions"]["placebo_sha256"]["placebo_x2_q2"]["CM"] = "CM" * 4
    elif case == "placebo_pair":
        loaded["whole"]["directions"]["placebo_sha256"]["placebo_x1"] = {g: f"{g}9" for g in WITH_EPI}
    elif case == "permutations_shared":
        loaded["2"]["directions"]["placebo_sha256"]["placebo_x2_q3"]["Epicardium"] = "Epicardium1"
    elif case == "placebo_missing":
        del loaded["2"]["directions"]["placebo_sha256"]["placebo_x2_q1"]
    elif case == "x2_is_x1":
        for row in loaded["1"]["rows"]:
            row["sha256"]["m1rc_gmom_x2"] = row["sha256"]["m1rc_gmom"]
    elif case == "placebo_is_x2":
        for row in loaded["chamber"]["rows"]:
            row["sha256"]["placebo_x2"] = row["sha256"]["m1rc_gmom_x2"]
    r = ds.t63_reading(loaded, committed, checks, source)
    assert r["outcome"] == "VOID" and r["reasons"] and r["follows"]["key"] == ds.VOID
    lines = ds.outcome_lines(r)
    assert lines[0].startswith("Outcome: VOID") and lines[1] == r["follows"]["line"]


def test_the_gates_within_the_tolerance_pass_and_name_the_largest_difference(ds):
    committed = _committed()
    committed["2"]["none_gmom"][4] += 5e-10
    committed["3"]["m1rc_gmom"][1] -= 5e-10
    assert _read(ds, _loaded(ds), committed)["outcome"] == "read"
    committed["3"]["m1rc_gmom"][1] -= 0.25
    r = _read(ds, _loaded(ds), committed)
    assert any("3: m1rc_gmom differs from row T62's committed m1rc_gmom by 0.25" in x for x in r["reasons"])


def test_committed_t62_reads_the_gated_arms(ds, tmp_path):
    (tmp_path / "lib1.json").write_text(json.dumps({"rows": [
        {"seed": s, "arms": {"none_gmom": {"score": 60.0 + s}, "m1rc_gmom": {"score": 61.0 + s},
                             "m1rc_gmom_epi_held": {"score": 0.0}}} for s in range(10)]}))
    got = ds.committed_t62("lib1.json", tmp_path)
    assert got == {"none_gmom": {s: 60.0 + s for s in range(10)}, "m1rc_gmom": {s: 61.0 + s for s in range(10)}}
    assert ds.T62_TASKS.parts[-2:] == ("t62_recomputed", "tasks")


def test_the_read_prints_no_score_and_writes_its_files_once(ds, tmp_path, monkeypatch, capsys):
    tries = tmp_path / "TRIES.md"
    tries.write_text("| T63 (one row) | T1 | x |\n")
    monkeypatch.setattr(ds, "TRIES", tries)
    out = tmp_path / "results" / "t63"
    (out / "tasks").mkdir(parents=True)
    loaded = _loaded(ds, chamber={"m1rc_gmom_x2": 1.234567})
    for k, f in ds.task_files().items():
        (out / "tasks" / f).write_text(json.dumps(loaded[k]))
    committed = _committed()
    monkeypatch.setattr(ds, "committed_t62", lambda f: committed[{v: k for k, v in ds.task_files().items()}[f]])
    monkeypatch.setattr(ds.rc, "relabel_checks", lambda **k: {"reasons": [], "passes": True})
    monkeypatch.setattr(ds.gv, "stamp", lambda started, files: {"commit": "e" * 40, "tree_clean": True, "utc": "x"})
    monkeypatch.setattr(ds.gv, "tree_state", lambda: {"commit": "e" * 40, "clean": True})
    rel = tmp_path / "reliability.json"
    rel.write_text(json.dumps({"groups": {"Epicardium": {"own_direction": True, "cells": {"E9.5": {"A": 6, "B": 5}},
                                                         "spearman_all": 0.35, "eb_factor_mean": 0.745,
                                                         "deciles": [{"spearman": 0.1}]}}}))
    monkeypatch.setattr(ds, "T62_RELIABILITY", rel)
    res = ds.read(out)
    printed = capsys.readouterr().out
    assert res["reading"]["S1"]["holds"] and res["tries"] == 45
    assert "1.23" not in printed and "51." not in printed and "+0." not in printed and "S1:" in printed
    assert "What follows" in printed
    text = (out / "t63.txt").read_text()
    assert "+1.23" in text and "J6 range" in text and "Mean scores" in text and "T1 368 to 413" in text
    assert "spearman_all" in text and "placebo_x2_q4 - none_gmom" in text
    with pytest.raises(SystemExit):
        ds.read(out)  # no second read into the same folder


def test_reliability_beside_reads_the_epicardium_group_only(ds, tmp_path):
    rel = tmp_path / "reliability.json"
    rel.write_text(json.dumps({"groups": {"Epicardium": {"own_direction": True, "cells": {}, "spearman_all": 0.35,
                                                         "eb_factor_mean": 0.745, "deciles": [{"spearman": 0.2}]},
                                          "CM": {"own_direction": True}}}))
    got = ds.reliability_beside(rel)
    assert got == {"own_direction": True, "cells": {}, "spearman_all": 0.35, "eb_factor_mean": 0.745,
                   "decile_spearman": [0.2]}
    assert ds.reliability_beside(tmp_path / "missing.json") is None
