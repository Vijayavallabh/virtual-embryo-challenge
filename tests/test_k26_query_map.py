"""Tests for experiments/t3_decomp/k26_query_map.py (row K26: the query-gene-conditioned response map's kill test).
Synthetic arrays, hand-made unit files and task files in temporary folders outside the repository; stubbed builders.
No competition file, no outside table and no score of any real fold is read (K9's committed x2 file is hashed for its
pin only)."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from vec import harness
from vec.adoption import SEEDS

REPO = Path(__file__).resolve().parents[1]
T3 = REPO / "experiments" / "t3_decomp"
T3X = REPO / "experiments" / "t3_ext"
for p in (T3, T3X):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))
import context_map as cmap  # noqa: E402
import k18_built_files as k18  # noqa: E402
import k24_wider_tables as k24  # noqa: E402
import k26_query_map as q  # noqa: E402
import t3x_tables as tb  # noqa: E402

COMMIT = "c" * 40
ROW = ("| K26 (one row, two scripts) | T3 | `experiments/t3_ext/context_map.py` and "
       "`experiments/t3_decomp/k26_query_map.py` |\n")


@pytest.fixture(autouse=True)
def clean_state():
    yield
    q._K26.clear()
    q._X2.clear()


def _sha(p: Path) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


# ------------------------------------------------------------------------------------------------ constants and row


def test_the_rows_constants():
    assert (q.ROW_LABEL, q.STAGE0_SEED, q.PERM_SEED, q.MIN_MOVED, q.MAX_CORR) == ("K26", 20261122, 20261113, 100, 0.9)
    assert (cmap.RANK, cmap.LAMBDA_PER, cmap.N_RANDOM, cmap.DETECTION_FLOOR, cmap.MIN_TYPE_CELLS) == (2, 10.0, 20,
                                                                                                    0.01, 100)
    assert q.QUERY == {"g": "Gata4", "x2": "Mab21l2"} and q.TRY_ARMS == ("q_map", "q_rand", "q_perm")
    assert q.CELLS == (("plain", "g1"), ("plain", "g2"), ("dm", "g1"), ("dm", "g2"), ("x2", "x2"))
    assert k24.S_STAR == 0.04949423311987974
    assert (q.MIN_GAIN, q.MIN_AHEAD, q.MAX_EXTRA_GATED, q.MAX_TRIES) == (0.3, 8, 1, 3)


def test_importing_the_script_registers_no_harness_candidate():
    assert not set(q.TRY_ARMS + (q.BASE_X,)) & set(harness.CANDIDATES)


def test_k9s_committed_x2_file_is_at_its_pin():
    assert _sha(q.K9_X2) == q.K9_X2_SHA256
    rows = q.k9_rows()
    assert [r["seed"] for r in rows] == list(SEEDS) and all(q.K9_X2_ARM in r["scores"] for r in rows)
    with pytest.raises(q.Refused, match="pinned sha256"):
        q.k9_rows(pin="0" * 64)


def test_the_row_must_be_k26_and_name_the_script(tmp_path, monkeypatch):
    f = tmp_path / "TRIES.md"
    monkeypatch.setattr(q, "TRIES_MD", f)
    f.write_text("| K25 | T3 | other |\n" + ROW)
    assert q.row_present() == "K26"
    for text in ("| K25 | T3 | nothing |\n", ROW.replace("K26", "K27"), ROW + ROW.replace("K26", "K27")):
        f.write_text(text)
        with pytest.raises(SystemExit) as e:
            q.row_present()
        assert e.value.code == q.EXIT_REFUSED


# ------------------------------------------------------------------------------------------------ labels


GENES = [f"g{i}" for i in range(12)]


def _ref():
    rng = np.random.default_rng(0)
    mu = {g: rng.dirichlet(np.ones(len(GENES))) for g in ("A", "B")}
    return {"genes": list(GENES), "groups": ["A", "B"], "shares": {"A": 0.6, "B": 0.4}, "mu": mu,
            "count": {"A": 6, "B": 4}}


def _unit(processed: Path, acc: str, gene: str, L: dict, refsha: str, genes=GENES) -> str:
    d = processed / acc
    d.mkdir(parents=True, exist_ok=True)
    p = d / f"{q.unit_name(acc, gene)}.unit.json"
    p.write_text(json.dumps({"accession": acc, "gene": gene, "genes": list(genes), "refstats_sha256": refsha,
                             "L": {g: [None if not np.isfinite(v) else float(v) for v in x] for g, x in L.items()}}))
    return _sha(p)


DATASETS = {"GSE1": {"status": "cleared"}, "GSE2": {"status": "held_phenocopy"},
            "GSE3": {"status": "excluded_design"}, "BULK": {"status": "cleared", "modality": "bulk"}}


def test_unit_tables_finish_each_unit_as_the_pipeline_and_refuse_held_sets_by_status(tmp_path):
    ref = _ref()
    rng = np.random.default_rng(1)
    La = {"A": rng.normal(0, 0.3, 12), "B": rng.normal(0, 0.3, 12)}
    La["A"][3] = np.nan
    Lb = {"__whole__": rng.normal(0, 0.3, 12)}
    man = {"refstats_sha256": "R", "unit_sha256": {"GSE1:Dnmt3a": _unit(tmp_path, "GSE1", "Dnmt3a", La, "R")}}
    out = q.unit_tables(man, "sc", tmp_path, ref, GENES, DATASETS, frozenset())
    np.testing.assert_allclose(out["Dnmt3a"][0], np.log(tb.finish_pooled(La, ref)[0]))
    bman = {"refstats_sha256": "R", "unit_sha256": {"BULK:Kmt2d": _unit(tmp_path, "BULK", "Kmt2d", Lb, "R")}}
    bout = q.unit_tables(bman, "bulk", tmp_path, ref, GENES, DATASETS, frozenset())
    np.testing.assert_allclose(bout["Kmt2d"][0], np.log(tb.finish(Lb["__whole__"], tb.whole_den(ref))))
    for acc in ("GSE2", "GSE3", "GSE9"):
        bad = {"refstats_sha256": "R", "unit_sha256": {f"{acc}:X": _unit(tmp_path, acc, "X", La, "R")}}
        with pytest.raises(q.Refused, match="status"):
            q.unit_tables(bad, "sc", tmp_path, ref, GENES, DATASETS, frozenset())
    with pytest.raises(q.Refused, match="status"):  # NEVER_FETCH refuses whatever the status says
        q.unit_tables(man, "sc", tmp_path, ref, GENES, DATASETS, frozenset({"GSE1"}))
    with pytest.raises(q.Refused, match="family"):  # a bulk dataset's unit under the single-cell manifest
        q.unit_tables(bman, "sc", tmp_path, ref, GENES, DATASETS, frozenset())
    wrong = {"refstats_sha256": "R", "unit_sha256": {"GSE1:Dnmt3a": "0" * 64}}
    with pytest.raises(q.Refused, match="sha256"):
        q.unit_tables(wrong, "sc", tmp_path, ref, GENES, DATASETS, frozenset())
    other_ref = {"refstats_sha256": "S", "unit_sha256": man["unit_sha256"]}
    with pytest.raises(q.Refused, match="another gene list, reference"):
        q.unit_tables(other_ref, "sc", tmp_path, ref, GENES, DATASETS, frozenset())


def test_outside_labels_average_units_inside_a_family_then_across_families(tmp_path, monkeypatch):
    ref = _ref()
    rng = np.random.default_rng(2)
    tables, processed = tmp_path / "tables", tmp_path / "processed"
    rs = tmp_path / "ref" / "refstats.npz"
    tb.save_ref_stats(ref, rs, "w")
    rsha = _sha(rs)
    L1 = {"A": rng.normal(0, 0.3, 12)}
    L2 = {"A": rng.normal(0, 0.3, 12), "B": rng.normal(0, 0.3, 12)}
    Lb = {"__whole__": rng.normal(0, 0.3, 12)}
    sc_units = {"GSE1:Dnmt3a": _unit(processed, "GSE1", "Dnmt3a", L1, rsha),
                "GSE1:Usb1": _unit(processed, "GSE1", "Usb1", L2, rsha)}
    # a second unit of Dnmt3a from another cleared dataset
    ds = dict(DATASETS, GSE4={"status": "cleared"})
    sc_units["GSE4:Dnmt3a"] = _unit(processed, "GSE4", "Dnmt3a", L2, rsha)
    bulk_units = {"BULK:Dnmt3a": _unit(processed, "BULK", "Dnmt3a", Lb, rsha)}
    recs = {}
    for name, units in (("all", sc_units), ("bulk_all", bulk_units)):
        (tables / name).mkdir(parents=True)
        mp = tables / name / f"{name}_manifest.json"
        mp.write_text(json.dumps({"unit_sha256": units, "refstats_sha256": rsha}))
        recs[name] = _sha(mp)
    loaded = {"listing_sha256": "L", "tables": {"sc3": (np.ones(12), {"manifest_sha256": recs["all"], "sha256": "a"}),
                                                "bulk": (np.ones(12), {"manifest_sha256": recs["bulk_all"],
                                                                       "sha256": "b"}),
                                                "two": (np.ones(12), {"sha256": "t"})}}
    monkeypatch.setattr(k24, "load_tables", lambda genes, env=None: loaded)
    monkeypatch.setattr(q._cm(), "load_datasets", lambda path=None: ds)
    env = {k24.ENV_TABLES: str(tables), k24.ENV_TWO: str(tables), q.ENV_PROCESSED: str(processed),
           q.ENV_REFSTATS: str(rs)}
    out, rec = q.outside_labels(GENES, env)
    sc_d = np.mean([np.log(tb.finish_pooled(L1, ref)[0]), np.log(tb.finish_pooled(L2, ref)[0])], 0)
    bk_d = np.log(tb.finish(Lb["__whole__"], tb.whole_den(ref)))
    np.testing.assert_allclose(out["Dnmt3a"], (sc_d + bk_d) / 2)
    np.testing.assert_allclose(out["Usb1"], np.log(tb.finish_pooled(L2, ref)[0]))
    assert rec["genes"] == ["Dnmt3a", "Usb1"]
    (tables / "all" / "all_manifest.json").write_text("{}")  # not the pinned table's manifest
    with pytest.raises(q.Refused, match="manifest"):
        q.outside_labels(GENES, env)


# ------------------------------------------------------------------------------------------------ features plan


def test_the_plan_keeps_outside_genes_with_a_context_and_draws_no_label_gene():
    det = np.linspace(0.02, 1.0, 400)
    names = [f"r{i}" for i in range(400)]
    names[10], names[200], names[399], names[50] = "Gata4", "Mab21l2", "Dnmt3a", "Usb1"
    dec = cmap.decile_of(det, names)
    plan = q.plan_features(["Dnmt3a", "Usb1", "Kmt2d", "Gata4"], dec)
    assert plan["outside_kept"] == ["Dnmt3a", "Usb1"]
    assert set(plan["outside_left_out"]) == {"Kmt2d", "Gata4"}
    drawn = {g for v in plan["stage0_random"].values() for g in v} | set(plan["q_rand"].values())
    assert not drawn & {"Dnmt3a", "Usb1", "Kmt2d", "Gata4", "Mab21l2"}
    assert all(len(v) == 20 and all(dec[g] == dec[h] for g in v) for h, v in plan["stage0_random"].items())
    assert dec[plan["q_rand"]["g"]] == dec["Gata4"] and dec[plan["q_rand"]["x2"]] == dec["Mab21l2"]
    assert plan == q.plan_features(["Dnmt3a", "Usb1", "Kmt2d", "Gata4"], dec)  # the seeds fix every draw
    assert set(plan["need"]) == {"Dnmt3a", "Usb1", "Gata4", "Mab21l2"} | drawn
    dec.pop("Mab21l2")
    with pytest.raises(cmap.ContextRefused, match="Mab21l2"):
        q.plan_features(["Dnmt3a"], dec)


# ------------------------------------------------------------------------------------------------ the pre-gate


def _pg_inputs(G=150, seed=3):
    rng = np.random.default_rng(seed)
    p = 60
    A = rng.normal(0, 1, (4, p))
    names = ["Gata4", "Mab21l2", "rg", "rx", *[f"o{i}" for i in range(8)]]
    X = rng.normal(0, 1, (len(names), 4)) @ A + 0.1 * rng.normal(0, 1, (len(names), p))
    ctx = {n: (x - x.mean()) / x.std() for n, x in zip(names, X, strict=True)}
    W = rng.normal(0, 1, (p, G)) / p
    outside = {f"o{i}": ctx[f"o{i}"] @ W for i in range(8)}
    frozen = {"Gata4": ctx["Gata4"] @ W, "Mab21l2": ctx["Mab21l2"] @ W}
    closed = {k: np.exp(rng.normal(0, 0.2, G)) for k in q.CLOSED}
    phi_m = np.exp(rng.normal(0, 0.3, G))
    return ctx, outside, frozen, closed, phi_m


def test_the_fold_training_leaves_out_the_folds_own_gene():
    _, outside, frozen, _, _ = _pg_inputs()
    assert sorted(q.fold_training("g", outside, frozen)) == sorted([*outside, "Mab21l2"])
    assert sorted(q.fold_training("x2", outside, frozen)) == sorted([*outside, "Gata4"])
    with pytest.raises(q.Refused, match="own gene"):
        q.fold_training("g", {**outside, "Gata4": frozen["Gata4"]}, frozen)


def test_the_pregate_keeps_a_fresh_table_sized_to_s_star():
    ctx, outside, frozen, closed, phi_m = _pg_inputs()
    pre = q.pregate_all(ctx, outside, frozen, {"g": "rg", "x2": "rx"}, closed, phi_m)
    assert pre["runs"] and pre["tries_if_run"] == 3 and not pre["dropped"]
    for f in q.FOLDS:
        v = pre["folds"][f]
        assert float(np.std(np.log(v["table"]))) == pytest.approx(k24.S_STAR, rel=1e-12)
        assert float(np.std(np.log(v["table_rand"]))) == pytest.approx(k24.S_STAR, rel=1e-12)
        assert v["figures"]["genes_ne_1"] >= q.MIN_MOVED and set(v["corr"]) == {"two", "sc3", "bulk", "mean",
                                                                                 "barred"}
    assert pre["folds"]["g"]["train"] == sorted([*outside, "Mab21l2"]) and pre["folds"]["g"]["rand_gene"] == "rg"


def test_the_pregate_drops_a_disguised_closed_table_and_with_it_the_whole_arm():
    ctx, outside, frozen, closed, phi_m = _pg_inputs()
    names = sorted(q.fold_training("g", outside, frozen))
    m = cmap.fit_map(np.vstack([ctx[g] for g in names]), np.vstack([q.fold_training("g", outside, frozen)[g]
                                                                    for g in names]))
    copy = np.exp(m.predict(ctx["Gata4"]) * 3.0)  # the same table at another size: corr 1
    for key in ("sc3", "barred", "mean"):
        c = dict(closed)
        pm = phi_m
        if key == "barred":
            pm = copy
        elif key == "mean":
            c = {k: copy for k in q.CLOSED}
        else:
            c["sc3"] = copy
        pre = q.pregate_all(ctx, outside, frozen, {"g": "rg", "x2": "rx"}, c, pm)
        assert not pre["runs"] and pre["tries_if_run"] == 0 and "g" in pre["dropped"]
        assert any(key in w for w in pre["dropped"]["g"])


def test_the_pregate_drops_a_table_that_moves_too_few_genes_or_has_no_size():
    ctx, outside, frozen, closed, phi_m = _pg_inputs(G=60)
    pre = q.pregate_all(ctx, outside, frozen, {"g": "rg", "x2": "rx"}, closed, phi_m)
    assert not pre["runs"] and any("fewer than 100" in w for w in pre["dropped"]["g"])
    ctx, outside, frozen, closed, phi_m = _pg_inputs()
    zero = {k: np.zeros_like(v) for k, v in outside.items()}
    fz = {k: np.zeros_like(v) for k, v in frozen.items()}
    pre = q.pregate_all(ctx, zero, fz, {"g": "rg", "x2": "rx"}, closed, phi_m)
    assert not pre["runs"] and any("no size exponent" in w for w in pre["dropped"]["x2"])


# ------------------------------------------------------------------------------------------------ the arms


def _load_tables(genes, rng):
    q._K26.update(genes=list(genes), pregate_sha256="P",
                  tables={f: {q.ARM: np.exp(rng.normal(0, 0.05, len(genes))),
                              q.RAND: np.exp(rng.normal(0, 0.05, len(genes)))} for f in q.FOLDS})


def test_fold_tables_follow_the_references_order_and_permute_q_map_within_deciles():
    rng = np.random.default_rng(4)
    genes = [f"g{i}" for i in range(40)]
    _load_tables(genes, rng)
    ref_genes = genes[::-1]
    basis = rng.uniform(0, 1, 40)
    t = q.fold_tables("x2", ref_genes, basis)
    np.testing.assert_array_equal(t[q.ARM], q._K26["tables"]["x2"][q.ARM][::-1])
    np.testing.assert_array_equal(t[q.PERM], harness.t3_decile_permutation(t[q.ARM], basis, 20261113, 10))
    assert not np.array_equal(t[q.PERM], t[q.ARM])


def _stub_bs(calls):
    t3d = SimpleNamespace(stable_rng=lambda seed, name: (seed, name))

    def prediction_inputs(ref, spec, extra, phi):
        return {"ratio": None, "alphas": None, "w": 1.0, "phi": phi, "f": "central"}, {}

    def construct(ref, n, rng, **kw):
        calls.append({"n": n, "rng": rng, **kw})
        return "built"

    return SimpleNamespace(t3d=t3d, prediction_inputs=prediction_inputs, construct=construct)


def test_the_g_arms_pass_the_barred_factor_times_their_table_on_the_one_stream():
    rng = np.random.default_rng(5)
    phi = np.exp(rng.normal(0, 0.3, 10))
    T = {a: np.exp(rng.normal(0, 0.05, 10)) for a in q.TRY_ARMS}
    F = {"barred_spec": {}, "phi_barred": phi, "k26": T}
    calls = []
    bs = _stub_bs(calls)
    for arm in q.G_ARMS:
        q.build_arm(bs, None, arm, 9, 4, F)
    assert all(c["rng"] == (4, q.STREAM) and c["n"] == 9 and c["f"] == "central" for c in calls)
    np.testing.assert_array_equal(calls[0]["phi"], phi)
    for c, a in zip(calls[1:], q.TRY_ARMS, strict=True):
        np.testing.assert_allclose(c["phi"], phi * T[a], rtol=0, atol=0)
    with pytest.raises(RuntimeError):
        q.build_arm(bs, None, "x_sc3", 9, 4, F)


def test_k18_arms_are_restored_after_the_context():
    before = (k18.ARMS, k18.build_arm, k18.fold_factors)
    with q.k18_arms(q.G_ARMS):
        assert k18.ARMS == q.G_ARMS and k18.build_arm is q.build_arm
    assert (k18.ARMS, k18.build_arm, k18.fold_factors) == before
    with pytest.raises(RuntimeError):
        with q.k18_arms((q.ARM,)):
            pass


def test_x2_base_is_k2as_build_and_a_try_arm_multiplies_k2as_factor(monkeypatch):
    rng = np.random.default_rng(6)
    genes = [f"g{i}" for i in range(20)]
    _load_tables(genes, rng)
    ref_genes = genes[::-1]
    phi = np.exp(rng.normal(0, 0.3, 20))
    basis = rng.uniform(0, 1, 20)
    calls = []
    bs = _stub_bs(calls)
    bs.assemble = lambda ref, built, gl, rng_: ("X", "coords", rng_)
    prep = {"phi": phi, "ref": SimpleNamespace(genes=ref_genes, V=SimpleNamespace(p="own"))}
    pool = SimpleNamespace(n_obs=33)
    ctx = SimpleNamespace(inputs={"wt": pool}, task=SimpleNamespace(reference_role="wt"), seed=2,
                          board=SimpleNamespace(genes=["a"], needs_coords=False))
    monkeypatch.setattr(harness, "_t3_factor_xko", lambda c: "k2a")
    monkeypatch.setattr(harness, "_k2_prepared", lambda c: prep)
    monkeypatch.setattr(harness, "_t3_head", lambda: bs)
    monkeypatch.setattr(harness, "count_for", lambda board, n: n + 1)
    monkeypatch.setattr(harness, "t3_mean_expression", lambda a: basis)
    monkeypatch.setattr(harness, "finish", lambda c, p, X, coords: (X, coords))
    assert q.build_x2(ctx, q.BASE_X) == "k2a"
    for arm in q.TRY_ARMS:
        assert q.build_x2(ctx, arm) == ("X", None)
    t = q.fold_tables("x2", ref_genes, basis)
    for c, a in zip(calls, q.TRY_ARMS, strict=True):
        assert c["rng"] == (2, harness.K2_STREAM) and c["n"] == 34 and c["f"] == "own"
        np.testing.assert_allclose(c["phi"], phi * t[a], rtol=0, atol=0)


# ------------------------------------------------------------------------------------------------ the gate and the line


def _rows(cell, arms, base_score, shift, seeds=SEEDS, de=0.01, gated=()):
    out = []
    for s in seeds:
        sc = {a: base_score + 0.01 * s + shift.get(a, 0.0) for a in arms}
        if cell[0] != "dm":
            sc[harness.FLOOR] = 50.0
        sk = {a: {"de_score": 0.1 + (de if a == q.ARM else 0.0), "de_direction": 0.2, "mmd_u": 0.3} for a in sc}
        t3 = {a: {"severity": {"state": ("beta_gate" if (a, s) in gated else "scored"), "beta": 1.0, "r2": 0.5}}
              for a in sc}
        out.append({"seed": s, "scores": sc, "skills": sk, "t3": t3})
    return out


def test_the_reproduction_gate_holds_within_k17s_tolerance_and_names_the_failed_base():
    refs = {"k18": {c: _rows(c, (q.BASE,), 55.0, {}) for c in q.CELLS if c[0] != "x2"},
            "k9_x2": [{"seed": r["seed"], "scores": {"raw": r["scores"][q.BASE_X]}}
                      for r in _rows(("x2", "x2"), (q.BASE_X,), 51.0, {})]}
    rows = {c: _rows(c, (q.base_of(c),), 51.0 if c[0] == "x2" else 55.0, {}) for c in q.CELLS}
    assert q.guard_eval(rows, refs)["passes"]
    rows[("x2", "x2")] = _rows(("x2", "x2"), (q.BASE_X,), 51.2, {})
    g = q.guard_eval(rows, refs)
    assert not g["passes"] and g["failed"] == [q.BASE_X]
    rows[("x2", "x2")] = _rows(("x2", "x2"), (q.BASE_X,), 51.0, {})
    rows[("dm", "g2")] = rows[("dm", "g2")][:9]  # a seed missing
    assert q.guard_eval(rows, refs)["failed"] == [q.BASE]


def _per(shifts: dict, de=0.01, gated=None):
    """summarise() per cell from rows where each arm is the base plus shifts[cell][arm]."""
    out = {}
    for c in q.CELLS:
        arms = q.X_ARMS if c[0] == "x2" else q.G_ARMS
        rows = _rows(c, arms, 50.0, shifts.get(c, {}), de=de, gated=(gated or {}).get(c, ()))
        full = (harness.FLOOR, *arms) if c[0] != "dm" else arms
        out[c] = q.summarise(rows, full, q.base_of(c))
    return out


GOOD = {c: {q.ARM: 1.0} for c in q.CELLS}


def test_the_pass_line_passes_when_every_clause_holds():
    d = q.decide(_per(GOOD))
    assert d["passes"] and all(d[k] for k in "abcde") and len(d["ahead_cells"]) == 5
    # a gain on x2 and one g-cell in each world is enough for (a)
    one = {("x2", "x2"): {q.ARM: 1.0}, ("plain", "g2"): {q.ARM: 1.0}, ("dm", "g1"): {q.ARM: 1.0}}
    assert q.decide(_per(one))["passes"]


def test_each_clause_of_the_pass_line_fails_alone():
    no_x2 = {c: v for c, v in GOOD.items() if c[0] != "x2"}
    d = q.decide(_per(no_x2))
    assert not d["passes"] and not d["a"]
    no_dm = {c: v for c, v in GOOD.items() if c[0] != "dm"}
    assert not q.decide(_per(no_dm))["a"]
    behind = {**GOOD, ("dm", "g2"): {q.ARM: -0.5}}
    d = q.decide(_per(behind))
    assert d["a"] and not d["b"] and not d["passes"]
    gated = {("plain", "g1"): {(q.ARM, s) for s in SEEDS[:2]}}
    d = q.decide(_per(GOOD, gated=gated))
    assert not d["c"] and not d["passes"]
    rand_close = {c: {q.ARM: 1.0, q.RAND: 0.8} for c in q.CELLS}
    d = q.decide(_per(rand_close))
    assert d["a"] and not d["d"] and not d["passes"]
    perm_close = {c: {q.ARM: 1.0, q.PERM: 0.75} for c in q.CELLS}
    assert not q.decide(_per(perm_close))["d"]
    d = q.decide(_per(GOOD, de=-0.01))
    assert not d["e"] and not d["passes"]


def test_ahead_needs_a_mean_above_three_tenths_and_eight_of_ten_seeds():
    v = {"mean": 0.31, "ahead": 8}
    assert q._beats(v) and not q._beats({"mean": 0.3, "ahead": 10}) and not q._beats({"mean": 2.0, "ahead": 7})


# ------------------------------------------------------------------------------------------------ the reading


def _prov(commit=COMMIT, machine="m", clean=True):
    return {"commit": commit, "tree_clean": clean, "machine": machine, "utc": "t"}


def _run_folder(tmp_path, shifts=GOOD, guard_pass=True, commit_of=None):
    R = tmp_path / "runs"
    R.mkdir(parents=True)
    k18.write_task_file(R / "features.json", {"plan": {"outside_left_out": {}}, "context_record": {},
                                               "provenance": _prov()})
    fsha = _sha(R / "features.json")
    k18.write_task_file(R / "stage0.json", {"stage0": {"passes": True, "per_gene": {}}, "features_sha256": fsha,
                                             "provenance": _prov()})
    k18.write_task_file(R / "pregate.json", {"pregate": {"runs": True, "dropped": {}, "folds": {}},
                                              "features_sha256": fsha, "stage0_sha256": _sha(R / "stage0.json"),
                                              "panel": [], "provenance": _prov()})
    psha = _sha(R / "pregate.json")
    refs = {"k18": {}, "k9_x2": []}

    def write(rel, cell, s, arms, shift, prefix=""):
        rows = _rows(cell, arms, 50.0, shift, seeds=[s])
        prov = _prov(commit=(commit_of or {}).get(rel, COMMIT))
        if cell[0] == "dm":
            body = rows[0] | {"task": cell[1], "provenance": prov, "pregate_sha256": psha}
        else:
            key = k18.TASKS[cell[1]].key if cell[0] == "plain" else q.X2_KEY
            body = {"task": {"key": key}, "rows": rows, "provenance": prov, "pregate_sha256": psha}
        p = R / prefix / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        k18.write_task_file(p, body)

    gsha = {}
    for f in q.guard_files():
        rel = f[len("guard/"):]
        cell, s = q.cell_of_file(rel)
        write(rel, cell, s, (q.base_of(cell),), {}, prefix="guard")
        gsha[f] = _sha(R / f)
    for c in q.CELLS:
        base_rows = _rows(c, (q.base_of(c),), 50.0, {})
        if c[0] == "x2":
            refs["k9_x2"] = [{"seed": r["seed"], "scores": {"raw": r["scores"][q.BASE_X]}} for r in base_rows]
        else:
            refs["k18"][c] = base_rows
    k18.write_task_file(R / "guard.json", {"cells": {}, "failed": [] if guard_pass else [q.BASE],
                                            "passes": guard_pass, "files_sha256": gsha, "pregate_sha256": psha,
                                            "provenance": _prov()})
    for f in q.score_files():
        cell, s = q.cell_of_file(f)
        write(f, cell, s, q.X_ARMS if cell[0] == "x2" else q.G_ARMS, shifts.get(cell, {}))
    return R, refs


def test_the_reading_passes_a_folder_that_meets_the_line(tmp_path):
    R, refs = _run_folder(tmp_path)
    res = q.reading(R, refs)
    assert not res["void"] and res["decision"]["passes"] and res["outcome"].startswith("PASS")
    assert res["tries_counted"] == 3 and res["stage0_counted"] == 1
    txt = q.text(res)
    assert "Try K26: PASS" in txt and "q_map - q_rand" in txt


def test_the_reading_fails_a_folder_without_a_gain(tmp_path):
    R, refs = _run_folder(tmp_path, shifts={})
    res = q.reading(R, refs)
    assert not res["void"] and not res["decision"]["passes"] and res["outcome"].startswith("FAIL")


def test_a_failed_reproduction_gate_is_void_and_counts_the_failed_base(tmp_path):
    R, refs = _run_folder(tmp_path, guard_pass=False)
    res = q.reading(R, refs)
    assert res["void"] and res["tries_counted"] == 1 and "reproduction gate failed" in res["why"][0]


def test_mixed_commits_or_a_base_off_the_guard_part_void_the_reading(tmp_path):
    R, refs = _run_folder(tmp_path / "a", commit_of={"x2/seed_3.json": "d" * 40})
    res = q.reading(R, refs)
    assert res["void"] and any("different commits" in w for w in res["why"])
    R, refs = _run_folder(tmp_path / "b", shifts={**GOOD, ("plain", "g1"): {q.ARM: 1.0, q.BASE: 0.2}})
    res = q.reading(R, refs)
    assert res["void"] and any("off the guard part's" in w for w in res["why"])


def test_a_missing_task_file_is_incomplete_not_a_reading(tmp_path):
    R, refs = _run_folder(tmp_path)
    (R / "dm" / "g2_seed_7.json").unlink()
    res = q.reading(R, refs)
    assert res["incomplete"] and res["missing"] == ["dm/g2_seed_7.json"]


def test_a_stage_file_made_from_another_features_file_voids_the_reading(tmp_path):
    R, refs = _run_folder(tmp_path)
    body = json.loads((R / "features.json").read_text())
    body["plan"]["outside_left_out"] = {"x": "changed"}
    (R / "features.json").write_text(json.dumps(body))
    res = q.reading(R, refs)
    assert res["void"] and any("another features.json" in w for w in res["why"])


def test_the_read_command_never_overwrites_a_reading(tmp_path, monkeypatch, capsys):
    R, refs = _run_folder(tmp_path)
    monkeypatch.setattr(q, "row_present", lambda path=None: "K26")
    monkeypatch.setattr(k18, "stamp", lambda: _prov())
    monkeypatch.setattr(q, "guard_refs", lambda *a, **k: refs)
    out = tmp_path / "reading"
    q.main(["read", "--runs", str(R), "--out", str(out), "--commit", COMMIT])
    assert (out / "reading.json").is_file() and (out / "reading.txt").is_file()
    printed = capsys.readouterr().out
    assert printed.strip().startswith("wrote ") and "reading.txt" in printed
    with pytest.raises(SystemExit) as e:
        q.main(["read", "--runs", str(R), "--out", str(out), "--commit", COMMIT])
    assert e.value.code == q.EXIT_EXISTS
    with pytest.raises(SystemExit) as e:  # another HEAD
        q.main(["read", "--runs", str(R), "--out", str(tmp_path / "r2"), "--commit", "d" * 40])
    assert e.value.code == q.EXIT_REFUSED


# ------------------------------------------------------------------------------------------------ step order


def test_steps_refuse_out_of_order(tmp_path, monkeypatch):
    monkeypatch.setattr(q, "preconditions", lambda commit: ("K26", _prov()))
    R = tmp_path / "runs"
    R.mkdir()
    k18.write_task_file(R / "features.json", {"plan": {}, "panel": [], "provenance": _prov()})
    # stage 0 failed: the pre-gate stops with exit 7
    k18.write_task_file(R / "stage0.json", {"stage0": {"passes": False}, "features_sha256": _sha(R / "features.json"),
                                             "provenance": _prov()})
    with pytest.raises(SystemExit) as e:
        q.main(["pregate", "--out", str(R), "--commit", COMMIT])
    assert e.value.code == q.EXIT_STAGE0
    with pytest.raises(SystemExit) as e:  # and a run refuses
        q.main(["run", "--world", "x2", "--out", str(R), "--commit", COMMIT])
    assert e.value.code == q.EXIT_REFUSED
    # stage 0 re-read from its file, never recomputed
    with pytest.raises(SystemExit) as e:
        q.main(["stage0", "--out", str(R), "--commit", COMMIT])
    assert e.value.code == q.EXIT_STAGE0
    # a passing stage 0, a pregate.json, but no guard.json: a scoring part refuses
    (R / "stage0.json").unlink()
    k18.write_task_file(R / "stage0.json", {"stage0": {"passes": True}, "features_sha256": _sha(R / "features.json"),
                                             "provenance": _prov()})
    k18.write_task_file(R / "pregate.json", {"pregate": {"runs": True, "folds": {}}, "panel": [],
                                              "provenance": _prov()})
    with pytest.raises(SystemExit) as e:
        q.main(["run", "--world", "plain", "--task", "g1", "--out", str(R), "--commit", COMMIT])
    assert e.value.code == q.EXIT_REFUSED
    # a step file of another commit is refused
    with pytest.raises(SystemExit) as e:
        q.load_step_file(R, "features.json", "d" * 40)
    assert e.value.code == q.EXIT_REFUSED


def test_a_rerun_of_a_done_target_free_step_reads_its_file_again(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(q, "preconditions", lambda commit: ("K26", _prov()))
    R = tmp_path / "runs"
    R.mkdir()
    k18.write_task_file(R / "features.json", {"provenance": _prov()})
    with pytest.raises(SystemExit) as e:
        q.main(["features", "--out", str(R), "--commit", COMMIT])
    assert e.value.code == 0 and capsys.readouterr().out == "features: present\n"
    with pytest.raises(SystemExit) as e:
        q.main(["features", "--out", str(R), "--commit", "d" * 40])
    assert e.value.code == q.EXIT_EXISTS
    k18.write_task_file(R / "pregate.json", {"pregate": {"runs": False}, "provenance": _prov()})
    with pytest.raises(SystemExit) as e:
        q.main(["pregate", "--out", str(R), "--commit", COMMIT])
    assert e.value.code == q.EXIT_DROPPED


def test_the_progress_lines_carry_no_figure(capsys):
    q.progress_log("score", "x2", "x2")("[x2] seed 7: floor 50.00, q_map 51.23")
    q.progress_log("guard", "plain", "g1")("whatever seed 3: 49.9")
    assert capsys.readouterr().out == "K26 x2 x2 seed 7: scored\nK26 guard plain g1 seed 3: scored\n"
