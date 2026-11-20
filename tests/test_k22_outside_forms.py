"""Tests for experiments/t3_decomp/k22_outside_forms.py (try K22: the barred factor times one outside-knockout table).
Offline, synthetic and fast, in the manner of tests/test_k21_factor_forms.py: K18's plain fold fixture with stubbed
construct, composition and assemble, hand-made tables, synthetic task files. No competition data and no outside data."""

from __future__ import annotations

import hashlib
import json
import re
import sys
from pathlib import Path
from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from vec import harness
from vec.adoption import SEEDS

REPO = Path(__file__).resolve().parents[1]
T3 = REPO / "experiments" / "t3_decomp"
if str(T3) not in sys.path:
    sys.path.insert(0, str(T3))
import k18_built_files as k18  # noqa: E402
import k21_factor_forms as k21  # noqa: E402
import k22_outside_forms as k22  # noqa: E402
import outside_forms as ext  # noqa: E402

G, N_POOL, N_DRAW = 14, 90, 150
GENES = [f"g{i}" for i in range(G)]
FILE_GENES = GENES[::-1] + ["extra"]
PINS = {"listing_sha256": "l" * 64, "tables": {"all": "a" * 64, "unrelated": "u" * 64}}


@pytest.fixture(autouse=True)
def pins(monkeypatch):
    monkeypatch.setattr(ext, "K22_LISTING_SHA256", PINS["listing_sha256"])
    monkeypatch.setattr(ext, "K22_TABLE_SHA256", dict(PINS["tables"]))
    yield
    k22._EXT.clear()


# ------------------------------------------------------------------------------------------------ constants


def test_arms_constants_and_the_shared_module():
    assert k22.ARMS == ("barred_form", "x_all", "x_all_perm", "x_unrel", "x_unrel_perm", "rep_all")
    assert k22.BASE == k18.BARRED and k22.STREAM == k18.STREAM == "t3_s2_form"
    assert k22.NAMEABLE == ("x_all", "x_unrel") and k22.CONTROL_OF == {"x_all": "x_all_perm",
                                                                       "x_unrel": "x_unrel_perm"}
    assert (k22.MIN_GAIN, k22.MIN_AHEAD, k22.MAX_EXTRA_GATED, k22.GUARD_MEAN, k22.GUARD_SEED) == (0.3, 8, 1, 0.1, 0.5)
    assert k22.ext is ext and k22.TASKS == k18.TASKS and k22.WORLDS == k18.WORLDS
    assert k22.K18_SHA256 == k21.K18_SHA256 and k22.RESULTS == T3 / "results" / "k22_outside_forms"
    assert {t for a, t in k22.TABLE_OF.items()} == set(ext.TRIES)


def test_the_run_refuses_without_its_row(tmp_path, capsys):
    p = tmp_path / "TRIES.md"
    p.write_text("| K21 (one row) | T3 | `experiments/t3_decomp/k21_factor_forms.py` |\n")
    with pytest.raises(SystemExit) as e:
        k22.row_present(p)
    assert e.value.code == 2 and "no row" in capsys.readouterr().err
    p.write_text(p.read_text() + "| K22 (one row) | T3 | `experiments/t3_decomp/k22_outside_forms.py` with ... |\n")
    assert k22.row_present(p) == "K22"
    p.write_text(p.read_text() + "| K23 (one row) | T3 | `experiments/t3_decomp/k22_outside_forms.py` again |\n")
    with pytest.raises(SystemExit) as e:
        k22.row_present(p)
    assert e.value.code == 2
    for argv in (["run", "--world", "plain", "--task", "g1", "--out", str(tmp_path / "o")],
                 ["pregate", "--out", str(tmp_path / "o")]):
        with pytest.MonkeyPatch.context() as m:
            m.setattr(k22, "TRIES_MD", tmp_path / "none.md")
            (tmp_path / "none.md").write_text("no rows\n")
            with pytest.raises(SystemExit) as e:
                k22.main(argv)
            assert e.value.code == 2
    assert not (tmp_path / "o").exists()


def test_a_run_refuses_unfilled_pins_an_untracked_listing_and_a_dirty_tree(tmp_path, monkeypatch, capsys):
    tries = tmp_path / "TRIES.md"
    tries.write_text("| K22 (one row) | T3 | `experiments/t3_decomp/k22_outside_forms.py` |\n")
    monkeypatch.setattr(k22, "TRIES_MD", tries)
    monkeypatch.setattr(k18, "stamp", lambda: {"commit": "abc", "tree_clean": True, "machine": "m1", "utc": "x"})
    monkeypatch.setattr(ext, "K22_LISTING_SHA256", None)
    argv = ["run", "--world", "plain", "--task", "g1", "--out", str(tmp_path / "o")]
    with pytest.raises(SystemExit) as e:
        k22.main(argv)
    assert e.value.code == 2 and "not filled" in capsys.readouterr().err
    monkeypatch.setattr(ext, "K22_LISTING_SHA256", "0" * 64)
    monkeypatch.setattr(ext, "TABLES_FILE", tmp_path / "TABLES.sha256")
    (tmp_path / "TABLES.sha256").write_text("")
    with pytest.raises(SystemExit) as e:
        k22.main(argv)
    assert e.value.code == 2 and "not tracked" in capsys.readouterr().err
    monkeypatch.setattr(k18, "stamp", lambda: {"commit": "abc", "tree_clean": False, "machine": "m1", "utc": "x"})
    with pytest.raises(SystemExit) as e:
        k22.main(argv)
    assert e.value.code == 2 and "not clean" in capsys.readouterr().err
    assert not (tmp_path / "o").exists()


def test_dropped_tables_drop_their_arms():
    assert k22.active_arms({"dropped": {}}) == k22.ARMS
    assert k22.active_arms({"dropped": {"all": "x"}}) == ("barred_form", "x_unrel", "x_unrel_perm")
    assert k22.active_arms({"dropped": {"unrelated": "x"}}) == ("barred_form", "x_all", "x_all_perm", "rep_all")
    assert k22.active_arms({"dropped": {"all": "x", "unrelated": "y"}}) == ("barred_form",)
    assert set(k22.dropped_arms({"dropped": {"all": "x"}})) == {"x_all", "x_all_perm", "rep_all"}


# ------------------------------------------------------------------------------------------------ the arms


@pytest.fixture
def fold(monkeypatch):
    """K18's plain fold fixture (as tests/test_k21_factor_forms.py): a synthetic pool and reference, the builder's
    construct, composition and assemble stubbed to log what they receive; the barred factor and the tables in a file
    gene order that differs from the reference's."""
    bs = harness._t3_head()
    rng = np.random.default_rng(5)
    C = rng.poisson(3.0, size=(N_POOL, G)).astype(np.int32)
    C[C.sum(1) == 0, 0] = 1
    lab = rng.integers(0, 3, N_POOL)
    P = (C / C.sum(1, keepdims=True)).astype(np.float32)
    V = SimpleNamespace(p=np.array([0.5, 0.3, 0.2]))
    ref = bs.Reference(preset="gata4", rules={}, genes=GENES, V=V, pool=bs.t3d.Pool(C, lab, P, 3),
                       coords=rng.normal(size=(N_POOL, 3)))
    pool = ad.AnnData(X=bs.t3d.lognorm(C).astype(np.float32), var=pd.DataFrame(index=GENES))
    group = np.array(["A", "B", "Unknown"])[rng.integers(0, 3, N_POOL)]
    phi = np.exp(rng.normal(0, 0.5, G))
    per = {"A": np.exp(rng.normal(0, 0.5, G)), "B": np.exp(rng.normal(0, 0.5, G))}
    prep = {"ref": ref, "phi": phi, "per": per, "group": group, "record": {}}
    f_s2 = np.array([0.6, 0.25, 0.15])
    barred_file = np.clip(np.exp(rng.normal(0, 0.8, G + 1)), 0.25, 4.0)
    t_all = np.where(rng.random(G + 1) < 0.6, np.exp(rng.normal(0, 0.5, G + 1)), 1.0)
    t_unrel = np.where(rng.random(G + 1) < 0.6, np.exp(rng.normal(0, 0.5, G + 1)), 1.0)
    spec, _, _ = k18.barred_spec()
    log = {"construct": [], "assemble": []}

    def construct(r, n, g, f=None, phi=None, w=1.0, ratio=None, alphas=None):
        assert r is ref and ratio is None and alphas is None
        rows = g.integers(0, N_POOL, n)
        log["construct"].append({"phi": None if phi is None else np.asarray(phi).copy(), "w": w, "row": rows.copy()})
        Cr = C[rows]
        if phi is not None:
            Cr = bs.t3d.apply_gene_factor(Cr, np.asarray(phi) ** w if w != 1.0 else np.asarray(phi), P[rows], g)
        return bs.Built(X=bs.t3d.lognorm(Cr), lab=lab[rows], source=np.full(n, "wt", dtype=object), row=rows)

    def assemble(r, built, panel, g):
        log["assemble"].append({"row": built.row.copy(), "draw": g.random()})
        return built.X, None, {}

    monkeypatch.setattr(k21.k7, "_prepared", lambda ctx: prep)
    monkeypatch.setattr(bs, "construct", construct)
    monkeypatch.setattr(bs, "composition", lambda v, s, params: (f_s2, {"states_missing": []}))
    monkeypatch.setattr(bs, "assemble", assemble)
    monkeypatch.setattr(harness, "count_for", lambda board, n: N_DRAW)
    monkeypatch.setattr(harness, "finish", lambda ctx, p, X, coords: X)
    monkeypatch.setattr(k18, "_BARRED", {"phi": barred_file, "genes": FILE_GENES, "spec": spec, "record": {}})
    k18._FOLD.clear()
    k22._EXT.clear()
    k22._EXT.update(genes=list(FILE_GENES), T={"all": t_all, "unrelated": t_unrel}, record=dict(PINS))

    def ctx(task: str, seed: int):
        return SimpleNamespace(task=SimpleNamespace(key=k18.TASKS[task].key, reference_role="wt"),
                               inputs={"wt": pool}, board=SimpleNamespace(genes=GENES, needs_coords=False), seed=seed)

    order = lambda v: k18.in_gene_order(v, FILE_GENES, GENES)  # noqa: E731
    yield SimpleNamespace(bs=bs, ctx=ctx, log=log, pool=pool, phi_barred=order(barred_file), t_all=order(t_all),
                          t_unrel=order(t_unrel), basis=harness.t3_mean_expression(pool))
    k18._FOLD.clear()


def test_barred_form_is_k18s_arm_value_for_value(fold):
    for seed in (0, 7):
        k18._FOLD.clear()
        want = k18._build_plain(fold.ctx("g1", seed), k18.BARRED)  # K18 as committed
        with k22.k18_arms(k22.ARMS):
            got = k18._build_plain(fold.ctx("g1", seed), k22.BASE)
        assert np.array_equal(got, want)
    assert k18.ARMS == k22._K18_ORIGINAL["ARMS"] and k18.build_arm is k22._K18_ORIGINAL["build_arm"]
    assert k18.fold_factors is k22._K18_ORIGINAL["fold_factors"]


def test_every_arm_draws_the_barred_s1s_cells_and_applies_its_factor(fold, monkeypatch):
    calls = []
    original = ext.product
    monkeypatch.setattr(ext, "product", lambda *a, **k: calls.append(a[1]) or original(*a, **k))
    want = {"barred_form": fold.phi_barred,
            "x_all": fold.phi_barred * fold.t_all,
            "x_all_perm": fold.phi_barred * harness.t3_decile_permutation(fold.t_all, fold.basis, 20261113, 10),
            "x_unrel": fold.phi_barred * fold.t_unrel,
            "x_unrel_perm": fold.phi_barred * harness.t3_decile_permutation(fold.t_unrel, fold.basis, 20261113, 10),
            "rep_all": fold.t_all}
    X = {}
    with k22.k18_arms(k22.ARMS):
        for arm in k22.ARMS:
            fold.log["construct"].clear()
            fold.log["assemble"].clear()
            X[arm] = k18._build_plain(fold.ctx("g2", 3), arm)
            (c,) = fold.log["construct"]
            assert c["w"] == 1.0 and np.array_equal(c["phi"], want[arm]), arm
            fold.log[arm] = (c["row"], fold.log["assemble"][0])
    assert calls == [1.0] * 4  # the tries and the controls go through the one shared product
    base_rows, base_asm = fold.log[k22.BASE]
    for arm in k22.ARMS:
        rows, asm = fold.log[arm]
        assert np.array_equal(rows, base_rows) and asm["draw"] == base_asm["draw"]
        if arm != k22.BASE:
            assert not np.array_equal(X[arm], X[k22.BASE])
    # the table in the reference's order: an absent gene (1) leaves phi_M as it is
    absent = fold.t_all == 1.0
    assert absent.any() and np.array_equal(want["x_all"][absent], fold.phi_barred[absent])


def test_the_practice_factor_equals_the_builder_s_at_x_1():
    """K22's arm factor and the builder's factor are outside_forms.product on the same vectors; at x = 1 equal."""
    rng = np.random.default_rng(2)
    phi_m = np.clip(np.exp(rng.normal(0, 0.6, 30)), 0.25, 4.0)
    T = np.where(rng.random(30) < 0.5, np.exp(rng.normal(0, 0.5, 30)), 1.0)
    F = {"phi_barred": phi_m, "ext": {"all": T, "unrelated": T}, "ext_perm": {"all": T, "unrelated": T}}
    bo = sys.modules.get("build_outside_form")
    if bo is None:
        import importlib.util

        spec = importlib.util.spec_from_file_location("build_outside_form", T3 / "build_outside_form.py")
        bo = importlib.util.module_from_spec(spec)
        sys.modules["build_outside_form"] = bo
        spec.loader.exec_module(bo)
    for x in (1.0, 1.5, 2.0):
        built, _ = bo.form_factor("xall", x, phi_m, T)
        assert np.array_equal(built, (phi_m if x == 1.0 else phi_m ** x) * T)
    assert np.array_equal(k22.arm_factor("x_all", F), bo.form_factor("xall", 1.0, phi_m, T)[0])
    assert bo.ext is k22.ext is ext


# ------------------------------------------------------------------------------------------------ the rule


def _rows(scores: dict, gated: dict | None = None, noise: dict | None = None) -> list[dict]:
    gated = gated or {}
    out = []
    for s in SEEDS:
        sc = {a: v + 0.1 * s + (noise[a][s] if noise and a in noise else 0.0) for a, v in scores.items()}
        out.append({"seed": s, "scores": sc,
                    "skills": {a: {"de_score": 0.5, "severity_slope": 0.6} for a in sc},
                    "t3": {a: {"severity": {"state": "r2_gate" if s < gated.get(a, 0) else "scored", "beta": 0.4,
                                            "r2": 0.05}} for a in sc}})
    return out


BASE_SCORE = 56.0


def _table(xa=None, xu=None, rep=None, ctl=None) -> dict:
    t = {harness.FLOOR: 50.0, **{a: BASE_SCORE for a in k22.ARMS}}
    for k, v in (("x_all", xa), ("x_unrel", xu), ("rep_all", rep)):
        if v is not None:
            t[k] = v
    t.update(ctl or {})
    return t


def _per(cells: dict, gated: dict | None = None, arms=k22.ARMS) -> dict:
    out: dict = {}
    for w in k22.WORLDS:
        a = (harness.FLOOR, *arms) if w == "plain" else arms
        out[w] = {n: k22._summarise(_rows({x: cells.get((w, n), _table())[x] for x in a},
                                          (gated or {}).get((w, n))), a) for n in k22.TASKS}
    return out


ALL_CELLS = [(w, n) for w in k22.WORLDS for n in k22.TASKS]


def test_qualified_needs_both_worlds_both_replicates_no_loss_the_gates_and_the_control():
    up = _table(xu=57.0)
    one_each = {("plain", "g1"): up, ("dm", "g2"): up}
    res = k22.decide(_per(one_each), k22.ARMS, 0.2)
    assert res["status"]["x_unrel"]["qualified"] and res["names"] == ["x_unrel"]
    assert not res["status"]["x_all"]["qualified"]
    for cells in ({("plain", "g1"): up, ("plain", "g2"): up}, {("plain", "g1"): up, ("dm", "g1"): up}):
        assert not k22.decide(_per(cells), k22.ARMS, 0.2)["status"]["x_unrel"]["qualified"]
    res = k22.decide(_per(one_each | {("dm", "g1"): _table(xu=55.5)}), k22.ARMS, 0.2)
    assert not res["status"]["x_unrel"]["qualified"] and res["status"]["x_unrel"]["behind_cells"] == ["dm g1"]
    res = k22.decide(_per(one_each, {("dm", "g1"): {"x_unrel": 2}}), k22.ARMS, 0.2)
    assert not res["status"]["x_unrel"]["c"] and not res["names"]
    # its control as high where it is ahead (plain): (d) fails
    res = k22.decide(_per({("plain", "g1"): _table(xu=57.0, ctl={"x_unrel_perm": 56.9}), ("dm", "g2"): up}),
                     k22.ARMS, 0.2)
    assert not res["status"]["x_unrel"]["d"] and not res["names"]
    # the dm cell where it is ahead does not beat the control: (d) fails there too
    res = k22.decide(_per({("plain", "g1"): up, ("dm", "g2"): _table(xu=57.0, ctl={"x_unrel_perm": 56.9})}),
                     k22.ARMS, 0.2)
    assert not res["status"]["x_unrel"]["d"]
    # a control or rep_all never qualifies or is named, however high
    res = k22.decide(_per({k: _table(rep=60.0, ctl={"x_all_perm": 60.0, "x_unrel_perm": 60.0}) for k in ALL_CELLS}),
                     k22.ARMS, 0.2)
    assert res["names"] == [] and set(res["status"]) == set(k22.NAMEABLE)


def test_the_ranking_the_file2_rule_and_at_most_two_names():
    all_up = {k: _table(xa=57.0, xu=58.0) for k in ALL_CELLS}
    res = k22.decide(_per(all_up), k22.ARMS, 0.5)
    assert res["ranked"] == ["x_unrel", "x_all"] and res["names"] == ["x_unrel", "x_all"]
    assert "file 1: x_unrel on base B" in res["outcome"] and "file 2: x_all on base B" in res["outcome"]
    assert "xunrel-<B>" in res["files"][0] and "xall-<B>" in res["files"][1]
    # file 2 only when corr(ln T_all, ln T_unrel) is below 0.9 (or known at all)
    for c in (0.9, 0.95, None):
        res = k22.decide(_per(all_up), k22.ARMS, c)
        assert res["names"] == ["x_unrel"] and res["file2_rule"]["applies"] and not res["file2_rule"]["named"]
    tie = {k: _table(xa=57.0, xu=57.0) for k in ALL_CELLS}
    assert k22.decide(_per(tie), k22.ARMS, 0.1)["names"] == ["x_all", "x_unrel"]
    arms = ("barred_form", "x_unrel", "x_unrel_perm")
    res = k22.decide(_per(tie, arms=arms), arms, None)
    assert res["status"]["x_all"] == {"ran": False, "qualified": False} and res["names"] == ["x_unrel"]


def test_the_doubt_statement_is_report_only():
    behind = {k: _table(rep=55.0) for k in ALL_CELLS}
    d = k22.doubt_statement(_per(behind), k22.ARMS)
    assert d["statement"] == "the planner's doubt holds on this fold" and len(d["behind_cells"]) == 4
    three = dict(behind)
    three[("dm", "g2")] = _table(rep=56.0)
    assert k22.doubt_statement(_per(three), k22.ARMS)["statement"] == "the doubt is not shown on this fold"
    assert not k22.doubt_statement(_per(three), ("barred_form", "x_unrel", "x_unrel_perm"))["ran"]
    res = k22.decide(_per({k: _table(rep=60.0) for k in ALL_CELLS}), k22.ARMS, 0.1)
    assert res["names"] == []


def test_on_pure_noise_a_name_is_rare():
    rng = np.random.default_rng(12)
    named = 0
    for _ in range(300):
        per = {}
        for w in k22.WORLDS:
            per[w] = {}
            common = {a: rng.normal(0, 1.0, len(SEEDS)) for a in k22.ARMS}
            for n in k22.TASKS:
                noise = {a: common[a] if w == "plain" else rng.normal(0, 1.0, len(SEEDS)) for a in k22.ARMS
                         if a != k22.BASE}
                per[w][n] = k22._summarise(_rows({a: BASE_SCORE for a in k22.ARMS}, noise=noise), k22.ARMS)
        named += bool(k22.decide(per, k22.ARMS, 0.1)["names"])
    assert named <= 6


# ------------------------------------------------------------------------------------------------ the reading


def _prov(commit="abc", clean=True, machine="m1"):
    return {"commit": commit, "tree_clean": clean, "machine": machine, "utc": "x"}


def _pre(dropped=None) -> dict:
    return {"n_genes": 500, "phi_m": {}, "tables": {}, "products": {}, "controls": {}, "dropped": dropped or {},
            "corr_ln_all_ln_unrelated": 0.3}


def _write(root: Path, cells=None, commit="abc", clean=True, machine="m1", shift=0.0, drop=None, dropped=None,
           ext_rec=None, bad_ext_file=None):
    cells = cells or {}
    pre = _pre(dropped)
    arms = k22.active_arms(pre)
    erec = ext_rec or {**PINS, "combine_sha256": "c" * 64}
    root.mkdir(parents=True)
    rec = {"dropped": k22.dropped_arms(pre), "arms": list(arms), "ext": erec}
    k18.write_task_file(root / "pregate.json", {"pregate": pre, "arms": list(arms), "ext": erec,
                                                "provenance": _prov(commit, clean, machine)})
    refs = {}
    for w in k22.WORLDS:
        for n in k22.TASKS:
            a = (harness.FLOOR, *arms) if w == "plain" else arms
            rows = _rows({x: cells.get((w, n), _table())[x] for x in a})
            refs[(w, n)] = [{"seed": r["seed"], "scores": {k22.BASE: r["scores"][k22.BASE] - shift}} for r in rows]
            for r in rows:
                s = r["seed"]
                f = root / (f"plain_{n}/seed_{s}.json" if w == "plain" else f"dm/{n}_seed_{s}.json")
                f.parent.mkdir(parents=True, exist_ok=True)
                if drop == str(f.relative_to(root)):
                    continue
                bad = bad_ext_file == str(f.relative_to(root))
                rr = {**rec, "ext": {**erec, "combine_sha256": "d"}} if bad else rec
                body = ({"task": {"key": k22.TASKS[n].key}, "rows": [r], "seed": s, "pregate": rr,
                         "provenance": _prov(commit, clean, machine)} if w == "plain" else
                        {**r, "check": {k18.DM_CHECK: 60.0 + s}, "target": k18.DM_TARGETS[n], "task": n,
                         "pregate": rr, "arms": list(arms), "provenance": _prov(commit, clean, machine)})
                k18.write_task_file(f, body)
    k18_dir = root.parent / f"k18_{root.name}"
    k18_dir.mkdir(parents=True, exist_ok=True)
    for n in k22.TASKS:
        (k18_dir / f"plain_{n}.json").write_text(json.dumps({"rows": refs[("plain", n)]}))
    (k18_dir / "dm.json").write_text(json.dumps({"rows": {n: refs[("dm", n)] for n in k22.TASKS}}))
    pins = {f: hashlib.sha256((k18_dir / f).read_bytes()).hexdigest() for f in k22.K18_SHA256}
    return root, k18_dir, pins


@pytest.fixture
def u3d(tmp_path):
    rows = [{"target": t, "seed": s, "world": "depth_matched", "variant": "dm:transfer:mab21l2_all", "t3": 60.0 + s}
            for t in ("rep1", "rep2") for s in SEEDS]
    p = tmp_path / "orth_records.csv"
    pd.DataFrame(rows).to_csv(p, index=False)
    return p


def test_the_reading_of_a_run_folder(tmp_path, u3d, monkeypatch):
    up = _table(xu=57.0)
    root, refs, pins = _write(tmp_path / "ok", {("plain", "g1"): up, ("dm", "g2"): up})
    monkeypatch.setattr(k22, "K18_SHA256", pins)
    res = k22.reading(root, refs, u3d)
    assert not res["void"] and res["names"] == ["x_unrel"] and res["tries_counted"] == 5
    assert all(d["max_abs_diff"] == 0.0 for d in res["guard"].values()) and len(res["files_sha256"]) == 41
    assert res["doubt"]["statement"] == "the doubt is not shown on this fold"
    txt = k22.text(res)
    assert "WORLD plain" in txt and "file 1: x_unrel on base B" in txt and "doubt statement" in txt
    assert "x_unrel - x_unrel_perm" in txt and "barred_form - floor" in txt and "pre-gate" in txt


@pytest.mark.parametrize("kw,why", [({"shift": 0.2}, "guard"), ({"clean": False}, "dirty tree"),
                                    ({"bad_ext_file": "dm/g1_seed_3.json"}, "tables differ"),
                                    ({"ext_rec": {"listing_sha256": "0" * 64, "tables": PINS["tables"],
                                                  "combine_sha256": "c"}}, "not the pinned ones")])
def test_void_readings(tmp_path, u3d, monkeypatch, kw, why):
    root, refs, pins = _write(tmp_path / "v", **kw)
    monkeypatch.setattr(k22, "K18_SHA256", pins)
    res = k22.reading(root, refs, u3d)
    assert res["void"] and any(why in w for w in res["why"]), res["why"]
    assert "no file" in res["outcome"] and "VOID" in k22.text(res)


def test_a_dropped_table_reads_with_fewer_tries(tmp_path, u3d, monkeypatch):
    root, refs, pins = _write(tmp_path / "d", dropped={"all": "test"})
    monkeypatch.setattr(k22, "K18_SHA256", pins)
    res = k22.reading(root, refs, u3d)
    assert not res["void"] and res["tries_counted"] == 2 and res["arms"] == ["barred_form", "x_unrel", "x_unrel_perm"]
    assert res["status"]["x_all"]["ran"] is False and not res["doubt"]["ran"]


def test_read_writes_once_records_its_code_and_never_overwrites(tmp_path, u3d, monkeypatch, capsys):
    root, refs, pins = _write(tmp_path / "r")
    monkeypatch.setattr(k22, "K18_SHA256", pins)
    monkeypatch.setattr(k22, "K18_TASKS", refs)
    monkeypatch.setattr(k18, "U3D_RECORDS", u3d)
    monkeypatch.setattr(k22, "reading", lambda runs, refs=refs, u3d=u3d, f=k22.reading: f(runs, refs, u3d))
    out = tmp_path / "reading"
    k22.main(["read", "--runs", str(root), "--out", str(out)])
    body = json.loads((out / "reading.json").read_text())
    assert body["reading_code"]["script_sha256"] == hashlib.sha256((T3 / "k22_outside_forms.py").read_bytes()
                                                                   ).hexdigest()
    assert body["reading_code"]["outside_forms_sha256"] == hashlib.sha256((T3 / "outside_forms.py").read_bytes()
                                                                          ).hexdigest()
    first = (out / "reading.json").read_bytes()
    with pytest.raises(SystemExit) as e:
        k22.main(["read", "--runs", str(root), "--out", str(out)])
    assert e.value.code == 4 and (out / "reading.json").read_bytes() == first
    capsys.readouterr()


def test_a_missing_file_is_not_a_reading(tmp_path, u3d, monkeypatch, capsys):
    root, refs, pins = _write(tmp_path / "i", drop="plain_g2/seed_7.json")
    monkeypatch.setattr(k22, "K18_SHA256", pins)
    (root / "dm" / "g1_seed_0.json").write_text("{ not JSON: never opened")
    res = k22.reading(root, refs, u3d)
    assert res["incomplete"] and res["missing"] == ["plain_g2/seed_7.json"]
    out = tmp_path / "reading"
    with pytest.raises(SystemExit) as e:
        k22.main(["read", "--runs", str(root), "--out", str(out)])
    assert e.value.code == 3 and not out.exists() and "not a reading" in capsys.readouterr().out


# ------------------------------------------------------------------------------------------------ the plain run


@pytest.fixture
def fake_run(monkeypatch, tmp_path):
    calls = {"seeds": []}
    lines = {}
    tries = tmp_path / "TRIES.md"
    tries.write_text("| K22 (one row) | T3 | `experiments/t3_decomp/k22_outside_forms.py` |\n")

    def run_task(t, cands, seeds, data_dir, log=print):
        calls["seeds"] += list(seeds)
        calls["names"] = [c.name for c in cands]
        for s in seeds:
            lines[s] = f"{t.key} seed {s}: floor 50.00, barred_form {60 - s / 9:.2f}, x_all {61 + s / 7:.2f}"
            log(lines[s])
        return {"task": {"key": t.key}, "rows": [{"seed": s} for s in seeds],
                "provenance": {"commit": "abc", "tree_clean": True}}

    phi = np.clip(np.exp(np.random.default_rng(1).normal(0, 0.6, 500)), 0.25, 4.0)

    def load(*a):
        k18._BARRED.clear()
        k18._BARRED.update(phi=phi, genes=[f"g{i}" for i in range(500)], spec={}, record={})
        return {}

    def load_ext(genes, phi_m, basis=None, den=None):
        pre = _pre()
        k22._EXT.clear()
        k22._EXT.update(genes=list(genes), T={}, record={**PINS, "combine_sha256": "c" * 64}, pre=pre)
        return pre

    monkeypatch.setattr(k22, "TRIES_MD", tries)
    monkeypatch.setattr(k22, "preconditions", lambda: ("K22", k18.stamp()))
    monkeypatch.setattr(k22, "load_ext", load_ext)
    monkeypatch.setattr(k18, "load_barred_for_plain", load)
    monkeypatch.setattr(harness, "run_task", run_task)
    monkeypatch.setattr(k21.k7, "_PREP", {"fold": {"record": {}}})
    monkeypatch.setattr(k18, "stamp", lambda: {"commit": "abc", "tree_clean": True, "machine": "m1", "utc": "x"})
    return SimpleNamespace(calls=calls, lines=lines)


def test_the_pregate_writes_its_file_before_any_score_and_exits_5_when_both_tries_drop(tmp_path, monkeypatch,
                                                                                         capsys):
    import vec.prep

    genes = [f"g{i}" for i in range(G)]
    C = np.random.default_rng(0).poisson(3.0, size=(30, G)).astype(np.int32)
    C[C.sum(1) == 0, 0] = 1
    wt = ad.AnnData(X=np.log1p(C).astype(np.float32), var=pd.DataFrame(index=genes), layers={"counts": C})
    phi = np.clip(np.exp(np.random.default_rng(1).normal(0, 0.6, G)), 0.25, 4.0)
    seen: dict = {}

    def load(*a):
        k18._BARRED.clear()
        k18._BARRED.update(phi=phi, genes=genes, spec={}, record={"file": "f"})
        return {"file": "f"}

    def load_ext(g, phi_m, basis=None, den=None):
        seen.update(basis=basis, den=den)
        pre = _pre(seen.get("dropped"))
        k22._EXT.clear()
        k22._EXT.update(genes=list(g), T={}, record={**PINS, "combine_sha256": "c" * 64}, table_records={}, pre=pre)
        return pre

    monkeypatch.setattr(k22, "preconditions", lambda: ("K22", _prov()))
    monkeypatch.setattr(k18, "load_barred_for_plain", load)
    monkeypatch.setattr(vec.prep, "load_stage", lambda p, drop_duplicates=False: wt)
    monkeypatch.setattr(k22, "load_ext", load_ext)
    k22.main(["pregate", "--out", str(tmp_path / "a")])
    body = json.loads((tmp_path / "a" / "pregate.json").read_text())
    assert body["arms"] == list(k22.ARMS) and body["row"] == "K22" and body["ext"]["listing_sha256"] == "l" * 64
    assert np.allclose(seen["den"], (C / C.sum(1, keepdims=True)).mean(0)) and seen["basis"] is not None
    assert "pre-gate (the tables and phi_M alone; no score)" in capsys.readouterr().out
    with pytest.raises(SystemExit) as e:
        k22.main(["pregate", "--out", str(tmp_path / "a")])
    assert e.value.code == 2
    seen["dropped"] = {"all": "test", "unrelated": "test"}
    with pytest.raises(SystemExit) as e:
        k22.main(["pregate", "--out", str(tmp_path / "b")])
    assert e.value.code == 5 and (tmp_path / "b" / "pregate.json").exists()
    assert "nothing runs" in capsys.readouterr().out


def test_the_plain_run_is_seed_by_seed_resumable_and_logs_progress_only(tmp_path, fake_run, capsys):
    out = tmp_path / "run"
    k22.main(["run", "--world", "plain", "--task", "g1", "--seeds", "0-1", "--out", str(out)])
    assert fake_run.calls["seeds"] == [0, 1] and fake_run.calls["names"] == [harness.FLOOR, *k22.ARMS]
    printed = capsys.readouterr()
    for x in printed.out.splitlines():
        assert re.fullmatch(r"K22 plain g1 seed \d: scored|wrote \S+/seed_\d\.json: commit abc, clean True, "
                            r"machine m1", x), x
    for line in fake_run.lines.values():
        for score in re.findall(r"\d+\.\d\d", line):
            assert score not in printed.out and score not in printed.err
    body = json.loads((out / "plain_g1" / "seed_0.json").read_text())
    assert body["pregate"]["arms"] == list(k22.ARMS) and body["pregate"]["ext"]["listing_sha256"] == "l" * 64
    k22.main(["run", "--world", "plain", "--task", "g1", "--seeds", "0-2", "--out", str(out)])
    assert fake_run.calls["seeds"] == [0, 1, 2] and "seed 0: present, skipped" in capsys.readouterr().out
    assert k18.ARMS == k22._K18_ORIGINAL["ARMS"] and k18.build_arm is k22._K18_ORIGINAL["build_arm"]


def test_a_crash_is_masked_and_writes_nothing(tmp_path, fake_run, monkeypatch, capsys):
    def boom(*a, **k):
        raise RuntimeError("score 61.234 went wrong")

    monkeypatch.setattr(harness, "run_task", boom)
    with pytest.raises(SystemExit) as e:
        k22.main(["run", "--world", "plain", "--task", "g1", "--seeds", "0", "--out", str(tmp_path / "r")])
    err = capsys.readouterr().err
    assert e.value.code == 1 and "61.234" not in err and "<num>" in err
    assert not (tmp_path / "r" / "plain_g1" / "seed_0.json").exists()
    assert k18.ARMS == k22._K18_ORIGINAL["ARMS"]
