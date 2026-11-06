"""experiments/t1_hygiene/n12_hygiene_size.py (TRIES.md row N12, R2's gate) on synthetic inputs only: the factor as
build_t1 applies it, the top lists and their ties, the gate at its boundaries, the decile-matched null and its seed,
the reported readings, and a synthetic end-to-end run that writes stamped outputs and no cell count.

No real gene list is typed: synthetic symbols are made up; where the run needs the committed extracts it reads them."""

from __future__ import annotations

import copy
import hashlib
import importlib.util
import json
import re
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse
from scipy.stats import spearmanr

from vec.provenance import read_stamp, stamp_path

ROOT = Path(__file__).resolve().parents[1]
T1_CELLS = 2357  # the synthetic E9.5 file's cell count: never printed or written
# synthetic Feng cells per label, (E9.5 halves A, B), (E14.5 halves A, B), each above any gene count of the run; blood
# stays below min_cells_ext
FENG_CELLS = {"ventricular_cm": ((1311, 1307), (1419, 1401)), "atrial_cm": ((1113, 1109), (1127, 1131)),
              "endocardial_ec": ((1173, 1179), (1181, 1191)), "vascular_ec": ((1137, 1139), (1149, 1151)),
              "epicardial": ((1157, 1163), (1167, 1197)), "blood": ((3, 2), (4, 1))}


@pytest.fixture(scope="module")
def n12():
    spec = importlib.util.spec_from_file_location("t1_hygiene_n12_hygiene_size",
                                                  ROOT / "experiments" / "t1_hygiene" / "n12_hygiene_size.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def cfg(n12):
    return n12.hygiene.load_config()


def _cfg(cfg, draws=None, seed=None):
    c = copy.deepcopy(cfg)
    if draws is not None:
        c["n12"]["decile_null"]["draws"] = draws
    if seed is not None:
        c["n12"]["decile_null"]["seed"] = seed
    return c


def _two_groups(n=1600, n_pos=700):
    """Synthetic T1 genes and two groups: A ranks s0000 highest and s0700 lowest; B is A with the sign flipped. The
    last n - 2 n_pos genes have no change."""
    genes = [f"s{i:04d}" for i in range(n)]
    d = np.zeros(n)
    d[:n_pos] = np.linspace(2.0, 0.01, n_pos)
    d[n_pos:2 * n_pos] = -np.linspace(2.0, 0.01, n_pos)
    basis = {g: np.linspace(1.0, 50.0, n) for g in ("A", "B")}
    return genes, {"A": d, "B": -d}, basis


def _modules(d=(), a1=(), n1=(), n2=(), t5=()):
    return {"D": list(d), "A1": list(a1), "N1": list(n1), "N2": list(n2), "D_table5": list(t5)}


# ---------------------------------------------------------------------------------------------------------------------
# The factor, the lists and the gate


def test_log_factor_is_build_t1s_zero_where_unmapped_and_held(n12):
    d = np.array([1.0, np.nan, -2.0, 3.0, 0.5])
    held = np.array([0, 0, 0, 1, 0], dtype=bool)
    np.testing.assert_allclose(n12.m1_log_factor(d, 0.4, held), [0.4, 0.0, -0.8, 0.0, 0.2])
    assert np.isnan(d[1]) and d[3] == 3.0  # the input is left as it was


def test_top_lists_break_ties_by_symbol(n12):
    genes = [f"s{i:04d}" for i in range(1000)]
    genes[398:403] = ["t_e", "t_c", "t_a", "t_d", "t_b"]  # five genes tied at the 400th place, names scrambled
    logf = np.zeros(1000)
    logf[:398] = np.linspace(3.0, 1.0, 398)
    logf[398:403] = 0.5
    logf[403:700] = np.linspace(0.4, 0.01, 297)
    logf[700:1000] = -np.linspace(3.0, 0.01, 300)
    r = n12.rankings(logf, genes)
    top = [genes[i] for i in r["up"][:400]]
    assert top[-2:] == ["t_a", "t_b"] and "t_c" not in top
    assert [genes[i] for i in r["up"][398:403]] == ["t_a", "t_b", "t_c", "t_d", "t_e"]
    flat = np.array([-1.0, -1.0, 2.0, 2.0, 0.0])
    r2 = n12.rankings(flat, ["zz", "aa", "yb", "ya", "mid"])
    assert r2["up"].tolist() == [3, 2] and r2["down"].tolist() == [1, 0]  # 0 is on neither side


@pytest.mark.parametrize("counts, verdict", [
    ({"A": {"up": 133, "down": 0}, "B": {"up": 0, "down": 0}}, "build"),
    ({"A": {"up": 134, "down": 0}, "B": {"up": 0, "down": 0}}, "captain"),
    ({"A": {"up": 9, "down": 9}, "B": {"up": 9, "down": 9}}, "not built"),
    ({"A": {"up": 9, "down": 9}, "B": {"up": 9, "down": 10}}, "build"),
    ({"A": {"up": 400, "down": 0}, "B": {"up": 0, "down": 0}}, "captain"),
])
def test_gate_boundaries(n12, cfg, counts, verdict):
    assert n12.gate(counts, cfg["gate"]["upper"], cfg["gate"]["floor"]) == verdict


def test_gate_needs_a_list(n12):
    with pytest.raises(ValueError, match="no top list"):
        n12.gate({}, 133, 10)


@pytest.mark.parametrize("u_up, u_down, verdict", [(133, 0, "build"), (134, 0, "captain"), (9, 9, "not built"),
                                                   (10, 0, "build"), (9, 10, "build")])
def test_n12_decides_the_gate_from_the_unions_counts(n12, cfg, u_up, u_down, verdict):
    genes, dirs, basis = _two_groups()
    d_genes = genes[:u_up] + genes[700:700 + u_down]  # A's top up and top down; B's top down and top up
    res = n12.n12(genes, dirs, basis, [], _modules(d=d_genes), _cfg(cfg, draws=5))
    assert res["verdict"] == verdict
    assert res["gate"]["counts"] == {"A": {"up": u_up, "down": u_down}, "B": {"up": u_down, "down": u_up}}
    assert res["groups"]["A"]["lists"]["up_400"]["length"] == 400
    assert res["gate"]["lists_above_upper"] == (["A up", "B down"] if u_up > 133 else [])


def test_union_follows_the_lines(n12, cfg):
    genes, dirs, basis = _two_groups()
    mods = _modules(d=genes[:5], a1=genes[5:25], n1=genes[25:30], n2=genes[30:32])
    on = n12.n12(genes, dirs, basis, [], mods, _cfg(cfg, draws=5))
    assert on["gate"]["counts"]["A"]["up"] == 32 and on["settings"]["union_modules"] == ["D", "A1", "N1", "N2"]
    off = _cfg(cfg, draws=5)
    off["lines"]["A1"] = False
    res = n12.n12(genes, dirs, basis, [], mods, off)
    assert res["gate"]["counts"]["A"]["up"] == 12 and res["modules"]["A1"]["in_union"] is False
    assert res["groups"]["A"]["lists"]["up_400"]["counts"]["A1"] == 20  # still reported


def test_unmapped_and_held_genes_stay_at_zero_and_out_of_the_null(n12, cfg):
    genes, dirs, basis = _two_groups()
    dirs = {g: d.copy() for g, d in dirs.items()}
    for d in dirs.values():
        d[1] = np.nan  # s0001: no Feng value
    held = ["s0000", "Absent-syn"]  # s0000 has the largest change of A: held flat
    mods = _modules(d=["s0002"], n2=["s0000", "s0001", "Nowhere-syn"])
    res = n12.n12(genes, dirs, basis, held, mods, _cfg(cfg, draws=20))
    a = res["groups"]["A"]
    assert all(row["counts"]["N2"] == 0 for row in a["lists"].values())
    assert a["lists"]["up_400"]["module_genes"]["D"] == ["s0002"]
    assert a["null"]["N2"]["mapped"] == 2 and a["null"]["N2"]["in_universe"] == 0
    assert all(s["null_mean"] == 0 for s in a["null"]["N2"]["lists"].values())
    n2 = res["modules"]["N2"]
    assert (n2["source"], n2["mapped"], n2["unmapped"]) == (3, 2, ["Nowhere-syn"])
    assert n2["mapped_on_held_flat"] == ["s0000"] and n2["mapped_without_feng_value"] == ["s0001"]
    assert res["settings"]["held_flat_genes_in_t1"] == 1
    assert a["genes_with_feng_value"] == 1599 and a["nonzero_log_f"] == 1398  # 1400 changed, one NaN, one held
    assert res["unions"]["U"]["on_held_flat"] == ["s0000"]


# ---------------------------------------------------------------------------------------------------------------------
# The decile-matched null


def test_decile_bins_are_k2bs(n12):
    basis = np.array([5.0, 1.0, 1.0, 3.0, 9.0, 2.0, 7.0, 0.5, 4.0, 8.0, 6.0])
    universe = np.ones(11, dtype=bool)
    universe[4] = False
    bins = n12.decile_bins(basis, universe, 4)
    assert [len(b) for b in bins] == [3, 3, 2, 2]
    assert bins[0].tolist() == [7, 1, 2]  # lowest first; the tie at 1.0 keeps gene order
    assert 4 not in np.concatenate(bins)


def test_null_draws_keep_the_modules_decile_composition(n12):
    rng = np.random.default_rng(0)
    bins = [np.arange(0, 10), np.arange(10, 20), np.arange(20, 30)]
    members = np.array([1, 2, 3, 12, 35])  # three in bin 0, one in bin 1, one outside every bin
    k, drawn = n12.null_counts(members, bins, {"bin0": bins[0], "bin1": bins[1], "none": np.array([30, 31])}, 40,
                               200, rng)
    assert k == [3, 1, 0]
    assert (drawn["bin0"] == 3).all() and (drawn["bin1"] == 1).all() and (drawn["none"] == 0).all()


def test_null_is_fixed_by_its_seed(n12, cfg):
    genes, dirs, basis = _two_groups()
    rng = np.random.default_rng(5)
    mods = _modules(d=list(rng.choice(genes[:1400], 30, replace=False)), n2=list(rng.choice(genes, 19, replace=False)),
                    t5=list(rng.choice(genes, 140, replace=False)))
    c = _cfg(cfg, draws=200)
    one, two = n12.n12(genes, dirs, basis, [], mods, c), n12.n12(genes, dirs, basis, [], mods, c)
    other = n12.n12(genes, dirs, basis, [], mods, _cfg(cfg, draws=200, seed=20261107))
    for g in ("A", "B"):
        assert one["groups"][g]["null"] == two["groups"][g]["null"]
        assert one["groups"][g]["null"] != other["groups"][g]["null"]
        assert one["groups"][g]["lists"] == other["groups"][g]["lists"]  # only the null moves with the seed
    logf = n12.m1_log_factor(dirs["A"], 0.4, np.zeros(len(genes), bool))
    bins = n12.decile_bins(basis["A"], logf != 0, 10)
    own = [len(set(genes[i] for i in b) & set(mods["D"])) for b in bins]
    assert one["groups"]["A"]["null"]["D"]["per_decile"] == own
    assert one["settings"]["null"] == {"draws": 200, "seed": 20261106, "bins": 10}


def test_the_registered_null_runs_its_1000_draws(n12, cfg):
    genes, dirs, basis = _two_groups(n=900, n_pos=420)
    res = n12.n12(genes, dirs, basis, [], _modules(d=genes[::50]), cfg)
    s = res["groups"]["A"]["null"]["D"]["lists"]["up_400"]
    assert res["settings"]["null"] == {"draws": 1000, "seed": 20261106, "bins": 10}
    assert 0 <= s["null_q05"] <= s["null_mean"] <= s["null_q95"] and 1 / 1001 <= s["p_at_or_above"] <= 1


def test_null_summary_is_exact_on_given_draws(n12):
    draws = np.array([0, 1, 2, 3, 4, 5, 3, 3, 1, 0])
    s = n12.null_summary(3, draws)
    assert s["p_at_or_above"] == (1 + 5) / (1 + 10)  # five draws at or above 3
    assert s["null_mean"] == 2.2 and s["null_sd"] == pytest.approx(np.std(draws, ddof=1))
    assert (s["null_q05"], s["null_q95"]) == (np.percentile(draws, 5), np.percentile(draws, 95))
    assert n12.null_summary(6, draws)["p_at_or_above"] == 1 / 11 and n12.null_summary(0, draws)["p_at_or_above"] == 1


def test_null_is_exact_where_every_draw_is_certain(n12, cfg):
    """With _two_groups' basis rising with the gene index, A's universe (s0000-s1399) falls into bins of 140 genes in
    index order, and A's top 400 up are s0000-s0399: a module in bins 0 and 1 puts every drawn gene in that list, so
    every one of the 1,000 draws counts the module's size, and none enters the top 400 down."""
    genes, dirs, basis = _two_groups()
    d = genes[0:5] + genes[150:155]
    res = n12.n12(genes, dirs, basis, [], _modules(d=d), cfg)
    null = res["groups"]["A"]["null"]["D"]
    assert null["per_decile"] == [5, 5, 0, 0, 0, 0, 0, 0, 0, 0] and null["in_universe"] == 10
    assert null["lists"]["up_400"] == {"observed": 10, "null_mean": 10.0, "null_sd": 0.0, "null_q05": 10.0,
                                       "null_q95": 10.0, "p_at_or_above": 1.0}
    assert null["lists"]["down_400"] == {"observed": 0, "null_mean": 0.0, "null_sd": 0.0, "null_q05": 0.0,
                                         "null_q95": 0.0, "p_at_or_above": 1.0}
    assert res["groups"]["B"]["null"]["D"]["lists"]["down_400"]["null_mean"] == 10.0  # B is A with the sign flipped


@pytest.mark.parametrize("d_size, verdict", [(140, "captain"), (10, "build"), (9, "not built")])
def test_short_top_lists_keep_the_absolute_thresholds_and_report_their_length(n12, cfg, d_size, verdict):
    """A group with fewer than 400 genes on a side: its list is shorter, the gate still reads 133 and 10 as counts
    (not shares), and the report prints each list's length."""
    genes, dirs, basis = _two_groups(n=1000, n_pos=250)
    res = n12.n12(genes, dirs, basis, [], _modules(d=genes[:d_size]), _cfg(cfg, draws=5))
    a = res["groups"]["A"]["lists"]
    assert (a["up_100"]["length"], a["up_400"]["length"], a["up_1000"]["length"], a["down_400"]["length"]) \
        == (100, 250, 250, 250)
    assert res["gate"]["counts"] == {"A": {"up": d_size, "down": 0}, "B": {"up": 0, "down": d_size}}
    assert res["verdict"] == verdict and res["settings"]["upper"] == 133 and res["settings"]["floor"] == 10
    text = n12.report_text(res, {}, {"commit": "x" * 40, "tree_clean": True, "utc": "2026-11-06T00:00:00Z"}, {})
    assert f"  A: up {d_size} of 250, down 0 of 250" in text and f"VERDICT: {verdict} -- " in text
    assert "top 400 up (250 genes)" in text


def test_check_settings_refuses_a_pseudocount_that_is_not_feng_dirs(n12, cfg):
    n12.check_settings(cfg)  # the committed config passes
    bad = copy.deepcopy(cfg)
    bad["direction"]["pseudocount_cp10k"] = 1.0
    with pytest.raises(SystemExit, match=r"pseudocount 1\.0 is not feng_dirs' 0\.5"):
        n12.check_settings(bad)


# ---------------------------------------------------------------------------------------------------------------------
# The reported readings


def test_reported_readings(n12, cfg):
    genes, dirs, basis = _two_groups()
    genes[3], genes[710], genes[20] = "Rps9syn", "Cox5syn", "mt-Syn1"  # made-up symbols matching the patterns
    genes[21], genes[22] = "Cox9syn", "Rpl-syn"  # and two that do not
    mods = _modules(d=[genes[3], genes[10]], n1=["mt-Syn1"], t5=[genes[3], genes[10], genes[710]])
    res = n12.n12(genes, dirs, basis, [], mods, _cfg(cfg, draws=5))
    a = res["groups"]["A"]
    logf = 0.4 * dirs["A"]
    up100 = np.argsort(-logf, kind="stable")[:100]
    w = np.abs(logf[up100])
    assert a["lists"]["up_100"]["abs_log_f_share"]["D"] == pytest.approx((abs(logf[3]) + abs(logf[10])) / w.sum())
    assert a["lists"]["up_100"]["abs_log_f_share"]["U"] == pytest.approx(
        (abs(logf[3]) + abs(logf[10]) + abs(logf[20])) / w.sum())
    u = np.zeros(len(genes), bool)
    u[[3, 10, 20]] = True
    assert a["held_share"]["U"] == pytest.approx(np.abs(logf[u]).sum() / np.abs(logf).sum())
    assert a["spearman_with_without"]["U"] == pytest.approx(spearmanr(logf, np.where(u, 0, logf)).statistic)
    up, down = a["lists"]["up_400"], a["lists"]["down_400"]
    assert up["patterns"] == {"ribosomal": {"genes": 1, "held_by_U": 1}, "oxphos": {"genes": 0, "held_by_U": 0},
                              "mt": {"genes": 1, "held_by_U": 1}}
    assert down["patterns"]["oxphos"] == {"genes": 1, "held_by_U": 0}  # Cox5syn is in Table 5, not in U
    assert down["counts"]["U_table5"] == 1 and down["counts"]["U"] == 0
    assert up["pattern_genes"]["ribosomal"] == ["Rps9syn"] and up["module_genes"]["N1"] == ["mt-Syn1"]
    assert "module_genes" not in a["lists"]["up_100"]  # by name in the registered top lists only
    t5 = res["gate"]["with_table5_in_place_of_D"]
    assert t5["counts"]["A"] == {"up": 3, "down": 1} and t5["decides"] is False
    assert res["unions"]["U_table5"]["genes"] == 4 and res["unions"]["U"]["per_module"] == {"D": 2, "N1": 1}


# ---------------------------------------------------------------------------------------------------------------------
# A synthetic end-to-end run


def _feng_stats(root: Path, genes: list[str], seed: int = 3) -> Path:
    rng = np.random.default_rng(seed)
    keys, n, cp = [], [], []
    for label, ends in FENG_CELLS.items():
        for stage, halves in zip(("E9.5", "E14.5"), ends, strict=True):
            prof = rng.gamma(0.6, 20.0, len(genes))
            for half, k in zip("AB", halves, strict=True):
                keys.append(f"{stage}|LV|{label}|{half}")
                n.append(k)
                cp.append(k * prof)
    d = root / "data" / "external" / "feng2022_gse193346_geo"
    d.mkdir(parents=True)
    path = d / "cd1_stats.npz"
    n = np.array(n, dtype=np.float64)
    np.savez(path, genes=np.array(genes), keys=np.array(keys), n=n, S_cp=np.array(cp),
             S_log=np.zeros((len(keys), len(genes))), S_nz=np.zeros((len(keys), len(genes))), depth_sum=n * 1000.0)
    sha = hashlib.sha256(path.read_bytes()).hexdigest()
    (d / "MANIFEST.json").write_text(json.dumps({"outputs": {"cd1_stats.npz": {"sha256": sha}}}))
    return path


@pytest.fixture
def run_inputs(n12, tmp_path, monkeypatch, own_repo):
    hyg = n12.hygiene
    d, n2, t5 = hyg.read_extract("D"), hyg.read_extract("N2"), hyg.read_extract("D_table5")
    t5_only = [s for s in t5 if s not in set(d)]
    syn = [f"s{i:04d}" for i in range(800)]
    patterned = ["Hba-syn1", "Hbb-syn2", "Hbb-syn3", "mt-Syn1", "mt-Syn2"]
    t1 = syn + d[:12] + n2[:6] + t5_only[:5] + patterned
    feng = syn + [f"x{i}" for i in range(40)] + d[:10] + n2[:6] + t5_only[:5] + patterned + ["Hbb-syn3"]
    data = tmp_path / "data"
    (data / "T1").mkdir(parents=True)
    X = sparse.random(T1_CELLS, len(t1), density=0.01, random_state=0, format="csr", dtype=np.float32)
    a = ad.AnnData(X=X)
    a.var_names, a.obs_names = t1, [f"c{i}" for i in range(T1_CELLS)]
    a.write_h5ad(data / "T1" / "E9.5_RNA.h5ad")
    feng_root = tmp_path / "feng"
    _feng_stats(feng_root, feng)
    monkeypatch.setattr(n12.feng_dirs, "REPO", feng_root)
    held = tmp_path / "held_flat_genes.csv"
    pd.DataFrame({"gene": ["s0003", d[0], "Absent-syn"]}).to_csv(held, index=False)
    monkeypatch.setattr(n12.t1_scope, "HELD_FLAT_GENES", held)
    unpinned = dict(hyg.load_config(), union_sha256=None)  # as before N12's run, also once the build commit pins it
    monkeypatch.setattr(hyg, "load_config", lambda *a, **k: copy.deepcopy(unpinned))
    return {"data": data, "out": tmp_path / "out", "t1": t1, "d": d, "n2": n2, "held": held}


def test_synthetic_run_writes_stamped_outputs_and_no_cell_count(n12, run_inputs, capsys):
    out, inp = run_inputs["out"], run_inputs
    n12.main(["--out", str(out), "--data", str(inp["data"])])
    printed = capsys.readouterr().out
    files = {k: out / v for k, v in n12.NAMES.items()}
    text, body, csv = files["text"].read_text(), json.loads(files["json"].read_text()), files["union"].read_text()
    assert printed == text and re.search(r"^VERDICT: (build|captain|not built) -- ", text, re.M)
    assert body["verdict"] in ("build", "captain", "not built") and f"VERDICT: {body['verdict']} " in text
    assert sorted(body["groups"]) == ["CM", "Endocardium", "Endothelium", "Epicardium"]  # blood has too few cells
    # the union: the modules mapped by exact symbol; A1 keeps Feng's repeated symbol, once, though it has no Feng value
    t1 = set(inp["t1"])
    want = ({s for s in inp["d"] if s in t1} | {s for s in inp["n2"] if s in t1}
            | {"Hba-syn1", "Hbb-syn2", "Hbb-syn3", "mt-Syn1", "mt-Syn2"})
    rows = [ln for ln in csv.splitlines() if not ln.startswith("#")]
    assert rows[0] == "gene,modules" and {r.split(",")[0] for r in rows[1:]} == want
    assert body["modules"]["A1"]["source"] == 3 and body["modules"]["D"]["mapped"] == 12
    assert body["modules"]["A1"]["mapped_without_feng_value"] == ["Hbb-syn3"]
    assert body["modules"]["D"]["mapped_without_feng_value"] == inp["d"][10:12]
    assert body["modules"]["D"]["mapped_on_held_flat"] == [inp["d"][0]]
    assert body["modules"]["D"]["genes"] == inp["d"][:12]  # each module's list and its sha256 (PLAN 4.1, item 4)
    def listed(xs):
        return hashlib.sha256("".join(f"{x}\n" for x in xs).encode()).hexdigest()

    assert body["modules"]["D"]["sha256"] == listed(inp["d"][:12])
    assert body["modules"]["D"]["source_sha256"] == listed(inp["d"])
    assert len(body["modules"]["D_table5"]["unmapped"]) == 140 - 12 - 5
    # provenance: inputs, outputs, commit, clean flag, time; each output stamped with its own sha256
    sha = {k: hashlib.sha256(p.read_bytes()).hexdigest() for k, p in files.items()}
    assert body["outputs"] == {n12.NAMES["union"]: sha["union"], n12.NAMES["text"]: sha["text"]}
    assert body["union_file"]["sha256"] == sha["union"] and f"hygiene_genes.csv sha256 {sha['union']}" in text
    assert body["inputs"]["T1 gene names"]["sha256"] == hashlib.sha256(
        (inp["data"] / "T1" / "E9.5_RNA.h5ad").read_bytes()).hexdigest()
    assert body["inputs"]["held-flat list"]["sha256"] == hashlib.sha256(inp["held"].read_bytes()).hexdigest()
    assert body["inputs"]["config"]["sha256"] == hashlib.sha256(
        (ROOT / "configs" / "t1_hygiene.yaml").read_bytes()).hexdigest()
    prov = body["provenance"]
    assert re.fullmatch(r"[0-9a-f]{40}", prov["commit"]) and prov["tree_clean"] is True
    assert re.fullmatch(r"\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z", prov["utc"])
    head = [ln for ln in csv.splitlines() if ln.startswith("# ")]
    assert any(ln.startswith("# input T1 gene names: T1/E9.5_RNA.h5ad sha256 ") for ln in head)
    assert head[-1] == f"# commit {prov['commit']}, tree clean True, utc {prov['utc']}"
    for k, p in files.items():
        s = read_stamp(stamp_path(p))
        assert s["sha256"] == sha[k] and s["clean"] and s["commit"] == prov["commit"]
    assert body["feng"]["groups"]["CM"] == {"feng_labels": ["ventricular_cm", "atrial_cm"]}  # no Feng cell count
    # no cell count anywhere: T1's, or Feng's per group
    feng_totals = {sum(h) for e in FENG_CELLS.values() for h in e}
    feng_totals |= {sum(FENG_CELLS[a][i]) + sum(FENG_CELLS[b][i]) for a, b in (("ventricular_cm", "atrial_cm"),
                                                                             ("endocardial_ec", "vascular_ec"))
                    for i in (0, 1)}
    everything = printed + text + csv + files["json"].read_text() + "".join(
        stamp_path(p).read_text() for p in files.values())
    for n in {T1_CELLS, *(t for t in feng_totals if t > 100)}:
        assert not re.search(rf"(?<![\w.,]){n}(?![\w.,])", everything), n
        assert not re.search(rf"(?<![\w.,]){n:,}(?![\w.,])", everything), n


def test_the_union_file_is_what_the_build_reads_once_pinned(n12, run_inputs):
    n12.main(["--out", str(run_inputs["out"]), "--data", str(run_inputs["data"])])
    csv = run_inputs["out"] / n12.NAMES["union"]
    cfg = copy.deepcopy(n12.hygiene.load_config())
    with pytest.raises(n12.hygiene.Refused, match="pins no union_sha256"):
        n12.hygiene.load_union(csv, cfg)
    cfg["union_sha256"] = hashlib.sha256(csv.read_bytes()).hexdigest()
    u = n12.hygiene.load_union(csv, cfg)
    body = json.loads((run_inputs["out"] / n12.NAMES["json"]).read_text())
    assert u["genes"] == sorted(body["union_members"]) and len(u["genes"]) == body["unions"]["U"]["genes"]
    assert u["modules"] == body["unions"]["U"]["per_module"]


def test_run_stops_when_its_factor_is_not_external_shifts(n12, run_inputs, monkeypatch):
    real = n12.m1_log_factor
    monkeypatch.setattr(n12, "m1_log_factor", lambda d, fraction, held: real(d, 2 * fraction, held))
    with pytest.raises(SystemExit, match="differs from external_shift's report"):
        n12.main(["--out", str(run_inputs["out"]), "--data", str(run_inputs["data"])])
    assert not run_inputs["out"].exists()


def test_run_refuses_a_config_whose_fraction_is_not_the_bracket(n12, run_inputs, monkeypatch):
    bad = copy.deepcopy(n12.hygiene.load_config())
    bad["direction"]["fraction"] = 0.5
    monkeypatch.setattr(n12.hygiene, "load_config", lambda *a, **k: bad)
    with pytest.raises(SystemExit, match="is not bracket_fraction"):
        n12.main(["--out", str(run_inputs["out"]), "--data", str(run_inputs["data"])])
    assert not run_inputs["out"].exists()


def test_t1_gene_names_reads_only_the_var_index(n12, tmp_path, monkeypatch):
    """Names-only: of the released file, only var's index element is read (no matrix, layer, obs or uns element)."""
    C = sparse.random(37, 5, density=0.3, random_state=1, format="csr", dtype=np.float32)
    a = ad.AnnData(X=C, layers={"counts": C.copy()}, uns={"note": "synthetic"},
                   obs=pd.DataFrame({"celltype": ["x"] * 37}, index=[f"c{i}" for i in range(37)]))
    a.var_names = ["Syn-a", "Syn-b", "Syn-c", "Syn-d", "Syn-e"]
    a.var["gene_ids"] = [f"id{i}" for i in range(5)]
    path = tmp_path / "E9.5_RNA.h5ad"
    a.write_h5ad(path)
    seen, real = [], n12.ad.io.read_elem
    monkeypatch.setattr(n12.ad.io, "read_elem", lambda elem, *x, **k: seen.append(elem.name) or real(elem, *x, **k))
    assert n12.t1_gene_names(path) == ["Syn-a", "Syn-b", "Syn-c", "Syn-d", "Syn-e"]
    assert seen == ["/var/_index"]
    with n12.h5py.File(path, "r+") as f:  # not anndata's dataframe encoding: refused, not read whole
        del f["var"].attrs["_index"]
    with pytest.raises(SystemExit, match="cannot be read alone"):
        n12.t1_gene_names(path)


def test_tries_rows_are_the_table_rows_ids(n12, tmp_path):
    md = tmp_path / "TRIES.md"
    md.write_text("# Tries\n| ID | Board |\n|---|---|\n| R2 (ruling, no run) | T1 |\n| K2b | T3 |\n| E7.25 | x |\n"
                  "| T1:test | x |\n| Heart | x |\nN13 outside a table\n")
    assert n12.tries_rows(md) == {"R2", "K2b"}
    assert {"R2", "N12"} <= n12.tries_rows()  # the committed ledger's rows


def test_n12_runs_once_unless_a_new_tries_row_registers_a_rerun(n12, run_inputs, monkeypatch, capsys):
    """R2: N12's result is not revisited after its numbers are seen (a change is a new row): a second run, or a run
    while configs/t1_hygiene.yaml pins union_sha256, is refused before anything is read, unless --rerun-row names a
    TRIES row other than N12; that row is then recorded in every output."""
    out, data = run_inputs["out"], run_inputs["data"]
    argv = ["--out", str(out), "--data", str(data)]
    n12.main(argv)
    first = {p.name: p.read_bytes() for p in out.iterdir()}
    assert len(first) == 6  # three outputs and their stamps
    names = n12.t1_gene_names
    monkeypatch.setattr(n12, "t1_gene_names", lambda *a, **k: pytest.fail("refused before anything is read"))
    with pytest.raises(SystemExit, match="outputs exist already"):
        n12.main(argv)
    assert {p.name: p.read_bytes() for p in out.iterdir()} == first  # nothing overwritten
    cfg = n12.hygiene.load_config()  # run_inputs' unpinned config
    monkeypatch.setattr(n12.hygiene, "load_config", lambda *a, **k: dict(cfg, union_sha256="ab" * 32))
    with pytest.raises(SystemExit, match="pins union_sha256"):
        n12.main(["--out", str(out.parent / "elsewhere"), "--data", str(data)])
    assert not (out.parent / "elsewhere").exists()
    for row, why in (("N12", "other than N12"), ("Z99", "not a row"), ("n13", "not a row")):
        with pytest.raises(SystemExit, match=why):
            n12.main([*argv, "--rerun-row", row])
    monkeypatch.setattr(n12, "t1_gene_names", names)
    monkeypatch.setattr(n12, "tries_rows", lambda *a, **k: {"N12", "N13"})  # N13 stands for a registered new row
    capsys.readouterr()
    n12.main([*argv, "--rerun-row", "N13"])  # also while the union is pinned
    body = json.loads((out / n12.NAMES["json"]).read_text())
    assert body["provenance"]["rerun_row"] == "N13"
    rerun = ", rerun under experiments/harness/TRIES.md row N13"
    assert (out / n12.NAMES["text"]).read_text().rstrip("\n").endswith(rerun)
    head = [ln for ln in (out / n12.NAMES["union"]).read_text().splitlines() if ln.startswith("# ")]
    assert head[-1].endswith(rerun)
