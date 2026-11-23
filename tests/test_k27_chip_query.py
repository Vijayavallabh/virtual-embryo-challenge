"""Tests for experiments/t3_decomp/k27_chip_query.py (try K27, T3 candidate CHIPQ-S: the practice query's binding table
on top of row K23's named form at two sizes, with a placebo, a decile permutation and a size twin per size, after a
guard part). Offline, synthetic and fast: the made-up ChIP database of tests/test_chipq.py (antigens P01..,
genes g000..),
K18's plain fold fixture with stubbed construct, composition and assemble (as tests/test_k24_wider_tables.py), hand-made
listings and pins, synthetic task files. No competition data and no outside data; the ledger and the register are not
read. K23's committed reading is read for `two`'s pins only."""

from __future__ import annotations

import contextlib
import hashlib
import importlib.util
import io
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
T3X = REPO / "experiments" / "t3_ext"
for p in (T3, T3X):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))
import k18_built_files as k18  # noqa: E402
import k21_factor_forms as k21  # noqa: E402
import k23_outside_exponents as k23  # noqa: E402
import k27_chip_query as k27  # noqa: E402
import outside_forms as ext  # noqa: E402

cq, cf = k27.cq, k27.cf


def _helpers():
    """tests/test_chipq.py's synthetic database (make_db, FakeNet, kept_files), loaded under its own name."""
    spec = importlib.util.spec_from_file_location("k27_chipq_helpers", REPO / "tests" / "test_chipq.py")
    m = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = m
    spec.loader.exec_module(m)
    return m


H = _helpers()
GENES = H.GENES
PLACEBO_NAMES = sorted(set(H.PLACEBO))
TWO_SHA = "9490f524e40e146480c775f85ae8829ec9c005607104ad8ace129c5c46ed6228"


def _sha(p) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


@pytest.fixture(autouse=True)
def clean_state():
    yield
    k27._K27.clear()
    k18._FOLD.clear()


@pytest.fixture(autouse=True)
def no_pause(monkeypatch):
    monkeypatch.setattr(cf, "_sleep", lambda s: None)


# the synthetic database's placebos with bound sets of distinct sizes (in tests/test_chipq.py two pairs of placebo
# stand-ins share a size and so finish to equal tables, which P2 refuses: see the test of that refusal)
DISTINCT = {**H.PLACEBO, "P10": (2, 0, 24), "P12": (2, 0, 22)}


@pytest.fixture
def db(tmp_path, monkeypatch):
    monkeypatch.setattr(H, "PLACEBO", DISTINCT)
    H.use_net(monkeypatch, H.FakeNet(H.make_db()))
    root = tmp_path / "chipatlas"
    cf.step_lists(root)
    cf.step_rows(root)
    cf.step_targets(root)
    monkeypatch.setattr(cf, "DATA_ROOT", root)
    return root


def _kept(tmp_path, kept=None):
    d = tmp_path / "kept"
    d.mkdir(exist_ok=True)
    yp, sp = H.kept_files(d, kept)
    return cq.load_kept(yp, sp, check_committed=False)


@pytest.fixture
def kept(tmp_path):
    return _kept(tmp_path)


# -------------------------------------------------------------------------------- (1) constants and pins

BYTE_IDENTICAL = {
    "experiments/t3_decomp/outside_forms.py": "4837687d47f134ff690177ecaaec7c712952c238d0e7c806c5c8ef2cecb3ceb5",
    "experiments/t3_decomp/k22_outside_forms.py": "6482585242ad6213eb550ceac54f3b9f911da6f8b0c2901216ea74ff24c7dc19",
    "experiments/t3_decomp/k23_outside_exponents.py":
        "98bdac5b9f2413bc6d11e470fdbd2202df5ab3dc22548c8940812ac1f455eb98",
    "experiments/t3_decomp/k18_built_files.py": "4355e6d21fec4c7f9b89eef8508cc06c0174d9a36f120e75028a98d2ecc6c6a8",
    "experiments/t3_decomp/k21_factor_forms.py": "97ce564ddea8b03ad017f0d029b27c7d453ce75615d287e54860b8680f23c8c6",
    "experiments/t3_ext/t3x_tables.py": "f5886fd220d5040031b0746c841309478d9f414eb606c856b5f35d7d80d76690",
}


def test_the_shared_modules_are_byte_identical_and_the_reading_records_chipqs_bytes():
    assert {f: _sha(REPO / f) for f in BYTE_IDENTICAL} == BYTE_IDENTICAL
    assert {k: _sha(T3 / k) for k in ext.FROZEN_SHA256} == ext.FROZEN_SHA256
    src = (T3 / "k27_chip_query.py").read_text()
    for k in ("chipq_sha256", "chipq_fetch_sha256", "k18_sha256", "k23_sha256"):
        assert f'"{k}"' in src


def test_arms_sizes_and_constants_are_the_rows():
    assert k27.ARMS == ("barred_form", "x_e2", "chip_q_s1", "chip_rand_s1", "chip_perm_s1", "chip_size_s1",
                        "chip_q_s4", "chip_rand_s4", "chip_perm_s4", "chip_size_s4")
    assert len(k27.TRY_ARMS) == 8 == k27.N_TRIES and k27.GUARDS == ("barred_form", "x_e2")
    assert k27.S_STAR == 0.04949423311987974 == cq.S_STAR
    assert k27.SIZES == {"s1": 0.04949423311987974, "s4": 4 * 0.04949423311987974}
    assert (k27.MIN_GAIN, k27.MIN_AHEAD, k27.MAX_EXTRA_GATED) == (0.3, 8, 1)
    assert (k27.GUARD_MEAN, k27.GUARD_SEED, k27.MAX_ABS_CORR, k27.REL_TOL) == (0.1, 0.5, 0.5, 1e-12)
    assert (k27.P2_MIN_GENES_QUERY, k27.P2_MIN_GENES_PLACEBO) == (10, 1)
    assert (cq.P1_MIN_EXPERIMENTS, cq.P1_MIN_BOUND, cq.POOL_MIN_SIZE, cq.PLACEBO_SEED_BASE) == (3, 15, 10, 20261124)
    assert (ext.PERM_SEED, ext.PERM_BINS) == (20261113, 10)
    assert k27.TABLE_NAMES == ("chipq_Gata4", *(f"chipq_rand_seed{s}" for s in range(10)))
    assert k27.ENV_TABLES == "VEC_T3_EXT_TABLES_K27" and k27.ENV_TWO == "VEC_T3_EXT_TABLES_TWO"
    assert k27.K24_BOUNDS == {"max_abs_ln_scaled": 0.5, "max_sd_ratio": 1.05, "max_den_shift": 0.05,
                              "max_at_clip": 45, "e_range": [0.125, 8.0]}


def test_the_pins_and_the_listing_are_empty_until_the_second_commit():
    assert k27.pins_state() == "empty" and k27.K27_LISTING_SHA256 is None
    assert set(k27.K27_TABLE_SHA256) == set(k27.TABLE_NAMES)
    assert k27.LISTING_FILE.read_bytes() == b""
    with pytest.raises(k27.Refused, match="not filled"):
        k27.read_listing()


def _listing(tmp_path, tables: dict) -> tuple[Path, str]:
    p = tmp_path / "TABLES_K27.sha256"
    p.write_text(k27.listing_text(tables))
    return p, _sha(p)


def _hex(i: int) -> str:
    return hashlib.sha256(str(i).encode()).hexdigest()


def test_the_listing_against_the_pins(tmp_path):
    tabs = {n: _hex(i) for i, n in enumerate(k27.TABLE_NAMES)}
    p, sha = _listing(tmp_path, tabs)
    listing, got = k27.read_listing(p, sha, tabs, check_tracked=False)
    assert listing == tabs and got == sha
    assert k27.pins_state(sha, tabs) == "filled"
    assert k27.pins_state(None, {**tabs, "chipq_Gata4": None}) == "mixed"
    assert k27.pins_state(sha, {n: None for n in tabs}) == "mixed"
    for bad in ({**tabs, "chipq_rand_seed3": tabs["chipq_Gata4"]},):  # a repeated sha256 among the pins
        with pytest.raises(k27.Refused, match="twice"):
            k27.read_listing(p, sha, bad, check_tracked=False)
    with pytest.raises(k27.Refused, match="twice"):  # and among the listing lines
        k27.listing_text({**tabs, "chipq_rand_seed3": tabs["chipq_Gata4"]})
    with pytest.raises(k27.Refused, match="pinned K27_LISTING_SHA256"):
        k27.read_listing(p, _hex(99), tabs, check_tracked=False)
    p.write_text(k27.listing_text(tabs).replace(tabs["chipq_rand_seed2"], _hex(77)))
    with pytest.raises(k27.Refused):
        k27.read_listing(p, _sha(p), tabs, check_tracked=False)
    p.unlink()
    with pytest.raises(k27.Refused, match="missing"):
        k27.read_listing(p, sha, tabs, check_tracked=False)
    pins = k27.pins_text(tabs)
    assert pins.splitlines()[0] == f'K27_LISTING_SHA256: str | None = "{sha}"'
    assert all(f'    "{n}": "{s}",' in pins for n, s in tabs.items())


def test_p2_refuses_filled_pins_and_a_filled_listing(tmp_path, monkeypatch):
    p = tmp_path / "TABLES_K27.sha256"
    p.write_bytes(b"")
    k27.listing_empty(p, check_tracked=False)
    p.write_text("x\n")
    with pytest.raises(k27.Refused, match="not empty"):
        k27.listing_empty(p, check_tracked=False)
    p.write_bytes(b"")
    monkeypatch.setattr(k27, "K27_LISTING_SHA256", _hex(1))
    with pytest.raises(k27.Refused, match="already filled"):
        k27.listing_empty(p, check_tracked=False)


# -------------------------------------------------------------------------------- (2) P2 on the database


def test_p2_passes_on_the_synthetic_database_and_records_the_barred_gene(db, kept):
    res = k27.p2_compute(GENES, db, kept)
    assert res["passes"] and res["why"] == []
    assert (res["query_kept"], res["query_bound"], res["query_K"], res["pool_size"]) == (3, 20, 20, 10)
    assert res["query_genes_changed"] == 20 and res["barred_gene_antigen"] == 0
    assert sorted(res["placebos"], key=int) == [str(s) for s in SEEDS]
    assert sorted(v["pool_index"] for v in res["placebos"].values()) == list(range(10))
    for n, e in res["exponents"].items():
        assert set(e) == {"s1", "s4"} and e["s4"] == pytest.approx(4 * e["s1"], rel=1e-12), n
    real = cf._lists
    monkeypatch = pytest.MonkeyPatch()
    try:
        monkeypatch.setattr(cf, "_lists", lambda r, m: ({**real(r, m)[0], "Mab21l2": "u"}, real(r, m)[1]))
        assert k27.p2_compute(GENES, db, kept)["barred_gene_antigen"] == 1
    finally:
        monkeypatch.undo()


def test_p2_fails_on_its_bars(db, tmp_path, monkeypatch):
    thin = _kept(tmp_path, {"Ctnnb1": H.KEPT["Ctnnb1"], "Gata4": H.KEPT["Gata4"][:2]})
    res = k27.p2_compute(GENES, db, thin)
    assert not res["passes"] and "P1's bars" in res["why"][0]
    kept = _kept(tmp_path)
    real_pool = cq.placebo_pool
    monkeypatch.setattr(cq, "placebo_pool", lambda *a, **k: real_pool(*a, **k)[:9])
    res = k27.p2_compute(GENES, db, kept)
    assert not res["passes"] and any("fewer than 10" in w for w in res["why"])
    monkeypatch.undo()
    real_build = cq.build

    def flat(a, genes, *rest):
        r = real_build(a, genes, *rest)
        return r | {"C": np.ones(len(genes))} if a == "P02" else r

    monkeypatch.setattr(cq, "build", flat)
    monkeypatch.setattr(cf, "DATA_ROOT", db)
    res = k27.p2_compute(GENES, db, kept)
    assert not res["passes"] and any("differs from 1 at no gene" in w for w in res["why"])


def test_p2_refuses_a_placebo_equal_to_the_query_and_an_sd_off_s(db, kept, monkeypatch):
    real_build = cq.build
    q = real_build("Gata4", GENES, db, kept)

    def same(a, genes, *rest):
        r = real_build(a, genes, *rest)
        return r | {"C": q["C"].copy()} if a == "P07" else r

    monkeypatch.setattr(cq, "build", same)
    with pytest.raises(k27.Refused, match="equals the query"):
        k27.p2_compute(GENES, db, kept)
    monkeypatch.setattr(cq, "build", real_build)
    real_e = cq.exponent
    monkeypatch.setattr(cq, "exponent", lambda C, S: None if real_e(C, S) is None else real_e(C, S) * (1 + 1e-9))
    with pytest.raises(k27.Refused, match="not S within"):
        k27.p2_compute(GENES, db, kept)


def test_p2_refuses_two_equal_placebo_tables(db, kept, tmp_path, monkeypatch):
    monkeypatch.setattr(H, "PLACEBO", dict(H.PLACEBO) | {"P10": (2, 0, 18), "P12": (2, 0, 20)})
    H.use_net(monkeypatch, H.FakeNet(H.make_db()))
    root = tmp_path / "chipatlas_same"
    cf.step_lists(root)
    cf.step_rows(root)
    cf.step_targets(root)
    with pytest.raises(k27.Refused, match="twice among the 11"):
        k27.p2_compute(GENES, root, kept)


def test_p2_writes_the_eleven_tables_and_refuses_a_byte_equal_placebo(db, kept, tmp_path, monkeypatch):
    d = tmp_path / "tables_k27"
    d.mkdir()
    got = k27.p2_write_tables(GENES, d, db, kept)
    assert sorted(got) == list(k27.TABLE_NAMES)
    assert all(_sha(d / "factors" / f"{n}.csv") == s for n, s in got.items())
    d2 = tmp_path / "tables_k27_b"
    d2.mkdir()
    real = cq.tables

    def copies(genes, out, *rest):
        w = real(genes, out, *rest)
        q = (Path(out) / "factors" / "chipq_Gata4.csv").read_bytes()
        (Path(out) / "factors" / "chipq_rand_seed4.csv").write_bytes(q)
        return [(f, hashlib.sha256(q).hexdigest() if f == "chipq_rand_seed4.csv" else s) for f, s in w]

    monkeypatch.setattr(cq, "tables", copies)
    with pytest.raises(k27.Refused, match="byte for byte"):
        k27.p2_write_tables(GENES, d2, db, kept)


def _p2_env(tmp_path, monkeypatch, db, kept):
    listing = tmp_path / "TABLES_K27.sha256"
    listing.write_bytes(b"")
    d = tmp_path / "tables_k27"
    d.mkdir()
    st = {"commit": "a" * 40, "tree_clean": True, "machine": "m1", "utc": "x"}
    monkeypatch.setattr(k27, "LISTING_FILE", listing)
    monkeypatch.setattr(k27.ext, "tracked", lambda p: True)
    monkeypatch.setattr(k27, "row_present", lambda: "K27")
    monkeypatch.setattr(k27.k18, "stamp", lambda: st)
    monkeypatch.setattr(cq, "load_kept", lambda *a, **k: kept)
    monkeypatch.setattr(cq, "barred_genes", lambda *a: list(GENES))
    monkeypatch.setenv("VEC_T3_EXT_TABLES_K27", str(d))
    return SimpleNamespace(listing=listing, d=d, st=st, runs=tmp_path / "runs" / "k27_aaaaaaa")


NUM = re.compile(r"(?<![0-9a-f])\d+(?:\.\d+)?(?![0-9a-f])")


def _clean_of_names_and_scores(text_: str) -> None:
    """No panel gene symbol, no placebo antigen name and no database score (the synthetic scores are 50 and up; the
    counts printed are below 50) in `text_`, sha256 hex removed first."""
    t = re.sub(r"[0-9a-f]{64}|sha256", "", text_, flags=re.IGNORECASE)
    assert not re.search(r"\bg\d{3}\b", t)
    assert not any(re.search(rf"\b{re.escape(a)}\b", t) for a in PLACEBO_NAMES)
    assert all(float(x) < 50 for x in NUM.findall(t)), [x for x in NUM.findall(t) if float(x) >= 50]


def test_the_p2_step_prints_counts_and_listing_lines_only_and_records_p2_json(db, kept, tmp_path, monkeypatch,
                                                                                capsys):
    E = _p2_env(tmp_path, monkeypatch, db, kept)
    k27.main(["p2", "--out", str(E.runs)])
    out, err = capsys.readouterr()
    _clean_of_names_and_scores(out + err)
    body = json.loads((E.runs / "p2.json").read_text())
    _clean_of_names_and_scores(json.dumps({k: v for k, v in body.items() if k not in ("exponents", "provenance",
                                                                                         "readings")}))
    assert body["passes"] and body["barred_gene_antigen"] == 0 and body["kept_sha256"] == kept.sha256
    lines = out.splitlines()
    assert lines[-12:-1] == [f"{body['tables'][n]}  {n}" for n in k27.TABLE_NAMES]
    assert lines[-1] == f"listing sha256 {body['listing_sha256']}"
    assert body["listing_sha256"] == hashlib.sha256(body["listing"].encode()).hexdigest()
    assert "the analysis list holds an antigen for the barred factor's gene: 0" in out
    # made once
    with pytest.raises(SystemExit) as ex:
        k27.main(["p2", "--out", str(E.runs)])
    assert ex.value.code == 2
    # the pins step prints the listing and the pin lines, from the tables on disk
    k27.main(["pins", "--runs", str(E.runs)])
    pins_out = capsys.readouterr().out
    assert body["listing"] in pins_out and k27.pins_text(body["tables"]) in pins_out
    _clean_of_names_and_scores(pins_out)
    t = E.d / "factors" / "chipq_rand_seed5.csv"
    t.write_bytes(t.read_bytes() + b"\n")
    with pytest.raises(SystemExit) as ex:
        k27.main(["pins", "--runs", str(E.runs)])
    assert ex.value.code == 2 and "not the one P2 wrote" in capsys.readouterr().err


def test_the_p2_step_exits_3_on_a_bar_and_2_on_filled_pins_or_no_kept_list(db, tmp_path, monkeypatch, capsys):
    thin = _kept(tmp_path, {"Ctnnb1": H.KEPT["Ctnnb1"], "Gata4": H.KEPT["Gata4"][:2]})
    E = _p2_env(tmp_path, monkeypatch, db, thin)
    with pytest.raises(SystemExit) as ex:
        k27.main(["p2", "--out", str(E.runs)])
    assert ex.value.code == 3
    out = capsys.readouterr().out
    assert "P2: fails" in out and not list((E.d).glob("factors/*.csv"))
    assert json.loads((E.runs / "p2.json").read_text())["passes"] is False
    _clean_of_names_and_scores(out)
    monkeypatch.setattr(k27, "K27_TABLE_SHA256", {n: _hex(i) for i, n in enumerate(k27.TABLE_NAMES)})
    with pytest.raises(SystemExit) as ex:
        k27.main(["p2", "--out", str(tmp_path / "other")])
    assert ex.value.code == 2 and "already filled" in capsys.readouterr().err
    monkeypatch.undo()
    (tmp_path / "x").mkdir()
    E = _p2_env(tmp_path / "x", monkeypatch, db, thin)

    def absent(*a, **k):
        raise cq.Refused("chipq_kept.yaml is absent (the main session commits it after the screen): refused")

    monkeypatch.setattr(cq, "load_kept", absent)
    with pytest.raises(SystemExit) as ex:
        k27.main(["p2", "--out", str(E.runs)])
    assert ex.value.code == 2 and "absent" in capsys.readouterr().err


# -------------------------------------------------------------------------------- (3) loading the tables


@pytest.fixture
def written(db, kept, tmp_path, monkeypatch):
    """The 11 tables written by P2 into a K27 folder, the listing and the pins made from them."""
    d = tmp_path / "tables_k27"
    d.mkdir()
    tabs = k27.p2_write_tables(GENES, d, db, kept)
    p, sha = _listing(tmp_path, tabs)
    monkeypatch.setenv("VEC_T3_EXT_TABLES_K27", str(d))
    return SimpleNamespace(d=d, tabs=tabs, listing=p, sha=sha, fetch=k27.fetch_manifest_sha256(db))


def _load(w, **kw):
    return k27.load_k27_tables(GENES, w.listing, w.sha, dict(w.tabs), w.d, w.fetch, check_tracked=False, **kw)


def test_the_tables_load_by_key_and_their_exponents_hit_s(written):
    loaded = _load(written)
    assert sorted(loaded["tables"]) == list(k27.TABLE_NAMES) and loaded["listing_sha256"] == written.sha
    e = k27.exponents_of(loaded["tables"])
    for n, (T, rec) in loaded["tables"].items():
        assert rec["sha256"] == written.tabs[n]
        for s, S in k27.SIZES.items():
            assert float(np.std(np.log(T ** e[n][s]))) == pytest.approx(S, rel=1e-12)
    assert sorted(r["pool_index"] for n, (_, r) in loaded["tables"].items() if r["seed"] is not None) == \
        list(range(10))


@pytest.mark.parametrize("case", ["bytes", "manifest", "factors_json", "order", "fetch", "pin", "missing"])
def test_a_changed_table_or_record_is_refused(written, case):
    f = written.d / "factors" / "chipq_rand_seed2.csv"
    kw = {}
    if case == "bytes":
        f.write_bytes(f.read_bytes().replace(b"1.0\n", b"1.00\n", 1))
    elif case == "manifest":
        m = written.d / "chipq_rand_seed2_manifest.json"
        body = json.loads(m.read_text())
        body["seed"] = 3
        m.write_text(json.dumps(body))
    elif case == "factors_json":
        m = written.d / "factors" / "factors.json"
        body = json.loads(m.read_text())
        body["files"]["chipq_rand_seed2"]["sha256"] = _hex(5)
        m.write_text(json.dumps(body))
    elif case == "order":
        with pytest.raises(k27.Refused, match="order"):
            k27.load_k27_tables(GENES[::-1], written.listing, written.sha, dict(written.tabs), written.d,
                                written.fetch, check_tracked=False)
        return
    elif case == "fetch":
        with pytest.raises(k27.Refused, match="fetch manifest"):
            k27.load_k27_tables(GENES, written.listing, written.sha, dict(written.tabs), written.d, _hex(3),
                                check_tracked=False)
        return
    elif case == "pin":
        with pytest.raises(k27.Refused):
            k27.load_k27_tables(GENES, written.listing, written.sha, {**written.tabs, "chipq_rand_seed2": _hex(8)},
                                written.d, written.fetch, check_tracked=False)
        return
    else:
        f.unlink()
    with pytest.raises(k27.Refused):
        _load(written, **kw)


def test_the_counted_steps_refuse_without_the_listing_with_empty_pins_or_a_changed_kept_list(tmp_path, monkeypatch,
                                                                                           capsys):
    st = {"commit": "b" * 40, "tree_clean": True, "machine": "m1", "utc": "x"}
    monkeypatch.setattr(k27, "row_present", lambda: "K27")
    monkeypatch.setattr(k27.k18, "stamp", lambda: st)
    # empty pins (the state at the registration commit)
    with pytest.raises(SystemExit) as ex:
        k27.preconditions()
    assert ex.value.code == 2 and "not filled" in capsys.readouterr().err
    tabs = {n: _hex(i) for i, n in enumerate(k27.TABLE_NAMES)}
    p, sha = _listing(tmp_path, tabs)
    monkeypatch.setattr(k27, "K27_LISTING_SHA256", sha)
    monkeypatch.setattr(k27, "K27_TABLE_SHA256", dict(tabs))
    monkeypatch.setattr(k27.ext, "tracked", lambda q: True)
    # no listing
    monkeypatch.setattr(k27, "LISTING_FILE", tmp_path / "absent.sha256")
    with pytest.raises(SystemExit) as ex:
        k27.preconditions()
    assert ex.value.code == 2 and "missing" in capsys.readouterr().err
    # a changed listing
    q = tmp_path / "changed.sha256"
    q.write_text(p.read_text() + "\n")
    monkeypatch.setattr(k27, "LISTING_FILE", q)
    with pytest.raises(SystemExit) as ex:
        k27.preconditions()
    assert ex.value.code == 2 and "not the pinned" in capsys.readouterr().err
    # the kept list changed after P2, or P2 made at this commit, or on another machine
    p2 = {"passes": True, "tables": dict(tabs), "listing_sha256": sha, "kept_sha256": "k" * 64,
          "provenance": {"commit": "a" * 40, "tree_clean": True, "machine": "m1"}}
    monkeypatch.setattr(k27, "is_ancestor", lambda a, b: True)
    k27.p2_matches(p2, st, "k" * 64)
    with pytest.raises(k27.Refused, match="kept list"):
        k27.p2_matches(p2, st, "j" * 64)
    with pytest.raises(k27.Refused, match="wrong commit"):
        k27.p2_matches({**p2, "provenance": {**p2["provenance"], "commit": st["commit"]}}, st, "k" * 64)
    with pytest.raises(k27.Refused, match="another machine"):
        k27.p2_matches({**p2, "provenance": {**p2["provenance"], "machine": "m2"}}, st, "k" * 64)
    with pytest.raises(k27.Refused, match="not the pinned"):
        k27.p2_matches({**p2, "tables": {**tabs, "chipq_Gata4": _hex(50)}}, st, "k" * 64)
    with pytest.raises(k27.Refused, match="did not pass"):
        k27.p2_matches({**p2, "passes": False}, st, "k" * 64)
    monkeypatch.setattr(k27, "is_ancestor", lambda a, b: False)
    with pytest.raises(k27.Refused, match="wrong commit"):
        k27.p2_matches(p2, st, "k" * 64)


# -------------------------------------------------------------------------------- (4) P3


G3 = 200


def _synth(seed=0, n=G3, corr_phi=False, corr_two=False, flat=False):
    rng = np.random.default_rng(seed)
    phi = np.clip(np.exp(rng.normal(0, 0.6, n)), 0.25, 4.0)
    two = np.where(rng.random(n) < 0.6, np.exp(rng.normal(0, 0.02, n)), 1.0)
    tabs = {}
    for i, name in enumerate(k27.TABLE_NAMES):
        w = np.zeros(n)
        pos = np.random.default_rng(seed + 100 + i).choice(n, 30, replace=False)
        if corr_phi and name == k27.QUERY_TABLE:
            pos = np.argsort(np.log(phi))[:30]   # the lowest ln phi_M get the weight: aligned with phi_M
        if corr_two and name == k27.QUERY_TABLE:
            pos = np.argsort(np.log(two))[:30]
        w[pos] = np.linspace(1.0, 0.5, 30)
        tabs[name] = np.ones(n) if (flat and name == k27.QUERY_TABLE) else cq.raw_table(w)
    if corr_phi:
        tabs[k27.QUERY_TABLE] = np.exp(0.9 * np.log(phi) + 0.01 * rng.normal(size=n))
    if corr_two:
        tabs[k27.QUERY_TABLE] = np.exp(40 * np.log(two) + 0.001 * rng.normal(size=n))
    e = {nm: {s: cq.exponent(T, S) for s, S in k27.SIZES.items()} for nm, T in tabs.items()} if not flat else None
    den = rng.random(n)
    den /= den.sum()
    bases = {"plain": {"genes": None, "basis": rng.random(n)}, "dm": {"genes": None, "basis": rng.random(n)}}
    return SimpleNamespace(phi=phi, two=two, tabs=tabs, e=e, den=den, bases=bases)


def test_the_pregate_figures_the_twin_from_the_whole_arm_and_each_worlds_permutation():
    s = _synth(1)
    pre = k27.pregate_k27(s.phi, s.two, s.tabs, s.e, s.bases, s.den,
                          pool_index={n: i for i, n in enumerate(k27.PLACEBO_TABLE.values())})
    assert pre["dropped"] is None
    x = s.phi * s.two ** 2
    for k, S in k27.SIZES.items():
        v = pre["arms"][k]
        Ce = s.tabs[k27.QUERY_TABLE] ** s.e[k27.QUERY_TABLE][k]
        P = x * Ce
        a = np.std(np.log(P)) / np.std(np.log(s.phi))
        assert v["a"] == pytest.approx(a, rel=1e-12)
        assert v["k"] == pytest.approx((s.den @ P) / (s.den @ s.phi ** a), rel=1e-12)
        assert v["table"]["sd_ln"] == pytest.approx(S, rel=1e-12)
        for w in ("plain", "dm"):
            perm = harness.t3_decile_permutation(Ce, s.bases[w]["basis"], 20261113, 10)
            want = ext.figures(s.phi, x * perm, s.den)
            assert v["controls"][w] == pytest.approx(want)
        assert pre["k24_bounds"][k]["applied"] is False
        assert set(pre["k24_bounds"][k]) >= {"max_abs_ln_scaled", "sd_ratio", "den_shift", "genes_at_clip", "e"}
    assert v["controls"]["plain"] != v["controls"]["dm"]
    assert sorted(pre["placebos"], key=int) == [str(i) for i in range(10)]
    assert pre["corr"]["ln_phi_m"] == pytest.approx(ext.corr_between(s.tabs[k27.QUERY_TABLE], s.phi))
    assert pre["corr"]["ln_two_sq"] == pytest.approx(ext.corr_between(s.tabs[k27.QUERY_TABLE], s.two ** 2))
    c = k27.constants_of(pre)
    assert c["twins"]["s4"] == {"a": pre["arms"]["s4"]["a"], "k": pre["arms"]["s4"]["k"]}


@pytest.mark.parametrize("case", ["phi", "two", "undefined"])
def test_the_pregate_drops_on_either_correlation(case):
    s = _synth(2, corr_phi=case == "phi", corr_two=case == "two")
    if case == "undefined":
        tabs = dict(s.tabs)
        tabs[k27.QUERY_TABLE] = np.ones(G3)
        e = {n: dict(v) for n, v in s.e.items()}
        pre = k27.pregate_k27(s.phi, s.two, tabs, e, s.bases, s.den)
        assert pre["corr"]["ln_phi_m"] is None or pre["corr"]["ln_two_sq"] is None
    else:
        pre = k27.pregate_k27(s.phi, s.two, s.tabs, s.e, s.bases, s.den)
        c = pre["corr"]["ln_phi_m" if case == "phi" else "ln_two_sq"]
        assert abs(c) >= 0.5
    assert k27.dropped(pre)
    assert "corr" in pre["dropped"]


def test_the_corr_bar_is_inclusive(monkeypatch):
    s = _synth(3)
    monkeypatch.setattr(k27.ext, "corr_between", lambda a, b: 0.5 if b is s.phi else 0.1)
    assert k27.dropped(k27.pregate_k27(s.phi, s.two, s.tabs, s.e, s.bases, s.den))
    monkeypatch.setattr(k27.ext, "corr_between", lambda a, b: 0.1 if b is s.phi else -0.5)
    assert k27.dropped(k27.pregate_k27(s.phi, s.two, s.tabs, s.e, s.bases, s.den))
    monkeypatch.setattr(k27.ext, "corr_between", lambda a, b: 0.4999)
    assert not k27.dropped(k27.pregate_k27(s.phi, s.two, s.tabs, s.e, s.bases, s.den))


def test_check_world_compares_the_runs_chip_perm_with_the_pregates():
    s = _synth(4)
    pre = k27.pregate_k27(s.phi, s.two, s.tabs, s.e, s.bases, s.den)
    k27._K27["consts"] = k27.constants_of(pre)
    basis = s.bases["plain"]["basis"]
    F = {"phi_barred": s.phi, "two": s.two,
         "scaled": {n: {k: T ** s.e[n][k] for k in k27.SIZE_KEYS} for n, T in s.tabs.items()}}
    F["q_perm"] = {k: ext.permuted(F["scaled"][k27.QUERY_TABLE][k], basis) for k in k27.SIZE_KEYS}
    k27.check_world(pre, "plain", F)
    with pytest.raises(k27._Mismatch):
        k27.check_world(pre, "dm", F)


# -------------------------------------------------------------------------------- (5) the arms in a fold

G, N_POOL, N_DRAW = 14, 90, 150
FG = [f"h{i}" for i in range(G)]
FILE_GENES = FG[::-1] + ["extra"]


@pytest.fixture
def fold(monkeypatch):
    """K18's plain fold fixture (as tests/test_k24_wider_tables.py), the barred factor, `two` and 11 tables in a file
    gene order that differs from the reference's."""
    bs = harness._t3_head()
    rng = np.random.default_rng(5)
    C = rng.poisson(3.0, size=(N_POOL, G)).astype(np.int32)
    C[C.sum(1) == 0, 0] = 1
    lab = rng.integers(0, 3, N_POOL)
    P = (C / C.sum(1, keepdims=True)).astype(np.float32)
    V = SimpleNamespace(p=np.array([0.5, 0.3, 0.2]))
    ref = bs.Reference(preset="gata4", rules={}, genes=FG, V=V, pool=bs.t3d.Pool(C, lab, P, 3),
                       coords=rng.normal(size=(N_POOL, 3)))
    pool = ad.AnnData(X=bs.t3d.lognorm(C).astype(np.float32), var=pd.DataFrame(index=FG))
    group = np.array(["A", "B", "Unknown"])[rng.integers(0, 3, N_POOL)]
    phi = np.exp(rng.normal(0, 0.5, G))
    per = {"A": np.exp(rng.normal(0, 0.5, G)), "B": np.exp(rng.normal(0, 0.5, G))}
    prep = {"ref": ref, "phi": phi, "per": per, "group": group, "record": {}}
    f_s2 = np.array([0.6, 0.25, 0.15])
    barred_file = np.clip(np.exp(rng.normal(0, 0.8, G + 1)), 0.25, 4.0)
    two = np.where(rng.random(G + 1) < 0.7, np.exp(rng.normal(0, 0.05, G + 1)), 1.0)
    tabs = {}
    for i, n in enumerate(k27.TABLE_NAMES):
        w = np.zeros(G + 1)
        w[np.random.default_rng(50 + i).choice(G + 1, 5, replace=False)] = np.linspace(1.0, 0.5, 5)
        tabs[n] = cq.raw_table(w)
    e = {n: {s: cq.exponent(T, S) for s, S in k27.SIZES.items()} for n, T in tabs.items()}
    consts = {"e": e, "twins": {"s1": {"a": 1.02, "k": 0.97}, "s4": {"a": 1.08, "k": 0.93}}}
    spec, _, _ = k18.barred_spec()
    log = {"construct": [], "assemble": []}

    def construct(r, n, g, f=None, phi=None, w=1.0, ratio=None, alphas=None):
        assert r is ref and ratio is None and alphas is None
        rows = g.integers(0, N_POOL, n)
        log["construct"].append({"phi": None if phi is None else np.asarray(phi).copy(), "w": w, "row": rows.copy()})
        Cr = C[rows]
        if phi is not None:
            Cr = bs.t3d.apply_gene_factor(Cr, np.asarray(phi), P[rows], g)
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
    k27._K27.clear()
    k27._K27.update(genes=list(FILE_GENES), T=dict(tabs), two=two, consts=consts)

    def ctx(task: str, seed: int):
        return SimpleNamespace(task=SimpleNamespace(key=k18.TASKS[task].key, reference_role="wt"),
                               inputs={"wt": pool}, board=SimpleNamespace(genes=FG, needs_coords=False), seed=seed)

    order = lambda v: k18.in_gene_order(v, FILE_GENES, FG)  # noqa: E731
    yield SimpleNamespace(bs=bs, ctx=ctx, log=log, phi_barred=order(barred_file), two=order(two),
                          T={k: order(v) for k, v in tabs.items()}, basis=harness.t3_mean_expression(pool), e=e,
                          consts=consts)
    k18._FOLD.clear()


def test_every_arm_is_the_rows_product_and_draws_the_barred_s1s_cells(fold):
    pb, two = fold.phi_barred, fold.two
    x = pb * two ** 2
    for seed in (2, 7):
        want = {"barred_form": pb, "x_e2": x}
        for s in k27.SIZE_KEYS:
            Ce = fold.T["chipq_Gata4"] ** fold.e["chipq_Gata4"][s]
            pl = k27.PLACEBO_TABLE[seed]
            tw = fold.consts["twins"][s]
            want |= {f"chip_q_{s}": x * Ce, f"chip_rand_{s}": x * fold.T[pl] ** fold.e[pl][s],
                     f"chip_perm_{s}": x * harness.t3_decile_permutation(Ce, fold.basis, 20261113, 10),
                     f"chip_size_{s}": tw["k"] * pb ** tw["a"]}
        X = {}
        with k27.k18_arms(k27.ARMS):
            for a in k27.ARMS:
                fold.log["construct"].clear()
                fold.log["assemble"].clear()
                X[a] = k18._build_plain(fold.ctx("g1", seed), a)
                (cc,) = fold.log["construct"]
                assert cc["w"] == 1.0 and np.allclose(cc["phi"], want[a], rtol=1e-14, atol=0), a
                fold.log[a] = (cc["row"], fold.log["assemble"][0])
        base_rows, base_asm = fold.log["barred_form"]
        for a in k27.ARMS:
            rows, asm = fold.log[a]
            assert np.array_equal(rows, base_rows) and asm["draw"] == base_asm["draw"]
    # a different placebo per seed: the ten seeds' chip_rand factors differ
    F = k18._FOLD.get("F")
    assert F is None  # restored on exit
    with k27.k18_arms(k27.ARMS):
        k18._build_plain(fold.ctx("g1", 0), "barred_form")
        F = k18._FOLD["F"]
        rands = [k27.arm_factor("chip_rand_s1", F, s) for s in SEEDS]
    assert all(not np.array_equal(rands[i], rands[j]) for i in range(10) for j in range(i + 1, 10))
    with pytest.raises(RuntimeError, match="seed"):
        k27.arm_factor("chip_rand_s1", F, None)


def test_barred_form_and_x_e2_are_k18s_and_k23s_value_for_value(fold):
    k18._FOLD.clear()
    want = k18._build_plain(fold.ctx("g2", 3), k18.BARRED)
    import k22_outside_forms as k22

    k22._EXT.clear()
    k22._EXT.update(genes=list(FILE_GENES), T={"all": k27._K27["two"]}, record={})
    k23._K23.update({"2": {"a": 1.0, "k": 1.0}, "4": {"a": 1.0, "k": 1.0}})
    try:
        with k23.k18_arms(("barred_form", "x_e2")):
            x23 = k18._build_plain(fold.ctx("g2", 3), "x_e2")
    finally:
        k22._EXT.clear()
        k23._K23.clear()
    with k27.k18_arms(k27.ARMS):
        got = k18._build_plain(fold.ctx("g2", 3), "barred_form")
        x27 = k18._build_plain(fold.ctx("g2", 3), "x_e2")
    assert np.array_equal(got, want) and np.array_equal(x27, x23)
    assert k18.ARMS == k27._K18_ORIGINAL["ARMS"] and k18.fold_factors is k27._K18_ORIGINAL["fold_factors"]
    with pytest.raises(ZeroDivisionError):
        with k27.k18_arms(k27.GUARDS):
            raise ZeroDivisionError
    assert k18.build_arm is k27._K18_ORIGINAL["build_arm"]


# -------------------------------------------------------------------------------- (6) the rule

BASE_SCORE = 56.0
CELLS = [(w, n) for w in k27.WORLDS for n in k27.TASKS]


def _rows(scores: dict, gated: dict | None = None, noise: dict | None = None) -> list[dict]:
    gated = gated or {}
    out = []
    for s in SEEDS:
        sc = {a: v + 0.1 * s + (noise[a][s] if noise and a in noise else 0.0) for a, v in scores.items()}
        out.append({"seed": s, "scores": sc,
                    "skills": {a: {"de_score": v / 200, "de_direction": 0.3 + v / 400, "severity_slope": 0.6}
                               for a, v in sc.items()},
                    "t3": {a: {"severity": {"state": "r2_gate" if s < gated.get(a, 0) else "scored", "beta": 0.4,
                                            "r2": 0.05}} for a in sc}})
    return out


def _table(**arms) -> dict:
    return {harness.FLOOR: 50.0, **{a: BASE_SCORE for a in k27.ARMS}, "barred_form": 55.5, **arms}


def _per(cells: dict, gated: dict | None = None) -> dict:
    out: dict = {}
    for w in k27.WORLDS:
        a = (harness.FLOOR, *k27.ARMS) if w == "plain" else k27.ARMS
        out[w] = {n: k27.summarise(_rows({x: cells.get((w, n), _table())[x] for x in a}, (gated or {}).get((w, n))),
                                   a) for n in k27.TASKS}
    return out


def _everywhere(**arms) -> dict:
    return {c: _table(**arms) for c in CELLS}


def test_a_pass_at_each_size_alone():
    r = k27.decide(_per(_everywhere(chip_q_s1=57.0)))
    assert r["passes"] and r["named_size"] == "S*" and r["qualified"] == ["S*"]
    assert r["outcome"] == k27.S_PASS.format(size="S*") and "PASS at size S*:" in r["outcome"]
    r = k27.decide(_per(_everywhere(chip_q_s4=57.0)))
    assert r["passes"] and r["named_size"] == "4 S*" and r["qualified"] == ["4 S*"]
    assert "at size 4 S* on base B" in r["outcome"]


def test_both_sizes_with_and_without_the_pivot():
    r = k27.decide(_per(_everywhere(chip_q_s1=57.0, chip_q_s4=57.5)))
    assert r["qualified"] == ["S*", "4 S*"] and r["pivot"]["holds"] and r["named_size"] == "4 S*"
    r = k27.decide(_per(_everywhere(chip_q_s1=57.0, chip_q_s4=57.2)))  # 4 S* - S* = 0.2: not AHEAD
    assert r["qualified"] == ["S*", "4 S*"] and not r["pivot"]["holds"] and r["named_size"] == "S*"
    cells = _everywhere(chip_q_s1=57.0, chip_q_s4=57.5)
    cells[("dm", "g2")] = _table(chip_q_s1=57.4, chip_q_s4=57.05)   # 4 S* - S* = -0.35 in one cell: no pivot
    r = k27.decide(_per(cells))
    assert r["qualified"] == ["S*", "4 S*"] and not r["pivot"]["holds"] and r["named_size"] == "S*"


def test_clause_g_fails_a_size_whose_control_is_ahead_of_the_named_form():
    for ctrl in ("chip_rand", "chip_perm", "chip_size"):
        # chip_q beats its control by 0.6 and the control is itself 0.4 ahead of x_e2
        r = k27.decide(_per(_everywhere(chip_q_s4=57.0, **{f"{ctrl}_s4": 56.4})))
        st = r["status"]["s4"]
        assert st["a"] and st["d"] and st["d2"] and st["e"] and not st["g"] and not st["qualified"], ctrl
        assert not r["passes"] and r["outcome"] == k27.S_FAIL
        # the same control ahead at the other size does not touch this size
        r = k27.decide(_per(_everywhere(chip_q_s4=57.0, **{f"{ctrl}_s1": 56.4})))
        assert r["named_size"] == "4 S*"


@pytest.mark.parametrize("letter,cells", [
    ("a", {c: _table(chip_q_s1=57.0 if c[0] == "plain" else 56.0) for c in CELLS}),   # no dm cell AHEAD
    ("b", {**_everywhere(chip_q_s1=57.0), ("dm", "g2"): _table(chip_q_s1=55.6)}),
    ("d", _everywhere(chip_q_s1=57.0, chip_rand_s1=56.9)),
    ("d2", _everywhere(chip_q_s1=57.0, chip_perm_s1=56.9)),
    ("e", _everywhere(chip_q_s1=57.0, chip_size_s1=56.9)),
])
def test_each_clause_failing(letter, cells):
    r = k27.decide(_per(cells))
    st = r["status"]["s1"]
    assert not st[letter] and not st["qualified"] and not r["passes"]


def test_clause_c_gated_seeds_and_the_gate_subclause_of_a():
    gated = {c: {"chip_q_s1": 2} for c in CELLS}   # x_e2 has 0: two more gated seeds
    r = k27.decide(_per(_everywhere(chip_q_s1=57.0), gated))
    assert not r["status"]["s1"]["c"] and not r["status"]["s1"]["a"]
    gated = {c: {"chip_q_s1": 1} for c in CELLS if c[0] == "plain"}
    r = k27.decide(_per(_everywhere(chip_q_s1=57.0), gated))
    assert r["status"]["s1"]["c"] and r["status"]["s1"]["qualified"]


def test_clause_f_needs_a_positive_de_change():
    per = _per(_everywhere(chip_q_s1=57.0))
    for w in k27.WORLDS:
        for n in k27.TASKS:
            per[w][n]["de_change"]["chip_q_s1"] = -0.001
    r = k27.decide(per)
    assert not r["status"]["s1"]["f"] and not r["passes"]


def test_nothing_qualifies_on_the_named_form_alone():
    r = k27.decide(_per(_everywhere()))
    assert not r["passes"] and r["outcome"] == k27.S_FAIL and r["named_size"] is None


def test_the_summary_pairs_and_the_de_change():
    v = _per(_everywhere(chip_q_s1=57.0))["plain"]["g1"]
    for s in k27.SIZE_KEYS:
        for c in k27.CONTROLS:
            assert f"chip_q_{s} - {c}_{s}" in v["paired"]
    assert "chip_q_s4 - chip_q_s1" in v["paired"] and "chip_q_s1 - x_e2" in v["paired"]
    assert v["paired"]["chip_q_s1 - x_e2"]["mean"] == pytest.approx(1.0)
    assert v["de_change"]["chip_q_s1"] == pytest.approx(1.0 / 200 + 1.0 / 400)


# -------------------------------------------------------------------------------- (7) the guard part


def _guard_refs():
    refs = {"k18": {}, "k23_rows": {}, "k23_cells": {}}
    for c in CELLS:
        rows = _rows({"barred_form": 55.5, "x_e2": BASE_SCORE})
        refs["k18"][c] = [{"seed": r["seed"], "scores": {"barred_form": r["scores"]["barred_form"]}} for r in rows]
        refs["k23_rows"][c] = [{"seed": r["seed"], "scores": {"x_e2": r["scores"]["x_e2"]}} for r in rows]
        refs["k23_cells"][c] = {"mean": float(np.mean([q["scores"]["x_e2"] for q in rows])), "severity_gated": 0}
    return refs


@pytest.mark.parametrize("shift,failed", [
    ({}, []), ({"barred_form": 0.11}, ["barred_form"]), ({"x_e2": -0.11}, ["x_e2"]),
    ({"barred_form": 0.2, "x_e2": 0.2}, ["barred_form", "x_e2"]), ({"barred_form": 0.09, "x_e2": 0.09}, []),
])
def test_the_guards_refuse_off_tolerance_rows(shift, failed):
    rows = {c: _rows({"barred_form": 55.5 + shift.get("barred_form", 0) * (c == ("dm", "g1")),
                      "x_e2": BASE_SCORE + shift.get("x_e2", 0) * (c == ("dm", "g1"))}) for c in CELLS}
    res = k27.guard_eval(rows, _guard_refs())
    assert res["failed"] == failed and res["passes"] == (not failed)
    rows = {c: _rows({"barred_form": 55.5, "x_e2": BASE_SCORE}) for c in CELLS}
    rows[("plain", "g2")][3]["scores"]["x_e2"] += 0.6   # one seed beyond 0.5
    assert k27.guard_eval(rows, _guard_refs())["failed"] == ["x_e2"]
    rows = {c: _rows({"barred_form": 55.5, "x_e2": BASE_SCORE}, {"x_e2": 1} if c == ("dm", "g2") else None)
            for c in CELLS}
    assert k27.guard_eval(rows, _guard_refs())["failed"] == ["x_e2"]   # a gated count off K23's


# -------------------------------------------------------------------------------- (8) the reading

EREC = {"listing_sha256": "f" * 64, "tables": {n: _hex(i) for i, n in enumerate(k27.TABLE_NAMES)},
        "table_manifests": {n: _hex(100 + i) for i, n in enumerate(k27.TABLE_NAMES)},
        "fetch_manifest_sha256": "a" * 64, "kept_sha256": "b" * 64,
        "two": {"listing_sha256": ext.K22_LISTING_SHA256, "combine_sha256": "c" * 64, "tables": {"all": TWO_SHA}}}


@pytest.fixture
def pins(monkeypatch):
    monkeypatch.setattr(k27, "K27_LISTING_SHA256", EREC["listing_sha256"])
    monkeypatch.setattr(k27, "K27_TABLE_SHA256", dict(EREC["tables"]))


def _prov(commit="abc", clean=True, machine="m1"):
    return {"commit": commit, "tree_clean": clean, "machine": machine, "utc": "x"}


def _pre(drop=False) -> dict:
    s = _synth(7, corr_phi=drop)
    return json.loads(json.dumps(k27.pregate_k27(s.phi, s.two, s.tabs, s.e, s.bases, s.den), default=float))


def _write_run(root: Path, cells=None, commit="abc", clean=True, machine="m1", guard=None, drop=None, bad_file=None,
               pre_drop=False, guard_only=False, p1_exit=None, p2_pass=True, p2_only=False):
    cells = cells or {}
    root.mkdir(parents=True)
    if p1_exit is not None:
        (root / "p1.out").write_text("4 22\n")
        (root / "p1.exit").write_text(f"{p1_exit}\n")
    p2 = {"query_kept": 3, "query_bound": 20, "query_K": 20, "query_genes_changed": 20, "pool_size": 10,
          "placebos": {str(s): {"pool_index": s, "kept": 2, "bound": 18, "K": 18, "genes_changed": 18} for s in SEEDS},
          "barred_gene_antigen": 0, "passes": p2_pass, "why": [] if p2_pass else ["the placebo pool holds fewer "
                                                                                    "than 10 antigens"],
          "tables": dict(EREC["tables"]) if p2_pass else {}, "provenance": _prov("reg", True, machine)}
    k18.write_task_file(root / "p2.json", p2)
    ns = SimpleNamespace(root=root, refs=root.parent / "no_refs", k23=root.parent / "no_k23",
                         k18pins=dict(k27.K18_SHA256))
    if p2_only or p1_exit == 3:
        return ns
    pre = _pre(pre_drop)
    consts = k27.constants_of(pre)
    rec = {"arms": list(k27.ARMS), "ext": EREC, "constants": consts, "sizes": dict(k27.SIZES)}
    k18.write_task_file(root / "pregate.json", {"pregate": pre, "arms": list(k27.ARMS), "constants": consts,
                                                "ext": EREC, "p2": {"sha256": _sha(root / "p2.json")},
                                                "provenance": _prov(commit, clean, machine)})
    refs, gsha = {}, {}
    if not pre_drop:
        for part, prefix, part_arms in (("guard", "guard/", k27.GUARDS), ("score", "", k27.ARMS)):
            if part == "score" and guard_only:
                continue
            for w in k27.WORLDS:
                for n in k27.TASKS:
                    a = (harness.FLOOR, *part_arms) if w == "plain" else part_arms
                    rows = _rows({x: cells.get((w, n), _table())[x] for x in a})
                    refs[(w, n)] = [{"seed": r["seed"], "scores": {"barred_form": r["scores"]["barred_form"]}}
                                    for r in rows]
                    for r in rows:
                        s = r["seed"]
                        rel = prefix + (f"plain_{n}/seed_{s}.json" if w == "plain" else f"dm/{n}_seed_{s}.json")
                        f = root / rel
                        f.parent.mkdir(parents=True, exist_ok=True)
                        if drop == rel:
                            continue
                        rr = rec if bad_file != rel else {**rec, "constants": {**consts, "twins": {}}}
                        body = ({"task": {"key": k27.TASKS[n].key}, "rows": [r], "seed": s, "pregate": rr,
                                 "provenance": _prov(commit, clean, machine)} if w == "plain" else
                                {**r, "check": {k18.DM_CHECK: 60.0 + s}, "target": k18.DM_TARGETS[n], "task": n,
                                 "pregate": rr, "arms": list(part_arms), "provenance": _prov(commit, clean, machine)})
                        k18.write_task_file(f, body)
                        if part == "guard":
                            gsha[rel] = _sha(f)
        g = guard or {"passes": True, "failed": []}
        k18.write_task_file(root / "guard.json", {**g, "guard1": {}, "guard2": {}, "files_sha256": gsha,
                                                  "provenance": _prov(commit, clean, machine)})
    k18_dir = root.parent / f"k18_{root.name}"
    k18_dir.mkdir()
    if refs:
        for n in k27.TASKS:
            (k18_dir / f"plain_{n}.json").write_text(json.dumps({"rows": refs[("plain", n)]}))
        (k18_dir / "dm.json").write_text(json.dumps({"rows": {n: refs[("dm", n)] for n in k27.TASKS}}))
    k23_dir = root.parent / f"k23_{root.name}"
    k23_dir.mkdir()
    (k23_dir / "pregate.json").write_text("{}")
    (k23_dir / "reading.json").write_text(json.dumps({"ext": EREC["two"]}))
    ns.refs, ns.k23 = k18_dir, k23_dir
    ns.k18pins = {f: _sha(k18_dir / f) for f in k27.K18_SHA256} if refs else dict(k27.K18_SHA256)
    return ns


@pytest.fixture
def u3d(tmp_path):
    rows = [{"target": t, "seed": s, "world": "depth_matched", "variant": "dm:transfer:mab21l2_all", "t3": 60.0 + s}
            for t in ("rep1", "rep2") for s in SEEDS]
    p = tmp_path / "orth_records.csv"
    pd.DataFrame(rows).to_csv(p, index=False)
    return p


def _read(w, u3d, monkeypatch):
    monkeypatch.setattr(k27, "K18_SHA256", w.k18pins)
    return k27.reading(w.root, w.refs, u3d, w.k23)


def test_a_pass_is_read_with_eight_tries(tmp_path, u3d, monkeypatch, pins):
    w = _write_run(tmp_path / "ok", _everywhere(chip_q_s4=57.0))
    res = _read(w, u3d, monkeypatch)
    assert not res["void"] and res["passes"] and res["named_size"] == "4 S*" and res["tries_counted"] == 8
    assert res["outcome"] == k27.S_PASS.format(size="4 S*")
    t = k27.text(res)
    assert t.splitlines()[0] == "OUTCOME " + res["outcome"]


def test_nothing_qualifies_is_read_with_eight_tries(tmp_path, u3d, monkeypatch, pins):
    res = _read(_write_run(tmp_path / "f"), u3d, monkeypatch)
    assert res["outcome"] == k27.S_FAIL and res["tries_counted"] == 8 and not res["void"]


def test_a_drop_is_read_with_eight_tries(tmp_path, u3d, monkeypatch, pins):
    res = _read(_write_run(tmp_path / "d", pre_drop=True), u3d, monkeypatch)
    assert res["outcome"] == k27.S_DROP and res["tries_counted"] == 8 and res["stop"] == "P3"


@pytest.mark.parametrize("failed,tries", [(["barred_form"], 9), (["x_e2"], 9), (["barred_form", "x_e2"], 10)])
def test_a_failed_guard_is_void_and_adds_each_failed_guard_arm(tmp_path, u3d, monkeypatch, pins, failed, tries):
    w = _write_run(tmp_path / "g", guard={"passes": False, "failed": failed}, guard_only=True)
    res = _read(w, u3d, monkeypatch)
    assert res["void"] and res["tries_counted"] == tries
    assert res["outcome"].startswith("Try K27: VOID: the guard part failed (")
    assert res["outcome"].endswith("; no file; tries counted: 8 (plus the failed guard arms).")


@pytest.mark.parametrize("kw,why", [
    ({"commit": None}, "different commits"),
    ({"clean": False}, "dirty tree"),
    ({"bad_file": "dm/g1_seed_3.json"}, "constants differ"),
])
def test_void_readings_count_eight(tmp_path, u3d, monkeypatch, pins, kw, why):
    if kw.get("commit") is None and "commit" in kw:
        w = _write_run(tmp_path / "v")
        f = w.root / "plain_g2" / "seed_4.json"
        body = json.loads(f.read_text())
        body["provenance"]["commit"] = "other"
        f.write_text(json.dumps(body))
    else:
        w = _write_run(tmp_path / "v", **kw)
    res = _read(w, u3d, monkeypatch)
    assert res["void"] and res["tries_counted"] == 8 and any(why in x for x in res["why"])
    assert res["outcome"].startswith("Try K27: VOID: ")


def test_tables_off_the_pins_void_the_reading(tmp_path, u3d, monkeypatch, pins):
    monkeypatch.setattr(k27, "K27_TABLE_SHA256", {**EREC["tables"], "chipq_Gata4": _hex(999)})
    res = _read(_write_run(tmp_path / "p", _everywhere(chip_q_s4=57.0)), u3d, monkeypatch)
    assert res["void"] and any("pinned" in x for x in res["why"])


def test_an_incomplete_run_is_not_a_reading_and_counts_nothing_new(tmp_path, u3d, monkeypatch, pins, capsys):
    w = _write_run(tmp_path / "i", drop="dm/g2_seed_9.json")
    res = _read(w, u3d, monkeypatch)
    assert res["incomplete"] and "tries_counted" not in res and res["outcome"].startswith("Try K27: not a reading")
    out = tmp_path / "out"
    monkeypatch.setattr(k27, "reading", lambda runs: res)
    with pytest.raises(SystemExit) as ex:
        k27.main(["read", "--runs", str(w.root), "--out", str(out)])
    assert ex.value.code == 3 and not out.exists()


def test_the_stops_at_p1_and_p2_count_nothing(tmp_path, u3d, monkeypatch):
    w = _write_run(tmp_path / "s1", p1_exit=3)
    res = _read(w, u3d, monkeypatch)
    assert res["outcome"] == k27.S_P1 and res["tries_counted"] == 0 and res["p1"] == {"exit": "3", "kept": 4,
                                                                                        "bound": 22}
    w = _write_run(tmp_path / "s2", p1_exit=0, p2_pass=False, p2_only=True)
    res = _read(w, u3d, monkeypatch)
    assert res["outcome"] == k27.S_P2 and res["tries_counted"] == 0
    assert "P2: practice query kept 3" in k27.text(res)
    w = _write_run(tmp_path / "s3", p1_exit=0, p2_only=True)
    assert _read(w, u3d, monkeypatch)["incomplete"]


def test_the_read_writes_once_and_records_its_code(tmp_path, u3d, monkeypatch, pins, capsys):
    w = _write_run(tmp_path / "r", _everywhere(chip_q_s1=57.0))
    monkeypatch.setattr(k27, "K18_SHA256", w.k18pins)
    real = k27.reading
    monkeypatch.setattr(k27, "reading", lambda runs: real(runs, w.refs, u3d, w.k23))
    out = tmp_path / "results"
    k27.main(["read", "--runs", str(w.root), "--out", str(out)])
    first = capsys.readouterr().out
    assert first.splitlines()[0] == "OUTCOME " + k27.S_PASS.format(size="S*")
    body = json.loads((out / "reading.json").read_text())
    assert body["reading_code"]["script"] == k27.SCRIPT and len(body["reading_code"]["chipq_sha256"]) == 64
    with pytest.raises(SystemExit) as ex:
        k27.main(["read", "--runs", str(w.root), "--out", str(out)])
    assert ex.value.code == 4


def test_a_pregate_with_empty_pins_refuses_the_read(tmp_path, u3d, capsys):
    w = _write_run(tmp_path / "e")
    with pytest.raises(k27.Refused, match="pins are not filled"):
        k27.reading(w.root, w.refs, u3d, w.k23)
    with pytest.raises(SystemExit) as ex:
        k27.main(["read", "--runs", str(w.root), "--out", str(tmp_path / "o")])
    assert ex.value.code == 2 and not (tmp_path / "o").exists()


def test_nothing_the_reading_writes_names_a_gene_or_a_placebo(tmp_path, u3d, monkeypatch, pins):
    res = _read(_write_run(tmp_path / "n", _everywhere(chip_q_s4=57.0)), u3d, monkeypatch)
    t = k27.text(res)
    assert not re.search(r"\bg\d{3}\b", t) and not any(re.search(rf"\b{a}\b", t) for a in PLACEBO_NAMES)


# -------------------------------------------------------------------------------- (9) the script's text


def test_the_script_opens_no_target_file_and_names_no_held_out_fact():
    src = (T3 / "k27_chip_query.py").read_text()
    # chipq.py is the only opener of a target file (tests/test_chipq.py: test_only_open_target_reads_a_target_file)
    assert all(w not in src for w in ("checked" "_bytes", "target" "_rel", ".5" ".tsv"))
    assert "rmtree" not in src and "unlink(" not in src
    with contextlib.redirect_stdout(io.StringIO()):
        assert k27.TRY_ARMS and k27.S_VOID.format(clause="x").endswith("(plus the failed guard arms).")


def test_genes_changed_per_arm_and_world_are_reported(tmp_path):
    runs = tmp_path / "runs"
    for s in SEEDS:
        f = runs / "plain_g1" / f"seed_{s}.json"
        f.parent.mkdir(parents=True, exist_ok=True)
        f.write_text(json.dumps({"world_figures": {"x_e2": {"genes_changed": 7},
                                                   "chip_rand_s1": {str(s): {"genes_changed": 10 + s}}}}))
    d = runs / "dm" / "g1_seed_0.json"
    d.parent.mkdir(parents=True)
    d.write_text(json.dumps({"world_figures": {"chip_q_s4": {"genes_changed": 9},
                                               "chip_rand_s4": {str(s): {"genes_changed": 20 + s} for s in SEEDS}}}))
    got = k27.genes_changed(runs)
    assert got["plain"]["x_e2"] == 7 and got["plain"]["chip_rand_s1"] == {str(s): 10 + s for s in SEEDS}
    assert got["dm"]["chip_q_s4"] == 9 and got["dm"]["chip_rand_s4"]["9"] == 29
