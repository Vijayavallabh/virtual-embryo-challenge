"""CHIPQ-S's database part (experiments/t3_ext/chipq_fetch.py and chipq.py) on synthetic files only, no network.

The lists, the target-gene files and the panel are made up here: antigens P01.. are placebo stand-ins, genes g000..
are made-up symbols, and the kept list is a temporary file. Nothing here reads the ledger, the register or any real
download."""

from __future__ import annotations

import csv
import http.client
import importlib.util
import io
import json
import re
import sys
import urllib.request
import urllib.response
import zipfile
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
EXT = ROOT / "experiments" / "t3_ext"
T3 = ROOT / "experiments" / "t3_decomp"
for p in (EXT, T3):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))


def _load(name):
    spec = importlib.util.spec_from_file_location(f"t3_ext_{name}", EXT / f"{name}.py")
    m = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = m
    spec.loader.exec_module(m)
    return m


cq = _load("chipq")
cf = cq.cf   # the very module chipq uses, so a patched DATA_ROOT reaches it

GENES = [f"g{i:03d}" for i in range(80)][::-1]   # the panel, in a file order that is not sorted
CLASS = "Unclassified"


# ---------------------------------------------------------------------------------------------------------------------
# Synthetic database


def _zip(name: str, text: str) -> bytes:
    b = io.BytesIO()
    with zipfile.ZipFile(b, "w") as z:
        z.writestr(name, text)
    return b.getvalue()


def target_file(antigen: str, ids: list[str], rows: dict[str, list[float]], avg="1.0", string="0") -> bytes:
    head = ["Target_genes", f"{antigen}|Average", *(f"{i}|{CLASS}" for i in ids), "STRING"]
    lines = ["\t".join(head)]
    for g, v in rows.items():
        lines.append("\t".join([g, avg if isinstance(avg, str) else avg(g), *(f"{x:g}" for x in v),
                                string if isinstance(string, str) else string(g)]))
    return ("\n".join(lines) + "\n").encode()


def bound_rows(n_bound: int, n_exp: int, n_half: int = 0, seed: int = 0, offset: int = 0) -> dict[str, list[float]]:
    """n_bound panel genes bound in every experiment with distinct means, n_half bound in one experiment only, and two
    genes that are not on the panel."""
    rng = np.random.default_rng(seed)
    rows = {}
    for k in range(n_bound):
        rows[GENES[(offset + k) % len(GENES)]] = list(100.0 + 10 * k + rng.integers(1, 5, n_exp))
    for k in range(n_half):
        rows[GENES[(offset + n_bound + k) % len(GENES)]] = [50.0] + [0.0] * (n_exp - 1)
    rows["NotOnPanelA"] = [7.0] * n_exp
    rows["NotOnPanelB"] = [3.0] * n_exp
    return rows


QUERY_IDS = {"Ctnnb1": ["SRX900", "SRX901", "SRX902", "SRX903"], "Gata4": ["SRX100", "SRX101", "SRX102", "SRX103"]}
KEPT = {"Ctnnb1": ["SRX900", "SRX901", "SRX902"], "Gata4": ["SRX100", "SRX101", "SRX102"]}
# placebo stand-ins: (kept-row count, dropped-row count, bound genes). Gata4's K is 20, so the band is [15, 25].
PLACEBO = {"P01": (2, 1, 15), "P02": (2, 0, 25), "P03": (3, 1, 20), "P04": (2, 0, 14), "P05": (2, 0, 26),
           "P06": (1, 2, 20), "P07": (2, 0, 18), "Gata6": (2, 0, 20),
           **{f"P{n:02d}": (2, 0, 8 + n) for n in range(8, 14)}}   # P08..P13 in the band: a pool of 10
POOL = ["P01", "P02", "P03", "P07", "P08", "P09", "P10", "P11", "P12", "P13"]
CLEAN = "wild type untreated ES cells"
DIRTY = ["Ctnnb1 shRNA knockdown", "treated with an inhibitor", "Mesp1-Cre; flox/flox", "CHIR99021 for 24h"]


def make_db():
    """{url: bytes} of a synthetic archive, plus the files' parts."""
    head = list(cf.ANALYSIS_COLUMNS) + ["Target Genes (TSS ±1k)", "Target Genes (TSS ±5k)",
                                        "Target Genes (TSS ±10k)", "Colocalization (TSV)", "Colocalization (GML)"]
    base = "http://dbarchive.biosciencedbc.jp/kyushu-u/mm10/target/"
    antigens = ["Ctnnb1", "Gata4", *PLACEBO]
    out = io.StringIO()
    w = csv.writer(out, lineterminator="\r\n")
    w.writerow(head)
    for a in antigens:
        w.writerow([a, "", "+", "mm10", f"{base}{a}.1.tsv", f"{base}{a}.5.tsv", f"{base}{a}.10.tsv", "x", "y"])
    w.writerow(["Ctnnb1", "", "+", "hg38", "u1", "u5", "u10", "x", "y"])   # another genome: never read
    analysis = _zip("chip_atlas_analysis_list.csv", out.getvalue())

    exp = io.StringIO()
    w = csv.writer(exp, lineterminator="\r\n")
    w.writerow(cf.EXPERIMENT_COLUMNS)

    def row(i, a, title, genome="mm10"):
        w.writerow([i, genome, "TFs and others", a, "Embryo", "ES cells", "NA", "1,2,3,4", title,
                    "source=" + title, *(["http://dbarchive.biosciencedbc.jp/kyushu-u/mm10/eachData/bw/x.bw"] * 7)])

    for q, ids in QUERY_IDS.items():
        for i in ids:
            row(i, q, CLEAN)
    row("SRX999", "Mab21l2", CLEAN)   # held in the experiment list but not in the analysis list
    n = 1000
    placebo_ids = {}
    for a, (k, d, _) in PLACEBO.items():
        placebo_ids[a] = []
        for _ in range(k):
            row(f"SRX{n}", a, CLEAN)
            placebo_ids[a].append(f"SRX{n}")
            n += 1
        for j in range(d):
            row(f"SRX{n}", a, DIRTY[j % len(DIRTY)])
            placebo_ids[a].append(f"SRX{n}")
            n += 1
    row("SRX5000", "P01", CLEAN, genome="mm9")
    experiments = _zip("chip_atlas_experiment_list.csv", exp.getvalue())

    files = {cf.LISTS["chip_atlas_analysis_list.zip"]: analysis,
             cf.LISTS["chip_atlas_experiment_list.zip"]: experiments}
    https = base.replace("http://", "https://")
    files[f"{https}Ctnnb1.5.tsv"] = target_file("Ctnnb1", QUERY_IDS["Ctnnb1"], bound_rows(16, 4, 3, seed=1))
    files[f"{https}Gata4.5.tsv"] = target_file("Gata4", QUERY_IDS["Gata4"], bound_rows(20, 4, 5, seed=2, offset=7))
    for a, (_, _, nb) in PLACEBO.items():
        files[f"{https}{a}.5.tsv"] = target_file(a, placebo_ids[a], bound_rows(nb, len(placebo_ids[a]), 0,
                                                                                seed=sum(map(ord, a)), offset=nb))
    return files


class FakeNet(urllib.request.HTTPSHandler):
    """A fake https handler inside the fetch's own opener (chipq_fetch._opener), so urllib's redirect and error
    processing run as they would on the network. `redirects` maps a URL to the Location of a 302 answer."""

    def __init__(self, files, redirects=None):
        super().__init__()
        self.files = files
        self.redirects = dict(redirects or {})
        self.calls = []

    def https_open(self, req):
        url = req.full_url
        self.calls.append((url, req.get_header("User-agent")))
        hdrs = http.client.HTTPMessage()
        if url in self.redirects:
            hdrs["Location"] = self.redirects[url]
            r = urllib.response.addinfourl(io.BytesIO(b"moved"), hdrs, url, 302)
            r.msg = "Found"
            return r
        if url not in self.files:
            raise AssertionError(f"unexpected request {url}")
        r = urllib.response.addinfourl(io.BytesIO(self.files[url]), hdrs, url, 200)
        r.msg = "OK"
        return r


_OPENER = cf._opener   # the module's own opener, before any test patches it


def use_net(monkeypatch, f):
    real = _OPENER
    monkeypatch.setattr(urllib.request, "urlopen", lambda *a, **k: pytest.fail("urlopen is never used"))
    monkeypatch.setattr(cf, "_opener", lambda *extra: real(f, *extra))
    return f


@pytest.fixture(autouse=True)
def no_pause(monkeypatch):
    pauses = []
    monkeypatch.setattr(cf, "_sleep", pauses.append)
    return pauses


@pytest.fixture
def net(monkeypatch):
    return use_net(monkeypatch, FakeNet(make_db()))


class _Broken(io.BytesIO):
    """A response body that fails half-way: the connection drops after the first block."""

    def read(self, n=-1):
        if self.tell() > 0:
            raise http.client.IncompleteRead(b"")
        return super().read(min(n, 10) if n and n > 0 else 10)


class FlakyNet(FakeNet):
    """FakeNet whose answers for `flaky` URLs fail: the first `fails` requests drop half-way (a .part is written)."""

    def __init__(self, files, flaky, fails):
        super().__init__(files)
        self.flaky, self.left = set(flaky), dict.fromkeys(flaky, fails)

    def https_open(self, req):
        r = super().https_open(req)
        u = req.full_url
        if u in self.flaky and self.left[u] > 0:
            self.left[u] -= 1
            r = urllib.response.addinfourl(_Broken(self.files[u]), r.headers, u, 200)
            r.msg = "OK"
        return r


@pytest.fixture
def db(tmp_path, net):
    root = tmp_path / "chipatlas"
    cf.step_lists(root)
    cf.step_rows(root)
    cf.step_targets(root)
    return root


def kept_files(tmp_path, kept=None, sha=None):
    import yaml

    yp, sp = tmp_path / "chipq_kept.yaml", tmp_path / "CHIPQ_KEPT.sha256"
    raw = yaml.safe_dump(KEPT if kept is None else kept, sort_keys=True).encode()
    yp.write_bytes(raw)
    sp.write_text(f"{cf.sha256_bytes(raw) if sha is None else sha}  chipq_kept.yaml\n")
    return yp, sp


@pytest.fixture
def kept(tmp_path):
    yp, sp = kept_files(tmp_path)
    return cq.load_kept(yp, sp, check_committed=False)


# ---------------------------------------------------------------------------------------------------------------------
# The fetch


def test_fetch_steps_hosts_agent_counts_and_manifest(tmp_path, net):
    root = tmp_path / "chipatlas"
    assert cf.step_lists(root) == {"files": 2}
    # the specification's three kinds of file only: the two lists and the target-gene files (no README page)
    assert set(cf.LISTS) == {"chip_atlas_analysis_list.zip", "chip_atlas_experiment_list.zip"}
    r = cf.step_rows(root)
    assert r["rows"] == {"Ctnnb1": 4, "Gata4": 4, "Mab21l2": 0}
    # P01 (2 kept), P02, P03 (3), P04, P05, P07, P08..P13, Gata6 keep at least two; P06 keeps one
    fetched = [a for a in PLACEBO if a != "P06"]
    assert r["placebo_at_least_two"] == len(fetched) == 13
    t = cf.step_targets(root)
    assert (t["query_files"], t["placebo_files"]) == (2, 13)
    assert t["total_bytes"] == sum((root / "targets" / f"{a}.5.tsv").stat().st_size
                                   for a in ["Ctnnb1", "Gata4", *fetched])
    assert not (root / "targets" / "P06.5.tsv").exists() and not (root / "screen" / "Mab21l2.rows.tsv").exists()
    for url, agent in net.calls:
        assert url.startswith("https://") and re.match(r"https://(dbarchive\.biosciencedbc\.jp|chip-atlas\.dbcls\.jp)/",
                                                         url)
        assert agent == cf.USER_AGENT
        assert not re.search(r"bigwig|\.bw$|\.bed|bigbed|\.bb$|\.1\.tsv|\.10\.tsv", url, re.IGNORECASE)
    assert not re.search(r"embryo|vec|mithril|team|challenge|claude", cf.USER_AGENT, re.IGNORECASE)
    kinds = re.compile(r"/chip_atlas_(analysis|experiment)_list\.zip$|/target/[^/]+\.5\.tsv$")
    assert all(kinds.search(u) for u, _ in net.calls)
    man = json.loads((root / "MANIFEST.json").read_text())
    for rel, rec in man["files"].items():
        assert {"url", "size", "sha256", "fetched"} <= set(rec)
        assert (root / rel).stat().st_size == rec["size"]
    assert "screen/placebo_kept.json" in man["files"] and "screen/Gata4.rows.tsv" in man["files"]
    assert cf.step_verify(root)
    # the readers' file holds the antigen's rows only, every field, tab-separated, with the header
    lines = (root / "screen" / "Gata4.rows.tsv").read_text().splitlines()
    assert lines[0].split("\t") == list(cf.EXPERIMENT_COLUMNS)
    assert sorted(x.split("\t")[0] for x in lines[1:]) == QUERY_IDS["Gata4"]
    # a rerun makes no request and changes nothing
    n = len(net.calls)
    cf.step_lists(root), cf.step_rows(root), cf.step_targets(root)
    assert len(net.calls) == n


def test_placebo_kept_is_the_keyword_filter_on_mm10_rows(db):
    pk = cf.placebo_kept(db)
    assert "Ctnnb1" not in pk and "Gata4" not in pk and "Mab21l2" not in pk
    assert len(pk["P01"]) == 2 and len(pk["P03"]) == 3 and len(pk["P06"]) == 1
    assert all(v == sorted(v) for v in pk.values())
    assert "SRX5000" not in pk["P01"]   # an mm9 row


def test_fetch_refuses_a_changed_file_and_overwrites_nothing(db, net):
    p = db / "lists" / "chip_atlas_experiment_list.zip"
    raw = p.read_bytes()
    p.write_bytes(raw[:-1] + bytes([raw[-1] ^ 1]))
    assert not cf.step_verify(db)
    with pytest.raises(cf.Refused):
        cf.step_rows(db)
    with pytest.raises(cf.Refused):
        cf.step_lists(db)
    assert p.read_bytes() != raw   # not re-fetched over the changed file
    p.write_bytes(raw)
    assert cf.step_verify(db)
    t = db / "targets" / "P02.5.tsv"
    t.write_bytes(t.read_bytes() + b"x")
    assert main_exit(cf.main, ["verify"], db) == 2


def main_exit(fn, argv, root, monkeypatch=None):
    old = cf.DATA_ROOT
    cf.DATA_ROOT = root
    try:
        return fn(argv)
    finally:
        cf.DATA_ROOT = old


def test_fetch_sets_leftover_parts_aside_and_refuses_unknown_files(tmp_path, net):
    root = tmp_path / "chipatlas"
    (root / "lists").mkdir(parents=True)
    (root / "lists" / "chip_atlas_analysis_list.zip.part").write_bytes(b"half")
    (root / "lists" / "chip_atlas_analysis_list.zip.failed.1").write_bytes(b"older")
    assert cf.step_lists(root) == {"files": 2}   # no removal by hand
    assert (root / "lists" / "chip_atlas_analysis_list.zip.failed.1").read_bytes() == b"older"
    assert (root / "lists" / "chip_atlas_analysis_list.zip.failed.2").read_bytes() == b"half"
    assert not (root / "lists" / "chip_atlas_analysis_list.zip.part").exists()
    assert cf.step_verify(root)   # the set-aside files are refused nowhere and never in the manifest
    assert not any(".failed." in rel for rel in cf.load_manifest(root)["files"])
    root2 = tmp_path / "chipatlas2"
    (root2 / "lists").mkdir(parents=True)
    (root2 / "lists" / "chip_atlas_analysis_list.zip").write_bytes(b"not from the manifest")
    with pytest.raises(cf.Refused):
        cf.step_lists(root2)
    assert (root2 / "lists" / "chip_atlas_analysis_list.zip").read_bytes() == b"not from the manifest"


P03_URL = "https://dbarchive.biosciencedbc.jp/kyushu-u/mm10/target/P03.5.tsv"


def test_a_failed_download_is_retried_and_its_part_set_aside(tmp_path, monkeypatch, no_pause):
    files = make_db()
    f = use_net(monkeypatch, FlakyNet(files, [P03_URL], fails=2))
    root = tmp_path / "chipatlas"
    cf.step_lists(root), cf.step_rows(root), cf.step_targets(root)
    t = root / "targets"
    assert (t / "P03.5.tsv").read_bytes() == files[P03_URL] and cf.step_verify(root)
    assert sorted(p.name for p in t.glob("P03.5.tsv.*")) == ["P03.5.tsv.failed.1", "P03.5.tsv.failed.2"]
    assert [u for u, _ in f.calls].count(P03_URL) == 3
    # a short pause before every request, a longer one before each retry
    assert no_pause.count(cf.RETRY_PAUSE) == 2 and no_pause.count(cf.REQUEST_PAUSE) == len(f.calls) - 2
    assert cf.REQUEST_PAUSE > 0 and cf.RETRY_PAUSE > 0 and cf.RETRIES == 3


def test_a_download_failing_every_try_is_refused_and_a_rerun_needs_no_hand(tmp_path, monkeypatch, capsys):
    files = make_db()
    f = use_net(monkeypatch, FlakyNet(files, [P03_URL], fails=4))
    root = tmp_path / "chipatlas"
    monkeypatch.setattr(cf, "DATA_ROOT", root)
    assert cf.main(["lists"]) == 0 and cf.main(["rows"]) == 0
    capsys.readouterr()
    assert cf.main(["targets"]) == 2
    err = capsys.readouterr().err
    assert "failed 4 times" in err and "IncompleteRead" in err and "P03" not in err
    assert [u for u, _ in f.calls].count(P03_URL) == 4
    t = root / "targets"
    assert len(list(t.glob("P03.5.tsv.failed.*"))) == 4 and not (t / "P03.5.tsv").exists()
    assert not list(t.glob("*.part")) and "targets/P03.5.tsv" not in cf.load_manifest(root)["files"]
    # the network is back: the rerun fetches the file and keeps every set-aside file as it was
    use_net(monkeypatch, FakeNet(files))
    assert cf.main(["targets"]) == 0 and cf.main(["verify"]) == 0
    assert (t / "P03.5.tsv").read_bytes() == files[P03_URL] and len(list(t.glob("P03.5.tsv.failed.*"))) == 4


def test_a_refusal_is_not_retried(tmp_path, monkeypatch):
    f = use_net(monkeypatch, FakeNet(make_db(), redirects={cf.LISTS["chip_atlas_analysis_list.zip"]: "/x"}))
    with pytest.raises(cf.Refused, match="redirect"):
        cf.step_lists(tmp_path / "chipatlas")
    assert len(f.calls) == 1
    # over the byte limit: refused at once, the partial file set aside (not left for removal by hand)
    f = use_net(monkeypatch, FakeNet(make_db()))
    monkeypatch.setattr(cf, "MAX_BYTES", {"list": 5, "target": 5})
    root = tmp_path / "chipatlas2"
    with pytest.raises(cf.Refused, match="byte limit"):
        cf.step_lists(root)
    assert len(f.calls) == 1 and not list(root.rglob("*.part"))
    assert [p.name for p in (root / "lists").iterdir()] == ["chip_atlas_analysis_list.zip.failed.1"]


def test_urls_are_refused_off_the_two_hosts_and_for_signal_files():
    ok = cf.https_url("http://dbarchive.biosciencedbc.jp/kyushu-u/mm10/target/P01.5.tsv")
    assert ok == "https://dbarchive.biosciencedbc.jp/kyushu-u/mm10/target/P01.5.tsv"
    cf.https_url("https://chip-atlas.dbcls.jp/data/mm10/target/P01.5.tsv")
    for bad in ("https://example.org/a.tsv", "https://dbarchive.biosciencedbc.jp:8080/a.tsv",
                "https://dbarchive.biosciencedbc.jp/a.tsv?x=1", "ftp://dbarchive.biosciencedbc.jp/a.tsv",
                "https://dbarchive.biosciencedbc.jp/kyushu-u/mm10/eachData/bw/SRX1.bw",
                "https://chip-atlas.dbcls.jp/data/mm10/eachData/bed05/SRX1.05.bed",
                "https://chip-atlas.dbcls.jp/data/mm10/eachData/bb05/SRX1.05.bb"):
        with pytest.raises(cf.Refused):
            cf.https_url(bad)


@pytest.mark.parametrize("location", ["https://example.org/x", "http://example.org/chip_atlas_analysis_list.zip",
                                      "https://dbarchive.biosciencedbc.jp/data/chip-atlas/other/x.zip", "/elsewhere"])
def test_a_redirect_is_refused_with_no_part_and_no_manifest_entry(tmp_path, monkeypatch, location):
    files = make_db()
    url = cf.LISTS["chip_atlas_analysis_list.zip"]
    f = use_net(monkeypatch, FakeNet(files, redirects={url: location}))
    root = tmp_path / "chipatlas"
    with pytest.raises(cf.Refused, match="redirect"):
        cf.step_lists(root)
    assert [u for u, _ in f.calls] == [url]   # the redirect was never followed
    assert not list(root.rglob("*.part")) and not (root / "lists" / "chip_atlas_analysis_list.zip").exists()
    assert not (root / "MANIFEST.json").exists() or cf.load_manifest(root)["files"] == {}
    assert main_exit(cf.main, ["lists"], root) == 2


def test_an_answer_from_another_url_is_refused(tmp_path, monkeypatch):
    files = make_db()

    class Moved(FakeNet):
        def https_open(self, req):
            r = super().https_open(req)
            r.url = "https://example.org/x"   # as if a handler had followed a redirect
            return r

    use_net(monkeypatch, Moved(files))
    root = tmp_path / "chipatlas"
    with pytest.raises(cf.Refused, match="another URL"):
        cf.step_lists(root)
    assert not list(root.rglob("*.part")) and not (root / "MANIFEST.json").exists()


def test_outputs_stay_outside_the_checkout():
    with pytest.raises(cf.Refused):
        cf.root_of(ROOT / "data" / "chipatlas")


def test_analysis_list_url_must_name_its_own_antigen(tmp_path):
    head = list(cf.ANALYSIS_COLUMNS) + ["a", "Target Genes (TSS ±5k)", "b", "c", "d"]
    base = "http://dbarchive.biosciencedbc.jp/kyushu-u/mm10/target/"
    out = io.StringIO()
    w = csv.writer(out)
    w.writerow(head)
    w.writerow(["P01", "", "+", "mm10", "x", f"{base}P02.5.tsv", "x", "x", "x"])
    with pytest.raises(cf.Refused):
        cf.analysis_rows(_zip("a.csv", out.getvalue()))


# ---------------------------------------------------------------------------------------------------------------------
# The keyword filter


DROP_EXAMPLES = {
    "mutant": "Ctcf mutant ES cells", "flox": "Smad4 floxed liver", "Cre": "Mesp1-Cre heart",
    "knockout": "Tet2 knock-out bone marrow", "knockdown": "Suz12 knockdown", "sh/si/shRNA": "shRNA against Ezh2",
    "CRISPR": "CRISPR-edited clone", "transgene": "transgenic line", "reporter": "GFP reporter mice",
    "overexpression": "Sox2 over-expression", "mutant (genotype notation)": "Dnmt1 -/- cells",
    "stabilised": "stabilized allele", "activated": "activated macrophages", "exon-3-deleted": "exon 3 deleted allele",
    "constitutively active": "constitutively active form", "Wnt-stimulated / Wnt3a / Wnt-ligand": "Wnt3a conditioned",
    "stimulated": "LPS stimulated", "GSK3-inhibitor": "GSK-3 blocked", "CHIR": "CHIR99021 3uM", "LiCl": "LiCl 20mM",
    "a medium naming a pathway drug (2i)": "ES cells in 2i/LIF", "drug": "drug exposure", "inhibitor": "Mek inhibitor",
    "agonist": "RAR agonist", "ligand": "Notch ligand", "growth-factor": "growth factor starved",
    "treatment": "treated with retinoic acid",
}


def test_every_drop_class_word_drops():
    words = [w for _, w, _ in cf.DROP_PATTERNS]
    assert sorted(words) == sorted(DROP_EXAMPLES)
    for w, text in DROP_EXAMPLES.items():
        assert f"{[c for c, ww, _ in cf.DROP_PATTERNS if ww == w][0]}:{w}" in cf.dropped_by(["NA", text]), w
    assert {c for c, _, _ in cf.DROP_PATTERNS} == {"i", "ii", "iii"}   # class iv is "the same list"


@pytest.mark.parametrize("text", ["untreated", "wild type", "Wild type", "wild-type", "WT untreated E14 cells",
                                  "Untreated wild type ES cells", "neural crest", "Six1 ChIP", "SHH ChIP", "Creb1",
                                  "Crebbp", "Crem", "CREB1", "Kdm5b", "Kdm1a", "Tokyo", "secreted", "+/+", "SKOV3"])
def test_plain_words_do_not_drop(text):
    assert cf.keep_row(["TFs and others", "Pax3", "Embryo", "ES cells", text, text]), cf.dropped_by([text])


# the reviewer's cases (blocker 1): inflections, the KO/KD abbreviation with a prefix or a plural, fl/fl and OE straight
# after a symbol's digit, the genotype notation with a Unicode minus or an en dash, Cre joined to a driver name
VARIANT_DROPS = {
    "knockout": ["knocked out", "cKO heart", "Gata4cKO", "dKO cells", "iKO", "KOs", "knock-outs", "SmadKO", "KO mice"],
    "knockdown": ["knocked down", "knocking down", "cKD", "Gata4KD", "KDs"],
    "flox": ["Ctnnb1fl/fl", "fl/fl", "Gata4 fl / fl", "fl/+"],
    "overexpression": ["Gata4OE", "OE line"],
    "mutant (genotype notation)": ["Gata4−/−", "Gata4 –/–", "Gata4+/−", "het +/-"],
    "Cre": ["AlbCre", "Mesp1Cre heart", "Vil1cre", "LysMcre", "iCre", "CreERT2", "CRE line", "Mesp1-Cre"],
}


@pytest.mark.parametrize("word,text", [(w, t) for w, ts in VARIANT_DROPS.items() for t in ts])
def test_a_doubt_drops_variants(word, text):
    assert f"i:{word}" in cf.dropped_by([text]), (text, cf.dropped_by([text]))


def test_a_doubt_drops_other_variants():
    for t in ("Non-targeting siRNA", "siGata", "sh-Ctrl", "Inhibition of", "pre-treated"):
        assert not cf.keep_row([t]), t


# the rule agent's A13 addendum (kwdecide): spellings, notations and abbreviations of the paragraph's own words, each
# hitting the pattern of the word it belongs to
KWDECIDE_DROPS = {
    "i:transgene": ["Tg(Myh6-cre) line", "Tg/+", "tg/+ embryos", "Tg"],
    "i:sh/si/shRNA": ["esiRNA", "esiRNA pool", "shctrl", "sictrl", "shscr", "sh_ctrl", "si-control", "shScrambled"],
    "i:flox": ["fl/-", "Gata4 fl/\u2212", "fl / \u2013", "f/f", "Gata4 f/f"],
    "i:mutant (genotype notation)": ["-/+", "Gata4 \u2212/+", "het \u2013 / +"],
    "ii:exon-3-deleted": ["ex3", "delta-ex3", "\u0394ex3", "dex3", "del(ex3)", "Ex3 allele", "ex-3"],
    "i:mutant": ["ENU-mutagenized line", "mutagenized", "S33Y", "p.V600E knock-in", "R270* allele", "G12D"],
    "i:CRISPR": ["Cas9", "dCas9-KRAB", "SpCas9 clone", "cas-9"],
}


@pytest.mark.parametrize("hit,text", [(h, t) for h, ts in KWDECIDE_DROPS.items() for t in ts])
def test_kwdecide_forms_drop_by_their_word(hit, text):
    assert hit in cf.dropped_by([text]), (text, cf.dropped_by([text]))


@pytest.mark.parametrize("text", ["sgRNA", "sgRNA against Pax3", "sgRNAs", "Tgfb1", "Tgif1", "E14Tg2a", "SHC1",
                                  "Sicily", "Shh", "Dex 3uM", "Exo3", "Ext1", "K562", "S33", "+/+", "fl"])
def test_kwdecide_not_matched_and_look_alikes_are_kept(text):
    # sgRNA alone is not a form of CRISPR (A13: a named agent without the challenger's word is kept); the rest are
    # symbols or plain names near a new form that spell no drop word
    assert cf.keep_row([text]), (text, cf.dropped_by([text]))


@pytest.mark.parametrize("text", ["treatment: none", "treatment=untreated", "Treatment: No treatment", "no treatment",
                                  "treatments = none", "treatment:none", "TREATMENT: NONE ",
                                  "source=liver|treatment=none|strain=C57BL/6", "strain: B6; treatment: untreated"])
def test_a_field_whose_whole_value_states_no_treatment_is_kept(text):
    assert cf.keep_row([text]), (text, cf.dropped_by([text]))


@pytest.mark.parametrize("text", ["treatment: vehicle", "treatment: DMSO", "treatment: control", "treatment: n/a",
                                  "treatment:", "treatment: ", "treatment=none, then RA", "pre-treatment: none",
                                  "treatment: none|agent: treated", "treatment: no treatment with RA",
                                  "no treatment then RA", "treatment: none (RA 2 days later)"])
def test_only_the_whole_no_treatment_value_is_exempt(text):
    assert "iii:treatment" in cf.dropped_by([text]), (text, cf.dropped_by([text]))


def test_a_named_agent_without_a_class_word_is_kept():
    # points_settled A13 and B8: the filter matches the challenger's words and their variants only; it adds no word
    for t in ("tamoxifen 1uM", "dox 2ug/ml", "Activin A", "4-OHT"):
        assert cf.keep_row([t]), (t, cf.dropped_by([t]))


def test_made_up_symbols_trip_no_pattern_unless_they_spell_a_drop_word():
    """An antigen's own symbol is a screened field: made-up symbols (a capital, one or two lowercase letters, an
    optional suffix) drop only when the symbol itself spells a drop word (Cre, Ko, Kd, Oe, Wnt, LiCl, and since the
    A13 addendum Tg, the abbreviation of transgene: Tg, Tg1, Tg2a)."""
    low = "abcdefghijklmnopqrstuvwxyz"
    syms = [c.upper() + a + b + s for c in low for a in ("", *low) for b in low for s in ("", "1", "2a", "l1")]
    dropped = sorted(s for s in syms if not cf.keep_row([s]))
    assert all(re.fullmatch(r"(?:Cre|Kd|Ko|Oe|Wnt|Licl|Tg)[0-9a-z]*", s) for s in dropped), dropped
    assert len(dropped) == 20


def test_the_filter_reads_the_metadata_fields_only():
    f = ["SRX1", "mm10", "TFs and others", "P01", "Embryo", "ES cells", "NA", "1,2,3,4", CLEAN, CLEAN,
         *(["http://dbarchive.biosciencedbc.jp/kyushu-u/mm10/eachData/bed05/x.bed"] * 7)]
    assert cf.keep_row(cf.screened_fields(f))
    f[8] = "shRNA"
    assert not cf.keep_row(cf.screened_fields(f))


# ---------------------------------------------------------------------------------------------------------------------
# The kept list


def test_kept_list_absent_disagreeing_or_uncommitted_is_refused(tmp_path):
    yp, sp = kept_files(tmp_path)
    assert cq.load_kept(yp, sp, check_committed=False)["Gata4"] == tuple(KEPT["Gata4"])
    with pytest.raises(cq.Refused, match="not committed"):
        cq.load_kept(yp, sp)
    with pytest.raises(cq.Refused, match="absent"):
        cq.load_kept(yp, tmp_path / "missing.sha256", check_committed=False)
    with pytest.raises(cq.Refused, match="absent"):
        cq.load_kept(tmp_path / "missing.yaml", sp, check_committed=False)
    (tmp_path / "x").mkdir()
    yp2, sp2 = kept_files(tmp_path / "x", sha="0" * 64)
    with pytest.raises(cq.Refused, match="does not hash"):
        cq.load_kept(yp2, sp2, check_committed=False)


def test_kept_list_form_is_checked(tmp_path):
    for n, bad in enumerate(({"Ctnnb1": ["SRX2", "SRX1"]}, {"P01": ["SRX1"]}, {"Gata4": ["notanid"]},
                             {"Gata4": "SRX1"})):
        d = tmp_path / f"case{n}"
        d.mkdir()
        yp, sp = kept_files(d, bad)
        with pytest.raises(cq.Refused):
            cq.load_kept(yp, sp, check_committed=False)


def test_committing_the_kept_list_changes_no_code():
    src = (EXT / "chipq.py").read_text()
    assert 'KEPT_YAML = HERE / "chipq_kept.yaml"' in src and 'KEPT_SHA = HERE / "CHIPQ_KEPT.sha256"' in src
    assert not re.search(r"KEPT_SHA256\s*=\s*\"[0-9a-f]{64}\"", src)


# ---------------------------------------------------------------------------------------------------------------------
# The one reader


def test_an_experiment_not_in_the_pinned_kept_list_is_refused(db, kept):
    with pytest.raises(cq.Refused, match="not in the pinned kept list"):
        cq.open_target("Gata4", GENES, db, kept, ["SRX103"])
    with pytest.raises(cq.Refused, match="not in the pinned kept list"):
        cq.open_target("Gata4", GENES, db, kept, ["SRX100", "SRX999"])
    with pytest.raises(cq.Refused, match="load_kept"):
        cq.open_target("Gata4", GENES, db, {"Gata4": ["SRX103"]})
    with pytest.raises(cq.Refused, match="placebo_kept.json only"):
        cq.open_target("P01", GENES, db, kept)
    E, ids = cq.open_target("Gata4", GENES, db, kept)
    assert ids == KEPT["Gata4"] and E.shape == (len(GENES), 3)


def test_a_kept_id_without_a_column_is_refused(db, tmp_path):
    (tmp_path / "k").mkdir()
    yp, sp = kept_files(tmp_path / "k", {"Gata4": ["SRX100", "SRX777"]})
    with pytest.raises(cq.Refused, match="no column"):
        cq.open_target("Gata4", GENES, db, cq.load_kept(yp, sp, check_committed=False))


def test_a_placebo_uses_its_kept_ids_that_have_a_column(db, kept):
    """Review blocker 5: a placebo's filter-kept id without a column in its target file is left out, not refused; the
    pool's two-experiment bar counts the columned ids; a query antigen keeps the refusal (test above)."""
    pk = cf.placebo_kept(db)
    k03, k07 = pk["P03"], pk["P07"]
    assert len(k03) == 3 and len(k07) == 2
    # P03 loses one kept column (2 remain: stays in the pool); P07 loses one (1 remains: leaves the pool)
    _rewrite_target(db, "P03", target_file("P03", k03[1:], bound_rows(20, 2, 0, seed=3, offset=20)))
    _rewrite_target(db, "P07", target_file("P07", k07[:1], bound_rows(18, 1, 0, seed=4, offset=18)))
    E, ids = cq.open_target("P03", GENES, db)
    assert ids == k03[1:] and E.shape == (len(GENES), 2)
    E, ids = cq.open_target("P07", GENES, db)
    assert ids == k07[:1] and E.shape == (len(GENES), 1)
    r = cq.build("P03", GENES, db)
    assert r["n_kept"] == 2 and r["n_bound"] == 20
    with pytest.raises(cq.Refused, match="no column"):   # asked for by name: still refused
        cq.open_target("P03", GENES, db, None, [k03[0]])
    rq = cq.build("Gata4", GENES, db, kept)
    pool = cq.placebo_pool("Gata4", rq["K"], GENES, db)
    assert "P03" in pool and "P07" not in pool
    # no kept column at all: E has no column, and the antigen is out of the pool
    _rewrite_target(db, "P03", target_file("P03", ["SRX4242"], bound_rows(20, 1, 0, seed=3, offset=20)))
    E, ids = cq.open_target("P03", GENES, db)
    assert ids == [] and E.shape == (len(GENES), 0)
    assert "P03" not in cq.placebo_pool("Gata4", rq["K"], GENES, db)


def _rewrite_target(root, antigen, raw):
    """Replace a target file and its manifest entry (the test stands in for a fresh fetch)."""
    p = root / "targets" / f"{antigen}.5.tsv"
    p.write_bytes(raw)
    man = json.loads((root / "MANIFEST.json").read_text())
    man["files"][f"targets/{antigen}.5.tsv"].update(size=len(raw), sha256=cf.sha256_bytes(raw))
    (root / "MANIFEST.json").write_text(json.dumps(man))


def test_poisoned_average_string_and_dropped_columns_change_nothing(db, kept):
    ref = cq.build("Gata4", GENES, db, kept)
    rows = bound_rows(20, 4, 5, seed=2, offset=7)
    for avg, string in (("NaN", "garbage"), ("-1e300", "xéy"), (lambda g: "9999", lambda g: "Ctnnb1")):
        _rewrite_target(db, "Gata4", target_file("Gata4", QUERY_IDS["Gata4"], rows, avg=avg, string=string))
        r = cq.build("Gata4", GENES, db, kept)
        assert np.array_equal(r["C"], ref["C"]) and r["n_bound"] == ref["n_bound"]
    # the dropped experiment (the fourth column) poisoned: not a number, or a huge score on every gene
    for v in ("nope", "1e9"):
        lines = target_file("Gata4", QUERY_IDS["Gata4"], rows).decode().splitlines()
        out = [lines[0]] + ["\t".join(f[:5] + [v] + f[6:]) for f in (x.split("\t") for x in lines[1:])]
        _rewrite_target(db, "Gata4", ("\n".join(out) + "\n").encode())
        r = cq.build("Gata4", GENES, db, kept)
        assert np.array_equal(r["C"], ref["C"])
    # a kept column poisoned is refused, and the message names no gene, score or column
    lines = target_file("Gata4", QUERY_IDS["Gata4"], rows).decode().splitlines()
    out = [lines[0]] + ["\t".join(f[:2] + ["oops"] + f[3:]) for f in (x.split("\t") for x in lines[1:])]
    _rewrite_target(db, "Gata4", ("\n".join(out) + "\n").encode())
    with pytest.raises(cq.Refused) as e:
        cq.build("Gata4", GENES, db, kept)
    assert not re.search(r"g\d{3}|SRX\d+|oops|\d{3}", str(e.value))


def test_the_reader_refuses_a_changed_file(db, kept):
    p = db / "targets" / "Gata4.5.tsv"
    p.write_bytes(p.read_bytes().replace(b"\t1.0\t", b"\t2.0\t", 1))
    with pytest.raises(cq.Refused, match="differ from the manifest"):
        cq.open_target("Gata4", GENES, db, kept)


# The read script (and its tests) may CALL chipq's functions (open_target, build, placebo_pool, ...), but nothing
# outside chipq.py opens a target-gene file itself: no target file name or path, no target_rel, no checked_bytes.
READ_CALLERS = ("experiments/t3_decomp/k27_chip_query.py", "tests/test_k27_chip_query.py")
OPENS_A_TARGET = re.compile(r"\.5\.tsv|target_rel|checked_bytes|[\"']targets/|/\s*[\"']targets[\"']")
NAMES_THE_READER = re.compile(r"chipatlas|open_target")


def target_access(rel: str, text: str) -> bool:
    """True when the file at `rel` (repository-relative) reaches a target-gene file in a way the rule forbids."""
    if OPENS_A_TARGET.search(text):
        return True
    return rel not in READ_CALLERS and bool(NAMES_THE_READER.search(text))


def test_the_static_rule_lets_the_read_script_call_chipq_but_not_open_a_file():
    calls = "import chipq as cq\nE, ids = cq.open_target('Gata4', genes)\nr = cq.build('Gata4', genes)\n"
    assert not target_access(READ_CALLERS[0], calls) and not target_access(READ_CALLERS[1], calls)
    assert target_access("experiments/t3_decomp/k24_wider_tables.py", calls)   # any other file: not even a call
    for opens in ("p = root / 'targets' / f'{a}.5.tsv'", "cq.cf.checked_bytes(root, cq.cf.target_rel(a))",
                  "open(root / \"targets/\" + a)", "raw = cf.checked_bytes(root, rel)"):
        assert target_access(READ_CALLERS[0], opens), opens


def test_only_open_target_reads_a_target_file():
    src = (EXT / "chipq.py").read_text()
    assert src.count("checked_bytes(") == 1
    i = src.index("checked_bytes(")
    assert src.rfind("def ", 0, i) == src.index("def open_target")
    hits = []
    for base in (ROOT / "experiments", ROOT / "src", ROOT / "tests"):
        for f in base.rglob("*.py"):
            if f.name in ("chipq.py", "chipq_fetch.py", "test_chipq.py"):
                continue
            if target_access(f.relative_to(ROOT).as_posix(), f.read_text(errors="replace")):
                hits.append(str(f.relative_to(ROOT)))
    assert hits == []
    assert "STRING" not in (EXT / "chipq_fetch.py").read_text()


# ---------------------------------------------------------------------------------------------------------------------
# Steps 3 to 7


def test_bound_set_weights_and_table():
    E = np.array([[5, 5, 0], [9, 0, 0], [1, 1, 1], [0, 0, 0], [3, 3, 0], [3, 3, 0]], dtype=float)
    s, bound, order = cq.bound_rank(E)
    assert bound.tolist() == [True, False, True, False, True, True]   # 2 * hits >= 3
    assert order.tolist() == [0, 4, 5, 2]                               # by s_g, ties in file order
    w, K = cq.weights(order, 6)
    assert K == 4 and w.tolist() == [1.0, 0.0, 0.5, 0.0, 1 - 0.5 / 3, 1 - 1.0 / 3]
    C = cq.raw_table(w)
    assert np.all(C[w == 0] == 1.0) and np.allclose(C[w > 0], np.exp(-w[w > 0]), rtol=0, atol=0)
    w1, K1 = cq.weights(np.array([3]), 6)
    assert K1 == 1 and w1[3] == 1.0 and w1.sum() == 1.0
    w0, K0 = cq.weights(np.array([], dtype=int), 6)
    assert K0 == 0 and np.all(cq.raw_table(w0) == 1.0) and cq.exponent(cq.raw_table(w0), cq.S_STAR) is None
    order = np.arange(60)
    w, K = cq.weights(order, 80)
    assert K == 40 and w[0] == 1.0 and w[39] == 0.5 and np.all(w[40:] == 0)


def test_the_table_is_exactly_one_at_unbound_genes(db, kept):
    r = cq.build("Gata4", GENES, db, kept)
    assert r["K"] == 20 and r["n_bound"] == 20 and r["n_kept"] == 3
    assert np.all(r["C"][r["w"] == 0] == 1.0) and (r["C"] != 1.0).sum() == 20
    assert np.all((r["C"] >= np.exp(-1)) & (r["C"] <= np.exp(-0.5)) | (r["C"] == 1.0))


def test_the_exponent_reaches_the_size_at_both_sizes(db, kept):
    import k24_wider_tables as k24

    assert cq.S_STAR == k24.S_STAR and cq.SIZES == {"S*": k24.S_STAR, "4S*": 4 * k24.S_STAR}
    for a, kp in (("Gata4", kept), ("P03", None)):
        r = cq.build(a, GENES, db, kp)
        for name, S in cq.SIZES.items():
            e = r["exponents"][name]
            assert abs(float(np.std(e * np.log(r["C"]))) - S) < 1e-12
            assert abs(float(np.std(np.log(r["C"] ** e))) - S) < 1e-12


# ---------------------------------------------------------------------------------------------------------------------
# The placebo pool


def test_placebo_pool_band_exclusions_and_reproducible_draw(db, kept):
    rq = cq.build("Gata4", GENES, db, kept)
    pool = cq.placebo_pool("Gata4", rq["K"], GENES, db)
    # band [15, 25] inclusive: P01 (15), P02 (25), P03 (20), P07 (18), P08..P13 (16..21); P04 (14) and P05 (26) out;
    # P06 keeps one experiment; Gata6 is excluded whatever its size
    assert pool == POOL
    assert cq.placebo_pool("P03", 20, GENES, db) == [a for a in POOL if a != "P03"]


def test_the_placebo_draw_is_sequential_without_replacement():
    """Settled point A1: seeds 0..9 in order; for seed s, remaining[default_rng(20261124 + s).integers(len(remaining))]
    over the sorted pool minus earlier draws; ten different antigens; a pool under 10 refuses."""
    from vec.adoption import SEEDS

    pool = [f"A{n:02d}" for n in range(14)]
    draws = cq.placebo_draws(pool)
    remaining, expect = sorted(pool), []
    for s in range(10):   # written out here, independently of the module
        a = remaining.pop(int(np.random.default_rng(20261124 + s).integers(len(remaining))))
        expect.append((s, a, sorted(pool).index(a)))
    assert draws == expect and [d[0] for d in draws] == list(SEEDS)
    assert len({a for _, a, _ in draws}) == 10
    for s, a, i in draws:
        assert cq.placebo_for_seed(pool, s) == (a, i) and pool[i] == a
    assert cq.placebo_draws(pool) == draws   # reproducible
    # exactly 10: every antigen drawn once
    assert sorted(a for _, a, _ in cq.placebo_draws(pool[:10])) == pool[:10]
    for bad in (pool[:9], [], pool[::-1], pool[:10] + ["A00"]):
        with pytest.raises(cq.Refused):
            cq.placebo_draws(bad)
    with pytest.raises(cq.Refused, match="fewer than 10"):
        cq.placebo_for_seed(pool[:9], 0)
    with pytest.raises(cq.Refused, match="SEEDS"):
        cq.placebo_for_seed(pool, 10)


# ---------------------------------------------------------------------------------------------------------------------
# Gene order through the loader K24's code uses, and the tables step


def _barred_run(tmp_path, genes):
    import build_submission as bs

    run = tmp_path / "run"
    (run / "factors").mkdir(parents=True)
    p = run / "factors" / f"{cq.BARRED_NAME}.csv"
    p.write_text("gene,factor\n" + "".join(f"{g},1.0\n" for g in genes))
    rec = {"files": {cq.BARRED_NAME: {"sha256": cf.sha256_file(p)}}}
    (run / "factors" / "factors.json").write_text(json.dumps(rec))
    return bs, run


def test_gene_order_is_the_barred_factors_through_the_builders_loader(db, kept, tmp_path, monkeypatch):
    bs, run = _barred_run(tmp_path, GENES)
    phi, rec = bs.load_factor(run, cq.BARRED_NAME, GENES, {})   # K18's barred_factor -> build_submission.load_factor
    import k18_built_files as k18

    def fake_load(data, run_dir, decision):
        k18._BARRED.clear()
        k18._BARRED.update(phi=phi, genes=list(GENES), spec={}, record=rec)
        return {**rec, "sha256": cq.PHI_M_SHA256}

    monkeypatch.setattr(k18, "load_barred_for_plain", fake_load)
    genes = cq.barred_genes(tmp_path, run, tmp_path)
    assert genes == GENES
    out = tmp_path / "tables"
    written = cq.tables(genes, out, seeds=(0, 1, 2), root=db, kept=kept)
    assert [w[0] for w in written] == ["chipq_Gata4.csv", "chipq_rand_seed0.csv", "chipq_rand_seed1.csv",
                                       "chipq_rand_seed2.csv"]
    for name, sha in written:
        stem = name[:-4]
        T, trec = bs.load_factor(out, stem, GENES, {})   # the builder's loader takes it in the panel's order
        assert trec["sha256"] == sha and len(T) == len(GENES)
        m = json.loads((out / f"{stem}_manifest.json").read_text())
        assert m["files"][stem]["sha256"] == sha
        assert not PLACEBO_NAME.search(json.dumps(m))
        if stem != "chipq_Gata4":
            s = int(stem.rsplit("seed", 1)[1])
            assert m["pool_size"] == 10 and POOL[m["pool_index"]] == cq.placebo_for_seed(POOL, s)[0]
    with pytest.raises(cq.Refused, match="nothing is overwritten"):
        cq.tables(genes, out, seeds=(0,), root=db, kept=kept)
    monkeypatch.setattr(k18, "load_barred_for_plain", lambda *a: {**rec, "sha256": "0" * 64})
    with pytest.raises(cq.Refused):
        cq.barred_genes(tmp_path, run, tmp_path)


# the placebo stand-ins' names (Gata6, named in the code's own exclusion list, is not a drawn placebo)
PLACEBO_NAME = re.compile("|".join(sorted(a for a in PLACEBO if a != "Gata6")))


def _main_tables(db, kept, monkeypatch, out):
    monkeypatch.setattr(cq, "load_kept", lambda *a, **k: kept)
    monkeypatch.setattr(cq, "barred_genes", lambda *a, **k: list(GENES))
    monkeypatch.setattr(cf, "DATA_ROOT", db)
    out.mkdir(parents=True, exist_ok=True)
    monkeypatch.setenv(cq.ENV_TABLES, str(out))
    return cq.main(["tables"])


def test_the_tables_step_writes_into_the_folder_its_variable_names(db, kept, tmp_path, monkeypatch, capsys):
    out = tmp_path / "tables_k27"
    assert _main_tables(db, kept, monkeypatch, out) == 0
    lines = capsys.readouterr().out.splitlines()
    names = ["chipq_Gata4.csv", *(f"chipq_rand_seed{s}.csv" for s in range(10))]
    assert [x.split("  ")[1] for x in lines] == names
    for x in lines:
        sha, name = x.split("  ")
        assert cf.sha256_file(out / "factors" / name) == sha
    with pytest.raises(SystemExit):   # no path from the command line
        cq.main(["tables", "--out", str(tmp_path / "elsewhere")])


def test_the_tables_folder_is_checked(tmp_path, monkeypatch):
    d = tmp_path / "k27"
    d.mkdir()
    assert cq.tables_folder({cq.ENV_TABLES: str(d)}) == d.resolve()
    (tmp_path / "link").symlink_to(d)
    (tmp_path / "file").write_text("x")
    for v, msg in ((None, "not set"), ("", "not set"), ("relative/k27", "absolute"),
                   (str(tmp_path / "link"), "symlink"), (str(tmp_path / "file"), "directory"),
                   (str(tmp_path / "missing"), "directory"),
                   (str(ROOT / "experiments"), "checkout|repository")):
        with pytest.raises(cq.Refused, match=msg):
            cq.tables_folder({} if v is None else {cq.ENV_TABLES: v})
    monkeypatch.delenv(cq.ENV_TABLES, raising=False)
    with pytest.raises(cq.Refused, match="not set"):
        cq.tables_folder()
    assert cq.ENV_TABLES == "VEC_T3_EXT_TABLES_K27"


@pytest.mark.parametrize("case", ["corrupted", "leftover_part"])
def test_no_message_names_a_placebo_antigen(db, kept, tmp_path, monkeypatch, capsys, case):
    """A corrupted placebo file, or a placebo file missing with a leftover partial download beside it: neither the
    tables step of chipq.py nor the targets step of chipq_fetch.py names the antigen in what it prints."""
    p = db / "targets" / "P07.5.tsv"
    raw = p.read_bytes()
    if case == "corrupted":
        p.write_bytes(raw + b"x")
    else:
        p.rename(p.with_name("P07.5.tsv.part"))
    assert _main_tables(db, kept, monkeypatch, tmp_path / "t") == 2
    err = capsys.readouterr()
    assert "REFUSED" in err.err and not PLACEBO_NAME.search(err.err + err.out), err
    rec = json.loads((db / "MANIFEST.json").read_text())["files"]["targets/P07.5.tsv"]
    assert rec["sha256"][:12] in err.err   # the non-identifying key
    code = cf.main(["targets"])
    err = capsys.readouterr()
    assert code == 2 and not PLACEBO_NAME.search(err.err + err.out), err
    # nothing was overwritten or deleted
    if case == "corrupted":
        assert p.read_bytes() == raw + b"x"
    else:
        assert (db / "targets" / "P07.5.tsv.part").read_bytes() == raw


def test_a_placebo_download_refused_mid_way_names_no_antigen(tmp_path, monkeypatch, capsys):
    files = make_db()
    use_net(monkeypatch, FakeNet(files, redirects={
        "https://dbarchive.biosciencedbc.jp/kyushu-u/mm10/target/P03.5.tsv": "https://example.org/P03.5.tsv"}))
    root = tmp_path / "chipatlas"
    monkeypatch.setattr(cf, "DATA_ROOT", root)
    assert cf.main(["lists"]) == 0 and cf.main(["rows"]) == 0
    capsys.readouterr()
    assert cf.main(["targets"]) == 2
    err = capsys.readouterr()
    assert "redirect" in err.err and not PLACEBO_NAME.search(err.err + err.out), err


def test_file_labels():
    assert cf.file_label("lists/chip_atlas_analysis_list.zip") == "lists/chip_atlas_analysis_list.zip"
    assert cf.file_label("targets/P01.5.tsv", {"sha256": "ab" * 32}) == \
        "a target-gene file (manifest entry sha256 abababababab)"
    assert "P01" not in cf.file_label("targets/P01.5.tsv")


def test_tables_refuse_a_pool_under_ten_and_write_nothing(db, kept, tmp_path):
    _rewrite_target(db, "P13", target_file("P13", cf.placebo_kept(db)["P13"], bound_rows(30, 2, 0, seed=5)))
    with pytest.raises(cq.Refused, match="fewer than 10"):
        cq.tables(GENES, tmp_path / "t", seeds=(0,), root=db, kept=kept)
    assert not (tmp_path / "t").exists()


def test_tables_never_open_or_write_the_test_query(db, kept, tmp_path, monkeypatch):
    opened = []
    real = cq.open_target
    monkeypatch.setattr(cq, "open_target", lambda a, *x, **k: (opened.append(a), real(a, *x, **k))[1])
    written = cq.tables(GENES, tmp_path / "t", seeds=range(10), root=db, kept=kept)
    assert "Ctnnb1" not in opened and not any("Ctnnb1" in n for n, _ in written) and len(written) == 11


# ---------------------------------------------------------------------------------------------------------------------
# P1


@pytest.mark.parametrize("n_bound,code", [(16, 0), (14, 3)])
def test_p1_prints_two_integers_and_nothing_else(db, kept, tmp_path, monkeypatch, capsys, n_bound, code):
    if n_bound != 16:
        _rewrite_target(db, "Ctnnb1", target_file("Ctnnb1", QUERY_IDS["Ctnnb1"], bound_rows(n_bound, 4, 3, seed=1)))
    monkeypatch.setattr(cq, "load_kept", lambda *a, **k: kept)
    monkeypatch.setattr(cq, "barred_genes", lambda *a, **k: (print("noise from a loader"), list(GENES))[1])
    monkeypatch.setattr(cf, "DATA_ROOT", db)
    assert cq.main(["p1"]) == code
    out = capsys.readouterr().out
    assert re.fullmatch(r"\d+ \d+\n", out) and out == f"3 {n_bound}\n"


def test_p1_with_no_kept_experiment_opens_nothing(db, tmp_path, monkeypatch, capsys):
    yp, sp = kept_files(tmp_path, {"Ctnnb1": [], "Gata4": KEPT["Gata4"]})
    k = cq.load_kept(yp, sp, check_committed=False)
    monkeypatch.setattr(cq, "open_target", lambda *a, **kw: pytest.fail("opened"))
    assert cq.p1(GENES, db, k) == (0, 0, False)
    monkeypatch.setattr(cq, "load_kept", lambda *a, **kw: k)
    monkeypatch.setattr(cq, "barred_genes", lambda *a, **kw: list(GENES))
    monkeypatch.setattr(cf, "DATA_ROOT", db)
    assert cq.main(["p1"]) == 3 and capsys.readouterr().out == "0 0\n"


def test_p1_refuses_with_exit_2_while_the_kept_list_is_absent(capsys, monkeypatch, tmp_path):
    monkeypatch.setattr(cq, "KEPT_YAML", tmp_path / "chipq_kept.yaml")
    monkeypatch.setattr(cq, "KEPT_SHA", tmp_path / "CHIPQ_KEPT.sha256")
    assert cq.main(["p1"]) == 2
    assert capsys.readouterr().out == ""
