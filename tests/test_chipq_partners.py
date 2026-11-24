"""Tests for experiments/t3_ext/chipq_partners.py (CHIPQ-C, ledger row K28: steps 1 to 3 of the specification's file
02, classes (v) and (vi) of its file 03, and gate G0a's counts) and for the every-row mode of chipq.open_target.

Hand-made files only. The archive is made up here: antigens Ant01.., cell-type classes with made-up names, gene symbols
gs00000.., experiment ids SRX10000.., one made-up score. It is installed through the fetch's own steps over the fake
network of tests/test_chipq.py, so the manifest, the keyword filter's record and the record of unfetched antigens are
the ones the real code writes. No network, no real download, no competition data; the ledger and the register are not
read; no byte built from released or outside data is pinned."""

from __future__ import annotations

import ast
import csv
import importlib.util
import io
import json
import re
import sys
from pathlib import Path
from types import SimpleNamespace
from typing import NamedTuple

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
EXT = ROOT / "experiments" / "t3_ext"
if str(EXT) not in sys.path:
    sys.path.insert(0, str(EXT))
import chipq_partners as cp  # noqa: E402

cq, cf = cp.cq, cp.cf


def _helpers():
    """tests/test_chipq.py's synthetic-archive helpers (FakeNet, use_net, target_file, turl, kept_files), loaded under
    their own name."""
    name = "k28_chipq_helpers"
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, ROOT / "tests" / "test_chipq.py")
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


H = _helpers()

# ---------------------------------------------------------------------------------------------------------------------
# The synthetic archive

CLEAN = "wild type untreated cells"
KA, LA, MU, NU, XI = "Kappa tissue", "Lambda tissue", "Mu tissue", "Nu tissue", "Xi tissue"  # made-up classes
SCORE = 7777.5  # the one made-up score (column j holds SCORE + j at its bound genes, 0 elsewhere)
N_ROWS = 260
SYMBOLS = [f"gs{k:05d}" for k in range(N_ROWS)][::-1]  # made-up gene symbols, in a file order that is not sorted
# every made-up name of the archive: none may reach a message, a printed line, a log or a written file
MADE_UP = re.compile(r"Ant\d\d|Bare\d|Gata6|Mab21l2|Kappa|Lambda|\bMu tissue|\bNu tissue|\bXi tissue|Unclassified"
                     r"|gs\d{5}|SRX\d+|7777")


class X(NamedTuple):
    """One row of the experiment list, and what its column of the antigen's target-gene file holds."""

    id: str
    antigen: str
    cls: str = KA
    text: str = CLEAN  # the title, and the metadata after "source="
    bound: int = 250  # genes with a score above 0 in its column
    genome: str = "mm10"
    column: bool = True  # a column in its antigen's target-gene file
    kept: bool = True  # a query's experiment only: on the readers' kept list
    cell_type: str = "cells"


def experiments_zip(exps: list[X]) -> bytes:
    out = io.StringIO()
    w = csv.writer(out, lineterminator="\r\n")
    w.writerow(cf.EXPERIMENT_COLUMNS)
    for x in exps:
        w.writerow([x.id, x.genome, "TFs and others", x.antigen, x.cls, x.cell_type, "NA", "1,2,3,4", x.text,
                    "source=" + x.text, *(["link"] * 7)])
    return H._zip("chip_atlas_experiment_list.csv", out.getvalue())


def analysis_zip(with_file: list[str], without_file: list[str]) -> bytes:
    head = [*cf.ANALYSIS_COLUMNS, "Target Genes (TSS ±1k)", "Target Genes (TSS ±5k)", "Target Genes (TSS ±10k)",
            "Colocalization (TSV)", "Colocalization (GML)"]
    out = io.StringIO()
    w = csv.writer(out, lineterminator="\r\n")
    w.writerow(head)
    for a in with_file:
        w.writerow([a, "", "+", "mm10", "u1", H.turl(a), "u10", "x", "y"])
    for a in without_file:  # held for mm10, no target-gene file recorded
        w.writerow([a, "", "-", "mm10", "", "", "", "x", "y"])
    w.writerow([with_file[0], "", "+", "hg38", "u1", "u5", "u10", "x", "y"])  # another genome: never read
    return H._zip("chip_atlas_analysis_list.csv", out.getvalue())


def target_bytes(antigen: str, cols: list[X], **kw) -> bytes:
    """The antigen's target-gene file: N_ROWS rows, column j bound (SCORE + j) at its first x.bound rows. Every column
    head names tests/test_chipq.py's class "Unclassified": an experiment's class is never taken from this file."""
    rows = {g: [SCORE + j if k < x.bound else 0.0 for j, x in enumerate(cols)] for k, g in enumerate(SYMBOLS)}
    return H.target_file(antigen, [x.id for x in cols], rows, **kw)


def columns(exps: list[X], antigen: str) -> list[X]:
    return [x for x in exps if x.antigen == antigen and x.genome == "mm10" and x.column]


def build_db(tmp_path: Path, monkeypatch, exps: list[X], name: str = "db", no_file=(), unserved=(),
             targets: bool = True) -> Path:
    """A pinned copy of the synthetic archive under tmp_path/name, through chipq_fetch's lists, rows and targets
    steps. `no_file`: antigens the analysis list records without a target-gene file. `unserved`: antigens whose file
    the archive does not serve (HTTP 404: recorded by the targets step)."""
    with_file = sorted({x.antigen for x in exps if x.genome == "mm10"} - set(no_file))
    files = {cf.LISTS[cf.ANALYSIS]: analysis_zip(with_file, sorted(no_file)),
             cf.LISTS[cf.EXPERIMENTS]: experiments_zip(exps)}
    for a in with_file:
        files[H.turl(a)] = target_bytes(a, columns(exps, a))
    H.use_net(monkeypatch, H.FakeNet(files, status={H.turl(a): 404 for a in unserved}))
    monkeypatch.setattr(cf, "_sleep", lambda s: None)
    root = tmp_path / name
    cf.step_lists(root)
    cf.step_rows(root)
    if targets:
        cf.step_targets(root)
    return root


def kept_ids(exps: list[X]) -> dict[str, list[str]]:
    return {q: sorted(x.id for x in exps if x.antigen == q and x.kept and x.genome == "mm10") for q in cp.QUERIES}


def kept_for(tmp_path: Path, exps: list[X] | None = None, name: str = "kept", ids: dict | None = None):
    """The pinned kept list (chipq.load_kept on a hand-made chipq_kept.yaml and its sha256 file)."""
    d = tmp_path / name
    d.mkdir(exist_ok=True)
    yp, sp = H.kept_files(d, kept_ids(exps) if ids is None else ids)
    return cq.load_kept(yp, sp, check_committed=False)


def base_exps() -> list[X]:
    """Twelve antigens with three usable experiments in each of two classes; the two queries, one class each; and one
    antigen for every way an experiment is not usable or an antigen is not a candidate."""
    n = iter(range(10000, 100000))

    def x(antigen, cls=KA, count=1, **kw):
        return [X(f"SRX{next(n)}", antigen, cls, **kw) for _ in range(count)]

    e = x("Ctnnb1", KA, 3) + x("Ctnnb1", LA, 1) + x("Ctnnb1", KA, 1, kept=False)  # class Kappa, 3 usable there
    e += x("Gata4", LA, 4) + x("Gata4", KA, 3)  # class Lambda, 4 usable there (and 3 in Kappa)
    for k in range(1, 13):
        e += x(f"Ant{k:02d}", KA, 3) + x(f"Ant{k:02d}", LA, 3)  # the 12 candidates of either class
    e += x("Ant01", KA, 1, text="GSK-3 blocked")  # one matching row among Kappa's 65: under 5%
    e += x("Gata6", KA, 3) + x("Gata6", LA, 3)  # never a candidate (named)
    e += x("Mab21l2", KA, 3)  # never a candidate (named; and it has no kept list)
    e += x("Ant13", KA, 2) + x("Ant13", KA, 1, text="hepatocellular carcinoma line")  # class (v) leaves two
    e += x("Ant14", KA, 2) + x("Ant14", KA, 1, bound=199)  # 199 genes: two usable
    e += x("Ant15", KA, 2) + x("Ant15", KA, 1, column=False)  # no column: two usable
    e += x("Ant16", MU, 3) + x("Ant16", MU, 1, text="Wnt3a conditioned")  # its class is excluded (1 row of 4)
    e += x("Ant17", "Unclassified", 3)  # its class is excluded by name
    e += x("Ant18", KA, 3)  # the archive serves no file for it (build_db's `unserved`)
    e += x("Ant19", KA, 1) + x("Ant19", KA, 2, text="Ant19 shRNA knockdown")  # one kept: the fetch takes no file
    e += x("Ant20", NU, 19) + x("Ant20", NU, 1, text="CHIR99021 3uM")  # 1 row of 20, exactly 5%: excluded
    e += x("Ant21", XI, 20) + x("Ant21", XI, 1, text="LiCl 20mM")  # 1 row of 23: not excluded, whatever follows
    e += x("Ant21", XI, 2, text="LPS stimulated")  # entry 16 of the filter: not one of class (vi)'s entries
    e += x("Ant21", XI, 5, text="Wnt3a conditioned", genome="mm9")  # another genome: not counted
    e += x("Bare1", XI, 5, text="ES cells in 2i/LIF")  # an antigen without a target-gene file: not counted
    return e


BASE_COUNTS = {"classes_excluded": 3,
               "queries": {"Ctnnb1": {"usable_in_class": 3, "candidates": 12, "a1": True, "a2": True},
                           "Gata4": {"usable_in_class": 4, "candidates": 12, "a1": True, "a2": True}},
               "passes": True}
BASE_OPENED = sorted(["Ctnnb1", "Gata4", *(f"Ant{k:02d}" for k in range(1, 13)), "Ant14", "Ant15"])
TWELVE = [f"Ant{k:02d}" for k in range(1, 13)]


def build_base(tmp_path: Path, monkeypatch, exps: list[X] | None = None, **kw) -> SimpleNamespace:
    exps = base_exps() if exps is None else exps
    root = build_db(tmp_path, monkeypatch, exps, no_file=("Bare1",), unserved=("Ant18",), **kw)
    return SimpleNamespace(root=root, exps=exps, kept=kept_for(tmp_path, exps), tmp=tmp_path)


def changed(exps: list[X], pick, **fields) -> list[X]:
    """The experiments with `fields` replaced on those `pick` selects."""
    return [x._replace(**fields) if pick(x) else x for x in exps]


@pytest.fixture(autouse=True)
def no_pause(monkeypatch):
    monkeypatch.setattr(cf, "_sleep", lambda s: None)


@pytest.fixture(scope="module")
def base(tmp_path_factory):
    """The base archive, built once and only read by the tests that take it."""
    mp = pytest.MonkeyPatch()
    try:
        return build_base(tmp_path_factory.mktemp("base"), mp)
    finally:
        mp.undo()


@pytest.fixture
def opened(monkeypatch):
    """Every call of the one reader: (antigen, the ids asked for, the mode)."""
    calls = []
    real = cq.open_target

    def spy(antigen, genes, root=None, kept=None, ids=None, man=None, **kw):
        calls.append((antigen, None if ids is None else list(ids), bool(kw.get("every_row")), genes))
        return real(antigen, genes, root, kept, ids, man, **kw)

    monkeypatch.setattr(cq, "open_target", spy)
    return calls


# ---------------------------------------------------------------------------------------------------------------------
# The constants


def test_the_constants_are_the_specifications():
    assert (cp.MIN_BOUND_GENES, cp.QUERY_MIN_USABLE, cp.CANDIDATE_MIN_USABLE, cp.MIN_CANDIDATES) == (200, 2, 3, 11)
    assert cp.QUERIES == ("Ctnnb1", "Gata4") == (cq.TEST_QUERY, cq.PRACTICE_QUERY)
    assert cp.CANDIDATE_EXCLUDE == ("Ctnnb1", "Gata4", "Gata6", "Mab21l2") == (*cq.POOL_EXCLUDE, "Mab21l2")
    assert set(cf.QUERY_ANTIGENS) <= set(cp.CANDIDATE_EXCLUDE)
    assert cp.CLASS_V_WORDS == ("carcinoma", "sarcoma", "lymphoma", "leukemia", "leukaemia", "myeloma", "melanoma",
                                "glioma", "blastoma", "adenoma", "tumor", "tumour", "cancer", "neoplas", "malignan",
                                "metasta", "transform", "immortal")
    assert cp.CLASS_VI_ENTRIES == (15, 17, 18, 19, 20) and cp.CLASS_VI_SHARE == (5, 100)
    assert cp.CLASS_VI_NAMED == ("Unclassified", "No description", "Others")
    assert cp.CLASS_FIELD == "Cell type class" == cf.EXPERIMENT_COLUMNS[4] and cp._CLASS_COL == 4
    assert 4 in cf.SCREENED_COLS  # the class is one of the seven screened fields
    assert cp.readings()["MIN_BOUND_GENES"] == 200 and cp.readings()["CLASS_VI_ENTRIES"] == [15, 17, 18, 19, 20]


def test_class_vi_entries_are_the_filters_class_ii_entries_other_than_the_five_words():
    words = [w for _, w, _ in cf.DROP_PATTERNS]
    assert [words[i] for i in cp.CLASS_VI_ENTRIES] == ["Wnt-stimulated / Wnt3a / Wnt-ligand", "GSK3-inhibitor", "CHIR",
                                                       "LiCl", "a medium naming a pathway drug (2i)"]
    ii = [i for i, (c, _, _) in enumerate(cf.DROP_PATTERNS) if c == "ii"]
    assert ii == list(range(11, 21))
    assert sorted(set(ii) - set(cp.CLASS_VI_ENTRIES)) == [11, 12, 13, 14, 16]
    assert [words[i] for i in (11, 12, 13, 14, 16)] == list(cp.CLASS_VI_OTHER_WORDS)
    pats = cp.class_vi_patterns()
    assert [p.pattern for p in pats] == [cf.DROP_PATTERNS[i][2] for i in cp.CLASS_VI_ENTRIES]
    assert all(p.flags & re.IGNORECASE for p in pats)


@pytest.mark.parametrize("change", ["inserted", "renamed", "removed"])
def test_a_changed_keyword_filter_is_never_read_by_position(monkeypatch, change):
    pats = list(cf.DROP_PATTERNS)
    if change == "inserted":
        pats.insert(3, ("i", "another word", r"another"))
    elif change == "renamed":
        pats[16] = ("ii", "stimulation", pats[16][2])
    else:
        del pats[18]
    monkeypatch.setattr(cf, "DROP_PATTERNS", tuple(pats))
    monkeypatch.setattr(cf, "_DROP", tuple((c, w, re.compile(p, re.IGNORECASE)) for c, w, p in pats))
    with pytest.raises(cp.Refused, match="keyword filter"):
        cp.class_vi_patterns()


# ---------------------------------------------------------------------------------------------------------------------
# Classes (v) and (vi), in code


@pytest.mark.parametrize("word", cp.CLASS_V_WORDS)
def test_every_class_v_word_drops_in_any_case_and_any_field(word):
    clean = ["TFs and others", "Ant01", KA, "cells", "NA", CLEAN, "source=" + CLEAN]
    assert not cp.dropped_v(clean)
    for form in (word, word.upper(), word.capitalize(), f"x{word}y"):
        for k in range(len(clean)):
            fields = list(clean)
            fields[k] = f"a {form} line"
            assert cp.dropped_v(fields), (form, k)


@pytest.mark.parametrize("text", ["immortalized fibroblasts", "Transformed B cells", "neoplastic", "metastasis",
                                  "MALIGNANT", "neuroblastoma", "adenocarcinoma", "Leukaemia", "HeLa cancer cells"])
def test_class_v_takes_stems_and_inflections(text):
    assert cp.dropped_v(["NA", text])


@pytest.mark.parametrize("text", [CLEAN, "ES cells", "embryonic fibroblasts", "liver", "wild type", "tumo", "cance r",
                                  "trans form", "", "bone marrow", "thymus"])
def test_plain_words_are_not_class_v(text):
    assert not cp.dropped_v(["TFs and others", "Ant01", KA, text, text])


def test_class_vi_counts_its_five_entries_and_no_other_word():
    ex = H.DROP_EXAMPLES
    counted = [w for i, (_, w, _) in enumerate(cf.DROP_PATTERNS) if i in cp.CLASS_VI_ENTRIES]
    for _, w, _ in cf.DROP_PATTERNS:
        assert cp.matches_vi(["NA", ex[w]]) == (w in counted), w
    assert cp.matches_vi([CLEAN, "cultured in 2i medium"]) and cp.matches_vi(["wnt3A", CLEAN])  # any field, any case
    assert not cp.matches_vi([CLEAN, "LPS stimulated", "activated", "stabilized allele", "constitutively active"])


def test_the_share_is_exact_at_five_percent():
    assert cp.share_excludes(1, 20) and cp.share_excludes(5, 100) and cp.share_excludes(3, 4)
    assert not cp.share_excludes(1, 21) and not cp.share_excludes(4, 100) and not cp.share_excludes(0, 7)
    assert not cp.share_excludes(0, 0)
    assert cp.share_excludes(50, 1000) and not cp.share_excludes(49, 1000)


# ---------------------------------------------------------------------------------------------------------------------
# The reader's every-row mode (chipq.open_target, row K28's one new mode)


def small_exps() -> list[X]:
    e = [X("SRX20000", "Ctnnb1", KA, bound=250), X("SRX20001", "Ctnnb1", KA, bound=210),
         X("SRX20002", "Ctnnb1", LA, bound=205), X("SRX20003", "Ctnnb1", KA, bound=260, kept=False),
         X("SRX20004", "Gata4", LA, bound=230), X("SRX20005", "Gata4", LA, bound=1)]
    return e + [X(f"SRX2001{k}", "Ant01", KA, bound=200 + k) for k in range(3)]


@pytest.fixture
def small(tmp_path, monkeypatch):
    exps = small_exps()
    root = build_db(tmp_path, monkeypatch, exps)
    return SimpleNamespace(root=root, exps=exps, kept=kept_for(tmp_path, exps), tmp=tmp_path)


def test_every_row_mode_gives_every_row_for_kept_columns_in_file_order(small):
    E, ids, symbols = cq.open_target("Ctnnb1", None, small.root, small.kept, every_row=True)
    assert ids == ["SRX20000", "SRX20001", "SRX20002"] == list(small.kept["Ctnnb1"])
    assert symbols == SYMBOLS and E.shape == (N_ROWS, 3) and E.dtype == np.float64
    assert (E > 0).sum(axis=0).tolist() == [250, 210, 205]
    assert E[0].tolist() == [SCORE, SCORE + 1, SCORE + 2] and not E[N_ROWS - 1].any()
    # the panel mode's scores, row for row, for any gene list
    genes = SYMBOLS[190:230]
    Ep, idp = cq.open_target("Ctnnb1", genes, small.root, small.kept)
    assert idp == ids and np.array_equal(Ep, E[190:230])
    # asked for by name, in the order asked
    E2, ids2, sym2 = cq.open_target("Ctnnb1", None, small.root, small.kept, ["SRX20002", "SRX20000"], every_row=True)
    assert ids2 == ["SRX20002", "SRX20000"] and sym2 == SYMBOLS and np.array_equal(E2, E[:, [2, 0]])
    # an antigen of the keyword filter's record: its kept ids, no kept list given
    E3, ids3, _ = cq.open_target("Ant01", None, small.root, None, every_row=True)
    assert ids3 == ["SRX20010", "SRX20011", "SRX20012"] and (E3 > 0).sum(axis=0).tolist() == [200, 201, 202]


def test_every_row_mode_keeps_the_readers_refusals(small):
    with pytest.raises(cq.Refused, match="not in the pinned kept list"):
        cq.open_target("Ctnnb1", None, small.root, small.kept, ["SRX20003"], every_row=True)  # not kept
    with pytest.raises(cq.Refused, match="not in the pinned kept list"):
        cq.open_target("Ant01", None, small.root, None, ["SRX20000"], every_row=True)
    with pytest.raises(cq.Refused, match="no gene list"):
        cq.open_target("Ctnnb1", SYMBOLS, small.root, small.kept, every_row=True)
    with pytest.raises(cq.Refused, match="placebo_kept.json only"):
        cq.open_target("Ant01", None, small.root, small.kept, every_row=True)
    with pytest.raises(cq.Refused, match="load_kept"):
        cq.open_target("Ctnnb1", None, small.root, {"Ctnnb1": ["SRX20000"]}, every_row=True)
    with pytest.raises(cq.Refused, match="asked for twice"):
        cq.open_target("Ctnnb1", None, small.root, small.kept, ["SRX20000", "SRX20000"], every_row=True)


def test_every_row_mode_leaves_out_an_experiment_without_a_column(tmp_path, monkeypatch):
    """CHIPQ-C's step 1: a column is a condition of a usable experiment, for a query too. The panel mode still refuses
    a query's kept id without one."""
    exps = changed(small_exps(), lambda x: x.id == "SRX20001", column=False)
    root = build_db(tmp_path, monkeypatch, exps)
    kept = kept_for(tmp_path, exps)
    assert list(kept["Ctnnb1"]) == ["SRX20000", "SRX20001", "SRX20002"]
    E, ids, symbols = cq.open_target("Ctnnb1", None, root, kept, every_row=True)
    assert ids == ["SRX20000", "SRX20002"] and E.shape == (N_ROWS, 2) and len(symbols) == N_ROWS
    E, ids, symbols = cq.open_target("Ctnnb1", None, root, kept, ["SRX20001"], every_row=True)
    assert ids == [] and E.shape == (N_ROWS, 0) and symbols == SYMBOLS
    with pytest.raises(cq.Refused, match="no column"):
        cq.open_target("Ctnnb1", SYMBOLS[:10], root, kept)


def _poisoned(raw: bytes, column: str, value: str) -> bytes:
    """The target file with `value` in every row of the column whose head begins with `column`."""
    lines = raw.decode().splitlines()
    head = lines[0].split("\t")
    k = next(i for i, h in enumerate(head) if h.split("|")[0] == column)
    out = [lines[0]] + ["\t".join(f[:k] + [value] + f[k + 1:]) for f in (ln.split("\t") for ln in lines[1:])]
    return ("\n".join(out) + "\n").encode()


def test_every_row_mode_never_reads_average_string_or_a_column_not_kept(small):
    ref, ids, symbols = cq.open_target("Ctnnb1", None, small.root, small.kept, every_row=True)
    cols = columns(small.exps, "Ctnnb1")
    for avg, string in (("NaN", "garbage"), ("-1e300", "xéy"), (lambda g: "oops", lambda g: "Ant01")):
        H._rewrite_target(small.root, "Ctnnb1", target_bytes("Ctnnb1", cols, avg=avg, string=string))
        E, i2, s2 = cq.open_target("Ctnnb1", None, small.root, small.kept, every_row=True)
        assert np.array_equal(E, ref) and i2 == ids and s2 == symbols
    for value in ("nope", "-5", "inf", ""):  # the experiment the readers did not keep
        H._rewrite_target(small.root, "Ctnnb1", _poisoned(target_bytes("Ctnnb1", cols), "SRX20003", value))
        assert np.array_equal(cq.open_target("Ctnnb1", None, small.root, small.kept, every_row=True)[0], ref)
    # a kept column that is not a finite, non-negative number is refused, and the message names nothing
    for value in ("oops", "-5", "inf", "nan", ""):
        H._rewrite_target(small.root, "Ctnnb1", _poisoned(target_bytes("Ctnnb1", cols), "SRX20001", value))
        with pytest.raises(cq.Refused) as e:
            cq.open_target("Ctnnb1", None, small.root, small.kept, every_row=True)
        assert not MADE_UP.search(str(e.value)) and "oops" not in str(e.value), str(e.value)
        # the columns asked for are the only ones parsed
        E, got, _ = cq.open_target("Ctnnb1", None, small.root, small.kept, ["SRX20000", "SRX20002"], every_row=True)
        assert got == ["SRX20000", "SRX20002"] and np.array_equal(E, ref[:, [0, 2]])


@pytest.mark.parametrize("case", ["two rows", "no symbol", "short row", "long row", "blank lines and CRLF"])
def test_every_row_mode_on_odd_rows(small, case):
    cols = columns(small.exps, "Ctnnb1")
    lines = target_bytes("Ctnnb1", cols).decode().splitlines()
    ref = cq.open_target("Ctnnb1", None, small.root, small.kept, every_row=True)[0]
    row = lines[7].split("\t")
    fresh = ["gs99999", *row[1:]]  # a further gene, so that only the fault named is in the file
    if case == "two rows":
        text, msg = "\n".join([*lines, lines[7]]) + "\n", "two rows"
    elif case == "no symbol":
        text, msg = "\n".join([*lines, "\t".join(["", *row[1:]])]) + "\n", "no gene symbol"
    elif case == "short row":
        text, msg = "\n".join([*lines, "\t".join(fresh[:-1])]) + "\n", "number of fields"
    elif case == "long row":
        text, msg = "\n".join([*lines[:9], "\t".join([*fresh, "0"]), *lines[9:]]) + "\n", "number of fields"
    else:
        text = "\r\n".join([lines[0], "", *lines[1:], "", ""]) + "\r\n"
    H._rewrite_target(small.root, "Ctnnb1", text.encode())
    if case == "blank lines and CRLF":
        E, _, symbols = cq.open_target("Ctnnb1", None, small.root, small.kept, every_row=True)
        assert np.array_equal(E, ref) and symbols == SYMBOLS
        return
    with pytest.raises(cq.Refused, match=msg) as e:
        cq.open_target("Ctnnb1", None, small.root, small.kept, every_row=True)
    assert not MADE_UP.search(str(e.value)), str(e.value)
    if case == "two rows":  # the panel mode's own rule stands: a second row of a gene off the panel is not its matter
        other = [g for g in SYMBOLS[:20] if g != row[0]]
        assert cq.open_target("Ctnnb1", other, small.root, small.kept)[0].shape == (len(other), 3)
    # the fault alone: the same file with a well-formed further row is read
    H._rewrite_target(small.root, "Ctnnb1", ("\n".join([*lines, "\t".join(fresh)]) + "\n").encode())
    E, _, symbols = cq.open_target("Ctnnb1", None, small.root, small.kept, every_row=True)
    assert symbols == [*SYMBOLS, "gs99999"] and np.array_equal(E[:-1], ref) and np.array_equal(E[-1], E[6])


# ---------------------------------------------------------------------------------------------------------------------
# The screen: classes (v) and (vi) before any count


def test_the_screen_excludes_classes_by_share_and_by_name(base, opened):
    sc = cp.screen(base.root, base.kept)
    # Mu: 1 row of 4; Nu: 1 of 20, exactly 5%; Unclassified: by name. Kappa (1 of 65) and Lambda stay. Xi stays: 1 of
    # its 23 rows counts; its rows of entry 16, of another genome and of an antigen without a file do not
    assert sc.excluded == frozenset({MU, NU, "Unclassified"})
    assert opened == []  # the screen reads the two lists and the records only
    el = sc.eligible
    assert set(el) == {x.antigen for x in base.exps} - {"Mab21l2", "Bare1"}
    assert [c for _, c in el["Ctnnb1"]] == [KA, KA, KA, LA] and [c for _, c in el["Gata4"]] == [LA] * 4 + [KA] * 3
    assert [i for i, _ in el["Ctnnb1"]] == list(base.kept["Ctnnb1"])
    assert len(el["Ant01"]) == 6 and len(el["Ant13"]) == 2 and len(el["Ant14"]) == 3 and len(el["Ant15"]) == 3
    assert el["Ant16"] == [] and el["Ant17"] == [] and el["Ant20"] == []  # their classes are excluded
    assert len(el["Ant18"]) == 3 and len(el["Ant19"]) == 1 and len(el["Ant21"]) == 20
    assert sc.no_file == frozenset({"Ant18", "Ant19"})
    assert all(i == sorted(i) for i in ([j for j, _ in v] for v in el.values()))


@pytest.mark.parametrize("wrong", ["entry 16 counted", "another genome counted", "an antigen without a file counted"])
def test_the_share_counts_only_what_class_vi_names(base, wrong):
    """Each wrong reading of class (vi) would exclude Xi in the base archive: the test of the test above."""
    exps = base.exps
    in_xi = [x for x in exps if x.cls == XI]
    with_file = [x for x in in_xi if x.genome == "mm10" and x.antigen != "Bare1"]
    assert len(with_file) == 23 and sum(1 for x in with_file if cp.matches_vi([x.text])) == 1
    if wrong == "entry 16 counted":
        hits, rows = sum(1 for x in with_file if not cf.keep_row([x.text])), 23
    elif wrong == "another genome counted":
        pool = [x for x in in_xi if x.antigen != "Bare1"]
        hits, rows = sum(1 for x in pool if cp.matches_vi([x.text])), len(pool)
    else:
        pool = [x for x in in_xi if x.genome == "mm10"]
        hits, rows = sum(1 for x in pool if cp.matches_vi([x.text])), len(pool)
    assert cp.share_excludes(hits, rows) and not cp.share_excludes(1, 23)


def test_a_query_row_is_screened_by_classes_v_and_vi_too(tmp_path, monkeypatch):
    exps = changed(base_exps(), lambda x: x.antigen == "Ctnnb1" and x.cls == LA, text="immortalized line")
    exps = changed(exps, lambda x: x.antigen == "Gata4" and x.cls == KA, cls=MU)
    b = build_base(tmp_path, monkeypatch, exps)
    sc = cp.screen(b.root, b.kept)
    assert [c for _, c in sc.eligible["Ctnnb1"]] == [KA] * 3  # class (v) dropped its Lambda experiment
    assert [c for _, c in sc.eligible["Gata4"]] == [LA] * 4  # class (vi) excludes the class of the other three
    assert len(b.kept["Ctnnb1"]) == 4 and len(b.kept["Gata4"]) == 7  # all of them on the readers' kept list


def test_the_screen_refuses_what_the_pinned_files_do_not_bear_out(base, tmp_path, monkeypatch):
    k = kept_ids(base.exps)
    ant01 = cf.placebo_kept(base.root)["Ant01"]
    cases = {
        "no mm10 row": ({**k, "Ctnnb1": sorted([*k["Ctnnb1"], "SRX99999"])}, "no mm10 row"),
        "another antigen": ({"Ctnnb1": sorted([*k["Ctnnb1"], k["Gata4"][0]]), "Gata4": k["Gata4"][1:]},
                            "names another antigen"),
        "kept twice": ({**k, "Ctnnb1": sorted([*k["Ctnnb1"], ant01[0]])}, "kept twice"),
        "no query": ({"Gata4": k["Gata4"]}, "no entry for a query"),
    }
    for n, (name, (ids, msg)) in enumerate(cases.items()):
        kept = kept_for(tmp_path, name=f"kept{n}", ids=ids)
        with pytest.raises(cp.Refused, match=msg) as e:
            cp.screen(base.root, kept)
        assert not MADE_UP.search(str(e.value)), (name, str(e.value))
    with pytest.raises(cp.Refused, match="load_kept"):
        cp.screen(base.root, {q: tuple(v) for q, v in k.items()})
    # another genome's row does not stand in for an mm10 row
    mm9 = next(x.id for x in base.exps if x.genome == "mm9")
    with pytest.raises(cp.Refused, match="1 kept experiment"):
        cp.screen(base.root, kept_for(tmp_path, name="kept_mm9", ids={**k, "Gata4": sorted([*k["Gata4"], mm9])}))
    # an mm10 id twice in the experiment list
    real = cf.experiment_rows
    monkeypatch.setattr(cf, "experiment_rows", lambda raw: [*real(raw), next(iter(real(raw)))])
    with pytest.raises(cp.Refused, match="appears twice"):
        cp.screen(base.root, base.kept)


def test_the_screen_refuses_before_a_complete_targets_pass(tmp_path, monkeypatch):
    exps = small_exps()
    root = build_db(tmp_path, monkeypatch, exps, targets=False)
    with pytest.raises(cp.Refused, match="incomplete"):
        cp.screen(root, kept_for(tmp_path, exps))
    with pytest.raises(cp.Refused, match="incomplete"):
        cp.g0a(root, kept_for(tmp_path, exps))


# Class (vi)'s two counts (reading R1, the ruling of 2026-11-24 on OPEN-4): over the rows of the antigens the analysis
# list records with a target-gene file, and over the rows of the antigens whose target-gene file is in the manifest


def two_counts(cls: str) -> list[X]:
    """Rows of one class that class (vi)'s two counts read apart, beside an antigen the fetch holds no file for.
    Mu: 1 matching row of 20 among antigens with a fetched file (Ant02), 1 of 23 once the 3 rows of Ant03 are added,
    whose file the archive does not serve (build_db's `unserved`): only the manifest count excludes Mu. Nu: no
    matching row among the 20 of Ant04, whose file is fetched, 2 of 23 once the rows of Ant05 are added, which keeps
    one experiment (the fetch takes no file for it): only the analysis-list count excludes Nu."""
    n = iter(range(40000 if cls == MU else 50000, 60000))

    def x(antigen, count, **kw):
        return [X(f"SRX{next(n)}", antigen, cls, **kw) for _ in range(count)]

    if cls == MU:
        return x("Ant02", 19) + x("Ant02", 1, text="CHIR99021 3uM") + x("Ant03", 3)
    return x("Ant04", 20) + x("Ant05", 1) + x("Ant05", 2, text="LiCl 20mM")


@pytest.mark.parametrize("cls,alone", [(MU, "manifest"), (NU, "analysis list")])
def test_a_class_that_one_count_alone_excludes_is_excluded(tmp_path, monkeypatch, cls, alone):
    exps = small_exps() + two_counts(cls)
    root = build_db(tmp_path, monkeypatch, exps, unserved=("Ant03",) if cls == MU else ())
    kept = kept_for(tmp_path, exps)
    sc = cp.screen(root, kept)
    # the archive is as the case says: the analysis list records a target-gene file for both antigens of the class;
    # the manifest holds the first one's only
    fetched, unfetched = ("Ant02", "Ant03") if cls == MU else ("Ant04", "Ant05")
    man = cf.load_manifest(root)
    urls = {rec.get("url") for rec in man["files"].values()}
    assert H.turl(fetched) in urls and H.turl(unfetched) not in urls and sc.no_file == frozenset({unfetched})
    assert {fetched, unfetched} <= set(cf._lists(root, man)[0])
    # the two shares, by hand: one count alone reaches 5%
    rows = [x for x in exps if x.cls == cls]
    rows_m = [x for x in rows if x.antigen == fetched]
    hits, hits_m = (sum(1 for x in r if cp.matches_vi([x.text])) for r in (rows, rows_m))
    assert (len(rows), len(rows_m)) == (23, 20)
    by_list, by_manifest = cp.share_excludes(hits, len(rows)), cp.share_excludes(hits_m, len(rows_m))
    assert (by_list, by_manifest) == ((False, True) if alone == "manifest" else (True, False))
    # the screen excludes the class whichever count excludes it, with every experiment in it
    assert sc.excluded == frozenset({cls})
    assert sc.eligible[fetched] == [] and sc.eligible[unfetched] == []
    assert cp.g0a(root, kept)["classes_excluded"] == 1


def test_the_number_of_classes_excluded_counts_a_class_either_count_excludes(tmp_path, monkeypatch):
    exps = small_exps() + two_counts(MU) + two_counts(NU)
    root = build_db(tmp_path, monkeypatch, exps, unserved=("Ant03",))
    kept = kept_for(tmp_path, exps)
    assert cp.screen(root, kept).excluded == frozenset({MU, NU})
    assert cp.g0a(root, kept)["classes_excluded"] == 2  # the number G0a prints
    # with the named classes too, where the list has them
    exps += [X("SRX60000", "Ant06", "Others"), X("SRX60001", "Ant06", "Others")]
    root = build_db(tmp_path, monkeypatch, exps, name="db2", unserved=("Ant03",))
    kept = kept_for(tmp_path, exps, name="kept2")
    assert cp.screen(root, kept).excluded == frozenset({MU, NU, "Others"})
    assert cp.g0a(root, kept)["classes_excluded"] == 3
    assert cp.readings()["CLASS_VI_ROWS"] == ["analysis list", "manifest"]
    assert cp.readings()["CLASS_VI_EXCLUDED_BY"] == "either count"


# ---------------------------------------------------------------------------------------------------------------------
# Step 1: usable experiments


def test_usable_experiments_need_a_column_and_two_hundred_bound_genes(base, opened):
    sc = cp.screen(base.root, base.kept)
    by = {a: [x for x in base.exps if x.antigen == a and x.genome == "mm10"] for a in ("Ant14", "Ant15", "Ant01")}
    u14, u15 = cp.usable(sc, "Ant14"), cp.usable(sc, "Ant15")
    assert sorted(u14) == sorted(x.id for x in by["Ant14"] if x.bound >= 200) and len(u14) == 2
    assert sorted(u15) == sorted(x.id for x in by["Ant15"] if x.column) and len(u15) == 2
    u01 = cp.usable(sc, "Ant01")
    assert u01 == {x.id: x.cls for x in by["Ant01"] if x.text == CLEAN} and sorted(u01.values()) == [KA] * 3 + [LA] * 3
    # each file opened once, in the every-row mode, for the kept experiments left after classes (v) and (vi) only
    assert cp.usable(sc, "Ant01") is u01 and [c[0] for c in opened] == ["Ant14", "Ant15", "Ant01"]
    assert all(mode and genes is None for _, _, mode, genes in opened)
    assert opened[0][1] == [i for i, _ in sc.eligible["Ant14"]] and len(opened[1][1]) == 3
    # no file is opened where nothing can be usable
    del opened[:]
    assert cp.usable(sc, "Ant18") == {} and cp.usable(sc, "Ant19") == {}  # the fetch holds no file
    assert cp.usable(sc, "Ant16") == {} and cp.usable(sc, "Ant17") == {} and cp.usable(sc, "Ant20") == {}
    assert opened == []
    assert len(cp.usable(sc, "Ant13")) == 2 and len(opened[0][1]) == 2  # the class (v) row's column is not asked for
    for a in ("Mab21l2", "Bare1", "Nobody"):
        with pytest.raises(cp.Refused, match="no kept list") as e:
            cp.usable(sc, a)
        assert not MADE_UP.search(str(e.value))


@pytest.mark.parametrize("bound,n", [(200, 3), (199, 0), (201, 3), (0, 0)])
def test_the_two_hundred_gene_bar_is_exact(tmp_path, monkeypatch, bound, n):
    exps = changed(small_exps(), lambda x: x.antigen == "Ant01", bound=bound)
    root = build_db(tmp_path, monkeypatch, exps)
    sc = cp.screen(root, kept_for(tmp_path, exps))
    assert len(cp.usable(sc, "Ant01")) == n
    assert sorted(cp.usable(sc, "Ctnnb1").values()) == [KA, KA, LA]  # 250, 210 and 205 genes
    assert list(cp.usable(sc, "Gata4").values()) == [LA]  # 230 genes; the other has 1


def test_an_antigen_missing_from_the_manifest_without_a_record_is_refused(tmp_path, monkeypatch):
    """Reading R6: only an antigen the targets pass recorded, or one with fewer than two kept experiments, is without a
    file; any other antigen without one is an incomplete fetch, refused by the reader, naming nothing."""
    b = build_base(tmp_path, monkeypatch)
    mf = b.root / "MANIFEST.json"
    man = json.loads(mf.read_bytes())
    rel = next(r for r, rec in man["files"].items() if rec.get("url") == H.turl("Ant05"))
    del man["files"][rel]
    mf.write_text(json.dumps(man))
    with pytest.raises(cp.Refused, match="not in the manifest") as e:
        cp.g0a(b.root, b.kept)
    assert not MADE_UP.search(str(e.value)), str(e.value)


# ---------------------------------------------------------------------------------------------------------------------
# Step 2: the query's class


def test_the_class_is_the_one_with_the_most_usable_experiments():
    assert cp.query_class({"a": KA, "b": LA, "c": LA}) == (LA, 2)
    assert cp.query_class({"a": KA}) == (KA, 1) and cp.query_class({}) == (None, 0)
    # a tie goes to the class name first in byte order, whatever the order of the experiments
    for us in ({"a": LA, "b": KA}, {"a": KA, "b": LA}, {"b": LA, "a": KA, "c": "Zeta", "d": "Zeta", "e": KA, "f": LA}):
        assert cp.query_class(us)[0] == KA
    assert cp.query_class({"a": "b", "b": "B"}) == ("B", 1)  # 0x42 before 0x62
    assert cp.query_class({"a": "é", "b": "z"}) == ("z", 1)  # 0x7a before 0xc3 0xa9
    assert cp.query_class({"a": "ab", "b": "a"}) == ("a", 1) and cp.query_class({"a": "a b", "b": "a-b"})[0] == "a b"


def test_a_tie_between_a_querys_classes_goes_by_byte_order(tmp_path, monkeypatch):
    lambda4 = [x.id for x in base_exps() if x.antigen == "Gata4" and x.cls == LA]
    exps = changed(base_exps(), lambda x: x.id == lambda4[0], kept=False)  # Gata4: 3 in Lambda, 3 in Kappa
    exps = changed(exps, lambda x: x.antigen == "Ant12" and x.cls == LA, bound=199)  # 11 candidates in Lambda only
    b = build_base(tmp_path, monkeypatch, exps)
    sc = cp.screen(b.root, b.kept)
    st = cp.structure(sc, "Gata4")
    assert st["class"] == KA and len(st["query_ids"]) == 3 and list(st["candidates"]) == TWELVE
    assert cp.query_class(cp.usable(sc, "Gata4")) == (KA, 3)


# ---------------------------------------------------------------------------------------------------------------------
# Step 3 and the structure


def test_the_structure_of_the_base_archive(base, opened):
    sc = cp.screen(base.root, base.kept)
    st = cp.structure(sc, "Ctnnb1")
    kappa3 = sorted(x.id for x in base.exps if x.antigen == "Ctnnb1" and x.cls == KA and x.kept)
    assert st["class"] == KA and st["query_ids"] == kappa3 and list(st["candidates"]) == TWELVE
    for a, ids in st["candidates"].items():
        assert ids == sorted(x.id for x in base.exps if x.antigen == a and x.cls == KA and x.text == CLEAN)
        assert len(ids) == 3
    sg = cp.structure(sc, "Gata4")
    assert sg["class"] == LA and len(sg["query_ids"]) == 4 and list(sg["candidates"]) == TWELVE
    assert all(len(v) == 3 for v in sg["candidates"].values())
    # reading R7: a file is opened once, and only for a query or an antigen that can reach the bar in a query's class
    assert sorted(c[0] for c in opened) == BASE_OPENED
    for q in ("Ant01", "P01", "Mab21l2", "Gata6"):
        with pytest.raises(cp.Refused, match="not a query"):
            cp.structure(sc, q)


def test_named_antigens_and_the_query_are_never_candidates(base):
    sc = cp.screen(base.root, base.kept)
    # in Kappa: Gata4 and Gata6 have three usable experiments each, Mab21l2 three rows; in Lambda, Gata6 has three
    assert len([c for c in cp.usable(sc, "Gata4").values() if c == KA]) == 3
    assert len([c for c in cp.usable(sc, "Gata6").values() if c == KA]) == 3
    for q, L in (("Ctnnb1", KA), ("Gata4", LA), ("Ctnnb1", LA), ("Gata4", KA)):
        assert list(cp.candidates(sc, q, L)) == TWELVE
    assert list(cp.candidates(sc, "Ctnnb1", XI)) == ["Ant21"] and len(cp.candidates(sc, "Ctnnb1", XI)["Ant21"]) == 20
    assert cp.candidates(sc, "Ctnnb1", MU) == {} and cp.candidates(sc, "Ctnnb1", "no such class") == {}


def test_candidates_come_in_byte_order_of_their_names(tmp_path, monkeypatch):
    exps = small_exps() + [X(f"SRX3000{k}", a, KA) for k, a in enumerate(["b1", "b1", "b1", "B2", "B2", "B2",
                                                                         "a-1", "a-1", "a-1", "a.1", "a.1", "a.1"])]
    root = build_db(tmp_path, monkeypatch, exps)
    sc = cp.screen(root, kept_for(tmp_path, exps))
    assert list(cp.candidates(sc, "Ctnnb1", KA)) == ["Ant01", "B2", "a-1", "a.1", "b1"]


# ---------------------------------------------------------------------------------------------------------------------
# Gate G0a's counts: each bar met and not met


def test_g0a_passes_on_the_base_archive_with_counts_only(base, opened, capsys):
    res = cp.g0a(base.root, base.kept)
    assert res == BASE_COUNTS
    assert sorted(c[0] for c in opened) == BASE_OPENED and len(opened) == len(BASE_OPENED)
    out = capsys.readouterr()
    assert out.out == "" and out.err == ""  # the module prints nothing
    flat = json.dumps(res)
    assert not MADE_UP.search(flat) and set(res["queries"]) == {"Ctnnb1", "Gata4"}
    for v in res["queries"].values():
        assert {type(x) for x in v.values()} == {int, bool} and set(v) == {"usable_in_class", "candidates", "a1", "a2"}
    assert type(res["classes_excluded"]) is int and type(res["passes"]) is bool


@pytest.mark.parametrize("fail,n,ok", [(("Ant12",), 11, True), (("Ant11", "Ant12"), 10, False), (TWELVE, 0, False)])
def test_bar_a2_eleven_candidates_in_the_class(tmp_path, monkeypatch, fail, n, ok):
    exps = changed(base_exps(), lambda x: x.antigen in fail and x.cls == KA, bound=199)
    b = build_base(tmp_path, monkeypatch, exps)
    res = cp.g0a(b.root, b.kept)
    assert res["queries"]["Ctnnb1"] == {"usable_in_class": 3, "candidates": n, "a1": True, "a2": ok}
    assert res["queries"]["Gata4"] == BASE_COUNTS["queries"]["Gata4"]  # its class, Lambda, is as it was
    assert res["passes"] is ok and res["classes_excluded"] == 3


@pytest.mark.parametrize("query,cls", [("Ctnnb1", KA), ("Gata4", LA)])
@pytest.mark.parametrize("left,ok", [(2, True), (1, False)])
def test_bar_a1_two_usable_experiments_in_the_class(tmp_path, monkeypatch, query, cls, left, ok):
    mine = [x.id for x in base_exps() if x.antigen == query and x.cls == cls and x.kept]
    other = [x.id for x in base_exps() if x.antigen == query and x.cls != cls and x.kept]
    # `left` of the query's experiments in its class stay usable; those of its other class fall under 200 genes too,
    # so that its class stays the one with the most
    drop = set(mine[left:]) | set(other)
    exps = changed(base_exps(), lambda x: x.id in drop, bound=150)
    b = build_base(tmp_path, monkeypatch, exps)
    res = cp.g0a(b.root, b.kept)
    assert cp.structure(cp.screen(b.root, b.kept), query)["class"] == cls
    assert res["queries"][query] == {"usable_in_class": left, "candidates": 12, "a1": ok, "a2": True}
    other_q = "Gata4" if query == "Ctnnb1" else "Ctnnb1"
    assert res["queries"][other_q] == BASE_COUNTS["queries"][other_q] and res["passes"] is ok


@pytest.mark.parametrize("how", ["too few genes", "class (v)", "excluded class", "no column", "empty kept list"])
def test_a_query_without_a_usable_experiment_meets_neither_bar(tmp_path, monkeypatch, opened, how):
    """A bar that cannot be computed is not met: no usable experiment, no class, no candidates."""
    q = lambda x: x.antigen == "Ctnnb1"  # noqa: E731
    exps = base_exps()
    ids = None
    if how == "too few genes":
        exps = changed(exps, q, bound=199)
    elif how == "class (v)":
        exps = changed(exps, q, cell_type="Sarcoma cells")
    elif how == "excluded class":
        exps = changed(exps, q, cls="Others")
    elif how == "no column":  # (the experiment the readers did not keep stays the file's one column)
        exps = changed(exps, lambda x: q(x) and x.kept, column=False)
    else:
        ids = {"Ctnnb1": [], "Gata4": kept_ids(exps)["Gata4"]}
    b = build_base(tmp_path, monkeypatch, exps)
    kept = b.kept if ids is None else kept_for(tmp_path, name="kept_empty", ids=ids)
    res = cp.g0a(b.root, kept)
    assert res["queries"]["Ctnnb1"] == {"usable_in_class": 0, "candidates": 0, "a1": False, "a2": False}
    assert res["queries"]["Gata4"] == BASE_COUNTS["queries"]["Gata4"] and res["passes"] is False
    assert res["classes_excluded"] == (4 if how == "excluded class" else 3)
    was_opened = "Ctnnb1" in [c[0] for c in opened]
    assert was_opened == (how in ("too few genes", "no column"))  # not opened when the lists already leave nothing
    st = cp.structure(cp.screen(b.root, kept), "Ctnnb1")
    assert st == {"class": None, "query_ids": [], "candidates": {}}


def test_both_bars_are_read_for_both_queries(tmp_path, monkeypatch):
    exps = changed(base_exps(), lambda x: x.antigen in ("Ant11", "Ant12") and x.cls == LA, bound=10)
    exps = changed(exps, lambda x: x.antigen == "Ctnnb1" and x.cls == KA, bound=10)
    b = build_base(tmp_path, monkeypatch, exps)
    res = cp.g0a(b.root, b.kept)
    assert res == {"classes_excluded": 3, "passes": False,
                   "queries": {"Ctnnb1": {"usable_in_class": 1, "candidates": 10, "a1": False, "a2": False},
                               "Gata4": {"usable_in_class": 4, "candidates": 10, "a1": True, "a2": False}}}


def test_bars_on_a_structure():
    def st(n, m, cls=KA):
        return {"class": cls, "query_ids": [f"SRX{i}" for i in range(n)],
                "candidates": {f"Ant{i:02d}": [] for i in range(m)}}

    assert cp.bars(st(2, 11)) == {"usable_in_class": 2, "candidates": 11, "a1": True, "a2": True}
    assert cp.bars(st(1, 11)) == {"usable_in_class": 1, "candidates": 11, "a1": False, "a2": True}
    assert cp.bars(st(2, 10)) == {"usable_in_class": 2, "candidates": 10, "a1": True, "a2": False}
    assert cp.bars(st(9, 99)) == {"usable_in_class": 9, "candidates": 99, "a1": True, "a2": True}
    assert cp.bars(st(0, 0, None)) == {"usable_in_class": 0, "candidates": 0, "a1": False, "a2": False}
    # without a class no bar can be computed, whatever else the structure holds: neither is met
    assert cp.bars(st(5, 20, None)) == {"usable_in_class": 5, "candidates": 20, "a1": False, "a2": False}


# ---------------------------------------------------------------------------------------------------------------------
# Nothing is printed, logged or written, and nothing but the one reader opens a target-gene file


def test_the_module_prints_logs_and_writes_nothing():
    tree = ast.parse((EXT / "chipq_partners.py").read_text())
    calls = [n.func for n in ast.walk(tree) if isinstance(n, ast.Call)]
    assert not {f.id for f in calls if isinstance(f, ast.Name)} & {"print", "open", "exec", "eval", "input"}
    assert not {f.attr for f in calls if isinstance(f, ast.Attribute)} & {
        "write", "write_text", "write_bytes", "mkdir", "info", "warning", "debug", "error", "exception", "run",
        "system", "read_bytes", "read_text", "urlopen"}
    mods = set()
    for n in ast.walk(tree):
        if isinstance(n, ast.Import):
            mods |= {a.name for a in n.names}
        elif isinstance(n, ast.ImportFrom):
            mods.add(n.module)
    assert mods == {"__future__", "re", "sys", "collections", "pathlib", "chipq"}


def test_no_refusal_and_no_count_carries_a_made_up_name(base, tmp_path, monkeypatch, caplog, capsys):
    """Every way the base archive's reading can end, with everything the module hands back or raises searched for the
    archive's antigens, classes, gene symbols, experiment ids and its score."""
    caplog.set_level(0)
    before = sorted(p.relative_to(base.root).as_posix() for p in base.root.rglob("*"))
    # the constants' record names the code's own exclusion list and class names, and nothing of the archive beyond
    own = re.sub(r"\b(?:Gata6|Mab21l2|Unclassified)\b", "", json.dumps(cp.readings()))
    seen = [json.dumps(cp.g0a(base.root, base.kept)), own]
    k = kept_ids(base.exps)
    for ids in ({**k, "Gata4": sorted([*k["Gata4"], "SRX99999"])}, {"Gata4": k["Gata4"]}):
        with pytest.raises(cp.Refused) as e:
            cp.g0a(base.root, kept_for(tmp_path, name=f"k{len(seen)}", ids=ids))
        seen.append(str(e.value))

    def refuses(*a, **kw):
        raise cq.Refused("a kept experiment's value is not a number: refused")

    with monkeypatch.context() as mp:
        mp.setattr(cq, "open_target", refuses)
        with pytest.raises(cp.Refused) as e:
            cp.g0a(base.root, base.kept)
        seen.append(str(e.value))
    out = capsys.readouterr()
    seen += [out.out, out.err, caplog.text]
    for text in seen:
        assert not MADE_UP.search(text), text
    assert out.out == "" and out.err == "" and caplog.text == ""
    assert sorted(p.relative_to(base.root).as_posix() for p in base.root.rglob("*")) == before  # nothing written
