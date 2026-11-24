"""Tests for experiments/t3_decomp/k28_chip_partners.py (row K28, T3 candidate CHIPQ-C: the command of the
database-only gate G0a).

Hand-made files only: the synthetic archive of tests/test_chipq_partners.py (antigens Ant01.., made-up classes, gene
symbols gs00000..), a hand-made ledger line and register line in temporary files, and a temporary git repository for
the tree and commit checks. No test reads the repository's own ledger or register, a real download or competition
data, and no byte built from released or outside data is pinned."""

from __future__ import annotations

import ast
import hashlib
import importlib.util
import json
import os
import re
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

REPO = Path(__file__).resolve().parents[1]
T3 = REPO / "experiments" / "t3_decomp"
if str(T3) not in sys.path:
    sys.path.insert(0, str(T3))
import k28_chip_partners as k28  # noqa: E402

cp, cq, cf = k28.cp, k28.cq, k28.cf


def _partners_helpers():
    """tests/test_chipq_partners.py's synthetic archive (build_base, base_exps, changed, MADE_UP), loaded under its
    own name."""
    name = "k28_partners_helpers"
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, REPO / "tests" / "test_chipq_partners.py")
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


P = _partners_helpers()
KA, LA = P.KA, P.LA
HEAD = "c" * 40
ROW_LINE = ("| K28 (one row, one script; the database-only gate G0, uncounted) | T3 | "
            "`experiments/t3_decomp/k28_chip_partners.py` | wt_identity |\n")
WORDS = "Reopened by the captain's order of 2026-11-24"
REGISTER_LINE = (f"| X018 | a generic database of public ChIP-seq | fetched 2026-11-23. {WORDS} for the database-only "
                 "gate of one practice row (K28, T3 candidate CHIPQ-C); no new fetch, no new kind of file. |\n")
LEDGER = "| id | task | code |\n|---|---|---|\n| K27 (one row, one script) | T3 | `k27_chip_query.py` |\n"
REGISTER = "| id | source | use |\n|---|---|---|\n| X017 | another source | other |\n"
BLOCK = re.compile(r"G0a: classes excluded by class \(vi\): (\d+)\n"
                   r"G0a Ctnnb1: usable experiments in its class: (\d+); candidates: (\d+)\n"
                   r"G0a Gata4: usable experiments in its class: (\d+); candidates: (\d+)\n")
BASE_BLOCK = ("G0a: classes excluded by class (vi): 3\n"
              "G0a Ctnnb1: usable experiments in its class: 3; candidates: 12\n"
              "G0a Gata4: usable experiments in its class: 4; candidates: 12\n")
S_G0 = "Try K28: STOP at G0: the database bars are not met; no fold read, no try, no file."


def make_env(tmp_path: Path, monkeypatch, exps=None, **kw) -> SimpleNamespace:
    """The step's surroundings, all under tmp_path: the synthetic archive as the fetch's data root, a ledger with row
    K28, a register whose row X018 carries the amendment's words, an empty results folder, a clean tree at HEAD."""
    b = P.build_base(tmp_path, monkeypatch, exps, **kw)
    ledger, register = tmp_path / "TRIES.md", tmp_path / "EXTERNAL_DATA.md"
    ledger.write_text(LEDGER + ROW_LINE)
    register.write_text(REGISTER + REGISTER_LINE)
    results = tmp_path / "results" / "k28_chip_partners"
    st = {"commit": HEAD, "tree_clean": True, "machine": "m1"}
    monkeypatch.setattr(k28, "TRIES_MD", ledger)
    monkeypatch.setattr(k28, "REGISTER", register)
    monkeypatch.setattr(k28, "RESULTS", results)
    monkeypatch.setattr(k28, "stamp", lambda: dict(st))
    monkeypatch.setattr(cq, "load_kept", lambda *a, **k: b.kept)
    monkeypatch.setattr(cf, "DATA_ROOT", b.root)
    monkeypatch.setattr(cf, "_opener", lambda *a: pytest.fail("a request: the step fetches nothing"))
    return SimpleNamespace(b=b, root=b.root, ledger=ledger, register=register, results=results, st=st, tmp=tmp_path)


@pytest.fixture
def env(tmp_path, monkeypatch):
    return make_env(tmp_path, monkeypatch)


def run(*argv) -> int:
    return k28.main(["g0a", "--commit", HEAD] if not argv else list(argv))


def tree(root: Path) -> dict[str, str]:
    """Every file under root: its path and the sha256 of its bytes."""
    return {p.relative_to(root).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(root.rglob("*")) if p.is_file()}


def never(what: str):
    def fail(*a, **k):
        pytest.fail(f"{what} was reached")

    return fail


# ---------------------------------------------------------------------------------------------------------------------
# The row's constants


def test_the_rows_constants():
    assert k28.S_G0 == S_G0
    assert (k28.EXIT_PASS, k28.EXIT_BAR, k28.EXIT_REFUSED) == (0, 3, 2)
    assert k28.ROW == "K28" and k28.REGISTER_ROW == "X018" and k28.REGISTER_WORDS == WORDS
    assert k28.SCRIPT == "experiments/t3_decomp/k28_chip_partners.py" and (REPO / k28.SCRIPT).is_file()
    # the results folder is named as row K27's is: the script's name, under experiments/t3_decomp/results/
    assert k28.RESULTS == T3 / "results" / "k28_chip_partners" == T3 / "results" / Path(k28.SCRIPT).stem
    assert k28.TRIES_MD == REPO / "experiments" / "harness" / "TRIES.md" and k28.REGISTER == REPO / "EXTERNAL_DATA.md"
    assert k28.QUERIES == ("Ctnnb1", "Gata4") and sorted(k28.STEPS) == ["g0a"]
    assert k28.READING_FILES == ("reading.json", "reading.txt")


def test_the_script_imports_the_database_part_only():
    """G0 reads database files only: the script imports no harness, no fold, no factor and no scoring code."""
    tree_ = ast.parse((REPO / k28.SCRIPT).read_text())
    mods = set()
    for n in ast.walk(tree_):
        if isinstance(n, ast.Import):
            mods |= {a.name for a in n.names}
        elif isinstance(n, ast.ImportFrom):
            mods.add(n.module)
    assert mods == {"__future__", "argparse", "contextlib", "hashlib", "json", "os", "re", "socket", "subprocess",
                    "sys", "pathlib", "chipq_partners"}


# ---------------------------------------------------------------------------------------------------------------------
# G0a: the printed block, the exit codes and each bar


def test_g0a_passes_with_three_lines_of_counts_and_nothing_else(env, capsys):
    assert run() == 0
    out = capsys.readouterr()
    assert out.out == BASE_BLOCK and BLOCK.fullmatch(out.out) and out.err == ""
    assert not env.results.exists() and not env.results.parent.exists()  # a pass writes nothing


def test_the_printed_lines_exact_form():
    res = {"classes_excluded": 7, "passes": False,
           "queries": {"Gata4": {"usable_in_class": 9, "candidates": 31, "a1": True, "a2": True},
                       "Ctnnb1": {"usable_in_class": 0, "candidates": 0, "a1": False, "a2": False}}}
    assert k28.g0a_lines(res) == ["G0a: classes excluded by class (vi): 7",
                                  "G0a Ctnnb1: usable experiments in its class: 0; candidates: 0",
                                  "G0a Gata4: usable experiments in its class: 9; candidates: 31"]
    m = BLOCK.fullmatch("".join(f"{ln}\n" for ln in k28.g0a_lines(res)))
    assert m and m.groups() == ("7", "0", "0", "9", "31")
    assert k28.LINE_CLASSES.count("{") == 1 and k28.LINE_QUERY.count("{") == 3  # counts, and the query's own name


def a2_not_met() -> list:
    """The base archive with ten candidates left in the test query's class (two fall under 200 genes there)."""
    return P.changed(P.base_exps(), lambda x: x.antigen in ("Ant11", "Ant12") and x.cls == KA, bound=199)


def a1_not_met(left: int) -> list:
    """The base archive with `left` usable experiments of the practice query, all in its class."""
    lam = [x.id for x in P.base_exps() if x.antigen == "Gata4" and x.cls == LA]
    return P.changed(P.base_exps(), lambda x: x.antigen == "Gata4" and (x.cls == KA or x.id in lam[left:]), bound=199)


@pytest.mark.parametrize("left", [1, 0])
def test_g0a_bar_a1_not_met_exits_3(tmp_path, monkeypatch, capsys, left):
    E = make_env(tmp_path, monkeypatch, a1_not_met(left))
    assert run() == 3
    out = capsys.readouterr()
    m = BLOCK.fullmatch(out.out)
    # with no usable experiment the practice query has no class: a bar that cannot be computed is not met
    assert m and m.groups() == ("3", "3", "12", str(left), "12" if left else "0") and out.err == ""
    rec = json.loads((E.results / "reading.json").read_text())
    assert rec["outcome"] == S_G0 and rec["g0a"]["queries"]["Gata4"] == {
        "usable_in_class": left, "candidates": 12 if left else 0, "a1": False, "a2": bool(left)}
    assert rec["g0a"]["queries"]["Ctnnb1"] == {"usable_in_class": 3, "candidates": 12, "a1": True, "a2": True}


def test_g0a_bar_a2_not_met_exits_3_and_makes_the_reading_once(tmp_path, monkeypatch, capsys):
    E = make_env(tmp_path, monkeypatch, a2_not_met())
    assert run() == 3
    out = capsys.readouterr()
    m = BLOCK.fullmatch(out.out)
    assert m and m.groups() == ("3", "3", "10", "4", "12") and out.err == ""
    # the reading: the fixed sentence, nothing counted, the counts, what was read and which code read it
    assert sorted(p.name for p in E.results.iterdir()) == ["reading.json", "reading.txt"]
    text = (E.results / "reading.txt").read_text()
    lines = text.splitlines()
    assert lines[0] == f"OUTCOME {S_G0}" and lines[1] == "tries counted: 0"
    assert "\n".join(lines[3:6]) + "\n" == out.out
    assert "bars met: Ctnnb1 a1 yes, a2 no; Gata4 a1 yes, a2 yes" in lines
    assert ("bars: a1, 2 or more usable experiments in the query's class; a2, 11 or more candidates in that class"
            in lines)
    rec = json.loads((E.results / "reading.json").read_text())
    assert rec["outcome"] == S_G0 and rec["passes"] is False and rec["tries_counted"] == 0
    assert (rec["stop"], rec["step"], rec["row"], rec["candidate"], rec["register_row"]) == (
        "G0", "g0a", "K28", "CHIPQ-C", "X018")
    assert rec["g0a"] == {"classes_excluded": 3, "passes": False,
                          "queries": {"Ctnnb1": {"usable_in_class": 3, "candidates": 10, "a1": True, "a2": False},
                                      "Gata4": {"usable_in_class": 4, "candidates": 12, "a1": True, "a2": True}}}
    assert rec["bars"] == {"a1_min_usable_in_class": 2, "a2_min_candidates": 11} and rec["readings"] == cp.readings()
    assert rec["kept_sha256"] == E.b.kept.sha256
    assert rec["fetch_manifest_sha256"] == hashlib.sha256((E.root / "MANIFEST.json").read_bytes()).hexdigest()
    code = rec["reading_code"]
    assert code["commit"] == HEAD and code["tree_clean"] is True and code["machine"] == "m1"
    assert code["script"] == k28.SCRIPT
    for key, path in (("script_sha256", REPO / k28.SCRIPT), ("chipq_partners_sha256", k28.T3X / "chipq_partners.py"),
                      ("chipq_sha256", k28.T3X / "chipq.py"), ("chipq_fetch_sha256", k28.T3X / "chipq_fetch.py")):
        assert code[key] == hashlib.sha256(path.read_bytes()).hexdigest()
    assert text == k28.reading_text(rec)
    # made once: with a reading there the candidate is closed; no step runs again and nothing is opened or rewritten
    before = tree(E.results)
    monkeypatch.setattr(cf, "load_manifest", never("a database file"))
    monkeypatch.setattr(cp, "g0a", never("the count"))
    assert run() == 2
    out = capsys.readouterr()
    assert out.out == "" and "never overwritten" in out.err and "closed" in out.err
    assert tree(E.results) == before


def test_the_reading_is_never_overwritten(tmp_path):
    out = tmp_path / "results" / "k28_chip_partners"
    st = {"commit": HEAD, "tree_clean": True, "machine": "m1"}
    counts = {"classes_excluded": 1, "passes": False,
              "queries": {q: {"usable_in_class": 1, "candidates": 2, "a1": False, "a2": False} for q in k28.QUERIES}}
    rec = k28.reading_record("g0a", counts, st, "a" * 64, "b" * 64)
    k28.write_reading(rec, out)
    first = tree(out)
    assert sorted(first) == ["reading.json", "reading.txt"]
    other = k28.reading_record("g0a", {**counts, "classes_excluded": 9}, st, "c" * 64, "d" * 64)
    with pytest.raises(k28.Refused, match="never overwritten"):
        k28.write_reading(other, out)
    assert tree(out) == first
    # either file alone, or a link under its name, stands for a reading
    for n, name in enumerate(k28.READING_FILES):
        for kind in ("file", "dangling link"):
            d = tmp_path / f"alone{n}{kind[0]}"
            d.mkdir()
            if kind == "file":
                (d / name).write_text("an earlier reading\n")
            else:
                (d / name).symlink_to(d / "nowhere")
            with pytest.raises(k28.Refused, match="never overwritten"):
                k28.write_reading(other, d)
            with pytest.raises(k28.Refused, match="never overwritten"):
                k28.reading_absent(d)
            assert sorted(p.name for p in d.iterdir()) == [name]
            if kind == "file":
                assert (d / name).read_text() == "an earlier reading\n"
    k28.reading_absent(tmp_path / "no such folder")


@pytest.mark.parametrize("name", ["reading.json", "reading.txt"])
def test_the_write_itself_never_replaces_a_file(tmp_path, monkeypatch, name):
    """Past the check too (a file that appears between the check and the write): each file is created new."""
    st = {"commit": HEAD, "tree_clean": True, "machine": "m1"}
    counts = {"classes_excluded": 1, "passes": False,
              "queries": {q: {"usable_in_class": 1, "candidates": 2, "a1": False, "a2": False} for q in k28.QUERIES}}
    out = tmp_path / "results"
    out.mkdir()
    (out / name).write_text("appeared after the check\n")
    monkeypatch.setattr(k28, "reading_absent", lambda out=None: None)
    with pytest.raises(k28.Refused, match="never overwritten"):
        k28.write_reading(k28.reading_record("g0a", counts, st, "a" * 64, "b" * 64), out)
    assert (out / name).read_text() == "appeared after the check\n"


@pytest.mark.parametrize("name", ["reading.json", "reading.txt"])
def test_a_reading_that_appears_while_the_step_counts_is_not_overwritten(tmp_path, monkeypatch, capsys, name):
    E = make_env(tmp_path, monkeypatch, a2_not_met())
    real = cp.g0a

    def late(*a, **k):
        res = real(*a, **k)
        E.results.mkdir(parents=True)
        (E.results / name).write_text("another run's reading\n")
        return res

    monkeypatch.setattr(cp, "g0a", late)
    assert run() == 2
    out = capsys.readouterr()
    assert out.out == "" and "never overwritten" in out.err
    assert {p.name: p.read_text() for p in E.results.iterdir()} == {name: "another run's reading\n"}


def test_the_reading_text_holds_the_sentence_and_counts_only():
    st = {"commit": "d" * 40, "tree_clean": True, "machine": "m2"}
    counts = {"classes_excluded": 4, "passes": False,
              "queries": {"Ctnnb1": {"usable_in_class": 1, "candidates": 12, "a1": False, "a2": True},
                          "Gata4": {"usable_in_class": 5, "candidates": 3, "a1": True, "a2": False}}}
    rec = k28.reading_record("g0a", counts, st, "a" * 64, "b" * 64)
    lines = k28.reading_text(rec).splitlines()
    assert lines[:7] == [
        f"OUTCOME {S_G0}", "tries counted: 0",
        "stopped at step g0a of gate G0 (database files only: no fold, no barred factor, no score)",
        "G0a: classes excluded by class (vi): 4",
        "G0a Ctnnb1: usable experiments in its class: 1; candidates: 12",
        "G0a Gata4: usable experiments in its class: 5; candidates: 3",
        "bars: a1, 2 or more usable experiments in the query's class; a2, 11 or more candidates in that class"]
    assert lines[7] == "bars met: Ctnnb1 a1 no, a2 yes; Gata4 a1 yes, a2 no"
    assert lines[8] == f"kept list {'a' * 64}; fetch manifest {'b' * 64}"
    assert lines[9].startswith(f"read by {k28.SCRIPT} ") and lines[9].endswith(
        f"at commit {'d' * 40}, clean True, machine m2") and len(lines) == 10
    assert json.loads(json.dumps(rec)) == rec  # plain JSON


# ---------------------------------------------------------------------------------------------------------------------
# Each refusal: exit 2, nothing printed on stdout, nothing written


def refused(E, capsys, *argv, match: str) -> str:
    before = tree(E.tmp)
    assert run(*argv) == 2
    out = capsys.readouterr()
    assert out.out == "" and out.err.startswith("refused: ") and re.search(match, out.err), out.err
    assert tree(E.tmp) == before and not E.results.exists()
    assert not P.MADE_UP.search(out.err), out.err
    return out.err


def test_the_step_refuses_before_any_database_file_without_its_row_register_tree_or_commit(env, monkeypatch, capsys):
    E = env
    monkeypatch.setattr(cf, "load_manifest", never("a database file"))
    monkeypatch.setattr(cf, "step_verify", never("the verification"))
    monkeypatch.setattr(cp, "g0a", never("the count"))
    monkeypatch.setattr(cq, "load_kept", never("the kept list"))
    # the ledger row
    for text, msg in ((LEDGER, "no row K28"), (LEDGER + ROW_LINE + ROW_LINE, "more than one row K28"),
                      (LEDGER + ROW_LINE.replace("| K28 ", "| K280 "), "no row K28"),
                      (LEDGER + ROW_LINE.replace("| K28 ", "| K28b "), "no row K28"),
                      (LEDGER + "the row K28 is described here, in a paragraph\n", "no row K28"),
                      (LEDGER + " " + ROW_LINE, "no row K28")):
        E.ledger.write_text(text)
        refused(E, capsys, match=msg)
    E.ledger.unlink()
    refused(E, capsys, match="TRIES.md is missing")
    E.ledger.write_text(LEDGER + ROW_LINE)
    # the register's wording
    plain = REGISTER_LINE.replace(WORDS, "Opened on 2026-11-23")
    for text, msg in ((REGISTER + plain, "opening words"),
                      (REGISTER + plain + f"| X019 | another | {WORDS} |\n", "opening words"),
                      (REGISTER + plain + f"\n{WORDS} for row X018.\n", "opening words"),
                      (REGISTER + REGISTER_LINE.replace("2026-11-24", "2026-11-23"), "opening words"),
                      (REGISTER + REGISTER_LINE.replace("captain's", "captains"), "opening words"),
                      (REGISTER, "no row X018"), (REGISTER + REGISTER_LINE + REGISTER_LINE, "more than one row X018")):
        E.register.write_text(text)
        refused(E, capsys, match=msg)
    E.register.unlink()
    refused(E, capsys, match="EXTERNAL_DATA.md is missing")
    E.register.write_text(REGISTER + REGISTER_LINE)
    # the tree and the commit
    E.st["tree_clean"] = False
    refused(E, capsys, match="the tree is not clean")
    E.st["tree_clean"] = True
    for commit in ("d" * 40, HEAD[:7], HEAD[:39], HEAD + "c", HEAD.upper(), "HEAD", ""):
        refused(E, capsys, "g0a", "--commit", commit, match="HEAD is not the commit given")
    E.st["commit"] = None
    refused(E, capsys, match="HEAD is not the commit given")


def test_the_step_passes_again_once_what_it_refused_for_is_in_place(env, capsys):
    env.ledger.write_text(LEDGER)
    refused(env, capsys, match="no row K28")
    env.ledger.write_text(LEDGER + ROW_LINE)
    assert run() == 0 and capsys.readouterr().out == BASE_BLOCK


def test_a_file_off_its_sha256_refuses_before_anything_is_counted(env, monkeypatch, capsys):
    E = env
    monkeypatch.setattr(cp, "g0a", never("the count"))
    monkeypatch.setattr(cq, "load_kept", never("the kept list"))
    man = json.loads((E.root / "MANIFEST.json").read_text())
    # a list; a target-gene file that the count itself would never open (its antigen's class is excluded); a record
    unopened = next(r for r, rec in man["files"].items() if rec.get("url") == P.H.turl("Ant16"))
    for rel in (f"lists/{cf.EXPERIMENTS}", f"lists/{cf.ANALYSIS}", unopened, "screen/placebo_kept.json", cf.UNFETCHED):
        p = E.root / rel
        raw = p.read_bytes()
        p.write_bytes(raw[:-1] + bytes([raw[-1] ^ 1]))
        err = refused(E, capsys, match="off its sha256")
        assert rel not in err
        p.write_bytes(raw)
        p.rename(p.with_name(p.name + ".away"))
        refused(E, capsys, match="off its sha256")
        p.with_name(p.name + ".away").rename(p)
    # an empty manifest, and none at all
    mf = E.root / "MANIFEST.json"
    raw = mf.read_bytes()
    mf.write_text(json.dumps({"root": "x", "files": {}}))
    refused(E, capsys, match="off its sha256")
    mf.rename(E.root / "MANIFEST.json.away")
    refused(E, capsys, match="off its sha256")
    (E.root / "MANIFEST.json.away").rename(mf)
    mf.write_bytes(raw)


def test_the_manifest_is_verified_first_then_the_kept_list_then_the_count(env, monkeypatch, capsys):
    order = []
    real_verify, real_g0a, kept = cf.step_verify, cp.g0a, env.b.kept
    monkeypatch.setattr(k28, "preconditions", lambda c: (order.append("preconditions"), dict(env.st))[1])
    monkeypatch.setattr(cf, "step_verify", lambda root: (order.append("verify"), real_verify(root))[1])
    monkeypatch.setattr(cq, "load_kept", lambda *a, **k: (order.append("kept"), kept)[1])
    monkeypatch.setattr(cp, "g0a", lambda *a, **k: (order.append("count"), real_g0a(*a, **k))[1])
    assert run() == 0 and order == ["preconditions", "verify", "kept", "count"]
    assert capsys.readouterr().out == BASE_BLOCK


def test_the_step_refuses_without_the_kept_list_or_a_complete_fetch(tmp_path, monkeypatch, capsys):
    E = make_env(tmp_path, monkeypatch)

    def absent(*a, **k):
        raise cq.Refused("chipq_kept.yaml is absent (the main session commits it after the screen): refused")

    with monkeypatch.context() as mp:
        mp.setattr(cq, "load_kept", absent)
        mp.setattr(cp, "g0a", never("the count"))
        refused(E, capsys, match="chipq_kept.yaml is absent")
    # a fetch whose targets pass never finished
    (tmp_path / "two").mkdir()
    E2 = make_env(tmp_path / "two", monkeypatch, P.small_exps(), targets=False)
    refused(E2, capsys, match="targets pass is incomplete")
    # a data root inside a checkout
    monkeypatch.setattr(cf, "DATA_ROOT", REPO / "data" / "k28_no_such_folder")
    refused(E2, capsys, match="inside a checkout")


def test_an_uncommitted_kept_list_is_refused(tmp_path, monkeypatch, capsys):
    """The step loads the kept list with chipq.load_kept's own check: tracked by git and unchanged against HEAD."""
    E = make_env(tmp_path, monkeypatch)
    yp, sp = P.H.kept_files(tmp_path / "kept", P.kept_ids(E.b.exps))
    monkeypatch.undo()  # chipq.load_kept is the real one again
    monkeypatch.setattr(cf, "_sleep", lambda s: None)
    for name, value in (("TRIES_MD", E.ledger), ("REGISTER", E.register), ("RESULTS", E.results)):
        monkeypatch.setattr(k28, name, value)
    monkeypatch.setattr(k28, "stamp", lambda: dict(E.st))
    monkeypatch.setattr(cf, "DATA_ROOT", E.root)
    monkeypatch.setattr(cq, "KEPT_YAML", yp)
    monkeypatch.setattr(cq, "KEPT_SHA", sp)
    monkeypatch.setattr(cp, "g0a", never("the count"))
    refused(E, capsys, match="is not committed")
    monkeypatch.setattr(cq, "KEPT_YAML", tmp_path / "missing.yaml")
    refused(E, capsys, match="is absent")


def test_a_crash_shows_its_type_and_never_its_message(env, monkeypatch, capsys):
    def crash(*a, **k):
        raise ValueError("could not convert 'Ant01 gs00017 SRX10042 7777.5 Kappa tissue'")

    monkeypatch.setattr(cp, "g0a", crash)
    err = refused(env, capsys, match="ValueError in step g0a")
    assert err == "refused: ValueError in step g0a\n"

    def killed(*a, **k):
        raise MemoryError

    monkeypatch.setattr(cp, "g0a", killed)
    assert refused(env, capsys, match="MemoryError") == "refused: MemoryError in step g0a\n"


def test_usage(env, capsys):
    for argv in (["g0a"], ["g0b", "--commit", HEAD], ["read", "--commit", HEAD], ["--commit", HEAD],
                 ["g0a", "--commit", HEAD, "--root", str(env.tmp)], ["g0a", "--commit", HEAD, "--out", str(env.tmp)]):
        with pytest.raises(SystemExit) as ex:
            k28.main(argv)
        assert ex.value.code == 2 and capsys.readouterr().out == ""  # no path is taken from the command line
    assert not env.results.exists()


def test_what_a_loader_prints_never_reaches_stdout(env, monkeypatch, capsys):
    kept = env.b.kept
    monkeypatch.setattr(cq, "load_kept", lambda *a, **k: (print("noise from a loader"), kept)[1])
    assert run() == 0
    out = capsys.readouterr()
    assert out.out == BASE_BLOCK and out.err == "noise from a loader\n"


# ---------------------------------------------------------------------------------------------------------------------
# The ledger row and the register's wording, as functions


def test_row_present_and_register_reopened(tmp_path):
    led, reg = tmp_path / "TRIES.md", tmp_path / "EXTERNAL_DATA.md"
    led.write_text(LEDGER + ROW_LINE + "| K29 (another row) | T3 |\n")
    reg.write_text(REGISTER + REGISTER_LINE)
    assert k28.row_present(led) == "K28" and k28.register_reopened(reg) is None
    led.write_text(LEDGER + "| K28 | T3 | the shortest form of the row |\n")
    assert k28.row_present(led) == "K28"
    led.write_text(LEDGER)
    with pytest.raises(k28.Refused, match="no row K28 in TRIES.md"):
        k28.row_present(led)
    (tmp_path / "link.md").symlink_to(reg)
    with pytest.raises(k28.Refused, match="missing or a symlink"):
        k28.register_reopened(tmp_path / "link.md")
    with pytest.raises(k28.Refused, match="missing or a symlink"):
        k28.row_present(tmp_path / "nothing.md")
    reg.write_text(REGISTER + "| X018 | a generic database | fetched 2026-11-23 |\n")
    with pytest.raises(k28.Refused, match="opening words"):
        k28.register_reopened(reg)


# ---------------------------------------------------------------------------------------------------------------------
# The tree and the commit, on a real temporary repository

GIT_VARS = ("GIT_DIR", "GIT_WORK_TREE", "GIT_INDEX_FILE", "GIT_OBJECT_DIRECTORY", "GIT_COMMON_DIR", "GIT_NAMESPACE",
            "GIT_ALTERNATE_OBJECT_DIRECTORIES", "GIT_CEILING_DIRECTORIES")


def _git(repo: Path, *args) -> str:
    env = {k: v for k, v in os.environ.items() if k not in GIT_VARS}
    env |= {"GIT_CONFIG_GLOBAL": os.devnull, "GIT_CONFIG_NOSYSTEM": "1", "GIT_CEILING_DIRECTORIES": str(repo.parent)}
    return subprocess.run(["git", "-C", str(repo), "-c", "user.name=t", "-c", "user.email=t@example.org",
                           "-c", "commit.gpgsign=false", *args], check=True, capture_output=True, text=True,
                          env=env).stdout.strip()


def test_the_steps_git_calls_never_follow_the_shells_git_variables():
    """A shell that exports GIT_DIR, GIT_WORK_TREE or GIT_INDEX_FILE does not point the step at another repository:
    the three are left out of every git call's environment (shown on a plain mapping; no variable is set here)."""
    shell = {"PATH": "/usr/bin", "HOME": "/home/x", "GIT_DIR": "/other/.git", "GIT_WORK_TREE": "/other",
             "GIT_INDEX_FILE": "/other/.git/index", "GIT_CEILING_DIRECTORIES": "/tmp"}
    assert k28._git_env(shell) == {"PATH": "/usr/bin", "HOME": "/home/x", "GIT_CEILING_DIRECTORIES": "/tmp",
                                   "GIT_OPTIONAL_LOCKS": "0"}
    assert not {"GIT_DIR", "GIT_WORK_TREE", "GIT_INDEX_FILE"} & set(k28._git_env())
    src = (REPO / k28.SCRIPT).read_text()
    assert src.count("subprocess.run(") == 1 and "env=_git_env()" in src and '"-C", str(REPO)' in src


def test_the_stamp_reads_the_checkouts_own_tree_and_commit(tmp_path, monkeypatch):
    repo, other = tmp_path / "repo", tmp_path / "other"
    for r, text in ((repo, "one\n"), (other, "two\n")):
        r.mkdir()
        _git(r, "init", "-q")
        (r / "tracked.txt").write_text(text)
        _git(r, "add", ".")
        _git(r, "commit", "-q", "-m", "a commit")
    head = _git(repo, "rev-parse", "HEAD")
    assert re.fullmatch(r"[0-9a-f]{40}", head) and head != _git(other, "rev-parse", "HEAD")
    monkeypatch.setattr(k28, "REPO", repo)
    monkeypatch.setenv("GIT_CEILING_DIRECTORIES", str(tmp_path))
    st = k28.stamp()
    assert st["commit"] == head and st["tree_clean"] is True and st["machine"]
    # an untracked file (the reading's own folder is one) leaves the tree clean; a changed tracked file does not
    (repo / "experiments").mkdir()
    (repo / "experiments" / "reading.txt").write_text("untracked\n")
    assert k28.stamp()["tree_clean"] is True
    led, reg = tmp_path / "TRIES.md", tmp_path / "EXTERNAL_DATA.md"
    led.write_text(LEDGER + ROW_LINE)
    reg.write_text(REGISTER + REGISTER_LINE)
    for name, value in (("TRIES_MD", led), ("REGISTER", reg), ("RESULTS", tmp_path / "results")):
        monkeypatch.setattr(k28, name, value)
    assert k28.preconditions(head) == st
    for commit in (head[:12], _git(other, "rev-parse", "HEAD"), "main", "HEAD"):
        with pytest.raises(k28.Refused, match="HEAD is not the commit given"):
            k28.preconditions(commit)
    (repo / "tracked.txt").write_text("changed\n")
    assert k28.stamp()["tree_clean"] is False
    with pytest.raises(k28.Refused, match="the tree is not clean"):
        k28.preconditions(head)
    _git(repo, "checkout", "--", "tracked.txt")
    assert k28.preconditions(head) == st
    # a staged change is a change; a folder that is no repository has no commit and is not clean
    (repo / "tracked.txt").write_text("staged\n")
    _git(repo, "add", "tracked.txt")
    assert k28.stamp()["tree_clean"] is False
    bare = tmp_path / "not_a_repository"
    bare.mkdir()
    monkeypatch.setattr(k28, "REPO", bare)
    assert k28.stamp() | {"machine": ""} == {"commit": None, "tree_clean": False, "machine": ""}
    with pytest.raises(k28.Refused, match="the tree is not clean"):
        k28.preconditions(head)


# ---------------------------------------------------------------------------------------------------------------------
# No name from a hand-made file reaches stdout, stderr, a log or a written file


def test_no_made_up_name_reaches_stdout_stderr_a_log_or_a_written_file(tmp_path, monkeypatch, capsys, caplog):
    """A pass, a stop at each bar with its reading, a second run over a reading, and refusals raised from inside the
    count, on the synthetic archive: what is printed, what is logged and every file written hold none of the archive's
    antigens, classes, gene symbols, experiment ids or its score. Nothing else is written, and the pinned copy, the
    kept list, the ledger and the register are left byte for byte as they were."""
    caplog.set_level(0)
    texts = []

    def collect() -> None:
        out = capsys.readouterr()
        texts.extend([out.out, out.err])

    reading = ["results/k28_chip_partners/reading.json", "results/k28_chip_partners/reading.txt"]
    for n, (exps, code) in enumerate(((P.base_exps(), 0), (a2_not_met(), 3), (a1_not_met(1), 3), (a1_not_met(0), 3))):
        d = tmp_path / f"case{n}"
        d.mkdir()
        E = make_env(d, monkeypatch, exps)
        before = tree(d)
        assert run() == code
        collect()
        after = tree(d)
        assert {k: after[k] for k in before} == before and sorted(set(after) - set(before)) == (reading if code else [])
        if code:
            rec = json.loads((d / reading[0]).read_text())
            assert rec.pop("readings") == cp.readings()  # the code's own constants (its exclusion list among them)
            texts.extend([json.dumps(rec), (d / reading[1]).read_text()])
            assert run() == 2  # the reading is there: refused
            collect()
            assert tree(d) == after
            continue
        # refusals raised from inside the count: a file that changed after the verification, and a kept value that is
        # not a number in a file the manifest holds
        man = json.loads((E.root / "MANIFEST.json").read_text())
        p = E.root / next(r for r, rec in man["files"].items() if rec.get("url") == P.H.turl("Ant03"))
        raw = p.read_bytes()
        monkeypatch.setattr(cf, "step_verify", lambda root: True)
        p.write_bytes(raw + b"x")
        assert run() == 2
        collect()
        P.H._rewrite_target(E.root, "Ant03", raw.replace(b"7777.5", b"oops", 1))
        assert run() == 2
        collect()
        assert not (d / "results").exists()
    assert caplog.text == ""
    for t in texts:
        clean = re.sub(r"[0-9a-f]{12,64}", "", t)  # sha256 values and the commit
        assert not P.MADE_UP.search(clean), clean
        assert not re.search(r"Kappa|Lambda|tissue|oops", clean), clean
    for said in ("STOP at G0", "differ from the manifest", "not a number", "never overwritten", "a2 no", "a1 no"):
        assert any(said in t for t in texts), said
