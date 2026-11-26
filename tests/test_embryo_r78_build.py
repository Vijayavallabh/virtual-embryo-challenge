"""Row R78's two files (the captain's order of 2026-11-26: the embryo's fourth refinement round on O, both sent at once
with no gate): the builder's one table EMBRYO_R78_TABLE (experiments/builds/build_test.py embryo) and the base checker
experiments/builds/embryo_r78_check.py, on synthetic stages and a synthetic test board whose limits are never shown (the
synthetic world of tests/test_mix_lean_build.py). Nothing is read from data/ or boards/; a synthetic TRIES.md stands in
for the ledger (two tests ask, by code and printing nothing, whether the repository's row R78 names the table as the
builder needs it).

Pinned: the table is the row's two pairs in its order, each name O's with its lean's exact fraction as its mark, and no
string or name of it is an earlier row's; the design's lean name later2_half and the struck file's string (series at
5/6) are not taken, and no lean gives t_mix 1 (the E8.0 twin, still declined); each lean passes its exact t_mix's
nearest double to vec.t2_build.interpolate_clouds and nothing else changes in O's call, the count of E7.5 rows being
interpolate_clouds' own int(round((1 - t_mix) n)), and every row keeps its bracket's released expression; each file's
record gains the r78 block and shows no count; each is taken only with --captain-ruling R78 on O's command with exactly
its string, at the default seed on the test board, while row R78's first cell names the string paired with the file and
the ruling, every refusal before anything is read; a lean only row R78 has (later_11_12, earlier) is refused under
every other ruling and with none; row R78's check runs first in build_embryo and, when it returns a file, row R75's,
row R59's and series_back's checks are skipped; the checker gives O's bytes from each file's own path with its t_mix put
back at 1/2 and fails when the route is not as named; and every other embryo command gives what the builder of BASE
(main with row R78 placed, before its route) gives, run in this same process from git's copy of build_test.py: the same
h5ad bytes and record, or the same refusal, but for the commands row R78 itself changes (its two files, now built; its
leans under any other ruling or none, and any other command under --captain-ruling R78, refused at BASE and now by row
R78's check). Built through experiments/builds/build_test.py."""

from __future__ import annotations

import hashlib
import inspect
import json
import re
import sys
import types
from fractions import Fraction
from pathlib import Path

import numpy as np
import pytest
from conftest import git
from test_embryo_r75_build import ARGS as R75_ARGS
from test_embryo_r75_build import ROW_R75
from test_mix_lean_build import R48_ROW, R59_ROW, R60_RS_ROW, SECRET_MAX, SECRET_TRUTH, _file, _shows, _stage, install

REPO = Path(__file__).resolve().parents[1]
BASE = "4d91239fc6537a57004a424fa0fb566732cef12d"  # main with rows R77, R78 and R79 placed, before row R78's route
FLOOR = ["--method", "fgw_shape0", "--positions", "floor_exact"]
OV = ["--shape-interp", "overlay", "--shape-interp-test-ok"]
O_ARGS = [*FLOOR, *OV]
O_CMD = "embryo " + " ".join(O_ARGS)
OT_ARGS = [*O_ARGS, "--mix-lean", "later", "--captain-ruling", "R59"]
SO_ARGS = [*FLOOR, "--expression", "series", *OV]
W_ARGS = [*FLOOR, "--expression", "series_back", *OV, "--captain-ruling", "R48"]
DRY = "T2:embryo:val_interp"
TABLE = [("T2embryo_s2x_ovt1112.h5ad", "--mix-lean later_11_12"),
         ("T2embryo_s2x_ovt13.h5ad", "--mix-lean earlier")]
OPTS = dict(TABLE)
NAMES = list(OPTS)
ARGS = {f: [*O_ARGS, *o.split(), "--captain-ruling", "R78"] for f, o in TABLE}
T_MIX = {"T2embryo_s2x_ovt1112.h5ad": Fraction(11, 12), "T2embryo_s2x_ovt13.h5ad": Fraction(1, 3)}
LEAN = {"T2embryo_s2x_ovt1112.h5ad": "later_11_12", "T2embryo_s2x_ovt13.h5ad": "earlier"}
# row R78's first cell as the repository's ledger holds it (tests below check the copy against it by code)
FIRST_R78 = (
    "R78 (the captain's order of 2026-11-26, ruling, no run: \"Yes to T2embryo_s2x_ovt1112.h5ad under R22, sent "
    "at once in one ruling row, with no gate: my named exception on the embryo's overlay form and occupancy "
    "term for this file.\" and \"Yes to T2embryo_s2x_ovt13.h5ad under R22 in the same row, sent at once, with no "
    "gate: my named exception on the embryo's overlay form and occupancy term for this file, and the lean-early "
    "mixing line declined on 2026-11-06 (e) reopened for this file alone.\" (their pick \"Both, with exceptions\" "
    "in a dialog with no recommended option, the AI assistant being recused on the embryo's occupancy term and "
    "shape forms; the option recorded the words above as the captain's); the embryo's fourth refinement round "
    "on O (EMBRYO_REF4_design.md with EMBRYO_REF4_challenge.md's fixes, under which both are sent at once; the "
    "design's second file, SOt56, struck as a re-draw in substance); the embryo builder is to take on `embryo "
    "--method fgw_shape0 --positions floor_exact --shape-interp overlay --shape-interp-test-ok` exactly these "
    "option strings, each only with `--captain-ruling R78`: `--mix-lean later_11_12` "
    "(T2embryo_s2x_ovt1112.h5ad), `--mix-lean earlier` (T2embryo_s2x_ovt13.h5ad); t_mix 11/12 and 1/3; both "
    "brackets' tissue kept, depth moved only through the shares, the E8.0 twin (t_mix 1) still declined; "
    "decision log 2026-11-26 (h); not a try)")
ROW_R78 = f"| {FIRST_R78} | T2 embryo | synthetic |"
LEDGER = (R48_ROW, R59_ROW, R60_RS_ROW, ROW_R75, ROW_R78)
RULINGS_OTHER = [None, "R59", "R48", "R71", "R60", "R75", "R76", "R77", "R79", "R780", "r78", "R78a", "R7", "E58"]
SOT56 = "--expression series --mix-lean later2"  # the design's second file, struck as a re-draw in substance


def _pairs(cell: str) -> list[tuple[str, str]]:
    return [(f, o) for o, f in re.findall(r"`([^`]*)` \((T2embryo_[A-Za-z0-9_]+\.h5ad)\)", cell)]


def _values(name: str) -> dict:
    return {"file": name, "option": OPTS[name], "expression": "mixture", "mix_lean": LEAN[name], "t_mix": T_MIX[name]}


@pytest.fixture(scope="module")
def stages():
    return {"E7.25.h5ad": _stage("E7.25", 11), "E7.5.h5ad": _stage("E7.5", 12), "E8.0.h5ad": _stage("E8.0", 13)}


@pytest.fixture(scope="module")
def bt():
    sys.path.insert(0, str(REPO / "experiments" / "builds"))
    import build_test

    return build_test


@pytest.fixture(scope="module")
def chk(bt):
    import embryo_r78_check

    assert embryo_r78_check.bt is bt  # the checker drives this same builder module
    return embryo_r78_check


@pytest.fixture
def env(bt, monkeypatch, own_repo, tmp_path, stages):
    """The synthetic stages, boards, size JSON and ledger (rows R48, R59, R60's synthetic ones and rows R75's and R78's
    first cells) installed in the builder; returns the list of stage reads."""
    reads: list[str] = []
    install(bt, monkeypatch.setattr, monkeypatch.setitem, tmp_path / "shared", reads, stages, ledger_rows=LEDGER)
    return reads


def _build(mod, tmp_path: Path, folder: str, *argv: str, name: str = "T2embryo.h5ad") -> tuple[Path, dict, str]:
    out = tmp_path / folder / name
    mod.main(["embryo", *argv, "--out", str(out)])
    text = out.with_suffix(".build.json").read_text()
    return out, json.loads(text), text


def _capture(bt, monkeypatch) -> list:
    """Record each interpolation the builder makes: its arguments and its result."""
    seen = []
    inner = bt.t2_build.interpolate_clouds

    def spy(*a, **kw):
        res = inner(*a, **kw)
        seen.append((a, kw, res))
        return res

    monkeypatch.setattr(bt.t2_build, "interpolate_clouds", spy)
    return seen


def _refused(bt, env, tmp_path, argv, match: str | None, name="T2embryo.h5ad") -> str:
    out = tmp_path / "r" / name
    with pytest.raises(SystemExit, match=match) as e:
        bt.main(["embryo", *argv, "--out", str(out)])
    assert not out.exists() and env == []  # refused before any stage was read
    assert not _shows(str(e.value), SECRET_MAX, SECRET_TRUTH)
    return str(e.value)


# ------------------------------------------------------------------------------------------------ the one table
def test_the_one_table_is_the_rows_two_pairs_in_its_order(bt):
    assert [tuple(p) for p in bt.EMBRYO_R78_TABLE] == TABLE == _pairs(FIRST_R78)
    assert re.findall(r"`([^`]*)`", FIRST_R78) == [O_CMD, "--captain-ruling R78", *OPTS.values()]
    assert bt.EMBRYO_R78_ROW == "R78" and list(bt.EMBRYO_R78_FILES) == NAMES
    for name, (option, expression, lean, t_mix) in bt.EMBRYO_R78_FILES.items():
        assert (option, expression, lean, t_mix) == (OPTS[name], "mixture", LEAN[name], T_MIX[name])
        assert type(t_mix) is Fraction  # exact fractions, not floats
        assert bt.EMBRYO_R78_COMMAND.format(option=option) == "embryo " + " ".join(ARGS[name])
        assert bt._r78_values(name, option) == ("mixture", LEAN[name], T_MIX[name])
    assert bt.EMBRYO_R78_LEANS == {"later_11_12": Fraction(11, 12), "earlier": Fraction(1, 3)}
    assert "t_mix 11/12 and 1/3;" in FIRST_R78 and [T_MIX[f] for f in NAMES] == [Fraction(11, 12), Fraction(1, 3)]
    assert [float(T_MIX[f]) for f in NAMES] == [11 / 12, 1 / 3]
    assert bt.EMBRYO_R78_RULES == {f: bt._r75_rule("mixture", T_MIX[f]) for f in NAMES}  # row R75's rule, its t_mix


def test_the_parser_takes_row_r78s_two_leans_and_no_other(bt, capsys):
    assert bt.EMBRYO_R78_NEW_LEANS == ("later_11_12", "earlier") and bt.MIX_LEANS == ("later",)
    action = {a.dest: a for a in bt.build_parser()._subparsers._group_actions[0].choices["embryo"]._actions}
    assert action["mix_lean"].choices == ("later", "later_half", "later2", "earlier_half", "later_11_12", "earlier")
    assert action["mix_lean"].default is None and action["captain_ruling"].default is None
    # later2_half: the design's name for 11/12, renamed by the challenge; no lean gives t_mix 1 (the E8.0 twin)
    for value in ("later2_half", "Earlier", "later_11/12", "later_1112", "11/12", "1/3", "earlier_third", "later_1",
                  "twin", "later_11_12 ", ""):
        with pytest.raises(SystemExit) as e:
            bt.build_parser().parse_args(["embryo", *O_ARGS, "--mix-lean", value, "--captain-ruling", "R78",
                                          "--out", "x.h5ad"])
        assert e.value.code == 2 and "invalid choice" in capsys.readouterr().err


def test_no_string_or_name_of_row_r78_is_an_earlier_rows(bt):
    earlier_files = {bt.EMBRYO_OT23_FILE, bt.EMBRYO_OVSB_FILE, bt.EMBRYO_RS_FILE, *bt.EMBRYO_BF_COMMANDS,
                     *bt.EMBRYO_R75_FILES, "T2embryo_s2x_ov.h5ad", "T2embryo_s2xs_ov.h5ad", "T2embryo_s2x_ovhi.h5ad",
                     "T2embryo_s2x.h5ad", "T2embryo_s2xs.h5ad", "T2embryo_s2xs_ovt56.h5ad"}
    assert not set(NAMES) & earlier_files
    commands = {bt.EMBRYO_R78_COMMAND.format(option=o) for o in OPTS.values()}
    assert not commands & {bt.EMBRYO_OT23_COMMAND, bt.EMBRYO_OVSB_COMMAND, *bt.EMBRYO_BF_COMMANDS.values(),
                           bt.EMBRYO_RS_COMMAND.format(row="R78"),
                           *(bt.EMBRYO_R75_COMMAND.format(option=o) for o, *_ in bt.EMBRYO_R75_FILES.values())}
    assert not set(OPTS.values()) & ({o for o, *_ in bt.EMBRYO_R75_FILES.values()} | {"--mix-lean later", SOT56})
    assert not set(bt.EMBRYO_R78_LEANS) & {*bt.MIX_LEANS, *bt.EMBRYO_R75_LEANS}
    assert set(bt.EMBRYO_R78_LEANS.values()).isdisjoint({*bt.EMBRYO_R75_LEANS.values(), Fraction(1, 2), Fraction(1)})
    assert "later2_half" not in bt.EMBRYO_R78_LEANS and SOT56 not in OPTS.values()  # renamed; struck


@pytest.mark.parametrize("name, option", [
    ("T2embryo_s2x_ovt1112.h5ad", "--mix-lean later2_half"),  # the design's name, renamed by the challenge
    ("T2embryo_s2xs_ovt56.h5ad", SOT56),  # the design's second file, struck
    ("T2embryo_s2x_ovt56.h5ad", "--mix-lean later2"),  # row R75's own pair: not row R78's form
    ("T2embryo_s2x_ovt512.h5ad", "--mix-lean earlier_half"),
    ("T2embryo_s2x_ovt23.h5ad", "--mix-lean later"),  # row R59's
    ("T2embryo_s2x_ovt11.h5ad", "--mix-lean later_1"),  # the E8.0 twin: no such lean
    ("T2embryo_s2x_ov.h5ad", ""),  # O's own: no string
    ("T2embryo_s2x_ovt1112.h5ad", "--expression mixture --mix-lean later_11_12"),  # O's default named: not the string
    ("T2embryo_s2x_ovt1112.h5ad", "--mix-lean  later_11_12"),  # two spaces
    ("T2embryo_s2x_ovt1112.h5ad", " --mix-lean later_11_12"),
    ("T2embryo_s2x_ovt1112.h5ad", "--mix-lean later_11_12 "),
    ("T2embryo_s2x_ovt1112.h5ad", "--mix-lean\tlater_11_12"),
    ("T2embryo_s2x_ovt1112.h5ad", "--mix-lean Later_11_12"),
    ("T2embryo_s2x_ovt1112.h5ad", "--mix-lean 11/12"),
    ("T2embryo_s2x_ovt1112.h5ad", "--mix-lean=later_11_12"),
    ("T2embryo_s2x_ovt1112.h5ad", "--mix-lean later_11_12 --mix-lean later_11_12"),
    ("T2embryo_s2x_ovt1112.h5ad", "--mix-lean later_11_12 --herd"),
    ("T2embryo_s2x_ovt1112.h5ad", "--mix-lean"),
    ("T2embryo_s2x_ovt1112.h5ad", "--mix-lean earlier"),  # the other file's string
    ("T2embryo_s2x_ovt13.h5ad", "--mix-lean later_11_12"),
    ("T2embryo_s2x_ovt13.h5ad", "--mix-lean Earlier"),
    ("T2embryo_s2x_ovt13.h5ad", "--mix-lean 1/3"),
    ("T2embryo_s2x_ovt13.h5ad", "--mix-lean earlier --expression series"),
    ("T2embryo_s2x_ovt2224.h5ad", "--mix-lean later_11_12"),  # 11/12 named unreduced
    ("T2embryo_s2x_ovt26.h5ad", "--mix-lean earlier"),  # 1/3 named unreduced
    ("T2embryo_s2x_ovt0.92.h5ad", "--mix-lean later_11_12"),
    ("T2embryo_s2x_ovt0.33.h5ad", "--mix-lean earlier"),
    ("T2embryo_s2xs_ovt13.h5ad", "--mix-lean earlier"),  # a series mark without series
    ("T2embryo_s2x_ovt13d.h5ad", "--mix-lean earlier"),  # a depth mark without depth
    ("T2embryo_s2x_ovt13", "--mix-lean earlier"),  # no suffix
    ("t2embryo_s2x_ovt13.h5ad", "--mix-lean earlier"),
    ("T2embryo_s2x_ovT13.h5ad", "--mix-lean earlier"),
    ("T2embryo_s2x_ovt1112.H5AD", "--mix-lean later_11_12"),
    ("T2embryo_s2x_ovt13.h5ad.h5ad", "--mix-lean earlier"),
])
def test_a_table_entry_off_the_rows_form_is_not_read(bt, name, option):
    with pytest.raises(ValueError):
        bt._r78_values(name, option)


def test_the_real_rows_first_cell_pairs_the_tables_strings(bt):
    """Row R78's first cell in the repository's ledger, read by code and printed nowhere: exactly one row R78, its first
    cell pairing the table's two strings with their files in the table's order, naming in backticks what the copy
    names, word for word the copy, and passing the builder's own check for each file."""
    text = (REPO / "experiments" / "harness" / "TRIES.md").read_text(encoding="utf-8")
    cells = [ln[2:].split("|", 1)[0].strip() for ln in text.splitlines() if ln.startswith("| R78 ")]
    ok = (len(cells) == 1 and _pairs(cells[0]) == TABLE
          and re.findall(r"`([^`]*)`", cells[0]) == re.findall(r"`([^`]*)`", FIRST_R78)
          and cells[0].startswith("R78 (the captain's order") and "ruling, no run" in cells[0]
          and all(bt._r78_row_names(f) for f in NAMES))
    assert ok, "the repository's row R78 does not name the table as the copy does (nothing of it printed)"
    same = len(cells) == 1 and cells[0] == FIRST_R78
    assert same, "the copy is not the row's first cell (nothing of it printed)"


def test_the_ruling_check_takes_only_a_single_row_r78_pairing_the_whole_string_with_its_file(bt, tmp_path):
    p = tmp_path / "TRIES.md"
    p.write_text(f"# ledger\n{ROW_R78}\n")
    assert all(bt._r78_row_names(f, p) for f in NAMES)
    for name, option in TABLE:
        other = [f for f in NAMES if f != name][0]
        for bad in (ROW_R78.replace(f"`{option}`", option),  # not backticked
                    ROW_R78.replace(f"`{option}`", f"`{option} --herd`"),  # inside a longer string
                    ROW_R78.replace(f"`{option}` ({name})", f"`{option}` ({other})"),  # paired with another file
                    ROW_R78.replace(f"`{option}` ({name})", f"`{option}`  ({name})")):  # not followed by " (<file>)"
            p.write_text(f"# ledger\n{bad}\n")
            assert [bt._r78_row_names(f, p) for f in NAMES] == [f != name for f in NAMES], (name, bad[-80:])
    for bad in (ROW_R78.replace("`--captain-ruling R78`", "--captain-ruling R78"),
                ROW_R78.replace("`--captain-ruling R78`", "`--captain-ruling R75`"),
                ROW_R78.replace("ruling, no run", "a try"), ROW_R78.replace("R78 (the captain's order", "R78 (a draft"),
                ROW_R78.replace("R78 (the captain's order", "R78 (the captain's ruling"),
                ROW_R78.replace("| R78 (", "| R780 ("), ROW_R78.replace("| R78 (", "| R78a ("),
                f"{ROW_R78}\n{ROW_R78}", ROW_R75, "| R59 (the captain's order |",
                "| R78 (the captain's order of X, ruling, no run) | T2 embryo | "
                + ", ".join(f"`{o}` ({f})" for f, o in TABLE) + " `--captain-ruling R78` |"):  # outside the first cell
        p.write_text(f"# ledger\n{bad}\n")
        assert not any(bt._r78_row_names(f, p) for f in NAMES), bad[:60]


def test_the_launchers_list_names_both_test_files():
    sh = (REPO / "experiments" / "external" / "n54_step.sh").read_text()
    step = sh.split("  tests)")[1].split(";;")[0]
    assert "tests/test_embryo_r78_build.py" in step and "tests/test_embryo_r78_launcher.py" in step


def test_build_embryo_runs_row_r78s_check_first(bt):
    """Row R78's check is build_embryo's first call; row R75's, row R59's and series_back's are each made only when it
    returns no file. The source is read by code and not shown."""
    src = inspect.getsource(bt.build_embryo)
    at = src.find("r78 = _r78_refusals(a)")
    calls = ("_r75_refusals(a)", "_mix_lean_refusals(a)", "_bayes_refusals(a)", "_series_back_refusals(a)",
             "_run_sorted_refusals(a)", "_copy_rows_refusals(a)", "_board(", "load_stage(")
    first = [0 <= at < src.find(call) for call in calls]
    assert first == [True] * len(calls)
    opening = src.split("\n")[1].strip().startswith("r78 = _r78_refusals(a)")  # the body's first line
    assert opening, "row R78's check is not build_embryo's first line (source not shown)"
    counts =(src.count("_r78_refusals("), src.count("_r75_refusals("), src.count("_mix_lean_refusals("),
              src.count("_series_back_refusals("))
    assert counts == (1, 1, 1, 1)


@pytest.mark.parametrize("name", NAMES)
def test_a_file_of_row_r78_skips_row_r75s_row_r59s_and_series_backs_checks(bt, env, tmp_path, monkeypatch, name):
    def never(check):
        def run(a):
            raise AssertionError(f"{check} ran on row R78's file")
        return run

    for check in ("_r75_refusals", "_mix_lean_refusals", "_series_back_refusals"):
        monkeypatch.setattr(bt, check, never(check))
    out, rec, _ = _build(bt, tmp_path, "f", *ARGS[name], name=name)
    assert out.exists() and rec["r78"]["file"] == name
    assert not {"r75", "mix_lean", "series_back"} & set(rec)


def test_the_other_checks_still_run_for_every_other_command(bt, env, tmp_path, monkeypatch):
    calls = []
    for check in ("_r75_refusals", "_mix_lean_refusals", "_series_back_refusals"):
        inner = getattr(bt, check)

        def spy(a, inner=inner, check=check):
            calls.append(check)
            return inner(a)

        monkeypatch.setattr(bt, check, spy)
    _build(bt, tmp_path, "o", *O_ARGS)
    assert calls == ["_r75_refusals", "_mix_lean_refusals", "_series_back_refusals"]


# ------------------------------------------------------------------------------------------------ routing and refusals
@pytest.mark.parametrize("name", NAMES)
def test_each_file_is_taken_with_row_r78_and_nothing_but_the_ledger_is_read(bt, env, name):
    for argv in (ARGS[name], [*ARGS[name], "--board", "T2:embryo:test_interp"], [*ARGS[name], "--seed", "1000001"]):
        a = bt.build_parser().parse_args(["embryo", *argv, "--out", "x.h5ad"])
        assert bt._r78_file(a) == name and bt._r78_refusals(a) == _values(name)
        assert bt._r75_file(a) is None and bt._r75_refusals(a) is None  # row R75's check passes it through
    a = bt.build_parser().parse_args(["embryo", *FLOOR, "--expression", "mixture", *OV, *OPTS[name].split(),
                                      "--captain-ruling", "R78", "--out", "x.h5ad"])
    assert bt._r78_refusals(a) == _values(name)  # the launcher's spelling with O's default named is the same command
    assert env == []


@pytest.mark.parametrize("name", NAMES)
@pytest.mark.parametrize("ruling", RULINGS_OTHER)
def test_each_files_string_under_any_other_ruling_or_none_is_refused(bt, env, tmp_path, name, ruling):
    argv = [*O_ARGS, *OPTS[name].split(), *([] if ruling is None else ["--captain-ruling", ruling])]
    a = bt.build_parser().parse_args(["embryo", *argv, "--out", "x.h5ad"])
    assert bt._r78_file(a) is None
    msg = _refused(bt, env, tmp_path, argv, f"--mix-lean {LEAN[name]}")
    want = ("needs --captain-ruling R78" if ruling is None
            else f"is refused with --captain-ruling {ruling}: only row R78 orders it")
    assert want in msg


@pytest.mark.parametrize("lean", ["later_11_12", "earlier"])
@pytest.mark.parametrize("argv", [
    [*FLOOR, "--expression", "series", *OV],  # SO's command with the lean
    [*FLOOR, "--expression", "depth", *OV],
    [*FLOOR, "--expression", "series_back", *OV],
    ["--method", "fgw"],
    FLOOR,
    [*FLOOR, "--shape-interp", "overlay_bayes", "--shape-interp-test-ok"],
])
@pytest.mark.parametrize("ruling", [None, "R59", "R75", "R48", "R71"])
def test_a_lean_only_row_r78_has_is_refused_on_any_command_under_any_other_ruling(bt, env, tmp_path, lean, argv,
                                                                                   ruling):
    full = [*argv, "--mix-lean", lean, *([] if ruling is None else ["--captain-ruling", ruling])]
    msg = _refused(bt, env, tmp_path, full, f"--mix-lean {lean}")
    assert "row R78" in msg


@pytest.mark.parametrize("name", NAMES)
@pytest.mark.parametrize("change", [
    ["--seed", "1000002"], ["--seed", "1000000"], ["--move-coupling", "aligned"],
    ["--row-order", "copy_parity", "--row-order-test-ok"], ["--neighbourhood", "run_sorted"], ["--herd"],
    ["--draw", "midpoint"], ["--copy-rows", "--copy-rows-test-ok"], ["--board", DRY],
    ["--expression", "series"], ["--expression", "depth"], ["--expression", "series_back"],
    ["--expression", "render_hi"], ["--expression", "blend"],
])
def test_every_other_use_of_a_files_string_under_r78_is_refused(bt, env, tmp_path, name, change):
    name_out = "T2embryo_dry.h5ad" if DRY in change else "T2embryo.h5ad"
    _refused(bt, env, tmp_path, [*ARGS[name], *change], "orders two files", name_out)


@pytest.mark.parametrize("argv", [
    [*O_ARGS, "--captain-ruling", "R78"],  # O itself
    [*OT_ARGS[:-2], "--captain-ruling", "R78"],  # row R59's file under row R78
    [*SO_ARGS, "--captain-ruling", "R78"],  # SO under row R78
    [*FLOOR, "--expression", "series_back", *OV, "--captain-ruling", "R78"],  # W under row R78
    [*FLOOR, *SOT56.split(), *OV, "--captain-ruling", "R78"],  # the struck file
    [*O_ARGS, *SOT56.split(), "--captain-ruling", "R78"],
    *[[*R75_ARGS[f][:-1], "R78"] for f in R75_ARGS],  # row R75's five strings under row R78
    [*FLOOR, "--mix-lean", "later_11_12", "--captain-ruling", "R78"],  # no form
    [*FLOOR, "--shape-interp", "overlay", "--mix-lean", "earlier", "--captain-ruling", "R78"],  # no yes flag
    [*FLOOR, "--shape-interp", "stretch_ratio", "--shape-interp-test-ok", "--mix-lean", "earlier",
     "--captain-ruling", "R78"],
    [*FLOOR, "--shape-interp", "overlay_stretch", "--shape-interp-test-ok", "--mix-lean", "later_11_12",
     "--captain-ruling", "R78"],
    [*FLOOR, "--shape-interp", "overlay_bayes", "--shape-interp-test-ok", "--mix-lean", "earlier",
     "--captain-ruling", "R78"],
    ["--method", "fgw", *OV[:2], "--mix-lean", "later_11_12", "--captain-ruling", "R78"],
    ["--method", "fgw_shape0", *OV, "--mix-lean", "earlier", "--captain-ruling", "R78"],  # positions moved
    ["--method", "fgw", "--mix-lean", "earlier", "--captain-ruling", "R78"],
])
def test_every_other_command_under_r78_is_refused(bt, env, tmp_path, argv):
    _refused(bt, env, tmp_path, argv, "orders two files")


@pytest.mark.parametrize("name", NAMES)
@pytest.mark.parametrize("which", ["none", "twice", "unpaired", "unticked", "draft", "r75_only"])
def test_a_file_whose_row_does_not_name_it_is_refused(bt, env, tmp_path, monkeypatch, name, which):
    option = OPTS[name]
    other = [f for f in NAMES if f != name][0]
    rows = {"none": (R59_ROW,), "twice": (ROW_R78, ROW_R78),
            "unpaired": (ROW_R78.replace(f"`{option}` ({name})", f"`{option}` ({other})"),),
            "unticked": (ROW_R78.replace("`--captain-ruling R78`", "--captain-ruling R78"),),
            "draft": (ROW_R78.replace("(the captain's order", "(the AI assistant's order"),),
            "r75_only": (ROW_R75,)}[which]
    p = tmp_path / "ledger" / "TRIES.md"
    p.parent.mkdir()
    p.write_text("| ID | what |\n|---|---|\n" + "".join(r + "\n" for r in rows))
    monkeypatch.setattr(bt, "TRIES_FILE", p)
    msg = _refused(bt, env, tmp_path, ARGS[name], f"row R78's file {name} is refused")
    assert f"`{option}` ({name}) and `--captain-ruling R78`" in msg


def test_the_help_names_row_r78s_files(bt):
    p = bt.build_parser()._subparsers._group_actions[0].choices["embryo"]
    helps = {a.dest: a.help for a in p._actions}
    for name, option in TABLE:
        assert f"`{option}` ({name})" in helps["mix_lean"]
    want = f"{bt.EMBRYO_R78_COMMAND.format(option='<OPTION>')}, t_mix 11/12, 1/3. Row R75's files"
    assert want in helps["mix_lean"]
    assert "with one of row R78's two option strings on O's command: R78; " in helps["captain_ruling"]
    assert helps["captain_ruling"].endswith("on O's command: R75; given once")  # row R75's words still close it
    assert p.format_help()  # the help still formats


# ------------------------------------------------------------------------------------------------ the builds
@pytest.mark.parametrize("name", NAMES)
def test_each_lean_passes_its_exact_t_mix_and_nothing_else_changes_in_os_call(bt, env, tmp_path, monkeypatch, capsys,
                                                                             stages, name):
    seen = _capture(bt, monkeypatch)
    o, rec_o, _ = _build(bt, tmp_path, "o", *O_ARGS)
    f, rec_f, text = _build(bt, tmp_path, "f", *ARGS[name], name=name)
    (a_o, kw_o, res_o), (a_f, kw_f, res_f) = seen
    tm = float(T_MIX[name])
    assert "t_mix" not in kw_o and kw_f["t_mix"] == tm and type(kw_f["t_mix"]) is float
    assert {k: v for k, v in kw_f.items() if k != "t_mix"} == kw_o  # form, positions, t_pos, seed, ... as O's
    assert a_f[4] == a_o[4] == 0.5 and a_f[5:] == a_o[5:]  # t = 1/2 still sets the moves; O's count and seed
    for x, y in zip(a_f[:4], a_o[:4], strict=True):
        np.testing.assert_array_equal(x, y)
    n = len(res_o.rows)
    n_b = int(round((1.0 - tm) * n))  # interpolate_clouds' own count (its rounding, not the exact fraction's)
    assert len(res_f.rows) == n and int((res_f.stage == 0).sum()) == n_b and int((res_f.stage == 1).sum()) == n - n_b
    assert 0 < n_b < n  # both brackets' tissue kept
    assert res_f.record["t_mix"] == tm and res_f.record["t"] == 0.5 and res_f.record["t_pos"] == 0.0
    w, n75 = res_f.pool_weights, stages["E7.5.h5ad"].n_obs  # the overlay's pool weights follow t_mix
    assert w[:n75].sum() == pytest.approx(1 - tm) and w[n75:].sum() == pytest.approx(tm)
    Xo, Co, no, vo = _file(o)
    Xf, Cf, nf, vf = _file(f)
    assert len(Xf) == len(Xo) == n and vf == vo and nf == no  # O's count, genes and names
    for C in (Co, Cf):  # the size rule, after the form
        rms = float(np.sqrt(((C.astype(np.float64) - C.mean(0)) ** 2).sum(1).mean()))
        assert rms == pytest.approx(6.0, rel=1e-3)
    assert not np.array_equal(Cf, Co) and not np.array_equal(Xf, Xo)
    # every row's expression is its bracket's released cell (the mixture): depth moves only through the shares
    Eb = stages["E7.5.h5ad"][:, list(vo)].X.toarray()
    Ea = stages["E8.0.h5ad"][:, list(vo)].X.toarray()
    np.testing.assert_array_equal(Xf[res_f.stage == 0], Eb[res_f.rows[res_f.stage == 0]])
    np.testing.assert_array_equal(Xf[res_f.stage == 1], Ea[res_f.rows[res_f.stage == 1]])
    r, ro = rec_f["record"], rec_o["record"]
    assert r["t_mix"] == tm and ro["t_mix"] == 0.5 and set(r) == set(ro)
    assert {k for k in r if r[k] != ro[k]} <= {"t_mix", "scale_from_canonical", "shape_interp_thin_sd"}
    assert r["rotation_deg"] == ro["rotation_deg"] and not set(bt.COUNT_FIELDS) & set(r)
    shown = text + capsys.readouterr().out
    assert not _shows(shown, SECRET_MAX, SECRET_TRUTH, n_b, n - n_b)
    f2, _, _ = _build(bt, tmp_path, "f2", *ARGS[name], name=name)  # built again: the same bytes
    assert f2.read_bytes() == f.read_bytes()


def test_the_two_files_and_row_r75s_leans_are_all_different_files(bt, env, tmp_path):
    built = {name: _build(bt, tmp_path, name.removesuffix(".h5ad"), *ARGS[name], name=name)[0].read_bytes()
             for name in NAMES}
    for name in ("T2embryo_s2x_ovt56.h5ad", "T2embryo_s2x_ovt512.h5ad"):
        built[name] = _build(bt, tmp_path, name.removesuffix(".h5ad"), *R75_ARGS[name], name=name)[0].read_bytes()
    built["O"] = _build(bt, tmp_path, "o", *O_ARGS)[0].read_bytes()
    assert len({hashlib.sha256(b).hexdigest() for b in built.values()}) == len(built)


@pytest.mark.parametrize("name", NAMES)
def test_each_files_record(bt, env, tmp_path, capsys, name):
    _, rec, text = _build(bt, tmp_path, "f", *ARGS[name], name=name)
    blk = rec["r78"]
    t = T_MIX[name]
    what = ("the captain's named exception on the embryo's overlay form and occupancy term for this file"
            + (", and the lean-early mixing line declined on 2026-11-06 (e) reopened for this file alone"
               if name == "T2embryo_s2x_ovt13.h5ad" else ""))
    words = {"T2embryo_s2x_ovt1112.h5ad": "Yes to T2embryo_s2x_ovt1112.h5ad under R22, sent at once in one ruling row, "
                                          "with no gate",
             "T2embryo_s2x_ovt13.h5ad": "Yes to T2embryo_s2x_ovt13.h5ad under R22 in the same row, sent at once, with "
                                        "no gate"}[name]
    assert words in FIRST_R78  # the captain's words as the row records them
    assert blk == {"option": OPTS[name], "file": name, "captain_ruling": "R78",
                   "order": bt.EMBRYO_R78_ORDER.format(option=OPTS[name], file=name, words=words, what=what),
                   "command": "embryo " + " ".join(ARGS[name]), "o_command": O_CMD,
                   "expression": "mixture", "mix_lean": LEAN[name], "t_mix": float(t),
                   "t_mix_exact": f"{t.numerator}/{t.denominator}", "rule": bt.EMBRYO_R78_RULES[name]}
    assert all(not isinstance(v, (dict, list)) for v in blk.values())
    assert f"`{OPTS[name]}` ({name}) and `--captain-ruling R78`" in blk["order"]
    assert "the E8.0 twin (t_mix 1) still declined" in blk["order"] and "decision log 2026-11-26 (h)" in blk["order"]
    assert blk["rule"].startswith(f"t_mix = {t.numerator}/{t.denominator} in vec.t2_build.")
    assert rec["n_cells"] == "max_cells" and not {"mix_lean", "series_back", "bayes", "r75"} & set(rec)
    assert not _shows(text + capsys.readouterr().out, SECRET_MAX, SECRET_TRUTH)


# ------------------------------------------------------------------------------------------------ the base checker
def _check(chk, tmp_path: Path, name: str, *extra: str) -> tuple[int, Path]:
    out = tmp_path / "half" / f"{name.removesuffix('.h5ad')}_half.h5ad"
    return chk.main(["half", "embryo", *ARGS[name], *extra, "--out", str(out)]), out


@pytest.mark.parametrize("name", NAMES)
def test_the_checker_gives_os_bytes_from_each_files_own_path(bt, chk, env, tmp_path, monkeypatch, name):
    seen = _capture(bt, monkeypatch)
    k, _, _ = _build(bt, tmp_path, "known", *O_ARGS)
    code, out = _check(chk, tmp_path, name)
    assert code == 0 and out.read_bytes() == k.read_bytes()
    assert seen[-1][1]["t_mix"] == 0.5  # the file's path, its t_mix put back at O's 1/2
    assert out.with_name(out.name + ".stamp.json").exists()  # stamped as any build: the launcher reads its fields
    f, _, _ = _build(bt, tmp_path, "file", *ARGS[name], name=name)
    assert f.read_bytes() != k.read_bytes()  # the file itself is not O


def test_the_checkers_mode_is_the_launchers(chk):
    sh = (REPO / "experiments" / "builds" / "embryo_r78_line.sh").read_text()
    assert chk.MODES == ("half",) and chk.HALF == Fraction(1, 2)
    assert 'build_base "$STEP" "${NAME}_half" half $ARGS || finish $?' in sh
    assert set(chk.EXPRESSIONS) == {"series_rendered_expression", "depth_matched_expression",
                                    "depth_rendered_expression", "series_back_expression", "blend_expression"}


@pytest.mark.parametrize("argv, match", [
    (O_ARGS, "checks row R78's files only"),
    (OT_ARGS, "checks row R78's files only"),
    (R75_ARGS["T2embryo_s2x_ovt56.h5ad"], "checks row R78's files only"),
    ([*ARGS["T2embryo_s2x_ovt13.h5ad"], "--seed", "1000002"], "orders two files"),
    ([*O_ARGS, "--mix-lean", "earlier", "--captain-ruling", "R59"], "only row R78 orders it"),
])
def test_the_checker_refuses_a_command_it_does_not_check(bt, chk, env, tmp_path, argv, match):
    out = tmp_path / "c" / "T2embryo_check.h5ad"
    with pytest.raises(SystemExit, match=match):
        chk.main(["half", "embryo", *argv, "--out", str(out)])
    assert not out.exists() and env == []


@pytest.mark.parametrize("argv", [[], ["half"], ["base", "embryo"], ["half", "t1"], ["Half", "embryo"],
                                  ["mixture", "embryo"], ["half", "heart"]])
def test_the_checkers_usage(chk, argv):
    with pytest.raises(SystemExit, match="usage: embryo_r78_check.py"):
        chk.main(argv)


def test_the_checker_fails_when_the_route_interpolates_more_than_once(bt, chk, env, tmp_path, monkeypatch):
    """A build that makes one interpolation more than the route's (here a small synthetic one before build_embryo's
    own) is not the route the checker names: it exits 1 though the file's own interpolation is as named."""
    code, _ = _check(chk, tmp_path / "once", "T2embryo_s2x_ovt1112.h5ad")
    assert code == 0
    inner_build = bt.build_embryo

    def twice(a):
        rng = np.random.default_rng(0)
        C, E = rng.normal(size=(60, 3)), rng.gamma(2.0, 1.0, size=(60, 12))
        bt.t2_build.interpolate_clouds(C, C + 1.0, E, E, 0.5, 20, 1, target_rms=1.0)
        return inner_build(a)

    monkeypatch.setattr(bt, "build_embryo", twice)
    code, _ = _check(chk, tmp_path / "twice", "T2embryo_s2x_ovt1112.h5ad")
    assert code == 1


@pytest.mark.parametrize("fn", ["series_rendered_expression", "depth_matched_expression", "blend_expression"])
def test_the_checker_fails_when_an_expression_function_runs(bt, chk, env, tmp_path, monkeypatch, fn):
    inner_build = bt.build_embryo

    def with_expression(a):
        try:
            getattr(bt.t2_build, fn)(None, None, None, None, 0.5, [1, 7])  # a call the route never makes
        except Exception:
            pass
        return inner_build(a)

    monkeypatch.setattr(bt, "build_embryo", with_expression)
    code, _ = _check(chk, tmp_path, "T2embryo_s2x_ovt13.h5ad")
    assert code == 1


def test_the_checker_leaves_the_builder_as_it_found_it(bt, chk, env, tmp_path):
    before = (bt._r78_refusals, bt.t2_build.interpolate_clouds,
              *(getattr(bt.t2_build, fn) for fn in chk.EXPRESSIONS))
    code, _ = _check(chk, tmp_path, "T2embryo_s2x_ovt13.h5ad")
    after = (bt._r78_refusals, bt.t2_build.interpolate_clouds, *(getattr(bt.t2_build, fn) for fn in chk.EXPRESSIONS))
    assert code == 0 and after == before


# ------------------------------------------------------------------------------------------------ every other command
def _module(name: str, rel: str) -> types.ModuleType:
    """A module made from BASE's text of `rel`, registered under `name` (dataclasses look their module up there)."""
    mod = types.ModuleType(name)
    mod.__file__ = str(REPO / rel)
    sys.modules[name] = mod
    exec(compile(git(REPO, "show", f"{BASE}:{rel}") + "\n", str(REPO / rel), "exec"), mod.__dict__)
    return mod


@pytest.fixture(scope="module")
def base_bt(bt):
    """BASE's builder: row R78's route changes build_test.py only, and every module is this checkout's in both
    builders, so a difference below is the route's."""
    mod = _module("_build_test_at_base_r78", "experiments/builds/build_test.py")
    assert not hasattr(mod, "EMBRYO_R78_TABLE") and hasattr(mod, "EMBRYO_R75_TABLE")
    return mod


# name: (argv, file name, what the builder of BASE does with it: "built", "refused" or None for either)
COMMANDS = {
    "O": (O_ARGS, "T2embryo.h5ad", "built"),
    "O_launcher": ([*FLOOR, "--expression", "mixture", *OV], "T2embryo.h5ad", "built"),
    "Ot23": (OT_ARGS, "T2embryo.h5ad", "built"),
    "W": (W_ARGS, "T2embryo.h5ad", "built"),
    "SO": (SO_ARGS, "T2embryo.h5ad", "built"),
    "O_hi": ([*FLOOR, "--expression", "render_hi", *OV], "T2embryo.h5ad", "built"),
    **{f"R75_{f.removesuffix('.h5ad')}": (R75_ARGS[f], f, "built") for f in R75_ARGS},  # row R75's five files
    "P": (FLOOR, "T2embryo.h5ad", "built"),
    "P_depth": ([*FLOOR, "--expression", "depth"], "T2embryo.h5ad", "built"),
    "s1": (["--method", "fgw"], "T2embryo.h5ad", "built"),
    "s2": (["--method", "fgw_shape0"], "T2embryo.h5ad", "built"),
    "series": (["--method", "fgw", "--expression", "series"], "T2embryo.h5ad", "built"),
    "seed": ([*O_ARGS, "--seed", "1000003"], "T2embryo.h5ad", "built"),
    "dry_O": ([*FLOOR, "--shape-interp", "overlay", "--board", DRY], "T2embryo_dry.h5ad", "built"),
    "RS": ([*O_ARGS, "--neighbourhood", "run_sorted", "--captain-ruling", "R60"], "T2embryo.h5ad", None),
    "BF": ([*FLOOR, "--shape-interp", "overlay_bayes", "--shape-interp-test-ok", "--captain-ruling", "R71"],
           "T2embryo_s2x_ovbf.h5ad", None),
    "copy_parity": ([*FLOOR, "--row-order", "copy_parity", "--row-order-test-ok"], "T2embryo.h5ad", None),
    "copy_rows": ([*FLOOR, "--copy-rows", "--copy-rows-test-ok"], "T2embryo.h5ad", None),
    # the refusals as they were, rows R48's, R59's and R75's among them
    "W_no_ruling": ([*FLOOR, "--expression", "series_back", *OV], "T2embryo.h5ad", "refused"),
    "Ot23_no_ruling": ([*O_ARGS, "--mix-lean", "later"], "T2embryo.h5ad", "refused"),
    "Ot23_R75": ([*O_ARGS, "--mix-lean", "later", "--captain-ruling", "R75"], "T2embryo.h5ad", "refused"),
    "Ot23_seed": ([*OT_ARGS, "--seed", "1000002"], "T2embryo.h5ad", "refused"),
    "R75_lean_no_ruling": ([*O_ARGS, "--mix-lean", "later_half"], "T2embryo.h5ad", "refused"),
    "R75_lean_R59": ([*O_ARGS, "--mix-lean", "earlier_half", "--captain-ruling", "R59"], "T2embryo.h5ad", "refused"),
    "R75_seed": ([*R75_ARGS["T2embryo_s2x_ovt56.h5ad"], "--seed", "1000002"], "T2embryo.h5ad", "refused"),
    "R75_depth_no_ruling": ([*FLOOR, "--expression", "depth", *OV], "T2embryo.h5ad", "refused"),
    "sot56_R75": ([*FLOOR, *SOT56.split(), *OV, "--captain-ruling", "R75"], "T2embryo.h5ad", "refused"),
    "later2_half_R78": ([*O_ARGS, "--mix-lean", "later2_half", "--captain-ruling", "R78"], "T2embryo.h5ad",
                        "refused"),  # the design's name: the parser's refusal at BASE and now
    "ruling_R48_alone": ([*O_ARGS, "--captain-ruling", "R48"], "T2embryo.h5ad", "refused"),
    "ruling_R59_alone": ([*O_ARGS, "--captain-ruling", "R59"], "T2embryo.h5ad", "refused"),
    "ruling_R75_alone": ([*O_ARGS, "--captain-ruling", "R75"], "T2embryo.h5ad", "refused"),
    "ruling_R77_alone": ([*O_ARGS, "--captain-ruling", "R77"], "T2embryo.h5ad", "refused"),
    "ruling_R79_alone": ([*O_ARGS, "--captain-ruling", "R79"], "T2embryo.h5ad", "refused"),
    "form_needs_ok": ([*FLOOR, "--shape-interp", "overlay"], "T2embryo.h5ad", "refused"),
    "seed_low": ([*O_ARGS, "--seed", "5"], "T2embryo.h5ad", "refused"),
    "dry_name": ([*FLOOR, "--board", DRY], "T2embryo.h5ad", "refused"),
}


def _outcome(mod, tmp: Path, stages, which: str, argv: list[str], name: str):
    """What one builder does with one command: the h5ad's sha256 and its record's text, or ("refused", its message)."""
    with pytest.MonkeyPatch.context() as mp:
        install(mod, mp.setattr, mp.setitem, tmp / "shared", [], stages, ledger_rows=LEDGER)
        out = tmp / which / name
        try:
            mod.main(["embryo", *argv, "--out", str(out)])
        except SystemExit as e:
            return "refused", str(e.code).replace(str(out), "<out>")
        except Exception as e:  # a step's own refusal (vec.run_sort's, say): the same exception from both builders
            return "refused", f"{type(e).__name__}: {e}"
        return hashlib.sha256(out.read_bytes()).hexdigest(), out.with_suffix(".build.json").read_text()


@pytest.mark.parametrize("key", list(COMMANDS))
def test_every_other_command_does_what_base_does(bt, base_bt, own_repo, tmp_path, stages, key):
    argv, name, kind = COMMANDS[key]
    want = _outcome(base_bt, tmp_path, stages, "base", argv, name)
    got = _outcome(bt, tmp_path, stages, "new", argv, name)
    base_did = "refused" if want[0] == "refused" else "built"
    assert kind is None or base_did == kind, (base_did, kind)
    assert got[0] == want[0] and got[1] == want[1]  # the same bytes and record, or the same refusal


@pytest.mark.parametrize("name", NAMES)
def test_row_r78s_two_commands_were_refused_at_base_and_are_built_now(bt, base_bt, own_repo, tmp_path, stages, name):
    want = _outcome(base_bt, tmp_path, stages, "base", ARGS[name], name)
    got = _outcome(bt, tmp_path, stages, "new", ARGS[name], name)
    assert want == ("refused", "2")  # BASE's parser knew neither lean
    assert got[0] != "refused" and json.loads(got[1])["r78"]["file"] == name


@pytest.mark.parametrize("name", NAMES)
@pytest.mark.parametrize("ruling", [None, "R59", "R75", "R48"])
def test_row_r78s_leans_elsewhere_were_refused_at_base_and_are_refused_now_by_row_r78(bt, base_bt, own_repo, tmp_path,
                                                                                      stages, name, ruling):
    argv = [*O_ARGS, *OPTS[name].split(), *([] if ruling is None else ["--captain-ruling", ruling])]
    want = _outcome(base_bt, tmp_path, stages, "base", argv, "T2embryo.h5ad")
    got = _outcome(bt, tmp_path, stages, "new", argv, "T2embryo.h5ad")
    assert want == ("refused", "2")  # BASE's parser
    assert got[0] == "refused" and "row R78" in got[1] and f"--mix-lean {LEAN[name]}" in got[1]


@pytest.mark.parametrize("argv", [
    [*O_ARGS, "--captain-ruling", "R78"],
    [*R75_ARGS["T2embryo_s2x_ovt56.h5ad"][:-1], "R78"],
    [*FLOOR, *SOT56.split(), *OV, "--captain-ruling", "R78"],
    [*FLOOR, "--expression", "depth", *OV, "--captain-ruling", "R78"],
])
def test_any_other_command_under_r78_was_refused_at_base_and_is_refused_now(bt, base_bt, own_repo, tmp_path, stages,
                                                                            argv):
    want = _outcome(base_bt, tmp_path, stages, "base", argv, "T2embryo.h5ad")
    got = _outcome(bt, tmp_path, stages, "new", argv, "T2embryo.h5ad")
    assert want[0] == "refused" and got[0] == "refused" and "orders two files" in got[1]
