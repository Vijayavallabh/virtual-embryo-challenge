"""Row R75's five files (the captain's order of 2026-11-26: the embryo's refinement round on O): the builder's one
table EMBRYO_R75_TABLE (experiments/builds/build_test.py embryo) and the base checker experiments/builds/
embryo_r75_check.py, on synthetic stages and a synthetic test board whose limits are never shown (the synthetic world
of tests/test_mix_lean_build.py). Nothing is read from data/ or boards/; a synthetic TRIES.md stands in for the ledger
(two tests ask, by code and printing nothing, whether the repository's row R75 names the table as the builder needs it).

Pinned: the table is the row's five pairs in its order, each name O's with its expression's mark and its lean's
fraction, and no string or name of it is an earlier row's or the struck lean's; a lean passes its exact t_mix's nearest
double to vec.t2_build.interpolate_clouds and nothing else changes in O's call, the count of E7.5 rows being
interpolate_clouds' own int(round((1 - t_mix) n)); the series file is SO's expression on row R59's rows (Ot23's
coordinates bit for bit, its E8.0 rows' expression Ot23's); the depth file keeps O's rows and coordinates bit for bit
and thins the deeper bracket to DEPTH_RULE's depth; each file's record gains the r75 block and shows no count; each is
taken only with --captain-ruling R75 on O's command with exactly its string, at the default seed on the test board,
while row R75's first cell names the string paired with the file and the ruling, every refusal before anything is
read; the checker gives O's, SO's and Ot23's bytes from the files' own paths with their one change put back and fails
when the route is not as named; and every other embryo command gives what the builder of BASE (main before row R75's
route) gives, run in this same process from git's copy of build_test.py: the same h5ad bytes and record, or the same
refusal, but for the commands row R75 itself changes (its five files, now built, and the depth expression with the
overlay form and no ruling, now refused). Built through experiments/builds/build_test.py."""

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
from test_mix_lean_build import R48_ROW, R59_ROW, R60_RS_ROW, SECRET_MAX, SECRET_TRUTH, _file, _shows, _stage, install

REPO = Path(__file__).resolve().parents[1]
BASE = "e0a106299712abbfc992acd8caaf493916a83241"  # main with rows R75 and R76 placed, before row R75's route
FLOOR = ["--method", "fgw_shape0", "--positions", "floor_exact"]
OV = ["--shape-interp", "overlay", "--shape-interp-test-ok"]
O_ARGS = [*FLOOR, *OV]
O_CMD = "embryo " + " ".join(O_ARGS)
OT_ARGS = [*O_ARGS, "--mix-lean", "later", "--captain-ruling", "R59"]
SO_ARGS = [*FLOOR, "--expression", "series", *OV]
W_ARGS = [*FLOOR, "--expression", "series_back", *OV, "--captain-ruling", "R48"]
DRY = "T2:embryo:val_interp"
TABLE = [("T2embryo_s2x_ovt712.h5ad", "--mix-lean later_half"),
         ("T2embryo_s2x_ovt56.h5ad", "--mix-lean later2"),
         ("T2embryo_s2x_ovt512.h5ad", "--mix-lean earlier_half"),
         ("T2embryo_s2xs_ovt23.h5ad", "--expression series --mix-lean later"),
         ("T2embryo_s2x_ovd.h5ad", "--expression depth")]
OPTS = dict(TABLE)
NAMES = list(OPTS)
ARGS = {f: [*O_ARGS, *o.split(), "--captain-ruling", "R75"] for f, o in TABLE}
T_MIX = {"T2embryo_s2x_ovt712.h5ad": Fraction(7, 12), "T2embryo_s2x_ovt56.h5ad": Fraction(5, 6),
         "T2embryo_s2x_ovt512.h5ad": Fraction(5, 12), "T2embryo_s2xs_ovt23.h5ad": Fraction(2, 3),
         "T2embryo_s2x_ovd.h5ad": None}
EXPRESSION = {"T2embryo_s2x_ovt712.h5ad": "mixture", "T2embryo_s2x_ovt56.h5ad": "mixture",
              "T2embryo_s2x_ovt512.h5ad": "mixture", "T2embryo_s2xs_ovt23.h5ad": "series",
              "T2embryo_s2x_ovd.h5ad": "depth"}
LEANS = [f for f in NAMES if EXPRESSION[f] == "mixture"]  # the three files that change t_mix alone
SERIES, DEPTH = "T2embryo_s2xs_ovt23.h5ad", "T2embryo_s2x_ovd.h5ad"
# row R75's first cell as the repository's ledger holds it (tests below check the copy against it by code)
FIRST_R75 = (
    "R75 (the captain's order of 2026-11-26, ruling, no run: \"Yes to T2embryo_s2x_ovt712.h5ad, "
    "T2embryo_s2x_ovt56.h5ad, T2embryo_s2x_ovt512.h5ad, T2embryo_s2xs_ovt23.h5ad and T2embryo_s2x_ovd.h5ad under R22, "
    "sent at once in one "
    "ruling row, with no gate: my named exception on the embryo's overlay form and occupancy term for the first four; "
    "the lean-early mixing line declined on 2026-11-06 (e) reopened for T2embryo_s2x_ovt512.h5ad alone; the embryo's "
    "depth line reopened for T2embryo_s2x_ovd.h5ad alone, beside O and never in its place.\" (his pick \"All five\" in "
    "a dialog with no recommended option, the AI assistant being recused on the embryo's occupancy term and shape "
    "forms; the option recorded the words above as the captain's); the embryo's refinement round on O "
    "(EMBRYO_REF3_design.md with EMBRYO_REF3_challenge.md's fixes; the design's fourth file struck as already struck "
    "in an earlier round); the embryo builder is to take on `embryo --method fgw_shape0 --positions floor_exact "
    "--shape-interp overlay --shape-interp-test-ok` exactly these option strings, each only with `--captain-ruling "
    "R75`: `--mix-lean later_half` (T2embryo_s2x_ovt712.h5ad), `--mix-lean later2` (T2embryo_s2x_ovt56.h5ad), "
    "`--mix-lean earlier_half` (T2embryo_s2x_ovt512.h5ad), `--expression series --mix-lean later` "
    "(T2embryo_s2xs_ovt23.h5ad), `--expression depth` (T2embryo_s2x_ovd.h5ad); t_mix 7/12, 5/6, 5/12 and 2/3; "
    "decision log 2026-11-26 (a); not a try)")
ROW_R75 = f"| {FIRST_R75} | T2 embryo | synthetic |"
LEDGER = (R48_ROW, R59_ROW, R60_RS_ROW, ROW_R75)
RULINGS_OTHER = [None, "R59", "R48", "R71", "R60", "R73", "R74", "R76", "R750", "r75", "R75a", "R7", "E58"]


def _pairs(cell: str) -> list[tuple[str, str]]:
    return [(f, o) for o, f in re.findall(r"`([^`]*)` \((T2embryo_[A-Za-z0-9_]+\.h5ad)\)", cell)]


def _values(name: str) -> dict:
    option = OPTS[name]
    lean = option.split()[-1] if "--mix-lean" in option else None
    return {"file": name, "option": option, "expression": EXPRESSION[name], "mix_lean": lean, "t_mix": T_MIX[name]}


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
    import embryo_r75_check

    assert embryo_r75_check.bt is bt  # the checker drives this same builder module
    return embryo_r75_check


@pytest.fixture
def env(bt, monkeypatch, own_repo, tmp_path, stages):
    """The synthetic stages, boards, size JSON and ledger (rows R48, R59, R60's synthetic ones and row R75's first
    cell) installed in the builder; returns the list of stage reads."""
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
def test_the_one_table_is_the_rows_five_pairs_in_its_order(bt):
    assert [tuple(p) for p in bt.EMBRYO_R75_TABLE] == TABLE == _pairs(FIRST_R75)
    assert re.findall(r"`([^`]*)`", FIRST_R75) == [O_CMD, "--captain-ruling R75", *OPTS.values()]
    assert bt.EMBRYO_R75_ROW == "R75" and list(bt.EMBRYO_R75_FILES) == NAMES
    for name, (option, expression, lean, t_mix) in bt.EMBRYO_R75_FILES.items():
        assert (option, expression, lean, t_mix) == tuple(_values(name)[k] for k in ("option", "expression",
                                                                                       "mix_lean", "t_mix"))
        assert t_mix is None or type(t_mix) is Fraction
        assert bt.EMBRYO_R75_COMMAND.format(option=option) == "embryo " + " ".join(ARGS[name])
        assert bt._r75_values(name, option) == (expression, lean, t_mix)
    assert bt.EMBRYO_R75_LEANS == {"later_half": Fraction(7, 12), "later2": Fraction(5, 6),
                                   "earlier_half": Fraction(5, 12), "later": Fraction(2, 3)}
    # the row's t_mix values in its order; later stays row R59's own float
    assert [T_MIX[f] for f in NAMES[:4]] == [Fraction(7, 12), Fraction(5, 6), Fraction(5, 12), Fraction(2, 3)]
    assert "t_mix 7/12, 5/6, 5/12 and 2/3;" in FIRST_R75
    assert float(bt.EMBRYO_R75_LEANS["later"]) == bt.EMBRYO_OT23_T_MIX == 2.0 / 3.0
    assert [float(T_MIX[f]) for f in LEANS] == [7 / 12, 5 / 6, 5 / 12]


def test_row_r59s_lean_stays_and_the_parser_takes_row_r75s_three_more(bt, capsys):
    assert bt.MIX_LEANS == ("later",) and bt.EMBRYO_R75_NEW_LEANS == ("later_half", "later2", "earlier_half")
    action = {a.dest: a for a in bt.build_parser()._subparsers._group_actions[0].choices["embryo"]._actions}
    assert action["mix_lean"].choices == ("later", "later_half", "later2", "earlier_half")
    assert action["mix_lean"].default is None and action["captain_ruling"].default is None
    for value in ("earlier", "Later_half", "later_halves", "7/12", "5/6", "later3", ""):  # earlier: the struck file
        with pytest.raises(SystemExit) as e:
            bt.build_parser().parse_args(["embryo", *O_ARGS, "--mix-lean", value, "--captain-ruling", "R75",
                                          "--out", "x.h5ad"])
        assert e.value.code == 2 and "invalid choice" in capsys.readouterr().err


def test_no_string_or_name_of_row_r75_is_an_earlier_rows_or_the_struck_leans(bt):
    earlier_files = {bt.EMBRYO_OT23_FILE, bt.EMBRYO_OVSB_FILE, bt.EMBRYO_RS_FILE, *bt.EMBRYO_BF_COMMANDS,
                     "T2embryo_s2x_ov.h5ad", "T2embryo_s2xs_ov.h5ad", "T2embryo_s2x_ovhi.h5ad", "T2embryo_s2x.h5ad",
                     "T2embryo_s2xs.h5ad", "T2embryo_s2x_ovt13.h5ad"}
    assert not set(NAMES) & earlier_files
    commands = {bt.EMBRYO_R75_COMMAND.format(option=o) for o in OPTS.values()}
    assert not commands & {bt.EMBRYO_OT23_COMMAND, bt.EMBRYO_OVSB_COMMAND, *bt.EMBRYO_BF_COMMANDS.values(),
                           bt.EMBRYO_RS_COMMAND.format(row="R75")}
    assert not {o for o in OPTS.values()} & {"--mix-lean later", "--expression series", "--expression series_back",
                                            "--expression render_hi", "--mix-lean earlier"}
    assert "earlier" not in bt.EMBRYO_R75_LEANS  # the design's fourth file, struck as in an earlier round


@pytest.mark.parametrize("name, option", [
    ("T2embryo_s2x_ovt13.h5ad", "--mix-lean earlier"),  # the struck lean: no value for it
    ("T2embryo_s2x_ovt712.h5ad", "--expression mixture --mix-lean later_half"),  # O's default named: not the string
    ("T2embryo_s2xs_ovt23.h5ad", "--mix-lean later --expression series"),  # the lean before the expression
    ("T2embryo_s2x_ovt712.h5ad", "--mix-lean  later_half"),  # two spaces
    ("T2embryo_s2x_ovt712.h5ad", " --mix-lean later_half"),
    ("T2embryo_s2x_ovt712.h5ad", "--mix-lean later_half "),
    ("T2embryo_s2x_ovt712.h5ad", "--mix-lean Later_half"),
    ("T2embryo_s2x_ovt712.h5ad", "--mix-lean 7/12"),
    ("T2embryo_s2x_ovt712.h5ad", "--mix-lean=later_half"),
    ("T2embryo_s2x_ovt712.h5ad", "--mix-lean later_half --mix-lean later_half"),
    ("T2embryo_s2x_ovt712.h5ad", "--mix-lean later_half --herd"),
    ("T2embryo_s2x_ovt712.h5ad", "--mix-lean"),
    ("T2embryo_s2x_ovt56.h5ad", "--mix-lean later_half"),  # another file's string
    ("T2embryo_s2x_ovt712.h5ad", "--mix-lean later2"),
    ("T2embryo_s2x_ovt1012.h5ad", "--mix-lean later2"),  # 5/6 named unreduced
    ("T2embryo_s2x_ovt0.83.h5ad", "--mix-lean later2"),
    ("T2embryo_s2xs_ovt512.h5ad", "--mix-lean earlier_half"),  # a series mark without series
    ("T2embryo_s2x_ovt23.h5ad", "--expression series --mix-lean later"),  # row R59's name on the series string
    ("T2embryo_s2x_ovt23s.h5ad", "--expression series --mix-lean later"),
    ("T2embryo_s2xs_ov.h5ad", "--expression series"),  # SO's own string: no lean
    ("T2embryo_s2x_ovt23.h5ad", "--mix-lean later"),  # well formed, but row R59's file: not in the table (below)
    ("T2embryo_s2x_ovd.h5ad", "--expression depth --mix-lean later"),  # depth with a lean
    ("T2embryo_s2x_ovt23d.h5ad", "--expression depth --mix-lean later"),
    ("T2embryo_s2x_ovsb.h5ad", "--expression series_back"),  # row R48's
    ("T2embryo_s2x_ovsbt23.h5ad", "--expression series_back --mix-lean later"),  # W with a lean
    ("T2embryo_s2x_ovhi.h5ad", "--expression render_hi"),
    ("T2embryo_s2x_ov.h5ad", ""),  # O's own: no string
    ("T2embryo_s2x_ovd.h5ad", "--expression  depth"),
    ("T2embryo_s2x_ovd.h5ad", "--expression=depth"),
    ("T2embryo_s2x_ovd.h5ad", "--expression Depth"),
    ("T2embryo_s2x_ovd", "--expression depth"),  # no suffix
    ("t2embryo_s2x_ovd.h5ad", "--expression depth"),
    ("T2embryo_s2x_ovD.h5ad", "--expression depth"),
    ("T2embryo_s2x_dov.h5ad", "--expression depth"),
])
def test_a_table_entry_off_the_rows_form_is_not_read(bt, name, option):
    if (name, option) == ("T2embryo_s2x_ovt23.h5ad", "--mix-lean later"):
        assert bt._r75_values(name, option) == ("mixture", "later", Fraction(2, 3))  # the form reads; the table has
        assert (name, option) not in TABLE and name == bt.EMBRYO_OT23_FILE  # no such entry: row R59's own file
        return
    with pytest.raises(ValueError):
        bt._r75_values(name, option)


def test_the_real_rows_first_cell_pairs_the_tables_strings(bt):
    """Row R75's first cell in the repository's ledger, read by code and printed nowhere: exactly one row R75, its first
    cell pairing the table's five strings with their files in the table's order, naming in backticks what the copy
    names, word for word the copy, and passing the builder's own check for each file."""
    text = (REPO / "experiments" / "harness" / "TRIES.md").read_text(encoding="utf-8")
    cells = [ln[2:].split("|", 1)[0].strip() for ln in text.splitlines() if ln.startswith("| R75 ")]
    ok = (len(cells) == 1 and _pairs(cells[0]) == TABLE
          and re.findall(r"`([^`]*)`", cells[0]) == re.findall(r"`([^`]*)`", FIRST_R75)
          and cells[0].startswith("R75 (the captain's order") and "ruling, no run" in cells[0]
          and all(bt._r75_row_names(f) for f in NAMES))
    assert ok, "the repository's row R75 does not name the table as the copy does (nothing of it printed)"
    assert cells[0] == FIRST_R75, "the copy is not the row's first cell (nothing of it printed)"


def test_the_ruling_check_takes_only_a_single_row_r75_pairing_the_whole_string_with_its_file(bt, tmp_path):
    p = tmp_path / "TRIES.md"
    p.write_text(f"# ledger\n{ROW_R75}\n")
    assert all(bt._r75_row_names(f, p) for f in NAMES)
    for name, option in TABLE:
        other = [f for f in NAMES if f != name][0]
        for bad in (ROW_R75.replace(f"`{option}`", option),  # not backticked
                    ROW_R75.replace(f"`{option}`", f"`{option} --herd`"),  # inside a longer string
                    ROW_R75.replace(f"`{option}` ({name})", f"`{option}` ({other})"),  # paired with another file
                    ROW_R75.replace(f"`{option}` ({name})", f"`{option}`  ({name})")):  # not followed by " (<file>)"
            p.write_text(f"# ledger\n{bad}\n")
            assert [bt._r75_row_names(f, p) for f in NAMES] == [f != name for f in NAMES], (name, bad[-80:])
    for bad in (ROW_R75.replace("`--captain-ruling R75`", "--captain-ruling R75"),
                ROW_R75.replace("`--captain-ruling R75`", "`--captain-ruling R59`"),
                ROW_R75.replace("ruling, no run", "a try"), ROW_R75.replace("R75 (the captain's order", "R75 (a draft"),
                ROW_R75.replace("R75 (the captain's order", "R75 (the captain's ruling"),
                ROW_R75.replace("| R75 (", "| R750 ("), ROW_R75.replace("| R75 (", "| R75a ("),
                f"{ROW_R75}\n{ROW_R75}", "| R59 (the captain's order |",
                "| R75 (the captain's order of X, ruling, no run) | T2 embryo | "
                + ", ".join(f"`{o}` ({f})" for f, o in TABLE) + " `--captain-ruling R75` |"):  # outside the first cell
        p.write_text(f"# ledger\n{bad}\n")
        assert not any(bt._r75_row_names(f, p) for f in NAMES), bad[:60]


def test_the_launchers_list_names_both_test_files():
    sh = (REPO / "experiments" / "external" / "n54_step.sh").read_text()
    step = sh.split("  tests)")[1].split(";;")[0]
    assert "tests/test_embryo_r75_build.py" in step and "tests/test_embryo_r75_launcher.py" in step


def test_build_embryo_runs_row_r75s_refusals_first(bt):
    src = inspect.getsource(bt.build_embryo)
    at = src.index("r75 = _r75_refusals(a)")
    for call in ("_mix_lean_refusals(a)", "_bayes_refusals(a)", "_series_back_refusals(a)", "_run_sorted_refusals(a)",
                 "_copy_rows_refusals(a)", "_board(", "load_stage("):
        assert at < src.index(call), call
    assert "mix_lean = _mix_lean_refusals(a) if r75 is None else False" in src
    assert "series_back = _series_back_refusals(a) if not (mix_lean or bayes or r75) else False" in src
    assert src.count("_r75_refusals(") == 1


# ------------------------------------------------------------------------------------------------ routing and refusals
@pytest.mark.parametrize("name", NAMES)
def test_each_file_is_taken_with_row_r75_and_nothing_but_the_ledger_is_read(bt, env, name):
    for argv in (ARGS[name], [*ARGS[name], "--board", "T2:embryo:test_interp"], [*ARGS[name], "--seed", "1000001"]):
        a = bt.build_parser().parse_args(["embryo", *argv, "--out", "x.h5ad"])
        assert bt._r75_file(a) == name and bt._r75_refusals(a) == _values(name)
    if EXPRESSION[name] == "mixture":  # the launcher's spelling with O's default named is the same command
        a = bt.build_parser().parse_args(["embryo", *FLOOR, "--expression", "mixture", *OV, *OPTS[name].split(),
                                          "--captain-ruling", "R75", "--out", "x.h5ad"])
        assert bt._r75_refusals(a) == _values(name)
    assert env == []


@pytest.mark.parametrize("name", NAMES)
@pytest.mark.parametrize("ruling", RULINGS_OTHER)
def test_each_files_string_under_any_other_ruling_or_none_is_refused(bt, env, tmp_path, name, ruling):
    argv = [*O_ARGS, *OPTS[name].split(), *([] if ruling is None else ["--captain-ruling", ruling])]
    a = bt.build_parser().parse_args(["embryo", *argv, "--out", "x.h5ad"])
    assert bt._r75_file(a) is None
    _refused(bt, env, tmp_path, argv, None)


@pytest.mark.parametrize("name", NAMES)
@pytest.mark.parametrize("change", [
    ["--seed", "1000002"], ["--seed", "1000000"], ["--move-coupling", "aligned"],
    ["--row-order", "copy_parity", "--row-order-test-ok"], ["--neighbourhood", "run_sorted"], ["--herd"],
    ["--draw", "midpoint"], ["--copy-rows", "--copy-rows-test-ok"], ["--board", DRY],
])
def test_every_other_use_of_a_files_string_under_r75_is_refused(bt, env, tmp_path, name, change):
    name_out = "T2embryo_dry.h5ad" if DRY in change else "T2embryo.h5ad"
    _refused(bt, env, tmp_path, [*ARGS[name], *change], "orders five files", name_out)


@pytest.mark.parametrize("argv", [
    [*O_ARGS, "--captain-ruling", "R75"],  # O itself
    [*OT_ARGS[:-2], "--captain-ruling", "R75"],  # row R59's file under row R75
    [*SO_ARGS, "--captain-ruling", "R75"],  # SO under row R75
    [*FLOOR, "--expression", "series_back", *OV, "--captain-ruling", "R75"],  # W under row R75
    [*FLOOR, "--expression", "series_back", *OV, "--mix-lean", "later_half", "--captain-ruling", "R75"],  # W, a lean
    [*FLOOR, "--expression", "render_hi", *OV, "--mix-lean", "later", "--captain-ruling", "R75"],
    [*FLOOR, "--expression", "depth", *OV, "--mix-lean", "later", "--captain-ruling", "R75"],  # depth with a lean
    [*FLOOR, "--expression", "depth", *OV, "--mix-lean", "later2", "--captain-ruling", "R75"],
    [*FLOOR, "--expression", "series", *OV, "--mix-lean", "later_half", "--captain-ruling", "R75"],  # another pair
    [*FLOOR, "--expression", "series", *OV, "--mix-lean", "later2", "--captain-ruling", "R75"],
    [*FLOOR, "--expression", "blend", *OV, "--mix-lean", "later_half", "--captain-ruling", "R75"],
    [*FLOOR, "--mix-lean", "later_half", "--captain-ruling", "R75"],  # no form
    [*FLOOR, "--shape-interp", "overlay", "--mix-lean", "later_half", "--captain-ruling", "R75"],  # no yes flag
    [*FLOOR, "--shape-interp", "stretch_ratio", "--shape-interp-test-ok", "--mix-lean", "later_half",
     "--captain-ruling", "R75"],
    [*FLOOR, "--shape-interp", "overlay_stretch", "--shape-interp-test-ok", "--expression", "depth",
     "--captain-ruling", "R75"],
    [*FLOOR, "--shape-interp", "overlay_bayes", "--shape-interp-test-ok", "--mix-lean", "later2",
     "--captain-ruling", "R75"],
    ["--method", "fgw", *OV[:2], "--mix-lean", "later_half", "--captain-ruling", "R75"],
    ["--method", "fgw_shape0", *OV, "--expression", "depth", "--captain-ruling", "R75"],  # positions moved
    ["--method", "fgw", "--expression", "depth", "--captain-ruling", "R75"],  # try E7's arm under row R75
])
def test_every_other_command_under_r75_is_refused(bt, env, tmp_path, argv):
    _refused(bt, env, tmp_path, argv, "orders five files")


@pytest.mark.parametrize("lean", ["later_half", "later2", "earlier_half"])
def test_a_lean_only_row_r75_has_needs_its_ruling(bt, env, tmp_path, lean):
    msg = _refused(bt, env, tmp_path, [*O_ARGS, "--mix-lean", lean], f"--mix-lean {lean}")
    assert "needs --captain-ruling R75" in msg
    for ruling in ("R59", "R48", "R71", "R60", "R76", "R750"):
        argv = [*O_ARGS, "--mix-lean", lean, "--captain-ruling", ruling]
        msg = _refused(bt, env, tmp_path, argv, f"--mix-lean {lean}")
        assert f"refused with --captain-ruling {ruling}: only row R75 orders it" in msg
    _refused(bt, env, tmp_path, ["--method", "fgw", "--mix-lean", lean], "needs --captain-ruling R75")


@pytest.mark.parametrize("extra", [[], ["--seed", "1000003"], ["--move-coupling", "aligned"], ["--board", DRY]])
def test_depth_with_the_overlay_form_needs_row_r75(bt, env, tmp_path, extra):
    name = "T2embryo_dry.h5ad" if DRY in extra else "T2embryo.h5ad"
    argv = [*FLOOR, "--expression", "depth", *OV, *extra]
    msg = _refused(bt, env, tmp_path, argv, "needs --captain-ruling R75", name)
    assert "row R75's file alone" in msg and bt.EMBRYO_R75_COMMAND.format(option="--expression depth") in msg


@pytest.mark.parametrize("name", NAMES)
@pytest.mark.parametrize("which", ["none", "twice", "unpaired", "unticked", "draft"])
def test_a_file_whose_row_does_not_name_it_is_refused(bt, env, tmp_path, monkeypatch, name, which):
    option = OPTS[name]
    other = [f for f in NAMES if f != name][0]
    rows = {"none": (R59_ROW,), "twice": (ROW_R75, ROW_R75),
            "unpaired": (ROW_R75.replace(f"`{option}` ({name})", f"`{option}` ({other})"),),
            "unticked": (ROW_R75.replace("`--captain-ruling R75`", "--captain-ruling R75"),),
            "draft": (ROW_R75.replace("(the captain's order", "(the AI assistant's order"),)}[which]
    p = tmp_path / "ledger" / "TRIES.md"
    p.parent.mkdir()
    p.write_text("| ID | what |\n|---|---|\n" + "".join(r + "\n" for r in rows))
    monkeypatch.setattr(bt, "TRIES_FILE", p)
    msg = _refused(bt, env, tmp_path, ARGS[name], f"row R75's file {name} is refused")
    assert f"`{option}` ({name}) and `--captain-ruling R75`" in msg


def test_the_help_names_row_r75s_files(bt):
    p = bt.build_parser()._subparsers._group_actions[0].choices["embryo"]
    helps = {a.dest: a.help for a in p._actions}
    for name, option in TABLE:
        assert f"`{option}` ({name})" in helps["mix_lean"]
    assert helps["mix_lean"].endswith(f"{bt.EMBRYO_R75_COMMAND.format(option='<OPTION>')}, t_mix 7/12, 5/6, 5/12, 2/3")
    assert helps["captain_ruling"].endswith("on O's command: R75; given once")
    assert "depth is row R75's file T2embryo_s2x_ovd.h5ad alone" in helps["expression"]
    assert p.format_help()  # the help still formats


# ------------------------------------------------------------------------------------------------ the builds
@pytest.mark.parametrize("name", LEANS)
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
    # every row's expression is its bracket's released cell (the mixture)
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


def test_the_series_file_is_sos_expression_on_row_r59s_rows(bt, env, tmp_path, monkeypatch, stages):
    seen = _capture(bt, monkeypatch)
    ot, _, _ = _build(bt, tmp_path, "ot", *OT_ARGS)
    s, rec, _ = _build(bt, tmp_path, "s", *ARGS[SERIES], name=SERIES)
    (a_t, kw_t, res_t), (a_s, kw_s, res_s) = seen
    assert kw_s == kw_t and kw_s["t_mix"] == bt.EMBRYO_OT23_T_MIX and a_s[4:] == a_t[4:]  # row R59's call exactly
    assert res_s.rows.tobytes() == res_t.rows.tobytes() and res_s.stage.tobytes() == res_t.stage.tobytes()
    Xt, Ct, nt_, vt = _file(ot)
    Xs, Cs, ns, vs = _file(s)
    assert Cs.tobytes() == Ct.tobytes() and ns == nt_ and vs == vt  # Ot23's coordinates bit for bit
    late = res_s.stage == 1
    np.testing.assert_array_equal(Xs[late], Xt[late])  # the E8.0 rows as released, as in Ot23
    assert not np.array_equal(Xs[~late], Xt[~late])  # the E7.5 rows rendered into E8.0's run profile and depth
    r = rec["record"]
    assert r["expression"] == "series" and r["series_rule"] == bt.t2_build.SERIES_RULE and r["t_mix"] == 2 / 3
    assert rec["r75"]["expression"] == "series" and rec["r75"]["t_mix_exact"] == "2/3"


def test_the_depth_file_keeps_os_rows_and_coordinates_bit_for_bit(bt, env, tmp_path, monkeypatch):
    seen = _capture(bt, monkeypatch)
    o, rec_o, _ = _build(bt, tmp_path, "o", *O_ARGS)
    d, rec_d, _ = _build(bt, tmp_path, "d", *ARGS[DEPTH], name=DEPTH)
    (a_o, kw_o, res_o), (a_d, kw_d, res_d) = seen
    assert kw_d == kw_o and "t_mix" not in kw_d and a_d[4:] == a_o[4:]  # interpolate_clouds called as O's
    Xo, Co, no, vo = _file(o)
    Xd, Cd, nd, vd = _file(d)
    assert Cd.tobytes() == Co.tobytes() and nd == no and vd == vo  # O's rows, order and coordinates bit for bit
    early = res_d.stage == 0
    np.testing.assert_array_equal(Xd[early], Xo[early])  # the shallower bracket is left as released
    assert not np.array_equal(Xd[~early], Xo[~early])  # the deeper one thinned to the predicted depth
    assert (Xd > 0).sum() < (Xo > 0).sum()
    r = rec_d["record"]
    assert r["expression"] == "depth" and r["depth_rule"] == bt.t2_build.DEPTH_RULE and r["t_mix"] == 0.5
    assert bt.t2_build.DEPTH_RULE in rec_d["r75"]["rule"] and rec_d["r75"]["t_mix"] is None
    assert {k: v for k, v in r.items() if not k.startswith("depth_") and k != "expression"} == rec_o["record"]


@pytest.mark.parametrize("name", NAMES)
def test_each_files_record(bt, env, tmp_path, capsys, name):
    _, rec, text = _build(bt, tmp_path, "f", *ARGS[name], name=name)
    blk = rec["r75"]
    t = T_MIX[name]
    what = ("the embryo's depth line reopened for this file alone, beside O and never in its place" if name == DEPTH
            else "the captain's named exception on the embryo's overlay form and occupancy term for this file"
            + (", and the lean-early mixing line declined on 2026-11-06 (e) reopened for this file alone"
               if name == "T2embryo_s2x_ovt512.h5ad" else ""))
    assert blk == {"option": OPTS[name], "file": name, "captain_ruling": "R75",
                   "order": bt.EMBRYO_R75_ORDER.format(option=OPTS[name], file=name, what=what),
                   "command": "embryo " + " ".join(ARGS[name]), "o_command": O_CMD,
                   "expression": EXPRESSION[name], "mix_lean": _values(name)["mix_lean"],
                   "t_mix": None if t is None else float(t),
                   "t_mix_exact": None if t is None else f"{t.numerator}/{t.denominator}",
                   "rule": bt.EMBRYO_R75_RULES[name]}
    assert all(not isinstance(v, (dict, list)) for v in blk.values())
    assert f"`{OPTS[name]}` ({name}) and `--captain-ruling R75`" in blk["order"]
    assert "under R22, sent at once in one ruling row, with no gate" in blk["order"]
    assert t is None or blk["rule"].startswith(f"t_mix = {t.numerator}/{t.denominator} in vec.t2_build.")
    assert rec["n_cells"] == "max_cells" and not {"mix_lean", "series_back", "bayes"} & set(rec)
    assert not _shows(text + capsys.readouterr().out, SECRET_MAX, SECRET_TRUTH)


# ------------------------------------------------------------------------------------------------ the base checker
BASES = [("T2embryo_s2x_ovt712.h5ad", "half", O_ARGS), ("T2embryo_s2x_ovt56.h5ad", "half", O_ARGS),
         ("T2embryo_s2x_ovt512.h5ad", "half", O_ARGS), (SERIES, "half", SO_ARGS), (SERIES, "mixture", OT_ARGS),
         (DEPTH, "mixture", O_ARGS)]


def _check(chk, tmp_path: Path, mode: str, name: str) -> tuple[int, Path]:
    out = tmp_path / f"{mode}" / f"{name.removesuffix('.h5ad')}_{mode}.h5ad"
    return chk.main([mode, "embryo", *ARGS[name], "--out", str(out)]), out


@pytest.mark.parametrize("name, mode, known", BASES)
def test_the_checker_gives_the_known_files_bytes_from_each_files_own_path(bt, chk, env, tmp_path, name, mode, known):
    k, _, _ = _build(bt, tmp_path, "known", *known)
    code, out = _check(chk, tmp_path, mode, name)
    assert code == 0 and out.read_bytes() == k.read_bytes()
    assert out.with_name(out.name + ".stamp.json").exists()  # stamped as any build: the launcher reads its fields
    f, _, _ = _build(bt, tmp_path, "file", *ARGS[name], name=name)
    assert f.read_bytes() != k.read_bytes()  # the file itself is none of the files it is checked against


def test_the_checkers_modes_are_the_launchers_and_each_files_change(bt, chk):
    sh = (REPO / "experiments" / "builds" / "embryo_r75_line.sh").read_text()
    table = dict(re.findall(r'^  (T2embryo_[a-z0-9_]+\.h5ad)\) OPT="[^"]+"; STEM=[a-z0-9]+; BASES="([^"]+)" ;;$', sh,
                            flags=re.M))
    keys = {tuple(O_ARGS): "O", tuple(SO_ARGS): "SO", tuple(OT_ARGS): "OT"}
    want = {}
    for name, mode, known in BASES:
        want.setdefault(name, []).append(f"{mode}:{keys[tuple(known)]}")
    assert {f: v.split() for f, v in table.items()} == want and chk.MODES == ("half", "mixture")
    assert chk.EXPRESSIONS == {"series": "series_rendered_expression", "depth": "depth_matched_expression"}
    for name in NAMES:  # half for each file with a lean, mixture for each with an expression
        modes = [m for f, m, _ in BASES if f == name]
        assert modes == ([] if T_MIX[name] is None else ["half"]) + ([] if EXPRESSION[name] == "mixture"
                                                                      else ["mixture"])


@pytest.mark.parametrize("mode, argv, match", [
    ("half", ARGS[DEPTH], "half takes a file with a lean"),
    *[("mixture", ARGS[f], "mixture takes a file with --expression") for f in LEANS],
    ("half", O_ARGS, "checks row R75's files only"),
    ("mixture", OT_ARGS, "checks row R75's files only"),
    ("half", [*ARGS["T2embryo_s2x_ovt712.h5ad"], "--seed", "1000002"], "orders five files"),
])
def test_the_checker_refuses_a_command_or_mode_it_does_not_check(bt, chk, env, tmp_path, mode, argv, match):
    out = tmp_path / "c" / "T2embryo_check.h5ad"
    with pytest.raises(SystemExit, match=match):
        chk.main([mode, "embryo", *argv, "--out", str(out)])
    assert not out.exists() and env == []


@pytest.mark.parametrize("argv", [[], ["half"], ["base", "embryo"], ["half", "t1"], ["Half", "embryo"],
                                  ["mixture", "heart"]])
def test_the_checkers_usage(chk, argv):
    with pytest.raises(SystemExit, match="usage: embryo_r75_check.py"):
        chk.main(argv)


@pytest.mark.parametrize("name, fn", [(SERIES, "series_rendered_expression"), (DEPTH, "depth_matched_expression")])
def test_the_checker_fails_when_the_expression_moves_the_interpolation(bt, chk, env, tmp_path, monkeypatch, name, fn):
    inner = getattr(bt.t2_build, fn)

    def moving(before, after, res, *a, **k):
        res.coords[0, 0] += 1e-9  # a step that touches the coordinates
        return inner(before, after, res, *a, **k)

    monkeypatch.setattr(bt.t2_build, fn, moving)
    code, _ = _check(chk, tmp_path, "mixture", name)
    assert code == 1


def test_the_checker_fails_when_the_route_interpolates_more_than_once(bt, chk, env, tmp_path, monkeypatch):
    """A build that makes one interpolation more than the route's (here a small synthetic one before build_embryo's
    own) is not the route the checker names: it exits 1 though the file's own interpolation is as named."""
    code, _ = _check(chk, tmp_path / "once", "half", "T2embryo_s2x_ovt712.h5ad")
    assert code == 0
    inner_build = bt.build_embryo

    def twice(a):
        rng = np.random.default_rng(0)
        C, E = rng.normal(size=(60, 3)), rng.gamma(2.0, 1.0, size=(60, 12))
        bt.t2_build.interpolate_clouds(C, C + 1.0, E, E, 0.5, 20, 1, target_rms=1.0)
        return inner_build(a)

    monkeypatch.setattr(bt, "build_embryo", twice)
    code, _ = _check(chk, tmp_path / "twice", "half", "T2embryo_s2x_ovt712.h5ad")
    assert code == 1


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
    """BASE's builder: row R75's route changes build_test.py only, and every module is this checkout's in both
    builders, so a difference below is the route's."""
    mod = _module("_build_test_at_base_r75", "experiments/builds/build_test.py")
    assert not hasattr(mod, "EMBRYO_R75_TABLE") and hasattr(mod, "EMBRYO_OT23_T_MIX")
    return mod


# name: (argv, file name, what the builder of BASE does with it: "built", "refused" or None for either)
COMMANDS = {
    "O": (O_ARGS, "T2embryo.h5ad", "built"),
    "O_launcher": ([*FLOOR, "--expression", "mixture", *OV], "T2embryo.h5ad", "built"),
    "Ot23": (OT_ARGS, "T2embryo.h5ad", "built"),
    "Ot23_launcher": ([*FLOOR, "--expression", "mixture", *OV, "--mix-lean", "later", "--captain-ruling", "R59"],
                      "T2embryo.h5ad", "built"),
    "W": (W_ARGS, "T2embryo.h5ad", "built"),
    "SO": (SO_ARGS, "T2embryo.h5ad", "built"),
    "O_hi": ([*FLOOR, "--expression", "render_hi", *OV], "T2embryo.h5ad", "built"),
    "P": (FLOOR, "T2embryo.h5ad", "built"),
    "P_depth": ([*FLOOR, "--expression", "depth"], "T2embryo.h5ad", "built"),  # depth without the overlay form
    "s1": (["--method", "fgw"], "T2embryo.h5ad", "built"),
    "s1_depth": (["--method", "fgw", "--expression", "depth"], "T2embryo.h5ad", "built"),  # try E7's arm
    "s2": (["--method", "fgw_shape0"], "T2embryo.h5ad", "built"),
    "series": (["--method", "fgw", "--expression", "series"], "T2embryo.h5ad", "built"),
    "midpoint": (["--method", "fgw", "--draw", "midpoint"], "T2embryo.h5ad", "built"),
    "seed": ([*O_ARGS, "--seed", "1000003"], "T2embryo.h5ad", "built"),
    "dry_O": ([*FLOOR, "--shape-interp", "overlay", "--board", DRY], "T2embryo_dry.h5ad", "built"),
    "dry_SO": ([*FLOOR, "--expression", "series", "--shape-interp", "overlay", "--board", DRY], "T2embryo_dry.h5ad",
               "built"),
    "depth_stretch": ([*FLOOR, "--expression", "depth", "--shape-interp", "stretch_ratio", "--shape-interp-test-ok"],
                      "T2embryo.h5ad", None),
    "RS": ([*O_ARGS, "--neighbourhood", "run_sorted", "--captain-ruling", "R60"], "T2embryo.h5ad", None),
    "BF": ([*FLOOR, "--shape-interp", "overlay_bayes", "--shape-interp-test-ok", "--captain-ruling", "R71"],
           "T2embryo_s2x_ovbf.h5ad", None),
    "copy_parity": ([*FLOOR, "--row-order", "copy_parity", "--row-order-test-ok"], "T2embryo.h5ad", None),
    "copy_rows": ([*FLOOR, "--copy-rows", "--copy-rows-test-ok"], "T2embryo.h5ad", None),
    # the refusals as they were, rows R48's and R59's among them
    "W_no_ruling": ([*FLOOR, "--expression", "series_back", *OV], "T2embryo.h5ad", "refused"),
    "W_lean": ([*W_ARGS, "--mix-lean", "later"], "T2embryo.h5ad", "refused"),
    "Ot23_no_ruling": ([*O_ARGS, "--mix-lean", "later"], "T2embryo.h5ad", "refused"),
    "Ot23_R48": ([*O_ARGS, "--mix-lean", "later", "--captain-ruling", "R48"], "T2embryo.h5ad", "refused"),
    "Ot23_seed": ([*OT_ARGS, "--seed", "1000002"], "T2embryo.h5ad", "refused"),
    "series_later_none": ([*SO_ARGS, "--mix-lean", "later"], "T2embryo.h5ad", "refused"),
    "series_later_R59": ([*SO_ARGS, "--mix-lean", "later", "--captain-ruling", "R59"], "T2embryo.h5ad", "refused"),
    "depth_R59": ([*FLOOR, "--expression", "depth", *OV, "--captain-ruling", "R59"], "T2embryo.h5ad", "refused"),
    "depth_R48": ([*FLOOR, "--expression", "depth", *OV, "--captain-ruling", "R48"], "T2embryo.h5ad", "refused"),
    "depth_R71": ([*FLOOR, "--expression", "depth", *OV, "--captain-ruling", "R71"], "T2embryo.h5ad", "refused"),
    "earlier": ([*O_ARGS, "--mix-lean", "earlier", "--captain-ruling", "R59"], "T2embryo.h5ad", "refused"),
    "ruling_R48_alone": ([*O_ARGS, "--captain-ruling", "R48"], "T2embryo.h5ad", "refused"),
    "ruling_R59_alone": ([*O_ARGS, "--captain-ruling", "R59"], "T2embryo.h5ad", "refused"),
    "ruling_R60_alone": ([*O_ARGS, "--captain-ruling", "R60"], "T2embryo.h5ad", "refused"),
    "ruling_R76_alone": ([*O_ARGS, "--captain-ruling", "R76"], "T2embryo.h5ad", "refused"),
    "form_needs_ok": ([*FLOOR, "--shape-interp", "overlay"], "T2embryo.h5ad", "refused"),
    "positions": (["--method", "fgw", "--positions", "floor_exact"], "T2embryo.h5ad", "refused"),
    "seed_low": ([*O_ARGS, "--seed", "5"], "T2embryo.h5ad", "refused"),
    "dry_name": ([*FLOOR, "--board", DRY], "T2embryo.h5ad", "refused"),
    "test_name": (FLOOR, "T2embryo_dry.h5ad", "refused"),
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
def test_row_r75s_five_commands_were_refused_at_base_and_are_built_now(bt, base_bt, own_repo, tmp_path, stages, name):
    want = _outcome(base_bt, tmp_path, stages, "base", ARGS[name], name)
    got = _outcome(bt, tmp_path, stages, "new", ARGS[name], name)
    assert want[0] == "refused" and got[0] != "refused"
    assert json.loads(got[1])["r75"]["file"] == name


@pytest.mark.parametrize("extra", [[], ["--seed", "1000003"]])
def test_depth_with_the_overlay_form_and_no_ruling_was_built_at_base_and_is_refused_now(bt, base_bt, own_repo, tmp_path,
                                                                                       stages, extra):
    argv = [*FLOOR, "--expression", "depth", *OV, *extra]
    want = _outcome(base_bt, tmp_path, stages, "base", argv, "T2embryo.h5ad")
    got = _outcome(bt, tmp_path, stages, "new", argv, "T2embryo.h5ad")
    assert want[0] != "refused" and got[0] == "refused" and "needs --captain-ruling R75" in got[1]


@pytest.mark.parametrize("lean", ["later_half", "later2", "earlier_half"])
def test_a_new_lean_without_its_ruling_is_refused_by_both(bt, base_bt, own_repo, tmp_path, stages, lean):
    argv = [*O_ARGS, "--mix-lean", lean]
    want = _outcome(base_bt, tmp_path, stages, "base", argv, "T2embryo.h5ad")
    got = _outcome(bt, tmp_path, stages, "new", argv, "T2embryo.h5ad")
    assert want == ("refused", "2")  # BASE's parser knew no such lean
    assert got[0] == "refused" and "needs --captain-ruling R75" in got[1]
