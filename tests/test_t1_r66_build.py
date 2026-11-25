"""`build_test.py t1` under row R66 of experiments/harness/TRIES.md (the captain's order of 2026-11-25, T1's third-round
step-on slots): thirty files on T1_s2rcMS20's command (`t1 --g-mom --base lib3 --feng-groups recomputed
--feng-fraction-mult 2.0 --direction-momentum all`, row R51's file), each taken only with --captain-ruling R66 through
_captain_order_names on row R66's first cell, each as its file's name gives it (T1_s2rcMS20 plus the suffixes in d, H,
m order):
  K1, K2, K3  two or three of line 12's steps together: --g-mom-damping 1.0 | 1.5 | 2.0 (d10, d15, d20),
              --held-momentum 0.5 | 1.0 (H05, H10), --direction-momentum-damping 1.0 | 2.0 (m10, m20)
  C1          --g-mom-damping 2.5 (T1_s2rcMS20d25)
  C3          --held-momentum 2.0 (T1_s2rcMS20H20)
Blood moves as in row R60's H files (R1 reversed for Blood, T64's close lifted) in every file with --held-momentum.

Held here on synthetic stages and directions only (the fixtures of tests/test_t1_r51_build.py), with temporary ledgers
holding a copy of row R66's real first cell (one test checks the copy against the repository's ledger, by hash,
printing neither): (1) the builder's one table is the row's thirty (file, option string) pairs in the row's order, each
name the one its string gives, and the values it routes on are read off the strings; (2) each file is taken with row
R66 only, every other use of it under R66 is refused before anything but the ledger is read, and each set under any
other ruling or none, and every other set under R66, is routed as before and refused there; (3) the builder at 07ca2f1
(read with `git show`) and this one route every argument set of a grid alike but the thirty sets under R66; (4) each
file's steps get the file's values and its record names the row, the slot, the base file and each change; Blood's
rows move only in files with Blood's step, every other row of such a file being the same command's without that step
byte for byte, and Blood's rows of every other file are T1_s2rcMS20's; with every step at T1_s2rcMS20's value
(experiments/builds/t1_r66_check.py base) each file is T1_s2rcMS20's bytes; a file builds the same bytes twice; MS20's
and line 12's files are the same bytes and records with the new route taken out. No released file, outside-atlas
file, board or score is read."""

from __future__ import annotations

import hashlib
import importlib.util
import itertools
import json
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
for _p in ("experiments/external", "experiments/t1_m1", "experiments/t1_direction", "experiments/t1_hygiene",
           "experiments/builds"):
    sys.path.insert(0, str(ROOT / _p))

import build_test as btm  # noqa: E402  (the module the fixture `bt` returns)
import feng_relabel as fr  # noqa: E402
import test_t1_l12_build as l12t  # noqa: E402
import test_t1_r51_build as r51t  # noqa: E402
from conftest import git_env  # noqa: E402
from test_hidden_counts import t1_env  # noqa: E402, F401  (the synthetic build_t1 fixture)
from test_t1_candidates_build import RC, _build, _no_read, _same_frame  # noqa: E402
from test_t1_fraction2_build import _record_sha  # noqa: E402

bt, ms4 = r51t.bt, r51t.ms4  # the builder module and the three-library recomputed setup

SEED = r51t.SEED
BASE_CMD = "t1 --g-mom --base lib3 --feng-groups recomputed"
MS20_OPT = "--feng-fraction-mult 2.0 --direction-momentum all"
MS20 = [*MS20_OPT.split(), "--captain-ruling", "R51"]  # row R51's file, the base of row R66's
STEP_D = {"d10": "1.0", "d15": "1.5", "d20": "2.0"}  # line 12's three values of each axis, in the row's words
STEP_H = {"H05": "0.5", "H10": "1.0"}
STEP_M = {"m10": "1.0", "m20": "2.0"}
FLAG = {"d": "--g-mom-damping", "H": "--held-momentum", "m": "--direction-momentum-damping"}


def _k_table() -> list[tuple[str, str]]:
    """K1, K2 and K3's twenty-eight (file, option string) pairs in the row's order: d with H, d with m, H with m, then
    all three, each axis's values in line 12's order."""
    out = []
    for axes in (("d", "H"), ("d", "m"), ("H", "m"), ("d", "H", "m")):
        values = {"d": STEP_D, "H": STEP_H, "m": STEP_M}
        for combo in itertools.product(*(values[x].items() for x in axes)):
            name = "T1_s2rcMS20" + "".join(suffix for suffix, _ in combo)
            out.append((name, " ".join([MS20_OPT, *(f"{FLAG[s[0]]} {v}" for s, v in combo)])))
    return out


# the thirty (file, option string) pairs, written here from row R66's words, in its order
TABLE = [*_k_table(), ("T1_s2rcMS20d25", f"{MS20_OPT} --g-mom-damping 2.5"),
         ("T1_s2rcMS20H20", f"{MS20_OPT} --held-momentum 2.0")]
OPTS = dict(TABLE)
ARGS = {f: [*o.split(), "--captain-ruling", "R66"] for f, o in TABLE}
# A copy of row R66's first cell as the repository's ledger holds it (test_the_copy_is_the_ledgers_own)
REAL_R66 = (
    "R66 (the captain's order of 2026-11-25, ruling, no run: \"Yes to K1, K2, K3, C1 and C3 of the round's T1 step-on "
    "design under R22, as gated (headlines only, read once by 26 Nov 06:00 UTC), in one ruling row; R1 reversed for "
    "Blood and T64's close lifted in the files with a Blood step, as in R60.\"; T1's third-round step-on slots on the "
    "banked file T1_s2rcMS20; the T1 builder is to take on `t1 --g-mom --base lib3 --feng-groups recomputed` exactly "
    "these thirty option strings, each only with `--captain-ruling R66` and each as its file's name gives it "
    "(`T1_s2rcMS20` plus the suffixes in d, H, m order): for K1, K2 and K3, two or three of the steps together, "
    + ", ".join(f"`{o}` ({f})" for f, o in TABLE[:28])
    + f"; for C1, `{OPTS['T1_s2rcMS20d25']}` (T1_s2rcMS20d25); for C3, `{OPTS['T1_s2rcMS20H20']}` (T1_s2rcMS20H20); "
    "decision log 2026-11-25 (hh); not a try)")
ROW_R66 = f"| {REAL_R66} | T1 | step-on |"
# Numbers a refusal may show: line 12's (tests/test_t1_l12_build.py) and row R66's ID and C1's and C3's values
NUMBERS_ALLOWED = l12t.NUMBERS_ALLOWED | {"66", "25", "2.5"}
RULINGS_OTHER = [None, "R51", "R60", "R46", "R22", "R23", "R44", "R45", "R58", "R660", "r66", "R66a", "R6"]
ORDER = ("experiments/harness/TRIES.md row R66, the captain's order of 2026-11-25 (\"ruling, no run\"; decision log "
         "2026-11-25 (hh), T1's third-round step-on {slot}) naming `{option}` and `--captain-ruling R66`: file {file}, "
         "T1_s2rcMS20's command (row R51's file: candidate MS's step exactly as row R23 fixes it, the atlas step at "
         "2.0 x the calendar fraction, g_mom's within-type damping 0.5, MS's damping harness.B0_DAMPING) with "
         "{changes}; nothing else changed{blood}")
BLOOD_RULING = ("; R1 reversed for Blood and T64's close lifted in this file, as row R66 orders for each file with a "
                "Blood step and as row R60 did in its H files")


def _steps(name: str) -> list[str]:
    """The name's suffixes in order (d15H05m20 -> d15, H05, m20)."""
    suffix = name.removeprefix("T1_s2rcMS20")
    steps = re.findall(r"[dHm]\d\d", suffix)
    assert "".join(steps) == suffix, name
    return steps


def _values(name: str) -> tuple:
    """(multiplier, g_mom's damping, Blood's damping or None, MS's damping) as the file's name gives them."""
    got = {s[0]: int(s[1:]) / 10 for s in _steps(name)}
    return 2.0, got.get("d", 0.5), got.get("H"), got.get("m", 0.5)


def _changes(name: str) -> list[str]:
    """The file's changes in row R60's words for each step (tests/test_t1_r51_build.py's _want writes them out)."""
    mult, damping, held, ms = _values(name)
    out = []
    for step in _steps(name):
        if step[0] == "d":
            out.append(f"g_mom's within-type damping {damping:.1f}")
        elif step[0] == "H":
            out.append(f"Blood's own released E9.5 -> E10.5 change at damping {held:.1f} on the drawn Blood cells, "
                       "after MS's step and before g_mom (R1 reversed for Blood in this one file, Macrophage and Early "
                       "VIC held)")
        else:
            out.append(f"MS's own damping {ms:.1f} in place of harness.B0_DAMPING")
    return out


def _slot(name: str) -> str:
    return {"T1_s2rcMS20d25": "slot C1", "T1_s2rcMS20H20": "slot C3"}.get(name, "slots K1, K2 and K3")


def _want(name: str) -> dict:
    """_t1_l12's values for a file, written out here from the row's words."""
    mult, damping, held, ms = _values(name)
    changes = _changes(name)
    count = {1: "one change", 2: "two changes", 3: "three changes"}[len(changes)]
    listed = changes[0] if len(changes) == 1 else ", ".join(changes[:-1]) + " and " + changes[-1]
    order = ORDER.format(slot=_slot(name), option=OPTS[name], file=name, changes=f"{count}, {listed}",
                         blood=BLOOD_RULING if held is not None else "")
    return {"direction_momentum": "all", "captain_ruling": "R66", "file": name, "base_file": "T1_s2rcMS20",
            "option": f"`{OPTS[name]}`", "feng_fraction_mult": mult, "g_mom_damping": damping, "held_momentum": held,
            "direction_momentum_damping": ms, "order": order,
            "command": f"{BASE_CMD} {OPTS[name]} --captain-ruling R66", "slot": _slot(name), "changes": changes}


def _ledger_text(r66: str | None = ROW_R66) -> str:
    return l12t._ledger_text() + (f"{r66}\n" if r66 else "")


@pytest.fixture
def ledger(bt, tmp_path, monkeypatch):
    """A temporary ledger: tests/test_t1_l12_build.py's (R20, R22, R23, R26, R37, synthetic R44, R45 and R51, and
    copies of rows R60's and R46's real first cells) with a copy of row R66's real first cell; the builder reads it in
    place of the repository's."""
    p = tmp_path / "TRIES.md"
    p.write_text(_ledger_text())
    monkeypatch.setattr(bt, "TRIES_FILE", p)
    return p


def _parse(argv: list[str]):
    return btm.build_parser().parse_args(["t1", *argv, "--out", "x.h5ad"])


def _refused(monkeypatch, tmp_path, argv: list[str], match: str | None = None) -> str:
    """Refused before anything but the ledger is read, showing no number but the allowed ones and the command's."""
    for obj, name in ((btm, "load_stage"), (btm, "load_boards"), (btm.hygiene, "load_union"),
                      (btm.feng_dirs, "feng_directions"), (fr, "committed_at_head"), (btm, "_t1_hm_step"),
                      (btm, "_t1_ms_step")):
        monkeypatch.setattr(obj, name, _no_read)
    name = "T1_dry.h5ad" if "--board" in argv else "T1.h5ad"
    with pytest.raises(SystemExit, match=match) as e:
        btm.main([*argv, "--out", str(tmp_path / name)])
    msg = str(e.value)
    given = set(re.findall(r"\d+(?:\.\d+)?", " ".join(argv)))
    assert set(re.findall(r"\d+(?:\.\d+)?", msg)) <= NUMBERS_ALLOWED | given, msg
    return msg


def _check():
    import t1_r66_check

    return t1_r66_check


# ------------------------------------------------------------------ the one table and the ruling check


def test_the_builders_one_table_is_the_rows_thirty_pairs_in_its_order():
    assert len(TABLE) == 30 == len(OPTS) == len(set(OPTS.values()))
    assert [tuple(p) for p in btm.T1_R66_TABLE] == TABLE
    assert list(btm.T1_R66_FILES) == list(OPTS) and {f: v[0] for f, v in btm.T1_R66_FILES.items()} == OPTS
    # the row pairs each string with its file, in this order, and names nothing else in backticks but the command,
    # the ruling and the stem
    assert [(f, o) for o, f in re.findall(r"`([^`]*)` \((T1_s2rcMS20[A-Za-z0-9]*)\)", REAL_R66)] == TABLE
    assert re.findall(r"`([^`]*)`", REAL_R66) == [BASE_CMD, "--captain-ruling R66", "T1_s2rcMS20", *OPTS.values()]
    assert btm.T1_R66_ROW == "R66" and btm.T1_R66_BASE_FILE == "T1_s2rcMS20"
    for _, o in TABLE:
        assert btm.T1_R66_COMMAND.format(option=o) == f"{BASE_CMD} {o} --captain-ruling R66"


@pytest.mark.parametrize("name", list(OPTS))
def test_each_name_is_the_one_its_option_string_gives_and_the_values_are_read_off_the_string(name):
    option, *values = btm.T1_R66_FILES[name]
    a = _parse([*RC, *option.split(), "--captain-ruling", "R66"])
    assert tuple(values) == (a.feng_fraction_mult, a.g_mom_damping, a.held_momentum, a.direction_momentum_damping)
    mult, damping, held, ms = _values(name)
    assert values[0] == mult and (values[1] if values[1] is not None else 0.5) == damping
    assert (None if values[2] is None else float(values[2])) == held
    assert (0.5 if values[3] is None else float(values[3])) == ms
    # the suffixes in d, H, m order, each the option's flag with its value, and nothing else in the string
    want = " ".join([MS20_OPT, *(f"{FLAG[s[0]]} {int(s[1:]) / 10:.1f}" for s in _steps(name))])
    assert option == want


def test_the_copy_is_the_ledgers_own():
    """The copy above is row R66's first cell in the repository's ledger (compared by hash only)."""
    got = []
    for line in (ROOT / "experiments" / "harness" / "TRIES.md").read_text(encoding="utf-8").splitlines():
        if line.startswith("| R66 "):
            got.append(hashlib.sha256(line[2:].split("|", 1)[0].strip().encode()).hexdigest())
    same = got == [hashlib.sha256(REAL_R66.encode()).hexdigest()]
    assert same, "the copy is not the ledger's first cell (hashes only)"


def test_the_three_rows_sets_are_disjoint_and_no_earlier_file_is_named_again():
    r66 = {tuple(v[1:]) for v in btm.T1_R66_FILES.values()}
    r60 = {tuple(v[1:]) for v in btm.T1_R60_FILES.values()}
    r51 = {tuple(v[1:]) for v in btm.T1_R51_FILES.values()}
    assert len(r66) == 30 and not r66 & r60 and not r66 & r51 and btm.T1_R46_VALUES not in r66
    earlier = {*btm.T1_R60_FILES, *btm.T1_R51_FILES, btm.T1_R46_FILE, btm.T1_R58_FILE, btm.T1_R54_FILE,
               btm.T1_MS25_FILE, *btm.T1_HM_FILES.values(), *btm.T1_LINE6_FILES.values(), btm.T1_LINE7_FILE}
    assert not set(btm.T1_R66_FILES) & earlier
    assert {f"T1_s2rcMS20{s}" for s in (*STEP_D, *STEP_H, *STEP_M)} == set(btm.T1_R60_FILES)  # line 12's steps


@pytest.mark.parametrize("bad", [
    MS20_OPT,  # no change
    f"{MS20_OPT} --held-momentum 0.5 --g-mom-damping 1.0",  # not in d, H, m order
    f"{MS20_OPT} --g-mom-damping 1.0 --g-mom-damping 1.5",  # a flag twice
    f"{MS20_OPT} --g-mom-damping",  # no value
    f"{MS20_OPT} --herd x",
    "--feng-fraction-mult 2.5 --direction-momentum all --g-mom-damping 1.0",
    "--direction-momentum all --feng-fraction-mult 2.0 --g-mom-damping 1.0",
])
def test_a_string_off_the_rows_form_is_not_read(bad):
    with pytest.raises(ValueError):
        btm._r66_values(bad)


def test_each_files_changes_and_slot_are_line_12s_words_for_its_steps():
    for name in OPTS:
        assert list(btm.T1_R66_CHANGES[name]) == _changes(name) and btm.T1_R66_SLOTS[name] == _slot(name), name
        for step, text in zip(_steps(name), btm.T1_R66_CHANGES[name], strict=True):
            line12 = f"T1_s2rcMS20{step}"
            if line12 in btm.T1_R60_CHANGES:  # K's steps: word for word row R60's change for that step
                assert text == btm.T1_R60_CHANGES[line12], name
    assert [f for f in OPTS if btm.T1_R66_SLOTS[f] != "slots K1, K2 and K3"] == ["T1_s2rcMS20d25", "T1_s2rcMS20H20"]


def test_row_r66s_first_cell_names_each_set_through_captain_order_names(bt, ledger, monkeypatch):
    seen = []
    inner = bt._captain_order_names
    monkeypatch.setattr(bt, "_captain_order_names", lambda row, opt, *a, **k: seen.append((row, opt)) or inner(
        row, opt, *a, **k))
    assert all(bt._l12_row_names("R66", f) for f in OPTS)
    assert seen == [("R66", f"`{o}`") for o in OPTS.values()]


def test_row_r66s_check_fails_on_a_cell_that_lacks_a_set(bt, tmp_path):
    p = tmp_path / "TRIES.md"
    for o in OPTS.values():
        p.write_text(f"# ledger\n| {REAL_R66.replace(f'`{o}`', f'`{o} --herd`')} | T1 |\n")
        assert [bt._captain_order_names("R66", f"`{x}`", p) for x in OPTS.values()] == [x != o for x in OPTS.values()]
    p.write_text(f"# ledger\n| {REAL_R66.replace('`--captain-ruling R66`', '`--captain-ruling R60`')} | T1 |\n")
    assert not any(bt._captain_order_names("R66", f"`{x}`", p) for x in OPTS.values())
    p.write_text(f"# ledger\n{ROW_R66}\n{ROW_R66}\n")  # two rows with one ID
    assert not any(bt._captain_order_names("R66", f"`{x}`", p) for x in OPTS.values())


# ------------------------------------------------------------------ routing and refusals


@pytest.mark.parametrize("name", list(OPTS))
def test_each_file_is_taken_with_row_r66_and_nothing_but_the_ledger_is_read(bt, ledger, monkeypatch, name):
    for obj, attr in ((bt, "load_stage"), (bt, "load_boards"), (fr, "committed_at_head")):
        monkeypatch.setattr(obj, attr, _no_read)
    a = _parse([*RC, *ARGS[name]])
    assert bt._l12_hit(a) == ("R66", name)
    assert bt._t1_l12(a) == _want(name) and bt._t1_r51(a) == _want(name)
    assert bt._l12_ordered(a) and bt._bet_refusals(a) is None and bt._t1_pooled(a) is None
    assert bt._t1_r54(a) is None and bt._t1_r58(a) is None


@pytest.mark.parametrize("name", list(OPTS))
@pytest.mark.parametrize("ruling", RULINGS_OTHER)
def test_each_set_under_any_other_ruling_is_refused(bt, ledger, monkeypatch, tmp_path, name, ruling):
    argv = ["t1", *RC, *OPTS[name].split(), *([] if ruling is None else ["--captain-ruling", ruling])]
    a = _parse(argv[1:])
    assert bt._t1_l12(a) is None and not bt._l12_ordered(a)  # routed as before, where it is refused
    _refused(monkeypatch, tmp_path, argv)


@pytest.mark.parametrize("extra", [
    MS20_OPT.split(),  # MS20's own set: row R51's file
    *[o.split() for o in l12t.OPTS.values()],  # row R60's seven sets
    l12t.PAIR_OPT.split(),  # row R46's pair
    [*MS20_OPT.split(), "--g-mom-damping", "2.5", "--held-momentum", "0.5"],  # C1's value with another step
    [*MS20_OPT.split(), "--g-mom-damping", "2.5", "--direction-momentum-damping", "1.0"],
    [*MS20_OPT.split(), "--g-mom-damping", "1.0", "--held-momentum", "2.0"],  # C3's value with another step
    [*MS20_OPT.split(), "--held-momentum", "2.0", "--direction-momentum-damping", "2.0"],
    [*MS20_OPT.split(), "--g-mom-damping", "0.5", "--held-momentum", "0.5"],  # the base's value given as a step
    [*MS20_OPT.split(), "--g-mom-damping", "3.0", "--held-momentum", "0.5"],
    [*MS20_OPT.split(), "--g-mom-damping", "1.0", "--held-momentum", "0.50"],  # the literal strings only
    [*MS20_OPT.split(), "--g-mom-damping", "1.0", "--direction-momentum-damping", "1"],
    [*MS20_OPT.split(), "--held-momentum", "0.5", "--direction-momentum-damping", "3.0"],
    [*MS20_OPT.split(), "--held-momentum", "1.5", "--direction-momentum-damping", "1.0"],
    ["--feng-fraction-mult", "1.5", "--direction-momentum", "all", "--g-mom-damping", "1.0", "--held-momentum", "0.5"],
    ["--feng-fraction-mult", "2.5", "--direction-momentum", "all", "--g-mom-damping", "1.0", "--held-momentum", "0.5"],
    ["--feng-fraction-mult", "3.0", "--direction-momentum", "all", "--held-momentum", "2.0"],
    ["--feng-fraction-mult", "2.0", "--g-mom-damping", "1.0", "--held-momentum", "0.5"],  # without MS
    ["--direction-momentum", "all", "--g-mom-damping", "1.0", "--held-momentum", "0.5"],  # without the multiplier
    ["--g-mom-damping", "2.5"],
    [],
])
def test_every_other_set_under_r66_is_refused(bt, ledger, monkeypatch, tmp_path, extra):
    argv = ["t1", *RC, *extra, "--captain-ruling", "R66"]
    a = _parse(argv[1:])
    assert bt._t1_l12(a) is None and not bt._l12_ordered(a)  # routed as before
    _refused(monkeypatch, tmp_path, argv)


@pytest.mark.parametrize("name", list(OPTS))
@pytest.mark.parametrize("cmd, extra, match", [
    (RC, ["--released-bet", "share_trend", "--bet-strength", "0.5"], "refused with --released-bet, --bet-strength"),
    (RC, ["--bet-strength", "0.5"], "refused with --bet-strength"),
    (RC, ["--no-hold-unbracketed"], "R1's hold"),
    (RC, ["--seed", "1000002"], "default seed on the test board"),
    (RC, ["--board", "T1:val"], "default seed on the test board"),
    (["--g-mom", "--base", "lib3", "--feng-groups", "none"], [], "default seed on the test board"),
])
def test_every_other_use_of_a_file_under_r66_is_refused(bt, ledger, monkeypatch, tmp_path, name, cmd, extra, match):
    msg = _refused(monkeypatch, tmp_path, ["t1", *cmd, *ARGS[name], *extra], match)
    assert f"row R66's file {name}" in msg


@pytest.mark.parametrize("name", list(OPTS))
def test_a_file_whose_row_does_not_name_it_is_refused(bt, ledger, monkeypatch, tmp_path, name):
    ledger.write_text(_ledger_text(r66=ROW_R66.replace(f"`{OPTS[name]}`", "")))
    _refused(monkeypatch, tmp_path, ["t1", *RC, *ARGS[name]], "has no single row R66")
    ledger.write_text(_ledger_text(r66=ROW_R66.replace("R66 (the captain's order", "R66 (a proposal")))
    _refused(monkeypatch, tmp_path, ["t1", *RC, *ARGS[name]], "has no single row R66")
    ledger.write_text(_ledger_text(r66=None))
    _refused(monkeypatch, tmp_path, ["t1", *RC, *ARGS[name]], "has no single row R66")


def test_the_help_names_row_r66s_files():
    import io
    from contextlib import redirect_stdout

    buf = io.StringIO()
    with redirect_stdout(buf), pytest.raises(SystemExit):
        btm.build_parser().parse_args(["t1", "--help"])
    text = re.sub(r"-\s+", "-", " ".join(buf.getvalue().split()))
    for word in ("With --captain-ruling R66 only: row R66's thirty files", *OPTS, *OPTS.values(),
                 "--captain-ruling R66", "T1_s2rcMS20H20"):
        assert word in text, word


# ------------------------------------------------------------------ the builder at 07ca2f1 and this one, side by side


OLD_COMMIT = "07ca2f1"  # main before row R66's routes
OLD_SHA = "86bf0813ae380a85df915ee0f016649c0df31f732cc3a32fa1c15c10c47d06eb"  # its experiments/builds/build_test.py
CMDS = [RC, ["--g-mom", "--base", "lib3"], ["--g-mom", "--base", "lib3", "--feng-groups", "none"]]
DAMPINGS = [None, "0.5", "1.0", "1.5", "2.0", "2.5"]
MULTS = [None, "1.5", "2.0", "2.5"]
HELD = [None, "0.5", "1.0", "2.0"]
DMS = [None, "all"]
DMD = [None, "1.0", "2.0"]
GRID_RULINGS = [None, "R51", "R60", "R46", "R66", "R22"]
GRID_SETS = 3 * 6 * 4 * 4 * 2 * 3 * 6


def _old_builder(tmp_path: Path):
    env = {**git_env(), "GIT_OPTIONAL_LOCKS": "0"}
    res = subprocess.run(["git", "-C", str(ROOT), "show", f"{OLD_COMMIT}:experiments/builds/build_test.py"],
                         capture_output=True, env=env)
    if res.returncode != 0:
        pytest.skip(f"commit {OLD_COMMIT} is not available in this checkout")
    assert hashlib.sha256(res.stdout).hexdigest() == OLD_SHA
    path = tmp_path / "old" / "build_test_before_r66.py"
    path.parent.mkdir()
    path.write_bytes(res.stdout)
    spec = importlib.util.spec_from_file_location("build_test_before_r66", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _grid() -> list[list[str]]:
    sets = []
    for cmd, damping, mult, held, dm, dmd, ruling in itertools.product(CMDS, DAMPINGS, MULTS, HELD, DMS, DMD,
                                                                        GRID_RULINGS):
        out = ["t1", *cmd]
        for flag, value in (("--g-mom-damping", damping), ("--feng-fraction-mult", mult), ("--held-momentum", held),
                            ("--direction-momentum", dm), ("--direction-momentum-damping", dmd),
                            ("--captain-ruling", ruling)):
            out += [] if value is None else [flag, value]
        sets.append(out)
    return sets


def test_the_builder_at_07ca2f1_and_this_one_route_every_set_alike_but_r66s(bt, tmp_path, monkeypatch):
    """Side by side on the grid, under the ledger with row R66 and without it: every set other than one of the thirty
    option sets under R66 is refused with the same message, or reaches the same first read with the same candidate,
    MS and R51-route values, in both builders; the old builder refuses each of the thirty, and this one takes each on
    the recomputed command while the ledger names it, refusing it otherwise. Nothing is built: every read is stopped
    where it would start."""
    old = _old_builder(tmp_path)
    assert not hasattr(old, "T1_R66_FILES") and hasattr(bt, "T1_R66_FILES")
    ledger = tmp_path / "TRIES.md"
    record = tmp_path / "relabel" / "MANIFEST.json"
    monkeypatch.setattr(fr, "RECORD", record)
    monkeypatch.setattr(fr, "committed_at_head",
                        lambda p: json.dumps({"route": "recomputed"}).encode() if Path(p) == record else None)
    monkeypatch.setattr(fr, "later_ruling", lambda *a, **k: "no later ruling (synthetic)")
    monkeypatch.setattr(bt.hygiene, "load_union", r51t._reach("union"))
    monkeypatch.setattr(bt.strong_null, "check_approved", r51t._reach("strong_null"))
    monkeypatch.setattr(bt.feng_dirs, "feng_directions", r51t._reach("feng"))
    spies, parsers = {}, {}
    for mod in (old, bt):
        parsers[mod] = mod.build_parser()
        monkeypatch.setattr(mod, "build_parser", lambda p=parsers[mod]: p)
        monkeypatch.setattr(mod, "TRIES_FILE", ledger)
        monkeypatch.setattr(mod, "load_stage", r51t._reach("stage"))
        monkeypatch.setattr(mod, "load_boards", r51t._reach("board"))
        for name in ("_t1_ms", "_t1_candidate", "_t1_pooled", "_t1_r51"):
            inner = getattr(mod, name)
            monkeypatch.setattr(mod, name, lambda a, inner=inner, name=name: spies.setdefault(name, inner(a)))
    out = tmp_path / "out"
    out.mkdir()

    def route(mod, argv):
        spies.clear()
        try:
            mod.main([*argv, "--out", str(out / "T1.h5ad")])
        except SystemExit as e:
            return ("refused", str(e.code))
        except r51t._Reached as e:
            return ("taken", e.args[0], spies.get("_t1_ms"), spies.get("_t1_candidate"), spies.get("_t1_pooled"),
                    spies.get("_t1_r51"))
        raise AssertionError(f"{argv} neither refused nor read")

    grid = _grid()
    assert len({tuple(x) for x in grid}) == len(grid) == GRID_SETS
    count = {"alike_taken": 0, "alike_refused": 0, "touch": 0}
    new = []
    with open(tmp_path / "stderr.txt", "w") as err, pytest.MonkeyPatch.context() as quiet:
        quiet.setattr(sys, "stderr", err)
        for text, named in ((_ledger_text(), True), (_ledger_text(r66=None), False)):
            ledger.write_text(text)
            for argv in grid:
                a = _parse(argv[1:])
                o, n = route(old, argv), route(bt, argv)
                hit = bt._l12_hit(a)
                if hit is None or hit[0] != "R66" or a.captain_ruling != "R66":  # not R66's set under R66: as before
                    assert n == o, argv
                    count["alike_taken" if o[0] == "taken" else "alike_refused"] += 1
                    continue
                assert o[0] == "refused", argv  # the old builder took none of them
                if named and argv[1:len(RC) + 1] == RC:
                    assert n[:2] == ("taken", "union") and n[-1] == _want(hit[1]), argv
                    new.append(hit[1])
                else:
                    assert n[0] == "refused", argv
                    count["touch"] += 1
    assert count["alike_taken"] and count["alike_refused"] and count["touch"]
    # the grid reaches each of the thirty files once on the recomputed command, under the ledger naming them
    assert sorted(new) == sorted(OPTS)


# ------------------------------------------------------------------ the thirty files, built on synthetic stages


def _without_blood(name: str) -> tuple[str, list[str]]:
    """The file's command without Blood's step: row R66's file, row R60's file or T1_s2rcMS20 (row R51's)."""
    rest = re.sub(r"H\d\d", "", name)
    option = re.sub(r" --held-momentum \S+", "", OPTS[name])
    if rest in OPTS:
        return rest, [*option.split(), "--captain-ruling", "R66"]
    if rest in l12t.OPTS:
        assert l12t.OPTS[rest] == option
        return rest, [*option.split(), "--captain-ruling", "R60"]
    assert rest == "T1_s2rcMS20" and option == MS20_OPT
    return rest, MS20


@pytest.mark.parametrize("name", list(OPTS))
def test_each_files_steps_record_and_blood_rows(bt, ms4, ledger, monkeypatch, tmp_path, name):
    """The steps get the file's values; the record names the row, the slot, the base file and each change; Blood's
    rows move in a file with Blood's step and nothing else of it differs from the same command without that step,
    while a file without it keeps T1_s2rcMS20's Blood rows byte for byte and moves other rows."""
    seen = r51t._calls(bt, monkeypatch)
    f, rec = _build(bt, tmp_path, name, ARGS[name])
    mult, damping, held, ms = _values(name)
    assert seen["ms"] == [{"fraction": mult * 0.4, "seed": SEED, "kw": {"damping": ms}}]
    assert len(seen["g_mom"]) == 1
    g = seen["g_mom"][0]
    assert g["damping"] == damping and g["fraction_mult"] == mult
    assert g["line"] == f"row R66's file {name}, T1's step-on {_slot(name)}, on T1_s2rcMS20's command"
    assert seen["hm"] == ([] if held is None else [{"damping": held, "seed": SEED, "kw": {"base_file": "T1_s2rcMS20"}}])
    cand = rec["t1_candidate"]
    assert {k: cand[k] for k in _want(name)} == _want(name)
    assert cand["applied_fraction"] == mult * cand["calendar_fraction"]
    if "d" in name.removeprefix("T1_s2rcMS20"):
        assert f"damped by {damping:g} in place of B0_DAMPING (row R66's file {name}" in rec["g_mom"]["rule"]
    else:
        assert rec["g_mom"]["damping"] == 0.5 and "damped by B0_DAMPING" in rec["g_mom"]["rule"]
    note = cand["momentum"].get("damping_option")
    if ms != 0.5:
        assert note.startswith(f"--direction-momentum-damping {ms:.1f} (row R66's file {name})") and "R51" not in note
        assert cand["momentum"]["damping"] == ms
    else:
        assert note is None and cand["momentum"]["damping"] == 0.5
    if held is None:
        assert "held_momentum_step" not in cand
        other, other_args = "T1_s2rcMS20", MS20
    else:
        assert cand["held_momentum_step"]["status"] == "shifted" and cand["held_momentum_step"]["damping"] == held
        assert "as T1_s2rcMS20 has it" in cand["held_momentum_step"]["rule"]
        other, other_args = _without_blood(name)
    f0, _ = _build(bt, tmp_path, f"{other}_0", other_args)
    x0, x1 = _same_frame(f0, f)
    labels = ms4[-1]["before"].obs["celltype"].astype(str).to_numpy()  # both builds draw the same cells at the seed
    blood = labels == "Blood"
    assert blood.any() and not blood.all()
    changed = (x0 != x1).any(axis=1)
    if held is None:  # Blood's rows are T1_s2rcMS20's; other rows move
        assert x0[blood].tobytes() == x1[blood].tobytes() and changed.any()
    else:  # every non-Blood row is the file's without Blood's step; a Blood row moves
        assert x0[~blood].tobytes() == x1[~blood].tobytes()
        assert changed.any() and not (changed & ~blood).any()


@pytest.mark.parametrize("name", ["T1_s2rcMS20d15H05", "T1_s2rcMS20d20m10", "T1_s2rcMS20H10m20",
                                  "T1_s2rcMS20d10H10m20", "T1_s2rcMS20d25", "T1_s2rcMS20H20"])
def test_a_file_builds_the_same_bytes_twice_and_is_not_ms20(bt, ms4, ledger, tmp_path, name):
    f0, _ = _build(bt, tmp_path, "T1_s2rcMS20", MS20)
    fa, ra = _build(bt, tmp_path, f"{name}a", ARGS[name])
    fb, rb = _build(bt, tmp_path, f"{name}b", ARGS[name])
    assert fa.read_bytes() == fb.read_bytes() and _record_sha(ra) == _record_sha(rb)
    x0, x1 = _same_frame(f0, fa)
    assert (x0 != x1).any()


def test_the_checker_gives_ms20s_bytes_for_every_file_with_its_steps_at_ms20s_values(bt, ms4, ledger, tmp_path):
    chk = _check()
    f0, _ = _build(bt, tmp_path, "T1_s2rcMS20", MS20)
    before = (bt._t1_l12, bt._t1_ms_step, bt._t1_g_mom, bt.t1_scope.momentum_shift)
    for name in OPTS:
        out = tmp_path / f"{name}_base.h5ad"
        assert chk.main(["base", "t1", *RC, *ARGS[name], "--out", str(out)]) == 0, name
        assert out.read_bytes() == f0.read_bytes(), name
        assert (bt._t1_l12, bt._t1_ms_step, bt._t1_g_mom, bt.t1_scope.momentum_shift) == before  # put back
        rec = json.loads(out.with_suffix(".build.json").read_text())
        assert rec["t1_candidate"]["file"] == name and rec["t1_candidate"]["captain_ruling"] == "R66"  # its own path


def test_the_checkers_settings_are_each_files_steps():
    chk = _check()
    keys = {"d": "g_mom_damping", "H": "held_momentum", "m": "direction_momentum_damping"}
    assert chk.SETTINGS == {f: tuple(keys[s[0]] for s in _steps(f)) for f in OPTS}
    assert chk.BASE_VALUES == {"g_mom_damping": 0.5, "direction_momentum_damping": 0.5}


@pytest.mark.parametrize("name, kept, differs", [
    ("T1_s2rcMS20d15m20", ("direction_momentum_damping",), True),  # g_mom's damping left at 1.5
    ("T1_s2rcMS20d15m20", ("g_mom_damping",), True),  # MS's damping left at 2.0
    ("T1_s2rcMS20H05m10", ("direction_momentum_damping",), False),  # Blood's step not named: its factors still at 1
])
def test_the_checker_fails_when_its_settings_forget_a_step(bt, ms4, ledger, tmp_path, monkeypatch, name, kept,
                                                           differs):
    """A checker whose settings forget one of the file's steps fails its call check; where the forgotten step is left
    at the file's value, the file it builds is not T1_s2rcMS20's bytes either (the launcher's check)."""
    chk = _check()
    f0, _ = _build(bt, tmp_path, "T1_s2rcMS20", MS20)
    monkeypatch.setitem(chk.SETTINGS, name, kept)
    out = tmp_path / "left.h5ad"
    assert chk.main(["base", "t1", *RC, *ARGS[name], "--out", str(out)]) == 1
    assert (out.read_bytes() != f0.read_bytes()) is differs


def test_the_checker_fails_off_row_r66s_route(bt, ms4, ledger, tmp_path):
    chk = _check()
    assert chk.main(["base", "t1", *RC, *MS20, "--out", str(tmp_path / "ms20.h5ad")]) == 1
    with pytest.raises(SystemExit, match="row R66's files only"):  # a line-12 file: its row's route, not R66's
        chk.main(["base", "t1", *RC, *l12t.ARGS["T1_s2rcMS20d10"], "--out", str(tmp_path / "d10.h5ad")])
    with pytest.raises(SystemExit):  # a row R66 set under another ruling: refused as before
        chk.main(["base", "t1", *RC, *ARGS["T1_s2rcMS20d25"][:-1], "R60", "--out", str(tmp_path / "x.h5ad")])
    with pytest.raises(SystemExit, match="usage"):
        chk.main(["msb", "t1"])


@pytest.mark.parametrize("base, args", [("T1_s2rcMS20", MS20), ("T1_s2rcMS20d10", l12t.ARGS["T1_s2rcMS20d10"]),
                                        ("T1_s2rcMS20H05", l12t.ARGS["T1_s2rcMS20H05"]),
                                        ("T1_s2rcMS20m20", l12t.ARGS["T1_s2rcMS20m20"])])
def test_without_the_new_route_ms20_and_line_12s_files_are_the_same_bytes_and_records(bt, ms4, ledger, tmp_path,
                                                                                      monkeypatch, base, args):
    f1, r1 = _build(bt, tmp_path, f"{base}_with", args)
    monkeypatch.setattr(bt, "T1_R66_FILES", {})
    f0, r0 = _build(bt, tmp_path, f"{base}_without", args)
    assert f1.read_bytes() == f0.read_bytes() and _record_sha(r1) == _record_sha(r0)


def test_the_checker_prints_and_deletes_nothing():
    src = (ROOT / "experiments" / "builds" / "t1_r66_check.py").read_text()
    code = "\n".join(line.split("#")[0] for line in src.splitlines())
    assert "print(" not in code
    assert not re.search(r"\b(unlink|rmtree|remove|rmdir)\b", code)
    assert np is not None
