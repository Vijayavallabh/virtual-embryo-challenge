"""`build_test.py t1` under row R73 of experiments/harness/TRIES.md (the captain's order of 2026-11-25, "All four": T1's
second refinement round): four files, each T1_s2rcMSBm10's command (row R51's file: `t1 --g-mom --base lib3
--feng-groups recomputed --feng-fraction-mult 1.5 --direction-momentum all --direction-momentum-damping 1.0`) with one
further change, each taken only with --captain-ruling R73 and each as row R73 names it, the string paired with its file:
  T1_s2rcMSBm10d10  ... --g-mom-damping 1.0 --direction-momentum-damping 1.0 --captain-ruling R73
  T1_s2rcMSBm10H05  ... --held-momentum 0.5 --direction-momentum-damping 1.0 --captain-ruling R73
  T1_s2rcMSBm10d15  ... --g-mom-damping 1.5 --direction-momentum-damping 1.0 --captain-ruling R73
  T1_s2rcMSBm10H10  ... --held-momentum 1.0 --direction-momentum-damping 1.0 --captain-ruling R73
Row R66's route (_t1_l12) on MSB's pair, the strings read by _r66_values with MSB's pair as their head.

Held here on synthetic stages and directions only (the fixtures of tests/test_t1_r51_build.py), with temporary ledgers
holding a synthetic row R73 in place of the repository's ledger; one test reads the repository's row R73 by code and
prints nothing of it: (1) the builder's one table is the row's four (file, option string) pairs in its order, each name
the one its string gives (near-miss names such as T1_s2rcMSBd10m10, T1_s2rcMSBm10m10 and T1_s2rcMSBm10d10H05 refused),
and no set is an earlier row's; (2) each file is taken with row R73 only, every other use under R73 is refused before
anything but the ledger is read, and each set under any other ruling or none is refused as before; (3) the builder at
f7029c5 (main before row R73's route, read with `git show`) and this one route every argument set of a grid alike but
R73's under R73; (4) each file's steps get the file's values, its record names the row, the base file T1_s2rcMSBm10 and
its change (Blood's step's record naming T1_s2rcMSBm10 and R73); an H file's Blood rows are those of row R51's H file at
its value and its other rows T1_s2rcMSBm10's, byte for byte, while a d file keeps T1_s2rcMSBm10's Blood rows; with its
one further step put back at T1_s2rcMSBm10's value (experiments/builds/t1_r73_check.py base) each file is
T1_s2rcMSBm10's bytes; each builds the same bytes twice; MSB's, MSBm10's, row R51's single-step files and row R66's
are the same bytes and records with the new route taken out. No released file, outside-atlas file, board or score is
read."""

from __future__ import annotations

import hashlib
import importlib.util
import itertools
import json
import re
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
for _p in ("experiments/external", "experiments/t1_m1", "experiments/t1_direction", "experiments/t1_hygiene",
           "experiments/builds"):
    sys.path.insert(0, str(ROOT / _p))

import build_test as btm  # noqa: E402  (the module the fixture `bt` returns)
import feng_relabel as fr  # noqa: E402
import test_t1_r51_build as r51t  # noqa: E402
import test_t1_r66_build as r66t  # noqa: E402
from conftest import git_env  # noqa: E402
from test_hidden_counts import t1_env  # noqa: E402, F401  (the synthetic build_t1 fixture)
from test_t1_candidates_build import RC, _build, _no_read, _same_frame  # noqa: E402
from test_t1_fraction2_build import _record_sha  # noqa: E402
from test_t1_lm_m20_build import ROW_R58  # noqa: E402
from test_t1_msb_build import MSB  # noqa: E402
from test_t1_r67_build import ROW_R67  # noqa: E402

bt, ms4 = r51t.bt, r51t.ms4  # the builder module and the three-library recomputed setup

SEED = r51t.SEED
ROW = "R73"
BASE_CMD = "t1 --g-mom --base lib3 --feng-groups recomputed"
MSB_OPT = "--feng-fraction-mult 1.5 --direction-momentum all"
M10_FILE = "T1_s2rcMSBm10"
# the four (file, option string) pairs, written here from row R73's words, in its order
TABLE = [("T1_s2rcMSBm10d10", f"{MSB_OPT} --g-mom-damping 1.0 --direction-momentum-damping 1.0"),
         ("T1_s2rcMSBm10H05", f"{MSB_OPT} --held-momentum 0.5 --direction-momentum-damping 1.0"),
         ("T1_s2rcMSBm10d15", f"{MSB_OPT} --g-mom-damping 1.5 --direction-momentum-damping 1.0"),
         ("T1_s2rcMSBm10H10", f"{MSB_OPT} --held-momentum 1.0 --direction-momentum-damping 1.0")]
OPTS = dict(TABLE)
NAMES = list(OPTS)
ARGS = {f: [*o.split(), "--captain-ruling", ROW] for f, o in TABLE}
# file -> (g_mom's damping, Blood's damping or None), MS's own damping 1.0 in each
VALUES = {"T1_s2rcMSBm10d10": (1.0, None), "T1_s2rcMSBm10H05": (0.5, 0.5), "T1_s2rcMSBm10d15": (1.5, None),
          "T1_s2rcMSBm10H10": (0.5, 1.0)}
# file -> row R51's single-step file at the same value (T1_s2rcMSB's command with that one step)
SINGLE = {"T1_s2rcMSBm10d10": "T1_s2rcMSBd10", "T1_s2rcMSBm10H05": "T1_s2rcMSBH05",
          "T1_s2rcMSBm10d15": "T1_s2rcMSBd15", "T1_s2rcMSBm10H10": "T1_s2rcMSBH10"}
H_NAMES = [f for f in NAMES if VALUES[f][1] is not None]
M10 = r51t.ARGS[M10_FILE]  # row R51's file, MS's own damping 1.0: the base of row R73's
ORDER = ("experiments/harness/TRIES.md row R73, the captain's order of 2026-11-25 (\"All four\", the four files sent "
         "at once; \"ruling, no run\"; decision log 2026-11-25 (bbbb), T1's second refinement round) naming `{option}` "
         "({file}) and `--captain-ruling R73`: file {file}, T1_s2rcMSBm10's command (row R51's file: candidate MS's "
         "step exactly as row R23 fixes it on candidate B's command, the atlas step at 1.5 x the calendar fraction, "
         "g_mom's within-type damping 0.5, MS's own damping 1.0 in place of harness.B0_DAMPING) with one further "
         "change, {change}, the change row R51's file {single} makes to T1_s2rcMSB's command; nothing else changed")
CHANGES = {"T1_s2rcMSBm10d10": "g_mom's within-type damping 1.0", "T1_s2rcMSBm10d15": "g_mom's within-type damping 1.5",
           **{f: (f"Blood's own released E9.5 -> E10.5 change at damping {VALUES[f][1]:.1f} on the drawn Blood cells, "
                  "after MS's step and before g_mom (R1 reversed for Blood in this one file, Macrophage and Early VIC "
                  "held)") for f in H_NAMES}}
NOTE = ("--direction-momentum-damping 1.0 (row R73's file {name}): MS's own damping in place of harness.B0_DAMPING "
        "(0.5); as row R51's file T1_s2rcMSBm10, this file's base, has it; the file's one further change as its order "
        "names it")
LINE = "row R73's file {name}, T1's second refinement round, on T1_s2rcMSBm10's command"
BLOOD_TAIL = ("; Blood's step is row R73's one further change to T1_s2rcMSBm10's command (row R51's file, MS's own "
              "damping 1.0), as row R51's H files make it to T1_s2rcMSB's")
# Row R73's first cell in the ledger's words, less one clause that names no option string (the real row's backticked
# strings and pairs are checked against this copy's by code in test_the_real_rows_first_cell_pairs_the_tables_strings)
FIRST_R73 = ("R73 (the captain's order of 2026-11-25, ruling, no run: \"All four\" (his pick in a dialog whose "
             "question named T1_s2rcMSBm10d10, T1_s2rcMSBm10d15, T1_s2rcMSBm10H05 and T1_s2rcMSBm10H10, all sent at "
             "once, and whose option said this row names the four strings, each only with `--captain-ruling R73`); the "
             "T1 builder is to take on `t1 --g-mom --base lib3 --feng-groups recomputed` exactly these four option "
             "strings, each only with `--captain-ruling R73`: " + ", ".join(f"`{o}` ({f})" for f, o in TABLE)
             + "; decision log 2026-11-25 (bbbb); not a try)")
ROW_R73 = f"| {FIRST_R73} | T1 | refinement |"
# Numbers a refusal may show: row R66's (tests/test_t1_r66_build.py), row R73's ID, the late mix's row, file and stages,
# and rows R58's and R67's IDs and values (their routes refuse a set of theirs under their own ruling)
NUMBERS_ALLOWED = r66t.NUMBERS_ALLOWED | {"73", "54", "16", "17", "58", "67", "0.75", "1.25", "075", "125"}
RULINGS_OTHER = [None, "R51", "R60", "R66", "R67", "R58", "R26", "R23", "R22", "R46", "R730", "r73", "R73a", "R7"]
NEAR_MISS_NAMES = ["T1_s2rcMSBd10m10", "T1_s2rcMSBm10m10", "T1_s2rcMSBm10d10H05"]


def _ledger_text(r73: str | None = ROW_R73) -> str:
    return r66t._ledger_text() + f"{ROW_R58}\n{ROW_R67}\n" + (f"{r73}\n" if r73 else "")


@pytest.fixture
def ledger(bt, tmp_path, monkeypatch):
    """A temporary ledger: tests/test_t1_r66_build.py's (R20, R22, R23, R26, R37, synthetic R44, R45 and R51, copies of
    rows R60's, R46's and R66's first cells) with row R58 as tests/test_t1_lm_m20_build.py copies it, the synthetic
    row R67 of tests/test_t1_r67_build.py and the synthetic row R73 above; the builder reads it in place of the
    repository's."""
    p = tmp_path / "TRIES.md"
    p.write_text(_ledger_text())
    monkeypatch.setattr(bt, "TRIES_FILE", p)
    return p


def _parse(argv: list[str]):
    return btm.build_parser().parse_args(["t1", *argv, "--out", "x.h5ad"])


def _want(name: str) -> dict:
    """_t1_l12's values for a file, written out here from the row's words."""
    damping, held = VALUES[name]
    return {"direction_momentum": "all", "captain_ruling": ROW, "file": name, "base_file": M10_FILE,
            "option": f"`{OPTS[name]}`", "feng_fraction_mult": 1.5, "g_mom_damping": damping, "held_momentum": held,
            "direction_momentum_damping": 1.0,
            "order": ORDER.format(option=OPTS[name], file=name, change=CHANGES[name], single=SINGLE[name]),
            "command": f"{BASE_CMD} {OPTS[name]} --captain-ruling {ROW}"}


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
    import t1_r73_check

    return t1_r73_check


def _pairs(cell: str) -> list[tuple[str, str]]:
    return [(f, o) for o, f in re.findall(r"`([^`]*)` \((T1_[A-Za-z0-9_]+)\)", cell)]


# ------------------------------------------------------------------ the one table and the ruling check


def test_the_builders_one_table_is_the_rows_four_pairs_in_its_order():
    assert [tuple(p) for p in btm.T1_R73_TABLE] == TABLE and _pairs(FIRST_R73) == TABLE
    assert list(btm.T1_R73_FILES) == NAMES and {f: v[0] for f, v in btm.T1_R73_FILES.items()} == OPTS
    assert re.findall(r"`([^`]*)`", FIRST_R73) == ["--captain-ruling R73", BASE_CMD, "--captain-ruling R73",
                                                   *OPTS.values()]
    assert btm.T1_R73_ROW == "R73" and btm.T1_R73_BASE_FILE == M10_FILE and btm.T1_R73_BASE_OPTION == MSB_OPT
    assert btm.T1_R73_SINGLE == SINGLE and btm.T1_R73_CHANGES == CHANGES
    for name, option in TABLE:
        assert btm.T1_R73_COMMAND.format(option=option) == f"{BASE_CMD} {option} --captain-ruling R73"
        damping, held = VALUES[name]
        # the values the route matches on, as the parser gives them for the string
        a = _parse([*RC, *ARGS[name]])
        assert tuple(btm.T1_R73_FILES[name][1:]) == (a.feng_fraction_mult, a.g_mom_damping, a.held_momentum,
                                                    a.direction_momentum_damping)
        assert btm.T1_R73_FILES[name][1:] == (1.5, None if damping == 0.5 else damping,
                                              None if held is None else f"{held:.1f}", "1.0")
        # MSBm10's string with the one step of row R51's file at the same value put in before MS's damping
        single = r51t.OPTS[SINGLE[name]]
        assert option == f"{single} --direction-momentum-damping 1.0" and single != r51t.OPTS[M10_FILE]
        assert btm.T1_R51_CHANGES[SINGLE[name]] == CHANGES[name]


def test_no_set_or_name_of_row_r73_is_an_earlier_rows():
    r73 = {tuple(v[1:]) for v in btm.T1_R73_FILES.values()}
    earlier_sets = ({tuple(v[1:]) for v in btm.T1_R66_FILES.values()}
                    | {tuple(v[1:]) for v in btm.T1_R60_FILES.values()}
                    | {tuple(v[1:]) for v in btm.T1_R51_FILES.values()} | {btm.T1_R46_VALUES})
    assert len(r73) == 4 and not r73 & earlier_sets
    earlier = {*btm.T1_R66_FILES, *btm.T1_R60_FILES, *btm.T1_R51_FILES, *btm.T1_R67_FILES, btm.T1_R46_FILE,
               btm.T1_R58_FILE, btm.T1_R54_FILE, btm.T1_MS25_FILE, *btm.T1_HM_FILES.values()}
    assert not set(NAMES) & earlier
    assert not {o for o in OPTS.values()} & {v[0] for v in (*btm.T1_R66_FILES.values(), *btm.T1_R60_FILES.values(),
                                                           *btm.T1_R51_FILES.values(), *btm.T1_R67_FILES.values())}


def test_the_real_rows_first_cell_pairs_the_tables_strings():
    """Row R73's first cell in the repository's ledger, read by code and printed nowhere: exactly one row R73, its
    first cell pairing the table's four strings with their files in the table's order, naming in backticks what the
    synthetic copy names, the copy its words but for one clause, and passing the builder's own check for each file."""
    text = (ROOT / "experiments" / "harness" / "TRIES.md").read_text(encoding="utf-8")
    cells = [ln[2:].split("|", 1)[0].strip() for ln in text.splitlines() if ln.startswith("| R73 ")]
    ok = (len(cells) == 1 and _pairs(cells[0]) == TABLE
          and re.findall(r"`([^`]*)`", cells[0]) == re.findall(r"`([^`]*)`", FIRST_R73)
          and cells[0].startswith("R73 (the captain's order") and "ruling, no run" in cells[0]
          and all(btm._r73_row_names(f) for f in NAMES))
    assert ok, "the repository's row R73 does not name the table as the copy does (nothing of it printed)"
    kept = [part for part in cells[0].split("; ") if not part.startswith("T1's second refinement round")]
    same = len(kept) == len(cells[0].split("; ")) - 1 and "; ".join(kept) == FIRST_R73
    assert same, "the copy is not the row's first cell less its one clause (nothing of it printed)"


@pytest.mark.parametrize("name, option", [
    *[(bad, OPTS["T1_s2rcMSBm10d10"]) for bad in NEAR_MISS_NAMES],  # the near-miss names on d10's string
    ("T1_s2rcMSBd10m10", OPTS["T1_s2rcMSBm10d10"]),  # row R66's d, H, m order of the same file
    ("T1_s2rcMSBm10m10", f"{MSB_OPT} --direction-momentum-damping 1.0"),  # MSBm10's own string, no further step
    ("T1_s2rcMSBm10d10H05", f"{MSB_OPT} --g-mom-damping 1.0 --held-momentum 0.5 --direction-momentum-damping 1.0"),
    ("T1_s2rcMSBm10", f"{MSB_OPT} --direction-momentum-damping 1.0"),
    ("T1_s2rcMSBm10d1", f"{MSB_OPT} --g-mom-damping 1.0 --direction-momentum-damping 1.0"),
    ("T1_s2rcMSBm10d10", f"{MSB_OPT} --g-mom-damping 1.00 --direction-momentum-damping 1.0"),  # the literal value only
    ("T1_s2rcMSBm10d10", f"{MSB_OPT} --g-mom-damping 1 --direction-momentum-damping 1.0"),
    ("T1_s2rcMSBm10H05", f"{MSB_OPT} --held-momentum .5 --direction-momentum-damping 1.0"),
    ("T1_s2rcMSBm10d10", f"{MSB_OPT} --g-mom-damping 1.0 --direction-momentum-damping 2.0"),  # another MS damping
    ("T1_s2rcMSBm10d10", f"{MSB_OPT} --g-mom-damping 1.0 --direction-momentum-damping 1"),
    ("T1_s2rcMSBm10d10", f"{MSB_OPT} --g-mom-damping 1.0"),  # row R51's file, no MS damping
    ("T1_s2rcMSBm10d10", f"{MSB_OPT} --direction-momentum-damping 1.0 --g-mom-damping 1.0"),  # m before d
    ("T1_s2rcMSBm10d10", "--feng-fraction-mult 2.0 --direction-momentum all --g-mom-damping 1.0 "
                         "--direction-momentum-damping 1.0"),  # row R66's file at 2.0
    ("T1_s2rcMSBm10d10", "--direction-momentum all --feng-fraction-mult 1.5 --g-mom-damping 1.0 "
                         "--direction-momentum-damping 1.0"),
    ("T1_s2rcMSBm10H05", f"{MSB_OPT} --held-momentum 0.5 --direction-momentum-damping 1.0 --herd"),
    ("T1_s2rcMSBm10d10", f"{MSB_OPT} --g-mom-damping 1.0 --direction-momentum-damping 1.0 --direction-momentum-damping "
                         "1.0"),
])
def test_a_table_entry_off_the_rows_form_is_not_read(name, option):
    with pytest.raises(ValueError):
        btm._r73_values(name, option)


def test_r66s_reader_reads_row_r66s_table_as_before_and_r73s_only_on_msbs_pair():
    for option, *values in btm.T1_R66_FILES.values():
        assert btm._r66_values(option) == tuple(values)
    with pytest.raises(ValueError, match="not a row R66 option string"):
        btm._r66_values(OPTS["T1_s2rcMSBm10d10"])  # MSB's pair is not row R66's head
    with pytest.raises(ValueError, match="not a row R73 option string"):
        btm._r66_values(r66t.OPTS["T1_s2rcMS20d10m10"], btm.T1_R73_BASE_OPTION, "R73")


def test_the_ruling_check_takes_only_a_single_row_r73_pairing_the_whole_string_with_its_file(tmp_path):
    p = tmp_path / "TRIES.md"
    p.write_text(f"# ledger\n{ROW_R73}\n")
    assert all(btm._r73_row_names(f, p) for f in NAMES)
    for name, option in TABLE:
        others = [f for f in NAMES if f != name]
        for bad in (ROW_R73.replace(f"`{option}`", option),  # not backticked
                    ROW_R73.replace(f"`{option}`", f"`{option} --herd`"),  # inside a longer string
                    ROW_R73.replace(f"`{option}` ({name})", f"`{option}` ({others[0]})")):  # paired with another file
            p.write_text(f"# ledger\n{bad}\n")
            assert [btm._r73_row_names(f, p) for f in NAMES] == [f != name for f in NAMES], (name, bad[-60:])
    for bad in (ROW_R73.replace("`--captain-ruling R73`", "--captain-ruling R73"),
                ROW_R73.replace("`--captain-ruling R73`", "`--captain-ruling R67`"),
                ROW_R73.replace("ruling, no run", "a try"), ROW_R73.replace("R73 (the captain's order", "R73 (a draft"),
                ROW_R73.replace("| R73 (", "| R730 ("), f"{ROW_R73}\n{ROW_R73}", "| R66 (the captain's order |",
                "| R73 (the captain's order of X, ruling, no run) | T1 | "
                + ", ".join(f"`{o}` ({f})" for f, o in TABLE) + " `--captain-ruling R73` |"):  # outside the first cell
        p.write_text(f"# ledger\n{bad}\n")
        assert not any(btm._r73_row_names(f, p) for f in NAMES), bad[:60]


def test_the_options_are_off_by_default(bt):
    a = _parse(RC)
    assert a.direction_momentum_damping is None and a.captain_ruling is None
    assert bt._t1_l12(a) is None and not bt._l12_ordered(a) and bt._l12_hit(a) is None


# ------------------------------------------------------------------ routing and refusals


@pytest.mark.parametrize("name", NAMES)
def test_each_file_is_taken_with_row_r73_and_nothing_but_the_ledger_is_read(bt, ledger, monkeypatch, name):
    for obj, attr in ((bt, "load_stage"), (bt, "load_boards"), (fr, "committed_at_head")):
        monkeypatch.setattr(obj, attr, _no_read)
    a = _parse([*RC, *ARGS[name]])
    assert bt._l12_hit(a) == ("R73", name)
    assert bt._t1_l12(a) == _want(name) and bt._t1_r51(a) == _want(name) and bt._l12_row_names("R73", name)
    assert bt._l12_ordered(a) and bt._bet_refusals(a) is None and bt._t1_pooled(a) is None
    assert bt._t1_r54(a) is None and bt._t1_r58(a) is None and bt._t1_r67(a) is None and bt._r51_file(a) is None
    for other in (bt._ruled_fraction_ordered, bt._ruled_damping_ordered, bt._ms_ordered, bt._msb_ordered,
                  bt._pooled_ordered, bt._hm_ordered, bt._ms25_ordered, bt._r51_ordered, bt._r54_ordered,
                  bt._r58_ordered, bt._r67_ordered):
        assert not other(a), other.__name__


@pytest.mark.parametrize("name", NAMES)
@pytest.mark.parametrize("ruling", RULINGS_OTHER)
def test_each_set_under_any_other_ruling_is_refused_as_before(bt, ledger, monkeypatch, tmp_path, name, ruling):
    argv = ["t1", *RC, *OPTS[name].split(), *([] if ruling is None else ["--captain-ruling", ruling])]
    a = _parse(argv[1:])
    assert bt._t1_l12(a) is None and not bt._l12_ordered(a)  # routed as before, where it is refused
    _refused(monkeypatch, tmp_path, argv)


@pytest.mark.parametrize("extra", [
    MSB_OPT.split(),  # MSB's own set (row R26's)
    [*MSB_OPT.split(), "--direction-momentum-damping", "1.0"],  # MSBm10's own set (row R51's)
    *[r51t.OPTS[f].split() for f in SINGLE.values()],  # row R51's single-step sets
    [*MSB_OPT.split(), "--g-mom-damping", "1.0", "--held-momentum", "0.5", "--direction-momentum-damping", "1.0"],
    [*MSB_OPT.split(), "--g-mom-damping", "2.0", "--direction-momentum-damping", "1.0"],  # values the row does not name
    [*MSB_OPT.split(), "--g-mom-damping", "2.5", "--direction-momentum-damping", "1.0"],
    [*MSB_OPT.split(), "--held-momentum", "2.0", "--direction-momentum-damping", "1.0"],
    [*MSB_OPT.split(), "--g-mom-damping", "1.0", "--direction-momentum-damping", "2.0"],
    [*MSB_OPT.split(), "--held-momentum", "0.5", "--direction-momentum-damping", "1.25"],
    [*MSB_OPT.split(), "--held-momentum", "0.50", "--direction-momentum-damping", "1.0"],  # the literal strings only
    [*MSB_OPT.split(), "--g-mom-damping", "1.0", "--direction-momentum-damping", "1"],
    [*MSB_OPT.split(), "--g-mom-damping", "1.0", "--direction-momentum-damping", "1.00"],
    [*MSB_OPT.split(), "--direction-momentum-damping", "1.25"],  # row R67's set
    r66t.OPTS["T1_s2rcMS20d10m10"].split(),  # row R66's file at 2.0
    ["--feng-fraction-mult", "1.5", "--g-mom-damping", "1.0", "--direction-momentum-damping", "1.0"],  # without MS
    ["--direction-momentum", "all", "--g-mom-damping", "1.0", "--direction-momentum-damping", "1.0"],
    [],
])
def test_every_other_set_under_r73_is_refused(bt, ledger, monkeypatch, tmp_path, extra):
    argv = ["t1", *RC, *extra, "--captain-ruling", "R73"]
    msg = _refused(monkeypatch, tmp_path, argv, "--captain-ruling R73 orders four files, each with exactly its own")
    assert all(f"{f} `{o}`" in msg for f, o in TABLE)


@pytest.mark.parametrize("name", NAMES)
@pytest.mark.parametrize("cmd, extra, match", [
    (RC, ["--released-bet", "share_trend", "--bet-strength", "0.5"],
     "refused with --released-bet, --bet-strength: it is T1_s2rcMSBm10's command with one change"),
    (RC, ["--bet-strength", "0.5"], "refused with --bet-strength"),
    (RC, ["--no-hold-unbracketed"], "R1's hold"),
    (RC, ["--seed", "1000002"], "changes T1_s2rcMSBm10 only: it needs the command"),
    (RC, ["--board", "T1:val"], "default seed on the test board"),
    (RC, ["--late-mix", "e15_17"], "--late-mix e15_17 is row R54's file T1_s2rcMSB_LM16 alone and needs "
                                   "--captain-ruling R54, not R73"),
    (RC, ["--cm-share", "0.6"], "--feng-groups recomputed builds one file"),
    (RC, ["--herd"], "--feng-groups recomputed builds one file"),
    (RC, ["--base", "pooled"], "--feng-groups recomputed builds one file"),
    (["--g-mom", "--base", "lib3"], [], "changes T1_s2rcMSBm10 only"),
    ([], [], "changes T1_s2rcMSBm10 only"),
    (["--g-mom", "--base", "lib3", "--feng-groups", "none"], [], "default seed on the test board"),
])
def test_every_other_use_of_a_file_under_r73_is_refused(bt, ledger, monkeypatch, tmp_path, name, cmd, extra, match):
    _refused(monkeypatch, tmp_path, ["t1", *cmd, *ARGS[name], *extra], match)


@pytest.mark.parametrize("name", NAMES)
def test_a_file_whose_row_does_not_name_it_is_refused(bt, ledger, monkeypatch, tmp_path, name):
    option, other = OPTS[name], [f for f in NAMES if f != name][0]
    for row in (None, ROW_R73.replace(f"`{option}`", option), ROW_R73.replace(f"`{option}`", f"`{option} --herd`"),
                ROW_R73.replace(f"`{option}` ({name})", f"`{option}` ({other})"),
                ROW_R73.replace("`--captain-ruling R73`", "--captain-ruling R73"),
                ROW_R73.replace("ruling, no run", "a try"), ROW_R73.replace("R73 (the captain's order", "R73 (a draft"),
                f"{ROW_R73}\n{ROW_R73}"):
        ledger.write_text(_ledger_text(r73=row))
        msg = _refused(monkeypatch, tmp_path, ["t1", *RC, *ARGS[name]],
                       f"row R73's file {name} is refused: experiments/harness/TRIES.md has no single row R73 whose "
                       "first cell is the captain's order")
        assert f"`{option}` ({name}) and `--captain-ruling R73`, each as a whole backticked string" in msg


def test_the_bet_refusal_lets_only_the_four_commands_through(bt, ledger):
    def ns(name, **k):
        damping, held = VALUES[name]
        base = dict(released_bet=None, bet_strength=None, captain_ruling=ROW, feng_fraction_mult=1.5,
                    g_mom_damping=None if damping == 0.5 else damping, direction_momentum="all",
                    held_momentum=None if held is None else f"{held:.1f}", direction_momentum_damping="1.0",
                    late_mix=None, feng_groups="recomputed", seed=SEED, board=None, base="lib3")
        return SimpleNamespace(**{**base, **k})

    for name in NAMES:
        assert bt._bet_refusals(ns(name)) is None and bt._l12_ordered(ns(name))
        for bad in (ns(name, captain_ruling="R66"), ns(name, bet_strength=0.5), ns(name, direction_momentum=None),
                    ns(name, feng_fraction_mult=2.0), ns(name, direction_momentum_damping="2.0"),
                    ns(name, direction_momentum_damping=None)):
            assert not bt._l12_ordered(bad)
    ledger.write_text(_ledger_text(r73=None))
    for name in NAMES:
        assert not bt._l12_ordered(ns(name))
        with pytest.raises(SystemExit, match="with --released-bet only"):
            bt._bet_refusals(ns(name))


def test_the_help_names_row_r73s_files():
    import io
    from contextlib import redirect_stdout

    buf = io.StringIO()
    with redirect_stdout(buf), pytest.raises(SystemExit):
        btm.build_parser().parse_args(["t1", "--help"])
    text = re.sub(r"-\s+", "-", " ".join(buf.getvalue().split()))
    for word in ("With --captain-ruling R73 only: row R73's four files", *NAMES, *OPTS.values(),
                 "--captain-ruling R73", "0.5 and 1.0 are row R73's files T1_s2rcMSBm10H05 and T1_s2rcMSBm10H10",
                 "Row R73's four files: the literal value 1.0", "With exactly one of row R73's four option sets"):
        assert word in text, word


# ------------------------------------------------------------------ the builder at f7029c5 and this one, side by side


OLD_COMMIT = "f7029c5"  # main before row R73's route
OLD_SHA = "c6fc34efebdae2ece921afef4db430fcc58d1cc744feb137571219eb5f6708b2"  # its experiments/builds/build_test.py
CMDS = [RC, ["--g-mom", "--base", "lib3"], []]
MULTS = [None, "1.5", "2.0"]
DMS = [None, "all"]
DAMPINGS = [None, "1.0", "1.5"]
HELD = [None, "0.5", "1.0"]
DMD = [None, "1.0", "2.0"]
GRID_RULINGS = [None, "R73", "R51", "R66", "R60", "R67", "R26"]
GRID_SETS = 3 * 3 * 2 * 3 * 3 * 3 * 7


def _old_builder(tmp_path: Path):
    """The builder at OLD_COMMIT, read with `git show` (skipped where that commit is not in the checkout)."""
    env = {**git_env(), "GIT_OPTIONAL_LOCKS": "0"}
    res = subprocess.run(["git", "-C", str(ROOT), "show", f"{OLD_COMMIT}:experiments/builds/build_test.py"],
                         capture_output=True, env=env)
    if res.returncode != 0:
        pytest.skip(f"commit {OLD_COMMIT} is not available in this checkout")
    assert hashlib.sha256(res.stdout).hexdigest() == OLD_SHA
    path = tmp_path / "old" / "build_test_before_r73.py"
    path.parent.mkdir()
    path.write_bytes(res.stdout)
    spec = importlib.util.spec_from_file_location("build_test_before_r73", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _grid() -> list[list[str]]:
    sets = []
    for cmd, mult, dm, damping, held, dmd, ruling in itertools.product(CMDS, MULTS, DMS, DAMPINGS, HELD, DMD,
                                                                        GRID_RULINGS):
        out = ["t1", *cmd]
        for flag, value in (("--feng-fraction-mult", mult), ("--direction-momentum", dm), ("--g-mom-damping", damping),
                            ("--held-momentum", held), ("--direction-momentum-damping", dmd),
                            ("--captain-ruling", ruling)):
            out += [] if value is None else [flag, value]
        sets.append(out)
    return sets


def test_the_builder_at_f7029c5_and_this_one_route_every_set_alike_but_r73s(bt, tmp_path, monkeypatch):
    """Side by side on the grid, under the ledger with row R73 and without it: every set that does not give
    --captain-ruling R73 is refused with the same message, or reaches the same first read with the same candidate, MS,
    pooled and R51-route values, in both builders; the old builder refuses every set under R73, and this one takes each
    of the four on the recomputed command while the ledger names it, refusing every other set under R73. Nothing is
    built: every read is stopped where it would start."""
    old = _old_builder(tmp_path)
    assert not hasattr(old, "T1_R73_FILES") and hasattr(bt, "T1_R73_FILES") and hasattr(old, "T1_R67_FILES")
    ledger = tmp_path / "TRIES.md"
    record = tmp_path / "relabel" / "MANIFEST.json"  # never written: the relabelling's record is the stand-in below
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
        for text, named in ((_ledger_text(), True), (_ledger_text(r73=None), False)):
            ledger.write_text(text)
            for argv in grid:
                a = _parse(argv[1:])
                o, n = route(old, argv), route(bt, argv)
                if a.captain_ruling != "R73":  # not under R73: as before
                    assert n == o, argv
                    count["alike_taken" if o[0] == "taken" else "alike_refused"] += 1
                    continue
                assert o[0] == "refused", argv  # the old builder took none of them
                hit = bt._l12_hit(a)
                if named and hit is not None and hit[0] == "R73" and argv[1:len(RC) + 1] == RC:
                    assert n[:2] == ("taken", "union") and n[2:5] == (None, None, None), argv
                    assert n[-1] == _want(hit[1]), argv
                    new.append(hit[1])
                else:
                    assert n[0] == "refused", argv
                    count["touch"] += 1
    assert count["alike_taken"] and count["alike_refused"] and count["touch"]
    # the grid reaches each of the four files once on the recomputed command, under the ledger naming them
    assert sorted(new) == sorted(NAMES)


# ------------------------------------------------------------------ the four files, built on synthetic stages


@pytest.mark.parametrize("name", NAMES)
def test_each_files_steps_and_record(bt, ms4, ledger, monkeypatch, tmp_path, name):
    """The steps get the file's values; the record names the row, the base file T1_s2rcMSBm10 and the change, in MS's
    note, g_mom's rule (the d files) and Blood's step's rule (the H files)."""
    seen = r51t._calls(bt, monkeypatch)
    _, rec = _build(bt, tmp_path, name, ARGS[name])
    damping, held = VALUES[name]
    assert seen["ms"] == [{"fraction": 1.5 * 0.4, "seed": SEED, "kw": {"damping": 1.0}}]
    assert len(seen["g_mom"]) == 1
    g = seen["g_mom"][0]
    assert g["damping"] == damping and g["fraction_mult"] == 1.5 and g["line"] == LINE.format(name=name)
    assert seen["hm"] == ([] if held is None else
                          [{"damping": held, "seed": SEED, "kw": {"base_file": M10_FILE, "row": "R73"}}])
    cand = rec["t1_candidate"]
    assert {k: cand[k] for k in _want(name)} == _want(name)
    assert cand["applied_fraction"] == 1.5 * cand["calendar_fraction"]
    assert cand["momentum"]["damping"] == 1.0 and cand["momentum"]["damping_option"] == NOTE.format(name=name)
    if held is None:
        assert f"damped by {damping:g} in place of B0_DAMPING ({LINE.format(name=name)})" in rec["g_mom"]["rule"]
        assert "held_momentum_step" not in cand
        return
    assert rec["g_mom"]["damping"] == 0.5 and "damped by B0_DAMPING" in rec["g_mom"]["rule"]
    step = cand["held_momentum_step"]
    assert step["status"] == "shifted" and step["damping"] == held
    # Blood's step's rule, compared by code and printed nowhere: it ends with the row's own words and its own text
    # names the base file T1_s2rcMSBm10, no other file and no ruling row
    rule = step["rule"]
    head = rule.removesuffix(BLOOD_TAIL)
    ok = (rule.endswith(BLOOD_TAIL) and f"as {M10_FILE} has it" in head
          and set(re.findall(r"T1_s2rc\w*", head)) == {M10_FILE} and not set(re.findall(r"\bR\d+\b", head)) - {"R1"})
    assert ok, "Blood's step's rule does not name row R73 and T1_s2rcMSBm10 alone (nothing of it printed)"


@pytest.mark.parametrize("name", NAMES)
def test_an_h_files_blood_rows_are_r51s_h_files_and_its_other_rows_msbm10s(bt, ms4, ledger, tmp_path, name):
    """On the same seed's draw: an H file's Blood rows are those of row R51's H file at its value and every other row
    is T1_s2rcMSBm10's, byte for byte, a Blood row having moved; a d file keeps T1_s2rcMSBm10's Blood rows and moves
    other rows."""
    f, _ = _build(bt, tmp_path, name, ARGS[name])
    labels = ms4[-1]["before"].obs["celltype"].astype(str).to_numpy()
    f0, _ = _build(bt, tmp_path, M10_FILE, M10)
    assert (ms4[-1]["before"].obs["celltype"].astype(str).to_numpy() == labels).all()  # the same cells drawn
    blood = labels == "Blood"
    assert blood.any() and not blood.all()
    x0, x1 = _same_frame(f0, f)
    changed = (x0 != x1).any(axis=1)
    if VALUES[name][1] is None:
        assert x0[blood].tobytes() == x1[blood].tobytes() and changed.any()
        return
    assert x0[~blood].tobytes() == x1[~blood].tobytes() and (changed & blood).any()
    fs, _ = _build(bt, tmp_path, SINGLE[name], r51t.ARGS[SINGLE[name]])
    xs, x2 = _same_frame(fs, f)
    assert xs[blood].tobytes() == x2[blood].tobytes() and (xs[~blood] != x2[~blood]).any()


def test_each_file_builds_the_same_bytes_twice_and_is_none_of_the_files_it_is_built_on(bt, ms4, ledger, tmp_path):
    shas = {}
    for name, args in (("T1_s2rcMSB", MSB), (M10_FILE, M10), *((s, r51t.ARGS[s]) for s in SINGLE.values())):
        f, _ = _build(bt, tmp_path, name, args)
        shas[name] = hashlib.sha256(f.read_bytes()).hexdigest()
    for name in NAMES:
        fa, ra = _build(bt, tmp_path, f"{name}a", ARGS[name])
        fb, rb = _build(bt, tmp_path, f"{name}b", ARGS[name])
        assert fa.read_bytes() == fb.read_bytes() and _record_sha(ra) == _record_sha(rb), name
        shas[name] = hashlib.sha256(fa.read_bytes()).hexdigest()
    assert len(set(shas.values())) == len(shas) == 10


def test_the_checker_gives_msbm10s_bytes_for_each_file_with_its_step_at_msbm10s_value(bt, ms4, ledger, tmp_path):
    chk = _check()
    f0, _ = _build(bt, tmp_path, M10_FILE, M10)
    before = (bt._t1_l12, bt._t1_ms_step, bt._t1_g_mom, bt.t1_scope.momentum_shift)
    for name in NAMES:
        out = tmp_path / f"{name}_base.h5ad"
        assert chk.main(["base", "t1", *RC, *ARGS[name], "--out", str(out)]) == 0, name
        assert out.read_bytes() == f0.read_bytes(), name
        assert (bt._t1_l12, bt._t1_ms_step, bt._t1_g_mom, bt.t1_scope.momentum_shift) == before  # put back
        cand = json.loads(out.with_suffix(".build.json").read_text())["t1_candidate"]
        assert cand["file"] == name and cand["captain_ruling"] == "R73"  # its own path
        assert cand["direction_momentum_damping"] == 1.0 and cand["g_mom_damping"] == 0.5


def test_the_checkers_settings_are_each_files_one_step():
    chk = _check()
    assert chk.SETTINGS == {f: ("held_momentum",) if VALUES[f][1] is not None else ("g_mom_damping",) for f in NAMES}
    assert chk.BASE_VALUES == {"g_mom_damping": 0.5, "direction_momentum_damping": 1.0}
    assert chk.BLOOD_WITHIN == {"blood": ["Blood"]}


@pytest.mark.parametrize("name, differs", [("T1_s2rcMSBm10d15", True), ("T1_s2rcMSBm10H05", False)])
def test_the_checker_fails_when_its_settings_forget_the_step(bt, ms4, ledger, tmp_path, monkeypatch, name, differs):
    """A checker whose settings forget the file's step fails its call check; for a d file the step is left at the
    file's value and the file it builds is not T1_s2rcMSBm10's bytes either (the launcher's check)."""
    chk = _check()
    f0, _ = _build(bt, tmp_path, M10_FILE, M10)
    monkeypatch.setitem(chk.SETTINGS, name, ())
    out = tmp_path / "left.h5ad"
    assert chk.main(["base", "t1", *RC, *ARGS[name], "--out", str(out)]) == 1
    assert (out.read_bytes() != f0.read_bytes()) is differs


def test_the_checker_fails_when_ms_momentum_does_not_run_at_its_damping(bt, ms4, ledger, tmp_path, monkeypatch):
    chk = _check()
    f0, _ = _build(bt, tmp_path, M10_FILE, M10)
    inner = bt._t1_ms_step

    def stubborn(ref, library, directions, fraction, groups, held, seed, **k):
        return inner(ref, library, directions, fraction, groups, held, seed, **{**k, "damping": 0.5})

    monkeypatch.setattr(bt, "_t1_ms_step", stubborn)
    out = tmp_path / "stubborn.h5ad"
    assert chk.main(["base", "t1", *RC, *ARGS["T1_s2rcMSBm10d10"], "--out", str(out)]) == 1
    assert out.read_bytes() != f0.read_bytes()


def test_the_checker_fails_off_row_r73s_route(bt, ms4, ledger, tmp_path):
    chk = _check()
    assert chk.main(["base", "t1", *RC, *M10, "--out", str(tmp_path / "m10.h5ad")]) == 1  # row R51's route
    assert chk.main(["base", "t1", *RC, *MSB, "--out", str(tmp_path / "msb.h5ad")]) == 1  # MSB's route
    with pytest.raises(SystemExit, match="row R73's files only"):  # a row R66 file: its row's route, not R73's
        chk.main(["base", "t1", *RC, *r66t.ARGS["T1_s2rcMS20d10m10"], "--out", str(tmp_path / "d10m10.h5ad")])
    with pytest.raises(SystemExit):  # a row R73 set under another ruling: refused as before
        chk.main(["base", "t1", *RC, *ARGS["T1_s2rcMSBm10d10"][:-1], "R66", "--out", str(tmp_path / "x.h5ad")])
    for argv in ([], ["base"], ["base", "heart"], ["m10", "t1"]):
        with pytest.raises(SystemExit, match="usage"):
            chk.main(argv)


@pytest.mark.parametrize("base, args", [("T1_s2rcMSB", MSB), (M10_FILE, M10),
                                        *((s, r51t.ARGS[s]) for s in SINGLE.values()),
                                        ("T1_s2rcMS20H05m10", r66t.ARGS["T1_s2rcMS20H05m10"])])
def test_without_the_new_route_the_earlier_files_are_the_same_bytes_and_records(bt, ms4, ledger, tmp_path, monkeypatch,
                                                                                base, args):
    f1, r1 = _build(bt, tmp_path, f"{base}_with", args)
    monkeypatch.setattr(bt, "T1_R73_FILES", {})
    f0, r0 = _build(bt, tmp_path, f"{base}_without", args)
    assert f1.read_bytes() == f0.read_bytes() and _record_sha(r1) == _record_sha(r0)
    assert "R73" not in json.dumps(r1)


def test_the_checker_prints_and_deletes_nothing():
    src = (ROOT / "experiments" / "builds" / "t1_r73_check.py").read_text()
    body = src.split('"""', 2)[2]
    code = "\n".join(line.split("#")[0] for line in body.splitlines())
    assert "print(" not in code
    assert not re.search(r"\b(unlink|rmtree|remove|rmdir)\b", code)
    for word in ("write_text", "write_bytes", "upload", "submit", ".build.json", "read_text"):
        assert word not in code, word
    assert np is not None
