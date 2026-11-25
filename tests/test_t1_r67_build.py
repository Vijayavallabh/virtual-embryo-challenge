"""`build_test.py t1` under row R67 of experiments/harness/TRIES.md (the captain's order of 2026-11-25): three files,
each T1_s2rcMSB's command (row R26's: `t1 --g-mom --base lib3 --feng-groups recomputed --feng-fraction-mult 1.5
--direction-momentum all`) with one change, MS's own damping, on row R58's route (T1_s2rcMSBm20, the same command at
2.0), each taken only with --captain-ruling R67 and each as row R67 names it:
  T1_s2rcMSBm125  ... --direction-momentum-damping 1.25 --captain-ruling R67
  T1_s2rcMSBm15   ... --direction-momentum-damping 1.5 --captain-ruling R67
  T1_s2rcMSBm075  ... --direction-momentum-damping 0.75 --captain-ruling R67

Held here on synthetic stages and directions only (the fixtures of tests/test_t1_r51_build.py), with temporary ledgers
holding a synthetic row R67 in place of the repository's ledger; one test reads the repository's row R67 by code and
prints nothing of it: (1) the builder's one table is the row's three (file, option string) pairs in its order, each name
the one its value gives, and the synthetic row names what the real one names; (2) each file is taken with row R67 only,
every other use under R67 is refused before anything but the ledger is read, and each set under any other ruling or
none is refused as before; (3) the builder at main e5d7e71 (read with `git show`) and this one route every argument
set of a grid alike but R67's under R67; (4) each file's steps are MSB's with MS's damping at the file's value, its
record is MSB's plus its candidate field with the file's own literal note; with MS's damping put back at B0's each file
is MSB's bytes (experiments/builds/t1_r67_check.py base); each builds the same bytes twice and the three are none of
MSB, MSBm10 and MSBm20. No released file, outside-atlas file, board or score is read."""

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

import pytest

ROOT = Path(__file__).resolve().parents[1]
for _p in ("experiments/external", "experiments/t1_m1", "experiments/t1_direction", "experiments/t1_hygiene",
           "experiments/builds", "experiments/t1_scope"):
    sys.path.insert(0, str(ROOT / _p))

import build_test as _bt  # noqa: E402  (the module itself, for the helpers)
import feng_relabel as fr  # noqa: E402
import test_t1_candidates_build as tcb  # noqa: E402
import test_t1_feng_groups_none  # noqa: E402, F401  (imported before any fixture, as tests/test_feng_relabel.py does)
import test_t1_r51_build as r51t  # noqa: E402
from conftest import git_env  # noqa: E402
from test_feng_relabel import _recomputed_env  # noqa: E402
from test_hidden_counts import t1_env  # noqa: E402, F401  (the synthetic build_t1 fixture)
from test_t1_candidates_build import RC, _build, _no_read, _spy  # noqa: E402
from test_t1_fraction2_build import _record_sha  # noqa: E402
from test_t1_lm_m20_build import ROW_R58  # noqa: E402
from test_t1_msb_build import MSB, _three_libraries  # noqa: E402

from vec import harness  # noqa: E402

bt = tcb.bt  # the builder module's fixture

SEED = harness.BUILD_SEED_OFFSET + 1
ROW = "R67"
PAIR = ["--feng-fraction-mult", "1.5", "--direction-momentum", "all"]
BASE_CMD = "t1 --g-mom --base lib3 --feng-groups recomputed"
MSB_OPT = "--feng-fraction-mult 1.5 --direction-momentum all"
# the three (file, option string) pairs, written here from row R67's words, in its order
TABLE = [("T1_s2rcMSBm125", f"{MSB_OPT} --direction-momentum-damping 1.25"),
         ("T1_s2rcMSBm15", f"{MSB_OPT} --direction-momentum-damping 1.5"),
         ("T1_s2rcMSBm075", f"{MSB_OPT} --direction-momentum-damping 0.75")]
OPTS = dict(TABLE)
NAMES = list(OPTS)
VALUES = {"T1_s2rcMSBm125": "1.25", "T1_s2rcMSBm15": "1.5", "T1_s2rcMSBm075": "0.75"}
ARGS = {f: [*o.split(), "--captain-ruling", ROW] for f, o in TABLE}
M10 = r51t.ARGS["T1_s2rcMSBm10"]  # row R51's file, MS's damping 1.0
M20 = [*PAIR, "--direction-momentum-damping", "2.0", "--captain-ruling", "R58"]  # row R58's file, MS's damping 2.0
# each file's note in MS's record, written out
NOTES = {
    "T1_s2rcMSBm125": ("--direction-momentum-damping 1.25 (row R67's file T1_s2rcMSBm125): MS's own damping in place "
                       "of harness.B0_DAMPING (0.5); nothing else changed"),
    "T1_s2rcMSBm15": ("--direction-momentum-damping 1.5 (row R67's file T1_s2rcMSBm15): MS's own damping in place of "
                      "harness.B0_DAMPING (0.5); nothing else changed"),
    "T1_s2rcMSBm075": ("--direction-momentum-damping 0.75 (row R67's file T1_s2rcMSBm075): MS's own damping in place "
                       "of harness.B0_DAMPING (0.5); nothing else changed"),
}
ORDER = ("experiments/harness/TRIES.md row R67, the captain's order of 2026-11-25 (\"Yes to T1_s2rcMSBm125, "
         "T1_s2rcMSBm15 and T1_s2rcMSBm075 under R22, sent at once, in one ruling row: MSB's command with MS's own "
         "damping at 1.25, 1.5 and 0.75.\"; decision log 2026-11-25 (kkk)) naming `{option}` and `--captain-ruling "
         "R67`: file {file}, T1_s2rcMSB's command with one change, MS's own damping {value} in place of "
         "harness.B0_DAMPING; nothing else changed")
# Row R67's first cell in the ledger's words, less one clause that names no option string (the real row's backticked
# strings and pairs are checked against this copy's by code in test_the_real_rows_first_cell_pairs_the_tables_strings)
FIRST_R67 = ("R67 (the captain's order of 2026-11-25, ruling, no run: \"Yes to T1_s2rcMSBm125, T1_s2rcMSBm15 and "
             "T1_s2rcMSBm075 under R22, sent at once, in one ruling row: MSB's command with MS's own damping at 1.25, "
             "1.5 and 0.75.\"; the T1 builder is to take on `t1 --g-mom --base lib3 --feng-groups recomputed` exactly "
             "these three option strings, each only with `--captain-ruling R67`: "
             + ", ".join(f"`{o}` ({f})" for f, o in TABLE) + "; decision log 2026-11-25 (kkk); not a try)")
ROW_R67 = f"| {FIRST_R67} | T1 | damping |"
# Numbers a refusal may show: the options' values, lib3, the rows' IDs, the files' names and the late mix's stages
NUMBERS_ALLOWED = {"0.5", "0.75", "1.0", "1.25", "1.5", "2.0", "1", "2", "3", "10", "15", "16", "17", "20", "075",
                   "125", "13", "26", "51", "54", "58", "67"}
RULINGS_OTHER = [None, "R51", "R58", "R26", "R23", "R54", "R60", "R66", "R22", "R670", "r67", "R67a", "R6"]


def _ledger_text(r67: str | None = ROW_R67) -> str:
    return r51t._ledger_text() + f"{ROW_R58}\n" + (f"{r67}\n" if r67 else "")


@pytest.fixture
def ledger(bt, tmp_path, monkeypatch):
    """A temporary ledger: tests/test_t1_r51_build.py's (R20, R22, R23, R26, R37, synthetic R44, R45 and R51) with row
    R58 as tests/test_t1_lm_m20_build.py copies it and the synthetic row R67 above; the builder reads it in place of
    the repository's."""
    p = tmp_path / "TRIES.md"
    p.write_text(_ledger_text())
    monkeypatch.setattr(bt, "TRIES_FILE", p)
    return p


@pytest.fixture
def ms4(bt, t1_env, monkeypatch):  # noqa: F811
    _recomputed_env(bt, t1_env, monkeypatch)
    _three_libraries(t1_env)
    return _spy(bt, monkeypatch)


def _sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def _parse(argv: list[str]):
    return _bt.build_parser().parse_args([*argv, "--out", "x.h5ad"])


def _want(name: str) -> dict:
    """_t1_r67's values for a file, written out here from the row's words."""
    return {"direction_momentum": "all", "captain_ruling": ROW, "file": name, "base_file": "T1_s2rcMSB",
            "option": f"`{OPTS[name]}`", "feng_fraction_mult": 1.5, "g_mom_damping": 0.5, "held_momentum": None,
            "direction_momentum_damping": float(VALUES[name]),
            "order": ORDER.format(option=OPTS[name], file=name, value=VALUES[name]),
            "command": f"{BASE_CMD} {OPTS[name]} --captain-ruling {ROW}"}


def _refused(monkeypatch, tmp_path, argv: list[str], match: str | None = None) -> str:
    """Refused before any stage, board, union, Feng file, relabelling record or step is read (the ledger only), with
    no number in the message but the allowed ones and what the command line said."""
    for obj, name in ((_bt, "load_stage"), (_bt, "load_boards"), (_bt.hygiene, "load_union"),
                      (_bt.feng_dirs, "feng_directions"), (fr, "committed_at_head"), (_bt, "_t1_hm_step"),
                      (_bt, "_t1_ms_step")):
        monkeypatch.setattr(obj, name, _no_read)
    name = "T1_dry.h5ad" if "--board" in argv else "T1.h5ad"
    with pytest.raises(SystemExit, match=match) as e:
        _bt.main([*argv, "--out", str(tmp_path / name)])
    msg = str(e.value)
    given = set(re.findall(r"\d+(?:\.\d+)?", " ".join(argv)))
    assert set(re.findall(r"\d+(?:\.\d+)?", msg)) <= NUMBERS_ALLOWED | given, msg
    return msg


def _check():
    import t1_r67_check

    return t1_r67_check


def _pairs(cell: str) -> list[tuple[str, str]]:
    return [(f, o) for o, f in re.findall(r"`([^`]*)` \((T1_[A-Za-z0-9_]+)\)", cell)]


# ------------------------------------------------------------------ the one table and the ruling check


def test_the_builders_one_table_is_the_rows_three_pairs_in_its_order():
    assert [tuple(p) for p in _bt.T1_R67_TABLE] == TABLE and _pairs(FIRST_R67) == TABLE
    assert {f: v for f, v in _bt.T1_R67_FILES.items()} == {f: (o, VALUES[f]) for f, o in TABLE}
    assert re.findall(r"`([^`]*)`", FIRST_R67) == [BASE_CMD, "--captain-ruling R67", *OPTS.values()]
    assert _bt.T1_R67_ROW == "R67" and _bt.R67_STACKED == _bt.R58_STACKED and _bt.T1_R67_NOTES == NOTES
    for name, option in TABLE:
        assert option.split()[-1] == VALUES[name] and name == "T1_s2rcMSBm" + VALUES[name].replace(".", "")
        assert _bt.T1_R67_COMMAND.format(option=option) == f"{BASE_CMD} {option} --captain-ruling R67"
    # no file of row R67 is an earlier file, and no value is one another route already takes
    earlier = {*_bt.T1_R51_FILES, _bt.T1_R58_FILE, _bt.T1_R54_FILE, *_bt.T1_R60_FILES, *_bt.T1_R66_FILES}
    assert not set(NAMES) & earlier and not {float(v) for v in VALUES.values()} & {0.5, 1.0, 2.0}


def test_the_real_rows_first_cell_pairs_the_tables_strings():
    """Row R67's first cell in the repository's ledger, read by code and printed nowhere: exactly one row R67, its
    first cell pairing the table's three strings with their files in the table's order, naming in backticks what the
    synthetic copy names, and passing the builder's own check for each file."""
    text = (ROOT / "experiments" / "harness" / "TRIES.md").read_text(encoding="utf-8")
    cells = [ln[2:].split("|", 1)[0] for ln in text.splitlines() if ln.startswith("| R67 ")]
    ok = (len(cells) == 1 and _pairs(cells[0]) == TABLE
          and re.findall(r"`([^`]*)`", cells[0]) == re.findall(r"`([^`]*)`", FIRST_R67)
          and cells[0].startswith("R67 (the captain's order") and "ruling, no run" in cells[0]
          and all(_bt._r67_row_names(f) for f in NAMES))
    assert ok, "the repository's row R67 does not name the table as the copy does (nothing of it printed)"


@pytest.mark.parametrize("name, option", [
    ("T1_s2rcMSBm125", f"{MSB_OPT} --direction-momentum-damping 1.250"),  # the name no longer the value's
    ("T1_s2rcMSBm12", f"{MSB_OPT} --direction-momentum-damping 1.25"),
    ("T1_s2rcMSBm10", f"{MSB_OPT} --direction-momentum-damping 1.0"),  # row R51's value
    ("T1_s2rcMSBm20", f"{MSB_OPT} --direction-momentum-damping 2.0"),  # row R58's value
    ("T1_s2rcMSBm05", f"{MSB_OPT} --direction-momentum-damping 0.5"),  # B0's value
    ("T1_s2rcMSBm1", f"{MSB_OPT} --direction-momentum-damping 1"),  # no point
    ("T1_s2rcMSBm75", f"{MSB_OPT} --direction-momentum-damping .75"),
    ("T1_s2rcMSBm125", f"{MSB_OPT} --direction-momentum-damping 1.25 --herd"),
    ("T1_s2rcMSBm125", "--feng-fraction-mult 2.0 --direction-momentum all --direction-momentum-damping 1.25"),
    ("T1_s2rcMSBm125", "--direction-momentum all --feng-fraction-mult 1.5 --direction-momentum-damping 1.25"),
    ("T1_s2rcMSBm125", f"{MSB_OPT} --g-mom-damping 1.25"),
])
def test_a_table_entry_off_the_rows_form_is_not_read(name, option):
    with pytest.raises(ValueError):
        _bt._r67_value(name, option)


def test_the_ruling_check_takes_only_a_single_row_r67_pairing_the_whole_string_with_its_file(tmp_path):
    p = tmp_path / "TRIES.md"
    p.write_text(f"# ledger\n{ROW_R67}\n")
    assert all(_bt._r67_row_names(f, p) for f in NAMES)
    for name, option in TABLE:
        others = [f for f in NAMES if f != name]
        for bad in (ROW_R67.replace(f"`{option}`", option),  # not backticked
                    ROW_R67.replace(f"`{option}`", f"`{option} --herd`"),  # inside a longer string
                    ROW_R67.replace(f"`{option}` ({name})", f"`{option}` ({others[0]})")):  # paired with another file
            p.write_text(f"# ledger\n{bad}\n")
            assert [_bt._r67_row_names(f, p) for f in NAMES] == [f != name for f in NAMES], (name, bad[-60:])
    for bad in (ROW_R67.replace("`--captain-ruling R67`", "--captain-ruling R67"),
                ROW_R67.replace("`--captain-ruling R67`", "`--captain-ruling R58`"),
                ROW_R67.replace("ruling, no run", "a try"), ROW_R67.replace("R67 (the captain's order", "R67 (a draft"),
                ROW_R67.replace("| R67 (", "| R670 ("), f"{ROW_R67}\n{ROW_R67}", "| R66 (the captain's order |",
                "| R67 (the captain's order of X, ruling, no run) | T1 | "
                + ", ".join(f"`{o}` ({f})" for f, o in TABLE) + " `--captain-ruling R67` |"):  # outside the first cell

        p.write_text(f"# ledger\n{bad}\n")
        assert not any(_bt._r67_row_names(f, p) for f in NAMES), bad[:60]


def test_the_options_are_off_by_default(bt):
    a = _parse(["t1", *RC])
    assert a.direction_momentum_damping is None and a.captain_ruling is None
    assert bt._t1_r67(a) is None and not bt._r67_ordered(a) and bt._r67_file(a) is None


# ------------------------------------------------------------------ routing and refusals


@pytest.mark.parametrize("name", NAMES)
def test_each_file_is_taken_with_its_row_and_routes_nowhere_else(bt, ledger, monkeypatch, name):
    for obj, attr in ((bt, "load_stage"), (bt, "load_boards")):
        monkeypatch.setattr(obj, attr, _no_read)
    monkeypatch.setattr(fr, "committed_at_head", lambda p: json.dumps({"route": "recomputed"}).encode())
    a = _parse(["t1", *RC, *ARGS[name]])
    assert bt._r67_file(a) == name and bt._t1_r67(a) == _want(name) and bt._r67_ordered(a)
    assert bt._t1_r58(a) is None and bt._t1_r54(a) is None and bt._t1_l12(a) is None
    assert bt._bet_refusals(a) is None and bt._feng_recomputed_refusals(a) is None and bt._t1_pooled(a) is None
    for other in (bt._ruled_fraction_ordered, bt._ruled_damping_ordered, bt._ms_ordered, bt._msb_ordered,
                  bt._pooled_ordered, bt._hm_ordered, bt._ms25_ordered, bt._r51_ordered, bt._r54_ordered,
                  bt._r58_ordered, bt._l12_ordered):
        assert not other(a), other.__name__


@pytest.mark.parametrize("name", NAMES)
@pytest.mark.parametrize("ruling", RULINGS_OTHER)
def test_each_set_under_any_other_ruling_is_refused_as_before(bt, ledger, monkeypatch, tmp_path, name, ruling):
    argv = ["t1", *RC, *OPTS[name].split(), *([] if ruling is None else ["--captain-ruling", ruling])]
    a = _parse(argv)
    assert bt._t1_r67(a) is None and not bt._r67_ordered(a)  # routed as before, where it is refused
    match = ("--captain-ruling R58 orders one file, T1_s2rcMSBm20" if ruling == "R58" else
             f"--direction-momentum-damping {VALUES[name]} is refused: it takes 1.0 only, and only in row R51's file")
    _refused(monkeypatch, tmp_path, argv, match)


@pytest.mark.parametrize("argv, match", [
    ([*PAIR, "--captain-ruling", "R67"], "--captain-ruling R67 orders three files, each with exactly its own option"),
    *[([*PAIR, "--direction-momentum-damping", v, "--captain-ruling", "R67"], "orders three files")
      for v in ("0.5", "1.0", "2.0", "1.250", "1.50", "0.750", ".75", "1.2", "0.8", "1", "1.75", "2.5")],
    (["--feng-fraction-mult", "2.0", "--direction-momentum", "all", "--direction-momentum-damping", "1.25",
      "--captain-ruling", "R67"], "orders three files"),
    (["--feng-fraction-mult", "1.5", "--direction-momentum-damping", "1.5", "--captain-ruling", "R67"],
     "orders three files"),
    (["--direction-momentum", "all", "--direction-momentum-damping", "0.75", "--captain-ruling", "R67"],
     "orders three files"),
    (["--captain-ruling", "R67"], "orders three files"),
    (["--g-mom-damping", "1.0", "--captain-ruling", "R67"], "orders three files"),
])
def test_every_other_set_under_r67_is_refused(bt, ledger, monkeypatch, tmp_path, argv, match):
    msg = _refused(monkeypatch, tmp_path, ["t1", *RC, *argv], match)
    assert all(f"{f} `{o}`" in msg for f, o in TABLE)


@pytest.mark.parametrize("name", NAMES)
@pytest.mark.parametrize("cmd, extra, match", [
    (RC, ["--g-mom-damping", "1.0"], "refused with --g-mom-damping: it is T1_s2rcMSB's command with one change"),
    (RC, ["--held-momentum", "0.5"], "refused with --held-momentum"),
    (RC, ["--released-bet", "share_trend", "--bet-strength", "0.5"], "refused with --released-bet, --bet-strength"),
    (RC, ["--bet-strength", "0.5"], "refused with --bet-strength"),
    (RC, ["--late-mix", "e15_17"], "--late-mix e15_17 is row R54's file T1_s2rcMSB_LM16 alone and needs "
                                   "--captain-ruling R54, not R67"),
    (RC, ["--no-hold-unbracketed"], "refused with --no-hold-unbracketed"),
    (RC, ["--seed", "1000002"], "changes T1_s2rcMSB only: it needs the command"),
    (RC, ["--board", "T1:val"], "changes T1_s2rcMSB only"),
    (RC, ["--cm-share", "0.6"], "--feng-groups recomputed builds one file"),
    (RC, ["--herd"], "--feng-groups recomputed builds one file"),
    (RC, ["--base", "pooled"], "--feng-groups recomputed builds one file"),
    (["--g-mom", "--base", "lib3"], [], "changes T1_s2rcMSB only"),
    ([], [], "changes T1_s2rcMSB only"),
    (["--g-mom", "--base", "lib3", "--feng-groups", "none"], [], "changes T1_s2rcMSB only"),
])
def test_every_other_use_of_a_file_under_r67_is_refused(bt, ledger, monkeypatch, tmp_path, name, cmd, extra, match):
    _refused(monkeypatch, tmp_path, ["t1", *cmd, *ARGS[name], *extra], match)


@pytest.mark.parametrize("name", NAMES)
def test_a_file_whose_row_does_not_name_it_is_refused(bt, ledger, monkeypatch, tmp_path, name):
    option, other = OPTS[name], [f for f in NAMES if f != name][0]
    for row in (None, ROW_R67.replace(f"`{option}`", option), ROW_R67.replace(f"`{option}`", f"`{option} --herd`"),
                ROW_R67.replace(f"`{option}` ({name})", f"`{option}` ({other})"),
                ROW_R67.replace("`--captain-ruling R67`", "--captain-ruling R67"),
                ROW_R67.replace("ruling, no run", "a try"), ROW_R67.replace("R67 (the captain's order", "R67 (a draft"),
                f"{ROW_R67}\n{ROW_R67}"):
        ledger.write_text(_ledger_text(r67=row))
        _refused(monkeypatch, tmp_path, ["t1", *RC, *ARGS[name]],
                 f"row R67's file {name} is refused: experiments/harness/TRIES.md has no single row R67 whose first "
                 "cell is the captain's order")


def test_the_bet_refusal_lets_only_the_three_commands_through(bt, ledger):
    def ns(name, **k):
        base = dict(released_bet=None, bet_strength=None, captain_ruling=ROW, feng_fraction_mult=1.5,
                    g_mom_damping=None, direction_momentum="all", held_momentum=None,
                    direction_momentum_damping=VALUES[name], late_mix=None, feng_groups="recomputed", seed=SEED,
                    board=None, base="lib3")
        return SimpleNamespace(**{**base, **k})

    for name in NAMES:
        assert bt._bet_refusals(ns(name)) is None and bt._r67_ordered(ns(name))
        for bad in (ns(name, captain_ruling="R58"), ns(name, bet_strength=0.5), ns(name, direction_momentum=None),
                    ns(name, feng_fraction_mult=2.0), ns(name, g_mom_damping=1.0), ns(name, held_momentum="0.5"),
                    ns(name, late_mix="e15_17"), ns(name, direction_momentum_damping="1.0")):
            assert not bt._r67_ordered(bad)
    ledger.write_text(_ledger_text(r67=None))
    for name in NAMES:
        assert not bt._r67_ordered(ns(name))
        with pytest.raises(SystemExit, match="with --released-bet only"):
            bt._bet_refusals(ns(name))


def test_the_help_names_row_r67s_files():
    import io
    from contextlib import redirect_stdout

    buf = io.StringIO()
    with redirect_stdout(buf), pytest.raises(SystemExit):
        _bt.build_parser().parse_args(["t1", "--help"])
    text = re.sub(r"-\s+", "-", " ".join(buf.getvalue().split()))
    for word in ("Row R67's three files", *NAMES, *OPTS.values(), "--captain-ruling R67",
                 *(f"{f} the literal value {v}" for f, v in VALUES.items())):
        assert word in text, word


# ------------------------------------------------------------------ the builder at e5d7e71 and this one, side by side


OLD_COMMIT = "e5d7e71"  # main before row R67's route
OLD_SHA = "a4eab7a61e40e032dbd2f437b297447b83ea5cbfd4f0e9e8f3c38e9c7c6a8331"  # its experiments/builds/build_test.py
CMDS = [RC, ["--g-mom", "--base", "lib3"], ["--g-mom", "--base", "pooled", "--feng-groups", "recomputed"], []]
MULTS = [None, "1.5", "2.0"]
DMS = [None, "all"]
DMD = [None, "0.75", "1.0", "1.25", "1.5", "2.0"]
GRID_RULINGS = [None, "R67", "R58", "R51", "R26", "R66", "R22"]
EXTRAS = [[], ["--seed", "1000002"], ["--g-mom-damping", "1.0"], ["--held-momentum", "0.5"]]
GRID_SETS = 4 * 3 * 2 * 6 * 7 * 4


def _old_builder(tmp_path: Path):
    """The builder at OLD_COMMIT, read with `git show` (skipped where that commit is not in the checkout)."""
    env = {**git_env(), "GIT_OPTIONAL_LOCKS": "0"}
    res = subprocess.run(["git", "-C", str(ROOT), "show", f"{OLD_COMMIT}:experiments/builds/build_test.py"],
                         capture_output=True, env=env)
    if res.returncode != 0:
        pytest.skip(f"commit {OLD_COMMIT} is not available in this checkout")
    assert hashlib.sha256(res.stdout).hexdigest() == OLD_SHA
    path = tmp_path / "old" / "build_test_before_r67.py"
    path.parent.mkdir()
    path.write_bytes(res.stdout)
    spec = importlib.util.spec_from_file_location("build_test_before_r67", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _grid() -> list[list[str]]:
    out = []
    for cmd, mult, dm, dmd, ruling, extra in itertools.product(CMDS, MULTS, DMS, DMD, GRID_RULINGS, EXTRAS):
        argv = ["t1", *cmd]
        for flag, value in (("--feng-fraction-mult", mult), ("--direction-momentum", dm),
                            ("--direction-momentum-damping", dmd), ("--captain-ruling", ruling)):
            argv += [] if value is None else [flag, value]
        out.append([*argv, *extra])
    return out


def _kind(argv: list[str], named: bool) -> str:
    """The file's name for exactly one of row R67's three commands under a ledger naming the row, "touch" for any
    other set giving --captain-ruling R67, "other" for every other set."""
    tail = argv[len(RC) + 1:] if argv[:len(RC) + 1] == ["t1", *RC] else None
    for name, args in ARGS.items():
        if tail == args and named:
            return name
    if "--captain-ruling" in argv and argv[argv.index("--captain-ruling") + 1] == ROW:
        return "touch"
    return "other"


def test_the_builder_at_e5d7e71_and_this_one_route_every_set_alike_but_r67s(bt, tmp_path, monkeypatch):
    """Side by side on the grid, under the ledger with row R67 and without it: every set that does not give
    --captain-ruling R67 is refused with the same message, or reaches the same first read with the same candidate, MS,
    pooled and R51-route values, in both builders, and this one's R67 route returns None for it, reading nothing; the
    old builder refuses every set under R67, and this one takes each of the three on the recomputed command while the
    ledger names it, refusing every other set under R67. Nothing is built: every read is stopped where it would
    start."""
    old = _old_builder(tmp_path)
    assert not hasattr(old, "_t1_r67") and hasattr(bt, "_t1_r67") and hasattr(old, "_t1_r58")
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
        new = ("_t1_r67",) if mod is bt else ()
        for name in ("_t1_ms", "_t1_candidate", "_t1_pooled", "_t1_r51", *new):
            inner = getattr(mod, name)
            monkeypatch.setattr(mod, name, lambda a, inner=inner, name=name: spies.setdefault(name, inner(a)))
    out = tmp_path / "out"
    out.mkdir()

    def route(mod, argv):
        spies.clear()
        name = "T1_dry.h5ad" if "--board" in argv else "T1.h5ad"
        try:
            mod.main([*argv, "--out", str(out / name)])
        except SystemExit as e:
            return ("refused", str(e.code)), dict(spies)
        except r51t._Reached as e:
            return (("taken", e.args[0], spies.get("_t1_ms"), spies.get("_t1_candidate"), spies.get("_t1_pooled"),
                     spies.get("_t1_r51")), dict(spies))
        raise AssertionError(f"{argv} neither refused nor read")

    grid = _grid()
    assert len({tuple(x) for x in grid}) == len(grid) == GRID_SETS
    count = {"other": 0, "touch": 0, "taken_alike": 0, "refused_alike": 0}
    new_taken = []
    with open(tmp_path / "stderr.txt", "w") as err, pytest.MonkeyPatch.context() as quiet:
        quiet.setattr(sys, "stderr", err)
        for text, named in ((_ledger_text(), True), (_ledger_text(r67=None), False)):
            ledger.write_text(text)
            for argv in grid:
                kind = _kind(argv, named)
                (o, _), (n, seen) = route(old, argv), route(bt, argv)
                if kind == "other":
                    count["other"] += 1
                    assert n == o, argv
                    assert seen.get("_t1_r67") is None, argv
                    count["taken_alike" if o[0] == "taken" else "refused_alike"] += 1
                    continue
                assert o[0] == "refused", argv  # the old builder took none of these
                if kind == "touch":
                    count["touch"] += 1
                    assert n[0] == "refused", argv
                    continue
                assert n == ("taken", "union", None, None, None, None), argv
                assert seen["_t1_r67"] == _want(kind), argv
                new_taken.append(kind)
    assert count["taken_alike"] and count["refused_alike"] and count["touch"]
    assert count["other"] + count["touch"] + len(new_taken) == 2 * GRID_SETS
    assert sorted(new_taken) == sorted(NAMES)  # each of the three once, on the recomputed command, under the row


# ------------------------------------------------------------------ the three files, built on synthetic stages


@pytest.mark.parametrize("name", NAMES)
def test_each_files_steps_are_msbs_but_ms_damping(bt, ms4, ledger, tmp_path, monkeypatch, name):
    seen = r51t._calls(bt, monkeypatch)
    _build(bt, tmp_path, "T1_s2rcMSB", MSB)
    _build(bt, tmp_path, name, ARGS[name])
    assert seen["ms"][0] == {"fraction": 1.5 * 0.4, "seed": SEED, "kw": {}}
    assert seen["ms"][1] == {"fraction": 1.5 * 0.4, "seed": SEED,
                             "kw": {"damping": float(VALUES[name]), "option_note": NOTES[name]}}
    assert seen["g_mom"][0] == seen["g_mom"][1] and seen["g_mom"][1]["damping"] == 0.5
    assert "line" not in seen["g_mom"][1] and seen["hm"] == []


@pytest.mark.parametrize("name", NAMES)
def test_each_files_record_is_msbs_plus_its_candidate_field_with_its_own_note(bt, ms4, ledger, tmp_path, name):
    _, r0 = _build(bt, tmp_path, "T1_s2rcMSB", MSB)
    _, r1 = _build(bt, tmp_path, name, ARGS[name])
    skip = {"t1_candidate", "g_mom"}
    assert {k: v for k, v in r1.items() if k not in skip} == {k: v for k, v in r0.items() if k not in skip}
    g0, g1 = dict(r0["g_mom"]), dict(r1["g_mom"])
    g0.pop("median_depth_scale"), g1.pop("median_depth_scale")
    assert g0 == g1
    cand = r1["t1_candidate"]
    assert {k: cand[k] for k in _want(name)} == _want(name)
    assert set(cand) == {*_want(name), "momentum", "calendar_fraction", "applied_fraction"}
    assert cand["applied_fraction"] == 1.5 * cand["calendar_fraction"]
    m0, m1 = dict(r0["t1_candidate"]["momentum"]), dict(cand["momentum"])
    assert m1.pop("damping") == float(VALUES[name]) and m0.pop("damping") == 0.5
    assert m1.pop("damping_option") == NOTES[name] and m1 == m0
    assert "R67" not in json.dumps(r0) and "R51" not in NOTES[name] and "late_mix" not in r1


def test_each_file_builds_the_same_bytes_twice_and_is_none_of_msb_m10_m20_or_another(bt, ms4, ledger, tmp_path):
    shas = {}
    for name, args in (("T1_s2rcMSB", MSB), ("T1_s2rcMSBm10", M10), ("T1_s2rcMSBm20", M20)):
        f, _ = _build(bt, tmp_path, name, args)
        shas[name] = _sha(f)
    for name in NAMES:
        fa, ra = _build(bt, tmp_path, f"{name}a", ARGS[name])
        fb, rb = _build(bt, tmp_path, f"{name}b", ARGS[name])
        assert fa.read_bytes() == fb.read_bytes() and _record_sha(ra) == _record_sha(rb), name
        shas[name] = _sha(fa)
    assert len(set(shas.values())) == 6


@pytest.mark.parametrize("name", NAMES)
def test_with_ms_damping_at_b0s_each_file_is_msbs_bytes_by_the_checker(bt, ms4, ledger, tmp_path, name):
    chk = _check()
    f0, _ = _build(bt, tmp_path, "T1_s2rcMSB", MSB)
    out = tmp_path / f"{name}_base.h5ad"
    before = (bt._t1_r67, bt._t1_ms_step, bt._t1_g_mom, bt.t1_scope.momentum_shift)
    assert chk.main(["base", "t1", *RC, *ARGS[name], "--out", str(out)]) == 0
    assert out.read_bytes() == f0.read_bytes()
    assert (bt._t1_r67, bt._t1_ms_step, bt._t1_g_mom, bt.t1_scope.momentum_shift) == before  # put back
    rec = json.loads(out.with_suffix(".build.json").read_text())
    cand = rec["t1_candidate"]
    assert cand["file"] == name and cand["captain_ruling"] == "R67"  # its own path
    assert cand["direction_momentum_damping"] == 0.5 and "damping_option" not in cand["momentum"]


def test_the_checker_fails_off_row_r67s_route(bt, ms4, ledger, tmp_path):
    chk = _check()
    assert chk.main(["base", "t1", *RC, *MSB, "--out", str(tmp_path / "msb.h5ad")]) == 1  # MSB's route
    assert chk.main(["base", "t1", *RC, *M20, "--out", str(tmp_path / "m20.h5ad")]) == 1  # row R58's route
    assert chk.main(["base", "t1", *RC, *M10, "--out", str(tmp_path / "m10.h5ad")]) == 1  # row R51's route
    with pytest.raises(SystemExit):  # a row R67 set under another ruling: refused as before
        chk.main(["base", "t1", *RC, *ARGS["T1_s2rcMSBm125"][:-1], "R58", "--out", str(tmp_path / "x.h5ad")])
    for argv in ([], ["base"], ["base", "heart"], ["msb", "t1"]):
        with pytest.raises(SystemExit, match="usage"):
            chk.main(argv)


def test_the_checker_fails_when_ms_momentum_does_not_run_at_b0s_damping(bt, ms4, ledger, tmp_path, monkeypatch):
    """A step that does not pass on the damping it is given (here: always the file's) fails the checker's momentum
    check, and the file it builds is not MSB's bytes either (the launcher's check)."""
    chk = _check()
    f0, _ = _build(bt, tmp_path, "T1_s2rcMSB", MSB)
    inner = bt._t1_ms_step

    def stubborn(ref, library, directions, fraction, groups, held, seed, **k):
        return inner(ref, library, directions, fraction, groups, held, seed, **{**k, "damping": 1.25})

    monkeypatch.setattr(bt, "_t1_ms_step", stubborn)
    out = tmp_path / "stubborn.h5ad"
    assert chk.main(["base", "t1", *RC, *ARGS["T1_s2rcMSBm125"], "--out", str(out)]) == 1
    assert out.read_bytes() != f0.read_bytes()


def test_the_checker_prints_and_deletes_nothing():
    src = (ROOT / "experiments" / "builds" / "t1_r67_check.py").read_text()
    body = src.split('"""', 2)[2]
    code = "\n".join(line.split("#")[0] for line in body.splitlines())
    assert "print(" not in code
    assert not re.search(r"\b(unlink|rmtree|remove|rmdir)\b", code)
    for word in ("write_text", "write_bytes", "upload", "submit", ".build.json", "read_text"):
        assert word not in code, word
