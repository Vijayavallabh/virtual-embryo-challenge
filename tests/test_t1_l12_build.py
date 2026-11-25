"""`build_test.py t1` under row R60 of experiments/harness/TRIES.md (the captain's order of 2026-11-25, T1's line 12):
seven files, each one change to T1_s2rcMS20's command (`t1 --g-mom --base lib3 --feng-groups recomputed
--feng-fraction-mult 2.0 --direction-momentum all`, row R51's file), each taken only with --captain-ruling R60 through
_captain_order_names on row R60's first cell; and row R46's pair T1_s2rcB25MSH05 (`--feng-fraction-mult 2.5
--direction-momentum all --held-momentum 0.5`, row R44's Blood step on T1_s2rcB25MS's command), taken only with
--captain-ruling R46 while row R46's first cell is the captain's and names the file:
  T1_s2rcMS20d10, d15, d20   ... --g-mom-damping 1.0 | 1.5 | 2.0
  T1_s2rcMS20H05, H10        ... --held-momentum 0.5 | 1.0 (Blood moves as in row R51's MSBH05 and MSBH10)
  T1_s2rcMS20m10, m20        ... --direction-momentum-damping 1.0 | 2.0

Held here on synthetic stages and directions only (the fixtures of tests/test_t1_r51_build.py), with temporary ledgers
holding copies of rows R60's and R46's real first cells (one test checks the copies against the repository's ledger,
by hash, printing neither): (1) the eight option sets are the brief's, and the real first cells pass the checks; (2)
each file is taken with its own row only, every other use of it under that row is refused before anything but the
ledger is read, and each set under any other ruling or none, and every other set under R60 or R46, is routed as
before and refused there; (3) the builder at main 0d829db (read with `git show`) and this one route every argument
set of a grid alike but the eight sets under their own rows; (4) each file
builds the same bytes twice, its steps get the file's values, and with its one setting at its base's value
(experiments/builds/t1_l12_check.py base) it is its base's bytes (T1_s2rcMS20's, T1_s2rcB25MS's); the H files and the
pair are their bases' but for Blood rows; MS20's and B25MS's bytes and records are the same with the new route taken
out. No released file, outside-atlas file, board or score is read."""

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
import test_t1_r51_build as r51t  # noqa: E402
from conftest import git_env  # noqa: E402
from test_hidden_counts import t1_env  # noqa: E402, F401  (the synthetic build_t1 fixture)
from test_t1_candidates_build import RC, _build, _no_read, _same_frame  # noqa: E402
from test_t1_fraction2_build import _record_sha  # noqa: E402
from test_t1_hm_build import MS25  # noqa: E402

bt, ms4 = r51t.bt, r51t.ms4  # the builder module and the three-library recomputed setup

SEED = r51t.SEED
BASE_CMD = "t1 --g-mom --base lib3 --feng-groups recomputed"
MS20_OPT = "--feng-fraction-mult 2.0 --direction-momentum all"
# file -> its option string, as the brief's table and row R60 have it
OPTS = {"T1_s2rcMS20d10": f"{MS20_OPT} --g-mom-damping 1.0", "T1_s2rcMS20H05": f"{MS20_OPT} --held-momentum 0.5",
        "T1_s2rcMS20d20": f"{MS20_OPT} --g-mom-damping 2.0",
        "T1_s2rcMS20m10": f"{MS20_OPT} --direction-momentum-damping 1.0",
        "T1_s2rcMS20d15": f"{MS20_OPT} --g-mom-damping 1.5", "T1_s2rcMS20H10": f"{MS20_OPT} --held-momentum 1.0",
        "T1_s2rcMS20m20": f"{MS20_OPT} --direction-momentum-damping 2.0"}
PAIR = "T1_s2rcB25MSH05"
PAIR_OPT = "--feng-fraction-mult 2.5 --direction-momentum all --held-momentum 0.5"
ROWS = {**{f: "R60" for f in OPTS}, PAIR: "R46"}
ALL_OPTS = {**OPTS, PAIR: PAIR_OPT}
ARGS = {f: [*o.split(), "--captain-ruling", ROWS[f]] for f, o in ALL_OPTS.items()}
# file -> (multiplier, g_mom's damping, Blood's damping or None, MS's damping)
VALUES = {"T1_s2rcMS20d10": (2.0, 1.0, None, 0.5), "T1_s2rcMS20H05": (2.0, 0.5, 0.5, 0.5),
          "T1_s2rcMS20d20": (2.0, 2.0, None, 0.5), "T1_s2rcMS20m10": (2.0, 0.5, None, 1.0),
          "T1_s2rcMS20d15": (2.0, 1.5, None, 0.5), "T1_s2rcMS20H10": (2.0, 0.5, 1.0, 0.5),
          "T1_s2rcMS20m20": (2.0, 0.5, None, 2.0), PAIR: (2.5, 0.5, 0.5, 0.5)}
H_FILES = ("T1_s2rcMS20H05", "T1_s2rcMS20H10", PAIR)
MS20 = [*MS20_OPT.split(), "--captain-ruling", "R51"]  # row R51's file, the base of row R60's
BASES = {**{f: ("T1_s2rcMS20", MS20) for f in OPTS}, PAIR: ("T1_s2rcB25MS", MS25)}
# Copies of the two rows' first cells as the repository's ledger holds them (test_the_copies_are_the_ledgers_own)
REAL_R60 = (
    "R60 (the captain's order of 2026-11-25, ruling, no run: \"Yes, T1_s2rcMS20d10, T1_s2rcMS20H05, T1_s2rcMS20d20 and "
    "T1_s2rcMS20m10 together, then T1_s2rcMS20d15, T1_s2rcMS20H10 and T1_s2rcMS20m20 as gated, under R22; lift T64's "
    "close and reverse R1 for Blood in MS20H05 and MS20H10; MS20m20 only on R58's G0 pass.\"; T1's line 12, seven "
    "files around the banked file T1_s2rcMS20; the T1 builder is to take "
    + ", ".join(f"`{o}`" for o in OPTS.values())
    + " on `t1 --g-mom --base lib3 --feng-groups recomputed`, each only with `--captain-ruling R60`; decision log "
    "2026-11-25 (i); not a try)")
REAL_R46 = ("R46 (the captain's order of 2026-11-24, ruling, no run: \"Yes, T1_s2rcB25MSH05 if both beat\"; the pair "
            "of T1's lines 8 and 9, named now and coded and built only if both parts beat the banked file; decision "
            "log 2026-11-24 (ss); not a try)")
ROW_R60 = f"| {REAL_R60} | T1 | line 12 |"
ROW_R46 = f"| {REAL_R46} | T1 | the pair |"
NUMBERS_ALLOWED = r51t.NUMBERS_ALLOWED | {"60", "46", "12"}
RULINGS_OTHER = [None, "R51", "R44", "R45", "R22", "R23", "R26", "R20", "R37", "R58", "R600", "r60", "R46a"]


def _sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def _ledger_text(r60: str | None = ROW_R60, r46: str | None = ROW_R46) -> str:
    return r51t._ledger_text() + "".join(f"{row}\n" for row in (r60, r46) if row)


@pytest.fixture
def ledger(bt, tmp_path, monkeypatch):
    """A temporary ledger: tests/test_t1_r51_build.py's (R20, R22, R23, R26, R37, synthetic R44, R45 and R51) with
    copies of rows R60's and R46's real first cells; the builder reads it in place of the repository's."""
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


def _want(name: str) -> dict:
    """_t1_l12's values for a file, written out here from the brief's words."""
    mult, damping, held, ms_damping = VALUES[name]
    if name == PAIR:
        return {"direction_momentum": "all", "captain_ruling": "R46", "file": PAIR, "base_file": "T1_s2rcB25MS",
                "option": f"`{PAIR_OPT}`", "feng_fraction_mult": 2.5, "g_mom_damping": 0.5, "held_momentum": 0.5,
                "direction_momentum_damping": 0.5, "order": btm.T1_R46_ORDER,
                "command": f"{BASE_CMD} {PAIR_OPT} --captain-ruling R46"}
    return {"direction_momentum": "all", "captain_ruling": "R60", "file": name, "base_file": "T1_s2rcMS20",
            "option": f"`{OPTS[name]}`", "feng_fraction_mult": mult, "g_mom_damping": damping, "held_momentum": held,
            "direction_momentum_damping": ms_damping,
            "order": btm.T1_R60_ORDER.format(option=OPTS[name], file=name, change=btm.T1_R60_CHANGES[name]),
            "command": f"{BASE_CMD} {OPTS[name]} --captain-ruling R60"}


def _check():
    import t1_l12_check

    return t1_l12_check


# ------------------------------------------------------------------ the sets and the ruling checks


def test_the_seven_sets_and_the_pair_are_the_briefs():
    assert {f: v[0] for f, v in btm.T1_R60_FILES.items()} == OPTS
    assert list(btm.T1_R60_FILES) == list(OPTS)  # the order of row R60's two waves
    assert btm.T1_R46_OPTION == PAIR_OPT and btm.T1_R46_FILE == PAIR and btm.T1_R46_BASE_FILE == "T1_s2rcB25MS"
    for f in OPTS:
        assert btm.T1_R60_COMMAND.format(option=OPTS[f]) == f"{BASE_CMD} {OPTS[f]} --captain-ruling R60"
    assert btm.T1_R46_COMMAND == f"{BASE_CMD} {PAIR_OPT} --captain-ruling R46"
    a = _parse(RC)
    assert a.direction_momentum_damping is None and a.held_momentum is None and a.captain_ruling is None
    assert btm._t1_l12(a) is None and not btm._l12_ordered(a)


def test_the_copies_are_the_ledgers_own():
    """The copies above are rows R60's and R46's first cells in the repository's ledger (compared by hash only)."""
    got = {}
    for line in (ROOT / "experiments" / "harness" / "TRIES.md").read_text(encoding="utf-8").splitlines():
        for row in ("R60", "R46"):
            if line.startswith(f"| {row} "):
                got.setdefault(row, []).append(hashlib.sha256(line[2:].split("|", 1)[0].strip().encode()).hexdigest())
    want = {row: [hashlib.sha256(cell.encode()).hexdigest()] for row, cell in (("R60", REAL_R60), ("R46", REAL_R46))}
    same = got == want
    assert same, "a copy is not the ledger's first cell (hashes only)"


def test_row_r60s_first_cell_names_each_set_through_captain_order_names(bt, ledger, monkeypatch):
    seen = []
    inner = bt._captain_order_names
    monkeypatch.setattr(bt, "_captain_order_names", lambda row, opt, *a, **k: seen.append((row, opt)) or inner(
        row, opt, *a, **k))
    assert all(bt._l12_row_names("R60", f) for f in OPTS)
    assert seen == [("R60", f"`{o}`") for o in OPTS.values()]


@pytest.mark.parametrize("cell, ok", [
    (REAL_R46, True),
    (REAL_R46.replace("R46 (the captain's order", "R46 (a proposal"), False),
    (REAL_R46.replace("ruling, no run", "a try"), False),
    (REAL_R46.replace("T1_s2rcB25MSH05", "T1_s2rcB25MSH05x"), False),
    (REAL_R46.replace("T1_s2rcB25MSH05", "xT1_s2rcB25MSH05"), False),
    (REAL_R46.replace("T1_s2rcB25MSH05", "T1_s2rcB25MSH10"), False),
    (REAL_R46.replace("R46 (", "R46b ("), False),
])
def test_row_r46s_check_takes_only_the_captains_first_cell_naming_the_file(bt, tmp_path, cell, ok):
    p = tmp_path / "TRIES.md"
    p.write_text(f"# ledger\n| {cell} | T1 |\n")
    assert bt._r46_row_names(p) is ok
    p.write_text(f"# ledger\n| {cell} | T1 |\n| {cell} | T1 |\n")  # two rows with one ID
    assert bt._r46_row_names(p) is False
    p.write_text("# ledger\n| R46 (the captain's order of X, ruling, no run) | T1_s2rcB25MSH05 |\n")
    assert bt._r46_row_names(p) is False  # named outside the first cell only


def test_row_r60s_check_fails_on_a_cell_that_lacks_a_set(bt, tmp_path):
    p = tmp_path / "TRIES.md"
    for o in OPTS.values():
        p.write_text(f"# ledger\n| {REAL_R60.replace(f'`{o}`', f'`{o} --herd`')} | T1 |\n")
        assert [bt._captain_order_names("R60", f"`{x}`", p) for x in OPTS.values()] == [x != o for x in OPTS.values()]


# ------------------------------------------------------------------ routing and refusals


@pytest.mark.parametrize("name", list(ALL_OPTS))
def test_each_file_is_taken_with_its_row_and_nothing_but_the_ledger_is_read(bt, ledger, monkeypatch, name):
    for obj, attr in ((bt, "load_stage"), (bt, "load_boards"), (fr, "committed_at_head")):
        monkeypatch.setattr(obj, attr, _no_read)
    a = _parse([*RC, *ARGS[name]])
    assert bt._t1_l12(a) == _want(name) and bt._t1_r51(a) == _want(name)
    assert bt._l12_ordered(a) and bt._bet_refusals(a) is None and bt._t1_pooled(a) is None


@pytest.mark.parametrize("name", list(ALL_OPTS))
@pytest.mark.parametrize("ruling", RULINGS_OTHER)
def test_each_set_under_any_other_ruling_is_refused(bt, ledger, monkeypatch, tmp_path, name, ruling):
    argv = ["t1", *RC, *ALL_OPTS[name].split(), *([] if ruling is None else ["--captain-ruling", ruling])]
    a = _parse(argv[1:])
    assert bt._t1_l12(a) is None and not bt._l12_ordered(a)  # routed as before, where it is refused
    _refused(monkeypatch, tmp_path, argv)


@pytest.mark.parametrize("name", list(ALL_OPTS))
def test_the_two_rows_are_not_each_others(bt, ledger, monkeypatch, tmp_path, name):
    other = "R46" if ROWS[name] == "R60" else "R60"
    argv = ["t1", *RC, *ALL_OPTS[name].split(), "--captain-ruling", other]
    assert bt._t1_l12(_parse(argv[1:])) is None  # routed as before, where it is refused
    _refused(monkeypatch, tmp_path, argv)


@pytest.mark.parametrize("ruling, extra", [
    ("R60", MS20_OPT.split()),  # MS20's own set: row R51's file, not R60's
    ("R60", [*MS20_OPT.split(), "--g-mom-damping", "1.0", "--held-momentum", "0.5"]),  # two changes
    ("R60", [*MS20_OPT.split(), "--g-mom-damping", "0.5"]),
    ("R60", [*MS20_OPT.split(), "--held-momentum", "2.0"]),
    ("R60", [*MS20_OPT.split(), "--held-momentum", "0.50"]),  # the literal strings only
    ("R60", [*MS20_OPT.split(), "--direction-momentum-damping", "1"]),
    ("R60", [*MS20_OPT.split(), "--direction-momentum-damping", "0.5"]),
    ("R60", ["--feng-fraction-mult", "1.5", "--direction-momentum", "all", "--g-mom-damping", "1.0"]),
    ("R60", ["--feng-fraction-mult", "2.0", "--g-mom-damping", "1.0"]),  # without MS
    ("R60", []),
    ("R46", ["--feng-fraction-mult", "2.5", "--direction-momentum", "all"]),  # B25MS's set
    ("R46", ["--feng-fraction-mult", "2.5", "--direction-momentum", "all", "--held-momentum", "1.0"]),
    ("R46", ["--feng-fraction-mult", "2.5", "--held-momentum", "0.5"]),  # B25H05's set
    ("R46", ["--feng-fraction-mult", "2.0", "--direction-momentum", "all", "--held-momentum", "0.5"]),
    ("R46", [*PAIR_OPT.split(), "--g-mom-damping", "1.0"]),
    ("R46", []),
])
def test_every_other_set_under_r60_or_r46_is_refused(bt, ledger, monkeypatch, tmp_path, ruling, extra):
    argv = ["t1", *RC, *extra, "--captain-ruling", ruling]
    assert bt._t1_l12(_parse(argv[1:])) is None and not bt._l12_ordered(_parse(argv[1:]))  # routed as before
    _refused(monkeypatch, tmp_path, argv)


@pytest.mark.parametrize("name", list(ALL_OPTS))
@pytest.mark.parametrize("cmd, extra, match", [
    (RC, ["--released-bet", "share_trend", "--bet-strength", "0.5"], "refused with --released-bet, --bet-strength"),
    (RC, ["--bet-strength", "0.5"], "refused with --bet-strength"),
    (RC, ["--no-hold-unbracketed"], "R1's hold"),
    (RC, ["--seed", "1000002"], "default seed on the test board"),
    (RC, ["--board", "T1:val"], "default seed on the test board"),
    (["--g-mom", "--base", "lib3", "--feng-groups", "none"], [], "default seed on the test board"),
])
def test_every_other_use_of_a_file_under_its_row_is_refused(bt, ledger, monkeypatch, tmp_path, name, cmd, extra,
                                                            match):
    _refused(monkeypatch, tmp_path, ["t1", *cmd, *ARGS[name], *extra], match)


@pytest.mark.parametrize("name", list(ALL_OPTS))
def test_a_file_whose_row_does_not_name_it_is_refused(bt, ledger, monkeypatch, tmp_path, name):
    if ROWS[name] == "R60":
        ledger.write_text(_ledger_text(r60=ROW_R60.replace(f"`{OPTS[name]}`", "")))
    else:
        ledger.write_text(_ledger_text(r46=ROW_R46.replace(PAIR, "T1_s2rcB25MS")))
    _refused(monkeypatch, tmp_path, ["t1", *RC, *ARGS[name]], "has no single row")
    ledger.write_text(_ledger_text(r60=None, r46=None))
    _refused(monkeypatch, tmp_path, ["t1", *RC, *ARGS[name]], "has no single row")


# ------------------------------------------------------------------ the builder at main and this one, side by side


OLD_COMMIT = "0d829db"  # main before row R60's and R46's routes
OLD_SHA = "b4476d41ad295edf5e11601ddcb6365a58f19002d895efb66450cb35949473d7"  # its experiments/builds/build_test.py
CMDS = [RC, ["--g-mom", "--base", "lib3"], ["--g-mom", "--base", "lib3", "--feng-groups", "none"]]
DAMPINGS = [None, "0.5", "1.0", "1.5", "2.0"]
MULTS = [None, "1.5", "2.0", "2.5"]
HELD = [None, "0.5", "1.0"]
DMS = [None, "all"]
DMD = [None, "1.0", "2.0"]
GRID_RULINGS = [None, "R51", "R44", "R45", "R22", "R26", "R60", "R46"]
GRID_SETS = 3 * 5 * 4 * 3 * 2 * 3 * 8


def _old_builder(tmp_path: Path):
    env = {**git_env(), "GIT_OPTIONAL_LOCKS": "0"}
    res = subprocess.run(["git", "-C", str(ROOT), "show", f"{OLD_COMMIT}:experiments/builds/build_test.py"],
                         capture_output=True, env=env)
    if res.returncode != 0:
        pytest.skip(f"commit {OLD_COMMIT} is not available in this checkout")
    assert hashlib.sha256(res.stdout).hexdigest() == OLD_SHA
    path = tmp_path / "old" / "build_test_before_l12.py"
    path.parent.mkdir()
    path.write_bytes(res.stdout)
    spec = importlib.util.spec_from_file_location("build_test_before_l12", path)
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


def test_the_builder_at_main_and_this_one_route_every_set_alike_but_r60s_and_r46s(bt, tmp_path, monkeypatch):
    """Side by side on the grid, under the ledger with rows R60 and R46 and without them: every set other than one
    of the eight option sets under its own row (R60's seven, R46's pair) is refused with the same
    message, or reaches the same first read with the same candidate, MS and R51-route values, in both builders; the
    old builder refuses each of those, and this one takes each on the recomputed command while the ledger names it,
    refusing it otherwise. Nothing is built: every read is stopped where it would start."""
    old = _old_builder(tmp_path)
    assert not hasattr(old, "_t1_l12") and hasattr(bt, "_t1_l12")
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
        for text, named in ((_ledger_text(), True), (_ledger_text(r60=None, r46=None), False)):
            ledger.write_text(text)
            for argv in grid:
                a = _parse(argv[1:])
                o, n = route(old, argv), route(bt, argv)
                hit = bt._l12_hit(a)
                if hit is None or hit[0] != a.captain_ruling:  # no file's set under its own row: as before
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
    # the grid reaches each of the eight files once on the recomputed command, under the ledger naming them
    assert sorted(new) == sorted(ALL_OPTS)


# ------------------------------------------------------------------ the eight files, built on synthetic stages


@pytest.mark.parametrize("name", list(ALL_OPTS))
def test_each_file_builds_the_same_bytes_twice_and_is_not_its_base(bt, ms4, ledger, tmp_path, name):
    base, base_args = BASES[name]
    f0, _ = _build(bt, tmp_path, base, base_args)
    fa, ra = _build(bt, tmp_path, f"{name}a", ARGS[name])
    fb, rb = _build(bt, tmp_path, f"{name}b", ARGS[name])
    assert fa.read_bytes() == fb.read_bytes() and _record_sha(ra) == _record_sha(rb)
    x0, x1 = _same_frame(f0, fa)  # the same cells, order, obs and var as the base's file
    assert (x0 != x1).any()
    assert ra["t1_candidate"]["file"] == name and ra["t1_candidate"]["captain_ruling"] == ROWS[name]


def test_the_eight_files_are_eight_files(bt, ms4, ledger, tmp_path):
    files = {name: _build(bt, tmp_path, name, ARGS[name])[0].read_bytes() for name in ALL_OPTS}
    assert len(set(files.values())) == 8


@pytest.mark.parametrize("name", list(ALL_OPTS))
def test_each_step_gets_the_files_values_and_nothing_else(bt, ms4, ledger, monkeypatch, tmp_path, name):
    seen = r51t._calls(bt, monkeypatch)
    _build(bt, tmp_path, name, ARGS[name])
    mult, damping, held, ms_damping = VALUES[name]
    assert seen["ms"] == [{"fraction": mult * 0.4, "seed": SEED, "kw": {"damping": ms_damping}}]
    assert len(seen["g_mom"]) == 1
    g = seen["g_mom"][0]
    assert g["damping"] == damping and g["fraction_mult"] == mult
    if ROWS[name] == "R60":
        assert g["line"] == f"row R60's file {name}, T1's line 12, on T1_s2rcMS20's command"
    if held is None:
        assert seen["hm"] == []
    else:
        assert seen["hm"] == [{"damping": held, "seed": SEED, "kw": {"base_file": BASES[name][0]}}]


@pytest.mark.parametrize("name", list(ALL_OPTS))
def test_the_checker_gives_the_bases_bytes_with_each_files_setting_at_the_bases_value(bt, ms4, ledger, tmp_path,
                                                                                         name):
    chk = _check()
    base, base_args = BASES[name]
    f0, _ = _build(bt, tmp_path, base, base_args)
    out = tmp_path / f"{name}_base.h5ad"
    before = (bt._t1_l12, bt._t1_ms_step, bt._t1_g_mom, bt.t1_scope.momentum_shift)
    assert chk.main(["base", "t1", *RC, *ARGS[name], "--out", str(out)]) == 0
    assert out.read_bytes() == f0.read_bytes()
    assert (bt._t1_l12, bt._t1_ms_step, bt._t1_g_mom, bt.t1_scope.momentum_shift) == before  # put back
    rec = json.loads(out.with_suffix(".build.json").read_text())
    assert rec["t1_candidate"]["file"] == name  # the file's own path, its row's route


def test_the_checker_fails_off_the_new_route(bt, ms4, ledger, tmp_path):
    chk = _check()
    assert chk.main(["base", "t1", *RC, *MS20, "--out", str(tmp_path / "ms20.h5ad")]) == 1
    with pytest.raises(SystemExit):
        chk.main(["base", "t1", *RC, *ARGS["T1_s2rcMS20d10"][:-1], "R51", "--out", str(tmp_path / "x.h5ad")])
    with pytest.raises(SystemExit, match="usage"):
        chk.main(["msb", "t1"])


@pytest.mark.parametrize("name", H_FILES)
def test_the_h_files_and_the_pair_are_their_base_but_for_blood_rows(bt, ms4, ledger, tmp_path, name):
    base, base_args = BASES[name]
    f0, _ = _build(bt, tmp_path, base, base_args)
    fh, rh = _build(bt, tmp_path, name, ARGS[name])
    x0, xh = _same_frame(f0, fh)
    labels = ms4[-1]["before"].obs["celltype"].astype(str).to_numpy()  # the H build's draw: its base's at the seed
    blood = labels == "Blood"
    assert blood.any() and not blood.all()
    assert x0[~blood].tobytes() == xh[~blood].tobytes()  # every non-Blood row: the base's bytes
    changed = (x0 != xh).any(axis=1)
    assert changed.any() and not (changed & ~blood).any()
    assert rh["t1_candidate"]["held_momentum_step"]["status"] == "shifted"
    assert f"as {base} has it" in rh["t1_candidate"]["held_momentum_step"]["rule"]


def test_the_m_files_records_name_row_r60s_file(bt, ms4, ledger, tmp_path):
    for name in ("T1_s2rcMS20m10", "T1_s2rcMS20m20"):
        _, rec = _build(bt, tmp_path, name, ARGS[name])
        note = rec["t1_candidate"]["momentum"]["damping_option"]
        assert note.startswith(f"--direction-momentum-damping {VALUES[name][3]:.1f} (row R60's file {name})")
        assert "R51" not in note and rec["t1_candidate"]["momentum"]["damping"] == VALUES[name][3]
    _, rec = _build(bt, tmp_path, "T1_s2rcMS20d20", ARGS["T1_s2rcMS20d20"])
    assert "row R60's file T1_s2rcMS20d20, T1's line 12" in rec["g_mom"]["rule"]


@pytest.mark.parametrize("base, args", [("T1_s2rcMS20", MS20), ("T1_s2rcB25MS", MS25),
                                        ("T1_s2rcMSBm10", r51t.ARGS["T1_s2rcMSBm10"]),
                                        ("T1_s2rcMSBd10", r51t.ARGS["T1_s2rcMSBd10"])])
def test_without_the_new_route_the_earlier_files_are_the_same_bytes_and_records(bt, ms4, ledger, tmp_path,
                                                                                monkeypatch, base, args):
    f1, r1 = _build(bt, tmp_path, f"{base}_with", args)
    monkeypatch.setattr(bt, "_t1_l12", lambda a: None)
    monkeypatch.setattr(bt, "_l12_ordered", lambda a: False)
    f0, r0 = _build(bt, tmp_path, f"{base}_without", args)
    assert f1.read_bytes() == f0.read_bytes() and _record_sha(r1) == _record_sha(r0)


def test_the_checker_prints_and_deletes_nothing():
    src = (ROOT / "experiments" / "builds" / "t1_l12_check.py").read_text()
    code = "\n".join(line.split("#")[0] for line in src.splitlines())
    assert "print(" not in code
    assert not re.search(r"\b(unlink|rmtree|remove|rmdir)\b", code)
    assert np is not None
