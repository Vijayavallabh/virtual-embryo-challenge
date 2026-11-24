"""`build_test.py t1` under row R51 of experiments/harness/TRIES.md (the captain's order of 2026-11-24, T1's line 10):
seven files, each one change to T1_s2rcMSB's command (row R26's: `t1 --g-mom --base lib3 --feng-groups recomputed
--feng-fraction-mult 1.5 --direction-momentum all`), each taken only with --captain-ruling R51 while the ledger holds
exactly one row R51 whose first cell is the captain's ruling naming the file's option string and `--captain-ruling R51`
as whole backticked strings:
  T1_s2rcMS20, T1_s2rcMS30      --feng-fraction-mult 2.0 | 3.0 --direction-momentum all
  T1_s2rcMSBd10, T1_s2rcMSBd15  --feng-fraction-mult 1.5 --direction-momentum all --g-mom-damping 1.0 | 1.5
  T1_s2rcMSBH05, T1_s2rcMSBH10  --feng-fraction-mult 1.5 --direction-momentum all --held-momentum 0.5 | 1.0
  T1_s2rcMSBm10                 --feng-fraction-mult 1.5 --direction-momentum all --direction-momentum-damping 1.0

Held here on synthetic stages and directions only (the fixtures of tests/test_t1_hm_build.py), with temporary ledgers
in place of the repository's (no test here reads the repository's own ledger): (1) nothing that builds today changes:
every existing command (R13, A, B, C, R20, R22's lines 1 to 7, R23, R26, R37, R44's three files and R45's) gives the
bytes and records pinned for it, R44's and R45's pins made by the builder at main 32c8c09 (before these options
existed), twice in separate processes with the same result; (2) with each new setting at MSB's value the new path gives
MSB's bytes (experiments/builds/t1_r51_check.py msb, and by hand); MSBH05's and MSBH10's non-Blood rows are MSB's at the
same seed byte for byte and a Blood row moves; with Blood's factors at 1 an H file is MSB's bytes; MSBm10 with MS's
damping at the default is MSB's bytes; with every momentum factor at 1 MS20, MS30, MSBd10 and MSBm10 are the no-MS
files of their commands (R20's, R22's 3.0, candidate C's, candidate B's); (3) each file is taken only with its row,
every other use is refused before anything but the ledger is read, and the new option is refused on every other
command; (4) the builder at main 32c8c09 (read with `git show`) and this one route every other argument set of a grid
exactly alike against synthetic ledgers; (5) the checker's two modes pass on these files and fail on altered ones,
printing booleans only. No released file, outside-atlas file, board or score is read."""

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

import anndata as ad
import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
for _p in ("experiments/external", "experiments/t1_m1", "experiments/t1_direction", "experiments/t1_hygiene",
           "experiments/builds"):
    sys.path.insert(0, str(ROOT / _p))

import feng_relabel as fr  # noqa: E402
import test_t1_candidates_build as tcb  # noqa: E402
import test_t1_feng_groups_none  # noqa: E402, F401  (imported before any fixture, as tests/test_feng_relabel.py does)
from conftest import git_env  # noqa: E402
from test_feng_relabel import _recomputed_env  # noqa: E402
from test_hidden_counts import t1_env  # noqa: E402, F401  (the synthetic build_t1 fixture)
from test_t1_candidates_build import RC, THREE, _build, _dense, _no_read, _same_frame, _spy  # noqa: E402
from test_t1_damping_line_build import DAMP as LINE2  # noqa: E402
from test_t1_damping_line_build import LINE1_BEFORE  # noqa: E402
from test_t1_fraction2_build import BEFORE, R20, _record_sha  # noqa: E402
from test_t1_fraction_line_build import R20_BEFORE  # noqa: E402
from test_t1_fraction_line_build import R22 as LINE1  # noqa: E402
from test_t1_hm_build import (  # noqa: E402
    BANKED,
    BANKED_MS4_BEFORE,
    DEFAULT_MS4_BEFORE,
    HM,
    LINE6_BEFORE,
    LINE7_BEFORE,
    MS25,
    POOLED_BEFORE,
    ROW_H,
    ROW_M,
)
from test_t1_hm_build import NAMES as HM_NAMES  # noqa: E402
from test_t1_l67pb_build import LINE5_BEFORE, LINE6, LINE7, MSB_BEFORE, ROW_R37  # noqa: E402
from test_t1_l67pb_build import ROW_R22 as ROW_R22_L67  # noqa: E402
from test_t1_line3_build import LINE2_BEFORE, LINE3, ROW_R20  # noqa: E402
from test_t1_line5_build import LINE5  # noqa: E402
from test_t1_ms_build import LINE3_BEFORE, ROW_R23  # noqa: E402
from test_t1_ms_build import MS as MS_ARGS  # noqa: E402
from test_t1_msb_build import MS_BEFORE, MSB, ROW_R26, _three_libraries  # noqa: E402

from vec import harness  # noqa: E402

bt, rc4 = tcb.bt, tcb.rc4  # the candidates' fixtures (the builder module; the two-library recomputed setup)

SEED = harness.BUILD_SEED_OFFSET + 1
ROW = "R51"
# file -> its option string, as row R51 names it (the brief's and the design's words)
OPTS = {"T1_s2rcMS20": "--feng-fraction-mult 2.0 --direction-momentum all",
        "T1_s2rcMS30": "--feng-fraction-mult 3.0 --direction-momentum all",
        "T1_s2rcMSBd10": "--feng-fraction-mult 1.5 --direction-momentum all --g-mom-damping 1.0",
        "T1_s2rcMSBd15": "--feng-fraction-mult 1.5 --direction-momentum all --g-mom-damping 1.5",
        "T1_s2rcMSBH05": "--feng-fraction-mult 1.5 --direction-momentum all --held-momentum 0.5",
        "T1_s2rcMSBH10": "--feng-fraction-mult 1.5 --direction-momentum all --held-momentum 1.0",
        "T1_s2rcMSBm10": "--feng-fraction-mult 1.5 --direction-momentum all --direction-momentum-damping 1.0"}
ARGS = {name: [*opt.split(), "--captain-ruling", ROW] for name, opt in OPTS.items()}
# file -> (multiplier, g_mom's damping, Blood's damping or None, MS's damping)
VALUES = {"T1_s2rcMS20": (2.0, 0.5, None, 0.5), "T1_s2rcMS30": (3.0, 0.5, None, 0.5),
          "T1_s2rcMSBd10": (1.5, 1.0, None, 0.5), "T1_s2rcMSBd15": (1.5, 1.5, None, 0.5),
          "T1_s2rcMSBH05": (1.5, 0.5, 0.5, 0.5), "T1_s2rcMSBH10": (1.5, 0.5, 1.0, 0.5),
          "T1_s2rcMSBm10": (1.5, 0.5, None, 1.0)}
H_FILES = ("T1_s2rcMSBH05", "T1_s2rcMSBH10")
BASE_CMD = "t1 --g-mom --base lib3 --feng-groups recomputed"
# The ruling row's first cell, in the words the ledger's first cells use, inside a synthetic row of that form
FIRST_CELL_R51 = ("R51 (the captain's order of a later day, ruling, no run: \"Yes, the seven files under R22 as "
                  "gated\"; T1's line 10, seven files around the combined-momentum file; the T1 builder is to take "
                  + ", ".join(f"`{o}`" for o in OPTS.values()) + f" on `{BASE_CMD}`, each only with "
                  "`--captain-ruling R51`; not a try)")
ROW_R51 = f"| {FIRST_CELL_R51} | T1 | line 10 |"
# Row R44's three files and row R45's file on the three-library fixture (ms4 below), built by the builder at main
# 32c8c09 (before row R51's options existed), twice in separate processes with the same result: the file's sha256 and
# its record's (feng_hygiene.commit set to "own_repo"), the same kind of synthetic-fixture pin as the earlier lines'.
HM_BEFORE = {0.5: ("442e304f9aaefceddfcc7c67c9c4fa81b3822022eb9dfea71cfde3996621572c",
                   "2874fe4fb57e4de306823c98e3ea94f722b7a842b021da82fd8b553b717f8765"),
             1.0: ("3331f5911be4bc290431cf5a7d0ddd18eda513a04c9acc463044595a6f0b8610",
                   "60cbed062be8f1687450574ebec56974d93f184ad1ce1d697d760b1d849d1025"),
             2.0: ("52b9e6e121d945f5299d562128e5659f322355217969c74667ceca96af1208dd",
                   "51773142ec1fcaed22d0e6bb3332116066a80d8ffbdfce095c634978fce5c80d")}
MS25_BEFORE = ("4c4ab9e3174855cb6ff5b4e34e4be942f0c1ec6916d040fae5b411fd71bf5cb1",
               "b82d4a693f70b432dc42939de7cec7f08acea92b08a297b69a334dfd635fa61e")
# Numbers a refusal may show: the options' values, lib3, the libraries, the rows' IDs and the files' names (T1_s2rcMS20
# reads as 1, 2 and 20; T1_s2rcMSBH05 as 1, 2 and 05), besides what the command line itself said
NUMBERS_ALLOWED = {"0.5", "1.0", "1.5", "2.0", "2.5", "3.0", "1", "2", "3", "05", "10", "15", "20", "25", "30", "13",
                   "22", "23", "26", "37", "44", "45", "51"}


def _sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def _ledger_text(r51: str | None = ROW_R51) -> str:
    rows = (ROW_R20, ROW_R22_L67, ROW_R23, ROW_R26, ROW_R37, ROW_H, ROW_M, r51)
    return "# ledger\n" + "".join(f"{row}\n" for row in rows if row)


@pytest.fixture
def ledger(bt, tmp_path, monkeypatch):
    """A temporary ledger holding the ruling rows R20, R22 (lines 1 to 7), R23, R26, R37, the synthetic rows R44 and
    R45 of tests/test_t1_hm_build.py and a synthetic row R51 naming the seven option sets; the builder reads it in
    place of the repository's."""
    p = tmp_path / "TRIES.md"
    p.write_text(_ledger_text())
    monkeypatch.setattr(bt, "TRIES_FILE", p)
    return p


@pytest.fixture
def ms4(bt, t1_env, monkeypatch):  # noqa: F811
    """The four-group recomputed setup on E10.5 names ending in _1, _2 and _3 in turn (MS's fit needs libraries 1 and
    2); returns build_submission's spy."""
    _recomputed_env(bt, t1_env, monkeypatch)
    _three_libraries(t1_env)
    return _spy(bt, monkeypatch)


@pytest.fixture
def ms3(bt, t1_env, monkeypatch):  # noqa: F811
    _recomputed_env(bt, t1_env, monkeypatch, groups=THREE)
    _three_libraries(t1_env)
    return _spy(bt, monkeypatch)


def _parse(bt, argv: list[str]):
    return bt.build_parser().parse_args([*argv, "--out", "x.h5ad"])


def _labels(pop) -> np.ndarray:
    return pop.obs["celltype"].astype(str).to_numpy()


def _want(bt, name: str) -> dict:
    """_t1_r51's values for a file, written out here from the brief's words."""
    mult, damping, held, ms_damping = VALUES[name]
    changes = {"T1_s2rcMS20": "the atlas step at 2.0 x the calendar fraction",
               "T1_s2rcMS30": "the atlas step at 3.0 x the calendar fraction",
               "T1_s2rcMSBd10": "g_mom's within-type damping 1.0", "T1_s2rcMSBd15": "g_mom's within-type damping 1.5",
               "T1_s2rcMSBm10": "MS's own damping 1.0 in place of harness.B0_DAMPING"}
    change = changes.get(name) or (f"Blood's own released E9.5 -> E10.5 change at damping {held:.1f} on the drawn "
                                   "Blood cells, after MS's step and before g_mom (R1 reversed for Blood in this one "
                                   "file, Macrophage and Early VIC held)")
    return {"direction_momentum": "all", "captain_ruling": ROW, "file": name, "base_file": "T1_s2rcMSB",
            "option": f"`{OPTS[name]}`", "feng_fraction_mult": mult, "g_mom_damping": damping,
            "held_momentum": held, "direction_momentum_damping": ms_damping,
            "order": bt.T1_R51_ORDER.format(option=OPTS[name], file=name, change=change),
            "command": f"{BASE_CMD} {OPTS[name]} --captain-ruling {ROW}"}


def _refused(bt, monkeypatch, tmp_path, argv: list[str], match: str) -> str:
    """The command is refused before anything but the ledger is read (no T1 file, Feng file, board, union or
    relabelling record; no step is set up), and its message shows no number but the allowed ones and what the command
    line said."""
    for obj, name in ((bt, "load_stage"), (bt, "load_boards"), (bt.hygiene, "load_union"),
                      (bt.feng_dirs, "feng_directions"), (fr, "committed_at_head"), (bt, "_t1_hm_step"),
                      (bt, "_t1_ms_step")):
        monkeypatch.setattr(obj, name, _no_read)
    name = "T1_dry.h5ad" if "--board" in argv else "T1.h5ad"
    with pytest.raises(SystemExit, match=match) as e:
        bt.main([*argv, "--out", str(tmp_path / name)])
    msg = str(e.value)
    given = set(re.findall(r"\d+(?:\.\d+)?", " ".join(argv)))
    assert set(re.findall(r"\d+(?:\.\d+)?", msg)) <= NUMBERS_ALLOWED | given, msg
    return msg


def _calls(bt, monkeypatch) -> dict[str, list]:
    """Around the builder's MS step, Blood's step and g_mom: per build, the arguments each gets."""
    seen = {"ms": [], "hm": [], "g_mom": []}
    ms, hm, gm = bt._t1_ms_step, bt._t1_hm_step, bt._t1_g_mom

    def ms_step(ref, library, directions, fraction, groups, held, seed, **k):
        seen["ms"].append({"fraction": fraction, "seed": seed, "kw": dict(k)})
        return ms(ref, library, directions, fraction, groups, held, seed, **k)

    def hm_step(ref, library, damping, held, seed, **k):
        seen["hm"].append({"damping": damping, "seed": seed, "kw": dict(k)})
        return hm(ref, library, damping, held, seed, **k)

    def g_mom(*a, **k):
        seen["g_mom"].append(dict(k))
        return gm(*a, **k)

    monkeypatch.setattr(bt, "_t1_ms_step", ms_step)
    monkeypatch.setattr(bt, "_t1_hm_step", hm_step)
    monkeypatch.setattr(bt, "_t1_g_mom", g_mom)
    return seen


def _factors_one(bt, monkeypatch, blood_only: bool = False) -> list:
    """momentum_shift with its damping forced to 0 (every factor exactly 1): every call, or Blood's alone."""
    calls = []
    inner = bt.t1_scope.momentum_shift

    def one(start, ref, target_like, damping, *a, **k):
        blood = k.get("within") == {"blood": ["Blood"]}
        calls.append(("blood" if blood else "other", damping))
        return inner(start, ref, target_like, 0.0 if (blood or not blood_only) else damping, *a, **k)

    monkeypatch.setattr(bt.t1_scope, "momentum_shift", one)
    return calls


def _check():
    import t1_r51_check

    return t1_r51_check


# ------------------------------------------------------------------ the constants, the help and the ruling check


def test_the_option_is_off_by_default_and_the_seven_sets_are_row_r51s(bt):
    a = _parse(bt, ["t1", *RC])
    assert a.direction_momentum_damping is None and a.held_momentum is None and a.captain_ruling is None
    assert bt.T1_R51_ROW == ROW and bt.T1_R51_BASE_FILE == "T1_s2rcMSB" and bt.T1_R51_M10_FILE == "T1_s2rcMSBm10"
    assert {f: v[0] for f, v in bt.T1_R51_FILES.items()} == OPTS
    for f, (opt, mult, damping, held, dmd) in bt.T1_R51_FILES.items():  # the values each option string gives
        a = _parse(bt, ["t1", *RC, *opt.split()])
        assert (a.feng_fraction_mult, a.g_mom_damping, a.held_momentum, a.direction_momentum_damping) == (
            mult, damping, held, dmd), f
        assert a.direction_momentum == "all"
        assert bt.T1_R51_COMMAND.format(option=opt) == " ".join(["t1", *RC, *ARGS[f]])
    assert bt.T1_R51_DMD_CHOICES == ("1.0",)
    assert set(bt.T1_R51_CHANGES) == set(OPTS)
    # the base file's values are MSB's, and the earlier options' constants are unchanged
    assert bt.T1_MSB_MULT == 1.5 and bt.G_MOM_DAMPINGS == (0.5, 1.0) and bt.harness.B0_DAMPING == 0.5
    assert bt.T1_MSB_COMMAND == " ".join(["t1", *RC, *MSB]) and bt.T1_MS_COMMAND == " ".join(["t1", *RC, *MS_ARGS])
    assert bt.T1_HM_CHOICES == ("0.5", "1.0", "2.0") and bt.T1_BANKED_FILE == "T1_s2rcB25"
    assert bt.T1_HM_STREAM == 6 and bt.T1_MS_STREAM == 1 and bt.T1_MS_FIT_LIBRARIES == ("1", "2")


def test_the_help_and_the_docstring_name_the_new_files(bt):
    import io
    from contextlib import redirect_stdout

    buf = io.StringIO()
    with redirect_stdout(buf), pytest.raises(SystemExit):
        bt.build_parser().parse_args(["t1", "--help"])
    text = re.sub(r"-\s+", "-", " ".join(buf.getvalue().split()))
    for word in ("--direction-momentum-damping {1.0}", "row R51's file T1_s2rcMSBm10", *OPTS, *OPTS.values(),
                 "With --captain-ruling R51: row R51's seven files"):
        assert word in text, word
    doc = " ".join(bt.__doc__.split())
    for word in ("Row R51 of experiments/harness/TRIES.md", *OPTS, "--direction-momentum-damping 1.0",
                 "as a whole backticked string", "library 3 never in the fit", "[seed, 6]", "held_momentum_step",
                 "taken or refused exactly as before"):
        assert word in doc, word


@pytest.mark.parametrize("text, ok", [
    (ROW_R51, True),
    (f"| R51 (the captain's yes of X; ruling, no run; `{OPTS['T1_s2rcMS20']}` `--captain-ruling R51`) |", True),
    (f"| R51 (the captain's order of X, ruling, no run: `{OPTS['T1_s2rcMS20']}`, each only with "
     "`--captain-ruling R51`; not a try) | T1 |", True),
    # the option inside a longer backticked string, or around a shorter one, is not the option
    (f"| R51 (the captain's order of X; ruling, no run; `{OPTS['T1_s2rcMS20']} --g-mom-damping 1.0` "
     "`--captain-ruling R51`) |", False),
    (f"| R51 (the captain's order of X; ruling, no run; `t1 --g-mom {OPTS['T1_s2rcMS20']}` `--captain-ruling R51`) |",
     False),
    ("| R51 (the captain's order of X; ruling, no run; `--feng-fraction-mult 2.0` `--direction-momentum all` "
     "`--captain-ruling R51`) |", False),
    (f"| R51 (the captain's order of X; ruling, no run; {OPTS['T1_s2rcMS20']} `--captain-ruling R51`) |", False),
    (f"| R51 (the captain's order of X; ruling, no run; `{OPTS['T1_s2rcMS20']}` --captain-ruling R51) |", False),
    (f"| R51 (the captain's order of X; ruling, no run; `{OPTS['T1_s2rcMS20']}` `--captain-ruling R510`) |", False),
    (f"| R51 (the captain's order of X; ruling, no run; `{OPTS['T1_s2rcMS20']}` `--captain-ruling R22`) |", False),
    (f"| R51 (the captain's order of X; ruling, no run; `{OPTS['T1_s2rcMS20']} --captain-ruling R51`) |", False),
    (f"| R51 (the captain's order of X; a try; `{OPTS['T1_s2rcMS20']}` `--captain-ruling R51`) |", False),
    (f"| R51 (a proposal; ruling, no run; `{OPTS['T1_s2rcMS20']}` `--captain-ruling R51`) |", False),
    (f"| R51 awaits (the captain's order; ruling, no run; `{OPTS['T1_s2rcMS20']}` `--captain-ruling R51`) |", False),
    (f"| R51 (the captain's order of X; ruling, no run) | T1 | `{OPTS['T1_s2rcMS20']}` `--captain-ruling R51` |",
     False),
    (f"| R510 (the captain's order of X; ruling, no run; `{OPTS['T1_s2rcMS20']}` `--captain-ruling R51`) |", False),
    (f"| R5 (the captain's order of X; ruling, no run; `{OPTS['T1_s2rcMS20']}` `--captain-ruling R51`) |", False),
])
def test_the_ruling_check_takes_only_a_single_row_r51_naming_the_whole_option_string(bt, tmp_path, text, ok):
    p = tmp_path / "TRIES.md"
    p.write_text(f"# ledger\n{text}\n")
    assert bt._r51_row_names(OPTS["T1_s2rcMS20"], p) is ok
    p.write_text(f"# ledger\n{text}\n{text}\n")  # two rows with one ID: never
    assert not bt._r51_row_names(OPTS["T1_s2rcMS20"], p)


def test_each_option_string_is_named_alone_and_no_string_passes_on_another(bt, ledger, tmp_path):
    for opt in OPTS.values():
        assert bt._r51_row_names(opt, ledger)
    for opt in ("--feng-fraction-mult 1.5 --direction-momentum all", "--direction-momentum all",
                "--feng-fraction-mult 2.5 --direction-momentum all", "--held-momentum 0.5", "--g-mom-damping 1.0",
                "--feng-fraction-mult 2.0 --direction-momentum all --g-mom-damping 1.0"):  # MSB's pair, parts, P1
        assert not bt._r51_row_names(opt, ledger), opt
    p = tmp_path / "one.md"
    for name, opt in OPTS.items():  # a row naming one set orders that set alone
        p.write_text(f"# ledger\n| R51 (the captain's order of X; ruling, no run; `{opt}` only with "
                     "`--captain-ruling R51`) |\n")
        assert [n for n, o in OPTS.items() if bt._r51_row_names(o, p)] == [name]
    # a longer string that holds d10's string names neither d10 nor MSB's pair
    p.write_text(f"# ledger\n| R51 (the captain's order of X; ruling, no run; `{OPTS['T1_s2rcMSBd10']} "
                 "--held-momentum 0.5` `--captain-ruling R51`) |\n")
    assert not any(bt._r51_row_names(o, p) for o in OPTS.values())


# ------------------------------------------------------------------ routing: each file under its row


@pytest.mark.parametrize("name", list(OPTS))
def test_each_file_is_taken_with_its_row_and_routes_nowhere_else(bt, ledger, monkeypatch, name):
    monkeypatch.setattr(fr, "committed_at_head", lambda p: json.dumps({"route": "recomputed"}).encode())
    a = _parse(bt, ["t1", *RC, *ARGS[name]])
    assert bt._r51_file(a) == name and bt._r51_ordered(a)
    assert bt._t1_r51(a) == _want(bt, name)
    assert bt._bet_refusals(a) is None and bt._feng_recomputed_refusals(a) is None and bt._t1_pooled(a) is None
    for other in (bt._ruled_fraction_ordered, bt._ruled_damping_ordered, bt._ruled_line3_ordered,
                  bt._ruled_line5_ordered, bt._ruled_line6_ordered, bt._ruled_line7_ordered, bt._ms_ordered,
                  bt._msb_ordered, bt._pooled_ordered, bt._hm_ordered, bt._ms25_ordered):
        assert not other(a), other.__name__  # no earlier row orders the command
    # without the row (or with a row naming another set only) the command is not row R51's
    ledger.write_text(_ledger_text(r51=None))
    assert bt._r51_file(a) is None
    if name != "T1_s2rcMSBm10":  # (the new option alone is refused by _t1_r51 itself: below)
        assert bt._t1_r51(a) is None
    other = next(o for n, o in OPTS.items() if n != name)
    ledger.write_text(_ledger_text(r51=f"| R51 (the captain's order of X; ruling, no run; `{other}` only with "
                                       "`--captain-ruling R51`) |"))
    assert bt._r51_file(a) is None and not bt._r51_ordered(a)


def test_without_the_row_each_file_is_refused_as_before(bt, ledger, monkeypatch, tmp_path):
    """Without row R51 the six commands that need no new option are refused exactly as before (the earlier paths'
    messages); the seventh is refused by the new option's own guard."""
    ledger.write_text(_ledger_text(r51=None))
    before = {"T1_s2rcMS20": "refused with --captain-ruling R51: only row R23 orders it",
              "T1_s2rcMS30": "refused with --captain-ruling R51: only row R23 orders it",
              "T1_s2rcMSBd10": "refused with --captain-ruling R51: only row R23 orders it",
              "T1_s2rcMSBd15": "refused with --captain-ruling R51: only row R23 orders it",
              "T1_s2rcMSBH05": "named but not coded", "T1_s2rcMSBH10": "named but not coded",
              "T1_s2rcMSBm10": "--direction-momentum-damping 1.0 is refused: it takes 1.0 only"}
    for name, match in before.items():
        _refused(bt, monkeypatch, tmp_path, ["t1", *RC, *ARGS[name]], match)


# ------------------------------------------------------------------ nothing existing changes


@pytest.mark.parametrize("name", list(BEFORE))
def test_the_existing_commands_give_the_same_bytes_and_records(bt, rc4, ledger, tmp_path, name):
    extra, file_sha, record_sha = BEFORE[name]
    f, rec = _build(bt, tmp_path, f"T1_rc{name}", extra)
    assert _sha(f) == file_sha and _record_sha(rec) == record_sha


def test_the_ruled_files_of_r20_and_r22s_seven_lines_and_r37_give_the_bytes_and_records_as_before(bt, rc4, ledger,
                                                                                                  tmp_path):
    files = [("T1_s2rcB2", R20, R20_BEFORE), ("T1_s2rcB25", LINE1[2.5], LINE1_BEFORE[2.5]),
             ("T1_s2rcB30", LINE1[3.0], LINE1_BEFORE[3.0]), ("T1_s2rcA15", LINE2[1.5], LINE2_BEFORE[1.5]),
             ("T1_s2rcA20", LINE2[2.0], LINE2_BEFORE[2.0]), ("T1_s2rcCm20", LINE3[2.0], LINE3_BEFORE[2.0]),
             ("T1_s2rcCm25", LINE3[2.5], LINE3_BEFORE[2.5]), ("T1_s2rcB2d15", LINE5[1.5], LINE5_BEFORE[1.5]),
             ("T1_s2rcB2d20", LINE5[2.0], LINE5_BEFORE[2.0]), ("T1_s2rcB25d15", LINE6[1.5], LINE6_BEFORE[1.5]),
             ("T1_s2rcB25d20", LINE6[2.0], LINE6_BEFORE[2.0]), ("T1_s2rcCm30", LINE7, LINE7_BEFORE),
             ("T1_s2rc_pooled", ["--base", "pooled", "--captain-ruling", "R37"], POOLED_BEFORE)]
    for name, argv, before in files:
        f, rec = _build(bt, tmp_path, name, argv)
        assert (_sha(f), _record_sha(rec)) == before, name
        assert "R51" not in json.dumps(rec), name


def test_ms_msb_the_banked_file_the_default_and_rows_r44s_and_r45s_files_give_the_same_bytes(bt, ms4, ledger,
                                                                                            tmp_path):
    files = [("T1_s2rc_MS", MS_ARGS, MS_BEFORE), ("T1_s2rcMSB", MSB, MSB_BEFORE),
             ("T1_s2rcB25", BANKED, BANKED_MS4_BEFORE), ("T1_s2rc", [], DEFAULT_MS4_BEFORE),
             *((HM_NAMES[v], HM[v], HM_BEFORE[v]) for v in HM_BEFORE), ("T1_s2rcB25MS", MS25, MS25_BEFORE)]
    for name, argv, before in files:
        f, rec = _build(bt, tmp_path, name, argv)
        assert (_sha(f), _record_sha(rec)) == before, name
        assert "R51" not in json.dumps(rec) and "damping_option" not in json.dumps(rec), name


def test_without_row_r51s_commands_its_steps_get_what_they_got_before(bt, ms4, ledger, monkeypatch, tmp_path):
    seen = _calls(bt, monkeypatch)
    for name, argv in (("T1_s2rcMSB", MSB), ("T1_s2rcB25MS", MS25), ("T1_s2rcB25H05", HM[0.5]), ("T1_s2rcB25", BANKED)):
        _build(bt, tmp_path, name, argv)
    assert [c["kw"] for c in seen["ms"]] == [{}, {}]  # MS's step as before: no damping given
    assert [c["kw"] for c in seen["hm"]] == [{}]  # Blood's step as before: no base file
    assert all("line" not in k for k in seen["g_mom"]) and len(seen["g_mom"]) == 4


# ------------------------------------------------------------------ the seven files


@pytest.mark.parametrize("name", list(OPTS))
def test_each_file_builds_the_same_bytes_twice_and_differs_from_msb(bt, ms4, ledger, tmp_path, name):
    f0, _ = _build(bt, tmp_path, "T1_s2rcMSB", MSB)
    fa, ra = _build(bt, tmp_path, f"{name}a", ARGS[name])
    fb, rb = _build(bt, tmp_path, f"{name}b", ARGS[name])
    assert fa.read_bytes() == fb.read_bytes() and _record_sha(ra) == _record_sha(rb)
    x0, x1 = _same_frame(f0, fa)  # the same cells, order, obs and var as MSB's file
    assert (x0 != x1).any()
    known = {MSB_BEFORE[0], MS_BEFORE[0], BANKED_MS4_BEFORE[0], DEFAULT_MS4_BEFORE[0], MS25_BEFORE[0],
             *(b[0] for b in HM_BEFORE.values())}
    assert _sha(fa) not in known


def test_the_seven_files_are_seven_files(bt, ms4, ledger, tmp_path):
    files = {name: _build(bt, tmp_path, name, ARGS[name])[0].read_bytes() for name in OPTS}
    assert len(set(files.values())) == 7


@pytest.mark.parametrize("name", list(OPTS))
def test_each_step_gets_the_files_values_and_nothing_else(bt, ms4, ledger, monkeypatch, tmp_path, name):
    seen = _calls(bt, monkeypatch)
    _build(bt, tmp_path, "T1_s2rcMSB", MSB)
    _build(bt, tmp_path, name, ARGS[name])
    mult, damping, held, ms_damping = VALUES[name]
    assert seen["ms"][0] == {"fraction": 1.5 * 0.4, "seed": SEED, "kw": {}}  # MSB's own call
    assert seen["ms"][1] == {"fraction": mult * 0.4, "seed": SEED, "kw": {"damping": ms_damping}}
    assert seen["g_mom"][0]["damping"] == 0.5 and seen["g_mom"][0]["fraction_mult"] == 1.5
    assert seen["g_mom"][1]["damping"] == damping and seen["g_mom"][1]["fraction_mult"] == mult
    assert seen["g_mom"][1]["line"] == f"row R51's file {name}, T1's line 10, on T1_s2rcMSB's command"
    assert {k: v for k, v in seen["g_mom"][1].items() if k not in ("damping", "fraction_mult", "line")} == {
        k: v for k, v in seen["g_mom"][0].items() if k not in ("damping", "fraction_mult")}
    if held is None:
        assert seen["hm"] == []
    else:
        assert seen["hm"] == [{"damping": held, "seed": SEED, "kw": {"base_file": "T1_s2rcMSB"}}]
    fractions = [e["fraction"] for e in ms4]  # MS takes the atlas step in its own step: build_submission none
    assert fractions == [None, None] and all(e["directions"] == [] for e in ms4)


@pytest.mark.parametrize("name", H_FILES)
def test_each_h_file_is_msbs_but_for_the_blood_rows(bt, ms4, ledger, tmp_path, name):
    f0, _ = _build(bt, tmp_path, "T1_s2rcMSB", MSB)
    fh, _ = _build(bt, tmp_path, name, ARGS[name])
    x0, xh = _same_frame(f0, fh)
    pop = ms4[-1]["before"]  # the population the H build's post_cells gets: MSB's draw at the same seed
    labels = _labels(pop)
    blood = labels == "Blood"
    assert blood.any() and not blood.all() and {"Macrophage", "Early VIC"} <= set(labels)
    assert x0[~blood].tobytes() == xh[~blood].tobytes()  # every non-Blood row: MSB's bytes
    changed = (x0 != xh).any(axis=1)
    assert changed.any() and not (changed & ~blood).any()  # only Blood rows move, and some do
    b0 = ms4[0]["before"]  # the draw is MSB's
    assert list(b0.obs_names) == list(pop.obs_names) and b0.obs.equals(pop.obs)
    assert np.array_equal(_dense(b0.layers["counts"]), _dense(pop.layers["counts"]))


def test_the_two_h_files_move_blood_further_at_the_larger_damping(bt, ms4, ledger, tmp_path):
    x0 = _dense(ad.read_h5ad(_build(bt, tmp_path, "T1_s2rcMSB", MSB)[0]).X)
    blood = _labels(ms4[0]["before"]) == "Blood"
    moved = [float(np.abs(_dense(ad.read_h5ad(_build(bt, tmp_path, n, ARGS[n])[0]).X)[blood] - x0[blood]).mean())
             for n in H_FILES]
    assert 0 < moved[0] < moved[1]


def test_the_h_step_is_row_r44s_on_msbs_population_after_ms(bt, ms4, t1_env, ledger, tmp_path,  # noqa: F811
                                                            monkeypatch):
    calls = []
    inner = bt.t1_scope.momentum_shift

    def spy(start, ref, target_like, damping, held_flat=(), **k):
        out, rep = inner(start, ref, target_like, damping, held_flat, **k)
        calls.append({"ref": list(ref.obs_names), "labels": sorted(set(_labels(target_like))), "damping": damping,
                      "kw": dict(k), "rows": rep.rows})
        return out, rep

    monkeypatch.setattr(bt.t1_scope, "momentum_shift", spy)
    _build(bt, tmp_path, "T1_s2rcMSBH10", ARGS["T1_s2rcMSBH10"])
    assert len(calls) == 2  # MS's momentum first, then Blood's
    ms, blood = calls
    assert ms["damping"] == 0.5 and "within" not in ms["kw"] and ms["kw"]["seed"] == [SEED, 1]
    assert "Blood" not in ms["labels"]
    lib = np.array([n.rsplit("_", 1)[1] for n in t1_env.ref.obs_names])
    fit = [n for n, x in zip(t1_env.ref.obs_names, lib, strict=True) if x in ("1", "2")]  # library 3 never in a fit
    assert ms["ref"] == fit == blood["ref"]
    spec = bt.t1_scope.load_spec()
    assert blood["damping"] == 1.0 and blood["labels"] == ["Blood"]
    assert blood["kw"] == {"label_key": "celltype", "eps": 0.5, "min_cells": 20, "keep_depth": True,
                           "within": {"blood": ["Blood"]}, "seed": [SEED, 6],
                           "start_fallback": spec.fallbacks("E9.5_RNA"), "ref_fallback": spec.fallbacks("E10.5_RNA")}
    assert [r["type"] for r in blood["rows"]] == ["blood"] and blood["rows"][0]["status"] == "shifted"


# ------------------------------------------------------------------ the checks at MSB's values


@pytest.mark.parametrize("name", list(OPTS))
def test_the_checker_builds_msbs_bytes_with_each_files_setting_at_msbs_value(bt, ms4, ledger, tmp_path, name):
    chk = _check()
    f0, _ = _build(bt, tmp_path, "T1_s2rcMSB", MSB)
    out = tmp_path / f"{name}_msb.h5ad"
    before = (bt._t1_r51, bt._t1_ms_step, bt._t1_g_mom, bt.t1_scope.momentum_shift)
    assert chk.main(["msb", "t1", *RC, *ARGS[name], "--out", str(out)]) == 0
    assert out.read_bytes() == f0.read_bytes()
    assert (bt._t1_r51, bt._t1_ms_step, bt._t1_g_mom, bt.t1_scope.momentum_shift) == before  # put back
    rec = json.loads(out.with_suffix(".build.json").read_text())
    assert rec["t1_candidate"]["file"] == name  # the file's own path, row R51's route


@pytest.mark.parametrize("name", H_FILES)
def test_with_bloods_factors_at_1_the_h_file_is_msbs_bytes(bt, ms4, ledger, tmp_path, monkeypatch, name):
    f0, _ = _build(bt, tmp_path, "T1_s2rcMSB", MSB)
    calls = _factors_one(bt, monkeypatch, blood_only=True)
    f1, r1 = _build(bt, tmp_path, f"{name}_one", ARGS[name])
    assert calls == [("other", 0.5), ("blood", VALUES[name][2])]  # MS's at B0's damping, then Blood's, forced to 0
    assert f1.read_bytes() == f0.read_bytes()
    assert r1["t1_candidate"]["held_momentum_step"]["status"] == "shifted"


def test_m10_with_ms_damping_at_the_default_is_msbs_bytes(bt, ms4, ledger, tmp_path, monkeypatch):
    f0, _ = _build(bt, tmp_path, "T1_s2rcMSB", MSB)
    inner = bt._t1_r51
    monkeypatch.setattr(bt, "_t1_r51", lambda a: {**inner(a), "direction_momentum_damping": 0.5})
    f1, r1 = _build(bt, tmp_path, "T1_s2rcMSBm10_default", ARGS["T1_s2rcMSBm10"])
    assert f1.read_bytes() == f0.read_bytes() and r1["t1_candidate"]["momentum"]["damping"] == 0.5
    monkeypatch.setattr(bt, "_t1_r51", inner)
    f2, _ = _build(bt, tmp_path, "T1_s2rcMSBm10", ARGS["T1_s2rcMSBm10"])
    assert f2.read_bytes() != f0.read_bytes()


@pytest.mark.parametrize("name, plain", [("T1_s2rcMS20", R20), ("T1_s2rcMS30", LINE1[3.0]),
                                         ("T1_s2rcMSBd10", ["--g-mom-damping", "1.0", "--feng-fraction-mult", "1.5"]),
                                         ("T1_s2rcMSBm10", ["--feng-fraction-mult", "1.5"]),
                                         ("T1_s2rcMSBH05", ["--feng-fraction-mult", "1.5"])])
def test_with_every_momentum_factor_at_1_a_file_is_its_commands_file_without_ms(bt, ms4, ledger, tmp_path, monkeypatch,
                                                                              name, plain):
    """MS's momentum (and Blood's) at factor 1 leaves the atlas step at the file's multiplier and g_mom at its damping:
    MS20 is row R20's 2.0 file, MS30 row R22's 3.0 file, MSBd10 candidate C, MSBm10 and MSBH05 candidate B."""
    f0, _ = _build(bt, tmp_path, "plain", plain)
    calls = _factors_one(bt, monkeypatch)
    f1, _ = _build(bt, tmp_path, f"{name}_all_one", ARGS[name])
    assert [c[0] for c in calls] == (["other", "blood"] if name in H_FILES else ["other"])
    assert f1.read_bytes() == f0.read_bytes()


# ------------------------------------------------------------------ the record


@pytest.mark.parametrize("name", list(OPTS))
def test_the_record_is_msbs_plus_its_candidate_field(bt, ms4, ledger, tmp_path, name):
    _, r0 = _build(bt, tmp_path, "T1_s2rcMSB", MSB)
    _, r1 = _build(bt, tmp_path, name, ARGS[name])
    mult, damping, held, ms_damping = VALUES[name]
    skip = {"t1_candidate", "fraction", "g_mom"}
    assert {k: v for k, v in r1.items() if k not in skip} == {k: v for k, v in r0.items()
                                                                if k not in skip}
    assert r1["fraction"] == mult * 0.4
    g0, g1 = dict(r0["g_mom"]), dict(r1["g_mom"])
    assert g1.pop("damping") == damping and g0.pop("damping") == 0.5
    rule1, rule0 = g1.pop("rule"), g0.pop("rule")
    g0.pop("median_depth_scale"), g1.pop("median_depth_scale")  # the edit's own outcome on the file's population
    assert g1 == g0
    if damping == 0.5:
        assert rule1 == rule0
    else:
        line = f"row R51's file {name}, T1's line 10, on T1_s2rcMSB's command"
        assert rule1 == rule0.replace("damped by B0_DAMPING", f"damped by {damping:g} in place of B0_DAMPING ({line})")
    cand = r1["t1_candidate"]
    assert {k: v for k, v in cand.items() if k not in ("momentum", "calendar_fraction", "applied_fraction",
                                                       "held_momentum_step")} == _want(bt, name)
    assert cand["calendar_fraction"] == 0.4 and cand["applied_fraction"] == mult * 0.4
    m0, m1 = dict(r0["t1_candidate"]["momentum"]), dict(cand["momentum"])
    assert m1.pop("damping") == ms_damping and m0.pop("damping") == 0.5
    option = m1.pop("damping_option", None)
    assert m1 == m0  # MS's record otherwise MSB's
    assert (option is not None) == (name == "T1_s2rcMSBm10")
    if option is not None:
        assert "--direction-momentum-damping 1.0" in option and "row R51's file T1_s2rcMSBm10" in option
    assert ("held_momentum_step" in cand) == (held is not None)
    if held is not None:
        h = cand["held_momentum_step"]
        assert h["fit"] == {"start": "E9.5_RNA", "later": "E10.5_RNA", "later_libraries": ["1", "2"],
                            "library_3_in_fit": False, "within": {"blood": ["Blood"]},
                            "start_cells": "own, restricted", "later_cells": "own, restricted"}
        assert (h["label"], h["status"], h["damping"], h["stream"]) == ("Blood", "shifted", held, [SEED, 6])
        assert "after MS's step" in h["rule"] and "as T1_s2rcMSB has it" in h["rule"]
        assert "banked file" not in h["rule"] and "after the atlas step" not in h["rule"]
    assert not re.search(r"n_target|n_start|n_ref|median_depth|\bn_cells\b", json.dumps(cand))  # no cell number


# ------------------------------------------------------------------ refusals, before anything but the ledger is read


@pytest.mark.parametrize("name", list(OPTS))
@pytest.mark.parametrize("extra, match", [
    (["--released-bet", "share_trend", "--bet-strength", "0.5"], "refused with --released-bet, --bet-strength: it is "
                                                                 "T1_s2rcMSB's command with one change"),
    (["--bet-strength", "0.5"], "refused with --bet-strength"),
    (["--no-hold-unbracketed"], "refused with --no-hold-unbracketed: the cells without a group keep R1's hold"),
    (["--seed", "1000002"], "changes T1_s2rcMSB only: it needs the command"),
    (["--board", "T1:val"], "changes T1_s2rcMSB only"),
    (["--cm-share", "0.6"], "--feng-groups recomputed builds one file"),
    (["--herd"], "--feng-groups recomputed builds one file"),
    (["--no-feng-hygiene"], "--feng-groups recomputed builds one file"),
    (["--strong-null"], "--feng-groups recomputed builds one file"),
    (["--base", "pooled"], "--feng-groups recomputed builds one file"),
    (["--feng-none-ruling", "x"], "--feng-groups recomputed builds one file"),
])
def test_every_other_use_of_a_file_under_its_row_is_refused(bt, ledger, monkeypatch, tmp_path, name, extra, match):
    _refused(bt, monkeypatch, tmp_path, ["t1", *RC, *ARGS[name], *extra], match)


@pytest.mark.parametrize("name", list(OPTS))
@pytest.mark.parametrize("cmd", [["--g-mom", "--base", "lib3"], [], ["--g-mom", "--base", "lib3", "--feng-groups",
                                                                       "none"]])
def test_a_file_is_refused_on_any_other_command(bt, ledger, monkeypatch, tmp_path, name, cmd):
    _refused(bt, monkeypatch, tmp_path, ["t1", *cmd, *ARGS[name]], f"row R51's file {name} changes T1_s2rcMSB only")


@pytest.mark.parametrize("argv, match", [
    # two of the changes together, the struck pair P1, and a value no file takes: the earlier paths' refusals
    ([*OPTS["T1_s2rcMS20"].split(), "--g-mom-damping", "1.0"], "refused with --captain-ruling R51: only row R23"),
    ([*OPTS["T1_s2rcMSBd10"].split(), "--held-momentum", "0.5"], "named but not coded"),
    ([*OPTS["T1_s2rcMSBH05"].split(), "--g-mom-damping", "1.5"], "named but not coded"),
    ([*OPTS["T1_s2rcMS30"].split(), "--held-momentum", "1.0"], "named but not coded"),
    (["--feng-fraction-mult", "1.5", "--direction-momentum", "all", "--held-momentum", "2.0"], "named but not coded"),
    (["--feng-fraction-mult", "1.5", "--direction-momentum", "all", "--g-mom-damping", "2.0"],
     "refused with --captain-ruling R51: only row R23"),
    (["--feng-fraction-mult", "2.5", "--direction-momentum", "all"], "refused with --captain-ruling R51: only row R23"),
    (["--feng-fraction-mult", "1.5", "--direction-momentum", "all"], "refused with --captain-ruling R51: only row R23"),
    (["--feng-fraction-mult", "2.0"], "refused with --captain-ruling R51: only row R20 orders it"),
    (["--g-mom-damping", "1.0", "--feng-fraction-mult", "1.5"], "go with --released-bet only"),
    (["--feng-fraction-mult", "2.5", "--held-momentum", "0.5"], "no single row R51 whose first cell"),
    ([*OPTS["T1_s2rcMSBH05"].split(), "--direction-momentum-damping", "1.0"], "--direction-momentum-damping 1.0 is "
                                                                              "refused: it takes 1.0 only"),
    ([*OPTS["T1_s2rcMSBd10"].split(), "--direction-momentum-damping", "1.0"], "it takes 1.0 only"),
    ([*OPTS["T1_s2rcMS20"].split(), "--direction-momentum-damping", "1.0"], "it takes 1.0 only"),
])
def test_combinations_the_row_does_not_name_are_refused(bt, ledger, monkeypatch, tmp_path, argv, match):
    _refused(bt, monkeypatch, tmp_path, ["t1", *RC, *argv, "--captain-ruling", ROW], match)


def test_the_struck_pair_stays_refused_even_where_a_row_names_it(bt, ledger, monkeypatch, tmp_path):
    p1 = "--feng-fraction-mult 2.0 --direction-momentum all --g-mom-damping 1.0"
    ledger.write_text(_ledger_text(r51=ROW_R51.replace("; not a try)", f", `{p1}`; not a try)")))
    assert bt._r51_row_names(p1)
    _refused(bt, monkeypatch, tmp_path, ["t1", *RC, *p1.split(), "--captain-ruling", ROW],
             "refused with --captain-ruling R51: only row R23 orders it")


@pytest.mark.parametrize("argv", [
    [*OPTS["T1_s2rcMSBm10"].split()[:-1], "0.5"], [*OPTS["T1_s2rcMSBm10"].split()[:-1], "1"],
    [*OPTS["T1_s2rcMSBm10"].split()[:-1], "1.00"], [*OPTS["T1_s2rcMSBm10"].split()[:-1], "2.0"],
    [*OPTS["T1_s2rcMSBm10"].split()[:-1], "abc"],
    ["--direction-momentum-damping", "1.0"], ["--feng-fraction-mult", "1.5", "--direction-momentum-damping", "1.0"],
    ["--direction-momentum", "all", "--direction-momentum-damping", "1.0"],
])
@pytest.mark.parametrize("ruling", [None, "R51", "R26", "R23", "R44"])
def test_the_new_option_is_refused_on_every_other_command(bt, ledger, monkeypatch, tmp_path, argv, ruling):
    tail = [] if ruling is None else ["--captain-ruling", ruling]
    _refused(bt, monkeypatch, tmp_path, ["t1", *RC, *argv, *tail], "--direction-momentum-damping .* is refused: it "
                                                                   "takes 1.0 only, and only in row R51's file "
                                                                   "T1_s2rcMSBm10")


@pytest.mark.parametrize("ruling", [None, "R26", "R23", "R22", "R45", "R510", "r51"])
def test_m10s_set_with_another_ruling_is_refused(bt, ledger, monkeypatch, tmp_path, ruling):
    tail = [] if ruling is None else ["--captain-ruling", ruling]
    _refused(bt, monkeypatch, tmp_path, ["t1", *RC, *OPTS["T1_s2rcMSBm10"].split(), *tail], "it takes 1.0 only")


@pytest.mark.parametrize("name", list(OPTS))
def test_a_file_whose_set_the_row_does_not_name_is_refused(bt, ledger, monkeypatch, tmp_path, name):
    rows = [f"| R51 (the captain's order of X; ruling, no run; `{OPTS[name]}`) |",  # no ruling flag
            f"| R51 (the captain's order of X; ruling, no run; `{OPTS[name]}` `--captain-ruling R22`) |",
            f"| R51 (the captain's order of X; a try; `{OPTS[name]}` `--captain-ruling R51`) |",
            f"| R51 (a proposal; ruling, no run; `{OPTS[name]}` `--captain-ruling R51`) |",
            f"| R51 (the captain's order of X; ruling, no run) | T1 | `{OPTS[name]}` `--captain-ruling R51` |",
            f"| R51 (the captain's order of X; ruling, no run; `{OPTS[name]} x` `--captain-ruling R51`) |",
            f"{ROW_R51}\n{ROW_R51}"]  # two rows with one ID
    for row in rows:
        ledger.write_text(_ledger_text(r51=row))
        _refused(bt, monkeypatch, tmp_path, ["t1", *RC, *ARGS[name]], "")
        assert bt._r51_file(_parse(bt, ["t1", *RC, *ARGS[name]])) is None


def test_the_bet_refusal_lets_only_the_seven_commands_through(bt, ledger):
    def ns(name, **k):
        mult, damping, held, dmd = bt.T1_R51_FILES[name][1:]
        base = dict(released_bet=None, bet_strength=None, captain_ruling=ROW, feng_fraction_mult=mult,
                    g_mom_damping=damping, direction_momentum="all", held_momentum=held,
                    direction_momentum_damping=dmd, feng_groups="recomputed", seed=SEED, board=None, base="lib3")
        return SimpleNamespace(**{**base, **k})

    for name in OPTS:
        assert bt._bet_refusals(ns(name)) is None and bt._r51_ordered(ns(name))
        for bad in (ns(name, captain_ruling="R26"), ns(name, bet_strength=0.5), ns(name, direction_momentum=None),
                    ns(name, feng_fraction_mult=2.5)):
            assert not bt._r51_ordered(bad)
    ledger.write_text(_ledger_text(r51=None))
    for name in OPTS:
        assert not bt._r51_ordered(ns(name))
        with pytest.raises(SystemExit, match="with --released-bet only"):
            bt._bet_refusals(ns(name))


# ------------------------------------------------------------------ the old builder against this one


OLD_COMMIT = "32c8c09"  # main before these options existed
OLD_SHA = "4aa1fc236d600b4666428a1ae98d64b22883e8bc2f70f19cd39383f04ef9c337"  # its experiments/builds/build_test.py


class _Reached(Exception):
    """The command passed every refusal and reached its first read (the argument names which)."""


def _reach(where: str):
    def read(*a, **k):
        raise _Reached(where)

    return read


def _old_builder(tmp_path: Path):
    """The builder at OLD_COMMIT, read with `git show` (skipped where that commit is not in the checkout)."""
    env = {**git_env(), "GIT_OPTIONAL_LOCKS": "0"}
    res = subprocess.run(["git", "-C", str(ROOT), "show", f"{OLD_COMMIT}:experiments/builds/build_test.py"],
                         capture_output=True, env=env)
    if res.returncode != 0:
        pytest.skip(f"commit {OLD_COMMIT} is not available in this checkout")
    assert hashlib.sha256(res.stdout).hexdigest() == OLD_SHA
    path = tmp_path / "old" / "build_test_before_r51.py"
    path.parent.mkdir()
    path.write_bytes(res.stdout)
    spec = importlib.util.spec_from_file_location("build_test_before_r51", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


CMDS = [RC, [], ["--g-mom", "--base", "lib3"], ["--g-mom", "--base", "pooled", "--feng-groups", "recomputed"],
        ["--base", "lib3", "--feng-groups", "recomputed"], ["--g-mom", "--base", "lib3", "--feng-groups", "none"]]
DAMPINGS = [None, "0.5", "1.0", "1.5"]
MULTS = [None, "1.5", "2.0", "2.5", "3.0"]
HELD = [None, "0.5", "1.0", "2.0"]
DMS = [None, "all"]
DMD = [None, "1.0"]
RULINGS = [None, "R51", "R20", "R22", "R23", "R26", "R37", "R44", "R45"]  # every ruling row a T1 file is built under
EXTRAS = [["--seed", "1000002"], ["--board", "T1:val"], ["--released-bet", "share_trend", "--bet-strength", "0.5"],
          ["--bet-strength", "0.5"], ["--no-hold-unbracketed"], ["--cm-share", "0.6"], ["--herd"],
          ["--no-feng-hygiene"], ["--strong-null"], ["--lib3-pair"], ["--feng-none-ruling", "x"],
          ["--blood-share", "0.1"], ["--direction-momentum-damping", "0.5"]]
GRID_SETS = 6 * 4 * 5 * 4 * 2 * 2 * 9 + 7 * 2 * 13
SETS_COMPARED = 3 * GRID_SETS


def _grid() -> list[list[str]]:
    def argv(cmd, damping, mult, held, dm, dmd, ruling, extra=()):
        out = ["t1", *cmd]
        for flag, value in (("--g-mom-damping", damping), ("--feng-fraction-mult", mult), ("--held-momentum", held),
                            ("--direction-momentum", dm), ("--direction-momentum-damping", dmd),
                            ("--captain-ruling", ruling)):
            out += [] if value is None else [flag, value]
        return [*out, *extra]

    sets = [argv(*x) for x in itertools.product(CMDS, DAMPINGS, MULTS, HELD, DMS, DMD, RULINGS)]
    for name, ruling, extra in itertools.product(OPTS, [None, ROW], EXTRAS):
        sets.append(["t1", *RC, *OPTS[name].split(), *([] if ruling is None else ["--captain-ruling", ruling]),
                     *extra])
    return sets


def _kind(bt, a, named: set[str]) -> tuple[str, str | None]:
    """("new", file) for exactly one of row R51's commands under a ledger naming its set (the recomputed command at
    the default seed on the test board, nothing else given); ("touch", file or None) for any other set holding one of
    the seven option sets with R51 while the row names it, or the new option; ("other", None) for every other set."""
    plain = (a.g_mom and a.feng_groups == "recomputed" and a.base == "lib3" and a.seed == SEED and a.board is None
             and a.released_bet is None and a.bet_strength is None and a.hold_unbracketed and a.feng_hygiene is None
             and not a.strong_null and a.cm_share is None and a.blood_share is None and not a.lib3_pair
             and not a.herd and a.feng_none_ruling is None)
    given = (a.feng_fraction_mult, a.g_mom_damping, a.held_momentum, a.direction_momentum_damping)
    name = next((n for n, v in bt.T1_R51_FILES.items() if tuple(v[1:]) == given), None)
    if name is not None and a.direction_momentum == "all" and a.captain_ruling == ROW and name in named:
        return ("new" if plain else "touch"), name
    return ("touch" if a.direction_momentum_damping is not None else "other"), None


def test_the_old_builder_and_this_one_route_every_other_argument_set_alike(bt, tmp_path, monkeypatch):
    """Side by side on the grid, under three synthetic ledgers (with row R51 naming every set; without it; with R51
    naming three sets only and a longer string holding a fourth): every set that is not one of row R51's commands under
    a row naming its set and holds no --direction-momentum-damping is refused with the same message, or reaches the
    same first read with the same candidate, MS and pooled values, in both builders, and this one's R51 route returns
    None for it; every such command is taken by this one with its values (the old builder refuses each), and every
    other set holding one of them or the new option is refused by both. Nothing is built: every read is stopped where
    it would start."""
    old = _old_builder(tmp_path)
    assert not hasattr(old, "_t1_r51") and hasattr(bt, "_t1_r51")
    ledger = tmp_path / "TRIES.md"
    record = tmp_path / "relabel" / "MANIFEST.json"  # never written: the relabelling's record is the stand-in below
    monkeypatch.setattr(fr, "RECORD", record)
    monkeypatch.setattr(fr, "committed_at_head",
                        lambda p: json.dumps({"route": "recomputed"}).encode() if Path(p) == record else None)
    monkeypatch.setattr(fr, "later_ruling", lambda *a, **k: "no later ruling (synthetic)")
    monkeypatch.setattr(bt.hygiene, "load_union", _reach("union"))
    monkeypatch.setattr(bt.strong_null, "check_approved", _reach("strong_null"))
    monkeypatch.setattr(bt.feng_dirs, "feng_directions", _reach("feng"))
    spies, parsers = {}, {}
    for mod in (old, bt):
        parsers[mod] = mod.build_parser()  # each builder's own parser, built once and reused (parsing keeps no state)
        monkeypatch.setattr(mod, "build_parser", lambda p=parsers[mod]: p)
        monkeypatch.setattr(mod, "TRIES_FILE", ledger)
        monkeypatch.setattr(mod, "load_stage", _reach("stage"))
        monkeypatch.setattr(mod, "load_boards", _reach("board"))
        for name in ("_t1_ms", "_t1_candidate", "_t1_pooled", *(("_t1_r51",) if mod is bt else ())):
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
            return ("refused", str(e.code)), spies.get("_t1_r51")
        except _Reached as e:
            return (("taken", e.args[0], spies.get("_t1_ms"), spies.get("_t1_candidate"), spies.get("_t1_pooled")),
                    spies.get("_t1_r51"))
        raise AssertionError(f"{argv} neither refused nor read")

    partial = ("| R51 (the captain's order of X; ruling, no run; " + ", ".join(
        f"`{OPTS[n]}`" for n in ("T1_s2rcMS20", "T1_s2rcMSBH10", "T1_s2rcMSBm10"))
        + f", `{OPTS['T1_s2rcMSBd10']} --g-mom-damping 1.5`, each only with `--captain-ruling R51`) |")
    ledgers = {"after": (_ledger_text(), set(OPTS)), "before": (_ledger_text(r51=None), set()),
               "partial": (_ledger_text(r51=partial), {"T1_s2rcMS20", "T1_s2rcMSBH10", "T1_s2rcMSBm10"})}
    grid = _grid()
    assert len({tuple(x) for x in grid}) == len(grid) == GRID_SETS
    count = {"other": 0, "touch": 0, "new": 0, "taken_alike": 0, "refused_alike": 0}
    new_files = []
    # the old parser's usage errors on --direction-momentum-damping go to a file, not the captured output
    with open(tmp_path / "stderr.txt", "w") as err, pytest.MonkeyPatch.context() as quiet:
        quiet.setattr(sys, "stderr", err)
        for text, named in ledgers.values():
            ledger.write_text(text)
            for argv in grid:
                a = _parse(bt, argv)
                kind, name = _kind(bt, a, named)
                (o, _), (n, r51) = route(old, argv), route(bt, argv)
                count[kind] += 1
                if kind == "other":
                    assert n == o and r51 is None, argv
                    count["taken_alike" if o[0] == "taken" else "refused_alike"] += 1
                    continue
                assert o[0] == "refused", argv  # the old builder took none of these
                if kind == "new":
                    assert n == ("taken", "union", None, None, None) and r51 == _want(bt, name), argv
                    new_files.append(name)
                else:
                    assert n[0] == "refused", argv
    assert sum(count[k] for k in ("other", "touch", "new")) == SETS_COMPARED
    assert count["taken_alike"] and count["refused_alike"] and count["touch"]
    # the new sets in the grid: all seven under "after" (the grid's values reach each once, on R13's command with
    # R51), MS20, MSBH10 and MSBm10 under "partial" (MSBd10 there only inside a longer string), none under "before"
    assert sorted(new_files) == sorted([*OPTS, "T1_s2rcMS20", "T1_s2rcMSBH10", "T1_s2rcMSBm10"])


# ------------------------------------------------------------------ experiments/builds/t1_r51_check.py


def _lines(capsys) -> dict[str, str]:
    out = capsys.readouterr().out.splitlines()
    assert all(re.fullmatch(r"[A-Za-z0-9_]+ (True|False)", line) for line in out), out  # booleans only
    return dict(line.split(" ") for line in out)


def test_the_checkers_msb_command_is_row_r26s_and_its_settings_are_the_files(bt):
    chk = _check()
    assert chk.MSB_ARGV == ["t1", *RC, *MSB]
    assert set(chk.SETTING) == set(OPTS)
    for name, key in chk.SETTING.items():
        mult, damping, held, ms_damping = VALUES[name]
        changed = {k for k, v, m in (("feng_fraction_mult", mult, 1.5), ("g_mom_damping", damping, 0.5),
                                     ("held_momentum", held, None), ("direction_momentum_damping", ms_damping, 0.5))
                   if v != m}
        assert changed == {key}, name  # the one setting the file changes from MSB's
    assert chk.MSB_VALUES == {"feng_fraction_mult": 1.5, "g_mom_damping": 0.5, "direction_momentum_damping": 0.5}


def test_the_checkers_msb_fails_off_row_r51s_route_and_passes_refusals_through(bt, ms4, ledger, tmp_path):
    chk = _check()
    assert chk.main(["msb", "t1", *RC, *MSB, "--out", str(tmp_path / "MSB_msb.h5ad")]) == 1  # R26's route, not R51's
    with pytest.raises(SystemExit, match="changes T1_s2rcMSB only"):
        chk.main(["msb", "t1", *RC, *ARGS["T1_s2rcMS20"], "--seed", "1000002", "--out", str(tmp_path / "x.h5ad")])
    for argv in ([], ["msb"], ["msb", "heart"], ["rows", "a"], ["one", "t1"], ["x"]):
        with pytest.raises(SystemExit, match="usage"):
            chk.main(argv)


@pytest.mark.parametrize("name, wrong", [("T1_s2rcMSBd10", "feng_fraction_mult"),
                                         ("T1_s2rcMSBm10", "g_mom_damping"), ("T1_s2rcMSBH05", "g_mom_damping"),
                                         ("T1_s2rcMS20", "direction_momentum_damping")])
def test_the_checkers_msb_fails_when_the_files_own_setting_is_left_as_it_is(bt, ms4, ledger, tmp_path, monkeypatch,
                                                                            name, wrong):
    """Forcing another setting than the file's own leaves the file's setting at its value, which the checker's calls
    then show: it exits 1."""
    chk = _check()
    monkeypatch.setitem(chk.SETTING, name, wrong)
    assert chk.main(["msb", "t1", *RC, *ARGS[name], "--out", str(tmp_path / f"{name}_wrong.h5ad")]) == 1


def _write_altered(src: Path, dst: Path, rows: np.ndarray) -> Path:
    a = ad.read_h5ad(src)
    x = _dense(a.X).copy()
    x[rows, 0] += 1.0
    a.X = x
    a.write_h5ad(dst, compression="gzip")
    return dst


def test_the_checkers_rows_pass_on_the_h_files_and_fail_on_altered_ones(bt, ms4, ledger, tmp_path, capsys,
                                                                       monkeypatch):
    import tempfile

    chk = _check()
    monkeypatch.setattr(tempfile, "tempdir", str(tmp_path / "tmp"))
    (tmp_path / "tmp").mkdir()
    f0, _ = _build(bt, tmp_path, "T1_s2rcMSB", MSB)
    hs = [_build(bt, tmp_path, name, ARGS[name])[0] for name in H_FILES]
    capsys.readouterr()
    assert chk.main(["rows", str(f0), *map(str, hs)]) == 0
    got = _lines(capsys)
    assert set(got.values()) == {"True"} and len(got) == 3 + 2 * 3
    blood = _labels(ms4[0]["before"]) == "Blood"
    pop = chk.msb_population()
    assert list(pop.obs_names) == list(ms4[0]["before"].obs_names)  # MSB's draw, the builder's own
    capsys.readouterr()
    bad_other = _write_altered(hs[0], tmp_path / "T1_s2rcMSBH05x.h5ad", np.flatnonzero(~blood)[:1])
    assert chk.main(["rows", str(f0), str(bad_other)]) == 1
    assert _lines(capsys)["T1_s2rcMSBH05x_other_rows_are_msbs"] == "False"
    assert chk.main(["rows", str(f0), str(f0)]) == 1  # MSB's file as an H file: no Blood row moved
    assert _lines(capsys)["T1_s2rcMSB_a_blood_row_moved"] == "False"
    bad_msb = _write_altered(f0, tmp_path / "T1_s2rcMSBx.h5ad", np.flatnonzero(blood)[:1])
    assert chk.main(["rows", str(bad_msb), str(hs[0])]) == 1  # the draw's Blood rows not where the file holds them
    assert _lines(capsys)["msb_blood_rows_at_drawn_counts"] == "False"
    assert [p for p in (tmp_path / "tmp").rglob("*") if p.is_file()] == []  # the draws wrote nothing


def test_the_checker_writes_and_deletes_nothing_and_prints_booleans_only():
    src = (ROOT / "experiments" / "builds" / "t1_r51_check.py").read_text()
    body = src.split('"""', 2)[2]  # after the module docstring
    code = "\n".join(line.split("#")[0] for line in body.splitlines())
    for word in ("unlink", "rmtree", "os.remove", "rmdir", "write_text", "write_h5ad", "write_bytes", "n_obs)",
                 "print(len", "print(x", "upload", "submit", ".build.json", "read_text"):
        assert word not in code, word
    assert code.count("print(") == 1 and "print(k, bool(v))" in code  # booleans only
