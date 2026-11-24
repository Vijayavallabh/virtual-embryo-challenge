"""`build_test.py t1` with the two steps of the T1 design round of 2026-11-24 on the banked file T1_s2rcB25 (R13's
command with --feng-fraction-mult 2.5 --captain-ruling R22), each under a ruling row named on the command line:
  --feng-fraction-mult 2.5 --held-momentum 0.5|1.0|2.0 --captain-ruling <ROW>   T1_s2rcB25H05, H10, H20 (F1, F7, F8)
  --feng-fraction-mult 2.5 --direction-momentum all --captain-ruling <ROW>      T1_s2rcB25MS (F2)
The held-momentum step, as the round's cold challenge fixes it (FIX (a)): on the drawn cells labelled Blood only,
Blood's own released change fitted on Blood-labelled cells at both ends (momentum_shift(..., within={"blood":
["Blood"]})) from E9.5 to E10.5's portion libraries 1 and 2 (library 3 never in the fit), with _t1_ms_step's constants
(eps 0.5, min_cells 20, the build's held genes, keep_depth), at damping d, on stream [seed, 6], after the atlas step and
before g_mom. MS at 2.5 is row R23's step exactly (_t1_ms_step) on the banked command.

Held here on synthetic stages and directions only (the build_t1 fixture of tests/test_hidden_counts.py, the recomputed
setup of tests/test_feng_relabel.py and the three-library names of tests/test_t1_msb_build.py), with temporary ledgers
in place of the repository's (no test here reads the repository's own ledger): (1) nothing that builds today changes:
every existing command (R13, A, B, C, R20, R22's lines 1 to 7, R23, R26, R37) gives the bytes and records pinned for
it, the pins below made by the builder at main ed6722a (before these options existed), twice in separate processes
with the same result; (2) every non-Blood row of an H file is the banked file's at the same seed byte for byte and a
Blood row moves; with Blood's factors forced to 1 the H file is the banked file byte for byte; with MS's factors
forced to 1 the MS file is the banked file byte for byte; the record shows Blood's fit (E9.5 to libraries 1 and 2,
Blood cells only, library 3 out, shifted); (3) each is taken only with a ruling row naming it, and every other use is
refused, before anything is read but the ledger (Blood and MS together stay refused); (4) the builder at main ed6722a
(read with `git show`) and this one route every other argument set of a grid exactly alike against synthetic ledgers;
(5) experiments/builds/t1_hm_check.py's three checks pass on these files and fail on altered ones, printing booleans
only. No released file, outside-atlas file, board or score is read."""

from __future__ import annotations

import hashlib
import importlib.util
import inspect
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
from test_feng_relabel import GROUPS4, _recomputed_env  # noqa: E402
from test_hidden_counts import E95, E105, t1_env  # noqa: E402, F401  (the synthetic build_t1 fixture)
from test_t1_candidates_build import RC, THREE, _build, _dense, _no_read, _same_frame, _spy  # noqa: E402
from test_t1_damping_line_build import DAMP as LINE2  # noqa: E402
from test_t1_damping_line_build import LINE1_BEFORE  # noqa: E402
from test_t1_fraction2_build import BEFORE, R20, _record_sha  # noqa: E402
from test_t1_fraction_line_build import R20_BEFORE  # noqa: E402
from test_t1_fraction_line_build import R22 as LINE1  # noqa: E402
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
ROW_H_ID, ROW_M_ID = "R44", "R45"  # the synthetic ruling rows' IDs below (not the ledger's own rows)
HM_OPT = {0.5: "`--held-momentum 0.5`", 1.0: "`--held-momentum 1.0`", 2.0: "`--held-momentum 2.0`"}
MS25_OPT = "`--feng-fraction-mult 2.5 --direction-momentum all`"
HM = {v: ["--feng-fraction-mult", "2.5", "--held-momentum", f"{v:.1f}", "--captain-ruling", ROW_H_ID] for v in HM_OPT}
MS25 = ["--feng-fraction-mult", "2.5", "--direction-momentum", "all", "--captain-ruling", ROW_M_ID]
BANKED = LINE1[2.5]  # the banked file T1_s2rcB25's options: --feng-fraction-mult 2.5 --captain-ruling R22
NAMES = {0.5: "T1_s2rcB25H05", 1.0: "T1_s2rcB25H10", 2.0: "T1_s2rcB25H20"}
# The ruling rows' first cells, in the words the ledger's first cells use, inside synthetic rows of that form
FIRST_CELL_H = ("R44 (the captain's yes of a later day, \"Yes, T1_s2rcB25H05: lift the close and reverse R1 for Blood "
                "in it\", with H10 and H20; ruling, no run, not a try; the T1 builder is to take "
                f"{HM_OPT[0.5]}, {HM_OPT[1.0]} and {HM_OPT[2.0]}, each only with `--captain-ruling R44`)")
ROW_H = f"| {FIRST_CELL_H} | T1 | Blood's momentum |"
FIRST_CELL_M = ("R45 (the captain's yes of a later day, \"Yes, T1_s2rcB25MS, MS at 2.5\"; ruling, no run, not a try; "
                f"the T1 builder is to take {MS25_OPT} only with `--captain-ruling R45`)")
ROW_M = f"| {FIRST_CELL_M} | T1 | MS at 2.5 |"
# Line 6's two files, line 7's file and row R37's pooled file (on the shared two-library fixture), and the banked file
# and the default file on the three-library fixture, built by the builder at main ed6722a (before these options
# existed), twice in separate processes with the same result: the file's sha256 and its record's (feng_hygiene.commit
# set to "own_repo"), the same kind of synthetic-fixture pin as those of the earlier lines.
LINE6_BEFORE = {1.5: ("bc4f7c5605d462c51a5cc240c44b211dad59b916a3cc17dcc649f0d4f04d2bf7",
                      "9e9fc9b9a733a6ad38b0eeb7d3eb36cb4c0b410eb7376038b68caebc4876b08e"),
                2.0: ("d1d222a12f0688c681ff933eeb08cdc575edbb9b21272c338096b4dde83f2894",
                      "def388d9b1d931465f209abf7ffb2f9b2b043403eb87e0c8a61eefb155a646d8")}
LINE7_BEFORE = ("75d7e3c6f5f98e6551ad4459a5ae0c4aed8a20608511120f33bbf48293eb3fba",
                "a4e492a9eda4d153f957b505283032cccf67be72c5973daada5fb8386273539d")
POOLED_BEFORE = ("2fd892aa67b8dd47a6a14f5127280fed89d5bfe6f0f3d7049e5b21846b068828",
                 "83c476fd1b7e7654537bdf5f7cf4b13dc60358921255943ad67d7b9af46ef3b9")
BANKED_MS4_BEFORE = ("4b7ee589fff225035b6ca8a10f8877a3a2848cf951da026a80ee8f0cb42be348",
                     "5dba02827804fd10f9f088fa6aadec52993adb31af1307903abb537d5f78335f")
DEFAULT_MS4_BEFORE = ("c89273f066e5e27a3648eb41de6a16347a3d95f3063a64b092eef20d8ede498b",
                      "b24d95268e48ab588feb19daa2133f92ac2dcdccbcb72f1a25ccd509808090ce")
# Numbers a refusal may show: the options' values, lib3's, the rows' IDs and the files' names (T1_s2rcB25H05 and
# T1_s2rcB25MSH05 read as 1, 2, 25 and 05), besides what the command line itself said
NUMBERS_ALLOWED = {"0.5", "1.0", "2.0", "2.5", "1.5", "3.0", "1", "2", "3", "05", "25", "13", "20", "22", "23", "26",
                   "37", "44", "45"}


def _sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def _ledger_text(h: str | None = ROW_H, m: str | None = ROW_M, r22: str | None = ROW_R22_L67) -> str:
    return "# ledger\n" + "".join(f"{row}\n" for row in (ROW_R20, r22, ROW_R23, ROW_R26, ROW_R37, h, m) if row)


@pytest.fixture
def ledger(bt, tmp_path, monkeypatch):
    """A temporary ledger holding the ruling rows R20, R22 (lines 1 to 7), R23, R26 and R37 and the synthetic rows R44
    (Blood's momentum at the three values) and R45 (MS at 2.5); the builder reads it in place of the repository's."""
    p = tmp_path / "TRIES.md"
    p.write_text(_ledger_text())
    monkeypatch.setattr(bt, "TRIES_FILE", p)
    return p


@pytest.fixture
def ms4(bt, t1_env, monkeypatch):  # noqa: F811
    """The four-group recomputed setup on E10.5 names ending in _1, _2 and _3 in turn (the fit needs libraries 1 and
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


def _labels(pop) -> np.ndarray:
    return pop.obs["celltype"].astype(str).to_numpy()


def _hm_steps(bt, monkeypatch) -> list[dict]:
    """Around the builder's held-momentum step: per build, its arguments, its record and what the step gets and
    returns."""
    seen = []
    inner = bt._t1_hm_step

    def capture(ref, library, damping, held, seed):
        entry = {"library": np.array(library), "damping": damping, "held": list(held), "seed": seed}
        step, record = inner(ref, library, damping, held, seed)
        entry["record"] = record

        def wrapped(pop):
            entry["in"] = pop.copy()
            out = step(pop)
            entry["out"] = out.copy()
            return out

        seen.append(entry)
        return wrapped, record

    monkeypatch.setattr(bt, "_t1_hm_step", capture)
    return seen


def _shifts(bt, monkeypatch) -> list[dict]:
    """Around t1_scope.momentum_shift: per call, its inputs' cell names and labels, its damping, keywords and report."""
    calls = []
    inner = bt.t1_scope.momentum_shift

    def spy(start, ref, target_like, damping, held_flat=(), **k):
        out, rep = inner(start, ref, target_like, damping, held_flat, **k)
        calls.append({"start": list(start.obs_names), "ref": list(ref.obs_names),
                      "target": list(target_like.obs_names), "target_labels": _labels(target_like).tolist(),
                      "damping": damping, "held": list(held_flat), "kw": dict(k), "report": rep.as_dict()})
        return out, rep

    monkeypatch.setattr(bt.t1_scope, "momentum_shift", spy)
    return calls


def _hooks(bt, monkeypatch) -> list[dict]:
    """Around the builder's g_mom: per build, its keyword arguments and the hook it returns."""
    seen = []
    inner = bt._t1_g_mom

    def capture(*a, **k):
        hook, record = inner(*a, **k)
        seen.append({"kw": dict(k), "hook": hook})
        return hook, record

    monkeypatch.setattr(bt, "_t1_g_mom", capture)
    return seen


def _factors_one(bt, monkeypatch) -> list:
    """Every momentum_shift call with its damping forced to 0: log f = 0, so every factor is exactly 1."""
    calls = []
    inner = bt.t1_scope.momentum_shift

    def one(start, ref, target_like, damping, *a, **k):
        calls.append(damping)
        return inner(start, ref, target_like, 0.0, *a, **k)

    monkeypatch.setattr(bt.t1_scope, "momentum_shift", one)
    return calls


def _known() -> set[str]:
    """Every file hash pinned for an existing command."""
    return {*(b[1] for b in BEFORE.values()), R20_BEFORE[0], *(b[0] for b in LINE1_BEFORE.values()),
            *(b[0] for b in LINE2_BEFORE.values()), *(b[0] for b in LINE3_BEFORE.values()),
            *(b[0] for b in LINE5_BEFORE.values()), *(b[0] for b in LINE6_BEFORE.values()), LINE7_BEFORE[0],
            POOLED_BEFORE[0], BANKED_MS4_BEFORE[0], DEFAULT_MS4_BEFORE[0], MS_BEFORE[0], MSB_BEFORE[0]}


# ------------------------------------------------------------------ the constants, the help and the ruling check


def test_the_options_are_off_by_default_and_have_their_own_constants(bt):
    a = _parse(bt, ["t1", *RC])
    assert a.held_momentum is None and a.direction_momentum is None and a.captain_ruling is None
    assert bt.T1_HM_CHOICES == ("0.5", "1.0", "2.0") and bt.T1_HM_OPTIONS == HM_OPT
    assert bt.T1_HM_FILES == NAMES and bt.T1_BANKED_FILE == "T1_s2rcB25"
    assert bt.T1_BANKED_MULT == 2.5 == bt.FENG_FRACTION_MULTS_R22[0]  # the banked file's multiplier, no new value
    assert (bt.T1_HM_LABEL, bt.T1_HM_CATEGORY) == ("Blood", "blood")
    assert bt.T1_HM_STREAM == 6 and bt.T1_HM_STREAM not in (bt.T1_MS_STREAM, bt.T1_MS_DIRECTION_STREAM, 5,
                                                            bt.T1_MID_STREAM)
    sig = inspect.signature(bt.t1_scope.momentum_shift).parameters  # _t1_ms_step's calls take these defaults
    assert (bt.T1_HM_EPS, bt.T1_HM_MIN_CELLS) == (sig["eps"].default, sig["min_cells"].default) == (0.5, 20)
    assert sig["keep_depth"].default is True
    assert bt.T1_HM_COMMAND.format(value="0.5", row="R44") == " ".join(["t1", *RC, *HM[0.5]])
    assert bt.T1_MS25_OPTION == MS25_OPT and bt.T1_MS25_FILE == "T1_s2rcB25MS"
    assert bt.T1_MS25_COMMAND.format(row="R45") == " ".join(["t1", *RC, *MS25])
    assert bt.T1_RULING_ROW_ID == r"R[0-9]+[a-z]?"
    # the earlier options' constants are unchanged
    assert bt.G_MOM_DAMPINGS == (0.5, 1.0) and bt.FENG_FRACTION_MULTS == (1.0, 1.5)
    assert bt.DIRECTION_MOMENTUM_ARMS == ("all",) and bt.T1_R23_OPTION == "`--direction-momentum all`"
    assert bt.T1_MS_COMMAND == " ".join(["t1", *RC, *MS_ARGS]) and bt.T1_MSB_COMMAND == " ".join(["t1", *RC, *MSB])
    assert bt.T1_MS_FIT_LIBRARIES == ("1", "2") and bt.T1_DISSECTED_LIBRARY == "3"
    # the labels the step reads: Blood and Macrophage share category blood (so the fit restricts to Blood), Early VIC
    # is endocardium's
    spec = bt.t1_scope.load_spec()
    assert spec.labels["Blood"] == spec.labels["Macrophage"] == "blood" and spec.labels["Early VIC"] == "endocardium"
    assert "blood" not in bt.G_MOM_EDIT_CATEGORIES  # g_mom leaves Blood alone


def test_the_help_and_the_docstring_name_the_new_files(bt):
    import io
    from contextlib import redirect_stdout

    buf = io.StringIO()
    with redirect_stdout(buf), pytest.raises(SystemExit):
        bt.build_parser().parse_args(["t1", "--help"])
    text = re.sub(r"-\s+", "-", " ".join(buf.getvalue().split()))  # the help's lines may break after a hyphen
    for word in ("--held-momentum {0.5,1.0,2.0}", "T1_s2rcB25H05", "T1_s2rcB25H10", "T1_s2rcB25H20", "T1_s2rcB25MS",
                 "the banked file T1_s2rcB25", "with --held-momentum, or with --feng-fraction-mult 2.5"):
        assert word in text, word
    doc = " ".join(bt.__doc__.split())
    for word in ("--held-momentum 0.5 | 1.0 | 2.0", "T1_s2rcB25H05", "T1_s2rcB25MS", 'within={"blood": ["Blood"]}',
                 "[seed, 6]", "library 3 never in the fit", "T1_s2rcB25MSH05, is named but not coded"):
        assert word in doc, word


@pytest.mark.parametrize("text, ok", [
    (ROW_H, True),
    (f"| R44 (the captain's order of X; ruling, no run; {HM_OPT[0.5]} only with `--captain-ruling R44`) | T1 |", True),
    (f"| R44 (the captain's named exception of X; ruling, no run; {HM_OPT[0.5]} `--captain-ruling R44`) |", True),
    (f"| R44 (the captain's words of X; ruling, no run; {HM_OPT[0.5]} `--captain-ruling R44`) |", True),
    ("| R44 (the captain's yes of X; ruling, no run; --held-momentum 0.5 `--captain-ruling R44`) |", False),
    ("| R44 (the captain's yes of X; ruling, no run; `--held-momentum 0.50` `--captain-ruling R44`) |", False),
    ("| R44 (the captain's yes of X; ruling, no run; `--feng-fraction-mult 2.5 --held-momentum 0.5` "
     "`--captain-ruling R44`) |", False),  # not the option string as the file's command gives it
    (f"| R44 (the captain's yes of X; ruling, no run; {HM_OPT[0.5]} `--captain-ruling R22`) |", False),
    (f"| R44 (the captain's yes of X; ruling, no run; {HM_OPT[0.5]}) |", False),
    (f"| R44 (the captain's yes of X; a try; {HM_OPT[0.5]} `--captain-ruling R44`) |", False),
    (f"| R44 awaits (the captain's yes; ruling, no run; {HM_OPT[0.5]} `--captain-ruling R44`) |", False),
    (f"| R44 (a proposal; ruling, no run; {HM_OPT[0.5]} `--captain-ruling R44`) |", False),
    (f"| R44 (the captain's yes of X; ruling, no run) | T1 | {HM_OPT[0.5]} `--captain-ruling R44` |", False),
    (f"| R440 (the captain's yes of X; ruling, no run; {HM_OPT[0.5]} `--captain-ruling R44`) |", False),
    (f"| R4 (the captain's yes of X; ruling, no run; {HM_OPT[0.5]} `--captain-ruling R44`) |", False),
])
def test_the_ruling_check_takes_only_a_single_ruling_row_naming_the_option(bt, tmp_path, text, ok):
    p = tmp_path / "TRIES.md"
    p.write_text(f"# ledger\n{text}\n")
    assert bt._captain_ruling_row("R44", HM_OPT[0.5], p) is ok
    p.write_text(f"# ledger\n{text}\n{text}\n")  # two rows with one ID: never
    assert not bt._captain_ruling_row("R44", HM_OPT[0.5], p)


def test_the_ruling_check_takes_only_a_ruling_rows_id_and_each_option_alone(bt, ledger, tmp_path):
    for v, opt in HM_OPT.items():
        assert bt._captain_ruling_row("R44", opt, ledger), v
        assert not bt._captain_ruling_row("R45", opt, ledger) and not bt._captain_ruling_row("R22", opt, ledger)
    assert bt._captain_ruling_row("R45", MS25_OPT, ledger) and not bt._captain_ruling_row("R44", MS25_OPT, ledger)
    for opt in ("`--direction-momentum all`", "`--feng-fraction-mult 2.5`"):  # the pair's parts are not the pair
        assert not bt._captain_ruling_row("R45", opt, ledger)
    assert not bt._captain_ruling_row("R23", MS25_OPT, ledger) and not bt._captain_ruling_row("R26", MS25_OPT, ledger)
    p = tmp_path / "ids.md"
    for rid in ("r44", "T64", "R44x1", "R", "44", "R44 ", " R44", "R-44", "R44ab", "R٤٤", "", None, 44):
        text = f"| {rid} (the captain's yes of X; ruling, no run; {HM_OPT[0.5]} `--captain-ruling {rid}`) |"
        p.write_text(f"# ledger\n{text}\n")
        assert not bt._captain_ruling_row(rid, HM_OPT[0.5], p), rid
    p.write_text(f"# ledger\n| R3b (the captain's yes of X; ruling, no run; {HM_OPT[0.5]} `--captain-ruling R3b`) |\n")
    assert bt._captain_ruling_row("R3b", HM_OPT[0.5], p)  # an ID with its letter, as R3b's
    # a row naming one value orders that value alone
    p.write_text(f"# ledger\n| R44 (the captain's yes of X; ruling, no run; {HM_OPT[0.5]} `--captain-ruling R44`) |\n")
    assert bt._captain_ruling_row("R44", HM_OPT[0.5], p)
    assert not bt._captain_ruling_row("R44", HM_OPT[1.0], p) and not bt._captain_ruling_row("R44", HM_OPT[2.0], p)


def test_the_bet_refusal_lets_only_the_ruled_commands_through(bt, ledger):
    def ns(**k):
        base = dict(released_bet=None, bet_strength=None, captain_ruling="R44", feng_fraction_mult=2.5,
                    g_mom_damping=None, direction_momentum=None, held_momentum="0.5", feng_groups="recomputed",
                    seed=SEED, board=None, base="lib3")
        return SimpleNamespace(**{**base, **k})

    ms = dict(held_momentum=None, direction_momentum="all", captain_ruling="R45")
    goods = (ns(), ns(held_momentum="1.0"), ns(held_momentum="2.0"), ns(**ms))
    for good in goods:
        assert bt._bet_refusals(good) is None
    assert all(bt._hm_ordered(g) for g in goods[:3]) and bt._ms25_ordered(goods[3])
    bad = (ns(captain_ruling="R45"), ns(held_momentum="0.75"), ns(held_momentum="1"), ns(bet_strength=0.5),
           ns(captain_ruling="R22"), ns(**{**ms, "captain_ruling": "R44"}), ns(**ms, g_mom_damping=0.5),
           ns(**{**ms, "feng_fraction_mult": 1.5}), ns(**{**ms, "held_momentum": "0.5"}))
    for b in bad:
        assert not (bt._hm_ordered(b) or bt._ms25_ordered(b)), b
    for b in bad[:4]:  # (with R22 the namespace is the banked command's, which _t1_hm refuses before this guard)
        with pytest.raises(SystemExit, match="with --released-bet only"):
            bt._bet_refusals(b)
    ledger.write_text(_ledger_text(h=None, m=None))  # the ledger without the two rows: the shortcuts close
    for good in goods:
        assert not (bt._hm_ordered(good) or bt._ms25_ordered(good))
        with pytest.raises(SystemExit, match="with --released-bet only"):
            bt._bet_refusals(good)


def test_each_new_command_is_taken_with_its_ruling_and_a_ledger_naming_it(bt, ledger, monkeypatch):
    monkeypatch.setattr(fr, "committed_at_head", lambda p: json.dumps({"route": "recomputed"}).encode())
    for v in HM_OPT:
        a = _parse(bt, ["t1", *RC, *HM[v]])
        hm = bt._t1_hm(a)
        assert hm == {"held_momentum": v, "captain_ruling": "R44", "g_mom_damping": 0.5, "feng_fraction_mult": 2.5,
                      "file": NAMES[v], "order": bt.T1_HM_ORDER.format(row="R44", option=HM_OPT[v], value=f"{v:.1f}"),
                      "command": " ".join(["t1", *RC, *HM[v]])}
        assert bt._t1_candidate(a) == {k: hm[k] for k in ("g_mom_damping", "feng_fraction_mult", "captain_ruling",
                                                          "order")}
        assert bt._t1_ms(a) is None and bt._t1_pooled(a) is None and bt._bet_refusals(a) is None
        assert bt._feng_recomputed_refusals(a) is None
        for word in ("row R44", "the captain's ruling naming", HM_OPT[v], "T1_s2rcB25", f"at damping {v:.1f}",
                     "R1 reversed for Blood in this one file", "Macrophage and Early VIC held", "nothing else changed"):
            assert word in hm["order"], word
    a = _parse(bt, ["t1", *RC, *MS25])
    ms = bt._t1_ms(a)
    assert ms == {"direction_momentum": "all", "captain_ruling": "R45", "g_mom_damping": 0.5, "feng_fraction_mult": 2.5,
                  "file": "T1_s2rcB25MS", "order": bt.T1_MS25_ORDER.format(row="R45"),
                  "command": " ".join(["t1", *RC, *MS25])}
    assert bt._t1_candidate(a) == {k: ms[k] for k in ("g_mom_damping", "feng_fraction_mult", "captain_ruling", "order")}
    assert bt._t1_hm(a) is None and bt._t1_pooled(a) is None and bt._bet_refusals(a) is None
    assert bt._feng_recomputed_refusals(a) is None
    for word in ("row R45", MS25_OPT, "candidate MS's step exactly as row R23 fixes it", "T1_s2rcB25",
                 "nothing else changed"):
        assert word in ms["order"], word
    # one row naming every option orders every file under one ID
    ledger.write_text(_ledger_text(h=ROW_H.replace("each only with `--captain-ruling R44`)",
                                                   f"and {MS25_OPT}, each only with `--captain-ruling R44`)"), m=None))
    a = _parse(bt, ["t1", *RC, *MS25[:-1], "R44"])
    assert bt._t1_ms(a)["captain_ruling"] == "R44"
    assert bt._t1_hm(_parse(bt, ["t1", *RC, *HM[2.0]]))["held_momentum"] == 2.0


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
        assert "held_momentum" not in json.dumps(rec) and "T1_s2rcB25MS" not in json.dumps(rec), name


def test_ms_the_combined_ms_file_the_banked_file_and_the_default_give_the_same_bytes_on_three_libraries(bt, ms4,
                                                                                                       ledger,
                                                                                                       tmp_path):
    for name, argv, before in (("T1_s2rc_MS", MS_ARGS, MS_BEFORE), ("T1_s2rcMSB", MSB, MSB_BEFORE),
                               ("T1_s2rcB25", BANKED, BANKED_MS4_BEFORE), ("T1_s2rc", [], DEFAULT_MS4_BEFORE)):
        f, rec = _build(bt, tmp_path, name, argv)
        assert (_sha(f), _record_sha(rec)) == before, name
    f, rec = _build(bt, tmp_path, "T1_s2rc_MS2", MS_ARGS)
    assert rec["t1_candidate"]["command"] == bt.T1_MS_COMMAND and "file" not in rec["t1_candidate"]


def test_without_the_options_no_new_step_runs(bt, ms4, ledger, monkeypatch, tmp_path):
    hm = _hm_steps(bt, monkeypatch)
    calls = _shifts(bt, monkeypatch)
    for name, argv in (("T1_s2rc", []), ("T1_s2rcB25", BANKED), ("T1_s2rcB25d15", LINE6[1.5])):
        _build(bt, tmp_path, name, argv)
    assert hm == [] and calls == []


# ------------------------------------------------------------------ Blood's momentum: the H files


@pytest.mark.parametrize("v", list(HM_OPT))
def test_each_h_file_is_the_banked_file_but_for_the_blood_rows_and_builds_the_same_bytes_twice(bt, ms4, ledger,
                                                                                               tmp_path, v):
    f0, r0 = _build(bt, tmp_path, "T1_s2rcB25", BANKED)
    fa, ra = _build(bt, tmp_path, f"{NAMES[v]}a", HM[v])
    fb, rb = _build(bt, tmp_path, f"{NAMES[v]}b", HM[v])
    assert fa.read_bytes() == fb.read_bytes() and _record_sha(ra) == _record_sha(rb)
    assert _sha(f0) == BANKED_MS4_BEFORE[0]
    x0, xh = _same_frame(f0, fa)  # the same cells, order, obs and var as the banked file
    pop = ms4[-1]["before"]  # the population the H build's post_cells gets: the banked file's at g_mom's entry
    labels = _labels(pop)
    blood = labels == "Blood"
    assert blood.any() and not blood.all() and {"Macrophage", "Early VIC"} <= set(labels)
    assert x0[~blood].tobytes() == xh[~blood].tobytes()  # every non-Blood row: the banked file's bytes
    changed = (x0 != xh).any(axis=1)
    assert changed.any() and not (changed & ~blood).any()  # only Blood rows move, and some do
    b0 = ms4[0]["before"]  # the draw and the atlas step are the banked file's
    assert list(b0.obs_names) == list(pop.obs_names) and b0.obs.equals(pop.obs)
    assert np.array_equal(_dense(b0.layers["counts"]), _dense(pop.layers["counts"]))
    assert np.array_equal(_dense(b0.X), _dense(pop.X))
    assert ms4[0]["fraction"] == ms4[-1]["fraction"] == 2.5 * 0.4 and ms4[0]["directions"] == ms4[-1]["directions"]
    assert _sha(fa) not in _known()


def test_the_three_values_give_three_files_and_the_banked_files_g_mom_and_fraction(bt, ms4, ledger, tmp_path,
                                                                                  monkeypatch):
    hooks = _hooks(bt, monkeypatch)
    f0, _ = _build(bt, tmp_path, "T1_s2rcB25", BANKED)
    files = {v: _build(bt, tmp_path, NAMES[v], HM[v])[0] for v in HM_OPT}
    assert len({f.read_bytes() for f in (f0, *files.values())}) == 4
    assert all(h["kw"] == hooks[0]["kw"] for h in hooks)  # g_mom: the banked file's damping 0.5 and multiplier
    assert hooks[0]["kw"]["damping"] == 0.5 and hooks[0]["kw"]["fraction_mult"] == 2.5
    x0 = _dense(ad.read_h5ad(f0).X)
    blood = _labels(ms4[0]["before"]) == "Blood"
    moved = [float(np.abs(_dense(ad.read_h5ad(files[v]).X)[blood] - x0[blood]).mean()) for v in HM_OPT]
    assert 0 < moved[0] < moved[1] < moved[2]  # a larger damping moves Blood further from its drawn expression


@pytest.mark.parametrize("v", list(HM_OPT))
def test_with_bloods_factors_at_1_the_h_file_is_the_banked_file(bt, ms4, ledger, tmp_path, monkeypatch, v):
    f0, _ = _build(bt, tmp_path, "T1_s2rcB25", BANKED)
    calls = _factors_one(bt, monkeypatch)
    f1, r1 = _build(bt, tmp_path, f"{NAMES[v]}_one", HM[v])
    assert calls == [v]  # one call, Blood's, at the file's damping, forced to 0 here
    assert f1.read_bytes() == f0.read_bytes()
    assert r1["t1_candidate"]["momentum"]["status"] == "shifted" and r1["t1_candidate"]["momentum"]["damping"] == v


def test_the_step_fits_on_blood_cells_of_e9_5_and_libraries_1_and_2_and_moves_the_drawn_blood_cells_alone(
        bt, ms4, t1_env, ledger, tmp_path, monkeypatch):  # noqa: F811
    steps = _hm_steps(bt, monkeypatch)
    calls = _shifts(bt, monkeypatch)
    _build(bt, tmp_path, "T1_s2rcB25H10", HM[1.0])
    assert len(calls) == 1 and len(steps) == 1
    c, s = calls[0], steps[0]
    assert c["damping"] == 1.0 and s["damping"] == 1.0 and s["seed"] == SEED
    assert c["start"] == list(t1_env.start.obs_names)  # E9.5, every cell (the fit restricts to Blood)
    lib = np.array([n.rsplit("_", 1)[1] for n in t1_env.ref.obs_names])
    assert np.array_equal(lib, s["library"])
    assert c["ref"] == [n for n, x in zip(t1_env.ref.obs_names, lib, strict=True) if x in ("1", "2")]  # no library 3
    pop = s["in"]
    blood = _labels(pop) == "Blood"
    assert c["target"] == list(pop.obs_names[blood]) and set(c["target_labels"]) == {"Blood"}
    spec = bt.t1_scope.load_spec()
    assert c["kw"] == {"label_key": "celltype", "eps": 0.5, "min_cells": 20, "keep_depth": True,
                       "within": {"blood": ["Blood"]}, "seed": [SEED, 6],
                       "start_fallback": spec.fallbacks("E9.5_RNA"), "ref_fallback": spec.fallbacks("E10.5_RNA")}
    assert c["held"] == t1_env.held_calls[-1] == s["held"]  # the build's held list: R2's union and the pre-registered
    rows = c["report"]["rows"]
    assert [r["type"] for r in rows] == ["blood"] and rows[0]["status"] == "shifted"
    assert (rows[0]["start_cells"], rows[0]["ref_cells"]) == ("own, restricted", "own, restricted")
    ref_labels = t1_env.ref.obs["celltype"].astype(str).to_numpy()
    assert rows[0]["n_start"] == E95["Blood"]  # E9.5's Blood cells
    assert rows[0]["n_ref"] == int(((ref_labels == "Blood") & np.isin(lib, ["1", "2"])).sum())  # no Macrophage
    assert int(((ref_labels == "Macrophage") & np.isin(lib, ["1", "2"])).sum()) > 0  # which the category holds
    # the step's own output: Blood rows moved by that call, every other row as it came
    out = s["out"]
    assert list(out.obs_names) == list(pop.obs_names) and out.obs.equals(pop.obs)
    assert np.array_equal(_dense(out.layers["counts"])[~blood], _dense(pop.layers["counts"])[~blood])
    assert np.array_equal(_dense(out.X)[~blood], _dense(pop.X)[~blood])
    assert not np.array_equal(_dense(out.layers["counts"])[blood], _dense(pop.layers["counts"])[blood])


def test_the_step_is_the_released_momentum_computed_directly_then_g_mom(bt, ms4, t1_env, ledger, tmp_path,  # noqa: F811
                                                                       monkeypatch):
    steps = _hm_steps(bt, monkeypatch)
    hooks = _hooks(bt, monkeypatch)
    _build(bt, tmp_path, "T1_s2rcB25H05", HM[0.5])
    s = steps[0]
    pop = s["in"]
    spec = bt.t1_scope.load_spec()
    start = harness._shallow(t1_env.start.copy())
    bt.t1_scope.classify(start, source="E9.5_RNA", spec=spec)
    later = t1_env.ref[np.flatnonzero(np.isin(s["library"], ["1", "2"]))].copy()
    bt.t1_scope.classify(later, source="E10.5_RNA", spec=spec)
    idx = np.flatnonzero(_labels(pop) == "Blood")
    want, rep = bt.t1_scope.momentum_shift(start, later, pop[idx].copy(), 0.5, s["held"], label_key="celltype",
                                           within={"blood": ["Blood"]}, seed=[SEED, 6],
                                           start_fallback=spec.fallbacks("E9.5_RNA"),
                                           ref_fallback=spec.fallbacks("E10.5_RNA"))
    got = s["out"][idx]
    assert rep.rows[0]["status"] == "shifted"
    assert np.array_equal(_dense(want.layers["counts"]), _dense(got.layers["counts"]))
    assert np.array_equal(_dense(want.X), _dense(got.X))
    # depth kept in expectation: the momentum's factors are scaled per cell to keep its expected total
    assert rep.rows[0]["median_depth_scale"] != 1.0 and rep.keep_depth is True
    # then g_mom, the banked file's hook, on the step's output; it leaves Blood's counts as the step left them
    final = ms4[-1]["after"]
    again = hooks[-1]["hook"](s["out"].copy())
    assert np.array_equal(_dense(again.layers["counts"]), _dense(final.layers["counts"]))
    assert np.array_equal(_dense(final.layers["counts"])[idx], _dense(got.layers["counts"]))


def test_the_h_record_is_the_banked_files_plus_its_candidate_field(bt, ms4, ledger, tmp_path):
    _, r0 = _build(bt, tmp_path, "T1_s2rcB25", BANKED)
    for v in HM_OPT:
        _, rh = _build(bt, tmp_path, NAMES[v], HM[v])
        assert {k: x for k, x in rh.items() if k != "t1_candidate"} == {k: x for k, x in r0.items()
                                                                       if k != "t1_candidate"}
        cand = rh["t1_candidate"]
        assert set(cand) == {"held_momentum", "captain_ruling", "g_mom_damping", "feng_fraction_mult", "file", "order",
                             "command", "momentum", "calendar_fraction", "applied_fraction"}
        assert (cand["held_momentum"], cand["captain_ruling"], cand["file"]) == (v, "R44", NAMES[v])
        assert (cand["g_mom_damping"], cand["feng_fraction_mult"]) == (0.5, 2.5) == (
            r0["t1_candidate"]["g_mom_damping"], r0["t1_candidate"]["feng_fraction_mult"])
        assert cand["command"] == " ".join(["t1", *RC, *HM[v]])
        assert cand["calendar_fraction"] == 0.4 and cand["applied_fraction"] == rh["fraction"] == r0["fraction"] == 1.0
        m = cand["momentum"]
        assert m["fit"] == {"start": "E9.5_RNA", "later": "E10.5_RNA", "later_libraries": ["1", "2"],
                            "library_3_in_fit": False, "within": {"blood": ["Blood"]},
                            "start_cells": "own, restricted", "later_cells": "own, restricted"}
        assert (m["label"], m["category"], m["status"], m["damping"]) == ("Blood", "blood", "shifted", v)
        assert (m["eps"], m["min_cells"], m["keep_depth"], m["stream"]) == (0.5, 20, True, [SEED, 6])
        assert "R2's union" in m["held"] and "pre-registered" in m["held"] and "T64" in m["practice_arm"]
        assert "before g_mom" in m["rule"] and "Macrophage and Early VIC keep R1's hold" in m["rule"]
        assert not re.search(r"n_target|n_start|n_ref|median_depth", json.dumps(cand))  # no cell number


def test_the_step_is_refused_before_the_draw_when_it_would_be_idle(bt, ms4, ledger, tmp_path, monkeypatch):
    out = tmp_path / "T1_s2rcB25H05.h5ad"
    monkeypatch.setattr(bt, "T1_HM_MIN_CELLS", 10 ** 6)  # Blood's own cells below min_cells at an end of the fit
    with pytest.raises(SystemExit, match="Blood's own cells at the fit's start stage do not reach min_cells"):
        bt.main(["t1", *RC, *HM[0.5], "--out", str(out)])
    assert ms4 == [] and not out.exists()  # before the draw: build_submission never ran


def test_the_step_is_refused_unless_the_portion_libraries_are_1_and_2(bt, rc4, ledger, tmp_path):
    """The shared fixture's E10.5 names end in _1 and _3 only: the fit's libraries are not there, and the build stops
    before the draw."""
    out = tmp_path / "T1_s2rcB25H05.h5ad"
    with pytest.raises(SystemExit, match=r"libraries other than 3 are \['1'\], not the fit's \['1', '2'\]"):
        bt.main(["t1", *RC, *HM[0.5], "--out", str(out)])
    assert rc4 == [] and not out.exists()


def test_the_step_never_moves_a_cell_the_banked_file_moves_when_epicardium_has_no_direction(bt, ms3, ledger, tmp_path):
    f0, _ = _build(bt, tmp_path, "T1_s2rcB25", BANKED)
    fh, _ = _build(bt, tmp_path, "T1_s2rcB25H20", HM[2.0])
    x0, xh = _same_frame(f0, fh)
    blood = _labels(ms3[0]["before"]) == "Blood"
    assert x0[~blood].tobytes() == xh[~blood].tobytes() and (x0[blood] != xh[blood]).any()


# ------------------------------------------------------------------ MS on the banked command: T1_s2rcB25MS


@pytest.mark.parametrize("fixture", ["ms4", "ms3"])
def test_the_ms25_file_is_ms_on_the_banked_command_and_builds_the_same_bytes_twice(bt, ledger, tmp_path, request,
                                                                                  monkeypatch, fixture):
    seen = request.getfixturevalue(fixture)
    with_groups = GROUPS4 if fixture == "ms4" else THREE
    hooks = _hooks(bt, monkeypatch)
    f0, r0 = _build(bt, tmp_path, "T1_s2rcB25", BANKED)
    fa, ra = _build(bt, tmp_path, "T1_s2rcB25MSa", MS25)
    fb, rb = _build(bt, tmp_path, "T1_s2rcB25MSb", MS25)
    assert fa.read_bytes() == fb.read_bytes() and _record_sha(ra) == _record_sha(rb)
    x0, xm = _same_frame(f0, fa)
    pop = seen[-1]["before"]
    labels = _labels(pop)
    groups = bt.feng_dirs.label_groups(sorted(set(labels)))
    covered = np.array([groups.get(lab) in with_groups for lab in labels])
    assert covered.any() and not covered.all()
    changed = (x0 != xm).any(axis=1)
    assert changed.any() and not (changed & ~covered).any()
    assert x0[~covered].tobytes() == xm[~covered].tobytes()  # Blood and every cell without a group: the banked bytes
    # MS takes the atlas step inside its own step, at the banked file's fraction; g_mom as the banked file's
    assert seen[0]["fraction"] == 2.5 * 0.4 and seen[-1]["fraction"] is None and seen[-1]["directions"] == []
    assert hooks[0]["kw"] == hooks[-1]["kw"] and hooks[-1]["kw"]["damping"] == 0.5
    assert _sha(fa) not in _known()


def test_with_ms_factors_at_1_the_ms25_file_is_the_banked_file(bt, ms4, ledger, tmp_path, monkeypatch):
    f0, _ = _build(bt, tmp_path, "T1_s2rcB25", BANKED)
    calls = _factors_one(bt, monkeypatch)
    f1, _ = _build(bt, tmp_path, "T1_s2rcB25MS_one", MS25)
    assert calls == [0.5] and f1.read_bytes() == f0.read_bytes()  # MS's call, B0's damping, forced to 0 here


def test_the_ms25_step_is_r23s_step_with_the_banked_fraction(bt, ms4, ledger, tmp_path, monkeypatch):
    seen = []
    inner = bt._t1_ms_step

    def capture(ref, library, directions, fraction, groups, held, seed):
        seen.append({"fraction": fraction, "seed": seed, "directions": sorted(directions)})
        return inner(ref, library, directions, fraction, groups, held, seed)

    monkeypatch.setattr(bt, "_t1_ms_step", capture)
    _build(bt, tmp_path, "T1_s2rc_MS", MS_ARGS)
    _build(bt, tmp_path, "T1_s2rcB25MS", MS25)
    assert [s["fraction"] for s in seen] == [0.4, 2.5 * 0.4] and all(s["seed"] == SEED for s in seen)
    assert seen[0]["directions"] == seen[1]["directions"] == sorted(GROUPS4)


def test_the_ms25_record_is_the_banked_files_plus_its_candidate_field(bt, ms4, ledger, tmp_path):
    _, r0 = _build(bt, tmp_path, "T1_s2rcB25", BANKED)
    _, rm = _build(bt, tmp_path, "T1_s2rcB25MS", MS25)
    _, r23 = _build(bt, tmp_path, "T1_s2rc_MS", MS_ARGS)
    assert {k: x for k, x in rm.items() if k != "t1_candidate"} == {k: x for k, x in r0.items() if k != "t1_candidate"}
    cand = rm["t1_candidate"]
    assert set(cand) == {"direction_momentum", "captain_ruling", "g_mom_damping", "feng_fraction_mult", "file", "order",
                         "command", "momentum", "calendar_fraction", "applied_fraction"}
    assert (cand["captain_ruling"], cand["file"], cand["g_mom_damping"], cand["feng_fraction_mult"]) == (
        "R45", "T1_s2rcB25MS", 0.5, 2.5)
    assert cand["command"] == " ".join(["t1", *RC, *MS25]) and cand["order"] == bt.T1_MS25_ORDER.format(row="R45")
    assert cand["applied_fraction"] == rm["fraction"] == 1.0 and cand["calendar_fraction"] == 0.4
    assert cand["momentum"] == r23["t1_candidate"]["momentum"]  # R23's step's record exactly


# ------------------------------------------------------------------ refusals, before anything but the ledger is read


@pytest.mark.parametrize("value", ["0.75", "1", "2", "0.50", "3.0", "0", "-0.5", "1.5", "0.25", "abc", "2.00", ".5"])
def test_any_other_held_momentum_value_is_refused(bt, ledger, monkeypatch, tmp_path, value):
    _refused(bt, monkeypatch, tmp_path, ["t1", *RC, "--feng-fraction-mult", "2.5", "--held-momentum", value,
                                         "--captain-ruling", "R44"], "it takes 0.5, 1.0, 2.0 only")


@pytest.mark.parametrize("v", list(HM_OPT))
@pytest.mark.parametrize("extra, match", [
    (["--feng-fraction-mult", "2.5"], "is the captain's ruling alone: it needs --captain-ruling <ROW>"),
    *[(["--feng-fraction-mult", "2.5", "--captain-ruling", r], f"no single row {r} whose first cell")
      for r in ("R22", "R23", "R26", "R37", "R45", "R20", "R40", "R440")],
    *[(["--feng-fraction-mult", "2.5", "--captain-ruling", r], "no single row") for r in ("r44", "T64", "R44x")],
    (["--feng-fraction-mult", "2.5", "--captain-ruling", "R44", "--direction-momentum", "all"], "named but not coded"),
    (["--feng-fraction-mult", "2.5", "--captain-ruling", "R45", "--direction-momentum", "all"], "named but not coded"),
    (["--direction-momentum", "all"], "named but not coded"),
    *[(["--feng-fraction-mult", "2.5", "--captain-ruling", "R44", "--g-mom-damping", d], "refused with --g-mom-damping")
      for d in ("0.5", "1.0", "1.5", "2.0")],
    (["--feng-fraction-mult", "2.5", "--captain-ruling", "R44", "--released-bet", "share_trend", "--bet-strength",
      "0.5"], "refused with --released-bet, --bet-strength"),
    (["--feng-fraction-mult", "2.5", "--captain-ruling", "R44", "--bet-strength", "0.5"],
     "refused with --bet-strength"),
    (["--captain-ruling", "R44"], "without --feng-fraction-mult 2.5 or with another multiplier"),
    *[(["--feng-fraction-mult", m, "--captain-ruling", "R44"], "without --feng-fraction-mult 2.5 or with another")
      for m in ("1.0", "1.5", "2.0", "3.0", "2.25")],
    (["--feng-fraction-mult", "2.5", "--captain-ruling", "R44", "--no-hold-unbracketed"],
     "refused with --no-hold-unbracketed"),
    (["--feng-fraction-mult", "2.5", "--captain-ruling", "R44", "--seed", "1000002"], "changes the banked file"),
    (["--feng-fraction-mult", "2.5", "--captain-ruling", "R44", "--board", "T1:val"], "changes the banked file"),
])
def test_every_other_use_of_the_held_momentum_is_refused(bt, ledger, monkeypatch, tmp_path, v, extra, match):
    _refused(bt, monkeypatch, tmp_path, ["t1", *RC, "--held-momentum", f"{v:.1f}", *extra], match)


@pytest.mark.parametrize("argv", [
    ["t1"],  # the authors' groups
    ["t1", "--g-mom", "--base", "lib3"],  # the first upload's command
    ["t1", "--g-mom", "--base", "lib3", "--feng-groups", "none"],  # the none file
    ["t1", "--g-mom", "--base", "mid", "--cm-share", "0.60"],  # T1's second upload
])
def test_the_held_momentum_is_refused_on_any_other_command(bt, ledger, monkeypatch, tmp_path, argv):
    _refused(bt, monkeypatch, tmp_path, [*argv, *HM[0.5]], "changes the banked file T1_s2rcB25 only")


@pytest.mark.parametrize("extra", [["--cm-share", "0.6"], ["--blood-share", "0.1"], ["--herd"], ["--no-feng-hygiene"],
                                   ["--lib3-pair"], ["--strong-null"], ["--base", "pooled"], ["--base", "mid"],
                                   ["--feng-none-ruling", "x"]])
def test_the_recomputed_commands_own_guards_hold_with_both_new_commands(bt, ledger, monkeypatch, tmp_path, extra):
    for new in (HM[0.5], HM[2.0], MS25):
        _refused(bt, monkeypatch, tmp_path, ["t1", *RC, *extra, *new],
                 "--feng-groups recomputed builds one file: the first upload's command")
    _refused(bt, monkeypatch, tmp_path, ["t1", "--base", "lib3", "--feng-groups", "recomputed", *HM[1.0]],
             "--feng-groups recomputed builds one file")  # without --g-mom


def test_a_held_momentum_value_is_refused_while_its_row_does_not_name_it(bt, ledger, monkeypatch, tmp_path):
    rows = [None,  # no row R44
            f"| R44 (the captain's yes of X; ruling, no run; {HM_OPT[0.5]}) |",  # no ruling flag
            f"| R44 (the captain's yes of X; ruling, no run; {HM_OPT[0.5]} `--captain-ruling R22`) |",
            f"| R44 (awaits the captain's yes; ruling, no run; {HM_OPT[0.5]} `--captain-ruling R44`) |",
            f"| R44 (the captain's yes of X; a try; {HM_OPT[0.5]} `--captain-ruling R44`) |",
            f"| R44 (the captain's yes of X; ruling, no run) | T1 | {HM_OPT[0.5]} `--captain-ruling R44` |",
            f"{ROW_H}\n{ROW_H}"]  # two rows with one ID
    for row in rows:
        ledger.write_text(_ledger_text(h=row))
        msg = _refused(bt, monkeypatch, tmp_path, ["t1", *RC, *HM[0.5]], "is refused: experiments/harness/TRIES.md "
                                                                         "has no single row R44")
        assert HM_OPT[0.5] in msg and "`--captain-ruling R44`" in msg and "captain's ruling" in msg
    # a row naming the first value orders it alone
    ledger.write_text(_ledger_text(h=f"| R44 (the captain's yes of X; ruling, no run; {HM_OPT[0.5]} only with "
                                     "`--captain-ruling R44`) |"))
    assert bt._t1_hm(_parse(bt, ["t1", *RC, *HM[0.5]]))["held_momentum"] == 0.5
    for v in (1.0, 2.0):
        _refused(bt, monkeypatch, tmp_path, ["t1", *RC, *HM[v]], re.escape(HM_OPT[v]))


@pytest.mark.parametrize("extra, match", [
    (["--g-mom-damping", "0.5"], "refused with --g-mom-damping: row R45's file is candidate MS's step on the banked"),
    (["--g-mom-damping", "1.0"], "refused with --g-mom-damping: row R45's file"),
    (["--g-mom-damping", "1.5"], "refused with --g-mom-damping: row R45's file"),
    (["--bet-strength", "0.5"], "refused with --bet-strength: row R45's file"),
    (["--released-bet", "share_trend", "--bet-strength", "0.5"], "refused with --released-bet, --bet-strength"),
    (["--no-hold-unbracketed"], "refused with --no-hold-unbracketed: the cells without a group keep R1's hold"),
    (["--seed", "1000002"], "changes the banked file T1_s2rcB25 only"),
    (["--board", "T1:val"], "changes the banked file T1_s2rcB25 only"),
    (["--held-momentum", "0.5"], "named but not coded"),
    (["--held-momentum", "2.0"], "named but not coded"),
])
def test_every_other_use_of_ms_at_2_5_under_its_row_is_refused(bt, ledger, monkeypatch, tmp_path, extra, match):
    _refused(bt, monkeypatch, tmp_path, ["t1", *RC, *MS25, *extra], match)


@pytest.mark.parametrize("ruling, match", [
    (None, "is the captain's order of row R23 alone and needs --captain-ruling R23"),
    ("R22", "refused with --captain-ruling R22: only row R23 orders it"),
    ("R23", "refused with --feng-fraction-mult: candidate MS is one change"),
    ("R26", "refused with --captain-ruling R26: only row R23 orders it"),
    ("R44", "refused with --captain-ruling R44: only row R23 orders it"),  # a row that names Blood's options only
    ("R46", "refused with --captain-ruling R46: only row R23 orders it"),
])
def test_ms_at_2_5_without_a_row_naming_it_is_refused_as_row_r23s_path_refuses_it(bt, ledger, monkeypatch, tmp_path,
                                                                                 ruling, match):
    argv = ["t1", *RC, "--feng-fraction-mult", "2.5", "--direction-momentum", "all"]
    _refused(bt, monkeypatch, tmp_path, argv if ruling is None else [*argv, "--captain-ruling", ruling], match)
    ledger.write_text(_ledger_text(m=None))  # without row R45 the pair with R45 is R23's to refuse, as before
    _refused(bt, monkeypatch, tmp_path, [*argv, "--captain-ruling", "R45"], "refused with --captain-ruling R45: only "
                                                                            "row R23 orders it")


def test_the_ruling_flag_alone_or_with_another_option_is_refused_as_before(bt, ledger, monkeypatch, tmp_path):
    for r in ("R44", "R45"):
        _refused(bt, monkeypatch, tmp_path, ["t1", *RC, "--captain-ruling", r], "with --released-bet only")
        _refused(bt, monkeypatch, tmp_path, ["t1", *RC, "--feng-fraction-mult", "2.5", "--captain-ruling", r],
                 f"refused with --captain-ruling {r}: only row R22 orders it")
        _refused(bt, monkeypatch, tmp_path, ["t1", *RC, "--direction-momentum", "all", "--captain-ruling", r],
                 f"refused with --captain-ruling {r}: only row R23 orders it")


# ------------------------------------------------------------------ the old builder against this one


OLD_COMMIT = "ed6722a"  # main before these options existed
OLD_SHA = "217b5f5ac5b12821fec08435d8f8c86b270b49e8539c20c17fc68868c4871256"  # its experiments/builds/build_test.py


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
    path = tmp_path / "old" / "build_test_before.py"
    path.parent.mkdir()
    path.write_bytes(res.stdout)
    spec = importlib.util.spec_from_file_location("build_test_before_hm", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


CMDS = [[], ["--g-mom", "--base", "lib3"], RC, ["--g-mom", "--base", "pooled", "--feng-groups", "recomputed"],
        ["--g-mom", "--feng-groups", "recomputed"], ["--base", "lib3", "--feng-groups", "recomputed"],
        ["--g-mom", "--base", "lib3", "--feng-groups", "none"], ["--g-mom", "--base", "mid", "--cm-share", "0.60"]]
DAMPINGS = [None, "0.5", "1.0", "1.5"]
MULTS = [None, "1.5", "2.5", "3.0"]
HELD = [None, "0.5", "2.0", "0.75"]
DMS = [None, "all"]
RULINGS = [None, "R22", "R23", "R26", "R44", "R45"]
EXTRAS = [["--seed", "1000002"], ["--board", "T1:val"], ["--released-bet", "share_trend", "--bet-strength", "0.5"],
          ["--bet-strength", "0.5"], ["--no-hold-unbracketed"], ["--cm-share", "0.6"], ["--herd"],
          ["--no-feng-hygiene"], ["--strong-null"], ["--lib3-pair"], ["--feng-none-ruling", "x"],
          ["--blood-share", "0.1"]]
EXTRA_SETS = [(None, "2.5", "0.5", None), (None, "2.5", "1.0", None), (None, "2.5", None, "all"),
              (None, "2.5", "2.0", "all")]
GRID_SETS = 8 * 4 * 4 * 4 * 2 * 6 + 4 * 3 * 12
SETS_COMPARED = 3 * GRID_SETS


def _grid() -> list[list[str]]:
    def argv(cmd, damping, mult, held, dm, ruling, extra=()):
        out = ["t1", *cmd]
        for flag, value in (("--g-mom-damping", damping), ("--feng-fraction-mult", mult), ("--held-momentum", held),
                            ("--direction-momentum", dm), ("--captain-ruling", ruling)):
            out += [] if value is None else [flag, value]
        return [*out, *extra]

    sets = [argv(*x) for x in itertools.product(CMDS, DAMPINGS, MULTS, HELD, DMS, RULINGS)]
    for (damping, mult, held, dm), ruling, extra in itertools.product(EXTRA_SETS, [None, "R44", "R45"], EXTRAS):
        sets.append(argv(RC, damping, mult, held, dm, ruling, extra))
    return sets


def _kind(a, h_row: str | None, h_named: set[str], m_row: str | None) -> str:
    """"new_h" or "new_m" for exactly one of the new commands under the ledger's rows (Blood's momentum at a value its
    row names, or MS at 2.5, each with its row's ID, on the banked command and nothing else); "touch" for any other
    set holding --held-momentum or MS at 2.5 with the row naming that pair; "other" for every other set."""
    plain = (a.g_mom and a.feng_groups == "recomputed" and a.base == "lib3" and a.seed == SEED and a.board is None
             and a.released_bet is None and a.bet_strength is None and a.hold_unbracketed and a.feng_hygiene is None
             and not a.strong_null and a.cm_share is None and a.blood_share is None and not a.lib3_pair
             and not a.herd and a.feng_none_ruling is None and a.g_mom_damping is None
             and a.feng_fraction_mult == 2.5)
    if a.held_momentum is not None:
        ok = (plain and a.direction_momentum is None and a.captain_ruling == h_row and a.held_momentum in h_named)
        return "new_h" if ok else "touch"
    if (a.direction_momentum == "all" and a.feng_fraction_mult == 2.5 and m_row is not None
            and a.captain_ruling == m_row):
        return "new_m" if plain else "touch"
    return "other"


def test_the_old_builder_and_this_one_route_every_other_argument_set_alike(bt, tmp_path, monkeypatch):
    """Side by side on the grid, under three synthetic ledgers (with rows R44 and R45; without them; with R44 naming
    the first value only and row R22 naming MS at 2.5): every set that holds no --held-momentum and is not MS at 2.5
    under a row naming that pair is refused with the same message, or reaches the same first read with the same
    candidate, MS and pooled values, in both builders; every set that holds one is refused by the old builder (which
    has no such option or refuses the pair), and by this one unless it is exactly a new command whose row names it,
    which this one takes with the banked file's values. Nothing is built: every read is stopped where it would start."""
    old = _old_builder(tmp_path)
    assert not hasattr(old, "_t1_hm") and hasattr(bt, "_t1_hm")
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
        for name in ("_t1_ms", "_t1_candidate", "_t1_pooled"):
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
            return ("refused", str(e.code))
        except _Reached as e:
            return ("taken", e.args[0], spies.get("_t1_ms"), spies.get("_t1_candidate"), spies.get("_t1_pooled"))
        raise AssertionError(f"{argv} neither refused nor read")

    def want_h(row, v):
        hm = bt._t1_hm(_parse(bt, ["t1", *RC, "--feng-fraction-mult", "2.5", "--held-momentum", v, "--captain-ruling",
                                   row]))
        return ("taken", "union", None, {k: hm[k] for k in ("g_mom_damping", "feng_fraction_mult", "captain_ruling",
                                                            "order")}, None)

    def want_m(row):
        ms = {"direction_momentum": "all", "captain_ruling": row, "g_mom_damping": 0.5, "feng_fraction_mult": 2.5,
              "file": "T1_s2rcB25MS", "order": bt.T1_MS25_ORDER.format(row=row),
              "command": bt.T1_MS25_COMMAND.format(row=row)}
        return ("taken", "union", ms, {k: ms[k] for k in ("g_mom_damping", "feng_fraction_mult", "captain_ruling",
                                                          "order")}, None)

    one_h = f"| R44 (the captain's yes of X; ruling, no run; {HM_OPT[0.5]} only with `--captain-ruling R44`) |"
    r22_m = ROW_R22_L67.replace("each only with `--captain-ruling R22`)",
                                f"and {MS25_OPT}, each only with `--captain-ruling R22`)")
    assert r22_m != ROW_R22_L67 and MS25_OPT in r22_m.split("|")[1]
    ledgers = {"after": (_ledger_text(), "R44", {"0.5", "1.0", "2.0"}, "R45"),
               "before": (_ledger_text(h=None, m=None), None, set(), None),
               "partial": (_ledger_text(h=one_h, m=None, r22=r22_m), "R44", {"0.5"}, "R22")}
    grid = _grid()
    assert len({tuple(x) for x in grid}) == len(grid) == GRID_SETS
    count = {"other": 0, "touch": 0, "new_h": 0, "new_m": 0, "taken_alike": 0, "refused_alike": 0, "taken_new": 0}
    # the old parser's usage errors on --held-momentum go to a file, not the captured output
    with open(tmp_path / "stderr.txt", "w") as err, pytest.MonkeyPatch.context() as quiet:
        quiet.setattr(sys, "stderr", err)
        for text, h_row, h_named, m_row in ledgers.values():
            ledger.write_text(text)
            for argv in grid:
                a = _parse(bt, argv)
                kind = _kind(a, h_row, h_named, m_row)
                o, n = route(old, argv), route(bt, argv)
                count[kind] += 1
                if kind == "other":
                    assert n == o, argv
                    count["taken_alike" if o[0] == "taken" else "refused_alike"] += 1
                    continue
                assert o[0] == "refused", argv  # the old builder took none of these
                if kind == "new_h":
                    assert n == want_h(h_row, a.held_momentum), argv
                    count["taken_new"] += 1
                elif kind == "new_m":
                    assert n == want_m(m_row), argv
                    count["taken_new"] += 1
                else:
                    assert n[0] == "refused", argv
    assert sum(count[k] for k in ("other", "touch", "new_h", "new_m")) == SETS_COMPARED
    assert count["taken_alike"] and count["refused_alike"] and count["touch"]
    # the new sets in the grid: R13's command with 2.5 and the values 0.5 and 2.0 under R44 (both named after, the first
    # alone partially), and the MS pair under R45 (after) and R22 (partial)
    assert (count["new_h"], count["new_m"]) == (2 + 1, 1 + 1) and count["taken_new"] == 5


# ------------------------------------------------------------------ experiments/builds/t1_hm_check.py


def _check():
    import t1_hm_check

    return t1_hm_check


def _lines(capsys) -> dict[str, str]:
    out = capsys.readouterr().out.splitlines()
    assert all(re.fullmatch(r"[A-Za-z0-9_]+ (True|False)", line) for line in out), out  # booleans only
    return dict(line.split(" ") for line in out)


def test_the_checkers_banked_command_is_the_launchers_and_its_draw_is_the_builders(bt, ms4, ledger, tmp_path,
                                                                                   monkeypatch):
    import tempfile

    chk = _check()
    assert chk.BANKED_ARGV == ["t1", *RC, *BANKED]
    _build(bt, tmp_path, "T1_s2rcB25", BANKED)
    monkeypatch.setattr(tempfile, "tempdir", str(tmp_path / "tmp"))
    (tmp_path / "tmp").mkdir()
    pop = chk.banked_population()
    b0 = ms4[0]["before"]  # the banked build's population at g_mom's entry
    assert list(pop.obs_names) == list(b0.obs_names) and pop.obs.equals(b0.obs)
    assert np.array_equal(_dense(pop.layers["counts"]), _dense(b0.layers["counts"]))
    assert len(ms4) == 2 and "after" not in ms4[1]  # stopped inside build_submission: g_mom never ran
    assert [p.name for p in (tmp_path / "tmp").rglob("*") if p.is_file()] == []  # nothing written


@pytest.mark.parametrize("new", ["H05", "H10", "H20", "MS25"])
def test_the_checkers_one_builds_the_banked_bytes_with_the_factors_at_1(bt, ms4, ledger, tmp_path, new):
    chk = _check()
    f0, _ = _build(bt, tmp_path, "T1_s2rcB25", BANKED)
    argv = {"H05": HM[0.5], "H10": HM[1.0], "H20": HM[2.0], "MS25": MS25}[new]
    out = tmp_path / f"{new}_one.h5ad"
    before = bt.t1_scope.momentum_shift
    assert chk.main(["one", "t1", *RC, *argv, "--out", str(out)]) == 0
    assert out.read_bytes() == f0.read_bytes()
    assert bt.t1_scope.momentum_shift is before  # the forced call is put back


def test_the_checkers_one_fails_without_a_momentum_step_and_on_a_refusal(bt, ms4, ledger, tmp_path):
    chk = _check()
    assert chk.main(["one", "t1", *RC, *BANKED, "--out", str(tmp_path / "B25_one.h5ad")]) == 1  # no step ran
    with pytest.raises(SystemExit, match="named but not coded"):
        chk.main(["one", "t1", *RC, *HM[0.5], "--direction-momentum", "all", "--out", str(tmp_path / "x.h5ad")])
    for argv in ([], ["one"], ["one", "heart"], ["rows", "a"], ["record", "a"], ["record", "a", "0.5", "b"], ["x"]):
        with pytest.raises(SystemExit, match="usage"):
            chk.main(argv)


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
    f0, _ = _build(bt, tmp_path, "T1_s2rcB25", BANKED)
    hs = [_build(bt, tmp_path, NAMES[v], HM[v])[0] for v in HM_OPT]
    capsys.readouterr()
    assert chk.main(["rows", str(f0), *map(str, hs)]) == 0
    got = _lines(capsys)
    assert set(got.values()) == {"True"} and len(got) == 3 + 3 * 3
    blood = _labels(ms4[0]["before"]) == "Blood"
    bad_other = _write_altered(hs[0], tmp_path / "T1_s2rcB25H05x.h5ad", np.flatnonzero(~blood)[:1])
    assert chk.main(["rows", str(f0), str(bad_other)]) == 1
    got = _lines(capsys)
    assert got["T1_s2rcB25H05x_other_rows_are_the_banked_files"] == "False"
    assert chk.main(["rows", str(f0), str(f0)]) == 1  # the banked file as an H file: no Blood row moved
    assert _lines(capsys)["T1_s2rcB25_a_blood_row_moved"] == "False"
    bad_banked = _write_altered(f0, tmp_path / "T1_s2rcB25x.h5ad", np.flatnonzero(blood)[:1])
    assert chk.main(["rows", str(bad_banked), str(hs[0])]) == 1  # the draw's Blood rows not where the file holds them
    assert _lines(capsys)["banked_blood_rows_at_drawn_counts"] == "False"
    assert [p for p in (tmp_path / "tmp").rglob("*") if p.is_file()] == []  # the draws wrote nothing


def test_the_checkers_record_passes_on_the_h_records_and_fails_on_altered_ones(bt, ms4, ledger, tmp_path, capsys):
    chk = _check()
    hs = {v: _build(bt, tmp_path, NAMES[v], HM[v])[0] for v in HM_OPT}
    capsys.readouterr()
    pairs = [x for v, f in hs.items() for x in (str(f), f"{v:.1f}")]
    assert chk.main(["record", *pairs]) == 0
    got = _lines(capsys)
    assert set(got.values()) == {"True"} and len(got) == 3 * 6
    assert chk.main(["record", str(hs[0.5]), "1.0"]) == 1  # another damping
    assert _lines(capsys)["T1_s2rcB25H05_its_damping"] == "False"
    side = hs[0.5].with_suffix(".build.json")
    rec = json.loads(side.read_text())
    for path, value, key in ((("fit", "library_3_in_fit"), True, "library_3_out"),
                             (("fit", "within"), {"blood": ["Blood", "Macrophage"]}, "blood_cells_only"),
                             (("fit", "later_libraries"), ["1", "2", "3"], "fit_e9_5_to_e10_5_libraries_1_and_2"),
                             (("status",), "skipped: fewer than 20 cells at a stage", "blood_shifted"),
                             (("stream",), [SEED, 1], "its_own_stream")):
        bad = json.loads(json.dumps(rec))
        node = bad["t1_candidate"]["momentum"]
        for k in path[:-1]:
            node = node[k]
        node[path[-1]] = value
        f = tmp_path / "alt" / "T1_s2rcB25H05.h5ad"
        f.parent.mkdir(exist_ok=True)
        f.with_suffix(".build.json").write_text(json.dumps(bad))
        assert chk.main(["record", str(f), "0.5"]) == 1, key
        assert _lines(capsys)[f"T1_s2rcB25H05_{key}"] == "False", key


def test_the_checker_writes_and_deletes_nothing_and_reads_only_its_fields():
    src = (ROOT / "experiments" / "builds" / "t1_hm_check.py").read_text()
    body = src.split('"""', 2)[2]  # after the module docstring
    code = "\n".join(line.split("#")[0] for line in body.splitlines())
    for word in ("unlink", "rmtree", "os.remove", "rmdir", "write_text", "write_h5ad", "write_bytes", "n_obs)",
                 "print(len", "print(x", "upload", "submit"):
        assert word not in code, word
    assert code.count("print(") == 1 and "print(k, bool(v))" in code  # booleans only
    assert ".build.json" in code and code.count("read_text()") == 1  # the record's fields, read by code
