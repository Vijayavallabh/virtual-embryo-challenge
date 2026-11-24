"""`build_test.py t1` with T1 lines 7 and 6 of row R22 and row R37's pooled base (experiments/harness/TRIES.md; the
captain's orders of 2026-11-24, ruling rows R39 and R40, and his named exception of 2026-11-24; docs/PLAN.md decision
log 2026-11-24 (bb) and (z)), each on R13's command (t1 --g-mom --base lib3 --feng-groups recomputed):
  line 7  --g-mom-damping 1.0 --feng-fraction-mult 3.0 --captain-ruling R22      T1_s2rcCm30 (file C's, step at 3.0)
  line 6  --g-mom-damping 1.5|2.0 --feng-fraction-mult 2.5 --captain-ruling R22  T1_s2rcB25d15, T1_s2rcB25d20
  R37     --base pooled --captain-ruling R37                                    T1_s2rc_pooled (every E10.5 cell)

Held here on synthetic stages and directions only (the build_t1 fixture of tests/test_hidden_counts.py and the
recomputed setup of tests/test_feng_relabel.py), with temporary ledgers in place of the repository's (no test here
reads the repository's own ledger or register): (1) nothing that builds today changes: the existing commands give the
bytes and records pinned by the earlier lines' tests, and the fifth line's two files and T1_s2rcMSB those the builder
at main 0470b57 (before these options existed) gave on these fixtures, pinned below; (2) each new file differs from its
base file only as its line says; (3) each new set is taken only with its ruling and only while that ruling row's first
cell names it, and every other use is refused before anything is read; the multiplier's routing still refuses 2.5 and
3.0 with any other damping; (4) the builder at main 0470b57 (read with `git show`) and this one route every other
argument set of a grid exactly alike against synthetic ledgers; (5) experiments/builds/t1_pool_depth_record.py reads
E10.5 alone and prints ratios of its pools' median library sizes. No released file, outside-atlas file, board or score
is read."""

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

import feng_relabel as fr  # noqa: E402
import test_t1_candidates_build as tcb  # noqa: E402
import test_t1_feng_groups_none  # noqa: E402, F401  (imported before any fixture, as tests/test_feng_relabel.py does)
from conftest import git_env  # noqa: E402
from test_feng_relabel import _recomputed_env  # noqa: E402
from test_hidden_counts import E105, _t1_stage, t1_env  # noqa: E402, F401  (the synthetic build_t1 fixture)
from test_t1_candidates_build import RC, THREE, _build, _no_read, _rows_changed, _same_frame, _spy  # noqa: E402
from test_t1_damping_line_build import DAMP as LINE2  # noqa: E402
from test_t1_damping_line_build import LINE1_BEFORE  # noqa: E402
from test_t1_fraction2_build import BEFORE, R20, _record_sha  # noqa: E402
from test_t1_fraction_line_build import R20_BEFORE  # noqa: E402
from test_t1_fraction_line_build import R22 as LINE1  # noqa: E402
from test_t1_line3_build import LINE2_BEFORE, LINE3, ROW_R20  # noqa: E402
from test_t1_line5_build import FIRST_CELL as FIRST_CELL_LINE5  # noqa: E402
from test_t1_line5_build import LINE5  # noqa: E402
from test_t1_ms_build import LINE3_BEFORE, ROW_R23  # noqa: E402
from test_t1_ms_build import MS as MS_ARGS  # noqa: E402
from test_t1_msb_build import MS_BEFORE, MSB, ROW_R26, _three_libraries  # noqa: E402

from vec import harness  # noqa: E402

bt, rc3, rc4 = tcb.bt, tcb.rc3, tcb.rc4  # the candidates' fixtures (the builder module; the recomputed setups)

OPT_CM30 = "`--g-mom-damping 1.0 --feng-fraction-mult 3.0`"
OPT_B25D15 = "`--g-mom-damping 1.5 --feng-fraction-mult 2.5`"
OPT_B25D20 = "`--g-mom-damping 2.0 --feng-fraction-mult 2.5`"
OPT_POOLED = "`--base pooled` with `--feng-groups recomputed`"
LINE7 = ["--g-mom-damping", "1.0", "--feng-fraction-mult", "3.0", "--captain-ruling", "R22"]
LINE6 = {1.5: ["--g-mom-damping", "1.5", "--feng-fraction-mult", "2.5", "--captain-ruling", "R22"],
         2.0: ["--g-mom-damping", "2.0", "--feng-fraction-mult", "2.5", "--captain-ruling", "R22"]}
NAMES6 = {1.5: "T1_s2rcB25d15", 2.0: "T1_s2rcB25d20"}
C = BEFORE["C"][0]  # file C's options: --g-mom-damping 1.0 --feng-fraction-mult 1.5
POOLED = ["t1", "--g-mom", "--base", "pooled", "--feng-groups", "recomputed", "--captain-ruling", "R37"]
# Row R22's sentence for lines 6 and 7 and row R37's first cell, in the words the ledger's first cells use, inside
# synthetic rows of that form (not the ledger's own rows)
SIX_SEVEN = ("; by T1's lines 6 and 7 added in decision log 2026-11-24 (bb), on the captain's own order of 2026-11-24 "
             f"(ruling rows R39 and R40), the T1 builder is also to take {OPT_CM30}, {OPT_B25D15} and {OPT_B25D20}, "
             "each only with `--captain-ruling R22`")
FIRST_CELL = FIRST_CELL_LINE5[:-1] + SIX_SEVEN + ")"
ROW_R22 = f"| {FIRST_CELL} | T1 | the lines |"
ROW_R22_BEFORE = f"| {FIRST_CELL_LINE5} | T1 | the lines |"  # the row before lines 6 and 7: neither named
FIRST_CELL_R37 = ("R37 (the captain's named exception of 2026-11-24, \"Yes, T1 pooled base, my exception\": T1's "
                  "depth, closed after the 4 November exposure, may change for one file, R13's command on the pooled "
                  "base; decision log 2026-11-24 (z) and (bb); ruling, no run, not a try; the T1 builder is to take "
                  f"{OPT_POOLED} only with `--captain-ruling R37`)")
ROW_R37 = f"| {FIRST_CELL_R37} | T1 | the pooled base |"
# The fifth line's two files (on the shared four-group fixture) and T1_s2rcMSB (on the three-library fixture of
# tests/test_t1_msb_build.py), built by the builder at main 0470b57 (before these options existed), twice in separate
# processes with the same result: the file's sha256 and its record's (feng_hygiene.commit set to "own_repo"), the same
# kind of synthetic-fixture pin as those of the earlier lines.
LINE5_BEFORE = {1.5: ("8d7e768577a5acf79b794c09aaf8d00ec56d1d48595d17cf5c5303e46fc675a2",
                      "5499eaf563e221cd2177421a5dc7b8e1690b84e05db377266da7c42399e6e95f"),
                2.0: ("d09dccb103ecdc181092265bc4459a775c7f27e16a336efd9a68f425cee464b8",
                      "c9bfb32e1f3965529ae9ff35e3d2cebb787768a7b40ce774a93bcff2a2fa565c")}
MSB_BEFORE = ("b36ffda8f4a2f6617adb9566e81d1d852496cdd7cbb79dddc321408f8431eb3f",
              "3399635651da5796ba461bfdebe4e275c0ba520510a17a65c48b8e94e2fb37da")
NUMBERS_ALLOWED = {"0.5", "1", "1.5", "2", "2.5", "3", "13", "20", "22", "23", "26", "37"}  # values, lib3, rows' IDs
SEED = harness.BUILD_SEED_OFFSET + 1


def _sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def _ledger_text(r22: str | None = ROW_R22, r37: str | None = ROW_R37) -> str:
    return "# ledger\n" + "".join(f"{row}\n" for row in (ROW_R20, r22, ROW_R23, ROW_R26, r37) if row)


@pytest.fixture
def ledger(bt, tmp_path, monkeypatch):
    """A temporary ledger holding the ruling rows R20, R22 (with lines 6 and 7), R23, R26 and R37; the builder reads it
    in place of the repository's."""
    p = tmp_path / "TRIES.md"
    p.write_text(_ledger_text())
    monkeypatch.setattr(bt, "TRIES_FILE", p)
    return p


@pytest.fixture
def ms4(bt, t1_env, monkeypatch):  # noqa: F811
    _recomputed_env(bt, t1_env, monkeypatch)
    _three_libraries(t1_env)
    return _spy(bt, monkeypatch)


def _refused(bt, monkeypatch, tmp_path, argv: list[str], match: str) -> str:
    """The command is refused before anything is read (no T1 file, Feng file, board, union or relabelling record), and
    its message shows no number but the options' values, lib3's, the rows' IDs and what the command line said."""
    for obj, name in ((bt, "load_stage"), (bt, "load_boards"), (bt.hygiene, "load_union"),
                      (bt.feng_dirs, "feng_directions"), (fr, "committed_at_head")):
        monkeypatch.setattr(obj, name, _no_read)
    name = "T1_dry.h5ad" if "--board" in argv else "T1.h5ad"
    with pytest.raises(SystemExit, match=match) as e:
        bt.main([*argv, "--out", str(tmp_path / name)])
    msg = str(e.value)
    given = set(re.findall(r"\d+(?:\.\d+)?", " ".join(argv)))
    assert set(re.findall(r"\d+(?:\.\d+)?", msg)) <= NUMBERS_ALLOWED | given, msg
    return msg


def _parse(bt, argv: list[str]):
    return bt.build_parser().parse_args([*argv, "--out", "x.h5ad"])


# ------------------------------------------------------------------ the constants and the ruling check


def test_the_new_sets_have_their_own_constants_and_the_earlier_ones_are_unchanged(bt):
    assert bt.G_MOM_DAMPINGS == (0.5, 1.0) and bt.FENG_FRACTION_MULTS == (1.0, 1.5)
    assert bt.FENG_FRACTION_MULTS_R22 == (2.5, 3.0) and set(bt.FENG_FRACTION_MULTS_RULED) == {2.0, 2.5, 3.0}
    assert bt.T1_R22_LINE3_OPTIONS == {2.0: "`--g-mom-damping 1.0 --feng-fraction-mult 2.0`",
                                       2.5: "`--g-mom-damping 1.0 --feng-fraction-mult 2.5`"}
    assert bt.T1_R22_LINE5_OPTIONS == {1.5: "`--g-mom-damping 1.5 --feng-fraction-mult 2.0`",
                                       2.0: "`--g-mom-damping 2.0 --feng-fraction-mult 2.0`"}
    assert bt.T1_R22_LINE6_OPTIONS == {1.5: OPT_B25D15, 2.0: OPT_B25D20} and bt.FENG_FRACTION_MULT_LINE6 == 2.5
    assert bt.G_MOM_DAMPINGS_LINE6 == (1.5, 2.0) and bt.T1_LINE6_FILES == NAMES6
    assert bt.T1_R22_LINE7_OPTION == OPT_CM30 and (bt.G_MOM_DAMPING_LINE7, bt.FENG_FRACTION_MULT_LINE7) == (1.0, 3.0)
    assert bt.T1_LINE7_FILE == "T1_s2rcCm30"
    assert (bt.T1_R37_ROW, bt.T1_R37_KIND, bt.T1_R37_OPTION) == ("R37", "named exception", OPT_POOLED)
    assert bt.T1_POOLED_COMMAND == " ".join(POOLED) and bt.T1_POOLED_FILE == "T1_s2rc_pooled.h5ad"
    assert bt.T1_BASES == ("pooled", "lib3", "mid") and bt.T1_CANDIDATE_COMMAND == " ".join(["t1", *RC])
    a = _parse(bt, ["t1", *RC])
    assert (a.g_mom_damping, a.feng_fraction_mult, a.captain_ruling, a.base) == (None, None, None, "lib3")


def test_the_help_and_the_docstring_name_lines_6_and_7_and_row_r37(bt):
    import io
    from contextlib import redirect_stdout

    buf = io.StringIO()
    with redirect_stdout(buf), pytest.raises(SystemExit):
        bt.build_parser().parse_args(["t1", "--help"])
    text = " ".join(buf.getvalue().split())
    assert "sixth search line (T1 line 6)" in text and "seventh search line (T1 line 7)" in text
    assert "sixth T1 search line" in text and "seventh T1 search line" in text
    assert "row R37's one file" in text and "the captain's named exception for R13's command on the pooled base" in text
    assert "fifth search line" in text and "fifth T1 search line" in text  # the fifth line's, as before
    doc = " ".join(bt.__doc__.split())
    for word in ("T1_s2rcB25d15", "T1_s2rcB25d20", "T1_s2rcCm30", "T1_s2rc_pooled", "T1 line 6", "T1 line 7",
                 "--captain-ruling R37", "t1_pool_depth_record.py"):
        assert word in doc, word


def test_the_order_check_reads_each_new_option_from_its_rows_first_cell(bt, ledger):
    for opt in (OPT_CM30, OPT_B25D15, OPT_B25D20):
        assert bt._captain_order_names("R22", opt, ledger)
        assert bt._captain_order_names("R22", opt, ledger, kind="order")  # the default kind, spelled out
        assert not bt._captain_order_names("R37", opt, ledger, kind="named exception")
    assert bt._captain_order_names("R37", OPT_POOLED, ledger, kind="named exception")
    assert not bt._captain_order_names("R37", OPT_POOLED, ledger)  # row R37 is no "captain's order"
    for opt in (*bt.T1_R22_LINE3_OPTIONS.values(), *bt.T1_R22_LINE5_OPTIONS.values(), "`--feng-fraction-mult 2.5`",
                "`--feng-fraction-mult 3.0`", "`--g-mom-damping 1.5`", "`--g-mom-damping 2.0`"):
        assert bt._captain_order_names("R22", opt, ledger)  # the earlier lines' options, as before
    ledger.write_text(_ledger_text(r22=ROW_R22_BEFORE, r37=None))  # the ledger before these rulings
    for opt in (OPT_CM30, OPT_B25D15, OPT_B25D20):
        assert not bt._captain_order_names("R22", opt, ledger)
    assert not bt._captain_order_names("R37", OPT_POOLED, ledger, kind="named exception")
    assert "`--feng-fraction-mult 3.0`" in ROW_R22_BEFORE and "`--g-mom-damping 1.5`" in ROW_R22_BEFORE


@pytest.mark.parametrize("text, ok", [
    (ROW_R22, True),
    (f"| R22 (the captain's order of X; ruling, no run; {OPT_CM30} with `--captain-ruling R22`) | T1 |", True),
    ("| R22 (the captain's order of X; ruling, no run; `--g-mom-damping 1.0` `--feng-fraction-mult 3.0` "
     "`--captain-ruling R22`) |", False),  # the parts, not the combined option
    ("| R22 (the captain's order of X; ruling, no run; --g-mom-damping 1.0 --feng-fraction-mult 3.0 "
     "`--captain-ruling R22`) |", False),  # not in backticks
    (f"| R22 (the captain's order of X; a try; {OPT_CM30} `--captain-ruling R22`) |", False),
    (f"| R22 (awaits the captain's order; ruling, no run; {OPT_CM30} `--captain-ruling R22`) |", False),
    (f"| R22 (the captain's order of X; ruling, no run) | T1 | {OPT_CM30} `--captain-ruling R22` |", False),
    (f"| R22 (the captain's order of X; ruling, no run; {OPT_CM30}) |", False),
    (f"| R22 (the captain's order of X; ruling, no run; {OPT_CM30} `--captain-ruling R39`) |", False),
    (f"| R39 (the captain's order of X; ruling, no run; {OPT_CM30} `--captain-ruling R22`) |", False),  # R39's row
])
def test_the_order_check_takes_only_a_single_ruling_row_naming_line_7(bt, tmp_path, text, ok):
    p = tmp_path / "TRIES.md"
    p.write_text(f"# ledger\n{text}\n")
    assert bt._captain_order_names("R22", OPT_CM30, p) is ok
    p.write_text(f"# ledger\n{text}\n{text}\n")  # two rows with one ID: never
    assert not bt._captain_order_names("R22", OPT_CM30, p)


@pytest.mark.parametrize("text, ok15, ok20", [
    (ROW_R22, True, True),
    (f"| R22 (the captain's order of X; ruling, no run; {OPT_B25D15} with `--captain-ruling R22`) |", True, False),
    ("| R22 (the captain's order of X; ruling, no run; `--g-mom-damping 1.5` `--g-mom-damping 2.0` "
     "`--feng-fraction-mult 2.5` `--captain-ruling R22`) |", False, False),  # the parts, not the combined options
    (f"| R22 (the captain's order of X; ruling, no run) | T1 | {OPT_B25D15} {OPT_B25D20} `--captain-ruling R22` |",
     False, False),  # outside the first cell
    (f"| R40 (the captain's order of X; ruling, no run; {OPT_B25D15} {OPT_B25D20} `--captain-ruling R22`) |", False,
     False),  # R40's row: the builder reads row R22's
])
def test_the_order_check_takes_only_a_single_ruling_row_naming_each_line_6_pair(bt, tmp_path, text, ok15, ok20):
    p = tmp_path / "TRIES.md"
    p.write_text(f"# ledger\n{text}\n")
    assert bt._captain_order_names("R22", OPT_B25D15, p) is ok15
    assert bt._captain_order_names("R22", OPT_B25D20, p) is ok20


@pytest.mark.parametrize("text, ok", [
    (ROW_R37, True),
    (f"| R37 (the captain's named exception of X; ruling, no run; {OPT_POOLED} only with `--captain-ruling R37`) |",
     True),
    (f"| R37 (the captain's order of X; ruling, no run; {OPT_POOLED} only with `--captain-ruling R37`) |", False),
    (f"| R37 (awaits the captain's named exception; ruling, no run; {OPT_POOLED} `--captain-ruling R37`) |", False),
    (f"| R37 (the captain's named exception of X; a try; {OPT_POOLED} `--captain-ruling R37`) |", False),
    ("| R37 (the captain's named exception of X; ruling, no run; `--base pooled` and `--feng-groups recomputed` "
     "`--captain-ruling R37`) |", False),  # the parts, not the pair as the row writes it
    ("| R37 (the captain's named exception of X; ruling, no run; --base pooled with --feng-groups recomputed "
     "`--captain-ruling R37`) |", False),  # not in backticks
    (f"| R37 (the captain's named exception of X; ruling, no run; {OPT_POOLED}) |", False),
    (f"| R37 (the captain's named exception of X; ruling, no run; {OPT_POOLED} `--captain-ruling R22`) |", False),
    (f"| R37 (the captain's named exception of X; ruling, no run) | T1 | {OPT_POOLED} `--captain-ruling R37` |",
     False),  # outside the first cell
    (f"| R36 (the captain's named exception of X; ruling, no run; {OPT_POOLED} `--captain-ruling R37`) |", False),
])
def test_the_order_check_takes_only_a_single_named_exception_row_r37_naming_the_pair(bt, tmp_path, text, ok):
    p = tmp_path / "TRIES.md"
    p.write_text(f"# ledger\n{text}\n")
    assert bt._captain_order_names("R37", OPT_POOLED, p, kind="named exception") is ok
    p.write_text(f"# ledger\n{text}\n{text}\n")  # two rows with one ID: never
    assert not bt._captain_order_names("R37", OPT_POOLED, p, kind="named exception")


# ------------------------------------------------------------------ nothing existing changes


@pytest.mark.parametrize("name", list(BEFORE))
def test_the_existing_commands_give_the_same_bytes_and_records(bt, rc4, ledger, tmp_path, name):
    extra, file_sha, record_sha = BEFORE[name]
    f, rec = _build(bt, tmp_path, f"T1_rc{name}", extra)
    assert _sha(f) == file_sha and _record_sha(rec) == record_sha
    assert "base_exception" not in rec


def test_the_earlier_ruled_files_give_the_bytes_and_records_as_before(bt, rc4, ledger, tmp_path):
    f, rec = _build(bt, tmp_path, "T1_s2rcB2", R20)
    assert (_sha(f), _record_sha(rec)) == R20_BEFORE
    files = [("T1_s2rcB25", LINE1[2.5], LINE1_BEFORE[2.5]), ("T1_s2rcB30", LINE1[3.0], LINE1_BEFORE[3.0]),
             ("T1_s2rcA15", LINE2[1.5], LINE2_BEFORE[1.5]), ("T1_s2rcA20", LINE2[2.0], LINE2_BEFORE[2.0]),
             ("T1_s2rcCm20", LINE3[2.0], LINE3_BEFORE[2.0]), ("T1_s2rcCm25", LINE3[2.5], LINE3_BEFORE[2.5]),
             ("T1_s2rcB2d15", LINE5[1.5], LINE5_BEFORE[1.5]), ("T1_s2rcB2d20", LINE5[2.0], LINE5_BEFORE[2.0])]
    for name, argv, before in files:
        f, rec = _build(bt, tmp_path, name, argv)
        assert (_sha(f), _record_sha(rec)) == before, name
        text = json.dumps(rec)
        assert "sixth" not in text and "seventh" not in text and "line 6" not in text and "line 7" not in text, name
        assert "base_exception" not in rec, name


def test_ms_and_the_combined_ms_file_give_the_bytes_and_records_as_before(bt, ms4, ledger, tmp_path):
    f, rec = _build(bt, tmp_path, "T1_s2rc_MS", MS_ARGS)
    assert (_sha(f), _record_sha(rec)) == MS_BEFORE
    f, rec = _build(bt, tmp_path, "T1_s2rcMSB", MSB)
    assert (_sha(f), _record_sha(rec)) == MSB_BEFORE


# ------------------------------------------------------------------ line 7: T1_s2rcCm30


def test_line_7_is_file_c_with_the_step_at_3_0_and_its_record_names_line_7(bt, rc4, ledger, tmp_path):
    fc, rc = _build(bt, tmp_path, "T1_s2rcC", C)
    f7, r7 = _build(bt, tmp_path, "T1_s2rcCm30", LINE7)
    calendar = bt.feng_dirs.bracket_fraction(10.5, 12.5)
    assert rc4[0]["fraction"] == rc["fraction"] == pytest.approx(1.5 * calendar, abs=1e-12)
    assert rc4[-1]["fraction"] == r7["fraction"] == 3.0 * calendar == pytest.approx(1.2, abs=1e-12)
    assert rc4[0]["groups"] == rc4[-1]["groups"] and rc4[0]["directions"] == rc4[-1]["directions"]
    b0, b7 = rc4[0]["before"], rc4[-1]["before"]
    assert list(b0.obs_names) == list(b7.obs_names) and b0.obs.equals(b7.obs)
    cand = r7["t1_candidate"]
    assert (cand["g_mom_damping"], cand["feng_fraction_mult"], cand["captain_ruling"]) == (1.0, 3.0, "R22")
    assert cand["order"] == bt.T1_LINE7_ORDER_R22.format(row="R22", v=3.0)
    for word in ("row R22", "the captain's order of 2026-11-24", "ruling row R39", "T1 line 7",
                 "seventh T1 search line", "file C's command", "damping held at 1.0",
                 "direction at 3.0 x the calendar fraction", "nothing else changed"):
        assert word in cand["order"], word
    assert "sixth" not in cand["order"] and "fifth" not in cand["order"] and "third T1" not in cand["order"]
    assert "lifting the third line's cap by one half step" in cand["order"]
    assert cand["calendar_fraction"] == calendar and cand["applied_fraction"] == r7["fraction"]
    assert cand["command"] == bt.T1_CANDIDATE_COMMAND
    assert set(cand) == {"g_mom_damping", "feng_fraction_mult", "captain_ruling", "order", "calendar_fraction",
                         "applied_fraction", "command"}  # the ruled files' fields, no other
    # the record is C's but for the step's size and the candidate field; g_mom's part (damping 1.0, candidate A's
    # text) is C's exactly
    strip = lambda r: {k: v for k, v in r.items() if k not in ("t1_candidate", "fraction")}  # noqa: E731
    assert strip(r7) == strip(rc) and r7["g_mom"] == rc["g_mom"] and r7["g_mom"]["damping"] == 1.0
    assert "candidate A" in r7["g_mom"]["rule"] and "line" not in r7["g_mom"]["rule"]
    assert "base_exception" not in r7
    known = {*(b[1] for b in BEFORE.values()), R20_BEFORE[0], *(b[0] for b in LINE1_BEFORE.values()),
             *(b[0] for b in LINE2_BEFORE.values()), *(b[0] for b in LINE3_BEFORE.values()),
             *(b[0] for b in LINE5_BEFORE.values())}
    assert _sha(f7) not in known


def test_line_7_changes_only_cells_with_a_direction_and_builds_the_same_bytes_twice(bt, rc3, ledger, tmp_path):
    fc, _ = _build(bt, tmp_path, "T1_s2rcCm25", LINE3[2.5])
    fa, _ = _build(bt, tmp_path, "T1_s2rcCm30a", LINE7)
    fb, _ = _build(bt, tmp_path, "T1_s2rcCm30b", LINE7)
    assert fa.read_bytes() == fb.read_bytes()
    x25, x30 = _same_frame(fc, fa)  # same cells, order, obs and var as T1_s2rcCm25 (the line's base, one step back)
    changed = _rows_changed(x25, x30)
    assert changed.any()
    pop = rc3[1]["before"]
    labels = pop.obs["celltype"].astype(str).to_numpy()
    with_dir = np.array([rc3[1]["groups"].get(lab) in THREE for lab in labels])
    assert not with_dir.all() and not (changed & ~with_dir).any()  # Epicardium and the rest keep their counts
    fs, _ = _build(bt, tmp_path, "T1_s2rcCm30s", ["--g-mom-damping", "1", "--feng-fraction-mult", "3",
                                                  "--captain-ruling", "R22"])
    assert fs.read_bytes() == fa.read_bytes()  # the values spelled 1 and 3: the same file


# ------------------------------------------------------------------ line 6: T1_s2rcB25d15 and T1_s2rcB25d20


@pytest.mark.parametrize("v", [1.5, 2.0])
def test_line_6_is_the_2_5_file_with_only_g_moms_damping_changed(bt, rc4, ledger, tmp_path, v):
    fx, rx = _build(bt, tmp_path, "T1_s2rcB25", LINE1[2.5])
    fv, rv = _build(bt, tmp_path, NAMES6[v], LINE6[v])
    fa, ra = _build(bt, tmp_path, f"A{v}", LINE2[v])  # the second line's file at the same damping, multiplier 1.0
    f5, r5 = _build(bt, tmp_path, f"B2d{v}", LINE5[v])  # the fifth line's file at the same damping, multiplier 2.0
    calendar = bt.feng_dirs.bracket_fraction(10.5, 12.5)
    assert rc4[0]["fraction"] == rc4[1]["fraction"] == rx["fraction"] == rv["fraction"] == 2.5 * calendar
    assert rc4[0]["groups"] == rc4[1]["groups"] and rc4[0]["directions"] == rc4[1]["directions"]
    b0, bv = rc4[0]["before"], rc4[1]["before"]  # the population g_mom's hook gets: the 2.5 file's
    assert list(b0.obs_names) == list(bv.obs_names) and b0.obs.equals(bv.obs)
    assert np.array_equal(tcb._dense(b0.layers["counts"]), tcb._dense(bv.layers["counts"]))
    cand = rv["t1_candidate"]
    assert (cand["g_mom_damping"], cand["feng_fraction_mult"], cand["captain_ruling"]) == (v, 2.5, "R22")
    assert cand["order"] == bt.T1_LINE6_ORDER_R22.format(row="R22", v=v)
    for word in ("row R22", "the captain's order of 2026-11-24", "ruling row R40", "T1 line 6", "sixth T1 search line",
                 "2.5 file's command", f"damping at {v:.1f}", "nothing else changed"):
        assert word in cand["order"], word
    assert "fifth" not in cand["order"] and "second" not in cand["order"] and "seventh" not in cand["order"]
    assert cand["calendar_fraction"] == calendar and cand["applied_fraction"] == rv["fraction"]
    assert cand["command"] == bt.T1_CANDIDATE_COMMAND
    assert set(cand) == {"g_mom_damping", "feng_fraction_mult", "captain_ruling", "order", "calendar_fraction",
                         "applied_fraction", "command"}
    # the record is the 2.5 file's but for the candidate field and g_mom's damping and text
    strip = lambda r: {k: x for k, x in r.items() if k not in ("t1_candidate", "g_mom")}  # noqa: E731
    assert strip(rv) == strip(rx) and "base_exception" not in rv
    gx, gv = rx["g_mom"], rv["g_mom"]
    assert gx["damping"] == 0.5 and gv["damping"] == v
    assert {k: x for k, x in gv.items() if k not in ("damping", "rule")} == \
        {k: x for k, x in gx.items() if k not in ("damping", "rule")}
    assert (f"damped by {v:g} in place of B0_DAMPING (row R22's sixth T1 search line, on the 2.5 file's command, "
            "docs/PLAN.md decision log 2026-11-24 (bb), T1 line 6)") in gv["rule"]
    assert gv["rule"].replace(gv["rule"][gv["rule"].index("damped by"):gv["rule"].index(", zeroed")],
                              "damped by B0_DAMPING") == gx["rule"]
    assert "second T1 search line" in ra["g_mom"]["rule"] and "fifth T1 search line" in r5["g_mom"]["rule"]
    assert _sha(fa) == LINE2_BEFORE[v][0] and _sha(f5) == LINE5_BEFORE[v][0]  # their files as before
    known = {*(b[1] for b in BEFORE.values()), R20_BEFORE[0], *(b[0] for b in LINE1_BEFORE.values()),
             *(b[0] for b in LINE2_BEFORE.values()), *(b[0] for b in LINE3_BEFORE.values()),
             *(b[0] for b in LINE5_BEFORE.values())}
    assert _sha(fv) not in known


@pytest.mark.parametrize("v", [1.5, 2.0])
def test_line_6_changes_only_cells_g_mom_edits_and_builds_the_same_bytes_twice(bt, rc3, ledger, tmp_path, v):
    fx, _ = _build(bt, tmp_path, "T1_s2rcB25", LINE1[2.5])
    fa, _ = _build(bt, tmp_path, f"{NAMES6[v]}a", LINE6[v])
    fb, _ = _build(bt, tmp_path, f"{NAMES6[v]}b", LINE6[v])
    assert fa.read_bytes() == fb.read_bytes()
    xx, xv = _same_frame(fx, fa)  # same cells, order, obs and var as the 2.5 file
    changed = _rows_changed(xx, xv)
    assert changed.any()
    pop = rc3[1]["before"]
    labels = pop.obs["celltype"].astype(str).to_numpy()
    cats = pop.obs["category"].astype(str).to_numpy()
    held = bt._labels_without_direction(np.unique(labels), list(THREE))
    editable = np.isin(cats, bt.G_MOM_EDIT_CATEGORIES) & ~np.isin(labels, sorted(held))
    assert not editable.all() and not (changed & ~editable).any()  # only cells g_mom edits change
    b0 = rc3[0]["before"]  # the population g_mom's hook gets (after the direction step) is the 2.5 file's
    assert list(b0.obs_names) == list(pop.obs_names) and b0.obs.equals(pop.obs)
    assert np.array_equal(tcb._dense(b0.X), tcb._dense(pop.X))


def test_the_two_line_6_values_give_two_files_and_neither_is_an_earlier_file(bt, rc4, ledger, tmp_path):
    f15, _ = _build(bt, tmp_path, NAMES6[1.5], LINE6[1.5])
    f20, r20 = _build(bt, tmp_path, NAMES6[2.0], LINE6[2.0])
    f2, r2 = _build(bt, tmp_path, "T1_s2rcB25d2", ["--g-mom-damping", "2", "--feng-fraction-mult", "2.5",
                                                   "--captain-ruling", "R22"])
    f7, _ = _build(bt, tmp_path, "T1_s2rcCm30", LINE7)
    assert f2.read_bytes() == f20.read_bytes() and r2["t1_candidate"] == r20["t1_candidate"]
    assert len({f15.read_bytes(), f20.read_bytes(), f7.read_bytes()}) == 3


# ------------------------------------------------------------------ row R37: T1_s2rc_pooled


def _spy_pools(bt, monkeypatch) -> list[dict]:
    """Around build_submission: per call, the pool's cell names, the shares and the count it is given."""
    seen = []
    inner = bt.t1_scope.build_submission

    def capture(*a, **k):
        seen.append({"names": [str(x) for x in a[0].obs_names], "shares": dict(a[2]), "n": a[3]})
        return inner(*a, **k)

    monkeypatch.setattr(bt.t1_scope, "build_submission", capture)
    return seen


def test_the_pooled_file_is_r13s_command_on_every_e105_cell_at_the_stages_own_mix(bt, rc4, t1_env, ledger,  # noqa: F811
                                                                                  monkeypatch, tmp_path):
    pools = _spy_pools(bt, monkeypatch)
    _, r13 = _build(bt, tmp_path, "T1_s2rc", [])
    fp, rp = _build(bt, tmp_path, "T1_s2rc_pooled", ["--base", "pooled", "--captain-ruling", "R37"])
    ref = t1_env.ref
    names = [str(x) for x in ref.obs_names]
    assert pools[1]["names"] == names  # every E10.5 cell, both libraries
    assert pools[0]["names"] == [x for x in names if x.endswith("_3")] != names  # R13's: library 3's alone
    own = bt.harness.t1_chamber_shares(bt.harness._shallow(ref), "E10.5_RNA")
    assert pools[1]["shares"] == own != pools[0]["shares"]  # the stage's own mix
    assert pools[0]["n"] == pools[1]["n"]  # the count rule is R13's
    # the direction step and g_mom's fit and mask are R13's
    assert rc4[0]["fraction"] == rc4[1]["fraction"] == r13["fraction"] == rp["fraction"]
    assert rc4[0]["groups"] == rc4[1]["groups"] and rc4[0]["directions"] == rc4[1]["directions"]
    assert rp["base"] == "pooled" and r13["base"] == "lib3"
    assert rp["base_exception"] == {"base": "pooled", "captain_ruling": "R37", "order": bt.T1_POOLED_ORDER,
                                    "command": bt.T1_POOLED_COMMAND}
    for word in ("row R37", "named exception of 2026-11-24", "Yes, T1 pooled base, my exception", "T1's depth",
                 "R13's command on the pooled base", "nothing else changed"):
        assert word in rp["base_exception"]["order"], word
    assert "t1_candidate" not in rp and "mid_base" not in rp and "lib3_pair" not in rp and "released_bet" not in rp
    composition = ("base", "base_exception", "shares", "cm_share_of_heart", "blood_share_of_chamber", "realised",
                   "g_mom")
    assert {k: v for k, v in rp.items() if k not in composition} == \
        {k: v for k, v in r13.items() if k not in composition}
    assert {k: v for k, v in rp["g_mom"].items() if k != "median_depth_scale"} == \
        {k: v for k, v in r13["g_mom"].items() if k != "median_depth_scale"}
    assert rp["shares"] == own and rp["feng_groups"] == r13["feng_groups"]


def test_the_pooled_file_builds_the_same_bytes_twice_and_is_no_earlier_file(bt, rc4, ledger, tmp_path):
    fa, ra = _build(bt, tmp_path, "T1_s2rc_pooled_a", ["--base", "pooled", "--captain-ruling", "R37"])
    fb, rb = _build(bt, tmp_path, "T1_s2rc_pooled_b", ["--base", "pooled", "--captain-ruling", "R37"])
    assert fa.read_bytes() == fb.read_bytes() and _record_sha(ra) == _record_sha(rb)
    out = tmp_path / "T1_s2rc_pooled_c.h5ad"  # the base given by its default (pooled): the same command and file
    bt.main(["t1", "--g-mom", "--feng-groups", "recomputed", "--captain-ruling", "R37", "--out", str(out)])
    assert out.read_bytes() == fa.read_bytes()
    known = {*(b[1] for b in BEFORE.values()), R20_BEFORE[0], *(b[0] for b in LINE1_BEFORE.values()),
             *(b[0] for b in LINE2_BEFORE.values()), *(b[0] for b in LINE3_BEFORE.values()),
             *(b[0] for b in LINE5_BEFORE.values())}
    assert _sha(fa) not in known


# ------------------------------------------------------------------ refusals, before anything is read


LINE_SETS = {"Cm30": LINE7[:4], "B25d15": LINE6[1.5][:4], "B25d20": LINE6[2.0][:4]}  # each pair, without its ruling


@pytest.mark.parametrize("pair", list(LINE_SETS))
@pytest.mark.parametrize("extra, match", [
    ([], "is the captain's order of row R22 alone and needs --captain-ruling R22"),
    *[(["--captain-ruling", r], f"refused with --captain-ruling {r}: only row R22")
      for r in ("R20", "R21", "R23", "R26", "R37", "R39", "R40", "r22", "R220", "T63")],
    (["--captain-ruling", "R22", "--bet-strength", "0.5"], "with --released-bet and --bet-strength"),
    (["--captain-ruling", "R22", "--released-bet", "share_trend", "--bet-strength", "0.5"],
     "with --released-bet and --bet-strength"),
    (["--captain-ruling", "R22", "--seed", "1000002"], "recomputed last file only"),
    (["--captain-ruling", "R22", "--board", "T1:val"], "recomputed last file only"),  # a dry run
    (["--captain-ruling", "R22", "--no-hold-unbracketed"], "--feng-groups recomputed builds one file"),
    (["--captain-ruling", "R22", "--direction-momentum", "all"], "refused"),  # candidate MS's option never stacks
])
def test_every_other_use_of_a_line_pair_is_refused(bt, ledger, monkeypatch, tmp_path, pair, extra, match):
    _refused(bt, monkeypatch, tmp_path, ["t1", *RC, *LINE_SETS[pair], *extra], match)


@pytest.mark.parametrize("pair", list(LINE_SETS))
@pytest.mark.parametrize("argv", [
    ["t1"],  # the authors' groups
    ["t1", "--g-mom", "--base", "lib3"],  # the first upload's command
    ["t1", "--g-mom", "--base", "lib3", "--feng-groups", "none"],  # the none file
    ["t1", "--g-mom", "--base", "mid", "--cm-share", "0.60"],  # T1's second upload
])
def test_a_line_pair_is_refused_on_any_other_command(bt, ledger, monkeypatch, tmp_path, pair, argv):
    _refused(bt, monkeypatch, tmp_path, [*argv, *LINE_SETS[pair], "--captain-ruling", "R22"],
             "recomputed last file only")


@pytest.mark.parametrize("extra", [["--cm-share", "0.6"], ["--blood-share", "0.1"], ["--herd"], ["--no-feng-hygiene"],
                                   ["--lib3-pair"], ["--strong-null"], ["--base", "pooled"], ["--base", "mid"]])
def test_the_recomputed_commands_own_guards_hold_with_a_line_pair(bt, ledger, monkeypatch, tmp_path, extra):
    for pair in LINE_SETS.values():
        _refused(bt, monkeypatch, tmp_path, ["t1", *RC, *extra, *pair, "--captain-ruling", "R22"],
                 "--feng-groups recomputed builds one file: the first upload's command")


def test_a_line_pair_is_refused_while_row_r22_does_not_name_it(bt, ledger, monkeypatch, tmp_path):
    ledger.write_text(_ledger_text(r22=ROW_R22_BEFORE))  # the row before lines 6 and 7
    for pair, opt in (("Cm30", OPT_CM30), ("B25d15", OPT_B25D15), ("B25d20", OPT_B25D20)):
        msg = _refused(bt, monkeypatch, tmp_path, ["t1", *RC, *LINE_SETS[pair], "--captain-ruling", "R22"],
                       "is refused: experiments/harness/TRIES.md has no single row R22")
        assert opt in msg and "`--captain-ruling R22`" in msg
    ledger.write_text(_ledger_text(r22=None))  # no row R22
    _refused(bt, monkeypatch, tmp_path, ["t1", *RC, *LINE7], "no single row R22")
    ledger.write_text(_ledger_text(r22="| R22 (the captain's order of X; ruling, no run; `--g-mom-damping 1.0` "
                                       "`--feng-fraction-mult 3.0` `--captain-ruling R22`) | T1 |"))  # the parts
    _refused(bt, monkeypatch, tmp_path, ["t1", *RC, *LINE7], "no single row R22")
    ledger.write_text(_ledger_text(r22=f"| R22 (awaits the captain's order; {OPT_CM30} {OPT_B25D15} {OPT_B25D20} "
                                       "`--captain-ruling R22`) | T1 |"))
    for argv in (LINE7, LINE6[1.5], LINE6[2.0]):
        _refused(bt, monkeypatch, tmp_path, ["t1", *RC, *argv], "no single row R22")
    # the orders in rows R39 and R40 alone: the builder reads row R22's first cell, as the rulings say
    ledger.write_text(_ledger_text(r22=ROW_R22_BEFORE) + f"| R39 (the captain's order of X; ruling, no run; {OPT_CM30} "
                      "`--captain-ruling R22`) |\n" + f"| R40 (the captain's order of X; ruling, no run; {OPT_B25D15} "
                      f"{OPT_B25D20} `--captain-ruling R22`) |\n")
    for argv in (LINE7, LINE6[1.5], LINE6[2.0]):
        _refused(bt, monkeypatch, tmp_path, ["t1", *RC, *argv], "no single row R22")
    # a row naming one line 6 pair orders that pair alone; line 7 alone orders line 7 alone
    ledger.write_text(_ledger_text(r22=f"| R22 (the captain's order of X; ruling, no run; {OPT_B25D15} with "
                                       "`--captain-ruling R22`) | T1 |"))
    _refused(bt, monkeypatch, tmp_path, ["t1", *RC, *LINE6[2.0]], f"no single row R22 .* {re.escape(OPT_B25D20)}")
    _refused(bt, monkeypatch, tmp_path, ["t1", *RC, *LINE7], f"no single row R22 .* {re.escape(OPT_CM30)}")
    assert bt._t1_candidate(_parse(bt, ["t1", *RC, *LINE6[1.5]]))["captain_ruling"] == "R22"
    ledger.write_text(_ledger_text(r22=f"| R22 (the captain's order of X; ruling, no run; {OPT_CM30} with "
                                       "`--captain-ruling R22`) | T1 |"))
    for v in (1.5, 2.0):
        _refused(bt, monkeypatch, tmp_path, ["t1", *RC, *LINE6[v]], "no single row R22")
    assert bt._t1_candidate(_parse(bt, ["t1", *RC, *LINE7]))["feng_fraction_mult"] == 3.0


@pytest.mark.parametrize("ruling", [[], ["--captain-ruling", "R22"], ["--captain-ruling", "R20"]])
@pytest.mark.parametrize("damping, mult", [
    *[(d, "2.5") for d in ("0.75", "1.25", "1.75", "2.5", "3.0", "0", "-1.0")],
    *[(d, "3.0") for d in ("0.75", "1.25", "1.5", "2.0", "2.5", "3.0", "0", "-1.0")],
])
def test_the_multipliers_routing_still_refuses_2_5_and_3_0_with_any_other_damping(bt, ledger, monkeypatch, tmp_path,
                                                                                  ruling, damping, mult):
    msg = _refused(bt, monkeypatch, tmp_path, ["t1", *RC, "--g-mom-damping", damping, "--feng-fraction-mult", mult,
                                               *ruling], r"refused|takes|needs")
    if ruling == ["--captain-ruling", "R22"]:
        assert f"--feng-fraction-mult {float(mult):g} is refused with --g-mom-damping other than 0.5" in msg


def test_the_multipliers_values_are_taken_with_exactly_the_dampings_of_their_lines(bt, ledger):
    """With row R22's ruling, 2.5 is taken at damping 0.5 (or none: the first line), 1.0 (the third), 1.5 and 2.0
    (the sixth), and 3.0 at 0.5 (the first) and 1.0 (the seventh); the order texts say which line."""
    lines = {("0.5", "2.5"): "first search line", (None, "2.5"): "first search line",
             ("1.0", "2.5"): "third T1 search line", ("1.5", "2.5"): "sixth T1 search line",
             ("2.0", "2.5"): "sixth T1 search line", ("0.5", "3.0"): "first search line",
             (None, "3.0"): "first search line", ("1.0", "3.0"): "seventh T1 search line"}
    for (damping, mult), line in lines.items():
        argv = ["t1", *RC, "--feng-fraction-mult", mult, "--captain-ruling", "R22"]
        a = _parse(bt, argv if damping is None else [*argv, "--g-mom-damping", damping])
        cand = bt._t1_candidate(a)
        assert cand["captain_ruling"] == "R22" and cand["feng_fraction_mult"] == float(mult)
        assert line in cand["order"], (damping, mult)
        assert bt._bet_refusals(a) is None and bt._t1_pooled(a) is None


# the pooled base: refusals


def test_the_pooled_base_without_row_r37s_ruling_is_refused_as_before(bt, ledger, monkeypatch, tmp_path):
    for argv in (POOLED[:-2], ["t1", "--g-mom", "--feng-groups", "recomputed"], ["t1", *RC, "--base", "pooled"]):
        _refused(bt, monkeypatch, tmp_path, argv, "--feng-groups recomputed builds one file: the first upload's "
                                                  "command .* refused with --base pooled or mid")
    for r in ("R22", "R13", "R20", "R23", "R26", "R36", "R38", "R39", "R40", "r37", "R370", "T37"):
        _refused(bt, monkeypatch, tmp_path, [*POOLED[:-1], r], "--captain-ruling and --bet-strength go with "
                                                               "--released-bet only")


def test_the_pooled_base_is_refused_while_row_r37_does_not_name_it(bt, ledger, monkeypatch, tmp_path):
    rows = [None,  # no row R37
            f"| R37 (the captain's order of X; ruling, no run; {OPT_POOLED} only with `--captain-ruling R37`) |",
            f"| R37 (awaits the captain's named exception; ruling, no run; {OPT_POOLED} `--captain-ruling R37`) |",
            "| R37 (the captain's named exception of X; ruling, no run; `--base pooled` `--feng-groups recomputed` "
            "`--captain-ruling R37`) |",
            f"| R37 (the captain's named exception of X; ruling, no run; {OPT_POOLED}) |",
            f"| R37 (the captain's named exception of X; a try; {OPT_POOLED} `--captain-ruling R37`) |",
            f"| R37 (the captain's named exception of X; ruling, no run) | T1 | {OPT_POOLED} `--captain-ruling R37` |",
            f"{ROW_R37}\n{ROW_R37}"]  # two rows with one ID
    for row in rows:
        ledger.write_text(_ledger_text(r37=row))
        msg = _refused(bt, monkeypatch, tmp_path, POOLED, "--base pooled with --feng-groups recomputed is refused: "
                                                          "experiments/harness/TRIES.md has no single row R37")
        assert "captain's named exception" in msg and OPT_POOLED in msg


@pytest.mark.parametrize("extra, match", [
    (["--g-mom-damping", "0.5"], "is refused with --g-mom-damping: row R37's file is R13's command on the pooled"),
    (["--g-mom-damping", "1.0"], "is refused with --g-mom-damping: row R37's file"),
    (["--feng-fraction-mult", "1.0"], "is refused with --feng-fraction-mult: row R37's file"),
    (["--feng-fraction-mult", "1.5"], "is refused with --feng-fraction-mult: row R37's file"),
    (["--g-mom-damping", "1.0", "--feng-fraction-mult", "1.5"], "with --g-mom-damping, --feng-fraction-mult: row R37"),
    (["--bet-strength", "0.5"], "is refused with --bet-strength: row R37's file"),
    (["--released-bet", "share_trend", "--bet-strength", "0.5"], "with --released-bet, --bet-strength: row R37's"),
    (["--direction-momentum", "all"], "refused with --captain-ruling R37: only row R23"),
    (["--g-mom-damping", "1.5"], "refused with --captain-ruling R37: only row R22"),
    (["--feng-fraction-mult", "2.0"], "refused with --captain-ruling R37: only row R20"),
    (["--feng-fraction-mult", "3.0"], "refused with --captain-ruling R37: only row R22"),
    (LINE7[:4], "refused with --captain-ruling R37: only row R22"),
    (LINE6[1.5][:4], "refused with --captain-ruling R37: only row R22"),
    (["--seed", "1000002"], "is row R37's one file: it needs the command"),
    (["--board", "T1:val"], "is row R37's one file: it needs the command"),  # a dry run
    *[(x, "--base pooled with --feng-groups recomputed builds row R37's one file") for x in (
        ["--cm-share", "0.6"], ["--blood-share", "0.1"], ["--lib3-pair"], ["--herd"], ["--no-hold-unbracketed"],
        ["--no-feng-hygiene"], ["--strong-null"], ["--feng-none-ruling", "x"])],
])
def test_every_other_use_of_the_pooled_base_with_r37_is_refused(bt, ledger, monkeypatch, tmp_path, extra, match):
    _refused(bt, monkeypatch, tmp_path, [*POOLED, *extra], match)


def test_the_pooled_base_with_r37_is_refused_without_g_mom_and_on_any_other_command(bt, ledger, monkeypatch,
                                                                                    tmp_path):
    _refused(bt, monkeypatch, tmp_path, [p for p in POOLED if p != "--g-mom"],
             "--base pooled with --feng-groups recomputed builds row R37's one file")
    for argv in (["t1", "--g-mom", "--base", "pooled"], ["t1", "--g-mom", "--base", "pooled", "--feng-groups", "none"],
                 ["t1", *RC], ["t1", "--g-mom", "--base", "mid", "--feng-groups", "recomputed"],
                 ["t1", "--g-mom", "--base", "mid", "--cm-share", "0.60"], ["t1"]):
        _refused(bt, monkeypatch, tmp_path, [*argv, "--captain-ruling", "R37"],
                 "--captain-ruling and --bet-strength go with --released-bet only")


def test_the_bet_refusal_lets_only_the_ruled_commands_through(bt, ledger):
    def ns(**k):
        base = dict(released_bet=None, bet_strength=None, captain_ruling="R22", feng_fraction_mult=3.0,
                    feng_groups="recomputed", seed=SEED, board=None, g_mom_damping=1.0, base="lib3",
                    direction_momentum=None)
        return SimpleNamespace(**{**base, **k})

    pooled = dict(captain_ruling="R37", base="pooled", g_mom_damping=None, feng_fraction_mult=None)
    goods = (ns(), ns(g_mom_damping=1.5, feng_fraction_mult=2.5), ns(g_mom_damping=2.0, feng_fraction_mult=2.5),
             ns(**pooled))
    for good in (*goods, ns(g_mom_damping=None), ns(g_mom_damping=1.5, feng_fraction_mult=2.0)):
        assert bt._bet_refusals(good) is None
    assert bt._ruled_line7_ordered(goods[0]) and bt._pooled_ordered(goods[3])
    assert bt._ruled_line6_ordered(goods[1]) and bt._ruled_line6_ordered(goods[2])
    refused = (ns(captain_ruling="R39"), ns(captain_ruling="R37"), ns(bet_strength=0.5),
               ns(**pooled, feng_groups="none"), ns(**{**pooled, "captain_ruling": "R22"}),
               ns(**{**pooled, "base": "mid"}), ns(**{**pooled, "g_mom_damping": 0.5}),
               ns(**pooled, direction_momentum="all"))
    for bad in (*refused, ns(g_mom_damping=1.25, feng_fraction_mult=2.5),
                ns(**pooled, released_bet="share_trend", bet_strength=0.5)):
        assert not (bt._ruled_line7_ordered(bad) or bt._ruled_line6_ordered(bad) or bt._pooled_ordered(bad))
    for bad in refused:
        with pytest.raises(SystemExit, match="with --released-bet only"):
            bt._bet_refusals(bad)
    ledger.write_text(_ledger_text(r22=ROW_R22_BEFORE, r37=None))  # without the orders the new shortcuts close
    for good in goods:
        assert not (bt._ruled_line7_ordered(good) or bt._ruled_line6_ordered(good) or bt._pooled_ordered(good))
    # (a line pair under the bare multiplier's shortcut is _t1_candidate's to refuse, and it runs first)
    with pytest.raises(SystemExit, match="with --released-bet only"):
        bt._bet_refusals(goods[3])


def test_each_new_set_is_taken_with_its_ruling_and_a_ledger_naming_it(bt, ledger, tmp_path, monkeypatch):
    monkeypatch.setattr(fr, "committed_at_head", lambda p: json.dumps({"route": "recomputed"}).encode())
    for argv, want in ((LINE7, (1.0, 3.0)), (LINE6[1.5], (1.5, 2.5)), (LINE6[2.0], (2.0, 2.5))):
        a = _parse(bt, ["t1", *RC, *argv])
        cand = bt._t1_candidate(a)
        assert (cand["g_mom_damping"], cand["feng_fraction_mult"], cand["captain_ruling"]) == (*want, "R22")
        assert bt._t1_ms(a) is None and bt._t1_pooled(a) is None and bt._bet_refusals(a) is None
        assert bt._feng_recomputed_refusals(a) is None
    a = _parse(bt, POOLED)
    assert bt._t1_ms(a) is None and bt._t1_candidate(a) is None and bt._bet_refusals(a) is None
    assert bt._t1_pooled(a) == {"base": "pooled", "captain_ruling": "R37", "order": bt.T1_POOLED_ORDER,
                                "command": bt.T1_POOLED_COMMAND}
    assert bt._feng_recomputed_refusals(a) is None  # the ruled exception: the pooled base passes R13's guards


# ------------------------------------------------------------------ the old builder against this one


OLD_COMMIT = "0470b57"  # main before these options existed
OLD_SHA = "6ecd24ede9cad0b5e02dcdb3a7e758fa5980db6c0eabeb4b5e4775b4c30f4586"  # its experiments/builds/build_test.py


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
    spec = importlib.util.spec_from_file_location("build_test_before_l67pb", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


CMDS = [[], ["--g-mom", "--base", "lib3"], RC, ["--g-mom", "--base", "pooled", "--feng-groups", "recomputed"],
        ["--base", "pooled", "--feng-groups", "recomputed"], ["--g-mom", "--feng-groups", "recomputed"],
        ["--g-mom", "--base", "mid", "--feng-groups", "recomputed"],
        ["--g-mom", "--base", "lib3", "--feng-groups", "none"],
        ["--g-mom", "--base", "pooled", "--feng-groups", "none"], ["--g-mom", "--base", "pooled"],
        ["--g-mom", "--base", "mid", "--cm-share", "0.60"]]
DAMPINGS = [None, "0.5", "1.0", "1.5", "2.0", "0.75", "2.5"]
MULTS = [None, "1.0", "1.5", "2.0", "2.5", "3.0", "2.25"]
RULINGS = [None, "R20", "R22", "R23", "R26", "R37", "R39", "R40"]
EXTRAS = [["--direction-momentum", "all"], ["--seed", "1000002"], ["--board", "T1:val"],
          ["--released-bet", "share_trend", "--bet-strength", "0.5"], ["--bet-strength", "0.5"],
          ["--no-hold-unbracketed"], ["--cm-share", "0.6"], ["--herd"], ["--no-feng-hygiene"], ["--strong-null"],
          ["--lib3-pair"], ["--feng-none-ruling", "x"]]
EXTRA_PAIRS = [(None, None), ("1.0", "3.0"), ("1.5", "2.5"), ("2.0", "2.5"), ("0.5", "2.5"), ("1.0", "2.5"),
               ("1.5", "2.0"), (None, "3.0"), ("1.0", None), (None, "1.5")]
NEW_PAIRS = {(1.0, 3.0): "line7", (1.5, 2.5): "line6", (2.0, 2.5): "line6"}
# the grid: every command, damping, multiplier and ruling; then R13's command and the pooled pair with each extra
# option, on the pairs that matter and three rulings. Each set is compared under each of the three ledgers below.
GRID_SETS = 4312 + 720
SETS_COMPARED = 3 * GRID_SETS


def _grid() -> list[list[str]]:
    def argv(cmd, damping, mult, ruling, extra=()):
        out = ["t1", *cmd]
        out += [] if damping is None else ["--g-mom-damping", damping]
        out += [] if mult is None else ["--feng-fraction-mult", mult]
        out += [] if ruling is None else ["--captain-ruling", ruling]
        return [*out, *extra]

    sets = [argv(*x) for x in itertools.product(CMDS, DAMPINGS, MULTS, RULINGS)]
    for cmd, (damping, mult), ruling, extra in itertools.product(CMDS[2:4], EXTRA_PAIRS, [None, "R22", "R37"], EXTRAS):
        sets.append(argv(cmd, damping, mult, ruling, extra))
    return sets


def _kind(a) -> str:
    """"new" for exactly one of the four new sets (R13's command with a line 6 or 7 pair and --captain-ruling R22, or
    R13's command on the pooled base with --captain-ruling R37, nothing else); "touch" for any other set holding a new
    pair or the pooled pair with R37; "other" for every other set."""
    pair = (a.g_mom_damping, a.feng_fraction_mult)
    pooled = a.base == "pooled" and a.feng_groups == "recomputed" and a.captain_ruling == "R37"
    if pair not in NEW_PAIRS and not pooled:
        return "other"
    plain = (a.g_mom and a.feng_groups == "recomputed" and a.seed == SEED and a.board is None
             and a.released_bet is None and a.bet_strength is None and a.direction_momentum is None
             and a.hold_unbracketed and a.feng_hygiene is None and not a.strong_null and a.cm_share is None
             and a.blood_share is None and not a.lib3_pair and not a.herd and a.feng_none_ruling is None)
    if plain and pair in NEW_PAIRS and a.base == "lib3" and a.captain_ruling == "R22":
        return "new"
    if plain and pooled and pair == (None, None):
        return "new"
    return "touch"


def test_the_old_builder_and_this_one_route_every_other_argument_set_alike(bt, tmp_path, monkeypatch):
    """Side by side on the grid, under three synthetic ledgers (the rows after these rulings; the rows before them;
    line 7 alone with row R37 written as an order, not a named exception): every set that holds none of the new pairs
    and not the pooled pair with R37 is refused with the same message, or reaches the same first read with the same
    candidate and MS values, in both builders; every set that holds one is refused by the old builder, and by this one
    unless it is exactly a new set whose ruling row names it, which this one takes with the line's values. Nothing is
    built: every read is stopped where it would start."""
    old = _old_builder(tmp_path)
    assert not hasattr(old, "_t1_pooled") and hasattr(bt, "_t1_pooled")
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
            if hasattr(mod, name):
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

    want = {"line7": ("taken", "union", None, {"g_mom_damping": 1.0, "feng_fraction_mult": 3.0,
                                               "captain_ruling": "R22",
                                               "order": bt.T1_LINE7_ORDER_R22.format(row="R22", v=3.0)}, None),
            "pooled": ("taken", "union", None, None, {"base": "pooled", "captain_ruling": "R37",
                                                      "order": bt.T1_POOLED_ORDER,
                                                      "command": bt.T1_POOLED_COMMAND})}
    for v in (1.5, 2.0):
        want[f"line6_{v}"] = ("taken", "union", None, {"g_mom_damping": v, "feng_fraction_mult": 2.5,
                                                       "captain_ruling": "R22",
                                                       "order": bt.T1_LINE6_ORDER_R22.format(row="R22", v=v)}, None)
    partial_r22 = f"| {FIRST_CELL_LINE5[:-1]}; and {OPT_CM30} only with `--captain-ruling R22`) | T1 | the lines |"
    as_order = f"| R37 (the captain's order of X; ruling, no run; {OPT_POOLED} only with `--captain-ruling R37`) |"
    ledgers = {"after": (_ledger_text(), {"line7", "line6", "pooled"}),
               "before": (_ledger_text(r22=ROW_R22_BEFORE, r37=None), set()),
               "partial": (_ledger_text(r22=partial_r22, r37=as_order), {"line7"})}
    grid = _grid()
    assert len({tuple(x) for x in grid}) == len(grid) == GRID_SETS
    count = {"other": 0, "touch": 0, "new": 0, "taken_alike": 0, "refused_alike": 0, "taken_new": 0}
    for text, named in ledgers.values():
        ledger.write_text(text)
        for argv in grid:
            a = _parse(bt, argv)
            kind = _kind(a)
            o, n = route(old, argv), route(bt, argv)
            count[kind] += 1
            if kind == "other":
                assert n == o, argv
                count["taken_alike" if o[0] == "taken" else "refused_alike"] += 1
                continue
            assert o[0] == "refused", argv  # the old builder took none of these
            line = "pooled" if a.base == "pooled" else NEW_PAIRS.get((a.g_mom_damping, a.feng_fraction_mult))
            if kind == "new" and line in named:
                key = line if line != "line6" else f"line6_{a.g_mom_damping}"
                assert n == want[key], argv
                count["taken_new"] += 1
            else:
                assert n[0] == "refused", argv
    assert sum(count[k] for k in ("other", "touch", "new")) == SETS_COMPARED
    assert count["taken_alike"] and count["refused_alike"] and count["touch"]
    # the exact new sets in the grid: each line pair with R22 on R13's command, the pooled pair with R37 written with
    # --base pooled or with the base left at its default; taken under the ledger after the rulings (4 + 1 + ... sets)
    assert count["new"] == 3 * 5 and count["taken_new"] == 5 + 1


# ------------------------------------------------------------------ the pooled base's depth record (f)(3)


def _depth_stage():
    """A synthetic E10.5 with released labels (chamber and other categories), library suffixes _1, _2 and _3, and
    library 3's cells twice as deep (counts doubled), so that the pools' medians differ."""
    from scipy import sparse

    ref = _t1_stage(E105, 2)
    libs = np.array(["1", "2", "3"])[np.arange(ref.n_obs) % 3]
    ref.obs_names = [f"c{i}_{lib}" for i, lib in enumerate(libs)]
    counts = sparse.csr_matrix(ref.layers["counts"]).astype(np.int64).toarray()
    counts[libs == "3"] *= 2
    ref.layers["counts"] = sparse.csr_matrix(counts)
    return ref, libs, counts


def test_the_pool_depth_record_reads_e105_alone_and_prints_ratios_of_median_library_sizes(monkeypatch, capsys):
    import t1_pool_depth_record as pdr

    from vec import boards, t1_scope

    ref, libs, counts = _depth_stage()
    loads = []

    def load(path, **kw):
        loads.append((Path(path), kw))
        return ref.copy()

    monkeypatch.setattr(pdr, "load_stage", load)
    monkeypatch.setattr(boards, "load_boards", _no_read)
    monkeypatch.setattr("anndata.read_h5ad", _no_read)  # no built file or other file is opened
    pdr.main([])
    assert loads == [(pdr.E105, {"lineage": False})] and pdr.E105.name == "E10.5_RNA.h5ad"
    out = json.loads(capsys.readouterr().out)
    spec = t1_scope.load_spec()
    cats = np.asarray(t1_scope.classify(ref.copy(), source="E10.5_RNA", spec=spec, key_added=None)).astype(str)
    chamber = np.isin(cats, spec.in_scope("chamber"))
    assert chamber.any() and not chamber.all()  # the fixture holds cells outside the chamber categories
    depth = counts.sum(axis=1)
    whole = np.median(depth[chamber])
    lib3 = np.median(depth[chamber & (libs == "3")]) / whole
    pools = out["pools"]
    assert set(pools) == {"pooled", "lib3"}
    assert pools["pooled"]["median_library_over_whole_chamber_pool"] == 1.0  # the pooled base's pool is the whole
    assert pools["lib3"]["median_library_over_whole_chamber_pool"] == pytest.approx(lib3, rel=1e-12)
    assert lib3 > 1.2 and out["pooled_over_lib3"] == pytest.approx(1.0 / lib3, rel=1e-12)
    assert set(out) == {"record", "source", "statistic", "denominator", "pools", "pooled_over_lib3"}
    numbers = [v for p in pools.values() for k, v in p.items() if not isinstance(v, str)] + [out["pooled_over_lib3"]]
    assert len(numbers) == 3 and all(isinstance(x, float) for x in numbers)  # ratios only, no cell number
    assert "(f)(3)" in out["record"] and "no board, test count, built file or build record is read" in out["source"]


def test_the_pool_depth_record_takes_no_argument_and_refuses_names_without_libraries(monkeypatch):
    import t1_pool_depth_record as pdr

    with pytest.raises(SystemExit, match="takes no arguments"):
        pdr.main(["--board", "T1:test"])
    ref, _, _ = _depth_stage()
    ref.obs_names = [f"c{i}" for i in range(ref.n_obs)]
    monkeypatch.setattr(pdr, "load_stage", lambda *a, **k: ref.copy())
    with pytest.raises(SystemExit, match="do not separate library 3"):
        pdr.main([])
    src = Path(pdr.__file__).read_text()
    assert "load_boards" not in src and "build.json" not in src and "read_h5ad" not in src
