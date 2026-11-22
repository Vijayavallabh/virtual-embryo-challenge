"""`build_test.py t1 --g-mom-damping 1.0 --feng-fraction-mult 2.0|2.5 --captain-ruling R22`: the third T1 search line
under row R22 of experiments/harness/TRIES.md (the captain's order of 2026-11-22; the line added in docs/PLAN.md
decision log 2026-11-22): file C's command (row R13's recomputed-groups file with g_mom's damping at 1.0 and the
outside atlas's direction step at 1.5 x the calendar fraction) with the step at 2.0 or 2.5 x the calendar fraction,
the damping held at 1.0 and nothing else changed.

Held here, on synthetic stages and directions only (the build_t1 fixture of tests/test_hidden_counts.py and the
recomputed setup of tests/test_feng_relabel.py), with a temporary ledger in place of the repository's: (1) nothing that
builds today changes: the default command and candidates A, B and C give the bytes and records pinned in
tests/test_t1_fraction2_build.py, row R20's file those pinned in tests/test_t1_fraction_line_build.py, the first line's
2.5 and 3.0 files those pinned in tests/test_t1_damping_line_build.py, and the second line's 1.5 and 2.0 files the
bytes and records the builder at main a6648c1 (before this option existed) gave on this fixture, pinned below (the same
kind of synthetic-fixture pin as those); (2) each new file differs from C only through the direction step's size; (3)
each pair is taken only with --captain-ruling R22 while row R22's first cell is the captain's order naming the combined
option, and every other use is refused before anything is read. No released file, outside-atlas file, board or score is
read."""

from __future__ import annotations

import hashlib
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
for _p in ("experiments/external", "experiments/t1_m1", "experiments/t1_direction", "experiments/t1_hygiene",
           "experiments/builds"):
    sys.path.insert(0, str(ROOT / _p))

import test_t1_candidates_build as tcb  # noqa: E402
import test_t1_feng_groups_none  # noqa: E402, F401  (imported before any fixture, as tests/test_feng_relabel.py does)
from test_hidden_counts import t1_env  # noqa: E402, F401  (the synthetic build_t1 fixture)
from test_t1_candidates_build import RC, THREE, _build, _rows_changed, _same_frame  # noqa: E402
from test_t1_damping_line_build import DAMP as LINE2  # noqa: E402
from test_t1_damping_line_build import LINE1_BEFORE  # noqa: E402
from test_t1_fraction2_build import BEFORE, R20, _record_sha, _refused  # noqa: E402
from test_t1_fraction_line_build import R20_BEFORE  # noqa: E402
from test_t1_fraction_line_build import R22 as LINE1  # noqa: E402

from vec import harness  # noqa: E402

bt, rc3, rc4 = tcb.bt, tcb.rc3, tcb.rc4  # the candidates' fixtures (the builder module; the recomputed setups)

C = BEFORE["C"][0]  # file C's options: --g-mom-damping 1.0 --feng-fraction-mult 1.5
LINE3 = {2.0: ["--g-mom-damping", "1.0", "--feng-fraction-mult", "2.0", "--captain-ruling", "R22"],
         2.5: ["--g-mom-damping", "1.0", "--feng-fraction-mult", "2.5", "--captain-ruling", "R22"]}
NAMES = {2.0: "T1_s2rcCm20", 2.5: "T1_s2rcCm25"}
R22 = ["--captain-ruling", "R22"]
# The second line's two files on this fixture, built by the builder at main a6648c1 (before the third line existed),
# twice in separate processes with the same result: the file's sha256 and its record's (feng_hygiene.commit set to
# "own_repo").
LINE2_BEFORE = {1.5: ("ccb2ed1e8638ce45b666ab6b2abffc1fa705d6eb864d4897d946eca105862558",
                      "dd1be3c5095e10a8933e7acb699e56a16e50418a6fd4155391e4c2971ff60388"),
                2.0: ("6f86bb0f36af012188d546f6f7288e51bd4786d30a71bac1011d197fa56b1fa6",
                      "0abb21f96390b763b8cad6dbbf09bde31a0a581c3de438715e039f2a994702cf")}
# Row R22's first cell, verbatim as decision log 2026-11-22 (nn) left it, inside a synthetic row of that form
FIRST_CELL = (
    "R22 (the captain's order of 2026-11-22: a score-guided search of round steps, decision log 2026-11-22 (j); the "
    "team's own bar on setting a number from a score is lifted for such searches and for nothing else; ruling, no run, "
    "not a try; the builder is to take `--feng-fraction-mult 2.5` and `--feng-fraction-mult 3.0` only with "
    "`--captain-ruling R22`; by the lines added in decision log 2026-11-22 the T1 builder is also to take "
    "`--g-mom-damping 1.5` and `--g-mom-damping 2.0`, each only with `--captain-ruling R22`; by the T3 line added in "
    "decision log 2026-11-22 (x) the T3 outside-form builder is to take `--table-exponent 4.0` on row K23's file for "
    "the barred s2, only with `--captain-ruling R22`; by the third T1 search line added in decision log 2026-11-22 "
    "(nn) the T1 builder is also to take `--g-mom-damping 1.0 --feng-fraction-mult 2.0` and "
    "`--g-mom-damping 1.0 --feng-fraction-mult 2.5`, each only with `--captain-ruling R22`)")
THIRD = ("; by the third T1 search line added in decision log 2026-11-22 (nn) the T1 builder is also to take "
         "`--g-mom-damping 1.0 --feng-fraction-mult 2.0` and `--g-mom-damping 1.0 --feng-fraction-mult 2.5`, each only "
         "with `--captain-ruling R22`")
ROW_R20 = ("| R20 (the captain's order of a day; ruling, no run, not a try; the builder is to take "
           "`--feng-fraction-mult 2.0` only with `--captain-ruling R20`) | T1 | a file |")
ROW_R22 = f"| {FIRST_CELL} | T1 | the lines |"
ROW_R22_BEFORE = f"| {FIRST_CELL.replace(THIRD, '')} | T1 | the lines |"  # the row before the third line was added


def _sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


@pytest.fixture
def ledger(bt, tmp_path, monkeypatch):
    """A temporary ledger holding the two ruling rows; the builder reads it in place of the repository's."""
    p = tmp_path / "TRIES.md"
    p.write_text(f"# ledger\n{ROW_R20}\n{ROW_R22}\n")
    monkeypatch.setattr(bt, "TRIES_FILE", p)
    return p


def _strip(r: dict) -> dict:
    """A record without what the step's size sets: the candidate field and the applied fraction."""
    return {k: v for k, v in r.items() if k not in ("t1_candidate", "fraction")}


# ------------------------------------------------------------------ the constants and the ruling check


def test_the_line_has_its_own_constants_and_the_earlier_ones_are_unchanged(bt):
    assert bt.G_MOM_DAMPINGS == (0.5, 1.0) and bt.FENG_FRACTION_MULTS == (1.0, 1.5)
    assert bt.G_MOM_DAMPING_LINE3 == 1.0 and bt.FENG_FRACTION_MULTS_LINE3 == (2.0, 2.5)
    assert bt.T1_R22_LINE3_OPTIONS == {2.0: "`--g-mom-damping 1.0 --feng-fraction-mult 2.0`",
                                       2.5: "`--g-mom-damping 1.0 --feng-fraction-mult 2.5`"}
    assert bt.FENG_FRACTION_MULTS_R22 == (2.5, 3.0) and set(bt.FENG_FRACTION_MULTS_RULED) == {2.0, 2.5, 3.0}
    assert bt.FENG_FRACTION_MULTS_RULED[2.0][0] == "R20" and bt.FENG_FRACTION_MULTS_RULED[2.5][0] == "R22"
    assert bt.G_MOM_DAMPINGS_R22 == (1.5, 2.0) and bt.T1_R22_ROW == "R22"
    a = bt.build_parser().parse_args(["t1", *RC, "--out", "x.h5ad"])
    assert a.g_mom_damping is None and a.feng_fraction_mult is None and a.captain_ruling is None


def test_the_order_check_reads_each_combined_option_from_row_r22s_first_cell(bt, ledger):
    for v in (2.0, 2.5):
        assert bt._captain_order_names("R22", bt.T1_R22_LINE3_OPTIONS[v], ledger)
        assert not bt._captain_order_names("R20", bt.T1_R22_LINE3_OPTIONS[v], ledger)
    for v in (2.5, 3.0):  # the first and second lines' options are still read from the same cell
        assert bt._captain_order_names("R22", bt.T1_R22_OPTIONS[v], ledger)
    for v in (1.5, 2.0):
        assert bt._captain_order_names("R22", bt.T1_R22_DAMPING_OPTIONS[v], ledger)
    assert not bt._captain_order_names("R22", "`--g-mom-damping 1.0 --feng-fraction-mult 3.0`", ledger)
    # the row as it stood before the third line: it names the bare `--feng-fraction-mult 2.5`, never the pair
    ledger.write_text(f"# ledger\n{ROW_R20}\n{ROW_R22_BEFORE}\n")
    assert "`--feng-fraction-mult 2.5`" in ROW_R22_BEFORE and "--g-mom-damping 1.0" not in ROW_R22_BEFORE
    for v in (2.0, 2.5):
        assert not bt._captain_order_names("R22", bt.T1_R22_LINE3_OPTIONS[v], ledger)
    assert bt._captain_order_names("R22", bt.T1_R22_OPTIONS[2.5], ledger)


@pytest.mark.parametrize("text, ok20, ok25", [
    (ROW_R22, True, True),
    (f"| R22 (the captain's order of X; ruling, no run{THIRD}) | T1 |", True, True),
    ("| R22 (the captain's order of X; ruling, no run; `--g-mom-damping 1.0 --feng-fraction-mult 2.0` with "
     "`--captain-ruling R22`) |", True, False),  # names only 2.0's pair
    ("| R22 (the captain's order of X; ruling, no run; `--g-mom-damping 1.0` `--feng-fraction-mult 2.0` "
     "`--feng-fraction-mult 2.5` `--captain-ruling R22`) |", False, False),  # the parts, not the combined option
    ("| R22 (the captain's order of X; ruling, no run; `--g-mom-damping 1 --feng-fraction-mult 2` "
     "`--g-mom-damping 1.0 --feng-fraction-mult 2.5` `--captain-ruling R22`) |", False, True),  # as the row writes it
    ("| R22 (the captain's order of X; ruling, no run; --g-mom-damping 1.0 --feng-fraction-mult 2.0 and "
     "--g-mom-damping 1.0 --feng-fraction-mult 2.5 `--captain-ruling R22`) |", False, False),  # not in backticks
    ("| R22 (the captain's order of X; a try; `--g-mom-damping 1.0 --feng-fraction-mult 2.0` "
     "`--g-mom-damping 1.0 --feng-fraction-mult 2.5` `--captain-ruling R22`) |", False, False),
    ("| R22 awaits (the captain's order of X; ruling, no run; `--g-mom-damping 1.0 --feng-fraction-mult 2.0` "
     "`--g-mom-damping 1.0 --feng-fraction-mult 2.5` `--captain-ruling R22`) |", False, False),
    ("| R22 (the captain's order of X; ruling, no run) | T1 | `--g-mom-damping 1.0 --feng-fraction-mult 2.0` "
     "`--g-mom-damping 1.0 --feng-fraction-mult 2.5` `--captain-ruling R22` |", False, False),  # outside the cell
    ("| R22 (the captain's order of X; ruling, no run; `--g-mom-damping 1.0 --feng-fraction-mult 2.0` "
     "`--g-mom-damping 1.0 --feng-fraction-mult 2.5`) |", False, False),
    ("| R22 (the captain's order of X; ruling, no run; `--g-mom-damping 1.0 --feng-fraction-mult 2.0` "
     "`--g-mom-damping 1.0 --feng-fraction-mult 2.5` `--captain-ruling R20`) |", False, False),
    ("| R23 (the captain's order of X; ruling, no run; `--g-mom-damping 1.0 --feng-fraction-mult 2.0` "
     "`--g-mom-damping 1.0 --feng-fraction-mult 2.5` `--captain-ruling R22`) |", False, False),
])
def test_the_order_check_takes_only_a_single_ruling_row_naming_each_combined_option(bt, tmp_path, text, ok20, ok25):
    p = tmp_path / "TRIES.md"
    p.write_text(f"# ledger\n{text}\n")
    assert bt._captain_order_names("R22", bt.T1_R22_LINE3_OPTIONS[2.0], p) is ok20
    assert bt._captain_order_names("R22", bt.T1_R22_LINE3_OPTIONS[2.5], p) is ok25
    p.write_text(f"# ledger\n{text}\n{text}\n")  # two rows with one ID: never
    assert not bt._captain_order_names("R22", bt.T1_R22_LINE3_OPTIONS[2.0], p)
    assert not bt._captain_order_names("R22", bt.T1_R22_LINE3_OPTIONS[2.5], p)


# ------------------------------------------------------------------ nothing existing changes


@pytest.mark.parametrize("name", list(BEFORE))
def test_the_existing_commands_give_the_same_bytes_and_records(bt, rc4, ledger, tmp_path, name):
    extra, file_sha, record_sha = BEFORE[name]
    f, rec = _build(bt, tmp_path, f"T1_rc{name}", extra)
    assert _sha(f) == file_sha and _record_sha(rec) == record_sha
    assert "captain_ruling" not in rec.get("t1_candidate", {})


def test_the_ruled_files_of_r20_and_the_first_two_lines_give_the_bytes_and_records_as_before(bt, rc4, ledger,
                                                                                            tmp_path):
    f, rec = _build(bt, tmp_path, "T1_s2rcB2", R20)
    assert (_sha(f), _record_sha(rec)) == R20_BEFORE
    for v, name in ((2.5, "T1_s2rcB25"), (3.0, "T1_s2rcB30")):
        f, rec = _build(bt, tmp_path, name, LINE1[v])
        assert (_sha(f), _record_sha(rec)) == LINE1_BEFORE[v]
        # an explicit damping at its default (0.5) changes nothing on the first line
        fd, recd = _build(bt, tmp_path, f"{name}_d", [*LINE1[v], "--g-mom-damping", "0.5"])
        assert (_sha(fd), _record_sha(recd)) == LINE1_BEFORE[v]
    for v, name in ((1.5, "T1_s2rcA15"), (2.0, "T1_s2rcA20")):
        f, rec = _build(bt, tmp_path, name, LINE2[v])
        assert (_sha(f), _record_sha(rec)) == LINE2_BEFORE[v]
    fd, recd = _build(bt, tmp_path, "T1_s2rcB2_d", [*R20, "--g-mom-damping", "0.5"])  # row R20's at damping 0.5
    assert (_sha(fd), _record_sha(recd)) == R20_BEFORE


def test_a_b_and_c_need_no_ruling(bt, rc4, ledger, tmp_path):
    ledger.write_text("# ledger\n")  # no ruling row at all: A, B and C build all the same
    for name in ("A", "B", "C"):
        extra, file_sha, record_sha = BEFORE[name]
        f, rec = _build(bt, tmp_path, f"T1_rc{name}_nr", extra)
        assert (_sha(f), _record_sha(rec)) == (file_sha, record_sha)
        assert "captain_ruling" not in rec["t1_candidate"]


# ------------------------------------------------------------------ the line's two files


@pytest.mark.parametrize("v", [2.0, 2.5])
def test_the_line_file_is_file_c_with_only_the_steps_size_changed(bt, rc4, ledger, tmp_path, v):
    fc, rc = _build(bt, tmp_path, "T1_s2rcC", C)
    fv, rv = _build(bt, tmp_path, NAMES[v], LINE3[v])
    calendar = bt.feng_dirs.bracket_fraction(10.5, 12.5)
    assert calendar == 0.4
    assert rc4[0]["fraction"] == rc["fraction"] == pytest.approx(1.5 * calendar, abs=1e-12)
    assert rc4[-1]["fraction"] == rv["fraction"] == v * calendar == pytest.approx(v * 0.4, abs=1e-12)
    # what build_submission is given besides the step's size is C's
    assert rc4[0]["groups"] == rc4[-1]["groups"] and rc4[0]["directions"] == rc4[-1]["directions"]
    b0, bv = rc4[0]["before"], rc4[-1]["before"]
    assert list(b0.obs_names) == list(bv.obs_names) and b0.obs.equals(bv.obs)
    cand = rv["t1_candidate"]
    assert cand["g_mom_damping"] == 1.0 and cand["feng_fraction_mult"] == v and cand["captain_ruling"] == "R22"
    assert "row R22" in cand["order"] and "decision log 2026-11-22" in cand["order"]
    assert "third T1 search line" in cand["order"] and f"direction at {v:.1f} x the calendar fraction" in cand["order"]
    assert "damping held at 1.0" in cand["order"] and "file C's command" in cand["order"]
    assert cand["calendar_fraction"] == calendar and cand["applied_fraction"] == rv["fraction"]
    assert cand["command"] == bt.T1_CANDIDATE_COMMAND
    assert set(cand) == {"g_mom_damping", "feng_fraction_mult", "captain_ruling", "order", "calendar_fraction",
                         "applied_fraction", "command"}  # the ruled files' fields, no other
    assert "released_bet" not in rv
    # the record is C's but for the step's size and the candidate field; g_mom's part (damping 1.0, candidate A's
    # text) is C's exactly
    assert _strip(rv) == _strip(rc) and rv["g_mom"] == rc["g_mom"] and rv["g_mom"]["damping"] == 1.0
    assert {k: rc["t1_candidate"][k] for k in ("g_mom_damping", "calendar_fraction", "command")} == \
        {k: cand[k] for k in ("g_mom_damping", "calendar_fraction", "command")}
    assert _sha(fv) not in {*(b[1] for b in BEFORE.values()), R20_BEFORE[0], *(b[0] for b in LINE1_BEFORE.values()),
                            *(b[0] for b in LINE2_BEFORE.values())}


@pytest.mark.parametrize("v", [2.0, 2.5])
def test_the_line_file_changes_only_cells_with_a_direction_and_builds_the_same_bytes_twice(bt, rc3, ledger, tmp_path,
                                                                                           v):
    fc, _ = _build(bt, tmp_path, "T1_s2rcC", C)
    fa, _ = _build(bt, tmp_path, f"{NAMES[v]}a", LINE3[v])
    fb, _ = _build(bt, tmp_path, f"{NAMES[v]}b", LINE3[v])
    assert fa.read_bytes() == fb.read_bytes()
    xc, xv = _same_frame(fc, fa)  # same cells, order, obs and var as C
    changed = _rows_changed(xc, xv)
    assert changed.any()
    pop = rc3[-1]["before"]
    labels = pop.obs["celltype"].astype(str).to_numpy()
    groups = rc3[-1]["groups"]
    with_dir = np.array([groups.get(lab) in THREE for lab in labels])
    assert not with_dir.all() and not (changed & ~with_dir).any()  # Epicardium and the rest keep C's counts
    b0 = rc3[0]["before"]  # the cells, their order and obs are C's
    assert list(b0.obs_names) == list(pop.obs_names) and b0.obs.equals(pop.obs)


def test_the_two_values_give_two_files_and_2_is_2_0_and_neither_is_an_earlier_file(bt, rc4, ledger, tmp_path):
    f20, r20 = _build(bt, tmp_path, NAMES[2.0], LINE3[2.0])
    f25, _ = _build(bt, tmp_path, NAMES[2.5], LINE3[2.5])
    f2, r2 = _build(bt, tmp_path, "T1_s2rcCm2", ["--g-mom-damping", "1", "--feng-fraction-mult", "2",
                                                  "--captain-ruling", "R22"])
    fx, _ = _build(bt, tmp_path, "T1_s2rcB2", R20)  # row R20's 2.0 file, damping 0.5
    fy, _ = _build(bt, tmp_path, "T1_s2rcB25", LINE1[2.5])  # the first line's 2.5 file, damping 0.5
    fc, _ = _build(bt, tmp_path, "T1_s2rcC", C)
    assert f2.read_bytes() == f20.read_bytes() and r2["t1_candidate"] == r20["t1_candidate"]
    assert len({f20.read_bytes(), f25.read_bytes(), fx.read_bytes(), fy.read_bytes(), fc.read_bytes()}) == 5


# ------------------------------------------------------------------ refusals, before anything is read


@pytest.mark.parametrize("v", ["2.0", "2.5"])
@pytest.mark.parametrize("extra, match", [
    ([], "is the captain's order of row R22 alone and needs --captain-ruling R22"),
    (["--captain-ruling", "R20"], "refused with --captain-ruling R20: only row R22"),
    (["--captain-ruling", "R21"], "refused with --captain-ruling R21: only row R22"),
    (["--captain-ruling", "R2"], "refused with --captain-ruling R2: only row R22"),
    (["--captain-ruling", "r22"], "refused with --captain-ruling r22: only row R22"),
    (["--captain-ruling", "R220"], "refused with --captain-ruling R220: only row R22"),
    (["--captain-ruling", "T63"], "refused with --captain-ruling T63: only row R22"),
    (["--captain-ruling", "R22", "--bet-strength", "0.5"], "with --released-bet and --bet-strength"),
    (["--captain-ruling", "R22", "--released-bet", "share_trend", "--bet-strength", "0.5"],
     "with --released-bet and --bet-strength"),
    (["--captain-ruling", "R22", "--seed", "1000002"], "recomputed last file only"),
    (["--captain-ruling", "R22", "--board", "T1:val"], "recomputed last file only"),  # a dry run
])
def test_every_other_use_of_a_line_pair_is_refused(bt, ledger, monkeypatch, tmp_path, v, extra, match):
    _refused(bt, monkeypatch, tmp_path, ["t1", *RC, "--g-mom-damping", "1.0", "--feng-fraction-mult", v, *extra],
             match)


@pytest.mark.parametrize("v", ["2.0", "2.5"])
@pytest.mark.parametrize("argv", [
    ["t1"],  # the authors' groups
    ["t1", "--g-mom", "--base", "lib3"],  # the first upload's command
    ["t1", "--g-mom", "--base", "lib3", "--feng-groups", "none"],  # the none file
    ["t1", "--g-mom", "--base", "mid", "--cm-share", "0.60"],  # T1's second upload
])
def test_a_line_pair_is_refused_on_any_other_command(bt, ledger, monkeypatch, tmp_path, v, argv):
    _refused(bt, monkeypatch, tmp_path, [*argv, "--g-mom-damping", "1.0", "--feng-fraction-mult", v,
                                         "--captain-ruling", "R22"], "recomputed last file only")


@pytest.mark.parametrize("extra", [["--cm-share", "0.6"], ["--herd"], ["--no-feng-hygiene"], ["--lib3-pair"],
                                   ["--base", "pooled"]])
def test_the_recomputed_commands_own_guards_hold_with_a_line_pair(bt, ledger, monkeypatch, tmp_path, extra):
    for v in ("2.0", "2.5"):
        _refused(bt, monkeypatch, tmp_path, ["t1", *RC, *extra, "--g-mom-damping", "1.0", "--feng-fraction-mult", v,
                                             "--captain-ruling", "R22"], "--feng-groups recomputed builds one file")


def test_a_line_pair_is_refused_when_row_r22_does_not_name_the_combined_option(bt, ledger, monkeypatch, tmp_path):
    ledger.write_text(f"# ledger\n{ROW_R20}\n{ROW_R22_BEFORE}\n")  # the bare `--feng-fraction-mult 2.5` only
    for v in ("2.0", "2.5"):
        _refused(bt, monkeypatch, tmp_path, ["t1", *RC, *LINE3[float(v)]],
                 f"no single row R22 .* `--g-mom-damping 1.0 --feng-fraction-mult {v}`")
    ledger.write_text(f"# ledger\n{ROW_R20}\n| R22 (the captain's order of X; ruling, no run; `--g-mom-damping 1.0` "
                      "`--feng-fraction-mult 2.0` `--feng-fraction-mult 2.5` `--captain-ruling R22`) | T1 |\n")
    _refused(bt, monkeypatch, tmp_path, ["t1", *RC, *LINE3[2.0]], "no single row R22")  # the parts, apart
    ledger.write_text(f"# ledger\n{ROW_R20}\n| R22 (awaits the captain's order; "
                      "`--g-mom-damping 1.0 --feng-fraction-mult 2.0` `--captain-ruling R22`) | T1 |\n")
    _refused(bt, monkeypatch, tmp_path, ["t1", *RC, *LINE3[2.0]], "no single row R22")
    ledger.write_text(f"# ledger\n{ROW_R20}\n")  # no row R22
    _refused(bt, monkeypatch, tmp_path, ["t1", *RC, *LINE3[2.5]], "no single row R22")
    # a row naming 2.0's pair alone orders that pair alone
    ledger.write_text(f"# ledger\n{ROW_R20}\n| R22 (the captain's order of X; ruling, no run; "
                      "`--g-mom-damping 1.0 --feng-fraction-mult 2.0` with `--captain-ruling R22`) | T1 |\n")
    _refused(bt, monkeypatch, tmp_path, ["t1", *RC, *LINE3[2.5]],
             "no single row R22 .* `--g-mom-damping 1.0 --feng-fraction-mult 2.5`")
    a = bt.build_parser().parse_args(["t1", *RC, *LINE3[2.0], "--out", str(tmp_path / "T1.h5ad")])
    assert bt._t1_candidate(a)["captain_ruling"] == "R22"


@pytest.mark.parametrize("extra, match", [
    # 2.0 at damping 0.5 is row R20's file as before, and refused with R22
    (["--feng-fraction-mult", "2.0", "--captain-ruling", "R22"], "refused with --captain-ruling R22: only row R20"),
    (["--g-mom-damping", "0.5", "--feng-fraction-mult", "2.0", "--captain-ruling", "R22"],
     "refused with --captain-ruling R22: only row R20"),
    (["--feng-fraction-mult", "2.0"], "--feng-fraction-mult takes 1 .* needs --captain-ruling R20"),
    # 2.5 and 3.0 at damping 0.5 under R20 or without a ruling, as before
    (["--feng-fraction-mult", "2.5", "--captain-ruling", "R20"], "refused with --captain-ruling R20: only row R22"),
    (["--g-mom-damping", "0.5", "--feng-fraction-mult", "3.0"], "--feng-fraction-mult takes 1 .* R22"),
    # 3.0 with damping 1.0: no line takes it
    (["--g-mom-damping", "1.0", "--feng-fraction-mult", "3.0", "--captain-ruling", "R22"],
     "with --g-mom-damping other than 0.5"),
    (["--g-mom-damping", "1.0", "--feng-fraction-mult", "3.0"], "needs --captain-ruling R22"),
    # A, B and C with a ruling: refused as before
    (["--g-mom-damping", "1.0", "--feng-fraction-mult", "1.5", "--captain-ruling", "R22"], "with --released-bet only"),
    (["--feng-fraction-mult", "1.5", "--captain-ruling", "R22"], "with --released-bet only"),
    (["--g-mom-damping", "1.0", "--captain-ruling", "R22"], "with --released-bet only"),
    # no value between the half steps and none above 2.5 on C's command (R22: --captain-ruling R22)
    (["--g-mom-damping", "1.0", "--feng-fraction-mult", "2.25", *R22], "--feng-fraction-mult takes"),
    (["--g-mom-damping", "1.0", "--feng-fraction-mult", "1.75", *R22], "--feng-fraction-mult takes"),
    (["--g-mom-damping", "1.0", "--feng-fraction-mult", "2.75", *R22], "--feng-fraction-mult takes"),
    (["--g-mom-damping", "1.0", "--feng-fraction-mult", "3.5", *R22], "--feng-fraction-mult takes"),
])
def test_the_neighbouring_pairs_are_taken_or_refused_as_before(bt, ledger, monkeypatch, tmp_path, extra, match):
    _refused(bt, monkeypatch, tmp_path, ["t1", *RC, *extra], match)


@pytest.mark.parametrize("ruling", [[], ["--captain-ruling", "R22"], ["--captain-ruling", "R20"]])
@pytest.mark.parametrize("damping", ["0.75", "1.25", "1.5", "2.0", "0", "-1.0"])
@pytest.mark.parametrize("mult", ["1.75", "2.0", "2.5", "3.0"])
def test_any_multiplier_above_1_5_with_a_damping_other_than_0_5_or_1_0_is_refused(bt, ledger, monkeypatch, tmp_path,
                                                                                  ruling, damping, mult):
    _refused(bt, monkeypatch, tmp_path, ["t1", *RC, "--g-mom-damping", damping, "--feng-fraction-mult", mult,
                                         *ruling], r"refused|takes|needs")


@pytest.mark.parametrize("damping", ["1.5", "2.0"])
@pytest.mark.parametrize("mult", ["0.8", "1.5", "2.0", "2.5", "3.0"])
def test_the_second_lines_dampings_with_any_multiplier_but_1_are_refused(bt, ledger, monkeypatch, tmp_path, damping,
                                                                        mult):
    _refused(bt, monkeypatch, tmp_path, ["t1", *RC, "--g-mom-damping", damping, "--feng-fraction-mult", mult,
                                         "--captain-ruling", "R22"], r"refused|takes")


def test_the_bet_refusal_lets_only_the_ruled_commands_through(bt, ledger):
    def ns(**k):
        base = dict(released_bet=None, bet_strength=None, captain_ruling="R22", feng_fraction_mult=2.0,
                    feng_groups="recomputed", seed=harness.BUILD_SEED_OFFSET + 1, board=None, g_mom_damping=1.0)
        return SimpleNamespace(**{**base, **k})

    for good in (ns(), ns(feng_fraction_mult=2.5), ns(g_mom_damping=None, feng_fraction_mult=2.5),
                 ns(g_mom_damping=None, captain_ruling="R20"), ns(g_mom_damping=1.5, feng_fraction_mult=None)):
        assert bt._bet_refusals(good) is None
    # (the damping beside row R20's 2.0 or the first line's 3.0 is _t1_candidate's to refuse, and it runs first)
    for bad in (ns(feng_fraction_mult=3.0, captain_ruling="R20"), ns(feng_fraction_mult=1.5),
                ns(feng_fraction_mult=None), ns(feng_fraction_mult=2.25), ns(captain_ruling="R21"),
                ns(g_mom_damping=0.75), ns(g_mom_damping=1.5), ns(feng_fraction_mult=2.5, captain_ruling="R20"),
                ns(bet_strength=0.5)):
        with pytest.raises(SystemExit, match="with --released-bet only"):
            bt._bet_refusals(bad)
    ledger.write_text(f"# ledger\n{ROW_R20}\n{ROW_R22_BEFORE}\n")  # without the combined options the shortcut closes
    with pytest.raises(SystemExit, match="with --released-bet only"):  # (2.5's pair: _t1_candidate refuses it first)
        bt._bet_refusals(ns())
    assert bt._bet_refusals(ns(g_mom_damping=None, feng_fraction_mult=2.5)) is None  # the first line's, as before
