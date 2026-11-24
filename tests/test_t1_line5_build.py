"""`build_test.py t1 --g-mom-damping 1.5|2.0 --feng-fraction-mult 2.0 --captain-ruling R22`: the fifth T1 search line
under row R22 of experiments/harness/TRIES.md (T1 line 5 of docs/PLAN.md's decision log): the 2.0 file's command (row
R13's recomputed-groups file with the outside atlas's direction step at 2.0 x the calendar fraction) with g_mom's
within-type damping at 1.5 or 2.0, nothing else changed (files T1_s2rcB2d15 and T1_s2rcB2d20).

Held here, on synthetic stages and directions only (the build_t1 fixture of tests/test_hidden_counts.py and the
recomputed setup of tests/test_feng_relabel.py), with a temporary ledger in place of the repository's: (1) nothing that
builds today changes: the default command and candidates A, B and C give the bytes and records pinned in
tests/test_t1_fraction2_build.py, row R20's file those pinned in tests/test_t1_fraction_line_build.py, the first line's
files those pinned in tests/test_t1_damping_line_build.py and the second line's files those pinned in
tests/test_t1_line3_build.py (their g_mom text still names the second line); (2) each new file differs from the 2.0 file
only through g_mom's damping; (3) each pair is taken only with --captain-ruling R22 while row R22's first cell is the
captain's order naming the combined option, and every other use is refused before anything is read. No released file,
outside-atlas file, board or score is read."""

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
from test_t1_fraction2_build import BEFORE, R20, _no_read, _record_sha, _refused  # noqa: E402
from test_t1_fraction_line_build import R20_BEFORE  # noqa: E402
from test_t1_fraction_line_build import R22 as LINE1  # noqa: E402
from test_t1_line3_build import FIRST_CELL as FIRST_CELL_LINE3  # noqa: E402
from test_t1_line3_build import LINE2_BEFORE, LINE3, ROW_R20  # noqa: E402

from vec import harness  # noqa: E402

bt, rc3, rc4 = tcb.bt, tcb.rc3, tcb.rc4  # the candidates' fixtures (the builder module; the recomputed setups)

LINE5 = {1.5: ["--g-mom-damping", "1.5", "--feng-fraction-mult", "2.0", "--captain-ruling", "R22"],
         2.0: ["--g-mom-damping", "2.0", "--feng-fraction-mult", "2.0", "--captain-ruling", "R22"]}
NAMES = {1.5: "T1_s2rcB2d15", 2.0: "T1_s2rcB2d20"}
OPT15 = "`--g-mom-damping 1.5 --feng-fraction-mult 2.0`"
OPT20 = "`--g-mom-damping 2.0 --feng-fraction-mult 2.0`"
FIFTH = (f"; by the fifth T1 search line the T1 builder is also to take {OPT15} and {OPT20}, each only with "
         "`--captain-ruling R22`")
FIRST_CELL = FIRST_CELL_LINE3[:-1] + FIFTH + ")"  # row R22's first cell with the fifth line's amendment
ROW_R22 = f"| {FIRST_CELL} | T1 | the lines |"
ROW_R22_BEFORE = f"| {FIRST_CELL_LINE3} | T1 | the lines |"  # the row before the fifth line: the bare options only


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
    """A record without what the damping sets: the candidate field and g_mom's part."""
    return {k: v for k, v in r.items() if k not in ("t1_candidate", "g_mom")}


# ------------------------------------------------------------------ the constants and the ruling check


def test_the_line_has_its_own_constants_and_the_earlier_ones_are_unchanged(bt):
    assert bt.G_MOM_DAMPINGS == (0.5, 1.0) and bt.FENG_FRACTION_MULTS == (1.0, 1.5)
    assert bt.FENG_FRACTION_MULT_LINE5 == 2.0 and bt.G_MOM_DAMPINGS_LINE5 == (1.5, 2.0)
    assert bt.T1_R22_LINE5_OPTIONS == {1.5: OPT15, 2.0: OPT20}
    assert bt.T1_R22_LINE3_OPTIONS == {2.0: "`--g-mom-damping 1.0 --feng-fraction-mult 2.0`",
                                       2.5: "`--g-mom-damping 1.0 --feng-fraction-mult 2.5`"}
    assert bt.T1_R22_DAMPING_OPTIONS == {1.5: "`--g-mom-damping 1.5`", 2.0: "`--g-mom-damping 2.0`"}
    assert set(bt.FENG_FRACTION_MULTS_RULED) == {2.0, 2.5, 3.0} and bt.FENG_FRACTION_MULTS_RULED[2.0][0] == "R20"
    assert bt.G_MOM_DAMPINGS_R22 == (1.5, 2.0) and bt.T1_R22_ROW == "R22"


def test_the_help_names_the_fifth_line(bt):
    import io
    from contextlib import redirect_stdout

    buf = io.StringIO()
    with redirect_stdout(buf), pytest.raises(SystemExit):
        bt.build_parser().parse_args(["t1", "--help"])
    text = " ".join(buf.getvalue().split())
    assert "fifth search line" in text and "fifth T1 search line" in text
    assert "fifth T1 search line" in bt.__doc__ and "T1_s2rcB2d15" in bt.__doc__


def test_the_order_check_reads_each_combined_option_from_row_r22s_first_cell(bt, ledger):
    for v in (1.5, 2.0):
        assert bt._captain_order_names("R22", bt.T1_R22_LINE5_OPTIONS[v], ledger)
        assert not bt._captain_order_names("R20", bt.T1_R22_LINE5_OPTIONS[v], ledger)
        assert bt._captain_order_names("R22", bt.T1_R22_DAMPING_OPTIONS[v], ledger)  # the second line's, as before
    for v in (2.0, 2.5):
        assert bt._captain_order_names("R22", bt.T1_R22_LINE3_OPTIONS[v], ledger)
    # the row as it stood before: it names the bare `--g-mom-damping 1.5` and `--g-mom-damping 2.0` and row R20 the
    # bare `--feng-fraction-mult 2.0`, never the pair
    ledger.write_text(f"# ledger\n{ROW_R20}\n{ROW_R22_BEFORE}\n")
    assert "`--g-mom-damping 1.5`" in ROW_R22_BEFORE and "`--feng-fraction-mult 2.0`" in ROW_R20
    for v in (1.5, 2.0):
        assert not bt._captain_order_names("R22", bt.T1_R22_LINE5_OPTIONS[v], ledger)
        assert bt._captain_order_names("R22", bt.T1_R22_DAMPING_OPTIONS[v], ledger)


@pytest.mark.parametrize("text, ok15, ok20", [
    (ROW_R22, True, True),
    (f"| R22 (the captain's order of X; ruling, no run{FIFTH}) | T1 |", True, True),
    (f"| R22 (the captain's order of X; ruling, no run; {OPT15} with `--captain-ruling R22`) |", True, False),
    ("| R22 (the captain's order of X; ruling, no run; `--g-mom-damping 1.5` `--g-mom-damping 2.0` "
     "`--feng-fraction-mult 2.0` `--captain-ruling R22`) |", False, False),  # the parts, not the combined option
    ("| R22 (the captain's order of X; ruling, no run; --g-mom-damping 1.5 --feng-fraction-mult 2.0 and "
     "--g-mom-damping 2.0 --feng-fraction-mult 2.0 `--captain-ruling R22`) |", False, False),  # not in backticks
    (f"| R22 (the captain's order of X; a try; {OPT15} {OPT20} `--captain-ruling R22`) |", False, False),
    (f"| R22 (the captain's order of X; ruling, no run) | T1 | {OPT15} {OPT20} `--captain-ruling R22` |", False,
     False),  # outside the first cell
    (f"| R22 (the captain's order of X; ruling, no run; {OPT15} {OPT20}) |", False, False),
    (f"| R22 (the captain's order of X; ruling, no run; {OPT15} {OPT20} `--captain-ruling R20`) |", False, False),
    (f"| R23 (the captain's order of X; ruling, no run; {OPT15} {OPT20} `--captain-ruling R22`) |", False, False),
])
def test_the_order_check_takes_only_a_single_ruling_row_naming_each_combined_option(bt, tmp_path, text, ok15, ok20):
    p = tmp_path / "TRIES.md"
    p.write_text(f"# ledger\n{text}\n")
    assert bt._captain_order_names("R22", OPT15, p) is ok15
    assert bt._captain_order_names("R22", OPT20, p) is ok20
    p.write_text(f"# ledger\n{text}\n{text}\n")  # two rows with one ID: never
    assert not bt._captain_order_names("R22", OPT15, p) and not bt._captain_order_names("R22", OPT20, p)


# ------------------------------------------------------------------ nothing existing changes


@pytest.mark.parametrize("name", list(BEFORE))
def test_the_existing_commands_give_the_same_bytes_and_records(bt, rc4, ledger, tmp_path, name):
    extra, file_sha, record_sha = BEFORE[name]
    f, rec = _build(bt, tmp_path, f"T1_rc{name}", extra)
    assert _sha(f) == file_sha and _record_sha(rec) == record_sha


def test_the_ruled_files_of_r20_and_the_first_two_lines_give_the_bytes_and_records_as_before(bt, rc4, ledger,
                                                                                            tmp_path):
    f, rec = _build(bt, tmp_path, "T1_s2rcB2", R20)
    assert (_sha(f), _record_sha(rec)) == R20_BEFORE
    for v, name in ((2.5, "T1_s2rcB25"), (3.0, "T1_s2rcB30")):
        f, rec = _build(bt, tmp_path, name, LINE1[v])
        assert (_sha(f), _record_sha(rec)) == LINE1_BEFORE[v]
    for v, name in ((1.5, "T1_s2rcA15"), (2.0, "T1_s2rcA20")):  # the second line: bytes, record and text as before
        f, rec = _build(bt, tmp_path, name, LINE2[v])
        assert (_sha(f), _record_sha(rec)) == LINE2_BEFORE[v]
        assert "second T1 search line" in rec["g_mom"]["rule"] and "fifth" not in rec["g_mom"]["rule"]
        assert "second T1 search line" in rec["t1_candidate"]["order"] and "fifth" not in str(rec)


# ------------------------------------------------------------------ the line's two files


@pytest.mark.parametrize("v", [1.5, 2.0])
def test_the_line_file_is_the_2_0_file_with_only_g_moms_damping_changed(bt, rc4, ledger, tmp_path, v):
    fx, rx = _build(bt, tmp_path, "T1_s2rcB2", R20)
    fv, rv = _build(bt, tmp_path, NAMES[v], LINE5[v])
    fa, ra = _build(bt, tmp_path, f"A{v}", LINE2[v])  # the second line's file at the same damping, multiplier 1.0
    calendar = bt.feng_dirs.bracket_fraction(10.5, 12.5)
    assert rc4[0]["fraction"] == rc4[1]["fraction"] == rx["fraction"] == rv["fraction"] == 2.0 * calendar
    # what build_submission is given is the 2.0 file's: the same direction step, groups and directions
    assert rc4[0]["groups"] == rc4[1]["groups"] and rc4[0]["directions"] == rc4[1]["directions"]
    b0, bv = rc4[0]["before"], rc4[1]["before"]
    assert list(b0.obs_names) == list(bv.obs_names) and b0.obs.equals(bv.obs)
    cand = rv["t1_candidate"]
    assert cand["g_mom_damping"] == v and cand["feng_fraction_mult"] == 2.0 and cand["captain_ruling"] == "R22"
    assert "row R22" in cand["order"] and "fifth T1 search line" in cand["order"]
    assert "2.0 file's command" in cand["order"] and f"damping at {v:.1f}" in cand["order"]
    assert "second" not in cand["order"] and "third" not in cand["order"]
    assert cand["calendar_fraction"] == calendar and cand["applied_fraction"] == rv["fraction"]
    assert cand["command"] == bt.T1_CANDIDATE_COMMAND
    assert set(cand) == {"g_mom_damping", "feng_fraction_mult", "captain_ruling", "order", "calendar_fraction",
                         "applied_fraction", "command"}  # the ruled files' fields, no other
    # the record is the 2.0 file's but for the candidate field and g_mom's damping and text
    assert _strip(rv) == _strip(rx)
    gx, gv, ga = rx["g_mom"], rv["g_mom"], ra["g_mom"]
    assert gx["damping"] == 0.5 and gv["damping"] == ga["damping"] == v
    assert {k: x for k, x in gv.items() if k not in ("damping", "rule")} == \
        {k: x for k, x in gx.items() if k not in ("damping", "rule")}
    assert f"damped by {v:g} in place of B0_DAMPING (row R22's fifth T1 search line" in gv["rule"]
    assert "second" not in gv["rule"] and "candidate A" not in gv["rule"]
    assert gv["rule"].replace(gv["rule"][gv["rule"].index("damped by"):gv["rule"].index(", zeroed")],
                              "damped by B0_DAMPING") == gx["rule"]
    assert "second T1 search line" in ga["rule"]  # the second line's text at multiplier 1.0
    known = {*(b[1] for b in BEFORE.values()), R20_BEFORE[0], *(b[0] for b in LINE1_BEFORE.values()),
             *(b[0] for b in LINE2_BEFORE.values())}
    assert _sha(fv) not in known and _sha(fa) == LINE2_BEFORE[v][0]


@pytest.mark.parametrize("v", [1.5, 2.0])
def test_the_line_file_changes_only_cells_with_a_direction_and_builds_the_same_bytes_twice(bt, rc3, ledger, tmp_path,
                                                                                           v):
    fx, _ = _build(bt, tmp_path, "T1_s2rcB2", R20)
    fa, _ = _build(bt, tmp_path, f"{NAMES[v]}a", LINE5[v])
    fb, _ = _build(bt, tmp_path, f"{NAMES[v]}b", LINE5[v])
    assert fa.read_bytes() == fb.read_bytes()
    xx, xv = _same_frame(fx, fa)  # same cells, order, obs and var as the 2.0 file
    changed = _rows_changed(xx, xv)
    assert changed.any()
    pop = rc3[1]["before"]
    labels = pop.obs["celltype"].astype(str).to_numpy()
    groups = rc3[1]["groups"]
    with_dir = np.array([groups.get(lab) in THREE for lab in labels])
    assert not with_dir.all() and not (changed & ~with_dir).any()  # Epicardium and the rest keep the 2.0 file's counts
    b0 = rc3[0]["before"]  # the population g_mom's hook gets (after the direction step) is the 2.0 file's
    assert list(b0.obs_names) == list(pop.obs_names) and b0.obs.equals(pop.obs)
    assert np.array_equal(tcb._dense(b0.X), tcb._dense(pop.X))


def test_the_two_values_give_two_files_and_2_is_2_0_and_neither_is_an_earlier_file(bt, rc4, ledger, tmp_path):
    f15, _ = _build(bt, tmp_path, NAMES[1.5], LINE5[1.5])
    f20, r20 = _build(bt, tmp_path, NAMES[2.0], LINE5[2.0])
    f2, r2 = _build(bt, tmp_path, "T1_s2rcB2d2", ["--g-mom-damping", "2", "--feng-fraction-mult", "2",
                                                  "--captain-ruling", "R22"])
    fx, _ = _build(bt, tmp_path, "T1_s2rcB2", R20)
    fc, _ = _build(bt, tmp_path, "T1_s2rcCm20", LINE3[2.0])  # the third line's 2.0 file, damping 1.0
    assert f2.read_bytes() == f20.read_bytes() and r2["t1_candidate"] == r20["t1_candidate"]
    assert len({f15.read_bytes(), f20.read_bytes(), fx.read_bytes(), fc.read_bytes()}) == 4


# ------------------------------------------------------------------ refusals, before anything is read


@pytest.mark.parametrize("v", ["1.5", "2.0"])
@pytest.mark.parametrize("extra, match", [
    ([], "is the captain's order of row R22 alone and needs --captain-ruling R22"),
    (["--captain-ruling", "R20"], "refused with --captain-ruling R20: only row R22"),
    (["--captain-ruling", "R23"], "refused with --captain-ruling R23: only row R22"),
    (["--captain-ruling", "R26"], "refused with --captain-ruling R26: only row R22"),
    (["--captain-ruling", "R21"], "refused with --captain-ruling R21: only row R22"),
    (["--captain-ruling", "r22"], "refused with --captain-ruling r22: only row R22"),
    (["--captain-ruling", "R220"], "refused with --captain-ruling R220: only row R22"),
    (["--captain-ruling", "T63"], "refused with --captain-ruling T63: only row R22"),
    (["--captain-ruling", "R22", "--bet-strength", "0.5"], "with --released-bet and --bet-strength"),
    (["--captain-ruling", "R22", "--released-bet", "share_trend", "--bet-strength", "0.5"],
     "with --released-bet and --bet-strength"),
    (["--captain-ruling", "R22", "--no-hold-unbracketed"], "refused"),
    (["--captain-ruling", "R22", "--seed", "1000002"], "recomputed last file only"),
    (["--captain-ruling", "R22", "--board", "T1:val"], "recomputed last file only"),  # a dry run
])
def test_every_other_use_of_a_line_pair_is_refused(bt, ledger, monkeypatch, tmp_path, v, extra, match):
    _refused(bt, monkeypatch, tmp_path, ["t1", *RC, "--g-mom-damping", v, "--feng-fraction-mult", "2.0", *extra],
             match)


@pytest.mark.parametrize("v", ["1.5", "2.0"])
@pytest.mark.parametrize("ruling", ["R22", "R23", "R26"])
def test_a_line_pair_is_refused_with_direction_momentum(bt, ledger, monkeypatch, tmp_path, v, ruling):
    """Candidate MS's option never stacks on the pair (its own refusal names row R23's numbers, so it is matched here
    without _refused's number check)."""
    monkeypatch.setattr(bt, "load_stage", _no_read)
    monkeypatch.setattr(bt.feng_dirs, "feng_directions", _no_read)
    with pytest.raises(SystemExit, match="refused"):
        bt.main(["t1", *RC, "--g-mom-damping", v, "--feng-fraction-mult", "2.0", "--direction-momentum", "all",
                 "--captain-ruling", ruling, "--out", str(tmp_path / "T1.h5ad")])


@pytest.mark.parametrize("v", ["1.5", "2.0"])
@pytest.mark.parametrize("argv", [
    ["t1"],  # the authors' groups
    ["t1", "--g-mom", "--base", "lib3"],  # the first upload's command
    ["t1", "--g-mom", "--base", "lib3", "--feng-groups", "none"],  # the none file
    ["t1", "--g-mom", "--base", "mid", "--cm-share", "0.60"],  # T1's second upload
])
def test_a_line_pair_is_refused_on_any_other_command(bt, ledger, monkeypatch, tmp_path, v, argv):
    _refused(bt, monkeypatch, tmp_path, [*argv, "--g-mom-damping", v, "--feng-fraction-mult", "2.0",
                                         "--captain-ruling", "R22"], "recomputed last file only")


@pytest.mark.parametrize("extra", [["--cm-share", "0.6"], ["--herd"], ["--no-feng-hygiene"], ["--lib3-pair"],
                                   ["--base", "pooled"]])
def test_the_recomputed_commands_own_guards_hold_with_a_line_pair(bt, ledger, monkeypatch, tmp_path, extra):
    for v in ("1.5", "2.0"):
        _refused(bt, monkeypatch, tmp_path, ["t1", *RC, *extra, *LINE5[float(v)]],
                 "--feng-groups recomputed builds one file")


def test_a_line_pair_is_refused_when_row_r22_does_not_name_the_combined_option(bt, ledger, monkeypatch, tmp_path):
    ledger.write_text(f"# ledger\n{ROW_R20}\n{ROW_R22_BEFORE}\n")  # the bare options only
    for v in ("1.5", "2.0"):
        _refused(bt, monkeypatch, tmp_path, ["t1", *RC, *LINE5[float(v)]],
                 f"no single row R22 .* `--g-mom-damping {v} --feng-fraction-mult 2.0`")
    ledger.write_text(f"# ledger\n{ROW_R20}\n| R22 (the captain's order of X; ruling, no run; `--g-mom-damping 1.5` "
                      "`--g-mom-damping 2.0` `--feng-fraction-mult 2.0` `--captain-ruling R22`) | T1 |\n")
    _refused(bt, monkeypatch, tmp_path, ["t1", *RC, *LINE5[1.5]], "no single row R22")  # the parts, apart
    ledger.write_text(f"# ledger\n{ROW_R20}\n| R22 (awaits the captain's order; {OPT15} `--captain-ruling R22`) |\n")
    _refused(bt, monkeypatch, tmp_path, ["t1", *RC, *LINE5[1.5]], "no single row R22")
    ledger.write_text(f"# ledger\n{ROW_R20}\n")  # no row R22
    _refused(bt, monkeypatch, tmp_path, ["t1", *RC, *LINE5[2.0]], "no single row R22")
    # a row naming 1.5's pair alone orders that pair alone
    ledger.write_text(f"# ledger\n{ROW_R20}\n| R22 (the captain's order of X; ruling, no run; {OPT15} with "
                      "`--captain-ruling R22`) | T1 |\n")
    _refused(bt, monkeypatch, tmp_path, ["t1", *RC, *LINE5[2.0]],
             "no single row R22 .* `--g-mom-damping 2.0 --feng-fraction-mult 2.0`")
    a = bt.build_parser().parse_args(["t1", *RC, *LINE5[1.5], "--out", str(tmp_path / "T1.h5ad")])
    assert bt._t1_candidate(a)["captain_ruling"] == "R22"


def test_the_pair_is_taken_under_r22_with_a_ledger_naming_the_combined_option(bt, ledger, tmp_path):
    for v in (1.5, 2.0):
        a = bt.build_parser().parse_args(["t1", *RC, *LINE5[v], "--out", str(tmp_path / "T1.h5ad")])
        cand = bt._t1_candidate(a)
        assert (cand["g_mom_damping"], cand["feng_fraction_mult"], cand["captain_ruling"]) == (v, 2.0, "R22")
        assert bt._bet_refusals(a) is None and bt._t1_ms(a) is None


@pytest.mark.parametrize("extra, match", [
    # 2.0 at damping 0.5 is row R20's file as before, and refused with R22
    (["--feng-fraction-mult", "2.0", "--captain-ruling", "R22"], "refused with --captain-ruling R22: only row R20"),
    (["--g-mom-damping", "0.5", "--feng-fraction-mult", "2.0", "--captain-ruling", "R22"],
     "refused with --captain-ruling R22: only row R20"),
    (["--feng-fraction-mult", "2.0"], "--feng-fraction-mult takes 1 .* needs --captain-ruling R20"),
    # the fifth line's dampings with any other multiplier above 1.5: refused as before
    (["--g-mom-damping", "1.5", "--feng-fraction-mult", "2.5", "--captain-ruling", "R22"],
     "with --g-mom-damping other than 0.5"),
    (["--g-mom-damping", "2.0", "--feng-fraction-mult", "3.0", "--captain-ruling", "R22"],
     "with --g-mom-damping other than 0.5"),
    # the second line's dampings with the multiplier 1.5: refused as before
    (["--g-mom-damping", "1.5", "--feng-fraction-mult", "1.5", "--captain-ruling", "R22"],
     "with --feng-fraction-mult other than 1"),
    # no value between the half steps on the 2.0 file's command
    (["--g-mom-damping", "1.25", "--feng-fraction-mult", "2.0", "--captain-ruling", "R22"],
     "refused with --captain-ruling R22: only row R20"),
    (["--g-mom-damping", "2.5", "--feng-fraction-mult", "2.0", "--captain-ruling", "R22"],
     "refused with --captain-ruling R22: only row R20"),
    (["--g-mom-damping", "0.75", "--feng-fraction-mult", "2.0", "--captain-ruling", "R20"],
     "with --g-mom-damping other than 0.5"),
    # the third line's pair: as before
    (["--g-mom-damping", "1.0", "--feng-fraction-mult", "2.0", "--captain-ruling", "R20"],
     "refused with --captain-ruling R20: only row R22"),
])
def test_the_neighbouring_pairs_are_taken_or_refused_as_before(bt, ledger, monkeypatch, tmp_path, extra, match):
    _refused(bt, monkeypatch, tmp_path, ["t1", *RC, *extra], match)


def test_the_second_and_third_lines_are_still_taken(bt, ledger, tmp_path):
    for argv, want in ((LINE2[1.5], (1.5, 1.0)), (LINE2[2.0], (2.0, 1.0)), (LINE3[2.0], (1.0, 2.0)),
                       (R20, (0.5, 2.0))):
        a = bt.build_parser().parse_args(["t1", *RC, *argv, "--out", str(tmp_path / "T1.h5ad")])
        cand = bt._t1_candidate(a)
        assert (cand["g_mom_damping"], cand["feng_fraction_mult"]) == want
        assert "fifth" not in cand["order"]


def test_the_bet_refusal_lets_only_the_ruled_commands_through(bt, ledger):
    def ns(**k):
        base = dict(released_bet=None, bet_strength=None, captain_ruling="R22", feng_fraction_mult=2.0,
                    feng_groups="recomputed", seed=harness.BUILD_SEED_OFFSET + 1, board=None, g_mom_damping=1.5)
        return SimpleNamespace(**{**base, **k})

    for good in (ns(), ns(g_mom_damping=2.0), ns(g_mom_damping=1.0), ns(g_mom_damping=None, captain_ruling="R20"),
                 ns(g_mom_damping=1.5, feng_fraction_mult=None)):
        assert bt._bet_refusals(good) is None
    assert bt._ruled_line5_ordered(ns()) and bt._ruled_line5_ordered(ns(g_mom_damping=2.0))
    for bad in (ns(captain_ruling="R20"), ns(captain_ruling="R23"), ns(g_mom_damping=1.25), ns(bet_strength=0.5)):
        assert not bt._ruled_line5_ordered(bad)
    for bad in (ns(captain_ruling="R21"), ns(g_mom_damping=1.25), ns(bet_strength=0.5)):
        with pytest.raises(SystemExit, match="with --released-bet only"):
            bt._bet_refusals(bad)
    ledger.write_text(f"# ledger\n{ROW_R20}\n{ROW_R22_BEFORE}\n")  # without the combined options the shortcut closes
    assert not bt._ruled_line5_ordered(ns()) and not bt._ruled_line5_ordered(ns(g_mom_damping=2.0))
    with pytest.raises(SystemExit, match="with --released-bet only"):
        bt._bet_refusals(ns())
