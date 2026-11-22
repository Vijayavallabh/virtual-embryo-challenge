"""`build_test.py t1 --g-mom-damping 1.5|2.0 --captain-ruling R22`: the second T1 search line under row R22 of
experiments/harness/TRIES.md (the captain's order of 2026-11-22; the lines added in docs/PLAN.md decision log
2026-11-22): the recomputed-groups file (row R13) with g_mom's within-type damping at 1.5 or 2.0 and nothing else
changed. --g-mom-damping 1.0 stays candidate A, with no ruling.

Held here, on synthetic stages and directions only (the build_t1 fixture of tests/test_hidden_counts.py and the
recomputed setup of tests/test_feng_relabel.py), with a temporary ledger in place of the repository's: (1) nothing that
builds today changes: the default command and candidates A, B and C give the bytes and records pinned in
tests/test_t1_fraction2_build.py, row R20's file those pinned in tests/test_t1_fraction_line_build.py, and the first
line's 2.5 and 3.0 files the bytes and records the builder at campaign/t1-fraction-line fc53487 (before this option
existed) gave on this fixture, pinned below (the same kind of pin as tests/test_t1_fraction2_build.py's BEFORE, a
synthetic-fixture pin that experiments/external/n54_step.sh's tests step runs on the DGX as well as here); (2) at 1.5
and 2.0 the only change against the default is g_mom's damping; (3) each value is taken only with --captain-ruling R22
while row R22's first cell is the captain's order naming that value's option, and every other use is refused before
anything is read. No released file, outside-atlas file, board or score is read."""

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
from test_t1_fraction2_build import BEFORE, R20, _record_sha, _refused  # noqa: E402
from test_t1_fraction_line_build import R20_BEFORE  # noqa: E402
from test_t1_fraction_line_build import R22 as LINE1  # noqa: E402

from vec import harness  # noqa: E402

bt, rc3, rc4 = tcb.bt, tcb.rc3, tcb.rc4  # the candidates' fixtures (the builder module; the recomputed setups)

DAMP = {1.5: ["--g-mom-damping", "1.5", "--captain-ruling", "R22"],
        2.0: ["--g-mom-damping", "2.0", "--captain-ruling", "R22"]}
NAMES = {1.5: "T1_s2rcA15", 2.0: "T1_s2rcA20"}
# The first line's two files on this fixture, built by the builder at campaign/t1-fraction-line fc53487 (before the
# damping values existed), twice in separate processes with the same result: the file's sha256 and its record's
# (feng_hygiene.commit set to "own_repo").
LINE1_BEFORE = {2.5: ("4cf67ec7d4533a23f342cabf5e9bafdfadf587863ec63663303c41fa32eeaa61",
                      "f80039201987b00c932192711272690b17550b3a73b03c1db6c8d4b57538aa97"),
                3.0: ("27875e5969ae6054ad9801a43f75d48184b4d187cc9bc0118ef67bb5a3dd75aa",
                      "e8af2cf7d901ea2ad3c2728329ef484e6e138dd4ffdc3b75b33ce25d3be02be9")}
# The sentence the main session adds to row R22's first cell, verbatim, inside a synthetic row of that form
SENTENCE = ("by the lines added in decision log 2026-11-22 the T1 builder is also to take `--g-mom-damping 1.5` and "
            "`--g-mom-damping 2.0`, and the T3 outside-form builder `--table-exponent 4.0`, each only with "
            "`--captain-ruling R22`")
ROW_R20 = ("| R20 (the captain's order of a day; ruling, no run, not a try; the builder is to take "
           "`--feng-fraction-mult 2.0` only with `--captain-ruling R20`) | T1 | a file |")
ROW_R22 = ("| R22 (the captain's order of a later day: a search; ruling, no run, not a try; the builder is to take "
           "`--feng-fraction-mult 2.5` and `--feng-fraction-mult 3.0` only with `--captain-ruling R22`; "
           f"{SENTENCE}) | T1 | the lines |")


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
    """A record without what the damping sets: the candidate field and g_mom's damping and rule text."""
    out = {k: v for k, v in r.items() if k != "t1_candidate"}
    out["g_mom"] = {k: v for k, v in r["g_mom"].items() if k not in ("damping", "rule")}
    return out


# ------------------------------------------------------------------ the constants and the ruling check


def test_the_line_has_its_own_constants_and_candidate_a_is_unchanged(bt):
    assert bt.G_MOM_DAMPINGS == (0.5, 1.0) and bt.FENG_FRACTION_MULTS == (1.0, 1.5)
    assert bt.G_MOM_DAMPINGS_R22 == (1.5, 2.0) and bt.T1_R22_ROW == "R22"
    assert bt.T1_R22_DAMPING_OPTIONS == {1.5: "`--g-mom-damping 1.5`", 2.0: "`--g-mom-damping 2.0`"}
    assert not set(bt.G_MOM_DAMPINGS_R22) & set(bt.G_MOM_DAMPINGS)
    assert bt.FENG_FRACTION_MULTS_R22 == (2.5, 3.0) and set(bt.FENG_FRACTION_MULTS_RULED) == {2.0, 2.5, 3.0}
    a = bt.build_parser().parse_args(["t1", *RC, "--out", "x.h5ad"])
    assert a.g_mom_damping is None and a.captain_ruling is None


def test_the_order_check_reads_each_option_from_row_r22s_first_cell(bt, ledger):
    for v in (1.5, 2.0):
        assert bt._captain_order_names("R22", bt.T1_R22_DAMPING_OPTIONS[v], ledger)
        assert not bt._captain_order_names("R20", bt.T1_R22_DAMPING_OPTIONS[v], ledger)
    for v in (2.5, 3.0):  # the first line's options are still read from the same cell
        assert bt._captain_order_names("R22", bt.T1_R22_OPTIONS[v], ledger)
    assert not bt._captain_order_names("R22", "`--g-mom-damping 1.0`", ledger)
    assert not bt._captain_ruling_names("R22", "share_trend_l0.5", ledger)


@pytest.mark.parametrize("text, ok15, ok20", [
    (ROW_R22, True, True),
    (f"| R22 (the captain's order of X; ruling, no run; {SENTENCE}) | T1 |", True, True),
    ("| R22 (the captain's order of X; ruling, no run; `--g-mom-damping 1.5` with `--captain-ruling R22`) |",
     True, False),  # names only 1.5: 2.0 is not ordered
    ("| R22 (the captain's order of X; ruling, no run; `--g-mom-damping 1.5` `--g-mom-damping 2` "
     "`--captain-ruling R22`) |", True, False),  # the option as the row must write it: 2.0
    ("| R22 (the captain's order of X; ruling, no run; --g-mom-damping 1.5 and --g-mom-damping 2.0 "
     "`--captain-ruling R22`) |", False, False),  # not in backticks
    ("| R22 (the captain's order of X; a try; `--g-mom-damping 1.5` `--g-mom-damping 2.0` `--captain-ruling R22`) |",
     False, False),
    ("| R22 (the captain's ruling of X; ruling, no run; `--g-mom-damping 1.5` `--g-mom-damping 2.0` "
     "`--captain-ruling R22`) |", False, False),
    ("| R22 awaits (the captain's order of X; ruling, no run; `--g-mom-damping 1.5` `--g-mom-damping 2.0` "
     "`--captain-ruling R22`) |", False, False),  # the first cell must start with the order
    ("| R22 (the captain's order of X; ruling, no run) | T1 | `--g-mom-damping 1.5` `--g-mom-damping 2.0` "
     "`--captain-ruling R22` |", False, False),  # the options named outside the first cell
    ("| R22 (the captain's order of X; ruling, no run; `--g-mom-damping 1.5` `--g-mom-damping 2.0`) |",
     False, False),
    ("| R22 (the captain's order of X; ruling, no run; `--g-mom-damping 1.5` `--g-mom-damping 2.0` "
     "`--captain-ruling R20`) |", False, False),
    ("| R23 (the captain's order of X; ruling, no run; `--g-mom-damping 1.5` `--g-mom-damping 2.0` "
     "`--captain-ruling R22`) |", False, False),
])
def test_the_order_check_takes_only_a_single_ruling_row_naming_each_option(bt, tmp_path, text, ok15, ok20):
    p = tmp_path / "TRIES.md"
    p.write_text(f"# ledger\n{text}\n")
    assert bt._captain_order_names("R22", bt.T1_R22_DAMPING_OPTIONS[1.5], p) is ok15
    assert bt._captain_order_names("R22", bt.T1_R22_DAMPING_OPTIONS[2.0], p) is ok20
    p.write_text(f"# ledger\n{text}\n{text}\n")  # two rows with one ID: never
    assert not bt._captain_order_names("R22", bt.T1_R22_DAMPING_OPTIONS[1.5], p)
    assert not bt._captain_order_names("R22", bt.T1_R22_DAMPING_OPTIONS[2.0], p)


# ------------------------------------------------------------------ nothing existing changes


@pytest.mark.parametrize("name", list(BEFORE))
def test_the_existing_commands_give_the_same_bytes_and_records(bt, rc4, ledger, tmp_path, name):
    extra, file_sha, record_sha = BEFORE[name]
    f, rec = _build(bt, tmp_path, f"T1_rc{name}", extra)
    assert _sha(f) == file_sha and _record_sha(rec) == record_sha
    assert "captain_ruling" not in rec.get("t1_candidate", {})


def test_the_r20_file_and_the_first_lines_files_give_the_bytes_and_records_as_before(bt, rc4, ledger, tmp_path):
    f, rec = _build(bt, tmp_path, "T1_s2rcB2", R20)
    assert (_sha(f), _record_sha(rec)) == R20_BEFORE
    for v, name in ((2.5, "T1_s2rcB25"), (3.0, "T1_s2rcB30")):
        f, rec = _build(bt, tmp_path, name, LINE1[v])
        assert (_sha(f), _record_sha(rec)) == LINE1_BEFORE[v]


def test_candidate_a_needs_no_ruling_and_keeps_its_record_text(bt, rc4, ledger, tmp_path):
    ledger.write_text("# ledger\n")  # no ruling row at all: A builds all the same
    f, rec = _build(bt, tmp_path, "T1_rcA", ["--g-mom-damping", "1.0"])
    assert (_sha(f), _record_sha(rec)) == BEFORE["A"][1:]
    assert "candidate A, docs/PLAN.md decision log 2026-11-20 (i)" in rec["g_mom"]["rule"]
    assert "R22" not in rec["g_mom"]["rule"] and "captain_ruling" not in rec["t1_candidate"]


# ------------------------------------------------------------------ the line's two files


@pytest.mark.parametrize("v", [1.5, 2.0])
def test_the_line_file_is_the_default_with_only_the_damping_changed(bt, rc4, ledger, tmp_path, v):
    f0, r0 = _build(bt, tmp_path, "T1_rc", [])
    fv, rv = _build(bt, tmp_path, NAMES[v], DAMP[v])
    calendar = bt.feng_dirs.bracket_fraction(10.5, 12.5)
    assert rc4[0]["fraction"] == rc4[-1]["fraction"] == r0["fraction"] == rv["fraction"] == calendar == 0.4
    assert rc4[0]["groups"] == rc4[-1]["groups"] and rc4[0]["directions"] == rc4[-1]["directions"]
    b0, bv = rc4[0]["before"], rc4[-1]["before"]
    assert list(b0.obs_names) == list(bv.obs_names) and b0.obs.equals(bv.obs)
    cand = rv["t1_candidate"]
    assert cand["captain_ruling"] == "R22" and "row R22" in cand["order"] and "2026-11-22" in cand["order"]
    assert "second T1 search line" in cand["order"] and f"damping at {v:.1f}" in cand["order"]
    assert cand["g_mom_damping"] == v and cand["feng_fraction_mult"] == 1.0
    assert cand["applied_fraction"] == cand["calendar_fraction"] == 0.4 and cand["command"] == bt.T1_CANDIDATE_COMMAND
    assert set(cand) == {"g_mom_damping", "feng_fraction_mult", "captain_ruling", "order", "calendar_fraction",
                         "applied_fraction", "command"}  # the ruled files' fields, no other
    assert rv["g_mom"]["damping"] == v and "released_bet" not in rv
    assert (f"damped by {v:g} in place of B0_DAMPING (row R22's second T1 search line, docs/PLAN.md decision log "
            "2026-11-22)") in rv["g_mom"]["rule"]
    assert "candidate A" not in rv["g_mom"]["rule"]
    # the record is the default's but for the candidate field and g_mom's damping, and candidate A's likewise
    assert _strip(rv) == _strip(r0)
    fa, ra = _build(bt, tmp_path, "T1_rcA", ["--g-mom-damping", "1.0"])
    assert _strip(rv) == _strip(ra)
    assert ra["g_mom"]["rule"].replace("1 in place of B0_DAMPING (candidate A, docs/PLAN.md decision log 2026-11-20 "
                                       "(i))", f"{v:g} in place of B0_DAMPING (row R22's second T1 search line, "
                                       "docs/PLAN.md decision log 2026-11-22)") == rv["g_mom"]["rule"]
    assert _sha(fv) not in {*(b[1] for b in BEFORE.values()), R20_BEFORE[0], *(b[0] for b in LINE1_BEFORE.values())}


@pytest.mark.parametrize("v", [1.5, 2.0])
def test_the_line_file_changes_only_cells_g_mom_edits_and_builds_the_same_bytes_twice(bt, rc3, ledger, tmp_path, v):
    f0, _ = _build(bt, tmp_path, "T1_rc", [])
    fa, _ = _build(bt, tmp_path, f"{NAMES[v]}a", DAMP[v])
    fb, _ = _build(bt, tmp_path, f"{NAMES[v]}b", DAMP[v])
    assert fa.read_bytes() == fb.read_bytes()
    x0, xv = _same_frame(f0, fa)
    changed = _rows_changed(x0, xv)
    assert changed.any()
    pop = rc3[-1]["before"]
    labels = pop.obs["celltype"].astype(str).to_numpy()
    groups = rc3[-1]["groups"]
    with_dir = np.array([groups.get(lab) in THREE for lab in labels])
    assert not with_dir.all() and not (changed & ~with_dir).any()  # Epicardium and the rest keep their counts
    b0 = rc3[0]["before"]
    assert list(b0.obs_names) == list(pop.obs_names) and b0.obs.equals(pop.obs)


def test_the_two_values_give_two_files_and_2_is_2_0_and_an_explicit_multiplier_1_changes_nothing(bt, rc4, ledger,
                                                                                                  tmp_path):
    fa, _ = _build(bt, tmp_path, "T1_rcA", ["--g-mom-damping", "1.0"])
    f15, _ = _build(bt, tmp_path, NAMES[1.5], DAMP[1.5])
    f20, r20 = _build(bt, tmp_path, NAMES[2.0], DAMP[2.0])
    f2, r2 = _build(bt, tmp_path, "T1_s2rcA2", ["--g-mom-damping", "2", "--captain-ruling", "R22"])
    fm, rm = _build(bt, tmp_path, "T1_s2rcA20_m", [*DAMP[2.0], "--feng-fraction-mult", "1.0"])
    assert len({fa.read_bytes(), f15.read_bytes(), f20.read_bytes()}) == 3
    assert f2.read_bytes() == f20.read_bytes() == fm.read_bytes()
    assert r2["t1_candidate"] == r20["t1_candidate"] == rm["t1_candidate"]


# ------------------------------------------------------------------ refusals, before anything is read


@pytest.mark.parametrize("v", ["1.5", "2.0"])
@pytest.mark.parametrize("extra, match", [
    ([], "--g-mom-damping takes 0.5 .* needs --captain-ruling R22"),
    (["--captain-ruling", "R20"], "refused with --captain-ruling R20: only row R22"),
    (["--captain-ruling", "R21"], "refused with --captain-ruling R21: only row R22"),
    (["--captain-ruling", "R2"], "refused with --captain-ruling R2: only row R22"),
    (["--captain-ruling", "r22"], "refused with --captain-ruling r22: only row R22"),
    (["--captain-ruling", "T63"], "refused with --captain-ruling T63: only row R22"),
    (["--captain-ruling", "R22", "--feng-fraction-mult", "1.5"], "with --feng-fraction-mult other than 1"),
    (["--captain-ruling", "R22", "--feng-fraction-mult", "0.8"], "with --feng-fraction-mult other than 1"),
    (["--captain-ruling", "R22", "--feng-fraction-mult", "2.5"], "with --g-mom-damping other than 0.5"),
    (["--captain-ruling", "R22", "--feng-fraction-mult", "3.0"], "with --g-mom-damping other than 0.5"),
    (["--captain-ruling", "R22", "--feng-fraction-mult", "2.0"], "refused with --captain-ruling R22: only row R20"),
    (["--captain-ruling", "R20", "--feng-fraction-mult", "2.0"], "with --g-mom-damping other than 0.5"),
    (["--captain-ruling", "R22", "--bet-strength", "0.5"], "with --released-bet and --bet-strength"),
    (["--captain-ruling", "R22", "--released-bet", "share_trend", "--bet-strength", "0.5"],
     "with --released-bet and --bet-strength"),
    (["--captain-ruling", "R22", "--seed", "1000002"], "recomputed last file only"),
    (["--captain-ruling", "R22", "--board", "T1:val"], "recomputed last file only"),
])
def test_every_other_use_of_a_line_value_is_refused(bt, ledger, monkeypatch, tmp_path, v, extra, match):
    _refused(bt, monkeypatch, tmp_path, ["t1", *RC, "--g-mom-damping", v, *extra], match)


@pytest.mark.parametrize("extra, match", [
    (["--g-mom-damping", "1.0", "--captain-ruling", "R22"], "with --released-bet only"),  # A needs no ruling
    (["--g-mom-damping", "0.5", "--captain-ruling", "R22"], "with --released-bet only"),
    (["--g-mom-damping", "1.25", "--captain-ruling", "R22"], "--g-mom-damping takes 0.5"),  # no value between
    (["--g-mom-damping", "1.75", "--captain-ruling", "R22"], "--g-mom-damping takes 0.5"),
    (["--g-mom-damping", "2.5", "--captain-ruling", "R22"], "--g-mom-damping takes 0.5"),  # none above 2.0
    (["--g-mom-damping", "3", "--captain-ruling", "R22"], "--g-mom-damping takes 0.5"),
    (["--g-mom-damping", "0.75", "--captain-ruling", "R22"], "--g-mom-damping takes 0.5"),
    (["--g-mom-damping", "0"], "--g-mom-damping takes 0.5"),
    (["--g-mom-damping", "-1.5", "--captain-ruling", "R22"], "--g-mom-damping takes 0.5"),
])
def test_the_row_takes_no_other_damping(bt, ledger, monkeypatch, tmp_path, extra, match):
    _refused(bt, monkeypatch, tmp_path, ["t1", *RC, *extra], match)


@pytest.mark.parametrize("v", ["1.5", "2.0"])
@pytest.mark.parametrize("argv", [
    ["t1"],  # the authors' groups
    ["t1", "--g-mom", "--base", "lib3"],  # the first upload's command
    ["t1", "--g-mom", "--base", "lib3", "--feng-groups", "none"],  # the none file
    ["t1", "--g-mom", "--base", "mid", "--cm-share", "0.60"],  # T1's second upload
])
def test_a_line_value_is_refused_on_any_other_command(bt, ledger, monkeypatch, tmp_path, v, argv):
    _refused(bt, monkeypatch, tmp_path, [*argv, "--g-mom-damping", v, "--captain-ruling", "R22"],
             "recomputed last file only")


@pytest.mark.parametrize("extra", [["--cm-share", "0.6"], ["--herd"], ["--no-feng-hygiene"], ["--lib3-pair"],
                                   ["--base", "pooled"]])
def test_the_recomputed_commands_own_guards_hold_with_a_line_value(bt, ledger, monkeypatch, tmp_path, extra):
    for v in ("1.5", "2.0"):
        _refused(bt, monkeypatch, tmp_path, ["t1", *RC, *extra, "--g-mom-damping", v, "--captain-ruling", "R22"],
                 "--feng-groups recomputed builds one file")


def test_a_line_value_is_refused_when_the_row_does_not_order_it(bt, ledger, monkeypatch, tmp_path):
    # row R22 as it stands before the sentence is added: the first line's options only
    ledger.write_text(f"# ledger\n{ROW_R20}\n| R22 (the captain's order of X; ruling, no run; "
                      "`--feng-fraction-mult 2.5` `--feng-fraction-mult 3.0` `--captain-ruling R22`) | T1 |\n")
    for v in (1.5, 2.0):
        _refused(bt, monkeypatch, tmp_path, ["t1", *RC, *DAMP[v]], "no single row R22 .* `--g-mom-damping")
    ledger.write_text(f"# ledger\n{ROW_R20}\n| R22 (awaits the captain's order; `--g-mom-damping 1.5` "
                      "`--g-mom-damping 2.0` `--captain-ruling R22`) | T1 |\n")
    _refused(bt, monkeypatch, tmp_path, ["t1", *RC, *DAMP[1.5]], "no single row R22")
    ledger.write_text(f"# ledger\n{ROW_R20}\n")  # no row R22
    _refused(bt, monkeypatch, tmp_path, ["t1", *RC, *DAMP[2.0]], "no single row R22")
    # a row naming 1.5 alone orders 1.5 alone
    ledger.write_text(f"# ledger\n{ROW_R20}\n| R22 (the captain's order of X; ruling, no run; "
                      "`--g-mom-damping 1.5` with `--captain-ruling R22`) | T1 |\n")
    _refused(bt, monkeypatch, tmp_path, ["t1", *RC, *DAMP[2.0]], "no single row R22 .* `--g-mom-damping 2.0`")
    a = bt.build_parser().parse_args(["t1", *RC, *DAMP[1.5], "--out", str(tmp_path / "T1.h5ad")])
    assert bt._t1_candidate(a)["captain_ruling"] == "R22"


def test_the_bet_refusal_lets_only_the_ruled_commands_through(bt, ledger):
    def ns(**k):
        base = dict(released_bet=None, bet_strength=None, captain_ruling="R22", feng_fraction_mult=None,
                    feng_groups="recomputed", seed=harness.BUILD_SEED_OFFSET + 1, board=None, g_mom_damping=1.5)
        return SimpleNamespace(**{**base, **k})

    for good in (ns(), ns(g_mom_damping=2.0), ns(feng_fraction_mult=1.0),
                 ns(g_mom_damping=None, feng_fraction_mult=2.5),
                 ns(g_mom_damping=None, feng_fraction_mult=2.0, captain_ruling="R20")):
        assert bt._bet_refusals(good) is None
    for bad in (ns(captain_ruling="R20"), ns(g_mom_damping=1.0), ns(g_mom_damping=0.5), ns(g_mom_damping=1.75),
                ns(feng_fraction_mult=1.5), ns(captain_ruling="R21"), ns(g_mom_damping=None)):
        with pytest.raises(SystemExit, match="with --released-bet only"):
            bt._bet_refusals(bad)
    ledger.write_text(f"# ledger\n{ROW_R20}\n")  # without row R22 the shortcut is closed
    with pytest.raises(SystemExit, match="with --released-bet only"):
        bt._bet_refusals(ns())
