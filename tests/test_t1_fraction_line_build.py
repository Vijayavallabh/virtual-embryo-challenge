"""`build_test.py t1 --feng-fraction-mult 2.5|3.0 --captain-ruling R22`: the first search line of row R22 of
experiments/harness/TRIES.md (the captain's order of 2026-11-22, docs/PLAN.md decision log 2026-11-22 (j), the
amendment of docs/S2_RULES.md "a score-guided search of round steps may run"): the recomputed-groups file (row R13) with
the outside atlas's direction step at 2.5 or 3.0 times the calendar fraction and nothing else changed.

Held here, on synthetic stages and directions only (the build_t1 fixture of tests/test_hidden_counts.py and the
recomputed setup of tests/test_feng_relabel.py), with a temporary ledger in place of the repository's: (1) nothing that
builds today changes: the default command and candidates A, B and C give the bytes and records pinned in
tests/test_t1_fraction2_build.py, and row R20's 2.0 file gives the bytes and record the builder at main a431e10 (before
this option existed) gave on this fixture, pinned below; (2) at 2.5 and 3.0 the only change against the default is the
step's size; (3) each value is taken only with --captain-ruling R22 while row R22 is the captain's order naming that
value's option, and every other use is refused before anything is read. No released file, outside-atlas file, board or
score is read."""

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

from vec import harness  # noqa: E402

bt, rc3, rc4 = tcb.bt, tcb.rc3, tcb.rc4  # the candidates' fixtures (the builder module; the recomputed setups)

R22 = {2.5: ["--feng-fraction-mult", "2.5", "--captain-ruling", "R22"],
       3.0: ["--feng-fraction-mult", "3.0", "--captain-ruling", "R22"]}
NAMES = {2.5: "T1_s2rcB25", 3.0: "T1_s2rcB30"}
# Row R20's file on this fixture, built by the builder at main a431e10 (before the R22 values existed), twice in
# separate processes with the same result: the file's sha256 and its record's (feng_hygiene.commit set to "own_repo").
R20_BEFORE = ("c2d8d6a2916f5442e18a8ba9fef44f26ddc64c6c1f4a3656deacea7788b3d7d2",
              "20c3e74da7a396c29498cbb2fd3c7a07f0abd20f773cde8023d8757e24c9d792")
# Synthetic ruling rows in the form the builder checks (first cell: the captain's order, "ruling, no run", the options)
ROW_R20 = ("| R20 (the captain's order of a day; ruling, no run, not a try; the builder is to take "
           "`--feng-fraction-mult 2.0` only with `--captain-ruling R20`) | T1 | a file |")
ROW_R22 = ("| R22 (the captain's order of a later day: a search line; ruling, no run, not a try; the builder is to "
           "take `--feng-fraction-mult 2.5` and `--feng-fraction-mult 3.0` only with `--captain-ruling R22`) | T1 | "
           "a line |")


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
    """A record without what the step's size sets: the candidate field, the fraction and g_mom's explicit damping."""
    out = {k: v for k, v in r.items() if k not in ("t1_candidate", "fraction", "realised")}
    out["g_mom"] = {k: v for k, v in r["g_mom"].items() if k != "damping"}
    return out


# ------------------------------------------------------------------ the constants and the ruling check


def test_the_line_has_its_own_constants_and_r20_is_unchanged(bt):
    assert bt.FENG_FRACTION_MULTS == (1.0, 1.5) and bt.G_MOM_DAMPINGS == (0.5, 1.0)
    assert bt.FENG_FRACTION_MULTS_R22 == (2.5, 3.0) and bt.T1_R22_ROW == "R22"
    assert bt.T1_R22_OPTIONS == {2.5: "`--feng-fraction-mult 2.5`", 3.0: "`--feng-fraction-mult 3.0`"}
    assert bt.FENG_FRACTION_MULT_R20 == 2.0 and bt.T1_R20_ROW == "R20"
    assert bt.T1_R20_OPTION == "`--feng-fraction-mult 2.0`"
    ruled = bt.FENG_FRACTION_MULTS_RULED
    assert set(ruled) == {2.0, 2.5, 3.0} and not set(ruled) & set(bt.FENG_FRACTION_MULTS)
    assert ruled[2.0][:2] == ("R20", "`--feng-fraction-mult 2.0`")
    assert ruled[2.0][2] == ("experiments/harness/TRIES.md row R20, the captain's order of 2026-11-21 (docs/PLAN.md "
                             "decision log 2026-11-21 (x)): the direction at 2.0 x the calendar fraction, nothing else "
                             "changed")
    for v in (2.5, 3.0):
        row, option, order = ruled[v]
        assert row == "R22" and option == bt.T1_R22_OPTIONS[v]
        assert "row R22" in order and "2026-11-22 (j)" in order and f"at {v:.1f} x the calendar fraction" in order


def test_the_order_check_reads_each_option_from_row_r22s_first_cell(bt, ledger):
    for v in (2.5, 3.0):
        assert bt._captain_order_names("R22", bt.T1_R22_OPTIONS[v], ledger)
        assert not bt._captain_order_names("R20", bt.T1_R22_OPTIONS[v], ledger)
    assert bt._captain_order_names("R20", bt.T1_R20_OPTION, ledger)
    assert not bt._captain_order_names("R22", bt.T1_R20_OPTION, ledger)
    assert not bt._captain_ruling_names("R22", "share_trend_l0.5", ledger)  # the bets' check takes no R22


@pytest.mark.parametrize("text, ok25, ok30", [
    (ROW_R22, True, True),
    ("| R22 (the captain's order of X; ruling, no run; `--feng-fraction-mult 2.5` with `--captain-ruling R22`) |",
     True, False),  # names only 2.5: 3.0 is not ordered
    ("| R22 (the captain's order of X; ruling, no run; `--feng-fraction-mult 3` with `--captain-ruling R22`) |",
     False, False),  # the option as the row must write it: 3.0
    ("| R22 (the captain's order of X; a try; `--feng-fraction-mult 2.5` `--feng-fraction-mult 3.0` "
     "`--captain-ruling R22`) |", False, False),
    ("| R22 (the captain's ruling of X; ruling, no run; `--feng-fraction-mult 2.5` `--feng-fraction-mult 3.0` "
     "`--captain-ruling R22`) |", False, False),
    ("| R22 (the captain's order of X; ruling, no run) | T1 | `--feng-fraction-mult 2.5` `--feng-fraction-mult 3.0` "
     "`--captain-ruling R22` |", False, False),  # the options named outside the first cell
    ("| R22 (the captain's order of X; ruling, no run; `--feng-fraction-mult 2.5` `--feng-fraction-mult 3.0`) |",
     False, False),
    ("| R22 (the captain's order of X; ruling, no run; `--feng-fraction-mult 2.5` `--feng-fraction-mult 3.0` "
     "`--captain-ruling R20`) |", False, False),
    ("| R23 (the captain's order of X; ruling, no run; `--feng-fraction-mult 2.5` `--feng-fraction-mult 3.0` "
     "`--captain-ruling R22`) |", False, False),
])
def test_the_order_check_takes_only_a_single_ruling_row_naming_each_option(bt, tmp_path, text, ok25, ok30):
    p = tmp_path / "TRIES.md"
    p.write_text(f"# ledger\n{text}\n")
    assert bt._captain_order_names("R22", bt.T1_R22_OPTIONS[2.5], p) is ok25
    assert bt._captain_order_names("R22", bt.T1_R22_OPTIONS[3.0], p) is ok30
    p.write_text(f"# ledger\n{text}\n{text}\n")  # two rows with one ID: never
    assert not bt._captain_order_names("R22", bt.T1_R22_OPTIONS[2.5], p)
    assert not bt._captain_order_names("R22", bt.T1_R22_OPTIONS[3.0], p)


# ------------------------------------------------------------------ nothing existing changes


@pytest.mark.parametrize("name", list(BEFORE))
def test_the_existing_commands_give_the_same_bytes_and_records(bt, rc4, ledger, tmp_path, name):
    extra, file_sha, record_sha = BEFORE[name]
    f, rec = _build(bt, tmp_path, f"T1_rc{name}", extra)
    assert _sha(f) == file_sha and _record_sha(rec) == record_sha


def test_the_r20_file_gives_the_bytes_and_record_main_gave(bt, rc4, ledger, tmp_path):
    f, rec = _build(bt, tmp_path, "T1_s2rcB2", R20)
    assert (_sha(f), _record_sha(rec)) == R20_BEFORE
    assert rec["t1_candidate"]["captain_ruling"] == "R20" and rec["t1_candidate"]["feng_fraction_mult"] == 2.0


# ------------------------------------------------------------------ the line's two files


@pytest.mark.parametrize("v", [2.5, 3.0])
def test_the_line_file_is_the_default_with_only_the_steps_size_changed(bt, rc4, ledger, tmp_path, v):
    f0, r0 = _build(bt, tmp_path, "T1_rc", [])
    fv, rv = _build(bt, tmp_path, NAMES[v], R22[v])
    calendar = bt.feng_dirs.bracket_fraction(10.5, 12.5)
    assert rc4[0]["fraction"] == r0["fraction"] == calendar == 0.4
    assert rc4[-1]["fraction"] == v * calendar == pytest.approx(v * 0.4, abs=1e-12)
    # what build_submission is given besides the step's size is the default's
    assert rc4[0]["groups"] == rc4[-1]["groups"] and rc4[0]["directions"] == rc4[-1]["directions"]
    b0, bv = rc4[0]["before"], rc4[-1]["before"]
    assert list(b0.obs_names) == list(bv.obs_names) and b0.obs.equals(bv.obs)
    cand = rv["t1_candidate"]
    assert cand["captain_ruling"] == "R22" and "row R22" in cand["order"] and "2026-11-22 (j)" in cand["order"]
    assert cand["feng_fraction_mult"] == v and cand["g_mom_damping"] == 0.5
    assert cand["applied_fraction"] == rv["fraction"] == rc4[-1]["fraction"]
    assert cand["calendar_fraction"] == 0.4 and cand["command"] == bt.T1_CANDIDATE_COMMAND
    assert set(cand) == {"g_mom_damping", "feng_fraction_mult", "captain_ruling", "order", "calendar_fraction",
                         "applied_fraction", "command"}  # the R20 file's fields, no other
    assert rv["g_mom"]["damping"] == 0.5 and "released_bet" not in rv
    # the record is the default's but for the step's size and the candidate field
    assert _strip(rv) == _strip(r0)
    _, rb = _build(bt, tmp_path, "T1_rcB", ["--feng-fraction-mult", "1.5"])
    assert _strip(rv) == _strip(rb)
    assert _sha(fv) not in {*(b[1] for b in BEFORE.values()), R20_BEFORE[0]}


@pytest.mark.parametrize("v", [2.5, 3.0])
def test_the_line_file_changes_only_cells_with_a_direction_and_builds_the_same_bytes_twice(bt, rc3, ledger, tmp_path,
                                                                                           v):
    f0, _ = _build(bt, tmp_path, "T1_rc", [])
    fa, _ = _build(bt, tmp_path, f"{NAMES[v]}a", R22[v])
    fb, _ = _build(bt, tmp_path, f"{NAMES[v]}b", R22[v])
    assert fa.read_bytes() == fb.read_bytes()
    x0, xv = _same_frame(f0, fa)
    changed = _rows_changed(x0, xv)
    assert changed.any()
    pop = rc3[-1]["before"]
    labels = pop.obs["celltype"].astype(str).to_numpy()
    groups = rc3[-1]["groups"]
    with_dir = np.array([groups.get(lab) in THREE for lab in labels])
    assert not with_dir.all() and not (changed & ~with_dir).any()  # Epicardium and the rest keep their counts
    b0 = rc3[0]["before"]  # the cells, their order and obs are the default's
    assert list(b0.obs_names) == list(pop.obs_names) and b0.obs.equals(pop.obs)


def test_the_two_values_give_two_files_and_3_is_3_0(bt, rc4, ledger, tmp_path):
    f25, _ = _build(bt, tmp_path, NAMES[2.5], R22[2.5])
    f30, r30 = _build(bt, tmp_path, NAMES[3.0], R22[3.0])
    f3, r3 = _build(bt, tmp_path, "T1_s2rcB3", ["--feng-fraction-mult", "3", "--captain-ruling", "R22"])
    fd, rd = _build(bt, tmp_path, "T1_s2rcB30_d", [*R22[3.0], "--g-mom-damping", "0.5"])
    assert f25.read_bytes() != f30.read_bytes()
    assert f3.read_bytes() == f30.read_bytes() == fd.read_bytes()
    assert r3["t1_candidate"] == r30["t1_candidate"] == rd["t1_candidate"]


# ------------------------------------------------------------------ refusals, before anything is read


@pytest.mark.parametrize("v", ["2.5", "3.0"])
@pytest.mark.parametrize("extra, match", [
    ([], "--feng-fraction-mult takes 1 .* needs --captain-ruling R22"),
    (["--captain-ruling", "R20"], "refused with --captain-ruling R20: only row R22"),
    (["--captain-ruling", "R21"], "refused with --captain-ruling R21: only row R22"),
    (["--captain-ruling", "R2"], "refused with --captain-ruling R2: only row R22"),
    (["--captain-ruling", "r22"], "refused with --captain-ruling r22: only row R22"),
    (["--captain-ruling", "T63"], "refused with --captain-ruling T63: only row R22"),
    # damping 1.0 with 2.5 is the third line's pair, refused here only for want of its order in this ledger (in
    # tests/test_t1_line3_build.py), and 3.0 with damping 1.0 is the seventh line's, 2.5 with 1.5 or 2.0 the sixth
    # line's (tests/test_t1_l67pb_build.py; 1.5 below)
    (["--captain-ruling", "R22", "--g-mom-damping", "0.75"], "with --g-mom-damping other than 0.5"),
    (["--captain-ruling", "R22", "--bet-strength", "0.5"], "with --released-bet and --bet-strength"),
    (["--captain-ruling", "R22", "--released-bet", "share_trend", "--bet-strength", "0.5"],
     "with --released-bet and --bet-strength"),
    (["--captain-ruling", "R22", "--seed", "1000002"], "recomputed last file only"),
    (["--captain-ruling", "R22", "--board", "T1:val"], "recomputed last file only"),
])
def test_every_other_use_of_a_line_value_is_refused(bt, ledger, monkeypatch, tmp_path, v, extra, match):
    _refused(bt, monkeypatch, tmp_path, ["t1", *RC, "--feng-fraction-mult", v, *extra], match)


def test_a_line_value_with_the_second_lines_damping_is_refused(bt, ledger, monkeypatch, tmp_path):
    """Damping 1.5 with 3.0: refused as before; with 2.5 it is the sixth line's pair (tests/test_t1_l67pb_build.py),
    refused here for want of its order in this ledger."""
    _refused(bt, monkeypatch, tmp_path, ["t1", *RC, *R22[3.0], "--g-mom-damping", "1.5"],
             "with --g-mom-damping other than 0.5")
    _refused(bt, monkeypatch, tmp_path, ["t1", *RC, *R22[2.5], "--g-mom-damping", "1.5"],
             "no single row R22 .* `--g-mom-damping 1.5 --feng-fraction-mult 2.5`")


@pytest.mark.parametrize("extra, match", [
    (["--feng-fraction-mult", "2.0", "--captain-ruling", "R22"], "refused with --captain-ruling R22: only row R20"),
    (["--feng-fraction-mult", "1.5", "--captain-ruling", "R22"], "with --released-bet only"),
    (["--feng-fraction-mult", "1.0", "--captain-ruling", "R22"], "with --released-bet only"),
    (["--captain-ruling", "R22"], "with --released-bet only"),
    (["--g-mom-damping", "1.0", "--captain-ruling", "R22"], "with --released-bet only"),
    (["--feng-fraction-mult", "2.25", "--captain-ruling", "R22"], "--feng-fraction-mult takes 1"),  # no value between
    (["--feng-fraction-mult", "2.75", "--captain-ruling", "R22"], "--feng-fraction-mult takes 1"),
    (["--feng-fraction-mult", "3.5", "--captain-ruling", "R22"], "--feng-fraction-mult takes 1"),  # the line ends at 3
    (["--feng-fraction-mult", "4", "--captain-ruling", "R22"], "--feng-fraction-mult takes 1"),
])
def test_the_row_takes_no_other_value(bt, ledger, monkeypatch, tmp_path, extra, match):
    _refused(bt, monkeypatch, tmp_path, ["t1", *RC, *extra], match)


@pytest.mark.parametrize("v", ["2.5", "3.0"])
@pytest.mark.parametrize("argv", [
    ["t1"],  # the authors' groups
    ["t1", "--g-mom", "--base", "lib3"],  # the first upload's command
    ["t1", "--g-mom", "--base", "lib3", "--feng-groups", "none"],  # the none file
    ["t1", "--g-mom", "--base", "mid", "--cm-share", "0.60"],  # T1's second upload
])
def test_a_line_value_is_refused_on_any_other_command(bt, ledger, monkeypatch, tmp_path, v, argv):
    _refused(bt, monkeypatch, tmp_path, [*argv, "--feng-fraction-mult", v, "--captain-ruling", "R22"],
             "recomputed last file only")


@pytest.mark.parametrize("extra", [["--cm-share", "0.6"], ["--herd"], ["--no-feng-hygiene"], ["--lib3-pair"],
                                   ["--base", "pooled"]])
def test_the_recomputed_commands_own_guards_hold_with_a_line_value(bt, ledger, monkeypatch, tmp_path, extra):
    for v in ("2.5", "3.0"):
        _refused(bt, monkeypatch, tmp_path, ["t1", *RC, *extra, "--feng-fraction-mult", v, "--captain-ruling", "R22"],
                 "--feng-groups recomputed builds one file")


def test_a_line_value_is_refused_when_the_row_does_not_order_it(bt, ledger, monkeypatch, tmp_path):
    ledger.write_text(f"# ledger\n{ROW_R20}\n| R22 (awaits the captain's order; `--feng-fraction-mult 2.5` "
                      "`--feng-fraction-mult 3.0` `--captain-ruling R22`) | T1 |\n")
    for v in ("2.5", "3.0"):
        _refused(bt, monkeypatch, tmp_path, ["t1", *RC, *R22[float(v)]], "no single row R22")
    ledger.write_text(f"# ledger\n{ROW_R20}\n")  # no row R22
    _refused(bt, monkeypatch, tmp_path, ["t1", *RC, *R22[2.5]], "no single row R22")
    # a row naming 2.5 alone orders 2.5 alone
    ledger.write_text(f"# ledger\n{ROW_R20}\n| R22 (the captain's order of X; ruling, no run; "
                      "`--feng-fraction-mult 2.5` with `--captain-ruling R22`) | T1 |\n")
    _refused(bt, monkeypatch, tmp_path, ["t1", *RC, *R22[3.0]], "no single row R22 .* `--feng-fraction-mult 3.0`")
    a = bt.build_parser().parse_args(["t1", *RC, *R22[2.5], "--out", str(tmp_path / "T1.h5ad")])
    assert bt._t1_candidate(a)["captain_ruling"] == "R22"


def test_the_bet_refusal_lets_only_the_ruled_commands_through(bt, ledger):
    def ns(**k):
        base = dict(released_bet=None, bet_strength=None, captain_ruling="R22", feng_fraction_mult=2.5,
                    feng_groups="recomputed", seed=harness.BUILD_SEED_OFFSET + 1, board=None, g_mom_damping=None)
        return SimpleNamespace(**{**base, **k})

    for good in (ns(), ns(feng_fraction_mult=3.0), ns(feng_fraction_mult=2.0, captain_ruling="R20")):
        assert bt._bet_refusals(good) is None
    for bad in (ns(captain_ruling="R20"), ns(feng_fraction_mult=2.0), ns(feng_fraction_mult=1.5),
                ns(feng_fraction_mult=None), ns(feng_fraction_mult=2.75), ns(captain_ruling="R21"),
                ns(feng_fraction_mult=3.0, captain_ruling="R20")):
        with pytest.raises(SystemExit, match="with --released-bet only"):
            bt._bet_refusals(bad)
    ledger.write_text(f"# ledger\n{ROW_R20}\n")  # without row R22 the shortcut is closed
    with pytest.raises(SystemExit, match="with --released-bet only"):
        bt._bet_refusals(ns())
