"""`build_test.py t1 --feng-fraction-mult 2.0 --captain-ruling R20`: row R20 of experiments/harness/TRIES.md (the
captain's order of 2026-11-21, docs/PLAN.md decision log 2026-11-21 (x), the last amendment of docs/S2_RULES.md): the
recomputed-groups file (row R13) with the outside atlas's direction step at 2.0 times the calendar fraction and nothing
else changed.

Two things are held here. (1) Nothing that builds today changes: the default command (R13's file), candidate A
(--g-mom-damping 1.0), B (--feng-fraction-mult 1.5) and C (both) give, on the synthetic fixture, the file bytes and the
records that the builder at main f7bd2f0 (before this option existed) gave, pinned below; FENG_FRACTION_MULTS stays
(1.0, 1.5). (2) The value 2.0 is taken only on row R20's one command, and every other use is refused before anything is
read. Synthetic stages and directions only (the build_t1 fixture of tests/test_hidden_counts.py and the recomputed setup
of tests/test_feng_relabel.py): no released file, outside-atlas file, board or score is read."""

from __future__ import annotations

import hashlib
import json
import re
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
for _p in ("experiments/external", "experiments/t1_m1", "experiments/t1_direction", "experiments/t1_hygiene",
           "experiments/builds"):
    sys.path.insert(0, str(ROOT / _p))

import test_t1_candidates_build as tcb  # noqa: E402
import test_t1_feng_groups_none  # noqa: E402, F401  (imported before any fixture, as tests/test_feng_relabel.py does)
from test_hidden_counts import t1_env  # noqa: E402, F401  (the synthetic build_t1 fixture)
from test_t1_candidates_build import RC, THREE, TODAY_RC_SHA, _build, _no_read, _rows_changed, _same_frame  # noqa: E402

from vec import harness  # noqa: E402

bt, rc3, rc4 = tcb.bt, tcb.rc3, tcb.rc4  # the candidates' fixtures (the builder module; the recomputed setups)

R20 = ["--feng-fraction-mult", "2.0", "--captain-ruling", "R20"]
# The synthetic fixture's four existing commands, built by the builder at main f7bd2f0 (before this option existed),
# twice in separate processes with the same result: each file's sha256 and the sha256 of its record dumped with sorted
# keys, the record's feng_hygiene.commit (the hash of the test's own throwaway repository, new on every run) set to
# "own_repo". The default's file hash is TODAY_RC_SHA of tests/test_t1_candidates_build.py.
BEFORE = {
    "D": ([], "a5899e2d0a4b97f8fcabca6f8781a613afdd070194c0a3c76e0b14d92249d132",
          "c6428f37fa92008735ceab9b6b29d806eb41544cc2d9c70eb30d92fe8290ff8a"),
    "A": (["--g-mom-damping", "1.0"], "5191a9f62d07946f1cd7d566e7c450ea1e0a607c0b8921c44300b8e9e27f9d5e",
          "fe29b18413da185bbe72edc03a4a1823be08cf164b7b448b16b52966eba042a2"),
    "B": (["--feng-fraction-mult", "1.5"], "36fc0355435d581ac6865582e563b82df0a2e1171042ef2f8553071e63a71a46",
          "e4a7cb8773e1568fe8ab73b199125faed91710c23cd73bbb4fd75cf60eb93b0e"),
    "C": (["--g-mom-damping", "1.0", "--feng-fraction-mult", "1.5"],
          "f22325c52725a9db33790d768b3c6af012b6ccef3f9e8d0a6a364a0a9531a4da",
          "810a4827528c71c0dd3f43e85f343949b40e63c665a5de092e8a95b9b5486e4f"),
}


def _sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def _record_sha(rec: dict) -> str:
    r = json.loads(json.dumps(rec))
    r["feng_hygiene"]["commit"] = "own_repo"
    return hashlib.sha256(json.dumps(r, sort_keys=True).encode()).hexdigest()


# ------------------------------------------------------------------ the constants and the ruling row


def test_the_existing_values_stay_and_the_ruled_value_has_its_own_constant(bt):
    assert bt.FENG_FRACTION_MULTS == (1.0, 1.5) and bt.G_MOM_DAMPINGS == (0.5, 1.0)
    assert bt.FENG_FRACTION_MULT_R20 == 2.0 and bt.FENG_FRACTION_MULT_R20 not in bt.FENG_FRACTION_MULTS
    assert bt.T1_R20_ROW == "R20" and bt.T1_R20_OPTION == "`--feng-fraction-mult 2.0`"
    a = bt.build_parser().parse_args(["t1", *RC, "--out", "x.h5ad"])
    assert a.feng_fraction_mult is None and a.captain_ruling is None


def test_row_r20_of_the_ledger_is_the_captains_order_naming_the_option(bt):
    assert bt._captain_order_names("R20", bt.T1_R20_OPTION)
    assert not bt._captain_order_names("R19", bt.T1_R20_OPTION)  # another ruling row
    assert not bt._captain_order_names("R20", "`--feng-fraction-mult 2.5`")
    # the old ruling check, for --released-bet, does not take R20 for any arm
    assert not bt._captain_ruling_names("R20", "share_trend_l0.5")
    rules = " ".join((ROOT / "docs" / "S2_RULES.md").read_text().split())
    assert ("build_test.py t1 --g-mom --base lib3 --feng-groups recomputed --feng-fraction-mult 2.0 --captain-ruling "
            "R20 --out runs/builds/T1_s2rcB2.h5ad") in rules


@pytest.mark.parametrize("text, ok", [
    ("| R20 (the captain's order of X; ruling, no run, not a try; `--feng-fraction-mult 2.0` with "
     "`--captain-ruling R20`) | T1 |", True),
    ("| R20 (the captain's order of X; a try; `--feng-fraction-mult 2.0` with `--captain-ruling R20`) | T1 |", False),
    ("| R20 (the captain's ruling of X; ruling, no run; `--feng-fraction-mult 2.0` with `--captain-ruling R20`) |",
     False),
    ("| R20 (the captain's order of X; ruling, no run) | T1 | `--feng-fraction-mult 2.0` `--captain-ruling R20` |",
     False),  # the option named outside the first cell
    ("| R20 (the captain's order of X; ruling, no run; `--feng-fraction-mult 2.0`) | T1 |", False),
    ("| R21 (the captain's order of X; ruling, no run; `--feng-fraction-mult 2.0` with `--captain-ruling R20`) |",
     False),
])
def test_the_order_check_takes_only_a_single_ruling_row_naming_the_option(bt, tmp_path, text, ok):
    p = tmp_path / "TRIES.md"
    p.write_text(f"# ledger\n{text}\n")
    assert bt._captain_order_names("R20", bt.T1_R20_OPTION, p) is ok
    p.write_text(f"# ledger\n{text}\n{text}\n")  # two rows with one ID: never
    assert bt._captain_order_names("R20", bt.T1_R20_OPTION, p) is False


# ------------------------------------------------------------------ nothing existing changes


@pytest.mark.parametrize("name", list(BEFORE))
def test_the_existing_commands_give_mains_bytes_and_records(bt, rc4, tmp_path, name):
    extra, file_sha, record_sha = BEFORE[name]
    f, rec = _build(bt, tmp_path, f"T1_rc{name}", extra)
    assert _sha(f) == file_sha
    assert _record_sha(rec) == record_sha
    assert "captain_ruling" not in rec.get("t1_candidate", {}) and "released_bet" not in rec
    if name == "D":
        assert file_sha == TODAY_RC_SHA and "t1_candidate" not in rec


# ------------------------------------------------------------------ row R20's file


def test_the_r20_file_is_the_direction_at_0_8_and_its_record_names_the_row(bt, rc4, tmp_path):
    _, r0 = _build(bt, tmp_path, "T1_rc", [])
    f2, r2 = _build(bt, tmp_path, "T1_s2rcB2", R20)
    assert [e["fraction"] for e in rc4][0] == 0.4
    assert rc4[-1]["fraction"] == pytest.approx(0.8, abs=1e-12)
    assert rc4[-1]["fraction"] == 2.0 * bt.feng_dirs.bracket_fraction(10.5, 12.5)
    cand = r2["t1_candidate"]
    assert cand["captain_ruling"] == "R20" and "row R20" in cand["order"] and "2026-11-21 (x)" in cand["order"]
    assert cand["feng_fraction_mult"] == 2.0 and cand["g_mom_damping"] == 0.5
    assert cand["applied_fraction"] == r2["fraction"] == rc4[-1]["fraction"]
    assert cand["calendar_fraction"] == r0["fraction"] == 0.4
    assert cand["command"] == bt.T1_CANDIDATE_COMMAND
    assert r2["g_mom"]["damping"] == 0.5 and "released_bet" not in r2
    # everything but the fraction and the candidate field is B's record
    _, rb = _build(bt, tmp_path, "T1_rcB", ["--feng-fraction-mult", "1.5"])
    strip = lambda r: {k: v for k, v in r.items() if k not in ("t1_candidate", "fraction", "realised")}  # noqa: E731
    assert strip(r2) == strip(rb)
    assert _sha(f2) not in {v[1] for v in BEFORE.values()}


def test_the_r20_file_changes_only_cells_with_a_direction_and_builds_the_same_bytes_twice(bt, rc3, tmp_path):
    f0, _ = _build(bt, tmp_path, "T1_rc", [])
    fa, _ = _build(bt, tmp_path, "T1_s2rcB2a", R20)
    fb, _ = _build(bt, tmp_path, "T1_s2rcB2b", R20)
    assert fa.read_bytes() == fb.read_bytes()
    x0, x2 = _same_frame(f0, fa)
    changed = _rows_changed(x0, x2)
    assert changed.any()
    pop = rc3[-1]["before"]
    labels = pop.obs["celltype"].astype(str).to_numpy()
    groups = rc3[-1]["groups"]
    with_dir = np.array([groups.get(lab) in THREE for lab in labels])
    assert not with_dir.all() and not (changed & ~with_dir).any()
    b0 = rc3[0]["before"]  # the cells, their order and obs are the default's
    assert list(b0.obs_names) == list(pop.obs_names) and b0.obs.equals(pop.obs)


def test_the_default_damping_given_explicitly_builds_the_same_r20_file(bt, rc4, tmp_path):
    fa, ra = _build(bt, tmp_path, "T1_s2rcB2", R20)
    fb, rb = _build(bt, tmp_path, "T1_s2rcB2_d", [*R20, "--g-mom-damping", "0.5"])
    assert fa.read_bytes() == fb.read_bytes() and ra["t1_candidate"] == rb["t1_candidate"]


# ------------------------------------------------------------------ refusals, before anything is read


NUMBERS_ALLOWED = {"0.5", "1", "1.5", "2", "3", "20"}  # the option's values, lib3 and the row's ID: no data


def _refused(bt, monkeypatch, tmp_path, argv: list[str], match: str) -> str:
    monkeypatch.setattr(bt, "load_stage", _no_read)
    monkeypatch.setattr(bt.feng_dirs, "feng_directions", _no_read)
    name = "T1_dry.h5ad" if "--board" in argv else "T1.h5ad"
    with pytest.raises(SystemExit, match=match) as e:
        bt.main([*argv, "--out", str(tmp_path / name)])
    msg = str(e.value)
    given = set(re.findall(r"\d+(?:\.\d+)?", " ".join(argv)))  # what the command line itself said
    assert set(re.findall(r"\d+(?:\.\d+)?", msg)) <= NUMBERS_ALLOWED | given, msg
    return msg


@pytest.mark.parametrize("extra, match", [
    (["--feng-fraction-mult", "2.0"], "--feng-fraction-mult takes 1 .* needs --captain-ruling R20"),
    (["--feng-fraction-mult", "2.0", "--captain-ruling", "R19"], "refused with --captain-ruling R19"),
    (["--feng-fraction-mult", "2.0", "--captain-ruling", "R2"], "refused with --captain-ruling R2"),
    (["--feng-fraction-mult", "2.0", "--captain-ruling", "r20"], "refused with --captain-ruling r20"),
    (["--feng-fraction-mult", "2.0", "--captain-ruling", "R99"], "refused with --captain-ruling R99"),
    (["--feng-fraction-mult", "2.0", "--captain-ruling", "T63"], "refused with --captain-ruling T63"),
    (["--feng-fraction-mult", "1.5", "--captain-ruling", "R20"], "with --released-bet only"),
    (["--feng-fraction-mult", "1.0", "--captain-ruling", "R20"], "with --released-bet only"),
    (["--captain-ruling", "R20"], "with --released-bet only"),
    (["--g-mom-damping", "1.0", "--captain-ruling", "R20"], "with --released-bet only"),
    ([*R20, "--g-mom-damping", "1.0"], "with --g-mom-damping other than 0.5"),
    ([*R20, "--g-mom-damping", "0.75"], "with --g-mom-damping other than 0.5"),
    ([*R20, "--bet-strength", "0.5"], "with --released-bet and --bet-strength"),
    ([*R20, "--released-bet", "share_trend", "--bet-strength", "0.5"], "with --released-bet and --bet-strength"),
    ([*R20, "--seed", "1000002"], "recomputed last file only"),
    ([*R20, "--board", "T1:val"], "recomputed last file only"),
    (["--feng-fraction-mult", "2.5", "--captain-ruling", "R20"], "--feng-fraction-mult takes 1"),
    (["--feng-fraction-mult", "3", "--captain-ruling", "R20"], "--feng-fraction-mult takes 1"),
    (["--feng-fraction-mult", "0.8"], "--feng-fraction-mult takes 1"),
])
def test_every_other_use_of_the_value_or_the_row_is_refused(bt, monkeypatch, tmp_path, extra, match):
    _refused(bt, monkeypatch, tmp_path, ["t1", *RC, *extra], match)


@pytest.mark.parametrize("argv", [
    ["t1"],  # the authors' groups
    ["t1", "--g-mom", "--base", "lib3"],  # the first upload's command
    ["t1", "--g-mom", "--base", "lib3", "--feng-groups", "none"],  # the none file
    ["t1", "--g-mom", "--base", "mid", "--cm-share", "0.60"],  # T1's second upload
])
def test_the_r20_value_is_refused_on_any_other_command(bt, monkeypatch, tmp_path, argv):
    _refused(bt, monkeypatch, tmp_path, [*argv, *R20], "recomputed last file only")


@pytest.mark.parametrize("extra", [["--cm-share", "0.6"], ["--herd"], ["--no-feng-hygiene"], ["--lib3-pair"]])
def test_the_recomputed_commands_own_guards_hold_with_the_r20_value(bt, monkeypatch, tmp_path, extra):
    _refused(bt, monkeypatch, tmp_path, ["t1", *RC, *extra, *R20], "--feng-groups recomputed builds one file")


def test_the_r20_value_is_refused_when_the_row_is_not_the_captains_order(bt, monkeypatch, tmp_path):
    ledger = tmp_path / "TRIES.md"
    ledger.write_text("| R20 (awaits the captain's order; `--feng-fraction-mult 2.0` `--captain-ruling R20`) | T1 |\n")
    monkeypatch.setattr(bt, "TRIES_FILE", ledger)
    _refused(bt, monkeypatch, tmp_path, ["t1", *RC, *R20], "no single row R20")
    ledger.write_text("# no such row\n")
    _refused(bt, monkeypatch, tmp_path, ["t1", *RC, *R20], "no single row R20")


def test_the_bet_refusal_lets_only_row_r20s_command_through(bt):
    """_bet_refusals still refuses --captain-ruling without --released-bet unless the command is row R20's (which
    _t1_candidate checked first): the shortcut needs the value, the row, no bet and the order in the ledger."""
    from types import SimpleNamespace

    def ns(**k):
        base = dict(released_bet=None, bet_strength=None, captain_ruling="R20", feng_fraction_mult=2.0,
                    feng_groups="recomputed", seed=harness.BUILD_SEED_OFFSET + 1, board=None, g_mom_damping=None)
        return SimpleNamespace(**{**base, **k})

    assert bt._bet_refusals(ns()) is None
    for bad in (ns(feng_fraction_mult=1.5), ns(feng_fraction_mult=None), ns(captain_ruling="R19"),
                ns(feng_fraction_mult=1.0)):
        with pytest.raises(SystemExit, match="with --released-bet only"):
            bt._bet_refusals(bad)
