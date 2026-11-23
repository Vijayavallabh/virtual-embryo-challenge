"""`build_test.py t1 --g-mom --base lib3 --feng-groups recomputed --direction-momentum all --captain-ruling R23`:
candidate MS under row R23 of experiments/harness/TRIES.md (the captain's order of 2026-11-23, docs/PLAN.md decision
log 2026-11-23 (c)). R13's file with one change: the drawn cells whose label's group has a recomputed direction also
get B0's per-category damped mean momentum, fitted from E9.5 to E10.5's portion libraries 1 and 2, before the atlas
step and g_mom.

Held here, on synthetic stages and directions only (the build_t1 fixture of tests/test_hidden_counts.py and the
recomputed setup of tests/test_feng_relabel.py), with a temporary ledger in place of the repository's: (1) nothing that
builds today changes: the default command, candidates A, B and C, row R20's file and the three search lines' files give
the bytes and records pinned for them (the third line's pinned below, built by the builder at main ba05554 before this
option existed); (2) the MS file differs from the default file only in cells whose label's group has a direction, the
others byte for byte the default's, and it builds the same bytes twice; (3) its values are the practice arm's
construction (experiments/t1_m1/ms_read.py's `ms_population`, then recomputed_read's `rc_hook`) for the same inputs;
(4) it is taken only with --captain-ruling R23 while row R23's first cell is the captain's order naming the option, and
every other use is refused. No released file, outside-atlas file, board or score is read."""

from __future__ import annotations

import hashlib
import inspect
import re
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
from test_feng_relabel import GROUPS4, _recomputed_env  # noqa: E402
from test_hidden_counts import t1_env  # noqa: E402, F401  (the synthetic build_t1 fixture)
from test_t1_candidates_build import RC, THREE, _build, _dense, _rows_changed, _same_frame, _spy  # noqa: E402
from test_t1_damping_line_build import DAMP as LINE2  # noqa: E402
from test_t1_damping_line_build import LINE1_BEFORE  # noqa: E402
from test_t1_fraction2_build import BEFORE, R20, _record_sha  # noqa: E402
from test_t1_fraction_line_build import R20_BEFORE  # noqa: E402
from test_t1_fraction_line_build import R22 as LINE1  # noqa: E402
from test_t1_line3_build import LINE2_BEFORE, LINE3, ROW_R20, ROW_R22  # noqa: E402

from vec import harness  # noqa: E402

bt = tcb.bt  # the builder module's fixture

MS = ["--direction-momentum", "all", "--captain-ruling", "R23"]
OPTION = "`--direction-momentum all`"
# The third search line's two files on this fixture, built by the builder at main ba05554 (before this option existed),
# twice in separate processes with the same result: the file's sha256 and its record's (feng_hygiene.commit set to
# "own_repo"), the same kind of synthetic-fixture pin as those of the other lines.
LINE3_BEFORE = {2.0: ("200f3e3bf25bcb9e41a94240a3d4a922580d017e405ab4401ebb3c7690cadd77",
                      "a8b2d8072cda3b00e7a7f9a175bb67a81e15cee727d8eecb1a26d577585f7faf"),
                2.5: ("079f6bd9ca36602e2861c1421a0cdcc40ca21fa95f101c25c940e6729a866a9b",
                      "042d0783163e1d35e95669ddbf498c440ff332748f439093d33821eb57903b93")}
# Row R23's first cell, verbatim as the captain's order of 2026-11-23 put it in the ledger, inside a synthetic row
FIRST_CELL = ("R23 (the captain's order of 2026-11-23: T1's candidate MS is named, decision log 2026-11-23 (c); the "
              "captain's ruling on M1's design for this one file, after row T65 read MS AHEAD on practice; a candidate "
              "fixed before any score of it by the T1 order rule of decision log 2026-11-21 (y); ruling, no run, not a "
              "try; the T1 builder is to take `--direction-momentum all` only with `--captain-ruling R23`)")
ROW_R23 = f"| {FIRST_CELL} | T1 | Candidate MS |"
NUMBERS_ALLOWED = {"0.5", "1", "1.5", "2", "2.5", "3", "13", "20", "22", "23"}  # values, lib3, R13 and the rows' IDs
LIBS = ("1", "2", "3")


def _sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


@pytest.fixture
def ledger(bt, tmp_path, monkeypatch):
    """A temporary ledger holding the three ruling rows; the builder reads it in place of the repository's."""
    p = tmp_path / "TRIES.md"
    p.write_text(f"# ledger\n{ROW_R20}\n{ROW_R22}\n{ROW_R23}\n")
    monkeypatch.setattr(bt, "TRIES_FILE", p)
    return p


def _three_libraries(env) -> None:
    """The fixture's E10.5 cell names ending in _1, _2 and _3 in turn, as the released file's three libraries (the
    shared fixture's names end in _1 and _3 only)."""
    env.ref.obs_names = [f"c2_{i}_{LIBS[i % 3]}" for i in range(env.ref.n_obs)]


@pytest.fixture
def ms4(bt, t1_env, monkeypatch):  # noqa: F811
    _recomputed_env(bt, t1_env, monkeypatch)
    _three_libraries(t1_env)
    return _spy(bt, monkeypatch)


@pytest.fixture
def ms3(bt, t1_env, monkeypatch):  # noqa: F811
    _recomputed_env(bt, t1_env, monkeypatch, groups=THREE)
    _three_libraries(t1_env)
    return _spy(bt, monkeypatch)


@pytest.fixture
def rc4(bt, t1_env, monkeypatch):  # noqa: F811
    _recomputed_env(bt, t1_env, monkeypatch)
    return _spy(bt, monkeypatch)


def _steps(bt, monkeypatch) -> list[dict]:
    """Around the builder's MS step: per build, its arguments and what the step gets and returns."""
    seen = []
    inner = bt._t1_ms_step

    def capture(ref, library, directions, fraction, groups, held, seed):
        entry = {"ref": ref.copy(), "library": np.array(library), "directions": dict(directions),
                 "fraction": fraction, "groups": dict(groups), "held": list(held), "seed": seed}
        step, record = inner(ref, library, directions, fraction, groups, held, seed)

        def wrapped(pop):
            entry["in"] = pop.copy()
            out = step(pop)
            entry["out"] = out.copy()
            return out

        seen.append(entry)
        return wrapped, record

    monkeypatch.setattr(bt, "_t1_ms_step", capture)
    return seen


def _no_read(*a, **k):
    pytest.fail("nothing is read before the refusal")


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


# ------------------------------------------------------------------ the constants and the ruling check


def test_the_option_is_off_by_default_and_has_its_own_constants(bt):
    a = bt.build_parser().parse_args(["t1", *RC, "--out", "x.h5ad"])
    assert a.direction_momentum is None and a.captain_ruling is None
    assert bt.DIRECTION_MOMENTUM_ARMS == ("all",) and bt.T1_R23_ROW == "R23" and bt.T1_R23_OPTION == OPTION
    assert bt.T1_MS_COMMAND == ("t1 --g-mom --base lib3 --feng-groups recomputed --direction-momentum all "
                                "--captain-ruling R23")
    assert bt.T1_MS_STREAM == 1 and bt.T1_MS_DIRECTION_STREAM == 2 and bt.T1_MS_FIT_LIBRARIES == ("1", "2")
    assert bt.T1_DISSECTED_LIBRARY == "3" and bt.T1_DISSECTED_LIBRARY not in bt.T1_MS_FIT_LIBRARIES
    assert harness.B0_DAMPING == 0.5
    # the earlier options' constants are unchanged
    assert bt.G_MOM_DAMPINGS == (0.5, 1.0) and bt.FENG_FRACTION_MULTS == (1.0, 1.5)
    assert bt._t1_ms(SimpleNamespace(direction_momentum=None)) is None and bt._t1_ms(SimpleNamespace()) is None


def test_the_practice_arm_uses_the_same_constants_and_streams(bt):
    import ms_read

    assert ms_read.DAMPING == harness.B0_DAMPING
    assert ms_read.MOMENTUM_STREAM == bt.T1_MS_STREAM and ms_read.DIRECTION_STREAM == bt.T1_MS_DIRECTION_STREAM
    assert ms_read.MOM == "m1rc_gmom_mom"
    sig = inspect.signature(bt.t1_scope.momentum_shift).parameters  # the record's constants are the call's defaults
    assert sig["eps"].default == 0.5 and sig["min_cells"].default == 20 and sig["keep_depth"].default is True


def test_the_order_check_passes_on_row_r23s_first_cell_and_fails_without_the_option(bt, ledger, tmp_path):
    assert bt._captain_order_names("R23", OPTION, ledger)
    assert not bt._captain_order_names("R22", OPTION, ledger)
    p = tmp_path / "other.md"
    for text in (
            ROW_R23.replace(" `--direction-momentum all`", ""),  # the option taken out of the first cell
            ROW_R23.replace("`--direction-momentum all`", "--direction-momentum all"),  # not in backticks
            ROW_R23.replace("`--captain-ruling R23`", "`--captain-ruling R22`"),
            ROW_R23.replace("ruling, no run", "a try"),
            ROW_R23.replace("R23 (the captain's order", "R23 awaits (the captain's order"),
            "| R23 (the captain's order of X; ruling, no run) | T1 | `--direction-momentum all` "
            "`--captain-ruling R23` |",  # named outside the first cell
            ROW_R23.replace("| R23 (", "| R24 ("),  # another row
            f"{ROW_R23}\n{ROW_R23}",  # two rows with one ID
    ):
        p.write_text(f"# ledger\n{text}\n")
        assert not bt._captain_order_names("R23", OPTION, p), text


def test_the_bet_refusal_lets_only_the_ruled_command_through(bt, ledger):
    def ns(**k):
        base = dict(released_bet=None, bet_strength=None, captain_ruling="R23", feng_fraction_mult=None,
                    g_mom_damping=None, direction_momentum="all", feng_groups="recomputed",
                    seed=harness.BUILD_SEED_OFFSET + 1, board=None)
        return SimpleNamespace(**{**base, **k})

    assert bt._bet_refusals(ns()) is None
    for bad in (ns(direction_momentum=None), ns(captain_ruling="R22"), ns(g_mom_damping=1.0),
                ns(feng_fraction_mult=1.5), ns(bet_strength=0.5)):
        with pytest.raises(SystemExit, match="with --released-bet only"):
            bt._bet_refusals(bad)
    ledger.write_text(f"# ledger\n{ROW_R20}\n{ROW_R22}\n")  # no row R23: the shortcut closes
    with pytest.raises(SystemExit, match="with --released-bet only"):
        bt._bet_refusals(ns())


# ------------------------------------------------------------------ nothing existing changes


@pytest.mark.parametrize("name", list(BEFORE))
def test_the_default_and_a_b_c_give_the_same_bytes_and_records(bt, rc4, ledger, tmp_path, name):
    extra, file_sha, record_sha = BEFORE[name]
    f, rec = _build(bt, tmp_path, f"T1_rc{name}", extra)
    assert _sha(f) == file_sha and _record_sha(rec) == record_sha


def test_row_r20s_file_and_the_three_search_lines_give_the_same_bytes_and_records(bt, rc4, ledger, tmp_path):
    f, rec = _build(bt, tmp_path, "T1_s2rcB2", R20)
    assert (_sha(f), _record_sha(rec)) == R20_BEFORE
    for v in (2.5, 3.0):
        f, rec = _build(bt, tmp_path, f"T1_s2rcB{int(v * 10)}", LINE1[v])
        assert (_sha(f), _record_sha(rec)) == LINE1_BEFORE[v]
    for v in (1.5, 2.0):
        f, rec = _build(bt, tmp_path, f"T1_s2rcA{int(v * 10)}", LINE2[v])
        assert (_sha(f), _record_sha(rec)) == LINE2_BEFORE[v]
    for v in (2.0, 2.5):
        f, rec = _build(bt, tmp_path, f"T1_s2rcCm{int(v * 10)}", LINE3[v])
        assert (_sha(f), _record_sha(rec)) == LINE3_BEFORE[v]


def test_without_the_option_no_ms_step_runs(bt, ms4, ledger, monkeypatch, tmp_path):
    steps = _steps(bt, monkeypatch)
    _, rec = _build(bt, tmp_path, "T1_rc", [])
    assert steps == [] and "t1_candidate" not in rec


# ------------------------------------------------------------------ the MS file


@pytest.mark.parametrize("fixture", ["ms4", "ms3"])
def test_the_ms_file_changes_only_cells_with_a_direction_and_builds_the_same_bytes_twice(bt, ledger, tmp_path,
                                                                                         request, fixture):
    seen = request.getfixturevalue(fixture)
    with_groups = GROUPS4 if fixture == "ms4" else THREE
    f0, r0 = _build(bt, tmp_path, "T1_s2rc", [])
    fa, ra = _build(bt, tmp_path, "T1_s2rc_MSa", MS)
    fb, rb = _build(bt, tmp_path, "T1_s2rc_MSb", MS)
    assert fa.read_bytes() == fb.read_bytes() and ra == rb
    x0, xm = _same_frame(f0, fa)  # the same cells, order, obs and var as the default file
    changed = _rows_changed(x0, xm)
    pop = seen[-1]["before"]  # the drawn population MS's step gets
    labels = pop.obs["celltype"].astype(str).to_numpy()
    groups = bt.feng_dirs.label_groups(sorted(set(labels)))
    covered = np.array([groups.get(lab) in with_groups for lab in labels])
    assert covered.any() and not covered.all()
    assert changed.any() and not (changed & ~covered).any()
    assert x0[~covered].tobytes() == xm[~covered].tobytes()  # cells without a group: the default's bytes
    # the draw is the default's: the same cells in the same order with the same obs
    b0 = seen[0]["before"]
    assert list(b0.obs_names) == list(pop.obs_names) and b0.obs.equals(pop.obs)
    # MS is called without the direction step in build_submission (its own step takes the direction), at the same
    # fraction, and the default with it
    assert seen[0]["fraction"] == 0.4 and seen[0]["directions"] == sorted(with_groups)
    assert seen[-1]["fraction"] is None and seen[-1]["directions"] == []
    assert _sha(fa) not in {*(b[1] for b in BEFORE.values()), R20_BEFORE[0], *(b[0] for b in LINE1_BEFORE.values()),
                            *(b[0] for b in LINE2_BEFORE.values()), *(b[0] for b in LINE3_BEFORE.values())}


def test_the_ms_record_is_the_defaults_plus_its_candidate_field(bt, ms4, ledger, tmp_path):
    _, r0 = _build(bt, tmp_path, "T1_s2rc", [])
    _, rm = _build(bt, tmp_path, "T1_s2rc_MS", MS)
    assert {k: v for k, v in rm.items() if k != "t1_candidate"} == r0
    cand = rm["t1_candidate"]
    assert set(cand) == {"direction_momentum", "captain_ruling", "order", "momentum", "calendar_fraction",
                         "applied_fraction", "command"}
    assert cand["direction_momentum"] == "all" and cand["captain_ruling"] == "R23"
    assert "row R23" in cand["order"] and "decision log 2026-11-23 (c)" in cand["order"]
    assert cand["command"] == bt.T1_MS_COMMAND
    assert cand["calendar_fraction"] == cand["applied_fraction"] == rm["fraction"] == 0.4
    m = cand["momentum"]
    assert (m["damping"], m["eps"], m["min_cells"], m["keep_depth"]) == (0.5, 0.5, 20, True)
    assert m["stream"] == [harness.BUILD_SEED_OFFSET + 1, 1] and m["direction_stream"] == [harness.BUILD_SEED_OFFSET + 1,
                                                                                            2]
    assert m["fit"] == {"start": "E9.5_RNA", "later": "E10.5_RNA", "later_libraries": ["1", "2"],
                        "library_3_in_fit": False}
    assert "R2's union" in m["held"] and "pre-registered" in m["held"]
    assert not re.search(r"n_target|n_start|n_ref|cells\"", str(cand))  # no cell number


def test_the_fit_reads_e9_5_and_libraries_1_and_2_only(bt, ms4, ledger, monkeypatch, tmp_path, t1_env):  # noqa: F811
    steps = _steps(bt, monkeypatch)
    calls = []
    inner = bt.t1_scope.momentum_shift

    def spy(start, ref, target_like, damping, held_flat=(), **k):
        calls.append({"start": list(start.obs_names), "ref": list(ref.obs_names), "damping": damping, "kw": k})
        return inner(start, ref, target_like, damping, held_flat, **k)

    monkeypatch.setattr(bt.t1_scope, "momentum_shift", spy)
    _build(bt, tmp_path, "T1_s2rc_MS", MS)
    assert len(calls) == 1 and calls[0]["damping"] == 0.5
    assert calls[0]["start"] == list(t1_env.start.obs_names)  # E9.5, every cell
    assert calls[0]["ref"] and all(n.rsplit("_", 1)[1] in ("1", "2") for n in calls[0]["ref"])
    assert calls[0]["ref"] == [n for n in t1_env.ref.obs_names if n.rsplit("_", 1)[1] in ("1", "2")]
    assert calls[0]["kw"]["seed"] == [harness.BUILD_SEED_OFFSET + 1, 1]
    assert set(calls[0]["kw"]) == {"label_key", "seed", "start_fallback", "ref_fallback"}  # eps, min_cells, keep_depth:
    assert steps[0]["held"] == sorted(steps[0]["held"])  # the defaults; the held list is the build's own


def test_the_values_are_the_practice_arms_construction_for_the_same_inputs(bt, ms4, ledger, monkeypatch, tmp_path,
                                                                        t1_env):  # noqa: F811
    """The builder's MS step on the drawn cells equals ms_read.ms_population (the practice arm m1rc_gmom_mom's step)
    given the same population, the same E9.5 and E10.5 libraries 1 and 2 (classified as the builder classifies them),
    the same directions, fraction, groups, held genes and seed, with the shift keywords of ms_read's cut arm (the
    chamber config's stand-ins for the two file stems). And the population g_mom's hook returns equals
    recomputed_read.rc_hook (the practice g_mom) on that step's output, with the factors bets_read.gmom_fits fits on
    the same two stages. The only difference from the practice is the fit's stages: there the task's start and
    reference (E8.5 -> E9.5), here E9.5 -> E10.5's portion libraries."""
    import bets_read as br
    import ms_read
    import recomputed_read as rcr

    steps = _steps(bt, monkeypatch)
    _build(bt, tmp_path, "T1_s2rc_MS", MS)
    s = steps[0]
    spec = bt.t1_scope.load_spec()
    start = harness._shallow(t1_env.start.copy())
    bt.t1_scope.classify(start, source="E9.5_RNA", spec=spec)
    lib = np.array([n.rsplit("_", 1)[1] for n in t1_env.ref.obs_names])
    assert np.array_equal(lib, s["library"])
    later = t1_env.ref[np.flatnonzero(np.isin(lib, ["1", "2"]))].copy()
    bt.t1_scope.classify(later, source="E10.5_RNA", spec=spec)
    shift_kw = {"start_fallback": spec.fallbacks("E9.5_RNA"), "ref_fallback": spec.fallbacks("E10.5_RNA")}
    want, rep = ms_read.ms_population(s["in"].copy(), start, later, s["directions"], s["fraction"], s["groups"],
                                      s["held"], s["seed"], shift_kw)
    got = s["out"]
    assert rep is not None and any(r["status"] == "shifted" for r in rep.as_dict()["rows"])
    assert list(want.obs_names) == list(got.obs_names) and want.obs.equals(got.obs)
    assert (_dense(want.layers["counts"]) == _dense(got.layers["counts"])).all()
    assert (_dense(want.X) == _dense(got.X)).all()
    assert s["fraction"] == 0.4 and s["seed"] == harness.BUILD_SEED_OFFSET + 1
    assert sorted(s["directions"]) == sorted(GROUPS4)
    # g_mom: the practice hook on the step's output gives what the builder's hook gave
    fits = br.gmom_fits(t1_env.ref.copy(), "E10.5_RNA", t1_env.start.copy(), "E9.5_RNA")
    labels = got.obs["celltype"].astype(str).unique()
    _, covered, with_direction = rcr.direction_kit(labels, {g: None for g in GROUPS4})
    held = bt.t1_scope.held_flat_genes()
    final = rcr.rc_hook(fits["file"], s["seed"], held, with_direction, covered)(got.copy())
    after = ms4[-1]["after"]
    assert (_dense(final.layers["counts"]) == _dense(after.layers["counts"])).all()
    assert (_dense(final.X) == _dense(after.X)).all()


def test_with_every_covered_category_skipped_the_ms_file_is_the_default(bt, ms4, ledger, monkeypatch, tmp_path):
    """The step's other calls are the default's: with the momentum's factors all at 1 (min_cells above any category)
    the MS file is the default file byte for byte, so the momentum is the one change."""
    f0, _ = _build(bt, tmp_path, "T1_s2rc", [])
    inner = bt.t1_scope.momentum_shift
    monkeypatch.setattr(bt.t1_scope, "momentum_shift", lambda *a, **k: inner(*a, **{**k, "min_cells": 10 ** 9}))
    fm, _ = _build(bt, tmp_path, "T1_s2rc_MS0", MS)
    assert fm.read_bytes() == f0.read_bytes()


# ------------------------------------------------------------------ refusals


@pytest.mark.parametrize("extra, match", [
    (["--direction-momentum", "all"], "is the captain's order of row R23 alone and needs --captain-ruling R23"),
    (["--direction-momentum", "all", "--captain-ruling", "R22"], "refused with --captain-ruling R22: only row R23"),
    (["--direction-momentum", "all", "--captain-ruling", "R20"], "refused with --captain-ruling R20: only row R23"),
    (["--direction-momentum", "all", "--captain-ruling", "r23"], "refused with --captain-ruling r23: only row R23"),
    (["--direction-momentum", "all", "--captain-ruling", "R230"], "refused with --captain-ruling R230: only row R23"),
    (["--direction-momentum", "all", "--captain-ruling", "T65"], "refused with --captain-ruling T65: only row R23"),
    ([*MS, "--g-mom-damping", "1.0"], "refused with --g-mom-damping: candidate MS is one change"),
    ([*MS, "--g-mom-damping", "0.5"], "refused with --g-mom-damping: candidate MS is one change"),
    ([*MS, "--feng-fraction-mult", "1.5"], "refused with --feng-fraction-mult: candidate MS is one change"),
    ([*MS, "--feng-fraction-mult", "1.0"], "refused with --feng-fraction-mult: candidate MS is one change"),
    ([*MS, "--feng-fraction-mult", "2.0"], "refused with --feng-fraction-mult: candidate MS is one change"),
    ([*MS, "--g-mom-damping", "1.0", "--feng-fraction-mult", "1.5"],
     "refused with --g-mom-damping, --feng-fraction-mult: candidate MS"),
    ([*MS, "--released-bet", "share_trend", "--bet-strength", "0.5"],
     "refused with --released-bet, --bet-strength: candidate MS"),
    ([*MS, "--bet-strength", "0.5"], "refused with --bet-strength: candidate MS"),
    ([*MS, "--no-hold-unbracketed"], "refused with --no-hold-unbracketed: the cells without a group keep R1's hold"),
    ([*MS, "--seed", "1000002"], "changes T1's recomputed last file only"),
    ([*MS, "--board", "T1:val"], "changes T1's recomputed last file only"),  # a dry run
])
def test_every_other_use_of_the_option_is_refused_before_anything_is_read(bt, ledger, monkeypatch, tmp_path, extra,
                                                                          match):
    _refused(bt, monkeypatch, tmp_path, ["t1", *RC, *extra], match)


@pytest.mark.parametrize("argv", [
    ["t1"],  # the authors' groups
    ["t1", "--g-mom", "--base", "lib3"],  # the first upload's command
    ["t1", "--g-mom", "--base", "lib3", "--feng-groups", "none"],  # the none file
    ["t1", "--g-mom", "--base", "mid", "--cm-share", "0.60"],  # T1's second upload
])
def test_the_option_is_refused_on_any_other_command(bt, ledger, monkeypatch, tmp_path, argv):
    _refused(bt, monkeypatch, tmp_path, [*argv, *MS], "changes T1's recomputed last file only")


@pytest.mark.parametrize("extra", [["--cm-share", "0.6"], ["--herd"], ["--no-feng-hygiene"], ["--lib3-pair"],
                                   ["--strong-null"], ["--blood-share", "0.1"]])
def test_the_recomputed_commands_own_guards_hold_with_the_option(bt, ledger, monkeypatch, tmp_path, extra):
    _refused(bt, monkeypatch, tmp_path, ["t1", *RC, *extra, *MS], "--feng-groups recomputed builds one file")


def test_without_g_mom_the_recomputed_guard_refuses_the_option(bt, ledger, monkeypatch, tmp_path):
    _refused(bt, monkeypatch, tmp_path, ["t1", "--base", "lib3", "--feng-groups", "recomputed", *MS],
             "--feng-groups recomputed builds one file")


def test_the_ruling_alone_is_refused(bt, ledger, monkeypatch, tmp_path):
    _refused(bt, monkeypatch, tmp_path, ["t1", *RC, "--captain-ruling", "R23"], "with --released-bet only")


def test_another_arm_is_refused_by_the_parser(bt, ledger, monkeypatch, tmp_path):
    monkeypatch.setattr(bt, "load_stage", _no_read)
    with pytest.raises(SystemExit) as e:
        bt.main(["t1", *RC, "--direction-momentum", "held", "--captain-ruling", "R23", "--out",
                 str(tmp_path / "T1.h5ad")])
    assert e.value.code == 2


def test_the_option_is_refused_when_row_r23_does_not_name_it(bt, ledger, monkeypatch, tmp_path):
    for text in (f"{ROW_R20}\n{ROW_R22}\n",  # no row R23
                 f"{ROW_R23.replace(' `--direction-momentum all`', '')}\n",  # the row without the option
                 f"{ROW_R23}\n{ROW_R23}\n",  # two rows with one ID
                 "| R23 (awaits the captain's order; `--direction-momentum all` `--captain-ruling R23`) | T1 |\n"):
        ledger.write_text(f"# ledger\n{text}")
        _refused(bt, monkeypatch, tmp_path, ["t1", *RC, *MS], "no single row R23 .* `--direction-momentum all`")
    ledger.write_text(f"# ledger\n{ROW_R23}\n")
    a = bt.build_parser().parse_args(["t1", *RC, *MS, "--out", str(tmp_path / "T1.h5ad")])
    assert bt._t1_ms(a)["captain_ruling"] == "R23"


def test_the_option_is_refused_unless_the_portion_libraries_are_1_and_2(bt, rc4, ledger, tmp_path):
    """The shared fixture's E10.5 names end in _1 and _3 only: the fit's libraries are not there, and the build stops
    before any file is written."""
    out = tmp_path / "T1_s2rc_MS.h5ad"
    with pytest.raises(SystemExit, match=r"libraries other than 3 are \['1'\], not the fit's \['1', '2'\]"):
        bt.main(["t1", *RC, *MS, "--out", str(out)])
    assert not out.exists()
