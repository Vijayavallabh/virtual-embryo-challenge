"""`build_test.py t1 --g-mom --base lib3 --feng-groups recomputed --feng-fraction-mult 1.5 --direction-momentum all
--captain-ruling R26`: the one file of the T1 combined MS list under row R26 of experiments/harness/TRIES.md (the
captain's order of 2026-11-23, "combine MS"; docs/PLAN.md decision log 2026-11-23 (j)). Candidate MS's step exactly as
row R23 fixes it (`_t1_ms_step`, unchanged), put on candidate B's command exactly as B has it (the atlas step at 1.5 x
the calendar fraction, g_mom's within-type damping 0.5), nothing else changed.

Held here, on synthetic stages and directions only (the fixtures of tests/test_t1_ms_build.py), with a temporary ledger
in place of the repository's: (1) nothing that builds today changes: the default command, A, B, C, row R20's file, the
three search lines' files and MS give the bytes and records pinned for them (MS's pinned below, built by the builder at
main 15c938b before this option existed); (2) the combined file is "B's command with MS's step" by construction: its
cells without a group are B's bytes, its population is B's g_mom hook on MS's step at 1.5 x the calendar fraction, with
every category's momentum skipped it is B's file byte for byte, and with the multiplier forced to 1.0 it is MS's file
byte for byte; (3) the atlas step gets the momentum's output on stream [seed, 2]; (4) keep_depth holds; (5) it builds
the same bytes twice; (6) it is taken only with --captain-ruling R26 while row R26's first cell is the captain's order
naming the pair, and every other use is refused before anything is read. No released file, outside-atlas file, board
or score is read."""

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
from test_t1_candidates_build import RC, THREE, _build, _dense, _same_frame, _spy  # noqa: E402
from test_t1_damping_line_build import DAMP as LINE2  # noqa: E402
from test_t1_damping_line_build import LINE1_BEFORE  # noqa: E402
from test_t1_fraction2_build import BEFORE, R20, _record_sha  # noqa: E402
from test_t1_fraction_line_build import R20_BEFORE  # noqa: E402
from test_t1_fraction_line_build import R22 as LINE1  # noqa: E402
from test_t1_line3_build import LINE2_BEFORE, LINE3, ROW_R20, ROW_R22  # noqa: E402
from test_t1_ms_build import LINE3_BEFORE, ROW_R23  # noqa: E402
from test_t1_ms_build import MS as MS_ARGS  # noqa: E402

from vec import harness  # noqa: E402

bt = tcb.bt  # the builder module's fixture

PAIR = ["--feng-fraction-mult", "1.5", "--direction-momentum", "all"]
MSB = [*PAIR, "--captain-ruling", "R26"]
B_ARGS = ["--feng-fraction-mult", "1.5"]
OPTION = "`--feng-fraction-mult 1.5 --direction-momentum all`"
RULING = "`--captain-ruling R26`"
# MS's file on this fixture (the three-library fixture below), built by the builder at main 15c938b (before this option
# existed), twice in separate processes with the same result: the file's sha256 and its record's (feng_hygiene.commit
# set to "own_repo"), the same kind of synthetic-fixture pin as those of the other lines.
MS_BEFORE = ("88d648554592137ae02abb3d3aa377ab38bcbffadd1a2a75902aa482d3c7d5e6",
             "866f5a8bddc40f8bd3619f2b898c8f01ab335c610120f0b50403fdbb95da6516")
# A synthetic first cell of the form the list fixes for row R26 (not the ledger's own text)
FIRST_CELL = ("R26 (the captain's order of 2026-11-23, \"combine MS\": the T1 combined MS list, decision log "
              f"2026-11-23; ruling, no run, not a try; the T1 builder is to take {OPTION} only with {RULING})")
ROW_R26 = f"| {FIRST_CELL} | T1 | Combined MS |"
NUMBERS_ALLOWED = {"0.5", "1", "1.5", "2", "2.5", "3", "13", "20", "22", "23", "26"}  # values, lib3 and the rows' IDs
LIBS = ("1", "2", "3")
SEED = harness.BUILD_SEED_OFFSET + 1


def _sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


@pytest.fixture
def ledger(bt, tmp_path, monkeypatch):
    """A temporary ledger holding the four ruling rows; the builder reads it in place of the repository's."""
    p = tmp_path / "TRIES.md"
    p.write_text(f"# ledger\n{ROW_R20}\n{ROW_R22}\n{ROW_R23}\n{ROW_R26}\n")
    monkeypatch.setattr(bt, "TRIES_FILE", p)
    return p


def _three_libraries(env) -> None:
    """The fixture's E10.5 cell names ending in _1, _2 and _3 in turn, as tests/test_t1_ms_build.py names them."""
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
        entry = {"directions": dict(directions), "fraction": fraction, "groups": dict(groups), "held": list(held),
                 "seed": seed}
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


def _no_read(*a, **k):
    pytest.fail("nothing is read before the refusal")


def _refused(bt, monkeypatch, tmp_path, argv: list[str], match: str) -> str:
    monkeypatch.setattr(bt, "load_stage", _no_read)
    monkeypatch.setattr(bt.feng_dirs, "feng_directions", _no_read)
    monkeypatch.setattr(bt, "_t1_ms_step", _no_read)
    name = "T1_dry.h5ad" if "--board" in argv else "T1.h5ad"
    with pytest.raises(SystemExit, match=match) as e:
        bt.main([*argv, "--out", str(tmp_path / name)])
    msg = str(e.value)
    given = set(re.findall(r"\d+(?:\.\d+)?", " ".join(argv)))  # what the command line itself said
    assert set(re.findall(r"\d+(?:\.\d+)?", msg)) <= NUMBERS_ALLOWED | given, msg
    return msg


def _covered(bt, pop, with_groups) -> np.ndarray:
    labels = pop.obs["celltype"].astype(str).to_numpy()
    groups = bt.feng_dirs.label_groups(sorted(set(labels)))
    return np.array([groups.get(lab) in with_groups for lab in labels])


# ------------------------------------------------------------------ the constants and the ruling check


def test_the_pair_is_off_by_default_and_has_its_own_constants(bt):
    a = bt.build_parser().parse_args(["t1", *RC, "--out", "x.h5ad"])
    assert a.direction_momentum is None and a.captain_ruling is None and a.feng_fraction_mult is None
    assert bt.T1_R26_ROW == "R26" and bt.T1_R26_OPTION == OPTION and bt.T1_MSB_FILE == "T1_s2rcMSB.h5ad"
    assert bt.T1_MSB_MULT == 1.5 == bt.FENG_FRACTION_MULTS[1]  # candidate B's multiplier, no new constant
    assert bt.T1_MSB_COMMAND == ("t1 --g-mom --base lib3 --feng-groups recomputed --feng-fraction-mult 1.5 "
                                 "--direction-momentum all --captain-ruling R26")
    assert "row R26" in bt.T1_MSB_ORDER and "combine MS" in bt.T1_MSB_ORDER and "2026-11-23" in bt.T1_MSB_ORDER
    # row R23's and the earlier options' constants are unchanged
    assert bt.T1_R23_ROW == "R23" and bt.T1_R23_OPTION == "`--direction-momentum all`"
    assert bt.T1_MS_COMMAND == ("t1 --g-mom --base lib3 --feng-groups recomputed --direction-momentum all "
                                "--captain-ruling R23")
    assert bt.T1_MS_STREAM == 1 and bt.T1_MS_DIRECTION_STREAM == 2 and bt.T1_MS_FIT_LIBRARIES == ("1", "2")
    assert bt.G_MOM_DAMPINGS == (0.5, 1.0) and bt.FENG_FRACTION_MULTS == (1.0, 1.5) and harness.B0_DAMPING == 0.5
    assert bt.DIRECTION_MOMENTUM_ARMS == ("all",)


def test_the_order_check_passes_on_row_r26s_first_cell_and_fails_otherwise(bt, ledger, tmp_path):
    assert bt._captain_order_names("R26", OPTION, ledger)
    assert not bt._captain_order_names("R23", OPTION, ledger)  # row R23 does not name the pair
    p = tmp_path / "other.md"
    for text in (
            ROW_R26.replace(f" {OPTION}", ""),  # the option string taken out
            ROW_R26.replace(OPTION, "`--direction-momentum all`"),  # only MS's option
            ROW_R26.replace(OPTION, "`--feng-fraction-mult 1.5`"),  # only B's option
            ROW_R26.replace(OPTION, OPTION.strip("`")),  # not in backticks
            ROW_R26.replace(RULING, "`--captain-ruling R23`"),
            ROW_R26.replace("ruling, no run", "a try"),
            ROW_R26.replace("R26 (the captain's order", "R26 awaits (the captain's order"),
            f"| R26 (the captain's order of X; ruling, no run) | T1 | {OPTION} {RULING} |",  # outside the first cell
            ROW_R26.replace("| R26 (", "| R23 ("),  # R23 in place of R26
            ROW_R26.replace("| R26 (", "| R27 ("),  # another row
            f"{ROW_R26}\n{ROW_R26}",  # two rows with one ID
    ):
        p.write_text(f"# ledger\n{ROW_R23}\n{text}\n")
        assert not bt._captain_order_names("R26", OPTION, p), text


def test_the_bet_refusal_lets_only_the_ruled_pair_through(bt, ledger):
    def ns(**k):
        base = dict(released_bet=None, bet_strength=None, captain_ruling="R26", feng_fraction_mult=1.5,
                    g_mom_damping=None, direction_momentum="all", feng_groups="recomputed", seed=SEED, board=None)
        return SimpleNamespace(**{**base, **k})

    assert bt._bet_refusals(ns()) is None
    for bad in (ns(direction_momentum=None), ns(feng_fraction_mult=None), ns(feng_fraction_mult=1.0),
                ns(feng_fraction_mult=2.0), ns(g_mom_damping=0.5), ns(g_mom_damping=1.0), ns(bet_strength=0.5),
                ns(captain_ruling="R23", feng_fraction_mult=1.5), ns(captain_ruling="R22")):
        with pytest.raises(SystemExit, match="with --released-bet only"):
            bt._bet_refusals(bad)
    ledger.write_text(f"# ledger\n{ROW_R20}\n{ROW_R22}\n{ROW_R23}\n")  # no row R26: the shortcut closes
    with pytest.raises(SystemExit, match="with --released-bet only"):
        bt._bet_refusals(ns())
    assert bt._bet_refusals(ns(captain_ruling="R23", feng_fraction_mult=None)) is None  # MS's as before


def test_the_pair_takes_candidate_bs_values_and_names_row_r26(bt, ledger, tmp_path):
    a = bt.build_parser().parse_args(["t1", *RC, *MSB, "--out", str(tmp_path / "T1.h5ad")])
    ms = bt._t1_ms(a)
    assert ms == {"direction_momentum": "all", "captain_ruling": "R26", "g_mom_damping": 0.5,
                  "feng_fraction_mult": 1.5, "list_file": "T1_s2rcMSB.h5ad", "order": bt.T1_MSB_ORDER}
    b = bt.build_parser().parse_args(["t1", *RC, *B_ARGS, "--out", str(tmp_path / "T1.h5ad")])
    cand = bt._t1_candidate(a)
    assert cand == bt._t1_candidate(b)  # B's values, as B's build has them: damping 0.5, multiplier 1.5
    assert (cand["g_mom_damping"], cand["feng_fraction_mult"]) == (0.5, 1.5)
    assert bt._bet_refusals(a) is None


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


def test_the_ms_file_gives_the_same_bytes_and_record(bt, ms4, ledger, tmp_path):
    f, rec = _build(bt, tmp_path, "T1_s2rc_MS", MS_ARGS)
    assert (_sha(f), _record_sha(rec)) == MS_BEFORE
    assert rec["t1_candidate"]["command"] == bt.T1_MS_COMMAND and rec["t1_candidate"]["captain_ruling"] == "R23"


# ------------------------------------------------------------------ the combined file


@pytest.mark.parametrize("fixture", ["ms4", "ms3"])
def test_the_combined_file_is_bs_outside_the_groups_and_builds_the_same_bytes_twice(bt, ledger, tmp_path, request,
                                                                                   fixture):
    seen = request.getfixturevalue(fixture)
    with_groups = GROUPS4 if fixture == "ms4" else THREE
    fb, _ = _build(bt, tmp_path, "T1_s2rcB", B_ARGS)
    fa, ra = _build(bt, tmp_path, "T1_s2rcMSBa", MSB)
    f2, r2 = _build(bt, tmp_path, "T1_s2rcMSBb", MSB)
    assert fa.read_bytes() == f2.read_bytes() and ra == r2
    xb, xm = _same_frame(fb, fa)  # the same cells, order, obs and var as B's file
    pop = seen[-1]["before"]
    covered = _covered(bt, pop, with_groups)
    assert covered.any() and not covered.all()
    changed = (xb != xm).any(axis=1)
    assert changed.any() and not (changed & ~covered).any()
    assert xb[~covered].tobytes() == xm[~covered].tobytes()  # cells without a group: B's bytes
    b0 = seen[0]["before"]  # B's draw is the combined file's draw
    assert list(b0.obs_names) == list(pop.obs_names) and b0.obs.equals(pop.obs)
    # B takes the direction inside build_submission, the combined file inside MS's step, at the same fraction
    assert seen[0]["fraction"] == 1.5 * 0.4 and seen[0]["directions"] == sorted(with_groups)
    assert seen[-1]["fraction"] is None and seen[-1]["directions"] == []
    known = {*(b[1] for b in BEFORE.values()), R20_BEFORE[0], *(b[0] for b in LINE1_BEFORE.values()),
             *(b[0] for b in LINE2_BEFORE.values()), *(b[0] for b in LINE3_BEFORE.values()), MS_BEFORE[0]}
    assert _sha(fa) not in known


def test_the_population_is_bs_g_mom_hook_on_ms_step_at_bs_fraction(bt, ms4, ledger, monkeypatch, tmp_path):
    steps = _steps(bt, monkeypatch)
    hooks = _hooks(bt, monkeypatch)
    _build(bt, tmp_path, "T1_s2rcB", B_ARGS)
    _build(bt, tmp_path, "T1_s2rcMSB", MSB)
    assert len(steps) == 1 and len(hooks) == 2  # MS's step runs in the combined build only
    assert hooks[0]["kw"]["damping"] == hooks[1]["kw"]["damping"] == 0.5  # passed as B's build passes it, not None
    assert hooks[0]["kw"]["with_direction"] == hooks[1]["kw"]["with_direction"]
    s = steps[0]
    assert s["fraction"] == 1.5 * bt.feng_dirs.bracket_fraction(10.5, 12.5) == ms4[0]["fraction"]
    assert s["seed"] == SEED and sorted(s["directions"]) == sorted(GROUPS4)
    assert ms4[-1]["before"].obs.equals(s["in"].obs)  # the step gets the draw
    want = hooks[0]["hook"](s["out"].copy())  # B's own g_mom hook on MS's step's output
    got = ms4[-1]["after"]
    assert list(want.obs_names) == list(got.obs_names) and want.obs.equals(got.obs)
    assert (_dense(want.layers["counts"]) == _dense(got.layers["counts"])).all()
    assert (_dense(want.X) == _dense(got.X)).all()


def test_the_atlas_step_gets_the_momentums_output_on_stream_2(bt, ms4, ledger, monkeypatch, tmp_path):
    calls = []
    mom_inner, ext_inner = bt.t1_scope.momentum_shift, bt.t1_scope.external_shift

    def mom(*a, **k):
        out, rep = mom_inner(*a, **k)
        calls.append(("momentum", k.get("seed"), out.copy(), rep))
        return out, rep

    def ext(target_like, directions, fraction, group_of, held_flat=(), **k):
        calls.append(("atlas", k.get("seed"), target_like.copy(), fraction))
        return ext_inner(target_like, directions, fraction, group_of, held_flat, **k)

    monkeypatch.setattr(bt.t1_scope, "momentum_shift", mom)
    monkeypatch.setattr(bt.t1_scope, "external_shift", ext)
    _build(bt, tmp_path, "T1_s2rcMSB", MSB)
    assert [c[0] for c in calls] == ["momentum", "atlas"]
    assert calls[0][1] == [SEED, 1] and calls[1][1] == [SEED, 2]
    assert calls[1][3] == 1.5 * 0.4
    m, x = calls[0][2], calls[1][2]
    assert list(m.obs_names) == list(x.obs_names)
    assert (_dense(m.layers["counts"]) == _dense(x.layers["counts"])).all()  # the momentum's output, unchanged


def test_keep_depth_holds_and_each_covered_cells_expected_depth_is_kept(bt, ms4, ledger, monkeypatch, tmp_path):
    """momentum_shift and external_shift are called without keep_depth (their default True); the momentum's
    apply_factors call gets keep_depth True, and its per-cell constant keeps each cell's expected total count."""
    sig = inspect.signature
    assert sig(bt.t1_scope.momentum_shift).parameters["keep_depth"].default is True
    assert sig(bt.t1_scope.external_shift).parameters["keep_depth"].default is True
    state = {"in_momentum": False}
    seen = []
    mom_inner, ext_inner, apply_inner = (bt.t1_scope.momentum_shift, bt.t1_scope.external_shift,
                                         bt.t1_scope.apply_factors)

    def mom(*a, **k):
        assert "keep_depth" not in k
        state["in_momentum"] = True
        try:
            return mom_inner(*a, **k)
        finally:
            state["in_momentum"] = False

    def ext(*a, **k):
        assert "keep_depth" not in k
        return ext_inner(*a, **k)

    def apply(C, F, tidx, held=None, keep_depth=True, seed=0):
        out = apply_inner(C, F, tidx, held, keep_depth=keep_depth, seed=seed)
        if state["in_momentum"]:
            seen.append({"C": C, "F": np.asarray(F), "tidx": np.asarray(tidx), "held": held,
                         "keep_depth": keep_depth, "scale": out[1]})
        return out

    monkeypatch.setattr(bt.t1_scope, "momentum_shift", mom)
    monkeypatch.setattr(bt.t1_scope, "external_shift", ext)
    monkeypatch.setattr(bt.t1_scope, "apply_factors", apply)
    _, rec = _build(bt, tmp_path, "T1_s2rcMSB", MSB)
    assert rec["t1_candidate"]["momentum"]["keep_depth"] is True
    assert len(seen) == 1 and seen[0]["keep_depth"] is True
    s = seen[0]
    C = _dense(s["C"]).astype(float)
    free = np.ones(C.shape[1], dtype=bool) if s["held"] is None else ~np.asarray(s["held"], dtype=bool)
    f = s["F"][s["tidx"]] * s["scale"][:, None]
    before = (C * free).sum(axis=1)
    expected = (C * f * free).sum(axis=1)
    assert np.allclose(expected, before)  # each covered cell's expected depth is unchanged by the momentum


def test_with_every_category_skipped_the_combined_file_is_bs(bt, ms4, ledger, monkeypatch, tmp_path):
    fb, _ = _build(bt, tmp_path, "T1_s2rcB", B_ARGS)
    inner = bt.t1_scope.momentum_shift
    monkeypatch.setattr(bt.t1_scope, "momentum_shift", lambda *a, **k: inner(*a, **{**k, "min_cells": 10 ** 9}))
    fm, _ = _build(bt, tmp_path, "T1_s2rcMSB0", MSB)
    assert fm.read_bytes() == fb.read_bytes()


def test_with_the_multiplier_at_its_default_the_combined_file_is_ms(bt, ms4, ledger, monkeypatch, tmp_path):
    """B's one extra option at its default (multiplier 1.0, forced here in the candidate values; g_mom's damping 0.5 is
    B0_DAMPING, MS's own): the combined file is MS's byte for byte, so the multiplier is the only difference."""
    fms, _ = _build(bt, tmp_path, "T1_s2rc_MS", MS_ARGS)
    assert _sha(fms) == MS_BEFORE[0]
    steps = _steps(bt, monkeypatch)
    inner = bt._t1_candidate
    monkeypatch.setattr(bt, "_t1_candidate", lambda a: {**inner(a), "feng_fraction_mult": 1.0})
    f1, _ = _build(bt, tmp_path, "T1_s2rcMSB1", MSB)
    assert steps[0]["fraction"] == 0.4
    assert f1.read_bytes() == fms.read_bytes()


def test_the_combined_record_is_bs_plus_its_own_candidate_field(bt, ms4, ledger, tmp_path):
    _, rb = _build(bt, tmp_path, "T1_s2rcB", B_ARGS)
    _, rm = _build(bt, tmp_path, "T1_s2rcMSB", MSB)
    assert {k: v for k, v in rm.items() if k != "t1_candidate"} == {k: v for k, v in rb.items() if k != "t1_candidate"}
    cand = rm["t1_candidate"]
    assert set(cand) == {"direction_momentum", "captain_ruling", "g_mom_damping", "feng_fraction_mult", "list_file",
                         "order", "momentum", "calendar_fraction", "applied_fraction", "command"}
    assert (cand["direction_momentum"], cand["captain_ruling"]) == ("all", "R26")
    assert (cand["g_mom_damping"], cand["feng_fraction_mult"]) == (0.5, 1.5)
    assert cand["list_file"] == "T1_s2rcMSB.h5ad" and cand["command"] == bt.T1_MSB_COMMAND
    assert "row R26" in cand["order"] and "combine MS" in cand["order"]
    assert cand["calendar_fraction"] == 0.4 and cand["applied_fraction"] == rm["fraction"] == 1.5 * 0.4
    assert rb["t1_candidate"]["applied_fraction"] == cand["applied_fraction"]
    m = cand["momentum"]  # MS's step's record, unchanged
    assert (m["damping"], m["eps"], m["min_cells"], m["keep_depth"]) == (0.5, 0.5, 20, True)
    assert m["stream"] == [SEED, 1] and m["direction_stream"] == [SEED, 2]
    assert m["fit"] == {"start": "E9.5_RNA", "later": "E10.5_RNA", "later_libraries": ["1", "2"],
                        "library_3_in_fit": False}
    assert not re.search(r"n_target|n_start|n_ref|cells\"", str(cand))  # no cell number


# ------------------------------------------------------------------ refusals


@pytest.mark.parametrize("extra, match", [
    # the pair without R26: row R23's path refuses it, as today
    (PAIR, "is the captain's order of row R23 alone and needs --captain-ruling R23"),
    ([*PAIR, "--captain-ruling", "R23"], "refused with --feng-fraction-mult: candidate MS is one change"),
    ([*PAIR, "--captain-ruling", "R22"], "refused with --captain-ruling R22: only row R23"),
    ([*PAIR, "--captain-ruling", "R20"], "refused with --captain-ruling R20: only row R23"),
    ([*PAIR, "--captain-ruling", "r26"], "refused with --captain-ruling r26: only row R23"),
    ([*PAIR, "--captain-ruling", "R260"], "refused with --captain-ruling R260: only row R23"),
    # the pair with R26 and anything else
    ([*MSB, "--g-mom-damping", "1.0"], "refused with --g-mom-damping: row R26's file"),
    ([*MSB, "--g-mom-damping", "0.5"], "refused with --g-mom-damping: row R26's file"),
    ([*MSB, "--g-mom-damping", "1.5"], "refused with --g-mom-damping: row R26's file"),
    ([*MSB, "--released-bet", "share_trend", "--bet-strength", "0.5"],
     "refused with --released-bet, --bet-strength: row R26's file"),
    ([*MSB, "--bet-strength", "0.5"], "refused with --bet-strength: row R26's file"),
    ([*MSB, "--no-hold-unbracketed"], "refused with --no-hold-unbracketed: the cells without a group keep R1's hold"),
    ([*MSB, "--seed", "1000002"], "changes T1's recomputed last file only: .* \\(row R26\\)"),
    ([*MSB, "--board", "T1:val"], "changes T1's recomputed last file only: .* \\(row R26\\)"),  # a dry run
    # the momentum with any other multiplier stays refused, with R26 or R23
    (["--feng-fraction-mult", "1.0", "--direction-momentum", "all", "--captain-ruling", "R26"],
     "refused with --captain-ruling R26: only row R23"),
    (["--feng-fraction-mult", "2.0", "--direction-momentum", "all", "--captain-ruling", "R26"],
     "refused with --captain-ruling R26: only row R23"),
    (["--feng-fraction-mult", "2.5", "--direction-momentum", "all", "--captain-ruling", "R26"],
     "refused with --captain-ruling R26: only row R23"),
    (["--feng-fraction-mult", "3.0", "--direction-momentum", "all", "--captain-ruling", "R26"],
     "refused with --captain-ruling R26: only row R23"),
    (["--feng-fraction-mult", "2.0", "--direction-momentum", "all", "--captain-ruling", "R23"],
     "refused with --feng-fraction-mult: candidate MS is one change"),
    # R26 with any other option set
    (["--direction-momentum", "all", "--captain-ruling", "R26"], "refused with --captain-ruling R26: only row R23"),
    (["--feng-fraction-mult", "1.5", "--captain-ruling", "R26"], "with --released-bet only"),
    (["--g-mom-damping", "1.0", "--feng-fraction-mult", "1.5", "--captain-ruling", "R26"], "with --released-bet only"),
    (["--g-mom-damping", "1.0", "--captain-ruling", "R26"], "with --released-bet only"),
    (["--g-mom-damping", "1.5", "--captain-ruling", "R26"], "refused with --captain-ruling R26: only row R22"),
    (["--feng-fraction-mult", "2.0", "--captain-ruling", "R26"], "refused with --captain-ruling R26: only row R20"),
    (["--captain-ruling", "R26"], "with --released-bet only"),
])
def test_every_other_use_is_refused_before_anything_is_read(bt, ledger, monkeypatch, tmp_path, extra, match):
    _refused(bt, monkeypatch, tmp_path, ["t1", *RC, *extra], match)


@pytest.mark.parametrize("argv", [
    ["t1"],  # the authors' groups
    ["t1", "--g-mom", "--base", "lib3"],  # the first upload's command
    ["t1", "--g-mom", "--base", "lib3", "--feng-groups", "none"],  # the none file
    ["t1", "--g-mom", "--base", "mid", "--cm-share", "0.60"],  # T1's second upload
])
def test_the_pair_is_refused_on_any_other_command(bt, ledger, monkeypatch, tmp_path, argv):
    _refused(bt, monkeypatch, tmp_path, [*argv, *MSB], "changes T1's recomputed last file only")


@pytest.mark.parametrize("extra", [["--cm-share", "0.6"], ["--herd"], ["--no-feng-hygiene"], ["--lib3-pair"],
                                   ["--strong-null"], ["--blood-share", "0.1"]])
def test_the_recomputed_commands_own_guards_hold_with_the_pair(bt, ledger, monkeypatch, tmp_path, extra):
    _refused(bt, monkeypatch, tmp_path, ["t1", *RC, *extra, *MSB], "--feng-groups recomputed builds one file")


def test_without_g_mom_the_recomputed_guard_refuses_the_pair(bt, ledger, monkeypatch, tmp_path):
    _refused(bt, monkeypatch, tmp_path, ["t1", "--base", "lib3", "--feng-groups", "recomputed", *MSB],
             "--feng-groups recomputed builds one file")


def test_the_pair_is_refused_when_row_r26_does_not_name_it(bt, ledger, monkeypatch, tmp_path):
    for text in (f"{ROW_R20}\n{ROW_R22}\n{ROW_R23}\n",  # no row R26
                 f"{ROW_R23}\n{ROW_R26.replace(f' {OPTION}', '')}\n",  # the row without the option string
                 f"{ROW_R23}\n{ROW_R26.replace(OPTION, '`--direction-momentum all`')}\n",  # MS's option only
                 f"{ROW_R23}\n{ROW_R26}\n{ROW_R26}\n",  # two rows with one ID
                 f"{ROW_R23}\n{ROW_R26.replace('| R26 (', '| R23 (')}\n",  # R23 in place of R26
                 f"| R26 (awaits the captain's order; {OPTION} {RULING}) | T1 |\n"):
        ledger.write_text(f"# ledger\n{text}")
        _refused(bt, monkeypatch, tmp_path, ["t1", *RC, *MSB],
                 "no single row R26 .* `--feng-fraction-mult 1.5 --direction-momentum all`")
    ledger.write_text(f"# ledger\n{ROW_R26}\n")
    a = bt.build_parser().parse_args(["t1", *RC, *MSB, "--out", str(tmp_path / "T1.h5ad")])
    assert bt._t1_ms(a)["captain_ruling"] == "R26"


def test_bare_ms_keeps_row_r23s_path_with_or_without_row_r26(bt, ledger, tmp_path):
    a = bt.build_parser().parse_args(["t1", *RC, *MS_ARGS, "--out", str(tmp_path / "T1.h5ad")])
    want = {"direction_momentum": "all", "captain_ruling": "R23", "order": bt.T1_MS_ORDER}
    assert bt._t1_ms(a) == want and bt._t1_candidate(a) is None
    ledger.write_text(f"# ledger\n{ROW_R23}\n")
    assert bt._t1_ms(a) == want


def test_the_pair_is_refused_unless_the_portion_libraries_are_1_and_2(bt, rc4, ledger, tmp_path):
    """The shared fixture's E10.5 names end in _1 and _3 only: MS's step refuses as for MS, before any file is
    written."""
    out = tmp_path / "T1_s2rcMSB.h5ad"
    with pytest.raises(SystemExit, match=r"libraries other than 3 are \['1'\], not the fit's \['1', '2'\]"):
        bt.main(["t1", *RC, *MSB, "--out", str(out)])
    assert not out.exists()
