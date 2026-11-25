"""`experiments/builds/build_test.py heart --entry organ_scope --trend-strength 1.25 --run-profile <sha256> --z-mode raw
--captain-ruling R56`: file OXP125 (row R56 of experiments/harness/TRIES.md, the captain's order of 2026-11-24, ruling,
no run), OTS125's command with E10.5's run profile removed (row H36's test profile, by its sha256).

On the synthetic stages and synthetic test board of tests/test_heart_trend_build.py, E10.5 as
tests/test_build_heart_organ_scope.py's `organ_stage` and that file's size stand-in, with a temporary ledger (rows R28,
R29, R49 and R52 as tests/test_heart_trend_line_g.py has them, and a synthetic first cell for row R56) and a temporary
profile file: a zero profile gives OTS125's bytes; any profile keeps OTS125's rows, coordinates, cell number and size
field and, after halve.py, its half's rows, and changes only the expression; the step is _trend_strength_step's at
1.25 followed by vec.t2_build.apply_run_profile, bit for bit; the same bytes twice; the record names the file, R56, its
first cell, the command and the sha256; every refusal, before anything is read (and the profile's after the boards,
before any stage); without --run-profile every command is refused or built as before (tests/test_heart_trend_line_g.py
pins their bytes and refusals); the R56 check read by code on the committed ledger. Nothing is read from data/."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from test_build_heart_iso_trim import _build  # noqa: E402
from test_build_heart_organ_scope import _use_stage, organ_stage  # noqa: E402
from test_heart_trend_build import GENES, install, sha  # noqa: E402
from test_heart_trend_line_d import _ctx_d, _organ_size  # noqa: E402
from test_heart_trend_line_f import _X, _coords  # noqa: E402
from test_heart_trend_line_g import LEDGER as G_LEDGER  # noqa: E402
from test_heart_trend_line_g import OTSF, OTSF_BEFORE  # noqa: E402

R56_FIRST = ("R56 (the captain's order of 2026-11-24, ruling, no run: \"... and yes to XPC\", a synthetic stand-in; "
             "the "
             "heart builder is to take `--entry organ_scope` with `--trend-strength 1.25` and `--run-profile` only "
             "with `--captain-ruling R56`; not a try)")
LEDGER = G_LEDGER + f"| {R56_FIRST} | heart | file OXP125 |\n"
OTS125 = OTSF[1.25]
ORG = ["--entry", "organ_scope"]


def oxp(sha256: str, row: str = "R56") -> list[str]:
    return [*ORG, "--trend-strength", "1.25", "--run-profile", sha256, "--z-mode", "raw", "--captain-ruling", row]


@pytest.fixture(scope="module")
def bt():
    sys.path.insert(0, str(REPO / "experiments" / "builds"))
    import build_test

    return build_test


@pytest.fixture
def ledger(bt, monkeypatch, tmp_path):
    path = tmp_path / "TRIES.md"
    path.write_text(LEDGER)
    monkeypatch.setattr(bt, "TRIES_FILE", path)
    return path


@pytest.fixture
def env(bt, monkeypatch, own_repo, ledger, tmp_path):
    """The synthetic stages and board, E10.5 as `organ_stage`, the organ size stand-in and a temporary profile path;
    the list records every stage read and every board read."""
    reads: list[str] = []
    install(bt, monkeypatch.setattr, reads)
    _use_stage(bt, monkeypatch, organ_stage(), reads)
    _organ_size(bt, monkeypatch)
    boards = bt.load_boards

    def load_boards():
        reads.append("boards")
        return boards()

    monkeypatch.setattr(bt, "load_boards", load_boards)
    monkeypatch.setattr(bt, "HEART_XPC_PROFILE", tmp_path / "profile" / "xpc_profile.json")
    return reads


def profile(bt, values, **over) -> str:
    """Write a test profile at bt.HEART_XPC_PROFILE (the form xpc_table.record writes) and return its sha256."""
    rec = {"axis": "heart-xpc-profile", "kind": "test", "passed": True, "genes": list(GENES),
           "delta": [float(x) for x in values], **over}
    text = json.dumps(rec, indent=1, sort_keys=True) + "\n"
    bt.HEART_XPC_PROFILE.parent.mkdir(parents=True, exist_ok=True)
    bt.HEART_XPC_PROFILE.write_text(text)
    return hashlib.sha256(text.encode()).hexdigest()


def delta(seed: int = 0) -> np.ndarray:
    return np.clip(np.random.default_rng(seed).normal(0, 0.6, len(GENES)), -1.0, 1.0)


def _why(bt, tmp_path, *argv) -> str:
    out = tmp_path / ("h_dry.h5ad" if "--board" in argv and "T2:heart:val" in " ".join(argv) else "h.h5ad")
    with pytest.raises(SystemExit) as e:
        bt.main(["heart", *argv, "--out", str(out)])
    assert not out.exists()
    return str(e.value)


# ------------------------------------------------------------------------------------------------ off by default
def test_the_option_is_off_by_default(bt):
    a = bt.build_parser().parse_args(["heart", *OTS125, "--out", "runs/x/h.h5ad"])
    assert a.run_profile is None
    assert bt.HEART_XPC_ROW == "R56" and bt.HEART_XPC_CHOICE == "1.25" and bt.HEART_XPC_FILE == "OXP125"
    assert bt.HEART_XPC_COMMAND.format(sha256="S") == "heart " + " ".join(oxp("S"))
    assert bt.HEART_XPC_NAMES == ("`--entry organ_scope`", "`--trend-strength 1.25`", "`--run-profile`")
    assert bt.HEART_XPC_PROFILE == REPO / "experiments" / "heart_profile" / "results" / "xpc_profile.json"
    assert bt.HEART_FURTHER_ROW == "R49"  # OTS125's own row stays R49's


def test_ots125_keeps_its_bytes_and_r56_is_refused_for_it_as_before(bt, env, tmp_path):
    out, rec, _ = _build(bt, tmp_path, "ots125", *OTS125)
    assert sha(out) == OTSF_BEFORE[1.25][0] and "run_profile" not in rec and "R56" not in json.dumps(rec)
    msg = _why(bt, tmp_path, *ORG, "--trend-strength", "1.25", "--z-mode", "raw", "--captain-ruling", "R56")
    assert msg == ("--entry organ_scope --trend-strength 1.25 is refused: experiments/harness/TRIES.md must hold "
                   "exactly one row R56 whose first cell is the captain's order (\"ruling, no run\") naming "
                   "`--entry organ_scope`, `--trend-strength 0.75`, `--trend-strength 1.25` and "
                   "`--captain-ruling R56`")


# ------------------------------------------------------------------------------------------------ the file
def test_a_zero_profile_gives_ots125s_bytes(bt, env, tmp_path):
    s = profile(bt, np.zeros(len(GENES)))
    out, rec, _ = _build(bt, tmp_path, "oxp0", *oxp(s))
    assert sha(out) == OTSF_BEFORE[1.25][0]
    ts, rp = rec["trend_strength"], rec["run_profile"]
    assert ts["file"] == "OXP125" and ts["captain_ruling"] == "R56" and ts["order_row"] == R56_FIRST
    assert ts["command"] == "heart " + " ".join(oxp(s)) and ts["strength"] == 1.25
    assert rp["sha256"] == s and rp["captain_ruling"] == "R56" and rp["stream"] == 37
    assert rp["profile"] == "experiments/heart_profile/results/xpc_profile.json"


def test_rows_coordinates_size_and_half_are_ots125s_and_only_the_expression_differs(bt, env, tmp_path):
    sys.path.insert(0, str(REPO / "experiments" / "builds"))
    import halve

    s = profile(bt, delta())
    o, o_rec, _ = _build(bt, tmp_path, "ots125", *OTS125)
    x, x_rec, _ = _build(bt, tmp_path, "oxp", *oxp(s))
    A, B = ad.read_h5ad(o), ad.read_h5ad(x)
    assert A.n_obs == B.n_obs and list(A.obs_names) == list(B.obs_names) and list(A.var_names) == list(B.var_names)
    np.testing.assert_array_equal(_coords(o), _coords(x))
    assert not np.array_equal(_X(o), _X(x))
    for k in ("n_cells", "seed", "count_decision", "trim", "size", "z_mode", "board", "dry_run", "depth_step",
              "scope_option"):
        assert x_rec[k] == o_rec[k], k
    assert x_rec["trend_strength"]["trend"] == o_rec["trend_strength"]["trend"]
    ho, hx = tmp_path / "ots125_half.h5ad", tmp_path / "oxp_half.h5ad"
    halve.halve(str(o), str(ho))
    halve.halve(str(x), str(hx))
    HA, HB = ad.read_h5ad(ho), ad.read_h5ad(hx)
    assert list(HA.obs_names) == list(HB.obs_names)
    np.testing.assert_array_equal(_coords(ho), _coords(hx))


def test_the_same_bytes_twice_and_the_test_boards_own_name(bt, env, tmp_path):
    s = profile(bt, delta(1))
    a, _, ta = _build(bt, tmp_path, "a", *oxp(s))
    b, _, tb = _build(bt, tmp_path, "b", *oxp(s))
    c, _, _ = _build(bt, tmp_path, "c", *oxp(s), "--board", "T2:heart:test_extrap")
    assert a.read_bytes() == b.read_bytes() == c.read_bytes() and ta == tb


def test_the_step_is_ots125s_trend_then_the_profile(bt):
    """A zero profile is OTS125's step bit for bit; any profile is row H28's trend step at 1.25 (the counts
    h28.trend_expression returns, which the practice arms of row H36 also start from) followed by
    vec.t2_build.apply_run_profile at stream [seed, 37] and harness.finish on the same draw."""
    import importlib.util

    from vec import harness, t2_build

    ctx = _ctx_d(bt)
    ots, rec = bt._trend_strength_step(ctx, 1.25)
    zero, zrec = bt._run_profile_step(ctx, 1.25, np.zeros(len(GENES)))
    np.testing.assert_array_equal(zero.X, ots.X)
    np.testing.assert_array_equal(zero.coords, ots.coords)
    assert zrec == rec
    d = delta(2)
    st, _ = bt._run_profile_step(ctx, 1.25, d)
    np.testing.assert_array_equal(st.coords, ots.coords)
    assert not np.array_equal(st.X, ots.X)
    spec = importlib.util.spec_from_file_location("h28_for_xpc",
                                                  REPO / "experiments/harness/h28_heart_trend_strength.py")
    h28 = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = h28
    spec.loader.exec_module(h28)
    _, counts, _, _ = h28.trend_expression(ctx, 1.25)
    new, X, _ = t2_build.apply_run_profile(counts, d, [ctx.seed, 37])
    pool, _, _, coords = harness.draw_copy(ctx)
    by_hand = harness.finish(ctx, pool, X, coords)
    np.testing.assert_array_equal(st.X, by_hand.X)
    lib0, lib1 = np.asarray(counts.sum(1)).ravel(), np.asarray(new.sum(1)).ravel()
    assert abs(np.median(lib1) / np.median(lib0) - 1) < 0.05


# ------------------------------------------------------------------------------------------------ refusals
GOOD = "a" * 64


@pytest.mark.parametrize("argv,start", [
    (oxp("A" * 64), "--run-profile AAAA"), (oxp("a" * 63), "--run-profile aaa"), (oxp("g" * 64), "--run-profile ggg"),
    (oxp(GOOD + "0"), "--run-profile aaa"), (oxp(""), "--run-profile  is refused"),
    ([*ORG, "--trend-strength", "1.5", "--run-profile", GOOD, "--z-mode", "raw", "--captain-ruling", "R56"],
     "--run-profile is refused except as"),
    ([*ORG, "--trend-strength", "1.250", "--run-profile", GOOD, "--z-mode", "raw", "--captain-ruling", "R56"],
     "--run-profile is refused except as"),
    ([*ORG, "--run-profile", GOOD, "--z-mode", "raw", "--captain-ruling", "R56"], "--run-profile is refused except as"),
    (["--entry", "core_trim", "--trend-strength", "1.25", "--run-profile", GOOD, "--z-mode", "raw", "--captain-ruling",
      "R56"], "--run-profile is refused except as"),
    ([*ORG, "--trend-strength", "1.25", "--run-profile", GOOD, "--captain-ruling", "R56"],
     "--run-profile is refused except as"),
    ([*oxp(GOOD), "--trend"], "--run-profile is refused except as"),
    ([*oxp(GOOD), "--depth"], "--run-profile is refused except as"),
    ([*oxp(GOOD), "--herd"], "--run-profile is refused except as"),
    ([*oxp(GOOD), "--cfm"], "--run-profile is refused except as"),
    ([*oxp(GOOD), "--expression", "copy"], "--run-profile is refused except as"),
    ([*oxp(GOOD), "--seed", "1000002"], "--run-profile is refused except as"),
    ([*oxp(GOOD), "--board", "T2:heart:val_extrap"], "--run-profile is refused except as"),
    ([*ORG, "--trend-strength", "1.25", "--run-profile", GOOD, "--z-mode", "raw"], "--run-profile is the captain's"),
    (oxp(GOOD, "R49"), "--run-profile is refused with --captain-ruling R49"),
    (oxp(GOOD, "R52"), "--run-profile is refused with --captain-ruling R52"),
    (oxp(GOOD, "r56"), "--run-profile is refused with --captain-ruling r56"),
])
def test_every_refusal_comes_before_anything_is_read(bt, env, tmp_path, argv, start):
    msg = _why(bt, tmp_path, *argv)
    assert msg.startswith(start), msg
    assert env == []


@pytest.mark.parametrize("text", [
    G_LEDGER,  # no row R56
    G_LEDGER + f"| {R56_FIRST.replace('`--run-profile`', '--run-profile')} | heart | x |\n",
    G_LEDGER + f"| {R56_FIRST.replace('`--trend-strength 1.25`', '`--trend-strength 1.5`')} | heart | x |\n",
    G_LEDGER + f"| {R56_FIRST.replace('ruling, no run', 'ruling')} | heart | x |\n",
    G_LEDGER + f"| {R56_FIRST.replace('order', 'note')} | heart | x |\n",
    G_LEDGER + f"| {R56_FIRST} | heart | x |\n| {R56_FIRST} | heart | y |\n",
    G_LEDGER + f"| R56 (the captain's order, ruling, no run) | heart | {R56_FIRST} |\n",
])
def test_the_r56_check_refuses_a_ledger_without_the_order(bt, env, ledger, tmp_path, text):
    ledger.write_text(text)
    msg = _why(bt, tmp_path, *oxp(GOOD))
    assert msg.startswith("--run-profile is refused: experiments/harness/TRIES.md must hold exactly one row R56")
    assert env == []


def test_the_profiles_refusals_come_after_the_boards_and_before_any_stage(bt, env, tmp_path):
    msg = _why(bt, tmp_path, *oxp(GOOD))
    assert "does not exist" in msg and env == ["boards"]
    s = profile(bt, delta())
    for sha256, over, part in ((GOOD, {}, "is not the file of sha256"),
                               (None, {"passed": False}, "is not a passed test profile"),
                               (None, {"kind": "practice"}, "is not a passed test profile"),
                               (None, {"genes": list(GENES)[::-1]}, "does not name the panel's genes"),
                               (None, {"delta": [2.0] * len(GENES)}, "within log 3")):
        env.clear()
        got = profile(bt, delta(), **over) if over else s
        msg = _why(bt, tmp_path, *oxp(sha256 or got))
        assert part in msg and env == ["boards"], (over, msg)


def test_the_committed_ledger_holds_r56_in_the_form_the_builder_takes(bt):
    """experiments/harness/TRIES.md as committed, read by code only: row R56's first cell is the captain's order naming
    the three strings and `--captain-ruling R56`; and it does not order OTS125's own path (R49's stays the one)."""
    assert bt.TRIES_FILE == REPO / "experiments" / "harness" / "TRIES.md"
    assert all(bt._captain_order_names("R56", s) for s in bt.HEART_XPC_NAMES)
    assert not bt._heart_further_ordered("R56") and bt._heart_further_ordered("R49")
