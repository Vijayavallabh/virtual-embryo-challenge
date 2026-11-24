"""`experiments/builds/build_test.py heart --entry organ_scope --trend-strength 0.75 | 1.25 --z-mode raw
--captain-ruling <ROW>`: file D's two further strengths, files OTS075 and OTS125 (the design round of 2026-11-24, its
FAST proposal for the heart with 0.25 dropped by the round's challenge; its ruling row "R-F" is the row given on the
command line): the organ cut D's command with --expression copy replaced by row H28's strength-scaled trend step, on
OTS05's draw, organ-cut pool, count rule, size procedure and seed.

On the synthetic stages and synthetic test board (its limit never shown) of tests/test_heart_trend_build.py, E10.5
relabelled as tests/test_build_heart_organ_scope.py's `organ_stage` and that file's size stand-in, with a temporary
ledger in place of TRIES.md (a synthetic first cell for row R49, the ruling row; the builder pins R49 since row R52's
change, campaign/heart-ots15, and tests/test_heart_trend_line_g.py holds that pin): the log factors and the damping the
step hands momentum_shift at each strength are exactly s times those at strength 1 (as row H34's G2 holds them); the
step is row H28's own at each strength; the rows, coordinates, cell number, size field and, after halve.py at seed
1000001, the half's rows are OTS05's and only the expression differs; the same bytes twice; the record names the
option, the strength, the file, the row, the order text, the row's first cell and the command, with 0.75 and 1.25
written as given; every refusal, before anything is read; the row check on temporary ledgers. Every heart path that
existed before (rows R28's and R29's files, files C and D, the default build and the other heart commands) gives the
bytes and record it gave before, and every refusal of an old path, 0.25's included (even under a row that names it), is
word for word the one it gave before (pinned at main ed6722a, before this option; since row R52's change 1.5 is that
row's value and left this set, and the stand-in row is R49: the set is pinned again at main b42442e, before that
change). Nothing is read from data/ and nothing is scored."""

from __future__ import annotations

import functools
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from test_build_heart_iso_trim import _build, canon  # noqa: E402
from test_build_heart_organ_scope import ORGAN, _use_stage, organ_stage  # noqa: E402
from test_heart_trend_build import SECRET_MAX, install, sha  # noqa: E402
from test_heart_trend_line import ARGV, BEFORE, R28_FIRST, TS  # noqa: E402
from test_heart_trend_line_d import OTS, R29_FIRST, TS_BEFORE, _ctx_d, _organ_size  # noqa: E402

ROW = "R49"  # the ruling row given on the command line, pinned since row R52's change (its first cell synthetic here)
PAIR = (0.75, 1.25)
VALUES = ("0.75", "1.25")
ORG = ["--entry", "organ_scope"]
CORE = ["--entry", "core_trim"]
OTSF = {s: [*ORG, "--trend-strength", v, "--z-mode", "raw", "--captain-ruling", ROW]
        for s, v in zip(PAIR, VALUES, strict=True)}
FILES = {0.75: "OTS075", 1.25: "OTS125"}
F_FIRST = ("R49 (the captain's order, a synthetic stand-in for the heart line's two further strengths on file D: "
           "ruling, no run, not a try; the heart builder is to take `--entry organ_scope` with "
           "`--trend-strength 0.75` and `--trend-strength 1.25` only with `--captain-ruling R49`)")
# the same row naming 0.25 as well (the design's first form): 0.25 stays refused even so
F_FIRST_025 = F_FIRST.replace("with `--trend-strength 0.75` and",
                              "with `--trend-strength 0.25`, `--trend-strength 0.75` and")
F_STRINGS = ("`--entry organ_scope`", "`--trend-strength 0.75`", "`--trend-strength 1.25`", "`--captain-ruling R49`")
LEDGER = (f"| R27 (other) | T1 | x |\n| {R28_FIRST} | heart | the line's two files |\n"
          f"| {R29_FIRST} | heart | file D's two files |\n| {F_FIRST} | heart | file D's two further files |\n"
          "| H28 (row) | heart | x |\n")
LEDGERS = {"full": LEDGER, "no_r29": LEDGER.replace(f"| {R29_FIRST} | heart | file D's two files |\n", ""),
           "no_r28": LEDGER.replace(f"| {R28_FIRST} | heart | the line's two files |\n", ""),
           "row_names_025": LEDGER.replace(F_FIRST, F_FIRST_025)}

# sha256 of (the .h5ad, its .build.json text, the record's `canon`) that build_test.py heart wrote for row R29's two
# commands on `organ_env`'s synthetic stages and size stand-in, with LEDGER (then with row R44 in R49's place; R29's
# record holds no ledger text), computed with the builder of main at ed6722a (before this option) by the same builds as
# below. A record may differ from its text only in the last bits of floats. Rows R28's files are pinned in TS_BEFORE
# (tests/test_heart_trend_line_d.py) and every other heart command in BEFORE (tests/test_heart_trend_line.py); both are
# rechecked here.
OTS_BEFORE = {
    0.5: ("3de4e3a2d91a837ba8accf923de2d433df4d1b47284db69b75b8ad702a7df10a",
          "26ad4f98351f9c63fb76d75759ae7b1d67f6c44b20375bb521e0794af97fa1ac",
          "0ca322c2c5175e653bf8362584f8061042f435490969c89ce409a0e48e761237"),
    1.0: ("a0485f4c6da2d6798e805066a59516ee2db9bc5dcc5b5fee7a41ebb65fdd0ca5",
          "a919c214b29e07c8bc612680d7638e6aefcfa8711b8c4edb76a714fedf4ff03d",
          "a5c76b781c7a7e3ab0f154de62eb6e3a04604a490825fe0d9f107b591f0ec879"),
}
# Every refusal of an old path: case -> (argv, ledger). None uses 0.75 or 1.25, the two values this option adds, or
# 1.5, row R52's value since campaign/heart-ots15; 0.25 is an old path (refused before, refused now). The messages were
# computed with the builder of main at ed6722a (before this option); with the stand-in row R44 renamed R49 and 1.5 taken
# out by row R52's change they were computed again with the builder of main at b42442e (before that change), and
# OLD_REFUSALS_SHA256 is the sha256 of json.dumps({case: message}, sort_keys=True) over all of them. OLD_MESSAGES spells
# out one case of each kind of refusal.
OLD_REFUSALS = {
    "ruling_alone_d": ([*ORGAN, "--captain-ruling", "R29"], "full"),
    "ruling_alone_row": ([*ORG, "--z-mode", "raw", "--captain-ruling", ROW], "full"),
    "ruling_alone_s1_trend": (["--entry", "s1", "--trend", "--captain-ruling", ROW], "full"),
    "ruling_alone_c": ([*CORE, "--expression", "copy", "--z-mode", "raw", "--captain-ruling", "R28"], "full"),
    **{f"value_{v}_{e[1]}_{r}": ([*e, "--trend-strength", v, "--z-mode", "raw", "--captain-ruling", r], "full")
       for v in ("0", "1", "0.50", "-1.0", "2.0", "0.2", "0.8", "1.2", "0.250", ".25", "1.25e0", "0.125", "",
                 "0.25", "0.750", "1.250", ".75")
       for e, r in ((ORG, "R29"), (ORG, ROW), (CORE, "R28"))},
    **{f"value_0.25_row_names_025_{e[1]}": (
        [*e, "--trend-strength", "0.25", "--z-mode", "raw", "--captain-ruling", ROW], "row_names_025")
       for e in (ORG, CORE)},
    **{f"zmode_{e[1]}_{s}": ([*e, "--trend-strength", s, "--z-mode", "recorded", "--captain-ruling", r], "full")
       for s in ("0.5", "1.0") for e, r in ((ORG, "R29"), (CORE, "R28"))},
    **{f"nozmode_{e[1]}_{s}": ([*e, "--trend-strength", s, "--captain-ruling", r], "full")
       for s in ("0.5", "1.0") for e, r in ((ORG, "R29"), (CORE, "R28"))},
    **{f"flag_{f[2:]}_{s}": ([*OTS[float(s)], *f.split()], "full")
       for f in ("--expression copy", "--trend", "--depth", "--herd", "--cfm") for s in ("0.5", "1.0")},
    **{f"entry_{e}_{s}": (["--entry", e, "--trend-strength", s, "--z-mode", "raw", "--captain-ruling", "R29"], "full")
       for e in ("s1", "s2", "s2_trim", "core_trim") for s in ("0.5", "1.0")},
    **{f"entry_{e}_{s}_r28": (["--entry", e, "--trend-strength", s, "--z-mode", "raw", "--captain-ruling", "R28"],
                              "full")
       for e in ("s1", "s2", "s2_trim", "organ_scope") for s in ("0.5", "1.0")},
    **{f"noruling_{e[1]}_{s}": ([*e, "--trend-strength", s, "--z-mode", "raw"], "full")
       for s in ("0.5", "1.0", "0.25") for e in (ORG, CORE)},
    **{f"ruling_{r}_{e[1]}_{s}": ([*e, "--trend-strength", s, "--z-mode", "raw", "--captain-ruling", r], "full")
       for s in ("0.5", "1.0") for e in (ORG, CORE) for r in ("R21", "R30", "r29", ROW, "R45")},
    **{f"board_{b}_{s}": ([*OTS[float(s)], "--board", b], "full")
       for b in ("T2:heart:val_extrap", "T2:heart:val_interp") for s in ("0.5", "1.0")},
    **{f"board_{b}_{s}_c": ([*TS[float(s)], "--board", b], "full")
       for b in ("T2:heart:val_extrap", "T2:heart:val_interp") for s in ("0.5", "1.0")},
    **{f"ledger_no_r29_{s}": (OTS[float(s)], "no_r29") for s in ("0.5", "1.0")},
    **{f"ledger_no_r28_{s}": (TS[float(s)], "no_r28") for s in ("0.5", "1.0")},
}
OLD_REFUSALS_SHA256 = "7d05fb1e624a04794c52e2f64d22b02703c9f5e406ec5ba3ea92cb1da31ebac2"
_COMMANDS = ("heart --entry core_trim --trend-strength 0.5 --z-mode raw --captain-ruling R28, heart --entry core_trim "
             "--trend-strength 1.0 --z-mode raw --captain-ruling R28, heart --entry organ_scope --trend-strength 0.5 "
             "--z-mode raw --captain-ruling R29 or heart --entry organ_scope --trend-strength 1.0 --z-mode raw "
             "--captain-ruling R29")
OLD_MESSAGES = {
    "ruling_alone_d": f"heart --captain-ruling goes with --trend-strength only ({_COMMANDS})",
    "ruling_alone_row": f"heart --captain-ruling goes with --trend-strength only ({_COMMANDS})",
    "value_2.0_organ_scope_R29": "--trend-strength 2.0 is refused: it takes 0.5 or 1.0 only",
    "value_0.25_organ_scope_R49": "--trend-strength 0.25 is refused: it takes 0.5 or 1.0 only",
    "value_0.25_row_names_025_organ_scope": "--trend-strength 0.25 is refused: it takes 0.5 or 1.0 only",
    "value_0.8_organ_scope_R49": "--trend-strength 0.8 is refused: it takes 0.5 or 1.0 only",
    "value_1.250_organ_scope_R49": "--trend-strength 1.250 is refused: it takes 0.5 or 1.0 only",
    "noruling_organ_scope_0.25": "--trend-strength 0.25 is refused: it takes 0.5 or 1.0 only",
    "zmode_organ_scope_0.5": (f"--trend-strength is refused except as {_COMMANDS}: --entry core_trim or organ_scope "
                              "and --z-mode raw only, and with --expression copy, --trend, --depth, --herd and --cfm"),
    "board_T2:heart:val_extrap_0.5": ("--trend-strength is refused on T2:heart:val_extrap: the heart search line "
                                      "builds on the test board only (no dry run)"),
    "noruling_organ_scope_0.5": ("--entry organ_scope --trend-strength 0.5 is the captain's order of row R29 alone and "
                                 "needs --captain-ruling R29"),
    "ruling_R49_organ_scope_0.5": ("--entry organ_scope --trend-strength 0.5 is refused with --captain-ruling R49: "
                                   "only row R29 orders it"),
    "ledger_no_r29_0.5": ("--entry organ_scope --trend-strength 0.5 is refused: experiments/harness/TRIES.md must hold "
                          "exactly one row R29 whose first cell is the captain's order (\"ruling, no run\") naming "
                          "`--entry organ_scope`, `--trend-strength 0.5`, `--trend-strength 1.0` and "
                          "`--captain-ruling R29`"),
    "noruling_core_trim_1.0": ("--trend-strength 1.0 is the captain's order of row R28 alone and needs "
                               "--captain-ruling R28"),
    "ruling_R49_core_trim_1.0": "--trend-strength 1.0 is refused with --captain-ruling R49: only row R28 orders it",
    "ledger_no_r28_1.0": ("--trend-strength 1.0 is refused: experiments/harness/TRIES.md must hold exactly one row R28 "
                          "whose first cell is the captain's order (\"ruling, no run\") naming `--trend-strength 0.5`, "
                          "`--trend-strength 1.0` and `--captain-ruling R28`"),
}


@pytest.fixture(scope="module")
def bt():
    sys.path.insert(0, str(REPO / "experiments" / "builds"))
    import build_test

    return build_test


@pytest.fixture(scope="module")
def h28():
    spec = importlib.util.spec_from_file_location("h28_for_heart_trend_line_f",
                                                  REPO / "experiments" / "harness" / "h28_heart_trend_strength.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture
def ledger(bt, monkeypatch, tmp_path):
    path = tmp_path / "TRIES.md"
    path.write_text(LEDGER)
    monkeypatch.setattr(bt, "TRIES_FILE", path)
    return path


@pytest.fixture
def organ_env(bt, monkeypatch, own_repo, ledger):
    """install()'s stages and board, E10.5 as `organ_stage` (file D's synthetic input) and the organ size stand-in;
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
    return reads


def _X(path):
    from vec.io import dense

    return dense(ad.read_h5ad(path).X).astype(np.float32)


def _coords(path):
    from vec.io import COORD_KEY

    return np.asarray(ad.read_h5ad(path).obsm[COORD_KEY])


def _recorder(bt, monkeypatch) -> list:
    """momentum_shift as it is (its signature too: vec.harness reads its keyword defaults), with each call's damping and
    log-factor mapping copied into the returned list."""
    calls: list = []
    own = bt.t1_scope.momentum_shift

    @functools.wraps(own)
    def recording(start, end, like, damping, held, **kw):
        calls.append((damping, {k: np.array(v, copy=True) for k, v in kw["logf_by_type"].items()}))
        return own(start, end, like, damping, held, **kw)

    monkeypatch.setattr(bt.t1_scope, "momentum_shift", recording)
    return calls


# ------------------------------------------------------------------------------------------------ off by default
def test_the_option_is_off_by_default_and_the_old_tables_stand(bt):
    a = bt.build_parser().parse_args(["heart", *ORGAN, "--out", "runs/x/h.h5ad"])
    assert a.trend_strength is None and a.captain_ruling is None
    assert bt.HEART_FURTHER_CHOICES == VALUES and bt.HEART_FURTHER_FILES == FILES
    assert {s: c.format(row=ROW) for s, c in bt.HEART_FURTHER_COMMANDS.items()} == {
        s: "heart " + " ".join(v) for s, v in OTSF.items()}
    # the strength is written as given: the old tables' one decimal would write 0.75 and 1.25 as 0.8 and 1.2
    assert [f"{s:.1f}" for s in PAIR] == ["0.8", "1.2"]
    assert all(f"--trend-strength {v} " in bt.HEART_FURTHER_COMMANDS[float(v)] for v in VALUES)
    assert bt.HEART_FURTHER_OPTIONS == {float(v): f"`--trend-strength {v}`" for v in VALUES}
    assert "0.25" not in bt.HEART_FURTHER_CHOICES and 0.25 not in bt.HEART_FURTHER_FILES
    # rows R28's and R29's tables and file D's command are as before
    assert bt.HEART_TS_CHOICES == ("0.5", "1.0") and bt.HEART_TS_ROW == "R28" and bt.HEART_OTS_ROW == "R29"
    assert bt.HEART_TS_COMMANDS == {s: "heart " + " ".join(v) for s, v in TS.items()}
    assert bt.HEART_OTS_COMMANDS == {s: "heart " + " ".join(v) for s, v in OTS.items()}
    assert bt.HEART_TS_FILES == {0.5: "TS05", 1.0: "TS10"} and bt.HEART_OTS_FILES == {0.5: "OTS05", 1.0: "OTS10"}
    assert bt.HEART_ORGAN_SCOPE == "heart " + " ".join(ORGAN)


# ------------------------------------------------------------------------------------------------ old paths unchanged
@pytest.mark.parametrize("name", sorted(BEFORE))
def test_every_heart_command_before_the_option_keeps_its_bytes_and_record(bt, own_repo, ledger, monkeypatch,
                                                                         tmp_path, name):
    install(bt, monkeypatch.setattr, [])
    if name.startswith("organ"):
        _use_stage(bt, monkeypatch, organ_stage(), [])
        _organ_size(bt, monkeypatch)
    out, rec, text = _build(bt, tmp_path, name, *ARGV[name])
    assert sha(out) == BEFORE[name][0]
    assert hashlib.sha256(text.encode()).hexdigest() == BEFORE[name][1] or canon(text) == BEFORE[name][2]
    assert "trend_strength" not in rec and ROW not in text and "OTS0" not in text and "OTS1" not in text


@pytest.mark.parametrize("s", [0.5, 1.0])
def test_row_r28s_two_files_keep_their_bytes_and_record(bt, own_repo, ledger, monkeypatch, tmp_path, s):
    install(bt, monkeypatch.setattr, [])
    out, rec, text = _build(bt, tmp_path, f"ts{s:g}", *TS[s])
    assert sha(out) == TS_BEFORE[s][0]
    assert hashlib.sha256(text.encode()).hexdigest() == TS_BEFORE[s][1] or canon(text) == TS_BEFORE[s][2]
    assert rec["trend_strength"]["captain_ruling"] == "R28" and "order_row" not in rec["trend_strength"]
    assert ROW not in text


@pytest.mark.parametrize("s", [0.5, 1.0])
def test_row_r29s_two_files_keep_their_bytes_and_record(bt, organ_env, tmp_path, s):
    out, rec, text = _build(bt, tmp_path, f"ots{s:g}", *OTS[s])
    assert sha(out) == OTS_BEFORE[s][0]
    assert hashlib.sha256(text.encode()).hexdigest() == OTS_BEFORE[s][1] or canon(text) == OTS_BEFORE[s][2]
    assert rec["trend_strength"]["captain_ruling"] == "R29" and "order_row" not in rec["trend_strength"]
    assert rec["trend_strength"]["file"] == {0.5: "OTS05", 1.0: "OTS10"}[s] and ROW not in text


def test_every_old_refusal_is_word_for_word_as_before(bt, organ_env, ledger, tmp_path):
    msgs = {}
    for name, (argv, led) in OLD_REFUSALS.items():
        ledger.write_text(LEDGERS[led])
        out = tmp_path / ("h_dry.h5ad" if "--board" in argv else "h.h5ad")
        with pytest.raises(SystemExit) as e:
            bt.main(["heart", *argv, "--out", str(out)])
        assert organ_env == [] and not out.exists(), name
        msgs[name] = str(e.value)
    assert len(msgs) == len(OLD_REFUSALS) == 129
    for name, msg in OLD_MESSAGES.items():
        assert msgs[name] == msg, name
    assert hashlib.sha256(json.dumps(msgs, sort_keys=True).encode()).hexdigest() == OLD_REFUSALS_SHA256


# ------------------------------------------------------------------------------------------------ the step
def test_the_factor_at_each_strength_is_s_times_the_strength_one_factor(bt, organ_env, monkeypatch, tmp_path):
    """Through the builder's own commands: the log-factor vectors and the damping the step hands momentum_shift at 0.75
    and 1.25 are exactly s times those at strength 1 (row R29's OTS10; as row H34's G2 holds them, no tolerance), over
    the same lineages, and the strength-1 vectors are vec.harness.b0_trend_factors' own on D's context."""
    from vec import harness

    calls = _recorder(bt, monkeypatch)
    _build(bt, tmp_path, "ots10", *OTS[1.0])
    for s in PAIR:
        _build(bt, tmp_path, f"f{s:g}", *OTSF[s])
    assert len(calls) == 3  # one step per build
    (d1, m1), *rest = calls
    assert d1 == harness.B0_DAMPING and m1 and max(float(np.abs(v).max()) for v in m1.values()) > 0.01
    for s, (d, m) in zip(PAIR, rest, strict=True):
        assert d == s * d1, s
        assert sorted(m) == sorted(m1), s
        for k in m1:
            assert m[k].dtype == m1[k].dtype == np.float64
            np.testing.assert_array_equal(m[k], s * m1[k])
            assert np.array_equal(m1[k] != 0, m[k] != 0) and (not np.any(m1[k]) or not np.array_equal(m[k], m1[k]))
    ctx = _ctx_d(bt)
    pool, rows, _, _ = harness.draw_copy(ctx)
    merged = np.asarray(pool.obs["merged"]).astype(str)[rows]
    types = sorted(set(np.where(merged == "Untracked", "Untracked@end", merged)))
    logf, _ = harness.b0_trend_factors(ctx, list(pool.var_names), types)
    assert sorted(logf) == sorted(m1)
    for k in m1:
        np.testing.assert_array_equal(m1[k], logf[k])


@pytest.mark.parametrize("s", PAIR)
def test_the_step_is_row_h28s_at_each_strength_on_ds_context(bt, h28, s):
    ctx = _ctx_d(bt)
    st, rec = bt._trend_strength_step(ctx, s)
    ref, _, _, ref_rec = h28.trend_expression(ctx, s)
    np.testing.assert_array_equal(st.X, ref.X)
    np.testing.assert_array_equal(st.coords, ref.coords)
    assert rec["budget_days"] == pytest.approx(0.5) and ref_rec["effective_budget_days"] == pytest.approx(0.5 * s)


def test_rows_and_coordinates_are_ots05s_and_only_the_expression_differs(bt, organ_env, tmp_path):
    """Rows, coordinates, cell number, size field, D's own record fields and, after halve.py at seed 1000001, the
    half's rows and coordinates are OTS05's (and so D's); the expression differs from OTS05's, OTS10's and the other
    strength's."""
    sys.path.insert(0, str(REPO / "experiments" / "builds"))
    import halve

    o, o_rec, _ = _build(bt, tmp_path, "ots05", *OTS[0.5])
    ten, _, _ = _build(bt, tmp_path, "ots10", *OTS[1.0])
    built = {s: _build(bt, tmp_path, f"f{s:g}", *OTSF[s]) for s in PAIR}
    A = ad.read_h5ad(o)
    assert halve.DEFAULT_SEED == 1000001
    ho = tmp_path / "ots05_half.h5ad"
    halve.halve(str(o), str(ho))
    HA = ad.read_h5ad(ho)
    for s, (t, t_rec, _) in built.items():
        B = ad.read_h5ad(t)
        assert A.n_obs == B.n_obs and list(A.obs_names) == list(B.obs_names), s
        assert list(A.var_names) == list(B.var_names), s
        np.testing.assert_array_equal(_coords(o), _coords(t))
        for k in ("n_cells", "seed", "count_decision", "trim", "size", "z_mode", "board", "dry_run", "depth_step",
                  "scope_option"):
            assert t_rec[k] == o_rec[k], (s, k)
        ht = tmp_path / f"f{s:g}_half.h5ad"
        halve.halve(str(t), str(ht))
        HB = ad.read_h5ad(ht)
        assert list(HA.obs_names) == list(HB.obs_names) and HA.n_obs == HB.n_obs, s
        np.testing.assert_array_equal(_coords(ho), _coords(ht))
        assert not np.array_equal(_X(ho), _X(ht)), s
    Xs = {"0.5": _X(o), "1.0": _X(ten), **{f"{s:g}": _X(t) for s, (t, _, _) in built.items()}}
    names = sorted(Xs)
    for i, p in enumerate(names):
        for q in names[i + 1:]:
            assert not np.array_equal(Xs[p], Xs[q]), (p, q)


@pytest.mark.parametrize("s", PAIR)
def test_the_same_bytes_twice(bt, organ_env, tmp_path, s):
    a, _, ta = _build(bt, tmp_path, "a", *OTSF[s])
    b, _, tb = _build(bt, tmp_path, "b", *OTSF[s])
    assert a.read_bytes() == b.read_bytes() and ta == tb


@pytest.mark.parametrize("s", PAIR)
def test_the_record_names_option_strength_file_row_order_and_command(bt, organ_env, tmp_path, capsys, s):
    _, rec, text = _build(bt, tmp_path, "t", *OTSF[s])
    ts = rec["trend_strength"]
    v = VALUES[PAIR.index(s)]
    assert ts["option"] == "--trend-strength" and ts["strength"] == s and ts["captain_ruling"] == ROW
    assert ts["file"] == FILES[s] and json.loads(text)["trend_strength"]["strength"] == float(v)
    assert f'"strength": {v},' in text and f'"strength": {s:.1f},' not in text  # 0.75 and 1.25, never 0.8 or 1.2
    assert ts["order"].startswith(f"experiments/harness/TRIES.md row {ROW}, the captain's order")
    assert all(f"`--trend-strength {u}`" in ts["order"] for u in VALUES) and f"`--captain-ruling {ROW}`" in ts["order"]
    assert "OTS075 and OTS125" in ts["order"] and "nothing else changed" in ts["order"]
    assert "0.25 dropped" in ts["order"] and "`--trend-strength 0.25`" not in ts["order"]
    assert ts["order_row"] == F_FIRST  # the row's first cell as the ledger holds it
    assert ts["command"] == "heart " + " ".join(OTSF[s]) and f"--trend-strength {v} " in ts["command"]
    assert f"--trend-strength {s:.1f} " not in text
    assert "D's draw" in ts["step"]
    assert ts["budget_days"] == pytest.approx(0.5) and ts["test_budget_days"] == pytest.approx(0.5 * s)
    assert ts["trend"]["rule"] == "vec.harness.B0_TREND_RULE"
    assert rec["scope_option"]["command"] == bt.HEART_ORGAN_SCOPE and "R21" in rec["scope_option"]["ruling"]
    assert rec["entry"] == f"trend strength {v} organ_scope (organ)" and "expression" not in rec
    assert rec["z_mode"] == "raw" and rec["n_cells"] == "max_cells" and "trend" not in rec
    assert "R28" not in text and "R29" not in text
    assert str(SECRET_MAX) not in text + capsys.readouterr().out


def test_the_reads_are_ots05s(bt, organ_env, tmp_path):
    """The two read what OTS05 reads, in the same order: the boards, E10.5, E9.5 and the three earlier stages."""
    _build(bt, tmp_path, "ots05", *OTS[0.5])
    want = list(organ_env)
    assert want == ["boards", "E10.5.h5ad", "E9.5.h5ad", "E8.75.h5ad", "E8.5.h5ad", "E8.25_late.h5ad", "boards"]
    for s in PAIR:
        organ_env.clear()
        _build(bt, tmp_path, f"f{s:g}", *OTSF[s])
        assert organ_env == want, s


# ------------------------------------------------------------------------------------------------ refusals
NEW_REFUSED = [
    *([*ORG, "--trend-strength", v, "--z-mode", "raw"] for v in VALUES),  # no ROW
    *([*ORG, "--trend-strength", "0.75", "--z-mode", "raw", "--captain-ruling", r]  # not a ruling row's id
      for r in ("r49", "R-F", "49", "R49 ", " R49", "R49;", "R", "RR49", "R49ab", "R4.9", "", "R49`")),
    *([*ORG, "--trend-strength", v, "--z-mode", "raw", "--captain-ruling", r]  # a row that does not order the line
      for v in VALUES for r in ("R45", "R29", "R28", "R27", "R4", "R490")),
    *(["--entry", e, "--trend-strength", v, "--z-mode", "raw", "--captain-ruling", ROW]
      for e in ("s1", "s2", "s2_trim", "core_trim") for v in VALUES),
    *([*ORG, "--trend-strength", v, "--captain-ruling", ROW] for v in VALUES),
    *([*ORG, "--trend-strength", v, "--z-mode", "recorded", "--captain-ruling", ROW] for v in VALUES),
    *([*OTSF[s], *f.split()] for s in PAIR for f in ("--expression copy", "--trend", "--depth", "--herd", "--cfm")),
    *([*ORG, "--trend-strength", v, "--z-mode", "raw", "--captain-ruling", ROW]  # values outside the two
      for v in ("0.25", "0.2", "0.8", "1.2", "0.750", "1.250", ".75", "1.25e0", "0.125", "1.5", "0.5 ", "2.5", "0.5",
                "1.0")),
    [*ORG, "--z-mode", "raw", "--captain-ruling", ROW], [*ORGAN, "--captain-ruling", ROW],
]


@pytest.mark.parametrize("argv", NEW_REFUSED)
def test_every_other_use_is_refused_before_anything_is_read(bt, organ_env, tmp_path, argv):
    with pytest.raises(SystemExit) as e:
        bt.main(["heart", *argv, "--out", str(tmp_path / "h.h5ad")])
    assert "refused" in str(e.value) or "needs --captain-ruling" in str(e.value) or "goes with" in str(e.value)
    assert organ_env == [] and not (tmp_path / "h.h5ad").exists()


def test_the_two_values_refusals_say_why(bt, organ_env, ledger, tmp_path):
    def why(*argv):
        with pytest.raises(SystemExit) as e:
            bt.main(["heart", *argv, "--out", str(tmp_path / "h.h5ad")])
        return str(e.value)

    for v in VALUES:
        assert why(*ORG, "--trend-strength", v, "--z-mode", "raw") == (
            f"--entry organ_scope --trend-strength {v} needs --captain-ruling <ROW>, the row of "
            "experiments/harness/TRIES.md that is the captain's order for file D's two further strengths")
        assert why(*ORG, "--trend-strength", v, "--z-mode", "raw", "--captain-ruling", "R45") == (
            f"--entry organ_scope --trend-strength {v} is refused: experiments/harness/TRIES.md must hold exactly one "
            "row R45 whose first cell is the captain's order (\"ruling, no run\") naming `--entry organ_scope`, "
            "`--trend-strength 0.75`, `--trend-strength 1.25` and `--captain-ruling R45`")
        assert why(*CORE, "--trend-strength", v, "--z-mode", "raw", "--captain-ruling", ROW) == (
            f"--trend-strength {v} is refused except as heart --entry organ_scope --trend-strength {v} --z-mode raw "
            "--captain-ruling <ROW>: --entry organ_scope and --z-mode raw only, and with --expression copy, --trend, "
            "--depth, --herd and --cfm")
    assert why(*ORG, "--trend-strength", "0.75", "--z-mode", "raw", "--captain-ruling", "r49") == (
        "--entry organ_scope --trend-strength 0.75 is refused with --captain-ruling r49: <ROW> is a ruling row's id, R "
        "and its number (e.g. R44)")
    assert why(*ORG, "--trend-strength", "1.25", "--z-mode", "raw", "--captain-ruling", "R29").startswith(
        "--entry organ_scope --trend-strength 1.25 is refused: experiments/harness/TRIES.md must hold exactly one row "
        "R29 ")
    # 0.25 is refused as before, word for word, even under a row that names it
    for led in ("full", "row_names_025"):
        ledger.write_text(LEDGERS[led])
        assert why(*ORG, "--trend-strength", "0.25", "--z-mode", "raw", "--captain-ruling", ROW) == (
            "--trend-strength 0.25 is refused: it takes 0.5 or 1.0 only")
    assert why(*ORG, "--trend-strength", "0.8", "--z-mode", "raw", "--captain-ruling", ROW) == (
        "--trend-strength 0.8 is refused: it takes 0.5 or 1.0 only")
    assert organ_env == []


@pytest.mark.parametrize("board", ["T2:heart:val_extrap", "T2:heart:val_interp", "T2:embryo:test_interp"])
@pytest.mark.parametrize("s", PAIR)
def test_any_other_board_is_refused_before_anything_is_read(bt, organ_env, tmp_path, board, s):
    with pytest.raises(SystemExit) as e:
        bt.main(["heart", *OTSF[s], "--board", board, "--out", str(tmp_path / "h_dry.h5ad")])
    # a heart validation board meets the line's own refusal; another task's board the builder's earlier one
    assert "test board only" in str(e.value) or (board.startswith("T2:embryo") and "takes only" in str(e.value))
    assert organ_env == [] and not (tmp_path / "h_dry.h5ad").exists()


def test_the_test_boards_own_name_is_its_default(bt, organ_env, tmp_path):
    out, rec, _ = _build(bt, tmp_path, "named", *OTSF[0.75], "--board", "T2:heart:test_extrap")
    default, _, _ = _build(bt, tmp_path, "default", *OTSF[0.75])
    assert out.read_bytes() == default.read_bytes() and rec["dry_run"] is False


# ------------------------------------------------------------------------------------------------ the row check
def _use(bt, monkeypatch, path: Path, text: str) -> None:
    path.write_text(text)
    monkeypatch.setattr(bt, "TRIES_FILE", path)


def test_the_row_check_passes_on_a_ledger_of_the_required_form(bt, monkeypatch, tmp_path):
    _use(bt, monkeypatch, tmp_path / "ok.md", LEDGER)
    assert bt._heart_further_ordered(ROW) and bt._heart_d_line_ordered() and bt._heart_line_ordered()
    _use(bt, monkeypatch, tmp_path / "only.md", f"| {F_FIRST} | heart | x |\n")
    assert bt._heart_further_ordered(ROW) and not bt._heart_d_line_ordered() and not bt._heart_line_ordered()
    # a row that also names 0.25 still orders the two (0.25 itself stays refused: test_the_two_values_refusals_say_why)
    _use(bt, monkeypatch, tmp_path / "with025.md", f"| {F_FIRST_025} | heart | x |\n")
    assert bt._heart_further_ordered(ROW)
    # the check reads the row given: another id of the form, in a row of its own (the builder then refuses any row but
    # R49 after this check, tests/test_heart_trend_line_g.py)
    _use(bt, monkeypatch, tmp_path / "r51.md", f"| {F_FIRST.replace('R49', 'R51')} | heart | x |\n")
    assert bt._heart_further_ordered("R51") and not bt._heart_further_ordered(ROW)
    _use(bt, monkeypatch, tmp_path / "r49a.md", f"| {F_FIRST.replace('R49', 'R49a')} | heart | x |\n")
    assert bt._heart_further_ordered("R49a") and not bt._heart_further_ordered(ROW)


@pytest.mark.parametrize("drop", [*F_STRINGS, "ruling, no run", "the captain's order"])
def test_the_row_check_fails_without_each_required_string(bt, monkeypatch, tmp_path, drop):
    first = F_FIRST.replace(drop, "something else")
    assert first != F_FIRST
    _use(bt, monkeypatch, tmp_path / "t.md", LEDGER.replace(F_FIRST, first))
    assert not bt._heart_further_ordered(ROW) and bt._heart_d_line_ordered() and bt._heart_line_ordered()


@pytest.mark.parametrize("moved", F_STRINGS)
def test_the_row_check_needs_the_strings_in_the_first_cell(bt, monkeypatch, tmp_path, moved):
    first = F_FIRST.replace(moved, "something else")
    _use(bt, monkeypatch, tmp_path / "t.md", f"| {first} | heart | {moved} | {' '.join(F_STRINGS)} |\n")
    assert not bt._heart_further_ordered(ROW)


def test_the_row_check_fails_with_two_rows_another_rows_number_or_another_form(bt, monkeypatch, tmp_path):
    _use(bt, monkeypatch, tmp_path / "two.md", LEDGER + f"| {F_FIRST} | heart | y |\n")
    assert not bt._heart_further_ordered(ROW)
    _use(bt, monkeypatch, tmp_path / "other.md", LEDGER.replace("| R49 (", "| R45 ("))
    assert not bt._heart_further_ordered(ROW) and not bt._heart_further_ordered("R45")  # R45's cell names R49's ruling
    _use(bt, monkeypatch, tmp_path / "longer.md", LEDGER.replace(F_FIRST, F_FIRST.replace("R49", "R490")))
    assert not bt._heart_further_ordered(ROW) and bt._heart_further_ordered("R490")
    _use(bt, monkeypatch, tmp_path / "ruled.md", LEDGER.replace("`--captain-ruling R49`", "`--captain-ruling R45`"))
    assert not bt._heart_further_ordered(ROW)
    _use(bt, monkeypatch, tmp_path / "lower.md", f"| {F_FIRST.replace('R49', 'r49')} | heart | x |\n")
    assert not bt._heart_further_ordered("r49")  # not a ruling row's id, whatever the ledger holds
    _use(bt, monkeypatch, tmp_path / "none.md", "| R27 (x) | T1 | x |\n")
    assert not bt._heart_further_ordered(ROW)
    for strength in ("`--trend-strength 0.75`", "`--trend-strength 1.25`"):
        _use(bt, monkeypatch, tmp_path / "wide.md", LEDGER.replace(strength, strength.replace("5`", "50`")))
        assert not bt._heart_further_ordered(ROW), strength  # `0.750` is not `0.75`
        rounded = f"`--trend-strength {float(strength[18:-1]):.1f}`"
        _use(bt, monkeypatch, tmp_path / "round.md", LEDGER.replace(strength, rounded))
        assert not bt._heart_further_ordered(ROW), strength  # `0.8` is not `0.75`


@pytest.mark.parametrize("text", [
    "| R27 (other) | T1 | x |\n",
    f"| {R29_FIRST} | heart | x |\n",
    LEDGER + f"| {F_FIRST} | heart | y |\n",
    LEDGER.replace(F_FIRST, F_FIRST.replace("`--entry organ_scope`", "`--entry core_trim`")),
    LEDGER.replace(F_FIRST, F_FIRST.replace("`--trend-strength 0.75`", "`--trend-strength 0.8`")),
    LEDGER.replace(F_FIRST, F_FIRST.replace("`--trend-strength 1.25`", "`--trend-strength 1.2`")),
    LEDGER.replace(F_FIRST, F_FIRST.replace("`--captain-ruling R49`", "`--captain-ruling R29`")),
    LEDGER.replace(F_FIRST, F_FIRST.replace("ruling, no run", "awaiting a ruling"))])
def test_the_build_is_refused_on_a_ledger_without_the_order(bt, organ_env, ledger, tmp_path, text):
    ledger.write_text(text)
    for s in PAIR:
        with pytest.raises(SystemExit) as e:
            bt.main(["heart", *OTSF[s], "--out", str(tmp_path / "h.h5ad")])
        assert f"must hold exactly one row {ROW}" in str(e.value)
    assert organ_env == [] and not (tmp_path / "h.h5ad").exists()


def test_rows_r28s_and_r29s_lines_and_the_rows_line_are_apart(bt, organ_env, ledger, tmp_path):
    """On a ledger with the new row alone, rows R28's and R29's files are refused as before; on one with R28 and R29
    alone, the two are refused; the row never orders 0.5 or 1.0."""
    ledger.write_text(f"| {F_FIRST} | heart | x |\n")
    for s in (0.5, 1.0):
        with pytest.raises(SystemExit) as e:
            bt.main(["heart", *OTS[s], "--out", str(tmp_path / "h.h5ad")])
        assert "must hold exactly one row R29" in str(e.value)
        with pytest.raises(SystemExit) as e:
            bt.main(["heart", *TS[s], "--out", str(tmp_path / "h.h5ad")])
        assert "must hold exactly one row R28" in str(e.value)
        with pytest.raises(SystemExit) as e:
            bt.main(["heart", *ORG, "--trend-strength", f"{s:.1f}", "--z-mode", "raw", "--captain-ruling", ROW,
                     "--out", str(tmp_path / "h.h5ad")])
        assert "only row R29 orders it" in str(e.value)
    ledger.write_text(LEDGERS["full"].replace(f"| {F_FIRST} | heart | file D's two further files |\n", ""))
    for s in PAIR:
        with pytest.raises(SystemExit) as e:
            bt.main(["heart", *OTSF[s], "--out", str(tmp_path / "h.h5ad")])
        assert f"must hold exactly one row {ROW}" in str(e.value)
    assert organ_env == [] and not (tmp_path / "h.h5ad").exists()
