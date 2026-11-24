"""`experiments/builds/build_test.py heart --entry organ_scope --trend-strength 1.5 --z-mode raw --captain-ruling R52`:
file D's strength 1.5, file OTS15 (row R52 of experiments/harness/TRIES.md, the captain's order of 2026-11-24, ruling,
no run: his named exception on the heart's scope and size for this file, with H34's stop and the line's cap lifted for
this file only): the organ cut D's command with --expression copy replaced by row H28's strength-scaled trend step at
1.5, on OTS05's draw, organ-cut pool, count rule, size procedure and seed; and the pin on file D's two further
strengths, 0.75 and 1.25 only with --captain-ruling R49.

On the synthetic stages and synthetic test board (its limit never shown) of tests/test_heart_trend_build.py, E10.5
relabelled as tests/test_build_heart_organ_scope.py's `organ_stage` and that file's size stand-in, with a temporary
ledger in place of TRIES.md (synthetic first cells for rows R28, R29, R49 and R52): the log factors and the damping the
step hands momentum_shift at 1.5 are exactly 1.5 times those at strength 1 (row R29's OTS10); the step is row H28's own
at 1.5; the rows, coordinates, cell number, size field and, after halve.py at seed 1000001, the half's rows are OTS05's
and only the expression differs; the same bytes twice; the record names the option, the strength, the file, the row,
the order text, the row's first cell and the command; every refusal, before anything is read; the R52 check on
temporary ledgers and, read by code, on the committed ledger. Every heart path that existed before (rows R28's, R29's
and R49's files, files C and D, the default build and the other heart commands) gives the bytes and record it gave
before, and every refusal of an old path is word for word the one it gave before (pinned at main b42442e, before this
change); row R49's two values are refused under any other row, even one whose first cell orders them. Nothing is read
from data/ and nothing is scored."""

from __future__ import annotations

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
from test_heart_trend_line_f import _COMMANDS, _X, OTS_BEFORE, _coords, _recorder  # noqa: E402

ORG = ["--entry", "organ_scope"]
CORE = ["--entry", "core_trim"]
OTS15 = [*ORG, "--trend-strength", "1.5", "--z-mode", "raw", "--captain-ruling", "R52"]
PAIR = (0.75, 1.25)
VALUES = ("0.75", "1.25")
OTSF = {s: [*ORG, "--trend-strength", v, "--z-mode", "raw", "--captain-ruling", "R49"]
        for s, v in zip(PAIR, VALUES, strict=True)}
R49_FIRST = ("R49 (the captain's order of 2026-11-24, ruling, no run: a synthetic stand-in; the heart builder is to "
             "take `--entry organ_scope` with `--trend-strength 0.75` and `--trend-strength 1.25` only with "
             "`--captain-ruling R49`; not a try)")
R52_FIRST = ("R52 (the captain's order of 2026-11-24, ruling, no run: \"Yes to OTS15\", a synthetic stand-in; the "
             "heart builder is to take `--entry organ_scope` with `--trend-strength 1.5` only with "
             "`--captain-ruling R52`; not a try)")
R52_STRINGS = ("`--entry organ_scope`", "`--trend-strength 1.5`", "`--captain-ruling R52`")
LEDGER = (f"| R27 (other) | T1 | x |\n| {R28_FIRST} | heart | the line's two files |\n"
          f"| {R29_FIRST} | heart | file D's two files |\n| {R49_FIRST} | heart | file D's two further files |\n"
          f"| {R52_FIRST} | heart | file D's strength 1.5 |\n| H28 (row) | heart | x |\n")
R44_FIRST = R49_FIRST.replace("R49", "R44")  # row R49's form under another id: main built the two under it
LEDGERS = {
    "full": LEDGER,
    "no_r28": LEDGER.replace(f"| {R28_FIRST} | heart | the line's two files |\n", ""),
    "no_r29": LEDGER.replace(f"| {R29_FIRST} | heart | file D's two files |\n", ""),
    "no_r49": LEDGER.replace(f"| {R49_FIRST} | heart | file D's two further files |\n", ""),
    "no_r52": LEDGER.replace(f"| {R52_FIRST} | heart | file D's strength 1.5 |\n", ""),
    "r49_rounded": LEDGER.replace(R49_FIRST, R49_FIRST.replace("`--trend-strength 0.75`", "`--trend-strength 0.8`")),
    "r44": LEDGER + f"| {R44_FIRST} | heart | x |\n",
    # row R52's first cell also naming row R49's two, and a row of R52's form under another id
    "r52_wide": LEDGER.replace(R52_FIRST, R52_FIRST.replace(
        "with `--trend-strength 1.5`",
        "with `--trend-strength 0.75`, `--trend-strength 1.25` and `--trend-strength 1.5`")),
    "r53": LEDGER + f"| {R52_FIRST.replace('R52', 'R53')} | heart | x |\n",
}

# sha256 of (the .h5ad, its .build.json text, the record's `canon`) that build_test.py heart wrote for row R49's two
# commands (--captain-ruling R49) on `organ_env`'s synthetic stages and size stand-in, with LEDGER, computed with the
# builder of main at b42442e (before this change) by the same builds as below. A record may differ from its text only in
# the last bits of floats. Rows R28's and R29's files are pinned in TS_BEFORE and OTS_BEFORE, and every other heart
# command in BEFORE; all are rechecked here.
OTSF_BEFORE = {
    0.75: ("4fe9da095523e97a3326b08521d7ee1333630a153d59d59b73aff0eee9d247a1",
           "ed5fc01dc55477e8d984bb80b8a00c8386fd65869d4a9ae792a26279f0329839",
           "14cf5f138855a003ee052c52267ae34c3bd67cb4b5e729490795cddd7e5eac39"),
    1.25: ("3dfaedcaa54f5b5184d834904cb4d6b8cf404480c4a0b570538d69fede7c1153",
           "5f9add5a8882a40a38ae742ef9d1fbaf4b281994906915be0857fbede06b1344",
           "8c3150b475983ce570da44b04c950a00311faae97180c09bf78baea45907d551"),
}

# Every refusal of an old path, before anything is read: case -> (argv, ledger). None uses --trend-strength 1.5, the
# value this change adds; every other spelling of it (1.50, 1.5e0, ...) is an old path. The messages were computed with
# the builder of main at b42442e (before this change), and OLD_REFUSALS_SHA256 is the sha256 of
# json.dumps({case: message}, sort_keys=True) over all of them. OLD_MESSAGES spells out one case of each kind.
RULED = ("0.5", "1.0", "0.75", "1.25")  # the values rows R28, R29 and R49 order
OTHER_VALUES = ("1.50", "1.5e0", "1.500", "01.5", "+1.5", "1.5 ", " 1.5", "1,5", "15", "1.05", "1.55", "1.4", "1.6",
                "1.75", "2.0", "3.0", "0", "1", "0.50", "-1.0", "-1.5", "0.2", "0.8", "1.2", "0.250", ".25", "1.25e0",
                "0.125", "", "0.25", "0.750", "1.250", ".75", ".5", "1.")
BAD_IDS = ("r49", "R-F", "49", "R49 ", " R49", "R49;", "R", "RR49", "R49ab", "R4.9", "", "R49`")
FLAGS = ("--expression copy", "--trend", "--depth", "--herd", "--cfm")
VAL_BOARDS = ("T2:heart:val_extrap", "T2:heart:val_interp")
OLD_REFUSALS = {
    **{f"ruling_alone_{n}": (argv, "full") for n, argv in (
        ("d_R52", [*ORGAN, "--captain-ruling", "R52"]), ("d_R49", [*ORGAN, "--captain-ruling", "R49"]),
        ("org_R52", [*ORG, "--z-mode", "raw", "--captain-ruling", "R52"]),
        ("s1_trend_R52", ["--entry", "s1", "--trend", "--captain-ruling", "R52"]),
        ("c_R52", [*CORE, "--expression", "copy", "--z-mode", "raw", "--captain-ruling", "R52"]))},
    **{f"value_{v}_{e[1]}_{r}": ([*e, "--trend-strength", v, "--z-mode", "raw", "--captain-ruling", r], "full")
       for v in OTHER_VALUES
       for e, r in ((ORG, "R52"), (ORG, "R49"), (ORG, "R29"), (CORE, "R28"), (CORE, "R52"),
                    (["--entry", "s1"], "R52"))},
    **{f"r52_{v}_{e[1]}": ([*e, "--trend-strength", v, "--z-mode", "raw", "--captain-ruling", "R52"], "full")
       for v in RULED for e in (ORG, CORE)},
    # rows R28's and R29's lines
    **{f"zmode_{e[1]}_{s}": ([*e, "--trend-strength", s, "--z-mode", "recorded", "--captain-ruling", r], "full")
       for s in ("0.5", "1.0") for e, r in ((ORG, "R29"), (CORE, "R28"))},
    **{f"nozmode_{e[1]}_{s}": ([*e, "--trend-strength", s, "--captain-ruling", r], "full")
       for s in ("0.5", "1.0") for e, r in ((ORG, "R29"), (CORE, "R28"))},
    **{f"flag_{f[2:]}_{s}_{w}": ([*line[float(s)], *f.split()], "full")
       for f in FLAGS for s in ("0.5", "1.0") for w, line in (("d", OTS), ("c", TS))},
    **{f"entry_{e}_{s}_r29": (["--entry", e, "--trend-strength", s, "--z-mode", "raw", "--captain-ruling", "R29"],
                              "full")
       for e in ("s1", "s2", "s2_trim", "core_trim") for s in ("0.5", "1.0")},
    **{f"entry_{e}_{s}_r28": (["--entry", e, "--trend-strength", s, "--z-mode", "raw", "--captain-ruling", "R28"],
                              "full")
       for e in ("s1", "s2", "s2_trim", "organ_scope") for s in ("0.5", "1.0")},
    **{f"noruling_{e[1]}_{s}": ([*e, "--trend-strength", s, "--z-mode", "raw"], "full")
       for s in ("0.5", "1.0", "0.25") for e in (ORG, CORE)},
    **{f"ruling_{r}_{e[1]}_{s}": ([*e, "--trend-strength", s, "--z-mode", "raw", "--captain-ruling", r], "full")
       for s in ("0.5", "1.0") for e in (ORG, CORE) for r in ("R21", "R30", "r29", "R44", "R45", "R49", "R52")},
    **{f"board_{b}_{s}_{w}": ([*line[float(s)], "--board", b], "full")
       for b in VAL_BOARDS for s in ("0.5", "1.0") for w, line in (("d", OTS), ("c", TS))},
    **{f"ledger_no_r29_{s}": (OTS[float(s)], "no_r29") for s in ("0.5", "1.0")},
    **{f"ledger_no_r28_{s}": (TS[float(s)], "no_r28") for s in ("0.5", "1.0")},
    # row R49's line (the row given on the command line, R49 the one that orders it)
    **{f"f_noruling_{v}": ([*ORG, "--trend-strength", v, "--z-mode", "raw"], "full") for v in VALUES},
    **{f"f_id_{r}_{v}": ([*ORG, "--trend-strength", v, "--z-mode", "raw", "--captain-ruling", r], "full")
       for v in VALUES for r in BAD_IDS},
    **{f"f_row_{r}_{v}": ([*ORG, "--trend-strength", v, "--z-mode", "raw", "--captain-ruling", r], "full")
       for v in VALUES for r in ("R45", "R29", "R28", "R27", "R4", "R490", "R52", "R44", "R50")},
    **{f"f_entry_{e}_{v}": (["--entry", e, "--trend-strength", v, "--z-mode", "raw", "--captain-ruling", "R49"],
                            "full")
       for e in ("s1", "s2", "s2_trim", "core_trim") for v in VALUES},
    **{f"f_zmode_{v}": ([*ORG, "--trend-strength", v, "--z-mode", "recorded", "--captain-ruling", "R49"], "full")
       for v in VALUES},
    **{f"f_nozmode_{v}": ([*ORG, "--trend-strength", v, "--captain-ruling", "R49"], "full") for v in VALUES},
    **{f"f_flag_{f[2:]}_{v}": ([*OTSF[float(v)], *f.split()], "full") for f in FLAGS for v in VALUES},
    **{f"f_board_{b}_{v}": ([*OTSF[float(v)], "--board", b], "full") for b in VAL_BOARDS for v in VALUES},
    **{f"f_ledger_{led}_{v}": (OTSF[float(v)], led) for led in ("no_r49", "r49_rounded") for v in VALUES},
}
OLD_REFUSALS_N = 387
OLD_REFUSALS_SHA256 = "d5407a37251334397bcb26ff51e29091d9d38cbd874bc8df4748fdce1f3bbe5a"
OLD_MESSAGES = {
    "ruling_alone_d_R52": f"heart --captain-ruling goes with --trend-strength only ({_COMMANDS})",
    "value_1.50_organ_scope_R52": "--trend-strength 1.50 is refused: it takes 0.5 or 1.0 only",
    "value_1.5e0_organ_scope_R52": "--trend-strength 1.5e0 is refused: it takes 0.5 or 1.0 only",
    "value_+1.5_core_trim_R52": "--trend-strength +1.5 is refused: it takes 0.5 or 1.0 only",
    "r52_0.5_organ_scope": ("--entry organ_scope --trend-strength 0.5 is refused with --captain-ruling R52: only row "
                            "R29 orders it"),
    "r52_1.0_core_trim": "--trend-strength 1.0 is refused with --captain-ruling R52: only row R28 orders it",
    "r52_0.75_organ_scope": ("--entry organ_scope --trend-strength 0.75 is refused: experiments/harness/TRIES.md must "
                             "hold exactly one row R52 whose first cell is the captain's order (\"ruling, no run\") "
                             "naming `--entry organ_scope`, `--trend-strength 0.75`, `--trend-strength 1.25` and "
                             "`--captain-ruling R52`"),
    "f_noruling_1.25": ("--entry organ_scope --trend-strength 1.25 needs --captain-ruling <ROW>, the row of "
                        "experiments/harness/TRIES.md that is the captain's order for file D's two further strengths"),
    "f_id_r49_0.75": ("--entry organ_scope --trend-strength 0.75 is refused with --captain-ruling r49: <ROW> is a "
                      "ruling row's id, R and its number (e.g. R44)"),
    "f_row_R44_0.75": ("--entry organ_scope --trend-strength 0.75 is refused: experiments/harness/TRIES.md must hold "
                       "exactly one row R44 whose first cell is the captain's order (\"ruling, no run\") naming "
                       "`--entry organ_scope`, `--trend-strength 0.75`, `--trend-strength 1.25` and "
                       "`--captain-ruling R44`"),
    "f_ledger_no_r49_1.25": ("--entry organ_scope --trend-strength 1.25 is refused: experiments/harness/TRIES.md must "
                             "hold exactly one row R49 whose first cell is the captain's order (\"ruling, no run\") "
                             "naming `--entry organ_scope`, `--trend-strength 0.75`, `--trend-strength 1.25` and "
                             "`--captain-ruling R49`"),
    "f_entry_core_trim_0.75": ("--trend-strength 0.75 is refused except as heart --entry organ_scope "
                               "--trend-strength 0.75 --z-mode raw --captain-ruling <ROW>: --entry organ_scope and "
                               "--z-mode raw only, and with --expression copy, --trend, --depth, --herd and --cfm"),
}
# The default build's refusals (no --trend-strength): after the boards are read and before any stage is, as before.
# Pinned the same way at main b42442e.
OLD_DEFAULT_REFUSALS = {
    "org_bare": [*ORG], "org_raw": [*ORG, "--z-mode", "raw"], "org_copy": [*ORG, "--expression", "copy"],
    "org_copy_raw_trend": [*ORGAN, "--trend"], "org_copy_raw_depth": [*ORGAN, "--depth"],
    "org_copy_raw_cfm": [*ORGAN, "--cfm"], "core_bare": [*CORE], "core_copy": [*CORE, "--expression", "copy"],
    "s2_trend": ["--entry", "s2", "--trend"], "s2_copy": ["--entry", "s2", "--expression", "copy"],
    "s1_copy_trend": ["--entry", "s1", "--expression", "copy", "--trend"],
    "s2_trim_copy": ["--entry", "s2_trim", "--expression", "copy"],
}
OLD_DEFAULT_REFUSALS_SHA256 = "1c0d0b4aa673ec47c2fafab665dde29dc45e1974c6493f168f2de51fdf7ab2cd"


@pytest.fixture(scope="module")
def bt():
    sys.path.insert(0, str(REPO / "experiments" / "builds"))
    import build_test

    return build_test


@pytest.fixture(scope="module")
def h28():
    spec = importlib.util.spec_from_file_location("h28_for_heart_trend_line_g",
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


def _why(bt, tmp_path, *argv) -> str:
    out = tmp_path / ("h_dry.h5ad" if "--board" in argv else "h.h5ad")
    with pytest.raises(SystemExit) as e:
        bt.main(["heart", *argv, "--out", str(out)])
    assert not out.exists()
    return str(e.value)


# ------------------------------------------------------------------------------------------------ off by default
def test_the_option_is_off_by_default_and_the_old_tables_stand(bt):
    a = bt.build_parser().parse_args(["heart", *ORGAN, "--out", "runs/x/h.h5ad"])
    assert a.trend_strength is None and a.captain_ruling is None
    assert bt.HEART_OTS15_ROW == "R52" and bt.HEART_OTS15_CHOICE == "1.5" and bt.HEART_OTS15_STRENGTH == 1.5
    assert bt.HEART_OTS15_OPTION == "`--trend-strength 1.5`" and bt.HEART_OTS15_FILE == "OTS15"
    assert bt.HEART_OTS15_COMMAND == "heart " + " ".join(OTS15)
    assert bt.HEART_FURTHER_ROW == "R49"
    assert {s: c.format(row=bt.HEART_FURTHER_ROW) for s, c in bt.HEART_FURTHER_COMMANDS.items()} == {
        s: "heart " + " ".join(v) for s, v in OTSF.items()}
    # the old tables are as they were: 1.5 is in none of them
    assert bt.HEART_TS_CHOICES == ("0.5", "1.0") and bt.HEART_FURTHER_CHOICES == VALUES
    assert bt.HEART_TS_ROW == "R28" and bt.HEART_OTS_ROW == "R29"
    assert bt.HEART_TS_COMMANDS == {s: "heart " + " ".join(v) for s, v in TS.items()}
    assert bt.HEART_OTS_COMMANDS == {s: "heart " + " ".join(v) for s, v in OTS.items()}
    assert all(1.5 not in t for t in (bt.HEART_TS_FILES, bt.HEART_OTS_FILES, bt.HEART_FURTHER_FILES))
    assert bt.HEART_ORGAN_SCOPE == "heart " + " ".join(ORGAN)


# ------------------------------------------------------------------------------------------------ old paths unchanged
@pytest.mark.parametrize("name", sorted(BEFORE))
def test_every_heart_command_before_keeps_its_bytes_and_record(bt, own_repo, ledger, monkeypatch, tmp_path, name):
    install(bt, monkeypatch.setattr, [])
    if name.startswith("organ"):
        _use_stage(bt, monkeypatch, organ_stage(), [])
        _organ_size(bt, monkeypatch)
    out, rec, text = _build(bt, tmp_path, name, *ARGV[name])
    assert sha(out) == BEFORE[name][0]
    assert hashlib.sha256(text.encode()).hexdigest() == BEFORE[name][1] or canon(text) == BEFORE[name][2]
    assert "trend_strength" not in rec and "R52" not in text and "OTS15" not in text


@pytest.mark.parametrize("s", [0.5, 1.0])
def test_row_r28s_two_files_keep_their_bytes_and_record(bt, own_repo, ledger, monkeypatch, tmp_path, s):
    install(bt, monkeypatch.setattr, [])
    out, rec, text = _build(bt, tmp_path, f"ts{s:g}", *TS[s])
    assert sha(out) == TS_BEFORE[s][0]
    assert hashlib.sha256(text.encode()).hexdigest() == TS_BEFORE[s][1] or canon(text) == TS_BEFORE[s][2]
    assert rec["trend_strength"]["captain_ruling"] == "R28" and "R52" not in text


@pytest.mark.parametrize("s", [0.5, 1.0])
def test_row_r29s_two_files_keep_their_bytes_and_record(bt, organ_env, tmp_path, s):
    out, rec, text = _build(bt, tmp_path, f"ots{s:g}", *OTS[s])
    assert sha(out) == OTS_BEFORE[s][0]
    assert hashlib.sha256(text.encode()).hexdigest() == OTS_BEFORE[s][1] or canon(text) == OTS_BEFORE[s][2]
    assert rec["trend_strength"]["captain_ruling"] == "R29" and "R52" not in text


@pytest.mark.parametrize("s", PAIR)
def test_row_r49s_two_files_keep_their_bytes_and_record(bt, organ_env, tmp_path, s):
    out, rec, text = _build(bt, tmp_path, f"f{s:g}", *OTSF[s])
    assert OTSF_BEFORE[s] is not None
    assert sha(out) == OTSF_BEFORE[s][0]
    assert hashlib.sha256(text.encode()).hexdigest() == OTSF_BEFORE[s][1] or canon(text) == OTSF_BEFORE[s][2]
    ts = rec["trend_strength"]
    assert ts["captain_ruling"] == "R49" and ts["order_row"] == R49_FIRST and ts["file"] == {0.75: "OTS075",
                                                                                            1.25: "OTS125"}[s]
    assert "R52" not in text and "OTS15" not in text


def test_every_old_refusal_is_word_for_word_as_before(bt, organ_env, ledger, tmp_path):
    msgs = {}
    for name, (argv, led) in OLD_REFUSALS.items():
        ledger.write_text(LEDGERS[led])
        msgs[name] = _why(bt, tmp_path, *argv)
        assert organ_env == [], name
    assert not any(a[i:i + 2] == ["--trend-strength", "1.5"] for a, _ in OLD_REFUSALS.values() for i in range(len(a)))
    assert OLD_REFUSALS_N is not None and len(msgs) == len(OLD_REFUSALS) == OLD_REFUSALS_N
    for name, msg in OLD_MESSAGES.items():
        assert msgs[name] == msg, name
    assert hashlib.sha256(json.dumps(msgs, sort_keys=True).encode()).hexdigest() == OLD_REFUSALS_SHA256


def test_every_old_refusal_of_the_default_build_is_word_for_word_as_before(bt, organ_env, tmp_path):
    msgs = {}
    for name, argv in OLD_DEFAULT_REFUSALS.items():
        organ_env.clear()
        msgs[name] = _why(bt, tmp_path, *argv)
        assert organ_env == ["boards"], name  # the boards only, as before: no stage is read
    assert msgs["org_bare"].startswith(f"--entry organ_scope is refused except as {bt.HEART_ORGAN_SCOPE} ")
    assert hashlib.sha256(json.dumps(msgs, sort_keys=True).encode()).hexdigest() == OLD_DEFAULT_REFUSALS_SHA256


# ------------------------------------------------------------------------------------------------ the step
def test_the_factor_at_1_5_is_1_5_times_the_strength_one_factor(bt, organ_env, monkeypatch, tmp_path):
    """Through the builder's own commands: the log-factor vectors and the damping the step hands momentum_shift at 1.5
    are exactly 1.5 times those at strength 1 (row R29's OTS10; no tolerance), over the same lineages, and the
    strength-1 vectors are vec.harness.b0_trend_factors' own on D's context."""
    from vec import harness

    calls = _recorder(bt, monkeypatch)
    _build(bt, tmp_path, "ots10", *OTS[1.0])
    _build(bt, tmp_path, "ots15", *OTS15)
    assert len(calls) == 2  # one step per build
    (d1, m1), (d, m) = calls
    assert d1 == harness.B0_DAMPING and m1 and max(float(np.abs(v).max()) for v in m1.values()) > 0.01
    assert d == 1.5 * d1 == 0.75
    assert sorted(m) == sorted(m1)
    for k in m1:
        assert m[k].dtype == m1[k].dtype == np.float64
        np.testing.assert_array_equal(m[k], 1.5 * m1[k])
        assert np.array_equal(m1[k] != 0, m[k] != 0) and (not np.any(m1[k]) or not np.array_equal(m[k], m1[k]))
    ctx = _ctx_d(bt)
    pool, rows, _, _ = harness.draw_copy(ctx)
    merged = np.asarray(pool.obs["merged"]).astype(str)[rows]
    types = sorted(set(np.where(merged == "Untracked", "Untracked@end", merged)))
    logf, _ = harness.b0_trend_factors(ctx, list(pool.var_names), types)
    assert sorted(logf) == sorted(m1)
    for k in m1:
        np.testing.assert_array_equal(m1[k], logf[k])


def test_the_step_is_row_h28s_at_1_5_on_ds_context(bt, h28):
    ctx = _ctx_d(bt)
    st, rec = bt._trend_strength_step(ctx, 1.5)
    ref, _, _, ref_rec = h28.trend_expression(ctx, 1.5)
    np.testing.assert_array_equal(st.X, ref.X)
    np.testing.assert_array_equal(st.coords, ref.coords)
    assert rec["budget_days"] == pytest.approx(0.5) and ref_rec["effective_budget_days"] == pytest.approx(0.75)


def test_rows_and_coordinates_are_ots05s_and_only_the_expression_differs(bt, organ_env, tmp_path):
    """Rows, coordinates, cell number, size field, D's own record fields and, after halve.py at seed 1000001, the
    half's rows and coordinates are OTS05's (and so D's); the expression differs from OTS05's, OTS10's, OTS075's and
    OTS125's."""
    sys.path.insert(0, str(REPO / "experiments" / "builds"))
    import halve

    o, o_rec, _ = _build(bt, tmp_path, "ots05", *OTS[0.5])
    t, t_rec, _ = _build(bt, tmp_path, "ots15", *OTS15)
    A, B = ad.read_h5ad(o), ad.read_h5ad(t)
    assert A.n_obs == B.n_obs and list(A.obs_names) == list(B.obs_names) and list(A.var_names) == list(B.var_names)
    np.testing.assert_array_equal(_coords(o), _coords(t))
    for k in ("n_cells", "seed", "count_decision", "trim", "size", "z_mode", "board", "dry_run", "depth_step",
              "scope_option"):
        assert t_rec[k] == o_rec[k], k
    assert halve.DEFAULT_SEED == 1000001
    ho, ht = tmp_path / "ots05_half.h5ad", tmp_path / "ots15_half.h5ad"
    halve.halve(str(o), str(ho))
    halve.halve(str(t), str(ht))
    HA, HB = ad.read_h5ad(ho), ad.read_h5ad(ht)
    assert list(HA.obs_names) == list(HB.obs_names) and HA.n_obs == HB.n_obs
    np.testing.assert_array_equal(_coords(ho), _coords(ht))
    assert not np.array_equal(_X(ho), _X(ht))
    others = {"0.5": o, "1.0": _build(bt, tmp_path, "ots10", *OTS[1.0])[0],
              **{f"{s:g}": _build(bt, tmp_path, f"f{s:g}", *OTSF[s])[0] for s in PAIR}}
    X15 = _X(t)
    for name, path in others.items():
        assert not np.array_equal(X15, _X(path)), name


def test_the_same_bytes_twice(bt, organ_env, tmp_path):
    a, _, ta = _build(bt, tmp_path, "a", *OTS15)
    b, _, tb = _build(bt, tmp_path, "b", *OTS15)
    assert a.read_bytes() == b.read_bytes() and ta == tb


def test_the_record_names_option_strength_file_row_order_and_command(bt, organ_env, tmp_path, capsys):
    _, rec, text = _build(bt, tmp_path, "t", *OTS15)
    ts = rec["trend_strength"]
    assert ts["option"] == "--trend-strength" and ts["strength"] == 1.5 and ts["captain_ruling"] == "R52"
    assert ts["file"] == "OTS15" and '"strength": 1.5,' in text
    assert ts["order"] == bt.HEART_OTS15_ORDER
    assert ts["order"].startswith("experiments/harness/TRIES.md row R52, the captain's order of 2026-11-24")
    assert "Yes to OTS15" in ts["order"] and "OTS15" in ts["order"] and "nothing else changed" in ts["order"]
    assert ts["order_row"] == R52_FIRST  # the row's first cell as the ledger holds it
    assert ts["command"] == "heart " + " ".join(OTS15) == bt.HEART_OTS15_COMMAND
    assert ts["step"] == bt.HEART_OTS15_STEP and "D's draw" in ts["step"]
    assert ts["budget_days"] == pytest.approx(0.5) and ts["test_budget_days"] == pytest.approx(0.75)
    assert ts["trend"]["rule"] == "vec.harness.B0_TREND_RULE"
    assert rec["scope_option"]["command"] == bt.HEART_ORGAN_SCOPE and "R21" in rec["scope_option"]["ruling"]
    assert rec["entry"] == "trend strength 1.5 organ_scope (organ)" and "expression" not in rec
    assert rec["z_mode"] == "raw" and rec["n_cells"] == "max_cells" and "trend" not in rec
    assert "R28" not in text and "R29" not in text and "R49" not in text
    assert str(SECRET_MAX) not in text + capsys.readouterr().out


def test_the_reads_are_ots05s(bt, organ_env, tmp_path):
    """OTS15 reads what OTS05 reads, in the same order: the boards, E10.5, E9.5 and the three earlier stages."""
    _build(bt, tmp_path, "ots05", *OTS[0.5])
    want = list(organ_env)
    assert want == ["boards", "E10.5.h5ad", "E9.5.h5ad", "E8.75.h5ad", "E8.5.h5ad", "E8.25_late.h5ad", "boards"]
    organ_env.clear()
    _build(bt, tmp_path, "ots15", *OTS15)
    assert organ_env == want


# ------------------------------------------------------------------------------------------------ refusals
NEW_REFUSED = [
    [*ORG, "--trend-strength", "1.5", "--z-mode", "raw"],  # no row
    *([*ORG, "--trend-strength", "1.5", "--z-mode", "raw", "--captain-ruling", r]  # a row that does not order it
      for r in ("R29", "R28", "R49", "R44", "R21", "R5", "R520", "R53", "r52", "R52 ", " R52", "R52a", "R-F", "")),
    *(["--entry", e, "--trend-strength", "1.5", "--z-mode", "raw", "--captain-ruling", "R52"]
      for e in ("s1", "s2", "s2_trim", "core_trim")),
    *(["--entry", e, "--trend-strength", "1.5", "--z-mode", "raw", "--captain-ruling", r]
      for e in ("s1", "core_trim") for r in ("R28", "R29")),
    [*ORG, "--trend-strength", "1.5", "--captain-ruling", "R52"],
    [*ORG, "--trend-strength", "1.5", "--z-mode", "recorded", "--captain-ruling", "R52"],
    *([*OTS15, *f.split()] for f in FLAGS),
]


@pytest.mark.parametrize("argv", NEW_REFUSED)
def test_every_other_use_of_1_5_is_refused_before_anything_is_read(bt, organ_env, tmp_path, argv):
    msg = _why(bt, tmp_path, *argv)
    assert "refused" in msg or "needs --captain-ruling R52" in msg
    assert organ_env == []


def test_the_new_values_refusals_say_why(bt, organ_env, ledger, tmp_path):
    cmd = "heart " + " ".join(OTS15)
    assert _why(bt, tmp_path, *ORG, "--trend-strength", "1.5", "--z-mode", "raw") == (
        "--entry organ_scope --trend-strength 1.5 is the captain's order of row R52 alone and needs --captain-ruling "
        "R52")
    for r in ("R29", "R49", "R44", "R53", "r52", "R52a", ""):
        assert _why(bt, tmp_path, *ORG, "--trend-strength", "1.5", "--z-mode", "raw", "--captain-ruling", r) == (
            f"--entry organ_scope --trend-strength 1.5 is refused with --captain-ruling {r}: only row R52 orders it")
    for argv in ([*CORE, "--trend-strength", "1.5", "--z-mode", "raw", "--captain-ruling", "R52"],
                 [*ORG, "--trend-strength", "1.5", "--captain-ruling", "R52"], [*OTS15, "--expression", "copy"],
                 [*OTS15, "--trend"]):
        assert _why(bt, tmp_path, *argv) == (
            f"--trend-strength 1.5 is refused except as {cmd}: --entry organ_scope and --z-mode raw only, and with "
            "--expression copy, --trend, --depth, --herd and --cfm")
    assert _why(bt, tmp_path, *OTS15, "--board", "T2:heart:val_extrap") == (
        "--trend-strength 1.5 is refused on T2:heart:val_extrap: the heart search line builds on the test board only "
        "(no dry run)")
    ledger.write_text(LEDGERS["no_r52"])
    assert _why(bt, tmp_path, *OTS15) == (
        "--entry organ_scope --trend-strength 1.5 is refused: experiments/harness/TRIES.md must hold exactly one row "
        "R52 whose first cell is the captain's order (\"ruling, no run\") naming `--entry organ_scope`, "
        "`--trend-strength 1.5` and `--captain-ruling R52`")
    # every other spelling of 1.5 is another value, refused as before, under row R52 too
    ledger.write_text(LEDGER)
    for v in ("1.50", "1.5e0", "1.500", "+1.5", "1.5 "):
        assert _why(bt, tmp_path, *ORG, "--trend-strength", v, "--z-mode", "raw", "--captain-ruling", "R52") == (
            f"--trend-strength {v} is refused: it takes 0.5 or 1.0 only")
    assert organ_env == []


def test_row_r49s_two_are_refused_under_any_other_row_and_r52s_value_under_any_other(bt, organ_env, ledger,
                                                                                    monkeypatch, tmp_path):
    """The pin: a row whose first cell is the captain's order for the two further strengths in row R49's form (main
    built 0.75 and 1.25 under it) passes the row check but is refused, after it, so every refusal of the line before
    the pin stands word for word; so is row R52 when its first cell also names the two. 1.5 is refused under a row of
    R52's form with another id, before the ledger is read for it."""
    ledger.write_text(LEDGERS["r44"])
    assert bt._heart_further_ordered("R44") and bt._heart_further_ordered("R49")
    for v in VALUES:
        assert _why(bt, tmp_path, *ORG, "--trend-strength", v, "--z-mode", "raw", "--captain-ruling", "R44") == (
            f"--entry organ_scope --trend-strength {v} is refused with --captain-ruling R44: only row R49 orders it")
    ledger.write_text(LEDGERS["r52_wide"])
    assert bt._heart_further_ordered("R52") and bt._heart_ots15_ordered()
    for v in VALUES:
        assert _why(bt, tmp_path, *ORG, "--trend-strength", v, "--z-mode", "raw", "--captain-ruling", "R52") == (
            f"--entry organ_scope --trend-strength {v} is refused with --captain-ruling R52: only row R49 orders it")
    ledger.write_text(LEDGERS["r53"])
    assert _why(bt, tmp_path, *ORG, "--trend-strength", "1.5", "--z-mode", "raw", "--captain-ruling", "R53") == (
        "--entry organ_scope --trend-strength 1.5 is refused with --captain-ruling R53: only row R52 orders it")
    assert organ_env == []


@pytest.mark.parametrize("board", ["T2:heart:val_extrap", "T2:heart:val_interp", "T2:embryo:test_interp"])
def test_any_other_board_is_refused_before_anything_is_read(bt, organ_env, tmp_path, board):
    msg = _why(bt, tmp_path, *OTS15, "--board", board)
    # a heart validation board meets the line's own refusal; another task's board the builder's earlier one
    assert "test board only" in msg or (board.startswith("T2:embryo") and "takes only" in msg)
    assert organ_env == []


def test_the_test_boards_own_name_is_its_default(bt, organ_env, tmp_path):
    out, rec, _ = _build(bt, tmp_path, "named", *OTS15, "--board", "T2:heart:test_extrap")
    default, _, _ = _build(bt, tmp_path, "default", *OTS15)
    assert out.read_bytes() == default.read_bytes() and rec["dry_run"] is False


# ------------------------------------------------------------------------------------------------ the R52 check
def _use(bt, monkeypatch, path: Path, text: str) -> None:
    path.write_text(text)
    monkeypatch.setattr(bt, "TRIES_FILE", path)


def test_the_r52_check_passes_on_a_ledger_of_the_required_form(bt, monkeypatch, tmp_path):
    _use(bt, monkeypatch, tmp_path / "ok.md", LEDGER)
    assert bt._heart_ots15_ordered() and bt._heart_further_ordered("R49")
    assert bt._heart_d_line_ordered() and bt._heart_line_ordered()
    _use(bt, monkeypatch, tmp_path / "only.md", f"| {R52_FIRST} | heart | x |\n")
    assert bt._heart_ots15_ordered() and not bt._heart_further_ordered("R49")
    assert not bt._heart_d_line_ordered() and not bt._heart_line_ordered()
    _use(bt, monkeypatch, tmp_path / "r49.md", f"| {R49_FIRST} | heart | x |\n")
    assert not bt._heart_ots15_ordered() and bt._heart_further_ordered("R49")


@pytest.mark.parametrize("drop", [*R52_STRINGS, "ruling, no run", "the captain's order"])
def test_the_r52_check_fails_without_each_required_string(bt, monkeypatch, tmp_path, drop):
    first = R52_FIRST.replace(drop, "something else")
    assert first != R52_FIRST
    _use(bt, monkeypatch, tmp_path / "t.md", LEDGER.replace(R52_FIRST, first))
    assert not bt._heart_ots15_ordered() and bt._heart_further_ordered("R49") and bt._heart_d_line_ordered()


@pytest.mark.parametrize("moved", R52_STRINGS)
def test_the_r52_check_needs_the_strings_in_the_first_cell(bt, monkeypatch, tmp_path, moved):
    first = R52_FIRST.replace(moved, "something else")
    _use(bt, monkeypatch, tmp_path / "t.md", f"| {first} | heart | {moved} | {' '.join(R52_STRINGS)} |\n")
    assert not bt._heart_ots15_ordered()


def test_the_r52_check_fails_with_two_rows_another_rows_number_or_another_form(bt, monkeypatch, tmp_path):
    _use(bt, monkeypatch, tmp_path / "two.md", LEDGER + f"| {R52_FIRST} | heart | y |\n")
    assert not bt._heart_ots15_ordered()
    _use(bt, monkeypatch, tmp_path / "other.md", LEDGER.replace("| R52 (", "| R53 ("))
    assert not bt._heart_ots15_ordered()
    _use(bt, monkeypatch, tmp_path / "longer.md", LEDGER.replace(R52_FIRST, R52_FIRST.replace("R52", "R520")))
    assert not bt._heart_ots15_ordered()
    _use(bt, monkeypatch, tmp_path / "ruled.md", LEDGER.replace("`--captain-ruling R52`", "`--captain-ruling R49`"))
    assert not bt._heart_ots15_ordered()
    for wrong in ("`--trend-strength 1.50`", "`--trend-strength 1.5e0`", "`--trend-strength 1.25`",
                  "`--trend-strength 15`", "--trend-strength 1.5"):
        _use(bt, monkeypatch, tmp_path / "wrong.md", LEDGER.replace("`--trend-strength 1.5`", wrong))
        assert not bt._heart_ots15_ordered(), wrong
    _use(bt, monkeypatch, tmp_path / "kind.md", LEDGER.replace("R52 (the captain's order", "R52 (the captain's words"))
    assert not bt._heart_ots15_ordered()
    _use(bt, monkeypatch, tmp_path / "none.md", "| R27 (x) | T1 | x |\n")
    assert not bt._heart_ots15_ordered()


@pytest.mark.parametrize("text", [
    "| R27 (other) | T1 | x |\n", LEDGERS["no_r52"], LEDGER + f"| {R52_FIRST} | heart | y |\n",
    LEDGER.replace(R52_FIRST, R52_FIRST.replace("`--entry organ_scope`", "`--entry core_trim`")),
    LEDGER.replace(R52_FIRST, R52_FIRST.replace("`--trend-strength 1.5`", "`--trend-strength 1.50`")),
    LEDGER.replace(R52_FIRST, R52_FIRST.replace("`--captain-ruling R52`", "`--captain-ruling R49`")),
    LEDGER.replace(R52_FIRST, R52_FIRST.replace("ruling, no run", "awaiting a ruling")),
    LEDGER.replace(R52_FIRST, R52_FIRST.replace("the captain's order", "the captain's named exception")),
    LEDGER.replace("| R52 (", "| R53 (")])
def test_the_build_is_refused_on_a_ledger_without_the_order(bt, organ_env, ledger, tmp_path, text):
    ledger.write_text(text)
    assert "must hold exactly one row R52" in _why(bt, tmp_path, *OTS15)
    assert organ_env == []


def test_rows_r28s_r29s_and_r49s_lines_are_refused_on_a_ledger_with_only_row_r52(bt, organ_env, ledger, tmp_path):
    ledger.write_text(f"| {R52_FIRST} | heart | x |\n")
    for s in (0.5, 1.0):
        assert "must hold exactly one row R29" in _why(bt, tmp_path, *OTS[s])
        assert "must hold exactly one row R28" in _why(bt, tmp_path, *TS[s])
    for s in PAIR:
        assert "must hold exactly one row R49" in _why(bt, tmp_path, *OTSF[s])
    assert organ_env == []


def test_the_committed_ledger_holds_the_orders_the_launcher_builds_under(bt):
    """experiments/harness/TRIES.md as committed, read by code only: row R52's first cell is the captain's order in the
    form the builder takes, and so are rows R29's and R49's, whose files are the launcher's identity steps."""
    assert bt.TRIES_FILE == REPO / "experiments" / "harness" / "TRIES.md"
    assert bt._heart_ots15_ordered()
    assert bt._heart_further_ordered(bt.HEART_FURTHER_ROW) and bt._heart_d_line_ordered()
