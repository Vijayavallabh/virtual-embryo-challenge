"""`experiments/builds/build_test.py heart --entry organ_scope --trend-strength <S> --trend-form category --z-mode raw
--captain-ruling R74`: row R74's four files (experiments/harness/TRIES.md row R74, the captain's order of 2026-11-25,
ruling, no run), OTK125's command (row R55's category file) with one change, the category form's strength S, taken from
the builder's one literal table HEART_R74_TABLE in the row's order: OTK10 (1.0), OTK15 (1.5), OTK075 (0.75) and OTK175
(1.75).

On the synthetic stages and synthetic test board (its limit never shown) of tests/test_heart_trend_line_h.py, with a
temporary ledger holding row R55's and row R74's orders: the table and its names; the strength trap of the round's
challenge (HEART_REF3_challenge.md section 6 (a)): each file's step is the category form at the table's strength, never
R55's constant 1.25, so its bytes are not OTK125's and the four differ from each other; its rows, coordinates, cell
number, size field and half rows at seed 1000001 are OTK125's; the same bytes twice; the record names the option string,
the strength string verbatim, the file, the row, its order, the row's first cell and the command, and holds no cell
number; every other use is refused before anything is read (near-miss strengths, forms, rows, entries, options, seeds
and boards); the row check on temporary ledgers and, read by code, on the committed ledger. Every earlier heart command
builds as before: R55's four forms keep the bytes and records pinned here at main 5622ff3 (before this change), no
earlier command reaches row R74's code, and each earlier refusal meets --captain-ruling R74 as it met any other row.
Nothing is read from data/ and nothing is scored."""

from __future__ import annotations

import hashlib
import json
import re
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from test_build_heart_iso_trim import _build, canon  # noqa: E402
from test_build_heart_organ_scope import ORGAN, _use_stage, organ_stage  # noqa: E402
from test_heart_trend_build import SECRET_MAX, STAGES, install, sha  # noqa: E402
from test_heart_trend_line_d import OTS, _organ_size  # noqa: E402
from test_heart_trend_line_f import _X, _coords  # noqa: E402
from test_heart_trend_line_g import OTS15, OTSF  # noqa: E402
from test_heart_trend_line_h import (  # noqa: E402
    FORMS,
    LEDGER_TF,
    ORG,
    TF,
    _reads,
    _why,
    relabel,
    tf_context,
    tf_organ_stage,
)

# row R74's table as the row names it, in its order: file -> (option string, strength string, file stem)
TABLE = {"OTK10": ("--trend-strength 1.0 --trend-form category", "1.0", "T2heart_organ_ts10_tk"),
         "OTK15": ("--trend-strength 1.5 --trend-form category", "1.5", "T2heart_organ_ts15_tk"),
         "OTK075": ("--trend-strength 0.75 --trend-form category", "0.75", "T2heart_organ_ts075_tk"),
         "OTK175": ("--trend-strength 1.75 --trend-form category", "1.75", "T2heart_organ_ts175_tk")}
NAMES = list(TABLE)
R74 = {n: [*ORG, *TABLE[n][0].split(), "--z-mode", "raw", "--captain-ruling", "R74"] for n in NAMES}
R74_FIRST = ("R74 (the captain's order of 2026-11-25, ruling, no run: \"Yes to OTK10, OTK15, OTK075 and OTK175\", a "
             "synthetic stand-in; the heart builder is to take on `build_test.py heart` exactly these strings, each "
             "only with `--captain-ruling R74`: `--entry organ_scope`, `--z-mode raw`, `--trend-strength 1.0 "
             "--trend-form category` (OTK10), `--trend-strength 1.5 --trend-form category` (OTK15), `--trend-strength "
             "0.75 --trend-form category` (OTK075), `--trend-strength 1.75 --trend-form category` (OTK175), "
             "`--captain-ruling R74`; each file then halved by `halve.py ... --seed 1000001`; not a try)")
R74_ROW = f"| {R74_FIRST} | heart | the refinement round's four files |\n"
LEDGER_R74 = LEDGER_TF + R74_ROW
# sha256 of (the .h5ad, its .build.json text, the record's `canon`) that build_test.py heart wrote for row R55's four
# forms on tests/test_heart_trend_line_h.py's trend-form environment, computed with the builder of main at 5622ff3
# (before row R74's code) by the same builds as below
R55_BEFORE = {
    "category": ("03b29917b11a4a4665d6b68e61f31264fb6de9962317f70357ef0bafb11726dc",
                 "2edbbd78ca29c9a3389059b2e5bcf8a52c819f83fe6ca6da7e2190c33120013a",
                 "6b5cca64b84312456c5889cc226ba128f367e76b5099394156711699b34f51b0"),
    "contrast": ("24516055217778e02af7b778865e37f56797c842030943838f7de14623fc3861",
                 "1b74a467a90c4c399c78438087480ad3cf5516b00527eb7ef430b163a3f464f6",
                 "cc911fb4bae1aa6aba980f96a0822a9aed7424466da64111bf490e08e855e3ad"),
    "whole": ("f45d0441fda5f344c6085f798e3ac0e39c1ee6cccf739b4a0c3fe4d27a10efe8",
              "14b80b597c8461f155deaba6da5a482da4a5688c68969ce380096183d94d3c3c",
              "bc0846aad040f826a6fe2133d847a33f0c3fa94ff1b0e276d723e52196869460"),
    "portion": ("218578c39d06456b010505d3232e4c48ed604a612f965dfec8d08d2d34d54b6c",
                "7c5b0986b3825c08d41c1dfd32bdec1faafdd7d5964e0532c68f0e27acb11c72",
                "94ebf30fa95f160a5234f35d54e9d0f2330244372762886770777d9cce5213cb"),
}


@pytest.fixture(scope="module")
def bt():
    sys.path.insert(0, str(REPO / "experiments" / "builds"))
    import build_test

    return build_test


@pytest.fixture(scope="module")
def tf():
    from vec import trend_forms

    return trend_forms


@pytest.fixture
def ledger(bt, monkeypatch, tmp_path):
    """The temporary ledger of tests/test_heart_trend_line_h.py with row R74's order added."""
    path = tmp_path / "TRIES.md"
    path.write_text(LEDGER_R74)
    monkeypatch.setattr(bt, "TRIES_FILE", path)
    return path


@pytest.fixture
def organ_env(bt, monkeypatch, own_repo, ledger):
    """tests/test_heart_trend_line_h.py's organ environment (the earlier lines'), with rows R55 and R74 in the
    ledger."""
    reads: list[str] = []
    install(bt, monkeypatch.setattr, reads)
    _use_stage(bt, monkeypatch, organ_stage(), reads)
    _organ_size(bt, monkeypatch)
    more = _reads(bt, monkeypatch)
    return reads, more


@pytest.fixture
def tf_env(bt, monkeypatch, own_repo, ledger):
    """tests/test_heart_trend_line_h.py's trend-form environment: install()'s board and stages, every stage relabelled,
    E10.5 as `tf_organ_stage`, the organ size stand-in; the list records every stage and board read, in order."""
    reads: list[str] = []
    install(bt, monkeypatch.setattr, reads)
    load = bt.load_stage

    def wrapped(p, **kw):
        name = Path(p).name
        if name == "E10.5.h5ad":
            reads.append(name)
            return tf_organ_stage()
        return relabel(load(p, **kw), STAGES[name][2])

    monkeypatch.setattr(bt, "load_stage", wrapped)
    _organ_size(bt, monkeypatch)
    boards = bt.load_boards

    def load_boards():
        reads.append("boards")
        return boards()

    monkeypatch.setattr(bt, "load_boards", load_boards)
    return reads


@pytest.fixture(scope="module")
def cw(bt):
    return tf_context(bt)


def _use(bt, monkeypatch, path: Path, text: str) -> None:
    path.write_text(text)
    monkeypatch.setattr(bt, "TRIES_FILE", path)


def _spy(monkeypatch) -> list:
    """Every (strength, form) vec.trend_forms.form_step is called with, the call passed through."""
    from vec import trend_forms

    seen: list = []
    own = trend_forms.form_step

    def spy(ctx, strength, form=None, whole=None, *a, **k):
        seen.append((strength, form))
        return own(ctx, strength, form, whole, *a, **k)
    monkeypatch.setattr(trend_forms, "form_step", spy)
    return seen


# ------------------------------------------------------------------------------------------------ the table
def test_the_table_is_the_rows_strings_in_its_order(bt):
    assert bt.HEART_R74_ROW == "R74" and bt.HEART_R74_FORM == "category" and bt.HEART_R74_BASE_FILE == "OTK125"
    assert [tuple(p) for p in bt.HEART_R74_TABLE] == [(n, v[0]) for n, v in TABLE.items()]
    assert bt.HEART_R74_FILES == {n: (v[0], v[1]) for n, v in TABLE.items()}
    assert bt.HEART_R74_STEMS == {n: v[2] for n, v in TABLE.items()}
    assert bt.HEART_R74_BY_OPTION == {v[0]: n for n, v in TABLE.items()}
    assert bt.HEART_R74_NAMED == ("--entry organ_scope", "--z-mode raw", "--captain-ruling R74")
    for n in NAMES:
        assert bt.HEART_R74_COMMAND.format(option=TABLE[n][0]) == "heart " + " ".join(R74[n])
        # each string, whole and backticked, in the row, followed by its file's name
        assert f"`{TABLE[n][0]}` ({n})" in R74_FIRST
    strengths = [float(v[1]) for v in TABLE.values()]
    assert len(set(strengths)) == 4 and bt.HEART_TF_STRENGTH not in strengths
    assert all(v[1] != bt.HEART_TF_STRENGTH_CHOICE for v in TABLE.values())
    # R55's constants as they were: its strength, row and commands
    assert bt.HEART_TF_STRENGTH == 1.25 and bt.HEART_TF_STRENGTH_CHOICE == "1.25" and bt.HEART_TF_ROW == "R55"
    assert {f: bt.HEART_TF_COMMANDS[f] for f in FORMS} == {f: "heart " + " ".join(TF[f]) for f in FORMS}
    assert bt.HEART_TF_FILES["category"] == "OTK125" and bt.HEART_TF_TAGS["category"] == "tk"


def test_the_option_is_off_by_default(bt):
    a = bt.build_parser().parse_args(["heart", *ORGAN, "--out", "runs/x/h.h5ad"])
    assert a.trend_form is None and a.captain_ruling is None and a.trend_strength is None
    assert bt._r74_file(a) is None
    a = bt.build_parser().parse_args(["heart", *TF["category"], "--out", "runs/x/h.h5ad"])
    assert bt._r74_file(a) is None  # R55's file: not row R74's


@pytest.mark.parametrize("name, option", [
    ("OTK10", "--trend-strength 1.5 --trend-form category"),  # another file's string
    ("OTK15", "--trend-strength 1.0 --trend-form category"),
    ("OTK1", "--trend-strength 1.0 --trend-form category"), ("OTK100", "--trend-strength 1.0 --trend-form category"),
    ("OTK10", "--trend-strength 1.00 --trend-form category"), ("OTK10", "--trend-strength 1 --trend-form category"),
    ("OTK75", "--trend-strength 0.75 --trend-form category"), ("OTK075", "--trend-strength .75 --trend-form category"),
    ("OTK125", "--trend-strength 1.25 --trend-form category"),  # R55's own string is not row R74's
    ("OTK10", "--trend-strength 1.0 --trend-form contrast"), ("OTK10", "--trend-strength 1.0 --trend-form Category"),
    ("OTK10", "--trend-form category --trend-strength 1.0"), ("OTK10", "--trend-strength 1.0  --trend-form category"),
    ("OTK10", " --trend-strength 1.0 --trend-form category"), ("OTK10", "--trend-strength 1.0 --trend-form category "),
    ("OTK10", "--trend-strength 1.0 --trend-form category --z-mode raw"), ("otk10", "--trend-strength 1.0 "
                                                                            "--trend-form category"),
    ("OTX10", "--trend-strength 1.0 --trend-form category"),
    ("OTK10.h5ad", "--trend-strength 1.0 --trend-form category"),
])
def test_a_name_must_match_its_option_string(bt, name, option):
    with pytest.raises(ValueError, match="does not match its option string"):
        bt._r74_strength(name, option)


def test_each_row_pair_gives_its_strength_string(bt):
    assert {n: bt._r74_strength(n, v[0]) for n, v in TABLE.items()} == {n: v[1] for n, v in TABLE.items()}


# ------------------------------------------------------------------------------------------------ the strength trap
@pytest.mark.parametrize("name", NAMES)
def test_each_file_steps_at_the_tables_strength_never_r55s(bt, tf_env, monkeypatch, tmp_path, name):
    """The trap of the round's challenge (a): R55's refusals return the constant HEART_TF_STRENGTH (1.25), so a branch
    that reused them would build every R74 file as OTK125. Here the step is called once, at the table's strength and
    the category form, and the file is not OTK125's bytes (pinned below)."""
    seen = _spy(monkeypatch)
    out, rec, _ = _build(bt, tmp_path, name, *R74[name])
    assert seen == [(float(TABLE[name][1]), "category")] and seen[0][0] != 1.25
    assert sha(out) != R55_BEFORE["category"][0]
    assert rec["trend_strength"]["strength"] == float(TABLE[name][1])
    seen.clear()
    b, _, _ = _build(bt, tmp_path, "otk125", *TF["category"])
    assert seen == [(1.25, "category")] and sha(b) == R55_BEFORE["category"][0]


def test_the_strength_is_the_tables_even_when_r55s_constant_moves(bt, tf_env, monkeypatch, tmp_path):
    """R55's constant moved to another value: row R74's files keep their bytes and step strength (nothing of R74 reads
    it), and R55's own build follows it."""
    first = {n: sha(_build(bt, tmp_path, f"{n}_a", *R74[n])[0]) for n in NAMES}
    monkeypatch.setattr(bt, "HEART_TF_STRENGTH", 9.0)
    seen = _spy(monkeypatch)
    again = {n: sha(_build(bt, tmp_path, f"{n}_b", *R74[n])[0]) for n in NAMES}
    assert again == first and seen == [(float(TABLE[n][1]), "category") for n in NAMES]


def test_the_table_strength_scales_otk125s_category_factors_and_damping_only(tf, cw):
    """At each table strength the step's arguments are OTK125's (the category form at 1.25) with the factors and the
    fallback damping scaled by s / 1.25, the draw and the typing unchanged: one change."""
    from vec import harness

    ctx, _ = cw
    base = tf.form_arguments(ctx, 1.25, "category")
    for _, raw, _ in TABLE.values():
        s = float(raw)
        got = tf.form_arguments(ctx, s, "category")
        np.testing.assert_array_equal(got.rows, base.rows)
        np.testing.assert_array_equal(got.like.obs["b0_type"], base.like.obs["b0_type"])
        assert got.damping == harness.B0_DAMPING * s and sorted(got.mapping) == sorted(base.mapping)
        for k, v in base.mapping.items():
            np.testing.assert_allclose(got.mapping[k], v * (s / 1.25), rtol=1e-12, atol=0)
        assert got.record["strength"] == s and got.record["split_types"] == base.record["split_types"]


# ------------------------------------------------------------------------------------------------ the files
@pytest.mark.parametrize("name", NAMES)
def test_each_file_keeps_otk125s_rows_coordinates_count_size_and_half(bt, tf_env, tmp_path, name):
    """Rows, coordinates, cell number, size field, D's own record fields and, after halve.py at seed 1000001, the
    half's rows and coordinates are OTK125's; the expression is not."""
    sys.path.insert(0, str(REPO / "experiments" / "builds"))
    import halve

    o, o_rec, _ = _build(bt, tmp_path, "otk125", *TF["category"])
    t, t_rec, _ = _build(bt, tmp_path, name, *R74[name])
    A, B = ad.read_h5ad(o), ad.read_h5ad(t)
    assert A.n_obs == B.n_obs and list(A.obs_names) == list(B.obs_names) and list(A.var_names) == list(B.var_names)
    np.testing.assert_array_equal(_coords(o), _coords(t))
    for k in ("n_cells", "seed", "count_decision", "trim", "size", "z_mode", "board", "dry_run", "depth_step",
              "scope_option"):
        assert t_rec[k] == o_rec[k], k
    ho, ht = tmp_path / "otk125_half.h5ad", tmp_path / f"{name}_half.h5ad"
    halve.halve(str(o), str(ho), 1000001)
    halve.halve(str(t), str(ht), 1000001)
    HA, HB = ad.read_h5ad(ho), ad.read_h5ad(ht)
    assert list(HA.obs_names) == list(HB.obs_names) and HA.n_obs == HB.n_obs
    np.testing.assert_array_equal(_coords(ho), _coords(ht))
    assert not np.array_equal(_X(o), _X(t)) and not np.array_equal(_X(ho), _X(ht))


def test_the_four_files_differ_from_each_other_and_from_otk125(bt, tf_env, tmp_path):
    X = {n: _X(_build(bt, tmp_path, n, *R74[n])[0]) for n in NAMES}
    X["OTK125"] = _X(_build(bt, tmp_path, "otk125", *TF["category"])[0])
    keys = list(X)
    for i, f in enumerate(keys):
        for g in keys[i + 1:]:
            assert not np.array_equal(X[f], X[g]), (f, g)


@pytest.mark.parametrize("name", NAMES)
def test_the_same_bytes_twice(bt, tf_env, tmp_path, name):
    a, _, ta = _build(bt, tmp_path, "a", *R74[name])
    b, _, tb = _build(bt, tmp_path, "b", *R74[name])
    assert a.read_bytes() == b.read_bytes() and ta == tb


@pytest.mark.parametrize("name", NAMES)
def test_the_record_names_the_file_its_strength_string_row_order_and_command_and_no_cell_number(bt, tf_env, tmp_path,
                                                                                                capsys, name):
    option, raw, stem = TABLE[name]
    _, rec, text = _build(bt, tmp_path, name, *R74[name])
    ts = rec["trend_strength"]
    assert ts["option"] == f"`{option}`" and ts["form"] == "category"
    assert ts["strength"] == float(raw) and ts["strength_option"] == raw  # the row's string verbatim
    assert ts["file"] == name and ts["file_stem"] == stem and ts["base_file"] == "OTK125"
    assert ts["captain_ruling"] == "R74" and ts["order_row"] == R74_FIRST  # the row's first cell as the ledger holds it
    assert ts["order"] == bt.HEART_R74_ORDER.format(option=option, file=name, strength=raw)
    assert ts["order"].startswith("experiments/harness/TRIES.md row R74, the captain's order of 2026-11-25")
    assert f"`{option}` ({name})" in ts["order"] and f"strength {raw} in place of 1.25" in ts["order"]
    assert ts["command"] == "heart " + " ".join(R74[name]) == bt.HEART_R74_COMMAND.format(option=option)
    assert "1.25" not in ts["command"] and "R55" not in ts["command"]
    assert ts["step"] == bt.HEART_R74_STEP.format(strength=raw) and f"strength {raw} (row R74's table)" in ts["step"]
    assert ts["budget_days"] == pytest.approx(0.5) and ts["test_budget_days"] == pytest.approx(float(raw) * 0.5)
    assert ts["trend"]["stages"] == ["e8_25_late", "e8_5", "e8_75", "start", "ref"]
    assert set(ts) >= set(bt.HEART_TF_RECORD["category"]) and {t["status"] for t in ts["types"]} <= {
        "shifted", "skipped: fewer than 20 cells at a stage"}
    assert not ({"tag", "trend_whole"} & set(ts))  # R55's form tag and the whole form's field are not this record's
    assert rec["entry"] == f"trend strength {float(raw):g} form category organ_scope (organ)"
    assert "expression" not in rec
    assert rec["scope_option"]["command"] == bt.HEART_ORGAN_SCOPE and rec["n_cells"] == "max_cells"
    assert "R52" not in text and "R49" not in text and text.count('"OTK125"') == 1  # base_file only

    def walk(x, path="trend_strength"):
        if isinstance(x, dict):
            for k, v in x.items():
                assert not k.startswith("n_") and not any(w in k.lower() for w in ("count", "min_")), k
                walk(v, f"{path}.{k}")
        elif isinstance(x, list):
            for i, v in enumerate(x):
                walk(v, f"{path}[{i}]")
        else:  # strings, floats and booleans only: an integer could be a cell number
            assert isinstance(x, (str, float, bool)) or x is None, (path, type(x))

    walk(ts)
    assert not re.search(rf"(?<![\d.]){SECRET_MAX}(?![\d.])", text + capsys.readouterr().out)
    assert json.loads(text) == rec


@pytest.mark.parametrize("name", NAMES)
def test_the_reads_are_otk125s(bt, tf_env, tmp_path, name):
    _build(bt, tmp_path, "otk125", *TF["category"])
    want = list(tf_env)
    assert want == ["boards", "E10.5.h5ad", "E9.5.h5ad", "E8.75.h5ad", "E8.5.h5ad", "E8.25_late.h5ad", "boards"]
    tf_env.clear()
    _build(bt, tmp_path, name, *R74[name])
    assert tf_env == want


def test_the_test_boards_own_name_is_its_default(bt, tf_env, tmp_path):
    out, rec, _ = _build(bt, tmp_path, "named", *R74["OTK075"], "--board", "T2:heart:test_extrap")
    default, _, _ = _build(bt, tmp_path, "default", *R74["OTK075"])
    assert out.read_bytes() == default.read_bytes() and rec["dry_run"] is False


# ------------------------------------------------------------------------------------------------ refusals
def _with(name: str, flag: str, value: str) -> list[str]:
    """Row R74's command for `name` with the value after `flag` replaced."""
    argv = list(R74[name])
    argv[argv.index(flag) + 1] = value
    return argv


NEW_REFUSED = [
    [*ORG, *TABLE["OTK10"][0].split(), "--z-mode", "raw"],  # no row
    *(_with("OTK10", "--captain-ruling", r) for r in ("R55", "R7", "R740", "r74", "R74 ", " R74", "R74a", "", "R73",
                                                     "R75", "R29", "R64")),
    *(_with("OTK15", "--trend-strength", v) for v in ("1", "1.00", "1.0 ", " 1.0", "+1.0", "1.0e0", "01.0", ".75",
                                                     "0.750", "0.7", "1.7", "1.8", "1.25", "0.5", "2.0", "1.2", "1.50",
                                                     "1.75 ", "1,5", "one")),
    *(_with("OTK075", "--trend-form", f) for f in ("contrast", "whole", "portion", "Category", "category ", "",
                                                   "banked", "category_portion", "categories")),
    [*ORG, "--trend-form", "category", "--z-mode", "raw", "--captain-ruling", "R74"],  # no strength
    *(_with("OTK175", "--entry", e) for e in ("s1", "s2", "s2_trim", "core_trim")),
    [*ORG, *TABLE["OTK175"][0].split(), "--captain-ruling", "R74"],  # the recorded z mode (the default)
    _with("OTK175", "--z-mode", "recorded"),
    *([*R74["OTK10"], *f.split()] for f in ("--expression copy", "--trend", "--depth", "--herd", "--cfm")),
    *([*R74[n], "--seed", s] for n in ("OTK10", "OTK15") for s in ("1000002", "1000000")),
    [*R74["OTK075"], "--board", "T2:heart:val_extrap"], [*R74["OTK075"], "--board", "T2:heart:val_interp"],
]


@pytest.mark.parametrize("argv", NEW_REFUSED)
def test_every_other_use_is_refused_before_anything_is_read(bt, tf_env, tmp_path, argv):
    msg = _why(bt, tmp_path, *argv)
    assert "refused" in msg or "needs --captain-ruling R55" in msg or "orders four files" in msg
    assert tf_env == []


def test_the_refusals_say_why(bt, tf_env, ledger, tmp_path):
    listed = "; ".join(f"{n} `{v[0]}`" for n, v in TABLE.items())
    four = ("--captain-ruling R74 orders four files, each with exactly its own strings on heart --entry organ_scope "
            f"<OPTION> --z-mode raw --captain-ruling R74: {listed}")
    assert _why(bt, tmp_path, *TF["category"][:-1], "R74") == four  # R55's strength under row R74
    assert _why(bt, tmp_path, *_with("OTK10", "--trend-strength", "1")) == four
    assert _why(bt, tmp_path, *_with("OTK10", "--trend-form", "whole")) == four
    cmd = "heart " + " ".join(R74["OTK10"])
    for extra in (["--expression", "copy"], ["--trend"], ["--cfm"]):
        assert _why(bt, tmp_path, *R74["OTK10"], *extra) == (
            f"row R74's file OTK10 is refused except as {cmd}: --entry organ_scope and --z-mode raw only, and with "
            "--expression copy, --trend, --depth, --herd and --cfm")
    assert _why(bt, tmp_path, *_with("OTK10", "--entry", "core_trim")).startswith(
        f"row R74's file OTK10 is refused except as {cmd}")
    for extra in (["--board", "T2:heart:val_extrap"], ["--seed", "1000002"]):
        assert _why(bt, tmp_path, *R74["OTK10"], *extra) == (
            f"row R74's file OTK10 is refused on another board or seed: it builds {cmd} at the default seed on the "
            "test board only (no dry run)")
    # R55's route as before: an R74 string under row R55, or under no row, is R55's refusal
    r55 = "heart " + " ".join(TF["category"])
    old = (f"--trend-form category is refused except as {r55}: OTS125's line only (--entry organ_scope, "
           "--trend-strength 1.25, --z-mode raw), and with --expression copy, --trend, --depth, --herd and --cfm")
    assert _why(bt, tmp_path, *_with("OTK10", "--captain-ruling", "R55")) == old
    assert _why(bt, tmp_path, *R74["OTK10"][:-2]) == old
    ledger.write_text(LEDGER_TF)  # row R74 gone
    opt = TABLE["OTK15"][0]
    assert _why(bt, tmp_path, *R74["OTK15"]) == (
        "row R74's file OTK15 is refused: experiments/harness/TRIES.md must hold exactly one row R74 whose first cell "
        "is the captain's order (\"ruling, no run\") naming `--entry organ_scope`, `--z-mode raw`, "
        f"`{opt}` and `--captain-ruling R74`, each as a whole backticked string, `{opt}` followed there by "
        "\" (OTK15)\"")
    assert tf_env == []


def test_every_earlier_value_and_row_meets_its_old_refusal_with_r74_in_the_ledger(bt, tf_env, tmp_path):
    """Without --trend-form, --captain-ruling R74 is another row to every earlier line, refused as before."""
    assert _why(bt, tmp_path, *OTSF[1.25][:-1], "R74") == (
        "--entry organ_scope --trend-strength 1.25 is refused: experiments/harness/TRIES.md must hold exactly one row "
        "R74 whose first cell is the captain's order (\"ruling, no run\") naming `--entry organ_scope`, "
        "`--trend-strength 0.75`, `--trend-strength 1.25` and `--captain-ruling R74`")
    assert _why(bt, tmp_path, *OTS15[:-1], "R74") == (
        "--entry organ_scope --trend-strength 1.5 is refused with --captain-ruling R74: only row R52 orders it")
    assert _why(bt, tmp_path, *OTS[1.0][:-1], "R74") == (
        "--entry organ_scope --trend-strength 1.0 is refused with --captain-ruling R74: only row R29 orders it")
    assert _why(bt, tmp_path, *ORG, "--trend-strength", "1.75", "--z-mode", "raw", "--captain-ruling", "R74") == (
        "--trend-strength 1.75 is refused: it takes 0.5 or 1.0 only")
    assert _why(bt, tmp_path, *ORGAN, "--captain-ruling", "R74").startswith(
        "heart --captain-ruling goes with --trend-strength only")
    assert tf_env == []


# ------------------------------------------------------------------------------------------------ the row check
def test_the_r74_check_on_ledgers_of_the_required_form(bt, monkeypatch, tmp_path):
    _use(bt, monkeypatch, tmp_path / "ok.md", LEDGER_R74)
    assert all(bt._heart_r74_ordered(n) for n in NAMES) and all(bt._heart_tf_ordered(f) for f in FORMS)
    _use(bt, monkeypatch, tmp_path / "only.md", R74_ROW)
    assert all(bt._heart_r74_ordered(n) for n in NAMES) and not bt._heart_tf_ordered("category")
    # a row naming two of the files orders those two only
    two = R74_FIRST.replace("`--trend-strength 1.0 --trend-form category` (OTK10), ", "").replace(
        "`--trend-strength 0.75 --trend-form category` (OTK075), ", "")
    _use(bt, monkeypatch, tmp_path / "two.md", f"| {two} | heart | x |\n")
    assert [n for n in NAMES if bt._heart_r74_ordered(n)] == ["OTK15", "OTK175"]
    _use(bt, monkeypatch, tmp_path / "ok2.md", LEDGER_R74)
    assert not any(bt._heart_r74_ordered(n) for n in ("OTK125", "OTK1", "", "otk10", "OTK10 "))


@pytest.mark.parametrize("drop", ["`--entry organ_scope`", "`--z-mode raw`", "`--captain-ruling R74`, `--entry",
                                  "ruling, no run", "the captain's order"])
def test_the_r74_check_fails_without_each_required_string(bt, monkeypatch, tmp_path, drop):
    if drop == "`--captain-ruling R74`, `--entry":  # both mentions of the ruling's string
        first = R74_FIRST.replace("`--captain-ruling R74`", "--captain-ruling R74")
    else:
        first = R74_FIRST.replace(drop, "something else")
    _use(bt, monkeypatch, tmp_path / "t.md", LEDGER_TF + f"| {first} | heart | x |\n")
    assert not any(bt._heart_r74_ordered(n) for n in NAMES)


@pytest.mark.parametrize("text", [
    LEDGER_TF, LEDGER_R74 + R74_ROW,  # none, two
    LEDGER_R74.replace("| R74 (", "| R740 ("), LEDGER_R74.replace("| R74 (", "| R75 ("),
    LEDGER_R74.replace("| R74 (", "| R74b ("),
    LEDGER_R74.replace("`--captain-ruling R74`", "`--captain-ruling R73`"),
    # the file's string inside a longer backticked string, without backticks, or paired with another file's name
    LEDGER_R74.replace("`--trend-strength 1.0 --trend-form category`",
                       "`--trend-strength 1.0 --trend-form category --z-mode raw`"),
    LEDGER_R74.replace("`--trend-strength 1.0 --trend-form category`",
                       "--trend-strength 1.0 --trend-form category"),
    LEDGER_R74.replace("`--trend-strength 1.0 --trend-form category` (OTK10)",
                       "`--trend-strength 1.0 --trend-form category` (OTK15)"),
    LEDGER_R74.replace("`--trend-strength 1.0 --trend-form category` (OTK10)",
                       "`--trend-strength 1.0 --trend-form category`, OTK10"),
    LEDGER_R74.replace("`--trend-strength 1.0 --trend-form category`",
                       "`--trend-strength 1.0` `--trend-form category`"),
    LEDGER_R74.replace("`--entry organ_scope`, `--z-mode raw`", "`--entry organ_scope --z-mode raw`"),
    LEDGER_R74.replace("ruling, no run", "awaiting a ruling"),
    LEDGER_R74.replace("R74 (the captain's order", "R74 (the captain's named exception"),
    # the strings outside the first cell do not count
    LEDGER_TF + R74_ROW.replace("`--trend-strength 1.0 --trend-form category` (OTK10), ", "").replace(
        "| heart |", "| `--trend-strength 1.0 --trend-form category` (OTK10) |")])
def test_the_otk10_build_is_refused_on_a_ledger_without_its_order(bt, tf_env, ledger, tmp_path, text):
    ledger.write_text(text)
    msg = _why(bt, tmp_path, *R74["OTK10"])
    assert "must hold exactly one row R74" in msg
    assert tf_env == []


def test_the_committed_ledger_holds_r74_in_the_builders_form():
    """experiments/harness/TRIES.md as committed, read by code only: row R74's first cell is the captain's order naming
    the four files' strings in the form the builder takes (and row R55's, the identity step's)."""
    sys.path.insert(0, str(REPO / "experiments" / "builds"))
    import build_test as bt

    assert bt.TRIES_FILE == REPO / "experiments" / "harness" / "TRIES.md"
    assert all(bt._heart_r74_ordered(n) for n in bt.HEART_R74_FILES) and list(bt.HEART_R74_FILES) == NAMES
    assert bt._heart_tf_ordered("category")


# ------------------------------------------------------------------------------------------------ earlier rows
@pytest.mark.parametrize("form", FORMS)
def test_row_r55s_forms_keep_their_bytes_and_records(bt, tf_env, tmp_path, form):
    out, rec, text = _build(bt, tmp_path, form, *TF[form])
    assert sha(out) == R55_BEFORE[form][0]
    assert hashlib.sha256(text.encode()).hexdigest() == R55_BEFORE[form][1] or canon(text) == R55_BEFORE[form][2]
    assert "R74" not in text and rec["trend_strength"]["captain_ruling"] == "R55"


def _broken(bt, monkeypatch) -> None:
    def broken(*a, **k):
        raise AssertionError("row R74's code was reached")
    for f in ("_r74_refusals", "_r74_record", "_heart_r74_ordered", "_r74_strength"):
        monkeypatch.setattr(bt, f, broken)


def test_no_earlier_form_reaches_row_r74s_code(bt, tf_env, monkeypatch, tmp_path):
    _broken(bt, monkeypatch)
    for f in FORMS:
        out, _, _ = _build(bt, tmp_path, f, *TF[f])
        assert sha(out) == R55_BEFORE[f][0]


def test_no_earlier_strength_command_reaches_row_r74s_code(bt, organ_env, monkeypatch, tmp_path):
    _broken(bt, monkeypatch)
    for name, argv in (("d", ORGAN), ("o5", OTS[0.5]), ("o10", OTS[1.0]), ("o075", OTSF[0.75]), ("o125", OTSF[1.25]),
                       ("o15", OTS15)):
        _, rec, text = _build(bt, tmp_path, name, *argv)
        assert "R74" not in text
