"""`experiments/builds/build_test.py heart --entry organ_scope --trend-strength <S> --trend-form category --z-mode raw
--captain-ruling R77`: row R77's two files (experiments/harness/TRIES.md row R77, the captain's order of 2026-11-26,
ruling, no run), OTK125's command (row R55's category file) with one change, the category form's strength S, taken from
the builder's one literal table HEART_R77_TABLE in the row's order: OTK05 (0.5) and OTK20 (2.0). Row R74's route, its
code copied beside row R74's and R74's own code left as it was.

On the synthetic stages and synthetic test board (its limit never shown) of tests/test_heart_trend_line_h.py, with a
temporary ledger holding rows R55's, R74's and R77's orders: the table and its names; the three traps of the round's
challenge (HEART_REF4_challenge.md section 6): the strength is the table's, never R55's constant 1.25 (so a file is not
OTK125's bytes, nor an R74 file's, and the two differ); the record is row R77's own, chosen before R55's form record
(else it names OTK125 with R55's command); build_heart's order row reads row R77's first cell (0.5 and 2.0 are not
R49's values). Its rows, coordinates, cell number, size field and half rows at seed 1000001 are OTK125's; the same
bytes twice; the record holds the strength string verbatim and no cell number; every other use is refused before
anything is read (near-miss strengths "2", "2.00", ".5" and "0.50" among them, R74's four strings under R77 and R77's
two under R74); "0.5" and "2.0" without --trend-form keep their plain routes (R29's file, and the old refusal); the row
check on temporary ledgers and, read by code, on the committed ledger. Old against new: the builder of main at
4d91239 (before this change, loaded from git) and this one give the same bytes and records for R55's four forms and
R74's four files (also pinned here) and the same route for every command but row R77's; no earlier command reaches row
R77's code. Nothing is read from data/ and nothing is scored."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import re
import subprocess
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from conftest import git_env  # noqa: E402
from test_build_heart_iso_trim import _build, canon  # noqa: E402
from test_build_heart_organ_scope import ORGAN, _use_stage, organ_stage  # noqa: E402
from test_heart_r74_build import LEDGER_R74, R55_BEFORE, R74, R74_FIRST  # noqa: E402
from test_heart_r74_build import NAMES as R74_NAMES  # noqa: E402
from test_heart_r74_build import TABLE as R74_TABLE  # noqa: E402
from test_heart_trend_build import SECRET_MAX, STAGES, install, sha  # noqa: E402
from test_heart_trend_line_d import OTS, _organ_size  # noqa: E402
from test_heart_trend_line_f import _X, OTS_BEFORE, _coords  # noqa: E402
from test_heart_trend_line_g import OTS15, OTSF, OTSF_BEFORE  # noqa: E402
from test_heart_trend_line_h import (  # noqa: E402
    FORMS,
    ORG,
    OTS15_BEFORE,
    TF,
    _reads,
    _why,
    relabel,
    tf_context,
    tf_organ_stage,
)

# row R77's table as the row names it, in its order: file -> (option string, strength string, file stem)
TABLE = {"OTK05": ("--trend-strength 0.5 --trend-form category", "0.5", "T2heart_organ_ts05_tk"),
         "OTK20": ("--trend-strength 2.0 --trend-form category", "2.0", "T2heart_organ_ts20_tk")}
NAMES = list(TABLE)
R77 = {n: [*ORG, *TABLE[n][0].split(), "--z-mode", "raw", "--captain-ruling", "R77"] for n in NAMES}
R77_FIRST = ("R77 (the captain's order of 2026-11-26, ruling, no run: \"Yes to OTK05 and OTK20\", a synthetic "
             "stand-in; the heart builder is to take on `build_test.py heart` exactly these strings, each only with "
             "`--captain-ruling R77`: `--entry organ_scope`, `--z-mode raw`, `--trend-strength 0.5 --trend-form "
             "category` (OTK05), `--trend-strength 2.0 --trend-form category` (OTK20), `--captain-ruling R77`; each "
             "file then halved by `halve.py ... --seed 1000001` and the half goes up; not a try)")
R77_ROW = f"| {R77_FIRST} | heart | refinement round 4's two files |\n"
LEDGER_R77 = LEDGER_R74 + R77_ROW  # rows R55, R74 (tests/test_heart_r74_build.py's stand-ins) and R77
# sha256 of (the .h5ad, its .build.json text, the record's `canon`) that build_test.py heart wrote for row R74's four
# files on tests/test_heart_trend_line_h.py's trend-form environment with tests/test_heart_r74_build.py's ledger,
# computed with the builder of main at 4d91239 (before row R77's code) by the same builds as below; row R55's four forms
# gave R55_BEFORE (tests/test_heart_r74_build.py) there, rechecked here.
R74_BEFORE = {
    "OTK10": ("be5d6e7f51f60e1d8186b768fb16bf249819228fe83942f42c11b9a8753b2a2e",
              "cb8a6c3ee103cee7547a28d25b192c7d42b16d521711bcc58842ec46f8a1e8a5",
              "4d1daf69239677951e213f382b80d6440565ad1334b9b36600dda7caf44e9c56"),
    "OTK15": ("c6bbdf29c55d1f2791be8f9642565158b00d918dab68d97290d63e29357aed9e",
              "206e3fd2b4524464c14cdd31fe3f7f75c12476a3367ce66df8e12b393bbd5d85",
              "000c810ff32dc4fcf0cd02353f2aae1404c43e4b259871f35ca8f0a83f549bec"),
    "OTK075": ("e74d65eb333a3b5cd60b8c131f169702113750a8b940887a93edf29fcf28e9fd",
               "f9a513c4ddca5cc3fb374b9d029d6ea47dc7cb015c927efeda179bf962e0dd8c",
               "ae1ab5793e1399a54515366c080890906cea232c494f96d4206cb41a88ec2e60"),
    "OTK175": ("fc43b40340854998743f4947095911b1cf405e5d5cb271e222aceae3840562ce",
               "f2567721bb89af49353f0a7d253189c78b75537dbcc5d609be77c600d42681b7",
               "8e632dbd130104bcbcb1c2859f7ca60841135ea3b3c20da628341c1c3f8cf5b0"),
}
BEFORE_COMMIT = "4d91239fc6537a57004a424fa0fb566732cef12d"  # main before row R77's code


@pytest.fixture(scope="module")
def bt():
    sys.path.insert(0, str(REPO / "experiments" / "builds"))
    import build_test

    return build_test


@pytest.fixture(scope="module")
def tf():
    from vec import trend_forms

    return trend_forms


@pytest.fixture(scope="module")
def old_bt(bt, tmp_path_factory):
    """The builder of main at BEFORE_COMMIT (before row R77's code), loaded from git beside this one."""
    out = subprocess.run(["git", "-C", str(REPO), "show", f"{BEFORE_COMMIT}:experiments/builds/build_test.py"],
                         capture_output=True, text=True, env=git_env())
    if out.returncode != 0:
        pytest.skip(f"git cannot show the builder at {BEFORE_COMMIT} here")
    assert "HEART_R77" not in out.stdout and "HEART_R74_TABLE" in out.stdout
    path = tmp_path_factory.mktemp("before_r77") / "build_test_before_r77.py"
    path.write_text(out.stdout)
    spec = importlib.util.spec_from_file_location("build_test_before_r77", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture
def ledger(bt, monkeypatch, tmp_path):
    """The temporary ledger: rows R55's and R74's orders (tests/test_heart_r74_build.py's) and row R77's."""
    path = tmp_path / "TRIES.md"
    path.write_text(LEDGER_R77)
    monkeypatch.setattr(bt, "TRIES_FILE", path)
    return path


@pytest.fixture
def organ_env(bt, monkeypatch, own_repo, ledger):
    """tests/test_heart_trend_line_h.py's organ environment (the earlier lines'), with rows R55, R74 and R77."""
    reads: list[str] = []
    install(bt, monkeypatch.setattr, reads)
    _use_stage(bt, monkeypatch, organ_stage(), reads)
    _organ_size(bt, monkeypatch)
    more = _reads(bt, monkeypatch)
    return reads, more


def _tf_install(m, monkeypatch, reads: list) -> None:
    """tests/test_heart_trend_line_h.py's trend-form environment on builder module `m`: install()'s board and stages,
    every stage relabelled, E10.5 as `tf_organ_stage`, the organ size stand-in; `reads` records every stage and board
    read, in order."""
    install(m, monkeypatch.setattr, reads)
    load = m.load_stage

    def wrapped(p, **kw):
        name = Path(p).name
        if name == "E10.5.h5ad":
            reads.append(name)
            return tf_organ_stage()
        return relabel(load(p, **kw), STAGES[name][2])

    monkeypatch.setattr(m, "load_stage", wrapped)
    _organ_size(m, monkeypatch)
    boards = m.load_boards

    def load_boards():
        reads.append("boards")
        return boards()

    monkeypatch.setattr(m, "load_boards", load_boards)


@pytest.fixture
def tf_env(bt, monkeypatch, own_repo, ledger):
    reads: list[str] = []
    _tf_install(bt, monkeypatch, reads)
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


def _broken(bt, monkeypatch, names: tuple[str, ...]) -> None:
    def broken(*a, **k):
        raise AssertionError("code that must not be reached was reached")
    for f in names:
        monkeypatch.setattr(bt, f, broken)


R77_CODE = ("_r77_refusals", "_r77_record", "_heart_r77_ordered", "_r77_strength")


# ------------------------------------------------------------------------------------------------ the table
def test_the_table_is_the_rows_strings_in_its_order(bt):
    assert bt.HEART_R77_ROW == "R77" and bt.HEART_R77_FORM == "category" and bt.HEART_R77_BASE_FILE == "OTK125"
    assert [tuple(p) for p in bt.HEART_R77_TABLE] == [(n, v[0]) for n, v in TABLE.items()]
    assert bt.HEART_R77_FILES == {n: (v[0], v[1]) for n, v in TABLE.items()}
    assert bt.HEART_R77_STEMS == {n: v[2] for n, v in TABLE.items()}
    assert bt.HEART_R77_BY_OPTION == {v[0]: n for n, v in TABLE.items()}
    assert bt.HEART_R77_NAMED == ("--entry organ_scope", "--z-mode raw", "--captain-ruling R77")
    for n in NAMES:
        assert bt.HEART_R77_COMMAND.format(option=TABLE[n][0]) == "heart " + " ".join(R77[n])
        assert f"`{TABLE[n][0]}` ({n})" in R77_FIRST  # each string, whole and backticked, then its file's name
    # no R77 value is R55's 1.25 or one of row R74's four, and no option string is row R74's
    strengths = {float(v[1]) for v in TABLE.values()}
    assert len(strengths) == 2 and not strengths & {1.25, *(float(v[1]) for v in R74_TABLE.values())}
    assert not set(bt.HEART_R77_BY_OPTION) & set(bt.HEART_R74_BY_OPTION)
    # row R74's and R55's constants as they were
    assert [tuple(p) for p in bt.HEART_R74_TABLE] == [(n, v[0]) for n, v in R74_TABLE.items()]
    assert bt.HEART_R74_ROW == "R74" and bt.HEART_R74_NAMED == ("--entry organ_scope", "--z-mode raw",
                                                                 "--captain-ruling R74")
    assert bt.HEART_TF_STRENGTH == 1.25 and bt.HEART_TF_STRENGTH_CHOICE == "1.25" and bt.HEART_TF_ROW == "R55"
    assert {f: bt.HEART_TF_COMMANDS[f] for f in FORMS} == {f: "heart " + " ".join(TF[f]) for f in FORMS}


def test_the_table_is_the_committed_rows_strings():
    """R77's first cell as committed in experiments/harness/TRIES.md, read by code and never shown: the pairs it
    names (each backticked string followed by its file's name) are the table's, in its order."""
    sys.path.insert(0, str(REPO / "experiments" / "builds"))
    import build_test as bt

    lines = [ln for ln in bt.TRIES_FILE.read_text(encoding="utf-8").splitlines() if ln.startswith("| R77 ")]
    assert len(lines) == 1
    first = lines[0][2:].split("|", 1)[0]
    pairs = re.findall(r"`(--trend-strength [^`]*)` \((OTK[0-9]+)\)", first)
    assert [(n, o) for o, n in pairs] == [(n, v[0]) for n, v in TABLE.items()]
    for s in ("--entry organ_scope", "--z-mode raw", "--captain-ruling R77"):
        assert f"`{s}`" in first


def test_the_option_is_off_by_default(bt):
    a = bt.build_parser().parse_args(["heart", *ORGAN, "--out", "runs/x/h.h5ad"])
    assert a.trend_form is None and a.captain_ruling is None and a.trend_strength is None
    assert bt._r77_file(a) is None
    for argv in (TF["category"], R74["OTK10"]):  # R55's file and an R74 file: not row R77's
        a = bt.build_parser().parse_args(["heart", *argv, "--out", "runs/x/h.h5ad"])
        assert bt._r77_file(a) is None


@pytest.mark.parametrize("name, option", [
    ("OTK05", "--trend-strength 2.0 --trend-form category"),  # the other file's string
    ("OTK20", "--trend-strength 0.5 --trend-form category"),
    ("OTK5", "--trend-strength 0.5 --trend-form category"), ("OTK050", "--trend-strength 0.5 --trend-form category"),
    ("OTK2", "--trend-strength 2.0 --trend-form category"), ("OTK200", "--trend-strength 2.0 --trend-form category"),
    ("OTK20", "--trend-strength 2 --trend-form category"), ("OTK20", "--trend-strength 2.00 --trend-form category"),
    ("OTK05", "--trend-strength .5 --trend-form category"), ("OTK05", "--trend-strength 0.50 --trend-form category"),
    ("OTK5", "--trend-strength .5 --trend-form category"), ("OTK2", "--trend-strength 2 --trend-form category"),
    ("OTK125", "--trend-strength 1.25 --trend-form category"),  # R55's own string
    ("OTK10", "--trend-strength 1.0 --trend-form category"),  # row R74's strings and values
    ("OTK15", "--trend-strength 1.5 --trend-form category"), ("OTK150", "--trend-strength 1.50 --trend-form category"),
    ("OTK075", "--trend-strength 0.75 --trend-form category"),
    ("OTK175", "--trend-strength 1.75 --trend-form category"),
    ("OTK05", "--trend-strength 0.5 --trend-form contrast"), ("OTK05", "--trend-strength 0.5 --trend-form Category"),
    ("OTK05", "--trend-form category --trend-strength 0.5"), ("OTK05", "--trend-strength 0.5  --trend-form category"),
    ("OTK05", " --trend-strength 0.5 --trend-form category"), ("OTK05", "--trend-strength 0.5 --trend-form category "),
    ("OTK20", "--trend-strength 2.0 --trend-form category --z-mode raw"),
    ("otk05", "--trend-strength 0.5 --trend-form category"), ("OTX05", "--trend-strength 0.5 --trend-form category"),
    ("OTK20.h5ad", "--trend-strength 2.0 --trend-form category"),
])
def test_a_name_must_match_its_option_string(bt, name, option):
    with pytest.raises(ValueError, match="does not match its option string"):
        bt._r77_strength(name, option)


def test_each_row_pair_gives_its_strength_string(bt):
    assert {n: bt._r77_strength(n, v[0]) for n, v in TABLE.items()} == {n: v[1] for n, v in TABLE.items()}


# ------------------------------------------------------------------------------------------------ trap (a): strength
@pytest.mark.parametrize("name", NAMES)
def test_each_file_steps_at_the_tables_strength_never_r55s(bt, tf_env, monkeypatch, tmp_path, name):
    """R55's refusals return the constant HEART_TF_STRENGTH (1.25), so a branch that reused them would build each R77
    file as OTK125. Here the step is called once, at the table's strength and the category form, and the file is not
    OTK125's bytes nor any R74 file's (all pinned)."""
    seen = _spy(monkeypatch)
    out, rec, _ = _build(bt, tmp_path, name, *R77[name])
    assert seen == [(float(TABLE[name][1]), "category")] and seen[0][0] != 1.25
    assert sha(out) not in {R55_BEFORE["category"][0], *(v[0] for v in R74_BEFORE.values())}
    assert rec["trend_strength"]["strength"] == float(TABLE[name][1])
    seen.clear()
    b, _, _ = _build(bt, tmp_path, "otk125", *TF["category"])
    assert seen == [(1.25, "category")] and sha(b) == R55_BEFORE["category"][0]


def test_the_strength_is_the_tables_even_when_r55s_constant_or_r74s_route_moves(bt, tf_env, monkeypatch, tmp_path):
    """R55's constant moved to another value and R55's and R74's refusals broken: row R77's files keep their bytes and
    step strength (nothing of R77 reads them)."""
    first = {n: sha(_build(bt, tmp_path, f"{n}_a", *R77[n])[0]) for n in NAMES}
    monkeypatch.setattr(bt, "HEART_TF_STRENGTH", 9.0)
    _broken(bt, monkeypatch, ("_trend_form_refusals", "_r74_refusals", "_heart_tf_ordered", "_heart_r74_ordered"))
    seen = _spy(monkeypatch)
    again = {n: sha(_build(bt, tmp_path, f"{n}_b", *R77[n])[0]) for n in NAMES}
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


# ------------------------------------------------------------------------------------------------ trap (a): record
@pytest.mark.parametrize("name", NAMES)
def test_the_record_is_row_r77s_own_never_r55s_form_record(bt, tf_env, monkeypatch, tmp_path, name):
    """The trap of the round's challenge (a): `_trend_strength_record` sends any other --trend-form build to R55's
    `_trend_form_record`, which names OTK125 with R55's command; row R77's branch comes first. With R55's and R74's
    record functions broken the files still build, with row R77's record; with R77's broken they do not."""
    _broken(bt, monkeypatch, ("_trend_form_record", "_r74_record"))
    _, rec, text = _build(bt, tmp_path, name, *R77[name])
    ts = rec["trend_strength"]
    assert ts["file"] == name and ts["captain_ruling"] == "R77" and ts["command"] == "heart " + " ".join(R77[name])
    assert "R55" not in ts["command"] and "1.25" not in ts["command"] and "tag" not in ts
    _broken(bt, monkeypatch, ("_r77_record",))
    with pytest.raises(AssertionError, match="must not be reached"):
        _build(bt, tmp_path, f"{name}_b", *R77[name])


def test_the_record_function_takes_row_r77s_branch_before_the_forms(bt, cw):
    """`_trend_strength_record` itself: with form "category" and r77 a file, row R77's record; with r77 None, R55's
    record as before (OTK125, R55's command); with r74 a file, R74's."""
    from vec import trend_forms

    ctx, ref = cw
    _, _, _, rec = trend_forms.form_step(ctx, 0.5, "category", ref)
    order = R77_FIRST
    got = bt._trend_strength_record(0.5, rec, "organ_scope", "R77", order, "category", None, "OTK05")
    assert got["file"] == "OTK05" and got["captain_ruling"] == "R77" and got["order_row"] == order
    assert got["command"] == bt.HEART_R77_COMMAND.format(option=TABLE["OTK05"][0])
    old = bt._trend_strength_record(0.5, rec, "organ_scope", "R77", order, "category", None)
    assert old["file"] == "OTK125" and old["command"] == bt.HEART_TF_COMMANDS["category"]  # the trap's record
    r74 = bt._trend_strength_record(1.0, rec, "organ_scope", "R74", R74_FIRST, "category", "OTK10")
    assert r74["file"] == "OTK10" and r74["captain_ruling"] == "R74"


# ------------------------------------------------------------------------------------------------ trap (b): order row
@pytest.mark.parametrize("name", NAMES)
def test_the_order_row_is_row_r77s_first_cell(bt, tf_env, monkeypatch, tmp_path, name):
    """build_heart's order-row line reads the ruling row's first cell only for R49's values, R52's and row R74's files;
    0.5 and 2.0 are none of these, so without row R77 in that line the record would hold no order row."""
    assert float(TABLE[name][1]) not in bt.HEART_FURTHER_FILES and float(TABLE[name][1]) != bt.HEART_OTS15_STRENGTH
    asked: list = []
    own = bt._ruling_first_cell

    def first_cell(row):
        asked.append(row)
        return own(row)
    monkeypatch.setattr(bt, "_ruling_first_cell", first_cell)
    _, rec, _ = _build(bt, tmp_path, name, *R77[name])
    assert asked == ["R77"] and rec["trend_strength"]["order_row"] == R77_FIRST


# ------------------------------------------------------------------------------------------------ the files
@pytest.mark.parametrize("name", NAMES)
def test_each_file_keeps_otk125s_rows_coordinates_count_size_and_half(bt, tf_env, tmp_path, name):
    """Rows, coordinates, cell number, size field, D's own record fields and, after halve.py at seed 1000001, the
    half's rows and coordinates are OTK125's; the expression is not."""
    sys.path.insert(0, str(REPO / "experiments" / "builds"))
    import halve

    o, o_rec, _ = _build(bt, tmp_path, "otk125", *TF["category"])
    t, t_rec, _ = _build(bt, tmp_path, name, *R77[name])
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


def test_the_two_files_differ_from_each_other_from_otk125_and_from_row_r74s(bt, tf_env, tmp_path):
    X = {n: _X(_build(bt, tmp_path, n, *R77[n])[0]) for n in NAMES}
    X["OTK125"] = _X(_build(bt, tmp_path, "otk125", *TF["category"])[0])
    for n in R74_NAMES:
        X[n] = _X(_build(bt, tmp_path, n, *R74[n])[0])
    keys = list(X)
    for i, f in enumerate(keys):
        for g in keys[i + 1:]:
            assert not np.array_equal(X[f], X[g]), (f, g)


@pytest.mark.parametrize("name", NAMES)
def test_the_same_bytes_twice(bt, tf_env, tmp_path, name):
    a, _, ta = _build(bt, tmp_path, "a", *R77[name])
    b, _, tb = _build(bt, tmp_path, "b", *R77[name])
    assert a.read_bytes() == b.read_bytes() and ta == tb


@pytest.mark.parametrize("name", NAMES)
def test_the_record_names_the_file_its_strength_string_row_order_and_command_and_no_cell_number(bt, tf_env, tmp_path,
                                                                                                capsys, name):
    option, raw, stem = TABLE[name]
    _, rec, text = _build(bt, tmp_path, name, *R77[name])
    ts = rec["trend_strength"]
    assert ts["option"] == f"`{option}`" and ts["form"] == "category"
    assert ts["strength"] == float(raw) and ts["strength_option"] == raw  # the row's string verbatim
    assert ts["file"] == name and ts["file_stem"] == stem and ts["base_file"] == "OTK125"
    assert ts["captain_ruling"] == "R77" and ts["order_row"] == R77_FIRST  # the row's first cell as the ledger holds it
    assert ts["order"] == bt.HEART_R77_ORDER.format(option=option, file=name, strength=raw)
    assert ts["order"].startswith("experiments/harness/TRIES.md row R77, the captain's order of 2026-11-26")
    assert f"`{option}` ({name})" in ts["order"] and f"strength {raw} in place of 1.25" in ts["order"]
    assert ts["command"] == "heart " + " ".join(R77[name]) == bt.HEART_R77_COMMAND.format(option=option)
    assert "1.25" not in ts["command"] and "R55" not in ts["command"] and "R74" not in ts["command"]
    assert ts["step"] == bt.HEART_R77_STEP.format(strength=raw) and f"strength {raw} (row R77's table)" in ts["step"]
    assert ts["budget_days"] == pytest.approx(0.5) and ts["test_budget_days"] == pytest.approx(float(raw) * 0.5)
    assert ts["trend"]["stages"] == ["e8_25_late", "e8_5", "e8_75", "start", "ref"]
    assert set(ts) >= set(bt.HEART_TF_RECORD["category"]) and {t["status"] for t in ts["types"]} <= {
        "shifted", "skipped: fewer than 20 cells at a stage"}
    assert not ({"tag", "trend_whole"} & set(ts))  # R55's form tag and the whole form's field are not this record's
    assert rec["entry"] == f"trend strength {float(raw):g} form category organ_scope (organ)"
    assert "expression" not in rec
    assert rec["scope_option"]["command"] == bt.HEART_ORGAN_SCOPE and rec["n_cells"] == "max_cells"
    assert "R52" not in text and "R49" not in text and "R74" not in text and text.count('"OTK125"') == 1  # base_file

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
    _build(bt, tmp_path, name, *R77[name])
    assert tf_env == want


def test_the_test_boards_own_name_is_its_default(bt, tf_env, tmp_path):
    out, rec, _ = _build(bt, tmp_path, "named", *R77["OTK20"], "--board", "T2:heart:test_extrap")
    default, _, _ = _build(bt, tmp_path, "default", *R77["OTK20"])
    assert out.read_bytes() == default.read_bytes() and rec["dry_run"] is False


# ------------------------------------------------------------------------------------------------ refusals
def _with(name: str, flag: str, value: str) -> list[str]:
    """Row R77's command for `name` with the value after `flag` replaced."""
    argv = list(R77[name])
    argv[argv.index(flag) + 1] = value
    return argv


NEAR_STRENGTHS = ("2", "2.00", ".5", "0.50", "2.0 ", " 2.0", "+2.0", "2.0e0", "02.0", "0.05", "0.5 ", " 0.5", "5e-1",
                  "00.5", "1.0", "1.5", "0.75", "1.75", "1.25", "2.5", "0.25", "1.2", "2,0", "0,5", "two", "")
NEW_REFUSED = [
    [*ORG, *TABLE["OTK05"][0].split(), "--z-mode", "raw"],  # no row
    *(_with("OTK05", "--captain-ruling", r) for r in ("R55", "R74", "R7", "R770", "r77", "R77 ", " R77", "R77a", "",
                                                     "R76", "R78", "R29", "R65", "R077")),
    *(_with("OTK20", "--trend-strength", v) for v in NEAR_STRENGTHS),
    *(_with("OTK05", "--trend-strength", v) for v in NEAR_STRENGTHS),
    *(_with("OTK05", "--trend-form", f) for f in ("contrast", "whole", "portion", "Category", "category ", "",
                                                  "banked", "category_portion", "categories", "insitu")),
    [*ORG, "--trend-form", "category", "--z-mode", "raw", "--captain-ruling", "R77"],  # no strength
    *(_with("OTK20", "--entry", e) for e in ("s1", "s2", "s2_trim", "core_trim")),
    [*ORG, *TABLE["OTK20"][0].split(), "--captain-ruling", "R77"],  # the recorded z mode (the default)
    _with("OTK20", "--z-mode", "recorded"),
    *([*R77["OTK05"], *f.split()] for f in ("--expression copy", "--trend", "--depth", "--herd", "--cfm",
                                            "--run-profile " + "0" * 64, "--insitu-chords " + "0" * 64)),
    *([*R77[n], "--seed", s] for n in NAMES for s in ("1000002", "1000000")),
    [*R77["OTK05"], "--board", "T2:heart:val_extrap"], [*R77["OTK20"], "--board", "T2:heart:val_interp"],
    # row R74's four strings under row R77, and row R77's two under rows R74 and R55
    *([*ORG, *R74_TABLE[n][0].split(), "--z-mode", "raw", "--captain-ruling", "R77"] for n in R74_NAMES),
    *(_with(n, "--captain-ruling", r) for n in NAMES for r in ("R74", "R55")),
]


@pytest.mark.parametrize("argv", NEW_REFUSED)
def test_every_other_use_is_refused_before_anything_is_read(bt, tf_env, tmp_path, argv):
    msg = _why(bt, tmp_path, *argv)
    assert ("refused" in msg or "needs --captain-ruling R55" in msg or "orders two files" in msg
            or "orders four files" in msg or "needs the sha256" in msg or "is the captain's order" in msg), msg
    assert tf_env == []


def test_the_refusals_say_why(bt, tf_env, ledger, tmp_path):
    listed = "; ".join(f"{n} `{v[0]}`" for n, v in TABLE.items())
    two = ("--captain-ruling R77 orders two files, each with exactly its own strings on heart --entry organ_scope "
           f"<OPTION> --z-mode raw --captain-ruling R77: {listed}")
    assert _why(bt, tmp_path, *TF["category"][:-1], "R77") == two  # R55's strength under row R77
    for v in ("2", "2.00", ".5", "0.50"):
        assert _why(bt, tmp_path, *_with("OTK20", "--trend-strength", v)) == two
    for n in R74_NAMES:  # row R74's strings under row R77
        argv = [*ORG, *R74_TABLE[n][0].split(), "--z-mode", "raw", "--captain-ruling", "R77"]
        assert _why(bt, tmp_path, *argv) == two
    assert _why(bt, tmp_path, *_with("OTK05", "--trend-form", "whole")) == two
    four = "; ".join(f"{n} `{v[0]}`" for n, v in R74_TABLE.items())
    for n in NAMES:  # row R77's strings under row R74
        assert _why(bt, tmp_path, *_with(n, "--captain-ruling", "R74")) == (
            "--captain-ruling R74 orders four files, each with exactly its own strings on heart --entry organ_scope "
            f"<OPTION> --z-mode raw --captain-ruling R74: {four}")
    cmd = "heart " + " ".join(R77["OTK05"])
    for extra in (["--expression", "copy"], ["--trend"], ["--cfm"]):
        assert _why(bt, tmp_path, *R77["OTK05"], *extra) == (
            f"row R77's file OTK05 is refused except as {cmd}: --entry organ_scope and --z-mode raw only, and with "
            "--expression copy, --trend, --depth, --herd and --cfm")
    assert _why(bt, tmp_path, *_with("OTK05", "--entry", "core_trim")).startswith(
        f"row R77's file OTK05 is refused except as {cmd}")
    for extra in (["--board", "T2:heart:val_extrap"], ["--seed", "1000002"]):
        assert _why(bt, tmp_path, *R77["OTK05"], *extra) == (
            f"row R77's file OTK05 is refused on another board or seed: it builds {cmd} at the default seed on the "
            "test board only (no dry run)")
    # R55's route as before: an R77 string under row R55, or under no row, is R55's refusal
    r55 = "heart " + " ".join(TF["category"])
    old = (f"--trend-form category is refused except as {r55}: OTS125's line only (--entry organ_scope, "
           "--trend-strength 1.25, --z-mode raw), and with --expression copy, --trend, --depth, --herd and --cfm")
    assert _why(bt, tmp_path, *_with("OTK20", "--captain-ruling", "R55")) == old
    assert _why(bt, tmp_path, *R77["OTK20"][:-2]) == old
    ledger.write_text(LEDGER_R74)  # row R77 gone
    opt = TABLE["OTK20"][0]
    assert _why(bt, tmp_path, *R77["OTK20"]) == (
        "row R77's file OTK20 is refused: experiments/harness/TRIES.md must hold exactly one row R77 whose first cell "
        "is the captain's order (\"ruling, no run\") naming `--entry organ_scope`, `--z-mode raw`, "
        f"`{opt}` and `--captain-ruling R77`, each as a whole backticked string, `{opt}` followed there by "
        "\" (OTK20)\"")
    assert tf_env == []


def test_the_plain_routes_of_05_and_20_are_as_before(bt, organ_env, monkeypatch, tmp_path):
    """Without --trend-form, 0.5 stays on R29's plain route (OTS05, its pinned bytes) and --captain-ruling R77 is
    another row to it; 2.0 is refused as before (the builder never took R65's plain 2.0)."""
    _broken(bt, monkeypatch, R77_CODE)
    out, rec, text = _build(bt, tmp_path, "ots05", *OTS[0.5])
    assert sha(out) == OTS_BEFORE[0.5][0] and "R77" not in text
    assert hashlib.sha256(text.encode()).hexdigest() == OTS_BEFORE[0.5][1] or canon(text) == OTS_BEFORE[0.5][2]
    assert _why(bt, tmp_path, *OTS[0.5][:-1], "R77") == (
        "--entry organ_scope --trend-strength 0.5 is refused with --captain-ruling R77: only row R29 orders it")
    for row in ("R77", "R65", "R29"):
        assert _why(bt, tmp_path, *ORG, "--trend-strength", "2.0", "--z-mode", "raw", "--captain-ruling", row) == (
            "--trend-strength 2.0 is refused: it takes 0.5 or 1.0 only")
    assert _why(bt, tmp_path, *ORGAN, "--captain-ruling", "R77").startswith(
        "heart --captain-ruling goes with --trend-strength only")


# ------------------------------------------------------------------------------------------------ the row check
def test_the_r77_check_on_ledgers_of_the_required_form(bt, monkeypatch, tmp_path):
    _use(bt, monkeypatch, tmp_path / "ok.md", LEDGER_R77)
    assert all(bt._heart_r77_ordered(n) for n in NAMES) and all(bt._heart_r74_ordered(n) for n in R74_NAMES)
    assert all(bt._heart_tf_ordered(f) for f in FORMS)
    _use(bt, monkeypatch, tmp_path / "only.md", R77_ROW)
    assert all(bt._heart_r77_ordered(n) for n in NAMES) and not any(bt._heart_r74_ordered(n) for n in R74_NAMES)
    # a row naming one of the files orders that one only
    one = R77_FIRST.replace("`--trend-strength 0.5 --trend-form category` (OTK05), ", "")
    _use(bt, monkeypatch, tmp_path / "one.md", f"| {one} | heart | x |\n")
    assert [n for n in NAMES if bt._heart_r77_ordered(n)] == ["OTK20"]
    _use(bt, monkeypatch, tmp_path / "ok2.md", LEDGER_R77)
    assert not any(bt._heart_r77_ordered(n) for n in ("OTK125", "OTK10", "OTK5", "", "otk05", "OTK20 ", *R74_NAMES))


@pytest.mark.parametrize("drop", ["`--entry organ_scope`", "`--z-mode raw`", "`--captain-ruling R77`, `--entry",
                                  "ruling, no run", "the captain's order"])
def test_the_r77_check_fails_without_each_required_string(bt, monkeypatch, tmp_path, drop):
    if drop == "`--captain-ruling R77`, `--entry":  # both mentions of the ruling's string
        first = R77_FIRST.replace("`--captain-ruling R77`", "--captain-ruling R77")
    else:
        first = R77_FIRST.replace(drop, "something else")
    _use(bt, monkeypatch, tmp_path / "t.md", LEDGER_R74 + f"| {first} | heart | x |\n")
    assert not any(bt._heart_r77_ordered(n) for n in NAMES)


@pytest.mark.parametrize("text", [
    LEDGER_R74, LEDGER_R77 + R77_ROW,  # none, two
    LEDGER_R77.replace("| R77 (", "| R770 ("), LEDGER_R77.replace("| R77 (", "| R78 ("),
    LEDGER_R77.replace("| R77 (", "| R77b ("),
    LEDGER_R77.replace("`--captain-ruling R77`", "`--captain-ruling R74`"),
    # the file's string inside a longer backticked string, without backticks, or paired with another file's name
    LEDGER_R77.replace("`--trend-strength 0.5 --trend-form category`",
                       "`--trend-strength 0.5 --trend-form category --z-mode raw`"),
    LEDGER_R77.replace("`--trend-strength 0.5 --trend-form category`", "--trend-strength 0.5 --trend-form category"),
    LEDGER_R77.replace("`--trend-strength 0.5 --trend-form category` (OTK05)",
                       "`--trend-strength 0.5 --trend-form category` (OTK20)"),
    LEDGER_R77.replace("`--trend-strength 0.5 --trend-form category` (OTK05)",
                       "`--trend-strength 0.5 --trend-form category`, OTK05"),
    LEDGER_R77.replace("`--trend-strength 0.5 --trend-form category`",
                       "`--trend-strength 0.5` `--trend-form category`"),
    LEDGER_R77.replace("`--trend-strength 0.5 --trend-form category`", "`--trend-strength 0.50 --trend-form category`"),
    LEDGER_R77.replace("`--entry organ_scope`, `--z-mode raw`, `--trend-strength 0.5",
                       "`--entry organ_scope --z-mode raw`, `--trend-strength 0.5"),
    LEDGER_R77.replace("ruling, no run: \"Yes to OTK05", "awaiting a ruling: \"Yes to OTK05"),
    LEDGER_R77.replace("R77 (the captain's order", "R77 (the captain's named exception"),
    # the strings outside the first cell do not count
    LEDGER_R74 + R77_ROW.replace("`--trend-strength 0.5 --trend-form category` (OTK05), ", "").replace(
        "| heart |", "| `--trend-strength 0.5 --trend-form category` (OTK05) |")])
def test_the_otk05_build_is_refused_on_a_ledger_without_its_order(bt, tf_env, ledger, tmp_path, text):
    ledger.write_text(text)
    msg = _why(bt, tmp_path, *R77["OTK05"])
    assert "must hold exactly one row R77" in msg
    assert tf_env == []


def test_the_committed_ledger_holds_r77_in_the_builders_form():
    """experiments/harness/TRIES.md as committed, read by code only: row R77's first cell is the captain's order naming
    the two files' strings in the form the builder takes (and rows R55's and R74's, as before)."""
    sys.path.insert(0, str(REPO / "experiments" / "builds"))
    import build_test as bt

    assert bt.TRIES_FILE == REPO / "experiments" / "harness" / "TRIES.md"
    assert all(bt._heart_r77_ordered(n) for n in bt.HEART_R77_FILES) and list(bt.HEART_R77_FILES) == NAMES
    assert all(bt._heart_r74_ordered(n) for n in bt.HEART_R74_FILES) and bt._heart_tf_ordered("category")


# ------------------------------------------------------------------------------------------------ old against new
@pytest.mark.parametrize("form", FORMS)
def test_row_r55s_forms_keep_their_bytes_and_records(bt, tf_env, tmp_path, form):
    out, rec, text = _build(bt, tmp_path, form, *TF[form])
    assert sha(out) == R55_BEFORE[form][0]
    assert hashlib.sha256(text.encode()).hexdigest() == R55_BEFORE[form][1] or canon(text) == R55_BEFORE[form][2]
    assert "R77" not in text and rec["trend_strength"]["captain_ruling"] == "R55"


@pytest.mark.parametrize("name", R74_NAMES)
def test_row_r74s_files_keep_their_bytes_and_records(bt, tf_env, tmp_path, name):
    out, rec, text = _build(bt, tmp_path, name, *R74[name])
    assert sha(out) == R74_BEFORE[name][0]
    assert hashlib.sha256(text.encode()).hexdigest() == R74_BEFORE[name][1] or canon(text) == R74_BEFORE[name][2]
    assert "R77" not in text and rec["trend_strength"]["captain_ruling"] == "R74"
    assert rec["trend_strength"]["order_row"] == R74_FIRST


def test_the_old_and_new_builders_give_r55s_and_r74s_files_alike(bt, old_bt, monkeypatch, own_repo, ledger, tmp_path):
    """The builder before row R77's code, loaded beside this one, on the same environment and ledger: the same bytes and
    the same record text for R55's four forms and R74's four files."""
    monkeypatch.setattr(old_bt, "TRIES_FILE", ledger)
    reads_new: list[str] = []
    reads_old: list[str] = []
    _tf_install(bt, monkeypatch, reads_new)
    _tf_install(old_bt, monkeypatch, reads_old)
    for key, argv in (*((f, TF[f]) for f in FORMS), *((n, R74[n]) for n in R74_NAMES)):
        new, _, new_text = _build(bt, tmp_path / "new", key, *argv)
        old, _, old_text = _build(old_bt, tmp_path / "old", key, *argv)
        assert new.read_bytes() == old.read_bytes() and new_text == old_text, key
        want = R55_BEFORE.get(key) or R74_BEFORE[key]
        assert sha(new) == want[0], key
    assert reads_new == reads_old


def _route(m, argv: list[str]) -> tuple:
    """What builder module `m` makes of a heart command before it reads anything: the in-situ or run-profile route,
    the argument parser's or the strength line's refusal, or the strength with row R74's file, row R77's file (None
    for the builder before R77) and the ruling row the order-row line reads."""
    try:
        a = m.build_parser().parse_args(["heart", *argv, "--out", "x.h5ad"])
    except SystemExit as e:
        return ("parser", str(e))
    if getattr(a, "insitu_chords", None) is not None or getattr(a, "trend_form", None) == m.HEART_OIC_FORM:
        return ("oic",)
    if getattr(a, "run_profile", None) is not None:
        return ("xpc",)
    try:
        s = m._trend_strength_refusals(a)
    except SystemExit as e:
        return ("refused", str(e))
    r74 = m._r74_file(a)
    r77 = m._r77_file(a) if hasattr(m, "_r77_file") else None
    row = (a.captain_ruling if r74 is not None or r77 is not None or s in m.HEART_FURTHER_FILES
           or s == m.HEART_OTS15_STRENGTH else None)
    return ("built", s, r74, r77, row)


ROUTE_ROWS = ("R77", "R74", "R55", "R29", "R28", "R49", "R52", "R65", "R70", "R56", "R7", "R770", "r77", "")
ROUTE_STRENGTHS = ("0.5", "2.0", "1.0", "1.5", "0.75", "1.75", "1.25", "2", "2.00", ".5", "0.50", "0.25", "3.0")
ROUTES = [
    *([*ORG, "--trend-strength", v, "--trend-form", f, "--z-mode", z, "--captain-ruling", r]
      for v in ROUTE_STRENGTHS for f in ("category", "contrast", "whole", "portion", "insitu", "Category")
      for z in ("raw", "recorded") for r in ROUTE_ROWS),
    *([*ORG, "--trend-strength", v, "--z-mode", "raw", "--captain-ruling", r] for v in ROUTE_STRENGTHS
      for r in ROUTE_ROWS),
    *([*ORG, "--trend-strength", v, "--z-mode", "raw"] for v in ROUTE_STRENGTHS),
    *([*ORG, "--trend-strength", v, "--trend-form", "category", "--z-mode", "raw"] for v in ROUTE_STRENGTHS),
    *([*e, "--trend-strength", v, "--trend-form", "category", "--z-mode", "raw", "--captain-ruling", r]
      for e in (["--entry", "core_trim"], ["--entry", "s1"]) for v in ("0.5", "2.0", "1.0") for r in ("R77", "R74")),
    *([*R77[n], *x] for n in NAMES for x in (["--seed", "1000002"], ["--board", "T2:heart:val_extrap"],
                                             ["--board", "T2:heart:test_extrap"], ["--trend"], ["--cfm"])),
    *([*R74[n], *x] for n in R74_NAMES for x in ([], ["--board", "T2:heart:test_extrap"], ["--seed", "1000002"])),
    ORGAN, [*ORGAN, "--captain-ruling", "R77"], ["--entry", "s1"], ["--entry", "s2_trim", "--z-mode", "raw"],
    [*ORG, "--trend-form", "category", "--z-mode", "raw", "--captain-ruling", "R77"],
]


def test_every_command_but_row_r77s_takes_the_old_route(bt, old_bt, monkeypatch, tmp_path):
    """Over every command above, read before anything is loaded: a command outside row R77's --trend-form route gives
    the same strength, R74 file, ruling row and refusal text in the old and new builders; under row R77 with
    --trend-form, the old builder gave R55's refusal, and the new one builds exactly R77's two strings at the
    table's strength (on the test board, named or by default) and refuses everything else."""
    path = tmp_path / "TRIES.md"
    path.write_text(LEDGER_R77)
    monkeypatch.setattr(bt, "TRIES_FILE", path)
    monkeypatch.setattr(old_bt, "TRIES_FILE", path)
    built = []
    for argv in ROUTES:
        new, old = _route(bt, argv), _route(old_bt, argv)
        ruling = argv[argv.index("--captain-ruling") + 1] if "--captain-ruling" in argv else None
        if ruling == "R77" and "--trend-form" in argv and new[0] not in ("oic", "xpc", "parser"):
            assert old[0] == "refused", argv  # R55's route refused it
            if new[0] == "built":
                built.append((argv, new))
            continue
        assert new == old, argv
    assert sorted({(r[1], r[3], r[4]) for _, r in built}) == [(0.5, "OTK05", "R77"), (2.0, "OTK20", "R77")]
    for argv, _ in built:
        option = " ".join(argv[argv.index("--trend-strength"):argv.index("--trend-strength") + 4])
        assert option in bt.HEART_R77_BY_OPTION and "--z-mode raw" in " ".join(argv)


def test_no_earlier_form_or_r74_file_reaches_row_r77s_code(bt, tf_env, monkeypatch, tmp_path):
    _broken(bt, monkeypatch, R77_CODE)
    for f in FORMS:
        out, _, _ = _build(bt, tmp_path, f, *TF[f])
        assert sha(out) == R55_BEFORE[f][0]
    for n in R74_NAMES:
        out, _, _ = _build(bt, tmp_path, n, *R74[n])
        assert sha(out) == R74_BEFORE[n][0]


def test_no_earlier_strength_command_reaches_row_r77s_code(bt, organ_env, monkeypatch, tmp_path):
    _broken(bt, monkeypatch, R77_CODE)
    pins = {"o5": OTS_BEFORE[0.5], "o10": OTS_BEFORE[1.0], "o075": OTSF_BEFORE[0.75], "o125": OTSF_BEFORE[1.25],
            "o15": OTS15_BEFORE}
    for name, argv in (("d", ORGAN), ("o5", OTS[0.5]), ("o10", OTS[1.0]), ("o075", OTSF[0.75]), ("o125", OTSF[1.25]),
                       ("o15", OTS15)):
        out, rec, text = _build(bt, tmp_path, name, *argv)
        assert "R77" not in text
        if name in pins:
            assert sha(out) == pins[name][0], name
