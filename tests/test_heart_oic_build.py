"""`experiments/builds/build_test.py heart --entry organ_scope --trend-strength 1.25 --trend-form insitu --insitu-chords
<sha256> --z-mode raw --captain-ruling R70`: file OIC125 (row R70 of experiments/harness/TRIES.md, the captain's order
of 2026-11-25, ruling, no run), OTK125's command with the in-situ chords of row H38's G0 added (vec.trend_forms'
insitu form; HEART_MOON1_design.md P1 with HEART_MOON1_challenge.md's FIX 1 and FIX 6).

On the synthetic stages and test board of tests/test_heart_trend_line_h.py (its `tf_env`), with a temporary ledger (row
R55 as that file has it, a synthetic row R70 and a synthetic row H38 whose outcome cell holds G0's PASS statement and
the read's PASS) and a temporary chord file: the in-situ step with its chord out is the category form bit for bit; with
a chord it adds P exactly where P exists (FIX 1) and nowhere else (FIX 6, to 1e-12); the file keeps OTK125's rows,
coordinates, cell number, size field and half rows and changes the expression only; the same bytes twice; the record
names R70, OIC125, the command and the sha256; every refusal (another ruling or option, no R70 order, no PASS in row
H38, a missing or other chord file, a chord a type cannot carry), before any file is written. Nothing is read from
data/ and nothing is scored."""

# ruff: noqa: F811 (the fixtures of tests/test_heart_trend_line_h.py are imported and taken by name)
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
sys.path.insert(0, str(REPO / "experiments" / "harness"))
import h38_row  # noqa: E402
from test_build_heart_iso_trim import _build  # noqa: E402
from test_heart_trend_line_f import _coords  # noqa: E402
from test_heart_trend_line_h import LEDGER_TF, TF, bt, ledger, tf, tf_context, tf_env  # noqa: E402,F401

R70_FIRST = ("R70 (the captain's order of 2026-11-25, ruling, no run: \"Both\", a synthetic stand-in for the heart's "
             "in-situ late-contrast file OIC125; built only after its gates read PASS, only with `--captain-ruling "
             "R70`; not a try)")
OIC = ["--entry", "organ_scope", "--trend-strength", "1.25", "--trend-form", "insitu", "--insitu-chords", "{sha}",
       "--z-mode", "raw", "--captain-ruling", "R70"]


def oic(sha: str, ruling: str = "R70") -> list[str]:
    return [a.format(sha=sha) for a in OIC[:-1]] + [ruling]


def chords_doc(genes, cm_idx=range(0, 30), endo_idx=range(10, 20), epi=None) -> dict:
    P = {"cm": [None] * len(genes), "endocardium": [None] * len(genes)}
    for i in cm_idx:
        P["cm"][i] = 0.3
    for i in endo_idx:
        P["endocardium"][i] = -0.2
    if epi is not None:
        P["epicardium"] = [0.25 if i < 15 else None for i in range(len(genes))]
    return {"row": "H38", "file": "OIC125", "ref": "E10.5", "genes": list(genes), "P": P}


@pytest.fixture
def oic_env(bt, tf_env, monkeypatch, tmp_path):
    """tf_env, a chord file of OIC125 on the synthetic genes, and a ledger with rows R55, R70 and H38 (PASS)."""
    from test_heart_trend_build import heart_test_board

    genes = list(heart_test_board().genes)
    path = tmp_path / "results" / "chords_E10.5.json"
    path.parent.mkdir()

    def write(doc, cell=None):
        path.write_text(json.dumps(doc))
        sha = hashlib.sha256(path.read_bytes()).hexdigest()
        cell = cell if cell is not None else (h38_row.g0_statement(sha, "b" * 64, "c" * 64) + " "
                                              + h38_row.outcome_statement("PASS"))
        tries = tmp_path / "TRIES_oic.md"
        tries.write_text(LEDGER_TF + f"| {R70_FIRST} | heart | file OIC125 |\n"
                         + f"| {h38_row.FIRST_CELL} | heart | " + " | ".join(["x"] * 8) + f" | {cell} |\n")
        monkeypatch.setattr(bt, "TRIES_FILE", tries)
        return sha

    monkeypatch.setattr(bt, "HEART_OIC_CHORDS", path)
    sha = write(chords_doc(genes))
    return {"reads": tf_env, "genes": genes, "sha": sha, "write": write, "path": path}


def _why(bt, tmp_path, *argv) -> str:
    out = tmp_path / "refused" / "h.h5ad"
    with pytest.raises(SystemExit) as e:
        bt.main(["heart", *argv, "--out", str(out)])
    assert not out.exists()
    return str(e.value)


# ---- the form ----------------------------------------------------------------------------------------------------
def test_zero_chords_are_the_category_form_bit_for_bit(bt, tf):
    ctx, _ = tf_context(bt)
    st0 = tf.form_step(ctx, 1.25, "category")[0]
    genes = list(ctx.inputs["ref"].var_names)
    for chords in (None, {}, {"cm": np.full(len(genes), np.nan)}):
        z = tf.insitu_step(ctx, 1.25, chords)[0]
        assert np.array_equal(z.X, st0.X) and np.array_equal(z.coords, st0.coords)
    some = {"cm": np.full(len(genes), 0.3)}
    z = tf.insitu_step(ctx, 1.25, some, zero=True)[0]
    assert np.array_equal(z.X, st0.X)


def test_the_chord_is_added_where_it_exists_and_nowhere_else(bt, tf):
    ctx, _ = tf_context(bt)
    genes = list(ctx.inputs["ref"].var_names)
    P = {"cm": np.full(len(genes), np.nan), "endocardium": np.full(len(genes), np.nan)}
    P["cm"][:30] = 0.3
    P["endocardium"][10:20] = -0.2
    base = tf.form_arguments(ctx, 1.25, "category")
    args = tf.insitu_arguments(ctx, 1.25, P)
    assert tf.insitu_checks(ctx, args, base, P) == []
    new, old = tf.applied(ctx, args), tf.applied(ctx, base)
    i, j = new["types"].index("CM"), old["types"].index("CM")
    np.testing.assert_allclose(new["logF"][i][:30] - old["logF"][j][:30], 0.3, atol=1e-12)
    np.testing.assert_allclose(new["logF"][i][30:], old["logF"][j][30:], atol=1e-12)
    st, counts, rows, rec = tf.insitu_step(ctx, 1.25, P)
    st0, c0, rows0, _ = tf.form_step(ctx, 1.25, "category")
    assert np.array_equal(rows, rows0) and np.array_equal(st.coords, st0.coords) and not np.array_equal(st.X, st0.X)
    assert rec["form"] == "insitu" and rec["carried"] == ["cm", "endocardium"]
    lib = np.asarray(counts.sum(1)).ravel().mean() / np.asarray(c0.sum(1)).ravel().mean()
    assert abs(lib - 1.0) < 0.05  # depth kept, in expectation
    anti = tf.insitu_step(ctx, 1.25, {k: -v for k, v in P.items()})[0]
    assert not np.array_equal(anti.X, st.X)


def test_a_misplaced_or_idle_chord_is_a_fault(bt, tf):
    ctx, _ = tf_context(bt)
    genes = list(ctx.inputs["ref"].var_names)
    P = {"cm": np.full(len(genes), np.nan)}
    P["cm"][:5] = 0.2
    base = tf.form_arguments(ctx, 1.25, "category")
    args = tf.insitu_arguments(ctx, 1.25, P)
    args.mapping["Endothelial/endocardium"] = args.mapping["Endothelial/endocardium"] + 1e-9
    assert any("off P's entries" in w for w in tf.insitu_checks(ctx, args, base, P))
    # synthetic epicardium cells carry a merged lineage of their own (not Untracked): the key finds none of them, so a
    # chord for epicardium would be idle, and the step refuses
    with pytest.raises(ValueError):
        tf.insitu_step(ctx, 1.25, {"epicardium": np.full(len(genes), 0.1)})
    with pytest.raises(ValueError):
        tf.insitu_arguments(ctx, 1.25, {"blood": np.zeros(len(genes))})
    with pytest.raises(ValueError):
        tf.insitu_arguments(ctx, 1.25, {"cm": np.zeros(3)})


def test_epicardium_under_its_own_key_carries_its_chord(bt, tf):
    """Released-style epicardium (merged null, so Untracked in OTK125): the synthetic Epi cells (lineage Tiny) made
    Untracked at every stage; the key types them, 0 + P_epi is applied where P exists, OTK125's 0 elsewhere."""
    from vec import harness

    ctx, _ = tf_context(bt)
    inputs = {}
    for r, a in ctx.inputs.items():
        m = np.asarray(a.obs["merged"]).astype(object)
        m[m == "Tiny"] = "Untracked"
        inputs[r] = harness._shallow(a, merged=m.astype(str))
    ctx = harness.BuildContext(task=ctx.task, board=ctx.board, inputs=inputs, seed=ctx.seed)
    genes = list(inputs["ref"].var_names)
    P = {"epicardium": np.full(len(genes), np.nan)}
    P["epicardium"][:12] = 0.25
    base = tf.form_arguments(ctx, 1.25, "category")
    args = tf.insitu_arguments(ctx, 1.25, P)
    assert tf.insitu_checks(ctx, args, base, P) == []
    new = tf.applied(ctx, args)
    assert {"type": tf.EPI_KEY, "status": "shifted"} in new["statuses"]
    row = new["logF"][new["types"].index(tf.EPI_KEY)]
    np.testing.assert_allclose(row[:12], 0.25, atol=1e-12)
    assert np.all(row[12:] == 0)


# ---- the builder -------------------------------------------------------------------------------------------------
def test_the_option_is_off_by_default(bt):
    a = bt.build_parser().parse_args(["heart", "--entry", "organ_scope", "--out", "runs/x/h.h5ad"])
    assert a.insitu_chords is None
    assert bt.HEART_OIC_ROW == "R70" and bt.HEART_OIC_FILE == "OIC125"
    assert bt.HEART_OIC_COMMAND.format(sha256="S") == "heart " + " ".join(oic("S"))


def test_oic125_keeps_otk125s_rows_coordinates_count_size_and_half(bt, oic_env, tmp_path):
    sys.path.insert(0, str(REPO / "experiments" / "builds"))
    import halve

    k, k_rec, _ = _build(bt, tmp_path, "otk125", *TF["category"])
    o, o_rec, text = _build(bt, tmp_path, "oic125", *oic(oic_env["sha"]))
    A, B = ad.read_h5ad(k), ad.read_h5ad(o)
    assert A.n_obs == B.n_obs and list(A.obs_names) == list(B.obs_names) and list(A.var_names) == list(B.var_names)
    np.testing.assert_array_equal(_coords(k), _coords(o))
    assert not np.array_equal(np.asarray(A.X.todense() if hasattr(A.X, "todense") else A.X),
                              np.asarray(B.X.todense() if hasattr(B.X, "todense") else B.X))
    for key in ("n_cells", "seed", "count_decision", "trim", "size", "z_mode", "board", "dry_run", "depth_step",
                "scope_option"):
        assert o_rec[key] == k_rec[key], key
    ts = o_rec["trend_strength"]
    assert ts["file"] == "OIC125" and ts["captain_ruling"] == "R70" and ts["base_file"] == "OTK125"
    assert ts["command"] == "heart " + " ".join(oic(oic_env["sha"])) and ts["carried"] == ["cm", "endocardium"]
    assert o_rec["insitu_chords"]["sha256"] == oic_env["sha"] and "R70" in ts["order_row"]
    hk, ho = tmp_path / "k_half.h5ad", tmp_path / "o_half.h5ad"
    halve.halve(str(k), str(hk))
    halve.halve(str(o), str(ho))
    assert list(ad.read_h5ad(hk).obs_names) == list(ad.read_h5ad(ho).obs_names)
    o2, _, _ = _build(bt, tmp_path, "oic125b", *oic(oic_env["sha"]))
    assert o.read_bytes() == o2.read_bytes()


def test_refusals_before_any_file(bt, oic_env, tmp_path):
    sha = oic_env["sha"]
    assert "64 lowercase" in _why(bt, tmp_path, *oic("xyz"))
    assert "R70 alone" in _why(bt, tmp_path, *oic(sha)[:-2])
    assert "only row R70" in _why(bt, tmp_path, *oic(sha, "R55"))
    assert "refused except as" in _why(bt, tmp_path, *oic(sha), "--cfm")
    assert "refused except as" in _why(bt, tmp_path, *oic(sha), "--seed", "1000003")
    assert "refused except as" in _why(bt, tmp_path, *[x if x != "1.25" else "1.0" for x in oic(sha)])
    no_sha = [x for x in oic(sha) if x not in ("--insitu-chords", sha)]
    assert "64 lowercase" in _why(bt, tmp_path, *no_sha)
    # the gates have not read PASS
    oic_env["write"](chords_doc(oic_env["genes"]), cell=h38_row.g0_statement(sha, "b" * 64, "c" * 64))
    assert "not read PASS" in _why(bt, tmp_path, *oic(sha))
    sha2 = oic_env["write"](chords_doc(oic_env["genes"], cm_idx=range(0, 31)))
    assert "not read PASS" in _why(bt, tmp_path, *oic(sha))  # the row names sha2, the command another file
    oic_env["path"].write_text("{}")
    assert "not the file of that sha256" in _why(bt, tmp_path, *oic(sha2))


def test_no_r70_order_refuses(bt, oic_env, monkeypatch, tmp_path):
    text = bt.TRIES_FILE.read_text().replace(f"| {R70_FIRST} |", "| R71 (another row) |")
    p = tmp_path / "TRIES_no_r70.md"
    p.write_text(text)
    monkeypatch.setattr(bt, "TRIES_FILE", p)
    assert "exactly one row R70" in _why(bt, tmp_path, *oic(oic_env["sha"]))


def test_a_chord_no_type_can_carry_refuses_the_build(bt, oic_env, tmp_path):
    sha = oic_env["write"](chords_doc(oic_env["genes"], epi=True))
    assert "refused" in _why(bt, tmp_path, *oic(sha))


def test_no_earlier_command_reaches_the_insitu_code(bt, oic_env, monkeypatch, tmp_path):
    from vec import trend_forms

    def broken(*a, **k):
        raise AssertionError("the in-situ code was reached")
    monkeypatch.setattr(trend_forms, "insitu_step", broken)
    monkeypatch.setattr(trend_forms, "insitu_arguments", broken)
    _build(bt, tmp_path, "otk", *TF["category"])
