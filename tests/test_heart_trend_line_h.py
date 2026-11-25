"""`experiments/builds/build_test.py heart --entry organ_scope --trend-strength 1.25 --trend-form <F> --z-mode raw
--captain-ruling R55`: the heart's trend-form line (row R55 of experiments/harness/TRIES.md, the captain's order of
2026-11-24: his named exception on the heart's scope and size for OTK125, OTX125, OTW125 and OTP125, bound to their
screen, row H35), OTS125's command with one change, the direction of the trend step (vec.trend_forms).

On the synthetic stages and synthetic test board (its limit never shown) of tests/test_heart_trend_build.py, every stage
labelled released-style (the Endothelial lineage half endocardium, half vessels, its endocardium with a trend of its
own, a third of E10.5's vessels a detached island the organ cut leaves out), with a temporary ledger for TRIES.md:
each form keeps OTS125's rows, coordinates, cell number, size field and half rows and changes the expression only; each
form with its change out is OTS125's step bit for bit; each form's factors are what its rule says (the split, the whole
stage's point, the shared change out with the norm held, the slope before the reference with the norm held); the same
bytes twice; the record names the option, the form, the file, the row, the order, the row's first cell and the command
and holds no cell number; every refusal, before anything is read; the R55 check on temporary ledgers and, read by code,
on the committed ledger. Every heart command before this change gives the bytes and record it gave before (OTS15 and
OTS125 pinned here at main 397b0e1, before this change; the others pinned in the earlier lines' tests). Nothing is read
from data/ and nothing is scored."""

from __future__ import annotations

import functools
import hashlib
import json
import re
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pytest
from scipy import sparse

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from test_build_heart_iso_trim import _build, canon  # noqa: E402
from test_build_heart_organ_scope import ORGAN, _use_stage, organ_stage  # noqa: E402
from test_heart_trend_build import SECRET_MAX, STAGES, heart_test_board, install, sha, synthetic_stage  # noqa: E402
from test_heart_trend_line_d import OTS, _organ_size  # noqa: E402
from test_heart_trend_line_f import _X, _coords  # noqa: E402
from test_heart_trend_line_g import LEDGER, OTS15, OTSF, OTSF_BEFORE  # noqa: E402

ORG = ["--entry", "organ_scope"]
FORMS = ("category", "contrast", "whole", "portion")
TF = {f: [*ORG, "--trend-strength", "1.25", "--trend-form", f, "--z-mode", "raw", "--captain-ruling", "R55"]
      for f in FORMS}
R55_FIRST = ("R55 (the captain's order of 2026-11-24, ruling, no run: \"Yes to OTK125, OTX125, OTW125 and OTP125\", a "
             "synthetic stand-in; the heart builder is to take `--entry organ_scope` with `--trend-strength 1.25` and "
             "`--trend-form category`, `--trend-form contrast`, `--trend-form whole` or `--trend-form portion` only "
             "with `--captain-ruling R55`; not a try)")
R55_ROW = f"| {R55_FIRST} | heart | the trend-form line's four files |\n"
LEDGER_TF = LEDGER + R55_ROW
# sha256 of (the .h5ad, its .build.json text, the record's `canon`) that build_test.py heart wrote for row R52's OTS15
# on the organ environment of tests/test_heart_trend_line_g.py, computed with the builder of main at 397b0e1 (before
# --trend-form) by the same build as below; OTS125's (and OTS075's) are that file's OTSF_BEFORE, rechecked here.
OTS15_BEFORE = ("98258f1b6c1c9379c022f2b379ccfa168a6e4de1ac661958c4722f125a695b8c",
                "ae8f1962385731d43a1ed681472ec4d585b36d64de9d805142ce88b1317a12b8",
                "5b4738c7e852fe9ea0764e6ac300dfcc98cd865989d6071ada4b44619a18d1ac")
ISLAND_SHIFT = 80.0
LABEL_OF = {"CM": "V-CM", "Tiny": "Epi", "Untracked": "Late-a-PHM"}


@pytest.fixture(scope="module")
def bt():
    sys.path.insert(0, str(REPO / "experiments" / "builds"))
    import build_test

    return build_test


@pytest.fixture(scope="module")
def tf():
    from vec import trend_forms

    return trend_forms


def relabel(a: ad.AnnData, shift: float) -> ad.AnnData:
    """Released-style labels by lineage (CM V-CM, the core; Tiny Epi; Untracked Late-a-PHM) and the Endothelial cells
    alternately Endo (endocardium) and BEC (vessels); the endocardium cells get extra counts on genes 40 to 49 that
    fall with the stage's `shift` (a trend the vessels do not have), the log1p CP10k matrix rebuilt from the counts."""
    merged = np.asarray(a.obs["merged"]).astype(str)
    labels = np.array([LABEL_OF.get(m, "x") for m in merged], dtype=object)
    endo = np.flatnonzero(merged == "Endothelial")
    labels[endo[0::2]] = "Endo"
    labels[endo[1::2]] = "BEC"
    C = np.asarray(sparse.csr_matrix(a.layers["counts"]).todense()).astype(np.int64)
    rng = np.random.default_rng(int(10 * shift) + 17)
    C[np.ix_(endo[0::2], np.arange(40, 50))] += rng.poisson(4.0 * (2.0 - shift), (len(endo[0::2]), 10))
    a.layers["counts"] = sparse.csr_matrix(C.astype(np.int32))
    a.X = sparse.csr_matrix(np.log1p(1e4 * C / C.sum(1, keepdims=True)).astype(np.float32))
    a.obs["celltype"] = labels.astype(str)
    a.obs["core"] = merged == "CM"
    return a


def tf_organ_stage() -> ad.AnnData:
    """The synthetic E10.5, relabelled, with a third of its vessel cells moved far from the core (a detached island
    the organ cut leaves out, so the whole stage's Endothelial point is not D's pool's)."""
    a = relabel(synthetic_stage(*STAGES["E10.5.h5ad"]), STAGES["E10.5.h5ad"][2])
    bec = np.flatnonzero(np.asarray(a.obs["celltype"]) == "BEC")
    C = np.asarray(a.obsm["spatial_3D"], dtype=np.float32).copy()
    C[bec[: len(bec) // 3], 0] += ISLAND_SHIFT
    a.obsm["spatial_3D"] = C
    return a


@pytest.fixture
def ledger(bt, monkeypatch, tmp_path):
    path = tmp_path / "TRIES.md"
    path.write_text(LEDGER_TF)
    monkeypatch.setattr(bt, "TRIES_FILE", path)
    return path


def _reads(bt, monkeypatch):
    boards = bt.load_boards
    reads = []

    def load_boards():
        reads.append("boards")
        return boards()

    monkeypatch.setattr(bt, "load_boards", load_boards)
    return reads


@pytest.fixture
def organ_env(bt, monkeypatch, own_repo, ledger):
    """The earlier lines' organ environment (tests/test_heart_trend_line_g.py), with row R55 in the ledger."""
    reads: list[str] = []
    install(bt, monkeypatch.setattr, reads)
    _use_stage(bt, monkeypatch, organ_stage(), reads)
    _organ_size(bt, monkeypatch)
    more = _reads(bt, monkeypatch)
    return reads, more


@pytest.fixture
def tf_env(bt, monkeypatch, own_repo, ledger):
    """install()'s board and stages, every stage relabelled (`relabel`), E10.5 as `tf_organ_stage`, the organ size
    stand-in; the list records every stage and board read, in order."""
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


def tf_context(bt):
    """The trend-form files' context as the builder builds it on `tf_env` (D's pool of `tf_organ_stage` the reference
    input, the other stages relabelled), and the whole E10.5."""
    from vec import harness, t2_build
    from vec.io import COORD_KEY

    stages = {f: relabel(synthetic_stage(*v), v[2]) for f, v in STAGES.items()}
    ref = tf_organ_stage()
    keep, _ = t2_build.organ_scope(np.asarray(ref.obsm[COORD_KEY], dtype=float)[:, :3],
                                   ref.obs[t2_build.T2_LABEL_KEY].to_numpy(), np.asarray(ref.obs["core"], dtype=bool))
    assert 0 < keep.sum() < ref.n_obs
    roles = {"ref": ("T2/E10.5.h5ad", 10.5), "start": ("T2/E9.5.h5ad", 9.5), **bt.HEART_TREND_STAGES}
    inputs = {"ref": ref[keep].copy(), "start": stages["E9.5.h5ad"],
              **{r: stages[Path(f).name] for r, (f, _) in bt.HEART_TREND_STAGES.items()}}
    view = bt._view("test:T2h", "T2", "extrap", "heart", heart_test_board().key, roles, "ref", 12.5)
    ctx = harness.BuildContext(task=view, board=harness._candidate_board(heart_test_board()), inputs=inputs,
                               seed=harness.BUILD_SEED_OFFSET + 1)
    return ctx, ref


@pytest.fixture(scope="module")
def cw(bt):
    return tf_context(bt)


def _why(bt, tmp_path, *argv) -> str:
    out = tmp_path / ("h_dry.h5ad" if "--board" in argv else "h.h5ad")
    with pytest.raises(SystemExit) as e:
        bt.main(["heart", *argv, "--out", str(out)])
    assert not out.exists()
    return str(e.value)


# ------------------------------------------------------------------------------------------------ off by default
def test_the_option_is_off_by_default_and_names_the_designs_files(bt, tf):
    a = bt.build_parser().parse_args(["heart", *ORGAN, "--out", "runs/x/h.h5ad"])
    assert a.trend_form is None
    assert bt.HEART_TF_ROW == "R55" and bt.HEART_TF_FORMS == tf.FORMS == FORMS
    # the design's table: OTK125 category, OTX125 contrast, OTW125 whole, OTP125 portion; files T2heart_organ_ts125_<t>
    assert bt.HEART_TF_FILES == tf.FILES == {"category": "OTK125", "contrast": "OTX125", "whole": "OTW125",
                                             "portion": "OTP125"}
    assert bt.HEART_TF_TAGS == tf.TAGS == {"category": "tk", "contrast": "tx", "whole": "tw", "portion": "tp"}
    assert {f: bt.HEART_TF_COMMANDS[f] for f in FORMS} == {f: "heart " + " ".join(TF[f]) for f in FORMS}
    assert bt.HEART_TF_STRENGTH == 1.25 and bt.HEART_TF_STRENGTH_CHOICE in bt.HEART_FURTHER_CHOICES
    assert bt.HEART_FURTHER_COMMANDS[1.25].format(row="R49") == "heart " + " ".join(OTSF[1.25])  # OTS125's own
    assert bt.HEART_FURTHER_ROW == "R49" and bt.HEART_OTS15_ROW == "R52"


# ------------------------------------------------------------------------------------------------ old paths unchanged
def test_ots15_keeps_its_bytes_and_record(bt, organ_env, tmp_path):
    out, rec, text = _build(bt, tmp_path, "ots15", *OTS15)
    assert sha(out) == OTS15_BEFORE[0]
    assert hashlib.sha256(text.encode()).hexdigest() == OTS15_BEFORE[1] or canon(text) == OTS15_BEFORE[2]
    assert "R55" not in text and "trend-form" not in text and "form" not in rec["trend_strength"]


@pytest.mark.parametrize("s", [0.75, 1.25])
def test_row_r49s_files_keep_their_bytes_and_record(bt, organ_env, tmp_path, s):
    out, rec, text = _build(bt, tmp_path, f"f{s:g}", *OTSF[s])
    assert sha(out) == OTSF_BEFORE[s][0]
    assert hashlib.sha256(text.encode()).hexdigest() == OTSF_BEFORE[s][1] or canon(text) == OTSF_BEFORE[s][2]
    assert rec["trend_strength"]["captain_ruling"] == "R49" and "R55" not in text


def test_no_earlier_command_reaches_the_forms_code(bt, organ_env, monkeypatch, tmp_path):
    """vec.trend_forms is never called without --trend-form: every earlier heart command builds with it broken."""
    from vec import trend_forms

    def broken(*a, **k):
        raise AssertionError("the forms' code was reached")
    monkeypatch.setattr(trend_forms, "form_step", broken)
    monkeypatch.setattr(trend_forms, "form_arguments", broken)
    for name, argv in (("d", ORGAN), ("o5", OTS[0.5]), ("o10", OTS[1.0]), ("o075", OTSF[0.75]), ("o125", OTSF[1.25]),
                       ("o15", OTS15)):
        _build(bt, tmp_path, name, *argv)


# ------------------------------------------------------------------------------------------------ the files
@pytest.mark.parametrize("form", FORMS)
def test_each_form_keeps_ots125s_rows_coordinates_count_size_and_half(bt, tf_env, tmp_path, form):
    """Rows, coordinates, cell number, size field, D's own record fields and, after halve.py at seed 1000001, the
    half's rows and coordinates are OTS125's; the expression is not."""
    sys.path.insert(0, str(REPO / "experiments" / "builds"))
    import halve

    o, o_rec, _ = _build(bt, tmp_path, "ots125", *OTSF[1.25])
    t, t_rec, _ = _build(bt, tmp_path, form, *TF[form])
    A, B = ad.read_h5ad(o), ad.read_h5ad(t)
    assert A.n_obs == B.n_obs and list(A.obs_names) == list(B.obs_names) and list(A.var_names) == list(B.var_names)
    np.testing.assert_array_equal(_coords(o), _coords(t))
    for k in ("n_cells", "seed", "count_decision", "trim", "size", "z_mode", "board", "dry_run", "depth_step",
              "scope_option"):
        assert t_rec[k] == o_rec[k], k
    ho, ht = tmp_path / "ots125_half.h5ad", tmp_path / f"{form}_half.h5ad"
    assert halve.DEFAULT_SEED == 1000001
    halve.halve(str(o), str(ho))
    halve.halve(str(t), str(ht))
    HA, HB = ad.read_h5ad(ho), ad.read_h5ad(ht)
    assert list(HA.obs_names) == list(HB.obs_names) and HA.n_obs == HB.n_obs
    np.testing.assert_array_equal(_coords(ho), _coords(ht))
    assert not np.array_equal(_X(o), _X(t)) and not np.array_equal(_X(ho), _X(ht))


def test_the_four_forms_differ_from_each_other(bt, tf_env, tmp_path):
    X = {f: _X(_build(bt, tmp_path, f, *TF[f])[0]) for f in FORMS}
    for i, f in enumerate(FORMS):
        for g in FORMS[i + 1:]:
            assert not np.array_equal(X[f], X[g]), (f, g)


@pytest.mark.parametrize("form", FORMS)
def test_the_same_bytes_twice(bt, tf_env, tmp_path, form):
    a, _, ta = _build(bt, tmp_path, "a", *TF[form])
    b, _, tb = _build(bt, tmp_path, "b", *TF[form])
    assert a.read_bytes() == b.read_bytes() and ta == tb


@pytest.mark.parametrize("form", FORMS)
def test_the_record_names_option_form_file_row_order_and_command_and_no_cell_number(bt, tf_env, tmp_path, capsys,
                                                                                    form):
    _, rec, text = _build(bt, tmp_path, form, *TF[form])
    ts = rec["trend_strength"]
    assert ts["option"] == "--trend-form" and ts["form"] == form and ts["strength"] == 1.25
    assert ts["file"] == bt.HEART_TF_FILES[form] and ts["tag"] == bt.HEART_TF_TAGS[form]
    assert ts["base_file"] == "OTS125" and ts["captain_ruling"] == "R55"
    assert ts["order"] == bt.HEART_TF_ORDER.format(file=bt.HEART_TF_FILES[form], form=form)
    assert ts["order"].startswith("experiments/harness/TRIES.md row R55, the captain's order of 2026-11-24")
    assert "bound to their screen" in ts["order"] and "row H35" in ts["order"] and "nothing else changed" in ts["order"]
    assert ts["order_row"] == R55_FIRST  # the row's first cell as the ledger holds it
    assert ts["command"] == "heart " + " ".join(TF[form]) == bt.HEART_TF_COMMANDS[form]
    assert ts["step"] == bt.HEART_TF_STEP and ts["budget_days"] == pytest.approx(0.5)
    assert ts["test_budget_days"] == pytest.approx(0.625)
    assert ts["trend"]["stages"] == ["e8_25_late", "e8_5", "e8_75", "start", "ref"]
    assert set(ts) >= set(bt.HEART_TF_RECORD[form]) and {t["status"] for t in ts["types"]} <= {
        "shifted", "skipped: fewer than 20 cells at a stage"}
    assert rec["entry"] == f"trend strength 1.25 form {form} organ_scope (organ)" and "expression" not in rec
    assert rec["scope_option"]["command"] == bt.HEART_ORGAN_SCOPE and rec["n_cells"] == "max_cells"
    assert "R52" not in text and text.count("R49") == 1 and "row R49's) with one change" in ts["order"]

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
    # the board's limit nowhere as a number (a float's digits may hold the same three digits)
    assert not re.search(rf"(?<![\d.]){SECRET_MAX}(?![\d.])", text + capsys.readouterr().out)


@pytest.mark.parametrize("form", FORMS)
def test_the_reads_are_ots125s(bt, tf_env, tmp_path, form):
    _build(bt, tmp_path, "ots125", *OTSF[1.25])
    want = list(tf_env)
    assert want == ["boards", "E10.5.h5ad", "E9.5.h5ad", "E8.75.h5ad", "E8.5.h5ad", "E8.25_late.h5ad", "boards"]
    tf_env.clear()
    _build(bt, tmp_path, form, *TF[form])
    assert tf_env == want


# ------------------------------------------------------------------------------------------------ the step
def test_the_banked_step_is_the_builders_and_at_strength_1_the_harness_b0_trend(bt, tf, cw):
    from vec import harness

    ctx, _ = cw
    st, rec = bt._trend_strength_step(ctx, 1.25)
    mine, _, _, mrec = tf.form_step(ctx, 1.25, None)
    np.testing.assert_array_equal(st.X, mine.X)
    np.testing.assert_array_equal(st.coords, mine.coords)
    assert mrec["form"] == "banked" and mrec["trend"]["lineages"] == rec["lineages"]
    one = tf.form_step(ctx, 1.0, None)[0]
    b0 = harness._b0_t2(ctx, trend=True)
    np.testing.assert_array_equal(one.X, b0.X)


@pytest.mark.parametrize("form", FORMS)
def test_each_form_with_its_change_out_is_ots125s_step_bit_for_bit(bt, tf, cw, form):
    ctx, whole = cw
    st, _ = bt._trend_strength_step(ctx, 1.25)
    z = tf.form_step(ctx, 1.25, form, whole, zero=True)[0]
    np.testing.assert_array_equal(z.X, st.X)
    np.testing.assert_array_equal(z.coords, st.coords)
    za, ba = tf.form_arguments(ctx, 1.25, form, whole, zero=True), tf.form_arguments(ctx, 1.25, None)
    assert za.damping == ba.damping and sorted(za.mapping) == sorted(ba.mapping)
    for k in ba.mapping:
        assert np.array_equal(za.mapping[k], ba.mapping[k]), k
    # and with the change in, the step is another
    on = tf.form_step(ctx, 1.25, form, whole)[0]
    assert not np.array_equal(on.X, st.X)
    np.testing.assert_array_equal(on.coords, st.coords)


def test_category_splits_the_endothelial_lineage_only(bt, tf, cw):
    from dataclasses import replace

    from vec import harness

    ctx, _ = cw
    ba, ka = tf.form_arguments(ctx, 1.25, None), tf.form_arguments(ctx, 1.25, "category")
    np.testing.assert_array_equal(ka.rows, ba.rows)
    split = {"Endothelial/endocardium", "Endothelial/vascular_endothelium"}
    assert set(ka.like.obs["b0_type"]) == (set(ba.like.obs["b0_type"]) - {"Endothelial"}) | split
    assert ka.record["split_types"] == sorted(split) and not any(ka.record["collapsed"].values())
    assert sorted(ka.mapping) == sorted((set(ba.mapping) - {"Endothelial"}) | split)
    for k in set(ba.mapping) - {"Endothelial"}:  # every other lineage's factor is OTS125's, exactly
        assert np.array_equal(ka.mapping[k], ba.mapping[k]), k
    # each part's factor is the trend over its own cells at every stage: keyed by label, the same arithmetic
    keyed = {r: harness._shallow(a, part=tf.split_types(a, None)) for r, a in ctx.inputs.items()}
    logf, _ = harness.b0_trend_factors(replace(ctx, inputs=keyed), list(ka.like.var_names), sorted(split), key="part")
    for k in split:
        np.testing.assert_array_equal(ka.mapping[k], 1.25 * logf[k])
        assert not np.allclose(ka.mapping[k], ba.mapping["Endothelial"])
    # the drawn Endothelial cells are typed by their own labels
    lab = np.asarray(ka.like.obs["celltype"]).astype(str)
    typ = np.asarray(ka.like.obs["b0_type"]).astype(str)
    assert set(typ[lab == "Endo"]) == {"Endothelial/endocardium"}
    assert set(typ[lab == "BEC"]) == {"Endothelial/vascular_endothelium"}
    assert set(typ[lab == "V-CM"]) == {"CM"} and set(typ[lab == "Epi"]) == {"Tiny"}
    assert set(typ[lab == "Late-a-PHM"]) <= {"Untracked@end"}  # outside the organ cut: none is drawn here


def test_whole_takes_the_reference_point_on_the_whole_stage(bt, tf, cw):
    from dataclasses import replace

    from vec import harness

    ctx, whole = cw
    ba, wa = tf.form_arguments(ctx, 1.25, None), tf.form_arguments(ctx, 1.25, "whole", whole)
    logf, _ = harness.b0_trend_factors(replace(ctx, inputs={**ctx.inputs, "ref": whole}),
                                       list(ba.like.var_names), sorted(set(ba.like.obs["b0_type"])))
    assert sorted(wa.mapping) == sorted(logf)
    for k in logf:
        np.testing.assert_array_equal(wa.mapping[k], 1.25 * logf[k])
    # the core is always in D's pool: CM's factor is OTS125's; the Endothelial point moves (the island)
    assert np.array_equal(wa.mapping["CM"], ba.mapping["CM"])
    assert not np.allclose(wa.mapping["Endothelial"], ba.mapping["Endothelial"])
    assert wa.like is ba.like or np.array_equal(wa.like.obs["b0_type"], ba.like.obs["b0_type"])
    assert wa.end.n_obs == ctx.inputs["ref"].n_obs < whole.n_obs  # momentum_shift stays on D's pool
    with pytest.raises(ValueError, match="needs the whole reference stage"):
        tf.form_arguments(ctx, 1.25, "whole")


def test_contrast_is_the_trend_less_the_whole_stages_shared_change_with_its_norm_held(bt, tf, cw):
    from dataclasses import replace

    from vec import harness

    ctx, whole = cw
    ba, xa = tf.form_arguments(ctx, 1.25, None), tf.form_arguments(ctx, 1.25, "contrast", whole)
    wtypes = sorted(set(tf.merged_types(whole, "end")))
    wlog, _ = harness.b0_trend_factors(replace(ctx, inputs={**ctx.inputs, "ref": whole}),
                                       list(ba.like.var_names), wtypes)
    assert sorted(wlog) == ["CM", "Endothelial", "Tiny"] and xa.record["shared_lineages"] == sorted(wlog)
    sbar = np.mean(np.vstack([wlog[k] for k in sorted(wlog)]), axis=0)
    assert xa.record["shared_norm"] == pytest.approx(float(np.linalg.norm(sbar)))
    assert sorted(xa.mapping) == sorted(ba.mapping)
    for k, v in ba.mapping.items():
        c = v / 1.25 - sbar
        got = xa.mapping[k]
        assert float(np.linalg.norm(got)) == pytest.approx(float(np.linalg.norm(v)), rel=1e-12)  # the norm held
        assert harness._cosine(got, c) == pytest.approx(1.0, abs=1e-12)  # the direction: the trend less s-bar
        assert harness._cosine(got, v) < 0.999
    with pytest.raises(ValueError, match="needs the whole reference stage"):
        tf.form_arguments(ctx, 1.25, "contrast")


def test_portion_is_the_slope_before_the_reference_with_its_norm_held(bt, tf, cw):
    from vec import harness
    from vec.core import counts_of, mean_cp10k

    ctx, _ = cw
    ba, pa = tf.form_arguments(ctx, 1.25, None), tf.form_arguments(ctx, 1.25, "portion")
    assert sorted(pa.mapping) == sorted(ba.mapping) and pa.record["banked_factor_kept"] == []
    genes = list(ba.like.var_names)
    by = {r["lineage"]: r for r in ba.record["trend"]["lineages"]}
    for k, v in ba.mapping.items():
        used = [(r, d) for r, d in zip(by[k]["stages"], by[k]["days"], strict=True) if r != "ref"]
        assert pa.record["lineages"][k]["stages"] == [r for r, _ in used] and len(used) >= 2
        y = []
        for r, _ in used:
            a = ctx.inputs[r]
            rows = np.flatnonzero(np.asarray(a.obs["merged"]).astype(str) == k)
            y.append(np.log(mean_cp10k(counts_of(a, genes)[rows]) + 0.5))
        p = harness.trend_log_factor(y, [d for _, d in used], 0.5)
        got = pa.mapping[k]
        assert float(np.linalg.norm(got)) == pytest.approx(float(np.linalg.norm(v)), rel=1e-12)
        assert harness._cosine(got, p) == pytest.approx(1.0, abs=1e-12)
    assert pa.record["lineages"]["Tiny"]["stages"] == ["e8_75", "start"]  # its two usable stages before E10.5


def test_the_damping_and_the_stream_are_ots125s_in_every_form(bt, tf, cw, monkeypatch):
    """momentum_shift gets OTS125's damping (B0_DAMPING x 1.25), key, held genes and stream in every form."""
    from vec import harness, t1_scope

    ctx, whole = cw
    seen = []
    own = t1_scope.momentum_shift

    @functools.wraps(own)  # b0_trend_factors reads momentum_shift's defaults from its signature
    def spy(*a, **k):
        # keep_depth and every other default are momentum_shift's own: each cell's expected depth is kept (the counts
        # of every moved cell are drawn afresh, so depth is held in expectation only)
        assert set(k) == {"key", "seed", "logf_by_type"} and len(a) == 5
        seen.append((a[3], a[4], k.get("key"), k.get("seed")))
        return own(*a, **k)
    monkeypatch.setattr(t1_scope, "momentum_shift", spy)
    bt._trend_strength_step(ctx, 1.25)
    for f in FORMS:
        tf.form_step(ctx, 1.25, f, whole)
    assert len(seen) == 5 and len(set(map(repr, seen))) == 1
    assert seen[0] == (harness.B0_DAMPING * 1.25, (), "b0_type", [ctx.seed, 3])


def test_a_direction_of_no_length_cannot_carry_the_norm(tf):
    v = np.array([0.3, -0.2, 0.1])
    np.testing.assert_array_equal(tf._held(v, v.copy(), "x"), v)  # c == v gives v back bit for bit
    with pytest.raises(ValueError, match="its norm cannot be held"):
        tf._held(v, np.zeros(3), "x")
    with pytest.raises(ValueError, match="its norm cannot be held"):
        tf._held(v, np.array([np.nan, 0.0, 0.0]), "x")
    with pytest.raises(ValueError, match="unknown trend form"):
        tf.form_arguments(None, 1.25, "banked2")


def test_an_unmapped_endothelial_label_refuses_the_category_form_only(bt, tf_env, monkeypatch, tmp_path):
    load = bt.load_stage

    def odd(p, **kw):
        a = load(p, **kw)
        if Path(p).name == "E8.5.h5ad":
            lab = np.asarray(a.obs["celltype"]).astype(object)
            lab[np.flatnonzero(lab == "Endo")[0]] = "Endo-unmapped"
            a.obs["celltype"] = lab.astype(str)
        return a
    monkeypatch.setattr(bt, "load_stage", odd)
    msg = _why(bt, tmp_path, *TF["category"])
    assert "without a T2 chamber category" in msg and "Endo-unmapped" in msg
    for f in ("contrast", "whole", "portion"):
        _build(bt, tmp_path, f, *TF[f])


# ------------------------------------------------------------------------------------------------ refusals
NEW_REFUSED = [
    [*ORG, "--trend-strength", "1.25", "--trend-form", "category", "--z-mode", "raw"],  # no row
    *([*ORG, "--trend-strength", "1.25", "--trend-form", "contrast", "--z-mode", "raw", "--captain-ruling", r]
      for r in ("R49", "R52", "R54", "R56", "R550", "r55", "R55 ", " R55", "R55a", "")),
    *([*ORG, "--trend-strength", v, "--trend-form", "whole", "--z-mode", "raw", "--captain-ruling", "R55"]
      for v in ("0.75", "1.5", "1.0", "0.5", "1.250", "1.25e0", "1.2", "+1.25")),
    [*ORG, "--trend-form", "whole", "--z-mode", "raw", "--captain-ruling", "R55"],  # no strength
    *(["--entry", e, "--trend-strength", "1.25", "--trend-form", "portion", "--z-mode", "raw", "--captain-ruling",
       "R55"] for e in ("s1", "s2", "s2_trim", "core_trim")),
    [*ORG, "--trend-strength", "1.25", "--trend-form", "portion", "--captain-ruling", "R55"],
    [*ORG, "--trend-strength", "1.25", "--trend-form", "portion", "--z-mode", "recorded", "--captain-ruling", "R55"],
    *([*TF["category"], *f.split()] for f in ("--expression copy", "--trend", "--depth", "--herd", "--cfm")),
    [*TF["contrast"], "--seed", "1000002"],
    *([*ORG, "--trend-strength", "1.25", "--trend-form", v, "--z-mode", "raw", "--captain-ruling", "R55"]
      for v in ("Category", "categories", "", "whole ", "contrast,whole", "banked", "ks", "count-kernel")),
]


@pytest.mark.parametrize("argv", NEW_REFUSED)
def test_every_other_use_is_refused_before_anything_is_read(bt, tf_env, tmp_path, argv):
    msg = _why(bt, tmp_path, *argv)
    assert "refused" in msg or "needs --captain-ruling R55" in msg
    assert tf_env == []


def test_the_refusals_say_why(bt, tf_env, ledger, tmp_path):
    cmd = "heart " + " ".join(TF["category"])
    assert _why(bt, tmp_path, *TF["category"][:-2]) == (
        "--trend-form category is the captain's order of row R55 alone and needs --captain-ruling R55")
    for r in ("R49", "R54", "r55", ""):
        assert _why(bt, tmp_path, *TF["category"][:-1], r) == (
            f"--trend-form category is refused with --captain-ruling {r}: only row R55 orders it")
    for argv in ([*ORG, "--trend-strength", "0.75", "--trend-form", "category", "--z-mode", "raw",
                  "--captain-ruling", "R55"], [*TF["category"], "--expression", "copy"], [*TF["category"], "--trend"]):
        assert _why(bt, tmp_path, *argv) == (
            f"--trend-form category is refused except as {cmd}: OTS125's line only (--entry organ_scope, "
            "--trend-strength 1.25, --z-mode raw), and with --expression copy, --trend, --depth, --herd and --cfm")
    for extra in (["--board", "T2:heart:val_extrap"], ["--seed", "1000002"]):
        assert _why(bt, tmp_path, *TF["category"], *extra) == (
            f"--trend-form category is refused on another board or seed: it builds {cmd} at the default seed on the "
            "test board only (no dry run)")
    assert _why(bt, tmp_path, *ORG, "--trend-strength", "1.25", "--trend-form", "ks", "--z-mode", "raw") == (
        "--trend-form ks is refused: it takes category, contrast, whole, portion only")
    ledger.write_text(LEDGER)
    assert _why(bt, tmp_path, *TF["whole"]) == (
        "--trend-form whole is refused: experiments/harness/TRIES.md must hold exactly one row R55 whose first cell is "
        "the captain's order (\"ruling, no run\") naming `--entry organ_scope`, `--trend-strength 1.25`, "
        "`--trend-form whole` and `--captain-ruling R55`, each as a whole backticked string")
    assert tf_env == []


def test_every_earlier_value_and_row_meets_its_old_refusal_with_r55_in_the_ledger(bt, tf_env, tmp_path):
    """Without --trend-form, --captain-ruling R55 is another row to every earlier line, refused as before."""
    assert _why(bt, tmp_path, *OTSF[1.25][:-1], "R55") == (
        "--entry organ_scope --trend-strength 1.25 is refused: experiments/harness/TRIES.md must hold exactly one row "
        "R55 whose first cell is the captain's order (\"ruling, no run\") naming `--entry organ_scope`, "
        "`--trend-strength 0.75`, `--trend-strength 1.25` and `--captain-ruling R55`")
    assert _why(bt, tmp_path, *OTS15[:-1], "R55") == (
        "--entry organ_scope --trend-strength 1.5 is refused with --captain-ruling R55: only row R52 orders it")
    assert _why(bt, tmp_path, *ORGAN, "--captain-ruling", "R55").startswith(
        "heart --captain-ruling goes with --trend-strength only")
    assert tf_env == []


def _use(bt, monkeypatch, path: Path, text: str) -> None:
    path.write_text(text)
    monkeypatch.setattr(bt, "TRIES_FILE", path)


def test_the_r55_check_on_ledgers_of_the_required_form(bt, monkeypatch, tmp_path):
    _use(bt, monkeypatch, tmp_path / "ok.md", LEDGER_TF)
    assert all(bt._heart_tf_ordered(f) for f in FORMS) and bt._heart_further_ordered("R49")
    _use(bt, monkeypatch, tmp_path / "only.md", R55_ROW)
    assert all(bt._heart_tf_ordered(f) for f in FORMS) and not bt._heart_further_ordered("R49")
    # a row naming two of the forms orders those two only
    two = R55_FIRST.replace("`--trend-form category`, `--trend-form contrast`, ", "")
    _use(bt, monkeypatch, tmp_path / "two.md", f"| {two} | heart | x |\n")
    assert [f for f in FORMS if bt._heart_tf_ordered(f)] == ["whole", "portion"]
    assert not bt._heart_tf_ordered("banked") and not bt._heart_tf_ordered("")


@pytest.mark.parametrize("drop", ["`--entry organ_scope`", "`--trend-strength 1.25`", "`--trend-form category`",
                                  "`--captain-ruling R55`", "ruling, no run", "the captain's order"])
def test_the_r55_check_fails_without_each_required_string(bt, monkeypatch, tmp_path, drop):
    first = R55_FIRST.replace(drop, "something else")
    _use(bt, monkeypatch, tmp_path / "t.md", LEDGER.replace("| R27 (other)", f"| {first} | heart | x |\n| R27 (other)"))
    assert not bt._heart_tf_ordered("category")
    if drop == "`--trend-form category`":
        assert bt._heart_tf_ordered("contrast")


@pytest.mark.parametrize("text", [
    "| R27 (other) | T1 | x |\n", LEDGER, LEDGER_TF + R55_ROW,  # none, none, two
    LEDGER_TF.replace("| R55 (", "| R550 ("), LEDGER_TF.replace("| R55 (", "| R56 ("),
    LEDGER_TF.replace("`--captain-ruling R55`", "`--captain-ruling R49`"),
    LEDGER_TF.replace("`--trend-strength 1.25`", "`--trend-strength 1.250`"),
    LEDGER_TF.replace("`--trend-strength 1.25`", "--trend-strength 1.25"),
    LEDGER_TF.replace("`--trend-form category`", "`--trend-form category2`"),
    LEDGER_TF.replace("`--trend-form category`", "`--trend-form category --z-mode raw`"),
    LEDGER_TF.replace("`--trend-form category`", "`x --trend-form category`"),
    LEDGER_TF.replace("`--entry organ_scope`", "`--entry organ_scope --z-mode raw`"),
    LEDGER_TF.replace("ruling, no run", "awaiting a ruling"),
    LEDGER_TF.replace("R55 (the captain's order", "R55 (the captain's named exception"),
    # the strings outside the first cell do not count
    LEDGER_TF.replace(R55_ROW, R55_ROW.replace("`--trend-form category`, ", "").replace(
        "| heart |", "| `--trend-form category` |"))])
def test_the_category_build_is_refused_on_a_ledger_without_the_order(bt, tf_env, ledger, tmp_path, text):
    ledger.write_text(text)
    msg = _why(bt, tmp_path, *TF["category"])
    assert "must hold exactly one row R55" in msg
    assert tf_env == []


def test_the_committed_ledger_holds_r55_in_the_builders_form():
    """experiments/harness/TRIES.md as committed, read by code only: row R55's first cell is the captain's order naming
    the four forms' files in the form the builder takes (rows R29's, R49's and R52's too, the launcher's identity
    steps)."""
    sys.path.insert(0, str(REPO / "experiments" / "builds"))
    import build_test as bt

    assert bt.TRIES_FILE == REPO / "experiments" / "harness" / "TRIES.md"
    assert all(bt._heart_tf_ordered(f) for f in FORMS)
    assert bt._heart_further_ordered(bt.HEART_FURTHER_ROW) and bt._heart_d_line_ordered() and bt._heart_ots15_ordered()


def test_the_test_boards_own_name_is_its_default(bt, tf_env, tmp_path):
    out, rec, _ = _build(bt, tmp_path, "named", *TF["portion"], "--board", "T2:heart:test_extrap")
    default, _, _ = _build(bt, tmp_path, "default", *TF["portion"])
    assert out.read_bytes() == default.read_bytes() and rec["dry_run"] is False


def test_the_records_json_text_round_trips(bt, tf_env, tmp_path):
    _, rec, text = _build(bt, tmp_path, "c", *TF["category"])
    assert json.loads(text) == rec and rec["trend_strength"]["collapsed"] == {
        "ref": False, "start": False, "e8_75": False, "e8_5": False, "e8_25_late": False, "draw": False}
