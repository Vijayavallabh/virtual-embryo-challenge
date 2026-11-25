"""Row H35 (experiments/harness/h35_heart_trend_forms.py): synthetic arrays, hand-made tables and temporary files, and
committed result files of rows H27 and H29b read for arithmetic only (G1's dry run). Nothing here loads a competition
file, scores a real cell or reads a figure of row H34; the repository's ledger is read by code only where a test says
so. The launcher's tests are in tests/test_run_h35_launcher.py."""

from __future__ import annotations

import ast
import copy
import hashlib
import importlib.util
import inspect
import itertools
import json
import re
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from vec import harness, prep, t2_build, trend_forms

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
SCRIPT = REPO_ROOT / "experiments" / "harness" / "h35_heart_trend_forms.py"
COMMIT = "c" * 40
KEY1, KEY2 = "T2h/<=E9.5->E10.5", "T2h/<=E8.75->E10.5"
FORMS = ("category", "contrast", "whole", "portion")
METRICS = ("de_score", "de_direction", "mmd_u", "variogram", "neighborhood_mmd", "d2_shape", "occupancy_dice",
           "scale_log_ratio")
GROUPS = ("expression_change", "cell_state", "shape_growth", "local_organisation")
FIRST = "H35 (the counted screen of row R55's trend forms: ten tries; G0 first, target-free, not a try)"
LABEL_OF = {"CM": "V-CM", "Untracked": "NCC", "Tiny": "Epi"}  # cm (the core), surrounding, epicardium


@pytest.fixture(scope="module")
def h35():
    spec = importlib.util.spec_from_file_location("h35_heart_trend_forms_module", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def _labelled(a):
    """T2 labels by lineage, the Endothelial cells alternately Endo (endocardium) and BEC (vessels); the core is CM."""
    merged = np.asarray(a.obs["merged"]).astype(str)
    labels = np.array([LABEL_OF.get(m, "x") for m in merged], dtype=object)
    endo = np.flatnonzero(merged == "Endothelial")
    labels[endo[0::2]], labels[endo[1::2]] = "Endo", "BEC"
    a.obs["celltype"] = labels.astype(str)
    a.obs["core"] = merged == "CM"
    return a


@pytest.fixture(scope="module")
def synth():
    """test_b0_trend's synthetic T2 stages, labelled. The tiny lineage has enough cells at `start` and at the reference
    only on the first transition: the trend mapping leaves it out there and it takes B0's own factor in every form."""
    spec = importlib.util.spec_from_file_location("b0_trend_tests_for_h35", REPO_ROOT / "tests" / "test_b0_trend.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    stages = {"T2/E9.5.h5ad": mod._stage(1600, seed=2, shift=0.0, tiny=160),
              "T2/E8.75.h5ad": mod._stage(1100, seed=3, shift=0.8, tiny=120),
              "T2/E8.5.h5ad": mod._stage(900, seed=5, shift=1.6, tiny=10),
              "T2/E8.25_late.h5ad": mod._stage(800, seed=9, shift=2.0, tiny=10),
              "T2/E10.5.h5ad": mod._stage(1400, seed=7, shift=-0.5, tiny=50)}
    for a in stages.values():
        _labelled(a)
    return mod, stages


ROLES = {KEY1: {"ref": "T2/E9.5.h5ad", "start": "T2/E8.75.h5ad", "earlier": "T2/E8.5.h5ad"},
         KEY2: {"ref": "T2/E8.75.h5ad", "start": "T2/E8.5.h5ad", "earlier": "T2/E8.25_late.h5ad"}}


def _arrays(a):
    return (np.asarray(a.obsm["spatial_3D"], float), a.obs["celltype"].to_numpy(), np.asarray(a.obs["core"], bool))


def _geoms(h34, stages, files):
    sp = h34.h27.sp
    return [sp.geometry_from_arrays(harness.stem(f), day, _arrays(stages[f])[0], _arrays(stages[f])[2])
            for f, day in files]


@pytest.fixture(scope="module")
def dform(h35, synth):
    """D's form on the synthetic stages of the first transition: the pool, the whole reference stage, D_f's size and a
    context maker (row H34's test's own construction)."""
    mod, stages = synth
    h34 = h35.h34
    inputs = {r: stages[f] for r, f in ROLES[KEY1].items()}
    whole = inputs["ref"]
    PC, labels, pcore = _arrays(whole)
    masks, _ = h34.h29b.organ_pools(PC, labels, pcore)
    ratios = {p: t2_build.organ_core_rms_ratio(PC, masks[p], pcore) for p in h34.h29b.POOLS}
    size = h34.h27.analogue_size_result(_geoms(h34, stages, [("T2/E8.5.h5ad", 8.5), ("T2/E8.75.h5ad", 8.75),
                                                              ("T2/E9.5.h5ad", 9.5)]), 10.5)
    ft = h34.h27.forecast_targets(size)
    ins = {**inputs, "ref": whole[masks["organ"]].copy()}

    def ctx(seed=4, view=None):
        return mod._ctx(ins, view=view, seed=harness.BUILD_SEED_OFFSET + seed)
    target = h34.h29b.d_targets(ft, {h34.h29b.CORE_FORM: float("nan")}, ratios)[h34.D_F]
    return SimpleNamespace(mod=mod, inputs=ins, whole=whole, masks=masks, ctx=ctx, target=target)


# ------------------------------------------------------------------------------------------------ constants, imports
def test_the_fixed_constants_and_names(h35):
    assert h35.ROW == "H35" and h35.RULING_ROW == "R55" and h35.FIRST_CELL == FIRST and h35.N_CELLS == 11
    assert h35.FORMS == FORMS == trend_forms.FORMS and h35.STRENGTH == 1.25
    assert h35.FILES == {"category": "OTK125", "contrast": "OTX125", "whole": "OTW125", "portion": "OTP125"}
    assert h35.TASK_KEYS == (KEY1, KEY2) and tuple(h35.SEEDS) == tuple(range(10))
    assert h35.WORLDS == (("own", ("own",)), ("CLEAN", ("radial_0", "thin_0.5", "thin_0")))
    # ten tries: the control and the four forms on the two transitions; the copy is a reference
    assert len((h35.CONTROL, *h35.FORMS)) * len(h35.TASK_KEYS) == 10 and h35.COPY == "copy"
    assert (h35.DROP_BAR, h35.MIN_SEEDS, h35.IDLE_BAR) == (0.30, 8, 0.10)
    assert (h35.G0_IDLE_SHARE, h35.G0_WILD_RATIO) == (0.02, 3.0)
    assert h35.G0_EXIT == {"PASS": 0, "STOP": 3, "VOID": 5} and h35.GATE_EXIT == {"PASS": 0, "FAIL": 3}
    assert h35.G0_ARGV == ("heart", "--entry", "organ_scope", "--trend-strength", "1.25", "--trend-form", "category",
                           "--z-mode", "raw", "--captain-ruling", "R55")
    assert h35.RESULTS == REPO_ROOT / "experiments" / "harness" / "results" / "h35_heart_trend_forms"
    assert h35.PRACTICE_BUDGET == {KEY1: 0.375, KEY2: 0.125} and h35.TOL == 1e-9


def test_row_h34s_code_is_imported_not_copied_and_no_h34_figure_is_read(h35):
    assert Path(h35.h34.__file__).name == "h34_heart_organ_trend.py"
    assert h35.h28 is h35.h34.h28 and h35.h29b is h35.h34.h29b and h35.TASK_KEYS is h35.h34.TASK_KEYS
    src = SCRIPT.read_text()
    tree = ast.parse(src)
    defined = {n.name for n in ast.walk(tree) if isinstance(n, ast.FunctionDef)}
    for theirs in ("size_and_half", "score_arms", "build_control", "draw_record", "record_differs", "pool_draw",
                   "means", "paired", "committed_ok", "load_h29b_reading", "load_h27_rows", "trend_expression"):
        assert theirs not in defined, theirs
    for used in ("h34.load_transition(", "h34.context(", "h34.size_and_half(", "h34.score_arms(", "h34.build_control(",
                 "h34.draw_record(", "h34.record_differs(", "h34.pool_draw(", "h34.means(", "h34.paired(",
                 "h34.committed_ok(", "h34.load_h29b_reading(", "h34.load_h27_rows(", "h34.compare(",
                 "trend_forms.form_step(", "trend_forms.form_arguments(", "trend_forms.applied(",
                 "trend_forms.g0_statistics("):
        assert used in src, used
    # row H34's figures: its results folder, parts, gate and reading are never named or loaded
    for bad in ("h34.RESULTS", "h34_heart_organ_trend/", "h34.load_parts", "h34.load_gate", "h34.reading(",
                "h34.gate_scores", "h34.json", "results/h34"):
        assert bad not in src, bad


def test_vec_is_this_checkouts_own(h35):
    body = ast.parse(SCRIPT.read_text()).body
    first_vec = next(i for i, n in enumerate(body) if isinstance(n, ast.ImportFrom) and (n.module or "")
                     .split(".")[0] == "vec")
    insert = next(i for i, n in enumerate(body) if isinstance(n, ast.Expr)
                  and ast.unparse(n).startswith("sys.path.insert(0, str(HERE / 'src'))"))
    assert insert < first_vec
    assert h35.HERE == REPO_ROOT and h35.vec_home() == REPO_ROOT


# ------------------------------------------------------------------------------------------------ the row
def _row(first=FIRST, outcome="pending", cells=11, board="heart"):
    mid = " | ".join(["x"] * (cells - 3))
    return f"| {first} | {board} | {mid} | {outcome} |\n"


def _led(*rows):
    return "# Tries\n\n| ID | Board |\n|---|---|\n| H34 (one row) | heart | x |\n" + "".join(rows)


def test_the_row_must_be_exactly_one_of_eleven_cells_with_the_fixed_first_cell(h35):
    ok = _led(_row())
    assert h35.ledger_row(ok)[0] == FIRST and len(h35.ledger_row(ok)) == 11
    assert h35.row_formed(ok) and h35.row_open(ok) and h35.screen_outcome(ok) is None
    bad = {"none": _led(), "two": _led(_row(), _row()), "ten": _led(_row(cells=10)), "twelve": _led(_row(cells=12)),
           "pipe_in_a_cell": _led(_row(outcome="||f|| pending")),
           "other_first": _led(_row(first=FIRST.replace("ten tries", "eight tries"))),
           "bold_first": _led(_row(first=f"**{FIRST}**")), "other_id": _led(_row(first=FIRST.replace("H35", "H36"))),
           "longer_id": _led(_row().replace("| H35 (", "| H350 (")), "unclosed": _led(_row().rstrip(" |\n") + "\n")}
    for name, text in bad.items():
        assert not h35.row_formed(text), name
        assert h35.screen_outcome(text) is None and not h35.row_open(text), name
    assert h35.ledger_row(bad["ten"]) is None and h35.ledger_row(bad["other_first"]) is not None
    # an outcome in the outcome cell closes the row to the steps
    done = _led(_row(outcome='Run ... "OUTCOME H35: category KEPT; contrast IDLE; whole DROP; portion IDLE."'))
    assert h35.row_formed(done) and not h35.row_open(done)


def test_every_outcome_line_reads_back_from_the_outcome_cell(h35):
    """All 81 word sets: the outcome line the read prints, quoted in the outcome cell as the main session quotes the
    read's fixed words, reads back by code as the same words, the KEPT forms allowed, STOP exactly when none is
    KEPT."""
    for words in itertools.product(("KEPT", "DROP", "IDLE"), repeat=4):
        w = dict(zip(FORMS, words, strict=True))
        verdict = h35.verdict_of(w)
        line = h35.outcome_line(verdict, w)
        assert h35.OUTCOME_RE.fullmatch(line), line
        got = h35.screen_outcome(_led(_row(outcome=f'**Run at x.** In the read\'s fixed words: "{line}" "A screen."')))
        assert got == {"outcome": verdict, "words": w, "allowed": [f for f in FORMS if w[f] == "KEPT"]}, line
        assert (verdict == "STOP") == ("KEPT" not in words) and line.startswith(
            "OUTCOME H35: STOP: " if verdict == "STOP" else "OUTCOME H35: category ")
    void = h35.outcome_line("VOID", why=["v3: no gate record", "v4: x"])
    assert void == "OUTCOME H35: VOID: v3: no gate record; v4: x; nothing is built."
    assert h35.screen_outcome(_led(_row(outcome=void))) == {"outcome": "VOID", "words": {}, "allowed": []}


@pytest.mark.parametrize("cell", [
    "pending", "", "OUTCOME H35:", "OUTCOME H35: category KEPT",
    "OUTCOME H35: category KEPT; contrast KEPT; whole KEPT; portion KEPT",  # no full stop
    "OUTCOME H35: category KEPT; contrast KEPT; portion KEPT; whole KEPT.",  # another order
    "OUTCOME H35: category kept; contrast KEPT; whole KEPT; portion KEPT.",
    "OUTCOME H35: category **KEPT**; contrast IDLE; whole IDLE; portion IDLE.",
    "OUTCOME H35: STOP: category KEPT; contrast IDLE; whole IDLE; portion IDLE.",  # STOP with a KEPT form
    "OUTCOME H35: category IDLE; contrast IDLE; whole DROP; portion IDLE.",  # none KEPT without STOP
    "OUTCOME H35: category KEPT; contrast IDLE; whole IDLE; portion IDLE. OUTCOME H35: VOID: x; nothing is built.",
    "OUTCOME H35: category KEPT; contrast IDLE; whole IDLE; portion IDLE. OUTCOME H35: category KEPT; contrast "
    "IDLE; whole IDLE; portion IDLE.",
    "OUTCOME H34: category KEPT; contrast IDLE; whole IDLE; portion IDLE.",
    "OUTCOME H35: category PASS; contrast IDLE; whole IDLE; portion IDLE.",
    "OUTCOME H35: STOP (category IDLE; contrast IDLE; whole IDLE; portion IDLE)."])
def test_a_malformed_outcome_cell_allows_nothing(h35, cell):
    assert h35.screen_outcome(_led(_row(outcome=cell))) is None


# ------------------------------------------------------------------------------------------------ the arms (synthetic)
def test_every_arm_is_ots125s_practice_form_with_one_change(h35, dform):
    out, extra = h35.build_arms(dform.ctx(), dform.whole, dform.target)
    assert list(out) == ["copy", "control", *FORMS] == list(extra["records"])
    c = extra["checks"]
    for k in ("same_draw", "pre_size_coords_equal", "post_size_coords_equal"):
        assert c[k] is True, k
    assert len(set(c["n_cells"].values())) == 1 and set(c["changed_counts"]) == {"control", *FORMS}
    for st in out.values():
        np.testing.assert_array_equal(st.coords, out["copy"].coords)
    # the control is the banked step at 1.25 sized and halved as D_f; each form its own step, sized and halved alike
    for name, form in (("control", None), *((f, f) for f in FORMS)):
        st = trend_forms.form_step(dform.ctx(), 1.25, form, dform.whole)[0]
        want, _ = h35.h34.size_and_half(dform.ctx(), st, dform.target, name)
        np.testing.assert_array_equal(out[name].X, want.X)
        rec = extra["records"][name]["trend"]
        assert rec["form"] == (form or "banked") and rec["strength"] == 1.25 and rec["budget_days"] == 0.375
        assert rec["ref_cells"] == dform.inputs["ref"].n_obs < dform.whole.n_obs  # the reference is D's pool
        assert "Tiny" in rec["shifted"]  # B0's own factor, the mapping's fallback, in every arm
    for f in FORMS:
        assert not np.array_equal(out[f].X, out["control"].X), f
    assert extra["records"]["copy"]["trend"] == {"form": None, "strength": None}
    # the control's step is the banked builder step: at strength 1 the harness's b0_trend
    ident = h35.identity(dform.ctx(), dform.whole)
    assert ident == {"s1_equals_harness_b0_trend": True, "zeroed_equals_banked": {f: True for f in FORMS}, "ok": True}


def test_the_identities_can_fail(h35, dform, monkeypatch):
    own = trend_forms.form_arguments

    def off(ctx, strength, form=None, whole=None, zero=False, spec=None):
        args = own(ctx, strength, form, whole, zero, spec)
        if form == "portion" and zero:
            args.mapping = {k: v * 2.0 for k, v in args.mapping.items()}
        return args
    monkeypatch.setattr(trend_forms, "form_arguments", off)
    ident = h35.identity(dform.ctx(), dform.whole)
    assert ident["ok"] is False and ident["zeroed_equals_banked"]["portion"] is False
    assert all(ident["zeroed_equals_banked"][f] for f in FORMS if f != "portion")


def test_g0_on_practice_contexts_idle_wild_fall_through_and_void(h35, dform, monkeypatch):
    """Each form's word and condition (what G0 writes: those two and nothing else) and the figures behind them (in
    memory only), on practice contexts."""
    ctx = dform.ctx()
    banked = trend_forms.form_arguments(ctx, 1.25, None)
    applied = trend_forms.applied(ctx, banked)
    assert set(h35.G0_CONDITIONS) == set(h35.G0_CONDITION_WORDS) == {"none", "wild", "not formed", "fall-through",
                                                                      "idle", "code fault"}
    # the forms on these stages move the step: the word follows the figures, which are not in what G0 writes
    for f in FORMS:
        rec, fig = h35.g0_assess(ctx, dform.whole, f, banked, applied)
        assert rec == h35.g0_form(ctx, dform.whole, f, banked, applied) and set(rec) == {"verdict", "condition"}
        want = ("none", "PASS") if fig["changed_share"] >= 0.02 else ("idle", "IDLE")
        assert (rec["condition"], rec["verdict"]) == want and fig["wild_ratio"] < 3, f
        assert {"type": "CM", "status": "shifted"} in fig["statuses"]
    assert {h35.g0_form(ctx, dform.whole, f, banked, applied)["verdict"] for f in ("contrast", "portion")} == {"PASS"}
    # a form whose change is out is idle: changed share 0
    own = trend_forms.form_arguments
    monkeypatch.setattr(trend_forms, "form_arguments",
                        lambda c, s, form=None, whole=None, zero=False, spec=None: own(c, s, form, whole, True, spec))
    rec, fig = h35.g0_assess(ctx, dform.whole, "contrast", banked, applied)
    assert rec == {"verdict": "IDLE", "condition": "idle"} and fig["changed_share"] == 0.0
    monkeypatch.setattr(trend_forms, "form_arguments", own)
    # wild: three times the banked mean |log f| is the edge (not wild), above it wild; idle under 0.02
    stats = trend_forms.g0_statistics
    for ratio, share, want in ((3.0, 0.5, ("PASS", "none")), (3.0 + 1e-9, 0.5, ("DROP", "wild")),
                               (3.0 + 1e-9, 0.0, ("DROP", "wild")), (1.0, 0.02, ("PASS", "none")),
                               (1.0, 0.02 - 1e-12, ("IDLE", "idle"))):
        monkeypatch.setattr(trend_forms, "g0_statistics", lambda f, b, ratio=ratio, share=share: {
            **stats(f, b), "wild_ratio": ratio, "changed_share": share})
        rec = h35.g0_form(ctx, dform.whole, "whole", banked, applied)
        assert (rec["verdict"], rec["condition"]) == want, (ratio, share)
    monkeypatch.setattr(trend_forms, "g0_statistics", stats)
    # fall-through: every Endothelial cell of one category at every stage
    one = {r: harness._shallow(a, celltype=np.where(np.asarray(a.obs["celltype"]) == "BEC", "Endo",
                                                     np.asarray(a.obs["celltype"])))
           for r, a in ctx.inputs.items()}
    from dataclasses import replace
    c1 = replace(ctx, inputs=one)
    b1 = trend_forms.form_arguments(c1, 1.25, None)
    rec, fig = h35.g0_assess(c1, dform.whole, "category", b1, trend_forms.applied(c1, b1))
    assert rec == {"verdict": "IDLE", "condition": "fall-through"} and fig["fall_through"] is True
    assert fig["changed_share"] == pytest.approx(0.0, abs=1e-20)
    # a form that cannot be formed is dropped; a draw that is not OTS125's is a code fault; neither has figures
    assert h35.g0_assess(ctx, None, "whole", banked, applied) == ({"verdict": "DROP", "condition": "not formed"}, {})
    monkeypatch.setattr(trend_forms, "form_arguments", lambda *a, **k: SimpleNamespace(
        **{**own(*a, **k).__dict__, "rows": banked.rows[::-1]}))
    assert h35.g0_assess(ctx, dform.whole, "portion", banked, applied) == (
        {"verdict": "VOID", "condition": "code fault"}, {})
    monkeypatch.setattr(trend_forms, "form_arguments", own)
    # the reading of every form: a code fault is VOID, with fixed words naming the form and no figure
    monkeypatch.setattr(h35, "g0_form", lambda c, w, f, b, a: (
        {"verdict": "VOID", "condition": "code fault"} if f == "whole" else {"verdict": "PASS", "condition": "none"}))
    res = h35.g0_reading(ctx, dform.whole)
    assert res["verdict"] == "VOID" and res["why"] == ["G0: whole: its draw or coordinates are not OTS125's (a code "
                                                       "fault)"]
    assert res["line"] == ("OUTCOME H35: VOID: G0: whole: its draw or coordinates are not OTS125's (a code fault); "
                           "nothing is built.")


def test_the_g0_statistics_by_hand(h35):
    Lb = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 2.0]])  # types: skipped, a, b
    Lf = np.array([[0.0, 0.0], [1.0, 1.0], [0.0, 2.0], [3.0, 0.0]])  # skipped, a, b, a new type c
    banked = {"logF": Lb, "tidx": np.array([0, 1, 1, 2, 2, 2])}
    form = {"logF": Lf, "tidx": np.array([0, 1, 3, 2, 2, 2])}
    st = trend_forms.g0_statistics(form, banked)
    # changed: one cell (1, 1) - (1, 0) = 1; one cell (3, 0) - (1, 0) = 4; banked: 2 x 1 + 3 x 4 = 14
    assert st["changed_share"] == pytest.approx(5 / 14)
    assert st["mean_abs_logf_banked"] == pytest.approx((0 + 0.5 + 0.5 + 1 + 1 + 1) / 6)
    assert st["mean_abs_logf"] == pytest.approx((0 + 1 + 1.5 + 1 + 1 + 1) / 6)
    assert st["wild_ratio"] == pytest.approx(5.5 / 4)
    nothing = {"logF": np.zeros((1, 2)), "tidx": np.zeros(6, dtype=int)}
    assert trend_forms.g0_statistics(nothing, nothing)["changed_share"] == 0.0
    assert trend_forms.g0_statistics(nothing, nothing)["wild_ratio"] == 0.0
    assert trend_forms.g0_statistics(banked, nothing)["changed_share"] == float("inf")
    assert trend_forms.g0_statistics(banked, nothing)["wild_ratio"] == float("inf")
    with pytest.raises(ValueError, match="drew different cells"):
        trend_forms.g0_statistics({"logF": Lf, "tidx": np.zeros(2, dtype=int)}, banked)


# ------------------------------------------------------------------------------------------------ G0 on the test inputs
@pytest.fixture
def g0_env(h35, monkeypatch, own_repo, tmp_path):
    """The heart builder on the synthetic test inputs and test board of tests/test_heart_trend_line_h.py (every stage
    relabelled, E10.5 with a detached island of vessels), the ledger there with row R55, a temporary home for G0's
    record: the builder is loaded as G0 loads it."""
    from test_heart_trend_build import install
    from test_heart_trend_line_d import _organ_size
    from test_heart_trend_line_h import LEDGER_TF, STAGES, relabel, tf_organ_stage

    bt = h35.builder()
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
    led = tmp_path / "TRIES.md"
    led.write_text(LEDGER_TF)
    monkeypatch.setattr(bt, "TRIES_FILE", led)
    monkeypatch.setattr(h35, "refuse_state", lambda commit, what="run": None)
    monkeypatch.setattr(h35, "git_state", lambda: {"commit": COMMIT, "tree_clean": True, "utc": "x"})
    return bt, reads


G0_KEYS = {"row", "axis", "commit", "provenance", "label", "verdict", "forms", "why", "line", "follows", "source"}


def _names_out(text: str) -> str:
    """The text without the names that hold digits (the commit, rows H35 and R55, the files, the step G0): what is
    left of a G0 output must hold no numeral."""
    for name in (COMMIT, "OTK125", "OTX125", "OTW125", "OTP125", "OTS125", "H35", "R55", "G0"):
        text = text.replace(name, "")
    return text


def _no_number(x, path="record"):
    """No number anywhere: strings, booleans and None only (a figure over the drawn cells gives back their number)."""
    if isinstance(x, dict):
        for k, v in x.items():
            _no_number(v, f"{path}.{k}")
    elif isinstance(x, list | tuple):
        for i, v in enumerate(x):
            _no_number(v, f"{path}[{i}]")
    else:
        assert isinstance(x, (str, bool)) or x is None, (path, type(x), x)


def test_g0_takes_the_builders_own_context_and_writes_each_forms_word_and_condition_only(h35, g0_env, tmp_path,
                                                                                         capsys):
    bt, reads = g0_env
    ctx, whole = h35.builder_context(bt=bt)
    assert reads == ["E10.5.h5ad", "E9.5.h5ad", "E8.75.h5ad", "E8.5.h5ad", "E8.25_late.h5ad"]
    assert ctx.seed == harness.BUILD_SEED_OFFSET + 1 and ctx.task.key == "test:T2h" and ctx.board.max_cells == 733
    assert ctx.inputs["ref"].n_obs < whole.n_obs and set(ctx.inputs) == {"ref", "start", "e8_75", "e8_5",
                                                                          "e8_25_late"}
    assert bt._trend_form_step.__name__ == "_trend_form_step"  # the builder is as it was
    # the builder's step on that context is OTS125's with the form, bit for bit
    st, _ = bt._trend_form_step(ctx, 1.25, "whole", whole)
    np.testing.assert_array_equal(st.X, trend_forms.form_step(ctx, 1.25, "whole", whole)[0].X)
    out = tmp_path / "g0"
    code = h35.run_g0(out, None, COMMIT, bt=bt)
    printed = capsys.readouterr().out.splitlines()
    g = json.loads((out / "g0.json").read_text())
    assert code == 0 and g["verdict"] == "PASS" and g["axis"] == h35.G0_AXIS and g["commit"] == COMMIT
    # the record: the run's envelope and, per form, its word and the condition that set it; no figure
    assert set(g) == G0_KEYS and g["why"] == [] and g["line"] is None and g["follows"] is None
    assert g["provenance"] == {"commit": COMMIT, "tree_clean": True, "utc": "x"}
    for f in FORMS:
        r = g["forms"][f]
        assert set(r) == {"verdict", "condition"} and h35.G0_CONDITIONS[r["condition"]] == r["verdict"], f
    _no_number(g)
    assert printed == [f"G0 {f} {g['forms'][f]['verdict']}" for f in FORMS] + ["G0 PASS"]
    assert h35.passed_forms(g) == [f for f in FORMS if g["forms"][f]["verdict"] == "PASS"] and h35.passed_forms(g)
    # the text: its fixed head, one line per form (word, condition and the condition's fixed words), the verdict
    assert (out / "g0.txt").read_text().splitlines() == [
        h35.G0_TEXT_HEAD.format(row="H35", commit=COMMIT, ruling="R55"),
        *(f"G0 {f} ({h35.FILES[f]}): {g['forms'][f]['verdict']} ({g['forms'][f]['condition']}: "
          f"{h35.G0_CONDITION_WORDS[g['forms'][f]['condition']]})." for f in FORMS), "G0 PASS."]
    for words in (h35.G0_TEXT_HEAD, *h35.G0_CONDITION_WORDS.values()):
        assert not re.search(r"\d", _names_out(words)), words  # no numeral in G0's fixed words but in names
    with pytest.raises(SystemExit, match="already holds a G0 record"):
        h35.run_g0(out, None, COMMIT, bt=bt)


def test_no_number_from_the_test_draw_reaches_any_g0_output(h35, g0_env, monkeypatch, tmp_path, capsys):
    """The review of aeb070e: a figure over the drawn cells gives back their number, the hidden count. With the test
    board's limit set to two different counts (two draws of two sizes from the same pool), G0's two files and its
    printed lines are the same bytes, while the figures it keeps in memory differ; with those figures written, as
    before the fix, the files differ. So nothing G0 writes or prints is a function of the draw beyond the words."""
    from dataclasses import replace as dc_replace

    from test_heart_trend_build import heart_test_board

    bt, _ = g0_env

    def at(n: int):
        board = dc_replace(heart_test_board(), max_cells=n)
        monkeypatch.setattr(bt, "load_boards", lambda: {board.key: board})

    def run(n: int, name: str):
        at(n)
        assert h35.run_g0(tmp_path / name, None, COMMIT, bt=bt) == 0
        printed = capsys.readouterr().out
        return {f: (tmp_path / name / f).read_bytes() for f in ("g0.json", "g0.txt")}, printed

    def figures(n: int):
        at(n)
        ctx, whole = h35.builder_context(bt=bt)
        banked = trend_forms.form_arguments(ctx, 1.25, None)
        applied = trend_forms.applied(ctx, banked)
        assert len(banked.rows) == n  # the draw's size is the board's limit
        return {f: h35.g0_assess(ctx, whole, f, banked, applied)[1] for f in FORMS}

    a, printed_a = run(733, "a")
    b, printed_b = run(641, "b")
    assert a == b and printed_a == printed_b  # every byte G0 writes or prints, at either count
    fa, fb = figures(733), figures(641)
    assert all(fa[f]["changed_share"] != fb[f]["changed_share"] for f in FORMS)
    assert all(fa[f]["mean_abs_logf_banked"] != fb[f]["mean_abs_logf_banked"] for f in FORMS)
    for files in (a, b):
        _no_number(json.loads(files["g0.json"]))
        assert not re.search(r"\d", _names_out(files["g0.txt"].decode()))  # no numeral but in names
    # the same check sees a figure written: G0 as it was before the fix (its figures beside each word)
    own = h35.g0_form
    monkeypatch.setattr(h35, "g0_form", lambda c, w, f, bk, ap: {
        **own(c, w, f, bk, ap), "changed_share": h35.g0_assess(c, w, f, bk, ap)[1].get("changed_share")})
    c, _ = run(733, "c")
    d, _ = run(641, "d")
    assert c["g0.json"] != d["g0.json"]


def test_g0_stops_when_no_form_passes_and_says_so(h35, g0_env, monkeypatch, tmp_path, capsys):
    bt, _ = g0_env
    stats = trend_forms.g0_statistics
    monkeypatch.setattr(trend_forms, "g0_statistics", lambda f, b: {**stats(f, b), "changed_share": 0.01})
    code = h35.run_g0(tmp_path / "g0", None, COMMIT, bt=bt)
    printed = capsys.readouterr().out.splitlines()
    assert code == 3 and printed[-3] == "G0 STOP"
    assert printed[-2] == "OUTCOME H35: STOP: category IDLE; contrast IDLE; whole IDLE; portion IDLE."
    assert printed[-1] == h35.follows_line("STOP") and "nothing is built under row R55" in printed[-1]
    assert h35.screen_outcome(_led(_row(outcome=printed[-2])))["allowed"] == []
    g = json.loads((tmp_path / "g0" / "g0.json").read_text())
    assert h35.g0_record_reasons(g, COMMIT) == ["v3: G0 said 'STOP' and passed no form"]


# ------------------------------------------------------------------------------------------------ hand-made tables
POOL_CELLS, DRAWN = 900, 700
SIZE_REC = {"sizing": "f", "target_rms_full": 3.0, "base_rms": 2.0, "ratio": 1.5, "log_rms_after_vs_target": 0.0}
DRAW = {"pool_cells": POOL_CELLS, "drawn": DRAWN, "short": False, "repeated_rows_at_least": 0}
G1_BASE = {"core_iso_half:f": 49.0, "half_iso_half:f": 49.4, "full_iso_half:f": 49.9}


def _off(v, q, s) -> float:
    if isinstance(v, dict):
        v = v.get("own" if q == "own" else "CLEAN", 0.0)
    return float(np.broadcast_to(np.asarray(v, dtype=float), (10,))[s])


def _base(h35, key, s, q) -> float:
    return 50.0 + 0.5 * h35.TASK_KEYS.index(key) + 0.01 * s + 0.001 * h35.TRUTHS.index(q)


def _trec(h35, key, name):
    if name == "copy":
        return {"form": None, "strength": None}
    return {"form": "banked" if name == "control" else name, "strength": 1.25,
            "budget_days": h35.PRACTICE_BUDGET[key], "ref_cells": POOL_CELLS, "shifted": ["CM", "Endothelial"]}


def _prow(h35, key, s, forms, score=None):
    """score: arm -> board points added to the control's score (the copy sits 0.4 below the control)."""
    score = score or {}
    built = ("copy", "control", *forms)
    sc, sk, gr = {}, {}, {}
    for q in h35.TRUTHS:
        base = _base(h35, key, s, q)
        sc[q] = {"floor": 50.0, "copy": base - 0.4, **{a: base + _off(score.get(a, 0.0), q, s)
                                                       for a in ("control", *forms)}}
        sk[q] = {a: {m: 0.5 for m in METRICS} for a in ("floor", *built)}
        gr[q] = {a: {g: 0.5 for g in GROUPS} for a in ("floor", *built)}
    checks = {"same_draw": True, "pre_size_coords_equal": True, "post_size_coords_equal": True,
              "changed_counts": {a: 9 for a in ("control", *forms)}, "n_cells": {a: 350 for a in built}}
    return {"seed": s, "build_seed": harness.BUILD_SEED_OFFSET + s, "n_cells": {a: 350 for a in built},
            "checks": checks, "records": {a: {"size": dict(SIZE_REC), "trend": _trec(h35, key, a)} for a in built},
            "draw": dict(DRAW), "scores": sc, "skills": sk, "groups": gr}


def _facts(h35):
    organ = {"rule": "x", "label_config": "c.yaml", "label_config_sha256": "0" * 64,
             "kept_categories": ["cm", "endocardium"], "link_percentile": 99.0, "link_length": 2.5, "kept_share": 0.5,
             "chamber_share": 0.6, "core_share_of_kept": 0.7, "components_kept": 2, "components_total": 5}
    return {"truth_info": {q: {"cells": 1000 + i, "core_share": 0.4, "xy_scope": 1.0}
                           for i, q in enumerate(h35.TRUTHS)},
            "pred_info": {"organ": {"pool_cells": POOL_CELLS, "core_share": 0.7, "xy_scope": 0.1}},
            "organ_record": {"reference": organ}, "ratios": {"organ": 1.5, "label": 1.8},
            "d_targets": {"organ:f": {"target": 3.0, "base": 2.0, "ratio": 1.5, "source": "x"}},
            "forecast_targets": {"core_iso_half": 2.0}, "v7": {"ok": True, "log_diff": 0.0}}


def _part(h35, key, seeds, forms=FORMS, commit=COMMIT, clean=True, **kw):
    return {"axis": h35.AXIS, "task": key, "forms": list(forms), **_facts(h35),
            "rows": [_prow(h35, key, s, forms, **kw) for s in seeds],
            "provenance": {"commit": commit, "tree_clean": clean, "utc": "x"}}


def _parts(h35, forms=FORMS, **kw):
    return {k: [_part(h35, k, range(0, 5), forms, **kw), _part(h35, k, range(5, 10), forms, **kw)]
            for k in h35.TASK_KEYS}


def _g0(h35, verdicts=None, verdict=None, commit=COMMIT, clean=True):
    verdicts = verdicts or {f: "PASS" for f in FORMS}
    v = verdict or ("PASS" if "PASS" in verdicts.values() else "STOP")
    return {"row": "H35", "axis": h35.G0_AXIS, "verdict": v, "commit": commit,
            "forms": {f: {"verdict": w, "why": []} for f, w in verdicts.items()}, "why": [],
            "provenance": {"commit": commit, "tree_clean": clean}}


def _gate_rows(h35, seeds=range(10)):
    rows = {}
    for k in h35.TASK_KEYS:
        rows[k] = []
        for s in seeds:
            sc = _prow(h35, k, s, ())["scores"]
            rows[k].append({"seed": s, "scores": {q: {a: sc[q][a] for a in ("floor", "copy")} for q in h35.TRUTHS},
                            "size_record": dict(SIZE_REC), "draw": dict(DRAW),
                            "identity": {"ok": True, "s1_equals_harness_b0_trend": True,
                                         "zeroed_equals_banked": {f: True for f in FORMS}}})
    return rows


def _gate(h35, verdict="PASS", commit=COMMIT, clean=True, seeds=range(10)):
    return {"row": "H35", "axis": h35.GATE_AXIS, "verdict": verdict, "commit": commit,
            "why": [] if verdict == "PASS" else ["G1 (M): x"], "found": {}, "rows": _gate_rows(h35, seeds),
            "transitions": {k: _facts(h35) for k in h35.TASK_KEYS}, "provenance": {"commit": commit,
                                                                                  "tree_clean": clean}}


def _read(h35, parts=None, g0="default", gate="default", commit=COMMIT, unreadable=(), **kw):
    return h35.reading(_parts(h35, **kw) if parts is None else parts, _g0(h35) if g0 == "default" else g0,
                       _gate(h35) if gate == "default" else gate, commit, unreadable)


KEPT_ALL = {f: {"own": 0.2, "CLEAN": 0.2} for f in FORMS}


def test_the_words_kept_drop_idle_and_the_outcome(h35):
    res = _read(h35, score=KEPT_ALL)
    assert res["void"] is False and res["outcome"] == "KEPT" and res["words"] == {f: "KEPT" for f in FORMS}
    assert res["line"] == "OUTCOME H35: category KEPT; contrast KEPT; whole KEPT; portion KEPT."
    assert res["follows"] == h35.follows_line("KEPT") and "built twice" in res["follows"]
    score = {"category": -0.5, "contrast": {"own": 0.05, "CLEAN": -0.08}, "whole": {"own": 0.4, "CLEAN": 0.0},
             "portion": {"own": 0.0, "CLEAN": -0.35}}
    res = _read(h35, score=score)
    assert res["words"] == {"category": "DROP", "contrast": "IDLE", "whole": "KEPT", "portion": "DROP"}
    assert res["line"] == "OUTCOME H35: category DROP; contrast IDLE; whole KEPT; portion DROP."
    res = _read(h35, score={"category": -0.5, "contrast": 0.0, "whole": 0.05, "portion": -0.05})
    assert res["outcome"] == "STOP" and res["line"] == (
        "OUTCOME H35: STOP: category DROP; contrast IDLE; whole IDLE; portion IDLE.")
    assert res["follows"] == h35.follows_line("STOP") and "nothing is built" in res["follows"]


def _word(h35, own, clean):
    tab = {w: {h35.FOLD: {h35.pname("whole", "control"): {"score": s}}} for w, s in (("own", own), ("CLEAN", clean))}
    return h35.form_word(tab, "whole")["word"]


def test_the_bars_at_their_edges(h35):
    ok = {"mean": 0.2, "behind": 0}
    assert _word(h35, {"mean": -0.3000001, "behind": 8}, ok) == "DROP"
    assert _word(h35, {"mean": -0.3, "behind": 10}, ok) == "KEPT"  # a mean of exactly -0.30 is not below it
    assert _word(h35, {"mean": -0.9, "behind": 7}, ok) == "KEPT"  # seven seeds below: not DROP
    assert _word(h35, ok, {"mean": -0.31, "behind": 8}) == "DROP"  # either world
    assert _word(h35, {"mean": 0.1, "behind": 0}, {"mean": -0.1, "behind": 9}) == "IDLE"  # within, in both
    assert _word(h35, {"mean": 0.1000001, "behind": 0}, {"mean": 0.0, "behind": 0}) == "KEPT"
    assert _word(h35, {"mean": 0.05, "behind": 0}, {"mean": -0.35, "behind": 8}) == "DROP"  # DROP before IDLE
    assert _word(h35, {"mean": 0.0, "behind": 0}, {"mean": 0.0, "behind": 0}) == "IDLE"  # no change reads IDLE

    def word(score):  # through the tables: a difference of exactly 0 in a seed is neither above nor below
        return _read(h35, score={"category": score})["words"]["category"]
    assert word(np.array([0.0] * 2 + [-0.6] * 8)) == "DROP"
    assert word(np.array([0.0] * 3 + [-0.6] * 7)) == "KEPT"
    assert word(np.array([-0.4] * 8 + [0.2, 0.2])) == "KEPT"  # mean -0.28 with eight below: not DROP


def test_g0_words_are_kept_for_the_forms_g0_did_not_pass(h35):
    g0 = _g0(h35, {"category": "IDLE", "contrast": "PASS", "whole": "DROP", "portion": "PASS"})
    res = _read(h35, g0=g0, forms=("contrast", "portion"), score={"contrast": 0.5, "portion": 0.0})
    assert res["void"] is False and res["forms"] == ["contrast", "portion"]
    assert res["words"] == {"category": "IDLE", "contrast": "KEPT", "whole": "DROP", "portion": "IDLE"}
    assert res["g0_words"] == {"category": "IDLE", "whole": "DROP"}
    assert "category (OTK125): IDLE at G0, not read." in h35.text(res)
    # the forms the parts hold must be G0's passing forms
    res = _read(h35, g0=g0, forms=FORMS)
    assert res["void"] is True and any("forms are not the ones G0 passed" in w for w in res["why"])


def test_every_void_condition(h35):
    def why(**kw):
        return _read(h35, **kw)["why"]
    assert why() == []
    assert why(g0=None) == ["v3: no G0 record"]
    assert why(g0=_g0(h35, commit="d" * 40)) == ["v3: the G0 record is not of the run commit"]
    assert why(g0=_g0(h35, {f: "IDLE" for f in FORMS})) == ["v3: G0 said 'STOP' and passed no form"]
    assert why(gate=None) == ["v3: no gate record"]
    assert why(gate=_gate(h35, "FAIL")) == ["v3: the gate said 'FAIL', not PASS"]
    assert why(gate=_gate(h35, clean=False)) == ["v3: the gate record is not of the run commit"]
    p = _parts(h35)
    del p[KEY2]
    assert why(parts=p) == [f"v4: missing transition {KEY2}"]
    p = _parts(h35)
    p[KEY1][1]["rows"] = p[KEY1][1]["rows"][:-1]
    assert why(parts=p) == [f"v4: {KEY1}: seeds {list(range(9))}, not {list(range(10))}"]
    p = _parts(h35)
    p[KEY1][1]["rows"][0]["seed"] = 4
    assert any("occur in more than one part" in w for w in why(parts=p))
    p = _parts(h35)
    del p[KEY1][0]["rows"][2]["scores"]["thin_0"]["whole"]
    assert why(parts=p) == [f"v4: {KEY1} seed 2: missing from scores: ['thin_0/whole']"]
    p = _parts(h35)
    p[KEY2][0]["rows"][1]["scores"]["own"]["portion"] = float("nan")
    assert why(parts=p) == [f"v4: {KEY2} seed 1: a score is not a finite number: ['own/portion']"]
    assert why(parts=_parts(h35, clean=False)) == ["v4: a part comes from an unclean tree"]
    assert why(parts=_parts(h35, commit="d" * 40)) == [
        f"v4: a part comes from a commit other than the run commit {COMMIT!r}"]
    assert why(unreadable=["x/h35.json"]) == ["v4: a part cannot be read: x/h35.json"]
    p = _parts(h35)
    p[KEY1][0]["rows"][0]["checks"]["same_draw"] = False
    assert why(parts=p) == [f"v5: {KEY1} seed 0: same draw fails"]
    p = _parts(h35)
    p[KEY1][0]["rows"][0]["skills"]["radial_0"]["whole"]["occupancy_dice"] += 2e-9
    assert why(parts=p) == [f"v5: {KEY1} seed 0 radial_0: whole's occupancy_dice is not the control's"]
    p = _parts(h35)
    p[KEY1][0]["rows"][0]["skills"]["own"]["contrast"]["d2_shape"] += 5e-10  # within the tolerance
    assert why(parts=p) == []
    p = _parts(h35)
    p[KEY2][0]["rows"][3]["records"]["whole"]["trend"]["strength"] = 1.0
    assert why(parts=p) == [f"v5: {KEY2} seed 3: whole's record does not hold form whole at strength 1.25"]
    p = _parts(h35)
    p[KEY2][0]["rows"][3]["records"]["control"]["trend"]["ref_cells"] = POOL_CELLS + 1
    assert why(parts=p) == [f"v5: {KEY2} seed 3: control's reference is not D's pool"]
    p = _parts(h35)
    p[KEY1][1]["rows"][0]["scores"]["own"]["copy"] += 1e-9
    assert why(parts=p) == [f"v5: {KEY1} seed 5 own: the run step's copy score is not the gate's"]
    res = _read(h35, gate=None)
    assert res["outcome"] == "VOID" and res["line"] == "OUTCOME H35: VOID: v3: no gate record; nothing is built."
    assert h35.screen_outcome(_led(_row(outcome=res["line"])))["outcome"] == "VOID"


# ------------------------------------------------------------------------------------------------ G1 and the gate
def _h27_rows(h35, seeds=range(10)):
    return {k: [{"seed": s, "scores": {q: {a: G1_BASE[a] + 0.3 * h35.TASK_KEYS.index(k) + 0.02 * ((7 * s) % 10)
                                           + 0.003 * h35.TRUTHS.index(q) for a in h35.h34.G1_ARMS}
                                       for q in h35.TRUTHS}} for s in seeds] for k in h35.TASK_KEYS}


def _committed(h35, seeds=range(10)):
    h29b = h35.h29b
    rows, others = _gate_rows(h35, seeds), _h27_rows(h35, seeds)
    f = _facts(h35)
    tables = {}
    for w, ts in h35.h34.G1_WORLDS:
        tables[w] = {}
        for u in h35.h34.G1_UNITS:
            c = h29b.series(rows, h29b.UNITS[u], ts, "copy")
            tables[w][u] = {"mean": {"organ:f": float(c.mean())},
                            "stated": {f"organ:f - {a}": h29b._stats(c - h29b.series(others, h29b.UNITS[u], ts, a))
                                       for a in h35.h34.G1_ARMS}}

    def per(v):
        return {k: copy.deepcopy(v) for k in h35.TASK_KEYS}
    short = {"organ": {"pool_cells": POOL_CELLS, "drawn": [DRAWN], "short_in_seeds": 0, "repeated_rows_at_least": 0}}
    return {"row": "H29b", "void": False, "outcome": "READ", "tables": tables, "truth_info": per(f["truth_info"]),
            "pool_info": per(f["pred_info"]), "organ_records": per(f["organ_record"]), "ratios": per(f["ratios"]),
            "d_targets": per(f["d_targets"]), "pools_short": per(short)}


def _gate_reasons(h35, rows=None, committed="default", h27="default"):
    return h35.gate_reasons(_gate_rows(h35) if rows is None else rows, {k: _facts(h35) for k in h35.TASK_KEYS},
                            _committed(h35) if committed == "default" else committed,
                            _h27_rows(h35) if h27 == "default" else h27)


def test_the_gate_passes_only_if_g1_and_every_identity_hold(h35):
    why, found = _gate_reasons(h35)
    assert why == [] and found["M"] == {"means": 18, "unequal": [], "largest_abs_diff": 0.0}
    assert found["P"] == {"numbers": 270, "unequal": []} and found["I"] == {"contexts": 20, "failed": 0}
    rows = _gate_rows(h35)
    rows[KEY2][7]["identity"]["ok"] = False
    assert _gate_reasons(h35, rows=rows)[0] == [f"I: 1 of 20 contexts fail an identity (first: {KEY2} seed 7)"]
    c = _committed(h35)
    c["tables"]["CLEAN"]["fold E10.5"]["mean"]["organ:f"] *= 1 + 1e-12
    assert _gate_reasons(h35, committed=c)[0] == [
        "G1 (M): 1 of 18 means are not row H29b's committed D_f means (first: CLEAN fold E10.5)"]
    c = _committed(h35)
    c["tables"]["own"][KEY1]["stated"]["organ:f - full_iso_half:f"]["se"] *= 1 + 1e-12
    assert _gate_reasons(h35, committed=c)[0] == [
        f"G1 (P): 1 of 270 paired statistics are not row H29b's committed (first: own {KEY1} organ:f - "
        "full_iso_half:f se)"]
    assert _gate_reasons(h35, committed=None)[0] == [
        "G1: row H29b's reading.json is not at hand as the file its row records (missing, another sha256, or void)"]
    assert _gate_reasons(h35, h27=None)[0] == ["G1 (P): row H27's committed parts are not at hand"]
    rows = _gate_rows(h35, range(9))
    assert "G1: the copy D_f was not built and scored on both transitions and every seed" in _gate_reasons(
        h35, rows=rows)[0]


def test_g1_is_row_h34s_g1_with_the_copy_in_its_controls_place(h35, monkeypatch):
    """On the same hand-made rows (the copy's scores under row H34's control name), row H34's own gate_reasons finds
    the same G1 failures, word for word but the arm's name."""
    h34 = h35.h34
    for change in (None, "M", "P"):
        c = _committed(h35)
        if change == "M":
            c["tables"]["own"][KEY2]["mean"]["organ:f"] += 1e-9
        if change == "P":
            c["tables"]["CLEAN"][KEY1]["stated"]["organ:f - core_iso_half:f"]["mean"] += 1e-9
        mine = h35.g1_reasons(_gate_rows(h35), {k: _facts(h35) for k in h35.TASK_KEYS}, c, _h27_rows(h35))
        theirs_rows = {k: [{**r, "scores": {q: {"floor": v["floor"], "control": v["copy"]}
                                            for q, v in r["scores"].items()},
                            "factors": {"ok": True, "lineages": []}} for r in rs]
                       for k, rs in _gate_rows(h35).items()}
        theirs = h34.gate_reasons(theirs_rows, {k: _facts(h35) for k in h35.TASK_KEYS}, c, _h27_rows(h35))
        assert mine[0] == [w.replace("the control's draw", "the copy's draw") for w in theirs[0]], change
        assert {k: v for k, v in mine[1].items()} == {k: v for k, v in theirs[1].items() if k != "G2"}, change


def test_g1_dry_runs_on_committed_results(h35):
    """Rule 9 J7: G1's arithmetic (row H34's means and paired, which g1_reasons calls) on committed results other than
    D_f's: row H27's committed parts of A_f, B_f and C_f give the means row H29b committed for those three controls
    exactly, and their paired differences row H27 committed (54 and 216 numbers). A copy with one seed's scores moved
    by one part in 1e12 fails both. Committed files of rows H27 and H29b only, read for arithmetic; nothing is built or
    scored, and no file of row H34 is read."""
    h34 = h35.h34
    committed = h34.load_h29b_reading()
    assert committed is not None and h34.committed_ok(committed)
    rows = h34.load_h27_rows()
    h27_reading = h35.h29b.load_h27_reading()
    letter = h34.h27.LETTER

    def m_off(parts):
        return [(a, w, u) for a in h34.G1_ARMS for w, _ in h34.G1_WORLDS for u in h34.G1_UNITS
                if h34.means(parts, a)[w][u] != committed["tables"][w][u]["mean"][a]]

    def p_off(parts):
        off = []
        for x, y in h34.h26.PAIRS:
            got = h34.paired(parts, f"{x}:f", parts, f"{y}:f")
            for w, _ in h34.G1_WORLDS:
                for u in h34.G1_UNITS:
                    want = h27_reading["tables"][w][u]["diff"]["f"][f"{letter[x]}-{letter[y]}"]
                    off += [(x, y, w, u, k) for k in ("mean", "lo", "hi", "ahead") if got[w][u][k] != want[k]]
        return off
    assert m_off(rows) == [] and p_off(rows) == []
    shifted = json.loads(json.dumps(rows))
    for q in h35.TRUTHS:
        for a in h34.G1_ARMS:
            shifted[KEY1][4]["scores"][q][a] *= 1 + 1e-12
    assert m_off(shifted) and p_off(shifted)


def test_g1_goes_through_row_h34s_means_and_paired(h35, monkeypatch):
    seen = []

    def recorded(name):
        own = getattr(h35.h34, name)

        def call(*a):
            seen.append((name, a[1], *a[3:]))
            return own(*a)
        return call
    for name in ("means", "paired"):
        monkeypatch.setattr(h35.h34, name, recorded(name))
    assert _gate_reasons(h35)[0] == []
    assert seen == [("means", "copy"), *(("paired", "copy", other) for other in h35.h34.G1_ARMS)]


# ------------------------------------------------------------------------------------------------ end to end
@pytest.fixture()
def fake_io(h35, synth, monkeypatch, tmp_path):
    """Every loader of rows H29b's, H34's and this row's scripts stood in by the synthetic stages, and the scorer by a
    recorder: nothing is read from disk and nothing is scored."""
    mod, stages = synth
    h34, h29b, h27 = h35.h34, h35.h29b, h35.h27

    def stage_of(path, *a, **k):
        return stages["T2/" + Path(path).name].copy()

    def task(key):
        return mod._heart_task(key=key, inputs=ROLES[key], reference=ROLES[key]["ref"])
    log = []

    def score(board, pred, truth, ref, seed=0, anchors="computed", build=None, protocol=None):
        log.append(1)
        return SimpleNamespace(score=50.0 + 1e-3 * len(log), skills={m: 0.5 for m in METRICS},
                               group_skill={g: 0.5 for g in GROUPS}, anchors={"of": [board.truth_cells, seed]})
    monkeypatch.setattr(h29b.payoff, "heart_task", task)
    monkeypatch.setattr(h27, "input_geometries", lambda t, data_dir: (_geoms(h34, stages, h27.input_files(t)), []))
    monkeypatch.setattr(prep, "add_lineage", lambda a, source=None: None)
    monkeypatch.setattr(harness, "task_board", lambda t, boards, n_target, n_ref: mod.BOARD)
    for m in (h34, h29b, h35):
        monkeypatch.setattr(m, "load_boards", lambda: {"T2:heart:val_extrap": mod.BOARD})
        monkeypatch.setattr(m, "load_stage", stage_of)
        if hasattr(m, "read_stage"):
            monkeypatch.setattr(m, "read_stage", stage_of)
        if hasattr(m, "score_board"):
            monkeypatch.setattr(m, "score_board", score)
    monkeypatch.setattr(h35, "refuse_state", lambda commit, what="run": None)
    return log


def _stub_scores(h35, offsets):
    def score_arms(tr, s, preds):
        sc, sk, gr = {}, {}, {}
        for q in h35.TRUTHS:
            base = _base(h35, tr.key, s, q)
            sc[q] = {"floor": 50.0, **{a: (base - 0.4 if a == "copy" else base + _off(offsets.get(a, 0.0), q, s))
                                       for a in preds}}
            sk[q] = {a: {m: 0.5 for m in METRICS} for a in ("floor", *preds)}
            gr[q] = {a: {g: 0.5 for g in GROUPS} for a in ("floor", *preds)}
        return sc, sk, gr
    return score_arms


def test_gate_run_and_read_end_to_end_on_stand_in_stages(h35, fake_io, monkeypatch, tmp_path, capsys):
    """The steps through the command line on stand-in stages, two seeds in place of ten (the scores come from a
    stand-in scorer; the builds, the factors, the records and every check are the script's own). G0 is stood in by a
    record that passes three forms."""
    seeds = (0, 1)
    monkeypatch.setattr(h35, "SEEDS", seeds)
    monkeypatch.setattr(h35.h34, "score_arms", _stub_scores(h35, {f: 0.5 for f in FORMS}))
    monkeypatch.setattr(h35, "git_state", lambda: {"commit": COMMIT, "tree_clean": True, "utc": "x"})
    h27 = tmp_path / "h27"
    for k, rows in _h27_rows(h35, seeds).items():
        (h27 / h35.h26.safe(k)).mkdir(parents=True)
        (h27 / h35.h26.safe(k) / "h27.json").write_text(json.dumps({"axis": h35.h27.AXIS, "task": k, "rows": rows}))
    committed = _committed(h35, seeds)
    for key in h35.TASK_KEYS:
        tr = h35.load_transition(key, tmp_path)
        d = h35.h34.build_control(h35.context(tr, 0), tr.target)[2]
        committed["pool_info"][key], committed["organ_records"][key] = tr.facts["pred_info"], tr.facts["organ_record"]
        committed["ratios"][key], committed["d_targets"][key] = tr.facts["ratios"], tr.facts["d_targets"]
        committed["pools_short"][key] = {"organ": {"pool_cells": d["pool_cells"], "drawn": [d["drawn"]],
                                                   "short_in_seeds": len(seeds) * d["short"],
                                                   "repeated_rows_at_least": d["repeated_rows_at_least"]}}
    reading = tmp_path / "h29b_reading.json"
    reading.write_text(json.dumps(committed, default=float))
    monkeypatch.setattr(h35.h34, "H29B_READING_SHA256", hashlib.sha256(reading.read_bytes()).hexdigest())
    out = tmp_path / "results"
    g0_dir, gate_dir, parts_dir = out / "g0", out / "gate", out / "parts"
    g0_dir.mkdir(parents=True)
    (g0_dir / "g0.json").write_text(json.dumps(_g0(h35, {"category": "PASS", "contrast": "PASS", "whole": "IDLE",
                                                         "portion": "PASS"})))
    common = ["--commit", COMMIT, "--data", str(tmp_path)]
    run = ["run", "--task", KEY1, "--seeds", "0-1", "--out", str(parts_dir / "e95" / "0-1"), "--g0", str(g0_dir),
           "--gate", str(gate_dir), *common]
    with pytest.raises(SystemExit, match="has not said PASS"):
        h35.main(run)
    with pytest.raises(SystemExit) as e:
        h35.main(["gate", "--out", str(gate_dir), "--g0", str(g0_dir), "--h29b", str(reading), "--h27", str(h27),
                  *common])
    printed = capsys.readouterr().out.splitlines()
    g = json.loads((gate_dir / "gate.json").read_text())
    assert e.value.code == 0 and printed[-1] == "GATE PASS" and g["why"] == [], g["why"]
    assert g["found"]["I"] == {"contexts": 4, "failed": 0}
    assert all(r["identity"]["ok"] for k in h35.TASK_KEYS for r in g["rows"][k])
    assert set(g["rows"][KEY1][0]["scores"]["own"]) == {"floor", "copy"}  # only the copy is built and scored
    h35.main(run)
    h35.main(["run", "--task", KEY2, "--seeds", "0-1", "--out", str(parts_dir / "e875" / "0-1"), "--g0", str(g0_dir),
              "--gate", str(gate_dir), *common])
    printed = capsys.readouterr().out.splitlines()
    assert [bool(re.fullmatch(r"T2h/[^ ]+ seed [01] \([0-9.]+ s\) done", ln)) for ln in printed] == \
        [True, True, False] * 2
    parts, unreadable = h35.load_parts(parts_dir)
    assert unreadable == [] and parts[KEY1][0]["forms"] == ["category", "contrast", "portion"]
    assert set(parts[KEY1][0]["rows"][0]["scores"]["own"]) == {"floor", "copy", "control", "category", "contrast",
                                                                "portion"}
    h35.main(["read", "--runs", str(parts_dir), "--g0", str(g0_dir), "--gate", str(gate_dir), "--out", str(out),
              "--commit", COMMIT])
    res = json.loads((out / "reading.json").read_text())
    assert res["void"] is False, res["why"]
    # two seeds cannot reach DROP's eight; +0.5 in both worlds: KEPT; whole keeps G0's IDLE
    assert res["words"] == {"category": "KEPT", "contrast": "KEPT", "whole": "IDLE", "portion": "KEPT"}
    line = "OUTCOME H35: category KEPT; contrast KEPT; whole IDLE; portion KEPT."
    assert capsys.readouterr().out == line + "\n" + h35.follows_line("KEPT") + "\n"
    assert (out / "reading.txt").read_text().startswith(line + "\n")
    assert res["records"][KEY1]["category"]["trend"]["split_types"] == ["Endothelial/endocardium",
                                                                       "Endothelial/vascular_endothelium"]
    with pytest.raises(SystemExit, match="already holds a reading"):
        h35.main(["read", "--runs", str(parts_dir), "--g0", str(g0_dir), "--gate", str(gate_dir), "--out", str(out),
                  "--commit", COMMIT])


# ------------------------------------------------------------------------------------------------ refusals and prints
def _calls(tmp_path, commit=COMMIT):
    return {"g0": ["g0", "--out", str(tmp_path / "z"), "--commit", commit],
            "gate": ["gate", "--out", str(tmp_path / "go"), "--g0", str(tmp_path / "g0"), "--commit", commit],
            "run": ["run", "--task", KEY1, "--seeds", "0-1", "--out", str(tmp_path / "o"), "--g0", str(tmp_path / "g0"),
                    "--gate", str(tmp_path / "gate"), "--commit", commit],
            "read": ["read", "--runs", str(tmp_path / "runs"), "--g0", str(tmp_path / "g0"), "--gate",
                     str(tmp_path / "gate"), "--out", str(tmp_path / "r"), "--commit", commit]}


def _asks(tmp_path, commit=COMMIT):
    return {"g0": ["refusals", "--step", "g0", "--commit", commit, "--out", str(tmp_path / "z")],
            "gate": ["refusals", "--step", "gate", "--commit", commit, "--out", str(tmp_path / "go"), "--g0",
                     str(tmp_path / "g0")],
            "run": ["refusals", "--step", "run", "--commit", commit, "--g0", str(tmp_path / "g0"), "--gate",
                    str(tmp_path / "gate")],
            "read": ["refusals", "--step", "read", "--commit", commit, "--out", str(tmp_path / "r")]}


def _nothing_loads(h35, monkeypatch):
    def boom(*a, **k):
        raise AssertionError("something loaded")
    for name in ("load_transition", "load_parts", "builder_context"):
        monkeypatch.setattr(h35, name, boom)
    for name in ("load_h29b_reading", "load_h27_rows"):
        monkeypatch.setattr(h35.h34, name, boom)


def _records(h35, tmp_path):
    for d, rec in (("g0", _g0(h35)), ("gate", _gate(h35))):
        (tmp_path / d).mkdir(exist_ok=True)
        (tmp_path / d / f"{d}.json").write_text(json.dumps(rec))


def test_every_step_refuses_without_a_formed_open_row_before_anything_loads(h35, monkeypatch, tmp_path, capsys):
    _nothing_loads(h35, monkeypatch)
    _records(h35, tmp_path)
    monkeypatch.setattr(h35, "git_state", lambda: {"commit": COMMIT, "tree_clean": True, "utc": "x"})
    calls, asks = _calls(tmp_path), _asks(tmp_path)
    led = tmp_path / "TRIES.md"
    monkeypatch.setattr(h35, "TRIES", led)
    for text, msg in ((_led(), "must hold exactly one row H35 of 11 cells"),
                      (_led(_row(), _row()), "must hold exactly one row H35 of 11 cells"),
                      (_led(_row(cells=10)), "must hold exactly one row H35 of 11 cells"),
                      (_led(_row(first=FIRST.replace("ten", "eight"))), "first cell is not the text"),
                      (_led(_row(outcome="OUTCOME H35: VOID: x; nothing is built.")), "already holds an outcome")):
        led.write_text(text)
        for step in calls:
            for args in (calls[step], asks[step]):
                with pytest.raises(SystemExit, match=msg):
                    h35.main(args)
    led.unlink()
    with pytest.raises(SystemExit, match="must hold exactly one row H35"):
        h35.main(calls["g0"])
    led.write_text(_led(_row()))
    for step in calls:
        assert h35.main(asks[step]) is None
        with pytest.raises(AssertionError, match="something loaded"):
            h35.main(calls[step])
    assert capsys.readouterr().out == "" and not any((tmp_path / d).exists() for d in ("z", "go", "o", "r"))


def test_the_other_refusals_of_every_step(h35, monkeypatch, tmp_path):
    _nothing_loads(h35, monkeypatch)
    _records(h35, tmp_path)
    led = tmp_path / "TRIES.md"
    led.write_text(_led(_row()))
    monkeypatch.setattr(h35, "TRIES", led)
    state = {"commit": COMMIT, "tree_clean": True, "utc": "x"}
    monkeypatch.setattr(h35, "git_state", lambda: dict(state))
    calls, asks = _calls(tmp_path), _asks(tmp_path)
    both = [*calls.values(), *asks.values()]
    state["tree_clean"] = False
    for args in both:
        with pytest.raises(SystemExit, match="H35 is refused: tracked files differ from the commit"):
            h35.main(args)
    state["tree_clean"], state["commit"] = True, "d" * 40
    for args in both:
        with pytest.raises(SystemExit, match="H35 is refused: HEAD is not the run commit"):
            h35.main(args)
    state["commit"] = COMMIT
    monkeypatch.setattr(h35, "vec_home", lambda: tmp_path)
    for args in both:
        with pytest.raises(SystemExit, match="vec is imported from"):
            h35.main(args)
    monkeypatch.setattr(h35, "vec_home", lambda: h35.HERE)
    # the gate and the run: G0 of this commit must have passed a form; the run: the gate of this commit said PASS
    g0 = tmp_path / "g0" / "g0.json"
    for rec in (_g0(h35, {f: "IDLE" for f in FORMS}), _g0(h35, commit="d" * 40), _g0(h35, clean=False),
                _g0(h35, verdict="VOID")):
        g0.write_text(json.dumps(rec))
        for args in (calls["gate"], asks["gate"], calls["run"], asks["run"]):
            with pytest.raises(SystemExit, match="has passed no form"):
                h35.main(args)
    g0.write_text(json.dumps(_g0(h35)))
    gate = tmp_path / "gate" / "gate.json"
    for rec in (_gate(h35, "FAIL"), _gate(h35, commit="d" * 40), _gate(h35, clean=False)):
        gate.write_text(json.dumps(rec))
        for args in (calls["run"], asks["run"]):
            with pytest.raises(SystemExit, match="has not said PASS"):
                h35.main(args)
    gate.write_text(json.dumps(_gate(h35)))
    (tmp_path / "o").mkdir()
    (tmp_path / "o" / "h35.json").write_text("{}")
    with pytest.raises(SystemExit, match="h35.json exists"):
        h35.main(calls["run"])
    with pytest.raises(SystemExit, match="--task must be one of"):
        h35.main(["run", "--task", "T2h/E8.25_late->E8.5", "--seeds", "0-1", "--out", str(tmp_path / "o2"), "--g0",
                  str(tmp_path / "g0"), "--gate", str(tmp_path / "gate"), "--commit", COMMIT])
    for d, f, step in (("z", "g0.json", "g0"), ("go", "gate.txt", "gate"), ("r", "reading.json", "read")):
        (tmp_path / d).mkdir(exist_ok=True)
        (tmp_path / d / f).write_text("x")
        for args in (calls[step], asks[step]):
            with pytest.raises(SystemExit, match="already holds"):
                h35.main(args)
    with pytest.raises(SystemExit, match="refusals --step gate needs --g0"):
        h35.main(["refusals", "--step", "gate", "--commit", COMMIT, "--out", str(tmp_path / "go2")])
    with pytest.raises(SystemExit, match="refusals --step run needs --g0 and --gate"):
        h35.main(["refusals", "--step", "run", "--commit", COMMIT, "--gate", str(tmp_path / "gate")])


def _prints(fn):
    tree = ast.parse(inspect.getsource(fn).lstrip())
    return [ast.unparse(n) for n in ast.walk(tree) if isinstance(n, ast.Call) and getattr(n.func, "id", "") == "print"]


def test_each_step_prints_only_its_lines(h35):
    run = sorted(_prints(h35.run_transition))
    assert len(run) == 2 and "s) done" in run[0] and run[1] == "print(f'wrote {path}', flush=True)"
    assert len(_prints(h35.run_gate)) == 4 and len(_prints(h35.run_g0)) == 4
    assert sorted(_prints(h35.cmd_read)) == ["print(res['follows'])", "print(res['line'])"]
    own = [fn for _, fn in inspect.getmembers(h35, inspect.isfunction) if fn.__module__ == h35.__name__]
    assert len(own) > 40
    for fn in own:
        if fn not in (h35.run_transition, h35.run_gate, h35.run_g0, h35.cmd_read):
            assert _prints(fn) == [], fn.__name__
    # the fixed words are fixed: no placeholder is left in them, and the ones that follow a read name what is built
    for words in (h35.SCREEN, h35.NOTHING, h35.BUILD, *h35.GATE_LINES.values(), h35.UNPRICED, h35.RULE):
        assert "{" not in words and "}" not in words, words
    assert "heart_trend_line_h.sh" in h35.BUILD and "row R55" in h35.BUILD and "nothing is built" in h35.NOTHING


def test_the_launchers_inline_outcome_reader_is_the_scripts(h35, tmp_path):
    """experiments/builds/heart_trend_line_h.sh reads row H35's outcome cell with an inline stdlib reader: on every
    ledger here it allows exactly the forms screen_outcome allows, and refuses where screen_outcome gives nothing or a
    VOID or STOP outcome."""
    import subprocess

    sh = (REPO_ROOT / "experiments" / "builds" / "heart_trend_line_h.sh").read_text()
    start = sh.index("python3 -I -B -c '\nimport re")
    code = sh[start + len("python3 -I -B -c '"):sh.index("' experiments/harness/TRIES.md", start)]
    cases = [_led(_row(outcome=h35.outcome_line(h35.verdict_of(dict(zip(FORMS, w, strict=True))),
                                                 dict(zip(FORMS, w, strict=True)))))
             for w in itertools.product(("KEPT", "DROP", "IDLE"), repeat=4)]
    cases += [_led(_row(outcome="OUTCOME H35: VOID: x; nothing is built.")), _led(_row()), _led(),
              _led(_row(), _row()), _led(_row(first=FIRST + " ")),
              _led(_row(outcome="OUTCOME H35: STOP: category KEPT; contrast IDLE; whole IDLE; portion IDLE.")),
              _led(_row(outcome="OUTCOME H35: category KEPT; contrast KEPT; whole KEPT; portion KEPT. OUTCOME H35: "
                                "category KEPT; contrast KEPT; whole KEPT; portion KEPT."))]
    seen = set()
    for i, text in enumerate(cases):
        p = tmp_path / f"t{i}.md"
        p.write_text(text)
        r = subprocess.run(["python3", "-I", "-B", "-c", code, str(p)], capture_output=True, text=True)
        want = h35.screen_outcome(text)
        if want is None or not want["allowed"]:
            assert r.returncode != 0 and r.stdout == "", (i, text)
        else:
            assert r.returncode == 0 and r.stdout.split() == want["allowed"], (i, text)
        seen.add(r.returncode)
    assert len(seen) >= 3  # allowed, nothing allowed, and a malformed or missing row are told apart


def test_the_text_of_a_read_holds_the_rule_the_unpriced_steps_and_the_tables(h35):
    res = _read(h35, score={"category": 0.5, "contrast": 0.0, "whole": -0.5, "portion": 0.2})
    txt = h35.text(res)
    lines = txt.splitlines()
    assert lines[0] == res["line"] and lines[1] == res["follows"]
    assert h35.RULE in txt and h35.UNPRICED in txt and "never set beside OTS125's headline" in txt
    assert "category (OTK125): KEPT; own +0.50" in txt and "whole (OTW125): DROP; own -0.50" in txt
    assert txt.count(h35.TAG) == 2 * 3 * (4 + 1 + 4 + 6)  # worlds x units x pairs
    void = h35.text(_read(h35, gate=None))
    assert "VOID, the conditions found, one by one:" in void and "v3: no gate record" in void
