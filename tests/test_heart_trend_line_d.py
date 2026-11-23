"""`experiments/builds/build_test.py heart --entry organ_scope --trend-strength 0.5 | 1.0 --z-mode raw --captain-ruling
R29`: the heart search line's step on file D (row R29 of experiments/harness/TRIES.md, the captain's order "i want it
too"), files OTS05 and OTS10: the organ cut D's command with --expression copy replaced by row H28's strength-scaled
trend step, before the size step on D's own draw, organ-cut pool, count rule, size procedure and seed.

On the synthetic stages and synthetic test board (its limit never shown) of tests/test_heart_trend_build.py, E10.5
relabelled as tests/test_build_heart_organ_scope.py's `organ_stage` and that file's size stand-in, with a temporary
ledger in place of TRIES.md: at strength 1.0 the expression is `vec.harness._b0_t2(ctx, trend=True)` on file D's
context bit for bit; at both strengths the step is row H28's own; the rows, coordinates, cell number, size field and,
after halve.py at seed 1000001, the half's rows are D's and only the expression differs; the same bytes twice; every
heart command that existed before this change, row R28's two included, gives the bytes and record it gave before
(pinned at main bf1dac7, before this option); every refusal, before anything is read; the R29 check on temporary
ledgers. Nothing is read from data/ and nothing is scored."""

from __future__ import annotations

import hashlib
import importlib.util
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from test_build_heart_iso_trim import _build, canon  # noqa: E402
from test_build_heart_organ_scope import ORGAN, _organ_mask, _use_stage, organ_stage  # noqa: E402
from test_heart_trend_build import SECRET_MAX, STAGES, heart_test_board, install, sha, synthetic_stage  # noqa: E402
from test_heart_trend_line import ARGV, BEFORE, R28_FIRST, TS  # noqa: E402

OTS = {0.5: ["--entry", "organ_scope", "--trend-strength", "0.5", "--z-mode", "raw", "--captain-ruling", "R29"],
       1.0: ["--entry", "organ_scope", "--trend-strength", "1.0", "--z-mode", "raw", "--captain-ruling", "R29"]}
R29_FIRST = ("R29 (the captain's order, \"i want it too\": a synthetic stand-in; ruling, no run, not a try; the heart "
             "builder is to take `--entry organ_scope` with `--trend-strength 0.5` and `--trend-strength 1.0` only "
             "with `--captain-ruling R29`)")
LEDGER = (f"| R27 (other) | T1 | x |\n| {R28_FIRST} | heart | the line's two files |\n"
          f"| {R29_FIRST} | heart | file D's two files |\n| H28 (row) | heart | x |\n")
R29_STRINGS = ("`--entry organ_scope`", "`--trend-strength 0.5`", "`--trend-strength 1.0`", "`--captain-ruling R29`")

# sha256 of (the .h5ad, its .build.json text, the record's `canon`) that build_test.py heart wrote for row R28's two
# commands on `install`'s synthetic stages and size stand-in, with LEDGER's row R28, computed with the builder of main
# at bf1dac7 (before this option) by the same builds as below. A record may differ from its text only in the last bits
# of floats. Every other heart command is pinned in BEFORE (tests/test_heart_trend_line.py) and rechecked here.
TS_BEFORE = {
    0.5: ("f3c4dc1de48e6fc5a140fa506add9446abd3ef569870d0d8c1eea07ccdc52b1d",
          "7d6a5d4af506161c94dd0ff70e6696c157cc36adfe4a33b2daf3581a53331a11",
          "f2b1eeb407d53a6550ad88ac9c3dcccc3a42c348cfa7b0854579189f4aae09d7"),
    1.0: ("c4230bae0315fb1bb488354b278ac68d0bc15ae9ce4e51cbef20dffd4fe5026b",
          "c393f55ec7c67fd9639e015591a57a5d5abdd98d5ededd56b7a32e66a303d4fa",
          "e50d2c30ac40828c127e5fc77362bc1a87d3f257605258e9d6857b7b5dac0fbb"),
}


@pytest.fixture(scope="module")
def bt():
    sys.path.insert(0, str(REPO / "experiments" / "builds"))
    import build_test

    return build_test


@pytest.fixture(scope="module")
def h28():
    spec = importlib.util.spec_from_file_location("h28_for_heart_trend_line_d",
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


def _organ_size(bt, monkeypatch) -> None:
    """tests/test_build_heart_organ_scope.py's size stand-in: scales by 1.25 and reports the core forecast's field."""
    from vec import t2_build

    def size(coords, size_json, hypothesis=None, **kw):
        C = np.asarray(coords, dtype=float)
        return C * 1.25, {"hypothesis": hypothesis, "factor": 1.25,
                          "target_rms_full": t2_build.rms_radius(C) * 1.25,
                          "target_source": "z_modes.raw.pred_rms_core", "warnings": []}

    monkeypatch.setattr(bt.t2_build, "size_or_hold", size)


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


def _ctx_d(bt):
    """File D's context with the trend stages, rebuilt here as the builder builds it: E10.5 is `organ_stage()` (what
    `organ_env` hands the builder), cut by vec.t2_build.organ_scope on its 3D coordinates, released-style labels and
    core (the builder's own call, the same mask as organ_scope_mask); "ref" is that organ-cut pool, "start" the
    synthetic E9.5, the three earlier stages HEART_TREND_STAGES' synthetic files; the test board's view and the
    builder's default seed. It is D's: its copy draw is D's expression (test below)."""
    from vec import harness, t2_build
    from vec.io import COORD_KEY

    stages = {f: synthetic_stage(*v) for f, v in STAGES.items()}
    ref = organ_stage()
    keep, _ = t2_build.organ_scope(np.asarray(ref.obsm[COORD_KEY], dtype=float)[:, :3],
                                   ref.obs[t2_build.T2_LABEL_KEY].to_numpy(), np.asarray(ref.obs["core"], dtype=bool))
    np.testing.assert_array_equal(keep, _organ_mask(ref))
    assert 0 < keep.sum() < ref.n_obs
    roles = {"ref": ("T2/E10.5.h5ad", 10.5), "start": ("T2/E9.5.h5ad", 9.5), **bt.HEART_TREND_STAGES}
    inputs = {"ref": ref[keep].copy(), "start": stages["E9.5.h5ad"],
              **{r: stages[Path(f).name] for r, (f, _) in bt.HEART_TREND_STAGES.items()}}
    view = bt._view("test:T2h", "T2", "extrap", "heart", heart_test_board().key, roles, "ref", 12.5)
    return harness.BuildContext(task=view, board=harness._candidate_board(heart_test_board()), inputs=inputs,
                                seed=harness.BUILD_SEED_OFFSET + 1)


def _X(path):
    from vec.io import dense

    return dense(ad.read_h5ad(path).X).astype(np.float32)


# ------------------------------------------------------------------------------------------------ off by default
def test_the_option_is_off_by_default(bt):
    a = bt.build_parser().parse_args(["heart", *ORGAN, "--out", "runs/x/h.h5ad"])
    assert a.trend_strength is None and a.captain_ruling is None
    assert bt.HEART_OTS_COMMANDS == {s: "heart " + " ".join(v) for s, v in OTS.items()}
    assert bt.HEART_OTS_ROW == "R29" and bt.HEART_OTS_FILES == {0.5: "OTS05", 1.0: "OTS10"}
    assert bt.HEART_ORGAN_SCOPE == "heart " + " ".join(ORGAN)  # D's command is untouched
    assert bt.HEART_TS_COMMANDS == {s: "heart " + " ".join(v) for s, v in TS.items()} and bt.HEART_TS_ROW == "R28"


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
    assert "trend_strength" not in rec and "R29" not in text and "OTS" not in text


@pytest.mark.parametrize("s", [0.5, 1.0])
def test_row_r28s_two_files_keep_their_bytes_and_record(bt, own_repo, ledger, monkeypatch, tmp_path, s):
    install(bt, monkeypatch.setattr, [])
    out, rec, text = _build(bt, tmp_path, f"ts{s:g}", *TS[s])
    assert sha(out) == TS_BEFORE[s][0]
    assert hashlib.sha256(text.encode()).hexdigest() == TS_BEFORE[s][1] or canon(text) == TS_BEFORE[s][2]
    assert rec["trend_strength"]["captain_ruling"] == "R28" and "R29" not in text and "OTS" not in text


# ------------------------------------------------------------------------------------------------ the step
def test_the_context_is_file_ds(bt, organ_env, tmp_path):
    """The rebuilt context's copy draw (vec.harness._copy, D's own step) is file D's expression bit for bit."""
    from vec import harness

    d, _, _ = _build(bt, tmp_path, "D", *ORGAN)
    np.testing.assert_array_equal(_X(d), harness._copy(_ctx_d(bt)).X)


def test_strength_one_is_the_harness_b0_trend_on_file_ds_context(bt, organ_env, tmp_path):
    """At strength 1.0 the expression is `harness._b0_t2(ctx, trend=True)` on file D's context bit for bit."""
    from vec import harness

    out, _, _ = _build(bt, tmp_path, "ots10", *OTS[1.0])
    want = harness._b0_t2(_ctx_d(bt), trend=True)
    X = _X(out)
    assert X.shape == want.X.shape
    np.testing.assert_array_equal(X, want.X)
    assert organ_env == ["boards", "E10.5.h5ad", "E9.5.h5ad", "E8.75.h5ad", "E8.5.h5ad", "E8.25_late.h5ad", "boards"]


@pytest.mark.parametrize("s", [0.5, 1.0])
def test_the_step_is_row_h28s_at_each_strength_on_ds_context(bt, h28, s):
    ctx = _ctx_d(bt)
    st, rec = bt._trend_strength_step(ctx, s)
    ref, _, _, ref_rec = h28.trend_expression(ctx, s)
    np.testing.assert_array_equal(st.X, ref.X)
    np.testing.assert_array_equal(st.coords, ref.coords)
    assert rec["budget_days"] == pytest.approx(0.5) and ref_rec["effective_budget_days"] == pytest.approx(0.5 * s)


@pytest.mark.parametrize("s", [0.5, 1.0])
def test_only_the_expression_differs_from_d(bt, organ_env, tmp_path, s):
    """Rows, coordinates, cell number, size field, D's own record fields and, after halve.py at seed 1000001, the
    half's rows are D's; the expression differs, and 0.5 differs from 1.0."""
    sys.path.insert(0, str(REPO / "experiments" / "builds"))
    import halve

    from vec.io import COORD_KEY

    d, d_rec, _ = _build(bt, tmp_path, "D", *ORGAN)
    t, t_rec, _ = _build(bt, tmp_path, f"ots{s:g}", *OTS[s])
    other, _, _ = _build(bt, tmp_path, "other", *OTS[1.5 - s])
    A, B = ad.read_h5ad(d), ad.read_h5ad(t)
    assert A.n_obs == B.n_obs and list(A.obs_names) == list(B.obs_names) and list(A.var_names) == list(B.var_names)
    np.testing.assert_array_equal(np.asarray(A.obsm[COORD_KEY]), np.asarray(B.obsm[COORD_KEY]))
    assert not np.array_equal(_X(d), _X(t)) and not np.array_equal(_X(t), _X(other))
    for k in ("n_cells", "seed", "count_decision", "trim", "size", "z_mode", "board", "dry_run", "depth_step",
              "scope_option"):
        assert t_rec[k] == d_rec[k], k
    assert t_rec["size"]["hypothesis"] == "organ" and t_rec["scope_option"]["command"] == "heart " + " ".join(ORGAN)
    hd, ht = tmp_path / "d_half.h5ad", tmp_path / "t_half.h5ad"
    halve.halve(str(d), str(hd))
    halve.halve(str(t), str(ht))
    HA, HB = ad.read_h5ad(hd), ad.read_h5ad(ht)
    assert halve.DEFAULT_SEED == 1000001
    assert list(HA.obs_names) == list(HB.obs_names) and HA.n_obs == HB.n_obs
    np.testing.assert_array_equal(np.asarray(HA.obsm[COORD_KEY]), np.asarray(HB.obsm[COORD_KEY]))
    assert not np.array_equal(_X(hd), _X(ht))


@pytest.mark.parametrize("s", [0.5, 1.0])
def test_the_same_bytes_twice(bt, organ_env, tmp_path, s):
    a, _, ta = _build(bt, tmp_path, "a", *OTS[s])
    b, _, tb = _build(bt, tmp_path, "b", *OTS[s])
    assert a.read_bytes() == b.read_bytes() and ta == tb


@pytest.mark.parametrize("s", [0.5, 1.0])
def test_the_record_names_option_strength_file_ruling_order_and_command(bt, organ_env, tmp_path, capsys, s):
    _, rec, text = _build(bt, tmp_path, "t", *OTS[s])
    ts = rec["trend_strength"]
    assert ts["option"] == "--trend-strength" and ts["strength"] == s and ts["captain_ruling"] == "R29"
    assert ts["file"] == {0.5: "OTS05", 1.0: "OTS10"}[s]
    assert "row R29" in ts["order"] and "i want it too" in ts["order"]
    assert ts["command"] == "heart " + " ".join(OTS[s]) and "D's draw" in ts["step"]
    assert ts["budget_days"] == pytest.approx(0.5) and ts["test_budget_days"] == pytest.approx(0.5 * s)
    assert ts["trend"]["rule"] == "vec.harness.B0_TREND_RULE"
    assert rec["scope_option"]["command"] == bt.HEART_ORGAN_SCOPE and "R21" in rec["scope_option"]["ruling"]
    assert rec["entry"] == f"trend strength {s:g} organ_scope (organ)" and "expression" not in rec
    assert rec["z_mode"] == "raw" and rec["n_cells"] == "max_cells" and "trend" not in rec and "R28" not in text
    assert str(SECRET_MAX) not in text + capsys.readouterr().out


# ------------------------------------------------------------------------------------------------ refusals
ORG = ["--entry", "organ_scope"]


@pytest.mark.parametrize("argv", [
    *([*ORG, "--trend-strength", v, "--z-mode", "raw", "--captain-ruling", "R29"] for v in ("0.25", "1", "0", "1.5",
                                                                                            "0.50", "-1.0")),
    [*ORG, "--trend-strength", "0.5", "--captain-ruling", "R29"],
    [*ORG, "--trend-strength", "1.0", "--z-mode", "recorded", "--captain-ruling", "R29"],
    [*OTS[0.5], "--expression", "copy"], [*OTS[1.0], "--expression", "copy"],
    [*OTS[0.5], "--trend"], [*OTS[0.5], "--depth"], [*OTS[0.5], "--herd"], [*OTS[0.5], "--cfm"],
    [*OTS[1.0], "--trend"], [*OTS[1.0], "--depth"], [*OTS[1.0], "--herd"], [*OTS[1.0], "--cfm"],
    *([*e, "--trend-strength", s, "--z-mode", "raw", "--captain-ruling", "R29"]
      for e in (["--entry", "s1"], ["--entry", "s2"], ["--entry", "s2_trim"], ["--entry", "core_trim"])
      for s in ("0.5", "1.0")),
    *([*ORG, "--trend-strength", s, "--z-mode", "raw", "--captain-ruling", r]
      for s in ("0.5", "1.0") for r in ("R28", "R21", "R30", "r29")),
    [*ORG, "--trend-strength", "0.5", "--z-mode", "raw"], [*ORG, "--trend-strength", "1.0", "--z-mode", "raw"],
    [*ORGAN, "--captain-ruling", "R29"], [*ORG, "--z-mode", "raw", "--captain-ruling", "R29"],
    ["--entry", "s1", "--trend", "--captain-ruling", "R29"]])
def test_every_other_use_is_refused_before_anything_is_read(bt, organ_env, tmp_path, argv):
    with pytest.raises(SystemExit) as e:
        bt.main(["heart", *argv, "--out", str(tmp_path / "h.h5ad")])
    assert "refused" in str(e.value) or "needs --captain-ruling" in str(e.value) or "goes with" in str(e.value)
    assert organ_env == [] and not (tmp_path / "h.h5ad").exists()


@pytest.mark.parametrize("board", ["T2:heart:val_extrap", "T2:heart:val_interp", "T2:embryo:test_interp"])
@pytest.mark.parametrize("s", [0.5, 1.0])
def test_any_other_board_is_refused_before_anything_is_read(bt, organ_env, tmp_path, board, s):
    with pytest.raises(SystemExit) as e:
        bt.main(["heart", *OTS[s], "--board", board, "--out", str(tmp_path / "h_dry.h5ad")])
    # a heart validation board meets the line's own refusal; another task's board the builder's earlier one
    assert "test board only" in str(e.value) or (board.startswith("T2:embryo") and "takes only" in str(e.value))
    assert organ_env == [] and not (tmp_path / "h_dry.h5ad").exists()


def test_the_test_boards_own_name_is_its_default(bt, organ_env, tmp_path):
    out, rec, _ = _build(bt, tmp_path, "named", *OTS[0.5], "--board", "T2:heart:test_extrap")
    default, _, _ = _build(bt, tmp_path, "default", *OTS[0.5])
    assert out.read_bytes() == default.read_bytes() and rec["dry_run"] is False


# ------------------------------------------------------------------------------------------------ the R29 check
def _use(bt, monkeypatch, path: Path, text: str) -> None:
    path.write_text(text)
    monkeypatch.setattr(bt, "TRIES_FILE", path)


def test_the_r29_check_passes_on_a_ledger_of_the_required_form(bt, monkeypatch, tmp_path):
    _use(bt, monkeypatch, tmp_path / "ok.md", LEDGER)
    assert bt._heart_d_line_ordered() and bt._heart_line_ordered()
    _use(bt, monkeypatch, tmp_path / "only.md", f"| {R29_FIRST} | heart | x |\n")
    assert bt._heart_d_line_ordered() and not bt._heart_line_ordered()  # row R28's own check is its own


@pytest.mark.parametrize("drop", [*R29_STRINGS, "ruling, no run", "the captain's order"])
def test_the_r29_check_fails_without_each_required_string(bt, monkeypatch, tmp_path, drop):
    first = R29_FIRST.replace(drop, "something else")
    assert first != R29_FIRST
    _use(bt, monkeypatch, tmp_path / "t.md", LEDGER.replace(R29_FIRST, first))
    assert not bt._heart_d_line_ordered() and bt._heart_line_ordered()


@pytest.mark.parametrize("moved", R29_STRINGS)
def test_the_r29_check_needs_the_strings_in_the_first_cell(bt, monkeypatch, tmp_path, moved):
    first = R29_FIRST.replace(moved, "something else")
    _use(bt, monkeypatch, tmp_path / "t.md", f"| {first} | heart | {moved} | {' '.join(R29_STRINGS)} |\n")
    assert not bt._heart_d_line_ordered()


def test_the_r29_check_fails_with_two_rows_or_another_rows_number(bt, monkeypatch, tmp_path):
    _use(bt, monkeypatch, tmp_path / "two.md", LEDGER + f"| {R29_FIRST} | heart | y |\n")
    assert not bt._heart_d_line_ordered()
    _use(bt, monkeypatch, tmp_path / "other.md", LEDGER.replace("| R29 (", "| R30 ("))
    assert not bt._heart_d_line_ordered()
    _use(bt, monkeypatch, tmp_path / "r28.md", LEDGER.replace(R29_FIRST, R29_FIRST.replace("R29", "R28")))
    assert not bt._heart_d_line_ordered()
    _use(bt, monkeypatch, tmp_path / "none.md", "| R27 (x) | T1 | x |\n")
    assert not bt._heart_d_line_ordered()


@pytest.mark.parametrize("text", [
    "| R27 (other) | T1 | x |\n",
    f"| {R28_FIRST} | heart | x |\n",
    LEDGER + f"| {R29_FIRST} | heart | y |\n",
    LEDGER.replace(R29_FIRST, R29_FIRST.replace("`--entry organ_scope`", "`--entry core_trim`")),
    LEDGER.replace(R29_FIRST, R29_FIRST.replace("`--trend-strength 1.0`", "`--trend-strength 1`")),
    LEDGER.replace(R29_FIRST, R29_FIRST.replace("`--captain-ruling R29`", "`--captain-ruling R28`"))])
def test_the_build_is_refused_on_a_ledger_without_the_order(bt, organ_env, ledger, tmp_path, text):
    ledger.write_text(text)
    for s in (0.5, 1.0):
        with pytest.raises(SystemExit) as e:
            bt.main(["heart", *OTS[s], "--out", str(tmp_path / "h.h5ad")])
        assert "must hold exactly one row R29" in str(e.value)
    assert organ_env == [] and not (tmp_path / "h.h5ad").exists()


def test_row_r28s_line_is_refused_on_a_ledger_with_only_row_r29(bt, organ_env, ledger, tmp_path):
    ledger.write_text(f"| {R29_FIRST} | heart | x |\n")
    for s in (0.5, 1.0):
        with pytest.raises(SystemExit) as e:
            bt.main(["heart", *TS[s], "--out", str(tmp_path / "h.h5ad")])
        assert "must hold exactly one row R28" in str(e.value)
    assert organ_env == [] and not (tmp_path / "h.h5ad").exists()
