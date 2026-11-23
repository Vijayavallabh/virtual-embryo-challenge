"""`experiments/builds/build_test.py heart --entry core_trim --trend-strength 0.5 | 1.0 --z-mode raw --captain-ruling
R28`: the heart search line's two files TS05 and TS10 (row R28 of experiments/harness/TRIES.md, the captain's order of
2026-11-23, "start the heart search"), the full trim C's command with row H28's trend step at that strength. On the
synthetic stages, synthetic test board (its limit never shown) and size stand-in of tests/test_heart_trend_build.py,
with a temporary ledger in place of TRIES.md: at strength 1.0 the expression is `vec.harness._b0_t2(ctx, trend=True)` on
the core-trim context bit for bit (row H28's V5), and at both strengths the step is row H28's own; the rows,
coordinates, cell number, size field and half's rows are C's and only the expression differs; the same bytes twice;
every heart command that existed before gives the bytes and record it gave before this option (pinned at main
f829663); every refusal, before anything is read; the ruling check on temporary ledgers. Nothing is read from data/
and nothing is scored."""

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
from test_build_heart_core_trim import CORE_TRIM  # noqa: E402
from test_build_heart_iso_trim import _build, canon  # noqa: E402
from test_build_heart_organ_scope import _use_stage, organ_stage  # noqa: E402
from test_heart_trend_build import SECRET_MAX, STAGES, heart_test_board, install, sha, synthetic_stage  # noqa: E402

TS = {0.5: ["--entry", "core_trim", "--trend-strength", "0.5", "--z-mode", "raw", "--captain-ruling", "R28"],
      1.0: ["--entry", "core_trim", "--trend-strength", "1.0", "--z-mode", "raw", "--captain-ruling", "R28"]}
R28_FIRST = ("R28 (the captain's order of 2026-11-23, \"start the heart search\": a synthetic stand-in; ruling, no "
             "run, not a try; the heart builder is to take `--trend-strength 0.5` and `--trend-strength 1.0` only with "
             "`--captain-ruling R28`)")
LEDGER = f"| R27 (other) | T1 | x |\n| {R28_FIRST} | heart | the line's two files |\n| H28 (row) | heart | x |\n"

# sha256 of (the .h5ad, its .build.json text, the record's `canon`) that build_test.py heart wrote for every heart
# command that existed before --trend-strength, on `install`'s synthetic stages and size stand-in (the organ cut on
# test_build_heart_organ_scope's organ_stage and its size stand-in), computed with the builder of main at f829663
# (before this option) by the same builds as below. A record may differ from its text only in the last bits of floats.
BEFORE = {
    "s1": ("cf887247be3d77a3fbbb6a99848a080137481f45a59bb97e4be3b3187e02e31e",
           "7772c7ccfcee8d33914afc0bc38ef25e38f861ca32505d78274176fd268bcd32",
           "2b0674350ca8655ff0ab235b8a7c8284d5a42eccaae165db62461575c32e27ef"),
    "s2": ("cf887247be3d77a3fbbb6a99848a080137481f45a59bb97e4be3b3187e02e31e",
           "0d160a305dca2bf52a7ac7fd43030f2c04753e302c2d2be62138afd7adc597ab",
           "4bc39f7170f8181e5a4d81fbe7378cb1ac79c4a409c38a886c7d656e37694f79"),
    "s2_trim": ("fc37df8923c87f39e6432eec45a2d8d49eae7c4e3b572fa057588bce949d5b9a",
                "167e4e407fb0079238c994b6e9c3270115898f2c7ee34b7e549c9909480ae6d1",
                "40accedc96b1df4335b4d3fb1af8c9a12d8076df28d923ed27da57729a419728"),
    "s1_trend": ("187673f9dabbe90d49d0c34af581fcb1431b6bfe2dc6fe01ff9f55d95698b0b4",
                 "0f3304eaced7e5cea546fefc0ba6c4ff87679a184d85a4d694589f046a889d2c",
                 "9b6d79152a5e3adec4c059ebd17f3b08947fb42bada7a1cf1a7c7e96b2779497"),
    "s2_trim_trend": ("4014fb74785059eea40fa7c6e9d7b22be05a3747af1162c1aa9fd1f26b24643b",
                      "dad93fdbfc321a6cf82d7941cc8a8607bea28084d0308e0c62cf80d03b6c6d90",
                      "788ecb4d87720fb113b8bec85653a1489051500cb32e4643132a1621ad4361ec"),
    "s1_depth": ("cf887247be3d77a3fbbb6a99848a080137481f45a59bb97e4be3b3187e02e31e",
                 "af9da2076c9d3d7f9cd2577237087752c11e7c30547a72aee5b3b322257bcf23",
                 "2b5322c5e6ee7549d1fa3184d6dd3b7a1242d36e198aa178f36cc602f6fd7ef0"),
    "s1_herd": ("74cbd5eca79ecf5a02af58209ecfac679bd66fd52ab2e12d42eb170ccad60037",
                "7772c7ccfcee8d33914afc0bc38ef25e38f861ca32505d78274176fd268bcd32",
                "2b0674350ca8655ff0ab235b8a7c8284d5a42eccaae165db62461575c32e27ef"),
    "s1_cfm": ("04f206cad866d09b108baf9c65796a9913e92c231c56f9c7fd689202012eed3f",
               "14d171a67e67c77c87f4b5ecbd308cf7b519685fe87dd3b722e5041c4a253214",
               "4f4b151451f97a95c524336dc37823fdb15cb65105cb1f2c8a9857af19448fdf"),
    "s1_copy": ("606ccc9479eb9f181cd333931de54c949150e732229f8394a1cf83785e988f3e",
                "4b3adaa11240e4acaac035d764e4e1e1973f1af8d72fba8385a1541c353141e7",
                "961e7625fc56c426b85f566b1114e3dad9997c61e1b3528bd18d8b082b1b20cd"),
    "s1_copy_raw": ("606ccc9479eb9f181cd333931de54c949150e732229f8394a1cf83785e988f3e",
                    "ee4c9f8b97e51ec839dffb499efe2da8fbd9a8729373c4f38a4069f24e574b57",
                    "b41e33a1b0a7491e79b0bd4d40ef6364509482b5c37cb242289a65f4a1da1ed8"),
    "s2_trim_copy_raw": ("6287272082da44e01a11ed5f705e78852e355dd655ec5ba257be55d31db67bad",
                         "7e62d6013348cb4c109d862944c883312f9e45f5c03c39140e112ca0e095d7c7",
                         "e59e7dc7446ba99a50e98980ff350ed56e90903bb0e401f9ea46ba189bc8e849"),
    "core_trim_copy_raw": ("73e57af3bc859644198fc5061fbfc200fcd3742605a79838cf03fc113819f254",
                           "cf0e1d0a81cf38c8e5fa542947cd40fdebe0a5d8b4572f0d006d9be3e77bd52b",
                           "53b26a922696642e7c0098d595ed5d712d1f6e738eab6da7afdd70a884c22c55"),
    "organ_scope_copy_raw": ("3a304c42a368cb3f0325cf525ab3e4c93f1b678c218dbbff77e8739211d9b802",
                             "b38020f208dd146da4e61223a4dada4122192a126ab7e31c0c5b5415512f5278",
                             "a971e29b0a1647ba3975e3c4441133d01d6954d4f29b9f9c66b7304c3ebf11eb"),
}
ARGV = {
    "s1": ["--entry", "s1"], "s2": ["--entry", "s2"], "s2_trim": ["--entry", "s2_trim"],
    "s1_trend": ["--entry", "s1", "--trend"], "s2_trim_trend": ["--entry", "s2_trim", "--trend"],
    "s1_depth": ["--entry", "s1", "--depth"], "s1_herd": ["--entry", "s1", "--herd"],
    "s1_cfm": ["--entry", "s1", "--cfm"],
    "s1_copy": ["--entry", "s1", "--expression", "copy"],
    "s1_copy_raw": ["--entry", "s1", "--expression", "copy", "--z-mode", "raw"],
    "s2_trim_copy_raw": ["--entry", "s2_trim", "--expression", "copy", "--z-mode", "raw"],
    "core_trim_copy_raw": CORE_TRIM,
    "organ_scope_copy_raw": ["--entry", "organ_scope", "--expression", "copy", "--z-mode", "raw"],
}


@pytest.fixture(scope="module")
def bt():
    sys.path.insert(0, str(REPO / "experiments" / "builds"))
    import build_test

    return build_test


@pytest.fixture(scope="module")
def h28():
    spec = importlib.util.spec_from_file_location("h28_for_heart_trend_line",
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
def heart_env(bt, monkeypatch, own_repo, ledger):
    """install()'s stages, board and size stand-in; `reads` records every stage read and every board read."""
    reads: list[str] = []
    install(bt, monkeypatch.setattr, reads)
    boards = bt.load_boards

    def load_boards():
        reads.append("boards")
        return boards()

    monkeypatch.setattr(bt, "load_boards", load_boards)
    return reads


def _ctx(bt):
    """The builder's core-trim context with the trend stages, rebuilt here from the synthetic stages."""
    from vec import harness, t2_build
    from vec.io import COORD_KEY

    stages = {f: synthetic_stage(*v) for f, v in STAGES.items()}
    ref = stages["E10.5.h5ad"]
    keep = t2_build.radial_trim_mask(np.asarray(ref.obsm[COORD_KEY], dtype=float)[:, :3],
                                     np.asarray(ref.obs["core"], dtype=bool), 0.0)
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
    a = bt.build_parser().parse_args(["heart", *CORE_TRIM, "--out", "runs/x/h.h5ad"])
    assert a.trend_strength is None and a.captain_ruling is None
    assert bt.HEART_TS_COMMANDS == {s: "heart " + " ".join(v) for s, v in TS.items()}
    assert bt.HEART_TS_ROW == "R28" and bt.HEART_TS_FILES == {0.5: "TS05", 1.0: "TS10"}
    assert bt.HEART_CORE_TRIM == "heart " + " ".join(CORE_TRIM)  # C's command is untouched


@pytest.mark.parametrize("name", sorted(BEFORE))
def test_every_heart_command_before_the_option_keeps_its_bytes_and_record(bt, own_repo, ledger, monkeypatch,
                                                                         tmp_path, name):
    from vec import t2_build

    install(bt, monkeypatch.setattr, [])
    if name.startswith("organ"):
        _use_stage(bt, monkeypatch, organ_stage(), [])

        def size(coords, size_json, hypothesis=None, **kw):
            C = np.asarray(coords, dtype=float)
            return C * 1.25, {"hypothesis": hypothesis, "factor": 1.25,
                              "target_rms_full": t2_build.rms_radius(C) * 1.25,
                              "target_source": "z_modes.raw.pred_rms_core", "warnings": []}

        monkeypatch.setattr(bt.t2_build, "size_or_hold", size)
    out, rec, text = _build(bt, tmp_path, name, *ARGV[name])
    assert sha(out) == BEFORE[name][0]
    assert hashlib.sha256(text.encode()).hexdigest() == BEFORE[name][1] or canon(text) == BEFORE[name][2]
    assert "trend_strength" not in rec and "R28" not in text


# ------------------------------------------------------------------------------------------------ the step
def test_strength_one_is_the_harness_b0_trend_on_the_core_trim_context(bt, heart_env, tmp_path):
    """Row H28's V5 on the file: at strength 1.0 the expression is `harness._b0_t2(ctx, trend=True)` on the core-trim
    context bit for bit, and the coordinates are its own through the size stand-in."""
    from vec import harness
    from vec.io import COORD_KEY

    out, _, _ = _build(bt, tmp_path, "ts10", *TS[1.0])
    want = harness._b0_t2(_ctx(bt), trend=True)
    X = _X(out)
    assert X.shape == want.X.shape
    np.testing.assert_array_equal(X, want.X)
    np.testing.assert_array_equal(np.asarray(ad.read_h5ad(out).obsm[COORD_KEY], dtype=np.float32)[:, :3],
                                  (np.asarray(want.coords, dtype=float) * 1.25).astype(np.float32))
    assert heart_env == ["boards", "E10.5.h5ad", "E9.5.h5ad", "E8.75.h5ad", "E8.5.h5ad", "E8.25_late.h5ad", "boards"]


@pytest.mark.parametrize("s", [0.5, 1.0])
def test_the_step_is_row_h28s_at_each_strength(bt, h28, s):
    """The builder's step equals row H28's practice step (trend_expression, placebo off) on the same context, X and
    coordinates, and its record's budget is B0's (0.5 day at the test), the strength scaling it."""
    ctx = _ctx(bt)
    st, rec = bt._trend_strength_step(ctx, s)
    ref, _, _, ref_rec = h28.trend_expression(ctx, s)
    np.testing.assert_array_equal(st.X, ref.X)
    np.testing.assert_array_equal(st.coords, ref.coords)
    assert rec["budget_days"] == pytest.approx(0.5) and ref_rec["effective_budget_days"] == pytest.approx(0.5 * s)


@pytest.mark.parametrize("s", [0.5, 1.0])
def test_only_the_expression_differs_from_c(bt, heart_env, tmp_path, s):
    """Rows, coordinates, cell number, size field and, after halve.py at seed 1000001, the half's rows are C's; the
    expression differs, and 0.5 differs from 1.0."""
    sys.path.insert(0, str(REPO / "experiments" / "builds"))
    import halve

    from vec.io import COORD_KEY

    c, c_rec, _ = _build(bt, tmp_path, "C", *CORE_TRIM)
    t, t_rec, _ = _build(bt, tmp_path, f"ts{s:g}", *TS[s])
    other, _, _ = _build(bt, tmp_path, "other", *TS[1.5 - s])
    A, B = ad.read_h5ad(c), ad.read_h5ad(t)
    assert A.n_obs == B.n_obs and list(A.obs_names) == list(B.obs_names) and list(A.var_names) == list(B.var_names)
    np.testing.assert_array_equal(np.asarray(A.obsm[COORD_KEY]), np.asarray(B.obsm[COORD_KEY]))
    assert not np.array_equal(_X(c), _X(t)) and not np.array_equal(_X(t), _X(other))
    for k in ("n_cells", "seed", "count_decision", "trim", "size", "z_mode", "board", "dry_run", "depth_step"):
        assert t_rec[k] == c_rec[k], k
    assert t_rec["scope_option"]["size_source_field"] == c_rec["scope_option"]["size_source_field"] \
        == "z_modes.raw.pred_rms_core"
    keys = set(t_rec["scope_option"]) - {"command", "ruling"}
    assert {k: t_rec["scope_option"][k] for k in keys} == {k: c_rec["scope_option"][k] for k in keys}
    hc, ht = tmp_path / "c_half.h5ad", tmp_path / "t_half.h5ad"
    halve.halve(str(c), str(hc))
    halve.halve(str(t), str(ht))
    HA, HB = ad.read_h5ad(hc), ad.read_h5ad(ht)
    assert halve.DEFAULT_SEED == 1000001
    assert list(HA.obs_names) == list(HB.obs_names) and HA.n_obs == HB.n_obs
    np.testing.assert_array_equal(np.asarray(HA.obsm[COORD_KEY]), np.asarray(HB.obsm[COORD_KEY]))
    assert not np.array_equal(_X(hc), _X(ht))


@pytest.mark.parametrize("s", [0.5, 1.0])
def test_the_same_bytes_twice(bt, heart_env, tmp_path, s):
    a, _, ta = _build(bt, tmp_path, "a", *TS[s])
    b, _, tb = _build(bt, tmp_path, "b", *TS[s])
    assert a.read_bytes() == b.read_bytes() and ta == tb


@pytest.mark.parametrize("s", [0.5, 1.0])
def test_the_record_names_option_strength_ruling_order_and_command(bt, heart_env, tmp_path, capsys, s):
    _, rec, text = _build(bt, tmp_path, "t", *TS[s])
    ts = rec["trend_strength"]
    assert ts["option"] == "--trend-strength" and ts["strength"] == s and ts["captain_ruling"] == "R28"
    assert ts["file"] == {0.5: "TS05", 1.0: "TS10"}[s]
    assert "row R28" in ts["order"] and "start the heart search" in ts["order"]
    assert ts["command"] == "heart " + " ".join(TS[s])
    assert ts["budget_days"] == pytest.approx(0.5) and ts["test_budget_days"] == pytest.approx(0.5 * s)
    assert ts["trend"]["rule"] == "vec.harness.B0_TREND_RULE"
    assert rec["scope_option"]["command"] == ts["command"] and rec["scope_option"]["ruling"] == ts["order"]
    assert rec["entry"] == f"trend strength {s:g} core_trim (core)" and "expression" not in rec
    assert rec["z_mode"] == "raw" and rec["n_cells"] == "max_cells" and "trend" not in rec
    assert str(SECRET_MAX) not in text + capsys.readouterr().out


# ------------------------------------------------------------------------------------------------ refusals
@pytest.mark.parametrize("argv", [
    [*CORE_TRIM[:2], "--trend-strength", "0.25", "--z-mode", "raw", "--captain-ruling", "R28"],
    [*CORE_TRIM[:2], "--trend-strength", "1", "--z-mode", "raw", "--captain-ruling", "R28"],
    [*CORE_TRIM[:2], "--trend-strength", "0", "--z-mode", "raw", "--captain-ruling", "R28"],
    [*CORE_TRIM[:2], "--trend-strength", "1.5", "--z-mode", "raw", "--captain-ruling", "R28"],
    ["--entry", "core_trim", "--trend-strength", "0.5", "--captain-ruling", "R28"],
    ["--entry", "core_trim", "--trend-strength", "0.5", "--z-mode", "recorded", "--captain-ruling", "R28"],
    [*TS[0.5], "--expression", "copy"], [*TS[1.0], "--expression", "copy"],
    [*TS[0.5], "--trend"], [*TS[0.5], "--depth"], [*TS[0.5], "--herd"], [*TS[0.5], "--cfm"],
    [*TS[1.0], "--trend"], [*TS[1.0], "--depth"], [*TS[1.0], "--herd"], [*TS[1.0], "--cfm"],
    *([*e, "--trend-strength", "0.5", "--z-mode", "raw", "--captain-ruling", "R28"]
      for e in (["--entry", "s1"], ["--entry", "s2"], ["--entry", "s2_trim"], ["--entry", "organ_scope"])),
    ["--entry", "organ_scope", "--expression", "copy", "--trend-strength", "1.0", "--z-mode", "raw",
     "--captain-ruling", "R28"],
    ["--entry", "core_trim", "--trend-strength", "0.5", "--z-mode", "raw"],
    ["--entry", "core_trim", "--trend-strength", "1.0", "--z-mode", "raw", "--captain-ruling", "R27"],
    ["--entry", "core_trim", "--trend-strength", "1.0", "--z-mode", "raw", "--captain-ruling", "R29"],
    [*CORE_TRIM, "--captain-ruling", "R28"], ["--entry", "s1", "--trend", "--captain-ruling", "R28"],
    ["--entry", "core_trim", "--z-mode", "raw", "--captain-ruling", "R28"]])
def test_every_other_use_is_refused_before_anything_is_read(bt, heart_env, tmp_path, argv):
    with pytest.raises(SystemExit) as e:
        bt.main(["heart", *argv, "--out", str(tmp_path / "h.h5ad")])
    assert "refused" in str(e.value) or "needs --captain-ruling" in str(e.value) or "goes with" in str(e.value)
    assert heart_env == [] and not (tmp_path / "h.h5ad").exists()


@pytest.mark.parametrize("board", ["T2:heart:val_extrap", "T2:heart:val_interp"])
def test_a_dry_run_is_refused_before_anything_is_read(bt, heart_env, tmp_path, board):
    with pytest.raises(SystemExit) as e:
        bt.main(["heart", *TS[0.5], "--board", board, "--out", str(tmp_path / "h_dry.h5ad")])
    assert "test board only" in str(e.value)
    assert heart_env == [] and not (tmp_path / "h_dry.h5ad").exists()


def test_the_test_boards_own_name_is_its_default(bt, heart_env, tmp_path):
    out, rec, _ = _build(bt, tmp_path, "named", *TS[0.5], "--board", "T2:heart:test_extrap")
    default, _, _ = _build(bt, tmp_path, "default", *TS[0.5])
    assert out.read_bytes() == default.read_bytes() and rec["dry_run"] is False


# ------------------------------------------------------------------------------------------------ the ruling check
def _ledger_with(path: Path, text: str) -> Path:
    path.write_text(text)
    return path


@pytest.mark.parametrize("drop", ["`--trend-strength 0.5`", "`--trend-strength 1.0`", "`--captain-ruling R28`",
                                  "ruling, no run", "the captain's order"])
def test_the_ruling_check_fails_without_each_required_string(bt, monkeypatch, tmp_path, drop):
    monkeypatch.setattr(bt, "TRIES_FILE", _ledger_with(tmp_path / "ok.md", LEDGER))
    assert bt._heart_line_ordered()
    monkeypatch.setattr(bt, "TRIES_FILE", _ledger_with(tmp_path / "t.md", LEDGER.replace(drop, "something else")))
    assert not bt._heart_line_ordered()


def test_the_ruling_check_needs_the_strings_in_the_first_cell(bt, monkeypatch, tmp_path):
    first = R28_FIRST.replace(" `--trend-strength 0.5` and", "")
    text = f"| {first} | heart | `--trend-strength 0.5` |\n"
    monkeypatch.setattr(bt, "TRIES_FILE", _ledger_with(tmp_path / "t.md", text))
    assert not bt._heart_line_ordered()


def test_the_ruling_check_fails_with_two_rows_or_another_rows_number(bt, monkeypatch, tmp_path):
    monkeypatch.setattr(bt, "TRIES_FILE", _ledger_with(tmp_path / "two.md", LEDGER + f"| {R28_FIRST} | heart | y |\n"))
    assert not bt._heart_line_ordered()
    other = LEDGER.replace("| R28 (", "| R29 (")
    monkeypatch.setattr(bt, "TRIES_FILE", _ledger_with(tmp_path / "other.md", other))
    assert not bt._heart_line_ordered()
    monkeypatch.setattr(bt, "TRIES_FILE", _ledger_with(tmp_path / "none.md", "| R27 (x) | T1 | x |\n"))
    assert not bt._heart_line_ordered()


@pytest.mark.parametrize("text", [
    "| R27 (other) | T1 | x |\n",
    LEDGER + f"| {R28_FIRST} | heart | y |\n",
    LEDGER.replace("`--trend-strength 1.0`", "`--trend-strength 1`"),
    LEDGER.replace("`--captain-ruling R28`", "`--captain-ruling R29`")])
def test_the_build_is_refused_on_a_ledger_without_the_order(bt, heart_env, ledger, tmp_path, text):
    ledger.write_text(text)
    for s in (0.5, 1.0):
        with pytest.raises(SystemExit) as e:
            bt.main(["heart", *TS[s], "--out", str(tmp_path / "h.h5ad")])
        assert "must hold exactly one row R28" in str(e.value)
    assert heart_env == [] and not (tmp_path / "h.h5ad").exists()
