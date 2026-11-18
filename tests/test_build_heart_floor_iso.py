"""`experiments/builds/build_test.py heart --expression copy --z-mode raw` (the captain's direction for the heart's
second upload, docs/PLAN.md decision log 2026-11-18 (j)): E10.5's cells with their expression unchanged, scaled by one
factor on all three axes. On synthetic stages, a synthetic test board whose limit is never shown and a synthetic size
JSON: the defaults route to the old calls and write the old bytes, the copy path returns the copy draw's expression
unchanged, the raw mode keeps the per-axis spread ratios and reaches the size record's target, and the refusals.
Nothing is read from data/ and nothing is scored."""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from test_heart_trend_build import (  # noqa: E402  (the synthetic heart stages, board and builder stand-ins)
    GENES,
    PRE_TREND_SHA256,
    SECRET_MAX,
    STAGES,
    heart_test_board,
    install,
    sha,
    synthetic_stage,
)

SEED = 1_000_001


@pytest.fixture(scope="module")
def bt():
    sys.path.insert(0, str(REPO / "experiments" / "builds"))
    import build_test

    return build_test


@pytest.fixture(scope="module")
def sp():
    """experiments/coords/size_procedure.py, loaded by path (as tests/test_t2_build.py does)."""
    spec = importlib.util.spec_from_file_location("coords_size_procedure_iso",
                                                  REPO / "experiments" / "coords" / "size_procedure.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture
def heart_env(bt, monkeypatch, own_repo):
    """The synthetic heart environment of tests/test_heart_trend_build.py, with the size stand-in also recording the
    keyword arguments each call passed."""
    reads: list[str] = []
    install(bt, monkeypatch.setattr, reads)
    stand_in = bt.t2_build.size_or_hold
    calls: list[dict] = []

    def size(coords, size_json, hypothesis=None, **kw):
        calls.append({"hypothesis": hypothesis, **kw})
        return stand_in(coords, size_json, hypothesis, **kw)

    monkeypatch.setattr(bt.t2_build, "size_or_hold", size)
    return reads, calls


def _build(bt, tmp_path, name: str, *argv: str) -> tuple[Path, dict, str]:
    out = tmp_path / name / "T2heart.h5ad"
    bt.main(["heart", *argv, "--out", str(out)])
    text = out.with_suffix(".build.json").read_text()
    return out, json.loads(text), text


def _ctx(ref: ad.AnnData, start: ad.AnnData):
    from vec import harness

    days = {"ref": 10.5, "start": 9.5}
    view = harness.TaskView(key="expected", task="T2", kind="extrap", setting="heart",
                            test_board=heart_test_board().key, roles=tuple(days), role_days=days,
                            role_stems={r: r for r in days}, reference_role="ref", target_day=float("nan"))
    return harness.BuildContext(task=view, board=harness._candidate_board(heart_test_board()),
                                inputs={"ref": ref, "start": start}, seed=SEED)


# ---------------------------------------------------------------------------------------------------- defaults
def test_the_new_options_are_off_by_default(bt):
    a = bt.build_parser().parse_args(["heart", "--entry", "s1", "--out", "runs/x/h.h5ad"])
    assert a.expression == "b0" and a.z_mode == "recorded"
    a = bt.build_parser().parse_args(["heart", "--entry", "s1", "--expression", "copy", "--z-mode", "raw",
                                      "--out", "runs/x/h.h5ad"])
    assert a.expression == "copy" and a.z_mode == "raw"


@pytest.mark.parametrize("entry", ["s1", "s2_trim"])
def test_defaults_write_the_old_bytes_and_make_the_old_size_call(bt, heart_env, tmp_path, capsys, entry):
    """With the defaults: the file is byte for byte the one pinned before --trend, the size step gets no z_mode
    keyword (z_mode=None, the recorded mode), and the record has neither new key."""
    reads, calls = heart_env
    out, rec, text = _build(bt, tmp_path, "plain", "--entry", entry)
    assert sha(out) == PRE_TREND_SHA256[entry]
    assert reads == ["E10.5.h5ad", "E9.5.h5ad"]
    assert len(calls) == 1 and "z_mode" not in calls[0]
    assert calls[0]["hypothesis"] == {"s1": "H1", "s2_trim": "H2"}[entry]
    assert "expression" not in rec and "z_mode" not in rec and rec["entry"].startswith("B0 ")
    assert str(SECRET_MAX) not in text + capsys.readouterr().out


def test_the_default_trend_build_is_unchanged(bt, heart_env, tmp_path):
    """--trend with the defaults gives the same file and record as --trend --expression b0 --z-mode recorded."""
    a, rec_a, text_a = _build(bt, tmp_path, "a", "--entry", "s1", "--trend")
    b, rec_b, text_b = _build(bt, tmp_path, "b", "--entry", "s1", "--trend", "--expression", "b0",
                              "--z-mode", "recorded")
    assert a.read_bytes() == b.read_bytes() and text_a == text_b
    assert "z_mode" not in heart_env[1][0] and "z_mode" not in heart_env[1][1]


# ---------------------------------------------------------------------------------------------------- copy
def test_copy_returns_the_copy_draw_expression_unchanged(bt, heart_env, tmp_path, capsys):
    """--expression copy writes the copy baseline: the cells draw_copy picks from the full E10.5 (the draw B0 starts
    from, same context, seed and count) with their released expression, no thinning and no top-up; the coordinates are
    the draw's, then the size step (stand-in: x 1.25)."""
    from vec import harness
    from vec.io import COORD_KEY, dense

    out, rec, text = _build(bt, tmp_path, "copy", "--entry", "s1", "--expression", "copy")
    stages = {f: synthetic_stage(*v) for f, v in STAGES.items()}
    ref, start = stages["E10.5.h5ad"], stages["E9.5.h5ad"]
    ctx = _ctx(ref, start)
    pool, rows, X, coords = harness.draw_copy(ctx)
    got = ad.read_h5ad(out)
    Xg = dense(got.X).astype(np.float32)
    np.testing.assert_array_equal(Xg, X)
    np.testing.assert_array_equal(Xg, dense(ref[:, GENES].X[rows]).astype(np.float32))  # the released values
    np.testing.assert_array_equal(np.asarray(got.obsm[COORD_KEY], dtype=np.float32)[:, :3],
                                  (np.asarray(coords, dtype=np.float32).astype(float) * 1.25).astype(np.float32))
    b0 = harness._b0_t2(ctx)
    assert Xg.shape == b0.X.shape and not np.array_equal(Xg, b0.X)
    np.testing.assert_array_equal(np.asarray(b0.coords), np.asarray(coords, dtype=np.float32))  # B0's own draw
    assert rec["expression"] == "copy" and rec["entry"] == "copy s1 (H1)" and "z_mode" not in rec
    assert rec["n_cells"] == "max_cells" and "trend" not in rec
    assert str(SECRET_MAX) not in text + capsys.readouterr().out
    assert "z_mode" not in heart_env[1][0]


@pytest.mark.parametrize("extra", [["--entry", "s2"], ["--entry", "s2_trim"], ["--entry", "s1", "--trend"],
                                   ["--entry", "s1", "--depth"], ["--entry", "s1", "--herd"],
                                   ["--entry", "s1", "--cfm"]])
def test_copy_is_refused_outside_s1_and_with_any_step_before_reading_anything(bt, heart_env, tmp_path, extra):
    with pytest.raises(SystemExit) as e:
        bt.main(["heart", *extra, "--expression", "copy", "--out", str(tmp_path / "h.h5ad")])
    assert "--expression copy is refused" in str(e.value)
    assert heart_env[0] == [] and not (tmp_path / "h.h5ad").exists()


# ---------------------------------------------------------------------------------------------------- raw z mode
def test_raw_mode_is_passed_to_the_size_step_and_recorded(bt, heart_env, tmp_path):
    _, rec, _ = _build(bt, tmp_path, "raw", "--entry", "s1", "--expression", "copy", "--z-mode", "raw")
    assert heart_env[1] == [{"hypothesis": "H1", "source_stage": "E10.5", "reference_cloud": heart_env[1][0][
        "reference_cloud"], "z_mode": "raw", "expect_board": "heart"}]
    assert rec["z_mode"] == "raw" and rec["expression"] == "copy"


def _stage_cloud(scale, sd, seed):
    """As tests/test_t2_build.py: a core blob inside a wider shell, z on 15-unit sections."""
    rng = np.random.default_rng(seed)
    core = rng.normal(size=(1500, 3)) * np.asarray(sd) * 0.6
    shell = rng.normal(size=(1500, 3)) * np.asarray(sd)
    C = np.vstack([core, shell]) * scale
    C[:, 2] = np.round(C[:, 2] / 15.0) * 15.0
    return C, np.r_[np.ones(1500, bool), np.zeros(1500, bool)]


def test_raw_mode_scales_all_axes_by_one_factor_to_the_target(bt, sp, monkeypatch, own_repo, tmp_path):
    """The real size step on a synthetic size JSON that recorded xy_grow_z_hold: with --z-mode raw the written cloud is
    the copy draw's cloud times one factor (per-axis standard-deviation ratios kept) at the RMS radius the JSON's raw
    mode forecasts for H1; with the default the z axis is held and the ratios move."""
    from vec import harness, t2_build
    from vec.io import COORD_KEY

    sd = (100.0, 80.0, 30.0)  # core thin along z: the procedure records xy_grow_z_hold
    rng = np.random.default_rng(0)
    scales = np.cumprod(1.0 + rng.uniform(0.05, 0.6, size=3))
    clouds = [(f"s{i}", t, *_stage_cloud(s, sd, i)) for i, (t, s) in enumerate(zip((8.75, 9.5, 10.5), scales,
                                                                                     strict=True))]
    geoms = [sp.geometry_from_arrays(*c) for c in clouds]
    rep = sp.heart_report(geoms, 12.5)
    size = {"board": "heart", "target_t": 12.5, "timestamp": "2026-10-19T00:00:00Z", "valid": rep["valid"],
            "preregistered_settings": dict(sp.PREREGISTERED),
            "provenance": {"git_commit": "synthetic", "git_uncommitted": [], "size_procedure_sha256": "synthetic"},
            "inputs": [], **rep}
    assert size["valid"] is True and size["preregistered"]["z_mode"] == "xy_grow_z_hold"
    path = tmp_path / "size_heart.json"
    path.write_text(sp.dumps(size))

    real_size = t2_build.size_or_hold  # taken before install() puts its stand-in in place
    reads: list[str] = []
    install(bt, monkeypatch.setattr, reads)
    monkeypatch.setattr(bt.t2_build, "size_or_hold", real_size)  # the real size step
    monkeypatch.setitem(bt.SIZES, "heart", path)
    ref = synthetic_stage(3000, 2, 0.0, 60)  # E10.5 in the JSON's measured cloud
    ref.obsm["spatial_3D"] = clouds[-1][2].astype(np.float32)
    ref.obs["core"] = clouds[-1][3]
    start = synthetic_stage(*STAGES["E9.5.h5ad"])
    load = bt.load_stage
    monkeypatch.setattr(bt, "load_stage", lambda p, **kw: ref.copy() if Path(p).name == "E10.5.h5ad"
                        else load(p, **kw))

    out, rec, _ = _build(bt, tmp_path, "raw", "--entry", "s1", "--expression", "copy", "--z-mode", "raw")
    _, coords_draw = harness.draw_copy(_ctx(ref, start))[2:]
    D = np.asarray(coords_draw, dtype=np.float32).astype(float)
    G = np.asarray(ad.read_h5ad(out).obsm[COORD_KEY], dtype=float)[:, :3]
    sd_draw, sd_out = D.std(0), G.std(0)
    np.testing.assert_allclose(sd_out / sd_out[0], sd_draw / sd_draw[0], rtol=1e-5)
    target = size["z_modes"]["raw"]["scope"]["H1"]["pred_rms_full"]
    assert t2_build.rms_radius(G) == pytest.approx(target, rel=1e-5)
    assert rec["size"]["target_rms_full"] == pytest.approx(target) and rec["size"]["z_mode_used"] == "raw"
    assert rec["size"]["z_mode_recorded"] == "xy_grow_z_hold" and rec["size"]["z_held"] is False
    assert rec["size"]["hypothesis"] == "H1" and rec["z_mode"] == "raw"

    held, rec_h, _ = _build(bt, tmp_path, "held", "--entry", "s1", "--expression", "copy")
    H = np.asarray(ad.read_h5ad(held).obsm[COORD_KEY], dtype=float)[:, :3]
    assert rec_h["size"]["z_mode_used"] == "xy_grow_z_hold" and "z_mode" not in rec_h
    np.testing.assert_allclose(H[:, 2], D[:, 2], rtol=1e-6)  # z held at the draw's
    assert not np.allclose(H.std(0) / H.std(0)[0], sd_draw / sd_draw[0], rtol=1e-3)
