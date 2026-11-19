"""`experiments/builds/build_test.py heart --entry core_trim --expression copy --z-mode raw` (the full trim, branch
campaign/heart-scope-options; no ruling), beside the half trim (`--entry s2_trim --expression copy --z-mode raw`).
On synthetic stages, a synthetic test board whose limit is never shown and synthetic size JSONs: each option's pool is
radial_trim_mask's at its fraction; the full trim is sized to the raw mode's own core forecast (pred_rms_core of
z_modes.raw) and refuses any other mode and the fallback; each option's record names its scope rule, the field of its
kept share, its size field and the halving, and no cell number; the refusals. The byte-for-byte pins of every heart
output that existed before (s1, s2, s2_trim, the ruled heart2 command) are in tests/test_build_heart_iso_trim.py and
tests/test_heart_trend_build.py. Nothing is read from data/ and nothing is scored."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from test_build_heart_floor_iso import _stage_cloud  # noqa: E402
from test_build_heart_iso_trim import ISO_TRIM, _build, _ctx  # noqa: E402
from test_heart_trend_build import (  # noqa: E402
    SECRET_MAX,
    STAGES,
    install,
    sha,
    synthetic_stage,
)

CORE_TRIM = ["--entry", "core_trim", "--expression", "copy", "--z-mode", "raw"]


@pytest.fixture(scope="module")
def bt():
    sys.path.insert(0, str(REPO / "experiments" / "builds"))
    import build_test

    return build_test


@pytest.fixture(scope="module")
def sp():
    """experiments/coords/size_procedure.py, loaded by path (as tests/test_t2_build.py does)."""
    import importlib.util

    spec = importlib.util.spec_from_file_location("coords_size_procedure_core_trim",
                                                  REPO / "experiments" / "coords" / "size_procedure.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture
def heart_env(bt, monkeypatch, own_repo):
    reads: list[str] = []
    install(bt, monkeypatch.setattr, reads)
    stand_in = bt.t2_build.size_or_hold
    calls: list[dict] = []

    def size(coords, size_json, hypothesis=None, **kw):
        calls.append({"hypothesis": hypothesis, **kw})
        return stand_in(coords, size_json, hypothesis, **kw)

    monkeypatch.setattr(bt.t2_build, "size_or_hold", size)
    return reads, calls


def _trimmed(ref: ad.AnnData, frac: float) -> ad.AnnData:
    from vec import t2_build
    from vec.io import COORD_KEY

    keep = t2_build.radial_trim_mask(np.asarray(ref.obsm[COORD_KEY], dtype=float)[:, :3],
                                     np.asarray(ref.obs["core"], dtype=bool), frac)
    assert 0 < keep.sum() < ref.n_obs
    return ref[keep].copy()


def _size_json(sp, tmp_path: Path, **over) -> tuple[Path, dict, list]:
    """A synthetic heart size JSON from the real procedure on three synthetic stages (recorded mode xy_grow_z_hold),
    and the stages' clouds."""
    sd = (100.0, 80.0, 30.0)
    rng = np.random.default_rng(0)
    scales = np.cumprod(1.0 + rng.uniform(0.05, 0.6, size=3))
    clouds = [(f"s{i}", t, *_stage_cloud(s, sd, i)) for i, (t, s) in enumerate(zip((8.75, 9.5, 10.5), scales,
                                                                                     strict=True))]
    rep = sp.heart_report([sp.geometry_from_arrays(*c) for c in clouds], 12.5)
    size = {"board": "heart", "target_t": 12.5, "timestamp": "2026-10-19T00:00:00Z", "valid": rep["valid"],
            "preregistered_settings": dict(sp.PREREGISTERED),
            "provenance": {"git_commit": "synthetic", "git_uncommitted": [], "size_procedure_sha256": "synthetic"},
            "inputs": [], **rep, **over}
    path = tmp_path / "size_heart.json"
    path.write_text(sp.dumps(size))
    return path, size, clouds


# ------------------------------------------------------------------------------------------------ off by default
def test_the_options_are_off_by_default(bt):
    for entry in ("s1", "s2_trim", "core_trim"):
        a = bt.build_parser().parse_args(["heart", "--entry", entry, "--out", "runs/x/h.h5ad"])
        assert a.expression == "b0" and a.z_mode == "recorded" and not a.trend
    assert bt.HEART_CORE_TRIM == "heart " + " ".join(CORE_TRIM)
    assert bt.HEART_ISO_TRIM == "heart " + " ".join(ISO_TRIM)
    assert {k: v["frac"] for k, v in bt.HEART_SCOPE_OPTIONS.items()} == {"s2_trim": 0.5, "core_trim": 0.0}


def test_the_ruled_and_frozen_heart_commands_carry_no_scope_option(bt, heart_env, tmp_path):
    """s1 (first upload), s2_trim --trend (the frozen second), and the ruled heart2 command build with no
    scope_option in their record (their bytes are pinned in test_build_heart_iso_trim.py)."""
    for name, argv in (("s1", ["--entry", "s1", "--trend"]), ("s2", ["--entry", "s2_trim", "--trend"]),
                       ("heart2", ["--entry", "s1", "--expression", "copy", "--z-mode", "raw"])):
        _, rec, _ = _build(bt, tmp_path, name, *argv)
        assert "scope_option" not in rec


# ------------------------------------------------------------------------------------------------ the pools
@pytest.mark.parametrize("argv, frac, hyp", [(CORE_TRIM, 0.0, "core"), (ISO_TRIM, 0.5, "H2")])
def test_each_option_draws_from_radial_trim_masks_pool_at_its_fraction(bt, heart_env, tmp_path, capsys, argv, frac,
                                                                       hyp):
    """The file is the copy baseline's draw (vec.harness.draw_copy, same seed and count) from E10.5 under
    radial_trim_mask at the option's fraction, with the released expression; the size step gets the option's
    hypothesis, the raw mode and the full E10.5 cloud as its reference."""
    from vec import harness
    from vec.io import COORD_KEY, dense

    reads, calls = heart_env
    out, rec, text = _build(bt, tmp_path, "opt", *argv)
    ref, start = synthetic_stage(*STAGES["E10.5.h5ad"]), synthetic_stage(*STAGES["E9.5.h5ad"])
    pool = _trimmed(ref, frac)
    _, rows, X, coords = harness.draw_copy(_ctx(pool, start))
    got = ad.read_h5ad(out)
    np.testing.assert_array_equal(dense(got.X).astype(np.float32), X)
    np.testing.assert_array_equal(np.asarray(got.obsm[COORD_KEY], dtype=np.float32)[:, :3],
                                  (np.asarray(coords, dtype=np.float32).astype(float) * 1.25).astype(np.float32))
    assert rec["trim"]["kept_share"] == pytest.approx(pool.n_obs / ref.n_obs)
    assert rec["trim"]["xy_scope_after"] <= frac * rec["trim"]["xy_scope_before"] + 1e-12
    assert reads == ["E10.5.h5ad", "E9.5.h5ad"]
    assert len(calls) == 1 and calls[0]["hypothesis"] == hyp and calls[0]["z_mode"] == "raw"
    assert calls[0]["source_stage"] == "E10.5" and calls[0]["reference_cloud"].n_obs == ref.n_obs
    assert got.n_obs == SECRET_MAX  # the count rule does not read the pool
    assert rec["n_cells"] == "max_cells" and str(SECRET_MAX) not in text + capsys.readouterr().out


def test_the_full_trim_keeps_every_core_cell_and_is_inside_the_half_trim():
    ref = synthetic_stage(*STAGES["E10.5.h5ad"])
    from vec import t2_build
    from vec.io import COORD_KEY

    C, core = np.asarray(ref.obsm[COORD_KEY], dtype=float)[:, :3], np.asarray(ref.obs["core"], dtype=bool)
    full, half = t2_build.radial_trim_mask(C, core, 0.0), t2_build.radial_trim_mask(C, core, 0.5)
    assert full[core].all() and not (full & ~half).any() and full.sum() < half.sum()


# ------------------------------------------------------------------------------------------------ the records
@pytest.mark.parametrize("argv, entry, source, hyp", [
    (CORE_TRIM, "core_trim", "z_modes.raw.pred_rms_core", "core"),
    (ISO_TRIM, "s2_trim", "z_modes.raw.scope.H2.pred_rms_full", "H2")])
def test_each_record_names_rule_kept_share_field_size_field_and_halving(bt, heart_env, tmp_path, argv, entry, source,
                                                                       hyp):
    _, rec, _ = _build(bt, tmp_path, entry, *argv)
    opt = rec["scope_option"]
    assert opt["command"] == "heart " + " ".join(argv) and opt["ruling"] is None
    assert opt["scope_rule"] == rec["trim"]["rule"]
    assert opt["kept_share_field"] == "trim.kept_share" and "kept_share" in rec["trim"]
    assert opt["size_source_field"] == source and opt["size_hypothesis"] == hyp
    assert "halve.py" in opt["halving"] and "heart2" in opt["halving"]
    assert rec["entry"] == f"copy {entry} ({hyp})" and rec["expression"] == "copy" and rec["z_mode"] == "raw"
    assert not any(isinstance(v, int) and not isinstance(v, bool) for v in opt.values())


def test_the_half_of_the_full_trim_file(bt, heart_env, tmp_path):
    sys.path.insert(0, str(REPO / "experiments" / "builds"))
    import halve

    out, _, _ = _build(bt, tmp_path, "core", *CORE_TRIM)
    half = tmp_path / "half.h5ad"
    digest = halve.halve(str(out), str(half))
    full, got = ad.read_h5ad(out), ad.read_h5ad(half)
    assert got.n_obs == full.n_obs // 2 and digest == sha(half)
    assert list(got.obs_names) == list(full.obs_names[halve.half_rows(full.n_obs, halve.DEFAULT_SEED)])


# ------------------------------------------------------------------------------------------------ the size
def test_the_core_forecast_is_the_full_forecast_at_scope_ratio_one(sp, tmp_path):
    """In the size JSON every raw-mode scope hypothesis is the core forecast times exp(its log scope): the core
    forecast is that product at log scope 0 (full-over-core ratio one), the size the full trim takes."""
    _, size, _ = _size_json(sp, tmp_path)
    raw = size["z_modes"]["raw"]
    for h in ("H1", "H0", "H2"):
        assert raw["scope"][h]["pred_rms_full"] == pytest.approx(raw["pred_rms_core"]
                                                                 * np.exp(raw["scope"][h]["pred_log_scope"]))


def test_the_full_trim_reaches_the_raw_core_forecast_with_one_factor(bt, sp, monkeypatch, own_repo, tmp_path):
    from vec import harness, t2_build
    from vec.io import COORD_KEY

    path, size, clouds = _size_json(sp, tmp_path)
    assert size["valid"] is True and size["preregistered"]["z_mode"] == "xy_grow_z_hold"
    real_size = t2_build.size_or_hold
    install(bt, monkeypatch.setattr, [])
    monkeypatch.setattr(bt.t2_build, "size_or_hold", real_size)
    monkeypatch.setitem(bt.SIZES, "heart", path)
    ref = synthetic_stage(3000, 2, 0.0, 60)
    ref.obsm["spatial_3D"] = clouds[-1][2].astype(np.float32)
    ref.obs["core"] = clouds[-1][3]
    start = synthetic_stage(*STAGES["E9.5.h5ad"])
    load = bt.load_stage
    monkeypatch.setattr(bt, "load_stage", lambda p, **kw: ref.copy() if Path(p).name == "E10.5.h5ad"
                        else load(p, **kw))

    out, rec, _ = _build(bt, tmp_path, "core", *CORE_TRIM)
    D = np.asarray(harness.draw_copy(_ctx(_trimmed(ref, 0.0), start))[3], dtype=np.float32).astype(float)
    G = np.asarray(ad.read_h5ad(out).obsm[COORD_KEY], dtype=float)[:, :3]
    np.testing.assert_allclose(G.std(0) / G.std(0)[0], D.std(0) / D.std(0)[0], rtol=1e-5)
    target = size["z_modes"]["raw"]["pred_rms_core"]
    assert target != pytest.approx(size["preregistered"]["pred_rms_core"])  # the raw mode's, not the recorded one's
    assert t2_build.rms_radius(G) == pytest.approx(target, rel=1e-5)
    s = rec["size"]
    assert s["target_rms_full"] == pytest.approx(target) and s["target_source"] == "z_modes.raw.pred_rms_core"
    assert s["hypothesis"] == "core" and s["z_mode_used"] == "raw" and s["z_held"] is False
    assert "scale_vs_reference_full" not in s


@pytest.mark.parametrize("z_mode", [None, "xy_grow_z_hold", "common_step"])
def test_the_core_scope_is_sized_in_the_raw_mode_only(sp, tmp_path, z_mode):
    from vec import t2_build

    path, _, clouds = _size_json(sp, tmp_path)
    with pytest.raises(ValueError, match="core's own scope is sized in the"):
        t2_build.apply_size_procedure(clouds[-1][2], path, "core", z_mode=z_mode, expect_board="heart")


def test_the_core_scope_gets_no_fallback_size(sp, tmp_path):
    """An invalid size JSON: approved rule 3's fallback keeps E10.5's size for H1, H0 and H2 only; the core's own scope
    is refused, its size needing a ruling."""
    from vec import t2_build

    path, size, clouds = _size_json(sp, tmp_path)
    size["valid"] = False
    path.write_text(sp.dumps(size))
    with pytest.raises(ValueError, match="needs the captain's ruling"):
        t2_build.size_or_hold(clouds[-1][2], path, "core", source_stage="E10.5", reference_cloud=clouds[-1][2],
                              expect_board="heart")


def test_the_scope_hypotheses_size_as_before(sp, tmp_path):
    """H1 and H2 in every applicable mode: the same target and the same record keys as the pre-registered path gives
    (scale_vs_reference_full kept, no target_source)."""
    from vec import t2_build

    path, size, clouds = _size_json(sp, tmp_path)
    for mode in (None, "raw", "common_step", "xy_grow_z_hold"):
        for h in ("H1", "H2"):
            _, rec = t2_build.apply_size_procedure(clouds[-1][2], path, h, z_mode=mode, expect_board="heart")
            res = size["preregistered"] if mode in (None, "xy_grow_z_hold") else size["z_modes"][mode]
            assert rec["target_rms_full"] == pytest.approx(res["scope"][h]["pred_rms_full"])
            assert "target_source" not in rec and "scale_vs_reference_full" in rec


# ------------------------------------------------------------------------------------------------ refusals
@pytest.mark.parametrize("extra", [
    ["--entry", "core_trim"], ["--entry", "core_trim", "--z-mode", "raw"],
    ["--entry", "core_trim", "--expression", "copy"],
    ["--entry", "core_trim", "--expression", "copy", "--z-mode", "recorded"],
    ["--entry", "core_trim", "--trend"], ["--entry", "core_trim", "--depth"],
    ["--entry", "core_trim", "--herd"], ["--entry", "core_trim", "--cfm"],
    [*CORE_TRIM, "--trend"], [*CORE_TRIM, "--depth"], [*CORE_TRIM, "--herd"], [*CORE_TRIM, "--cfm"]])
def test_other_core_trim_combinations_are_refused_before_reading_anything(bt, heart_env, tmp_path, extra):
    with pytest.raises(SystemExit) as e:
        bt.main(["heart", *extra, "--out", str(tmp_path / "h.h5ad")])
    assert "refused" in str(e.value)
    assert heart_env[0] == [] and not (tmp_path / "h.h5ad").exists()


def test_the_core_trim_refusal_names_the_one_command(bt, heart_env, tmp_path):
    with pytest.raises(SystemExit) as e:
        bt.main(["heart", "--entry", "core_trim", "--out", str(tmp_path / "h.h5ad")])
    assert bt.HEART_CORE_TRIM in str(e.value)


def test_the_record_round_trips_as_json(bt, heart_env, tmp_path):
    _, rec, text = _build(bt, tmp_path, "core", *CORE_TRIM)
    assert json.loads(text) == rec
