"""`experiments/builds/build_test.py heart --entry organ_scope --expression copy --z-mode raw` (the heart's organ cut,
file D of the heart's order rule; row R21 of experiments/harness/TRIES.md, the captain's order of 2026-11-22 to code and
read it). On synthetic stages with released-style labels, a synthetic test board whose limit is never shown and
synthetic size JSONs: the option is off by default and refuses every other combination before reading anything; the
file is the copy draw from E10.5 under vec.t2_build.organ_scope_mask; it is sized to the raw-mode core forecast times
rms(organ pool) / rms(core) on E10.5, and gets no fallback; an unmapped label refuses the build; the record names the
scope, the label config's sha256, L, the pool's E10.5 rows and the size field, and no test-board count; the half
follows. Every heart output that existed before the option keeps its bytes (the eight commands pinned in
tests/test_build_heart_iso_trim.py, and here the half and full trims' copy files). Nothing is read from data/ and
nothing is scored."""

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
from test_build_heart_core_trim import _size_json  # noqa: E402
from test_build_heart_iso_trim import _build, _ctx, canon  # noqa: E402
from test_heart_trend_build import (  # noqa: E402
    SECRET_MAX,
    STAGES,
    install,
    sha,
    synthetic_stage,
)

ORGAN = ["--entry", "organ_scope", "--expression", "copy", "--z-mode", "raw"]
# sha256 of (the .h5ad, its .build.json, the record's `canon`) that build_test.py heart wrote for the half and full
# trims' copy files on the synthetic stages with the size stand-in of `install`, computed at main e395e37 (before this
# option) by the same builds as below. The other eight heart commands are pinned in test_build_heart_iso_trim.py.
BEFORE = {
    "s2_trim_copy_raw": ("6287272082da44e01a11ed5f705e78852e355dd655ec5ba257be55d31db67bad",
                         "7e62d6013348cb4c109d862944c883312f9e45f5c03c39140e112ca0e095d7c7",
                         "e59e7dc7446ba99a50e98980ff350ed56e90903bb0e401f9ea46ba189bc8e849"),
    "core_trim_copy_raw": ("73e57af3bc859644198fc5061fbfc200fcd3742605a79838cf03fc113819f254",
                           "cf0e1d0a81cf38c8e5fa542947cd40fdebe0a5d8b4572f0d006d9be3e77bd52b",
                           "53b26a922696642e7c0098d595ed5d712d1f6e738eab6da7afdd70a884c22c55"),
}
ARGV = {"s2_trim_copy_raw": ["--entry", "s2_trim", "--expression", "copy", "--z-mode", "raw"],
        "core_trim_copy_raw": ["--entry", "core_trim", "--expression", "copy", "--z-mode", "raw"]}
ISLAND_SHIFT = 80.0  # synthetic: vessels moved this far from the core form a detached island


@pytest.fixture(scope="module")
def bt():
    sys.path.insert(0, str(REPO / "experiments" / "builds"))
    import build_test

    return build_test


@pytest.fixture(scope="module")
def sp():
    """experiments/coords/size_procedure.py, loaded by path (as tests/test_t2_build.py does)."""
    import importlib.util

    spec = importlib.util.spec_from_file_location("coords_size_procedure_organ",
                                                  REPO / "experiments" / "coords" / "size_procedure.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def organ_stage() -> ad.AnnData:
    """The synthetic E10.5 of `install` with released-style labels: CM (the core) as V-CM, half the Endothelial cells
    as BEC beside the core and a sixth of them moved far away (a detached island), the rest as Dorsal-PAM (not
    chamber), Untracked as Late-a-PHM and the tiny lineage as Epi."""
    a = synthetic_stage(*STAGES["E10.5.h5ad"])
    merged = np.asarray(a.obs["merged"])
    labels = np.full(a.n_obs, "Late-a-PHM", dtype=object)
    labels[merged == "CM"] = "V-CM"
    labels[merged == "Tiny"] = "Epi"
    endo = np.flatnonzero(merged == "Endothelial")
    labels[endo[: len(endo) // 2]] = "BEC"
    labels[endo[len(endo) // 2:]] = "Dorsal-PAM"
    island = endo[: len(endo) // 6]
    C = np.asarray(a.obsm["spatial_3D"], dtype=np.float32).copy()
    C[island, 0] += ISLAND_SHIFT
    a.obsm["spatial_3D"] = C
    a.obs["celltype"] = labels.astype(str)
    a.obs["core"] = merged == "CM"
    return a


def _organ_mask(ref: ad.AnnData) -> np.ndarray:
    from vec import t2_build

    return t2_build.organ_scope_mask(np.asarray(ref.obsm["spatial_3D"], dtype=float)[:, :3],
                                     np.asarray(ref.obs["celltype"]), np.asarray(ref.obs["core"], dtype=bool))


def _use_stage(bt, monkeypatch, ref: ad.AnnData, reads: list) -> None:
    load = bt.load_stage

    def wrapped(p, **kw):
        if Path(p).name == "E10.5.h5ad":
            reads.append("E10.5.h5ad")
            return ref.copy()
        return load(p, **kw)

    monkeypatch.setattr(bt, "load_stage", wrapped)


@pytest.fixture
def organ_env(bt, monkeypatch, own_repo):
    """install()'s synthetic stages and board, E10.5 relabelled (organ_stage), and a size stand-in that scales by
    1.25 and reports its target as the core forecast (the real size step is tested below)."""
    from vec import t2_build

    reads: list[str] = []
    install(bt, monkeypatch.setattr, reads)
    _use_stage(bt, monkeypatch, organ_stage(), reads)
    calls: list[dict] = []

    def size(coords, size_json, hypothesis=None, **kw):
        calls.append({"hypothesis": hypothesis, **kw})
        C = np.asarray(coords, dtype=float)
        return C * 1.25, {"hypothesis": hypothesis, "factor": 1.25,
                          "target_rms_full": t2_build.rms_radius(C) * 1.25,
                          "target_source": "z_modes.raw.pred_rms_core", "warnings": []}

    monkeypatch.setattr(bt.t2_build, "size_or_hold", size)
    return reads, calls


# ------------------------------------------------------------------------------------------------ nothing else moves
def test_the_option_is_off_by_default(bt):
    a = bt.build_parser().parse_args(["heart", "--entry", "organ_scope", "--out", "runs/x/h.h5ad"])
    assert a.expression == "b0" and a.z_mode == "recorded" and not a.trend
    assert bt.HEART_ORGAN_SCOPE == "heart " + " ".join(ORGAN)
    assert set(bt.HEART_SCOPE_OPTIONS) == {"s2_trim", "core_trim"}  # the trims' table is untouched


@pytest.mark.parametrize("name", sorted(BEFORE))
def test_the_trims_copy_files_keep_their_bytes_and_record(bt, own_repo, monkeypatch, tmp_path, name):
    """The half and full trims' copy commands write the same .h5ad and build record, byte for byte, as main at e395e37
    wrote them on the same synthetic stages; a record may differ from that text only in the last bits of its floats."""
    install(bt, monkeypatch.setattr, [])
    out, rec, text = _build(bt, tmp_path, name, *ARGV[name])
    assert sha(out) == BEFORE[name][0]
    assert hashlib.sha256(text.encode()).hexdigest() == BEFORE[name][1] or canon(text) == BEFORE[name][2]
    assert "organ" not in text


# ------------------------------------------------------------------------------------------------ the new file
def test_the_file_is_the_copy_draw_from_the_organ_pool(bt, organ_env, tmp_path, capsys):
    from vec import harness, t2_build
    from vec.io import COORD_KEY, dense

    reads, calls = organ_env
    out, rec, text = _build(bt, tmp_path, "organ", *ORGAN)
    ref, start = organ_stage(), synthetic_stage(*STAGES["E9.5.h5ad"])
    keep = _organ_mask(ref)
    island = np.flatnonzero(np.asarray(ref.obsm["spatial_3D"])[:, 0] > ISLAND_SHIFT / 2)
    labels = np.asarray(ref.obs["celltype"])
    assert len(island) and not keep[island].any()  # the detached vessels are cut
    assert not keep[np.isin(labels, ["Dorsal-PAM", "Late-a-PHM"])].any() and keep[np.asarray(ref.obs["core"])].all()
    pool = ref[keep].copy()
    _, rows, X, coords = harness.draw_copy(_ctx(pool, start))
    got = ad.read_h5ad(out)
    np.testing.assert_array_equal(dense(got.X).astype(np.float32), X)
    C = np.asarray(ref.obsm["spatial_3D"], dtype=float)
    ratio = t2_build.organ_core_rms_ratio(C, keep, np.asarray(ref.obs["core"], dtype=bool))
    D = np.asarray(coords, dtype=np.float32).astype(float)
    want = t2_build.rescale_to_rms(D, t2_build.rms_radius(D) * 1.25 * ratio)
    np.testing.assert_allclose(np.asarray(got.obsm[COORD_KEY], dtype=float)[:, :3], want, rtol=1e-5, atol=1e-5)
    assert reads == ["E10.5.h5ad", "E9.5.h5ad"]
    assert len(calls) == 1 and calls[0]["hypothesis"] == "core" and calls[0]["z_mode"] == "raw"
    assert calls[0]["source_stage"] == "E10.5" and calls[0]["reference_cloud"].n_obs == ref.n_obs
    assert got.n_obs == SECRET_MAX  # the count rule does not read the pool
    assert rec["n_cells"] == "max_cells" and str(SECRET_MAX) not in text + capsys.readouterr().out


def test_the_record_names_scope_config_l_pool_rows_and_size_field(bt, organ_env, tmp_path):
    from vec import t2_build

    _, rec, text = _build(bt, tmp_path, "organ", *ORGAN)
    ref = organ_stage()
    keep = _organ_mask(ref)
    C, core = np.asarray(ref.obsm["spatial_3D"], dtype=float), np.asarray(ref.obs["core"], dtype=bool)
    spec = t2_build.load_t2_chamber_spec()
    opt, trim, size = rec["scope_option"], rec["trim"], rec["size"]
    assert opt["command"] == "heart " + " ".join(ORGAN) and "R21" in opt["ruling"] and opt["scope"] == "organ_scope"
    assert opt["label_config"] == "configs/t2_chamber_classes.yaml" and opt["label_config_sha256"] == spec.sha256
    assert opt["link_percentile"] == 99.0
    assert opt["link_length"] == pytest.approx(t2_build.organ_link_length(C, core))
    assert opt["kept_share_field"] == "trim.kept_share" and opt["pool_rows_field"] == "trim.pool_rows"
    assert trim["pool_rows"] == int(keep.sum()) and trim["kept_share"] == pytest.approx(keep.mean())
    assert opt["size_hypothesis"] == "organ" and opt["size_source_field"] == t2_build.SIZE_ORGAN_SOURCE
    assert "halve.py" in opt["halving"]
    assert rec["entry"] == "copy organ_scope (organ)" and rec["expression"] == "copy" and rec["z_mode"] == "raw"
    assert size["hypothesis"] == "organ" and size["target_source"] == t2_build.SIZE_ORGAN_SOURCE
    assert size["organ_core_rms_ratio"] == pytest.approx(t2_build.organ_core_rms_ratio(C, keep, core))
    assert trim["organ_core_rms_ratio"] == size["organ_core_rms_ratio"]
    assert size["core_target_source"] == "z_modes.raw.pred_rms_core" and size["ratio_stage"] == "E10.5"
    assert json.loads(text) == rec


def test_an_unmapped_label_refuses_the_build(bt, organ_env, monkeypatch, tmp_path):
    ref = organ_stage()
    labels = np.asarray(ref.obs["celltype"]).copy()
    labels[5] = "Mystery"
    ref.obs["celltype"] = labels
    _use_stage(bt, monkeypatch, ref, [])
    with pytest.raises(SystemExit) as e:
        bt.main(["heart", *ORGAN, "--out", str(tmp_path / "h.h5ad")])
    assert "Mystery" in str(e.value) and "refused" in str(e.value)
    assert not (tmp_path / "h.h5ad").exists()


def test_the_half_of_the_organ_file(bt, organ_env, tmp_path):
    sys.path.insert(0, str(REPO / "experiments" / "builds"))
    import halve

    out, _, _ = _build(bt, tmp_path, "organ", *ORGAN)
    half = tmp_path / "half.h5ad"
    digest = halve.halve(str(out), str(half))
    full, got = ad.read_h5ad(out), ad.read_h5ad(half)
    assert halve.DEFAULT_SEED == 1_000_001
    assert got.n_obs == full.n_obs // 2 and digest == sha(half)
    assert list(got.obs_names) == list(full.obs_names[halve.half_rows(full.n_obs, halve.DEFAULT_SEED)])


def test_two_builds_give_the_same_bytes(bt, organ_env, tmp_path):
    a, _, ta = _build(bt, tmp_path, "one", *ORGAN)
    b, _, tb = _build(bt, tmp_path, "two", *ORGAN)
    assert sha(a) == sha(b) and ta == tb


# ------------------------------------------------------------------------------------------------ the real size step
def _real_size_env(bt, sp, monkeypatch, tmp_path, **over):
    from vec import t2_build

    path, size, clouds = _size_json(sp, tmp_path, **over)
    real_size = t2_build.size_or_hold
    install(bt, monkeypatch.setattr, [])
    monkeypatch.setattr(bt.t2_build, "size_or_hold", real_size)
    monkeypatch.setitem(bt.SIZES, "heart", path)
    ref = synthetic_stage(3000, 2, 0.0, 60)
    C, core = clouds[-1][2], clouds[-1][3]
    ref.obsm["spatial_3D"] = C.astype(np.float32)
    ref.obs["core"] = core
    shell = np.flatnonzero(~core)
    labels = np.where(core, "V-CM", "Dorsal-PAM").astype(object)
    labels[shell[::2]] = "BEC"
    ref.obs["celltype"] = labels.astype(str)
    _use_stage(bt, monkeypatch, ref, [])
    return path, size, ref


def test_the_organ_file_reaches_the_raw_core_forecast_times_the_e10_5_ratio(bt, sp, monkeypatch, own_repo, tmp_path):
    from vec import harness, t2_build
    from vec.io import COORD_KEY

    _, size, ref = _real_size_env(bt, sp, monkeypatch, tmp_path)
    assert size["valid"] is True and size["preregistered"]["z_mode"] == "xy_grow_z_hold"
    out, rec, _ = _build(bt, tmp_path, "organ", *ORGAN)
    keep = _organ_mask(ref)
    C, core = np.asarray(ref.obsm["spatial_3D"], dtype=float), np.asarray(ref.obs["core"], dtype=bool)
    ratio = t2_build.organ_core_rms_ratio(C, keep, core)
    assert ratio != pytest.approx(1.0) and 0 < keep.sum() < ref.n_obs
    start = synthetic_stage(*STAGES["E9.5.h5ad"])
    D = np.asarray(harness.draw_copy(_ctx(ref[keep].copy(), start))[3], dtype=np.float32).astype(float)
    G = np.asarray(ad.read_h5ad(out).obsm[COORD_KEY], dtype=float)[:, :3]
    np.testing.assert_allclose(G.std(0) / G.std(0)[0], D.std(0) / D.std(0)[0], rtol=1e-5)  # one factor
    core_target = size["z_modes"]["raw"]["pred_rms_core"]
    assert core_target != pytest.approx(size["preregistered"]["pred_rms_core"])
    assert t2_build.rms_radius(G) == pytest.approx(core_target * ratio, rel=1e-5)
    s = rec["size"]
    assert s["core_target_rms"] == pytest.approx(core_target) and s["core_target_source"] == "z_modes.raw.pred_rms_core"
    assert s["target_rms_full"] == pytest.approx(core_target * ratio) and s["hypothesis"] == "organ"
    assert s["z_mode_used"] == "raw" and s["z_held"] is False and "scale_vs_reference_full" not in s


def test_the_organ_size_gets_no_fallback(bt, sp, monkeypatch, own_repo, tmp_path):
    from vec import t2_build

    path, size, ref = _real_size_env(bt, sp, monkeypatch, tmp_path, valid=False)
    with pytest.raises(ValueError, match="needs the captain's ruling"):
        t2_build.size_organ_scope(np.asarray(ref.obsm["spatial_3D"], dtype=float), path, 1.5, source_stage="E10.5",
                                  reference_cloud=ref, expect_board="heart")
    with pytest.raises(ValueError):
        bt.main(["heart", *ORGAN, "--out", str(tmp_path / "h.h5ad")])
    assert not (tmp_path / "h.h5ad").exists()


def test_the_organ_size_refuses_cells_from_another_stage(sp, tmp_path):
    from vec import t2_build

    path, _, clouds = _size_json(sp, tmp_path)
    with pytest.raises(ValueError, match="not E10.5"):
        t2_build.size_organ_scope(clouds[-1][2], path, 1.5, source_stage="E9.5", reference_cloud=clouds[-1][2])


# ------------------------------------------------------------------------------------------------ refusals
@pytest.mark.parametrize("extra", [
    ["--entry", "organ_scope"], ["--entry", "organ_scope", "--z-mode", "raw"],
    ["--entry", "organ_scope", "--expression", "copy"],
    ["--entry", "organ_scope", "--expression", "copy", "--z-mode", "recorded"],
    ["--entry", "organ_scope", "--trend"], ["--entry", "organ_scope", "--depth"],
    ["--entry", "organ_scope", "--herd"], ["--entry", "organ_scope", "--cfm"],
    [*ORGAN, "--trend"], [*ORGAN, "--depth"], [*ORGAN, "--herd"], [*ORGAN, "--cfm"]])
def test_other_organ_combinations_are_refused_before_reading_anything(bt, organ_env, tmp_path, extra):
    with pytest.raises(SystemExit) as e:
        bt.main(["heart", *extra, "--out", str(tmp_path / "h.h5ad")])
    assert "refused" in str(e.value) and bt.HEART_ORGAN_SCOPE in str(e.value)
    assert organ_env[0] == [] and not (tmp_path / "h.h5ad").exists()
