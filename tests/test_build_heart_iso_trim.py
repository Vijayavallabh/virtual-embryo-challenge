"""`experiments/builds/build_test.py heart --entry s2_trim --expression copy --z-mode raw` (the heart's scope bet,
branch campaign/heart-iso-trim; no ruling): the ruled second file's two fixes (the copy's expression, the size step's
raw mode) on the half trim, H2 of docs/S2_RULES.md. On synthetic stages, a synthetic test board whose limit is never
shown and a synthetic size JSON: every heart output that existed before the option keeps its bytes and its record; the
new file is the copy draw from the trimmed pool with its released expression, scaled by one factor to the raw mode's H2
target; the count does not read the pool; the half follows; the refusals. Nothing is read from data/ and nothing is
scored."""

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
from test_build_heart_floor_iso import _stage_cloud  # noqa: E402
from test_heart_trend_build import (  # noqa: E402
    GENES,
    SECRET_MAX,
    STAGES,
    heart_test_board,
    install,
    sha,
    synthetic_stage,
)

SEED = 1_000_001
# sha256 of (the .h5ad, its .build.json) that build_test.py heart wrote on the synthetic stages with the size stand-in
# of `install`, computed at main a4d58a6 (before this option) by the same builds as below
BEFORE = {
    "s1": ("cf887247be3d77a3fbbb6a99848a080137481f45a59bb97e4be3b3187e02e31e",
           "7772c7ccfcee8d33914afc0bc38ef25e38f861ca32505d78274176fd268bcd32"),
    "s2_trim": ("fc37df8923c87f39e6432eec45a2d8d49eae7c4e3b572fa057588bce949d5b9a",
                "167e4e407fb0079238c994b6e9c3270115898f2c7ee34b7e549c9909480ae6d1"),
    "s1_trend": ("187673f9dabbe90d49d0c34af581fcb1431b6bfe2dc6fe01ff9f55d95698b0b4",
                 "0f3304eaced7e5cea546fefc0ba6c4ff87679a184d85a4d694589f046a889d2c"),
    "s2_trim_trend": ("4014fb74785059eea40fa7c6e9d7b22be05a3747af1162c1aa9fd1f26b24643b",
                      "dad93fdbfc321a6cf82d7941cc8a8607bea28084d0308e0c62cf80d03b6c6d90"),
    "s1_copy": ("606ccc9479eb9f181cd333931de54c949150e732229f8394a1cf83785e988f3e",
                "4b3adaa11240e4acaac035d764e4e1e1973f1af8d72fba8385a1541c353141e7"),
    "s1_copy_raw": ("606ccc9479eb9f181cd333931de54c949150e732229f8394a1cf83785e988f3e",
                    "ee4c9f8b97e51ec839dffb499efe2da8fbd9a8729373c4f38a4069f24e574b57"),
    "s2_trim_raw": ("fc37df8923c87f39e6432eec45a2d8d49eae7c4e3b572fa057588bce949d5b9a",
                    "7c2a496054102ad577fa2acb9d664f7f628196779ef82c3fa338933751c05a1e"),
    "s2_raw": ("cf887247be3d77a3fbbb6a99848a080137481f45a59bb97e4be3b3187e02e31e",
               "bc3ab44995e7a5b0091e56e07ca9b87a3439a3f23f70680107badc948b9a77a1"),
}

# The same records with every float rounded to 12 significant digits (`canon`). The trend block's summaries (cosine,
# norm ratio, mean log factor) differ in their last bits between machines: relative 1e-16 between the local machine,
# where BEFORE was computed, and the DGX on 2026-11-19, where the exact text hash of two trend records failed while
# their .h5ad hashes matched. Computed from the records whose exact hashes are BEFORE's.
BEFORE_CANON = {
    "s1": "2b0674350ca8655ff0ab235b8a7c8284d5a42eccaae165db62461575c32e27ef",
    "s2_trim": "40accedc96b1df4335b4d3fb1af8c9a12d8076df28d923ed27da57729a419728",
    "s1_trend": "9b6d79152a5e3adec4c059ebd17f3b08947fb42bada7a1cf1a7c7e96b2779497",
    "s2_trim_trend": "788ecb4d87720fb113b8bec85653a1489051500cb32e4643132a1621ad4361ec",
    "s1_copy": "961e7625fc56c426b85f566b1114e3dad9997c61e1b3528bd18d8b082b1b20cd",
    "s1_copy_raw": "b41e33a1b0a7491e79b0bd4d40ef6364509482b5c37cb242289a65f4a1da1ed8",
    "s2_trim_raw": "64cbe46998344ecbccb1520abccff88ccd79324efb1c57a0c5e9f4e462c1cea9",
    "s2_raw": "9129317d538bfa53d58401c54c563e80c21c80fcd3da7651edac2d529c8589fe",
}


def canon(text: str) -> str:
    def rounded(o):
        if isinstance(o, float):
            return float(f"{o:.12g}")
        if isinstance(o, dict):
            return {k: rounded(v) for k, v in o.items()}
        if isinstance(o, list):
            return [rounded(v) for v in o]
        return o

    return hashlib.sha256(json.dumps(rounded(json.loads(text)), sort_keys=True).encode()).hexdigest()
ARGV = {"s1": ["--entry", "s1"], "s2_trim": ["--entry", "s2_trim"],
        "s1_trend": ["--entry", "s1", "--trend"], "s2_trim_trend": ["--entry", "s2_trim", "--trend"],
        "s1_copy": ["--entry", "s1", "--expression", "copy"],
        "s1_copy_raw": ["--entry", "s1", "--expression", "copy", "--z-mode", "raw"],
        "s2_trim_raw": ["--entry", "s2_trim", "--z-mode", "raw"],
        "s2_raw": ["--entry", "s2", "--z-mode", "raw"]}
ISO_TRIM = ["--entry", "s2_trim", "--expression", "copy", "--z-mode", "raw"]


@pytest.fixture(scope="module")
def bt():
    sys.path.insert(0, str(REPO / "experiments" / "builds"))
    import build_test

    return build_test


@pytest.fixture(scope="module")
def sp():
    """experiments/coords/size_procedure.py, loaded by path (as tests/test_t2_build.py does)."""
    import importlib.util

    spec = importlib.util.spec_from_file_location("coords_size_procedure_iso_trim",
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


def _trimmed(ref: ad.AnnData) -> ad.AnnData:
    from vec import t2_build
    from vec.io import COORD_KEY

    keep = t2_build.radial_trim_mask(np.asarray(ref.obsm[COORD_KEY], dtype=float)[:, :3],
                                     np.asarray(ref.obs["core"], dtype=bool), 0.5)
    assert 0 < keep.sum() < ref.n_obs
    return ref[keep].copy()


# ------------------------------------------------------------------------------------------------ nothing else moves
def test_the_option_is_off_by_default(bt):
    a = bt.build_parser().parse_args(["heart", "--entry", "s2_trim", "--out", "runs/x/h.h5ad"])
    assert a.expression == "b0" and a.z_mode == "recorded"
    assert bt.HEART_ISO_TRIM == "heart " + " ".join(ISO_TRIM)


@pytest.mark.parametrize("name", sorted(BEFORE))
def test_every_heart_output_before_the_option_keeps_its_bytes_and_record(bt, heart_env, tmp_path, name):
    """Each heart command that built before the option writes the same .h5ad and the same build record, byte for byte,
    as main at a4d58a6 wrote them on the same synthetic stages; a record may differ from that text only in the last
    bits of its floats (`canon`)."""
    out, _, text = _build(bt, tmp_path, name, *ARGV[name])
    assert sha(out) == BEFORE[name][0]
    assert hashlib.sha256(text.encode()).hexdigest() == BEFORE[name][1] or canon(text) == BEFORE_CANON[name]


# ------------------------------------------------------------------------------------------------ the new file
def test_iso_trim_is_the_copy_draw_from_the_trimmed_pool(bt, heart_env, tmp_path, capsys):
    """The file is the copy baseline's draw (`vec.harness.draw_copy`, same seed and count) from the half-trimmed E10.5
    pool, with the released expression of the drawn cells and no momentum; the size step gets H2 and the raw mode, with
    the full E10.5 cloud as its reference; the record says copy, s2_trim, H2, the raw mode and the trim, and holds no
    cell number."""
    from vec import harness
    from vec.io import COORD_KEY, dense

    reads, calls = heart_env
    out, rec, text = _build(bt, tmp_path, "iso_trim", *ISO_TRIM)
    stages = {f: synthetic_stage(*v) for f, v in STAGES.items()}
    ref, start = stages["E10.5.h5ad"], stages["E9.5.h5ad"]
    pool = _trimmed(ref)
    _, rows, X, coords = harness.draw_copy(_ctx(pool, start))
    got = ad.read_h5ad(out)
    Xg = dense(got.X).astype(np.float32)
    np.testing.assert_array_equal(Xg, X)
    np.testing.assert_array_equal(Xg, dense(pool[:, GENES].X[rows]).astype(np.float32))  # the released values
    np.testing.assert_array_equal(np.asarray(got.obsm[COORD_KEY], dtype=np.float32)[:, :3],
                                  (np.asarray(coords, dtype=np.float32).astype(float) * 1.25).astype(np.float32))
    assert set(np.asarray(pool.obs_names)[rows]) <= set(pool.obs_names)  # only trimmed-pool cells
    assert reads == ["E10.5.h5ad", "E9.5.h5ad"]
    assert len(calls) == 1 and calls[0]["hypothesis"] == "H2" and calls[0]["z_mode"] == "raw"
    assert calls[0]["source_stage"] == "E10.5"
    assert calls[0]["reference_cloud"].n_obs == ref.n_obs  # the full E10.5 cloud, not the trimmed pool
    assert rec["entry"] == "copy s2_trim (H2)" and rec["expression"] == "copy" and rec["z_mode"] == "raw"
    assert rec["trim"]["kept_share"] == pytest.approx(pool.n_obs / ref.n_obs)
    assert rec["trim"]["xy_scope_after"] <= 0.5 * rec["trim"]["xy_scope_before"] + 1e-12
    assert rec["n_cells"] == "max_cells" and "trend" not in rec and rec["depth_step"] is False
    assert str(SECRET_MAX) not in text + capsys.readouterr().out


def test_iso_trim_differs_from_the_ruled_file_by_scope_only(bt, heart_env, tmp_path):
    """Against the ruled second file (s1 copy raw): the same count, a pool that is a strict subset of E10.5, and the
    hypothesis H2 in place of H1 in the size call."""
    a, rec_a, _ = _build(bt, tmp_path, "iso", "--entry", "s1", "--expression", "copy", "--z-mode", "raw")
    b, rec_b, _ = _build(bt, tmp_path, "iso_trim", *ISO_TRIM)
    ga, gb = ad.read_h5ad(a), ad.read_h5ad(b)
    assert ga.n_obs == gb.n_obs == SECRET_MAX
    assert [c["hypothesis"] for c in heart_env[1]] == ["H1", "H2"]
    assert rec_a["trim"] is None and rec_b["trim"] is not None
    assert sha(a) != sha(b)


def test_the_count_does_not_read_the_pool():
    """The count rule (vec.t2_build.n_cells, through harness.count_for) gives the same count for any pool size: the
    trimmed pool gets s1's count, which is why the half follows as for the ruled file."""
    from vec import harness

    board = heart_test_board()
    counts = {harness.count_for(board, n) for n in (50, 400, SECRET_MAX, 2600, 10**6)}
    assert counts == {SECRET_MAX}


def test_the_half_of_the_iso_trim_file(bt, heart_env, tmp_path):
    """halve.py on the new file: a uniform random half of its rows (seed 1000001, order kept), nothing else changed."""
    sys.path.insert(0, str(REPO / "experiments" / "builds"))
    import halve

    out, _, _ = _build(bt, tmp_path, "iso_trim", *ISO_TRIM)
    half = tmp_path / "half.h5ad"
    digest = halve.halve(str(out), str(half))
    full, got = ad.read_h5ad(out), ad.read_h5ad(half)
    idx = halve.half_rows(full.n_obs, halve.DEFAULT_SEED)
    assert got.n_obs == full.n_obs // 2 and digest == sha(half)
    assert list(got.obs_names) == list(full.obs_names[idx])
    np.testing.assert_array_equal(np.asarray(got.X.todense() if hasattr(got.X, "todense") else got.X),
                                  np.asarray((full.X.todense() if hasattr(full.X, "todense") else full.X)[idx]))


def test_iso_trim_reaches_the_raw_mode_h2_target_with_one_factor(bt, sp, monkeypatch, own_repo, tmp_path):
    """The real size step on a synthetic size JSON that recorded xy_grow_z_hold: the trimmed copy draw is scaled about
    its centroid by one factor (per-axis spread ratios kept) to the raw mode's own H2 forecast (z_modes.raw.scope.H2),
    not to the recorded mode's H2 target and not to H1's."""
    from vec import harness, t2_build
    from vec.io import COORD_KEY

    sd = (100.0, 80.0, 30.0)
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
    raw = size["z_modes"]["raw"]["scope"]
    assert raw["H2"]["pred_rms_full"] != pytest.approx(raw["H1"]["pred_rms_full"])
    path = tmp_path / "size_heart.json"
    path.write_text(sp.dumps(size))

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

    out, rec, _ = _build(bt, tmp_path, "iso_trim", *ISO_TRIM)
    _, coords_draw = harness.draw_copy(_ctx(_trimmed(ref), start))[2:]
    D = np.asarray(coords_draw, dtype=np.float32).astype(float)
    G = np.asarray(ad.read_h5ad(out).obsm[COORD_KEY], dtype=float)[:, :3]
    np.testing.assert_allclose(G.std(0) / G.std(0)[0], D.std(0) / D.std(0)[0], rtol=1e-5)
    target = raw["H2"]["pred_rms_full"]
    assert t2_build.rms_radius(G) == pytest.approx(target, rel=1e-5)
    assert rec["size"]["target_rms_full"] == pytest.approx(target) and rec["size"]["z_mode_used"] == "raw"
    assert rec["size"]["hypothesis"] == "H2" and rec["size"]["z_held"] is False


# ------------------------------------------------------------------------------------------------ refusals
@pytest.mark.parametrize("extra", [["--entry", "s2_trim"], ["--entry", "s2_trim", "--z-mode", "recorded"],
                                   ["--entry", "s2", "--z-mode", "raw"],
                                   ["--entry", "s2_trim", "--z-mode", "raw", "--trend"],
                                   ["--entry", "s2_trim", "--z-mode", "raw", "--depth"],
                                   ["--entry", "s2_trim", "--z-mode", "raw", "--herd"],
                                   ["--entry", "s2_trim", "--z-mode", "raw", "--cfm"]])
def test_other_copy_combinations_are_refused_before_reading_anything(bt, heart_env, tmp_path, extra):
    with pytest.raises(SystemExit) as e:
        bt.main(["heart", *extra, "--expression", "copy", "--out", str(tmp_path / "h.h5ad")])
    assert "--expression copy is refused" in str(e.value)
    assert heart_env[0] == [] and not (tmp_path / "h.h5ad").exists()
