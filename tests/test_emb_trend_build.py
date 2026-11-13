"""`experiments/builds/build_test.py embryo --expression trend` (port/emb-s1-trend: fgw_trend_render as a first-upload
option, no entry until the captain rules), on synthetic stages and a synthetic test board whose limit is never shown:
the option is accepted and every existing embryo build is byte for byte what main's builder wrote; with every gene
common the builder's expression is the harness candidate's; a sample lacking two panel genes leaves those two columns
as fgw has them and renders the rest exactly as on the panel without them, and the record says 2; the refusals; a
sample with no shared extraembryonic label still raises the existing error. Nothing is read from data/, nothing is
scored, and no prediction is compared with any target."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from vec.boards import Board

REPO = Path(__file__).resolve().parents[1]
BUILDS = REPO / "experiments" / "builds"
PANEL = [f"g{i}" for i in range(24)]
OUTSIDE = [f"x{i}" for i in range(3)]
LACKING = ("g3", "g9")  # the two panel genes the synthetic E7.25 lacks in the missing-gene case
SECRET_MAX = 437  # stands for the test board's max_cells: never printed, never in a record
# label -> (lineage group, centre of its cells): one extraembryonic label every stage holds, two embryonic ones
LABELS = {"EXE-Endoderm": ("Endoderm", (0.0, 0.0, 0.0)), "SOM": ("ParaxialMesoderm", (6.0, 0.0, 0.0)),
          "PS": ("Mesoderm", (0.0, 6.0, 0.0))}
# file -> (cells per label, median depth, run factor on panel genes 0-5, development switch on genes 12-17 in SOM, seed)
STAGES = {"E6.75.h5ad": (130, 30.0, 0.7, 1.0, 21), "E7.25.h5ad": (140, 35.0, 0.85, 1.5, 22),
          "E7.5.h5ad": (180, 40.0, 1.0, 2.0, 23), "E8.0.h5ad": (220, 160.0, 2.2, 5.0, 24)}
# sha256 of the .h5ad build_test.py embryo wrote on `emb_env` at its default seed, computed with the builder, harness
# and t2_build of main at 4f39a11 before --expression trend was added (OPENBLAS/OMP/MKL threads 4)
PRE_TREND_SHA256 = {
    "fgw": "e93a4a08b2c60d58f731c5ec5fdd48b0b32e562b8a67fcbc44d396c370f76e3d",
    "fgw_shape0": "464ba84558b44e6f81882b02c775bb137a0ec64263cde3a27ea86335c7c57964",
    "fgw_shape1": "842f36d8a76710678a881ebf1aef47fd24b7b21364823d2697c8b752e56d4048",
    "fgw_aligned_series": "857b66becfd622b9a5ba2d1268dd1938de6c7f0347b16783cc83e40e54b5286b",
    "fgw_series": "0405193e09c7d0ba5626d560f3f09573b415bd14b910be1948315dbbe1d4db4b",
    "fgw_render_hi": "c22afcd9e6dbb19437d17dc6f7e457c7350ee2d36e57231fbba972fd6025354a",
    "fgw_depth": "797fe28941ff9c60d186b54131b17caa48d6f395846b828b3ce9efec388e4454",
    "fgw_blend": "1d18ff997e35e0eee3fd5a23e86fca8f3afeb8b3491e4dced4708835f9bd62e5",
    "fgw_herd": "59996b4479427aeddd1359d2351a52d9ae0fc42f173978aa017062e860cce019",
    "fgw_midpoint": "487169d8a49008b4f92bc8a9a0a69fec9a74c57a49a74b2e98bcdbc659fd0408",
    "pseudobulk_shift": "eac4f937a6964298126cf2ce2fc771e376420fcc7d01e8f3af6aaaed6a0b9451",
}
ARGV = {"fgw": ["--method", "fgw"], "fgw_shape0": ["--method", "fgw_shape0"], "fgw_shape1": ["--method", "fgw_shape1"],
        "fgw_aligned_series": ["--method", "fgw", "--move-coupling", "aligned", "--expression", "series"],
        "fgw_series": ["--method", "fgw", "--expression", "series"],
        "fgw_render_hi": ["--method", "fgw", "--expression", "render_hi"],
        "fgw_depth": ["--method", "fgw", "--expression", "depth"],
        "fgw_blend": ["--method", "fgw", "--expression", "blend"],
        "fgw_herd": ["--method", "fgw", "--herd"],
        "fgw_midpoint": ["--method", "fgw", "--draw", "midpoint"],
        "pseudobulk_shift": ["--method", "pseudobulk_shift"]}


def synthetic_stage(name: str, lacking=(), exe_prefix: str = "EXE-") -> ad.AnnData:
    """A released T2 embryo stage as read from its file (vec.prep.load_stage without a gene list): panel and outside
    genes, counts layer, X = log(1 + 1e4 c / N) with N over every gene, obs celltype and group, 3D coordinates.
    `lacking`: panel genes the file does not hold. `exe_prefix` renames the extraembryonic label."""
    per, depth, run, dev, seed = STAGES[name]
    rng = np.random.default_rng(seed)
    base = np.random.default_rng(5).gamma(2.0, 1.0, len(PANEL) + len(OUTSIDE))
    parts, obs, xyz = [], [], []
    for lab, (grp, centre) in LABELS.items():
        prof = base.copy()
        prof[:6] *= run
        if lab == "SOM":
            prof[12:18] *= dev
        prof /= prof.sum()
        d = np.maximum(np.round(np.exp(rng.normal(np.log(depth), 0.3, per))), 5)[:, None]
        C = rng.poisson(prof * d)
        C[C[:, :len(PANEL)].sum(1) == 0, 0] += 1
        parts.append(C)
        obs += [(lab.replace("EXE-", exe_prefix), grp)] * per
        xyz.append(np.asarray(centre) + rng.normal(size=(per, 3)) * 1.5)
    C = np.vstack(parts)
    N = C.sum(1, keepdims=True)
    a = ad.AnnData(X=sparse.csr_matrix(np.log1p(1e4 * C / N).astype(np.float32)),
                   obs=pd.DataFrame(obs, columns=["celltype", "group"],
                                    index=[f"{name[:-5]}_{i}" for i in range(len(C))]),
                   var=pd.DataFrame(index=PANEL + OUTSIDE))
    a.layers["counts"] = sparse.csr_matrix(C.astype(np.int32))
    a.obsm["spatial_3D"] = np.vstack(xyz).astype(np.float32)
    keep = [g for g in a.var_names if g not in set(lacking)]
    return a[:, keep].copy()


def embryo_test_board() -> Board:
    return Board(key="T2:embryo:test_interp", task="T2", genes=list(PANEL), needs_coords=True, split="test_interp",
                 setting="embryo", min_cells=100, max_cells=SECRET_MAX)


def harness_inputs(stages: dict) -> dict:
    """The four stages as the builder holds them after loading: the brackets on the panel (load_stage's gene list), the
    two earlier samples restricted to the panel genes they hold, in day order."""
    def held(a):
        return a[:, [g for g in PANEL if g in a.var_names]].copy()

    return {"e6_75": held(stages["E6.75.h5ad"]), "e7_25": held(stages["E7.25.h5ad"]),
            "before": stages["E7.5.h5ad"][:, PANEL].copy(), "after": stages["E8.0.h5ad"][:, PANEL].copy()}


def candidate_ctx(stages: dict, seed: int):
    """A harness context on the test geometry with the four stages (fgw_trend_render's own route)."""
    from vec import harness

    inputs = harness_inputs(stages)
    days = {"e6_75": 6.75, "e7_25": 7.25, "before": 7.5, "after": 8.0}
    view = harness.TaskView(key="expected", task="T2", kind="interp", setting="embryo",
                            test_board=embryo_test_board().key, roles=tuple(days), role_days=days,
                            role_stems={r: r for r in days}, reference_role="before", target_day=7.75)
    return harness.BuildContext(task=view, board=harness._candidate_board(embryo_test_board()), inputs=inputs,
                                seed=seed)


def install(bt, setattr_, reads: list, size_path: Path, lacking=(), exe_prefix: str = "EXE-") -> dict:
    """Point build_test's embryo builder at the synthetic stages and board: `setattr_(obj, name, value)` is
    monkeypatch.setattr in the tests (and plain setattr where the pre-change checksums were computed). The count
    decision is a stand-in, the size JSON a raw-mode result whose target RMS is the harness's size rule on the four
    synthetic stages, and pseudobulk_shift's size step a fixed scale (their own tests cover them). `lacking`: panel
    genes the E7.25 file lacks. Returns the stages."""
    from vec import harness

    stages = {f: synthetic_stage(f, lacking if f == "E7.25.h5ad" else (),
                                 exe_prefix if f == "E6.75.h5ad" else "EXE-") for f in STAGES}
    board = embryo_test_board()
    full = {f: synthetic_stage(f) for f in STAGES}
    rms = harness.target_rms(candidate_ctx(full, harness.BUILD_SEED_OFFSET + 1))
    size_path.write_text(json.dumps({"board": "embryo", "valid": True,
                                     "preregistered": {"z_mode": "raw", "pred_rms_full": rms}}))

    def load(path, genes=None, lineage=True, drop_duplicates=False):
        name = Path(path).name
        assert Path(path).parent.name == "T2" and drop_duplicates and lineage, path
        reads.append((name, genes is None))
        a = stages[name]
        if genes is None:
            return a.copy()
        missing = [g for g in genes if g not in a.var_names]
        if missing:  # as vec.prep.load_stage
            raise KeyError(f"{path}: {len(missing)} panel genes missing, e.g. {missing[:5]}")
        return a[:, list(genes)].copy()

    def size(coords, size_json, **kw):
        return np.asarray(coords, dtype=float) * 1.25, {"factor": 1.25}

    setattr_(bt, "load_stage", load)
    setattr_(bt, "load_boards", lambda: {board.key: board})
    setattr_(bt, "_t2_count", lambda b, n: {"decision": "synthetic"})
    setattr_(bt, "SIZES", {**bt.SIZES, "embryo": size_path})
    setattr_(bt.t2_build, "apply_size_procedure", size)
    return stages


def sha(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


# ---------------------------------------------------------------------------------------------------------------------
@pytest.fixture(scope="module")
def bt():
    sys.path.insert(0, str(BUILDS))
    import build_test

    return build_test


@pytest.fixture
def emb_env(bt, monkeypatch, own_repo, tmp_path):
    reads: list = []
    stages = install(bt, monkeypatch.setattr, reads, tmp_path / "size_embryo.json")
    return reads, stages


def _build(bt, tmp_path, name: str, *argv: str) -> tuple[Path, dict, str]:
    out = tmp_path / name / "T2embryo.h5ad"
    bt.main(["embryo", *argv, "--out", str(out)])
    text = out.with_suffix(".build.json").read_text()
    return out, json.loads(text), text


def _file_X(path: Path) -> np.ndarray:
    from vec.io import dense

    return dense(ad.read_h5ad(path).X).astype(np.float32)


def test_the_trend_option_is_accepted_and_off_by_default(bt):
    ap = bt.build_parser()
    base = ["embryo", "--method", "fgw", "--out", "runs/x/T2embryo_s1.h5ad"]
    assert ap.parse_args(base).expression == "mixture"
    assert ap.parse_args([*base, "--expression", "trend"]).expression == "trend"
    assert bt.EMBRYO_TREND_STAGES == {"e6_75": ("T2/E6.75.h5ad", 6.75), "e7_25": ("T2/E7.25.h5ad", 7.25)}
    help_text = " ".join(ap._subparsers._group_actions[0].choices["embryo"].format_help().split())
    assert "not an entry until the captain rules" in help_text
    # build_all.sh still builds the embryo s1 with the released cells
    run = [ln.split("#")[0] for ln in (BUILDS / "build_all.sh").read_text().splitlines()]
    assert "EMBRYO_EXPRESSION=mixture  " in run and not any("trend" in ln for ln in run if "embryo" in ln)
    assert sum("EMBRYO_EXPRESSION=" in ln for ln in run) == 1


@pytest.mark.parametrize("name", sorted(PRE_TREND_SHA256))
def test_every_existing_embryo_build_is_todays(bt, emb_env, tmp_path, capsys, name):
    """Without --expression trend each embryo build is byte for byte the file main's builder wrote before the option
    (pinned), it reads the brackets on the panel (and E7.25 for pseudobulk_shift) and nothing else."""
    reads, _ = emb_env
    out, rec, text = _build(bt, tmp_path, name, *ARGV[name])
    assert sha(out) == PRE_TREND_SHA256[name]
    want = [("E7.5.h5ad", False), ("E8.0.h5ad", False)]
    assert reads == (want + [("E7.25.h5ad", False)] if name == "pseudobulk_shift" else want)
    assert rec["n_cells"] == "max_cells" and str(SECRET_MAX) not in text + capsys.readouterr().out
    assert not any(k.startswith("trend") for k in rec["record"])


def test_with_every_gene_common_the_builder_writes_the_candidates_expression(bt, emb_env, tmp_path, capsys):
    """The file --expression trend writes holds fgw_trend_render's expression (`vec.harness` candidate on the four
    stages, the build's seed, the same size target) and fgw's positions; the build reads the two brackets on the panel
    and E6.75 and E7.25 without a gene list; the record holds the rule, the target day and 0 unrendered genes, and no
    cell number."""
    from vec import harness, t2_build

    reads, stages = emb_env
    out, rec, text = _build(bt, tmp_path, "trend", "--method", "fgw", "--expression", "trend")
    plain, _, _ = _build(bt, tmp_path, "plain", "--method", "fgw")
    assert reads[:4] == [("E7.5.h5ad", False), ("E8.0.h5ad", False), ("E6.75.h5ad", True), ("E7.25.h5ad", True)]
    want = harness.CANDIDATES["fgw_trend_render"].build(candidate_ctx(stages, harness.BUILD_SEED_OFFSET + 1))
    X = _file_X(out)
    np.testing.assert_array_equal(X, want.X)
    got = ad.read_h5ad(out)
    np.testing.assert_array_equal(np.asarray(got.obsm["spatial_3D"]),
                                  np.asarray(ad.read_h5ad(plain).obsm["spatial_3D"]))
    assert not np.array_equal(X, _file_X(plain))
    r = rec["record"]
    assert rec["entry"] == "fgw" and r["expression"] == "trend" and r["trend_panel_genes_unrendered"] == 0
    assert r["trend_rule"] == t2_build.TREND_COMMON_RULE and r["trend_target_day"] == 7.75
    assert rec["n_cells"] == "max_cells" and str(SECRET_MAX) not in text + capsys.readouterr().out
    assert not ({"n", "n_before", "n_after", "n_jittered"} & set(r)) and not any("moved" in k for k in r)


def _fgw_interp(bt, stages, board, seed):
    """fgw's own interpolation as the builder makes it (t = 1/2, the size JSON's target RMS)."""
    from vec import harness, t2_build
    from vec.io import dense

    target = json.loads(Path(bt.SIZES["embryo"]).read_text())["preregistered"]["pred_rms_full"]
    before, after = stages["E7.5.h5ad"][:, PANEL].copy(), stages["E8.0.h5ad"][:, PANEL].copy()
    Eb, Ea = (dense(a.X).astype(np.float32) for a in (before, after))
    res = t2_build.interpolate_clouds(np.asarray(before.obsm["spatial_3D"], float)[:, :3],
                                      np.asarray(after.obsm["spatial_3D"], float)[:, :3], Eb, Ea, 0.5,
                                      harness.count_for(board, before.n_obs), seed, target_rms=target)
    return res, before, after, Eb, Ea


def test_a_sample_lacking_two_panel_genes(bt, monkeypatch, own_repo, tmp_path, capsys):
    """E7.25 lacks two panel genes (load_stage with the panel would raise): the build goes through, the two genes'
    columns are fgw's mixture expression exactly, the other columns are trend_rendered_expression on the panel
    without the two genes (same draws, seed stream [seed, 7]), and the record says 2 and names no gene."""
    from vec import harness, t2_build

    reads: list = []
    stages = install(bt, monkeypatch.setattr, reads, tmp_path / "size.json", lacking=LACKING)
    with pytest.raises(KeyError, match="panel genes missing"):  # what the stage's file does to load_stage
        bt.load_stage(Path("T2/E7.25.h5ad"), genes=PANEL, drop_duplicates=True)
    reads.clear()
    out, rec, text = _build(bt, tmp_path, "lack", "--method", "fgw", "--expression", "trend")
    seed = harness.BUILD_SEED_OFFSET + 1
    res, before, after, Eb, Ea = _fgw_interp(bt, stages, embryo_test_board(), seed)
    X = _file_X(out)
    lack = [PANEL.index(g) for g in LACKING]
    common = [g for g in PANEL if g not in LACKING]
    mix = np.empty_like(X)
    mix[res.stage == 0], mix[res.stage == 1] = Eb[res.rows[res.stage == 0]], Ea[res.rows[res.stage == 1]]
    np.testing.assert_array_equal(X[:, lack], mix[:, lack])
    e675 = stages["E6.75.h5ad"][:, PANEL].copy()
    e725 = stages["E7.25.h5ad"][:, common].copy()
    sub, sub_rec = t2_build.trend_rendered_expression(before, after, res, common,
                                                      [(e675, 6.75), (e725, 7.25), (before, 7.5), (after, 8.0)],
                                                      7.75, [seed, 7])
    np.testing.assert_array_equal(X[:, [PANEL.index(g) for g in common]], sub)
    assert sub_rec["moved"]["before"] > 0 and sub_rec["moved"]["after"] > 0  # the common genes are rendered
    assert not np.array_equal(X[:, [PANEL.index(g) for g in common]], mix[:, [PANEL.index(g) for g in common]])
    r = rec["record"]
    assert r["trend_panel_genes_unrendered"] == 2 and r["trend_rule"] == t2_build.TREND_COMMON_RULE
    printed = text + capsys.readouterr().out
    assert all(f'"{g}"' not in printed for g in LACKING) and str(SECRET_MAX) not in printed


@pytest.mark.parametrize("extra", [["--method", "pseudobulk_shift"], ["--method", "fgw", "--herd"],
                                   ["--method", "fgw_shape0", "--herd"]])
def test_trend_is_refused_with_pseudobulk_shift_and_herd_before_reading_a_stage(bt, emb_env, tmp_path, extra):
    reads, _ = emb_env
    with pytest.raises(SystemExit) as e:
        bt.main(["embryo", *extra, "--expression", "trend", "--out", str(tmp_path / "e.h5ad")])
    assert "--expression trend is refused" in str(e.value)
    assert reads == [] and not (tmp_path / "e.h5ad").exists()


def test_a_sample_without_the_shared_extraembryonic_label_still_raises(bt, monkeypatch, own_repo, tmp_path):
    """E6.75's extraembryonic label renamed: no label is held by every sample, and the build stops with the existing
    error (no label, sample or factor is substituted) and writes no file."""
    reads: list = []
    install(bt, monkeypatch.setattr, reads, tmp_path / "size.json", exe_prefix="EMB-")
    with pytest.raises(ValueError, match="no extraembryonic label held by every sample"):
        bt.main(["embryo", "--method", "fgw", "--expression", "trend", "--out", str(tmp_path / "x" / "e.h5ad")])
    assert not (tmp_path / "x" / "e.h5ad").exists()
