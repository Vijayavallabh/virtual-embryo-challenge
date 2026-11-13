"""`experiments/builds/build_test.py heart --trend` (b0_trend in both heart files: the captain's ruling of 2026-11-13,
docs/PLAN.md decision log (w); ledger row H11), on synthetic stages and a synthetic test board whose limit is never
shown: the refusals, a repeatable build, today's build without the flag unchanged, a record with no cell number, and
build_all.sh's two heart lines. Nothing is read from data/ and nothing is scored."""

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
GENES = [f"g{i}" for i in range(120)]
LINEAGES = np.array(["CM", "Endothelial", "Untracked", "Tiny"])
SECRET_MAX = 733  # stands for the test board's max_cells: never printed, never in a record
# file -> (cells, seed, shift of the CM state module, cells of the lineage too small to be shifted)
STAGES = {"E10.5.h5ad": (2600, 2, 0.0, 60), "E9.5.h5ad": (1500, 3, 0.6, 40), "E8.75.h5ad": (1300, 5, 1.1, 30),
          "E8.5.h5ad": (1100, 7, 1.5, 10), "E8.25_late.h5ad": (900, 11, 1.8, 10)}
# sha256 of the .h5ad build_test.py heart wrote on `heart_env` at its default seed, without --trend, computed with the
# builder and harness of main at 7b04d82 (before --trend was added)
PRE_TREND_SHA256 = {"s1": "cf887247be3d77a3fbbb6a99848a080137481f45a59bb97e4be3b3187e02e31e",
                    "s2_trim": "fc37df8923c87f39e6432eec45a2d8d49eae7c4e3b572fa057588bce949d5b9a"}


def synthetic_stage(n: int, seed: int, shift: float, tiny: int) -> ad.AnnData:
    """A T2 heart stage as vec.prep.load_stage returns it (panel genes, counts layer, obs["merged"], obs["core"],
    3D coordinates): CM and Endothelial lineages with a within-lineage state module, Untracked cells and `tiny` cells of
    a lineage below min_cells. `shift` moves the CM module between stages."""
    rng = np.random.default_rng(seed)
    types = rng.choice(3, size=n, p=[0.5, 0.35, 0.15])
    types[:tiny] = 3
    pi = np.random.default_rng(1).gamma(1.0, 1.0, len(GENES))
    pi /= pi.sum()
    a = np.zeros(len(GENES))
    a[:30] = np.random.default_rng(4).choice([-1.0, 1.0], 30) * 1.5
    s = rng.normal(size=n)
    lam = pi * np.exp(np.outer(s, a) + np.outer(shift * s * (types == 0), a) + 0.3 * (types == 1)[:, None]
                      + 0.2 * shift * (types == 1)[:, None] * np.linspace(-1, 1, len(GENES)))
    depth = np.maximum(np.round(np.exp(rng.normal(np.log(80.0), 0.4, n))), 5)[:, None]
    C = rng.poisson(lam / lam.sum(1, keepdims=True) * depth).astype(np.int32)
    C[C.sum(1) == 0, 0] = 1
    X = np.log1p(1e4 * C / C.sum(1, keepdims=True)).astype(np.float32)
    centres = np.array([[0.0, 0.0, 0.0], [5.0, 0.0, 0.0], [0.0, 5.0, 0.0], [5.0, 5.0, 0.0]])
    out = ad.AnnData(X=sparse.csr_matrix(X), obs=pd.DataFrame(
        {"celltype": [f"type{t}" for t in types], "merged": LINEAGES[types], "core": types == 0},
        index=[f"c{seed}_{i}" for i in range(n)]), var=pd.DataFrame(index=GENES))
    out.layers["counts"] = sparse.csr_matrix(C)
    out.obsm["spatial_3D"] = (centres[types] + rng.normal(size=(n, 3)) * 1.5).astype(np.float32)
    return out


def heart_test_board() -> Board:
    return Board(key="T2:heart:test_extrap", task="T2", genes=list(GENES), needs_coords=True, split="test_extrap",
                 setting="heart", min_cells=100, max_cells=SECRET_MAX)


def install(bt, setattr_, reads: list) -> None:
    """Point build_test's heart builder at the synthetic stages and board: `setattr_(obj, name, value)` is
    monkeypatch.setattr in the tests (and plain setattr where the pre-change checksum was computed). The count
    decision and the size step are stand-ins (their own tests cover them); everything else is the builder's."""
    stages = {f: synthetic_stage(*v) for f, v in STAGES.items()}
    board = heart_test_board()

    def load(path, genes=None, drop_duplicates=False):
        name = Path(path).name
        assert Path(path).parent.name == "T2" and genes == list(GENES) and drop_duplicates, path
        reads.append(name)
        return stages[name].copy()

    def size(coords, size_json, hypothesis=None, **kw):
        return np.asarray(coords, dtype=float) * 1.25, {"hypothesis": hypothesis, "factor": 1.25}

    setattr_(bt, "load_stage", load)
    setattr_(bt, "load_boards", lambda: {board.key: board})
    setattr_(bt, "_t2_count", lambda b, n: {"decision": "synthetic"})
    setattr_(bt.t2_build, "size_or_hold", size)


def sha(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


# ---------------------------------------------------------------------------------------------------------------------
@pytest.fixture(scope="module")
def bt():
    sys.path.insert(0, str(BUILDS))
    import build_test

    return build_test


@pytest.fixture
def heart_env(bt, monkeypatch, own_repo):
    reads: list[str] = []
    install(bt, monkeypatch.setattr, reads)
    return reads


def _build(bt, tmp_path, name: str, *argv: str) -> tuple[Path, dict, str]:
    out = tmp_path / name / "T2heart.h5ad"
    bt.main(["heart", *argv, "--out", str(out)])
    text = out.with_suffix(".build.json").read_text()
    return out, json.loads(text), text


@pytest.mark.parametrize("extra", [["--entry", "s2"], ["--entry", "s1", "--depth"], ["--entry", "s1", "--herd"],
                                   ["--entry", "s1", "--cfm"], ["--entry", "s2_trim", "--cfm"]])
def test_trend_is_refused_with_s2_depth_herd_and_cfm_before_reading_anything(bt, heart_env, tmp_path, extra):
    with pytest.raises(SystemExit) as e:
        bt.main(["heart", *extra, "--trend", "--out", str(tmp_path / "h.h5ad")])
    assert "--trend is refused" in str(e.value)
    assert heart_env == [] and not (tmp_path / "h.h5ad").exists()


def test_trend_flag_is_off_by_default(bt):
    ap = bt.build_parser()
    assert ap.parse_args(["heart", "--entry", "s1", "--out", "runs/x/h.h5ad"]).trend is False
    assert ap.parse_args(["heart", "--entry", "s2_trim", "--trend", "--out", "runs/x/h.h5ad"]).trend is True
    assert bt.HEART_TREND_STAGES == {"e8_75": ("T2/E8.75.h5ad", 8.75), "e8_5": ("T2/E8.5.h5ad", 8.5),
                                     "e8_25_late": ("T2/E8.25_late.h5ad", 8.3)}


@pytest.mark.parametrize("entry", ["s1", "s2_trim"])
def test_a_build_without_trend_is_todays(bt, heart_env, tmp_path, capsys, entry):
    """The file without --trend is byte for byte the one main's builder wrote before --trend (pinned), it reads E10.5
    and E9.5 only, and its record has no "trend" key."""
    out, rec, text = _build(bt, tmp_path, "plain", "--entry", entry)
    assert sha(out) == PRE_TREND_SHA256[entry]
    assert heart_env == ["E10.5.h5ad", "E9.5.h5ad"]
    assert "trend" not in rec and rec["n_cells"] == "max_cells"
    assert str(SECRET_MAX) not in text + capsys.readouterr().out


@pytest.mark.parametrize("entry", ["s1", "s2_trim"])
def test_a_trend_build_repeats_byte_for_byte_and_differs_from_b0(bt, heart_env, tmp_path, entry):
    a, rec_a, text_a = _build(bt, tmp_path, "a", "--entry", entry, "--trend")
    b, rec_b, text_b = _build(bt, tmp_path, "b", "--entry", entry, "--trend")
    plain, _, _ = _build(bt, tmp_path, "plain", "--entry", entry)
    assert a.read_bytes() == b.read_bytes() and text_a == text_b
    assert sha(a) != sha(plain)
    assert heart_env[:5] == ["E10.5.h5ad", "E9.5.h5ad", "E8.75.h5ad", "E8.5.h5ad", "E8.25_late.h5ad"]
    assert rec_a["entry"] == f"B0 {entry} ({'H1' if entry == 's1' else 'H2'})"  # the size hypothesis is the entry's
    assert (rec_a["trim"] is None) == (entry == "s1")


def test_the_trend_record_holds_no_cell_number(bt, heart_env, tmp_path, capsys):
    """The "trend" block: the rule, the five stage roles and days (the reference last, start before it), the transition
    weights, per lineage the named fields only, the medians; no key that names a count and no integer anywhere (so no
    stage's, pool's or file's cell number), and the board's limit nowhere in the record or the printout."""
    _, rec, text = _build(bt, tmp_path, "t", "--entry", "s2_trim", "--trend")
    tr = rec["trend"]
    assert set(tr) == {"rule", "stages", "days", "weights", "budget_days", "lineages", "median_cosine",
                       "median_norm_ratio"}
    assert tr["rule"] == "vec.harness.B0_TREND_RULE"
    assert tr["stages"] == ["e8_25_late", "e8_5", "e8_75", "start", "ref"]
    assert tr["days"] == [8.3, 8.5, 8.75, 9.5, 10.5] and tr["budget_days"] == pytest.approx(0.5)
    np.testing.assert_allclose(np.round(tr["weights"], 4), [0.0500, 0.1095, 0.4118, 0.4287])
    assert {r["lineage"] for r in tr["lineages"]} == {"CM", "Endothelial"}
    for r in tr["lineages"]:
        assert set(r) == set(bt.TREND_LINEAGE_FIELDS)
        assert r["stages"] == tr["stages"] and not r["own_factor"]
    assert np.isfinite([tr["median_cosine"], tr["median_norm_ratio"]]).all()

    def walk(x, path="trend"):
        if isinstance(x, dict):
            for k, v in x.items():
                assert not k.startswith("n_") and not any(w in k.lower() for w in ("cell", "count", "size", "min_")), k
                walk(v, f"{path}.{k}")
        elif isinstance(x, list):
            for i, v in enumerate(x):
                walk(v, f"{path}[{i}]")
        else:  # strings, floats and the own_factor flags only: an integer could be a cell number
            assert isinstance(x, (str, float, bool)), (path, type(x))

    walk(tr)
    printed = capsys.readouterr().out
    assert str(SECRET_MAX) not in text + printed


def test_build_all_puts_trend_in_both_heart_lines():
    """The captain's ruling of 2026-11-13 (docs/PLAN.md decision log (w)): both heart files are b0_trend."""
    text = (BUILDS / "build_all.sh").read_text()
    assert "build_test.py heart --entry s1 --trend --out $O/builds/T2heart_s1.h5ad" in text
    assert "build_test.py heart --entry s2_trim --trend --out $O/builds/T2heart_s2.h5ad" in text
    assert text.count("build_test.py heart ") == 2
