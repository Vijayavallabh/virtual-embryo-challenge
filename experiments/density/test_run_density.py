"""Fast checks of the density experiment's bookkeeping on synthetic row indices (no competition data).

Run: uv run pytest -q experiments/density/test_run_density.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parent))
import run_density as rd  # noqa: E402

N_T, N_R, TC = 31671, 13295, 1583  # embryo_E8.0 sizes: target, reference, truth half


def pools(seed: int = 0):
    rng = np.random.default_rng(seed)
    wc = np.sort(rng.choice(N_T, int(round(0.1 * N_T)), replace=False))
    rwc = np.sort(rng.choice(N_R, int(round(0.1 * N_R)), replace=False))
    return np.setdiff1d(np.arange(N_T), wc), np.setdiff1d(np.arange(N_R), rwc)


def conds(seed: int, max_cells: int = 5000):
    pool, ref_pool = pools()
    return pool, ref_pool, {lab: (meta, src, idx, pert)
                            for lab, meta, src, idx, pert in rd.conditions(pool, ref_pool, TC, max_cells, seed)}


def test_rule_max_cells():
    assert rd.rule_max_cells(1583, 1330) == 5000  # floor binds
    assert rd.rule_max_cells(2687, 2483) == 3 * 2687
    assert rd.rule_max_cells(1241, 5872) == 3 * 5872  # reference working copy larger than the truth half


def test_labels_unique_and_complete():
    pool, ref_pool = pools()
    labels = [lab for lab, *_ in rd.conditions(pool, ref_pool, TC, 5000, 0)]
    assert len(labels) == len(set(labels)) == 66
    for pat in ("random", "shuffle50", "refstage"):
        assert {f"{pat}_{r}x" for r in rd.RATIOS} | {f"{pat}_{f:g}max" for f in rd.MAX_FRACS} <= set(labels)


def test_random_draws_keep_the_first_run_recipe():
    pool, ref_pool, c = conds(seed=3)
    perm = np.random.default_rng([3, 1]).permutation(len(pool))
    ref_perm = np.random.default_rng([3, 4]).permutation(len(ref_pool))
    for r in rd.RATIOS:
        n = min(int(round(r * TC)), len(pool))
        assert np.array_equal(c[f"random_{r}x"][2], pool[np.sort(perm[:n])])
        assert np.array_equal(c[f"head_{r}x"][2], pool[:n])
        n_ref = min(int(round(r * TC)), len(ref_pool))
        assert np.array_equal(c[f"refstage_{r}x"][2], ref_pool[np.sort(ref_perm[:n_ref])])
    # max_cells levels are prefixes of the same permutation (nested with the truth multiples)
    assert np.array_equal(c["random_1max"][2], pool[np.sort(perm[:5000])])
    assert c["random_0.5max"][0]["basis"] == "max" and c["random_0.5max"][0]["ratio"] == pytest.approx(2500 / TC)


def test_blocks_fit_and_centres_vary():
    starts = {}
    for seed in range(5):
        pool, ref_pool, c = conds(seed)
        for lab, (meta, src, idx, _) in c.items():
            if meta["pattern"] not in ("block", "refblock"):
                continue
            p = pool if src == "target" else ref_pool
            s = meta["block_start"]
            assert meta["block_draw"] == rd.BLOCK_DRAW
            assert 0 <= s <= len(p) - len(idx) and np.array_equal(idx, p[s:s + len(idx)])
            starts.setdefault(lab, set()).add(s)
    assert len(starts["block_1x"]) == 5 and len(starts["refblock_1x"]) == 5
    assert starts["refblock_8x"] == {0}  # 8x exceeds this reference pool: the block is the whole pool


def test_shuffle50_permutes_half_of_the_rows_of_the_random_cells():
    pool, ref_pool, c = conds(seed=1)
    meta, src, idx, pert = c["shuffle50_2x"]
    assert np.array_equal(idx, c["random_2x"][2]) and pert == ("shuffle", 0.5, 3)
    assert c["shuffle50_1max"][3][2] == 100 + len(rd.RATIOS) + 3
    n = len(idx)
    x = np.arange(n, dtype=np.float32)[:, None]
    pred = rd.Stage(X=x, coords=np.zeros((n, 3), np.float32))
    out, rms = rd.perturb(pred, pert, 1, 1.0, 1.0)
    moved = int((out.X[:, 0] != x[:, 0]).sum())
    assert rms == 0.0 and sorted(out.X[:, 0]) == sorted(x[:, 0]) and 0.45 * n < moved <= 0.5 * n + 1


def test_resolve_and_skip(tmp_path):
    assert rd.resolve(str(tmp_path / "E10.5.h5ad"))[0] is None
    (tmp_path / "E10.5_val.h5ad").write_bytes(b"")
    p, note = rd.resolve(str(tmp_path / "E10.5.h5ad"))
    assert p == tmp_path / "E10.5_val.h5ad" and "using" in note
    (tmp_path / "E10.5_other.h5ad").write_bytes(b"")
    assert rd.resolve(str(tmp_path / "E10.5.h5ad"))[0] is None
    assert rd.resolve(str(tmp_path / "E7.5.h5ad"))[0] is None  # E7.5* must not match E7.75 or E10.5 files


def test_old_block_rows_are_refused(tmp_path):
    f = tmp_path / "t_s0.jsonl"
    f.write_text(json.dumps(dict(condition="random_1x", pattern="random")) + "\n")
    assert rd.done_conditions(f) == {"random_1x"}
    f.write_text(f.read_text() + json.dumps(dict(condition="block_1x", pattern="block", block_start=5)) + "\n")
    with pytest.raises(SystemExit):
        rd.done_conditions(f)
