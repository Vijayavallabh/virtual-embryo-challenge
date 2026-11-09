"""experiments/harness/diag_e12_halves.py (row E12d of experiments/harness/TRIES.md), on synthetic stages only: the
halves are alternating sections, each world's truth is half A alone, the predictions are the harness's own builders on
half B, the reading follows its row, and nothing written or printed is a cell count."""

from __future__ import annotations

import importlib.util
import json
import re
import sys
import types
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from vec import adoption, harness
from vec.boards import Board
from vec.evaluate import Stage
from vec.harness import BUILD_SEED_OFFSET, CANDIDATES
from vec.io import COORD_KEY, dense

REPO = Path(__file__).resolve().parents[1]
GENES = [f"g{i}" for i in range(18)]
LEVELS = 12  # z levels of a synthetic stage
SIZES = {"E7.5": 413, "E8.0": 527}  # the synthetic stages' cells; never printed or written
SECRET = {"min_cells": 777, "max_cells": 4321, "truth_cells": 2468, "ref_cells": 1357}  # stand for the test board's


@pytest.fixture(scope="module")
def dg():
    spec = importlib.util.spec_from_file_location("harness_diag_e12_halves",
                                                  REPO / "experiments" / "harness" / "diag_e12_halves.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def _stage(name: str, seed: int) -> ad.AnnData:
    """A stage as load_stage returns it: z on whole section levels, x the thinnest axis, expression left-right symmetric
    (a function of |x|, y and the section's z), so the pair has FGW's mirror ambiguity. E8.0 carries a strong change on
    six genes, so a mixture of the two brackets has DE genes against E7.5 and the scorer's anchors are complete."""
    n = SIZES[name]
    rng = np.random.default_rng(seed)
    z_level = rng.integers(0, LEVELS, n)
    C = np.c_[rng.normal(size=n) * 0.5, rng.normal(size=n) * 1.6, z_level * 15.0]  # z exactly on the section levels
    centres = np.random.default_rng(0).normal(size=(len(GENES), 3)) * [0.4, 1.2, 40.0] + [0.0, 0.0, 80.0]
    F = np.c_[np.abs(C[:, 0]), C[:, 1], C[:, 2]]
    E = np.exp(-((F[:, None, :] - centres[None]) ** 2 / np.array([0.3, 1.0, 1200.0])).sum(-1))
    rate = 40.0 * E + 0.3
    if name == "E8.0":
        rate[:, :6] += 25.0
    counts = rng.poisson(rate)
    lib = counts.sum(1, keepdims=True).clip(min=1)
    a = ad.AnnData(X=np.log1p(1e4 * counts / lib).astype(np.float32),
                   obs=pd.DataFrame({"celltype": rng.choice(["a", "b"], n)}, index=[f"c{i}" for i in range(n)]),
                   var=pd.DataFrame(index=GENES))
    a.obsm[COORD_KEY] = C
    return a


def _patch(dg, monkeypatch, tmp_path, seed=3) -> Path:
    """Synthetic boards and stages in place of boards/ and data/, and a small scorer window; returns the data folder."""
    boards = {"T2:embryo:val_interp": Board(key="T2:embryo:val_interp", task="T2", genes=[], needs_coords=True,
                                            setting="embryo", split="val_interp"),
              "T2:embryo:test_interp": Board(key="T2:embryo:test_interp", task="T2", genes=GENES, needs_coords=True,
                                             setting="embryo", split="test_interp", **SECRET)}
    stages = {n: _stage(n, seed + i) for i, n in enumerate(SIZES)}
    monkeypatch.setattr(dg, "load_boards", lambda: boards)

    def load(path, genes=None, drop_duplicates=False):
        return stages[Path(path).stem][:, genes].copy()

    monkeypatch.setattr(dg, "load_stage", load)
    # a small scorer window and pseudo-board counts for a stage of a few hundred cells (the real ones read
    # SUBSAMPLE_FRAC and the released files' sizes)
    monkeypatch.setattr(dg, "FRAC", 1.0)

    def pseudo_counts(n_target, n_ref):
        return {"truth_cells": n_target // 2, "ref_cells": n_ref, "min_cells": 100, "max_cells": 300}

    monkeypatch.setattr(dg.harness, "pseudo_counts", pseudo_counts)
    data = tmp_path / "data"
    (data / "T2").mkdir(parents=True)
    for n in SIZES:
        (data / "T2" / f"{n}.h5ad").write_bytes(b"synthetic stage " + n.encode())
    return data


def _ctx(dg, a_stage, b_stage, board, seed):
    view = dg.task_view()
    return harness.BuildContext(task=view, board=board, inputs={"before": a_stage, "after": b_stage}, seed=seed)


def test_section_halves_alternate_snapped_z_levels(dg):
    z = np.array([0.0, 15.0, 30.0, 45.0, 60.0, 75.0, 0.0, 15.0, 30.0])
    coords = np.c_[np.zeros(9), np.zeros(9), z]
    a, b = dg.section_halves(coords)
    assert list(np.flatnonzero(a)) == [0, 2, 4, 6, 8]  # levels 0, 2, 4 (even) in ascending z
    assert list(np.flatnonzero(b)) == [1, 3, 5, 7]
    # the level is the snapped one: float noise on a level does not move a cell to the other half
    noisy = coords.copy()
    noisy[3, 2] += 1e-9
    assert list(np.flatnonzero(dg.section_halves(noisy)[0])) == [0, 2, 4, 6, 8]


def test_each_half_is_a_section_set_and_the_two_partition_a_stage(dg):
    a = _stage("E7.5", 11)
    half_a, half_b = dg.split(a)
    assert half_a.n_obs + half_b.n_obs == a.n_obs
    assert set(half_a.obs_names).isdisjoint(half_b.obs_names)
    assert set(half_a.obs_names) | set(half_b.obs_names) == set(a.obs_names)
    za, zb = (np.unique(np.round(h.obsm[COORD_KEY][:, 2]).astype(int)) for h in (half_a, half_b))
    assert set(za).isdisjoint(zb)


def test_predict_is_exactly_the_harness_builders(dg):
    b_b = _stage("E7.5", 21)
    b_a = _stage("E8.0", 22)
    board = Board(key="pseudo:E12d", task="T2", genes=GENES, needs_coords=True, setting="embryo", min_cells=100,
                  max_cells=350)
    cboard = harness._candidate_board(board)
    ctx = harness.BuildContext(task=dg.task_view(), board=cboard, inputs={"before": b_b, "after": b_a},
                               seed=BUILD_SEED_OFFSET + 3)
    for name, coupling in dg.CANDIDATES.items():
        stage, res = dg.predict(ctx, coupling)
        built = CANDIDATES[name].build(ctx)
        np.testing.assert_array_equal(stage.coords, built.coords)
        np.testing.assert_array_equal(stage.X, built.X)
        assert res.coupling is not None and res.record["t_pos"] == 0.0
    assert dg.CANDIDATES == {"fgw_shape0": "fgw", "fgw_shape0_aligned": "aligned"}


def test_world_truth_is_half_a_alone_in_the_earlier_brackets_geometry(dg):
    a_b, a_a = _stage("E7.5", 31), _stage("E8.0", 32)
    n = 2 * min(a_b.n_obs, a_a.n_obs)
    truth, res = dg.world_truth(a_b, a_a, GENES, n, seed=2_000_004, move_coupling="fgw")
    Eb, Ea = dense(a_b.X), dense(a_a.X)
    # half of the truth's cells are the earlier stage's real cells, half the later stage's, none repeated
    assert truth.n == n and (res.stage == 0).sum() == n // 2 == (res.stage == 1).sum()
    assert len(np.unique(res.rows[res.stage == 0])) == n // 2 and len(np.unique(res.rows[res.stage == 1])) == n // 2
    np.testing.assert_array_equal(truth.X[res.stage == 0], Eb[res.rows[res.stage == 0]].astype(np.float32))
    np.testing.assert_array_equal(truth.X[res.stage == 1], Ea[res.rows[res.stage == 1]].astype(np.float32))
    # E7.5's geometry: t_pos = 0 leaves the earlier stage's cells where they are, up to the size rule's similarity
    Cb = np.asarray(a_b.obsm[COORD_KEY], float)[:, :3]
    P = truth.coords[res.stage == 0].astype(float)
    Q = Cb[res.rows[res.stage == 0]]
    P, Q = P - P.mean(0), Q - Q.mean(0)
    s = np.linalg.norm(P) / np.linalg.norm(Q)
    np.testing.assert_allclose(P, Q * s, atol=1e-3)
    assert res.record["t"] == 0.5 and res.record["t_pos"] == 0.0
    # the aligned world differs only in where the mapped later half lands
    aligned, res_a = dg.world_truth(a_b, a_a, GENES, n, seed=2_000_004, move_coupling="aligned")
    np.testing.assert_array_equal(aligned.X, truth.X)
    assert not np.allclose(aligned.coords[res.stage == 1], truth.coords[res.stage == 1])
    assert res_a.record["move_coupling"] == "aligned"


def _summary(gain_a: float, ahead_a: int, gain_f: float, usable: int = 10) -> dict:
    return {"aligned_placed": {"gain_mean": gain_a, "seeds_ahead": ahead_a, "n_seeds": usable},
            "fgw_placed": {"gain_mean": gain_f, "n_seeds": usable}}


@pytest.mark.parametrize("gain_a, ahead_a, gain_f, full, reading", [
    (0.5, 9, 0.0, True, "supports"),
    (0.5, 9, -0.25, True, "supports"),
    (0.5, 9, -0.4, True, "no clear read"),  # below -0.3 in the FGW-placed world: not supports
    (0.5, 7, 0.0, True, "no clear read"),  # ahead in too few seeds of the aligned-placed world
    (0.25, 9, 0.0, True, "no clear read"),  # mean not above 0.3 in the aligned-placed world
    (-0.4, 2, 0.5, True, "against"),  # loses in the world that favours it
    (-0.2, 3, -0.2, True, "no clear read"),
    (0.5, 9, 0.0, False, "partial run: no reading"),
])
def test_reading_follows_the_row(dg, gain_a, ahead_a, gain_f, full, reading):
    assert dg.reading(_summary(gain_a, ahead_a, gain_f), full) == reading


def test_reading_bars_are_the_adoption_rules(dg):
    assert adoption.MIN_GAIN == 0.3 and adoption.MIN_AHEAD == 8 and adoption.MAX_LOSS == 0.3
    assert dg.MIN_USABLE_SEEDS == adoption.MIN_AHEAD
    # exactly at the bars: not above +0.3 is no support; not below -0.3 in the FGW-placed world is no veto
    assert dg.reading(_summary(0.3, 9, 0.0), True) == "no clear read"
    assert dg.reading(_summary(0.5, 8, -0.3), True) == "supports"
    assert dg.reading(_summary(-0.3, 0, 0.0), True) == "no clear read"


def test_a_world_with_too_few_usable_seeds_gives_no_reading(dg):
    assert dg.reading(_summary(0.5, 9, 0.0, usable=9), True) == "supports"
    assert dg.reading(_summary(0.5, 8, 0.0, usable=8), True) == "supports"
    assert dg.reading(_summary(0.5, 7, 0.0, usable=7), True) == "no clear read"
    assert dg.reading(_summary(-0.5, 2, 0.0, usable=7), True) == "no clear read"  # not "against" on so few seeds


def _world_row(entry_score: float, cand_score: float, degenerate: bool = False) -> dict:
    raw = {m: 0.1 for m in ("d2_shape", "occupancy_dice", "scale_log_ratio", "neighborhood_mmd")}
    row = {"scores": {}, "skills": {},
           "raw": {"fgw_shape0": dict(raw), "fgw_shape0_aligned": {**raw, "d2_shape": 0.08}}}
    if degenerate:
        row["degenerate"] = {"inverted": ["d2_shape"], "missing": []}
    else:
        skills = {m: 0.5 for m in ("d2_shape", "occupancy_dice", "scale_log_ratio", "neighborhood_mmd", "de_score",
                                   "de_direction", "mmd_u", "variogram")}
        row["scores"] = {"fgw_shape0": entry_score, "fgw_shape0_aligned": cand_score}
        row["skills"] = {"fgw_shape0": skills, "fgw_shape0_aligned": {**skills, "d2_shape": 0.6}}
    return row


def test_world_summary_scores_usable_seeds_only_and_keeps_raw_changes_for_all(dg):
    rows = [_world_row(60.0, 60.5), _world_row(61.0, 60.9), _world_row(0.0, 0.0, degenerate=True)]
    s = dg.world_summary(rows)
    assert s["n_seeds"] == 2 and s["n_seeds_total"] == 3 and len(s["degenerate"]) == 1
    assert s["gain_mean"] == pytest.approx(0.2) and s["seeds_ahead"] == 1
    assert s["score_mean"]["fgw_shape0_aligned"] == pytest.approx(60.7)
    assert s["skill_change_mean"]["d2_shape"] == pytest.approx(0.1) and s["skill_change_mean"]["mmd_u"] == 0.0
    # raw changes are over every seed, degenerate included
    assert s["raw_change_mean"]["d2_shape"] == pytest.approx(-0.02) and s["raw_change_mean"]["occupancy_dice"] == 0.0
    # no usable seed: no scores and no reading, nothing to divide by
    none = dg.world_summary([_world_row(0.0, 0.0, degenerate=True)])
    assert none["n_seeds"] == 0 and none["gain_mean"] is None and none["score_mean"] is None
    assert none["seeds_ahead"] == 0


def test_score_world_marks_a_degenerate_seed_and_refuses_a_broken_floor(dg, monkeypatch):
    metrics = {m: 0.1 for m in ("d2_shape", "occupancy_dice", "scale_log_ratio", "neighborhood_mmd")}
    complete = {m: {"floor": 1.0, "ceiling": 0.0} for m in dg.harness.RANKED_METRICS["T2"]}

    def fake_score_board(extras_inverted=(), anchors=complete, score=50.0):
        return types.SimpleNamespace(score=score, skills={}, metrics=metrics,
                                     extras={"inverted_anchors": list(extras_inverted)}, anchors=anchors)

    ref = Stage(X=np.zeros((4, 2)), coords=np.zeros((4, 3)))
    preds = {"fgw_shape0": ref, "fgw_shape0_aligned": ref}
    monkeypatch.setattr(dg, "subsample", lambda stage, frac, seed: stage)
    monkeypatch.setattr(dg, "working_copies", lambda *a: (ref, ref, ref))
    monkeypatch.setattr(dg, "panel_metrics", lambda *a, **kw: metrics)
    # an inverted shape anchor: no scores, the raw metrics kept, no exception
    monkeypatch.setattr(dg, "score_board", lambda *a, **kw: fake_score_board(extras_inverted=("d2_shape",)))
    out = dg.score_world(None, ref, ref, preds, seed=0)
    assert out["degenerate"]["inverted"] == ["d2_shape"] and out["scores"] == {} and out["skills"] == {}
    assert set(out["raw"]) == set(preds) and out["raw"]["fgw_shape0"]["d2_shape"] == 0.1
    # a missing anchor is degenerate too
    monkeypatch.setattr(dg, "score_board", lambda *a, **kw: fake_score_board(anchors={}))
    assert dg.score_world(None, ref, ref, preds, seed=0)["degenerate"]["missing"]
    # complete, sound anchors with a floor that does not score 50: a bug, refused
    monkeypatch.setattr(dg, "score_board", lambda *a, **kw: fake_score_board(score=49.0))
    with pytest.raises(ValueError, match="the floor scores"):
        dg.score_world(None, ref, ref, preds, seed=0)


def _leaves(x, keys, ints, texts):
    if isinstance(x, dict):
        keys.update(x)
        for v in x.values():
            _leaves(v, keys, ints, texts)
    elif isinstance(x, list):
        for v in x:
            _leaves(v, keys, ints, texts)
    elif isinstance(x, bool):
        pass
    elif isinstance(x, int):
        ints.add(x)
    elif isinstance(x, str):
        texts.append(x)


COUNT_KEYS = {"cells", "n", "n_obs", "n_cells", "n_before", "n_after", "n_jittered", "n_anchor", "coupling_nonzero",
              "truth_cells", "ref_cells", "min_cells", "max_cells"}


def test_main_runs_end_to_end_writes_no_cell_count_and_refuses_a_rerun(dg, tmp_path, monkeypatch, own_repo, capsys):
    data = _patch(dg, monkeypatch, tmp_path)
    full_b = {name: _stage(name, 3 + i) for i, name in enumerate(SIZES)}  # the stages _patch serves (same seeds)
    half_b_cells = {name: set(dg.split(a[:, GENES].copy())[1].obs_names) for name, a in full_b.items()}
    out = tmp_path / "out"
    seen = {"truth": 0, "predict": 0}
    real_truth, real_predict = dg.world_truth, dg.predict

    def world_truth(*a, **kw):
        seen["truth"] += 1
        return real_truth(*a, **kw)

    def predict(ctx, coupling):
        assert ctx.seed == BUILD_SEED_OFFSET and ctx.board.truth_cells is None  # the candidate's view of the board
        before, after = ctx.inputs["before"], ctx.inputs["after"]
        # the predictions are built from half B alone: its even-numbered sections are never among them
        assert set(before.obs_names) == half_b_cells["E7.5"] and set(after.obs_names) == half_b_cells["E8.0"]
        seen["predict"] += 1
        return real_predict(ctx, coupling)

    monkeypatch.setattr(dg, "world_truth", world_truth)
    monkeypatch.setattr(dg, "predict", predict)
    dg.main(["--out", str(out), "--data", str(data), "--seeds", "0"])
    printed = capsys.readouterr().out
    assert seen == {"truth": 2, "predict": 2}  # two worlds and two candidates per seed
    rep = json.loads((out / "e12d_halves.json").read_text())
    txt = (out / "e12d_halves.txt").read_text()
    assert printed.endswith(txt) and txt.endswith(f"UTC {rep['stamp']['utc']}\n")
    assert rep["reading"] == "partial run: no reading" and rep["n_seeds"] == 1 and rep["seeds"] == [0]
    assert set(rep["summary"]) == {"fgw_placed", "aligned_placed"} and set(rep["worlds"]) == set(rep["summary"])
    s = rep["summary"]["fgw_placed"]
    assert set(s) == {"score_mean", "gain_mean", "gain_sd", "seeds_ahead", "n_seeds", "n_seeds_total", "degenerate",
                      "skill_change_mean", "raw_change_mean"}
    assert s["n_seeds"] == 1 and s["degenerate"] == []  # this synthetic pair's anchors are sound at seed 0
    assert s["gain_mean"] == pytest.approx(rep["per_seed"][0]["worlds"]["fgw_placed"]["gain"])
    assert set(s["raw_change_mean"]) == {"d2_shape", "occupancy_dice", "scale_log_ratio", "neighborhood_mmd"}
    # the floor scored 50 (checked inside the run) and both candidates carry the ranked skills
    row = rep["per_seed"][0]
    assert set(row["worlds"]["aligned_placed"]["scores"]) == {"fgw_shape0", "fgw_shape0_aligned"}
    skills = row["worlds"]["fgw_placed"]["skills"]["fgw_shape0"]
    assert set(skills) >= {"d2_shape", "occupancy_dice", "neighborhood_mmd"}
    assert row["build_seed"] == BUILD_SEED_OFFSET and row["truth_seed"] == dg.TRUTH_SEED_OFFSET
    # the mechanism read is there for both halves and both couplings, target-free
    assert set(rep["mechanism"]["B"]) == {"fgw", "aligned"}
    assert set(rep["mechanism"]["A"]) == {"fgw_placed", "aligned_placed"}
    assert len(rep["mechanism"]["B"]["fgw"]["xcorr_min"]) == 1
    # no cell count, count field, limit or anchor is written or printed; whole numbers are seeds and seed counts only
    keys, ints, texts = set(), set(), []
    _leaves(rep, keys, ints, texts)
    assert not keys & COUNT_KEYS
    assert ints <= {0, 1, BUILD_SEED_OFFSET, dg.TRUTH_SEED_OFFSET}  # seeds, the seed's offsets, n_seeds and seeds_ahead
    forbidden = {*SECRET.values(), *SIZES.values()}
    for text in (*texts, txt, printed):
        masked = re.sub(r"\b[0-9a-f]{7,64}\b", "<hex>", text)  # commits and sha256s hold digit runs by chance
        masked = re.sub(r"\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z?", "<utc>", masked)
        masked = re.sub(r"/\S*", "<path>", masked)  # temporary paths
        found = {int(d) for d in re.findall(r"\d+", masked)}
        assert not found & forbidden, masked
    # outputs stamped: commit, clean tree, utc, the inputs' hashes and a build stamp beside each file
    assert rep["stamp"]["commit"] and rep["stamp"]["tree_clean"] is True and rep["stamp"]["utc"]
    assert len(rep["inputs"]) == 3 and all(len(v) == 64 for v in rep["inputs"].values())
    for name in dg.NAMES:
        assert (out / f"{name}.stamp.json").exists()
    with pytest.raises(SystemExit, match="a rerun needs a new TRIES row"):
        dg.main(["--out", str(out), "--data", str(data), "--seeds", "0"])


def test_main_gives_a_reading_only_on_a_full_run(dg, tmp_path, monkeypatch, own_repo, capsys):
    data = _patch(dg, monkeypatch, tmp_path)
    monkeypatch.setattr(dg, "SEEDS", (0, 1))  # the configured seeds of this synthetic run
    out = tmp_path / "out"
    dg.main(["--out", str(out), "--data", str(data), "--seeds", "0", "1"])
    printed = capsys.readouterr().out
    rep = json.loads((out / "e12d_halves.json").read_text())
    assert rep["n_seeds"] == 2 and len(rep["per_seed"]) == 2
    # a full run of this synthetic pair: tiny worlds, so too few usable seeds for a reading whatever the gains
    assert rep["reading"] == "no clear read"
    for world in dg.WORLDS:
        s = rep["summary"][world]
        usable = [r["worlds"][world] for r in rep["per_seed"] if "degenerate" not in r["worlds"][world]]
        degenerate = [r["worlds"][world] for r in rep["per_seed"] if "degenerate" in r["worlds"][world]]
        assert s["n_seeds"] == len(usable) and s["n_seeds_total"] == 2 and len(s["degenerate"]) == len(degenerate)
        for row in degenerate:  # no scores, but both candidates' raw shape metrics
            assert row["scores"] == {} and "gain" not in row
            assert set(row["raw"]) == {"fgw_shape0", "fgw_shape0_aligned"}
        if usable:
            gains = [r["gain"] for r in usable]
            assert s["gain_mean"] == pytest.approx(np.mean(gains)) and s["seeds_ahead"] == sum(g > 0 for g in gains)
    # the degenerate seeds are named in the report, and nothing else about them
    assert "degenerate" in (out / "e12d_halves.txt").read_text() or all(
        s["degenerate"] == [] for s in rep["summary"].values())
    assert "seed 1" in printed
    # the mechanism read keeps one value per seed for both halves and both couplings
    assert all(len(v["xcorr_min"]) == 2 for half in rep["mechanism"].values() for v in half.values())
