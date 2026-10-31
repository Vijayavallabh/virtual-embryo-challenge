"""vec.harness (feat-008) and vec.registry: pseudo-task rules, boards, paired scoring, decisions, on synthetic data."""

from __future__ import annotations

import csv
from dataclasses import replace

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from conftest import GENES, simulate_stage

from vec import harness, registry
from vec.boards import Board, load_boards
from vec.evaluate import Stage, score_board, subsample
from vec.harness import CANDIDATES, LeakageError, PseudoTask, check_task, load_tasks

# released row counts (docs/DATA_NOTES.md) of each validation board's target and reference
VALIDATION_ROWS = {"T1:val": (24979, 17057), "T2:embryo:val_interp": (11665, 13295),
                   "T2:heart:val_extrap": (167864, 53742), "T2:heart:val_interp": (23710, 58716),
                   "T3:gata4": (46396, 24826)}


def test_config_loads_and_every_task_passes_the_rules():
    tasks = load_tasks()
    assert set(tasks) == {"T1:test", "T2:embryo:test_interp", "T2:heart:test_extrap", "T3:bcatenin"}
    boards = load_boards()
    for test_board, ts in tasks.items():
        assert ts and all(t.test_board == test_board for t in ts)
        for t in ts:
            assert t.panel_board in boards and t.target not in t.inputs.values()
    keys = [t.key for ts in tasks.values() for t in ts]
    assert len(keys) == len(set(keys))


def _task(**over):
    base = PseudoTask(key="x", test_board="T2:heart:test_extrap", task="T2", setting="heart", kind="extrap",
                      target="T2/E10.5.h5ad", reference="T2/E9.5.h5ad",
                      inputs={"ref": "T2/E9.5.h5ad", "start": "T2/E8.75.h5ad"}, twin_board="T2:heart:val_extrap")
    return replace(base, **over)


def test_leakage_rules_refuse_what_could_reveal_the_target():
    check_task(_task())
    with pytest.raises(LeakageError, match="target is an input"):
        check_task(_task(inputs={"ref": "T2/E9.5.h5ad", "peek": "T2/E10.5.h5ad"}))
    with pytest.raises(LeakageError, match="not earlier"):
        check_task(_task(target="T2/E9.5.h5ad", reference="T2/E8.75.h5ad",
                         inputs={"ref": "T2/E8.75.h5ad", "later": "T2/E10.5.h5ad"}))
    with pytest.raises(LeakageError, match="forbidden"):
        check_task(_task(forbidden=("T2/E8.75.h5ad",)))
    interp = dict(task="T2", setting="embryo", kind="interp", target="T2/E7.5.h5ad", reference="T2/E7.25.h5ad",
                  twin_board="T2:embryo:val_interp")
    check_task(_task(**interp, inputs={"before": "T2/E7.25.h5ad", "after": "T2/E8.0.h5ad"}))
    with pytest.raises(LeakageError, match="each side"):
        check_task(_task(**interp, inputs={"before": "T2/E7.25.h5ad", "earlier": "T2/E6.75.h5ad"}))
    with pytest.raises(ValueError, match="earlier bracket"):
        check_task(_task(**interp | {"reference": "T2/E6.75.h5ad"},
                         inputs={"before": "T2/E7.25.h5ad", "after": "T2/E8.0.h5ad", "earlier": "T2/E6.75.h5ad"}))
    with pytest.raises(ValueError, match="unknown files"):
        check_task(_task(inputs={"ref": "T2/E9.5.h5ad", "x": "T2/E12.5.h5ad"}))
    with pytest.raises(ValueError, match="reference must be one of the inputs"):
        check_task(_task(inputs={"start": "T2/E8.75.h5ad"}))
    with pytest.raises(ValueError, match="validation board"):
        check_task(_task(twin_board="T2:heart:test_extrap"))
    with pytest.raises(ValueError, match="knockout"):
        check_task(_task(task="T3", setting=None, kind="ko", target="T2/E9.5.h5ad", reference="T2/E8.75.h5ad",
                         inputs={"wt": "T2/E8.75.h5ad"}, twin_board="T3:gata4"))


def test_count_rule_reproduces_every_validation_board():
    boards = load_boards()
    for key, (n_target, n_ref) in VALIDATION_ROWS.items():
        b = boards[key]
        counts = harness.pseudo_counts(n_target, n_ref)
        assert counts == {"truth_cells": b.truth_cells, "ref_cells": b.ref_cells, "min_cells": b.min_cells,
                          "max_cells": b.max_cells}, key


def test_task_board_keeps_validation_specs_and_builds_synthetic_ones():
    boards = load_boards()
    tasks = {t.key: t for ts in load_tasks().values() for t in ts}
    val = harness.task_board(tasks["T2h/<=E9.5->E10.5"], boards, *VALIDATION_ROWS["T2:heart:val_extrap"])
    assert val.key == "T2:heart:val_extrap" and val.anchors == {}
    with pytest.raises(ValueError, match="published counts"):
        harness.task_board(tasks["T2h/<=E9.5->E10.5"], boards, 1000, 53742)
    syn = harness.task_board(tasks["T2h/<=E8.75->E10.5"], boards, 167864, 24826)
    assert syn.key == "pseudo:T2h/<=E8.75->E10.5" and syn.genes == boards["T2:heart:val_extrap"].genes
    assert (syn.truth_cells, syn.ref_cells, syn.max_cells, syn.min_cells) == (8393, 2483, 25179, 1000)
    leaky = dict(boards) | {"T3:bcatenin": replace(boards["T3:gata4"], key="T3:bcatenin",
                                                   anchors={"mmd_u": {"floor": 0.1, "ceiling": 0.01}})}
    with pytest.raises(LeakageError):
        harness.check_index(leaky)


def test_candidates_apply_only_where_their_inputs_exist():
    tasks = {t.key: t for ts in load_tasks().values() for t in ts}
    pbs, st = CANDIDATES["pseudobulk_shift"], CANDIDATES["shift_transfer"]
    assert pbs.applies(tasks["T2h/<=E9.5->E10.5"]) and pbs.applies(tasks["T1/E9.5->E10.5"])
    assert pbs.applies(tasks["T2e/E7.25+E8.0->E7.5"])  # earlier bracket shifted from E6.75, as the organisers do
    assert not pbs.applies(tasks["T2e/E6.75+E7.5->E7.25"]) and not pbs.applies(tasks["T1/E8.5->E9.5"])
    assert not pbs.applies(tasks["T3/Gata4rep1->rep2"]) and st.applies(tasks["T3/Gata4->Mab21l2"])
    assert all(CANDIDATES["copy"].applies(t) for t in tasks.values())


# ------------------------------------------------------------------------------------------------ synthetic scoring
def _pool(n, rng, means, props, shift=0.0, scale=1.0, stretch=(1.0, 1.0, 1.0)):
    X, ct, coords = _sim(n, rng, means, props, shift, scale)
    coords = (coords * np.asarray(stretch, np.float32)).astype(np.float32)
    a = ad.AnnData(X=X, obs=pd.DataFrame({"celltype": ct}, index=[f"c{i}" for i in range(n)]),
                   var=pd.DataFrame(index=GENES))
    a.obsm["spatial_3D"] = coords
    return a


def _sim(n, rng, means, props, shift, scale):
    centres = np.array([[0, 0, 0], [6, 0, 0], [0, 6, 0]], float) * scale + shift
    st = simulate_stage(n, means, props, rng, coords_scale=1.5 * scale, coords_shift=centres)
    return st.X, st.celltype, st.coords


@pytest.fixture(scope="module")
def synthetic_t2():
    rng = np.random.default_rng(11)
    base = rng.gamma(2.0, 2.0, size=(3, len(GENES)))
    later = base.copy()
    later[:, :6] *= 3.0
    start = _pool(3000, rng, base, [0.5, 0.3, 0.2])
    ref = _pool(4000, rng, base * 1.2, [0.4, 0.35, 0.25], scale=1.1)
    target = _pool(6000, rng, later, [0.3, 0.4, 0.3], scale=1.3, stretch=(1.0, 1.0, 3.0))
    board = Board(key="pseudo:synthetic", task="T2", genes=list(GENES), needs_coords=True, setting="heart",
                  min_cells=100, max_cells=900, truth_cells=300, ref_cells=400)
    task = PseudoTask(key="syn", test_board="T2:heart:test_extrap", task="T2", setting="heart", kind="extrap",
                      target="T2/E10.5.h5ad", reference="T2/E9.5.h5ad",
                      inputs={"ref": "T2/E9.5.h5ad", "start": "T2/E8.75.h5ad"}, twin_board="T2:heart:val_extrap")

    def st(a):
        return Stage(X=np.asarray(a.X, np.float32), coords=np.asarray(a.obsm["spatial_3D"], np.float32),
                     celltype=np.asarray(a.obs["celltype"]))

    return task, board, st(target), st(ref), {"ref": ref, "start": start}


def test_an_anchor_map_scores_exactly_like_computed_anchors(synthetic_t2):
    task, board, target, ref, inputs = synthetic_t2
    pred = subsample(target, 0.2, 5)
    a = score_board(board, pred, target, ref, seed=3, anchors="computed", protocol="bake", build="main")
    b = score_board(board, pred, target, ref, seed=3, anchors=a.anchors, protocol="bake", build="main")
    assert a.score == b.score and a.skills == b.skills


def test_paired_scoring_floor_common_draws_and_build_seeds(synthetic_t2):
    task, board, target, ref, inputs = synthetic_t2
    twin = replace(CANDIDATES["copy"], name="copy_twin")
    cands = [CANDIDATES["floor"], CANDIDATES["copy"], twin, CANDIDATES["pseudobulk_shift"]]
    rows = harness.score_seeds(task, board, target, ref, inputs, cands, seeds=(0, 1, 2))
    for row in rows:
        assert row["scores"]["floor"] == 50.0
        assert row["scores"]["copy"] == row["scores"]["copy_twin"]  # same build seed, same draws, same score
        assert row["n_cells"]["copy"] == row["n_cells"]["pseudobulk_shift"] == board.max_cells
        assert row["build_seed"] == harness.BUILD_SEED_OFFSET + row["seed"]
    # the candidates' rows never reproduce the scorer's reference draw at the same seed
    ctx = harness.BuildContext(task=task, board=board, inputs=inputs, seed=harness.BUILD_SEED_OFFSET)
    _, rows0, _, _ = harness.draw_copy(ctx)
    floor_rows = np.random.default_rng(0).choice(inputs["ref"].n_obs, 400, replace=False)
    assert len(set(floor_rows) & set(rows0.tolist())) < 0.6 * len(floor_rows)
    # why the offset exists: numpy's [s, 0] stream (sample_cells' rows) is the scoring stream s itself
    assert np.array_equal(np.random.default_rng([0, 0]).random(8), np.random.default_rng(0).random(8))
    for s in harness.SEEDS:
        b = harness.BUILD_SEED_OFFSET + s
        assert b not in harness.SEEDS
        for k in harness.SEEDS:
            assert not np.array_equal(np.random.default_rng([b, 0]).random(8), np.random.default_rng(k).random(8))


def test_decide_pairs_only_tasks_where_both_candidates_ran():
    gains = [1.0, 0.9, 1.2, 0.8, 1.1, 0.7, 1.3, 0.9, 1.0, 0.6]
    r1 = {"task": {"key": "a"}, "candidates": ["floor", "copy"],
          "rows": [{"seed": s, "scores": {"floor": 50.0, "copy": 50.0 + g}} for s, g in enumerate(gains)]}
    r2 = {"task": {"key": "b"}, "candidates": ["floor"], "rows": [{"seed": s, "scores": {"floor": 50.0}}
                                                                   for s in range(10)]}
    d = harness.decide([r1, r2], "floor", "copy")
    assert d["pseudo_tasks"] == ["a"] and d["keep"] and d["need"] == 1
    assert not harness.decide([r2], "floor", "copy")["keep"]


def test_run_task_gives_candidates_only_the_listed_inputs(synthetic_t2, monkeypatch, tmp_path):
    task, board, target, ref, inputs = synthetic_t2
    loaded_inputs, scored = [], []
    monkeypatch.setattr(harness, "read_stage", lambda path, genes: scored.append(str(path)) or (
        ad.AnnData(X=np.asarray(target.X), var=pd.DataFrame(index=GENES), obsm={"spatial_3D": target.coords})
        if str(path).endswith(task.target) else
        ad.AnnData(X=np.asarray(ref.X), var=pd.DataFrame(index=GENES), obsm={"spatial_3D": ref.coords})))
    by_file = {f: inputs[r] for r, f in task.inputs.items()}
    monkeypatch.setattr(harness, "load_stage", lambda path, **kw: loaded_inputs.append(str(path)) or
                        by_file[str(path).split(str(tmp_path) + "/")[1]])
    monkeypatch.setattr(harness, "_sha256", lambda path: "0" * 64)
    monkeypatch.setattr(harness, "task_board", lambda t, boards, n_t, n_r: board)
    boards = {"T2:heart:val_extrap": replace(board, key="T2:heart:val_extrap")}
    res = harness.run_task(task, [CANDIDATES["floor"], CANDIDATES["copy"]], seeds=(0,), data_dir=tmp_path,
                           boards=boards, log=None)
    assert all(not p.endswith(task.target) for p in loaded_inputs)
    assert sorted(p.split(str(tmp_path) + "/")[1] for p in loaded_inputs) == sorted(task.inputs.values())
    assert res["candidates"] == ["floor", "copy"] and res["protocol"] == "bake" and res["build"] == "main"


def test_registry_header_and_refusals(tmp_path):
    path = tmp_path / "registry.csv"
    row = dict.fromkeys(registry.HEADER, "") | {"candidate": "copy", "protocol": "bake"}
    assert registry.append([row, row], path=path) == 2
    assert registry.append([row], path=path) == 1
    lines = list(csv.reader(path.open()))
    assert tuple(lines[0]) == registry.HEADER and len(lines) == 4
    with pytest.raises(ValueError, match="differ from the header"):
        registry.append([{"candidate": "x"}], path=path)
    other = tmp_path / "old.csv"
    other.write_text("utc_time,commit,board\n")
    with pytest.raises(ValueError, match="another header"):
        registry.append([row], path=other)
