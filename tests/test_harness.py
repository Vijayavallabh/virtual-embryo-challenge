"""vec.harness (feat-008) and vec.registry: pseudo-task rules, boards, paired scoring, the read guard, decisions, on
synthetic data. Revised after the 2026-11-01 four-lens review: every safeguard has a test that fails without it."""

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


# ------------------------------------------------------------------------------------------------ config and rules
def test_config_loads_and_every_task_passes_the_rules():
    tasks = load_tasks()
    assert set(tasks) == {"T1:test", "T2:embryo:test_interp", "T2:heart:test_extrap", "T3:bcatenin"}
    boards = load_boards()
    for test_board, ts in tasks.items():
        assert ts and all(t.test_board == test_board for t in ts)
        for t in ts:
            assert check_task(t) is t and t.panel_board in boards and t.target not in t.inputs.values()
    keys = [t.key for ts in tasks.values() for t in ts]
    assert len(keys) == len(set(keys))
    t3 = {t.key: t for t in tasks["T3:bcatenin"]}
    assert t3["T3/Gata4rep1->rep2"].replicate and not t3["T3/Gata4rep1->rep2"].counted


def test_load_tasks_refuses_a_leaky_or_duplicated_config(tmp_path):
    entry = ("  - {key: k, task: T2, setting: heart, kind: extrap, target: T2/E10.5.h5ad, reference: T2/E9.5.h5ad, "
             "inputs: {ref: T2/E9.5.h5ad, peek: %s}, twin_board: 'T2:heart:val_extrap'}\n")
    leaky = tmp_path / "leaky.yaml"
    leaky.write_text("T2:heart:test_extrap:\n" + entry % "T2/E10.5.h5ad")
    with pytest.raises(LeakageError):
        load_tasks(leaky)
    dup = tmp_path / "dup.yaml"
    dup.write_text("T2:heart:test_extrap:\n" + entry % "T2/E8.75.h5ad" + entry % "T2/E8.5.h5ad")
    with pytest.raises(ValueError, match="duplicate"):
        load_tasks(dup)


def _task(**over):
    base = PseudoTask(key="x", test_board="T2:heart:test_extrap", task="T2", setting="heart", kind="extrap",
                      target="T2/E10.5.h5ad", reference="T2/E9.5.h5ad",
                      inputs={"ref": "T2/E9.5.h5ad", "start": "T2/E8.75.h5ad"}, twin_board="T2:heart:val_extrap")
    return replace(base, **over)


INTERP = dict(task="T2", setting="embryo", kind="interp", target="T2/E7.5.h5ad", reference="T2/E7.25.h5ad",
              twin_board="T2:embryo:val_interp")
KO = dict(task="T3", setting=None, kind="ko", target="T3/E8.75_gata4_ko_2.h5ad", reference="T2/E8.75.h5ad",
          twin_board="T3:gata4")


@pytest.mark.parametrize("over, err, match", [
    (dict(inputs={"ref": "T2/E9.5.h5ad", "peek": "T2/E10.5.h5ad"}), LeakageError, "target is an input"),
    (dict(target="T2/E9.5.h5ad", reference="T2/E8.75.h5ad", inputs={"ref": "T2/E8.75.h5ad", "x": "T2/E10.5.h5ad"}),
     LeakageError, "later than the reference"),
    # another modality or genotype on the target's own day is not earlier than the target
    (dict(target="T2/E9.5.h5ad", reference="T2/E8.75.h5ad",
          inputs={"ref": "T2/E8.75.h5ad", "ko": "T3/E9.5_mab21l2_ko.h5ad"}), LeakageError, "later than the reference"),
    # a skipped stage cannot come back through another file: the reference is the latest stage allowed
    (dict(reference="T2/E8.75.h5ad", inputs={"ref": "T2/E8.75.h5ad", "t1": "T1/E9.5_RNA.h5ad"}), LeakageError,
     "later than the reference"),
    (dict(forbidden=("T2/E8.75.h5ad",)), LeakageError, "forbidden"),
    (dict(inputs={"ref": "T2/E9.5.h5ad", "x": "T2/E12.5.h5ad"}), ValueError, "unknown files"),
    (dict(inputs={"start": "T2/E8.75.h5ad"}), ValueError, "reference must be one of the inputs"),
    (dict(twin_board="T2:heart:test_extrap"), ValueError, "validation board"),
    (dict(validation_board="T2:heart:test_extrap", twin_board=None), ValueError, "validation board"),
    (dict(kind="extrapolation"), ValueError, "kind one of"),
    (dict(task="T4"), ValueError, "task must be one of"),
    (dict(setting=None), ValueError, "setting"),
    (dict(**INTERP, inputs={"before": "T2/E7.25.h5ad", "earlier": "T2/E6.75.h5ad"}), ValueError, "'after'"),
    (dict(**INTERP, inputs={"before": "T2/E7.25.h5ad", "after": "T2/E8.0.h5ad", "x": "T3/E8.75_gata4_ko_1.h5ad"}),
     None, None),
    (dict(**INTERP | {"target": "T2/E8.75.h5ad", "reference": "T2/E8.5.h5ad"},
          inputs={"before": "T2/E8.5.h5ad", "after": "T2/E9.5.h5ad", "ko": "T3/E8.75_gata4_ko_1.h5ad"}),
     LeakageError, "inside the bracket"),
    (dict(**INTERP | {"reference": "T2/E6.75.h5ad"},
          inputs={"before": "T2/E7.25.h5ad", "after": "T2/E8.0.h5ad", "earlier": "T2/E6.75.h5ad"}),
     ValueError, "'before' \\(the reference\\)"),
    (dict(**KO, inputs={"wt": "T2/E8.75.h5ad", "ko_train": "T3/E8.75_gata4_ko_2.h5ad"}), LeakageError,
     "target is an input"),
    (dict(**KO, inputs={"wt": "T2/E8.75.h5ad", "ko_train": "T3/E8.75_gata4_ko_1.h5ad"}), LeakageError,
     "own genotype"),
    (dict(**KO | {"reference": "T3/E9.5_mab21l2_ko.h5ad"}, inputs={"wt": "T3/E9.5_mab21l2_ko.h5ad"}), ValueError,
     "wild-type reference"),
    (dict(replicate=True), ValueError, "replicate fold"),
    (dict(target_filter={"keep": "chamber"}), ValueError, "target filter"),
])
def test_leakage_rules(over, err, match):
    t = _task(**over)
    if err is None:
        assert check_task(t) is t
    else:
        with pytest.raises(err, match=match):
            check_task(t)


def test_a_replicate_fold_may_train_on_the_sibling_and_never_counts():
    t = _task(**KO, inputs={"wt": "T2/E8.75.h5ad", "ko_train": "T3/E8.75_gata4_ko_1.h5ad"}, replicate=True)
    assert check_task(t) is t and not t.counted


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
    assert syn.anchors == {}
    mab = harness.task_board(tasks["T3/Gata4->Mab21l2"], boards, 50294, 53742)
    assert mab.target_ko == "mab21l2" and mab.anchors == {}
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
    assert not pbs.applies(_task(inputs={"ref": "T2/E9.5.h5ad", "earlier": "T2/E8.5.h5ad"}))  # extrap needs start
    assert not pbs.applies(tasks["T3/Gata4rep1->rep2"]) and st.applies(tasks["T3/Gata4->Mab21l2"])
    assert all(CANDIDATES["copy"].applies(t) for t in tasks.values())


def test_external_data_is_refused_inside_the_window():
    t = _task()  # E10.5 from E9.5: window (8.5, 11.5]
    ok = replace(CANDIDATES["copy"], externals=(("X002", 14.5, 18.5), ("X003", 8.0, 8.5)))
    harness.check_externals(ok, t)
    bad = replace(CANDIDATES["copy"], externals=(("X009", 11.0, 12.0),))
    with pytest.raises(LeakageError, match="window"):
        harness.check_externals(bad, t)
    interp = _task(**INTERP, inputs={"before": "T2/E7.25.h5ad", "after": "T2/E8.0.h5ad"})
    with pytest.raises(LeakageError):
        harness.check_externals(replace(CANDIDATES["copy"], externals=(("Y", 7.5, 7.5),)), interp)
    harness.check_externals(replace(CANDIDATES["copy"], externals=(("Y", 8.0, 9.0),)), interp)


# ------------------------------------------------------------------------------------------------ synthetic scoring
def _pool(n, rng, means, props, shift=0.0, scale=1.0, stretch=(1.0, 1.0, 1.0)):
    centres = np.array([[0, 0, 0], [6, 0, 0], [0, 6, 0]], float) * scale + shift
    st = simulate_stage(n, means, props, rng, coords_scale=1.5 * scale, coords_shift=centres)
    a = ad.AnnData(X=st.X, obs=pd.DataFrame({"celltype": st.celltype}, index=[f"c{i}" for i in range(n)]),
                   var=pd.DataFrame(index=GENES))
    a.obsm["spatial_3D"] = (st.coords * np.asarray(stretch, np.float32)).astype(np.float32)
    return a


def _stage(a):
    return Stage(X=np.asarray(a.X, np.float32), coords=np.asarray(a.obsm["spatial_3D"], np.float32),
                 celltype=np.asarray(a.obs["celltype"]))


@pytest.fixture(scope="module")
def synthetic_t2():
    rng = np.random.default_rng(11)
    base = rng.gamma(2.0, 2.0, size=(3, len(GENES)))
    later = base.copy()
    later[:, :6] *= 3.0
    start = _pool(3000, rng, base, [0.5, 0.3, 0.2])
    ref = _pool(12000, rng, base * 1.2, [0.4, 0.35, 0.25], scale=1.1)  # > 10,000 rows, like the real pools
    target = _pool(6000, rng, later, [0.3, 0.4, 0.3], scale=1.3, stretch=(1.0, 1.0, 3.0))
    board = Board(key="pseudo:synthetic", task="T2", genes=list(GENES), needs_coords=True, setting="heart",
                  min_cells=100, max_cells=3600, truth_cells=300, ref_cells=1200)
    task = _task(key="syn")
    return task, board, _stage(target), _stage(ref), {"ref": ref, "start": start}


def test_an_anchor_map_scores_exactly_like_computed_anchors_and_must_be_complete(synthetic_t2):
    task, board, target, ref, inputs = synthetic_t2
    pred = subsample(target, 0.2, 5)
    a = score_board(board, pred, target, ref, seed=3, anchors="computed", protocol="bake", build="main")
    b = score_board(board, pred, target, ref, seed=3, anchors=a.anchors, protocol="bake", build="main")
    assert a.score == b.score and a.skills == b.skills
    with pytest.raises(ValueError, match="anchor map"):
        score_board(board, pred, target, ref, seed=3, anchors={"mmd_u": {"floor": float("nan"), "ceiling": 0.1}},
                    protocol="bake")


def _spy(log):
    def build(ctx):
        log.append({"seed": ctx.seed, "roles": set(ctx.inputs), "board_truth": ctx.board.truth_cells,
                    "view_has_paths": any("/" in str(v) and ".h5ad" in str(v) for v in vars(ctx.task).values())})
        pool, rows, X, coords = harness.draw_copy(ctx)
        log[-1]["rows"] = rows
        pred = harness.finish(ctx, pool, X, coords)
        log[-1]["pred"] = pred
        return pred
    return harness.Candidate("spy", frozenset(harness.TASKS), harness._ANY, build)


def test_paired_scoring_on_shared_draws_with_offset_build_seeds(synthetic_t2):
    task, board, target, ref, inputs = synthetic_t2
    log = []
    twin = replace(CANDIDATES["copy"], name="copy_twin")
    cands = [CANDIDATES["floor"], CANDIDATES["copy"], twin, CANDIDATES["pseudobulk_shift"], _spy(log)]
    seeds = (0, 1, 2)
    rows = harness.score_seeds(task, board, target, ref, inputs, cands, seeds=seeds)
    assert [e["seed"] for e in log] == [harness.BUILD_SEED_OFFSET + s for s in seeds]
    for row, e in zip(rows, log, strict=True):
        s = row["seed"]
        assert row["scores"]["floor"] == 50.0
        assert row["scores"]["copy"] == row["scores"]["copy_twin"] == row["scores"]["spy"]  # common random numbers
        assert row["n_cells"]["copy"] == row["n_cells"]["pseudobulk_shift"] == board.max_cells
        # the harness score is the score the calibrated emulator gives this file on its own
        direct = score_board(board, e["pred"], target, ref, seed=s, anchors="computed", protocol="bake", build="main")
        assert row["scores"]["spy"] == direct.score
        # the row keeps each candidate's metric panel, de_score's audit trail included, and the floor's no-change zero
        assert row["metrics"]["spy"] == direct.metrics and set(row["metrics"]) == set(row["scores"])
        assert row["metrics"]["floor"]["de_score"] == 0.0 and "_de_chance" in row["metrics"]["copy"]
        # candidates see their inputs only, no truth count, no file paths
        assert e["roles"] == set(task.inputs) and e["board_truth"] is None and not e["view_has_paths"]
        # the build rows are not the scorer's reference draw (at the scoring seed they would contain it)
        floor_rows = np.random.default_rng(s).choice(12000, 1200, replace=False)
        assert len(set(floor_rows) & set(e["rows"].tolist())) < 0.6 * len(floor_rows)
        same = harness.draw_copy(harness.BuildContext(task=harness.TaskView.of(task), board=board, inputs=inputs,
                                                      seed=s))[1]
        assert set(floor_rows) <= set(same.tolist())
    assert not {harness.BUILD_SEED_OFFSET + k for k in harness.SEEDS} & set(harness.SEEDS)


def test_the_floor_check_refuses_a_floor_that_does_not_score_fifty(synthetic_t2, monkeypatch):
    task, board, target, ref, inputs = synthetic_t2
    real = harness.score_board

    def off(board_, pred, *a, **kw):
        r = real(board_, pred, *a, **kw)
        return replace(r, score=49.0) if kw.get("anchors") == "computed" else r

    monkeypatch.setattr(harness, "score_board", off)
    with pytest.raises(ValueError, match="the floor scores"):
        harness.score_seeds(task, board, target, ref, inputs, [CANDIDATES["copy"]], seeds=(0,))


def test_a_floor_working_copy_built_as_a_candidate_scores_fifty(synthetic_t2):
    task, board, target, ref, inputs = synthetic_t2
    cand = harness.Candidate("floor_copy", frozenset(harness.TASKS), harness._ANY,
                             lambda ctx: subsample(ref, 0.1, ctx.seed - harness.BUILD_SEED_OFFSET))
    rows = harness.score_seeds(task, board, target, ref, inputs, [cand], seeds=(4,))
    assert rows[0]["scores"]["floor_copy"] == 50.0


def test_the_read_guard_refuses_target_derived_files(synthetic_t2, tmp_path):
    task, board, target, ref, inputs = synthetic_t2
    (tmp_path / "experiments" / "coords" / "results").mkdir(parents=True)
    (tmp_path / "experiments" / "coords" / "results" / "size.json").write_text("{}")
    (tmp_path / "configs").mkdir()
    (tmp_path / "configs" / "ok.yaml").write_text("a: 1\n")
    data = tmp_path / "data" / "raw"
    data.mkdir(parents=True)
    (data / "T2").mkdir()
    ad.AnnData(X=np.zeros((2, 2), np.float32)).write_h5ad(data / "T2" / "E10.5.h5ad")

    def reads(path):
        def build(ctx):
            if str(path).endswith(".h5ad"):
                ad.read_h5ad(path)
            else:
                open(path).read()
            return harness._copy(ctx)
        return harness.Candidate("reader", frozenset(harness.TASKS), harness._ANY, build)

    kw = dict(seeds=(0,), root=tmp_path, data_dir=data)
    with pytest.raises(PermissionError, match="refused"):
        harness.score_seeds(task, board, target, ref, inputs,
                            [reads(tmp_path / "experiments" / "coords" / "results" / "size.json")], **kw)
    with pytest.raises(PermissionError, match="refused"):
        harness.score_seeds(task, board, target, ref, inputs, [reads(data / "T2" / "E10.5.h5ad")], **kw)
    rows = harness.score_seeds(task, board, target, ref, inputs, [reads(tmp_path / "configs" / "ok.yaml")], **kw)
    assert rows[0]["reads"]["reader"] == ["configs/ok.yaml"]
    open(tmp_path / "experiments" / "coords" / "results" / "size.json").read()  # the guard is off outside builds


def _fake_io(monkeypatch, tmp_path, task, target, ref, inputs, calls):
    def read(path, genes):
        calls["scored"].append(str(path).split(str(tmp_path) + "/")[1])
        st = target if str(path).endswith(task.target) else ref
        return ad.AnnData(X=np.asarray(st.X), var=pd.DataFrame(index=GENES),
                          obsm={"spatial_3D": st.coords} if st.coords is not None else {})

    by_file = {f: inputs[r] for r, f in task.inputs.items()}

    def load(path, **kw):
        calls["inputs"].append((str(path).split(str(tmp_path) + "/")[1], kw))
        return by_file[str(path).split(str(tmp_path) + "/")[1]]

    monkeypatch.setattr(harness, "read_stage", read)
    monkeypatch.setattr(harness, "load_stage", load)
    monkeypatch.setattr(harness, "_sha256", lambda path: "0" * 64)


def test_run_task_loads_target_and_reference_as_released_and_inputs_only(synthetic_t2, monkeypatch, tmp_path):
    task, board, target, ref, inputs = synthetic_t2
    calls = {"scored": [], "inputs": []}
    _fake_io(monkeypatch, tmp_path, task, target, ref, inputs, calls)
    monkeypatch.setattr(harness, "task_board", lambda t, boards, n_t, n_r: board)
    boards = {"T2:heart:val_extrap": replace(board, key="T2:heart:val_extrap")}
    res = harness.run_task(task, [CANDIDATES["copy"]], seeds=(0,), data_dir=tmp_path, boards=boards, log=None)
    assert calls["scored"] == [task.target, task.reference]
    assert sorted(f for f, _ in calls["inputs"]) == sorted(task.inputs.values())
    assert all(kw.get("drop_duplicates") is True for _, kw in calls["inputs"])  # T2: duplicates dropped from inputs
    assert res["candidates"] == ["floor", "copy"] and res["protocol"] == "bake" and res["build"] == "main"
    # an anchored non-validation board stops the run before anything is read
    calls = {"scored": [], "inputs": []}
    _fake_io(monkeypatch, tmp_path, task, target, ref, inputs, calls)
    leaky = boards | {"T2:heart:test_extrap": replace(board, key="T2:heart:test_extrap",
                                                      anchors={"mmd_u": {"floor": 0.1, "ceiling": 0.01}})}
    with pytest.raises(LeakageError):
        harness.run_task(task, [CANDIDATES["copy"]], seeds=(0,), data_dir=tmp_path, boards=leaky, log=None)
    assert calls == {"scored": [], "inputs": []}
    # a leaky task is refused before anything is read
    with pytest.raises(LeakageError):
        harness.run_task(replace(task, inputs={"ref": task.reference, "x": task.target}), [CANDIDATES["copy"]],
                         seeds=(0,), data_dir=tmp_path, boards=boards, log=None)
    assert calls == {"scored": [], "inputs": []}


def test_run_task_on_t1_loads_inputs_as_released(monkeypatch, tmp_path):
    rng = np.random.default_rng(5)
    base = rng.gamma(2.0, 2.0, size=(3, len(GENES)))

    def t1(n, m):
        st = simulate_stage(n, m, [0.5, 0.3, 0.2], rng)
        return ad.AnnData(X=st.X, obs=pd.DataFrame({"celltype": st.celltype}, index=[f"c{i}" for i in range(n)]),
                          var=pd.DataFrame(index=GENES))

    later, earlier = base.copy(), base.copy()
    later[:, :8] *= 4.0  # gene-specific changes: a uniform scale cancels in log1p CP10k and leaves no DE genes
    later[:, 8:16] *= 0.25
    earlier[:, 16:24] *= 3.0
    ref_a, start_a, tgt_a = t1(3000, base), t1(2500, earlier), t1(4000, later)
    task = PseudoTask(key="t1syn", test_board="T1:test", task="T1", kind="extrap", target="T1/E10.5_RNA.h5ad",
                      reference="T1/E9.5_RNA.h5ad", inputs={"ref": "T1/E9.5_RNA.h5ad", "start": "T1/E8.5_RNA.h5ad"},
                      twin_board="T1:val")
    board = Board(key="pseudo:t1syn", task="T1", genes=list(GENES), needs_coords=False, min_cells=100, max_cells=600,
                  truth_cells=200, ref_cells=300)
    calls = {"scored": [], "inputs": []}
    st = {"target": Stage(X=np.asarray(tgt_a.X)), "ref": Stage(X=np.asarray(ref_a.X))}
    _fake_io(monkeypatch, tmp_path, task, st["target"], st["ref"], {"ref": ref_a, "start": start_a}, calls)
    monkeypatch.setattr(harness, "task_board", lambda t, boards, n_t, n_r: board)
    res = harness.run_task(task, [CANDIDATES["copy"], CANDIDATES["pseudobulk_shift"]], seeds=(0,),
                           data_dir=tmp_path, boards={"T1:val": replace(board, key="T1:val")}, log=None)
    assert calls["scored"] == [task.target, task.reference]
    assert all(kw == {"lineage": False} for _, kw in calls["inputs"])
    assert res["rows"][0]["scores"]["floor"] == 50.0
    assert set(res["candidates"]) == {"floor", "copy", "pseudobulk_shift"}


def test_the_held_flat_list_is_handed_over_only_where_its_stages_are_inputs(monkeypatch):
    from vec import t1_scope

    d = harness.DERIVED["t1_held_flat"]
    assert d.path == t1_scope.HELD_FLAT_GENES and d.sources() == ("T1/E8.5_RNA.h5ad", "T1/E9.5_RNA.h5ad")
    tasks = {t.key: t for ts in load_tasks().values() for t in ts}
    b0 = CANDIDATES["b0"]
    assert b0.applies(tasks["T1/E9.5->E10.5"]) and b0.applies(tasks["T1/E9.5->E10.5:chamber"])
    assert not b0.applies(tasks["T1/E8.5->E9.5"])
    assert b0.applies(tasks["T2h/<=E9.5->E10.5"]) and harness.load_derived(b0, tasks["T2h/<=E9.5->E10.5"]) == {}
    held = harness.load_derived(b0, tasks["T1/E9.5->E10.5"])["t1_held_flat"]
    assert list(held) == pd.read_csv(t1_scope.HELD_FLAT_GENES)["gene"].astype(str).tolist() and "Xist" in held
    with harness.guarded_reads(), pytest.raises(PermissionError, match="refused"):
        t1_scope.held_flat_genes()  # a candidate reading the committed file itself is refused
    # a list computed on the target's stage would never be handed over: the candidate does not apply there
    monkeypatch.setattr(harness, "DERIVED", {"t1_held_flat": replace(
        d, sources=lambda: ("T1/E9.5_RNA.h5ad", "T1/E10.5_RNA.h5ad"))})
    assert not b0.applies(tasks["T1/E9.5->E10.5"]) and b0.applies(tasks["T2h/<=E9.5->E10.5"])
    with pytest.raises(LeakageError, match="derived input"):
        harness.load_derived(b0, tasks["T1/E9.5->E10.5"])


def test_b0_on_t1_builds_inside_the_read_guard_with_the_list_handed_over(monkeypatch, tmp_path):
    from scipy import sparse

    from vec import t1_scope
    from vec.prep import log1p_cp10k

    rng = np.random.default_rng(7)
    base = rng.gamma(2.0, 2.0, size=(3, len(GENES)))

    def t1(n, m, labels):  # released T1 labels, so the chamber config classifies them (E9.5 and E8.5 per-file rules)
        types = rng.choice(3, size=n, p=[0.5, 0.3, 0.2])
        C = rng.poisson(m[types] * 20).astype(np.int32)
        C[np.arange(n), rng.integers(0, len(GENES), n)] += 1
        a = ad.AnnData(X=log1p_cp10k(C), obs=pd.DataFrame({"celltype": np.asarray(labels)[types]},
                                                          index=[f"c{i}" for i in range(n)]),
                       var=pd.DataFrame(index=GENES))
        a.layers["counts"] = sparse.csr_matrix(C)
        return a

    later, earlier = base.copy(), base.copy()
    later[:, :8] *= 4.0
    later[:, 8:16] *= 0.25
    earlier[:, 16:24] *= 3.0
    ref_a = t1(3000, base, ["V-CM", "Endocardium", "Foregut"])
    start_a = t1(2500, earlier, ["V-CM", "Endothelium", "Foregut"])
    tgt_a = t1(4000, later, ["V-CM", "Endocardium", "Foregut"])
    task = PseudoTask(key="t1b0", test_board="T1:test", task="T1", kind="extrap", target="T1/E10.5_RNA.h5ad",
                      reference="T1/E9.5_RNA.h5ad", inputs={"ref": "T1/E9.5_RNA.h5ad", "start": "T1/E8.5_RNA.h5ad"},
                      twin_board="T1:val")
    board = Board(key="pseudo:t1b0", task="T1", genes=list(GENES), needs_coords=False, min_cells=100, max_cells=600,
                  truth_cells=200, ref_cells=300)
    calls = {"scored": [], "inputs": []}
    st = {"target": Stage(X=tgt_a.X.toarray()), "ref": Stage(X=ref_a.X.toarray())}
    _fake_io(monkeypatch, tmp_path, task, st["target"], st["ref"], {"ref": ref_a, "start": start_a}, calls)
    monkeypatch.setattr(harness, "task_board", lambda t, boards, n_t, n_r: board)
    held_seen, real = [], t1_scope.build_submission

    def spy(*a, **kw):
        held_seen.append(kw.get("held_flat"))
        return real(*a, **kw)

    monkeypatch.setattr(t1_scope, "build_submission", spy)
    res = harness.run_task(task, [CANDIDATES["b0"], CANDIDATES["b0_scope"]], seeds=(0,), data_dir=tmp_path,
                           boards={"T1:val": replace(board, key="T1:val")}, log=None)
    assert res["candidates"] == ["floor", "b0", "b0_scope"] and res["rows"][0]["scores"]["floor"] == 50.0
    want = pd.read_csv(t1_scope.HELD_FLAT_GENES)["gene"].astype(str).tolist()
    assert held_seen == [want, None]  # b0 shifts with the handed-over list; b0_scope does not shift
    assert not any("results" in p for p in res["rows"][0]["reads"]["b0"])
    assert res["params"]["b0"]["derived"] == {"t1_held_flat": {
        "file": "experiments/t1_direction/results/held_flat_genes.csv",
        "sources": ["T1/E8.5_RNA.h5ad", "T1/E9.5_RNA.h5ad"]}}
    assert "derived" not in res["params"]["b0_scope"] and res["provenance"]["derived_sha256"] == {
        "t1_held_flat": "0" * 64}


# ------------------------------------------------------------------------------------------------ deciding
GAINS = [1.0, 0.9, 1.2, 0.8, 1.1, 0.7, 1.3, 0.9, 1.0, 0.6]


def _result(key, cands=("floor", "copy"), gains=GAINS, **prov):
    rows = [{"seed": s, "build_seed": harness.BUILD_SEED_OFFSET + s,
             "scores": {"floor": 50.0, **({"copy": 50.0 + g} if "copy" in cands else {})}}
            for s, g in zip(harness.SEEDS, gains, strict=True)]
    p = {"commit": "abc", "tree_clean": True, "config_sha256": "c" * 64, "files_sha256": {"T2/E9.5.h5ad": "f"}} | prov
    return {"task": {"key": key, "test_board": "T2:heart:test_extrap"}, "candidates": list(cands), "rows": rows,
            "protocol": "bake", "build": "main", "anchors": "computed", "frac": 0.1, "seeds": list(harness.SEEDS),
            "provenance": p}


CONFIGURED = [_task(key="a"), _task(key="b")]


def test_decide_counts_the_configured_pseudo_tasks():
    d = harness.decide([_result("a"), _result("b")], CONFIGURED, "floor", "copy")
    assert d["keep"] and d["counted"] == ["a", "b"] and d["need"] == 2
    missing = harness.decide([_result("a")], CONFIGURED, "floor", "copy")  # b has no result: incomplete
    assert not missing["keep"] and "incomplete" in missing["reason"]
    unscored = harness.decide([_result("a"), _result("b", cands=("floor",))], CONFIGURED, "floor", "copy")
    assert not unscored["keep"]
    diag = harness.decide([_result("a")], [CONFIGURED[0], replace(CONFIGURED[1], adoption=False)], "floor", "copy")
    assert diag["keep"] and diag["excluded"] == {"b": "diagnostic"} and diag["decided_on"].startswith("1 of 2")
    tuned = harness.decide([_result("a"), _result("b")], CONFIGURED, "floor", "copy",
                           candidates=CANDIDATES | {"copy": replace(CANDIDATES["copy"], tuned_on=frozenset({"b"}))})
    assert tuned["counted"] == ["a"] and tuned["excluded"] == {"b": "change was tuned on it"}
    # the captain's reading (2026-11-01): a pseudo-task the base cannot run on is not counted, and the record says so
    narrow = replace(CANDIDATES["copy"], name="narrow", requires=lambda t: t.key == "a")
    rb = _result("b", cands=("floor", "copy"))
    ra = _result("a", cands=("floor", "copy"))
    for r in (ra, rb):
        for row in r["rows"]:
            row["scores"]["wide"] = row["scores"]["copy"] + 1.0
        r["candidates"] = [*r["candidates"], "wide"]
    for row in ra["rows"]:
        row["scores"]["narrow"] = row["scores"]["copy"]
    ra["candidates"].append("narrow")
    wide = replace(CANDIDATES["copy"], name="wide")
    d = harness.decide([ra, rb], CONFIGURED, "narrow", "wide", candidates=CANDIDATES | {"narrow": narrow, "wide": wide})
    assert d["keep"] and d["counted"] == ["a"] and d["excluded"] == {"b": "base does not apply"}
    assert d["decided_on"].startswith("1 of 2")


@pytest.mark.parametrize("bad, match", [
    (dict(commit="other"), "commits"), (dict(tree_clean=False), "dirty"), (dict(config_sha256="d" * 64), "config"),
    (dict(files_sha256={"T2/E9.5.h5ad": "g"}), "two checksums"),
])
def test_decide_refuses_results_that_do_not_belong_together(bad, match):
    with pytest.raises(harness.InconsistentResults, match=match):
        harness.decide([_result("a"), _result("b", **bad)], CONFIGURED, "floor", "copy")


def test_decide_refuses_unknown_repeated_or_off_protocol_results():
    with pytest.raises(harness.InconsistentResults, match="not in the config"):
        harness.decide([_result("a"), _result("zzz")], CONFIGURED, "floor", "copy")
    with pytest.raises(harness.InconsistentResults, match="repeated"):
        harness.decide([_result("a"), _result("a")], CONFIGURED, "floor", "copy")
    off = _result("b") | {"build": "prefix"}
    with pytest.raises(harness.InconsistentResults, match="build"):
        harness.decide([_result("a"), off], CONFIGURED, "floor", "copy")
    shifted = _result("b")
    shifted["rows"][3]["build_seed"] = 3
    with pytest.raises(harness.InconsistentResults, match="build seeds"):
        harness.decide([_result("a"), shifted], CONFIGURED, "floor", "copy")
    with pytest.raises(harness.InconsistentResults, match="config"):
        harness.decide([_result("a"), _result("b")], CONFIGURED, "floor", "copy", config_sha256="e" * 64)
    with pytest.raises(harness.InconsistentResults, match="t1_held_flat: two checksums"):
        harness.decide([_result("a", derived_sha256={"t1_held_flat": "x"}),
                        _result("b", derived_sha256={"t1_held_flat": "y"})], CONFIGURED, "floor", "copy")


# ------------------------------------------------------------------------------------------------ registry
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
    shuffled = tmp_path / "shuffled.csv"
    shuffled.write_text(",".join(reversed(registry.HEADER)) + "\n")
    with pytest.raises(ValueError, match="another header"):
        registry.append([row], path=shuffled)


def test_t1_scope_rule_judges_tissue_mix_changes_on_the_chamber_cut_truth():
    """The captain's 2026-11-01 reading: same-scope T1 candidates are judged on every counted T1 pseudo-task; a
    candidate of another tissue mix (B0 against the floor, decision point D2) only on the chamber-cut truth."""
    tasks = {t.key: t for t in load_tasks()["T1:test"]}
    whole, cut, early = tasks["T1/E9.5->E10.5"], tasks["T1/E9.5->E10.5:chamber"], tasks["T1/E8.5->E9.5"]
    assert cut.counted and cut.target_filter["keep"] == "chamber"

    def res(key, cands):
        r = _result(key, cands=tuple("copy" if c == "b0" else c for c in cands))
        r["task"]["test_board"], r["candidates"] = "T1:test", list(cands)
        if "b0" in cands:
            for row in r["rows"]:
                row["scores"]["b0"] = row["scores"].pop("copy")
        return r

    b0 = replace(CANDIDATES["copy"], name="b0", params={"scope": "chamber"},
                 requires=lambda t: "start" in t.inputs)
    cands = CANDIDATES | {"b0": b0}
    results = [res(k, ("floor", "b0")) for k in (whole.key, cut.key)] + [res(early.key, ("floor",))]
    d = harness.decide(results, [whole, cut, early], "floor", "b0", candidates=cands)
    assert d["counted"] == [cut.key] and d["keep"]
    assert d["excluded"][whole.key].startswith("different tissue mixes")
    same = harness.decide([res(k, ("floor", "copy")) for k in (whole.key, cut.key, early.key)],
                          [whole, cut, early], "floor", "copy")
    assert same["counted"] == [whole.key, cut.key, early.key]


# ------------------------------------------------------------------------------------------------ B0 (feat-009)
def _counted_pool(n, rng, means, props, groups, scale=1.0):
    types = rng.choice(len(props), size=n, p=props)
    counts = rng.poisson(means[types]).astype(np.int32)
    lib = counts.sum(1, keepdims=True).clip(min=1)
    X = np.log1p(1e4 * counts / lib).astype(np.float32)
    centres = np.array([[0, 0, 0], [6, 0, 0], [0, 6, 0]], float) * scale
    coords = (centres[types] + rng.normal(size=(n, 3)) * 1.5 * scale).astype(np.float32)
    a = ad.AnnData(X=X, obs=pd.DataFrame({"celltype": [f"type{k}" for k in types],
                                          "merged": [groups[k] for k in types]}, index=[f"c{i}" for i in range(n)]),
                   var=pd.DataFrame(index=GENES))
    a.layers["counts"] = counts
    a.obsm["spatial_3D"] = coords
    return a


@pytest.fixture(scope="module")
def counted_t2():
    rng = np.random.default_rng(23)
    base = rng.gamma(2.0, 2.0, size=(3, len(GENES)))
    later = base.copy()
    later[:, :8] *= 3.0
    groups = ["CM", "Endothelial", "Untracked"]
    start = _counted_pool(3000, rng, base, [0.5, 0.3, 0.2], groups)
    ref = _counted_pool(4000, rng, later, [0.4, 0.4, 0.2], groups, scale=1.1)
    board = Board(key="pseudo:b0", task="T2", genes=list(GENES), needs_coords=True, setting="heart",
                  min_cells=100, max_cells=1500)
    return board, {"ref": ref, "start": start}


def test_b0_on_t2_moves_the_copied_cells_in_count_space(counted_t2):
    board, inputs = counted_t2
    view = harness.TaskView.of(_task(key="b0syn"))
    ctx = harness.BuildContext(task=view, board=board, inputs=inputs, seed=harness.BUILD_SEED_OFFSET + 2)
    pool, rows, X0, coords0 = harness.draw_copy(ctx)
    pred = CANDIDATES["b0"].build(ctx)
    copy = CANDIDATES["copy"].build(ctx)
    assert np.array_equal(pred.coords, copy.coords) and pred.n == copy.n == board.max_cells  # the copy's own cells
    assert not np.allclose(pred.X, copy.X)
    assert np.all(pred.X[copy.X == 0] == 0)  # count-space edits keep zeros
    untracked = np.asarray(pool.obs["merged"])[rows] == "Untracked"
    assert untracked.any() and np.allclose(pred.X[untracked], copy.X[untracked])  # held flat


def test_b0_interpolates_by_the_target_fraction_of_the_bracket(counted_t2, monkeypatch):
    board, inputs = counted_t2
    task = _task(**INTERP, inputs={"before": "T2/E7.25.h5ad", "after": "T2/E8.0.h5ad"})
    seen = {}
    real = harness.t2_build

    def spy(start, end, like, damping, held, **kw):
        seen["damping"] = damping
        from vec import t1_scope
        return t1_scope.__dict__["_momentum_shift_real"](start, end, like, damping, held, **kw)

    from vec import t1_scope
    monkeypatch.setitem(t1_scope.__dict__, "_momentum_shift_real", t1_scope.momentum_shift)
    monkeypatch.setattr(t1_scope, "momentum_shift", spy)
    ctx = harness.BuildContext(task=harness.TaskView.of(task), board=replace(board, setting="embryo"),
                               inputs={"before": inputs["start"], "after": inputs["ref"]},
                               seed=harness.BUILD_SEED_OFFSET)
    CANDIDATES["b0"].build(ctx)
    assert seen["damping"] == pytest.approx((7.5 - 7.25) / (8.0 - 7.25)) and real is harness.t2_build
    assert CANDIDATES["b0"].applies(task) and not CANDIDATES["b0"].applies(_task(**KO, inputs={"wt": KO["reference"]}))


# ------------------------------------------------------------------------------------------------ T2 embryo (feat-010)
def test_embryo_candidates_resize_to_the_size_rule_and_interpolate_inside_the_guard():
    rng = np.random.default_rng(21)
    base = rng.gamma(2.0, 2.0, size=(3, len(GENES)))
    earlier = _pool(1500, rng, base * 0.8, [0.6, 0.3, 0.1], scale=0.7)
    before = _pool(2500, rng, base, [0.5, 0.3, 0.2])
    after = _pool(4000, rng, base * 1.5, [0.3, 0.4, 0.3], scale=1.8, stretch=(1.0, 1.4, 1.0))
    board = Board(key="pseudo:synthetic_embryo", task="T2", genes=list(GENES), needs_coords=True, setting="embryo",
                  min_cells=100, max_cells=1500)
    task = _task(**INTERP, key="syn_embryo",
                 inputs={"before": "T2/E7.25.h5ad", "after": "T2/E8.0.h5ad", "earlier": "T2/E6.75.h5ad"})
    tasks = {t.key: t for ts in load_tasks().values() for t in ts}
    assert CANDIDATES["fgw"].applies(task) and CANDIDATES["psb_size"].applies(task)
    assert CANDIDATES["fgw"].applies(tasks["T2e/E6.75+E7.5->E7.25"]) and not CANDIDATES["psb_size"].applies(
        tasks["T2e/E6.75+E7.5->E7.25"])  # no stage before E6.75: the shift, and so its resized version, cannot run
    assert not any(CANDIDATES[c].applies(tasks["T2h/<=E9.5->E10.5"]) for c in ("fgw", "psb_size"))
    ctx = harness.BuildContext(task=harness.TaskView.of(task), board=board,
                               inputs={"before": before, "after": after, "earlier": earlier}, seed=7)
    want = t2_rms = harness.target_rms(ctx)
    radii = [harness.t2_build.rms_radius(a.obsm["spatial_3D"]) for a in (earlier, before, after)]
    days = [6.75, 7.25, 8.0]
    assert want == pytest.approx(harness.t2_build.pchip_rms(days, radii, 7.5)) and radii[1] < want < radii[2]
    psb, sized = CANDIDATES["pseudobulk_shift"].build(ctx), CANDIDATES["psb_size"].build(ctx)
    np.testing.assert_array_equal(psb.X, sized.X)  # the same cells and expression: only the size step differs
    pool = np.asarray(before.obsm["spatial_3D"], float)
    f = t2_rms / harness.t2_build.rms_radius(pool)
    np.testing.assert_allclose(sized.coords, pool.mean(0) + (psb.coords - pool.mean(0)) * f, rtol=1e-5, atol=1e-3)
    fgw = CANDIDATES["fgw"].build(ctx)
    assert fgw.n == board.max_cells and harness.t2_build.rms_radius(fgw.coords) == pytest.approx(want, rel=0.01)
    assert fgw.X.shape == (board.max_cells, len(GENES))
    # the embryo s2 shape variants: fgw's cells and expression, positions from one bracket's geometry, same size rule
    for name in ("fgw_shape0", "fgw_shape1"):
        assert CANDIDATES[name].applies(task) and CANDIDATES[name].applies(tasks["T2e/E6.75+E7.5->E7.25"])
        assert not CANDIDATES[name].applies(tasks["T2h/<=E9.5->E10.5"])
        v = CANDIDATES[name].build(ctx)
        np.testing.assert_array_equal(v.X, fgw.X)
        assert harness.t2_build.rms_radius(v.coords) == pytest.approx(want, rel=0.01)
        assert not np.allclose(v.coords, fgw.coords)


def test_depth_tries_keep_fgw_geometry_and_change_only_the_expression():
    """Tries E7 and E8 (experiments/harness/TRIES.md): fgw's cells and positions, expression at the predicted depth."""
    from scipy import sparse as sp

    rng = np.random.default_rng(22)
    base = rng.gamma(2.0, 2.0, size=(3, len(GENES)))
    pools = {}
    for role, (n, f, sc) in {"before": (2500, 1.0, 1.0), "after": (4000, 4.0, 1.8)}.items():
        types = rng.choice(3, size=n, p=[0.4, 0.3, 0.3])  # the later stage four times deeper
        C = rng.poisson(base[types] * f)
        C[C.sum(1) == 0, 0] = 1
        X = np.log1p(1e4 * C / C.sum(1, keepdims=True)).astype(np.float32)
        centres = np.array([[0, 0, 0], [6, 0, 0], [0, 6, 0]], float) * sc
        a = ad.AnnData(X=sp.csr_matrix(X), obs=pd.DataFrame({"celltype": [f"type{t}" for t in types]},
                                                             index=[f"c{i}" for i in range(n)]),
                       var=pd.DataFrame(index=GENES))
        a.layers["counts"] = sp.csr_matrix(C.astype(np.int32))
        a.obsm["spatial_3D"] = (centres[types] + rng.normal(size=(n, 3)) * 1.5 * sc).astype(np.float32)
        pools[role] = a
    board = Board(key="pseudo:synthetic_embryo", task="T2", genes=list(GENES), needs_coords=True, setting="embryo",
                  min_cells=100, max_cells=1500)
    task = _task(**INTERP, key="syn_embryo_depth", inputs={"before": "T2/E7.25.h5ad", "after": "T2/E8.0.h5ad"})
    ctx = harness.BuildContext(task=harness.TaskView.of(task), board=board, inputs=pools, seed=7)
    fgw = CANDIDATES["fgw"].build(ctx)
    tasks = {t.key: t for ts in load_tasks().values() for t in ts}
    for name in ("fgw_depth", "fgw_blend"):
        c = CANDIDATES[name]
        assert c.applies(task) and c.applies(tasks["T2e/E6.75+E7.5->E7.25"])
        assert not c.applies(tasks["T2h/<=E9.5->E10.5"])
        v = c.build(ctx)
        np.testing.assert_array_equal(v.coords, fgw.coords)  # positions and draws are fgw's
        assert v.X.shape == fgw.X.shape and not np.array_equal(v.X, fgw.X)
        assert ((v.X > 0).sum(1) >= 1).all()
