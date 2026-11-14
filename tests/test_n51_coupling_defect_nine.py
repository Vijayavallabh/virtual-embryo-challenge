"""Row N51 (experiments/harness/n51_coupling_defect_nine.py): the coupling-defect flag on the nine practice tasks, on
a tiny synthetic task set whose stage files are written to a temporary folder. Nothing is read from data/ or boards/,
and nothing is scored."""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from vec import harness
from vec.io import COORD_KEY

HARNESS = Path(__file__).resolve().parents[1] / "experiments" / "harness"
GENES = [f"g{i}" for i in range(24)]
E75, E725 = "T2e/E7.25+E8.0->E7.5", "T2e/E6.75+E7.5->E7.25"


@pytest.fixture(scope="module")
def n51():
    sys.path.insert(0, str(HARNESS))
    spec = importlib.util.spec_from_file_location("n51_coupling_defect_nine", HARNESS / "n51_coupling_defect_nine.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _task(key, before, after, target):
    return harness.check_task(harness.PseudoTask(
        key=key, test_board="T2:embryo:test_interp", task="T2", kind="interp", target=f"T2/{target}.h5ad",
        reference=f"T2/{before}.h5ad", inputs={"before": f"T2/{before}.h5ad", "after": f"T2/{after}.h5ad"},
        twin_board="T2:embryo:val_interp", setting="embryo", adoption=False))


TASKS = {"e75": _task(E75, "E7.25", "E8.0", "E7.5"), "e725": _task(E725, "E6.75", "E7.5", "E7.25"),
         "x1": _task("T2e/E8.0+E8.75->E8.25_late", "E8.0", "E8.75", "E8.25_late")}


def _stage(seed: int, n: int = 90) -> ad.AnnData:
    rng = np.random.default_rng(seed)
    X = rng.poisson(2.0, size=(n, len(GENES))).astype(np.float32)
    obs = pd.DataFrame({"merged": rng.choice(["a", "b", "c"], n)}, index=[f"c{i}" for i in range(n)])
    a = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=GENES))
    a.obsm[COORD_KEY] = rng.normal(size=(n, 3)) * np.array([3.0, 2.0, 1.0])
    return a


def fake_measure(inputs, frac, seed):
    """A measure that depends on the draw only: seeds 0-2 flagged by xcorr, 3-4 by thin, 5 by both, 6-9 not."""
    s = seed - harness.BUILD_SEED_OFFSET
    x = 0.1 if s in (0, 1, 2, 5) else 0.8
    th = 0.5 if s in (3, 4, 5) else 0.95
    return {"xcorr_min": x + frac * 1e-3, "thin_ratio": th}


def committed_for(tasks, shift=0.0):
    """A committed diagnostic matching fake_measure on the two backtests, one value moved by `shift`."""
    pairs = {}
    for name in ("e75", "e725"):
        t = tasks[name]
        frac = harness.TaskView.of(t)
        frac = (frac.target_day - frac.role_days["before"]) / (frac.role_days["after"] - frac.role_days["before"])
        draws = [{"seed": s, "fgw": fake_measure(None, frac, harness.BUILD_SEED_OFFSET + s)} for s in range(10)]
        pairs[t.key] = {"before": harness.stem(t.inputs["before"]), "after": harness.stem(t.inputs["after"]),
                        "draws": draws}
    pairs[E725]["draws"][4]["fgw"]["thin_ratio"] += shift
    return {"commit": "c0", "pairs": pairs}


@pytest.fixture()
def world(n51, tmp_path, monkeypatch):
    data = tmp_path / "data"
    (data / "T2").mkdir(parents=True)
    for i, s in enumerate(["E6.75", "E7.25", "E7.5", "E8.0", "E8.75", "E8.25_late"]):
        _stage(i).write_h5ad(data / "T2" / f"{s}.h5ad")
    opened = []

    def read(path, genes):
        opened.append(Path(path))
        return ad.read_h5ad(path)

    tries = tmp_path / "TRIES.md"
    tries.write_text("| Row | what |\n|---|---|\n| N51 (diagnostic) | the flag on nine tasks |\n")
    committed = tmp_path / "diag.json"
    committed.write_text(json.dumps(committed_for(TASKS)))
    monkeypatch.setattr(n51, "nine", lambda: dict(TASKS))
    monkeypatch.setattr(n51, "_read", read)
    monkeypatch.setattr(n51, "_genes", lambda t: list(GENES))
    monkeypatch.setattr(n51, "_state", lambda: ("a" * 40, True, "2026-11-14T00:00:00Z"))
    monkeypatch.setattr(n51, "TRIES", tries)
    monkeypatch.setattr(n51, "COMMITTED", committed)
    monkeypatch.setattr(n51, "REPO", tmp_path)
    monkeypatch.setattr(n51, "measure_draw", fake_measure)
    return {"data": data, "out": tmp_path / "out", "opened": opened, "committed": committed, "tries": tries}


def run(n51, world, *extra):
    n51.main(["run", "--out", str(world["out"]), "--data", str(world["data"]), *extra])


def test_classify_follows_the_committed_thresholds(n51):
    assert n51.dc.DEFECT_XCORR == 0.3 and n51.dc.DEFECT_THIN == 0.9
    assert n51.classify(0.29, 0.95) == (True, "xcorr")
    assert n51.classify(0.3, 0.95) == (False, None)
    assert n51.classify(0.5, 0.89) == (True, "thin")
    assert n51.classify(0.5, 0.9) == (False, None)
    assert n51.classify(-0.4, 0.6) == (True, "both")


def test_full_run_flags_and_merges(n51, world):
    run(n51, world)
    rep = json.loads((world["out"] / "n51.json").read_text())
    assert set(rep["tasks"]) == set(TASKS)
    by = [d["by"] for d in rep["tasks"]["x1"]["draws"]]
    assert by == ["xcorr"] * 3 + ["thin"] * 2 + ["both"] + [None] * 4
    assert [d["defect"] for d in rep["tasks"]["e75"]["draws"]] == [True] * 6 + [False] * 4
    assert rep["tasks"]["e75"]["n_defect"] == 6
    assert rep["tasks"]["e725"]["fraction"] == pytest.approx(2 / 3)
    assert rep["gate"]["max_abs_difference"] == 0.0
    assert rep["clean"] is True and rep["commit"] == "a" * 40
    txt = (world["out"] / "n51.txt").read_text()
    assert "x x x t t b . . . ." in txt


def test_no_target_path_is_opened(n51, world):
    """Each task alone: every path its run asks the loader for is one of its two brackets, never its target (a
    stage can be one task's target and another's bracket, so the check is per task)."""
    for name, t in TASKS.items():
        world["opened"].clear()
        run(n51, world, "--tasks", name)
        brackets = {world["data"] / t.inputs["before"], world["data"] / t.inputs["after"]}
        assert set(world["opened"]) == brackets
        assert world["data"] / t.target not in world["opened"]


def test_a_bracket_equal_to_the_target_is_refused(n51, tmp_path):
    t = TASKS["x1"]
    bad = harness.PseudoTask(**{**t.__dict__, "target": t.inputs["after"]})  # unchecked: built past check_task
    with pytest.raises(harness.LeakageError):
        n51.bracket_paths(bad, tmp_path)


def test_gate_stops_the_run_at_1e5(n51, world):
    world["committed"].write_text(json.dumps(committed_for(TASKS, shift=1e-5)))
    with pytest.raises(SystemExit, match="reproduction gate"):
        run(n51, world)
    assert not world["out"].exists() or not any(world["out"].iterdir())  # no part file, not even e75's


def test_a_failed_gate_leaves_no_part_file_whatever_the_task_order(n51, world, monkeypatch):
    """The nine come in all_tasks() order (e75 and e725 after four others): the gate is checked before any part file
    is written, so a rerun after a failed gate is not refused as an overwrite."""
    monkeypatch.setattr(n51, "nine", lambda: {"x1": TASKS["x1"], "e75": TASKS["e75"], "e725": TASKS["e725"]})
    world["committed"].write_text(json.dumps(committed_for(TASKS, shift=1e-5)))
    with pytest.raises(SystemExit, match="reproduction gate"):
        run(n51, world)
    assert not world["out"].exists() or not any(world["out"].iterdir())
    world["committed"].write_text(json.dumps(committed_for(TASKS)))
    run(n51, world)
    rep = json.loads((world["out"] / "n51.json").read_text())
    assert list(rep["tasks"]) == ["x1", "e75", "e725"]  # the merged report keeps the task list's order


def test_gate_passes_below_tolerance(n51, world):
    world["committed"].write_text(json.dumps(committed_for(TASKS, shift=5e-7)))
    run(n51, world)
    rep = json.loads((world["out"] / "n51.json").read_text())
    assert rep["gate"]["max_abs_difference"] == pytest.approx(5e-7, rel=1e-3)


def test_the_row_must_exist(n51, world):
    world["tries"].write_text("| N50 (diagnostic) | another row |\n")
    with pytest.raises(SystemExit, match="N51 has no row"):
        run(n51, world)
    assert not world["out"].exists()


def test_refuses_a_dirty_tree(n51, world, monkeypatch):
    monkeypatch.setattr(n51, "_state", lambda: ("a" * 40, False, "t"))
    with pytest.raises(SystemExit, match="not clean"):
        run(n51, world)
    assert not world["out"].exists()


def test_refuses_an_existing_output(n51, world):
    world["out"].mkdir()
    (world["out"] / "n51_x1.json").write_text("{}")
    with pytest.raises(SystemExit, match="refusing to overwrite"):
        run(n51, world, "--tasks", "x1")
    run(n51, world, "--tasks", "e75")  # another task's part file is not in the way
    with pytest.raises(SystemExit, match="refusing to overwrite"):
        run(n51, world, "--tasks", "e75")


def test_split_runs_then_merge(n51, world):
    run(n51, world, "--tasks", "e75", "x1")
    assert not (world["out"] / "n51.json").exists()
    with pytest.raises(SystemExit, match="no part file"):
        n51.main(["merge", "--out", str(world["out"])])
    run(n51, world, "--tasks", "e725")
    n51.main(["merge", "--out", str(world["out"])])
    rep = json.loads((world["out"] / "n51.json").read_text())
    assert set(rep["tasks"]) == set(TASKS) and rep["gate"]["max_abs_difference"] == 0.0
    with pytest.raises(SystemExit, match="refusing to overwrite"):
        n51.main(["merge", "--out", str(world["out"])])


def test_merge_applies_the_gate(n51, world):
    run(n51, world, "--tasks", "e75", "e725", "x1")
    world["committed"].write_text(json.dumps(committed_for(TASKS, shift=1e-5)))
    with pytest.raises(SystemExit, match="reproduction gate"):
        n51.main(["merge", "--out", str(world["out"])])
    assert not (world["out"] / "n51.json").exists()


def test_merge_refuses_parts_of_another_commit(n51, world, monkeypatch):
    run(n51, world, "--tasks", "e75", "e725", "x1")
    monkeypatch.setattr(n51, "_state", lambda: ("b" * 40, True, "t"))
    with pytest.raises(SystemExit, match="commits"):
        n51.main(["merge", "--out", str(world["out"])])


def test_outputs_hold_no_cell_count(n51, world):
    run(n51, world)
    rep = json.loads((world["out"] / "n51.json").read_text())
    keys = {k for r in rep["tasks"].values() for d in r["draws"] for k in d} | {k for r in rep["tasks"].values()
                                                                             for k in r}
    assert keys == {"seed", "build_seed", "xcorr_min", "thin_ratio", "defect", "by",
                    "key", "brackets", "panel_board", "fraction", "n_defect", "draws"}


def test_measure_task_equals_a_harness_style_call(n51, tmp_path, monkeypatch):
    """measure_task (its real measure_draw) against a direct call of the harness's interpolate_clouds and the
    diagnostic's _measure, on clouds larger than the coupling's anchor count, so that the anchors depend on the seed.
    The direct call uses another draw size and target_rms: neither changes the coupling or the measures. The task's
    fraction is 2/3, so passing 1 - fraction anywhere would show."""
    from vec import t2_build

    n_cells = t2_build.INTERP_N_ANCHOR + 40
    data = tmp_path / "data"
    (data / "T2").mkdir(parents=True)
    for i, s in enumerate(["E6.75", "E7.5"]):
        _stage(20 + i, n_cells).write_h5ad(data / "T2" / f"{s}.h5ad")
    monkeypatch.setattr(n51, "_read", lambda path, genes: ad.read_h5ad(path))
    monkeypatch.setattr(n51, "_genes", lambda t: list(GENES))
    monkeypatch.setattr(n51, "SEEDS", (0, 3))
    t = TASKS["e725"]
    rec = n51.measure_task(t, data, log=lambda *a, **k: None)
    frac = rec["fraction"]
    assert frac == pytest.approx(2 / 3)
    before, after = (ad.read_h5ad(data / t.inputs[r]) for r in ("before", "after"))
    Cb, Ca, Eb, Ea, lab_b, lab_a = n51.dc._inputs(before, after, GENES)
    direct = []
    for s in (0, 3):
        seed = harness.BUILD_SEED_OFFSET + s
        res = t2_build.interpolate_clouds(Cb, Ca, Eb, Ea, frac, 1234, seed, target_rms=3.7, with_coupling=True)
        m = n51.dc._measure(res, frac, lab_b, lab_a, n51.dc._nn_rows(seed, len(Cb), len(Ca)))
        direct.append({k: float(m[k]) for k in n51.MEASURES})
    got = [{k: d[k] for k in n51.MEASURES} for d in rec["draws"]]
    assert got == direct
    assert direct[0] != direct[1]  # the seed reaches the anchors: the comparison is not blind to the seed


def test_read_drops_duplicate_records(n51, tmp_path):
    """_read loads as run_task does, duplicate records dropped (a copied cell a level down, as vec.prep's tests)."""
    from vec.prep import log1p_cp10k

    rng = np.random.default_rng(3)
    n, g = 300, len(GENES)
    counts = rng.poisson(4.0, size=(n, g)).astype(np.float32)
    coords = np.c_[rng.uniform(-100, 100, size=(n, 2)), np.repeat([30.0, 15.0, 0.0], n // 3)]
    counts = np.vstack([counts, counts[[5]]])
    coords = np.vstack([coords, coords[5] - [0, 0, 15.0]])
    a = ad.AnnData(X=log1p_cp10k(counts), var=pd.DataFrame(index=GENES))
    a.layers["counts"] = sparse.csr_matrix(counts)
    a.obsm[COORD_KEY] = coords
    a.obs["celltype"] = "V-CM"
    a.write_h5ad(tmp_path / "E7.5.h5ad")
    out = n51._read(tmp_path / "E7.5.h5ad", GENES[:5])
    assert out.n_obs == n and list(out.var_names) == GENES[:5]


def test_the_real_measure_runs_on_synthetic_clouds(n51):
    """The diagnostic's own measure (not faked) on two synthetic clouds: finite values, reproducible per seed."""
    inputs = n51.dc._inputs(_stage(10, 120), _stage(11, 110), GENES)
    x = n51.measure_draw(inputs, 0.5, harness.BUILD_SEED_OFFSET + 1)
    y = n51.measure_draw(inputs, 0.5, harness.BUILD_SEED_OFFSET + 1)
    assert set(x) == {"xcorr_min", "thin_ratio"} and all(np.isfinite(v) for v in x.values())
    assert x == y
