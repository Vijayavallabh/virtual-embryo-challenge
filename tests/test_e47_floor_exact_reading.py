"""experiments/harness/e47_floor_exact.py (tries E47a and E47b of experiments/harness/TRIES.md), part B: the guard's
arithmetic on the committed task files of tries E24, E25 and E28 (read, never recomputed), and the reading's rules on
hand-made task results. Nothing here builds a prediction or loads a competition file."""

from __future__ import annotations

import copy
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

from vec.adoption import SEEDS

REPO = Path(__file__).resolve().parents[1]
METRICS = ("de_score", "de_direction", "mmd_u", "variogram", "d2_shape", "occupancy_dice", "scale_log_ratio",
           "neighborhood_mmd")


@pytest.fixture(scope="module")
def e47():
    spec = importlib.util.spec_from_file_location("harness_e47_floor_exact",
                                                  REPO / "experiments" / "harness" / "e47_floor_exact.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


# ------------------------------------------------------------------------------------------------ the guard
def test_the_committed_files_hold_every_reference_arm_on_every_task(e47):
    committed, why = e47.committed_scores()
    assert why == []  # the committed runs agree among themselves (nine_tasks' own check)
    keys = [t.key for t in e47.tasks().values()]
    assert len(keys) == 9
    for key in keys:
        for arm in e47.REFERENCE:
            assert committed[key][arm].shape == (len(SEEDS),) and np.isfinite(committed[key][arm]).all()


def test_the_guard_on_the_committed_files(e47):
    committed, _ = e47.committed_scores()
    keys = [t.key for t in e47.tasks().values()]
    worst, why = e47.guard(committed, committed, keys)
    assert why == [] and worst == {"mean": 0.0, "seed": 0.0}  # the committed scores against themselves
    key, arm = keys[4], "fgw_shape0"

    def shifted(delta):
        run = copy.deepcopy(committed)
        run[key][arm] = run[key][arm] + delta
        return e47.guard(run, committed, keys)

    assert shifted(np.full(10, 0.09))[1] == []  # a mean 0.09 off: inside
    worst, why = shifted(np.full(10, 0.11))  # a mean 0.11 off: void
    assert len(why) == 1 and "mean" in why[0] and worst["mean"] == pytest.approx(0.11)
    one = np.zeros(10)
    one[3] = 0.49
    assert shifted(one)[1] == []  # one seed 0.49 off (mean 0.049): inside
    one[3] = 0.51
    worst, why = shifted(one)  # one seed 0.51 off (mean 0.051): void on the seed
    assert len(why) == 1 and "in a seed" in why[0] and worst["seed"] == pytest.approx(0.51)
    run = copy.deepcopy(committed)
    del run[key]["fgw"]
    assert any("no recomputed" in w for w in e47.guard(run, committed, keys)[1])


def test_the_guard_on_noise_at_the_seed_level(e47):
    """Pure noise: a recomputed arm whose seeds scatter by SD 0.3 around the committed ones keeps its mean within 0.1
    most of the time but crosses 0.5 in some seed of nine tasks and four arms: the guard then voids the run. The
    reproductions on record differ by at most 0.004 per seed (decision log 2026-11-18 (n))."""
    committed, _ = e47.committed_scores()
    keys = [t.key for t in e47.tasks().values()]
    rng = np.random.default_rng(0)
    voids = 0
    for _ in range(50):
        run = {k: {a: v + rng.normal(0.0, 0.3, v.shape) for a, v in arms.items()} for k, arms in committed.items()}
        voids += bool(e47.guard(run, committed, keys)[1])
    assert voids == 50
    run = {k: {a: v + rng.normal(0.0, 0.004, v.shape) for a, v in arms.items()} for k, arms in committed.items()}
    assert e47.guard(run, committed, keys)[1] == []


# ------------------------------------------------------------------------------------------------ the reading
def _results(e47, rng, change=None, commit="c" * 40):
    """Hand-made task results for every arm: reference and control arms random, each new arm its control plus
    `change(task_name, arm) -> (board delta per seed, {metric: skill delta per seed})` (default: equal)."""
    out = {}
    for name, t in e47.tasks().items():
        rows = []
        base = {a: rng.normal(60.0, 3.0, len(SEEDS)) for a in e47.ARMS}
        sk = {a: {m: rng.uniform(0.3, 0.9, len(SEEDS)) for m in METRICS} for a in e47.ARMS}
        for arm, ctl in e47.NEW.items():
            db, dm = change(name, arm) if change else (np.zeros(len(SEEDS)), {})
            base[arm] = base[ctl] + db
            sk[arm] = {m: sk[ctl][m] + dm.get(m, 0.0) for m in METRICS}
        for i, s in enumerate(SEEDS):
            rows.append({"seed": s, "scores": {a: float(base[a][i]) for a in e47.ARMS},
                         "skills": {a: {m: float(sk[a][m][i]) for m in METRICS} for a in e47.ARMS}})
        out[t.key] = {"task": {"key": t.key}, "rows": rows, "provenance": {"commit": commit, "tree_clean": True},
                      "floor_exact_fallback": {str(s): False for s in e47.build_seeds()}}
    return out


def _committed_of(e47, results):
    return ({k: {a: np.array([row["scores"][a] for row in sorted(r["rows"], key=lambda x: x["seed"])])
                 for a in e47.REFERENCE} for k, r in results.items()}, [])


def _read(e47, results):
    return e47.reading(results, committed=_committed_of(e47, results))


def test_truths_merge_the_tasks_that_share_a_target(e47):
    groups = e47.truths()
    assert sorted(map(sorted, groups.values())) == sorted(
        [["e75"], ["e725"], ["half", "late"], ["new1", "x1"], ["new2", "x2", "x3"]])
    v = {n: 0.0 for n in e47.tasks()}
    v["new2"] = 0.9
    assert e47.merged(v, groups) == pytest.approx(0.9 / 3 / 5)  # one task of a three-task truth, one truth of five


def test_an_arm_equal_to_its_control_replaces_it_and_pure_noise_does_too(e47):
    """Rule (ii) is a non-inferiority bar: an arm equal to its control, or differing by zero-mean seed noise, replaces
    it in the question. A shift of every task's d2_shape by -0.05 also passes; nothing in (ii) asks for a gain."""
    rng = np.random.default_rng(1)
    res = _read(e47, _results(e47, rng))
    assert not res["void"] and all(g["replaces"] for g in res["rule2"].values())
    noise = _read(e47, _results(e47, rng, lambda n, a: (rng.normal(0, 0.3, 10),
                                                        {"d2_shape": rng.normal(0, 0.01, 10),
                                                         "occupancy_dice": rng.normal(0, 0.01, 10),
                                                         "neighborhood_mmd": rng.normal(0, 0.005, 10)})))
    assert all(g["replaces"] for g in noise["rule2"].values())
    assert 0 < noise["noise"]["fgw_shape0x"]["merged_board_se"] < 0.3
    small = _read(e47, _results(e47, rng, lambda n, a: (np.zeros(10), {"d2_shape": np.full(10, -0.05)})))
    assert all(g["replaces"] for g in small["rule2"].values())


@pytest.mark.parametrize(("which", "change"), [
    ("a_board", lambda n, a: (np.full(10, -0.31), {})),  # every task 0.31 lower: the merged mean too
    ("b_nbhd", lambda n, a: (np.zeros(10), {"neighborhood_mmd": np.full(10, -0.021)})),
    ("c_shape", lambda n, a: (np.zeros(10), {"occupancy_dice": np.full(10, -0.11 if n == "x2" else 0.0)})),
    ("c_shape", lambda n, a: (np.zeros(10), {"d2_shape": np.full(10, -0.101 if n == "e75" else 0.2)})),
])
def test_each_clause_of_rule_two_keeps_the_control(e47, which, change):
    res = _read(e47, _results(e47, np.random.default_rng(2), change))
    for arm, g in res["rule2"].items():
        assert not g[which] and not g["replaces"] and g["outcome"] == f"{e47.NEW[arm]} stays"


def test_the_board_clause_reads_the_merged_mean_and_its_boundary(e47):
    # -0.9 on new2 alone: merged -0.06, unmerged -0.1; both inside the bar
    res = _read(e47, _results(e47, np.random.default_rng(3),
                              lambda n, a: (np.full(10, -0.9 if n == "new2" else 0.0), {})))
    g = res["rule2"]["fgw_shape0x"]
    assert g["board_merged"] == pytest.approx(-0.06) and g["a_board"] and g["replaces"]
    # -1.5 on both one-task truths: merged -0.6, beyond the bar
    res = _read(e47, _results(e47, np.random.default_rng(3),
                              lambda n, a: (np.full(10, -1.5 if n in ("e75", "e725") else 0.0), {})))
    assert res["rule2"]["fgw_shape0x"]["board_merged"] == pytest.approx(-0.6) and not res["rule2"]["fgw_shape0x"][
        "a_board"]


def test_the_two_arms_are_read_apart(e47):
    def change(n, a):
        return (np.full(10, -1.0 if a == "fgw_shape0x_series" else 0.0), {})

    res = _read(e47, _results(e47, np.random.default_rng(4), change))
    assert res["rule2"]["fgw_shape0x"]["replaces"] and not res["rule2"]["fgw_shape0x_series"]["replaces"]
    assert res["tasks"]["half"]["paired"]["fgw_shape0x_series"]["board"]["mean"] == pytest.approx(-1.0)
    assert res["tasks"]["half"]["paired"]["fgw_shape0x_series"]["board"]["ahead"] == 0


def test_the_pricing_table_is_against_fgw_both_ways(e47):
    rng = np.random.default_rng(5)
    results = _results(e47, rng)
    names = list(e47.tasks())
    for i, t in enumerate(e47.tasks().values()):
        for row in results[t.key]["rows"]:
            row["scores"]["fgw_aligned_series"] = row["scores"]["fgw"] + (1.0 if i % 2 == 0 else -2.0)
    p = _read(e47, results)["pricing"]["fgw_aligned_series"]
    gain = [n for i, n in enumerate(names) if i % 2 == 0]
    assert p["gain_tasks"] == len(gain) and p["loss_tasks"] == 9 - len(gain)
    assert p["mean_gain"] == pytest.approx(1.0) and p["mean_loss"] == pytest.approx(-2.0)
    assert p["tasks"][names[0]]["pair_value"] == pytest.approx(1.0)  # the better of two beside fgw
    assert p["tasks"][names[1]]["pair_value"] == pytest.approx(0.0)
    assert p["worst"][1] == pytest.approx(-2.0) and p["best"][1] == pytest.approx(1.0)
    assert set(_read(e47, results)["pricing"]) == set(e47.PRICED)


def test_void_conditions(e47):
    rng = np.random.default_rng(6)
    good = _results(e47, rng)
    assert not _read(e47, good)["void"]
    # a new arm's expression metric differs from its control's: they hold the same cells, so the run is void
    bad = copy.deepcopy(good)
    bad[next(iter(bad))]["rows"][0]["skills"]["fgw_shape0x"]["mmd_u"] += 1e-6
    res = _read(e47, bad)
    assert res["void"] and "mmd_u" in res["why"][0]
    # a reference arm off its committed scores
    res = e47.reading(good, committed=({k: {a: v + 0.2 for a, v in arms.items()}
                                        for k, arms in _committed_of(e47, good)[0].items()}, []))
    assert res["void"] and any("mean" in w for w in res["why"])
    # a missing task, a missing arm, two commits, a dirty tree
    miss = copy.deepcopy(good)
    del miss[next(iter(miss))]
    assert _read(e47, miss)["void"]
    arm = copy.deepcopy(good)
    for row in arm[next(iter(arm))]["rows"]:
        del row["scores"]["fgw_shape0x_series"]
    assert e47.reading(arm, committed=_committed_of(e47, good))["void"]
    two = copy.deepcopy(good)
    two[next(iter(two))]["provenance"]["commit"] = "d" * 40
    assert _read(e47, two)["void"]
    dirty = copy.deepcopy(good)
    dirty[next(iter(dirty))]["provenance"]["tree_clean"] = False
    assert _read(e47, dirty)["void"]
    assert e47.text(_read(e47, dirty)).startswith("VOID")
    assert "rule (ii)" in e47.text(_read(e47, good))
    # a task file without the floor-exact fallback record for every build seed
    nofb = copy.deepcopy(good)
    del nofb[next(iter(nofb))]["floor_exact_fallback"]
    assert _read(e47, nofb)["void"]


def test_a_fallback_is_listed_and_voids_nothing(e47):
    good = _results(e47, np.random.default_rng(7))
    res = _read(e47, good)
    assert res["fallback_tasks"] == [] and "fallback (" in e47.text(res) and "on no task" in e47.text(res)
    names = list(e47.tasks())
    fb = copy.deepcopy(good)
    fb[e47.tasks()["e725"].key]["floor_exact_fallback"][str(e47.build_seeds()[0])] = True
    res = _read(e47, fb)
    assert not res["void"] and res["fallback_tasks"] == ["e725"] and "e725" in names
    assert "e725" in e47.text(res).splitlines()[1]


def test_the_run_records_the_fallback_seen_by_the_builds(e47):
    """floor_exact_records sees what _floor_exact_positions reports, per build seed, and restores the function."""
    from test_e47_floor_exact import _two_types

    from vec import t2_build as tb

    orig = tb._floor_exact_positions
    Cb, Ca, Eb, Ea, _, _ = _two_types(np.random.default_rng(6), n_b=100, n_a=300)
    with e47.floor_exact_records() as seen:
        for seed, n in ((5, 80), (6, 130)):
            tb.interpolate_clouds(Cb, Ca, Eb, Ea, 0.5, n, seed=seed, target_rms=1.0, n_anchor=60, k=3, t_pos=0.0,
                                  positions=tb.POSITIONS[1])
        tb.interpolate_clouds(Cb, Ca, Eb, Ea, 0.5, 130, seed=7, target_rms=1.0, n_anchor=60, k=3, t_pos=0.0)
    assert seen == {5: False, 6: True}  # the default position mode never calls it
    assert tb._floor_exact_positions is orig
    assert e47.build_seeds() == [1_000_000 + s for s in SEEDS]


# ------------------------------------------------------------------------------------ before the row (J5, J7)
def test_the_guards_dry_run_on_committed_reproductions(e47):
    """J7: try E25's reference arms against try E24's committed scores (two commits, six tasks) sit within the guard,
    and the expression-identity guard holds on every committed pair that shares its cells. Both guards can fail on
    such data: a shift of 0.11 or 1e-6 voids, as test_the_guard_on_the_committed_files and test_void_conditions show."""
    d = e47.dry_run_guards()
    ref = d["reference"]
    assert ref["tasks"] == 6 and ref["arms"] == ["fgw", "fgw_shape0", "fgw_shape0_series"]
    assert ref["commits"] == ["798908a", "ccd33b5"] and ref["why"] == []
    assert ref["worst"]["mean"] <= e47.GUARD_MEAN and ref["worst"]["seed"] <= e47.GUARD_SEED
    expr = d["expression"]
    assert expr["files"] == 15 and expr["why"] == [] and expr["worst"] <= e47.EXPR_TOL


def test_the_committed_noise_is_spreads_only(e47):
    """J5 and the noise of rule (ii) from committed scores: spreads and chances, no mean of any arm."""
    n = e47.committed_noise()
    j5 = n["j5"]
    assert j5["noise_floor"] == max(j5["seed_se"], j5["between_truth_sd"]) > 0
    assert set(j5["task_se"]) == set(e47.tasks())
    p = n["rule2_proxy"]
    assert p["pair"] == ["fgw", "fgw_shape0"] and set(p["skill_sd"]) == set(e47.tasks())
    for v in (p["p_fail_a"], p["p_fail_b"], p["p_fail_c"]):
        assert 0.0 <= v <= 1.0


def test_the_run_refuses_without_its_row(e47, monkeypatch, tmp_path, capsys):
    monkeypatch.setattr(e47.n12, "tries_rows", lambda path: {"E28a"})
    with pytest.raises(SystemExit, match="E47a has no row"):
        e47.main(["run", "--task", "half", "--out", str(tmp_path)])
    assert not any(tmp_path.iterdir())
    assert e47.ARMS == ("floor", "fgw", "fgw_shape0", "fgw_shape0_series", "fgw_aligned_series", "fgw_shape0x",
                        "fgw_shape0x_series")
