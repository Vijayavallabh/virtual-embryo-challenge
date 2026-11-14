"""experiments/t1_m1/lib_step_read.py (try T46 of experiments/harness/TRIES.md) on synthetic data: the library truths,
the scaled arms and the reading. No released file, statistics file, board file or score is read."""

from __future__ import annotations

import importlib.util
import json
import sys

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from test_t1_step_read import _arms_world

from vec import harness
from vec.boards import REPO_ROOT, Board
from vec.evaluate import Stage

SCRIPT = REPO_ROOT / "experiments" / "t1_m1" / "lib_step_read.py"
SCALED = ("m1_r_gmom_x15", "m1_r_gmom_x2", "m1_r_gmom_x3", "m1_r_gmom_x4")
MULT = {"m1_r_gmom_x15": 1.5, "m1_r_gmom_x2": 2.0, "m1_r_gmom_x3": 3.0, "m1_r_gmom_x4": 4.0}


@pytest.fixture
def ls(monkeypatch, tmp_path):
    spec = importlib.util.spec_from_file_location("lib_step_read_module", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    monkeypatch.setattr(mod, "RESULTS", tmp_path / "results")
    return mod


# ---------------------------------------------------------------------------------------------------------------------
# The truths


GENES = [f"g{i}" for i in range(8)]
LABELS = ["V-CM", "Endocardium", "Blood", "Foregut", "Neural Tube", "aPHM", "Epicardium"]


def _target(n=240, seed=3):
    """A small E10.5-like target: chamber, surrounding and pharyngeal labels, cell names ending in _1, _2 or _3 in an
    interleaved order."""
    rng = np.random.default_rng(seed)
    X = rng.poisson(2.0, (n, len(GENES))).astype(np.float32)
    lib = rng.choice(["1", "2", "3"], n)
    obs = pd.DataFrame({"celltype": rng.choice(LABELS, n)}, index=[f"AC{i:04d}-x_{lib[i]}" for i in range(n)])
    return ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=GENES))


@pytest.fixture
def chamber(monkeypatch):
    a = _target()
    monkeypatch.setattr(harness, "read_stage", lambda path, genes: a[:, list(genes)].copy())
    t = {x.key: x for x in harness.load_tasks()["T1:test"]}["T1/E9.5->E10.5:chamber"]
    return t, a


def test_the_filtered_target_and_its_names_are_filter_targets_rows(ls, chamber, tmp_path):
    t, a = chamber
    direct, rec = harness.filter_target(t, tmp_path, GENES)
    assert rec["contamination"] > 0  # the filter draws surrounding cells here, so a library keeps its share of them
    got, names = ls.filtered_target(t, tmp_path, GENES)
    assert type(got) is Stage and harness.Stage is Stage
    assert np.array_equal(got.X, direct.X) and np.array_equal(got.celltype, direct.celltype)
    rows = [a.obs_names.get_loc(n) for n in names]
    assert rows == sorted(rows)
    assert np.array_equal(np.asarray(a.X)[rows], got.X)


def test_a_library_keeps_exactly_its_rows_in_order_and_the_three_partition_the_target(ls, chamber, tmp_path):
    t, a = chamber
    target, names = ls.filtered_target(t, tmp_path, GENES)
    parts = {}
    for truth, tag in (("lib1", "1"), ("lib2", "2"), ("lib3", "3")):
        want = [i for i, n in enumerate(names) if n.endswith(f"_{tag}")]
        assert list(ls.library_rows(names, truth)) == want
        lib = ls.load_truth(t, tmp_path, GENES, truth)
        assert np.array_equal(lib.X, target.X[want]) and np.array_equal(lib.celltype, target.celltype[want])
        parts[truth] = want
    joined = sorted(i for w in parts.values() for i in w)
    assert joined == list(range(target.n))
    assert sum(len(w) for w in parts.values()) == target.n


def test_pooled_hands_the_scorer_the_very_target_t44s_run_does(ls, chamber, tmp_path):
    t, _ = chamber
    t44 = harness.filter_target(t, tmp_path, GENES)[0]
    pooled = ls.load_truth(t, tmp_path, GENES, "pooled")
    assert pooled.X.tobytes() == t44.X.tobytes() and pooled.X.dtype == t44.X.dtype
    assert np.array_equal(pooled.celltype, t44.celltype) and pooled.coords is None and t44.coords is None


def test_the_board_alone_carries_the_truth_into_a_build():
    """The candidate's board is the only route from the truth to a prediction, and its count is max_cells: a truth too
    small to move max_cells leaves the count as the pooled truth's (min_cells may still change: a check, not a count),
    a larger one changes it."""
    twin = Board(key="T1:val", task="T1", genes=GENES, needs_coords=False, split="val", min_cells=1000, max_cells=9000)
    t = harness.PseudoTask(key="T1/x", test_board="T1:test", task="T1", kind="extrap", target="T1/E10.5_RNA.h5ad",
                           reference="T1/E9.5_RNA.h5ad", inputs={"ref": "T1/E9.5_RNA.h5ad"}, twin_board="T1:val")
    cb = {n: harness._candidate_board(harness.task_board(t, {"T1:val": twin}, n, 20000)) for n in (30000, 9000, 80000)}
    count = {n: harness.count_for(b, 2000) for n, b in cb.items()}
    assert cb[30000].max_cells == cb[9000].max_cells and count[30000] == count[9000]
    assert cb[30000].min_cells != cb[9000].min_cells
    assert cb[80000].max_cells != cb[30000].max_cells and count[80000] != count[30000]
    assert all(b.truth_cells is None and b.ref_cells is None for b in cb.values())


# ---------------------------------------------------------------------------------------------------------------------
# The arms


def test_the_arms_fractions_are_point_three_to_point_eight(ls):
    assert ls.MULTIPLIERS == MULT
    fr = ls.arm_fractions(0.2)
    assert fr["m1_r_gmom"] == 0.2
    assert [fr[a] for a in SCALED] == pytest.approx([0.3, 0.4, 0.6, 0.8])


def test_an_arm_differs_from_the_base_build_in_the_fraction_only(ls, monkeypatch):
    calls = []

    def spy(ctx, directions, groups, fraction, held_m1, held_factors):
        calls.append({"ctx": ctx, "directions": directions, "groups": groups, "fraction": fraction,
                      "held_m1": held_m1, "held_factors": held_factors})
        return {"m1_r": ("m1_r", fraction), "m1_r_gmom": ("m1_r_gmom", fraction)}

    monkeypatch.setattr(ls.sr.gc, "m1_r_arms", spy)
    ctx, dirs, groups, held_m1, held_f = object(), {"CM": np.zeros(3)}, {"V-CM": "CM"}, ["a", "b"], ["a"]
    arms = ls.scaled_arms(ctx, dirs, groups, 0.2, held_m1, held_f)
    assert set(arms) == {"m1_r_gmom", *SCALED}
    assert len(calls) == 5
    base = calls[0]
    assert base["fraction"] == 0.2
    assert [c["fraction"] for c in calls[1:]] == pytest.approx([0.3, 0.4, 0.6, 0.8], abs=1e-15)
    for c in calls:
        for key in ("ctx", "directions", "groups", "held_m1", "held_factors"):
            assert c[key] is base[key], f"an arm differs from the base build in {key}"
    assert arms["m1_r_gmom_x4"] == ("m1_r_gmom", 0.2 * 4)


def test_a_multiplier_of_one_is_the_base_bit_for_bit_and_the_larger_steps_build(ls, monkeypatch):
    """The idle case: at a multiplier of 1 the arm is the base build bit for bit; fractions 0.6 and 0.8 build without
    a refusal and move the cells."""
    ctx, groups, dirs, held = _arms_world()
    monkeypatch.setattr(ls, "MULTIPLIERS", {"idle": 1.0, **MULT})
    arms = ls.scaled_arms(ctx, dirs, groups, 0.2, held, held)
    direct = ls.gc.m1_r_arms(ctx, dirs, groups, 0.2, held, held)["m1_r_gmom"]
    assert arms["m1_r_gmom"].X.tobytes() == direct.X.tobytes()
    assert arms["idle"].X.tobytes() == arms["m1_r_gmom"].X.tobytes()
    assert ls.prediction_sha256(arms["idle"]) == ls.prediction_sha256(arms["m1_r_gmom"])
    for a in SCALED:
        assert arms[a].X.shape == arms["m1_r_gmom"].X.shape and np.isfinite(arms[a].X).all()
        assert not np.array_equal(arms[a].X, arms["m1_r_gmom"].X)
        assert ls.prediction_sha256(arms[a]) != ls.prediction_sha256(arms["m1_r_gmom"])


# ---------------------------------------------------------------------------------------------------------------------
# The reading


B0, BASE = 50.0, 53.0
W = {"de_score": 0.25, "de_direction": 0.25, "mmd_u": 0.30, "variogram": 0.20}


def _arm(score, de_part=0.5, mmd_u=0.1):
    """An arm whose skills give exactly its score: `de_part` of the points above BASE on the two de terms, the rest on
    mmd_u and the variogram (split by weight)."""
    gain = score - BASE
    sk = {"de_score": 0.4 + de_part * gain / 2 / (100 * W["de_score"]),
          "de_direction": 0.4 + de_part * gain / 2 / (100 * W["de_direction"]),
          "mmd_u": 0.3 + (1 - de_part) * gain * 0.6 / (100 * W["mmd_u"]),
          "variogram": 0.2 + (1 - de_part) * gain * 0.4 / (100 * W["variogram"])}
    return {"score": score, "groups": {"distribution": score / 100}, "metrics": {"mmd_u": mmd_u, "variogram": 0.01},
            "skills": sk}


def _task(truth, gains=None, per_seed=None, seeds=range(10), commit="c" * 40, clean=True, fraction=0.2, fractions=None,
          base_shift=0.0, identical=(), rows=None, de_part=0.5, sha=None):
    """A task file. gains: arm -> gain over m1_r_gmom on every seed; per_seed: arm -> ten gains (overrides gains); a
    scaled arm's raw metric differs from the base's by 1e-9 unless named in `identical`."""
    gains = gains or {}
    out = []
    for i, s in enumerate(seeds):
        arms = {"b0": _arm(B0), "m1_r_gmom": _arm(BASE + base_shift)}
        for a in SCALED:
            g = per_seed[a][i] if per_seed and a in per_seed else gains.get(a, 0.0)
            arms[a] = (dict(arms["m1_r_gmom"]) if a in identical else
                       _arm(BASE + base_shift + g, de_part=de_part, mmd_u=0.1 + 1e-9))
        out.append({"seed": s, "floor": 50.0, "arms": arms,
                    "sha256": {a: (sha or {}).get(a, f"{a}-{s}") for a in ("b0", "m1_r_gmom", *SCALED)}})
    fr = fractions or {"m1_r_gmom": fraction, **{a: fraction * m for a, m in MULT.items()}}
    return {"label": "synthetic", "task": "T1/E9.5->E10.5:chamber", "truth": truth, "truth_rows": rows,
            "fraction": fraction, "fractions": fr, "multipliers": dict(MULT), "seeds": list(seeds), "rows": out,
            "stamp": {"commit": commit, "tree_clean": clean, "utc": "x", "inputs": {}}}


ROWS = {"pooled": 7919, "lib1": 1777, "lib2": 2711, "lib3": 3431}


def _write(ls, tmp_path, truths, t44_shift=0.0):
    """truths: truth -> task record. T42's and T44's committed chamber files, synthetic, beside them."""
    out = tmp_path / "out"
    (out / "tasks").mkdir(parents=True)
    for k, rec in truths.items():
        (out / "tasks" / f"{k}.json").write_text(json.dumps(rec))
    t42, t44 = tmp_path / "t42.json", tmp_path / "t44.json"
    t42.write_text(json.dumps({"rows": [{"seed": s, "arms": {"b0": _arm(B0), "m1_r_gmom": _arm(BASE)}}
                                        for s in range(10)]}))
    pooled = truths.get("pooled")
    t44_rows = []
    for s in range(10):
        arms = {"b0": _arm(B0), "m1_r_gmom": _arm(BASE)}
        if pooled is not None:
            row = next((r for r in pooled["rows"] if r["seed"] == s), None)
            if row is not None:
                for a in ("m1_r_gmom_x15", "m1_r_gmom_x2"):
                    arms[a] = {**row["arms"][a], "score": row["arms"][a]["score"] + t44_shift}
        t44_rows.append({"seed": s, "arms": arms})
    t44.write_text(json.dumps({"rows": t44_rows}))
    return out, t42, t44


def _read(ls, tmp_path, by_truth=None, over=None, t44_shift=0.0, **kw):
    """by_truth: truth -> gains; over: truth -> task-file keyword overrides (None drops the truth)."""
    by_truth, over = by_truth or {}, over or {}
    truths = {}
    for k in ls.TRUTHS:
        if k in over and over[k] is None:
            continue
        args = {"gains": by_truth.get(k, {}), "rows": ROWS[k], **kw, **over.get(k, {})}
        truths[k] = _task(k, **args)
    out, t42, t44 = _write(ls, tmp_path, truths, t44_shift)
    return ls.read(out, t42, t44)["reading"]


UP = {"m1_r_gmom_x2": 1.0}


@pytest.mark.parametrize("by_truth, statement", [
    ({"lib3": UP}, "the larger step's gain is the kind's"),
    ({"lib1": UP, "lib2": UP}, "the larger step's gain is not the kind's"),
    ({"lib1": UP, "lib2": UP, "lib3": UP}, "the larger step's gain is not the kind's"),
    ({"lib1": UP, "lib3": UP}, "mixed"),
    ({"lib2": UP}, "mixed"),
    ({}, "mixed"),
    ({"pooled": UP, "lib3": {"m1_r_gmom_x15": 1.0}}, "mixed"),
])
def test_each_kind_statement_is_reached_by_a_table_built_for_it(ls, tmp_path, by_truth, statement):
    r = _read(ls, tmp_path, by_truth)
    assert r["valid"] and r["outcome"] == statement
    for k in ("lib1", "lib2", "lib3"):
        want = by_truth.get(k, {}).get("m1_r_gmom_x2", 0.0)
        g = r["kind"][k]
        assert g["gain"] == pytest.approx(want)
        assert g["de"] == pytest.approx(want / 2) and g["mmd_u_variogram"] == pytest.approx(want / 2)
        assert g["test"] == ("ahead" if want > 0.3 else "neither")
    txt = (tmp_path / "out" / "t46.txt").read_text()
    assert f"On the kind (twice the step): {statement}" in txt
    for k in ("lib1", "lib2", "lib3"):
        assert f"m1_r_gmom_x2 on {k}: " in txt
    assert not any(str(n) in txt for n in ROWS.values())  # no cell count is printed


def test_the_parts_split_the_gain_by_the_boards_weights(ls, tmp_path):
    r = _read(ls, tmp_path, {"lib1": {"m1_r_gmom_x2": 2.0}}, de_part=0.75)
    c = r["comparisons"]["lib1"]["m1_r_gmom_x2"]
    assert c["points"]["de"] == pytest.approx(1.5) and c["points"]["mmd_u_variogram"] == pytest.approx(0.5)
    assert c["points"]["de"] + c["points"]["mmd_u_variogram"] == pytest.approx(c["mean"])
    assert set(c["skills"]) == set(W)


def test_a_gain_of_exactly_point_three_is_not_ahead_and_a_loss_is_behind(ls, tmp_path):
    r = _read(ls, tmp_path, {"lib3": {"m1_r_gmom_x2": 0.3}})
    assert r["valid"] and r["comparisons"]["lib3"]["m1_r_gmom_x2"]["test"] == "neither" and r["outcome"] == ls.MIXED
    r = _read(ls, tmp_path / "exact", {"lib3": {"m1_r_gmom_x2": 0.3}}, over={"lib3": {"base_shift": -BASE}})
    c = r["comparisons"]["lib3"]["m1_r_gmom_x2"]
    assert r["valid"] and c["mean"] == pytest.approx(0.3, abs=1e-15) and c["seeds_ahead"] == 10
    assert c["test"] == "neither"
    r = _read(ls, tmp_path / "b", {"lib3": {"m1_r_gmom_x2": 0.31}})
    assert r["comparisons"]["lib3"]["m1_r_gmom_x2"]["test"] == "ahead" and r["outcome"] == ls.KIND
    r = _read(ls, tmp_path / "c", {"lib1": {"m1_r_gmom_x2": -0.3, "m1_r_gmom_x3": -0.31}})
    assert r["comparisons"]["lib1"]["m1_r_gmom_x2"]["test"] == "neither"
    assert r["comparisons"]["lib1"]["m1_r_gmom_x3"]["test"] == "behind"


def test_ahead_needs_eight_seeds(ls, tmp_path):
    seven = [1.0] * 7 + [-0.1] * 3
    eight = [1.0] * 8 + [-0.1] * 2
    r = _read(ls, tmp_path, over={"lib3": {"per_seed": {"m1_r_gmom_x2": seven}}})
    c = r["comparisons"]["lib3"]["m1_r_gmom_x2"]
    assert c["seeds_ahead"] == 7 and c["mean"] > 0.3 and c["test"] == "neither"
    assert c["range"] == pytest.approx([-0.1, 1.0])
    r = _read(ls, tmp_path / "b", over={"lib3": {"per_seed": {"m1_r_gmom_x2": eight}}})
    assert r["comparisons"]["lib3"]["m1_r_gmom_x2"]["test"] == "ahead" and r["outcome"] == ls.KIND


@pytest.mark.parametrize("gains, best, statement", [
    ({"m1_r_gmom_x15": 0.5, "m1_r_gmom_x2": 1.0, "m1_r_gmom_x3": 1.5, "m1_r_gmom_x4": 2.0}, 4.0, "still rising at 4"),
    ({"m1_r_gmom_x15": 0.5, "m1_r_gmom_x2": 1.0, "m1_r_gmom_x3": 0.4, "m1_r_gmom_x4": -1.0}, 2.0, "turns at 2"),
    ({"m1_r_gmom_x15": -0.5, "m1_r_gmom_x2": -1.0, "m1_r_gmom_x3": -2.0, "m1_r_gmom_x4": -3.0}, 1.0, "turns at 1"),
    ({"m1_r_gmom_x15": 0.5, "m1_r_gmom_x2": 0.5, "m1_r_gmom_x3": 0.2, "m1_r_gmom_x4": 0.1}, 1.5, "turns at 1.5"),
])
def test_where_the_response_turns(ls, tmp_path, gains, best, statement):
    r = _read(ls, tmp_path, {k: gains for k in ls.TRUTHS})
    for k in ls.TRUTHS:
        assert r["turn"][k]["best"] == best and r["turn"][k]["statement"] == statement
        assert r["mean_scores"][k]["m1_r_gmom"] == pytest.approx(BASE)


def test_the_predictions_sha256_agreement_is_reported(ls, tmp_path):
    same = {a: "same" for a in ("b0", "m1_r_gmom", *SCALED)}
    r = _read(ls, tmp_path, sha=same)
    assert r["valid"] and r["sha256_across_truths"]["all_agree"]
    r = _read(ls, tmp_path / "b", sha=same, over={"lib2": {"sha": {**same, "m1_r_gmom_x4": "other"}}})
    sa = r["sha256_across_truths"]
    assert r["valid"] and not sa["all_agree"] and sa["arms"]["m1_r_gmom_x4"]["seeds_differing"] == list(range(10))
    assert sa["arms"]["b0"]["agree"]


@pytest.mark.parametrize("over, kw, reason", [
    ({"lib2": None}, {}, "missing truths"),
    ({"lib1": {"seeds": range(9)}}, {}, "seeds"),
    ({"pooled": {"base_shift": 1e-6}}, {}, "do not reproduce"),
    ({}, {"t44_shift": 1e-6}, "do not reproduce"),
    ({"lib3": {"fraction": 0.4}}, {}, "fraction 0.4"),
    ({"lib3": {"fractions": {"m1_r_gmom": 0.2, "m1_r_gmom_x15": 0.3, "m1_r_gmom_x2": 0.4, "m1_r_gmom_x3": 0.6,
                             "m1_r_gmom_x4": 0.6}}}, {}, "wrong fraction"),
    ({"lib1": {"commit": "d" * 40}}, {}, "2 commits"),
    ({"lib2": {"clean": False}}, {}, "not clean"),
    ({"lib2": {"identical": ("m1_r_gmom_x3",)}}, {}, "exactly as m1_r_gmom"),
    ({"lib3": {"rows": 3430}}, {}, "do not add up"),
    ({"lib3": {"rows": None}}, {}, "do not add up"),
])
def test_the_void_cases(ls, tmp_path, over, kw, reason):
    t44_shift = kw.get("t44_shift", 0.0)
    r = _read(ls, tmp_path, {"lib3": UP}, over=over, t44_shift=t44_shift)
    assert not r["valid"] and r["outcome"] == "VOID"
    assert any(reason in x for x in r["validity"]["reasons"]), r["validity"]["reasons"]
    txt = (tmp_path / "out" / "t46.txt").read_text()
    assert not any(str(n) in txt for n in (*ROWS.values(), 3430))


def test_a_rescored_base_within_tolerance_stays_valid_and_an_arm_identical_on_some_seeds_too(ls, tmp_path):
    r = _read(ls, tmp_path, {"lib3": UP}, over={"pooled": {"base_shift": 1e-12}})
    assert r["valid"]
    per = {"m1_r_gmom_x3": [0.0] * 9 + [0.5]}
    r = _read(ls, tmp_path / "b", {"lib3": UP}, over={"lib1": {"per_seed": per}})
    assert r["valid"]


def test_a_missing_committed_file_voids(ls, tmp_path):
    truths = {k: _task(k, rows=ROWS[k]) for k in ls.TRUTHS}
    out, t42, t44 = _write(ls, tmp_path, truths)
    t44.unlink()
    r = ls.read(out, t42, t44)["reading"]
    assert not r["valid"] and any("do not reproduce" in x for x in r["validity"]["reasons"])


def test_the_reading_refuses_to_overwrite(ls, tmp_path):
    truths = {k: _task(k, rows=ROWS[k]) for k in ls.TRUTHS}
    out, t42, t44 = _write(ls, tmp_path, truths)
    ls.read(out, t42, t44)
    assert (out / "t46.json").exists() and (out / "t46.txt").exists()
    with pytest.raises(SystemExit, match="already holds"):
        ls.read(out, t42, t44)


# ---------------------------------------------------------------------------------------------------------------------
# The run's refusals


def test_the_run_needs_its_row(ls, monkeypatch):
    monkeypatch.setattr(ls.n12, "tries_rows", lambda: {"T44", "T45"})
    monkeypatch.setattr(ls.gc, "n22_passed", lambda path: True)
    with pytest.raises(SystemExit, match="T46 has no row"):
        ls.main(["run", "--truth", "lib3", "--out", "t46"])


def test_the_run_needs_n22(ls, monkeypatch):
    monkeypatch.setattr(ls.n12, "tries_rows", lambda: {"T46"})
    with pytest.raises(SystemExit, match="does not hold a passing N22"):
        ls.main(["run", "--truth", "pooled", "--out", "t46", "--n22", "n22x"])


@pytest.mark.parametrize("truth", ["pooled", "lib1", "lib2", "lib3"])
def test_the_run_refuses_an_existing_task_file(ls, monkeypatch, truth):
    monkeypatch.setattr(ls.n12, "tries_rows", lambda: {"T46"})
    monkeypatch.setattr(ls.gc, "n22_passed", lambda path: True)
    tasks = ls.RESULTS / "t46" / "tasks"
    tasks.mkdir(parents=True)
    (tasks / f"{truth}.json").write_text("{}")
    with pytest.raises(SystemExit, match="already holds"):
        ls.main(["run", "--truth", truth, "--out", "t46"])


def test_the_run_refuses_an_unknown_truth(ls):
    with pytest.raises(SystemExit):
        ls.main(["run", "--truth", "lib4", "--out", "t46"])

